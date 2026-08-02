!------------------------------------------------------------------------------
!
! MODULE: my_mpi
!
! DESCRIPTION:
!> @file my_mpi.F90
!> @brief Module providing MPI wrapper routines for neighbour exchange,
!> directional gather/scatter operations, and global solution assembly.
!>
!> @details
!> This module groups MPI-oriented helper procedures used throughout the
!> project to:
!> - identify neighbouring processes in the logical process grid,
!> - exchange overlapping spline-coefficient blocks between neighbouring
!>   subdomains,
!> - gather and scatter directional right-hand-side slices on
!>   one-dimensional communicator fibres,
!> - assemble the full global solution from distributed local pieces,
!> - convert between linear and two-dimensional storage layouts used by
!>   communication kernels.
!>
!> The module is tightly coupled with:
!> - \ref parallelism for process-grid metadata and decomposition tools,
!> - \ref communicators for logical process-cube rank tables,
!> - the ADS workflow for redistribution and directional solves.
!
!------------------------------------------------------------------------------
module my_mpi

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes the linear MPI rank of a neighbouring process in the
!> logical process grid.
!>
!> @details
!> This function adds the directional offset \p d to the current
!> process-grid coordinates and uses the global process-cube table
!> \ref processors from module \ref communicators to retrieve the
!> corresponding linear rank.
!>
!> The input offset is interpreted componentwise as a displacement in the
!> \f$x\f$-, \f$y\f$-, and \f$z\f$-directions relative to the current
!> process.
!
! Input:
! ------
!> @param[in] d
!> Integer displacement vector identifying the neighbouring process.
!
! Output:
! -------
!> @return idx
!> Linear MPI rank of the selected neighbouring process.
!
! Notes:
! ------
!> @warning
!> The routine assumes that the requested neighbour lies inside the valid
!> logical process grid.
!
!---------------------------------------------------------------------------
function neighbour(d) result(idx)
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ
   use communicators, ONLY: processors
   implicit none
!> @brief Directional displacement from the current process.
   integer(kind=4), dimension(:), intent(in) :: d
!> @brief Linear rank of the neighbouring process.
   integer(kind=4) :: idx
!> @brief Shifted logical process-grid coordinates.
   integer(kind=4) :: ix, iy, iz

   ix = MYRANKX + d(1) + 1
   iy = MYRANKY + d(2) + 1
   iz = MYRANKZ + d(3) + 1
   idx = processors(ix, iy, iz)

end function neighbour

!---------------------------------------------------------------------------
!> @brief Exchanges exactly the coefficient boxes required for reconstruction.
!>
!> The communication plan is prepared once by `AllocateADSdata`.  It contains
!> intersections between the trial-coefficient box owned by every source and
!> the arbitrarily wide coefficient halo required by every destination.
!> Consequently, this routine has no nearest-neighbour restriction and sends
!> neither complete remote blocks nor a globally replicated solution.
!---------------------------------------------------------------------------
subroutine DistributeSpline(part, ads_trial, ads_data)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK, NRPROC
   use mpi
   implicit none
!> @brief Locally owned trial coefficients in normalized X-Y-Z order.
   real(kind=8), dimension(:, :), intent(in) :: part
!> @brief Trial-space ownership metadata.
   type(ADS_setup), intent(in) :: ads_trial
!> @brief Halo plan, reusable communication buffers, and destination halo.
   type(ADS_compute_data), intent(inout) :: ads_data
   integer(kind=4), parameter :: HALO_TAG = 731
   integer(kind=4) :: peer, count, pos, column
   integer(kind=4) :: x, y, z, nreq, ierr
   integer(kind=4), dimension(3) :: owned_begin

   owned_begin = ads_trial%ibeg - 1

   ! Pack each source/destination intersection from the local coefficient
   ! block.  X is the fastest-changing coordinate in both layouts.
   do peer = 1, NRPROC
      count = ads_data%halo_send_count(peer)
      if (count == 0) cycle
      pos = ads_data%halo_send_displ(peer) + 1
      do z = ads_data%halo_send_begin(3, peer), ads_data%halo_send_end(3, peer)
         do y = ads_data%halo_send_begin(2, peer), ads_data%halo_send_end(2, peer)
            column = y - owned_begin(2) + &
                     (z - owned_begin(3))*ads_trial%s(2) + 1
            do x = ads_data%halo_send_begin(1, peer), ads_data%halo_send_end(1, peer)
               ads_data%halo_send_buffer(pos) = part(x - owned_begin(1) + 1, column)
               pos = pos + 1
            end do
         end do
      end do
   end do

   nreq = 0

   ! Post all receives before sends.  The self-intersection is copied without
   ! entering MPI; every other peer is independent of its distance in the
   ! logical process grid.
   do peer = 1, NRPROC
      count = ads_data%halo_recv_count(peer)
      if (count == 0 .or. peer == MYRANK + 1) cycle
      nreq = nreq + 1
      pos = ads_data%halo_recv_displ(peer) + 1
      call mpi_irecv(ads_data%halo_recv_buffer(pos:pos + count - 1), count, &
                     MPI_DOUBLE_PRECISION, peer - 1, HALO_TAG, MPI_COMM_WORLD, &
                     ads_data%halo_requests(nreq), ierr)
   end do

   do peer = 1, NRPROC
      count = ads_data%halo_send_count(peer)
      if (count == 0 .or. peer == MYRANK + 1) cycle
      nreq = nreq + 1
      pos = ads_data%halo_send_displ(peer) + 1
      call mpi_isend(ads_data%halo_send_buffer(pos:pos + count - 1), count, &
                     MPI_DOUBLE_PRECISION, peer - 1, HALO_TAG, MPI_COMM_WORLD, &
                     ads_data%halo_requests(nreq), ierr)
   end do

   peer = MYRANK + 1
   count = ads_data%halo_recv_count(peer)
   if (count > 0) then
      ads_data%halo_recv_buffer( &
         ads_data%halo_recv_displ(peer) + 1:ads_data%halo_recv_displ(peer) + count) = &
      ads_data%halo_send_buffer( &
         ads_data%halo_send_displ(peer) + 1:ads_data%halo_send_displ(peer) + count)
   end if

   if (nreq > 0) call mpi_waitall(nreq, ads_data%halo_requests, &
                                  ads_data%halo_statuses, ierr)

   ! Unpack the received intersections into one compact global-indexed halo.
   do peer = 1, NRPROC
      count = ads_data%halo_recv_count(peer)
      if (count == 0) cycle
      pos = ads_data%halo_recv_displ(peer) + 1
      do z = ads_data%halo_recv_begin(3, peer), ads_data%halo_recv_end(3, peer)
         do y = ads_data%halo_recv_begin(2, peer), ads_data%halo_recv_end(2, peer)
            do x = ads_data%halo_recv_begin(1, peer), ads_data%halo_recv_end(1, peer)
               ads_data%R(x - ads_data%halo_begin(1) + 1, &
                          y - ads_data%halo_begin(2) + 1, &
                          z - ads_data%halo_begin(3) + 1, 1) = &
                  ads_data%halo_recv_buffer(pos)
               pos = pos + 1
            end do
         end do
      end do
   end do

end subroutine DistributeSpline

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes the number of coefficients owned by a process at a
!> specified logical process-grid coordinate.
!>
!> @details
!> This function determines the owned basis-function range in each
!> coordinate direction for the process identified by \p point and then
!> returns the product of the corresponding local sizes.
!>
!> The ownership ranges are computed using \ref ComputeEndpoints from
!> module \ref parallelism.
!
! Input:
! ------
!> @param[in] point
!> Logical process-grid coordinates of the selected process.
!>
!> @param[in] n
!> Numbers of basis functions minus one in the three directions.
!>
!> @param[in] p
!> Polynomial degrees in the three directions.
!
! Output:
! -------
!> @return s
!> Number of coefficients stored by the selected process.
!
!---------------------------------------------------------------------------
function SizeOfPiece(point, n, p) result(s)
   use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints
   implicit none
!> @brief Logical process-grid coordinates of the selected process.
   integer(kind=4), intent(in), dimension(:) :: point
!> @brief Numbers of basis functions minus one.
   integer(kind=4), intent(in), dimension(:) :: n
!> @brief Polynomial degrees.
   integer(kind=4), intent(in), dimension(:) :: p
!> @brief Number of coefficients owned by the selected process.
   integer(kind=4) :: s
!> @brief Owned sizes in the three directions.
   integer(kind=4) :: sx, sy, sz
!> @brief Auxiliary ownership data returned by `ComputeEndpoints`.
   integer(kind=4) :: nrcpp, ibeg, iend
!> @brief Auxiliary element-range data returned by `ComputeEndpoints`.
   integer(kind=4) :: mine, maxe

   call ComputeEndpoints(point(1), NRPROCX, n(1), p(1), nrcpp, ibeg, iend, mine, maxe)
   sx = iend - ibeg + 1
   call ComputeEndpoints(point(2), NRPROCY, n(2), p(2), nrcpp, ibeg, iend, mine, maxe)
   sy = iend - ibeg + 1
   call ComputeEndpoints(point(3), NRPROCZ, n(3), p(3), nrcpp, ibeg, iend, mine, maxe)
   sz = iend - ibeg + 1

   s = sx*sy*sz

end function SizeOfPiece

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Gathers the full distributed solution to one selected root
!> process and stores it as a three-dimensional array.
!>
!> @details
!> This routine assembles the global spline coefficient field from the
!> local solution parts \p part owned by all processes.
!>
!> The procedure operates in two stages:
!> - first, all local pieces are collected consecutively into a linear
!>   buffer by `mpi_gatherv`,
!> - second, on the root process \p at, the gathered linear buffer is
!>   reordered into the three-dimensional array \p full according to the
!>   logical process decomposition and local ownership sizes.
!>
!> The local piece sizes are obtained from \ref SizeOfPiece, whereas the
!> ownership ranges required for final placement are computed through
!> \ref ComputeEndpoints.
!
! Input:
! ------
!> @param[in] at
!> Root process rank on which the full solution is assembled.
!>
!> @param[in] part
!> Local part of the distributed solution.
!>
!> @param[in] n
!> Numbers of basis functions minus one in the three directions.
!>
!> @param[in] p
!> Polynomial degrees in the three directions.
!>
!> @param[in] s
!> Owned sizes of the local solution block in the three directions.
!
! Output:
! -------
!> @param[out] full
!> Full gathered solution stored as a three-dimensional tensor on the
!> root process.
!
! Notes:
! ------
!> @note
!> Only the root process allocates the full three-dimensional output
!> array with nonzero size.
!
!> @note
!> The routine assumes the specific tensor-product storage layout
!> documented in the source comments: both process blocks and entries
!> within each block are ordered with the first direction changing
!> fastest.
!
!---------------------------------------------------------------------------
subroutine GatherFullSolution(at, part, full, n, p, s)
   use parallelism, ONLY: MYRANK, LINEARINDEX, NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints
   use mpi
   implicit none
   !> @brief Root process rank.
   integer(kind=4), intent(in) :: at
   !> @brief Numbers of basis functions minus one.
   integer(kind=4), dimension(:), intent(in) :: n
   !> @brief Polynomial degrees.
   integer(kind=4), dimension(:), intent(in) :: p
   !> @brief Owned sizes of the local solution block.
   integer(kind=4), dimension(:), intent(in) :: s
   !> @brief Local piece of the distributed solution.
   real(kind=8), intent(in) :: part(:, :)
   !> @brief Full gathered solution on the root process.
   real(kind=8), allocatable, dimension(:, :, :), intent(out) :: full
   !> @brief Linear receive buffer used by `mpi_gatherv`.
   real(kind=8), allocatable :: buffer(:)
   !> @brief Numbers of coefficients received from all processes.
   integer(kind=4), dimension(0:NRPROCX*NRPROCY*NRPROCZ - 1) :: recvcounts
   !> @brief Displacements of gathered coefficient blocks.
   integer(kind=4), dimension(0:NRPROCX*NRPROCY*NRPROCZ - 1) :: displs
   !> @brief Running offset and local piece size.
   integer(kind=4) :: offset, msize
   !> @brief Minimum number of coefficients for OpenMP root-side unpacking.
   integer(kind=4), parameter :: GATHER_FULL_OMP_THRESHOLD = 262144
   !> @brief Offset and size of one gathered process block.
   integer(kind=4) :: block_offset, block_size
   !> @brief MPI return code.
   integer(kind=4) :: ierr
   !> @brief Total number of gathered coefficients on the root process.
   integer(kind=4) :: array_size
   !> @brief Ownership bounds of a selected process block.
   integer(kind=4), dimension(3) :: begs, ends
   !> @brief Auxiliary element-range data returned by `ComputeEndpoints`.
   integer(kind=4) :: mine, maxe
   !> @brief Auxiliary nominal columns-per-processor value.
   integer(kind=4) :: nrcpp
   !> @brief Local sizes of a selected process block.
   integer(kind=4), dimension(3) :: ss
   !> @brief Loop counters inside one process block.
   integer(kind=4) :: xx, yy, zz
   !> @brief Logical process-grid coordinates.
   integer(kind=4) :: x, y, z
   !> @brief Linear rank and local linearized index.
   integer(kind=4) :: rank_idx, local_idx
   !> @brief Temporary process-grid coordinate triplet.
   integer(kind=4), dimension(3) :: tmp

! Only the root process needs buffer, but passing unallocated array
! is illegal in Fortran, hence we allocate it as array of size 0
! in other processes.
   if (MYRANK == at) then
      array_size = (n(1) + 1)*(n(2) + 1)*(n(3) + 1)
      allocate (full(0:n(1), 0:n(2), 0:n(3)))
   else
      array_size = 0
   end if

   allocate (buffer(0:array_size - 1))

! Just grab all the pieces and put it in the array one after another,
! reordering will be done later at the root.
   offset = 0
   do x = 0, NRPROCX - 1
      do y = 0, NRPROCY - 1
         do z = 0, NRPROCZ - 1
            rank_idx = LinearIndex(x, y, z)
            tmp = (/x, y, z/)
            msize = SizeOfPiece(tmp, n, p)
            recvcounts(rank_idx) = msize
            displs(rank_idx) = offset
            offset = offset + msize
         end do
      end do
   end do

   call mpi_gatherv(part, s(1)*s(2)*s(3), MPI_DOUBLE_PRECISION, buffer, &
                     recvcounts, displs, MPI_DOUBLE_PRECISION, at, MPI_COMM_WORLD, ierr)

! Reordering of the array at root
   if (MYRANK == at) then
      do x = 0, NRPROCX - 1
         do y = 0, NRPROCY - 1
            do z = 0, NRPROCZ - 1
               rank_idx = LinearIndex(x, y, z)
               block_offset = displs(rank_idx)
               block_size = recvcounts(rank_idx)
               call ComputeEndpoints(x, NRPROCX, n(1), p(1), nrcpp, begs(1), ends(1), mine, maxe)
               call ComputeEndpoints(y, NRPROCY, n(2), p(2), nrcpp, begs(2), ends(2), mine, maxe)
               call ComputeEndpoints(z, NRPROCZ, n(3), p(3), nrcpp, begs(3), ends(3), mine, maxe)
               ss(1) = ends(1) - begs(1) + 1
               ss(2) = ends(2) - begs(2) + 1
               ss(3) = ends(3) - begs(3) + 1

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(xx,yy,zz,local_idx) SCHEDULE(STATIC) &
!$OMP IF(block_size > GATHER_FULL_OMP_THRESHOLD)
               do zz = 0, ss(3) - 1
                  do yy = 0, ss(2) - 1
                     do xx = 0, ss(1) - 1
                        local_idx = (zz*ss(2) + yy)*ss(1) + xx
                        full(begs(1) - 1 + xx, begs(2) - 1 + yy, begs(3) - 1 + zz) = &
                              buffer(block_offset + local_idx)
                     end do
                  end do
               end do
!$OMP END PARALLEL DO
            end do
         end do
      end do
   end if

   deallocate (buffer)

end subroutine GatherFullSolution

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a linearized one-dimensional buffer into a
!> two-dimensional array.
!>
!> @details
!> This routine reconstructs a rank-two array \p F from the linear
!> storage buffer \p F_lin by assigning consecutive blocks of length
!> \p stride to successive rows of \p F.
!>
!> The mapping corresponds to the inverse of \ref Linearize.
!
! Input:
! ------
!> @param[in] F_lin
!> Input one-dimensional array containing linearized row data.
!>
!> @param[in] elems
!> First dimension of the output array.
!>
!> @param[in] stride
!> Second dimension of the output array.
!
! Output:
! -------
!> @param[out] F
!> Two-dimensional reconstructed array.
!
!---------------------------------------------------------------------------
subroutine Delinearize(F_lin, F, elems, stride)
   implicit none
!> @brief Dimensions of the output array.
   integer(kind=4), intent(in) :: elems, stride
!> @brief Input linearized data buffer.
   real(kind=8), dimension(:), intent(in) :: F_lin
!> @brief Reconstructed two-dimensional array.
   real(kind=8), dimension(:, :), intent(out) :: F
!> @brief Loop counter and segment bounds inside the linear buffer.
   integer(kind=4) :: i, a, b

   do i = 1, elems
      a = (i - 1)*stride + 1
      b = i*stride
      F(i, :) = F_lin(a:b)
   end do

end subroutine Delinearize

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Gathers distributed data along one logical axis onto the face
!> processes associated with the corresponding communicator.
!>
!> @details
!> This routine linearizes the local two-dimensional input array \p F,
!> gathers all contributions in communicator \p comm by `mpi_gatherv`,
!> and reconstructs the gathered result into the two-dimensional array
!> \p F_out.
!>
!> The vectors \p dims and \p shifts provide the receive counts and
!> displacements describing the distribution of directional slices among
!> processes.
!
! Input:
! ------
!> @param[in] F
!> Local input array to be gathered.
!>
!> @param[in] n
!> Global problem size in the gathered direction, minus one.
!>
!> @param[in] elems
!> Number of locally owned slices in the gathered direction.
!>
!> @param[in] stride
!> Total size of one transverse slice.
!>
!> @param[in] dims
!> Receive counts for all processes in the communicator.
!>
!> @param[in] shifts
!> Receive displacements for all processes in the communicator.
!>
!> @param[in] comm
!> Communicator spanning the selected logical process-grid fibre.
!
! Output:
! -------
!> @param[out] F_out
!> Gathered two-dimensional array.
!>
!> @param[out] ierr
!> MPI return code.
!
!---------------------------------------------------------------------------
subroutine Gather(F, F_out, n, elems, stride, dims, shifts, comm, ierr)
   use mpi
   implicit none
!> @brief Global size, local size, slice stride, and communicator handle.
   integer(kind=4), intent(in) :: n, elems, stride, comm
!> @brief Local input array.
   real(kind=8), dimension(:, :), intent(in) :: F
!> @brief Receive counts and displacements for communicator-local gathering.
   integer(kind=4), dimension(:), intent(in) :: dims, shifts
!> @brief Gathered output array.
   real(kind=8), dimension(:, :), intent(out) :: F_out
!> @brief MPI return code.
   integer(kind=4), intent(out) :: ierr
!> @brief Linearized input and gathered output buffers.
   real(kind=8), dimension(:), allocatable :: F_lin, F_out_lin

   allocate (F_lin(elems*stride))
   allocate (F_out_lin((n + 1)*stride))
   call Linearize(F, F_lin, elems, stride)

   call mpi_gatherv(F_lin, &
                     elems*stride, &
                     MPI_DOUBLE_PRECISION, &
                     F_out_lin, &
                     dims, shifts, &
                     MPI_DOUBLE_PRECISION, &
                     0, comm, ierr)

   call Delinearize(F_out_lin, F_out, n + 1, stride)

   if (allocated(F_lin)) deallocate (F_lin)
   if (allocated(F_out_lin)) deallocate (F_out_lin)

end subroutine Gather

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Scatters directional data from face processes back to all
!> processes of a communicator fibre.
!>
!> @details
!> This routine linearizes the input array \p F on the root process of
!> communicator \p comm, scatters the data using `mpi_scatterv`, and
!> reconstructs the received local piece into the two-dimensional array
!> \p F_out.
!>
!> The vectors \p dims and \p shifts provide the send counts and
!> displacements describing the directional partition of the global data.
!
! Input:
! ------
!> @param[in] F
!> Directional array to be scattered.
!>
!> @param[in] n
!> Global problem size in the scattered direction, minus one.
!>
!> @param[in] elems
!> Number of locally received slices.
!>
!> @param[in] stride
!> Total size of one transverse slice.
!>
!> @param[in] dims
!> Send counts for all processes in the communicator.
!>
!> @param[in] shifts
!> Send displacements for all processes in the communicator.
!>
!> @param[in] comm
!> Communicator spanning the selected logical process-grid fibre.
!
! Output:
! -------
!> @param[out] F_out
!> Local received array after scatter.
!>
!> @param[out] ierr
!> MPI return code.
!
!---------------------------------------------------------------------------
subroutine Scatter(F, F_out, n, elems, stride, dims, shifts, comm, ierr)
   use mpi
   implicit none
!> @brief Global size, local size, slice stride, and communicator handle.
   integer(kind=4), intent(in) :: n, elems, stride, comm
!> @brief Global input array on the communicator root.
   real(kind=8), dimension(:, :), intent(in) :: F
!> @brief Send counts and displacements for communicator-local scattering.
   integer(kind=4), dimension(:), intent(in) :: dims, shifts
!> @brief Local output array after scatter.
   real(kind=8), dimension(:, :), intent(out) :: F_out
!> @brief MPI return code.
   integer(kind=4), intent(out) :: ierr
!> @brief Linearized received buffer.
   real(kind=8), allocatable, dimension(:) :: F_out_lin
!> @brief Linearized global input buffer.
   real(kind=8), allocatable, dimension(:) :: F_lin

   allocate (F_out_lin(elems*stride))
   allocate (F_lin((n + 1)*stride))
   call Linearize(F, F_lin, n + 1, stride)

   call mpi_scatterv(F_lin, &
                     dims, shifts, &
                     MPI_DOUBLE_PRECISION, &
                     F_out_lin, &
                     elems*stride, &
                     MPI_DOUBLE_PRECISION, &
                     0, comm, ierr)

   call Delinearize(F_out_lin, F_out, elems, stride)

   if (allocated(F_lin)) deallocate (F_lin)
   if (allocated(F_out_lin)) deallocate (F_out_lin)

end subroutine Scatter

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Performs an all-gather of directional data along one logical
!> communicator fibre.
!>
!> @details
!> This routine linearizes the local two-dimensional input array \p F,
!> gathers all local pieces from communicator \p comm by
!> `mpi_allgatherv`, and reconstructs the resulting full directional
!> array into \p F_out.
!>
!> Unlike \ref Gather, the assembled result is returned on every process
!> of the communicator.
!
! Input:
! ------
!> @param[in] F
!> Local input array.
!>
!> @param[in] n
!> Global problem size in the gathered direction, minus one.
!>
!> @param[in] elems
!> Number of locally owned slices.
!>
!> @param[in] stride
!> Total size of one transverse slice.
!>
!> @param[in] dims
!> Receive counts for all processes in the communicator.
!>
!> @param[in] shifts
!> Receive displacements for all processes in the communicator.
!>
!> @param[in] comm
!> Communicator spanning the selected logical process-grid fibre.
!
! Output:
! -------
!> @param[out] F_out
!> Full gathered array available on every process in the communicator.
!
!---------------------------------------------------------------------------
subroutine AllGather(F, F_out, n, elems, stride, dims, shifts, comm)
   use mpi
   implicit none
!> @brief Global size, local size, slice stride, and communicator handle.
   integer(kind=4), intent(in) :: n, elems, stride, comm
!> @brief Local input array.
   real(kind=8), dimension(:, :), intent(in) :: F
!> @brief Full gathered array.
   real(kind=8), dimension(:, :), intent(out) :: F_out
!> @brief Receive counts and displacements for communicator-local all-gather.
   integer(kind=4), dimension(:) :: dims, shifts
!> @brief Linearized local input buffer.
   real(kind=8), dimension(elems*stride) :: F_lin
!> @brief Linearized gathered output buffer.
   real(kind=8), dimension((n + 1)*stride) :: F_out_lin
!> @brief MPI return code.
   integer(kind=4) :: ierr

   call Linearize(F, F_lin, elems, stride)

   call mpi_allgatherv(F_lin, &
                        elems*stride, &
                        MPI_DOUBLE_PRECISION, &
                        F_out_lin, &
                        dims, shifts, &
                        MPI_DOUBLE_PRECISION, &
                        comm, ierr)

   call Delinearize(F_out_lin, F_out, n + 1, stride)

end subroutine AllGather

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Linearizes a two-dimensional array into a one-dimensional
!> communication buffer.
!>
!> @details
!> This routine stores successive rows of the two-dimensional input array
!> \p F as consecutive segments inside the one-dimensional output buffer
!> \p F_lin.
!>
!> The mapping is the inverse of \ref Delinearize and is used by the MPI
!> gather/scatter wrappers implemented in this module.
!
! Input:
! ------
!> @param[in] F
!> Input two-dimensional array.
!>
!> @param[in] elems
!> First dimension of the input array.
!>
!> @param[in] stride
!> Second dimension of the input array.
!
! Output:
! -------
!> @param[out] F_lin
!> Linearized one-dimensional communication buffer.
!
!---------------------------------------------------------------------------
subroutine Linearize(F, F_lin, elems, stride)
   implicit none
!> @brief Dimensions of the input array.
   integer(kind=4), intent(in) :: elems, stride
!> @brief Input two-dimensional array.
   real(kind=8), dimension(:, :), intent(in) :: F
!> @brief Linearized output buffer.
   real(kind=8), dimension(:), intent(out) :: F_lin
!> @brief Loop counter and segment bounds inside the output buffer.
   integer(kind=4) :: i, a, b

   do i = 1, elems
      a = (i - 1)*stride + 1
      b = i*stride
      F_lin(a:b) = F(i, :)
   end do

end subroutine Linearize

end module my_mpi
