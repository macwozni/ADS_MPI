module my_mpi

contains


   !!! liczy pozycje sasiedniego procesora
   function neighbour(d) result(idx)
      use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ
      use communicators, ONLY: processors
      implicit none
      integer(kind = 4), intent(in), dimension(3) :: d
      integer(kind = 4) :: idx
      integer(kind = 4) :: ix, iy, iz

      ix = MYRANKX + d(1) + 1
      iy = MYRANKY + d(2) + 1
      iz = MYRANKZ + d(3) + 1
      idx = processors(ix, iy, iz)

   end function


   !!! przesyla cala kostke
   subroutine send_piece(items, dst, req, nrcpp)
      use mpi
      implicit none
      real (kind = 8), intent(in) :: items(:)
      integer, intent(in) :: dst, req
      integer(kind = 4), intent(in), dimension(3) :: nrcpp
      integer :: ierr

      call mpi_isend(items, nrcpp(3) * nrcpp(1) * nrcpp(2), &
      MPI_DOUBLE_PRECISION, dst, 0, MPI_COMM_WORLD, req, ierr)

   end subroutine


   !! odbiera cala kostke
   subroutine recv_piece(items, src, req, nrcpp)
      use mpi
      implicit none
      real (kind = 8) :: items(:)
      integer(kind = 4) :: src, req
      integer(kind = 4), dimension(3) :: nrcpp
      integer(kind = 4) :: ierr

      call mpi_irecv(items, nrcpp(3) * nrcpp(1) * nrcpp(2), &
      MPI_DOUBLE_PRECISION, src, 0, MPI_COMM_WORLD, req, ierr)

   end subroutine



   ! -------------------------------------------------------------------
   ! Distributes spline (e.g. solution from previous timestep, or parametrization
   ! approximation for Jacobian calculation) to neighbouring processes. It is
   ! essential that each process possess enough of the solution to be able to
   ! calculate values near the boundary, hence overlapping supports of B-splines
   ! necessitate partial sharing.
   ! -------------------------------------------------------------------
   subroutine DistributeSpline(spline, nrcpp, R)
      use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, NRPROCX, NRPROCY, NRPROCZ
      use mpi
      implicit none
      integer(kind = 4), intent(in), dimension(3) :: nrcpp
      real (kind = 8) :: spline(:,:,:,:)
      real (kind = 8), allocatable :: R(:,:,:,:)
      integer(kind = 4) :: i, j, k, s
      integer(kind = 4) :: request(3 * 3 * 3 * 2), stat(MPI_STATUS_SIZE)
      integer(kind = 4) :: ierr(3 * 3 * 3 * 2)
      integer(kind = 4) :: fierr
      integer(kind = 4) :: dst, src
      integer(kind = 4) :: temp(3)

      s = 1

      ! Right
      if (MYRANKX < NRPROCX - 1) then
         temp = (/1, 0, 0/)
         dst = neighbour(temp)
         call send_piece(spline(:, 2, 2, 2), dst, request(s), nrcpp)
         s = s + 1
      endif
      if (MYRANKX > 0) then
         temp = (/-1, 0, 0/)
         src = neighbour(temp)
         call recv_piece(spline(:, 1, 2, 2), src, request(s), nrcpp)
         s = s + 1
      endif

      do i = 1, s - 1
         call mpi_wait(request(i), stat, fierr)
      enddo
      s = 1

      ! Up
      if (MYRANKY > 0) then
         temp = (/0, -1, 0/)
         dst = neighbour(temp)
         call send_piece(spline(:, 2, 2, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 1, 2, 2), dst, request(s), nrcpp)
         s = s + 1
      endif
      if (MYRANKY < NRPROCY - 1) then
         temp = (/0, 1, 0/)
         src = neighbour(temp)
         call recv_piece(spline(:, 2, 3, 2), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 1, 3, 2), src, request(s), nrcpp)
         s = s + 1
      endif

      do i = 1, s - 1
         call mpi_wait(request(i), stat, fierr)
      enddo
      s = 1

      ! Left
      if (MYRANKX > 0) then
         temp = (/-1, 0, 0/)
         dst = neighbour(temp)
         call send_piece(R(:, 2, 2, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(R(:, 2, 3, 2), dst, request(s), nrcpp)
         s = s + 1
      endif
      if (MYRANKX < NRPROCX - 1) then
         temp = (/1, 0, 0/)
         src = neighbour(temp)
         call recv_piece(R(:, 3, 2, 2), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(R(:, 3, 3, 2), src, request(s), nrcpp)
         s = s + 1
      endif

      do i = 1, s - 1
         call mpi_wait(request(i), stat, fierr)
      enddo
      s = 1

      ! Above
      if (MYRANKZ < NRPROCZ - 1) then
         temp = (/0, 0, 1/)
         dst = neighbour(temp)
         call send_piece(spline(:, 2, 2, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 1, 2, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 1, 3, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 2, 3, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 3, 3, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 3, 2, 2), dst, request(s), nrcpp)
         s = s + 1
      endif
      if (MYRANKZ > 0) then
         temp = (/0, 0, -1/)
         src = neighbour(temp)
         call recv_piece(spline(:, 2, 2, 1), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 1, 2, 1), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 1, 3, 1), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 2, 3, 1), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 3, 3, 1), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 3, 2, 1), src, request(s), nrcpp)
         s = s + 1
      endif

      do i = 1, s - 1
         call mpi_wait(request(i), stat, fierr)
      enddo
      s = 1

      ! Down
      if (MYRANKY < NRPROCY - 1) then
         temp = (/0, 1, 0/)
         dst = neighbour(temp)
         call send_piece(spline(:, 2, 2, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 1, 2, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 3, 2, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 1, 2, 1), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 2, 2, 1), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 3, 2, 1), dst, request(s), nrcpp)
         s = s + 1
      endif
      if (MYRANKY > 0) then
         temp = (/0, -1, 0/)
         src = neighbour(temp)
         call recv_piece(spline(:, 2, 1, 2), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 1, 1, 2), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 3, 1, 2), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 1, 1, 1), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 2, 1, 1), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 3, 1, 1), src, request(s), nrcpp)
         s = s + 1
      endif

      do i = 1, s - 1
         call mpi_wait(request(i), stat, fierr)
      enddo
      s = 1

      ! Below
      if (MYRANKZ > 0) then
         temp = (/0, 0, -1/)
         dst = neighbour(temp)
         call send_piece(spline(:, 1, 1, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 1, 2, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 1, 3, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 2, 1, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 2, 2, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 2, 3, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 3, 1, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 3, 2, 2), dst, request(s), nrcpp)
         s = s + 1
         call send_piece(spline(:, 3, 3, 2), dst, request(s), nrcpp)
         s = s + 1
      endif
      if (MYRANKZ < NRPROCZ - 1) then
         temp = (/0, 0, 1/)
         src = neighbour(temp)
         call recv_piece(spline(:, 1, 1, 3), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 1, 2, 3), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 1, 3, 3), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 2, 1, 3), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 2, 2, 3), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 2, 3, 3), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 3, 1, 3), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 3, 2, 3), src, request(s), nrcpp)
         s = s + 1
         call recv_piece(spline(:, 3, 3, 3), src, request(s), nrcpp)
         s = s + 1
      endif

      do i = 1, s - 1
         call mpi_wait(request(i), stat, fierr)
      enddo
      s = 1


      call mpi_barrier(MPI_COMM_WORLD, fierr)

   end subroutine


   ! -------------------------------------------------------------------
   ! Calculates size of the piece corresponding to process with
   ! specified coordinates. Concretly, number of coefficients.
   !
   ! x, y, z    - coordinates
   ! -------------------------------------------------------------------
   function SizeOfPiece(point, n, p) result (s)
      use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints
      implicit none
      integer(kind = 4), intent(in), dimension(3) :: point
      integer(kind = 4), intent(in), dimension(3) :: n
      integer(kind = 4), intent(in), dimension(3) :: p
      integer(kind = 4) :: s
      integer(kind = 4) :: sx, sy, sz
      integer(kind = 4) :: nrcpp, ibeg, iend
      integer(kind = 4) :: mine, maxe

      call ComputeEndpoints(point(1), NRPROCX, n(1), p(1), nrcpp, ibeg, iend, mine, maxe)
      sx = iend - ibeg + 1
      call ComputeEndpoints(point(2), NRPROCY, n(2), p(2), nrcpp, ibeg, iend, mine, maxe)
      sy = iend - ibeg + 1
      call ComputeEndpoints(point(3), NRPROCZ, n(3), p(3), nrcpp, ibeg, iend, mine, maxe)
      sz = iend - ibeg + 1

      s = sx * sy * sz

   end function


   ! -------------------------------------------------------------------
   ! Gathers full solution at the specified process. It is stored in
   ! 3D array.
   !
   ! at    - process where to gather the solution
   ! part  - part of the solution of each process
   ! full  - full solution, combined from parts
   !
   ! The procedure relies crucially on specific data layout inside
   ! pieces at the end of each iteration.
   !
   ! Expected decomposition structure layout:
   !   Of the pieces: process coordinates = (z, y, x)
   !   Inside pieces: (z, y, x), i.e. x changes fastest
   ! -------------------------------------------------------------------
   subroutine GatherFullSolution(at, part, full, n, p, s)
      use parallelism, ONLY: MYRANK, LINEARINDEX, NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints
      use mpi
      implicit none
      integer(kind = 4), intent(in) :: at
      integer(kind = 4), intent(in), dimension(3) :: n
      integer(kind = 4), intent(in), dimension(3) :: p
      integer(kind = 4), intent(in), dimension(3) :: s
      real (kind = 8), intent(in) :: part(:,:)
      real (kind = 8), intent(out), allocatable :: full(:,:,:)
      real (kind = 8), allocatable :: buffer(:)
      integer(kind = 4) :: recvcounts(0:NRPROCX * NRPROCY * NRPROCZ - 1)
      integer(kind = 4) :: displs(0:NRPROCX * NRPROCY * NRPROCZ - 1)
      integer(kind = 4) :: offset, size
      integer(kind = 4) :: ierr
      integer(kind = 4) :: array_size
      integer(kind = 4), dimension(3) :: begs, ends
      integer(kind = 4) :: mine, maxe
      integer(kind = 4) :: nrcpp
      integer(kind = 4), dimension(3) :: ss
      integer(kind = 4) :: xx, yy, zz
      integer(kind = 4) :: x, y, z
      integer(kind = 4), dimension(3) :: i
      integer(kind = 4) :: idx
      integer(kind = 4), dimension(3) :: tmp

      ! Only the root process needs buffer, but passing unallocated array
      ! is illegal in Fortran, hence we allocate it as array of size 0
      ! in other processes.
      if (MYRANK == at) then
         array_size = (n(1) + 1)*(n(2) + 1)*(n(3) + 1)
         allocate(full(0:n(1), 0:n(2), 0:n(3)))
      else
         array_size = 0
      endif

      allocate(buffer(0:array_size - 1))

      ! Just grab all the pieces and put it in the array one after another,
      ! reordering will be done later at the root.
      offset = 0
      do x = 0, NRPROCX - 1
         do y = 0, NRPROCY - 1
            do z = 0, NRPROCZ - 1
               idx = LinearIndex(x, y, z)
               tmp = (/x, y, z/)
               size = SizeOfPiece(tmp, n, p)
               recvcounts(idx) = size
               displs(idx) = offset
               offset = offset + size
            enddo
         enddo
      enddo

      call mpi_gatherv(part, s(1) * s(2) * s(3), MPI_DOUBLE_PRECISION, buffer, &
      recvcounts, displs, MPI_DOUBLE_PRECISION, at, MPI_COMM_WORLD, ierr)

      ! Reordering of the array at root
      if (MYRANK == at) then
         offset = 0
         do x = 0, NRPROCX - 1
            do y = 0, NRPROCY - 1
               do z = 0, NRPROCZ - 1
                  call ComputeEndpoints(x, NRPROCX, n(1), p(1), nrcpp, begs(1), ends(1), mine, maxe)
                  call ComputeEndpoints(y, NRPROCY, n(2), p(2), nrcpp, begs(2), ends(2), mine, maxe)
                  call ComputeEndpoints(z, NRPROCZ, n(3), p(3), nrcpp, begs(3), ends(3), mine, maxe)
                  ss(1) = ends(1) - begs(1) + 1
                  ss(2) = ends(2) - begs(2) + 1
                  ss(3) = ends(3) - begs(3) + 1

                  do xx = 0, ss(1) - 1
                     do yy = 0, ss(2) - 1
                        do zz = 0, ss(3) - 1
                           i(1) = begs(1) - 1 + xx ! beg_ starts from 1, hence -1
                           i(2) = begs(2) - 1 + yy
                           i(3) = begs(3) - 1 + zz
                           idx = (zz * ss(2) + yy) * ss(1) + xx

                           full(i(1), i(2), i(3)) = buffer(offset + idx)
                        enddo
                     enddo
                  enddo

                  tmp = (/x, y, z/)
                  offset = offset + SizeOfPiece(tmp, n, p)
               enddo
            enddo
         enddo
      endif

      deallocate(buffer)

   end subroutine



   ! -------------------------------------------------------------------
   ! Delinearizes an array
   !
   ! F_lin    - input rank-1 array
   ! F        - output rank-2 array
   ! elems    - first dimension of F
   ! stride   - second dimension of F
   !
   ! Input array like this (s = stride, N = elems):
   !
   !     x11 x12 x13 ... x1s x21 x22 ... x2s... xN1 xN2 ... xNs
   !
   ! is reshaped to the following form:
   !
   !     x11 x12 x13 ... x1s
   !     x21 x22 x23 ... x2s
   !      .   .   .  .    .
   !      .   .   .    .  .
   !     xN1 xN2 xN3 ... xNs
   ! -------------------------------------------------------------------
   subroutine Delinearize(F_lin, F, elems, stride)
      implicit none
      integer(kind = 4), intent(in) :: elems, stride
      real (kind = 8), intent(in) :: F_lin(elems * stride)
      real (kind = 8), intent(out) :: F(elems, stride)
      integer(kind = 4) :: i, a, b

      do i = 1, elems
         a = (i - 1) * stride + 1
         b = i * stride
         F(i,:) = F_lin(a:b)
      enddo

   end subroutine



   ! -------------------------------------------------------------------
   ! Gathers data along one axis to the processes on the corresponding
   ! face of the cube.
   !
   ! F         - input array
   ! F_out     - output array
   ! n         - problem size (total length of the axis)
   ! elems     - length of slice assigned to this process
   ! stride    - total size of each slice layer
   ! dims      - sizes of slices for all processors
   ! shifts    - offsets (linearized) of slices for all processors
   ! comm      - communicator of the axis
   ! ierr      - error code output
   ! -------------------------------------------------------------------
   subroutine Gather(F, F_out, n, elems, stride, dims, shifts, comm, ierr)
      use mpi
      implicit none
      integer(kind = 4), intent(in) :: n, elems, stride, comm
      real (kind = 8), intent(in) :: F(elems, stride)
      integer(kind = 4), intent(in) :: dims(:), shifts(:)
      real (kind = 8), intent(out) :: F_out(n + 1, stride)
      integer(kind = 4), intent(out) :: ierr
      real (kind = 8) :: F_lin(elems * stride), F_out_lin((n + 1) * stride)

      call Linearize(F, F_lin, elems, stride)

      call mpi_gatherv(F_lin, &
      elems * stride, &
      MPI_DOUBLE_PRECISION, &
      F_out_lin, &
      dims, shifts, &
      MPI_DOUBLE_PRECISION, &
      0, comm, ierr)

      call Delinearize(F_out_lin, F_out, n + 1, stride)

   end subroutine



   ! -------------------------------------------------------------------
   ! Scatters computed partial solution along one axis, but does not
   ! delinearize the output (used at the very end of computation)
   !
   ! F         - data to scatter
   ! F_out     - buffer to receive data
   ! n         - problem size
   ! elems     - length of received slice
   ! stride    - total size of each slice layer
   ! dims      - sizes of slices for all processors
   ! shifts    - offsets (linearized) of slices for all processors
   ! comm      - communicator of the axis
   ! ierr      - error code output
   ! -------------------------------------------------------------------
   subroutine Scatter2(F, F_out, n, elems, stride, dims, shifts, comm, ierr)
      use mpi
      implicit none
      integer(kind = 4), intent(in) :: n, elems, stride, comm
      real (kind = 8), intent(in) :: F(n + 1, stride)
      integer(kind = 4), intent(in) :: dims(:), shifts(:)
      real (kind = 8), intent(out) :: F_out(elems * stride)
      integer(kind = 4), intent(out) :: ierr
      real (kind = 8) :: F_lin((n + 1) * stride)

      call Linearize(F, F_lin, n + 1, stride)

      call mpi_scatterv(F_lin, &
      dims, shifts, &
      MPI_DOUBLE_PRECISION, &
      F_out, &
      elems * stride, &
      MPI_DOUBLE_PRECISION, &
      0, comm, ierr)

   end subroutine


   ! -------------------------------------------------------------------
   ! Scatters computed partial solution along one axis
   !
   ! F         - data to scatter
   ! F_out     - buffer to receive data
   ! n         - problem size
   ! elems     - length of received slice
   ! stride    - total size of each slice layer
   ! dims      - sizes of slices for all processors
   ! shifts    - offsets (linearized) of slices for all processors
   ! comm      - communicator of the axis
   ! ierr      - error code output
   ! -------------------------------------------------------------------
   subroutine Scatter(F, F_out, n, elems, stride, dims, shifts, comm, ierr)
      implicit none
      integer(kind = 4), intent(in) :: n, elems, stride, comm
      real (kind = 8), intent(in) :: F(n + 1, stride)
      integer(kind = 4), intent(in) :: dims(:), shifts(:)
      real (kind = 8), intent(out) :: F_out(elems, stride)
      integer(kind = 4), intent(out) :: ierr
      real (kind = 8) :: F_out_lin(elems * stride)
      
      call Scatter2(F, F_out_lin, n, elems, stride, dims, shifts, comm, ierr)
      call Delinearize(F_out_lin, F_out, elems, stride)

   end subroutine


   ! -------------------------------------------------------------------
   ! Broadcasts computed partial solution along one axis
   !
   ! Input:
   ! ------
   ! F        - data to distribute
   ! n        - problem size
   ! elems    - length of received slice
   ! stride   - total size of each slice layer
   ! dims     - sizes of slices for all processors
   ! shifts   - offsets (linearized) of slices for all processors
   ! comm     - communicator of the axis
   !
   ! Output:
   ! ------
   ! F_out    - buffer to receive data
   !
   ! -------------------------------------------------------------------
   subroutine AllGather(F, F_out, n, elems, stride, dims, shifts, comm)
      use mpi
      implicit none
      integer(kind = 4), intent(in) :: n, elems, stride, comm
      real (kind = 8), intent(in) :: F(elems, stride)
      real (kind = 8), intent(out) :: F_out(n + 1, stride)
      integer(kind = 4) :: dims(:), shifts(:)
      real (kind = 8) :: F_lin(elems * stride), F_out_lin((n + 1) * stride)
      integer(kind = 4) :: ierr

      call Linearize(F, F_lin, elems, stride)

      call mpi_allgatherv(F_lin, &
      elems * stride, &
      MPI_DOUBLE_PRECISION, &
      F_out_lin, &
      dims, shifts, &
      MPI_DOUBLE_PRECISION, &
      comm, ierr)

      call Delinearize(F_out_lin, F_out, n + 1, stride)

   end subroutine



   ! -------------------------------------------------------------------
   ! Linearizes and transposes an array.
   !
   ! Input:
   ! ------
   ! F        - input rank-2 array
   ! elems    - first dimension of F
   ! stride   - second dimension of F
   !
   ! Output:
   ! ------
   ! F_lin    - output rank-1 array
   !
   ! For input array like this, where columns are consecutive in memory
   ! (Fortran uses column-major layout), with s = stride, N = elems:
   !
   !     x11 x12 x13 ... x1s
   !     x21 x22 x23 ... x2s
   !      .   .   .  .    .
   !      .   .   .    .  .
   !     xN1 xN2 xN3 ... xNs
   !
   ! output array has the form
   !
   !     x11 x12 x13 ... x1s x21 x22 ... x2s... xN1 xN2 ... xNs
   ! -------------------------------------------------------------------
   subroutine Linearize(F, F_lin, elems, stride)
      implicit none
      integer(kind = 4), intent(in) :: elems, stride
      real (kind = 8), intent(in) :: F(elems, stride)
      real (kind = 8), intent(out) :: F_lin(elems * stride)
      integer(kind = 4) :: i, a, b

      do i = 1, elems
         a = (i - 1) * stride + 1
         b = i * stride
         F_lin(a:b) = F(i,:)
      enddo

   end subroutine





end module
