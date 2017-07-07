module my_mpi
   
contains


!!! liczy pozycje sasiedniego procesora
function neighbour(dx, dy, dz) result(idx)
use parallelism, ONLY : MYRANKX,MYRANKY,MYRANKZ
use communicators, ONLY : processors
implicit none
integer(kind=4), intent(in) :: dx, dy, dz
integer(kind=4) :: idx
integer(kind=4) :: ix, iy, iz

  ix = MYRANKX + dx + 1
  iy = MYRANKY + dy + 1
  iz = MYRANKZ + dz + 1
  idx = processors(ix, iy, iz)

end function


!!! przesyla cala kostke
subroutine send_piece(items, dst, req, nrcppx,nrcppy,nrcppz)
implicit none
include "mpif.h"
real (kind=8), intent(in) :: items(:)
integer, intent(in) :: dst, req
integer(kind=4), intent(in) :: nrcppx,nrcppy,nrcppz
integer :: ierr

  call mpi_isend(items, nrcppz*nrcppx*nrcppy, &
    MPI_DOUBLE_PRECISION, dst, 0, MPI_COMM_WORLD, req, ierr)

end subroutine


! odbiera cala kostke
subroutine recv_piece(items,src,req,nrcppz,nrcppx,nrcppy)
implicit none
include "mpif.h"
real   (kind=8) :: items(:)
integer(kind=4) :: src, req
integer(kind=4) :: nrcppz,nrcppx,nrcppy
integer(kind=4) :: ierr

  call mpi_irecv(items, nrcppz*nrcppx*nrcppy, &
    MPI_DOUBLE_PRECISION, src, 0, MPI_COMM_WORLD, req, ierr)

end subroutine



! -------------------------------------------------------------------
! Distributes spline (e.g. solution from previous timestep, or parametrization
! approximation for Jacobian calculation) to neighbouring processes. It is
! essential that each process possess enough of the solution to be able to
! calculate values near the boundary, hence overlapping supports of B-splines
! necessitate partial sharing.
! -------------------------------------------------------------------
subroutine DistributeSpline(spline,nrcppz,nrcppx,nrcppy,R)
use parallelism, ONLY : MYRANKX,MYRANKY,MYRANKZ,NRPROCX,NRPROCY,NRPROCZ
implicit none
include "mpif.h"
integer(kind=4), intent(in) :: nrcppx,nrcppy,nrcppz
real   (kind=8) :: spline(:,:,:,:)
real (kind=8), allocatable :: R(:,:,:,:)
integer(kind=4) :: i, j, k, s
integer(kind=4) :: request(3*3*3*2), stat(MPI_STATUS_SIZE)
integer(kind=4) :: ierr(3*3*3*2)
integer(kind=4) :: dst, src

  s = 1

  ! Right
  if (MYRANKX < NRPROCX - 1) then
    dst = neighbour(1, 0, 0)
    call send_piece(spline(:,2,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif
  if (MYRANKX > 0) then
    src = neighbour(-1, 0, 0)
    call recv_piece(spline(:,1,2,2), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Up
  if (MYRANKY > 0) then
    dst = neighbour(0, -1, 0)
    call send_piece(spline(:,2,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,1,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif
  if (MYRANKY < NRPROCY - 1) then
    src = neighbour(0, 1, 0)
    call recv_piece(spline(:,2,3,2), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,1,3,2), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Left
  if (MYRANKX > 0) then
    dst = neighbour(-1, 0, 0)
    call send_piece(R(:,2,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(R(:,2,3,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif
  if (MYRANKX < NRPROCX - 1) then
    src = neighbour(1, 0, 0)
    call recv_piece(R(:,3,2,2), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(R(:,3,3,2), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Above
  if (MYRANKZ < NRPROCZ - 1) then
    dst = neighbour(0, 0, 1)
    call send_piece(spline(:,2,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,1,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,1,3,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,2,3,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,3,3,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,3,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif
  if (MYRANKZ > 0) then
    src = neighbour(0, 0, -1)
    call recv_piece(spline(:,2,2,1), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,1,2,1), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,1,3,1), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,2,3,1), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,3,3,1), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,3,2,1), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Down
  if (MYRANKY < NRPROCY - 1) then
    dst = neighbour(0, 1, 0)
    call send_piece(spline(:,2,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,1,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,3,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,1,2,1), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,2,2,1), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,3,2,1), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif
  if (MYRANKY > 0) then
    src = neighbour(0, -1, 0)
    call recv_piece(spline(:,2,1,2), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,1,1,2), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,3,1,2), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,1,1,1), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,2,1,1), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,3,1,1), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Below
  if (MYRANKZ > 0) then
    dst = neighbour(0, 0, -1)
    call send_piece(spline(:,1,1,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,1,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,1,3,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,2,1,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,2,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,2,3,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,3,1,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,3,2,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call send_piece(spline(:,3,3,2), dst, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif
  if (MYRANKZ < NRPROCZ - 1) then
    src = neighbour(0, 0, 1)
    call recv_piece(spline(:,1,1,3), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,1,2,3), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,1,3,3), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,2,1,3), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,2,2,3), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,2,3,3), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,3,1,3), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,3,2,3), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
    call recv_piece(spline(:,3,3,3), src, request(s),nrcppx,nrcppy,nrcppz)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1


  call mpi_barrier(MPI_COMM_WORLD,ierr)

end subroutine


! -------------------------------------------------------------------
! Calculates size of the piece corresponding to process with
! specified coordinates. Concretly, number of coefficients.
!
! x, y, z    - coordinates
! -------------------------------------------------------------------
function SizeOfPiece(x,y,z,nx,ny,nz,px,py,pz) result (s)
use parallelism, ONLY : NRPROCX,NRPROCY,NRPROCZ
use utils, ONLY : ComputeEndpoints
implicit none
integer(kind=4), intent(in) :: x, y, z
integer(kind=4), intent(in) :: nx, ny, nz
integer(kind=4), intent(in) :: px, py, pz
integer(kind=4) :: s
integer(kind=4) :: sx, sy, sz
integer(kind=4) :: nrcpp, ibeg, iend
integer(kind=4) :: mine, maxe

  call ComputeEndpoints(x, NRPROCX, nx, px, nrcpp, ibeg, iend, mine, maxe)
  sx = iend - ibeg + 1
  call ComputeEndpoints(y, NRPROCY, ny, py, nrcpp, ibeg, iend, mine, maxe)
  sy = iend - ibeg + 1
  call ComputeEndpoints(z, NRPROCZ, nz, pz, nrcpp, ibeg, iend, mine, maxe)
  sz = iend - ibeg + 1

  s = sx * sy * sz

end function


!!!!! przeniesc do my_mpi
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
subroutine GatherFullSolution(at,part,full,nx,ny,nz,px,py,pz,sx,sy,sz)
use parallelism, ONLY : MYRANK,LINEARINDEX,NRPROCX,NRPROCY,NRPROCZ
use utils, ONLY : ComputeEndpoints
implicit none
include "mpif.h"
integer(kind=4), intent(in) :: at
integer(kind=4), intent(in) :: nx,ny,nz
integer(kind=4), intent(in) :: px,py,pz
integer(kind=4), intent(in) :: sx,sy,sz
real   (kind=8), intent(in) :: part(:,:)
real   (kind=8), intent(out), allocatable :: full(:,:,:)
real   (kind=8), allocatable :: buffer(:)
integer(kind=4) :: recvcounts(0:NRPROCX*NRPROCY*NRPROCZ-1)
integer(kind=4) :: displs(0:NRPROCX*NRPROCY*NRPROCZ-1)
integer(kind=4) :: x, y, z
integer(kind=4) :: offset, size
integer(kind=4) :: ierr
integer(kind=4) :: array_size
integer(kind=4) :: begx, begy, begz, endx, endy, endz
integer(kind=4) :: mine, maxe
integer(kind=4) :: nrcpp
integer(kind=4) :: ssx, ssy, ssz
integer(kind=4) :: xx, yy, zz
integer(kind=4) :: ix, iy, iz, idx

  ! Only the root process needs buffer, but passing unallocated array
  ! is illegal in Fortran, hence we allocate it as array of size 0
  ! in other processes.
  if (MYRANK == at) then
    array_size = (nx+1)*(ny+1)*(nz+1)
    allocate(full(0:nx,0:ny,0:nz))
  else
    array_size = 0
  endif

  allocate(buffer(0:array_size-1))

  ! Just grab all the pieces and put it in the array one after another,
  ! reordering will be done later at the root.
  offset = 0
  do x = 0, NRPROCX-1
    do y = 0, NRPROCY-1
      do z = 0, NRPROCZ-1
        idx = LinearIndex(x, y, z)
        size = SizeOfPiece(x,y,z,nx,ny,nz,px,py,pz)
        recvcounts(idx) = size
        displs(idx) = offset
        offset = offset + size
      enddo
    enddo
  enddo

  call mpi_gatherv(part, sx*sy*sz, MPI_DOUBLE_PRECISION, buffer, &
    recvcounts, displs, MPI_DOUBLE_PRECISION, at, MPI_COMM_WORLD, ierr)

  ! Reordering of the array at root
  if (MYRANK == at) then
    offset = 0
    do x = 0, NRPROCX-1
      do y = 0, NRPROCY-1
        do z = 0, NRPROCZ-1
          call ComputeEndpoints(x, NRPROCX, nx, px, nrcpp, begx, endx, mine, maxe)
          call ComputeEndpoints(y, NRPROCY, ny, py, nrcpp, begy, endy, mine, maxe)
          call ComputeEndpoints(z, NRPROCZ, nz, pz, nrcpp, begz, endz, mine, maxe)
          ssx = endx - begx + 1
          ssy = endy - begy + 1
          ssz = endz - begz + 1

          do xx = 0, ssx-1
            do yy = 0, ssy-1
              do zz = 0, ssz-1
                ix = begx - 1 + xx   ! beg_ starts from 1, hence -1
                iy = begy - 1 + yy
                iz = begz - 1 + zz
                idx = (zz * ssy + yy) * ssx + xx

                full(ix, iy, iz) = buffer(offset + idx)
              enddo
            enddo
          enddo

          offset = offset + SizeOfPiece(x,y,z,nx,ny,nz,px,py,pz)
        enddo
      enddo
    enddo
  endif

  deallocate(buffer)

end subroutine


end module
