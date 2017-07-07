module utils

implicit none

contains


!!!! parallelism
! -------------------------------------------------------------------
! Calculates the range of the direction that is assigned to processor
! with specified rank.
!
! Input:
! ------
! rank    - rank of the current process in this direction
! nrproc  - # of processors for this direction
! n       - size of the problem
! p       - order of polynomial basis
!
! Output:
! -------
! nrcpp   - # of columns per processor
! ibeg    - index of first assigned slice
! iend    - index of last assigned slice
! mine    - index of first element corresponding to the assigned slice
! maxe    - index of last element corresponding to the assigned slice
! -------------------------------------------------------------------
subroutine ComputeEndpoints(rank, nrproc, n, p, nrcpp, ibeg, iend, mine, maxe)
implicit none
integer(kind=4), intent(in) :: rank, nrproc, n, p
integer(kind=4), intent(out):: nrcpp, ibeg, iend, mine, maxe
integer(kind=4) :: elems

elems = n + 1 - p
nrcpp = (n+1 + nrproc-1) / nrproc
ibeg = nrcpp * rank + 1

if(rank == nrproc-1)then
  iend = n+1
else
  iend = nrcpp * (rank + 1)
endif

mine = max(ibeg - p - 1, 1)
maxe = min(iend, elems)

end subroutine


!!!!! parallelism
! -------------------------------------------------------------------
! Calculates sizes and ranges of slices for each processor in
! given direction.
!
! dims     - sizes of 'slices' of dimension
! shifts   - offsets of slices
! nrcpp    - # of columns per processor
! stride   - size of full, 3d slice
! n        - size of the problem
! nrproc   - # of processors for this direction
! -------------------------------------------------------------------
subroutine FillDimVector(dims, shifts, nrcpp, stride, n, nrproc)
implicit none
integer(kind=4), intent(in) :: nrcpp, stride, n, nrproc
integer(kind=4), allocatable, intent(out) :: dims(:), shifts(:)
integer(kind=4) :: i

  allocate(dims(nrproc))
  allocate(shifts(nrproc))

  shifts = 0
  dims = 0

  do i = 1,nrproc-1
    dims(i) = nrcpp * stride
    if (i > 1) shifts(i) = shifts(i-1) + dims(i-1)
  enddo

  if (nrproc > 1) then
     dims(nrproc) = ((n+1) - nrcpp * (nrproc-1)) * stride
     shifts(nrproc) = shifts(nrproc-1) + dims(nrproc-1)
  else
     dims(1) = (n+1) * stride
     shifts(1) = 0
  endif
  
end subroutine



!!!!!! my_mpi
! -------------------------------------------------------------------
! Linearizes and transposes an array.
!
! F        - input rank-2 array
! F_lin    - output rank-1 array
! elems    - first dimension of F
! stride   - second dimension of F
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
integer(kind=4), intent(in) :: elems, stride
real   (kind=8), intent(in) :: F(elems, stride)
real   (kind=8), intent(out):: F_lin(elems * stride)
integer(kind=4) :: i, a, b

do i = 1,elems
  a = (i-1) * stride + 1
  b = i * stride
  F_lin(a:b) = F(i,:)
enddo

end subroutine



!!!!!! my_mpi
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
integer(kind=4), intent(in) :: elems, stride
real   (kind=8), intent(in)  :: F_lin(elems * stride)
real   (kind=8), intent(out) :: F(elems, stride)
integer(kind=4) :: i, a, b

do i = 1,elems
  a = (i-1) * stride + 1
  b = i * stride
  F(i,:) = F_lin(a:b)
enddo

end subroutine



!!!!!! my_mpi
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
implicit none
include "mpif.h"
integer(kind=4), intent(in) :: n, elems, stride, comm
real   (kind=8), intent(in) :: F(elems, stride)
integer(kind=4), intent(in) :: dims(:), shifts(:)
real   (kind=8), intent(out):: F_out(n+1, stride)
integer(kind=4), intent(out):: ierr
real   (kind=8) :: F_lin(elems * stride), F_out_lin((n+1) * stride)

call Linearize(F,F_lin,elems,stride)

call mpi_gatherv(F_lin, &
  elems * stride,       &
  MPI_DOUBLE_PRECISION, &
  F_out_lin,            &
  dims, shifts,         &
  MPI_DOUBLE_PRECISION, &
  0, comm, ierr)

call Delinearize(F_out_lin,F_out,n+1,stride)

end subroutine



!!!!!! my_mpi
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
implicit none
include "mpif.h"
integer(kind=4), intent(in) :: n, elems, stride, comm
real   (kind=8), intent(in) :: F(n+1, stride)
integer(kind=4), intent(in) :: dims(:), shifts(:)
real   (kind=8), intent(out):: F_out(elems * stride)
integer(kind=4), intent(out):: ierr
real   (kind=8) :: F_lin((n+1) * stride)

call Linearize(F, F_lin, n+1, stride)

call mpi_scatterv(F_lin, &
  dims, shifts,          &
  MPI_DOUBLE_PRECISION,  &
  F_out,                 &
  elems * stride,        &
  MPI_DOUBLE_PRECISION,  &
  0, comm, ierr)

end subroutine



!!!!!! my_mpi
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
integer(kind=4), intent(in) :: n, elems, stride, comm
real   (kind=8), intent(in) :: F(n+1, stride)
integer(kind=4), intent(in) :: dims(:), shifts(:)
real   (kind=8), intent(out):: F_out(elems, stride)
integer(kind=4), intent(out) :: ierr
real   (kind=8) :: F_out_lin(elems * stride)

call Scatter2(F,F_out_lin,n,elems,stride,dims,shifts,comm,ierr)
call Delinearize(F_out_lin, F_out, elems, stride)

end subroutine



!!!!!! my_mpi
! -------------------------------------------------------------------
! Broadcasts computed partial solution along one axis
!
! F        - data to distribute
! F_out    - buffer to receive data
! n        - problem size
! elems    - length of received slice
! stride   - total size of each slice layer
! dims     - sizes of slices for all processors
! shifts   - offsets (linearized) of slices for all processors
! comm     - communicator of the axis
! ierr     - error code output
! -------------------------------------------------------------------
subroutine AllGather(F, F_out, n, elems, stride, dims, shifts, comm)
implicit none
include "mpif.h"
integer(kind=4), intent(in) :: n, elems, stride, comm
real   (kind=8), intent(in) :: F(elems, stride)
real   (kind=8), intent(out) :: F_out(n+1, stride)
integer(kind=4) :: dims(:), shifts(:)
real   (kind=8) :: F_lin(elems*stride), F_out_lin((n+1)*stride)
integer(kind=4) :: ierr

call Linearize(F,F_lin,elems,stride)

call mpi_allgatherv(F_lin, &
  elems * stride,          &
  MPI_DOUBLE_PRECISION,    &
  F_out_lin,               &
  dims, shifts,            &
  MPI_DOUBLE_PRECISION,    &
  comm, ierr)

call Delinearize(F_out_lin,F_out,n+1,stride)

end subroutine


end module

