
module utils

implicit none

include "mpif.h"
contains


! -------------------------------------------------------------------
! Calculates the range of the direction that is assigned to processor
! with specified rank.
!
! rank    - rank of the current process
! nrproc  - # of processors for this direction
! n       - size of the problem
! nrcpp   - # of columns per processor
! ibeg    - index of first assigned slice
! iend    - index of last assigned slice
! -------------------------------------------------------------------
subroutine ComputeEndpoints(rank,nrproc,n,nrcpp,ibeg,iend)
integer :: rank, nrproc, n
integer, intent(out) :: nrcpp, ibeg, iend

  nrcpp = (n+1+1) / nrproc
  if(rank == nrproc-1)then
    ibeg = nrcpp*(nrproc-1)+1
    iend = n+1
  else
    ! ibeg = nrcpp*(nrproc-2)+1
    ! iend = nrcpp*(nrproc-1)
    ibeg = nrcpp * rank + 1
    iend = nrcpp * (rank + 1)
  endif

end subroutine


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
subroutine FillDimVector(dims,shifts,nrcpp,stride,n,nrproc)
integer, allocatable :: dims(:), shifts(:)
integer :: nrcpp, stride, n, nrproc
integer:: i

  allocate(dims(nrproc))
  allocate(shifts(nrproc))

  shifts = 0
  dims = 0

  do i = 1,nrproc-1
    dims(i) = nrcpp * stride
    if (i > 1) shifts(i) = shifts(i-1) + dims(i-1)
  enddo

  dims(nrproc) = ((n+1) - nrcpp * (nrproc-1)) * stride
  shifts(nrproc) = shifts(nrproc-1) + dims(nrproc-1)

end subroutine


! -------------------------------------------------------------------
! Linearizes and transposes an array
!
! F        - input rank-2 array
! F_lin    - output rank-1 array
! elems    - first dimension of F
! stride   - second dimension of F
! -------------------------------------------------------------------
subroutine Linearize(F,F_lin,elems,stride)
integer :: elems, stride
real (kind=8), intent(in)  :: F(elems, stride)
real (kind=8), intent(out) :: F_lin(elems * stride)
integer :: i, a, b

  do i = 1,elems
    a = (i-1) * stride + 1
    b = i * stride
    F_lin(a:b) = F(i,:)
  enddo

end subroutine


! -------------------------------------------------------------------
! Delinearizes an array
!
! F_lin    - input rank-1 array
! F        - output rank-2 array
! elems    - first dimension of F
! stride   - second dimension of F
! -------------------------------------------------------------------
subroutine Delinearize(F_lin,F,elems,stride)
integer :: elems, stride
real (kind=8), intent(in)  :: F_lin(elems * stride)
real (kind=8), intent(out) :: F(elems, stride)
integer :: i, a, b

  do i = 1,elems
    a = (i-1) * stride + 1
    b = i * stride
    F(i,:) = F_lin(a:b) 
  enddo

end subroutine


! -------------------------------------------------------------------
! Gathers data along one axis
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
subroutine Gather(F,F_out,n,elems,stride,dims,shifts,comm,ierr)
integer :: n, elems, stride, comm
real (kind=8), intent(in)  :: F(elems, stride)
real (kind=8), intent(out) :: F_out(n+1, stride)
integer :: dims(:), shifts(:)
integer, intent(out) :: ierr
real (kind=8) :: F_lin(elems * stride), F_out_lin((n+1) * stride)

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
subroutine Scatter2(F,F_out,n,elems,stride,dims,shifts,comm,ierr)
integer :: n, elems, stride, comm
real (kind=8), intent(in) :: F(n+1, stride)
real (kind=8), intent(out) :: F_out(elems * stride)
integer :: dims(:), shifts(:)
integer, intent(out) :: ierr
real (kind=8) :: F_lin((n+1) * stride)

  call Linearize(F, F_lin, n+1, stride)

  call mpi_scatterv(F_lin, &
    dims, shifts,          &
    MPI_DOUBLE_PRECISION,  &
    F_out,                 &
    elems * stride,        &
    MPI_DOUBLE_PRECISION,  &
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
subroutine Scatter(F,F_out,n,elems,stride,dims,shifts,comm,ierr)
integer :: n, elems, stride, comm
real (kind=8), intent(in) :: F(n+1, stride)
real (kind=8), intent(out) :: F_out(elems, stride)
integer :: dims(:), shifts(:)
integer, intent(out) :: ierr
real (kind=8) :: F_out_lin(elems * stride)

  call Scatter2(F,F_out_lin,n,elems,stride,dims,shifts,comm,ierr)
  call Delinearize(F_out_lin, F_out, elems, stride)

end subroutine


! -------------------------------------------------------------------
! Broadcasts computed partial solution along one axisutation)
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
subroutine AllGather(F,F_out,n,elems,stride,dims,shifts,comm)
integer :: n, elems, stride, comm
real (kind=8), intent(in) :: F(elems, stride)
real (kind=8), intent(out) :: F_out(n+1, stride)
integer :: dims(:), shifts(:)
real (kind=8) :: F_lin(elems*stride), F_out_lin((n+1)*stride)
integer :: ierr

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

