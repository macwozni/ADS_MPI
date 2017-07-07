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








end module

