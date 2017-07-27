module parallelism

implicit none

save

! Rank of this processor
integer(kind=4) :: MYRANK

! Number of processors
integer(kind=4) :: NRPROC

! Integer coordinates of processor along X, Y and Z
integer(kind=4) :: MYRANKX, MYRANKY, MYRANKZ

! Total number of processors along X, Y and Z
integer(kind=4) :: NRPROCX
integer(kind=4) :: NRPROCY
integer(kind=4) :: NRPROCZ

! Rank of this processor converted to a string
character(len=7) :: PRINTRANK

contains

     
! -------------------------------------------------------------------
! Initializes MPI communicators and global variables of this module.
! -------------------------------------------------------------------
subroutine InitializeParallelism(procx,procy,procz,ierr)
USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT ! access computing environment
implicit none
include "mpif.h"
integer(kind=4), intent(in)  :: procx,procy,procz
integer(kind=4), intent(out) :: ierr
character(4)    :: buffer
integer(kind=4) :: i1, i2, i3

NRPROCX = procx
NRPROCY = procy
NRPROCZ = procz

! Initialize MPI
call mpi_init(i1)
call mpi_comm_size(MPI_COMM_WORLD,NRPROC,i2)
call mpi_comm_rank(MPI_COMM_WORLD,MYRANK,i3)

if ((i1+i2+i3) /= 0) then
  write(ERROR_UNIT,*)MYRANK,': main: error initializing MPI!'
  call abort
endif

call Decompose(MYRANK,MYRANKX,MYRANKY,MYRANKZ)
call int2str(MYRANK, buffer)

if (MYRANK < 10) then
  PRINTRANK = "0000"//buffer
else if (MYRANK < 100) then
  PRINTRANK = "000"//buffer
else if (MYRANK < 1000) then
  PRINTRANK = "00"//buffer
else if (MYRANK < 10000) then
  PRINTRANK = "0"//buffer
else
  PRINTRANK = buffer
  stop
endif

ierr=0

end subroutine 


! -------------------------------------------------------------------
! Based on the linear index (process rank) computes its coordinates
! in 3D cube (NRPROCX x NRPROCY x NRPROCZ).
!
! rank    - linear rank of the process
! rankx   - x coordinate
! ranky   - y coordinate
! rankz   - z coordinate
!
! Order of components (from slowest changing): Z, Y, X
!   Rank    Coords
!    0    (0, 0, 0)
!    1    (0, 0, 1)
! etc.
! -------------------------------------------------------------------
subroutine Decompose(rank,rankx,ranky,rankz)
implicit none
integer(kind=4), intent(in)  :: rank
integer(kind=4), intent(out) :: rankx, ranky, rankz

rankz = rank / (NRPROCX*NRPROCY)
ranky = (rank - rankz*(NRPROCX*NRPROCY)) / NRPROCX
rankx = rank - rankz*(NRPROCX*NRPROCY) 
rankx = rankx - ranky*NRPROCX 

end subroutine 


! -------------------------------------------------------------------
! Computes global, linear index, based on coordinates in 3D cube.
!
! rankx   - x coordinate
! ranky   - y coordinate
! rankz   - z coordinate
!
! returns: linear index based on (Z, Y, X) lexicographic order
!          (Z changes slowest)
! -------------------------------------------------------------------
function LinearIndex(rankx,ranky,rankz) result (rank)
implicit none
integer(kind=4), intent(in) :: rankx, ranky, rankz
integer(kind=4) :: rank

rank = rankz
rank = rank * NRPROCY + ranky
rank = rank * NRPROCX + rankx

end function



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
subroutine ComputeEndpoints(rank,nrproc,n,p,nrcpp,ibeg,iend,mine,maxe)
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


subroutine CleanParallelism(ierr)
implicit none
integer(kind=4), intent(out) :: ierr
ierr=0
end subroutine



end module 
