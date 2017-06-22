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
subroutine InitializeParallelism
USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT ! access computing environment
implicit none
include "mpif.h"
character(4)    :: buffer
integer(kind=4) :: i1, i2, i3

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
function LinearIndex(rankx, ranky, rankz) result (rank)
implicit none
integer(kind=4), intent(in) :: rankx, ranky, rankz
integer(kind=4) :: rank

rank = rankz
rank = rank * NRPROCY + ranky
rank = rank * NRPROCX + rankx

end function


end module 
