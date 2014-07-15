      module parallelism
c
      save

c  ...rank of this processor
      integer :: MYRANK
      integer :: NRPROC
c  ...integer coordinates of processor along X, Y and Z
      integer :: MYRANKX,MYRANKY,MYRANKZ
c  ...rank of this processor converted to a string
      character(len=6) :: PRINTRANK
c  ...total number of processors along X, Y and Z
      integer, parameter :: NRPROCX = 3
      integer, parameter :: NRPROCY = 2
      integer, parameter :: NRPROCZ = 2

c  ...total parallel execution time
      double precision :: DTIME_PARALLEL
      integer :: ITIME_PARALLEL
c
      contains
c
c     subroutine sets NRPROC, MYRANK(X,Y,Z)
c     
      subroutine initialize_parallelism
c
      include "mpif.h"c
c
      character(4) :: temp_printrank
c  ...initialize mpi
      call mpi_init(i1)
      call mpi_comm_size(MPI_COMM_WORLD,NRPROC,i2)
      call mpi_comm_rank(MPI_COMM_WORLD,MYRANK,i3)
      if ((i1+i2+i3)/=0) then
        write(*,*)MYRANK,':',
     .  'main: error initializing MPI!'
        call abort
      endif
      call decompose(MYRANK,MYRANKX,MYRANKY,MYRANKZ)
      call int2str(MYRANK, temp_printrank)
      if(MYRANK.lt.10)then
        PRINTRANK = "000"//temp_printrank
      elseif(MYRANK.lt.100)then
        PRINTRANK = "00"//temp_printrank
      elseif(MYRANK.lt.1000)then
        PRINTRANK = "0"//temp_printrank
      elseif(MYRANK.lt.10000)then
        PRINTRANK = temp_printrank
      else
        write(*,*)'more then 10000 processors'
        stop
      endif
      end subroutine initialize_parallelism
c
      subroutine decompose(myrank,myrankx,myranky,myrankz)
      integer :: myrank,myrankx,myranky,myrankz
c
      myrankz = myrank / (NRPROCX*NRPROCY)
      myranky = ( myrank - myrankz*(NRPROCX*NRPROCY) ) / nrprocx
      myrankx = myrank - myrankz*(NRPROCX*NRPROCY) 
      myrankx = myrankx - myranky*NRPROCX 

      end subroutine decompose

      end module parallelism
