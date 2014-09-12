      program main
c    
      use parallelism
      use projection_engine
      use communicators
      use utils
c
      implicit none
c
      integer(kind=4) :: n, p, nelem
      real   (kind=8), allocatable, dimension(:) :: U
      real   (kind=8), allocatable, dimension(:,:) :: M
c-> DEBUG
      real (kind=8), allocatable, dimension(:,:) :: MFULL
      real (kind=8), allocatable, dimension(:) :: FFULL
c<-DEBUG
      real (kind=8), allocatable, dimension(:,:) :: F,F2,F3
      real (kind=8), allocatable, dimension(:,:) :: F_out
      real (kind=8), allocatable, dimension(:,:) :: F2_out
      real (kind=8), allocatable, dimension(:,:) :: F3_out
      real (kind=8), allocatable, dimension(:,:) :: Result
      real (kind=8), allocatable :: R(:,:,:,:)
      integer :: request(3*3*3*2), stat(MPI_STATUS_SIZE)

      integer, allocatable, dimension(:) :: dimensionsX
      integer, allocatable, dimension(:) :: dimensionsY
      integer, allocatable, dimension(:) :: dimensionsZ
      integer, allocatable, dimension(:) :: shiftsX
      integer, allocatable, dimension(:) :: shiftsY
      integer, allocatable, dimension(:) :: shiftsZ
      integer :: Fsize,Fsize_out
c-> DEBUG
      integer, allocatable, dimension(:) ::testF
      integer, allocatable, dimension(:) ::testF_out
c<-
      integer, allocatable, dimension(:) :: IPIV
      integer :: KL, KU
      integer :: nrcppx,nrcppy,nrcppz
      integer :: ibegx,iendx,ibegy,iendy,ibegz,iendz
      integer :: sx,sy,sz
      integer :: i,j,k,s,nzero,ind
      integer :: iret,ierr
      integer :: idebug,iprint,iinfo
      integer :: iclock,iclock_init
      integer :: iclock_gather1,iclock_gather2,iclock_gather3
      integer :: iclock_solve1,iclock_solve2,iclock_solve3
      integer :: iclock_scatter1,iclock_scatter2,iclock_scatter3
      integer :: iclock_i1,iclock_i2,iclock_i3,iclock_i4
      real (kind=8) ::dtime,dtime_init
      real (kind=8) ::dtime_gather1,dtime_gather2,dtime_gather3
      real (kind=8) ::dtime_scatter1,dtime_scatter2,dtime_scatter3
      real (kind=8) ::dtime_solve1,dtime_solve2,dtime_solve3
      real (kind=8) ::dtime_i1,dtime_i2,dtime_i3,dtime_i4

      integer :: iter = 0, steps = 10
      real (kind=8) :: t = 0, Dt = 1.d-3
      integer :: obegx,oendx,obegy,oendy,obegz,oendz
      integer, dimension(3) :: ibegsx,iendsx,ibegsy,iendsy,ibegsz,iendsz

      idebug=0
      iprint=0
      iinfo=1
c
      call start_clock(iclock)
      call start_clock(iclock_init)
c
      if(MYRANK.eq.0)then
      call start_clock(iclock_i1)
      endif
c
      call initialize_parameters
      call initialize_parallelism
      call create_communicators

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(MYRANK.eq.0)then
      call stop_clock(dtime_i1,iclock_i1)
      write(*,*)'create_communicators:',dtime_i1
      call start_clock(iclock_i2)
      endif
c
c      call create_grids
c
c      if(MYRANK.eq.0)iprint=1
      if(MYRANK.eq.0)iinfo=1
c
      if(iinfo.eq.1)write(*,*)PRINTRANK,'INITIALIZATION'
c
c  ...prepare the problem dimensions
      p=ORDER !order
      n=SIZE !intervals

      if(iinfo.eq.1)then
        write(*,*)'p',p,'n',n,'size of U',n+p+2
      endif
c
      if(SIZE<NRPROCX.or.SIZE<NRPROCY.or.SIZE<NRPROCZ)then
        write(*,*)'Number of elements smaller than number of processors'
        stop
      endif
c
      allocate(U(n+p+2)) !knot vector
      KL=p; KU=p
c PARALLEL - distributed
c number of columns per processors
      call ComputeEndpoints(MYRANKX,NRPROCX,n,nrcppx,ibegx,iendx)
      call ComputeEndpoints(MYRANKY,NRPROCY,n,nrcppy,ibegy,iendy)
      call ComputeEndpoints(MYRANKZ,NRPROCZ,n,nrcppz,ibegz,iendz)

      sx = iendx - ibegx + 1
      sy = iendy - ibegy + 1
      sz = iendz - ibegz + 1

      if(iinfo.eq.1)then
        write(*,*)'Number of cols per processor:',nrcppx,nrcppy,nrcppz
        write(*,*)'ibegx,iendx',ibegx,iendx
        write(*,*)'ibegy,iendy',ibegy,iendy
        write(*,*)'ibegz,iendz',ibegz,iendz
      endif
c
c ....prepare dimensions vectors
      call FillDimVector(dimensionsX,shiftsX,nrcppx,sy*sz,n,NRPROCX)
      call FillDimVector(dimensionsY,shiftsY,nrcppy,sx*sz,n,NRPROCY)
      call FillDimVector(dimensionsZ,shiftsZ,nrcppz,sx*sy,n,NRPROCZ)

c  ...check
      if(idebug.eq.1)then
      k=0
      do i=1,NRPROCX
        k=k+dimensionsX(i)
      enddo
      if(k.ne.(n+1)*sy*sz)then
        write(*,*)PRINTRANK,'problem with dimensionsX',dimensionsX
        write(*,*)PRINTRANK,'n+1',n+1
        write(*,*)PRINTRANK,'sy',sy
        write(*,*)PRINTRANK,'sz',sz
        write(*,*)PRINTRANK,'nrcppx',nrcppx
        stop
      endif
      k=0
      do i=1,NRPROCY
        k=k+dimensionsY(i)
      enddo
      if(k.ne.(n+1)*sx*sz)then
        write(*,*)PRINTRANK,'problem with dimensionsY',dimensionsY
        write(*,*)PRINTRANK,'n+1',n+1
        write(*,*)PRINTRANK,'sx',sx
        write(*,*)PRINTRANK,'sz',sz
        stop
      endif
      k=0
      do i=1,NRPROCZ
        k=k+dimensionsZ(i)
      enddo
      if(k.ne.(n+1)*sx*sy)then
        write(*,*)PRINTRANK,'problem with dimensionsZ',dimensionsZ
        write(*,*)PRINTRANK,'n+1',n+1
        write(*,*)PRINTRANK,'sx',sx
        write(*,*)PRINTRANK,'sy',sy
        stop
      endif
      endif     
c

      if(iprint.eq.1)then
        write(*,*)PRINTRANK,
     .     'MYRANKX,MYRANKY,MYRANKZ',MYRANKX,MYRANKY,MYRANKZ
        write(*,*)PRINTRANK,
     .    'NRPROCX,NRPROCY,NRPROCZ',NRPROCX,NRPROCY,NRPROCZ
        write(*,*)PRINTRANK,'n+1',n+1
        write(*,*)PRINTRANK,'nrcppx,nrcppy,nrcppz',nrcppx,nrcppy,nrcppz
        write(*,*)PRINTRANK,'ibegx,iendx',ibegx,iendx
        write(*,*)PRINTRANK,'ibegy,iendy',ibegy,iendy
        write(*,*)PRINTRANK,'ibegz,iendz',ibegz,iendz
      endif
c Now, we distribute RHS matrices

      allocate(M(2*KL+KU+1,n+1))
c OLD: MP start with system fully generated along X
c      allocate( F((n+1),(sy)*(sz))) !x,y,z
      allocate( F(sx,sy*sz)) !x,y,z
      allocate(F2(sy,sx*sz)) !y,x,z
      allocate(F3(sz,sx*sy)) !z,x,y
      allocate(Result(sz,sx*sy)) !z,x,y
      if(MYRANKX.eq.0.or.MYRANKY.eq.0.or.MYRANKZ.eq.0)then
        allocate(IPIV(n+1))
      endif
      allocate(R((nrcppz+p-2)*(nrcppx+p-2)*(nrcppy+p-2),3,3,3))

c      allocate(M(2*KL+KU+1,n+1))
c      allocate(F((n+1),(n+1)*(n+1))) !x,y,z
c      allocate(F2((n+1),(n+1)*(n+1))) !y,x,z
c      allocate(F3((n+1),(n+1)*(n+1))) !z,x,y
c-> DEBUG
      if(idebug.eq.1)then
        allocate(MFULL(n+1,n+1))
        allocate(FFULL((n+1)*(n+1)*(n+1)))
      endif

! Compute indices for neighbours
      ibegsx=-1
      iendsx=-1
      ibegsy=-1
      iendsy=-1
      ibegsz=-1
      iendsz=-1
      do i = max(MYRANKX-1,0)+1,min(MYRANKX+1,NRPROCX-1)+1
        call ComputeEndpoints(i-1,NRPROCX,n,nrcppx,
     .    ibegsx(i-MYRANKX+1),iendsx(i-MYRANKX+1))
      enddo
      do i = max(MYRANKY-1,0)+1,min(MYRANKY+1,NRPROCY-1)+1
        call ComputeEndpoints(i-1,NRPROCY,n,nrcppy,
     .    ibegsy(i-MYRANKY+1),iendsy(i-MYRANKY+1))
      enddo
      do i = max(MYRANKZ-1,0)+1,min(MYRANKZ+1,NRPROCZ-1)+1
        call ComputeEndpoints(i-1,NRPROCZ,n,nrcppz,
     .    ibegsz(i-MYRANKZ+1),iendsz(i-MYRANKZ+1))
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(MYRANK.eq.0)then
        call stop_clock(dtime_i2,iclock_i2)
        write(*,*)'allocations:',dtime_i2
        call start_clock(iclock_i3)
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      U(1:p+1)=0.d0
      U(n+2:n+p+2)=1.d0
      do i=p+2,n+1
        U(i)=real(i-p-1)/real(n-p+1)
      enddo
c
      nelem = CountSpans(n,p,U)
c      if(iprint.eq.1)then
      if(iinfo.eq.1)then
        write(*,*)'n,p,nelem',n,p,nelem
        write(*,*)'U',U
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Iterations 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do iter = 1,steps

      write(*,*)'Iteration',iter,'/',steps
      write(*,*)'t = ',t
c
c  ...generate the 1D matrix
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      M=0
      call Form1DMassMatrix(KL,KU,U,p,n,nelem,M)
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'M'
        do i=1,2*KL+KU+1
          write(*,*)PRINTRANK,M(i,1:n+1)
        enddo
      endif
      !!$ U is the knot vector
      !!$ p is the polynomial order
      !!$ n is the index of the last control point
      !!$ nelem you get from running CountSpans
      !!$ M is the dense matrix
c-> DEBUG
      if(idebug.eq.1)then
        MFULL=0
        call Form1DMassMatrixFULL(U,p,n,nelem,MFULL) 
        if(iprint.eq.1)then
          write(*,*)'MFULL='
          do i=1,n+1
            write(*,*)i,'row=',MFULL(i,1:n+1)
          enddo
          call pause
        endif
      endif
c<- DEBUG

c
c  ...generate the RHS vectors
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(MYRANK.eq.0)then
      call stop_clock(dtime_i3,iclock_i3)
      write(*,*)'Form 1D Mass Matrix:',dtime_i3
      call start_clock(iclock_i4)
      endif
      call Form3DRHS
     .   (U,p,n,nelem,nrcppx,
     .    U,p,n,nelem,nrcppy,
     .    U,p,n,nelem,nrcppz,
     .    ibegx,iendx,MYRANKX,NRPROCX,
     .    ibegy,iendy,MYRANKY,NRPROCY,
     .    ibegz,iendz,MYRANKZ,NRPROCZ,
     .    ibegsx,iendsx,ibegsy,iendsy,ibegsz,iendsz,
     .    F,R,t)
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'F'
        do i=1,sx
          write(*,*)PRINTRANK,F(i,1:sy*sz)
        enddo
      endif
c-> DEBUG
      if(idebug.eq.1)then
        call Form3DRHSFULL
     .   (U,p,n,nelem,U,p,n,nelem,U,p,n,nelem,FFULL)
        if(iprint.eq.1)then
          write(*,*)'FFULL=',FFULL
          call pause
        endif
      endif
c<- DEBUG
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(MYRANK.eq.0)then
      call stop_clock(dtime_i4,iclock_i4)
      write(*,*)'Form 3D RHS:',dtime_i4
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ...solve the first problem
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c
      if(MYRANK.eq.0)then
      call stop_clock(dtime_init,iclock_init)
      write(*,*)dtime_init
      call start_clock(iclock_gather1)
      endif
c
      if(iinfo.eq.1)write(*,*)PRINTRANK,'1a) GATHER'
c
c
c-> MP start with system fully generated along X
c-> DEBUG
c      if(idebug.eq.1)then
c        allocate(testF(1))
c        allocate(testF_out(NRPROCX))
c        testF=MYRANKX
c        testF_out=0
c        call mpi_gather(testF,
c     .   1,
c     .   MPI_INTEGER,
c     .   testF_out,NRPROCX,MPI_INTEGER,
c     .   nzero,COMMX,ierr)
c      endif
c<-
      allocate(F_out((n+1),sy*sz))
      call Gather(F,F_out,n,sx,sy*sz,dimensionsX,shiftsX,COMMX,ierr)
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'after call mpi_gather'
        write(*,*)PRINTRANK,'ierr',ierr
        write(*,*)PRINTRANK,'F_out:'
        do i=1,n+1
          write(*,*)PRINTRANK,i,'row=',F_out(i,1:sy*sz)
        enddo
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c<-

c       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
c       .. Scalar Arguments ..
c       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
c       .. Array Arguments ..
c       INTEGER            IPIV( * )
c       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )

      if(MYRANK.eq.0)then
        call stop_clock(dtime_gather1,iclock_gather1)
        write(*,*)dtime_gather1
        call start_clock(iclock_solve1)
      endif
c
      if(MYRANKX.eq.0)then
        IPIV=0
c
        if(iinfo.eq.1)write(*,*)PRINTRANK,'1b) SOLVE THE FIRST PROBLEM'
c
        if(iprint.eq.1)then
          write(*,*)'CALL DGBSV'
          write(*,*)'N=',n+1
          write(*,*)'KL=',KL
          write(*,*)'KU=',KU
          write(*,*)'NRHS=',sy*sz
          write(*,*)'AB='
          do i=1,2*KL+KU+1
            write(*,*)i,'row=',M(i,1:n+1)
          enddo
          write(*,*)'LDAB=',2*KL+KU+1
          write(*,*)'IPIV=',IPIV
          write(*,*)'B='
          do i=1,n+1
          write(*,*)i,'row=',F_out(i,1:sy*sz)
          enddo
          write(*,*)'LDB=',n+1
        endif
        call DGBSV(n+1,KL,KU,sy*sz,M,2*KL+KU+1,IPIV,F_out,n+1,iret)
        if(iprint.eq.1)then
          write(*,*)'iret=',iret
          write(*,*)'Solution='
          do i=1,n+1
            write(*,*)i,'row=',F_out(i,1:sy*sz)
          enddo
        endif
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(MYRANK.eq.0)then
        call stop_clock(dtime_solve1,iclock_solve1)
        write(*,*)dtime_solve1
        call start_clock(iclock_scatter1)
      endif

      if(iinfo.eq.1)write(*,*)PRINTRANK,'1c) SCATTER'
      allocate(F2_out(sx,sy*sz)) 
      call Scatter(F_out,F2_out,n,sx,sy*sz,dimensionsX,shiftsX,COMMX,
     .  ierr)
      deallocate(F_out)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(MYRANK.eq.0)then
        call stop_clock(dtime_scatter1,iclock_scatter1)
        write(*,*)dtime_scatter1
        call start_clock(iclock_gather2)
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-> DEBUG
      if(MYRANK.eq.0)iprint=1
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'call ReorderRHSForY'
        write(*,*)PRINTRANK,'U',U
        write(*,*)PRINTRANK,'p',p
        write(*,*)PRINTRANK,'n',n
        write(*,*)PRINTRANK,'ibegx,iendx,MYRANKX,NRPROCX',
     .    ibegx,iendx,MYRANKX,NRPROCX
        write(*,*)PRINTRANK,'ibegy,iendy,MYRANKY,NRPROCY',
     .    ibegy,iendy,MYRANKY,NRPROCY
        write(*,*)PRINTRANK,'ibegz,iendz,MYRANKZ,NRPROCZ',
     .    ibegz,iendz,MYRANKZ,NRPROCZ
        write(*,*)PRINTRANK,'F2_out',F2_out
c        write(*,*)PRINTRANK,'F2',F2
      endif
c<-

      call ReorderRHSForY
     .   (U,p,n,nelem,U,p,n,nelem,U,p,n,nelem, 
     .    ibegx,iendx,MYRANKX,NRPROCX,
     .    ibegy,iendy,MYRANKY,NRPROCY,
     .    ibegz,iendz,MYRANKZ,NRPROCZ,F2_out,F2)
      deallocate(F2_out)

      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'after ReorderRHSForY'
        write(*,*)PRINTRANK,'F2:'
        do i=1,sy
          write(*,*)PRINTRANK,i,'row=',F2(i,1:sx*sz)
        enddo
      endif
      iprint=0      
c  ...solve the second problem
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(iinfo.eq.1)write(*,*)PRINTRANK,'2a) GATHER'

      allocate(F2_out((n+1),sx*sz))
      call Gather(F2,F2_out,n,sy,sx*sz,dimensionsY,shiftsY,COMMY,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(MYRANK.eq.0)then
        call stop_clock(dtime_gather2,iclock_gather2)
        write(*,*)dtime_gather2
        call start_clock(iclock_solve2)
      endif
c
c       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
c       .. Scalar Arguments ..
c       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
c       .. Array Arguments ..
c       INTEGER            IPIV( * )
c       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )

      if(MYRANKY.eq.0)then
        IPIV=0
c
        if(iinfo.eq.1)write(*,*)PRINTRANK,'2b) SOLVE THE SECOND PROBLEM'
c
        M=0
        call Form1DMassMatrix(KL,KU,U,p,n,nelem,M)
        if(iprint.eq.1)then
          write(*,*)PRINTRANK,'M'
          do i=1,2*KL+KU+1
            write(*,*)PRINTRANK,M(i,1:n+1)
          enddo
        endif
c
        if(iprint.eq.1)then
          write(*,*)'CALL DGBSV'
          write(*,*)'N=',n+1
          write(*,*)'KL=',KL
          write(*,*)'KU=',KU
          write(*,*)'NRHS=',(n+1)*(n+1)
          write(*,*)'AB='
          do i=1,2*KL+KU+1
            write(*,*)i,'row=',M(i,1:n+1)
          enddo
          write(*,*)'LDAB=',2*KL+KU+1
          write(*,*)'IPIV=',IPIV
          write(*,*)'B='
          do i=1,n+1
          write(*,*)i,'row=',F2_out(i,1:sx*sz)
          enddo
          write(*,*)'LDB=',n+1
        endif
        call DGBSV(n+1,KL,KU,sx*sz,M,2*KL+KU+1,IPIV,F2_out,n+1,iret)
        if(iprint.eq.1)then
          write(*,*)'iret=',iret
          write(*,*)'Solution='
          do i=1,n+1
          write(*,*)i,'row=',F2_out(i,1:sx*sz)
          enddo
        endif
      endif
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(MYRANK.eq.0)then
      call stop_clock(dtime_solve2,iclock_solve2)
      write(*,*)dtime_solve2
      call start_clock(iclock_scatter2)
      endif

      if(iinfo.eq.1)write(*,*)PRINTRANK,'2c) SCATHER'

c-> CORRECTION
      allocate(F3_out(sy,sx*sz)) 
      call Scatter(F2_out,F3_out,n,sy,sx*sz,dimensionsY,shiftsY,COMMY,
     .  ierr)
      deallocate(F2_out)
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'after call mpi_scatterv'
        write(*,*)PRINTRANK,'ierr',ierr
        write(*,*)PRINTRANK,'F3_out:'
        do i=1,sy
          write(*,*)PRINTRANK,i,'row=',F3_out(i,1:sx*sz)
        enddo
      endif
c<- test
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(MYRANK.eq.0)then
      call stop_clock(dtime_scatter2,iclock_scatter2)
      write(*,*)dtime_scatter2
      call start_clock(iclock_gather3)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c-> DEBUG
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'call ReorderRHSForZ'
        write(*,*)PRINTRANK,'U',U
        write(*,*)PRINTRANK,'p',p
        write(*,*)PRINTRANK,'n',n
        write(*,*)PRINTRANK,'ibegx,iendx,MYRANKX,NRPROCX',
     .    ibegx,iendx,MYRANKX,NRPROCX
        write(*,*)PRINTRANK,'ibegy,iendy,MYRANKY,NRPROCY',
     .    ibegy,iendy,MYRANKY,NRPROCY
        write(*,*)PRINTRANK,'ibegz,iendz,MYRANKZ,NRPROCZ',
     .    ibegz,iendz,MYRANKZ,NRPROCZ
        write(*,*)PRINTRANK,'F2_out',F3_out
c        write(*,*)PRINTRANK,'F3',F3
      endif
      if(MYRANK.eq.0)iprint=1
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'before ReorderRHSForZ'
        write(*,*)PRINTRANK,'F3_out:'
        do i=1,sz
          write(*,*)PRINTRANK,i,'row=',F3_out(i,1:sx*sy)
        enddo
      endif
c<-
c  ...reorder right hand sides
      call ReorderRHSForZ
     .   (U,p,n,nelem,U,p,n,nelem,U,p,n,nelem, 
     .    ibegx,iendx,MYRANKX,NRPROCX,
     .    ibegy,iendy,MYRANKY,NRPROCY,
     .    ibegz,iendz,MYRANKZ,NRPROCZ,F3_out,F3)
      deallocate(F3_out)

      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'after ReorderRHSForZ'
        write(*,*)PRINTRANK,'F3:'
        do i=1,sz
          write(*,*)PRINTRANK,i,'row=',F3(i,1:sx*sy)
        enddo
      endif
      iprint=0
c  ...solve the third problem
c
      if(iinfo.eq.1)write(*,*)PRINTRANK,'3a) GATHER'
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      allocate(F3_out((n+1),sx*sy))
      call Gather(F3,F3_out,n,sz,sx*sy,dimensionsZ,shiftsZ,COMMZ,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(MYRANK.eq.0)then
      call stop_clock(dtime_gather3,iclock_gather3)
      write(*,*)dtime_gather3
      call start_clock(iclock_solve3)
      endif
c
c
c       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
c       .. Scalar Arguments ..
c       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
c       .. Array Arguments ..
c       INTEGER            IPIV( * )
c       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )

      if(MYRANKZ.eq.0)then
        IPIV=0
c
      if(iinfo.eq.1)write(*,*)PRINTRANK,'3b) SOLVE THE THIRD PROBLEM'
c
        M=0
        call Form1DMassMatrix(KL,KU,U,p,n,nelem,M)
        if(iprint.eq.1)then
          write(*,*)PRINTRANK,'M'
          do i=1,2*KL+KU+1
            write(*,*)PRINTRANK,M(i,1:n+1)
          enddo
        endif

        if(iprint.eq.1)then
          write(*,*)'CALL DGBSV'
          write(*,*)'N=',n+1
          write(*,*)'KL=',KL
          write(*,*)'KU=',KU
          write(*,*)'NRHS=',(n+1)*(n+1)
          write(*,*)'AB='
          do i=1,2*KL+KU+1
            write(*,*)i,'row=',M(i,1:n+1)
          enddo
          write(*,*)'LDAB=',2*KL+KU+1
          write(*,*)'IPIV=',IPIV
          write(*,*)'B='
          do i=1,n+1
          write(*,*)i,'row=',F3_out(i,1:sx*sy)
          enddo
          write(*,*)'LDB=',n+1
        endif
        call DGBSV(n+1,KL,KU,sx*sy,M,2*KL+KU+1,IPIV,F3_out,n+1,iret)
        if(iprint.eq.1)then
          write(*,*)'iret=',iret
          write(*,*)'Solution='
          do i=1,n+1
          write(*,*)i,'row=',F3_out(i,1:sx*sy)
          enddo
        endif
      endif
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c
      if(MYRANK.eq.0)then
      call stop_clock(dtime_solve3,iclock_solve3)
      write(*,*)dtime_solve3
      call start_clock(iclock_scatter3)
      endif

c
      if(iinfo.eq.1)write(*,*)PRINTRANK,'3c) SCATTER'
c
c-> CORRECTION
      call Scatter2(F3_out,Result,n,sz,sx*sy,dimensionsZ,shiftsZ,COMMZ,
     .  ierr)
      deallocate(F3_out)

      !do i=ibegx,iendx
      !  do j=ibegy,iendy
      !    do k=ibegz,iendz
      !       Result(k-ibegz+1,(j-ibegy)*sx + i-ibegx+1)=(10*i+j)*10 + k 
      !    enddo
      !  enddo
      !enddo

      s = 1
      do i = max(MYRANKX-1,0)+1,min(MYRANKX+1,NRPROCX-1)+1
        do j = max(MYRANKY-1,0)+1,min(MYRANKY+1,NRPROCY-1)+1
          do k = max(MYRANKZ-1,0)+1,min(MYRANKZ+1,NRPROCZ-1)+1
            call mpi_isend(Result,
     .        sx*sy*sz,
     .        MPI_DOUBLE_PRECISION,
     .        processors(i,j,k),
     .        0,
     .        MPI_COMM_WORLD,
     .        request(s),ierr)
            s = s + 1
            call mpi_irecv(R(:,i-MYRANKX+1,j-MYRANKY+1,k-MYRANKZ+1),
     .        nrcppx*nrcppy*nrcppz,
     .        MPI_DOUBLE_PRECISION,
     .        processors(i,j,k),
     .        0,
     .        MPI_COMM_WORLD,
     .        request(s),ierr)
            s = s + 1
          enddo
        enddo
      enddo

      do i = 1,s-1
        call mpi_wait(request(i),stat,ierr)
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if (MYRANK == 0) then
        write(*,*)PRINTRANK,'R:'
        do i = max(MYRANKX-1,0)+1,min(MYRANKX+1,NRPROCX-1)+1
        do j = max(MYRANKY-1,0)+1,min(MYRANKY+1,NRPROCY-1)+1
        do k = max(MYRANKZ-1,0)+1,min(MYRANKZ+1,NRPROCZ-1)+1
          write(*,*)'(i,j,k)=',i,j,k
          call ComputeEndpoints(i-1,NRPROCX,n,nrcppx,obegx,oendx)
          call ComputeEndpoints(j-1,NRPROCY,n,nrcppy,obegy,oendy)
          call ComputeEndpoints(k-1,NRPROCZ,n,nrcppz,obegz,oendz)
          write(*,*) reshape(
     .      R(:,i-MYRANKX+1,j-MYRANKY+1,k-MYRANKZ+1),
     .      (/ oendz-obegz+1,oendx-obegx+1,oendy-obegy+1 /))
        enddo
        enddo
        enddo
        write(*,*)'----'
      endif

      if(MYRANK.eq.0)iprint=1
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'Result:'
        do i=1,sz
          write(*,*)PRINTRANK,i,'row=',Result(i,:)
        enddo
      endif

      if(MYRANK.eq.0)then
      call stop_clock(dtime_scatter3,iclock_scatter3)
      write(*,*)dtime_scatter3
      endif


      call mpi_barrier(MPI_COMM_WORLD,ierr)

      t = t + Dt

c End of iterations
      enddo

      deallocate(shiftsX)
      deallocate(shiftsY)
      deallocate(shiftsZ)
      deallocate(dimensionsX)
      deallocate(dimensionsY)
      deallocate(dimensionsZ)

      if(allocated(IPIV))deallocate(IPIV)
      deallocate(U)
      deallocate(M)
      deallocate(F)

      if(idebug.eq.1)then
        deallocate(MFULL)
        deallocate(FFULL)
      endif

      call mpi_finalize(ierr)
      if(iinfo.eq.1)write(*,*)PRINTRANK,"Exiting..."

      if(MYRANK.eq.0)then
      call stop_clock(dtime,iclock)
      write(*,*)dtime
      endif

      end

