      program main
c    
      use parallelism
      use projection_engine
      use communicators
c
      implicit none
c
      include "mpif.h"
c
      integer(kind=4) :: n, p, nelem
      real   (kind=8), allocatable, dimension(:) :: U
      double precision, allocatable, dimension(:,:) :: M
c-> DEBUG
      double precision, allocatable, dimension(:,:) :: MFULL
      double precision, allocatable, dimension(:) :: FFULL
c<-DEBUG
      double precision, allocatable, dimension(:,:) :: F
      double precision, allocatable, dimension(:) :: F_in_lin
      double precision, allocatable, dimension(:) :: F_out_lin
      double precision, allocatable, dimension(:,:) :: F_out
      double precision, allocatable, dimension(:,:) :: F2
      double precision, allocatable, dimension(:) :: F2_in_lin
      double precision, allocatable, dimension(:) :: F2_out_lin
      double precision, allocatable, dimension(:,:) :: F2_out
      double precision, allocatable, dimension(:,:) :: F3
      double precision, allocatable, dimension(:) :: F3_in_lin
      double precision, allocatable, dimension(:) :: F3_out_lin
      double precision, allocatable, dimension(:) :: F3_out_lin2
      double precision, allocatable, dimension(:,:) :: F3_out
      double precision, allocatable, dimension(:,:) :: Result
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
      integer :: i,j,k,nzero,ind
      integer :: iret,ierr
      integer :: idebug,iprint,iinfo
      integer :: iclock,iclock_init
      integer :: iclock_gather1,iclock_gather2,iclock_gather3
      integer :: iclock_solve1,iclock_solve2,iclock_solve3
      integer :: iclock_scatter1,iclock_scatter2,iclock_scatter3
      integer :: iclock_i1,iclock_i2,iclock_i3,iclock_i4
      double precision ::dtime,dtime_init
      double precision ::dtime_gather1,dtime_gather2,dtime_gather3
      double precision ::dtime_scatter1,dtime_scatter2,dtime_scatter3
      double precision ::dtime_solve1,dtime_solve2,dtime_solve3
      double precision ::dtime_i1,dtime_i2,dtime_i3,dtime_i4

      integer :: iter = 0, steps = 10
      real (kind=8) :: t, Dt

      idebug=0
      iprint=0
      iinfo=0
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
      nrcppx = (n+1+1)/NRPROCX
      if(MYRANKX.eq.NRPROCX-1)then
        ibegx = nrcppx*(NRPROCX-1)+1
        iendx = n+1
      else
        ibegx = nrcppx*(NRPROCX-2)+1
        iendx = nrcppx*(NRPROCX-1)
      endif
      nrcppy = (n+1+1)/NRPROCY
      if(MYRANKY.eq.NRPROCY-1)then
        ibegy = nrcppy*(NRPROCY-1)+1
        iendy = n+1
      else
        ibegy = nrcppy*(NRPROCY-2)+1
        iendy = nrcppy*(NRPROCY-1)
      endif
      nrcppz = (n+1+1)/NRPROCZ
      if(MYRANKZ.eq.NRPROCZ-1)then
        ibegz = nrcppz*(NRPROCZ-1)+1
        iendz = n+1
      else
        ibegz = nrcppz*(NRPROCZ-2)+1
        iendz = nrcppz*(NRPROCZ-1)
      endif

      sx = iendx - ibegx + 1
      sy = iendy - ibegy + 1
      sz = iendz - ibegz + 1

      if(iinfo.eq.1)then
        write(*,*)'Number of columns per processor:',MYRANK
        write(*,*)'ibegx,iendx',ibegx,iendx
        write(*,*)'ibegy,iendy',ibegy,iendy
        write(*,*)'ibegz,iendz',ibegz,iendz
      endif
c
c ....prepare dimensions vectors
      allocate(dimensionsX(NRPROCX))
      allocate(shiftsX(NRPROCX))
      shiftsX=0;dimensionsX=0
      do i=1,NRPROCX-1
        dimensionsX(i)=nrcppx*sy*sz
        if(i.gt.1)shiftsX(i)=shiftsX(i-1)+dimensionsX(i-1)
      enddo
      dimensionsX(NRPROCX)=((n+1)-nrcppx*(NRPROCX-1))*sy*sz
      shiftsX(NRPROCX)=shiftsX(NRPROCX-1)+dimensionsX(NRPROCX-1)

      allocate(dimensionsY(NRPROCY))
      allocate(shiftsY(NRPROCY))
      shiftsY=0;dimensionsY=0
      do i=1,NRPROCY-1
        dimensionsY(i)=nrcppy*sx*sz
        if(i.gt.1)shiftsY(i)=shiftsY(i-1)+dimensionsY(i-1)
      enddo
      dimensionsY(NRPROCY)=((n+1)-nrcppy*(NRPROCY-1))*sx*sz
      shiftsY(NRPROCY)=shiftsY(NRPROCY-1)+dimensionsY(NRPROCY-1)

      allocate(dimensionsZ(NRPROCZ))
      allocate(shiftsZ(NRPROCZ))
      shiftsZ=0;dimensionsZ=0
      do i=1,NRPROCZ-1
        dimensionsZ(i)=nrcppy*sx*sy
        if(i.gt.1)shiftsZ(i)=shiftsZ(i-1)+dimensionsZ(i-1)
      enddo
      dimensionsZ(NRPROCZ)=((n+1)-nrcppz*(NRPROCZ-1))*sx*sy
      shiftsZ(NRPROCZ)=shiftsZ(NRPROCZ-1)+dimensionsZ(NRPROCZ-1)

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
c      allocate(M(2*KL+KU+1,n+1))
c      allocate(F((n+1),(n+1)*(n+1))) !x,y,z
c      allocate(F2((n+1),(n+1)*(n+1))) !y,x,z
c      allocate(F3((n+1),(n+1)*(n+1))) !z,x,y
c-> DEBUG
      if(idebug.eq.1)then
        allocate(MFULL(n+1,n+1))
        allocate(FFULL((n+1)*(n+1)*(n+1)))
      endif

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
c        if(iprint.eq.1)then
c          write(*,*)'i-p-1,n,(i-p-1)/n',i-p-1,n,real(i-p-1)/real(n)
c        endif
        U(i)=real(i-p-1)/real(n)
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
      
      F=0
      call Form3DRHS
     .   (U,p,n,nelem,U,p,n,nelem,U,p,n,nelem,
     .    ibegx,iendx,MYRANKX,NRPROCX,
     .    ibegy,iendy,MYRANKY,NRPROCY,
     .    ibegz,iendz,MYRANKZ,NRPROCZ,F)
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
      allocate(F_in_lin(sx*sy*sz))
      k=1
      do i=ibegx,iendx
        F_in_lin((k-1)*sy*sz+1:k*sy*sz) = F(k,1:sy*sz)
        k=k+1
      enddo

      nzero=0
c      if(MYRANKX.eq.0)then
        allocate(F_out_lin((n+1)*sy*sz))
        allocate(F_out((n+1),sy*sz))
c      endif
c
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'call mpi_gatherv'
        write(*,*)PRINTRANK,'F_in_lin:'
        do i=1,sx
          write(*,*)PRINTRANK,i,'row=',F_in_lin(1+(i-1)*sy*sz:i*sy*sz)
        enddo
        write(*,*)PRINTRANK,'size in:',sx*sy*sz
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'size out:dimensionsX',dimensionsX
        write(*,*)PRINTRANK,'size out:shiftsX',shiftsX
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'receiver',nzero    
        write(*,*)PRINTRANK,'COMMX',COMMX    
      endif

      call mpi_gatherv(F_in_lin,
     .   sx*sy*sz,
     .   MPI_DOUBLE_PRECISION,
     .   F_out_lin,
     .   dimensionsX,
     .   shiftsX,
     .   MPI_DOUBLE_PRECISION,
     .   nzero,COMMX,ierr)

      do i=1,n+1
        F_out(i,1:sy*sz)=F_out_lin((i-1)*sy*sz+1:i*sy*sz)
      enddo
      deallocate(F_in_lin)
      deallocate(F_out_lin)
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'after call mpi_gather'
        write(*,*)PRINTRANK,'ierr',ierr
        write(*,*)PRINTRANK,'F_out:'
        do i=1,n+1
          write(*,*)PRINTRANK,i,
     .      'row=',F_out(i,1:sy*sz)
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

      allocate(F_out_lin((n+1)*sy*sz))
      F_out_lin=0
      if(MYRANKX.eq.0)then
c      do i=1,n+1
c        F_out_lin(
c     .   (i-1)*(sy)*(sz)+1
c     .   :
c     .       i*(sy)*(sz)
c     .   ) =
c     .  F_out(i,1:(sy)*(sz))
c      enddo      
        ind=1
        do i=1,n+1
          do k=ibegz,iendz
            do j=ibegy,iendy
              F_out_lin(ind) = F_out(i,(j-ibegy+1)+(k-ibegz)*sy)
c            if(MYRANK.eq.0)
c     .        write(*,*)'i,j,k',i,j,k,'->',F_out_lin(ind)
              ind=ind+1
            enddo
          enddo
        enddo
      endif      
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'call mpi_scatterv'
        write(*,*)PRINTRANK,'F_out_lin:',F_out_lin
c        do i=1,n+1
c          write(*,*)PRINTRANK,i,
c     .     'row=',F_out_lin(1+(i-1)*(sy)*(sz):
c     .                       i*(sy)*(sz))
c        enddo
        write(*,*)PRINTRANK,'size in:dimensionsX',dimensionsX
        write(*,*)PRINTRANK,'size in:shiftsX',shiftsX
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'size out:',sx*sy*sz
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'sender',nzero    
        write(*,*)PRINTRANK,'COMMX',COMMX    
      endif
      allocate(F2_out(sx,sy*sz)) 
      allocate(F2_out_lin(sx*sy*sz)) 
      deallocate(F_out)
c<- test
      call mpi_scatterv(F_out_lin,
     .   dimensionsX,shiftsX,
     .   MPI_DOUBLE_PRECISION,
     .   F2_out_lin, !test: F2_out -> F2_out_lin
     .   sx*sy*sz,
     .   MPI_DOUBLE_PRECISION,
     .   nzero,COMMX,ierr)
c-> test
      do i=1,sx
        F2_out(i,1:sy*sz) = F2_out_lin((i-1)*sy*sz+1:i*sy*sz)
      enddo
      deallocate(F2_out_lin)
c<- test      
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'after call mpi_scatterv'
        write(*,*)PRINTRANK,'ierr',ierr
        write(*,*)PRINTRANK,'F_out:'
        do i=1,sx
          write(*,*)PRINTRANK,i,
     .      'row=',F2_out(i,1:sy*sz)
        enddo
      endif
      deallocate(F_out_lin)
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

      allocate(F2_in_lin(sx*(n+1)*sz))
      k=1
      do i=ibegy,iendy
        F2_in_lin((k-1)*sx*sz+1:k*sx*sz)=F2(k,1:sx*sz)
        k=k+1
      enddo
      !deallocate(F2)
      nzero=0
c      if(MYRANKX.eq.0)then
        allocate(F2_out_lin(sx*(n+1)*sz))
        allocate(F2_out((n+1),sx*sz))
c      endif

      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'call mpi_gatherv'
        write(*,*)PRINTRANK,'F_in_lin:'
        do i=1,sy
          write(*,*)PRINTRANK,i,'row=',F2_in_lin(1+(i-1)*sx*sz:i*sx*sz)
        enddo
        write(*,*)PRINTRANK,'size in:',sx*sy*sz
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'size out:dimensionsX',dimensionsY
        write(*,*)PRINTRANK,'size out:shiftsX',shiftsY
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'receiver',nzero    
        write(*,*)PRINTRANK,'COMMY',COMMY 
      endif

      call mpi_gatherv(F2_in_lin,
     .   sx*sy*sz,
     .   MPI_DOUBLE_PRECISION,
     .   F2_out_lin,
     .   dimensionsY,
     .   shiftsY,
     .   MPI_DOUBLE_PRECISION,
     .   nzero,COMMY,ierr)

      do i=1,n+1
        F2_out(i,1:sx*sz)=F2_out_lin((i-1)*sx*sz+1:i*sx*sz)
      enddo

      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'after call mpi_gather'
        write(*,*)PRINTRANK,'ierr',ierr
        write(*,*)PRINTRANK,'F2_out:'
        do i=1,n+1
          write(*,*)PRINTRANK,i,'row=',F2_out(i,1:sx*sz)
        enddo
      endif
      deallocate(F2_out_lin)
      deallocate(F2_in_lin)

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
        !if(.not.allocated(IPIV))allocate(IPIV(n+1))
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
      allocate(F2_out_lin((n+1)*sx*sz))
      F2_out_lin=0
      if(MYRANKY.eq.0)then
        ind=1
        do j=1,n+1
          do k=ibegz,iendz
            do i=ibegx,iendx
              F2_out_lin(ind)=F2_out(j,(i-ibegx+1)+(k-ibegz)*sx)
c            if(MYRANK.eq.0)
c     .        write(*,*)'i,j,k',i,j,k,'->',F2_out_lin(ind)
              ind=ind+1
            enddo
          enddo
        enddo
      endif      
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'call mpi_scatterv'
        write(*,*)PRINTRANK,'F2_out_lin:',F2_out_lin
        write(*,*)PRINTRANK,'size in:dimensionsY',dimensionsY
        write(*,*)PRINTRANK,'size in:shiftsY',shiftsY
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'size out:',sx*sy*sz
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'sender',nzero    
        write(*,*)PRINTRANK,'COMMY',COMMY 
      endif
c<- CORRECTION
      deallocate(F2_out)
      allocate(F3_out(sy,sx*sz)) 
      allocate(F3_out_lin(sx*sy*sz)) 
c<- test      
      call mpi_scatterv(F2_out_lin,
     .   dimensionsY,shiftsY,
     .   MPI_DOUBLE_PRECISION,
     .   F3_out_lin, !test F3_out -> F3_out_lin
     .   sx*sy*sz,
     .   MPI_DOUBLE_PRECISION,
     .   nzero,COMMY,ierr)
c-> test
      do i=1,sy
        F3_out(i,1:sx*sz) = F3_out_lin((i-1)*sx*sz+1:i*sx*sz)
      enddo
      deallocate(F2_out_lin)
      deallocate(F3_out_lin)
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

      allocate(F3_in_lin(sx*sy*sz))
      k=1
      do i=ibegz,iendz
        F3_in_lin((k-1)*sx*sy+1:k*sx*sy) = F3(k,1:sx*sy)
        k=k+1
      enddo
      nzero=0
c      if(MYRANKX.eq.0)then
        allocate(F3_out_lin2(sx*sy*(n+1)))
        allocate(F3_out((n+1),sx*sy))
c      endif
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'call mpi_gatherv'
        write(*,*)PRINTRANK,'F3_in_lin:'
        do i=1,sz
          write(*,*)PRINTRANK,i,
     .     'row=',F3_in_lin(1+(i-1)*sx*sy:i*sx*sy)
        enddo
        write(*,*)PRINTRANK,'size in:',sx*sy*sz
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'size out:dimensionsX',dimensionsZ
        write(*,*)PRINTRANK,'size out:shiftsX',shiftsZ
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'receiver',nzero    
        write(*,*)PRINTRANK,'COMMZ',COMMZ
      endif

      call mpi_gatherv(F3_in_lin,
     .   sx*sy*sz,
     .   MPI_DOUBLE_PRECISION,
     .   F3_out_lin2,
     .   dimensionsZ,
     .   shiftsZ,
     .   MPI_DOUBLE_PRECISION,
     .   nzero,COMMZ,ierr)

      do i=1,n+1
        F3_out(i,1:sx*sy)=F3_out_lin2((i-1)*sx*sy+1:i*sx*sy)
      enddo
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'after call mpi_gatherv'
        write(*,*)PRINTRANK,'ierr',ierr
        write(*,*)PRINTRANK,'F3_out:'
        do i=1,n+1
          write(*,*)PRINTRANK,i,'row=',F3_out(i,1:sx*sy)
        enddo
      endif
      deallocate(F3_out_lin2)
      deallocate(F3_in_lin)

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
      allocate(F3_out_lin((n+1)*sx*sy))
      F3_out_lin=0
      if(MYRANKZ.eq.0)then
        ind=1
        do i=1,n+1
          do k=ibegy,iendy
            do j=ibegx,iendx
              F3_out_lin(ind) = F3_out(i,(j-ibegx+1)+(k-ibegy)*sx)
c            if(MYRANK.eq.0)
c     .        write(*,*)'i,j,k',i,j,k,'->',F2_out_lin(ind)
              ind=ind+1
            enddo
          enddo
        enddo
      endif      
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'call mpi_scatterv'
        write(*,*)PRINTRANK,'F3_out_lin:',F3_out_lin
        write(*,*)PRINTRANK,'size in:dimensionsZ',dimensionsZ
        write(*,*)PRINTRANK,'size in:shiftsZ',shiftsZ
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'size out:',sx*sy*sz
        write(*,*)PRINTRANK,'MPI_DOUBLE_PRECISION'
        write(*,*)PRINTRANK,'sender',nzero    
        write(*,*)PRINTRANK,'COMMZ',COMMZ 
      endif
c<- CORRECTION
c
      call mpi_scatterv(F3_out_lin,
     .   dimensionsZ,shiftsZ,
     .   MPI_DOUBLE_PRECISION,
     .   Result,
     .   sx*sy*sz,
     .   MPI_DOUBLE_PRECISION,
     .   nzero,COMMZ,ierr)

      deallocate(F3_out_lin)
      deallocate(F3_out)

      call mpi_barrier(MPI_COMM_WORLD,ierr)
c
      if(MYRANK.eq.0)then
      call stop_clock(dtime_scatter3,iclock_scatter3)
      write(*,*)dtime_scatter3
      endif
c
      if(MYRANK.eq.0)iprint=1
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'Result:'
        do i=1,sz
          write(*,*)PRINTRANK,i,'row=',Result(i,1:sx*sy)
        enddo
      endif
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      enddo
      deallocate(shiftsX)
      deallocate(shiftsY)
      deallocate(shiftsZ)
      deallocate(dimensionsX)
      deallocate(dimensionsY)
      deallocate(dimensionsZ)
c
      if(allocated(IPIV))deallocate(IPIV)
      deallocate(U)
      deallocate(M)
      deallocate(F)
c-> DEBUG
      if(idebug.eq.1)then
        deallocate(MFULL)
        deallocate(FFULL)
      endif
c<- DEBUG
c

      call mpi_finalize(ierr)
      if(iinfo.eq.1)write(*,*)PRINTRANK,"Exiting..."
c
      if(MYRANK.eq.0)then
      call stop_clock(dtime,iclock)
      write(*,*)dtime
      endif
c
      end




