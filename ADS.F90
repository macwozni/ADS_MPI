module ADSS

type ADS_setup
   ! Number of functions in each dimension minus one
   integer(kind=4), dimension(3) :: n

   ! Degree of approximation
   integer(kind=4), dimension(3) :: p

   ! Knot vector
   real (kind=8), allocatable, dimension(:) :: Ux
   real (kind=8), allocatable, dimension(:) :: Uy
   real (kind=8), allocatable, dimension(:) :: Uz

   ! Mass matrix
   real (kind=8), allocatable, dimension(:,:) :: Mx
   real (kind=8), allocatable, dimension(:,:) :: My
   real (kind=8), allocatable, dimension(:,:) :: Mz

   real (kind=8), allocatable, dimension(:,:) :: F, F2, F3
   real (kind=8), allocatable, dimension(:,:) :: F_out, F2_out, F3_out

   ! Buffer for coefficients of solution corresponding to neighbouring
   ! parts of the domain. It is (Nx*Ny*Nz) x 3 x 3 x 3 array, where
   ! Nx*Ny*Nz is the size of part of solution for one fragment of domain.
   real (kind=8), allocatable :: R(:,:,:,:)

   ! Number of subintervals (currently n - p + 1)
   integer(kind=4) :: nelemx
   integer(kind=4) :: nelemy
   integer(kind=4) :: nelemz

   ! Size of slices of domain in each dimension
   integer(kind=4), allocatable, dimension(:) :: dimensionsX
   integer(kind=4), allocatable, dimension(:) :: dimensionsY
   integer(kind=4), allocatable, dimension(:) :: dimensionsZ

   ! Offsets of slices of domain in each direction
   integer(kind=4), allocatable, dimension(:) :: shiftsX
   integer(kind=4), allocatable, dimension(:) :: shiftsY
   integer(kind=4), allocatable, dimension(:) :: shiftsZ

   ! Pivot array for these processes that need to solve systems
   integer(kind=4), allocatable, dimension(:) :: IPIVx
   integer(kind=4), allocatable, dimension(:) :: IPIVy
   integer(kind=4), allocatable, dimension(:) :: IPIVz

   ! Number of lower and upper diagonal entries in mass matrix
   integer(kind=4) :: KLx, KUx
   integer(kind=4) :: KLy, KUy
   integer(kind=4) :: KLz, KUz

   ! Number of columns (average) per processor
   integer(kind=4) :: nrcppx,nrcppy,nrcppz

   ! Range of piece of domain assigned to this process
   integer(kind=4) :: ibegx,iendx
   integer(kind=4) :: ibegy,iendy
   integer(kind=4) :: ibegz,iendz

   ! Size of piece of domain assigned to this process
   integer(kind=4) :: sx,sy,sz

   ! Ranges of pieces of domain around the one assigned to this process
   integer(kind=4), dimension(3) :: ibegsx,iendsx
   integer(kind=4), dimension(3) :: ibegsy,iendsy
   integer(kind=4), dimension(3) :: ibegsz,iendsz

   ! Range of elements associated with basis functions assigned to this process
   integer(kind=4) :: minex, maxex
   integer(kind=4) :: miney, maxey
   integer(kind=4) :: minez, maxez
end type
   
contains




! -------------------------------------------------------------------
! Initialization of clocks and MPI
! -------------------------------------------------------------------
subroutine initialize (nx,ny,nz,px,py,pz,ads)
use parallelism, ONLY : NRPROCX,NRPROCY,NRPROCZ
use parallelism, ONLY : PRINTRANK
use knot_vector, ONLY : PrepareKnot
implicit none
include "mpif.h"
integer(kind=4), intent(in) :: nx,ny,nz
integer(kind=4), intent(in) :: px,py,pz
type(ADS_setup), intent(out) :: ads
integer(kind=4) :: ierr

  ads%p(1) = px ! order
  ads%p(2) = py ! order
  ads%p(3) = pz ! order
  ads%n(1) = nx  ! intervals
  ads%n(2) = ny  ! intervals
  ads%n(3) = nz  ! intervals

  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'INITIALIZATION'
    write(*,*)'px',px,'py',py,'pz',pz, &
    'nx',nx,'ny',ny,'nz',nz, &
    'size of Ux',nx+px+2,'size of Uy',ny+py+2,'size of Uz',nz+pz+2
#endif

  if (nx<NRPROCX .or. ny<NRPROCY .or. nz<NRPROCZ) then
    write(*,*)'Number of elements smaller than number of processors'
    stop
  endif
  
  ads%KLx = px
  ads%KUx = px
  ads%KLy = py
  ads%KUy = py
  ads%KLz = pz
  ads%KUz = pz

  call ComputeDecomposition( &
      ads%n, &
      ads%p, &
      ads%sx,ads%sy,ads%sz, &
      ads%ibegx,ads%ibegy,ads%ibegz, &
      ads%iendx,ads%iendy,ads%iendz, &
      ads%ibegsx,ads%ibegsy,ads%ibegsz, &
      ads%iendsx,ads%iendsy,ads%iendsz, &
      ads%minex,ads%miney,ads%minez, &
      ads%maxex,ads%maxey,ads%maxez, &
      ads%nrcppx,ads%nrcppy,ads%nrcppz, &
      ads%dimensionsX,ads%dimensionsY,ads%dimensionsZ, &
      ads%shiftsX,ads%shiftsY,ads%shiftsZ)
  
#ifdef IDEBUG
    call ValidateDimensions(&
      ads%n, &
      ads%sx,ads%sy,ads%sz, &
      ads%nrcppx,ads%nrcppy,ads%nrcppz, &
      ads%dimensionsX,ads%dimensionsY,ads%dimensionsZ)
#endif

#ifdef IPRINT
    call PrintDecompositionInfo(&
      ads%n, &
      ads%nrcppx,ads%nrcppy,ads%nrcppz, &
      ads%ibegx,ads%ibegy,ads%ibegz, &
      ads%iendx,ads%iendy,ads%iendz)
#endif
  
  call AllocateArrays(&
      nx,ny,nz, &
      ads%sx,ads%sy,ads%sz, &
      ads%nrcppx,ads%nrcppy,ads%nrcppz, &
      ads%Klx,ads%Kly,ads%Klz, &
      ads%KUx,ads%KUy,ads%KUz, &
      ads%Mx,ads%My,ads%Mz, &
      ads%F,ads%F2,ads%F3, &
      ads%IPIVx,ads%IPIVy,ads%IPIVz, &
      ads%R)
  
  call PrepareKnot(ads%Ux,nx,px,ads%nelemx)
  call PrepareKnot(ads%Uy,ny,py,ads%nelemy)
  call PrepareKnot(ads%Uz,nz,pz,ads%nelemz)
  
end subroutine


! -------------------------------------------------------------------
! Establishes decomposition of the domain. Calculates size and location
! of the piece for current process.
! -------------------------------------------------------------------
subroutine ComputeDecomposition( &
   n, &
   p, &
   sx,sy,sz, &
   ibegx,ibegy,ibegz, &
   iendx,iendy,iendz, &
   ibegsx,ibegsy,ibegsz, &
   iendsx,iendsy,iendsz, &
   minex,miney,minez, &
   maxex,maxey,maxez, &
   nrcppx,nrcppy,nrcppz, &
   dimensionsX,dimensionsY,dimensionsZ, &
   shiftsX,shiftsY,shiftsZ)
use parallelism, ONLY : MYRANKX,MYRANKY,MYRANKZ,PRINTRANK, &
   NRPROCX,NRPROCY,NRPROCZ,ComputeEndpoints,FillDimVector
implicit none
integer(kind=4), intent(in), dimension(3) :: n
integer(kind=4), intent(in), dimension(3) :: p
integer(kind=4), intent(out) :: ibegx,iendx
integer(kind=4), intent(out) :: ibegy,iendy
integer(kind=4), intent(out) :: ibegz,iendz
integer(kind=4), intent(out) :: minex, maxex
integer(kind=4), intent(out) :: miney, maxey
integer(kind=4), intent(out) :: minez, maxez
integer(kind=4), intent(out) :: nrcppx,nrcppy,nrcppz
integer(kind=4), intent(out) :: sx,sy,sz
integer(kind=4), intent(out), allocatable, dimension(:) :: dimensionsX
integer(kind=4), intent(out), allocatable, dimension(:) :: dimensionsY
integer(kind=4), intent(out), allocatable, dimension(:) :: dimensionsZ
integer(kind=4), intent(out), allocatable, dimension(:) :: shiftsX
integer(kind=4), intent(out), allocatable, dimension(:) :: shiftsY
integer(kind=4), intent(out), allocatable, dimension(:) :: shiftsZ
integer(kind=4), intent(out), dimension(3) :: ibegsx,iendsx
integer(kind=4), intent(out), dimension(3) :: ibegsy,iendsy
integer(kind=4), intent(out), dimension(3) :: ibegsz,iendsz
   
integer(kind=4) :: i
integer(kind=4) :: ix, iy, iz
integer(kind=4) :: mine, maxe

  ! number of columns per processors
  call ComputeEndpoints(MYRANKX, NRPROCX, n(1), p(1), nrcppx, ibegx, iendx, minex, maxex)
  call ComputeEndpoints(MYRANKY, NRPROCY, n(2), p(2), nrcppy, ibegy, iendy, miney, maxey)
  call ComputeEndpoints(MYRANKZ, NRPROCZ, n(3), p(3), nrcppz, ibegz, iendz, minez, maxez)

  sx = iendx - ibegx + 1
  sy = iendy - ibegy + 1
  sz = iendz - ibegz + 1

#ifdef IINFO
    write(*,*)PRINTRANK,'Number of cols per processor:',nrcppx,nrcppy,nrcppz
    write(*,*)PRINTRANK,'ibegx,iendx',ibegx,iendx
    write(*,*)PRINTRANK,'ibegy,iendy',ibegy,iendy
    write(*,*)PRINTRANK,'ibegz,iendz',ibegz,iendz
#endif

  ! prepare dimensions vectors
  call FillDimVector(dimensionsX, shiftsX, nrcppx, sy*sz, n(1), NRPROCX)
  call FillDimVector(dimensionsY, shiftsY, nrcppy, sx*sz, n(2), NRPROCY)
  call FillDimVector(dimensionsZ, shiftsZ, nrcppz, sx*sy, n(3), NRPROCZ)

  ! Compute indices for neighbours
  ibegsx = -1
  iendsx = -1
  ibegsy = -1
  iendsy = -1
  ibegsz = -1
  iendsz = -1

  do i = max(MYRANKX-1,0)+1, min(MYRANKX+1,NRPROCX-1)+1
    ix = i-MYRANKX+1
    call ComputeEndpoints(i-1, NRPROCX, n(1), p(1), nrcppx, ibegsx(ix), iendsx(ix), mine, maxe)
  enddo
  do i = max(MYRANKY-1,0)+1, min(MYRANKY+1,NRPROCY-1)+1
    iy = i-MYRANKY+1
    call ComputeEndpoints(i-1,NRPROCY, n(2), p(2), nrcppy, ibegsy(iy), iendsy(iy), mine, maxe)
  enddo
  do i = max(MYRANKZ-1,0)+1, min(MYRANKZ+1,NRPROCZ-1)+1
    iz = i-MYRANKZ+1
    call ComputeEndpoints(i-1,NRPROCZ, n(3), p(3), nrcppz, ibegsz(iz), iendsz(iz), mine, maxe)
  enddo

end subroutine


! -------------------------------------------------------------------
! Allocates most of the 'static' arrays
! -------------------------------------------------------------------
subroutine AllocateArrays(&
      nx,ny,nz, &
      sx,sy,sz, &
      nrcppx,nrcppy,nrcppz, &
      Klx,Kly,Klz, &
      KUx,KUy,KUz, &
      Mx,My,Mz, &
      F,F2,F3, &
      IPIVx,IPIVy,IPIVz, &
      R)
use parallelism, ONLY : MYRANKX,MYRANKY,MYRANKZ
implicit none
include "mpif.h"
integer(kind=4), intent(in) :: nx,ny,nz
integer(kind=4), intent(in) :: sx,sy,sz
integer(kind=4), intent(in) :: nrcppx,nrcppy,nrcppz
integer(kind=4), intent(in) :: KLx, KUx
integer(kind=4), intent(in) :: KLy, KUy
integer(kind=4), intent(in) :: KLz, KUz
real (kind=8), intent(out), allocatable, dimension(:,:) :: Mx
real (kind=8), intent(out), allocatable, dimension(:,:) :: My
real (kind=8), intent(out), allocatable, dimension(:,:) :: Mz
real (kind=8), intent(out), allocatable, dimension(:,:) :: F, F2, F3
integer(kind=4), intent(out), allocatable, dimension(:) :: IPIVx
integer(kind=4), intent(out), allocatable, dimension(:) :: IPIVy
integer(kind=4), intent(out), allocatable, dimension(:) :: IPIVz
real (kind=8), intent(out), allocatable :: R(:,:,:,:)
integer :: ierr

  allocate(Mx(2*KLx+KUx+1,nx+1))
  allocate(My(2*KLy+KUy+1,ny+1))
  allocate(Mz(2*KLz+KUz+1,nz+1))

  ! OLD: MP start with system fully generated along X
  ! allocate( F((n+1),(sy)*(sz))) !x,y,z
  allocate( F(sx,sy*sz)) !x,y,z
  allocate(F2(sy,sx*sz)) !y,x,z
  allocate(F3(sz,sx*sy)) !z,x,y


  ! Processes on the border need pivot vector for LAPACK call
  if (MYRANKX == 0 .or. MYRANKY == 0 .or. MYRANKZ == 0) then
    allocate(IPIVx(nx+1))
    allocate(IPIVy(ny+1))
    allocate(IPIVz(nz+1))
  endif

  allocate(R(nrcppz*nrcppx*nrcppy,3,3,3))
  R=0.d0

  call mpi_barrier(MPI_COMM_WORLD,ierr)

end subroutine



!!!!! przeniesc do debug
! -------------------------------------------------------------------
! Prints debugging information about results of distributing
! data to neighbouring processes.
! -------------------------------------------------------------------
subroutine PrintDistributedData(ads)
use parallelism, ONLY : MYRANKX,MYRANKY,MYRANKZ,PRINTRANK, &
   NRPROCX,NRPROCY,NRPROCZ,ComputeEndpoints
implicit none
type   (ADS_setup) :: ads
integer(kind=4) :: i, j, k
integer(kind=4) :: obegx,oendx,obegy,oendy,obegz,oendz
integer(kind=4) :: mine, maxe

  write(*,*)PRINTRANK,'R:'

  do i = max(MYRANKX-1,0)+1, min(MYRANKX+1,NRPROCX-1)+1
    do j = max(MYRANKY-1,0)+1, min(MYRANKY+1,NRPROCY-1)+1
      do k = max(MYRANKZ-1,0)+1, min(MYRANKZ+1,NRPROCZ-1)+1
        write(*,*)'(i,j,k)=',i+1,j+1,k+1

        call ComputeEndpoints(i-1,NRPROCX,ads%n(1),ads%p(1),ads%nrcppx,obegx,oendx,mine,maxe)
        call ComputeEndpoints(j-1,NRPROCY,ads%n(2),ads%p(2),ads%nrcppy,obegy,oendy,mine,maxe)
        call ComputeEndpoints(k-1,NRPROCZ,ads%n(3),ads%p(3),ads%nrcppz,obegz,oendz,mine,maxe)

        write(*,*) reshape(                                 &
          ads%R(:,i-MYRANKX+1,j-MYRANKY+1,k-MYRANKZ+1),         &
          (/ oendz-obegz+1,oendx-obegx+1,oendy-obegy+1 /))
      enddo
    enddo
  enddo

  write(*,*)'----'

end subroutine


!!!! przeniesc do solver
! -------------------------------------------------------------------
! Solves 1D linear system (one of 3 steps of solving the whole),
! using DGBSV.
!
! RHS   - vector of right-hand sides, of dimension (n+1) x eqnum
! eqnum - number of right-hand sides (equations)
! -------------------------------------------------------------------
subroutine SolveOneDirection(RHS, eqnum,n,KL,KU,p,M,IPIV)
implicit none
type   (ADS_setup) :: ads
real   (kind=8) :: RHS(:,:)
integer :: KL, KU
integer, dimension(:) :: IPIV
real (kind=8), dimension(:,:) :: M
integer :: n,p
integer(kind=4) :: eqnum
integer(kind=4) :: i, iret

  IPIV = 0

#ifdef IPRINT
    write(*,*)'CALL DGBSV'
    write(*,*)'N=',n+1
    write(*,*)'KL=',KL
    write(*,*)'KU=',KU
    write(*,*)'NRHS=',eqnum
    write(*,*)'AB='
    do i = 1,2*KL+KU+1
      write(*,*)i,'row=',M(i,1:n+1)
    enddo
    write(*,*)'LDAB=',2*KL+KU+1
    write(*,*)'IPIV=',IPIV
    write(*,*)'B='
    do i = 1,n+1
      write(*,*)i,'row=',RHS(i,1:eqnum)
    enddo
    write(*,*)'LDB=',n+1
#endif

  ! SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
  ! .. Scalar Arguments ..
  ! INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
  ! .. Array Arguments ..
  ! INTEGER            IPIV( * )
  ! DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )

  call DGBSV(n+1,KL,KU,eqnum,M,2*KL+KU+1,IPIV,RHS,n+1,iret)

#ifdef IPRINT
    write(*,*)'iret=',iret
    write(*,*)'Solution='
    do i = 1,n+1
      write(*,*)i,'row=',RHS(i,1:eqnum)
    enddo
#endif

end subroutine



!!!! podzielic na wraper i czesc wlasciwa
! przeniesc czesc do solver
! -------------------------------------------------------------------
! Performs one step of the simulation
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
subroutine Step(iter,RHS_fun,ads)
use parallelism, ONLY :PRINTRANK,MYRANKX,MYRANKY,MYRANKZ
use communicators, ONLY : COMMX,COMMY,COMMZ
use reorderRHS, ONLY : ReorderRHSForX,ReorderRHSForY,ReorderRHSForZ
use projection_engine, ONLY : Form3DRHS, ComputeMassMatrix
use my_mpi, ONLY : DistributeSpline,Gather,Scatter
implicit none
include "mpif.h"
interface
  subroutine RHS_fun( &
         Xx,Xy,Xz, &
         kx,ky,kz, &
         ex,ey,ez, &
         px,py,pz, &
         nelemx,nelemy,nelemz, &
         ax,ay,az, &
         bx,by,bz, &
         NNx,NNy,NNz, &
         dux,duy,duz, &
         ibegx,ibegy,ibegz, &
         iendx,iendy,iendz, &
         minex,miney,minez, &
         maxex,maxey,maxez, &
         Uval,J,W,F)
      implicit none
      integer(kind=4), intent(in)  :: px,py,pz
      real   (kind=8), intent(in)  :: Xx(px+1,nelemx)
      real   (kind=8), intent(in)  :: Xy(py+1,nelemy)
      real   (kind=8), intent(in)  :: Xz(pz+1,nelemz)
      integer(kind=4), intent(in)  :: kx,ky,kz
      integer(kind=4), intent(in)  :: ex,ey,ez
      integer(kind=4), intent(in)  :: nelemx,nelemy,nelemz
      real   (kind=8), intent(in)  :: Uval
      integer(kind=4), intent(in)  :: ibegx,ibegy,ibegz
      integer(kind=4), intent(in)  :: iendx,iendy,iendz
      integer(kind=4), intent(in)  :: maxex,maxey,maxez
      integer(kind=4), intent(in)  :: minex,miney,minez
      integer(kind=4), intent(in)  :: ax,ay,az
      integer(kind=4), intent(in)  :: bx,by,bz
      real   (kind=8), intent(in)  :: dux,duy,duz
      real   (kind=8), intent(in)  :: J,W
      real   (kind=8), intent(in)  :: NNx(0:px-1,0:px,px+1,nelemx), &
                         NNy(0:py-1,0:py,py+1,nelemy), &
                         NNz(0:pz-1,0:pz,pz+1,nelemz)
      real   (kind=8), intent(out) :: F
  end subroutine
end interface
type   (ADS_setup) :: ads
integer(kind=4) :: iter
integer(kind=4) :: i
integer(kind=4) :: iret, ierr

  ! generate the RHS vectors
    call Form3DRHS                                    &
      (ads%Ux,ads%p(1),ads%n(1),ads%nelemx,ads%nrcppx,  &
       ads%Uy,ads%p(2),ads%n(2),ads%nelemy,ads%nrcppy,  &
       ads%Uz,ads%p(3),ads%n(3),ads%nelemz,ads%nrcppz,  &
       ads%ibegx,ads%iendx,                             &
       ads%ibegy,ads%iendy,                             &
       ads%ibegz,ads%iendz,                             &
       ads%ibegsx,ads%iendsx,ads%ibegsy,     &
       ads%iendsy,ads%ibegsz,ads%iendsz,   &
       ads%minex,ads%maxex,ads%miney,              &
       ads%maxey,ads%minez,ads%maxez,              &
       ads%R,ads%F,RHS_fun)

#ifdef IPRINT
    write(*,*)PRINTRANK,'F'
    do i = 1,ads%sx
      write(*,*)PRINTRANK,ads%F(i,1:ads%sy*ads%sz)
    enddo
#endif

  !--------------------------------------------------------------------
  ! Solve the first problem
  !--------------------------------------------------------------------
  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'1a) GATHER'
#endif

  allocate(ads%F_out((ads%n(1)+1),ads%sy*ads%sz))

  call Gather(ads%F,ads%F_out,ads%n(1),ads%sx,ads%sy*ads%sz,ads%dimensionsX,ads%shiftsX,COMMX,ierr)

#ifdef IPRINT
    write(*,*)PRINTRANK,'after call mpi_gather'
    write(*,*)PRINTRANK,'ierr',ierr
    write(*,*)PRINTRANK,'F_out:'
    do i=1,ads%n(1)+1
      write(*,*)PRINTRANK,i,'row=',ads%F_out(i,1:ads%sy*ads%sz)
    enddo
#endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (MYRANKX == 0) then
#ifdef IINFO
     write(*,*)PRINTRANK,'1b) SOLVE THE FIRST PROBLEM'
#endif

    call ComputeMassMatrix(ads%KLx,ads%KUx,ads%Ux,ads%p(1),ads%n(1),ads%nelemx,ads%Mx)
    call SolveOneDirection(ads%F_out, ads%sy*ads%sz,ads%n(1),ads%KLx,ads%KUx,ads%p(1),ads%Mx,ads%IPIVx)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'1c) SCATTER'
#endif
  allocate(ads%F2_out(ads%sx,ads%sy*ads%sz))
  call Scatter(ads%F_out,ads%F2_out,ads%n(1),ads%sx,ads%sy*ads%sz,ads%dimensionsX,ads%shiftsX,COMMX,ierr)
  deallocate(ads%F_out)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'1d) REORDER'
#endif
  call ReorderRHSForY(ads%ibegx,ads%iendx,ads%ibegy,ads%iendy,ads%ibegz,ads%iendz,ads%F2_out,ads%F2)
  deallocate(ads%F2_out)

#ifdef IPRINT
    write(*,*)PRINTRANK,'after ReorderRHSForY'
    write(*,*)PRINTRANK,'F2:'
    do i = 1,ads%sy
      write(*,*)PRINTRANK,i,'row=',ads%F2(i,1:ads%sx*ads%sz)
    enddo
#endif

  !--------------------------------------------------------------------
  ! Solve the second problem
  !--------------------------------------------------------------------
  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'2a) GATHER'
#endif

  allocate(ads%F2_out((ads%n(2)+1),ads%sx*ads%sz))
  call Gather(ads%F2,ads%F2_out,ads%n(2),ads%sy,ads%sx*ads%sz,ads%dimensionsY,ads%shiftsY,COMMY,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (MYRANKY == 0) then
#ifdef IINFO
     write(*,*)PRINTRANK,'2b) SOLVE THE SECOND PROBLEM'
#endif

    call ComputeMassMatrix(ads%KLy,ads%KUy,ads%Uy,ads%p(2),ads%n(2),ads%nelemy,ads%My)
    call SolveOneDirection(ads%F2_out, ads%sx*ads%sz,ads%n(2),ads%KLy,ads%KUy,ads%p(2),ads%My,ads%IPIVy)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'2c) SCATHER'
#endif

  ! CORRECTION
  allocate(ads%F3_out(ads%sy,ads%sx*ads%sz))
  call Scatter(ads%F2_out,ads%F3_out,ads%n(2),ads%sy,ads%sx*ads%sz,ads%dimensionsY,ads%shiftsY,COMMY,ierr)
  deallocate(ads%F2_out)

#ifdef IPRINT
    write(*,*)PRINTRANK,'after call mpi_scatterv'
    write(*,*)PRINTRANK,'ierr',ierr
    write(*,*)PRINTRANK,'F3_out:'
    do i = 1,ads%sy
      write(*,*)PRINTRANK,i,'row=',ads%F3_out(i,1:ads%sx*ads%sz)
    enddo
#endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'2d) REORDER'
#endif
  ! Reorder right hand sides
  call ReorderRHSForZ(ads%ibegx,ads%iendx,ads%ibegy,ads%iendy,ads%ibegz,ads%iendz,ads%F3_out,ads%F3)
  deallocate(ads%F3_out)

#ifdef IPRINT
    write(*,*)PRINTRANK,'after ReorderRHSForZ'
    write(*,*)PRINTRANK,'F3:'
    do i = 1,ads%sz
      write(*,*)PRINTRANK,i,'row=',ads%F3(i,1:ads%sx*ads%sy)
    enddo
#endif

  !--------------------------------------------------------------------
  ! Solve the third problem
  !--------------------------------------------------------------------
#ifdef IINFO
    write(*,*)PRINTRANK,'3a) GATHER'
#endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  allocate(ads%F3_out((ads%n(3)+1),ads%sx*ads%sy))
  call Gather(ads%F3,ads%F3_out,ads%n(3),ads%sz,ads%sx*ads%sy,ads%dimensionsZ,ads%shiftsZ,COMMZ,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (MYRANKZ == 0) then
#ifdef IINFO
     write(*,*)PRINTRANK,'3b) SOLVE THE THIRD PROBLEM'
#endif

    call ComputeMassMatrix(ads%KLz,ads%KUz,ads%Uz,ads%p(3),ads%n(3),ads%nelemz,ads%Mz)
    call SolveOneDirection(ads%F3_out, ads%sx*ads%sy,ads%n(3),ads%KLz,ads%KUz,ads%p(3),ads%Mz,ads%IPIVz)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'3c) SCATTER'
#endif

  ! CORRECTION
  allocate(ads%F_out(ads%sz,ads%sx*ads%sy))
  call Scatter(ads%F3_out,ads%F_out,ads%n(3),ads%sz,ads%sx*ads%sy,ads%dimensionsZ,ads%shiftsZ,COMMZ,ierr)
  deallocate(ads%F3_out)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'3d) REORDER'
#endif
  ! Reorder right hand sides
  call ReorderRHSForX(ads%ibegx,ads%iendx,ads%ibegy,ads%iendy,ads%ibegz,ads%iendz,ads%F_out,ads%F)
  deallocate(ads%F_out)

#ifdef IPRINT
    write(*,*)PRINTRANK,'after ReorderRHSForX'
    write(*,*)PRINTRANK,'F:'
    do i = 1,ads%sx
      write(*,*)PRINTRANK,i,'row=',ads%F(i,1:ads%sy*ads%sz)
    enddo
#endif

#ifdef IINFO
    write(*,*)PRINTRANK,'3e) DISTRIBUTE SOLUTION'
#endif
  ads%R(1:ads%sx*ads%sy*ads%sz,2,2,2) = reshape(ads%F, [ads%sx*ads%sy*ads%sz])
  call DistributeSpline(ads%R,ads%nrcppz,ads%nrcppx,ads%nrcppy,ads%R)

#ifdef IPRINT
    write(*,*)PRINTRANK,'Result:'
    do i = 1,ads%sz
      write(*,*)PRINTRANK,i,'row=',ads%F(i,:)
    enddo
#endif


  call mpi_barrier(MPI_COMM_WORLD,ierr)

end subroutine




! -------------------------------------------------------------------
! Deallocates all the resources and finalizes MPI.
! -------------------------------------------------------------------
subroutine Cleanup(ads)
use parallelism, ONLY : PRINTRANK
implicit none
type   (ADS_setup) :: ads
integer(kind=4) :: ierr

  if (allocated(ads%shiftsX)) deallocate(ads%shiftsX)
  if (allocated(ads%shiftsY)) deallocate(ads%shiftsY)
  if (allocated(ads%shiftsZ)) deallocate(ads%shiftsZ)
  
  if (allocated(ads%dimensionsX)) deallocate(ads%dimensionsX)
  if (allocated(ads%dimensionsX)) deallocate(ads%dimensionsY)
  if (allocated(ads%dimensionsZ)) deallocate(ads%dimensionsZ)

  if (allocated(ads%IPIVx)) deallocate(ads%IPIVx)
  if (allocated(ads%IPIVy)) deallocate(ads%IPIVy)
  if (allocated(ads%IPIVz)) deallocate(ads%IPIVz)
  
  if (allocated(ads%Ux)) deallocate(ads%Ux)
  if (allocated(ads%Uy)) deallocate(ads%Uy)
  if (allocated(ads%Uz)) deallocate(ads%Uz)
  
  if (allocated(ads%Mx)) deallocate(ads%Mx)
  if (allocated(ads%My)) deallocate(ads%My)
  if (allocated(ads%Mz)) deallocate(ads%Mz)
  
  if (allocated(ads%F))  deallocate(ads%F)
  if (allocated(ads%F2)) deallocate(ads%F2)
  if (allocated(ads%F3)) deallocate(ads%F3)

!!!!!! wyciac
  call mpi_finalize(ierr)
#ifdef IINFO
  write(*,*)PRINTRANK,"Exiting..."
#endif

end subroutine



!!!! przeniesc do debug
! -------------------------------------------------------------------
! Sanity-check of dimensions vector
! -------------------------------------------------------------------
subroutine ValidateDimensions(&
      n, &
      sx,sy,sz, &
      nrcppx,nrcppy,nrcppz, &
      dimensionsX,dimensionsY,dimensionsZ)
use parallelism, ONLY : NRPROCX,NRPROCY,NRPROCZ,PRINTRANK
implicit none
include "mpif.h"
integer(kind=4), intent(in), dimension(3) :: n
integer(kind=4), intent(in) :: sx,sy,sz
integer(kind=4), intent(in) :: nrcppx,nrcppy,nrcppz
integer(kind=4), intent(in), allocatable, dimension(:) :: dimensionsX
integer(kind=4), intent(in), allocatable, dimension(:) :: dimensionsY
integer(kind=4), intent(in), allocatable, dimension(:) :: dimensionsZ
   
integer(kind=4) :: i, k

  k = 0
  do i = 1,NRPROCX
    k = k + dimensionsX(i)
  enddo
  if (k /= (n(1)+1)*sy*sz) then
    write(*,*)PRINTRANK,'problem with dimensionsX',dimensionsX
    write(*,*)PRINTRANK,'nx+1',n(1)+1
    write(*,*)PRINTRANK,'sy',sy
    write(*,*)PRINTRANK,'sz',sz
    write(*,*)PRINTRANK,'nrcppx',nrcppx
    stop
  endif

  k = 0
  do i = 1,NRPROCY
    k = k + dimensionsY(i)
  enddo
  if (k /= (n(2)+1)*sx*sz) then
    write(*,*)PRINTRANK,'problem with dimensionsY',dimensionsY
    write(*,*)PRINTRANK,'n+1',n(2)+1
    write(*,*)PRINTRANK,'sx',sx
    write(*,*)PRINTRANK,'sz',sz
    stop
  endif

  k = 0
  do i = 1,NRPROCZ
    k = k + dimensionsZ(i)
  enddo
  if (k /= (n(3)+1)*sx*sy) then
    write(*,*)PRINTRANK,'problem with dimensionsZ',dimensionsZ
    write(*,*)PRINTRANK,'n+1',n(3)+1
    write(*,*)PRINTRANK,'sx',sx
    write(*,*)PRINTRANK,'sy',sy
    stop
  endif

end subroutine



!!!!! przeniesc do debug
! -------------------------------------------------------------------
! Displays computed domain decomposition, for debugging.
! -------------------------------------------------------------------
subroutine PrintDecompositionInfo(&
      n, &
      nrcppx,nrcppy,nrcppz, &
      ibegx,ibegy,ibegz, &
      iendx,iendy,iendz)
use parallelism, ONLY : NRPROCX,NRPROCY,NRPROCZ,PRINTRANK, &
MYRANKX,MYRANKY,MYRANKZ
implicit none
integer(kind=4), intent(in), dimension(3) :: n
integer(kind=4), intent(in) :: nrcppx,nrcppy,nrcppz
integer(kind=4), intent(in) :: ibegx,iendx
integer(kind=4), intent(in) :: ibegy,iendy
integer(kind=4), intent(in) :: ibegz,iendz

  write(*,*)PRINTRANK,'MYRANKX,MYRANKY,MYRANKZ',MYRANKX,MYRANKY,MYRANKZ
  write(*,*)PRINTRANK,'NRPROCX,NRPROCY,NRPROCZ',NRPROCX,NRPROCY,NRPROCZ
  write(*,*)PRINTRANK,'nx+1',n(1)+1
  write(*,*)PRINTRANK,'ny+1',n(2)+1
  write(*,*)PRINTRANK,'nz+1',n(3)+1
  write(*,*)PRINTRANK,'nrcppx,nrcppy,nrcppz',nrcppx,nrcppy,nrcppz
  write(*,*)PRINTRANK,'ibegx,iendx',ibegx,iendx
  write(*,*)PRINTRANK,'ibegy,iendy',ibegy,iendy
  write(*,*)PRINTRANK,'ibegz,iendz',ibegz,iendz

end subroutine


! -------------------------------------------------------------------
! Gathers full solution and plots it
! -------------------------------------------------------------------
subroutine PrintSolution(iter, t, ads)
use parallelism, ONLY : MYRANK
use plot, ONLY : SaveSplinePlot,PlotParams
use vtk, ONLY : VtkOutput
use my_mpi, ONLY : GatherFullSolution
implicit none
type   (ADS_setup), intent(in) :: ads
integer(kind=4), intent(in) :: iter
real   (kind=8), intent(in) :: t
real   (kind=8), allocatable :: solution(:,:,:)
type (PlotParams) :: params
character(len=20) :: filename

  call GatherFullSolution(0,ads%F,solution, &
         ads%n(1),ads%n(2),ads%n(3),ads%p(1),ads%p(2),ads%p(3),ads%sx,ads%sy,ads%sz)

  if (MYRANK == 0) then
    write(filename,'(I10)') iter
    filename = 'step' // adjustl(filename)
    ! filename = trim(filename) // '_'

    params = PlotParams(0,1,0,1,0,1,31,31,31)
    call SaveSplinePlot(trim(filename), &
      ads%Ux,ads%p(1),ads%n(1),ads%nelemx, &
      ads%Uy,ads%p(2),ads%n(2),ads%nelemy, &
      ads%Uz,ads%p(3),ads%n(3),ads%nelemz, &
      ! solution, GnuPlotOutput, params)
      solution, VtkOutput, params)

    ! call SavePlot(trim(filename), ftest, GnuPlotOutput, params)
  endif

end subroutine


   
end module
