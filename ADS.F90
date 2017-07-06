module ADSS

type ADS_setup
   ! Number of functions in each dimension minus one
   integer(kind=4) :: nx
   integer(kind=4) :: ny
   integer(kind=4) :: nz

   ! Degree of approximation
   integer(kind=4) :: px
   integer(kind=4) :: py
   integer(kind=4) :: pz

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
implicit none
include "mpif.h"
integer(kind=4), intent(in) :: nx,ny,nz
integer(kind=4), intent(in) :: px,py,pz
type(ADS_setup), intent(out) :: ads
integer(kind=4) :: ierr

  ads%px = px ! order
  ads%py = py ! order
  ads%pz = pz ! order
  ads%nx = nx  ! intervals
  ads%ny = ny  ! intervals
  ads%nz = nz  ! intervals

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
      nx,ny,nz, &
      px,py,pz, &
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
      nx,ny,nz, &
      ads%sx,ads%sy,ads%sz, &
      ads%nrcppx,ads%nrcppy,ads%nrcppz, &
      ads%dimensionsX,ads%dimensionsY,ads%dimensionsZ)
#endif

#ifdef IPRINT
    call PrintDecompositionInfo(&
      nx,ny,nz, &
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
   nx,ny,nz, &
   px,py,pz, &
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
   NRPROCX,NRPROCY,NRPROCZ
use utils, ONLY : ComputeEndpoints,FillDimVector
implicit none
integer(kind=4), intent(in) :: nx,ny,nz
integer(kind=4), intent(in) :: px,py,pz
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
  call ComputeEndpoints(MYRANKX, NRPROCX, nx, px, nrcppx, ibegx, iendx, minex, maxex)
  call ComputeEndpoints(MYRANKY, NRPROCY, ny, py, nrcppy, ibegy, iendy, miney, maxey)
  call ComputeEndpoints(MYRANKZ, NRPROCZ, nz, pz, nrcppz, ibegz, iendz, minez, maxez)

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
  call FillDimVector(dimensionsX, shiftsX, nrcppx, sy*sz, nx, NRPROCX)
  call FillDimVector(dimensionsY, shiftsY, nrcppy, sx*sz, ny, NRPROCY)
  call FillDimVector(dimensionsZ, shiftsZ, nrcppz, sx*sy, nz, NRPROCZ)

  ! Compute indices for neighbours
  ibegsx = -1
  iendsx = -1
  ibegsy = -1
  iendsy = -1
  ibegsz = -1
  iendsz = -1

  do i = max(MYRANKX-1,0)+1, min(MYRANKX+1,NRPROCX-1)+1
    ix = i-MYRANKX+1
    call ComputeEndpoints(i-1, NRPROCX, nx, px, nrcppx, ibegsx(ix), iendsx(ix), mine, maxe)
  enddo
  do i = max(MYRANKY-1,0)+1, min(MYRANKY+1,NRPROCY-1)+1
    iy = i-MYRANKY+1
    call ComputeEndpoints(i-1,NRPROCY, ny, py, nrcppy, ibegsy(iy), iendsy(iy), mine, maxe)
  enddo
  do i = max(MYRANKZ-1,0)+1, min(MYRANKZ+1,NRPROCZ-1)+1
    iz = i-MYRANKZ+1
    call ComputeEndpoints(i-1,NRPROCZ, nz, pz, nrcppz, ibegsz(iz), iendsz(iz), mine, maxe)
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


! -------------------------------------------------------------------
! Allocates and fills the knot vector
! -------------------------------------------------------------------
subroutine PrepareKnot(U,n,p,nelem)
use utils, ONLY : FillOpenKnot
use basis, ONLY : CountSpans
implicit none
integer(kind=4), intent(in) :: n,p
real   (kind=8), allocatable, dimension(:), intent(out) :: U
integer(kind=4), intent(out) :: nelem


  allocate(U(n+p+2))
  call FillOpenKnot(U, n, p)
  nelem = CountSpans(n,p,U)

#ifdef IINFO
    write(*,*)'n,p,nelem',n,p,nelem
    write(*,*)'U',U
#endif

end subroutine


! -------------------------------------------------------------------
! Calculates mass matrix M
! -------------------------------------------------------------------
subroutine ComputeMassMatrix(KL,KU,U,p,n,nelem,M)
use parallelism, ONLY : PRINTRANK
use projection_engine, ONLY : Form1DMassMatrix
implicit none
integer(kind=4), intent(in)  :: KL,KU
integer(kind=4), intent(in)  :: n, p, nelem
real   (kind=8), intent(in)  :: U(0:n+p+1)
real   (kind=8), intent(out) :: M(0:(2*KL+KU),0:n)
integer :: i

  call Form1DMassMatrix(KL,KU,U,p,n,nelem,M)
#ifdef IPRINT
    write(*,*)PRINTRANK,'M'
    do i = 1,2*KL+KU+1
      write(*,*)PRINTRANK,M(i,1:n+1)
    enddo
#endif

end subroutine




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


subroutine send_piece(items, dst, req)
implicit none
include "mpif.h"
real (kind=8) :: items(:)
type   (ADS_setup) :: ads
integer :: dst, req
integer :: ierr

  call mpi_isend(items, ads%nrcppz*ads%nrcppx*ads%nrcppy, &
    MPI_DOUBLE_PRECISION, dst, 0, MPI_COMM_WORLD, req, ierr)

end subroutine


subroutine recv_piece(items, src, req,ads)
implicit none
include "mpif.h"
real   (kind=8) :: items(:)
type   (ADS_setup) :: ads
integer(kind=4) :: src, req
integer(kind=4) :: ierr

  call mpi_irecv(items, ads%nrcppz*ads%nrcppx*ads%nrcppy, &
    MPI_DOUBLE_PRECISION, src, 0, MPI_COMM_WORLD, req, ierr)

end subroutine


! -------------------------------------------------------------------
! Distributes spline (e.g. solution from previous timestep, or parametrization
! approximation for Jacobian calculation) to neighbouring processes. It is
! essential that each process possess enough of the solution to be able to
! calculate values near the boundary, hence overlapping supports of B-splines
! necessitate partial sharing.
! -------------------------------------------------------------------
subroutine DistributeSpline(spline,ads)
use parallelism, ONLY : MYRANKX,MYRANKY,MYRANKZ,NRPROCX,NRPROCY,NRPROCZ
implicit none
include "mpif.h"
type   (ADS_setup) :: ads
real   (kind=8) :: spline(:,:,:,:)
integer(kind=4) :: i, j, k, s
integer(kind=4) :: request(3*3*3*2), stat(MPI_STATUS_SIZE)
integer(kind=4) :: ierr(3*3*3*2)
integer(kind=4) :: dst, src

  s = 1

  ! Right
  if (MYRANKX < NRPROCX - 1) then
    dst = neighbour(1, 0, 0)
    call send_piece(spline(:,2,2,2), dst, request(s))
    s = s + 1
  endif
  if (MYRANKX > 0) then
    src = neighbour(-1, 0, 0)
    call recv_piece(spline(:,1,2,2), src, request(s),ads)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Up
  if (MYRANKY > 0) then
    dst = neighbour(0, -1, 0)
    call send_piece(spline(:,2,2,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,1,2,2), dst, request(s))
    s = s + 1
  endif
  if (MYRANKY < NRPROCY - 1) then
    src = neighbour(0, 1, 0)
    call recv_piece(spline(:,2,3,2), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,1,3,2), src, request(s),ads)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Left
  if (MYRANKX > 0) then
    dst = neighbour(-1, 0, 0)
    call send_piece(ads%R(:,2,2,2), dst, request(s))
    s = s + 1
    call send_piece(ads%R(:,2,3,2), dst, request(s))
    s = s + 1
  endif
  if (MYRANKX < NRPROCX - 1) then
    src = neighbour(1, 0, 0)
    call recv_piece(ads%R(:,3,2,2), src, request(s),ads)
    s = s + 1
    call recv_piece(ads%R(:,3,3,2), src, request(s),ads)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Above
  if (MYRANKZ < NRPROCZ - 1) then
    dst = neighbour(0, 0, 1)
    call send_piece(spline(:,2,2,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,1,2,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,1,3,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,2,3,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,3,3,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,3,2,2), dst, request(s))
    s = s + 1
  endif
  if (MYRANKZ > 0) then
    src = neighbour(0, 0, -1)
    call recv_piece(spline(:,2,2,1), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,1,2,1), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,1,3,1), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,2,3,1), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,3,3,1), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,3,2,1), src, request(s),ads)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Down
  if (MYRANKY < NRPROCY - 1) then
    dst = neighbour(0, 1, 0)
    call send_piece(spline(:,2,2,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,1,2,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,3,2,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,1,2,1), dst, request(s))
    s = s + 1
    call send_piece(spline(:,2,2,1), dst, request(s))
    s = s + 1
    call send_piece(spline(:,3,2,1), dst, request(s))
    s = s + 1
  endif
  if (MYRANKY > 0) then
    src = neighbour(0, -1, 0)
    call recv_piece(spline(:,2,1,2), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,1,1,2), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,3,1,2), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,1,1,1), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,2,1,1), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,3,1,1), src, request(s),ads)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Below
  if (MYRANKZ > 0) then
    dst = neighbour(0, 0, -1)
    call send_piece(spline(:,1,1,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,1,2,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,1,3,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,2,1,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,2,2,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,2,3,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,3,1,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,3,2,2), dst, request(s))
    s = s + 1
    call send_piece(spline(:,3,3,2), dst, request(s))
    s = s + 1
  endif
  if (MYRANKZ < NRPROCZ - 1) then
    src = neighbour(0, 0, 1)
    call recv_piece(spline(:,1,1,3), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,1,2,3), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,1,3,3), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,2,1,3), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,2,2,3), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,2,3,3), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,3,1,3), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,3,2,3), src, request(s),ads)
    s = s + 1
    call recv_piece(spline(:,3,3,3), src, request(s),ads)
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1


  call mpi_barrier(MPI_COMM_WORLD,ierr)

end subroutine


! -------------------------------------------------------------------
! Prints debugging information about results of distributing
! data to neighbouring processes.
! -------------------------------------------------------------------
subroutine PrintDistributedData(ads)
use parallelism, ONLY : MYRANKX,MYRANKY,MYRANKZ,PRINTRANK, &
   NRPROCX,NRPROCY,NRPROCZ
use utils, ONLY : ComputeEndpoints
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

        call ComputeEndpoints(i-1,NRPROCX,ads%nx,ads%px,ads%nrcppx,obegx,oendx,mine,maxe)
        call ComputeEndpoints(j-1,NRPROCY,ads%ny,ads%py,ads%nrcppy,obegy,oendy,mine,maxe)
        call ComputeEndpoints(k-1,NRPROCZ,ads%nz,ads%pz,ads%nrcppz,obegz,oendz,mine,maxe)

        write(*,*) reshape(                                 &
          ads%R(:,i-MYRANKX+1,j-MYRANKY+1,k-MYRANKZ+1),         &
          (/ oendz-obegz+1,oendx-obegx+1,oendy-obegy+1 /))
      enddo
    enddo
  enddo

  write(*,*)'----'

end subroutine


! -------------------------------------------------------------------
! Calculates force (load)
!
! t - current time
! -------------------------------------------------------------------
subroutine ComputeRHS(iter,RHS_fun,ads)
use parallelism, ONLY : MYRANK,PRINTRANK
use projection_engine, ONLY : Form3DRHS
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
integer(kind=4), intent(in) :: iter
integer(kind=4) :: ierr,i

  call Form3DRHS                                    &
      (ads%Ux,ads%px,ads%nx,ads%nelemx,ads%nrcppx,  &
       ads%Uy,ads%py,ads%ny,ads%nelemy,ads%nrcppy,  &
       ads%Uz,ads%pz,ads%nz,ads%nelemz,ads%nrcppz,  &
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

  call mpi_barrier(MPI_COMM_WORLD,ierr)

end subroutine




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


! -------------------------------------------------------------------
! Performs one step of the simulation
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
subroutine Step(iter,RHS_fun,ads)
use parallelism, ONLY :PRINTRANK,MYRANKX,MYRANKY,MYRANKZ
use communicators, ONLY : COMMX,COMMY,COMMZ
use utils, ONLY : Gather,Scatter
use reorderRHS, ONLY : ReorderRHSForX,ReorderRHSForY,ReorderRHSForZ
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
  call ComputeRHS(iter,RHS_fun,ads)

  !--------------------------------------------------------------------
  ! Solve the first problem
  !--------------------------------------------------------------------
  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'1a) GATHER'
#endif

  allocate(ads%F_out((ads%nx+1),ads%sy*ads%sz))

  call Gather(ads%F,ads%F_out,ads%nx,ads%sx,ads%sy*ads%sz,ads%dimensionsX,ads%shiftsX,COMMX,ierr)

#ifdef IPRINT
    write(*,*)PRINTRANK,'after call mpi_gather'
    write(*,*)PRINTRANK,'ierr',ierr
    write(*,*)PRINTRANK,'F_out:'
    do i=1,ads%nx+1
      write(*,*)PRINTRANK,i,'row=',ads%F_out(i,1:ads%sy*ads%sz)
    enddo
#endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (MYRANKX == 0) then
#ifdef IINFO
     write(*,*)PRINTRANK,'1b) SOLVE THE FIRST PROBLEM'
#endif

    call ComputeMassMatrix(ads%KLx,ads%KUx,ads%Ux,ads%px,ads%nx,ads%nelemx,ads%Mx)
    call SolveOneDirection(ads%F_out, ads%sy*ads%sz,ads%nx,ads%KLx,ads%KUx,ads%px,ads%Mx,ads%IPIVx)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'1c) SCATTER'
#endif
  allocate(ads%F2_out(ads%sx,ads%sy*ads%sz))
  call Scatter(ads%F_out,ads%F2_out,ads%nx,ads%sx,ads%sy*ads%sz,ads%dimensionsX,ads%shiftsX,COMMX,ierr)
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

  allocate(ads%F2_out((ads%ny+1),ads%sx*ads%sz))
  call Gather(ads%F2,ads%F2_out,ads%ny,ads%sy,ads%sx*ads%sz,ads%dimensionsY,ads%shiftsY,COMMY,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (MYRANKY == 0) then
#ifdef IINFO
     write(*,*)PRINTRANK,'2b) SOLVE THE SECOND PROBLEM'
#endif

    call ComputeMassMatrix(ads%KLy,ads%KUy,ads%Uy,ads%py,ads%ny,ads%nelemy,ads%My)
    call SolveOneDirection(ads%F2_out, ads%sx*ads%sz,ads%ny,ads%KLy,ads%KUy,ads%py,ads%My,ads%IPIVy)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'2c) SCATHER'
#endif

  ! CORRECTION
  allocate(ads%F3_out(ads%sy,ads%sx*ads%sz))
  call Scatter(ads%F2_out,ads%F3_out,ads%ny,ads%sy,ads%sx*ads%sz,ads%dimensionsY,ads%shiftsY,COMMY,ierr)
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
  allocate(ads%F3_out((ads%nz+1),ads%sx*ads%sy))
  call Gather(ads%F3,ads%F3_out,ads%nz,ads%sz,ads%sx*ads%sy,ads%dimensionsZ,ads%shiftsZ,COMMZ,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (MYRANKZ == 0) then
#ifdef IINFO
     write(*,*)PRINTRANK,'3b) SOLVE THE THIRD PROBLEM'
#endif

    call ComputeMassMatrix(ads%KLz,ads%KUz,ads%Uz,ads%pz,ads%nz,ads%nelemz,ads%Mz)
    call SolveOneDirection(ads%F3_out, ads%sx*ads%sy,ads%nz,ads%KLz,ads%KUz,ads%pz,ads%Mz,ads%IPIVz)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifdef IINFO
  write(*,*)PRINTRANK,'3c) SCATTER'
#endif

  ! CORRECTION
  allocate(ads%F_out(ads%sz,ads%sx*ads%sy))
  call Scatter(ads%F3_out,ads%F_out,ads%nz,ads%sz,ads%sx*ads%sy,ads%dimensionsZ,ads%shiftsZ,COMMZ,ierr)
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
  call DistributeSpline(ads%R,ads)

#ifdef IPRINT
    write(*,*)PRINTRANK,'Result:'
    do i = 1,ads%sz
      write(*,*)PRINTRANK,i,'row=',ads%F(i,:)
    enddo
#endif


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
subroutine GatherFullSolution(at, part, full, ads)
use parallelism, ONLY : MYRANK,LINEARINDEX,NRPROCX,NRPROCY,NRPROCZ
use utils, ONLY : ComputeEndpoints
implicit none
include "mpif.h"
type   (ADS_setup) :: ads
integer(kind=4), intent(in) :: at
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
    array_size = (ads%nx+1)*(ads%ny+1)*(ads%nz+1)
    allocate(full(0:ads%nx,0:ads%ny,0:ads%nz))
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
        size = SizeOfPiece(x,y,z,ads%nx,ads%ny,ads%nz,ads%px,ads%py,ads%pz)
        recvcounts(idx) = size
        displs(idx) = offset
        offset = offset + size
      enddo
    enddo
  enddo

  call mpi_gatherv(part, ads%sx*ads%sy*ads%sz, MPI_DOUBLE_PRECISION, buffer, &
    recvcounts, displs, MPI_DOUBLE_PRECISION, at, MPI_COMM_WORLD, ierr)

  ! Reordering of the array at root
  if (MYRANK == at) then
    offset = 0
    do x = 0, NRPROCX-1
      do y = 0, NRPROCY-1
        do z = 0, NRPROCZ-1
          call ComputeEndpoints(x, NRPROCX, ads%nx, ads%px, nrcpp, begx, endx, mine, maxe)
          call ComputeEndpoints(y, NRPROCY, ads%ny, ads%py, nrcpp, begy, endy, mine, maxe)
          call ComputeEndpoints(z, NRPROCZ, ads%nz, ads%pz, nrcpp, begz, endz, mine, maxe)
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

          offset = offset + SizeOfPiece(x,y,z,ads%nx,ads%ny,ads%nz,ads%px,ads%py,ads%pz)
        enddo
      enddo
    enddo
  endif

  deallocate(buffer)

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

  call mpi_finalize(ierr)
#ifdef IINFO
  write(*,*)PRINTRANK,"Exiting..."
#endif

end subroutine


! -------------------------------------------------------------------
! Sanity-check of dimensions vector
! -------------------------------------------------------------------
subroutine ValidateDimensions(&
      nx,ny,nz, &
      sx,sy,sz, &
      nrcppx,nrcppy,nrcppz, &
      dimensionsX,dimensionsY,dimensionsZ)
use parallelism, ONLY : NRPROCX,NRPROCY,NRPROCZ,PRINTRANK
implicit none
include "mpif.h"
integer(kind=4), intent(in) :: nx,ny,nz
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
  if (k /= (nx+1)*sy*sz) then
    write(*,*)PRINTRANK,'problem with dimensionsX',dimensionsX
    write(*,*)PRINTRANK,'nx+1',nx+1
    write(*,*)PRINTRANK,'sy',sy
    write(*,*)PRINTRANK,'sz',sz
    write(*,*)PRINTRANK,'nrcppx',nrcppx
    stop
  endif

  k = 0
  do i = 1,NRPROCY
    k = k + dimensionsY(i)
  enddo
  if (k /= (ny+1)*sx*sz) then
    write(*,*)PRINTRANK,'problem with dimensionsY',dimensionsY
    write(*,*)PRINTRANK,'n+1',ny+1
    write(*,*)PRINTRANK,'sx',sx
    write(*,*)PRINTRANK,'sz',sz
    stop
  endif

  k = 0
  do i = 1,NRPROCZ
    k = k + dimensionsZ(i)
  enddo
  if (k /= (nz+1)*sx*sy) then
    write(*,*)PRINTRANK,'problem with dimensionsZ',dimensionsZ
    write(*,*)PRINTRANK,'n+1',nz+1
    write(*,*)PRINTRANK,'sx',sx
    write(*,*)PRINTRANK,'sy',sy
    stop
  endif

end subroutine


! -------------------------------------------------------------------
! Displays computed domain decomposition, for debugging.
! -------------------------------------------------------------------
subroutine PrintDecompositionInfo(&
      nx,ny,nz, &
      nrcppx,nrcppy,nrcppz, &
      ibegx,ibegy,ibegz, &
      iendx,iendy,iendz)
use parallelism, ONLY : NRPROCX,NRPROCY,NRPROCZ,PRINTRANK, &
MYRANKX,MYRANKY,MYRANKZ
implicit none
integer(kind=4), intent(in) :: nx
integer(kind=4), intent(in) :: ny
integer(kind=4), intent(in) :: nz
integer(kind=4), intent(in) :: nrcppx,nrcppy,nrcppz
integer(kind=4), intent(in) :: ibegx,iendx
integer(kind=4), intent(in) :: ibegy,iendy
integer(kind=4), intent(in) :: ibegz,iendz

  write(*,*)PRINTRANK,'MYRANKX,MYRANKY,MYRANKZ',MYRANKX,MYRANKY,MYRANKZ
  write(*,*)PRINTRANK,'NRPROCX,NRPROCY,NRPROCZ',NRPROCX,NRPROCY,NRPROCZ
  write(*,*)PRINTRANK,'nx+1',nx+1
  write(*,*)PRINTRANK,'ny+1',ny+1
  write(*,*)PRINTRANK,'nz+1',nz+1
  write(*,*)PRINTRANK,'nrcppx,nrcppy,nrcppz',nrcppx,nrcppy,nrcppz
  write(*,*)PRINTRANK,'ibegx,iendx',ibegx,iendx
  write(*,*)PRINTRANK,'ibegy,iendy',ibegy,iendy
  write(*,*)PRINTRANK,'ibegz,iendz',ibegz,iendz

end subroutine


function ftest(x, y, z) result(val)
real (kind=8) :: x, y, z, val

  val = x**2 + y**2 + z**2

end function


! -------------------------------------------------------------------
! Gathers full solution and plots it
! -------------------------------------------------------------------
subroutine PrintSolution(iter, t, ads)
use parallelism, ONLY : MYRANK
use plot, ONLY : SaveSplinePlot,PlotParams
use vtk, ONLY : VtkOutput
implicit none
type   (ADS_setup), intent(in) :: ads
integer(kind=4), intent(in) :: iter
real   (kind=8), intent(in) :: t
real   (kind=8), allocatable :: solution(:,:,:)
type (PlotParams) :: params
character(len=20) :: filename

  call GatherFullSolution(0, ads%F, solution,ads)

  if (MYRANK == 0) then
    write(filename,'(I10)') iter
    filename = 'step' // adjustl(filename)
    ! filename = trim(filename) // '_'

    params = PlotParams(0,1,0,1,0,1,31,31,31)
    call SaveSplinePlot(trim(filename), &
      ads%Ux,ads%px,ads%nx,ads%nelemx, &
      ads%Uy,ads%py,ads%ny,ads%nelemy, &
      ads%Uz,ads%pz,ads%nz,ads%nelemz, &
      ! solution, GnuPlotOutput, params)
      solution, VtkOutput, params)

    ! call SavePlot(trim(filename), ftest, GnuPlotOutput, params)
  endif

end subroutine


   
end module
