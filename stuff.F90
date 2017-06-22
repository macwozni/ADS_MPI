module stuff

use parallelism
use projection_engine
use communicators
use utils
use debug
use plot
use gnuplot
use vtk
use basis
use reorderRHS
use input_data

implicit none

! Number of functions in each dimension minus one
integer(kind=4) :: nx
integer(kind=4) :: ny
integer(kind=4) :: nz

! Degree of approximation
integer(kind=4) :: px
integer(kind=4) :: py
integer(kind=4) :: pz

! Number of iterations
integer :: steps

! Time and timestep
real (kind=8) :: Dt


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

! Statistics computed during the simulation
real (kind=8) :: drained = 0, pollution = 0

! Buffer for coefficients of solution corresponding to neighbouring
! parts of the domain. It is (Nx*Ny*Nz) x 3 x 3 x 3 array, where
! Nx*Ny*Nz is the size of part of solution for one fragment of domain.
real (kind=8), allocatable :: R(:,:,:,:)

! Buffer for values of permeability function
real (kind=8), allocatable :: Kqvals(:,:,:,:,:,:)

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

! order of approximations
integer(kind=4) :: ORDER

! number of elements in one dimension
integer(kind=4) :: SIZE

contains


! -------------------------------------------------------------------
! Sets values of parameters (order and size)
! -------------------------------------------------------------------
subroutine InitializeParameters
implicit none
character(100) :: input

  ! ./l2 <size> <procx> <procy> <procz> <nsteps> <dt>
  ORDER = 2

  call getarg(1,input)
  read(input,*) SIZE
  call getarg(2,input)
  read(input,*) NRPROCX
  call getarg(3,input)
  read(input,*) NRPROCY
  call getarg(4,input)
  read(input,*) NRPROCZ
  call getarg(5,input)
  read(input,*) steps
  call getarg(6,input)
  read(input,*) Dt

end subroutine


! -------------------------------------------------------------------
! Initialization of clocks and MPI
! -------------------------------------------------------------------
subroutine Initialize
implicit none
include "mpif.h"
integer(kind=4) :: ierr

  call InitializeParameters
  call InitializeParallelism
  call CreateCommunicators

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iinfo == 1) write(*,*)PRINTRANK,'INITIALIZATION'

end subroutine


! -------------------------------------------------------------------
! Establishes decomposition of the domain. Calculates size and location
! of the piece for current process.
! -------------------------------------------------------------------
subroutine ComputeDecomposition
implicit none
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

  if (iinfo == 1) then
    write(*,*)PRINTRANK,'Number of cols per processor:',nrcppx,nrcppy,nrcppz
    write(*,*)PRINTRANK,'ibegx,iendx',ibegx,iendx
    write(*,*)PRINTRANK,'ibegy,iendy',ibegy,iendy
    write(*,*)PRINTRANK,'ibegz,iendz',ibegz,iendz
  endif

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
subroutine AllocateArrays
implicit none
include "mpif.h"
integer :: ierr

  allocate(Mx(2*KLx+KUx+1,nx+1))
  allocate(My(2*KLy+KUy+1,ny+1))
  allocate(Mz(2*KLz+KUz+1,nz+1))

  ! OLD: MP start with system fully generated along X
  ! allocate( F((n+1),(sy)*(sz))) !x,y,z
  allocate( F(sx,sy*sz)) !x,y,z
  allocate(F2(sy,sx*sz)) !y,x,z
  allocate(F3(sz,sx*sy)) !z,x,y

  allocate(Kqvals(px+1,py+1,pz+1,maxex-minex+1,maxey-miney+1,maxez-minez+1))

  ! Processes on the border need pivot vector for LAPACK call
  if (MYRANKX == 0 .or. MYRANKY == 0 .or. MYRANKZ == 0) then
    allocate(IPIVx(nx+1))
    allocate(IPIVy(ny+1))
    allocate(IPIVz(nz+1))
  endif

  allocate(R(nrcppz*nrcppx*nrcppy,3,3,3))

  call mpi_barrier(MPI_COMM_WORLD,ierr)

end subroutine


! -------------------------------------------------------------------
! Allocates and fills the knot vector
! -------------------------------------------------------------------
subroutine PrepareKnot(U,n,p)
implicit none
integer(kind=4), intent(in) :: n,p
real   (kind=8), allocatable, dimension(:), intent(out) :: U
integer(kind=4) :: nelem


  allocate(U(n+p+2))
  call FillOpenKnot(U, n, p)
  nelem = CountSpans(n,p,U)

  if (iinfo == 1) then
    write(*,*)'n,p,nelem',n,p,nelem
    write(*,*)'U',U
  endif

end subroutine


! -------------------------------------------------------------------
! Calculates mass matrix M
! -------------------------------------------------------------------
subroutine ComputeMassMatrix(KL,KU,U,p,n,nelem,M)
implicit none
integer(kind=4), intent(in)  :: KL,KU
integer(kind=4), intent(in)  :: n, p, nelem
real   (kind=8), intent(in)  :: U(0:n+p+1)
real   (kind=8), intent(out) :: M(0:(2*KL+KU),0:n)
integer :: i

  call Form1DMassMatrix(KL,KU,U,p,n,nelem,M)
  if (iprint == 1) then
    write(*,*)PRINTRANK,'M'
    do i = 1,2*KL+KU+1
      write(*,*)PRINTRANK,M(i,1:n+1)
    enddo
  endif

end subroutine


subroutine PrecomputeKq
implicit none

  call CacheKqValues                                &
      (Ux,px,nx,minex,maxex,nelemx,                     &
       Uy,py,ny,miney,maxey,nelemy,                     &
       Uz,pz,nz,minez,maxez,nelemz,                     &
       Kqvals)

end subroutine


function neighbour(dx, dy, dz) result(idx)
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
integer :: dst, req
integer :: ierr

  call mpi_isend(items, nrcppz*nrcppx*nrcppy, &
    MPI_DOUBLE_PRECISION, dst, 0, MPI_COMM_WORLD, req, ierr)

end subroutine


subroutine recv_piece(items, src, req)
implicit none
include "mpif.h"
real   (kind=8) :: items(:)
integer(kind=4) :: src, req
integer(kind=4) :: ierr

  call mpi_irecv(items, nrcppz*nrcppx*nrcppy, &
    MPI_DOUBLE_PRECISION, src, 0, MPI_COMM_WORLD, req, ierr)

end subroutine


! -------------------------------------------------------------------
! Distributes spline (e.g. solution from previous timestep, or parametrization
! approximation for Jacobian calculation) to neighbouring processes. It is
! essential that each process possess enough of the solution to be able to
! calculate values near the boundary, hence overlapping supports of B-splines
! necessitate partial sharing.
! -------------------------------------------------------------------
subroutine DistributeSpline(spline)
implicit none
include "mpif.h"
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
    call recv_piece(spline(:,1,2,2), src, request(s))
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
    call recv_piece(spline(:,2,3,2), src, request(s))
    s = s + 1
    call recv_piece(spline(:,1,3,2), src, request(s))
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1

  ! Left
  if (MYRANKX > 0) then
    dst = neighbour(-1, 0, 0)
    call send_piece(R(:,2,2,2), dst, request(s))
    s = s + 1
    call send_piece(R(:,2,3,2), dst, request(s))
    s = s + 1
  endif
  if (MYRANKX < NRPROCX - 1) then
    src = neighbour(1, 0, 0)
    call recv_piece(R(:,3,2,2), src, request(s))
    s = s + 1
    call recv_piece(R(:,3,3,2), src, request(s))
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
    call recv_piece(spline(:,2,2,1), src, request(s))
    s = s + 1
    call recv_piece(spline(:,1,2,1), src, request(s))
    s = s + 1
    call recv_piece(spline(:,1,3,1), src, request(s))
    s = s + 1
    call recv_piece(spline(:,2,3,1), src, request(s))
    s = s + 1
    call recv_piece(spline(:,3,3,1), src, request(s))
    s = s + 1
    call recv_piece(spline(:,3,2,1), src, request(s))
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
    call recv_piece(spline(:,2,1,2), src, request(s))
    s = s + 1
    call recv_piece(spline(:,1,1,2), src, request(s))
    s = s + 1
    call recv_piece(spline(:,3,1,2), src, request(s))
    s = s + 1
    call recv_piece(spline(:,1,1,1), src, request(s))
    s = s + 1
    call recv_piece(spline(:,2,1,1), src, request(s))
    s = s + 1
    call recv_piece(spline(:,3,1,1), src, request(s))
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
    call recv_piece(spline(:,1,1,3), src, request(s))
    s = s + 1
    call recv_piece(spline(:,1,2,3), src, request(s))
    s = s + 1
    call recv_piece(spline(:,1,3,3), src, request(s))
    s = s + 1
    call recv_piece(spline(:,2,1,3), src, request(s))
    s = s + 1
    call recv_piece(spline(:,2,2,3), src, request(s))
    s = s + 1
    call recv_piece(spline(:,2,3,3), src, request(s))
    s = s + 1
    call recv_piece(spline(:,3,1,3), src, request(s))
    s = s + 1
    call recv_piece(spline(:,3,2,3), src, request(s))
    s = s + 1
    call recv_piece(spline(:,3,3,3), src, request(s))
    s = s + 1
  endif

  do i = 1,s-1
    call mpi_wait(request(i),stat,ierr)
  enddo
  s = 1


  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !if (MYRANK == 0) then
  !  call PrintDistributedData
  !endif

  !do i=ibegx,iendx
  !  do j=ibegy,iendy
  !    do k=ibegz,iendz
  !       Result(k-ibegz+1,(j-ibegy)*sx + i-ibegx+1)=(10*i+j)*10 + k
  !    enddo
  !  enddo
  !enddo

end subroutine


! -------------------------------------------------------------------
! Prints debugging information about results of distributing
! data to neighbouring processes.
! -------------------------------------------------------------------
subroutine PrintDistributedData
implicit none
integer(kind=4) :: i, j, k
integer(kind=4) :: obegx,oendx,obegy,oendy,obegz,oendz
integer(kind=4) :: mine, maxe

  write(*,*)PRINTRANK,'R:'

  do i = max(MYRANKX-1,0)+1, min(MYRANKX+1,NRPROCX-1)+1
    do j = max(MYRANKY-1,0)+1, min(MYRANKY+1,NRPROCY-1)+1
      do k = max(MYRANKZ-1,0)+1, min(MYRANKZ+1,NRPROCZ-1)+1
        write(*,*)'(i,j,k)=',i+1,j+1,k+1

        call ComputeEndpoints(i-1,NRPROCX,nx,px,nrcppx,obegx,oendx,mine,maxe)
        call ComputeEndpoints(j-1,NRPROCY,ny,py,nrcppy,obegy,oendy,mine,maxe)
        call ComputeEndpoints(k-1,NRPROCZ,nz,pz,nrcppz,obegz,oendz,mine,maxe)

        write(*,*) reshape(                                 &
          R(:,i-MYRANKX+1,j-MYRANKY+1,k-MYRANKZ+1),         &
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
subroutine ComputeRHS(iter, t)
implicit none
include "mpif.h"
real   (kind=8) :: t
integer(kind=4) :: iter, i
integer(kind=4) :: ierr
real   (kind=8) :: l2norm, fullnorm

  call Form3DRHS                                    &
      (Ux,px,nx,nelemx,nrcppx,                          &
       Uy,py,ny,nelemy,nrcppy,                          &
       Uz,pz,nz,nelemz,nrcppz,                          &
       ibegx,iendx,                                 &
       ibegy,iendy,                                 &
       ibegz,iendz,                                 &
       ibegsx,iendsx,ibegsy,iendsy,ibegsz,iendsz,   &
       minex,maxex,miney,maxey,minez,maxez,         &
       Kqvals,Dt,t,R,F,drained,l2norm)

  if (iprint == 1) then
    write(*,*)PRINTRANK,'F'
    do i = 1,sx
      write(*,*)PRINTRANK,F(i,1:sy*sz)
    enddo
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call MPI_Reduce(l2norm, fullnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (MYRANK == 0) then
    write(*,*)iter, 'L2 norm:', fullnorm
  endif

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
real   (kind=8) :: RHS(:,:)
integer :: KL, KU
integer, dimension(:) :: IPIV
real (kind=8), dimension(:,:) :: M
integer :: n,p
integer(kind=4) :: eqnum
integer(kind=4) :: i, iret

  IPIV = 0

  if (iprint == 1) then
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
  endif

  ! SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
  ! .. Scalar Arguments ..
  ! INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
  ! .. Array Arguments ..
  ! INTEGER            IPIV( * )
  ! DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )

  call DGBSV(n+1,KL,KU,eqnum,M,2*KL+KU+1,IPIV,RHS,n+1,iret)

  if (iprint == 1) then
    write(*,*)'iret=',iret
    write(*,*)'Solution='
    do i = 1,n+1
      write(*,*)i,'row=',RHS(i,1:eqnum)
    enddo
  endif

end subroutine


! -------------------------------------------------------------------
! Performs one step of the simulation
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
subroutine Step(iter, t)
implicit none
include "mpif.h"
integer(kind=4) :: iter
real   (kind=8) :: t
integer(kind=4) :: i
integer(kind=4) :: iret, ierr

  ! generate the RHS vectors
  call ComputeRHS(iter, t)

  !--------------------------------------------------------------------
  ! Solve the first problem
  !--------------------------------------------------------------------
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iinfo == 1) write(*,*)PRINTRANK,'1a) GATHER'

  allocate(F_out((nx+1),sy*sz))

  call Gather(F,F_out,nx,sx,sy*sz,dimensionsX,shiftsX,COMMX,ierr)

  if (iprint == 1) then
    write(*,*)PRINTRANK,'after call mpi_gather'
    write(*,*)PRINTRANK,'ierr',ierr
    write(*,*)PRINTRANK,'F_out:'
    do i=1,nx+1
      write(*,*)PRINTRANK,i,'row=',F_out(i,1:sy*sz)
    enddo
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (MYRANKX == 0) then
    if (iinfo == 1) write(*,*)PRINTRANK,'1b) SOLVE THE FIRST PROBLEM'

    call ComputeMassMatrix(KLx,KUx,Ux,px,nx,nelemx,Mx)
    call SolveOneDirection(F_out, sy*sz,nx,KLx,KUx,px,Mx,IPIVx)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iinfo == 1) write(*,*)PRINTRANK,'1c) SCATTER'
  allocate(F2_out(sx,sy*sz))
  call Scatter(F_out,F2_out,nx,sx,sy*sz,dimensionsX,shiftsX,COMMX,ierr)
  deallocate(F_out)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iinfo == 1) write(*,*)PRINTRANK,'1d) REORDER'
  call ReorderRHSForY(ibegx,iendx,ibegy,iendy,ibegz,iendz,F2_out,F2)
  deallocate(F2_out)

  if (iprint == 1) then
    write(*,*)PRINTRANK,'after ReorderRHSForY'
    write(*,*)PRINTRANK,'F2:'
    do i = 1,sy
      write(*,*)PRINTRANK,i,'row=',F2(i,1:sx*sz)
    enddo
  endif
  iprint = 0

  !--------------------------------------------------------------------
  ! Solve the second problem
  !--------------------------------------------------------------------
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iinfo == 1) write(*,*)PRINTRANK,'2a) GATHER'

  allocate(F2_out((ny+1),sx*sz))
  call Gather(F2,F2_out,ny,sy,sx*sz,dimensionsY,shiftsY,COMMY,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (MYRANKY == 0) then
    if (iinfo == 1) write(*,*)PRINTRANK,'2b) SOLVE THE SECOND PROBLEM'

    call ComputeMassMatrix(KLy,KUy,Uy,py,ny,nelemy,My)
    call SolveOneDirection(F2_out, sx*sz,ny,KLy,KUy,py,My,IPIVy)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iinfo == 1) write(*,*)PRINTRANK,'2c) SCATHER'

  ! CORRECTION
  allocate(F3_out(sy,sx*sz))
  call Scatter(F2_out,F3_out,ny,sy,sx*sz,dimensionsY,shiftsY,COMMY,ierr)
  deallocate(F2_out)

  if(iprint == 1)then
    write(*,*)PRINTRANK,'after call mpi_scatterv'
    write(*,*)PRINTRANK,'ierr',ierr
    write(*,*)PRINTRANK,'F3_out:'
    do i = 1,sy
      write(*,*)PRINTRANK,i,'row=',F3_out(i,1:sx*sz)
    enddo
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iinfo == 1) write(*,*)PRINTRANK,'2d) REORDER'
  ! Reorder right hand sides
  call ReorderRHSForZ(ibegx,iendx,ibegy,iendy,ibegz,iendz,F3_out,F3)
  deallocate(F3_out)

  if (iprint == 1) then
    write(*,*)PRINTRANK,'after ReorderRHSForZ'
    write(*,*)PRINTRANK,'F3:'
    do i = 1,sz
      write(*,*)PRINTRANK,i,'row=',F3(i,1:sx*sy)
    enddo
  endif
  iprint = 0

  !--------------------------------------------------------------------
  ! Solve the third problem
  !--------------------------------------------------------------------
  if (iinfo == 1) write(*,*)PRINTRANK,'3a) GATHER'

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  allocate(F3_out((nz+1),sx*sy))
  call Gather(F3,F3_out,nz,sz,sx*sy,dimensionsZ,shiftsZ,COMMZ,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (MYRANKZ == 0) then
    if (iinfo == 1) write(*,*)PRINTRANK,'3b) SOLVE THE THIRD PROBLEM'

    call ComputeMassMatrix(KLz,KUz,Uz,pz,nz,nelemz,Mz)
    call SolveOneDirection(F3_out, sx*sy,nz,KLz,KUz,pz,Mz,IPIVz)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iinfo == 1) write(*,*)PRINTRANK,'3c) SCATTER'

  ! CORRECTION
  allocate(F_out(sz,sx*sy))
  call Scatter(F3_out,F_out,nz,sz,sx*sy,dimensionsZ,shiftsZ,COMMZ,ierr)
  deallocate(F3_out)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iinfo == 1) write(*,*)PRINTRANK,'3d) REORDER'
  ! Reorder right hand sides
  call ReorderRHSForX(ibegx,iendx,ibegy,iendy,ibegz,iendz,F_out,F)
  deallocate(F_out)

  if (iprint == 1) then
    write(*,*)PRINTRANK,'after ReorderRHSForX'
    write(*,*)PRINTRANK,'F:'
    do i = 1,sx
      write(*,*)PRINTRANK,i,'row=',F(i,1:sy*sz)
    enddo
  endif
  iprint = 0

  if (iinfo == 1) write(*,*)PRINTRANK,'3e) DISTRIBUTE SOLUTION'
  R(1:sx*sy*sz,2,2,2) = reshape(F, [sx*sy*sz])
  call DistributeSpline(R)

  ! if (MYRANK == 0) iprint=1
  if (iprint == 1) then
    write(*,*)PRINTRANK,'Result:'
    do i = 1,sz
      write(*,*)PRINTRANK,i,'row=',F(i,:)
    enddo
  endif


  call mpi_barrier(MPI_COMM_WORLD,ierr)

end subroutine


! -------------------------------------------------------------------
! Calculates size of the piece corresponding to process with
! specified coordinates. Concretly, number of coefficients.
!
! x, y, z    - coordinates
! -------------------------------------------------------------------
function SizeOfPiece(x, y, z) result (s)
implicit none
integer(kind=4), intent(in) :: x, y, z
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
subroutine GatherFullSolution(at, part, full)
implicit none
include "mpif.h"
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
    array_size = (nx+1)*(ny+1)*(nz+1)
    allocate(full(0:nx,0:ny,0:nz))
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
        size = SizeOfPiece(x, y, z)
        recvcounts(idx) = size
        displs(idx) = offset
        offset = offset + size
      enddo
    enddo
  enddo

  call mpi_gatherv(part, sx*sy*sz, MPI_DOUBLE_PRECISION, buffer, &
    recvcounts, displs, MPI_DOUBLE_PRECISION, at, MPI_COMM_WORLD, ierr)

  ! Reordering of the array at root
  if (MYRANK == at) then
    offset = 0
    do x = 0, NRPROCX-1
      do y = 0, NRPROCY-1
        do z = 0, NRPROCZ-1
          call ComputeEndpoints(x, NRPROCX, nx, px, nrcpp, begx, endx, mine, maxe)
          call ComputeEndpoints(y, NRPROCY, ny, py, nrcpp, begy, endy, mine, maxe)
          call ComputeEndpoints(z, NRPROCZ, nz, pz, nrcpp, begz, endz, mine, maxe)
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

          offset = offset + SizeOfPiece(x, y, z)
        enddo
      enddo
    enddo
  endif

  deallocate(buffer)

end subroutine


! -------------------------------------------------------------------
! Deallocates all the resources and finalizes MPI.
! -------------------------------------------------------------------
subroutine Cleanup
implicit none
integer(kind=4) :: ierr

  deallocate(shiftsX)
  deallocate(shiftsY)
  deallocate(shiftsZ)
  deallocate(dimensionsX)
  deallocate(dimensionsY)
  deallocate(dimensionsZ)

  if (allocated(IPIVx)) deallocate(IPIVx)
  if (allocated(IPIVy)) deallocate(IPIVy)
  if (allocated(IPIVz)) deallocate(IPIVz)
  deallocate(Ux)
  deallocate(Uy)
  deallocate(Uz)
  deallocate(Mx)
  deallocate(My)
  deallocate(Mz)
  deallocate(F)

  call mpi_finalize(ierr)
  if (iinfo == 1) write(*,*)PRINTRANK,"Exiting..."

end subroutine


! -------------------------------------------------------------------
! Sanity-check of dimensions vector
! -------------------------------------------------------------------
subroutine ValidateDimensions
implicit none
include "mpif.h"
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
subroutine PrintDecompositionInfo
implicit none

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
subroutine PrintSolution(iter, t)
implicit none
integer(kind=4) :: iter
real   (kind=8) :: t
real   (kind=8), allocatable :: solution(:,:,:)
type (PlotParams) :: params
character(len=20) :: filename

  call GatherFullSolution(0, F, solution)

  if (MYRANK == 0) then
    write(filename,'(I10)') iter
    filename = 'step' // adjustl(filename)
    ! filename = trim(filename) // '_'

    params = PlotParams(0,1,0,1,0,1,31,31,31)
    call SaveSplinePlot(trim(filename), &
      Ux,px,nx,nelemx, &
      Uy,py,ny,nelemy, &
      Uz,pz,nz,nelemz, &
      ! solution, GnuPlotOutput, params)
      solution, VtkOutput, params)

    ! call SavePlot(trim(filename), ftest, GnuPlotOutput, params)
  endif

end subroutine


subroutine ComputeResults()
implicit none
include "mpif.h"
real   (kind=8) :: fulldrained
integer(kind=4) :: ierr


  call MPI_Reduce(drained, fulldrained, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (MYRANK == 0) then
    write(*,*) fulldrained
  endif

end subroutine


end module

 
