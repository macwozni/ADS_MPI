
module time_data

implicit none

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

end module


module stuff

use parallelism
use projection_engine
use communicators
use utils
use stopwatch
use time_data
use debug
use plot
use gnuplot
use vtk
use basis
use reorderRHS
use input_data

implicit none

! Number of functions in each dimension minus one
integer(kind=4) :: n

! Degree of approximation
integer(kind=4) :: p

! Number of iterations
integer :: steps

! Time and timestep
real (kind=8) :: Dt


! Knot vector
real (kind=8), allocatable, dimension(:) :: U

! Mass matrix
real (kind=8), allocatable, dimension(:,:) :: M

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
integer(kind=4) :: nelem

! Size of slices of domain in each dimension
integer, allocatable, dimension(:) :: dimensionsX
integer, allocatable, dimension(:) :: dimensionsY
integer, allocatable, dimension(:) :: dimensionsZ

! Offsets of slices of domain in each direction
integer, allocatable, dimension(:) :: shiftsX
integer, allocatable, dimension(:) :: shiftsY
integer, allocatable, dimension(:) :: shiftsZ

! Pivot array for these processes that need to solve systems
integer, allocatable, dimension(:) :: IPIV

! Number of lower and upper diagonal entries in mass matrix
integer :: KL, KU

! Number of columns (average) per processor
integer :: nrcppx,nrcppy,nrcppz

! Range of piece of domain assigned to this process
integer :: ibegx,iendx
integer :: ibegy,iendy
integer :: ibegz,iendz

! Size of piece of domain assigned to this process
integer :: sx,sy,sz

! Ranges of pieces of domain around the one assigned to this process
integer, dimension(3) :: ibegsx,iendsx
integer, dimension(3) :: ibegsy,iendsy
integer, dimension(3) :: ibegsz,iendsz

! Range of elements associated with basis functions assigned to this process
integer :: minex, maxex
integer :: miney, maxey
integer :: minez, maxez

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
integer :: ierr

  call start_clock(iclock)
  call start_clock(iclock_init)

  if (MYRANK == 0) then
    call start_clock(iclock_i1)
  endif

  call InitializeParameters
  call InitializeParallelism
  call CreateCommunicators

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (iinfo == 1) then
    call stop_clock(dtime_i1,iclock_i1)
    write(*,*)'create_communicators:',dtime_i1
    call start_clock(iclock_i2)
  endif

  ! if (MYRANK == 0) iinfo=1
  if (iinfo == 1) write(*,*)PRINTRANK,'INITIALIZATION'

end subroutine


! -------------------------------------------------------------------
! Establishes decomposition of the domain. Calculates size and location
! of the piece for current process.
! -------------------------------------------------------------------
subroutine ComputeDecomposition
implicit none
integer :: i
integer :: ix, iy, iz
integer :: mine, maxe

  ! number of columns per processors
  call ComputeEndpoints(MYRANKX, NRPROCX, n, p, nrcppx, ibegx, iendx, minex, maxex)
  call ComputeEndpoints(MYRANKY, NRPROCY, n, p, nrcppy, ibegy, iendy, miney, maxey)
  call ComputeEndpoints(MYRANKZ, NRPROCZ, n, p, nrcppz, ibegz, iendz, minez, maxez)

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
  call FillDimVector(dimensionsX, shiftsX, nrcppx, sy*sz, n, NRPROCX)
  call FillDimVector(dimensionsY, shiftsY, nrcppy, sx*sz, n, NRPROCY)
  call FillDimVector(dimensionsZ, shiftsZ, nrcppz, sx*sy, n, NRPROCZ)

  ! Compute indices for neighbours
  ibegsx = -1
  iendsx = -1
  ibegsy = -1
  iendsy = -1
  ibegsz = -1
  iendsz = -1

  do i = max(MYRANKX-1,0)+1, min(MYRANKX+1,NRPROCX-1)+1
    ix = i-MYRANKX+1
    call ComputeEndpoints(i-1, NRPROCX, n, p, nrcppx, ibegsx(ix), iendsx(ix), mine, maxe)
  enddo
  do i = max(MYRANKY-1,0)+1, min(MYRANKY+1,NRPROCY-1)+1
    iy = i-MYRANKY+1
    call ComputeEndpoints(i-1,NRPROCY, n, p, nrcppy, ibegsy(iy), iendsy(iy), mine, maxe)
  enddo
  do i = max(MYRANKZ-1,0)+1, min(MYRANKZ+1,NRPROCZ-1)+1
    iz = i-MYRANKZ+1
    call ComputeEndpoints(i-1,NRPROCZ, n, p, nrcppz, ibegsz(iz), iendsz(iz), mine, maxe)
  enddo

end subroutine


! -------------------------------------------------------------------
! Allocates most of the 'static' arrays
! -------------------------------------------------------------------
subroutine AllocateArrays
implicit none
integer :: ierr

  allocate(M(2*KL+KU+1,n+1))

  ! OLD: MP start with system fully generated along X
  ! allocate( F((n+1),(sy)*(sz))) !x,y,z
  allocate( F(sx,sy*sz)) !x,y,z
  allocate(F2(sy,sx*sz)) !y,x,z
  allocate(F3(sz,sx*sy)) !z,x,y

  allocate(Kqvals(p+1,p+1,p+1,maxex-minex+1,maxey-miney+1,maxez-minez+1))

  ! Processes on the border need pivot vector for LAPACK call
  if (MYRANKX == 0 .or. MYRANKY == 0 .or. MYRANKZ == 0) then
    allocate(IPIV(n+1))
  endif

  allocate(R(nrcppz*nrcppx*nrcppy,3,3,3))

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (iinfo == 1) then
    call stop_clock(dtime_i2,iclock_i2)
    write(*,*)'allocations:',dtime_i2
    call start_clock(iclock_i3)
  endif

end subroutine


! -------------------------------------------------------------------
! Allocates and fills the knot vector
! -------------------------------------------------------------------
subroutine PrepareKnot
implicit none

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
subroutine ComputeMassMatrix
implicit none
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
      (U,p,n,minex,maxex,nelem,                     &
       U,p,n,miney,maxey,nelem,                     &
       U,p,n,minez,maxez,nelem,                     &
       Kqvals)

end subroutine


function neighbour(dx, dy, dz) result(idx)
implicit none
integer, intent(in) :: dx, dy, dz
integer :: idx
integer :: ix, iy, iz

  ix = MYRANKX + dx + 1
  iy = MYRANKY + dy + 1
  iz = MYRANKZ + dz + 1
  idx = processors(ix, iy, iz)

end function


subroutine send_piece(items, dst, req)
implicit none
real (kind=8) :: items(:)
integer :: dst, req
integer :: ierr

  call mpi_isend(items, nrcppz*nrcppx*nrcppy, &
    MPI_DOUBLE_PRECISION, dst, 0, MPI_COMM_WORLD, req, ierr)

end subroutine


subroutine recv_piece(items, src, req)
implicit none
real (kind=8) :: items(:)
integer :: src, req
integer :: ierr

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
real (kind=8) :: spline(:,:,:,:)
integer :: i, j, k, s
integer :: request(3*3*3*2), stat(MPI_STATUS_SIZE)
integer :: ierr(3*3*3*2)
integer :: dst, src

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
integer :: i, j, k
integer :: obegx,oendx,obegy,oendy,obegz,oendz
integer :: mine, maxe

  write(*,*)PRINTRANK,'R:'

  do i = max(MYRANKX-1,0)+1, min(MYRANKX+1,NRPROCX-1)+1
    do j = max(MYRANKY-1,0)+1, min(MYRANKY+1,NRPROCY-1)+1
      do k = max(MYRANKZ-1,0)+1, min(MYRANKZ+1,NRPROCZ-1)+1
        write(*,*)'(i,j,k)=',i+1,j+1,k+1

        call ComputeEndpoints(i-1,NRPROCX,n,p,nrcppx,obegx,oendx,mine,maxe)
        call ComputeEndpoints(j-1,NRPROCY,n,p,nrcppy,obegy,oendy,mine,maxe)
        call ComputeEndpoints(k-1,NRPROCZ,n,p,nrcppz,obegz,oendz,mine,maxe)

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
real (kind=8) :: t
integer :: iter, i
integer :: ierr
real (kind=8) :: l2norm, fullnorm

  call Form3DRHS                                    &
      (U,p,n,nelem,nrcppx,                          &
       U,p,n,nelem,nrcppy,                          &
       U,p,n,nelem,nrcppz,                          &
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

  if (iinfo == 1) then
    call stop_clock(dtime_i4,iclock_i4)
    write(*,*)'Form 3D RHS:',dtime_i4
  endif

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
subroutine SolveOneDirection(RHS, eqnum)
implicit none
real (kind=8) :: RHS(:,:)
integer :: eqnum
integer :: i, iret

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
integer :: iter
real (kind=8) :: t
integer :: i
integer :: iret, ierr

  ! generate the RHS vectors
  call ComputeRHS(iter, t)

  !--------------------------------------------------------------------
  ! Solve the first problem
  !--------------------------------------------------------------------
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iprint == 1) then
    call stop_clock(dtime_init,iclock_init)
    write(*,*)dtime_init
    call start_clock(iclock_gather1)
  endif

  if (iinfo == 1) write(*,*)PRINTRANK,'1a) GATHER'

  allocate(F_out((n+1),sy*sz))

  call Gather(F,F_out,n,sx,sy*sz,dimensionsX,shiftsX,COMMX,ierr)

  if (iprint == 1) then
    write(*,*)PRINTRANK,'after call mpi_gather'
    write(*,*)PRINTRANK,'ierr',ierr
    write(*,*)PRINTRANK,'F_out:'
    do i=1,n+1
      write(*,*)PRINTRANK,i,'row=',F_out(i,1:sy*sz)
    enddo
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iprint == 1) then
    call stop_clock(dtime_gather1,iclock_gather1)
    write(*,*)dtime_gather1
    call start_clock(iclock_solve1)
  endif

  if (MYRANKX == 0) then
    if (iinfo == 1) write(*,*)PRINTRANK,'1b) SOLVE THE FIRST PROBLEM'

    call ComputeMassMatrix
    call SolveOneDirection(F_out, sy*sz)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iprint == 1) then
    call stop_clock(dtime_solve1,iclock_solve1)
    write(*,*)dtime_solve1
    call start_clock(iclock_scatter1)
  endif

  if (iinfo == 1) write(*,*)PRINTRANK,'1c) SCATTER'
  allocate(F2_out(sx,sy*sz))
  call Scatter(F_out,F2_out,n,sx,sy*sz,dimensionsX,shiftsX,COMMX,ierr)
  deallocate(F_out)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iprint == 1) then
    call stop_clock(dtime_scatter1,iclock_scatter1)
    write(*,*)dtime_scatter1
    call start_clock(iclock_gather2)
  endif

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

  allocate(F2_out((n+1),sx*sz))
  call Gather(F2,F2_out,n,sy,sx*sz,dimensionsY,shiftsY,COMMY,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iprint == 1) then
    call stop_clock(dtime_gather2,iclock_gather2)
    write(*,*)dtime_gather2
    call start_clock(iclock_solve2)
  endif


  if (MYRANKY == 0) then
    if (iinfo == 1) write(*,*)PRINTRANK,'2b) SOLVE THE SECOND PROBLEM'

    call ComputeMassMatrix
    call SolveOneDirection(F2_out, sx*sz)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iprint == 1) then
    call stop_clock(dtime_solve2,iclock_solve2)
    write(*,*)dtime_solve2
    call start_clock(iclock_scatter2)
  endif

  if (iinfo == 1) write(*,*)PRINTRANK,'2c) SCATHER'

  ! CORRECTION
  allocate(F3_out(sy,sx*sz))
  call Scatter(F2_out,F3_out,n,sy,sx*sz,dimensionsY,shiftsY,COMMY,ierr)
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

  if (iprint == 1) then
    call stop_clock(dtime_scatter2,iclock_scatter2)
    write(*,*)dtime_scatter2
    call start_clock(iclock_gather3)
  endif

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
  allocate(F3_out((n+1),sx*sy))
  call Gather(F3,F3_out,n,sz,sx*sy,dimensionsZ,shiftsZ,COMMZ,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iprint == 1) then
    call stop_clock(dtime_gather3,iclock_gather3)
    write(*,*)dtime_gather3
    call start_clock(iclock_solve3)
  endif

  if (MYRANKZ == 0) then
    if (iinfo == 1) write(*,*)PRINTRANK,'3b) SOLVE THE THIRD PROBLEM'

    call ComputeMassMatrix
    call SolveOneDirection(F3_out, sx*sy)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iprint == 1) then
    call stop_clock(dtime_solve3,iclock_solve3)
    write(*,*)dtime_solve3
    call start_clock(iclock_scatter3)
  endif

  if (iinfo == 1) write(*,*)PRINTRANK,'3c) SCATTER'

  ! CORRECTION
  allocate(F_out(sz,sx*sy))
  call Scatter(F3_out,F_out,n,sz,sx*sy,dimensionsZ,shiftsZ,COMMZ,ierr)
  deallocate(F3_out)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if (iprint == 1) then
    call stop_clock(dtime_scatter3,iclock_scatter3)
    write(*,*)dtime_scatter3
  endif

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
integer, intent(in) :: x, y, z
integer :: s
integer :: sx, sy, sz
integer :: nrcpp, ibeg, iend
integer :: mine, maxe

  call ComputeEndpoints(x, NRPROCX, n, p, nrcpp, ibeg, iend, mine, maxe)
  sx = iend - ibeg + 1
  call ComputeEndpoints(y, NRPROCY, n, p, nrcpp, ibeg, iend, mine, maxe)
  sy = iend - ibeg + 1
  call ComputeEndpoints(z, NRPROCZ, n, p, nrcpp, ibeg, iend, mine, maxe)
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
integer, intent(in) :: at
real (kind=8), intent(in) :: part(:,:)
real (kind=8), intent(out), allocatable :: full(:,:,:)
real (kind=8), allocatable :: buffer(:)
integer :: recvcounts(0:NRPROCX*NRPROCY*NRPROCZ-1)
integer :: displs(0:NRPROCX*NRPROCY*NRPROCZ-1)
integer :: x, y, z
integer :: offset, size
integer :: ierr
integer :: array_size
integer :: begx, begy, begz, endx, endy, endz
integer :: mine, maxe
integer :: nrcpp
integer :: ssx, ssy, ssz
integer :: xx, yy, zz
integer :: ix, iy, iz, idx

  ! Only the root process needs buffer, but passing unallocated array
  ! is illegal in Fortran, hence we allocate it as array of size 0
  ! in other processes.
  if (MYRANK == at) then
    array_size = (n+1)*(n+1)*(n+1)
    allocate(full(0:n,0:n,0:n))
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
          call ComputeEndpoints(x, NRPROCX, n, p, nrcpp, begx, endx, mine, maxe)
          call ComputeEndpoints(y, NRPROCY, n, p, nrcpp, begy, endy, mine, maxe)
          call ComputeEndpoints(z, NRPROCZ, n, p, nrcpp, begz, endz, mine, maxe)
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
integer :: ierr

  deallocate(shiftsX)
  deallocate(shiftsY)
  deallocate(shiftsZ)
  deallocate(dimensionsX)
  deallocate(dimensionsY)
  deallocate(dimensionsZ)

  if (allocated(IPIV)) deallocate(IPIV)
  deallocate(U)
  deallocate(M)
  deallocate(F)

  call mpi_finalize(ierr)
  if (iinfo == 1) write(*,*)PRINTRANK,"Exiting..."

  if (iprint == 1) then
    call stop_clock(dtime,iclock)
    write(*,*)dtime
  endif
end subroutine


! -------------------------------------------------------------------
! Sanity-check of dimensions vector
! -------------------------------------------------------------------
subroutine ValidateDimensions
implicit none
integer :: i, k

  k = 0
  do i = 1,NRPROCX
    k = k + dimensionsX(i)
  enddo
  if (k /= (n+1)*sy*sz) then
    write(*,*)PRINTRANK,'problem with dimensionsX',dimensionsX
    write(*,*)PRINTRANK,'n+1',n+1
    write(*,*)PRINTRANK,'sy',sy
    write(*,*)PRINTRANK,'sz',sz
    write(*,*)PRINTRANK,'nrcppx',nrcppx
    stop
  endif

  k = 0
  do i = 1,NRPROCY
    k = k + dimensionsY(i)
  enddo
  if (k /= (n+1)*sx*sz) then
    write(*,*)PRINTRANK,'problem with dimensionsY',dimensionsY
    write(*,*)PRINTRANK,'n+1',n+1
    write(*,*)PRINTRANK,'sx',sx
    write(*,*)PRINTRANK,'sz',sz
    stop
  endif

  k = 0
  do i = 1,NRPROCZ
    k = k + dimensionsZ(i)
  enddo
  if (k /= (n+1)*sx*sy) then
    write(*,*)PRINTRANK,'problem with dimensionsZ',dimensionsZ
    write(*,*)PRINTRANK,'n+1',n+1
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
  write(*,*)PRINTRANK,'n+1',n+1
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
integer :: iter
real (kind=8) :: t
real (kind=8), allocatable :: solution(:,:,:)
type (PlotParams) :: params
character(len=20) :: filename

  call GatherFullSolution(0, F, solution)

  if (MYRANK == 0) then
    write(filename,'(I10)') iter
    filename = 'step' // adjustl(filename)
    ! filename = trim(filename) // '_'

    params = PlotParams(0,1,0,1,0,1,31,31,31)
    call SaveSplinePlot(trim(filename), &
      U,p,n,nelem, &
      U,p,n,nelem, &
      U,p,n,nelem, &
      ! solution, GnuPlotOutput, params)
      solution, VtkOutput, params)

    ! call SavePlot(trim(filename), ftest, GnuPlotOutput, params)
  endif

end subroutine


subroutine ComputeResults()
implicit none
real (kind=8) :: fullpollution, fulldrained
integer :: ierr

  call Contamination                                &
      (U,p,n,nelem,nrcppx,                          &
       U,p,n,nelem,nrcppy,                          &
       U,p,n,nelem,nrcppz,                          &
       ibegx,iendx,MYRANKX,NRPROCX,                 &
       ibegy,iendy,MYRANKY,NRPROCY,                 &
       ibegz,iendz,MYRANKZ,NRPROCZ,                 &
       ibegsx,iendsx,ibegsy,iendsy,ibegsz,iendsz,   &
       minex,maxex,miney,maxey,minez,maxez,         &
       R, pollution)

  call MPI_Reduce(pollution, fullpollution, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(drained, fulldrained, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (MYRANK == 0) then
    write(*,*) fulldrained
    write(*,*) fullpollution
  endif

end subroutine


end module



program main

use stuff
implicit none

! Iteration counter
integer :: iter = 0

! Current time
real (kind=8) :: t = 0

! -------------------------------------------------------------------
! Code
! -------------------------------------------------------------------

  call Initialize

  ! prepare the problem dimensions
  p = ORDER ! order
  n = SIZE  ! intervals

  if (iinfo == 1) then
    write(*,*)'p',p,'n',n,'size of U',n+p+2
  endif

  if (SIZE<NRPROCX .or. SIZE<NRPROCY .or. SIZE<NRPROCZ) then
    write(*,*)'Number of elements smaller than number of processors'
    stop
  endif

  KL = p
  KU = p

  call ComputeDecomposition

  if (idebug == 1) then
    call ValidateDimensions
  endif

  if (iprint == 1) then
    call PrintDecompositionInfo
  endif

  call AllocateArrays
  call PrepareKnot
  call InitInputData
  call PrecomputeKq

  ! Iterations
  do iter = 0,steps

    ! if (MYRANK == 0) then
    !   write(*,*)'Iteration',iter,'/',steps
    !   write(*,*)'t = ',t
    ! endif

    call Step(iter, t)
    t = t + Dt

    ! if (mod(iter, 100) == 0) then
    !   call PrintSolution(iter, t)
    ! endif
  enddo

  call ComputeResults
  call Cleanup

end

