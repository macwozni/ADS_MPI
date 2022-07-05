module ADSS

contains

! -------------------------------------------------------------------
! Initialization of clocks and MPI
! -------------------------------------------------------------------
   subroutine initialize(nelem, p_test, p_trial, continuity, ads_test, ads_trial, ads_data, mierr)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use knot_vector, ONLY: PrepareKnot
      use basis, ONLY: BasisData
      use mpi
      implicit none
      integer(kind=4), intent(in), dimension(3) :: nelem
      integer(kind=4), intent(in), dimension(3) :: p_test, p_trial
      integer(kind=4), intent(in), dimension(3) :: continuity
      type(ADS_setup), intent(out) :: ads_trial
      type(ADS_setup), intent(out) :: ads_test
      type(ADS_compute_data), intent(out) :: ads_data
      integer(kind=4), intent(out) :: mierr
      integer(kind=4) :: ierr
      integer(kind=4), dimension(3) :: n1,n2
      integer(kind=4), dimension(3) :: ads_trialng, ads_testng

      n1 = nelem+p_test-1
      n2 = nelem+p_trial-1

      ads_testng = p_test+1
      ads_trialng = p_trial+1

      if (p_test(1).GT.p_trial(1)) ads_trialng(1)=p_test(1)+1
      if (p_test(2).GT.p_trial(2)) ads_trialng(2)=p_test(2)+1
      if (p_test(3).GT.p_trial(3)) ads_trialng(3)=p_test(3)+1

      call initialize_setup(n1, p_test, continuity, ads_testng, ads_test, mierr)
      call initialize_setup(n2, p_trial, continuity, ads_trialng, ads_trial, mierr)

      call AllocateADSdata(ads_test, ads_trial, ads_data)

   end subroutine initialize

! -------------------------------------------------------------------
! Initialization of clocks and MPI
! -------------------------------------------------------------------
   subroutine initialize_setup(n, p, continuity, ng, ads, mierr)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ
      use knot_vector, ONLY: PrepareKnot
      use basis, ONLY: BasisData
      use mpi
      implicit none
      integer(kind=4), intent(in), dimension(3) :: n
      integer(kind=4), intent(in), dimension(3) :: p
      integer(kind=4), intent(in), dimension(3) :: continuity
      integer(kind=4), intent(in), dimension(3) :: ng
      type(ADS_setup), intent(out) :: ads
      integer(kind=4), intent(out) :: mierr
      integer(kind=4) :: ierr
      integer(kind=4), dimension(3) :: nelem
      real(kind=8), allocatable, dimension(:) :: Ux
      real(kind=8), allocatable, dimension(:) :: Uy
      real(kind=8), allocatable, dimension(:) :: Uz

      call PrepareKnot(n(1), p(1), Ux, nelem(1))
      call PrepareKnot(n(2), p(2), Uy, nelem(2))
      call PrepareKnot(n(3), p(3), Uz, nelem(3))

      call AllocateADS(n, nelem, p, ng, ads)

      call move_alloc(Ux, ads%Ux)
      call move_alloc(Uy, ads%Uy)
      call move_alloc(Uz, ads%Uz)

      ads%p = p ! order
      ads%n = n ! intervals

      call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
      write (*, *) PRINTRANK, 'INITIALIZATION'
      write (*, *) 'px', p1(1), 'py', p1(2), 'pz', p1(3), &
         'nx', n(1), 'ny', n(2), 'nz', n(3), &
         'size of Ux', n(1) + p1(1) + 2, 'size of Uy', n(2) + p1(2) + 2, 'size of Uz', n(3) + p1(3) + 2
      write (*, *) 'px', p2(1), 'py', p2(2), 'pz', p2(3), &
         'nx', n(1), 'ny', n(2), 'nz', n(3), &
         'size of Ux', n(1) + p2(1) + 2, 'size of Uy', n(2) + p2(2) + 2, 'size of Uz', n(3) + p2(3) + 2
#endif

      if (n(1) < NRPROCX .or. n(2) < NRPROCY .or. n(3) < NRPROCZ) then
         write (*, *) 'Number of elements smaller than number of processors'
         stop
      end if

      call ComputeDecomposition(ads)

#ifdef IDEBUG
      call ValidateDimensions( &
         ads%n, &
         ads%s, &
         ads%nrcpp, &
         ads%dimensionsX, ads%dimensionsY, ads%dimensionsZ)
#endif

#ifdef IPRINT
      call PrintDecompositionInfo( &
         ads%n, &
         ads%nrcpp, &
         ads%ibeg, &
         ads%iend)
#endif

      ads%nelem = nelem
      mierr = 0

      ads%m = ads%n+ads%p+1

      ads%ng = ng

      call BasisData(ads%p(1), ads%m(1), ads%Ux, 1, ads%ng(1), &
                     ads%nelem(1), ads%Ox, ads%Jx, ads%Wx, ads%Xx, ads%NNx)
      call BasisData(ads%p(2), ads%m(2), ads%Uy, 1, ads%ng(2), &
                     ads%nelem(2), ads%Oy, ads%Jy, ads%Wy, ads%Xy, ads%NNy)
      call BasisData(ads%p(3), ads%m(3), ads%Uz, 1, ads%ng(3), &
                     ads%nelem(3), ads%Oz, ads%Jz, ads%Wz, ads%Xz, ads%NNz)

      ads%lnelem = ads%maxe - ads%mine + 1

#ifdef IPRINT
      write (*, *) PRINTRANK, 'ex:', ads%mine(1), ads%maxe(1)
      write (*, *) PRINTRANK, 'ey:', ads%mine(2), ads%maxe(2)
      write (*, *) PRINTRANK, 'ez:', ads%mine(3), ads%maxe(3)
      write (*, *) PRINTRANK, 'ibegx,iendx', ads%ibeg(1), ads%iend(1)
      write (*, *) PRINTRANK, 'ibegy,iendy', ads%ibeg(2), ads%iend(2)
      write (*, *) PRINTRANK, 'ibegz,iendz', ads%ibeg(3), ads%iend(3)
#endif

   end subroutine initialize_setup

! -------------------------------------------------------------------
! Establishes decomposition of the domain. Calculates size and location
! of the piece for current process.
! -------------------------------------------------------------------
   subroutine ComputeDecomposition(ads)
      use Setup, ONLY: ADS_Setup
      use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, &
                             NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints, FillDimVector
      implicit none
      type(ADS_setup), intent(inout) :: ads
      integer(kind=4) :: i
      integer(kind=4) :: ix, iy, iz
      integer(kind=4) :: imine, imaxe

      ! number of columns per processors
      call ComputeEndpoints(MYRANKX, NRPROCX, ads%n(1), ads%p(1), ads%nrcpp(1), ads%ibeg(1), &
                            ads%iend(1), ads%mine(1), ads%maxe(1))
      call ComputeEndpoints(MYRANKY, NRPROCY, ads%n(2), ads%p(2), ads%nrcpp(2), ads%ibeg(2), &
                            ads%iend(2), ads%mine(2), ads%maxe(2))
      call ComputeEndpoints(MYRANKZ, NRPROCZ, ads%n(3), ads%p(3), ads%nrcpp(3), ads%ibeg(3), &
                            ads%iend(3), ads%mine(3), ads%maxe(3))

      ads%s(1) = ads%iend(1) - ads%ibeg(1) + 1
      ads%s(2) = ads%iend(2) - ads%ibeg(2) + 1
      ads%s(3) = ads%iend(3) - ads%ibeg(3) + 1

#ifdef IINFO
      write (*, *) PRINTRANK, 'Number of cols per processor:', ads%nrcpp(1), ads%nrcpp(2), ads%nrcpp(3)
      write (*, *) PRINTRANK, 'ibegx,iendx', ads%ibeg(1), ads%iend(1)
      write (*, *) PRINTRANK, 'ibegy,iendy', ads%ibeg(2), ads%iend(2)
      write (*, *) PRINTRANK, 'ibegz,iendz', ads%ibeg(3), ads%iend(3)
#endif

      ! prepare dimensions vectors
      call FillDimVector(ads%dimensionsX, ads%shiftsX, ads%nrcpp(1), ads%s(2)*ads%s(3), ads%n(1), NRPROCX)
      call FillDimVector(ads%dimensionsY, ads%shiftsY, ads%nrcpp(2), ads%s(1)*ads%s(3), ads%n(2), NRPROCY)
      call FillDimVector(ads%dimensionsZ, ads%shiftsZ, ads%nrcpp(3), ads%s(1)*ads%s(2), ads%n(3), NRPROCZ)

      ! Compute indices for neighbours
      ads%ibegsx = -1
      ads%iendsx = -1
      ads%ibegsy = -1
      ads%iendsy = -1
      ads%ibegsz = -1
      ads%iendsz = -1

      do i = max(MYRANKX - 1, 0) + 1, min(MYRANKX + 1, NRPROCX - 1) + 1
         ix = i - MYRANKX + 1
         call ComputeEndpoints(i - 1, NRPROCX, ads%n(1), ads%p(1), ads%nrcpp(1), ads%ibegsx(ix), &
                               ads%iendsx(ix), imine, imaxe)
      end do
      do i = max(MYRANKY - 1, 0) + 1, min(MYRANKY + 1, NRPROCY - 1) + 1
         iy = i - MYRANKY + 1
         call ComputeEndpoints(i - 1, NRPROCY, ads%n(2), ads%p(2), ads%nrcpp(2), ads%ibegsy(iy), &
                               ads%iendsy(iy), imine, imaxe)
      end do
      do i = max(MYRANKZ - 1, 0) + 1, min(MYRANKZ + 1, NRPROCZ - 1) + 1
         iz = i - MYRANKZ + 1
         call ComputeEndpoints(i - 1, NRPROCZ, ads%n(3), ads%p(3), ads%nrcpp(3), ads%ibegsz(iz), &
                               ads%iendsz(iz), imine, imaxe)
      end do

   end subroutine ComputeDecomposition

! -------------------------------------------------------------------
! Allocates most of the 'static' arrays
! -------------------------------------------------------------------
   subroutine AllocateADSdata(ads_test, ads_trial, ads_data)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      ! use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ
      use mpi
      implicit none
      type(ADS_setup), intent(in) :: ads_test, ads_trial
      type(ADS_compute_data), intent(out) :: ads_data
      integer :: ierr

      allocate (ads_data%Un(ads_trial%lnelem(1), ads_trial%lnelem(2), ads_trial%lnelem(3), &
                            ads_trial%ng(1), ads_trial%ng(2), ads_trial%ng(3)))
      allocate (ads_data%Un13(ads_trial%lnelem(1), ads_trial%lnelem(2), ads_trial%lnelem(3), &
                              ads_trial%ng(1), ads_trial%ng(2), ads_trial%ng(3)))
      allocate (ads_data%Un23(ads_trial%lnelem(1), ads_trial%lnelem(2), ads_trial%lnelem(3), &
                              ads_trial%ng(1), ads_trial%ng(2), ads_trial%ng(3)))
      allocate (ads_data%dUn(ads_trial%lnelem(1), ads_trial%lnelem(2), ads_trial%lnelem(3), &
                             ads_trial%ng(1), ads_trial%ng(2), ads_trial%ng(3), 3))

      ! OLD: MP start with system fully generated along X
      ! allocate( F((n+1),(sy)*(sz))) !x,y,z
      !allocate( ads_data % F_test(ads % s(1), ads % s(2) * ads % s(3))) !x,y,z


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO CHANGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate (ads_data%R(ads_trial%nrcpp(3)*ads_trial%nrcpp(1)*ads_trial%nrcpp(2), 3, 3, 3))
      ads_data%R = 0.d0

      call mpi_barrier(MPI_COMM_WORLD, ierr)

   end subroutine AllocateADSdata

! -------------------------------------------------------------------
! Allocates most of the 'static' arrays
! -------------------------------------------------------------------
   subroutine AllocateADS(n, nelem, p, ng, ads)
      use Setup, ONLY: ADS_Setup
      use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ
      use mpi
      implicit none
      integer(kind=4), dimension(3), intent(in) :: n, nelem, p, ng
      type(ADS_setup), intent(out) :: ads
      ! integer :: ierr

      allocate (ads%Ox(nelem(1)))
      allocate (ads%Oy(nelem(2)))
      allocate (ads%Oz(nelem(3)))

      allocate (ads%Jx(nelem(1)))
      allocate (ads%Jy(nelem(2)))
      allocate (ads%Jz(nelem(3)))

      allocate (ads%Xx(ng(1), nelem(1)))
      allocate (ads%Xy(ng(2), nelem(2)))
      allocate (ads%Xz(ng(3), nelem(3)))

      allocate (ads%NNx(0:1, 0:p(1), ng(1), nelem(1)))
      allocate (ads%NNy(0:1, 0:p(2), ng(2), nelem(2)))
      allocate (ads%NNz(0:1, 0:p(3), ng(3), nelem(3)))

      allocate (ads%Wx(ng(1)))
      allocate (ads%Wy(ng(2)))
      allocate (ads%Wz(ng(3)))

      ! Processes on the border need pivot vector for LAPACK call
      if (MYRANKX == 0 .or. MYRANKY == 0 .or. MYRANKZ == 0) then
         allocate (ads%IPIVx(n(1) + 1))
         allocate (ads%IPIVy(n(2) + 1))
         allocate (ads%IPIVz(n(3) + 1))
      end if
   end subroutine AllocateADS

!!!! przeniesc do solver
! -------------------------------------------------------------------
! Solves 1D linear system (one of 3 steps of solving the whole),
! using DGBSV.
!
! RHS   - vector of right-hand sides, of dimension (n+1) x eqnum
! eqnum - number of right-hand sides (equations)
!
!
!!!!!!!!!!!!!
! Transforms all data to MUMPS format and calls MUMPS solver
!
! -------------------------------------------------------------------
   subroutine SolveOneDirection(RHS, eqnum, n, p, sprsmtrx)
      use sparse
      use mpi
      implicit none
      include 'dmumps_struc.h'
      real(kind=8) :: RHS(:, :)
      integer :: n, p
      integer(kind=4) :: eqnum
      integer(kind=4) :: i!, iret
      type(sparse_matrix), pointer, intent(in) :: sprsmtrx
      type(dmumps_struc) :: mumps_par

!  initialize MUMPS
      mumps_par%comm = MPI_COMM_SELF
      mumps_par%job = -1
      mumps_par%par = 1
      mumps_par%N = n + 1
      call dmumps(mumps_par)

!  convert LHS and RHS to MUMPS format
      call to_mumps_format(sprsmtrx, mumps_par)
      allocate (mumps_par%rhs(mumps_par%n))

!  set MUMPS parameters
!  error output stream (non-positive to suppress)
      mumps_par%icntl(1) = 1 !1
!  diagnostic, statistics and warnings
      mumps_par%icntl(2) = 0! 1 !1
!  global information
      mumps_par%icntl(3) = 0!6 !6
!  printing level
      mumps_par%icntl(4) = 0!3 !3
!  input matrix in assembled or element format
      mumps_par%icntl(5) = 0
!  column permutation for zero-free diagonal (automatic)
!  mumps_par%icntl(6)  = 7
!  pivot order (automatic)
      mumps_par%icntl(7) = 4 !1 enforce ordering, 5 metis, 0 AMD, 7 auto
!  scaling (automatic)
!  mumps_par%icntl(8)  = 7
!  no transpose
!  mumps_par%icntl(9)  = 1
!  max steps for iterative refinement
!  mumps_par%icntl(10) = 0
!  statistics info
      mumps_par%icntl(11) = 2
!  controls parallelism
      mumps_par%icntl(12) = 0
!  use ScaLAPACK for root node
      mumps_par%icntl(13) = 1 !0 use 1 do not use
!  percentage increase in estimated workspace
      mumps_par%icntl(14) = 50
!  matrix distribution for assembled input
      mumps_par%icntl(18) = 0 !distributed
!  nonzero for Schur complement
      mumps_par%icntl(19) = 0
!  distribution of RHS (centralized on host)
      mumps_par%icntl(20) = 0
!  mumps_par%icntl(32) = 1

!  start MUMPS
      mumps_par%job = 1
      call dmumps(mumps_par)

      mumps_par%job = 2
      call dmumps(mumps_par)
      if (mumps_par%info(1) .ne. 0) then
         write (*, *) 'mumps_par%job=', mumps_par%job
         write (*, *) 'mumps_par%info=', mumps_par%info
         stop 1
      end if

      do i = 1, eqnum
         mumps_par%rhs(1:n + 1) = rhs(1:n + 1, i)
         mumps_par%job = 3
         call dmumps(mumps_par)
         rhs(1:n + 1, i) = mumps_par%rhs(1:n + 1)
      end do

!  stop MUMPS
      mumps_par%job = -2
      call dmumps(mumps_par)

      deallocate (mumps_par%irn)
      deallocate (mumps_par%jcn)
      deallocate (mumps_par%a)
      deallocate (mumps_par%rhs)

   end subroutine SolveOneDirection

! -------------------------------------------------------------------
! Performs one step of the simulation with multiple substeps
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
   subroutine MultiStep(iter, mix, RHS_fun, ads_test, ads_trial, ads_data, n, alpha_step, mierr)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use projection_engine, ONLY: FormUn
      use Interfaces, ONLY: forcing_fun
      implicit none
      integer(kind=4), intent(in) :: iter
      real(kind=8), intent(in) :: mix(4, 3)
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in), dimension(7, 3) :: alpha_step
      procedure(forcing_fun) :: RHS_fun
      type(ADS_setup), intent(in) :: ads_test, ads_trial
      type(ADS_compute_data), intent(inout) :: ads_data
      integer(kind=4), intent(out) :: mierr
      real(kind=8) :: mmix(4)
      integer(kind=4), dimension(3) :: direction
      integer(kind=4) :: substep
      integer(kind=4), dimension(3, 3) :: abc


      if (allocated(ads_data%R)) deallocate(ads_data%R)
      allocate (ads_data%R(ads_trial%nrcpp(3)*ads_trial%nrcpp(1)*ads_trial%nrcpp(2), 3, 3, 3))
      ads_data%R = 0.d0

      allocate (ads_data%F (ads_trial%s(1), ads_trial%s(2)*ads_trial%s(3))) !x,y,z
      allocate (ads_data%F2(ads_trial%s(2), ads_trial%s(3)*ads_trial%s(1))) !y,z,x
      allocate (ads_data%F3(ads_trial%s(3), ads_trial%s(1)*ads_trial%s(2))) !z,x,y
      allocate (ads_data%Ft (ads_test%s(1), ads_trial%s(2)*ads_trial%s(3))) !x,y,z
      allocate (ads_data%Ft2(ads_trial%s(2), ads_trial%s(3)*ads_test%s(1))) !y,z,x
      allocate (ads_data%Ft3(ads_trial%s(3), ads_test%s(1)*ads_trial%s(2))) !z,x,y

      mmix = mix(:, 1)
      direction = (/1, 0, 0/) ! x
      abc(:, 1) = (/1, 2, 3/) ! x y z
      abc(:, 2) = (/2, 3, 1/) ! y z x
      abc(:, 3) = (/3, 1, 2/) ! z x y
      substep = 1
      call FormUn(substep, ads_trial, ads_data)
      call Sub_Step(ads_test, ads_trial, iter, mmix, direction, substep, abc, &
                    n, alpha_step, RHS_fun, ads_data, mierr)
      if (allocated(ads_data%Ft)) deallocate(ads_data%Ft)
      if (allocated(ads_data%Ft2)) deallocate(ads_data%Ft2)
      if (allocated(ads_data%Ft3)) deallocate(ads_data%Ft3)
      if (allocated(ads_data%F)) deallocate(ads_data%F)
      if (allocated(ads_data%F2)) deallocate(ads_data%F2)
      if (allocated(ads_data%F3)) deallocate(ads_data%F3)

      if (allocated(ads_data%R)) deallocate(ads_data%R)
      allocate (ads_data%R(ads_trial%nrcpp(2)*ads_trial%nrcpp(3)*ads_trial%nrcpp(1), 3, 3, 3))
      ads_data%R = 0.d0

      allocate (ads_data%F (ads_trial%s(2), ads_trial%s(3)*ads_trial%s(1))) !y,z,x
      allocate (ads_data%F2(ads_trial%s(3), ads_trial%s(1)*ads_trial%s(2))) !z,x,y
      allocate (ads_data%F3(ads_trial%s(1), ads_trial%s(2)*ads_trial%s(3))) !x,y,z
      allocate (ads_data%Ft (ads_test%s(2), ads_trial%s(3)*ads_trial%s(1))) !y,z,x
      allocate (ads_data%Ft2(ads_trial%s(3), ads_trial%s(1)*ads_test%s(2))) !z.x.y
      allocate (ads_data%Ft3(ads_trial%s(1), ads_test%s(2)*ads_trial%s(3))) !x,y,z

      mmix = mix(:, 2)
      direction = (/0, 1, 0/) ! y
      abc(:, 1) = (/2, 3, 1/) ! y z x
      abc(:, 2) = (/3, 1, 2/) ! z x y
      abc(:, 3) = (/1, 2, 3/) ! x y z
      substep = 2
      call FormUn(substep, ads_trial, ads_data)
      call Sub_Step(ads_test, ads_trial, iter, mmix, direction, substep, abc, &
                    n, alpha_step, RHS_fun, ads_data, mierr)
      if (allocated(ads_data%Ft)) deallocate(ads_data%Ft)
      if (allocated(ads_data%Ft2)) deallocate(ads_data%Ft2)
      if (allocated(ads_data%Ft3)) deallocate(ads_data%Ft3)
      if (allocated(ads_data%F)) deallocate(ads_data%F)
      if (allocated(ads_data%F2)) deallocate(ads_data%F2)
      if (allocated(ads_data%F3)) deallocate(ads_data%F3)

      if (allocated(ads_data%R)) deallocate(ads_data%R)
      allocate (ads_data%R(ads_trial%nrcpp(1)*ads_trial%nrcpp(2)*ads_trial%nrcpp(3), 3, 3, 3))
      ads_data%R = 0.d0

      allocate (ads_data%F (ads_trial%s(3), ads_trial%s(1)*ads_trial%s(2))) !z,x,y
      allocate (ads_data%F2(ads_trial%s(1), ads_trial%s(2)*ads_trial%s(3))) !x,y,z
      allocate (ads_data%F3(ads_trial%s(2), ads_trial%s(3)*ads_trial%s(1))) !y,z,x
      allocate (ads_data%Ft (ads_test%s(3), ads_trial%s(1)*ads_trial%s(2))) !z,x,y
      allocate (ads_data%Ft2(ads_trial%s(1), ads_trial%s(2)*ads_test%s(3))) !x,y,z
      allocate (ads_data%Ft3(ads_trial%s(2), ads_test%s(3)*ads_trial%s(1))) !y,z,x

      mmix = mix(:, 3)
      direction = (/0, 0, 1/) ! z
      abc(:, 1) = (/3, 1, 2/) ! z y y
      abc(:, 2) = (/1, 2, 3/) ! x y z
      abc(:, 3) = (/2, 3, 1/) ! y z x
      substep = 3
      call FormUn(substep, ads_trial, ads_data)
      call Sub_Step(ads_test, ads_trial, iter, mmix, direction, substep, abc, &
                    n, alpha_step, RHS_fun, ads_data, mierr)
      if (allocated(ads_data%Ft)) deallocate(ads_data%Ft)
      if (allocated(ads_data%Ft2)) deallocate(ads_data%Ft2)
      if (allocated(ads_data%Ft3)) deallocate(ads_data%Ft3)
      if (allocated(ads_data%F)) deallocate(ads_data%F)
      if (allocated(ads_data%F2)) deallocate(ads_data%F2)
      if (allocated(ads_data%F3)) deallocate(ads_data%F3)

   end subroutine MultiStep

! -------------------------------------------------------------------
! Performs one step of the simulation with single substep
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
   subroutine Step(iter, RHS_fun, ads, ads_data, mierr)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use projection_engine, ONLY: FormUn
      use Interfaces, ONLY: forcing_fun
      implicit none
      integer(kind=4), intent(in) :: iter
      procedure(forcing_fun) :: RHS_fun
      type(ADS_setup), intent(in) :: ads
      type(ADS_compute_data), intent(inout) :: ads_data
      integer(kind=4), intent(out) :: mierr
      real(kind=8) :: mix(4)
      integer(kind=4), dimension(3) :: direction
      integer(kind=4) :: substep
      integer(kind=4), dimension(3, 3) :: abc
      real(kind=8), dimension(7, 3) :: alpha_step

      mix = (/1.d0, 0.d0, 0.d0, 0.d0/)
      direction = (/0, 0, 0/)
      abc(:, 1) = (/1, 2, 3/) ! x y z
      abc(:, 2) = (/2, 1, 3/) ! y x z
      abc(:, 3) = (/3, 1, 2/) ! z x y
      substep = 1
      ads_data%un13 = 0.d0
      ads_data%un23 = 0.d0
      call FormUn(1, ads, ads_data)
      alpha_step = 1.d0


      if (allocated(ads_data%R)) deallocate(ads_data%R)
      allocate (ads_data%R(ads%nrcpp(3)*ads%nrcpp(1)*ads%nrcpp(2), 3, 3, 3))
      ads_data%R = 0.d0

      allocate (ads_data%F (ads%s(1), ads%s(2)*ads%s(3))) !x,y,z
      allocate (ads_data%F2(ads%s(2), ads%s(3)*ads%s(1))) !y,x,z
      allocate (ads_data%F3(ads%s(3), ads%s(1)*ads%s(2))) !z,x,y

      !call Sub_Step(ads, ads, iter, mix,direction,substep,abc,RHS_fun,ads_data, mierr)
      call Sub_Step(ads, ads, iter, mix, direction, substep, abc, &
                    1, alpha_step, RHS_fun, ads_data, mierr)
      if (allocated(ads_data%F)) deallocate(ads_data%F)
      if (allocated(ads_data%F2)) deallocate(ads_data%F2)
      if (allocated(ads_data%F3)) deallocate(ads_data%F3)

   end subroutine Step

!!!! podzielic na wraper i czesc wlasciwa
! przeniesc czesc do solver
!---------------------------------------------------------------------------
!> @author Maciej Wozniak
!>
!> @brief
!> Calculates single substep for different time integration schemes.
!>
!
! Input:
! ------
!> @param[in] ads_test      - knot vector
!> @param[in] ads_trial
!> @param[in] iter         - number of the iteration
!> @param[in] mix          - mixing values for MKBBT matrices
!> @param[in] directon     - direction in which we enrich test space
!> @param[in] substep      - number of substep
!> @param[in] abc
!> @param[in] RHS_fun
!
! Input/Output:
! -------
!> @param[inout] ads_data  - data structures for ADS
!
! Output:
! -------
!> @param[out] mierr
!
! -------------------------------------------------------------------
   subroutine Sub_Step(ads_test, ads_trial, iter, mix, direction, substep, abc, &
                       n, alpha_step, &
                       RHS_fun, ads_data, mierr)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      ! use parallelism, ONLY: PRINTRANK, MYRANKX, MYRANKY, MYRANKZ
      ! use communicators, ONLY: COMMX, COMMY, COMMZ
      use reorderRHS, ONLY: ReorderRHSForX, ReorderRHSForY, ReorderRHSForZ
      use projection_engine, ONLY: Form3DRHS, ComputeMatrix
      use my_mpi, ONLY: DistributeSpline, Gather, Scatter
      use Interfaces, ONLY: forcing_fun
      use mpi
      use sparse
      implicit none
      type(ADS_setup), intent(in) :: ads_test
      type(ADS_setup), intent(in) :: ads_trial
      integer(kind=4), intent(in) :: iter
      real(kind=8), intent(in) :: mix(4)
      integer(kind=4), intent(in), dimension(3) :: direction
      integer(kind=4), intent(in) :: substep
      integer(kind=4), dimension(3, 3), intent(in) :: abc
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in), dimension(7, 3) :: alpha_step
      procedure(forcing_fun) :: RHS_fun
      type(ADS_compute_data), intent(inout) :: ads_data
      integer(kind=4), intent(out) :: mierr
      integer(kind=4) :: i,a,b,c
      integer(kind=4) :: ierr!, iret
      integer(kind=4), dimension(3) :: nrcpp
      ! real(kind=8) :: time1, time2
      logical :: igrm

#ifdef PERFORMANCE
      time1 = MPI_Wtime()
#endif

      if (allocated(ads_data%F)) ads_data%F = 0.d0
      if (allocated(ads_data%Ft)) ads_data%Ft = 0.d0
! generate the RHS vectors
      call Form3DRHS(ads_test, ads_trial, ads_data, direction, n, substep, alpha_step, RHS_fun, igrm)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write (*, *) "Form 3D RHS: ", time2 - time1
#endif

#ifdef IPRINT
      write (*, *) PRINTRANK, 'F'
      do i = 1, ads_trial%s(1)
         write (*, *) PRINTRANK, ads_trial%F(i, 1:ads_trial%s(2)*ads_trial%s(3))
      end do
#endif

      a=abc(1, 1)
      b=abc(2, 1)
      c=abc(3, 1)

!--------------------------------------------------------------------
! Solve the first problem
!--------------------------------------------------------------------
      call solve_problem(ads_test, ads_trial, abc(1, 1), abc(2, 1), abc(3, 1), &
                         mix, mix, mix, direction, igrm, ads_data%F, ads_data%F2, ads_data%Ft, ads_data%Ft2, ierr)

!--------------------------------------------------------------------
! Solve the second problem
!--------------------------------------------------------------------
      call solve_problem(ads_test, ads_trial, abc(1, 2), abc(2, 2), abc(3, 2), &
                         mix, mix, mix, direction, igrm, ads_data%F2, ads_data%F3, ads_data%Ft2, ads_data%Ft3, ierr)

!--------------------------------------------------------------------
! Solve the third problem
!--------------------------------------------------------------------
      call solve_problem(ads_test, ads_trial, abc(1, 3), abc(2, 3), abc(3, 3), &
                         mix, mix, mix, direction, igrm, ads_data%F3, ads_data%F, ads_data%Ft3, ads_data%Ft, ierr)

#ifdef IINFO
      write (*, *) PRINTRANK, '3e) DISTRIBUTE SOLUTION'
#endif
!  copy results to proper buffer
      do i = 1, ads_trial%s(b)*ads_trial%s(c)
         ads_data%R((i - 1)*ads_trial%s(a) + 1:i*ads_trial%s(a), 2, 2, 2) = ads_data%F(:, i)
      end do
!  nrcpp - number of columns (average) per processor
      nrcpp = (/ads_trial%nrcpp(c), ads_trial%nrcpp(a), ads_trial%nrcpp(b)/)
      call DistributeSpline(ads_data%R, nrcpp, ads_data%R)

#ifdef IPRINT
      write (*, *) PRINTRANK, 'Result:'
      do i = 1, ads_trial%s(3)
         write (*, *) PRINTRANK, i, 'row=', ads_trial%F(i, :)
      end do
#endif

      call mpi_barrier(MPI_COMM_WORLD, ierr)

      mierr = 0
   end subroutine Sub_Step

! -------------------------------------------------------------------
! Deallocates all the resources and finalizes MPI.
! -------------------------------------------------------------------
   subroutine Cleanup(ads_test, ads_trial, ads_data, mierr)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      ! use parallelism, ONLY: PRINTRANK
      implicit none
      type(ADS_setup), intent(inout) :: ads_test, ads_trial
      type(ADS_compute_data), intent(inout) :: ads_data
      integer(kind=4), intent(out) :: mierr
      integer(kind=4) :: ierr

      call Cleanup_ADS(ads_test, mierr)
      call Cleanup_ADS(ads_trial, mierr)

      if (allocated(ads_data%F)) deallocate (ads_data%F)
      if (allocated(ads_data%F2)) deallocate (ads_data%F2)
      if (allocated(ads_data%F3)) deallocate (ads_data%F3)

      if (allocated(ads_data%Un)) deallocate (ads_data%Un)
      if (allocated(ads_data%Un13)) deallocate (ads_data%Un13)
      if (allocated(ads_data%Un23)) deallocate (ads_data%Un23)
      if (allocated(ads_data%dUn)) deallocate (ads_data%dUn)
      !!!!!! wyciac
      call mpi_finalize(ierr)
#ifdef IINFO
      write (*, *) PRINTRANK, "Exiting..."
#endif

      mierr = 0

   end subroutine Cleanup

! -------------------------------------------------------------------
! Deallocates all the resources and finalizes MPI.
! -------------------------------------------------------------------
   subroutine Cleanup_ADS(ads, mierr)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      ! use parallelism, ONLY: PRINTRANK
      implicit none
      type(ADS_setup), intent(inout) :: ads
      integer(kind=4), intent(out) :: mierr
      ! integer(kind=4) :: ierr

      if (allocated(ads%shiftsX)) deallocate (ads%shiftsX)
      if (allocated(ads%shiftsY)) deallocate (ads%shiftsY)
      if (allocated(ads%shiftsZ)) deallocate (ads%shiftsZ)

      if (allocated(ads%dimensionsX)) deallocate (ads%dimensionsX)
      if (allocated(ads%dimensionsX)) deallocate (ads%dimensionsY)
      if (allocated(ads%dimensionsZ)) deallocate (ads%dimensionsZ)

      if (allocated(ads%IPIVx)) deallocate (ads%IPIVx)
      if (allocated(ads%IPIVy)) deallocate (ads%IPIVy)
      if (allocated(ads%IPIVz)) deallocate (ads%IPIVz)

      if (allocated(ads%Ux)) deallocate (ads%Ux)
      if (allocated(ads%Uy)) deallocate (ads%Uy)
      if (allocated(ads%Uz)) deallocate (ads%Uz)

      if (allocated(ads%Ox)) deallocate (ads%Ox)
      if (allocated(ads%Oy)) deallocate (ads%Oy)
      if (allocated(ads%Oz)) deallocate (ads%Oz)

      if (allocated(ads%Jx)) deallocate (ads%Jx)
      if (allocated(ads%Jy)) deallocate (ads%Jy)
      if (allocated(ads%Jz)) deallocate (ads%Jz)

      if (allocated(ads%Xx)) deallocate (ads%Xx)
      if (allocated(ads%Xy)) deallocate (ads%Xy)
      if (allocated(ads%Xz)) deallocate (ads%Xz)

      if (allocated(ads%NNx)) deallocate (ads%NNx)
      if (allocated(ads%NNy)) deallocate (ads%NNy)
      if (allocated(ads%NNz)) deallocate (ads%NNz)

      if (allocated(ads%Wx)) deallocate (ads%Wx)
      if (allocated(ads%Wy)) deallocate (ads%Wy)
      if (allocated(ads%Wz)) deallocate (ads%Wz)
      mierr = 0
   end subroutine Cleanup_ADS

! -------------------------------------------------------------------
! Gathers full solution and plots it
! -------------------------------------------------------------------
   subroutine PrintSolution(iter, t, ads, ads_data)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use parallelism, ONLY: MYRANK
      use plot, ONLY: SaveSplinePlot, PlotParams
      use vtk, ONLY: VtkOutput
      use my_mpi, ONLY: GatherFullSolution
      implicit none
      type(ADS_compute_data), intent(in) :: ads_data
      type(ADS_setup), intent(in) :: ads
      integer(kind=4), intent(in) :: iter
      real(kind=8), intent(in) :: t
      real(kind=8), allocatable :: solution(:, :, :)
      type(PlotParams) :: params
      character(len=20) :: filename

      call GatherFullSolution(0, ads_data%F, solution, &
                              ads%n, ads%p, ads%s)

      if (MYRANK == 0) then
         write (filename, '(I10)') iter
         filename = 'step'//adjustl(filename)
         ! filename = trim(filename) // '_'

         params = PlotParams(0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 31, 31, 31)
         call SaveSplinePlot(trim(filename), &
                             ads%Ux, ads%p(1), ads%n(1), ads%nelem(1), &
                             ads%Uy, ads%p(2), ads%n(2), ads%nelem(2), &
                             ads%Uz, ads%p(3), ads%n(3), ads%nelem(3), &
                             ! solution, GnuPlotOutput, params)
                             solution, VtkOutput, params)

         ! call SavePlot(trim(filename), ftest, GnuPlotOutput, params)
      end if

   end subroutine PrintSolution

!---------------------------------------------------------------------------
!> @author Maciej Wozniak
!>
!> @brief
!> Solves one sub-direction of single time sub-step.
!> Works both with and without iGRM.
!> Takes care about gathering and scattering result (and residuum).
!> Automagically takes care about transfering both parts of RHS.
!
! Input:
! ------
!> @param[in] ads_test  -
!> @param[in] ads_trial -
!> @param[in] a         -
!> @param[in] b         -
!> @param[in] c         -
!> @param[in] mix       - mixing proportions of M, K, B and BT matrices
!> @param[in] direction -
!> @param[in] igrm      -
!> @param[in] mixBT     -
!
! Output:
! -------
!> @param[out] F        -
!> @param[out] F2       -
!> @param[out] ierr     -
!
! -------------------------------------------------------------------
   subroutine solve_problem(ads_test, ads_trial, a, b, c, mixA, mixB, mixBT, direction, igrm, F, F2, Ft, Ft2, ierr)
      use Setup, ONLY: ADS_Setup
      use sparse
      use mpi
      use my_mpi, ONLY: Gather, Scatter
      use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, NRPROCX, NRPROCY, NRPROCZ, FillDimVector
      use communicators, ONLY: COMMX, COMMY, COMMZ
      use reorderRHS, ONLY: ReorderRHSForX, ReorderRHSForY, ReorderRHSForZ
      use projection_engine, ONLY: ComputeMatrix
      implicit none
      type(ADS_setup), intent(in) :: ads_test, ads_trial
      integer(kind=4), intent(in) :: a, b, c
      real(kind=8), dimension(4), intent(in) :: mixA, mixB, mixBT
      integer(kind=4), dimension(3), intent(in) :: direction
      logical, intent(in) :: igrm
      real(kind=8), allocatable, dimension(:, :) :: F, F2
      real(kind=8), allocatable, dimension(:, :) :: Ft, Ft2
      integer(kind=4), intent(out) :: ierr
      integer(kind=4), dimension(:), allocatable :: dimensions_test, dimensions_trial ! Size of slices of domain in each dimension
      integer(kind=4), allocatable, dimension(:) :: shifts_test, shifts_trial
      integer(kind=4) :: comm
      real(kind=8), allocatable, dimension(:) :: U_trial, U_test
      integer(kind=4) :: myrankdim ! Integer coordinates of processor along X, Y or Z
      type(sparse_matrix), pointer :: sprsmtrx
      real(kind=8), allocatable, dimension(:, :) :: Fs ! F-solve
      real(kind=8), allocatable, dimension(:, :) :: F_out, F2_out ! F-trial
      real(kind=8), allocatable, dimension(:, :) :: Ft_out, Ft2_out ! F-test
      integer(kind=4), dimension(3) :: ibeg, iend
      ! real(kind=8) :: time1, time2
      logical :: equ

!  we have identical test and trial spaces if equ=true in given direction
      equ = .TRUE.
!  if we have enriched test space in given direction
      if (direction(a) .EQ. 1) equ = .FALSE.

!  set proper paremeters depending on which direction we solve
!  we solve in x directon
      if (a .EQ. 1) then
         comm = COMMX
         myrankdim = MYRANKX
         U_trial = ads_trial%ux
         U_test = ads_test%ux
         shifts_trial = ads_trial%shiftsX
         shifts_test = ads_test%shiftsX
         dimensions_trial = ads_trial%dimensionsX
         dimensions_test = ads_test%dimensionsX
         ! prepare dimensions vectors
         call FillDimVector(dimensions_test, shifts_test,ads_test%nrcpp(1),&
            (direction(2)*ads_test%s(2)+(1-direction(2))*ads_trial%s(2))*&
            (direction(3)*ads_test%s(3)+(1-direction(3))*ads_trial%s(3)),&
            (direction(1)*ads_test%n(1)+(1-direction(1))*ads_trial%n(1)),&
            NRPROCX)
!  we solve in y directon
      else if (a .EQ. 2) then
         comm = COMMY
         myrankdim = MYRANKY
         U_trial = ads_trial%uy
         U_test = ads_test%uy
         shifts_trial = ads_trial%shiftsY
         shifts_test = ads_test%shiftsY
         dimensions_trial = ads_trial%dimensionsY
         dimensions_test = ads_test%dimensionsY
         ! prepare dimensions vectors
         call FillDimVector(dimensions_test, shifts_test,ads_test%nrcpp(2),&
            (direction(1)*ads_test%s(1)+(1-direction(1))*ads_trial%s(1))*&
            (direction(3)*ads_test%s(3)+(1-direction(3))*ads_trial%s(3)),&
            (direction(2)*ads_test%n(2)+(1-direction(2))*ads_trial%n(2)),&
            NRPROCY)
!  we solve in z directon
      else ! a.EQ.3
         comm = COMMZ
         myrankdim = MYRANKZ
         U_trial = ads_trial%uz
         U_test = ads_test%uz
         shifts_trial = ads_trial%shiftsZ
         shifts_test = ads_test%shiftsZ
         dimensions_trial = ads_trial%dimensionsZ
         dimensions_test = ads_test%dimensionsZ
         ! prepare dimensions vectors
         call FillDimVector(dimensions_test, shifts_test,ads_test%nrcpp(3),&
            (direction(1)*ads_test%s(1)+(1-direction(1))*ads_trial%s(1))*&
            (direction(2)*ads_test%s(2)+(1-direction(2))*ads_trial%s(2)),&
            (direction(3)*ads_test%n(3)+(1-direction(3))*ads_trial%n(3)),&
            NRPROCZ)
      end if

      ibeg = ads_trial%ibeg
      iend = ads_trial%iend
      if (direction(1) .EQ. 1) then
            ibeg(1) = ads_test%ibeg(1)
            iend(1) = ads_test%iend(1)
      endif
      if (direction(2) .EQ. 1) then
            ibeg(2) = ads_test%ibeg(2)
            iend(2) = ads_test%iend(2)
      endif
      if (direction(3) .EQ. 1) then
            ibeg(3) = ads_test%ibeg(3)
            iend(3) = ads_test%iend(3)
      endif

      call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
      write (*, *) PRINTRANK, a, 'a) GATHER'
      call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

!  allocate result buffer
      allocate (Fs((ads_trial%n(a) + 1 + direction(a)*(ads_test%n(a) + 1)), &
                   (ads_trial%s(b) + direction(b)*ads_test%s(b))*(ads_trial%s(c) + direction(c)*ads_test%s(c))))
      allocate (F_out((ads_trial%n(a) + 1), &
                      (ads_trial%s(b)*ads_trial%s(c))))
#ifdef PERFORMANCE
      time1 = MPI_Wtime()
#endif
!  gather onto the face of processors
      call Gather(F, F_out, ads_trial%n(a), &
                  ads_trial%s(a), &
                  ads_trial%s(b)*ads_trial%s(c), &
                  dimensions_trial, shifts_trial, comm, ierr)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write (*, *) "Gather", a, " : ", time2 - time1
#endif

      if (igrm) then
!  allocate result buffer
         allocate (Ft_out(((1 - direction(a))*ads_trial%s(a) + direction(a)*ads_test%s(a)), &
                          ((1 - direction(b))*ads_trial%s(b) + direction(b)*ads_test%s(b))* &
                          ((1 - direction(c))*ads_trial%s(c) + direction(c)*ads_test%s(c))))
#ifdef PERFORMANCE
         time1 = MPI_Wtime()
#endif
!  gather onto the face of processors
         call Gather(Ft, Ft_out, (1 - direction(a))*ads_trial%n(a) + direction(a)*ads_test%n(a), &
                     (1 - direction(a))*ads_trial%s(a) + direction(a)*ads_test%s(a), &
                     ((1 - direction(b))*ads_trial%s(b) + direction(b)*ads_test%s(b))* &
                     ((1 - direction(c))*ads_trial%s(c) + direction(c)*ads_test%s(c)), &
                     dimensions_test, shifts_test, comm, ierr)
#ifdef PERFORMANCE
         time2 = MPI_Wtime()
         write (*, *) "Gather", a, " : ", time2 - time1
#endif
         if (equ) then
            Fs(:, 1:ads_trial%s(b)*ads_trial%s(c)) = F_out
            Fs(:, (ads_trial%s(b)*ads_trial%s(c) + 1): &
               (ads_trial%s(b) + direction(b)*ads_test%s(b))*(ads_trial%s(c) + direction(c)*ads_test%s(c))) = Ft_out
         else
            Fs(1:ads_trial%n(a) + 1, :) = F_out
            Fs(ads_trial%n(a) + 2:ads_trial%n(a) + 1 + direction(a)*(ads_test%n(a) + 1), :) = Ft_out
         end if
      else
         Fs = F_out
      end if
      call mpi_barrier(MPI_COMM_WORLD, ierr)

!  performed only on face of processors
      if (myrankdim == 0) then
#ifdef IINFO
         write (*, *) PRINTRANK, a, 'b) SOLVE THE ', a, ' PROBLEM'
#endif

#ifdef PERFORMANCE
         time1 = MPI_Wtime()
#endif
!     compute LHS matrix
         call ComputeMatrix(U_test, ads_test%p(a), ads_test%n(a), ads_test%nelem(a), &
                            U_trial, ads_trial%p(a), ads_trial%n(a), ads_trial%nelem(a), &
                            mixA, mixB, mixBT, equ, sprsmtrx)
#ifdef PERFORMANCE
         time2 = MPI_Wtime()
         write (*, *) "Mass matrix", a, ": ", time2 - time1
         time1 = MPI_Wtime()
#endif
!     perform real solver
         call SolveOneDirection(Fs, (ads_trial%s(b) + direction(b)*ads_test%s(b)) &
                                *(ads_trial%s(c) + direction(c)*ads_test%s(c)), &
                                (ads_trial%n(a) + direction(a)*ads_test%n(a)), &
                                (ads_trial%n(a) + direction(a)*ads_test%n(a)), sprsmtrx)
!     clean buffers
         call clear_matrix(sprsmtrx)
#ifdef PERFORMANCE
         time2 = MPI_Wtime()
         write (*, *) "Solve ", a, ": ", time2 - time1
#endif
      end if

      call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
      write (*, *) PRINTRANK, a, 'c) SCATTER'
#endif
      if (igrm) then
         if (equ) then
            F_out = Fs(:, 1:ads_trial%s(b)*ads_trial%s(c))
            Ft_out = Fs(:, (ads_trial%s(b)*ads_trial%s(c) + 1): &
                        (ads_trial%s(b) + direction(b)*ads_test%s(b))*(ads_trial%s(c) + direction(c)*ads_test%s(c)))
         else
            F_out = Fs(1:ads_trial%n(a) + 1, :)
            Ft_out = Fs(ads_trial%n(a) + 2:ads_trial%n(a) + 1 + direction(a)*(ads_test%n(a) + 1), :)
         end if
!  allocate buffers
         allocate (Ft2_out(((1 - direction(a))*ads_trial%s(a) + direction(a)*ads_test%s(a)), &
                           ((1 - direction(b))*ads_trial%s(b) + direction(b)*ads_test%s(b))* &
                           ((1 - direction(c))*ads_trial%s(c) + direction(c)*ads_test%s(c))))
#ifdef PERFORMANCE
         time1 = MPI_Wtime()
#endif
!  scatter back onto the cube of processors
         call Scatter(Ft_out, Ft2_out, (1 - direction(a))*ads_trial%n(a) + direction(a)*ads_test%n(a), &
                      (1 - direction(a))*ads_trial%s(a) + direction(a)*ads_test%s(a), &
                      ((1 - direction(b))*ads_trial%s(b) + direction(b)*ads_test%s(b))* &
                      ((1 - direction(c))*ads_trial%s(c) + direction(c)*ads_test%s(c)), &
                      dimensions_test, shifts_test, comm, ierr)
#ifdef PERFORMANCE
         time2 = MPI_Wtime()
         write (*, *) "Scatter ", a, ": ", time2 - time1
#endif
      else
         F_out = Fs
      end if

!  allocate buffers
      allocate (F2_out(ads_trial%n(a) + 1, &
                       ads_trial%s(b)*ads_trial%s(c)))
#ifdef PERFORMANCE
      time1 = MPI_Wtime()
#endif
!  scatter back onto the cube of processors
      call Scatter(F_out, F2_out, ads_trial%n(a), &
                   ads_trial%s(a), &
                   ads_trial%s(b)*ads_trial%s(c), &
                   dimensions_trial, shifts_trial, comm, ierr)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write (*, *) "Scatter ", a, ": ", time2 - time1
#endif
!  cleanup
      if (allocated(F_out)) deallocate (F_out)
      if (allocated(Ft_out)) deallocate (Ft_out)

      call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
      write (*, *) PRINTRANK, a, 'd) REORDER'
#endif
      ! Reorder right hand sides
      if (a .EQ. 1) call ReorderRHSForY(ads_trial%ibeg, ads_trial%iend, F2_out, F2)
      if (a .EQ. 2) call ReorderRHSForZ(ads_trial%ibeg, ads_trial%iend, F2_out, F2)
      if (a .EQ. 3) call ReorderRHSForX(ads_trial%ibeg, ads_trial%iend, F2_out, F2)

      if (igrm) then
         if (a .EQ. 1) call ReorderRHSForY(ibeg, iend, Ft2_out, Ft2)
         if (a .EQ. 2) call ReorderRHSForZ(ibeg, iend, Ft2_out, Ft2)
         if (a .EQ. 3) call ReorderRHSForX(ibeg, iend, Ft2_out, Ft2)
      end if
!  cleanup
      if (allocated(F2_out)) deallocate (F2_out)
      if (allocated(Ft2_out)) deallocate (Ft2_out)

#ifdef IPRINT
      write (*, *) PRINTRANK, 'after ReorderRHS'
      write (*, *) PRINTRANK, 'F:'
      do i = 1, ads_trial%s(1)
         write (*, *) PRINTRANK, i, 'row=', F2(i, 1:ads_trial%s(2)*ads_trial%s(3))
      end do
#endif

      ierr = 0

   end subroutine solve_problem

end module ADSS
