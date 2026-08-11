!------------------------------------------------------------------------------
!
! MODULE: ads_lifecycle
!
! DESCRIPTION:
!> @file ads_lifecycle.F90
!> @brief ADS setup, allocation, decomposition, and cleanup routines.
!>
!> @details
!> This module owns the setup-time and teardown-time parts of the ADS
!> workflow. It is responsible for creating test and trial spline spaces,
!> preparing decomposition metadata, allocating basis and runtime buffers,
!> and releasing the allocated data at the end of a problem run.
!>
!> The public routines collected here used to live in the high-level
!> \ref ADSS module. Moving them out keeps \ref ADSS focused on time-step
!> orchestration while preserving the historical public entry points
!> re-exported by that module.
!>
!> The setup performed here is intentionally one-time work. Expensive
!> operations such as knot-vector preparation, basis-tabulation, and
!> iGRM-compatible buffer allocation are completed before the time loop
!> starts and are then reused by the solver path.
!
!------------------------------------------------------------------------------
module ads_lifecycle

   implicit none

   private
   public :: initialize
   public :: initialize_setup
   public :: ComputeDecomposition
   public :: AllocateADSdata
   public :: AllocateADS
   public :: Cleanup_data
   public :: Cleanup_ADS

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Initializes test and trial ADS spaces together with runtime
!> working buffers.
!>
!> @details
!> This routine computes the corresponding numbers of basis functions for
!> the test and trial spaces from the supplied element counts and
!> polynomial degrees, determines the quadrature sizes required by both
!> spaces, initializes the two setup structures through
!> \ref initialize_setup, and allocates the principal runtime buffers by
!> calling \ref AllocateADSdata.
!>
!> The trial-space quadrature order is increased where necessary so that
!> it is compatible with the test-space degree in each direction.
!
! Input:
! ------
!> @param[in] nelem
!> Numbers of elements in the three parametric directions.
!>
!> @param[in] p_test
!> Polynomial degrees of the test space.
!>
!> @param[in] p_trial
!> Polynomial degrees of the trial space.
!>
!> @param[in] continuity
!> Continuity-control parameter vector passed through the initialization
!> interface.
!
! Output:
! -------
!> @param[out] ads_test
!> Initialized setup structure for the test space.
!>
!> @param[out] ads_trial
!> Initialized setup structure for the trial space.
!>
!> @param[out] ads_data
!> Allocated runtime-data container.
!>
!> @param[out] mierr
!> Error/status code returned by the initialization workflow.
!
!---------------------------------------------------------------------------
subroutine initialize(nelem, p_test, p_trial, continuity, ads_test, ads_trial, ads_data, mierr)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use knot_vector, ONLY: PrepareKnot
      use basis, ONLY: BasisData
      use mpi
      implicit none
!> @brief Numbers of elements in the three directions.
      integer(kind=4), intent(in), dimension(3) :: nelem
!> @brief Polynomial degrees of the test and trial spaces.
      integer(kind=4), intent(in), dimension(3) :: p_test, p_trial
!> @brief Continuity-control vector passed into setup initialization.
      integer(kind=4), intent(in), dimension(3) :: continuity
!> @brief Output setup structure for the trial space.
      type(ADS_setup), intent(out) :: ads_trial
!> @brief Output setup structure for the test space.
      type(ADS_setup), intent(out) :: ads_test
!> @brief Output runtime-data container.
      type(ADS_compute_data), intent(out) :: ads_data
!> @brief Returned status code.
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

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Initializes a single ADS setup structure, including knot
!> vectors, decomposition metadata, and basis tables.
!>
!> @details
!> This routine prepares three open knot vectors, allocates the static
!> storage associated with the ADS setup, stores spline-space metadata,
!> computes the process-local domain decomposition, and precomputes
!> basis-function values and first derivatives at quadrature points by
!> calling \ref BasisData in each direction.
!>
!> The resulting setup structure contains:
!> - knot vectors,
!> - polynomial degrees and basis sizes,
!> - decomposition bounds and neighboring ranges,
!> - quadrature points and weights,
!> - nonzero basis-function tables on all elements.
!
! Input:
! ------
!> @param[in] n
!> Numbers of basis functions minus one in the three directions.
!>
!> @param[in] p
!> Polynomial degrees in the three directions.
!>
!> @param[in] continuity
!> Continuity-control vector passed through the interface.
!>
!> @param[in] ng
!> Numbers of quadrature points in the three directions.
!
! Output:
! -------
!> @param[out] ads
!> Initialized ADS setup structure.
!>
!> @param[out] mierr
!> Error/status code returned by the initialization workflow.
!
! Notes:
! ------
!> @note
!> The current implementation accepts \p continuity in the interface but
!> does not use it explicitly in the body of the routine.
!
!> @warning
!> The routine stops execution if the number of degrees of freedom in any
!> direction is smaller than the number of MPI processes assigned to that
!> direction. Empty ownership ranges are not supported.
!
!---------------------------------------------------------------------------
subroutine initialize_setup(n, p, continuity, ng, ads, mierr)
      use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use parallelism, ONLY: MYRANK, NRPROCX, NRPROCY, NRPROCZ
      use knot_vector, ONLY: PrepareKnot
      use basis, ONLY: BasisData
      use mpi
      implicit none
!> @brief Numbers of basis functions minus one.
      integer(kind=4), intent(in), dimension(3) :: n
!> @brief Polynomial degrees in the three directions.
      integer(kind=4), intent(in), dimension(3) :: p
!> @brief Continuity-control vector passed through the interface.
      integer(kind=4), intent(in), dimension(3) :: continuity
!> @brief Numbers of quadrature points in the three directions.
      integer(kind=4), intent(in), dimension(3) :: ng
!> @brief Output ADS setup structure.
      type(ADS_setup), intent(out) :: ads
!> @brief Returned status code.
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
      ads%tau = 0.d0

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

      if (n(1) < NRPROCX - 1 .or. n(2) < NRPROCY - 1 .or. &
          n(3) < NRPROCZ - 1) then
            if (MYRANK == 0) then
                  write (ERROR_UNIT, *) &
                     'Number of degrees of freedom smaller than number of processors'
                  flush (ERROR_UNIT)
            end if
            ! Every rank evaluates the same global size condition.  Complete
            ! rank zero's diagnostic before another rank can terminate the
            ! MPI job, otherwise MPI_Abort may discard the message.
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
            STOP 4
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

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes the process-local domain decomposition and neighbor
!> overlap metadata.
!>
!> @details
!> This routine determines, for the current MPI rank in each logical
!> direction:
!> - the owned basis-function range,
!> - the corresponding element range,
!> - the number of locally stored coefficients,
!> - dimension and shift vectors used by gather/scatter communication,
!> - neighboring overlap ranges required for stencil-like exchange.
!>
!> The calculations are delegated to helper routines from module
!> `parallelism`, in particular `ComputeEndpoints` and `FillDimVector`.
!
! Input/Output:
! -------------
!> @param[inout] ads
!> ADS setup structure updated with decomposition information.
!
!---------------------------------------------------------------------------
subroutine ComputeDecomposition(ads)
      use Setup, ONLY: ADS_Setup
      use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, &
                              NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints, FillDimVector
      implicit none
!> @brief ADS setup structure updated in place.
      type(ADS_setup), intent(inout) :: ads
      integer(kind=4) :: i
      integer(kind=4) :: ix, iy, iz
      integer(kind=4) :: dummy_nrcpp, imine, imaxe

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
            call ComputeEndpoints(i - 1, NRPROCX, ads%n(1), ads%p(1), dummy_nrcpp, ads%ibegsx(ix), &
                                    ads%iendsx(ix), imine, imaxe)
      end do
      do i = max(MYRANKY - 1, 0) + 1, min(MYRANKY + 1, NRPROCY - 1) + 1
            iy = i - MYRANKY + 1
            call ComputeEndpoints(i - 1, NRPROCY, ads%n(2), ads%p(2), dummy_nrcpp, ads%ibegsy(iy), &
                                    ads%iendsy(iy), imine, imaxe)
      end do
      do i = max(MYRANKZ - 1, 0) + 1, min(MYRANKZ + 1, NRPROCZ - 1) + 1
            iz = i - MYRANKZ + 1
            call ComputeEndpoints(i - 1, NRPROCZ, ads%n(3), ads%p(3), dummy_nrcpp, ads%ibegsz(iz), &
                                    ads%iendsz(iz), imine, imaxe)
      end do

end subroutine ComputeDecomposition

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Allocates the principal runtime arrays used during ADS
!> computations.
!>
!> @details
!> This routine allocates solution buffers sampled on element/quadrature
!> grids, the compact coefficient halo \p R, and a reusable sparse MPI
!> exchange plan used during local reconstruction.
!>
!> The halo bounds are derived from the union of test/trial element ranges
!> and the actual trial-basis support tables.  Ownership and demand boxes
!> are gathered once so later time steps communicate only their pairwise
!> intersections.
!
! Input:
! ------
!> @param[in] ads_test
!> Setup structure of the test space.
!>
!> @param[in] ads_trial
!> Setup structure of the trial space.
!
! Output:
! -------
!> @param[out] ads_data
!> Runtime-data structure with allocated core buffers.
!
!---------------------------------------------------------------------------
subroutine AllocateADSdata(ads_test, ads_trial, ads_data)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use parallelism, ONLY: NRPROC
      use mpi
      implicit none
!> @brief Test and trial setup structures.
      type(ADS_setup), intent(in) :: ads_test, ads_trial
!> @brief Output runtime-data container.
      type(ADS_compute_data), intent(out) :: ads_data
      integer :: ierr
      integer(kind=4), dimension(3) :: state_ng
      integer(kind=4), dimension(3) :: send_extent, recv_extent
      integer(kind=4), dimension(12) :: local_boxes
      integer(kind=4), allocatable, dimension(:, :) :: all_boxes
      integer(kind=4) :: peer, total_send, total_recv

      ads_data%state_mine = min(ads_test%mine, ads_trial%mine)
      ads_data%state_maxe = max(ads_test%maxe, ads_trial%maxe)
      ads_data%state_lnelem = ads_data%state_maxe - ads_data%state_mine + 1
      state_ng = max(ads_test%ng, ads_trial%ng)

      allocate (ads_data%Un(ads_data%state_lnelem(1), ads_data%state_lnelem(2), ads_data%state_lnelem(3), &
                              state_ng(1), state_ng(2), state_ng(3)))
      allocate (ads_data%Un13(ads_data%state_lnelem(1), ads_data%state_lnelem(2), ads_data%state_lnelem(3), &
                              state_ng(1), state_ng(2), state_ng(3)))
      allocate (ads_data%Un23(ads_data%state_lnelem(1), ads_data%state_lnelem(2), ads_data%state_lnelem(3), &
                              state_ng(1), state_ng(2), state_ng(3)))
      allocate (ads_data%dUn(ads_data%state_lnelem(1), ads_data%state_lnelem(2), ads_data%state_lnelem(3), &
                                state_ng(1), state_ng(2), state_ng(3), 3))
      allocate (ads_data%dUn0(ads_data%state_lnelem(1), ads_data%state_lnelem(2), ads_data%state_lnelem(3), &
                                state_ng(1), state_ng(2), state_ng(3), 3))
      allocate (ads_data%dUn13(ads_data%state_lnelem(1), ads_data%state_lnelem(2), ads_data%state_lnelem(3), &
                                state_ng(1), state_ng(2), state_ng(3), 3))
      allocate (ads_data%dUn23(ads_data%state_lnelem(1), ads_data%state_lnelem(2), ads_data%state_lnelem(3), &
                                state_ng(1), state_ng(2), state_ng(3), 3))
      ads_data%Un = 0.d0
      ads_data%Un13 = 0.d0
      ads_data%Un23 = 0.d0
      ads_data%dUn = 0.d0
      ads_data%dUn0 = 0.d0
      ads_data%dUn13 = 0.d0
      ads_data%dUn23 = 0.d0
      ads_data%rhs_du_state = 0

      ! The state buffers cover the union of the test and trial element
      ! ranges.  Reconstructing the trial solution on that union requires
      ! every trial coefficient supported on those elements.
      ads_data%halo_begin = (/ &
            ads_trial%Ox(ads_data%state_mine(1)), &
            ads_trial%Oy(ads_data%state_mine(2)), &
            ads_trial%Oz(ads_data%state_mine(3))/)
      ads_data%halo_end = (/ &
            ads_trial%Ox(ads_data%state_maxe(1)) + ads_trial%p(1), &
            ads_trial%Oy(ads_data%state_maxe(2)) + ads_trial%p(2), &
            ads_trial%Oz(ads_data%state_maxe(3)) + ads_trial%p(3)/)

      ! Gather ownership and demand boxes once.  The resulting peer plan is
      ! reused by every time step; no global metadata exchange is needed in
      ! the reconstruction hot path.
      local_boxes(1:3) = ads_trial%ibeg - 1
      local_boxes(4:6) = ads_trial%iend - 1
      local_boxes(7:9) = ads_data%halo_begin
      local_boxes(10:12) = ads_data%halo_end
      allocate (all_boxes(12, NRPROC))
      call mpi_allgather(local_boxes, 12, MPI_INTEGER, all_boxes, 12, &
                         MPI_INTEGER, MPI_COMM_WORLD, ierr)

      allocate (ads_data%halo_send_begin(3, NRPROC), ads_data%halo_send_end(3, NRPROC))
      allocate (ads_data%halo_recv_begin(3, NRPROC), ads_data%halo_recv_end(3, NRPROC))
      allocate (ads_data%halo_send_count(NRPROC), ads_data%halo_send_displ(NRPROC))
      allocate (ads_data%halo_recv_count(NRPROC), ads_data%halo_recv_displ(NRPROC))

      do peer = 1, NRPROC
            ads_data%halo_send_begin(:, peer) = max(local_boxes(1:3), all_boxes(7:9, peer))
            ads_data%halo_send_end(:, peer) = min(local_boxes(4:6), all_boxes(10:12, peer))
            send_extent = max(ads_data%halo_send_end(:, peer) - &
                              ads_data%halo_send_begin(:, peer) + 1, 0)
            ads_data%halo_send_count(peer) = product(send_extent)

            ads_data%halo_recv_begin(:, peer) = max(all_boxes(1:3, peer), local_boxes(7:9))
            ads_data%halo_recv_end(:, peer) = min(all_boxes(4:6, peer), local_boxes(10:12))
            recv_extent = max(ads_data%halo_recv_end(:, peer) - &
                              ads_data%halo_recv_begin(:, peer) + 1, 0)
            ads_data%halo_recv_count(peer) = product(recv_extent)
      end do

      ads_data%halo_send_displ(1) = 0
      ads_data%halo_recv_displ(1) = 0
      do peer = 2, NRPROC
            ads_data%halo_send_displ(peer) = ads_data%halo_send_displ(peer - 1) + &
                                             ads_data%halo_send_count(peer - 1)
            ads_data%halo_recv_displ(peer) = ads_data%halo_recv_displ(peer - 1) + &
                                             ads_data%halo_recv_count(peer - 1)
      end do
      total_send = ads_data%halo_send_displ(NRPROC) + ads_data%halo_send_count(NRPROC)
      total_recv = ads_data%halo_recv_displ(NRPROC) + ads_data%halo_recv_count(NRPROC)
      allocate (ads_data%halo_send_buffer(max(total_send, 1)))
      allocate (ads_data%halo_recv_buffer(max(total_recv, 1)))
      allocate (ads_data%halo_requests(max(2*(NRPROC - 1), 1)))
      allocate (ads_data%halo_statuses(MPI_STATUS_SIZE, max(2*(NRPROC - 1), 1)))
      ads_data%halo_send_buffer = 0.d0
      ads_data%halo_recv_buffer = 0.d0
      deallocate (all_boxes)

      allocate (ads_data%R( &
            ads_data%halo_end(1) - ads_data%halo_begin(1) + 1, &
            ads_data%halo_end(2) - ads_data%halo_begin(2) + 1, &
            ads_data%halo_end(3) - ads_data%halo_begin(3) + 1, 1))
      ads_data%R = 0.d0

end subroutine AllocateADSdata

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Allocates the mostly static arrays stored inside an ADS setup
!> structure.
!>
!> @details
!> This routine allocates basis-data arrays, quadrature arrays, element
!> Jacobians, and optionally pivot vectors used by linear solvers on
!> boundary-aligned ranks.
!>
!> The routine does not fill these arrays. Population is performed later
!> during the setup phase.
!
! Input:
! ------
!> @param[in] n
!> Numbers of basis functions minus one in the three directions.
!>
!> @param[in] nelem
!> Numbers of elements in the three directions.
!>
!> @param[in] p
!> Polynomial degrees in the three directions.
!>
!> @param[in] ng
!> Numbers of quadrature points in the three directions.
!
! Output:
! -------
!> @param[out] ads
!> ADS setup structure with allocated storage.
!
!---------------------------------------------------------------------------
subroutine AllocateADS(n, nelem, p, ng, ads)
      use Setup, ONLY: ADS_Setup
      use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ
      use mpi
      implicit none
!> @brief Basis sizes, element counts, degrees, and quadrature sizes.
      integer(kind=4), dimension(3), intent(in) :: n, nelem, p, ng
!> @brief Output ADS setup structure.
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

      allocate (ads%IPIVx(n(1) + 1))
      allocate (ads%IPIVy(n(2) + 1))
      allocate (ads%IPIVz(n(3) + 1))
      
end subroutine AllocateADS

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Deallocates runtime ADS data buffers.
!>
!> @details
!> This routine releases all allocatable arrays stored in the runtime
!> data structure, including right-hand-side buffers, solution-history
!> arrays, derivative arrays, and redistributed coefficient blocks.
!
! Input/Output:
! -------------
!> @param[inout] ads_data
!> Runtime-data structure to be cleaned.
!
! Output:
! -------
!> @param[out] mierr
!> Returned status code.
!
!---------------------------------------------------------------------------
subroutine Cleanup_data(ads_data, mierr)
      use Setup, ONLY: ADS_compute_data
      ! use parallelism, ONLY: PRINTRANK
      implicit none
!> @brief Runtime-data structure cleaned in place.
      type(ADS_compute_data), intent(inout) :: ads_data
!> @brief Returned status code.
      integer(kind=4), intent(out) :: mierr
      integer(kind=4) :: ierr

      if (allocated(ads_data%F)) deallocate (ads_data%F)
      if (allocated(ads_data%FF)) deallocate (ads_data%FF)
      if (allocated(ads_data%F2)) deallocate (ads_data%F2)
      if (allocated(ads_data%F3)) deallocate (ads_data%F3)

      if (allocated(ads_data%Ft)) deallocate (ads_data%Ft)
      if (allocated(ads_data%FFt)) deallocate (ads_data%FFt)
      if (allocated(ads_data%Ft2)) deallocate (ads_data%Ft2)
      if (allocated(ads_data%Ft3)) deallocate (ads_data%Ft3)

      if (allocated(ads_data%Un)) deallocate (ads_data%Un)
      if (allocated(ads_data%Un13)) deallocate (ads_data%Un13)
      if (allocated(ads_data%Un23)) deallocate (ads_data%Un23)
      if (allocated(ads_data%dUn)) deallocate (ads_data%dUn)
      if (allocated(ads_data%dUn0)) deallocate (ads_data%dUn0)
      if (allocated(ads_data%dUn13)) deallocate (ads_data%dUn13)
      if (allocated(ads_data%dUn23)) deallocate (ads_data%dUn23)

      if (allocated(ads_data%R)) deallocate (ads_data%R)
      if (allocated(ads_data%halo_send_begin)) deallocate (ads_data%halo_send_begin)
      if (allocated(ads_data%halo_send_end)) deallocate (ads_data%halo_send_end)
      if (allocated(ads_data%halo_recv_begin)) deallocate (ads_data%halo_recv_begin)
      if (allocated(ads_data%halo_recv_end)) deallocate (ads_data%halo_recv_end)
      if (allocated(ads_data%halo_send_count)) deallocate (ads_data%halo_send_count)
      if (allocated(ads_data%halo_send_displ)) deallocate (ads_data%halo_send_displ)
      if (allocated(ads_data%halo_recv_count)) deallocate (ads_data%halo_recv_count)
      if (allocated(ads_data%halo_recv_displ)) deallocate (ads_data%halo_recv_displ)
      if (allocated(ads_data%halo_send_buffer)) deallocate (ads_data%halo_send_buffer)
      if (allocated(ads_data%halo_recv_buffer)) deallocate (ads_data%halo_recv_buffer)
      if (allocated(ads_data%halo_requests)) deallocate (ads_data%halo_requests)
      if (allocated(ads_data%halo_statuses)) deallocate (ads_data%halo_statuses)
#ifdef IINFO
      write (*, *) PRINTRANK, "Exiting..."
#endif

      mierr = 0

end subroutine Cleanup_data

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Deallocates the arrays stored inside an ADS setup structure.
!>
!> @details
!> This routine releases knot vectors, decomposition metadata, pivot
!> vectors, basis tables, quadrature arrays, and related storage owned by
!> the setup structure.
!
! Input/Output:
! -------------
!> @param[inout] ads
!> ADS setup structure to be cleaned.
!
! Output:
! -------
!> @param[out] mierr
!> Returned status code.
!
!---------------------------------------------------------------------------
subroutine Cleanup_ADS(ads, mierr)
      use Setup, ONLY: ADS_Setup
      ! use parallelism, ONLY: PRINTRANK
      implicit none
!> @brief ADS setup structure cleaned in place.
      type(ADS_setup), intent(inout) :: ads
!> @brief Returned status code.
      integer(kind=4), intent(out) :: mierr
      ! integer(kind=4) :: ierr

      if (allocated(ads%Ux)) deallocate (ads%Ux)
      if (allocated(ads%Uy)) deallocate (ads%Uy)
      if (allocated(ads%Uz)) deallocate (ads%Uz)

      if (allocated(ads%dimensionsX)) deallocate (ads%dimensionsX)
      if (allocated(ads%dimensionsY)) deallocate (ads%dimensionsY)
      if (allocated(ads%dimensionsZ)) deallocate (ads%dimensionsZ)

      if (allocated(ads%shiftsX)) deallocate (ads%shiftsX)
      if (allocated(ads%shiftsY)) deallocate (ads%shiftsY)
      if (allocated(ads%shiftsZ)) deallocate (ads%shiftsZ)

      if (allocated(ads%IPIVx)) deallocate (ads%IPIVx)
      if (allocated(ads%IPIVy)) deallocate (ads%IPIVy)
      if (allocated(ads%IPIVz)) deallocate (ads%IPIVz)

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

end module ads_lifecycle
