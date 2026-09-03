module igrm_pollution_solver_test_data

   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

   implicit none
   private

   real(kind=8), parameter, public :: POLLUTION_CENTER(3) = &
      (/3000.d0, 2000.d0, 2000.d0/)
   real(kind=8), parameter, public :: POLLUTION_RADIUS = 25.d0

   public :: constant_source, affine_source, zero_source, pollution_source
   public :: zero_wind, fixed_wind, physical_wind, nan_wind

contains

   function constant_source(un, du, point) result(value)
      real(kind=8), intent(in) :: un, du(3), point(3)
      real(kind=8) :: value

      value = 1.d0 + 0.d0*(un + sum(du) + sum(point))
   end function constant_source

   function affine_source(un, du, point) result(value)
      real(kind=8), intent(in) :: un, du(3), point(3)
      real(kind=8) :: value

      value = 1.d0 + point(1) + 2.d0*point(2) + 3.d0*point(3) + &
              0.d0*(un + sum(du))
   end function affine_source

   function zero_source(un, du, point) result(value)
      real(kind=8), intent(in) :: un, du(3), point(3)
      real(kind=8) :: value

      value = 0.d0*(un + sum(du) + sum(point))
   end function zero_source

   function pollution_source(un, du, point) result(value)
      real(kind=8), intent(in) :: un, du(3), point(3)
      real(kind=8) :: value, normalized(3), radius_squared

      normalized = (point - POLLUTION_CENTER)/POLLUTION_RADIUS
      radius_squared = min(sum(normalized*normalized), 1.d0)
      value = (radius_squared - 1.d0)**2*(radius_squared + 1.d0)**2
      value = value + 0.d0*(un + sum(du))
   end function pollution_source

   function zero_wind(time) result(value)
      real(kind=8), intent(in) :: time
      real(kind=8) :: value(3)

      value = 0.d0*time
   end function zero_wind

   function fixed_wind(time) result(value)
      real(kind=8), intent(in) :: time
      real(kind=8) :: value(3)

      value = (/0.4d0, -0.5d0, 0.6d0/) + 0.d0*time
   end function fixed_wind

   function physical_wind(time) result(value)
      real(kind=8), intent(in) :: time
      real(kind=8) :: value(3), angle, phase, scaled_time
      real(kind=8), parameter :: PI = &
         3.141592653589793238462643383279502884197d0

      scaled_time = time/150.d0
      phase = sin(scaled_time) + 0.5d0*sin(2.3d0*scaled_time)
      angle = PI/3.d0*phase + 3.d0*PI/8.d0
      value(1) = 5.d0*cos(angle)
      value(2) = 5.d0*sin(angle)
      value(3) = 5.d0/3.d0*sin(scaled_time)
   end function physical_wind

   function nan_wind(time) result(value)
      real(kind=8), intent(in) :: time
      real(kind=8) :: value(3)

      value = 0.d0*time
      value(2) = ieee_value(0.d0, ieee_quiet_nan)
   end function nan_wind

end module igrm_pollution_solver_test_data


program test_igrm_pollution_solver
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
      ieee_quiet_nan
   use mpi
   use parallelism, only: MYRANK, InitializeParallelism, Cleanup_Parallelism
   use pollution_dpg_solver, only: PollutionDPGSolution, PollutionDPGStats, &
      POLLUTION_DPG_SUCCESS, POLLUTION_DPG_INVALID_CONFIGURATION, &
      POLLUTION_DPG_SYSTEM_TOO_LARGE, POLLUTION_DPG_IO_ERROR, &
      POLLUTION_DPG_SOURCE_ERROR, SolvePollutionDPG3D, &
      EvaluatePollutionDPG, WritePollutionDPGVTI, &
      CleanupPollutionDPGSolution
   use igrm_pollution_solver_test_data, only: POLLUTION_CENTER, &
      POLLUTION_RADIUS, constant_source, affine_source, zero_source, &
      pollution_source, zero_wind, fixed_wind, physical_wind, nan_wind

   implicit none

   real(kind=8), parameter :: UNIT_MIN(3) = (/0.d0, 0.d0, 0.d0/)
   real(kind=8), parameter :: UNIT_MAX(3) = (/1.d0, 1.d0, 1.d0/)
   real(kind=8), parameter :: ZERO_DIFFUSION(3) = (/0.d0, 0.d0, 0.d0/)
   real(kind=8), parameter :: PHYSICAL_DIFFUSION(3) = &
      (/50.d0, 50.d0, 0.5d0/)
   real(kind=8), parameter :: PI = &
      3.141592653589793238462643383279502884197d0
   real(kind=8), parameter :: ROUND_OFF_TOLERANCE = 2.d-10

   integer(kind=4) :: checks, failures, global_failures
   integer(kind=4) :: ierr, cleanup_ierr, reduce_ierr, process_count
   character(len=32) :: process_count_argument

   checks = 0
   failures = 0
   process_count = 1
   if (command_argument_count() == 1) then
      call get_command_argument(1, process_count_argument)
      read (process_count_argument, *, iostat=ierr) process_count
      if (ierr /= 0) process_count = 1
   end if

   call InitializeParallelism(process_count, 1, 1, ierr)
   call check('MPI initialization succeeds', ierr == 0)
   if (ierr == 0) then
      call test_invalid_configuration_and_errors
      call test_directional_substep_oracle
      call test_affine_operator_oracle
      call test_enriched_test_space
      call test_physical_pollution_source
      call test_adapted_mesh
   end if

   call MPI_Allreduce(failures, global_failures, 1, MPI_INTEGER, MPI_SUM, &
                      MPI_COMM_WORLD, reduce_ierr)
   if (reduce_ierr /= MPI_SUCCESS) global_failures = failures + 1
   if (MYRANK == 0) then
      if (global_failures == 0) then
         write (*, '(A,I0,A)') &
            'OK (', checks, ' iGRM pollution solver checks per rank)'
      else
         write (*, '(A,I0,A)') &
            'FAILED (', global_failures, ' iGRM pollution solver failures)'
      end if
   end if

   call Cleanup_Parallelism(cleanup_ierr)
   if (global_failures /= 0) stop 1

contains

   subroutine test_invalid_configuration_and_errors
      type(PollutionDPGSolution) :: solution
      type(PollutionDPGStats) :: stats
      real(kind=8) :: reversed_min(3), reversed_max(3), invalid_diffusion(3)
      real(kind=8) :: invalid_source_min(3), nan_value
      integer(kind=4) :: status

      nan_value = ieee_value(0.d0, ieee_quiet_nan)

      call SolvePollutionDPG3D(1, .false., 1, 0, 1, 0, 0, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('zero time steps are rejected before solution allocation', &
         status == POLLUTION_DPG_INVALID_CONFIGURATION .and. &
         .not. allocated(solution%coefficients))

      call SolvePollutionDPG3D(2, .false., 2, 0, 1, 0, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('a smaller test space is rejected before allocation', &
         status == POLLUTION_DPG_INVALID_CONFIGURATION .and. &
         .not. allocated(solution%coefficients))

      call SolvePollutionDPG3D(2, .false., 1, -1, 2, 1, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('a broken trial space without DG fluxes is rejected', &
         status == POLLUTION_DPG_INVALID_CONFIGURATION .and. &
         .not. allocated(solution%coefficients))

      call SolvePollutionDPG3D(2, .false., 1, 0, 2, -1, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('a broken test space without DG fluxes is rejected', &
         status == POLLUTION_DPG_INVALID_CONFIGURATION .and. &
         .not. allocated(solution%coefficients))

      call SolvePollutionDPG3D(1, .false., 1, 0, 1, 0, 1, nan_value, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('a NaN time step is rejected before allocation', &
         status == POLLUTION_DPG_INVALID_CONFIGURATION .and. &
         .not. allocated(solution%coefficients))

      invalid_diffusion = ZERO_DIFFUSION
      invalid_diffusion(2) = nan_value
      call SolvePollutionDPG3D(1, .false., 1, 0, 1, 0, 1, 0.3d0, &
         1.d0, invalid_diffusion, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('a NaN diffusion coefficient is rejected before allocation', &
         status == POLLUTION_DPG_INVALID_CONFIGURATION .and. &
         .not. allocated(solution%coefficients))

      call SolvePollutionDPG3D(400, .false., 1, 0, 1, 0, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('an oversized tensor problem is rejected without allocation', &
         status == POLLUTION_DPG_SYSTEM_TOO_LARGE .and. &
         .not. allocated(solution%coefficients))

      call SolvePollutionDPG3D(200, .false., 1, 0, 9, 0, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('enriched source tensor is capped while trial cube fits', &
         status == POLLUTION_DPG_SYSTEM_TOO_LARGE .and. &
         .not. allocated(solution%coefficients))

      call SolvePollutionDPG3D(1, .false., 1, 0, 1, 0, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 1000, .false., solution, stats, status)
      call check('oversized output grid is rejected without allocation', &
         status == POLLUTION_DPG_SYSTEM_TOO_LARGE .and. &
         .not. allocated(solution%coefficients))

      reversed_min = (/0.8d0, 0.d0, 0.d0/)
      reversed_max = (/0.2d0, 1.d0, 1.d0/)
      call SolvePollutionDPG3D(1, .false., 1, 0, 1, 0, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, reversed_min, &
         reversed_max, 2, .false., solution, stats, status)
      call check('reversed source bounds are rejected as configuration', &
         status == POLLUTION_DPG_INVALID_CONFIGURATION)
      call CleanupPollutionDPGSolution(solution)

      invalid_source_min = UNIT_MIN
      invalid_source_min(3) = nan_value
      call SolvePollutionDPG3D(1, .false., 1, 0, 1, 0, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, &
         invalid_source_min, UNIT_MAX, 2, .false., solution, stats, status)
      call check('a NaN source bound is rejected before space allocation', &
         status == POLLUTION_DPG_INVALID_CONFIGURATION .and. &
         .not. allocated(solution%coefficients))

      call SolvePollutionDPG3D(1, .false., 1, 0, 1, 0, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, zero_source, UNIT_MIN, UNIT_MAX, &
         2, .false., solution, stats, status)
      call check('a source with zero total emission is rejected', &
         status == POLLUTION_DPG_SOURCE_ERROR)
      call CleanupPollutionDPGSolution(solution)

      call SolvePollutionDPG3D(1, .false., 1, 0, 1, 0, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, nan_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('a non-finite wind is rejected collectively', &
         status == POLLUTION_DPG_INVALID_CONFIGURATION)
      call CleanupPollutionDPGSolution(solution)

      call WritePollutionDPGVTI(solution, 'unused.vti', 2, status)
      call check('VTI output rejects an unallocated solution', &
         status == POLLUTION_DPG_IO_ERROR)
   end subroutine test_invalid_configuration_and_errors


   subroutine test_directional_substep_oracle
      type(PollutionDPGSolution) :: solution
      type(PollutionDPGStats) :: stats
      real(kind=8) :: expected, point(3), value
      integer(kind=4) :: status

      call SolvePollutionDPG3D(2, .false., 1, 0, 1, 0, 2, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 4, .false., solution, stats, status)
      call check('equal-space constant-source solve succeeds', &
         status == POLLUTION_DPG_SUCCESS)
      if (status /= POLLUTION_DPG_SUCCESS) return

      expected = 0.6d0
      call check('three directional substeps advance by dt/3 each', &
         maxval(abs(solution%coefficients - expected)) <= &
         ROUND_OFF_TOLERANCE)
      call check('constant-source diagnostics preserve exact mass and L2', &
         abs(stats%l2_norm - expected) <= ROUND_OFF_TOLERANCE .and. &
         abs(stats%total_mass - expected) <= ROUND_OFF_TOLERANCE .and. &
         abs(stats%maximum_abs - expected) <= ROUND_OFF_TOLERANCE)
      call check('constant-source run reports exact step, time and integral', &
         stats%completed_steps == 2 .and. &
         abs(stats%time - 0.6d0) <= ROUND_OFF_TOLERANCE .and. &
         abs(stats%source_integral - 1.d0) <= ROUND_OFF_TOLERANCE)
      call check('directional normal solves have a small residual', &
         ieee_is_finite(stats%maximum_relative_residual) .and. &
         stats%maximum_relative_residual <= 1.d-10)

      point = (/-1.d0, 2.d0, 0.37d0/)
      call EvaluatePollutionDPG(solution, point, value)
      call check('evaluation clamps coordinates and preserves constants', &
         abs(value - expected) <= ROUND_OFF_TOLERANCE)
      call check_replicated(solution, &
         'equal-space solution is replicated on every MPI rank')
      call test_vti_output(solution, expected)
      call test_vti_axis_order(solution)
      call test_vti_direct_validation(solution)

      call CleanupPollutionDPGSolution(solution)
      call check('solution cleanup releases fields and resets metadata', &
         .not. allocated(solution%coefficients) .and. &
         .not. allocated(solution%Ux) .and. &
         .not. allocated(solution%Uy) .and. &
         .not. allocated(solution%Uz) .and. all(solution%n == -1))
   end subroutine test_directional_substep_oracle


   subroutine test_affine_operator_oracle
      type(PollutionDPGSolution) :: solution
      type(PollutionDPGStats) :: stats
      real(kind=8), parameter :: ORACLE_DIFFUSION(3) = &
         (/0.1d0, 0.2d0, 0.3d0/)
      real(kind=8) :: expected(2, 2, 2)
      integer(kind=4) :: status

      ! Independent p=1 oracle obtained from the exact one-element matrices
      ! M=[[1/3,1/6],[1/6,1/3]], K=[[1,-1],[-1,1]] and
      ! A=[[-1/2,1/2],[-1/2,1/2]] for source 1+x+2y+3z.
      expected(1, 1, 1) = 0.52270262943887402d0
      expected(2, 1, 1) = 0.77984548658172825d0
      expected(1, 2, 1) = 0.96786391976145147d0
      expected(2, 2, 1) = 1.2250067769043116d0
      expected(1, 1, 2) = 1.1050555706153429d0
      expected(2, 1, 2) = 1.3621984277582018d0
      expected(1, 2, 2) = 1.5502168609379243d0
      expected(2, 2, 2) = 1.8073597180807792d0

      call SolvePollutionDPG3D(1, .false., 1, 0, 1, 0, 1, 0.3d0, &
         1.d0, ORACLE_DIFFUSION, fixed_wind, affine_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('affine advection-diffusion oracle solve succeeds', &
         status == POLLUTION_DPG_SUCCESS)
      if (status == POLLUTION_DPG_SUCCESS) then
         call check('affine oracle verifies signs and XYZ directional sweeps', &
            maxval(abs(solution%coefficients - expected)) <= 2.d-10)
         call check('affine oracle preserves its exact source integral', &
            abs(stats%source_integral - 4.d0) <= ROUND_OFF_TOLERANCE)
         call check('affine oracle records the fixed wind used by the step', &
            maxval(abs(stats%wind - fixed_wind(0.d0))) <= &
            ROUND_OFF_TOLERANCE)
         call check_replicated(solution, &
            'affine oracle solution is replicated on every MPI rank')
      end if
      call CleanupPollutionDPGSolution(solution)
   end subroutine test_affine_operator_oracle


   subroutine test_enriched_test_space
      type(PollutionDPGSolution) :: solution
      type(PollutionDPGStats) :: stats
      integer(kind=4) :: status

      call SolvePollutionDPG3D(2, .false., 1, 0, 2, 1, 1, 0.3d0, &
         1.d0, ZERO_DIFFUSION, zero_wind, constant_source, UNIT_MIN, &
         UNIT_MAX, 2, .false., solution, stats, status)
      call check('genuinely enriched p2/C1 test-space solve succeeds', &
         status == POLLUTION_DPG_SUCCESS)
      if (status == POLLUTION_DPG_SUCCESS) then
         call check('enriched DPG projection exactly reproduces a constant', &
            maxval(abs(solution%coefficients - 0.3d0)) <= &
            ROUND_OFF_TOLERANCE)
         call check('enriched solve returns the trial-space dimensions', &
            all(solution%n == (/2, 2, 2/)) .and. &
            all(shape(solution%coefficients) == (/3, 3, 3/)))
         call check_replicated(solution, &
            'enriched solution is replicated on every MPI rank')
      end if
      call CleanupPollutionDPGSolution(solution)
   end subroutine test_enriched_test_space


   subroutine test_physical_pollution_source
      type(PollutionDPGSolution) :: solution
      type(PollutionDPGStats) :: stats
      real(kind=8) :: source_min(3), source_max(3), exact_integral
      real(kind=8) :: expected_wind(3), value
      integer(kind=4) :: status

      source_min = POLLUTION_CENTER - POLLUTION_RADIUS
      source_max = POLLUTION_CENTER + POLLUTION_RADIUS
      exact_integral = 2.d6*PI/231.d0
      call SolvePollutionDPG3D(2, .false., 1, 0, 2, 1, 1, 1.8d0, &
         5000.d0, PHYSICAL_DIFFUSION, physical_wind, pollution_source, &
         source_min, source_max, 2, .false., solution, stats, status)
      call check('physical compact-source pollution solve succeeds', &
         status == POLLUTION_DPG_SUCCESS)
      if (status /= POLLUTION_DPG_SUCCESS) return

      expected_wind = physical_wind(0.d0)
      call check('source-aware quadrature resolves the narrow emission', &
         abs(stats%source_integral - exact_integral) <= &
         5.d-6*exact_integral)
      call check('physical run reports finite nonzero diagnostics', &
         ieee_is_finite(stats%l2_norm) .and. stats%l2_norm > 0.d0 .and. &
         ieee_is_finite(stats%total_mass) .and. &
         abs(stats%total_mass) > 0.d0 .and. &
         ieee_is_finite(stats%maximum_abs) .and. &
         stats%maximum_abs > 0.d0 .and. &
         all(ieee_is_finite(solution%coefficients)))
      call check('physical run records the wind used at the first step', &
         maxval(abs(stats%wind - expected_wind)) <= &
         ROUND_OFF_TOLERANCE)
      call check('physical solve residual remains below solver tolerance', &
         ieee_is_finite(stats%maximum_relative_residual) .and. &
         stats%maximum_relative_residual <= 1.d-8)
      call EvaluatePollutionDPG(solution, POLLUTION_CENTER, value)
      call check('physical solution is finite and nonzero at source center', &
         ieee_is_finite(value) .and. abs(value) > 0.d0)
      call check_replicated(solution, &
         'physical pollution solution is replicated on every MPI rank')
      call CleanupPollutionDPGSolution(solution)
   end subroutine test_physical_pollution_source


   subroutine test_adapted_mesh
      type(PollutionDPGSolution) :: solution
      type(PollutionDPGStats) :: stats
      real(kind=8), parameter :: DOMAIN = 5000.d0
      real(kind=8) :: source_min(3), source_max(3), expected_knots(3)
      integer(kind=4) :: status

      source_min = 0.d0
      source_max = DOMAIN
      call SolvePollutionDPG3D(4, .true., 1, 0, 2, 1, 1, 0.1d0, &
         DOMAIN, ZERO_DIFFUSION, zero_wind, constant_source, source_min, &
         source_max, 2, .false., solution, stats, status)
      call check('adapted x-mesh solve succeeds', &
         status == POLLUTION_DPG_SUCCESS)
      if (status /= POLLUTION_DPG_SUCCESS) return

      expected_knots = (/2475.d0, 4950.d0, 4975.d0/)
      call check('x knots use the upstream 99/1 percent adaptation', &
         maxval(abs(solution%Ux(3:5) - expected_knots)) <= 1.d-10)
      expected_knots = (/1250.d0, 2500.d0, 3750.d0/)
      call check('y and z meshes remain uniform', &
         maxval(abs(solution%Uy(3:5) - expected_knots)) <= 1.d-10 .and. &
         maxval(abs(solution%Uz(3:5) - expected_knots)) <= 1.d-10)
      call check('adaptation preserves the constant-source oracle', &
         maxval(abs(solution%coefficients - 0.1d0)) <= 5.d-10)
      call check_replicated(solution, &
         'adapted solution is replicated on every MPI rank')
      call CleanupPollutionDPGSolution(solution)
   end subroutine test_adapted_mesh


   subroutine test_vti_output(solution, expected_value)
      type(PollutionDPGSolution), intent(in) :: solution
      real(kind=8), intent(in) :: expected_value
      integer(kind=4) :: status
      logical :: metadata_ok, values_ok

      if (MYRANK /= 0) then
         call check('VTI file ownership is restricted to rank zero', .true.)
         return
      end if

      call WritePollutionDPGVTI(solution, 'test_igrm_pollution_result.vti', &
                                4, status)
      call check('VTI output succeeds on every replicated solution', &
         status == POLLUTION_DPG_SUCCESS)
      if (status == POLLUTION_DPG_SUCCESS) then
         call inspect_vti('test_igrm_pollution_result.vti', expected_value, &
                          metadata_ok, values_ok)
         call check('VTI metadata uses the requested extent and spacing', &
            metadata_ok)
         call check('VTI stores exactly 5 cubed evaluated point values', &
            values_ok)
      end if

      call WritePollutionDPGVTI(solution, 'unused.vti', 0, status)
      call check('VTI output rejects a non-positive resolution', &
         status == POLLUTION_DPG_IO_ERROR)
      call WritePollutionDPGVTI(solution, &
         'missing_igrm_pollution_directory/result.vti', 2, status)
      call check('VTI output reports an open failure', &
         status == POLLUTION_DPG_IO_ERROR)
   end subroutine test_vti_output


   subroutine test_vti_axis_order(solution)
      type(PollutionDPGSolution), intent(in) :: solution
      type(PollutionDPGSolution) :: linear_solution
      integer(kind=4) :: ix, iy, iz, status, mpi_status
      logical :: local_order_ok, global_order_ok

      linear_solution = solution
      do iz = 1, linear_solution%n(3) + 1
         do iy = 1, linear_solution%n(2) + 1
            do ix = 1, linear_solution%n(1) + 1
               linear_solution%coefficients(ix, iy, iz) = &
                  real(ix - 1, kind=8)/real(linear_solution%n(1), kind=8) + &
                  10.d0*real(iy - 1, kind=8)/ &
                     real(linear_solution%n(2), kind=8) + &
                  100.d0*real(iz - 1, kind=8)/ &
                     real(linear_solution%n(3), kind=8)
            end do
         end do
      end do

      local_order_ok = .true.
      if (MYRANK == 0) then
         call WritePollutionDPGVTI(linear_solution, &
            'test_igrm_pollution_axes.vti', 4, status)
         local_order_ok = status == POLLUTION_DPG_SUCCESS
         if (local_order_ok) then
            call inspect_vti_axis_order('test_igrm_pollution_axes.vti', &
                                        local_order_ok)
         end if
      end if
      call MPI_Bcast(local_order_ok, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, &
                     mpi_status)
      global_order_ok = local_order_ok .and. mpi_status == MPI_SUCCESS
      call check('VTI values preserve the X-fastest XYZ axis order', &
                 global_order_ok)
      call CleanupPollutionDPGSolution(linear_solution)
   end subroutine test_vti_axis_order


   subroutine test_vti_direct_validation(solution)
      type(PollutionDPGSolution), intent(in) :: solution
      type(PollutionDPGSolution) :: invalid_solution
      integer(kind=4) :: status

      call WritePollutionDPGVTI(solution, 'unused.vti', 1000, status)
      call check('direct VTI output enforces its output-point limit', &
         status == POLLUTION_DPG_IO_ERROR)

      invalid_solution = solution
      invalid_solution%domain_length = &
         ieee_value(0.d0, ieee_quiet_nan)
      call WritePollutionDPGVTI(invalid_solution, 'unused.vti', 2, status)
      call check('direct VTI output rejects a non-finite domain', &
         status == POLLUTION_DPG_IO_ERROR)
      call CleanupPollutionDPGSolution(invalid_solution)
   end subroutine test_vti_direct_validation


   subroutine inspect_vti_axis_order(path, order_ok)
      character(len=*), intent(in) :: path
      logical, intent(out) :: order_ok
      character(len=512) :: line
      integer(kind=4), parameter :: RESOLUTION = 4
      integer(kind=4) :: unit_number, open_status, read_status, value_status
      integer(kind=4) :: value_count, ix, iy, iz
      real(kind=8) :: value, expected
      logical :: in_array

      order_ok = .false.
      in_array = .false.
      value_count = 0
      unit_number = 92
      open (unit=unit_number, file=path, status='old', action='read', &
            iostat=open_status)
      if (open_status /= 0) return

      do
         read (unit_number, '(A)', iostat=read_status) line
         if (read_status /= 0) exit
         if (index(line, '<DataArray Name="Result"') > 0) then
            in_array = .true.
         else if (index(line, '</DataArray>') > 0) then
            in_array = .false.
         else if (in_array) then
            read (line, *, iostat=value_status) value
            if (value_status /= 0) exit
            ix = mod(value_count, RESOLUTION + 1)
            iy = mod(value_count/(RESOLUTION + 1), RESOLUTION + 1)
            iz = value_count/((RESOLUTION + 1)*(RESOLUTION + 1))
            expected = real(ix, kind=8)/real(RESOLUTION, kind=8) + &
               10.d0*real(iy, kind=8)/real(RESOLUTION, kind=8) + &
               100.d0*real(iz, kind=8)/real(RESOLUTION, kind=8)
            if (abs(value - expected) > ROUND_OFF_TOLERANCE) exit
            value_count = value_count + 1
         end if
      end do
      close (unit_number, status='delete')
      order_ok = value_count == (RESOLUTION + 1)**3
   end subroutine inspect_vti_axis_order


   subroutine inspect_vti(path, expected_value, metadata_ok, values_ok)
      character(len=*), intent(in) :: path
      real(kind=8), intent(in) :: expected_value
      logical, intent(out) :: metadata_ok, values_ok
      character(len=512) :: line
      integer(kind=4) :: unit_number, open_status, read_status, value_status
      integer(kind=4) :: value_count
      real(kind=8) :: value
      logical :: extent_found, origin_found, spacing_found
      logical :: array_found, in_array, all_values_match

      metadata_ok = .false.
      values_ok = .false.
      extent_found = .false.
      origin_found = .false.
      spacing_found = .false.
      array_found = .false.
      in_array = .false.
      all_values_match = .true.
      value_count = 0
      unit_number = 91
      open (unit=unit_number, file=path, status='old', action='read', &
            iostat=open_status)
      if (open_status /= 0) return

      do
         read (unit_number, '(A)', iostat=read_status) line
         if (read_status /= 0) exit
         if (index(line, 'WholeExtent="0 4 0 4 0 4"') > 0) &
            extent_found = .true.
         if (index(line, 'Origin="0 0 0"') > 0) origin_found = .true.
         if (index(line, '2.5000000000000000E-01') > 0) &
            spacing_found = .true.
         if (index(line, '<DataArray Name="Result"') > 0) then
            array_found = .true.
            in_array = .true.
         else if (index(line, '</DataArray>') > 0) then
            in_array = .false.
         else if (in_array) then
            read (line, *, iostat=value_status) value
            if (value_status /= 0) then
               all_values_match = .false.
            else
               value_count = value_count + 1
               all_values_match = all_values_match .and. &
                  abs(value - expected_value) <= ROUND_OFF_TOLERANCE
            end if
         end if
      end do
      close (unit_number, status='delete')

      metadata_ok = extent_found .and. origin_found .and. spacing_found .and. &
                    array_found
      values_ok = value_count == 125 .and. all_values_match
   end subroutine inspect_vti


   subroutine check_replicated(solution, label)
      type(PollutionDPGSolution), intent(in) :: solution
      character(len=*), intent(in) :: label
      real(kind=8) :: local_checksum, minimum_checksum, maximum_checksum
      integer(kind=4) :: mpi_status

      local_checksum = sum(solution%coefficients)
      call MPI_Allreduce(local_checksum, minimum_checksum, 1, &
         MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpi_status)
      if (mpi_status == MPI_SUCCESS) then
         call MPI_Allreduce(local_checksum, maximum_checksum, 1, &
            MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_status)
      end if
      call check(label, mpi_status == MPI_SUCCESS .and. &
         abs(maximum_checksum - minimum_checksum) <= &
         ROUND_OFF_TOLERANCE*max(1.d0, abs(maximum_checksum)))
   end subroutine check_replicated


   subroutine check(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      checks = checks + 1
      if (.not. condition) then
         failures = failures + 1
         write (*, '(A,I0,2A)') 'rank ', MYRANK, ': FAIL: ', trim(label)
      end if
   end subroutine check

end program test_igrm_pollution_solver
