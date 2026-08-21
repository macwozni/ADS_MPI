module igrm_stokes_solver_test_data

   implicit none
   private
   public :: nonzero_forcing, zero_forcing, zero_velocity, &
             constant_x_velocity, zero_scalar

contains

   function nonzero_forcing(un, du, point) result(value)
      real(kind=8), intent(in) :: un, du(3), point(3)
      real(kind=8) :: value
      value = 1.d0 + point(1) + 2.d0*point(2) - point(3) + &
              0.d0*un + 0.d0*sum(du)
   end function nonzero_forcing

   function zero_forcing(un, du, point) result(value)
      real(kind=8), intent(in) :: un, du(3), point(3)
      real(kind=8) :: value
      value = 0.d0*(un + sum(du) + sum(point))
   end function zero_forcing

   function zero_velocity(point) result(value)
      real(kind=8), intent(in) :: point(3)
      real(kind=8) :: value(3)
      value = 0.d0*point
   end function zero_velocity

   function constant_x_velocity(point) result(value)
      real(kind=8), intent(in) :: point(3)
      real(kind=8) :: value(3)
      value = (/1.d0, 0.d0, 0.d0/) + 0.d0*point
   end function constant_x_velocity

   function zero_scalar(point) result(value)
      real(kind=8), intent(in) :: point(3)
      real(kind=8) :: value
      value = 0.d0*sum(point)
   end function zero_scalar

end module igrm_stokes_solver_test_data


program test_igrm_stokes_solver
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
      ieee_quiet_nan
   use parallelism, only: MYRANK, InitializeParallelism, Cleanup_Parallelism
   use sparse, only: sparse_matrix, clear_matrix, to_dense_matrix
   use igrm_stokes_solver, only: IGRMStokesSolution, IGRMStokesStats, &
      IGRM_STOKES_SUCCESS, IGRM_STOKES_INVALID_CONFIGURATION, &
      IGRM_STOKES_RESIDUAL_TOO_LARGE, IGRM_STOKES_IO_ERROR, &
      AssembleIGRMStokes3DSystem, SolveIGRMStokes3D, EvaluateIGRMStokes, &
      ComputeIGRMStokesErrors, WriteIGRMStokesVTI, &
      CleanupIGRMStokesSolution
   use igrm_stokes_solver_test_data, only: nonzero_forcing, zero_forcing, &
      zero_velocity, constant_x_velocity, zero_scalar

   implicit none

   integer(kind=4), parameter :: NELEM(3) = (/2, 2, 2/)
   integer(kind=4), parameter :: DEGREE(3) = (/1, 1, 1/)
   integer(kind=4), parameter :: TEST_SCALAR = 64
   integer(kind=4), parameter :: TRIAL_SCALAR = 27
   integer(kind=4), parameter :: SYSTEM_SIZE = 365
   integer(kind=4), parameter :: TEST_PRESSURE_BEGIN = 3*TEST_SCALAR
   integer(kind=4), parameter :: TRIAL_BEGIN = 4*TEST_SCALAR
   integer(kind=4), parameter :: TRIAL_PRESSURE_BEGIN = &
      TRIAL_BEGIN + 3*TRIAL_SCALAR
   integer(kind=4), parameter :: GAUGE_INDEX = SYSTEM_SIZE - 1
   real(kind=8), parameter :: TOLERANCE = 2.d-11

   integer(kind=4) :: checks, failures, ierr, cleanup_ierr, process_count
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
      call test_invalid_configuration
      call test_matrix_identities
      call test_boundary_penalty_uses_trial_degree
      call test_upstream_trial_continuity
      call test_direct_solve_and_output
   end if

   if (MYRANK == 0) then
      if (failures == 0) then
         write (*, '(A,I0,A)') 'OK (', checks, ' iGRM Stokes solver checks)'
      else
         write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                    ' iGRM Stokes solver checks)'
      end if
   end if
   call Cleanup_Parallelism(cleanup_ierr)
   if (failures /= 0) stop 1

contains

   subroutine test_invalid_configuration
      type(sparse_matrix), pointer :: matrix
      real(kind=8), allocatable :: rhs(:)
      type(IGRMStokesSolution) :: solution
      type(IGRMStokesStats) :: stats
      integer(kind=4) :: status, invalid_degree(3)
      real(kind=8) :: nan_value

      nullify (matrix)
      nan_value = ieee_value(0.d0, ieee_quiet_nan)
      call AssembleIGRMStokes3DSystem(NELEM, DEGREE, DEGREE, -1.d0, 1.d0, &
         nonzero_forcing, nonzero_forcing, nonzero_forcing, zero_velocity, &
         matrix, rhs, status)
      call check('assembly rejects non-positive viscosity', &
         status == IGRM_STOKES_INVALID_CONFIGURATION .and. &
         .not. associated(matrix) .and. .not. allocated(rhs))

      call AssembleIGRMStokes3DSystem(NELEM, DEGREE, DEGREE, 1.d0, &
         nan_value, nonzero_forcing, nonzero_forcing, nonzero_forcing, &
         zero_velocity, matrix, rhs, status)
      call check('assembly rejects NaN penalty factor', &
         status == IGRM_STOKES_INVALID_CONFIGURATION .and. &
         .not. associated(matrix) .and. .not. allocated(rhs))

      invalid_degree = (/1, 1, 1/)
      invalid_degree(2) = 2
      call AssembleIGRMStokes3DSystem(NELEM, DEGREE, invalid_degree, 1.d0, &
         1.d0, nonzero_forcing, nonzero_forcing, nonzero_forcing, &
         zero_velocity, matrix, rhs, status)
      call check('assembly rejects trial degree above test degree', &
         status == IGRM_STOKES_INVALID_CONFIGURATION .and. &
         .not. associated(matrix) .and. .not. allocated(rhs))

      invalid_degree = (/10, 1, 1/)
      call AssembleIGRMStokes3DSystem(NELEM, invalid_degree, DEGREE, 1.d0, &
         1.d0, nonzero_forcing, nonzero_forcing, nonzero_forcing, &
         zero_velocity, matrix, rhs, status)
      call check('assembly rejects degree beyond the Gauss table', &
         status == IGRM_STOKES_INVALID_CONFIGURATION .and. &
         .not. associated(matrix) .and. .not. allocated(rhs))

      call SolveIGRMStokes3D(NELEM, DEGREE, DEGREE, 1.d0, 1.d0, nan_value, &
         nonzero_forcing, nonzero_forcing, nonzero_forcing, zero_velocity, &
         solution, stats, status)
      call check('solve rejects NaN residual tolerance before allocation', &
         status == IGRM_STOKES_INVALID_CONFIGURATION .and. &
         .not. allocated(solution%coefficients))

      call SolveIGRMStokes3D(NELEM, DEGREE, DEGREE, 1.d0, 1.d0, 1.d-30, &
         nonzero_forcing, nonzero_forcing, nonzero_forcing, zero_velocity, &
         solution, stats, status)
      call check('solve rejects a residual above the requested tolerance', &
         status == IGRM_STOKES_RESIDUAL_TOO_LARGE .and. &
         .not. allocated(solution%coefficients))
   end subroutine test_invalid_configuration

   subroutine test_matrix_identities
      type(sparse_matrix), pointer :: matrix, changed_matrix
      real(kind=8), allocatable :: rhs(:), changed_rhs(:)
      real(kind=8), allocatable :: dense(:, :), changed_dense(:, :)
      real(kind=8) :: constant_pressure_defect, symmetry_defect
      integer(kind=4) :: status

      nullify (matrix, changed_matrix)
      call AssembleIGRMStokes3DSystem(NELEM, DEGREE, DEGREE, 1.d0, 1.d0, &
         nonzero_forcing, nonzero_forcing, nonzero_forcing, zero_velocity, &
         matrix, rhs, status)
      call check('C^-1/C0 system with interior facets assembles', &
         status == IGRM_STOKES_SUCCESS .and. associated(matrix))
      if (status /= IGRM_STOKES_SUCCESS .or. .not. associated(matrix)) return
      call check('small facet system has expected augmented dimensions', &
         matrix%x == SYSTEM_SIZE .and. matrix%y == SYSTEM_SIZE .and. &
         size(rhs) == SYSTEM_SIZE)
      call check('assembled system is genuinely sparse and nonempty', &
         matrix%total_entries > 0_8 .and. &
         matrix%total_entries < int(SYSTEM_SIZE, kind=8)**2)

      allocate (dense(0:SYSTEM_SIZE - 1, 0:SYSTEM_SIZE - 1))
      call to_dense_matrix(matrix, dense)
      symmetry_defect = maxval(abs(dense - transpose(dense)))
      call check('complete Gram/B/B^T/gauge matrix is symmetric', &
         symmetry_defect <= TOLERANCE)

      constant_pressure_defect = maxval(abs(sum( &
         dense(0:3*TEST_SCALAR - 1, &
               TRIAL_PRESSURE_BEGIN:TRIAL_PRESSURE_BEGIN + TRIAL_SCALAR - 1), &
         dim=2)))
      call check('constant trial pressure is annihilated by velocity rows', &
         constant_pressure_defect <= TOLERANCE)
      call check('mean-pressure gauge integrates the constant exactly', &
         abs(sum(dense(GAUGE_INDEX, &
            TRIAL_PRESSURE_BEGIN:TRIAL_PRESSURE_BEGIN + TRIAL_SCALAR - 1)) &
             - 1.d0) <= TOLERANCE)
      call check('gauge has zero diagonal and symmetric pressure coupling', &
         abs(dense(GAUGE_INDEX, GAUGE_INDEX)) <= TOLERANCE .and. &
         maxval(abs(dense(GAUGE_INDEX, :) - dense(:, GAUGE_INDEX))) <= &
         TOLERANCE)
      call check('pressure-test and all lower RHS blocks are zero', &
         maxval(abs(rhs(TEST_PRESSURE_BEGIN + 1:SYSTEM_SIZE))) <= TOLERANCE)

      call AssembleIGRMStokes3DSystem(NELEM, DEGREE, DEGREE, 1.d0, 2.d0, &
         nonzero_forcing, nonzero_forcing, nonzero_forcing, zero_velocity, &
         changed_matrix, changed_rhs, status)
      call check('second penalty-factor system assembles', &
         status == IGRM_STOKES_SUCCESS .and. associated(changed_matrix))
      if (associated(changed_matrix)) then
         allocate (changed_dense(0:SYSTEM_SIZE - 1, 0:SYSTEM_SIZE - 1))
         call to_dense_matrix(changed_matrix, changed_dense)
         call check('penalty factor changes velocity B facet blocks', &
            maxval(abs(changed_dense(0:3*TEST_SCALAR - 1, &
              TRIAL_BEGIN:TRIAL_BEGIN + 3*TRIAL_SCALAR - 1) - &
              dense(0:3*TEST_SCALAR - 1, &
              TRIAL_BEGIN:TRIAL_BEGIN + 3*TRIAL_SCALAR - 1))) > 1.d-8)
         deallocate (changed_dense, changed_rhs)
         call clear_matrix(changed_matrix)
      end if

      deallocate (dense, rhs)
      call clear_matrix(matrix)
      call check('sparse matrices are released cleanly', &
         .not. associated(matrix) .and. .not. associated(changed_matrix))
   end subroutine test_matrix_identities

   subroutine test_boundary_penalty_uses_trial_degree
      integer(kind=4), parameter :: ONE_ELEMENT(3) = (/1, 1, 1/)
      integer(kind=4), parameter :: QUADRATIC_TEST(3) = (/2, 2, 2/)
      integer(kind=4), parameter :: LINEAR_TRIAL(3) = (/1, 1, 1/)
      integer(kind=4), parameter :: QUADRATIC_TEST_SCALAR = 27
      type(sparse_matrix), pointer :: matrix_one, matrix_two
      real(kind=8), allocatable :: rhs_one(:), rhs_two(:)
      real(kind=8) :: observed_increment, expected_increment
      integer(kind=4) :: status

      nullify (matrix_one, matrix_two)
      call AssembleIGRMStokes3DSystem(ONE_ELEMENT, QUADRATIC_TEST, &
         LINEAR_TRIAL, 1.d0, 1.d0, zero_forcing, zero_forcing, zero_forcing, &
         constant_x_velocity, matrix_one, rhs_one, status)
      call check('p_test=2/p_trial=1 boundary system assembles', &
         status == IGRM_STOKES_SUCCESS .and. associated(matrix_one))
      if (status /= IGRM_STOKES_SUCCESS .or. &
          .not. associated(matrix_one)) return
      call AssembleIGRMStokes3DSystem(ONE_ELEMENT, QUADRATIC_TEST, &
         LINEAR_TRIAL, 1.d0, 2.d0, zero_forcing, zero_forcing, zero_forcing, &
         constant_x_velocity, matrix_two, rhs_two, status)
      call check('changed-factor p_test=2/p_trial=1 system assembles', &
         status == IGRM_STOKES_SUCCESS .and. associated(matrix_two))
      if (status == IGRM_STOKES_SUCCESS .and. associated(matrix_two)) then
         ! Six unit-square faces, h=sqrt(2), and eta increment
         ! (p_trial+1)(p_trial+2)=6.  Partition of unity turns the sum of
         ! velocity-x test rows into exactly 6*6/sqrt(2).
         observed_increment = sum(rhs_two(1:QUADRATIC_TEST_SCALAR) - &
                                  rhs_one(1:QUADRATIC_TEST_SCALAR))
         expected_increment = 36.d0/sqrt(2.d0)
         call check('Nitsche RHS eta is derived from p_trial, not p_test', &
            abs(observed_increment - expected_increment) <= 2.d-10)
         call check('constant x data add no y/z penalty increment', &
            maxval(abs(rhs_two(QUADRATIC_TEST_SCALAR + 1: &
                               3*QUADRATIC_TEST_SCALAR) - &
                       rhs_one(QUADRATIC_TEST_SCALAR + 1: &
                               3*QUADRATIC_TEST_SCALAR))) <= TOLERANCE)
         deallocate (rhs_two)
         call clear_matrix(matrix_two)
      end if
      deallocate (rhs_one)
      call clear_matrix(matrix_one)
   end subroutine test_boundary_penalty_uses_trial_degree

   subroutine test_upstream_trial_continuity
      integer(kind=4), parameter :: LOCAL_NELEM(3) = (/2, 1, 1/)
      integer(kind=4), parameter :: LOCAL_DEGREE(3) = (/3, 1, 1/)
      integer(kind=4), parameter :: EXPECTED_TEST_SCALAR = 32
      integer(kind=4), parameter :: EXPECTED_TRIAL_SCALAR = 24
      integer(kind=4), parameter :: EXPECTED_SYSTEM_SIZE = &
         4*EXPECTED_TEST_SCALAR + 4*EXPECTED_TRIAL_SCALAR + 1
      type(sparse_matrix), pointer :: matrix
      real(kind=8), allocatable :: rhs(:)
      integer(kind=4) :: status

      nullify (matrix)
      call AssembleIGRMStokes3DSystem(LOCAL_NELEM, LOCAL_DEGREE, &
         LOCAL_DEGREE, 1.d0, 1.d0, zero_forcing, zero_forcing, zero_forcing, &
         zero_velocity, matrix, rhs, status)
      call check('anisotropic cubic C1 trial system assembles', &
         status == IGRM_STOKES_SUCCESS .and. associated(matrix))
      if (status == IGRM_STOKES_SUCCESS .and. associated(matrix)) then
         call check('trial-space dimensions preserve upstream C1 continuity', &
            matrix%x == EXPECTED_SYSTEM_SIZE .and. &
            matrix%y == EXPECTED_SYSTEM_SIZE .and. &
            size(rhs) == EXPECTED_SYSTEM_SIZE)
         deallocate (rhs)
         call clear_matrix(matrix)
      end if
   end subroutine test_upstream_trial_continuity

   subroutine test_direct_solve_and_output
      type(IGRMStokesSolution) :: solution
      type(IGRMStokesStats) :: stats
      real(kind=8) :: point(3), velocity(3), pressure, divergence
      real(kind=8) :: velocity_l2, pressure_l2, divergence_l2
      integer(kind=4) :: status, unit_number, delete_status, read_status
      character(len=512) :: line
      logical :: exists, extent_found

      call SolveIGRMStokes3D(NELEM, DEGREE, DEGREE, 1.d0, 1.d0, 1.d-9, &
         nonzero_forcing, nonzero_forcing, nonzero_forcing, zero_velocity, &
         solution, stats, status)
      call check('real augmented MUMPS solve succeeds', &
         status == IGRM_STOKES_SUCCESS)
      if (status /= IGRM_STOKES_SUCCESS) return
      call check('solve returns all four replicated trial fields', &
         allocated(solution%coefficients) .and. &
         size(solution%coefficients, 1) == TRIAL_SCALAR .and. &
         size(solution%coefficients, 2) == 4)
      call check('reported residual is finite and below tolerance', &
         ieee_is_finite(stats%residual_rms) .and. &
         ieee_is_finite(stats%relative_residual) .and. &
         stats%relative_residual <= 1.d-9)
      call check('all solved coefficients are finite', &
         all(ieee_is_finite(solution%coefficients)))

      point = (/0.37d0, 0.41d0, 0.59d0/)
      call EvaluateIGRMStokes(solution, point, velocity, pressure, divergence)
      call check('field evaluation returns finite velocity, pressure and div', &
         all(ieee_is_finite(velocity)) .and. ieee_is_finite(pressure) .and. &
         ieee_is_finite(divergence))
      call ComputeIGRMStokesErrors(solution, zero_velocity, zero_scalar, &
         velocity_l2, pressure_l2, divergence_l2)
      call check('L2 postprocessing returns finite nonnegative diagnostics', &
         ieee_is_finite(velocity_l2) .and. ieee_is_finite(pressure_l2) .and. &
         ieee_is_finite(divergence_l2) .and. velocity_l2 >= 0.d0 .and. &
         pressure_l2 >= 0.d0 .and. divergence_l2 >= 0.d0)

      call WriteIGRMStokesVTI(solution, 'test_igrm_stokes_result.vti', status)
      call check('collective rank-zero VTI output succeeds', &
         status == IGRM_STOKES_SUCCESS)
      if (MYRANK == 0) then
         inquire (file='test_igrm_stokes_result.vti', exist=exists)
         call check('VTI output file exists on rank zero', exists)
         if (exists) then
            unit_number = 91
            open (unit=unit_number, file='test_igrm_stokes_result.vti', &
                  status='old', iostat=delete_status)
            extent_found = .false.
            if (delete_status == 0) then
               do
                  read (unit_number, '(A)', iostat=read_status) line
                  if (read_status /= 0) exit
                  if (index(line, &
                     'WholeExtent="0 50 0 50 0 50"') > 0) &
                     extent_found = .true.
               end do
               close (unit_number, status='delete')
            end if
            call check('VTI uses the upstream 51-cubed sampling grid', &
               extent_found)
         end if
      else
         call check('VTI rank-zero ownership check is not applicable', .true.)
      end if
      call WriteIGRMStokesVTI(solution, &
         'test_igrm_stokes_solver/result.vti', status)
      call check('VTI open failure is reported collectively', &
         status == IGRM_STOKES_IO_ERROR)
      call CleanupIGRMStokesSolution(solution)
      call check('solution cleanup releases coefficients and knots', &
         .not. allocated(solution%coefficients) .and. &
         .not. allocated(solution%Ux) .and. .not. allocated(solution%Uy) &
         .and. .not. allocated(solution%Uz))
   end subroutine test_direct_solve_and_output

   subroutine check(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      checks = checks + 1
      if (.not. condition) then
         failures = failures + 1
         write (*, '(A,I0,2A)') 'rank ', MYRANK, ': FAIL: ', trim(label)
      end if
   end subroutine check

end program test_igrm_stokes_solver
