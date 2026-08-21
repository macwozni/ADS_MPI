module igrm_eirksson_solver_test_forcing

   implicit none
   private

   public :: positive_forcing

contains

   real(kind=8) function positive_forcing(un, du, point) result(value)
      real(kind=8), intent(in) :: un, du(3), point(3)

      value = 1.d0 + point(1) + point(2) + point(3) + &
              0.d0*un + 0.d0*sum(du)
   end function positive_forcing

end module igrm_eirksson_solver_test_forcing


program test_igrm_eirksson_solver
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf
   use Setup, only: ADS_Setup, ADS_compute_data
   use parallelism, only: MYRANK, InitializeParallelism, Cleanup_Parallelism
   use ads_lifecycle, only: initialize, Cleanup_ADS, Cleanup_data
   use sparse, only: sparse_matrix, clear_matrix, to_dense_matrix
   use igrm_mumps_solver, only: IGRMMumpsStats, &
      IGRM_MUMPS_SUCCESS, IGRM_MUMPS_INVALID_CONFIGURATION, &
      AssembleIGRMMumps3DSystem, SolveIGRMMumps3D
   use igrm_eirksson_solver_test_forcing, only: positive_forcing

   implicit none

   integer(kind=4), parameter :: NELEM(3) = (/2, 2, 2/)
   integer(kind=4), parameter :: P_TEST(3) = (/2, 2, 2/)
   integer(kind=4), parameter :: P_TRIAL(3) = (/1, 1, 1/)
   integer(kind=4), parameter :: TEST_SIZE = 64
   integer(kind=4), parameter :: TRIAL_SIZE = 27
   integer(kind=4), parameter :: SYSTEM_SIZE = 91
   integer(kind=8), parameter :: EXPECTED_NNZ = 370_8
   real(kind=8), parameter :: EPSILON_VALUE = 1.d-2
   real(kind=8), parameter :: BETA(3) = (/1.d0, 1.d0, 1.d0/)
   real(kind=8), parameter :: STRUCTURE_TOLERANCE = 1.d-12
   real(kind=8), parameter :: SOLVER_TOLERANCE = 1.d-9

   type(ADS_Setup) :: ads_test, ads_trial
   type(ADS_compute_data) :: ads_data
   integer(kind=4) :: checks, failures, ierr, cleanup_ierr

   checks = 0
   failures = 0

   call InitializeParallelism(1, 1, 1, ierr)
   call check('one-rank MPI initialization succeeds', ierr == 0)

   call initialize(NELEM, P_TEST, P_TRIAL, P_TRIAL - 1, ads_test, &
                   ads_trial, ads_data, ierr)
   call check('real ADS test and trial spaces initialize', ierr == 0)

   if (ierr == 0) then
      call test_invalid_configuration
      call test_assembled_system
      call test_direct_solve
   end if

   call Cleanup_ADS(ads_test, cleanup_ierr)
   call check('test-space cleanup succeeds', cleanup_ierr == 0)
   call Cleanup_ADS(ads_trial, cleanup_ierr)
   call check('trial-space cleanup succeeds', cleanup_ierr == 0)
   call Cleanup_data(ads_data, cleanup_ierr)
   call check('runtime-data cleanup succeeds', cleanup_ierr == 0)

   if (MYRANK == 0) then
      if (failures == 0) then
         write (*, '(A,I0,A)') 'OK (', checks, &
                               ' iGRM Eriksson solver checks)'
      else
         write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                    ' iGRM Eriksson solver checks)'
      end if
   end if

   call Cleanup_Parallelism(cleanup_ierr)
   if (failures /= 0) stop 1

contains

   subroutine test_invalid_configuration
      type(sparse_matrix), pointer :: matrix
      real(kind=8), allocatable :: right_hand_side(:)
      real(kind=8), allocatable :: solution(:, :)
      type(ADS_Setup) :: invalid_trial
      type(IGRMMumpsStats) :: stats
      real(kind=8) :: invalid_beta(3), nan_value, positive_infinity
      integer(kind=4) :: status

      nan_value = ieee_value(0.d0, ieee_quiet_nan)
      positive_infinity = ieee_value(0.d0, ieee_positive_inf)
      nullify (matrix)
      call AssembleIGRMMumps3DSystem(ads_test, ads_trial, -EPSILON_VALUE, &
                                     BETA, positive_forcing, matrix, &
                                     right_hand_side, status)

      call check('assembly rejects a non-positive diffusion coefficient', &
                 status == IGRM_MUMPS_INVALID_CONFIGURATION)
      call check('rejected assembly does not allocate a sparse matrix', &
                 .not. associated(matrix))
      call check('rejected assembly does not allocate a right-hand side', &
                 .not. allocated(right_hand_side))

      call AssembleIGRMMumps3DSystem(ads_test, ads_trial, positive_infinity, &
                                     BETA, positive_forcing, matrix, &
                                     right_hand_side, status)
      call check('assembly rejects an infinite diffusion coefficient', &
                 status == IGRM_MUMPS_INVALID_CONFIGURATION .and. &
                 .not. associated(matrix) .and. &
                 .not. allocated(right_hand_side))

      invalid_beta = BETA
      invalid_beta(2) = nan_value
      call AssembleIGRMMumps3DSystem(ads_test, ads_trial, EPSILON_VALUE, &
                                     invalid_beta, positive_forcing, matrix, &
                                     right_hand_side, status)
      call check('assembly rejects a NaN advection coefficient', &
                 status == IGRM_MUMPS_INVALID_CONFIGURATION .and. &
                 .not. associated(matrix) .and. &
                 .not. allocated(right_hand_side))

      call SolveIGRMMumps3D(ads_test, ads_trial, EPSILON_VALUE, BETA, &
                            nan_value, positive_forcing, solution, stats, status)
      call check('solve rejects a NaN residual tolerance before allocation', &
                 status == IGRM_MUMPS_INVALID_CONFIGURATION .and. &
                 .not. allocated(solution))

      invalid_trial = ads_trial
      invalid_trial%p(1) = 0
      call AssembleIGRMMumps3DSystem(ads_test, invalid_trial, EPSILON_VALUE, &
                                     BETA, positive_forcing, matrix, &
                                     right_hand_side, status)
      call check('public assembly rejects a zero trial degree', &
                 status == IGRM_MUMPS_INVALID_CONFIGURATION .and. &
                 .not. associated(matrix) .and. &
                 .not. allocated(right_hand_side))

      invalid_trial = ads_trial
      invalid_trial%ng(3) = 0
      call AssembleIGRMMumps3DSystem(ads_test, invalid_trial, EPSILON_VALUE, &
                                     BETA, positive_forcing, matrix, &
                                     right_hand_side, status)
      call check('public assembly rejects a non-positive quadrature size', &
                 status == IGRM_MUMPS_INVALID_CONFIGURATION .and. &
                 .not. associated(matrix) .and. &
                 .not. allocated(right_hand_side))
   end subroutine test_invalid_configuration


   subroutine test_assembled_system
      type(sparse_matrix), pointer :: matrix
      real(kind=8), allocatable :: right_hand_side(:)
      real(kind=8), allocatable :: dense(:, :)
      integer(kind=4) :: status
      logical :: boundary_rows_ok, gram_ok, saddle_ok, lower_rhs_ok

      nullify (matrix)
      call AssembleIGRMMumps3DSystem(ads_test, ads_trial, EPSILON_VALUE, &
                                     BETA, positive_forcing, matrix, &
                                     right_hand_side, status)
      call check('full 3D saddle system assembles', &
                 status == IGRM_MUMPS_SUCCESS .and. associated(matrix))
      if (status /= IGRM_MUMPS_SUCCESS .or. .not. associated(matrix)) return

      call check('nelem=2, ptest=2, ptrial=1 gives a 91 by 91 system', &
                 matrix%x == SYSTEM_SIZE .and. matrix%y == SYSTEM_SIZE .and. &
                 size(right_hand_side) == SYSTEM_SIZE)
      call check('the reference small system has exactly 370 nonzeros', &
                 matrix%total_entries == EXPECTED_NNZ)

      allocate (dense(0:SYSTEM_SIZE - 1, 0:SYSTEM_SIZE - 1))
      call to_dense_matrix(matrix, dense)

      boundary_rows_ok = boundary_rows_are_identity(dense, ads_test%n, &
                                                    ads_trial%n)
      call check('all strong-Dirichlet boundary rows are identity rows', &
                 boundary_rows_ok)

      gram_ok = internal_gram_block_is_symmetric(dense, ads_test%n)
      call check('the internal G block is symmetric with positive diagonal', &
                 gram_ok)

      saddle_ok = internal_saddle_blocks_match(dense, ads_test%n, ads_trial%n)
      call check('internal off-diagonal blocks have the -B/B^T signs', &
                 saddle_ok)

      lower_rhs_ok = maxval(abs(right_hand_side(TEST_SIZE + 1:SYSTEM_SIZE))) &
                     <= STRUCTURE_TOLERANCE
      call check('the complete trial-space right-hand-side block is zero', &
                 lower_rhs_ok)

      deallocate (dense, right_hand_side)
      call clear_matrix(matrix)
      call check('assembled matrix cleanup nullifies its pointer', &
                 .not. associated(matrix))
   end subroutine test_assembled_system


   subroutine test_direct_solve
      type(IGRMMumpsStats) :: stats
      real(kind=8), allocatable :: solution(:, :)
      integer(kind=4) :: status
      logical :: finite_solution, zero_boundary, nonzero_interior

      call SolveIGRMMumps3D(ads_test, ads_trial, EPSILON_VALUE, BETA, &
                            SOLVER_TOLERANCE, positive_forcing, solution, &
                            stats, status)
      call check('real MUMPS factorization and solve succeeds', &
                 status == IGRM_MUMPS_SUCCESS)
      if (status /= IGRM_MUMPS_SUCCESS) return

      call check('solve reports the assembled size and nonzero count', &
                 stats%system_size == SYSTEM_SIZE .and. &
                 stats%nonzeros == EXPECTED_NNZ)
      call check('solve reports a finite small algebraic residual', &
                 ieee_is_finite(stats%residual_rms) .and. &
                 ieee_is_finite(stats%relative_residual) .and. &
                 stats%relative_residual <= SOLVER_TOLERANCE)

      call inspect_solution(solution, finite_solution, zero_boundary, &
                            nonzero_interior)
      call check('all returned trial coefficients are finite', finite_solution)
      call check('all six trial-space boundary faces remain zero', zero_boundary)
      call check('the solved trial field has a nonzero interior coefficient', &
                 nonzero_interior)

      if (allocated(solution)) deallocate (solution)
   end subroutine test_direct_solve


   logical function boundary_rows_are_identity(matrix, test_n, trial_n) &
      result(matches)
      real(kind=8), intent(in) :: matrix(0:, 0:)
      integer(kind=4), intent(in) :: test_n(3), trial_n(3)
      integer(kind=4) :: ix, iy, iz, row

      matches = .true.
      do iz = 0, test_n(3)
         do iy = 0, test_n(2)
            do ix = 0, test_n(1)
               if (.not. is_boundary(ix, iy, iz, test_n)) cycle
               row = linear_index(ix, iy, iz, test_n)
               matches = matches .and. row_is_identity(matrix, row)
            end do
         end do
      end do

      do iz = 0, trial_n(3)
         do iy = 0, trial_n(2)
            do ix = 0, trial_n(1)
               if (.not. is_boundary(ix, iy, iz, trial_n)) cycle
               row = TEST_SIZE + linear_index(ix, iy, iz, trial_n)
               matches = matches .and. row_is_identity(matrix, row)
            end do
         end do
      end do
   end function boundary_rows_are_identity


   logical function row_is_identity(matrix, row) result(matches)
      real(kind=8), intent(in) :: matrix(0:, 0:)
      integer(kind=4), intent(in) :: row
      integer(kind=4) :: column

      matches = abs(matrix(row, row) - 1.d0) <= STRUCTURE_TOLERANCE
      do column = 0, ubound(matrix, 2)
         if (column == row) cycle
         matches = matches .and. &
                   abs(matrix(row, column)) <= STRUCTURE_TOLERANCE
      end do
   end function row_is_identity


   logical function internal_gram_block_is_symmetric(matrix, test_n) &
      result(matches)
      real(kind=8), intent(in) :: matrix(0:, 0:)
      integer(kind=4), intent(in) :: test_n(3)
      integer(kind=4) :: ix, iy, iz, jx, jy, jz, row, column

      matches = .true.
      do iz = 1, test_n(3) - 1
         do iy = 1, test_n(2) - 1
            do ix = 1, test_n(1) - 1
               row = linear_index(ix, iy, iz, test_n)
               matches = matches .and. matrix(row, row) > 0.d0
               do jz = 1, test_n(3) - 1
                  do jy = 1, test_n(2) - 1
                     do jx = 1, test_n(1) - 1
                        column = linear_index(jx, jy, jz, test_n)
                        matches = matches .and. &
                           abs(matrix(row, column) - matrix(column, row)) &
                           <= STRUCTURE_TOLERANCE
                     end do
                  end do
               end do
            end do
         end do
      end do
   end function internal_gram_block_is_symmetric


   logical function internal_saddle_blocks_match(matrix, test_n, trial_n) &
      result(matches)
      real(kind=8), intent(in) :: matrix(0:, 0:)
      integer(kind=4), intent(in) :: test_n(3), trial_n(3)
      integer(kind=4) :: ix, iy, iz, test_row, trial_row
      logical :: has_positive_b, has_negative_b
      real(kind=8) :: lower_value

      trial_row = TEST_SIZE + linear_index(1, 1, 1, trial_n)
      matches = maxval(abs(matrix(trial_row, TEST_SIZE:SYSTEM_SIZE - 1))) &
                <= STRUCTURE_TOLERANCE
      has_positive_b = .false.
      has_negative_b = .false.

      do iz = 1, test_n(3) - 1
         do iy = 1, test_n(2) - 1
            do ix = 1, test_n(1) - 1
               test_row = linear_index(ix, iy, iz, test_n)
               lower_value = matrix(trial_row, test_row)
               matches = matches .and. &
                  abs(matrix(test_row, trial_row) + lower_value) &
                  <= STRUCTURE_TOLERANCE
               has_positive_b = has_positive_b .or. &
                                lower_value > STRUCTURE_TOLERANCE
               has_negative_b = has_negative_b .or. &
                                lower_value < -STRUCTURE_TOLERANCE
            end do
         end do
      end do

      matches = matches .and. has_positive_b .and. has_negative_b
   end function internal_saddle_blocks_match


   subroutine inspect_solution(solution, finite_values, zero_boundary, &
                               nonzero_interior)
      real(kind=8), intent(in) :: solution(:, :)
      logical, intent(out) :: finite_values, zero_boundary, nonzero_interior
      integer(kind=4) :: lx, ly, lz, gx, gy, gz, local_column
      real(kind=8) :: coefficient

      finite_values = all(ieee_is_finite(solution))
      zero_boundary = .true.
      nonzero_interior = .false.

      do lz = 1, ads_trial%s(3)
         gz = ads_trial%ibeg(3) + lz - 2
         do ly = 1, ads_trial%s(2)
            gy = ads_trial%ibeg(2) + ly - 2
            local_column = ly + (lz - 1)*ads_trial%s(2)
            do lx = 1, ads_trial%s(1)
               gx = ads_trial%ibeg(1) + lx - 2
               coefficient = solution(lx, local_column)
               if (is_boundary(gx, gy, gz, ads_trial%n)) then
                  zero_boundary = zero_boundary .and. &
                                  abs(coefficient) <= SOLVER_TOLERANCE
               else
                  nonzero_interior = nonzero_interior .or. &
                                     abs(coefficient) > SOLVER_TOLERANCE
               end if
            end do
         end do
      end do
   end subroutine inspect_solution


   pure integer(kind=4) function linear_index(ix, iy, iz, n) result(index_value)
      integer(kind=4), intent(in) :: ix, iy, iz, n(3)

      index_value = ix + (n(1) + 1)*(iy + (n(2) + 1)*iz)
   end function linear_index


   pure logical function is_boundary(ix, iy, iz, n) result(boundary)
      integer(kind=4), intent(in) :: ix, iy, iz, n(3)

      boundary = ix == 0 .or. ix == n(1) .or. &
                 iy == 0 .or. iy == n(2) .or. &
                 iz == 0 .or. iz == n(3)
   end function is_boundary

   subroutine check(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      checks = checks + 1
      if (.not. condition) then
         failures = failures + 1
         if (MYRANK == 0) write (*, '(A,A)') 'FAIL: ', trim(label)
      end if
   end subroutine check

end program test_igrm_eirksson_solver
