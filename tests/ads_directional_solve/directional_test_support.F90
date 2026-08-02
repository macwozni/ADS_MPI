module directional_test_support
   use Setup, only: ADS_Setup
   implicit none

   integer(kind=4) :: compute_calls = 0
   integer(kind=4) :: solve_calls = 0
   integer(kind=4) :: expected_p_test = 0
   integer(kind=4) :: expected_n_test = 0
   integer(kind=4) :: expected_nelem_test = 0
   integer(kind=4) :: expected_p_trial = 0
   integer(kind=4) :: expected_n_trial = 0
   integer(kind=4) :: expected_nelem_trial = 0
   integer(kind=4) :: expected_matrix_size = 0
   integer(kind=4) :: injected_solver_rank = -1
   integer(kind=4) :: injected_solver_status = 0
   integer(kind=4) :: injected_solver_rank_2 = -1
   integer(kind=4) :: injected_solver_status_2 = 0
   integer(kind=4) :: injected_matrix_size_delta = 0
   logical :: expected_equ = .true.
   logical :: compute_contract_ok = .true.
   logical :: solve_contract_ok = .true.
   real(kind=8) :: expected_mixA(4) = 0.d0
   real(kind=8) :: expected_mixB(4) = 0.d0
   real(kind=8) :: expected_mixBT(4) = 0.d0
   real(kind=8), allocatable :: expected_u_test(:)
   real(kind=8), allocatable :: expected_u_trial(:)
   real(kind=8), allocatable :: expected_packed_rhs(:, :)

contains

   subroutine configure_directional_spies(ads_test, ads_trial, axis, &
                                           mixA, mixB, mixBT, equ, packed_rhs)
      type(ADS_Setup), intent(in) :: ads_test, ads_trial
      integer(kind=4), intent(in) :: axis
      real(kind=8), intent(in) :: mixA(4), mixB(4), mixBT(4)
      logical, intent(in) :: equ
      real(kind=8), intent(in) :: packed_rhs(:, :)

      compute_calls = 0
      solve_calls = 0
      compute_contract_ok = .true.
      solve_contract_ok = .true.
      injected_solver_rank = -1
      injected_solver_status = 0
      injected_solver_rank_2 = -1
      injected_solver_status_2 = 0
      injected_matrix_size_delta = 0

      expected_p_test = ads_test%p(axis)
      expected_n_test = ads_test%n(axis)
      expected_nelem_test = ads_test%nelem(axis)
      expected_p_trial = ads_trial%p(axis)
      expected_n_trial = ads_trial%n(axis)
      expected_nelem_trial = ads_trial%nelem(axis)
      expected_matrix_size = size(packed_rhs, 1)
      expected_equ = equ
      expected_mixA = mixA
      expected_mixB = mixB
      expected_mixBT = mixBT

      if (allocated(expected_u_test)) deallocate(expected_u_test)
      if (allocated(expected_u_trial)) deallocate(expected_u_trial)
      if (allocated(expected_packed_rhs)) deallocate(expected_packed_rhs)

      select case (axis)
      case (1)
         allocate(expected_u_test(0:size(ads_test%Ux) - 1))
         allocate(expected_u_trial(0:size(ads_trial%Ux) - 1))
         expected_u_test = ads_test%Ux
         expected_u_trial = ads_trial%Ux
      case (2)
         allocate(expected_u_test(0:size(ads_test%Uy) - 1))
         allocate(expected_u_trial(0:size(ads_trial%Uy) - 1))
         expected_u_test = ads_test%Uy
         expected_u_trial = ads_trial%Uy
      case (3)
         allocate(expected_u_test(0:size(ads_test%Uz) - 1))
         allocate(expected_u_trial(0:size(ads_trial%Uz) - 1))
         expected_u_test = ads_test%Uz
         expected_u_trial = ads_trial%Uz
      case default
         stop 71
      end select

      allocate(expected_packed_rhs(size(packed_rhs, 1), size(packed_rhs, 2)))
      expected_packed_rhs = packed_rhs
   end subroutine configure_directional_spies


   subroutine inject_solver_failure(world_rank, status)
      integer(kind=4), intent(in) :: world_rank, status

      injected_solver_rank = world_rank
      injected_solver_status = status
   end subroutine inject_solver_failure


   subroutine inject_two_solver_failures(world_rank_1, status_1, &
                                         world_rank_2, status_2)
      integer(kind=4), intent(in) :: world_rank_1, status_1
      integer(kind=4), intent(in) :: world_rank_2, status_2

      injected_solver_rank = world_rank_1
      injected_solver_status = status_1
      injected_solver_rank_2 = world_rank_2
      injected_solver_status_2 = status_2
   end subroutine inject_two_solver_failures


   subroutine inject_matrix_size_mismatch(delta)
      integer(kind=4), intent(in) :: delta

      injected_matrix_size_delta = delta
   end subroutine inject_matrix_size_mismatch


   subroutine clear_solver_failure()
      injected_solver_rank = -1
      injected_solver_status = 0
      injected_solver_rank_2 = -1
      injected_solver_status_2 = 0
   end subroutine clear_solver_failure


   logical function exact_matrix(actual, expected) result(matches)
      real(kind=8), intent(in) :: actual(:, :), expected(:, :)

      matches = all(shape(actual) == shape(expected))
      if (.not. matches) return
      matches = all(actual == expected)
   end function exact_matrix

end module directional_test_support
