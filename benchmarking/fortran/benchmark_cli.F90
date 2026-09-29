module benchmark_cli

   use benchmark_contract, ONLY: BenchmarkConfiguration

   implicit none

   private
   public :: ReadBenchmarkConfiguration
   public :: WriteBenchmarkUsage

contains

   subroutine ReadBenchmarkConfiguration(config, status)
      use, intrinsic :: ieee_arithmetic, ONLY: ieee_is_finite
      type(BenchmarkConfiguration), intent(out) :: config
      integer(kind=4), intent(out) :: status
      character(len=64) :: argument
      integer(kind=4) :: argument_count, axis, write_samples_integer

      status = 0
      argument_count = command_argument_count()
      if (argument_count /= 17 .and. argument_count /= 18) then
         status = 5
         return
      end if

      call get_command_argument(1, config%scheme)
      config%scheme = adjustl(config%scheme)
      call ReadRealArgument(2, config%final_time, status)
      call ReadIntegerArgument(3, config%steps, status)
      do axis = 1, 3
         call ReadIntegerArgument(3 + axis, config%nelem(axis), status)
         call ReadIntegerArgument(6 + axis, config%p_test(axis), status)
         call ReadIntegerArgument(9 + axis, config%p_trial(axis), status)
         call ReadIntegerArgument(12 + axis, config%process_grid(axis), status)
      end do
      call ReadIntegerArgument(16, config%sample_points, status)
      call ReadIntegerArgument(17, write_samples_integer, status)
      if (argument_count == 18) then
         call ReadStringArgument(18, config%exact_case, status)
      else
         config%exact_case = 'temporal-polynomial'
      end if
      if (status /= 0) return

      if (trim(config%scheme) /= 'dg' .and. &
          trim(config%scheme) /= 'pr' .and. &
          trim(config%scheme) /= 'be') status = 5
      if (.not. ieee_is_finite(config%final_time) .or. &
          config%final_time <= 0.d0) status = 5
      if (config%steps <= 0) status = 5
      if (any(config%nelem <= 0)) status = 5
      if (any(config%p_test <= 0) .or. any(config%p_test > 9)) status = 5
      if (any(config%p_trial <= 0) .or. any(config%p_trial > 9)) status = 5
      select case (trim(config%exact_case))
      case ('temporal-polynomial')
         if (any(config%p_trial < 3)) status = 5
      case ('spatial-cosine')
         continue
      case default
         status = 5
      end select
      if (any(config%p_test <= config%p_trial)) status = 5
      if (any(config%process_grid <= 0) .or. &
          any(config%process_grid > 128)) status = 5
      if (config%sample_points < 2 .or. config%sample_points > 257) status = 5
      if (write_samples_integer /= 0 .and. write_samples_integer /= 1) status = 5
      config%write_samples = write_samples_integer == 1
      if (status /= 0) return

      config%dt = config%final_time/real(config%steps, kind=8)
      if (.not. ieee_is_finite(config%dt) .or. config%dt <= 0.d0) status = 5

      ! Keep the local variable in the interface deliberately small and
      ! ensure command_argument implementations cannot silently truncate a
      ! future wider numeric field.
      call get_command_argument(2, argument)
      if (len_trim(argument) == len(argument)) status = 5

   end subroutine ReadBenchmarkConfiguration

   subroutine ReadIntegerArgument(position, value, status)
      integer(kind=4), intent(in) :: position
      integer(kind=4), intent(out) :: value
      integer(kind=4), intent(inout) :: status
      character(len=64) :: argument
      integer(kind=4) :: io_status

      if (status /= 0) return
      call get_command_argument(position, argument)
      if (len_trim(argument) == len(argument)) then
         status = 5
         return
      end if
      read(argument, *, iostat=io_status) value
      if (io_status /= 0) status = 5

   end subroutine ReadIntegerArgument

   subroutine ReadStringArgument(position, value, status)
      integer(kind=4), intent(in) :: position
      character(len=*), intent(out) :: value
      integer(kind=4), intent(inout) :: status
      character(len=64) :: argument
      integer(kind=4) :: argument_length, argument_status

      if (status /= 0) return
      call get_command_argument( &
         position, argument, length=argument_length, status=argument_status)
      if (argument_status /= 0 .or. argument_length <= 0 .or. &
          argument_length > len(value) .or. argument_length >= len(argument)) then
         status = 5
         return
      end if
      value = adjustl(argument(:argument_length))

   end subroutine ReadStringArgument

   subroutine ReadRealArgument(position, value, status)
      integer(kind=4), intent(in) :: position
      real(kind=8), intent(out) :: value
      integer(kind=4), intent(inout) :: status
      character(len=64) :: argument
      integer(kind=4) :: io_status

      if (status /= 0) return
      call get_command_argument(position, argument)
      if (len_trim(argument) == len(argument)) then
         status = 5
         return
      end if
      read(argument, *, iostat=io_status) value
      if (io_status /= 0) status = 5

   end subroutine ReadRealArgument

   subroutine WriteBenchmarkUsage
      use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT

      write(ERROR_UNIT, '(A)') &
         'usage: <scheme> <T> <steps> <nx> <ny> <nz> ' // &
         '<ptest-x> <ptest-y> <ptest-z> <ptrial-x> <ptrial-y> ' // &
         '<ptrial-z> <proc-x> <proc-y> <proc-z> ' // &
         '<sample-points> <write-samples:0|1> [exact-case]'

   end subroutine WriteBenchmarkUsage

end module benchmark_cli
