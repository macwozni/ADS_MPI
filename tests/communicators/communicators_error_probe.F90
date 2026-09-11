program communicators_error_probe
   use communicators, only: COMMX, COMMY, COMMZ, Cleanup_Communicators, &
      CreateCommunicators
   use mpi, only: BARRIER_1_SENTINEL, BARRIER_2_SENTINEL, &
      BARRIER_3_SENTINEL, MPI_COMM_NULL, barrier_calls, comm_create_calls, &
      comm_free_calls, configure_failure, group_free_calls, group_incl_calls
   implicit none

   character(len=32) :: mode
   integer(kind=4) :: ierr

   if (command_argument_count() /= 1) stop 90
   call get_command_argument(1, mode)

   select case (trim(mode))
   case ('comm-group')
      call configure_failure(1)
   case ('group-z')
      call configure_failure(2)
   case ('group-y')
      call configure_failure(3)
   case ('group-x')
      call configure_failure(4)
   case ('world-group-free')
      call configure_failure(5)
   case ('comm-z')
      call configure_failure(6)
   case ('comm-y')
      call configure_failure(7)
   case ('comm-x')
      call configure_failure(8)
   case ('barrier-1')
      call check_barrier_failure(9, BARRIER_1_SENTINEL, 1, 0, 0)
      stop
   case ('barrier-2')
      call check_barrier_failure(10, BARRIER_2_SENTINEL, 2, 26, 0)
      stop
   case ('barrier-3')
      call check_barrier_failure(11, BARRIER_3_SENTINEL, 3, 26, 26)
      stop
   case ('cleanup')
      call configure_failure(0)
      call CreateCommunicators(ierr)
      if (ierr /= 0) stop 91

      call configure_failure(20)
      call Cleanup_Communicators(ierr)
      if (ierr /= 201) stop 92
      if (comm_free_calls /= 26) stop 93
      if (group_free_calls /= 26) stop 94
      if (COMMX /= MPI_COMM_NULL .or. COMMY /= MPI_COMM_NULL .or. &
          COMMZ /= MPI_COMM_NULL) stop 95
      write (*, '(A)') 'SUCCESS cleanup retained first error and continued'
      stop
   case default
      stop 90
   end select

   call CreateCommunicators(ierr)
   write (*, '(A)') 'UNEXPECTED SUCCESS'
   stop 98

contains

   subroutine check_barrier_failure(failure, sentinel, expected_barriers, &
      expected_groups, expected_comms)
      integer(kind=4), intent(in) :: failure, sentinel, expected_barriers
      integer(kind=4), intent(in) :: expected_groups, expected_comms
      integer(kind=4) :: cleanup_ierr, comm_frees_after_cleanup
      integer(kind=4) :: failures, group_frees_after_cleanup

      failures = 0
      call configure_failure(failure)
      call CreateCommunicators(ierr)

      call check(ierr == sentinel, &
         'CreateCommunicators returns the exact MPI_Barrier error', failures)
      call check(barrier_calls == expected_barriers, &
         'no later MPI_Barrier is entered', failures)
      call check(group_incl_calls == expected_groups, &
         'no later group-creation phase is entered', failures)
      call check(comm_create_calls == expected_comms, &
         'no later communicator-creation phase is entered', failures)
      call check(COMMX == MPI_COMM_NULL .and. COMMY == MPI_COMM_NULL .and. &
         COMMZ == MPI_COMM_NULL, &
         'failed creation does not publish local communicators', failures)

      call Cleanup_Communicators(cleanup_ierr)
      call check(cleanup_ierr == 0, &
         'partial communicator state can be cleaned', failures)
      call check(COMMX == MPI_COMM_NULL .and. COMMY == MPI_COMM_NULL .and. &
         COMMZ == MPI_COMM_NULL, &
         'cleanup nulls public handles after partial creation', failures)
      call check(group_free_calls == expected_groups + 1 .and. &
         comm_free_calls == expected_comms, &
         'cleanup releases exactly the objects created before the barrier', &
         failures)

      comm_frees_after_cleanup = comm_free_calls
      group_frees_after_cleanup = group_free_calls
      call Cleanup_Communicators(cleanup_ierr)
      call check(cleanup_ierr == 0 .and. &
         comm_free_calls == comm_frees_after_cleanup .and. &
         group_free_calls == group_frees_after_cleanup, &
         'cleanup after partial creation is repeatable', failures)

      if (failures /= 0) then
         write (*, '(A,I0,A)') 'FAILED (', failures, &
            ' MPI_Barrier fault-injection checks)'
         stop 1
      end if

      write (*, '(A,I0,A)') 'SUCCESS barrier ', expected_barriers, &
         ' failure stopped communicator creation'
   end subroutine check_barrier_failure


   subroutine check(condition, label, failures)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      integer(kind=4), intent(inout) :: failures

      if (condition) then
         write (*, '(A,A)') 'PASS ', label
      else
         write (*, '(A,A)') 'FAIL ', label
         failures = failures + 1
      end if
   end subroutine check
end program communicators_error_probe
