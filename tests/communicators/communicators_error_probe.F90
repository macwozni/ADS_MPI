program communicators_error_probe
   use communicators, only: COMMX, COMMY, COMMZ, Cleanup_Communicators, &
      CreateCommunicators
   use mpi, only: MPI_COMM_NULL, comm_free_calls, configure_failure, &
      group_free_calls
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
end program communicators_error_probe
