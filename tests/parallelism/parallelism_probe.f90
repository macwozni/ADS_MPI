program parallelism_probe
   use parallelism, only: InitializeParallelism, Cleanup_Parallelism, &
      AbortOnError
   implicit none

   integer(kind=4), dimension(3) :: process_grid
   integer(kind=4) :: ierr, injected_status

   call read_process_grid(process_grid)
   call InitializeParallelism(process_grid(1), process_grid(2), process_grid(3), ierr)
   if (command_argument_count() == 4) then
      call read_status(injected_status)
      call AbortOnError(injected_status, 'parallelism probe')
   end if
   call Cleanup_Parallelism(ierr)
   write (*, '(A)') 'UNEXPECTED SUCCESS'

contains

   subroutine read_process_grid(grid)
      integer(kind=4), dimension(3), intent(out) :: grid
      character(len=32) :: argument
      integer(kind=4) :: i, read_status

      if (command_argument_count() /= 3 .and. &
          command_argument_count() /= 4) stop 2

      do i = 1, 3
         call get_command_argument(i, argument)
         read (argument, *, iostat=read_status) grid(i)
         if (read_status /= 0) stop 2
      end do
   end subroutine read_process_grid


   subroutine read_status(status)
      integer(kind=4), intent(out) :: status
      character(len=32) :: argument
      integer(kind=4) :: read_status_code

      call get_command_argument(4, argument)
      read (argument, *, iostat=read_status_code) status
      if (read_status_code /= 0) stop 2
   end subroutine read_status

end program parallelism_probe
