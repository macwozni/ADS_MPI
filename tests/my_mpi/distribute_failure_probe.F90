program distribute_failure_probe
   use mpi, only: MPI_SUCCESS, FAIL_IRECV, FAIL_ISEND, FAIL_WAITALL, &
      IRECV_SENTINEL, ISEND_SENTINEL, WAITALL_SENTINEL, &
      configure_mpi_failure, first_mpi_error, irecv_calls, isend_calls, &
      waitall_calls, transport_posts_after_failure
   use Setup, only: ADS_setup, ADS_compute_data
   use my_mpi, only: DistributeSpline
   implicit none

   real(kind=8), parameter :: R_SENTINEL = -777.0d0
   real(kind=8), parameter :: RECV_SENTINEL = -333.0d0
   integer(kind=4), parameter :: REMOTE_PEERS = 2
   type(ADS_setup) :: trial_space
   type(ADS_compute_data) :: data
   real(kind=8), dimension(3, 1) :: part
   character(len=16) :: mode_name
   integer(kind=4) :: mode, expected_status, workflow_status
   integer(kind=4) :: failures, checks

   call get_command_argument(1, mode_name)
   select case (trim(mode_name))
   case ('irecv')
      mode = FAIL_IRECV
      expected_status = IRECV_SENTINEL
   case ('isend')
      mode = FAIL_ISEND
      expected_status = ISEND_SENTINEL
   case ('waitall')
      mode = FAIL_WAITALL
      expected_status = WAITALL_SENTINEL
   case default
      write(*, '(A)') 'usage: distribute_failure_probe {irecv|isend|waitall}'
      stop 2
   end select

   call prepare_exchange(trial_space, data, part)
   call configure_mpi_failure(mode)
   workflow_status = MPI_SUCCESS

   call distribute_workflow(part, trial_space, data, workflow_status)

   failures = 0
   checks = 0
   call check('fault injector emitted the exact sentinel status', &
      first_mpi_error == expected_status, failures, checks)
#ifdef DISTRIBUTE_ERROR_STATUS_API
   call check('DistributeSpline exposes the first MPI error to its caller', &
      workflow_status == expected_status, failures, checks)
#endif
   call check('no later receive or send is posted after the first MPI error', &
      transport_posts_after_failure == 0, failures, checks)
   call check('the failed exchange does not publish any coefficient into R', &
      all(data%R == R_SENTINEL), failures, checks)

   select case (mode)
   case (FAIL_IRECV)
      call check('MPI_Irecv failure stops receive posting immediately', &
         irecv_calls == 1 .and. isend_calls == 0 .and. waitall_calls == 0, &
         failures, checks)
      call check('MPI_Irecv failure skips later receive-buffer updates', &
         all(data%halo_recv_buffer == RECV_SENTINEL), failures, checks)
   case (FAIL_ISEND)
      call check('MPI_Isend failure skips later sends', &
         irecv_calls == REMOTE_PEERS .and. isend_calls == 1, &
         failures, checks)
      call check('MPI_Isend failure skips the later local self-copy', &
         data%halo_recv_buffer(1) == RECV_SENTINEL, failures, checks)
   case (FAIL_WAITALL)
      call check('MPI_Waitall is reached once after all requests are posted', &
         irecv_calls == REMOTE_PEERS .and. isend_calls == REMOTE_PEERS .and. &
         waitall_calls == 1, failures, checks)
   end select

   if (failures == 0) then
      write(*, '(A,A,A,I0,A)') 'OK: ', trim(mode_name), ' (', checks, ' checks)'
   else
      write(*, '(A,A,A,I0,A,I0,A)') 'FAILED: ', trim(mode_name), ' (', &
         failures, ' of ', checks, ' checks)'
      stop 1
   end if

contains

   subroutine distribute_workflow(local_part, setup, compute_data, status)
      real(kind=8), dimension(:, :), intent(in) :: local_part
      type(ADS_setup), intent(in) :: setup
      type(ADS_compute_data), intent(inout) :: compute_data
      integer(kind=4), intent(inout) :: status

#ifdef DISTRIBUTE_ERROR_STATUS_API
      call DistributeSpline(local_part, setup, compute_data, status)
#else
      call DistributeSpline(local_part, setup, compute_data)
#endif
   end subroutine distribute_workflow


   subroutine prepare_exchange(setup, compute_data, local_part)
      type(ADS_setup), intent(out) :: setup
      type(ADS_compute_data), intent(out) :: compute_data
      real(kind=8), dimension(:, :), intent(out) :: local_part
      integer(kind=4) :: peer

      setup%ibeg = (/1, 1, 1/)
      setup%s = (/3, 1, 1/)
      local_part(:, 1) = (/11.0d0, 22.0d0, 33.0d0/)

      compute_data%halo_begin = (/0, 0, 0/)
      compute_data%halo_end = (/2, 0, 0/)
      allocate(compute_data%halo_send_begin(3, 3))
      allocate(compute_data%halo_send_end(3, 3))
      allocate(compute_data%halo_recv_begin(3, 3))
      allocate(compute_data%halo_recv_end(3, 3))
      allocate(compute_data%halo_send_count(3))
      allocate(compute_data%halo_send_displ(3))
      allocate(compute_data%halo_recv_count(3))
      allocate(compute_data%halo_recv_displ(3))
      allocate(compute_data%halo_send_buffer(3))
      allocate(compute_data%halo_recv_buffer(3))
      allocate(compute_data%halo_requests(4))
      allocate(compute_data%halo_statuses(3, 4))
      allocate(compute_data%R(3, 1, 1, 1))

      do peer = 1, 3
         compute_data%halo_send_begin(:, peer) = (/peer - 1, 0, 0/)
         compute_data%halo_send_end(:, peer) = (/peer - 1, 0, 0/)
         compute_data%halo_recv_begin(:, peer) = (/peer - 1, 0, 0/)
         compute_data%halo_recv_end(:, peer) = (/peer - 1, 0, 0/)
      end do
      compute_data%halo_send_count = 1
      compute_data%halo_recv_count = 1
      compute_data%halo_send_displ = (/0, 1, 2/)
      compute_data%halo_recv_displ = (/0, 1, 2/)
      compute_data%halo_send_buffer = -222.0d0
      compute_data%halo_recv_buffer = RECV_SENTINEL
      compute_data%halo_requests = -1
      compute_data%halo_statuses = -1
      compute_data%R = R_SENTINEL
   end subroutine prepare_exchange


   subroutine check(label, condition, failure_count, check_count)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      integer(kind=4), intent(inout) :: failure_count, check_count

      check_count = check_count + 1
      if (condition) then
         write(*, '(A,A)') 'PASS: ', trim(label)
      else
         write(*, '(A,A)') 'FAIL: ', trim(label)
         failure_count = failure_count + 1
      end if
   end subroutine check

end program distribute_failure_probe
