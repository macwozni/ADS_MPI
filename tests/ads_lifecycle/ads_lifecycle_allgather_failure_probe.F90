program ads_lifecycle_allgather_failure_probe
   use mpi, only: ALLGATHER_FAILURE_SENTINEL, allgather_calls, &
                  reset_mpi_fault_stub
   use Setup, only: ADS_Setup, ADS_compute_data
   use ads_lifecycle, only: initialize, Cleanup_data, Cleanup_ADS
   implicit none

   type(ADS_Setup) :: test_space, trial_space
   type(ADS_compute_data) :: data
   integer(kind=4), parameter :: nelem(3) = (/2, 2, 2/)
   integer(kind=4), parameter :: degree(3) = (/1, 1, 1/)
   integer(kind=4), parameter :: continuity(3) = (/0, 0, 0/)
   integer(kind=4) :: checks, failures, initialize_ierr, cleanup_ierr
   logical :: cleanup_ok

   checks = 0
   failures = 0
   call reset_mpi_fault_stub()

   call initialize(nelem, degree, degree, continuity, test_space, &
                   trial_space, data, initialize_ierr)

   call assert_true('initialize preserves the MPI_Allgather failure code', &
                    initialize_ierr == ALLGATHER_FAILURE_SENTINEL)
   call assert_true('AllocateADSdata attempts MPI_Allgather exactly once', &
                    allgather_calls == 1)
   call assert_true('failed MPI_Allgather does not publish a halo plan', &
                    halo_plan_is_absent(data))

   call Cleanup_data(data, cleanup_ierr)
   cleanup_ok = cleanup_ierr == 0 .and. data_allocations_are_absent(data)
   call Cleanup_ADS(test_space, cleanup_ierr)
   cleanup_ok = cleanup_ok .and. cleanup_ierr == 0 .and. &
                setup_allocations_are_absent(test_space)
   call Cleanup_ADS(trial_space, cleanup_ierr)
   cleanup_ok = cleanup_ok .and. cleanup_ierr == 0 .and. &
                setup_allocations_are_absent(trial_space)
   call assert_true('partial initialization is safe to clean up', cleanup_ok)

   call Cleanup_data(data, cleanup_ierr)
   cleanup_ok = cleanup_ierr == 0 .and. data_allocations_are_absent(data)
   call Cleanup_ADS(test_space, cleanup_ierr)
   cleanup_ok = cleanup_ok .and. cleanup_ierr == 0
   call Cleanup_ADS(trial_space, cleanup_ierr)
   cleanup_ok = cleanup_ok .and. cleanup_ierr == 0
   call assert_true('cleanup after the injected failure is repeatable', &
                    cleanup_ok)

   if (failures == 0) then
      write (*, '(A,I0,A)') 'OK (', checks, &
                            ' ADS lifecycle MPI failure checks)'
   else
      write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                 ' ADS lifecycle MPI failure checks)'
      stop 1
   end if

contains

   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      checks = checks + 1
      if (condition) then
         write (*, '(A,A)') 'PASS ', trim(label)
      else
         failures = failures + 1
         write (*, '(A,A)') 'FAIL ', trim(label)
      end if
   end subroutine assert_true


   logical function halo_plan_is_absent(value) result(absent)
      type(ADS_compute_data), intent(in) :: value

      absent = .not. allocated(value%halo_send_begin) .and. &
               .not. allocated(value%halo_send_end) .and. &
               .not. allocated(value%halo_recv_begin) .and. &
               .not. allocated(value%halo_recv_end) .and. &
               .not. allocated(value%halo_send_count) .and. &
               .not. allocated(value%halo_send_displ) .and. &
               .not. allocated(value%halo_recv_count) .and. &
               .not. allocated(value%halo_recv_displ) .and. &
               .not. allocated(value%halo_send_buffer) .and. &
               .not. allocated(value%halo_recv_buffer) .and. &
               .not. allocated(value%halo_requests) .and. &
               .not. allocated(value%halo_statuses) .and. &
               .not. allocated(value%R)
   end function halo_plan_is_absent


   logical function data_allocations_are_absent(value) result(absent)
      type(ADS_compute_data), intent(in) :: value

      absent = halo_plan_is_absent(value) .and. &
               .not. allocated(value%F) .and. &
               .not. allocated(value%F2) .and. &
               .not. allocated(value%F3) .and. &
               .not. allocated(value%FF) .and. &
               .not. allocated(value%Ft) .and. &
               .not. allocated(value%Ft2) .and. &
               .not. allocated(value%Ft3) .and. &
               .not. allocated(value%FFt) .and. &
               .not. allocated(value%Un) .and. &
               .not. allocated(value%Un13) .and. &
               .not. allocated(value%Un23) .and. &
               .not. allocated(value%dUn) .and. &
               .not. allocated(value%dUn0) .and. &
               .not. allocated(value%dUn13) .and. &
               .not. allocated(value%dUn23)
   end function data_allocations_are_absent


   logical function setup_allocations_are_absent(value) result(absent)
      type(ADS_Setup), intent(in) :: value

      absent = .not. allocated(value%Ux) .and. &
               .not. allocated(value%Uy) .and. &
               .not. allocated(value%Uz) .and. &
               .not. allocated(value%dimensionsX) .and. &
               .not. allocated(value%dimensionsY) .and. &
               .not. allocated(value%dimensionsZ) .and. &
               .not. allocated(value%shiftsX) .and. &
               .not. allocated(value%shiftsY) .and. &
               .not. allocated(value%shiftsZ) .and. &
               .not. allocated(value%IPIVx) .and. &
               .not. allocated(value%IPIVy) .and. &
               .not. allocated(value%IPIVz) .and. &
               .not. allocated(value%Ox) .and. &
               .not. allocated(value%Oy) .and. &
               .not. allocated(value%Oz) .and. &
               .not. allocated(value%Jx) .and. &
               .not. allocated(value%Jy) .and. &
               .not. allocated(value%Jz) .and. &
               .not. allocated(value%Xx) .and. &
               .not. allocated(value%Xy) .and. &
               .not. allocated(value%Xz) .and. &
               .not. allocated(value%NNx) .and. &
               .not. allocated(value%NNy) .and. &
               .not. allocated(value%NNz) .and. &
               .not. allocated(value%Wx) .and. &
               .not. allocated(value%Wy) .and. &
               .not. allocated(value%Wz)
   end function setup_allocations_are_absent

end program ads_lifecycle_allgather_failure_probe
