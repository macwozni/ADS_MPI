module mumps_test_support

   implicit none

   integer(kind=4), parameter :: MAX_RECORDED_JOBS = 32
   integer(kind=4) :: recorded_jobs(MAX_RECORDED_JOBS) = 0
   integer(kind=4) :: recorded_job_count = 0
   integer(kind=4) :: failure_job = 999
   integer(kind=4) :: failure_occurrence = 1
   integer(kind=4) :: failure_code = 0
   integer(kind=4) :: finalization_code = 0
   integer(kind=4) :: job_occurrences(-2:3) = 0
   logical :: failure_uses_infog = .false.

contains

   subroutine configure_mumps_stub(job, code, occurrence, use_infog, finalize_code)
      integer(kind=4), intent(in) :: job, code
      integer(kind=4), intent(in), optional :: occurrence, finalize_code
      logical, intent(in), optional :: use_infog

      recorded_jobs = 0
      recorded_job_count = 0
      job_occurrences = 0
      failure_job = job
      failure_code = code
      failure_occurrence = 1
      finalization_code = 0
      failure_uses_infog = .false.
      if (present(occurrence)) failure_occurrence = occurrence
      if (present(use_infog)) failure_uses_infog = use_infog
      if (present(finalize_code)) finalization_code = finalize_code
   end subroutine configure_mumps_stub


   subroutine record_mumps_job(job, occurrence)
      integer(kind=4), intent(in) :: job
      integer(kind=4), intent(out) :: occurrence

      recorded_job_count = recorded_job_count + 1
      if (recorded_job_count <= MAX_RECORDED_JOBS) then
         recorded_jobs(recorded_job_count) = job
      end if
      if (job >= lbound(job_occurrences, 1) .and. &
          job <= ubound(job_occurrences, 1)) then
         job_occurrences(job) = job_occurrences(job) + 1
         occurrence = job_occurrences(job)
      else
         occurrence = 0
      end if
   end subroutine record_mumps_job


   logical function expected_jobs(actual) result(matches)
      integer(kind=4), intent(in) :: actual(:)

      matches = recorded_job_count == size(actual)
      if (.not. matches) return
      matches = all(recorded_jobs(1:recorded_job_count) == actual)
   end function expected_jobs

end module mumps_test_support
