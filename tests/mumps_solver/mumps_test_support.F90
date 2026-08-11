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
   logical :: initialization_starts_instance = .true.
   logical :: contract_recorded = .false.
   integer(kind=4) :: recorded_comm = 0
   integer(kind=4) :: recorded_sym = 0
   integer(kind=4) :: recorded_par = 0
   integer(kind=4) :: recorded_n = 0
   integer(kind=4) :: recorded_nz = 0
   integer(kind=4) :: recorded_icntl(40) = 0
   integer(kind=4) :: recorded_irn(64) = 0
   integer(kind=4) :: recorded_jcn(64) = 0
   real(kind=8) :: recorded_a(64) = 0.d0

contains

   subroutine configure_mumps_stub(job, code, occurrence, use_infog, finalize_code, &
                                   start_instance)
      integer(kind=4), intent(in) :: job, code
      integer(kind=4), intent(in), optional :: occurrence, finalize_code
      logical, intent(in), optional :: use_infog, start_instance

      recorded_jobs = 0
      recorded_job_count = 0
      job_occurrences = 0
      failure_job = job
      failure_code = code
      failure_occurrence = 1
      finalization_code = 0
      failure_uses_infog = .false.
      initialization_starts_instance = .true.
      contract_recorded = .false.
      recorded_comm = 0
      recorded_sym = 0
      recorded_par = 0
      recorded_n = 0
      recorded_nz = 0
      recorded_icntl = 0
      recorded_irn = 0
      recorded_jcn = 0
      recorded_a = 0.d0
      if (present(occurrence)) failure_occurrence = occurrence
      if (present(use_infog)) failure_uses_infog = use_infog
      if (present(finalize_code)) finalization_code = finalize_code
      if (present(start_instance)) initialization_starts_instance = start_instance
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


   subroutine record_mumps_contract(comm, sym, par, n, nz, irn, jcn, values, icntl)
      integer(kind=4), intent(in) :: comm, sym, par, n, nz
      integer(kind=4), intent(in) :: irn(:), jcn(:), icntl(:)
      real(kind=8), intent(in) :: values(:)
      integer(kind=4) :: entry_count, control_count

      contract_recorded = .true.
      recorded_comm = comm
      recorded_sym = sym
      recorded_par = par
      recorded_n = n
      recorded_nz = nz
      entry_count = min(nz, size(recorded_irn), size(irn), size(jcn), size(values))
      if (entry_count > 0) then
         recorded_irn(1:entry_count) = irn(1:entry_count)
         recorded_jcn(1:entry_count) = jcn(1:entry_count)
         recorded_a(1:entry_count) = values(1:entry_count)
      end if
      control_count = min(size(recorded_icntl), size(icntl))
      recorded_icntl(1:control_count) = icntl(1:control_count)
   end subroutine record_mumps_contract


   logical function expected_jobs(actual) result(matches)
      integer(kind=4), intent(in) :: actual(:)

      matches = recorded_job_count == size(actual)
      if (.not. matches) return
      matches = all(recorded_jobs(1:recorded_job_count) == actual)
   end function expected_jobs

end module mumps_test_support
