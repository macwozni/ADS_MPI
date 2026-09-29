module benchmark_harness

   use benchmark_contract, ONLY: BenchmarkAdapter, BenchmarkConfiguration, &
                                 BenchmarkMetrics, BenchmarkState, &
                                 ValidateAdapter

   implicit none

   private
   public :: RunManufacturedBenchmark

contains

   subroutine RunManufacturedBenchmark(adapter)
      use benchmark_cli, ONLY: ReadBenchmarkConfiguration, WriteBenchmarkUsage
      use manufactured_solution, ONLY: ActivateManufacturedCase
      use parallelism, ONLY: MYRANK
      use, intrinsic :: ieee_arithmetic, ONLY: ieee_is_finite
      use mpi
      type(BenchmarkAdapter), intent(in) :: adapter
      type(BenchmarkConfiguration) :: config
      type(BenchmarkMetrics) :: final_metrics, initial_metrics
      type(BenchmarkState) :: state
      integer(kind=4) :: cleanup_status, mpi_ierr
      integer(kind=4) :: status, step
      logical :: rank_zero
      real(kind=8) :: actual_final_time, expected_initial_norm
      real(kind=8) :: local_elapsed, start_time, time_tolerance

      call ValidateAdapter(adapter, status)
      if (status /= 0) then
         write(*, '(A)') 'invalid benchmark adapter registration'
         stop 5
      end if
      call ActivateManufacturedCase(adapter%exact_case, status)
      if (status /= 0) then
         write(*, '(A)') 'invalid manufactured-case registration'
         stop 5
      end if
      call ReadBenchmarkConfiguration(config, status)
      if (status /= 0) then
         call WriteBenchmarkUsage
         stop 5
      end if

      call adapter%initialize(config, state, status)
      call RequireCollectiveSuccess( &
         adapter, state, status, 'benchmark initialization failed')
      rank_zero = MYRANK == 0

      call adapter%project_initial(config, state, status)
      call RequireCollectiveSuccess( &
         adapter, state, status, 'initial projection failed')
      call adapter%measure(config, state, 0.d0, .false., &
                           initial_metrics, status)
      call RequireCollectiveSuccess( &
         adapter, state, status, 'initial-state measurement failed')

      expected_initial_norm = adapter%exact_case%initial_l2_norm
      status = 0
      if (.not. MetricsAreFinite(initial_metrics) .or. &
          initial_metrics%l2_error > 1.d-10 .or. &
          initial_metrics%linf_error > 1.d-10 .or. &
          abs(initial_metrics%solution_l2_norm - expected_initial_norm) > &
             1.d-10) status = 7
      call RequireCollectiveSuccess( &
         adapter, state, status, &
         'initial manufactured state failed its oracle')

      call MPI_Barrier(MPI_COMM_WORLD, mpi_ierr)
      if (mpi_ierr /= 0) call FailWithCleanup( &
         adapter, state, mpi_ierr, 'pre-step MPI barrier failed')
      start_time = MPI_Wtime()
      do step = 1, config%steps
         call adapter%advance(config, state, step, status)
         call RequireCollectiveSuccess( &
            adapter, state, status, 'physical ADS step failed')
      end do
      call MPI_Barrier(MPI_COMM_WORLD, mpi_ierr)
      if (mpi_ierr /= 0) call FailWithCleanup( &
         adapter, state, mpi_ierr, 'post-step MPI barrier failed')
      local_elapsed = MPI_Wtime() - start_time
      call MPI_Allreduce(local_elapsed, state%step_wall_seconds, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, &
                         mpi_ierr)
      if (mpi_ierr /= 0) call FailWithCleanup( &
         adapter, state, mpi_ierr, 'step-timing reduction failed')

      actual_final_time = real(config%steps, kind=8)*config%dt
      time_tolerance = 64.d0*epsilon(1.d0)*max(1.d0, config%final_time)
      status = 0
      if (abs(actual_final_time - config%final_time) > time_tolerance .or. &
          abs(state%ads_data%t - actual_final_time) > time_tolerance) status = 7
      call RequireCollectiveSuccess( &
         adapter, state, status, 'physical loop did not terminate at T')

      call adapter%measure(config, state, actual_final_time, &
                           config%write_samples, final_metrics, status)
      call RequireCollectiveSuccess( &
         adapter, state, status, 'final-state measurement failed')
      status = 0
      if (.not. MetricsAreFinite(final_metrics) .or. &
          final_metrics%l2_error < 0.d0 .or. &
          final_metrics%linf_error < 0.d0 .or. &
          final_metrics%solution_l2_norm < 0.d0 .or. &
          .not. ieee_is_finite(state%step_wall_seconds) .or. &
          state%step_wall_seconds < 0.d0) status = 7
      call RequireCollectiveSuccess( &
         adapter, state, status, 'non-finite manufactured result')

      call adapter%cleanup(state, cleanup_status)
      if (cleanup_status /= 0) then
         if (rank_zero) write(*, '(A,I0)') 'benchmark cleanup failed: ', &
                                             cleanup_status
         stop 1
      end if
      if (rank_zero) call WriteResult( &
         adapter, config, initial_metrics, final_metrics, actual_final_time, &
         state%step_wall_seconds)

   end subroutine RunManufacturedBenchmark

   logical function MetricsAreFinite(metrics)
      use, intrinsic :: ieee_arithmetic, ONLY: ieee_is_finite
      type(BenchmarkMetrics), intent(in) :: metrics

      MetricsAreFinite = &
         ieee_is_finite(metrics%l2_error) .and. &
         ieee_is_finite(metrics%linf_error) .and. &
         ieee_is_finite(metrics%solution_l2_norm) .and. &
         ieee_is_finite(metrics%field_checksum)

   end function MetricsAreFinite

   subroutine RequireCollectiveSuccess(adapter, state, local_status, message)
      use mpi
      type(BenchmarkAdapter), intent(in) :: adapter
      type(BenchmarkState), intent(inout) :: state
      integer(kind=4), intent(in) :: local_status
      character(len=*), intent(in) :: message
      integer(kind=4) :: any_failure, local_failure, mpi_ierr

      if (.not. state%parallel_active) then
         if (local_status /= 0) call FailWithCleanup( &
            adapter, state, local_status, message)
         return
      end if

      local_failure = 0
      if (local_status /= 0) local_failure = 1
      call MPI_Allreduce(local_failure, any_failure, 1, MPI_INTEGER, MPI_MAX, &
                         MPI_COMM_WORLD, mpi_ierr)
      if (mpi_ierr /= 0) call FailWithCleanup( &
         adapter, state, mpi_ierr, 'benchmark status reduction failed')
      if (any_failure /= 0) call FailWithCleanup( &
         adapter, state, 1, message)

   end subroutine RequireCollectiveSuccess

   subroutine FailWithCleanup(adapter, state, failure_status, message)
      use parallelism, ONLY: MYRANK
      type(BenchmarkAdapter), intent(in) :: adapter
      type(BenchmarkState), intent(inout) :: state
      integer(kind=4), intent(in) :: failure_status
      character(len=*), intent(in) :: message
      integer(kind=4) :: cleanup_status
      logical :: report

      report = .true.
      if (state%parallel_active) report = MYRANK == 0
      call adapter%cleanup(state, cleanup_status)
      if (report) then
         write(*, '(A,A,A,A,I0)') trim(message), '; adapter=', &
                                  trim(adapter%name), '; status=', failure_status
         if (cleanup_status /= 0) write(*, '(A,I0)') &
            'additional cleanup status=', cleanup_status
      end if
      stop 1

   end subroutine FailWithCleanup

   subroutine WriteResult(adapter, config, initial_metrics, final_metrics, &
                          actual_final_time, step_wall_seconds)
      type(BenchmarkAdapter), intent(in) :: adapter
      type(BenchmarkConfiguration), intent(in) :: config
      type(BenchmarkMetrics), intent(in) :: initial_metrics, final_metrics
      real(kind=8), intent(in) :: actual_final_time, step_wall_seconds
      character(len=4096) :: record
      character(len=5) :: wrote_samples

      if (config%write_samples) then
         wrote_samples = 'true '
      else
         wrote_samples = 'false'
      end if
      record = &
         'ADS_BENCHMARK_RESULT {' // &
         '"schema_version":1,' // &
         '"kind":"ads-manufactured-transient-result",' // &
         '"exact_case":"' // trim(adapter%exact_case%name) // '",' // &
         '"problem":"' // trim(adapter%name) // '",' // &
         '"scheme":"' // trim(config%scheme) // '",' // &
         '"requested_final_time":' // trim(JsonReal(config%final_time)) // ',' // &
         '"actual_final_time":' // trim(JsonReal(actual_final_time)) // ',' // &
         '"time_step":' // trim(JsonReal(config%dt)) // ',' // &
         '"steps":' // trim(JsonInteger(config%steps)) // ',' // &
         '"initial_l2_error":' // &
            trim(JsonReal(initial_metrics%l2_error)) // ',' // &
         '"initial_linf_error":' // &
            trim(JsonReal(initial_metrics%linf_error)) // ',' // &
         '"l2_error":' // trim(JsonReal(final_metrics%l2_error)) // ',' // &
         '"linf_error":' // trim(JsonReal(final_metrics%linf_error)) // ',' // &
         '"solution_l2_norm":' // &
            trim(JsonReal(final_metrics%solution_l2_norm)) // ',' // &
         '"field_checksum":' // &
            trim(JsonReal(final_metrics%field_checksum)) // ',' // &
         '"sample_points_per_axis":' // &
            trim(JsonInteger(config%sample_points)) // ',' // &
         '"field_samples_written":' // trim(wrote_samples) // ',' // &
         '"physical_step_wall_seconds":' // &
            trim(JsonReal(step_wall_seconds)) // ',' // &
         '"solver_status":0}'
      write(*, '(A)') trim(record)

   end subroutine WriteResult

   function JsonReal(value) result(text)
      real(kind=8), intent(in) :: value
      character(len=32) :: text

      write(text, '(ES24.16E3)') value
      text = adjustl(text)

   end function JsonReal

   function JsonInteger(value) result(text)
      integer(kind=4), intent(in) :: value
      character(len=16) :: text

      write(text, '(I0)') value
      text = adjustl(text)

   end function JsonInteger

end module benchmark_harness
