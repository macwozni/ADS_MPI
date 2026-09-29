module igrm_l2_benchmark_adapter

   use benchmark_contract, ONLY: BenchmarkAdapter, BenchmarkConfiguration, &
                                 BenchmarkState
   use benchmark_runtime, ONLY: AdvanceWithCallback, CleanupCommon, &
                                InitializeCommon, MeasureCommon, &
                                ProjectInitialCommon

   implicit none

   private
   public :: RegisterIgrmL2Adapter

contains

   subroutine RegisterIgrmL2Adapter(adapter)
      use manufactured_solution, ONLY: RegisterTemporalPolynomialCase
      type(BenchmarkAdapter), intent(out) :: adapter

      adapter%name = 'igrm_l2'
      adapter%capabilities = &
         'anisotropic test/trial degrees; repeated physical ADS steps'
      call RegisterTemporalPolynomialCase(adapter%exact_case)
      adapter%initialize => InitializeCommon
      adapter%project_initial => ProjectInitialCommon
      adapter%advance => AdvanceIgrmL2
      adapter%measure => MeasureCommon
      adapter%cleanup => CleanupCommon

   end subroutine RegisterIgrmL2Adapter

   subroutine AdvanceIgrmL2(config, state, step, status)
      use manufactured_solution, ONLY: ManufacturedRhsPoint
      type(BenchmarkConfiguration), intent(in) :: config
      type(BenchmarkState), intent(inout) :: state
      integer(kind=4), intent(in) :: step
      integer(kind=4), intent(out) :: status

      call AdvanceWithCallback(config, state, step, ManufacturedRhsPoint, status)

   end subroutine AdvanceIgrmL2

end module igrm_l2_benchmark_adapter
