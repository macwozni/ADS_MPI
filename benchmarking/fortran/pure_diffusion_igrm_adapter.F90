module pure_diffusion_igrm_benchmark_adapter

   use benchmark_contract, ONLY: BenchmarkAdapter, BenchmarkConfiguration, &
                                 BenchmarkState
   use benchmark_runtime, ONLY: AdvanceWithCallback, CleanupCommon, &
                                InitializeCommon, MeasureCommon, &
                                ProjectInitialCommon

   implicit none

   private
   public :: RegisterPureDiffusionIgrmAdapter

contains

   subroutine RegisterPureDiffusionIgrmAdapter(adapter)
      use manufactured_solution, ONLY: RegisterTemporalPolynomialCase
      type(BenchmarkAdapter), intent(out) :: adapter

      adapter%name = 'pure_diffusion_igrm'
      adapter%capabilities = &
         'anisotropic benchmark degrees; explicit mass projection'
      call RegisterTemporalPolynomialCase(adapter%exact_case)
      adapter%initialize => InitializeCommon
      adapter%project_initial => ProjectInitialCommon
      adapter%advance => AdvancePureDiffusionIgrm
      adapter%measure => MeasureCommon
      adapter%cleanup => CleanupCommon

   end subroutine RegisterPureDiffusionIgrmAdapter

   subroutine AdvancePureDiffusionIgrm(config, state, step, status)
      use manufactured_solution, ONLY: ManufacturedRhsPoint
      type(BenchmarkConfiguration), intent(in) :: config
      type(BenchmarkState), intent(inout) :: state
      integer(kind=4), intent(in) :: step
      integer(kind=4), intent(out) :: status

      call AdvanceWithCallback(config, state, step, ManufacturedRhsPoint, status)

   end subroutine AdvancePureDiffusionIgrm

end module pure_diffusion_igrm_benchmark_adapter
