module igrm_heat_benchmark_adapter

   use benchmark_contract, ONLY: BenchmarkAdapter, BenchmarkConfiguration, &
                                 BenchmarkState
   use benchmark_runtime, ONLY: AdvanceWithCallback, CleanupCommon, &
                                InitializeCommon, MeasureCommon, &
                                ProjectInitialCommon

   implicit none

   private
   public :: RegisterIgrmHeatAdapter

contains

   subroutine RegisterIgrmHeatAdapter(adapter)
      use manufactured_solution, ONLY: RegisterTemporalPolynomialCase
      type(BenchmarkAdapter), intent(out) :: adapter

      adapter%name = 'igrm_heat'
      adapter%capabilities = &
         'heat full-RHS callback; mass projection; no VTI in timed steps'
      call RegisterTemporalPolynomialCase(adapter%exact_case)
      adapter%initialize => InitializeCommon
      adapter%project_initial => ProjectInitialCommon
      adapter%advance => AdvanceIgrmHeat
      adapter%measure => MeasureCommon
      adapter%cleanup => CleanupCommon

   end subroutine RegisterIgrmHeatAdapter

   subroutine AdvanceIgrmHeat(config, state, step, status)
      use input_data, ONLY: t
      type(BenchmarkConfiguration), intent(in) :: config
      type(BenchmarkState), intent(inout) :: state
      integer(kind=4), intent(in) :: step
      integer(kind=4), intent(out) :: status

      ! The production heat callback uses t only to distinguish its initial
      ! projection branch.  Stage time for the manufactured source is set by
      ! HeatManufacturedRhsPoint below.
      t = real(step, kind=8)*config%dt
      call AdvanceWithCallback( &
         config, state, step, HeatManufacturedRhsPoint, status)

   end subroutine AdvanceIgrmHeat

   subroutine HeatManufacturedRhsPoint( &
      ads, point, quadrature, element, local_basis, gradient, history, &
      un11, un13, un23, ads_data, jacobian, weight, direction, substep, &
      alpha_step, forcing_cb, value)
      use Interfaces, ONLY: forcing_fun
      use manufactured_solution, ONLY: ManufacturedScalarSource, &
                                         SetCallbackSourceTime
      use RHS_fun, ONLY: heat_igrm_rhs_point
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      type(ADS_Setup), intent(in) :: ads
      real(kind=8), dimension(3), intent(in) :: point, gradient
      integer(kind=4), dimension(3), intent(in) :: quadrature, element
      integer(kind=4), dimension(3), intent(in) :: local_basis, direction
      integer(kind=4), intent(in) :: history, substep
      real(kind=8), intent(in) :: un11, un13, un23, jacobian, weight
      type(ADS_compute_data), intent(in) :: ads_data
      real(kind=8), dimension(7, 3), intent(in) :: alpha_step
      procedure(forcing_fun) :: forcing_cb
      real(kind=8), intent(out) :: value

      call SetCallbackSourceTime(substep)
      call heat_igrm_rhs_point( &
         ads, point, quadrature, element, local_basis, gradient, history, &
         un11, un13, un23, ads_data, jacobian, weight, direction, substep, &
         alpha_step, ManufacturedScalarSource, value)

   end subroutine HeatManufacturedRhsPoint

end module igrm_heat_benchmark_adapter
