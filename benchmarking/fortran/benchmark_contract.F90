! Public contract between the manufactured benchmark harness and a problem
! adapter.  The harness owns the lifecycle; adapters select the concrete ADS
! initialization, step/callback path, measurement, and cleanup procedures.
module benchmark_contract

   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use time_scheme, ONLY: TimeScheme3D
   use, intrinsic :: ieee_arithmetic, ONLY: ieee_is_finite

   implicit none

   private
   public :: BenchmarkAdapter
   public :: BenchmarkConfiguration
   public :: BenchmarkMetrics
   public :: BenchmarkState
   public :: ManufacturedCase
   public :: ValidateAdapter

   abstract interface
      pure function manufactured_value_fun(physical_time, point) result(value)
         real(kind=8), intent(in) :: physical_time
         real(kind=8), dimension(3), intent(in) :: point
         real(kind=8) :: value
      end function manufactured_value_fun

      pure function manufactured_source_fun(physical_time, point) result(value)
         real(kind=8), intent(in) :: physical_time
         real(kind=8), dimension(3), intent(in) :: point
         real(kind=8) :: value
      end function manufactured_source_fun
   end interface

   type ManufacturedCase
      character(len=32) :: name = ''
      real(kind=8) :: initial_l2_norm = 0.d0
      procedure(manufactured_value_fun), pointer, nopass :: value => null()
      procedure(manufactured_source_fun), pointer, nopass :: source => null()
   end type ManufacturedCase

   type BenchmarkConfiguration
      character(len=16) :: scheme = ''
      real(kind=8) :: final_time = 0.d0
      real(kind=8) :: dt = 0.d0
      integer(kind=4) :: steps = 0
      integer(kind=4), dimension(3) :: nelem = 0
      integer(kind=4), dimension(3) :: p_test = 0
      integer(kind=4), dimension(3) :: p_trial = 0
      integer(kind=4), dimension(3) :: process_grid = 0
      integer(kind=4) :: sample_points = 0
      logical :: write_samples = .false.
   end type BenchmarkConfiguration

   type BenchmarkState
      type(ADS_Setup) :: ads_test
      type(ADS_Setup) :: ads_trial
      type(ADS_compute_data) :: ads_data
      type(TimeScheme3D) :: initial_scheme
      type(TimeScheme3D) :: evolution_scheme
      real(kind=8) :: step_wall_seconds = 0.d0
      logical :: parallel_active = .false.
      logical :: communicators_active = .false.
      logical :: ads_active = .false.
   end type BenchmarkState

   type BenchmarkMetrics
      real(kind=8) :: l2_error = 0.d0
      real(kind=8) :: linf_error = 0.d0
      real(kind=8) :: solution_l2_norm = 0.d0
      real(kind=8) :: field_checksum = 0.d0
   end type BenchmarkMetrics

   abstract interface
      subroutine adapter_initialize(config, state, status)
         import :: BenchmarkConfiguration, BenchmarkState
         type(BenchmarkConfiguration), intent(in) :: config
         type(BenchmarkState), intent(inout) :: state
         integer(kind=4), intent(out) :: status
      end subroutine adapter_initialize

      subroutine adapter_project_initial(config, state, status)
         import :: BenchmarkConfiguration, BenchmarkState
         type(BenchmarkConfiguration), intent(in) :: config
         type(BenchmarkState), intent(inout) :: state
         integer(kind=4), intent(out) :: status
      end subroutine adapter_project_initial

      subroutine adapter_advance(config, state, step, status)
         import :: BenchmarkConfiguration, BenchmarkState
         type(BenchmarkConfiguration), intent(in) :: config
         type(BenchmarkState), intent(inout) :: state
         integer(kind=4), intent(in) :: step
         integer(kind=4), intent(out) :: status
      end subroutine adapter_advance

      subroutine adapter_measure(config, state, physical_time, write_field, &
                                 metrics, status)
         import :: BenchmarkConfiguration, BenchmarkMetrics, BenchmarkState
         type(BenchmarkConfiguration), intent(in) :: config
         type(BenchmarkState), intent(in) :: state
         real(kind=8), intent(in) :: physical_time
         logical, intent(in) :: write_field
         type(BenchmarkMetrics), intent(out) :: metrics
         integer(kind=4), intent(out) :: status
      end subroutine adapter_measure

      subroutine adapter_cleanup(state, status)
         import :: BenchmarkState
         type(BenchmarkState), intent(inout) :: state
         integer(kind=4), intent(out) :: status
      end subroutine adapter_cleanup
   end interface

   type BenchmarkAdapter
      character(len=32) :: name = ''
      character(len=128) :: capabilities = ''
      type(ManufacturedCase) :: exact_case
      procedure(adapter_initialize), pointer, nopass :: initialize => null()
      procedure(adapter_project_initial), pointer, nopass :: project_initial => null()
      procedure(adapter_advance), pointer, nopass :: advance => null()
      procedure(adapter_measure), pointer, nopass :: measure => null()
      procedure(adapter_cleanup), pointer, nopass :: cleanup => null()
   end type BenchmarkAdapter

contains

   subroutine ValidateAdapter(adapter, status)
      type(BenchmarkAdapter), intent(in) :: adapter
      integer(kind=4), intent(out) :: status

      status = 0
      if (len_trim(adapter%name) == 0) status = 1
      if (len_trim(adapter%capabilities) == 0) status = 1
      if (len_trim(adapter%exact_case%name) == 0) status = 1
      if (.not. ieee_is_finite(adapter%exact_case%initial_l2_norm) .or. &
          adapter%exact_case%initial_l2_norm <= 0.d0) status = 1
      if (.not. associated(adapter%exact_case%value)) status = 1
      if (.not. associated(adapter%exact_case%source)) status = 1
      if (.not. associated(adapter%initialize)) status = 1
      if (.not. associated(adapter%project_initial)) status = 1
      if (.not. associated(adapter%advance)) status = 1
      if (.not. associated(adapter%measure)) status = 1
      if (.not. associated(adapter%cleanup)) status = 1

   end subroutine ValidateAdapter

end module benchmark_contract
