! Shared ADS lifecycle, stepping, sampling, and cleanup implementation.
module benchmark_runtime

   use benchmark_contract, ONLY: BenchmarkConfiguration, BenchmarkMetrics, &
                                 BenchmarkState

   implicit none

   private
   public :: AdvanceWithCallback
   public :: CleanupCommon
   public :: InitializeCommon
   public :: MeasureCommon
   public :: ProjectInitialCommon

contains

   subroutine InitializeCommon(config, state, status)
      use ADSS, ONLY: Initialize
      use communicators, ONLY: CreateCommunicators
      use parallelism, ONLY: InitializeParallelism
      use time_scheme, ONLY: ConfigureBackwardEuler3DTimeScheme, &
                             ConfigureDouglasGunn3DTimeScheme, &
                             ConfigureMassOnly3DTimeScheme, &
                             ConfigurePeacemanRachford3DTimeScheme, &
                             ValidateSpaces
      type(BenchmarkConfiguration), intent(in) :: config
      type(BenchmarkState), intent(inout) :: state
      integer(kind=4), intent(out) :: status

      status = 0
      call InitializeParallelism(config%process_grid(1), config%process_grid(2), &
                                 config%process_grid(3), status)
      if (status /= 0) return
      state%parallel_active = .true.

      call CreateCommunicators(status)
      state%communicators_active = .true.
      if (status /= 0) return

      state%ads_active = .true.
      call Initialize(config%nelem, config%p_test, config%p_trial, &
                      config%p_trial - 1, state%ads_test, state%ads_trial, &
                      state%ads_data, status)
      if (status /= 0) return
      state%ads_data%t = 0.d0
      call ValidateSpaces(state%ads_test, state%ads_trial)

      call ConfigureMassOnly3DTimeScheme(state%initial_scheme)
      select case (trim(config%scheme))
      case ('dg')
         call ConfigureDouglasGunn3DTimeScheme( &
            config%dt, state%evolution_scheme, include_transport=.false.)
      case ('pr')
         call ConfigurePeacemanRachford3DTimeScheme( &
            config%dt, state%evolution_scheme, include_transport=.false.)
      case ('be')
         call ConfigureBackwardEuler3DTimeScheme( &
            config%dt, state%evolution_scheme, include_transport=.false.)
      case default
         status = 5
      end select
   end subroutine InitializeCommon

   subroutine ProjectInitialCommon(config, state, status)
      use manufactured_solution, ONLY: InitialProjectionRhs, &
                                         ManufacturedScalarSource
      use time_scheme, ONLY: BackwardEuler3DStep, DouglasGunn3DStep, &
                             PeacemanRachford3DStep
      type(BenchmarkConfiguration), intent(in) :: config
      type(BenchmarkState), intent(inout) :: state
      integer(kind=4), intent(out) :: status

      state%ads_test%tau = 1.d0
      state%ads_trial%tau = 1.d0
      state%ads_data%t = 0.d0
      select case (trim(config%scheme))
      case ('dg')
         call DouglasGunn3DStep( &
            state%initial_scheme, 0, ManufacturedScalarSource, &
            state%ads_test, state%ads_trial, state%ads_data, 1, status, &
            InitialProjectionRhs)
      case ('pr')
         call PeacemanRachford3DStep( &
            state%initial_scheme, 0, ManufacturedScalarSource, &
            state%ads_test, state%ads_trial, state%ads_data, 1, status, &
            InitialProjectionRhs)
      case ('be')
         call BackwardEuler3DStep( &
            state%initial_scheme, 0, ManufacturedScalarSource, &
            state%ads_test, state%ads_trial, state%ads_data, 1, status, &
            InitialProjectionRhs)
      case default
         status = 5
      end select

   end subroutine ProjectInitialCommon

   subroutine AdvanceWithCallback(config, state, step, rhs_point, status)
      use Interfaces, ONLY: rhs_point_fun
      use manufactured_solution, ONLY: ManufacturedScalarSource, SetStepClock
      use time_scheme, ONLY: BackwardEuler3DStep, DouglasGunn3DStep, &
                             PeacemanRachford3DStep
      type(BenchmarkConfiguration), intent(in) :: config
      type(BenchmarkState), intent(inout) :: state
      integer(kind=4), intent(in) :: step
      procedure(rhs_point_fun) :: rhs_point
      integer(kind=4), intent(out) :: status
      real(kind=8) :: start_time

      start_time = real(step - 1, kind=8)*config%dt
      state%ads_data%t = start_time
      state%ads_test%tau = config%dt
      state%ads_trial%tau = config%dt
      call SetStepClock(config%scheme, start_time, config%dt)

      select case (trim(config%scheme))
      case ('dg')
         call DouglasGunn3DStep( &
            state%evolution_scheme, step, ManufacturedScalarSource, &
            state%ads_test, state%ads_trial, state%ads_data, 1, status, &
            rhs_point)
      case ('pr')
         call PeacemanRachford3DStep( &
            state%evolution_scheme, step, ManufacturedScalarSource, &
            state%ads_test, state%ads_trial, state%ads_data, 1, status, &
            rhs_point)
      case ('be')
         call BackwardEuler3DStep( &
            state%evolution_scheme, step, ManufacturedScalarSource, &
            state%ads_test, state%ads_trial, state%ads_data, 1, status, &
            rhs_point)
      case default
         status = 5
      end select
      if (status == 0) state%ads_data%t = real(step, kind=8)*config%dt

   end subroutine AdvanceWithCallback

   subroutine MeasureCommon(config, state, physical_time, write_field, &
                            metrics, status)
      type(BenchmarkConfiguration), intent(in) :: config
      type(BenchmarkState), intent(in) :: state
      real(kind=8), intent(in) :: physical_time
      logical, intent(in) :: write_field
      type(BenchmarkMetrics), intent(out) :: metrics
      integer(kind=4), intent(out) :: status

      call MeasureL2Metrics( &
         state, physical_time, metrics%l2_error, metrics%solution_l2_norm, &
         status)
      if (status /= 0) return
      call SampleField(config, state, physical_time, write_field, &
                       metrics%linf_error, metrics%field_checksum, status)

   end subroutine MeasureCommon

   subroutine MeasureL2Metrics( &
      state, physical_time, l2_error, solution_l2_norm, status)
      use basis, ONLY: EvalSpline
      use gauss, ONLY: GaussRule
      use manufactured_solution, ONLY: ExactSolution
      use my_mpi, ONLY: GatherFullSolution
      use parallelism, ONLY: MYRANK
      use mpi
      type(BenchmarkState), intent(in) :: state
      real(kind=8), intent(in) :: physical_time
      real(kind=8), intent(out) :: l2_error, solution_l2_norm
      integer(kind=4), intent(out) :: status
      real(kind=8), allocatable, dimension(:, :, :) :: coefficients
      real(kind=8), dimension(0:9) :: gauss_x, gauss_y, gauss_z
      real(kind=8), dimension(0:9) :: weight_x, weight_y, weight_z
      real(kind=8), dimension(2) :: values
      real(kind=8) :: error_compensation, error_sum, error_term
      real(kind=8) :: exact_value, jacobian, solution_compensation
      real(kind=8) :: solution_sum, solution_term, spline_value
      real(kind=8) :: update, x, x_lower, y, y_lower, z, z_lower
      integer(kind=4) :: ex, ey, ez, ix, iy, iz, mpi_ierr
      integer(kind=4) :: quadrature_x, quadrature_y, quadrature_z

      ! The solver's trial-space quadrature is deliberately not reused here.
      ! Reusing it aliases a non-polynomial exact field against the projection
      ! (for example p=1 can report a zero cosine error).  Three extra points
      ! per axis provide an independent rule; ten is the library table limit
      ! and also integrates the square of a degree-nine spline exactly.
      quadrature_x = min(10, state%ads_trial%p(1) + 3)
      quadrature_y = min(10, state%ads_trial%p(2) + 3)
      quadrature_z = min(10, state%ads_trial%p(3) + 3)
      call ValidateUnitCubeElementMap(state, status)
      if (status /= 0) return
      call GaussRule( &
         quadrature_x, gauss_x(0:quadrature_x - 1), &
         weight_x(0:quadrature_x - 1))
      call GaussRule( &
         quadrature_y, gauss_y(0:quadrature_y - 1), &
         weight_y(0:quadrature_y - 1))
      call GaussRule( &
         quadrature_z, gauss_z(0:quadrature_z - 1), &
         weight_z(0:quadrature_z - 1))
      call GatherFullSolution( &
         0, state%ads_data%FF, coefficients, state%ads_trial%n, &
         state%ads_trial%p, state%ads_trial%s)

      values = 0.d0
      if (MYRANK == 0) then
         error_sum = 0.d0
         error_compensation = 0.d0
         solution_sum = 0.d0
         solution_compensation = 0.d0
         do ez = 1, state%ads_trial%nelem(3)
            ! Knot vectors are allocated with lower bound one, while
            ! BasisData remaps them to its explicit lower bound zero.  Adding
            ! lbound(Uz,1) reproduces BasisData's U(O+p) element boundary.
            z_lower = state%ads_trial%Uz( &
               lbound(state%ads_trial%Uz, 1) + state%ads_trial%Oz(ez) + &
               state%ads_trial%p(3))
            do ey = 1, state%ads_trial%nelem(2)
               y_lower = state%ads_trial%Uy( &
                  lbound(state%ads_trial%Uy, 1) + state%ads_trial%Oy(ey) + &
                  state%ads_trial%p(2))
               do ex = 1, state%ads_trial%nelem(1)
                  x_lower = state%ads_trial%Ux( &
                     lbound(state%ads_trial%Ux, 1) + &
                     state%ads_trial%Ox(ex) + state%ads_trial%p(1))
                  jacobian = state%ads_trial%Jx(ex)* &
                             state%ads_trial%Jy(ey)* &
                             state%ads_trial%Jz(ez)
                  do iz = 0, quadrature_z - 1
                     z = z_lower + state%ads_trial%Jz(ez)*(gauss_z(iz) + 1.d0)
                     do iy = 0, quadrature_y - 1
                        y = y_lower + &
                            state%ads_trial%Jy(ey)*(gauss_y(iy) + 1.d0)
                        do ix = 0, quadrature_x - 1
                           x = x_lower + &
                               state%ads_trial%Jx(ex)*(gauss_x(ix) + 1.d0)
                           spline_value = EvalSpline(0, &
                              state%ads_trial%Ux, state%ads_trial%p(1), &
                              state%ads_trial%n(1), &
                              state%ads_trial%nelem(1), &
                              state%ads_trial%Uy, state%ads_trial%p(2), &
                              state%ads_trial%n(2), &
                              state%ads_trial%nelem(2), &
                              state%ads_trial%Uz, state%ads_trial%p(3), &
                              state%ads_trial%n(3), &
                              state%ads_trial%nelem(3), coefficients, x, y, z)
                           exact_value = ExactSolution( &
                              physical_time, (/x, y, z/))
                           error_term = &
                              (spline_value - exact_value)**2*jacobian* &
                              weight_x(ix)*weight_y(iy)*weight_z(iz) - &
                              error_compensation
                           update = error_sum + error_term
                           error_compensation = &
                              (update - error_sum) - error_term
                           error_sum = update
                           solution_term = spline_value*spline_value* &
                              jacobian*weight_x(ix)*weight_y(iy)* &
                              weight_z(iz) - solution_compensation
                           update = solution_sum + solution_term
                           solution_compensation = &
                              (update - solution_sum) - solution_term
                           solution_sum = update
                        end do
                     end do
                  end do
               end do
            end do
         end do
         values = (/sqrt(max(error_sum, 0.d0)), &
                    sqrt(max(solution_sum, 0.d0))/)
      end if

      call MPI_Bcast(values, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, &
                     mpi_ierr)
      if (mpi_ierr /= 0) status = mpi_ierr
      l2_error = values(1)
      solution_l2_norm = values(2)
      if (allocated(coefficients)) deallocate(coefficients)

   end subroutine MeasureL2Metrics

   subroutine ValidateUnitCubeElementMap(state, status)
      type(BenchmarkState), intent(in) :: state
      integer(kind=4), intent(out) :: status
      real(kind=8) :: tolerance

      tolerance = 64.d0*epsilon(1.d0)
      status = 0
      call ValidateAxisElementMap( &
         state%ads_trial%Ux, state%ads_trial%Ox, state%ads_trial%Jx, &
         state%ads_trial%p(1), state%ads_trial%nelem(1), tolerance, status)
      call ValidateAxisElementMap( &
         state%ads_trial%Uy, state%ads_trial%Oy, state%ads_trial%Jy, &
         state%ads_trial%p(2), state%ads_trial%nelem(2), tolerance, status)
      call ValidateAxisElementMap( &
         state%ads_trial%Uz, state%ads_trial%Oz, state%ads_trial%Jz, &
         state%ads_trial%p(3), state%ads_trial%nelem(3), tolerance, status)

   end subroutine ValidateUnitCubeElementMap

   subroutine ValidateAxisElementMap( &
      knots, offsets, jacobians, degree, elements, tolerance, status)
      real(kind=8), dimension(:), intent(in) :: knots, jacobians
      integer(kind=4), dimension(:), intent(in) :: offsets
      integer(kind=4), intent(in) :: degree, elements
      real(kind=8), intent(in) :: tolerance
      integer(kind=4), intent(inout) :: status
      real(kind=8) :: element_lower, element_upper, previous_upper
      integer(kind=4) :: element, knot_index

      if (status /= 0) return
      previous_upper = 0.d0
      do element = 1, elements
         knot_index = lbound(knots, 1) + offsets(element) + degree
         element_lower = knots(knot_index)
         element_upper = element_lower + 2.d0*jacobians(element)
         if (abs(element_lower - previous_upper) > tolerance .or. &
             element_upper <= element_lower .or. &
             element_lower < -tolerance .or. &
             element_upper > 1.d0 + tolerance) then
            status = 7
            return
         end if
         previous_upper = element_upper
      end do
      if (abs(previous_upper - 1.d0) > tolerance) status = 7

   end subroutine ValidateAxisElementMap

   subroutine SampleField(config, state, physical_time, write_field, &
                          linf_error, checksum, status)
      use basis, ONLY: EvalSpline
      use manufactured_solution, ONLY: ExactSolution
      use my_mpi, ONLY: GatherFullSolution
      use parallelism, ONLY: MYRANK
      use mpi
      type(BenchmarkConfiguration), intent(in) :: config
      type(BenchmarkState), intent(in) :: state
      real(kind=8), intent(in) :: physical_time
      logical, intent(in) :: write_field
      real(kind=8), intent(out) :: linf_error, checksum
      integer(kind=4), intent(out) :: status
      real(kind=8), allocatable, dimension(:, :, :) :: coefficients
      real(kind=8), dimension(2) :: values
      real(kind=8) :: compensation, exact_value, sample_value
      real(kind=8) :: weighted_compensation, weighted_sum, x, y, z
      real(kind=8) :: sum_value, update, corrected
      integer(kind=4) :: index, ix, iy, iz, io_status, mpi_ierr
      integer(kind=4), parameter :: sample_unit = 91

      call GatherFullSolution(0, state%ads_data%FF, coefficients, &
                              state%ads_trial%n, state%ads_trial%p, &
                              state%ads_trial%s)
      status = 0
      io_status = 0
      values = 0.d0
      if (MYRANK == 0) then
         if (write_field) then
            open(unit=sample_unit, file='field_samples.csv', status='replace', &
                 action='write', iostat=io_status)
            if (io_status /= 0) status = 6
            if (status == 0) write(sample_unit, '(A)', iostat=io_status) &
               'x,y,z,numerical,exact,error'
            if (io_status /= 0) status = 6
         end if

         linf_error = 0.d0
         sum_value = 0.d0
         compensation = 0.d0
         weighted_sum = 0.d0
         weighted_compensation = 0.d0
         index = 0
         do iz = 0, config%sample_points - 1
            z = real(iz, kind=8)/real(config%sample_points - 1, kind=8)
            do iy = 0, config%sample_points - 1
               y = real(iy, kind=8)/real(config%sample_points - 1, kind=8)
               do ix = 0, config%sample_points - 1
                  x = real(ix, kind=8)/real(config%sample_points - 1, kind=8)
                  index = index + 1
                  sample_value = EvalSpline(0, &
                     state%ads_trial%Ux, state%ads_trial%p(1), &
                     state%ads_trial%n(1), state%ads_trial%nelem(1), &
                     state%ads_trial%Uy, state%ads_trial%p(2), &
                     state%ads_trial%n(2), state%ads_trial%nelem(2), &
                     state%ads_trial%Uz, state%ads_trial%p(3), &
                     state%ads_trial%n(3), state%ads_trial%nelem(3), &
                     coefficients, x, y, z)
                  exact_value = ExactSolution(physical_time, (/x, y, z/))
                  linf_error = max(linf_error, abs(sample_value - exact_value))

                  corrected = sample_value - compensation
                  update = sum_value + corrected
                  compensation = (update - sum_value) - corrected
                  sum_value = update
                  corrected = real(index, kind=8)*sample_value - &
                              weighted_compensation
                  update = weighted_sum + corrected
                  weighted_compensation = (update - weighted_sum) - corrected
                  weighted_sum = update

                  if (write_field .and. status == 0) then
                     write(sample_unit, '(6(ES24.16E3,:,","))', &
                           iostat=io_status) x, y, z, sample_value, exact_value, &
                           sample_value - exact_value
                     if (io_status /= 0) status = 6
                  end if
               end do
            end do
         end do
         checksum = sum_value + weighted_sum/real(index + 1, kind=8)
         if (write_field) close(sample_unit, iostat=io_status)
         if (write_field .and. io_status /= 0) status = 6
         values = (/linf_error, checksum/)
      end if

      call MPI_Bcast(status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
      if (status == 0 .and. mpi_ierr /= 0) status = mpi_ierr
      call MPI_Bcast(values, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, &
                     mpi_ierr)
      if (status == 0 .and. mpi_ierr /= 0) status = mpi_ierr
      linf_error = values(1)
      checksum = values(2)
      if (allocated(coefficients)) deallocate(coefficients)

   end subroutine SampleField

   subroutine CleanupCommon(state, status)
      use ADSS, ONLY: Cleanup_ADS, Cleanup_data
      use communicators, ONLY: Cleanup_Communicators
      use parallelism, ONLY: Cleanup_Parallelism
      use mpi
      type(BenchmarkState), intent(inout) :: state
      integer(kind=4), intent(out) :: status
      integer(kind=4) :: any_failure, cleanup_status, local_failure, mpi_ierr

      status = 0
      if (state%ads_active) then
         call Cleanup_ADS(state%ads_test, cleanup_status)
         if (status == 0 .and. cleanup_status /= 0) status = cleanup_status
         call Cleanup_ADS(state%ads_trial, cleanup_status)
         if (status == 0 .and. cleanup_status /= 0) status = cleanup_status
         call Cleanup_data(state%ads_data, cleanup_status)
         if (status == 0 .and. cleanup_status /= 0) status = cleanup_status
         state%ads_active = .false.
      end if
      if (state%communicators_active) then
         call Cleanup_Communicators(cleanup_status)
         if (status == 0 .and. cleanup_status /= 0) status = cleanup_status
         state%communicators_active = .false.
      end if
      if (state%parallel_active) then
         local_failure = 0
         if (status /= 0) local_failure = 1
         call MPI_Allreduce(local_failure, any_failure, 1, MPI_INTEGER, &
                            MPI_MAX, MPI_COMM_WORLD, mpi_ierr)
         if (mpi_ierr /= 0 .and. status == 0) status = mpi_ierr
         if (mpi_ierr == 0 .and. any_failure /= 0 .and. status == 0) status = 1
         call Cleanup_Parallelism(cleanup_status)
         if (status == 0 .and. cleanup_status /= 0) status = cleanup_status
         state%parallel_active = .false.
      end if

   end subroutine CleanupCommon

end module benchmark_runtime
