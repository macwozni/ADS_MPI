!------------------------------------------------------------------------------
!> @file pollution_dpg_solver.F90
!> @brief Corrected, problem-local 3D directional DPG pollution solver.
!------------------------------------------------------------------------------
module pollution_dpg_solver

   implicit none
   private

   integer(kind=4), parameter, public :: POLLUTION_DPG_SUCCESS = 0
   integer(kind=4), parameter, public :: POLLUTION_DPG_INVALID_CONFIGURATION = -7601
   integer(kind=4), parameter, public :: POLLUTION_DPG_SYSTEM_TOO_LARGE = -7602
   integer(kind=4), parameter, public :: POLLUTION_DPG_LINEAR_SOLVE_ERROR = -7603
   integer(kind=4), parameter, public :: POLLUTION_DPG_RESIDUAL_TOO_LARGE = -7604
   integer(kind=4), parameter, public :: POLLUTION_DPG_IO_ERROR = -7605
   integer(kind=4), parameter, public :: POLLUTION_DPG_MPI_ERROR = -7606
   integer(kind=4), parameter, public :: POLLUTION_DPG_SOURCE_ERROR = -7607

   integer(kind=4), parameter :: MAX_POLYNOMIAL_DEGREE = 9
   integer(kind=4), parameter :: SOURCE_QUADRATURE = 8
   integer(kind=4), parameter :: SOURCE_SUBDIVISIONS = 8
   integer(kind=8), parameter :: MAX_TENSOR_ENTRIES = 50000000_8
   integer(kind=8), parameter :: MAX_WORK_ENTRIES = 50000000_8
   integer(kind=8), parameter :: MAX_OUTPUT_POINTS = 10000000_8
   real(kind=8), parameter :: RESIDUAL_TOLERANCE = 1.d-8

   abstract interface
      function WindFunction(time) result(wind)
         real(kind=8), intent(in) :: time
         real(kind=8) :: wind(3)
      end function WindFunction
   end interface

   interface
      subroutine dgetrf(m, n, a, lda, ipiv, info)
         integer, intent(in) :: m, n, lda
         integer, intent(out) :: ipiv(*), info
         double precision, intent(inout) :: a(lda, *)
      end subroutine dgetrf

      subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         character(len=1), intent(in) :: trans
         integer, intent(in) :: n, nrhs, lda, ldb, ipiv(*)
         integer, intent(out) :: info
         double precision, intent(in) :: a(lda, *)
         double precision, intent(inout) :: b(ldb, *)
      end subroutine dgetrs
   end interface

   type, public :: PollutionDPGStats
      integer(kind=4) :: completed_steps = 0
      real(kind=8) :: time = 0.d0
      real(kind=8) :: l2_norm = 0.d0
      real(kind=8) :: total_mass = 0.d0
      real(kind=8) :: maximum_abs = 0.d0
      real(kind=8) :: source_integral = 0.d0
      real(kind=8) :: maximum_relative_residual = 0.d0
      real(kind=8) :: wind(3) = 0.d0
   end type PollutionDPGStats

   type, public :: PollutionDPGSolution
      integer(kind=4) :: n(3) = -1
      integer(kind=4) :: p(3) = -1
      integer(kind=4) :: nelem(3) = 0
      real(kind=8) :: domain_length = 0.d0
      real(kind=8), allocatable :: Ux(:)
      real(kind=8), allocatable :: Uy(:)
      real(kind=8), allocatable :: Uz(:)
      real(kind=8), allocatable :: coefficients(:, :, :)
   end type PollutionDPGSolution

   type :: Space1D
      integer(kind=4) :: n = -1
      integer(kind=4) :: p = -1
      integer(kind=4) :: nelem = 0
      integer(kind=4) :: ng = 0
      real(kind=8), allocatable :: knots(:)
      integer(kind=4), allocatable :: offset(:)
      real(kind=8), allocatable :: jacobian(:)
      real(kind=8), allocatable :: weights(:)
      real(kind=8), allocatable :: points(:, :)
      real(kind=8), allocatable :: basis(:, :, :, :)
   end type Space1D

   type :: MatrixSet1D
      real(kind=8), allocatable :: mass(:, :)
      real(kind=8), allocatable :: stiffness(:, :)
      real(kind=8), allocatable :: advection(:, :)
   end type MatrixSet1D

   type :: DenseFactor
      integer(kind=4) :: n = 0
      real(kind=8), allocatable :: lu(:, :)
      integer(kind=4), allocatable :: pivots(:)
      logical :: ready = .false.
   end type DenseFactor

   type :: DirectionalOperator
      real(kind=8), allocatable :: B(:, :)
      real(kind=8), allocatable :: normal(:, :)
      type(DenseFactor) :: normal_factor
   end type DirectionalOperator

   type :: TensorData
      real(kind=8), allocatable :: values(:, :, :)
   end type TensorData

   public :: WindFunction
   public :: SolvePollutionDPG3D
   public :: EvaluatePollutionDPG
   public :: WritePollutionDPGVTI
   public :: CleanupPollutionDPGSolution

contains

!------------------------------------------------------------------------------
!> Run the corrected three-sweep DPG scheme on rank zero and replicate the
!> final trial coefficients.  The existing library API remains unchanged.
!------------------------------------------------------------------------------
subroutine SolvePollutionDPG3D(nelem, adapt_mesh, p_trial, c_trial, p_test, &
                               c_test, steps, dt, domain_length, diffusion, &
                               wind_function, source, source_min, source_max, &
                               output_resolution, write_output, solution, &
                               stats, ierr)
   use Interfaces, only: forcing_fun
   use parallelism, only: MYRANK
   use mpi

   implicit none

   integer(kind=4), intent(in) :: nelem, p_trial, c_trial, p_test, c_test
   integer(kind=4), intent(in) :: steps, output_resolution
   logical, intent(in) :: adapt_mesh, write_output
   real(kind=8), intent(in) :: dt, domain_length, diffusion(3)
   real(kind=8), intent(in) :: source_min(3), source_max(3)
   procedure(WindFunction) :: wind_function
   procedure(forcing_fun) :: source
   type(PollutionDPGSolution), intent(out) :: solution
   type(PollutionDPGStats), intent(out) :: stats
   integer(kind=4), intent(out) :: ierr

   integer(kind=4) :: mpi_ierr

   call CleanupPollutionDPGSolution(solution)
   stats = PollutionDPGStats()
   ierr = POLLUTION_DPG_SUCCESS

   if (MYRANK == 0) then
      call SolveOnRoot(nelem, adapt_mesh, p_trial, c_trial, p_test, c_test, &
                       steps, dt, domain_length, diffusion, wind_function, &
                       source, source_min, source_max, output_resolution, &
                       write_output, solution, stats, ierr)
   end if

   call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) then
      ierr = POLLUTION_DPG_MPI_ERROR
      return
   end if
   if (ierr /= POLLUTION_DPG_SUCCESS) return

   call BroadcastStats(stats, ierr)
   if (ierr /= POLLUTION_DPG_SUCCESS) return
   call BroadcastSolution(solution, ierr)

end subroutine SolvePollutionDPG3D

subroutine SolveOnRoot(nelem, adapt_mesh, p_trial, c_trial, p_test, c_test, &
                       steps, dt, domain_length, diffusion, wind_function, &
                       source, source_min, source_max, output_resolution, &
                       write_output, solution, stats, ierr)
   use Interfaces, only: forcing_fun
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none

   integer(kind=4), intent(in) :: nelem, p_trial, c_trial, p_test, c_test
   integer(kind=4), intent(in) :: steps, output_resolution
   logical, intent(in) :: adapt_mesh, write_output
   real(kind=8), intent(in) :: dt, domain_length, diffusion(3)
   real(kind=8), intent(in) :: source_min(3), source_max(3)
   procedure(WindFunction) :: wind_function
   procedure(forcing_fun) :: source
   type(PollutionDPGSolution), intent(out) :: solution
   type(PollutionDPGStats), intent(out) :: stats
   integer(kind=4), intent(out) :: ierr

   type(Space1D) :: trial(3), test(3)
   type(MatrixSet1D) :: trial_matrix(3), test_matrix(3), mixed_matrix(3)
   type(DenseFactor) :: mass_factor(3), gram_factor(3)
   type(DirectionalOperator) :: directional(3)
   type(TensorData) :: source_tensor(3)
   real(kind=8), allocatable :: gram(:, :), rhs(:, :, :)
   real(kind=8) :: wind(3), h, residual
   integer(kind=4) :: axis, step
   character(len=64) :: filename

   stats = PollutionDPGStats()
   call ValidateConfiguration(nelem, p_trial, c_trial, p_test, c_test, &
                              steps, dt, domain_length, diffusion, &
                              output_resolution, ierr)
   if (ierr /= POLLUTION_DPG_SUCCESS) return
   if (.not. all(ieee_is_finite(source_min)) .or. &
       .not. all(ieee_is_finite(source_max)) .or. &
       any(source_max <= source_min)) then
      ierr = POLLUTION_DPG_INVALID_CONFIGURATION
      return
   end if

   do axis = 1, 3
      call BuildSpace1D(nelem, p_trial, c_trial, &
                        max(p_trial, p_test) + 1, domain_length, &
                        adapt_mesh .and. axis == 1, trial(axis), ierr)
      if (ierr /= POLLUTION_DPG_SUCCESS) return
      call BuildSpace1D(nelem, p_test, c_test, &
                        max(p_trial, p_test) + 1, domain_length, &
                        adapt_mesh .and. axis == 1, test(axis), ierr)
      if (ierr /= POLLUTION_DPG_SUCCESS) return

      call AssembleMatrixSet(trial(axis), trial(axis), trial_matrix(axis))
      call AssembleMatrixSet(test(axis), test(axis), test_matrix(axis))
      call AssembleMatrixSet(test(axis), trial(axis), mixed_matrix(axis))

      call FactorDense(trial_matrix(axis)%mass, mass_factor(axis), ierr)
      if (ierr /= POLLUTION_DPG_SUCCESS) return
      allocate (gram(size(test_matrix(axis)%mass, 1), &
                     size(test_matrix(axis)%mass, 2)))
      gram = test_matrix(axis)%mass + test_matrix(axis)%stiffness
      call FactorDense(gram, gram_factor(axis), ierr)
      deallocate (gram)
      if (ierr /= POLLUTION_DPG_SUCCESS) return
   end do

   call InitializeSolution(trial, domain_length, solution)
   call AssembleSourceTensor(test(1), trial(2), trial(3), source, &
                             source_min, source_max, source_tensor(1)%values, &
                             ierr)
   if (ierr /= POLLUTION_DPG_SUCCESS) return
   call AssembleSourceTensor(trial(1), test(2), trial(3), source, &
                             source_min, source_max, source_tensor(2)%values, &
                             ierr)
   if (ierr /= POLLUTION_DPG_SUCCESS) return
   call AssembleSourceTensor(trial(1), trial(2), test(3), source, &
                             source_min, source_max, source_tensor(3)%values, &
                             ierr)
   if (ierr /= POLLUTION_DPG_SUCCESS) return

   stats%source_integral = sum(source_tensor(1)%values)
   if (.not. ieee_is_finite(stats%source_integral) .or. &
       stats%source_integral <= 0.d0) then
      ierr = POLLUTION_DPG_SOURCE_ERROR
      return
   end if

   if (write_output) then
      call WritePollutionDPGVTI(solution, 'out_0.vti', output_resolution, ierr)
      if (ierr /= POLLUTION_DPG_SUCCESS) return
   end if

   h = dt/3.d0
   do step = 1, steps
      stats%time = real(step - 1, kind=8)*dt
      wind = wind_function(stats%time)
      if (.not. all(ieee_is_finite(wind))) then
         ierr = POLLUTION_DPG_INVALID_CONFIGURATION
         return
      end if
      stats%wind = wind

      do axis = 1, 3
         call BuildDirectionalOperator(mixed_matrix(axis), gram_factor(axis), &
                                       h, diffusion(axis), wind(axis), &
                                       directional(axis), ierr)
         if (ierr /= POLLUTION_DPG_SUCCESS) return
      end do

      do axis = 1, 3
         call BuildStageRHS(axis, solution%coefficients, h, diffusion, wind, &
                            trial_matrix, mixed_matrix, &
                            source_tensor(axis)%values, rhs)
         call AdvanceDirection(axis, rhs, directional(axis), &
                               gram_factor(axis), mass_factor, &
                               trial_matrix, solution%coefficients, &
                               residual, ierr)
         if (allocated(rhs)) deallocate (rhs)
         if (ierr /= POLLUTION_DPG_SUCCESS) return
         if (.not. ieee_is_finite(residual)) then
            ierr = POLLUTION_DPG_RESIDUAL_TOO_LARGE
            return
         end if
         stats%maximum_relative_residual = &
            max(stats%maximum_relative_residual, residual)
      end do

      if (.not. all(ieee_is_finite(solution%coefficients))) then
         ierr = POLLUTION_DPG_LINEAR_SOLVE_ERROR
         return
      end if
      stats%completed_steps = step
      stats%time = real(step, kind=8)*dt
      call ComputeDiagnostics(solution%coefficients, trial_matrix, &
                              stats%l2_norm, stats%total_mass, &
                              stats%maximum_abs)
      if (.not. ieee_is_finite(stats%l2_norm) .or. &
          .not. ieee_is_finite(stats%total_mass) .or. &
          .not. ieee_is_finite(stats%maximum_abs) .or. &
          .not. ieee_is_finite(stats%maximum_relative_residual) .or. &
          stats%maximum_relative_residual > RESIDUAL_TOLERANCE) then
         ierr = POLLUTION_DPG_RESIDUAL_TOO_LARGE
         return
      end if

      if (write_output) then
         write (filename, '("out_",I0,".vti")') step
         call WritePollutionDPGVTI(solution, trim(filename), &
                                   output_resolution, ierr)
         if (ierr /= POLLUTION_DPG_SUCCESS) return
      end if
   end do

end subroutine SolveOnRoot

subroutine ValidateConfiguration(nelem, p_trial, c_trial, p_test, c_test, &
                                 steps, dt, domain_length, diffusion, &
                                 output_resolution, ierr)
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none

   integer(kind=4), intent(in) :: nelem, p_trial, c_trial, p_test, c_test
   integer(kind=4), intent(in) :: steps, output_resolution
   real(kind=8), intent(in) :: dt, domain_length, diffusion(3)
   integer(kind=4), intent(out) :: ierr
   integer(kind=8) :: trial_dofs, test_dofs, trial_entries, mixed_entries
   integer(kind=8) :: output_points, estimated_work

   ierr = POLLUTION_DPG_INVALID_CONFIGURATION
   if (nelem <= 0 .or. p_trial < 1 .or. p_test < 1) return
   if (p_trial > MAX_POLYNOMIAL_DEGREE .or. &
       p_test > MAX_POLYNOMIAL_DEGREE) return
   ! This volume-only formulation is conforming.  Broken C^-1 spaces need
   ! facet fluxes and trace terms which the upstream prototype does not have.
   if (c_trial < 0 .or. c_trial > p_trial - 1) return
   if (c_test < 0 .or. c_test > p_test - 1) return
   if (steps <= 0 .or. output_resolution <= 0) return
   if (.not. ieee_is_finite(dt) .or. .not. ieee_is_finite(domain_length)) return
   if (.not. all(ieee_is_finite(diffusion))) return
   if (dt <= 0.d0 .or. domain_length <= 0.d0) return
   if (any(diffusion < 0.d0)) return

   trial_dofs = int(p_trial + 1, kind=8) + &
      int(nelem - 1, kind=8)*int(p_trial - c_trial, kind=8)
   test_dofs = int(p_test + 1, kind=8) + &
      int(nelem - 1, kind=8)*int(p_test - c_test, kind=8)
   if (test_dofs < trial_dofs .or. trial_dofs <= 0_8) return
   if (trial_dofs > huge(0_4) .or. test_dofs > huge(0_4)) then
      ierr = POLLUTION_DPG_SYSTEM_TOO_LARGE
      return
   end if
   if (test_dofs > MAX_TENSOR_ENTRIES/test_dofs) then
      ierr = POLLUTION_DPG_SYSTEM_TOO_LARGE
      return
   end if
   if (trial_dofs > MAX_TENSOR_ENTRIES/trial_dofs) then
      ierr = POLLUTION_DPG_SYSTEM_TOO_LARGE
      return
   end if
   trial_entries = trial_dofs*trial_dofs
   if (trial_entries > MAX_TENSOR_ENTRIES/trial_dofs) then
      ierr = POLLUTION_DPG_SYSTEM_TOO_LARGE
      return
   end if
   trial_entries = trial_entries*trial_dofs
   if (trial_entries > MAX_TENSOR_ENTRIES) then
      ierr = POLLUTION_DPG_SYSTEM_TOO_LARGE
      return
   end if
   if (test_dofs > MAX_TENSOR_ENTRIES/trial_dofs) then
      ierr = POLLUTION_DPG_SYSTEM_TOO_LARGE
      return
   end if
   mixed_entries = test_dofs*trial_dofs
   if (mixed_entries > MAX_TENSOR_ENTRIES/trial_dofs) then
      ierr = POLLUTION_DPG_SYSTEM_TOO_LARGE
      return
   end if
   mixed_entries = mixed_entries*trial_dofs
   estimated_work = 5_8*mixed_entries + 4_8*trial_entries
   if (estimated_work > MAX_WORK_ENTRIES) then
      ierr = POLLUTION_DPG_SYSTEM_TOO_LARGE
      return
   end if
   output_points = int(output_resolution, kind=8) + 1_8
   if (output_points > MAX_OUTPUT_POINTS/output_points) then
      ierr = POLLUTION_DPG_SYSTEM_TOO_LARGE
      return
   end if
   trial_entries = output_points*output_points
   if (trial_entries > MAX_OUTPUT_POINTS/output_points) then
      ierr = POLLUTION_DPG_SYSTEM_TOO_LARGE
      return
   end if
   ierr = POLLUTION_DPG_SUCCESS

end subroutine ValidateConfiguration

subroutine BuildSpace1D(nelem, degree, continuity, ng, domain_length, &
                        adapted, space, ierr)
   use knot_vector, only: PrepareKnot
   use basis, only: BasisData

   implicit none

   integer(kind=4), intent(in) :: nelem, degree, continuity, ng
   real(kind=8), intent(in) :: domain_length
   logical, intent(in) :: adapted
   type(Space1D), intent(out) :: space
   integer(kind=4), intent(out) :: ierr
   integer(kind=4) :: base_n, generated_nelem, index_value, m
   real(kind=8) :: unit_coordinate

   ierr = POLLUTION_DPG_SUCCESS
   base_n = nelem + degree - 1
   space%n = base_n
   call PrepareKnot(space%n, degree, 1, continuity, space%knots, &
                    generated_nelem)
   if (generated_nelem /= nelem) then
      ierr = POLLUTION_DPG_INVALID_CONFIGURATION
      return
   end if
   space%p = degree
   space%nelem = nelem
   space%ng = ng

   do index_value = 1, size(space%knots)
      unit_coordinate = space%knots(index_value)
      if (adapted) then
         space%knots(index_value) = domain_length*AdaptedUnitCoordinate(unit_coordinate)
      else
         space%knots(index_value) = domain_length*unit_coordinate
      end if
   end do

   allocate (space%offset(nelem), space%jacobian(nelem))
   allocate (space%weights(ng), space%points(ng, nelem))
   allocate (space%basis(0:1, 0:degree, ng, nelem))
   m = space%n + degree + 1
   call BasisData(degree, m, space%knots, 1, ng, nelem, space%offset, &
                  space%jacobian, space%weights, space%points, space%basis)

end subroutine BuildSpace1D

pure real(kind=8) function AdaptedUnitCoordinate(t) result(mapped)
   implicit none
   real(kind=8), intent(in) :: t

   if (t < 0.5d0) then
      mapped = 1.98d0*t
   else
      mapped = 0.99d0 + 0.02d0*(t - 0.5d0)
   end if

end function AdaptedUnitCoordinate

subroutine AssembleMatrixSet(row_space, column_space, matrices)
   implicit none

   type(Space1D), intent(in) :: row_space, column_space
   type(MatrixSet1D), intent(out) :: matrices
   integer(kind=4) :: element, q, a, b, row, column
   real(kind=8) :: scale, row_value, row_derivative
   real(kind=8) :: column_value, column_derivative

   allocate (matrices%mass(row_space%n + 1, column_space%n + 1))
   allocate (matrices%stiffness(row_space%n + 1, column_space%n + 1))
   allocate (matrices%advection(row_space%n + 1, column_space%n + 1))
   matrices%mass = 0.d0
   matrices%stiffness = 0.d0
   matrices%advection = 0.d0

   do element = 1, row_space%nelem
      do q = 1, row_space%ng
         scale = row_space%jacobian(element)*row_space%weights(q)
         do a = 0, row_space%p
            row = row_space%offset(element) + a + 1
            row_value = row_space%basis(0, a, q, element)
            row_derivative = row_space%basis(1, a, q, element)
            do b = 0, column_space%p
               column = column_space%offset(element) + b + 1
               column_value = column_space%basis(0, b, q, element)
               column_derivative = column_space%basis(1, b, q, element)
               matrices%mass(row, column) = matrices%mass(row, column) + &
                  row_value*column_value*scale
               matrices%stiffness(row, column) = &
                  matrices%stiffness(row, column) + &
                  row_derivative*column_derivative*scale
               matrices%advection(row, column) = &
                  matrices%advection(row, column) + &
                  row_value*column_derivative*scale
            end do
         end do
      end do
   end do

end subroutine AssembleMatrixSet

subroutine BuildDirectionalOperator(mixed, gram_factor, h, diffusion, wind, &
                                    directional, ierr)
   implicit none

   type(MatrixSet1D), intent(in) :: mixed
   type(DenseFactor), intent(in) :: gram_factor
   real(kind=8), intent(in) :: h, diffusion, wind
   type(DirectionalOperator), intent(inout) :: directional
   integer(kind=4), intent(out) :: ierr
   real(kind=8), allocatable :: transformed(:, :)

   if (allocated(directional%B)) deallocate (directional%B)
   if (allocated(directional%normal)) deallocate (directional%normal)
   allocate (directional%B(size(mixed%mass, 1), size(mixed%mass, 2)))
   directional%B = mixed%mass + h*(diffusion*mixed%stiffness + &
                                    wind*mixed%advection)
   allocate (transformed(size(directional%B, 1), size(directional%B, 2)))
   transformed = directional%B
   call SolveFactorMatrix(gram_factor, transformed, size(transformed, 2), ierr)
   if (ierr /= POLLUTION_DPG_SUCCESS) then
      deallocate (transformed)
      return
   end if
   allocate (directional%normal(size(directional%B, 2), &
                                size(directional%B, 2)))
   directional%normal = matmul(transpose(directional%B), transformed)
   directional%normal = 0.5d0*(directional%normal + &
                                transpose(directional%normal))
   deallocate (transformed)
   call FactorDense(directional%normal, directional%normal_factor, ierr)

end subroutine BuildDirectionalOperator

subroutine BuildStageRHS(axis, coefficients, h, diffusion, wind, trial, &
                         mixed, source_values, rhs)
   implicit none

   integer(kind=4), intent(in) :: axis
   real(kind=8), intent(in) :: coefficients(:, :, :), h, diffusion(3), wind(3)
   type(MatrixSet1D), intent(in) :: trial(3), mixed(3)
   real(kind=8), intent(in) :: source_values(:, :, :)
   real(kind=8), allocatable, intent(out) :: rhs(:, :, :)
   real(kind=8), allocatable :: term(:, :, :)

   select case (axis)
   case (1)
      call ApplyTensor(mixed(1)%mass, trial(2)%mass, trial(3)%mass, &
                       coefficients, rhs)
      call ApplyTensor(mixed(1)%mass, trial(2)%stiffness, trial(3)%mass, &
                       coefficients, term)
      rhs = rhs - h*diffusion(2)*term
      deallocate (term)
      call ApplyTensor(mixed(1)%mass, trial(2)%advection, trial(3)%mass, &
                       coefficients, term)
      rhs = rhs - h*wind(2)*term
      deallocate (term)
      call ApplyTensor(mixed(1)%mass, trial(2)%mass, trial(3)%stiffness, &
                       coefficients, term)
      rhs = rhs - h*diffusion(3)*term
      deallocate (term)
      call ApplyTensor(mixed(1)%mass, trial(2)%mass, trial(3)%advection, &
                       coefficients, term)
      rhs = rhs - h*wind(3)*term
   case (2)
      call ApplyTensor(trial(1)%mass, mixed(2)%mass, trial(3)%mass, &
                       coefficients, rhs)
      call ApplyTensor(trial(1)%stiffness, mixed(2)%mass, trial(3)%mass, &
                       coefficients, term)
      rhs = rhs - h*diffusion(1)*term
      deallocate (term)
      call ApplyTensor(trial(1)%advection, mixed(2)%mass, trial(3)%mass, &
                       coefficients, term)
      rhs = rhs - h*wind(1)*term
      deallocate (term)
      call ApplyTensor(trial(1)%mass, mixed(2)%mass, trial(3)%stiffness, &
                       coefficients, term)
      rhs = rhs - h*diffusion(3)*term
      deallocate (term)
      call ApplyTensor(trial(1)%mass, mixed(2)%mass, trial(3)%advection, &
                       coefficients, term)
      rhs = rhs - h*wind(3)*term
   case (3)
      call ApplyTensor(trial(1)%mass, trial(2)%mass, mixed(3)%mass, &
                       coefficients, rhs)
      call ApplyTensor(trial(1)%stiffness, trial(2)%mass, mixed(3)%mass, &
                       coefficients, term)
      rhs = rhs - h*diffusion(1)*term
      deallocate (term)
      call ApplyTensor(trial(1)%advection, trial(2)%mass, mixed(3)%mass, &
                       coefficients, term)
      rhs = rhs - h*wind(1)*term
      deallocate (term)
      call ApplyTensor(trial(1)%mass, trial(2)%stiffness, mixed(3)%mass, &
                       coefficients, term)
      rhs = rhs - h*diffusion(2)*term
      deallocate (term)
      call ApplyTensor(trial(1)%mass, trial(2)%advection, mixed(3)%mass, &
                       coefficients, term)
      rhs = rhs - h*wind(2)*term
   end select
   rhs = rhs + h*source_values
   if (allocated(term)) deallocate (term)

end subroutine BuildStageRHS

subroutine AdvanceDirection(axis, rhs, directional, gram_factor, &
                            mass_factor, trial, coefficients, &
                            relative_residual, ierr)
   implicit none

   integer(kind=4), intent(in) :: axis
   real(kind=8), intent(in) :: rhs(:, :, :)
   type(DirectionalOperator), intent(in) :: directional
   type(DenseFactor), intent(in) :: gram_factor, mass_factor(3)
   type(MatrixSet1D), intent(in) :: trial(3)
   real(kind=8), allocatable, intent(inout) :: coefficients(:, :, :)
   real(kind=8), intent(out) :: relative_residual
   integer(kind=4), intent(out) :: ierr
   real(kind=8), allocatable :: representer(:, :, :), normal_rhs(:, :, :)
   real(kind=8), allocatable :: product_values(:, :, :)
   real(kind=8) :: rhs_norm, residual_norm

   allocate (representer(size(rhs, 1), size(rhs, 2), size(rhs, 3)))
   representer = rhs
   call SolveAlongAxis(gram_factor, representer, axis, ierr)
   if (ierr /= POLLUTION_DPG_SUCCESS) return
   call ApplyAlongAxis(transpose(directional%B), representer, axis, normal_rhs)
   deallocate (representer)

   if (allocated(coefficients)) deallocate (coefficients)
   allocate (coefficients(size(normal_rhs, 1), size(normal_rhs, 2), &
                          size(normal_rhs, 3)))
   coefficients = normal_rhs
   call SolveAlongAxis(directional%normal_factor, coefficients, axis, ierr)
   if (ierr /= POLLUTION_DPG_SUCCESS) return
   select case (axis)
   case (1)
      call SolveAlongAxis(mass_factor(2), coefficients, 2, ierr)
      if (ierr == POLLUTION_DPG_SUCCESS) &
         call SolveAlongAxis(mass_factor(3), coefficients, 3, ierr)
      if (ierr == POLLUTION_DPG_SUCCESS) &
         call ApplyTensor(directional%normal, trial(2)%mass, &
                          trial(3)%mass, coefficients, product_values)
   case (2)
      call SolveAlongAxis(mass_factor(1), coefficients, 1, ierr)
      if (ierr == POLLUTION_DPG_SUCCESS) &
         call SolveAlongAxis(mass_factor(3), coefficients, 3, ierr)
      if (ierr == POLLUTION_DPG_SUCCESS) &
         call ApplyTensor(trial(1)%mass, directional%normal, &
                          trial(3)%mass, coefficients, product_values)
   case (3)
      call SolveAlongAxis(mass_factor(1), coefficients, 1, ierr)
      if (ierr == POLLUTION_DPG_SUCCESS) &
         call SolveAlongAxis(mass_factor(2), coefficients, 2, ierr)
      if (ierr == POLLUTION_DPG_SUCCESS) &
         call ApplyTensor(trial(1)%mass, trial(2)%mass, &
                          directional%normal, coefficients, product_values)
   end select
   if (ierr /= POLLUTION_DPG_SUCCESS) return

   product_values = product_values - normal_rhs
   residual_norm = StableTensorNorm(product_values)
   rhs_norm = StableTensorNorm(normal_rhs)
   if (rhs_norm > 0.d0) then
      relative_residual = residual_norm/rhs_norm
   else if (residual_norm == 0.d0) then
      relative_residual = 0.d0
   else
      relative_residual = huge(0.d0)
   end if
   deallocate (normal_rhs, product_values)

end subroutine AdvanceDirection

subroutine FactorDense(matrix, factor, ierr)
   implicit none

   real(kind=8), intent(in) :: matrix(:, :)
   type(DenseFactor), intent(inout) :: factor
   integer(kind=4), intent(out) :: ierr
   integer(kind=4) :: lapack_info

   ierr = POLLUTION_DPG_LINEAR_SOLVE_ERROR
   factor%ready = .false.
   if (size(matrix, 1) /= size(matrix, 2) .or. size(matrix, 1) <= 0) return
   factor%n = size(matrix, 1)
   if (allocated(factor%lu)) deallocate (factor%lu)
   if (allocated(factor%pivots)) deallocate (factor%pivots)
   allocate (factor%lu(factor%n, factor%n), factor%pivots(factor%n))
   factor%lu = matrix
   call dgetrf(factor%n, factor%n, factor%lu, factor%n, factor%pivots, &
               lapack_info)
   if (lapack_info /= 0) return
   factor%ready = .true.
   ierr = POLLUTION_DPG_SUCCESS

end subroutine FactorDense

subroutine SolveFactorMatrix(factor, rhs, column_count, ierr)
   implicit none

   type(DenseFactor), intent(in) :: factor
   integer(kind=4), intent(in) :: column_count
   real(kind=8), intent(inout) :: rhs(factor%n, column_count)
   integer(kind=4), intent(out) :: ierr
   integer(kind=4) :: lapack_info

   ierr = POLLUTION_DPG_LINEAR_SOLVE_ERROR
   if (.not. factor%ready .or. column_count <= 0) return
   call dgetrs('N', factor%n, column_count, factor%lu, factor%n, &
               factor%pivots, rhs, factor%n, lapack_info)
   if (lapack_info /= 0) return
   ierr = POLLUTION_DPG_SUCCESS

end subroutine SolveFactorMatrix

subroutine SolveAlongAxis(factor, values, axis, ierr)
   implicit none

   type(DenseFactor), intent(in) :: factor
   real(kind=8), intent(inout) :: values(:, :, :)
   integer(kind=4), intent(in) :: axis
   integer(kind=4), intent(out) :: ierr
   real(kind=8), allocatable :: rhs(:, :)
   integer(kind=4) :: ix, iy, iz, column

   select case (axis)
   case (1)
      if (size(values, 1) /= factor%n) then
         ierr = POLLUTION_DPG_LINEAR_SOLVE_ERROR
         return
      end if
      allocate (rhs(factor%n, size(values, 2)*size(values, 3)))
      rhs = reshape(values, shape(rhs))
      call SolveFactorMatrix(factor, rhs, size(rhs, 2), ierr)
      if (ierr == POLLUTION_DPG_SUCCESS) values = reshape(rhs, shape(values))
   case (2)
      if (size(values, 2) /= factor%n) then
         ierr = POLLUTION_DPG_LINEAR_SOLVE_ERROR
         return
      end if
      allocate (rhs(factor%n, size(values, 1)*size(values, 3)))
      column = 0
      do iz = 1, size(values, 3)
         do ix = 1, size(values, 1)
            column = column + 1
            rhs(:, column) = values(ix, :, iz)
         end do
      end do
      call SolveFactorMatrix(factor, rhs, size(rhs, 2), ierr)
      if (ierr == POLLUTION_DPG_SUCCESS) then
         column = 0
         do iz = 1, size(values, 3)
            do ix = 1, size(values, 1)
               column = column + 1
               values(ix, :, iz) = rhs(:, column)
            end do
         end do
      end if
   case (3)
      if (size(values, 3) /= factor%n) then
         ierr = POLLUTION_DPG_LINEAR_SOLVE_ERROR
         return
      end if
      allocate (rhs(factor%n, size(values, 1)*size(values, 2)))
      column = 0
      do iy = 1, size(values, 2)
         do ix = 1, size(values, 1)
            column = column + 1
            rhs(:, column) = values(ix, iy, :)
         end do
      end do
      call SolveFactorMatrix(factor, rhs, size(rhs, 2), ierr)
      if (ierr == POLLUTION_DPG_SUCCESS) then
         column = 0
         do iy = 1, size(values, 2)
            do ix = 1, size(values, 1)
               column = column + 1
               values(ix, iy, :) = rhs(:, column)
            end do
         end do
      end if
   case default
      ierr = POLLUTION_DPG_LINEAR_SOLVE_ERROR
      return
   end select
   if (allocated(rhs)) deallocate (rhs)

end subroutine SolveAlongAxis

pure function StableTensorNorm(values) result(norm)
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none

   real(kind=8), intent(in) :: values(:, :, :)
   real(kind=8) :: norm, scale, sum_squares, magnitude
   integer(kind=4) :: ix, iy, iz

   scale = 0.d0
   sum_squares = 1.d0
   do iz = 1, size(values, 3)
      do iy = 1, size(values, 2)
         do ix = 1, size(values, 1)
            magnitude = abs(values(ix, iy, iz))
            if (.not. ieee_is_finite(magnitude)) then
               norm = magnitude
               return
            end if
            if (magnitude == 0.d0) cycle
            if (scale < magnitude) then
               sum_squares = 1.d0 + sum_squares*(scale/magnitude)**2
               scale = magnitude
            else
               sum_squares = sum_squares + (magnitude/scale)**2
            end if
         end do
      end do
   end do
   if (scale == 0.d0) then
      norm = 0.d0
   else
      norm = scale*sqrt(sum_squares)
   end if

end function StableTensorNorm

subroutine ApplyAlongAxis(matrix, input_values, axis, output_values)
   implicit none

   real(kind=8), intent(in) :: matrix(:, :), input_values(:, :, :)
   integer(kind=4), intent(in) :: axis
   real(kind=8), allocatable, intent(out) :: output_values(:, :, :)
   integer(kind=4) :: ix, iy, iz

   select case (axis)
   case (1)
      allocate (output_values(size(matrix, 1), size(input_values, 2), &
                              size(input_values, 3)))
      do iz = 1, size(input_values, 3)
         do iy = 1, size(input_values, 2)
            output_values(:, iy, iz) = matmul(matrix, input_values(:, iy, iz))
         end do
      end do
   case (2)
      allocate (output_values(size(input_values, 1), size(matrix, 1), &
                              size(input_values, 3)))
      do iz = 1, size(input_values, 3)
         do ix = 1, size(input_values, 1)
            output_values(ix, :, iz) = matmul(matrix, input_values(ix, :, iz))
         end do
      end do
   case (3)
      allocate (output_values(size(input_values, 1), size(input_values, 2), &
                              size(matrix, 1)))
      do iy = 1, size(input_values, 2)
         do ix = 1, size(input_values, 1)
            output_values(ix, iy, :) = matmul(matrix, input_values(ix, iy, :))
         end do
      end do
   end select

end subroutine ApplyAlongAxis

subroutine ApplyTensor(matrix_x, matrix_y, matrix_z, input_values, &
                       output_values)
   implicit none

   real(kind=8), intent(in) :: matrix_x(:, :), matrix_y(:, :), matrix_z(:, :)
   real(kind=8), intent(in) :: input_values(:, :, :)
   real(kind=8), allocatable, intent(out) :: output_values(:, :, :)
   real(kind=8), allocatable :: temporary_x(:, :, :), temporary_y(:, :, :)

   call ApplyAlongAxis(matrix_x, input_values, 1, temporary_x)
   call ApplyAlongAxis(matrix_y, temporary_x, 2, temporary_y)
   deallocate (temporary_x)
   call ApplyAlongAxis(matrix_z, temporary_y, 3, output_values)
   deallocate (temporary_y)

end subroutine ApplyTensor

subroutine AssembleSourceTensor(space_x, space_y, space_z, source, &
                                source_min, source_max, values, ierr)
   use Interfaces, only: forcing_fun
   use gauss, only: GaussRule
   use basis, only: FindSpan, DersBasisFuns

   implicit none

   type(Space1D), intent(in) :: space_x, space_y, space_z
   procedure(forcing_fun) :: source
   real(kind=8), intent(in) :: source_min(3), source_max(3)
   real(kind=8), allocatable, intent(out) :: values(:, :, :)
   integer(kind=4), intent(out) :: ierr
   real(kind=8), allocatable :: breaks_x(:), breaks_y(:), breaks_z(:)
   real(kind=8), allocatable :: bx(:, :), by(:, :), bz(:, :)
   real(kind=8) :: gauss_x(0:SOURCE_QUADRATURE - 1)
   real(kind=8) :: gauss_w(0:SOURCE_QUADRATURE - 1)
   real(kind=8) :: X(3), gradient(3), value, weight
   real(kind=8) :: midpoint_x, midpoint_y, midpoint_z
   real(kind=8) :: half_x, half_y, half_z
   integer(kind=4) :: cell_x, cell_y, cell_z, qx, qy, qz
   integer(kind=4) :: span_x, span_y, span_z, first_x, first_y, first_z
   integer(kind=4) :: ax, ay, az, ix, iy, iz

   ierr = POLLUTION_DPG_SOURCE_ERROR
   allocate (values(space_x%n + 1, space_y%n + 1, space_z%n + 1))
   values = 0.d0
   if (any(source_max <= source_min)) return
   call BuildIntegrationBreaks(space_x, source_min(1), source_max(1), breaks_x)
   call BuildIntegrationBreaks(space_y, source_min(2), source_max(2), breaks_y)
   call BuildIntegrationBreaks(space_z, source_min(3), source_max(3), breaks_z)
   if (size(breaks_x) < 2 .or. size(breaks_y) < 2 .or. &
       size(breaks_z) < 2) return

   call GaussRule(SOURCE_QUADRATURE, gauss_x, gauss_w)
   allocate (bx(0:space_x%p, 0:0), by(0:space_y%p, 0:0), &
             bz(0:space_z%p, 0:0))
   gradient = 0.d0
   do cell_z = 1, size(breaks_z) - 1
      midpoint_z = 0.5d0*(breaks_z(cell_z) + breaks_z(cell_z + 1))
      half_z = 0.5d0*(breaks_z(cell_z + 1) - breaks_z(cell_z))
      do qz = 0, SOURCE_QUADRATURE - 1
         X(3) = midpoint_z + half_z*gauss_x(qz)
         span_z = FindSpan(space_z%n, space_z%p, X(3), space_z%knots)
         call DersBasisFuns(span_z, X(3), space_z%p, 0, space_z%knots, bz)
         first_z = span_z - space_z%p
         do cell_y = 1, size(breaks_y) - 1
            midpoint_y = 0.5d0*(breaks_y(cell_y) + breaks_y(cell_y + 1))
            half_y = 0.5d0*(breaks_y(cell_y + 1) - breaks_y(cell_y))
            do qy = 0, SOURCE_QUADRATURE - 1
               X(2) = midpoint_y + half_y*gauss_x(qy)
               span_y = FindSpan(space_y%n, space_y%p, X(2), space_y%knots)
               call DersBasisFuns(span_y, X(2), space_y%p, 0, &
                                  space_y%knots, by)
               first_y = span_y - space_y%p
               do cell_x = 1, size(breaks_x) - 1
                  midpoint_x = 0.5d0*(breaks_x(cell_x) + breaks_x(cell_x + 1))
                  half_x = 0.5d0*(breaks_x(cell_x + 1) - breaks_x(cell_x))
                  do qx = 0, SOURCE_QUADRATURE - 1
                     X(1) = midpoint_x + half_x*gauss_x(qx)
                     span_x = FindSpan(space_x%n, space_x%p, X(1), &
                                       space_x%knots)
                     call DersBasisFuns(span_x, X(1), space_x%p, 0, &
                                        space_x%knots, bx)
                     first_x = span_x - space_x%p
                     value = source(0.d0, gradient, X)
                     if (value == 0.d0) cycle
                     weight = value*half_x*half_y*half_z*gauss_w(qx)* &
                              gauss_w(qy)*gauss_w(qz)
                     do az = 0, space_z%p
                        iz = first_z + az + 1
                        do ay = 0, space_y%p
                           iy = first_y + ay + 1
                           do ax = 0, space_x%p
                              ix = first_x + ax + 1
                              values(ix, iy, iz) = values(ix, iy, iz) + &
                                 weight*bx(ax, 0)*by(ay, 0)*bz(az, 0)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do
   deallocate (breaks_x, breaks_y, breaks_z, bx, by, bz)
   ierr = POLLUTION_DPG_SUCCESS

end subroutine AssembleSourceTensor

subroutine BuildIntegrationBreaks(space, requested_lower, requested_upper, &
                                  breaks)
   implicit none

   type(Space1D), intent(in) :: space
   real(kind=8), intent(in) :: requested_lower, requested_upper
   real(kind=8), allocatable, intent(out) :: breaks(:)
   real(kind=8), allocatable :: candidates(:)
   real(kind=8) :: lower, upper, value, temporary
   integer(kind=4) :: count, index_value, other, unique_count, subdivision

   lower = max(space%knots(1), requested_lower)
   upper = min(space%knots(size(space%knots)), requested_upper)
   if (upper <= lower) then
      allocate (breaks(0))
      return
   end if
   allocate (candidates(size(space%knots) + SOURCE_SUBDIVISIONS + 1))
   count = 2
   candidates(1:2) = (/ lower, upper /)
   do index_value = 1, size(space%knots)
      value = space%knots(index_value)
      if (value > lower .and. value < upper) then
         count = count + 1
         candidates(count) = value
      end if
   end do
   do subdivision = 1, SOURCE_SUBDIVISIONS - 1
      count = count + 1
      candidates(count) = lower + (upper - lower)* &
         real(subdivision, kind=8)/real(SOURCE_SUBDIVISIONS, kind=8)
   end do

   do index_value = 2, count
      temporary = candidates(index_value)
      other = index_value - 1
      do while (other >= 1)
         if (candidates(other) <= temporary) exit
         candidates(other + 1) = candidates(other)
         other = other - 1
      end do
      candidates(other + 1) = temporary
   end do
   unique_count = 1
   do index_value = 2, count
      if (candidates(index_value) > candidates(unique_count)) then
         unique_count = unique_count + 1
         candidates(unique_count) = candidates(index_value)
      end if
   end do
   allocate (breaks(unique_count))
   breaks = candidates(1:unique_count)
   deallocate (candidates)

end subroutine BuildIntegrationBreaks

subroutine InitializeSolution(trial, domain_length, solution)
   implicit none

   type(Space1D), intent(in) :: trial(3)
   real(kind=8), intent(in) :: domain_length
   type(PollutionDPGSolution), intent(out) :: solution
   integer(kind=4) :: axis

   do axis = 1, 3
      solution%n(axis) = trial(axis)%n
      solution%p(axis) = trial(axis)%p
      solution%nelem(axis) = trial(axis)%nelem
   end do
   solution%domain_length = domain_length
   allocate (solution%Ux(size(trial(1)%knots)))
   allocate (solution%Uy(size(trial(2)%knots)))
   allocate (solution%Uz(size(trial(3)%knots)))
   solution%Ux = trial(1)%knots
   solution%Uy = trial(2)%knots
   solution%Uz = trial(3)%knots
   allocate (solution%coefficients(trial(1)%n + 1, trial(2)%n + 1, &
                                   trial(3)%n + 1))
   solution%coefficients = 0.d0

end subroutine InitializeSolution

subroutine ComputeDiagnostics(coefficients, matrices, l2_norm, total_mass, &
                              maximum_abs)
   implicit none

   real(kind=8), intent(in) :: coefficients(:, :, :)
   type(MatrixSet1D), intent(in) :: matrices(3)
   real(kind=8), intent(out) :: l2_norm, total_mass, maximum_abs
   real(kind=8), allocatable :: mass_times_solution(:, :, :)

   call ApplyTensor(matrices(1)%mass, matrices(2)%mass, matrices(3)%mass, &
                    coefficients, mass_times_solution)
   l2_norm = sqrt(max(sum(coefficients*mass_times_solution), 0.d0))
   total_mass = sum(mass_times_solution)
   maximum_abs = maxval(abs(coefficients))
   deallocate (mass_times_solution)

end subroutine ComputeDiagnostics

subroutine EvaluatePollutionDPG(solution, point, value)
   use basis, only: FindSpan, DersBasisFuns

   implicit none

   type(PollutionDPGSolution), intent(in) :: solution
   real(kind=8), intent(in) :: point(3)
   real(kind=8), intent(out) :: value
   real(kind=8) :: bx(0:max(0, min(MAX_POLYNOMIAL_DEGREE, solution%p(1))), 0:0)
   real(kind=8) :: by(0:max(0, min(MAX_POLYNOMIAL_DEGREE, solution%p(2))), 0:0)
   real(kind=8) :: bz(0:max(0, min(MAX_POLYNOMIAL_DEGREE, solution%p(3))), 0:0)
   real(kind=8) :: coordinate(3)
   integer(kind=4) :: span_x, span_y, span_z, first_x, first_y, first_z
   integer(kind=4) :: ax, ay, az, ix, iy, iz

   value = 0.d0
   if (.not. allocated(solution%coefficients)) return
   if (.not. allocated(solution%Ux) .or. .not. allocated(solution%Uy) .or. &
       .not. allocated(solution%Uz)) return
   if (any(solution%p < 0) .or. &
       any(solution%p > MAX_POLYNOMIAL_DEGREE)) return
   if (any(solution%n < solution%p)) return
   coordinate(1) = max(solution%Ux(1), &
                       min(solution%Ux(size(solution%Ux)), point(1)))
   coordinate(2) = max(solution%Uy(1), &
                       min(solution%Uy(size(solution%Uy)), point(2)))
   coordinate(3) = max(solution%Uz(1), &
                       min(solution%Uz(size(solution%Uz)), point(3)))
   span_x = FindSpan(solution%n(1), solution%p(1), coordinate(1), solution%Ux)
   span_y = FindSpan(solution%n(2), solution%p(2), coordinate(2), solution%Uy)
   span_z = FindSpan(solution%n(3), solution%p(3), coordinate(3), solution%Uz)
   call DersBasisFuns(span_x, coordinate(1), solution%p(1), 0, solution%Ux, bx)
   call DersBasisFuns(span_y, coordinate(2), solution%p(2), 0, solution%Uy, by)
   call DersBasisFuns(span_z, coordinate(3), solution%p(3), 0, solution%Uz, bz)
   first_x = span_x - solution%p(1)
   first_y = span_y - solution%p(2)
   first_z = span_z - solution%p(3)
   do az = 0, solution%p(3)
      iz = first_z + az + 1
      do ay = 0, solution%p(2)
         iy = first_y + ay + 1
         do ax = 0, solution%p(1)
            ix = first_x + ax + 1
            value = value + solution%coefficients(ix, iy, iz)* &
               bx(ax, 0)*by(ay, 0)*bz(az, 0)
         end do
      end do
   end do
end subroutine EvaluatePollutionDPG

subroutine PrepareEvaluationAxis(n, p, knots, resolution, spacing, &
                                 basis_values, first_dof)
   use basis, only: FindSpan, DersBasisFuns

   implicit none

   integer(kind=4), intent(in) :: n, p, resolution
   real(kind=8), intent(in) :: knots(0:n + p + 1), spacing
   real(kind=8), allocatable, intent(out) :: basis_values(:, :)
   integer(kind=4), allocatable, intent(out) :: first_dof(:)
   real(kind=8), allocatable :: work(:, :)
   real(kind=8) :: coordinate
   integer(kind=4) :: index_value, span

   allocate (basis_values(0:p, 0:resolution), first_dof(0:resolution))
   allocate (work(0:p, 0:0))
   do index_value = 0, resolution
      coordinate = spacing*real(index_value, kind=8)
      span = FindSpan(n, p, coordinate, knots)
      call DersBasisFuns(span, coordinate, p, 0, knots, work)
      basis_values(:, index_value) = work(:, 0)
      first_dof(index_value) = span - p
   end do
   deallocate (work)

end subroutine PrepareEvaluationAxis

subroutine WritePollutionDPGVTI(solution, path, resolution, ierr)
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none

   type(PollutionDPGSolution), intent(in) :: solution
   character(len=*), intent(in) :: path
   integer(kind=4), intent(in) :: resolution
   integer(kind=4), intent(out) :: ierr
   integer(kind=4) :: ix, iy, iz, ax, ay, az, coefficient_x
   integer(kind=4) :: coefficient_y, coefficient_z
   integer(kind=4) :: unit_number, io_status, close_status, candidate
   integer(kind=8) :: output_points, output_square
   integer(kind=4), allocatable :: first_x(:), first_y(:), first_z(:)
   real(kind=8), allocatable :: basis_x(:, :), basis_y(:, :), basis_z(:, :)
   real(kind=8) :: value, spacing
   logical :: opened

   ierr = POLLUTION_DPG_IO_ERROR
   if (resolution <= 0 .or. .not. allocated(solution%coefficients)) return
   if (.not. ieee_is_finite(solution%domain_length) .or. &
       solution%domain_length <= 0.d0) return
   output_points = int(resolution, kind=8) + 1_8
   if (output_points > MAX_OUTPUT_POINTS/output_points) return
   output_square = output_points*output_points
   if (output_square > MAX_OUTPUT_POINTS/output_points) return
   if (.not. allocated(solution%Ux) .or. .not. allocated(solution%Uy) .or. &
       .not. allocated(solution%Uz)) return
   if (any(solution%p < 0) .or. &
       any(solution%p > MAX_POLYNOMIAL_DEGREE)) return
   if (any(solution%n < solution%p)) return
   if (size(solution%Ux) /= solution%n(1) + solution%p(1) + 2 .or. &
       size(solution%Uy) /= solution%n(2) + solution%p(2) + 2 .or. &
       size(solution%Uz) /= solution%n(3) + solution%p(3) + 2) return
   if (any(shape(solution%coefficients) /= solution%n + 1)) return
   if (.not. all(ieee_is_finite(solution%Ux)) .or. &
       .not. all(ieee_is_finite(solution%Uy)) .or. &
       .not. all(ieee_is_finite(solution%Uz)) .or. &
       .not. all(ieee_is_finite(solution%coefficients))) return
   if (any(solution%Ux(2:) < solution%Ux(:size(solution%Ux) - 1)) .or. &
       any(solution%Uy(2:) < solution%Uy(:size(solution%Uy) - 1)) .or. &
       any(solution%Uz(2:) < solution%Uz(:size(solution%Uz) - 1))) return
   spacing = solution%domain_length/real(resolution, kind=8)
   call PrepareEvaluationAxis(solution%n(1), solution%p(1), solution%Ux, &
                              resolution, spacing, basis_x, first_x)
   call PrepareEvaluationAxis(solution%n(2), solution%p(2), solution%Uy, &
                              resolution, spacing, basis_y, first_y)
   call PrepareEvaluationAxis(solution%n(3), solution%p(3), solution%Uz, &
                              resolution, spacing, basis_z, first_z)
   unit_number = -1
   do candidate = 20, 99
      inquire (unit=candidate, opened=opened)
      if (.not. opened) then
         unit_number = candidate
         exit
      end if
   end do
   if (unit_number < 0) return
   open (unit=unit_number, file=path, status='replace', action='write', &
         form='formatted', iostat=io_status)
   if (io_status /= 0) return

   write (unit_number, '(A)', iostat=io_status) '<?xml version="1.0"?>'
   if (io_status /= 0) go to 900
   write (unit_number, '(A)', iostat=io_status) &
      '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'
   if (io_status /= 0) go to 900
   write (unit_number, '(A,I0,A,I0,A,I0,A,3(ES24.16,1X),A)', &
          iostat=io_status) &
      '  <ImageData WholeExtent="0 ', resolution, ' 0 ', resolution, &
      ' 0 ', resolution, '" Origin="0 0 0" Spacing="', &
      spacing, spacing, spacing, '">'
   if (io_status /= 0) go to 900
   write (unit_number, '(A,I0,A,I0,A,I0,A)', iostat=io_status) &
      '    <Piece Extent="0 ', resolution, ' 0 ', resolution, &
      ' 0 ', resolution, '">'
   if (io_status /= 0) go to 900
   write (unit_number, '(A)', iostat=io_status) &
      '      <PointData Scalars="Result">'
   if (io_status /= 0) go to 900
   write (unit_number, '(A)', iostat=io_status) &
      '        <DataArray Name="Result" type="Float64" ' // &
      'NumberOfComponents="1" format="ascii">'
   if (io_status /= 0) go to 900
   do iz = 0, resolution
      do iy = 0, resolution
         do ix = 0, resolution
            value = 0.d0
            do az = 0, solution%p(3)
               coefficient_z = first_z(iz) + az + 1
               do ay = 0, solution%p(2)
                  coefficient_y = first_y(iy) + ay + 1
                  do ax = 0, solution%p(1)
                     coefficient_x = first_x(ix) + ax + 1
                     value = value + solution%coefficients( &
                        coefficient_x, coefficient_y, coefficient_z)* &
                        basis_x(ax, ix)*basis_y(ay, iy)*basis_z(az, iz)
                  end do
               end do
            end do
            write (unit_number, '(ES24.16)', iostat=io_status) value
            if (io_status /= 0) go to 900
         end do
      end do
   end do
   write (unit_number, '(A)', iostat=io_status) '        </DataArray>'
   if (io_status /= 0) go to 900
   write (unit_number, '(A)', iostat=io_status) '      </PointData>'
   if (io_status /= 0) go to 900
   write (unit_number, '(A)', iostat=io_status) '      <CellData/>'
   if (io_status /= 0) go to 900
   write (unit_number, '(A)', iostat=io_status) '    </Piece>'
   if (io_status /= 0) go to 900
   write (unit_number, '(A)', iostat=io_status) '  </ImageData>'
   if (io_status /= 0) go to 900
   write (unit_number, '(A)', iostat=io_status) '</VTKFile>'
900 continue
   close (unit_number, iostat=close_status)
   if (io_status == 0 .and. close_status == 0) ierr = POLLUTION_DPG_SUCCESS

end subroutine WritePollutionDPGVTI

subroutine BroadcastStats(stats, ierr)
   use mpi
   implicit none

   type(PollutionDPGStats), intent(inout) :: stats
   integer(kind=4), intent(out) :: ierr
   integer(kind=4) :: mpi_ierr, integer_values(1)
   real(kind=8) :: real_values(10)

   integer_values(1) = stats%completed_steps
   real_values = 0.d0
   real_values(1) = stats%time
   real_values(2) = stats%l2_norm
   real_values(3) = stats%total_mass
   real_values(4) = stats%maximum_abs
   real_values(5) = stats%source_integral
   real_values(6) = stats%maximum_relative_residual
   real_values(7:9) = stats%wind
   call MPI_Bcast(integer_values, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr == MPI_SUCCESS) &
      call MPI_Bcast(real_values, 10, MPI_DOUBLE_PRECISION, 0, &
                     MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) then
      ierr = POLLUTION_DPG_MPI_ERROR
      return
   end if
   stats%completed_steps = integer_values(1)
   stats%time = real_values(1)
   stats%l2_norm = real_values(2)
   stats%total_mass = real_values(3)
   stats%maximum_abs = real_values(4)
   stats%source_integral = real_values(5)
   stats%maximum_relative_residual = real_values(6)
   stats%wind = real_values(7:9)
   ierr = POLLUTION_DPG_SUCCESS

end subroutine BroadcastStats

subroutine BroadcastSolution(solution, ierr)
   use parallelism, only: MYRANK
   use mpi
   implicit none

   type(PollutionDPGSolution), intent(inout) :: solution
   integer(kind=4), intent(out) :: ierr
   integer(kind=4) :: mpi_ierr, metadata(9), coefficient_count
   real(kind=8) :: domain_length

   if (MYRANK == 0) then
      metadata(1:3) = solution%n
      metadata(4:6) = solution%p
      metadata(7:9) = solution%nelem
      domain_length = solution%domain_length
   end if
   call MPI_Bcast(metadata, 9, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr == MPI_SUCCESS) &
      call MPI_Bcast(domain_length, 1, MPI_DOUBLE_PRECISION, 0, &
                     MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) then
      ierr = POLLUTION_DPG_MPI_ERROR
      return
   end if
   if (MYRANK /= 0) then
      solution%n = metadata(1:3)
      solution%p = metadata(4:6)
      solution%nelem = metadata(7:9)
      solution%domain_length = domain_length
      allocate (solution%Ux(solution%n(1) + solution%p(1) + 2))
      allocate (solution%Uy(solution%n(2) + solution%p(2) + 2))
      allocate (solution%Uz(solution%n(3) + solution%p(3) + 2))
      allocate (solution%coefficients(solution%n(1) + 1, &
                                      solution%n(2) + 1, &
                                      solution%n(3) + 1))
   end if
   call MPI_Bcast(solution%Ux, size(solution%Ux), MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr == MPI_SUCCESS) &
      call MPI_Bcast(solution%Uy, size(solution%Uy), MPI_DOUBLE_PRECISION, 0, &
                     MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr == MPI_SUCCESS) &
      call MPI_Bcast(solution%Uz, size(solution%Uz), MPI_DOUBLE_PRECISION, 0, &
                     MPI_COMM_WORLD, mpi_ierr)
   coefficient_count = size(solution%coefficients)
   if (mpi_ierr == MPI_SUCCESS) &
      call MPI_Bcast(solution%coefficients, coefficient_count, &
                     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) then
      ierr = POLLUTION_DPG_MPI_ERROR
   else
      ierr = POLLUTION_DPG_SUCCESS
   end if

end subroutine BroadcastSolution

subroutine CleanupPollutionDPGSolution(solution)
   implicit none
   type(PollutionDPGSolution), intent(inout) :: solution

   if (allocated(solution%Ux)) deallocate (solution%Ux)
   if (allocated(solution%Uy)) deallocate (solution%Uy)
   if (allocated(solution%Uz)) deallocate (solution%Uz)
   if (allocated(solution%coefficients)) deallocate (solution%coefficients)
   solution%n = -1
   solution%p = -1
   solution%nelem = 0
   solution%domain_length = 0.d0

end subroutine CleanupPollutionDPGSolution

end module pollution_dpg_solver
