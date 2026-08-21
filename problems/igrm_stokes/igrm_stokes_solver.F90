!------------------------------------------------------------------------------
!> @file igrm_stokes_solver.F90
!> @brief Problem-local direct DGiGRM solver for stationary 3D Stokes flow.
!------------------------------------------------------------------------------
module igrm_stokes_solver

   implicit none
   private

   integer(kind=4), parameter, public :: IGRM_STOKES_SUCCESS = 0
   integer(kind=4), parameter, public :: IGRM_STOKES_INVALID_CONFIGURATION = -7401
   integer(kind=4), parameter, public :: IGRM_STOKES_SYSTEM_TOO_LARGE = -7402
   integer(kind=4), parameter, public :: IGRM_STOKES_RESIDUAL_TOO_LARGE = -7403
   integer(kind=4), parameter, public :: IGRM_STOKES_IO_ERROR = -7404
   integer(kind=4), parameter, public :: IGRM_STOKES_MPI_ERROR = -7405

   integer(kind=4), parameter :: VELOCITY_X = 1
   integer(kind=4), parameter :: VELOCITY_Y = 2
   integer(kind=4), parameter :: VELOCITY_Z = 3
   integer(kind=4), parameter :: PRESSURE_FIELD = 4

   abstract interface
      function VectorField3D(point) result(value)
         real(kind=8), intent(in) :: point(3)
         real(kind=8) :: value(3)
      end function VectorField3D

      function ScalarField3D(point) result(value)
         real(kind=8), intent(in) :: point(3)
         real(kind=8) :: value
      end function ScalarField3D
   end interface

   public :: VectorField3D, ScalarField3D

   type, public :: IGRMStokesStats
      integer(kind=4) :: system_size = 0
      integer(kind=8) :: nonzeros = 0_8
      real(kind=8) :: residual_rms = huge(0.d0)
      real(kind=8) :: relative_residual = huge(0.d0)
   end type IGRMStokesStats

   type, public :: IGRMStokesSolution
      integer(kind=4) :: n(3) = -1
      integer(kind=4) :: p(3) = -1
      integer(kind=4) :: nelem(3) = 0
      real(kind=8), allocatable :: Ux(:)
      real(kind=8), allocatable :: Uy(:)
      real(kind=8), allocatable :: Uz(:)
      ! Columns 1:3 are velocity components; column 4 is pressure.
      real(kind=8), allocatable :: coefficients(:, :)
   end type IGRMStokesSolution

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

   type :: Space3D
      type(Space1D) :: direction(3)
      integer(kind=4) :: n(3) = -1
      integer(kind=4) :: p(3) = -1
      integer(kind=4) :: nelem(3) = 0
   end type Space3D

   type :: Matrix1D
      real(kind=8), allocatable :: mass(:, :)
      real(kind=8), allocatable :: stiffness(:, :)
      real(kind=8), allocatable :: derivative_row(:, :)
      real(kind=8), allocatable :: derivative_column(:, :)
      real(kind=8), allocatable :: element_mass(:, :, :)
   end type Matrix1D

   public :: AssembleIGRMStokes3DSystem
   public :: SolveIGRMStokes3D
   public :: EvaluateIGRMStokes
   public :: ComputeIGRMStokesErrors
   public :: WriteIGRMStokesVTI
   public :: CleanupIGRMStokesSolution

contains

!------------------------------------------------------------------------------
!> Solve the symmetric DGiGRM saddle system on rank zero and replicate the
!> four trial fields.  A final Lagrange multiplier enforces integral(p)=0.
!------------------------------------------------------------------------------
subroutine SolveIGRMStokes3D(nelem, p_test, p_trial, viscosity, &
                             penalty_factor, residual_tolerance, forcing_x, &
                             forcing_y, forcing_z, boundary_velocity, &
                             solution, stats, ierr)
   use Interfaces, only: forcing_fun
   use mumps_solver, only: SolveOneDirection
   use parallelism, only: MYRANK
   use sparse, only: sparse_matrix, clear_matrix
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mpi

   implicit none

   integer(kind=4), intent(in) :: nelem(3), p_test(3), p_trial(3)
   real(kind=8), intent(in) :: viscosity, penalty_factor, residual_tolerance
   procedure(forcing_fun) :: forcing_x, forcing_y, forcing_z
   procedure(VectorField3D) :: boundary_velocity
   type(IGRMStokesSolution), intent(out) :: solution
   type(IGRMStokesStats), intent(out) :: stats
   integer(kind=4), intent(out) :: ierr

   type(sparse_matrix), pointer :: matrix
   real(kind=8), allocatable :: rhs(:), buffer(:, :), global_fields(:, :)
   integer(kind=8) :: test_scalar_64, trial_scalar_64, system_size_64
   integer(kind=4) :: test_scalar, trial_scalar, system_size
   integer(kind=4) :: solver_ierr, mpi_ierr, field

   nullify (matrix)
   stats = IGRMStokesStats()
   call CleanupIGRMStokesSolution(solution)

   call ValidateConfiguration(nelem, p_test, p_trial, viscosity, &
                              penalty_factor, residual_tolerance, &
                              test_scalar_64, trial_scalar_64, &
                              system_size_64, ierr)
   if (ierr /= IGRM_STOKES_SUCCESS) return

   test_scalar = int(test_scalar_64, kind=4)
   trial_scalar = int(trial_scalar_64, kind=4)
   system_size = int(system_size_64, kind=4)
   stats%system_size = system_size
   allocate (global_fields(trial_scalar, 4))
   global_fields = 0.d0

   if (MYRANK == 0) then
      call AssembleSystem(nelem, p_test, p_trial, viscosity, penalty_factor, &
                          forcing_x, forcing_y, forcing_z, boundary_velocity, &
                          matrix, rhs, ierr)
      if (ierr == IGRM_STOKES_SUCCESS) then
         stats%nonzeros = matrix%total_entries
         allocate (buffer(system_size, 1))
         buffer(:, 1) = rhs
         call SolveOneDirection(buffer, 1, system_size - 1, 0, matrix, &
                                solver_ierr)
         if (solver_ierr /= 0) then
            ierr = solver_ierr
         else
            call ComputeResidual(matrix, rhs, buffer(:, 1), &
                                 stats%residual_rms, stats%relative_residual)
            if (.not. ieee_is_finite(stats%residual_rms) .or. &
                .not. ieee_is_finite(stats%relative_residual) .or. &
                stats%relative_residual > residual_tolerance) then
               ierr = IGRM_STOKES_RESIDUAL_TOO_LARGE
            else
               do field = 1, 4
                  global_fields(:, field) = buffer( &
                     4*test_scalar + (field - 1)*trial_scalar + 1: &
                     4*test_scalar + field*trial_scalar, 1)
               end do
            end if
         end if
      end if

      if (associated(matrix)) call clear_matrix(matrix)
      if (allocated(rhs)) deallocate (rhs)
      if (allocated(buffer)) deallocate (buffer)
   end if

   call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) ierr = IGRM_STOKES_MPI_ERROR
   if (ierr /= IGRM_STOKES_SUCCESS) then
      deallocate (global_fields)
      return
   end if

   call MPI_Bcast(stats%system_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) ierr = IGRM_STOKES_MPI_ERROR
   call MPI_Bcast(stats%nonzeros, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) ierr = IGRM_STOKES_MPI_ERROR
   call MPI_Bcast(stats%residual_rms, 1, MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) ierr = IGRM_STOKES_MPI_ERROR
   call MPI_Bcast(stats%relative_residual, 1, MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) ierr = IGRM_STOKES_MPI_ERROR
   call MPI_Bcast(global_fields, 4*trial_scalar, MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) ierr = IGRM_STOKES_MPI_ERROR
   if (ierr /= IGRM_STOKES_SUCCESS) then
      deallocate (global_fields)
      return
   end if

   call InitializeSolution(nelem, p_trial, global_fields, solution, ierr)
   deallocate (global_fields)

end subroutine SolveIGRMStokes3D

!------------------------------------------------------------------------------
!> Public assembly entry for focused matrix and facet tests.
!------------------------------------------------------------------------------
subroutine AssembleIGRMStokes3DSystem(nelem, p_test, p_trial, viscosity, &
                                      penalty_factor, forcing_x, forcing_y, &
                                      forcing_z, boundary_velocity, matrix, &
                                      rhs, ierr)
   use Interfaces, only: forcing_fun
   use sparse, only: sparse_matrix

   implicit none

   integer(kind=4), intent(in) :: nelem(3), p_test(3), p_trial(3)
   real(kind=8), intent(in) :: viscosity, penalty_factor
   procedure(forcing_fun) :: forcing_x, forcing_y, forcing_z
   procedure(VectorField3D) :: boundary_velocity
   type(sparse_matrix), pointer, intent(out) :: matrix
   real(kind=8), allocatable, intent(out) :: rhs(:)
   integer(kind=4), intent(out) :: ierr

   integer(kind=8) :: test_scalar, trial_scalar, system_size

   nullify (matrix)
   call ValidateConfiguration(nelem, p_test, p_trial, viscosity, &
                              penalty_factor, 1.d0, test_scalar, &
                              trial_scalar, system_size, ierr)
   if (ierr == IGRM_STOKES_SUCCESS) then
      call AssembleSystem(nelem, p_test, p_trial, viscosity, penalty_factor, &
                          forcing_x, forcing_y, forcing_z, boundary_velocity, &
                          matrix, rhs, ierr)
   end if

end subroutine AssembleIGRMStokes3DSystem

!------------------------------------------------------------------------------
!> Build both tensor spaces and all volume/facet contributions.
!------------------------------------------------------------------------------
subroutine AssembleSystem(nelem, p_test, p_trial, viscosity, penalty_factor, &
                          forcing_x, forcing_y, forcing_z, boundary_velocity, &
                          matrix, rhs, ierr)
   use Interfaces, only: forcing_fun
   use sparse, only: sparse_matrix, initialize_sparse

   implicit none

   integer(kind=4), intent(in) :: nelem(3), p_test(3), p_trial(3)
   real(kind=8), intent(in) :: viscosity, penalty_factor
   procedure(forcing_fun) :: forcing_x, forcing_y, forcing_z
   procedure(VectorField3D) :: boundary_velocity
   type(sparse_matrix), pointer, intent(out) :: matrix
   real(kind=8), allocatable, intent(out) :: rhs(:)
   integer(kind=4), intent(out) :: ierr

   type(Space3D) :: test_space, trial_space
   type(Matrix1D) :: gram_x, gram_y, gram_z
   type(Matrix1D) :: mixed_x, mixed_y, mixed_z
   integer(kind=4) :: ng, test_scalar, trial_scalar, system_size

   nullify (matrix)
   ierr = IGRM_STOKES_SUCCESS
   ng = max(2, max(maxval(p_test), maxval(p_trial)) + 1)
   call BuildSpace3D(nelem, p_test, -1, ng, test_space, ierr)
   if (ierr /= IGRM_STOKES_SUCCESS) return
   call BuildSpace3D(nelem, p_trial, 1, ng, trial_space, ierr)
   if (ierr /= IGRM_STOKES_SUCCESS) return

   test_scalar = TensorSize(test_space%n)
   trial_scalar = TensorSize(trial_space%n)
   system_size = 4*test_scalar + 4*trial_scalar + 1
   call initialize_sparse(system_size, system_size, matrix)
   allocate (rhs(system_size))
   rhs = 0.d0
   call AssembleMatrix1D(test_space%direction(1), &
                         test_space%direction(1), gram_x)
   call AssembleMatrix1D(test_space%direction(2), &
                         test_space%direction(2), gram_y)
   call AssembleMatrix1D(test_space%direction(3), &
                         test_space%direction(3), gram_z)
   call AssembleMatrix1D(test_space%direction(1), &
                         trial_space%direction(1), mixed_x)
   call AssembleMatrix1D(test_space%direction(2), &
                         trial_space%direction(2), mixed_y)
   call AssembleMatrix1D(test_space%direction(3), &
                         trial_space%direction(3), mixed_z)

   call AssembleVolumeBlocks(test_space, trial_space, gram_x, gram_y, &
                             gram_z, mixed_x, mixed_y, mixed_z, viscosity, &
                             matrix)
   call AssembleFacetBlocks(test_space, trial_space, gram_x, gram_y, &
                            gram_z, mixed_x, mixed_y, mixed_z, viscosity, &
                            penalty_factor, matrix)
   call AssembleVolumeLoad(test_space, forcing_x, forcing_y, forcing_z, rhs)
   call AssembleBoundaryLoad(test_space, p_trial, viscosity, penalty_factor, &
                             boundary_velocity, rhs)
   call AssemblePressureGauge(trial_space, matrix)

end subroutine AssembleSystem

!------------------------------------------------------------------------------
!> Validate sizes without overflowing the default integer used by sparse/MUMPS.
!------------------------------------------------------------------------------
subroutine ValidateConfiguration(nelem, p_test, p_trial, viscosity, &
                                 penalty_factor, residual_tolerance, &
                                 test_scalar, trial_scalar, system_size, ierr)
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none

   integer(kind=4), intent(in) :: nelem(3), p_test(3), p_trial(3)
   real(kind=8), intent(in) :: viscosity, penalty_factor, residual_tolerance
   integer(kind=8), intent(out) :: test_scalar, trial_scalar, system_size
   integer(kind=4), intent(out) :: ierr
   integer(kind=8) :: test_extent(3), trial_extent(3), trial_continuity(3)
   logical :: valid

   test_scalar = 0_8
   trial_scalar = 0_8
   system_size = 0_8
   ierr = IGRM_STOKES_INVALID_CONFIGURATION

   if (any(nelem <= 0) .or. any(p_test < 1) .or. any(p_trial < 1)) return
   if (any(p_test < p_trial) .or. any(p_test > 9) .or. any(p_trial > 9)) return
   if (.not. ieee_is_finite(viscosity) .or. &
       .not. ieee_is_finite(penalty_factor) .or. &
       .not. ieee_is_finite(residual_tolerance)) return
   if (viscosity <= 0.d0 .or. penalty_factor <= 0.d0 .or. &
       residual_tolerance <= 0.d0) return

   test_extent = int(nelem, kind=8)*(int(p_test, kind=8) + 1_8)
   trial_continuity = min(1_8, int(p_trial, kind=8) - 1_8)
   trial_extent = int(nelem, kind=8)* &
      (int(p_trial, kind=8) - trial_continuity) + &
      trial_continuity + 1_8
   call SafeProduct(test_extent, test_scalar, valid)
   if (.not. valid) then
      ierr = IGRM_STOKES_SYSTEM_TOO_LARGE
      return
   end if
   call SafeProduct(trial_extent, trial_scalar, valid)
   if (.not. valid) then
      ierr = IGRM_STOKES_SYSTEM_TOO_LARGE
      return
   end if

   if (test_scalar > (huge(0_8) - 1_8)/4_8 .or. &
       trial_scalar > (huge(0_8) - 1_8 - 4_8*test_scalar)/4_8) then
      ierr = IGRM_STOKES_SYSTEM_TOO_LARGE
      return
   end if
   system_size = 4_8*test_scalar + 4_8*trial_scalar + 1_8
   if (system_size > int(huge(0_4), kind=8)) then
      ierr = IGRM_STOKES_SYSTEM_TOO_LARGE
      return
   end if
   ierr = IGRM_STOKES_SUCCESS

end subroutine ValidateConfiguration

subroutine SafeProduct(extents, product_value, valid)
   implicit none
   integer(kind=8), intent(in) :: extents(3)
   integer(kind=8), intent(out) :: product_value
   logical, intent(out) :: valid
   integer(kind=4) :: direction

   product_value = 1_8
   valid = .true.
   do direction = 1, 3
      if (extents(direction) <= 0_8 .or. &
          product_value > huge(product_value)/extents(direction)) then
         product_value = 0_8
         valid = .false.
         return
      end if
      product_value = product_value*extents(direction)
   end do
end subroutine SafeProduct

!------------------------------------------------------------------------------
!> Construct C^-1 test spaces or trial spaces with selected continuity.
!------------------------------------------------------------------------------
subroutine BuildSpace3D(nelem, degree, continuity_selector, ng, space, ierr)
   implicit none
   integer(kind=4), intent(in) :: nelem(3), degree(3)
   integer(kind=4), intent(in) :: continuity_selector, ng
   type(Space3D), intent(out) :: space
   integer(kind=4), intent(out) :: ierr
   integer(kind=4) :: direction, continuity

   ierr = IGRM_STOKES_SUCCESS
   space%nelem = nelem
   space%p = degree
   do direction = 1, 3
      if (continuity_selector < 0) then
         continuity = -1
      else
         continuity = min(continuity_selector, degree(direction) - 1)
      end if
      call BuildSpace1D(nelem(direction), degree(direction), continuity, ng, &
                        space%direction(direction), ierr)
      if (ierr /= IGRM_STOKES_SUCCESS) return
      space%n(direction) = space%direction(direction)%n
   end do
end subroutine BuildSpace3D

subroutine BuildSpace1D(nelem, degree, continuity, ng, space, ierr)
   use knot_vector, only: PrepareKnot
   use basis, only: BasisData

   implicit none
   integer(kind=4), intent(in) :: nelem, degree, continuity, ng
   type(Space1D), intent(out) :: space
   integer(kind=4), intent(out) :: ierr
   integer(kind=4) :: generated_nelem, block_size, m

   ierr = IGRM_STOKES_SUCCESS
   space%n = nelem + degree - 1
   space%p = degree
   space%ng = ng
   if (continuity < degree - 1) then
      block_size = 1
   else
      block_size = nelem
   end if
   call PrepareKnot(space%n, degree, block_size, continuity, space%knots, &
                    generated_nelem)
   if (generated_nelem /= nelem) then
      ierr = IGRM_STOKES_INVALID_CONFIGURATION
      return
   end if
   space%nelem = generated_nelem
   m = space%n + degree + 1
   allocate (space%offset(nelem), space%jacobian(nelem))
   allocate (space%weights(ng), space%points(ng, nelem))
   allocate (space%basis(0:1, 0:degree, ng, nelem))
   call BasisData(degree, m, space%knots, 1, ng, nelem, space%offset, &
                  space%jacobian, space%weights, space%points, space%basis)

end subroutine BuildSpace1D

!------------------------------------------------------------------------------
!> Assemble global and element-local 1D products for a row/column space pair.
!------------------------------------------------------------------------------
subroutine AssembleMatrix1D(row_space, column_space, matrix)
   implicit none
   type(Space1D), intent(in) :: row_space, column_space
   type(Matrix1D), intent(out) :: matrix
   integer(kind=4) :: element, q, a, b, row, column
   real(kind=8) :: scale, row_value, column_value

   allocate (matrix%mass(0:row_space%n, 0:column_space%n))
   allocate (matrix%stiffness(0:row_space%n, 0:column_space%n))
   allocate (matrix%derivative_row(0:row_space%n, 0:column_space%n))
   allocate (matrix%derivative_column(0:row_space%n, 0:column_space%n))
   allocate (matrix%element_mass(0:row_space%n, 0:column_space%n, &
                                 row_space%nelem))
   matrix%mass = 0.d0
   matrix%stiffness = 0.d0
   matrix%derivative_row = 0.d0
   matrix%derivative_column = 0.d0
   matrix%element_mass = 0.d0

   do element = 1, row_space%nelem
      do q = 1, row_space%ng
         scale = row_space%jacobian(element)*row_space%weights(q)
         do a = 0, row_space%p
            row = row_space%offset(element) + a
            row_value = row_space%basis(0, a, q, element)
            do b = 0, column_space%p
               column = column_space%offset(element) + b
               column_value = column_space%basis(0, b, q, element)
               matrix%mass(row, column) = matrix%mass(row, column) + &
                  row_value*column_value*scale
               matrix%element_mass(row, column, element) = &
                  matrix%element_mass(row, column, element) + &
                  row_value*column_value*scale
               matrix%stiffness(row, column) = matrix%stiffness(row, column) + &
                  row_space%basis(1, a, q, element)* &
                  column_space%basis(1, b, q, element)*scale
               matrix%derivative_row(row, column) = &
                  matrix%derivative_row(row, column) + &
                  row_space%basis(1, a, q, element)*column_value*scale
               matrix%derivative_column(row, column) = &
                  matrix%derivative_column(row, column) + &
                  row_value*column_space%basis(1, b, q, element)*scale
            end do
         end do
      end do
   end do

end subroutine AssembleMatrix1D

!------------------------------------------------------------------------------
!> Volume Gram and Stokes operator blocks.
!------------------------------------------------------------------------------
subroutine AssembleVolumeBlocks(test_space, trial_space, gx, gy, gz, mx, my, &
                                mz, viscosity, matrix)
   use sparse, only: sparse_matrix

   implicit none
   type(Space3D), intent(in) :: test_space, trial_space
   type(Matrix1D), intent(in) :: gx, gy, gz, mx, my, mz
   real(kind=8), intent(in) :: viscosity
   type(sparse_matrix), pointer, intent(inout) :: matrix
   integer(kind=4) :: ix, iy, iz, jx, jy, jz, field
   integer(kind=4) :: row_scalar, column_scalar, nt, nu
   real(kind=8) :: mass_value, laplace_value, coupling

   nt = TensorSize(test_space%n)
   nu = TensorSize(trial_space%n)

   do iz = 0, test_space%n(3)
      do iy = 0, test_space%n(2)
         do ix = 0, test_space%n(1)
            row_scalar = LinearIndex3D(ix, iy, iz, test_space%n)
            do jz = 0, test_space%n(3)
               do jy = 0, test_space%n(2)
                  do jx = 0, test_space%n(1)
                     column_scalar = LinearIndex3D(jx, jy, jz, test_space%n)
                     mass_value = gx%mass(ix, jx)*gy%mass(iy, jy)* &
                                  gz%mass(iz, jz)
                     laplace_value = &
                        gx%stiffness(ix, jx)*gy%mass(iy, jy)*gz%mass(iz, jz) + &
                        gx%mass(ix, jx)*gy%stiffness(iy, jy)*gz%mass(iz, jz) + &
                        gx%mass(ix, jx)*gy%mass(iy, jy)*gz%stiffness(iz, jz)
                     do field = VELOCITY_X, VELOCITY_Z
                        call AddMatrixEntry(matrix, &
                           TestIndex(field, row_scalar, nt), &
                           TestIndex(field, column_scalar, nt), laplace_value)
                     end do
                     call AddMatrixEntry(matrix, &
                        TestIndex(PRESSURE_FIELD, row_scalar, nt), &
                        TestIndex(PRESSURE_FIELD, column_scalar, nt), mass_value)
                  end do
               end do
            end do
         end do
      end do
   end do

   do iz = 0, test_space%n(3)
      do iy = 0, test_space%n(2)
         do ix = 0, test_space%n(1)
            row_scalar = LinearIndex3D(ix, iy, iz, test_space%n)
            do jz = 0, trial_space%n(3)
               do jy = 0, trial_space%n(2)
                  do jx = 0, trial_space%n(1)
                     column_scalar = LinearIndex3D(jx, jy, jz, trial_space%n)
                     laplace_value = viscosity*( &
                        mx%stiffness(ix, jx)*my%mass(iy, jy)*mz%mass(iz, jz) + &
                        mx%mass(ix, jx)*my%stiffness(iy, jy)*mz%mass(iz, jz) + &
                        mx%mass(ix, jx)*my%mass(iy, jy)*mz%stiffness(iz, jz))
                     do field = VELOCITY_X, VELOCITY_Z
                        call AddSymmetricB(matrix, &
                           TestIndex(field, row_scalar, nt), &
                           TrialIndex(field, column_scalar, nt, nu), &
                           laplace_value)
                     end do

                     coupling = -mx%derivative_row(ix, jx)* &
                                 my%mass(iy, jy)*mz%mass(iz, jz)
                     call AddSymmetricB(matrix, &
                        TestIndex(VELOCITY_X, row_scalar, nt), &
                        TrialIndex(PRESSURE_FIELD, column_scalar, nt, nu), coupling)
                     coupling = -mx%mass(ix, jx)* &
                                 my%derivative_row(iy, jy)*mz%mass(iz, jz)
                     call AddSymmetricB(matrix, &
                        TestIndex(VELOCITY_Y, row_scalar, nt), &
                        TrialIndex(PRESSURE_FIELD, column_scalar, nt, nu), coupling)
                     coupling = -mx%mass(ix, jx)*my%mass(iy, jy)* &
                                 mz%derivative_row(iz, jz)
                     call AddSymmetricB(matrix, &
                        TestIndex(VELOCITY_Z, row_scalar, nt), &
                        TrialIndex(PRESSURE_FIELD, column_scalar, nt, nu), coupling)

                     coupling = mx%derivative_column(ix, jx)* &
                                my%mass(iy, jy)*mz%mass(iz, jz)
                     call AddSymmetricB(matrix, &
                        TestIndex(PRESSURE_FIELD, row_scalar, nt), &
                        TrialIndex(VELOCITY_X, column_scalar, nt, nu), coupling)
                     coupling = mx%mass(ix, jx)* &
                                my%derivative_column(iy, jy)*mz%mass(iz, jz)
                     call AddSymmetricB(matrix, &
                        TestIndex(PRESSURE_FIELD, row_scalar, nt), &
                        TrialIndex(VELOCITY_Y, column_scalar, nt, nu), coupling)
                     coupling = mx%mass(ix, jx)*my%mass(iy, jy)* &
                                mz%derivative_column(iz, jz)
                     call AddSymmetricB(matrix, &
                        TestIndex(PRESSURE_FIELD, row_scalar, nt), &
                        TrialIndex(VELOCITY_Z, column_scalar, nt, nu), coupling)
                  end do
               end do
            end do
         end do
      end do
   end do

end subroutine AssembleVolumeBlocks

!------------------------------------------------------------------------------
!> Assemble all interior and boundary facet forms used by upstream
!> DGiGRM_stokes_3D.  jump(v)=n(v_left-v_right), and boundary jumps equal v.
!------------------------------------------------------------------------------
subroutine AssembleFacetBlocks(test_space, trial_space, gx, gy, gz, mx, my, &
                               mz, viscosity, penalty_factor, matrix)
   use sparse, only: sparse_matrix

   implicit none
   type(Space3D), intent(in) :: test_space, trial_space
   type(Matrix1D), intent(in) :: gx, gy, gz, mx, my, mz
   real(kind=8), intent(in) :: viscosity, penalty_factor
   type(sparse_matrix), pointer, intent(inout) :: matrix

   real(kind=8), allocatable :: test_avg(:), test_avg_derivative(:)
   real(kind=8), allocatable :: test_jump(:)
   real(kind=8), allocatable :: trial_avg(:), trial_avg_derivative(:)
   real(kind=8), allocatable :: trial_jump(:)
   integer(kind=4) :: direction, tangent1, tangent2, face, element1, element2
   integer(kind=4) :: inormal, jnormal, a1, a2, b1, b2, field
   integer(kind=4) :: row_coordinates(3), column_coordinates(3)
   integer(kind=4) :: row_scalar, column_scalar, nt, nu
   integer(kind=4) :: row1, row2, column1, column2
   real(kind=8) :: normal, h, eta, tangential, value
   logical :: interior

   nt = TensorSize(test_space%n)
   nu = TensorSize(trial_space%n)

   do direction = 1, 3
      call TangentialDirections(direction, tangent1, tangent2)
      allocate (test_avg(0:test_space%n(direction)))
      allocate (test_avg_derivative(0:test_space%n(direction)))
      allocate (test_jump(0:test_space%n(direction)))
      allocate (trial_avg(0:trial_space%n(direction)))
      allocate (trial_avg_derivative(0:trial_space%n(direction)))
      allocate (trial_jump(0:trial_space%n(direction)))

      eta = penalty_factor*real((trial_space%p(direction) + 1)* &
                                (trial_space%p(direction) + 2), kind=8)
      do face = 0, test_space%nelem(direction)
         if (face == 0) then
            normal = -1.d0
         else
            normal = 1.d0
         end if
         interior = face > 0 .and. face < test_space%nelem(direction)
         call ComputeTrace1D(test_space%direction(direction), face, normal, &
                             test_avg, test_avg_derivative, test_jump)
         call ComputeTrace1D(trial_space%direction(direction), face, normal, &
                             trial_avg, trial_avg_derivative, trial_jump)

         do element2 = 1, test_space%nelem(tangent2)
            do element1 = 1, test_space%nelem(tangent1)
               h = sqrt((2.d0*test_space%direction(tangent1)% &
                         jacobian(element1))**2 + &
                        (2.d0*test_space%direction(tangent2)% &
                         jacobian(element2))**2)

               ! Facet parts of G for three velocity residuals and pressure.
               do inormal = 0, test_space%n(direction)
                  if (test_jump(inormal) == 0.d0) cycle
                  do jnormal = 0, test_space%n(direction)
                     if (test_jump(jnormal) == 0.d0) cycle
                     do a2 = 0, test_space%p(tangent2)
                        row2 = test_space%direction(tangent2)% &
                               offset(element2) + a2
                        do b2 = 0, test_space%p(tangent2)
                           column2 = test_space%direction(tangent2)% &
                                    offset(element2) + b2
                           do a1 = 0, test_space%p(tangent1)
                              row1 = test_space%direction(tangent1)% &
                                     offset(element1) + a1
                              do b1 = 0, test_space%p(tangent1)
                                 column1 = test_space%direction(tangent1)% &
                                          offset(element1) + b1
                                 tangential = ElementMass(gx, gy, gz, tangent1, &
                                                          row1, column1, element1)* &
                                              ElementMass(gx, gy, gz, tangent2, &
                                                          row2, column2, element2)
                                 if (tangential == 0.d0) cycle
                                 row_coordinates(direction) = inormal
                                 row_coordinates(tangent1) = row1
                                 row_coordinates(tangent2) = row2
                                 column_coordinates(direction) = jnormal
                                 column_coordinates(tangent1) = column1
                                 column_coordinates(tangent2) = column2
                                 row_scalar = LinearIndex3D(row_coordinates(1), &
                                    row_coordinates(2), row_coordinates(3), &
                                    test_space%n)
                                 column_scalar = LinearIndex3D( &
                                    column_coordinates(1), column_coordinates(2), &
                                    column_coordinates(3), test_space%n)
                                 value = test_jump(inormal)*test_jump(jnormal)* &
                                         tangential/h
                                 do field = VELOCITY_X, VELOCITY_Z
                                    call AddMatrixEntry(matrix, &
                                       TestIndex(field, row_scalar, nt), &
                                       TestIndex(field, column_scalar, nt), value)
                                 end do
                                 if (interior) then
                                    value = h*test_jump(inormal)* &
                                            test_jump(jnormal)*tangential
                                    call AddMatrixEntry(matrix, &
                                       TestIndex(PRESSURE_FIELD, row_scalar, nt), &
                                       TestIndex(PRESSURE_FIELD, column_scalar, nt), &
                                       value)
                                 end if
                              end do
                           end do
                        end do
                     end do
                  end do
               end do

               ! Facet parts of B, with row in C^-1 test and column in trial.
               do inormal = 0, test_space%n(direction)
                  if (test_jump(inormal) == 0.d0 .and. &
                      test_avg(inormal) == 0.d0 .and. &
                      test_avg_derivative(inormal) == 0.d0) cycle
                  do jnormal = 0, trial_space%n(direction)
                     if (trial_jump(jnormal) == 0.d0 .and. &
                         trial_avg(jnormal) == 0.d0 .and. &
                         trial_avg_derivative(jnormal) == 0.d0) cycle
                     do a2 = 0, test_space%p(tangent2)
                        row2 = test_space%direction(tangent2)% &
                               offset(element2) + a2
                        do b2 = 0, trial_space%p(tangent2)
                           column2 = trial_space%direction(tangent2)% &
                                    offset(element2) + b2
                           do a1 = 0, test_space%p(tangent1)
                              row1 = test_space%direction(tangent1)% &
                                     offset(element1) + a1
                              do b1 = 0, trial_space%p(tangent1)
                                 column1 = trial_space%direction(tangent1)% &
                                          offset(element1) + b1
                                 tangential = ElementMass(mx, my, mz, tangent1, &
                                                          row1, column1, element1)* &
                                              ElementMass(mx, my, mz, tangent2, &
                                                          row2, column2, element2)
                                 if (tangential == 0.d0) cycle
                                 row_coordinates(direction) = inormal
                                 row_coordinates(tangent1) = row1
                                 row_coordinates(tangent2) = row2
                                 column_coordinates(direction) = jnormal
                                 column_coordinates(tangent1) = column1
                                 column_coordinates(tangent2) = column2
                                 row_scalar = LinearIndex3D(row_coordinates(1), &
                                    row_coordinates(2), row_coordinates(3), &
                                    test_space%n)
                                 column_scalar = LinearIndex3D( &
                                    column_coordinates(1), column_coordinates(2), &
                                    column_coordinates(3), trial_space%n)

                                 value = viscosity*( &
                                    -normal*test_avg_derivative(inormal)* &
                                               trial_jump(jnormal) &
                                    -normal*trial_avg_derivative(jnormal)* &
                                               test_jump(inormal) &
                                    +eta/h*trial_jump(jnormal)* &
                                               test_jump(inormal))*tangential
                                 do field = VELOCITY_X, VELOCITY_Z
                                    call AddSymmetricB(matrix, &
                                       TestIndex(field, row_scalar, nt), &
                                       TrialIndex(field, column_scalar, nt, nu), &
                                       value)
                                 end do

                                 value = trial_avg(jnormal)*normal* &
                                         test_jump(inormal)*tangential
                                 call AddSymmetricB(matrix, &
                                    TestIndex(direction, row_scalar, nt), &
                                    TrialIndex(PRESSURE_FIELD, column_scalar, nt, nu), &
                                    value)
                                 value = -normal*trial_jump(jnormal)* &
                                         test_avg(inormal)*tangential
                                 call AddSymmetricB(matrix, &
                                    TestIndex(PRESSURE_FIELD, row_scalar, nt), &
                                    TrialIndex(direction, column_scalar, nt, nu), &
                                    value)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      deallocate (test_avg, test_avg_derivative, test_jump)
      deallocate (trial_avg, trial_avg_derivative, trial_jump)
   end do

end subroutine AssembleFacetBlocks

!------------------------------------------------------------------------------
!> Evaluate average, normal-signed jump, and averaged physical derivative.
!------------------------------------------------------------------------------
subroutine ComputeTrace1D(space, face, normal, average, average_derivative, jump)
   use basis, only: DersBasisFuns

   implicit none
   type(Space1D), intent(in) :: space
   integer(kind=4), intent(in) :: face
   real(kind=8), intent(in) :: normal
   real(kind=8), intent(out) :: average(0:), average_derivative(0:), jump(0:)
   real(kind=8), allocatable :: left_value(:), right_value(:)
   real(kind=8), allocatable :: left_derivative(:), right_derivative(:)
   real(kind=8) :: basis_values(0:space%p, 0:1), coordinate
   integer(kind=4) :: element, span, a, index_value
   logical :: interior

   allocate (left_value(0:space%n), right_value(0:space%n))
   allocate (left_derivative(0:space%n), right_derivative(0:space%n))
   left_value = 0.d0
   right_value = 0.d0
   left_derivative = 0.d0
   right_derivative = 0.d0
   coordinate = real(face, kind=8)/real(space%nelem, kind=8)
   interior = face > 0 .and. face < space%nelem

   if (face > 0) then
      element = face
      span = space%offset(element) + space%p
      call DersBasisFuns(span, coordinate, space%p, 1, space%knots, &
                         basis_values)
      do a = 0, space%p
         index_value = space%offset(element) + a
         left_value(index_value) = basis_values(a, 0)
         left_derivative(index_value) = basis_values(a, 1)
      end do
   end if
   if (face < space%nelem) then
      element = face + 1
      span = space%offset(element) + space%p
      call DersBasisFuns(span, coordinate, space%p, 1, space%knots, &
                         basis_values)
      do a = 0, space%p
         index_value = space%offset(element) + a
         right_value(index_value) = basis_values(a, 0)
         right_derivative(index_value) = basis_values(a, 1)
      end do
   end if

   if (interior) then
      average = 0.5d0*(left_value + right_value)
      average_derivative = 0.5d0*(left_derivative + right_derivative)
   else
      average = left_value + right_value
      average_derivative = left_derivative + right_derivative
   end if
   jump = normal*(left_value - right_value)
   deallocate (left_value, right_value, left_derivative, right_derivative)

end subroutine ComputeTrace1D

pure subroutine TangentialDirections(direction, tangent1, tangent2)
   implicit none
   integer(kind=4), intent(in) :: direction
   integer(kind=4), intent(out) :: tangent1, tangent2

   select case (direction)
   case (1)
      tangent1 = 2
      tangent2 = 3
   case (2)
      tangent1 = 1
      tangent2 = 3
   case default
      tangent1 = 1
      tangent2 = 2
   end select
end subroutine TangentialDirections

pure function ElementMass(mx, my, mz, direction, row, column, element) &
   result(value)
   implicit none
   type(Matrix1D), intent(in) :: mx, my, mz
   integer(kind=4), intent(in) :: direction, row, column, element
   real(kind=8) :: value

   select case (direction)
   case (1)
      value = mx%element_mass(row, column, element)
   case (2)
      value = my%element_mass(row, column, element)
   case default
      value = mz%element_mass(row, column, element)
   end select
end function ElementMass

!------------------------------------------------------------------------------
!> Volume forcing in the three velocity-test blocks.
!------------------------------------------------------------------------------
subroutine AssembleVolumeLoad(test_space, forcing_x, forcing_y, forcing_z, rhs)
   use Interfaces, only: forcing_fun

   implicit none
   type(Space3D), intent(in) :: test_space
   procedure(forcing_fun) :: forcing_x, forcing_y, forcing_z
   real(kind=8), intent(inout) :: rhs(:)
   integer(kind=4) :: ex, ey, ez, qx, qy, qz, ax, ay, az
   integer(kind=4) :: ix, iy, iz, scalar_index, nt
   real(kind=8) :: point(3), derivative(3), source(3), scale, basis_value

   nt = TensorSize(test_space%n)
   derivative = 0.d0
   do ez = 1, test_space%nelem(3)
      do ey = 1, test_space%nelem(2)
         do ex = 1, test_space%nelem(1)
            do qz = 1, test_space%direction(3)%ng
               do qy = 1, test_space%direction(2)%ng
                  do qx = 1, test_space%direction(1)%ng
                     point = (/test_space%direction(1)%points(qx, ex), &
                               test_space%direction(2)%points(qy, ey), &
                               test_space%direction(3)%points(qz, ez)/)
                     source(1) = forcing_x(0.d0, derivative, point)
                     source(2) = forcing_y(0.d0, derivative, point)
                     source(3) = forcing_z(0.d0, derivative, point)
                     scale = test_space%direction(1)%jacobian(ex)* &
                             test_space%direction(1)%weights(qx)* &
                             test_space%direction(2)%jacobian(ey)* &
                             test_space%direction(2)%weights(qy)* &
                             test_space%direction(3)%jacobian(ez)* &
                             test_space%direction(3)%weights(qz)
                     do az = 0, test_space%p(3)
                        iz = test_space%direction(3)%offset(ez) + az
                        do ay = 0, test_space%p(2)
                           iy = test_space%direction(2)%offset(ey) + ay
                           do ax = 0, test_space%p(1)
                              ix = test_space%direction(1)%offset(ex) + ax
                              scalar_index = LinearIndex3D(ix, iy, iz, &
                                                           test_space%n)
                              basis_value = &
                                 test_space%direction(1)%basis(0, ax, qx, ex)* &
                                 test_space%direction(2)%basis(0, ay, qy, ey)* &
                                 test_space%direction(3)%basis(0, az, qz, ez)
                              rhs(TestIndex(VELOCITY_X, scalar_index, nt) + 1) = &
                                 rhs(TestIndex(VELOCITY_X, scalar_index, nt) + 1) + &
                                 source(1)*basis_value*scale
                              rhs(TestIndex(VELOCITY_Y, scalar_index, nt) + 1) = &
                                 rhs(TestIndex(VELOCITY_Y, scalar_index, nt) + 1) + &
                                 source(2)*basis_value*scale
                              rhs(TestIndex(VELOCITY_Z, scalar_index, nt) + 1) = &
                                 rhs(TestIndex(VELOCITY_Z, scalar_index, nt) + 1) + &
                                 source(3)*basis_value*scale
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do

end subroutine AssembleVolumeLoad

!------------------------------------------------------------------------------
!> Nitsche data on all six boundary faces.
!------------------------------------------------------------------------------
subroutine AssembleBoundaryLoad(test_space, trial_degree, viscosity, &
                                penalty_factor, &
                                boundary_velocity, rhs)
   implicit none
   type(Space3D), intent(in) :: test_space
   integer(kind=4), intent(in) :: trial_degree(3)
   real(kind=8), intent(in) :: viscosity, penalty_factor
   procedure(VectorField3D) :: boundary_velocity
   real(kind=8), intent(inout) :: rhs(:)
   real(kind=8), allocatable :: average(:), average_derivative(:), jump(:)
   integer(kind=4) :: direction, tangent1, tangent2, side, face
   integer(kind=4) :: element1, element2, q1, q2, a1, a2, inormal, field
   integer(kind=4) :: coordinates(3), scalar_index, row1, row2, nt
   real(kind=8) :: normal, h, eta, point(3), velocity(3)
   real(kind=8) :: scale, tangential_basis, test_value, normal_derivative

   nt = TensorSize(test_space%n)
   do direction = 1, 3
      call TangentialDirections(direction, tangent1, tangent2)
      allocate (average(0:test_space%n(direction)))
      allocate (average_derivative(0:test_space%n(direction)))
      allocate (jump(0:test_space%n(direction)))
      eta = penalty_factor*real((trial_degree(direction) + 1)* &
                                (trial_degree(direction) + 2), kind=8)
      do side = 0, 1
         face = side*test_space%nelem(direction)
         if (side == 0) then
            normal = -1.d0
         else
            normal = 1.d0
         end if
         call ComputeTrace1D(test_space%direction(direction), face, normal, &
                             average, average_derivative, jump)
         do element2 = 1, test_space%nelem(tangent2)
            do element1 = 1, test_space%nelem(tangent1)
               h = sqrt((2.d0*test_space%direction(tangent1)% &
                         jacobian(element1))**2 + &
                        (2.d0*test_space%direction(tangent2)% &
                         jacobian(element2))**2)
               do q2 = 1, test_space%direction(tangent2)%ng
                  do q1 = 1, test_space%direction(tangent1)%ng
                     point(direction) = real(face, kind=8)/ &
                                        real(test_space%nelem(direction), kind=8)
                     point(tangent1) = test_space%direction(tangent1)% &
                                       points(q1, element1)
                     point(tangent2) = test_space%direction(tangent2)% &
                                       points(q2, element2)
                     velocity = boundary_velocity(point)
                     scale = test_space%direction(tangent1)% &
                             jacobian(element1)* &
                             test_space%direction(tangent1)%weights(q1)* &
                             test_space%direction(tangent2)% &
                             jacobian(element2)* &
                             test_space%direction(tangent2)%weights(q2)
                     do inormal = 0, test_space%n(direction)
                        if (average(inormal) == 0.d0 .and. &
                            average_derivative(inormal) == 0.d0) cycle
                        do a2 = 0, test_space%p(tangent2)
                           row2 = test_space%direction(tangent2)% &
                                  offset(element2) + a2
                           do a1 = 0, test_space%p(tangent1)
                              row1 = test_space%direction(tangent1)% &
                                     offset(element1) + a1
                              tangential_basis = &
                                 test_space%direction(tangent1)% &
                                    basis(0, a1, q1, element1)* &
                                 test_space%direction(tangent2)% &
                                    basis(0, a2, q2, element2)
                              coordinates(direction) = inormal
                              coordinates(tangent1) = row1
                              coordinates(tangent2) = row2
                              scalar_index = LinearIndex3D(coordinates(1), &
                                 coordinates(2), coordinates(3), test_space%n)
                              test_value = average(inormal)*tangential_basis
                              normal_derivative = normal* &
                                 average_derivative(inormal)*tangential_basis
                              do field = VELOCITY_X, VELOCITY_Z
                                 rhs(TestIndex(field, scalar_index, nt) + 1) = &
                                    rhs(TestIndex(field, scalar_index, nt) + 1) + &
                                    viscosity*(-normal_derivative*velocity(field) + &
                                      eta/h*velocity(field)*test_value)*scale
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      deallocate (average, average_derivative, jump)
   end do

end subroutine AssembleBoundaryLoad

!------------------------------------------------------------------------------
!> Add a symmetric mean-zero pressure constraint with one scalar multiplier.
!------------------------------------------------------------------------------
subroutine AssemblePressureGauge(trial_space, matrix)
   use sparse, only: sparse_matrix

   implicit none
   type(Space3D), intent(in) :: trial_space
   type(sparse_matrix), pointer, intent(inout) :: matrix
   real(kind=8), allocatable :: integral_x(:), integral_y(:), integral_z(:)
   integer(kind=4) :: ix, iy, iz, scalar_index, nt, nu, gauge_index, column
   real(kind=8) :: value

   nu = TensorSize(trial_space%n)
   nt = (matrix%x - 1 - 4*nu)/4
   gauge_index = matrix%x - 1
   call BasisIntegrals(trial_space%direction(1), integral_x)
   call BasisIntegrals(trial_space%direction(2), integral_y)
   call BasisIntegrals(trial_space%direction(3), integral_z)
   do iz = 0, trial_space%n(3)
      do iy = 0, trial_space%n(2)
         do ix = 0, trial_space%n(1)
            scalar_index = LinearIndex3D(ix, iy, iz, trial_space%n)
            column = TrialIndex(PRESSURE_FIELD, scalar_index, nt, nu)
            value = integral_x(ix)*integral_y(iy)*integral_z(iz)
            call AddMatrixEntry(matrix, gauge_index, column, value)
            call AddMatrixEntry(matrix, column, gauge_index, value)
         end do
      end do
   end do
   deallocate (integral_x, integral_y, integral_z)

end subroutine AssemblePressureGauge

subroutine BasisIntegrals(space, integrals)
   implicit none
   type(Space1D), intent(in) :: space
   real(kind=8), allocatable, intent(out) :: integrals(:)
   integer(kind=4) :: element, q, a, index_value

   allocate (integrals(0:space%n))
   integrals = 0.d0
   do element = 1, space%nelem
      do q = 1, space%ng
         do a = 0, space%p
            index_value = space%offset(element) + a
            integrals(index_value) = integrals(index_value) + &
               space%basis(0, a, q, element)*space%jacobian(element)* &
               space%weights(q)
         end do
      end do
   end do
end subroutine BasisIntegrals

!------------------------------------------------------------------------------
!> Initialize the replicated trial-field result and retain its knot vectors.
!------------------------------------------------------------------------------
subroutine InitializeSolution(nelem, degree, fields, solution, ierr)
   implicit none
   integer(kind=4), intent(in) :: nelem(3), degree(3)
   real(kind=8), intent(in) :: fields(:, :)
   type(IGRMStokesSolution), intent(out) :: solution
   integer(kind=4), intent(out) :: ierr
   type(Space3D) :: trial_space
   integer(kind=4) :: ng

   ng = max(2, maxval(degree) + 1)
   call BuildSpace3D(nelem, degree, 1, ng, trial_space, ierr)
   if (ierr /= IGRM_STOKES_SUCCESS) return
   solution%n = trial_space%n
   solution%p = degree
   solution%nelem = nelem
   allocate (solution%Ux(size(trial_space%direction(1)%knots)))
   allocate (solution%Uy(size(trial_space%direction(2)%knots)))
   allocate (solution%Uz(size(trial_space%direction(3)%knots)))
   solution%Ux = trial_space%direction(1)%knots
   solution%Uy = trial_space%direction(2)%knots
   solution%Uz = trial_space%direction(3)%knots
   allocate (solution%coefficients(size(fields, 1), 4))
   solution%coefficients = fields

end subroutine InitializeSolution

subroutine CleanupIGRMStokesSolution(solution)
   implicit none
   type(IGRMStokesSolution), intent(inout) :: solution

   if (allocated(solution%Ux)) deallocate (solution%Ux)
   if (allocated(solution%Uy)) deallocate (solution%Uy)
   if (allocated(solution%Uz)) deallocate (solution%Uz)
   if (allocated(solution%coefficients)) deallocate (solution%coefficients)
   solution%n = -1
   solution%p = -1
   solution%nelem = 0
end subroutine CleanupIGRMStokesSolution

!------------------------------------------------------------------------------
!> Evaluate velocity, pressure, and velocity divergence at one point.
!------------------------------------------------------------------------------
subroutine EvaluateIGRMStokes(solution, point, velocity, pressure, divergence)
   use basis, only: FindSpan, DersBasisFuns

   implicit none
   type(IGRMStokesSolution), intent(in) :: solution
   real(kind=8), intent(in) :: point(3)
   real(kind=8), intent(out) :: velocity(3), pressure, divergence
   real(kind=8), allocatable :: bx(:, :), by(:, :), bz(:, :)
   integer(kind=4) :: span_x, span_y, span_z, x0, y0, z0
   integer(kind=4) :: ax, ay, az, ix, iy, iz, index_value, field
   real(kind=8) :: basis_value, derivative_x, derivative_y, derivative_z
   real(kind=8) :: coordinate(3)

   velocity = 0.d0
   pressure = 0.d0
   divergence = 0.d0
   if (.not. allocated(solution%coefficients)) return
   coordinate = max(0.d0, min(1.d0, point))
   allocate (bx(0:solution%p(1), 0:1))
   allocate (by(0:solution%p(2), 0:1))
   allocate (bz(0:solution%p(3), 0:1))
   span_x = FindSpan(solution%n(1), solution%p(1), coordinate(1), solution%Ux)
   span_y = FindSpan(solution%n(2), solution%p(2), coordinate(2), solution%Uy)
   span_z = FindSpan(solution%n(3), solution%p(3), coordinate(3), solution%Uz)
   call DersBasisFuns(span_x, coordinate(1), solution%p(1), 1, solution%Ux, bx)
   call DersBasisFuns(span_y, coordinate(2), solution%p(2), 1, solution%Uy, by)
   call DersBasisFuns(span_z, coordinate(3), solution%p(3), 1, solution%Uz, bz)
   x0 = span_x - solution%p(1)
   y0 = span_y - solution%p(2)
   z0 = span_z - solution%p(3)

   do az = 0, solution%p(3)
      iz = z0 + az
      do ay = 0, solution%p(2)
         iy = y0 + ay
         do ax = 0, solution%p(1)
            ix = x0 + ax
            index_value = LinearIndex3D(ix, iy, iz, solution%n) + 1
            basis_value = bx(ax, 0)*by(ay, 0)*bz(az, 0)
            derivative_x = bx(ax, 1)*by(ay, 0)*bz(az, 0)
            derivative_y = bx(ax, 0)*by(ay, 1)*bz(az, 0)
            derivative_z = bx(ax, 0)*by(ay, 0)*bz(az, 1)
            do field = VELOCITY_X, VELOCITY_Z
               velocity(field) = velocity(field) + &
                  solution%coefficients(index_value, field)*basis_value
            end do
            pressure = pressure + &
               solution%coefficients(index_value, PRESSURE_FIELD)*basis_value
            divergence = divergence + &
               solution%coefficients(index_value, VELOCITY_X)*derivative_x + &
               solution%coefficients(index_value, VELOCITY_Y)*derivative_y + &
               solution%coefficients(index_value, VELOCITY_Z)*derivative_z
         end do
      end do
   end do
   deallocate (bx, by, bz)

end subroutine EvaluateIGRMStokes

!------------------------------------------------------------------------------
!> Compute global L2 velocity/pressure errors and L2 divergence by quadrature.
!------------------------------------------------------------------------------
subroutine ComputeIGRMStokesErrors(solution, exact_velocity, exact_pressure, &
                                   velocity_l2, pressure_l2, divergence_l2)
   use gauss, only: GaussRule, MAX_GAUSS_POINTS

   implicit none
   type(IGRMStokesSolution), intent(in) :: solution
   procedure(VectorField3D) :: exact_velocity
   procedure(ScalarField3D) :: exact_pressure
   real(kind=8), intent(out) :: velocity_l2, pressure_l2, divergence_l2
   integer(kind=4) :: ng, ex, ey, ez, qx, qy, qz
   real(kind=8), allocatable :: gauss_points(:), gauss_weights(:)
   real(kind=8) :: point(3), velocity(3), expected_velocity(3)
   real(kind=8) :: pressure, expected_pressure, divergence, scale

   velocity_l2 = 0.d0
   pressure_l2 = 0.d0
   divergence_l2 = 0.d0
   if (.not. allocated(solution%coefficients)) return
   ng = min(MAX_GAUSS_POINTS, max(4, maxval(solution%p) + 2))
   allocate (gauss_points(ng), gauss_weights(ng))
   call GaussRule(ng, gauss_points, gauss_weights)
   scale = 1.d0/(8.d0*real(product(solution%nelem), kind=8))

   do ez = 1, solution%nelem(3)
      do ey = 1, solution%nelem(2)
         do ex = 1, solution%nelem(1)
            do qz = 1, ng
               point(3) = (real(ez - 1, kind=8) + &
                  0.5d0*(gauss_points(qz) + 1.d0))/ &
                  real(solution%nelem(3), kind=8)
               do qy = 1, ng
                  point(2) = (real(ey - 1, kind=8) + &
                     0.5d0*(gauss_points(qy) + 1.d0))/ &
                     real(solution%nelem(2), kind=8)
                  do qx = 1, ng
                     point(1) = (real(ex - 1, kind=8) + &
                        0.5d0*(gauss_points(qx) + 1.d0))/ &
                        real(solution%nelem(1), kind=8)
                     call EvaluateIGRMStokes(solution, point, velocity, &
                                             pressure, divergence)
                     expected_velocity = exact_velocity(point)
                     expected_pressure = exact_pressure(point)
                     velocity_l2 = velocity_l2 + &
                        sum((velocity - expected_velocity)**2)* &
                        gauss_weights(qx)*gauss_weights(qy)* &
                        gauss_weights(qz)*scale
                     pressure_l2 = pressure_l2 + &
                        (pressure - expected_pressure)**2* &
                        gauss_weights(qx)*gauss_weights(qy)* &
                        gauss_weights(qz)*scale
                     divergence_l2 = divergence_l2 + divergence**2* &
                        gauss_weights(qx)*gauss_weights(qy)* &
                        gauss_weights(qz)*scale
                  end do
               end do
            end do
         end do
      end do
   end do
   velocity_l2 = sqrt(max(0.d0, velocity_l2))
   pressure_l2 = sqrt(max(0.d0, pressure_l2))
   divergence_l2 = sqrt(max(0.d0, divergence_l2))
   deallocate (gauss_points, gauss_weights)

end subroutine ComputeIGRMStokesErrors

!------------------------------------------------------------------------------
!> Write one VTI containing a 3-component Velocity and scalar Pressure array.
!> The routine is collective; only rank zero performs file I/O.
!------------------------------------------------------------------------------
subroutine WriteIGRMStokesVTI(solution, path, ierr)
   use parallelism, only: MYRANK
   use mpi

   implicit none
   type(IGRMStokesSolution), intent(in) :: solution
   character(len=*), intent(in) :: path
   integer(kind=4), intent(out) :: ierr
   integer(kind=4), parameter :: RESOLUTION = 50
   integer(kind=4) :: unit_number, io_status, mpi_ierr, ix, iy, iz
   real(kind=8) :: point(3), velocity(3), pressure, divergence, spacing
   logical :: opened

   ierr = IGRM_STOKES_SUCCESS
   if (MYRANK == 0) then
      unit_number = -1
      do ix = 20, 99
         inquire (unit=ix, opened=opened)
         if (.not. opened) then
            unit_number = ix
            exit
         end if
      end do
      if (unit_number < 0) then
         ierr = IGRM_STOKES_IO_ERROR
      else
         open (unit=unit_number, file=path, status='replace', action='write', &
               form='formatted', iostat=io_status)
         if (io_status /= 0) then
            ierr = IGRM_STOKES_IO_ERROR
         else
            spacing = 1.d0/real(RESOLUTION, kind=8)
            write (unit_number, '(A)', iostat=io_status) &
               '<?xml version="1.0"?>'
            if (io_status /= 0) go to 900
            write (unit_number, '(A)', iostat=io_status) &
               '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'
            if (io_status /= 0) go to 900
            write (unit_number, '(A,I0,A,I0,A,I0,A,3(ES16.8,1X),A)', &
                   iostat=io_status) &
               '  <ImageData WholeExtent="0 ', RESOLUTION, ' 0 ', RESOLUTION, &
               ' 0 ', RESOLUTION, '" Origin="0 0 0" Spacing="', &
               spacing, spacing, spacing, '">'
            if (io_status /= 0) go to 900
            write (unit_number, '(A,I0,A,I0,A,I0,A)', iostat=io_status) &
               '    <Piece Extent="0 ', RESOLUTION, ' 0 ', RESOLUTION, &
               ' 0 ', RESOLUTION, '">'
            if (io_status /= 0) go to 900
            write (unit_number, '(A)', iostat=io_status) &
               '      <PointData Scalars="Pressure" Vectors="Velocity">'
            if (io_status /= 0) go to 900
            write (unit_number, '(A)', iostat=io_status) &
               '        <DataArray Name="Velocity" type="Float64" ' // &
               'NumberOfComponents="3" format="ascii">'
            if (io_status /= 0) go to 900
            do iz = 0, RESOLUTION
               point(3) = real(iz, kind=8)/real(RESOLUTION, kind=8)
               do iy = 0, RESOLUTION
                  point(2) = real(iy, kind=8)/real(RESOLUTION, kind=8)
                  do ix = 0, RESOLUTION
                     point(1) = real(ix, kind=8)/real(RESOLUTION, kind=8)
                     call EvaluateIGRMStokes(solution, point, velocity, &
                                             pressure, divergence)
                     write (unit_number, '(3(ES24.16,1X))', iostat=io_status) &
                        velocity
                     if (io_status /= 0) go to 900
                  end do
               end do
            end do
            write (unit_number, '(A)', iostat=io_status) '        </DataArray>'
            if (io_status /= 0) go to 900
            write (unit_number, '(A)', iostat=io_status) &
               '        <DataArray Name="Pressure" type="Float64" ' // &
               'NumberOfComponents="1" format="ascii">'
            if (io_status /= 0) go to 900
            do iz = 0, RESOLUTION
               point(3) = real(iz, kind=8)/real(RESOLUTION, kind=8)
               do iy = 0, RESOLUTION
                  point(2) = real(iy, kind=8)/real(RESOLUTION, kind=8)
                  do ix = 0, RESOLUTION
                     point(1) = real(ix, kind=8)/real(RESOLUTION, kind=8)
                     call EvaluateIGRMStokes(solution, point, velocity, &
                                             pressure, divergence)
                     write (unit_number, '(ES24.16)', iostat=io_status) pressure
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
900         continue
            if (io_status /= 0) ierr = IGRM_STOKES_IO_ERROR
            close (unit_number, iostat=io_status)
            if (io_status /= 0) ierr = IGRM_STOKES_IO_ERROR
         end if
      end if
   end if

   call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
   if (mpi_ierr /= MPI_SUCCESS) ierr = IGRM_STOKES_MPI_ERROR

end subroutine WriteIGRMStokesVTI

!------------------------------------------------------------------------------
!> Sparse insertion helpers. Repeated contributions are accumulated by sparse.
!------------------------------------------------------------------------------
subroutine AddMatrixEntry(matrix, row, column, value)
   use sparse, only: sparse_matrix, add

   implicit none
   type(sparse_matrix), pointer, intent(inout) :: matrix
   integer(kind=4), intent(in) :: row, column
   real(kind=8), intent(in) :: value

   if (value == 0.d0) return
   call add(matrix, row, column, value)
end subroutine AddMatrixEntry

subroutine AddSymmetricB(matrix, row, column, value)
   use sparse, only: sparse_matrix

   implicit none
   type(sparse_matrix), pointer, intent(inout) :: matrix
   integer(kind=4), intent(in) :: row, column
   real(kind=8), intent(in) :: value

   call AddMatrixEntry(matrix, row, column, value)
   call AddMatrixEntry(matrix, column, row, value)
end subroutine AddSymmetricB

subroutine ComputeResidual(matrix, rhs, solution, residual_rms, &
                           relative_residual)
   use sparse, only: sparse_matrix

   implicit none
   type(sparse_matrix), pointer, intent(in) :: matrix
   real(kind=8), intent(in) :: rhs(:), solution(:)
   real(kind=8), intent(out) :: residual_rms, relative_residual
   real(kind=8), allocatable :: product_vector(:)
   real(kind=8) :: residual_norm, rhs_norm
   integer(kind=4) :: row, entry

   allocate (product_vector(size(rhs)))
   product_vector = 0.d0
   do row = 0, matrix%x - 1
      do entry = 1, matrix%rows(row)%nnz
         product_vector(row + 1) = product_vector(row + 1) + &
            matrix%rows(row)%val(entry)* &
            solution(matrix%rows(row)%col(entry) + 1)
      end do
   end do
   residual_norm = sqrt(sum((rhs - product_vector)**2))
   rhs_norm = sqrt(sum(rhs**2))
   residual_rms = residual_norm/sqrt(real(size(rhs), kind=8))
   relative_residual = residual_norm/max(rhs_norm, tiny(0.d0))
   deallocate (product_vector)
end subroutine ComputeResidual

pure integer(kind=4) function TensorSize(n) result(value)
   implicit none
   integer(kind=4), intent(in) :: n(3)
   value = (n(1) + 1)*(n(2) + 1)*(n(3) + 1)
end function TensorSize

pure integer(kind=4) function LinearIndex3D(ix, iy, iz, n) result(value)
   implicit none
   integer(kind=4), intent(in) :: ix, iy, iz, n(3)
   value = ix + (n(1) + 1)*(iy + (n(2) + 1)*iz)
end function LinearIndex3D

pure integer(kind=4) function TestIndex(field, scalar_index, nt) result(value)
   implicit none
   integer(kind=4), intent(in) :: field, scalar_index, nt
   value = (field - 1)*nt + scalar_index
end function TestIndex

pure integer(kind=4) function TrialIndex(field, scalar_index, nt, nu) &
   result(value)
   implicit none
   integer(kind=4), intent(in) :: field, scalar_index, nt, nu
   value = 4*nt + (field - 1)*nu + scalar_index
end function TrialIndex

end module igrm_stokes_solver
