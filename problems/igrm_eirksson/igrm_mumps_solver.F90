!------------------------------------------------------------------------------
!> @file igrm_mumps_solver.F90
!> @brief Full-tensor direct iGRM solver for the three-dimensional Eriksson case.
!------------------------------------------------------------------------------
module igrm_mumps_solver

   implicit none

   private

   integer(kind=4), parameter, public :: IGRM_MUMPS_SUCCESS = 0
   integer(kind=4), parameter, public :: IGRM_MUMPS_INVALID_CONFIGURATION = -7301
   integer(kind=4), parameter, public :: IGRM_MUMPS_SYSTEM_TOO_LARGE = -7302
   integer(kind=4), parameter, public :: IGRM_MUMPS_RESIDUAL_TOO_LARGE = -7303

   type, public :: IGRMMumpsStats
      integer(kind=4) :: system_size = 0
      integer(kind=8) :: nonzeros = 0_8
      real(kind=8) :: residual_rms = huge(0.d0)
      real(kind=8) :: relative_residual = huge(0.d0)
   end type IGRMMumpsStats

   type :: Matrix1D
      real(kind=8), allocatable :: mass(:, :)
      real(kind=8), allocatable :: stiffness(:, :)
      real(kind=8), allocatable :: advection(:, :)
   end type Matrix1D

   type :: Adjacency1D
      integer(kind=4), allocatable :: count(:)
      integer(kind=4), allocatable :: column(:, :)
   end type Adjacency1D

   ! A compact duplicate of the assembled contributions is retained only
   ! until the solve finishes.  The library sparse-row implementation is
   ! intentionally opaque; this problem-local representation lets us verify
   ! the direct solution without changing the core library API.
   type :: MatrixEntries
      integer(kind=4) :: count = 0
      integer(kind=4) :: capacity = 0
      integer(kind=4), allocatable :: row(:)
      integer(kind=4), allocatable :: column(:)
      real(kind=8), allocatable :: value(:)
   end type MatrixEntries

   public :: AssembleIGRMMumps3DSystem
   public :: SolveIGRMMumps3D

contains

!------------------------------------------------------------------------------
!> Assemble and solve [G,-B;B^T,0][r,u]^T=[-l,0]^T once with MUMPS.
!>
!> The solve is performed on rank zero through the existing sequential
!> `MPI_COMM_SELF` MUMPS wrapper.  Trial coefficients are then broadcast and
!> sliced into the canonical ADS ownership layout on every rank.
!------------------------------------------------------------------------------
subroutine SolveIGRMMumps3D(ads_test, ads_trial, epsilon, beta, &
                            residual_tolerance, forcing, local_solution, &
                            stats, ierr)
   use Interfaces, only: forcing_fun
   use Setup, only: ADS_Setup
   use mumps_solver, only: SolveOneDirection
   use parallelism, only: MYRANK
   use sparse, only: sparse_matrix, clear_matrix
   use, intrinsic :: iso_fortran_env, only: error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mpi

   implicit none

   type(ADS_Setup), intent(in) :: ads_test, ads_trial
   real(kind=8), intent(in) :: epsilon
   real(kind=8), intent(in), dimension(3) :: beta
   real(kind=8), intent(in) :: residual_tolerance
   procedure(forcing_fun) :: forcing
   real(kind=8), allocatable, intent(out) :: local_solution(:, :)
   type(IGRMMumpsStats), intent(out) :: stats
   integer(kind=4), intent(out) :: ierr

   type(sparse_matrix), pointer :: system_matrix
   type(MatrixEntries) :: entries
   real(kind=8), allocatable :: right_hand_side(:)
   real(kind=8), allocatable :: solve_buffer(:, :)
   real(kind=8), allocatable :: global_solution(:)
   integer(kind=8) :: test_size_64, trial_size_64, system_size_64
   integer(kind=4) :: test_size, trial_size, system_size
   integer(kind=4) :: solver_ierr, mpi_ierr

   nullify (system_matrix)
   stats = IGRMMumpsStats()
   ierr = IGRM_MUMPS_SUCCESS

   call ValidateConfiguration(ads_test, ads_trial, epsilon, beta, &
                              residual_tolerance, test_size_64, &
                              trial_size_64, ierr)
   if (ierr /= IGRM_MUMPS_SUCCESS) return

   if (test_size_64 > huge(0_8) - trial_size_64) then
      ierr = IGRM_MUMPS_SYSTEM_TOO_LARGE
      return
   end if
   system_size_64 = test_size_64 + trial_size_64
   if (system_size_64 > int(huge(0_4), kind=8)) then
      ierr = IGRM_MUMPS_SYSTEM_TOO_LARGE
      return
   end if

   test_size = int(test_size_64, kind=4)
   trial_size = int(trial_size_64, kind=4)
   system_size = int(system_size_64, kind=4)
   stats%system_size = system_size

   allocate (local_solution(ads_trial%s(1), ads_trial%s(2)*ads_trial%s(3)))
   allocate (global_solution(trial_size))
   local_solution = 0.d0
   global_solution = 0.d0

   if (MYRANK == 0) then
      call AssembleTrackedSystem(ads_test, ads_trial, epsilon, beta, forcing, &
                                 system_matrix, right_hand_side, entries, ierr)

      if (ierr == IGRM_MUMPS_SUCCESS) then
         stats%nonzeros = system_matrix%total_entries
         allocate (solve_buffer(system_size, 1))
         solve_buffer(:, 1) = right_hand_side

         call SolveOneDirection(solve_buffer, 1, system_size - 1, 0, &
                                system_matrix, solver_ierr)
         if (solver_ierr /= 0) then
            ierr = solver_ierr
         else
            call ComputeResidual(entries, right_hand_side, solve_buffer(:, 1), &
                                 stats%residual_rms, stats%relative_residual)
            if (.not. ieee_is_finite(stats%residual_rms) .or. &
                .not. ieee_is_finite(stats%relative_residual) .or. &
                stats%relative_residual > residual_tolerance) then
               write (error_unit, '(A,1X,ES24.16,1X,A,1X,ES24.16)') &
                  'iGRM-MUMPS relative residual', stats%relative_residual, &
                  'exceeds tolerance', residual_tolerance
               flush (error_unit)
               ierr = IGRM_MUMPS_RESIDUAL_TOO_LARGE
            else
               global_solution = solve_buffer(test_size + 1:system_size, 1)
            end if
         end if
      end if

      if (associated(system_matrix)) call clear_matrix(system_matrix)
      if (allocated(right_hand_side)) deallocate (right_hand_side)
      if (allocated(solve_buffer)) deallocate (solve_buffer)
      call ReleaseEntries(entries)
   end if

   call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
   if (ierr /= IGRM_MUMPS_SUCCESS) then
      deallocate (local_solution, global_solution)
      return
   end if

   call MPI_Bcast(stats%system_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
   call MPI_Bcast(stats%nonzeros, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpi_ierr)
   call MPI_Bcast(stats%residual_rms, 1, MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, mpi_ierr)
   call MPI_Bcast(stats%relative_residual, 1, MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, mpi_ierr)
   call MPI_Bcast(global_solution, trial_size, MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, mpi_ierr)

   call CopyOwnedCoefficients(ads_trial, global_solution, local_solution)
   deallocate (global_solution)

end subroutine SolveIGRMMumps3D

!------------------------------------------------------------------------------
!> Public assembly entry used by focused matrix-structure tests.
!------------------------------------------------------------------------------
subroutine AssembleIGRMMumps3DSystem(ads_test, ads_trial, epsilon, beta, &
                                     forcing, system_matrix, &
                                     right_hand_side, ierr)
   use Interfaces, only: forcing_fun
   use Setup, only: ADS_Setup
   use sparse, only: sparse_matrix

   implicit none

   type(ADS_Setup), intent(in) :: ads_test, ads_trial
   real(kind=8), intent(in) :: epsilon
   real(kind=8), intent(in), dimension(3) :: beta
   procedure(forcing_fun) :: forcing
   type(sparse_matrix), pointer, intent(out) :: system_matrix
   real(kind=8), allocatable, intent(out) :: right_hand_side(:)
   integer(kind=4), intent(out) :: ierr

   type(MatrixEntries) :: entries
   integer(kind=8) :: test_size, trial_size, system_size

   nullify (system_matrix)
   call ValidateConfiguration(ads_test, ads_trial, epsilon, beta, 1.d0, &
                              test_size, trial_size, ierr)
   if (test_size > huge(0_8) - trial_size) then
      ierr = IGRM_MUMPS_SYSTEM_TOO_LARGE
      call ReleaseEntries(entries)
      return
   end if
   system_size = test_size + trial_size
   if (ierr == IGRM_MUMPS_SUCCESS) then
      if (system_size > int(huge(0_4), kind=8)) then
         ierr = IGRM_MUMPS_SYSTEM_TOO_LARGE
      end if
   end if
   if (ierr == IGRM_MUMPS_SUCCESS) then
      call AssembleTrackedSystem(ads_test, ads_trial, epsilon, beta, forcing, &
                                 system_matrix, right_hand_side, entries, ierr)
   end if
   call ReleaseEntries(entries)

end subroutine AssembleIGRMMumps3DSystem

!------------------------------------------------------------------------------
!> Assemble the full three-dimensional iGRM saddle matrix and load vector.
!------------------------------------------------------------------------------
subroutine AssembleTrackedSystem(ads_test, ads_trial, epsilon, beta, forcing, &
                                 system_matrix, right_hand_side, entries, ierr)
   use Interfaces, only: forcing_fun
   use Setup, only: ADS_Setup
   use sparse, only: sparse_matrix, initialize_sparse

   implicit none

   type(ADS_Setup), intent(in) :: ads_test, ads_trial
   real(kind=8), intent(in) :: epsilon
   real(kind=8), intent(in), dimension(3) :: beta
   procedure(forcing_fun) :: forcing
   type(sparse_matrix), pointer, intent(out) :: system_matrix
   real(kind=8), allocatable, intent(out) :: right_hand_side(:)
   type(MatrixEntries), intent(out) :: entries
   integer(kind=4), intent(out) :: ierr

   type(Matrix1D) :: test_x, test_y, test_z
   type(Matrix1D) :: mixed_x, mixed_y, mixed_z
   type(Adjacency1D) :: test_adj_x, test_adj_y, test_adj_z
   type(Adjacency1D) :: mixed_adj_x, mixed_adj_y, mixed_adj_z
   type(Adjacency1D) :: trans_adj_x, trans_adj_y, trans_adj_z
   integer(kind=4) :: test_size, trial_size, system_size
   integer(kind=4) :: ix, iy, iz, jx, jy, jz
   integer(kind=4) :: ax, ay, az, row, column
   real(kind=8) :: h, h2, value

   nullify (system_matrix)
   entries%count = 0
   entries%capacity = 0
   ierr = IGRM_MUMPS_SUCCESS

   test_size = TensorSize(ads_test%n)
   trial_size = TensorSize(ads_trial%n)
   system_size = test_size + trial_size

   call AssembleOneSpace1D(ads_test%n(1), ads_test%p(1), &
                           ads_test%nelem(1), ads_test%ng(1), ads_test%Ox, &
                           ads_test%Jx, ads_test%Wx, ads_test%NNx, test_x)
   call AssembleOneSpace1D(ads_test%n(2), ads_test%p(2), &
                           ads_test%nelem(2), ads_test%ng(2), ads_test%Oy, &
                           ads_test%Jy, ads_test%Wy, ads_test%NNy, test_y)
   call AssembleOneSpace1D(ads_test%n(3), ads_test%p(3), &
                           ads_test%nelem(3), ads_test%ng(3), ads_test%Oz, &
                           ads_test%Jz, ads_test%Wz, ads_test%NNz, test_z)

   call AssembleMixed1D(ads_test%n(1), ads_test%p(1), ads_trial%n(1), &
                        ads_trial%p(1), ads_test%nelem(1), ads_test%ng(1), &
                        ads_test%Ox, ads_trial%Ox, ads_test%Jx, ads_test%Wx, &
                        ads_test%NNx, ads_trial%NNx, mixed_x)
   call AssembleMixed1D(ads_test%n(2), ads_test%p(2), ads_trial%n(2), &
                        ads_trial%p(2), ads_test%nelem(2), ads_test%ng(2), &
                        ads_test%Oy, ads_trial%Oy, ads_test%Jy, ads_test%Wy, &
                        ads_test%NNy, ads_trial%NNy, mixed_y)
   call AssembleMixed1D(ads_test%n(3), ads_test%p(3), ads_trial%n(3), &
                        ads_trial%p(3), ads_test%nelem(3), ads_test%ng(3), &
                        ads_test%Oz, ads_trial%Oz, ads_test%Jz, ads_test%Wz, &
                        ads_test%NNz, ads_trial%NNz, mixed_z)

   call BuildAdjacency(test_x%mass, .false., test_adj_x)
   call BuildAdjacency(test_y%mass, .false., test_adj_y)
   call BuildAdjacency(test_z%mass, .false., test_adj_z)
   call BuildAdjacency(mixed_x%mass, .false., mixed_adj_x)
   call BuildAdjacency(mixed_y%mass, .false., mixed_adj_y)
   call BuildAdjacency(mixed_z%mass, .false., mixed_adj_z)
   call BuildAdjacency(mixed_x%mass, .true., trans_adj_x)
   call BuildAdjacency(mixed_y%mass, .true., trans_adj_y)
   call BuildAdjacency(mixed_z%mass, .true., trans_adj_z)

   h = (2.d0*maxval(ads_trial%Jx) * 2.d0*maxval(ads_trial%Jy) * &
        2.d0*maxval(ads_trial%Jz))**(1.d0/3.d0)
   h2 = h*h

   call initialize_sparse(system_size, system_size, system_matrix)
   call InitializeEntries(entries)
   allocate (right_hand_side(system_size))
   right_hand_side = 0.d0

   ! Test rows: the H1-like Gram block G and the negative PDE block -B.
   do iz = 0, ads_test%n(3)
      do iy = 0, ads_test%n(2)
         do ix = 0, ads_test%n(1)
            row = LinearIndex3D(ix, iy, iz, ads_test%n)
            if (IsBoundary3D(ix, iy, iz, ads_test%n)) then
               call AddTracked(system_matrix, entries, row, row, 1.d0)
               cycle
            end if

            do az = 1, test_adj_z%count(iz)
               jz = test_adj_z%column(az, iz)
               do ay = 1, test_adj_y%count(iy)
                  jy = test_adj_y%column(ay, iy)
                  do ax = 1, test_adj_x%count(ix)
                     jx = test_adj_x%column(ax, ix)
                     if (IsBoundary3D(jx, jy, jz, ads_test%n)) cycle
                     column = LinearIndex3D(jx, jy, jz, ads_test%n)
                     value = GramTensorValue(test_x, test_y, test_z, &
                                             ix, iy, iz, jx, jy, jz, h2)
                     call AddTracked(system_matrix, entries, row, column, value)
                  end do
               end do
            end do

            do az = 1, mixed_adj_z%count(iz)
               jz = mixed_adj_z%column(az, iz)
               do ay = 1, mixed_adj_y%count(iy)
                  jy = mixed_adj_y%column(ay, iy)
                  do ax = 1, mixed_adj_x%count(ix)
                     jx = mixed_adj_x%column(ax, ix)
                     column = test_size + &
                              LinearIndex3D(jx, jy, jz, ads_trial%n)
                     value = OperatorTensorValue(mixed_x, mixed_y, mixed_z, &
                                                 ix, iy, iz, jx, jy, jz, &
                                                 epsilon, beta)
                     call AddTracked(system_matrix, entries, row, column, -value)
                  end do
               end do
            end do
         end do
      end do
   end do

   ! Trial rows: B^T for internal trial DOFs and homogeneous Dirichlet rows.
   do iz = 0, ads_trial%n(3)
      do iy = 0, ads_trial%n(2)
         do ix = 0, ads_trial%n(1)
            row = test_size + LinearIndex3D(ix, iy, iz, ads_trial%n)
            if (IsBoundary3D(ix, iy, iz, ads_trial%n)) then
               call AddTracked(system_matrix, entries, row, row, 1.d0)
               cycle
            end if

            do az = 1, trans_adj_z%count(iz)
               jz = trans_adj_z%column(az, iz)
               do ay = 1, trans_adj_y%count(iy)
                  jy = trans_adj_y%column(ay, iy)
                  do ax = 1, trans_adj_x%count(ix)
                     jx = trans_adj_x%column(ax, ix)
                     if (IsBoundary3D(jx, jy, jz, ads_test%n)) cycle
                     column = LinearIndex3D(jx, jy, jz, ads_test%n)
                     value = OperatorTensorValue(mixed_x, mixed_y, mixed_z, &
                                                 jx, jy, jz, ix, iy, iz, &
                                                 epsilon, beta)
                     call AddTracked(system_matrix, entries, row, column, value)
                  end do
               end do
            end do
         end do
      end do
   end do

   call AssembleLoadVector(ads_test, forcing, right_hand_side(1:test_size))

end subroutine AssembleTrackedSystem

!------------------------------------------------------------------------------
!> Validate the direct-solver contract and compute tensor sizes safely.
!------------------------------------------------------------------------------
subroutine ValidateConfiguration(ads_test, ads_trial, epsilon, beta, &
                                 residual_tolerance, test_size, trial_size, ierr)
   use Setup, only: ADS_Setup
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none

   type(ADS_Setup), intent(in) :: ads_test, ads_trial
   real(kind=8), intent(in) :: epsilon, beta(3), residual_tolerance
   integer(kind=8), intent(out) :: test_size, trial_size
   integer(kind=4), intent(out) :: ierr
   logical :: test_size_valid, trial_size_valid

   ierr = IGRM_MUMPS_SUCCESS
   if (any(ads_test%n < 0) .or. any(ads_trial%n < 0)) then
      test_size = 0_8
      trial_size = 0_8
      ierr = IGRM_MUMPS_INVALID_CONFIGURATION
      return
   end if
   call SafeTensorSize(ads_test%n, test_size, test_size_valid)
   call SafeTensorSize(ads_trial%n, trial_size, trial_size_valid)

   if (.not. test_size_valid .or. .not. trial_size_valid) then
      ierr = IGRM_MUMPS_SYSTEM_TOO_LARGE
   else if (.not. ieee_is_finite(epsilon) .or. &
            .not. ieee_is_finite(residual_tolerance) .or. &
            .not. all(ieee_is_finite(beta))) then
      ierr = IGRM_MUMPS_INVALID_CONFIGURATION
   else if (epsilon <= 0.d0 .or. residual_tolerance <= 0.d0) then
      ierr = IGRM_MUMPS_INVALID_CONFIGURATION
   else if (any(ads_test%nelem <= 0) .or. &
            any(ads_trial%nelem <= 0) .or. &
            any(ads_test%ng <= 0) .or. any(ads_trial%ng <= 0)) then
      ierr = IGRM_MUMPS_INVALID_CONFIGURATION
   else if (any(ads_test%nelem /= ads_trial%nelem)) then
      ierr = IGRM_MUMPS_INVALID_CONFIGURATION
   else if (any(ads_test%ng /= ads_trial%ng)) then
      ierr = IGRM_MUMPS_INVALID_CONFIGURATION
   else if (any(ads_trial%p < 1) .or. &
            any(ads_test%p <= ads_trial%p)) then
      ierr = IGRM_MUMPS_INVALID_CONFIGURATION
   else if (any(ads_test%n < ads_test%p) .or. &
            any(ads_trial%n < ads_trial%p)) then
      ierr = IGRM_MUMPS_INVALID_CONFIGURATION
   else if (any(ads_test%n - ads_test%p + 1 /= ads_test%nelem) .or. &
            any(ads_trial%n - ads_trial%p + 1 /= ads_trial%nelem)) then
      ierr = IGRM_MUMPS_INVALID_CONFIGURATION
   else if (test_size <= 0_8 .or. trial_size <= 0_8) then
      ierr = IGRM_MUMPS_INVALID_CONFIGURATION
   end if

end subroutine ValidateConfiguration

!------------------------------------------------------------------------------
!> Compute (n_x+1)(n_y+1)(n_z+1) without signed-integer overflow.
!------------------------------------------------------------------------------
subroutine SafeTensorSize(n, size_value, valid)
   implicit none

   integer(kind=4), intent(in), dimension(3) :: n
   integer(kind=8), intent(out) :: size_value
   logical, intent(out) :: valid
   integer(kind=8) :: extent
   integer(kind=4) :: direction

   size_value = 1_8
   valid = .true.
   do direction = 1, 3
      extent = int(n(direction), kind=8) + 1_8
      if (extent <= 0_8 .or. size_value > huge(size_value)/extent) then
         size_value = 0_8
         valid = .false.
         return
      end if
      size_value = size_value*extent
   end do

end subroutine SafeTensorSize

!------------------------------------------------------------------------------
!> Assemble one-dimensional test-space mass and stiffness matrices.
!------------------------------------------------------------------------------
subroutine AssembleOneSpace1D(n, p, nelem, ng, offsets, jacobian, weights, &
                              basis_values, matrix)
   implicit none

   integer(kind=4), intent(in) :: n, p, nelem, ng
   integer(kind=4), intent(in) :: offsets(:)
   real(kind=8), intent(in) :: jacobian(:), weights(:)
   real(kind=8), intent(in) :: basis_values(0:, 0:, :, :)
   type(Matrix1D), intent(out) :: matrix
   integer(kind=4) :: element, q, a, b, ia, ib
   real(kind=8) :: scale

   allocate (matrix%mass(0:n, 0:n), matrix%stiffness(0:n, 0:n))
   allocate (matrix%advection(0:n, 0:n))
   matrix%mass = 0.d0
   matrix%stiffness = 0.d0
   matrix%advection = 0.d0

   do element = 1, nelem
      do q = 1, ng
         scale = jacobian(element)*weights(q)
         do a = 0, p
            ia = offsets(element) + a
            do b = 0, p
               ib = offsets(element) + b
               matrix%mass(ia, ib) = matrix%mass(ia, ib) + &
                  basis_values(0, a, q, element)* &
                  basis_values(0, b, q, element)*scale
               matrix%stiffness(ia, ib) = matrix%stiffness(ia, ib) + &
                  basis_values(1, a, q, element)* &
                  basis_values(1, b, q, element)*scale
            end do
         end do
      end do
   end do

end subroutine AssembleOneSpace1D

!------------------------------------------------------------------------------
!> Assemble test-row/trial-column mass, stiffness and advection matrices.
!------------------------------------------------------------------------------
subroutine AssembleMixed1D(n_test, p_test, n_trial, p_trial, nelem, ng, &
                           offsets_test, offsets_trial, jacobian, weights, &
                           test_values, trial_values, matrix)
   implicit none

   integer(kind=4), intent(in) :: n_test, p_test, n_trial, p_trial, nelem, ng
   integer(kind=4), intent(in) :: offsets_test(:), offsets_trial(:)
   real(kind=8), intent(in) :: jacobian(:), weights(:)
   real(kind=8), intent(in) :: test_values(0:, 0:, :, :)
   real(kind=8), intent(in) :: trial_values(0:, 0:, :, :)
   type(Matrix1D), intent(out) :: matrix
   integer(kind=4) :: element, q, a, b, ia, ib
   real(kind=8) :: scale

   allocate (matrix%mass(0:n_test, 0:n_trial))
   allocate (matrix%stiffness(0:n_test, 0:n_trial))
   allocate (matrix%advection(0:n_test, 0:n_trial))
   matrix%mass = 0.d0
   matrix%stiffness = 0.d0
   matrix%advection = 0.d0

   do element = 1, nelem
      do q = 1, ng
         scale = jacobian(element)*weights(q)
         do a = 0, p_test
            ia = offsets_test(element) + a
            do b = 0, p_trial
               ib = offsets_trial(element) + b
               matrix%mass(ia, ib) = matrix%mass(ia, ib) + &
                  test_values(0, a, q, element)* &
                  trial_values(0, b, q, element)*scale
               matrix%stiffness(ia, ib) = matrix%stiffness(ia, ib) + &
                  test_values(1, a, q, element)* &
                  trial_values(1, b, q, element)*scale
               matrix%advection(ia, ib) = matrix%advection(ia, ib) + &
                  test_values(0, a, q, element)* &
                  trial_values(1, b, q, element)*scale
            end do
         end do
      end do
   end do

end subroutine AssembleMixed1D

!------------------------------------------------------------------------------
!> Build overlap lists from a one-dimensional mixed mass matrix.
!------------------------------------------------------------------------------
subroutine BuildAdjacency(mass, transposed, adjacency)
   implicit none

   real(kind=8), intent(in) :: mass(0:, 0:)
   logical, intent(in) :: transposed
   type(Adjacency1D), intent(out) :: adjacency
   integer(kind=4) :: rows, columns, row, column, entry_count

   if (transposed) then
      rows = ubound(mass, 2) + 1
      columns = ubound(mass, 1) + 1
   else
      rows = ubound(mass, 1) + 1
      columns = ubound(mass, 2) + 1
   end if

   allocate (adjacency%count(0:rows - 1))
   allocate (adjacency%column(1:max(1, columns), 0:rows - 1))
   adjacency%count = 0
   adjacency%column = -1

   do row = 0, rows - 1
      entry_count = 0
      do column = 0, columns - 1
         if (transposed) then
            if (mass(column, row) <= 0.d0) cycle
         else
            if (mass(row, column) <= 0.d0) cycle
         end if
         entry_count = entry_count + 1
         adjacency%column(entry_count, row) = column
      end do
      adjacency%count(row) = entry_count
   end do

end subroutine BuildAdjacency

!------------------------------------------------------------------------------
!> Tensor value of G=(v,w)+h^2(grad v,grad w).
!------------------------------------------------------------------------------
pure function GramTensorValue(mx, my, mz, ix, iy, iz, jx, jy, jz, h2) &
   result(value)
   implicit none

   type(Matrix1D), intent(in) :: mx, my, mz
   integer(kind=4), intent(in) :: ix, iy, iz, jx, jy, jz
   real(kind=8), intent(in) :: h2
   real(kind=8) :: value

   value = mx%mass(ix, jx)*my%mass(iy, jy)*mz%mass(iz, jz)
   value = value + h2*( &
      mx%stiffness(ix, jx)*my%mass(iy, jy)*mz%mass(iz, jz) + &
      mx%mass(ix, jx)*my%stiffness(iy, jy)*mz%mass(iz, jz) + &
      mx%mass(ix, jx)*my%mass(iy, jy)*mz%stiffness(iz, jz))

end function GramTensorValue

!------------------------------------------------------------------------------
!> Tensor value of B(v,u)=eps(grad v,grad u)+(v,beta.grad u).
!------------------------------------------------------------------------------
pure function OperatorTensorValue(mx, my, mz, ix, iy, iz, jx, jy, jz, &
                                  epsilon, beta) result(value)
   implicit none

   type(Matrix1D), intent(in) :: mx, my, mz
   integer(kind=4), intent(in) :: ix, iy, iz, jx, jy, jz
   real(kind=8), intent(in) :: epsilon
   real(kind=8), intent(in), dimension(3) :: beta
   real(kind=8) :: value

   value = epsilon*( &
      mx%stiffness(ix, jx)*my%mass(iy, jy)*mz%mass(iz, jz) + &
      mx%mass(ix, jx)*my%stiffness(iy, jy)*mz%mass(iz, jz) + &
      mx%mass(ix, jx)*my%mass(iy, jy)*mz%stiffness(iz, jz))
   value = value + &
      beta(1)*mx%advection(ix, jx)*my%mass(iy, jy)*mz%mass(iz, jz)
   value = value + &
      beta(2)*mx%mass(ix, jx)*my%advection(iy, jy)*mz%mass(iz, jz)
   value = value + &
      beta(3)*mx%mass(ix, jx)*my%mass(iy, jy)*mz%advection(iz, jz)

end function OperatorTensorValue

!------------------------------------------------------------------------------
!> Assemble the signed load block [-l,0].
!------------------------------------------------------------------------------
subroutine AssembleLoadVector(ads, forcing, load)
   use Interfaces, only: forcing_fun
   use Setup, only: ADS_Setup

   implicit none

   type(ADS_Setup), intent(in) :: ads
   procedure(forcing_fun) :: forcing
   real(kind=8), intent(inout) :: load(:)
   integer(kind=4) :: ex, ey, ez, qx, qy, qz, ax, ay, az
   integer(kind=4) :: ix, iy, iz, row
   real(kind=8), dimension(3) :: point, zero_derivative
   real(kind=8) :: scale, source, basis_value

   zero_derivative = 0.d0

   do ez = 1, ads%nelem(3)
      do ey = 1, ads%nelem(2)
         do ex = 1, ads%nelem(1)
            do qz = 1, ads%ng(3)
               do qy = 1, ads%ng(2)
                  do qx = 1, ads%ng(1)
                     point = (/ads%Xx(qx, ex), ads%Xy(qy, ey), &
                               ads%Xz(qz, ez)/)
                     source = forcing(0.d0, zero_derivative, point)
                     scale = ads%Jx(ex)*ads%Jy(ey)*ads%Jz(ez) * &
                             ads%Wx(qx)*ads%Wy(qy)*ads%Wz(qz)
                     do az = 0, ads%p(3)
                        iz = ads%Oz(ez) + az
                        do ay = 0, ads%p(2)
                           iy = ads%Oy(ey) + ay
                           do ax = 0, ads%p(1)
                              ix = ads%Ox(ex) + ax
                              row = LinearIndex3D(ix, iy, iz, ads%n)
                              basis_value = ads%NNx(0, ax, qx, ex) * &
                                            ads%NNy(0, ay, qy, ey) * &
                                            ads%NNz(0, az, qz, ez)
                              load(row + 1) = load(row + 1) - &
                                              source*basis_value*scale
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do

end subroutine AssembleLoadVector

!------------------------------------------------------------------------------
!> Add one contribution to both the library matrix and residual audit list.
!------------------------------------------------------------------------------
subroutine AddTracked(matrix, entries, row, column, value)
   use sparse, only: sparse_matrix, add

   implicit none

   type(sparse_matrix), pointer, intent(inout) :: matrix
   type(MatrixEntries), intent(inout) :: entries
   integer(kind=4), intent(in) :: row, column
   real(kind=8), intent(in) :: value

   if (value == 0.d0) return
   call add(matrix, row, column, value)
   call AppendEntry(entries, row, column, value)

end subroutine AddTracked

subroutine InitializeEntries(entries)
   implicit none
   type(MatrixEntries), intent(inout) :: entries

   entries%count = 0
   entries%capacity = 1024
   allocate (entries%row(entries%capacity))
   allocate (entries%column(entries%capacity))
   allocate (entries%value(entries%capacity))

end subroutine InitializeEntries

subroutine AppendEntry(entries, row, column, value)
   implicit none

   type(MatrixEntries), intent(inout) :: entries
   integer(kind=4), intent(in) :: row, column
   real(kind=8), intent(in) :: value
   integer(kind=4), allocatable :: new_row(:), new_column(:)
   real(kind=8), allocatable :: new_value(:)
   integer(kind=4) :: new_capacity

   if (entries%count == entries%capacity) then
      new_capacity = 2*entries%capacity
      allocate (new_row(new_capacity), new_column(new_capacity))
      allocate (new_value(new_capacity))
      new_row(1:entries%count) = entries%row(1:entries%count)
      new_column(1:entries%count) = entries%column(1:entries%count)
      new_value(1:entries%count) = entries%value(1:entries%count)
      call move_alloc(new_row, entries%row)
      call move_alloc(new_column, entries%column)
      call move_alloc(new_value, entries%value)
      entries%capacity = new_capacity
   end if

   entries%count = entries%count + 1
   entries%row(entries%count) = row
   entries%column(entries%count) = column
   entries%value(entries%count) = value

end subroutine AppendEntry

subroutine ReleaseEntries(entries)
   implicit none
   type(MatrixEntries), intent(inout) :: entries

   if (allocated(entries%row)) deallocate (entries%row)
   if (allocated(entries%column)) deallocate (entries%column)
   if (allocated(entries%value)) deallocate (entries%value)
   entries%count = 0
   entries%capacity = 0

end subroutine ReleaseEntries

!------------------------------------------------------------------------------
!> Compute absolute RMS and relative Euclidean residuals from tracked entries.
!------------------------------------------------------------------------------
subroutine ComputeResidual(entries, right_hand_side, solution, residual_rms, &
                           relative_residual)
   implicit none

   type(MatrixEntries), intent(in) :: entries
   real(kind=8), intent(in) :: right_hand_side(:), solution(:)
   real(kind=8), intent(out) :: residual_rms, relative_residual
   real(kind=8), allocatable :: product(:)
   real(kind=8) :: residual_norm, rhs_norm
   integer(kind=4) :: entry

   allocate (product(size(right_hand_side)))
   product = 0.d0
   do entry = 1, entries%count
      product(entries%row(entry) + 1) = &
         product(entries%row(entry) + 1) + &
         entries%value(entry)*solution(entries%column(entry) + 1)
   end do

   residual_norm = sqrt(sum((right_hand_side - product)**2))
   rhs_norm = sqrt(sum(right_hand_side**2))
   residual_rms = residual_norm/sqrt(real(size(right_hand_side), kind=8))
   relative_residual = residual_norm/max(rhs_norm, tiny(0.d0))
   deallocate (product)

end subroutine ComputeResidual

!------------------------------------------------------------------------------
!> Copy a replicated global coefficient tensor into canonical local layout.
!------------------------------------------------------------------------------
subroutine CopyOwnedCoefficients(ads, global_solution, local_solution)
   use Setup, only: ADS_Setup

   implicit none

   type(ADS_Setup), intent(in) :: ads
   real(kind=8), intent(in) :: global_solution(:)
   real(kind=8), intent(out) :: local_solution(:, :)
   integer(kind=4) :: lx, ly, lz, gx, gy, gz
   integer(kind=4) :: local_column, global_index

   do lz = 1, ads%s(3)
      gz = ads%ibeg(3) + lz - 2
      do ly = 1, ads%s(2)
         gy = ads%ibeg(2) + ly - 2
         local_column = ly + (lz - 1)*ads%s(2)
         do lx = 1, ads%s(1)
            gx = ads%ibeg(1) + lx - 2
            global_index = LinearIndex3D(gx, gy, gz, ads%n)
            local_solution(lx, local_column) = global_solution(global_index + 1)
         end do
      end do
   end do

end subroutine CopyOwnedCoefficients

pure integer(kind=4) function TensorSize(n) result(size_value)
   implicit none
   integer(kind=4), intent(in), dimension(3) :: n
   size_value = (n(1) + 1)*(n(2) + 1)*(n(3) + 1)
end function TensorSize

pure integer(kind=4) function LinearIndex3D(ix, iy, iz, n) result(index_value)
   implicit none
   integer(kind=4), intent(in) :: ix, iy, iz
   integer(kind=4), intent(in), dimension(3) :: n
   index_value = ix + (n(1) + 1)*(iy + (n(2) + 1)*iz)
end function LinearIndex3D

pure logical function IsBoundary3D(ix, iy, iz, n) result(boundary)
   implicit none
   integer(kind=4), intent(in) :: ix, iy, iz
   integer(kind=4), intent(in), dimension(3) :: n
   boundary = ix == 0 .or. ix == n(1) .or. &
              iy == 0 .or. iy == n(2) .or. &
              iz == 0 .or. iz == n(3)
end function IsBoundary3D

end module igrm_mumps_solver
