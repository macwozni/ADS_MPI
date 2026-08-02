module workflow_test_support
   implicit none

   integer(kind=4), parameter :: MAX_SOLVES = 12
   integer(kind=4) :: solve_call_count = 0
   integer(kind=4) :: solve_status(MAX_SOLVES) = 0
   integer(kind=4) :: form_un_call_count = 0
   integer(kind=4) :: form_un_substeps(3) = 0
   integer(kind=4) :: form_rhs_call_count = 0
   integer(kind=4) :: normalize_call_count = 0
   integer(kind=4) :: distribute_call_count = 0

contains

   subroutine reset_workflow_spy()
      solve_call_count = 0
      solve_status = 0
      form_un_call_count = 0
      form_un_substeps = 0
      form_rhs_call_count = 0
      normalize_call_count = 0
      distribute_call_count = 0
   end subroutine reset_workflow_spy


   subroutine fail_solve(call_number, status)
      integer(kind=4), intent(in) :: call_number, status

      if (call_number < 1 .or. call_number > MAX_SOLVES) then
         error stop 'invalid injected solve call number'
      end if
      solve_status(call_number) = status
   end subroutine fail_solve


   real(kind=8) function forcing(un, du, X) result(value)
      real(kind=8), intent(in) :: un, du(3), X(3)

      value = un + sum(du) + sum(X)
   end function forcing

end module workflow_test_support


module mpi
   implicit none
end module mpi


module sparse
   implicit none
end module sparse


module ads_lifecycle
   implicit none
   integer(kind=4), parameter :: initialize = 0
   integer(kind=4), parameter :: initialize_setup = 0
   integer(kind=4), parameter :: ComputeDecomposition = 0
   integer(kind=4), parameter :: AllocateADSdata = 0
   integer(kind=4), parameter :: AllocateADS = 0
   integer(kind=4), parameter :: Cleanup_data = 0
   integer(kind=4), parameter :: Cleanup_ADS = 0
end module ads_lifecycle


module mumps_solver
   implicit none
   integer(kind=4), parameter :: SolveOneDirection = 0
end module mumps_solver


module solution_output
   implicit none
   integer(kind=4), parameter :: PrintSolution = 0
end module solution_output


module projection_engine
   use Setup, only: ADS_Setup, ADS_compute_data
   use Interfaces, only: forcing_fun, rhs_point_fun
   use workflow_test_support, only: form_un_call_count, form_un_substeps, &
      form_rhs_call_count
   implicit none
   integer(kind=4), parameter :: ComputeMatrix = 0

contains

   subroutine FormUn(substep, ads, ads_data)
      integer(kind=4), intent(in) :: substep
      type(ADS_Setup), intent(in) :: ads
      type(ADS_compute_data), intent(inout) :: ads_data

      form_un_call_count = form_un_call_count + 1
      if (form_un_call_count <= size(form_un_substeps)) then
         form_un_substeps(form_un_call_count) = substep
      end if
   end subroutine FormUn


   subroutine Form3DRHS(ads_test, ads_trial, ads_data, direction, n, substep, &
                        alpha_step, forcing, igrm, rhs_point)
      type(ADS_Setup), intent(in) :: ads_test, ads_trial
      type(ADS_compute_data), intent(inout) :: ads_data
      integer(kind=4), intent(in) :: direction(3), n, substep
      real(kind=8), intent(in) :: alpha_step(7, 3)
      procedure(forcing_fun) :: forcing
      logical, intent(out) :: igrm
      procedure(rhs_point_fun), optional :: rhs_point

      form_rhs_call_count = form_rhs_call_count + 1
      igrm = any(direction == 1)
      if (allocated(ads_data%F)) ads_data%F = -10.0d0*form_rhs_call_count
      if (allocated(ads_data%Ft)) ads_data%Ft = -20.0d0*form_rhs_call_count
   end subroutine Form3DRHS

end module projection_engine


module ads_directional_solve
   use Setup, only: ADS_Setup
   use workflow_test_support, only: solve_call_count, solve_status
   implicit none

contains

   subroutine solve_problem(ads_test, ads_trial, a, b, c, mixA, mixB, mixBT, &
                            direction, igrm, F, F2, Ft, Ft2, ierr)
      type(ADS_Setup), intent(in) :: ads_test, ads_trial
      integer(kind=4), intent(in) :: a, b, c, direction(3)
      real(kind=8), intent(in) :: mixA(4), mixB(4), mixBT(4)
      logical, intent(in) :: igrm
      real(kind=8), allocatable, intent(inout) :: F(:, :), F2(:, :)
      real(kind=8), allocatable, intent(inout) :: Ft(:, :), Ft2(:, :)
      integer(kind=4), intent(out) :: ierr

      solve_call_count = solve_call_count + 1
      ierr = solve_status(solve_call_count)
      if (ierr /= 0) return

      if (allocated(F2)) F2 = real(solve_call_count, kind=8)
      if (igrm .and. allocated(Ft2)) Ft2 = real(solve_call_count, kind=8)
   end subroutine solve_problem

end module ads_directional_solve


module reorderRHS
   use Setup, only: ADS_Setup
   use workflow_test_support, only: normalize_call_count
   implicit none
   integer(kind=4), parameter :: ReorderRHSForX = 0
   integer(kind=4), parameter :: ReorderRHSForY = 0
   integer(kind=4), parameter :: ReorderRHSForZ = 0

contains

   subroutine NormalizeTrialBufferToXYZ(ads, order, F)
      type(ADS_Setup), intent(in) :: ads
      integer(kind=4), intent(in) :: order(3)
      real(kind=8), allocatable, intent(inout) :: F(:, :)

      normalize_call_count = normalize_call_count + 1
   end subroutine NormalizeTrialBufferToXYZ

end module reorderRHS


module my_mpi
   use Setup, only: ADS_Setup, ADS_compute_data
   use workflow_test_support, only: distribute_call_count
   implicit none
   integer(kind=4), parameter :: Gather = 0
   integer(kind=4), parameter :: Scatter = 0

contains

   subroutine DistributeSpline(part, ads_trial, ads_data)
      real(kind=8), intent(in) :: part(:, :)
      type(ADS_Setup), intent(in) :: ads_trial
      type(ADS_compute_data), intent(inout) :: ads_data

      distribute_call_count = distribute_call_count + 1
   end subroutine DistributeSpline

end module my_mpi
