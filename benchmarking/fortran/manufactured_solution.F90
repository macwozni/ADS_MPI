! Manufactured transients shared by every benchmark adapter.
!
! The temporal-polynomial case is exactly representable by every supported
! trial space with degree at least three.  The spatial-cosine case is not
! polynomial and is reserved for mesh/degree convergence measurements.
!
! q(s) = s^2(3-2s), Q(x,y,z) = q(x)q(y)q(z),
! u(x,y,z,t) = exp(-t)Q(x,y,z).
!
! Since q'(0)=q'(1)=0, the normal flux is zero on the unit cube.  With the
! library convention M u_t + K u = f, K representing -Delta, the source is
!
!   f = exp(-t) [-Q + (12x-6)q(y)q(z)
!                    + (12y-6)q(x)q(z)
!                    + (12z-6)q(x)q(y)].
!
! R(x,y,z) = cos(pi*x)cos(pi*y)cos(pi*z),
! u(x,y,z,t) = exp(-t)R(x,y,z).
!
! Its normal derivative also vanishes on every face and
!
!   f = (3*pi^2 - 1) exp(-t) R.
module manufactured_solution

   use benchmark_contract, ONLY: ManufacturedCase

   implicit none

   private
   public :: ExactAtEvaluationTime
   public :: ExactSolution
   public :: InitialProjectionRhs
   public :: ManufacturedRhsPoint
   public :: ManufacturedScalarSource
   public :: ManufacturedSource
   public :: ActivateManufacturedCase
   public :: RegisterManufacturedCase
   public :: RegisterSpatialCosineCase
   public :: RegisterTemporalPolynomialCase
   public :: SetCallbackSourceTime
   public :: SetEvaluationTime
   public :: SetStepClock
   public :: SourceEvaluationTime

   real(kind=8) :: evaluation_time = 0.d0
   real(kind=8) :: step_start_time = 0.d0
   real(kind=8) :: step_size = 0.d0
   real(kind=8) :: callback_source_time = 0.d0
   character(len=16) :: active_scheme = ''
   type(ManufacturedCase), save :: active_case
   real(kind=8), parameter :: PI = &
      3.141592653589793238462643383279502884197d0
!$omp threadprivate(callback_source_time)

contains

   pure function ShapeQ(s) result(value)
      real(kind=8), intent(in) :: s
      real(kind=8) :: value

      value = s*s*(3.d0 - 2.d0*s)

   end function ShapeQ

   pure function TemporalPolynomialValue(physical_time, point) result(value)
      real(kind=8), intent(in) :: physical_time
      real(kind=8), dimension(3), intent(in) :: point
      real(kind=8) :: value

      value = exp(-physical_time)*ShapeQ(point(1))*ShapeQ(point(2))* &
              ShapeQ(point(3))

   end function TemporalPolynomialValue

   pure function TemporalPolynomialSource(physical_time, point) result(value)
      real(kind=8), intent(in) :: physical_time
      real(kind=8), dimension(3), intent(in) :: point
      real(kind=8) :: value
      real(kind=8) :: qx, qy, qz

      qx = ShapeQ(point(1))
      qy = ShapeQ(point(2))
      qz = ShapeQ(point(3))
      value = exp(-physical_time)*( &
         -qx*qy*qz + &
         (12.d0*point(1) - 6.d0)*qy*qz + &
         (12.d0*point(2) - 6.d0)*qx*qz + &
         (12.d0*point(3) - 6.d0)*qx*qy)

   end function TemporalPolynomialSource

   subroutine RegisterTemporalPolynomialCase(case_descriptor)
      type(ManufacturedCase), intent(out) :: case_descriptor

      case_descriptor%name = 'temporal-polynomial'
      case_descriptor%initial_l2_norm = (13.d0/35.d0)**1.5d0
      case_descriptor%exact_initial_projection = .true.
      case_descriptor%value => TemporalPolynomialValue
      case_descriptor%source => TemporalPolynomialSource

   end subroutine RegisterTemporalPolynomialCase

   pure function SpatialCosineValue(physical_time, point) result(value)
      real(kind=8), intent(in) :: physical_time
      real(kind=8), dimension(3), intent(in) :: point
      real(kind=8) :: value

      value = exp(-physical_time)*cos(PI*point(1))*cos(PI*point(2))* &
              cos(PI*point(3))

   end function SpatialCosineValue

   pure function SpatialCosineSource(physical_time, point) result(value)
      real(kind=8), intent(in) :: physical_time
      real(kind=8), dimension(3), intent(in) :: point
      real(kind=8) :: value

      value = (3.d0*PI*PI - 1.d0)*SpatialCosineValue(physical_time, point)

   end function SpatialCosineSource

   subroutine RegisterSpatialCosineCase(case_descriptor)
      type(ManufacturedCase), intent(out) :: case_descriptor

      case_descriptor%name = 'spatial-cosine'
      case_descriptor%initial_l2_norm = 1.d0/sqrt(8.d0)
      case_descriptor%exact_initial_projection = .false.
      case_descriptor%value => SpatialCosineValue
      case_descriptor%source => SpatialCosineSource

   end subroutine RegisterSpatialCosineCase

   subroutine RegisterManufacturedCase(case_name, case_descriptor, status)
      character(len=*), intent(in) :: case_name
      type(ManufacturedCase), intent(out) :: case_descriptor
      integer(kind=4), intent(out) :: status

      status = 0
      select case (trim(case_name))
      case ('temporal-polynomial')
         call RegisterTemporalPolynomialCase(case_descriptor)
      case ('spatial-cosine')
         call RegisterSpatialCosineCase(case_descriptor)
      case default
         status = 5
      end select

   end subroutine RegisterManufacturedCase

   subroutine ActivateManufacturedCase(case_descriptor, status)
      type(ManufacturedCase), intent(in) :: case_descriptor
      integer(kind=4), intent(out) :: status

      status = 0
      if (len_trim(case_descriptor%name) == 0 .or. &
          .not. associated(case_descriptor%value) .or. &
          .not. associated(case_descriptor%source)) then
         status = 5
         return
      end if
      active_case%name = case_descriptor%name
      active_case%initial_l2_norm = case_descriptor%initial_l2_norm
      active_case%exact_initial_projection = &
         case_descriptor%exact_initial_projection
      active_case%value => case_descriptor%value
      active_case%source => case_descriptor%source

   end subroutine ActivateManufacturedCase

   function ExactSolution(physical_time, point) result(value)
      real(kind=8), intent(in) :: physical_time
      real(kind=8), dimension(3), intent(in) :: point
      real(kind=8) :: value

      value = active_case%value(physical_time, point)

   end function ExactSolution

   function ManufacturedSource(physical_time, point) result(value)
      real(kind=8), intent(in) :: physical_time
      real(kind=8), dimension(3), intent(in) :: point
      real(kind=8) :: value

      value = active_case%source(physical_time, point)

   end function ManufacturedSource

   subroutine SetEvaluationTime(physical_time)
      real(kind=8), intent(in) :: physical_time

      evaluation_time = physical_time

   end subroutine SetEvaluationTime

   function ExactAtEvaluationTime(point) result(value)
      real(kind=8), dimension(3), intent(in) :: point
      real(kind=8) :: value

      value = ExactSolution(evaluation_time, point)

   end function ExactAtEvaluationTime

   subroutine SetStepClock(scheme, physical_start, dt)
      character(len=*), intent(in) :: scheme
      real(kind=8), intent(in) :: physical_start, dt

      active_scheme = trim(scheme)
      step_start_time = physical_start
      step_size = dt

   end subroutine SetStepClock

   function SourceEvaluationTime(substep) result(physical_time)
      integer(kind=4), intent(in) :: substep
      real(kind=8) :: physical_time

      select case (trim(active_scheme))
      case ('dg')
         physical_time = step_start_time + 0.5d0*step_size
      case ('be')
         physical_time = step_start_time + step_size
      case ('pr')
         physical_time = step_start_time + &
                         (real(substep, kind=8) - 0.5d0)*step_size/3.d0
      case default
         physical_time = -huge(1.d0)
      end select

   end function SourceEvaluationTime

   function ManufacturedScalarSource(un, gradient, point) result(value)
      real(kind=8), intent(in) :: un
      real(kind=8), dimension(3), intent(in) :: gradient, point
      real(kind=8) :: value

      value = ManufacturedSource(callback_source_time, point)

   end function ManufacturedScalarSource

   subroutine SetCallbackSourceTime(substep)
      integer(kind=4), intent(in) :: substep

      callback_source_time = SourceEvaluationTime(substep)

   end subroutine SetCallbackSourceTime

   subroutine InitialProjectionRhs( &
      ads, point, quadrature, element, local_basis, gradient, history, &
      un11, un13, un23, ads_data, jacobian, weight, direction, substep, &
      alpha_step, forcing_cb, value)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use Interfaces, ONLY: forcing_fun
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
      real(kind=8) :: test_value

      test_value = ads%NNx(0, local_basis(1), quadrature(1), element(1))* &
                   ads%NNy(0, local_basis(2), quadrature(2), element(2))* &
                   ads%NNz(0, local_basis(3), quadrature(3), element(3))
      value = jacobian*weight*test_value*ExactSolution(0.d0, point)

   end subroutine InitialProjectionRhs

   subroutine ManufacturedRhsPoint( &
      ads, point, quadrature, element, local_basis, gradient, history, &
      un11, un13, un23, ads_data, jacobian, weight, direction, substep, &
      alpha_step, forcing_cb, value)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use Interfaces, ONLY: forcing_fun
      use RHS_eq, ONLY: ComputePointForRHS
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

      ! Form3DRHS invokes this routine inside an OpenMP region.  The source
      ! time is therefore thread-private, while the step clock is immutable
      ! for the duration of the complete MultiStep call.
      call SetCallbackSourceTime(substep)
      call ComputePointForRHS( &
         ads, point, quadrature, element, local_basis, gradient, history, &
         un11, un13, un23, ads_data, jacobian, weight, direction, substep, &
         alpha_step, ManufacturedScalarSource, value)

   end subroutine ManufacturedRhsPoint

end module manufactured_solution
