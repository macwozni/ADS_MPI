!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Parameters and manufactured data for the stationary 3D Eriksson
!> iGRM problem.
!
!> @details The problem posed on the unit cube is
!> \f$-\varepsilon\Delta u + \boldsymbol\beta\cdot\nabla u=f\f$, with
!> \f$\varepsilon=10^{-2}\f$, \f$\boldsymbol\beta=(1,1,1)\f$, and
!> homogeneous Dirichlet data.  The manufactured field is a tensor product
!> of one-dimensional outflow-layer profiles.  Its exponential is evaluated
!> as \f$\exp(-(1-x)/\varepsilon)\f$ to avoid the overflow-prone quotient
!> \f$\exp(x/\varepsilon)/\exp(1/\varepsilon)\f$.
!
!------------------------------------------------------------------------------
module input_data

   use argument_parser, only: ReadIntegerArgument, RequireNonnegativeInteger, &
                              RequirePositiveInteger, RequireSafeSplineDimensions

   implicit none
   private

   real(kind=8), parameter, public :: EPSILON_VALUE = 1.d-2
   real(kind=8), parameter, public :: BETA(3) = (/ 1.d0, 1.d0, 1.d0 /)
   real(kind=8), parameter, public :: RESIDUAL_TOLERANCE = 1.d-8

   integer(kind=4), public :: nelem(3)
   integer(kind=4), public :: p_test(3)
   integer(kind=4), public :: p_trial(3)
   integer(kind=4), public :: procx, procy, procz

   public :: InitializeParameters
   public :: eriksson_profile
   public :: eriksson_profile_derivative
   public :: eriksson_profile_second_derivative
   public :: exact_solution
   public :: exact_forcing

contains

!> @brief Read three element counts, three test degrees, three trial degrees,
!> and three MPI process-grid dimensions from the command line.
   subroutine InitializeParameters
      implicit none
      integer(kind=4) :: direction

      if (command_argument_count() /= 12) then
         write(*, '(A)') "proper usage with arguments: " // &
            "<nelem_x> <nelem_y> <nelem_z> " // &
            "<ptest_x> <ptest_y> <ptest_z> " // &
            "<ptrial_x> <ptrial_y> <ptrial_z> " // &
            "<procx> <procy> <procz>"
         stop 5
      end if

      do direction = 1, 3
         call ReadIntegerArgument(direction, nelem(direction))
         call ReadIntegerArgument(direction + 3, p_test(direction))
         call ReadIntegerArgument(direction + 6, p_trial(direction))
      end do
      call ReadIntegerArgument(10, procx)
      call ReadIntegerArgument(11, procy)
      call ReadIntegerArgument(12, procz)

      do direction = 1, 3
         call RequirePositiveInteger(nelem(direction), "number of elements")
         call RequireNonnegativeInteger(p_test(direction), "test polynomial degree")
         call RequirePositiveInteger(p_trial(direction), "trial polynomial degree")
      end do

      if (any(p_test <= p_trial)) then
         write(*, '(A)') &
            "test polynomial degree must be greater than trial polynomial degree"
         stop 5
      end if

      do direction = 1, 3
         call RequireSafeSplineDimensions(nelem(direction), p_test(direction))
         call RequireSafeSplineDimensions(nelem(direction), p_trial(direction))
      end do
      call RequirePositiveInteger(procx, "process-grid dimension")
      call RequirePositiveInteger(procy, "process-grid dimension")
      call RequirePositiveInteger(procz, "process-grid dimension")

   end subroutine InitializeParameters

!> @brief Stable one-dimensional Eriksson outflow-layer profile.
   pure function eriksson_profile(coordinate) result(value)
      implicit none
      real(kind=8), intent(in) :: coordinate
      real(kind=8) :: value
      real(kind=8), parameter :: LEFT_EXPONENTIAL = exp(-1.d0/EPSILON_VALUE)
      real(kind=8), parameter :: NORMALIZER = 1.d0 - LEFT_EXPONENTIAL
      real(kind=8) :: outflow_exponential

      outflow_exponential = exp(-(1.d0 - coordinate)/EPSILON_VALUE)
      value = coordinate - (outflow_exponential - LEFT_EXPONENTIAL)/NORMALIZER

   end function eriksson_profile

!> @brief First derivative of \ref eriksson_profile.
   pure function eriksson_profile_derivative(coordinate) result(value)
      implicit none
      real(kind=8), intent(in) :: coordinate
      real(kind=8) :: value
      real(kind=8), parameter :: NORMALIZER = 1.d0 - exp(-1.d0/EPSILON_VALUE)

      value = 1.d0 - exp(-(1.d0 - coordinate)/EPSILON_VALUE) / &
                        (EPSILON_VALUE*NORMALIZER)

   end function eriksson_profile_derivative

!> @brief Second derivative of \ref eriksson_profile.
   pure function eriksson_profile_second_derivative(coordinate) result(value)
      implicit none
      real(kind=8), intent(in) :: coordinate
      real(kind=8) :: value
      real(kind=8), parameter :: NORMALIZER = 1.d0 - exp(-1.d0/EPSILON_VALUE)

      value = -exp(-(1.d0 - coordinate)/EPSILON_VALUE) / &
               (EPSILON_VALUE*EPSILON_VALUE*NORMALIZER)

   end function eriksson_profile_second_derivative

!> @brief Exact tensor-product solution at a point in the unit cube.
   pure function exact_solution(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value

      value = eriksson_profile(X(1))*eriksson_profile(X(2))* &
              eriksson_profile(X(3))

   end function exact_solution

!> @brief Cancellation-free manufactured volume forcing.
!>
!> Since \f$-\varepsilon g''+g'=1\f$, applying the 3D operator to the
!> tensor product gives the sum of its three two-dimensional factors.
   pure function exact_forcing(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value
      real(kind=8) :: profile(3)

      profile(1) = eriksson_profile(X(1))
      profile(2) = eriksson_profile(X(2))
      profile(3) = eriksson_profile(X(3))
      value = profile(2)*profile(3) + profile(1)*profile(3) + &
              profile(1)*profile(2)

   end function exact_forcing

end module input_data
