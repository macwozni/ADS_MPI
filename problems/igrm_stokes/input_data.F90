!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Parameters and manufactured data for the stationary 3D Stokes
!> iGRM problem.
!
!> @details The manufactured solution is the `stokes3_type2` field from the
!> iga-ads DGiGRM example.  On the unit cube it satisfies
!> \f$-\nu\Delta\boldsymbol u+\nabla p=\boldsymbol f\f$ and
!> \f$\nabla\cdot\boldsymbol u=0\f$, with \f$\nu=1\f$.
!
!------------------------------------------------------------------------------
module input_data

   use argument_parser, only: ReadIntegerArgument, RequirePositiveInteger, &
                              RequireSafeSplineDimensions

   implicit none
   private

   real(kind=8), parameter, public :: VISCOSITY_VALUE = 1.d0
   real(kind=8), parameter, public :: PENALTY_FACTOR = 1.d0
   real(kind=8), parameter, public :: RESIDUAL_TOLERANCE = 1.d-8
   integer(kind=4), parameter, public :: MAX_POLYNOMIAL_DEGREE = 9

   integer(kind=4), public :: nelem(3)
   integer(kind=4), public :: p_test(3)
   integer(kind=4), public :: p_trial(3)
   integer(kind=4), public :: procx, procy, procz

   public :: InitializeParameters
   public :: stokes_profile
   public :: stokes_profile_derivative
   public :: stokes_profile_second_derivative
   public :: stokes_profile_third_derivative
   public :: exact_velocity
   public :: exact_velocity_divergence
   public :: exact_vx
   public :: exact_vy
   public :: exact_vz
   public :: exact_pressure
   public :: exact_pressure_gradient
   public :: exact_velocity_laplacian
   public :: body_force
   public :: exact_force

contains

!> @brief Read element counts, test and trial degrees, and the MPI grid.
   subroutine InitializeParameters
      implicit none
      integer(kind=4) :: direction
      integer(kind=8) :: discontinuous_extent, trial_extent, trial_continuity

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
         call RequirePositiveInteger(p_test(direction), &
                                     "test polynomial degree")
         call RequirePositiveInteger(p_trial(direction), &
                                     "trial polynomial degree")
      end do

      if (any(p_test > MAX_POLYNOMIAL_DEGREE)) then
         write(*, '(A,I0)') &
            "unsupported polynomial degree in iGRM test space; " // &
            "maximum supported degree is ", MAX_POLYNOMIAL_DEGREE
         stop 5
      end if

      if (any(p_trial > MAX_POLYNOMIAL_DEGREE)) then
         write(*, '(A,I0)') &
            "unsupported polynomial degree in iGRM trial space; " // &
            "maximum supported degree is ", MAX_POLYNOMIAL_DEGREE
         stop 5
      end if

      if (any(p_test < p_trial)) then
         write(*, '(A)') "test polynomial degree must be greater than " // &
                         "or equal to trial polynomial degree"
         stop 5
      end if

      do direction = 1, 3
         call RequireSafeSplineDimensions(nelem(direction), p_test(direction))
         call RequireSafeSplineDimensions(nelem(direction), p_trial(direction))
         discontinuous_extent = int(nelem(direction), kind=8)* &
            (int(p_test(direction), kind=8) + 1_8) - 1_8
         if (discontinuous_extent > &
             int(huge(nelem(direction)), kind=8)) then
            write(*, '(A)') "discontinuous test-space dimensions exceed " // &
                            "supported integer range"
            stop 5
         end if
         trial_continuity = min(1_8, int(p_trial(direction), kind=8) - 1_8)
         trial_extent = int(nelem(direction), kind=8)* &
            (int(p_trial(direction), kind=8) - trial_continuity) + &
            trial_continuity
         if (trial_extent > int(huge(nelem(direction)), kind=8)) then
            write(*, '(A)') "C1 trial-space dimensions exceed supported " // &
                            "integer range"
            stop 5
         end if
      end do

      call RequirePositiveInteger(procx, "process-grid dimension")
      call RequirePositiveInteger(procy, "process-grid dimension")
      call RequirePositiveInteger(procz, "process-grid dimension")

   end subroutine InitializeParameters

!> @brief Polynomial \f$a(t)=t^2(1-t)^2\f$ used by the exact velocity.
   pure function stokes_profile(coordinate) result(value)
      implicit none
      real(kind=8), intent(in) :: coordinate
      real(kind=8) :: value

      value = coordinate*coordinate*(1.d0 - coordinate)* &
              (1.d0 - coordinate)

   end function stokes_profile

!> @brief First derivative of \ref stokes_profile.
   pure function stokes_profile_derivative(coordinate) result(value)
      implicit none
      real(kind=8), intent(in) :: coordinate
      real(kind=8) :: value

      value = 2.d0*coordinate - 6.d0*coordinate*coordinate + &
              4.d0*coordinate*coordinate*coordinate

   end function stokes_profile_derivative

!> @brief Second derivative of \ref stokes_profile.
   pure function stokes_profile_second_derivative(coordinate) result(value)
      implicit none
      real(kind=8), intent(in) :: coordinate
      real(kind=8) :: value

      value = 2.d0 - 12.d0*coordinate + 12.d0*coordinate*coordinate

   end function stokes_profile_second_derivative

!> @brief Third derivative of \ref stokes_profile.
   pure function stokes_profile_third_derivative(coordinate) result(value)
      implicit none
      real(kind=8), intent(in) :: coordinate
      real(kind=8) :: value

      value = -12.d0 + 24.d0*coordinate

   end function stokes_profile_third_derivative

!> @brief Divergence-free manufactured velocity.
   pure function exact_velocity(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value(3)
      real(kind=8) :: profile(3), derivative(3)
      integer(kind=4) :: direction

      do direction = 1, 3
         profile(direction) = stokes_profile(X(direction))
         derivative(direction) = stokes_profile_derivative(X(direction))
      end do

      value(1) = profile(1)*(derivative(2) - derivative(3))
      value(2) = profile(2)*(derivative(3) - derivative(1))
      value(3) = profile(3)*(derivative(1) - derivative(2))

   end function exact_velocity

!> @brief Analytic divergence of \ref exact_velocity (identically zero).
   pure function exact_velocity_divergence(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value
      real(kind=8) :: derivative(3)
      integer(kind=4) :: direction

      do direction = 1, 3
         derivative(direction) = stokes_profile_derivative(X(direction))
      end do

      value = derivative(1)*(derivative(2) - derivative(3)) + &
              derivative(2)*(derivative(3) - derivative(1)) + &
              derivative(3)*(derivative(1) - derivative(2))

   end function exact_velocity_divergence

!> @brief Scalar x component of \ref exact_velocity.
   pure function exact_vx(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value

      value = stokes_profile(X(1))* &
              (stokes_profile_derivative(X(2)) - &
               stokes_profile_derivative(X(3)))

   end function exact_vx

!> @brief Scalar y component of \ref exact_velocity.
   pure function exact_vy(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value

      value = stokes_profile(X(2))* &
              (stokes_profile_derivative(X(3)) - &
               stokes_profile_derivative(X(1)))

   end function exact_vy

!> @brief Scalar z component of \ref exact_velocity.
   pure function exact_vz(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value

      value = stokes_profile(X(3))* &
              (stokes_profile_derivative(X(1)) - &
               stokes_profile_derivative(X(2)))

   end function exact_vz

!> @brief Mean-zero manufactured pressure.
   pure function exact_pressure(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value
      real(kind=8) :: y_factor, diagonal

      y_factor = 2.d0*X(2) - 1.d0
      diagonal = X(1) + X(2) - X(3)
      value = X(1)*X(1)*y_factor*y_factor*diagonal*diagonal - 5.d0/54.d0

   end function exact_pressure

!> @brief Analytic gradient of \ref exact_pressure.
   pure function exact_pressure_gradient(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value(3)
      real(kind=8) :: x_squared, y_factor, diagonal

      x_squared = X(1)*X(1)
      y_factor = 2.d0*X(2) - 1.d0
      diagonal = X(1) + X(2) - X(3)

      value(1) = 2.d0*X(1)*y_factor*y_factor*diagonal*diagonal + &
                 2.d0*x_squared*y_factor*y_factor*diagonal
      value(2) = 4.d0*x_squared*y_factor*diagonal*diagonal + &
                 2.d0*x_squared*y_factor*y_factor*diagonal
      value(3) = -2.d0*x_squared*y_factor*y_factor*diagonal

   end function exact_pressure_gradient

!> @brief Component-wise Laplacian of \ref exact_velocity.
   pure function exact_velocity_laplacian(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value(3)
      real(kind=8) :: profile(3), derivative(3)
      real(kind=8) :: second_derivative(3), third_derivative(3)
      integer(kind=4) :: direction

      do direction = 1, 3
         profile(direction) = stokes_profile(X(direction))
         derivative(direction) = stokes_profile_derivative(X(direction))
         second_derivative(direction) = &
            stokes_profile_second_derivative(X(direction))
         third_derivative(direction) = &
            stokes_profile_third_derivative(X(direction))
      end do

      value(1) = second_derivative(1)*(derivative(2) - derivative(3)) + &
                 profile(1)*(third_derivative(2) - third_derivative(3))
      value(2) = second_derivative(2)*(derivative(3) - derivative(1)) + &
                 profile(2)*(third_derivative(3) - third_derivative(1))
      value(3) = second_derivative(3)*(derivative(1) - derivative(2)) + &
                 profile(3)*(third_derivative(1) - third_derivative(2))

   end function exact_velocity_laplacian

!> @brief Body force \f$-\nu\Delta\boldsymbol u+\nabla p\f$.
   pure function body_force(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value(3)

      value = -VISCOSITY_VALUE*exact_velocity_laplacian(X) + &
              exact_pressure_gradient(X)

   end function body_force

!> @brief Backward-compatible descriptive alias for \ref body_force.
   pure function exact_force(X) result(value)
      implicit none
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value(3)

      value = body_force(X)

   end function exact_force

end module input_data
