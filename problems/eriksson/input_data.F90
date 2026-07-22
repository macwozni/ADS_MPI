!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Input data for the Eriksson example.
!>
!> @details
!> This module stores the spatial discretization, time-step parameters,
!> and MPI process-grid dimensions used by the Eriksson driver. It also
!> defines the compactly supported initial state. The pointwise forcing
!> callback expected by the current ADS API is defined in \ref RHS_fun.
!
!------------------------------------------------------------------------------
module input_data

   use argument_parser, ONLY: ReadIntegerArgument, ReadRealArgument, &
                              RequireNonnegativeInteger, RequirePositiveInteger, RequirePositiveReal, &
                              RequireSafeSplineDimensions

   implicit none

!> @brief Current physical time.
   real (kind = 8) :: t

!> @brief Time-step size.
   real (kind = 8) :: Dt

!> @brief Number of time iterations.
   integer :: steps

!> @brief Polynomial order of the approximation space.
   integer(kind = 4) :: ORDER

!> @brief Number of elements in each parametric direction.
   integer(kind = 4) :: SIZE

!> @brief Numbers of MPI processes in the three process-grid directions.
   integer(kind = 4) :: procx, procy, procz



contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads discretization, time, and process-grid parameters.
!>
!> @details
!> The expected argument list is:
!> `<size> <order> <steps> <dt> <procx> <procy> <procz>`.
!
!---------------------------------------------------------------------------
   subroutine InitializeParameters
      implicit none

      ! ./l2 <size> <order> <procx> <procy> <procz>

      if (COMMAND_ARGUMENT_COUNT() .NE. 7) then
         write(*,*) "proper usage with arguments: ", &
         "<size> <order> <steps> <dt> <procx> <procy> <procz>"
         STOP 5
      end if

      call ReadIntegerArgument(1, SIZE)
      call ReadIntegerArgument(2, ORDER)
      call ReadIntegerArgument(3, steps)
      call ReadRealArgument(4, Dt)
      call ReadIntegerArgument(5, procx)
      call ReadIntegerArgument(6, procy)
      call ReadIntegerArgument(7, procz)

      call RequirePositiveInteger(SIZE, "number of elements")
      call RequireNonnegativeInteger(ORDER, "polynomial order")
      call RequireSafeSplineDimensions(SIZE, ORDER)
      call RequireNonnegativeInteger(steps, "number of time steps")
      call RequirePositiveReal(Dt, "time step")
      call RequirePositiveInteger(procx, "process-grid dimension")
      call RequirePositiveInteger(procy, "process-grid dimension")
      call RequirePositiveInteger(procz, "process-grid dimension")

   end subroutine InitializeParameters

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the initial state \f$u(0)\f$.
!>
!> @details
!> The current profile is a scaled compact-support bump centered inside
!> the unit cube.
!
! Input:
! ------
!> @param[in] x
!> First physical coordinate.
!>
!> @param[in] y
!> Second physical coordinate.
!>
!> @param[in] z
!> Third physical coordinate.
!
! Output:
! -------
!> @return val
!> Initial-state value at the supplied point.
!
!---------------------------------------------------------------------------
   function initial_state(x, y, z) result (val)
      use math, ONLY: falloff, bump3d, lerp
      implicit none
      real (kind = 8), intent(in) :: x, y, z
      real (kind = 8) :: dist, val
      real (kind = 8), dimension(3) :: p1

      p1 = (/ x, y, z/)
      !dist = sqrt(dist_from_curves(p1, cx, cy, cz, cN, cL))
      dist = 0 !!!!!!
      !val = 1.d0 * lerp(falloff(0.d0, 0.2d0, dist), 0.d0, 1.d0) * bump3d(0.2d0, 0.6d0, x, y, z)
      val = 2.d0  * bump3d(0.05d0, 0.4d0, x, y, z)

   end function initial_state

end module input_data
