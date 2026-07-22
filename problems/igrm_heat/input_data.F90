!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Input parameters for the iGRM heat example.
!>
!> @details
!> This module stores command-line grid dimensions, test/trial polynomial
!> orders, MPI process-grid dimensions, time-loop parameters, and the
!> selected iGRM time scheme. It also defines the compactly supported
!> initial state used by the heat problem.
!
!------------------------------------------------------------------------------
module input_data

   use argument_parser, ONLY: ARG_TIME_SCHEME_BE => TIME_SCHEME_BE, &
                              ARG_TIME_SCHEME_DG => TIME_SCHEME_DG, &
                              ARG_TIME_SCHEME_PR => TIME_SCHEME_PR, &
                              ReadIntegerArgument, ReadRealArgument, &
                              ReadTimeSchemeArgument, RequireNonnegativeInteger, &
                              RequirePositiveInteger, RequirePositiveReal, &
                              RequireSafeSplineDimensions, SelectTimeScheme

   implicit none

!> @brief Numbers of elements in the three parametric directions.
   integer(kind = 4), dimension(3) :: nelem

!> @brief Polynomial orders for the enriched iGRM test space.
   integer(kind = 4), dimension(3) :: p_test

!> @brief Polynomial orders for the trial space.
   integer(kind = 4), dimension(3) :: p_trial

!> @brief Numbers of MPI processes in the three process-grid directions.
   integer(kind = 4) :: procx, procy, procz

!> @brief Current physical time.
   real(kind = 8) :: t

!> @brief Time-step size.
   real(kind = 8) :: Dt

!> @brief Number of physical time steps after the initial projection.
   integer(kind = 4) :: steps

!> @brief Supported iGRM time-scheme identifiers.
   integer(kind = 4), parameter :: TIME_SCHEME_DG = ARG_TIME_SCHEME_DG
   integer(kind = 4), parameter :: TIME_SCHEME_PR = ARG_TIME_SCHEME_PR
   integer(kind = 4), parameter :: TIME_SCHEME_BE = ARG_TIME_SCHEME_BE

!> @brief Time scheme selected for the iGRM multi-step driver.
   character(len = 32) :: time_scheme = "dg"

!> @brief Normalized iGRM time-scheme identifier selected from arguments.
   integer(kind = 4) :: selected_time_scheme = TIME_SCHEME_DG

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads iGRM heat parameters from command-line arguments.
!>
!> @details
!> The accepted argument lists are:
!> `<nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz> <steps> <dt>`,
!> or the same list followed by `<scheme>`.
!
!---------------------------------------------------------------------------
   subroutine InitializeParameters
      implicit none
      integer(kind = 4) :: nargs

      nargs = COMMAND_ARGUMENT_COUNT()

      if (nargs .NE. 14 .AND. nargs .NE. 15) then
         write(*,*) "proper usage with arguments: ", &
         "<nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> ", &
         "<ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz> ", &
         "<steps> <dt> [dg|pr|be]"
         STOP 5
      end if

      call ReadIntegerArgument(1, nelem(1))
      call ReadIntegerArgument(2, nelem(2))
      call ReadIntegerArgument(3, nelem(3))
      call ReadIntegerArgument(4, p_test(1))
      call ReadIntegerArgument(5, p_test(2))
      call ReadIntegerArgument(6, p_test(3))
      call ReadIntegerArgument(7, p_trial(1))
      call ReadIntegerArgument(8, p_trial(2))
      call ReadIntegerArgument(9, p_trial(3))
      call ReadIntegerArgument(10, procx)
      call ReadIntegerArgument(11, procy)
      call ReadIntegerArgument(12, procz)
      call ReadIntegerArgument(13, steps)
      call ReadRealArgument(14, Dt)
      call SelectTimeScheme("dg", selected_time_scheme, time_scheme)
      if (nargs .EQ. 15) call ReadTimeSchemeArgument(15, selected_time_scheme, time_scheme)

      call RequirePositiveInteger(nelem(1), "number of elements")
      call RequirePositiveInteger(nelem(2), "number of elements")
      call RequirePositiveInteger(nelem(3), "number of elements")
      call RequireNonnegativeInteger(p_test(1), "test polynomial degree")
      call RequireNonnegativeInteger(p_test(2), "test polynomial degree")
      call RequireNonnegativeInteger(p_test(3), "test polynomial degree")
      call RequireNonnegativeInteger(p_trial(1), "trial polynomial degree")
      call RequireNonnegativeInteger(p_trial(2), "trial polynomial degree")
      call RequireNonnegativeInteger(p_trial(3), "trial polynomial degree")
      if (any(p_test <= p_trial)) then
         write(*, '(A)') "test polynomial degree must be greater than trial polynomial degree"
         STOP 5
      end if
      call RequireSafeSplineDimensions(nelem(1), p_test(1))
      call RequireSafeSplineDimensions(nelem(2), p_test(2))
      call RequireSafeSplineDimensions(nelem(3), p_test(3))
      call RequireSafeSplineDimensions(nelem(1), p_trial(1))
      call RequireSafeSplineDimensions(nelem(2), p_trial(2))
      call RequireSafeSplineDimensions(nelem(3), p_trial(3))
      call RequirePositiveInteger(procx, "process-grid dimension")
      call RequirePositiveInteger(procy, "process-grid dimension")
      call RequirePositiveInteger(procz, "process-grid dimension")
      call RequireNonnegativeInteger(steps, "number of time steps")
      call RequirePositiveReal(Dt, "time step")

   end subroutine InitializeParameters

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the initial heat state \f$u(0)\f$.
!>
!> @details
!> The profile matches the standard heat example: a compact-support bump
!> centered inside the unit cube.
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
      use math, ONLY: bump3d
      implicit none
      real(kind = 8), intent(in) :: x, y, z
      real(kind = 8) :: val

      val = 2.d0*bump3d(0.05d0, 0.4d0, x, y, z)

   end function initial_state

end module input_data
