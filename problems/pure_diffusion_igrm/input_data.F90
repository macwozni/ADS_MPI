!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Input data and placeholder coefficients for the pure-diffusion
!> iGRM example.
!>
!> @details
!> This module stores discretization, time-step, and process-grid
!> parameters for the pure-diffusion iGRM driver. The pointwise forcing
!> callback required by the current ADS API is defined in \ref RHS_fun.
!> Placeholder functions for boundary/source data are retained from the older
!> problem formulation.
!
!------------------------------------------------------------------------------
module input_data

   use argument_parser, ONLY: ReadIntegerArgument, ReadRealArgument, ReadTimeSchemeArgument, &
                              RequireNonnegativeInteger, RequirePositiveInteger, RequirePositiveReal, &
                              RequireSafeSplineDimensions

   implicit none

!> @brief Number of generated curves and number of segments per curve
!> retained from the template.
   integer(kind = 4), parameter :: cN = 30, cL = 16
!> @brief Radius and source/sink strengths retained from the template
!> problem.
   real (kind = 8), parameter :: radius = 0.15, pumping_strength = 1, draining_strength = 1

!> @brief Curve-coordinate buffers retained from the template problem.
   real (kind = 8) :: cx(cN * cL), cy(cN * cL), cz(cN * cL)
!> @brief Nonlinear coefficient retained from the template problem.
   real (kind = 8) :: mi = 10.d0
!> @brief Ground-level parameter retained from the template problem.
   real (kind = 8) :: GROUND = 0.2
!> @brief Minimum and maximum coefficient values retained from the template.
   real (kind = 8), parameter :: Kqmin = 1.d0, Kqmax = 1000.d0
!> @brief Numbers of pump and drain points retained from the template
!> problem.
   integer(kind = 4) :: npumps, ndrains
!> @brief Pump and drain coordinate buffers retained from the template
!> problem.
   real (kind = 8), allocatable, dimension(:,:) :: pumps, drains

!> @brief Permeability cache retained from the template problem.
   real (kind = 8), allocatable :: Kqvals(:,:,:,:,:,:)

!> @brief Current physical time.
   real (kind = 8) :: t

!> @brief Time-step size.
   real (kind = 8) :: Dt

!> @brief Number of time iterations.
   integer :: steps

!> @brief Accumulated statistic retained from the template problem.
   real (kind = 8) :: pollution = 0

!> @brief Polynomial order of the approximation space.
   integer(kind = 4) :: ORDER

!> @brief Number of elements in each parametric direction.
   integer(kind = 4) :: SIZE

!> @brief Numbers of MPI processes in the three process-grid directions.
   integer(kind = 4) :: procx, procy, procz

!> @brief Time scheme selected for the iGRM multi-step driver.
   character(len = 32) :: time_scheme = "dg"

!> @brief Drained quantity retained from the template problem.
   real (kind = 8) :: drained = 0


contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads discretization, process-grid, and time parameters.
!>
!> @details
!> The expected argument list is:
!> `<size> <order> <procx> <procy> <procz> <steps> <dt> [scheme]`.
!
!---------------------------------------------------------------------------
   subroutine InitializeParameters
      implicit none
      integer(kind = 4) :: selected_time_scheme

      ! ./l2 <size> <procx> <procy> <procz> <nsteps> <dt>
      ORDER = 2

      if (COMMAND_ARGUMENT_COUNT() .NE. 7 .AND. COMMAND_ARGUMENT_COUNT() .NE. 8) then
         write(*,*) "proper usage with arguments: ", &
         "<size> <order> <procx> <procy> <procz> <steps> <dt> [dg|pr|be]"
         STOP 5
      end if

      call ReadIntegerArgument(1, SIZE)
      call ReadIntegerArgument(2, ORDER)
      call ReadIntegerArgument(3, procx)
      call ReadIntegerArgument(4, procy)
      call ReadIntegerArgument(5, procz)
      call ReadIntegerArgument(6, steps)
      call ReadRealArgument(7, Dt)
      time_scheme = "dg"
      if (COMMAND_ARGUMENT_COUNT() .EQ. 8) then
         call ReadTimeSchemeArgument(8, selected_time_scheme, time_scheme)
      end if

      call RequirePositiveInteger(SIZE, "number of elements")
      call RequireNonnegativeInteger(ORDER, "polynomial order")
      call RequireSafeSplineDimensions(SIZE, ORDER, 1)
      call RequireNonnegativeInteger(steps, "number of time steps")
      call RequirePositiveReal(Dt, "time step")
      call RequirePositiveInteger(procx, "process-grid dimension")
      call RequirePositiveInteger(procy, "process-grid dimension")
      call RequirePositiveInteger(procz, "process-grid dimension")

   end subroutine InitializeParameters

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Placeholder boundary source function.
!>
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
!> @return fval
!> Boundary source value.
!
!---------------------------------------------------------------------------
   function g(x, y, z) result (fval)
      implicit none
      real (kind = 8) :: x, y, z
      real (kind = 8) :: fval

      fval = 0.d0

   end function g

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Placeholder advective boundary coefficient.
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
!> @return fval
!> Boundary coefficient value.
!
!---------------------------------------------------------------------------
   function b(x, y, z) result (fval)
      implicit none
      real (kind = 8) :: x, y, z
      real (kind = 8) :: fval

      fval = 0.d0

   end function b

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Placeholder boundary normal factor.
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
!> @return fval
!> Boundary normal factor.
!
!---------------------------------------------------------------------------
   function n(x, y, z) result (fval)
      implicit none
      real (kind = 8) :: x, y, z
      real (kind = 8) :: fval

      fval = 0.d0

   end function n


end module input_data
