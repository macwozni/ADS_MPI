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
      character(100) :: input

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
      if (COMMAND_ARGUMENT_COUNT() .EQ. 8) call ReadStringArgument(8, time_scheme)

   contains

      subroutine ReadArgument(arg, input)
         implicit none
         integer(kind = 4), intent(in) :: arg
         character(*), intent(out) :: input
         integer(kind = 4) :: length
         integer(kind = 4) :: status

         call GET_COMMAND_ARGUMENT(arg, input, length, status)
         if (status /= 0) then
            write(*,*) "invalid command argument: ", arg
            STOP 5
         end if
      end subroutine ReadArgument

      subroutine ReadIntegerArgument(arg, value)
         implicit none
         integer(kind = 4), intent(in) :: arg
         integer(kind = 4), intent(out) :: value
         character(100) :: input
         integer(kind = 4) :: read_status

         call ReadArgument(arg, input)
         read(input, *, iostat = read_status) value
         if (read_status /= 0) then
            write(*,*) "invalid integer argument: ", arg
            STOP 5
         end if
      end subroutine ReadIntegerArgument

      subroutine ReadRealArgument(arg, value)
         implicit none
         integer(kind = 4), intent(in) :: arg
         real(kind = 8), intent(out) :: value
         character(100) :: input
         integer(kind = 4) :: read_status

         call ReadArgument(arg, input)
         read(input, *, iostat = read_status) value
         if (read_status /= 0) then
            write(*,*) "invalid real argument: ", arg
            STOP 5
         end if
      end subroutine ReadRealArgument

      subroutine ReadStringArgument(arg, value)
         implicit none
         integer(kind = 4), intent(in) :: arg
         character(*), intent(out) :: value
         character(100) :: input

         call ReadArgument(arg, input)
         value = adjustl(input)
         call Lowercase(value)
      end subroutine ReadStringArgument

      subroutine Lowercase(value)
         implicit none
         character(*), intent(inout) :: value
         integer(kind = 4) :: i, code

         do i = 1, len_trim(value)
            code = iachar(value(i:i))
            if (code >= iachar('A') .AND. code <= iachar('Z')) value(i:i) = achar(code + 32)
         end do
      end subroutine Lowercase

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
