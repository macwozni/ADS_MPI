!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Input data and source callback for the Eriksson example.
!>
!> @details
!> This module stores the spatial discretization, time-step parameters,
!> and MPI process-grid dimensions used by the Eriksson driver. It also
!> defines the compactly supported initial state and the pointwise forcing
!> callback expected by the current ADS API.
!
!------------------------------------------------------------------------------
module input_data

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
      character(100) :: input
      integer(kind = 4) :: length
      integer(kind = 4) :: status

      ! ./l2 <size> <order> <procx> <procy> <procz>

      call GET_COMMAND_ARGUMENT(1, input, length, status)
      read(input, *) SIZE
      call GET_COMMAND_ARGUMENT(2, input, length, status)
      read(input, *) ORDER
      call GET_COMMAND_ARGUMENT(3, input, length, status)
      read(input, *) steps
      call GET_COMMAND_ARGUMENT(4, input, length, status)
      read(input, *) Dt
      call GET_COMMAND_ARGUMENT(5, input, length, status)
      read(input, *) procx
      call GET_COMMAND_ARGUMENT(6, input, length, status)
      read(input, *) procy
      call GET_COMMAND_ARGUMENT(7, input, length, status)
      read(input, *) procz

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


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the pointwise forcing for the Eriksson example.
!>
!> @details
!> At the initial pseudo-step the callback injects the initial state. For
!> later time steps it returns zero, so the subsequent evolution is driven
!> by the ADS operator and stored state.
!
! Input:
! ------
!> @param[in] un
!> Previous solution value at the quadrature point.
!>
!> @param[in] du
!> Previous solution gradient at the quadrature point.
!>
!> @param[in] X
!> Physical coordinates of the quadrature point.
!
! Output:
! -------
!> @return ret
!> Pointwise forcing value.
!
!---------------------------------------------------------------------------
   function forcing(un, du, X) result(ret)
      implicit none
      real(kind = 8), intent(in) :: un
      real(kind = 8), intent(in), dimension(3) :: du
      real(kind = 8), intent(in), dimension(3) :: X
      real(kind = 8) :: ret

      if (t > 0.d0) then
         ret = 0.d0
      else
         ret = initial_state(X(1), X(2), X(3))
      endif

   end function forcing


end module input_data
