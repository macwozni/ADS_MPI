!------------------------------------------------------------------------------
! AGH University of Science and Technology, Krakow, Poland
!------------------------------------------------------------------------------
!
! MODULE: math
!
!> @author
!> Maciej Wozniak
!
! DESCRIPTION:
!> This module contains mathematical helpers.
!
! REVISION HISTORY:
! 21 11 2017 - Initial Version
! 
!------------------------------------------------------------------------------


module math

implicit none

!> Definition of \f$ \Pi \f$
real (kind = 8), parameter :: PI = 4.d0 * datan(1.d0)


contains

!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> Linear interpolation between two values.
!
! Input:
! ------
!> @param[in] t      - interpolation parameter \f$ t \in [0,1] \f$, point between x and y
!> @param[in] x, y   - points to interpolate between \f$ (t=0 \rightarrow x, t=1 \rightarrow y) \f$
!
! Output:
! -------
!> @return val       - interpolated value
! -------------------------------------------------------------------
function lerp(t, x, y) result (val)
   implicit none
   real (kind = 8), intent(in) :: t, x, y
   real (kind = 8) :: val

   val = (1.d0 - t) * x + t * y

end function lerp


!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> Smoothly decaying function. \f$ C^1 \f$ on \f$ [r, Rr] \f$. Function and its 1-st derivative 
!> vanish at the endpoints.
!> \f$ t \in [-\infty, r] \rightarrow 1 \f$,
!> \f$ t \in [r, R] \rightarrow [1,0] \f$,
!> \f$ t \in [R, \infty] \rightarrow 0 \f$
!
! Input:
! ------
!> @param[in] r   - left boundary of decaing part
!> @param[in] Rr  - right boundary of decaing part
!> @param[in] t   - point at which value should be computed
!
! Output:
! -------
!> @return fval   - function value at point t
! -------------------------------------------------------------------
function falloff(r, Rr, t) result (fval)
   implicit none
   real (kind = 8), intent(in) :: r, Rr, t
   real (kind = 8) :: h, fval

   if (t < r) then
      fval = 1.d0
   else if (t > Rr) then
      fval = 0.d0
   else
      h = (t - r) / (Rr - r)
      fval = ((h - 1.d0) * (h + 1.d0)) ** 2.d0
   endif

end function falloff



!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> \f$ C^1 \f$ bump function on \f$ [-1, 1] \f$. Function and its 1-st derivative 
!> vanish at the endpoints, maximum is attained for \f$ x = 0 \f$ and its
!> value is \f$ 1 \f$.
!
! Input:
! ------
!> @param[in] x - point at which value should be computed
!
! Output:
! -------
!> @return val  - function value at point t
! -------------------------------------------------------------------
function bump(x) result (val)
   implicit none
   real (kind = 8), intent(in) :: x
   real (kind = 8) :: val

   if (x > -1.d0 .and. x < 1.d0) then
      val = (x + 1.d0)**2.d0 * (x - 1.d0)**2.d0
   else
      val = 0.d0
   endif

end function bump


!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> \f$ C^1 \f$ bump function on \f$ [0, 1] \f$. Function and its 1st derivative 
!> vanish at the endpoints, maximum is attained for \f$ x = 0.5 \f$ and its
!> value is \f$ 1 \f$.
!> It's a rescaled version of bump subroutine.
!
! Input:
! ------
!> @param[in] x  - point at which value should be computed
!
! Output:
! -------
!> @return val   - function value at point t
! -------------------------------------------------------------------
function bump01(x) result (val)
   implicit none
   real (kind = 8), intent(in) :: x
   real (kind = 8) :: val

      val = bump(2.d0 * (x - 0.5d0))

end function bump01


!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> \f$ C^1 \f$ bump function on \f$ [r, Rr] \times [r,Rr] \times [r,Rr] \f$.
!> Function and its 1st derivative vanish at the endpoints, maximum is attained for \f$ || x,y,x || <= r \f$ and its
!> value is \f$ 1 \f$.
!> It's a 3D version of bump subroutine.
!> \f$ || x,y,x || \in [0, r] \rightarrow 1 \f$,
!> \f$ || x,y,x || \in [r, R] \rightarrow [1,0] \f$,
!> \f$ || x,y,x || \in [R, \infty] \rightarrow 0 \f$.
!
! Input:
! ------
!> @param[in] r    - left boundary of decaing part
!> @param[in] Rr   - right boundary of decaing part
!> @param[in] x    - point at which value should be computed
!> @param[in] y    - point at which value should be computed
!> @param[in] z    - point at which value should be computed
!
! Output:
! -------
!> @return val     - function value at point t
!--------------------------------------------------------------------------- 
function bump3d(r, Rr, x, y, z) result (val)
   implicit none
   real (kind = 8), intent(in) :: r, Rr, x, y, z
   real (kind = 8) :: val, t

   t = norm2((/x, y, z/) - 0.5d0)
   val = falloff(r/2.d0, Rr/2.d0, t)

end function bump3d

!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> Calculates the Euclidean vector norm (\f$ L_2 \f$ norm) of of vector x. 
!> Implementation of Fortran2008 function with the same name.
!
! Input:
! ------
!> @param[in] x   - vector
!
! Output:
! -------
!> @return norm2  - norm of vector x
!--------------------------------------------------------------------------- 
function norm2(x)
   implicit none
   intrinsic :: dot_product, sqrt
   real (kind = 8), intent(in) :: x(:)
   real (kind = 8) :: norm2
   norm2 = sqrt(dot_product(x,x))
end function norm2
   
end module math

