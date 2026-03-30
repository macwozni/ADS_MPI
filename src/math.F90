!------------------------------------------------------------------------------
!
! MODULE: math
!
! DESCRIPTION:
!> @file math.F90
!> @brief Module providing elementary mathematical helper routines and
!> constants used across the project.
!>
!> @details
!> This module groups lightweight numerical utilities that are reused by
!> higher-level components of the code base.
!>
!> The provided functionality includes:
!> - the mathematical constant \ref PI,
!> - linear interpolation through \ref lerp,
!> - a smooth radial falloff function through \ref falloff,
!> - one-dimensional bump functions through \ref bump and \ref bump01,
!> - a three-dimensional compact-support profile through \ref bump3d,
!> - Euclidean norm evaluation through \ref norm2.
!>
!> These routines are used primarily by plotting, test-function
!> construction, and auxiliary numerical workflows.
!
!------------------------------------------------------------------------------
module math

   implicit none

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Numerical approximation of the constant \f$ \pi \f$ in double
!> precision.
!>
!> @details
!> The value is constructed from the identity
!> \f$ \pi = 4 \arctan(1) \f$
!> using the intrinsic double-precision arctangent function.
!
!---------------------------------------------------------------------------
   real(kind=8), parameter :: PI = 4.d0*datan(1.d0)

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Performs linear interpolation between two scalar values.
!>
!> @details
!> This function returns the convex combination
!>
!> \f[
!> (1-t)x + ty,
!> \f]
!>
!> which corresponds to the value obtained by linearly interpolating
!> between \p x and \p y with interpolation parameter \p t.
!>
!> In particular:
!> - \f$ t = 0 \f$ returns \p x,
!> - \f$ t = 1 \f$ returns \p y,
!> - intermediate values of \p t produce points on the segment joining
!>   \p x and \p y.
!
! Input:
! ------
!> @param[in] t
!> Interpolation parameter, typically taken from the interval
!> \f$[0,1]\f$.
!>
!> @param[in] x
!> Left endpoint value.
!>
!> @param[in] y
!> Right endpoint value.
!
! Output:
! -------
!> @return val
!> Interpolated scalar value.
!
!---------------------------------------------------------------------------
function lerp(t, x, y) result(val)
   implicit none
!> @brief Interpolation parameter and endpoint values.
   real(kind=8), intent(in) :: t, x, y
!> @brief Interpolated value.
   real(kind=8) :: val

   val = (1.d0 - t)*x + t*y

end function lerp

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates a smooth one-dimensional decay profile with compact
!> transition region.
!>
!> @details
!> This function defines a piecewise profile that is constant and equal
!> to one to the left of \p r, constant and equal to zero to the right of
!> \p Rr, and smoothly decays on the interval \f$[r,Rr]\f$.
!>
!> The transition profile is based on the quartic polynomial
!>
!> \f[
!> \left((h-1)(h+1)\right)^2,
!> \qquad
!> h = \frac{t-r}{Rr-r},
!> \f]
!>
!> which yields a function with vanishing first derivative at the
!> transition endpoints.
!
! Input:
! ------
!> @param[in] r
!> Left boundary of the decay interval.
!>
!> @param[in] Rr
!> Right boundary of the decay interval.
!>
!> @param[in] t
!> Evaluation point.
!
! Output:
! -------
!> @return fval
!> Value of the falloff function at \p t.
!
! Notes:
! ------
!> @note
!> For \f$ t < r \f$, the function returns \f$1\f$.
!>
!> @note
!> For \f$ t > Rr \f$, the function returns \f$0\f$.
!
!> @warning
!> The implementation assumes \p Rr and \p r are distinct. If
!> \f$Rr=r\f$, the transition formula is singular.
!
!---------------------------------------------------------------------------
function falloff(r, Rr, t) result(fval)
   implicit none
!> @brief Transition interval endpoints and evaluation point.
   real(kind=8), intent(in) :: r, Rr, t
!> @brief Normalized transition coordinate and function value.
   real(kind=8) :: h, fval

   if (t < r) then
      fval = 1.d0
   else if (t > Rr) then
      fval = 0.d0
   else
      h = (t - r)/(Rr - r)
      fval = ((h - 1.d0)*(h + 1.d0))**2.d0
   end if

end function falloff

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates a one-dimensional \f$C^1\f$ bump function on the
!> interval \f$[-1,1]\f$.
!>
!> @details
!> This function is equal to zero outside the interval \f$[-1,1]\f$ and
!> inside that interval is given by the polynomial
!>
!> \f[
!> (x+1)^2(x-1)^2.
!> \f]
!>
!> The function vanishes together with its first derivative at the
!> endpoints \f$-1\f$ and \f$1\f$, and attains its maximum value
!> \f$1\f$ at \f$x=0\f$.
!
! Input:
! ------
!> @param[in] x
!> Evaluation point.
!
! Output:
! -------
!> @return val
!> Value of the bump function at \p x.
!
!---------------------------------------------------------------------------
function bump(x) result(val)
   implicit none
!> @brief Evaluation point.
   real(kind=8), intent(in) :: x
!> @brief Function value.
   real(kind=8) :: val

   if (x > -1.d0 .and. x < 1.d0) then
      val = (x + 1.d0)**2.d0*(x - 1.d0)**2.d0
   else
      val = 0.d0
   end if

end function bump

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates a rescaled \f$C^1\f$ bump function on the interval
!> \f$[0,1]\f$.
!>
!> @details
!> This function is a shifted and rescaled version of \ref bump. It maps
!> the interval \f$[0,1]\f$ to \f$[-1,1]\f$ via the transformation
!>
!> \f[
!> \xi = 2(x-0.5),
!> \f]
!>
!> and then evaluates \ref bump at \f$\xi\f$.
!>
!> The function vanishes together with its first derivative at
!> \f$0\f$ and \f$1\f$, and attains its maximum value \f$1\f$ at
!> \f$x=0.5\f$.
!
! Input:
! ------
!> @param[in] x
!> Evaluation point.
!
! Output:
! -------
!> @return val
!> Value of the rescaled bump function at \p x.
!
!---------------------------------------------------------------------------
function bump01(x) result(val)
   implicit none
!> @brief Evaluation point.
   real(kind=8), intent(in) :: x
!> @brief Function value.
   real(kind=8) :: val

   val = bump(2.d0*(x - 0.5d0))

end function bump01

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates a three-dimensional radially symmetric bump-like
!> profile centered at \f$(0.5,0.5,0.5)\f$.
!>
!> @details
!> This function computes the Euclidean distance of the point
!> \f$(x,y,z)\f$ from the center of the unit cube and applies the
!> one-dimensional falloff profile \ref falloff to that distance.
!>
!> More precisely, it evaluates
!>
!> \f[
!> t = \left\| (x,y,z) - (0.5,0.5,0.5) \right\|_2
!> \f]
!>
!> and returns a smooth radial profile transitioning from one to zero
!> between radii \f$r/2\f$ and \f$Rr/2\f$.
!>
!> This yields a compact-support-like three-dimensional shape convenient
!> for test forcing terms and localized profiles.
!
! Input:
! ------
!> @param[in] r
!> Inner diameter-like scale controlling the constant core.
!>
!> @param[in] Rr
!> Outer diameter-like scale controlling the end of the decay region.
!>
!> @param[in] x
!> First coordinate of the evaluation point.
!>
!> @param[in] y
!> Second coordinate of the evaluation point.
!>
!> @param[in] z
!> Third coordinate of the evaluation point.
!
! Output:
! -------
!> @return val
!> Value of the three-dimensional radial profile.
!
! Notes:
! ------
!> @note
!> The radial distance is evaluated with the local helper function
!> \ref norm2.
!
!---------------------------------------------------------------------------
function bump3d(r, Rr, x, y, z) result(val)
   implicit none
!> @brief Profile scales and Cartesian evaluation coordinates.
   real(kind=8), intent(in) :: r, Rr, x, y, z
!> @brief Function value and radial distance from the cube center.
   real(kind=8) :: val, t

   t = norm2((/x, y, z/) - 0.5d0)
   val = falloff(r/2.d0, Rr/2.d0, t)

end function bump3d

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes the Euclidean norm of a real vector.
!>
!> @details
!> This function evaluates the standard \f$L_2\f$ norm
!>
!> \f[
!> \|x\|_2 = \sqrt{x \cdot x},
!> \f]
!>
!> using the intrinsic procedures `dot_product` and `sqrt`.
!>
!> The routine serves as a project-local implementation corresponding to
!> the Fortran 2008 intrinsic function with the same name.
!
! Input:
! ------
!> @param[in] x
!> Input vector.
!
! Output:
! -------
!> @return norm2
!> Euclidean norm of the vector \p x.
!
!---------------------------------------------------------------------------
function norm2(x)
   implicit none
   intrinsic :: dot_product, sqrt
!> @brief Input vector.
   real(kind=8), dimension(:), intent(in) :: x
!> @brief Euclidean norm of the input vector.
   real(kind=8) :: norm2

   norm2 = sqrt(dot_product(x, x))
end function norm2

end module math