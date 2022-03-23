!------------------------------------------------------------------------------
! AGH University of Science and Technology, Krakow, Poland
!------------------------------------------------------------------------------
!
! MODULE: basis
!
!> @author
!> Maciej Wozniak
!
! DESCRIPTION
!> This module contains all functionality associated to basis functions
!
! REVISION HISTORY:
! 21 11 2017 - Initial Version
!
!------------------------------------------------------------------------------

module basis

   implicit none

contains

!---------------------------------------------------------------------------
!> @author Maciej Wozniak
!>
!> @brief
!> Computes all the information about basis B-spline functions needed to perform numerical integration.
!>
!> Computes all the information about basis B-spline functions needed to perform numerical integration.
!> Array of values - \f$ N \f$ - is indexed by:
!>
!>  - derivative
!>  - numer of nonzero function
!>  - index of Gauss quadrature point
!>  - index of element
!>
!> Second index is local - for each element, there are \f$ (p+1) \f$ nonzero
!> basis functions, indexed \f$ 0 \dots p \f$. It does not refer to global index
!> of basis functions on the whole interval.
!>
!> To translate between local and global indexing, use \f$ O \f$ array. \f$ O(i) \f$ is
!> the global index of first nonzero function on i-th element, hence
!>
!>     global_index = O(i) + local_index
!>
!> where local_index is local for i-th element.
!>
!> Value of \f$ q \f$ (number of quadrature points) is taken to be \f$ p+1 \f$, since
!> order of Gaussian quadrature is \f$ 2n - 1 \f$, hence using \f$ q = p+1 \f$ we get
!> exact results for polynomials of degree up to \f$ 2q - 1 = 2p + 1 \f$.
!> This is nice, as we integrate mostly bilinear forms in B-spline
!> basis, consisting of dogree p polynomials.
!
!> Input:
! ------
!> @param[in] p  - polynomial order
!> @param[in] m  - index of the last node in knot vector (number of nodes - 1)
!> @param[in] U  - knot vector
!> @param[in] d  - order of highest derivatives we want to compute
!> @param[in] q  - \f$ q = p+1 \f$ number of Gauss quadrature points
!> @param[in] r  - number of elements (subintervals)
!
!> Output:
! -------
!> @param[out] O - indexes of first nonzero functions on each element
!> @param[out] J - values of the Jacobian of elements
!> @param[out] W - weights of Gauss quadrature points
!> @param[out] X - points of Gauss quadrature
!> @param[out] N - values of \f$ (p+1) \f$ nonzero basis functions and their derivatives at points of Gauss quadrature
! -------------------------------------------------------------------
   subroutine BasisData(p, m, U, d, q, r, O, J, W, X, N)
      use gauss, ONLY: GaussRule
      implicit none
      integer(kind=4), intent(in) :: p, m
      real(kind=8), dimension(0:m), intent(in) :: U
      integer(kind=4), intent(in) :: d, q, r
      integer(kind=4), dimension(r), intent(out) :: O
      real(kind=8), dimension(r), intent(out) :: J
      real(kind=8), dimension(q), intent(out) :: W
      real(kind=8), dimension(q, r), intent(out) :: X
      real(kind=8), dimension(0:d, 0:p, q, r), intent(out) :: N
      integer(kind=4) :: i, iq, ir
      real(kind=8) :: uu
      real(kind=8), dimension(q) :: Xg
      real(kind=8), dimension(0:p, 0:d) :: basis

! Calculates first nonzero basis function for each element
      ir = 1
      do i = p, m - p
         if (U(i) /= U(i + 1)) then
            O(ir) = i - p
            ir = ir + 1
         end if
      end do

      call GaussRule(q, Xg, W)

      do ir = 1, r
         i = O(ir) + p
         J(ir) = (U(i + 1) - U(i))/2.0
         X(:, ir) = (Xg + 1.0)*J(ir) + U(i) ! translate Gauss [-1,1] -> [0,1]
         do iq = 1, q
            uu = X(iq, ir)
            call DersBasisFuns(i, uu, p, d, U, basis)
            N(:, :, iq, ir) = transpose(basis)
         end do
      end do

   end subroutine BasisData

!---------------------------------------------------------------------------
!> @author  Maciej Wozniak
!>
!> @brief
!> Computes values of basis functions and their derivatives at the
!> specified point.
!
! Input:
! ------
!> @param[in] i       - index of (last?) basis function to compute
!> @param[in] uu      - coordinate of the point (argument)
!> @param[in] p       - order of B-splines
!> @param[in] d       - order of highest derivatives to compute
!> @param[in] U       - knot vector
!
! Output:
! ------
!> @param[out] ders   - computed values of basis functions and derivatives
!>
!> Output array is indexed by local index of nonzero basis functions
!> on the element, and order of derivative.
! -------------------------------------------------------------------
   subroutine DersBasisFuns(i, uu, p, d, U, ders)
      implicit none
      integer(kind=4), intent(in) :: i, p, d
      real(kind=8), intent(in) :: uu
      real(kind=8), dimension(0:i + p), intent(in) :: U
      real(kind=8), dimension(0:p, 0:d), intent(out) :: ders
      integer(kind=4) :: j, k, r, s1, s2, rk, pk, j1, j2
      real(kind=8) :: saved, temp, der
      real(kind=8), dimension(p) :: left, right
      real(kind=8), dimension(0:p, 0:p) :: ndu
      real(kind=8), dimension(0:1, 0:p) :: a

      ndu(0, 0) = 1.d0
      do j = 1, p
         left(j) = uu - U(i + 1 - j)
         right(j) = U(i + j) - uu
         saved = 0.0
         do r = 0, j - 1
            ndu(j, r) = right(r + 1) + left(j - r)
            temp = ndu(r, j - 1)/ndu(j, r)
            ndu(r, j) = saved + right(r + 1)*temp
            saved = left(j - r)*temp
         end do
         ndu(j, j) = saved
      end do

      ders(:, 0) = ndu(:, p)

      do r = 0, p
         s1 = 0; s2 = 1; 
         a(0, 0) = 1.0
         do k = 1, d
            der = 0.0
            rk = r - k; pk = p - k; 
            if (r >= k) then
               a(s2, 0) = a(s1, 0)/ndu(pk + 1, rk)
               der = a(s2, 0)*ndu(rk, pk)
            end if
            if (rk > -1) then
               j1 = 1
            else
               j1 = -rk
            end if
            if (r - 1 <= pk) then
               j2 = k - 1
            else
               j2 = p - r
            end if
            do j = j1, j2
               a(s2, j) = (a(s1, j) - a(s1, j - 1))/ndu(pk + 1, rk + j)
               der = der + a(s2, j)*ndu(rk + j, pk)
            end do
            if (r <= pk) then
               a(s2, k) = -a(s1, k - 1)/ndu(pk + 1, r)
               der = der + a(s2, k)*ndu(r, pk)
            end if
            ders(r, k) = der
            j = s1; s1 = s2; s2 = j; 
         end do
      end do
      r = p
      do k = 1, d
         ders(:, k) = ders(:, k)*r
         r = r*(p - k)
      end do

   end subroutine DersBasisFuns

!---------------------------------------------------------------------------
!> @author  Maciej Wozniak
!>
!> @brief
!> Finds the element that contains specified point. It is given as
!> the index of knot u_i such that \f$  [u_i, u_{i+1}] \f$ contains the point.
!> In case the point lies on the edge, index of innermost edge node
!> is returned.
!
! Input:
! ------
!> @param[in] n  - number of functions (control points) minus 1
!> @param[in] p  - order of basis functions
!> @param[in] uu - coordinate of the point
!> @param[in] U  - knot vector
!
! Output:
! ------
!> @return span  - index of element
!
! -------------------------------------------------------------------
   function FindSpan(n, p, uu, U) result(span)
      implicit none
      integer(kind=4), intent(in) :: n, p
      real(kind=8), intent(in) :: uu
      real(kind=8), dimension(0:n + p + 1), intent(in) :: U
      integer(kind=4) :: span
      integer(kind=4) :: low, high

      ! check edge cases
      if (uu >= U(n + 1)) then
         span = n
         return
      end if

      if (uu <= U(p)) then
         span = p
         return
      end if

      ! Binary search for uu
      low = p
      high = n + 1
      span = (low + high)/2

      do while (uu < U(span) .or. uu >= U(span + 1))
         if (uu < U(span)) then
            high = span
         else
            low = span
         end if
         span = (low + high)/2
      end do

   end function FindSpan

!---------------------------------------------------------------------------
!> @author  Maciej Wozniak
!>
!> @brief
!> Calculates number of elements (nonempty subintervals) in the knot vector.
!
! Input:
! ------
!> @param[in] n      - number of functions (control points) minus 1
!> @param[in] p      - order of basis functions
!> @param[in] U      - knot vector
!
! Output:
! ------
!> @param[out] nelem - number of elements
!
! -------------------------------------------------------------------
   function CountSpans(n, p, U) result(nelem)
      implicit none
      integer(kind=4), intent(in) :: n, p
      real(kind=8), dimension(0:n + p + 1), intent(in) :: U
      integer(kind=4) :: i, nelem

      nelem = 0
      i = p
      do while (i <= n)
         ! skip multiple knots
         do while (i < n .and. U(i) == U(i + 1))
            i = i + 1
         end do
#ifdef IPRINT
         write (*, *) 'CountSpans:i,n,U(i),U(i+1)', i, n, U(i), U(i + 1)
#endif
         nelem = nelem + 1
         i = i + 1
      end do

   end function CountSpans

!---------------------------------------------------------------------------
!> @author  Maciej Wozniak
!>
!> @brief
!> Evaluates linear combination of basis spline functions in specified point.
!
! Input:
! ------
!> @param[in] d          - derivative to evaluate
!> @param[in] U_         - knot points
!> @param[in] p_         - degree of splines
!> @param[in] n_         - nubmer of functions minus one
!> @param[in] nelem_     - number of elements
!> @param[in] coeffs     - 3D array of coefficients of basis functions (0-based)
!> @param[in] x, y, z    - point to evaluate at
!
! Output:
! ------
!> @param[out] val        - value in specified point
!
! -------------------------------------------------------------------
   function EvalSpline(d, &
                       Ux, px, nx, nelemx, &
                       Uy, py, ny, nelemy, &
                       Uz, pz, nz, nelemz, &
                       coeffs, x, y, z) result(val)
      implicit none
      integer(kind=4), intent(in) :: d
      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real(kind=8), dimension(0:nx + px + 1), intent(in) :: Ux
      real(kind=8), dimension(0:ny + py + 1), intent(in) :: Uy
      real(kind=8), dimension(0:nz + pz + 1), intent(in) :: Uz
      real(kind=8), dimension(0:nx, 0:ny, 0:nz), intent(in) :: coeffs
      real(kind=8), intent(in) :: x, y, z
      real(kind=8) :: val
      real(kind=8), dimension(0:px, 0:d) :: bx
      real(kind=8), dimension(0:py, 0:d) :: by
      real(kind=8), dimension(0:pz, 0:d) :: bz
      real(kind=8) :: b
      integer(kind=4) :: xspan, yspan, zspan
      integer(kind=4) :: ix, iy, iz, x0, y0, z0

      xspan = FindSpan(nx, px, x, Ux)
      yspan = FindSpan(ny, py, y, Uy)
      zspan = FindSpan(nz, pz, z, Uz)

      call DersBasisFuns(xspan, x, px, 0, ux, bx)
      call DersBasisFuns(yspan, y, py, 0, uy, by)
      call DersBasisFuns(zspan, z, pz, 0, uz, bz)

      x0 = xspan - px
      y0 = yspan - py
      z0 = zspan - pz

      val = 0

      do ix = 0, px
         do iy = 0, py
            do iz = 0, pz
               b = bx(ix, d)*by(iy, d)*bz(iz, d)
               val = val + b*coeffs(x0 + ix, y0 + iy, z0 + iz)
            end do
         end do
      end do

   end function EvalSpline

end module basis

