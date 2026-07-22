!------------------------------------------------------------------------------
!
! MODULE: basis
!
! DESCRIPTION:
!> @file basis.F90
!> @brief Module providing utilities related to B-spline basis functions.
!>
!> @details
!> This module groups procedures used to construct, differentiate,
!> query, and evaluate B-spline basis functions and spline fields.
!>
!> The provided functionality includes:
!> - preparation of quadrature-related basis data through \ref BasisData,
!> - evaluation of basis functions and derivatives through
!>   \ref DersBasisFuns,
!> - span lookup in a knot vector through \ref FindSpan,
!> - counting nonempty knot spans through \ref CountSpans,
!> - evaluation of a tensor-product spline field through \ref EvalSpline.
!>
!> These procedures are intended to support numerical integration,
!> assembly, and pointwise post-processing in spline-based solvers.
!
!------------------------------------------------------------------------------
module basis

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes quadrature-related data for nonzero B-spline basis
!> functions on all elements.
!>
!> @details
!> This routine prepares all information required to perform numerical
!> integration over a one-dimensional B-spline basis defined by the
!> knot vector \p U and polynomial order \p p.
!>
!> For each nonempty knot span, the procedure:
!> - identifies the first globally indexed nonzero basis function,
!> - computes the element Jacobian,
!> - maps Gauss points from the reference interval \f$[-1,1]\f$ to the
!>   physical knot span,
!> - evaluates the \f$(p+1)\f$ locally nonzero basis functions and
!>   their derivatives up to order \p d at all quadrature points.
!>
!> The output array \p N is indexed by:
!> - derivative order,
!> - local index of a nonzero basis function,
!> - quadrature-point index,
!> - element index.
!>
!> The basis-function index is local to a given element. For each
!> element there are exactly \f$(p+1)\f$ nonzero basis functions,
!> indexed from \f$0\f$ to \f$p\f$. Conversion to the global basis
!> numbering is obtained through the array \p O according to
!>
!> \f[
!> \text{global\_index} = O(i) + \text{local\_index}.
!> \f]
!>
!> In typical use, the number of quadrature points is chosen as
!> \f$q=p+1\f$, which yields exact Gauss-Legendre integration for
!> polynomials up to degree \f$2p+1\f$.
!
! Input:
! ------
!> @param[in] p
!> Polynomial order of the B-spline basis.
!>
!> @param[in] m
!> Index of the last entry in the knot vector, i.e. number of knot
!> entries minus one.
!>
!> @param[in] U
!> Knot vector defining the spline space.
!>
!> @param[in] d
!> Highest derivative order to be evaluated.
!>
!> @param[in] q
!> Number of Gauss quadrature points.
!>
!> @param[in] r
!> Number of nonempty knot spans, i.e. elements.
!
! Output:
! -------
!> @param[out] O
!> Global indices of the first nonzero basis functions on each element.
!>
!> @param[out] J
!> Element Jacobians for the mapping from the reference interval to
!> physical knot spans.
!>
!> @param[out] W
!> Weights of the Gauss quadrature rule on \f$[-1,1]\f$.
!>
!> @param[out] X
!> Physical coordinates of quadrature points for each element.
!>
!> @param[out] N
!> Values of locally nonzero basis functions and their derivatives at
!> the quadrature points.
!
! Notes:
! ------
!> @note
!> The procedure calls \ref GaussRule to obtain quadrature points and
!> weights and \ref DersBasisFuns to evaluate basis values.
!
!> @warning
!> The argument \p r is assumed to match the number of nonempty knot
!> spans in the knot vector.
!
!---------------------------------------------------------------------------
subroutine BasisData(p, m, U, d, q, r, O, J, W, X, N)
   use gauss, ONLY: GaussRule
   implicit none
!> @brief Polynomial order and last knot index.
   integer(kind=4), intent(in) :: p, m
!> @brief Knot vector.
   real(kind=8), dimension(0:m), intent(in) :: U
!> @brief Highest derivative order, quadrature size, and number of elements.
   integer(kind=4), intent(in) :: d, q, r
!> @brief First global nonzero basis-function index on each element.
   integer(kind=4), dimension(r), intent(out) :: O
!> @brief Jacobians of the element mappings.
   real(kind=8), dimension(r), intent(out) :: J
!> @brief Gauss weights on the reference interval.
   real(kind=8), dimension(q), intent(out) :: W
!> @brief Physical quadrature points for all elements.
   real(kind=8), dimension(q, r), intent(out) :: X
!> @brief Basis-function values and derivatives at quadrature points.
   real(kind=8), dimension(0:d, 0:p, q, r), intent(out) :: N

   integer(kind=4) :: i, iq, ir
   !> @brief Loop counters.
   real(kind=8) :: uu
   !> @brief Physical coordinate of a quadrature point.
   real(kind=8), dimension(q) :: Xg
   !> @brief Gauss points on the reference interval.
   real(kind=8), dimension(0:p, 0:d) :: basis
   !> @brief Local basis values and derivatives at one point.

! Calculates first nonzero basis function for each element
   ir = 1
   do i = p, m - p - 1
      if (U(i) /= U(i + 1)) then
         O(ir) = i - p
         ir = ir + 1
      end if
   end do

   call GaussRule(q, Xg, W)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ir,i,iq,uu,basis) SCHEDULE(STATIC)
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
!$OMP END PARALLEL DO

end subroutine BasisData


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes values of nonzero B-spline basis functions and their
!> derivatives at a specified point.
!>
!> @details
!> This routine evaluates the \f$(p+1)\f$ basis functions that are
!> nonzero at the coordinate \p uu, together with their derivatives up
!> to order \p d. The parameter \p i denotes the knot-span index used
!> to identify the active basis functions.
!>
!> The implementation follows the standard tabular algorithm based on
!> the `ndu` triangular array and the auxiliary array `a` for
!> derivative accumulation.
!>
!> The output array \p ders is indexed by:
!> - local basis-function index from \f$0\f$ to \f$p\f$,
!> - derivative order from \f$0\f$ to \f$d\f$.
!
! Input:
! ------
!> @param[in] i
!> Index of the knot span or, equivalently, the last active basis
!> function associated with the current point.
!>
!> @param[in] uu
!> Coordinate of the evaluation point.
!>
!> @param[in] p
!> Polynomial order of the B-spline basis.
!>
!> @param[in] d
!> Highest derivative order to be computed.
!>
!> @param[in] U
!> Knot vector.
!
! Output:
! -------
!> @param[out] ders
!> Values of active basis functions and their derivatives.
!
! Notes:
! ------
!> @note
!> The routine returns only locally nonzero basis functions. Global
!> indexing must be reconstructed externally from the span index.
!
!> @warning
!> The routine assumes that the point \p uu lies in a valid span
!> compatible with \p i and the knot vector \p U.
!
!---------------------------------------------------------------------------
subroutine DersBasisFuns(i, uu, p, d, U, ders)
   implicit none
!> @brief Span index, polynomial order, and highest derivative order.
   integer(kind=4), intent(in) :: i, p, d
!> @brief Evaluation coordinate.
   real(kind=8), intent(in) :: uu
!> @brief Knot vector.
   real(kind=8), dimension(0:i + p), intent(in) :: U
!> @brief Basis-function values and derivatives.
   real(kind=8), dimension(0:p, 0:d), intent(out) :: ders

!> @brief Integer work variables and loop counters.
   integer(kind=4) :: j, k, r, s1, s2, rk, pk, j1, j2
!> @brief Scalar work variables used by the recurrence.
   real(kind=8) :: saved, temp, der
!> @brief Distances from the evaluation point to surrounding knots.
   real(kind=8), dimension(p) :: left, right
!> @brief Triangular table of basis recursion data.
   real(kind=8), dimension(0:p, 0:p) :: ndu
!> @brief Auxiliary array used in derivative evaluation.
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
!
! DESCRIPTION:
!> @brief Finds the knot span containing a specified point.
!>
!> @details
!> This function returns the index \p span such that the interval
!> \f$[U(\text{span}), U(\text{span}+1))\f$ contains the point \p uu.
!>
!> If the point lies at or beyond the right boundary of the spline
!> domain, the last admissible span index \p n is returned. If the
!> point lies at or before the left admissible boundary, the index
!> \p p is returned.
!>
!> For interior points, the span is found by a binary search over the
!> knot vector.
!
! Input:
! ------
!> @param[in] n
!> Number of basis functions minus one.
!>
!> @param[in] p
!> Polynomial order of the basis.
!>
!> @param[in] uu
!> Coordinate of the query point.
!>
!> @param[in] U
!> Knot vector.
!
! Output:
! -------
!> @return span
!> Index of the knot span containing the point.
!
! Notes:
! ------
!> @note
!> The edge handling is designed for open knot vectors commonly used in
!> spline discretizations.
!
!---------------------------------------------------------------------------
function FindSpan(n, p, uu, U) result(span)
   implicit none
!> @brief Number of basis functions minus one and polynomial order.
   integer(kind=4), intent(in) :: n, p
!> @brief Coordinate of the evaluation point.
   real(kind=8), intent(in) :: uu
!> @brief Knot vector.
   real(kind=8), dimension(0:n + p + 1), intent(in) :: U
!> @brief Index of the containing span.
   integer(kind=4) :: span
!> @brief Bounds used by the binary search.
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
!
! DESCRIPTION:
!> @brief Counts the number of nonempty knot spans in a knot vector.
!>
!> @details
!> This function scans the admissible part of the knot vector and
!> counts all subintervals with nonzero length. Repeated knot values are
!> skipped so that only actual elements are counted.
!>
!> In the present context, the returned value is the number of elements
!> induced by the knot vector.
!
! Input:
! ------
!> @param[in] n
!> Number of basis functions minus one.
!>
!> @param[in] p
!> Polynomial order of the basis.
!>
!> @param[in] U
!> Knot vector.
!
! Output:
! -------
!> @return nelem
!> Number of nonempty knot spans.
!
! Notes:
! ------
!> @note
!> The function treats consecutive repeated knots as belonging to the
!> same span boundary and does not count zero-length intervals.
!
!---------------------------------------------------------------------------
function CountSpans(n, p, U) result(nelem)
   implicit none
   !> @brief Number of basis functions minus one and polynomial order.
   integer(kind=4), intent(in) :: n, p
   !> @brief Knot vector.
   real(kind=8), dimension(0:n + p + 1), intent(in) :: U
   !> @brief Loop counter and returned number of elements.
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
!
! DESCRIPTION:
!> @brief Evaluates a three-dimensional tensor-product spline or its
!> derivative at a specified point.
!>
!> @details
!> This function evaluates a scalar spline field defined by the
!> coefficient array \p coeffs and three one-dimensional B-spline bases
!> associated with the knot vectors \p Ux, \p Uy, and \p Uz.
!>
!> The containing span is determined independently in each coordinate
!> direction by \ref FindSpan. The active one-dimensional basis
!> functions are then evaluated by \ref DersBasisFuns and combined in a
!> tensor-product sum over the local support region.
!>
!> The returned value corresponds to the derivative order \p d in each
!> parametric direction as assembled by the expression
!> \f$bx(ix,d)\,by(iy,d)\,bz(iz,d)\f$.
!
! Input:
! ------
!> @param[in] d
!> Derivative order selected in each coordinate direction.
!>
!> @param[in] Ux
!> Knot vector in the first direction.
!>
!> @param[in] px
!> Polynomial order in the first direction.
!>
!> @param[in] nx
!> Number of basis functions minus one in the first direction.
!>
!> @param[in] nelemx
!> Number of elements in the first direction.
!>
!> @param[in] Uy
!> Knot vector in the second direction.
!>
!> @param[in] py
!> Polynomial order in the second direction.
!>
!> @param[in] ny
!> Number of basis functions minus one in the second direction.
!>
!> @param[in] nelemy
!> Number of elements in the second direction.
!>
!> @param[in] Uz
!> Knot vector in the third direction.
!>
!> @param[in] pz
!> Polynomial order in the third direction.
!>
!> @param[in] nz
!> Number of basis functions minus one in the third direction.
!>
!> @param[in] nelemz
!> Number of elements in the third direction.
!>
!> @param[in] coeffs
!> Tensor-product spline coefficients stored on a three-dimensional
!> control grid.
!>
!> @param[in] x
!> Evaluation coordinate in the first direction.
!>
!> @param[in] y
!> Evaluation coordinate in the second direction.
!>
!> @param[in] z
!> Evaluation coordinate in the third direction.
!
! Output:
! -------
!> @return val
!> Evaluated spline value or derivative value at the point
!> \f$(x,y,z)\f$.
!
! Notes:
! ------
!> @note
!> The arguments \p nelemx, \p nelemy, and \p nelemz are currently
!> passed through the interface but are not used in the implementation.
!
!> @warning
!> The procedure assumes that the query point belongs to the spline
!> domain in all three directions.
!
!---------------------------------------------------------------------------
function EvalSpline(d, &
                     Ux, px, nx, nelemx, &
                     Uy, py, ny, nelemy, &
                     Uz, pz, nz, nelemz, &
                     coeffs, x, y, z) result(val)
   implicit none
!> @brief Derivative order selected in each direction.
   integer(kind=4), intent(in) :: d
!> @brief Basis size, polynomial order, and element count in the first direction.
   integer(kind=4), intent(in) :: nx, px, nelemx
!> @brief Basis size, polynomial order, and element count in the second direction.
   integer(kind=4), intent(in) :: ny, py, nelemy
!> @brief Basis size, polynomial order, and element count in the third direction.
   integer(kind=4), intent(in) :: nz, pz, nelemz
!> @brief Knot vector in the first direction.
   real(kind=8), dimension(0:nx + px + 1), intent(in) :: Ux
!> @brief Knot vector in the second direction.
   real(kind=8), dimension(0:ny + py + 1), intent(in) :: Uy
!> @brief Knot vector in the third direction.
   real(kind=8), dimension(0:nz + pz + 1), intent(in) :: Uz
!> @brief Tensor-product spline coefficients.
   real(kind=8), dimension(0:nx, 0:ny, 0:nz), intent(in) :: coeffs
!> @brief Evaluation coordinates.
   real(kind=8), intent(in) :: x, y, z
!> @brief Evaluated spline value.
   real(kind=8) :: val
!> @brief Basis values or derivatives in the first direction.
   real(kind=8), dimension(0:px, 0:d) :: bx
!> @brief Basis values or derivatives in the second direction.
   real(kind=8), dimension(0:py, 0:d) :: by
!> @brief Basis values or derivatives in the third direction.
   real(kind=8), dimension(0:pz, 0:d) :: bz
!> @brief Temporary tensor-product basis value.
   real(kind=8) :: b
!> @brief Containing spans in the three directions.
   integer(kind=4) :: xspan, yspan, zspan
!> @brief Local loop counters and first active global indices.
   integer(kind=4) :: ix, iy, iz, x0, y0, z0

   xspan = FindSpan(nx, px, x, Ux)
   yspan = FindSpan(ny, py, y, Uy)
   zspan = FindSpan(nz, pz, z, Uz)

   call DersBasisFuns(xspan, x, px, d, ux, bx)
   call DersBasisFuns(yspan, y, py, d, uy, by)
   call DersBasisFuns(zspan, z, pz, d, uz, bz)

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

