!------------------------------------------------------------------------------
!
! MODULE: knot_vector
!
! DESCRIPTION:
!> @file knot_vector.F90
!> @brief Module providing procedures for creation and modification of
!> spline knot vectors.
!>
!> @details
!> This module groups utilities used to allocate, fill, and post-process
!> one-dimensional knot vectors defined on the interval \f$[0,1]\f$.
!>
!> The provided functionality includes:
!> - preparation of an open knot vector through \ref PrepareKnot,
!> - direct filling of an open knot vector through \ref FillOpenKnot,
!> - counting nonempty knot spans through \ref CountSpans,
!> - construction of knot vectors with repeated internal knots through
!>   \ref repeatedKnot.
!>
!> The routines are intended for use in spline-space setup and mesh-like
!> parametrization stages preceding assembly or evaluation procedures.
!
!------------------------------------------------------------------------------
module knot_vector

   implicit none

contains


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Allocates and fills an open knot vector on the interval
!> \f$[0,1]\f$.
!>
!> @details
!> This routine allocates the knot vector \p U for a spline space with
!> \f$n+1\f$ basis functions and polynomial degree \p p, fills it with
!> a standard open uniform knot sequence, and computes the number of
!> nonempty knot spans.
!>
!> In this convention, the number of subintervals is
!> \f$N = n - p + 1\f$. The endpoints \f$0\f$ and \f$1\f$ are repeated
!> \f$p+1\f$ times, which yields an open knot vector.
!>
!> The actual filling of the vector is delegated to \ref FillOpenKnot,
!> while the number of elements is obtained through \ref CountSpans.
!
! Input:
! ------
!> @param[in] n
!> Number of basis functions minus one.
!>
!> @param[in] p
!> Polynomial degree of the spline basis.
!
! Output:
! -------
!> @param[out] U
!> Allocated and filled knot vector.
!>
!> @param[out] nelem
!> Number of nonempty knot spans in the resulting knot vector.
!
! Notes:
! ------
!> @note
!> The allocated size of \p U is \f$n + p + 2\f$.
!
!---------------------------------------------------------------------------
subroutine PrepareKnot(n, p, U, nelem)
   implicit none
!> @brief Number of basis functions minus one and polynomial degree.
   integer(kind=4), intent(in) :: n, p
!> @brief Allocated knot vector.
   real(kind=8), allocatable, dimension(:), intent(out) :: U
!> @brief Number of nonempty knot spans.
   integer(kind=4), intent(out) :: nelem

   allocate (U(n + p + 2))
   call FillOpenKnot(n, p, U)
   nelem = CountSpans(n, p, U)

#ifdef IINFO
   write (*, *) 'n,p,nelem', n, p, nelem
   write (*, *) 'U', U
#endif

end subroutine PrepareKnot


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Fills an existing array with an open knot vector on the
!> interval \f$[0,1]\f$.
!>
!> @details
!> This routine writes a standard open uniform knot vector into the
!> output array \p U. The first and last knot values are repeated
!> \f$p+1\f$ times, and the internal knots are distributed uniformly
!> between \f$0\f$ and \f$1\f$.
!>
!> The number of subintervals induced by the knot vector is
!> \f$N = n - p + 1\f$.
!
! Input:
! ------
!> @param[in] n
!> Number of basis functions minus one.
!>
!> @param[in] p
!> Polynomial degree of the spline basis.
!
! Output:
! -------
!> @param[out] U
!> Knot vector filled with the open uniform distribution.
!
! Notes:
! ------
!> @note
!> The output array \p U must already be allocated with length
!> \f$n + p + 2\f$.
!
!---------------------------------------------------------------------------
subroutine FillOpenKnot(n, p, U)
   implicit none
!> @brief Number of basis functions minus one and polynomial degree.
   integer(kind=4), intent(in) :: n, p
!> @brief Output knot vector.
   real(kind=8), dimension(1:n + p + 2), intent(out) :: U
!> @brief Loop counter.
   integer(kind=4) :: i

   U(1:p + 1) = 0.d0
   U(n + 2:n + p + 2) = 1.d0

   do i = p + 2, n + 1
      U(i) = real(i - p - 1)/real(n - p + 1)
   end do

end subroutine FillOpenKnot


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Counts the number of nonempty knot spans in a knot vector.
!>
!> @details
!> This function scans the admissible part of the knot vector and
!> counts all intervals of positive length. Consecutive repeated knot
!> values are skipped, so only actual nonzero-length spans contribute
!> to the result.
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
!> Polynomial order of the basis functions.
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
!> Zero-length spans created by repeated knots are not counted.
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
!> @brief Constructs an open knot vector with repeated internal knots
!> at block interfaces.
!>
!> @details
!> This routine starts from the standard open uniform knot vector on
!> \f$[0,1]\f$ and repeats internal knots located at block interfaces,
!> i.e. after every \p iblock elements.
!>
!> The continuity at a block interface is prescribed by
!> \p c_continuity. For a spline basis of degree \p p and internal knot
!> multiplicity \f$m\f$, the continuity is \f$C^{p-m}\f$. Therefore, the
!> required multiplicity is \f$m = p - c\_continuity\f$.
!>
!> On exit, both the knot vector \p U and the value of \p n are updated
!> so that they are consistent with the inserted knot repetitions.
!
! Input:
! ------
!> @param[in,out] n
!> Number of basis functions minus one. On input, it defines the base
!> open uniform space. On output, it is updated to reflect the enlarged
!> space after knot repetitions have been inserted.
!>
!> @param[in] p
!> Polynomial degree of the spline basis.
!>
!> @param[in] iblock
!> Block length measured in elements. A repeated internal knot is
!> inserted after every \p iblock elements.
!>
!> @param[in] c_continuity
!> Desired continuity at block interfaces.
!
! Output:
! -------
!> @param[out] U
!> Allocated knot vector with repeated internal knots at block
!> interfaces.
!>
!> @param[out] nelem
!> Number of nonempty knot spans in the base open uniform
!> discretization.
!
! Notes:
! ------
!> @note
!> For \f$c\_continuity = p - 1\f$, the routine reproduces the standard
!> open knot vector.
!>
!> @note
!> For \f$c\_continuity = 0\f$, the routine enforces \f$C^0\f$
!> continuity at block interfaces.
!>
!> @warning
!> For \f$c\_continuity = -1\f$, the resulting spline space is
!> discontinuous across block interfaces.
!>
!> @warning
!> The routine assumes an open uniform parametrization on \f$[0,1]\f$
!> and that \p iblock is compatible with the element partition.
!
!---------------------------------------------------------------------------
subroutine repeatedKnot(n, p, iblock, c_continuity, U, nelem)

   implicit none

!> @brief Number of basis functions minus one, updated on exit.
   integer(kind=4), intent(inout) :: n
!> @brief Polynomial degree, block size, and target continuity.
   integer(kind=4), intent(in) :: p, iblock, c_continuity
!> @brief Allocated knot vector with repeated internal knots.
   real(kind=8), allocatable, dimension(:), intent(out) :: U
!> @brief Number of nonempty spans of the base open knot vector.
   integer(kind=4), intent(out) :: nelem

   integer(kind=4) :: n_old, mult_target, extra_mult, nbreaks, n_new
!> @brief Auxiliary dimensions, multiplicities, loop counters, and indices.
   integer(kind=4) :: e, k, pos
!> @brief Current internal knot value.
   real(kind=8) :: xi

   n_old = n

   nelem = n_old - p + 1

   ! continuity = p - multiplicity  =>  multiplicity = p - continuity
   mult_target = p - c_continuity

   ! additional copies relative to the standard open knot vector
   extra_mult = mult_target - 1

   ! number of block interfaces: iblock, 2*iblock, ...
   nbreaks = (nelem - 1) / iblock

   ! updated number of basis functions minus one
   n_new = n_old + nbreaks*extra_mult

   allocate (U(n_new + p + 2))

   pos = 1

   ! left boundary: p+1 zeros
   do k = 1, p + 1
      U(pos) = 0.0d0
      pos = pos + 1
   end do

   ! internal knots
   do e = 1, nelem - 1

      xi = real(e, kind=8) / real(nelem, kind=8)

      if (mod(e, iblock) == 0) then
         do k = 1, mult_target
            U(pos) = xi
            pos = pos + 1
         end do
      else
         U(pos) = xi
         pos = pos + 1
      end if

   end do

   ! right boundary: p+1 ones
   do k = 1, p + 1
      U(pos) = 1.0d0
      pos = pos + 1
   end do

   n = n_new

end subroutine repeatedKnot

end module knot_vector
