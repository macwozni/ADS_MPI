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
!> @brief Constructs a knot vector with repeated internal knots grouped
!> by blocks.
!>
!> @details
!> This routine first builds a reference open knot vector, determines
!> the number of nonempty spans, and then reallocates the knot vector so
!> that selected internal knots are repeated according to the block size
!> \p iblock and polynomial degree \p p.
!>
!> The procedure currently provides specialized logic for the cases
!> \f$p=2\f$, \f$p=3\f$, and \f$p=4\f$. For each supported degree, the
!> routine modifies both the size of the knot vector and the value of
!> \p n so that the resulting spline space reflects the introduced knot
!> repetitions.
!>
!> This type of construction is intended for block-wise repeated-knot
!> layouts, for example in reduced-continuity spline spaces.
!
! Input:
! ------
!> @param[in,out] n
!> Number of basis functions minus one. On exit, this value is updated
!> to reflect the modified knot-vector structure.
!>
!> @param[in] p
!> Polynomial degree of the spline basis.
!>
!> @param[in] iblock
!> Block length controlling the placement of repeated internal knots.
!
! Output:
! -------
!> @param[out] U
!> Allocated knot vector with repeated internal knots.
!>
!> @param[out] nelem
!> Number of nonempty spans in the initially constructed open knot
!> vector used as a reference for the repeated-knot layout.
!
! Notes:
! ------
!> @note
!> The implementation currently handles only degrees \f$p=2\f$,
!> \f$p=3\f$, and \f$p=4\f$.
!
!> @warning
!> No fallback branch is provided for unsupported polynomial degrees.
!>
!> @warning
!> The routine assumes that \p iblock is compatible with the chosen
!> degree and element count.
!
!---------------------------------------------------------------------------
subroutine repeatedKnot(n, p, iblock, U, nelem)
   implicit none
   integer(kind=4), intent(inout) :: n
!> @brief Number of basis functions minus one, updated on exit.
   integer(kind=4), intent(in) :: p, iblock
!> @brief Polynomial degree and block size controlling repetitions.
   real(kind=8), allocatable, dimension(:), intent(out) :: U
!> @brief Allocated knot vector with repeated internal knots.
   integer(kind=4), intent(out) :: nelem
!> @brief Number of nonempty spans of the initial open knot vector.
   integer(kind=4) :: i, j, k, l
!> @brief Loop counters and auxiliary indexing variables.

   !---------------------------------------------------------------------------
   if (p .EQ. 2) then

      allocate (U(n + p + 2)) !knot vector
      U(1:p + 1) = 0.d0
      U(n + 2:n + p + 2) = 1.d0
      do i = p + 2, n + 1
         U(i) = real(i - p - 1)/real(n - 1)
      end do

      nelem = CountSpans(n, p, U)

      deallocate (U)
      allocate (U(n + p + 2 + (nelem/iblock)*(p - 1) - 1))
      n = n + (nelem/iblock)*(p - 1) - 1
      U(1:p + 1) = 0.d0
      i = p + 2; l = 1
      do j = 1, nelem/IBLOCK
         do k = 0, IBLOCK - p
            U(i + k) = real(l + k)/real(nelem)
         end do
         if (j .lt. nelem/IBLOCK) then
            do k = IBLOCK - p + 1, IBLOCK
               U(i + k) = real(l + IBLOCK - p + 1)/real(nelem)
            end do
            i = i + IBLOCK + 1; l = l + IBLOCK - p + 2
         else
            i = i + IBLOCK + 1 - p
         end if
      end do
      U(i:i + p) = 1.d0
   end if
   !-------------------------------------------------------------------------
   if (p .EQ. 3) then
      allocate (U(n + p + 2)) !knot vector
      U(1:p + 1) = 0.d0
      U(n + 2:n + p + 2) = 1.d0
      do i = p + 2, n + 1
         U(i) = real(i - p - 1)/real(n - p + 1)
      end do
      !
      nelem = CountSpans(n, p, U)

      deallocate (U)
      allocate (U(n + p + 2 + (nelem/iblock)*(p - 1) - 2))
      n = n + (nelem/iblock)*(p - 1) - 2
      U(1:p + 1) = 0.d0
      i = p + 2; l = 1
      do j = 1, nelem/IBLOCK
         do k = 0, IBLOCK - p + 1
            U(i + k) = real(l + k)/real(nelem)
         end do
         if (j .lt. nelem/IBLOCK) then
            do k = IBLOCK - p + 1, IBLOCK
               U(i + k + 1) = real(l + IBLOCK - p + 2)/real(nelem)
            end do
            i = i + IBLOCK + 1 - p + 4; l = l + IBLOCK - p + 3
         else
            i = i + IBLOCK + 1 - p + 1
         end if
      end do
      U(i:i + p) = 1.d0
   end if
   !----------------------------------------------------------------------
   if (p .EQ. 4) then
      allocate (U(n + p + 2)) !knot vector
      U(1:p + 1) = 0.d0
      U(n + 2:n + p + 2) = 1.d0
      do i = p + 2, n + 1
         U(i) = real(i - p - 1)/real(n - p + 1)
      end do

      nelem = CountSpans(n, p, U)

      deallocate (U)
      allocate (U(n + p + 2 + (nelem/iblock)*(p - 1) - 3))
      n = n + (nelem/iblock)*(p - 1) - 3
      U(1:p + 1) = 0.d0
      i = p + 2; l = 1
      do j = 1, nelem/IBLOCK
         do k = 0, IBLOCK - p + 2
            U(i + k) = real(l + k)/real(nelem)
         end do
         if (j .lt. nelem/IBLOCK) then
            do k = IBLOCK - p + 1, IBLOCK
               U(i + k + 2) = real(l + IBLOCK - p + 3)/real(nelem)
            end do
            i = i + IBLOCK + 1 - p + 6; l = l + IBLOCK - p + 4
         else
            i = i + IBLOCK + 1 - p + 2
         end if
      end do
   end if

   !<- rIGA
end subroutine repeatedKnot


!> @brief Builds an open uniform knot vector with repeated internal knots
!!        at block interfaces.
!!
!! This routine constructs a knot vector on the interval [0,1] for a
!! B-spline space of degree @p p. Starting from the standard open uniform
!! knot vector, it increases the multiplicity of internal knots located at
!! block interfaces, i.e. after every @p iblock elements.
!!
!! The target continuity at each block interface is given by
!! @p c_continuity. For a B-spline of degree @p p and an internal knot
!! multiplicity @f$m@f$, the continuity is
!! @f[
!!   C^{p-m}.
!! @f]
!! Therefore the multiplicity required to obtain continuity
!! @p c_continuity is
!! @f[
!!   m = p - c\_continuity.
!! @f]
!!
!! In particular:
!! - `c_continuity = p-1` gives the standard open knot vector,
!! - `c_continuity = 0` gives `C^0` continuity at block interfaces,
!! - `c_continuity = -1` gives a fully discontinuous interface (`C^{-1}`).
!!
!! On output, @p n is updated consistently with the new knot vector size.
!!
!! @param[in,out] n
!!   Number of basis functions minus one. On input, it defines the base
!!   open uniform space. On output, it is updated to reflect the enlarged
!!   space after knot repetitions have been inserted.
!!
!! @param[in] p
!!   Polynomial degree of the B-spline basis.
!!
!! @param[in] iblock
!!   Block length measured in elements. A repeated knot is inserted after
!!   every `iblock` elements, then after `2*iblock`, etc.
!!
!! @param[in] c_continuity
!!   Desired continuity at block interfaces.
!!
!! @param[out] U
!!   Allocated knot vector with repeated knots at block interfaces.
!!
!! @param[out] nelem
!!   Number of nonzero knot spans of the base open uniform discretization,
!!   equal to `n_old - p + 1`, where `n_old` is the input value of @p n.
!!
!! @note
!!   This routine assumes an open uniform parametrization on [0,1].
!!
!! @warning
!!   For `c_continuity = -1`, the resulting spline space is discontinuous
!!   across block interfaces. Any downstream assembly or connectivity logic
!!   must be compatible with this.
!!
subroutine repeatedKnot2(n, p, iblock, c_continuity, U, nelem)

   implicit none

   integer(kind=4), intent(inout)         :: n
   integer(kind=4), intent(in)            :: p
   integer(kind=4), intent(in)            :: iblock
   integer(kind=4), intent(in)            :: c_continuity
   real(kind=8), allocatable, intent(out) :: U(:)
   integer(kind=4), intent(out)           :: nelem

   integer(kind=4) :: n_old
   integer(kind=4) :: mult_target
   integer(kind=4) :: extra_mult
   integer(kind=4) :: nbreaks
   integer(kind=4) :: n_new
   integer(kind=4) :: e
   integer(kind=4) :: k
   integer(kind=4) :: pos
   real(kind=8)    :: xi

   n_old = n

   nelem = n_old - p + 1

   ! continuity = p - multiplicity  =>  multiplicity = p - continuity
   mult_target = p - c_continuity

   ! additional copies relative to the standard open knot vector
   extra_mult = mult_target - 1

   ! number of block interfaces: iblock, 2*iblock, ...
   nbreaks = (nelem - 1) / iblock

   ! updated number of basis functions minus one
   n_new = n_old + nbreaks * extra_mult

   allocate(U(n_new + p + 2))

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

end subroutine repeatedKnot2

end module knot_vector
