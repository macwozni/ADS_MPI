!------------------------------------------------------------------------------
!
! MODULE: knot_vector
!
! DESCRIPTION:
!> This module contains all functionality associated to creating and filling knot vectors.
!
!------------------------------------------------------------------------------

module knot_vector

   implicit none

contains

!---------------------------------------------------------------------------
!> @brief
!> Allocates and fills the knot vector on \f$  [0, 1] \f$.
!> Number of subintervals is \f$ N = n-p+1 \f$.
!> \f$ 0 \f$ and \f$ 1 \f$ are repeated \f$ (p+1) \f$ times.
!
! Input:
! ------
!> @param[in] n      - number of functions on the knot minus one
!> @param[in] p      - degree of polynomial
!
! Output:
! -------
!> @param[out] U     - array to fill with points
!> @param[out] nelem - number of elements
!---------------------------------------------------------------------------
   subroutine PrepareKnot(n, p, U, nelem)
      implicit none
      integer(kind=4), intent(in) :: n, p
      real(kind=8), allocatable, dimension(:), intent(out) :: U
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
!> @brief
!> Fills knot vector on \f$ [0, 1] \f$.
!> Number of subintervals is \f$ N = n-p+1 \f$.
!> \f$ 0 \f$ and \f$ 1 \f$ are repeated \f$ (p+1) \f$ times.
!
! Input:
! ------
!> @param[in] n  - number of functions on the knot minus one
!> @param[in] p  - degree of polynomial
!
! Output:
! -------
!> @param[out] U  - array to fill with points
!---------------------------------------------------------------------------
   subroutine FillOpenKnot(n, p, U)
      implicit none
      integer(kind=4), intent(in) :: n, p
      real(kind=8), dimension(1:n + p + 2), intent(out) :: U
      integer(kind=4) :: i

      U(1:p + 1) = 0.d0
      U(n + 2:n + p + 2) = 1.d0

      do i = p + 2, n + 1
         U(i) = real(i - p - 1)/real(n - p + 1)
      end do

   end subroutine FillOpenKnot

!---------------------------------------------------------------------------
!> @brief
!> Calculates number of elements (nonempty subintervals) in the knot vector.
!
! Input:
! ------
!> @param[in] n  - number of functions (control points) minus 1
!> @param[in] p  - order of basis functions
!> @param[in] U  - knot vector
!
! Output:
! -------
!> @return nelem - number of elements
!---------------------------------------------------------------------------
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

   subroutine repeatedKnot(n, p, iblock, U, nelem)
      implicit none
      integer(kind=4), intent(inout) :: n
      integer(kind=4), intent(in) :: p, iblock
      real(kind=8), allocatable, dimension(:), intent(out) :: U
      integer(kind=4), intent(out) :: nelem
      integer(kind=4) :: i, j, k, l

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

end module knot_vector
