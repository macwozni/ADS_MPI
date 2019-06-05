!------------------------------------------------------------------------------
! AGH University of Science and Technology, Krakow, Poland
!------------------------------------------------------------------------------
!
! MODULE: knot_vector
!
!> @author
!> Maciej Wozniak
!
! DESCRIPTION: 
!> This module contains all functionality associated to creating and filling knot vectors.
!
! REVISION HISTORY:
! 21 11 2017 - Initial Version
! 
!------------------------------------------------------------------------------

module knot_vector


contains

!---------------------------------------------------------------------------  
!> @author 
!> Maciej Wozniak
!
! DESCRIPTION: 
!> Allocates and fills the knot vector on \f$  [0, 1] \f$.
!>Number of subintervals is \f$ N = n-p+1 \f$.
!> \f$ 0 \f$ and \f$ 1 \f$ are repeated \f$ (p+1) \f$ times.
!
! REVISION HISTORY:
! 21 11 2017 - Initial Version
!
!> @param[in] n  - number of functions on the knot minus one
!> @param[in] p  - degree of polynomial
!> @param[out] U  - array to fill with points
!> @param[out] nelem - number of elements
!---------------------------------------------------------------------------  
subroutine PrepareKnot(U, n, p, nelem)
   implicit none
   integer(kind = 4), intent(in) :: n, p
   real (kind = 8), allocatable, dimension(:), intent(out) :: U
   integer(kind = 4), intent(out) :: nelem


   allocate(U(n + p + 2))
   call FillOpenKnot(U, n, p)
   nelem = CountSpans(n, p, U)

#ifdef IINFO
   write(*, *) 'n,p,nelem', n, p, nelem
   write(*, *) 'U', U
#endif

end subroutine PrepareKnot


!---------------------------------------------------------------------------  
!> @author 
!> Maciej Wozniak
!
! DESCRIPTION: 
!> Fills knot vector on [0, 1].
!> Number of subintervals is \f$ N = n-p+1 \f$.
!> \f$ 0 \f$ and \f$ 1 \f$ are repeated \f$ (p+1) \f$ times.
!
! REVISION HISTORY:
! 21 11 2017 - Initial Version
!
!> @param[in] n  - number of functions on the knot minus one
!> @param[in] p  - degree of polynomial
!> @param[out] U  - array to fill with points
!--------------------------------------------------------------------------- 
subroutine FillOpenKnot(U, n, p)
   implicit none
   integer(kind = 4), intent(in) :: n, p
   real (kind = 8), dimension(1:n + p + 2), intent(out) :: U
   integer(kind = 4) :: i

   U(1: p + 1) = 0.d0
   U(n + 2: n + p + 2) = 1.d0

   do i = p + 2, n + 1
      U(i) = real(i - p - 1) / real(n - p + 1)
   enddo

end subroutine FillOpenKnot

!---------------------------------------------------------------------------  
!> @author 
!> Maciej Wozniak
!
! DESCRIPTION: 
!> Calculates number of elements (nonempty subintervals) in the knot vector.
!
! REVISION HISTORY:
! 21 11 2017 - Initial Version
!
!> @param[in] n - number of functions (control points) minus 1
!> @param[in] p - order of basis functions
!> @param[in] U - knot vector
!> @return nelem - number of elements
!---------------------------------------------------------------------------  
function CountSpans(n, p, U) result (nelem)
   implicit none
   integer(kind = 4), intent(in) :: n, p
   real (kind = 8), dimension(0:n + p + 1), intent(in) :: U
   integer(kind = 4) :: i, nelem

   nelem = 0
   i = p
   do while (i <= n)
      ! skip multiple knots
      do while (i < n .and. U(i) == U(i + 1))
         i = i + 1
      enddo
#ifdef IPRINT
      write(*, *) 'CountSpans:i,n,U(i),U(i+1)', i, n, U(i), U(i + 1)
#endif
      nelem = nelem + 1
      i = i + 1
   enddo

end function CountSpans




end module knot_vector
