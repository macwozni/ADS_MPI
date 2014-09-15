

module core

use debug
!use parallelism
!use projection_engine
!use communicators
!use utils

implicit none
private

public FillOpenKnot


contains

! -------------------------------------------------------------------
! Fills knot vector on [0, 1] 
! 
! U  - array to fill with points
! n  - number of functions on the knot minus one
! p  - degree of polynomial
!
! Number of subintervals is N = n-p+1.
! 0 and 1 are repeated (p+1) times.
! -------------------------------------------------------------------
subroutine FillOpenKnot(U, n, p)
integer :: n, p
real (kind=8), intent(out) :: U(1:n+p+2)
integer :: i

  U(1 : p+1) = 0
  U(n+2 : n+p+2) = 1

  do i = p+2, n+1 
    U(i) = (i-p-1) / real(n-p+1)
  enddo

end subroutine


end module
