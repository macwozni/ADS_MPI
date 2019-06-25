!------------------------------------------------------------------------------
! AGH University of Science and Technology, Krakow, Poland
!------------------------------------------------------------------------------
!
! MODULE: gauss
!
!> @author
!> Maciej Wozniak
!
! DESCRIPTION:
!> This module contains all functionality associated to creating gauss integration points and weights
!
! REVISION HISTORY:
! 21 11 2017 - Initial Version
! 
!------------------------------------------------------------------------------


module gauss

implicit none

contains


!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> Fills supplied arrays with points and weights for Gauss quadrature.
!
! Input:
! ------
!> @param[in] n  - number of quadrature points
!
! Output:
! ------
!> @param[out] X - coordinates of quadrature points
!> @param[out] W - weights
!
!> Note: \f$ n \f$- point Gauss quadrature yields exact results for polynomials
!> of degree up to \f$ 2n - 1 \f$.
! -------------------------------------------------------------------
subroutine GaussRule(n, X, W)
   implicit none
   integer(kind = 4), intent(in) :: n
   real (kind = 8), intent(out) :: X(0:n - 1)
   real (kind = 8), intent(out) :: W(0:n - 1)

   select case (n)
   case (1) ! p = 1
      X(0) = 0.0d0
      W(0) = 2.0d0
   case (2) ! p = 3
      X(0) = - 1.0d0 / sqrt(3.0d0)
      X(1) = -X(0)
      W(0) = 1.0d0
      W(1) = W(0)
   case (3) ! p = 5
      X(0) = - sqrt(3.0d0 / 5.0d0)
      X(1) = 0.0d0
      X(2) = -X(0)
      W(0) = 5.0d0/9.0d0
      W(1) = 8.0d0/9.0d0
      W(2) = W(0)
   case (4) ! p = 7
      X(0) = -sqrt((3.0d0+2.0d0*sqrt(6.0d0/5.0d0))/7.0d0)
      X(1) = -sqrt((3.0d0-2.0d0*sqrt(6.0d0/5.0d0))/7.0d0)
      X(2) = -X(1)
      X(3) = -X(0)
      W(0) = (18.0d0-sqrt(30.0d0))/36.0d0
      W(1) = (18.0d0+sqrt(30.0d0))/36.0d0
      W(2) = W(1)
      W(3) = W(0)
   case (5) ! p = 9
      X(0) = -1.0d0/3.0d0*sqrt(5.0d0+2.0d0*sqrt(10.0d0/7.0d0))
      X(1) = -1.0d0/3.0d0*sqrt(5.0d0-2.0d0*sqrt(10.0d0/7.0d0))
      X(2) = 0.0d0
      X(3) = -X(1)
      X(4) = -X(0)
      W(0) = (322.0d0-13.0d0*sqrt(70.0d0))/900.0d0
      W(1) = (322.0d0+13.0d0*sqrt(70.0d0))/900.0d0
      W(2) = 128.0d0/225.0d0
      W(3) = W(1)
      W(4) = W(0)
   case (6)
      X(0) = -0.9324695142031520d0
      X(1) = -0.6612093864662645d0
      X(2) = -0.2386191860831969d0
      X(3) = -X(2)
      X(4) = -X(1)
      X(5) = -X(0)
      W(0) = 0.1713244923791703d0
      W(1) = 0.3607615730481386d0
      W(2) = 0.4679139345726911d0
      W(3) = W(2)
      W(4) = W(1)
      W(5) = W(0)
   case (7)
      X(0) = -0.9491079123427585d0
      X(1) = -0.7415311855993944d0
      X(2) = -0.4058451513773972d0
      X(3) = 0.0
      X(4) = -X(2)
      X(5) = -X(1)
      X(6) = -X(0)
      W(0) = 0.1294849661688697d0
      W(1) = 0.2797053914892767d0
      W(2) = 0.3818300505051189d0
      W(3) = 0.4179591836734694d0
      W(4) = W(2)
      W(5) = W(1)
      W(6) = W(0)
   case (8)
      X(0) = -0.9602898564975362d0
      X(1) = -0.7966664774136267d0
      X(2) = -0.5255324099163290d0
      X(3) = -0.1834346424956498d0
      X(4) = -X(3)
      X(5) = -X(2)
      X(6) = -X(1)
      X(7) = -X(0)
      W(0) = 0.1012285362903763d0
      W(1) = 0.2223810344533745d0
      W(2) = 0.3137066458778873d0
      W(3) = 0.3626837833783620d0
      W(4) = W(3)
      W(5) = W(2)
      W(6) = W(1)
      W(7) = W(0)
   case (9)
      X(0) = -0.9681602395076261d0
      X(1) = -0.8360311073266358d0
      X(2) = -0.6133714327005904d0
      X(3) = -0.3242534234038089d0
      X(4) = 0.0
      X(5) = -X(3)
      X(6) = -X(2)
      X(7) = -X(1)
      X(8) = -X(0)
      W(0) = 0.0812743883615744d0
      W(1) = 0.1806481606948574d0
      W(2) = 0.2606106964029354d0
      W(3) = 0.3123470770400029d0
      W(4) = 0.3302393550012598d0
      W(5) = W(3)
      W(6) = W(2)
      W(7) = W(1)
      W(8) = W(0)
   case (10)
      X(0) = -0.973906528517172d0
      X(1) = -0.865063366688985d0
      X(2) = -0.679409568299024d0
      X(3) = -0.433395394129247d0
      X(4) = -0.148874338981631d0
      X(5) = -X(4)
      X(6) = -X(3)
      X(7) = -X(2)
      X(8) = -X(1)
      X(9) = -X(0)
      W(0) = 0.066671344308688d0
      W(1) = 0.149451349150581d0
      W(2) = 0.219086362515982d0
      W(3) = 0.269266719309996d0
      W(4) = 0.295524224714753d0
      W(5) = W(4)
      W(6) = W(3)
      W(7) = W(2)
      W(8) = W(1)
      W(9) = W(0)
   case default
      X = 0.0d0
      W = 0.0d0
   end select

end subroutine GaussRule

end module gauss


