!------------------------------------------------------------------------------
!
! MODULE: gauss
!
! DESCRIPTION:
!> This module contains all functionality associated to creating gauss integration points and weights
!
!------------------------------------------------------------------------------

module gauss

   implicit none

!> Checks if private varaibles are initialized
   logical :: initialized = .FALSE.

   real(kind=8), dimension(0:0) :: X1
   real(kind=8), dimension(0:0) :: W1

   real(kind=8), dimension(0:1) :: X2
   real(kind=8), dimension(0:1) :: W2

   real(kind=8), dimension(0:2) :: X3
   real(kind=8), dimension(0:2) :: W3

   real(kind=8), dimension(0:3) :: X4
   real(kind=8), dimension(0:3) :: W4

   real(kind=8), dimension(0:4) :: X5
   real(kind=8), dimension(0:4) :: W5

   real(kind=8), dimension(0:5) :: X6
   real(kind=8), dimension(0:5) :: W6

   real(kind=8), dimension(0:6) :: X7
   real(kind=8), dimension(0:6) :: W7

   real(kind=8), dimension(0:7) :: X8
   real(kind=8), dimension(0:7) :: W8

   real(kind=8), dimension(0:8) :: X9
   real(kind=8), dimension(0:8) :: W9

   real(kind=8), dimension(0:9) :: X10
   real(kind=8), dimension(0:9) :: W10

   PRIVATE :: X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
   PRIVATE :: W1, W2, W3, W4, W5, W6, W7, W8, W9, W10
   PROTECTED :: initialized

contains

!---------------------------------------------------------------------------
!> @brief
!> Precomputes all Gaussian quadrature points and weights
! -------------------------------------------------------------------
   subroutine initialize()
      initialized = .TRUE.
      ! p = 1
      X1(0) = 0.0d0
      W1(0) = 2.0d0
      ! p = 3
      X2(0) = -1.0d0/sqrt(3.0d0)
      X2(1) = -X2(0)
      W2(0) = 1.0d0
      W2(1) = W2(0)
      ! p = 5
      X3(0) = -sqrt(3.0d0/5.0d0)
      X3(1) = 0.0d0
      X3(2) = -X3(0)
      W3(0) = 5.0d0/9.0d0
      W3(1) = 8.0d0/9.0d0
      W3(2) = W3(0)
      ! p = 7
      X4(0) = -sqrt((3.0d0 + 2.0d0*sqrt(6.0d0/5.0d0))/7.0d0)
      X4(1) = -sqrt((3.0d0 - 2.0d0*sqrt(6.0d0/5.0d0))/7.0d0)
      X4(2) = -X4(1)
      X4(3) = -X4(0)
      W4(0) = (18.0d0 - sqrt(30.0d0))/36.0d0
      W4(1) = (18.0d0 + sqrt(30.0d0))/36.0d0
      W4(2) = W4(1)
      W4(3) = W4(0)
      ! p = 9
      X5(0) = -1.0d0/3.0d0*sqrt(5.0d0 + 2.0d0*sqrt(10.0d0/7.0d0))
      X5(1) = -1.0d0/3.0d0*sqrt(5.0d0 - 2.0d0*sqrt(10.0d0/7.0d0))
      X5(2) = 0.0d0
      X5(3) = -X5(1)
      X5(4) = -X5(0)
      W5(0) = (322.0d0 - 13.0d0*sqrt(70.0d0))/900.0d0
      W5(1) = (322.0d0 + 13.0d0*sqrt(70.0d0))/900.0d0
      W5(2) = 128.0d0/225.0d0
      W5(3) = W5(1)
      W5(4) = W5(0)

      X6(0) = -0.9324695142031520d0
      X6(1) = -0.6612093864662645d0
      X6(2) = -0.2386191860831969d0
      X6(3) = -X6(2)
      X6(4) = -X6(1)
      X6(5) = -X6(0)
      W6(0) = 0.1713244923791703d0
      W6(1) = 0.3607615730481386d0
      W6(2) = 0.4679139345726911d0
      W6(3) = W6(2)
      W6(4) = W6(1)
      W6(5) = W6(0)

      X7(0) = -0.9491079123427585d0
      X7(1) = -0.7415311855993944d0
      X7(2) = -0.4058451513773972d0
      X7(3) = 0.0d0
      X7(4) = -X7(2)
      X7(5) = -X7(1)
      X7(6) = -X7(0)
      W7(0) = 0.1294849661688697d0
      W7(1) = 0.2797053914892767d0
      W7(2) = 0.3818300505051189d0
      W7(3) = 0.4179591836734694d0
      W7(4) = W7(2)
      W7(5) = W7(1)
      W7(6) = W7(0)

      X8(0) = -0.9602898564975362d0
      X8(1) = -0.7966664774136267d0
      X8(2) = -0.5255324099163290d0
      X8(3) = -0.1834346424956498d0
      X8(4) = -X8(3)
      X8(5) = -X8(2)
      X8(6) = -X8(1)
      X8(7) = -X8(0)
      W8(0) = 0.1012285362903763d0
      W8(1) = 0.2223810344533745d0
      W8(2) = 0.3137066458778873d0
      W8(3) = 0.3626837833783620d0
      W8(4) = W8(3)
      W8(5) = W8(2)
      W8(6) = W8(1)
      W8(7) = W8(0)

      X9(0) = -0.9681602395076261d0
      X9(1) = -0.8360311073266358d0
      X9(2) = -0.6133714327005904d0
      X9(3) = -0.3242534234038089d0
      X9(4) = 0.d0
      X9(5) = -X9(3)
      X9(6) = -X9(2)
      X9(7) = -X9(1)
      X9(8) = -X9(0)
      W9(0) = 0.0812743883615744d0
      W9(1) = 0.1806481606948574d0
      W9(2) = 0.2606106964029354d0
      W9(3) = 0.3123470770400029d0
      W9(4) = 0.3302393550012598d0
      W9(5) = W9(3)
      W9(6) = W9(2)
      W9(7) = W9(1)
      W9(8) = W9(0)

      X10(0) = -0.973906528517172d0
      X10(1) = -0.865063366688985d0
      X10(2) = -0.679409568299024d0
      X10(3) = -0.433395394129247d0
      X10(4) = -0.148874338981631d0
      X10(5) = -X10(4)
      X10(6) = -X10(3)
      X10(7) = -X10(2)
      X10(8) = -X10(1)
      X10(9) = -X10(0)
      W10(0) = 0.066671344308688d0
      W10(1) = 0.149451349150581d0
      W10(2) = 0.219086362515982d0
      W10(3) = 0.269266719309996d0
      W10(4) = 0.295524224714753d0
      W10(5) = W10(4)
      W10(6) = W10(3)
      W10(7) = W10(2)
      W10(8) = W10(1)
      W10(9) = W10(0)
   end subroutine initialize

!---------------------------------------------------------------------------
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
      integer(kind=4), intent(in) :: n
      real(kind=8), dimension(0:n - 1), intent(out) :: X
      real(kind=8), dimension(0:n - 1), intent(out) :: W

      if (.NOT. initialized) call initialize()

      select case (n)
      case (1) ! p = 1
         X = X1
         W = W1
      case (2) ! p = 3
         X = X2
         W = W2
      case (3) ! p = 5
         X = X3
         W = W3
      case (4) ! p = 7
         X = X4
         W = W4
      case (5) ! p = 9
         X = X5
         W = W5
      case (6)
         X = X6
         W = W6
      case (7)
         X = X7
         W = W7
      case (8)
         X = X8
         W = W8
      case (9)
         X = X9
         W = W9
      case (10)
         X = X10
         W = W10
      case default
         X = 0.0d0
         W = 0.0d0
      end select

   end subroutine GaussRule

end module gauss

