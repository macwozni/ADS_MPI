!------------------------------------------------------------------------------
!
! MODULE: gauss
!
! DESCRIPTION:
!> @file gauss.F90
!> @brief Module providing Gauss-Legendre quadrature points and weights.
!>
!> @details
!> This module stores precomputed or explicitly constructed sets of
!> quadrature points and weights for the Gauss-Legendre rule on the
!> reference interval \f$[-1,1]\f$ for numbers of integration points
!> from 1 to 10.
!>
!> The public interface is limited to the procedure \ref GaussRule,
!> which returns the appropriate set of points and weights for a given
!> quadrature order. The internal data are initialized lazily on the
!> first call.
!>
!> @note
!> For \f$n\f$ quadrature points, the rule is exact for polynomials
!> of degree not greater than \f$2n-1\f$.
!>
!> @warning
!> If the argument \p n is outside the range from 1 to 10, the
!> procedure \ref GaussRule returns arrays filled with zeros.
!
!------------------------------------------------------------------------------
module gauss

   implicit none

!> @brief Maximum number of points available in the Gauss rule table.
   integer(kind=4), parameter, public :: MAX_GAUSS_POINTS = 10

!> @brief Initialization status flag for module data.
!
!> @details
!> This variable indicates whether all private arrays containing
!> quadrature points and weights have already been prepared by
!> the procedure \ref initialize. It is set to `.TRUE.` after
!> successful completion of the initialization step.
   logical :: initialized = .FALSE.

!> @brief Gauss-Legendre quadrature point(s) for the 1-point rule.
   real(kind=8), dimension(0:0) :: X1
!> @brief Gauss-Legendre quadrature weight(s) for the 1-point rule.
   real(kind=8), dimension(0:0) :: W1

!> @brief Gauss-Legendre quadrature points for the 2-point rule.
   real(kind=8), dimension(0:1) :: X2
!> @brief Gauss-Legendre quadrature weights for the 2-point rule.
   real(kind=8), dimension(0:1) :: W2

!> @brief Gauss-Legendre quadrature points for the 3-point rule.
   real(kind=8), dimension(0:2) :: X3
!> @brief Gauss-Legendre quadrature weights for the 3-point rule.
   real(kind=8), dimension(0:2) :: W3

!> @brief Gauss-Legendre quadrature points for the 4-point rule.
   real(kind=8), dimension(0:3) :: X4
!> @brief Gauss-Legendre quadrature weights for the 4-point rule.
   real(kind=8), dimension(0:3) :: W4

!> @brief Gauss-Legendre quadrature points for the 5-point rule.
   real(kind=8), dimension(0:4) :: X5
!> @brief Gauss-Legendre quadrature weights for the 5-point rule.
   real(kind=8), dimension(0:4) :: W5

!> @brief Gauss-Legendre quadrature points for the 6-point rule.
   real(kind=8), dimension(0:5) :: X6
!> @brief Gauss-Legendre quadrature weights for the 6-point rule.
   real(kind=8), dimension(0:5) :: W6

!> @brief Gauss-Legendre quadrature points for the 7-point rule.
   real(kind=8), dimension(0:6) :: X7
!> @brief Gauss-Legendre quadrature weights for the 7-point rule.
   real(kind=8), dimension(0:6) :: W7

!> @brief Gauss-Legendre quadrature points for the 8-point rule.
   real(kind=8), dimension(0:7) :: X8
!> @brief Gauss-Legendre quadrature weights for the 8-point rule.
   real(kind=8), dimension(0:7) :: W8

!> @brief Gauss-Legendre quadrature points for the 9-point rule.
   real(kind=8), dimension(0:8) :: X9
!> @brief Gauss-Legendre quadrature weights for the 9-point rule.
   real(kind=8), dimension(0:8) :: W9

!> @brief Gauss-Legendre quadrature points for the 10-point rule.
   real(kind=8), dimension(0:9) :: X10
!> @brief Gauss-Legendre quadrature weights for the 10-point rule.
   real(kind=8), dimension(0:9) :: W10

   PRIVATE :: X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
   PRIVATE :: W1, W2, W3, W4, W5, W6, W7, W8, W9, W10
   PROTECTED :: initialized

contains


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Initializes all quadrature point and weight sets.
!> Prepares the private module arrays containing Gauss-Legendre
!> quadrature points and weights for rules with 1 to 10 integration
!> points.
!>
!> For lower orders, explicit analytical formulas are used. For higher
!> orders, the values are constructed by algebraic or iterative
!> expressions implemented directly in this procedure.
!>
!> After successful completion, the procedure sets the flag
!> \ref initialized to `.TRUE.`, preventing repeated initialization
!> during subsequent calls to \ref GaussRule.
!
! Notes:
! ------
!> @note
!> This procedure is internal to the module and is not intended to be
!> the primary user-level interface.
!
!> @warning
!> This procedure assumes floating-point arithmetic compatible with
!> `real(kind=8)`.
!
!---------------------------------------------------------------------------
   subroutine initialize()
      real(kind=8) :: PI, A, D0, D1, Y0, Y1, Y2, Y3, S, M, EPS
      integer(kind=4) :: I, J
      real(kind=8), dimension(0:4) :: WP, Z
      initialized = .TRUE.

      PI = 4.0d0*atan(1.0d0)

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

      ! p = 11
      A = acos(sqrt(15.0d0)/35.0d0)/3.0d0
      Y0 = 5.0d0/11.0d0 + 4.0d0*sqrt(5.0d0/3.0d0)*cos(A)/11.0d0
      Y1 = 5.0d0/11.0d0 + 4.0d0*sqrt(5.0d0/3.0d0)*cos(A + 4.0d0*PI/3.0d0)/11.0d0
      Y2 = 5.0d0/11.0d0 + 4.0d0*sqrt(5.0d0/3.0d0)*cos(A + 2.0d0*PI/3.0d0)/11.0d0
      X6(0) = -sqrt(Y0)
      X6(1) = -sqrt(Y1)
      X6(2) = -sqrt(Y2)
      X6(3) = -X6(2)
      X6(4) = -X6(1)
      X6(5) = -X6(0)
      W6(0) = 128.0d0/(441.0d0*Y0*(1.0d0 - Y0)*(33.0d0*Y0*Y0 - 30.0d0*Y0 + 5.0d0)**2)
      W6(1) = 128.0d0/(441.0d0*Y1*(1.0d0 - Y1)*(33.0d0*Y1*Y1 - 30.0d0*Y1 + 5.0d0)**2)
      W6(2) = 128.0d0/(441.0d0*Y2*(1.0d0 - Y2)*(33.0d0*Y2*Y2 - 30.0d0*Y2 + 5.0d0)**2)
      W6(3) = W6(2)
      W6(4) = W6(1)
      W6(5) = W6(0)

      ! p = 13
      A = acos(-sqrt(231.0d0)/189.0d0)/3.0d0
      Y0 = 7.0d0/13.0d0 + 4.0d0*sqrt(231.0d0)*cos(A)/143.0d0
      Y1 = 7.0d0/13.0d0 + 4.0d0*sqrt(231.0d0)*cos(A + 4.0d0*PI/3.0d0)/143.0d0
      Y2 = 7.0d0/13.0d0 + 4.0d0*sqrt(231.0d0)*cos(A + 2.0d0*PI/3.0d0)/143.0d0
      X7(0) = -sqrt(Y0)
      X7(1) = -sqrt(Y1)
      X7(2) = -sqrt(Y2)
      X7(3) = 0.0d0
      X7(4) = -X7(2)
      X7(5) = -X7(1)
      X7(6) = -X7(0)
      W7(0) = 128.0d0/(9.0d0*Y0*Y0*(1.0d0 - Y0)*(429.0d0*Y0*Y0 - 462.0d0*Y0 + 105.0d0)**2)
      W7(1) = 128.0d0/(9.0d0*Y1*Y1*(1.0d0 - Y1)*(429.0d0*Y1*Y1 - 462.0d0*Y1 + 105.0d0)**2)
      W7(2) = 128.0d0/(9.0d0*Y2*Y2*(1.0d0 - Y2)*(429.0d0*Y2*Y2 - 462.0d0*Y2 + 105.0d0)**2)
      W7(3) = 512.0d0/1225.0d0
      W7(4) = W7(2)
      W7(5) = W7(1)
      W7(6) = W7(0)

      ! p = 15
      A = acos(sqrt(2310.0d0)/55.0d0)/3.0d0
      M = 448.0d0/2925.0d0 + 32.0d0*sqrt(2310.0d0)*cos(A)/6435.0d0
      S = sqrt(M)
      D0 = 448.0d0/975.0d0 - M - 3584.0d0/(482625.0d0*S)
      D1 = 448.0d0/975.0d0 - M + 3584.0d0/(482625.0d0*S)
      Y0 = 7.0d0/15.0d0 - 0.5d0*S - 0.5d0*sqrt(D0)
      Y1 = 7.0d0/15.0d0 - 0.5d0*S + 0.5d0*sqrt(D0)
      Y2 = 7.0d0/15.0d0 + 0.5d0*S - 0.5d0*sqrt(D1)
      Y3 = 7.0d0/15.0d0 + 0.5d0*S + 0.5d0*sqrt(D1)
      X8(0) = -sqrt(Y3)
      X8(1) = -sqrt(Y2)
      X8(2) = -sqrt(Y1)
      X8(3) = -sqrt(Y0)
      X8(4) = -X8(3)
      X8(5) = -X8(2)
      X8(6) = -X8(1)
      X8(7) = -X8(0)
      W8(0) = 512.0d0/((1.0d0 - Y3)*Y3*(6435.0d0*Y3*Y3*Y3 - 9009.0d0*Y3*Y3 + 3465.0d0*Y3 - 315.0d0)**2)
      W8(1) = 512.0d0/((1.0d0 - Y2)*Y2*(6435.0d0*Y2*Y2*Y2 - 9009.0d0*Y2*Y2 + 3465.0d0*Y2 - 315.0d0)**2)
      W8(2) = 512.0d0/((1.0d0 - Y1)*Y1*(6435.0d0*Y1*Y1*Y1 - 9009.0d0*Y1*Y1 + 3465.0d0*Y1 - 315.0d0)**2)
      W8(3) = 512.0d0/((1.0d0 - Y0)*Y0*(6435.0d0*Y0*Y0*Y0 - 9009.0d0*Y0*Y0 + 3465.0d0*Y0 - 315.0d0)**2)
      W8(4) = W8(3)
      W8(5) = W8(2)
      W8(6) = W8(1)
      W8(7) = W8(0)

      ! p = 17
      A = acos(sqrt(6006.0d0)/91.0d0)/3.0d0
      M = -48.0d0/1445.0d0 + 16.0d0*sqrt(6006.0d0)*cos(A)/12155.0d0
      S = sqrt(2.0d0*M + 288.0d0/1445.0d0)
      D0 = -2.0d0*M + 288.0d0/1445.0d0 + 1536.0d0/(319345.0d0*S)
      D1 = -2.0d0*M + 288.0d0/1445.0d0 - 1536.0d0/(319345.0d0*S)
      Y0 = 9.0d0/17.0d0 - 0.5d0*S - 0.5d0*sqrt(D0)
      Y1 = 9.0d0/17.0d0 - 0.5d0*S + 0.5d0*sqrt(D0)
      Y2 = 9.0d0/17.0d0 + 0.5d0*S - 0.5d0*sqrt(D1)
      Y3 = 9.0d0/17.0d0 + 0.5d0*S + 0.5d0*sqrt(D1)
      X9(0) = -sqrt(Y3)
      X9(1) = -sqrt(Y2)
      X9(2) = -sqrt(Y1)
      X9(3) = -sqrt(Y0)
      X9(4) = 0.0d0
      X9(5) = -X9(3)
      X9(6) = -X9(2)
      X9(7) = -X9(1)
      X9(8) = -X9(0)
      W9(0) = 512.0d0/(81.0d0*(1.0d0 - Y3)*(715.0d0*Y3*Y3*Y3 - 1001.0d0*Y3*Y3 + 385.0d0*Y3 - 35.0d0)**2)
      W9(1) = 512.0d0/(81.0d0*(1.0d0 - Y2)*(715.0d0*Y2*Y2*Y2 - 1001.0d0*Y2*Y2 + 385.0d0*Y2 - 35.0d0)**2)
      W9(2) = 512.0d0/(81.0d0*(1.0d0 - Y1)*(715.0d0*Y1*Y1*Y1 - 1001.0d0*Y1*Y1 + 385.0d0*Y1 - 35.0d0)**2)
      W9(3) = 512.0d0/(81.0d0*(1.0d0 - Y0)*(715.0d0*Y0*Y0*Y0 - 1001.0d0*Y0*Y0 + 385.0d0*Y0 - 35.0d0)**2)
      W9(4) = 32768.0d0/99225.0d0
      W9(5) = W9(3)
      W9(6) = W9(2)
      W9(7) = W9(1)
      W9(8) = W9(0)

      ! p = 19
      EPS = 1.0d-15
      do I = 0, 4
         Z(I) = cos((4.0d0*dble(I + 1) - 1.0d0)*acos(-1.0d0)/42.0d0)
         do J = 1, 50
            Y0 = Z(I)*Z(I)
            M = (((((46189.0d0*Y0 - 109395.0d0)*Y0 + 90090.0d0)*Y0 - 30030.0d0)*Y0 + 3465.0d0)*Y0 - 63.0d0)/256.0d0
            S = Z(I)*(230945.0d0*Y0*Y0*Y0*Y0 - 437580.0d0*Y0*Y0*Y0 + 270270.0d0*Y0*Y0 - 60060.0d0*Y0 + 3465.0d0)/128.0d0
            A = Z(I)
            Z(I) = Z(I) - M/S
            if (abs(Z(I) - A) < EPS) exit
         end do
         Y0 = Z(I)*Z(I)
         S = Z(I)*(230945.0d0*Y0*Y0*Y0*Y0 - 437580.0d0*Y0*Y0*Y0 + 270270.0d0*Y0*Y0 - 60060.0d0*Y0 + 3465.0d0)/128.0d0
         WP(I) = 2.0d0/((1.0d0 - Y0)*S*S)
      end do
      X10(0) = -Z(0)
      X10(1) = -Z(1)
      X10(2) = -Z(2)
      X10(3) = -Z(3)
      X10(4) = -Z(4)
      X10(5) = -X10(4)
      X10(6) = -X10(3)
      X10(7) = -X10(2)
      X10(8) = -X10(1)
      X10(9) = -X10(0)
      W10(0) = WP(0)
      W10(1) = WP(1)
      W10(2) = WP(2)
      W10(3) = WP(3)
      W10(4) = WP(4)
      W10(5) = W10(4)
      W10(6) = W10(3)
      W10(7) = W10(2)
      W10(8) = W10(1)
      W10(9) = W10(0)
   end subroutine initialize


!---------------------------------------------------------------------------
!> @brief Returns Gauss-Legendre quadrature points and weights
!>        for a given number of integration points.
!
!> @details
!> This procedure fills the arrays \p X and \p W with the quadrature
!> points and weights corresponding to the Gauss-Legendre rule defined
!> by the requested number of integration points \p n.
!>
!> If the module data have not been initialized yet, the procedure
!> automatically calls \ref initialize.
!>
!> An \f$n\f$-point Gauss-Legendre quadrature rule is exact for
!> polynomials of degree not greater than \f$2n-1\f$.
!
! Input:
! ------
!> @param[in]  n
!> Number of quadrature points. Valid range: from 1 to 10.
!
! Output:
! -------
!> @param[out] X
!> Array of quadrature points on the reference interval \f$[-1,1]\f$.
!> Its size must be consistent with the declaration
!> `dimension(0:n-1)`.
!
!> @param[out] W
!> Array of quadrature weights corresponding to the points stored
!> in \p X. Its size must be consistent with the declaration
!> `dimension(0:n-1)`.
!
! Notes:
! ------
!> @note
!> The points and weights are returned in ascending order with respect
!> to the quadrature points.
!
!> @warning
!> If \p n is outside the range from 1 to 10, the procedure does not
!> raise an error; instead, it fills both \p X and \p W with zeros.
!---------------------------------------------------------------------------
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
      case (6) ! p = 11
         X = X6
         W = W6
      case (7) ! p = 13
         X = X7
         W = W7
      case (8) ! p =15
         X = X8
         W = W8
      case (9) ! p = 17
         X = X9
         W = W9
      case (10) ! p = 19
         X = X10
         W = W10
      case default
         X = 0.0d0
         W = 0.0d0
      end select

   end subroutine GaussRule

end module gauss
