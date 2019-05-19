module Interfaces

   interface
! -------------------------------------------------------------------
! Right-hand side of the equation.
!
! Input:
! ------
! ads             - ADS setup structure
! X_              - quadrature points
! k_              - indexes for quadrature points
! e_              - indexes for elements
! a_              - indexes of basis functions
! du_             - value of derivative from previous time step
! Uval            - previous solution coefficient at given point
! ads_data        - data structures for ADS
! J               - jacobian
! W               - weight for quadratures
!
! Output:
! -------
! ret             - value of RHS function at given point
!
! -------------------------------------------------------------------
      subroutine RHS_fun_int(&
         ads, &
         X, &
         k, &
         e, &
         a, &
         du, &
         Uval, &
         ads_data, J, W, l2norm, ret)
         use Setup
         implicit none
         type (ADS_setup), intent(in) :: ads
         real (kind = 8), intent(in), dimension(3) :: X
         integer(kind = 4), intent(in), dimension(3) :: k
         integer(kind = 4), intent(in), dimension(3) :: e
         integer(kind = 4), intent(in), dimension(3) :: a
         real (kind = 8), intent(in), dimension(3) :: du
         real (kind = 8), intent(in) :: Uval
         type (ADS_compute_data), intent(in) :: ads_data
         real (kind = 8), intent(in) :: J, W
         real (kind = 8), intent(out) :: l2norm
         real (kind = 8), intent(out) :: ret
      end subroutine

      
! -------------------------------------------------------------------
! Create 1D Matrix from Mass-Matrix, Stifness-Matrix and Advection-Matrix
!
! Input:
! ------
! KL     - number of lower diagonals of the resulting matrix
! KU     - number of upper diagonals of the resulting matrix
! n      - number of control points minus one
! M      - mass matrix
! K      - stifness matrix
! A      - advection matrix
!
! Output:
! -------
! O      - combination of M, K, A
!
! -------------------------------------------------------------------
      subroutine Form1DMatrix (KL,KU,n,M,K,A, O)
         implicit none
         integer(kind = 4), intent(in) :: KL, KU
         integer(kind = 4), intent(in) :: n
         real (kind = 8), intent(in) :: M(0:(2 * KL + KU), 0:n)
         real (kind = 8), intent(in) :: K(0:(2 * KL + KU), 0:n)
         real (kind = 8), intent(in) :: A(0:(2 * KL + KU), 0:n)
         real (kind = 8), intent(out) :: O(0:(2 * KL + KU), 0:n)
      end subroutine Form1DMatrix
   end interface 

end module Interfaces