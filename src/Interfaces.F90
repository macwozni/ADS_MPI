module Interfaces

   interface
      function forcing_fun (un, X) result (ret)
         implicit none
         real (kind = 8), intent(in), dimension(:) :: un
         real (kind = 8), intent(in), dimension(3) :: X
         real (kind = 8) :: ret
      end function

   end interface


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
! n               - nuber of previous time steps
! Un              - U_n, previous solution coefficient at given point
! Un13            - U_n+1/3
! Un23            - U_n+2/3
! ads_data        - data structures for ADS
! J               - jacobian
! W               - weight for quadratures
! directon        - direction for the substep
! substep         - number of substep
!
! Output:
! -------
! ret             - value of RHS function at given point
! l2norm          -
!
! -------------------------------------------------------------------
      subroutine RHS_fun_int(&
         ads, &
         X, &
         k, &
         e, &
         a, &
         du, &
         n, &
         un, &
         un13, &
         un23, &
         ads_data, J, W, direction, substep, l2norm, ret)
         use Setup
         implicit none
         type (ADS_setup), intent(in) :: ads
         real (kind = 8), intent(in), dimension(3) :: X
         integer(kind = 4), intent(in), dimension(3) :: k
         integer(kind = 4), intent(in), dimension(3) :: e
         integer(kind = 4), intent(in), dimension(3) :: a
         real (kind = 8), intent(in), dimension(3) :: du
         integer (kind = 4), intent(in) :: n
         real (kind = 8), intent(in), dimension(n)  :: un
         real (kind = 8), intent(in) :: un13,un23
         type (ADS_compute_data), intent(in) :: ads_data
         real (kind = 8), intent(in) :: J, W
         integer (kind=4), intent(in) :: direction
         integer (kind=4), intent(in) :: substep
         real (kind = 8), intent(out) :: l2norm
         real (kind = 8), intent(out) :: ret
      end subroutine

   end interface

end module Interfaces