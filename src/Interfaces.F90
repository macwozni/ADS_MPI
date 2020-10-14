!------------------------------------------------------------------------------
! AGH University of Science and Technology, Krakow, Poland
!------------------------------------------------------------------------------
!
! MODULE: Interfaces
!
!> @author
!> Maciej Wozniak
!
! DESCRIPTION:
!> This module contains interfaces for subroutines to be implemented for each problem.
!
! REVISION HISTORY:
! 13 10 2020 - Initial Version
! 
!------------------------------------------------------------------------------



module Interfaces

   interface
!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief Interface for Right-hand-side of the equation.
!> Subroutine that is implemented as RHS-equation should comply to this interface
!
! Input:
! ------
!> @param[in] ads         - ADS setup structure
!> @param[in] X_          - quadrature points
!> @param[in] k_          - indexes for quadrature points
!> @param[in] e_          - indexes for elements
!> @param[in] a_          - indexes of basis functions
!> @param[in] du_         - value of derivative from previous time step
!> @param[in] n           - nuber of previous time steps
!> @param[in] Un          - \f$ U_n \f$, previous solution coefficient at given point
!> @param[in] Un13        - \f$ U_{n+1/3} \f$
!> @param[in] Un23        - \f$ U_{n+2/3} \f$
!> @param[in] ads_data    - data structures for ADS
!> @param[in] J           - jacobian
!> @param[in] W           - weight for quadratures
!> @param[in] directon    - direction for the substep
!> @param[in] substep     - number of substep
!
! Output:
! -------
!> @param[out] ret        - value of RHS function at given point
!> @param[out] l2norm     -
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