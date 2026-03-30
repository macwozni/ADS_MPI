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
         real(kind=8), intent(in) :: un
         real(kind=8), intent(in), dimension(3) :: du
         real(kind=8), intent(in), dimension(3) :: X
         real(kind=8) :: ret
      end function

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Abstract interface for scalar forcing functions used in the
!> right-hand-side assembly process.
!>
!> @details
!> This interface specifies the signature of callback procedures that
!> evaluate a scalar forcing term at a given point, using:
!> - the reconstructed scalar solution value \p un,
!> - the reconstructed directional derivative vector \p du,
!> - the physical coordinates \p X of the evaluation point.
!>
!> The returned value is a scalar contribution later used by the ADS
!> pointwise assembly layer.
!
! Input:
! ------
!> @param[in] un
!> Reconstructed scalar solution value at the current quadrature point.
!>
!> @param[in] du
!> Reconstructed directional derivatives of the solution.
!>
!> @param[in] X
!> Physical coordinates of the current quadrature point.
!
! Output:
! -------
!> @return ret
!> Scalar forcing value evaluated at the supplied state and point.
!
! Notes:
! ------
!> @note
!> The interface is abstract and does not provide an implementation.
!> Concrete problem-dependent forcing functions must conform to this
!> signature.
!
!---------------------------------------------------------------------------
interface
   function forcing_fun(un, du, X) result(ret)
      implicit none
!> @brief Reconstructed scalar solution value at the evaluation point.
      real(kind=8), intent(in) :: un
!> @brief Directional derivatives of the solution.
      real(kind=8), intent(in), dimension(3) :: du
!> @brief Physical coordinates of the evaluation point.
      real(kind=8), intent(in), dimension(3) :: X
!> @brief Scalar forcing value returned by the callback.
      real(kind=8) :: ret
   end function

end interface

end module Interfaces