!------------------------------------------------------------------------------
!
! MODULE: time_scheme
!
! DESCRIPTION:
!> @file time_scheme.F90
!> @brief Time-scheme setup helpers for ADS/iGRM multi-step workflows.
!>
!> @details
!> This module prepares the coefficient tables consumed by \ref MultiStep
!> and performs one-time compatibility checks required by iGRM time
!> schemes. It intentionally does not execute ADS substeps; execution
!> remains in \ref ADSS.
!
!------------------------------------------------------------------------------
module time_scheme

   implicit none

   private

   public :: ValidateIGRMTimeSchemeSpaces
   public :: ConfigureDouglasGunn
   public :: ConfigureDouglasGunn3D
   public :: ConfigurePeacemanRachford3D
   public :: ConfigureBackwardEuler3D

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Resolves whether time-scheme transport terms should be active.
!>
!> @details
!> The time-scheme configurators are diffusion-only by default. Callers
!> can pass \p include_transport to also mirror each directional
!> diffusion coefficient into the first-derivative transport row.
!
! Input:
! ------
!> @param[in] include_transport
!> Optional logical flag selecting transport-term inclusion.
!
! Output:
! -------
!> @return UseTransport
!> `.TRUE.` when transport terms should be included.
!
!---------------------------------------------------------------------------
logical function UseTransport(include_transport)
      implicit none
!> @brief Optional logical flag selecting transport-term inclusion.
      logical, intent(in), optional :: include_transport

      UseTransport = .FALSE.
      if (present(include_transport)) UseTransport = include_transport

end function UseTransport

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Initializes mass-only coefficient tables for \ref MultiStep.
!>
!> @details
!> The returned \p mix table keeps legacy substep behavior mass-only, while
!> \p lhs_mix prepares a physical-direction/substep table whose default is
!> also mass-only. Scheme configurators then overwrite selected axes.
!
! Output:
! -------
!> @param[out] mix
!> Legacy substep matrix-mixing table.
!>
!> @param[out] lhs_mix
!> Directional LHS matrix-mixing table indexed by direction and substep.
!
!---------------------------------------------------------------------------
subroutine ConfigureMassTables(mix, lhs_mix)
      implicit none
!> @brief Legacy substep matrix-mixing table.
      real(kind=8), intent(out), dimension(4, 3) :: mix
!> @brief Directional LHS matrix-mixing table.
      real(kind=8), intent(out), dimension(4, 3, 3) :: lhs_mix

      mix = 0.d0
      mix(1, :) = 1.d0
      lhs_mix = 0.d0
      lhs_mix(1, :, :) = 1.d0

end subroutine ConfigureMassTables

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Sets one implicit directional operator in a time-scheme LHS table.
!>
!> @details
!> The selected \p axis and \p substep receive a mass-plus-directional
!> operator. Diffusion uses row 2 of the mixing vector; optional transport
!> uses row 4, matching the existing matrix assembly convention.
!
! Input:
! ------
!> @param[in] axis
!> Physical direction selected for the implicit operator.
!>
!> @param[in] substep
!> Multi-step substep where the operator is active.
!>
!> @param[in] theta
!> Operator coefficient used in the LHS matrix.
!>
!> @param[in] include_transport
!> Optional logical flag selecting transport-term inclusion.
!
! Input/Output:
! -------------
!> @param[inout] lhs_mix
!> Directional LHS matrix-mixing table updated in place.
!
!---------------------------------------------------------------------------
subroutine SetImplicitAxis(lhs_mix, axis, substep, theta, include_transport)
      implicit none
!> @brief Directional LHS matrix-mixing table updated in place.
      real(kind=8), intent(inout), dimension(4, 3, 3) :: lhs_mix
!> @brief Physical direction and multi-step substep being configured.
      integer(kind=4), intent(in) :: axis, substep
!> @brief Operator coefficient used in the LHS matrix.
      real(kind=8), intent(in) :: theta
!> @brief Optional logical flag selecting transport-term inclusion.
      logical, intent(in), optional :: include_transport

      lhs_mix(:, axis, substep) = 0.d0
      lhs_mix(1, axis, substep) = 1.d0
      lhs_mix(2, axis, substep) = theta
      if (UseTransport(include_transport)) lhs_mix(4, axis, substep) = theta

end subroutine SetImplicitAxis

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Adds one explicit directional operator to an RHS coefficient table.
!>
!> @details
!> The routine maps physical axis numbers to the RHS coefficient rows used
!> by \ref ComputePointForRHS. It also stores which reconstructed
!> derivative buffer should be used for that explicit contribution.
!
! Input:
! ------
!> @param[in] axis
!> Physical direction selected for the explicit operator.
!>
!> @param[in] substep
!> Multi-step substep where the RHS contribution is assembled.
!>
!> @param[in] coeff
!> Coefficient added to the RHS table.
!>
!> @param[in] state
!> Derivative-state selector consumed by \ref ComputePointForRHS.
!>
!> @param[in] include_transport
!> Optional logical flag selecting transport-term inclusion.
!
! Input/Output:
! -------------
!> @param[inout] alpha_step
!> RHS coefficient table updated in place.
!>
!> @param[inout] rhs_du_state
!> RHS derivative-state selector table updated in place.
!
!---------------------------------------------------------------------------
subroutine AddExplicitAxis(alpha_step, rhs_du_state, axis, substep, coeff, state, include_transport)
      implicit none
!> @brief RHS coefficient table updated in place.
      real(kind=8), intent(inout), dimension(7, 3) :: alpha_step
!> @brief RHS derivative-state selector table updated in place.
      integer(kind=4), intent(inout), dimension(6, 3) :: rhs_du_state
!> @brief Physical direction, substep, and derivative-state selector.
      integer(kind=4), intent(in) :: axis, substep, state
!> @brief Coefficient added to the RHS table.
      real(kind=8), intent(in) :: coeff
!> @brief Optional logical flag selecting transport-term inclusion.
      logical, intent(in), optional :: include_transport
!> @brief RHS rows associated with diffusion and transport terms.
      integer(kind=4) :: diffusion_row, transport_row

      select case (axis)
      case (1)
            diffusion_row = 1
            transport_row = 2
      case (2)
            diffusion_row = 3
            transport_row = 4
      case (3)
            diffusion_row = 5
            transport_row = 6
      case default
            stop "wrong time-scheme axis"
      end select

      alpha_step(diffusion_row, substep) = alpha_step(diffusion_row, substep) + coeff
      rhs_du_state(diffusion_row, substep) = state
      if (UseTransport(include_transport)) then
            alpha_step(transport_row, substep) = alpha_step(transport_row, substep) + coeff
            rhs_du_state(transport_row, substep) = state
      end if

end subroutine AddExplicitAxis

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Validates test/trial mesh compatibility for iGRM time schemes.
!>
!> @details
!> All multi-step iGRM schemes require the enriched test space and trial
!> space to live on identical knot locations in every direction. Call this
!> once after the spaces are initialized, before the time loop starts. The
!> detailed one-dimensional check is delegated to \ref ValidateIGRMMesh.
!
! Input:
! ------
!> @param[in] ads_test
!> Enriched iGRM test-space setup.
!>
!> @param[in] ads_trial
!> Trial-space setup.
!
!---------------------------------------------------------------------------
subroutine ValidateIGRMTimeSchemeSpaces(ads_test, ads_trial)
      use Setup, ONLY: ADS_Setup
      use projection_engine, ONLY: ValidateIGRMMesh
      implicit none
!> @brief Enriched iGRM test-space and trial-space setup structures.
      type(ADS_setup), intent(in) :: ads_test, ads_trial

      call ValidateIGRMMesh(ads_test%Ux, ads_test%p(1), ads_test%n(1), ads_test%nelem(1), &
                            ads_trial%Ux, ads_trial%p(1), ads_trial%n(1), ads_trial%nelem(1))
      call ValidateIGRMMesh(ads_test%Uy, ads_test%p(2), ads_test%n(2), ads_test%nelem(2), &
                            ads_trial%Uy, ads_trial%p(2), ads_trial%n(2), ads_trial%nelem(2))
      call ValidateIGRMMesh(ads_test%Uz, ads_test%p(3), ads_test%n(3), ads_test%nelem(3), &
                            ads_trial%Uz, ads_trial%p(3), ads_trial%n(3), ads_trial%nelem(3))

end subroutine ValidateIGRMTimeSchemeSpaces

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Builds coefficient tables for the 3D Douglas-Gunn scheme.
!>
!> @details
!> The coefficients follow the 3D Douglas-Gunn formulation:
!> the x, y, and z implicit axes use half-step LHS operators, while
!> explicit RHS terms are selected from either the previous state or the
!> active intermediate state through \p rhs_du_state.
!
! Input:
! ------
!> @param[in] tau
!> Time-step length.
!>
!> @param[in] include_transport
!> Optional logical flag selecting transport-term inclusion.
!
! Output:
! -------
!> @param[out] mix
!> Legacy substep matrix-mixing table.
!>
!> @param[out] alpha_step
!> RHS coefficient table for the three substeps.
!>
!> @param[out] lhs_mix
!> Directional LHS matrix-mixing table.
!>
!> @param[out] rhs_du_state
!> RHS derivative-state selector table.
!
!---------------------------------------------------------------------------
subroutine ConfigureDouglasGunn3D(tau, mix, alpha_step, lhs_mix, rhs_du_state, include_transport)
      implicit none
!> @brief Time-step length.
      real(kind=8), intent(in) :: tau
!> @brief Legacy substep matrix-mixing table.
      real(kind=8), intent(out), dimension(4, 3) :: mix
!> @brief RHS coefficient table for the three substeps.
      real(kind=8), intent(out), dimension(7, 3) :: alpha_step
!> @brief Directional LHS matrix-mixing table.
      real(kind=8), intent(out), dimension(4, 3, 3) :: lhs_mix
!> @brief RHS derivative-state selector table.
      integer(kind=4), intent(out), dimension(6, 3) :: rhs_du_state
!> @brief Optional logical flag selecting transport-term inclusion.
      logical, intent(in), optional :: include_transport

      call ConfigureMassTables(mix, lhs_mix)
      alpha_step = 0.d0
      rhs_du_state = 0

      call SetImplicitAxis(lhs_mix, 1, 1, 0.5d0*tau, include_transport)
      call SetImplicitAxis(lhs_mix, 2, 2, 0.5d0*tau, include_transport)
      call SetImplicitAxis(lhs_mix, 3, 3, 0.5d0*tau, include_transport)

      call AddExplicitAxis(alpha_step, rhs_du_state, 1, 1, -0.5d0, 1, include_transport)
      call AddExplicitAxis(alpha_step, rhs_du_state, 2, 1, -1.d0, 1, include_transport)
      call AddExplicitAxis(alpha_step, rhs_du_state, 3, 1, -1.d0, 1, include_transport)
      alpha_step(7, 1) = 1.d0

      call AddExplicitAxis(alpha_step, rhs_du_state, 2, 2, 0.5d0, 1, include_transport)
      call AddExplicitAxis(alpha_step, rhs_du_state, 3, 3, 0.5d0, 1, include_transport)

end subroutine ConfigureDouglasGunn3D

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Builds coefficient tables for a cyclic 3D Peaceman-Rachford step.
!>
!> @details
!> This is the three-dimensional extension used by the code path: each
!> substep treats one physical direction implicitly and the remaining
!> directions explicitly with half-step coefficients. The forcing is split
!> evenly across the three directional substeps.
!
! Input:
! ------
!> @param[in] tau
!> Time-step length.
!>
!> @param[in] include_transport
!> Optional logical flag selecting transport-term inclusion.
!
! Output:
! -------
!> @param[out] mix
!> Legacy substep matrix-mixing table.
!>
!> @param[out] alpha_step
!> RHS coefficient table for the three substeps.
!>
!> @param[out] lhs_mix
!> Directional LHS matrix-mixing table.
!>
!> @param[out] rhs_du_state
!> RHS derivative-state selector table.
!
!---------------------------------------------------------------------------
subroutine ConfigurePeacemanRachford3D(tau, mix, alpha_step, lhs_mix, rhs_du_state, include_transport)
      implicit none
!> @brief Time-step length.
      real(kind=8), intent(in) :: tau
!> @brief Legacy substep matrix-mixing table.
      real(kind=8), intent(out), dimension(4, 3) :: mix
!> @brief RHS coefficient table for the three substeps.
      real(kind=8), intent(out), dimension(7, 3) :: alpha_step
!> @brief Directional LHS matrix-mixing table.
      real(kind=8), intent(out), dimension(4, 3, 3) :: lhs_mix
!> @brief RHS derivative-state selector table.
      integer(kind=4), intent(out), dimension(6, 3) :: rhs_du_state
!> @brief Optional logical flag selecting transport-term inclusion.
      logical, intent(in), optional :: include_transport

      call ConfigureMassTables(mix, lhs_mix)
      alpha_step = 0.d0
      rhs_du_state = 0

      call SetImplicitAxis(lhs_mix, 1, 1, 0.5d0*tau, include_transport)
      call SetImplicitAxis(lhs_mix, 2, 2, 0.5d0*tau, include_transport)
      call SetImplicitAxis(lhs_mix, 3, 3, 0.5d0*tau, include_transport)

      call AddExplicitAxis(alpha_step, rhs_du_state, 2, 1, -0.5d0, 1, include_transport)
      call AddExplicitAxis(alpha_step, rhs_du_state, 3, 1, -0.5d0, 1, include_transport)
      alpha_step(7, 1) = 1.d0/3.d0

      call AddExplicitAxis(alpha_step, rhs_du_state, 1, 2, -0.5d0, 0, include_transport)
      call AddExplicitAxis(alpha_step, rhs_du_state, 3, 2, -0.5d0, 0, include_transport)
      alpha_step(7, 2) = 1.d0/3.d0

      call AddExplicitAxis(alpha_step, rhs_du_state, 1, 3, -0.5d0, 0, include_transport)
      call AddExplicitAxis(alpha_step, rhs_du_state, 2, 3, -0.5d0, 0, include_transport)
      alpha_step(7, 3) = 1.d0/3.d0

end subroutine ConfigurePeacemanRachford3D

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Builds coefficient tables for a split 3D Backward Euler step.
!>
!> @details
!> The scheme is executed through \ref MultiStep. Each substep treats one
!> physical direction implicitly. The forcing term is distributed evenly
!> over the three directional substeps.
!
! Input:
! ------
!> @param[in] tau
!> Time-step length.
!>
!> @param[in] include_transport
!> Optional logical flag selecting transport-term inclusion.
!
! Output:
! -------
!> @param[out] mix
!> Legacy substep matrix-mixing table.
!>
!> @param[out] alpha_step
!> RHS coefficient table for the three substeps.
!>
!> @param[out] lhs_mix
!> Directional LHS matrix-mixing table.
!>
!> @param[out] rhs_du_state
!> RHS derivative-state selector table.
!
!---------------------------------------------------------------------------
subroutine ConfigureBackwardEuler3D(tau, mix, alpha_step, lhs_mix, rhs_du_state, include_transport)
      implicit none
!> @brief Time-step length.
      real(kind=8), intent(in) :: tau
!> @brief Legacy substep matrix-mixing table.
      real(kind=8), intent(out), dimension(4, 3) :: mix
!> @brief RHS coefficient table for the three substeps.
      real(kind=8), intent(out), dimension(7, 3) :: alpha_step
!> @brief Directional LHS matrix-mixing table.
      real(kind=8), intent(out), dimension(4, 3, 3) :: lhs_mix
!> @brief RHS derivative-state selector table.
      integer(kind=4), intent(out), dimension(6, 3) :: rhs_du_state
!> @brief Optional logical flag selecting transport-term inclusion.
      logical, intent(in), optional :: include_transport
!> @brief Directional substep length used for implicit axis operators.
      real(kind=8) :: sub_tau

      call ConfigureMassTables(mix, lhs_mix)
      alpha_step = 0.d0
      rhs_du_state = 0
      sub_tau = 0.5d0*tau

      call SetImplicitAxis(lhs_mix, 1, 1, sub_tau, include_transport)
      call SetImplicitAxis(lhs_mix, 2, 2, sub_tau, include_transport)
      call SetImplicitAxis(lhs_mix, 3, 3, sub_tau, include_transport)
      alpha_step(7, :) = 1.d0/3.d0

end subroutine ConfigureBackwardEuler3D

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Builds legacy mass-only coefficient tables for iGRM multi-step.
!>
!> @details
!> The routine only prepares the coefficients consumed by \ref MultiStep.
!> It does not call \ref Step, because \ref Step is the Forward Euler path
!> and intentionally keeps `alpha_step = 1`.
!>
!> The matrix mixing convention is:
!> - `mix(1,:)` for mass matrices,
!> - `mix(2,:)` for stiffness matrices,
!> - `mix(3,:)` and `mix(4,:)` for first-derivative blocks.
!>
!> The RHS coefficient convention is:
!> - rows 1,3,5 for diffusion-like derivative terms,
!> - rows 2,4,6 for transport-like first-derivative terms,
!> - row 7 for the scalar forcing callback.
!>
!> This routine is retained for older callers. The actual Douglas-Gunn
!> wrapper uses \ref ConfigureDouglasGunn3D so it can pass directional LHS
!> matrices and RHS derivative-state selectors into \ref MultiStep.
!
! Input:
! ------
!> @param[in] tau
!> Time-step length retained in the interface for future schemes.
!
! Output:
! -------
!> @param[out] mix
!> Matrix mixing coefficients for the three iGRM substeps.
!>
!> @param[out] alpha_step
!> RHS coefficients for the three iGRM substeps.
!
!---------------------------------------------------------------------------
subroutine ConfigureDouglasGunn(tau, mix, alpha_step)
      implicit none
!> @brief Time-step length.
      real(kind=8), intent(in) :: tau
!> @brief Matrix mixing coefficients for \ref MultiStep.
      real(kind=8), intent(out), dimension(4, 3) :: mix
!> @brief RHS coefficient table for \ref MultiStep.
      real(kind=8), intent(out), dimension(7, 3) :: alpha_step

      mix = 0.d0
      mix(1, :) = 1.d0
      alpha_step = 1.d0

end subroutine ConfigureDouglasGunn

end module time_scheme
