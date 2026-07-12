!------------------------------------------------------------------------------
!
! MODULE: igrm_space
!
! DESCRIPTION:
!> @file igrm_space.F90
!> @brief iGRM test/trial spline-space compatibility checks.
!>
!> @details
!> This module groups checks related to the pair of spline spaces used by
!> the iGRM path. The central requirement is that the test and trial
!> spaces share the same geometric mesh while the test space is enriched
!> relative to the trial space.
!>
!> Knot multiplicity is not treated as a separate mesh location. This is
!> important for iGRM setups where the number of degrees of freedom may be
!> increased by lowering spline continuity rather than only by raising the
!> polynomial degree.
!>
!> The public validation routine is intended to be called during problem
!> setup, before the time loop starts, so the compatibility checks are not
!> repeated at every time step.
!
!------------------------------------------------------------------------------
module igrm_space

   implicit none

   private
   public :: ValidateIGRMMesh

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Checks whether two spline spaces use the same geometric mesh.
!>
!> @details
!> The comparison walks through both knot vectors and compares only the
!> ordered sequence of distinct knot locations. Consecutive knots whose
!> values differ by less than the local tolerance are treated as one
!> location. This makes the routine insensitive to multiplicity changes
!> introduced to lower spline continuity.
!>
!> This distinction matters for iGRM: the test and trial spaces must
!> describe the same geometric element partition, but the test space may
!> still be enriched by changing continuity and therefore by repeating
!> knots. Such repeated knots increase the number of basis functions
!> without creating additional element boundaries.
!>
!> The routine returns `.true.` only when both knot vectors contain the
!> same distinct locations in the same order and both vectors are
!> exhausted at the same time.
!
! Input:
! ------
!> @param[in] U1
!> First knot vector.
!>
!> @param[in] p1
!> Polynomial degree associated with \p U1.
!>
!> @param[in] n1
!> Number of control points minus one for the first space.
!>
!> @param[in] U2
!> Second knot vector.
!>
!> @param[in] p2
!> Polynomial degree associated with \p U2.
!>
!> @param[in] n2
!> Number of control points minus one for the second space.
!
! Output:
! -------
!> @return
!> `.true.` when both vectors define the same unique knot-location
!> sequence; `.false.` otherwise.
!
! Notes:
! ------
!> @note
!> The routine is intentionally local and deterministic. It performs no
!> MPI communication and assumes that both compared vectors are already
!> available in the calling rank's setup data.
!
!> @warning
!> The comparison tolerance is fixed at \f$10^{-12}\f$. Knot vectors whose
!> distinct locations are closer than this tolerance may be interpreted as
!> having the same geometric location.
!
!---------------------------------------------------------------------------
logical function SameKnotLocations(U1, p1, n1, U2, p2, n2) result(same)
   implicit none
!> @brief Polynomial degrees and basis sizes minus one.
   integer(kind=4), intent(in) :: p1, n1, p2, n2
!> @brief Knot vectors to compare.
   real(kind=8), intent(in) :: U1(0:n1 + p1 + 1)
   real(kind=8), intent(in) :: U2(0:n2 + p2 + 1)
   real(kind=8), parameter :: knot_tol = 1.d-12
   integer(kind=4) :: i1, i2, last1, last2
   real(kind=8) :: k1, k2

   same = .FALSE.
   i1 = 0
   i2 = 0
   last1 = n1 + p1 + 1
   last2 = n2 + p2 + 1

   do while (i1 <= last1 .and. i2 <= last2)
      k1 = U1(i1)
      k2 = U2(i2)
      if (abs(k1 - k2) > knot_tol) return

      do
         if (i1 > last1) exit
         if (abs(U1(i1) - k1) > knot_tol) exit
         i1 = i1 + 1
      end do
      do
         if (i2 > last2) exit
         if (abs(U2(i2) - k2) > knot_tol) exit
         i2 = i2 + 1
      end do
   end do

   same = (i1 > last1 .and. i2 > last2)

end function SameKnotLocations

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Validates the iGRM mixed-space assumptions before matrix assembly.
!>
!> @details
!> This procedure checks the compatibility conditions between one
!> test-space direction and the corresponding trial-space direction used
!> by the iGRM path. It is intended to be called once while constructing
!> the problem setup, before operator assembly and before entering the
!> time loop.
!>
!> The validation currently enforces three conditions:
!> - the test polynomial degree must be greater than the trial degree,
!> - the test and trial knot vectors must contain the same distinct knot
!>   locations,
!> - the stored element counts must agree.
!>
!> The knot-location check is delegated to \ref SameKnotLocations, so
!> knot multiplicities are ignored when comparing the geometric mesh.
!> This prevents repeated knots used for continuity control from being
!> mistaken for a different element partition.
!>
!> If any condition fails, the routine writes a diagnostic to
!> `ERROR_UNIT` and stops execution with exit code `1`. This fail-fast
!> behaviour avoids assembling mixed iGRM operators for incompatible
!> spaces.
!
! Input:
! ------
!> @param[in] U_test
!> Test-space knot vector for one parametric direction.
!>
!> @param[in] p_test
!> Test-space polynomial degree.
!>
!> @param[in] n_test
!> Number of test-space control points minus one.
!>
!> @param[in] nelem_test
!> Number of geometric elements recorded for the test space.
!>
!> @param[in] U_trial
!> Trial-space knot vector for the same parametric direction.
!>
!> @param[in] p_trial
!> Trial-space polynomial degree.
!>
!> @param[in] n_trial
!> Number of trial-space control points minus one.
!>
!> @param[in] nelem_trial
!> Number of geometric elements recorded for the trial space.
!
! Notes:
! ------
!> @note
!> The routine validates one direction at a time. A complete 3D iGRM
!> setup should call it independently for the x, y, and z directions.
!
!> @note
!> The procedure checks compatibility metadata only. It does not allocate
!> basis tables, assemble operators, or inspect the numerical contents of
!> solution vectors.
!
!> @warning
!> The current degree check is stricter than a general "test space richer
!> than trial space" condition because it requires \p p_test > \p p_trial.
!> Setups enriched purely through reduced continuity should be handled by
!> a future validation rule based directly on degrees of freedom and
!> unique knot locations.
!
!---------------------------------------------------------------------------
subroutine ValidateIGRMMesh(U_test, p_test, n_test, nelem_test, &
                            U_trial, p_trial, n_trial, nelem_trial)
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
   implicit none
!> @brief Test-space metadata.
   integer(kind=4), intent(in) :: p_test, n_test, nelem_test
!> @brief Trial-space metadata.
   integer(kind=4), intent(in) :: p_trial, n_trial, nelem_trial
!> @brief Test and trial knot vectors.
   real(kind=8), intent(in) :: U_test(0:n_test + p_test + 1)
   real(kind=8), intent(in) :: U_trial(0:n_trial + p_trial + 1)

   if (p_test <= p_trial) then
      write (ERROR_UNIT, *) 'iGRM requires test degree greater than trial degree:', &
         p_test, p_trial
      stop 1
   end if

   if (.not. SameKnotLocations(U_test, p_test, n_test, U_trial, p_trial, n_trial)) then
      write (ERROR_UNIT, *) 'iGRM requires identical test/trial knot locations'
      stop 1
   end if

   if (nelem_test /= nelem_trial) then
      write (ERROR_UNIT, *) 'iGRM mesh metadata mismatch:', nelem_test, nelem_trial
      stop 1
   end if

end subroutine ValidateIGRMMesh

end module igrm_space
