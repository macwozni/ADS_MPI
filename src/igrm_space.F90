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
!> The comparison ignores knot multiplicities and only checks the sequence of
!> distinct knot locations. This matches the iGRM requirement that test and
!> trial spaces share the same element partition, while still allowing
!> repeated knots used to lower continuity.
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
