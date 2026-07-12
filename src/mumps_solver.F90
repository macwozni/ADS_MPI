!------------------------------------------------------------------------------
!
! MODULE: mumps_solver
!
! DESCRIPTION:
!> @file mumps_solver.F90
!> @brief MUMPS-based local one-dimensional solve wrapper.
!>
!> @details
!> This module isolates the direct linear solve stage used by directional
!> ADS subproblems. The sparse one-dimensional operator is converted to
!> MUMPS triplet storage and solved for all right-hand-side columns owned
!> by the calling processor face.
!>
!> The solver is configured on `MPI_COMM_SELF`, matching the current
!> sequential-MUMPS linking model used by each MPI rank. The surrounding
!> ADS computation may run with many MPI ranks, but each rank calls MUMPS
!> as a local sequential direct solver.
!>
!> Keeping this wrapper separate makes the MUMPS-specific control
!> parameters, error checks, and triplet conversion boundary explicit
!> without mixing them into the time-stepping code.
!
!------------------------------------------------------------------------------
module mumps_solver

   implicit none

   private
   public :: SolveOneDirection

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Solves a one-dimensional sparse linear system with multiple
!> right-hand sides using MUMPS.
!>
!> @details
!> This routine converts the sparse matrix \p sprsmtrx to MUMPS format,
!> configures selected MUMPS control parameters, performs analysis and
!> factorization, and then solves the system separately for each column
!> of the right-hand-side array \p RHS.
!>
!> The solve is executed in `MPI_COMM_SELF`, which means each calling
!> process performs its own local/direct solve once the data have been
!> gathered to the solving face.
!
! Input:
! ------
!> @param[inout] RHS
!> Array of right-hand sides overwritten by the computed solutions.
!>
!> @param[in] eqnum
!> Number of right-hand sides.
!>
!> @param[in] n
!> Number of unknowns minus one.
!>
!> @param[in] p
!> Polynomial degree or bandwidth-related parameter passed through the
!> interface.
!>
!> @param[in] sprsmtrx
!> Sparse matrix describing the one-dimensional operator.
!
! Notes:
! ------
!> @note
!> The argument \p p is accepted by the interface but is not used
!> directly in the current implementation.
!
!> @warning
!> The routine stops execution if any MUMPS phase returns a negative
!> primary or global error code.
!
!---------------------------------------------------------------------------
subroutine SolveOneDirection(RHS, eqnum, n, p, sprsmtrx)
use sparse
      use mpi
      implicit none
      include 'dmumps_struc.h'
!> @brief Right-hand-side matrix overwritten by the solution.
      real(kind=8) :: RHS(:, :)
!> @brief System size minus one and auxiliary interface parameter.
      integer :: n, p
!> @brief Number of right-hand sides.
      integer(kind=4) :: eqnum
      integer(kind=4) :: i!, iret
!> @brief Sparse matrix of the one-dimensional operator.
      type(sparse_matrix), pointer, intent(in) :: sprsmtrx
      type(dmumps_struc) :: mumps_par

      !  initialize MUMPS
      mumps_par%sym = 0
      mumps_par%comm = MPI_COMM_SELF
      mumps_par%job = -1
      mumps_par%par = 1
      mumps_par%N = n + 1
      call dmumps(mumps_par)
      call CheckMumpsStatus('initialization')

      !  convert LHS and RHS to MUMPS format
      call to_mumps_format(sprsmtrx, mumps_par)
      allocate (mumps_par%rhs(mumps_par%n))

!  set MUMPS parameters
!  error output stream (non-positive to suppress)
      mumps_par%icntl(1) = 1 !1
!  diagnostic, statistics and warnings
      mumps_par%icntl(2) = 0! 1 !1
!  global information
      mumps_par%icntl(3) = 0!6 !6
!  printing level
      mumps_par%icntl(4) = 0!3 !3
!  input matrix in assembled or element format
      mumps_par%icntl(5) = 0
!  column permutation for zero-free diagonal (automatic)
!  mumps_par%icntl(6)  = 7
!  pivot order (automatic)
      mumps_par%icntl(7) = 7 !1 enforce ordering, 5 metis, 0 AMD, 7 auto
!  scaling (automatic)
!  mumps_par%icntl(8)  = 7
!  no transpose
!  mumps_par%icntl(9)  = 1
!  max steps for iterative refinement
!  mumps_par%icntl(10) = 0
!  statistics info
      mumps_par%icntl(11) = 2
!  controls parallelism
      mumps_par%icntl(12) = 0
!  use ScaLAPACK for root node
      mumps_par%icntl(13) = 1 !0 use 1 do not use
!  percentage increase in estimated workspace
      mumps_par%icntl(14) = 50
!  matrix distribution for assembled input
      mumps_par%icntl(18) = 0 !distributed
!  nonzero for Schur complement
      mumps_par%icntl(19) = 0
!  distribution of RHS (centralized on host)
      mumps_par%icntl(20) = 0
!  mumps_par%icntl(32) = 1

!  start MUMPS
      mumps_par%job = 1
      call dmumps(mumps_par)
      call CheckMumpsStatus('analysis')

      mumps_par%job = 2
      call dmumps(mumps_par)
      call CheckMumpsStatus('factorization')

      do i = 1, eqnum
            mumps_par%rhs(1:n + 1) = rhs(1:n + 1, i)
            mumps_par%job = 3
            call dmumps(mumps_par)
            call CheckMumpsStatus('solve')
            rhs(1:n + 1, i) = mumps_par%rhs(1:n + 1)
      end do

!  stop MUMPS
      mumps_par%job = -2
      call dmumps(mumps_par)
      call CheckMumpsStatus('finalization')

      deallocate (mumps_par%irn)
      deallocate (mumps_par%jcn)
      deallocate (mumps_par%a)
      deallocate (mumps_par%rhs)

contains

      subroutine CheckMumpsStatus(phase)
            implicit none
            character(*), intent(in) :: phase

            if (mumps_par%info(1) .lt. 0 .or. mumps_par%infog(1) .lt. 0) then
                  write (*, *) 'MUMPS failed during ', trim(phase)
                  write (*, *) 'mumps_par%job=', mumps_par%job
                  write (*, *) 'mumps_par%info=', mumps_par%info
                  write (*, *) 'mumps_par%infog=', mumps_par%infog
                  stop 1
            end if
      end subroutine CheckMumpsStatus

end subroutine SolveOneDirection

end module mumps_solver
