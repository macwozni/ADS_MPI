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
! Output:
! -------
!> @param[out] ierr
!> Zero after a successful solve, otherwise the first negative MUMPS
!> status reported by initialization, analysis, factorization, solve, or
!> finalization. On failure, \p RHS must not be used as a complete solution.
!> The argument is optional for source compatibility. If an older caller
!> omits it and MUMPS fails, the complete MPI job is aborted instead of
!> stopping only the calling rank.
!
!---------------------------------------------------------------------------
subroutine SolveOneDirection(RHS, eqnum, n, p, sprsmtrx, ierr)
      use sparse
      use mpi
      implicit none
      include 'dmumps_struc.h'
!> @brief Right-hand-side matrix overwritten by the solution.
      real(kind=8), intent(inout) :: RHS(:, :)
!> @brief System size minus one and auxiliary interface parameter.
      integer, intent(in) :: n, p
!> @brief Number of right-hand sides.
      integer(kind=4), intent(in) :: eqnum
!> @brief Sparse matrix of the one-dimensional operator.
      type(sparse_matrix), pointer, intent(in) :: sprsmtrx
!> @brief First MUMPS error code, or zero after a successful solve.
      integer(kind=4), intent(out), optional :: ierr
      integer(kind=4) :: i
      integer(kind=4) :: phase_ierr, first_error, abort_ierr
      logical :: mumps_instance_started
      type(dmumps_struc) :: mumps_par

      first_error = 0
      mumps_instance_started = .false.
      nullify(mumps_par%irn, mumps_par%jcn, mumps_par%a, mumps_par%rhs)
      nullify(mumps_par%intr_encoding)

      !  initialize MUMPS
      mumps_par%sym = 0
      mumps_par%comm = MPI_COMM_SELF
      mumps_par%job = -1
      mumps_par%par = 1
      mumps_par%N = n + 1
      call dmumps(mumps_par)
      ! MUMPS creates INTR_ENCODING before entering the initialization
      ! driver.  Its presence means that JOB=-2 is valid even when the
      ! initialization driver itself reported an error.
      mumps_instance_started = associated(mumps_par%intr_encoding)
      call CheckMumpsStatus('initialization', phase_ierr)
      call RecordFirstError(first_error, phase_ierr)
      if (first_error /= 0) go to 900

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
      call CheckMumpsStatus('analysis', phase_ierr)
      call RecordFirstError(first_error, phase_ierr)
      if (first_error /= 0) go to 900

      mumps_par%job = 2
      call dmumps(mumps_par)
      call CheckMumpsStatus('factorization', phase_ierr)
      call RecordFirstError(first_error, phase_ierr)
      if (first_error /= 0) go to 900

      do i = 1, eqnum
            mumps_par%rhs(1:n + 1) = rhs(1:n + 1, i)
            mumps_par%job = 3
            call dmumps(mumps_par)
            call CheckMumpsStatus('solve', phase_ierr)
            call RecordFirstError(first_error, phase_ierr)
            if (first_error /= 0) go to 900
            rhs(1:n + 1, i) = mumps_par%rhs(1:n + 1)
      end do

900   continue

      ! Finalize every instance which MUMPS actually started, including a
      ! partially initialized instance. A successful finalization must not
      ! overwrite an earlier phase error.
      if (mumps_instance_started) then
            mumps_par%job = -2
            call dmumps(mumps_par)
            call CheckMumpsStatus('finalization', phase_ierr)
            call RecordFirstError(first_error, phase_ierr)
      end if

      if (associated(mumps_par%irn)) deallocate (mumps_par%irn)
      if (associated(mumps_par%jcn)) deallocate (mumps_par%jcn)
      if (associated(mumps_par%a)) deallocate (mumps_par%a)
      if (associated(mumps_par%rhs)) deallocate (mumps_par%rhs)

      if (present(ierr)) then
            ierr = first_error
      else if (first_error /= 0) then
            call MPI_Abort(MPI_COMM_WORLD, first_error, abort_ierr)
            stop 5
      end if

contains

      subroutine CheckMumpsStatus(phase, status)
            implicit none
            character(*), intent(in) :: phase
            integer(kind=4), intent(out) :: status

            status = 0
            if (mumps_par%infog(1) < 0) then
                  status = mumps_par%infog(1)
            else if (mumps_par%info(1) < 0) then
                  status = mumps_par%info(1)
            end if

            if (status /= 0) then
                  write (*, *) 'MUMPS failed during ', trim(phase)
                  write (*, *) 'mumps_par%job=', mumps_par%job
                  write (*, *) 'mumps_par%info(1:2)=', mumps_par%info(1:2)
                  write (*, *) 'mumps_par%infog(1:2)=', mumps_par%infog(1:2)
            end if
      end subroutine CheckMumpsStatus

      subroutine RecordFirstError(first_error, candidate_error)
            implicit none
            integer(kind=4), intent(inout) :: first_error
            integer(kind=4), intent(in) :: candidate_error

            if (first_error == 0 .and. candidate_error /= 0) then
                  first_error = candidate_error
            end if
      end subroutine RecordFirstError

end subroutine SolveOneDirection

end module mumps_solver
