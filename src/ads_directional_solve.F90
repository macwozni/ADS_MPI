!------------------------------------------------------------------------------
!
! MODULE: ads_directional_solve
!
! DESCRIPTION:
!> @file ads_directional_solve.F90
!> @brief Directional ADS gather/solve/scatter/reorder workflow.
!>
!> @details
!> This module contains the one-direction solve orchestration used by
!> \ref Sub_Step. A directional solve gathers the active tensor slice to
!> the appropriate processor face, builds the one-dimensional operator,
!> calls the local MUMPS solve wrapper, scatters the solution back to the
!> full process grid, and reorders the result for the next directional
!> stage.
!>
!> The implementation deliberately preserves the existing numerical solve
!> order and storage permutations. The refactor only separates the
!> direction-specific orchestration from \ref ADSS so that the time-step
!> code can stay readable.
!>
!> Both standard ADS and iGRM enriched solves are handled here. The iGRM
!> path gathers and scatters the auxiliary test-space buffers alongside
!> the trial-space RHS and therefore shares the same MPI communication
!> conventions as the legacy implementation.
!
!------------------------------------------------------------------------------
module ads_directional_solve

   implicit none

   integer(kind=4), parameter :: MATRIX_SIZE_MISMATCH_ERROR = -1001

   private
   public :: solve_problem

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Solves one directional subproblem, including gather, matrix
!> assembly, linear solve, scatter, and RHS reordering.
!>
!> @details
!> This routine manages one directional one-dimensional solve embedded in
!> the tensor-product ADS workflow. Depending on the active direction and
!> the iGRM flag, it:
!> - selects the proper communicator and directional metadata,
!> - gathers the relevant RHS data onto the solving face of processors,
!> - optionally gathers the enriched test-space contribution,
!> - assembles the operator matrix by calling \ref ComputeMatrix,
!> - solves the resulting linear systems through \ref SolveOneDirection,
!> - scatters the results back to the full process grid,
!> - reorders the RHS arrays for the next directional pass.
!>
!> The routine supports both equal-space and mixed-space solves through
!> the logical flag \p equ, determined from \p direction and the current
!> solving axis \p a.
!
! Input:
! ------
!> @param[in] ads_test
!> Test-space setup structure.
!>
!> @param[in] ads_trial
!> Trial-space setup structure.
!>
!> @param[in] a
!> Direction currently being solved.
!>
!> @param[in] b
!> First transverse direction used in tensor reordering.
!>
!> @param[in] c
!> Second transverse direction used in tensor reordering.
!>
!> @param[in] mixA
!> Mixing coefficients for the diagonal or equal-space operator block.
!>
!> @param[in] mixB
!> Mixing coefficients for the mixed off-diagonal block.
!>
!> @param[in] mixBT
!> Mixing coefficients for the transposed mixed block.
!>
!> @param[in] direction
!> Directional enrichment selector.
!>
!> @param[in] igrm
!> Logical flag indicating whether iGRM mode is active.
!
! Input/Output:
! -------------
!> @param[inout] F
!> Trial-space RHS/solution buffer for the current solve.
!>
!> @param[inout] F2
!> Trial-space output buffer after scatter and reordering.
!>
!> @param[inout] Ft
!> Test-space RHS/solution buffer for the current solve in iGRM mode.
!>
!> @param[inout] Ft2
!> Test-space output buffer after scatter and reordering in iGRM mode.
!
! Output:
! -------
!> @param[out] ierr
!> Returned status code.
!
! Notes:
! ------
!> @note
!> The routine uses direction-specific reorderers `ReorderRHSForX`,
!> `ReorderRHSForY`, and `ReorderRHSForZ`.
!
!---------------------------------------------------------------------------
subroutine solve_problem(ads_test, ads_trial, a, b, c, mixA, mixB, mixBT, direction, igrm, F, F2, Ft, Ft2, ierr)
      use Setup, ONLY: ADS_Setup
      use sparse
      use mumps_solver, ONLY: SolveOneDirection
      use mpi
      use my_mpi, ONLY: Gather, Scatter
      use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, NRPROCX, NRPROCY, NRPROCZ, FillDimVector
      use communicators, ONLY: COMMX, COMMY, COMMZ
      use reorderRHS, ONLY: ReorderRHSForX, ReorderRHSForY, ReorderRHSForZ
      use projection_engine, ONLY: ComputeMatrix
      implicit none
!> @brief Test and trial setup structures.
      type(ADS_setup), intent(in) :: ads_test, ads_trial
!> @brief Solving direction and transverse-direction permutation.
      integer(kind=4), intent(in) :: a, b, c
!> @brief Mixing vectors used for matrix assembly.
      real(kind=8), dimension(4), intent(in) :: mixA, mixB, mixBT
!> @brief Directional enrichment selector.
      integer(kind=4), dimension(3), intent(in) :: direction
!> @brief Flag indicating whether iGRM mode is active.
      logical, intent(in) :: igrm
!> @brief Trial-space buffers for input and reordered output.
      real(kind=8), allocatable, dimension(:, :) :: F, F2
!> @brief Test-space buffers for input and reordered output.
      real(kind=8), allocatable, dimension(:, :) :: Ft, Ft2
!> @brief Returned status code.
      integer(kind=4), intent(out) :: ierr
      integer(kind=4), dimension(:), allocatable :: dimensions_test, dimensions_trial ! Size of slices of domain in each dimension
      integer(kind=4), allocatable, dimension(:) :: shifts_test, shifts_trial
      integer(kind=4) :: comm
      real(kind=8), allocatable, dimension(:) :: U_trial, U_test
      integer(kind=4) :: myrankdim ! Integer coordinates of processor along X, Y or Z
      type(sparse_matrix), pointer :: sprsmtrx
      real(kind=8), allocatable, dimension(:, :) :: Fs ! F-solve
      real(kind=8), allocatable, dimension(:, :) :: F_out, F2_out ! F-trial
      real(kind=8), allocatable, dimension(:, :) :: Ft_out, Ft2_out ! F-test
      integer(kind=4), dimension(3) :: ibeg, iend
      integer(kind=4) :: system_n
      integer(kind=4) :: trial_owned_a, mixed_owned_a
      integer(kind=4) :: mpi_ierr, solver_ierr, synchronized_ierr
      ! real(kind=8) :: time1, time2
      logical :: equ

      ierr = 0
      nullify(sprsmtrx)

      !  we have identical test and trial spaces if equ=true in given direction
      equ = .TRUE.
      !  if we have enriched test space in given direction
      if (direction(a) .EQ. 1) equ = .FALSE.

      !  set proper paremeters depending on which direction we solve
      !  we solve in x directon
      if (a .EQ. 1) then
            comm = COMMX
            myrankdim = MYRANKX
            U_trial = ads_trial%ux
            U_test = ads_test%ux
            shifts_trial = ads_trial%shiftsX
            shifts_test = ads_test%shiftsX
            dimensions_trial = ads_trial%dimensionsX
            dimensions_test = ads_test%dimensionsX
            ! prepare dimensions vectors
            call FillDimVector(dimensions_test, shifts_test,&
            (1-direction(1))*ads_trial%nrcpp(1) + direction(1)*ads_test%nrcpp(1),&
            (direction(2)*ads_test%s(2)+(1-direction(2))*ads_trial%s(2))*&
            (direction(3)*ads_test%s(3)+(1-direction(3))*ads_trial%s(3)),&
            (direction(1)*ads_test%n(1)+(1-direction(1))*ads_trial%n(1)),&
            NRPROCX)
      !  we solve in y directon
      else if (a .EQ. 2) then
            comm = COMMY
            myrankdim = MYRANKY
            U_trial = ads_trial%uy
            U_test = ads_test%uy
            shifts_trial = ads_trial%shiftsY
            shifts_test = ads_test%shiftsY
            dimensions_trial = ads_trial%dimensionsY
            dimensions_test = ads_test%dimensionsY
            ! prepare dimensions vectors
            call FillDimVector(dimensions_test, shifts_test,&
            (1-direction(2))*ads_trial%nrcpp(2) + direction(2)*ads_test%nrcpp(2),&
            (direction(1)*ads_test%s(1)+(1-direction(1))*ads_trial%s(1))*&
            (direction(3)*ads_test%s(3)+(1-direction(3))*ads_trial%s(3)),&
            (direction(2)*ads_test%n(2)+(1-direction(2))*ads_trial%n(2)),&
            NRPROCY)
      !  we solve in z directon
      else ! a.EQ.3
            comm = COMMZ
            myrankdim = MYRANKZ
            U_trial = ads_trial%uz
            U_test = ads_test%uz
            shifts_trial = ads_trial%shiftsZ
            shifts_test = ads_test%shiftsZ
            dimensions_trial = ads_trial%dimensionsZ
            dimensions_test = ads_test%dimensionsZ
            ! prepare dimensions vectors
            call FillDimVector(dimensions_test, shifts_test,&
            (1-direction(3))*ads_trial%nrcpp(3) + direction(3)*ads_test%nrcpp(3),&
            (direction(1)*ads_test%s(1)+(1-direction(1))*ads_trial%s(1))*&
            (direction(2)*ads_test%s(2)+(1-direction(2))*ads_trial%s(2)),&
            (direction(3)*ads_test%n(3)+(1-direction(3))*ads_trial%n(3)),&
            NRPROCZ)
      end if

      ibeg = ads_trial%ibeg
      iend = ads_trial%iend
      if (direction(1) .EQ. 1) then
            ibeg(1) = ads_test%ibeg(1)
            iend(1) = ads_test%iend(1)
      endif
      if (direction(2) .EQ. 1) then
            ibeg(2) = ads_test%ibeg(2)
            iend(2) = ads_test%iend(2)
      endif
      if (direction(3) .EQ. 1) then
            ibeg(3) = ads_test%ibeg(3)
            iend(3) = ads_test%iend(3)
      endif
      trial_owned_a = ads_trial%s(a)
      mixed_owned_a = (1 - direction(a))*ads_trial%s(a) + direction(a)*ads_test%s(a)

      ! This collective also provides the synchronization previously supplied
      ! by the world barrier, while making every rank start with one status.
      call SynchronizeFirstError(ierr, synchronized_ierr)
      ierr = synchronized_ierr
      if (ierr /= 0) return

#ifdef IINFO
      write (*, *) PRINTRANK, a, 'a) GATHER'
      call SynchronizeFirstError(ierr, synchronized_ierr)
      ierr = synchronized_ierr
      if (ierr /= 0) return
#endif

      !  allocate result buffer
      allocate (Fs((ads_trial%n(a) + 1 + direction(a)*(ads_test%n(a) + 1)), &
                        (ads_trial%s(b) + direction(b)*ads_test%s(b))*(ads_trial%s(c) + direction(c)*ads_test%s(c))))
      allocate (F_out((ads_trial%n(a) + 1), &
                  (ads_trial%s(b)*ads_trial%s(c))))
#ifdef PERFORMANCE
      time1 = MPI_Wtime()
#endif
      !  gather onto the face of processors
      call Gather(F, F_out, ads_trial%n(a), &
                  ads_trial%s(a), &
                  ads_trial%s(b)*ads_trial%s(c), &
            dimensions_trial, shifts_trial, comm, mpi_ierr)
      call RecordFirstError(ierr, mpi_ierr)
      call SynchronizeFirstError(ierr, synchronized_ierr)
      ierr = synchronized_ierr
      if (ierr /= 0) go to 900
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write (*, *) "Gather", a, " : ", time2 - time1
#endif

      if (igrm) then
      !  allocate result buffer
            allocate (Ft_out(((1 - direction(a))*ads_trial%n(a) + direction(a)*ads_test%n(a) + 1), &
                              ((1 - direction(b))*ads_trial%s(b) + direction(b)*ads_test%s(b))* &
                        ((1 - direction(c))*ads_trial%s(c) + direction(c)*ads_test%s(c))))
#ifdef PERFORMANCE
            time1 = MPI_Wtime()
#endif
      !  gather onto the face of processors
            call Gather(Ft, Ft_out, (1 - direction(a))*ads_trial%n(a) + direction(a)*ads_test%n(a), &
                        (1 - direction(a))*ads_trial%s(a) + direction(a)*ads_test%s(a), &
                        ((1 - direction(b))*ads_trial%s(b) + direction(b)*ads_test%s(b))* &
                        ((1 - direction(c))*ads_trial%s(c) + direction(c)*ads_test%s(c)), &
                        dimensions_test, shifts_test, comm, mpi_ierr)
            call RecordFirstError(ierr, mpi_ierr)
            call SynchronizeFirstError(ierr, synchronized_ierr)
            ierr = synchronized_ierr
            if (ierr /= 0) go to 900
#ifdef PERFORMANCE
            time2 = MPI_Wtime()
            write (*, *) "Gather", a, " : ", time2 - time1
#endif
            if (equ) then
            Fs(:, 1:ads_trial%s(b)*ads_trial%s(c)) = F_out
            Fs(:, (ads_trial%s(b)*ads_trial%s(c) + 1): &
                  (ads_trial%s(b) + direction(b)*ads_test%s(b))*(ads_trial%s(c) + direction(c)*ads_test%s(c))) = Ft_out
            else
            Fs(1:ads_trial%n(a) + 1, :) = F_out
            Fs(ads_trial%n(a) + 2:ads_trial%n(a) + 1 + direction(a)*(ads_test%n(a) + 1), :) = Ft_out
            end if
      else
            Fs = F_out
      end if

      !  performed only on face of processors
      solver_ierr = 0
      if (myrankdim == 0) then
#ifdef IINFO
            write (*, *) PRINTRANK, a, 'b) SOLVE THE ', a, ' PROBLEM'
#endif

#ifdef PERFORMANCE
            time1 = MPI_Wtime()
#endif
      !     compute LHS matrix
            call ComputeMatrix(U_test, ads_test%p(a), ads_test%n(a), ads_test%nelem(a), &
                        U_trial, ads_trial%p(a), ads_trial%n(a), ads_trial%nelem(a), &
                        mixA, mixB, mixBT, equ, sprsmtrx)
#ifdef PERFORMANCE
            time2 = MPI_Wtime()
            write (*, *) "Mass matrix", a, ": ", time2 - time1
            time1 = MPI_Wtime()
#endif
      !     perform real solver
            system_n = size(Fs, 1) - 1
            if (sprsmtrx%x /= size(Fs, 1) .or. sprsmtrx%y /= size(Fs, 1)) then
                  write (*, *) 'matrix/Fs size mismatch', sprsmtrx%x, sprsmtrx%y, size(Fs, 1)
                  solver_ierr = MATRIX_SIZE_MISMATCH_ERROR
            else
                  call SolveOneDirection(Fs, size(Fs, 2), system_n, system_n, sprsmtrx, solver_ierr)
            end if
      !     clean buffers
            call clear_matrix(sprsmtrx)
#ifdef PERFORMANCE
            time2 = MPI_Wtime()
            write (*, *) "Solve ", a, ": ", time2 - time1
#endif
      end if

      ! Only the processor face executes MUMPS. Propagate its status before
      ! any rank enters the scatter phase so all ranks take the same path.
      call SynchronizeFirstError(solver_ierr, synchronized_ierr)
      call RecordFirstError(ierr, synchronized_ierr)
      if (ierr /= 0) go to 900

#ifdef IINFO
      write (*, *) PRINTRANK, a, 'c) SCATTER'
#endif
      if (igrm) then
            if (equ) then
            F_out = Fs(:, 1:ads_trial%s(b)*ads_trial%s(c))
            Ft_out = Fs(:, (ads_trial%s(b)*ads_trial%s(c) + 1): &
                        (ads_trial%s(b) + direction(b)*ads_test%s(b))*(ads_trial%s(c) + direction(c)*ads_test%s(c)))
            else
            F_out = Fs(1:ads_trial%n(a) + 1, :)
            Ft_out = Fs(ads_trial%n(a) + 2:ads_trial%n(a) + 1 + direction(a)*(ads_test%n(a) + 1), :)
            end if
      !  allocate buffers
            allocate (Ft2_out(mixed_owned_a, &
                        ((1 - direction(b))*ads_trial%s(b) + direction(b)*ads_test%s(b))* &
                        ((1 - direction(c))*ads_trial%s(c) + direction(c)*ads_test%s(c))))
#ifdef PERFORMANCE
            time1 = MPI_Wtime()
#endif
      !  scatter back onto the cube of processors
            call Scatter(Ft_out, Ft2_out, (1 - direction(a))*ads_trial%n(a) + direction(a)*ads_test%n(a), &
                        mixed_owned_a, &
                        ((1 - direction(b))*ads_trial%s(b) + direction(b)*ads_test%s(b))* &
                        ((1 - direction(c))*ads_trial%s(c) + direction(c)*ads_test%s(c)), &
                  dimensions_test, shifts_test, comm, mpi_ierr)
            call RecordFirstError(ierr, mpi_ierr)
            call SynchronizeFirstError(ierr, synchronized_ierr)
            ierr = synchronized_ierr
            if (ierr /= 0) go to 900
#ifdef PERFORMANCE
            time2 = MPI_Wtime()
            write (*, *) "Scatter ", a, ": ", time2 - time1
#endif
      else
            F_out = Fs
      end if

      !  allocate buffers
      allocate (F2_out(trial_owned_a, &
                        ads_trial%s(b)*ads_trial%s(c)))
#ifdef PERFORMANCE
      time1 = MPI_Wtime()
#endif
      !  scatter back onto the cube of processors
      call Scatter(F_out, F2_out, ads_trial%n(a), &
                        trial_owned_a, &
                  ads_trial%s(b)*ads_trial%s(c), &
                  dimensions_trial, shifts_trial, comm, mpi_ierr)
      call RecordFirstError(ierr, mpi_ierr)
      call SynchronizeFirstError(ierr, synchronized_ierr)
      ierr = synchronized_ierr
      if (ierr /= 0) go to 900
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write (*, *) "Scatter ", a, ": ", time2 - time1
#endif
      !  cleanup
      if (allocated(F_out)) deallocate (F_out)
      if (allocated(Ft_out)) deallocate (Ft_out)

#ifdef IINFO
      write (*, *) PRINTRANK, a, 'd) REORDER'
#endif
      ! Reorder right hand sides
      if (a .EQ. 1) call ReorderRHSForY(ads_trial%ibeg, ads_trial%iend, F2_out, F2)
      if (a .EQ. 2) call ReorderRHSForZ(ads_trial%ibeg, ads_trial%iend, F2_out, F2)
      if (a .EQ. 3) call ReorderRHSForX(ads_trial%ibeg, ads_trial%iend, F2_out, F2)

      if (igrm) then
            if (a .EQ. 1) call ReorderRHSForY(ibeg, iend, Ft2_out, Ft2)
            if (a .EQ. 2) call ReorderRHSForZ(ibeg, iend, Ft2_out, Ft2)
            if (a .EQ. 3) call ReorderRHSForX(ibeg, iend, Ft2_out, Ft2)
      end if
      !  cleanup
      if (allocated(F2_out)) deallocate (F2_out)
      if (allocated(Ft2_out)) deallocate (Ft2_out)

#ifdef IPRINT
      write (*, *) PRINTRANK, 'after ReorderRHS'
      write (*, *) PRINTRANK, 'F:'
      do i = 1, ads_trial%s(1)
            write (*, *) PRINTRANK, i, 'row=', F2(i, 1:ads_trial%s(2)*ads_trial%s(3))
      end do
#endif

900   continue

      if (associated(sprsmtrx)) call clear_matrix(sprsmtrx)

end subroutine solve_problem

!---------------------------------------------------------------------------
!> @brief Preserves the first nonzero status produced by a staged workflow.
!---------------------------------------------------------------------------
subroutine RecordFirstError(first_error, candidate_error)
      implicit none
      integer(kind=4), intent(inout) :: first_error
      integer(kind=4), intent(in) :: candidate_error

      if (first_error == 0 .and. candidate_error /= 0) then
            first_error = candidate_error
      end if
end subroutine RecordFirstError

!---------------------------------------------------------------------------
!> @brief Makes a local error status identical on every world rank.
!>
!> The lowest failing world rank supplies the returned code. With no local
!> failures, the all-reduction acts as the synchronization point required by
!> the directional gather/solve/scatter workflow.
!---------------------------------------------------------------------------
subroutine SynchronizeFirstError(local_error, global_error)
      use mpi
      use parallelism, ONLY: MYRANK
      implicit none
      integer(kind=4), intent(in) :: local_error
      integer(kind=4), intent(out) :: global_error
      integer(kind=4) :: candidate(2), selected(2)
      integer(kind=4) :: mpi_ierr, abort_ierr

      candidate(1) = 1
      if (local_error /= 0) candidate(1) = 0
      candidate(2) = MYRANK

      call MPI_Allreduce(candidate, selected, 1, MPI_2INTEGER, MPI_MINLOC, &
                         MPI_COMM_WORLD, mpi_ierr)
      if (mpi_ierr /= MPI_SUCCESS) then
            global_error = local_error
            if (global_error == 0) global_error = mpi_ierr
            ! A failed world collective cannot be recovered by another
            ! collective without risking divergent rank paths.
            call MPI_Abort(MPI_COMM_WORLD, mpi_ierr, abort_ierr)
            return
      end if

      if (selected(1) /= 0) then
            global_error = 0
            return
      end if

      global_error = 0
      if (MYRANK == selected(2)) global_error = local_error
      call MPI_Bcast(global_error, 1, MPI_INTEGER, selected(2), &
                     MPI_COMM_WORLD, mpi_ierr)
      if (mpi_ierr /= MPI_SUCCESS .and. global_error == 0) then
            global_error = mpi_ierr
      end if
      if (mpi_ierr /= MPI_SUCCESS) then
            call MPI_Abort(MPI_COMM_WORLD, mpi_ierr, abort_ierr)
      end if
end subroutine SynchronizeFirstError

end module ads_directional_solve
