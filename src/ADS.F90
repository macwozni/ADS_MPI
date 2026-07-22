!------------------------------------------------------------------------------
!
! MODULE: ADSS
!
! DESCRIPTION:
!> @file ADS.F90
!> @brief High-level ADS workflow and compatibility API module.
!>
!> @details
!> This module owns the time-stepping orchestration routines used by the
!> ADS workflow and re-exports the historical public entry points that
!> were split into focused implementation modules.
!>
!> The procedures implemented directly here are:
!> - \ref Step, the single-step Forward Euler path,
!> - \ref MultiStep, the three-substep directional ADS/iGRM path,
!> - \ref Sub_Step, the shared per-substep orchestration layer.
!>
!> Setup, cleanup, directional solve details, buffer normalization, MUMPS
!> interaction, and solution output now live in separate modules and are
!> re-exported here for compatibility with existing problem drivers.
!>
!> The Forward Euler convention in \ref Step intentionally keeps
!> `alpha_step = 1`. Douglas-Gunn, Peaceman-Rachford, Backward Euler, and
!> related ADI/iGRM schemes should enter through \ref MultiStep via the
!> wrapper layer in \ref time_scheme.
!
!------------------------------------------------------------------------------
module ADSS

      use ads_lifecycle, ONLY: initialize, initialize_setup, ComputeDecomposition, &
                              AllocateADSdata, AllocateADS, Cleanup_data, Cleanup_ADS
      use ads_directional_solve, ONLY: solve_problem
      use mumps_solver, ONLY: SolveOneDirection
      use reorderRHS, ONLY: NormalizeTrialBufferToXYZ
      use solution_output, ONLY: PrintSolution

      implicit none

contains







!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Performs one full simulation step composed of three directional
!> substeps.
!>
!> @details
!> This routine executes the three-stage ADS workflow corresponding to
!> the x-, y-, and z-oriented directional solves. For each substep it:
!> - allocates directional RHS/work buffers,
!> - reconstructs the current solution at quadrature points by calling
!>   \ref FormUn,
!> - invokes \ref Sub_Step with the proper directional metadata,
!> - stores intermediate arrays through `move_alloc`,
!> - reinitializes the redistribution buffer \p ads_data%R for the next
!>   stage.
!>
!> The routine is responsible for directional permutations encoded in
!> \p abc and for selecting the appropriate mixing vector from \p mix.
!
! Input:
! ------
!> @param[in] iter
!> Number of the outer iteration or time step.
!>
!> @param[in] mix
!> Mixing coefficients for the three directional substeps.
!>
!> @param[in] RHS_fun
!> Callback used for pointwise RHS evaluation.
!>
!> @param[in] ads_test
!> Test-space setup structure.
!>
!> @param[in] ads_trial
!> Trial-space setup structure.
!>
!> @param[in] n
!> Time-integration or history index passed to lower-level routines.
!>
!> @param[in] alpha_step
!> Coefficient table used by the time-integration scheme.
!
! Input/Output:
! -------------
!> @param[inout] ads_data
!> Runtime data structure updated throughout all substeps.
!
! Output:
! -------
!> @param[out] mierr
!> Returned status code.
!
!---------------------------------------------------------------------------
subroutine MultiStep(iter, mix, RHS_fun, ads_test, ads_trial, ads_data, n, alpha_step, mierr, RHS_point, &
                  lhs_mix, rhs_du_state)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use projection_engine, ONLY: FormUn
      use Interfaces, ONLY: forcing_fun, rhs_point_fun
      implicit none
!> @brief Outer iteration index.
      integer(kind=4), intent(in) :: iter
!> @brief Mixing vectors for the three directional substeps.
      real(kind=8), intent(in) :: mix(4, 3)
!> @brief Time-integration or history index.
      integer(kind=4), intent(in) :: n
!> @brief Coefficient table for the time-integration scheme.
      real(kind=8), intent(in), dimension(7, 3) :: alpha_step
!> @brief Callback used for pointwise RHS evaluation.
      procedure(forcing_fun) :: RHS_fun
!> @brief Optional callback overriding the default pointwise RHS integrand.
      procedure(rhs_point_fun), optional :: RHS_point
!> @brief Optional LHS mixing vectors indexed by physical direction and substep.
      real(kind=8), intent(in), optional :: lhs_mix(4, 3, 3)
!> @brief Optional RHS derivative-state selector indexed by coefficient row and substep.
      integer(kind=4), intent(in), optional :: rhs_du_state(6, 3)
!> @brief Test and trial setup structures.
      type(ADS_setup), intent(in) :: ads_test, ads_trial
!> @brief Runtime data updated in place.
      type(ADS_compute_data), intent(inout) :: ads_data
!> @brief Returned status code.
      integer(kind=4), intent(out) :: mierr
      real(kind=8) :: mmix(4)
      integer(kind=4), dimension(3) :: direction
      integer(kind=4) :: substep
      integer(kind=4), dimension(3, 3) :: abc

      ads_data%rhs_du_state = 0
      if (present(rhs_du_state)) ads_data%rhs_du_state = rhs_du_state

      allocate (ads_data%F (ads_trial%s(1), ads_trial%s(2)*ads_trial%s(3))) !x,y,z
      allocate (ads_data%F2(ads_trial%s(2), ads_trial%s(3)*ads_trial%s(1))) !y,z,x
      allocate (ads_data%F3(ads_trial%s(3), ads_trial%s(1)*ads_trial%s(2))) !z,x,y
      allocate (ads_data%Ft (ads_test%s(1), ads_trial%s(2)*ads_trial%s(3))) !x,y,z
      allocate (ads_data%Ft2(ads_trial%s(2), ads_trial%s(3)*ads_test%s(1))) !y,z,x
      allocate (ads_data%Ft3(ads_trial%s(3), ads_test%s(1)*ads_trial%s(2))) !z,x,y

      mmix = mix(:, 1)
      direction = (/1, 0, 0/) ! x
      abc(:, 1) = (/1, 2, 3/) ! x y z
      abc(:, 2) = (/2, 1, 3/) ! y x z
      abc(:, 3) = (/3, 1, 2/) ! z x y
      substep = 1
      call FormUn(substep, ads_trial, ads_data)
      if (allocated(ads_data%R)) deallocate(ads_data%R)
      allocate (ads_data%R(ads_trial%nrcpp(3)*ads_trial%nrcpp(1)*ads_trial%nrcpp(2), 3, 3, 3))
      ads_data%R = 0.d0
      if (present(lhs_mix)) then
            call Sub_Step(ads_test, ads_trial, iter, mmix, direction, substep, abc, &
                        n, alpha_step, RHS_fun, ads_data, mierr, RHS_point, lhs_mix(:, :, substep))
      else
            call Sub_Step(ads_test, ads_trial, iter, mmix, direction, substep, abc, &
                        n, alpha_step, RHS_fun, ads_data, mierr, RHS_point)
      end if
      if (allocated(ads_data%FFt)) deallocate(ads_data%FFt)
      if (allocated(ads_data%Ft2)) deallocate(ads_data%Ft2)
      if (allocated(ads_data%Ft3)) deallocate(ads_data%Ft3)
      if (allocated(ads_data%FF)) deallocate(ads_data%FF)
      if (allocated(ads_data%F2)) deallocate(ads_data%F2)
      if (allocated(ads_data%F3)) deallocate(ads_data%F3)

      call move_alloc(ads_data%F,ads_data%FF)
      call move_alloc(ads_data%Ft,ads_data%FFt)
      allocate (ads_data%F (ads_trial%s(2), ads_trial%s(1)*ads_trial%s(3))) !y,x,z
      allocate (ads_data%F2(ads_trial%s(3), ads_trial%s(1)*ads_trial%s(2))) !z,x,y
      allocate (ads_data%F3(ads_trial%s(1), ads_trial%s(2)*ads_trial%s(3))) !x,y,z
      allocate (ads_data%Ft (ads_test%s(2), ads_trial%s(1)*ads_trial%s(3))) !y,x,z
      allocate (ads_data%Ft2(ads_trial%s(3), ads_trial%s(1)*ads_test%s(2))) !z.x.y
      allocate (ads_data%Ft3(ads_trial%s(1), ads_test%s(2)*ads_trial%s(3))) !x,y,z

      mmix = mix(:, 2)
      direction = (/0, 1, 0/) ! y
      abc(:, 1) = (/2, 1, 3/) ! y x z
      abc(:, 2) = (/3, 1, 2/) ! z x y
      abc(:, 3) = (/1, 2, 3/) ! x y z
      substep = 2
      call FormUn(substep, ads_trial, ads_data)
      if (allocated(ads_data%R)) deallocate(ads_data%R)
      allocate (ads_data%R(ads_trial%nrcpp(2)*ads_trial%nrcpp(3)*ads_trial%nrcpp(1), 3, 3, 3))
      ads_data%R = 0.d0
      if (present(lhs_mix)) then
            call Sub_Step(ads_test, ads_trial, iter, mmix, direction, substep, abc, &
                        n, alpha_step, RHS_fun, ads_data, mierr, RHS_point, lhs_mix(:, :, substep))
      else
            call Sub_Step(ads_test, ads_trial, iter, mmix, direction, substep, abc, &
                        n, alpha_step, RHS_fun, ads_data, mierr, RHS_point)
      end if
      if (allocated(ads_data%FFt)) deallocate(ads_data%FFt)
      if (allocated(ads_data%Ft2)) deallocate(ads_data%Ft2)
      if (allocated(ads_data%Ft3)) deallocate(ads_data%Ft3)
      if (allocated(ads_data%FF)) deallocate(ads_data%FF)
      if (allocated(ads_data%F2)) deallocate(ads_data%F2)
      if (allocated(ads_data%F3)) deallocate(ads_data%F3)

      call move_alloc(ads_data%F,ads_data%FF)
      call move_alloc(ads_data%Ft,ads_data%FFt)
      allocate (ads_data%F (ads_trial%s(3), ads_trial%s(1)*ads_trial%s(2))) !z,x,y
      allocate (ads_data%F2(ads_trial%s(1), ads_trial%s(2)*ads_trial%s(3))) !x,y,z
      allocate (ads_data%F3(ads_trial%s(2), ads_trial%s(3)*ads_trial%s(1))) !y,z,x
      allocate (ads_data%Ft (ads_test%s(3), ads_trial%s(1)*ads_trial%s(2))) !z,x,y
      allocate (ads_data%Ft2(ads_trial%s(1), ads_trial%s(2)*ads_test%s(3))) !x,y,z
      allocate (ads_data%Ft3(ads_trial%s(2), ads_test%s(3)*ads_trial%s(1))) !y,z,x

      mmix = mix(:, 3)
      direction = (/0, 0, 1/) ! z
      abc(:, 1) = (/3, 1, 2/) ! z x y
      abc(:, 2) = (/1, 2, 3/) ! x y z
      abc(:, 3) = (/2, 1, 3/) ! y x z
      substep = 3
      call FormUn(substep, ads_trial, ads_data)
      if (allocated(ads_data%R)) deallocate(ads_data%R)
      allocate (ads_data%R(ads_trial%nrcpp(1)*ads_trial%nrcpp(2)*ads_trial%nrcpp(3), 3, 3, 3))
      ads_data%R = 0.d0
      if (present(lhs_mix)) then
            call Sub_Step(ads_test, ads_trial, iter, mmix, direction, substep, abc, &
                        n, alpha_step, RHS_fun, ads_data, mierr, RHS_point, lhs_mix(:, :, substep))
      else
            call Sub_Step(ads_test, ads_trial, iter, mmix, direction, substep, abc, &
                        n, alpha_step, RHS_fun, ads_data, mierr, RHS_point)
      end if
      if (allocated(ads_data%FFt)) deallocate(ads_data%FFt)
      if (allocated(ads_data%Ft2)) deallocate(ads_data%Ft2)
      if (allocated(ads_data%Ft3)) deallocate(ads_data%Ft3)
      if (allocated(ads_data%FF)) deallocate(ads_data%FF)
      if (allocated(ads_data%F2)) deallocate(ads_data%F2)
      if (allocated(ads_data%F3)) deallocate(ads_data%F3)

      call move_alloc(ads_data%F,ads_data%FF)
      call move_alloc(ads_data%Ft,ads_data%FFt)
      ads_data%rhs_du_state = 0

end subroutine MultiStep

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Performs a Forward Euler single-step solve without directional
!> enrichment.
!>
!> @details
!> This routine executes one single-substep workflow in which:
!> - the mixing vector is fixed to a pure mass contribution,
!> - `alpha_step` is intentionally fixed to one for the Forward Euler RHS,
!> - no enriched direction is activated,
!> - the current solution is reconstructed through \ref FormUn,
!> - \ref Sub_Step is called once with identical test and trial spaces.
!>
!> It is intended for a simpler or baseline execution path compared with
!> \ref MultiStep.
!
! Input:
! ------
!> @param[in] iter
!> Outer iteration index.
!>
!> @param[in] RHS_fun
!> Callback used for pointwise RHS evaluation.
!>
!> @param[in] ads
!> Setup structure used simultaneously as test and trial space.
!
! Input/Output:
! -------------
!> @param[inout] ads_data
!> Runtime data structure updated during the solve.
!
! Output:
! -------
!> @param[out] mierr
!> Returned status code.
!
!---------------------------------------------------------------------------
subroutine Step(iter, RHS_fun, ads, ads_data, mierr, RHS_point)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use projection_engine, ONLY: FormUn
      use Interfaces, ONLY: forcing_fun, rhs_point_fun
      implicit none
!> @brief Outer iteration index.
      integer(kind=4), intent(in) :: iter
!> @brief Callback used for pointwise RHS evaluation.
      procedure(forcing_fun) :: RHS_fun
!> @brief Optional callback overriding the default pointwise RHS integrand.
      procedure(rhs_point_fun), optional :: RHS_point
!> @brief Setup structure used as both test and trial space.
      type(ADS_setup), intent(in) :: ads
!> @brief Runtime data updated during the solve.
      type(ADS_compute_data), intent(inout) :: ads_data
!> @brief Returned status code.
      integer(kind=4), intent(out) :: mierr
      real(kind=8) :: mix(4)
      integer(kind=4), dimension(3) :: direction
      integer(kind=4) :: substep
      integer(kind=4), dimension(3, 3) :: abc
      real(kind=8), dimension(7, 3) :: alpha_step

      mix = (/1.d0, 0.d0, 0.d0, 0.d0/)
      direction = (/0, 0, 0/)
      abc(:, 1) = (/1, 2, 3/) ! x y z
      abc(:, 2) = (/2, 1, 3/) ! y x z
      abc(:, 3) = (/3, 1, 2/) ! z x y
      substep = 1
      ads_data%un13 = 0.d0
      ads_data%un23 = 0.d0
      call FormUn(1, ads, ads_data)
      alpha_step = 1.d0


      if (allocated(ads_data%R)) deallocate(ads_data%R)
      allocate (ads_data%R(ads%nrcpp(3)*ads%nrcpp(1)*ads%nrcpp(2), 3, 3, 3))
      ads_data%R = 0.d0

      allocate (ads_data%F (ads%s(1), ads%s(2)*ads%s(3))) !x,y,z
      allocate (ads_data%F2(ads%s(2), ads%s(3)*ads%s(1))) !y,x,z
      allocate (ads_data%F3(ads%s(3), ads%s(1)*ads%s(2))) !z,x,y

      !call Sub_Step(ads, ads, iter, mix,direction,substep,abc,RHS_fun,ads_data, mierr)
      call Sub_Step(ads, ads, iter, mix, direction, substep, abc, &
                        1, alpha_step, RHS_fun, ads_data, mierr, RHS_point)
      if (allocated(ads_data%FF)) deallocate(ads_data%FF)
      if (allocated(ads_data%F2)) deallocate(ads_data%F2)
      if (allocated(ads_data%F3)) deallocate(ads_data%F3)

      call move_alloc(ads_data%F,ads_data%FF)

end subroutine Step

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Executes one ADS substep, including RHS formation, directional
!> solves, and redistribution of the resulting spline coefficients.
!>
!> @details
!> This routine is the central orchestration layer for a single ADS
!> substep. It:
!> - resets current RHS buffers,
!> - forms the three-dimensional right-hand side by calling
!>   \ref Form3DRHS,
!> - solves three one-dimensional directional problems through
!>   \ref solve_problem using the permutation table \p abc,
!> - injects the resulting coefficients into the central block of
!>   \p ads_data%R,
!> - redistributes the updated coefficient block with `DistributeSpline`.
!>
!> The routine operates both in the standard and iGRM-enriched modes.
!
! Input:
! ------
!> @param[in] ads_test
!> Test-space setup structure.
!>
!> @param[in] ads_trial
!> Trial-space setup structure.
!>
!> @param[in] iter
!> Outer iteration index.
!>
!> @param[in] mix
!> Mixing vector for the operator matrices.
!>
!> @param[in] direction
!> Directional enrichment selector.
!>
!> @param[in] substep
!> Current substep number.
!>
!> @param[in] abc
!> Permutation table controlling the sequence of directional solves.
!>
!> @param[in] n
!> Time-integration or history index.
!>
!> @param[in] alpha_step
!> Coefficient table of the time-integration scheme.
!>
!> @param[in] RHS_fun
!> Callback used for pointwise RHS evaluation.
!
! Input/Output:
! -------------
!> @param[inout] ads_data
!> Runtime data structure updated throughout the substep.
!
! Output:
! -------
!> @param[out] mierr
!> Returned status code.
!
!---------------------------------------------------------------------------
subroutine Sub_Step(ads_test, ads_trial, iter, mix, direction, substep, abc, &
                  n, alpha_step, &
                  RHS_fun, ads_data, mierr, RHS_point, lhs_mix)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      ! use parallelism, ONLY: PRINTRANK, MYRANKX, MYRANKY, MYRANKZ
      ! use communicators, ONLY: COMMX, COMMY, COMMZ
      use reorderRHS, ONLY: ReorderRHSForX, ReorderRHSForY, ReorderRHSForZ
      use projection_engine, ONLY: Form3DRHS, ComputeMatrix
      use my_mpi, ONLY: DistributeSpline, Gather, Scatter
      use Interfaces, ONLY: forcing_fun, rhs_point_fun
      use mpi
      use sparse
      implicit none
!> @brief Test-space setup structure.
      type(ADS_setup), intent(in) :: ads_test
!> @brief Trial-space setup structure.
      type(ADS_setup), intent(in) :: ads_trial
!> @brief Outer iteration index.
      integer(kind=4), intent(in) :: iter
!> @brief Mixing vector for operator assembly.
      real(kind=8), intent(in) :: mix(4)
!> @brief Directional enrichment selector.
      integer(kind=4), intent(in), dimension(3) :: direction
!> @brief Current substep number.
      integer(kind=4), intent(in) :: substep
!> @brief Permutation table defining directional solve order.
      integer(kind=4), dimension(3, 3), intent(in) :: abc
!> @brief Time-integration or history index.
      integer(kind=4), intent(in) :: n
!> @brief Coefficient table of the time-integration scheme.
      real(kind=8), intent(in), dimension(7, 3) :: alpha_step
!> @brief Callback used for pointwise RHS evaluation.
      procedure(forcing_fun) :: RHS_fun
!> @brief Optional callback overriding the default pointwise RHS integrand.
      procedure(rhs_point_fun), optional :: RHS_point
!> @brief Optional LHS mixing vectors indexed by physical direction.
      real(kind=8), intent(in), optional :: lhs_mix(4, 3)
!> @brief Runtime data updated in place.
      type(ADS_compute_data), intent(inout) :: ads_data
!> @brief Returned status code.
      integer(kind=4), intent(out) :: mierr
!> @brief Minimum number of copied values for OpenMP redistribution-buffer fill.
      integer(kind=4), parameter :: COPY_TO_R_OMP_THRESHOLD = 262144
      integer(kind=4) :: i,a,b,c
      integer(kind=4) :: copy_columns, copy_size, copy_width
      integer(kind=4) :: diridx
      integer(kind=4) :: ierr!, iret
      integer(kind=4), dimension(3) :: nrcpp
      real(kind=8) :: mass_mix(4)
      real(kind=8) :: solve_mix(4, 3)
      real(kind=8) :: solve_mixA(4), solve_mixB(4), solve_mixBT(4)
      ! real(kind=8) :: time1, time2
      logical :: igrm

      mass_mix = (/1.d0, 0.d0, 0.d0, 0.d0/)
      if (present(lhs_mix)) then
            solve_mix = lhs_mix
      else
            do diridx = 1, 3
                  solve_mix(:, diridx) = mix
            end do
      end if

#ifdef PERFORMANCE
      time1 = MPI_Wtime()
#endif

      if (allocated(ads_data%F)) ads_data%F = 0.d0
      if (allocated(ads_data%Ft)) ads_data%Ft = 0.d0
      ! generate the RHS vectors
      call Form3DRHS(ads_test, ads_trial, ads_data, direction, n, substep, &
                     alpha_step, RHS_fun, igrm, RHS_point)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write (*, *) "Form 3D RHS: ", time2 - time1
#endif

#ifdef IPRINT
      write (*, *) PRINTRANK, 'F'
      do i = 1, ads_trial%s(1)
            write (*, *) PRINTRANK, ads_trial%F(i, 1:ads_trial%s(2)*ads_trial%s(3))
      end do
#endif

      a=abc(1, 1)
      b=abc(2, 1)
      c=abc(3, 1)

!--------------------------------------------------------------------
! Solve the first problem
!--------------------------------------------------------------------
      if (direction(a) .EQ. 1) then
            solve_mixA = mass_mix
      else
            solve_mixA = solve_mix(:, a)
      end if
      solve_mixB = solve_mix(:, a)
      solve_mixBT = solve_mix(:, a)
      call solve_problem(ads_test, ads_trial, abc(1, 1), abc(2, 1), abc(3, 1), &
                        solve_mixA, solve_mixB, solve_mixBT, direction, igrm, &
                        ads_data%F, ads_data%F2, ads_data%Ft, ads_data%Ft2, ierr)

!--------------------------------------------------------------------
! Solve the second problem
!--------------------------------------------------------------------
      a=abc(1, 2)
      if (direction(a) .EQ. 1) then
            solve_mixA = mass_mix
      else
            solve_mixA = solve_mix(:, a)
      end if
      solve_mixB = solve_mix(:, a)
      solve_mixBT = solve_mix(:, a)
      call solve_problem(ads_test, ads_trial, abc(1, 2), abc(2, 2), abc(3, 2), &
                        solve_mixA, solve_mixB, solve_mixBT, direction, igrm, &
                        ads_data%F2, ads_data%F3, ads_data%Ft2, ads_data%Ft3, ierr)

!--------------------------------------------------------------------
! Solve the third problem
!--------------------------------------------------------------------
      a=abc(1, 3)
      if (direction(a) .EQ. 1) then
            solve_mixA = mass_mix
      else
            solve_mixA = solve_mix(:, a)
      end if
      solve_mixB = solve_mix(:, a)
      solve_mixBT = solve_mix(:, a)
      call solve_problem(ads_test, ads_trial, abc(1, 3), abc(2, 3), abc(3, 3), &
                        solve_mixA, solve_mixB, solve_mixBT, direction, igrm, &
                        ads_data%F3, ads_data%F, ads_data%Ft3, ads_data%Ft, ierr)

#ifdef IINFO
      write (*, *) PRINTRANK, '3e) DISTRIBUTE SOLUTION'
#endif
      call NormalizeTrialBufferToXYZ(ads_trial, abc(:, 1), ads_data%F)
      a = 1
      b = 2
      c = 3

      !  copy results to proper buffer
      copy_width = ads_trial%s(a)
      copy_columns = ads_trial%s(b)*ads_trial%s(c)
      copy_size = copy_width*copy_columns
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC) &
!$OMP IF(copy_size > COPY_TO_R_OMP_THRESHOLD)
      do i = 1, copy_columns
            ads_data%R((i - 1)*copy_width + 1:i*copy_width, 2, 2, 2) = ads_data%F(:, i)
      end do
!$OMP END PARALLEL DO
      !  nrcpp - number of columns (average) per processor
      nrcpp = (/ads_trial%nrcpp(c), ads_trial%nrcpp(a), ads_trial%nrcpp(b)/)
      call DistributeSpline(ads_data%R, nrcpp)

#ifdef IPRINT
      write (*, *) PRINTRANK, 'Result:'
      do i = 1, ads_trial%s(3)
            write (*, *) PRINTRANK, i, 'row=', ads_trial%F(i, :)
      end do
#endif

      call mpi_barrier(MPI_COMM_WORLD, ierr)

      mierr = 0
end subroutine Sub_Step





end module ADSS
