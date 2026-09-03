!------------------------------------------------------------------------------
!> @file main.F90
!> @brief Driver for the corrected three-dimensional DPG pollution problem.
!------------------------------------------------------------------------------
program main

   use input_data, only: InitializeParameters, nelem, adapt_mesh, p_trial, &
                         c_trial, p_test, c_test, steps, procx, procy, procz, &
                         output_resolution, TIME_STEP, DOMAIN_LENGTH, &
                         DIFFUSION, SOURCE_CENTER, SOURCE_RADIUS, wind_at_time
   use parallelism, only: MYRANK, InitializeParallelism, Cleanup_Parallelism, &
                          AbortOnError
   use pollution_dpg_solver, only: PollutionDPGSolution, PollutionDPGStats, &
                                   SolvePollutionDPG3D, &
                                   CleanupPollutionDPGSolution
   use RHS_fun, only: forcing

   implicit none

   type(PollutionDPGSolution) :: solution
   type(PollutionDPGStats) :: stats
   real(kind=8) :: source_min(3), source_max(3)
   integer(kind=4) :: ierr, cleanup_ierr

   call InitializeParameters
   call InitializeParallelism(procx, procy, procz, ierr)
   call AbortOnError(ierr, 'parallel initialization')

   source_min = SOURCE_CENTER - SOURCE_RADIUS
   source_max = SOURCE_CENTER + SOURCE_RADIUS
   call SolvePollutionDPG3D(nelem, adapt_mesh, p_trial, c_trial, p_test, &
                            c_test, steps, TIME_STEP, DOMAIN_LENGTH, &
                            DIFFUSION, wind_at_time, forcing, source_min, &
                            source_max, output_resolution, .true., solution, &
                            stats, ierr)
   call AbortOnError(ierr, 'iGRM pollution solve')

   if (MYRANK == 0) then
      write (*, '(A,I0)') 'iGRM pollution completed steps: ', &
                           stats%completed_steps
      write (*, '(A,ES24.16)') 'iGRM pollution time: ', stats%time
      write (*, '(A,ES24.16)') 'iGRM pollution source integral: ', &
                               stats%source_integral
      write (*, '(A,ES24.16)') 'iGRM pollution L2 norm: ', stats%l2_norm
      write (*, '(A,ES24.16)') 'iGRM pollution total mass: ', &
                               stats%total_mass
      write (*, '(A,ES24.16)') 'iGRM pollution maximum coefficient abs: ', &
                               stats%maximum_abs
      write (*, '(A,ES24.16)') 'iGRM pollution relative residual: ', &
                               stats%maximum_relative_residual
      write (*, '(A,3(ES24.16,1X))') 'iGRM pollution wind: ', stats%wind
   end if

   call CleanupPollutionDPGSolution(solution)
   call Cleanup_Parallelism(cleanup_ierr)

end program main
