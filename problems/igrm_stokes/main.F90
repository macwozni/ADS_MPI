!------------------------------------------------------------------------------
!> @file main.F90
!> @brief Stationary three-dimensional DGiGRM Stokes driver.
!------------------------------------------------------------------------------
program main

   use parallelism, only: MYRANK, InitializeParallelism, Cleanup_Parallelism, &
                          AbortOnError
   use communicators, only: CreateCommunicators, Cleanup_Communicators
   use input_data, only: nelem, p_test, p_trial, procx, procy, procz, &
                         VISCOSITY_VALUE, PENALTY_FACTOR, &
                         RESIDUAL_TOLERANCE, InitializeParameters, &
                         exact_velocity, exact_pressure
   use RHS_fun, only: forcing_x, forcing_y, forcing_z
   use igrm_stokes_solver, only: IGRMStokesSolution, IGRMStokesStats, &
                                  SolveIGRMStokes3D, &
                                  ComputeIGRMStokesErrors, &
                                  WriteIGRMStokesVTI

   implicit none

   type(IGRMStokesSolution) :: solution
   type(IGRMStokesStats) :: stats
   real(kind=8) :: velocity_l2, pressure_l2, divergence_l2
   integer(kind=4) :: ierr

   call InitializeParameters

   call InitializeParallelism(procx, procy, procz, ierr)
   call AbortOnError(ierr, 'parallel initialization')
   call CreateCommunicators(ierr)
   call AbortOnError(ierr, 'communicator initialization')

   call SolveIGRMStokes3D(nelem, p_test, p_trial, VISCOSITY_VALUE, &
                          PENALTY_FACTOR, &
                          RESIDUAL_TOLERANCE, forcing_x, forcing_y, forcing_z, &
                          exact_velocity, solution, stats, ierr)
   call AbortOnError(ierr, 'iGRM-Stokes direct solve')

   call ComputeIGRMStokesErrors(solution, exact_velocity, exact_pressure, &
                                velocity_l2, pressure_l2, divergence_l2)

   if (MYRANK == 0) then
      write (*, '(A,1X,I0)') 'iGRM-Stokes system size:', stats%system_size
      write (*, '(A,1X,I0)') 'iGRM-Stokes nonzeros:', stats%nonzeros
      write (*, '(A,1X,ES24.16)') 'iGRM-Stokes residual RMS:', &
                                  stats%residual_rms
      write (*, '(A,1X,ES24.16)') 'iGRM-Stokes relative residual:', &
                                  stats%relative_residual
      write (*, '(A,1X,ES24.16)') 'L2 velocity error:', velocity_l2
      write (*, '(A,1X,ES24.16)') 'L2 pressure error:', pressure_l2
      write (*, '(A,1X,ES24.16)') 'L2 divergence:', divergence_l2
   end if

   call WriteIGRMStokesVTI(solution, 'result.vti', ierr)
   call AbortOnError(ierr, 'iGRM-Stokes VTI output')

   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

end program main
