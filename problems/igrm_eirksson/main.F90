!------------------------------------------------------------------------------
!> @file main.F90
!> @brief Stationary three-dimensional Eriksson iGRM-MUMPS driver.
!------------------------------------------------------------------------------
program main

   use Setup, only: ADS_Setup, ADS_compute_data
   use parallelism, only: MYRANK, InitializeParallelism, Cleanup_Parallelism, &
                          AbortOnError
   use communicators, only: CreateCommunicators, Cleanup_Communicators
   use time_scheme, only: ValidateSpaces
   use ADSS, only: Initialize, Cleanup_ADS, Cleanup_data, PrintSolution
   use utils, only: NormL2
   use input_data, only: nelem, p_test, p_trial, procx, procy, procz, &
                         EPSILON_VALUE, BETA, RESIDUAL_TOLERANCE, &
                         InitializeParameters, exact_solution
   use RHS_fun, only: forcing
   use igrm_mumps_solver, only: IGRMMumpsStats, SolveIGRMMumps3D

   implicit none

   type(ADS_Setup) :: ads_test, ads_trial
   type(ADS_compute_data) :: ads_data
   type(IGRMMumpsStats) :: stats
   real(kind=8), allocatable :: solution(:, :)
   real(kind=8) :: l2_error
   integer(kind=4) :: ierr

   call InitializeParameters

   call InitializeParallelism(procx, procy, procz, ierr)
   call AbortOnError(ierr, 'parallel initialization')
   call CreateCommunicators(ierr)
   call AbortOnError(ierr, 'communicator initialization')

   call Initialize(nelem, p_test, p_trial, p_trial - 1, ads_test, ads_trial, &
                   ads_data, ierr)
   call AbortOnError(ierr, 'iGRM space initialization')
   call ValidateSpaces(ads_test, ads_trial)

   call SolveIGRMMumps3D(ads_test, ads_trial, EPSILON_VALUE, BETA, &
                         RESIDUAL_TOLERANCE, forcing, solution, stats, ierr)
   call AbortOnError(ierr, 'iGRM-MUMPS direct solve')

   call NormL2(ads_trial, solution, l2_error, exact_solution)

   if (MYRANK == 0) then
      write (*, '(A,1X,I0)') 'iGRM-MUMPS system size:', stats%system_size
      write (*, '(A,1X,I0)') 'iGRM-MUMPS nonzeros:', stats%nonzeros
      write (*, '(A,1X,ES24.16)') 'iGRM-MUMPS residual RMS:', &
                                  stats%residual_rms
      write (*, '(A,1X,ES24.16)') 'iGRM-MUMPS relative residual:', &
                                  stats%relative_residual
      write (*, '(A,1X,ES24.16)') 'L2 error:', l2_error
   end if

   call PrintSolution(0, ads_trial, solution)

   if (allocated(solution)) deallocate (solution)
   call Cleanup_ADS(ads_test, ierr)
   call Cleanup_ADS(ads_trial, ierr)
   call Cleanup_data(ads_data, ierr)
   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

end program main
