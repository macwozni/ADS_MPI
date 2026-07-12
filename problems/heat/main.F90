!------------------------------------------------------------------------------
!
! PROGRAM: main
!
! DESCRIPTION:
!> @file main.F90
!> @brief Driver for the heat transient example.
!>
!> @details
!> The program initializes the ADS runtime, advances the solution for the
!> requested number of steps using the forcing from \ref input_data, prints
!> iteration numbers on rank zero, and exports the solution after each step.
!
!------------------------------------------------------------------------------
program main

   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK
   use parallelism, ONLY: PRINTRANK, InitializeParallelism, Cleanup_Parallelism
   use communicators, ONLY: CreateCommunicators, Cleanup_Communicators
   use time_scheme, ONLY: ForwardEuler3DStep
   use ADSS
   use RHS_fun
   use input_data
   use mpi
   use plot

   implicit none

   ! Iteration counter
   integer :: iter = 0

   integer(kind = 4) :: ierr
   integer(kind = 4), dimension(3) :: nelem, p

   type (ADS_setup) :: ads_test, ads_trial
   type (ADS_compute_data) :: ads_data

#ifdef DEBUG
   write (*, *) 'debug'
#endif

   t = 0

   ! -------------------------------------------------------------------
   ! Code
   ! -------------------------------------------------------------------

   call InitializeParameters

   ! prepare the problem dimensions

   call InitializeParallelism(procx, procy, procz, ierr)
   call CreateCommunicators(ierr)
   nelem = (/ SIZE, SIZE, SIZE /)
   p = (/ ORDER, ORDER, ORDER /)
   call Initialize(nelem, p, p, p - 1, ads_test, ads_trial, ads_data, ierr)
   
   ! Iterations
   do iter = 0, steps

      if (t > 0.d0) then
         ads_trial%tau = Dt
      else
         ads_trial%tau = 1.d0
      endif
      call ForwardEuler3DStep(iter, forcing, ads_trial, ads_data, ierr, heat_rhs_point)
      if (MYRANK == 0) then
         write(*, *) iter
      endif
      t = t + Dt
      call PrintSolution(iter, ads_trial, ads_data%FF)

   enddo

   call Cleanup_ADS(ads_test, ierr)
   call Cleanup_ADS(ads_trial, ierr)
   call Cleanup_data(ads_data, ierr)
   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

end program main
