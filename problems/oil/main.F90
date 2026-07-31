!------------------------------------------------------------------------------
!
! PROGRAM: main
!
! DESCRIPTION:
!> @file main.F90
!> @brief Driver for the oil transport example.
!>
!> @details
!> The program initializes the ADS runtime, prepares the oil-problem input
!> data and permeability cache, advances the solution for the requested
!> number of steps, and prints the accumulated drained quantity.
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

   implicit none


   ! Iteration counter
   integer :: iter = 0

   integer(kind = 4) :: ierr
   integer(kind = 4), dimension(3) :: nelem, p

   type (ADS_setup) :: ads_test, ads_trial
   type (ADS_compute_data) :: ads_data


   integer :: values(1:8), k
   integer, dimension(:), allocatable :: seed
   real(8) :: r

   !call date_and_time(values = values)
   !values = (/ 0.d8, 0.d7, 0.d6, 0.d5, 0.d4, 0.d3, 0.d2, 0.d1 /)
   !call random_seed(size = k)
   !allocate(seed(1:k))
   !seed(:) = values(3)
   !call random_seed(put = seed)


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
   nelem = (/ ELEMENTS_PER_DIRECTION, ELEMENTS_PER_DIRECTION, ELEMENTS_PER_DIRECTION /)
   p = (/ ORDER, ORDER, ORDER /)
   call Initialize(nelem, p, p, p - 1, ads_test, ads_trial, ads_data, ierr)
   allocate(Kqvals(ads_trial % p(1) + 1, ads_trial % p(2) + 1, ads_trial % p(3) + 1, &
   ads_trial % maxe(1) - ads_trial % mine(1) + 1, &
   ads_trial % maxe(2) - ads_trial % mine(2) + 1, ads_trial % maxe(3) - ads_trial % mine(3) + 1))
   call InitializeDrainedAccumulator(ads_trial)
   call InitInputData
   call PrecomputeKq(ads_trial)

   ! Iterations
   do iter = 0, steps

      if (t > 0.d0) then
         ads_trial%tau = Dt
      else
         ads_trial%tau = 1.d0
      endif
      call ForwardEuler3DStep(iter, forcing, ads_trial, ads_data, ierr, oil_rhs_point)
      if (MYRANK == 0) then
         write(*, *) iter
      endif
      t = t + Dt

   enddo

   call ComputeResults
   call CleanupDrainedAccumulator
   call Cleanup_ADS(ads_test, ierr)
   call Cleanup_ADS(ads_trial, ierr)
   call Cleanup_data(ads_data, ierr)
   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

end program main
