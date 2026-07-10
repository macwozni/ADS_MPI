!------------------------------------------------------------------------------
!
! PROGRAM: main
!
! DESCRIPTION:
!> @file main.F90
!> @brief Driver for the pure-diffusion iGRM example.
!>
!> @details
!> The program initializes the ADS runtime, advances the zero-forcing
!> pure-diffusion example for the requested number of steps, and cleans up
!> the ADS, communicator, and MPI resources.
!
!------------------------------------------------------------------------------
program main

   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK
   use parallelism, ONLY: PRINTRANK, InitializeParallelism, Cleanup_Parallelism
   use communicators, ONLY: CreateCommunicators, Cleanup_Communicators
   use time_scheme, ONLY: BackwardEuler3DStep, DouglasGunn3DStep, PeacemanRachford3DStep, &
                          ValidateIGRMTimeSchemeSpaces
   use ADSS
   use RHS_fun
   use input_data
   use mpi

   implicit none


   ! Iteration counter
   integer :: iter = 0

   integer(kind = 4) :: ierr
   integer(kind = 4), dimension(3) :: nelem, p_test, p_trial
   integer(kind = 4) :: nn

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
   nelem = (/ SIZE, SIZE, SIZE /)
   p_trial = (/ ORDER, ORDER, ORDER /)
   p_test = p_trial + 1
   call Initialize(nelem, p_test, p_trial, p_trial - 1, ads_test, ads_trial, ads_data, ierr)
   call ValidateIGRMTimeSchemeSpaces(ads_test, ads_trial)
   nn = 1
   ! Iterations
   do iter = 0, steps

      if (t > 0.d0) then
         ads_trial%tau = Dt
         ads_test%tau = Dt
      else
         ads_trial%tau = 1.d0
         ads_test%tau = 1.d0
      endif
      select case (trim(time_scheme))
      case ("dg", "douglas-gunn", "douglas_gunn")
         call DouglasGunn3DStep(iter, forcing, ads_test, ads_trial, ads_data, nn, ierr, &
                                include_transport=.false.)
      case ("pr", "peaceman-rachford", "peaceman_rachford")
         call PeacemanRachford3DStep(iter, forcing, ads_test, ads_trial, ads_data, nn, ierr, &
                                     include_transport=.false.)
      case ("be", "backward-euler", "backward_euler", "backwardeuler")
         call BackwardEuler3DStep(iter, forcing, ads_test, ads_trial, ads_data, nn, ierr, &
                                  include_transport=.false.)
      case default
         if (MYRANK == 0) write(*, *) "unknown time scheme: ", trim(time_scheme)
         stop 5
      end select
      if (MYRANK == 0) then
         write(*, *) iter
      endif
      t = t + Dt

   enddo

   call Cleanup_ADS(ads_test, ierr)
   call Cleanup_ADS(ads_trial, ierr)
   call Cleanup_data(ads_data, ierr)
   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

end program main
