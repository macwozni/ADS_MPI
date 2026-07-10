!------------------------------------------------------------------------------
!
! PROGRAM: main
!
! DESCRIPTION:
!> @file main.F90
!> @brief Driver for the iGRM L2 projection example.
!>
!> @details
!> The program initializes separate test and trial ADS spaces, configures
!> a selected iGRM time scheme through the time-scheme wrapper layer, and
!> prints the resulting trial-space solution.
!
!------------------------------------------------------------------------------
program main

   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK
   use parallelism, ONLY: PRINTRANK, InitializeParallelism, Cleanup_Parallelism
   use communicators, ONLY: CreateCommunicators, Cleanup_Communicators
   use time_scheme, ONLY: BackwardEuler3DStep, ConfigureBackwardEuler3DTimeScheme, &
                          ConfigureDouglasGunn3DTimeScheme, ConfigurePeacemanRachford3DTimeScheme, &
                          DouglasGunn3DStep, PeacemanRachford3DStep, TimeScheme3D, &
                          ValidateIGRMTimeSchemeSpaces
   use RHS_fun
   use ADSS
   use input_data
   use mpi

   implicit none


   ! Iteration counter
   integer(kind = 4) :: iter = 0
   integer(kind = 4) :: i, j

   integer(kind = 4) :: ierr
   integer(kind = 4), dimension(3) :: nelem, p1, p2
   integer(kind = 4), parameter :: SCHEME_DG = 1, SCHEME_PR = 2, SCHEME_BE = 3
   integer(kind = 4) :: selected_scheme

   logical :: prnt = .FALSE.
   logical :: ok = .TRUE.

   type (ADS_setup) :: ads_test, ads_trial
   type (ADS_compute_data) :: ads_data
   type (TimeScheme3D) :: scheme

   real (kind = 8) :: epsilon = 1.E-10

   real (kind = 8) :: l2norm, fullnorm

   real (kind=8) :: tau
   integer (kind = 4) :: nn

#ifdef DEBUG
   write (*, *) 'debug'
#endif


   ! -------------------------------------------------------------------
   ! Code
   ! -------------------------------------------------------------------

   call InitializeParameters

   ! prepare the problem dimensions

   call InitializeParallelism(procx, procy, procz, ierr)
   call CreateCommunicators(ierr)
   !nelem = (/ isizex, isizey, isizez /)
   !p1 = (/ order, order, order /)
   !p2 = (/ order, order, order /)
   nelem = (/ 2,2,2 /)
   p1 = (/3,3,3/)
   p2 = (/1,1,1/)
   call Initialize(nelem, p1, p2, p2-1, ads_test, ads_trial, ads_data, ierr)
   call ValidateIGRMTimeSchemeSpaces(ads_test, ads_trial)
   tau = scheme_tau
   ads_test%tau = tau
   ads_trial%tau = tau
   select case (trim(time_scheme))
   case ("dg", "douglas-gunn", "douglas_gunn")
      selected_scheme = SCHEME_DG
      call ConfigureDouglasGunn3DTimeScheme(tau, scheme, include_transport=.false.)
   case ("pr", "peaceman-rachford", "peaceman_rachford")
      selected_scheme = SCHEME_PR
      call ConfigurePeacemanRachford3DTimeScheme(tau, scheme, include_transport=.false.)
   case ("be", "backward-euler", "backward_euler", "backwardeuler")
      selected_scheme = SCHEME_BE
      call ConfigureBackwardEuler3DTimeScheme(tau, scheme, include_transport=.false.)
   case ("fe", "forward-euler", "forward_euler", "forwardeuler")
      if (MYRANK == 0) write(*, *) "forward euler is not an iGRM L2 time-scheme option"
      stop 5
   case default
      if (MYRANK == 0) write(*, *) "unknown time scheme: ", trim(time_scheme)
      stop 5
   end select

   fullnorm = 0.d0
   iter = 0
   l2norm = 0.d0

   nn = 1
   select case (selected_scheme)
   case (SCHEME_DG)
      call DouglasGunn3DStep(scheme, iter, forcing, ads_test, ads_trial, ads_data, nn, ierr)
   case (SCHEME_PR)
      call PeacemanRachford3DStep(scheme, iter, forcing, ads_test, ads_trial, ads_data, nn, ierr)
   case (SCHEME_BE)
      call BackwardEuler3DStep(scheme, iter, forcing, ads_test, ads_trial, ads_data, nn, ierr)
   end select
   call PrintSolution(iter, ads_trial, ads_data%FF)

   call Cleanup_ADS(ads_test, ierr)
   call Cleanup_ADS(ads_trial, ierr)
   call Cleanup_data(ads_data, ierr)
   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

end program main
 
