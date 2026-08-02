!------------------------------------------------------------------------------
!
! PROGRAM: main
!
! DESCRIPTION:
!> @file main.F90
!> @brief Driver for the iGRM heat transient example.
!>
!> @details
!> The program initializes separate test and trial ADS spaces, projects
!> the heat initial condition with a mass-only iGRM setup, advances the
!> solution with the selected iGRM time scheme, and exports the
!> trial-space solution after each step.
!
!------------------------------------------------------------------------------
program main

   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK
   use parallelism, ONLY: PRINTRANK, InitializeParallelism, &
      Cleanup_Parallelism, AbortOnError
   use communicators, ONLY: CreateCommunicators, Cleanup_Communicators
   use time_scheme, ONLY: BackwardEuler3DStep, ConfigureBackwardEuler3DTimeScheme, &
                          ConfigureDouglasGunn3DTimeScheme, ConfigureMassOnly3DTimeScheme, &
                          ConfigurePeacemanRachford3DTimeScheme, DouglasGunn3DStep, &
                          PeacemanRachford3DStep, TimeScheme3D, ValidateSpaces
   use RHS_fun
   use ADSS
   use input_data
   use mpi
   use plot

   implicit none

!> @brief Iteration counter.
   integer(kind = 4) :: iter = 0
!> @brief MPI return code.
   integer(kind = 4) :: ierr
!> @brief Time-history index used by the iGRM step wrappers.
   integer(kind = 4) :: nn
!> @brief Test and trial setup structures.
   type(ADS_setup) :: ads_test, ads_trial
!> @brief Runtime data buffers.
   type(ADS_compute_data) :: ads_data
!> @brief Initial projection and physical time-step schemes.
   type(TimeScheme3D) :: init_scheme, heat_scheme

#ifdef DEBUG
   write (*, *) 'debug'
#endif

   t = 0.d0

   call InitializeParameters

   call InitializeParallelism(procx, procy, procz, ierr)
   call CreateCommunicators(ierr)
   call Initialize(nelem, p_test, p_trial, p_trial - 1, ads_test, ads_trial, ads_data, ierr)
   call ValidateSpaces(ads_test, ads_trial)

   call ConfigureMassOnly3DTimeScheme(init_scheme)
   select case (selected_time_scheme)
   case (TIME_SCHEME_DG)
      call ConfigureDouglasGunn3DTimeScheme(Dt, heat_scheme, include_transport=.false.)
   case (TIME_SCHEME_PR)
      call ConfigurePeacemanRachford3DTimeScheme(Dt, heat_scheme, include_transport=.false.)
   case (TIME_SCHEME_BE)
      call ConfigureBackwardEuler3DTimeScheme(Dt, heat_scheme, include_transport=.false.)
   end select

   nn = 1
   do iter = 0, steps
      if (t > 0.d0) then
         ads_test%tau = Dt
         ads_trial%tau = Dt
         call RunHeatStep(heat_scheme, iter, ads_test, ads_trial, ads_data, nn, ierr)
      else
         ads_test%tau = 1.d0
         ads_trial%tau = 1.d0
         call RunHeatStep(init_scheme, iter, ads_test, ads_trial, ads_data, nn, ierr)
      end if
      call AbortOnError(ierr, 'iGRM heat step')
      if (MYRANK == 0) then
         write(*, *) iter
      endif
      t = t + Dt
      call PrintSolution(iter, ads_trial, ads_data%FF)
   end do

   call Cleanup_ADS(ads_test, ierr)
   call Cleanup_ADS(ads_trial, ierr)
   call Cleanup_data(ads_data, ierr)
   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Dispatches one iGRM heat step through the selected scheme wrapper.
!>
!> @details
!> The same dispatcher is used for the mass-only initialization scheme and
!> the physical heat time-step scheme. The wrapper chosen by the
!> command-line scheme name is reused so that data movement follows the
!> same path for initialization and time stepping.
!
!---------------------------------------------------------------------------
   subroutine RunHeatStep(scheme, iter, ads_test, ads_trial, ads_data, nn, ierr)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use time_scheme, ONLY: BackwardEuler3DStep, DouglasGunn3DStep, PeacemanRachford3DStep, TimeScheme3D
      implicit none
!> @brief Persistent iGRM scheme configuration.
      type(TimeScheme3D), intent(in) :: scheme
!> @brief Iteration index.
      integer(kind = 4), intent(in) :: iter
!> @brief Time-history index used by the iGRM wrappers.
      integer(kind = 4), intent(in) :: nn
!> @brief Test and trial spaces.
      type(ADS_Setup), intent(in) :: ads_test, ads_trial
!> @brief Runtime data updated by the selected wrapper.
      type(ADS_compute_data), intent(inout) :: ads_data
!> @brief MPI/status return code.
      integer(kind = 4), intent(out) :: ierr

      select case (selected_time_scheme)
      case (TIME_SCHEME_DG)
         call DouglasGunn3DStep(scheme, iter, forcing, ads_test, ads_trial, ads_data, nn, ierr, &
                                heat_igrm_rhs_point)
      case (TIME_SCHEME_PR)
         call PeacemanRachford3DStep(scheme, iter, forcing, ads_test, ads_trial, ads_data, nn, ierr, &
                                     heat_igrm_rhs_point)
      case (TIME_SCHEME_BE)
         call BackwardEuler3DStep(scheme, iter, forcing, ads_test, ads_trial, ads_data, nn, ierr, &
                                  heat_igrm_rhs_point)
      end select

   end subroutine RunHeatStep

end program main
