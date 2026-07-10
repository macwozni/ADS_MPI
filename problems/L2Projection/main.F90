!------------------------------------------------------------------------------
!
! PROGRAM: main
!
! DESCRIPTION:
!> @file main.F90
!> @brief Driver for the scalar L2 projection regression problem.
!>
!> @details
!> The program initializes the ADS runtime, projects the constant forcing
!> supplied by \ref input_data, and checks on rank zero that the resulting
!> coefficient field remains equal to one within a small tolerance.
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
   integer(kind = 4) :: iter = 0
   integer(kind = 4) :: i, j

   integer(kind = 4) :: ierr
   integer(kind = 4), dimension(3) :: nelem, p
   
   logical :: prnt = .FALSE.
   logical :: ok = .TRUE.

   type (ADS_setup) :: ads_test, ads_trial
   type (ADS_compute_data) :: ads_data
   
   real (kind = 8) :: epsilon = 1.E-10

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
   nelem = (/ isizex, isizey, isizez /)
   p = (/ order, order, order /)
   call Initialize(nelem, p, p, p - 1, ads_test, ads_trial, ads_data, ierr)
   ads_trial%tau = 1.d0

   iter = 0

   call ForwardEuler3DStep(iter, forcing, ads_trial, ads_data, ierr)

   if (MYRANK == 0) then
      if (prnt) then
         write(*, *) 'Result:'
         do i = 1, ads_trial % s(3)
            write(*, *) i, 'row='
            do j = 1, ads_trial % s(1) * ads_trial % s(2)
               write(*, *) ads_data % FF(i,j)
            enddo
         enddo
      endif
      do i = 1, ads_trial % s(1)
         do j = 1, ads_trial % s(2) * ads_trial % s(3)
            if (abs(ads_data % FF(i,j) - 1.d0) > epsilon) then
               ok = .FALSE.
            endif
         enddo
      enddo
      if (ok .eqv. .FALSE.) then
         write(*,*) 'not OK'
      endif
   endif

   call Cleanup_ADS(ads_test, ierr)
   call Cleanup_ADS(ads_trial, ierr)
   call Cleanup_data(ads_data, ierr)
   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

end program main
 
