

program main

   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK
   use parallelism, ONLY: PRINTRANK, InitializeParallelism, Cleanup_Parallelism
   use communicators, ONLY: CreateCommunicators
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
   
   logical :: prnt = .FALSE.
   logical :: ok = .TRUE.

   type (ADS_setup) :: ads_test, ads_trial
   type (ADS_compute_data) :: ads_data
   
   real (kind = 8) :: epsilon = 1.E-10

   real (kind = 8) :: l2norm, fullnorm

   real(kind=8) :: mix(4,3)
   real (kind=8), dimension(7,3) :: alpha_step
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

   fullnorm = 0.d0
   iter = 0
   l2norm = 0.d0
   
   mix(:, 1) = (/1.d0, 0.d0, 0.d0, 0.d0/)
   mix(:, 2) = (/1.d0, 0.d0, 0.d0, 0.d0/)
   mix(:, 3) = (/1.d0, 0.d0, 0.d0, 0.d0/)
   alpha_step=1.d0
   nn=1.d0
   call MultiStep(iter, mix, forcing, ads_test, ads_trial, ads_data,nn,alpha_step, ierr)
   call PrintSolution(iter, ads_trial, ads_data%FF)

   call Cleanup_ADS(ads_test, ierr)
   call Cleanup_ADS(ads_trial, ierr)
   call Cleanup_data(ads_data, ierr)
   call Cleanup_Parallelism(ierr)

end program main
 
