

program main

   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK
   use parallelism, ONLY: PRINTRANK, InitializeParallelism, CleanParallelism
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
   integer(kind = 4), dimension(3) :: n, p1, p2
   
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
   n = (/ isizex, isizey, isizez /)
   p1 = (/ order, order, order /)
   p2 = (/ order, order, order /)
   call Initialize(n, p1, p2, p1-1, ads_test, ads_trial, ads_data, ierr)

   fullnorm = 0.d0
   iter = 0
   l2norm = 0.d0
   
   mix=1.d0
   alpha_step=1.d0
   nn=1.d0
   ! call MultiStep(iter, mix, forcing, ads_test, ads_trial, ads_data,nn,alpha_step, l2norm, ierr)


   call Step(iter, forcing, ads_trial, ads_data, l2norm, ierr)
   call MPI_Reduce(l2norm, fullnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   
   ! if (MYRANK == 0) then
   !    write(*, *)'L2 norm:', fullnorm
   !    write(*,*) (abs(fullnorm - 1.d0) < epsilon)
   !    if (prnt) then
   !       write(*, *) 'Result:'
   !       do i = 1, ads_tial % s(3)
   !          write(*, *) i, 'row='
   !          do j = 1, ads_tial % s(1) * ads_tial % s(2)
   !             write(*, *) ads_data % F(i,j)
   !          enddo
   !       enddo
   !    endif
   !    do i = 1, ads % s(1)
   !       do j = 1, ads_tial % s(2) * ads_tial % s(3)
   !          if (abs(ads_data % F(i,j) - 1.d0) > epsilon) then
   !             ok = .FALSE.
   !          endif
   !       enddo
   !    enddo
   !    if (ok .eqv. .FALSE.) then
   !       write(*,*) 'not OK'
   !    endif
   ! endif

   call Cleanup(ads_test, ads_trial, ads_data, ierr)
   call CleanParallelism(ierr)

end program main
 
