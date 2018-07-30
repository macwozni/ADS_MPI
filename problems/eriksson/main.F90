

program main

   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK
   use parallelism, ONLY: PRINTRANK, InitializeParallelism, CleanParallelism
   use communicators, ONLY: CreateCommunicators
   use RHS_eq
   use ADSS
   use input_data
   use mpi
   use plot

   implicit none

   ! Iteration counter
   integer :: iter = 0

   integer(kind = 4) :: ierr
   integer(kind = 4), dimension(3) :: p1, p2

   type (ADS_setup) :: ads
   type (ADS_compute_data) :: ads_data
   
   real (kind = 8) :: l2norm, fullnorm

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
   p1 = (/ SIZE, SIZE, SIZE /)
   p2 = (/ ORDER, ORDER, ORDER /)
   call Initialize(p1, p2, ads, ads_data, ierr)
   
   ! Iterations
   do iter = 0, steps

      l2norm = 0.d0
      call Step(iter, ComputePointForRHS, ads, ads_data, l2norm, ierr)
      call MPI_Reduce(l2norm, fullnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (MYRANK == 0) then
         write(*, *) iter, 'L2 norm:', fullnorm
      endif
      t = t + Dt
      call PrintSolution(iter, t, ads, ads_data)

   enddo

   call Cleanup(ads, ads_data, ierr)
   call CleanParallelism(ierr)

end

