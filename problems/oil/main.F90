

program main

   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK
   use parallelism, ONLY: PRINTRANK, InitializeParallelism, CleanParallelism
   use communicators, ONLY: CreateCommunicators
   use RHS_eq
   use ADSS
   use input_data
   use mpi

   implicit none


   ! Iteration counter
   integer :: iter = 0

   integer(kind = 4) :: ierr
   integer(kind = 4), dimension(3) :: p1, p2

   type (ADS_setup) :: ads
   type (ADS_compute_data) :: ads_data


   integer :: values(1:8), k
   integer, dimension(:), allocatable :: seed
   real(8) :: r

   call date_and_time(values = values)
   values = (/ 0.d8, 0.d7, 0.d6, 0.d5, 0.d4, 0.d3, 0.d2, 0.d1 /)
   call random_seed(size = k)
   allocate(seed(1:k))
   seed(:) = values(3)
   call random_seed(put = seed)


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
   allocate(Kqvals(ads % p(1) + 1, ads % p(2) + 1, ads % p(3) + 1, ads % maxe(1) - ads % mine(1) + 1, &
   ads % maxe(2) - ads % mine(2) + 1, ads % maxe(3) - ads % mine(3) + 1))
   call InitInputData
   call PrecomputeKq(ads)

   ! Iterations
   do iter = 0, steps

      l2norm = 0
      call Step(iter, ComputePointForRHS, ads, ads_data, ierr)
      call MPI_Reduce(l2norm, fullnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (MYRANK == 0) then
         write(*, *) iter, 'L2 norm:', fullnorm
      endif
      t = t + Dt

   enddo

   call ComputeResults
   call Cleanup(ads, ads_data, ierr)
   call CleanParallelism(ierr)

end
