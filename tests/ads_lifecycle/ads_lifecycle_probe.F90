program ads_lifecycle_probe
   use Setup, only: ADS_Setup
   use ads_lifecycle, only: initialize_setup
   use parallelism, only: InitializeParallelism
   implicit none

   type(ADS_Setup) :: ads
   integer(kind=4), parameter :: n(3) = (/0, 0, 0/)
   integer(kind=4), parameter :: p(3) = (/0, 0, 0/)
   integer(kind=4), parameter :: continuity(3) = (/0, 0, 0/)
   integer(kind=4), parameter :: ng(3) = (/1, 1, 1/)
   integer(kind=4) :: ierr

   call InitializeParallelism(2, 1, 1, ierr)
   call initialize_setup(n, p, continuity, ng, ads, ierr)

   write (*, '(A)') 'UNEXPECTED SUCCESS'
   stop 98
end program ads_lifecycle_probe
