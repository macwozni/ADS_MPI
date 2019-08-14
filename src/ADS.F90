module ADSS

contains




! -------------------------------------------------------------------
! Initialization of clocks and MPI
! -------------------------------------------------------------------
subroutine initialize(n, p, ads, ads_data, mierr)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ
   use parallelism, ONLY: PRINTRANK
   use knot_vector, ONLY: PrepareKnot
   use basis, ONLY: BasisData
   use mpi
   implicit none
   integer(kind = 4), intent(in), dimension(3) :: n
   integer(kind = 4), intent(in), dimension(3) :: p
   type(ADS_setup), intent(out) :: ads
   type (ADS_compute_data), intent(out) :: ads_data
   integer(kind = 4), intent(out) :: mierr
   integer(kind = 4) :: ierr
   integer(kind = 4), dimension(3) :: nelem
   real (kind = 8), allocatable, dimension(:) :: Ux
   real (kind = 8), allocatable, dimension(:) :: Uy
   real (kind = 8), allocatable, dimension(:) :: Uz

   call PrepareKnot(n(1), p(1), Ux, nelem(1))
   call PrepareKnot(n(2), p(2), Uy, nelem(2))
   call PrepareKnot(n(3), p(3), Uz, nelem(3))
   
   call AllocateADS(n,nelem,p,ads)
   
   call move_alloc(Ux, ads % Ux)
   call move_alloc(Uy, ads % Uy)
   call move_alloc(Uz, ads % Uz)
   
   ads % p = p ! order
   ads % n = n ! intervals

   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, 'INITIALIZATION'
   write(*, *) 'px', p(1), 'py', p(2), 'pz', p(3), &
   'nx', n(1), 'ny', n(2), 'nz', n(3), &
   'size of Ux', n(1) + p(1) + 2, 'size of Uy', n(2) + p(2) + 2, 'size of Uz', n(3) + p(3) + 2
#endif

   if (n(1) < NRPROCX .or. n(2) < NRPROCY .or. n(3) < NRPROCZ) then
      write(*, *) 'Number of elements smaller than number of processors'
      stop
   endif

   call ComputeDecomposition(ads)

#ifdef IDEBUG
   call ValidateDimensions(&
   ads % n, &
   ads % s, &
   ads % nrcpp, &
   ads % dimensionsX, ads % dimensionsY, ads % dimensionsZ)
#endif
   
#ifdef IPRINT
   call PrintDecompositionInfo(&
   ads % n, &
   ads % nrcpp, &
   ads % ibeg, &
   ads % iend)
#endif

   ads % nelem = nelem
   mierr = 0
   
   ads % m(1) = ads % n(1) + ads % p(1) + 1
   ads % ng(1) = ads % p(1) + 1
   ads % m(2) = ads % n(2) + ads % p(2) + 1
   ads % ng(2) = ads % p(2) + 1
   ads % m(3) = ads % n(3) + ads % p(3) + 1
   ads % ng(3) = ads % p(3) + 1
   
   
   call BasisData(ads % p(1), ads % m(1), ads % Ux, 1, ads % ng(1), &
   ads % nelem(1), ads % Ox, ads % Jx, ads % Wx, ads % Xx, ads % NNx)
   call BasisData(ads % p(2), ads % m(2), ads % Uy, 1, ads % ng(2), &
   ads % nelem(2), ads % Oy, ads % Jy, ads % Wy, ads % Xy, ads % NNy)
   call BasisData(ads % p(3), ads % m(3), ads % Uz, 1, ads % ng(3), &
   ads % nelem(3), ads % Oz, ads % Jz, ads % Wz, ads % Xz, ads % NNz)

   ads % lnelem(1) = ads % maxe(1) - ads % mine(1) + 1
   ads % lnelem(2) = ads % maxe(2) - ads % mine(2) + 1
   ads % lnelem(3) = ads % maxe(3) - ads % mine(3) + 1
   call AllocateADSdata(ads, ads_data)
   
#ifdef IPRINT
      write(*, *) PRINTRANK, 'ex:', ads % mine(1), ads % maxe(1)
      write(*, *) PRINTRANK, 'ey:', ads % mine(2), ads % maxe(2)
      write(*, *) PRINTRANK, 'ez:', ads % mine(3), ads % maxe(3)
      write(*, *) PRINTRANK, 'ibegx,iendx', ads % ibeg(1), ads % iend(1)
      write(*, *) PRINTRANK, 'ibegy,iendy', ads % ibeg(2), ads % iend(2)
      write(*, *) PRINTRANK, 'ibegz,iendz', ads % ibeg(3), ads % iend(3)
#endif
      
end subroutine initialize


! -------------------------------------------------------------------
! Establishes decomposition of the domain. Calculates size and location
! of the piece for current process.
! -------------------------------------------------------------------
subroutine ComputeDecomposition(ads)
   use Setup, ONLY: ADS_Setup
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, PRINTRANK, &
   NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints, FillDimVector
   implicit none
   type(ADS_setup), intent(inout) :: ads
   integer(kind = 4) :: i
   integer(kind = 4) :: ix, iy, iz
   integer(kind = 4) :: imine, imaxe

   ! number of columns per processors
   call ComputeEndpoints(MYRANKX, NRPROCX, ads % n(1), ads % p(1), ads % nrcpp(1), ads % ibeg(1), &
   ads % iend(1), ads % mine(1), ads % maxe(1))
   call ComputeEndpoints(MYRANKY, NRPROCY, ads % n(2), ads % p(2), ads % nrcpp(2), ads % ibeg(2), &
   ads % iend(2), ads % mine(2), ads % maxe(2))
   call ComputeEndpoints(MYRANKZ, NRPROCZ, ads % n(3), ads % p(3), ads % nrcpp(3), ads % ibeg(3), &
   ads % iend(3), ads % mine(3), ads % maxe(3))

   ads % s(1) = ads % iend(1) - ads % ibeg(1) + 1
   ads % s(2) = ads % iend(2) - ads % ibeg(2) + 1
   ads % s(3) = ads % iend(3) - ads % ibeg(3) + 1

#ifdef IINFO
   write(*, *) PRINTRANK, 'Number of cols per processor:', ads % nrcpp(1), ads % nrcpp(2), ads % nrcpp(3)
   write(*, *) PRINTRANK, 'ibegx,iendx', ads % ibeg(1), ads % iend(1)
   write(*, *) PRINTRANK, 'ibegy,iendy', ads % ibeg(2), ads % iend(2)
   write(*, *) PRINTRANK, 'ibegz,iendz', ads % ibeg(3), ads % iend(3)
#endif
   
   ! prepare dimensions vectors
   call FillDimVector(ads % dimensionsX, ads % shiftsX, ads % nrcpp(1), ads % s(2) * ads % s(3), ads % n(1), NRPROCX)
   call FillDimVector(ads % dimensionsY, ads % shiftsY, ads % nrcpp(2), ads % s(1) * ads % s(3), ads % n(2), NRPROCY)
   call FillDimVector(ads % dimensionsZ, ads % shiftsZ, ads % nrcpp(3), ads % s(1) * ads % s(2), ads % n(3), NRPROCZ)

   ! Compute indices for neighbours
   ads % ibegsx = -1
   ads % iendsx = -1
   ads % ibegsy = -1
   ads % iendsy = -1
   ads % ibegsz = -1
   ads % iendsz = -1

   do i = max(MYRANKX - 1, 0) + 1, min(MYRANKX + 1, NRPROCX - 1) + 1
      ix = i - MYRANKX + 1
      call ComputeEndpoints(i - 1, NRPROCX, ads % n(1), ads % p(1), ads % nrcpp(1), ads % ibegsx(ix), &
      ads % iendsx(ix), imine, imaxe)
   enddo
   do i = max(MYRANKY - 1, 0) + 1, min(MYRANKY + 1, NRPROCY - 1) + 1
      iy = i - MYRANKY + 1
      call ComputeEndpoints(i - 1, NRPROCY, ads % n(2), ads % p(2), ads % nrcpp(2), ads % ibegsy(iy), &
      ads % iendsy(iy), imine, imaxe)
   enddo
   do i = max(MYRANKZ - 1, 0) + 1, min(MYRANKZ + 1, NRPROCZ - 1) + 1
      iz = i - MYRANKZ + 1
      call ComputeEndpoints(i - 1, NRPROCZ, ads % n(3), ads % p(3), ads % nrcpp(3), ads % ibegsz(iz), &
      ads % iendsz(iz), imine, imaxe)
   enddo

end subroutine ComputeDecomposition


! -------------------------------------------------------------------
! Allocates most of the 'static' arrays
! -------------------------------------------------------------------
subroutine AllocateADSdata(ads, ads_data)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ
   use mpi
   implicit none
   type(ADS_setup), intent(in) :: ads
   type (ADS_compute_data), intent(out) :: ads_data
   integer :: ierr

   allocate(ads_data % Un(ads%lnelem(1),ads%lnelem(2),ads % lnelem(3),ads%ng(1),ads%ng(2),ads%ng(3)))
   allocate(ads_data % Un13(ads%lnelem(1),ads%lnelem(2),ads % lnelem(3),ads%ng(1),ads%ng(2),ads%ng(3)))
   allocate(ads_data % Un23(ads%lnelem(1),ads%lnelem(2),ads % lnelem(3),ads%ng(1),ads%ng(2),ads%ng(3)))
   allocate(ads_data % dUn(ads%lnelem(1),ads%lnelem(2),ads % lnelem(3),ads%ng(1),ads%ng(2),ads%ng(3),3))
   
   ! OLD: MP start with system fully generated along X
   ! allocate( F((n+1),(sy)*(sz))) !x,y,z
   allocate( ads_data % F_test(ads % s(1), ads % s(2) * ads % s(3))) !x,y,z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO CHANGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   allocate( ads_data % F(ads % s(1), ads % s(2) * ads % s(3))) !x,y,z
   allocate(ads_data % F2(ads % s(2), ads % s(1) * ads % s(3))) !y,x,z
   allocate(ads_data % F3(ads % s(3), ads % s(1) * ads % s(2))) !z,x,y

   allocate(ads_data % R(ads % nrcpp(3) * ads % nrcpp(1) * ads % nrcpp(2), 3, 3, 3))
   ads_data % R = 0.d0

   call mpi_barrier(MPI_COMM_WORLD, ierr)

end subroutine AllocateADSdata


! -------------------------------------------------------------------
! Allocates most of the 'static' arrays
! -------------------------------------------------------------------
subroutine AllocateADS(n,nelem,p,ads)
   use Setup, ONLY: ADS_Setup
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ
   use mpi
   implicit none
   integer (kind=4), dimension(3), intent(in) :: n,nelem,p
   type(ADS_setup), intent(out) :: ads
   integer :: ierr

   allocate(ads % Ox(nelem(1)))
   allocate(ads % Oy(nelem(2)))
   allocate(ads % Oz(nelem(3)))
      
   allocate(ads % Jx(nelem(1)))
   allocate(ads % Jy(nelem(2)))
   allocate(ads % Jz(nelem(3)))

   allocate(ads % Xx(p(1) + 1, nelem(1)))
   allocate(ads % Xy(p(2) + 1, nelem(2)))
   allocate(ads % Xz(p(3) + 1, nelem(3)))

   allocate(ads % NNx(0:1, 0:p(1), p(1) + 1, nelem(1)))
   allocate(ads % NNy(0:1, 0:p(2), p(2) + 1, nelem(2)))
   allocate(ads % NNz(0:1, 0:p(3), p(3) + 1, nelem(3)))

   allocate(ads % Wx(p(1) + 1))
   allocate(ads % Wy(p(2) + 1))
   allocate(ads % Wz(p(3) + 1))

   ! Processes on the border need pivot vector for LAPACK call
   if (MYRANKX == 0 .or. MYRANKY == 0 .or. MYRANKZ == 0) then
      allocate(ads % IPIVx(n(1) + 1))
      allocate(ads % IPIVy(n(2) + 1))
      allocate(ads % IPIVz(n(3) + 1))
   endif
end subroutine AllocateADS


!!!!! przeniesc do debug
! -------------------------------------------------------------------
! Prints debugging information about results of distributing
! data to neighbouring processes.
! -------------------------------------------------------------------
subroutine PrintDistributedData(ads, ads_data)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, PRINTRANK, &
   NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints
   implicit none
   type (ADS_setup), intent(inout) :: ads
   type (ADS_compute_data), intent(in) :: ads_data
   integer(kind = 4) :: i, j, k
   integer(kind = 4) :: obegx, oendx, obegy, oendy, obegz, oendz
   integer(kind = 4) :: mine, maxe

   write(*, *) PRINTRANK, 'R:'

   do i = max(MYRANKX - 1, 0) + 1, min(MYRANKX + 1, NRPROCX - 1) + 1
      do j = max(MYRANKY - 1, 0) + 1, min(MYRANKY + 1, NRPROCY - 1) + 1
         do k = max(MYRANKZ - 1, 0) + 1, min(MYRANKZ + 1, NRPROCZ - 1) + 1
            write(*, *) '(i,j,k)=', i + 1, j + 1, k + 1

            call ComputeEndpoints(i - 1, NRPROCX, ads % n(1), ads % p(1), ads % nrcpp(1), obegx, oendx, mine, maxe)
            call ComputeEndpoints(j - 1, NRPROCY, ads % n(2), ads % p(2), ads % nrcpp(2), obegy, oendy, mine, maxe)
            call ComputeEndpoints(k - 1, NRPROCZ, ads % n(3), ads % p(3), ads % nrcpp(3), obegz, oendz, mine, maxe)

            write(*, *) reshape(&
            ads_data % R(:, i - MYRANKX + 1, j - MYRANKY + 1, k - MYRANKZ + 1), &
            (/ oendz - obegz + 1, oendx - obegx + 1, oendy - obegy + 1 /))
         enddo
      enddo
   enddo

   write(*, *) '----'

end subroutine PrintDistributedData


!!!! przeniesc do solver
! -------------------------------------------------------------------
! Solves 1D linear system (one of 3 steps of solving the whole),
! using DGBSV.
!
! RHS   - vector of right-hand sides, of dimension (n+1) x eqnum
! eqnum - number of right-hand sides (equations)
! -------------------------------------------------------------------
subroutine SolveOneDirection(RHS, eqnum, n, p, sprsmtrx)
   use sparse
   use mpi
   implicit none
   include 'dmumps_struc.h'
   real (kind = 8) :: RHS(:,:)
   integer :: n, p
   integer(kind = 4) :: eqnum
   integer(kind = 4) :: i, iret
   type(sparse_matrix), pointer, intent(in) :: sprsmtrx
    type(dmumps_struc) :: mumps_par

   mumps_par % comm = MPI_COMM_SELF
   mumps_par%job = -1
   mumps_par%par = 1
   mumps_par%N = n+1
   call dmumps(mumps_par)

   call to_mumps_format(sprsmtrx, mumps_par)
   allocate(mumps_par%rhs(mumps_par%n))


!     error output stream (non-positive to suppress)
      mumps_par%icntl(1)  = 1 !1
!     diagnostic, statistics and warnings
      mumps_par%icntl(2)  =0! 1 !1
!     global information
      mumps_par%icntl(3)  = 0!6 !6
!     printing level
      mumps_par%icntl(4)  = 0!3 !3
!     input matrix in assembled or element format
      mumps_par%icntl(5)  = 0
!     column permutation for zero-free diagonal (automatic)
!     mumps_par%icntl(6)  = 7
!     pivot order (automatic)
      mumps_par%icntl(7)  = 4 !1 enforce ordering, 5 metis, 0 AMD, 7 auto
!     scaling (automatic)
!     mumps_par%icntl(8)  = 7
!     no transpose
!     mumps_par%icntl(9)  = 1
!     max steps for iterative refinement
!     mumps_par%icntl(10) = 0
!     statistics info
      mumps_par%icntl(11) = 2
!     controls parallelism
      mumps_par%icntl(12) = 0
!     use ScaLAPACK for root node
      mumps_par%icntl(13) = 1 !0 use 1 do not use
!     percentage increase in estimated workspace
      mumps_par%icntl(14) = 50
!     matrix distribution for assembled input
      mumps_par%icntl(18) = 0 !distributed
!     nonzero for Schur complement
      mumps_par%icntl(19) = 0
!     distribution of RHS (centralized on host)
      mumps_par%icntl(20) = 0
!     mumps_par%icntl(32) = 1

      mumps_par%job = 1
      call dmumps(mumps_par)
      
      mumps_par%job = 2
      call dmumps(mumps_par)
      if (mumps_par%info(1).ne.0) then
        write (*,*) 'mumps_par%job=',mumps_par%job
        write (*,*) 'mumps_par%info=',mumps_par%info
        stop 1
      endif

   do i = 1, eqnum
      mumps_par%rhs(1:n+1) = rhs(1:n+1,i)
      mumps_par%job = 3
      call dmumps(mumps_par)
      rhs(1:n+1,i) = mumps_par%rhs(1:n+1)
   enddo
      
      mumps_par%job = -2
      call dmumps(mumps_par)


      deallocate(mumps_par%irn)
      deallocate(mumps_par%jcn)
      deallocate(mumps_par%a)
      deallocate(mumps_par%rhs)
   
end subroutine SolveOneDirection





! -------------------------------------------------------------------
! Performs one step of the simulation with multiple substeps
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
subroutine MultiStep(iter, mix, RHS_fun, ads, ads_data, l2norm, mierr)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use projection_engine, ONLY: FormUn
   use Interfaces, ONLY: RHS_fun_int
   implicit none
   integer(kind = 4), intent(in) :: iter
   real(kind=8), intent(in) :: mix(4,3)
   procedure(RHS_fun_int) :: RHS_fun
   type (ADS_setup), intent(in) :: ads
   type (ADS_compute_data), intent(inout) :: ads_data
   real (kind = 8), intent(out) :: l2norm
   integer(kind = 4), intent(out) :: mierr
   real(kind=8) :: mmix(4)
   integer (kind=4) :: direction
   integer (kind=4) :: substep
   
   mmix = mix(:,1)
   direction = 1
   substep = 1
   call FormUn(1, ads, ads_data)
   ads_data % un13 = 0.d0
   ads_data % un23 = 0.d0
   call Sub_Step(ads, ads, iter, mmix,direction,substep,RHS_fun,ads_data, l2norm, mierr)
   
   mmix = mix(:,2)
   direction = 2
   substep = 2
   call FormUn(2, ads, ads_data)
   ads_data % un23 = 0.d0
   call Sub_Step(ads, ads, iter, mmix,direction,substep,RHS_fun,ads_data, l2norm, mierr)
   
   mmix = mix(:,3)
   direction = 3
   substep = 3
   call FormUn(3, ads, ads_data)
   call Sub_Step(ads, ads, iter, mmix,direction,substep,RHS_fun,ads_data, l2norm, mierr)
   
   
end subroutine MultiStep

! -------------------------------------------------------------------
! Performs one step of the simulation with single substep
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
subroutine Step(iter, RHS_fun, ads, ads_data, l2norm, mierr)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use projection_engine, ONLY: FormUn
   use Interfaces, ONLY: RHS_fun_int
   implicit none
   integer(kind = 4), intent(in) :: iter
   procedure(RHS_fun_int) :: RHS_fun
   type (ADS_setup), intent(in) :: ads
   type (ADS_compute_data), intent(inout) :: ads_data
   real (kind = 8), intent(out) :: l2norm
   integer(kind = 4), intent(out) :: mierr
   real(kind=8) :: mix(4)
   integer (kind=4) :: direction
   integer (kind=4) :: substep
   
   mix = (/ 1.d0, 0.d0, 0.d0, 0.d0 /)
   direction = 1
   substep = 1
   ads_data % un13 = 0.d0
   ads_data % un23 = 0.d0
   call FormUn(1, ads, ads_data)
   
   call Sub_Step(ads, ads, iter, mix,direction,substep,RHS_fun,ads_data, l2norm, mierr)
   
end subroutine Step
   
!!!! podzielic na wraper i czesc wlasciwa
! przeniesc czesc do solver
! -------------------------------------------------------------------
! Performs one substep of the simulation
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
subroutine Sub_Step(ads, ads_trial, iter, mix,direction,substep,RHS_fun,ads_data, l2norm, mierr)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY:PRINTRANK, MYRANKX, MYRANKY, MYRANKZ
   use communicators, ONLY: COMMX, COMMY, COMMZ
   use reorderRHS, ONLY: ReorderRHSForX, ReorderRHSForY, ReorderRHSForZ
   use projection_engine, ONLY: Form3DRHS, ComputeMatrix
   use my_mpi, ONLY: DistributeSpline, Gather, Scatter
   use Interfaces, ONLY: RHS_fun_int
   use mpi
   use sparse
   implicit none
   type (ADS_setup), intent(in) :: ads
   type (ADS_setup), intent(in) :: ads_trial
   integer(kind = 4), intent(in) :: iter
   real(kind=8), intent(in) :: mix(4)
   integer (kind=4), intent(in) :: direction
   integer (kind=4), intent(in) :: substep
   procedure(RHS_fun_int) :: RHS_fun
   type (ADS_compute_data), intent(inout) :: ads_data
   real (kind = 8), intent(out) :: l2norm
   integer(kind = 4), intent(out) :: mierr
   integer(kind = 4) :: i
   integer(kind = 4) :: iret, ierr
   integer(kind = 4), dimension(3) :: nrcpp
   real(kind = 8) :: time1, time2
   type(sparse_matrix), pointer :: sprsmtrx

#ifdef PERFORMANCE
   time1 = MPI_Wtime()
#endif
   ! generate the RHS vectors
   call Form3DRHS(ads, ads_trial, ads_data, direction, substep, RHS_fun, l2norm)
#ifdef PERFORMANCE
   time2 = MPI_Wtime()
   write(*,*) "Form 3D RHS: ", time2 - time1
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO CHANGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ads_data%F = 0.d0
   ads_data%F = ads_data%F_test

#ifdef IPRINT
   write(*, *) PRINTRANK, 'F'
   do i = 1, ads_trial % s(1)
      write(*, *) PRINTRANK, ads % F(i, 1:ads_trial % s(2) * ads_trial % s(3))
   enddo
#endif

   !--------------------------------------------------------------------
   ! Solve the first problem
   !--------------------------------------------------------------------

   
   call solve_problem(ads, ads_trial, 1, 2, 3, &
   mix, sprsmtrx, ads_data % F, ads_data % F2, ierr)

   !--------------------------------------------------------------------
   ! Solve the second problem
   !--------------------------------------------------------------------
   
   call solve_problem(ads, ads_trial, 2, 1, 3, &
   mix, sprsmtrx, ads_data % F2, ads_data % F3, ierr)

   !--------------------------------------------------------------------
   ! Solve the third problem
   !--------------------------------------------------------------------

   call solve_problem(ads, ads_trial, 3, 1, 2, &
   mix, sprsmtrx, ads_data % F3, ads_data % F, ierr)

      
#ifdef IINFO
   write(*, *) PRINTRANK, '3e) DISTRIBUTE SOLUTION'
#endif
   do i=1,ads % s(2) * ads % s(3)
      ads_data % R((i-1)*ads % s(1)+1:i*ads % s(1), 2, 2, 2) = ads_data % F(:,i)
   enddo
   nrcpp = (/ ads % nrcpp(3), ads % nrcpp(1), ads % nrcpp(2) /)
   call DistributeSpline(ads_data % R, nrcpp, ads_data % R)

#ifdef IPRINT
   write(*, *) PRINTRANK, 'Result:'
   do i = 1, ads % s(3)
      write(*, *) PRINTRANK, i, 'row=', ads % F(i,:)
   enddo
#endif
      

   call mpi_barrier(MPI_COMM_WORLD, ierr)

   mierr = 0
end subroutine Sub_Step




! -------------------------------------------------------------------
! Deallocates all the resources and finalizes MPI.
! -------------------------------------------------------------------
subroutine Cleanup(ads, ads_data, mierr)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: PRINTRANK
   implicit none
   type (ADS_setup), intent(inout) :: ads
   type (ADS_compute_data), intent(inout) :: ads_data
   integer(kind = 4), intent(out) :: mierr
   integer(kind = 4) :: ierr

   call Cleanup_ADS(ads, mierr)

   if (allocated(ads_data % F)) deallocate(ads_data % F)
   if (allocated(ads_data % F2)) deallocate(ads_data % F2)
   if (allocated(ads_data % F3)) deallocate(ads_data % F3)

   if (allocated(ads_data % Un)) deallocate(ads_data % Un)
   if (allocated(ads_data % Un13)) deallocate(ads_data % Un13)
   if (allocated(ads_data % Un23)) deallocate(ads_data % Un23)
   if (allocated(ads_data % dUn)) deallocate(ads_data % dUn)
   !!!!!! wyciac
   call mpi_finalize(ierr)
#ifdef IINFO
   write(*, *) PRINTRANK, "Exiting..."
#endif

   mierr = 0

end subroutine Cleanup


! -------------------------------------------------------------------
! Deallocates all the resources and finalizes MPI.
! -------------------------------------------------------------------
subroutine Cleanup_ADS(ads, mierr)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: PRINTRANK
   implicit none
   type (ADS_setup), intent(inout) :: ads
   integer(kind = 4), intent(out) :: mierr
   integer(kind = 4) :: ierr

   if (allocated(ads % shiftsX)) deallocate(ads % shiftsX)
   if (allocated(ads % shiftsY)) deallocate(ads % shiftsY)
   if (allocated(ads % shiftsZ)) deallocate(ads % shiftsZ)

   if (allocated(ads % dimensionsX)) deallocate(ads % dimensionsX)
   if (allocated(ads % dimensionsX)) deallocate(ads % dimensionsY)
   if (allocated(ads % dimensionsZ)) deallocate(ads % dimensionsZ)

   if (allocated(ads % IPIVx)) deallocate(ads % IPIVx)
   if (allocated(ads % IPIVy)) deallocate(ads % IPIVy)
   if (allocated(ads % IPIVz)) deallocate(ads % IPIVz)

   if (allocated(ads % Ux)) deallocate(ads % Ux)
   if (allocated(ads % Uy)) deallocate(ads % Uy)
   if (allocated(ads % Uz)) deallocate(ads % Uz)

   if (allocated(ads % Ox)) deallocate(ads % Ox)
   if (allocated(ads % Oy)) deallocate(ads % Oy)
   if (allocated(ads % Oz)) deallocate(ads % Oz)

   if (allocated(ads % Jx)) deallocate(ads % Jx)
   if (allocated(ads % Jy)) deallocate(ads % Jy)
   if (allocated(ads % Jz)) deallocate(ads % Jz)

   if (allocated(ads % Xx)) deallocate(ads % Xx)
   if (allocated(ads % Xy)) deallocate(ads % Xy)
   if (allocated(ads % Xz)) deallocate(ads % Xz)

   if (allocated(ads % NNx)) deallocate(ads % NNx)
   if (allocated(ads % NNy)) deallocate(ads % NNy)
   if (allocated(ads % NNz)) deallocate(ads % NNz)

   if (allocated(ads % Wx)) deallocate(ads % Wx)
   if (allocated(ads % Wy)) deallocate(ads % Wy)
   if (allocated(ads % Wz)) deallocate(ads % Wz)
   mierr = 0

end subroutine Cleanup_ADS


!!!! przeniesc do debug
! -------------------------------------------------------------------
! Sanity-check of dimensions vector
! -------------------------------------------------------------------
subroutine ValidateDimensions(n, s, nrcpp, &
   dimensionsX, dimensionsY, dimensionsZ)
   use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ, PRINTRANK
   use mpi
   implicit none
   integer(kind = 4), intent(in), dimension(3) :: n
   integer(kind = 4), intent(in), dimension(3) :: s
   integer(kind = 4), intent(in), dimension(3) :: nrcpp
   integer(kind = 4), intent(in), allocatable, dimension(:) :: dimensionsX
   integer(kind = 4), intent(in), allocatable, dimension(:) :: dimensionsY
   integer(kind = 4), intent(in), allocatable, dimension(:) :: dimensionsZ

   integer(kind = 4) :: i, k

   k = 0
   do i = 1, NRPROCX
      k = k + dimensionsX(i)
   enddo
   if (k /= (n(1) + 1) * s(2) * s(3)) then
      write(*, *) PRINTRANK, 'problem with dimensionsX', dimensionsX
      write(*, *) PRINTRANK, 'nx+1', n(1) + 1
      write(*, *) PRINTRANK, 'sy', s(2)
      write(*, *) PRINTRANK, 'sz', s(3)
      write(*, *) PRINTRANK, 'nrcppx', nrcpp(1)
      stop
   endif

   k = 0
   do i = 1, NRPROCY
      k = k + dimensionsY(i)
   enddo
   if (k /= (n(2) + 1) * s(1) * s(3)) then
      write(*, *) PRINTRANK, 'problem with dimensionsY', dimensionsY
      write(*, *) PRINTRANK, 'n+1', n(2) + 1
      write(*, *) PRINTRANK, 'sx', s(1)
      write(*, *) PRINTRANK, 'sz', s(3)
      stop
   endif

   k = 0
   do i = 1, NRPROCZ
      k = k + dimensionsZ(i)
   enddo
   if (k /= (n(3) + 1) * s(1) * s(2)) then
      write(*, *) PRINTRANK, 'problem with dimensionsZ', dimensionsZ
      write(*, *) PRINTRANK, 'n+1', n(3) + 1
      write(*, *) PRINTRANK, 'sx', s(1)
      write(*, *) PRINTRANK, 'sy', s(2)
      stop
   endif

end subroutine ValidateDimensions



!!!!! przeniesc do debug
! -------------------------------------------------------------------
! Displays computed domain decomposition, for debugging.
! -------------------------------------------------------------------
subroutine PrintDecompositionInfo(n, nrcpp, ibeg, iend)
   use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ, PRINTRANK, &
   MYRANKX, MYRANKY, MYRANKZ
   implicit none
   integer(kind = 4), intent(in), dimension(3) :: n
   integer(kind = 4), intent(in), dimension(3) :: nrcpp
   integer(kind = 4), intent(in), dimension(3) :: ibeg
   integer(kind = 4), intent(in), dimension(3) :: iend

   write(*, *) PRINTRANK, 'MYRANKX,MYRANKY,MYRANKZ', MYRANKX, MYRANKY, MYRANKZ
   write(*, *) PRINTRANK, 'NRPROCX,NRPROCY,NRPROCZ', NRPROCX, NRPROCY, NRPROCZ
   write(*, *) PRINTRANK, 'nx+1', n(1) + 1
   write(*, *) PRINTRANK, 'ny+1', n(2) + 1
   write(*, *) PRINTRANK, 'nz+1', n(3) + 1
   write(*, *) PRINTRANK, 'nrcppx,nrcppy,nrcppz', nrcpp(1), nrcpp(2), nrcpp(3)
   write(*, *) PRINTRANK, 'ibegx,iendx', ibeg(1), iend(1)
   write(*, *) PRINTRANK, 'ibegy,iendy', ibeg(2), iend(2)
   write(*, *) PRINTRANK, 'ibegz,iendz', ibeg(3), iend(3)

end subroutine PrintDecompositionInfo


! -------------------------------------------------------------------
! Gathers full solution and plots it
! -------------------------------------------------------------------
subroutine PrintSolution(iter, t, ads, ads_data)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK
   use plot, ONLY: SaveSplinePlot, PlotParams
   use vtk, ONLY: VtkOutput
   use my_mpi, ONLY: GatherFullSolution
   implicit none
   type (ADS_compute_data), intent(in) :: ads_data
   type (ADS_setup), intent(in) :: ads
   integer(kind = 4), intent(in) :: iter
   real (kind = 8), intent(in) :: t
   real (kind = 8), allocatable :: solution(:,:,:)
   type (PlotParams) :: params
   character(len = 20) :: filename

   call GatherFullSolution(0, ads_data % F, solution, &
   ads % n, ads % p, ads % s)

   if (MYRANK == 0) then
      write(filename, '(I10)') iter
      filename = 'step' // adjustl(filename)
      ! filename = trim(filename) // '_'

      params = PlotParams(0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 31, 31, 31)
      call SaveSplinePlot(trim(filename), &
      ads % Ux, ads % p(1), ads % n(1), ads % nelem(1), &
      ads % Uy, ads % p(2), ads % n(2), ads % nelem(2), &
      ads % Uz, ads % p(3), ads % n(3), ads % nelem(3), &
      ! solution, GnuPlotOutput, params)
      solution, VtkOutput, params)

      ! call SavePlot(trim(filename), ftest, GnuPlotOutput, params)
   endif

end subroutine PrintSolution







subroutine solve_problem(ads, ads_trial, a, b, c, mix, sprsmtrx, F, F2, ierr)
   use Setup, ONLY: ADS_Setup
   use sparse
   use mpi
   use my_mpi, ONLY: Gather, Scatter
   use parallelism, ONLY:PRINTRANK, MYRANKX, MYRANKY, MYRANKZ
   use communicators, ONLY: COMMX, COMMY, COMMZ
   use reorderRHS, ONLY: ReorderRHSForX, ReorderRHSForY, ReorderRHSForZ
   use projection_engine, ONLY: ComputeMatrix
   implicit none
   type (ADS_setup), intent(in) :: ads,ads_trial
   integer(kind=4), intent(in) :: a,b,c
   integer(kind = 4), dimension(:), allocatable :: dimensions ! Size of slices of domain in each dimension
   integer(kind = 4), allocatable, dimension(:) :: shifts
   real (kind = 8), allocatable, dimension(:,:) :: F,F2
   integer(kind = 4) :: comm
   real(kind=8), intent(in) :: mix(4)
   real (kind = 8), allocatable, dimension(:) :: U
   integer(kind = 4) :: myrankdim ! Integer coordinates of processor along X, Y or Z
   type(sparse_matrix), pointer, intent(inout) :: sprsmtrx
   integer(kind=4), intent(out) :: ierr
   real (kind = 8), allocatable, dimension(:,:) :: F_out, F2_out
   real(kind = 8) :: time1, time2

   if (a.EQ.1) then
     comm = COMMX
     myrankdim = MYRANKX
     u = ads % ux
     shifts = ads % shiftsX
     dimensions = ads % dimensionsX
   else if (a .EQ. 2) then
     comm = COMMY
     myrankdim = MYRANKY
     u = ads % uy
     shifts = ads % shiftsY
     dimensions = ads % dimensionsY
   else
     comm = COMMZ
     myrankdim = MYRANKZ
     u = ads % uz
     shifts = ads % shiftsZ
     dimensions = ads % dimensionsZ
   endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, a,'a) GATHER'
#endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)
   allocate(F_out((ads % n(a) + 1), ads % s(b) * ads % s(c)))
#ifdef PERFORMANCE
   time1 = MPI_Wtime()
#endif
   call Gather(F, F_out, ads % n(a), ads % s(a), ads % s(b) * ads % s(c), &
   dimensions, shifts, comm, ierr)
#ifdef PERFORMANCE
   time2 = MPI_Wtime()
   write(*,*) "Gather", a, " : ", time2 - time1
#endif
   call mpi_barrier(MPI_COMM_WORLD, ierr)

   if (myrankdim == 0) then
#ifdef IINFO
      write(*, *) PRINTRANK, a,'b) SOLVE THE ',a,' PROBLEM'
#endif

#ifdef PERFORMANCE
      time1 = MPI_Wtime()
#endif
      call ComputeMatrix(U, ads % p(a), ads % n(a), ads % nelem(a),&
      U, ads % p(a), ads % n(a), ads % nelem(a), &
      mix, sprsmtrx)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write(*,*) "Mass matrix", a, ": ", time2 - time1
      time1 = MPI_Wtime()
#endif
      call SolveOneDirection(F_out, ads % s(b) * ads % s(c), ads % n(a), ads % p(a), sprsmtrx)
      call clear_matrix(sprsmtrx)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write(*,*) "Solve ", a,": ", time2 - time1
#endif
   endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, a,'c) SCATTER'
#endif

   allocate(F2_out(ads % s(a), ads % s(b) * ads % s(c)))
#ifdef PERFORMANCE
   time1 = MPI_Wtime()
#endif
   call Scatter(F_out, F2_out, ads % n(a), ads % s(a), ads % s(b) * ads % s(c), &
   dimensions, shifts, comm, ierr)
#ifdef PERFORMANCE
   time2 = MPI_Wtime()
   write(*,*) "Scatter ", a,": ", time2 - time1
#endif
   deallocate(F_out)

   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, a,'d) REORDER'
#endif
   ! Reorder right hand sides
   if (a .EQ. 1) call ReorderRHSForY(ads % ibeg, ads % iend, F2_out, F2)
   if (a .EQ. 2) call ReorderRHSForZ(ads % ibeg, ads % iend, F2_out, F2)
   if (a .EQ. 3) call ReorderRHSForX(ads % ibeg, ads % iend, F2_out, F2)
   deallocate(F2_out)

#ifdef IPRINT
   write(*, *) PRINTRANK, 'after ReorderRHS'
   write(*, *) PRINTRANK, 'F:'
   do i = 1, ads % s(1)
      write(*, *) PRINTRANK, i, 'row=', F2(i, 1:ads % s(2) * ads % s(3))
   enddo
#endif
      
end subroutine solve_problem 



end module ADSS
