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

   ads % KL = p
   ads % KU = p

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

   call AllocateArrays(ads, ads_data)

   call PrepareKnot(ads % Ux, n(1), p(1), ads % nelem(1))
   call PrepareKnot(ads % Uy, n(2), p(2), ads % nelem(2))
   call PrepareKnot(ads % Uz, n(3), p(3), ads % nelem(3))

   mierr = 0
   
   ads % m(1) = ads % n(1) + ads % p(1) + 1
   ads % ng(1) = ads % p(1) + 1
   ads % m(2) = ads % n(2) + ads % p(2) + 1
   ads % ng(2) = ads % p(2) + 1
   ads % m(3) = ads % n(3) + ads % p(3) + 1
   ads % ng(3) = ads % p(3) + 1
   
   call AllocateStatic(ads, ads_data)
   
   call BasisData(ads % p(1), ads % m(1), ads % Ux, 1, ads % ng(1), &
   ads % nelem(1), ads % Ox, ads % Jx, ads % Wx, ads % Xx, ads % NNx)
   call BasisData(ads % p(2), ads % m(2), ads % Uy, 1, ads % ng(2), &
   ads % nelem(2), ads % Oy, ads % Jy, ads % Wy, ads % Xy, ads % NNy)
   call BasisData(ads % p(3), ads % m(3), ads % Uz, 1, ads % ng(3), &
   ads % nelem(3), ads % Oz, ads % Jz, ads % Wz, ads % Xz, ads % NNz)

   ads % lnelem(1) = ads % maxe(1) - ads % mine(1) + 1
   ads % lnelem(2) = ads % maxe(2) - ads % mine(2) + 1
   ads % lnelem(3) = ads % maxe(3) - ads % mine(3) + 1
   
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
subroutine AllocateArrays(ads, ads_data)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ
   use mpi
   implicit none
   type(ADS_setup), intent(inout) :: ads
   type (ADS_compute_data), intent(inout) :: ads_data
   integer :: ierr

   allocate(ads_data % Mx(2 * ads % KL(1) + ads % KU(1) + 1, ads % n(1) + 1))
   allocate(ads_data % My(2 * ads % KL(2) + ads % KU(2) + 1, ads % n(2) + 1))
   allocate(ads_data % Mz(2 * ads % KL(3) + ads % KU(3) + 1, ads % n(3) + 1))

   ! OLD: MP start with system fully generated along X
   ! allocate( F((n+1),(sy)*(sz))) !x,y,z
   allocate( ads_data % F(ads % s(1), ads % s(2) * ads % s(3))) !x,y,z
   allocate(ads_data % F2(ads % s(2), ads % s(1) * ads % s(3))) !y,x,z
   allocate(ads_data % F3(ads % s(3), ads % s(1) * ads % s(2))) !z,x,y

   ! Processes on the border need pivot vector for LAPACK call
   if (MYRANKX == 0 .or. MYRANKY == 0 .or. MYRANKZ == 0) then
      allocate(ads % IPIVx(ads % n(1) + 1))
      allocate(ads % IPIVy(ads % n(2) + 1))
      allocate(ads % IPIVz(ads % n(3) + 1))
   endif

   allocate(ads_data % R(ads % nrcpp(3) * ads % nrcpp(1) * ads % nrcpp(2), 3, 3, 3))
   ads_data % R = 0.d0

   call mpi_barrier(MPI_COMM_WORLD, ierr)

end subroutine AllocateArrays


! -------------------------------------------------------------------
! Allocates most of the 'static' arrays
! -------------------------------------------------------------------
subroutine AllocateStatic(ads, ads_data)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ
   use mpi
   implicit none
   type(ADS_setup), intent(inout) :: ads
   type (ADS_compute_data), intent(inout) :: ads_data
   integer :: ierr

   allocate(ads % Ox(ads % nelem(1)))
   allocate(ads % Oy(ads % nelem(2)))
   allocate(ads % Oz(ads % nelem(3)))
      
   allocate(ads % Jx(ads % nelem(1)))
   allocate(ads % Jy(ads % nelem(2)))
   allocate(ads % Jz(ads % nelem(3)))

   allocate(ads % Xx(ads % p(1) + 1, ads % nelem(1)))
   allocate(ads % Xy(ads % p(2) + 1, ads % nelem(2)))
   allocate(ads % Xz(ads % p(3) + 1, ads % nelem(3)))

   allocate(ads % NNx(0:1, 0:ads % p(1), ads % p(1) + 1, ads % nelem(1)))
   allocate(ads % NNy(0:1, 0:ads % p(2), ads % p(2) + 1, ads % nelem(2)))
   allocate(ads % NNz(0:1, 0:ads % p(3), ads % p(3) + 1, ads % nelem(3)))

   allocate(ads % Wx(ads % p(1) + 1))
   allocate(ads % Wy(ads % p(2) + 1))
   allocate(ads % Wz(ads % p(3) + 1))

end subroutine AllocateStatic


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
subroutine SolveOneDirection(RHS, eqnum, n, KL, KU, p, M, IPIV)
   implicit none
   real (kind = 8) :: RHS(:,:)
   integer :: KL, KU
   integer, dimension(:) :: IPIV
   real (kind = 8), dimension(:,:) :: M
   integer :: n, p
   integer(kind = 4) :: eqnum
   integer(kind = 4) :: i, iret

   IPIV = 0

#ifdef IPRINT
   write(*, *) 'CALL DGBSV'
   write(*, *) 'N=', n + 1
   write(*, *) 'KL=', KL
   write(*, *) 'KU=', KU
   write(*, *) 'NRHS=', eqnum
   write(*, *) 'AB='
   do i = 1, 2 * KL + KU + 1
      write(*, *) i, 'row=', M(i, 1:n + 1)
   enddo
   write(*, *) 'LDAB=', 2 * KL + KU + 1
   write(*, *) 'IPIV=', IPIV
   write(*, *) 'B='
   do i = 1, n + 1
      write(*, *) i, 'row=', RHS(i, 1:eqnum)
   enddo
   write(*, *) 'LDB=', n + 1
#endif
   
   ! SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
   ! .. Scalar Arguments ..
   ! INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
   ! .. Array Arguments ..
   ! INTEGER            IPIV( * )
   ! DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )

   call DGBSV(n + 1, KL, KU, eqnum, M, 2 * KL + KU + 1, IPIV, RHS, n + 1, iret)

#ifdef IPRINT
   write(*, *) 'iret=', iret
   write(*, *) 'Solution='
   do i = 1, n + 1
      write(*, *) i, 'row=', RHS(i, 1:eqnum)
   enddo
#endif
   
end subroutine SolveOneDirection





! -------------------------------------------------------------------
! Performs one step of the simulation with multiple substeps
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
subroutine MultiStep(iter, mix, RHS_fun, ads, ads_data, l2norm, mierr)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY:PRINTRANK, MYRANKX, MYRANKY, MYRANKZ
   use communicators, ONLY: COMMX, COMMY, COMMZ
   use reorderRHS, ONLY: ReorderRHSForX, ReorderRHSForY, ReorderRHSForZ
   use projection_engine, ONLY: Form3DRHS, ComputeMatrix
   use my_mpi, ONLY: DistributeSpline, Gather, Scatter
   use Interfaces, ONLY: RHS_fun_int
   use mpi
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
   real (kind = 8) :: un13,un23
   
   mmix = mix(:,1)
   direction = 1
   substep = 1
   un13 = 0.d0
   un23 = 0.d0
   call Sub_Step(iter, mmix,direction,substep,un13,un23, RHS_fun, ads, ads_data, l2norm, mierr)
   
   mmix = mix(:,2)
   direction = 2
   substep = 2
   un13 = 0.d0
   un23 = 0.d0
   call Sub_Step(iter, mmix,direction,substep,un13,un23, RHS_fun, ads, ads_data, l2norm, mierr)
   
   mmix = mix(:,3)
   direction = 3
   substep = 3
   un13 = 0.d0
   un23 = 0.d0
   call Sub_Step(iter, mmix,direction,substep,un13,un23, RHS_fun, ads, ads_data, l2norm, mierr)
end subroutine MultiStep

! -------------------------------------------------------------------
! Performs one step of the simulation with single substep
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
subroutine Step(iter, RHS_fun, ads, ads_data, l2norm, mierr)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY:PRINTRANK, MYRANKX, MYRANKY, MYRANKZ
   use communicators, ONLY: COMMX, COMMY, COMMZ
   use reorderRHS, ONLY: ReorderRHSForX, ReorderRHSForY, ReorderRHSForZ
   use projection_engine, ONLY: Form3DRHS, ComputeMatrix
   use my_mpi, ONLY: DistributeSpline, Gather, Scatter
   use Interfaces, ONLY: RHS_fun_int
   use mpi
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
   real (kind = 8) :: un13,un23
   
   mix = (/ 1.d0, 0.d0, 0.d0, 0.d0 /)
   direction = 1
   substep = 1
   un13 = 0.d0
   un23 = 0.d0
   call Sub_Step(iter, mix,direction,substep,un13,un23, RHS_fun, ads, ads_data, l2norm, mierr)
end subroutine Step
   
!!!! podzielic na wraper i czesc wlasciwa
! przeniesc czesc do solver
! -------------------------------------------------------------------
! Performs one substep of the simulation
!
! iter - number of the iteration
! t    - time at the beginning of step
! -------------------------------------------------------------------
subroutine Sub_Step(iter, mix,direction,substep,un13,un23,RHS_fun, ads, ads_data, l2norm, mierr)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY:PRINTRANK, MYRANKX, MYRANKY, MYRANKZ
   use communicators, ONLY: COMMX, COMMY, COMMZ
   use reorderRHS, ONLY: ReorderRHSForX, ReorderRHSForY, ReorderRHSForZ
   use projection_engine, ONLY: Form3DRHS, ComputeMatrix
   use my_mpi, ONLY: DistributeSpline, Gather, Scatter
   use Interfaces, ONLY: RHS_fun_int
   use mpi
   implicit none
   integer(kind = 4), intent(in) :: iter
   real(kind=8), intent(in) :: mix(4)
   integer (kind=4), intent(in) :: direction
   integer (kind=4), intent(in) :: substep
   real (kind = 8), intent(in) :: un13,un23
   procedure(RHS_fun_int) :: RHS_fun
   type (ADS_setup), intent(in) :: ads
   type (ADS_compute_data), intent(inout) :: ads_data
   real (kind = 8), intent(out) :: l2norm
   integer(kind = 4), intent(out) :: mierr
   integer(kind = 4) :: i
   integer(kind = 4) :: iret, ierr
   integer(kind = 4), dimension(3) :: nrcpp
   real(kind = 8) :: time1, time2

#ifdef PERFORMANCE
   time1 = MPI_Wtime()
#endif
   ! generate the RHS vectors
   call Form3DRHS(ads, ads_data, direction, substep,un13,un23,RHS_fun,l2norm)
#ifdef PERFORMANCE
   time2 = MPI_Wtime()
   write(*,*) "Form 3D RHS: ", time2 - time1
#endif

#ifdef IPRINT
   write(*, *) PRINTRANK, 'F'
   do i = 1, ads % s(1)
      write(*, *) PRINTRANK, ads % F(i, 1:ads % s(2) * ads % s(3))
   enddo
#endif

   !--------------------------------------------------------------------
   ! Solve the first problem
   !--------------------------------------------------------------------
   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, '1a) GATHER'
#endif

   allocate(ads_data % F_out((ads % n(1) + 1), ads % s(2) * ads % s(3)))

#ifdef PERFORMANCE
   time1 = MPI_Wtime()
#endif
   call Gather(ads_data % F, ads_data % F_out, ads % n(1), ads % s(1), ads % s(2) &
   *ads % s(3), ads % dimensionsX, ads % shiftsX, COMMX, ierr)
#ifdef PERFORMANCE
   time2 = MPI_Wtime()
   write(*,*) "Gather 1: ", time2 - time1
#endif

#ifdef IPRINT
   write(*, *) PRINTRANK, 'after call mpi_gather'
   write(*, *) PRINTRANK, 'ierr', ierr
   write(*, *) PRINTRANK, 'F_out:'
   do i = 1, ads % n(1) + 1
      write(*, *) PRINTRANK, i, 'row=', ads % F_out(i, 1:ads % s(2) * ads % s(3))
   enddo
#endif
   call mpi_barrier(MPI_COMM_WORLD, ierr)

   if (MYRANKX == 0) then
#ifdef IINFO
      write(*, *) PRINTRANK, '1b) SOLVE THE FIRST PROBLEM'
#endif

#ifdef PERFORMANCE
      time1 = MPI_Wtime()
#endif
      call ComputeMatrix(ads % KL(1), ads % KU(1), ads % Ux, ads % p(1), &
      ads % n(1), ads % nelem(1), mix, ads_data % Mx)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write(*,*) "Mass matrix 1: ", time2 - time1
      time1 = MPI_Wtime()
#endif
      call SolveOneDirection(ads_data % F_out, ads % s(2) * ads % s(3), ads % n(1), &
      ads % KL(1), ads % KU(1), ads % p(1), ads_data % Mx, ads % IPIVx)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write(*,*) "Solve 1: ", time2 - time1
#endif
   endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, '1c) SCATTER'
#endif
   allocate(ads_data % F2_out(ads % s(1), ads % s(2) * ads % s(3)))
#ifdef PERFORMANCE
   time1 = MPI_Wtime()
#endif
   call Scatter(ads_data % F_out, ads_data % F2_out, ads % n(1), ads % s(1), ads % s(2) * &
   ads % s(3), ads % dimensionsX, ads % shiftsX, COMMX, ierr)
#ifdef PERFORMANCE
   time2 = MPI_Wtime()
   write(*,*) "Scatter 1: ", time2 - time1
#endif
   deallocate(ads_data % F_out)

   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, '1d) REORDER'
#endif
   call ReorderRHSForY(ads % ibeg, ads % iend, ads_data % F2_out, ads_data % F2)
   deallocate(ads_data % F2_out)

#ifdef IPRINT
   write(*, *) PRINTRANK, 'after ReorderRHSForY'
   write(*, *) PRINTRANK, 'F2:'
   do i = 1, ads % s(2)
      write(*, *) PRINTRANK, i, 'row=', ads_data % F2(i, 1:ads % s(1) * ads % s(3))
   enddo
#endif

   !--------------------------------------------------------------------
   ! Solve the second problem
   !--------------------------------------------------------------------
   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, '2a) GATHER'
#endif

   allocate(ads_data % F2_out((ads % n(2) + 1), ads % s(1) * ads % s(3)))
#ifdef PERFORMANCE
   time1 = MPI_Wtime()
#endif
   call Gather(ads_data % F2, ads_data % F2_out, ads % n(2), ads % s(2), ads % s(1) * ads % s(3), &
   ads % dimensionsY, ads % shiftsY, COMMY, ierr)
#ifdef PERFORMANCE
   time2 = MPI_Wtime()
   write(*,*) "Gather 2: ", time2 - time1
#endif), intent

   call mpi_barrier(MPI_COMM_WORLD, ierr)

   if (MYRANKY == 0) then
#ifdef IINFO
      write(*, *) PRINTRANK, '2b) SOLVE THE SECOND PROBLEM'
#endif

#ifdef PERFORMANCE
      time1 = MPI_Wtime()
#endif
      call ComputeMatrix(ads % KL(2), ads % KU(2), ads % Uy, ads % p(2), ads % n(2), &
      ads % nelem(2), mix, ads_data % My)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write(*,*) "Mass matrix 2: ", time2 - time1
      time1 = MPI_Wtime()
#endif
      call SolveOneDirection(ads_data % F2_out, ads % s(1) * ads % s(3), ads % n(2), ads % KL(2), &
      ads % KU(2), ads % p(2), ads_data % My, ads % IPIVy)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write(*,*) "Solve 2: ", time2 - time1
#endif

   endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, '2c) SCATHER'
#endif

   ! CORRECTION
   allocate(ads_data % F3_out(ads % s(2), ads % s(1) * ads % s(3)))
#ifdef PERFORMANCE
   time1 = MPI_Wtime()
#endif
   call Scatter(ads_data % F2_out, ads_data % F3_out, ads % n(2), ads % s(2), ads % s(1) * ads % s(3), &
   ads % dimensionsY, ads % shiftsY, COMMY, ierr)
#ifdef PERFORMANCE
   time2 = MPI_Wtime()
   write(*,*) "Scatter 2: ", time2 - time1
#endif
   deallocate(ads_data % F2_out)

#ifdef IPRINT
   write(*, *) PRINTRANK, 'after call mpi_scatterv'
   write(*, *) PRINTRANK, 'ierr', ierr
   write(*, *) PRINTRANK, 'F3_out:'
   do i = 1, ads % s(2)
      write(*, *) PRINTRANK, i, 'row=', ads % F3_out(i, 1:ads % s(1) * ads % s(3))
   enddo
#endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, '2d) REORDER'
#endif
   ! Reorder right hand sides
   call ReorderRHSForZ(ads % ibeg, ads % iend, ads_data % F3_out, ads_data % F3)
   deallocate(ads_data % F3_out)

#ifdef IPRINT
   write(*, *) PRINTRANK, 'after ReorderRHSForZ'
   write(*, *) PRINTRANK, 'F3:'
   do i = 1, ads % s(3)
      write(*, *) PRINTRANK, i, 'row=', ads_data % F3(i, 1:ads % s(1) * ads % s(2))
   enddo
#endif

   !--------------------------------------------------------------------
   ! Solve the third problem
   !--------------------------------------------------------------------
#ifdef IINFO
   write(*, *) PRINTRANK, '3a) GATHER'
#endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)
   allocate(ads_data % F3_out((ads % n(3) + 1), ads % s(1) * ads % s(2)))
#ifdef PERFORMANCE
   time1 = MPI_Wtime()
#endif
   call Gather(ads_data % F3, ads_data % F3_out, ads % n(3), ads % s(3), ads % s(1) * ads % s(2), &
   ads % dimensionsZ, ads % shiftsZ, COMMZ, ierr)
#ifdef PERFORMANCE
   time2 = MPI_Wtime()
   write(*,*) "Gather 3: ", time2 - time1
#endif
   call mpi_barrier(MPI_COMM_WORLD, ierr)

   if (MYRANKZ == 0) then
#ifdef IINFO
      write(*, *) PRINTRANK, '3b) SOLVE THE THIRD PROBLEM'
#endif

#ifdef PERFORMANCE
      time1 = MPI_Wtime()
#endif
      call ComputeMatrix(ads % KL(3), ads % KU(3), ads % Uz, ads % p(3), ads % n(3), &
      ads % nelem(3), mix, ads_data % Mz)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write(*,*) "Mass matrix 3: ", time2 - time1
      time1 = MPI_Wtime()
#endif
      call SolveOneDirection(ads_data % F3_out, ads % s(1) * ads % s(2), ads % n(3), ads % KL(3), &
      ads % KU(3), ads % p(3), ads_data % Mz, ads % IPIVz)
#ifdef PERFORMANCE
      time2 = MPI_Wtime()
      write(*,*) "Solve 3: ", time2 - time1
#endif
   endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, '3c) SCATTER'
#endif

   ! CORRECTION
   allocate(ads_data % F_out(ads % s(3), ads % s(1) * ads % s(2)))
#ifdef PERFORMANCE
   time1 = MPI_Wtime()
#endif
   call Scatter(ads_data % F3_out, ads_data % F_out, ads % n(3), ads % s(3), ads % s(1) * ads % s(2), &
   ads % dimensionsZ, ads % shiftsZ, COMMZ, ierr)
#ifdef PERFORMANCE
   time2 = MPI_Wtime()
   write(*,*) "Scatter 3: ", time2 - time1
#endif
   deallocate(ads_data % F3_out)

   call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
   write(*, *) PRINTRANK, '3d) REORDER'
#endif
   ! Reorder right hand sides
   call ReorderRHSForX(ads % ibeg, ads % iend, ads_data % F_out, ads_data % F)
   deallocate(ads_data % F_out)

#ifdef IPRINT
   write(*, *) PRINTRANK, 'after ReorderRHSForX'
   write(*, *) PRINTRANK, 'F:'
   do i = 1, ads % s(1)
      write(*, *) PRINTRANK, i, 'row=', ads % F(i, 1:ads % s(2) * ads % s(3))
   enddo
#endif
      
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

   if (allocated(ads_data % Mx)) deallocate(ads_data % Mx)
   if (allocated(ads_data % My)) deallocate(ads_data % My)
   if (allocated(ads_data % Mz)) deallocate(ads_data % Mz)

   if (allocated(ads_data % F)) deallocate(ads_data % F)
   if (allocated(ads_data % F2)) deallocate(ads_data % F2)
   if (allocated(ads_data % F3)) deallocate(ads_data % F3)

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

   !!!!!! wyciac
   call mpi_finalize(ierr)
#ifdef IINFO
   write(*, *) PRINTRANK, "Exiting..."
#endif

   mierr = 0

end subroutine Cleanup



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
   integer (kind = 4), dimension(3) :: tmp1, tmp2, tmp3

   tmp1 = (/ads % n(1), ads % n(2), ads % n(3)/)
   tmp2 = (/ads % p(1), ads % p(2), ads % p(3)/)
   tmp3 = (/ads % s(1), ads % s(2), ads % s(3)/)

   call GatherFullSolution(0, ads_data % F, solution, &
   tmp1, tmp2, tmp3)

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

end module ADSS
