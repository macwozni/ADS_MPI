module ADSS

   type ADS_setup
      ! Number of functions in each dimension minus one
      integer(kind = 4), dimension(3) :: n

      ! Degree of approximation
      integer(kind = 4), dimension(3) :: p

      ! Knot vector
      real (kind = 8), allocatable, dimension(:) :: Ux
      real (kind = 8), allocatable, dimension(:) :: Uy
      real (kind = 8), allocatable, dimension(:) :: Uz

      ! Mass matrix
      real (kind = 8), allocatable, dimension(:,:) :: Mx
      real (kind = 8), allocatable, dimension(:,:) :: My
      real (kind = 8), allocatable, dimension(:,:) :: Mz

      real (kind = 8), allocatable, dimension(:,:) :: F, F2, F3
      real (kind = 8), allocatable, dimension(:,:) :: F_out, F2_out, F3_out

      ! Buffer for coefficients of solution corresponding to neighbouring
      ! parts of the domain. It is (Nx*Ny*Nz) x 3 x 3 x 3 array, where
      ! Nx*Ny*Nz is the size of part of solution for one fragment of domain.
      real (kind = 8), allocatable :: R(:,:,:,:)

      ! Number of subintervals (currently n - p + 1)
      integer(kind = 4), dimension(3) :: nelem

      ! Size of slices of domain in each dimension
      integer(kind = 4), allocatable, dimension(:) :: dimensionsX
      integer(kind = 4), allocatable, dimension(:) :: dimensionsY
      integer(kind = 4), allocatable, dimension(:) :: dimensionsZ

      ! Offsets of slices of domain in each direction
      integer(kind = 4), allocatable, dimension(:) :: shiftsX
      integer(kind = 4), allocatable, dimension(:) :: shiftsY
      integer(kind = 4), allocatable, dimension(:) :: shiftsZ

      ! Pivot array for these processes that need to solve systems
      integer(kind = 4), allocatable, dimension(:) :: IPIVx
      integer(kind = 4), allocatable, dimension(:) :: IPIVy
      integer(kind = 4), allocatable, dimension(:) :: IPIVz

      ! Number of lower and upper diagonal entries in mass matrix
      integer(kind = 4), dimension(3) :: KL, KU

      ! Number of columns (average) per processor
      integer(kind = 4), dimension(3) :: nrcpp

      ! Range of piece of domain assigned to this process
      integer(kind = 4), dimension(3) :: ibeg, iend

      ! Size of piece of domain assigned to this process
      integer(kind = 4), dimension(3) :: s

      ! Ranges of pieces of domain around the one assigned to this process
      integer(kind = 4), dimension(3) :: ibegsx, iendsx
      integer(kind = 4), dimension(3) :: ibegsy, iendsy
      integer(kind = 4), dimension(3) :: ibegsz, iendsz

      ! Range of elements associated with basis functions assigned to this process
      integer(kind = 4), dimension(3) :: mine, maxe
   end type

contains




   ! -------------------------------------------------------------------
   ! Initialization of clocks and MPI
   ! -------------------------------------------------------------------
   subroutine initialize(n, p, ads, mierr)
      use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ
      use parallelism, ONLY: PRINTRANK
      use knot_vector, ONLY: PrepareKnot
      use mpi
      implicit none
      integer(kind = 4), intent(in), dimension(3) :: n
      integer(kind = 4), intent(in), dimension(3) :: p
      type(ADS_setup), intent(out) :: ads
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

   call ComputeDecomposition(&
   ads % n, &
   ads % p, &
   ads % s, &
   ads % ibeg, &
   ads % iend, &
   ads % ibegsx, ads % ibegsy, ads % ibegsz, &
   ads % iendsx, ads % iendsy, ads % iendsz, &
   ads % mine, &
   ads % maxe, &
   ads % nrcpp, &
   ads % dimensionsX, ads % dimensionsY, ads % dimensionsZ, &
   ads % shiftsX, ads % shiftsY, ads % shiftsZ)

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

   call AllocateArrays(&
   ads % n, &
   ads % s, &
   ads % nrcpp, &
   ads % Kl, &
   ads % KU, &
   ads % Mx, ads % My, ads % Mz, &
   ads % F, ads % F2, ads % F3, &
   ads % IPIVx, ads % IPIVy, ads % IPIVz, &
   ads % R)

   call PrepareKnot(ads % Ux, n(1), p(1), ads % nelem(1))
   call PrepareKnot(ads % Uy, n(2), p(2), ads % nelem(2))
   call PrepareKnot(ads % Uz, n(3), p(3), ads % nelem(3))

   mierr = 0
end subroutine


! -------------------------------------------------------------------
! Establishes decomposition of the domain. Calculates size and location
! of the piece for current process.
! -------------------------------------------------------------------
subroutine ComputeDecomposition(&
   n, &
   p, &
   s, &
   ibeg, &
   iend, &
   ibegsx, ibegsy, ibegsz, &
   iendsx, iendsy, iendsz, &
   mine, &
   maxe, &
   nrcpp, &
   dimensionsX, dimensionsY, dimensionsZ, &
   shiftsX, shiftsY, shiftsZ)
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, PRINTRANK, &
   NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints, FillDimVector
   implicit none
   integer(kind = 4), intent(in), dimension(3) :: n
   integer(kind = 4), intent(in), dimension(3) :: p
   integer(kind = 4), intent(out), dimension(3) :: ibeg, iend
   integer(kind = 4), intent(out), dimension(3) :: mine, maxe
   integer(kind = 4), intent(out), dimension(3) :: nrcpp
   integer(kind = 4), intent(out), dimension(3) :: s
   integer(kind = 4), intent(out), allocatable, dimension(:) :: dimensionsX
   integer(kind = 4), intent(out), allocatable, dimension(:) :: dimensionsY
   integer(kind = 4), intent(out), allocatable, dimension(:) :: dimensionsZ
   integer(kind = 4), intent(out), allocatable, dimension(:) :: shiftsX
   integer(kind = 4), intent(out), allocatable, dimension(:) :: shiftsY
   integer(kind = 4), intent(out), allocatable, dimension(:) :: shiftsZ
   integer(kind = 4), intent(out), dimension(3) :: ibegsx, iendsx
   integer(kind = 4), intent(out), dimension(3) :: ibegsy, iendsy
   integer(kind = 4), intent(out), dimension(3) :: ibegsz, iendsz

   integer(kind = 4) :: i
   integer(kind = 4) :: ix, iy, iz
   integer(kind = 4) :: imine, imaxe

   ! number of columns per processors
   call ComputeEndpoints(MYRANKX, NRPROCX, n(1), p(1), nrcpp(1), ibeg(1), iend(1), mine(1), maxe(1))
   call ComputeEndpoints(MYRANKY, NRPROCY, n(2), p(2), nrcpp(2), ibeg(2), iend(2), mine(2), maxe(2))
   call ComputeEndpoints(MYRANKZ, NRPROCZ, n(3), p(3), nrcpp(3), ibeg(3), iend(3), mine(3), maxe(3))

   s(1) = iend(1) - ibeg(1) + 1
   s(2) = iend(2) - ibeg(2) + 1
   s(3) = iend(3) - ibeg(3) + 1

#ifdef IINFO
   write(*, *) PRINTRANK, 'Number of cols per processor:', nrcpp(1), nrcpp(2), nrcpp(3)
   write(*, *) PRINTRANK, 'ibegx,iendx', ibeg(1), iend(1)
   write(*, *) PRINTRANK, 'ibegy,iendy', ibeg(2), iend(2)
   write(*, *) PRINTRANK, 'ibegz,iendz', ibeg(3), iend(3)
#endif
   
      ! prepare dimensions vectors
   call FillDimVector(dimensionsX, shiftsX, nrcpp(1), s(2) * s(3), n(1), NRPROCX)
   call FillDimVector(dimensionsY, shiftsY, nrcpp(2), s(1) * s(3), n(2), NRPROCY)
   call FillDimVector(dimensionsZ, shiftsZ, nrcpp(3), s(1) * s(2), n(3), NRPROCZ)

   ! Compute indices for neighbours
   ibegsx = -1
   iendsx = -1
   ibegsy = -1
   iendsy = -1
   ibegsz = -1
   iendsz = -1

   do i = max(MYRANKX - 1, 0) + 1, min(MYRANKX + 1, NRPROCX - 1) + 1
      ix = i - MYRANKX + 1
      call ComputeEndpoints(i - 1, NRPROCX, n(1), p(1), nrcpp(1), ibegsx(ix), iendsx(ix), imine, imaxe)
   enddo
   do i = max(MYRANKY - 1, 0) + 1, min(MYRANKY + 1, NRPROCY - 1) + 1
      iy = i - MYRANKY + 1
      call ComputeEndpoints(i - 1, NRPROCY, n(2), p(2), nrcpp(2), ibegsy(iy), iendsy(iy), imine, imaxe)
   enddo
   do i = max(MYRANKZ - 1, 0) + 1, min(MYRANKZ + 1, NRPROCZ - 1) + 1
      iz = i - MYRANKZ + 1
      call ComputeEndpoints(i - 1, NRPROCZ, n(3), p(3), nrcpp(3), ibegsz(iz), iendsz(iz), imine, imaxe)
   enddo

end subroutine


! -------------------------------------------------------------------
! Allocates most of the 'static' arrays
! -------------------------------------------------------------------
subroutine AllocateArrays(&
   n, &
   s, &
   nrcpp, &
   Kl, &
   KU, &
   Mx, My, Mz, &
   F, F2, F3, &
   IPIVx, IPIVy, IPIVz, &
   R)
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ
   use mpi
   implicit none
   integer(kind = 4), intent(in), dimension(3) :: n
   integer(kind = 4), intent(in), dimension(3) :: s
   integer(kind = 4), intent(in), dimension(3) :: nrcpp
   integer(kind = 4), intent(in), dimension(3) :: KL, KU
   real (kind = 8), intent(out), allocatable, dimension(:,:) :: Mx
   real (kind = 8), intent(out), allocatable, dimension(:,:) :: My
   real (kind = 8), intent(out), allocatable, dimension(:,:) :: Mz
   real (kind = 8), intent(out), allocatable, dimension(:,:) :: F, F2, F3
   integer(kind = 4), intent(out), allocatable, dimension(:) :: IPIVx
   integer(kind = 4), intent(out), allocatable, dimension(:) :: IPIVy
   integer(kind = 4), intent(out), allocatable, dimension(:) :: IPIVz
   real (kind = 8), intent(out), allocatable :: R(:,:,:,:)
   integer :: ierr

   allocate(Mx(2 * KL(1) + KU(1) + 1, n(1) + 1))
   allocate(My(2 * KL(2) + KU(2) + 1, n(2) + 1))
   allocate(Mz(2 * KL(3) + KU(3) + 1, n(3) + 1))

   ! OLD: MP start with system fully generated along X
   ! allocate( F((n+1),(sy)*(sz))) !x,y,z
   allocate( F(s(1), s(2) * s(3))) !x,y,z
   allocate(F2(s(2), s(1) * s(3))) !y,x,z
   allocate(F3(s(3), s(1) * s(2))) !z,x,y


   ! Processes on the border need pivot vector for LAPACK call
   if (MYRANKX == 0 .or. MYRANKY == 0 .or. MYRANKZ == 0) then
      allocate(IPIVx(n(1) + 1))
      allocate(IPIVy(n(2) + 1))
      allocate(IPIVz(n(3) + 1))
   endif

   allocate(R(nrcpp(3) * nrcpp(1) * nrcpp(2), 3, 3, 3))
   R = 0.d0

   call mpi_barrier(MPI_COMM_WORLD, ierr)

end subroutine



!!!!! przeniesc do debug
! -------------------------------------------------------------------
! Prints debugging information about results of distributing
! data to neighbouring processes.
! -------------------------------------------------------------------
subroutine PrintDistributedData(ads)
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, PRINTRANK, &
   NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints
   implicit none
   type (ADS_setup) :: ads
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
            ads % R(:, i - MYRANKX + 1, j - MYRANKY + 1, k - MYRANKZ + 1), &
            (/ oendz - obegz + 1, oendx - obegx + 1, oendy - obegy + 1 /))
         enddo
      enddo
   enddo

   write(*, *) '----'

end subroutine


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
   type (ADS_setup) :: ads
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
   
   end subroutine



   !!!! podzielic na wraper i czesc wlasciwa
   ! przeniesc czesc do solver
   ! -------------------------------------------------------------------
   ! Performs one step of the simulation
   !
   ! iter - number of the iteration
   ! t    - time at the beginning of step
   ! -------------------------------------------------------------------
   subroutine Step(iter, RHS_fun, ads, mierr)
      use parallelism, ONLY:PRINTRANK, MYRANKX, MYRANKY, MYRANKZ
      use communicators, ONLY: COMMX, COMMY, COMMZ
      use reorderRHS, ONLY: ReorderRHSForX, ReorderRHSForY, ReorderRHSForZ
      use projection_engine, ONLY: Form3DRHS, ComputeMassMatrix
      use my_mpi, ONLY: DistributeSpline, Gather, Scatter
      use mpi
      implicit none
      interface
         subroutine RHS_fun(&
            X, &
            k, &
            e, &
            p, &
            nelem, &
            a, &
            b, &
            du, &
            ibeg, &
            iend, &
            mine, &
            maxe, &
            NNx, NNy, NNz, &
            Uval, J, W, F)
            implicit none
            integer(kind = 4), intent(in), dimension(3) :: p
            real (kind = 8), intent(in), dimension(3) :: X
            integer(kind = 4), intent(in), dimension(3) :: k
            integer(kind = 4), intent(in), dimension(3) :: e
            integer(kind = 4), intent(in), dimension(3) :: nelem
            integer(kind = 4), intent(in), dimension(3) :: ibeg
            integer(kind = 4), intent(in), dimension(3) :: iend
            integer(kind = 4), intent(in), dimension(3) :: maxe
            integer(kind = 4), intent(in), dimension(3) :: mine
            integer(kind = 4), intent(in), dimension(3) :: a
            integer(kind = 4), intent(in), dimension(3) :: b
            real (kind = 8), intent(in), dimension(3) :: du
            real (kind = 8), intent(in) :: Uval
            real (kind = 8), intent(in) :: J, W
            real (kind = 8), intent(in) :: NNx(0:p(1) - 1, 0:p(1), p(1) + 1, nelem(1)), &
            NNy(0:p(2) - 1, 0:p(2), p(2) + 1, nelem(2)), &
            NNz(0:p(3) - 1, 0:p(3), p(3) + 1, nelem(3))
            real (kind = 8), intent(out) :: F
         end subroutine
      end interface
      type (ADS_setup) :: ads
      integer(kind = 4), intent(out) :: mierr
      integer(kind = 4) :: iter
      integer(kind = 4) :: i
      integer(kind = 4) :: iret, ierr
      integer(kind = 4), dimension(3) :: nrcpp

      ! generate the RHS vectors
      call Form3DRHS(&
      ads % p, ads % n, ads % nelem, ads % nrcpp, &
      ads % ibeg, ads % iend, &
      ads % mine, ads % maxe, &
      ads % Ux, ads % Uy, ads % Uz, &
      ads % ibegsx, ads % iendsx, &
      ads % ibegsy, ads % iendsy, &
      ads % ibegsz, ads % iendsz, &
      ads % R, ads % F, RHS_fun)

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
      
      allocate(ads % F_out((ads % n(1) + 1), ads % s(2) * ads % s(3)))

      call Gather(ads % F, ads % F_out, ads % n(1), ads % s(1), ads % s(2) &
      *ads % s(3), ads % dimensionsX, ads % shiftsX, COMMX, ierr)

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
         
         call ComputeMassMatrix(ads % KL(1), ads % KU(1), ads % Ux, ads % p(1), &
         ads % n(1), ads % nelem(1), ads % Mx)
         call SolveOneDirection(ads % F_out, ads % s(2) * ads % s(3), ads % n(1), &
         ads % KL(1), ads % KU(1), ads % p(1), ads % Mx, ads % IPIVx)
      endif

      call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
      write(*, *) PRINTRANK, '1c) SCATTER'
#endif
      allocate(ads % F2_out(ads % s(1), ads % s(2) * ads % s(3)))
      call Scatter(ads % F_out, ads % F2_out, ads % n(1), ads % s(1), ads % s(2) * &
      ads % s(3), ads % dimensionsX, ads % shiftsX, COMMX, ierr)
      deallocate(ads % F_out)

      call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
      write(*, *) PRINTRANK, '1d) REORDER'
#endif
      call ReorderRHSForY(ads % ibeg, ads % iend, ads % F2_out, ads % F2)
      deallocate(ads % F2_out)

#ifdef IPRINT
      write(*, *) PRINTRANK, 'after ReorderRHSForY'
      write(*, *) PRINTRANK, 'F2:'
      do i = 1, ads % s(2)
         write(*, *) PRINTRANK, i, 'row=', ads % F2(i, 1:ads % s(1) * ads % s(3))
      enddo
#endif
      
      !--------------------------------------------------------------------
      ! Solve the second problem
      !--------------------------------------------------------------------
      call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
      write(*, *) PRINTRANK, '2a) GATHER'
#endif
      
      allocate(ads % F2_out((ads % n(2) + 1), ads % s(1) * ads % s(3)))
      call Gather(ads % F2, ads % F2_out, ads % n(2), ads % s(2), ads % s(1) * ads % s(3), &
      ads % dimensionsY, ads % shiftsY, COMMY, ierr)

      call mpi_barrier(MPI_COMM_WORLD, ierr)

      if (MYRANKY == 0) then
#ifdef IINFO
         write(*, *) PRINTRANK, '2b) SOLVE THE SECOND PROBLEM'
#endif
         
         call ComputeMassMatrix(ads % KL(2), ads % KU(2), ads % Uy, ads % p(2), ads % n(2), &
         ads % nelem(2), ads % My)
         call SolveOneDirection(ads % F2_out, ads % s(1) * ads % s(3), ads % n(2), ads % KL(2), &
         ads % KU(2), ads % p(2), ads % My, ads % IPIVy)
      endif

      call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
      write(*, *) PRINTRANK, '2c) SCATHER'
#endif
      
      ! CORRECTION
      allocate(ads % F3_out(ads % s(3), ads % s(1) * ads % s(3)))
      call Scatter(ads % F2_out, ads % F3_out, ads % n(2), ads % s(2), ads % s(1) * ads % s(3), &
      ads % dimensionsY, ads % shiftsY, COMMY, ierr)
      deallocate(ads % F2_out)

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
      call ReorderRHSForZ(ads % ibeg, ads % iend, ads % F3_out, ads % F3)
      deallocate(ads % F3_out)

#ifdef IPRINT
      write(*, *) PRINTRANK, 'after ReorderRHSForZ'
      write(*, *) PRINTRANK, 'F3:'
      do i = 1, ads % s(3)
         write(*, *) PRINTRANK, i, 'row=', ads % F3(i, 1:ads % s(1) * ads % s(2))
      enddo
#endif
      
      !--------------------------------------------------------------------
      ! Solve the third problem
      !--------------------------------------------------------------------
#ifdef IINFO
      write(*, *) PRINTRANK, '3a) GATHER'
#endif
      
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      allocate(ads % F3_out((ads % n(3) + 1), ads % s(1) * ads % s(2)))
      call Gather(ads % F3, ads % F3_out, ads % n(3), ads % s(3), ads % s(1) * ads % s(2), &
      ads % dimensionsZ, ads % shiftsZ, COMMZ, ierr)

      call mpi_barrier(MPI_COMM_WORLD, ierr)

      if (MYRANKZ == 0) then
#ifdef IINFO
         write(*, *) PRINTRANK, '3b) SOLVE THE THIRD PROBLEM'
#endif
         
         call ComputeMassMatrix(ads % KL(3), ads % KU(3), ads % Uz, ads % p(3), ads % n(3), &
         ads % nelem(3), ads % Mz)
         call SolveOneDirection(ads % F3_out, ads % s(1) * ads % s(2), ads % n(3), ads % KL(3), &
         ads % KU(3), ads % p(3), ads % Mz, ads % IPIVz)
      endif

      call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
      write(*, *) PRINTRANK, '3c) SCATTER'
#endif
      
      ! CORRECTION
      allocate(ads % F_out(ads % s(3), ads % s(1) * ads % s(2)))
      call Scatter(ads % F3_out, ads % F_out, ads % n(3), ads % s(3), ads % s(1) * ads % s(2), &
      ads % dimensionsZ, ads % shiftsZ, COMMZ, ierr)
      deallocate(ads % F3_out)

      call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef IINFO
      write(*, *) PRINTRANK, '3d) REORDER'
#endif
      ! Reorder right hand sides
      call ReorderRHSForX(ads % ibeg, ads % iend, ads % F_out, ads % F)
      deallocate(ads % F_out)

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
      ads % R(1:ads % s(1) * ads % s(2) * ads % s(3), 2, 2, 2) = reshape(ads % F, &
      (/ ads % s(1) * ads % s(2) * ads % s(3) /) )
      nrcpp = (/ ads % nrcpp(3), ads % nrcpp(1), ads % nrcpp(2) /)
      call DistributeSpline(ads % R, nrcpp, ads % R)

#ifdef IPRINT
      write(*, *) PRINTRANK, 'Result:'
      do i = 1, ads % s(3)
         write(*, *) PRINTRANK, i, 'row=', ads % F(i,:)
      enddo
#endif
      

      call mpi_barrier(MPI_COMM_WORLD, ierr)

      mierr = 0
   end subroutine




   ! -------------------------------------------------------------------
   ! Deallocates all the resources and finalizes MPI.
   ! -------------------------------------------------------------------
   subroutine Cleanup(ads, mierr)
      use parallelism, ONLY: PRINTRANK
      implicit none
      type (ADS_setup) :: ads
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

      if (allocated(ads % Mx)) deallocate(ads % Mx)
      if (allocated(ads % My)) deallocate(ads % My)
      if (allocated(ads % Mz)) deallocate(ads % Mz)

      if (allocated(ads % F)) deallocate(ads % F)
      if (allocated(ads % F2)) deallocate(ads % F2)
      if (allocated(ads % F3)) deallocate(ads % F3)

      !!!!!! wyciac
      call mpi_finalize(ierr)
#ifdef IINFO
      write(*, *) PRINTRANK, "Exiting..."
#endif

      mierr = 0

   end subroutine



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

   end subroutine



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

   end subroutine


   ! -------------------------------------------------------------------
   ! Gathers full solution and plots it
   ! -------------------------------------------------------------------
   subroutine PrintSolution(iter, t, ads)
      use parallelism, ONLY: MYRANK
      use plot, ONLY: SaveSplinePlot, PlotParams
      use vtk, ONLY: VtkOutput
      use my_mpi, ONLY: GatherFullSolution
      implicit none
      type (ADS_setup), intent(in) :: ads
      integer(kind = 4), intent(in) :: iter
      real (kind = 8), intent(in) :: t
      real (kind = 8), allocatable :: solution(:,:,:)
      type (PlotParams) :: params
      character(len = 20) :: filename

      call GatherFullSolution(0, ads % F, solution, &
      [ads % n(1), ads % n(2), ads % n(3)], [ads % p(1), ads % p(2), ads % p(3)], [ads % s(1), ads % s(2), ads % s(3)])

      if (MYRANK == 0) then
         write(filename, '(I10)') iter
         filename = 'step' // adjustl(filename)
         ! filename = trim(filename) // '_'

         params = PlotParams(0, 1, 0, 1, 0, 1, 31, 31, 31)
         call SaveSplinePlot(trim(filename), &
         ads % Ux, ads % p(1), ads % n(1), ads % nelem(1), &
         ads % Uy, ads % p(2), ads % n(2), ads % nelem(2), &
         ads % Uz, ads % p(3), ads % n(3), ads % nelem(3), &
         ! solution, GnuPlotOutput, params)
         solution, VtkOutput, params)

         ! call SavePlot(trim(filename), ftest, GnuPlotOutput, params)
      endif

   end subroutine



end module
