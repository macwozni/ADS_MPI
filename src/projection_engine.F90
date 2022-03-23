!------------------------------------------------------------------------------
! AGH University of Science and Technology, Krakow, Poland
!------------------------------------------------------------------------------
!
! MODULE: projection_engine
!
!> @author
!> Maciej Wozniak
!
! DESCRIPTION:
!> This module contains all functionality associated to projection.
!
! REVISION HISTORY:
! 21 11 2017 - Initial Version
!
!------------------------------------------------------------------------------

module projection_engine

   implicit none

contains

!---------------------------------------------------------------------------
!> @author Maciej Wozniak
!>
!> @brief
!> Calculates matrices M, K, B and BT.
!>
!> Calculates:
!>
!>  - the mass matrix M
!>  - the stifness matrix K
!>  - the advection matrix B
!>  - the advection matrix transposed BT.
!
! Input:
! ------
!> @param[in] U1      - knot vector
!> @param[in] p1      - degree of approximation
!> @param[in] n1      - number of control points minus one
!> @param[in] nelem1  - number of subintervals in knot
!> @param[in] U2      - knot vector
!> @param[in] p2      - degree of approximation
!> @param[in] n2      - number of control points minus one
!> @param[in] nelem2  - number of subintervals in knot
!> @param[in] mix     - mixing proportions of M, K, B and BT matrices
!
! Output:
! -------
!> @param[out] sprsmtrx - sparse matrix, logically \f$ (n+1) \times (n+1) \f$, combination of M, K, B and BT matrices
!>
!> Values in the matrix are stored in the band format, i.e. while \f$ M \f$
!> is \f$ (n+1) \times (n+1) \f$, it is stored as \f$ (2 KL + KU + 1) \times n \f$, and the
!> index correspondence is given by:
!>
!>     A(i, j) = M(KL + KU + 1 + i - j, j)
!>
!> \f$ M = u*v \f$
! -------------------------------------------------------------------
   subroutine MKBBT_large(nelem, U1, p1, n1, U2, p2, n2, mixA, mixB, mixBT, sprsmtrx)
      use basis, ONLY: BasisData
      use omp_lib
      use sparse
      implicit none
      integer(kind=4), intent(in) :: nelem
      integer(kind=4), intent(in) :: n1, p1
      integer(kind=4), intent(in) :: n2, p2
      real(kind=8), intent(in) :: U1(0:n1 + p1 + 1)
      real(kind=8), intent(in) :: U2(0:n2 + p2 + 1)
      real(kind=8), dimension(4), intent(in) :: mixA, mixB, mixBT
      real(kind=8), dimension(nelem) :: J ! values of the Jacobian of elements
      real(kind=8), dimension(p1 + 1) :: W ! weights of Gauss quadrature points
      real(kind=8), dimension(p1 + 1, nelem) :: X ! points of Gauss quadrature
      real(kind=8), dimension(0:1, 0:p1, p1 + 1, nelem) :: NN1 ! values of (p1+1) nonzero basis functions and their derivatives at points of Gauss quadrature
      real(kind=8), dimension(0:1, 0:p2, p1 + 1, nelem) :: NN2 ! values of (p2+1) nonzero basis functions and their derivatives at points of Gauss quadrature
      integer(kind=4) :: dd ! order of highest derivatives we want to compute
      integer(kind=4) :: ia, ib
      integer(kind=4) :: mm1, mm2
      integer(kind=4) :: ng ! number of Gauss quadrature points
      integer(kind=4) :: e, i, c, d
      integer(kind=4) :: O1(nelem) ! indexes of first nonzero functions on each element
      integer(kind=4) :: O2(nelem) ! indexes of first nonzero functions on each element
      type(sparse_matrix), pointer, intent(out) :: sprsmtrx
      real(kind=8) :: val
      real(kind=8) :: M, K, B, BT

      mm1 = n1 + p1 + 1
      ng = p1 + 1
      dd = 1
      mm2 = n2 + p2 + 1

! test
      call BasisData(p1, mm1, U1, dd, ng, nelem, O1, J, W, X, NN1)
! trial
      call BasisData(p2, mm2, U2, dd, ng, nelem, O2, J, W, X, NN2)

      call initialize_sparse(n1 + n2 + 2, n1 + n2 + 2, sprsmtrx)

      ! total_size = (nelem1)*(ng1)*(p1 + 1)*(p1 + 1)
! submatrix A
! new parallel loop
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(PRIVATE) &
! !$OMP PRIVATE(d,c,i,e,ia,ib,tmp) &
! !$OMP SHARED(nelem,ng,p,O,NN,W,J,total_size) &
! !$OMP REDUCTION(+:M) &
! !$OMP REDUCTION(+:K) &
! !$OMP REDUCTION(+:B) &
! !$OMP REDUCTION(+:BT)
! loop over elements
      do e = 1, nelem
! loop over Gauss points
         do i = 1, ng
! loop over shape functions over elements (p+1 functions)
            do c = 0, p1
               ! loop over shape functions over elements (p+1 functions)
               do d = 0, p1
                  ! O(e) + c = first dof of element + 1st local shape function index
                  ! O(e) + d = first dof of element + 2nd local shape function index
                  ! NN(0,c,i,e) = value of shape function c at Gauss point i over element e
                  ! NN(0,d,i,e) = value of shape function d at Gauss point i over element e
                  ! W(i) weight for Gauss point i
                  ! J(e) jacobian for element e
                  ia = O1(e) + c
                  ib = O1(e) + d
                  ! M = u*v
                  M = NN1(0, c, i, e)*NN1(0, d, i, e)*J(e)*W(i)
                  K = NN1(1, c, i, e)*NN1(1, d, i, e)*J(e)*W(i)
                  B = NN1(1, c, i, e)*NN1(0, d, i, e)*J(e)*W(i)
                  BT = NN1(0, c, i, e)*NN1(1, d, i, e)*J(e)*W(i)
                  val = mixA(1)*M + mixA(2)*K + mixA(3)*B + mixA(4)*BT
                  call add(sprsmtrx, ia, ib, val)
               end do
            end do
         end do
      end do
! !$OMP END PARALLEL DO

      ! total_size = (nelem1)*(ng1)*(p1 + 1)*(p2 + 1)
! submatrix B
! new parallel loop
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(PRIVATE) &
! !$OMP PRIVATE(d,c,i,e,ia,ib,tmp) &
! !$OMP SHARED(nelem,ng,p,O,NN,W,J,total_size) &
! !$OMP REDUCTION(+:M) &
! !$OMP REDUCTION(+:K) &
! !$OMP REDUCTION(+:B) &
! !$OMP REDUCTION(+:BT)
! loop over elements
      do e = 1, nelem
! loop over Gauss points
         do i = 1, ng
! loop over shape functions over elements (p+1 functions)
            do c = 0, p2
               ! loop over shape functions over elements (p+1 functions)
               do d = 0, p1
                  ! O(e) + c = first dof of element + 1st local shape function index
                  ! O(e) + d = first dof of element + 2nd local shape function index
                  ! NN(0,c,i,e) = value of shape function c at Gauss point i over element e
                  ! NN(0,d,i,e) = value of shape function d at Gauss point i over element e
                  ! W(i) weight for Gauss point i
                  ! J(e) jacobian for element e
                  ia = O2(e) + c
                  ib = O1(e) + d + n1 + 1
                  ! M = u*v
                  M = NN2(0, c, i, e)*NN1(0, d, i, e)*J(e)*W(i)
                  K = NN2(1, c, i, e)*NN1(1, d, i, e)*J(e)*W(i)
                  B = NN2(1, c, i, e)*NN1(0, d, i, e)*J(e)*W(i)
                  BT = NN2(0, c, i, e)*NN1(1, d, i, e)*J(e)*W(i)
                  val = mixB(1)*M + mixB(2)*K + mixB(3)*B + mixB(4)*BT
                  call add(sprsmtrx, ia, ib, val)
               end do
            end do
         end do
      end do
! !$OMP END PARALLEL DO

      ! total_size = (nelem2)*(ng2)*(p2 + 1)*(p1 + 1)
! submatrix BT
! new parallel loop
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(PRIVATE) &
! !$OMP PRIVATE(d,c,i,e,ia,ib,tmp) &
! !$OMP SHARED(nelem,ng,p,O,NN,W,J,total_size) &
! !$OMP REDUCTION(+:M) &
! !$OMP REDUCTION(+:K) &
! !$OMP REDUCTION(+:B) &
! !$OMP REDUCTION(+:BT)
      do e = 1, nelem
! loop over Gauss points
         do i = 1, ng
! loop over shape functions over elements (p+1 functions)
            do c = 0, p1
! loop over shape functions over elements (p+1 functions)
               do d = 0, p2
                  ! O(e) + c = first dof of element + 1st local shape function index
                  ! O(e) + d = first dof of element + 2nd local shape function index
                  ! NN(0,c,i,e) = value of shape function c at Gauss(i) weight for Gauss point i
                  ! J(e) jacobian for element e
                  ia = O1(e) + c + n2 + 1
                  ib = O2(e) + d
                  ! M = u*v
                  M = NN1(0, c, i, e)*NN2(0, d, i, e)*J(e)*W(i)
                  K = NN1(1, c, i, e)*NN2(1, d, i, e)*J(e)*W(i)
                  B = NN1(1, c, i, e)*NN2(0, d, i, e)*J(e)*W(i)
                  BT = NN1(0, c, i, e)*NN2(1, d, i, e)*J(e)*W(i)
                  val = mixBT(1)*M + mixBT(2)*K + mixBT(3)*B + mixBT(4)*BT
                  call add(sprsmtrx, ia, ib, val)
               end do
            end do
         end do
      end do
! !$OMP END PARALLEL DO

   end subroutine MKBBT_large

!---------------------------------------------------------------------------
!> @author Maciej Wozniak
!>
!> @brief
!> Calculates matrices M, K, B and BT.
!>
!> Calculates:
!>
!>  - the mass matrix M
!>  - the stifness matrix K
!>  - the advection matrix B
!>  - the advection matrix transposed BT.
!
! Input:
! ------
!> @param[in] U1      - knot vector
!> @param[in] p1      - degree of approximation
!> @param[in] n1      - number of control points minus one
!> @param[in] nelem1  - number of subintervals in knot
!> @param[in] U2      - knot vector
!> @param[in] p2      - degree of approximation
!> @param[in] n2      - number of control points minus one
!> @param[in] nelem2  - number of subintervals in knot
!> @param[in] mix     - mixing proportions of M, K, B and BT matrices
!
! Output:
! -------
!> @param[out] sprsmtrx - sparse matrix, logically \f$ (n+1) \times (n+1) \f$, combination of M, K, B and BT matrices
!>
!> Values in the matrix are stored in the band format, i.e. while \f$ M \f$
!> is \f$ (n+1) \times (n+1) \f$, it is stored as \f$ (2 KL + KU + 1) \times n \f$, and the
!> index correspondence is given by:
!>
!>     A(i, j) = M(KL + KU + 1 + i - j, j)
!>
!> \f$ M = u*v \f$
! -------------------------------------------------------------------
   subroutine MKBBT_small(nelem, U, p, n, mix, sprsmtrx)
      use basis, ONLY: BasisData
      use omp_lib
      use sparse
      implicit none
      integer(kind=4), intent(in) :: n, p, nelem
      real(kind=8), intent(in) :: U(0:n + p + 1)
      real(kind=8), dimension(4), intent(in) :: mix
      real(kind=8), dimension(nelem) :: J ! values of the Jacobian of elements
      real(kind=8), dimension(p + 1) :: W ! weights of Gauss quadrature points
      real(kind=8), dimension(p + 1, nelem) :: X ! points of Gauss quadrature
      real(kind=8), dimension(0:1, 0:p, p + 1, nelem) :: NN ! values of (p1+1) nonzero basis functions and their derivatives at points of Gauss quadrature
      integer(kind=4) :: dd ! order of highest derivatives we want to compute
      integer(kind=4) :: ia, ib
      integer(kind=4) :: mm
      integer(kind=4) :: ng ! number of Gauss quadrature points
      integer(kind=4) :: e, i, c, d
      integer(kind=4) :: O(nelem) ! indexes of first nonzero functions on each element
      type(sparse_matrix), pointer, intent(out) :: sprsmtrx
      real(kind=8) :: val
      real(kind=8) :: M, K, B, BT

      mm = n + p + 1
      ng = p + 1
      dd = 1

      call BasisData(p, mm, U, dd, ng, nelem, O, J, W, X, NN)

      call initialize_sparse(n + 1, n + 1, sprsmtrx)

! submatrix A
! new parallel loop
!!$OMP PARALLEL DO &
!!$OMP DEFAULT(PRIVATE) &
!!$OMP PRIVATE(d,c,i,e,ia,ib,tmp) &
!!$OMP SHARED(nelem,ng,p,O,NN,W,J,total_size) &
!!$OMP REDUCTION(+:M) &
!!$OMP REDUCTION(+:K) &
!!$OMP REDUCTION(+:B) &
!!$OMP REDUCTION(+:BT)
! loop over elements
      do e = 1, nelem
! loop over Gauss points
         do i = 1, ng
            ! loop over shape functions over elements (p+1 functions)
            do c = 0, p
               ! loop over shape functions over elements (p+1 functions)
               do d = 0, p
                  ! O(e) + c = first dof of element + 1st local shape function index
                  ! O(e) + d = first dof of element + 2nd local shape function index
                  ! NN(0,c,i,e) = value of shape function c at Gauss point i over element e
                  ! NN(0,d,i,e) = value of shape function d at Gauss point i over element e
                  ! W(i) weight for Gauss point i
                  ! J(e) jacobian for element e
                  ia = O(e) + c
                  ib = O(e) + d
                  ! M = u*v
                  M = NN(0, c, i, e)*NN(0, d, i, e)*J(e)*W(i)
                  K = NN(1, c, i, e)*NN(1, d, i, e)*J(e)*W(i)
                  B = NN(1, c, i, e)*NN(0, d, i, e)*J(e)*W(i)
                  BT = NN(0, c, i, e)*NN(1, d, i, e)*J(e)*W(i)
                  val = mix(1)*M + mix(2)*K + mix(3)*B + mix(4)*BT
                  call add(sprsmtrx, ia, ib, val)
               end do
            end do
         end do
      end do
!!$OMP END PARALLEL DO

   end subroutine MKBBT_small

!---------------------------------------------------------------------------
!> @author Maciej Wozniak
!>
!> @brief
!> Calculate right-hand side of the equation.
!
! Input:
! ------
!> @param[in] ads          - ADS setup structure
!> @param[in] directon     - direction for the substep
!> @param[in] substep      - number of substep
! n               - nuber of previous time steps
!
! Input/Output:
! -------
!> @param[inout] ads_data  - data structures for ADS
!
! Output:
! -------
!> @param[out] l2norm      -
! -------------------------------------------------------------------
   subroutine Form3DRHS(ads_test, ads_trial, ads_data, direction, n, substep, alpha_step, forcing, igrm, l2norm)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      ! use parallelism, ONLY: PRINTRANK
      use Interfaces, ONLY: forcing_fun
      USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
      use omp_lib
      use RHS_eq
      implicit none
      procedure(forcing_fun) :: forcing
      type(ADS_setup), intent(in) :: ads_test, ads_trial
      integer(kind=4), dimension(3), intent(in) :: direction
      integer(kind=4), intent(in) :: substep
      type(ADS_compute_data), intent(inout) :: ads_data
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in), dimension(7, 3) :: alpha_step
      real(kind=8), intent(out) :: l2norm
      integer(kind=4) :: kx, ky, kz, ax, ay, az, ex, ey, ez!, exx,eyy,ezz
      real(kind=8) :: J, W
      integer(kind=4) :: ind, ind1, ind23, indx, indy, indz
      real(kind=8) :: resvalue
      real(kind=8), dimension(3) :: X
      integer(kind=4), dimension(3) :: k, e, a
      ! integer (kind = 4) :: tmp, all
      ! integer (kind = 4) :: total_size
      real(kind=8), dimension(3)  :: du
      ! integer(kind = 4) :: indbx, indby, indbz
      real(kind=8) :: Uval
      real(kind=8) :: Uval13
      real(kind=8) :: Uval23
      real(kind=8), dimension(:, :, :), allocatable :: elarr
      real(kind=8) :: l2normtmp
      type(ADS_setup) :: ads
      logical, intent(out) :: igrm

!  copy default space as trial space
      ads = ads_trial
      igrm = .FALSE.

!  if we have enriched one direction, then modify default space
      if (direction(1) .EQ. 1) then
         ads%n(1) = ads_test%n(1)
         ads%p(1) = ads_test%p(1)
         ads%Ux = ads_test%Ux
         ads%nelem(1) = ads_test%nelem(1)
         ads%dimensionsX = ads_test%dimensionsX
         ads%shiftsX = ads_test%shiftsX
         ads%IPIVx = ads_test%IPIVx
         ads%nrcpp(1) = ads_test%nrcpp(1)
         ads%ibeg(1) = ads_test%ibeg(1)
         ads%iend(1) = ads_test%iend(1)
         ads%s(1) = ads_test%s(1)
         ads%ibegsx(1) = ads_test%ibegsx(1)
         ads%iendsx(1) = ads_test%iendsx(1)
         ads%ibegsy(1) = ads_test%ibegsy(1)
         ads%iendsy(1) = ads_test%iendsy(1)
         ads%ibegsz(1) = ads_test%ibegsz(1)
         ads%iendsz(1) = ads_test%iendsz(1)
         ads%mine(1) = ads_test%mine(1)
         ads%maxe(1) = ads_test%maxe(1)
         ads%lnelem(1) = ads_test%lnelem(1)
         ads%m(1) = ads_test%m(1)
         ads%ng(1) = ads_test%ng(1)
         ads%Ox = ads_test%Ox
         ads%Jx = ads_test%Jx
         ads%Xx = ads_test%Xx
         ads%NNx = ads_test%NNx
         ads%Wx = ads_test%Wx
         igrm = .TRUE.
      end if
      if (direction(2) .EQ. 1) then
         ads%n(2) = ads_test%n(2)
         ads%p(2) = ads_test%p(2)
         ads%Uy = ads_test%Uy
         ads%nelem(2) = ads_test%nelem(2)
         ads%dimensionsY = ads_test%dimensionsY
         ads%shiftsY = ads_test%shiftsY
         ads%IPIVy = ads_test%IPIVy
         ads%nrcpp(2) = ads_test%nrcpp(2)
         ads%ibeg(2) = ads_test%ibeg(2)
         ads%iend(2) = ads_test%iend(2)
         ads%s(2) = ads_test%s(2)
         ads%ibegsx(2) = ads_test%ibegsx(2)
         ads%iendsx(2) = ads_test%iendsx(2)
         ads%ibegsy(2) = ads_test%ibegsy(2)
         ads%iendsy(2) = ads_test%iendsy(2)
         ads%ibegsz(2) = ads_test%ibegsz(2)
         ads%iendsz(2) = ads_test%iendsz(2)
         ads%mine(2) = ads_test%mine(2)
         ads%maxe(2) = ads_test%maxe(2)
         ads%lnelem(2) = ads_test%lnelem(2)
         ads%m(2) = ads_test%m(2)
         ads%ng(2) = ads_test%ng(2)
         ads%Oy = ads_test%Oy
         ads%Jy = ads_test%Jy
         ads%Xy = ads_test%Xy
         ads%NNy = ads_test%NNy
         ads%Wy = ads_test%Wy
         igrm = .TRUE.
      end if
      if (direction(3) .EQ. 1) then
         ads%n(3) = ads_test%n(3)
         ads%p(3) = ads_test%p(3)
         ads%Uz = ads_test%Uz
         ads%nelem(3) = ads_test%nelem(3)
         ads%dimensionsZ = ads_test%dimensionsZ
         ads%shiftsZ = ads_test%shiftsZ
         ads%IPIVz = ads_test%IPIVz
         ads%nrcpp(3) = ads_test%nrcpp(3)
         ads%ibeg(3) = ads_test%ibeg(3)
         ads%iend(3) = ads_test%iend(3)
         ads%s(3) = ads_test%s(3)
         ads%ibegsx(3) = ads_test%ibegsx(3)
         ads%iendsx(3) = ads_test%iendsx(3)
         ads%ibegsy(3) = ads_test%ibegsy(3)
         ads%iendsy(3) = ads_test%iendsy(3)
         ads%ibegsz(3) = ads_test%ibegsz(3)
         ads%iendsz(3) = ads_test%iendsz(3)
         ads%mine(3) = ads_test%mine(3)
         ads%maxe(3) = ads_test%maxe(3)
         ads%lnelem(3) = ads_test%lnelem(3)
         ads%m(3) = ads_test%m(3)
         ads%ng(3) = ads_test%ng(3)
         ads%Oz = ads_test%Oz
         ads%Jz = ads_test%Jz
         ads%Xz = ads_test%Xz
         ads%NNz = ads_test%NNz
         ads%Wz = ads_test%Wz
         igrm = .TRUE.
      end if

      allocate (elarr(0:ads%p(1), 0:ads%p(2), 0:ads%p(3)))
      ! total_size = ads % lnelem(1) * ads % lnelem(2) * ads % lnelem(3)

      l2norm = 0.d0
!   if (allocated(ads_data%F)) ads_data%F = 0.d0
!   if (allocated(ads_data%Ft)) ads_data%Ft = 0.d0

!      loop over points
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(SHARED) &
! !$OMP SHARED(ads,ads_data,total_size) &
! !$OMP PRIVATE(tmp,ex,ey,ez,e,kx,ky,kz,k,W,ax,ay,az,a,ind,indx,indy,indz,ind1,ind23,J) &
! !$OMP PRIVATE(X,du,resvalue) &
! !$OMP PRIVATE(indbx,indby,indbz,Uval,elarr,l2normtmp,Uval_m,Uval13,Uval23) &
! !$OMP REDUCTION(+:l2norm)
      ! do all = 1, total_size
! translate coefficients to local
      ! ez = modulo(all - 1, ads % lnelem(3))
      ! tmp = (all - ez)/ads % lnelem(3) + 1
      ! ey = modulo(tmp - 1, ads % lnelem(2))
      ! ex = (tmp - ey)/ads % lnelem(2)
      ! write(*,*) size(ads%Jx) , ads % lnelem(1), ads % mine(1)
      ! write(*,*) size(ads%Jy) , ads % lnelem(2), ads % mine(2)
      ! write(*,*) size(ads%Jz) , ads % lnelem(3), ads % mine(3)
      ! do exx=1,ads % lnelem(1)
      ! do eyy=1,ads % lnelem(2)
      ! do ezz=1,ads % lnelem(3)
      do ex = ads%mine(1), ads%maxe(1)
         do ey = ads%mine(2), ads%maxe(2)
            do ez = ads%mine(3), ads%maxe(3)
! fix distributed part
               ! ex = exx + ads % mine(1)
               ! ey = eyy + ads % mine(2)
               ! ez = ezz + ads % mine(3)
! Jacobian
               J = ads%Jx(ex)*ads%Jy(ey)*ads%Jz(ez)
               e = (/ex, ey, ez/)
               elarr = 0.d0
! loop over quadrature points
               do kx = 1, ads%ng(1)
                  do ky = 1, ads%ng(2)
                     do kz = 1, ads%ng(3)
                        k = (/kx, ky, kz/)
! weigths
                        W = ads%Wx(kx)*ads%Wy(ky)*ads%Wz(kz)
                        Uval = ads_data%Un(ex, ey, ez, kx, ky, kz)
                        Uval13 = ads_data%Un13(ex, ey, ez, kx, ky, kz)
                        Uval23 = ads_data%Un23(ex, ey, ez, kx, ky, kz)
                        du = ads_data%dUn(ex, ey, ez, kx, ky, kz, :)

!                 loop over degrees of freedom
                        do ax = 0, ads%p(1)
                           do ay = 0, ads%p(2)
                              do az = 0, ads%p(3)
                                 indx = (ads%Ox(ex) + ax)
                                 indy = (ads%Oy(ey) + ay)
                                 indz = (ads%Oz(ez) + az)
                                 ind = indx + (indy + indz*(ads%n(2) + 1))*(ads%n(1) + 1)

                                 if ((indx < ads%ibeg(1) - 1) .or. (indx > ads%iend(1) - 1) .or. &
                                     (indy < ads%ibeg(2) - 1) .or. (indy > ads%iend(2) - 1) .or. &
                                     (indz < ads%ibeg(3) - 1) .or. (indz > ads%iend(3) - 1)) then
                                 else
                                    ind1 = indx - ads%ibeg(1) + 1
                                    ind23 = (indy - ads%ibeg(2) + 1) + &
                                            (indz - ads%ibeg(3) + 1)*(ads%iend(2) - ads%ibeg(2) + 1)
                                    X = (/ads%Xx(kx, ex), ads%Xy(ky, ey), ads%Xz(kz, ez)/)
                                    a = (/ax, ay, az/)

                                    ! call RHS_fun(&
                                    ! ads, &
                                    ! X, &
                                    ! k, &
                                    ! e, &
                                    ! a, &
                                    ! du, &
                                    ! 1, Uval_m, Uval13,Uval23, &
                                    ! ads_data, J, W, direction, substep, l2normtmp, resvalue)

                                    call ComputePointForRHS( &
                                       ads, &
                                       X, &
                                       k, &
                                       e, &
                                       a, &
                                       du, &
                                       n, &
                                       Uval, &
                                       Uval13, &
                                       Uval23, &
                                       ads_data, J, W, direction, substep, &
                                       alpha_step, &
                                       forcing, &
                                       l2normtmp, resvalue)

                                    elarr(ax, ay, az) = elarr(ax, ay, az) + resvalue
                                    l2norm = l2norm + l2normtmp
                                 end if
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
! moving results from temporary array to main one
! !$OMP CRITICAL
               do ax = 0, ads%p(1)
                  do ay = 0, ads%p(2)
                     do az = 0, ads%p(3)
                        indx = (ads%Ox(ex) + ax)
                        indy = (ads%Oy(ey) + ay)
                        indz = (ads%Oz(ez) + az)
                        ind1 = indx - ads%ibeg(1) + 1
                        ind23 = (indy - ads%ibeg(2) + 1) + &
                                (indz - ads%ibeg(3) + 1)*(ads%iend(2) - ads%ibeg(2) + 1)
                        if ((indx < ads%ibeg(1) - 1) .or. (indx > ads%iend(1) - 1) .or. &
                            (indy < ads%ibeg(2) - 1) .or. (indy > ads%iend(2) - 1) .or. &
                            (indz < ads%ibeg(3) - 1) .or. (indz > ads%iend(3) - 1)) then
                        else
                           if (igrm) then
                              ads_data%Ft(ind1 + 1, ind23 + 1) = &
                                 ads_data%Ft(ind1 + 1, ind23 + 1) &
                                 + elarr(ax, ay, az)
                           else
                              ads_data%F(ind1 + 1, ind23 + 1) = &
                                 ads_data%F(ind1 + 1, ind23 + 1) &
                                 + elarr(ax, ay, az)
                           end if
                        end if
                     end do
                  end do
               end do
! !$OMP END CRITICAL
            end do
         end do
      end do
! !$OMP END PARALLEL DO

      deallocate (elarr)

   end subroutine Form3DRHS

!---------------------------------------------------------------------------
!> @author Maciej Wozniak
!>
!> @brief
!> Calculates value of derivative from previous time step - dUn.
!> Calculates previous solution coefficient - Un.
!
! Input:
! ------
!> @param[in] ads      - ADS setup structure
!> @param[in] ads_data - data structures for ADS
!
! Output:
! -------
!> @param[out] Un      -
!> @param[out] dUn     -
! -------------------------------------------------------------------
   subroutine FormUn(subun, ads, ads_data)
      USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      ! use parallelism, ONLY: PRINTRANK
      use Interfaces, ONLY: forcing_fun
      USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
      use omp_lib
      implicit none
      integer(kind=4), intent(in) :: subun
      type(ADS_setup), intent(in) :: ads
      type(ADS_compute_data), intent(inout) :: ads_data
      integer(kind=4) :: kx, ky, kz, ex, ey, ez
      integer(kind=4) :: ind
      integer(kind=4) :: tmp, all
      integer(kind=4) :: total_size
      integer(kind=4) :: rx, ry, rz, ix, iy, iz, sx, sy, sz
      integer(kind=4) :: bx, by, bz
      real(kind=8) :: dux, duy, duz
      ! real   (kind=8), dimension(3)  :: du
      integer(kind=4) :: indbx, indby, indbz
      real(kind=8) :: Uval, ucoeff
      real(kind=8) :: dvx, dvy, dvz, v

      ads_data%un = 0.d0
      ads_data%un13 = 0.d0
      ads_data%un23 = 0.d0

      if (subun .EQ. 1) then
         ads_data%Un = 0.d0
      else if (subun .EQ. 2) then
         ads_data%Un13 = 0.d0
      else if (subun .EQ. 3) then
         ads_data%Un23 = 0.d0
      else
         write (ERROR_UNIT, *) "wrong substep"
      end if
      ads_data%dUn = 0.d0
      total_size = ads%lnelem(1)*ads%lnelem(2)*ads%lnelem(3)

!      loop over points
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(SHARED) &
! !$OMP SHARED(ads,ads_data,total_size) &
! !$OMP PRIVATE(tmp,ex,ey,ez,kx,ky,kz,ind) &
! !$OMP PRIVATE(bx,by,bz,rx,ry,rz,ix,iy,iz,sx,sy,sz,Ucoeff,dvx,dvy,dvz,du) &
! !$OMP PRIVATE(indbx,indby,indbz,Uval,dux,duy,duz,v)
      do all = 1, total_size
!        translate coefficients to local
         ez = modulo(all - 1, ads%lnelem(3))
         tmp = (all - ez)/ads%lnelem(3) + 1
         ey = modulo(tmp - 1, ads%lnelem(2))
         ex = (tmp - ey)/ads%lnelem(2)
!        fix distributed part
         ex = ex + ads%mine(1)
         ey = ey + ads%mine(2)
         ez = ez + ads%mine(3)
!        loop over quadrature points
         do kx = 1, ads%ng(1)
            do ky = 1, ads%ng(2)
               do kz = 1, ads%ng(3)
                  Uval = 0.d0
                  dux = 0.d0
                  duy = 0.d0
                  duz = 0.d0
!                 compute value of derivative from previous time step - du
!                 compute previous solution coefficient at given point - Uval
                  do bx = 0, ads%p(1)
                     do by = 0, ads%p(2)
                        do bz = 0, ads%p(3)
                           indbx = (ads%Ox(ex) + bx)
                           indby = (ads%Oy(ey) + by)
                           indbz = (ads%Oz(ez) + bz)
                           ind = indbx + (indby + indbz*(ads%n(2) + 1))*(ads%n(1) + 1)

                           rx = 2
                           ry = 2
                           rz = 2
                           if (indbx < ads%ibeg(1) - 1) rx = 1
                           if (indbx > ads%iend(1) - 1) rx = 3
                           if (indby < ads%ibeg(2) - 1) ry = 1
                           if (indby > ads%iend(2) - 1) ry = 3
                           if (indbz < ads%ibeg(3) - 1) rz = 1
                           if (indbz > ads%iend(3) - 1) rz = 3

                           ix = indbx - ads%ibegsx(rx) + 1
                           iy = indby - ads%ibegsy(ry) + 1
                           iz = indbz - ads%ibegsz(rz) + 1
                           sx = ads%iendsx(rx) - ads%ibegsx(rx) + 1
                           sy = ads%iendsy(ry) - ads%ibegsy(ry) + 1
                           sz = ads%iendsz(rz) - ads%ibegsz(rz) + 1
                           ind = ix + sx*(iy + sy*iz)

#ifdef IDEBUG
                           if (ind < 0 .or. ind > ads%nrcpp(3)*ads%nrcpp(1)*ads%nrcpp(2) - 1) then
                              write (ERROR_UNIT, *) PRINTRANK, 'Oh crap', ix, iy, iz
                              write (ERROR_UNIT, *) PRINTRANK, 'r', rx, ry, rz
                              write (ERROR_UNIT, *) PRINTRANK, 'x', ads%ibeg(1), ads%iend(1)
                              write (ERROR_UNIT, *) PRINTRANK, 'y', ads%ibeg(2), ads%iend(2)
                              write (ERROR_UNIT, *) PRINTRANK, 'z', ads%ibeg(3), ads%iend(3)
                              write (ERROR_UNIT, *) PRINTRANK, 'sizes=', sx, sy, sz
                              write (ERROR_UNIT, *) PRINTRANK, 'begsx=', ads%ibegsx
                              write (ERROR_UNIT, *) PRINTRANK, 'endsx=', ads%iendsx
                              write (ERROR_UNIT, *) PRINTRANK, 'begsy=', ads%ibegsy
                              write (ERROR_UNIT, *) PRINTRANK, 'endsy=', ads%iendsy
                              write (ERROR_UNIT, *) PRINTRANK, 'begsz=', ads%ibegsz
                              write (ERROR_UNIT, *) PRINTRANK, 'endsz=', ads%iendsz
                           end if
#endif

                           Ucoeff = ads_data%R(ind + 1, rx, ry, rz)
                           v = ads%NNx(0, bx, kx, ex)*ads%NNy(0, by, ky, ey)*ads%NNz(0, bz, kz, ez)
                           dvx = ads%NNx(1, bx, kx, ex)*ads%NNy(0, by, ky, ey)*ads%NNz(0, bz, kz, ez)
                           dvy = ads%NNx(0, bx, kx, ex)*ads%NNy(1, by, ky, ey)*ads%NNz(0, bz, kz, ez)
                           dvz = ads%NNx(0, bx, kx, ex)*ads%NNy(0, by, ky, ey)*ads%NNz(1, bz, kz, ez)

                           Uval = Uval + Ucoeff*v
                           dux = dux + Ucoeff*dvx
                           duy = duy + Ucoeff*dvy
                           duz = duz + Ucoeff*dvz
                        end do
                     end do
                  end do
                  ads_data%dUn(ex, ey, ez, kx, ky, kz, :) = (/dux, duy, duz/)
                  if (subun .EQ. 1) then
                     ads_data%Un(ex, ey, ez, kx, ky, kz) = Uval
                  else if (subun .EQ. 2) then
                     ads_data%Un13(ex, ey, ez, kx, ky, kz) = Uval
                  else if (subun .EQ. 3) then
                     ads_data%Un23(ex, ey, ez, kx, ky, kz) = Uval
                  else
                     write (ERROR_UNIT, *) "wrong substep"
                  end if
               end do
            end do
         end do
      end do
! !$OMP END PARALLEL DO

   end subroutine FormUn

!!!!! to nie tu
!---------------------------------------------------------------------------
!> @author Maciej Wozniak
!>
!> @brief
!> Translates global linearized index given by
!>
!>     ind = (z * (ny+1) + y) * (nx+1) + x
!>
!> to triple \f$ (x,y,z) \f$.
!
! Input:
! ------
!> @param[in] ind        - linearized index
!> @param[in] nx,ny,nz   - sizes of the solution cube minus one
!
! Output:
! -------
!> @param[out] x,y,z      - output coordinates
! -------------------------------------------------------------------
   subroutine global2local(ind, n, x, y, z)
      implicit none
      integer(kind=4), intent(in) :: ind
      integer(kind=4), dimension(3), intent(in) :: n
      integer(kind=4), intent(out) :: x, y, z
      integer(kind=4) :: tmp

      z = ind/((n(1) + 1)*(n(2) + 1))
      tmp = ind - z*(n(1) + 1)*(n(2) + 1)
      y = tmp/(n(1) + 1)
      x = tmp - y*(n(1) + 1)

   end subroutine global2local

!---------------------------------------------------------------------------
!> @author Maciej Wozniak
!>
!> @brief
!> Calculates mass matrix M
!
! Input:
! ------
!> @param[in] U      - knot vector
!> @param[in] p      - degree of approximation
!> @param[in] n      - number of control points minus one
!> @param[in] nelem  - number of subintervals in knot
!> @param[in] mix    - mixing values for MKBBT matrices
!
! Output:
! -------
!> @param[out] sprsmtrx     - sparse matrix, logically \f$ (n+1) \times (n+1) \f$
!
! -------------------------------------------------------------------
   subroutine ComputeMatrix(U1, p1, n1, nelem1, U2, p2, n2, nelem2, mixA, mixB, mixBT, equ, sprsmtrx)
      ! use parallelism, ONLY: PRINTRANK
      use sparse
      implicit none
      integer(kind=4), intent(in) :: n1, p1, nelem1
      real(kind=8), dimension(0:n1 + p1 + 1), intent(in) :: U1
      integer(kind=4), intent(in) :: n2, p2, nelem2
      real(kind=8), dimension(0:n2 + p2 + 1), intent(in) :: U2
      real(kind=8), dimension(4), intent(in) :: mixA, mixB, mixBT
      logical, intent(in) :: equ
      type(sparse_matrix), pointer, intent(out) :: sprsmtrx
      ! integer :: i

      if (equ) then
         call MKBBT_small(nelem2, U2, p2, n2, mixA, sprsmtrx)
      else
         call MKBBT_large(nelem2, U1, p1, n1, U2, p2, n2, mixA, mixB, mixBT, sprsmtrx)
      end if

   end subroutine ComputeMatrix

end module projection_engine

