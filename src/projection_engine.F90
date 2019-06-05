module projection_engine

   implicit none

contains


   ! -------------------------------------------------------------------
   ! Calculates the mass matrix M. 
   !
   ! Input:
   ! ------
   ! KL     - number of lower diagonals of the resulting matrix
   ! KU     - number of upper diagonals of the resulting matrix
   ! U      - knot vector
   ! p      - degree of approximation
   ! n      - number of control points minus one
   ! nelem  - number of subintervals in knot
   !
   ! Output:
   ! -------
   ! M      - mass matrix, logically (n+1) x (n+1)
   !
   ! Values in the matrix are stored in the band format, i.e. while M
   ! is (n+1) x (n+1), it is stored as (2 KL + KU + 1) x n, and the
   ! index correspondence is given by:
   !
   !     A(i, j) = M(KL + KU + 1 + i - j, j)
   !
   ! M = u*v
   !
   ! -------------------------------------------------------------------
   subroutine Form1DMassMatrix(KL, KU, U, p, n, nelem, M)
      use basis, ONLY: BasisData
      use omp_lib
      implicit none
      integer(kind = 4), intent(in) :: KL, KU
      integer(kind = 4), intent(in) :: n, p, nelem
      real (kind = 8), intent(in) :: U(0:n + p + 1)
      real (kind = 8), intent(out) :: M(0:(2 * KL + KU), 0:n)
      real (kind = 8) :: J(nelem)
      real (kind = 8) :: W(p + 1)
      real (kind = 8) :: X(p + 1, nelem)
      real (kind = 8) :: NN(0:0, 0:p, p + 1, nelem)
      integer(kind = 4) :: d
      integer(kind = 4) :: ia, ib
      integer(kind = 4) :: mm, ng, e, k, a, b
      integer(kind = 4) :: O(nelem)
      integer(kind = 4) :: all, tmp, total_size

      mm = n + p + 1
      ng = p + 1
      d = 0
      M = 0

      call BasisData(p, mm, U, d, ng, nelem, O, J, W, X, NN)

      total_size = nelem * ng * (p + 1)*(p + 1)

! new parallel loop
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) &
!$OMP PRIVATE(b,a,k,e,ia,ib,tmp) &
!$OMP SHARED(nelem,ng,p,O,KL,KU,NN,W,J,total_size) &
!$OMP REDUCTION(+:M)
      do all = 1, total_size
! loop over shape functions over elements (p+1 functions)
         b = modulo(all - 1, p + 1)
         tmp = (all - b) / (p + 1)
! loop over shape functions over elements (p+1 functions)
         a = modulo(tmp, p + 1)
         tmp = (tmp - a) / (p + 1)
! loop over Gauss points
         k = modulo(tmp, ng) + 1
! loop over elements
         e = (tmp - k + 1) / ng + 1
         ! O(e) + a = first dof of element + 1st local shape function index
         ! O(e) + b = first dof of element + 2nd local shape function index
         ! NN(0,a,k,e) = value of shape function a at Gauss point k over element e
         ! NN(0,b,k,e) = value of shape function b at Gauss point k over element e
         ! W(k) weight for Gauss point k
         ! J(e) jacobian for element e
         ia = O(e) + a
         ib = O(e) + b
         ! M = u*v
         M(KL + KU + ia - ib, ib) = M(KL + KU + ia - ib, ib) + NN(0, a, k, e) * NN(0, b, k, e) * J(e) * W(k)


      enddo
!$OMP END PARALLEL DO

   end subroutine Form1DMassMatrix


   !!!!!!!!!!!!!!! TODO !!!!!!!!!!!!!!!!!!!!
   ! -------------------------------------------------------------------
   ! Calculates the stifness matrix M. 
   !
   ! Input:
   ! ------
   ! KL     - number of lower diagonals of the resulting matrix
   ! KU     - number of upper diagonals of the resulting matrix
   ! U      - knot vector
   ! p      - degree of approximation
   ! n      - number of control points minus one
   ! nelem  - number of subintervals in knot
   !
   ! Output:
   ! -------
   ! M      - stifness matrix, logically (n+1) x (n+1)
   !
   ! Values in the matrix are stored in the band format, i.e. while M
   ! is (n+1) x (n+1), it is stored as (2 KL + KU + 1) x n, and the
   ! index correspondence is given by:
   !
   !     A(i, j) = M(KL + KU + 1 + i - j, j)
   !
   ! M = du/dx*dv/dx
   !
   ! -------------------------------------------------------------------
   subroutine Form1DStifnessMatrix(KL, KU, U, p, n, nelem, M)
      use basis, ONLY: BasisData
      use omp_lib
      implicit none
      integer(kind = 4), intent(in) :: KL, KU
      integer(kind = 4), intent(in) :: n, p, nelem
      real (kind = 8), intent(in) :: U(0:n + p + 1)
      real (kind = 8), intent(out) :: M(0:(2 * KL + KU), 0:n)
      real (kind = 8) :: J(nelem)
      real (kind = 8) :: W(p + 1)
      real (kind = 8) :: X(p + 1, nelem)
      real (kind = 8) :: NN(0:1, 0:p, p + 1, nelem)
      integer(kind = 4) :: d
      integer(kind = 4) :: ia, ib
      integer(kind = 4) :: mm, ng, e, k, a, b
      integer(kind = 4) :: O(nelem)
      integer(kind = 4) :: all, tmp, total_size

      mm = n + p + 1
      ng = p + 1
      d = 1
      M = 0

      call BasisData(p, mm, U, d, ng, nelem, O, J, W, X, NN)

      total_size = nelem * ng * (p + 1)*(p + 1)

! new parallel loop
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) &
!$OMP PRIVATE(b,a,k,e,ia,ib,tmp) &
!$OMP SHARED(nelem,ng,p,O,KL,KU,NN,W,J,total_size) &
!$OMP REDUCTION(+:M)
      do all = 1, total_size
! loop over shape functions over elements (p+1 functions)
         b = modulo(all - 1, p + 1)
         tmp = (all - b) / (p + 1)
! loop over shape functions over elements (p+1 functions)
         a = modulo(tmp, p + 1)
         tmp = (tmp - a) / (p + 1)
! loop over Gauss points
         k = modulo(tmp, ng) + 1
! loop over elements
         e = (tmp - k + 1) / ng + 1
         ! O(e) + a = first dof of element + 1st local shape function index
         ! O(e) + b = first dof of element + 2nd local shape function index
         ! NN(0,a,k,e) = value of shape function a at Gauss point k over element e
         ! NN(0,b,k,e) = value of shape function b at Gauss point k over element e
         ! W(k) weight for Gauss point k
         ! J(e) jacobian for element e
         ia = O(e) + a
         ib = O(e) + b
         ! M = du/dx*dv/dx
         M(KL + KU + ia - ib, ib) = M(KL + KU + ia - ib, ib) + NN(1, a, k, e) * NN(1, b, k, e) * J(e) * W(k)


      enddo
!$OMP END PARALLEL DO

   end subroutine Form1DStifnessMatrix
   
   
   !!!!!!!!!!!!!!! TODO !!!!!!!!!!!!!!!!!!!!
   ! -------------------------------------------------------------------
   ! Calculates the advection matrix M. 
   !
   ! Input:
   ! ------
   ! KL     - number of lower diagonals of the resulting matrix
   ! KU     - number of upper diagonals of the resulting matrix
   ! U      - knot vector
   ! p      - degree of approximation
   ! n      - number of control points minus one
   ! nelem  - number of subintervals in knot
   !
   ! Output:
   ! -------
   ! M      - advection matrix, logically (n+1) x (n+1)
   !
   ! Values in the matrix are stored in the band format, i.e. while M
   ! is (n+1) x (n+1), it is stored as (2 KL + KU + 1) x n, and the
   ! index correspondence is given by:
   !
   !     A(i, j) = M(KL + KU + 1 + i - j, j)
   !
   ! M = du/dx*v
   !
   ! -------------------------------------------------------------------
   subroutine Form1DAdvectionMatrix(KL, KU, U, p, n, nelem, M)
      use basis, ONLY: BasisData
      use omp_lib
      implicit none
      integer(kind = 4), intent(in) :: KL, KU
      integer(kind = 4), intent(in) :: n, p, nelem
      real (kind = 8), intent(in) :: U(0:n + p + 1)
      real (kind = 8), intent(out) :: M(0:(2 * KL + KU), 0:n)
      real (kind = 8) :: J(nelem)
      real (kind = 8) :: W(p + 1)
      real (kind = 8) :: X(p + 1, nelem)
      real (kind = 8) :: NN(0:1, 0:p, p + 1, nelem)
      integer(kind = 4) :: d
      integer(kind = 4) :: ia, ib
      integer(kind = 4) :: mm, ng, e, k, a, b
      integer(kind = 4) :: O(nelem)
      integer(kind = 4) :: all, tmp, total_size

      mm = n + p + 1
      ng = p + 1
      d = 1
      M = 0

      call BasisData(p, mm, U, d, ng, nelem, O, J, W, X, NN)

      total_size = nelem * ng * (p + 1)*(p + 1)

! new parallel loop
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) &
!$OMP PRIVATE(b,a,k,e,ia,ib,tmp) &
!$OMP SHARED(nelem,ng,p,O,KL,KU,NN,W,J,total_size) &
!$OMP REDUCTION(+:M)
      do all = 1, total_size
! loop over shape functions over elements (p+1 functions)
         b = modulo(all - 1, p + 1)
         tmp = (all - b) / (p + 1)
! loop over shape functions over elements (p+1 functions)
         a = modulo(tmp, p + 1)
         tmp = (tmp - a) / (p + 1)
! loop over Gauss points
         k = modulo(tmp, ng) + 1
! loop over elements
         e = (tmp - k + 1) / ng + 1
         ! O(e) + a = first dof of element + 1st local shape function index
         ! O(e) + b = first dof of element + 2nd local shape function index
         ! NN(0,a,k,e) = value of shape function a at Gauss point k over element e
         ! NN(0,b,k,e) = value of shape function b at Gauss point k over element e
         ! W(k) weight for Gauss point k
         ! J(e) jacobian for element e
         ia = O(e) + a
         ib = O(e) + b
         ! M = du/dx*v
         M(KL + KU + ia - ib, ib) = M(KL + KU + ia - ib, ib) + NN(1, a, k, e) * NN(0, b, k, e) * J(e) * W(k)


      enddo
!$OMP END PARALLEL DO

   end subroutine Form1DAdvectionMatrix
   
   
   !!!!!!!!!!!!!!! TODO !!!!!!!!!!!!!!!!!!!!
   ! -------------------------------------------------------------------
   ! Calculates the advection matrix M transposed. 
   !
   ! Input:
   ! ------
   ! KL     - number of lower diagonals of the resulting matrix
   ! KU     - number of upper diagonals of the resulting matrix
   ! U      - knot vector
   ! p      - degree of approximation
   ! n      - number of control points minus one
   ! nelem  - number of subintervals in knot
   !
   ! Output:
   ! -------
   ! M      - advection matrix transposed, logically (n+1) x (n+1)
   !
   ! Values in the matrix are stored in the band format, i.e. while M
   ! is (n+1) x (n+1), it is stored as (2 KL + KU + 1) x n, and the
   ! index correspondence is given by:
   !
   !     A(i, j) = M(KL + KU + 1 + i - j, j)
   !
   ! M = du/dx*v
   !
   ! -------------------------------------------------------------------
   subroutine Form1DAdvectionMatrixT(KL, KU, U, p, n, nelem, M)
      use basis, ONLY: BasisData
      use omp_lib
      implicit none
      integer(kind = 4), intent(in) :: KL, KU
      integer(kind = 4), intent(in) :: n, p, nelem
      real (kind = 8), intent(in) :: U(0:n + p + 1)
      real (kind = 8), intent(out) :: M(0:(2 * KL + KU), 0:n)
      real (kind = 8) :: J(nelem)
      real (kind = 8) :: W(p + 1)
      real (kind = 8) :: X(p + 1, nelem)
      real (kind = 8) :: NN(0:1, 0:p, p + 1, nelem)
      integer(kind = 4) :: d
      integer(kind = 4) :: ia, ib
      integer(kind = 4) :: mm, ng, e, k, a, b
      integer(kind = 4) :: O(nelem)
      integer(kind = 4) :: all, tmp, total_size

      mm = n + p + 1
      ng = p + 1
      d = 1
      M = 0

      call BasisData(p, mm, U, d, ng, nelem, O, J, W, X, NN)

      total_size = nelem * ng * (p + 1)*(p + 1)

! new parallel loop
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) &
!$OMP PRIVATE(b,a,k,e,ia,ib,tmp) &
!$OMP SHARED(nelem,ng,p,O,KL,KU,NN,W,J,total_size) &
!$OMP REDUCTION(+:M)
      do all = 1, total_size
! loop over shape functions over elements (p+1 functions)
         b = modulo(all - 1, p + 1)
         tmp = (all - b) / (p + 1)
! loop over shape functions over elements (p+1 functions)
         a = modulo(tmp, p + 1)
         tmp = (tmp - a) / (p + 1)
! loop over Gauss points
         k = modulo(tmp, ng) + 1
! loop over elements
         e = (tmp - k + 1) / ng + 1
         ! O(e) + a = first dof of element + 1st local shape function index
         ! O(e) + b = first dof of element + 2nd local shape function index
         ! NN(0,a,k,e) = value of shape function a at Gauss point k over element e
         ! NN(0,b,k,e) = value of shape function b at Gauss point k over element e
         ! W(k) weight for Gauss point k
         ! J(e) jacobian for element e
         ia = O(e) + a
         ib = O(e) + b
         ! M = du/dx*v
         M(KL + KU + ia - ib, ib) = M(KL + KU + ia - ib, ib) + NN(0, a, k, e) * NN(1, b, k, e) * J(e) * W(k)


      enddo
!$OMP END PARALLEL DO

   end subroutine Form1DAdvectionMatrixT
   
   
   ! -------------------------------------------------------------------
   ! Calculate right-hand side of the equation.
   !
   ! Input:
   ! ------
   ! U_              - knot vectors
   ! p_              - degrees of approximation
   ! n_              - numbers of functions minus one
   ! nelem_          - number of subintervals
   ! nrcpp_          - number of basis functions per process
   ! ibeg_, iend_    - piece of domain associated with this process
   ! ibegs_, iends_  - pieces of domain surrounding this process' piece
   ! mine_, maxe_    - indices of first and last elements in each direction
   ! R               - previous solution coefficients
   !
   ! Output:
   ! -------
   ! F               - output rhs (multiple vectors)
   !
   ! R has two 'kinds' of dimensions - it's orgainzed as 3x3x3 array of
   ! domain pieces.
   ! -------------------------------------------------------------------
   subroutine Form3DRHS(ads, ads_data, RHS_fun,l2norm)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use parallelism, ONLY: PRINTRANK
      use Interfaces, ONLY: RHS_fun_int
      USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
      use omp_lib
      implicit none
      procedure(RHS_fun_int) :: RHS_fun
      type (ADS_setup), intent(in) :: ads
      type (ADS_compute_data), intent(inout) :: ads_data
      real (kind = 8), intent(out) :: l2norm
      integer(kind = 4) :: kx, ky, kz, ax, ay, az, ex, ey, ez
      real (kind = 8) :: J, W
      integer(kind = 4) :: ind, ind1, ind23, indx, indy, indz
      real (kind = 8) :: resvalue
      real (kind = 8), dimension(3) :: X
      integer(kind = 4), dimension(3) :: k, e, a
      integer (kind = 4) :: tmp, all
      integer (kind = 4) :: total_size
      integer(kind = 4) :: rx, ry, rz, ix, iy, iz, sx, sy, sz
      integer(kind = 4) :: bx, by, bz
      real   (kind=8), dimension(3)  :: du
      real (kind = 8) :: dux, duy, duz
      integer(kind = 4) :: mx, my, mz, ngx, ngy, ngz
      integer(kind = 4) :: indbx, indby, indbz
      real (kind = 8) :: Uval, ucoeff
      real   (kind=8) :: dvx,dvy,dvz,v
      integer(kind = 4) :: threadcnt,threadid
      real (kind = 8) :: elarr(0:ads % p(1),0:ads % p(2),0:ads % p(3))
      real (kind = 8) :: l2normtmp

      total_size = ads % lnelem(1) * ads % lnelem(2) * ads % lnelem(3)

      l2norm=0.d0
      ads_data % F = 0.d0
      
!      loop over points
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP SHARED(ads,ads_data,total_size) &
!$OMP PRIVATE(tmp,ex,ey,ez,e,kx,ky,kz,k,W,ax,ay,az,a,ind,indx,indy,indz,ind1,ind23,J) &
!$OMP PRIVATE(bx,by,bz,rx,ry,rz,ix,iy,iz,sx,sy,sz,Ucoeff,dvx,dvy,dvz,X,du,resvalue) &
!$OMP PRIVATE(indbx,indby,indbz,Uval,dux,duy,duz,v,threadid,elarr,l2normtmp) &
!$OMP REDUCTION(+:l2norm)
      do all = 1, total_size
!        translate coefficients to local
         ez = modulo(all - 1, ads % lnelem(3))
         tmp = (all - ez)/ads % lnelem(3) + 1
         ey = modulo(tmp - 1, ads % lnelem(2))
         ex = (tmp - ey)/ads % lnelem(2)
!        fix distributed part
         ex = ex + ads % mine(1)
         ey = ey + ads % mine(2)
         ez = ez + ads % mine(3)
!        Jacobian
         J = ads % Jx(ex) * ads % Jy(ey) * ads % Jz(ez)
         e = (/ ex, ey, ez /)
         elarr = 0.d0
!        loop over quadrature points
         do kx = 1, ads % ng(1)
            do ky = 1, ads % ng(2)
               do kz = 1, ads % ng(3)
                  k = (/ kx, ky, kz /)
!                 weigths
                  W = ads % Wx(kx) * ads % Wy(ky) * ads % Wz(kz)
                  Uval = 0
                  dux = 0
                  duy = 0
                  duz = 0
!                 compute value of derivative from previous time step - du
!                 compute previous solution coefficient at given point - Uval
                  do bx = 0, ads % p(1)
                     do by = 0, ads % p(2)
                        do bz = 0, ads % p(3)
                           indbx = (ads % Ox(ex) + bx)
                           indby = (ads % Oy(ey) + by)
                           indbz = (ads % Oz(ez) + bz)
                           ind = indbx + (indby + indbz*(ads % n(2) + 1))*(ads % n(1) + 1)

                           rx = 2
                           ry = 2
                           rz = 2
                           if (indbx < ads % ibeg(1) - 1) rx = 1
                           if (indbx > ads % iend(1) - 1) rx = 3
                           if (indby < ads % ibeg(2) - 1) ry = 1
                           if (indby > ads % iend(2) - 1) ry = 3
                           if (indbz < ads % ibeg(3) - 1) rz = 1
                           if (indbz > ads % iend(3) - 1) rz = 3

                           ix = indbx - ads % ibegsx(rx) + 1
                           iy = indby - ads % ibegsy(ry) + 1
                           iz = indbz - ads % ibegsz(rz) + 1
                           sx = ads % iendsx(rx) - ads % ibegsx(rx) + 1
                           sy = ads % iendsy(ry) - ads % ibegsy(ry) + 1
                           sz = ads % iendsz(rz) - ads % ibegsz(rz) + 1
                           ind = ix + sx * (iy + sy * iz)

#ifdef IDEBUG
                           if (ind < 0 .or. ind > ads % nrcpp(3) * ads % nrcpp(1) * ads % nrcpp(2) - 1) then
                              write(ERROR_UNIT, *) PRINTRANK, 'Oh crap', ix, iy, iz
                              write(ERROR_UNIT, *) PRINTRANK, 'r', rx, ry, rz
                              write(ERROR_UNIT, *) PRINTRANK, 'x', ads % ibeg(1), ads % iend(1)
                              write(ERROR_UNIT, *) PRINTRANK, 'y', ads % ibeg(2), ads % iend(2)
                              write(ERROR_UNIT, *) PRINTRANK, 'z', ads % ibeg(3), ads % iend(3)
                              write(ERROR_UNIT, *) PRINTRANK, 'idxs', indx, indy, indz
                              write(ERROR_UNIT, *) PRINTRANK, 'sizes=', sx, sy, sz
                              write(ERROR_UNIT, *) PRINTRANK, 'begsx=', ads % ibegsx
                              write(ERROR_UNIT, *) PRINTRANK, 'endsx=', ads % iendsx
                              write(ERROR_UNIT, *) PRINTRANK, 'begsy=', ads % ibegsy
                              write(ERROR_UNIT, *) PRINTRANK, 'endsy=', ads % iendsy
                              write(ERROR_UNIT, *) PRINTRANK, 'begsz=', ads % ibegsz
                              write(ERROR_UNIT, *) PRINTRANK, 'endsz=', ads % iendsz
                           endif
#endif

                           Ucoeff = ads_data % R(ind + 1, rx, ry, rz)
                           v = ads % NNx(0, bx, k(1), e(1)) * ads % NNy(0, by, k(2), e(2)) * ads % NNz(0, bz, k(3), e(3))
                           dvx = ads % NNx(1, bx, k(1), e(1)) * ads % NNy(0, by, k(2), e(2)) * ads % NNz(0, bz, k(3), e(3))
                           dvy = ads % NNx(0, bx, k(1), e(1)) * ads % NNy(1, by, k(2), e(2)) * ads % NNz(0, bz, k(3), e(3))
                           dvz = ads % NNx(0, bx, k(1), e(1)) * ads % NNy(0, by, k(2), e(2)) * ads % NNz(1, bz, k(3), e(3))

                           Uval = Uval + Ucoeff * v
                           dux = dux + Ucoeff * dvx
                           duy = duy + Ucoeff * dvy
                           duz = duz + Ucoeff * dvz
                        enddo
                     enddo
                  enddo
                  du = (/ dux, duy, duz /)

!                 loop over degrees of freedom
                  do ax = 0, ads % p(1)
                     do ay = 0, ads % p(2)
                        do az = 0, ads % p(3)
                           indx = (ads % Ox(ex) + ax)
                           indy = (ads % Oy(ey) + ay)
                           indz = (ads % Oz(ez) + az)
                           ind = indx + (indy + indz*(ads % n(2) + 1))*(ads % n(1) + 1)
                           
                           if ((indx < ads % ibeg(1) - 1) .or. (indx > ads % iend(1) - 1) .or. &
                           (indy < ads % ibeg(2) - 1) .or. (indy > ads % iend(2) - 1) .or. &
                           (indz < ads % ibeg(3) - 1) .or. (indz > ads % iend(3) - 1)) then
                           else 
                              ind1 = indx - ads % ibeg(1) + 1
                              ind23 = (indy - ads % ibeg(2) + 1) + &
                              (indz - ads % ibeg(3) + 1)*(ads % iend(2) - ads % ibeg(2) + 1)
                              X = (/ ads % Xx(kx, ex), ads % Xy(ky, ey), ads % Xz(kz, ez) /)
                              a = (/ ax, ay, az /)
                              call RHS_fun(&
                              ads, &
                              X, &
                              k, &
                              e, &
                              a, &
                              du, &
                              1, (/Uval/), 0.0d0, 0.0d0, &
                              ads_data, J, W, 1, 1, l2normtmp, resvalue)
                              elarr(ax,ay,az) = elarr(ax,ay,az) + resvalue
                              l2norm = l2norm + l2normtmp
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
!        moving results from temporary array to main one
!$OMP CRITICAL
         do ax = 0, ads % p(1)
            do ay = 0, ads % p(2)
               do az = 0, ads % p(3)
                  indx = (ads % Ox(ex) + ax)
                  indy = (ads % Oy(ey) + ay)
                  indz = (ads % Oz(ez) + az)
                  ind1 = indx - ads % ibeg(1) + 1
                  ind23 = (indy - ads % ibeg(2) + 1) + &
                  (indz - ads % ibeg(3) + 1)*(ads % iend(2) - ads % ibeg(2) + 1)
                  if ((indx < ads % ibeg(1) - 1) .or. (indx > ads % iend(1) - 1) .or. &
                  (indy < ads % ibeg(2) - 1) .or. (indy > ads % iend(2) - 1) .or. &
                  (indz < ads % ibeg(3) - 1) .or. (indz > ads % iend(3) - 1)) then
                  else 
                     ads_data % F(ind1 + 1, ind23 + 1) = ads_data % F(ind1 + 1, ind23 + 1) + elarr(ax,ay,az)
                  endif
               enddo
            enddo
         enddo
!$OMP END CRITICAL
      enddo
!$OMP END PARALLEL DO

   end subroutine Form3DRHS



!!!!!!! debug?????
!!!!!!! nie jest uzywany
! -------------------------------------------------------------------
! Checks whether the index is within range.
!
! idnx,indy,indz   - 3D index
! ibeg(1),iend(1)      - x range
! ibeg(2),iend(2)      - y range
! ibeg(3),iend(3)      - z range
! -------------------------------------------------------------------
logical function IndexInRange(ind, ibeg, iend)
   implicit none
   integer(kind = 4), intent(in), dimension(3) :: ind
   integer(kind = 4), intent(in), dimension(3) :: ibeg, iend

   IndexInRange = .true.
   if (ind(1) < ibeg(1) - 1 .or. ind(1) > iend(1) - 1) IndexInRange = .false.
   if (ind(2) < ibeg(2) - 1 .or. ind(2) > iend(2) - 1) IndexInRange = .false.
   if (ind(3) < ibeg(3) - 1 .or. ind(3) > iend(3) - 1) IndexInRange = .false.

end function IndexInRange


!!!!! to nie tu
! -------------------------------------------------------------------
! Translates global linearized index given by
!
!     ind = (z * (ny+1) + y) * (nx+1) + x
!
! to triple (x,y,z).
!
! ind        - linearized index
! nx,ny,nz   - sizes of the solution cube minus one
! x,y,z      - output coordinates
! -------------------------------------------------------------------
subroutine global2local(ind, n, x, y, z)
   implicit none
   integer(kind = 4), intent(in) :: ind
   integer(kind = 4), intent(in), dimension(3) :: n
   integer(kind = 4), intent(out) :: x, y, z
   integer(kind = 4) :: tmp

   z = ind / ((n(1) + 1)*(n(2) + 1))
   tmp = ind - z * (n(1) + 1)*(n(2) + 1)
   y = tmp / (n(1) + 1)
   x = tmp - y * (n(1) + 1)

end subroutine global2local




! -------------------------------------------------------------------
! Calculates mass matrix M
!
! Input:
! ------
! KL     - number of lower diagonals of the resulting matrix
! KU     - number of upper diagonals of the resulting matrix
! U      - knot vector
! p      - degree of approximation
! n      - number of control points minus one
! nelem  - number of subintervals in knot
! MKAT    - logical values determining if M-mass, K-stifness, A-advection matrices should be computed
! mix    - mixing values for MKAT matrices
!
! Output:
! -------
! O      - matrix, logically (n+1) x (n+1)
!
! -------------------------------------------------------------------
subroutine ComputeMatrix(KL, KU, U, p, n, nelem, MKAT, mix, O)
   use parallelism, ONLY: PRINTRANK
   implicit none
   integer(kind = 4), intent(in) :: KL, KU
   integer(kind = 4), intent(in) :: n, p, nelem
   real (kind = 8), intent(in) :: U(0:n + p + 1)
   logical, intent(in) :: MKAT(4)
   real (kind = 8), intent(in) :: mix(4)
   real (kind = 8), intent(out) :: O(0:(2 * KL + KU), 0:n)
   real (kind = 8) :: M(0:(2 * KL + KU), 0:n)
   real (kind = 8) :: K(0:(2 * KL + KU), 0:n)
   real (kind = 8) :: B(0:(2 * KL + KU), 0:n)
   real (kind = 8) :: BT(0:(2 * KL + KU), 0:n)
   integer :: i

   M = 0.0d0
   K = 0.0d0
   B = 0.0d0
   BT = 0.0d0
   
   if (MKAT(1)) call Form1DMassMatrix(KL, KU, U, p, n, nelem, M)
   if (MKAT(2)) call Form1DStifnessMatrix(KL, KU, U, p, n, nelem, K)
   if (MKAT(3)) call Form1DAdvectionMatrix(KL, KU, U, p, n, nelem, B)
   if (MKAT(4)) call Form1DAdvectionMatrixT(KL, KU, U, p, n, nelem, BT)
#ifdef IPRINT
   write(*, *) PRINTRANK, 'M'
   do i = 1, 2 * KL + KU !+ 1
      write(*, *) PRINTRANK, M(i, 1:n)
   enddo
   write(*, *) PRINTRANK, 'K'
   do i = 1, 2 * KL + KU !+ 1
      write(*, *) PRINTRANK, K(i, 1:n)
   enddo
   write(*, *) PRINTRANK, 'B'
   do i = 1, 2 * KL + KU !+ 1
      write(*, *) PRINTRANK, B(i, 1:n)
   enddo
   write(*, *) PRINTRANK, 'BT'
   do i = 1, 2 * KL + KU !+ 1
      write(*, *) PRINTRANK, BT(i, 1:n)
   enddo
#endif
   
   O = mix(1)*M + mix(2)*K + mix(3)*B + mix(4)*BT
   
#ifdef IPRINT
   write(*, *) PRINTRANK, 'O'
   do i = 1, 2 * KL + KU !+ 1
      write(*, *) PRINTRANK, O(i, 1:n)
   enddo
#endif

end subroutine ComputeMatrix


end module projection_engine

