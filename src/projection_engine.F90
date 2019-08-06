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
!> @param[in] KL1     - number of lower diagonals of the resulting matrix
!> @param[in] KU1     - number of upper diagonals of the resulting matrix
!> @param[in] U1      - knot vector
!> @param[in] p1      - degree of approximation
!> @param[in] n1      - number of control points minus one
!> @param[in] nelem1  - number of subintervals in knot
!> @param[in] KL2     - number of lower diagonals of the resulting matrix
!> @param[in] KU2     - number of upper diagonals of the resulting matrix
!> @param[in] U2      - knot vector
!> @param[in] p2      - degree of approximation
!> @param[in] n2      - number of control points minus one
!> @param[in] nelem2  - number of subintervals in knot
!
! Output:
! -------
!> @param[out] M     - mass matrix, logically \f$ (n+1) \times (n+1) \f$
!> @param[out] K     - stifness matrix, logically \f$ (n+1) \times (n+1) \f$
!> @param[out] B     - advection matrix, logically \f$ (n+1) \times (n+1) \f$
!> @param[out] BT    - advection matrix transposed, logically \f$ (n+1) \times (n+1) \f$
!>
!> Values in the matrix are stored in the band format, i.e. while \f$ M \f$
!> is \f$ (n+1) \times (n+1) \f$, it is stored as \f$ (2 KL + KU + 1) \times n \f$, and the
!> index correspondence is given by:
!>
!>     A(i, j) = M(KL + KU + 1 + i - j, j)
!>
!> \f$ M = u*v \f$
! -------------------------------------------------------------------
subroutine MKBBT(KL1, KU1, U1, p1, n1, nelem1,&
    KL2, KU2, U2, p2, n2,n elem2&
    M,K,B,BT)
   use basis, ONLY: BasisData
   use omp_lib
   implicit none
   integer(kind = 4), intent(in) :: KL1, KU1
   integer(kind = 4), intent(in) :: n1, p1, nelem1
   real (kind = 8), intent(in) :: U1(0:n1 + p1 + 1)
   integer(kind = 4), intent(in) :: KL2, KU2
   integer(kind = 4), intent(in) :: n2, p2, nelem2
   real (kind = 8), intent(in) :: U2(0:n2 + p2 + 1)
   real (kind = 8), dimension(0:(2 * KL + KU), 0:n), intent(out) :: M,K,B,BT
   real (kind = 8), dimension(nelem1) :: J1
   real (kind = 8), dimension(p1 + 1) :: W1
   real (kind = 8), dimension(p1 + 1, nelem1) :: X1
   real (kind = 8), dimension(0:1, 0:p1, p1 + 1, nelem1) :: NN1
   real (kind = 8), dimension(nelem2) :: J2
   real (kind = 8), dimension(p2 + 1) :: W2
   real (kind = 8), dimension(p2 + 1, nelem2) :: X2
   real (kind = 8), dimension(0:1, 0:p2, p2 + 1, nelem2) :: NN2
   integer(kind = 4) :: dd1, dd2
   integer(kind = 4) :: ia, ib
   integer(kind = 4) :: mm1, mm2, ng1, ng2, e, i, c, d
   integer(kind = 4) :: O1(nelem1)
   integer(kind = 4) :: O2(nelem2)
   integer(kind = 4) :: all, tmp, total_size

   mm = n + p + 1
   ng = p + 1
   dd = 1
   M = 0.d0
   K = 0.d0
   B = 0.d0
   BT = 0.d0

   call BasisData(p1, mm1, U1, dd1, ng1, nelem1, O1, J1, W1, X1, NN1)
   call BasisData(p2, mm2, U2, dd2, ng2, nelem2, O2, J2, W2, X2, NN2)

   total_size = nelem * ng * (p + 1)*(p + 1)

! new parallel loop
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) &
!$OMP PRIVATE(d,c,i,e,ia,ib,tmp) &
!$OMP SHARED(nelem,ng,p,O,KL,KU,NN,W,J,total_size) &
!$OMP REDUCTION(+:M) &
!$OMP REDUCTION(+:K) &
!$OMP REDUCTION(+:B) &
!$OMP REDUCTION(+:BT)
   do all = 1, total_size
! loop over shape functions over elements (p+1 functions)
      d = modulo(all - 1, p + 1)
      tmp = (all - d) / (p + 1)
! loop over shape functions over elements (p+1 functions)
      c = modulo(tmp, p + 1)
      tmp = (tmp - c) / (p + 1)
! loop over Gauss points
      i = modulo(tmp, ng) + 1
! loop over elements
      e = (tmp - i + 1) / ng + 1
      ! O(e) + c = first dof of element + 1st local shape function index
      ! O(e) + d = first dof of element + 2nd local shape function index
      ! NN(0,c,i,e) = value of shape function c at Gauss point i over element e
      ! NN(0,d,i,e) = value of shape function d at Gauss point i over element e
      ! W(i) weight for Gauss point i
      ! J(e) jacobian for element e
      ia = O(e) + c
      ib = O(e) + d
      ! M = u*v
      M(KL + KU + ia - ib, ib) = M(KL + KU + ia - ib, ib) + NN(0, c, i, e) * NN(0, d, i, e) * J(e) * W(i)
      K(KL + KU + ia - ib, ib) = K(KL + KU + ia - ib, ib) + NN(1, c, i, e) * NN(1, d, i, e) * J(e) * W(i)
      B(KL + KU + ia - ib, ib) = B(KL + KU + ia - ib, ib) + NN(1, c, i, e) * NN(0, d, i, e) * J(e) * W(i)
      BT(KL + KU + ia - ib, ib) = BT(KL + KU + ia - ib, ib) + NN(0, c, i, e) * NN(1, d, i, e) * J(e) * W(i)


   enddo
!$OMP END PARALLEL DO

end subroutine MKBBT


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
!
! Input/Output
!> @param[inout] ads_data  - data structures for ADS
!
! Output:
! -------
!> @param[out] l2norm      -
! -------------------------------------------------------------------
subroutine Form3DRHS(ads, ads_data, direction, substep,RHS_fun,l2norm)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: PRINTRANK
   use Interfaces, ONLY: RHS_fun_int
   USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
   use omp_lib
   implicit none
   procedure(RHS_fun_int) :: RHS_fun
   type (ADS_setup), intent(in) :: ads
   integer (kind=4), intent(in) :: direction
   integer (kind=4), intent(in) :: substep
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
   real   (kind=8), dimension(3)  :: du
   integer(kind = 4) :: indbx, indby, indbz
   real (kind = 8) :: Uval
   real (kind = 8) :: Uval13
   real (kind = 8) :: Uval23
   real (kind = 8) :: Uval_m(1)
   real (kind = 8), dimension(0:ads % p(1),0:ads % p(2),0:ads % p(3)) :: elarr
   real (kind = 8) :: l2normtmp

   total_size = ads % lnelem(1) * ads % lnelem(2) * ads % lnelem(3)

   l2norm=0.d0
   ads_data % F = 0.d0

!      loop over points
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP SHARED(ads,ads_data,total_size) &
!$OMP PRIVATE(tmp,ex,ey,ez,e,kx,ky,kz,k,W,ax,ay,az,a,ind,indx,indy,indz,ind1,ind23,J) &
!$OMP PRIVATE(X,du,resvalue) &
!$OMP PRIVATE(indbx,indby,indbz,Uval,elarr,l2normtmp,Uval_m,Uval13,Uval23) &
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
               Uval = ads_data % Un(ex,ey,ez,kx,ky,kz)
               Uval13 = ads_data % Un13(ex,ey,ez,kx,ky,kz)
               Uval23 = ads_data % Un23(ex,ey,ez,kx,ky,kz)
               du = ads_data % dUn(ex,ey,ez,kx,ky,kz,:)

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
                           Uval_m = (/Uval/)
                           call RHS_fun(&
                           ads, &
                           X, &
                           k, &
                           e, &
                           a, &
                           du, &
                           1, Uval_m, Uval13,Uval23, &
                           ads_data, J, W, direction, substep, l2normtmp, resvalue)
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
   use parallelism, ONLY: PRINTRANK
   use Interfaces, ONLY: RHS_fun_int
   USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
   use omp_lib
   implicit none
   integer (kind = 4), intent(in) :: subun
   type (ADS_setup), intent(in) :: ads
   type (ADS_compute_data), intent(inout) :: ads_data
   integer(kind = 4) :: kx, ky, kz, ex, ey, ez
   integer(kind = 4) :: ind
   integer (kind = 4) :: tmp, all
   integer (kind = 4) :: total_size
   integer(kind = 4) :: rx, ry, rz, ix, iy, iz, sx, sy, sz
   integer(kind = 4) :: bx, by, bz
   real (kind = 8) :: dux, duy, duz
   real   (kind=8), dimension(3)  :: du
   integer(kind = 4) :: indbx, indby, indbz
   real (kind = 8) :: Uval, ucoeff
   real   (kind=8) :: dvx,dvy,dvz,v


   if (subun .EQ. 1) then
      ads_data % Un = 0.d0
   else if (subun .EQ. 2) then
      ads_data % Un13 = 0.d0
   else if (subun .EQ. 3) then
      ads_data % Un23 = 0.d0
   else
      write(ERROR_UNIT, *) "wrong substep"
   end if
   ads_data % dUn = 0.d0
   total_size = ads % lnelem(1) * ads % lnelem(2) * ads % lnelem(3)

!      loop over points
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP SHARED(ads,ads_data,total_size) &
!$OMP PRIVATE(tmp,ex,ey,ez,kx,ky,kz,ind) &
!$OMP PRIVATE(bx,by,bz,rx,ry,rz,ix,iy,iz,sx,sy,sz,Ucoeff,dvx,dvy,dvz,du) &
!$OMP PRIVATE(indbx,indby,indbz,Uval,dux,duy,duz,v)
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
!        loop over quadrature points
      do kx = 1, ads % ng(1)
         do ky = 1, ads % ng(2)
            do kz = 1, ads % ng(3)
               Uval = 0.d0
               dux = 0.d0
               duy = 0.d0
               duz = 0.d0
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
                        v = ads % NNx(0, bx, kx, ex) * ads % NNy(0, by, ky, ey) * ads % NNz(0, bz, kz, ez)
                        dvx = ads % NNx(1, bx, kx, ex) * ads % NNy(0, by, ky, ey) * ads % NNz(0, bz, kz, ez)
                        dvy = ads % NNx(0, bx, kx, ex) * ads % NNy(1, by, ky, ey) * ads % NNz(0, bz, kz, ez)
                        dvz = ads % NNx(0, bx, kx, ex) * ads % NNy(0, by, ky, ey) * ads % NNz(1, bz, kz, ez)

                        Uval = Uval + Ucoeff * v
                        dux = dux + Ucoeff * dvx
                        duy = duy + Ucoeff * dvy
                        duz = duz + Ucoeff * dvz
                     enddo
                  enddo
               enddo
               ads_data % dUn(ex,ey,ez,kx,ky,kz,:) = (/ dux, duy, duz /)
               if (subun .EQ. 1) then
                  ads_data % Un(ex,ey,ez,kx,ky,kz) = Uval
               else if (subun .EQ. 2) then
                  ads_data % Un13(ex,ey,ez,kx,ky,kz) = Uval
               else if (subun .EQ. 3) then
                  ads_data % Un23(ex,ey,ez,kx,ky,kz) = Uval
               else
                  write(ERROR_UNIT, *) "wrong substep"
               end if
            enddo
         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

end subroutine FormUn


!!!!!!! debug?????
!!!!!!! nie jest uzywany
!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> Checks whether the index is within range.
!
! Input:
! ------
!> @param[in] idn       - 3D index
!> @param[in] ibeg      - 3D range beginning
!> @param[in] iend      - 3D range end
!
! Output:
! -------
!> @return IndexInRange - in 3D range?
! -------------------------------------------------------------------
logical function IndexInRange(ind, ibeg, iend)
   implicit none
   integer(kind = 4), dimension(3), intent(in) :: ind
   integer(kind = 4), dimension(3), intent(in) :: ibeg, iend

   IndexInRange = .true.
   if (ind(1) < ibeg(1) - 1 .or. ind(1) > iend(1) - 1) IndexInRange = .false.
   if (ind(2) < ibeg(2) - 1 .or. ind(2) > iend(2) - 1) IndexInRange = .false.
   if (ind(3) < ibeg(3) - 1 .or. ind(3) > iend(3) - 1) IndexInRange = .false.

end function IndexInRange


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
   integer(kind = 4), intent(in) :: ind
   integer(kind = 4), dimension(3), intent(in) :: n
   integer(kind = 4), intent(out) :: x, y, z
   integer(kind = 4) :: tmp

   z = ind / ((n(1) + 1)*(n(2) + 1))
   tmp = ind - z * (n(1) + 1)*(n(2) + 1)
   y = tmp / (n(1) + 1)
   x = tmp - y * (n(1) + 1)

end subroutine global2local




!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> Calculates mass matrix M
!
! Input:
! ------
!> @param[in] KL     - number of lower diagonals of the resulting matrix
!> @param[in] KU     - number of upper diagonals of the resulting matrix
!> @param[in] U      - knot vector
!> @param[in] p      - degree of approximation
!> @param[in] n      - number of control points minus one
!> @param[in] nelem  - number of subintervals in knot
!> @param[in] mix    - mixing values for MKBBT matrices
!
! Output:
! -------
!> @param[out] O     - matrix, logically \f$ (n+1) \times (n+1) \f$
!
! -------------------------------------------------------------------
subroutine ComputeMatrix(KL, KU, U, p, n, nelem, mix, O)
   use parallelism, ONLY: PRINTRANK
   implicit none
   integer(kind = 4), intent(in) :: KL, KU
   integer(kind = 4), intent(in) :: n, p, nelem
   real (kind = 8), dimension(0:n + p + 1), intent(in) :: U
   real (kind = 8), dimension(4), intent(in) :: mix
   real (kind = 8), dimension(0:(2 * KL + KU), 0:n), intent(out) :: O
   real (kind = 8), dimension(0:(2 * KL + KU), 0:n) :: M,K,B,BT
   integer :: i

   M = 0.0d0
   K = 0.0d0
   B = 0.0d0
   BT = 0.0d0
   
   call MKBBT(KL, KU, U, p, n, nelem, M,K,B,BT)
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

