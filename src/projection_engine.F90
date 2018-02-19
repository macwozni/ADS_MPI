module projection_engine

   implicit none

contains


   ! -------------------------------------------------------------------
   ! Calculates the mass matrix. 
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
   ! -------------------------------------------------------------------
   subroutine Form1DMassMatrix(KL, KU, U, p, n, nelem, M)
      use basis, ONLY: BasisData
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

      mm = n + p + 1
      ng = p + 1
      d = 0
      M = 0

      call BasisData(p, mm, U, d, ng, nelem, O, J, W, X, NN)

      ! loop over elements
      do e = 1, nelem
         ! loop over Gauss points
         do k = 1, ng
            ! loop over shape functions over elements (p+1 functions)
            do a = 0, p
               ! loop over shape functions over elements (p+1 functions)
               do b = 0, p
                  ! O(e) + a = first dof of element + 1st local shape function index
                  ! O(e) + b = first dof of element + 2nd local shape function index
                  ! NN(0,a,k,e) = value of shape function a at Gauss point k over element e
                  ! NN(0,b,k,e) = value of shape function b at Gauss point k over element e
                  ! W(k) weight for Gauss point k
                  ! J(e) jacobian for element e
                  ia = O(e) + a
                  ib = O(e) + b
                  M(KL + KU + ia - ib, ib) = M(KL + KU + ia - ib, ib) + NN(0, a, k, e) * NN(0, b, k, e) * J(e) * W(k)

               enddo
            enddo
         enddo
      enddo

   end subroutine


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
   subroutine Form3DRHS(ads, ads_data, RHS_fun)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use parallelism, ONLY: PRINTRANK
      USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
      use basis, ONLY: BasisData
      implicit none
      interface
         subroutine RHS_fun(&
            ads, &
            X, &
            k, &
            e, &
            a, &
            b, &
            du, &
            NNx, NNy, NNz, &
            Uval, J, W, ret)
            use Setup, ONLY: ADS_Setup
            implicit none
            type (ADS_setup), intent(in) :: ads
            real (kind = 8), intent(in), dimension(3) :: X
            integer(kind = 4), intent(in), dimension(3) :: k
            integer(kind = 4), intent(in), dimension(3) :: e
            integer(kind = 4), intent(in), dimension(3) :: a
            integer(kind = 4), intent(in), dimension(3) :: b
            real (kind = 8), intent(in), dimension(3) :: du
            real (kind = 8), intent(in) :: Uval
            real (kind = 8), intent(in) :: J, W
            real (kind = 8), intent(in) :: &
            NNx(0:1, 0:ads % p(1), ads % p(1) + 1, ads % nelem(1)), &
            NNy(0:1, 0:ads % p(2), ads % p(2) + 1, ads % nelem(2)), &
            NNz(0:1, 0:ads % p(3), ads % p(3) + 1, ads % nelem(3))
            real (kind = 8), intent(out) :: ret
         end subroutine
      end interface
      type (ADS_setup), intent(in) :: ads
      type (ADS_compute_data), intent(inout) :: ads_data
      integer(kind = 4) :: mx, my, mz, ngx, ngy, ngz, ex, ey, ez
      integer(kind = 4) :: kx, ky, kz, ax, ay, az, bx, by, bz, d
      integer(kind = 4) :: Ox(ads % nelem(1)), Oy(ads % nelem(2)), Oz(ads % nelem(3))
      real (kind = 8) :: Jx(ads % nelem(1)), Jy(ads % nelem(2)), Jz(ads % nelem(3))
      real (kind = 8) :: Wx(ads % p(1) + 1), Wy(ads % p(2) + 1), Wz(ads % p(3) + 1)
      real (kind = 8) :: Xx(ads % p(1) + 1, ads % nelem(1))
      real (kind = 8) :: Xy(ads % p(2) + 1, ads % nelem(2))
      real (kind = 8) :: Xz(ads % p(3) + 1, ads % nelem(3))
      real (kind = 8) :: &
      NNx(0:1, 0:ads % p(1), ads % p(1) + 1, ads % nelem(1)), &
      NNy(0:1, 0:ads % p(2), ads % p(2) + 1, ads % nelem(2)), &
      NNz(0:1, 0:ads % p(3), ads % p(3) + 1, ads % nelem(3))
      real (kind = 8) :: J, W, Uval, ucoeff
      real (kind = 8) :: v, rhs
      real (kind = 8) :: dux, duy, duz, dvx, dvy, dvz
      integer(kind = 4) :: nreppx, nreppy, nreppz !# elements per proc along x,y,z
      integer(kind = 4) :: ind, ind1, ind23, indx, indy, indz
      integer(kind = 4) :: indbx, indby, indbz
      integer(kind = 4) :: rx, ry, rz, ix, iy, iz, sx, sy, sz
      real (kind = 8) :: resvalue
      real (kind = 8), dimension(3) :: X, du
      integer(kind = 4), dimension(3) :: k, e, a, b

      d = 0
      mx = ads % n(1) + ads % p(1) + 1
      ngx = ads % p(1) + 1
      my = ads % n(2) + ads % p(2) + 1
      ngy = ads % p(2) + 1
      mz = ads % n(3) + ads % p(3) + 1
      ngz = ads % p(3) + 1

      call BasisData(ads % p(1), mx, ads % Ux, 1, ngx, ads % nelem(1), Ox, Jx, Wx, Xx, NNx)
      call BasisData(ads % p(2), my, ads % Uy, 1, ngy, ads % nelem(2), Oy, Jy, Wy, Xy, NNy)
      call BasisData(ads % p(3), mz, ads % Uz, 1, ngz, ads % nelem(3), Oz, Jz, Wz, Xz, NNz)

#ifdef IPRINT
      write(*, *) PRINTRANK, 'ex:', ads % mine(1), ads % maxe(1)
      write(*, *) PRINTRANK, 'ey:', ads % mine(2), ads % maxe(2)
      write(*, *) PRINTRANK, 'ez:', ads % mine(3), ads % maxe(3)
      write(*, *) PRINTRANK, 'ibegx,iendx', ads % ibeg(1), ads % iend(1)
      write(*, *) PRINTRANK, 'ibegy,iendy', ads % ibeg(2), ads % iend(2)
      write(*, *) PRINTRANK, 'ibegz,iendz', ads % ibeg(3), ads % iend(3)
#endif
      
      ads_data % F = 0

      do ex = ads % mine(1), ads % maxe(1)
         do ey = ads % mine(2), ads % maxe(2)
            do ez = ads % mine(3), ads % maxe(3)
               J = Jx(ex) * Jy(ey) * Jz(ez)
               do kx = 1, ngx
                  do ky = 1, ngy
                     do kz = 1, ngz
                        W = Wx(kx) * Wy(ky) * Wz(kz)
                        do ax = 0, ads % p(1)
                           do ay = 0, ads % p(2)
                              do az = 0, ads % p(3)
                                 ind = (Ox(ex) + ax) + (Oy(ey) + ay)*(ads % n(1) + 1) + (Oz(ez) + az)* &
                                 (ads % n(2) + 1)*(ads % n(1) + 1)
                                 call global2local(ind, ads % n, indx, indy, indz)

                                 if (indx < ads % ibeg(1) - 1 .or. indx > ads % iend(1) - 1) cycle
                                 if (indy < ads % ibeg(2) - 1 .or. indy > ads % iend(2) - 1) cycle
                                 if (indz < ads % ibeg(3) - 1 .or. indz > ads % iend(3) - 1) cycle

                                 ind1 = indx - ads % ibeg(1) + 1
                                 ind23 = (indy - ads % ibeg(2) + 1) + (indz - ads % ibeg(3) + 1)*(ads % iend(2)- &
                                 ads % ibeg(2) + 1)

                                 Uval = 0
                                 dux = 0
                                 duy = 0
                                 duz = 0
                                 do bx = 0, ads % p(1)
                                    do by = 0, ads % p(2)
                                       do bz = 0, ads % p(3)
                                          ind = (Ox(ex) + bx) + (Oy(ey) + by)*(ads % n(1) + 1) + (Oz(ez) + bz)* &
                                          (ads % n(2) + 1)*(ads % n(1) + 1)
                                          call global2local(ind, ads % n, indbx, indby, indbz)

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
                                          write(ERROR_UNIT, *) PRINTRANK, 'idxs', indbx, indby, indbz
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
                                       v = NNx(0, bx, kx, ex) * NNy(0, by, ky, ey) * NNz(0, bz, kz, ez)
                                       dvx = NNx(1, bx, kx, ex) * NNy(0, by, ky, ey) * NNz(0, bz, kz, ez)
                                       dvy = NNx(0, bx, kx, ex) * NNy(1, by, ky, ey) * NNz(0, bz, kz, ez)
                                       dvz = NNx(0, bx, kx, ex) * NNy(0, by, ky, ey) * NNz(1, bz, kz, ez)

                                       Uval = Uval + Ucoeff * v
                                       dux = dux + Ucoeff * dvx
                                       duy = duy + Ucoeff * dvy
                                       duz = duz + Ucoeff * dvz
                                    enddo
                                 enddo
                              enddo
                              X = (/ Xx(kx, ex), Xy(ky, ey), Xz(kz, ez) /)
                              k = (/ kx, ky, kz /)
                              e = (/ ex, ey, ez /)
                              a = (/ ax, ay, az /)
                              b = (/ bx, by, bz/)
                              du = (/ dux, duy, duz /)
                              call RHS_fun(&
                              ads, &
                              X, &
                              k, &
                              e, &
                              a, &
                              b, &
                              du, &
                              NNx, NNy, NNz, &
                              Uval, J, W, resvalue)

                              ads_data % F(ind1 + 1, ind23 + 1) = ads_data % F(ind1 + 1, ind23 + 1) + resvalue

                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

end subroutine



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

end function


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

end subroutine




! -------------------------------------------------------------------
! Calculates mass matrix M
! -------------------------------------------------------------------
subroutine ComputeMassMatrix(KL, KU, U, p, n, nelem, M)
   use parallelism, ONLY: PRINTRANK
   implicit none
   integer(kind = 4), intent(in) :: KL, KU
   integer(kind = 4), intent(in) :: n, p, nelem
   real (kind = 8), intent(in) :: U(0:n + p + 1)
   real (kind = 8), intent(out) :: M(0:(2 * KL + KU), 0:n)
   integer :: i

   call Form1DMassMatrix(KL, KU, U, p, n, nelem, M)
#ifdef IPRINT
   write(*, *) PRINTRANK, 'M'
   do i = 1, 2 * KL + KU !+ 1
      write(*, *) PRINTRANK, M(i, 1:n)
   enddo
#endif
   
end subroutine


end module

