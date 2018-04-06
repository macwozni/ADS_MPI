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
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(PRIVATE) &
! !$OMP PRIVATE(b,a,k,e,ia,ib,tmp) &
! !$OMP SHARED(nelem,ng,p,O,KL,KU,NN,W,J,total_size) &
! !$OMP REDUCTION(+:M)
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
         M(KL + KU + ia - ib, ib) = M(KL + KU + ia - ib, ib) + NN(0, a, k, e) * NN(0, b, k, e) * J(e) * W(k)


      enddo
! !$OMP END PARALLEL DO

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
      use omp_lib
      implicit none
      interface
         subroutine RHS_fun(&
            ads, &
            X, &
            k, &
            e, &
            a, & 
            ads_data, J, W, ret)
            use Setup, ONLY: ADS_Setup,ADS_compute_data
            implicit none
            type (ADS_setup), intent(in) :: ads
            real (kind = 8), intent(in), dimension(3) :: X
            integer(kind = 4), intent(in), dimension(3) :: k
            integer(kind = 4), intent(in), dimension(3) :: e
            integer(kind = 4), intent(in), dimension(3) :: a
            type (ADS_compute_data), intent(in) :: ads_data
            real (kind = 8), intent(in) :: J, W
            real (kind = 8), intent(out) :: ret
         end subroutine
      end interface
      type (ADS_setup), intent(in) :: ads
      type (ADS_compute_data), intent(inout) :: ads_data
      integer(kind = 4) :: kx, ky, kz, ax, ay, az, ex, ey, ez
      real (kind = 8) :: J, W
      integer(kind = 4) :: ind, ind1, ind23, indx, indy, indz
      real (kind = 8) :: resvalue
      real (kind = 8), dimension(3) :: X
      integer(kind = 4), dimension(3) :: k, e, a
      integer (kind = 4) :: tmp, all
      integer (kind = 4) :: total_size
      real (kind = 8) :: F(ads % s(1), ads % s(2) * ads % s(3))
!      real (kind = 8),allocatable,dimension(:,:) :: F

!      allocate(F(ads % s(1), ads % s(2) * ads % s(3)))

      F = 0.d0

      total_size = ads % lnelem(1) * ads % lnelem(2) * ads % lnelem(3)

!      loop over points
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(SHARED) &
! !$OMP SHARED(ads,ads_data,total_size) &
! !$OMP PRIVATE(tmp,ex,ey,ez,e,kx,ky,kz,k,W,ax,ay,az,a,ind,indx,indy,indz,ind1,ind23,J) &
! !$OMP PRIVATE(X,resvalue) &
! !$OMP REDUCTION(+:F)      
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
!        loop over quadrature points
         do kx = 1, ads % ng(1)
            do ky = 1, ads % ng(2)
               do kz = 1, ads % ng(3)
!                 weigths
                  W = ads % Wx(kx) * ads % Wy(ky) * ads % Wz(kz)
!                 loop over degrees of freedom
                  do ax = 0, ads % p(1)
                     do ay = 0, ads % p(2)
                        do az = 0, ads % p(3)
                           indx = (ads % Ox(ex) + ax)
                           indy = (ads % Oy(ey) + ay)
                           indz = (ads % Oz(ez) + az)
                           ind = indx + (indy + indz*(ads % n(2) + 1))*(ads % n(1) + 1)
                           
                           if (indx < ads % ibeg(1) - 1 .or. indx > ads % iend(1) - 1) cycle
                           if (indy < ads % ibeg(2) - 1 .or. indy > ads % iend(2) - 1) cycle
                           if (indz < ads % ibeg(3) - 1 .or. indz > ads % iend(3) - 1) cycle

                           ind1 = indx - ads % ibeg(1) + 1
                           ind23 = (indy - ads % ibeg(2) + 1) + &
                           (indz - ads % ibeg(3) + 1)*(ads % iend(2) - &
                           ads % ibeg(2) + 1)
                           X = (/ ads % Xx(kx, ex), ads % Xy(ky, ey), ads % Xz(kz, ez) /)
                           k = (/ kx, ky, kz /)
                           e = (/ ex, ey, ez /)
                           a = (/ ax, ay, az /)
                           call RHS_fun(&
                           ads, &
                           X, &
                           k, &
                           e, &
                           a, &
                           ads_data, J, W, resvalue)
!!$OMP FLUSH(F)
                           F(ind1 + 1, ind23 + 1) = F(ind1 + 1, ind23 + 1) + resvalue

                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
! !$OMP END PARALLEL DO

      ads_data % F = F

!if (allocated(F)) deallocate (F)

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

