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
subroutine Form1DMassMatrix(KL,KU,U,p,n,nelem,M)
use basis, ONLY : BasisData
implicit none
integer(kind=4), intent(in)  :: KL,KU
integer(kind=4), intent(in)  :: n, p, nelem
real   (kind=8), intent(in)  :: U(0:n+p+1)
real   (kind=8), intent(out) :: M(0:(2*KL+KU),0:n)
real   (kind=8) :: J(nelem)
real   (kind=8) :: W(p+1)
real   (kind=8) :: X(p+1,nelem)
real   (kind=8) :: NN(0:0,0:p,p+1,nelem)
integer(kind=4) :: d
integer(kind=4) :: ia, ib
integer(kind=4) :: mm,ng,e,k,a,b
integer(kind=4) :: O(nelem)

  mm = n+p+1
  ng = p+1
  d = 0
  M = 0

  call BasisData(p,mm,U,d,ng,nelem,O,J,W,X,NN) 

  ! loop over elements
  do e = 1,nelem
    ! loop over Gauss points
    do k = 1,ng
      ! loop over shape functions over elements (p+1 functions)
      do a = 0,p
        ! loop over shape functions over elements (p+1 functions)
        do b = 0,p
          ! O(e) + a = first dof of element + 1st local shape function index
          ! O(e) + b = first dof of element + 2nd local shape function index
          ! NN(0,a,k,e) = value of shape function a at Gauss point k over element e
          ! NN(0,b,k,e) = value of shape function b at Gauss point k over element e
          ! W(k) weight for Gauss point k
          ! J(e) jacobian for element e
          ia = O(e) + a
          ib = O(e) + b
          M(KL+KU+ia-ib,ib) = M(KL+KU+ia-ib,ib) + NN(0,a,k,e)*NN(0,b,k,e)*J(e)*W(k)

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
subroutine Form3DRHS(          &
   p,n,nelem,nrcpp,     &
   ibeg,iend,                &
   mine,maxe,                &
   Ux,Uy,Uz,     &
   ibegsx,iendsx,              &
   ibegsy,iendsy,              &
   ibegsz,iendsz,              &
   R,F,RHS_fun)
use parallelism, ONLY : PRINTRANK
USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT ! access computing environment
use basis, ONLY : BasisData
implicit none
interface
  subroutine RHS_fun( &
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
         NNx,NNy,NNz, &
         Uval,J,W,F)
      implicit none
      integer(kind=4), intent(in), dimension(3)  :: p
      real   (kind=8), intent(in), dimension(3)  :: X
      integer(kind=4), intent(in), dimension(3)  :: k
      integer(kind=4), intent(in), dimension(3)  :: e
      integer(kind=4), intent(in), dimension(3)  :: nelem
      integer(kind=4), intent(in), dimension(3)  :: ibeg
      integer(kind=4), intent(in), dimension(3)  :: iend
      integer(kind=4), intent(in), dimension(3)  :: maxe
      integer(kind=4), intent(in), dimension(3)  :: mine
      integer(kind=4), intent(in), dimension(3)  :: a
      integer(kind=4), intent(in), dimension(3)  :: b
      real   (kind=8), intent(in), dimension(3)  :: du
      real   (kind=8), intent(in)  :: Uval
      real   (kind=8), intent(in)  :: J,W
      real   (kind=8), intent(in)  :: NNx(0:p(1)-1,0:p(1),p(1)+1,nelem(1)), &
                           NNy(0:p(2)-1,0:p(2),p(2)+1,nelem(2)), &
                           NNz(0:p(3)-1,0:p(3),p(3)+1,nelem(3))
      real   (kind=8), intent(out) :: F
  end subroutine
end interface
integer(kind=4), intent(in), dimension(3)  :: n, p, nelem, nrcpp
integer(kind=4), intent(in), dimension(3)  :: mine,maxe
integer(kind=4), intent(in), dimension(3)  :: ibeg,iend
real   (kind=8), intent(in)  :: Ux(0:n(1)+p(1)+1)
real   (kind=8), intent(in)  :: Uy(0:n(2)+p(2)+1)
real   (kind=8), intent(in)  :: Uz(0:n(3)+p(3)+1)
real   (kind=8), intent(in)  :: R(0:nrcpp(3)*nrcpp(1)*nrcpp(2)-1,3,3,3)
integer(kind=4), intent(in)  :: ibegsx(3),iendsx(3),ibegsy(3),iendsy(3),ibegsz(3),iendsz(3)
                               
real   (kind=8), intent(out) :: F(0:(iend(1)-ibeg(1)+1)-1, &
  0:(iend(2)-ibeg(2)+1)*(iend(3)-ibeg(3)+1)-1)
integer(kind=4) :: mx,my,mz,ngx,ngy,ngz,ex,ey,ez
integer(kind=4) :: kx,ky,kz,ax,ay,az,bx,by,bz,d
integer(kind=4) :: Ox(nelem(1)),Oy(nelem(2)),Oz(nelem(3))
real   (kind=8) :: Jx(nelem(1)),Jy(nelem(2)),Jz(nelem(3))
real   (kind=8) :: Wx(p(1)+1),Wy(p(2)+1),Wz(p(3)+1)
real   (kind=8) :: Xx(p(1)+1,nelem(1))
real   (kind=8) :: Xy(p(2)+1,nelem(2))
real   (kind=8) :: Xz(p(3)+1,nelem(3))
real   (kind=8) :: NNx(0:p(1)-1,0:p(1),p(1)+1,nelem(1)), &
                   NNy(0:p(2)-1,0:p(2),p(2)+1,nelem(2)), &
                   NNz(0:p(3)-1,0:p(3),p(3)+1,nelem(3))
real   (kind=8) :: J,W,Uval,ucoeff
real   (kind=8) :: v, rhs
real   (kind=8) :: dux,duy,duz,dvx,dvy,dvz
integer(kind=4) :: nreppx,nreppy,nreppz !# elements per proc along x,y,z
integer(kind=4) :: ind,ind1,ind23,indx,indy,indz
integer(kind=4) :: indbx,indby,indbz
integer(kind=4) :: rx,ry,rz, ix,iy,iz, sx,sy,sz
real   (kind=8) :: resvalue

  d = 0
  mx  = n(1) + p(1) + 1
  ngx = p(1) + 1
  my  = n(2) + p(2) + 1
  ngy = p(2) + 1
  mz  = n(3) + p(3) + 1
  ngz = p(3) + 1

  call BasisData(p(1),mx,Ux,p(1)-1,ngx,nelem(1),Ox,Jx,Wx,Xx,NNx)
  call BasisData(p(2),my,Uy,p(2)-1,ngy,nelem(2),Oy,Jy,Wy,Xy,NNy)
  call BasisData(p(3),mz,Uz,p(3)-1,ngz,nelem(3),Oz,Jz,Wz,Xz,NNz)

#ifdef IPRINT
    write(*,*)PRINTRANK,'ex:',mine(1),maxe(1)
    write(*,*)PRINTRANK,'ey:',mine(2),maxe(2)
    write(*,*)PRINTRANK,'ez:',mine(3),maxe(3)
    write(*,*)PRINTRANK,'ibegx,iendx',ibeg(1),iend(1)
    write(*,*)PRINTRANK,'ibegy,iendy',ibeg(2),iend(2)
    write(*,*)PRINTRANK,'ibegz,iendz',ibeg(3),iend(3)
#endif

  F = 0

  do ex = mine(1), maxe(1)
  do ey = mine(2), maxe(2)
  do ez = mine(3), maxe(3)
    J = Jx(ex)*Jy(ey)*Jz(ez)
    do kx = 1,ngx
    do ky = 1,ngy
    do kz = 1,ngz
      W = Wx(kx)*Wy(ky)*Wz(kz)  
      do ax = 0,p(1)
      do ay = 0,p(2)
      do az = 0,p(3)
        ind = (Ox(ex)+ax) + (Oy(ey)+ay)*(n(1)+1) + (Oz(ez)+az)*(n(2)+1)*(n(1)+1)
        call global2local(ind,n,indx,indy,indz)

        if (indx < ibeg(1)-1 .or. indx > iend(1)-1) cycle
        if (indy < ibeg(2)-1 .or. indy > iend(2)-1) cycle
        if (indz < ibeg(3)-1 .or. indz > iend(3)-1) cycle

        ind1 = indx-ibeg(1)+1
        ind23 = (indy-ibeg(2)+1) + (indz-ibeg(3)+1)*(iend(2)-ibeg(2)+1)

        Uval = 0
        dux = 0
        duy = 0
        duz = 0
        do bx = 0,p(1)
        do by = 0,p(2)
        do bz = 0,p(3)
          ind = (Ox(ex)+bx) + (Oy(ey)+by)*(n(1)+1) + (Oz(ez)+bz)*(n(2)+1)*(n(1)+1)
          call global2local(ind,n,indbx,indby,indbz)

          rx = 2
          ry = 2
          rz = 2
          if (indbx < ibeg(1)-1) rx = 1
          if (indbx > iend(1)-1) rx = 3
          if (indby < ibeg(2)-1) ry = 1
          if (indby > iend(2)-1) ry = 3
          if (indbz < ibeg(3)-1) rz = 1
          if (indbz > iend(3)-1) rz = 3

          ix = indbx - ibegsx(rx) + 1
          iy = indby - ibegsy(ry) + 1
          iz = indbz - ibegsz(rz) + 1
          sx = iendsx(rx) - ibegsx(rx) + 1
          sy = iendsy(ry) - ibegsy(ry) + 1
          sz = iendsz(rz) - ibegsz(rz) + 1
          ind = ix + sx * (iy + sy * iz)

#ifdef IDEBUG
          if (ind < 0 .or. ind > nrcpp(3)*nrcpp(1)*nrcpp(2)-1) then
            write(ERROR_UNIT,*)PRINTRANK,'Oh crap',ix,iy,iz
            write(ERROR_UNIT,*)PRINTRANK,'r',rx,ry,rz
            write(ERROR_UNIT,*)PRINTRANK,'x',ibeg(1),iend(1)
            write(ERROR_UNIT,*)PRINTRANK,'y',ibeg(2),iend(2)
            write(ERROR_UNIT,*)PRINTRANK,'z',ibeg(3),iend(3)
            write(ERROR_UNIT,*)PRINTRANK,'idxs',indbx,indby,indbz
            write(ERROR_UNIT,*)PRINTRANK,'sizes=',sx,sy,sz
            write(ERROR_UNIT,*)PRINTRANK,'begsx=',ibegsx
            write(ERROR_UNIT,*)PRINTRANK,'endsx=',iendsx
            write(ERROR_UNIT,*)PRINTRANK,'begsy=',ibegsy
            write(ERROR_UNIT,*)PRINTRANK,'endsy=',iendsy
            write(ERROR_UNIT,*)PRINTRANK,'begsz=',ibegsz
            write(ERROR_UNIT,*)PRINTRANK,'endsz=',iendsz
          endif
#endif
          
          Ucoeff = R(ind,rx,ry,rz)
          v   = NNx(0,bx,kx,ex) * NNy(0,by,ky,ey) * NNz(0,bz,kz,ez)
          dvx = NNx(1,bx,kx,ex) * NNy(0,by,ky,ey) * NNz(0,bz,kz,ez) 
          dvy = NNx(0,bx,kx,ex) * NNy(1,by,ky,ey) * NNz(0,bz,kz,ez) 
          dvz = NNx(0,bx,kx,ex) * NNy(0,by,ky,ey) * NNz(1,bz,kz,ez) 

          Uval = Uval + Ucoeff * v
          dux = dux + Ucoeff * dvx
          duy = duy + Ucoeff * dvy
          duz = duz + Ucoeff * dvz
        enddo
        enddo
        enddo
          call RHS_fun ( &
               [Xx(kx,ex),Xy(ky,ey),Xz(kz,ez)], &
               [kx,ky,kz], &
               [ex,ey,ez], &
               p, &
               nelem, &
               [ax,ay,az], &
               [bx,by,bz], &
               [dux,duy,duz], &
               [ibeg(1),ibeg(2),ibeg(3)], &
               [iend(1),iend(2),iend(3)], &
               [mine(1),mine(2),mine(3)], &
               [maxe(1),maxe(2),maxe(3)], &
               NNx,NNy,NNz, &
               Uval,J,W,resvalue)
   
          F(ind1,ind23) = F(ind1,ind23) + resvalue
   
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
logical function IndexInRange(ind,ibeg,iend)
implicit none
integer(kind=4), intent(in), dimension(3) :: ind
integer(kind=4), intent(in), dimension(3) :: ibeg,iend

IndexInRange = .true.
if (ind(1) < ibeg(1)-1 .or. ind(1) > iend(1)-1) IndexInRange = .false.
if (ind(2) < ibeg(2)-1 .or. ind(2) > iend(2)-1) IndexInRange = .false.
if (ind(3) < ibeg(3)-1 .or. ind(3) > iend(3)-1) IndexInRange = .false.

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
subroutine global2local(ind,n,x,y,z)
implicit none
integer(kind=4), intent(in)  :: ind
integer(kind=4), intent(in), dimension(3)  :: n
integer(kind=4), intent(out) :: x,y,z
integer(kind=4) :: tmp

z = ind / ((n(1)+1)*(n(2)+1))
tmp = ind - z*(n(1)+1)*(n(2)+1)
y = tmp / (n(1)+1)
x = tmp - y*(n(1)+1)

end subroutine




! -------------------------------------------------------------------
! Calculates mass matrix M
! -------------------------------------------------------------------
subroutine ComputeMassMatrix(KL,KU,U,p,n,nelem,M)
use parallelism, ONLY : PRINTRANK
implicit none
integer(kind=4), intent(in)  :: KL,KU
integer(kind=4), intent(in)  :: n, p, nelem
real   (kind=8), intent(in)  :: U(0:n+p+1)
real   (kind=8), intent(out) :: M(0:(2*KL+KU),0:n)
integer :: i

  call Form1DMassMatrix(KL,KU,U,p,n,nelem,M)
#ifdef IPRINT
    write(*,*)PRINTRANK,'M'
    do i = 1,2*KL+KU+1
      write(*,*)PRINTRANK,M(i,1:n+1)
    enddo
#endif

end subroutine


end module

