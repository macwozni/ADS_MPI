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
   Ux,px,nx,nelemx,nrcppx,     &
   Uy,py,ny,nelemy,nrcppy,     &
   Uz,pz,nz,nelemz,nrcppz,     &
   ibegx,iendx,                &
   ibegy,iendy,                &
   ibegz,iendz,                &
   ibegsx,iendsx,              &
   ibegsy,iendsy,              &
   ibegsz,iendsz,              &
   minex,maxex,                &
   miney,maxey,                &
   minez,maxez,                &
   R,F,RHS_fun)
use parallelism, ONLY : PRINTRANK
USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT ! access computing environment
use basis, ONLY : BasisData
implicit none
interface
  subroutine RHS_fun( &
         Xx,Xy,Xz, &
         kx,ky,kz, &
         ex,ey,ez, &
         px,py,pz, &
         nelemx,nelemy,nelemz, &
         ax,ay,az, &
         bx,by,bz, &
         NNx,NNy,NNz, &
         dux,duy,duz, &
         ibegx,ibegy,ibegz, &
         iendx,iendy,iendz, &
         minex,miney,minez, &
         maxex,maxey,maxez, &
         Uval,J,W,F)
      implicit none
      integer(kind=4), intent(in)  :: px,py,pz
      real   (kind=8), intent(in)  :: Xx(px+1,nelemx)
      real   (kind=8), intent(in)  :: Xy(py+1,nelemy)
      real   (kind=8), intent(in)  :: Xz(pz+1,nelemz)
      integer(kind=4), intent(in)  :: kx,ky,kz
      integer(kind=4), intent(in)  :: ex,ey,ez
      integer(kind=4), intent(in)  :: nelemx,nelemy,nelemz
      real   (kind=8), intent(in)  :: Uval
      integer(kind=4), intent(in)  :: ibegx,ibegy,ibegz
      integer(kind=4), intent(in)  :: iendx,iendy,iendz
      integer(kind=4), intent(in)  :: maxex,maxey,maxez
      integer(kind=4), intent(in)  :: minex,miney,minez
      integer(kind=4), intent(in)  :: ax,ay,az
      integer(kind=4), intent(in)  :: bx,by,bz
      real   (kind=8), intent(in)  :: dux,duy,duz
      real   (kind=8), intent(in)  :: J,W
      real   (kind=8), intent(in)  :: NNx(0:px-1,0:px,px+1,nelemx), &
                         NNy(0:py-1,0:py,py+1,nelemy), &
                         NNz(0:pz-1,0:pz,pz+1,nelemz)
      real   (kind=8), intent(out) :: F
  end subroutine
end interface
integer(kind=4), intent(in)  :: nx, px, nelemx, nrcppx
integer(kind=4), intent(in)  :: ny, py, nelemy, nrcppy
integer(kind=4), intent(in)  :: nz, pz, nelemz, nrcppz
integer(kind=4), intent(in)  :: minex,maxex,miney,maxey,minez,maxez
real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
real   (kind=8), intent(in)  :: Uz(0:nz+pz+1)
real   (kind=8), intent(in)  :: R(0:nrcppz*nrcppx*nrcppy-1,3,3,3)
integer(kind=4), intent(in) :: ibegsx(3),iendsx(3),ibegsy(3),iendsy(3),ibegsz(3),iendsz(3)
integer(kind=4), intent(in)  :: ibegx,ibegy,ibegz
integer(kind=4), intent(in)  :: iendx,iendy,iendz
                               
real   (kind=8), intent(out) :: F(0:(iendx-ibegx+1)-1, &
  0:(iendy-ibegy+1)*(iendz-ibegz+1)-1)
integer(kind=4) :: mx,my,mz,ngx,ngy,ngz,ex,ey,ez
integer(kind=4) :: kx,ky,kz,ax,ay,az,bx,by,bz,d
integer(kind=4) :: Ox(nelemx),Oy(nelemy),Oz(nelemz)
real   (kind=8) :: Jx(nelemx),Jy(nelemy),Jz(nelemz)
real   (kind=8) :: Wx(px+1),Wy(py+1),Wz(pz+1)
real   (kind=8) :: Xx(px+1,nelemx)
real   (kind=8) :: Xy(py+1,nelemy)
real   (kind=8) :: Xz(pz+1,nelemz)
real   (kind=8) :: NNx(0:px-1,0:px,px+1,nelemx), &
                   NNy(0:py-1,0:py,py+1,nelemy), &
                   NNz(0:pz-1,0:pz,pz+1,nelemz)
real   (kind=8) :: J,W,Uval,ucoeff
real   (kind=8) :: v, rhs
real   (kind=8) :: dux,duy,duz,dvx,dvy,dvz
integer(kind=4) :: nreppx,nreppy,nreppz !# elements per proc along x,y,z
integer(kind=4) :: ind,ind1,ind23,indx,indy,indz
integer(kind=4) :: indbx,indby,indbz
integer(kind=4) :: rx,ry,rz, ix,iy,iz, sx,sy,sz
real   (kind=8) :: resvalue

  d = 0
  mx  = nx + px + 1
  ngx = px + 1
  my  = ny + py + 1
  ngy = py + 1
  mz  = nz + pz + 1
  ngz = pz + 1

  call BasisData(px,mx,Ux,px-1,ngx,nelemx,Ox,Jx,Wx,Xx,NNx)
  call BasisData(py,my,Uy,py-1,ngy,nelemy,Oy,Jy,Wy,Xy,NNy)
  call BasisData(pz,mz,Uz,pz-1,ngz,nelemz,Oz,Jz,Wz,Xz,NNz)

#ifdef IPRINT
    write(*,*)PRINTRANK,'ex:',minex,maxex
    write(*,*)PRINTRANK,'ey:',miney,maxey
    write(*,*)PRINTRANK,'ez:',minez,maxez
    write(*,*)PRINTRANK,'ibegx,iendx',ibegx,iendx
    write(*,*)PRINTRANK,'ibegy,iendy',ibegy,iendy
    write(*,*)PRINTRANK,'ibegz,iendz',ibegz,iendz
#endif

  F = 0

  do ex = minex, maxex
  do ey = miney, maxey
  do ez = minez, maxez
    J = Jx(ex)*Jy(ey)*Jz(ez)
    do kx = 1,ngx
    do ky = 1,ngy
    do kz = 1,ngz
      W = Wx(kx)*Wy(ky)*Wz(kz)  
      do ax = 0,px
      do ay = 0,py
      do az = 0,pz
        ind = (Ox(ex)+ax) + (Oy(ey)+ay)*(nx+1) + (Oz(ez)+az)*(ny+1)*(nx+1)
        call global2local(ind,nx,ny,nz,indx,indy,indz)

        if (indx < ibegx-1 .or. indx > iendx-1) cycle
        if (indy < ibegy-1 .or. indy > iendy-1) cycle
        if (indz < ibegz-1 .or. indz > iendz-1) cycle

        ind1 = indx-ibegx+1
        ind23 = (indy-ibegy+1) + (indz-ibegz+1)*(iendy-ibegy+1)

        Uval = 0
        dux = 0
        duy = 0
        duz = 0
        do bx = 0,px
        do by = 0,py
        do bz = 0,pz
          ind = (Ox(ex)+bx) + (Oy(ey)+by)*(nx+1) + (Oz(ez)+bz)*(ny+1)*(nx+1)
          call global2local(ind,nx,ny,nz,indbx,indby,indbz)

          rx = 2
          ry = 2
          rz = 2
          if (indbx < ibegx-1) rx = 1
          if (indbx > iendx-1) rx = 3
          if (indby < ibegy-1) ry = 1
          if (indby > iendy-1) ry = 3
          if (indbz < ibegz-1) rz = 1
          if (indbz > iendz-1) rz = 3

          ix = indbx - ibegsx(rx) + 1
          iy = indby - ibegsy(ry) + 1
          iz = indbz - ibegsz(rz) + 1
          sx = iendsx(rx) - ibegsx(rx) + 1
          sy = iendsy(ry) - ibegsy(ry) + 1
          sz = iendsz(rz) - ibegsz(rz) + 1
          ind = ix + sx * (iy + sy * iz)

!!!!! debug
          if (ind < 0 .or. ind > nrcppz*nrcppx*nrcppy-1) then
            write(ERROR_UNIT,*)PRINTRANK,'Oh crap',ix,iy,iz
            write(ERROR_UNIT,*)PRINTRANK,'r',rx,ry,rz
            write(ERROR_UNIT,*)PRINTRANK,'x',ibegx,iendx
            write(ERROR_UNIT,*)PRINTRANK,'y',ibegy,iendy
            write(ERROR_UNIT,*)PRINTRANK,'z',ibegz,iendz
            write(ERROR_UNIT,*)PRINTRANK,'idxs',indbx,indby,indbz
            write(ERROR_UNIT,*)PRINTRANK,'sizes=',sx,sy,sz
            write(ERROR_UNIT,*)PRINTRANK,'begsx=',ibegsx
            write(ERROR_UNIT,*)PRINTRANK,'endsx=',iendsx
            write(ERROR_UNIT,*)PRINTRANK,'begsy=',ibegsy
            write(ERROR_UNIT,*)PRINTRANK,'endsy=',iendsy
            write(ERROR_UNIT,*)PRINTRANK,'begsz=',ibegsz
            write(ERROR_UNIT,*)PRINTRANK,'endsz=',iendsz
          endif
          
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
               Xx,Xy,Xz, &
               kx,ky,kz, &
               ex,ey,ez, &
               px,py,pz, &
               nelemx,nelemy,nelemz, &
               ax,ay,az, &
               bx,by,bz, &
               NNx,NNy,NNz, &
               dux,duy,duz, &
               ibegx,ibegy,ibegz, &
               iendx,iendy,iendz, &
               minex,miney,minez, &
               maxex,maxey,maxez, &
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
! ibegx,iendx      - x range
! ibegy,iendy      - y range
! ibegz,iendz      - z range
! -------------------------------------------------------------------
logical function IndexInRange(indx,indy,indz,ibegx,iendx,ibegy,iendy,ibegz,iendz)
implicit none
integer(kind=4) :: indx,indy,indz
integer(kind=4) :: ibegx,iendx,ibegy,iendy,ibegz,iendz

IndexInRange = .true.
if (indx < ibegx-1 .or. indx > iendx-1) IndexInRange = .false.
if (indy < ibegy-1 .or. indy > iendy-1) IndexInRange = .false.
if (indz < ibegz-1 .or. indz > iendz-1) IndexInRange = .false.

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
subroutine global2local(ind,nx,ny,nz,x,y,z)
implicit none
integer(kind=4), intent(in)  :: ind
integer(kind=4), intent(in)  :: nx,ny,nz
integer(kind=4), intent(out) :: x,y,z
integer(kind=4) :: tmp

z = ind / ((nx+1)*(ny+1))
tmp = ind - z*(nx+1)*(ny+1)
y = tmp / (nx+1)
x = tmp - y*(nx+1)

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

