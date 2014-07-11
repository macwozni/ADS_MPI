module projection_engine

use gauss
use basis
use parallelism
use reorderRHS


! order of approximations
integer :: ORDER
! number of elements in one dimension
integer :: SIZE

contains

subroutine initialize_parameters
  ORDER = 2
  SIZE = 2
end subroutine initialize_parameters

 
!!$ KL = number of subdiagonals, KU = number of superdiagonals
!!$ d is dimension
!!$ U is the knot vector
!!$ p is the polynomial order
!!$ n is the index of the last control point
!!$ nelem you get from running CountSpans
!!$ M is the dense matrix
subroutine Form1DMassMatrix(KL,KU,U,p,n,nelem,M)
implicit none
integer :: KL,KU
integer(kind=4), intent(in) :: n, p, nelem
real   (kind=8), intent(in) :: U(0:n+p+1)
! M is a band storage matrix A(n,n)->M(2*KL+KU+1,n)
double precision, intent(out) :: M(0:(2*KL+KU),0:n)
! double precision, intent(out) :: M(0:n,0:n)
integer(kind=4) :: mm,ng,e,k,a,b
integer(kind=4) :: O(nelem)
real   (kind=8) :: J(nelem)
real   (kind=8) :: W(p+1)
real   (kind=8) :: X(p+1,nelem)
! real   (kind=8) :: NN(0:d,0:p,p+1,nelem)
real   (kind=8) :: NN(0:0,0:p,p+1,nelem)
integer :: d
integer :: nrepp !# elements per processor
integer :: iprint

  iprint=0
  ! if(MYRANK.eq.0)iprint=1

  mm = n+p+1
  ng = p+1
  d = 0
  call BasisData(p,mm,U,d,ng,nelem,O,J,W,X,NN) 
  M = 0
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
          ! J(e) jacobian ? for element e

          ! M is a band storage matrix A(i,j)->M(KL+KU+1+i-j,j)
          !  write(*,*)'i,j',O(e)+a,O(e)+b,'->',
          !    KL+KU+O(e)+a-(O(e)+b),O(e)+b
          M(KL+KU+O(e)+a-(O(e)+b),O(e)+b) = &
            M(KL+KU+O(e)+a-(O(e)+b),O(e)+b) + NN(0,a,k,e)*NN(0,b,k,e)*J(e)*W(k)
          !  M(O(e)+a,O(e)+b) =
          !    M(O(e)+a,O(e)+b) + NN(0,a,k,e)*NN(0,b,k,e)*J(e)*W(k)
          !      NN(0,a,k,e)*NN(0,b,k,e)*J(e)*W(k)
        end do
      end do
    end do
  end do

end subroutine Form1DMassMatrix

!-> DEBUG
subroutine Form1DMassMatrixFULL(U,p,n,nelem,M) 
!!$ d is dimension
!!$ U is the knot vector
!!$ p is the polynomial order
!!$ n is the index of the last control point
!!$ nelem you get from running CountSpans
!!$ M is the dense matrix
implicit none
integer(kind=4), intent(in) :: n, p, nelem
real   (kind=8), intent(in) :: U(0:n+p+1)
double precision, intent(out) :: M(0:n,0:n)
integer(kind=4) :: mm,ng,e,k,a,b
integer(kind=4) :: O(nelem)
real   (kind=8) :: J(nelem)
real   (kind=8) :: W(p+1)
real   (kind=8) :: X(p+1,nelem)
!      real   (kind=8) :: NN(0:d,0:p,p+1,nelem)
real   (kind=8) :: NN(0:0,0:p,p+1,nelem)
integer :: d

  mm  = n+p+1
  ng  = p+1
  d = 0
  call BasisData(p,mm,U,d,ng,nelem,O,J,W,X,NN) 

  M = 0.d0
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
! J(e) jacobian ? for element e
          M(O(e)+a,O(e)+b) = &
             M(O(e)+a,O(e)+b) + NN(0,a,k,e)*NN(0,b,k,e)*J(e)*W(k)
        end do
      end do
    end do
  end do

end subroutine Form1DMassMatrixFULL
!<- DEBUG


! Ux,Uy,Uz - knots vectors
! px,py,py - orders
! nx,ny,nz - number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
! nelemx,nelemy,nelemz - number of elements
! F output rhs (multiple vectors)
subroutine Form3DRHS(          &
   Ux,px,nx,nelemx,            &
   Uy,py,ny,nelemy,            &
   Uz,pz,nz,nelemz,            &
   ibegx,iendx,nrankx,nrpx,    &
   ibegy,iendy,nranky,nrpy,    &
   ibegz,iendz,nrankz,nrpz,F)
implicit none
integer(kind=4), intent(in) :: nx, px, nelemx
integer(kind=4), intent(in) :: ny, py, nelemy
integer(kind=4), intent(in) :: nz, pz, nelemz
real   (kind=8), intent(in) :: Ux(0:nx+px+1)
real   (kind=8), intent(in) :: Uy(0:ny+py+1)
real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
double precision, intent(out) :: F(0:(iendx-ibegx+1)-1, &
  0:(iendy-ibegy+1)*(iendz-ibegz+1)-1)
integer(kind=4) :: mx,my,mz,ngx,ngy,ngz,ex,ey,ez
integer(kind=4) :: kx,ky,kz,ax,ay,az,d
integer(kind=4) :: Ox(nelemx),Oy(nelemy),Oz(nelemz)
real   (kind=8) :: Jx(nelemx),Jy(nelemy),Jz(nelemz)
real   (kind=8) :: Wx(px+1),Wy(py+1),Wz(pz+1)
real   (kind=8) :: Xx(px+1,nelemx)
real   (kind=8) :: Xy(py+1,nelemy)
real   (kind=8) :: Xz(pz+1,nelemz)
real   (kind=8) :: NNx(0:0,0:px,px+1,nelemx), &
                   NNy(0:0,0:py,py+1,nelemy), &
                   NNz(0:0,0:pz,pz+1,nelemz)
real   (kind=8) :: J,W,value
integer, intent(in) :: ibegx,ibegy,ibegz
integer, intent(in) :: iendx,iendy,iendz
integer, intent(in) :: nrankx,nranky,nrankz
integer, intent(in) :: nrpx,nrpy,nrpz
integer :: nreppx,nreppy,nreppz !# elements per proc along x,y,z
integer :: ind,ind1,ind23,ind23a,indx,indy,indz
integer :: iprint

  iprint=0
  ! if(MYRANK.eq.2)iprint=1

  d=0
  mx  = nx+px+1
  ngx = px+1
  my  = ny+py+1
  ngy = py+1
  mz  = nz+pz+1
  ngz = pz+1

  call BasisData(px,mx,Ux,0,ngx,nelemx,Ox,Jx,Wx,Xx,NNx)
  call BasisData(py,my,Uy,0,ngy,nelemy,Oy,Jy,Wy,Xy,NNy)
  call BasisData(pz,mz,Uz,0,ngz,nelemz,Oz,Jz,Wz,Xz,NNz)

  !-> parallel number of elements per processors
  nreppx = nelemx/nrpx
  nreppy = nelemy/nrpy
  nreppz = nelemz/nrpz

  if(iprint == 1)then
    write(*,*)PRINTRANK,'ex:',max(nreppx*nrankx-px+1,1), min(nelemx,nreppx*(nrankx+1)+px)
    write(*,*)PRINTRANK,'ey:',max(nreppy*nranky-py+1,1), min(nelemy,nreppy*(nranky+1)+py)
    write(*,*)PRINTRANK,'ez:',max(nreppz*nrankz-pz+1,1), min(nelemz,nreppz*(nrankz+1)+pz)
    write(*,*)PRINTRANK,'ibegx,iendx',ibegx,iendx
    write(*,*)PRINTRANK,'ibegy,iendy',ibegy,iendy
    write(*,*)PRINTRANK,'ibegz,iendz',ibegz,iendz
  endif

  F = 0.d0

  do ex = max(nreppx*nrankx-px+1,1), min(nelemx,nreppx*(nrankx+1)+px)
  do ey = max(nreppy*nranky-py+1,1), min(nelemy,nreppy*(nranky+1)+py)
  do ez = max(nreppz*nrankz-pz+1,1), min(nelemz,nreppz*(nrankz+1)+pz)
    J = Jx(ex)*Jy(ey)*Jz(ez)
    do kx = 1,ngx
    do ky = 1,ngy
    do kz = 1,ngz
      W = Wx(kx)*Wy(ky)*Wz(kz)
      value = fvalue(Xx(kx,ex),Xy(ky,ey),Xz(kz,ez))
      do ax = 0,px
      do ay = 0,py
      do az = 0,pz
        ind = (Ox(ex)+ax) + (Oy(ey)+ay)*(nx+1) + (Oz(ez)+az)*(ny+1)*(nx+1)
        call global2local(ind,indx,indy,indz,nx,ny,nz)

        if (indx < ibegx-1 .or. indx > iendx-1) cycle
        if (indy < ibegy-1 .or. indy > iendy-1) cycle
        if (indz < ibegz-1 .or. indz > iendz-1) cycle

        ind1 = indx-ibegx+1
        ind23 = (indy-ibegy+1) + (indz-ibegz+1)*(iendy-ibegy+1)

        if (iprint == 1) then
          write(*,*)PRINTRANK,'ind->x,y,z', ind,indx,indy,indz,'->',ind1,ind23
        endif

        if (ind1 < 0 .or. ind1 > (iendx-ibegx)) cycle
        if (ind23 < 0 .or. ind23 > (iendy-ibegy+1)*(iendz-ibegz+1)-1) cycle

        F(ind1,ind23) = F(ind1,ind23) + &
            NNx(0,ax,kx,ex) * NNy(0,ay,ky,ey) * NNz(0,az,kz,ez)*J*W*value

        if (iprint == 1) then
          write(*,*)PRINTRANK, 'ind',ind,'->',ind1,ind23, ex,ey,ez
        endif
      end do
      end do
      end do
    end do
    end do
    end do         
  end do
  end do
  end do
  
end subroutine Form3DRHS

!-> DEBUG
! Ux,Uy,Uz - knots vectors
! px,py,py - orders
! nx,ny,nz - number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
! nelemx,nelemy,nelemz - number of elements
! F output rhs (multiple vectors)
subroutine Form3DRHSFULL(Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz,F)
implicit none
integer(kind=4), intent(in) :: nx, px, nelemx
integer(kind=4), intent(in) :: ny, py, nelemy
integer(kind=4), intent(in) :: nz, pz, nelemz
real   (kind=8), intent(in) :: Ux(0:nx+px+1)
real   (kind=8), intent(in) :: Uy(0:ny+py+1)
real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
double precision, intent(out) :: F(0:(nx+1)*(ny+1)*(nz+1)-1)
integer(kind=4) :: mx,my,mz,ngx,ngy,ngz,ex,ey,ez
integer(kind=4) :: kx,ky,kz,ax,ay,az,d
integer(kind=4) :: Ox(nelemx),Oy(nelemy),Oz(nelemz)
real   (kind=8) :: Jx(nelemx),Jy(nelemy),Jz(nelemz)
real   (kind=8) :: Wx(px+1),Wy(py+1),Wz(pz+1)
real   (kind=8) :: Xx(px+1,nelemx)
real   (kind=8) :: Xy(py+1,nelemy)
real   (kind=8) :: Xz(pz+1,nelemz)
real   (kind=8) :: NNx(0:0,0:px,px+1,nelemx), &
                   NNy(0:0,0:py,py+1,nelemy), &
                   NNz(0:0,0:pz,pz+1,nelemz)
real   (kind=8) :: J,W,value
integer :: ind

  d=0
  mx  = nx+px+1
  ngx = px+1
  my  = ny+py+1
  ngy = py+1
  mz  = nz+pz+1
  ngz = pz+1

  call BasisData(px,mx,Ux,0,ngx,nelemx,Ox,Jx,Wx,Xx,NNx)
  call BasisData(py,my,Uy,0,ngy,nelemy,Oy,Jy,Wy,Xy,NNy)
  call BasisData(pz,mz,Uz,0,ngz,nelemz,Oz,Jz,Wz,Xz,NNz)

  F = 0.d0
  do ex = 1,nelemx
  do ey = 1,nelemy
  do ez = 1,nelemz
    J = Jx(ex)*Jy(ey)*Jz(ez)
    do kx = 1,ngx
    do ky = 1,ngy
    do kz = 1,ngz
      W = Wx(kx)*Wy(ky)*Wz(kz)
      value = fvalue(Xx(kx,ex),Xy(ky,ey),Xz(kz,ez))
      do ax = 0,px
      do ay = 0,py
      do az = 0,pz
        ind = (Oz(ez)+az)*(ny+1)*(nx+1)
        ind = ind + (Oy(ey)+ay)*(nx+1)
        ind = ind + (Ox(ex)+ax)

        if(ind.gt.(nx+1)*(ny+1)*(nz+1))then
          write(*,*)'Form3DRHSFULL:indexing problem'
          stop
        endif

        F(ind) = F(ind) + &
          NNx(0,ax,kx,ex) * NNy(0,ay,ky,ey)* NNz(0,az,kz,ez)*J*W*value
      end do
      end do
      end do

    end do
    end do
    end do         

  end do
  end do
  end do

end subroutine Form3DRHSFULL


function fvalue(x,y,z) result (fval)
implicit none
real   (kind=8) :: x,y,z
real   (kind=8) :: fval

  fval = 1.d0

end function fvalue


subroutine global2local(ind,indx,indy,indz,nx,ny,nz)
integer,intent(in) :: ind
integer :: ind_temp
integer,intent(in) :: nx,ny,nz
integer,intent(out) :: indx,indy,indz

  indz = ind / ((nx+1)*(ny+1))
  ind_temp = ind - indz*(nx+1)*(ny+1)
  indy = ind_temp / (nx+1)
  ind_temp = ind_temp - indy*(nx+1)
  indx = ind_temp 

end subroutine global2local

end module projection_engine

