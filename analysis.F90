module analysis

implicit none

contains


! Ux,Uy,Uz  knots vectors
! px,py,py  orders
! nx,ny,nz  number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
! nelemx,nelemy,nelemz - number of elements
! F coefficients of the function
subroutine NormL2(              &
  Ux,px,nx,nelemx,              &
  Uy,py,ny,nelemy,              &
  Uz,pz,nz,nelemz,              &
  ibegx,iendx,nrankx,nrpx,      &
  ibegy,iendy,nranky,nrpy,      &
  ibegz,iendz,nrankz,nrpz,F)
use basis, ONLY : BasisData
use projection_engine, ONLY : global2local
implicit none
integer(kind=4), intent(in)  :: nx, px, nelemx
integer(kind=4), intent(in)  :: ny, py, nelemy
integer(kind=4), intent(in)  :: nz, pz, nelemz
real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
real   (kind=8), intent(in)  :: Uz(0:nz+pz+1)
real   (kind=8), intent(out) :: F(0:(iendx-ibegx+1)-1, &
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
integer(kind=4), intent(in) :: ibegx,ibegy,ibegz
integer(kind=4), intent(in) :: iendx,iendy,iendz
integer(kind=4), intent(in) :: nrankx,nranky,nrankz
integer(kind=4), intent(in) :: nrpx,nrpy,nrpz
integer(kind=4) :: nreppx,nreppy,nreppz !# elements per proc along x,y,z
integer(kind=4) :: ind,ind1,ind23,ind23a,indx,indy,indz
integer(kind=4) :: iprint

  iprint=0

  d = 0
  mx  = nx+px+1
  ngx = px+1
  my  = ny+py+1
  ngy = py+1
  mz  = nz+pz+1
  ngz = pz+1

  call BasisData(px,mx,Ux,0,ngx,nelemx,Ox,Jx,Wx,Xx,NNx)
  call BasisData(py,my,Uy,0,ngy,nelemy,Oy,Jy,Wy,Xy,NNy)
  call BasisData(pz,mz,Uz,0,ngz,nelemz,Oz,Jz,Wz,Xz,NNz)

  ! parallel number of elements per processors
  nreppx = nelemx/nrpx
  nreppy = nelemy/nrpy
  nreppz = nelemz/nrpz
  F = 0
  do ex = max(nreppx*nrankx-px+1,1),min(nelemx,nreppx*(nrankx+1)+px)
  do ey = max(nreppy*nranky-py+1,1),min(nelemy,nreppy*(nranky+1)+py)
  do ez = max(nreppz*nrankz-pz+1,1),min(nelemz,nreppz*(nrankz+1)+pz)
    J = Jx(ex)*Jy(ey)*Jz(ez)
    do kx = 1,ngx
    do ky = 1,ngy
    do kz = 1,ngz
      W = Wx(kx)*Wy(ky)*Wz(kz)
      do ax = 0,px
      do ay = 0,py
      do az = 0,pz
        d = (Ox(ex)+ax)+(Oy(ey)+ay)*(nx+1)+(Oz(ez)+az)*(ny+1)*(nx+1)
        call global2local(ind,nx,ny,nz,indx,indy,indz)
        if (indx < ibegx-1 .or. indx > iendx-1) cycle
        if (indy < ibegy-1 .or. indy > iendy-1) cycle
        if (indz < ibegz-1 .or. indz > iendz-1) cycle
        ind1 = indx-ibegx+1
        ind23 = (indy-ibegy+1)+(indz-ibegz+1)*(iendy-ibegy+1)

        if (ind1 < 0 .or. ind1 > (iendx-ibegx)) cycle
        if (ind23 < 0 .or. ind23 > (iendy-ibegy+1)*(iendz-ibegz+1)-1) cycle

        ! parallel
        F(ind1,ind23) = F(ind1,ind23) + &
          NNx(0,ax,kx,ex)*NNy(0,ay,ky,ey)*NNz(0,az,kz,ez)*J*W*value
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

end module

