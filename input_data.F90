module input_data

use math
use basis

implicit none


integer, parameter :: pN = 80
real   (kind=8) :: s = 0.05d0
real   (kind=8) :: px(pN), py(pN), pz(pN)

integer, parameter :: lN = 30
real   (kind=8) :: lx(lN), ly(lN), lz(lN)

! Curve number and lentgh
integer, parameter :: cN = 20, cL=16
real   (kind=8) :: cx(cN*cL), cy(cN*cL), cz(cN*cL)

real   (kind=8) :: mi = 10.d0

contains

subroutine InitInputData()
integer       :: i,j
real (kind=8) :: t(3), x(3), dx(3), ddx(3)
real (kind=8) :: f = 0.1d0, step = 0.05d0

  call random_number(px)
  call random_number(py)
  call random_number(pz)

  call random_number(lx)
  call random_number(ly)
  call random_number(lz)

  do i = 0,cN-1
    call random_number(x)
    cx(i*cL + 1) = x(1)
    cy(i*cL + 1) = x(2)
    cz(i*cL + 1) = x(3)
    call random_number(dx)

    do j = 2,cL
      call random_number(t)
      call random_number(ddx)
      ddx = 2 * (ddx - 0.5d0) - 0.4d0 * dx * sum(dx**2)
      dx = dx + 0.4d0 * ddx
      cx(i*cL + j) = cx(i*cL + j - 1) + step * dx(1)
      cy(i*cL + j) = cy(i*cL + j - 1) + step * dx(2)
      cz(i*cL + j) = cz(i*cL + j - 1) + step * dx(3)
    end do
  end do

end subroutine InitInputData

function dist_from_segment(x,y,z,ax,ay,az,bx,by,bz) result (d)
real   (kind=8) :: x,y,z,ax,ay,az,bx,by,bz
real   (kind=8) :: dx,dy,dz,cx,cy,cz,xx,yy,zz
real   (kind=8) :: dot, len2, proj, d

 dx = bx - ax
 dy = by - ay
 dz = bz - az

 cx = x - ax
 cy = y - ay
 cz = z - az

 dot = dx * cx + dy * cy + dz * cz
 len2 = dx ** 2 + dy **2 + dz ** 2
 proj = dot / len2

 if (proj < 0) then
   xx = ax
   yy = ay
   zz = az
 else if (proj > 1) then
   xx = bx
   yy = by
   zz = bz
 else
   xx = ax + proj * dx
   yy = ay + proj * dy
   zz = az + proj * dz
 end if

 d = (x - xx)**2 + (y - yy)**2 + (z - zz)**2

end function


function dist_from_curves(x, y, z) result (fval)
real   (kind=8) :: x, y, z
real   (kind=8) :: ax,ay,az,bx,by,bz,fval
integer         :: i, j

 fval = 1e3
 do i = 0,cN-1
   do j = 2,cL
     ax = cx(i*cL + j - 1)
     bx = cx(i*cL + j)
     ay = cy(i*cL + j - 1)
     by = cy(i*cL + j)
     az = cz(i*cL + j - 1)
     bz = cz(i*cL + j)

     fval = min(fval, dist_from_segment(x,y,z,ax,ay,az,bx,by,bz))
   end do
 end do

end function

! Force function
! x, y, z - point in space
function fvalue(x, y, z) result (fval)
real   (kind=8) :: x, y, z
real   (kind=8) :: fval

  fval = 0.d0 ! 1.d0 + sin(2*PI*x) * sin(2*PI*y) * sin(2*PI*z)

end function


function bump3d(x, y, z, x0, y0, z0, sx, sy, sz) result (val)
real   (kind=8) :: x, y, z, x0, y0, z0, sx, sy, sz
real   (kind=8) :: val

  val = bump((x-x0)/sx) * bump((y-y0)/sy) * bump((z-z0)/sz)

end function


function points(x, y, z, sx, sy, sz, n, XX, YY, ZZ) result (val)
real   (kind=8) :: x, y, z, sx, sy, sz
integer :: n
real   (kind=8) :: val, XX(n), YY(n), ZZ(n)
integer :: i

  val = 0

  do i = 1,n
    val = val+bump3d(x,y,z,XX(i),YY(i),ZZ(i),sx,sy,sz)
  end do
  !val = log(1+3*val) / log(4.d0)

end function


function kq(x, y, z) result (val)
real   (kind=8) :: x, y, z
real   (kind=8) :: val

  val = points(x, y, z, s, s, s, pN, px, py, pz)
  val = val + points(x,y,z,0.04d0,0.4d0,0.04d0,lN,lx,ly,lz)
  !val = points(x,y,z,0.04d0,0.04d0,0.04d0,cN*cL,cx,cy,cz)
  val = val + bump(20.d0 * sqrt(dist_from_curves(x,y,z)))
  val = min(1.d0, log(1+3*val) / log(4.d0))

end function


function bq(x, y, z, u) result (val)
real   (kind=8) :: x, y, z, u
real   (kind=8) :: val

  val = exp(mi * u)

end function


function h(x, y, z) result (val)
real   (kind=8) :: x, y, z
real   (kind=8) :: val

  val = 1 + sin(2*PI*x) * sin(2*PI*y) * sin(2*PI*z)

end function


! Initial state of the system - u(0)
function initial_state(x, y, z) result (val)
real   (kind=8) :: x, y, z
real   (kind=8) :: val

  val = 0.1d0 * kq(x, y, z)

end function



subroutine CacheKqValues(Ux,px,nx,minex,maxex,nelemx,Uy,py,ny,miney,maxey,nelemy,Uz,pz,nz,minez,maxez,nelemz,Kq_vals)
implicit none
integer(kind=4), intent(in) :: nx, px, minex, maxex, nelemx
integer(kind=4), intent(in) :: ny, py, miney, maxey, nelemy
integer(kind=4), intent(in) :: nz, pz, minez, maxez, nelemz
real   (kind=8), intent(in) :: Ux(0:nx+px+1)
real   (kind=8), intent(in) :: Uy(0:ny+py+1)
real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
real   (kind=8), intent(out) :: Kq_vals(px+1,py+1,pz+1,maxex-minex+1,maxey-miney+1,maxez-minez+1)
integer(kind=4) :: mx,my,mz,ngx,ngy,ngz,ex,ey,ez
integer(kind=4) :: kx,ky,kz,ax,ay,az,d
integer(kind=4) :: Ox(nelemx),Oy(nelemy),Oz(nelemz)
real   (kind=8) :: Jx(nelemx),Jy(nelemy),Jz(nelemz)
real   (kind=8) :: Wx(px+1),Wy(py+1),Wz(pz+1)
real   (kind=8) :: Xx(px+1,nelemx)
real   (kind=8) :: Xy(py+1,nelemy)
real   (kind=8) :: Xz(pz+1,nelemz)
real   (kind=8) :: NNx(0:2,0:px,px+1,nelemx), NNy(0:2,0:py,py+1,nelemy), NNz(0:2,0:pz,pz+1,nelemz)
d=2
mx  = nx+px+1
ngx = px+1
my  = ny+py+1
ngy = py+1
mz  = nz+pz+1
ngz = pz+1

call BasisData(px,mx,Ux,d,ngx,nelemx,Ox,Jx,Wx,Xx,NNx)
call BasisData(py,my,Uy,d,ngy,nelemy,Oy,Jy,Wy,Xy,NNy)
call BasisData(pz,mz,Uz,d,ngz,nelemz,Oz,Jz,Wz,Xz,NNz)

do ex = minex,maxex
do ey = miney,maxey
do ez = minez,maxez
  do kx = 1,ngx
  do ky = 1,ngy
  do kz = 1,ngz
    Kq_vals(kx,ky,kz,ex-minex+1,ey-miney+1,ez-minez+1) = kq(Xx(kx,ex),Xy(ky,ey),Xz(kz,ez))
  end do
  end do
  end do
end do
end do
end do
end subroutine


end module input_data
