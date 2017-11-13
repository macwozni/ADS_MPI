module input_data

implicit none

! Curve number and lentgh
integer(kind=4), parameter :: cN = 30, cL=16
real   (kind=8), parameter :: radius = 0.15, pumping_strength = 1, draining_strength = 1

!!! krzywe dane przez segmenty
real   (kind=8) :: cx(cN*cL), cy(cN*cL), cz(cN*cL)
real   (kind=8) :: mi = 10.d0
real   (kind=8) :: GROUND = 0.2
real   (kind=8), parameter :: Kqmin = 1.d0, Kqmax = 1000.d0
integer(kind=4) :: npumps, ndrains
real   (kind=8), allocatable, dimension(:,:) :: pumps, drains

! Buffer for values of permeability function
real (kind=8), allocatable :: Kqvals(:,:,:,:,:,:)

! Current time
real   (kind=8) :: t

! Time and timestep
real (kind=8) :: Dt

! Number of iterations
integer :: steps

! Statistics computed during the simulation
real (kind=8) ::  pollution = 0

! order of approximations
integer(kind=4) :: ORDER

! number of elements in one dimension
integer(kind=4) :: SIZE

integer(kind=4) :: procx,procy,procz

real (kind=8) :: drained = 0
real (kind=8) :: l2norm,fullnorm


contains


! -------------------------------------------------------------------
! Sets values of parameters (order and size)
! -------------------------------------------------------------------
subroutine InitializeParameters
implicit none
character(100) :: input

  ! ./l2 <size> <procx> <procy> <procz> <nsteps> <dt>
  ORDER = 2

  call getarg(1,input)
  read(input,*) SIZE
  call getarg(2,input)
  read(input,*) procx
  call getarg(3,input)
  read(input,*) procy
  call getarg(4,input)
  read(input,*) procz
  call getarg(5,input)
  read(input,*) steps
  call getarg(6,input)
  read(input,*) Dt

end subroutine






subroutine ComputeResults()
use parallelism, ONLY : MYRANK
implicit none
include "mpif.h"
real   (kind=8) :: fulldrained
integer(kind=4) :: ierr


  call MPI_Reduce(drained, fulldrained, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (MYRANK == 0) then
    write(*,*) fulldrained
  endif

end subroutine




subroutine InitPumps()
implicit none
character(100)  :: input
integer(kind=4) :: i, arg = 7 ! First argument after "technical" ones

call getarg(arg,input)
read(input,*) npumps
arg = arg + 1
allocate(pumps(3,npumps))

do i = 1,npumps
  call getarg(arg, input)
  read(input,*) pumps(1,i)
  call getarg(arg + 1, input)
  read(input,*) pumps(2,i)
  call getarg(arg + 2, input)
  read(input,*) pumps(3,i)
  arg = arg + 3
enddo

call getarg(arg,input)
read(input,*) ndrains
arg = arg + 1
allocate(drains(3,ndrains))

do i = 1,ndrains
  call getarg(arg, input)
  read(input,*) drains(1,i)
  call getarg(arg + 1, input)
  read(input,*) drains(2,i)
  call getarg(arg + 2, input)
  read(input,*) drains(3,i)
  arg = arg + 3
enddo

end subroutine


subroutine InitInputData()
implicit none
integer(kind=4) :: i,j
real   (kind=8) :: t(3), x(3), dx(3), ddx(3)
real   (kind=8) :: f = 0.1d0, step = 0.05d0

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

call InitPumps()

end subroutine



!!!! liczy odleglosc punktu od odcinka
!!!!! przeniesc gdzies indziej
function dist_from_segment(point,ibeg,iend) result (d)
implicit none
real   (kind=8), intent(in), dimension(3) :: point,ibeg,iend
real   (kind=8) :: dx,dy,dz,cx,cy,cz,xx,yy,zz
real   (kind=8) :: dot, len2, proj, d

dx = iend(1) - ibeg(1)
dy = iend(2) - ibeg(2)
dz = iend(3) - ibeg(3)

cx = point(1) - ibeg(1)
cy = point(2) - ibeg(2)
cz = point(3) - ibeg(3)

dot = dx * cx + dy * cy + dz * cz
len2 = dx ** 2 + dy **2 + dz ** 2
proj = dot / len2

if (proj < 0) then
  xx = ibeg(1)
  yy = ibeg(2)
  zz = ibeg(3)
else if (proj > 1) then
  xx = iend(1)
  yy = iend(2)
  zz = iend(3)
else
  xx = ibeg(1) + proj * dx
  yy = ibeg(2) + proj * dy
  zz = ibeg(3) + proj * dz
end if

d = (point(1) - xx)**2 + (point(2) - yy)**2 + (point(3) - zz)**2

end function



!!!!! liczy odleglosc punktu od krzywych
!!!!! krzywe jako zmienne globalne
!!!!! przeniesc do gdziesc indziej
function dist_from_curves(point,cx,cy,cz,cN,cL) result (fval)
implicit none
real   (kind=8), intent(in), dimension(3) :: point
integer(kind=4), intent(in) :: cN, cL
real   (kind=8), intent(in) :: cx(cN*cL), cy(cN*cL), cz(cN*cL)
real   (kind=8) :: ax,ay,az,bx,by,bz,fval
integer(kind=4) :: i, j

fval = 1e3
do i = 0,cN-1
  do j = 2,cL
    ax = cx(i*cL + j - 1)
    bx = cx(i*cL + j)
    ay = cy(i*cL + j - 1)
    by = cy(i*cL + j)
    az = cz(i*cL + j - 1)
    bz = cz(i*cL + j)

    fval = min(fval, dist_from_segment(point,[ax,ay,az],[bx,by,bz]))
  end do
end do

end function



! Pumping
! x, y, z - point in space
function pumping(x, y, z) result (fval)
use math, ONLY : falloff
implicit none
real   (kind=8) :: x, y, z
real   (kind=8) :: fval
integer(kind=4) :: i

fval = 0.d0
do i = 1,npumps
  fval = fval + pumping_strength * falloff(0.d0, radius, norm2(pumps(:,i) - [x, y, z]))
enddo

end function




! Draining
! x, y, z - point in space
function draining(u, x, y, z) result (fval)
use math, ONLY : falloff
implicit none
real   (kind=8) :: u, x, y, z
real   (kind=8) :: fval
integer :: i

fval = 0.d0
do i = 1,ndrains
  fval = fval + draining_strength * falloff(0.d0, radius, norm2(drains(:,i) - [x, y, z]))
enddo
fval = fval * u

end function


function kq(x, y, z) result (val)
use math, ONLY : falloff,lerp
implicit none
real   (kind=8),intent(in) :: x, y, z
real   (kind=8) :: val, dist

dist = sqrt(dist_from_curves([x,y,z],cx,cy,cz,cN,cL))
val = lerp(falloff(0.d0, 0.06d0, dist), Kqmin, Kqmax)

end function


function bq(x, y, z, u) result (val)
real   (kind=8) :: x, y, z, u
real   (kind=8) :: val

val = exp(mi * u)

end function



! Initial state of the system - u(0)
function initial_state(x, y, z) result (val)
use math, ONLY : falloff,bump3d,lerp
implicit none
real   (kind=8), intent(in) :: x, y, z
real   (kind=8) :: dist, val

dist = sqrt(dist_from_curves([x,y,z],cx,cy,cz,cN,cL))
val = 0.1d0 * lerp(falloff(0.d0, 0.1d0, dist), 0.d0, 1.d0) * bump3d(0.2d0, 0.6d0, x, y, z)

end function



subroutine CacheKqValues(Ux,px,nx,minex,maxex,nelemx,Uy,py,ny,miney,maxey,nelemy,Uz,pz,nz,minez,maxez,nelemz,Kq_vals)
use basis, ONLY : BasisData
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



subroutine PrecomputeKq(ads)
use ADSS, ONLY : ADS_setup
implicit none
type   (ADS_setup) :: ads

  call CacheKqValues                                        &
      (ads%Ux,ads%p(1),ads%n(1),ads%mine(1),ads%maxe(1),ads%nelem(1), &
       ads%Uy,ads%p(2),ads%n(2),ads%mine(2),ads%maxe(2),ads%nelem(2), &
       ads%Uz,ads%p(3),ads%n(3),ads%mine(3),ads%maxe(3),ads%nelem(3), &
       Kqvals)
 
end subroutine


end module