      module projection_engine

      use reorderRHS
      use input_data
      
      contains

      subroutine ComputeProjection
     .   (f,Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz,U0)
      implicit none

      interface
        function f(x, y, z) result (val)
        real (kind=8) :: x,y,z,val
        end function f
      end interface

      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real   (kind=8), intent(in) :: Ux(0:nx+px+1)
      real   (kind=8), intent(in) :: Uy(0:ny+py+1)
      real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
      real   (kind=8), intent(out) :: U0(0:(nx+1)-1,0:(ny+1)*(nz+1)-1)
      integer(kind=4) :: mx,my,mz,ngx,ngy,ngz,ex,ey,ez
      integer(kind=4) :: kx,ky,kz,ax,ay,az,d
      integer(kind=4) :: Ox(nelemx),Oy(nelemy),Oz(nelemz)
      real   (kind=8) :: Jx(nelemx),Jy(nelemy),Jz(nelemz)
      real   (kind=8) :: Wx(px+1),Wy(py+1),Wz(pz+1)
      real   (kind=8) :: Xx(px+1,nelemx)
      real   (kind=8) :: Xy(py+1,nelemy)
      real   (kind=8) :: Xz(pz+1,nelemz)
      real   (kind=8) :: NNx(0:2,0:px,px+1,nelemx),
     .                   NNy(0:2,0:py,py+1,nelemy),
     .                   NNz(0:2,0:pz,pz+1,nelemz)
      real   (kind=8) :: J,W,value,Utvalue
      integer(kind=4) :: bx,by,bz
      integer :: ind,ind1,ind23,inda,ind1a,ind23a
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

      U0 = 0.d0
      do ex = 1,nelemx !loops over elements in x,y,z
         do ey = 1,nelemy
            do ez = 1,nelemz
               J = Jx(ex)*Jy(ey)*Jz(ez)
               do kx = 1,ngx !loops over p+1 Gaussian points (= number of B-spline functions over element)
                  do ky = 1,ngy
                     do kz = 1,ngz
                        W = Wx(kx)*Wy(ky)*Wz(kz)
                        value = f(Xx(kx,ex),Xy(ky,ey),Xz(kz,ez))
                        do ax = 0,px
                           do ay = 0,py
                              do az = 0,pz
                                 ind1 = (Ox(ex)+ax) !along x
                                 ind23 = (Oy(ey)+ay) !along y
                                 ind23 = ind23 + (Oz(ez)+az)*(ny+1) !along z

                                 U0(ind1,ind23) = U0(ind1,ind23) +
     .                               NNx(0,ax,kx,ex)*
     .                               NNy(0,ay,ky,ey)*
     .                               NNz(0,az,kz,ez)*J*W*value
                              end  do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      end subroutine ComputeProjection


      subroutine FillWithInitialState
     .   (Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz,U0)
      implicit none
      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real   (kind=8), intent(in) :: Ux(0:nx+px+1)
      real   (kind=8), intent(in) :: Uy(0:ny+py+1)
      real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
      real   (kind=8), intent(out) :: U0(0:(nx+1)-1,0:(ny+1)*(nz+1)-1)

      call InitInputData()
      call ComputeProjection(initial_state,
     .  Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz,U0)

      end subroutine FillWithInitialState
 
c KL - number of subdiagonals
c KU - number of superdiagonals
c d - derivatives
c U - knot vector
c p - polynomial order
c n - index of the last control point
c nelem - you get from running CountSpans
c M - dense matrix, output
      subroutine Form1DMassMatrix(KL,KU,U,p,n,nelem,M) 
      implicit none
      integer :: KL,KU
      integer(kind=4), intent(in) :: n, p, nelem
      real   (kind=8), intent(in) :: U(0:n+p+1)
c M is a band storage matrix A(n,n)->M(2*KL+KU+1,n)
      double precision, intent(out) :: M(0:(2*KL+KU),0:n)
c      double precision, intent(out) :: M(0:n,0:n)
      integer(kind=4) :: mm,ng,e,k,a,b
      integer(kind=4) :: O(nelem)
      real   (kind=8) :: J(nelem)
      real   (kind=8) :: W(p+1)
      real   (kind=8) :: X(p+1,nelem)
c      real   (kind=8) :: NN(0:d,0:p,p+1,nelem)
      real   (kind=8) :: NN(0:0,0:p,p+1,nelem)
      integer :: d
c
      mm  = n+p+1
      ng  = p+1
      d = 0
      !write(*,*)'Form1DMassMatrix:p,n,mm',p,n,mm
      !write(*,*)'Form1DMassMatrix:U',U
      call BasisData(p,mm,U,d,ng,nelem,O,J,W,X,NN) 

      M = 0.d0
c loop over elements
      do e = 1,nelem
c loop over Gauss points
         do k = 1,ng
c loop over shape functions over elements (p+1 functions)
           do a = 0,p
c loop over shape functions over elements (p+1 functions)
              do b = 0,p
c O(e) + a = first dof of element + 1st local shape function index
c O(e) + b = first dof of element + 2nd local shape function index
c NN(0,a,k,e) = value of shape function a at Gauss point k over element e
c NN(0,b,k,e) = value of shape function b at Gauss point k over element e
c W(k) weight for Gauss point k
c J(e) jacobian ? for element e

c M is a band storage matrix A(i,j)->M(KL+KU+1+i-j,j)
c                 write(*,*)'i,j',O(e)+a,O(e)+b,'->',
c     .             KL+KU+O(e)+a-(O(e)+b),O(e)+b
                 M(KL+KU+O(e)+a-(O(e)+b),O(e)+b) = 
     .             M(KL+KU+O(e)+a-(O(e)+b),O(e)+b) +
     .               NN(0,a,k,e)*NN(0,b,k,e)*J(e)*W(k)
c                 M(O(e)+a,O(e)+b) = 
c     .             M(O(e)+a,O(e)+b) + NN(0,a,k,e)*NN(0,b,k,e)*J(e)*W(k)
              end do
           end do
        end do
      end do

      end subroutine Form1DMassMatrix


c Ux,Uy,Uz - knots vectors
c px,py,py - orders
c nx,ny,nz - number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
c nelemx,nelemy,nelemz - number of elements
c F output rhs (multiple vectors)
      subroutine Form3DRHS
     .   (Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz,F,t,Deltat,Ut,
     .    Kq)
      implicit none
      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real   (kind=8), intent(in) :: Ux(0:nx+px+1)
      real   (kind=8), intent(in) :: Uy(0:ny+py+1)
      real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
      double precision, intent(out) :: F(0:(nx+1)-1,0:(ny+1)*(nz+1)-1)
      double precision :: t,Deltat
      double precision, intent(in) :: Ut(0:(nx+1)-1,0:(ny+1)*(nz+1)-1)
      real   (kind=8) :: Kq(0:(nx+1)-1, 0:(ny+1)*(nz+1)-1)
      integer(kind=4) :: mx,my,mz,ngx,ngy,ngz,ex,ey,ez
      integer(kind=4) :: kx,ky,kz,ax,ay,az,d
      integer(kind=4) :: Ox(nelemx),Oy(nelemy),Oz(nelemz)
      real   (kind=8) :: Jx(nelemx),Jy(nelemy),Jz(nelemz)
      real   (kind=8) :: Wx(px+1),Wy(py+1),Wz(pz+1)
      real   (kind=8) :: Xx(px+1,nelemx)
      real   (kind=8) :: Xy(py+1,nelemy)
      real   (kind=8) :: Xz(pz+1,nelemz)
      real   (kind=8) :: NNx(0:2,0:px,px+1,nelemx),
     .                   NNy(0:2,0:py,py+1,nelemy),
     .                   NNz(0:2,0:pz,pz+1,nelemz)
      real   (kind=8) :: J,W,ftvalue,Utvalue,Laplasjanvalue,kqval,mi,rhs
      real   (kind=8) :: dxu, dyu, dzu, dxkq, dykq, dzkq, change
      integer(kind=4) :: bx,by,bz,cx,cy,cz
      integer :: ind,ind1,ind23,inda,ind1a,ind23a
      integer :: idebug
      idebug=0      
      d=2
      mx  = nx+px+1
      ngx = px+1
      my  = ny+py+1
      ngy = py+1
      mz  = nz+pz+1
      ngz = pz+1

      change = 0.d0

      call BasisData(px,mx,Ux,d,ngx,nelemx,Ox,Jx,Wx,Xx,NNx)
      call BasisData(py,my,Uy,d,ngy,nelemy,Oy,Jy,Wy,Xy,NNy)
      call BasisData(pz,mz,Uz,d,ngz,nelemz,Oz,Jz,Wz,Xz,NNz)

      F = 0.d0
      do ex = 1,nelemx !loops over elements in x,y,z
         do ey = 1,nelemy
            do ez = 1,nelemz
               J = Jx(ex)*Jy(ey)*Jz(ez)
               do kx = 1,ngx !loops over p+1 Gaussian points (= number of B-spline functions over element
                  do ky = 1,ngy
                     do kz = 1,ngz
                        W = Wx(kx)*Wy(ky)*Wz(kz)
                        ftvalue = fvalue(Xx(kx,ex),Xy(ky,ey),Xz(kz,ez))
                        do ax = 0,px ! px+1 functions over element
                           do ay = 0,py
                              do az = 0,pz
c                                 ind = (Oz(ez)+az)*(ny+1)*(nx+1)
c                                 ind = ind + (Oy(ey)+ay)*(nx+1)
c                                 ind = ind + (Ox(ex)+ax)
c                                 F(ind) = F(ind) +
c F is a multiple columns vector
                                 ind1 = (Ox(ex)+ax) !along x
                                 ind23 = (Oy(ey)+ay) !along y
                                 ind23 = ind23 + (Oz(ez)+az)*(ny+1) !along z
c                          write(*,*)'ind',ind,'->(',ind1,',',ind23,')'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c STANDART PROJECTION
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                 F(ind1,ind23) = F(ind1,ind23) +
c     .                               NNx(0,ax,kx,ex)*
c     .                               NNy(0,ay,ky,ey)*
c     .                               NNz(0,az,kz,ez)*J*W*ftvalue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c TIME DEPENDENT HEAT EQUATION
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Ut+1=Ut + Deltat*Laplasjan Ut + Deltat Ut
c  Int_Omega Ut*Vdx + Deltat Int_Omega Laplasjan Ut*Vdx + Deltat Int_Omega ft*V dx
c        (A)                     (B)                          (C)
c
c                   go to 13
c Computations of (A)
                                Utvalue=0.d0
                                kqval = 0.d0
                                dxu=0.d0
                                dyu=0.d0
                                dzu=0.d0
                                dxkq=0.d0
                                dykq=0.d0
                                dzkq=0.d0

                                do bx = 0,px ! px+1 functions over element
                                   do by = 0,py
                                      do bz = 0,pz
                                 ind1a = (Ox(ex)+bx) !along x
                                 ind23a = (Oy(ey)+by) !along y
                                 ind23a = ind23a + (Oz(ez)+bz)*(ny+1) !along z

                                 ! U
                                 Utvalue = Utvalue +
     .                             NNx(0,bx,kx,ex)*
     .                             NNy(0,by,ky,ey)*
     .                             NNz(0,bz,kz,ez)*Ut(ind1a,ind23a)

                                 dxu = dxu +
     .                             NNx(1,bx,kx,ex)*
     .                             NNy(0,by,ky,ey)*
     .                             NNz(0,bz,kz,ez)*Ut(ind1a,ind23a)

                                 dyu = dyu +
     .                             NNx(0,bx,kx,ex)*
     .                             NNy(1,by,ky,ey)*
     .                             NNz(0,bz,kz,ez)*Ut(ind1a,ind23a)

                                 dzu = dzu +
     .                             NNx(0,bx,kx,ex)*
     .                             NNy(0,by,ky,ey)*
     .                             NNz(1,bz,kz,ez)*Ut(ind1a,ind23a)

                                 Laplasjanvalue = Laplasjanvalue +
     .                               (NNx(2,bx,kx,ex)* !Laplasjan of basis function
     .                                NNy(0,by,ky,ey)*
     .                                NNz(0,bz,kz,ez)+
     .                                NNx(0,bx,kx,ex)*
     .                                NNy(2,by,ky,ey)*
     .                                NNz(0,bz,kz,ez)+
     .                                NNx(0,bx,kx,ex)*
     .                                NNy(0,by,ky,ey)*
     .                                NNz(2,bz,kz,ez))*Ut(ind1a,ind23a) !coefficients

                                 ! Kq
                                 kqval = kqval +
     .                             NNx(0,bx,kx,ex)*
     .                             NNy(0,by,ky,ey)*
     .                             NNz(0,bz,kz,ez)*Kq(ind1a,ind23a)

                                 dxkq = dxkq +
     .                             NNx(1,bx,kx,ex)*
     .                             NNy(0,by,ky,ey)*
     .                             NNz(0,bz,kz,ez)*Kq(ind1a,ind23a)

                                 dykq = dykq +
     .                             NNx(0,bx,kx,ex)*
     .                             NNy(1,by,ky,ey)*
     .                             NNz(0,bz,kz,ez)*Kq(ind1a,ind23a)

                                 dzkq = dzkq +
     .                             NNx(0,bx,kx,ex)*
     .                             NNy(0,by,ky,ey)*
     .                             NNz(1,bz,kz,ez)*Kq(ind1a,ind23a)

                                       enddo
                                    enddo
                                 enddo

             rhs = Deltat * exp(mi * Utvalue) * (
     .          mi * kqval  * (dxu * dxu + dyu * dyu + dzu * dzu) +     ! mi kq |\/u|^2
     .          dxkq * dxu + dykq * dyu + dzkq * dzu +                  ! \/kq * \/u
     .          kqval * Laplasjanvalue) +                               ! kq /\u
     .          Deltat * h(Xx(kx,ex),Xy(ky,ey),Xz(kz,ez))

             change = change + F(ind1,ind23) +
     .           NNx(0,ax,kx,ex)*
     .           NNy(0,ay,ky,ey)*
     .           NNz(0,az,kz,ez)*J*W*rhs**2

             rhs = rhs + Utvalue

                                 F(ind1,ind23) = F(ind1,ind23) +
     .                               NNx(0,ax,kx,ex)* !represents Bi, Bj, Bk for v
     .                               NNy(0,ay,ky,ey)*
     .                               NNz(0,az,kz,ez)*J*W*rhs

c Boundary compensation
c       u = u* + ub
c       u(t+1)  - u(t)  = dt ( /\u(t) + f )
c       u*(t+1) - u*(t) = dt ( /\u*(t) + /\ub + f)
c                                       ------
c                                       only change
c
c       /\ub = 0 for constant boundary conditions, so no change is
c       necessary to impose it
c
        if(idebug.eq.1)then
        write(*,*)'ftvalue,Utvalue,Laplasjanvalue',
     .    ftvalue,Utvalue,Laplasjanvalue
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 13     continue
                              end  do
                           end do
                        end do

                     end do
                  end do
               end do         

            end do
         end do
      end do

      write(*,*)'L2 DIFFERENCE = ', change
  
      end subroutine Form3DRHS

      end module projection_engine

