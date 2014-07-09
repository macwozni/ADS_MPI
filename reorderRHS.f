      module reorderRHS

      use basis
          
      contains
      
c Ux,Uy,Uz - knots vectors
c px,py,py - orders
c nx,ny,nz - number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
c nelemx,nelemy,nelemz - number of elements
c F input rhs (multiple vectors) 
c F2 output rhs
      subroutine ReorderRHSForY
     .   (Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz,F,F2) 
      implicit none
      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real   (kind=8), intent(in) :: Ux(0:nx+px+1)
      real   (kind=8), intent(in) :: Uy(0:ny+py+1)
      real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
      double precision, intent(in) :: F(0:(nx+1)-1,0:(ny+1)*(nz+1)-1)
      double precision, intent(out) :: F2(0:(ny+1)-1,0:(nx+1)*(nz+1)-1)
      integer(kind=4) :: mx,my,mz,ngx,ngy,ngz,ex,ey,ez
      integer(kind=4) :: kx,ky,kz,ax,ay,az,d
      integer(kind=4) :: Ox(nelemx),Oy(nelemy),Oz(nelemz)
      real   (kind=8) :: Jx(nelemx),Jy(nelemy),Jz(nelemz)
      real   (kind=8) :: Wx(px+1),Wy(py+1),Wz(pz+1)
      real   (kind=8) :: Xx(px+1,nelemx)
      real   (kind=8) :: Xy(py+1,nelemy)
      real   (kind=8) :: Xz(pz+1,nelemz)
      real   (kind=8) :: NNx(0:0,0:px,px+1,nelemx),
     .                   NNy(0:0,0:py,py+1,nelemy),
     .                   NNz(0:0,0:pz,pz+1,nelemz)
c      real   (kind=8) :: NNx(0:d,0:px,px+1,nelemx),
c     .                   NNy(0:d,0:py,py+1,nelemy),
c     .                   NNz(0:d,0:pz,pz+1,nelemz)
      real   (kind=8) :: J,W,value
      integer :: ind,ind1,ind23,ind2,ind13
      d=0
      mx  = nx+px+1
      ngx = px+1
      my  = ny+py+1
      ngy = py+1
      mz  = nz+pz+1
      ngz = pz+1

      call BasisData(px,mx,Ux,d,ngx,nelemx,Ox,Jx,Wx,Xx,NNx)
      call BasisData(py,my,Uy,d,ngy,nelemy,Oy,Jy,Wy,Xy,NNy)
      call BasisData(pz,mz,Uz,d,ngz,nelemz,Oz,Jz,Wz,Xz,NNz)

      F2 = 0.d0
      do ex = 1,nelemx
         do ey = 1,nelemy
            do ez = 1,nelemz
                do ax = 0,px
                   do ay = 0,py
                      do az = 0,pz
c                        new ordering
                         ind2 = (Oy(ey)+ay) !along y
                         ind13 = (Ox(ex)+ax) !along x
                         ind13 = ind13 + (Oz(ez)+az)*(nx+1) !along z
c                        old ordering
                         ind1 = (Ox(ex)+ax) !along x
                         ind23 = (Oy(ey)+ay) !along y
                         ind23 = ind23 + (Oz(ez)+az)*(ny+1) !along z
                         F2(ind2,ind13) = F(ind1,ind23) 
                      end  do
                   end do
                end do
            end do
         end do
      end do
  
      end subroutine ReorderRHSforY

c Ux,Uy,Uz - knots vectors
c px,py,py - orders
c nx,ny,nz - number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
c nelemx,nelemy,nelemz - number of elements
c F input rhs (multiple vectors) 
c F2 output rhs
      subroutine ReorderRHSForZ
     .   (Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz,F2,F3) 
      implicit none
      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real   (kind=8), intent(in) :: Ux(0:nx+px+1)
      real   (kind=8), intent(in) :: Uy(0:ny+py+1)
      real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
      double precision, intent(in) :: F2(0:(ny+1)-1,0:(nx+1)*(nz+1)-1)
      double precision, intent(out) :: F3(0:(nz+1)-1,0:(nx+1)*(ny+1)-1)
      integer(kind=4) :: mx,my,mz,ngx,ngy,ngz,ex,ey,ez
      integer(kind=4) :: kx,ky,kz,ax,ay,az,d
      integer(kind=4) :: Ox(nelemx),Oy(nelemy),Oz(nelemz)
      real   (kind=8) :: Jx(nelemx),Jy(nelemy),Jz(nelemz)
      real   (kind=8) :: Wx(px+1),Wy(py+1),Wz(pz+1)
      real   (kind=8) :: Xx(px+1,nelemx)
      real   (kind=8) :: Xy(py+1,nelemy)
      real   (kind=8) :: Xz(pz+1,nelemz)
      real   (kind=8) :: NNx(0:0,0:px,px+1,nelemx),
     .                   NNy(0:0,0:py,py+1,nelemy),
     .                   NNz(0:0,0:pz,pz+1,nelemz)
c      real   (kind=8) :: NNx(0:d,0:px,px+1,nelemx),
c     .                   NNy(0:d,0:py,py+1,nelemy),
c     .                   NNz(0:d,0:pz,pz+1,nelemz)
      real   (kind=8) :: J,W,value
      integer :: ind3,ind12,ind2,ind13
      d=0
      mx  = nx+px+1
      ngx = px+1
      my  = ny+py+1
      ngy = py+1
      mz  = nz+pz+1
      ngz = pz+1

      call BasisData(px,mx,Ux,d,ngx,nelemx,Ox,Jx,Wx,Xx,NNx)
      call BasisData(py,my,Uy,d,ngy,nelemy,Oy,Jy,Wy,Xy,NNy)
      call BasisData(pz,mz,Uz,d,ngz,nelemz,Oz,Jz,Wz,Xz,NNz)

      F3 = 0.d0
      do ex = 1,nelemx
         do ey = 1,nelemy
            do ez = 1,nelemz
                do ax = 0,px
                   do ay = 0,py
                      do az = 0,pz
c                        new ordering
                         ind3 = (Oz(ez)+az) !along z
                         ind12 = (Ox(ex)+ax) !along x
                         ind12 = ind12 + (Oy(ey)+ay)*(nx+1) !along y
c                        old ordering
                         ind2 = (Oy(ey)+ay) !along y
                         ind13 = (Ox(ex)+ax) !along x
                         ind13 = ind13 + (Oz(ez)+az)*(nx+1) !along z
                         F3(ind3,ind12) = F2(ind2,ind13) 
                      end  do
                   end do
                end do
            end do
         end do
      end do
  
      end subroutine ReorderRHSforZ
      
      end module reorderRHS
