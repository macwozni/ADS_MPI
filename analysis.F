      module analysis


      use gauss
      use basis
      use projection_engine

      contains


c Ux,Uy,Uz  knots vectors
c px,py,py  orders
c nx,ny,nz  number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
c nelemx,nelemy,nelemz - number of elements
c F coefficients of the function
      subroutine NormL2
     .   (Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz, 
     .    ibegx,iendx,nrankx,nrpx,
     .    ibegy,iendy,nranky,nrpy,
     .    ibegz,iendz,nrankz,nrpz,F)
      implicit none
      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real   (kind=8), intent(in) :: Ux(0:nx+px+1)
      real   (kind=8), intent(in) :: Uy(0:ny+py+1)
      real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
      double precision, intent(out) :: F(0:(iendx-ibegx+1)-1,
     .     0:(iendy-ibegy+1)*(iendz-ibegz+1)-1)
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
      real   (kind=8) :: J,W,value
      integer, intent(in) :: ibegx,ibegy,ibegz
      integer, intent(in) :: iendx,iendy,iendz
      integer, intent(in) :: nrankx,nranky,nrankz
      integer, intent(in) :: nrpx,nrpy,nrpz
      integer :: nreppx,nreppy,nreppz !# elements per proc along x,y,z
      integer :: ind,ind1,ind23,ind23a,indx,indy,indz
      integer :: iprint

      iprint=0
c      if(MYRANK.eq.2)iprint=1

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

c-> parallel number of elements per processors
      nreppx = nelemx/nrpx
      nreppy = nelemy/nrpy
      nreppz = nelemz/nrpz
c<- parallel
      F = 0.d0
c-> parallel
c      do ex = 1,nelemx
       do ex = max(nreppx*nrankx-px+1,1),
     .         min(nelemx,nreppx*(nrankx+1)+px)
c         do ey = 1,nelemy
         do ey = max(nreppy*nranky-py+1,1),
     .           min(nelemy,nreppy*(nranky+1)+py)
c            do ez = 1,nelemz
            do ez = max(nreppz*nrankz-pz+1,1),
     .              min(nelemz,nreppz*(nrankz+1)+pz)
c<- parallel
               J = Jx(ex)*Jy(ey)*Jz(ez)
               do kx = 1,ngx
                  do ky = 1,ngy
                     do kz = 1,ngz
                        W = Wx(kx)*Wy(ky)*Wz(kz)
                        !value = fvalue(Xx(kx,ex),Xy(ky,ey),Xz(kz,ez))
                        do ax = 0,px
                           do ay = 0,py
                              do az = 0,pz
                                 ind = (Ox(ex)+ax)+(Oy(ey)+ay)*(nx+1)+
     .                               (Oz(ez)+az)*(ny+1)*(nx+1)
                                 call global2local(ind,indx,indy,indz,
     .                              nx,ny,nz)
c                                 if(indx.ne.(Ox(ex)+ax))stop
c                                 if(indy.ne.(Oy(ey)+ay))stop
c                                 if(indz.ne.(Oz(ez)+az))stop
                                 if(indx.lt.ibegx-1.or.
     .                              indx.gt.iendx-1)cycle
                                 if(indy.lt.ibegy-1.or.
     .                              indy.gt.iendy-1)cycle
                                 if(indz.lt.ibegz-1.or.
     .                              indz.gt.iendz-1)cycle
                                 ind1 = indx-ibegx+1
                                 ind23 = (indy-ibegy+1)+
     .                             (indz-ibegz+1)*(iendy-ibegy+1)
c-> OLD
cc F is a multiple columns vector
cc-> parallel now we have distributed rhs
cc                                 ind1 = (Ox(ex)+ax) !along x
cc                                 ind23 = (Oy(ey)+ay) !along y
cc                                 ind23 = ind23 + (Oz(ez)+az)*(ny+1) !along z
c                                 ind1 = (Ox(ex)+ax-ibegx+1) !along x
c                                 ind23 = (Oy(ey)+ay-ibegy+1) !along y
c                                 ind23a = (Oz(ez)+az-ibegz+1)
c                                 ind23 =ind23+ind23a*(iendy-ibegy+1) !along z
c<- OLD
                        if(ind1.lt.0.or.ind1.gt.(iendx-ibegx))cycle
                        if(ind23.lt.0.or.ind23.gt.
     .                     (iendy-ibegy+1)*(iendz-ibegz+1)-1)cycle
c<- parallel
                                 F(ind1,ind23) = F(ind1,ind23) + 
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
  
      end subroutine NormL2


      end module analysis

