      module reorderRHS

      use basis
      use parallelism
          
      contains
      
      
c Ux,Uy,Uz - knots vectors
c px,py,py - orders
c nx,ny,nz - number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
c nelemx,nelemy,nelemz - number of elements
c F input rhs (multiple vectors) 
c F2 output rhs
      subroutine ReorderRHSForY
     .   (Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz, 
     .    ibegx,iendx,nrankx,nrpx,
     .    ibegy,iendy,nranky,nrpy,
     .    ibegz,iendz,nrankz,nrpz,F,F2)
      implicit none
      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real   (kind=8), intent(in) :: Ux(0:nx+px+1)
      real   (kind=8), intent(in) :: Uy(0:ny+py+1)
      real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
      double precision, intent(in) :: F(0:(iendx-ibegx+1)-1,
     .   0:(iendy-ibegy+1)*(iendz-ibegz+1)-1)
      double precision, intent(out) :: F2(0:(iendy-ibegy+1)-1,
     .   0:(iendx-ibegx+1)*(iendz-ibegz+1)-1)
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
      integer, intent(in) :: ibegx,ibegy,ibegz
      integer, intent(in) :: iendx,iendy,iendz
      integer, intent(in) :: nrankx,nranky,nrankz
      integer, intent(in) :: nrpx,nrpy,nrpz
      integer :: nreppx,nreppy,nreppz !# elements per proc along x,y,z
      real   (kind=8) :: J,W,value
      integer :: ind,ind1,ind23,ind2,ind13,indx,indy,indz
      integer :: i
      integer :: iprint
c
      iprint=0
c      if(MYRANK.eq.2)iprint=1
c
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'On entry to ReorderRHSForY'
        write(*,*)PRINTRANK,'Ux',Ux
        write(*,*)PRINTRANK,'px',px
        write(*,*)PRINTRANK,'nx',nx
        write(*,*)PRINTRANK,'Uy',Uy
        write(*,*)PRINTRANK,'py',py
        write(*,*)PRINTRANK,'ny',ny
        write(*,*)PRINTRANK,'Uz',Uz
        write(*,*)PRINTRANK,'pz',pz
        write(*,*)PRINTRANK,'nz',nz
        write(*,*)PRINTRANK,'ibegx,iendx,MYRANKX,NRPROCX',
     .    ibegx,iendx,MYRANKX,NRPROCX
        write(*,*)PRINTRANK,'ibegy,iendy,MYRANKY,NRPROCY',
     .    ibegy,iendy,MYRANKY,NRPROCY
        write(*,*)PRINTRANK,'ibegz,iendz,MYRANKZ,NRPROCZ',
     .    ibegz,iendz,MYRANKZ,NRPROCZ
        write(*,*)PRINTRANK,'F_out',F
        write(*,*)PRINTRANK,'F2',F2
c        do i=0,iendx-ibegx
c          write(*,*)i,'row=',F(i,0:(iendy-ibegy+1)*(iendz-ibegz+1)-1)
c        enddo
      endif
c
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
      if(iprint.eq.1)then
         write(*,*)PRINTRANK,'ibegx,iendx,ibegy,iendy,ibegz,iendz',
     .     ibegx,iendx,ibegy,iendy,ibegz,iendz
         write(*,*)PRINTRANK,'nelemx,nelemy,nelemz',nelemx,nelemy,nelemz
         write(*,*)PRINTRANK,'nrpx,nrpy,nrpz',nrpx,nrpy,nrpz
         write(*,*)PRINTRANK,'nreppx,nreppy,nreppz',nreppx,nreppy,nreppz
         write(*,*)PRINTRANK,'ex=',
     .         max(nreppx*nrankx-px+1,1),
     .         min(nelemx,nreppx*(nrankx+1)+px)
         write(*,*)PRINTRANK,'ey=',
     .         max(nreppy*nranky-py+1,1),
     .         min(nelemy,nreppy*(nranky+1)+py)
         write(*,*)PRINTRANK,'ez=',
     .         max(nreppz*nrankz-pz+1,1),
     .         min(nelemz,nreppz*(nrankz+1)+pz)
      endif
      F2 = 0.d0
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
                do ax = 0,px
                   do ay = 0,py
                      do az = 0,pz
c                         ind = (Oy(ey)+ay)+(Ox(ex)+ax)*(ny+1)+
c     .                         (Oz(ez)+az)*(nx+1)*(ny+1)
c                         call global2local(ind,indx,indy,indz,
c     .                        nx,ny,nz)
c                         if(iprint.eq.1)then
c                            write(*,*)PRINTRANK,'NEW ind',ind,'->x,y,z',
c     .                        indx,indy,indz,'->'
c                         endif
                         indx=(Ox(ex)+ax)
                         indy=(Oy(ey)+ay)
                         indz=(Oz(ez)+az)
                         if(indx.lt.ibegx-1.or.
     .                      indx.gt.iendx-1)cycle
                         if(indy.lt.ibegy-1.or.
     .                      indy.gt.iendy-1)cycle
                         if(indz.lt.ibegz-1.or.
     .                      indz.gt.iendz-1)cycle
                         ind2 = indy-ibegy+1
                         ind13 = (indx-ibegx+1)+
     .                           (indz-ibegz+1)*(iendx-ibegx+1)
c                         if(iprint.eq.1)then
c                            write(*,*)PRINTRANK,'NEW ind->',ind2,ind13
c                         endif
c                        old ordering
c                         ind = (Ox(ex)+ax)+(Oy(ey)+ay)*(nx+1)+
c     .                         (Oz(ez)+az)*(nx+1)*(ny+1)
c                         call global2local(ind,indx,indy,indz,
c     .                        nx,ny,nz)
c                         indx=(Ox(ex)+ax))
c                         indy=(Oy(ey)+ay))
c                         indz=(Oz(ez)+az))
c                         if(indx.ne.(Ox(ex)+ax))stop
c                         if(indy.ne.(Oy(ey)+ay))stop
c                         if(indz.ne.(Oz(ez)+az))stop
                         if(iprint.eq.1)then
                            write(*,*)PRINTRANK,'OLD ind->x,y,z',
     .                        indx,indy,indz,'->'
                         endif
c                         if(indx.lt.ibegx-1.or.
c     .                      indx.gt.iendx-1)cycle
c                         if(indy.lt.ibegy-1.or.
c     .                      indy.gt.iendy-1)cycle
c                         if(indz.lt.ibegz-1.or.
c     .                      indz.gt.iendz-1)cycle
                         ind1 = indx-ibegx+1
                         ind23 = (indy-ibegy+1)+
     .                           (indz-ibegz+1)*(iendy-ibegy+1)
                         if(iprint.eq.1)then
                            write(*,*)PRINTRANK,'OLD ->',ind1,ind23
                         endif
                         F2(ind2,ind13) = F(ind1,ind23) 
                      end  do
                   end do
                end do
            end do
         end do
      end do

      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'On exit from ReorderRHSForY'
        do i=0,iendy-ibegy
          write(*,*)i,'row=',F2(i,0:(iendx-ibegx+1)*(iendz-ibegz+1)-1)
        enddo
      endif

      end subroutine ReorderRHSforY


c Ux,Uy,Uz - knots vectors
c px,py,py - orders
c nx,ny,nz - number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
c nelemx,nelemy,nelemz - number of elements
c F input rhs (multiple vectors) 
c F2 output rhs
      subroutine ReorderRHSForZ
     .   (Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz, 
     .    ibegx,iendx,nrankx,nrpx,
     .    ibegy,iendy,nranky,nrpy,
     .    ibegz,iendz,nrankz,nrpz,F2,F3)
      implicit none
      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real   (kind=8), intent(in) :: Ux(0:nx+px+1)
      real   (kind=8), intent(in) :: Uy(0:ny+py+1)
      real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
      double precision, intent(in) :: F2(0:(iendy-ibegy+1)-1,
     .    0:(iendx-ibegx+1)*(iendz-ibegz+1)-1)
      double precision, intent(out) :: F3(0:(iendz-ibegz+1)-1,
     .    0:(iendx-ibegx+1)*(iendy-ibegy+1)-1)
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
      integer, intent(in) :: ibegx,ibegy,ibegz
      integer, intent(in) :: iendx,iendy,iendz
      integer, intent(in) :: nrankx,nranky,nrankz
      integer, intent(in) :: nrpx,nrpy,nrpz
      integer :: nreppx,nreppy,nreppz !# elements per proc along x,y,z
      real   (kind=8) :: J,W,value
      integer :: ind,ind3,ind12,ind2,ind13,indx,indy,indz
      integer :: iprint
      integer :: i
c
      iprint=0
c      if(MYRANK.eq.0)iprint=1
c
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'On entry to ReorderRHSForZ'
        write(*,*)PRINTRANK,'Ux',Ux
        write(*,*)PRINTRANK,'px',px
        write(*,*)PRINTRANK,'nx',nx
        write(*,*)PRINTRANK,'Uy',Uy
        write(*,*)PRINTRANK,'py',py
        write(*,*)PRINTRANK,'ny',ny
        write(*,*)PRINTRANK,'Uz',Uz
        write(*,*)PRINTRANK,'pz',pz
        write(*,*)PRINTRANK,'nz',nz
        write(*,*)PRINTRANK,'ibegx,iendx,MYRANKX,NRPROCX',
     .    ibegx,iendx,MYRANKX,NRPROCX
        write(*,*)PRINTRANK,'ibegy,iendy,MYRANKY,NRPROCY',
     .    ibegy,iendy,MYRANKY,NRPROCY
        write(*,*)PRINTRANK,'ibegz,iendz,MYRANKZ,NRPROCZ',
     .    ibegz,iendz,MYRANKZ,NRPROCZ
        write(*,*)PRINTRANK,'F2',F2
c        write(*,*)PRINTRANK,'F2',F3
c        do i=0,iendx-ibegx
c          write(*,*)i,'row=',F(i,0:(iendy-ibegy+1)*(iendz-ibegz+1)-1)
c        enddo
      endif
c
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
      F3 = 0.d0
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
                do ax = 0,px
                   do ay = 0,py
                      do az = 0,pz
c                        new ordering
c                         ind = (Oz(ez)+az)+(Ox(ex)+ax)*(nz+1)+
c     .                         (Oz(ey)+ay)*(nx+1)*(nz+1)
c                         call global2local(ind,indx,indy,indz,
c     .                        nx,ny,nz)
                         indx=(Ox(ex)+ax)
                         indy=(Oy(ey)+ay)
                         indz=(Oz(ez)+az)
                         if(indx.lt.ibegx-1.or.
     .                      indx.gt.iendx-1)cycle
                         if(indy.lt.ibegy-1.or.
     .                      indy.gt.iendy-1)cycle
                         if(indz.lt.ibegz-1.or.
     .                      indz.gt.iendz-1)cycle
                         ind3 = indz-ibegz+1
                         ind12 = (indx-ibegx+1)+
     .                           (indy-ibegy+1)*(iendx-ibegx+1)
c                         if(iprint.eq.1)then
c                            write(*,*)PRINTRANK,'ind->x,y,z',
c     .                        ind,indx,indy,indz,'->',ind3,ind12
c                         endif
c                        old ordering
c                         ind = (Oy(ey)+ay)+(Ox(ex)+ax)*(ny+1)+
c     .                         (Oz(ez)+az)*(nx+1)*(ny+1)
c                         call global2local(ind,indx,indy,indz,
c     .                        nx,ny,nz)
c                         if(indx.lt.ibegx-1.or.
c     .                      indx.gt.iendx-1)cycle
c                         if(indy.lt.ibegy-1.or.
c     .                      indy.gt.iendy-1)cycle
c                         if(indz.lt.ibegz-1.or.
c     .                      indz.gt.iendz-1)cycle
                         ind2 = indy-ibegy+1
                         ind13 = (indx-ibegx+1)+
     .                           (indz-ibegz+1)*(iendx-ibegx+1)
c                         if(iprint.eq.1)then
c                            write(*,*)PRINTRANK,'ind->x,y,z',
c     .                        ind,indx,indy,indz,'->',ind2,ind13
c                         endif
                         F3(ind3,ind12) = F2(ind2,ind13) 
                      end  do
                   end do
                end do
            end do
         end do
      end do

      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'On exit from ReorderRHSForZ'
       do i=0,iendz-ibegz
          write(*,*)i,'row=',F3(i,0:(iendx-ibegx+1)*(iendy-ibegy+1)-1)
        enddo
      endif
  
      end subroutine ReorderRHSforZ


c Ux,Uy,Uz - knots vectors
c px,py,py - orders
c nx,ny,nz - number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
c nelemx,nelemy,nelemz - number of elements
c F input rhs (multiple vectors) 
c F2 output rhs
      subroutine ReorderRHSForX
     .   (Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz, 
     .    ibegx,iendx,nrankx,nrpx,
     .    ibegy,iendy,nranky,nrpy,
     .    ibegz,iendz,nrankz,nrpz,F,F2)
      implicit none
      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real   (kind=8), intent(in) :: Ux(0:nx+px+1)
      real   (kind=8), intent(in) :: Uy(0:ny+py+1)
      real   (kind=8), intent(in) :: Uz(0:nz+pz+1)
      double precision, intent(in) :: F(0:(iendz-ibegz+1)-1,
     .   0:(iendx-ibegx+1)*(iendy-ibegy+1)-1)
      double precision, intent(out) :: F2(0:(iendx-ibegx+1)-1,
     .   0:(iendy-ibegy+1)*(iendz-ibegz+1)-1)
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
      integer, intent(in) :: ibegx,ibegy,ibegz
      integer, intent(in) :: iendx,iendy,iendz
      integer, intent(in) :: nrankx,nranky,nrankz
      integer, intent(in) :: nrpx,nrpy,nrpz
      integer :: nreppx,nreppy,nreppz !# elements per proc along x,y,z
      real   (kind=8) :: J,W,value
      integer :: ind,ind3,ind12,ind1,ind23,indx,indy,indz
      integer :: i
      integer :: iprint
c
      iprint=0
c
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
       do ex = max(nreppx*nrankx-px+1,1),
     .         min(nelemx,nreppx*(nrankx+1)+px)
         do ey = max(nreppy*nranky-py+1,1),
     .           min(nelemy,nreppy*(nranky+1)+py)
            do ez = max(nreppz*nrankz-pz+1,1),
     .              min(nelemz,nreppz*(nrankz+1)+pz)
                do ax = 0,px
                   do ay = 0,py
                      do az = 0,pz
                         indx=(Ox(ex)+ax)
                         indy=(Oy(ey)+ay)
                         indz=(Oz(ez)+az)
                         if(indx.lt.ibegx-1.or.
     .                      indx.gt.iendx-1)cycle
                         if(indy.lt.ibegy-1.or.
     .                      indy.gt.iendy-1)cycle
                         if(indz.lt.ibegz-1.or.
     .                      indz.gt.iendz-1)cycle
                         ind1 = indx-ibegx+1
                         ind23 = (indy-ibegy+1)+
     .                           (indz-ibegz+1)*(iendy-ibegy+1)
                         ind3 = indz-ibegz+1
                         ind12 = (indx-ibegx+1)+
     .                           (indy-ibegy+1)*(iendx-ibegx+1)
                         if(iprint.eq.1)then
                            write(*,*)PRINTRANK,'OLD ->',ind1,ind23
                         endif
                         F2(ind1,ind23) = F(ind3,ind12) 
                      end  do
                   end do
                end do
            end do
         end do
      end do

      end subroutine ReorderRHSforX

      end module reorderRHS
