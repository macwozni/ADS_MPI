      module projection_engine
c
      use gauss
      use basis
      use parallelism
c

c  ...order of approximations
      integer :: ORDER
c  ...number of elements in one dimension
      integer :: SIZE

      contains

      subroutine initialize_parameters
      ORDER=2
      SIZE=2
      end subroutine initialize_parameters

 
      subroutine Form1DMassMatrix(KL,KU,U,p,n,nelem,M)
      !!$ KL = number of subdiagonals, KU = number of superdiagonals
      !!$ d is dimension
      !!$ U is the knot vector
      !!$ p is the polynomial order
      !!$ n is the index of the last control point
      !!$ nelem you get from running CountSpans
      !!$ M is the dense matrix
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
      integer :: nrepp !# elements per processor
      integer :: iprint
c
      iprint=0
c      if(MYRANK.eq.0)iprint=1
c
      mm  = n+p+1
      ng  = p+1
      d = 0
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
c     .               NN(0,a,k,e)*NN(0,b,k,e)*J(e)*W(k)
              end do
           end do
        end do
      end do

      end subroutine Form1DMassMatrix

c-> DEBUG
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
c      real   (kind=8) :: NN(0:d,0:p,p+1,nelem)
      real   (kind=8) :: NN(0:0,0:p,p+1,nelem)
      integer :: d
c
      mm  = n+p+1
      ng  = p+1
      d = 0
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
                 M(O(e)+a,O(e)+b) = 
     .             M(O(e)+a,O(e)+b) + NN(0,a,k,e)*NN(0,b,k,e)*J(e)*W(k)
              end do
           end do
        end do
      end do

      end subroutine Form1DMassMatrixFULL
c<- DEBUG


c Ux,Uy,Uz - knots vectors
c px,py,py - orders
c nx,ny,nz - number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
c nelemx,nelemy,nelemz - number of elements
c F output rhs (multiple vectors)
      subroutine Form3DRHS
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
c      real   (kind=8) :: NNx(0:d,0:px,px+1,nelemx),
c     .                   NNy(0:d,0:py,py+1,nelemy),
c     .                   NNz(0:d,0:pz,pz+1,nelemz)
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
      if(iprint.eq.1)then
        write(*,*)PRINTRANK,'ex:',max(nreppx*nrankx-px+1,1),
     .         min(nelemx,nreppx*(nrankx+1)+px)
        write(*,*)PRINTRANK,'ey:',max(nreppy*nranky-py+1,1),
     .         min(nelemy,nreppy*(nranky+1)+py)
        write(*,*)PRINTRANK,'ez:',max(nreppz*nrankz-pz+1,1),
     .         min(nelemz,nreppz*(nrankz+1)+pz)
        write(*,*)PRINTRANK,'ibegx,iendx',ibegx,iendx
        write(*,*)PRINTRANK,'ibegy,iendy',ibegy,iendy
        write(*,*)PRINTRANK,'ibegz,iendz',ibegz,iendz
      endif
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
                        value = fvalue(Xx(kx,ex),Xy(ky,ey),Xz(kz,ez))
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
                                 if(iprint.eq.1)then
                                   write(*,*)PRINTRANK,'ind->x,y,z',
     .                                ind,indx,indy,indz,'->',ind1,ind23
                                 endif
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
                        if(iprint.eq.1)then
                           write(*,*)PRINTRANK,
     .                      'ind',ind,'->',ind1,ind23,
     .                      ex,ey,ez
                        endif
                              end  do
                           end do
                        end do

                     end do
                  end do
               end do         

            end do
         end do
      end do
  
      end subroutine Form3DRHS

c-> DEBUG
c Ux,Uy,Uz - knots vectors
c px,py,py - orders
c nx,ny,nz - number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
c nelemx,nelemy,nelemz - number of elements
c F output rhs (multiple vectors)
      subroutine Form3DRHSFULL
     .   (Ux,px,nx,nelemx,Uy,py,ny,nelemy,Uz,pz,nz,nelemz,F) 
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
      real   (kind=8) :: NNx(0:0,0:px,px+1,nelemx),
     .                   NNy(0:0,0:py,py+1,nelemy),
     .                   NNz(0:0,0:pz,pz+1,nelemz)
c      real   (kind=8) :: NNx(0:d,0:px,px+1,nelemx),
c     .                   NNy(0:d,0:py,py+1,nelemy),
c     .                   NNz(0:d,0:pz,pz+1,nelemz)
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
                                 F(ind) = F(ind) + 
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
  
      end subroutine Form3DRHSFULL
c<- DEBUG

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
c
      indz = ind / ((nx+1)*(ny+1))
      ind_temp = ind - indz*(nx+1)*(ny+1)
      indy = ind_temp / (nx+1)
      ind_temp = ind_temp - indy*(nx+1)
      indx = ind_temp 
c
      end subroutine global2local

      end module projection_engine

