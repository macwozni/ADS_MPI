module RHS_eq

implicit none
real (kind=8) :: drained = 0
real (kind=8) :: l2norm

contains


subroutine ComputePointForRHS(Xx,Xy,Xz,px,py,pz,kx,ky,kz, &
   ex,ey,ez,nelemx,nelemy,nelemz,Uval,v,ax,ay,az,NNx,NNy,NNz, &
   minex,miney,minez,Kq,maxex,maxey,maxez,t,Dt,mi,dux,duy,duz, &
   ibegx,ibegy,ibegz,iendx,iendy,iendz,F,J,W)
use input_data, ONLY : pumping, draining, initial_state
implicit none
integer(kind=4), intent(in)  :: px,py,pz
real   (kind=8), intent(in)  :: Xx(px+1,nelemx)
real   (kind=8), intent(in)  :: Xy(py+1,nelemy)
real   (kind=8), intent(in)  :: Xz(pz+1,nelemz)
integer(kind=4), intent(in)  :: kx,ky,kz
integer(kind=4), intent(in)  :: ex,ey,ez
integer(kind=4), intent(in)  :: nelemx,nelemy,nelemz
real   (kind=8), intent(in)  :: Uval
real   (kind=8), intent(inout)  :: v
integer(kind=4), intent(in)  :: ibegx,ibegy,ibegz
integer(kind=4), intent(in)  :: iendx,iendy,iendz
real   (kind=8), intent(out) :: F
integer(kind=4), intent(in)  :: maxex,maxey,maxez
real   (kind=8), intent(in)  :: Kq(px+1,py+1,pz+1,maxex-minex+1,maxey-miney+1,maxez-minez+1)
integer(kind=4), intent(in)  :: minex,miney,minez
integer(kind=4), intent(in)  :: ax,ay,az
real   (kind=8), intent(in)  :: dux,duy,duz
real   (kind=8), intent(in)  :: t,Dt,mi,J,W
real   (kind=8), intent(in)  :: NNx(0:1,0:px,px+1,nelemx), &
                   NNy(0:1,0:py,py+1,nelemy), &
                   NNz(0:1,0:pz,pz+1,nelemz)
real   (kind=8) :: fval,vpump,vdrain,kqval
real   (kind=8) :: Umax = -1d10, Umin = 1d10
real   (kind=8) :: dvx,dvy,dvz,rhs

        vpump = pumping(Xx(kx,ex),Xy(ky,ey),Xz(kz,ez))    
        Umax = max(Umax, Uval)
        Umin = min(Umin, Uval)

        v   = NNx(0,ax,kx,ex) * NNy(0,ay,ky,ey) * NNz(0,az,kz,ez)
        dvx = NNx(1,ax,kx,ex) * NNy(0,ay,ky,ey) * NNz(0,az,kz,ez) 
        dvy = NNx(0,ax,kx,ex) * NNy(1,ay,ky,ey) * NNz(0,az,kz,ez) 
        dvz = NNx(0,ax,kx,ex) * NNy(0,ay,ky,ey) * NNz(1,az,kz,ez) 

        kqval = Kq(kx,ky,kz,ex-minex+1,ey-miney+1,ez-minez+1)
        vdrain = max(0.d0, draining(Uval, Xx(kx,ex),Xy(ky,ey),Xz(kz,ez)))
        fval = vpump - vdrain
        !--- Real
        if (t > 0.0) then
          rhs = Dt * ( - kqval * exp(mi * Uval) * (dux*dvx + duy*dvy + duz*dvz) + v * fval)
          F = J*W*(v * Uval + rhs)

          drained = drained + J*W*v*Dt*vdrain
          l2norm = l2norm + J*W*v*Uval*Uval
        else
          fval = initial_state(Xx(kx,ex),Xy(ky,ey),Xz(kz,ez))
          F= J*W*v*fval
          l2norm = l2norm + J*W*v*fval*fval
        endif

end subroutine


end module