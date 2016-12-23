module reorderRHS

use basis
use parallelism

implicit none

contains

subroutine ReorderRHSForY(ibegx,iendx,ibegy,iendy,ibegz,iendz,F,F2)
integer, intent(in) :: ibegx,ibegy,ibegz
integer, intent(in) :: iendx,iendy,iendz
real   (kind=8), intent(in) :: F(0:(iendx-ibegx+1)-1,    &
  0:(iendy-ibegy+1)*(iendz-ibegz+1)-1)
real   (kind=8), intent(out) :: F2(0:(iendy-ibegy+1)-1,  &
  0:(iendx-ibegx+1)*(iendz-ibegz+1)-1)
integer :: ix, iy, iz
integer :: ind2, ind13, ind1, ind23

  do ix = ibegx,iendx
  do iy = ibegy,iendy
  do iz = ibegz,iendz

    ind2 = iy-ibegy
    ind13 = (ix-ibegx)+(iz-ibegz)*(iendx-ibegx+1)

    ind1 = ix-ibegx
    ind23 = (iy-ibegy)+(iz-ibegz)*(iendy-ibegy+1)

    F2(ind2,ind13) = F(ind1,ind23)
  enddo
  enddo
  enddo

end subroutine


subroutine ReorderRHSForZ(ibegx,iendx,ibegy,iendy,ibegz,iendz,F,F2)
integer, intent(in) :: ibegx,ibegy,ibegz
integer, intent(in) :: iendx,iendy,iendz
real   (kind=8), intent(in) :: F(0:(iendy-ibegy+1)-1,  &
  0:(iendx-ibegx+1)*(iendz-ibegz+1)-1)
real   (kind=8), intent(out) :: F2(0:(iendz-ibegz+1)-1, &
  0:(iendx-ibegx+1)*(iendy-ibegy+1)-1)
integer :: ix, iy, iz
integer :: ind3, ind12, ind2, ind13

  do ix = ibegx,iendx
  do iy = ibegy,iendy
  do iz = ibegz,iendz

    ind3 = iz-ibegz
    ind12 = (ix-ibegx)+(iy-ibegy)*(iendx-ibegx+1)

    ind2 = iy-ibegy
    ind13 = (ix-ibegx)+(iz-ibegz)*(iendx-ibegx+1)

    F2(ind3,ind12) = F(ind2,ind13)

  enddo
  enddo
  enddo

end subroutine


subroutine ReorderRHSForX(ibegx,iendx,ibegy,iendy,ibegz,iendz,F,F2)
integer, intent(in) :: ibegx,ibegy,ibegz
integer, intent(in) :: iendx,iendy,iendz
real   (kind=8), intent(in) :: F(0:(iendz-ibegz+1)-1,    &
  0:(iendx-ibegx+1)*(iendy-ibegy+1)-1)
real   (kind=8), intent(out) :: F2(0:(iendx-ibegx+1)-1,  &
  0:(iendy-ibegy+1)*(iendz-ibegz+1)-1)
integer :: ix, iy, iz
integer :: ind1, ind23, ind3, ind12

  do ix = ibegx,iendx
  do iy = ibegy,iendy
  do iz = ibegz,iendz

    ind1 = ix-ibegx
    ind23 = (iy-ibegy)+(iz-ibegz)*(iendy-ibegy+1)

    ind3 = iz-ibegz
    ind12 = (ix-ibegx)+(iy-ibegy)*(iendx-ibegx+1)

    F2(ind1,ind23) = F(ind3,ind12)

  enddo
  enddo
  enddo

end subroutine

end module
