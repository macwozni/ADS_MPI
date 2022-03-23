module reorderRHS

   implicit none

contains

   subroutine ReorderRHSForY(ibeg, iend, F, F2)
      implicit none
      integer(kind = 4), intent(in), dimension(3) :: ibeg
      integer(kind = 4), intent(in), dimension(3) :: iend
      real (kind = 8), intent(in) :: F(0:(iend(1) - ibeg(1) + 1) - 1, &
      0:(iend(2) - ibeg(2) + 1)*(iend(3) - ibeg(3) + 1) - 1)
      real (kind = 8), intent(out) :: F2(0:(iend(2) - ibeg(2) + 1) - 1, &
      0:(iend(1) - ibeg(1) + 1)*(iend(3) - ibeg(3) + 1) - 1)
      integer(kind = 4) :: ix, iy, iz
      integer(kind = 4) :: ind2, ind13, ind1, ind23

      do ix = ibeg(1), iend(1)
         do iy = ibeg(2), iend(2)
            do iz = ibeg(3), iend(3)

               ind2 = iy - ibeg(2)
               ind13 = (ix - ibeg(1))+(iz - ibeg(3))*(iend(1) - ibeg(1) + 1)

               ind1 = ix - ibeg(1)
               ind23 = (iy - ibeg(2))+(iz - ibeg(3))*(iend(2) - ibeg(2) + 1)

               F2(ind2, ind13) = F(ind1, ind23)
            enddo
         enddo
      enddo

   end subroutine


   subroutine ReorderRHSForZ(ibeg, iend, F, F2)
      implicit none
      integer(kind = 4), intent(in), dimension(3) :: ibeg
      integer(kind = 4), intent(in), dimension(3) :: iend
      real (kind = 8), intent(in) :: F(0:(iend(2) - ibeg(2) + 1) - 1, &
      0:(iend(1) - ibeg(1) + 1)*(iend(3) - ibeg(3) + 1) - 1)
      real (kind = 8), intent(out) :: F2(0:(iend(3) - ibeg(3) + 1) - 1, &
      0:(iend(1) - ibeg(1) + 1)*(iend(2) - ibeg(2) + 1) - 1)
      integer(kind = 4) :: ix, iy, iz
      integer(kind = 4) :: ind3, ind12, ind2, ind13

      do ix = ibeg(1), iend(1)
         do iy = ibeg(2), iend(2)
            do iz = ibeg(3), iend(3)

               ind3 = iz - ibeg(3)
               ind12 = (ix - ibeg(1))+(iy - ibeg(2))*(iend(1) - ibeg(1) + 1)

               ind2 = iy - ibeg(2)
               ind13 = (ix - ibeg(1))+(iz - ibeg(3))*(iend(1) - ibeg(1) + 1)

               F2(ind3, ind12) = F(ind2, ind13)

            enddo
         enddo
      enddo

   end subroutine


   subroutine ReorderRHSForX(ibeg, iend, F, F2)
      implicit none
      integer(kind = 4), intent(in), dimension(3) :: ibeg
      integer(kind = 4), intent(in), dimension(3) :: iend
      real (kind = 8), intent(in) :: F(0:(iend(3) - ibeg(3) + 1) - 1, &
      0:(iend(1) - ibeg(1) + 1)*(iend(2) - ibeg(2) + 1) - 1)
      real (kind = 8), intent(out) :: F2(0:(iend(1) - ibeg(1) + 1) - 1, &
      0:(iend(2) - ibeg(2) + 1)*(iend(3) - ibeg(3) + 1) - 1)
      integer(kind = 4) :: ix, iy, iz
      integer(kind = 4) :: ind1, ind23, ind3, ind12

      do ix = ibeg(1), iend(1)
         do iy = ibeg(2), iend(2)
            do iz = ibeg(3), iend(3)

               ind1 = ix - ibeg(1)
               ind23 = (iy - ibeg(2))+(iz - ibeg(3))*(iend(2) - ibeg(2) + 1)

               ind3 = iz - ibeg(3)
               ind12 = (ix - ibeg(1))+(iy - ibeg(2))*(iend(1) - ibeg(1) + 1)

               F2(ind1, ind23) = F(ind3, ind12)

            enddo
         enddo
      enddo

   end subroutine

end module