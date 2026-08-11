!------------------------------------------------------------------------------
!
! MODULE: reorderRHS
!
! DESCRIPTION:
!> @file reorderRHS.F90
!> @brief Module providing directional reordering routines for
!> right-hand-side arrays in tensor-product storage.
!>
!> @details
!> This module groups helper procedures used to permute local
!> two-dimensional right-hand-side arrays between the storage layouts
!> associated with different leading directions of the tensor-product
!> decomposition.
!>
!> The provided functionality includes:
!> - reordering from an \f$(x, yz)\f$ layout to a \f$(y, xz)\f$ layout
!>   through \ref ReorderRHSForY,
!> - reordering from a \f$(y, xz)\f$ layout to a \f$(z, xy)\f$ layout
!>   through \ref ReorderRHSForZ,
!> - reordering from a \f$(z, xy)\f$ layout back to an \f$(x, yz)\f$
!>   layout through \ref ReorderRHSForX,
!> - normalizing any axis-order permutation to canonical \f$(x, yz)\f$
!>   storage through \ref NormalizeTrialBufferToXYZ.
!>
!> These permutations are used by the ADS workflow to rotate directional
!> right-hand-side data between successive one-dimensional solve stages.
!
!------------------------------------------------------------------------------
module reorderRHS

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reorders a local right-hand-side array from \f$(x, yz)\f$
!> storage to \f$(y, xz)\f$ storage.
!>
!> @details
!> The input array \p F is interpreted as a two-dimensional view of a
!> three-dimensional tensor-product block where:
!> - the first array index corresponds to the first physical direction,
!> - the second array index is a flattened combination of the second and
!>   third physical directions.
!>
!> This routine traverses all locally owned tensor-product indices and
!> rewrites the data into \p F2 so that:
!> - the first array index corresponds to the second physical direction,
!> - the second array index is a flattened combination of the first and
!>   third physical directions.
!>
!> In symbolic form, the mapping is:
!>
!> \f[ F(i_x,\, i_y + n_y i_z) \longrightarrow
!> F2(i_y,\, i_x + n_x i_z), \f]
!>
!> where the local extents \f$n_x\f$ and \f$n_y\f$ are implied by the
!> ownership bounds \p ibeg and \p iend.
!
! Input:
! ------
!> @param[in] ibeg
!> Lower bounds of the locally owned index range in the three
!> directions.
!>
!> @param[in] iend
!> Upper bounds of the locally owned index range in the three
!> directions.
!>
!> @param[in] F
!> Input right-hand-side array stored in \f$(x, yz)\f$ layout.
!
! Output:
! -------
!> @param[out] F2
!> Reordered right-hand-side array stored in \f$(y, xz)\f$ layout.
!
! Notes:
! ------
!> @note
!> The routine performs a pure permutation of local data and does not
!> modify numerical values.
!
!---------------------------------------------------------------------------
subroutine ReorderRHSForY(ibeg, iend, F, F2)
   implicit none
!> @brief Lower bounds of the local ownership range.
   integer(kind = 4), intent(in), dimension(3) :: ibeg
!> @brief Upper bounds of the local ownership range.
   integer(kind = 4), intent(in), dimension(3) :: iend
!> @brief Input array in \f$(x, yz)\f$ storage.
   real (kind = 8), intent(in) :: F(0:(iend(1) - ibeg(1) + 1) - 1, &
         0:(iend(2) - ibeg(2) + 1)*(iend(3) - ibeg(3) + 1) - 1)
!> @brief Output array in \f$(y, xz)\f$ storage.
   real (kind = 8), intent(out) :: F2(0:(iend(2) - ibeg(2) + 1) - 1, &
         0:(iend(1) - ibeg(1) + 1)*(iend(3) - ibeg(3) + 1) - 1)
!> @brief Loop counters over local tensor-product indices.
   integer(kind = 4) :: ix, iy, iz
!> @brief Flattened local indices in the input and output layouts.
   integer(kind = 4) :: ind2, ind13, ind1, ind23
!> @brief Local extents in the first and second physical directions.
   integer(kind = 4) :: n1, n2

   n1 = iend(1) - ibeg(1) + 1
   n2 = iend(2) - ibeg(2) + 1

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
!$OMP PRIVATE(ix,iy,iz,ind2,ind13,ind1,ind23) SCHEDULE(STATIC)
   do iz = ibeg(3), iend(3)
      do ix = ibeg(1), iend(1)
         do iy = ibeg(2), iend(2)

            ind2 = iy - ibeg(2)
            ind13 = (ix - ibeg(1))+(iz - ibeg(3))*n1

            ind1 = ix - ibeg(1)
            ind23 = (iy - ibeg(2))+(iz - ibeg(3))*n2

            F2(ind2, ind13) = F(ind1, ind23)
         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

end subroutine ReorderRHSForY

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reorders a local right-hand-side array from \f$(y, xz)\f$
!> storage to \f$(z, xy)\f$ storage.
!>
!> @details
!> The input array \p F is interpreted as a two-dimensional view of a
!> three-dimensional tensor-product block where:
!> - the first array index corresponds to the second physical direction,
!> - the second array index is a flattened combination of the first and
!>   third physical directions.
!>
!> This routine rewrites the data into \p F2 so that:
!> - the first array index corresponds to the third physical direction,
!> - the second array index is a flattened combination of the first and
!>   second physical directions.
!>
!> In symbolic form, the mapping is:
!>
!> \f[ F(i_y,\, i_x + n_x i_z) \longrightarrow
!> F2(i_z,\, i_x + n_x i_y), \f]
!>
!> with local extents derived from \p ibeg and \p iend.
!
! Input:
! ------
!> @param[in] ibeg
!> Lower bounds of the locally owned index range in the three
!> directions.
!>
!> @param[in] iend
!> Upper bounds of the locally owned index range in the three
!> directions.
!>
!> @param[in] F
!> Input right-hand-side array stored in \f$(y, xz)\f$ layout.
!
! Output:
! -------
!> @param[out] F2
!> Reordered right-hand-side array stored in \f$(z, xy)\f$ layout.
!
!---------------------------------------------------------------------------
subroutine ReorderRHSForZ(ibeg, iend, F, F2)
   implicit none
!> @brief Lower bounds of the local ownership range.
   integer(kind = 4), intent(in), dimension(3) :: ibeg
!> @brief Upper bounds of the local ownership range.
   integer(kind = 4), intent(in), dimension(3) :: iend
!> @brief Input array in \f$(y, xz)\f$ storage.
   real (kind = 8), intent(in) :: F(0:(iend(2) - ibeg(2) + 1) - 1, &
         0:(iend(1) - ibeg(1) + 1)*(iend(3) - ibeg(3) + 1) - 1)
!> @brief Output array in \f$(z, xy)\f$ storage.
   real (kind = 8), intent(out) :: F2(0:(iend(3) - ibeg(3) + 1) - 1, &
         0:(iend(1) - ibeg(1) + 1)*(iend(2) - ibeg(2) + 1) - 1)
!> @brief Loop counters over local tensor-product indices.
   integer(kind = 4) :: ix, iy, iz
!> @brief Flattened local indices in the input and output layouts.
   integer(kind = 4) :: ind3, ind12, ind2, ind13
!> @brief Local extent in the first physical direction.
   integer(kind = 4) :: n1

   n1 = iend(1) - ibeg(1) + 1

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
!$OMP PRIVATE(ix,iy,iz,ind3,ind12,ind2,ind13) SCHEDULE(STATIC)
   do iz = ibeg(3), iend(3)
      do ix = ibeg(1), iend(1)
         do iy = ibeg(2), iend(2)

            ind3 = iz - ibeg(3)
            ind12 = (ix - ibeg(1))+(iy - ibeg(2))*n1

            ind2 = iy - ibeg(2)
            ind13 = (ix - ibeg(1))+(iz - ibeg(3))*n1

            F2(ind3, ind12) = F(ind2, ind13)

         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

end subroutine ReorderRHSForZ

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reorders a local right-hand-side array from \f$(z, xy)\f$
!> storage to \f$(x, yz)\f$ storage.
!>
!> @details
!> The input array \p F is interpreted as a two-dimensional view of a
!> three-dimensional tensor-product block where:
!> - the first array index corresponds to the third physical direction,
!> - the second array index is a flattened combination of the first and
!>   second physical directions.
!>
!> This routine rewrites the data into \p F2 so that:
!> - the first array index corresponds to the first physical direction,
!> - the second array index is a flattened combination of the second and
!>   third physical directions.
!>
!> In symbolic form, the mapping is:
!>
!> \f[ F(i_z,\, i_x + n_x i_y) \longrightarrow
!> F2(i_x,\, i_y + n_y i_z). \f]
!>
!> This permutation closes the cyclic sequence of storage reorders used
!> between successive directional solves.
!
! Input:
! ------
!> @param[in] ibeg
!> Lower bounds of the locally owned index range in the three
!> directions.
!>
!> @param[in] iend
!> Upper bounds of the locally owned index range in the three
!> directions.
!>
!> @param[in] F
!> Input right-hand-side array stored in \f$(z, xy)\f$ layout.
!
! Output:
! -------
!> @param[out] F2
!> Reordered right-hand-side array stored in \f$(x, yz)\f$ layout.
!
!---------------------------------------------------------------------------
subroutine ReorderRHSForX(ibeg, iend, F, F2)
   implicit none
!> @brief Lower bounds of the local ownership range.
   integer(kind = 4), intent(in), dimension(3) :: ibeg
!> @brief Upper bounds of the local ownership range.
   integer(kind = 4), intent(in), dimension(3) :: iend
!> @brief Input array in \f$(z, xy)\f$ storage.
   real (kind = 8), intent(in) :: F(0:(iend(3) - ibeg(3) + 1) - 1, &
         0:(iend(1) - ibeg(1) + 1)*(iend(2) - ibeg(2) + 1) - 1)
!> @brief Output array in \f$(x, yz)\f$ storage.
   real (kind = 8), intent(out) :: F2(0:(iend(1) - ibeg(1) + 1) - 1, &
         0:(iend(2) - ibeg(2) + 1)*(iend(3) - ibeg(3) + 1) - 1)
!> @brief Loop counters over local tensor-product indices.
   integer(kind = 4) :: ix, iy, iz
!> @brief Flattened local indices in the input and output layouts.
   integer(kind = 4) :: ind1, ind23, ind3, ind12
!> @brief Local extents in the first and second physical directions.
   integer(kind = 4) :: n1, n2

   n1 = iend(1) - ibeg(1) + 1
   n2 = iend(2) - ibeg(2) + 1

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
!$OMP PRIVATE(ix,iy,iz,ind1,ind23,ind3,ind12) SCHEDULE(STATIC)
   do iz = ibeg(3), iend(3)
      do iy = ibeg(2), iend(2)
         do ix = ibeg(1), iend(1)

            ind1 = ix - ibeg(1)
            ind23 = (iy - ibeg(2))+(iz - ibeg(3))*n2

            ind3 = iz - ibeg(3)
            ind12 = (ix - ibeg(1))+(iy - ibeg(2))*n1

            F2(ind1, ind23) = F(ind3, ind12)

         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

end subroutine ReorderRHSForX

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a local trial-space buffer back to canonical x-y-z
!> storage.
!>
!> @details
!> Multi-step ADS substeps may start from y-x-z or z-x-y storage so that
!> the existing directional reorder kernels can be reused. This helper
!> rewrites the final two-dimensional buffer to the canonical
!> `(x, y + z*s_y)` layout expected by history reconstruction and output.
!>
!> The entries of \p order identify the physical axes stored as the first,
!> second-fastest, and third-fastest input indices. For example, the order
!> `(2, 1, 3)` describes a \f$(y, xz)\f$ input buffer.
!
!---------------------------------------------------------------------------
subroutine NormalizeTrialBufferToXYZ(ads, order, F)
   use Setup, ONLY: ADS_Setup
   implicit none
!> @brief Trial-space setup containing the three local buffer extents.
   type(ADS_setup), intent(in) :: ads
!> @brief Physical-axis order represented by the input buffer.
   integer(kind=4), dimension(3), intent(in) :: order
!> @brief Input buffer, replaced with the canonical \f$(x, yz)\f$ layout.
   real(kind=8), allocatable, dimension(:, :), intent(inout) :: F
!> @brief Temporary canonical buffer.
   real(kind=8), allocatable, dimension(:, :) :: Fxyz
!> @brief Current zero-based indices in canonical physical-axis order.
   integer(kind=4), dimension(3) :: idx
!> @brief Loop indices in canonical physical-axis order.
   integer(kind=4) :: ix, iy, iz
!> @brief Flattened input and output indices.
   integer(kind=4) :: in1, in23, out23

   if (order(1) == 1 .and. order(2) == 2 .and. order(3) == 3) return

   allocate (Fxyz(ads%s(1), ads%s(2)*ads%s(3)))

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
!$OMP PRIVATE(ix,iy,iz,idx,in1,in23,out23) SCHEDULE(STATIC)
   do iz = 0, ads%s(3) - 1
      do iy = 0, ads%s(2) - 1
         do ix = 0, ads%s(1) - 1
            idx = (/ix, iy, iz/)
            in1 = idx(order(1)) + 1
            in23 = idx(order(2)) + 1 + idx(order(3))*ads%s(order(2))
            out23 = iy + 1 + iz*ads%s(2)
            Fxyz(ix + 1, out23) = F(in1, in23)
         end do
      end do
   end do
!$OMP END PARALLEL DO

   call move_alloc(Fxyz, F)

end subroutine NormalizeTrialBufferToXYZ

end module reorderRHS
