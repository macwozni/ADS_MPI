!------------------------------------------------------------------------------
!
! MODULE: utils
!
! DESCRIPTION:
!> @file utils.F90
!> @brief Module providing auxiliary utility procedures.
!>
!> @details
!> This module groups a small set of helper procedures used across the
!> code base. The available routines currently include:
!> - a placeholder procedure \ref NormL2,
!> - a distributed assembly-oriented routine \ref Norm_L2,
!> - an integer-to-string conversion utility \ref int2str.
!>
!> The procedures operate on numerical data structures and basic
!> character buffers and are intended to support higher-level
!> computational kernels.
!
!------------------------------------------------------------------------------
module utils

   implicit none

contains


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Placeholder routine for \f$L_2\f$-norm-related setup handling.
!>
!> @details
!> This procedure currently imports selected symbols from module
!> `Setup`, but it does not perform any computational work yet.
!> It may serve as a future entry point for norm evaluation logic
!> or setup-dependent pre-processing.
!
! Notes:
! ------
!> @note
!> At present, this routine is an empty stub.
!
!> @warning
!> The procedure has no effect in its current form.
!
!---------------------------------------------------------------------------
subroutine NormL2 ()
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   implicit none

end subroutine NormL2

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Assembles local contributions associated with a distributed
!> \f$L_2\f$-norm-related quantity.
!>
!> @details
!> This procedure prepares one-dimensional basis-function data for
!> the three parametric directions, traverses the local element range
!> assigned to the current process, and accumulates contributions into
!> the output array \p F.
!>
!> The routine uses knot vectors, spline orders, problem dimensions,
!> local index extents, and process-grid metadata to map global basis
!> indices to the local storage layout. Quadrature-related weights and
!> Jacobians are combined with basis values during the assembly stage.
!>
!> The storage layout of \p F corresponds to a tensor-product block
!> distributed over the first direction explicitly and over the second
!> and third directions through a flattened index.
!
! Input:
! ------
!> @param[in] Ux
!> Knot vector in the first parametric direction.
!>
!> @param[in] Uy
!> Knot vector in the second parametric direction.
!>
!> @param[in] Uz
!> Knot vector in the third parametric direction.
!>
!> @param[in] p
!> Polynomial orders in the three parametric directions.
!>
!> @param[in] n
!> Numbers of basis-function intervals in the three parametric
!> directions. The global problem size is based on
!> \f$(n_x+1)(n_y+1)(n_z+1)\f$.
!>
!> @param[in] nelem
!> Numbers of elements in the three parametric directions.
!>
!> @param[in] ibeg
!> Lower bounds of the locally owned index range in each direction.
!>
!> @param[in] iend
!> Upper bounds of the locally owned index range in each direction.
!>
!> @param[in] nrank
!> Coordinates of the current process in the logical process grid.
!>
!> @param[in] nrp
!> Numbers of processes in the three parametric directions.
!
! Output:
! -------
!> @param[out] F
!> Local array containing accumulated contributions for the owned
!> degrees of freedom.
!
! Notes:
! ------
!> @note
!> The procedure relies on external routines \ref BasisData and
!> \ref global2local to prepare basis information and convert index
!> numbering, respectively.
!
!> @warning
!> The procedure assumes that all array extents and local ownership
!> ranges are mutually consistent.
!
!---------------------------------------------------------------------------
subroutine Norm_L2( &
   Ux, Uy, Uz, &
   p, n, nelem, &
   ibeg, iend, nrank, nrp, &
   F)
   use basis, ONLY: BasisData
   use projection_engine, ONLY: global2local
   implicit none
!> @brief Numbers of basis-function intervals in each direction.
   integer(kind=4), dimension(3), intent(in) :: n
!> @brief Polynomial orders in each direction.
   integer(kind=4), dimension(3), intent(in) :: p
!> @brief Numbers of elements in each direction.
   integer(kind=4), dimension(3), intent(in) :: nelem
!> @brief Local starting indices in each direction.
   integer(kind=4), dimension(3), intent(in) :: ibeg
!> @brief Local ending indices in each direction.
   integer(kind=4), dimension(3), intent(in) :: iend
!> @brief Coordinates of the current process in the process grid.
   integer(kind=4), dimension(3), intent(in) :: nrank
!> @brief Numbers of processes in the process grid directions.
   integer(kind=4), dimension(3), intent(in) :: nrp
!> @brief Knot vector in the first parametric direction.
   real(kind=8), dimension(0:n(1) + p(1) + 1), intent(in) :: Ux
!> @brief Knot vector in the second parametric direction.
   real(kind=8), dimension(0:n(2) + p(2) + 1), intent(in) :: Uy
!> @brief Knot vector in the third parametric direction.
   real(kind=8), dimension(0:n(3) + p(3) + 1), intent(in) :: Uz
!> @brief Output array of locally assembled contributions.
   real(kind=8), intent(out) :: F(0:(iend(1) - ibeg(1) + 1) - 1, &
                                    0:(iend(2) - ibeg(2) + 1)*(iend(3) - ibeg(3) + 1) - 1)

   integer(kind=4) :: mx, my, mz, ngx, ngy, ngz, ex, ey, ez
   integer(kind=4) :: kx, ky, kz, ax, ay, az, d
   integer(kind=4), dimension(nelem(1)) :: Ox
   integer(kind=4), dimension(nelem(2)) :: Oy
   integer(kind=4), dimension(nelem(3)) :: Oz
   real(kind=8), dimension(nelem(1)) :: Jx
   real(kind=8), dimension(nelem(2)) :: Jy
   real(kind=8), dimension(nelem(3)) :: Jz
   real(kind=8), dimension(p(1) + 1) :: Wx
   real(kind=8), dimension(p(2) + 1) :: Wy
   real(kind=8), dimension(p(3) + 1) :: Wz
   real(kind=8), dimension(p(1) + 1, nelem(1)) :: Xx
   real(kind=8), dimension(p(2) + 1, nelem(2)) :: Xy
   real(kind=8), dimension(p(3) + 1, nelem(3)) :: Xz
   real(kind=8), dimension(0:0, 0:p(1), p(1) + 1, nelem(1)) :: NNx
   real(kind=8), dimension(0:0, 0:p(2), p(2) + 1, nelem(2)) :: NNy
   real(kind=8), dimension(0:0, 0:p(3), p(3) + 1, nelem(3)) :: NNz
   real(kind=8) :: J, W, value
   integer(kind=4) :: nreppx, nreppy, nreppz !# elements per proc along x,y,z
   integer(kind=4) :: ind, ind1, ind23, indx, indy, indz
   integer(kind=4) :: iprint

   iprint = 0

   d = 0
   mx = n(1) + p(1) + 1
   ngx = p(1) + 1
   my = n(2) + p(2) + 1
   ngy = p(2) + 1
   mz = n(3) + p(3) + 1
   ngz = p(3) + 1

   call BasisData(p(1), mx, Ux, 0, ngx, nelem(1), Ox, Jx, Wx, Xx, NNx)
   call BasisData(p(2), my, Uy, 0, ngy, nelem(2), Oy, Jy, Wy, Xy, NNy)
   call BasisData(p(3), mz, Uz, 0, ngz, nelem(3), Oz, Jz, Wz, Xz, NNz)

   ! parallel number of elements per processors
   nreppx = nelem(1)/nrp(1)
   nreppy = nelem(2)/nrp(2)
   nreppz = nelem(3)/nrp(3)
   F = 0
   do ex = max(nreppx*nrank(1) - p(1) + 1, 1), min(nelem(1), nreppx*(nrank(1) + 1) + p(1))
      do ey = max(nreppy*nrank(2) - p(2) + 1, 1), min(nelem(2), nreppy*(nrank(2) + 1) + p(2))
         do ez = max(nreppz*nrank(3) - p(3) + 1, 1), min(nelem(3), nreppz*(nrank(3) + 1) + p(3))
            J = Jx(ex)*Jy(ey)*Jz(ez)
            do kx = 1, ngx
               do ky = 1, ngy
                  do kz = 1, ngz
                     W = Wx(kx)*Wy(ky)*Wz(kz)
                     do ax = 0, p(1)
                        do ay = 0, p(2)
                           do az = 0, p(3)
                              d = (Ox(ex) + ax) + (Oy(ey) + ay)*(n(1) + 1) + (Oz(ez) + az)*(n(2) + 1)*(n(1) + 1)
                              call global2local(ind, [n(1), n(2), n(3)], indx, indy, indz)
                              if (indx < ibeg(1) - 1 .or. indx > iend(1) - 1) cycle
                              if (indy < ibeg(2) - 1 .or. indy > iend(2) - 1) cycle
                              if (indz < ibeg(3) - 1 .or. indz > iend(3) - 1) cycle
                              ind1 = indx - ibeg(1) + 1
                              ind23 = (indy - ibeg(2) + 1) + (indz - ibeg(3) + 1)*(iend(2) - ibeg(2) + 1)

                              if (ind1 < 0 .or. ind1 > (iend(1) - ibeg(1))) cycle
                              if (ind23 < 0 .or. ind23 > (iend(2) - ibeg(2) + 1)*(iend(3) - ibeg(3) + 1) - 1) cycle

                              ! parallel
                              F(ind1, ind23) = F(ind1, ind23) + &
                                                NNx(0, ax, kx, ex)*NNy(0, ay, ky, ey)*NNz(0, az, kz, ez)*J*W*value
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do

end subroutine Norm_L2


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts an integer value to a character string.
!>
!> @details
!> The procedure writes the input integer to an internal character
!> buffer, left-adjusts the result, and copies it to the output
!> argument \p str.
!>
!> If the target buffer is shorter than the trimmed textual
!> representation of the integer, a warning message is written to
!> `ERROR_UNIT`. The output string is still assigned, subject to
!> normal Fortran character truncation rules.
!
! Input:
! ------
!> @param[in] n
!> Integer value to be converted.
!
! Output:
! -------
!> @param[out] str
!> Character variable receiving the textual representation of \p n.
!
! Notes:
! ------
!> @note
!> The internal temporary buffer is sized for default 32-bit integers.
!
!> @warning
!> When \p str is too short, the textual representation may be
!> truncated after assignment.
!
!---------------------------------------------------------------------------
subroutine int2str(n, str)
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
   implicit none

!> @brief Integer to be converted.
   integer(kind=4), intent(in) :: n
!> @brief String containing the converted integer.
   character(len=*), intent(out) :: str
!> @brief Temporary buffer sized for the longest signed 32-bit integer. Longest is -2147483647
   character(len=11) :: longstr

   write (longstr, '(I11)') n
   longstr = adjustl(longstr)

   if (len_trim(longstr) > len(str)) then
      write (ERROR_UNIT, '(A,I3,A)') 'int2str: WARNING: can''t fit '//trim(longstr)// &
         ' into a ', len(str), '-character variable'
   end if

   str = longstr

end subroutine int2str

end module utils

