!------------------------------------------------------------------------------
!
! MODULE: solution_reconstruction
!
! DESCRIPTION:
!> @file solution_reconstruction.F90
!> @brief Quadrature-point solution reconstruction helpers.
!>
!> @details
!> This module reconstructs solution values and first derivatives at
!> quadrature points from the distributed spline coefficient buffers
!> exchanged between neighbouring MPI ranks.
!>
!> The reconstructed fields populate the time-history arrays stored in
!> \ref ADS_compute_data, including the base state and the intermediate
!> states used by multi-substep iGRM schemes. The active state is selected
!> by the caller and later consumed by RHS assembly.
!>
!> The module also contains the tensor-product index conversion helper
!> historically used by utilities and legacy problem-local sources.
!
!------------------------------------------------------------------------------
module solution_reconstruction

   implicit none

   private
   public :: FormUn
   public :: global2local

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reconstructs solution values and first derivatives at
!> quadrature points from neighbouring-domain coefficient blocks.
!>
!> @details
!> This routine evaluates:
!> - the scalar solution value at each quadrature point,
!> - the first derivatives with respect to the three coordinates,
!> - the corresponding storage buffers in \p ads_data for one of the
!>   solution states selected by \p subun.
!>
!> For each local element and quadrature point, the procedure loops over
!> the locally supported tensor-product basis functions, retrieves the
!> appropriate coefficient from the neighbouring-domain block array
!> \p ads_data%R, and forms the value and derivatives using the basis
!> tables stored in \p ads.
!>
!> The computed value is written to one of:
!> - \p ads_data%Un and \p ads_data%dUn0 for \p subun equal to 1,
!> - \p ads_data%Un13 and \p ads_data%dUn13 for \p subun equal to 2,
!> - \p ads_data%Un23 and \p ads_data%dUn23 for \p subun equal to 3,
!>
!> while the active derivative vector is always also written to
!> \p ads_data%dUn. The selected value and history-derivative buffers,
!> together with \p ads_data%dUn, are cleared before reconstruction.
!
! Input:
! ------
!> @param[in] subun
!> Selector identifying which solution buffer is updated.
!>
!> @param[in] ads
!> ADS setup structure containing basis tables and ownership metadata.
!>
!> @param[in,out] ads_data
!> Working data structure containing coefficient blocks and output
!> arrays.
!
! Notes:
! ------
!> @note
!> The neighbouring-domain block indices are encoded by the triplet
!> \f$(r_x,r_y,r_z)\f$ taking values in \f$\{1,2,3\}\f$ for each
!> direction.
!
!> @warning
!> The procedure assumes that \p ads_data%R and all output arrays have
!> been allocated consistently with the setup structure.
!>
!> @warning
!> The supported basis functions may reach only the current MPI rank and
!> its immediate neighbours represented by the 3-by-3-by-3 block stencil
!> in \p ads_data%R.
!
!---------------------------------------------------------------------------
subroutine FormUn(subun, ads, ads_data)
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
   use Setup, ONLY: ADS_Setup, ADS_compute_data
#ifdef IDEBUG
   use parallelism, ONLY: PRINTRANK
#endif
   use Interfaces, ONLY: forcing_fun
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
   use omp_lib
   implicit none
!> @brief Selector of the solution buffer to be updated.
   integer(kind=4), intent(in) :: subun
!> @brief Setup structure with basis tables and decomposition data.
   type(ADS_setup), intent(in) :: ads
!> @brief Working data structure updated in place.
   type(ADS_compute_data), intent(inout) :: ads_data
!> @brief Loop counters over quadrature points and elements.
   integer(kind=4) :: kx, ky, kz, ex, ey, ez, exx, eyy, ezz
!> @brief Linearized global or local index.
   integer(kind=4) :: ind
!> @brief Auxiliary index variables retained for alternate traversal strategies.
   integer(kind=4) :: tmp, all
!> @brief Total number of element blocks in the local partition.
   integer(kind=4) :: total_size
!> @brief Neighbour-block selectors and local coordinates within these blocks.
   integer(kind=4) :: rx, ry, rz, ix, iy, iz, sx, sy, sz
!> @brief Local basis-function indices in the tensor-product basis.
   integer(kind=4) :: bx, by, bz
!> @brief Reconstructed first derivatives of the solution.
   real(kind=8) :: dux, duy, duz
   ! real   (kind=8), dimension(3)  :: du
!> @brief Global basis-function coordinates in the three directions.
   integer(kind=4) :: indbx, indby, indbz
!> @brief Solution value and current coefficient value.
   real(kind=8) :: Uval, ucoeff
!> @brief Basis products for value and directional derivatives.
   real(kind=8) :: dvx, dvy, dvz, v
!> @brief Local indices into the state buffers stored in ads_data.
   integer(kind=4) :: statex, statey, statez

   select case (subun)
   case (1)
      ads_data%Un = 0.d0
      ads_data%dUn0 = 0.d0
   case (2)
      ads_data%Un13 = 0.d0
      ads_data%dUn13 = 0.d0
   case (3)
      ads_data%Un23 = 0.d0
      ads_data%dUn23 = 0.d0
   case default
      write (ERROR_UNIT, *) "wrong substep"
      stop 1
   end select
   ads_data%dUn = 0.d0
   !total_size = ads%lnelem(1)*ads%lnelem(2)*ads%lnelem(3)

!      loop over points
!$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) &
!$OMP PRIVATE(ex,ey,ez,kx,ky,kz,ind,statex,statey,statez) &
!$OMP PRIVATE(bx,by,bz,rx,ry,rz,ix,iy,iz,sx,sy,sz,Ucoeff,dvx,dvy,dvz) &
!$OMP PRIVATE(indbx,indby,indbz,Uval,dux,duy,duz,v) &
!$OMP SCHEDULE(STATIC)
   !do all = 1, total_size
!        translate coefficients to local
   !ez = modulo(all - 1, ads%lnelem(3))
   !tmp = (all - ez)/ads%lnelem(3) + 1
   !ey = modulo(tmp - 1, ads%lnelem(2))
   !ex = (tmp - ey)/ads%lnelem(2)
   do ex = ads_data%state_mine(1), ads_data%state_maxe(1)
      do ey = ads_data%state_mine(2), ads_data%state_maxe(2)
         do ez = ads_data%state_mine(3), ads_data%state_maxe(3)
            statex = ex - ads_data%state_mine(1) + 1
            statey = ey - ads_data%state_mine(2) + 1
            statez = ez - ads_data%state_mine(3) + 1
!        loop over quadrature points
            do kx = 1, ads%ng(1)
               do ky = 1, ads%ng(2)
                  do kz = 1, ads%ng(3)
                     Uval = 0.d0
                     dux = 0.d0
                     duy = 0.d0
                     duz = 0.d0
!                 compute value of derivative from previous time step - du
!                 compute previous solution coefficient at given point - Uval
                     do bx = 0, ads%p(1)
                        do by = 0, ads%p(2)
                           do bz = 0, ads%p(3)
                              indbx = (ads%Ox(ex) + bx)
                              indby = (ads%Oy(ey) + by)
                              indbz = (ads%Oz(ez) + bz)
                              ind = indbx + (indby + indbz*(ads%n(2) + 1))*(ads%n(1) + 1)

                              rx = 2
                              ry = 2
                              rz = 2
                              if (indbx < ads%ibeg(1) - 1) rx = 1
                              if (indbx > ads%iend(1) - 1) rx = 3
                              if (indby < ads%ibeg(2) - 1) ry = 1
                              if (indby > ads%iend(2) - 1) ry = 3
                              if (indbz < ads%ibeg(3) - 1) rz = 1
                              if (indbz > ads%iend(3) - 1) rz = 3

                              ix = indbx - ads%ibegsx(rx) + 1
                              iy = indby - ads%ibegsy(ry) + 1
                              iz = indbz - ads%ibegsz(rz) + 1
                              sx = ads%iendsx(rx) - ads%ibegsx(rx) + 1
                              sy = ads%iendsy(ry) - ads%ibegsy(ry) + 1
                              sz = ads%iendsz(rz) - ads%ibegsz(rz) + 1
                              ind = ix + sx*(iy + sy*iz)

#ifdef IDEBUG
                              if (ind < 0 .or. ind > ads%nrcpp(3)*ads%nrcpp(1)*ads%nrcpp(2) - 1) then
                                 write (ERROR_UNIT, *) PRINTRANK, 'Oh crap', ix, iy, iz
                                 write (ERROR_UNIT, *) PRINTRANK, 'r', rx, ry, rz
                                 write (ERROR_UNIT, *) PRINTRANK, 'x', ads%ibeg(1), ads%iend(1)
                                 write (ERROR_UNIT, *) PRINTRANK, 'y', ads%ibeg(2), ads%iend(2)
                                 write (ERROR_UNIT, *) PRINTRANK, 'z', ads%ibeg(3), ads%iend(3)
                                 write (ERROR_UNIT, *) PRINTRANK, 'sizes=', sx, sy, sz
                                 write (ERROR_UNIT, *) PRINTRANK, 'begsx=', ads%ibegsx
                                 write (ERROR_UNIT, *) PRINTRANK, 'endsx=', ads%iendsx
                                 write (ERROR_UNIT, *) PRINTRANK, 'begsy=', ads%ibegsy
                                 write (ERROR_UNIT, *) PRINTRANK, 'endsy=', ads%iendsy
                                 write (ERROR_UNIT, *) PRINTRANK, 'begsz=', ads%ibegsz
                                 write (ERROR_UNIT, *) PRINTRANK, 'endsz=', ads%iendsz
                              end if
#endif

                              Ucoeff = ads_data%R(ind + 1, rx, ry, rz)
                              v = ads%NNx(0, bx, kx, ex)*ads%NNy(0, by, ky, ey)*ads%NNz(0, bz, kz, ez)
                              dvx = ads%NNx(1, bx, kx, ex)*ads%NNy(0, by, ky, ey)*ads%NNz(0, bz, kz, ez)
                              dvy = ads%NNx(0, bx, kx, ex)*ads%NNy(1, by, ky, ey)*ads%NNz(0, bz, kz, ez)
                              dvz = ads%NNx(0, bx, kx, ex)*ads%NNy(0, by, ky, ey)*ads%NNz(1, bz, kz, ez)

                              Uval = Uval + Ucoeff*v
                              dux = dux + Ucoeff*dvx
                              duy = duy + Ucoeff*dvy
                              duz = duz + Ucoeff*dvz
                           end do
                        end do
                     end do
                     ads_data%dUn(statex, statey, statez, kx, ky, kz, :) = (/dux, duy, duz/)
                     if (subun .EQ. 1) then
                        ads_data%Un(statex, statey, statez, kx, ky, kz) = Uval
                        ads_data%dUn0(statex, statey, statez, kx, ky, kz, :) = (/dux, duy, duz/)
                     else if (subun .EQ. 2) then
                        ads_data%Un13(statex, statey, statez, kx, ky, kz) = Uval
                        ads_data%dUn13(statex, statey, statez, kx, ky, kz, :) = (/dux, duy, duz/)
                     else if (subun .EQ. 3) then
                        ads_data%Un23(statex, statey, statez, kx, ky, kz) = Uval
                        ads_data%dUn23(statex, statey, statez, kx, ky, kz, :) = (/dux, duy, duz/)
                     else
                        write (ERROR_UNIT, *) "wrong substep"
                     end if
                  end do
               end do
            end do
         end do
      end do
   end do
!$OMP END PARALLEL DO

end subroutine FormUn

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a linearized tensor-product index into Cartesian
!> coordinates.
!>
!> @details
!> This routine maps the one-dimensional index
!>
!> \f[
!> \text{ind} = (z\,(n_y+1) + y)\,(n_x+1) + x
!> \f]
!>
!> to the corresponding coordinate triple \f$(x,y,z)\f$ in a structured
!> tensor-product grid of logical size
!> \f$(n_x+1) \times (n_y+1) \times (n_z+1)\f$.
!
! Input:
! ------
!> @param[in] ind
!> Linearized global index.
!>
!> @param[in] n
!> Sizes of the structured grid minus one in the three directions.
!
! Output:
! -------
!> @param[out] x
!> Coordinate in the first direction.
!>
!> @param[out] y
!> Coordinate in the second direction.
!>
!> @param[out] z
!> Coordinate in the third direction.
!
!---------------------------------------------------------------------------
subroutine global2local(ind, n, x, y, z)
   implicit none
!> @brief Linearized tensor-product index.
   integer(kind=4), intent(in) :: ind
!> @brief Grid sizes minus one in the three directions.
   integer(kind=4), dimension(3), intent(in) :: n
!> @brief Output Cartesian coordinates.
   integer(kind=4), intent(out) :: x, y, z
!> @brief Auxiliary remainder used during index splitting.
   integer(kind=4) :: tmp

   z = ind/((n(1) + 1)*(n(2) + 1))
   tmp = ind - z*(n(1) + 1)*(n(2) + 1)
   y = tmp/(n(1) + 1)
   x = tmp - y*(n(1) + 1)

end subroutine global2local

end module solution_reconstruction
