!------------------------------------------------------------------------------
!
! MODULE: rhs_assembly
!
! DESCRIPTION:
!> @file rhs_assembly.F90
!> @brief Three-dimensional ADS right-hand-side assembly routines.
!>
!> @details
!> This module owns the three-dimensional quadrature assembly used to form
!> ADS right-hand-side buffers before the directional solves. It
!> reconstructs the active mixed space, evaluates pointwise contributions,
!> and accumulates deterministic element-local RHS data into the storage
!> layout expected by the subsequent solve stage.
!>
!> The module is shared by the standard ADS path and the iGRM path. When a
!> directional enrichment is active, \ref create_mixed_space selects the
!> enriched test-space metadata for that direction while preserving the
!> trial-space data in the remaining directions.
!>
!> Problem-specific weak-form behavior enters through the callback
!> interfaces from \ref Interfaces; the default pointwise implementation
!> remains in \ref RHS_eq.
!
!------------------------------------------------------------------------------
module rhs_assembly

   implicit none

   private
   public :: Form3DRHS
   public :: create_mixed_space

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Forms the local three-dimensional right-hand-side array for an
!> ADS substep.
!>
!> @details
!> This procedure assembles contributions to the right-hand side of a
!> three-dimensional tensor-product problem. Depending on the enrichment
!> encoded in \p direction, the routine first constructs a mixed test-
!> trial space by calling \ref create_mixed_space.
!>
!> The assembly then proceeds by:
!> - iterating over local elements and quadrature points,
!> - extracting solution values and derivatives from previously computed
!>   buffers,
!> - evaluating the pointwise contribution through
!>   `ComputePointForRHS`,
!> - accumulating the result into a temporary element-local array,
!> - scattering the local contributions into either \p ads_data%F or
!>   \p ads_data%Ft depending on whether iGRM mode is active.
!>
!> The routine uses the quadrature and basis tables already stored in
!> the setup structure and assumes they are allocated consistently.
!
! Input:
! ------
!> @param[in] ads_test
!> Setup structure describing the test space.
!>
!> @param[in] ads_trial
!> Setup structure describing the trial space.
!>
!> @param[in,out] ads_data
!> Working data structure holding buffers used during assembly.
!>
!> @param[in] direction
!> Directional indicator of the currently enriched dimension.
!>
!> @param[in] n
!> Index or number associated with previous time-step data.
!>
!> @param[in] substep
!> Number of the current substep.
!>
!> @param[in] alpha_step
!> Coefficient table used by the substep formula.
!>
!> @param[in] forcing
!> External forcing callback used during pointwise evaluation.
!
! Output:
! -------
!> @param[out] igrm
!> Logical flag indicating whether the current substep uses iGRM.
!
! Notes:
! ------
!> @note
!> The routine allocates a temporary array `elarr` storing the
!> element-local accumulated contributions.
!
!> @warning
!> The procedure assumes that all arrays inside \p ads_test,
!> \p ads_trial, and \p ads_data are dimensionally consistent.
!
!---------------------------------------------------------------------------
subroutine Form3DRHS(ads_test, ads_trial, ads_data, direction, n, substep, &
   alpha_step, forcing, igrm, rhs_point)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   ! use parallelism, ONLY: PRINTRANK
   use Interfaces, ONLY: forcing_fun, rhs_point_fun
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
   use omp_lib
   use RHS_eq
   implicit none
!> @brief Forcing term callback used in pointwise RHS evaluation.
   procedure(forcing_fun) :: forcing
!> @brief Optional callback overriding the default pointwise RHS integrand.
   procedure(rhs_point_fun), optional :: rhs_point
!> @brief Setup structures of the test and trial spaces.
   type(ADS_setup), intent(in) :: ads_test, ads_trial
!> @brief Enrichment indicator for the current substep.
   integer(kind=4), dimension(3), intent(in) :: direction
!> @brief Current substep number.
   integer(kind=4), intent(in) :: substep
!> @brief Working data buffers updated during assembly.
   type(ADS_compute_data), intent(inout) :: ads_data
!> @brief Index associated with previous time-step data.
   integer(kind=4), intent(in) :: n
!> @brief Coefficients of the time-stepping substep formula.
   real(kind=8), intent(in), dimension(7, 3) :: alpha_step
!> @brief Loop counters over quadrature points, local basis indices, and elements.
   integer(kind=4) :: kx, ky, kz, ax, ay, az, ex, ey, ez!, exx,eyy,ezz
!> @brief Linearized local element and element-basis indices.
   integer(kind=4) :: elidx, ldof
!> @brief Numbers of local elements and local basis functions.
   integer(kind=4) :: local_elements, local_dofs
!> @brief Element Jacobian product and quadrature weight product.
   real(kind=8) :: J, W
!> @brief Global and local index variables.
   integer(kind=4) :: ind, ind1, ind23, indx, indy, indz
!> @brief Local indices into the state buffers stored in ads_data.
   integer(kind=4) :: statex, statey, statez
!> @brief Pointwise contribution returned by `ComputePointForRHS`.
   real(kind=8) :: resvalue
!> @brief Physical coordinates of the current quadrature point.
   real(kind=8), dimension(3) :: X
!> @brief Direction-dependent quadrature, element, and local basis indices.
   integer(kind=4), dimension(3) :: k, e, a, indb
!> @brief Local element counts in each direction.
   integer(kind=4), dimension(3) :: lnelem
   ! integer (kind = 4) :: tmp, all
   ! integer (kind = 4) :: total_size
!> @brief Values of the gradient-like quantity from the previous step.
   real(kind=8), dimension(3)  :: du
   ! integer(kind = 4) :: indbx, indby, indbz
!> @brief Solution value at the current quadrature point.
   real(kind=8) :: Uval
!> @brief Auxiliary intermediate solution value.
   real(kind=8) :: Uval13
!> @brief Auxiliary intermediate solution value.
   real(kind=8) :: Uval23
!> @brief Element-local temporary accumulation array.
   real(kind=8), dimension(:, :, :), allocatable :: elarr
!> @brief Element-wise RHS contributions before deterministic global scatter.
   real(kind=8), dimension(:, :), allocatable :: element_rhs
!> @brief Mixed-space setup used for the current assembly pass.
   type(ADS_setup) :: ads
!> @brief Flag indicating use of the iGRM path.
   logical, intent(out) :: igrm
!> @brief Permutation of coordinate directions used internally.
   integer(kind=4) :: dira,dirb,dirc

   call create_mixed_space(ads_test, ads_trial, direction,&
   ads, dira, dirb, dirc, igrm)

   ! Each OpenMP thread needs a private element-local accumulator.
   lnelem = ads%maxe - ads%mine + 1
   local_elements = lnelem(1)*lnelem(2)*lnelem(3)
   local_dofs = (ads%p(dira) + 1)*(ads%p(dirb) + 1)*(ads%p(dirc) + 1)
   allocate (element_rhs(local_dofs, local_elements))
   element_rhs = 0.d0

!   if (allocated(ads_data%F)) ads_data%F = 0.d0
!   if (allocated(ads_data%Ft)) ads_data%Ft = 0.d0

!      loop over points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ex,ey,ez,e,kx,ky,kz,k,W,ax,ay,az,a,ind,indx,indy,indz,ind1,ind23,J) &
!$OMP PRIVATE(statex,statey,statez,X,du,resvalue,indb,Uval,Uval13,Uval23,elarr,elidx,ldof)
   allocate (elarr(0:ads%p(dira), 0:ads%p(dirb), 0:ads%p(dirc)))
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
   ! do all = 1, total_size
! translate coefficients to local
   ! ez = modulo(all - 1, ads % lnelem(3))
   ! tmp = (all - ez)/ads % lnelem(3) + 1
   ! ey = modulo(tmp - 1, ads % lnelem(2))
   ! ex = (tmp - ey)/ads % lnelem(2)
   ! write(*,*) size(ads%Jx) , ads % lnelem(1), ads % mine(1)
   ! write(*,*) size(ads%Jy) , ads % lnelem(2), ads % mine(2)
   ! write(*,*) size(ads%Jz) , ads % lnelem(3), ads % mine(3)
   ! do exx=1,ads % lnelem(1)
   ! do eyy=1,ads % lnelem(2)
   ! do ezz=1,ads % lnelem(3)
   do ex = ads%mine(1), ads%maxe(1)
      do ey = ads%mine(2), ads%maxe(2)
         do ez = ads%mine(3), ads%maxe(3)
! fix distributed part
            ! ex = exx + ads % mine(1)
            ! ey = eyy + ads % mine(2)
            ! ez = ezz + ads % mine(3)
! Jacobian
            J = ads%Jx(ex)*ads%Jy(ey)*ads%Jz(ez)
            e = (/ex, ey, ez/)
            elidx = ((ex - ads%mine(1))*lnelem(2) + &
               (ey - ads%mine(2)))*lnelem(3) + &
               (ez - ads%mine(3)) + 1
            elarr = 0.d0
! loop over quadrature points
            do kx = 1, ads%ng(dira)
               do ky = 1, ads%ng(dirb)
                  do kz = 1, ads%ng(dirc)
                     k(dira) = kx
                     k(dirb) = ky
                     k(dirc) = kz
! weigths
                     W = ads%Wx(k(1))*ads%Wy(k(2))*ads%Wz(k(3))
                     statex = ex - ads_data%state_mine(1) + 1
                     statey = ey - ads_data%state_mine(2) + 1
                     statez = ez - ads_data%state_mine(3) + 1
                     Uval = ads_data%Un(statex, statey, statez, k(1), k(2), k(3))
                     Uval13 = ads_data%Un13(statex, statey, statez, k(1), k(2), k(3))
                     Uval23 = ads_data%Un23(statex, statey, statez, k(1), k(2), k(3))
                     du = ads_data%dUn(statex, statey, statez, k(1), k(2), k(3), :)

!                 loop over degrees of freedom
                     do ax = 0, ads%p(dira)
                        do ay = 0, ads%p(dirb)
                           do az = 0, ads%p(dirc)
                              a(dira) = ax
                              a(dirb) = ay
                              a(dirc) = az

                              indb(1) = ads%Ox(ex) + a(1)
                              indb(2) = ads%Oy(ey) + a(2)
                              indb(3) = ads%Oz(ez) + a(3)

                              indx = indb(1)
                              indy = indb(2)
                              indz = indb(3)
                              ind = indx + (indy + indz*(ads%n(2) + 1))*(ads%n(1) + 1)

                              if ((indb(1) < ads%ibeg(1) - 1) .or. (indb(1) > ads%iend(1) - 1) .or. &
                                 (indb(2) < ads%ibeg(2) - 1) .or. (indb(2) > ads%iend(2) - 1) .or. &
                                 (indb(3) < ads%ibeg(3) - 1) .or. (indb(3) > ads%iend(3) - 1)) then
                              else
                                 ind1 = indb(dira) - ads%ibeg(dira) + 1
                                 ind23 = (indb(dirb) - ads%ibeg(dirb) + 1) + &
                                    (indb(dirc) - ads%ibeg(dirc) + 1)*(ads%iend(dirb) - ads%ibeg(dirb) + 1)

                                 X(1) = ads%Xx(k(1), ex)
                                 X(2) = ads%Xy(k(2), ey)
                                 X(3) = ads%Xz(k(3), ez)


                                 ! call RHS_fun(&
                                 ! ads, &
                                 ! X, &
                                 ! k, &
                                 ! e, &
                                 ! a, &
                                 ! du, &
                                 ! 1, Uval_m, Uval13,Uval23, &
                                 ! ads_data, J, W, direction, substep, resvalue)

                                 if (present(rhs_point)) then
                                    call rhs_point( &
                                       ads, &
                                       X, &
                                       k, &
                                       e, &
                                       a, &
                                       du, &
                                       n, &
                                       Uval, &
                                       Uval13, &
                                       Uval23, &
                                       ads_data, J, W, direction, substep, &
                                       alpha_step, &
                                       forcing, &
                                       resvalue)
                                 else
                                    call ComputePointForRHS( &
                                       ads, &
                                       X, &
                                       k, &
                                       e, &
                                       a, &
                                       du, &
                                       n, &
                                       Uval, &
                                       Uval13, &
                                       Uval23, &
                                       ads_data, J, W, direction, substep, &
                                       alpha_step, &
                                       forcing, &
                                       resvalue)
                                 end if

                                 elarr(ax, ay, az) = elarr(ax, ay, az) + resvalue
                              end if
                           end do
                        end do
                     end do
                  end do
               end do
            end do
            do ax = 0, ads%p(dira)
               do ay = 0, ads%p(dirb)
                  do az = 0, ads%p(dirc)
                     ldof = ax + (ads%p(dira) + 1)*(ay + (ads%p(dirb) + 1)*az) + 1
                     element_rhs(ldof, elidx) = elarr(ax, ay, az)
                  end do
               end do
            end do
         end do
      end do
   end do
!$OMP END DO
   deallocate (elarr)
!$OMP END PARALLEL

! moving results from temporary array to main one in deterministic element order
   do ex = ads%mine(1), ads%maxe(1)
      do ey = ads%mine(2), ads%maxe(2)
         do ez = ads%mine(3), ads%maxe(3)
            elidx = ((ex - ads%mine(1))*lnelem(2) + &
               (ey - ads%mine(2)))*lnelem(3) + &
               (ez - ads%mine(3)) + 1
            do ax = 0, ads%p(dira)
               do ay = 0, ads%p(dirb)
                  do az = 0, ads%p(dirc)
                     a(dira) = ax
                     a(dirb) = ay
                     a(dirc) = az

                     indb(1) = ads%Ox(ex) + a(1)
                     indb(2) = ads%Oy(ey) + a(2)
                     indb(3) = ads%Oz(ez) + a(3)

                     if ((indb(1) < ads%ibeg(1) - 1) .or. (indb(1) > ads%iend(1) - 1) .or. &
                        (indb(2) < ads%ibeg(2) - 1) .or. (indb(2) > ads%iend(2) - 1) .or. &
                        (indb(3) < ads%ibeg(3) - 1) .or. (indb(3) > ads%iend(3) - 1)) then
                     else
                        ind1 = indb(dira) - ads%ibeg(dira) + 1
                        ind23 = (indb(dirb) - ads%ibeg(dirb) + 1) + &
                           (indb(dirc) - ads%ibeg(dirc) + 1)*(ads%iend(dirb) - ads%ibeg(dirb) + 1)
                        ldof = ax + (ads%p(dira) + 1)*(ay + (ads%p(dirb) + 1)*az) + 1

                        if (igrm) then
                           ads_data%Ft(ind1 + 1, ind23 + 1) = &
                              ads_data%Ft(ind1 + 1, ind23 + 1) &
                              + element_rhs(ldof, elidx)
                        else
                           ads_data%F(ind1 + 1, ind23 + 1) = &
                              ads_data%F(ind1 + 1, ind23 + 1) &
                              + element_rhs(ldof, elidx)
                        end if
                     end if
                  end do
               end do
            end do
         end do
      end do
   end do
   deallocate (element_rhs)

end subroutine Form3DRHS

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Creates the mixed space used for a directional iGRM substep.
!>
!> @details
!> This routine starts from the trial-space setup and selectively
!> replaces one coordinate direction by the corresponding data from the
!> test-space setup. The selected direction is identified by the vector
!> \p direction.
!>
!> Depending on the enriched direction, the routine also sets:
!> - the permutation \p dira, \p dirb, \p dirc controlling the ordering
!>   of tensor loops,
!> - the flag \p igrm indicating whether an iGRM configuration is active.
!>
!> If no direction is marked for enrichment, the output setup remains
!> identical to the trial space and \p igrm is returned as `.FALSE.`.
!
! Input:
! ------
!> @param[in] ads_test
!> Setup structure of the test space.
!>
!> @param[in] ads_trial
!> Setup structure of the trial space.
!>
!> @param[in] direction
!> Directional enrichment selector.
!
! Output:
! -------
!> @param[out] ads
!> Resulting mixed-space setup.
!>
!> @param[out] dira
!> First direction used in tensor traversal.
!>
!> @param[out] dirb
!> Second direction used in tensor traversal.
!>
!> @param[out] dirc
!> Third direction used in tensor traversal.
!>
!> @param[out] igrm
!> Logical flag indicating whether iGRM is active.
!
!---------------------------------------------------------------------------
subroutine create_mixed_space(ads_test, ads_trial, direction,&
   ads, dira, dirb, dirc, igrm)
   use Setup, ONLY: ADS_Setup
   implicit none
!> @brief Setup structures of the test and trial spaces.
   type(ADS_setup), intent(in) :: ads_test, ads_trial
!> @brief Directional enrichment selector.
   integer(kind=4), dimension(3), intent(in) :: direction
!> @brief Output mixed-space setup.
   type(ADS_setup), intent(out) :: ads
!> @brief Internal permutation of traversal directions.
   integer(kind=4), intent(out) :: dira,dirb,dirc
!> @brief Flag indicating use of the iGRM configuration.
   logical, intent(out) :: igrm

   dira=1
   dirb=2
   dirc=3

!  copy default space as trial space
   ads = ads_trial
   igrm = .FALSE.

!  if we have enriched one direction, then modify default space
   if (direction(1) .EQ. 1) then
      ads%n(1) = ads_test%n(1)
      ads%p(1) = ads_test%p(1)
      ads%Ux = ads_test%Ux
      ads%nelem(1) = ads_test%nelem(1)
      ads%dimensionsX = ads_test%dimensionsX
      ads%shiftsX = ads_test%shiftsX
      ads%IPIVx = ads_test%IPIVx
      ads%nrcpp(1) = ads_test%nrcpp(1)
      ads%ibeg(1) = ads_test%ibeg(1)
      ads%iend(1) = ads_test%iend(1)
      ads%s(1) = ads_test%s(1)
      ads%ibegsx = ads_test%ibegsx
      ads%iendsx = ads_test%iendsx
      ads%mine(1) = ads_test%mine(1)
      ads%maxe(1) = ads_test%maxe(1)
      ads%lnelem(1) = ads_test%lnelem(1)
      ads%m(1) = ads_test%m(1)
      ads%ng(1) = ads_test%ng(1)
      ads%Ox = ads_test%Ox
      ads%Jx = ads_test%Jx
      ads%Xx = ads_test%Xx
      ads%NNx = ads_test%NNx
      ads%Wx = ads_test%Wx
      igrm = .TRUE.
      dira=1
      dirb=2
      dirc=3
   end if
   if (direction(2) .EQ. 1) then
      ads%n(2) = ads_test%n(2)
      ads%p(2) = ads_test%p(2)
      ads%Uy = ads_test%Uy
      ads%nelem(2) = ads_test%nelem(2)
      ads%dimensionsY = ads_test%dimensionsY
      ads%shiftsY = ads_test%shiftsY
      ads%IPIVy = ads_test%IPIVy
      ads%nrcpp(2) = ads_test%nrcpp(2)
      ads%ibeg(2) = ads_test%ibeg(2)
      ads%iend(2) = ads_test%iend(2)
      ads%s(2) = ads_test%s(2)
      ads%ibegsy = ads_test%ibegsy
      ads%iendsy = ads_test%iendsy
      ads%mine(2) = ads_test%mine(2)
      ads%maxe(2) = ads_test%maxe(2)
      ads%lnelem(2) = ads_test%lnelem(2)
      ads%m(2) = ads_test%m(2)
      ads%ng(2) = ads_test%ng(2)
      ads%Oy = ads_test%Oy
      ads%Jy = ads_test%Jy
      ads%Xy = ads_test%Xy
      ads%NNy = ads_test%NNy
      ads%Wy = ads_test%Wy
      igrm = .TRUE.
      dira=2
      dirb=1
      dirc=3
   end if
   if (direction(3) .EQ. 1) then
      ads%n(3) = ads_test%n(3)
      ads%p(3) = ads_test%p(3)
      ads%Uz = ads_test%Uz
      ads%nelem(3) = ads_test%nelem(3)
      ads%dimensionsZ = ads_test%dimensionsZ
      ads%shiftsZ = ads_test%shiftsZ
      ads%IPIVz = ads_test%IPIVz
      ads%nrcpp(3) = ads_test%nrcpp(3)
      ads%ibeg(3) = ads_test%ibeg(3)
      ads%iend(3) = ads_test%iend(3)
      ads%s(3) = ads_test%s(3)
      ads%ibegsz = ads_test%ibegsz
      ads%iendsz = ads_test%iendsz
      ads%mine(3) = ads_test%mine(3)
      ads%maxe(3) = ads_test%maxe(3)
      ads%lnelem(3) = ads_test%lnelem(3)
      ads%m(3) = ads_test%m(3)
      ads%ng(3) = ads_test%ng(3)
      ads%Oz = ads_test%Oz
      ads%Jz = ads_test%Jz
      ads%Xz = ads_test%Xz
      ads%NNz = ads_test%NNz
      ads%Wz = ads_test%Wz
      igrm = .TRUE.
      dira=3
      dirb=1
      dirc=2
   end if

end subroutine create_mixed_space

end module rhs_assembly
