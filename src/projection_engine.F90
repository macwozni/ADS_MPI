!------------------------------------------------------------------------------
!
! MODULE: projection_engine
!
! DESCRIPTION:
!> @file projection_engine.F90
!> @brief Module providing matrix-assembly and right-hand-side formation
!> routines for projection-related operators.
!>
!> @details
!> This module groups procedures used to:
!> - assemble sparse one-dimensional operator matrices composed of mass,
!>   stiffness, and advection-like contributions,
!> - form three-dimensional right-hand-side arrays for ADS-based
!>   substeps,
!> - construct mixed test-trial spaces used by the iGRM workflow,
!> - evaluate solution values and derivatives at quadrature points from
!>   stored coefficient blocks,
!> - translate linearized tensor-product indices into Cartesian
!>   coordinates.
!>
!> The implementation relies on spline-basis data prepared externally
!> and on sparse-matrix helper routines supplied by module `sparse`.
!>
!------------------------------------------------------------------------------
module projection_engine

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Checks whether two spline spaces use the same geometric mesh.
!>
!> @details
!> The comparison ignores knot multiplicities and only checks the sequence of
!> distinct knot locations. This matches the iGRM requirement that test and
!> trial spaces share the same element partition, while still allowing
!> repeated knots used to lower continuity.
!
!---------------------------------------------------------------------------
logical function SameKnotLocations(U1, p1, n1, U2, p2, n2) result(same)
   implicit none
!> @brief Polynomial degrees and basis sizes minus one.
   integer(kind=4), intent(in) :: p1, n1, p2, n2
!> @brief Knot vectors to compare.
   real(kind=8), intent(in) :: U1(0:n1 + p1 + 1)
   real(kind=8), intent(in) :: U2(0:n2 + p2 + 1)
   real(kind=8), parameter :: knot_tol = 1.d-12
   integer(kind=4) :: i1, i2, last1, last2
   real(kind=8) :: k1, k2

   same = .FALSE.
   i1 = 0
   i2 = 0
   last1 = n1 + p1 + 1
   last2 = n2 + p2 + 1

   do while (i1 <= last1 .and. i2 <= last2)
      k1 = U1(i1)
      k2 = U2(i2)
      if (abs(k1 - k2) > knot_tol) return

      do
         if (i1 > last1) exit
         if (abs(U1(i1) - k1) > knot_tol) exit
         i1 = i1 + 1
      end do
      do
         if (i2 > last2) exit
         if (abs(U2(i2) - k2) > knot_tol) exit
         i2 = i2 + 1
      end do
   end do

   same = (i1 > last1 .and. i2 > last2)

end function SameKnotLocations

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Validates the iGRM mixed-space assumptions before matrix assembly.
!
!---------------------------------------------------------------------------
subroutine ValidateIGRMMesh(U_test, p_test, n_test, nelem_test, &
                            U_trial, p_trial, n_trial, nelem_trial)
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
   implicit none
!> @brief Test-space metadata.
   integer(kind=4), intent(in) :: p_test, n_test, nelem_test
!> @brief Trial-space metadata.
   integer(kind=4), intent(in) :: p_trial, n_trial, nelem_trial
!> @brief Test and trial knot vectors.
   real(kind=8), intent(in) :: U_test(0:n_test + p_test + 1)
   real(kind=8), intent(in) :: U_trial(0:n_trial + p_trial + 1)

   if (p_test <= p_trial) then
      write (ERROR_UNIT, *) 'iGRM requires test degree greater than trial degree:', &
         p_test, p_trial
      stop 1
   end if

   if (.not. SameKnotLocations(U_test, p_test, n_test, U_trial, p_trial, n_trial)) then
      write (ERROR_UNIT, *) 'iGRM requires identical test/trial knot locations'
      stop 1
   end if

   if (nelem_test /= nelem_trial) then
      write (ERROR_UNIT, *) 'iGRM mesh metadata mismatch:', nelem_test, nelem_trial
      stop 1
   end if

end subroutine ValidateIGRMMesh

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Assembles a block sparse matrix composed of three mixed-space
!> submatrices.
!>
!> @details
!> This procedure constructs a sparse matrix associated with two
!> one-dimensional spline spaces. The resulting matrix contains:
!> - a test-test block assembled from \p mixA,
!> - a trial-test coupling block assembled from \p mixB,
!> - a test-trial coupling block assembled from \p mixBT.
!>
!> For each contribution, the routine evaluates combinations of:
!> - the mass term \f$M = u\,v\f$,
!> - the stiffness term \f$K = u'\,v'\f$,
!> - the advection-like term \f$B = u'\,v\f$,
!> - the transposed advection-like term \f$B^T = u\,v'\f$.
!>
!> The basis values and first derivatives are computed by \ref BasisData
!> on the test and trial spaces separately. The assembled entries are
!> inserted into a sparse matrix through routine `add`.
!>
!> The total matrix size is \f$(n_1+n_2+2) \times (n_1+n_2+2)\f$.
!
! Input:
! ------
!> @param[in] nelem
!> Number of elements used during assembly.
!>
!> @param[in] U1
!> Knot vector of the first spline space.
!>
!> @param[in] p1
!> Polynomial degree of the first spline space.
!>
!> @param[in] n1
!> Number of control points minus one for the first spline space.
!>
!> @param[in] U2
!> Knot vector of the second spline space.
!>
!> @param[in] p2
!> Polynomial degree of the second spline space.
!>
!> @param[in] n2
!> Number of control points minus one for the second spline space.
!>
!> @param[in] mixA
!> Mixing coefficients for the first diagonal block.
!>
!> @param[in] mixB
!> Mixing coefficients for the off-diagonal coupling block.
!>
!> @param[in] mixBT
!> Mixing coefficients for the transposed off-diagonal coupling block.
!
! Output:
! -------
!> @param[out] sprsmtrx
!> Sparse matrix containing the assembled block operator.
!
! Notes:
! ------
!> @note
!> The routine requests first derivatives only, therefore the derivative
!> order passed to \ref BasisData is equal to one.
!
!> @warning
!> The argument \p nelem is used for both spline spaces and is assumed to
!> be consistent with both knot vectors.
!
!---------------------------------------------------------------------------
subroutine MKBBT_large(nelem, U1, p1, n1, U2, p2, n2, mixA, mixB, mixBT, sprsmtrx)
   use basis, ONLY: BasisData
   use omp_lib
   use sparse
   implicit none
   !> @brief Number of elements used during assembly.
   integer(kind=4), intent(in) :: nelem
   !> @brief Basis size minus one and polynomial degree for the first space.
   integer(kind=4), intent(in) :: n1, p1
   !> @brief Basis size minus one and polynomial degree for the second space.
   integer(kind=4), intent(in) :: n2, p2
   !> @brief Knot vector of the first space.
   real(kind=8), intent(in) :: U1(0:n1 + p1 + 1)
!> @brief Knot vector of the second space.
   real(kind=8), intent(in) :: U2(0:n2 + p2 + 1)
!> @brief Mixing coefficients for the assembled operator terms.
   real(kind=8), dimension(4), intent(in) :: mixA, mixB, mixBT
!> @brief Element Jacobians.
   real(kind=8), dimension(nelem) :: J
!> @brief Gauss weights for the quadrature rule.
   real(kind=8), dimension(p1 + 1) :: W
!> @brief Physical Gauss-point coordinates on all elements.
   real(kind=8), dimension(p1 + 1, nelem) :: X
!> @brief Basis values and first derivatives for the first space.
   real(kind=8), dimension(0:1, 0:p1, p1 + 1, nelem) :: NN1
!> @brief Basis values and first derivatives for the second space.
   real(kind=8), dimension(0:1, 0:p2, p1 + 1, nelem) :: NN2
!> @brief Highest derivative order requested from basis evaluation.
   integer(kind=4) :: dd
!> @brief Row and column indices of the assembled sparse entry.
   integer(kind=4) :: ia, ib
!> @brief Last valid knot indices in the two knot vectors.
   integer(kind=4) :: mm1, mm2
!> @brief Number of Gauss quadrature points.
   integer(kind=4) :: ng
!> @brief Loop counters over elements, points, and local basis indices.
   integer(kind=4) :: e, i, c, d
!> @brief First nonzero basis-function indices for the first space.
   integer(kind=4) :: O1(nelem)
!> @brief First nonzero basis-function indices for the second space.
   integer(kind=4) :: O2(nelem)
!> @brief Output sparse matrix.
   type(sparse_matrix), pointer, intent(out) :: sprsmtrx
!> @brief Assembled entry inserted into the sparse matrix.
   real(kind=8) :: val
!> @brief Elementary mass, stiffness, and advection-like terms.
   real(kind=8) :: M, K, B, BT

   mm1 = n1 + p1 + 1
   ng = p1 + 1
   dd = 1
   mm2 = n2 + p2 + 1

! test
   call BasisData(p1, mm1, U1, dd, ng, nelem, O1, J, W, X, NN1)
! trial
   call BasisData(p2, mm2, U2, dd, ng, nelem, O2, J, W, X, NN2)

   call initialize_sparse(n1 + n2 + 2, n1 + n2 + 2, sprsmtrx)

   ! total_size = (nelem1)*(ng1)*(p1 + 1)*(p1 + 1)
! submatrix A: trial-trial block matching the leading rows of Fs
! new parallel loop
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(PRIVATE) &
! !$OMP PRIVATE(d,c,i,e,ia,ib,tmp) &
! !$OMP SHARED(nelem,ng,p,O,NN,W,J,total_size) &
! !$OMP REDUCTION(+:M) &
! !$OMP REDUCTION(+:K) &
! !$OMP REDUCTION(+:B) &
! !$OMP REDUCTION(+:BT)
! loop over elements
   do e = 1, nelem
! loop over Gauss points
      do i = 1, ng
! loop over shape functions over elements (p+1 functions)
         do c = 0, p2
            ! loop over shape functions over elements (p+1 functions)
            do d = 0, p2
               ! O(e) + c = first dof of element + 1st local shape function index
               ! O(e) + d = first dof of element + 2nd local shape function index
               ! NN(0,c,i,e) = value of shape function c at Gauss point i over element e
               ! NN(0,d,i,e) = value of shape function d at Gauss point i over element e
               ! W(i) weight for Gauss point i
               ! J(e) jacobian for element e
               ia = O2(e) + c
               ib = O2(e) + d
               ! M = u*v
               M = NN2(0, c, i, e)*NN2(0, d, i, e)*J(e)*W(i)
               K = NN2(1, c, i, e)*NN2(1, d, i, e)*J(e)*W(i)
               B = NN2(1, c, i, e)*NN2(0, d, i, e)*J(e)*W(i)
               BT = NN2(0, c, i, e)*NN2(1, d, i, e)*J(e)*W(i)
               val = mixA(1)*M + mixA(2)*K + mixA(3)*B + mixA(4)*BT
               call add(sprsmtrx, ia, ib, val)
            end do
         end do
      end do
   end do
! !$OMP END PARALLEL DO

   ! total_size = (nelem1)*(ng1)*(p1 + 1)*(p2 + 1)
! submatrix B
! new parallel loop
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(PRIVATE) &
! !$OMP PRIVATE(d,c,i,e,ia,ib,tmp) &
! !$OMP SHARED(nelem,ng,p,O,NN,W,J,total_size) &
! !$OMP REDUCTION(+:M) &
! !$OMP REDUCTION(+:K) &
! !$OMP REDUCTION(+:B) &
! !$OMP REDUCTION(+:BT)
! loop over elements
   do e = 1, nelem
! loop over Gauss points
      do i = 1, ng
! loop over shape functions over elements (p+1 functions)
         do c = 0, p2
            ! loop over shape functions over elements (p+1 functions)
            do d = 0, p1
               ! O(e) + c = first dof of element + 1st local shape function index
               ! O(e) + d = first dof of element + 2nd local shape function index
               ! NN(0,c,i,e) = value of shape function c at Gauss point i over element e
               ! NN(0,d,i,e) = value of shape function d at Gauss point i over element e
               ! W(i) weight for Gauss point i
               ! J(e) jacobian for element e
               ia = O2(e) + c
               ib = O1(e) + d + n2 + 1
               ! M = u*v
               M = NN2(0, c, i, e)*NN1(0, d, i, e)*J(e)*W(i)
               K = NN2(1, c, i, e)*NN1(1, d, i, e)*J(e)*W(i)
               B = NN2(1, c, i, e)*NN1(0, d, i, e)*J(e)*W(i)
               BT = NN2(0, c, i, e)*NN1(1, d, i, e)*J(e)*W(i)
               val = mixB(1)*M + mixB(2)*K + mixB(3)*B + mixB(4)*BT
               call add(sprsmtrx, ia, ib, val)
            end do
         end do
      end do
   end do
! !$OMP END PARALLEL DO

   ! total_size = (nelem2)*(ng2)*(p2 + 1)*(p1 + 1)
! submatrix BT
! new parallel loop
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(PRIVATE) &
! !$OMP PRIVATE(d,c,i,e,ia,ib,tmp) &
! !$OMP SHARED(nelem,ng,p,O,NN,W,J,total_size) &
! !$OMP REDUCTION(+:M) &
! !$OMP REDUCTION(+:K) &
! !$OMP REDUCTION(+:B) &
! !$OMP REDUCTION(+:BT)
   do e = 1, nelem
! loop over Gauss points
      do i = 1, ng
! loop over shape functions over elements (p+1 functions)
         do c = 0, p1
! loop over shape functions over elements (p+1 functions)
            do d = 0, p2
               ! O(e) + c = first dof of element + 1st local shape function index
               ! O(e) + d = first dof of element + 2nd local shape function index
               ! NN(0,c,i,e) = value of shape function c at Gauss(i) weight for Gauss point i
               ! J(e) jacobian for element e
               ia = O1(e) + c + n2 + 1
               ib = O2(e) + d
               ! M = u*v
               M = NN1(0, c, i, e)*NN2(0, d, i, e)*J(e)*W(i)
               K = NN1(1, c, i, e)*NN2(1, d, i, e)*J(e)*W(i)
               B = NN1(1, c, i, e)*NN2(0, d, i, e)*J(e)*W(i)
               BT = NN1(0, c, i, e)*NN2(1, d, i, e)*J(e)*W(i)
               val = mixBT(1)*M + mixBT(2)*K + mixBT(3)*B + mixBT(4)*BT
               call add(sprsmtrx, ia, ib, val)
            end do
         end do
      end do
   end do
! !$OMP END PARALLEL DO

end subroutine MKBBT_large

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Assembles a sparse matrix on a single spline space from a
!> weighted combination of standard bilinear-form terms.
!>
!> @details
!> This routine constructs a sparse operator matrix on one
!> one-dimensional spline space. For each element and Gauss point it
!> combines:
!> - the mass term \f$M = u\,v\f$,
!> - the stiffness term \f$K = u'\,v'\f$,
!> - the advection-like term \f$B = u'\,v\f$,
!> - the transposed advection-like term \f$B^T = u\,v'\f$.
!>
!> The final entry added to the sparse matrix is obtained from the
!> mixing vector \p mix.
!
! Input:
! ------
!> @param[in] nelem
!> Number of elements used during assembly.
!>
!> @param[in] U
!> Knot vector of the spline space.
!>
!> @param[in] p
!> Polynomial degree of the spline basis.
!>
!> @param[in] n
!> Number of control points minus one.
!>
!> @param[in] mix
!> Mixing coefficients for the elementary terms \f$M\f$, \f$K\f$,
!> \f$B\f$, and \f$B^T\f$.
!
! Output:
! -------
!> @param[out] sprsmtrx
!> Sparse matrix containing the assembled one-space operator.
!
! Notes:
! ------
!> @note
!> The output matrix is initialized with logical dimensions
!> \f$(n+1) \times (n+1)\f$.
!
!---------------------------------------------------------------------------
subroutine MKBBT_small(nelem, U, p, n, mix, sprsmtrx)
   use basis, ONLY: BasisData
   use omp_lib
   use sparse
   implicit none
!> @brief Basis size minus one, polynomial degree, and number of elements.
   integer(kind=4), intent(in) :: n, p, nelem
!> @brief Knot vector of the spline space.
   real(kind=8), intent(in) :: U(0:n + p + 1)
!> @brief Mixing coefficients for the elementary bilinear forms.
   real(kind=8), dimension(4), intent(in) :: mix
!> @brief Element Jacobians.
   real(kind=8), dimension(nelem) :: J
!> @brief Gauss quadrature weights.
   real(kind=8), dimension(p + 1) :: W
!> @brief Physical Gauss-point coordinates.
   real(kind=8), dimension(p + 1, nelem) :: X
!> @brief Basis values and first derivatives at quadrature points.
   real(kind=8), dimension(0:1, 0:p, p + 1, nelem) :: NN
!> @brief Highest derivative order requested from basis evaluation.
   integer(kind=4) :: dd
!> @brief Row and column indices of the assembled sparse entry.
   integer(kind=4) :: ia, ib
!> @brief Last valid knot index.
   integer(kind=4) :: mm
!> @brief Number of Gauss points.
   integer(kind=4) :: ng
!> @brief Loop counters.
   integer(kind=4) :: e, i, c, d
!> @brief First nonzero basis-function index on each element.
   integer(kind=4) :: O(nelem)
!> @brief Output sparse matrix.
   type(sparse_matrix), pointer, intent(out) :: sprsmtrx
!> @brief Assembled matrix entry.
   real(kind=8) :: val
!> @brief Elementary mass, stiffness, and advection-like terms.
   real(kind=8) :: M, K, B, BT

   mm = n + p + 1
   ng = p + 1
   dd = 1

   call BasisData(p, mm, U, dd, ng, nelem, O, J, W, X, NN)

   call initialize_sparse(n + 1, n + 1, sprsmtrx)

! submatrix A
! new parallel loop
!!$OMP PARALLEL DO &
!!$OMP DEFAULT(PRIVATE) &
!!$OMP PRIVATE(d,c,i,e,ia,ib,tmp) &
!!$OMP SHARED(nelem,ng,p,O,NN,W,J,total_size) &
!!$OMP REDUCTION(+:M) &
!!$OMP REDUCTION(+:K) &
!!$OMP REDUCTION(+:B) &
!!$OMP REDUCTION(+:BT)
! loop over elements
   do e = 1, nelem
! loop over Gauss points
      do i = 1, ng
         ! loop over shape functions over elements (p+1 functions)
         do c = 0, p
            ! loop over shape functions over elements (p+1 functions)
            do d = 0, p
               ! O(e) + c = first dof of element + 1st local shape function index
               ! O(e) + d = first dof of element + 2nd local shape function index
               ! NN(0,c,i,e) = value of shape function c at Gauss point i over element e
               ! NN(0,d,i,e) = value of shape function d at Gauss point i over element e
               ! W(i) weight for Gauss point i
               ! J(e) jacobian for element e
               ia = O(e) + c
               ib = O(e) + d
               ! M = u*v
               M = NN(0, c, i, e)*NN(0, d, i, e)*J(e)*W(i)
               K = NN(1, c, i, e)*NN(1, d, i, e)*J(e)*W(i)
               B = NN(1, c, i, e)*NN(0, d, i, e)*J(e)*W(i)
               BT = NN(0, c, i, e)*NN(1, d, i, e)*J(e)*W(i)
               val = mix(1)*M + mix(2)*K + mix(3)*B + mix(4)*BT
               call add(sprsmtrx, ia, ib, val)
            end do
         end do
      end do
   end do
!!$OMP END PARALLEL DO

end subroutine MKBBT_small

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
!> @brief Element Jacobian product and quadrature weight product.
   real(kind=8) :: J, W
!> @brief Global and local index variables.
   integer(kind=4) :: ind, ind1, ind23, indx, indy, indz
!> @brief Pointwise contribution returned by `ComputePointForRHS`.
   real(kind=8) :: resvalue
!> @brief Physical coordinates of the current quadrature point.
   real(kind=8), dimension(3) :: X
!> @brief Direction-dependent quadrature, element, and local basis indices.
   integer(kind=4), dimension(3) :: k, e, a, indb
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
!> @brief Mixed-space setup used for the current assembly pass.
   type(ADS_setup) :: ads
!> @brief Flag indicating use of the iGRM path.
   logical, intent(out) :: igrm
!> @brief Permutation of coordinate directions used internally.
   integer(kind=4) :: dira,dirb,dirc

   call create_mixed_space(ads_test, ads_trial, direction,&
   ads, dira, dirb, dirc, igrm)

   allocate (elarr(0:ads%p(dira), 0:ads%p(dirb), 0:ads%p(dirc)))
   ! total_size = ads % lnelem(1) * ads % lnelem(2) * ads % lnelem(3)

!   if (allocated(ads_data%F)) ads_data%F = 0.d0
!   if (allocated(ads_data%Ft)) ads_data%Ft = 0.d0

!      loop over points
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(SHARED) &
! !$OMP SHARED(ads,ads_data,total_size) &
! !$OMP PRIVATE(tmp,ex,ey,ez,e,kx,ky,kz,k,W,ax,ay,az,a,ind,indx,indy,indz,ind1,ind23,J) &
! !$OMP PRIVATE(X,du,resvalue) &
! !$OMP PRIVATE(indbx,indby,indbz,Uval,elarr,,Uval_m,Uval13,Uval23)
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
                     Uval = ads_data%Un(ex, ey, ez, k(1), k(2), k(3))
                     Uval13 = ads_data%Un13(ex, ey, ez, k(1), k(2), k(3))
                     Uval23 = ads_data%Un23(ex, ey, ez, k(1), k(2), k(3))
                     du = ads_data%dUn(ex, ey, ez, k(1), k(2), k(3), :)

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
! moving results from temporary array to main one
! !$OMP CRITICAL
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

                        if (igrm) then
                           ads_data%Ft(ind1 + 1, ind23 + 1) = &
                              ads_data%Ft(ind1 + 1, ind23 + 1) &
                              + elarr(ax, ay, az)
                        else
                           ads_data%F(ind1 + 1, ind23 + 1) = &
                              ads_data%F(ind1 + 1, ind23 + 1) &
                              + elarr(ax, ay, az)
                        end if
                     end if
                  end do
               end do
            end do
! !$OMP END CRITICAL
         end do
      end do
   end do
! !$OMP END PARALLEL DO

   deallocate (elarr)

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
!> - \p ads_data%Un,
!> - \p ads_data%Un13,
!> - \p ads_data%Un23,
!>
!> according to the selector \p subun. The derivative vector is written
!> to \p ads_data%dUn.
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
!
!---------------------------------------------------------------------------
subroutine FormUn(subun, ads, ads_data)
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   ! use parallelism, ONLY: PRINTRANK
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

   select case (subun)
   case (1)
      ads_data%Un = 0.d0
   case (2)
      ads_data%Un13 = 0.d0
   case (3)
      ads_data%Un23 = 0.d0
   case default
      write (ERROR_UNIT, *) "wrong substep"
      stop 1
   end select
   ads_data%dUn = 0.d0
   !total_size = ads%lnelem(1)*ads%lnelem(2)*ads%lnelem(3)

!      loop over points
! !$OMP PARALLEL DO &
! !$OMP DEFAULT(SHARED) &
! !$OMP SHARED(ads,ads_data,total_size) &
! !$OMP PRIVATE(tmp,ex,ey,ez,kx,ky,kz,ind) &
! !$OMP PRIVATE(bx,by,bz,rx,ry,rz,ix,iy,iz,sx,sy,sz,Ucoeff,dvx,dvy,dvz,du) &
! !$OMP PRIVATE(indbx,indby,indbz,Uval,dux,duy,duz,v)
   !do all = 1, total_size
!        translate coefficients to local
   !ez = modulo(all - 1, ads%lnelem(3))
   !tmp = (all - ez)/ads%lnelem(3) + 1
   !ey = modulo(tmp - 1, ads%lnelem(2))
   !ex = (tmp - ey)/ads%lnelem(2)
   do exx=1,ads%lnelem(1)
      do eyy=1,ads%lnelem(2)
         do ezz=1,ads%lnelem(3)
!        fix distributed part
            ex = exx + ads%mine(1)-1
            ey = eyy + ads%mine(2)-1
            ez = ezz + ads%mine(3)-1
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
                     ads_data%dUn(ex, ey, ez, kx, ky, kz, :) = (/dux, duy, duz/)
                     if (subun .EQ. 1) then
                        ads_data%Un(ex, ey, ez, kx, ky, kz) = Uval
                     else if (subun .EQ. 2) then
                        ads_data%Un13(ex, ey, ez, kx, ky, kz) = Uval
                     else if (subun .EQ. 3) then
                        ads_data%Un23(ex, ey, ez, kx, ky, kz) = Uval
                     else
                        write (ERROR_UNIT, *) "wrong substep"
                     end if
                  end do
               end do
            end do
         end do
      end do
   end do
! !$OMP END PARALLEL DO

end subroutine FormUn

!!!!! to nie tu
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

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Selects and assembles the appropriate sparse operator matrix.
!>
!> @details
!> This routine dispatches matrix assembly to one of two lower-level
!> procedures:
!> - \ref MKBBT_small for the equal-space case,
!> - \ref MKBBT_large for the mixed-space case.
!>
!> The choice is controlled by the logical flag \p equ. In the equal-space
!> case, only one spline space and the mixing vector \p mixA are used.
!> Otherwise, a mixed block matrix is assembled from the two spline
!> spaces and the three mixing vectors \p mixA, \p mixB, and \p mixBT.
!
! Input:
! ------
!> @param[in] U1
!> Knot vector of the first spline space.
!>
!> @param[in] p1
!> Polynomial degree of the first spline space.
!>
!> @param[in] n1
!> Number of control points minus one for the first spline space.
!>
!> @param[in] nelem1
!> Number of elements in the first spline space.
!>
!> @param[in] U2
!> Knot vector of the second spline space.
!>
!> @param[in] p2
!> Polynomial degree of the second spline space.
!>
!> @param[in] n2
!> Number of control points minus one for the second spline space.
!>
!> @param[in] nelem2
!> Number of elements in the second spline space.
!>
!> @param[in] mixA
!> Mixing coefficients for the diagonal block or the equal-space case.
!>
!> @param[in] mixB
!> Mixing coefficients for one mixed off-diagonal block.
!>
!> @param[in] mixBT
!> Mixing coefficients for the transposed mixed off-diagonal block.
!>
!> @param[in] equ
!> Logical flag indicating whether both spaces are treated as equal.
!
! Output:
! -------
!> @param[out] sprsmtrx
!> Sparse matrix returned by the selected assembly routine.
!
!---------------------------------------------------------------------------
subroutine ComputeMatrix(U1, p1, n1, nelem1, U2, p2, n2, nelem2, mixA, mixB, mixBT, equ, sprsmtrx)
   ! use parallelism, ONLY: PRINTRANK
   use sparse
   implicit none
!> @brief Basis size minus one, polynomial degree, and element count for the first space.
   integer(kind=4), intent(in) :: n1, p1, nelem1
!> @brief Knot vector of the first space.
   real(kind=8), dimension(0:n1 + p1 + 1), intent(in) :: U1
!> @brief Knot vector of the second space.
   integer(kind=4), intent(in) :: n2, p2, nelem2
!> @brief Knot vector of the second space.
   real(kind=8), dimension(0:n2 + p2 + 1), intent(in) :: U2
!> @brief Mixing coefficients used by the selected assembly routine.
   real(kind=8), dimension(4), intent(in) :: mixA, mixB, mixBT
!> @brief Equality flag selecting equal-space or mixed-space assembly.
   logical, intent(in) :: equ
!> @brief Output sparse matrix.
   type(sparse_matrix), pointer, intent(out) :: sprsmtrx
   ! integer :: i

   if (equ) then
      call MKBBT_small(nelem2, U2, p2, n2, mixA, sprsmtrx)
   else
      call ValidateIGRMMesh(U1, p1, n1, nelem1, U2, p2, n2, nelem2)
      call MKBBT_large(nelem2, U1, p1, n1, U2, p2, n2, mixA, mixB, mixBT, sprsmtrx)
   end if

end subroutine ComputeMatrix

end module projection_engine
