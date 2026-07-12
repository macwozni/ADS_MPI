!------------------------------------------------------------------------------
!
! MODULE: operator_assembly
!
! DESCRIPTION:
!> @file operator_assembly.F90
!> @brief One-dimensional ADS/iGRM operator matrix assembly routines.
!>
!> @details
!> This module owns the one-dimensional operator builders used by the
!> directional ADS solve path. It assembles mass, stiffness, and
!> first-derivative coupling contributions into the sparse format later
!> converted for MUMPS.
!>
!> The routines support both equal-space ADS operators and mixed
!> test/trial iGRM block operators. They are kept together because they
!> share basis evaluation, element coloring, sparse-pattern construction,
!> and deterministic sparse insertion conventions.
!>
!> The assembly routines do not solve the resulting systems. They only
!> build the local sparse operator selected by \ref ComputeMatrix.
!
!------------------------------------------------------------------------------
module operator_assembly

   implicit none

   private
   public :: MKBBT_large
   public :: MKBBT_small
   public :: ComputeMatrix

contains

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
!> @brief Element color counts used to avoid concurrent writes to the same rows.
   integer(kind=4) :: color, ncolors1, ncolors2
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
!> @brief Element-local matrices accumulated before sparse scatter.
   real(kind=8), dimension(0:p2, 0:p2) :: elmatA
   real(kind=8), dimension(0:p2, 0:p1) :: elmatB
   real(kind=8), dimension(0:p1, 0:p2) :: elmatBT

   mm1 = n1 + p1 + 1
   ng = p1 + 1
   dd = 1
   mm2 = n2 + p2 + 1

! test
   call BasisData(p1, mm1, U1, dd, ng, nelem, O1, J, W, X, NN1)
! trial
   call BasisData(p2, mm2, U2, dd, ng, nelem, O2, J, W, X, NN2)

   call initialize_sparse(n1 + n2 + 2, n1 + n2 + 2, sprsmtrx)

! Build the sparse pattern once for all three mixed blocks. The numerical
! phase can then update existing entries without reallocating row storage.
   do e = 1, nelem
      do c = 0, p2
         ia = O2(e) + c
         do d = 0, p2
            ib = O2(e) + d
            call add(sprsmtrx, ia, ib, 0.d0)
         end do
         do d = 0, p1
            ib = O1(e) + d + n2 + 1
            call add(sprsmtrx, ia, ib, 0.d0)
         end do
      end do
      do c = 0, p1
         ia = O1(e) + c + n2 + 1
         do d = 0, p2
            ib = O2(e) + d
            call add(sprsmtrx, ia, ib, 0.d0)
         end do
      end do
   end do

   ncolors2 = p2 + 1
   do color = 1, ncolors2
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(e,i,c,d,ia,ib,M,K,B,BT,val,elmatA,elmatB) SCHEDULE(STATIC)
      do e = color, nelem, ncolors2
         elmatA = 0.d0
         elmatB = 0.d0
! loop over Gauss points
         do i = 1, ng
! submatrix A: trial-trial block matching the leading rows of Fs
            do c = 0, p2
               do d = 0, p2
                  M = NN2(0, c, i, e)*NN2(0, d, i, e)*J(e)*W(i)
                  K = NN2(1, c, i, e)*NN2(1, d, i, e)*J(e)*W(i)
                  B = NN2(1, c, i, e)*NN2(0, d, i, e)*J(e)*W(i)
                  BT = NN2(0, c, i, e)*NN2(1, d, i, e)*J(e)*W(i)
                  val = mixA(1)*M + mixA(2)*K + mixA(3)*B + mixA(4)*BT
                  elmatA(c, d) = elmatA(c, d) + val
               end do
            end do

! submatrix B
            do c = 0, p2
               do d = 0, p1
                  M = NN2(0, c, i, e)*NN1(0, d, i, e)*J(e)*W(i)
                  K = NN2(1, c, i, e)*NN1(1, d, i, e)*J(e)*W(i)
                  B = NN2(1, c, i, e)*NN1(0, d, i, e)*J(e)*W(i)
                  BT = NN2(0, c, i, e)*NN1(1, d, i, e)*J(e)*W(i)
                  val = mixB(1)*M + mixB(2)*K + mixB(3)*B + mixB(4)*BT
                  elmatB(c, d) = elmatB(c, d) + val
               end do
            end do
         end do

         do c = 0, p2
            ia = O2(e) + c
            do d = 0, p2
               ib = O2(e) + d
               call add_existing(sprsmtrx, ia, ib, elmatA(c, d))
            end do
            do d = 0, p1
               ib = O1(e) + d + n2 + 1
               call add_existing(sprsmtrx, ia, ib, elmatB(c, d))
            end do
         end do
      end do
!$OMP END PARALLEL DO
   end do

   ncolors1 = p1 + 1
   do color = 1, ncolors1
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(e,i,c,d,ia,ib,M,K,B,BT,val,elmatBT) SCHEDULE(STATIC)
      do e = color, nelem, ncolors1
         elmatBT = 0.d0
! loop over Gauss points
         do i = 1, ng
! submatrix BT
            do c = 0, p1
               do d = 0, p2
                  M = NN1(0, c, i, e)*NN2(0, d, i, e)*J(e)*W(i)
                  K = NN1(1, c, i, e)*NN2(1, d, i, e)*J(e)*W(i)
                  B = NN1(1, c, i, e)*NN2(0, d, i, e)*J(e)*W(i)
                  BT = NN1(0, c, i, e)*NN2(1, d, i, e)*J(e)*W(i)
                  val = mixBT(1)*M + mixBT(2)*K + mixBT(3)*B + mixBT(4)*BT
                  elmatBT(c, d) = elmatBT(c, d) + val
               end do
            end do
         end do

         do c = 0, p1
            ia = O1(e) + c + n2 + 1
            do d = 0, p2
               ib = O2(e) + d
               call add_existing(sprsmtrx, ia, ib, elmatBT(c, d))
            end do
         end do
      end do
!$OMP END PARALLEL DO
   end do

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
!> @brief Element color used to avoid concurrent writes to the same rows.
   integer(kind=4) :: color, ncolors
!> @brief First nonzero basis-function index on each element.
   integer(kind=4) :: O(nelem)
!> @brief Output sparse matrix.
   type(sparse_matrix), pointer, intent(out) :: sprsmtrx
!> @brief Assembled matrix entry.
   real(kind=8) :: val
!> @brief Elementary mass, stiffness, and advection-like terms.
   real(kind=8) :: M, K, B, BT
!> @brief Element-local matrix accumulated before sparse scatter.
   real(kind=8), dimension(0:p, 0:p) :: elmat

   mm = n + p + 1
   ng = p + 1
   dd = 1

   call BasisData(p, mm, U, dd, ng, nelem, O, J, W, X, NN)

   call initialize_sparse(n + 1, n + 1, sprsmtrx)

! Build the sparse pattern once. The numerical phase can then update
! existing entries without reallocating row storage.
   do e = 1, nelem
      do c = 0, p
         ia = O(e) + c
         do d = 0, p
            ib = O(e) + d
            call add(sprsmtrx, ia, ib, 0.d0)
         end do
      end do
   end do

   ncolors = p + 1
   do color = 1, ncolors
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(e,i,c,d,ia,ib,M,K,B,BT,val,elmat) SCHEDULE(STATIC)
      do e = color, nelem, ncolors
         elmat = 0.d0
! loop over Gauss points
         do i = 1, ng
            ! loop over shape functions over elements (p+1 functions)
            do c = 0, p
               ! loop over shape functions over elements (p+1 functions)
               do d = 0, p
                  ! M = u*v
                  M = NN(0, c, i, e)*NN(0, d, i, e)*J(e)*W(i)
                  K = NN(1, c, i, e)*NN(1, d, i, e)*J(e)*W(i)
                  B = NN(1, c, i, e)*NN(0, d, i, e)*J(e)*W(i)
                  BT = NN(0, c, i, e)*NN(1, d, i, e)*J(e)*W(i)
                  val = mix(1)*M + mix(2)*K + mix(3)*B + mix(4)*BT
                  elmat(c, d) = elmat(c, d) + val
               end do
            end do
         end do

         do c = 0, p
            ia = O(e) + c
            do d = 0, p
               ib = O(e) + d
               call add_existing(sprsmtrx, ia, ib, elmat(c, d))
            end do
         end do
      end do
!$OMP END PARALLEL DO
   end do

end subroutine MKBBT_small

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
      call MKBBT_large(nelem2, U1, p1, n1, U2, p2, n2, mixA, mixB, mixBT, sprsmtrx)
   end if

end subroutine ComputeMatrix

end module operator_assembly
