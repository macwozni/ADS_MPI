module basis

use gauss

implicit none

contains
     

! -------------------------------------------------------------------
! Computes all the information about basis B-spline functions needed
! to perform numerical integration.
!
! Input:
! ------
! p       - polynomial order
! m       - index of the last node in knot vector (number of nodes - 1)
! U       - knot vector
! d       - order of highest derivatives we want to compute
! q = p+1 - number of Gauss quadrature points
! r       - number of elements (subintervals)
! 
! Output:
! -------
! O       - indexes of first nonzero functions on each element
! J       - values of the Jacobian of elements
! W       - weights of Gauss quadrature points
! X       - points of Gauss quadrature
! N       - values of (p+1) nonzero basis functions and their derivatives
!           at points of Gauss quadrature
!
! Array of values - N - is indexed by:
!
!  - derivative
!  - numer of nonzero function
!  - index of Gauss quadrature point
!  - index of element
!
! Second index is local - for each element, there are (p+1) nonzero
! basis functions, indexed 0...p. It does not refer to global index
! of basis functions on the whole interval. 
!
! To translate between local and global indexing, use O array. O(i) is 
! the global index of first nonzero function on i-th element, hence
!
!     global_index = O(i) + local_index
!
! where local_index is local for i-th element.
!
! Value of q (number of quadrature points) is taken to be p+1, since
! order of Gaussian quadrature is 2n - 1, hence using q = p+1 we get
! exact results for polynomials of degree up to 2q - 1 = 2p + 1.
! This is nice, as we integrate mostly bilinear forms in B-spline
! basis, consisting of dogree p polynomials.
! -------------------------------------------------------------------
subroutine BasisData(p, m, U, d, q, r, O, J, W, X, N)
integer (kind=4), intent(in)  :: p, m
real    (kind=8), intent(in)  :: U(0:m)
integer (kind=4), intent(in)  :: d, q, r
integer (kind=4), intent(out) :: O(r)
real    (kind=8), intent(out) :: J(r)
real    (kind=8), intent(out) :: W(q)
real    (kind=8), intent(out) :: X(q,r)
real    (kind=8), intent(out) :: N(0:d,0:p,q,r)

integer (kind=4) i, iq, ir
real    (kind=8) uu, Xg(q)
real    (kind=8) basis(0:p,0:d)

  ! Calculates first nonzero basis function for each element
  ir = 1
  do i = p, m-p
     if (U(i) /= U(i+1)) then
        O(ir) = i - p
        ir = ir + 1
     endif
  enddo

  call GaussRule(q, Xg, W)

  do ir = 1, r
     i = O(ir) + p
     J(ir) = (U(i+1) - U(i)) / 2.0
     X(:,ir) = (Xg + 1.0) * J(ir) + U(i) ! translate Gauss [-1,1] -> [0,1]
     do iq = 1, q
        uu = X(iq,ir)
        call DersBasisFuns(i,uu,p,d,U,basis)
        N(:,:,iq,ir) = transpose(basis)
     enddo
  enddo

end subroutine


! -------------------------------------------------------------------
! Computes values of basis functions and their derivatives at the
! specified point.
!
! Input:
! ------
! i      - index of (last?) basis function to compute
! uu     - coordinate of the point (argument)
! p      - order of B-splines
! d      - order of highest derivatives to compute
! U      - knot vector
!
! Output:
! ------
! ders   - computed values of basis functions and derivatives
!
! Output array is indexed by local index of nonzero basis functions
! on the element, and order of derivative.
! -------------------------------------------------------------------
subroutine DersBasisFuns(i, uu, p, d, U, ders)
integer(kind=4), intent(in) :: i, p, d
real   (kind=8), intent(in) :: uu, U(0:i+p)
real   (kind=8), intent(out):: ders(0:p,0:d)
integer(kind=4) :: j, k, r, s1, s2, rk, pk, j1, j2
real   (kind=8) :: saved, temp, der
real   (kind=8) :: left(p), right(p)
real   (kind=8) :: ndu(0:p,0:p), a(0:1,0:p)

  ndu(0,0) = 1.0
  do j = 1, p
     left(j)  = uu - U(i+1-j)
     right(j) = U(i+j) - uu
     saved = 0.0
     do r = 0, j-1
        ndu(j,r) = right(r+1) + left(j-r)
        temp = ndu(r,j-1) / ndu(j,r)
        ndu(r,j) = saved + right(r+1) * temp
        saved = left(j-r) * temp
     enddo
     ndu(j,j) = saved
  enddo

  ders(:,0) = ndu(:,p)

  do r = 0, p
     s1 = 0; s2 = 1;
     a(0,0) = 1.0
     do k = 1, d
        der = 0.0
        rk = r-k; pk = p-k;
        if (r >= k) then
           a(s2,0) = a(s1,0) / ndu(pk+1,rk)
           der =  a(s2,0) * ndu(rk,pk)
        endif
        if (rk > -1) then
           j1 = 1
        else
           j1 = -rk
        endif
        if (r-1 <= pk) then
           j2 = k-1
        else
           j2 = p-r
        endif
        do j = j1, j2
           a(s2,j) = (a(s1,j) - a(s1,j-1)) / ndu(pk+1,rk+j)
           der =  der + a(s2,j) * ndu(rk+j,pk)
        enddo
        if (r <= pk) then
           a(s2,k) = - a(s1,k-1) / ndu(pk+1,r)
           der =  der + a(s2,k) * ndu(r,pk)
        endif
        ders(r,k) = der
        j = s1; s1 = s2; s2 = j;
     enddo
  enddo
  r = p
  do k = 1, d
     ders(:,k) = ders(:,k) * r
     r = r * (p-k)
  enddo

end subroutine

      
! -------------------------------------------------------------------
! Finds the element that contains specified point. It is given as
! the index of knot u_i such that [u_i, u_(i+1)) contains the point.
! In case the point lies on the edge, index of innermost edge node
! is returned.
! 
! n    - number of functions (control points) minus 1
! p    - order of basis functions
! uu   - coordinate of the point
! U    - knot vector
! -------------------------------------------------------------------
function FindSpan(n, p, uu, U) result (span)
integer(kind=4), intent(in) :: n, p
real   (kind=8), intent(in) :: uu, U(0:n+p+1)
integer(kind=4)  :: span
integer(kind=4) low, high

  ! check edge cases
  if (uu >= U(n+1)) then
     span = n
     return
  endif

  if (uu <= U(p)) then
     span = p
     return
  endif

  ! Binary search for uu
  low  = p
  high = n+1
  span = (low + high) / 2

  do while (uu < U(span) .or. uu >= U(span+1))
     if (uu < U(span)) then
        high = span
     else
        low  = span
     endif
     span = (low + high) / 2
  enddo

end function
 

! -------------------------------------------------------------------
! Calculates number of elements (nonempty subintervals) in the knot 
! vector.
!
! n   - number of functions (control points) minus 1
! p   - order of basis functions
! U   - knot vector
! -------------------------------------------------------------------
function CountSpans(n, p, U) result (nelem)
integer(kind=4), intent(in) :: n, p
real   (kind=8), intent(in) :: U(0:n+p+1)
integer(kind=4) :: i, nelem  

  nelem = 0
  i = p
  do while (i <= n)
     ! skip multiple knots
     do while (i < n .and. U(i) == U(i+1))
        i = i + 1
     enddo
     write(*,*)'CountSpans:i,n,U(i),U(i+1)',i,n,U(i),U(i+1)
     nelem = nelem + 1
     i = i + 1
  enddo

end function
      
end module

