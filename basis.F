      module basis

      use gauss
      
      contains
       
c p - polynomial order
c m - 
c U - knot vector
c d - derivatives
c q=p+1 - order for Gauss integration points
c r - number of elements
c O,J,W,X,N - output
      subroutine BasisData(p,m,U,d,q,r,O,J,W,X,N)
      implicit none
      integer(kind=4), intent(in)  :: p, m
      real   (kind=8), intent(in)  :: U(0:m)
      integer(kind=4), intent(in)  :: d, q, r
      integer(kind=4), intent(out) :: O(r)
      real   (kind=8), intent(out) :: J(r)
      real   (kind=8), intent(out) :: W(q)
      real   (kind=8), intent(out) :: X(q,r)
      real   (kind=8), intent(out) :: N(0:d,0:p,q,r)  

      integer(kind=4) i, iq, ir
      real   (kind=8) uu, Xg(q)
      real   (kind=8) basis(0:p,0:d)

c      write(*,*)'BasisData:p,r,U',p,r,U
      ir = 1
      do i = p, m-p
         if (U(i) /= U(i+1)) then
            O(ir) = i - p
            ir = ir + 1
         end if
      end do
c      write(*,*)'O',O
      call GaussRule(q,Xg,W)
c      write(*,*)'W',W
      do ir = 1, r
         i = O(ir) + p
         J(ir) = (U(i+1)-U(i))/2.0
         X(:,ir) = (Xg + 1.0) * J(ir) + U(i)
         do iq = 1, q
            uu = X(iq,ir)
            call DersBasisFuns(i,uu,p,d,U,basis)
            N(:,:,iq,ir) = transpose(basis)
         end do
      end do
c      write(*,*)'J',J
c      write(*,*)'X',X
c      write(*,*)'N',N
c      call pause
      end subroutine BasisData
 
c i - function index
c uu - coordinates
c p - order
c d - derivatives
c U - knot vector
c ders - output rhs
      subroutine DersBasisFuns(i,uu,p,d,U,ders)
      implicit none
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
         end do
         ndu(j,j) = saved
      end do
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
            end if
            if (rk > -1) then
               j1 = 1
            else
               j1 = -rk
            end if
            if (r-1 <= pk) then
               j2 = k-1
            else
               j2 = p-r
            end if
            do j = j1, j2
               a(s2,j) = (a(s1,j) - a(s1,j-1)) / ndu(pk+1,rk+j)
               der =  der + a(s2,j) * ndu(rk+j,pk)
            end do
            if (r <= pk) then
               a(s2,k) = - a(s1,k-1) / ndu(pk+1,r)
               der =  der + a(s2,k) * ndu(r,pk)
            end if
            ders(r,k) = der
            j = s1; s1 = s2; s2 = j;
         end do
      end do
      r = p
      do k = 1, d
         ders(:,k) = ders(:,k) * r
         r = r * (p-k)
      end do
      end subroutine DersBasisFuns   

      
c n - 
c p - order
c uu -coordinates
c U - knot vector
      function FindSpan(n,p,uu,U) result (span)
      implicit none
      integer(kind=4), intent(in) :: n, p
      real   (kind=8), intent(in) :: uu, U(0:n+p+1)
      integer(kind=4)             :: span
      integer(kind=4) low, high
      if (uu >= U(n+1)) then
         span = n
         return
      end if
      if (uu <= U(p)) then
         span = p
         return
      end if
      low  = p
      high = n+1
      span = (low + high) / 2
      do while (uu < U(span) .or. uu >= U(span+1))
         if (uu < U(span)) then
            high = span
         else
            low  = span
         end if
         span = (low + high) / 2
      end do
      end function FindSpan
 
      function CountSpans(n,p,U) result (nelem)
      implicit none
      integer(kind=4), intent(in) :: n, p
      real   (kind=8), intent(in) :: U(0:n+p+1)
      integer(kind=4) :: i,nelem  
 
      nelem = 0
      i=p
      do while (i.le.n)
         do while ((i.le.n) .and. (U(i) == U(i+1)))
            i = i + 1
         end do
         write(*,*)'CountSpans:i,n,U(i),U(i+1)',i,n,U(i),U(i+1)
         nelem = nelem + 1
         i = i + 1
      end do
      end function CountSpans 
      
      end module basis
