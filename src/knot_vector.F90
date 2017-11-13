module knot_vector

contains


   ! -------------------------------------------------------------------
   ! Allocates and fills the knot vector on [0, 1]
   !
   ! Input:
   ! ------
   ! U  - array to fill with points
   ! n  - number of functions on the knot minus one
   ! p  - degree of polynomial
   !
   ! Output:
   ! ------
   ! nelem - number of elements
   !
   ! Number of subintervals is N = n-p+1.
   ! 0 and 1 are repeated (p+1) times.
   ! -------------------------------------------------------------------
   subroutine PrepareKnot(U, n, p, nelem)
      implicit none
      integer(kind = 4), intent(in) :: n, p
      real (kind = 8), allocatable, dimension(:), intent(out) :: U
      integer(kind = 4), intent(out) :: nelem


      allocate(U(n + p + 2))
      call FillOpenKnot(U, n, p)
      nelem = CountSpans(n, p, U)

#ifdef IINFO
      write(*, *) 'n,p,nelem', n, p, nelem
      write(*, *) 'U', U
#endif
      
end subroutine



      ! -------------------------------------------------------------------
      ! Fills knot vector on [0, 1]
      !
      ! Input:
      ! ------
      ! n  - number of functions on the knot minus one
      ! p  - degree of polynomial
      !
      ! Output:
      ! ------
      ! U  - array to fill with points
      !
      ! Number of subintervals is N = n-p+1.
      ! 0 and 1 are repeated (p+1) times.
      ! -------------------------------------------------------------------
      subroutine FillOpenKnot(U, n, p)
         implicit none
         integer(kind = 4), intent(in) :: n, p
         real (kind = 8), intent(out) :: U(1:n + p + 2)
         integer(kind = 4) :: i

         U(1: p + 1) = 0.d0
         U(n + 2: n + p + 2) = 1.d0

         do i = p + 2, n + 1
            U(i) = real(i - p - 1) / real(n - p + 1)
         enddo

      end subroutine



      ! -------------------------------------------------------------------
      ! Calculates number of elements (nonempty subintervals) in the knot 
      ! vector.
      !
      ! Input:
      ! ------
      ! n     - number of functions (control points) minus 1
      ! p     - order of basis functions
      ! U     - knot vector
      !
      ! Output:
      ! ------
      ! nelem - number of elements
      !
      ! -------------------------------------------------------------------
      function CountSpans(n, p, U) result (nelem)
         implicit none
         integer(kind = 4), intent(in) :: n, p
         real (kind = 8), intent(in) :: U(0:n + p + 1)
         integer(kind = 4) :: i, nelem

         nelem = 0
         i = p
         do while (i <= n)
            ! skip multiple knots
            do while (i < n .and. U(i) == U(i + 1))
               i = i + 1
            enddo
#ifdef IPRINT
            write(*, *) 'CountSpans:i,n,U(i),U(i+1)', i, n, U(i), U(i + 1)
#endif
            nelem = nelem + 1
            i = i + 1
         enddo

      end function




   end module
