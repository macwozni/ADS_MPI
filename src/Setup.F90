module Setup

   implicit none

   type ADS_setup
      ! Number of functions in each dimension minus one
      integer(kind = 4), dimension(3) :: n

      ! Degree of approximation
      integer(kind = 4), dimension(3) :: p

      ! Knot vector
      real (kind = 8), allocatable, dimension(:) :: Ux
      real (kind = 8), allocatable, dimension(:) :: Uy
      real (kind = 8), allocatable, dimension(:) :: Uz

      ! Number of subintervals (currently n - p + 1)
      integer(kind = 4), dimension(3) :: nelem

      ! Size of slices of domain in each dimension
      integer(kind = 4), allocatable, dimension(:) :: dimensionsX
      integer(kind = 4), allocatable, dimension(:) :: dimensionsY
      integer(kind = 4), allocatable, dimension(:) :: dimensionsZ

      ! Offsets of slices of domain in each direction
      integer(kind = 4), allocatable, dimension(:) :: shiftsX
      integer(kind = 4), allocatable, dimension(:) :: shiftsY
      integer(kind = 4), allocatable, dimension(:) :: shiftsZ

      ! Pivot array for these processes that need to solve systems
      integer(kind = 4), allocatable, dimension(:) :: IPIVx
      integer(kind = 4), allocatable, dimension(:) :: IPIVy
      integer(kind = 4), allocatable, dimension(:) :: IPIVz

      ! Number of lower and upper diagonal entries in mass matrix
      integer(kind = 4), dimension(3) :: KL, KU

      ! Number of columns (average) per processor
      integer(kind = 4), dimension(3) :: nrcpp

      ! Range of piece of domain assigned to this process
      integer(kind = 4), dimension(3) :: ibeg, iend

      ! Size of piece of domain assigned to this process
      integer(kind = 4), dimension(3) :: s

      ! Ranges of pieces of domain around the one assigned to this process
      integer(kind = 4), dimension(3) :: ibegsx, iendsx
      integer(kind = 4), dimension(3) :: ibegsy, iendsy
      integer(kind = 4), dimension(3) :: ibegsz, iendsz

      ! Range of elements associated with basis functions assigned to this process
      integer(kind = 4), dimension(3) :: mine, maxe
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      integer(kind = 4), dimension(3) :: m
      
      integer(kind = 4), dimension(3) :: ng
      
      integer(kind = 4), allocatable, dimension(:) :: Ox, Oy, Oz
      
      real (kind = 8), allocatable, dimension(:) :: Jx, Jy, Jz
      
      real (kind = 8), allocatable, dimension(:,:) :: Xx, Xy, Xz
      
      real (kind = 8), allocatable, dimension(:,:,:,:) :: NNx, NNy, NNz
      
      real (kind = 8), allocatable, dimension(:) :: Wx, Wy, Wz
      
      integer (kind = 4), dimension(3) :: lnelem
      
      ! time step length
      real (kind = 8) :: tau
   end type ADS_setup


   type ADS_compute_data
      ! Mass matrix
      real (kind = 8), allocatable, dimension(:,:) :: Mx
      real (kind = 8), allocatable, dimension(:,:) :: My
      real (kind = 8), allocatable, dimension(:,:) :: Mz

      real (kind = 8), allocatable, dimension(:,:) :: F, F2, F3
      real (kind = 8), allocatable, dimension(:,:) :: F_out, F2_out, F3_out

      ! Buffer for coefficients of solution corresponding to neighbouring
      ! parts of the domain. It is (Nx*Ny*Nz) x 3 x 3 x 3 array, where
      ! Nx*Ny*Nz is the size of part of solution for one fragment of domain.
      real (kind = 8), allocatable, dimension(:,:,:,:) :: R
      
      ! current time step
      real (kind = 8) :: t
   end type ADS_compute_data
            
end module Setup