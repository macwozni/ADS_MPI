!------------------------------------------------------------------------------
! AGH University of Science and Technology, Krakow, Poland
!------------------------------------------------------------------------------
!
! MODULE: Setup
!
!> @author
!> Maciej Wozniak
!
! DESCRIPTION:
!> This module contains structures that hold configuration of ADS_MPI and current state of computations.
!
! REVISION HISTORY:
! 13 10 2020 - Initial Version
! 
!------------------------------------------------------------------------------

module Setup

implicit none


!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief Structure containing configuration of the solver.
!> Should not be modified by a user.
!---------------------------------------------------------------------------  
type ADS_setup
!> Number of functions in each dimension minus one
   integer(kind = 4), dimension(3) :: n

!> Degree of approximation
   integer(kind = 4), dimension(3) :: p

!> Knot vector in x direction
   real (kind = 8), allocatable, dimension(:) :: Ux
!> Knot vector in y direction
   real (kind = 8), allocatable, dimension(:) :: Uy
!> Knot vector in z direction
   real (kind = 8), allocatable, dimension(:) :: Uz

!> Number of subintervals (currently \f$ n - p + 1\f$)
   integer(kind = 4), dimension(3) :: nelem

!> Size of slices of domain in x dimension
   integer(kind = 4), allocatable, dimension(:) :: dimensionsX
!> Size of slices of domain in y dimension
   integer(kind = 4), allocatable, dimension(:) :: dimensionsY
!> Size of slices of domain in z dimension
   integer(kind = 4), allocatable, dimension(:) :: dimensionsZ

!> Offsets of slices of domain in x direction
   integer(kind = 4), allocatable, dimension(:) :: shiftsX
!> Offsets of slices of domain in y direction
   integer(kind = 4), allocatable, dimension(:) :: shiftsY
!> Offsets of slices of domain in z direction
   integer(kind = 4), allocatable, dimension(:) :: shiftsZ

!> Pivot array for these processes that need to solve systems (x direction)
   integer(kind = 4), allocatable, dimension(:) :: IPIVx
!> Pivot array for these processes that need to solve systems (y direction)
   integer(kind = 4), allocatable, dimension(:) :: IPIVy
!> Pivot array for these processes that need to solve systems (z direction)
   integer(kind = 4), allocatable, dimension(:) :: IPIVz

!> Number of columns (average) per processor
   integer(kind = 4), dimension(3) :: nrcpp

!> Range of piece of domain assigned to this process (beginning)
   integer(kind = 4), dimension(3) :: ibeg
!> Range of piece of domain assigned to this process (end)
   integer(kind = 4), dimension(3) :: iend

!> Size of piece of domain assigned to this process
   integer(kind = 4), dimension(3) :: s

!> Ranges of pieces of domain around the one assigned to this process
   integer(kind = 4), dimension(3) :: ibegsx, iendsx
   integer(kind = 4), dimension(3) :: ibegsy, iendsy
   integer(kind = 4), dimension(3) :: ibegsz, iendsz

!> Range of elements associated with basis functions assigned to this process
   integer(kind = 4), dimension(3) :: mine, maxe

!> Index of the last node in knot vector (number of nodes - 1)
   integer(kind = 4), dimension(3) :: m

!> Order of highest derivatives we want to compute
   integer(kind = 4), dimension(3) :: ng

!> Indexes of first nonzero functions on each element
   integer(kind = 4), allocatable, dimension(:) :: Ox, Oy, Oz

!> Values of the Jacobian of elements
   real (kind = 8), allocatable, dimension(:) :: Jx, Jy, Jz

!> Points of Gauss quadrature
   real (kind = 8), allocatable, dimension(:,:) :: Xx, Xy, Xz

!> Values of \f$ (p+1) \f$ nonzero basis functions and their derivatives at points of Gauss quadrature
   real (kind = 8), allocatable, dimension(:,:,:,:) :: NNx, NNy, NNz

!> Weights of Gauss quadrature points
   real (kind = 8), allocatable, dimension(:) :: Wx, Wy, Wz

!> Local number of subintervals (currently \f$ maxe - mine + 1\f$)
   integer (kind = 4), dimension(3) :: lnelem

! Time step length
   real (kind = 8) :: tau
end type ADS_setup

!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief Structure containing current tate of simulation.
!> Should not be modified by a user.
!--------------------------------------------------------------------------- 
type ADS_compute_data
   real (kind = 8), allocatable, dimension(:,:) :: F, F2, F3
   real (kind = 8), allocatable, dimension(:,:) :: F_out, F2_out, F3_out

!> Buffer for coefficients of solution corresponding to neighbouring
!> parts of the domain. It is \f$ (N_x \times N_y \times N_z) x 3 x 3 x 3 \f$ array, where
!> \f$ N_x \times N_y \times N_z \f$ is the size of part of solution for one fragment of domain.
   real (kind = 8), allocatable, dimension(:,:,:,:) :: R

!> current time step
   real (kind = 8) :: t
   
!> previous solution coefficient at given point
   real (kind = 8), allocatable, dimension(:,:,:,:,:,:) :: Un
!> previous solution coefficient at given point   
   real (kind = 8), allocatable, dimension(:,:,:,:,:,:) :: Un13
!> previous solution coefficient at given point
   real (kind = 8), allocatable, dimension(:,:,:,:,:,:) :: Un23
!> value of derivative from previous time step
   real (kind = 8), allocatable, dimension(:,:,:,:,:,:,:) :: dUn
end type ADS_compute_data
   
end module Setup