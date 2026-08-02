!------------------------------------------------------------------------------
!
! MODULE: Setup
!
! DESCRIPTION:
!> @file Setup.F90
!> @brief Module defining setup and runtime data containers for ADS-based
!> computations.
!>
!> @details
!> This module provides derived types used to store:
!> - discretization metadata,
!> - spline-space definitions,
!> - domain-decomposition data,
!> - quadrature and basis-function tables,
!> - temporary buffers and state variables required during computation.
!>
!> The module is intended as a central repository of configuration and
!> runtime arrays shared by higher-level assembly and solver procedures.
!>
!> The principal derived types are:
!> - \ref ADS_setup, storing static or semi-static problem setup data,
!> - \ref ADS_compute_data, storing time-dependent working arrays and
!>   solution-related buffers.
!
!------------------------------------------------------------------------------
module Setup

   implicit none

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Container holding discretization, quadrature, and
!> domain-decomposition setup data.
!>
!> @details
!> This derived type stores the structural information required to define
!> the approximation space and the local computational environment for a
!> given process.
!>
!> The stored data include:
!> - spline-space sizes and polynomial degrees,
!> - knot vectors in all three parametric directions,
!> - element counts and quadrature metadata,
!> - decomposition of the domain into local and neighbouring subranges,
!> - indexing and pivot data used by local linear-algebra operations,
!> - precomputed basis-function values and quadrature coordinates.
!>
!> Most allocatable components are expected to be initialized during the
!> setup phase before the computational kernels are executed.
!
! Notes:
! ------
!> @note
!> The type does not enforce allocation or consistency constraints by
!> itself. These are assumed to be managed externally by the calling code.
!
!---------------------------------------------------------------------------
type ADS_setup

!> @brief Numbers of basis functions in each direction minus one.
      integer(kind=4), dimension(3) :: n

!> @brief Polynomial degrees of approximation in the three directions.
      integer(kind=4), dimension(3) :: p

!> @brief Knot vector in the first parametric direction.
      real(kind=8), allocatable, dimension(:) :: Ux
!> @brief Knot vector in the second parametric direction.
      real(kind=8), allocatable, dimension(:) :: Uy
!> @brief Knot vector in the third parametric direction.
      real(kind=8), allocatable, dimension(:) :: Uz

!> @brief Numbers of elements in the three directions.
!>
!> @details
!> In the present convention this is currently associated with
!> \f$ n - p + 1 \f$ in each coordinate direction.
      integer(kind=4), dimension(3) :: nelem

!> @brief Sizes of domain slices in the first direction.
      integer(kind=4), allocatable, dimension(:) :: dimensionsX
!> @brief Sizes of domain slices in the second direction.
      integer(kind=4), allocatable, dimension(:) :: dimensionsY
!> @brief Sizes of domain slices in the third direction.
      integer(kind=4), allocatable, dimension(:) :: dimensionsZ

!> @brief Offsets of domain slices in the first direction.
      integer(kind=4), allocatable, dimension(:) :: shiftsX
!> @brief Offsets of domain slices in the second direction.
      integer(kind=4), allocatable, dimension(:) :: shiftsY
!> @brief Offsets of domain slices in the third direction.
      integer(kind=4), allocatable, dimension(:) :: shiftsZ

!> @brief Pivot array for processes solving systems associated with the
!> first direction.
      integer(kind=4), allocatable, dimension(:) :: IPIVx
!> @brief Pivot array for processes solving systems associated with the
!> second direction.
      integer(kind=4), allocatable, dimension(:) :: IPIVy
!> @brief Pivot array for processes solving systems associated with the
!> third direction.
      integer(kind=4), allocatable, dimension(:) :: IPIVz

!> @brief Average numbers of columns per processor in each direction.
      integer(kind=4), dimension(3) :: nrcpp

!> @brief Lower and upper bounds of the domain fragment assigned to the
!> current process.
      integer(kind=4), dimension(3) :: ibeg, iend

!> @brief Sizes of the domain fragment assigned to the current process.
      integer(kind=4), dimension(3) :: s

!> @brief Bounds of neighbouring domain fragments in the first direction.
      integer(kind=4), dimension(3) :: ibegsx, iendsx
!> @brief Bounds of neighbouring domain fragments in the second direction.
      integer(kind=4), dimension(3) :: ibegsy, iendsy
!> @brief Bounds of neighbouring domain fragments in the third direction.
      integer(kind=4), dimension(3) :: ibegsz, iendsz

!> @brief Ranges of elements associated with basis functions assigned to
!> the current process.
      integer(kind=4), dimension(3) :: mine, maxe, lnelem

!> @brief Last valid indices in the knot vectors.
!>
!> @details
!> Equivalently, these are the numbers of knot entries minus one.
      integer(kind=4), dimension(3) :: m

!> @brief Numbers of Gauss quadrature points in the three directions.
      integer(kind=4), dimension(3) :: ng

!> @brief Indices of the first nonzero basis functions on elements in the
!> first direction.
      integer(kind=4), allocatable, dimension(:) :: Ox
!> @brief Indices of the first nonzero basis functions on elements in the
!> second direction.
      integer(kind=4), allocatable, dimension(:) :: Oy
!> @brief Indices of the first nonzero basis functions on elements in the
!> third direction.
      integer(kind=4), allocatable, dimension(:) :: Oz

!> @brief Jacobians of element mappings in the first direction.
      real(kind=8), allocatable, dimension(:) :: Jx
!> @brief Jacobians of element mappings in the second direction.
      real(kind=8), allocatable, dimension(:) :: Jy
!> @brief Jacobians of element mappings in the third direction.
      real(kind=8), allocatable, dimension(:) :: Jz

!> @brief Gauss quadrature points mapped to physical elements in the
!> first direction.
      real(kind=8), allocatable, dimension(:, :) :: Xx
!> @brief Gauss quadrature points mapped to physical elements in the
!> second direction.
      real(kind=8), allocatable, dimension(:, :) :: Xy
!> @brief Gauss quadrature points mapped to physical elements in the
!> third direction.
      real(kind=8), allocatable, dimension(:, :) :: Xz

!> @brief Values of the nonzero basis functions and their derivatives at
!> Gauss points in the first direction.
      real(kind=8), allocatable, dimension(:, :, :, :) :: NNx
!> @brief Values of the nonzero basis functions and their derivatives at
!> Gauss points in the second direction.
      real(kind=8), allocatable, dimension(:, :, :, :) :: NNy
!> @brief Values of the nonzero basis functions and their derivatives at
!> Gauss points in the third direction.
      real(kind=8), allocatable, dimension(:, :, :, :) :: NNz

!> @brief Weights of Gauss quadrature points in the first direction.
      real(kind=8), allocatable, dimension(:) :: Wx
!> @brief Weights of Gauss quadrature points in the second direction.
      real(kind=8), allocatable, dimension(:) :: Wy
!> @brief Weights of Gauss quadrature points in the third direction.
      real(kind=8), allocatable, dimension(:) :: Wz

!> @brief Time-step length.
      real(kind=8) :: tau

end type ADS_setup


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Container holding temporary arrays and time-dependent data used
!> during ADS computations.
!>
!> @details
!> This derived type stores working right-hand-side arrays, auxiliary
!> buffers, neighbouring-domain coefficient blocks, and solution data
!> from current or previous time levels.
!>
!> The components are intended to support tensor-directional assembly,
!> residual handling, and storage of historical solution values and
!> derivatives needed by time-stepping procedures.
!
! Notes:
! ------
!> @note
!> Several arrays are high-dimensional and their exact interpretation
!> depends on the calling algorithm and decomposition layout.
!
!---------------------------------------------------------------------------
type ADS_compute_data

!> @brief Small right-hand-side arrays associated with test functions.
      real(kind=8), allocatable, dimension(:, :) :: Ft, Ft2, Ft3

!> @brief Small right-hand-side arrays associated with trial functions.
      real(kind=8), allocatable, dimension(:, :) :: F, F2, F3

!> @brief Auxiliary right-hand-side remnant buffers.
      real(kind=8), allocatable, dimension(:, :) :: FF, FFt

!> @brief Compact coefficient halo required by the local element range.
!>
!> @details
!> The first three dimensions correspond to the global coefficient box
!> `halo_begin:halo_end`; the final singleton dimension preserves the
!> historical four-dimensional declaration of this component.
      real(kind=8), allocatable, dimension(:, :, :, :) :: R

!> @brief Inclusive global coefficient bounds stored in \ref R.
      integer(kind=4), dimension(3) :: halo_begin, halo_end

!> @brief Inclusive global boxes sent to and received from every MPI rank.
!>
!> @details
!> Column `rank + 1` stores the exact intersection used for that peer.  An
!> empty intersection is represented by an end smaller than its begin.
      integer(kind=4), allocatable, dimension(:, :) :: halo_send_begin, halo_send_end
      integer(kind=4), allocatable, dimension(:, :) :: halo_recv_begin, halo_recv_end

!> @brief Packed-message sizes and zero-based displacements for halo peers.
      integer(kind=4), allocatable, dimension(:) :: halo_send_count, halo_send_displ
      integer(kind=4), allocatable, dimension(:) :: halo_recv_count, halo_recv_displ

!> @brief Reusable packed communication buffers for coefficient exchange.
      real(kind=8), allocatable, dimension(:) :: halo_send_buffer, halo_recv_buffer

!> @brief Reusable MPI request and status storage for halo exchange.
      integer(kind=4), allocatable, dimension(:) :: halo_requests
      integer(kind=4), allocatable, dimension(:, :) :: halo_statuses

!> @brief Current time value or time-step position.
      real(kind=8) :: t

!> @brief Global element range covered by the solution-value buffers.
      integer(kind=4), dimension(3) :: state_mine, state_maxe, state_lnelem

!> @brief Solution coefficients from the previous time level.
      real(kind=8), allocatable, dimension(:, :, :, :, :, :) :: Un

!> @brief Auxiliary solution-coefficient buffer associated with a
!> first/third-directional intermediate state.
      real(kind=8), allocatable, dimension(:, :, :, :, :, :) :: Un13

!> @brief Auxiliary solution-coefficient buffer associated with a
!> second/third-directional intermediate state.
      real(kind=8), allocatable, dimension(:, :, :, :, :, :) :: Un23

!> @brief Derivative values inherited from the previous time step.
      real(kind=8), allocatable, dimension(:, :, :, :, :, :, :) :: dUn

!> @brief Derivative values reconstructed from the previous time level.
      real(kind=8), allocatable, dimension(:, :, :, :, :, :, :) :: dUn0

!> @brief Derivative values reconstructed from the first intermediate state.
      real(kind=8), allocatable, dimension(:, :, :, :, :, :, :) :: dUn13

!> @brief Derivative values reconstructed from the second intermediate state.
      real(kind=8), allocatable, dimension(:, :, :, :, :, :, :) :: dUn23

!> @brief RHS derivative-state selector for time schemes.
!>
!> @details
!> Rows 1:6 match RHS derivative coefficients. Value 0 means "use the
!> currently active derivative buffer", while 1, 2, and 3 select
!> dUn0, dUn13, and dUn23, respectively.
      integer(kind=4), dimension(6, 3) :: rhs_du_state

end type ADS_compute_data

end module Setup
