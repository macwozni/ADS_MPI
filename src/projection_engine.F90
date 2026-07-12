!------------------------------------------------------------------------------
!
! MODULE: projection_engine
!
! DESCRIPTION:
!> @file projection_engine.F90
!> @brief Compatibility facade for projection, assembly, and reconstruction helpers.
!>
!> @details
!> The projection and assembly implementation is split across focused
!> modules:
!> - \ref igrm_space for iGRM space compatibility checks,
!> - \ref operator_assembly for one-dimensional sparse operators,
!> - \ref rhs_assembly for three-dimensional RHS assembly,
!> - \ref solution_reconstruction for quadrature-point reconstruction.
!>
!> This facade preserves the historical `projection_engine` API used by
!> ADS, utilities, and legacy problem-local sources. New code may import
!> directly from the focused modules, but existing call sites can continue
!> to use this module without changing semantics.
!
!------------------------------------------------------------------------------
module projection_engine

   use igrm_space, ONLY: ValidateIGRMMesh
   use operator_assembly, ONLY: MKBBT_large, MKBBT_small, ComputeMatrix
   use rhs_assembly, ONLY: Form3DRHS, create_mixed_space
   use solution_reconstruction, ONLY: FormUn, global2local

   implicit none

   public :: ValidateIGRMMesh
   public :: MKBBT_large, MKBBT_small, ComputeMatrix
   public :: Form3DRHS, create_mixed_space
   public :: FormUn, global2local

end module projection_engine
