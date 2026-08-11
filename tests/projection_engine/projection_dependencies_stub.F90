module igrm_space
   implicit none
contains
   subroutine ValidateIGRMMesh(marker)
      integer, intent(inout) :: marker
      marker = marker + 1
   end subroutine ValidateIGRMMesh
end module igrm_space

module operator_assembly
   implicit none
contains
   subroutine MKBBT_large(marker)
      integer, intent(inout) :: marker
      marker = marker + 2
   end subroutine MKBBT_large

   subroutine MKBBT_small(marker)
      integer, intent(inout) :: marker
      marker = marker + 4
   end subroutine MKBBT_small

   subroutine ComputeMatrix(marker)
      integer, intent(inout) :: marker
      marker = marker + 8
   end subroutine ComputeMatrix
end module operator_assembly

module rhs_assembly
   implicit none
contains
   subroutine Form3DRHS(marker)
      integer, intent(inout) :: marker
      marker = marker + 16
   end subroutine Form3DRHS

   subroutine create_mixed_space(marker)
      integer, intent(inout) :: marker
      marker = marker + 32
   end subroutine create_mixed_space
end module rhs_assembly

module solution_reconstruction
   implicit none
contains
   subroutine FormUn(marker)
      integer, intent(inout) :: marker
      marker = marker + 64
   end subroutine FormUn

   subroutine global2local(marker)
      integer, intent(inout) :: marker
      marker = marker + 128
   end subroutine global2local
end module solution_reconstruction
