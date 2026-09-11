module reorderRHS
   implicit none

   integer(kind=4) :: reorder_calls = 0

contains

   subroutine reset_reorder_spy()
      reorder_calls = 0
   end subroutine reset_reorder_spy


   subroutine ReorderRHSForX(ibeg, iend, input, output)
      integer(kind=4), intent(in) :: ibeg(3), iend(3)
      real(kind=8), intent(in) :: input(:, :)
      real(kind=8), allocatable, intent(inout) :: output(:, :)

      reorder_calls = reorder_calls + 1
      output = input
   end subroutine ReorderRHSForX


   subroutine ReorderRHSForY(ibeg, iend, input, output)
      integer(kind=4), intent(in) :: ibeg(3), iend(3)
      real(kind=8), intent(in) :: input(:, :)
      real(kind=8), allocatable, intent(inout) :: output(:, :)

      reorder_calls = reorder_calls + 1
      output = input
   end subroutine ReorderRHSForY


   subroutine ReorderRHSForZ(ibeg, iend, input, output)
      integer(kind=4), intent(in) :: ibeg(3), iend(3)
      real(kind=8), intent(in) :: input(:, :)
      real(kind=8), allocatable, intent(inout) :: output(:, :)

      reorder_calls = reorder_calls + 1
      output = input
   end subroutine ReorderRHSForZ

end module reorderRHS
