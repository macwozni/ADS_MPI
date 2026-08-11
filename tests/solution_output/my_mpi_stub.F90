module my_mpi
   use Setup, only: ADS_Setup
   use parallelism, only: MYRANK
   implicit none

   integer(kind=4) :: gather_calls = 0
   integer(kind=4) :: captured_root = -1
   integer(kind=4), dimension(3) :: captured_n = 0
   integer(kind=4), dimension(3) :: captured_p = 0
   integer(kind=4), dimension(3) :: captured_s = 0
   real(kind=8), allocatable :: captured_part(:, :)

contains

   subroutine reset_my_mpi_stub()
      gather_calls = 0
      captured_root = -1
      captured_n = 0
      captured_p = 0
      captured_s = 0
      if (allocated(captured_part)) deallocate (captured_part)
   end subroutine reset_my_mpi_stub


   subroutine GatherFullSolution(at, part, full, n, p, s)
      integer(kind=4), intent(in) :: at
      real(kind=8), intent(in) :: part(:, :)
      real(kind=8), allocatable, intent(out) :: full(:, :, :)
      integer(kind=4), intent(in) :: n(:), p(:), s(:)
      integer :: ix, iy, iz

      gather_calls = gather_calls + 1
      captured_root = at
      captured_n = n
      captured_p = p
      captured_s = s
      if (allocated(captured_part)) deallocate (captured_part)
      allocate (captured_part(size(part, 1), size(part, 2)))
      captured_part = part

      if (MYRANK == at) then
         allocate (full(0:n(1), 0:n(2), 0:n(3)))
         do iz = 0, n(3)
            do iy = 0, n(2)
               do ix = 0, n(1)
                  full(ix, iy, iz) = coefficient_value(ix, iy, iz)
               end do
            end do
         end do
      end if
   end subroutine GatherFullSolution


   function coefficient_value(ix, iy, iz) result(value)
      integer, intent(in) :: ix, iy, iz
      real(kind=8) :: value

      value = real(100*ix + 10*iy + iz, kind=8)
   end function coefficient_value

end module my_mpi
