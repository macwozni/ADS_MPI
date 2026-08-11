module mpi
   use parallelism, only: MYRANK
   implicit none

   integer(kind=4), parameter :: MPI_DOUBLE_PRECISION = 1
   integer(kind=4), parameter :: MPI_CHARACTER = 2
   integer(kind=4), parameter :: MPI_COMM_WORLD = 3

   integer(kind=4) :: real_bcast_calls = 0
   integer(kind=4) :: character_bcast_calls = 0
   logical :: bcast_contract_ok = .true.
   character(len=20) :: broadcast_filename = ''

   interface MPI_Bcast
      module procedure bcast_real_3d
      module procedure bcast_character
   end interface MPI_Bcast

contains

   subroutine reset_mpi_stub()
      real_bcast_calls = 0
      character_bcast_calls = 0
      bcast_contract_ok = .true.
      broadcast_filename = ''
   end subroutine reset_mpi_stub


   subroutine bcast_real_3d(buffer, count, datatype, root, comm, ierr)
      real(kind=8), intent(inout) :: buffer(:, :, :)
      integer(kind=4), intent(in) :: count, datatype, root, comm
      integer(kind=4), intent(out) :: ierr
      integer :: ix, iy, iz

      real_bcast_calls = real_bcast_calls + 1
      bcast_contract_ok = bcast_contract_ok .and. &
                          count == size(buffer) .and. &
                          datatype == MPI_DOUBLE_PRECISION .and. &
                          root == 0 .and. comm == MPI_COMM_WORLD
      if (MYRANK /= root) then
         do iz = 1, size(buffer, 3)
            do iy = 1, size(buffer, 2)
               do ix = 1, size(buffer, 1)
                  buffer(ix, iy, iz) = coefficient_value(ix - 1, iy - 1, iz - 1)
               end do
            end do
         end do
      end if
      ierr = 0
   end subroutine bcast_real_3d


   subroutine bcast_character(buffer, count, datatype, root, comm, ierr)
      character(len=*), intent(inout) :: buffer
      integer(kind=4), intent(in) :: count, datatype, root, comm
      integer(kind=4), intent(out) :: ierr

      character_bcast_calls = character_bcast_calls + 1
      bcast_contract_ok = bcast_contract_ok .and. &
                          count == len(buffer) .and. &
                          datatype == MPI_CHARACTER .and. &
                          root == 0 .and. comm == MPI_COMM_WORLD
      if (MYRANK /= root) buffer = broadcast_filename
      ierr = 0
   end subroutine bcast_character


   function coefficient_value(ix, iy, iz) result(value)
      integer, intent(in) :: ix, iy, iz
      real(kind=8) :: value

      value = real(100*ix + 10*iy + iz, kind=8)
   end function coefficient_value

end module mpi
