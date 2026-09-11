module mpi
   implicit none

   integer(kind=4), parameter :: MPI_SUCCESS = 0
   integer(kind=4), parameter :: MPI_DOUBLE_PRECISION = 1
   integer(kind=4), parameter :: MPI_COMM_WORLD = 2
   integer(kind=4), parameter :: MPI_STATUS_SIZE = 3

   integer(kind=4), parameter :: FAIL_NONE = 0
   integer(kind=4), parameter :: FAIL_IRECV = 1
   integer(kind=4), parameter :: FAIL_ISEND = 2
   integer(kind=4), parameter :: FAIL_WAITALL = 3

   integer(kind=4), parameter :: IRECV_SENTINEL = 9101
   integer(kind=4), parameter :: ISEND_SENTINEL = 9202
   integer(kind=4), parameter :: WAITALL_SENTINEL = 9303

   integer(kind=4) :: failure_mode = FAIL_NONE
   integer(kind=4) :: first_mpi_error = MPI_SUCCESS
   integer(kind=4) :: irecv_calls = 0
   integer(kind=4) :: isend_calls = 0
   integer(kind=4) :: waitall_calls = 0
   integer(kind=4) :: transport_posts_after_failure = 0
   logical :: failure_seen = .false.

   interface mpi_gatherv
      module procedure mpi_gatherv_rank_one
      module procedure mpi_gatherv_rank_two
   end interface mpi_gatherv

contains

   subroutine configure_mpi_failure(mode)
      integer(kind=4), intent(in) :: mode

      failure_mode = mode
      first_mpi_error = MPI_SUCCESS
      irecv_calls = 0
      isend_calls = 0
      waitall_calls = 0
      transport_posts_after_failure = 0
      failure_seen = .false.
   end subroutine configure_mpi_failure


   integer(kind=4) function injected_sentinel() result(status)
      select case (failure_mode)
      case (FAIL_IRECV)
         status = IRECV_SENTINEL
      case (FAIL_ISEND)
         status = ISEND_SENTINEL
      case (FAIL_WAITALL)
         status = WAITALL_SENTINEL
      case default
         status = MPI_SUCCESS
      end select
   end function injected_sentinel


   subroutine note_post_after_failure()
      if (failure_seen) then
         transport_posts_after_failure = transport_posts_after_failure + 1
      end if
   end subroutine note_post_after_failure


   subroutine inject_failure(status)
      integer(kind=4), intent(out) :: status

      status = injected_sentinel()
      first_mpi_error = status
      failure_seen = .true.
   end subroutine inject_failure


   subroutine mpi_irecv(buffer, count, datatype, source, tag, comm, request, ierr)
      real(kind=8), dimension(:), intent(inout) :: buffer
      integer(kind=4), intent(in) :: count, datatype, source, tag, comm
      integer(kind=4), intent(out) :: request, ierr

      call note_post_after_failure()
      irecv_calls = irecv_calls + 1
      request = 100 + irecv_calls
      if (failure_mode == FAIL_IRECV .and. irecv_calls == 1) then
         call inject_failure(ierr)
         return
      end if

      ierr = MPI_SUCCESS
      if (count > 0) buffer(1:count) = 100.0d0 + dble(source)
      if (datatype + tag + comm == -huge(0)) request = -1
   end subroutine mpi_irecv


   subroutine mpi_isend(buffer, count, datatype, destination, tag, comm, request, ierr)
      real(kind=8), dimension(:), intent(in) :: buffer
      integer(kind=4), intent(in) :: count, datatype, destination, tag, comm
      integer(kind=4), intent(out) :: request, ierr

      call note_post_after_failure()
      isend_calls = isend_calls + 1
      request = 200 + isend_calls
      if (failure_mode == FAIL_ISEND .and. isend_calls == 1) then
         call inject_failure(ierr)
         return
      end if

      ierr = MPI_SUCCESS
      if (count + datatype + destination + tag + comm == -huge(0)) request = nint(buffer(1))
   end subroutine mpi_isend


   subroutine mpi_waitall(count, requests, statuses, ierr)
      integer(kind=4), intent(in) :: count
      integer(kind=4), dimension(:), intent(inout) :: requests
      integer(kind=4), dimension(:, :), intent(out) :: statuses
      integer(kind=4), intent(out) :: ierr

      waitall_calls = waitall_calls + 1
      statuses = 0
      if (failure_mode == FAIL_WAITALL .and. waitall_calls == 1) then
         call inject_failure(ierr)
         return
      end if

      ierr = MPI_SUCCESS
      if (count == -huge(0)) requests(1) = -1
   end subroutine mpi_waitall


   subroutine mpi_gatherv_rank_one(send_buffer, send_count, send_type, recv_buffer, &
                                   recv_counts, displacements, recv_type, root, comm, ierr)
      real(kind=8), dimension(:), intent(in) :: send_buffer
      integer(kind=4), intent(in) :: send_count, send_type, recv_type, root, comm
      real(kind=8), dimension(:), intent(out) :: recv_buffer
      integer(kind=4), dimension(:), intent(in) :: recv_counts, displacements
      integer(kind=4), intent(out) :: ierr

      ierr = MPI_SUCCESS
      if (size(recv_buffer) > 0 .and. size(send_buffer) > 0) recv_buffer(1) = send_buffer(1)
      if (send_count + send_type + recv_type + root + comm + &
          size(recv_counts) + size(displacements) == -huge(0)) ierr = 1
   end subroutine mpi_gatherv_rank_one


   subroutine mpi_gatherv_rank_two(send_buffer, send_count, send_type, recv_buffer, &
                                   recv_counts, displacements, recv_type, root, comm, ierr)
      real(kind=8), dimension(:, :), intent(in) :: send_buffer
      integer(kind=4), intent(in) :: send_count, send_type, recv_type, root, comm
      real(kind=8), dimension(:), intent(out) :: recv_buffer
      integer(kind=4), dimension(:), intent(in) :: recv_counts, displacements
      integer(kind=4), intent(out) :: ierr

      ierr = MPI_SUCCESS
      if (size(recv_buffer) > 0 .and. size(send_buffer) > 0) recv_buffer(1) = send_buffer(1, 1)
      if (send_count + send_type + recv_type + root + comm + &
          size(recv_counts) + size(displacements) == -huge(0)) ierr = 1
   end subroutine mpi_gatherv_rank_two


   subroutine mpi_scatterv(send_buffer, send_counts, displacements, send_type, &
                           recv_buffer, recv_count, recv_type, root, comm, ierr)
      real(kind=8), dimension(:), intent(in) :: send_buffer
      integer(kind=4), dimension(:), intent(in) :: send_counts, displacements
      integer(kind=4), intent(in) :: send_type, recv_count, recv_type, root, comm
      real(kind=8), dimension(:), intent(out) :: recv_buffer
      integer(kind=4), intent(out) :: ierr

      ierr = MPI_SUCCESS
      if (size(recv_buffer) > 0 .and. size(send_buffer) > 0) recv_buffer(1) = send_buffer(1)
      if (size(send_counts) + size(displacements) + send_type + recv_count + &
          recv_type + root + comm == -huge(0)) ierr = 1
   end subroutine mpi_scatterv


   subroutine mpi_allgatherv(send_buffer, send_count, send_type, recv_buffer, &
                             recv_counts, displacements, recv_type, comm, ierr)
      real(kind=8), dimension(:), intent(in) :: send_buffer
      integer(kind=4), intent(in) :: send_count, send_type, recv_type, comm
      real(kind=8), dimension(:), intent(out) :: recv_buffer
      integer(kind=4), dimension(:), intent(in) :: recv_counts, displacements
      integer(kind=4), intent(out) :: ierr

      ierr = MPI_SUCCESS
      if (size(recv_buffer) > 0 .and. size(send_buffer) > 0) recv_buffer(1) = send_buffer(1)
      if (send_count + send_type + recv_type + comm + &
          size(recv_counts) + size(displacements) == -huge(0)) ierr = 1
   end subroutine mpi_allgatherv

end module mpi
