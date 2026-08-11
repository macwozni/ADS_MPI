module mpi
   implicit none

   integer(kind=4), parameter :: MPI_SUCCESS = 0
   integer(kind=4), parameter :: MPI_COMM_NULL = 0
   integer(kind=4), parameter :: MPI_GROUP_NULL = 0
   integer(kind=4), parameter :: MPI_COMM_WORLD = 1

   integer(kind=4) :: failure_mode = 0
   integer(kind=4) :: comm_free_calls = 0
   integer(kind=4) :: group_free_calls = 0
   integer(kind=4) :: group_incl_calls = 0
   integer(kind=4) :: comm_create_calls = 0

contains

   subroutine configure_failure(mode)
      integer(kind=4), intent(in) :: mode

      failure_mode = mode
      comm_free_calls = 0
      group_free_calls = 0
      group_incl_calls = 0
      comm_create_calls = 0
   end subroutine configure_failure


   subroutine mpi_comm_group(comm, group, ierr)
      integer(kind=4), intent(in) :: comm
      integer(kind=4), intent(out) :: group, ierr

      group = 900 + 0*comm
      ierr = MPI_SUCCESS
      if (failure_mode == 1) ierr = 101
   end subroutine mpi_comm_group


   subroutine mpi_barrier(comm, ierr)
      integer(kind=4), intent(in) :: comm
      integer(kind=4), intent(out) :: ierr

      ierr = MPI_SUCCESS + 0*comm
   end subroutine mpi_barrier


   subroutine mpi_group_incl(parent, n, ranks, group, ierr)
      integer(kind=4), intent(in) :: parent, n
      integer(kind=4), intent(in) :: ranks(:)
      integer(kind=4), intent(out) :: group, ierr

      group_incl_calls = group_incl_calls + 1
      group = 100*n + group_incl_calls + 0*parent + 0*sum(ranks)
      ierr = MPI_SUCCESS
      if (failure_mode == 2 .and. n == 4) ierr = 102
      if (failure_mode == 3 .and. n == 3) ierr = 103
      if (failure_mode == 4 .and. n == 2) ierr = 104
   end subroutine mpi_group_incl


   subroutine mpi_group_free(group, ierr)
      integer(kind=4), intent(inout) :: group
      integer(kind=4), intent(out) :: ierr

      group_free_calls = group_free_calls + 1
      ierr = MPI_SUCCESS
      if (failure_mode == 5 .and. group == 900) ierr = 105
      if (failure_mode == 20 .and. group_free_calls == 1) ierr = 202
      group = MPI_GROUP_NULL
   end subroutine mpi_group_free


   subroutine mpi_comm_create(parent, group, comm, ierr)
      integer(kind=4), intent(in) :: parent, group
      integer(kind=4), intent(out) :: comm, ierr
      integer(kind=4) :: axis_size

      comm_create_calls = comm_create_calls + 1
      comm = 1000 + comm_create_calls + 0*parent
      axis_size = group/100
      ierr = MPI_SUCCESS
      if (failure_mode == 6 .and. axis_size == 4) ierr = 106
      if (failure_mode == 7 .and. axis_size == 3) ierr = 107
      if (failure_mode == 8 .and. axis_size == 2) ierr = 108
   end subroutine mpi_comm_create


   subroutine mpi_comm_free(comm, ierr)
      integer(kind=4), intent(inout) :: comm
      integer(kind=4), intent(out) :: ierr

      comm_free_calls = comm_free_calls + 1
      ierr = MPI_SUCCESS
      if (failure_mode == 20 .and. comm_free_calls == 1) ierr = 201
      comm = MPI_COMM_NULL
   end subroutine mpi_comm_free

end module mpi
