module mpi
   implicit none

   integer(kind=4), parameter :: MPI_DOUBLE_PRECISION = 17
   integer(kind=4), parameter :: MPI_SUM = 23
   integer(kind=4), parameter :: MPI_COMM_WORLD = 29
   integer(kind=4) :: reduce_calls = 0
   integer(kind=4) :: recorded_count = 0
   integer(kind=4) :: recorded_datatype = 0
   integer(kind=4) :: recorded_operation = 0
   integer(kind=4) :: recorded_root = -1
   integer(kind=4) :: recorded_communicator = 0
   real(kind=8) :: recorded_send_value = 0.d0

contains

   subroutine MPI_Reduce(sendbuf, recvbuf, count, datatype, operation, root, communicator, ierr)
      real(kind=8), intent(in) :: sendbuf
      real(kind=8), intent(out) :: recvbuf
      integer(kind=4), intent(in) :: count, datatype, operation, root, communicator
      integer(kind=4), intent(out) :: ierr

      reduce_calls = reduce_calls + 1
      recorded_send_value = sendbuf
      recorded_count = count
      recorded_datatype = datatype
      recorded_operation = operation
      recorded_root = root
      recorded_communicator = communicator
      recvbuf = sendbuf
      ierr = 0
   end subroutine MPI_Reduce

end module mpi
