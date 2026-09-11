module mpi
   implicit none

   integer(kind=4), parameter :: MPI_SUCCESS = 0
   integer(kind=4), parameter :: MPI_COMM_WORLD = 700
   integer(kind=4), parameter :: MPI_2INTEGER = 701
   integer(kind=4), parameter :: MPI_MINLOC = 702
   integer(kind=4), parameter :: MPI_INTEGER = 703
   integer(kind=4), parameter :: FAIL_ON_LOCAL_ERROR = -1

   integer(kind=4) :: allreduce_calls = 0
   integer(kind=4) :: bcast_calls = 0
   integer(kind=4) :: abort_calls = 0
   integer(kind=4) :: failing_allreduce_call = 0
   integer(kind=4) :: allreduce_error = MPI_SUCCESS
   integer(kind=4) :: bcast_error = MPI_SUCCESS
   integer(kind=4) :: abort_code = MPI_SUCCESS
   integer(kind=4) :: abort_comm = -1
   integer(kind=4) :: bcast_root = -1
   logical :: collective_contract_ok = .true.

contains

   subroutine configure_collective_failures(allreduce_call, reduction_status, &
                                             broadcast_status)
      integer(kind=4), intent(in) :: allreduce_call
      integer(kind=4), intent(in) :: reduction_status, broadcast_status

      allreduce_calls = 0
      bcast_calls = 0
      abort_calls = 0
      failing_allreduce_call = allreduce_call
      allreduce_error = reduction_status
      bcast_error = broadcast_status
      abort_code = MPI_SUCCESS
      abort_comm = -1
      bcast_root = -1
      collective_contract_ok = .true.
   end subroutine configure_collective_failures


   subroutine MPI_Allreduce(sendbuf, recvbuf, count, datatype, operation, &
                            comm, ierr)
      integer(kind=4), intent(in) :: sendbuf(2)
      integer(kind=4), intent(out) :: recvbuf(2)
      integer(kind=4), intent(in) :: count, datatype, operation, comm
      integer(kind=4), intent(out) :: ierr

      allreduce_calls = allreduce_calls + 1
      collective_contract_ok = collective_contract_ok .and. count == 1
      collective_contract_ok = collective_contract_ok .and. &
                               datatype == MPI_2INTEGER
      collective_contract_ok = collective_contract_ok .and. &
                               operation == MPI_MINLOC
      collective_contract_ok = collective_contract_ok .and. &
                               comm == MPI_COMM_WORLD

      ierr = MPI_SUCCESS
      if (allreduce_calls == failing_allreduce_call .or. &
          (failing_allreduce_call == FAIL_ON_LOCAL_ERROR .and. &
           sendbuf(1) == 0)) then
         ierr = allreduce_error
         return
      end if
      recvbuf = sendbuf
   end subroutine MPI_Allreduce


   subroutine MPI_Bcast(buffer, count, datatype, root, comm, ierr)
      integer(kind=4), intent(inout) :: buffer
      integer(kind=4), intent(in) :: count, datatype, root, comm
      integer(kind=4), intent(out) :: ierr

      bcast_calls = bcast_calls + 1
      bcast_root = root
      collective_contract_ok = collective_contract_ok .and. count == 1
      collective_contract_ok = collective_contract_ok .and. &
                               datatype == MPI_INTEGER
      collective_contract_ok = collective_contract_ok .and. root == 0
      collective_contract_ok = collective_contract_ok .and. &
                               comm == MPI_COMM_WORLD
      ierr = bcast_error
   end subroutine MPI_Bcast


   subroutine MPI_Abort(comm, errorcode, ierr)
      integer(kind=4), intent(in) :: comm, errorcode
      integer(kind=4), intent(out) :: ierr

      abort_calls = abort_calls + 1
      abort_comm = comm
      abort_code = errorcode
      ierr = MPI_SUCCESS
   end subroutine MPI_Abort

end module mpi
