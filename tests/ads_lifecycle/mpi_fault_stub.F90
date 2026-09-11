module mpi
   implicit none

   integer, parameter :: MPI_COMM_WORLD = 91
   integer, parameter :: MPI_INTEGER = 92
   integer, parameter :: MPI_STATUS_SIZE = 5
   integer, parameter :: ALLGATHER_FAILURE_SENTINEL = 7319

   integer :: allgather_calls = 0

contains

   subroutine reset_mpi_fault_stub()
      allgather_calls = 0
   end subroutine reset_mpi_fault_stub


   subroutine mpi_barrier(comm, ierr)
      integer, intent(in) :: comm
      integer, intent(out) :: ierr

      if (comm == MPI_COMM_WORLD) continue
      ierr = 0
   end subroutine mpi_barrier


   subroutine mpi_allgather(sendbuf, sendcount, sendtype, recvbuf, &
                            recvcount, recvtype, comm, ierr)
      integer, intent(in) :: sendbuf(*)
      integer, intent(in) :: sendcount, sendtype, recvcount, recvtype, comm
      integer, intent(out) :: recvbuf(*)
      integer, intent(out) :: ierr

      allgather_calls = allgather_calls + 1

      ! A failed collective does not promise a valid receive buffer.  Keep
      ! the injected contents bounded so a broken caller can be diagnosed
      ! without turning the regression test into an accidental huge
      ! allocation.
      recvbuf(1:3) = 9
      recvbuf(4:6) = 8
      recvbuf(7:9) = 9
      recvbuf(10:12) = 8

      if (sendcount == 12 .and. recvcount == 12 .and. &
          sendtype == MPI_INTEGER .and. recvtype == MPI_INTEGER .and. &
          comm == MPI_COMM_WORLD .and. sendbuf(1) == sendbuf(1)) continue
      ierr = ALLGATHER_FAILURE_SENTINEL
   end subroutine mpi_allgather


   subroutine mpi_abort(comm, errorcode, ierr)
      integer, intent(in) :: comm, errorcode
      integer, intent(out) :: ierr

      if (comm == MPI_COMM_WORLD .and. errorcode /= 0) continue
      ierr = 0
   end subroutine mpi_abort

end module mpi
