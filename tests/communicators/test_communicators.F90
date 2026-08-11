program test_communicators
   use communicators, only: COMMX, COMMY, COMMZ, Cleanup_Communicators, &
      CreateCommunicators, processors
   use mpi
   use parallelism, only: Cleanup_Parallelism, InitializeParallelism, &
      MYRANK, MYRANKX, MYRANKY, MYRANKZ, NRPROC, NRPROCX, NRPROCY, NRPROCZ
   implicit none

   integer(kind=4), dimension(3) :: process_grid
   integer(kind=4) :: checks, failures, ierr

   call read_process_grid(process_grid)
   call InitializeParallelism(process_grid(1), process_grid(2), &
      process_grid(3), ierr)

   checks = 0
   failures = 0
   call assert_global('InitializeParallelism succeeds', &
      ierr == MPI_SUCCESS, failures)

   call Cleanup_Communicators(ierr)
   call assert_global('Cleanup_Communicators is safe before first creation', &
      ierr == MPI_SUCCESS .and. COMMX == MPI_COMM_NULL .and. &
      COMMY == MPI_COMM_NULL .and. COMMZ == MPI_COMM_NULL, failures)

   call CreateCommunicators(ierr)
   call assert_global('first CreateCommunicators succeeds', &
      ierr == MPI_SUCCESS, failures)
   call validate_communicators('first creation', failures)

   call Cleanup_Communicators(ierr)
   call assert_global('first Cleanup_Communicators succeeds', &
      ierr == MPI_SUCCESS, failures)
   call assert_global('first cleanup nulls public handles', &
      COMMX == MPI_COMM_NULL .and. COMMY == MPI_COMM_NULL .and. &
      COMMZ == MPI_COMM_NULL, failures)

   call CreateCommunicators(ierr)
   call assert_global('second CreateCommunicators succeeds', &
      ierr == MPI_SUCCESS, failures)
   call validate_communicators('second creation', failures)

   call Cleanup_Communicators(ierr)
   call assert_global('second Cleanup_Communicators succeeds', &
      ierr == MPI_SUCCESS, failures)
   call assert_global('second cleanup nulls public handles', &
      COMMX == MPI_COMM_NULL .and. COMMY == MPI_COMM_NULL .and. &
      COMMZ == MPI_COMM_NULL, failures)

   call Cleanup_Communicators(ierr)
   call assert_global('Cleanup_Communicators is idempotent', &
      ierr == MPI_SUCCESS .and. COMMX == MPI_COMM_NULL .and. &
      COMMY == MPI_COMM_NULL .and. COMMZ == MPI_COMM_NULL, failures)

   if (MYRANK == 0) then
      if (failures == 0) then
         write(*, '(A,I0,A)') 'OK (', checks, ' communicator checks)'
      else
         write(*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
            ' communicator checks)'
      end if
   end if

   call Cleanup_Parallelism(ierr)
   if (failures /= 0) stop 1

contains

   subroutine validate_communicators(stage, failure_count)
      character(len=*), intent(in) :: stage
      integer(kind=4), intent(inout) :: failure_count
      logical :: handles_ready

      call test_processors(stage, failure_count)

      call assert_global(trim(stage)//': local handles are non-null', &
         COMMX /= MPI_COMM_NULL .and. COMMY /= MPI_COMM_NULL .and. &
         COMMZ /= MPI_COMM_NULL, failure_count, handles_ready)

      if (.not. handles_ready) return

      call test_axis(stage, 'X', 1, COMMX, NRPROCX, MYRANKX, failure_count)
      call test_axis(stage, 'Y', 2, COMMY, NRPROCY, MYRANKY, failure_count)
      call test_axis(stage, 'Z', 3, COMMZ, NRPROCZ, MYRANKZ, failure_count)
   end subroutine validate_communicators


   subroutine test_processors(stage, failure_count)
      character(len=*), intent(in) :: stage
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4) :: expected_rank, i, j, k
      logical :: local_ok

      local_ok = .true.
      do k = 1, NRPROCZ
         do j = 1, NRPROCY
            do i = 1, NRPROCX
               expected_rank = (i - 1) + (j - 1)*NRPROCX + &
                  (k - 1)*NRPROCX*NRPROCY
               local_ok = local_ok .and. &
                  processors(i, j, k) == expected_rank
            end do
         end do
      end do

      call assert_global(trim(stage)//': processors table is exact', &
         local_ok, failure_count)
   end subroutine test_processors


   subroutine test_axis(stage, axis_name, axis, comm, expected_size, &
      expected_rank, failure_count)
      character(len=*), intent(in) :: stage, axis_name
      integer(kind=4), intent(in) :: axis, comm, expected_size, expected_rank
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4), allocatable :: gathered(:), expected(:)
      integer(kind=4) :: comm_rank, comm_size, comm_ierr, i
      logical :: local_ok

      call mpi_comm_size(comm, comm_size, comm_ierr)
      local_ok = comm_ierr == MPI_SUCCESS .and. comm_size == expected_size
      call mpi_comm_rank(comm, comm_rank, comm_ierr)
      local_ok = local_ok .and. comm_ierr == MPI_SUCCESS .and. &
         comm_rank == expected_rank
      call assert_global(trim(stage)//': COMM'//axis_name// &
         ' size and rank match axis', local_ok, failure_count)

      allocate(gathered(expected_size), expected(expected_size))
      call mpi_allgather(MYRANK, 1, MPI_INTEGER, gathered, 1, MPI_INTEGER, &
         comm, comm_ierr)

      do i = 0, expected_size - 1
         select case (axis)
         case (1)
            expected(i + 1) = i + MYRANKY*NRPROCX + &
               MYRANKZ*NRPROCX*NRPROCY
         case (2)
            expected(i + 1) = MYRANKX + i*NRPROCX + &
               MYRANKZ*NRPROCX*NRPROCY
         case (3)
            expected(i + 1) = MYRANKX + MYRANKY*NRPROCX + &
               i*NRPROCX*NRPROCY
         end select
      end do

      local_ok = comm_ierr == MPI_SUCCESS .and. all(gathered == expected)
      call assert_global(trim(stage)//': COMM'//axis_name// &
         ' gathers exact world-rank fibre', local_ok, failure_count)

      deallocate(gathered, expected)
   end subroutine test_axis


   subroutine assert_global(label, local_ok, failure_count, result)
      character(len=*), intent(in) :: label
      logical, intent(in) :: local_ok
      integer(kind=4), intent(inout) :: failure_count
      logical, intent(out), optional :: result
      integer(kind=4) :: global_flag, local_flag, reduce_ierr
      logical :: global_ok

      local_flag = merge(1, 0, local_ok)
      call mpi_allreduce(local_flag, global_flag, 1, MPI_INTEGER, MPI_MIN, &
         MPI_COMM_WORLD, reduce_ierr)
      global_ok = reduce_ierr == MPI_SUCCESS .and. global_flag == 1
      if (present(result)) result = global_ok

      checks = checks + 1
      if (.not. global_ok) failure_count = failure_count + 1

      if (MYRANK == 0) then
         if (global_ok) then
            write(*, '(A,A)') 'PASS ', trim(label)
         else
            write(*, '(A,A)') 'FAIL ', trim(label)
         end if
      end if
   end subroutine assert_global


   subroutine read_process_grid(grid)
      integer(kind=4), intent(out), dimension(3) :: grid
      character(len=32) :: argument
      integer(kind=4) :: i, read_status

      if (command_argument_count() /= 3) then
         write(*, '(A)') &
            'Usage: test_communicators <procx> <procy> <procz>'
         stop 2
      end if

      do i = 1, 3
         call get_command_argument(i, argument)
         read(argument, *, iostat=read_status) grid(i)
         if (read_status /= 0 .or. grid(i) <= 0) then
            write(*, '(A,I0,A,A)') 'Invalid process-grid argument ', i, &
               ': ', trim(argument)
            stop 2
         end if
      end do
   end subroutine read_process_grid

end program test_communicators
