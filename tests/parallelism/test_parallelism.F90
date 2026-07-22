program test_parallelism
   use mpi
   use parallelism, only: MYRANK, MYRANKX, MYRANKY, MYRANKZ, NRPROC, &
      NRPROCX, NRPROCY, NRPROCZ, PRINTRANK, InitializeParallelism, &
      Cleanup_Parallelism, Decompose, LinearIndex, ComputeEndpoints, &
      FillDimVector
   implicit none

   integer(kind=4), dimension(3) :: process_grid
   integer(kind=4) :: checks, failures, ierr, finalize_ierr
   logical :: finalized, cleanup_ok

   call read_process_grid(process_grid)
   call InitializeParallelism(process_grid(1), process_grid(2), process_grid(3), ierr)

   checks = 0
   failures = 0

   call test_initialized_metadata(process_grid, failures)
   call test_rank_mapping(failures)

   call check_partition('single-rank partition', 6, 3, 1, failures)
   call check_partition('uneven three-way partition', 9, 2, 3, failures)
   call check_partition('X-axis partition', 4, 1, NRPROCX, failures)
   call check_partition('Y-axis partition', 6, 2, NRPROCY, failures)
   call check_partition('Z-axis partition', 8, 3, NRPROCZ, failures)

   call check_dim_vectors('single-rank dimension vectors', 6, 1, 4, failures)
   call check_dim_vectors('uneven three-way vectors', 9, 3, 7, failures)
   call check_dim_vectors('X-axis dimension vectors', 4, NRPROCX, 2, failures)
   call check_dim_vectors('Y-axis dimension vectors', 6, NRPROCY, 3, failures)
   call check_dim_vectors('Z-axis dimension vectors', 8, NRPROCZ, 5, failures)

   finalized = .false.
   call Cleanup_Parallelism(ierr)
   call MPI_Finalized(finalized, finalize_ierr)
   cleanup_ok = ierr == MPI_SUCCESS .and. finalize_ierr == MPI_SUCCESS .and. &
      finalized
   checks = checks + 1
   if (.not. cleanup_ok) failures = failures + 1

   if (MYRANK == 0) then
      if (cleanup_ok) then
         write (*, '(A)') 'PASS Cleanup_Parallelism finalizes MPI'
      else
         write (*, '(A)') 'FAIL Cleanup_Parallelism finalizes MPI'
      end if

      if (failures == 0) then
         write (*, '(A,I0,A)') 'OK (', checks, ' parallelism checks)'
      else
         write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
            ' parallelism checks)'
      end if
   end if

   if (failures /= 0) stop 1

contains

   subroutine test_initialized_metadata(grid, failure_count)
      integer(kind=4), dimension(3), intent(in) :: grid
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4) :: world_rank, world_size, mpi_ierr
      character(len=7) :: expected_print_rank
      logical :: local_ok

      call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, mpi_ierr)
      local_ok = mpi_ierr == MPI_SUCCESS
      call MPI_Comm_size(MPI_COMM_WORLD, world_size, mpi_ierr)
      local_ok = local_ok .and. mpi_ierr == MPI_SUCCESS

      local_ok = local_ok .and. ierr == 0
      local_ok = local_ok .and. MYRANK == world_rank
      local_ok = local_ok .and. NRPROC == world_size
      local_ok = local_ok .and. NRPROC == product(grid)
      local_ok = local_ok .and. NRPROCX == grid(1)
      local_ok = local_ok .and. NRPROCY == grid(2)
      local_ok = local_ok .and. NRPROCZ == grid(3)
      local_ok = local_ok .and. MYRANKX >= 0 .and. MYRANKX < NRPROCX
      local_ok = local_ok .and. MYRANKY >= 0 .and. MYRANKY < NRPROCY
      local_ok = local_ok .and. MYRANKZ >= 0 .and. MYRANKZ < NRPROCZ
      local_ok = local_ok .and. LinearIndex(MYRANKX, MYRANKY, MYRANKZ) == MYRANK

      call assert_global('InitializeParallelism metadata', local_ok, failure_count)

      write (expected_print_rank, '(I5.5)') MYRANK
      call assert_global('PRINTRANK is a zero-padded rank prefix', &
         PRINTRANK == expected_print_rank, failure_count)
   end subroutine test_initialized_metadata


   subroutine test_rank_mapping(failure_count)
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4) :: rank, rankx, ranky, rankz
      integer(kind=4) :: expected_rankx, expected_ranky, expected_rankz
      logical :: local_ok

      local_ok = .true.
      do rank = 0, NRPROC - 1
         call Decompose(rank, rankx, ranky, rankz)
         expected_rankx = mod(rank, NRPROCX)
         expected_ranky = mod(rank/NRPROCX, NRPROCY)
         expected_rankz = rank/(NRPROCX*NRPROCY)
         local_ok = local_ok .and. rankx >= 0 .and. rankx < NRPROCX
         local_ok = local_ok .and. ranky >= 0 .and. ranky < NRPROCY
         local_ok = local_ok .and. rankz >= 0 .and. rankz < NRPROCZ
         local_ok = local_ok .and. rankx == expected_rankx
         local_ok = local_ok .and. ranky == expected_ranky
         local_ok = local_ok .and. rankz == expected_rankz
         local_ok = local_ok .and. LinearIndex(rankx, ranky, rankz) == rank
      end do

      call assert_global('rank mapping is X-fastest and exactly invertible', &
         local_ok, failure_count)
   end subroutine test_rank_mapping


   subroutine check_partition(label, n, p, nrproc, failure_count)
      character(len=*), intent(in) :: label
      integer(kind=4), intent(in) :: n, p, nrproc
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4), allocatable, dimension(:) :: ownership
      integer(kind=4) :: rank, nrcpp, ibeg, iend, mine, maxe
      integer(kind=4) :: expected_nrcpp, expected_ibeg, expected_iend
      integer(kind=4) :: expected_mine, expected_maxe, elems
      integer(kind=4) :: previous_end
      logical :: local_ok

      allocate (ownership(n + 1))
      ownership = 0
      elems = n + 1 - p
      expected_nrcpp = (n + nrproc)/nrproc
      previous_end = 0
      local_ok = n >= p .and. nrproc > 0

      do rank = 0, nrproc - 1
         call ComputeEndpoints(rank, nrproc, n, p, nrcpp, ibeg, iend, mine, maxe)

         expected_ibeg = expected_nrcpp*rank + 1
         expected_iend = min(expected_nrcpp*(rank + 1), n + 1)
         expected_mine = max(expected_ibeg - p - 1, 1)
         expected_maxe = min(expected_iend, elems)

         local_ok = local_ok .and. nrcpp == expected_nrcpp
         local_ok = local_ok .and. ibeg == expected_ibeg
         local_ok = local_ok .and. iend == expected_iend
         local_ok = local_ok .and. ibeg == previous_end + 1
         local_ok = local_ok .and. iend >= ibeg
         local_ok = local_ok .and. mine == expected_mine
         local_ok = local_ok .and. maxe == expected_maxe
         local_ok = local_ok .and. mine >= 1 .and. maxe <= elems
         local_ok = local_ok .and. mine <= maxe

         if (ibeg >= 1 .and. iend <= n + 1 .and. ibeg <= iend) then
            ownership(ibeg:iend) = ownership(ibeg:iend) + 1
         else
            local_ok = .false.
         end if
         previous_end = iend
      end do

      local_ok = local_ok .and. previous_end == n + 1
      local_ok = local_ok .and. all(ownership == 1)
      call assert_global(label, local_ok, failure_count)
      deallocate (ownership)
   end subroutine check_partition


   subroutine check_dim_vectors(label, n, nrproc, stride, failure_count)
      character(len=*), intent(in) :: label
      integer(kind=4), intent(in) :: n, nrproc, stride
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4), allocatable, dimension(:) :: dims, shifts
      integer(kind=4) :: i, nrcpp, expected_dim
      logical :: local_ok

      nrcpp = (n + nrproc)/nrproc
      call FillDimVector(dims, shifts, nrcpp, stride, n, nrproc)

      local_ok = allocated(dims) .and. allocated(shifts)
      local_ok = local_ok .and. size(dims) == nrproc
      local_ok = local_ok .and. size(shifts) == nrproc
      local_ok = local_ok .and. shifts(1) == 0
      local_ok = local_ok .and. sum(dims) == (n + 1)*stride
      local_ok = local_ok .and. all(dims > 0)
      local_ok = local_ok .and. all(mod(dims, stride) == 0)

      do i = 1, nrproc
         expected_dim = min(nrcpp, n + 1 - nrcpp*(i - 1))*stride
         local_ok = local_ok .and. dims(i) == expected_dim
         if (i > 1) then
            local_ok = local_ok .and. shifts(i) == shifts(i - 1) + dims(i - 1)
         end if
      end do

      local_ok = local_ok .and. shifts(nrproc) + dims(nrproc) == (n + 1)*stride
      call assert_global(label, local_ok, failure_count)
      deallocate (dims, shifts)
   end subroutine check_dim_vectors


   subroutine assert_global(label, local_ok, failure_count)
      character(len=*), intent(in) :: label
      logical, intent(in) :: local_ok
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4) :: global_flag, local_flag, reduce_ierr
      logical :: global_ok

      local_flag = merge(1, 0, local_ok)
      call MPI_Allreduce(local_flag, global_flag, 1, MPI_INTEGER, MPI_MIN, &
         MPI_COMM_WORLD, reduce_ierr)
      global_ok = reduce_ierr == MPI_SUCCESS .and. global_flag == 1

      checks = checks + 1
      if (.not. global_ok) failure_count = failure_count + 1

      if (MYRANK == 0) then
         if (global_ok) then
            write (*, '(A,A)') 'PASS ', trim(label)
         else
            write (*, '(A,A)') 'FAIL ', trim(label)
         end if
      end if
   end subroutine assert_global


   subroutine read_process_grid(grid)
      integer(kind=4), dimension(3), intent(out) :: grid
      character(len=32) :: argument
      integer(kind=4) :: i, read_status

      if (command_argument_count() /= 3) then
         write (*, '(A)') 'Usage: test_parallelism <procx> <procy> <procz>'
         stop 2
      end if

      do i = 1, 3
         call get_command_argument(i, argument)
         read (argument, *, iostat=read_status) grid(i)
         if (read_status /= 0 .or. grid(i) <= 0) then
            write (*, '(A,I0,A,A)') 'Invalid process-grid argument ', i, ': ', &
               trim(argument)
            stop 2
         end if
      end do
   end subroutine read_process_grid

end program test_parallelism
