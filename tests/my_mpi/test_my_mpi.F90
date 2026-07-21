program test_my_mpi
   use communicators, only: COMMX, COMMY, COMMZ, CreateCommunicators, &
      Cleanup_Communicators, processors
   use mpi
   use my_mpi, only: AllGather, Delinearize, DistributeSpline, Gather, &
      GatherFullSolution, Linearize, Scatter
   use parallelism, only: MYRANK, MYRANKX, MYRANKY, MYRANKZ, NRPROC, &
      NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints, FillDimVector, &
      InitializeParallelism, Cleanup_Parallelism
   implicit none

   integer(kind=4), dimension(3) :: process_grid
   integer(kind=4) :: checks, failures, ierr

   call read_process_grid(process_grid)
   call InitializeParallelism(process_grid(1), process_grid(2), process_grid(3), ierr)
   call CreateCommunicators(ierr)

   checks = 0
   failures = 0
   call test_linear_storage(failures)
   call test_directional_collectives(1, failures)
   call test_directional_collectives(2, failures)
   call test_directional_collectives(3, failures)
   call test_full_solution_gather(0, failures)
   if (NRPROC > 1) call test_full_solution_gather(NRPROC - 1, failures)
   call test_spline_distribution(failures)

   if (MYRANK == 0) then
      if (failures == 0) then
         write(*, '(A,I0,A)') 'OK (', checks, ' my_mpi checks)'
      else
         write(*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, ' my_mpi checks)'
      end if
   end if

   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

   if (failures /= 0) stop 1

contains

   subroutine test_linear_storage(failure_count)
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8), dimension(4, 3) :: original, reconstructed
      real(kind=8), dimension(12) :: linear
      integer(kind=4) :: i, j

      do i = 1, size(original, 1)
         do j = 1, size(original, 2)
            original(i, j) = storage_value(MYRANK, i, j)
         end do
      end do

      call Linearize(original, linear, size(original, 1), size(original, 2))
      call assert_global('Linearize uses row-major segments', &
         all(linear == reshape(transpose(original), shape(linear))), failure_count)

      reconstructed = -huge(0.0d0)
      call Delinearize(linear, reconstructed, size(original, 1), size(original, 2))
      call assert_global('Delinearize reverses Linearize', &
         all(reconstructed == original), failure_count)
   end subroutine test_linear_storage


   subroutine test_directional_collectives(axis, failure_count)
      integer(kind=4), intent(in) :: axis
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4), parameter :: n = 4
      integer(kind=4), parameter :: p = 1
      integer(kind=4), parameter :: stride = 4
      integer(kind=4), allocatable :: dims(:), shifts(:)
      real(kind=8), allocatable :: local(:,:), gathered(:,:), scattered(:,:), &
         all_gathered(:,:), global_values(:,:)
      integer(kind=4) :: axis_rank, axis_size, comm
      integer(kind=4) :: nrcpp, ibeg, iend, mine, maxe, owned
      integer(kind=4) :: local_row, global_row, lane
      logical :: local_ok
      character(len=96) :: label

      call select_axis(axis, axis_rank, axis_size, comm)
      call ComputeEndpoints(axis_rank, axis_size, n, p, nrcpp, ibeg, iend, mine, maxe)
      owned = iend - ibeg + 1
      call FillDimVector(dims, shifts, nrcpp, stride, n, axis_size)

      allocate(local(owned, stride))
      allocate(gathered(n + 1, stride))
      allocate(scattered(owned, stride))
      allocate(all_gathered(n + 1, stride))
      allocate(global_values(n + 1, stride))

      do global_row = 1, n + 1
         do lane = 1, stride
            global_values(global_row, lane) = directional_value(axis, global_row, lane)
         end do
      end do

      do local_row = 1, owned
         global_row = ibeg + local_row - 1
         local(local_row, :) = global_values(global_row, :)
      end do

      gathered = -huge(0.0d0)
      call Gather(local, gathered, n, owned, stride, dims, shifts, comm, ierr)
      local_ok = .true.
      if (axis_rank == 0) local_ok = all(gathered == global_values)
      write(label, '(A,I0)') 'Gather exact on axis ', axis
      call assert_global(trim(label), local_ok, failure_count)

      scattered = -huge(0.0d0)
      call Scatter(global_values, scattered, n, owned, stride, dims, shifts, comm, ierr)
      local_ok = all(scattered == local)
      write(label, '(A,I0)') 'Scatter exact on axis ', axis
      call assert_global(trim(label), local_ok, failure_count)

      all_gathered = -huge(0.0d0)
      call AllGather(local, all_gathered, n, owned, stride, dims, shifts, comm)
      local_ok = all(all_gathered == global_values)
      write(label, '(A,I0)') 'AllGather exact on axis ', axis
      call assert_global(trim(label), local_ok, failure_count)

      deallocate(dims, shifts)
      deallocate(local, gathered, scattered, all_gathered, global_values)
   end subroutine test_directional_collectives


   subroutine test_full_solution_gather(root, failure_count)
      integer(kind=4), intent(in) :: root
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4), dimension(3), parameter :: n = (/4, 5, 6/)
      integer(kind=4), dimension(3), parameter :: p = (/1, 2, 3/)
      integer(kind=4), dimension(3) :: begs, ends, owned
      real(kind=8), allocatable :: part(:,:), full(:,:,:)
      integer(kind=4) :: nrcpp, mine, maxe
      integer(kind=4) :: lx, ly, lz, local_yz, gx, gy, gz
      logical :: local_ok
      character(len=96) :: label

      call ComputeEndpoints(MYRANKX, NRPROCX, n(1), p(1), nrcpp, &
         begs(1), ends(1), mine, maxe)
      call ComputeEndpoints(MYRANKY, NRPROCY, n(2), p(2), nrcpp, &
         begs(2), ends(2), mine, maxe)
      call ComputeEndpoints(MYRANKZ, NRPROCZ, n(3), p(3), nrcpp, &
         begs(3), ends(3), mine, maxe)
      owned = ends - begs + 1

      allocate(part(owned(1), owned(2)*owned(3)))
      do lz = 1, owned(3)
         gz = begs(3) + lz - 2
         do ly = 1, owned(2)
            gy = begs(2) + ly - 2
            local_yz = ly + (lz - 1)*owned(2)
            do lx = 1, owned(1)
               gx = begs(1) + lx - 2
               part(lx, local_yz) = tensor_value(gx, gy, gz)
            end do
         end do
      end do

      call GatherFullSolution(root, part, full, n, p, owned)
      local_ok = .true.
      if (MYRANK == root) then
         local_ok = allocated(full)
         if (local_ok) then
            local_ok = all(shape(full) == n + 1)
            do gz = 0, n(3)
               do gy = 0, n(2)
                  do gx = 0, n(1)
                     local_ok = local_ok .and. full(gx, gy, gz) == tensor_value(gx, gy, gz)
                  end do
               end do
            end do
         end if
      else
         local_ok = .not. allocated(full)
      end if

      write(label, '(A,I0)') 'GatherFullSolution exact at root ', root
      call assert_global(trim(label), local_ok, failure_count)
      deallocate(part)
      if (allocated(full)) deallocate(full)
   end subroutine test_full_solution_gather


   subroutine test_spline_distribution(failure_count)
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4), dimension(3), parameter :: nrcpp = (/2, 2, 2/)
      real(kind=8), parameter :: sentinel = -987654321.25d0
      real(kind=8), allocatable :: spline(:,:,:,:)
      integer(kind=4) :: block_size, i, ix, iy, iz, nx, ny, nz, source_rank
      logical :: local_ok, neighbour_exists

      block_size = product(nrcpp)
      allocate(spline(block_size, 3, 3, 3))
      spline = sentinel
      do i = 1, block_size
         spline(i, 2, 2, 2) = neighbour_value(MYRANK, i)
      end do

      call DistributeSpline(spline, nrcpp)

      local_ok = .true.
      do iz = 1, 3
         nz = MYRANKZ + iz - 2
         do iy = 1, 3
            ny = MYRANKY + iy - 2
            do ix = 1, 3
               nx = MYRANKX + ix - 2
               neighbour_exists = nx >= 0 .and. nx < NRPROCX .and. &
                  ny >= 0 .and. ny < NRPROCY .and. nz >= 0 .and. nz < NRPROCZ
               if (neighbour_exists) then
                  source_rank = processors(nx + 1, ny + 1, nz + 1)
                  do i = 1, block_size
                     local_ok = local_ok .and. &
                        spline(i, ix, iy, iz) == neighbour_value(source_rank, i)
                  end do
               else
                  local_ok = local_ok .and. all(spline(:, ix, iy, iz) == sentinel)
               end if
            end do
         end do
      end do

      call assert_global('DistributeSpline fills valid neighbour blocks', &
         local_ok, failure_count)
      deallocate(spline)
   end subroutine test_spline_distribution


   subroutine select_axis(axis, axis_rank, axis_size, comm)
      integer(kind=4), intent(in) :: axis
      integer(kind=4), intent(out) :: axis_rank, axis_size, comm

      select case (axis)
      case (1)
         axis_rank = MYRANKX
         axis_size = NRPROCX
         comm = COMMX
      case (2)
         axis_rank = MYRANKY
         axis_size = NRPROCY
         comm = COMMY
      case (3)
         axis_rank = MYRANKZ
         axis_size = NRPROCZ
         comm = COMMZ
      case default
         write(*, '(A,I0)') 'Invalid test axis: ', axis
         stop 2
      end select
   end subroutine select_axis


   subroutine assert_global(label, local_ok, failure_count)
      character(len=*), intent(in) :: label
      logical, intent(in) :: local_ok
      integer(kind=4), intent(inout) :: failure_count
      logical :: global_ok
      integer(kind=4) :: global_flag, local_flag, reduce_ierr

      local_flag = merge(1, 0, local_ok)
      call MPI_Allreduce(local_flag, global_flag, 1, MPI_INTEGER, MPI_MIN, &
         MPI_COMM_WORLD, reduce_ierr)
      global_ok = global_flag == 1
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


   pure real(kind=8) function storage_value(rank, row, column) result(value)
      integer(kind=4), intent(in) :: rank, row, column

      value = 10000.0d0*rank + 101.0d0*row - 7.0d0*column + 0.125d0
   end function storage_value


   real(kind=8) function directional_value(axis, row, lane) result(value)
      integer(kind=4), intent(in) :: axis, row, lane
      integer(kind=4) :: fixed_one, fixed_two

      select case (axis)
      case (1)
         fixed_one = MYRANKY
         fixed_two = MYRANKZ
      case (2)
         fixed_one = MYRANKX
         fixed_two = MYRANKZ
      case (3)
         fixed_one = MYRANKX
         fixed_two = MYRANKY
      end select

      value = 1000000.0d0*axis + 100000.0d0*fixed_two + &
         10000.0d0*fixed_one + 101.0d0*row - 13.0d0*lane + 0.375d0
   end function directional_value


   pure real(kind=8) function tensor_value(x, y, z) result(value)
      integer(kind=4), intent(in) :: x, y, z

      value = 10000.0d0*z + 100.0d0*y + 3.0d0*x + 0.625d0
   end function tensor_value


   pure real(kind=8) function neighbour_value(rank, index) result(value)
      integer(kind=4), intent(in) :: rank, index

      value = 1000.0d0*rank + 11.0d0*index + 0.75d0
   end function neighbour_value


   subroutine read_process_grid(grid)
      integer(kind=4), intent(out), dimension(3) :: grid
      character(len=32) :: argument
      integer(kind=4) :: i, read_status

      if (command_argument_count() /= 3) then
         write(*, '(A)') 'Usage: test_my_mpi <procx> <procy> <procz>'
         stop 2
      end if

      do i = 1, 3
         call get_command_argument(i, argument)
         read(argument, *, iostat=read_status) grid(i)
         if (read_status /= 0 .or. grid(i) <= 0) then
            write(*, '(A,I0,A,A)') 'Invalid process-grid argument ', i, ': ', trim(argument)
            stop 2
         end if
      end do
   end subroutine read_process_grid

end program test_my_mpi
