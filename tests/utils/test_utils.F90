module norm_l2_test_support
   use Setup, only: ADS_Setup
   implicit none

contains

   pure real(kind=8) function affine_field(X) result(value)
      real(kind=8), intent(in), dimension(3) :: X

      value = 1.0d0 + 2.0d0*X(1) - 3.0d0*X(2) + 0.5d0*X(3)
   end function affine_field


   pure real(kind=8) function affine_field_plus_one(X) result(value)
      real(kind=8), intent(in), dimension(3) :: X

      value = affine_field(X) + 1.0d0
   end function affine_field_plus_one


   pure real(kind=8) function zero_field(X) result(value)
      real(kind=8), intent(in), dimension(3) :: X

      value = 0.0d0*sum(X)
   end function zero_field


   pure real(kind=8) function greville(U, p, basis_index) result(value)
      real(kind=8), intent(in), dimension(:) :: U
      integer(kind=4), intent(in) :: p, basis_index
      integer(kind=4) :: first

      first = lbound(U, 1) + basis_index + 1
      value = sum(U(first:first + p - 1))/real(p, kind=8)
   end function greville


   subroutine fill_affine_coefficients(ads, part)
      type(ADS_Setup), intent(in) :: ads
      real(kind=8), intent(out), dimension(:, :) :: part
      real(kind=8), dimension(3) :: X
      integer(kind=4) :: global_x, global_y, global_z
      integer(kind=4) :: local_x, local_y, local_z, local_yz

      do local_z = 1, ads%s(3)
         global_z = ads%ibeg(3) + local_z - 2
         X(3) = greville(ads%Uz, ads%p(3), global_z)
         do local_y = 1, ads%s(2)
            global_y = ads%ibeg(2) + local_y - 2
            X(2) = greville(ads%Uy, ads%p(2), global_y)
            local_yz = local_y + (local_z - 1)*ads%s(2)
            do local_x = 1, ads%s(1)
               global_x = ads%ibeg(1) + local_x - 2
               X(1) = greville(ads%Ux, ads%p(1), global_x)
               part(local_x, local_yz) = affine_field(X)
            end do
         end do
      end do
   end subroutine fill_affine_coefficients


   subroutine remap_geometry(ads)
      use basis, only: BasisData
      type(ADS_Setup), intent(inout) :: ads
      real(kind=8) :: t
      integer(kind=4) :: i

      do i = lbound(ads%Ux, 1), ubound(ads%Ux, 1)
         t = ads%Ux(i)
         ads%Ux(i) = -1.0d0 + 2.0d0*t*t
      end do
      do i = lbound(ads%Uy, 1), ubound(ads%Uy, 1)
         t = ads%Uy(i)
         ads%Uy(i) = 2.0d0*(0.25d0*t + 0.75d0*t*t)
      end do
      do i = lbound(ads%Uz, 1), ubound(ads%Uz, 1)
         t = ads%Uz(i)
         ads%Uz(i) = 1.0d0 + 3.0d0*(0.5d0*t + 0.5d0*t*t*t)
      end do

      call BasisData(ads%p(1), ads%m(1), ads%Ux, 1, ads%ng(1), &
         ads%nelem(1), ads%Ox, ads%Jx, ads%Wx, ads%Xx, ads%NNx)
      call BasisData(ads%p(2), ads%m(2), ads%Uy, 1, ads%ng(2), &
         ads%nelem(2), ads%Oy, ads%Jy, ads%Wy, ads%Xy, ads%NNy)
      call BasisData(ads%p(3), ads%m(3), ads%Uz, 1, ads%ng(3), &
         ads%nelem(3), ads%Oz, ads%Jz, ads%Wz, ads%Xz, ads%NNz)
   end subroutine remap_geometry


   subroutine check_norm(label, actual, expected, tolerance, failures)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use mpi
      use parallelism, only: MYRANK
      character(len=*), intent(in) :: label
      real(kind=8), intent(in) :: actual, expected, tolerance
      integer(kind=4), intent(inout) :: failures
      real(kind=8) :: global_min, global_max
      logical :: value_ok, ranks_agree
      integer(kind=4) :: ierr

      call MPI_Allreduce(actual, global_min, 1, MPI_DOUBLE_PRECISION, &
         MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(actual, global_max, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, MPI_COMM_WORLD, ierr)

      value_ok = ieee_is_finite(actual) .and. abs(actual - expected) <= tolerance
      ranks_agree = ieee_is_finite(global_min) .and. ieee_is_finite(global_max) .and. &
         abs(global_max - global_min) <= tolerance

      if (.not. value_ok) failures = failures + 1
      if (.not. ranks_agree) failures = failures + 1

      if (MYRANK == 0) then
         if (value_ok .and. ranks_agree) then
            write(*, '(A,A,A,ES14.6)') 'PASS ', trim(label), ': ', actual
         else
            write(*, '(A,A)') 'FAIL ', trim(label)
            write(*, '(A,3(1X,ES14.6))') '  actual/expected/tolerance:', &
               actual, expected, tolerance
            write(*, '(A,2(1X,ES14.6))') '  rank min/max:', global_min, global_max
         end if
      end if
   end subroutine check_norm

end module norm_l2_test_support


program test_utils
   use Setup, only: ADS_Setup
   use ads_lifecycle, only: initialize_setup, Cleanup_ADS
   use communicators, only: CreateCommunicators, Cleanup_Communicators
   use norm_l2_test_support
   use parallelism, only: MYRANK, NRPROC, InitializeParallelism, Cleanup_Parallelism
   use projection_engine, only: global2local_call_count, reset_global2local_spy
   use utils, only: NormL2, Norm_L2, int2str
   use mpi
   implicit none

   type(ADS_Setup) :: ads
   real(kind=8), allocatable, dimension(:, :) :: part
   real(kind=8) :: norm
   integer(kind=4), dimension(3) :: n, p, continuity, ng, process_grid
   integer(kind=4) :: failures, global_failures, ierr

   call read_process_grid(process_grid)
   call InitializeParallelism(process_grid(1), process_grid(2), process_grid(3), ierr)
   call CreateCommunicators(ierr)

   n = (/4, 6, 8/)
   p = (/1, 2, 3/)
   continuity = p - 1
   ng = p + 1
   call initialize_setup(n, p, continuity, ng, ads, ierr)
   call remap_geometry(ads)

   allocate(part(ads%s(1), ads%s(2)*ads%s(3)))
   failures = 0
   call check_int2str(failures)
   call check_legacy_norm_l2(failures)
   call reset_global2local_spy()

   call fill_affine_coefficients(ads, part)
   call NormL2(ads, part, norm)
   call check_norm('solution affine', norm, sqrt(61.0d0), 1.0d-11, failures)

   call NormL2(ads, part, norm, affine_field)
   call check_norm('matched affine error', norm, 0.0d0, 1.0d-12, failures)

   call NormL2(ads, part, norm, affine_field_plus_one)
   call check_norm('shifted affine error', norm, sqrt(12.0d0), 1.0d-11, failures)

   call NormL2(ads, part, norm, zero_field)
   call check_norm('zero-reference affine error', norm, sqrt(61.0d0), 1.0d-11, failures)

   part = 2.0d0
   call NormL2(ads, part, norm)
   call check_norm('constant solution', norm, sqrt(48.0d0), 1.0d-11, failures)

   call NormL2(ads, part, norm, affine_field)
   call check_norm('constant-versus-affine error', norm, sqrt(145.0d0), 1.0d-11, failures)

   part = 0.0d0
   call NormL2(ads, part, norm)
   call check_norm('zero solution', norm, 0.0d0, 1.0d-12, failures)

   call NormL2(ads, part, norm, zero_field)
   call check_norm('zero error', norm, 0.0d0, 1.0d-12, failures)

   if (global2local_call_count /= 0) failures = failures + 1

   call MPI_Allreduce(failures, global_failures, 1, MPI_INTEGER, &
      MPI_SUM, MPI_COMM_WORLD, ierr)

   if (MYRANK == 0) then
      if (global_failures == 0) then
         write(*, '(A,I0,A,I0,A)') 'OK (', 8, &
            ' NormL2 checks plus legacy assembly and int2str, ', NRPROC, ' MPI ranks)'
      else
         write(*, '(A,I0,A)') 'FAILED (', global_failures, ' accumulated failures)'
      end if
   end if

   deallocate(part)
   call Cleanup_ADS(ads, ierr)
   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

   if (global_failures /= 0) stop 1

contains

   subroutine check_legacy_norm_l2(local_failures)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      integer(kind=4), intent(inout) :: local_failures
      integer(kind=4), parameter :: legacy_n = 2, legacy_p = 1
      integer(kind=4), dimension(3) :: legacy_shape, legacy_degree
      integer(kind=4), dimension(3) :: legacy_elements, legacy_begin, legacy_end
      integer(kind=4), dimension(3) :: legacy_rank, legacy_processes
      real(kind=8), dimension(0:legacy_n + legacy_p + 1) :: knots
      real(kind=8), dimension(0:legacy_n, 0:(legacy_n + 1)*(legacy_n + 1) - 1) :: assembled
      real(kind=8), dimension(0:0, 0:0) :: centre
      real(kind=8), dimension(0:legacy_n) :: one_dimensional_integral
      real(kind=8) :: expected
      integer(kind=4) :: x, y, z, yz

      call reset_global2local_spy()
      legacy_shape = legacy_n
      legacy_degree = legacy_p
      legacy_elements = 2
      legacy_begin = 1
      legacy_end = legacy_n + 1
      legacy_rank = 0
      legacy_processes = 1
      knots = (/0.0d0, 0.0d0, 0.5d0, 1.0d0, 1.0d0/)
      one_dimensional_integral = (/0.25d0, 0.5d0, 0.25d0/)

      call Norm_L2(knots, knots, knots, legacy_degree, legacy_shape, &
         legacy_elements, legacy_begin, legacy_end, legacy_rank, &
         legacy_processes, assembled)

      do z = 0, legacy_n
         do y = 0, legacy_n
            yz = y + z*(legacy_n + 1)
            do x = 0, legacy_n
               expected = one_dimensional_integral(x)*one_dimensional_integral(y)* &
                  one_dimensional_integral(z)
               if (.not. (ieee_is_finite(assembled(x, yz)) .and. &
                          abs(assembled(x, yz) - expected) <= 1.0d-13)) then
                  local_failures = local_failures + 1
               end if
            end do
         end do
      end do

      legacy_begin = 2
      legacy_end = 2
      call Norm_L2(knots, knots, knots, legacy_degree, legacy_shape, &
         legacy_elements, legacy_begin, legacy_end, legacy_rank, &
         legacy_processes, centre)
      if (.not. (ieee_is_finite(centre(0, 0)) .and. &
                 abs(centre(0, 0) - 0.125d0) <= 1.0d-13)) then
         local_failures = local_failures + 1
      end if
      if (global2local_call_count == 0) local_failures = local_failures + 1

      call check_legacy_norm_l2_uneven_partition(local_failures)
   end subroutine check_legacy_norm_l2


   subroutine check_legacy_norm_l2_uneven_partition(local_failures)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      integer(kind=4), intent(inout) :: local_failures
      integer(kind=4), parameter :: uneven_n = 8, transverse_n = 2, degree = 1
      integer(kind=4), parameter :: virtual_ranks = 3
      integer(kind=4), dimension(3) :: shape, degrees, elements
      integer(kind=4), dimension(3) :: local_begin, local_end
      integer(kind=4), dimension(3) :: virtual_rank, process_grid
      real(kind=8), dimension(0:uneven_n + degree + 1) :: knots_x
      real(kind=8), dimension(0:transverse_n + degree + 1) :: knots_yz
      real(kind=8), dimension(0:2, 0:8) :: local_assembly
      real(kind=8), dimension(0:uneven_n) :: integral_x
      real(kind=8), dimension(0:transverse_n) :: integral_yz
      real(kind=8) :: expected
      integer(kind=4) :: global_x, local_x, rank_x, y, yz, z

      shape = (/uneven_n, transverse_n, transverse_n/)
      degrees = degree
      elements = (/uneven_n, transverse_n, transverse_n/)
      process_grid = (/virtual_ranks, 1, 1/)
      local_begin = 1
      local_end = (/3, transverse_n + 1, transverse_n + 1/)
      virtual_rank = 0
      knots_x = (/0.0d0, 0.0d0, 0.125d0, 0.25d0, 0.375d0, &
         0.5d0, 0.625d0, 0.75d0, 0.875d0, 1.0d0, 1.0d0/)
      knots_yz = (/0.0d0, 0.0d0, 0.5d0, 1.0d0, 1.0d0/)
      integral_x = 0.125d0
      integral_x(0) = 0.0625d0
      integral_x(uneven_n) = 0.0625d0
      integral_yz = (/0.25d0, 0.5d0, 0.25d0/)

      do rank_x = 0, virtual_ranks - 1
         virtual_rank(1) = rank_x
         local_begin(1) = 3*rank_x + 1
         local_end(1) = 3*rank_x + 3
         call Norm_L2(knots_x, knots_yz, knots_yz, degrees, shape, &
            elements, local_begin, local_end, virtual_rank, process_grid, &
            local_assembly)

         do z = 0, transverse_n
            do y = 0, transverse_n
               yz = y + z*(transverse_n + 1)
               do local_x = 0, 2
                  global_x = 3*rank_x + local_x
                  expected = integral_x(global_x)*integral_yz(y)*integral_yz(z)
                  if (.not. (ieee_is_finite(local_assembly(local_x, yz)) .and. &
                             abs(local_assembly(local_x, yz) - expected) <= &
                             1.0d-13)) then
                     local_failures = local_failures + 1
                  end if
               end do
            end do
         end do
      end do
   end subroutine check_legacy_norm_l2_uneven_partition

   subroutine check_int2str(local_failures)
      integer(kind=4), intent(inout) :: local_failures
      integer(kind=4), parameter :: max_i4 = huge(0_4)
      integer(kind=4) :: min_i4
      character(len=4) :: zero
      character(len=6) :: positive, negative
      character(len=10) :: maximum
      character(len=11) :: minimum

      min_i4 = -huge(0_4)
      min_i4 = min_i4 - 1_4

      call int2str(0_4, zero)
      call int2str(42_4, positive)
      call int2str(-42_4, negative)
      call int2str(max_i4, maximum)
      call int2str(min_i4, minimum)

      if (zero /= '0   ') local_failures = local_failures + 1
      if (positive /= '42    ') local_failures = local_failures + 1
      if (negative /= '-42   ') local_failures = local_failures + 1
      if (maximum /= '2147483647') local_failures = local_failures + 1
      if (minimum /= '-2147483648') local_failures = local_failures + 1
   end subroutine check_int2str

   subroutine read_process_grid(grid)
      integer(kind=4), intent(out), dimension(3) :: grid
      character(len=32) :: argument
      integer(kind=4) :: i, read_status

      if (command_argument_count() /= 3) then
         write(*, '(A)') 'Usage: test_norm_l2_mpi <procx> <procy> <procz>'
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

end program test_utils
