program test_solution_output
   use Setup, only: ADS_Setup
   use parallelism, only: MYRANK
   use mpi, only: MPI_COMM_WORLD, real_bcast_calls, character_bcast_calls, &
                  bcast_contract_ok, broadcast_filename, reset_mpi_stub
   use my_mpi, only: gather_calls, captured_gather_root => captured_root, &
                     captured_gather_n => captured_n, &
                     captured_gather_p => captured_p, &
                     captured_gather_s => captured_s, captured_part, &
                     reset_my_mpi_stub
   use plot, only: captured_plot_root => captured_root, &
                   captured_plot_comm => captured_comm, captured_degree, &
                   captured_plot_n => captured_n, captured_nelem, &
                   captured_filename, captured_params, captured_ux, &
                   captured_uy, captured_uz, captured_coeffs, &
                   spline_plot_calls, reset_plot_stub
   use vtk, only: vtk_calls, vtk_filename, vtk_params, vtk_values, &
                  reset_vtk_stub
   use solution_output, only: PrintSolution
   implicit none

   integer(kind=4) :: checks, failures
   type(ADS_Setup) :: ads
   real(kind=8) :: part(2, 4)

   checks = 0
   failures = 0
   call prepare_inputs(ads, part)

   call test_root_orchestration(ads, part)
   call test_nonroot_broadcast_path(ads, part)

   if (failures == 0) then
      write (*, '(A,I0,A)') 'OK (', checks, ' solution_output checks)'
   else
      write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                 ' solution_output checks)'
      stop 1
   end if

contains

   subroutine prepare_inputs(space, local_part)
      type(ADS_Setup), intent(out) :: space
      real(kind=8), intent(out) :: local_part(2, 4)
      integer :: i

      space%n = (/1, 1, 1/)
      space%p = (/1, 1, 1/)
      space%s = (/2, 2, 2/)
      space%nelem = (/1, 1, 1/)
      allocate (space%Ux(0:3), space%Uy(0:3), space%Uz(0:3))
      do i = 0, 3
         space%Ux(i) = real(100 + i, kind=8)
         space%Uy(i) = real(200 + i, kind=8)
         space%Uz(i) = real(300 + i, kind=8)
      end do
      local_part = reshape((/1.d0, 2.d0, 3.d0, 4.d0, &
                             5.d0, 6.d0, 7.d0, 8.d0/), shape(local_part))
   end subroutine prepare_inputs


   subroutine test_root_orchestration(space, local_part)
      type(ADS_Setup), intent(in) :: space
      real(kind=8), intent(in) :: local_part(:, :)
      logical :: matches

      MYRANK = 0
      call reset_stubs()
      call PrintSolution(42, space, local_part)

      matches = gather_calls == 1 .and. captured_gather_root == 0
      matches = matches .and. all(captured_gather_n == space%n)
      matches = matches .and. all(captured_gather_p == space%p)
      matches = matches .and. all(captured_gather_s == space%s)
      matches = matches .and. allocated(captured_part)
      if (matches) matches = all(captured_part == local_part)
      call assert_true('PrintSolution gathers the exact local block on root', matches)

      matches = real_bcast_calls == 1 .and. character_bcast_calls == 1 .and. &
                bcast_contract_ok
      call assert_true('PrintSolution broadcasts coefficients and filename exactly once', &
                       matches)

      matches = spline_plot_calls == 1
      matches = matches .and. trim(captured_filename) == 'step42'
      matches = matches .and. captured_plot_root == 0
      matches = matches .and. captured_plot_comm == MPI_COMM_WORLD
      matches = matches .and. all(captured_degree == space%p)
      matches = matches .and. all(captured_plot_n == space%n)
      matches = matches .and. all(captured_nelem == space%nelem)
      matches = matches .and. params_are_default()
      matches = matches .and. all(captured_ux == space%Ux)
      matches = matches .and. all(captured_uy == space%Uy)
      matches = matches .and. all(captured_uz == space%Uz)
      call assert_true('PrintSolution forwards exact spline metadata and step filename', &
                       matches)

      call assert_true('PrintSolution forwards gathered coefficients in XYZ order', &
                       coefficients_match())

      matches = vtk_calls == 1 .and. trim(vtk_filename) == 'step42'
      matches = matches .and. allocated(vtk_values)
      if (matches) then
         matches = all(shape(vtk_values) == (/31, 31, 31/))
         matches = matches .and. vtk_values(1, 1, 1) == 10101.d0
         matches = matches .and. vtk_values(2, 1, 1) == 20101.d0
         matches = matches .and. vtk_values(1, 2, 1) == 10201.d0
         matches = matches .and. vtk_values(1, 1, 2) == 10102.d0
         matches = matches .and. vtk_values(31, 31, 31) == 313131.d0
      end if
      call assert_true('PrintSolution invokes VTK with exact sampled callback order', &
                       matches)
   end subroutine test_root_orchestration


   subroutine test_nonroot_broadcast_path(space, local_part)
      type(ADS_Setup), intent(in) :: space
      real(kind=8), intent(in) :: local_part(:, :)
      logical :: matches

      MYRANK = 1
      call reset_stubs()
      broadcast_filename = 'step7'
      call PrintSolution(7, space, local_part)

      matches = gather_calls == 1 .and. real_bcast_calls == 1 .and. &
                character_bcast_calls == 1 .and. bcast_contract_ok
      matches = matches .and. spline_plot_calls == 1
      matches = matches .and. trim(captured_filename) == 'step7'
      matches = matches .and. coefficients_match()
      call assert_true('nonroot receives exact coefficient and filename broadcasts', &
                       matches)

      call assert_true('nonroot never invokes the VTK callback', vtk_calls == 0)
      MYRANK = 0
   end subroutine test_nonroot_broadcast_path


   subroutine reset_stubs()
      call reset_mpi_stub()
      call reset_my_mpi_stub()
      call reset_plot_stub()
      call reset_vtk_stub()
   end subroutine reset_stubs


   function coefficients_match() result(matches)
      logical :: matches
      integer :: ix, iy, iz

      matches = allocated(captured_coeffs)
      if (.not. matches) return
      matches = all(shape(captured_coeffs) == (/2, 2, 2/))
      if (.not. matches) return
      do iz = 0, 1
         do iy = 0, 1
            do ix = 0, 1
               matches = matches .and. captured_coeffs(ix, iy, iz) == &
                         real(100*ix + 10*iy + iz, kind=8)
            end do
         end do
      end do
   end function coefficients_match


   function params_are_default() result(matches)
      logical :: matches

      matches = captured_params%startx == 0.d0 .and. &
                captured_params%endx == 1.d0 .and. &
                captured_params%starty == 0.d0 .and. &
                captured_params%endy == 1.d0 .and. &
                captured_params%startz == 0.d0 .and. &
                captured_params%endz == 1.d0 .and. &
                captured_params%resx == 31 .and. &
                captured_params%resy == 31 .and. &
                captured_params%resz == 31
   end function params_are_default


   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      checks = checks + 1
      if (condition) then
         write (*, '(A)') 'PASS '//trim(label)
      else
         failures = failures + 1
         write (*, '(A)') 'FAIL '//trim(label)
      end if
   end subroutine assert_true

end program test_solution_output
