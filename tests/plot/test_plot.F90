program test_plot
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mpi
   use plot, only: PlotParams, SavePlot, SaveSplinePlot, SaveSplinePlotMPI
   use plot_test_callbacks, only: callback_calls, captured_filename, &
      captured_params, captured_values, scalar_function, capture_output, &
      reset_capture
   implicit none

   integer(kind=4) :: checks, failures, global_failures
   integer(kind=4) :: ierr, myrank, nranks
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

   checks = 0
   failures = 0
   callback_calls = 0

   call test_save_plot()
   call test_save_plot_singleton()
   call test_save_spline_plot()
   call test_save_spline_plot_singleton()
   call test_save_spline_plot_mpi()
   call test_mpi_singleton_resolution()

   call MPI_Allreduce(failures, global_failures, 1, MPI_INTEGER, MPI_SUM, &
                      MPI_COMM_WORLD, ierr)
   if (myrank == 0) then
      if (global_failures == 0) then
         write (*, '(A,I0,A,I0,A)') 'OK (', checks, ' plot checks on ', &
                                    nranks, ' rank(s))'
      else
         write (*, '(A,I0,A)') 'FAILED (global failure count ', &
                                global_failures, ')'
      end if
   end if

   call MPI_Finalize(ierr)
   if (global_failures /= 0) stop 1

contains

   subroutine test_save_plot()
      type(PlotParams) :: params
      integer :: ix, iy, iz
      real(kind=8) :: x, y, z
      logical :: matches

      params = PlotParams(-1.d0, 1.d0, 2.d0, 4.d0, 10.d0, 14.d0, 3, 2, 3)
      call reset_capture()
      call SavePlot('scalar_callback', scalar_function, capture_output, params)

      matches = callback_calls == 1
      matches = matches .and. trim(captured_filename) == 'scalar_callback'
      matches = matches .and. params_match(captured_params, params)
      matches = matches .and. allocated(captured_values)
      if (matches) matches = all(shape(captured_values) == (/3, 2, 3/))
      if (matches) then
         do iz = 1, params%resz
            z = coordinate(iz, params%resz, params%startz, params%endz)
            do iy = 1, params%resy
               y = coordinate(iy, params%resy, params%starty, params%endy)
               do ix = 1, params%resx
                  x = coordinate(ix, params%resx, params%startx, params%endx)
                  matches = matches .and. &
                     captured_values(ix, iy, iz) == scalar_function(x, y, z)
               end do
            end do
         end do
      end if
      call assert_true('SavePlot forwards exact filename, grid, and callback data', &
                       matches)
   end subroutine test_save_plot


   subroutine test_save_plot_singleton()
      type(PlotParams) :: params
      real(kind=8) :: expected
      logical :: matches

      params = PlotParams(2.d0, 9.d0, -3.d0, 7.d0, 4.d0, 11.d0, 1, 1, 1)
      call reset_capture()
      call SavePlot('scalar_singleton', scalar_function, capture_output, params)

      expected = scalar_function(params%startx, params%starty, params%startz)
      matches = callback_calls == 1
      matches = matches .and. trim(captured_filename) == 'scalar_singleton'
      matches = matches .and. allocated(captured_values)
      if (matches) then
         matches = all(shape(captured_values) == (/1, 1, 1/))
         matches = matches .and. ieee_is_finite(captured_values(1, 1, 1))
         matches = matches .and. captured_values(1, 1, 1) == expected
      end if
      call assert_true('SavePlot uses lower bounds for singleton resolution', matches)
   end subroutine test_save_plot_singleton


   subroutine test_save_spline_plot()
      type(PlotParams) :: params
      real(kind=8) :: ux(0:1), uy(0:1), uz(0:1), coeffs(0:0, 0:0, 0:0)
      integer :: ix, iy, iz
      real(kind=8) :: x, y, z, expected
      logical :: matches

      ux = (/0.d0, 1.d0/)
      uy = (/0.d0, 1.d0/)
      uz = (/0.d0, 1.d0/)
      coeffs = 7.d0
      params = PlotParams(0.d0, 2.d0, -1.d0, 1.d0, 3.d0, 5.d0, 2, 2, 2)
      call reset_capture()
      call SaveSplinePlot('serial_spline', &
                          ux, 0, 0, 1, uy, 0, 0, 1, uz, 0, 0, 1, &
                          coeffs, capture_output, params)

      matches = callback_calls == 1
      matches = matches .and. trim(captured_filename) == 'serial_spline'
      matches = matches .and. params_match(captured_params, params)
      matches = matches .and. allocated(captured_values)
      if (matches) then
         do iz = 1, params%resz
            z = coordinate(iz, params%resz, params%startz, params%endz)
            do iy = 1, params%resy
               y = coordinate(iy, params%resy, params%starty, params%endy)
               do ix = 1, params%resx
                  x = coordinate(ix, params%resx, params%startx, params%endx)
                  expected = coeffs(0, 0, 0) + 100.d0*x + 10.d0*y + z
                  matches = matches .and. captured_values(ix, iy, iz) == expected
               end do
            end do
         end do
      end if
      call assert_true('SaveSplinePlot samples the exact tensor-grid callback data', &
                       matches)
   end subroutine test_save_spline_plot


   subroutine test_save_spline_plot_singleton()
      type(PlotParams) :: params
      real(kind=8) :: ux(0:1), uy(0:1), uz(0:1), coeffs(0:0, 0:0, 0:0)
      real(kind=8) :: expected
      logical :: matches

      ux = (/0.d0, 1.d0/)
      uy = (/0.d0, 1.d0/)
      uz = (/0.d0, 1.d0/)
      coeffs = 5.d0
      params = PlotParams(2.d0, 9.d0, -3.d0, 7.d0, 4.d0, 11.d0, 1, 1, 1)
      call reset_capture()
      call SaveSplinePlot('spline_singleton', &
                          ux, 0, 0, 1, uy, 0, 0, 1, uz, 0, 0, 1, &
                          coeffs, capture_output, params)

      expected = coeffs(0, 0, 0) + 100.d0*params%startx + &
                 10.d0*params%starty + params%startz
      matches = callback_calls == 1
      matches = matches .and. trim(captured_filename) == 'spline_singleton'
      matches = matches .and. allocated(captured_values)
      if (matches) then
         matches = all(shape(captured_values) == (/1, 1, 1/))
         matches = matches .and. ieee_is_finite(captured_values(1, 1, 1))
         matches = matches .and. captured_values(1, 1, 1) == expected
      end if
      call assert_true('SaveSplinePlot uses lower bounds for singleton resolution', &
                       matches)
   end subroutine test_save_spline_plot_singleton


   subroutine test_save_spline_plot_mpi()
      type(PlotParams) :: params
      real(kind=8) :: ux(0:1), uy(0:1), uz(0:1), coeffs(0:0, 0:0, 0:0)
      integer :: ix, iy, iz, root
      real(kind=8) :: x, y, z, expected
      logical :: matches

      ux = (/0.d0, 1.d0/)
      uy = (/0.d0, 1.d0/)
      uz = (/0.d0, 1.d0/)
      coeffs = -2.d0
      params = PlotParams(-2.d0, 2.d0, 1.d0, 3.d0, 4.d0, 8.d0, 3, 2, 2)
      root = nranks - 1
      call reset_capture()
      call SaveSplinePlotMPI('mpi_spline', &
                             ux, 0, 0, 1, uy, 0, 0, 1, uz, 0, 0, 1, &
                             coeffs, capture_output, params, root, MPI_COMM_WORLD)

      matches = callback_calls == merge(1, 0, myrank == root)
      if (myrank == root) then
         matches = matches .and. trim(captured_filename) == 'mpi_spline'
         matches = matches .and. params_match(captured_params, params)
         matches = matches .and. allocated(captured_values)
         if (matches) then
            do iz = 1, params%resz
               z = coordinate(iz, params%resz, params%startz, params%endz)
               do iy = 1, params%resy
                  y = coordinate(iy, params%resy, params%starty, params%endy)
                  do ix = 1, params%resx
                     x = coordinate(ix, params%resx, params%startx, params%endx)
                     expected = coeffs(0, 0, 0) + 100.d0*x + 10.d0*y + z
                     matches = matches .and. &
                               captured_values(ix, iy, iz) == expected
                  end do
               end do
            end do
         end if
      end if
      call assert_true('SaveSplinePlotMPI reconstructs exact global data order on root', &
                       matches)
   end subroutine test_save_spline_plot_mpi


   subroutine test_mpi_singleton_resolution()
      type(PlotParams) :: params
      real(kind=8) :: ux(0:1), uy(0:1), uz(0:1), coeffs(0:0, 0:0, 0:0)
      real(kind=8) :: expected
      integer :: root
      logical :: matches

      ux = (/0.d0, 1.d0/)
      uy = (/0.d0, 1.d0/)
      uz = (/0.d0, 1.d0/)
      coeffs = 5.d0
      params = PlotParams(2.d0, 9.d0, -3.d0, 7.d0, 4.d0, 11.d0, 1, 1, 1)
      root = nranks - 1
      call reset_capture()
      call SaveSplinePlotMPI('singleton', &
                             ux, 0, 0, 1, uy, 0, 0, 1, uz, 0, 0, 1, &
                             coeffs, capture_output, params, root, MPI_COMM_WORLD)

      matches = callback_calls == merge(1, 0, myrank == root)
      if (myrank == root) then
         expected = 5.d0 + 100.d0*2.d0 + 10.d0*(-3.d0) + 4.d0
         matches = matches .and. allocated(captured_values)
         if (matches) then
            matches = all(shape(captured_values) == (/1, 1, 1/))
            matches = matches .and. ieee_is_finite(captured_values(1, 1, 1))
            matches = matches .and. captured_values(1, 1, 1) == expected
         end if
      end if
      call assert_true('SaveSplinePlotMPI uses lower bounds for singleton resolution', &
                       matches)
   end subroutine test_mpi_singleton_resolution


   function coordinate(index, resolution, first, last) result(value)
      integer, intent(in) :: index, resolution
      real(kind=8), intent(in) :: first, last
      real(kind=8) :: value

      if (resolution > 1) then
         value = (1.d0 - dble(index - 1)/dble(resolution - 1))*first + &
                 dble(index - 1)/dble(resolution - 1)*last
      else
         value = first
      end if
   end function coordinate


   function params_match(actual, expected) result(matches)
      type(PlotParams), intent(in) :: actual, expected
      logical :: matches

      matches = actual%startx == expected%startx .and. &
                actual%endx == expected%endx .and. &
                actual%starty == expected%starty .and. &
                actual%endy == expected%endy .and. &
                actual%startz == expected%startz .and. &
                actual%endz == expected%endz .and. &
                actual%resx == expected%resx .and. &
                actual%resy == expected%resy .and. &
                actual%resz == expected%resz
   end function params_match


   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      checks = checks + 1
      if (.not. condition) failures = failures + 1
      if (myrank == 0) then
         if (condition) then
            write (*, '(A)') 'PASS '//trim(label)
         else
            write (*, '(A)') 'FAIL '//trim(label)
         end if
      end if
   end subroutine assert_true

end program test_plot
