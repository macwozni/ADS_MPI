program test_gnuplot
   use gnuplot, only: BuildFileName, GnuPlotOutput, OutputLayer
   use plot, only: PlotParams
   implicit none

   integer(kind=4) :: checks, failures

   checks = 0
   failures = 0

   call test_build_file_name()
   call test_singleton_layer()
   call test_layer_stack_order()

   if (failures == 0) then
      write (*, '(A,I0,A)') 'OK (', checks, ' gnuplot checks)'
   else
      write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                 ' gnuplot checks)'
      stop 1
   end if

contains

   subroutine test_build_file_name()
      character(len=80) :: filename

      call BuildFileName('sample_', 12, filename)
      call assert_true('BuildFileName appends the layer and extension', &
                       trim(filename) == 'sample_12.plot')

      call BuildFileName('negative_', -3, filename)
      call assert_true('BuildFileName preserves a negative layer', &
                       trim(filename) == 'negative_-3.plot')
   end subroutine test_build_file_name


   subroutine test_singleton_layer()
      type(PlotParams) :: params
      real(kind=8) :: vals(1, 1)
      character(len=256) :: actual, expected
      integer :: unit, ios

      params = PlotParams(-2.d0, 9.d0, 3.d0, 8.d0, 4.d0, 7.d0, 1, 1, 1)
      vals(1, 1) = 7.25d0
      call delete_file('gnuplot_single_12.plot')

      call OutputLayer('gnuplot_single_', 12, vals, params)
      call assert_file_exists('OutputLayer creates the exact layer filename', &
                              'gnuplot_single_12.plot')

      unit = 81
      open (unit=unit, file='gnuplot_single_12.plot', status='old', &
            action='read', form='formatted')
      read (unit, '(A)', iostat=ios) actual
      write (expected, '(3(ES24.16,1X))') -2.d0, 3.d0, 7.25d0
      call assert_true('singleton layer uses the lower-bound coordinates', &
                       ios == 0 .and. trim(actual) == trim(expected))
      read (unit, '(A)', iostat=ios) actual
      call assert_true('singleton layer terminates its X sweep with a blank line', &
                       ios == 0 .and. len_trim(actual) == 0)
      read (unit, '(A)', iostat=ios) actual
      call assert_true('singleton layer contains no extra records', ios < 0)
      close (unit)

      call delete_file('gnuplot_single_12.plot')
   end subroutine test_singleton_layer


   subroutine test_layer_stack_order()
      type(PlotParams) :: params
      real(kind=8) :: vals(2, 2, 2)
      character(len=80) :: filename
      integer :: ix, iy, iz

      params = PlotParams(0.d0, 2.d0, 10.d0, 20.d0, -1.d0, 1.d0, 2, 2, 2)
      do iz = 1, 2
         do iy = 1, 2
            do ix = 1, 2
               vals(ix, iy, iz) = real(100*iz + 10*ix + iy, kind=8)
            end do
         end do
      end do

      call delete_file('gnuplot_stack_0.plot')
      call delete_file('gnuplot_stack_1.plot')
      call GnuPlotOutput('gnuplot_stack_', vals, params)

      do iz = 1, 2
         call BuildFileName('gnuplot_stack_', iz - 1, filename)
         call assert_file_exists('GnuPlotOutput creates every zero-based layer', &
                                 trim(filename))
         call check_stack_file(trim(filename), iz, vals(:, :, iz), params)
         call delete_file(trim(filename))
      end do
   end subroutine test_layer_stack_order


   subroutine check_stack_file(filename, layer, vals, params)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: layer
      real(kind=8), intent(in) :: vals(2, 2)
      type(PlotParams), intent(in) :: params
      character(len=256) :: actual, expected
      real(kind=8) :: x, y
      integer :: unit, ios, ix, iy
      logical :: matches

      unit = 82
      open (unit=unit, file=filename, status='old', action='read', form='formatted')
      matches = .true.
      do ix = 1, params%resx
         x = params%startx + (params%endx - params%startx)* &
             dble(ix - 1)/dble(params%resx - 1)
         do iy = 1, params%resy
            y = params%starty + (params%endy - params%starty)* &
                dble(iy - 1)/dble(params%resy - 1)
            read (unit, '(A)', iostat=ios) actual
            write (expected, '(3(ES24.16,1X))') x, y, vals(ix, iy)
            matches = matches .and. ios == 0 .and. &
                      trim(actual) == trim(expected)
         end do
         read (unit, '(A)', iostat=ios) actual
         matches = matches .and. ios == 0 .and. len_trim(actual) == 0
      end do
      read (unit, '(A)', iostat=ios) actual
      matches = matches .and. ios < 0
      close (unit)

      call assert_true('layer data order is X sweeps containing Y rows, layer '// &
                       digit(layer - 1), matches)
   end subroutine check_stack_file


   function digit(value) result(text)
      integer, intent(in) :: value
      character(len=16) :: text

      write (text, '(I0)') value
   end function digit


   subroutine assert_file_exists(label, filename)
      character(len=*), intent(in) :: label, filename
      logical :: exists

      inquire (file=filename, exist=exists)
      call assert_true(label, exists)
   end subroutine assert_file_exists


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


   subroutine delete_file(filename)
      character(len=*), intent(in) :: filename
      logical :: exists
      integer :: unit, ios

      inquire (file=filename, exist=exists)
      if (.not. exists) return
      unit = 89
      open (unit=unit, file=filename, status='old', iostat=ios)
      if (ios == 0) close (unit, status='delete')
   end subroutine delete_file

end program test_gnuplot
