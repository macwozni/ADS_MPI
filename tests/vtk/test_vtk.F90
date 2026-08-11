program test_vtk
   use plot, only: PlotParams
   use vtk, only: VtkOutput, VtkStructuredGridOutput
   implicit none

   integer(kind=4) :: checks, failures

   checks = 0
   failures = 0

   call test_image_data_order()
   call test_image_data_singleton()
   call test_structured_grid_order()

   if (failures == 0) then
      write (*, '(A,I0,A)') 'OK (', checks, ' VTK checks)'
   else
      write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                 ' VTK checks)'
      stop 1
   end if

contains

   subroutine test_image_data_order()
      type(PlotParams) :: params
      real(kind=8) :: vals(2, 2, 2)
      character(len=256) :: line, extent, origin, spacing
      integer :: unit, ios, ix, iy, iz
      logical :: matches

      params = PlotParams(1.d0, 3.d0, -2.d0, 2.d0, 4.d0, 10.d0, 2, 2, 2)
      do iz = 1, 2
         do iy = 1, 2
            do ix = 1, 2
               vals(ix, iy, iz) = real(100*iz + 10*iy + ix, kind=8)
            end do
         end do
      end do

      call delete_file('vtk_image.vti')
      call VtkOutput('vtk_image', vals, params)
      call assert_file_exists('VtkOutput appends the .vti extension', &
                              'vtk_image.vti')

      unit = 81
      open (unit=unit, file='vtk_image.vti', status='old', action='read', &
            form='formatted')
      matches = .true.
      call expect_line(unit, '<?xml version="1.0"?>', matches)
      call expect_line(unit, '<VTKFile type="ImageData" version="0.1">', matches)
      write (extent, '("0 ",(I5)," 0 ",(I5)," 0 ",(I5))') 1, 1, 1
      write (origin, '(3(ES24.16,1X))') 1.d0, -2.d0, 4.d0
      write (spacing, '(3(ES24.16,1X))') 2.d0, 4.d0, 6.d0
      call expect_line(unit, '  <ImageData WholeExtent="'//trim(extent)// &
                       '" Origin="'//trim(adjustl(origin))//'" Spacing="'// &
                       trim(adjustl(spacing))//'">', matches)
      call expect_line(unit, '    <Piece Extent="'//trim(extent)//'">', matches)
      call expect_line(unit, '      <PointData Scalars="Result">', matches)
      call expect_line(unit, &
                       '        <DataArray Name="Result" type="Float64" format="ascii">', &
                       matches)
      do iz = 1, 2
         do iy = 1, 2
            do ix = 1, 2
               write (line, '(F30.10)') vals(ix, iy, iz)
               call expect_line(unit, '          '//trim(adjustl(line)), matches)
            end do
         end do
      end do
      call expect_line(unit, '        </DataArray>', matches)
      call expect_line(unit, '      </PointData>', matches)
      call expect_line(unit, '    </Piece>', matches)
      call expect_line(unit, '  </ImageData>', matches)
      call expect_line(unit, '</VTKFile>', matches)
      read (unit, '(A)', iostat=ios) line
      matches = matches .and. ios < 0
      close (unit)

      call assert_true('VtkOutput emits exact XML and X-fastest point order', matches)
      call delete_file('vtk_image.vti')
   end subroutine test_image_data_order


   subroutine test_image_data_singleton()
      type(PlotParams) :: params
      real(kind=8) :: vals(1, 1, 1)
      character(len=256) :: line, extent, origin, spacing
      integer :: unit, ios, record
      logical :: matches

      params = PlotParams(-1.d0, 8.d0, 2.d0, 9.d0, 5.d0, 12.d0, 1, 1, 1)
      vals = -3.5d0
      call delete_file('vtk_single.vti')
      call VtkOutput('vtk_single', vals, params)

      write (extent, '("0 ",(I5)," 0 ",(I5)," 0 ",(I5))') 0, 0, 0
      write (origin, '(3(ES24.16,1X))') -1.d0, 2.d0, 5.d0
      write (spacing, '(3(ES24.16,1X))') 0.d0, 0.d0, 0.d0
      unit = 82
      open (unit=unit, file='vtk_single.vti', status='old', action='read', &
            form='formatted')
      matches = .true.
      do record = 1, 2
         read (unit, '(A)', iostat=ios) line
         matches = matches .and. ios == 0
      end do
      call expect_line(unit, '  <ImageData WholeExtent="'//trim(extent)// &
                       '" Origin="'//trim(adjustl(origin))//'" Spacing="'// &
                       trim(adjustl(spacing))//'">', matches)
      close (unit)

      call assert_true('VtkOutput represents singleton axes with zero spacing', &
                       matches)
      call delete_file('vtk_single.vti')
   end subroutine test_image_data_singleton


   subroutine test_structured_grid_order()
      type(PlotParams) :: params
      real(kind=8) :: vals(2, 1, 2), x(2, 1, 2), y(2, 1, 2), z(2, 1, 2)
      character(len=256) :: line, extent
      integer :: unit, ios, ix, iy, iz
      logical :: matches

      params = PlotParams(0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 2, 1, 2)
      do ix = 1, 2
         do iy = 1, 1
            do iz = 1, 2
               vals(ix, iy, iz) = real(100*ix + 10*iy + iz, kind=8)
               x(ix, iy, iz) = real(ix, kind=8) + 0.1d0*iz
               y(ix, iy, iz) = real(iy, kind=8) + 0.2d0*ix
               z(ix, iy, iz) = real(iz, kind=8) + 0.3d0*iy
            end do
         end do
      end do

      call delete_file('vtk_structured.vts')
      call VtkStructuredGridOutput('vtk_structured', vals, x, y, z, params)
      call assert_file_exists('VtkStructuredGridOutput appends the .vts extension', &
                              'vtk_structured.vts')

      unit = 83
      open (unit=unit, file='vtk_structured.vts', status='old', &
            action='read', form='formatted')
      matches = .true.
      call expect_line(unit, '<?xml version="1.0"?>', matches)
      call expect_line(unit, '<VTKFile type="StructuredGrid" version="0.1">', matches)
      write (extent, '("0 ",(I5)," 0 ",(I5)," 0 ",(I5))') 1, 0, 1
      call expect_line(unit, '  <StructuredGrid WholeExtent="'//trim(extent)// &
                       '" origin="0 0 0" spacing="1 1 1">', matches)
      call expect_line(unit, '    <Piece Extent="'//trim(extent)//'">', matches)
      call expect_line(unit, '      <PointData Scalars="Result">', matches)
      call expect_line(unit, '        <DataArray Name="Result" type="Float32">', &
                       matches)
      do ix = 1, 2
         do iy = 1, 1
            do iz = 1, 2
               write (line, '(F30.10)') vals(ix, iy, iz)
               call expect_line(unit, '          '//trim(adjustl(line)), matches)
            end do
         end do
      end do
      call expect_line(unit, '        </DataArray>', matches)
      call expect_line(unit, '      </PointData>', matches)
      call expect_line(unit, '      <Points>', matches)
      call expect_line(unit, &
                       '        <DataArray type="Float32" NumberOfComponents="3">', &
                       matches)
      do ix = 1, 2
         do iy = 1, 1
            do iz = 1, 2
               write (line, '(3F30.10)') x(ix, iy, iz), y(ix, iy, iz), &
                                          z(ix, iy, iz)
               call expect_line(unit, '          '//trim(adjustl(line)), matches)
            end do
         end do
      end do
      call expect_line(unit, '        </DataArray>', matches)
      call expect_line(unit, '      </Points>', matches)
      call expect_line(unit, '    </Piece>', matches)
      call expect_line(unit, '  </StructuredGrid>', matches)
      call expect_line(unit, '</VTKFile>', matches)
      read (unit, '(A)', iostat=ios) line
      matches = matches .and. ios < 0
      close (unit)

      call assert_true('structured VTK data and coordinates use Z-fastest order', &
                       matches)
      call delete_file('vtk_structured.vts')
   end subroutine test_structured_grid_order


   subroutine expect_line(unit, expected, matches)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: expected
      logical, intent(inout) :: matches
      character(len=512) :: actual
      integer :: ios

      read (unit, '(A)', iostat=ios) actual
      matches = matches .and. ios == 0 .and. trim(actual) == trim(expected)
   end subroutine expect_line


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

end program test_vtk
