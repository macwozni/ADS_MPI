!------------------------------------------------------------------------------
!
! MODULE: vtk
!
! DESCRIPTION:
!> @file vtk.F90
!> @brief Module providing routines for exporting scalar data to VTK XML
!> formats.
!>
!> @details
!> This module groups output procedures used to write numerical results
!> in VTK-compatible XML files for later visualization.
!>
!> The provided functionality includes:
!> - export to image-data format through \ref VtkOutput,
!> - export to structured-grid format through
!>   \ref VtkStructuredGridOutput.
!>
!> The routines serialize scalar fields sampled on regular or structured
!> three-dimensional grids and store them in files that can be processed
!> by standard visualization tools supporting VTK XML formats.
!
!------------------------------------------------------------------------------
module vtk

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Writes a scalar field to a VTK XML image-data file.
!>
!> @details
!> This procedure exports the array \p vals to a file in VTK XML image
!> format with extension `.vti`. The output is interpreted as point data
!> on a regular Cartesian grid whose dimensions are specified by the
!> plotting-parameter object \p params.
!>
!> The routine:
!> - opens an output file named `filename.vti`,
!> - constructs the VTK `WholeExtent` descriptor from the grid
!>   resolutions,
!> - writes the scalar field under the `PointData` section using the
!>   array name `Result`,
!> - closes all XML nodes and the file handle.
!>
!> The coordinate origin and spacing are derived from the sampling bounds
!> and resolution stored in \p params.
!
! Input:
! ------
!> @param[in] filename
!> Base name of the output file, without extension.
!>
!> @param[in] vals
!> Scalar values to be written as VTK point data.
!>
!> @param[in] params
!> Plot-parameter structure defining the grid resolution.
!
! Notes:
! ------
!> @note
!> The file is written in ASCII formatted XML form.
!
!> @note
!> The VTK data type is declared as `Float64`, matching the
!> `real(kind=8)` input data written in ASCII form.
!
!---------------------------------------------------------------------------
subroutine VtkOutput(filename, vals, params)
   use plot, ONLY: PlotParams
   implicit none
!> @brief Base name of the output file.
   character(len=*), intent(in) :: filename
!> @brief Plot parameters defining the grid resolution.
   type(PlotParams), intent(in) :: params
!> @brief Scalar field values written to the VTK file.
   real(kind=8), intent(in) :: vals(params%resx, params%resy, params%resz)
!> @brief Loop counters over the grid directions.
   integer(kind=4) :: ix, iy, iz
!> @brief Temporary buffers for formatted numeric and extent strings.
   character(len=200) :: temp, extent, origin, spacing
!> @brief Physical grid spacing in each coordinate direction.
   real(kind=8) :: hx, hy, hz

!> @brief Output unit number used for file writing.
   integer :: outFile = 57

#ifdef IPRINT
   write (*, *) 'Starting VTK output...'
#endif

   open (unit=outFile, file=trim(filename)//'.vti', &
         form='formatted', access='sequential', status='replace', &
         action='write')

   ! XML version/root
   write (outFile, '(A)') '<?xml version="1.0"?>'
   write (outFile, '(A)') '<VTKFile type="ImageData" version="0.1">'

   ! Prepare extent (count of parts in each dimension) for later use
   ! Format: "x1 x2 y1 y2 z1 z2"
   write (extent, '("0 ",(I5)," 0 ",(I5)," 0 ",(I5))') &
      params%resx - 1, params%resy - 1, params%resz - 1
   if (params%resx > 1) then
      hx = (params%endx - params%startx)/dble(params%resx - 1)
   else
      hx = 0.d0
   end if
   if (params%resy > 1) then
      hy = (params%endy - params%starty)/dble(params%resy - 1)
   else
      hy = 0.d0
   end if
   if (params%resz > 1) then
      hz = (params%endz - params%startz)/dble(params%resz - 1)
   else
      hz = 0.d0
   end if
   write (origin, '(3(ES24.16,1X))') params%startx, params%starty, params%startz
   write (spacing, '(3(ES24.16,1X))') hx, hy, hz

   ! Init ImageData structure for whole region
   temp = '  <ImageData WholeExtent="'//trim(extent)
   write (outFile, '(A)') trim(temp)//'" Origin="'//trim(adjustl(origin))// &
      '" Spacing="'//trim(adjustl(spacing))//'">'

   ! Region consists of one piece in one file
   write (outFile, '(A)') '    <Piece Extent="'//trim(extent)//'">'
   write (outFile, '(A)') '      <PointData Scalars="Result">'
   write (outFile, '(A)') '        <DataArray Name="Result" type="Float64" format="ascii">'

   ! outFile result values by X, Y and Z axis
#ifdef IPRINT
   write (*, *) 'Printing data'
#endif

   do iz = 1, params%resz
      do iy = 1, params%resy
         do ix = 1, params%resx
            write (temp, '(F30.10)') vals(ix, iy, iz)
            write (outFile, '(A)') '          '//trim(adjustl(temp))
         end do
      end do
   end do

   ! Close all XML nodes
   write (outFile, '(A)') '        </DataArray>'
   write (outFile, '(A)') '      </PointData>'
   write (outFile, '(A)') '    </Piece>'
   write (outFile, '(A)') '  </ImageData>'
   write (outFile, '(A)') '</VTKFile>'

   close (outFile)
#ifdef IPRINT
   write (*, *) 'Done with output.'
#endif

end subroutine VtkOutput

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Writes a scalar field and point coordinates to a VTK XML
!> structured-grid file.
!>
!> @details
!> This procedure exports the array \p vals together with the point
!> coordinates \p X, \p Y, and \p Z to a VTK XML structured-grid file
!> with extension `.vts`.
!>
!> The routine:
!> - opens an output file named `filename.vts`,
!> - constructs the VTK `WholeExtent` descriptor from the grid
!>   resolutions,
!> - writes the scalar field under the `PointData` section using the
!>   array name `Result`,
!> - writes the physical coordinates of all grid points under the
!>   `Points` section as a three-component data array,
!> - closes all XML nodes and the file handle.
!>
!> This format is suitable for grids whose topology is Cartesian but
!> whose point coordinates are not necessarily uniformly spaced.
!
! Input:
! ------
!> @param[in] filename
!> Base name of the output file, without extension.
!>
!> @param[in] vals
!> Scalar values to be written as VTK point data.
!>
!> @param[in] X
!> First Cartesian coordinates of the structured-grid points.
!>
!> @param[in] Y
!> Second Cartesian coordinates of the structured-grid points.
!>
!> @param[in] Z
!> Third Cartesian coordinates of the structured-grid points.
!>
!> @param[in] params
!> Plot-parameter structure defining the grid resolution.
!
! Notes:
! ------
!> @note
!> The file is written in ASCII formatted XML form.
!
!> @warning
!> The VTK data type declared in the XML output is `Float32`, whereas the
!> input arrays are stored as `real(kind=8)`.
!
!---------------------------------------------------------------------------
subroutine VtkStructuredGridOutput(filename, vals, X, Y, Z, params)
   use plot, ONLY: PlotParams
   implicit none
   !> @brief Base name of the output file.
   character(len=*), intent(in) :: filename
   !> @brief Plot parameters defining the grid resolution.
   type(PlotParams), intent(in) :: params
   !> @brief Structured-grid coordinates and scalar values.
   real(kind=8), dimension(params%resx, params%resy, params%resz), &
      intent(in) :: X, Y, Z, vals
   !> @brief Loop counters over the grid directions.
   integer(kind=4) :: ix, iy, iz
   !> @brief Temporary buffers for formatted numeric and extent strings.
   character(len=200) :: temp, extent

   !> @brief Output unit number used for file writing.
   integer :: outFile = 57

#ifdef IPRINT
   write (*, *) 'Starting VTK output...'
#endif

   open (unit=outFile, file=trim(filename)//'.vts', &
         form='formatted', access='sequential', status='unknown')

   ! XML version/root
   write (outFile, '(A)') '<?xml version="1.0"?>'
   write (outFile, '(A)') '<VTKFile type="StructuredGrid" version="0.1">'

   ! Prepare extent (count of parts in each dimension) for later use
   ! Rormat: "x1 x2 y1 y2 z1 z2"
   write (extent, '("0 ",(I5)," 0 ",(I5)," 0 ",(I5))') &
      params%resx - 1, params%resy - 1, params%resz - 1

   ! Init ImageData structure for whole region
   temp = '  <StructuredGrid WholeExtent="'//trim(extent)
   write (outFile, '(A)') trim(temp)//'" origin="0 0 0" spacing="1 1 1">'

   ! Region consists of one piece in one file
   write (outFile, '(A)') '    <Piece Extent="'//trim(extent)//'">'
   write (outFile, '(A)') '      <PointData Scalars="Result">'
   write (outFile, '(A)') '        <DataArray Name="Result" type="Float32">'

   ! output result values by X, Y and Z axis
#ifdef IPRINT
   write (*, *) 'Printing data'
#endif

   do ix = 1, params%resx
      do iy = 1, params%resy
         do iz = 1, params%resz
            write (temp, '(F30.10)') vals(ix, iy, iz)
            write (outFile, '(A)') '          '//trim(adjustl(temp))
         end do
      end do
   end do

   write (outFile, '(A)') '        </DataArray>'
   write (outFile, '(A)') '      </PointData>'
   write (outFile, '(A)') '      <Points>'
   write (outFile, '(A)') '        <DataArray type="Float32" NumberOfComponents="3">'

   ! outFile result values by X, Y and Z axis
   do ix = 1, params%resx
      do iy = 1, params%resy
         do iz = 1, params%resz
            write (temp, '(3F30.10)') X(ix, iy, iz), Y(ix, iy, iz), Z(ix, iy, iz)
            write (outFile, '(A)') '          '//trim(adjustl(temp))
         end do
      end do
   end do

   write (outFile, '(A)') '        </DataArray>'
   write (outFile, '(A)') '      </Points>'
   write (outFile, '(A)') '    </Piece>'
   write (outFile, '(A)') '  </StructuredGrid>'
   write (outFile, '(A)') '</VTKFile>'

   close (outFile)
#ifdef IPRINT
   write (*, *) 'Done with output.'
#endif

end subroutine VtkStructuredGridOutput

end module vtk
