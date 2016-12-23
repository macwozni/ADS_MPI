module vtk

use plot
use debug

implicit none

contains


! -------------------------------------------------------------------
! Data output function for VTK
!
! filename    - name of the output file
! vals        - values to plot
! params      - plot parameters
! -------------------------------------------------------------------
subroutine VtkOutput(filename, vals, params)
character(len=*), intent(in) :: filename
type (PlotParams), intent(in) :: params
real (kind=8), intent(in) :: vals(params%resx,params%resy,params%resz)
integer :: ix, iy, iz
character (len=200) :: temp, extent

integer :: outFile = 57 ! random value, Grothendieck's prime

  if (iprint == 1) then
    write(*,*) 'Starting VTK output...'
  endif

  open(unit=outFile, file=trim(filename) // '.vti', &
    form='formatted', access='sequential', status='unknown')

  ! XML version/root
  write(outFile,'(A)') '<?xml version="1.0"?>'
  write(outFile,'(A)') '<VTKFile type="ImageData" version="0.1">'

  ! Prepare extent (count of parts in each dimension) for later use
  ! Rormat: "x1 x2 y1 y2 z1 z2"
  write(extent, '("0 ",(I5)," 0 ",(I5)," 0 ",(I5))') &
    params%resx-1, params%resy-1, params%resz-1

  ! Init ImageData structure for whole region
  temp = '  <ImageData WholeExtent="' // trim(extent)
  write(outFile,'(A)') trim(temp) // '" origin="0 0 0" spacing="1 1 1">'

  ! Region consists of one piece in one file
  write(outFile,'(A)')'    <Piece Extent="' // trim(extent) // '">'
  write(outFile,'(A)')'      <PointData Scalars="Result">'
  write(outFile,'(A)')'        <DataArray Name="Result" type="Float32">'

  ! outFile result values by X, Y and Z axis
  if (iprint == 1) then
    write(*,*)'Printing data'
  endif

  do iz = 1, params%resz
    do iy = 1, params%resy
      do ix = 1, params%resx
        write(temp,'(F30.10)') vals(ix, iy, iz)
        write(outFile,'(A)') '          ' // trim(adjustl(temp))
      enddo
    enddo
  enddo

  ! Close all XML nodes
  write(outFile,'(A)') '        </DataArray>'
  write(outFile,'(A)') '      </PointData>'
  write(outFile,'(A)') '    </Piece>'
  write(outFile,'(A)') '  </ImageData>'
  write(outFile,'(A)') '</VTKFile>'

  close(outFile)
  if (iprint == 1) then
    write(*,*) 'Done with output.'
  endif

end subroutine


! -------------------------------------------------------------------
! Data output function for VTK, version for structured grid
!
! filename    - name of the output file
! vals        - values to plot
! X, Y, Z     - coordinates of grid points
! params      - plot parameters
! -------------------------------------------------------------------
subroutine VtkStructuredGridOutput(filename, vals, X, Y, Z, params)
character(len=*), intent(in) :: filename
type (PlotParams), intent(in) :: params
real (kind=8), dimension(params%resx,params%resy,params%resz), &
  intent(in) :: X, Y, Z, vals
integer :: ix, iy, iz
character(len=200) :: temp, extent

integer :: outFile = 57 ! random value, Grothendieck's prime

  if (iprint == 1) then
    write(*,*) 'Starting VTK output...'
  endif

  open(unit=outFile, file=trim(filename) // '.vts', &
    form='formatted', access='sequential', status='unknown')

  ! XML version/root
  write(outFile,'(A)') '<?xml version="1.0"?>'
  write(outFile,'(A)') '<VTKFile type="StructuredGrid" version="0.1">'

  ! Prepare extent (count of parts in each dimension) for later use
  ! Rormat: "x1 x2 y1 y2 z1 z2"
  write(extent, '("0 ",(I5)," 0 ",(I5)," 0 ",(I5))') &
    params%resx-1, params%resy-1, params%resz-1

  ! Init ImageData structure for whole region
  temp = '  <StructuredGrid WholeExtent="' // trim(extent)
  write(outFile,'(A)') trim(temp) // '" origin="0 0 0" spacing="1 1 1">'

  ! Region consists of one piece in one file
  write(outFile,'(A)')'    <Piece Extent="' // trim(extent) // '">'
  write(outFile,'(A)')'      <PointData Scalars="Result">'
  write(outFile,'(A)')'        <DataArray Name="Result" type="Float32">'

  ! output result values by X, Y and Z axis
  if (iprint == 1) then
    write(*,*)'Printing data'
  endif

  do ix = 1, params%resx
    do iy = 1, params%resy
      do iz = 1, params%resz
        write(temp,'(F30.10)') vals(ix, iy, iz)
        write(outFile,'(A)') '          ' // trim(adjustl(temp))
      enddo
    enddo
  enddo

  write(outFile,'(A)') '        </DataArray>'
  write(outFile,'(A)') '      </PointData>'
  write(outFile,'(A)')'      <Points>'
  write(outFile,'(A)') '        <DataArray type="Float32" NumberOfComponents="3">'

  ! outFile result values by X, Y and Z axis
  do ix = 1, params%resx
    do iy = 1, params%resy
      do iz = 1, params%resz
        write(temp,'(3F30.10)') X(ix, iy, iz), Y(ix, iy, iz), Z(ix, iy, iz)
        write(outFile,'(A)') '          ' // trim(adjustl(temp))
      enddo
    enddo
  enddo

  write(outFile,'(A)') '        </DataArray>'
  write(outFile,'(A)') '      </Points>'
  write(outFile,'(A)') '    </Piece>'
  write(outFile,'(A)') '  </StructuredGrid>'
  write(outFile,'(A)') '</VTKFile>'

  close(outFile)
  if (iprint == 1) then
    write(*,*) 'Done with output.'
  endif

end subroutine


end module

