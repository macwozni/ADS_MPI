module vtk
   use plot, only: PlotParams
   implicit none

   integer(kind=4) :: vtk_calls = 0
   character(len=80) :: vtk_filename = ''
   type(PlotParams) :: vtk_params
   real(kind=8), allocatable :: vtk_values(:, :, :)

contains

   subroutine reset_vtk_stub()
      vtk_calls = 0
      vtk_filename = ''
      if (allocated(vtk_values)) deallocate (vtk_values)
   end subroutine reset_vtk_stub


   subroutine VtkOutput(filename, vals, params)
      character(len=*), intent(in) :: filename
      type(PlotParams), intent(in) :: params
      real(kind=8), intent(in) :: vals(params%resx, params%resy, params%resz)

      vtk_calls = vtk_calls + 1
      vtk_filename = filename
      vtk_params = params
      allocate (vtk_values(params%resx, params%resy, params%resz))
      vtk_values = vals
   end subroutine VtkOutput

end module vtk
