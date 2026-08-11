module plot_test_callbacks
   use plot, only: PlotParams
   implicit none

   integer(kind=4) :: callback_calls = 0
   character(len=80) :: captured_filename = ''
   type(PlotParams) :: captured_params
   real(kind=8), allocatable :: captured_values(:, :, :)

contains

   function scalar_function(x, y, z) result(value)
      real(kind=8) :: x, y, z, value

      value = 100.d0*x + 10.d0*y + z
   end function scalar_function


   subroutine capture_output(filename, vals, params)
      character(len=*), intent(in) :: filename
      type(PlotParams), intent(in) :: params
      real(kind=8), intent(in) :: vals(params%resx, params%resy, params%resz)

      callback_calls = callback_calls + 1
      captured_filename = filename
      captured_params = params
      if (allocated(captured_values)) deallocate (captured_values)
      allocate (captured_values(params%resx, params%resy, params%resz))
      captured_values = vals
   end subroutine capture_output


   subroutine reset_capture()
      callback_calls = 0
      captured_filename = ''
      if (allocated(captured_values)) deallocate (captured_values)
   end subroutine reset_capture

end module plot_test_callbacks
