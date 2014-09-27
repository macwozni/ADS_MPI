
module plot

implicit none


! -------------------------------------------------------------------
! Structure describing spatial plot parameters, like range and
! resolution.
! -------------------------------------------------------------------
type PlotParams
  ! Argument range
  real (kind=8) :: startx, endx
  real (kind=8) :: starty, endy
  real (kind=8) :: startz, endz

  ! Number of points to sample
  integer :: resx, resy, resz
end type


contains


! -------------------------------------------------------------------
! Linear interpolation between two values.
!
! t      - interpolation parameter
! x, y   - points to interpolate between (t=0 -> x, t=1 -> y)
! -------------------------------------------------------------------
function lerp(t, x, y) result (val)
real (kind=8) :: t, x, y, val

  val = (1 - t) * x + t * y

end function


! -------------------------------------------------------------------
! Generic printing routine. Computes values of the function at
! points specified by passed PlotParams object, and writes them
! to file, using specified output function.
!
! filename   - name of the output file / pattern (passed to specified 
!              output function)
! f          - function to plot
! output     - function used to output calculated values
! params     - various options governing plot parameters
! -------------------------------------------------------------------
subroutine SavePlot(filename, f, output, params)
interface
  function f(x, y, z) result (val)
    real (kind=8) :: x, y, z, val
  end function
  
  subroutine output(filename, vals, params)
    import PlotParams
    character(len=*), intent(in) :: filename
    type (PlotParams) :: params
    real (kind=8) :: vals(params%resx,params%resy,params%resz)
  end subroutine
end interface

character(len=*), intent(in) :: filename
type (PlotParams) :: params
real (kind=8) :: vals(params%resx, params%resy, params%resz)
integer :: ix, iy, iz
real (kind=8) :: x, y, z, t

  ! Compute function values
  do ix = 1, params%resx
    t = (ix - 1) / dble(params%resx-1)
    x = lerp(t, params%startx, params%endx)
    do iy = 1, params%resy
      t = (iy - 1) / dble(params%resy-1)
      y = lerp(t, params%starty, params%endy)
      do iz = 1, params%resz
        t = (iz - 1) / dble(params%resz-1)
        z = lerp(t, params%startz, params%endz)
        vals(ix, iy, iz) = f(x, y, z)
      enddo
    enddo
  enddo

  call output(filename, vals, params)

end subroutine


end module

