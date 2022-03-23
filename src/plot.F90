! -------------------------------------------------------------------
! Module for ploting graphics.
! -------------------------------------------------------------------

module plot

   implicit none

   ! -------------------------------------------------------------------
   ! Structure describing spatial plot parameters, like range and
   ! resolution.
   ! -------------------------------------------------------------------
   type PlotParams
      ! Argument range
      real(kind=8) :: startx, endx
      real(kind=8) :: starty, endy
      real(kind=8) :: startz, endz

      ! Number of points to sample
      integer(kind=4) :: resx, resy, resz
   end type PlotParams

contains

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
      use math, ONLY: lerp
      implicit none
      interface
         function f(x, y, z) result(val)
            real(kind=8) :: x, y, z, val
         end function

         subroutine output(filename, vals, params)
            import PlotParams
            character(len=*), intent(in) :: filename
            type(PlotParams) :: params
            real(kind=8) :: vals(params%resx, params%resy, params%resz)
         end subroutine
      end interface

      character(len=*), intent(in) :: filename
      type(PlotParams) :: params
      real(kind=8) :: vals(params%resx, params%resy, params%resz)
      integer(kind=4) :: ix, iy, iz
      real(kind=8) :: x, y, z, t

      ! Compute function values
      do ix = 1, params%resx
         t = (ix - 1)/dble(params%resx - 1)
         x = lerp(t, params%startx, params%endx)
         do iy = 1, params%resy
            t = (iy - 1)/dble(params%resy - 1)
            y = lerp(t, params%starty, params%endy)
            do iz = 1, params%resz
               t = (iz - 1)/dble(params%resz - 1)
               z = lerp(t, params%startz, params%endz)
               vals(ix, iy, iz) = f(x, y, z)
            end do
         end do
      end do

      call output(filename, vals, params)

   end subroutine SavePlot

   ! -------------------------------------------------------------------
   ! Spline function plotting routine. Computes values of the function
   ! given by knot specification and coefficients of basis spline
   ! function at points specified by PlotParams object, and writes them
   ! to file, using specified output function.
   !
   ! filename   - name of the output file / pattern (passed to specified
   !              output function)
   ! U_         - knot points
   ! p_         - degree of splines
   ! n_         - nubmer of functions minus one
   ! nelem_     - number of elements
   ! coeffs     - 3D array of coefficients of basis functions (0-based)
   ! output     - function used to output calculated values
   ! params     - various options governing plot parameters
   ! -------------------------------------------------------------------
   subroutine SaveSplinePlot(filename, &
                             Ux, px, nx, nelemx, &
                             Uy, py, ny, nelemy, &
                             Uz, pz, nz, nelemz, &
                             coeffs, output, params)
      use math, ONLY: lerp
      use basis, ONLY: EvalSpline
      interface
         subroutine output(filename, vals, params)
            import PlotParams
            character(len=*), intent(in) :: filename
            type(PlotParams), intent(in) :: params
            real(kind=8), intent(in) :: vals(params%resx, params%resy, params%resz)
         end subroutine
      end interface

      character(len=*), intent(in) :: filename
      integer(kind=4), intent(in) :: nx, px, nelemx
      integer(kind=4), intent(in) :: ny, py, nelemy
      integer(kind=4), intent(in) :: nz, pz, nelemz
      real(kind=8), intent(in) :: Ux(0:nx + px + 1)
      real(kind=8), intent(in) :: Uy(0:ny + py + 1)
      real(kind=8), intent(in) :: Uz(0:nz + pz + 1)
      real(kind=8), intent(in) :: coeffs(0:nx, 0:ny, 0:nz)
      type(PlotParams) :: params
      real(kind=8) :: vals(params%resx, params%resy, params%resz)
      real(kind=8) :: x, y, z, t
      integer(kind=4) :: ix, iy, iz

      ! Compute function values
      do ix = 1, params%resx
         t = (ix - 1)/dble(params%resx - 1)
         x = lerp(t, params%startx, params%endx)
         do iy = 1, params%resy
            t = (iy - 1)/dble(params%resy - 1)
            y = lerp(t, params%starty, params%endy)
            do iz = 1, params%resz
               t = (iz - 1)/dble(params%resz - 1)
               z = lerp(t, params%startz, params%endz)
               !vals(ix, iy, iz) = f(x, y, z)
               vals(ix, iy, iz) = EvalSpline(0, Ux, px, nx, nelemx, Uy, py, ny, nelemy, &
                                             Uz, pz, nz, nelemz, coeffs, x, y, z)
            end do
         end do
      end do

      call output(filename, vals, params)

   end subroutine SaveSplinePlot

end module plot

