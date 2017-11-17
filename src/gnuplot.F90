module gnuplot

   implicit none

contains


   ! -------------------------------------------------------------------
   ! Plot output function for gnuplot.
   !
   ! Input:
   ! ------
   ! filename    - pattern for layer names
   ! vals        - values to plot
   ! params      - plot parameters
   ! -------------------------------------------------------------------
   subroutine GnuPlotOutput(filename, vals, params)
      use plot, ONLY: PlotParams
      implicit none
      character(len = *), intent(in) :: filename
      type (PlotParams), intent(in) :: params
      real (kind = 8), intent(in) :: vals(params % resx, params % resy, params % resz)
      integer :: z

      write(*, *) 'Starting GNUPLOT output...'
      do z = 1, params % resz
         call OutputLayer(filename, z - 1, vals(:,:, z), params)
      enddo
      write(*, *) 'Done with output.'

   end subroutine


   ! -------------------------------------------------------------------
   ! Auxilary function to create file name from pattern and layer number.
   ! Created name is: pattern{layer}.plot.
   !
   ! Input:
   ! ------
   ! pattern     - pattern of the name
   ! layer       - layer number
   ! filename    - output buffer to write name
   ! -------------------------------------------------------------------
   subroutine BuildFileName(pattern, layer, filename)
      implicit none
      character(len = *), intent(in) :: pattern
      integer :: layer
      character(len = *), intent(out) :: filename
      character(len = 10) :: buffer

      write(buffer, '(I10)') layer
      buffer = adjustl(buffer)
      filename = trim(pattern // buffer) // '.plot'

   end subroutine


   ! -------------------------------------------------------------------
   ! Writes one layer as gnuplot file. Its name is: pattern{layer}.plot
   !
   ! Input:
   ! ------
   ! pattern     - pattern for layer names
   ! zlayer      - number of layer
   ! vals        - values to plot (only one layer)
   ! params      - plot parameters
   ! -------------------------------------------------------------------
   subroutine OutputLayer(pattern, zlayer, vals, params)
      use plot, ONLY: PlotParams
      implicit none
      character(len = *), intent(in) :: pattern
      integer, intent(in) :: zlayer
      type (PlotParams), intent(in) :: params
      real (kind = 8), intent(in) :: vals(params % resx, params % resy)

      character(len = 50) :: filename, buf
      integer :: outFile = 57 ! random value, Grothendieck's prime
      integer :: ix, iy

      write(*, *) 'Layer', zlayer
      call BuildFileName(pattern, zlayer, filename)

      open(unit = outFile, file = filename, &
      form = 'formatted', access = 'sequential', status = 'unknown')

      do ix = 1, params % resx
         do iy = 1, params % resy
            write(buf, '(F30.10)') vals(ix, iy)
            write(outFile, *) ix - 1, ' ', iy - 1, ' ', buf
         enddo
         write(outFile, *) ! blank line
      enddo

      close(outFile)

   end subroutine


end module

