module plot
   implicit none

   type PlotParams
      real(kind=8) :: startx, endx
      real(kind=8) :: starty, endy
      real(kind=8) :: startz, endz
      integer(kind=4) :: resx, resy, resz
   end type PlotParams

end module plot
