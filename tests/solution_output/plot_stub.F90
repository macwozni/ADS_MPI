module plot
   use parallelism, only: MYRANK
   implicit none

   type PlotParams
      real(kind=8) :: startx, endx
      real(kind=8) :: starty, endy
      real(kind=8) :: startz, endz
      integer(kind=4) :: resx, resy, resz
   end type PlotParams

   integer(kind=4) :: spline_plot_calls = 0
   integer(kind=4) :: captured_root = -1
   integer(kind=4) :: captured_comm = -1
   integer(kind=4), dimension(3) :: captured_degree = 0
   integer(kind=4), dimension(3) :: captured_n = 0
   integer(kind=4), dimension(3) :: captured_nelem = 0
   character(len=80) :: captured_filename = ''
   type(PlotParams) :: captured_params
   real(kind=8), allocatable :: captured_ux(:), captured_uy(:), captured_uz(:)
   real(kind=8), allocatable :: captured_coeffs(:, :, :)

contains

   subroutine reset_plot_stub()
      spline_plot_calls = 0
      captured_root = -1
      captured_comm = -1
      captured_degree = 0
      captured_n = 0
      captured_nelem = 0
      captured_filename = ''
      if (allocated(captured_ux)) deallocate (captured_ux)
      if (allocated(captured_uy)) deallocate (captured_uy)
      if (allocated(captured_uz)) deallocate (captured_uz)
      if (allocated(captured_coeffs)) deallocate (captured_coeffs)
   end subroutine reset_plot_stub


   subroutine SaveSplinePlotMPI(filename, &
                                Ux, px, nx, nelemx, &
                                Uy, py, ny, nelemy, &
                                Uz, pz, nz, nelemz, &
                                coeffs, output, params, root, comm)
      interface
         subroutine output(filename, vals, params)
            import PlotParams
            character(len=*), intent(in) :: filename
            type(PlotParams), intent(in) :: params
            real(kind=8), intent(in) :: vals(params%resx, params%resy, params%resz)
         end subroutine output
      end interface
      character(len=*), intent(in) :: filename
      integer(kind=4), intent(in) :: px, nx, nelemx
      integer(kind=4), intent(in) :: py, ny, nelemy
      integer(kind=4), intent(in) :: pz, nz, nelemz
      real(kind=8), intent(in) :: Ux(0:nx + px + 1)
      real(kind=8), intent(in) :: Uy(0:ny + py + 1)
      real(kind=8), intent(in) :: Uz(0:nz + pz + 1)
      real(kind=8), intent(in) :: coeffs(0:nx, 0:ny, 0:nz)
      type(PlotParams), intent(in) :: params
      integer(kind=4), intent(in) :: root, comm
      real(kind=8), allocatable :: sampled(:, :, :)
      integer :: ix, iy, iz

      spline_plot_calls = spline_plot_calls + 1
      captured_filename = filename
      captured_root = root
      captured_comm = comm
      captured_degree = (/px, py, pz/)
      captured_n = (/nx, ny, nz/)
      captured_nelem = (/nelemx, nelemy, nelemz/)
      captured_params = params

      allocate (captured_ux(0:nx + px + 1))
      allocate (captured_uy(0:ny + py + 1))
      allocate (captured_uz(0:nz + pz + 1))
      allocate (captured_coeffs(0:nx, 0:ny, 0:nz))
      captured_ux = Ux
      captured_uy = Uy
      captured_uz = Uz
      captured_coeffs = coeffs

      if (MYRANK == root) then
         allocate (sampled(params%resx, params%resy, params%resz))
         do iz = 1, params%resz
            do iy = 1, params%resy
               do ix = 1, params%resx
                  sampled(ix, iy, iz) = coeffs(0, 0, 0) + &
                                        real(10000*ix + 100*iy + iz, kind=8)
               end do
            end do
         end do
         call output(filename, sampled, params)
         deallocate (sampled)
      end if
   end subroutine SaveSplinePlotMPI

end module plot
