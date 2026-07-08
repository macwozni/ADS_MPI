!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Input data, material coefficients, and forcing callback for the
!> oil transport example.
!>
!> @details
!> This module stores the oil-problem configuration, including randomly
!> generated channel curves, pump and drain locations, time-step data, and
!> precomputed permeability samples. It also provides the pointwise forcing
!> callback used by the current ADS API.
!
!------------------------------------------------------------------------------
module input_data

   implicit none

!> @brief Number of generated curves and number of segments per curve.
   integer(kind = 4), parameter :: cN = 30, cL = 16
!> @brief Pump/drain support radius and source/sink strengths.
   real (kind = 8), parameter :: radius = 0.15, pumping_strength = 1, draining_strength = 1

!> @brief Coordinates of piecewise-linear channel curves.
   real (kind = 8) :: cx(cN * cL), cy(cN * cL), cz(cN * cL)
!> @brief Nonlinear mobility exponent used by the legacy oil model.
   real (kind = 8) :: mi = 10.d0
!> @brief Ground-level parameter retained from the original problem data.
   real (kind = 8) :: GROUND = 0.2
!> @brief Minimum and maximum permeability values.
   real (kind = 8), parameter :: Kqmin = 1.d0, Kqmax = 1000.d0
!> @brief Numbers of pump and drain points.
   integer(kind = 4) :: npumps, ndrains
!> @brief Pump and drain coordinates stored by column.
   real (kind = 8), allocatable, dimension(:,:) :: pumps, drains

!> @brief Cached permeability values at local quadrature points.
   real (kind = 8), allocatable :: Kqvals(:,:,:,:,:,:)

!> @brief Current physical time.
   real (kind = 8) :: t

!> @brief Time-step size.
   real (kind = 8) :: Dt

!> @brief Number of time iterations.
   integer :: steps

!> @brief Accumulated pollution statistic retained for the problem model.
   real (kind = 8) :: pollution = 0

!> @brief Polynomial order of the approximation space.
   integer(kind = 4) :: ORDER

!> @brief Number of elements in each parametric direction.
   integer(kind = 4) :: SIZE

!> @brief Numbers of MPI processes in the three process-grid directions.
   integer(kind = 4) :: procx, procy, procz

!> @brief Amount of material drained over the simulation.
   real (kind = 8) :: drained = 0


contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads the oil-problem discretization and time parameters.
!>
!> @details
!> The expected leading argument list is:
!> `<size> <order> <procx> <procy> <procz> <steps> <dt>`.
!> Pump and drain coordinates are read later by \ref InitPumps.
!
!---------------------------------------------------------------------------
   subroutine InitializeParameters
      implicit none
      character(100) :: input
      integer(kind = 4) :: length
      integer(kind = 4) :: status

      ! ./l2 <size> <procx> <procy> <procz> <nsteps> <dt>
      ORDER = 2

      call GET_COMMAND_ARGUMENT(1, input, length, status)
      read(input, *) SIZE
      call GET_COMMAND_ARGUMENT(2, input, length, status)
      read(input, *) ORDER
      call GET_COMMAND_ARGUMENT(3, input, length, status)
      read(input, *) procx
      call GET_COMMAND_ARGUMENT(4, input, length, status)
      read(input, *) procy
      call GET_COMMAND_ARGUMENT(5, input, length, status)
      read(input, *) procz
      call GET_COMMAND_ARGUMENT(6, input, length, status)
      read(input, *) steps
      call GET_COMMAND_ARGUMENT(7, input, length, status)
      read(input, *) Dt

   end subroutine InitializeParameters

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reduces and prints the total drained quantity.
!>
!> @details
!> The process-local value \ref drained is summed over `MPI_COMM_WORLD`.
!> Rank zero writes the global value to standard output.
!
!---------------------------------------------------------------------------
   subroutine ComputeResults()
      use parallelism, ONLY: MYRANK
      use mpi
      implicit none
      real (kind = 8) :: fulldrained
      integer(kind = 4) :: ierr


      call MPI_Reduce(drained, fulldrained, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      if (MYRANK == 0) then
         write(*, *) fulldrained
      endif

   end subroutine ComputeResults

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads pump and drain coordinates from command-line arguments.
!>
!> @details
!> After the seven technical arguments consumed by
!> \ref InitializeParameters, the routine expects the number of pumps,
!> three coordinates per pump, the number of drains, and three coordinates
!> per drain. Coordinates are stored column-wise in \ref pumps and
!> \ref drains.
!
!---------------------------------------------------------------------------
   subroutine InitPumps()
      implicit none
      character(100) :: input
      integer(kind = 4) :: i, arg = 8 ! First argument after "technical" ones
      integer(kind = 4) :: length
      integer(kind = 4) :: status

      call GET_COMMAND_ARGUMENT(arg, input, length, status)
      read(input, *) npumps
      arg = arg + 1
      allocate(pumps(3, npumps))

      do i = 1, npumps
         call GET_COMMAND_ARGUMENT(arg, input, length, status)
         read(input, *) pumps(1, i)
         call GET_COMMAND_ARGUMENT(arg + 1, input, length, status)
         read(input, *) pumps(2, i)
         call GET_COMMAND_ARGUMENT(arg + 2, input, length, status)
         read(input, *) pumps(3, i)
         arg = arg + 3
      enddo

      call GET_COMMAND_ARGUMENT(arg, input, length, status)
      read(input, *) ndrains
      arg = arg + 1
      allocate(drains(3, ndrains))

      do i = 1, ndrains
         call GET_COMMAND_ARGUMENT(arg, input, length, status)
         read(input, *) drains(1, i)
         call GET_COMMAND_ARGUMENT(arg + 1, input, length, status)
         read(input, *) drains(2, i)
         call GET_COMMAND_ARGUMENT(arg + 2, input, length, status)
         read(input, *) drains(3, i)
         arg = arg + 3
      enddo

   end subroutine InitPumps

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Generates random channel curves and initializes pumps/drains.
!>
!> @details
!> The routine fills the global curve coordinate arrays with short random
!> walks and then calls \ref InitPumps to read external source/sink data.
!
!---------------------------------------------------------------------------
   subroutine InitInputData()
      implicit none
      integer(kind = 4) :: i, j
      real (kind = 8) :: t(3), x(3), dx(3), ddx(3)
      real (kind = 8) :: f = 0.1d0, step = 0.05d0

      do i = 0, cN - 1
         call random_number(x)
         cx(i * cL + 1) = x(1)
         cy(i * cL + 1) = x(2)
         cz(i * cL + 1) = x(3)
         call random_number(dx)

         do j = 2, cL
            call random_number(t)
            call random_number(ddx)
            ddx = 2 * (ddx - 0.5d0) - 0.4d0 * dx * sum(dx**2)
            dx = dx + 0.4d0 * ddx
            cx(i * cL + j) = cx(i * cL + j - 1) + step * dx(1)
            cy(i * cL + j) = cy(i * cL + j - 1) + step * dx(2)
            cz(i * cL + j) = cz(i * cL + j - 1) + step * dx(3)
         end do
      end do

      call InitPumps()

   end subroutine InitInputData

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes the squared distance from a point to a segment.
!>
!> @details
!> The point is orthogonally projected onto the segment whenever the
!> projection lies inside the segment interval. Otherwise the nearest
!> endpoint is used.
!
! Input:
! ------
!> @param[in] point
!> Point whose distance from the segment is evaluated.
!>
!> @param[in] ibeg
!> First segment endpoint.
!>
!> @param[in] iend
!> Second segment endpoint.
!
! Output:
! -------
!> @return d
!> Squared Euclidean distance from \p point to the segment.
!
!---------------------------------------------------------------------------
   function dist_from_segment(point, ibeg, iend) result (d)
      implicit none
      real (kind = 8), intent(in), dimension(3) :: point, ibeg, iend
      real (kind = 8) :: dx, dy, dz, cx, cy, cz, xx, yy, zz
      real (kind = 8) :: dot, len2, proj, d

      dx = iend(1) - ibeg(1)
      dy = iend(2) - ibeg(2)
      dz = iend(3) - ibeg(3)

      cx = point(1) - ibeg(1)
      cy = point(2) - ibeg(2)
      cz = point(3) - ibeg(3)

      dot = dx * cx + dy * cy + dz * cz
      len2 = dx ** 2 + dy **2 + dz ** 2
      proj = dot / len2

      if (proj < 0) then
         xx = ibeg(1)
         yy = ibeg(2)
         zz = ibeg(3)
      else if (proj > 1) then
         xx = iend(1)
         yy = iend(2)
         zz = iend(3)
      else
         xx = ibeg(1) + proj * dx
         yy = ibeg(2) + proj * dy
         zz = ibeg(3) + proj * dz
      end if

      d = (point(1) - xx)**2 + (point(2) - yy)**2 + (point(3) - zz)**2

   end function dist_from_segment

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes the squared distance from a point to a set of curves.
!>
!> @details
!> Each curve is represented as consecutive line segments. The returned
!> value is the minimum squared point-to-segment distance over all curve
!> segments.
!
! Input:
! ------
!> @param[in] point
!> Point whose distance from the curve set is evaluated.
!>
!> @param[in] cx
!> First coordinates of all curve nodes.
!>
!> @param[in] cy
!> Second coordinates of all curve nodes.
!>
!> @param[in] cz
!> Third coordinates of all curve nodes.
!>
!> @param[in] cN
!> Number of curves.
!>
!> @param[in] cL
!> Number of nodes stored per curve.
!
! Output:
! -------
!> @return fval
!> Minimum squared distance from \p point to the curve set.
!
!---------------------------------------------------------------------------
   function dist_from_curves(point, cx, cy, cz, cN, cL) result (fval)
      implicit none
      real (kind = 8), intent(in), dimension(3) :: point
      integer(kind = 4), intent(in) :: cN, cL
      real (kind = 8), intent(in) :: cx(cN * cL), cy(cN * cL), cz(cN * cL)
      real (kind = 8) :: ax, ay, az, bx, by, bz, fval
      real (kind = 8), dimension(3) :: p1, p2
      integer(kind = 4) :: i, j

      fval = 1e3
      do i = 0, cN - 1
         do j = 2, cL
            ax = cx(i * cL + j - 1)
            bx = cx(i * cL + j)
            ay = cy(i * cL + j - 1)
            by = cy(i * cL + j)
            az = cz(i * cL + j - 1)
            bz = cz(i * cL + j)

            p1 = (/ ax, ay, az/)
            p2 = (/ bx, by, bz/)

            fval = min(fval, dist_from_segment(point, p1, p2))
         end do
      end do

   end function dist_from_curves
!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the pump source contribution at a point.
!>
!> @details
!> The value is the sum of compactly supported radial falloff profiles
!> centered at all configured pump points.
!
! Input:
! ------
!> @param[in] x
!> First physical coordinate.
!>
!> @param[in] y
!> Second physical coordinate.
!>
!> @param[in] z
!> Third physical coordinate.
!
! Output:
! -------
!> @return fval
!> Pumping source value at the supplied point.
!
!---------------------------------------------------------------------------
   function pumping(x, y, z) result (fval)
      use math, ONLY: falloff, norm_2
      implicit none
      real (kind = 8) :: x, y, z
      real (kind = 8) :: fval
      integer(kind = 4) :: i
      real (kind = 8), dimension(3) :: p1


      fval = 0.d0
      do i = 1, npumps
         p1 = (/ x, y, z/)
         fval = fval + pumping_strength * falloff(0.d0, radius, norm_2(pumps(:, i) - p1))
      enddo

   end function pumping
!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the drain sink contribution at a point.
!>
!> @details
!> The value is the sum of compactly supported radial falloff profiles
!> centered at all drain points, multiplied by the current solution value
!> \p u.
!
! Input:
! ------
!> @param[in] u
!> Local solution value multiplying the drain strength.
!>
!> @param[in] x
!> First physical coordinate.
!>
!> @param[in] y
!> Second physical coordinate.
!>
!> @param[in] z
!> Third physical coordinate.
!
! Output:
! -------
!> @return fval
!> Draining sink value at the supplied point.
!
!---------------------------------------------------------------------------
   function draining(u, x, y, z) result (fval)
      use math, ONLY: falloff, norm_2
      implicit none
      real (kind = 8) :: u, x, y, z
      real (kind = 8) :: fval
      integer :: i
      real (kind = 8), dimension(3) :: p1

      fval = 0.d0
      do i = 1, ndrains
         p1 = (/ x, y, z/)
         fval = fval + draining_strength * falloff(0.d0, radius, norm_2(drains(:, i) - p1))
      enddo
      fval = fval * u

   end function draining

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the spatial permeability coefficient.
!>
!> @details
!> The coefficient interpolates between \ref Kqmin and \ref Kqmax according
!> to the distance from the generated channel curves.
!
! Input:
! ------
!> @param[in] x
!> First physical coordinate.
!>
!> @param[in] y
!> Second physical coordinate.
!>
!> @param[in] z
!> Third physical coordinate.
!
! Output:
! -------
!> @return val
!> Permeability value at the supplied point.
!
!---------------------------------------------------------------------------
   function kq(x, y, z) result (val)
      use math, ONLY: falloff, lerp
      implicit none
      real (kind = 8), intent(in) :: x, y, z
      real (kind = 8) :: val, dist
      real (kind = 8), dimension(3) :: p1


      p1 = (/ x, y, z/)

      dist = sqrt(dist_from_curves(p1, cx, cy, cz, cN, cL))
      val = lerp(falloff(0.d0, 0.06d0, dist), Kqmin, Kqmax)

   end function kq

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the nonlinear mobility factor.
!>
!> @details
!> The current oil model uses an exponential response \f$\exp(\mu u)\f$.
!
! Input:
! ------
!> @param[in] x
!> First physical coordinate, retained for interface symmetry.
!>
!> @param[in] y
!> Second physical coordinate, retained for interface symmetry.
!>
!> @param[in] z
!> Third physical coordinate, retained for interface symmetry.
!>
!> @param[in] u
!> Local solution value.
!
! Output:
! -------
!> @return val
!> Nonlinear mobility value.
!
!---------------------------------------------------------------------------
   function bq(x, y, z, u) result (val)
      real (kind = 8) :: x, y, z, u
      real (kind = 8) :: val

      val = exp(mi * u)
 
   end function bq
!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the initial oil concentration profile.
!>
!> @details
!> The profile is localized by a compact-support bump and modulated by the
!> distance from the generated channel curves.
!
! Input:
! ------
!> @param[in] x
!> First physical coordinate.
!>
!> @param[in] y
!> Second physical coordinate.
!>
!> @param[in] z
!> Third physical coordinate.
!
! Output:
! -------
!> @return val
!> Initial concentration value.
!
!---------------------------------------------------------------------------
   function initial_state(x, y, z) result (val)
      use math, ONLY: falloff, bump3d, lerp
      implicit none
      real (kind = 8), intent(in) :: x, y, z
      real (kind = 8) :: dist, val
      real (kind = 8), dimension(3) :: p1


      p1 = (/ x, y, z/)
      dist = sqrt(dist_from_curves(p1, cx, cy, cz, cN, cL))
      val = 0.1d0 * lerp(falloff(0.d0, 0.1d0, dist), 0.d0, 1.d0) * bump3d(0.2d0, 0.6d0, x, y, z)

   end function initial_state


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the pointwise source term for the current ADS API.
!>
!> @details
!> At the initial pseudo-step the callback returns the initial state. For
!> later time steps it returns pump injection minus the nonnegative drain
!> contribution.
!
! Input:
! ------
!> @param[in] un
!> Previous solution value at the quadrature point.
!>
!> @param[in] du
!> Previous solution gradient at the quadrature point.
!>
!> @param[in] X
!> Physical coordinates of the quadrature point.
!
! Output:
! -------
!> @return ret
!> Pointwise source value.
!
!---------------------------------------------------------------------------
   function forcing(un, du, X) result(ret)
      implicit none
      real(kind = 8), intent(in) :: un
      real(kind = 8), intent(in), dimension(3) :: du
      real(kind = 8), intent(in), dimension(3) :: X
      real(kind = 8) :: ret

      if (t > 0.d0) then
         ret = pumping(X(1), X(2), X(3)) - &
               max(0.d0, draining(un, X(1), X(2), X(3)))
      else
         ret = initial_state(X(1), X(2), X(3))
      endif

   end function forcing


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Samples the permeability coefficient at local quadrature points.
!>
!> @details
!> The routine reconstructs one-dimensional quadrature data for the three
!> spline directions and fills \p Kq_vals for the locally owned element
!> range.
!
! Input:
! ------
!> @param[in] Ux
!> Knot vector in the first direction.
!>
!> @param[in] px
!> Polynomial degree in the first direction.
!>
!> @param[in] nx
!> Number of basis functions minus one in the first direction.
!>
!> @param[in] minex
!> First local element in the first direction.
!>
!> @param[in] maxex
!> Last local element in the first direction.
!>
!> @param[in] nelemx
!> Number of elements in the first direction.
!>
!> @param[in] Uy
!> Knot vector in the second direction.
!>
!> @param[in] py
!> Polynomial degree in the second direction.
!>
!> @param[in] ny
!> Number of basis functions minus one in the second direction.
!>
!> @param[in] miney
!> First local element in the second direction.
!>
!> @param[in] maxey
!> Last local element in the second direction.
!>
!> @param[in] nelemy
!> Number of elements in the second direction.
!>
!> @param[in] Uz
!> Knot vector in the third direction.
!>
!> @param[in] pz
!> Polynomial degree in the third direction.
!>
!> @param[in] nz
!> Number of basis functions minus one in the third direction.
!>
!> @param[in] minez
!> First local element in the third direction.
!>
!> @param[in] maxez
!> Last local element in the third direction.
!>
!> @param[in] nelemz
!> Number of elements in the third direction.
!
! Output:
! -------
!> @param[out] Kq_vals
!> Permeability values sampled at local quadrature points.
!
!---------------------------------------------------------------------------
   subroutine CacheKqValues(Ux, px, nx, minex, maxex, nelemx, Uy, py, ny, miney, maxey, nelemy, Uz, pz, nz, minez, &
      maxez, nelemz, Kq_vals)
      use basis, ONLY: BasisData
      implicit none
      integer(kind = 4), intent(in) :: nx, px, minex, maxex, nelemx
      integer(kind = 4), intent(in) :: ny, py, miney, maxey, nelemy
      integer(kind = 4), intent(in) :: nz, pz, minez, maxez, nelemz
      real (kind = 8), intent(in) :: Ux(0:nx + px + 1)
      real (kind = 8), intent(in) :: Uy(0:ny + py + 1)
      real (kind = 8), intent(in) :: Uz(0:nz + pz + 1)
      real (kind = 8), intent(out) :: Kq_vals(px + 1, py + 1, pz + 1, maxex - minex + 1, maxey - miney + 1, maxez - minez + 1)
      integer(kind = 4) :: mx, my, mz, ngx, ngy, ngz, ex, ey, ez
      integer(kind = 4) :: kx, ky, kz, ax, ay, az, d
      integer(kind = 4) :: Ox(nelemx), Oy(nelemy), Oz(nelemz)
      real (kind = 8) :: Jx(nelemx), Jy(nelemy), Jz(nelemz)
      real (kind = 8) :: Wx(px + 1), Wy(py + 1), Wz(pz + 1)
      real (kind = 8) :: Xx(px + 1, nelemx)
      real (kind = 8) :: Xy(py + 1, nelemy)
      real (kind = 8) :: Xz(pz + 1, nelemz)
      real (kind = 8) :: NNx(0:2, 0:px, px + 1, nelemx), NNy(0:2, 0:py, py + 1, nelemy), NNz(0:2, 0:pz, pz + 1, nelemz)

      d = 2
      mx = nx + px + 1
      ngx = px + 1
      my = ny + py + 1
      ngy = py + 1
      mz = nz + pz + 1
      ngz = pz + 1

      call BasisData(px, mx, Ux, d, ngx, nelemx, Ox, Jx, Wx, Xx, NNx)
      call BasisData(py, my, Uy, d, ngy, nelemy, Oy, Jy, Wy, Xy, NNy)
      call BasisData(pz, mz, Uz, d, ngz, nelemz, Oz, Jz, Wz, Xz, NNz)

      do ex = minex, maxex
         do ey = miney, maxey
            do ez = minez, maxez
               do kx = 1, ngx
                  do ky = 1, ngy
                     do kz = 1, ngz
                        Kq_vals(kx, ky, kz, ex - minex + 1, ey - miney + 1, ez - minez + 1) = kq(Xx(kx, ex), Xy(ky, ey), Xz(kz, ez))
                     end do
                  end do
               end do
            end do
         end do
      end do

   end subroutine CacheKqValues
!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Fills the module-level permeability cache for an ADS setup.
!>
!> @details
!> The routine extracts knot vectors, polynomial degrees, basis sizes, and
!> local element ranges from \p ads, then delegates the sampling work to
!> \ref CacheKqValues.
!
! Input:
! ------
!> @param[in] ads
!> ADS setup structure defining the local spline space.
!
!---------------------------------------------------------------------------
   subroutine PrecomputeKq(ads)
      use Setup, ONLY: ADS_Setup
      implicit none
      type (ADS_setup) :: ads

      call CacheKqValues &
      (ads % Ux, ads % p(1), ads % n(1), ads % mine(1), ads % maxe(1), ads % nelem(1), &
      ads % Uy, ads % p(2), ads % n(2), ads % mine(2), ads % maxe(2), ads % nelem(2), &
      ads % Uz, ads % p(3), ads % n(3), ads % mine(3), ads % maxe(3), ads % nelem(3), &
      Kqvals)

   end subroutine PrecomputeKq


end module input_data
