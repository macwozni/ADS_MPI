module input_data

   implicit none

   ! Curve number and lentgh
   integer(kind = 4), parameter :: cN = 30, cL = 16
   real (kind = 8), parameter :: radius = 0.15, pumping_strength = 1, draining_strength = 1

   !!! krzywe dane przez segmenty
   real (kind = 8) :: cx(cN * cL), cy(cN * cL), cz(cN * cL)
   real (kind = 8) :: mi = 10.d0
   real (kind = 8) :: GROUND = 0.2
   real (kind = 8), parameter :: Kqmin = 1.d0, Kqmax = 1000.d0
   integer(kind = 4) :: npumps, ndrains
   real (kind = 8), allocatable, dimension(:,:) :: pumps, drains

   ! Buffer for values of permeability function
   real (kind = 8), allocatable :: Kqvals(:,:,:,:,:,:)

   ! Current time
   real (kind = 8) :: t

   ! Time and timestep
   real (kind = 8) :: Dt

   ! Number of iterations
   integer :: steps

   ! Statistics computed during the simulation
   real (kind = 8) :: pollution = 0

   ! order of approximations
   integer(kind = 4) :: ORDER

   ! number of elements in one dimension
   integer(kind = 4) :: SIZE

   integer(kind = 4) :: procx, procy, procz

   real (kind = 8) :: drained = 0


contains


   ! -------------------------------------------------------------------
   ! Sets values of parameters (order and size)
   ! -------------------------------------------------------------------
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





   ! forcing
   ! x, y, z - point in space
   function forcing(x, y, z) result (fval)
      implicit none
      real (kind = 8) :: x, y, z
      real (kind = 8) :: fval

   end function forcing



   ! 
   ! x, y, z - point in space
   function g(x, y, z) result (fval)
      implicit none
      real (kind = 8) :: x, y, z
      real (kind = 8) :: fval

   end function g




   ! 
   ! x, y, z - point in space
   function b(x, y, z) result (fval)
      implicit none
      real (kind = 8) :: x, y, z
      real (kind = 8) :: fval

   end function b




   ! normal to boundary
   ! x, y, z - point in space
   function n(x, y, z) result (fval)
      implicit none
      real (kind = 8) :: x, y, z
      real (kind = 8) :: fval

   end function n


end module input_data
