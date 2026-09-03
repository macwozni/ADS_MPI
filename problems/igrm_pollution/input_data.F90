!------------------------------------------------------------------------------
!> @file input_data.F90
!> @brief Parameters and physical data for the corrected 3D pollution DPG case.
!------------------------------------------------------------------------------
module input_data

   use argument_parser, only: ReadIntegerArgument, RequirePositiveInteger, &
                              RequireSafeSplineDimensions

   implicit none
   private

   real(kind=8), parameter, public :: DOMAIN_MIN = 0.d0
   real(kind=8), parameter, public :: DOMAIN_LENGTH = 5000.d0
   real(kind=8), parameter, public :: DOMAIN_MAX = DOMAIN_MIN + DOMAIN_LENGTH
   real(kind=8), parameter, public :: TIME_STEP = 1.8d0
   real(kind=8), parameter, public :: DT = TIME_STEP
   real(kind=8), parameter, public :: DIFFUSION(3) = &
      (/ 50.d0, 50.d0, 0.5d0 /)
   real(kind=8), parameter, public :: WIND_SPEED = 5.d0
   real(kind=8), parameter, public :: SOURCE_CENTER(3) = &
      (/ 3000.d0, 2000.d0, 2000.d0 /)
   real(kind=8), parameter, public :: SOURCE_RADIUS = 25.d0
   integer(kind=4), parameter, public :: MAX_POLYNOMIAL_DEGREE = 9
   integer(kind=4), parameter, public :: DEFAULT_OUTPUT_RESOLUTION = 100

   integer(kind=4), public :: nelem
   logical, public :: adapt_mesh
   integer(kind=4), public :: p_trial, c_trial
   integer(kind=4), public :: p_test, c_test
   integer(kind=4), public :: steps
   integer(kind=4), public :: procx, procy, procz
   integer(kind=4), public :: output_resolution = DEFAULT_OUTPUT_RESOLUTION

   public :: InitializeParameters
   public :: SplineDofs
   public :: adapted_coordinate
   public :: emission
   public :: wind_at_time

contains

!------------------------------------------------------------------------------
!> Read the upstream seven problem arguments followed by the MPI process grid.
!------------------------------------------------------------------------------
subroutine InitializeParameters
   implicit none

   integer(kind=4) :: adapt_value
   integer(kind=8) :: test_dofs, trial_dofs

   if (command_argument_count() /= 10) then
      write (*, '(A)') 'proper usage with arguments: ' // &
         '<N> <adapt:0|1> <p_trial> <C_trial> <p_test> <C_test> ' // &
         '<steps> <procx> <procy> <procz>'
      stop 5
   end if

   call ReadIntegerArgument(1, nelem)
   call ReadIntegerArgument(2, adapt_value)
   call ReadIntegerArgument(3, p_trial)
   call ReadIntegerArgument(4, c_trial)
   call ReadIntegerArgument(5, p_test)
   call ReadIntegerArgument(6, c_test)
   call ReadIntegerArgument(7, steps)
   call ReadIntegerArgument(8, procx)
   call ReadIntegerArgument(9, procy)
   call ReadIntegerArgument(10, procz)

   call RequirePositiveInteger(nelem, 'number of elements')
   if (adapt_value /= 0 .and. adapt_value /= 1) then
      write (*, '(A)') 'adapt flag must be zero or one'
      stop 5
   end if
   adapt_mesh = adapt_value == 1

   call RequirePositiveInteger(p_trial, 'trial polynomial degree')
   call RequirePositiveInteger(p_test, 'test polynomial degree')
   if (p_trial > MAX_POLYNOMIAL_DEGREE) then
      write (*, '(A,I0)') &
         'unsupported polynomial degree in iGRM trial space; ' // &
         'maximum supported degree is ', MAX_POLYNOMIAL_DEGREE
      stop 5
   end if
   if (p_test > MAX_POLYNOMIAL_DEGREE) then
      write (*, '(A,I0)') &
         'unsupported polynomial degree in iGRM test space; ' // &
         'maximum supported degree is ', MAX_POLYNOMIAL_DEGREE
      stop 5
   end if

   call ValidateContinuity(c_trial, p_trial, 'trial')
   call ValidateContinuity(c_test, p_test, 'test')
   call RequirePositiveInteger(steps, 'number of time steps')
   call RequirePositiveInteger(procx, 'process-grid dimension')
   call RequirePositiveInteger(procy, 'process-grid dimension')
   call RequirePositiveInteger(procz, 'process-grid dimension')

   call RequireSafeSplineDimensions(nelem, p_trial)
   call RequireSafeSplineDimensions(nelem, p_test)
   trial_dofs = SplineDofs(nelem, p_trial, c_trial)
   test_dofs = SplineDofs(nelem, p_test, c_test)
   call ValidateSplineExtent(trial_dofs, p_trial)
   call ValidateSplineExtent(test_dofs, p_test)
   if (test_dofs < trial_dofs) then
      write (*, '(A)') &
         'test-space dimension must be greater than or equal to trial-space dimension'
      stop 5
   end if

   call ReadOutputResolution

end subroutine InitializeParameters

!------------------------------------------------------------------------------
!> Number of basis functions for N elements, degree p and continuity C.
!------------------------------------------------------------------------------
pure function SplineDofs(element_count, degree, continuity) result(dofs)
   implicit none

   integer(kind=4), intent(in) :: element_count, degree, continuity
   integer(kind=8) :: dofs

   dofs = int(degree, kind=8) + 1_8 + &
      (int(element_count, kind=8) - 1_8)* &
      int(degree - continuity, kind=8)

end function SplineDofs

!------------------------------------------------------------------------------
!> Map a uniform unit coordinate to the upstream x-adapted physical mesh.
!------------------------------------------------------------------------------
pure function adapted_coordinate(unit_coordinate) result(coordinate)
   implicit none

   real(kind=8), intent(in) :: unit_coordinate
   real(kind=8) :: coordinate, mapped

   if (unit_coordinate < 0.5d0) then
      mapped = unit_coordinate/0.5d0*0.99d0
   else
      mapped = (unit_coordinate - 0.5d0)/0.5d0*0.01d0 + 0.99d0
   end if
   coordinate = DOMAIN_MIN + DOMAIN_LENGTH*mapped

end function adapted_coordinate

!------------------------------------------------------------------------------
!> Compact polynomial source centered at (3000,2000,2000), radius 25.
!------------------------------------------------------------------------------
pure function emission(X) result(value)
   implicit none

   real(kind=8), intent(in) :: X(3)
   real(kind=8) :: value
   real(kind=8) :: normalized(3), radius_squared

   normalized = (X - SOURCE_CENTER)/SOURCE_RADIUS
   radius_squared = min(sum(normalized*normalized), 1.d0)
   value = (radius_squared - 1.d0)**2*(radius_squared + 1.d0)**2

end function emission

!------------------------------------------------------------------------------
!> Corrected wind law used uniformly from t=0 onward.
!------------------------------------------------------------------------------
pure function wind_at_time(time) result(wind)
   implicit none

   real(kind=8), intent(in) :: time
   real(kind=8) :: wind(3)
   real(kind=8), parameter :: PI = &
      3.141592653589793238462643383279502884197d0
   real(kind=8) :: angle, phase, scaled_time

   scaled_time = time/150.d0
   phase = sin(scaled_time) + 0.5d0*sin(2.3d0*scaled_time)
   angle = PI/3.d0*phase + 3.d0*PI/8.d0
   wind(1) = WIND_SPEED*cos(angle)
   wind(2) = WIND_SPEED*sin(angle)
   wind(3) = WIND_SPEED/3.d0*sin(scaled_time)

end function wind_at_time

subroutine ValidateContinuity(continuity, degree, space_name)
   implicit none

   integer(kind=4), intent(in) :: continuity, degree
   character(*), intent(in) :: space_name

   if (continuity < 0 .or. continuity > degree - 1) then
      write (*, '(A)') trim(space_name) // &
         ' continuity must be between zero and polynomial degree minus one'
      stop 5
   end if

end subroutine ValidateContinuity

subroutine ValidateSplineExtent(dofs, degree)
   implicit none

   integer(kind=8), intent(in) :: dofs
   integer(kind=4), intent(in) :: degree
   integer(kind=8) :: max_default_integer

   max_default_integer = int(huge(degree), kind=8)
   if (dofs <= 0_8 .or. &
       dofs > max_default_integer - int(degree, kind=8) - 1_8) then
      write (*, '(A)') 'spline dimensions exceed supported integer range'
      stop 5
   end if

end subroutine ValidateSplineExtent

subroutine ReadOutputResolution
   implicit none

   character(len=128) :: raw_value
   integer(kind=4) :: env_length, env_status, read_status, value

   output_resolution = DEFAULT_OUTPUT_RESOLUTION
   raw_value = ''
   call get_environment_variable('ADS_POLLUTION_OUTPUT_RESOLUTION', &
                                 raw_value, env_length, env_status, &
                                 trim_name=.true.)
   if (env_status == 1) return
   if (env_status /= 0 .or. env_length <= 0 .or. &
       env_length > len(raw_value)) then
      call InvalidOutputResolution
   end if
   if (.not. IsIntegerToken(raw_value(1:env_length))) then
      call InvalidOutputResolution
   end if
   read (raw_value(1:env_length), *, iostat=read_status) value
   if (read_status /= 0 .or. value <= 0) call InvalidOutputResolution
   output_resolution = value

end subroutine ReadOutputResolution

logical function IsIntegerToken(input) result(valid)
   implicit none

   character(*), intent(in) :: input
   character(len=len(input)) :: token
   integer(kind=4) :: index, token_length

   token = adjustl(input)
   token_length = len_trim(token)
   valid = .false.
   if (token_length == 0) return

   index = 1
   if (token(index:index) == '+' .or. token(index:index) == '-') &
      index = index + 1
   if (index > token_length) return
   do while (index <= token_length)
      if (token(index:index) < '0' .or. token(index:index) > '9') return
      index = index + 1
   end do
   valid = .true.

end function IsIntegerToken

subroutine InvalidOutputResolution
   implicit none

   write (*, '(A)') &
      'ADS_POLLUTION_OUTPUT_RESOLUTION must be a positive integer'
   stop 5

end subroutine InvalidOutputResolution

end module input_data
