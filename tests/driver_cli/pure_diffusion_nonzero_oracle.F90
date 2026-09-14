! Test-only manufactured solution for the real MPI/iGRM diffusion path.
!
! The exactly representable steady manufactured solution
!
!   u(x,y,z) = q(x)*q(y)*q(z),  q(s) = s^2*(3-2*s),
!
! has a zero normal derivative on every face.  Its source f=-Laplace(u)
! exercises all three directional diffusion operators.  A mass-only iGRM
! projection creates the initial coefficients without touching private ADS
! buffers; a correct diffusion step must then preserve this steady state.
module pure_diffusion_oracle_rhs

   implicit none

contains

   pure function q(s) result(value)
      implicit none
      real(kind=8), intent(in) :: s
      real(kind=8) :: value

      value = s*s*(3.0d0 - 2.0d0*s)

   end function q

   pure function exact_field(x) result(value)
      implicit none
      real(kind=8), dimension(3), intent(in) :: x
      real(kind=8) :: value

      value = q(x(1))*q(x(2))*q(x(3))

   end function exact_field

   function manufactured_source(un, du, x) result(value)
      implicit none
      real(kind=8), intent(in) :: un
      real(kind=8), dimension(3), intent(in) :: du, x
      real(kind=8) :: value

      value = (12.0d0*x(1) - 6.0d0)*q(x(2))*q(x(3)) + &
              (12.0d0*x(2) - 6.0d0)*q(x(1))*q(x(3)) + &
              (12.0d0*x(3) - 6.0d0)*q(x(1))*q(x(2))

   end function manufactured_source

   subroutine initial_projection_rhs( &
      ads, x, k, e, a, du, n, un11, un13, un23, ads_data, j, w, &
      direction, substep, alpha_step, forcing_cb, value)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use Interfaces, ONLY: forcing_fun
      implicit none
      type(ADS_Setup), intent(in) :: ads
      real(kind=8), dimension(3), intent(in) :: x, du
      integer(kind=4), dimension(3), intent(in) :: k, e, a, direction
      integer(kind=4), intent(in) :: n, substep
      real(kind=8), intent(in) :: un11, un13, un23, j, w
      type(ADS_compute_data), intent(in) :: ads_data
      real(kind=8), dimension(7, 3), intent(in) :: alpha_step
      procedure(forcing_fun) :: forcing_cb
      real(kind=8), intent(out) :: value
      real(kind=8) :: test_value

      test_value = ads%NNx(0, a(1), k(1), e(1))* &
                   ads%NNy(0, a(2), k(2), e(2))* &
                   ads%NNz(0, a(3), k(3), e(3))
      value = j*w*test_value*exact_field(x)

   end subroutine initial_projection_rhs

end module pure_diffusion_oracle_rhs

program pure_diffusion_nonzero_oracle

   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANK, InitializeParallelism, Cleanup_Parallelism, &
                          AbortOnError
   use communicators, ONLY: CreateCommunicators, Cleanup_Communicators
   use time_scheme, ONLY: BackwardEuler3DStep, &
                          ConfigureBackwardEuler3DTimeScheme, &
                          ConfigureDouglasGunn3DTimeScheme, &
                          ConfigureMassOnly3DTimeScheme, &
                          ConfigurePeacemanRachford3DTimeScheme, &
                          DouglasGunn3DStep, PeacemanRachford3DStep, &
                          TimeScheme3D, ValidateSpaces
   use ADSS, ONLY: Initialize, Cleanup_ADS, Cleanup_data, PrintSolution
   use pure_diffusion_oracle_rhs, ONLY: exact_field, &
                                        initial_projection_rhs, &
                                        manufactured_source
   use, intrinsic :: ieee_arithmetic, ONLY: ieee_is_finite

   implicit none

   integer(kind=4), dimension(3), parameter :: NELEM = (/3, 3, 3/)
   integer(kind=4), dimension(3), parameter :: P_TEST = (/4, 4, 4/)
   integer(kind=4), dimension(3), parameter :: P_TRIAL = (/3, 3, 3/)

   type(ADS_Setup) :: ads_test, ads_trial
   type(ADS_compute_data) :: ads_data
   type(TimeScheme3D) :: initial_scheme, scheme
   character(len=32) :: argument, scheme_name
   real(kind=8) :: dt
   integer(kind=4) :: ierr, step, steps

   if (command_argument_count() /= 3) then
      write(*, '(A)') &
         'usage: pure_diffusion_nonzero_oracle <dg|pr|be> <dt> <steps>'
      stop 5
   end if
   call get_command_argument(1, scheme_name)
   call get_command_argument(2, argument)
   read(argument, *, iostat=ierr) dt
   if (ierr /= 0) then
      write(*, '(A,A)') 'invalid time step: ', trim(argument)
      stop 5
   end if
   if (.not. ieee_is_finite(dt) .or. dt <= 0.0d0) then
      write(*, '(A,A)') 'invalid time step: ', trim(argument)
      stop 5
   end if
   call get_command_argument(3, argument)
   read(argument, *, iostat=ierr) steps
   if (ierr /= 0) then
      write(*, '(A,A)') 'invalid step count: ', trim(argument)
      stop 5
   end if
   if (steps <= 0) then
      write(*, '(A,A)') 'invalid step count: ', trim(argument)
      stop 5
   end if

   call InitializeParallelism(2, 1, 1, ierr)
   call AbortOnError(ierr, 'manufactured diffusion parallel initialization')
   call CreateCommunicators(ierr)
   call AbortOnError(ierr, 'manufactured diffusion communicators')
   call Initialize(NELEM, P_TEST, P_TRIAL, P_TRIAL - 1, &
                   ads_test, ads_trial, ads_data, ierr)
   call AbortOnError(ierr, 'manufactured diffusion initialization')
   call ValidateSpaces(ads_test, ads_trial)

   call ConfigureMassOnly3DTimeScheme(initial_scheme)
   select case (trim(scheme_name))
   case ('dg')
      call ConfigureDouglasGunn3DTimeScheme(dt, scheme, &
                                             include_transport=.false.)
   case ('pr')
      call ConfigurePeacemanRachford3DTimeScheme(dt, scheme, &
                                                  include_transport=.false.)
   case ('be')
      call ConfigureBackwardEuler3DTimeScheme(dt, scheme, &
                                               include_transport=.false.)
   case default
      if (MYRANK == 0) write(*, '(A,A)') 'unknown time scheme: ', &
                                             trim(scheme_name)
      stop 5
   end select

   ads_test%tau = 1.0d0
   ads_trial%tau = 1.0d0
   call DouglasGunn3DStep(initial_scheme, 0, manufactured_source, &
                          ads_test, ads_trial, ads_data, 1, ierr, &
                          initial_projection_rhs)
   call AbortOnError(ierr, 'manufactured diffusion initial projection')
   call ReportError(0, ads_trial, ads_data)
   call PrintSolution(0, ads_trial, ads_data%FF)

   ads_test%tau = dt
   ads_trial%tau = dt
   do step = 1, steps
      call RunStep(step, scheme_name, scheme, ads_test, ads_trial, &
                   ads_data, ierr)
      call AbortOnError(ierr, 'manufactured diffusion step')
      call ReportError(step, ads_trial, ads_data)
      call PrintSolution(step, ads_trial, ads_data%FF)
   end do

   call Cleanup_ADS(ads_test, ierr)
   call Cleanup_ADS(ads_trial, ierr)
   call Cleanup_data(ads_data, ierr)
   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

contains

   subroutine RunStep(step, name, configured_scheme, test_space, trial_space, &
                      data, status)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use time_scheme, ONLY: BackwardEuler3DStep, DouglasGunn3DStep, &
                             PeacemanRachford3DStep, TimeScheme3D
      implicit none
      integer(kind=4), intent(in) :: step
      character(len=*), intent(in) :: name
      type(TimeScheme3D), intent(in) :: configured_scheme
      type(ADS_Setup), intent(in) :: test_space, trial_space
      type(ADS_compute_data), intent(inout) :: data
      integer(kind=4), intent(out) :: status

      select case (trim(name))
      case ('dg')
         call DouglasGunn3DStep(configured_scheme, step, manufactured_source, &
                                test_space, trial_space, data, 1, status)
      case ('pr')
         call PeacemanRachford3DStep(configured_scheme, step, manufactured_source, &
                                     test_space, trial_space, data, 1, status)
      case ('be')
         call BackwardEuler3DStep(configured_scheme, step, manufactured_source, &
                                  test_space, trial_space, data, 1, status)
      end select

   end subroutine RunStep

   subroutine ReportError(step, trial_space, data)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use utils, ONLY: NormL2
      implicit none
      integer(kind=4), intent(in) :: step
      type(ADS_Setup), intent(in) :: trial_space
      type(ADS_compute_data), intent(in) :: data
      real(kind=8) :: error

      call NormL2(trial_space, data%FF, error, exact_field)
      if (MYRANK == 0) then
         write(*, '(A,I0,A,ES24.16)') 'manufactured L2 error step ', step, &
                                      ': ', error
      end if

   end subroutine ReportError

end program pure_diffusion_nonzero_oracle
