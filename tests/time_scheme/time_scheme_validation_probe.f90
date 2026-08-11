program time_scheme_validation_probe
   use Setup, only: ADS_Setup
   use time_scheme, only: ValidateSpaces
   implicit none

   type(ADS_Setup) :: ads, ads_test, ads_trial
   character(len=16) :: mode, field, token
   integer(kind=4) :: axis, value

   if (command_argument_count() /= 4) stop 64
   call get_command_argument(1, mode)
   call get_command_argument(2, token)
   read (token, *, err=900) axis
   call get_command_argument(3, field)
   call get_command_argument(4, token)
   read (token, *, err=900) value
   if (axis < 1 .or. axis > 3) stop 64

   ads%p = 1
   ads%ng = 2
   ads_test%p = 2
   ads_test%ng = 3
   ads_trial%p = 1
   ads_trial%ng = 2

   select case (trim(mode))
   case ('standard')
      call set_value(ads, axis, field, value)
      call ValidateSpaces(ads)
   case ('test')
      call set_value(ads_test, axis, field, value)
      call ValidateSpaces(ads_test, ads_trial)
   case ('trial')
      call set_value(ads_trial, axis, field, value)
      call ValidateSpaces(ads_test, ads_trial)
   case default
      stop 64
   end select

   stop 0

900 stop 64

contains

   subroutine set_value(space, selected_axis, selected_field, selected_value)
      type(ADS_Setup), intent(inout) :: space
      integer(kind=4), intent(in) :: selected_axis, selected_value
      character(len=*), intent(in) :: selected_field

      select case (trim(selected_field))
      case ('p')
         space%p(selected_axis) = selected_value
      case ('ng')
         space%ng(selected_axis) = selected_value
      case default
         stop 64
      end select
   end subroutine set_value

end program time_scheme_validation_probe
