program solution_reconstruction_probe
   use ISO_FORTRAN_ENV, only: ERROR_UNIT
   use Setup, only: ADS_Setup, ADS_compute_data
   use solution_reconstruction, only: FormUn
   implicit none

   character(len=32) :: mode
   integer(kind=4) :: selector
   type(ADS_Setup) :: ads
   type(ADS_compute_data) :: ads_data

   if (command_argument_count() /= 1) stop 90
   call get_command_argument(1, mode)

   select case (trim(mode))
   case ("zero")
      selector = 0
   case ("four")
      selector = 4
   case default
      write(ERROR_UNIT, '(A)') "unknown probe mode"
      stop 90
   end select

   call FormUn(selector, ads, ads_data)

   write(ERROR_UNIT, '(A)') "validation unexpectedly succeeded"
   stop 98

end program solution_reconstruction_probe
