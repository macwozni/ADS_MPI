program int2str_probe
   use utils, only: int2str
   implicit none

   character(len=32) :: integer_arg, length_arg
   character(len=:), allocatable :: value
   integer(kind=4) :: number, target_length
   integer(kind=4) :: read_status

   if (command_argument_count() /= 2) stop 90

   call get_command_argument(1, integer_arg)
   call get_command_argument(2, length_arg)

   read (integer_arg, *, iostat=read_status) number
   if (read_status /= 0) stop 91

   read (length_arg, *, iostat=read_status) target_length
   if (read_status /= 0 .or. target_length < 0) stop 92

   allocate(character(len=target_length) :: value)
   call int2str(number, value)

   write (*, '("[",A,"]")') value

end program int2str_probe
