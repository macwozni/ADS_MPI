module input_data

   implicit none

   ! order of approximations
   integer(kind = 4) :: ORDER

   ! number of elements in one dimension
   integer(kind = 4) :: SIZE

   integer(kind = 4) :: procx, procy, procz

   real (kind = 8) :: l2norm, fullnorm


contains


   ! -------------------------------------------------------------------
   ! Sets values of parameters (order and size)
   ! -------------------------------------------------------------------
   subroutine InitializeParameters
      implicit none
      character(100) :: input
      integer(kind = 4) :: length
      integer(kind = 4) :: status

      ! ./l2 <size> <procx> <procy> <procz>
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

   end subroutine



end module