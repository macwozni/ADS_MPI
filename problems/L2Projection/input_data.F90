module input_data

   implicit none

   ! order of approximations
   integer(kind = 4) :: order

   ! number of elements in one dimension
   integer(kind = 4) :: isizex, isizey, isizez

   integer(kind = 4) :: procx, procy, procz


contains


   ! -------------------------------------------------------------------
   ! Sets values of parameters (order and size)
   ! -------------------------------------------------------------------
   subroutine InitializeParameters
      implicit none
      character(100) :: input
      integer(kind = 4) :: length
      integer(kind = 4) :: status

      ! ./l2 <size> <order> <procx> <procy> <procz>

      call GET_COMMAND_ARGUMENT(1, input, length, status)
      read(input, *) isizex
      call GET_COMMAND_ARGUMENT(2, input, length, status)
      read(input, *) isizey
      call GET_COMMAND_ARGUMENT(3, input, length, status)
      read(input, *) isizez
      call GET_COMMAND_ARGUMENT(4, input, length, status)
      read(input, *) order
      call GET_COMMAND_ARGUMENT(5, input, length, status)
      read(input, *) procx
      call GET_COMMAND_ARGUMENT(6, input, length, status)
      read(input, *) procy
      call GET_COMMAND_ARGUMENT(7, input, length, status)
      read(input, *) procz
      
   end subroutine



end module
