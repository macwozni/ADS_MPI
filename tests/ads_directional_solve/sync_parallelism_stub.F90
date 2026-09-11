module parallelism
   implicit none

   integer(kind=4) :: MYRANK = 0
   integer(kind=4) :: MYRANKX = 0
   integer(kind=4) :: MYRANKY = 0
   integer(kind=4) :: MYRANKZ = 0
   integer(kind=4) :: NRPROCX = 1
   integer(kind=4) :: NRPROCY = 1
   integer(kind=4) :: NRPROCZ = 1

contains

   subroutine FillDimVector(dims, shifts, nrcpp, stride, n, nrproc)
      integer(kind=4), allocatable, intent(out) :: dims(:), shifts(:)
      integer(kind=4), intent(in) :: nrcpp, stride, n, nrproc

      allocate(dims(nrproc), shifts(nrproc))
      dims = (n + 1)*stride
      shifts = 0
   end subroutine FillDimVector

end module parallelism
