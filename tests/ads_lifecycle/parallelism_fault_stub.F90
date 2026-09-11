module parallelism
   implicit none

   integer(kind=4) :: MYRANK = 0
   integer(kind=4) :: MYRANKX = 0
   integer(kind=4) :: MYRANKY = 0
   integer(kind=4) :: MYRANKZ = 0
   integer(kind=4) :: NRPROC = 1
   integer(kind=4) :: NRPROCX = 1
   integer(kind=4) :: NRPROCY = 1
   integer(kind=4) :: NRPROCZ = 1

contains

   subroutine ComputeEndpoints(rank, nrproc, n, p, nrcpp, ibeg, iend, &
                               mine, maxe)
      integer(kind=4), intent(in) :: rank, nrproc, n, p
      integer(kind=4), intent(out) :: nrcpp, ibeg, iend, mine, maxe
      integer(kind=4) :: block_size, remainder

      block_size = (n + 1)/nrproc
      remainder = mod(n + 1, nrproc)
      nrcpp = (n + 1 + nrproc - 1)/nrproc
      ibeg = rank*block_size + min(rank, remainder) + 1
      iend = ibeg + block_size - 1
      if (rank < remainder) iend = iend + 1
      mine = max(ibeg - p - 1, 1)
      maxe = min(iend, n + 1 - p)
   end subroutine ComputeEndpoints


   subroutine FillDimVector(dims, shifts, nrcpp, stride, n, nrproc)
      integer(kind=4), allocatable, intent(out) :: dims(:), shifts(:)
      integer(kind=4), intent(in) :: nrcpp, stride, n, nrproc
      integer(kind=4) :: i, local_nrcpp, ibeg, iend, mine, maxe

      allocate (dims(nrproc), shifts(nrproc))
      do i = 1, nrproc
         call ComputeEndpoints(i - 1, nrproc, n, 0, local_nrcpp, ibeg, &
                               iend, mine, maxe)
         dims(i) = (iend - ibeg + 1)*stride
         shifts(i) = (ibeg - 1)*stride
      end do

      if (nrcpp < 0) continue
   end subroutine FillDimVector

end module parallelism
