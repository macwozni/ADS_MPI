module parallelism
   implicit none

   integer(kind=4) :: MYRANK = 0
   integer(kind=4) :: MYRANKX = 0
   integer(kind=4) :: MYRANKY = 0
   integer(kind=4) :: MYRANKZ = 0
   integer(kind=4) :: NRPROC = 3
   integer(kind=4) :: NRPROCX = 3
   integer(kind=4) :: NRPROCY = 1
   integer(kind=4) :: NRPROCZ = 1

contains

   integer(kind=4) function LINEARINDEX(x, y, z) result(rank)
      integer(kind=4), intent(in) :: x, y, z

      rank = x + y*NRPROCX + z*NRPROCX*NRPROCY
   end function LINEARINDEX


   subroutine ComputeEndpoints(coordinate, process_count, n, p, nrcpp, &
                               ibeg, iend, mine, maxe)
      integer(kind=4), intent(in) :: coordinate, process_count, n, p
      integer(kind=4), intent(out) :: nrcpp, ibeg, iend, mine, maxe
      integer(kind=4) :: owned, remainder

      owned = (n + 1)/process_count
      remainder = mod(n + 1, process_count)
      nrcpp = owned
      if (remainder > 0) nrcpp = nrcpp + 1
      ibeg = coordinate*owned + min(coordinate, remainder) + 1
      iend = ibeg + owned - 1
      if (coordinate < remainder) iend = iend + 1
      mine = max(1, ibeg - p)
      maxe = min(n - p + 1, iend)
   end subroutine ComputeEndpoints

end module parallelism


module communicators
   implicit none

   integer(kind=4), dimension(3, 1, 1) :: processors = &
      reshape((/0, 1, 2/), (/3, 1, 1/))

end module communicators
