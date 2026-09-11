module Setup
   implicit none

   type ADS_Setup
      integer(kind=4) :: n(3) = 0
      integer(kind=4) :: p(3) = 0
      integer(kind=4) :: s(3) = 0
      integer(kind=4) :: nelem(3) = 0
      integer(kind=4) :: nrcpp(3) = 0
      integer(kind=4) :: ibeg(3) = 0
      integer(kind=4) :: iend(3) = 0
      real(kind=8), allocatable :: Ux(:), Uy(:), Uz(:)
      integer(kind=4), allocatable :: dimensionsX(:), dimensionsY(:), &
                                      dimensionsZ(:)
      integer(kind=4), allocatable :: shiftsX(:), shiftsY(:), shiftsZ(:)
   end type ADS_Setup

end module Setup
