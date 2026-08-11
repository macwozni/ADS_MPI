module Setup
   implicit none

   type ADS_Setup
      integer(kind=4), dimension(3) :: n, p, s, nelem
      real(kind=8), allocatable :: Ux(:), Uy(:), Uz(:)
   end type ADS_Setup

end module Setup
