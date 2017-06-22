module math

implicit none

! No predefined constant anywhere else
real (kind=8), parameter :: PI = 4.d0 * datan(1.d0) 


contains

! -------------------------------------------------------------------
! Linear interpolation between two values.
!
! t      - interpolation parameter
! x, y   - points to interpolate between (t=0 -> x, t=1 -> y)
! -------------------------------------------------------------------
function lerp(t, x, y) result (val)
real (kind=8) :: t, x, y, val

val = (1 - t) * x + t * y

end function


! -------------------------------------------------------------------
! Smoothly decaying function
! t in [0, r] -> 1
! t in [r, R] -> from 1 to 0
! t in [R, inf] -> 0
! -------------------------------------------------------------------
function falloff(r, Rr, t) result (fval)
real (kind=8) :: r, Rr, t
real (kind=8) :: h, fval

if (t < r) then
  fval = 1.d0
else if (t > Rr) then
  fval = 0.d0
else
  h = (t - r) / (Rr - r)
  fval = ((h - 1) * (h + 1)) ** 2
endif

end function


! -------------------------------------------------------------------
! C^1 bump function on [0, 1]. Function and its 1st derivative 
! vanish at the endpoints, maximum is attained for x = 0.5 and its
! value is 1.
!
! x  - argument
! -------------------------------------------------------------------
function bump(x) result (val)
real (kind=8), intent(in) :: x
real (kind=8) :: val

if (x > -1 .and. x < 1) then
  val = (x + 1)**2 * (x - 1)**2
else
  val = 0
endif

end function


! -------------------------------------------------------------------
! C^1 bump function on [0, 1]. Function and its 1st derivative 
! vanish at the endpoints, maximum is attained for x = 0.5 and its
! value is 1.
!
! x  - argument
! -------------------------------------------------------------------
function bump01(x) result (val)
real (kind=8), intent(in) :: x
real (kind=8) :: val

  val = bump(4 * (x - 0.5))

end function


function bump3d(r, Rr, x, y, z) result (val)
real (kind=8), intent(in) :: r, Rr, x, y, z
real (kind=8) :: val, t

t = norm2([x, y, z] - 0.5d0)
val = falloff(r/2, Rr/2, t)

end function

end module

