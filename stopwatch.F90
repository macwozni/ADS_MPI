
module stopwatch


contains

! -------------------------------------------------------------------
! Writes current system clock to specified location.
!
! Icount    - value of the clock
! -------------------------------------------------------------------
subroutine start_clock(Icount)
#include "syscom.blk"

  call system_clock(Icount,ir,im)

end subroutine


! -------------------------------------------------------------------
! Computes time difference between specified point in time and
! current time.
!
! Dtime     - difference between current and previous clock values
! Icount    - input:  previous clock value
!             output: current clock value
! -------------------------------------------------------------------
subroutine stop_clock(Dtime,Icount)
#include "syscom.blk"

  call system_clock(ic,ir,im)

  if (ic > Icount) then
    Dtime = dble(ic-Icount)/dble(ir)
  elseif (ic < Icount) then
    Dtime = dble(im+ic-Icount)/dble(ir)
  else
    Dtime = 0.d0
  endif

end subroutine


! -------------------------------------------------------------------
! Stops the clock, and immediately starts it again.
!
! Dtime     - difference between current and previous clock values
! Icount    - input:  previous clock value
!             output: current clock value
! -------------------------------------------------------------------
subroutine reset_clock(Dtime,Icount)
#include "syscom.blk"

  call stop_clock(Dtime,Icount)
  call start_clock(Icount)

end subroutine

end module

