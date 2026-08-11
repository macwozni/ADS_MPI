program test_setup
   use Setup, only: ADS_Setup, ADS_compute_data
   implicit none

   type(ADS_Setup) :: original, copied
   type(ADS_compute_data) :: work
   integer :: checks, failures

   checks = 0
   failures = 0

   call check(.not. allocated(original%Ux), 'new setup has no knot storage')
   call check(.not. allocated(work%F), 'new compute data has no directional buffers')
   call check(.not. allocated(work%R), 'new compute data has no halo storage')

   original%n = (/3, 4, 5/)
   original%p = (/1, 2, 3/)
   original%tau = 0.125d0
   allocate(original%Ux(0:3), original%dimensionsX(0:1))
   original%Ux = (/0.d0, 0.d0, 1.d0, 1.d0/)
   original%dimensionsX = (/2, 3/)

   copied = original
   original%Ux(1) = 9.d0
   original%dimensionsX(0) = 99

   call check(all(copied%n == (/3, 4, 5/)), 'setup scalar metadata survives assignment')
   call check(lbound(copied%Ux, 1) == 0 .and. ubound(copied%Ux, 1) == 3, &
              'setup assignment preserves knot bounds')
   call check(copied%Ux(1) == 0.d0, 'setup assignment deep-copies allocatables')
   call check(copied%dimensionsX(0) == 2, 'integer allocatables are independent after assignment')

   allocate(work%F(0:1, 0:2))
   allocate(work%R(-1:1, 2:3, 5:5, 1:1))
   work%F = 4.d0
   work%R = 7.d0
   call check(all(shape(work%F) == (/2, 3/)), 'directional buffer rank and shape match the API')
   call check(all(lbound(work%R) == (/-1, 2, 5, 1/)), 'halo storage preserves global lower bounds')
   call check(all(ubound(work%R) == (/1, 3, 5, 1/)), 'halo storage preserves global upper bounds')

   if (failures /= 0) then
      write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, ' Setup checks)'
      stop 1
   end if
   write (*, '(A,I0,A)') 'OK (', checks, ' Setup checks)'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      checks = checks + 1
      if (.not. condition) then
         failures = failures + 1
         write (*, '(A,A)') 'FAIL: ', trim(label)
      end if
   end subroutine check

end program test_setup
