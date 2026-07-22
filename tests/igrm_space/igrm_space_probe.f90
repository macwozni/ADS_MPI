program igrm_space_probe
   use ISO_FORTRAN_ENV, only: ERROR_UNIT
   use igrm_space, only: ValidateIGRMMesh
   implicit none

   character(len=32) :: mode
   integer(kind=4) :: nelem_test, nelem_trial, n_test, n_trial, p_test, p_trial
   real(kind=8), allocatable :: U_test(:), U_trial(:)

   if (command_argument_count() /= 1) stop 90
   call get_command_argument(1, mode)

   select case (trim(mode))
   case ("equal-degree")
      p_test = 1
      n_test = 3
      nelem_test = 2
      p_trial = 1
      n_trial = 2
      nelem_trial = 2
      allocate(U_test(0:5), U_trial(0:4))
      U_test = (/ 0.d0, 0.d0, 0.5d0, 0.5d0, 1.d0, 1.d0 /)
      U_trial = (/ 0.d0, 0.d0, 0.5d0, 1.d0, 1.d0 /)

   case ("lower-degree")
      p_test = 1
      n_test = 2
      nelem_test = 2
      p_trial = 2
      n_trial = 3
      nelem_trial = 2
      allocate(U_test(0:4), U_trial(0:6))
      U_test = (/ 0.d0, 0.d0, 0.5d0, 1.d0, 1.d0 /)
      U_trial = (/ 0.d0, 0.d0, 0.d0, 0.5d0, 1.d0, 1.d0, 1.d0 /)

   case ("shifted-location")
      p_test = 2
      n_test = 3
      nelem_test = 2
      p_trial = 1
      n_trial = 2
      nelem_trial = 2
      allocate(U_test(0:6), U_trial(0:4))
      U_test = (/ 0.d0, 0.d0, 0.d0, 0.5d0 + 2.d-12, 1.d0, 1.d0, 1.d0 /)
      U_trial = (/ 0.d0, 0.d0, 0.5d0, 1.d0, 1.d0 /)

   case ("extra-test-location")
      p_test = 2
      n_test = 3
      nelem_test = 2
      p_trial = 1
      n_trial = 1
      nelem_trial = 2
      allocate(U_test(0:6), U_trial(0:3))
      U_test = (/ 0.d0, 0.d0, 0.d0, 0.5d0, 1.d0, 1.d0, 1.d0 /)
      U_trial = (/ 0.d0, 0.d0, 1.d0, 1.d0 /)

   case ("extra-trial-location")
      p_test = 2
      n_test = 2
      nelem_test = 2
      p_trial = 1
      n_trial = 2
      nelem_trial = 2
      allocate(U_test(0:5), U_trial(0:4))
      U_test = (/ 0.d0, 0.d0, 0.d0, 1.d0, 1.d0, 1.d0 /)
      U_trial = (/ 0.d0, 0.d0, 0.5d0, 1.d0, 1.d0 /)

   case ("metadata-mismatch")
      p_test = 2
      n_test = 3
      nelem_test = 2
      p_trial = 1
      n_trial = 2
      nelem_trial = 3
      allocate(U_test(0:6), U_trial(0:4))
      U_test = (/ 0.d0, 0.d0, 0.d0, 0.5d0, 1.d0, 1.d0, 1.d0 /)
      U_trial = (/ 0.d0, 0.d0, 0.5d0, 1.d0, 1.d0 /)

   case default
      write(ERROR_UNIT, '(A)') "unknown probe mode"
      stop 90
   end select

   call ValidateIGRMMesh(U_test, p_test, n_test, nelem_test, &
                         U_trial, p_trial, n_trial, nelem_trial)

   write(ERROR_UNIT, '(A)') "validation unexpectedly succeeded"
   stop 98

end program igrm_space_probe
