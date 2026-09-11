program solution_output_error_contract_probe
   use Setup, only: ADS_Setup
   use parallelism, only: MYRANK
   use mpi, only: real_bcast_calls, character_bcast_calls, &
                  real_bcast_error, character_bcast_error, reset_mpi_stub
   use my_mpi, only: reset_my_mpi_stub
   use plot, only: spline_plot_calls, reset_plot_stub
   use vtk, only: vtk_calls, reset_vtk_stub
   use solution_output, only: PrintSolution
   implicit none

   integer(kind=4), parameter :: COEFFICIENT_BCAST_ERROR = 7101
   integer(kind=4), parameter :: FILENAME_BCAST_ERROR = 7102
   type(ADS_Setup) :: ads
   real(kind=8) :: part(2, 4)
   integer(kind=4) :: checks, failures, status

   checks = 0
   failures = 0
   MYRANK = 0
   call prepare_inputs(ads, part)

   call reset_stubs()
   real_bcast_error = COEFFICIENT_BCAST_ERROR
   status = -1
   call PrintSolution(61, ads, part, status)
   call check('coefficient Bcast status reaches the caller exactly', &
              status == COEFFICIENT_BCAST_ERROR)
   call check('coefficient Bcast status stops all later output work', &
              real_bcast_calls == 1 .and. character_bcast_calls == 0 .and. &
              spline_plot_calls == 0 .and. vtk_calls == 0)

   call reset_stubs()
   character_bcast_error = FILENAME_BCAST_ERROR
   status = -1
   call PrintSolution(62, ads, part, status)
   call check('filename Bcast status reaches the caller exactly', &
              status == FILENAME_BCAST_ERROR)
   call check('filename Bcast status prevents output with partial metadata', &
              real_bcast_calls == 1 .and. character_bcast_calls == 1 .and. &
              spline_plot_calls == 0 .and. vtk_calls == 0)

   if (failures == 0) then
      write (*, '(A,I0,A)') 'OK (', checks, ' solution-output error-contract checks)'
   else
      write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                 ' solution-output error-contract checks)'
      stop 1
   end if

contains

   subroutine prepare_inputs(space, local_part)
      type(ADS_Setup), intent(out) :: space
      real(kind=8), intent(out) :: local_part(2, 4)
      integer :: i

      space%n = (/1, 1, 1/)
      space%p = (/1, 1, 1/)
      space%s = (/2, 2, 2/)
      space%nelem = (/1, 1, 1/)
      allocate(space%Ux(0:3), space%Uy(0:3), space%Uz(0:3))
      do i = 0, 3
         space%Ux(i) = real(100 + i, kind=8)
         space%Uy(i) = real(200 + i, kind=8)
         space%Uz(i) = real(300 + i, kind=8)
      end do
      local_part = reshape((/1.d0, 2.d0, 3.d0, 4.d0, &
                             5.d0, 6.d0, 7.d0, 8.d0/), shape(local_part))
   end subroutine prepare_inputs


   subroutine reset_stubs()
      call reset_mpi_stub()
      call reset_my_mpi_stub()
      call reset_plot_stub()
      call reset_vtk_stub()
   end subroutine reset_stubs


   subroutine check(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      checks = checks + 1
      if (condition) then
         write (*, '(A)') 'PASS '//trim(label)
      else
         failures = failures + 1
         write (*, '(A)') 'FAIL '//trim(label)
      end if
   end subroutine check

end program solution_output_error_contract_probe
