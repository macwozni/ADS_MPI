subroutine dmumps(mumps_par)
   use mumps_test_support, only: failure_job, failure_occurrence, failure_code, &
      finalization_code, failure_uses_infog, record_mumps_job
   use mumps_test_support, only: record_mumps_contract
   implicit none
   include 'dmumps_struc.h'
   type(dmumps_struc), intent(inout) :: mumps_par
   integer(kind=4) :: occurrence
   logical :: fail_this_call

   call record_mumps_job(mumps_par%job, occurrence)
   if (mumps_par%job == 1) then
      call record_mumps_contract(mumps_par%comm, mumps_par%sym, mumps_par%par, &
         mumps_par%n, mumps_par%nz, mumps_par%irn, mumps_par%jcn, &
         mumps_par%a, mumps_par%icntl)
   end if
   mumps_par%info = 0
   mumps_par%infog = 0

   if (mumps_par%job == -1) then
      allocate(mumps_par%intr_encoding(1))
   else if (mumps_par%job == -2) then
      if (associated(mumps_par%intr_encoding)) then
         deallocate(mumps_par%intr_encoding)
      end if
   end if

   fail_this_call = mumps_par%job == failure_job .and. &
                    occurrence == failure_occurrence
   if (fail_this_call) then
      if (failure_uses_infog) then
         mumps_par%infog(1) = failure_code
      else
         mumps_par%info(1) = failure_code
      end if
   else if (mumps_par%job == -2 .and. finalization_code /= 0) then
      mumps_par%info(1) = finalization_code
   else if (mumps_par%job == 3) then
      mumps_par%rhs = 2.d0*mumps_par%rhs
   end if
end subroutine dmumps
