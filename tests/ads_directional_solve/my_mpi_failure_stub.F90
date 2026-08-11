module my_mpi
   implicit none

   integer(kind=4) :: transport_operation_count = 0
   integer(kind=4) :: gather_call_count = 0
   integer(kind=4) :: scatter_call_count = 0
   integer(kind=4) :: failing_operation = 0
   integer(kind=4) :: injected_transport_status = 0

contains

   subroutine configure_transport_failure(operation, status)
      integer(kind=4), intent(in) :: operation, status

      transport_operation_count = 0
      gather_call_count = 0
      scatter_call_count = 0
      failing_operation = operation
      injected_transport_status = status
   end subroutine configure_transport_failure


   subroutine Gather(F, F_out, n, elems, stride, dims, shifts, comm, ierr)
      real(kind=8), intent(in) :: F(:, :)
      real(kind=8), intent(out) :: F_out(:, :)
      integer(kind=4), intent(in) :: n, elems, stride, comm
      integer(kind=4), intent(in) :: dims(:), shifts(:)
      integer(kind=4), intent(out) :: ierr

      transport_operation_count = transport_operation_count + 1
      gather_call_count = gather_call_count + 1
      ierr = 0
      if (transport_operation_count == failing_operation) then
         ierr = injected_transport_status
         return
      end if

      if (any(shape(F_out) /= shape(F))) then
         stop 91
      end if
      F_out = F
   end subroutine Gather


   subroutine Scatter(F, F_out, n, elems, stride, dims, shifts, comm, ierr)
      real(kind=8), intent(in) :: F(:, :)
      real(kind=8), intent(out) :: F_out(:, :)
      integer(kind=4), intent(in) :: n, elems, stride, comm
      integer(kind=4), intent(in) :: dims(:), shifts(:)
      integer(kind=4), intent(out) :: ierr

      transport_operation_count = transport_operation_count + 1
      scatter_call_count = scatter_call_count + 1
      ierr = 0
      if (transport_operation_count == failing_operation) then
         ierr = injected_transport_status
         return
      end if

      if (any(shape(F_out) /= shape(F))) then
         stop 92
      end if
      F_out = F
   end subroutine Scatter

end module my_mpi
