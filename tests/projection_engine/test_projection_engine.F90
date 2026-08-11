program test_projection_engine
   use projection_engine, only: ValidateIGRMMesh, MKBBT_large, MKBBT_small, &
                                ComputeMatrix, Form3DRHS, create_mixed_space, &
                                FormUn, global2local
   implicit none

   integer :: marker

   marker = 0
   call ValidateIGRMMesh(marker)
   call MKBBT_large(marker)
   call MKBBT_small(marker)
   call ComputeMatrix(marker)
   call Form3DRHS(marker)
   call create_mixed_space(marker)
   call FormUn(marker)
   call global2local(marker)

   if (marker /= 255) then
      write (*, '(A,I0)') 'FAIL: projection_engine did not re-export every API entry, marker=', marker
      stop 1
   end if

   write (*, '(A)') 'OK (projection_engine compatibility exports)'
end program test_projection_engine
