program igrm_pollution_input_probe
   use input_data, only: InitializeParameters, nelem, adapt_mesh, p_trial, &
                         c_trial, p_test, c_test, steps, procx, procy, &
                         procz, output_resolution
   implicit none

   integer(kind=4) :: adapt_value

   call InitializeParameters
   if (adapt_mesh) then
      adapt_value = 1
   else
      adapt_value = 0
   end if
   write (*, '(11(I0,1X))') nelem, adapt_value, p_trial, c_trial, p_test, &
      c_test, steps, procx, procy, procz, output_resolution

end program igrm_pollution_input_probe
