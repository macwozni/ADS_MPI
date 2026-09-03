program igrm_heat_input_probe
   use input_data, only: InitializeParameters, nelem, p_test, p_trial, &
                         procx, procy, procz, steps, Dt, &
                         selected_time_scheme, time_scheme
   implicit none

   call InitializeParameters
   write(*, '(15(I0,1X),A)') nelem, p_test, p_trial, procx, procy, procz, &
      steps, nint(Dt*1.d6), selected_time_scheme, trim(time_scheme)

end program igrm_heat_input_probe
