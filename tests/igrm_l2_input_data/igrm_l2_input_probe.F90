program igrm_l2_input_probe
   use input_data, only: InitializeParameters, nelem, p_test, p_trial, &
                         procx, procy, procz, scheme_tau, time_scheme, &
                         selected_time_scheme
   implicit none

   call InitializeParameters
   write(*, '(12(I0,1X),I0,1X,A,1X,F0.6)') nelem, p_test, p_trial, &
      procx, procy, procz, selected_time_scheme, trim(time_scheme), scheme_tau

end program igrm_l2_input_probe
