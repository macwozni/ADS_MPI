program igrm_stokes_input_probe
   use input_data, only: InitializeParameters, nelem, p_test, p_trial, &
                         procx, procy, procz
   implicit none

   call InitializeParameters
   write(*, '(I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)') &
      nelem, p_test, p_trial, procx, procy, procz

end program igrm_stokes_input_probe
