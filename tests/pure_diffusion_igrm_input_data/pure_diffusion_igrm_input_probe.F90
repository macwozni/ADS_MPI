program pure_diffusion_igrm_input_probe
   use input_data, only: InitializeParameters, SIZE, ORDER, procx, procy, &
                         procz, steps, Dt, time_scheme
   implicit none

   call InitializeParameters
   write (*, '(6(I0,1X),I0,1X,A)') SIZE, ORDER, procx, procy, procz, &
      steps, nint(1000000.d0*Dt), trim(time_scheme)

end program pure_diffusion_igrm_input_probe
