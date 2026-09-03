program eriksson_input_probe
   use input_data, only: InitializeParameters, SIZE, ORDER, steps, Dt, &
                         procx, procy, procz
   implicit none

   call InitializeParameters
   write (*, '(7(I0,1X))') SIZE, ORDER, steps, nint(1000000.d0*Dt), &
      procx, procy, procz

end program eriksson_input_probe
