program l2_input_probe
   use input_data, only: InitializeParameters, order, isizex, isizey, isizez, &
                         procx, procy, procz
   implicit none

   call InitializeParameters
   write(*, '(I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)') &
      isizex, isizey, isizez, order, procx, procy, procz

end program l2_input_probe
