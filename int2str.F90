
subroutine int2str(n, str)
integer, intent(in) :: n             ! Integer to be converted
character(len=*), intent(out) :: str ! String containing number
character(len=11) :: longstr         ! Longest is -2147483647

  write(longstr,'(I11)') n
  longstr = adjustl(longstr)

  if (len_trim(longstr) > len(str)) then
    write(*,'(A,I3,A)') 'int2str: WARNING: can''t fit '//trim(longstr)// &
      ' into a ',len(str),'-character variable'
  endif

  str = longstr

end subroutine
