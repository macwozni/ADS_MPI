! -------------------------------------------------------------------
! Converts integer to string
!
! Input:
! ------
! n          - integer to be converted
!
! Output:
! ------
! str        - output string
!
! -------------------------------------------------------------------
subroutine int2str(n, str)
USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT ! access computing environment
integer (kind=4), intent(in)  :: n     ! Integer to be converted
character(len=*), intent(out) :: str   ! String containing number
character(len=11) :: longstr           ! Longest is -2147483647

write(longstr,'(I11)') n
longstr = adjustl(longstr)

if (len_trim(longstr) > len(str)) then
  write(ERROR_UNIT,'(A,I3,A)') 'int2str: WARNING: can''t fit '//trim(longstr)// &
    ' into a ',len(str),'-character variable'
endif

str = longstr

end subroutine
