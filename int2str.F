c
      subroutine int2str(int, str)
      implicit none
c     In
      integer, intent(IN) :: int           ! Integer to be converted
c     Out
      character(len=*), intent(OUT) :: str ! String containing number
c     Local
      character(len=11) :: longstr         ! Longest is -2147483647
c
      write(longstr,'(I11)') int
      longstr = adjustl(longstr)
      if (len_trim(longstr) > len(str)) then
        write(*,7001) 'int2str: WARNING: can''t fit '//trim(longstr)//
     .                ' into a ',len(str),'-character variable'
 7001   format(A,I3,A)
      endif
      str = longstr
      end subroutine int2str
