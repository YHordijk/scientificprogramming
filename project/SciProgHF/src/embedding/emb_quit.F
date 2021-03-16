      SUBROUTINE fde_quit(text)
         character(len=60) :: text
         write (*,'(a60)') text
         stop 
      END

