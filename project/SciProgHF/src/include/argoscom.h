!    -- FILE: argoscom.h --
!    This is for the ECP routine
      INTEGER IGRPAG, ECPNUC, ECPNONT(1000),                            &
     &       ndptag, MAXANG(200), irowmax(200,200), ireparg(200,200,100) 
      COMMON /argoscom/ IGRPAG, ECPNUC, ECPNONT,                        & 
     &                  ndptag, MAXANG, irowmax, ireparg           
!
      INTEGER         ECPINT_PRINT, ECPINT_DEBUG
      COMMON /ECPINT/ ECPINT_PRINT, ECPINT_DEBUG
!
      INTEGER ndptags(0:7)
      data ndptags/ 0,0,0,0,1,1,1,7/
!    -- end of argoscom.h --
