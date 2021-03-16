C
C     Fine structure constant; taken from DIRAC
C
       REAL*8 CVEL,ALPHAC,ALPHA2
CMI parameter (speed2=CVEL*CVEL), parameter (speed4=speed2*speed2) 
       REAL*8 speed2,speed4
       COMMON /AMFICB/ CVEL,ALPHAC,ALPHA2,speed2,speed4
