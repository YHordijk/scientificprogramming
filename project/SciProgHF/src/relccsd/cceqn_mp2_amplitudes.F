      SUBROUTINE amplitude_equation_mp2 (EPS,FVO,T1,S1,T2,S2)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Solves MP2 amplitude or lambda equations.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C     Intermediaries etc.
C
      REAL*8 EPS(*),FVO(*),T1(*),S1(*),T2(*),S2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "results.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 DIAG
      LOGICAL TEQ
C
C---------------Executable code--------------------------------------
C
C
C  INITIALIZE T2 AND T1 AMPLITUDES
C
      CALL XCOPY (NDIMT1,FVO,1,T1,1)
      CALL GETVVOO (T2)
      CALL DENOM (EPS,T1,T2,T1,T2)
      CALL ZCORE (T1,T2)
      !mlv
      CALL FCORE(T1,T2,EPS)
C
C  EVALUATE CORRELATION ENERGY
C
         CALL MP2ENG (T1,T2,FVO,S2,EMP2)
         CALL T1DIAG (DIAG,T1,S1)
C
         WRITE(IW,6071)
         WRITE(IW,6072) ESCF
         WRITE(IW,6073) EMP2
         ETOT = ESCF + EMP2
         WRITE(IW,6075) ESCF + EMP2
         WRITE(IW,6077) DIAG
C
 6071 FORMAT(//,'  MP2 results',/)
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 7000 FORMAT(//' Timing of routine MP2EQN :'
     &/' Before MP2EQN :',T30,F12.3,' seconds'
     &/' Energy, T1-Diag :',T30,F12.3,' seconds'
     &/' Untimed parts :',T30,F12.3,' seconds'
     &/' Total time in MP2EQN :',T30,F12.3,' seconds')
      END
