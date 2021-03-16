      SUBROUTINE VIBNOR(HESS,AMASS,DKIN,HESSMW,
     &                  EVAL,EVEC,EVECS,WRK,LWRK,
     &                  NCORD,NINTRN,NNEG,NOTALL,IPRINT)
#include "implicit.h"
#include "codata.h"
C
C     Find normal coordinates in DNORCR and normal frequences in FREQ
C     for projected Hessian eigenvalue problem
C     (HESS - EVAL * DKIN) * EVEC = 0
C
C     Modification KR 29.10.96 - Sort only according to frequency, not symmetry
C     Modificaton TUH 26.12.86 - Sort according to symmetry
C     Modified TUH 21.09.89 - Explicit Abelian symmetry
C
C     Input:   HESS,AMASS
C     Output:  EVAL,EVEC,EVECS,NNEG
C
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 ,
     &            THRSH1 = 1.0D-14, THRSH2 = 1.0D-3, THRSH3 = 1.0D-7)
      LOGICAL NOTALL
C
      DIMENSION HESS(NCORD,NCORD),
     &          AMASS(*), DKIN(NCORD), HESSMW(NCORD*(NCORD+1)/2),
     &          EVAL(NCORD), EVEC(NCORD,NCORD), EVECS(NCORD,NCORD),
     &          WRK(LWRK)
C
#include "abainf.h"
#include "symmet.h"
C
C     *******************************************
C     ***** Construct mass-weighted Hessian *****
C     *******************************************
C
      DO 100 I = 1, NCORD
         DKIN(I) = D1/SQRT(XFAMU*AMASS((I+2)/3))
 100  CONTINUE
      IJ = 0
      DO 200 I = 1, NCORD
         DO 210 J = 1, I
            IJ = IJ + 1
            HESSMW(IJ) = DKIN(J)*HESS(J,I)*DKIN(I)
 210     CONTINUE
 200  CONTINUE
C
C     *********************************************
C     ***** Diagonalize mass-weighted Hessian *****
C     *********************************************
C
      KWRK  = 1
      KIWRK = KWRK + NCORD
      CALL DUNIT(EVEC,NCORD)
      CALL JACO(HESSMW,EVEC,NCORD,NCORD,NCORD,WRK(KWRK),WRK(KIWRK))
C
C     **********************************************
C     ***** Eigenvalues and normal coordinates *****
C     **********************************************
C
      II = 0
      DO 300 I = 1,NCORD
         II = II + I
         EVAL(I) = HESSMW(II)
         DO 310 J = 1,NCORD
            EVEC(J,I) = DKIN(J)*EVEC(J,I)
 310     CONTINUE
 300  CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Eigenvalues of mass-weighted Hessian',-1)
         CALL OUTPUT(EVAL,1,1,1,NCORD,1,NCORD,1,LUPRI)
         CALL HEADER('Normal coordinates in Cartesian basis',-1)
         CALL OUTPUT(EVEC,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
      END IF
C
C     *******************************************
C     ***** Eliminate non-vibrational modes *****
C     *******************************************
C
      NZER = 0
      NNEG = 0
C/djw/ VALTHR set to 1.0D-10
C     VALTHR = THRSH1
C     IF (NUMHES) VALTHR = 1.0D-10
      VALTHR = 1.0D-10
      DO 400 I = 1, NCORD
         IF (ABS(EVAL(I)) .LE.  VALTHR) NZER = NZER + 1
         IF (    EVAL(I)  .LT. -VALTHR) NNEG = NNEG + 1
 400  CONTINUE
      CALL ORDER2(EVEC,EVAL,NCORD,NCORD)
      IZERO  = NCORD - NNEG - NZER + 1
      DO 500 INEG = NCORD - NNEG + 1, NCORD
         CALL DSWAP(NCORD,EVEC(1,INEG),1,EVEC(1,IZERO),1)
         CALL DSWAP(1,EVAL(INEG),1,EVAL(IZERO),1)
         IZERO = IZERO + 1
 500  CONTINUE
      NINTRN = NCORD - NZER
      IF (.NOT.NOTALL .AND. (NZER .NE. 5 .AND. NZER .NE. 6)) THEN
         WRITE (LUPRI,'(/A,I2,/,A)') ' The number of translational'//
     &      ' and rotational coordinates in VIBNOR',NZER,' is'//
     &      ' not equal to 5 or 6.'
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER
     &      ('Non-zero eigenvalues of mass-weighted Hessian',-1)
         CALL OUTPUT(EVAL,1,1,1,NINTRN,1,NCORD,1,LUPRI)
         CALL HEADER
     &      ('Non-zero normal coordinates in Cartesian basis',-1)
         CALL OUTPUT(EVEC,1,NCORD,1,NINTRN,NCORD,NCORD,1,LUPRI)
      END IF
C
C     **************************************
C     ***** Transform to symmetry basis ****
C     **************************************
C
      KCSTRA = 1
      KSCTRA = KCSTRA + NCORD*NCORD
      KLAST  = KSCTRA + NCORD*NCORD
      IF (KLAST .GT. LWRK) CALL STOPIT('VIBNOR','TRACOR',KLAST,LWRK)
      CALL TRACOR(WRK(KCSTRA),WRK(KSCTRA),1,NCORD,0)
C
      CALL DGEMM('T','N',NCORD,NCORD,NCORD,1.D0,
     &           WRK(KSCTRA),NCORD,
     &           EVEC,NCORD,0.D0,
     &           EVECS,NCORD)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Normal coordinates in symmetry basis',-1)
         CALL OUTPUT(EVECS,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
      END IF
C
C     **********************************************************
C     ***** Sort normal modes according to their frequency *****
C     **********************************************************
C
      DO 600 ICOOR = 1, NINTRN
         DO 610 JCOOR = 1, NINTRN
            IF (EVAL(ICOOR) .GT. EVAL(JCOOR)) THEN
               CALL DSWAP(1,EVAL(ICOOR),1,EVAL(JCOOR),1)
               CALL DSWAP(NCORD,EVEC(1,ICOOR),1,EVEC(1,JCOOR),1)
               IF (MAXREP .GT. 1) THEN
                  CALL DSWAP(NCORD,EVECS(1,ICOOR),1,EVECS(1,JCOOR),1)
               END IF
            END IF
 610     CONTINUE 
 600  CONTINUE 
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Sorted eigenvalues',-1)
         CALL OUTPUT(EVAL,1,1,1,NINTRN,1,NCORD,1,LUPRI)
         CALL HEADER('Sorted normal coordinates in Cartesian basis',-1)
         CALL OUTPUT(EVEC,1,NCORD,1,NINTRN,NCORD,NCORD,1,LUPRI)
         CALL HEADER('Sorted normal coordinates in symmetry basis',-1)
         CALL OUTPUT(EVECS,1,NCORD,1,NINTRN,NCORD,NCORD,1,LUPRI)
      END IF
      END
