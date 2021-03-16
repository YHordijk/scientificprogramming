      SUBROUTINE ORDER2(EVEC,EVAL,N,NEVEC)
C
C Purpose: order the N values in EVAL and their associated vectors
C          in EVEC so EVAL(i+1) .le. EVAL(i)
C          (this is opposite order of "ORDER")
C
C Revisions:
C   29-Jul-1992 hjaaj (only dswap if nevec .gt. 0)
C    5-Aug-1985 hjaaj (first version, based on ORDER)
C
#include "implicit.h"
      DIMENSION EVEC(*),EVAL(*)
      IF (N.LE.1) RETURN
      IN = 1
      DO 10 I=1,N-1
         EMAX = EVAL(I)
         IMAX = I
         DO 20 J=I+1,N
            IF (EVAL(J) .GT. EMAX) THEN
               EMAX = EVAL(J)
!radovan: this file is under safeflags because imax gets strange values
!         here with optimized ifort 12.0.0 32bit
               IMAX = J
            ENDIF
   20    CONTINUE
         IF (IMAX.NE.I) THEN
            EVAL(IMAX)=EVAL(I)
            EVAL(I)=EMAX
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMAX-1)*NEVEC+1),1)
            ENDIF
         ENDIF
         IN = IN + NEVEC
   10 CONTINUE
      END
