C
C
C  /* Deck dspsol */
      SUBROUTINE DSPSOL (N,NSIM,AP,B,KPVT,INFO)
C
C Written 8-Feb-1985 Hans Joergen Aa. Jensen
C No revisions.
C
C Purpose:
C  Solve the NSIM simultaneous eqautions:
C
C     B(n,nsim) := A(n,n) inverse * B(n,nsim)
C
C  AP is A in packed format.
C  KPVT is a scratch array of length N.
C
#include "implicit.h"
      DIMENSION AP(*),B(N,*),KPVT(*)
C
#if !defined (VAR_SBLAS)
      CALL DSPFA (AP,N,KPVT,INFO)
#else
      CALL SSPFA (AP,N,KPVT,INFO)
#endif
      IF (INFO.NE.0) RETURN
C
      DO 100 J = 1,NSIM
#if !defined (VAR_SBLAS)
        CALL DSPSL (AP,N,KPVT,B(1,J))
#else
        CALL SSPSL (AP,N,KPVT,B(1,J))
#endif
  100 CONTINUE
C
      RETURN
      END
C  /* Deck dspsli */
      SUBROUTINE DSPSLI (N,NSIM,AP,B,KPVT,INFO,DET,INERT)
C
C Written 24-Feb-1989 Hans Joergen Aa. Jensen, based on DSPSOL.
C No revisions.
C
C Purpose:
C  Solve the NSIM simultaneous eqautions:
C
C     B(n,nsim) := A(n,n) inverse * B(n,nsim)
C
C  AP is A in packed format.
C  KPVT is a scratch array of length N.
C
#include "implicit.h"
      DIMENSION AP(*),B(N,*),KPVT(*),DET(2),INERT(3)
C
#if !defined (VAR_SBLAS)
      CALL DSPFA (AP,N,KPVT,INFO)
#else
      CALL SSPFA (AP,N,KPVT,INFO)
#endif
      IF (INFO.NE.0) RETURN
C
#if !defined (VAR_SBLAS)
      CALL DSPDI(AP,N,KPVT,DET,INERT,DUMMY,110)
C     CALL DSPDI(AP,N,KPVT,DET,INERT,WORK,JOB)
#else
      CALL SSPDI(AP,N,KPVT,DET,INERT,DUMMY,110)
C     CALL SSPDI(AP,N,KPVT,DET,INERT,WORK,JOB)
#endif
C
      DO 100 J = 1,NSIM
#if !defined (VAR_SBLAS)
        CALL DSPSL (AP,N,KPVT,B(1,J))
#else
        CALL SSPSL (AP,N,KPVT,B(1,J))
#endif
  100 CONTINUE
C
      RETURN
      END
C  /* Deck dsplinpck */
#if !defined (VAR_NOLINPACK)
      SUBROUTINE DSPCO(AP,N,KPVT,RCOND,Z)
      INTEGER N,KPVT(*)
      DOUBLE PRECISION AP(*),Z(*)
      DOUBLE PRECISION RCOND
C
C     DSPCO FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX STORED IN
C     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING AND ESTIMATES
C     THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DSPFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DSPCO BY DSPSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPCO BY DSPSL.
C     TO COMPUTE  INVERSE(A) , FOLLOW DSPCO BY DSPDI.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPCO BY DSPDI.
C     TO COMPUTE  INERTIA(A), FOLLOW DSPCO BY DSPDI.
C
C     ON ENTRY
C
C        AP      DOUBLE PRECISION (N*(N+1)/2)
C                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
C                SEE COMMENTS BELOW FOR DETAILS.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     OUTPUT
C
C        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     PACKED STORAGE
C
C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
C          TRIANGLE OF A SYMMETRIC MATRIX.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K) = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DSPFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN DABS,DMAX1,IABS,DSIGN
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,EK,T
      DOUBLE PRECISION ANORM,S,DASUM,YNORM
      INTEGER I,IJ,IK,IKM1,IKP1,INFO,J,JM1,J1
      INTEGER K,KK,KM1K,KM1KM1,KP,KPS,KS
C
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
      J1 = 1
      DO 30 J = 1, N
         Z(J) = DASUM(J,AP(J1),1)
         IJ = J1
         J1 = J1 + J
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + DABS(AP(IJ))
            IJ = IJ + 1
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
C
C     FACTOR
C
      CALL DSPFA(AP,N,KPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE U*D*W = E
C
      EK = 1.0D0
      DO 50 J = 1, N
         Z(J) = 0.0D0
   50 CONTINUE
      K = N
      IK = (N*(N - 1))/2
   60 IF (K .EQ. 0) GO TO 120
         KK = IK + K
         IKM1 = IK - (K - 1)
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = IABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .EQ. KPS) GO TO 70
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
   70    CONTINUE
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL DAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
         IF (KS .EQ. 1) GO TO 80
            IF (Z(K-1) .NE. 0.0D0) EK = DSIGN(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL DAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
   80    CONTINUE
         IF (KS .EQ. 2) GO TO 100
            IF (DABS(Z(K)) .LE. DABS(AP(KK))) GO TO 90
               S = DABS(AP(KK))/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   90       CONTINUE
            IF (AP(KK) .NE. 0.0D0) Z(K) = Z(K)/AP(KK)
            IF (AP(KK) .EQ. 0.0D0) Z(K) = 1.0D0
         GO TO 110
  100    CONTINUE
            KM1K = IK + K - 1
            KM1KM1 = IKM1 + K - 1
            AK = AP(KK)/AP(KM1K)
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = Z(K)/AP(KM1K)
            BKM1 = Z(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
         K = K - KS
         IK = IK - K
         IF (KS .EQ. 2) IK = IK - (K + 1)
      GO TO 60
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(U)*Y = W
C
      K = 1
      IK = 0
  130 IF (K .GT. N) GO TO 160
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 150
            Z(K) = Z(K) + DDOT(K-1,AP(IK+1),1,Z(1),1)
            IKP1 = IK + K
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + DDOT(K-1,AP(IKP1+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 140
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  140       CONTINUE
  150    CONTINUE
         IK = IK + K
         IF (KS .EQ. 2) IK = IK + (K + 1)
         K = K + KS
      GO TO 130
  160 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE U*D*V = Y
C
      K = N
      IK = N*(N - 1)/2
  170 IF (K .EQ. 0) GO TO 230
         KK = IK + K
         IKM1 = IK - (K - 1)
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. KS) GO TO 190
            KP = IABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .EQ. KPS) GO TO 180
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
  180       CONTINUE
            CALL DAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
            IF (KS .EQ. 2) CALL DAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
  190    CONTINUE
         IF (KS .EQ. 2) GO TO 210
            IF (DABS(Z(K)) .LE. DABS(AP(KK))) GO TO 200
               S = DABS(AP(KK))/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  200       CONTINUE
            IF (AP(KK) .NE. 0.0D0) Z(K) = Z(K)/AP(KK)
            IF (AP(KK) .EQ. 0.0D0) Z(K) = 1.0D0
         GO TO 220
  210    CONTINUE
            KM1K = IK + K - 1
            KM1KM1 = IKM1 + K - 1
            AK = AP(KK)/AP(KM1K)
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = Z(K)/AP(KM1K)
            BKM1 = Z(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
         K = K - KS
         IK = IK - K
         IF (KS .EQ. 2) IK = IK - (K + 1)
      GO TO 170
  230 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE TRANS(U)*Z = V
C
      K = 1
      IK = 0
  240 IF (K .GT. N) GO TO 270
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 260
            Z(K) = Z(K) + DDOT(K-1,AP(IK+1),1,Z(1),1)
            IKP1 = IK + K
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + DDOT(K-1,AP(IKP1+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 250
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  250       CONTINUE
  260    CONTINUE
         IK = IK + K
         IF (KS .EQ. 2) IK = IK + (K + 1)
         K = K + KS
      GO TO 240
  270 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
      SUBROUTINE DSPFA(AP,N,KPVT,INFO)
      INTEGER N,KPVT(*),INFO
      DOUBLE PRECISION AP(*)
C
C     DSPFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX STORED IN
C     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.
C
C     TO SOLVE  A*X = B , FOLLOW DSPFA BY DSPSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPFA BY DSPSL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPFA BY DSPDI.
C     TO COMPUTE  INERTIA(A) , FOLLOW DSPFA BY DSPDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DSPFA BY DSPDI.
C
C     ON ENTRY
C
C        AP      DOUBLE PRECISION (N*(N+1)/2)
C                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
C                SEE COMMENTS BELOW FOR DETAILS.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     OUTPUT
C
C        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
C                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
C                     BUT IT DOES INDICATE THAT DSPSL OR DSPDI MAY
C                     DIVIDE BY ZERO IF CALLED.
C
C     PACKED STORAGE
C
C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
C          TRIANGLE OF A SYMMETRIC MATRIX.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K)  = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSWAP,IDAMAX
C     FORTRAN DABS,DMAX1,DSQRT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IDAMAX,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
      LOGICAL SWAP
C
C
C     INITIALIZE
C
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0
C
      INFO = 0
C
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
C
      K = N
      IK = (N*(N - 1))/2
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
C     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (AP(1) .EQ. 0.0D0) INFO = 1
C     ......EXIT
            GO TO 200
   20    CONTINUE
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
C        REQUIRED.
C
         KM1 = K - 1
         KK = IK + K
         ABSAKK = DABS(AP(KK))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = IDAMAX(K-1,AP(IK+1),1)
         IMK = IK + IMAX
         COLMAX = DABS(AP(IMK))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
C
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C           ROW IMAX.
C
            ROWMAX = 0.0D0
            IMAXP1 = IMAX + 1
            IM = IMAX*(IMAX - 1)/2
            IMJ = IM + 2*IMAX
            DO 40 J = IMAXP1, K
               ROWMAX = DMAX1(ROWMAX,DABS(AP(IMJ)))
               IMJ = IMJ + J
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = IDAMAX(IMAX-1,AP(IM+1),1)
               JMIM = JMAX + IM
               ROWMAX = DMAX1(ROWMAX,DABS(AP(JMIM)))
   50       CONTINUE
            IMIM = IMAX + IM
            IF (DABS(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100
C
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
C
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
C
C           1 X 1 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 120
C
C              PERFORM AN INTERCHANGE.
C
               CALL DSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)
               IMJ = IK + IMAX
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  JK = IK + J
                  T = AP(JK)
                  AP(JK) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  110          CONTINUE
  120       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            IJ = IK - (K - 1)
            DO 130 JJ = 1, KM1
               J = K - JJ
               JK = IK + J
               MULK = -AP(JK)/AP(KK)
               T = MULK
               CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
               IJJ = IJ + J
               AP(JK) = MULK
               IJ = IJ - (J - 1)
  130       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            KM1K = IK + K - 1
            IKM1 = IK - (K - 1)
            IF (.NOT.SWAP) GO TO 160
C
C              PERFORM AN INTERCHANGE.
C
               CALL DSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)
               IMJ = IKM1 + IMAX
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  JKM1 = IKM1 + J
                  T = AP(JKM1)
                  AP(JKM1) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  150          CONTINUE
               T = AP(KM1K)
               AP(KM1K) = AP(IMK)
               AP(IMK) = T
  160       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = AP(KK)/AP(KM1K)
               KM1KM1 = IKM1 + K - 1
               AKM1 = AP(KM1KM1)/AP(KM1K)
               DENOM = 1.0D0 - AK*AKM1
               IJ = IK - (K - 1) - (K - 2)
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  JK = IK + J
                  BK = AP(JK)/AP(KM1K)
                  JKM1 = IKM1 + J
                  BKM1 = AP(JKM1)/AP(KM1K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
                  T = MULKM1
                  CALL DAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)
                  AP(JK) = MULK
                  AP(JKM1) = MULKM1
                  IJJ = IJ + J
                  IJ = IJ - (J - 1)
  170          CONTINUE
  180       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         IK = IK - (K - 1)
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
      SUBROUTINE DSPSL(AP,N,KPVT,B)
      INTEGER N,KPVT(*)
      DOUBLE PRECISION AP(*),B(*)
C
C     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM
C     A * X = B
C     USING THE FACTORS COMPUTED BY DSPFA.
C
C     ON ENTRY
C
C        AP      DOUBLE PRECISION(N*(N+1)/2)
C                THE OUTPUT FROM DSPFA.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DSPFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO MAY OCCUR IF  DSPCO  HAS SET RCOND .EQ. 0.0
C        OR  DSPFA  HAS SET INFO .NE. 0  .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DSPFA(AP,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DSPSL(AP,N,KPVT,C(1,J))
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C     FORTRAN IABS
C
C     INTERNAL VARIABLES.
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP
      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
      K = N
      IK = (N*(N - 1))/2
   10 IF (K .EQ. 0) GO TO 80
         KK = IK + K
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL DAXPY(K-1,B(K),AP(IK+1),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/AP(KK)
            K = K - 1
            IK = IK - K
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IKM1 = IK - (K - 1)
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 INTERCHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL DAXPY(K-2,B(K),AP(IK+1),1,B(1),1)
               CALL DAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            KM1K = IK + K - 1
            KK = IK + K
            AK = AP(KK)/AP(KM1K)
            KM1KM1 = IKM1 + K - 1
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = B(K)/AP(KM1K)
            BKM1 = B(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0D0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
            IK = IK - (K + 1) - K
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
      IK = 0
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            IK = IK + K
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 140
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)
               IKP1 = IK + K
               B(K+1) = B(K+1) + DDOT(K-1,AP(IKP1+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            IK = IK + K + K + 1
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
      SUBROUTINE DSPDI(AP,N,KPVT,DET,INERT,WORK,JOB)
      INTEGER N,JOB
      DOUBLE PRECISION AP(*),WORK(*)
      DOUBLE PRECISION DET(2)
      INTEGER KPVT(*),INERT(3)
C
C     DSPDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE
C     OF A DOUBLE PRECISION SYMMETRIC MATRIX USING THE FACTORS FROM
C     DSPFA, WHERE THE MATRIX IS STORED IN PACKED FORM.
C
C     ON ENTRY
C
C        AP      DOUBLE PRECISION (N*(N+1)/2)
C                THE OUTPUT FROM DSPFA.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX A.
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DSPFA.
C
C        WORK    DOUBLE PRECISION(N)
C                WORK VECTOR.  CONTENTS IGNORED.
C
C        JOB     INTEGER
C                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE
C                   IF  C .NE. 0, THE INVERSE IS COMPUTED,
C                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,
C                   IF  A .NE. 0, THE INERTIA IS COMPUTED.
C
C                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.
C
C     ON RETURN
C
C        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.
C
C        AP     CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF
C               THE ORIGINAL MATRIX, STORED IN PACKED FORM.
C               THE COLUMNS OF THE UPPER TRIANGLE ARE STORED
C               SEQUENTIALLY IN A ONE-DIMENSIONAL ARRAY.
C
C        DET    DOUBLE PRECISION(2)
C               DETERMINANT OF ORIGINAL MATRIX.
C               DETERMINANT = DET(1) * 10.0**DET(2)
C               WITH 1.0 .LE. DABS(DET(1)) .LT. 10.0
C               OR DET(1) = 0.0.
C
C        INERT  INTEGER(3)
C               THE INERTIA OF THE ORIGINAL MATRIX.
C               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.
C               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.
C               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INVERSE IS REQUESTED
C        AND  DSPCO  HAS SET RCOND .EQ. 0.0
C        OR  DSPFA  HAS SET  INFO .NE. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DCOPY,DDOT,DSWAP
C     FORTRAN DABS,IABS,MOD
C
C     INTERNAL VARIABLES.
C
      DOUBLE PRECISION AKKP1,DDOT,TEMP
      DOUBLE PRECISION TEN,D,T,AK,AKP1
      INTEGER IJ,IK,IKP1,IKS,J,JB,JK,JKP1
      INTEGER K,KK,KKP1,KM1,KS,KSJ,KSKP1,KSTEP
      LOGICAL NOINV,NODET,NOERT
C
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
      NOERT = MOD(JOB,1000)/100 .EQ. 0
C
      IF (NODET .AND. NOERT) GO TO 140
         IF (NOERT) GO TO 10
            INERT(1) = 0
            INERT(2) = 0
            INERT(3) = 0
   10    CONTINUE
         IF (NODET) GO TO 20
            DET(1) = 1.0D0
            DET(2) = 0.0D0
            TEN = 10.0D0
   20    CONTINUE
         T = 0.0D0
         IK = 0
         DO 130 K = 1, N
            KK = IK + K
            D = AP(KK)
C
C           CHECK IF 1 BY 1
C
            IF (KPVT(K) .GT. 0) GO TO 50
C
C              2 BY 2 BLOCK
C              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S)
C                      (S  C)
C              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
C              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
C
               IF (T .NE. 0.0D0) GO TO 30
                  IKP1 = IK + K
                  KKP1 = IKP1 + K
                  T = DABS(AP(KKP1))
                  D = (D/T)*AP(KKP1+1) - T
               GO TO 40
   30          CONTINUE
                  D = T
                  T = 0.0D0
   40          CONTINUE
   50       CONTINUE
C
            IF (NOERT) GO TO 60
               IF (D .GT. 0.0D0) INERT(1) = INERT(1) + 1
               IF (D .LT. 0.0D0) INERT(2) = INERT(2) + 1
               IF (D .EQ. 0.0D0) INERT(3) = INERT(3) + 1
   60       CONTINUE
C
            IF (NODET) GO TO 120
               DET(1) = D*DET(1)
               IF (DET(1) .EQ. 0.0D0) GO TO 110
   70             IF (DABS(DET(1)) .GE. 1.0D0) GO TO 80
                     DET(1) = TEN*DET(1)
                     DET(2) = DET(2) - 1.0D0
                  GO TO 70
   80             CONTINUE
   90             IF (DABS(DET(1)) .LT. TEN) GO TO 100
                     DET(1) = DET(1)/TEN
                     DET(2) = DET(2) + 1.0D0
                  GO TO 90
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
            IK = IK + K
  130    CONTINUE
  140 CONTINUE
C
C     COMPUTE INVERSE(A)
C
      IF (NOINV) GO TO 270
         K = 1
         IK = 0
  150    IF (K .GT. N) GO TO 260
            KM1 = K - 1
            KK = IK + K
            IKP1 = IK + K
            KKP1 = IKP1 + K
            IF (KPVT(K) .LT. 0) GO TO 180
C
C              1 BY 1
C
               AP(KK) = 1.0D0/AP(KK)
               IF (KM1 .LT. 1) GO TO 170
                  CALL DCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 160 J = 1, KM1
                     JK = IK + J
                     AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  160             CONTINUE
                  AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)
  170          CONTINUE
               KSTEP = 1
            GO TO 220
  180       CONTINUE
C
C              2 BY 2
C
               T = DABS(AP(KKP1))
               AK = AP(KK)/T
               AKP1 = AP(KKP1+1)/T
               AKKP1 = AP(KKP1)/T
               D = T*(AK*AKP1 - 1.0D0)
               AP(KK) = AKP1/D
               AP(KKP1+1) = AK/D
               AP(KKP1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 210
                  CALL DCOPY(KM1,AP(IKP1+1),1,WORK,1)
                  IJ = 0
                  DO 190 J = 1, KM1
                     JKP1 = IKP1 + J
                     AP(JKP1) = DDOT(J,AP(IJ+1),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IKP1+1),1)
                     IJ = IJ + J
  190             CONTINUE
                  AP(KKP1+1) = AP(KKP1+1)
     *                         + DDOT(KM1,WORK,1,AP(IKP1+1),1)
                  AP(KKP1) = AP(KKP1)
     *                       + DDOT(KM1,AP(IK+1),1,AP(IKP1+1),1)
                  CALL DCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 200 J = 1, KM1
                     JK = IK + J
                     AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  200             CONTINUE
                  AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)
  210          CONTINUE
               KSTEP = 2
  220       CONTINUE
C
C           SWAP
C
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 250
               IKS = (KS*(KS - 1))/2
               CALL DSWAP(KS,AP(IKS+1),1,AP(IK+1),1)
               KSJ = IK + KS
               DO 230 JB = KS, K
                  J = K + KS - JB
                  JK = IK + J
                  TEMP = AP(JK)
                  AP(JK) = AP(KSJ)
                  AP(KSJ) = TEMP
                  KSJ = KSJ - (J - 1)
  230          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 240
                  KSKP1 = IKP1 + KS
                  TEMP = AP(KSKP1)
                  AP(KSKP1) = AP(KKP1)
                  AP(KKP1) = TEMP
  240          CONTINUE
  250       CONTINUE
            IK = IK + K
            IF (KSTEP .EQ. 2) IK = IK + K + 1
            K = K + KSTEP
         GO TO 150
  260    CONTINUE
  270 CONTINUE
      RETURN
      END
#endif
