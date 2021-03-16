!      Copyright (c) 2019 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org

C
C
Cjth - copied from dalton/sirius/sirnexkap.F
C
C  /* Deck ortvec */
      SUBROUTINE ORTVEC(NOLD,NNEW,NVDIM,THRLDP,VEC,LINDEP)
C
C 15-Mar-1985 hjaaj
C l.r. 4-May-1994 hjaaj (only elim. new vector if norm < THRLDP)
C
#include "implicit.h"
      DIMENSION VEC(NVDIM,*)
      LOGICAL   LINDEP(*)
C
      PARAMETER (THRRND=1.D-2, THRTT=1.D-4, D1=1.0D0)
C
C     Used from common blocks:
C     PRIUNIT: LUPRI
C
#include "priunit.h"
C
      CALL QENTER('ORTVEC')
C
      IF (NNEW.LE.0) GO TO 9999
      TLINDP = SQRT(THRLDP)
C
C     Normalize NNEW new vectors VEC(*,NOLD+1) TO VEC(*,NOLD+NNEW)
C
      IVEC = NOLD
      DO 200 INEW = 1,NNEW
         IVEC = IVEC + 1
         TT   = DNORM2(NVDIM,VEC(1,IVEC),1)
         IF (TT.LE.THRLDP) THEN
            LINDEP(INEW) = .TRUE.
            WRITE (LUPRI,8100) INEW,TT
         ELSE
            LINDEP(INEW) = .FALSE.
            IF (TT.LT.THRTT) THEN
               CALL DSCAL (NVDIM,(D1/TT),VEC(1,IVEC),1)
               TT = DNORM2(NVDIM,VEC(1,IVEC),1)
            END IF
            CALL DSCAL (NVDIM,(D1/TT),VEC(1,IVEC),1)
         END IF
  200 CONTINUE
C
C
C
      IROUND = 0
      ITURN  = 0
 1500 ITURN  = ITURN + 1
C
C     Orthogonalize new vectors against previous vectors
C
      DO 2000 K=1,NOLD
         DO 1900 J = NOLD+1,NOLD+NNEW
            TT = - DDOT(NVDIM,VEC(1,K),1,VEC(1,J),1)
            CALL DAXPY(NVDIM,TT,VEC(1,K),1,VEC(1,J),1)
 1900    CONTINUE
 2000 CONTINUE
C
C     Orthogonalize new vectors against each other
C     and normalization.
C
      DO 2400 INEW = 1,NNEW
         IF (.NOT.LINDEP(INEW)) THEN
C           ... orthogonalize using prev. vectors are normalized
            IVEC = NOLD + INEW
            DO 2300 JNEW = 1,(INEW-1)
               IF (.NOT.LINDEP(JNEW)) THEN
                  JVEC = NOLD + JNEW
                  TT = - DDOT(NVDIM,VEC(1,JVEC),1,VEC(1,IVEC),1)
                  CALL DAXPY(NVDIM,TT,VEC(1,JVEC),1,VEC(1,IVEC),1)
               END IF
 2300       CONTINUE
            TT = DNORM2(NVDIM,VEC(1,IVEC),1)
            IF (TT .LE. THRLDP) THEN
               LINDEP(INEW) = .TRUE.
               WRITE (LUPRI,8100) INEW,TT
            ELSE
               IF (TT .LT. THRRND) IROUND = IROUND+1
               IF (TT .LT. THRTT) THEN
                  CALL DSCAL(NVDIM,(D1/TT),VEC(1,IVEC),1)
                  TT = DNORM2(NVDIM,VEC(1,IVEC),1)
               END IF
               CALL DSCAL(NVDIM,(D1/TT),VEC(1,IVEC),1)
            END IF
         END IF
 2400 CONTINUE
C
C
      IF (IROUND.GT.0 .AND. ITURN.EQ.1) GO TO 1500
C
C
      JNEW = 0
      DO 4400 INEW = 1,NNEW
         IF (.NOT.LINDEP(INEW)) THEN
            JNEW = JNEW + 1
            IF (JNEW .LT. INEW) THEN
               CALL DCOPY (NVDIM,VEC(1,NOLD+INEW),1,VEC(1,NOLD+JNEW),1)
            END IF
         END IF
 4400 CONTINUE
      IF (JNEW .LT. NNEW) THEN
         WRITE (LUPRI,8200) NNEW,JNEW
         NNEW = JNEW
      END IF
C
C
C
 8100 FORMAT(/' ORTVEC: New vector no.',I3,
     *        ' is removed because of linear dependence;'
     *       /' norm of vector after Gram-Schmidt''s orthogonalization',
     *        1P,D9.2)
 8200 FORMAT(/' (ORTVEC) NNEW reduced from',I3,' to',I3)
C
C *** End of subroutine ORTVEC
C
 9999 CALL QEXIT('ORTVEC')
      RETURN
      END

