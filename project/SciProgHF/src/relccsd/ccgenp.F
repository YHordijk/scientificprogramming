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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION SUMT(A,N)
      REAL*8 SUMT,A(*)
      INTEGER N

#include "param.inc"

      SUMT=D0
      DO I=1,N
        SUMT=SUMT+A(I)
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FLIN(A,IDIM,IN,IM,DET)
      implicit none
C
C     LINEAR SIMULTANEOUS EQUATION
C
C     A(IN*IN) * X(IN*IM) = B(IN*IM)
C
C     A & B SHOULD BE STORED ON A(IN*(IN+IM))
C     SOLUTION X WILL BE STORED ON B PART IN DIMENSION A.
C
      integer idim, in, im, i, imax, j, k, m, n, npi, nr,jj, jmax
      real*8 a, aval, d, det, sign, temp, dott, stor
      DIMENSION A(IDIM,*)
      REAL ZERO, ONE
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      N=IN
      NR=IM
      JMAX=N+NR
      SIGN=ONE
C M IS THE STAGE OF ELIMINATION
      DO 49 M=1,N
      TEMP=ZERO
      DO 41 I=M,N
      IF(M.GT.1)A(I,M)=A(I,M)-DOTT(A(I,1),IDIM,A(1,M),1,M-1)
      AVAL=A(I,M)
      IF(ABS(AVAL).LE.TEMP)GOTO 41
      TEMP=ABS(AVAL)
      IMAX=I
 41   CONTINUE
      IF(TEMP.LE.ZERO)GOTO 999
      IF(IMAX.EQ.M)GOTO 45
      SIGN=-SIGN
      DO 44 J=1,JMAX
      STOR=A(M,J)
      A(M,J)=A(IMAX,J)
      A(IMAX,J)=STOR
 44   CONTINUE
 45   CONTINUE
      JJ=M+1
      IF(JJ.GT.JMAX)GOTO 49
      IF(M.GT.1)GOTO 47
      DO 46 J=JJ,JMAX
      A(1,J)=A(1,J)/A(1,1)
 46   CONTINUE
      D=A(1,1)
      GOTO 49
 47   CONTINUE
      DO 48 J=JJ,JMAX
      A(M,J)=(A(M,J)-DOTT(A(M,1),IDIM,A(1,J),1,M-1))/A(M,M)
 48   CONTINUE
      D=D*A(M,M)
 49   CONTINUE
      IF(NR.EQ.0) RETURN
      DO 59 I=1,NR
      NPI=N+I
      DO 58 K=2,N
      J=N+1-K
      A(J,NPI)=A(J,NPI)-DOTT(A(J,J+1),IDIM,A(J+1,NPI),1,K-1)
 58   CONTINUE
 59   CONTINUE
C***  IF(ABS(D).GE.1.0D+36) D=1.0D+36
C***  IF(ABS(D).LE.1.0D-36) D=1.0D-36
      DET=D*SIGN
      RETURN
C ON ZERO PIVOT, SET DET=0.AND RETURN TO CALLING PROGRAM NOV 1972
 999  DET=ZERO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION DOTT(A,NA,B,NB,N)
      implicit none
      real*8 a, b, dott
      integer nb,n,na, i, iapt, ibpt
      DIMENSION A(*),B(*)
      real*8 zero, one, d
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      IAPT=1
      IBPT=1
      D   =ZERO
      DO 10 I=1,N
        D   = D+A(IAPT)*B(IBPT)
        IAPT= IAPT+NA
        IBPT= IBPT+NB
 10   CONTINUE
      DOTT=D
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QUITC
      CALL QUIT ("Error in one of the C I/O c-routines")
      RETURN
      END
