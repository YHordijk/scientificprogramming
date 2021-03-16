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

C  /* Deck ran1 */
      FUNCTION RAN1(IDUM)
#include "implicit.h"
      PARAMETER (IA=16807,IM=2147483647,AM=1.0D0/IM,
     +           IQ=127773,IR=2386,NTAB=32,NDIV=1+(IM-1)/NTAB,
     +           EPS=1.0D-7,RNMX=1.0D0-EPS)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     "Minimal" random number generator of Park and Miller with
C     Bays-Durham shuffle and added safeguards. Return a uniform random
C     deviate between 0.0 and 1.0 (exclusive of the end point values).
C     Call with IDUM a negative integer to initialize; thereafter,
C     do not alter IDUM between succesive deviates in a sequence. RNMX
C     should approximate the largest floating value that is less than 1.
C
C     References:
C       Press,Teukolsky,Vetterling & Flannery: 
C         "Numerical Recipes in FORTRAN",
C         Cambridge University Press 2nd edition 1992
C       Park & Miller (1988), Communications of the ACM, vol.31,pp.1192-1201
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      INTEGER IV(NTAB)
      SAVE IV,IY
      DATA IV/NTAB*0/, IY /0/
C
C     Initialization:
C        Be sure to prevent IDUM = 0
C        Load the shuffle table(after 8 warm-ups)
C
      IF(IDUM.LE.0.OR.IY.EQ.0) THEN
        IDUM = MAX(-IDUM,1)
        DO 10 J = NTAB+8,1,-1
          K = IDUM/IQ
          IDUM = IA*(IDUM-K*IQ)-IR*K
          IF(IDUM.LT.0) IDUM = IDUM+IM
          IF(J.LE.NTAB) IV(J) = IDUM
   10   CONTINUE
        IY = IV(1)
      ENDIF
C
C     Start here when not initializaing
C
      K = IDUM/IQ
C
C     Compute IDUM = MOD(IA*IDUM,IM) without overflows by
C     Schrage's method
C
      IDUM = IA*(IDUM-K*IQ)-IR*K
      IF(IDUM.LT.0) IDUM = IDUM+IM
      J = 1+IY/NDIV
      IY = IV(J)
      IV(J) = IDUM
      RAN1 = MIN(AM*IY,RNMX)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ranelm */
      FUNCTION RANELM(RMIN,RMAX,IDUM)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Pick a random number in the range <RMIN,RMAX>.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
      R = RAN1(IDUM)
      RANELM = RMIN + (RMAX-RMIN)*R
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ranmat */
      SUBROUTINE RANMAT(TYP,A,LDA,N,RMIN,RMAX,IDUM)
C**********************************************************************
C
C	Generation of test matrices with random elements with
C	values in the range <RMIN,RMAX>.
C	  TYP = S	Symmetric matrix	
C	  TYP = A	Antisymmetric matrix
C	  TYP = G	General matrix
C	On initialization of random number generation IDUM should
C	be a negative number
C**********************************************************************
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
      CHARACTER TYP*1
      DIMENSION A(LDA,N)
      IF(TYP.EQ.'G') THEN
        DO 10 J = 1,N
          DO 20 I = 1,N
            A(I,J) = RANELM(RMIN,RMAX,IDUM)
   20     CONTINUE
   10   CONTINUE
      ELSEIF(TYP.EQ.'S') THEN
        DO 30 J = 1,N
          DO 40 I = (J+1),N
            A(I,J) = RANELM(RMIN,RMAX,IDUM)
            A(J,I) = A(I,J)
   40     CONTINUE
          A(J,J) = RANELM(RMIN,RMAX,IDUM)
   30   CONTINUE
      ELSEIF(TYP.EQ.'A') THEN
        DO 50 J = 1,N
          DO 60 I = (J+1),N
            A(I,J) = RANELM(RMIN,RMAX,IDUM)
            A(J,I) = -A(I,J)
   60     CONTINUE
        A(J,J) = D0
   50   CONTINUE
      ELSE
        CALL QUIT('RANMAT: UNKNOWN MATRIX TYPE !')
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fulmat */
      SUBROUTINE FULMAT(MATTYP,LDA,NDIM,A)
#include "implicit.h"
#include "priunit.h"
      CHARACTER MATTYP*1
      DIMENSION A(LDA,NDIM)
C*SYMM  : Real symmetric matrix
C*ASYM  : Real antisymmetric matrix
C
C   ON INPUT: A - lower triangular matrix
C   ON OUTPUT:A - full matrix
      IF(MATTYP.EQ.'S') THEN
      DO 10 J = 1,NDIM
        DO 20 I = J,NDIM
          A(I,J) = A(I,J) + A(J,I)
          A(J,I) = A(I,J)
   20     CONTINUE
   10   CONTINUE
      ELSEIF(MATTYP.EQ.'A') THEN
      DO 30 J = 1,NDIM
        DO 40 I = J,NDIM
          A(I,J) = A(I,J) - A(J,I)
          A(J,I) = -A(I,J)
   40     CONTINUE
   30   CONTINUE
      ELSE
        WRITE(LUPRI,'(A,A6,A)') '** ERROR in FULMAT **  : Keyword ',
     +    MATTYP,' not recognized!'
      ENDIF
      END
C  /* Deck fulmat */
      SUBROUTINE FULMAT2(MATTYP,LDA,NDIM,A)
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DFAC=0.5D0)
      CHARACTER MATTYP*1
      DIMENSION A(LDA,NDIM)
C*SYMM  : Real symmetric matrix
C*ASYM  : Real antisymmetric matrix
C
C   ON INPUT: A - lower triangular matrix
C   ON OUTPUT:A - full matrix
      IF(MATTYP.EQ.'S') THEN
      DO 10 J = 1,NDIM
        DO 20 I = J,NDIM
          A(I,J) = A(I,J) + A(J,I)
          A(J,I) = A(I,J)
   20     CONTINUE
   10   CONTINUE
      ELSEIF(MATTYP.EQ.'A') THEN
      DO 30 J = 1,NDIM
        DO 40 I = J,NDIM
          A(I,J) = A(I,J) - A(J,I)
          A(J,I) = -A(I,J)
   40     CONTINUE
   30   CONTINUE
      ELSE
        WRITE(LUPRI,'(A,A6,A)') '** ERROR in FULMAT **  : Keyword ',
     +    MATTYP,' not recognized!'
      ENDIF
      NN=LDA*NDIM
      CALL DSCAL(NN,DFAC,A,1)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck offtst */
      FUNCTION OFFTST(A,N,IZ)
C
C	Checks the absolute value of off-diagonal elements
C	of a matrix
C
#include "implicit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
      DIMENSION A(N,N,IZ)
      DIA = D0
      ELM = D0
      DIATST = D0
      OFFTST = D0
      DO 10 K = 1,IZ
        DO 20 J = 1,N
          DO 30 I = 1,N
            ELM = ELM + D1
            OFFTST = ((ELM-D1)*OFFTST + ABS(A(I,J,K)))/ELM
   30     CONTINUE
   20   CONTINUE
        DO 40 J = 1,N
          DIA = DIA + D1
          DIATST = ((DIA-D1)*DIATST + ABS(A(J,J,K)))/DIA
   40   CONTINUE
   10 CONTINUE
      OFFTST = (ELM*OFFTST - DIA*DIATST)/(ELM-DIA)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mxform */
      FUNCTION MXFORM(X,NDIG)
C***********************************************************************
C
C	Find the format that gives optimum precision for real
C	variable X
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
      CHARACTER MXFORM*6
C
      MXFORM = '      '
      IF(X.EQ.D0) THEN
        NDEC = NDIG
      ELSE
        NDEC = INT(LOG10(ABS(X)))
      ENDIF
      IF((NDEC.LT.0).OR.(NDEC.GE.(NDIG-3))) THEN
        MXFORM(1:1) = 'E'
        NDEC = NDIG - 7
      ELSE
        MXFORM(1:1) = 'F'
        NDEC = NDIG - 3 - NDEC
      ENDIF
      IF(NDEC.LE.0) THEN
        WRITE(LUPRI,*) 'MXFORM out of boundaries...'
        WRITE(LUPRI,*) 'X=',X
        WRITE(LUPRI,*) 'NDEC,NDIG',NDEC,NDIG
        CALL QUIT('MXFORM: NDEC.LE.0')
      ENDIF
      IF(NDIG.LT.10) THEN
        IDIG = ICHAR('0') + NDIG
        MXFORM(2:2) = CHAR(IDIG)
        IOFF= 0
      ELSE
        IDIG = ICHAR('0') + INT(NDIG/10)
        MXFORM(2:2) = CHAR(IDIG)
        IDIG = ICHAR('0') + MOD(NDIG,10)
        MXFORM(3:3) = CHAR(IDIG)
        IOFF = 1
      ENDIF
      IND = 3+IOFF
      MXFORM(IND:IND) = '.'
      IF(NDEC.LT.10) THEN
        IDIG = ICHAR('0') + NDEC
        IND  = 4+IOFF
        MXFORM(IND:IND) = CHAR(IDIG)
        IOFF = IOFF+ 0
      ELSE
        IDIG = ICHAR('0') + 1
        IND = 4 + IOFF
        MXFORM(IND:IND) = CHAR(IDIG)
        IND = IND + 1
        IDIG = ICHAR('0') + MOD(NDEC,10)
        MXFORM(IND:IND) = CHAR(IDIG)
        IOFF = IOFF + 1
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck opnfil */
      SUBROUTINE OPNFIL(IUNIT,NAME,STATUS,SUB)
C***********************************************************************
C
C	Open sequential, unformatted file.
C	Based on OPENDX
C       T.Saue, University of Tromsoe, Norway, Feb.-94
C       Revised Aug. 2015 hjaaj to find IUNIT if not OK number in input
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) NAME,STATUS,SUB
      LOGICAL       IUNIT_USED

      IF (IUNIT .LE. 0) THEN
         IUNIT = 21
    5    CONTINUE
            INQUIRE( UNIT=IUNIT, OPENED=IUNIT_USED )
            IF (IUNIT_USED) THEN
               IUNIT = IUNIT + 1
               GO TO 5
            END IF
      END IF
C
      IF (STATUS .EQ. 'NEW') GO TO 10
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 100
C
C      Open old file
C
      OPEN(IUNIT,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',
     +     ACCESS='SEQUENTIAL',ERR=10)
C
C     Just in case....
      REWIND IUNIT
      GOTO 20
   10 CONTINUE
      IF(STATUS.EQ.'OLD') GOTO 110
C
C      Open new file
C
      OPEN(IUNIT,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     +     ACCESS='SEQUENTIAL')
   20 CONTINUE
      RETURN
C
C      Error messages
C
  100 CONTINUE
      WRITE(LUPRI,'(6A,I2,A)') SUB,': Invalid STATUS keyword ',STATUS,
     +     ' for file ',NAME,', (unit ',IUNIT,').'
      CALL QUIT('OPNFIL:Invalid STATUS keyword')
  110 CONTINUE
      WRITE(LUPRI,'(4A,I2,A)') SUB,': Old file ',NAME,' (unit',IUNIT,
     +    ') not found.'
      CALL QUIT('OPNFIL: Old file not found')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck putltr */
      SUBROUTINE PUTLTR(NM,N,ATRI,B)
C***********************************************************************
C
C	Copy elements from a rowwise lower triangularly packed matrix 
C       into a square matrix
C	
C	ATRI 	- rowwise lower triangularly packed matrix
C	B    	- square matrix
C	NM	- leading dimension of B
C	N	- order of matrix B
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION ATRI(N*(N+1)/2),B(NM,N)
C
      DO 10 I = 1,N
        IOFF = I*(I-1)/2
        DO 20 J = 1,I
          B(I,J) = ATRI(IOFF+J)
   20   CONTINUE
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck putsqr */
      SUBROUTINE PUTSQR(NM,NROW,NCOL,ASQR,B)
C***********************************************************************
C
C	Copy elements from columnwise array into rectangular matrix
C	
C	ASQR 	- rowwise lower triangularly packed matrix
C	B    	- rectangular matrix
C	NM	- leading dimension of B
C	NROW	- number of rows in B
C	NCOL    - number of columns in B
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION ASQR(NROW*NCOL),B(NM,NCOL)
C
      DO 10 J = 1,NCOL
        IOFF = (J-1)*NROW
        DO 20 I = 1,NROW
          B(I,J) = ASQR(IOFF+I)
   20   CONTINUE
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck unfld */
      SUBROUTINE UNFLD(NDIM,AR,AI)
C**********************************************************************
C	A Kramer symmetric matrix has the following structure:
C	
C  		-------------
C		|  A  |  B  |
C		-------------
C               | -B* |  A* |
C		-------------
C
C	in which A is a complex Hermitian matrix and B a complex 
C       antisymmetric matrix.
C
C	This subroutines assumes on input that only the lower triangle 
C       of blocks A and B are given, abd will on output give the full 
C       blocks A and B.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION AR(NDIM,NDIM/2),AI(NDIM,NDIM/2)
      IF(MOD(NDIM,2).EQ.1) THEN
        WRITE(LUPRI,'(A,I10,A)') 'Error in UNFLD: NDIM = ',NDIM,
     +  ' is not an even number.'
        CALL QUIT('Uneven dimension in UNFLD')
      ENDIF
      NHALF = NDIM/2
      IOFF = NHALF + 1
C** Unfold real part:
      CALL UNFLDR('SYMM',NDIM,NHALF,AR)
      CALL UNFLDR('ASYM',NDIM,NHALF,AR(IOFF,1))
C** Unfold imaginary part:
      CALL UNFLDR('ASYM',NDIM,NHALF,AI)
      CALL UNFLDR('ASYM',NDIM,NHALF,AI(IOFF,1))
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck unfldr */
      SUBROUTINE UNFLDR(MATTYP,LDA,NDIM,A)
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
      CHARACTER*4 MATTYP
      DIMENSION A(LDA,NDIM)
C*SYMM  : Real symmetric matrix
C*ASYM  : Real antisymmetric matrix
C
C   ON INPUT: A - lower triangular matrix
C   ON OUTPUT:A - full matrix
      IF(MATTYP.EQ.'SYMM') THEN
      DO 10 J = 1,NDIM
        DO 20 I = (J+1),NDIM
          A(J,I) = A(I,J)
   20     CONTINUE
   10   CONTINUE
      ELSEIF(MATTYP.EQ.'ASYM') THEN
      DO 30 J = 1,NDIM
        DO 40 I = (J+1),NDIM
          A(J,I) = -A(I,J)
   40     CONTINUE
        A(J,J) = D0
   30   CONTINUE
      ELSE
        WRITE(LUPRI,'(A,A6,A)') '** ERROR in UNFLD **  : Keyword ',
     +    MATTYP,' not recognized!'
      ENDIF
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sectid */
      CHARACTER*12 FUNCTION SECTID(TID)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Given time in seconds this function return a character string
C     with time in suitable units
C     Revised 4-aug-2017 hjaaj: input TID is not modified any more.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      IMPLICIT NONE
      REAL*8   TID, TID_rest
      INTEGER  ITID, IMIN, IHOUR, IDAY
C
      ITID = TID
      IF(ITID.GE.86400) GOTO 3
      IF(ITID.GE.3600)  GOTO 2
      IF(ITID.GE.60)    GOTO 1
C
C     Format : Seconds
C
      IF (TID .LT. 0.0D0) THEN
         WRITE(SECTID,'(F11.3,A1)') TID,'s'
      ELSE
         WRITE(SECTID,'(F11.8,A1)') TID,'s'
      END IF
      RETURN
C
C     Format: Minutes and seconds
C
    1 CONTINUE
      IMIN = ITID/60
      TID_rest  = TID - IMIN*60
      WRITE(SECTID,'(I2,A3,F6.3,A1)') IMIN,'min',TID_rest,'s'
      RETURN
C
C     Format: Hours, minutes and seconds
C
    2 CONTINUE
      IHOUR = ITID/3600
      ITID  = MOD(ITID,3600)
      IMIN  = ITID/60
      ITID  = 3600*IHOUR + 60*IMIN
      ITID   = TID - ITID
      WRITE(SECTID,'(I3,A1,I2,A3,I2,A1)')
     &     IHOUR,'h',IMIN,'min',ITID,'s'
      RETURN
C
C     Format: Days, hours and minutes
C
    3 CONTINUE
      IDAY  = ITID/86400
      ITID  = MOD(ITID,86400)
      IHOUR = ITID/3600
      ITID  = MOD(ITID,3600)
      IMIN  = ITID/60
      IF (IDAY .LE. 999) THEN
         WRITE(SECTID,'(I3,A1,I2,A1,I2,A3)')
     &      IDAY,'d',IHOUR,'h',IMIN,'min'
      ELSE
         WRITE(SECTID,'(I5,A1,I2,A1,I2,A1)')
     &      IDAY,'d',IHOUR,'h',IMIN,'m'
      END IF
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck logset */
      SUBROUTINE LOGSET(N,LVAL,LVEC)
C***********************************************************************
C     This routine will initialize a logical array
C     Written by T.Saue , January 1995
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL LVAL,LVEC(*)
      DO 10 I = 1,N
        LVEC(I) = LVAL
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lbit */
C##############lbit##############################################
      LOGICAL FUNCTION LBIT(I,N)
C	Based on an analogous DISCO routine(Jan Almloef)
#include "implicit.h"
#include "ibtfun.h"
      LBIT = .FALSE.
      IF (IBTAND(IBTSHR(I,N-1),1).EQ.1) LBIT = .TRUE.
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lset */
      SUBROUTINE LSET(N,LVAL,LVEC)
#include "implicit.h"
      LOGICAL LVEC(N),LVAL
      DO 10 I = 1,N
        LVEC(I) = LVAL
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if defined (SYS_CRAY) || defined (SYS_T3D) 
CSK|| defined (INT_STAR8)
      REAL FUNCTION GAMMP(A,X)
#else
      DOUBLE PRECISION FUNCTION GAMMP(A,X)
#endif
#include "implicit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C****************************************************************C
C     SUBFUNCTION TO EVALUATE THE NORMALIZED INCOMPLETE GAMMA    C
C     FUNCTION USING SERIES EXPANSION FOR SMALL X AND CONTINUED  C
C     FRACTION REPRESENTATION FOR 'LARGE' X (I.E.  X > A+1)      C
C                                                                C
C     FROM THE NUMERICAL RECIPES LIBRARY, BUT MODIFIED BY AN     C
C     ADDITIONAL FACTOR OF GAMMA(A) SO THAT GAMMAP RETURNS THE   C
C     UN-NORMALIZED GAMMA-FUNCTION (6.5.2), AND NOT THE          C
C     P-FUNCTION (6.5.1). EQUATION NUMBERS AND DEFINITIONS ARE   C
C     FROM ABRAMOWITZ AND STEGUN                                 C
C                                                                C
C****************************************************************C
      IF(X.LT.D0.OR.A.LE.D0) THEN
         GAMMP=D0
         CALL QUIT('GAMMP: Negative A or X..')
      ELSEIF(X.LT.A+D1) THEN
         CALL GSER(GAMMP,A,X,GLN)
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP = D1 - GAMMCF
      ENDIF
C
C     RETURN THE UNNORMALIZED INCOMPLETE GAMMA FUNCTION
C
      GAMMP = GAMMP * EXP(GLN)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GSER(GAMSER,A,X,GLN)
#include "implicit.h"
C******************************************************************C
C     SERIES EXPANSION REPRESENTATION                              C
C******************************************************************C
      PARAMETER (ITMAX=100,EPS=1.D-12,D1 = 1.0D0)
      GLN=GAMMLN(A)
      AP=A
      SUM=D1/A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+D1
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      CALL QUIT('GSER: ITMAX Too Small')
C
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
#include "implicit.h"
C*******************************************************************C
C     CONTINUED FRACTION REPRESENTATION                             C
C*******************************************************************C
      PARAMETER (ITMAX=100,EPS=1.D-12,D0 = 0.0D0,D1 = 1.0D0)
      GLN=GAMMLN(A)
      GOLD=D0
      A0=D1
      A1=X
      B0=D0
      B1=D1
      FAC=D1
      DO 11 N=1,ITMAX
        AN=dble(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.D0)THEN
          FAC=D1/A1
          G=B1*FAC
          IF(DABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      CALL QUIT('GCF: ITMAX Too Small')
C
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*G
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DFACTR
#include "implicit.h"
      PARAMETER(D1 = 1.0D0)
C***************************************************************C
C     DFACTR INITIALIZES A VECTOR OF DOUBLE-FACTORIAL FUNCTIONS C
C     SUBJECT TO THE CONVENTIONS                                C
C                                                               C
C     DFACT(I) = (I-2)!!                                        C
C     (-1)!! = 1                                                C
C     0!!    = 1                                                C
C     N!!    = N*(N-2)!!                                        C
C                                                               C
C***************************************************************C
      COMMON/FACTRL/DFACT(30)
C
      DFACT(1)=D1
      DFACT(2)=D1
      DO 10 M=3,30
      DFACT(M)=dble(M-2)*DFACT(M-2)
10    CONTINUE
      RETURN
      END
      subroutine init_factorial_common_block
!     radovan: adaptation of dfactr for usual factorials

      implicit none

      real(8) :: fact
      integer :: i

      common /factrl/ fact(30)

      fact(1) = 1.0d0
     
      do i = 2, 30
        fact(i) = float(i)*fact(i - 1)
      end do

      end subroutine



C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DAMATR(N,AMAT,LDM,ATRI)
C***********************************************************************
C
C     Row-pack lower triangle of matrix AMAT into ATRI
C
C     Written by T.Saue - Aug 16 1996
C     Last revision: Aug 16 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
      DIMENSION AMAT(LDM,N), ATRI(*)
C
#include "nosync.h"
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J
            ATRI(JOFF+I) = AMAT(I,J)
  100    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dunit2 */
      SUBROUTINE DUNIT2(A,N,LR,LC,NZ)
C***********************************************************************
C
C     Subroutine DUNIT2 sets a general square matrix A equal
C     to a unit matrix.
C     /VER 2/ 14-Sep-1983 hjaaj
C     Include leading dimension - T.Saue Aug 16 1996
C
C***********************************************************************
#include "implicit.h"
      DIMENSION A(LR,LC,NZ)
      PARAMETER (D1=1.0D00, D0=0.0D00)
C
      IF (N.GT.LR .OR. N.GT.LC) THEN
        CALL QUIT('Inconsistent N,LR,LC in DUNIT2 call')
      ENDIF
      DO IZ = 1,NZ
        DO J = 1,N
          DO I = 1,N
            A(I,J,IZ) = D0
          ENDDO
        ENDDO
      ENDDO
      DO J = 1,N
        A(J,J,1) = D1
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck jaco2 */
      SUBROUTINE JACO2(F,V,NB,NMAX,NROWV,NEVEC,LDM,BIG,JBIG)
C
C Revisions:
C
C  26-Aug-1984 tsaue  Introduced a new parameter LDM, which
C                     is the leading dimension of V and must
C                     be larger or equal to NROWV.
C   2-Nov-1984 hjaaj (new parameter NROWV such that
C                     dim(V) = (NROWV,NMAX). This makes
C                     it possible to solve eigenproblem
C                     in a reduced basis but get the
C                     eigenvectors in the original full
C                     basis, e.g. less mo's than ao's)
C  23-Feb-1989 hjaaj  Note that if NROWV = 0 then only
C                     eigenvalues will be calculated,
C                     V matrix will not be referenced.
C  27-Jul-1990 hjaaj  Changed -CX,+SSX transformation to +CX,-SSX
C                     transformation; probably the -CX,+SSX
C                     transformation was responsible for that
C                     the eigenvectors easily changed sign.
C                     Changed initial test on NB. Changed SD.
C                     Optimized IR loop.
C     Jun-1992 ov     Parameters for 0.5, 1.5, ... (for Cray)
C  20-Jul-1992 hjaaj  Changed C1,C2 to THRZER
C  30-oct-1992 hjaaj  zero f(ab) to avoid round-off errors
C                     absolute conv.threshold SD=C1
C
C jun 2005,MI : Added NEVEC to decide whether to calculate eigenvectors
C
#include "implicit.h"
      DIMENSION F(*),V(*)
      DIMENSION BIG(*) ,JBIG(*)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, ROOT2 = 0.707106781186548D0)
      PARAMETER(DP5 = 0.5D0, D1P5 = 1.5D0, D1P375 = 1.375D0,
     *          D3P875 = 3.875D0, DP25 = 0.25D0)
#include "thrzer.h"
C      DATA C1,C2,C3,C4,C5,C6/THRZER,THRZER,1.D-20,1.D-14,1.D-9,1.D-5/
COLD   DATA C1,C2,C3,C4,C5,C6/1.D-12,1.D-12,1.D-20,1.D-14,1.D-9,1.D-5/
      DATA C1,C2,C3,C4,C5,C6/1.D-15,1.D-15,1.D-20,1.D-14,1.D-9,1.D-5/
      IF (NB.LE.1 .OR. NMAX.LE.0) RETURN
Cwas: IF (NB.EQ.1) RETURN !900727-hjaaj
      IF(LDM.LT.NROWV) CALL QUIT('JACO2: LDM.LT.NROWV')
      DO 190 I=1,NB
         JBIGI=0
         J=MIN(I-1,NMAX)
         IF (J .GT. 0) THEN
            II = (I*I-I)/2
            ABIGI=D0
            DO 18 K=1,J
            IF (ABIGI .GE. ABS(F(II+K))) GO TO  18
               ABIGI=ABS(F(II+K))
               JBIGI=K
   18       CONTINUE
         END IF
         IF (JBIGI .GT. 0) THEN
            JBIG(I) = JBIGI
            BIG(I)  = F(II+JBIGI)
         ELSE
            JBIG(I) = 0
            BIG(I)  = D0
         END IF
  190 CONTINUE
C
#if defined (VAR_OLDCODE)
C 900727-hjaaj:
C SD calculation was done in every Jacobi iteration.
C Now the largest absolute element in F is found once and
C the SD based on that value is used in every iteration.
  410 SD=1.05D 00
      DO 220 J=1,NMAX
         DAB=ABS(F(J*(J+1)/2))
CHJ-861103: commented out next line, it seems to make the loop
C           meaningless (setting SD equal to J=NMAX value always!)
C        IF (SD .GT. DAB) SD=DAB
  220    SD=MAX(SD,DAB)
      SD=MAX(C1,C2*SD)
#else
C 921030-hjaaj: SD = C1 now
      NNB = (NB*NB+NB)/2
C     SD = 1.05D0
C     DO 220 J = 1,NNB
C        SD = MAX(SD, ABS(F(J)) )
C 220 CONTINUE
C     SD=MAX(C1,C2*SD)
      SD=C1
C
      MXITJA = 50*NNB
      ITJACO = 0
  410 ITJACO = ITJACO + 1
      IF (ITJACO .GT. MXITJA) THEN
         CALL QUIT('ERROR: JACO did not converge ...')
      END IF
#endif
      T = D0
      DO 230 I=2,NB
      IF (T .GE.  ABS(BIG(I))) GO TO 230
         T = ABS(BIG(I))
         IB= I
  230 CONTINUE
      IF(T.LT.SD) GO TO 420
         IA =JBIG(IB)
         IAA=IA*(IA-1)/2
         IBB=IB*(IB-1)/2
         DIF=F(IAA+IA)-F(IBB+IB)
         IF( ABS(DIF) .GT. C3) GO TO 271
            SSX=ROOT2
            CX=ROOT2
         GO TO 270
  271       T2X2 =BIG(IB)/DIF
            T2X25=T2X2*T2X2
         IF(T2X25 .GT. C4) GO TO 240
            CX=1.D 00
            SSX=T2X2
         GO TO 270
  240    IF(T2X25 .GT. C5) GO TO 250
            SSX=T2X2*(D1 - D1P5*T2X25)
            CX=D1 - DP5*T2X25
         GO TO 270
  250    IF(T2X25 . GT . C6) GO TO 260
            CX=D1+T2X25*(T2X25*D1P375 - DP5 )
            SSX= T2X2*(D1 + T2X25*(T2X25*D3P875 - D1P5))
         GO TO 270
  260       T = DP25  / SQRT(DP25   + T2X25)
            CX= SQRT(DP5   + T)
            SSX= SIGN( SQRT(DP5 - T),T2X2)
  270    CONTINUE
#include "ivdep.h"
         DO 275 IR=1,IA
            T        = F(IAA+IR)*SSX
            F(IAA+IR)= F(IAA+IR)*CX+F(IBB+IR)*SSX
            F(IBB+IR)=-T           +F(IBB+IR)*CX
  275    CONTINUE
         IEAA=IAA+IA
         IEAB=IBB+IA
         TT  =F(IEAB)
         F(IEAB)=BIG(IB)
         IF (JBIG(IA) .EQ. 0) THEN
            IRST = IA   + 1
            IEAR = IEAA + IA
            IEBR = IEAB + 1
         ELSE
            IRST = IA
            IEAR = IEAA
            IEBR = IEAB
         END IF
         DO 390 IR = IRST,NB
#if !defined (VAR_OLDCODE)
            IF (IR .EQ. IA) GO TO 360
C              ... we have checked above that JBIG(IA) .ne. 0
#else
            IF (IR .EQ. IA) THEN
               GO TO 360
C              ... we have checked above that JBIG(IA) .ne. 0
C              IF(JBIG(IR)) 360,380,360
            END IF
#endif
            T      = F(IEAR)*SSX
            F(IEAR)= F(IEAR)*CX+F(IEBR)*SSX
            F(IEBR)=-T         +F(IEBR)*CX
            T   =F(IEAR)
            IT  =IA
            IF(IR-IB) 340,310,320
  310          F(IEAA)=F(IEAA)*CX+F(IEAB)*SSX
C              921030+hjaaj: zero f(ab) to avoid round-off errors
C              F(IEAB)=     TT*CX+F(IEBR)*SSX
               F(IEAB)=     D0
               F(IEBR)=    -TT*SSX+F(IEBR)*CX
            GO TO 360
  320       IF(ABS(T) .GE.  ABS(F(IEBR))) GO TO 340
               T   =F(IEBR)
               IT  =IB
  340       IF(ABS(T) .LT.  ABS(BIG(IR))) GO TO 350
               BIG(IR)  = T
               JBIG(IR) = IT
            GO TO 380
  350       IF(IA .NE. JBIG(IR) .AND. IB .NE. JBIG(IR))  GO TO 380
  360          K= IEAR - IA
               JBIGI = 0
               IR1=MIN (IR-1,NMAX)
               IF (IR1 .GT. 0) THEN
                  ABIGI = D0
                  DO 370 I=1,IR1
                  IF(ABIGI .GE. ABS(F(K+I)))  GO TO 370
                     ABIGI = ABS(F(K+I))
                     JBIGI =I
  370             CONTINUE
               END IF
               IF (JBIGI .GT. 0) THEN
                  JBIG(IR) = JBIGI
                  BIG(IR)  = F(K+JBIGI)
               ELSE
                  JBIG(IR) = 0
                  BIG(IR)  = D0
               END IF
  380          CONTINUE
               IEAR = IEAR + IR
               IF (IR .GE. IB) THEN
                  IEBR = IEBR + IR
               ELSE
                  IEBR = IEBR + 1
               END IF
  390       CONTINUE
Cwas         JAA=(IA-1)*NROWV
Cwas         JBB=(IB-1)*NROWV
         JAA=(IA-1)*LDM
         JBB=(IB-1)*LDM
#include "ivdep.h"
CMI ... ensure that no mem-error happens when we do not want 
CMI     to calculate eigenvectors!!
       IF (NEVEC.GT.0) THEN
         DO 400 I=1,NROWV
            T=V(JBB+I)*SSX
            V(JBB+I)=-V(JAA+I)*SSX + V(JBB+I)*CX
  400       V(JAA+I)= V(JAA+I)*CX + T
       ENDIF
      GO TO 410
  420 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qorder3 */
      SUBROUTINE QORDER3(N,EVAL,QEVEC,LVR,LVC,NEVEC,NZ,IJOB)
C***********************************************************************
C
C Purpose: order the N values in EVAL and
C          their associated quaternion vectors in QEVEC:
C     IJOB.EQ. 2  -  absolute ascending order  
C                    ABS(EVAL(i+1)) .ge. ABS(EVAL(i))
C     IJOB.EQ. 1  -  ascending order   EVAL(i+1) .ge. EVAL(i)
C     IJOB.EQ.-1  -  descending order   EVAL(i+1) .le. EVAL(i)
C     IJOB.EQ.-2  -  absolute descending order  
C                    ABS(EVAL(i+1)).le. ABS(EVAL(i))
C     IJOB.EQ.0  -  return
C          
C          Generalization of routine ORDER3 by T.Saue Aug 1996
C          to quaternion eigenvectors by H.J.Aa.Jensen Feb 2008
C
C***********************************************************************
#include "implicit.h"
      DIMENSION EVAL(N), QEVEC(LVR,LVC,NZ)
C     active dimension: QEVEC(1:NEVEC,1:N,NZ)
      IF (N.LE.1.OR.IJOB.EQ.0) RETURN
      II = IJOB + 3
      GOTO (1,2,3,4,5) II
1     CONTINUE
C***********************************************************************
C IJOB.EQ.-2  -  absolute descending order 
C                ABS(EVAL(i+1)).le. ABS(EVAL(i))
C***********************************************************************
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (ABS(EVAL(J)) .GE. ABS(EMIN)) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)   =EMIN
            IF (NEVEC .GT. 0) THEN
               DO IZ = 1,NZ
                  CALL DSWAP(NEVEC,QEVEC(1,I,IZ),1,QEVEC(1,IMIN,IZ),1)
               END DO
            ENDIF
         ENDIF
      ENDDO
      RETURN
2     CONTINUE
C***********************************************************************
C      IJOB.EQ.-1  -  descending order   EVAL(i+1) .le. EVAL(i)
C***********************************************************************
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (EVAL(J) .GE. EMIN) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)=EMIN
            IF (NEVEC .GT. 0) THEN
               DO IZ = 1,NZ
                  CALL DSWAP(NEVEC,QEVEC(1,I,IZ),1,QEVEC(1,IMIN,IZ),1)
               END DO
            ENDIF
         ENDIF
      ENDDO
      RETURN
3     CONTINUE
C***********************************************************************
C   IJOB.EQ. 0  -  return
C***********************************************************************
      RETURN
4     CONTINUE
C***********************************************************************
C   IJOB.EQ. 1  -  ascending order   
C                  EVAL(i+1) .ge. EVAL(i)
C***********************************************************************
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (EVAL(J) .LT. EMIN) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)=EMIN
            IF (NEVEC .GT. 0) THEN
               DO IZ = 1,NZ
                  CALL DSWAP(NEVEC,QEVEC(1,I,IZ),1,QEVEC(1,IMIN,IZ),1)
               END DO
            ENDIF
         ENDIF
      ENDDO
      RETURN
5     CONTINUE
C***********************************************************************
C   IJOB.EQ. 2  -  absolute ascending order  
C                  ABS(EVAL(i+1)) .ge. ABS(EVAL(i))
C***********************************************************************
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (ABS(EVAL(J)).LT.ABS(EMIN)) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)=EMIN
            IF (NEVEC .GT. 0) THEN
               DO IZ = 1,NZ
                  CALL DSWAP(NEVEC,QEVEC(1,I,IZ),1,QEVEC(1,IMIN,IZ),1)
               END DO
            ENDIF
         ENDIF
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck order3 */
      SUBROUTINE ORDER3(EVEC,EVAL,LDM,N,NEVEC,IJOB)
C***********************************************************************
C
C Purpose: order the N values in EVAL and their associated vectors 
C     in EVEC:
C     IJOB.EQ. 2  -  absolute ascending order  
C                    ABS(EVAL(i+1)) .ge. ABS(EVAL(i))
C     IJOB.EQ. 1  -  ascending order   EVAL(i+1) .ge. EVAL(i)
C     IJOB.EQ.-1  -  descending order   EVAL(i+1) .le. EVAL(i)
C     IJOB.EQ.-2  -  absolute descending order  
C                    ABS(EVAL(i+1)).le. ABS(EVAL(i))
C     IJOB.EQ.0  -  return
C          
C          Generalization of routines ORDER/ORDER2
C          Written by T.Saue - Aug 26 1996
C          Last revision: Aug 26 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
      DIMENSION EVEC(*),EVAL(*)
      IF (N.LE.1.OR.IJOB.EQ.0) RETURN
      II = IJOB + 3
      GOTO (1,2,3,4,5) II
1     CONTINUE
C***********************************************************************
C IJOB.EQ.-2  -  absolute descending order 
C                ABS(EVAL(i+1)).le. ABS(EVAL(i))
C***********************************************************************
      IN = 1
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (ABS(EVAL(J)) .GE. ABS(EMIN)) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)   =EMIN
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*LDM+1),1)
            ENDIF
         ENDIF
         IN = IN + LDM
      ENDDO
      RETURN
2     CONTINUE
C***********************************************************************
C      IJOB.EQ.-1  -  descending order   EVAL(i+1) .le. EVAL(i)
C***********************************************************************
      IN = 1
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (EVAL(J) .GE. EMIN) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)=EMIN
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*LDM+1),1)
            ENDIF
         ENDIF
         IN = IN + LDM
      ENDDO
      RETURN
3     CONTINUE
C***********************************************************************
C   IJOB.EQ. 0  -  return
C***********************************************************************
      RETURN
4     CONTINUE
C***********************************************************************
C   IJOB.EQ. 1  -  ascending order   
C                  EVAL(i+1) .ge. EVAL(i)
C***********************************************************************
      IN = 1
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (EVAL(J) .LT. EMIN) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)=EMIN
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*LDM+1),1)
            ENDIF
         ENDIF
         IN = IN + LDM
      ENDDO
      RETURN
5     CONTINUE
C***********************************************************************
C   IJOB.EQ. 2  -  absolute ascending order  
C                  ABS(EVAL(i+1)) .ge. ABS(EVAL(i))
C***********************************************************************
      IN = 1
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (ABS(EVAL(J)).LT.ABS(EMIN)) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)=EMIN
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*LDM+1),1)
            ENDIF
         ENDIF
         IN = IN + LDM
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gatpak */
      SUBROUTINE GATPAK(NA,A,NB,B,IND,IOFF)
C***********************************************************************
C
C     This routine will gather elements of the same "kin" in the
C     rowpacked lower triangular matrix A and place them in the
C     rowpacked lower triangular matrix B
C
C     "Kinship" is determined by the integer array IND pointing from
C     elements of A to elements of B
C
C     Written by T.Saue  -  July 23 1994
C     LAST VERSION: July 23 1994
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION A(*),B(*),IND(*)
C
      KA = 0
      DO 10 IA = 1,NA
        IB = IND(IA) - IOFF
        DO 20 JA = 1,IA
          JB = IND(JA) - IOFF
          KA = KA + 1
          KB = (IB*(IB-1)/2)+JB
          B(KB) = B(KB) + A(KA)
   20   CONTINUE
   10 CONTINUE
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck igath */
      SUBROUTINE IGATH(IV1,N1,IV2,N2,N,IND)
C***********************************************************************
C
C     Reorder an integer array
C
C     Written by T..Saue November 1994
C     Last revision: Nov 16 1994 - tsaue
C***********************************************************************
#include "implicit.h"
      DIMENSION IV1(N1),IV2(N2),IND(*)
      DO 10 I = 1,N
        IV2(I) = IV1(IND(I))
   10 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck privec */
      SUBROUTINE PRIVEC(NAME,IVEC,NVEC)
C***********************************************************************
C
C     Print integer array with indices in compact form
C
C     Written by T.Saue, November 1994
C     Last revised: Nov 24 1994 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER NAME*6
      DIMENSION IVEC(NVEC)
      NROW = (NVEC+15)/16
      J1 = 1
      DO 10 J = 1,NROW
        J2 = MIN(J1+15,NVEC)
        WRITE(LUPRI,'(3X,A3,4X,20I6)') 'IND',(JJ, JJ = J1,J2)
        WRITE(LUPRI,'(3X,A6,1X,20I6)')
     &            NAME,(IVEC(JJ),JJ = J1,J2)
        J1 = J1 + 16
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtcomp */
      SUBROUTINE MTCOMP(COMP,AMAT,NDIM,NZ,CMAT)
C***********************************************************************
C
C     This subroutine will pick out symmetric('S') or antisymmetric('A')
C     component of a quadratic matrix AMAT and place it in CMAT.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DP5 = 0.5D0)
      DIMENSION AMAT(NDIM,NDIM,NZ),CMAT(NDIM,NDIM,NZ)
      CHARACTER COMP*1
      N2DIMQ = NDIM*NDIM*NZ
C
C     Symmetric component
C
      IF    (COMP.EQ.'S') THEN
        DO 10 IZ = 1,NZ
          DO 20 J = 1,NDIM
            DO 30 I = 1,NDIM
              CMAT(I,J,IZ) = AMAT(I,J,IZ)+AMAT(J,I,IZ)
 30         CONTINUE
 20       CONTINUE
 10     CONTINUE
        CALL DSCAL(N2DIMQ,DP5,CMAT,1)
      ELSEIF(COMP.EQ.'A') THEN
        DO 40 IZ = 1,NZ
          DO 50 J = 1,NDIM
            DO 60 I = 1,NDIM
              CMAT(I,J,IZ) = AMAT(I,J,IZ)-AMAT(J,I,IZ)
 60         CONTINUE
 50       CONTINUE
 40     CONTINUE
        CALL DSCAL(N2DIMQ,DP5,CMAT,1)
      ELSE
        CALL QUIT('MTCOMP:Unknown component '//COMP)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck difmmh */
      SUBROUTINE DIFMMH(AMAT,NDIM,NZ,LRA,LCA)
C***********************************************************************
C
C     Take difference between a general quaternionic matrix AMAT and
C     its Hermitian conjugate; return in AMAT. Note that this gives
C     an anti-Hermitian result !
C
C     Written by T.Saue, August 23 1995
C     Last revision: tsaue - aug23-95
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION AMAT(LRA,LCA,NZ)
C
C     Take difference
C
      DO J = 1,NDIM
        DO I = 1,(J-1)
          AMAT(I,J,1) = AMAT(I,J,1)-AMAT(J,I,1)
          AMAT(J,I,1) = -AMAT(I,J,1)
        ENDDO
        AMAT(J,J,1) = D0
      ENDDO
      DO IZ = 2,NZ
        DO J = 1,NDIM
          DO I = 1,(J-1)
            AMAT(I,J,IZ) = AMAT(I,J,IZ)+AMAT(J,I,IZ)
            AMAT(J,I,IZ) = AMAT(I,J,IZ)
          ENDDO
          AMAT(J,J,IZ) = AMAT(J,J,IZ)+AMAT(J,J,IZ)
        ENDDO
      ENDDO
C
      RETURN
      END

      SUBROUTINE QGETHE(AMAT,NDIM,NZ,LRA,LCA)
C***********************************************************************
C
C     Quaternian GEneral matrix To HErmitian matrix:
C     generate the Hermitian component of the general matrix AMAT;
C     return in AMAT.
C
C     Written by H.J.Aa. Jensen, June 2002
C     (based on DIFMMH)
C
C***********************************************************************
#include "implicit.h"
      PARAMETER(D0=0.0D0, DP5=0.5D0)
      DIMENSION AMAT(LRA,LCA,NZ)
C
C     Generate A + A(dagger)
C
      DO J = 1,NDIM
        DO I = 1,(J-1)
          AMAT(I,J,1) = DP5*(AMAT(I,J,1)+AMAT(J,I,1))
          AMAT(J,I,1) = AMAT(I,J,1)
        ENDDO
      ENDDO
      DO IZ = 2,NZ
        DO J = 1,NDIM
          DO I = 1,(J-1)
            AMAT(I,J,IZ) = DP5*(AMAT(I,J,IZ)-AMAT(J,I,IZ))
            AMAT(J,I,IZ) = -AMAT(I,J,IZ)
          ENDDO
          AMAT(J,J,IZ) = D0
        ENDDO
      ENDDO
C
      RETURN
      END

      SUBROUTINE QGETHCGE(AMAT,NDIM,NZ,LRA,LCA)
C*****************************************************************************
C
C     Quaternion GEneral matrix To 
C                HErmitian Conjugate GEneral matrix:
C
C     generate the Hermitian conjugate matrix 
C     of the general matrix AMAT;
C
C     return the new general matrix in AMAT 
C       (so AMAT is overwritten !)
C
C     Written by Miro Ilias, July 2003
C     (based on QGETHE)
C
C*****************************************************************************
#include "implicit.h"
      DIMENSION AMAT(LRA,LCA,NZ)
C
C  Generate A(dagger), where A(dagger)ij = Aji^*
C
      DO J = 1,NDIM
        DO I = 1, (J-1)
          AMAT(J,I,1) = AMAT(I,J,1)
        ENDDO
      ENDDO
      DO IZ = 2,NZ
        DO J = 1,NDIM
          DO I = 1,(J-1)
            AMAT(J,I,IZ) = -AMAT(I,J,IZ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END

      SUBROUTINE QAHM(AMAT,NDIM,NZ,LRA,LCA)
C*****************************************************************************
C
C   Make -A^(dagger) from A general quanternion matrix
C
C   Written by Miro Ilias, August 2005
C
C*****************************************************************************
#include "implicit.h"
      DIMENSION AMAT(LRA,LCA,NZ)
C
C  Generate -A(dagger), where A(dagger)ij = Aji^*
C
      DO J = 1, NDIM
        DO I = 1, (J-1)
          ZP  =  AMAT(J,I,1)
          AMAT(J,I,1) = -AMAT(I,J,1)
          AMAT(I,J,1) = -ZP
        ENDDO
        AMAT(J,J,1)=-AMAT(J,J,1)
      ENDDO
C  ... quaternion imaginary parts - 
C  -   only transpose the matrix and preserve the signs
      DO IZ = 2, NZ
        DO J = 1,NDIM
          DO I = 1,(J-1)
            ZP = AMAT(J,I,IZ)
            AMAT(J,I,IZ) = AMAT(I,J,IZ)
            AMAT(I,J,IZ) = ZP
          ENDDO
          AMAT(J,J,IZ)= AMAT(J,J,IZ)
        ENDDO
      ENDDO
C
      RETURN
      END

      SUBROUTINE QGETAM(AMAT,NDIM,NZ,LRA,LCA,IQSYM)
C*****************************************************************************
C
C     Quaternion GEneral matrix To Antihermitian matrix:
C     generate the ANTIHermitian component of the general matrix AMAT
C     return in AMAT.
C
C     Do A + A^dagger, where A is general matrix constructed from
C     antuhermitian operator - connection matrix T_x/y/z. 
C
C     Written by MI, september 2003, Prievidza, Slovakia.
C
C*****************************************************************************
#include "implicit.h"
      DIMENSION AMAT(LRA,LCA,NZ),IQSYM(4)
C
C     Generate A - A(dagger)
C

      DO IZ = 1, NZ
       IF (IQSYM(IZ).EQ.2) THEN ! A
         DO J = 1, NDIM
          DO I = 1, J
           AMAT(I,J,IZ) = AMAT(I,J,IZ)-AMAT(J,I,IZ)
           AMAT(J,I,IZ) = -AMAT(I,J,IZ)
          ENDDO
         ENDDO
       ELSEIF (IQSYM(IZ).EQ.1) THEN ! S
        DO J = 1, NDIM
         DO I = 1, J
           AMAT(I,J,IZ) = (AMAT(I,J,IZ)+AMAT(J,I,IZ))
           AMAT(J,I,IZ) = AMAT(I,J,IZ)
         ENDDO
        ENDDO
       ELSE
        CALL QUIT('QGETAM: IQSYM(IZ)<>1,2!')
       ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrtdac */
      SUBROUTINE WRTDAC (IT,N,X,IREC)
#include "implicit.h"
      DIMENSION X(N)
      WRITE (IT, REC=IREC) X
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck readac */
      SUBROUTINE READAC (IT,N,X,IREC)
#include "implicit.h"
#include "priunit.h"
      DIMENSION X(N)
      READ (IT, REC=IREC, IOSTAT=IOS, ERR=10) X
      RETURN
 10   CONTINUE
CMI   ...adding more on error report
      WRITE(LUPRI,'(/,6X,A)')  '*** Error in READAC ***'
      write(lupri,'(2X,A,I4)') 'file number IT=',IT
      write(lupri,'(2X,A,I8)') 'number of elements to be read N=',N
      write(lupri,'(2X,A,I4)') ' RECORD #: IREC=',IREC
      write(lupri,'(2X,A,I4)') ' IOSTAT=',IOS
      CALL QUIT('READAC: Error reading file')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prmutc */
      SUBROUTINE PRMUTC(AMAT,NCOL,NCTL,LUPRI)
C***********************************************************************
C
C  Based on OUTPAK by Hans Jorgen Aa. Jensen.
C
C PRMUTC prints a real symmetric matrix stored in column-packed upper
C triangular form (see diagram below) in formatted form with numbered
C rows and columns.  The input is as follows:
C
C        AMAT(*).............packed matrix
C        NCOL................number of columns to be output
C        NCTL................carriage control flag: 1 for single space,
C                                                   2 for double space,
C                                                   3 for triple space.
C
C The matrix elements are arranged in storage as follows:
C
C      1   2   4   7  11
C          3   5   8  12
C              6   9  13
C                 10  14
C                     15
C       and so on.
C
C PRMUTC is set up to handle 6 columns/page with a 6F20.14 format
C for the columns.  If a different number of columns is required, change
C formats 1000 and 2000, and initialize kcol with the new number of
C columns.
#include "implicit.h"
      DIMENSION AMAT(*)
      INTEGER BEGIN
      CHARACTER*1 ASA(3),BLANK,CTL,FX,FPRI
      CHARACTER   PFMT*27, COLUMN*8, LFMT*13,XFMT*6,RFMT*6
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
C     First check whether the matrix is all zero
C
      J = NCOL*(NCOL+1)/2
      AMAX = ZERO
      AMIN = ZERO
      DO 5 I=1,J
         AMAX = MAX(AMAX,ABS(AMAT(I)))
         AMIN = MIN(AMIN,ABS(AMAT(I)))
    5 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 200
      END IF
C
C     Determine output format
C
      IF (AMIN.GT.FFMIN.AND.AMAX.LE.FFMAX) THEN
C        use F output format
         LFMT = '(   A1,I7,2X,'
         XFMT = '(15X),'
         RFMT = 'F15.8)'
      ELSE
C        use 1PD output format
         LFMT = '(1P,A1,I7,2X,'
         XFMT = '(15X),'
         RFMT = 'E15.8)'
      END IF
C
C LAST is the last column number in the row currently being printed
C
      LAST = MIN(NCOL,KCOL)
C
C BEGIN is the first column number in the row currently being printed.
C
C.....BEGIN NON STANDARD DO LOOP.
      BEGIN= 1
 1050 CONTINUE
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 40 IROW = 1,LAST
            ICOL = MAX(IROW,BEGIN)
            IOFF = (ICOL*(ICOL-1))/2 + IROW
            NPRI = LAST - ICOL + 1
            DO 10 J = 0,(NPRI-1)
               IF (AMAT(IOFF+ICOL*J+(J*(J-1))/2).NE.ZERO) GO TO 20
   10       CONTINUE
            GO TO 40
   20       CONTINUE
            IDIG = ICHAR('0') + NPRI
            FPRI = CHAR(IDIG)
            NX   = ICOL - BEGIN
            IDIG = ICHAR('0') + NX
            FX   = CHAR(IDIG)
            IF(NX.GT.0) THEN
              PFMT = LFMT//FX//XFMT//FPRI//RFMT
            ELSE
              PFMT = LFMT//'      '//FPRI//RFMT
            ENDIF
            WRITE (LUPRI,PFMT) CTL,IROW,
     &           (AMAT(ICOL*J+(J*(J-1))/2+IOFF),J=0,(NPRI-1))
   40    CONTINUE
         LAST = MIN(LAST+KCOL,NCOL)
         BEGIN= BEGIN + KCOL
      IF (BEGIN.LE.NCOL) GO TO 1050
  200 CONTINUE
      RETURN
C
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
C2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtinfo */
      SUBROUTINE GTINFO(TEXT)
C
C Written by T.Saue Sept 1 1995
C Based on TSTAMP by Hans Joergen Aa. Jensen 9-Jul-1990
C
#include "implicit.h"
      CHARACTER TEXT*(*)
C
#if defined (SYS_CRAY) || defined (SYS_T3D)
      TEXT = ' '
      WRITE(TEXT(1:8),'(A8)') DATE()
      WRITE(TEXT(12:19),'(A8)') CLOCK()
#else
      CHARACTER*(24) FDATE
      TEXT = FDATE()
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prsymb */
      SUBROUTINE PRSYMB(IUNIT,SYMB,N,IX)
C***********************************************************************
C
C     Print a symbol/character N times on the same line
C
C     Written by T.Saue Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER SYMB*1,FMT*9
      FMT = '         '
      IF(IX.EQ.0) THEN
        WRITE(FMT(1:5),'(I3,A1)') N,'A1'
      ELSE
        WRITE(FMT,'(I3,A2,I3,A1)') IX,'X,',N,'A1'
      ENDIF
      WRITE(IUNIT,'('//FMT//')') (SYMB,I=1,N)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck readns */
      SUBROUTINE REDSYM(IRREP,IGROUP,KSOP,NELM)
C***********************************************************************
C
C     Determine point group of reduced symmetry for a non-symmetric
C     operator
C
C     Written by T.Saue - May 15 1996
C     Last revision: May 15 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "pgroup.h"
C
      DIMENSION KSOP(0:7),KBUF(0:7)
      DO I = 0,7
        KBUF(I) = 0
      ENDDO
      NELM = -1
C
C     Go through characters for the irrep to find the
C     totally symmetric operations
C     KSOP is a pointer to those....
C
      DO I = 0,MAXREP
      II = JSOP(I)
      IF(IXVAL(IRREP,II).EQ.1) THEN
        NELM = NELM + 1
        KSOP(NELM) = II
        KBUF(II)  = 1
      ENDIF
      ENDDO
      NELM = NELM+1
C
C     Sum up number of rotations, reflections and inversion to
C     determine the reduced group
C
      KROTS = KBUF(3)+KBUF(5)+KBUF(6)
      KREFL = KBUF(4)+KBUF(2)+KBUF(1)
      KINVC = KBUF(7)
      IGROUP = MIN(7,NINT((4*KROTS+8*KINVC+6*KREFL)/3.0))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gatmat */
      SUBROUTINE GATMAT(IOPT,LDA,N,A,B,IND,NIND)
C***********************************************************************
C
C     This routine will collect columns elements of square matrix A
C     in groups based on pointer index IND.
C     IOPT discerns two options:
C	IOPT = 0	pick out the elements of largest absolute value
C	IOPT = 1        sum up elements
C
C	Written by T.Saue, June 1994
C	Last revision: October 9 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION A(LDA,N),B(NIND,NIND),IND(N)
C
      CALL QENTER('GATMAT')
C
CTROND No initialization !!!!      CALL DZERO(B,NIND*NIND)
C
C     Pick out largest elements for each group
C     ========================================
C
      IF    (IOPT.EQ.0) THEN
        DO 10 J = 1,N
          JJ = IND(J)
          DO 20 I = 1,N
            II = IND(I)
            IF(ABS(A(I,J)).GT.ABS(B(II,JJ))) B(II,JJ) = ABS(A(I,J))
   20     CONTINUE
   10   CONTINUE
C
C     Sum up elements within each group
C     =================================
C
      ELSEIF(IOPT.EQ.1) THEN
        DO 30 J = 1,N
          JJ = IND(J)
          DO 40 I = 1,N
            II = IND(I)
            B(II,JJ) = B(II,JJ) + A(I,J)
   40     CONTINUE
   30   CONTINUE
      ELSE
        CALL QUIT('GATMAT:Unknown option IOPT!')
      ENDIF
C
      CALL QEXIT('GATMAT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gatvcc */
      SUBROUTINE GATVCC(IOPT,NVEC,A,N,B,IND,NIND)
C***********************************************************************
C
C     This routine will collect columns elements of vector A
C     in groups based on pointer index IND.
C     IOPT discerns two options:
C	IOPT = 0	pick out the elements of largest absolute value
C	IOPT = 1        sum up elements
C
C	Written by T.Saue, June 1994
C	Last revision: October 9 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION A(N,NVEC),B(NIND,NVEC),IND(N)
CTROND No initialization !!!!      CALL DZERO(B,NIND*NIND)
C
C     Pick out largest elements for each group
C     ========================================
C
      IF    (IOPT.EQ.0) THEN
        DO J = 1,NVEC
          DO I = 1,N
            II = IND(I)
            B(II,J) = MAX(B(II,J),ABS(A(I,J)))
          ENDDO
        ENDDO
C
C     Sum up elements within each group
C     =================================
C
      ELSEIF(IOPT.EQ.1) THEN
        DO J = 1,N
          DO I = 1,N
            II = IND(I)
            B(II,J) = B(II,J) + A(I,J)
          ENDDO
        ENDDO
      ELSE
        CALL QUIT('GATVCC:Unknown option IOPT!')
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qaxpy */
      SUBROUTINE QAXPY(N,FA,TA,A,IQPA,NZA,
     &                 X,NX,INCX,IQPX,NZX,Y,NY,INCY,IQPY,NZY)
C***********************************************************************
C     Computes a constant times a vector plus av vector
C
C        Y := A*X + Y
C
C     The vectors and constant may be real(NZ=1),complex(NZ=2),
C     or quaternionic(NZ = 4).
C     NX and NY are the total lengths of vectors X and Y,
C     whereas INCX and INCY are increments.
C
C     op(A) is determined by the corresponding character variable
C     FA:
C       FA = 'N'        Normal 
C       FA = 'C'        Complex conjugate
C     and TX
C       TA = 'N'        No transform    a+bj
C       TA = 'I'        i-transform   -i(a+bj)i = a - bj
C       TA = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TA = 'K'        k-transform   -k(a+bj)k = a*-b*j
C     Calls BLAS routine DAXPY
C
C     Written by T.Saue, University of Oslo, Apr 7 1997
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
C     Global variables
C
      DIMENSION A(NZA),X(NX,NZX),Y(NY,NZY),
     &          IQPA(NZA),IQPX(NZX),IQPY(NZY)
      CHARACTER FA*1,TA*1
C
C     Local variables
C
      DIMENSION IPA(4)
C
      IF (N.LE.0)       RETURN
C
      IF(N*INCX.GT.NX) CALL QUIT('QAXPY:N*INCX.GT.NX')
      IF(N*INCY.GT.NY) CALL QUIT('QAXPY:N*INCY.GT.NY')
C
C     Make pointer IPA
C     ================
C
      CALL IZERO(IPA,4)
      DO IZA = 1,NZA
        IPA(IQPA(IZA)) = IZA
      ENDDO
C
C     Constant A
C

C     Determine FX:
C     =============
C     Normal number...
      IF    (FA.EQ.'N') THEN
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        IFA = 2
      ELSE
        CALL QUIT('QAXPY:Unknown FA '//FA//' of scalar A!')
      ENDIF
C
C     Determine TX:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QAXPY:Unknown TA '//TA//' of scalar A!')
      ENDIF
C
      DO 10 IZY = 1,NZY
        IQY = IQPY(IZY)
        DO 20 IZX = 1,NZX
          IQX = IQPX(IZX)
          IQA = IQMULT(IQX,IQY,1)
          IZA = IPA(IQA)
          IF(IZA.EQ.0) GOTO 20
          IF(A(IZA).EQ.D0) GOTO 20
          AA = A(IZA)*(IQSIGN(IQA,IFA,ITA)*IQPHASE(IQA,IQX,1))
          CALL DAXPY(N,AA,X(1,IZX),INCX,Y(1,IZY),INCY)
   20   CONTINUE
   10 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qxapy */
      SUBROUTINE QXAPY(N,FA,TA,A,IQPA,NZA,
     &                         X,NX,INCX,IQPX,NZX,
     &                         Y,NY,INCY,IQPY,NZY)
C***********************************************************************
C     Computes a vector times a constant plus av vector
C
C        Y := X*A + Y
C
C     The vectors and constant may be real(NZ=1),complex(NZ=2),
C     or quaternionic(NZ = 4).
C     NX and NY are the total lengths of vectors X and Y,
C     whereas INCX and INCY are increments.
C
C     op(A) is determined by the corresponding character variable
C     FA:
C       FA = 'N'        Normal 
C       FA = 'C'        Complex conjugate
C     and TX
C       TA = 'N'        No transform    a+bj
C       TA = 'I'        i-transform   -i(a+bj)i = a - bj
C       TA = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TA = 'K'        k-transform   -k(a+bj)k = a*-b*j
C     Calls BLAS routine DAXPY
C
C     Written by T.Saue, University of Oslo, Apr 7 1997
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
C     Global variables
C
      DIMENSION A(NZA),X(NX,NZX),Y(NY,NZY),
     &          IQPA(NZA),IQPX(NZX),IQPY(NZY)
      CHARACTER FA*1,TA*1
C
C     Local variables
C
      DIMENSION IPA(4)
C
      IF (N.LE.0)       RETURN
C
      IF(N*INCX.GT.NX) CALL QUIT('QXAPY:N*INCX.GT.NX')
      IF(N*INCY.GT.NY) CALL QUIT('QXAPY:N*INCY.GT.NY')
C
C     Make pointer IPA
C     ================
C
      CALL IZERO(IPA,4)
      DO IZA = 1,NZA
        IPA(IQPA(IZA)) = IZA
      ENDDO
C
C     Constant A
C

C     Determine FX:
C     =============
C     Normal number...
      IF    (FA.EQ.'N') THEN
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        IFA = 2
      ELSE
        CALL QUIT('QXAPY:Unknown FA '//FA//' of scalar A!')
      ENDIF
C
C     Determine TX:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QXAPY:Unknown TA '//TA//' of scalar A!')
      ENDIF
C
      DO 10 IZY = 1,NZY
        IQY = IQPY(IZY)
        DO 20 IZX = 1,NZX
          IQX = IQPX(IZX)
          IQA = IQMULT(IQX,IQY,1)
          IZA = IPA(IQA)
          IF(IZA.EQ.0) GOTO 20
          IF(A(IZA).EQ.D0) GOTO 20
          AA = A(IZA)*(IQSIGN(IQA,IFA,ITA)*IQPHASE(IQX,IQA,1))
          CALL DAXPY(N,AA,X(1,IZX),INCX,Y(1,IZY),INCY)
   20   CONTINUE
   10 CONTINUE
      RETURN
C
      END
C  /* Deck fndamx */
      SUBROUTINE FNDAMX(VEC,NDIM,MXELMN,IPLACE,NELMN,IPRT,THRES)
C
C 890205 hjaaj. Based on FNDMN3 by Jeppe Olsen.
C Last revision 900405 hjaaj.
C
C FIND MXELMN/NELMN ELEMENTS IN VEC with highest absolute value.
C NELMN IS THE LARGEST NUMBER LOWER THAN MXELMN THAT DOES NOT
C SPLIT DEGENERATE PAIRS.
C ORIGINAL PLACES OF THE LOWEST ELEMENTS ARE STORED IN IPLACE.
C
#include "implicit.h"
      DIMENSION VEC(NDIM),IPLACE(*)
C
      PARAMETER ( D1 = 1.0D0 )
C
C. FIRST OCCURANCE OF LOWEST ELEMENT AND LARGEST ELEMENT
      XMIN = ABS(VEC(1))
      XMAX = ABS(VEC(1))
      IMIN = 1
      IMAX = 1
      DO 100 I = 2,NDIM
        AVECI = ABS(VEC(I))
        IF( AVECI .GT. XMAX) THEN
           XMAX = AVECI
           IMAX = I
        ELSE IF( AVECI .LT. XMIN) THEN
           XMIN = AVECI
           IMIN = I
        END IF
  100 CONTINUE
C
      IF (IPRT.GT.0) WRITE(6,*) ' >> output from FNDAMX <<'
      IF(IPRT .GE. 5 ) THEN
         WRITE(6,*) ' LOWEST VALUE AND PLACE ',XMIN,IMIN
         WRITE(6,*) ' HIGHST VALUE AND PLACE ',XMAX,IMAX
      END IF
C
      IPLACE(1) = IMAX
      NDEG = 1
      IF(IPRT .GE. 5 ) WRITE(6,'(A,I8,1P,D15.5,I8,I4)')
     &    'IELMNT value IMAX NDEG ', 1,VEC(IMAX),IMAX,NDEG
C
      ITOP = MIN(NDIM,MXELMN+1)
      DO 200 IELMNT = 2,ITOP
        XMAXPR = XMAX
        IMAXPR = IMAX
        XMAX   = -D1
        IMAX   = -1
        DO 150 I = 1,NDIM
          AVECI = ABS(VEC(I))
          IF (AVECI .GT. XMAX) THEN
            IF (AVECI .LT. XMAXPR) THEN
               IMAX = I
               XMAX = AVECI
            ELSE IF (AVECI .EQ. XMAXPR .AND. I .GT. IMAXPR) THEN
               IMAX = I
               XMAX = AVECI
               GO TO 151
            END IF
          END IF
  150   CONTINUE
  151   CONTINUE
        IF(XMAXPR-XMAX .LE. THRES )THEN
          NDEG = NDEG + 1
        ELSE
          NDEG = 1
        END IF
C
        IF( IELMNT .LE. MXELMN ) IPLACE(IELMNT) = IMAX
        IF( IPRT .GE. 5 ) WRITE(6,'(A,I8,1P,D15.5,I8,I4)')
     &    ' IELMNT value IMAX NDEG ',IELMNT,VEC(IMAX),IMAX,NDEG
  200 CONTINUE
C
C CHECK DEGENERACY ON LAST VALUE
      IF(MXELMN .LT. NDIM ) THEN
        NELMN = MXELMN+1-NDEG
      ELSE
        NELMN = NDIM
      END IF
      IF (IPRT .GT. 0)
     &   WRITE(6,*) ' NUMBER OF ELEMENTS OBTAINED IN FNDAMX ',NELMN
C
C
      IF( IPRT  .GE. 3 ) THEN
        WRITE(6,*) ' FROM FNDAMX : '
        WRITE(6,*) '   PLACES OF ELEMENTS of highest absolut value'
        CALL IWRTMA(IPLACE,1,NELMN ,1,NELMN )
      END IF
C
      IF( IPRT .GE. 1 ) THEN
       WRITE(6,*)
     & ' MAX AND MIN IN SELECTED SUPSPACE ',
     &   VEC(IPLACE(1)),VEC(IPLACE(NELMN))
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Extelm */
      SUBROUTINE EXTELM(FA,A,LDA,NA,B,LDB,NBR,NBC,IBR,IBC)
C***********************************************************************
C
C     Extract elements of matrix A(LDA,NA) and
C     insert into matrix B(LDB,NBC) using 
C     row array IBR and column array IBC.
C
C     Written by T.Saue Mar 24 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER FA*1
      DIMENSION A(LDA,NA),B(LDB,NBC),IBR(NBR),IBC(NBC)
C
C     Extract from normal matrix
C
      IF(FA.EQ.'N') THEN
        DO J = 1,NBC
          JA = IBC(J)
          DO I = 1,NBR
            B(I,J) = A(IBR(I),JA)
          ENDDO
        ENDDO
C
C     Extract from transposed matrix
C
      ELSEIF(FA.EQ.'T') THEN
        DO J = 1,NBC
          JA = IBC(J)
          DO I = 1,NBR
            B(I,J) = A(JA,IBR(I))
          ENDDO
        ENDDO
      ELSE
        CALL QUIT('EXTELM:Unknown FA '//FA//' of matrix A!')
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Extvec */
      SUBROUTINE EXTVEC(AVEC,LDA,BVEC,LDB,NB,NVEC,IBR)
C***********************************************************************
C
C     Extract elements of vectors A(LDA,NVEC) and
C     insert into vector BVEC(LDB,NVEC) using 
C     row array IBR.
C
C     Written by T.Saue Sep 2013
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AVEC(LDA,NVEC),BVEC(LDB,NVEC),IBR(NB)
C
      DO J = 1,NVEC
        DO I = 1,NB
          BVEC(I,J) = AVEC(IBR(I),J)
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lmnrep */
      SUBROUTINE LMNREP(IREP,NHKTA,KHKTA,LVALUE,MVALUE,NVALUE)
C***********************************************************************
C
C     Get out xyz-powers for a given irrep IREP.
C     Based on LMNVAL by tuh March 87
C
C     Written by T.Saue - Mar 24 1998
C***********************************************************************
#include "implicit.h"
#include "maxmom.h"
#include "xyzpow.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION LVALUE(KHKTA), MVALUE(KHKTA), NVALUE(KHKTA)
#include "symmet.h"
C
      ICOMP = 0
      DO 100 I = 1, KHKTA
      IF(ISYMAO(NHKTA,I).EQ.IREP) THEN
        ICOMP = ICOMP + 1
        LVALUE(ICOMP) = NHKTA - ISTEP(I)
        MVALUE(ICOMP) = MVAL(I)
        NVALUE(ICOMP) = NVAL(I)
      ENDIF
  100 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck qchold */
      SUBROUTINE QCHOLD(A,N,NZ,LRA,LCA,D,STOL,NEFF,JOB,IPIVOT)
C***********************************************************************
C
C     This routine performs upper triangle Cholesky decomposition of
C         * real symmetric matrices                     NZ = 1
C         * complex Hermitian matrices                  NZ = 2
C         * quaternionic Hermitian matrices             NZ = 4
C                                                (dagger)
C     Upper triangle Cholesky decomposition: A=UU
C     where U is an upper triangular matrix.  
C     On input/output A contains the lower triangle and diagonal of
C     the Hermitian matrix to be decomposed.
C     On output the upper triangle of the Cholesky matrix U is stored
C     in the upper triangle of A, and the diagonal in D.
C     With pivoting(JOB.ne.0) the Cholesky matrix has been scrambled
C     according to the ordering indicated by IPIVOT.
C     
C     Written by T.Saue March 25 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION A(LRA,LCA,NZ),D(N),IPIVOT(*),P(4),IP(4)

C
      CALL DZERO(P,4)
C
C     Initialization
C
      DO IZ = 1,NZ
        IP(IZ) = IZ
      ENDDO
      NA = LRA*LCA
      NEFF = 0
      J = 1
      KMAX = J-1
C
C     Decomposition without pivoting
C
      IF(JOB.EQ.0) THEN
        IF(A(J,J,1).LE.STOL) RETURN
        NEFF = J
        D(J) = SQRT(A(J,J,1))
        DO I = 2,N
          A(J,I,1) = A(I,J,1)/D(1)
          DO IZ = 2,NZ
            A(J,I,IZ) = -A(I,J,IZ)/D(1)
          ENDDO
        ENDDO
        DO J = 2,N
          KMAX = J-1
          CALL QDOT(D(J),IP,1,KMAX,
     &               'H',IP,A(1,J,1),NA,NZ,1,
     &               'N',IP,A(1,J,1),NA,NZ,1)
          D(J) = A(J,J,1) - D(J)
          IF(D(J).LE.STOL) RETURN
          NEFF = J
          D(J) = SQRT(D(J))
          DO I = (J+1),N
            CALL QDOT(P,IP,NZ,KMAX,
     &               'H',IP,A(1,I,1),NA,NZ,1,
     &               'N',IP,A(1,J,1),NA,NZ,1)
            A(J,I,1) = (A(I,J,1)-P(1))/D(J)
            DO IZ = 2,NZ
              A(J,I,IZ) = (P(IZ)-A(I,J,IZ))/D(J)
            ENDDO
          ENDDO
        ENDDO
      ELSE
C
C     Cholesky decomposition with full pivoting
C
C...intialize pivot array
        DO I = 1,N
          IPIVOT(I) = I
        ENDDO
C...scan diagonal  of matrix A for largest element
        D(J) = A(J,J,1)
        JJ   = J
        DO K = 2,N
          IF(A(K,K,1).GT.D(J)) THEN
            D(J) = A(K,K,1)
            JJ   = K
          ENDIF
        ENDDO
        IF(D(J).LE.STOL) RETURN
        IF(JJ.NE.J) THEN
          NEFF = IPIVOT(JJ)
          IPIVOT(JJ) = IPIVOT(J)
          IPIVOT(J)  = NEFF
          JJ = NEFF
        ENDIF
        NEFF = J
C...generate elements of Cholesky matrix
        D(J) = SQRT(D(J))
        DO I = 2,N
          II = IPIVOT(I)
          IF(II.GE.JJ) THEN
            A(J,I,1) = A(II,JJ,1)/D(J)
            DO IZ = 2,NZ
              A(J,I,IZ) = -A(II,JJ,IZ)/D(J)
            ENDDO
          ELSE
            A(J,I,1) = A(JJ,II,1)/D(J)
            DO IZ = 2,NZ
              A(J,I,IZ) = A(JJ,II,IZ)/D(J)
            ENDDO
          ENDIF
        ENDDO
C...iterate
        DO J = 2,N
          KMAX = J-1
C...      generate diagonal elements of Cholesky matrix
C...      and find largest element
          JJ = IPIVOT(J)
          CALL QDOT(D(J),IP,1,KMAX,'H',IP,A(1,J,1),NA,NZ,1,
     &                             'N',IP,A(1,J,1),NA,NZ,1)
          D(J) = A(JJ,JJ,1) - D(J)
          DO K = (J+1),N
            KK = IPIVOT(K)
            CALL QDOT(D(K),IP,1,KMAX,'H',IP,A(1,K,1),NA,NZ,1,
     &                               'N',IP,A(1,K,1),NA,NZ,1)
            D(K) = A(KK,KK,1) - D(K)
            IF(D(K).GT.D(J)) THEN
              D(J) = D(K)
              JJ   = K
            ENDIF
          ENDDO
          IF(D(J).LE.STOL) RETURN
C..       if pivoting, interchange columns of previous rows of Cholesky matrix
          IF(JJ.NE.IPIVOT(J)) THEN
            DO IZ = 1,NZ
              CALL DSWAP(KMAX,A(1,J,IZ),1,A(1,JJ,IZ),1)
            ENDDO
            NEFF = IPIVOT(JJ)
            IPIVOT(JJ) = IPIVOT(J)
            IPIVOT(J)  = NEFF
            JJ = NEFF
          ENDIF
          NEFF = J
C...      generate elements of Cholesky matrix
          D(J) = SQRT(D(J))
          DO I = (J+1),N
            II = IPIVOT(I)
            CALL QDOT(P,IP,NZ,KMAX,
     &               'H',IP,A(1,I,1),NA,NZ,1,
     &               'N',IP,A(1,J,1),NA,NZ,1)
            IF(II.GE.JJ) THEN
              A(J,I,1) = (A(II,JJ,1)-P(1))/D(J)
              DO IZ = 2,NZ
                A(J,I,IZ) = (P(IZ)-A(II,JJ,IZ))/D(J)
              ENDDO
            ELSE
              A(J,I,1) = (A(JJ,II,1)-P(1))/D(J)
              DO IZ = 2,NZ
                A(J,I,IZ) = (P(IZ)+A(JJ,II,IZ))/D(J)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qdot */
      SUBROUTINE QDOT(P,IQPP,NZP,N,FX,IQPX,X,LRX,NZX,INCX,
     &                             FY,IQPY,Y,LRY,NZY,INCY)
C***********************************************************************
C
C     Performs generally real,complex or quaternion dot product
C     (determined by NZ)
C     FX/FY = 'N' : Normal vectors
C     FX/FY = 'H' : Hermitian conjugate
C
C     Written by T.Saue March 25 1999
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION X(LRX,NZX),Y(LRY,NZY),P(NZP),
     &          IQPX(NZX),IQPY(NZY),IQPP(NZP)
      CHARACTER  FX*1,FY*1
      DIMENSION IPY(4)
C
C
C     Initialization
C
      IFX = 1
      IF(FX.EQ.'H') IFX = 2
      IFY = 1
      IF(FY.EQ.'H') IFY = 2
C
C     Make pointer IPY
C     ================
C
      CALL IZERO(IPY,4)
      DO IZY = 1,NZY
        IPY(IQPY(IZY)) = IZY
      ENDDO
C
      DO IZP = 1,NZP
        IQP = IQPP(IZP)
        P(IZP) = D0
        DO 10 IZX = 1,NZX
          IQX = IQPX(IZX)
          IQY = IQMULT(IQX,IQP,1)
          IZY = IPY(IQY)
          IF(IZY.EQ.0) GOTO 10
          FAC = IQSIGN(IQX,IFX,1)*IQSIGN(IQY,IFY,1)*IQPHASE(IQX,IQY,1)
          P(IZP) = P(IZP) + FAC*DDOT(N,X(1,IZX),INCX,Y(1,IZY),INCY)
 10     CONTINUE
      ENDDO
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qdotm */
      SUBROUTINE QDOTM(P,IQPP,NZP,FX,IQPX,X,NX,LRX,LCX,NZX,INCX,
     &                            FY,IQPY,Y,NY,LRY,LCY,NZY,INCY,S,LDS)
C***********************************************************************
C
C     Performs generally real,complex or quaternion dot product
C     (determined by NZ)
C     FX/FY = 'N' : Normal vectors
C     FX/FY = 'H' : Hermitian conjugate
C     with REAL metric S
C
C     Written by T.Saue Nov 3 2000
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION X(LRX,LCX,NZX),Y(LRY,LCY,NZY),S(LDS,NY),P(NZP),
     &          IQPP(NZP),IQPX(NZX),IQPY(NZY)
      CHARACTER  FX*1,FY*1
      DIMENSION IPY(4)
C
C
C     Found out INC[XY] is not used! /hjaaj 4-Nov-2001
C
      IF (INCX .NE. 1 .OR. INCY .NE. 1) THEN
         CALL QUIT('QDOTM not implemented with INCX/Y .ne. 1')
      END IF
C
C     Initialization
C
      IFX = 1
      IF(FX.EQ.'H') IFX = 2
      IFY = 1
      IF(FY.EQ.'H') IFY = 2
C
C
C     Make pointer IPY
C     ================
C
      CALL IZERO(IPY,4)
      DO IZY = 1,NZY
        IPY(IQPY(IZY)) = IZY
      ENDDO
      DO IZP = 1,NZP
        IQP = IQPP(IZP)
        P(IZP) = D0
        DO 10 IZX = 1,NZX
          IQX = IQPX(IZX)
          IQY = IQMULT(IQX,IQP,1)
          IZY = IPY(IQY)
          IF(IZY.EQ.0) GOTO 10
          FAC = IQSIGN(IQX,IFX,1)*IQSIGN(IQY,IFY,1)*IQPHASE(IQX,IQY,1)
          P(IZP) = P(IZP)+FAC*DDOTM(NX,NY,LDS,X(1,1,IZX),S,Y(1,1,IZY))
 10     CONTINUE
      ENDDO
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qnrm2 */
      FUNCTION QNRM2(N,NZ,X,NX,INCX)
C***********************************************************************
C
C     Gives squared norm of real,complex or quaternion vector
C     (determined by NZ)
C
C     Written by T.Saue March 25 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION X(NX,NZ)
C
      QNRM2 = D0
      DO IZ = 1,NZ
        QNRM2 = QNRM2 + DDOT(N,X(1,IZ),INCX,X(1,IZ),INCX)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck qchols */
      SUBROUTINE QCHOLS(A,N,NEFF,NSIM,NZ,LRA,LCA,D,
     &                  B,LRB,LCB,X,LRX,LCX,JOB,IPIVOT,W)
C***********************************************************************
C
C     This routine solves the NSIM simultaneous equations
C         A(n,n,nz)*x(n,nsim,nz) = b(n,nsim,nz)
C     where A and D is output from the Cholesky decomposition
C     routine QCHOLD. The upper triagonal of A contains the 
C     upper triangular Cholesky matrix and the diagonal
C     and lower triangle is the original matrix A
C         * a real symmetric matrix                     NZ = 1
C         * a complex Hermitian matrix                  NZ = 2
C         * a quaternion Hermitian matrix               NZ = 4
C     The vector D contains the diagonal elements of the
C     Cholesky matrix.
C
C     The system of linear equations is written as
C
C        Ax = U(dagger)Ux = U(dagger)y = b ; Ux = y
C
C     JOB.NE.0 means that full pivoting should be used....
C     W(N) is a work array needed for sorting in case of pivot.
C     
C     Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0, D0=0.0D0)
      DIMENSION A(LRA,LCA,NZ),D(N),B(LRB,LCB,NZ),X(LRX,LCX,NZ),
     &          IPIVOT(*),W(*)
      DIMENSION P(4),IQA(4),IQB(4)
      CALL DZERO(P,4)
C
C     Initialization
C
      DO IZ = 1,NZ
        IQA(IZ) = IZ
        IQB(IZ) = IZ
      ENDDO
      NA = LRA*LCA
c     NB = LRB*LCB
      NX = LRX*LCX
C
C     First solve U(dagger)y = b by forward substitution
C
      DO J = 1,NSIM
C
C       In case of pivot, reorder coefficient array b
C
        IF(JOB.NE.0) THEN
          DO IZ = 1,NZ
            CALL DCOPY(N,B(1,J,IZ),1,W,1)
            DO I = 1,NEFF
              B(I,J,IZ) = W(IPIVOT(I))
            ENDDO
          ENDDO
        ENDIF
        DO IZ = 1,NZ
          X(1,J,IZ) = B(1,J,IZ)/D(1)
        ENDDO
        DO I = 2,NEFF
          KMAX = I-1
            CALL QDOT(P,IQA,NZ,KMAX,
     &               'H',IQA,A(1,I,1),NA,NZ,1,
     &               'N',IQA,X(1,J,1),NX,NZ,1)
          DO IZ = 1,NZ
            X(I,J,IZ) = (B(I,J,IZ)-P(IZ))/D(I)
          ENDDO
        ENDDO
      ENDDO
C
C     Solve Ux = y by backsubstitution
C
      DO J = 1,NSIM
        DO IZ = 1,NZ
          X(NEFF,J,IZ) = X(NEFF,J,IZ)/D(NEFF)
        ENDDO
      ENDDO
      DO I = (NEFF-1),1,-1
        P(2) = D1/D(I)
        P(1) = -P(2)
        NELM = NEFF-I
        KINT = I+1
        CALL QGEMM(1,NSIM,NELM,P(1),
     &             'N','N',IQA,A(I,KINT,1),LRA,LCA,NZ,
     &             'N','N',IQA,X(KINT,1,1),LRX,LCX,NZ,
     &                P(2),IQB,X(I,1,1),LRX,LCX,NZ)
      ENDDO
C
C     In case of pivot, reorder solution vectors
C
      IF(JOB.NE.0) THEN
        DO J = 1,NSIM
          DO IZ = 1,NZ
            CALL DCOPY(NEFF,X(1,J,IZ),1,W,1)
            DO I = 1,NEFF
              X(IPIVOT(I),J,IZ) = W(I)
            ENDDO
            DO I = (NEFF+1),N
              X(IPIVOT(I),J,IZ) = D0
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* qhmrst */
      SUBROUTINE QHMRST(A,N,NZ,LRA,LCA)
C***********************************************************************
C
C     In input is given a (real/complex/quaternion) Hermitian
C     matrix where the upper diagonal has been destroyed. Restore.
C
C     Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(LRA,LCA,NZ)
C
      DO J = 1,N
        DO I = 1,(J-1)
          A(I,J,1) = A(J,I,1)
        ENDDO
      ENDDO
      DO IZ = 2,NZ
        DO J = 1,N
          DO I = 1,(J-1)
            A(I,J,IZ) = -A(J,I,IZ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QSWAP(N,DX,ISYMX,FAC,DAGGER)
C***********************************************************************
C
C     Swap upper and lower part of quaternionic vector and scale
C     with a factor FAC.
C     If DAGGER then also complex conjugate vector.
C
C     Calls BLAS routine DSWAP.
C
C     Written by panor 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
C
      PARAMETER ( DM1 = -1.0D0 )
C
      LOGICAL DAGGER
      DIMENSION DX(N,NZ)
C
      IF (N.LE.0) RETURN
C
      DO IZ=1,NZ
         CALL DSWAP(N/2,DX(1,IZ),1,DX(N/2+1,IZ),1)
         IF (DAGGER .AND. IPQTOQ(IZ,ISYMX-1).GT.1) 
     &        CALL DSCAL(N,DM1,DX(1,IZ),1)
      END DO
      CALL DSCAL(N*NZ,FAC,DX,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck curve */
      SUBROUTINE CURVE(NORDER,NPOINTS,A,B,X,Y,CHISQ,KVPT)
C***********************************************************************
C
C     Least square polynomial fit solving normal equations
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
      PARAMETER (D0=0.0D0)
      DIMENSION A(*),B(*),KVPT(*)
      DIMENSION X(NPOINTS),Y(NPOINTS)
      NDIM = NORDER + 1
C* Stabilize by subtracting mean from Y-vector:
      YMEAN = D0
      DO I = 1,NPOINTS
        YMEAN = YMEAN + Y(I)
      ENDDO
      YMEAN = YMEAN/NPOINTS
      DO I = 1,NPOINTS
        Y(I) = Y(I) - YMEAN
      ENDDO
C* Calculate upper triangle column packed matrix A :
      IOFF = 0
      DO 10 J = 1,NDIM
        DO 20 I = 1,J
          IOFF = IOFF + 1
          IORDER = I + J - 2
          A(IOFF) = XN(IORDER,NPOINTS,X)
   20 CONTINUE
   10 CONTINUE  
C* Calculate column vector b:
      DO 30 I = 1,NDIM
        IORDER = I - 1
        B(I) = XNY(IORDER,NPOINTS,X,Y)
   30 CONTINUE  
C* Solve linear equations
      CALL DSPSOL(NDIM,1,A,B,KVPT,INFO)      
      DO I = 1,NPOINTS
        Y(I) = Y(I) + YMEAN
      ENDDO
      B(1) = B(1) + YMEAN
C*    Calculate chisquare
      CHISQ = D0
      DO I = 1,NPOINTS
        DIFF  = Y(I) - POLVAL(NORDER,B,X(I))
        CHISQ = CHISQ + DIFF*DIFF
      ENDDO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Polsvd */
      SUBROUTINE POLSVD(ND,NP,A,B,X,Y,W,RV,CHISQ,ISKIP)
C***********************************************************************
C
C     Least square polynomial fit using singular value decomposition
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER(TOL=1.0D-16)
      DIMENSION A(NP,ND),B(NP),X(NP),Y(NP),W(ND),RV(ND)
C
      NORDER = ND-1
C* Stabilize by subtracting mean from Y-vector:
      YMEAN = D0
      DO I = 1,NP
        YMEAN = YMEAN + Y(I)
      ENDDO
      YMEAN = YMEAN/NP
      DO I = 1,NP
        Y(I) = Y(I) - YMEAN
      ENDDO
C     Set up design matrix A 
      CALL DCOPY(NP,D1,0,A(1,1),1)
      CALL DCOPY(NP,X,1,A(1,2),1)
      DO J = 3,ND
        DO I = 1,NP
          A(I,J)=A(I,J-1)*A(I,2)
        ENDDO
      ENDDO
C     Set up response vector B
      CALL DCOPY(NP,Y,1,B,1)
C     Perform singular value decomposition
      CALL MINFIT(NP,NP,ND,A,W,1,B,IERR,RV)
C     Edit the singular values: Note that this is done to B !!
      WA = D0
      DO I = 1,ND
        IF(W(I).GT.WA) WA=W(I)
      ENDDO
      WA = TOL*WA
      ISKIP = 0
      DO I = 1,ND
        IF(W(I).LT.WA) THEN
          B(I)=D0
          ISKIP = ISKIP + 1
        ELSE
          B(I)=B(I)/W(I)
        ENDIF
      ENDDO
C     Backsubsitution to get parameters      
      CALL DGEMV('N',ND,ND,D1,A,NP,B,1,D0,W,1)
      DO I = 1,NP
        Y(I) = Y(I) + YMEAN
      ENDDO
      W(1) = W(1) + YMEAN
C     Calculate chisquare
      CHISQ = D0
      DO I = 1,NP
        DIFF  = Y(I) - POLVAL(NORDER,W,X(I))
        CHISQ = CHISQ + DIFF*DIFF
      ENDDO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE NEWRAP(NORDER,DTOL,B,XVAL,NITER,IERR)
#include "implicit.h"
#include "priunit.h"
      DIMENSION B(NORDER + 1)
      IERR = 0
      NSTEP = 0
   10 CONTINUE
      NSTEP = NSTEP + 1
      XOLD = XVAL  
      XFAC = POLDER(NORDER,B,XOLD)/POL2DER(NORDER,B,XOLD)
      XVAL = XOLD - XFAC
      IF((ABS((XVAL-XOLD)).GT.DTOL).AND.(NSTEP.LE.NITER)) GOTO 10
      IF(NSTEP.GT.NITER) THEN
        IERR = 1
      ENDIF
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION XN(NORDER,NPOINTS,X)
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
      DIMENSION X(NPOINTS)
      XN = D0
      DO 10 I = 1,NPOINTS
        XN = XN + X(I)**NORDER
   10 CONTINUE    
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION XNY(NORDER,NPOINTS,X,Y)
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
      DIMENSION X(NPOINTS),Y(NPOINTS)
      XNY = D0
      FAC = NORDER
      DO 10 I = 1,NPOINTS
        XNYBUF = Y(I)
        XBUF = X(I)
        DO 20 J = 1,NORDER
          XNYBUF = XNYBUF*XBUF
   20   CONTINUE    
C      XNY = XNY + XNYBUF
      XNY = XNY + Y(I)*(X(I)**FAC)
   10 CONTINUE    
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POLVAL(NORDER,B,XVAL)
#include "implicit.h"
      DIMENSION B(NORDER+1)
      POLVAL = B(1)
      XBUF = 1
      DO 10 I = 2,(NORDER+1)
        XBUF = XBUF*XVAL
        POLVAL = POLVAL + XBUF*B(I)
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POLDER(NORDER,B,XVAL)
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
      DIMENSION B(NORDER+1)
      IF(NORDER.LT.1) THEN
        POLDER = D0
        RETURN
      ENDIF
      POLDER = B(2)
      DO 10 I = 2,NORDER
        POLDER = POLDER + (XVAL**(I-1))*I*B(I+1)
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POL2DER(NORDER,B,XVAL)
#include "implicit.h"
      PARAMETER(D0 = 0.0D0,D2 = 2.0D0)
      DIMENSION B(NORDER+1)
      IF(NORDER.LT.2) THEN
        POL2DER = D0
        RETURN
      ENDIF
      POL2DER = D2*B(3)
      DO 10 I = 3,NORDER
        POL2DER = POL2DER + (XVAL**(I-2))*I*(I-1)*B(I+1)
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POL3DER(NORDER,B,XVAL)
#include "implicit.h"
      PARAMETER(D0 = 0.0D0,D6 = 6.0D0)
      DIMENSION B(NORDER+1)
      IF(NORDER.LT.3) THEN
        POL3DER = D0
        RETURN
      ENDIF
      POL3DER = D6*B(4)
      DO 10 I = 4,NORDER
        POL3DER = POL3DER + (XVAL**(I-3))*(I-2)*(I-1)*I*B(I+1)
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POLNDER(NORDER,B,XVAL,N)
#include "implicit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION B(NORDER+1)
      IF(NORDER.LT.N) THEN
        POLNDER = D0
        RETURN
      ENDIF
      FAC = D1
      DO J = 2,N
        FAC = FAC*J
      ENDDO
      POLNDER = FAC*B(N+1)
      DO 10 I = N+1,NORDER
        FAC = I
        DO J = 1,N-1
          FAC = FAC*(I-J)
        ENDDO
        POLNDER = POLNDER + (XVAL**(I-N))*B(I+1)*FAC
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION MIND(NPOINTS,Y)
#include "implicit.h"
      DIMENSION Y(NPOINTS)
      YVAL = Y(1)
      DO 10 I = 2,NPOINTS
        IF(Y(I).LE.YVAL) THEN
          YVAL = Y(I)
          MIND = I
        ENDIF
   10 CONTINUE
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qtrans */
      SUBROUTINE QTRANS(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,
     &                  FAO,LRAO,LCAO,NZAO,IQAO,
     &                  FMO,LRMO,LCMO,NZMO,IQMO,
     &                  TM1,LR1,LC1,NZTM1,IQTM1,
     &                  TM2,LR2,LC2,NZTM2,IQTM2,
     &                  BUF,LBUF,IPRINT)
C***********************************************************************
C
C     This routine performs unitary transformations indicated by TYP:
C
C       TYP = AOMO : AO-to-MO-transformation  FMO = (C+)FAOC + FADD*FMO
C       TYP = MOAO : MO-to-AO-transformation  FAO = CFMO(C+) + FADD*FAO
C
C     FADD gives the possibility of adding the results
C
C     TREV indicates symmetry of F-matrix under timereversal
C
C       TREV = 'S' - symmetric     (t =  1)
C       TREV = 'A' - anti-symmtric (t = -1)
C
C     The unitary C-matrix has a time-symmetric structure.
C
C     The AO-to-MO-transformation can be set up quaternionically as
C
C       FMO = [(Ca+)-t(CbT)j][(FAOa)+(FAOb)j][(Ca)+(Cb)j]
C
C       F(NBAS,NBAS,NZ) --->  F(NORB,NORB,NZ)
C
C     and the MO-to-Ao-transformation can be set up as:
C
C       FAO = [(Ca)+t(Cb)j][(FMOa)+(FMOb)j][(Ca+)-(Cb)j]
C
C       F(NORB,NORB,NZ) --->  F(NBAS,NBAS,NZ)
C
C     Written by T.Saue,January 1995
C     Last revision: 29 January 1997 : Luuk Visscher
C**********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
      CHARACTER TREV*1,TRA*1,TYP*4
      DIMENSION FAO(LRAO,LCAO,*),IQAO(*),
     &          FMO(LRMO,LCMO,*),IQMO(*),
     &          TM1(LR1,LC1,*),IQTM1(*),
     &          TM2(LR2,LC2,*),IQTM2(*),
     &          BUF(LBUF)
      DIMENSION IQW(4)
C
      CALL QENTER('QTRANS')
C
C     Add the results on to the contents of the MO matrix ?
C
!#define DEBUG
#ifdef DEBUG
      print*, 'QTRANS: The AO matrix'
      CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IQAO,LUPRI)
      print*, 'QTRANS: The transformation matrix TM1'
      CALL PRQMAT(TM1,NRAO,NCMO,LR1,LC1,NZTM1,IQTM1,LUPRI)
      print*, 'QTRANS: The transformation matrix TM2'
      CALL PRQMAT(TM2,NRAO,NCMO,LR2,LC2,NZTM2,IQTM2,LUPRI)
#undef DEBUG
#endif
C
C     Time reversal symmetry of F-matrix
C
      IF    (TREV.EQ.'S') THEN
        TRA = 'N'
      ELSEIF(TREV.EQ.'A') THEN
        TRA = 'I'
      ELSE
        CALL QUIT('QTRANS: Unknown keyword TREV '//TREV)
      ENDIF
C
C     *******************************************************
C     *****  AO - to - MO -transformation : F' = (C+)FC *****
C     *******************************************************
C
      IF(TYP.EQ.'AOMO') THEN
C
C       Look for the most efficient transformation
C       ------------------------------------------
C
        IF (NCAO.GE.NRAO) THEN
C
C         Perform the two steps:
C           1. W  = FC
C           2. F' = (C+)W
C
          CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
          NBUF = NRAO*NCMO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
             CALL QUIT('ERROR: need more work space')
          ENDIF
         
          IF (IPRINT.GE.12) THEN
            CALL HEADER('QTRANS: AO-to-MO-transformation',-1)

            WRITE(LUPRI,'(2X,A,4I2)')
     &      'QTRANS: The AO matrix; IQAO=',(IQAO(IZ),IZ=1,NZAO)
            WRITE(LUPRI,'(2X,A,A1)') ' Time rev=',TREV
            CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IQAO,LUPRI)

            WRITE(LUPRI,'(/2X,A,4I2)')
     &      'QTRANS: The transformation matrix TM1; IQTM1=',
     &      (IQTM1(IZ),IZ=1,NZTM1)
            CALL PRQMAT(TM1,NRAO,NCMO,LR1,LC1,NZTM1,IQTM1,LUPRI)

            WRITE(LUPRI,'(/2X,A,4I2)')
     &      'QTRANS: The transformation matrix TM2; IQTM2=',
     &      (IQTM2(IZ),IZ=1,NZTM2)
            CALL PRQMAT(TM2,NRAO,NCMO,LR2,LC2,NZTM2,IQTM2,LUPRI)
          ENDIF
C
C         First part of transformation: W = Fao C
C         ------------------------------------
C
          CALL QGEMM(NRAO,NCMO,NCAO,D1,
     &              'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &              'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &               D0,IQW,BUF,NRAO,NCMO,NZW)

          IF (IPRINT.GE.13) THEN
            WRITE(LUPRI,'(/2X,A,4I2)')
     &      'QTRANS: First part transformation matrix W, W=FC;'//
     &      ' IQW=',(IQW(IZ),IZ=1,NZW)
            CALL PRQMAT(BUF,NRAO,NCMO,NRAO,NCMO,NZW,IQW,LUPRI)
          ENDIF
C
C         Second part of transformation: Fmo = (C+)W
C         ----------------------------------------
C
          CALL QGEMM(NRMO,NCMO,NRAO,D1,
     &               'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,
     &               'N','N',IQW,BUF,NRAO,NCMO,NZW,
     &               FADD,IQMO,FMO,LRMO,LCMO,NZMO)

          IF (IPRINT.GE.13) THEN
            WRITE(LUPRI,'(/2X,A,4I2)')
     &     'QTRANS: Second part of transformation H,H=(C+).W; '//
     &     'IQMO=',(IQMO(IZ),IZ=1,NZMO)
            CALL PRQMAT(FMO,NRMO,NCMO,LRMO,LCMO,NZMO,IQMO,LUPRI)
          ENDIF

        ELSE
C
C         Perform the two steps:
C           1. W  = (C+)F
C           2. F' = WC
C
          CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
          NBUF = NRMO*NCAO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
             CALL QUIT('ERROR: need more work space')
          ENDIF

C
C         First part of transformation: W = (C+)F
C         ----------------------------------------
C
          CALL QGEMM(NRMO,NCAO,NRAO,D1,
     &       'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,
     &       'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &        D0,IQW,BUF,NRMO,NCAO,NZW)

	  IF(IPRINT.GE.15) THEN
        CALL HEADER('QTRANS:first part of transformation,W = (C+)F ',-1)
        CALL PRQMAT(BUF,NRMO,NCAO,NRMO,NCAO,NZW,IQW,LUPRI)       
	  ENDIF
	 
C 
C         Second part of transformation: H = WC
C         -------------------------------------
C
          CALL QGEMM(NRMO,NCMO,NCAO,D1,
     &       'N','N',IQW ,BUF,NRMO,NCAO,NZW,
     &       'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
C
        ENDIF
		
        IF(IPRINT.GE.11) THEN
          CALL HEADER('QTRANS: Final AO-to-MO-transformed matrix',-1)
          CALL PRQMAT(FMO,NRMO,NCMO,LRMO,LCMO,NZMO,IQMO,LUPRI)
        ENDIF
C
C     *******************************************************
C     *****  MO - to - AO -transformation : F' = CF(C+) *****
C     *******************************************************
C
      ELSEIF(TYP.EQ.'MOAO') THEN
C
#ifdef DEBUG
      print*, 'QTRANS: The MO matrix'
      CALL PRQMAT(FMO,NRMO,NCMO,LRMO,LCMO,NZMO,IQMO,LUPRI)
      print*, 'QTRANS: The transformation matrix TM1'
      CALL PRQMAT(TM1,NRAO,NCMO,LR1,LC1,NZTM1,IQTM1,LUPRI)
      print*, 'QTRANS: The transformation matrix TM2'
      CALL PRQMAT(TM2,NRAO,NCMO,LR2,LC2,NZTM2,IQTM2,LUPRI)
#endif
C       Look for the most efficient transformation
C       ------------------------------------------
C
        IF(NCMO.GE.NRMO) THEN
          CALL IQPACK(IQMO,NZMO,IQTM2,NZTM2,IQW,NZW)
          NBUF = NRMO*NCAO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
             CALL QUIT('ERROR: need more work space')
          ENDIF
C
C         First part of transformation: W = Fmo(C+)
C
          CALL QGEMM(NRMO,NCAO,NCMO,D1,
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       D0,IQW,BUF,NRMO,NCAO,NZW)
C
C         Second part of transformation: Fao = CW
C
          CALL QGEMM(NRAO,NCAO,NRMO,D1,
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQW ,BUF,NRMO,NCAO,NZW,
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
        ELSE
          CALL IQPACK(IQMO,NZMO,IQTM1,NZTM1,IQW,NZW)
          NBUF = NRAO*NCMO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
             CALL QUIT('ERROR: need more work space')
          ENDIF
C
C         First part of transformation: W = CFmo
C
          CALL QGEMM(NRAO,NCMO,NRMO,D1,
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,
     &       D0,IQW,BUF,NRAO,NCMO,NZW)
C
C         Second part of transformation: Fao = W(C+)
C
          CALL QGEMM(NRAO,NCAO,NCMO,D1,
     &       'N','N',IQW ,BUF,NRAO,NCMO,NZW,
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
        ENDIF
       IF(IPRINT.GE.11) THEN
          CALL HEADER('QTRANS: MO-to-AO-transformed Fock matrix',-1)
          CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IQAO,LUPRI)
        ENDIF
      ELSE
        CALL QUIT('FTRANS: Unknown keyword TYP '//TYP)
      ENDIF
C
      CALL QEXIT('QTRANS')
      RETURN
C
      END

      SUBROUTINE QTRANS90(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,
     &                  FAO,LRAO,LCAO,NZAO,IQAO,
     &                  FMO,LRMO,LCMO,NZMO,IQMO,
     &                  TM1,LR1,LC1,NZTM1,IQTM1,
     &                  TM2,LR2,LC2,NZTM2,IQTM2,IPRINT)
c     Like QTRANS, but uses fortran 90 memory allocation internally
      use memory_allocator
#include "implicit.h"
      CHARACTER TREV*1,TYP*4
      real(8) :: FAO(LRAO,LCAO,NZAO)
      real(8) :: FMO(LRMO,LCMO,NZMO)
      real(8) :: TM1(LR1 ,LC1 ,NZTM1)
      real(8) :: TM2(LR2 ,LC2 ,NZTM2)
      integer :: IQAO(nzao),IQMO(nzmo),IQTM1(nztm1),IQTM2(nztm2)
      real(8), allocatable :: work(:)
      integer IQW(4)
C     find out how much buffer is needed in qtrans
      IF(TYP.EQ.'AOMO') THEN
         IF (NCAO.GE.NRAO) THEN
            CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
            NBUF = NRAO*NCMO*NZW
         else
            CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
            NBUF = NRMO*NCAO*NZW
         endif
      else
         IF(NCMO.GE.NRMO) THEN
            CALL IQPACK(IQMO,NZMO,IQTM2,NZTM2,IQW,NZW)
            NBUF = NRMO*NCAO*NZW
         else
            CALL IQPACK(IQMO,NZMO,IQTM1,NZTM1,IQW,NZW)
            NBUF = NRAO*NCMO*NZW
         endif
      endif

      if(nbuf.le.0)then
        call quit(' *** error in QTRANS90: buffer space allocation not
     & possible for length <= 0.***')
      end if

      call alloc(work,nbuf, id='QTRANS90-buffer-space')
!     allocate( work(nbuf), stat=ierr )
!     if (ierr.ne.0) call quit('QTRANS90: error in allocation !')
      work = 0
      call QTRANS(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,
     &                  FAO,LRAO,LCAO,NZAO,IQAO,
     &                  FMO,LRMO,LCMO,NZMO,IQMO,
     &                  TM1,LR1,LC1,NZTM1,IQTM1,
     &                  TM2,LR2,LC2,NZTM2,IQTM2,
     &                  work,nbuf,IPRINT)
      call dealloc(work, id='QTRANS90-buffer-space')
CMI    ... necessary to eliminate gfortran memory bug...
!     deallocate(work,  stat=ierr)
!     if (ierr.ne.0) call quit('QTRANS90: error in dallocation !')
      
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Sqstra */
      SUBROUTINE SQSTRA(IOPT,A,N,LRA,LCA,NZA,WORK,LWORK)
C***********************************************************************
C
C     IOPT = -1:
C     Generate the matrix  Inv{Sqrt(S)} to be used in 
C     symmetric orthonormalization
C   
C     IOPT = +1:
C     Generate the matrix Sqrt(S)
C
C     IOPT = +2:
CMI   Generate the matrix Inv(S)     
C
C     NOTE that matrix A gets destroyed in the process !
C
C     Written by T. Saue Oct 10 2004
C     Last changes: M.Ilias, Sept 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(LRA,LCA,NZA),WORK(LWORK)
C
C     Memory allocation
C
#include "memint.h"

      IF (.NOT.(IOPT.EQ.1.OR.IOPT.EQ.-1.OR.IOPT.EQ.2))
     &   CALL QUIT('SQSMAT: Invalid IOPT! Must be 1 or -1 or 2 !')

      N2DIM = LRA*N*NZA
      CALL MEMGET('REAL',KEVEC,N2DIM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVAL,N,WORK,KFREE,LFREE)
C
      CALL SQSTR2(IOPT,A,N,LRA,LCA,NZA,WORK(KEVEC),WORK(KEVAL),
     &            WORK,KFREE,LFREE)
C
      CALL MEMREL('SQSTRA',WORK,KWORK,KWORK,KFREE,LFREE)      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Sqstr2 */
      SUBROUTINE SQSTR2(IOPT,A,N,LRA,LCA,NZA,EVEC,EVAL,WORK,KFREE,LFREE)
C***********************************************************************
C
C     IOPT = -1:
C     Generate the matrix  Inv{Sqrt(S)} to be used in 
C     symmetric orthonormalization
C   
C     IOPT = +1:
C     Generate the matrix Sqrt(S)
C
C     IOPT = +2:
CMI   Generate the matrix Inv(S)     
C
C     NOTE that matrix A gets destroyed in the process !
C
C     Written by T. Saue Oct 10 2004
C     Last changes: M.Ilias, Sept 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0,D0=0.0D0)
      DIMENSION A(LRA,LCA,NZA),EVEC(LRA,N,NZA),EVAL(N),WORK(*)
#include "dgroup.h"
      KFRSAV = KFREE
C
C     First diagonalize the matrix and save eigenvectors 
C
      CALL QDIAG(NZA,N,A,LRA,LCA,EVAL,1,EVEC,LRA,N,
     &           WORK(KFREE),LFREE,IERR)
      N2DIM = LRA*N*NZA
      CALL MEMGET('REAL',KBUF,N2DIM,WORK,KFREE,LFREE)
      CALL DCOPY(N2DIM,EVEC,1,WORK(KBUF),1)
C
C     Then form Lowdin matrix for canonical orthonormalization
C     (here one should perhaps introduce a check on linear
C     dependency)

      IF(IOPT.EQ.-1) THEN
        DO J = 1,N
          FAC = D1/SQRT(EVAL(J))
          DO JZ = 1,NZA
            CALL DSCAL(N,FAC,EVEC(1,J,JZ),1)
          ENDDO
        ENDDO
      ELSE IF(IOPT.EQ.1) THEN
        DO J = 1,N
          FAC = SQRT(EVAL(J))
          DO JZ = 1,NZA
            CALL DSCAL(N,FAC,EVEC(1,J,JZ),1)
          ENDDO
        ENDDO
      ELSE IF(IOPT.EQ.2) THEN
CMI     ... this is the pure inversion matrix ...
        DO J = 1,N
          FAC = D1/(EVAL(J))
          DO JZ = 1,NZA
            CALL DSCAL(N,FAC,EVEC(1,J,JZ),1)
          ENDDO
        ENDDO
      ELSE
       CALL QUIT('SQSTR2: wrong value of IOPT !')       
      ENDIF
C

C     Finally get the full matrix
C
      CALL QGEMM(N,N,N,D1,'N','N',IQDEF,WORK(KBUF),LRA,N,NZA,
     &                    'H','N',IQDEF,EVEC      ,LRA,N,NZA,
     &                         D0,IQDEF,A,LRA,LCA,NZA)
C
      CALL MEMREL('SQSTR2',WORK,KFRSAV,KFRSAV,KFREE,LFREE)              
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prdiag */
      SUBROUTINE PRDIAG(A,LRA,N)
C***********************************************************************
C
C     Print diagonal elements of a matrix
C
C     Written by T. Saue Oct 25 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(LRA,N)
C
      WRITE(LUPRI,'(A,I5)') 'PRDIAG output...',N
      DO I = 1,N
        WRITE(LUPRI,*) I,A(I,I)
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck absind */
      SUBROUTINE ABSIND(N,ARR,INDX)
C***********************************************************************
C
C     This subroutine goes through array ARR and returns the ordering
C     index array INDX giving indices from smallest to largest
C     absolute value
C     Based on Luuk's routine INDEXX
C
C***********************************************************************
      INTEGER N,INDX(N),M,NSTACK
      REAL*8 ARR(N)
      PARAMETER (M=7,NSTACK=50)
      INTEGER I,INDXT,IR,ITEMP,J,JSTACK,K,L,ISTACK(NSTACK)
      REAL*8 A
      DO J=1,N
         INDX(J)=J
      END DO
      JSTACK=0
      L=1
      IR=N
 1    IF(IR-L.GE.M)THEN
         K=(L+IR)/2
         ITEMP=INDX(K)
         INDX(K)=INDX(L+1)
         INDX(L+1)=ITEMP
         IF(ABS(ARR(INDX(L+1))).GT.ABS(ARR(INDX(IR))))THEN
            ITEMP=INDX(L+1)
            INDX(L+1)=INDX(IR)
            INDX(IR)=ITEMP
         ENDIF
         IF(ABS(ARR(INDX(L))).GT.ABS(ARR(INDX(IR))))THEN
            ITEMP=INDX(L)
            INDX(L)=INDX(IR)
            INDX(IR)=ITEMP
         ENDIF
         IF(ABS(ARR(INDX(L+1))).GT.ABS(ARR(INDX(L))))THEN
            ITEMP=INDX(L+1)
            INDX(L+1)=INDX(L)
            INDX(L)=ITEMP
         ENDIF
         I=L+1
         J=IR
         INDXT=INDX(L)
         A=ABS(ARR(INDXT))
 3       CONTINUE
         I=I+1
         IF(ABS(ARR(INDX(I))).LT.A)GOTO 3
 4       CONTINUE
         J=J-1
         IF(ABS(ARR(INDX(J))).GT.A)GOTO 4
         IF(J.LT.I)GOTO 5
         ITEMP=INDX(I)
         INDX(I)=INDX(J)
         INDX(J)=ITEMP
         GOTO 3
 5       INDX(L)=INDX(J)
         INDX(J)=INDXT
         JSTACK=JSTACK+2
         IF(JSTACK.GT.NSTACK)
     $        CALL QUIT('NSTACK TOO SMALL IN INDEXX')
         IF(IR-I+1.GE.J-L)THEN
            ISTACK(JSTACK)=IR
            ISTACK(JSTACK-1)=I
            IR=J-1
         ELSE
            ISTACK(JSTACK)=J-1
            ISTACK(JSTACK-1)=L
            L=I
         ENDIF
      ELSE
         DO J=L+1,IR
            INDXT=INDX(J)
            A=ABS(ARR(INDXT))
            DO I=J-1,1,-1
               IF(ABS(ARR(INDX(I))).LE.A)GOTO 2
               INDX(I+1)=INDX(I)
            END DO
            I=0
 2          INDX(I+1)=INDXT
         END DO
         IF(JSTACK.EQ.0)RETURN
         IR=ISTACK(JSTACK)
         L=ISTACK(JSTACK-1)
         JSTACK=JSTACK-2
      ENDIF
      GOTO 1
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Factln */
      FUNCTION FACTLN(N)
C***********************************************************************
C
C     Returns ln(n!)
C
C     From "Numerical Recipies"
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(NDIM=100,DM1=-1.0D0,D1=1.0D0,D0=0.0D0)
      DIMENSION A(NDIM)
      SAVE A
      DATA A/NDIM*DM1/
      IF(N.LT.0) CALL QUIT('FACTLN: Negative factorial !')
      IF(N.LE.99) THEN
         IF(A(N+1).LT.D0) A(N+1)=GAMMLN(N+D1)
         FACTLN=A(N+1)
      ELSE
         FACTLN=GAMMLN(N+D1)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Gammln */
      FUNCTION GAMMLN(XX)
C***********************************************************************
C
C     Returns the value ln[Gamma(xx)] for xx > 0
C
C     From Numerical Recipies
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
      DIMENSION COF(6)
      SAVE COF,STP
      DATA COF,STP/76.18009172947146D0,-86.50532032941677D0,
     &             24.01409824083091D0,-1.231739572450155D0,
     &             0.1208650973866179D-2,-0.5395239384953D-5,
     &             2.5066282746310005D0/
      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO J = 1,6
         Y=Y+D1
         SER=SER+COF(J)/Y
      ENDDO
      GAMMLN=TMP+LOG(STP*SER/X)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Bico */
      FUNCTION BICO(N,K)
C***********************************************************************
C
C     Returns the binomial coefficient (n k) as a floating-point number
C
C     From "Numerical Recipies"
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      BICO=NINT(EXP(FACTLN(N)-FACTLN(K)-FACTLN(N-K)))      
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck orblo3 */
      SUBROUTINE OCCUPY(ISEL,N)
C***********************************************************************
C
C     Written by S.Dubillard Jan 21 2004
C     Last revision:
C
C     Construction of the array ISEL for the selection of the occupied
C     electronic orbital coefficients
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION ISEL(N)
      CALL QENTER('OCCUPY')
C
      DO I = 1,N
         ISEL(I) = I
      ENDDO
C
      CALL QEXIT('OCCUPY')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck matrix_inversion */
      subroutine matrix_inversion ( A, NP, IPIV, INDXR, INDXC )
C     !-------------------------------------------------------------------------
C     !
C     !	      Taken from "Numeric recipes".  The original program was
C     !       GAUSSJ which solves linear equations by the Gauss_Jordon
C     !       elimination method.  Only the parts required to invert
C     !	      matrices have been retained.
C     !
C     !	      J.P. Griffith  6/88
C     !
C     !  Added by Miro Ilias, Strasbourg, 2005
C     !
C     !  Called from SPHGEN (dirac/dirtra.F)
C     !-------------------------------------------------------------------------
#include "implicit.h" 
#include "priunit.h"  

CMI     PARAMETER (NMAX=50)

      DIMENSION A(NP,NP), IPIV(*), INDXR(*), INDXC(*)

      n = np

      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
CMI             PAUSE 'Singular matrix'
                CALL QUIT('matrix_inversion: Singular matrix!')
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE

        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
CMI     IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        IF (A(ICOL,ICOL).EQ.0.) CALL QUIT('Singular matrix.')
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MINV(AB,N,LDA,SCRATCH,DET,EPS,M,MODE)
!--------------------------------------------------------------------------------
!C
!C     A subroutine that calculates the determinant and inverse of
!C          a matrix, as well as solving systems of linear equations.
!C     Martin J. McBride.  11/25/85.
!C     General Electric CRD, Information System Operation.
!C
!C  THE COMMENTS ABOVE REFER TO AN OLD "SCIPORT" VERSION OF THIS ROUTINE,
!C  PRESUMABLY WRITTEN BY MR. MC BRIDE.  UNFORTUNATELY, IT SEEMS THAT 
!C  HE IS A BAD PROGRAMMER AND THIS ROUTINE DID NOT WORK.  MINV IS NOW
!C  A SIMPLE INTERFACE TO LINPACK ROUTINES WHICH *WORK*
!C
!CEND
!
! miro: called from routine DIIS, hsmrcc/high_sect_fscc_external/ccsd1.F;
!      routine formerly part of spin-adapted hihger sector FS CC codes
!
!--------------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AB(LDA,1),SCRATCH(1),X(2)
C
      IOFF=1+N
      CALL DGECO(AB,LDA,N,SCRATCH,X,SCRATCH(IOFF))
      IF(M.EQ.1)CALL DGESL(AB,LDA,N,SCRATCH,AB(1,N+1),0)
      IF(M.GT.1)THEN
       WRITE(6,1000)
1000   FORMAT(T3,'@MINV-F, Only one system can be solved at a time.')
       CALL QUIT('MINV: Only one system can be solved at a time.')
      ENDIF
!              ... library routine pdpack/dge.F
      IF(M.EQ.0)CALL DGEDI(AB,LDA,N,SCRATCH,X,SCRATCH(IOFF),11)
!C
!C THIS IS WRONG, BUT WHO CARES?
!C
      DET=X(1)      
C
      RETURN
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C*/ Deck Totmass */
      FUNCTION TOTMASS()
C***********************************************************************
C
C  Calculate total mass of molecule
C
C  Rip-off from LNROUT in DALTON.
C  However, somewhat simplified: As a first go, this routine
C  provides the mass of the molecule using the most abundant isotopes.
C
C  Written by T. Saue Feb 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
#include "mxcent.h"
C
#include "nuclei.h"
      TOTMASS = D0
      DO IATOM = 1, NUCIND
        NCHARG = IZATOM(IATOM)
        IF (NCHARG .NE. 0 .AND. .NOT. NOORBT(IATOM)) THEN
          TOTMASS = TOTMASS
     &            + NUCDEG(IATOM)*DISOTP(NCHARG,ISOTOP(IATOM),'MASS')
        ENDIF
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MODEHARM(OMEGA,V2,XOPT,UM,C,NP,NORDER,ILOGG)
C***********************************************************************
C
C     Given a polynomial fit, find harmonic constant
C
C     Trond Saue, Nov 17 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
C
      PARAMETER(DM1 = -1.0D0,D1 = 1.0D0,D2 = 2.0D0)
      DIMENSION C(NP)
C
      V2 = POL2DER(NORDER,C,XOPT)
      OMEGA = SQRT(V2/UM)
      TEMP  = V2*XTJ/XTANGM10/XTANGM10
      WRITE(LUPRI,'(A,1P,E18.4,A)') 
     &   '* Force constant  :',TEMP, ' N/m'
      WRITE(ILOGG,'(A,1P,E18.4,A)') 
     &   '* Force constant  :',TEMP, ' N/m'
      WRITE(LUPRI,'(A,1P,E12.5,A)') 
     &   '* Frequency       :',OMEGA*XTHZ,' Hz',
     &   '                   ',OMEGA*XTKAYS,' cm-1'
      WRITE(ILOGG,'(A,1P,E12.5,A)') 
     &   '* Frequency       :',OMEGA*XTHZ,' Hz',
     &   '                   ',OMEGA*XTKAYS,' cm-1'
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MODEANHARM(WXE,OMEGA,XOPT,UM,C,NP,NORDER,ILOGG)
C***********************************************************************
C
C     Given a polynomial fit, find anharmonic constant
C
C     Trond Saue, Nov 17 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
C
      PARAMETER(DM1 = -1.0D0,D1 = 1.0D0,D2 = 2.0D0,D4=4.0D0)
C
      DIMENSION C(NP)
C
      V2 = POLNDER(NORDER,C,XOPT,2)
      V3 = POLNDER(NORDER,C,XOPT,3)
      V4 = POLNDER(NORDER,C,XOPT,4)

      TEMP = 7.0D0*V3*V3/V2/9.0D0
      T2 = 100.0D0*(TEMP - V4)/V4
      WRITE(ILOGG,'(A,D20.10)')
     &   '* Third derivative  (au)             : ', V3,
     &   '* Fourth derivative (au)             : ', V4
      WRITE(ILOGG,'(A,D20.10, F9.2,A)')
     &   '* Fourth derivative (Morse estimate) : ', TEMP,T2,' \\%'
      TEMP = D1/UM/OMEGA
      T2   = TEMP*TEMP
      WXE = (5.0D0*TEMP*T2*V3*V3/48.0D0/OMEGA)-(V4*T2/16.0D0)
      WRITE(LUPRI,'(A,1P,E12.5,A)') 
     &   '* Omega*x_e       :',WXE*XTKAYS,' cm-1'
      WRITE(ILOGG,'(A,1P,E12.5,A)') 
     &   '* Omega*x_e       :',WXE*XTKAYS,' cm-1'
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MODEDISP(X,NP,XOPT,UM,OMEGA,ILOGG)
C***********************************************************************
C
C     Given a reduced mass and harmonic frequency, check 
C     for points outside mean displacement in
C     harmonic ground state
C
C     Trond Saue, Nov 17 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
C
      PARAMETER(DM1 = -1.0D0,D1 = 1.0D0,D2=2.0D0)
      DIMENSION X(NP)
C
      XDEV=D1/SQRT(UM*OMEGA)
      XMIN=XOPT-XDEV
      XMAX=XOPT+XDEV
      TEMP=XDEV*XTANG
      WRITE(LUPRI,'(A,F10.4,A/A,A1,F10.4,A1,F10.4,A1)') 
     &   '* Mean displacement in harmonic ground state:',TEMP,
     &   ' Angstroms',
     &   '  corresponding to interval: ',
     &   '[',XMIN*XTANG,',',XMAX*XTANG,']'
      WRITE(ILOGG,'(A,F10.4,A/A,A1,F10.4,A1,F10.4,A1)') 
     &   '* Mean displacement in harmonic ground state:',TEMP,
     &   ' Angstroms',
     &   '  corresponding to interval: ',
     &   '[',XMIN*XTANG,',',XMAX*XTANG,']'
      IOUT=0
      DO I = 1,NP
        IF(X(I).LT.XMIN.OR.X(I).GT.XMAX) IOUT=IOUT+1
      ENDDO
      IF(IOUT.GT.0) THEN
        WRITE(LUPRI,'(A,I5,A/A)')
     &   '* WARNING : ',IOUT,' points lie outside interval...',
     &   'This may reduce the quality of your spectroscopic constants.'
        WRITE(ILOGG,'(A,I5,A/A)')
     &   '* WARNING : ',IOUT,' points lie outside interval...',
     &   'This may reduce the quality of your spectroscopic constants.'
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MODEDISS(OMEGA,WXE,ILOGG)
C***********************************************************************
C
C     For a polynomial fit, estimate dissociation energy
C     from a Morse model
C
C     Trond Saue, Nov 17 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
C
      PARAMETER(D2=2.0D0,D4=4.0D0)
C
      DE = OMEGA*OMEGA/WXE/D4
      WRITE(LUPRI,'(A,2X,E10.4,A)') '* Estimated D_e:',DE*XTEV,   ' eV'
      WRITE(LUPRI,'(A,2X,E10.4,A)')
     &  '* Estimated D_e:',DE*XKJMOL, ' kJ/mol'
      WRITE(LUPRI,'(A,2X,E10.4,A)')
     & '* Estimated D_e:',DE*XTKAYS, ' cm-1'
      WRITE(ILOGG,'(A,2X,E10.4,A)') '* Estimated D_e:',DE*XTEV, ' eV'
      D0 = DE-(OMEGA/D2)+WXE/D4
      WRITE(LUPRI,'(A,2X,E10.4,A)') '* Estimated dissociation energy:',
     &       D0*XTEV, ' eV'
      WRITE(LUPRI,'(A,2X,E10.4,A)') '* Estimated dissociation energy:',
     &       D0*XKJMOL, ' kJ/mol'
      WRITE(LUPRI,'(A,2X,E10.4,A)') '* Estimated dissociation energy:',
     &       D0*XTKAYS  , ' cm-1'
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MODEFIT(NP,X,Y,A,B,C,D,MORDER,NORDER,ILOGG)
C***********************************************************************
C
C     From NP points along a normal coordinate, do a polynomial fit.
C     MORDER = maximum order of polynomial fit
C     NORDER = actual order of polynomial fit
C
C     Written by T. Saue Nov 17 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
C
      DIMENSION X(NP),Y(NP),A(NP,*),B(NP),C(NP),D(NP)
C
      WRITE(LUPRI,'(A)') 
     &  'Polynomial fit: Give order of polynomial '
      READ(LUSTDIN,*) NORDER
      IF(NORDER.LT.2) THEN
        WRITE(LUPRI,'(A)') 
     &    '* This program requires at least a quadratic fit.'
        STOP 'Too low order fit !'
      ENDIF
      IF(NORDER.LT.4) THEN
        WRITE(LUPRI,'(A)') 
     &    '* WARNING: a quartic fit is needed for full functionality.'
      ENDIF
      IF(NP.LE.NORDER) THEN
        WRITE(LUPRI,'(A,I3)') 
     +    'Too few points for a polynomial fit of order',NORDER
        STOP 'Too few points !'
      ENDIF
      IF(NORDER.GT.MORDER) THEN
        STOP 'Order to large ....'
      ENDIF
      NDIM = NORDER + 1
      CALL POLSVD(NDIM,NP,A,B,X,Y,C,D,CHISQ,ISKIP)
C* Estimate fit:
      WRITE(LUPRI,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(LUPRI,'(3X,A)') '* Coefficients:'
      WRITE(LUPRI,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &      C(I),I=1,NDIM)  
C
      WRITE(ILOGG,'(72A1)') ('=',I=1,72)
      WRITE(ILOGG,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(ILOGG,'(3X,A)') '* Coefficients:'
      WRITE(ILOGG,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &      C(I),I=1,NDIM)  
      WRITE(LUPRI,'(8X,A,4X,2A)')
     &   'X',' Predicted Y',' Relative error'
      DO I = 1,NP
        YP = POLVAL(NORDER,C,X(I))
        DEV = (YP-Y(I))/Y(I)
        WRITE(LUPRI,'(3X,F6.3,4X,1P,E20.12,E12.4)') XTANG*X(I),YP,DEV
      ENDDO
      WRITE(LUPRI,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(LUPRI,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP
      WRITE(LUPRI,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP
      WRITE(ILOGG,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MODEMIN(XOPT,NP,X,Y,C,NORDER,NITER,ILOGG)
C***********************************************************************
C
C     Given a polynomial fit, find local minimum,
C     using Newton-Raphsons method
C
C     Trond Saue, Nov 17 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
C
      PARAMETER(DTOL = 1.0D-4)
C
      DIMENSION X(NP),Y(NP),C(NP)
C
C* Find smallest y-value and take corresponding x-value as start value
C* for Newton-Raphsons method      
      IYOPT = IDMIN(NP,Y,1)
      XOPT = X(IYOPT)
      CALL NEWRAP(NORDER,DTOL,C,XOPT,NITER,IERR)
      IF(IERR.EQ.1) THEN
        WRITE(ILOGG,'(A)') 'Newtons method failed'
        STOP 'Did not find minimum. Newtons method failed...'
      ELSE
        YOPT = POLVAL(NORDER,C,XOPT)
        TEMP = XOPT*XTANG
        WRITE(LUPRI,'(A,F18.5,A,F18.5,A)') 
     &    '* Local minimum   :',TEMP,' Angstroms = ',XOPT,' Bohrs'
        WRITE(ILOGG,'(A,F18.5,A)') 
     &    '* Local minimum   :',TEMP,' Angstroms'
        WRITE(LUPRI,'(A,1P,E18.10,A)') 
     &    '* Expected energy :',YOPT,' Hartrees'
        WRITE(ILOGG,'(A,1P,E18.10,A)') 
     &    '* Expected energy :',YOPT,' Hartrees'
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INDEXFILL(INDEX,N,ISTART,INC)
C***********************************************************************
C
C     Fill an index according to first element and increment, see below
C
C***********************************************************************
      DIMENSION INDEX(N)
      II = ISTART
      DO I = 1,N
        INDEX(I)=II
        II=II+INC
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION SYMCHECK(A,LR,LC,N,NZ)
C***********************************************************************
C
C     Check symmetry aboyut the diagonal of a matrix; 
C     uses absolute values  in order to work for both 
C     symmetric and anti-symmetric matrices
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(LR,LC,NZ)
      SYMCHECK = 0.0D0
      DO IZ = 1,NZ
        DO J = 1,N-1
          DO I = J+1,N
            SYMCHECK = SYMCHECK +
     &        (DABS(A(I,J,IZ))-DABS(A(J,I,IZ)))
          ENDDO
        ENDDO
      ENDDO
      SYMCHECK=SYMCHECK/DFLOAT(NZ*N*(N-1)/2)
      RETURN
      END
C  /* Deck indexi */
      SUBROUTINE INDEXI(N,IARR,INDX)
C***********************************************************************
C
C     This subroutine goes through an integer array IARR and returns the ordering
C     index array INDX giving indices from smallest to largest
C
C     The array IARR is not touched.
C
C***********************************************************************
      INTEGER N,INDX(N),M,NSTACK,IARR(N)
      PARAMETER (M=7,NSTACK=50)
      INTEGER I,IA,INDXT,IR,ITEMP,J,JSTACK,K,L,ISTACK(NSTACK)
      DO J=1,N
         INDX(J)=J
      END DO
      JSTACK=0
      L=1
      IR=N
 1    IF(IR-L.GE.M)THEN
         K=(L+IR)/2
         ITEMP=INDX(K)
         INDX(K)=INDX(L+1)
         INDX(L+1)=ITEMP
         IF(IARR(INDX(L+1)).GT.IARR(INDX(IR)))THEN
            ITEMP=INDX(L+1)
            INDX(L+1)=INDX(IR)
            INDX(IR)=ITEMP
         ENDIF
         IF(IARR(INDX(L)).GT.IARR(INDX(IR)))THEN
            ITEMP=INDX(L)
            INDX(L)=INDX(IR)
            INDX(IR)=ITEMP
         ENDIF
         IF(IARR(INDX(L+1)).GT.IARR(INDX(L)))THEN
            ITEMP=INDX(L+1)
            INDX(L+1)=INDX(L)
            INDX(L)=ITEMP
         ENDIF
         I=L+1
         J=IR
         INDXT=INDX(L)
         IA=IARR(INDXT)
 3       CONTINUE
         I=I+1
         IF(IARR(INDX(I)).LT.IA)GOTO 3
 4       CONTINUE
         J=J-1
         IF(IARR(INDX(J)).GT.IA)GOTO 4
         IF(J.LT.I)GOTO 5
         ITEMP=INDX(I)
         INDX(I)=INDX(J)
         INDX(J)=ITEMP
         GOTO 3
 5       INDX(L)=INDX(J)
         INDX(J)=INDXT
         JSTACK=JSTACK+2
         IF(JSTACK.GT.NSTACK)
     $        CALL QUIT('NSTACK TOO SMALL IN INDEXI')
         IF(IR-I+1.GE.J-L)THEN
            ISTACK(JSTACK)=IR
            ISTACK(JSTACK-1)=I
            IR=J-1
         ELSE
            ISTACK(JSTACK)=J-1
            ISTACK(JSTACK-1)=L
            L=I
         ENDIF
      ELSE
         DO J=L+1,IR
            INDXT=INDX(J)
            IA=IARR(INDXT)
            DO I=J-1,1,-1
               IF(IARR(INDX(I)).LE.IA)GOTO 2
               INDX(I+1)=INDX(I)
            END DO
            I=0
 2          INDX(I+1)=INDXT
         END DO
         IF(JSTACK.EQ.0)RETURN
         IR=ISTACK(JSTACK)
         L=ISTACK(JSTACK-1)
         JSTACK=JSTACK-2
      ENDIF
      GOTO 1
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck cleanmat */
      SUBROUTINE CLEANMAT(A,DTOL,N,LR,LC,NZ)
C***********************************************************************
C
C     Subroutine CLEANMAT will set matrix elements below the given 
C     threshold to zero
C     Written by T. Saue Oct. 2 2013
C***********************************************************************
#include "implicit.h"
      DIMENSION A(LR,LC,NZ)
      PARAMETER (D1=1.0D00, D0=0.0D00)
      DO IZ = 1,NZ
        DO J = 1,N
          DO I = 1,N
            IF(ABS(A(I,J,IZ)).LT.DTOL) A(I,J,IZ) = D0
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Absolute_Overlap */
      SUBROUTINE ABSOLUTE_OVERLAP(CMO,NBAS,NVEC,ISEL,NSEL,
     &                            CREF,NCOL,IREF,SMAT,AMAT,
     &                            IPQ,NZ,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate overlap between vector set; 
C     for readability report square of overlap
C
C     Written by T. Saue Oct 3 2013
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
      DIMENSION CMO(NBAS,NVEC,NZ),CREF(NBAS,NCOL,NZ),
     &          SMAT(NBAS,NBAS),AMAT(NSEL,NSEL,NZ),IPQ(NZ),WORK(*)
C.....Calculate overlap between reference and selected orbitals
      CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NSEL,NSEL,
     &         SMAT,NBAS,NBAS,1,IPQ,
     &         AMAT,NSEL,NSEL,NZ,IPQ,
     &         CREF(1,IREF,1),NBAS,NCOL,NZ,IPQ,
     &         CMO(1,ISEL,1),NBAS,NVEC,NZ,IPQ,
     &         WORK(KFREE),LFREE,IPRINT)
      DO J = 1, NSEL
        DO I = 1,NSEL
           AMAT(I,J,1) = AMAT(I,J,1)*AMAT(I,J,1)
        ENDDO
      ENDDO
      DO IZ = 2,NZ
        DO J = 1, NSEL
          DO I = 1,NSEL
             AMAT(I,J,1) = AMAT(I,J,1)+AMAT(I,J,IZ)*AMAT(I,J,IZ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END

      SUBROUTINE NZ_MTX_AMPLITUDES(AMAT,AMAT1,NTOT,NZ)
!------------------------------------------------------------------------
! Get amplitudes of the full quaternion matrix AMAT into AMAT1
!
! AMAT1(i,j,1)=sqrt(AMAT(i,j,1)+...+AMAT(i,j,NZ))
!
! Written by Miro Ilias, GSI, 2016
!------------------------------------------------------------------------
      integer, intent(in)  :: NTOT,NZ
      real*8,  intent(in)  :: AMAT(NTOT,NTOT,NZ)
      real*8,  intent(out) :: AMAT1(NTOT,NTOT,1)
      integer              :: I,J,IZ

      DO J = 1, NTOT
        DO I = 1,NTOT
!         ...accumulate the real part
           AMAT1(I,J,1) = (AMAT(I,J,1)*AMAT(I,J,1))
        ENDDO
      ENDDO

      DO IZ = 2,NZ
        DO J = 1, NTOT
          DO I = 1,NTOT
!           ... accumulate quaternion-imaginary parts
             AMAT1(I,J,1) = AMAT1(I,J,1)+(AMAT(I,J,IZ)*AMAT(I,J,IZ))
          ENDDO
        ENDDO
      ENDDO
!    finally, get the amplitude
      DO J = 1, NTOT
        DO I = 1,NTOT
           AMAT1(I,J,1) = DSQRT(AMAT1(I,J,1))
        ENDDO
      ENDDO
      END

      SUBROUTINE SUM_NZ_MATRIX(AMAT,NROW,NCOL,LRQ,LCQ,NZ,
     &                         AMAT_SUM_NZ,AMAT_SUM)
!----------------------------------------------------------------------------
! Sum all elements in the entering quaternion matrix AMAT into
! quaternion AMAT_SUM_NZ and real amplitude AMAT_SUM.
!
! Written by Miro Ilias, GSI, 2016
!----------------------------------------------------------------------------
      integer, intent(in)    :: NROW,NCOL,LRQ,LCQ,NZ
      real*8,  intent(in)    :: AMAT(LRQ,LCQ,NZ)
      real*8,  intent(out)   :: AMAT_SUM_NZ(1:4),AMAT_SUM
      integer                :: I,J,IZ
      AMAT_SUM_NZ(1)=0.0D0; AMAT_SUM_NZ(2)=0.0D0
      AMAT_SUM_NZ(3)=0.0D0; AMAT_SUM_NZ(4)=0.0D0 
      DO IZ=1, NZ
        DO I=1, NROW
          DO J=1, NCOL
            AMAT_SUM_NZ(IZ)=AMAT_SUM_NZ(IZ)+AMAT(I,J,IZ)
          ENDDO
        ENDDO
      ENDDO
! get the amplitude ...
      AMAT_SUM =
     &      DSQRT((AMAT_SUM_NZ(1)*AMAT_SUM_NZ(1)) + 
     &            (AMAT_SUM_NZ(2)*AMAT_SUM_NZ(2)) +
     &            (AMAT_SUM_NZ(3)*AMAT_SUM_NZ(3)) +
     &            (AMAT_SUM_NZ(4)*AMAT_SUM_NZ(4)))
      END

      SUBROUTINE SUM_RE_MATRIX(AMAT,NROW,NCOL,LRQ,LCQ,AMAT_SUM)
!----------------------------------------------------------------------------
! Sum all elements in the entering real matrix AMAT, containing real 
! amplitutes, into AMAT_SUM.
!
! Written by Miro Ilias, GSI, 2016
!----------------------------------------------------------------------------
      integer, intent(in)  :: NROW,NCOL,LRQ,LCQ
      real*8,  intent(in)  :: AMAT(LRQ,LCQ,1)
      real*8,  intent(out) :: AMAT_SUM
      integer              :: I,J
      AMAT_SUM=0.0D0
      DO I=1, NROW
      DO J=1, NCOL
         AMAT_SUM=AMAT_SUM+AMAT(I,J,1)
      ENDDO
      ENDDO
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck iicopy */
      SUBROUTINE IICOPY(NDIM,IA,IA1,IASTEP,IB,IBSTEP,IB1)
C***********************************************************************
C
C     WRapper routine allowing use of ICOPY for WORK
C
C***********************************************************************
      implicit none
      integer, intent(in) :: IA(*),NDIM,IA1,IASTEP,IBSTEP,IB1
      integer, intent(out):: IB(*)
      CALL ICOPY(NDIM,IA(IA1),IASTEP,IB(IB1),IBSTEP)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CHECK_FOR_OPEN_FILES
C***********************************************************************
C
C     Check for open files within some range of unit numbers
C
C***********************************************************************
      INTEGER I
      LOGICAL ITSOPEN
      DO I = 1,300
         INQUIRE(UNIT=I,OPENED=ITSOPEN)
         IF(ITSOPEN) THEN
           WRITE(6,*) 'Open file. Unit = ',I
         ENDIF
      ENDDO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LEVI_CEVITA(I,J,K,IJK)
C***********************************************************************
C
C     Given indices i and j, give k and e_ijk
C
C     Written by Trond Saue Oct 25 2018      
C***********************************************************************
      implicit none
      integer, intent(in)  :: I,J
      integer, intent(out) :: K,IJK
      INTEGER :: IJKTAB(3,3),KTAB(3,3)
      SAVE IJKTAB,KTAB
      DATA IJKTAB/ 0,-1, 1, 1, 0,-1,-1, 1, 0/
      DATA   KTAB/ 0, 3, 2, 3, 0, 1, 2, 1, 0/
      IF(I.LT.1.OR.I.GT.3) CALL QUIT('LEVI-CEVITA: Not valid I!')
      IF(J.LT.1.OR.J.GT.3) CALL QUIT('LEVI-CEVITA: Not valid J!')
      K   = KTAB(I,J)
      IJK = IJKTAB(I,J)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE NOISE(A,N,FACTOR)
C***********************************************************************
C
C     Input is an array A of length N. On output noise has been added.
C     FACTOR controls the magnitude of the noise ("noise level")
C     Since RANDOM_NUMBER returns numbers in the interval (0,1],
C     we extend to (-1,1].      
C
C     Written by Trond Saue Nov 7 2019
C
C***********************************************************************
      implicit none
      integer, intent(in)   :: N
      real*8, intent(in)    :: factor
      real*8, intent(inout) :: A(N)
      integer               :: I
      real*8                :: R

      DO I = 1,N
        CALL RANDOM_NUMBER(r)
        A(I) = A(I) + FACTOR*(2.0D00*R-1.0D00)
      ENDDO
      RETURN
      END
C/* Deck Usym */
      SUBROUTINE USYM(A,LDA,N)
C***********************************************************************
C
C     Symmetrize an upper triangular matrix
C
C***********************************************************************
       INTEGER LDA,N,I,J
       REAL*8 A(LDA,N)
       DO I=2,N
          DO J=1,I
            A(I,J) = A(J,I)
          ENDDO
       ENDDO
       END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Fdiag */
      SUBROUTINE FDIAG(A,LDA,N)
C***********************************************************************
C
C     Force diagonal form by zeroing out all non-diagonal elements
C
C***********************************************************************
       INTEGER LDA,N,I,J
       REAL*8 A(LDA,N)
       DO I=1,N
          DO J=1,N
            IF(I.NE.J) A(I,J) = 0.0D0
          ENDDO
       ENDDO
       END

