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
C  /* Deck flshfo */
      SUBROUTINE FLSHFO (IUNIT)
C
C *** THIS SUBROUTINE IS SYSTEM DEPENDENT ***
C
C     Flush formatted output unit (empty buffers).
C     If no flush utility, this is achieved by
C     CLose and reOPen Formatted Output
C
C Written 21-Nov-1983 by Hans Jorgen Aa. Jensen in Uppsala, Sweden.
C Last revision 16-Jul-1984 hjaaj / 30-Oct-1984 hjaaj (extendsize)
C 10-Feb-1989 hjaaj, renamed CLOPFO to FLSHFO
C
C Calls to this subroutine makes it possible to read the output
C up to the moment of the last call while the program continues
C executing (provided the computer allows shared access).
C This subroutine may be a dummy routine.
C
C     Force transfer of all buffered output to the file or device
C     associated with logical unit IUNIT.
C
C 21-02-2007  LV : Removed outdated machines, after 24 years (!) the
C                  whole routine can now probably be deleted (we always call FLUSH)
      CALL FLUSH(IUNIT)
      RETURN
      END
C  /* Deck ordrss */
      SUBROUTINE ORDRSS(EVEC,EVAL,ISS,N,NEVEC)
C
C 920729-hjaaj (based on ORDER)
C
C Purpose: order the N values in EVAL and their associated vectors
C          in EVEC so EVAL(i+1) .ge. EVAL(i),
C          but only within the class of vectors having the
C          same value in the ISS array (which could be the
C          supersymmetry of the orbital).
C
#include "implicit.h"
      DIMENSION EVEC(*),EVAL(*),ISS(*)
      IF (N.LE.1) RETURN
      IN = 1
      DO 10 I=1,N-1
        EMIN = EVAL(I)
        IMIN = I
        ISSI = ISS(I)
        DO 20 J=I+1,N
          IF (ISS(J) .NE. ISSI) GO TO 20
          IF (EVAL(J) .LT. EMIN) THEN
            EMIN = EVAL(J)
            IMIN = J
          ENDIF
   20   CONTINUE
        IF (IMIN.NE.I) THEN
          EVAL(IMIN)=EVAL(I)
          EVAL(I)=EMIN
          IF (NEVEC .GT. 0) THEN
            CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*NEVEC+1),1)
          ENDIF
        ENDIF
        IN = IN + NEVEC
   10 CONTINUE
      RETURN
      END
C  /* Deck ord2ss */
      SUBROUTINE ORD2SS(EVEC,EVAL,ISS,N,NEVEC)
C
C Purpose: order the N values in EVAL and their associated vectors
C          in EVEC so EVAL(i+1) .le. EVAL(i) using the infomation
C          in ISS
C          (this is opposite order of "ORDRSS")
C
#include "implicit.h"
      DIMENSION EVEC(*),EVAL(*),ISS(*)
      IF (N.LE.1) RETURN
      IN = 1
      DO 10 I=1,N-1
         EMAX = EVAL(I)
         IMAX = I
         ISSI = ISS(I)
         DO 20 J=I+1,N
            IF (ISS(J) .NE. ISSI) GO TO 20
            IF (EVAL(J) .GT. EMAX) THEN
               EMAX = EVAL(J)
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
      RETURN
      END
C  /* Deck order */
      SUBROUTINE ORDER(EVEC,EVAL,N,NEVEC)
C
C Purpose: order the N values in EVAL and their associated vectors
C          in EVEC so EVAL(i+1) .ge. EVAL(i)
C
C Revisions:
C   29-Jul-1992 hjaaj (only dswap if nevec .gt. 0)
C    2-Nov-1984 hjaaj (new parameter NEVEC, EVEC(1:NEVEC,1:N))
C   27-Oct-1984 hjaaj (reduced number of swaps)
C
#include "implicit.h"
      DIMENSION EVEC(*),EVAL(*)
      IF (N.LE.1) RETURN
      IN = 1
      DO 10 I=1,N-1
        EMIN = EVAL(I)
        IMIN = I
        DO 20 J=I+1,N
          IF (EVAL(J) .LT. EMIN) THEN
            EMIN = EVAL(J)
            IMIN = J
          ENDIF
   20   CONTINUE
        IF (IMIN.NE.I) THEN
          EVAL(IMIN)=EVAL(I)
          EVAL(I)=EMIN
          IF (NEVEC .GT. 0) THEN
            CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*NEVEC+1),1)
          ENDIF
        ENDIF
        IN = IN + NEVEC
   10 CONTINUE
      RETURN
      END
C  /* Deck trace */
      SUBROUTINE TRACE
C
C Written 4-Dec-1983 hjaaj
C
C Major clean-up LV (feb. 2007)
C
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL QTRACE(6)
      CALL TRBK()
#endif
      CALL QTRACE(6)
      WRITE(*,*) 'This OS has no obvious system traceback facility.'
      RETURN
      END
C  /* Deck canon */
      SUBROUTINE CANON(I,J,K,L)
C
C     Reorder IJKL to canonical order for real 2-electron entities:
C       I .ge. J, K .ge. L, IJ .ge. KL
C
#include "implicit.h"
      IP=MAX(I,J)
      JP=I+J-IP
      KP=MAX(K,L)
      LP=K+L-KP
      IF (IP.GT.KP) THEN
         I=IP
         J=JP
         K=KP
         L=LP
      ELSE
         I=KP
         J=LP
         K=IP
         L=JP
         IF(I.NE.K)RETURN
         IF(J.GT.L)RETURN
         J=JP
         L=LP
      END IF
      RETURN
      END
C  /* Deck dunit */
      SUBROUTINE DUNIT(A,N)
C
C  SUBROUTINE DUNIT SETS THE REAL SQUARE MATRIX A EQUAL
C  TO A UNIT MATRIX.
C  /VER 2/ 14-Sep-1983 hjaaj
C
#include "implicit.h"
      DIMENSION A(*)
      PARAMETER (D1=1.0D00, D0=0.0D00)
C
      NN = N*N
      DO 100 I = 1,NN
         A(I) = D0
  100 CONTINUE
      N1 = N + 1
      DO 200 I = 1,NN,N1
         A(I) = D1
  200 CONTINUE
      RETURN
      END
C  /* Deck dzero */
      SUBROUTINE DZERO(DX,LENGTH)
#include "implicit.h"
C
C Last revision 5-May-1984 by Hans Jorgen Aa. Jensen
C
C   Subroutine DZERO sets a real array of length *LENGTH*
C   to zero.
C...................................................................
C     DIMENSION DX(*)
      DIMENSION DX(LENGTH)
C
      IF (LENGTH.LE.0) RETURN
C
      DO 100 I = 1,LENGTH
  100    DX(I) = 0.0D00
C
      RETURN
      END
C  /* Deck dzero8 */
      SUBROUTINE DZERO8(DX,LENGTH8)
#include "implicit.h"
C   Last revision 23-Nov-2007 by S. Knecht
C
C   Subroutine DZERO8 sets a real array of length *LENGTH8*
C   to zero.
C...................................................................
      INTEGER*8 LENGTH8
      DIMENSION DX(LENGTH8)
C
      IF (LENGTH8.LE.0) RETURN
C
      DO 100 I = 1,LENGTH8
  100    DX(I) = 0.0D00
C
      RETURN
      END
C  /* Deck isum */
      FUNCTION ISUM(N,IA,INCA)
C
C     8-Feb-1987 hjaaj
C     Sums elements of a integer vector.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     -- based on DDOT from LINPACK
C     DDOT: JACK DONGARRA, LINPACK, 3/11/78.
C
      INTEGER ISUM,  IA(*), ITEMP
      INTEGER I,INCA,JA,M,MP1,N
C
      ISUM  = 0
      IF(N.LE.0)RETURN
      ITEMP = 0
      IF(INCA.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      JA = 1
      IF(INCA.LT.0)JA = (-N+1)*INCA + 1
      DO 10 I = 1,N
        ITEMP = ITEMP + IA(JA)
        JA = JA + INCA
   10 CONTINUE
      ISUM = ITEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        ITEMP = ITEMP + IA(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         ITEMP = ITEMP     + IA(I)     + IA(I + 1)
     *         + IA(I + 2) + IA(I + 3) + IA(I + 4)
   50 CONTINUE
   60 ISUM = ITEMP
      RETURN
      END
C  /* Deck izero */
      SUBROUTINE IZERO(INTI,LENGTH)
C...................................................................
C Written 5-May-1984 by Hans Jorgen Aa. Jensen
C
C   Subroutine IZERO sets an integer array of length *LENGTH*
C   to zero.
C...................................................................
      INTEGER INTI(*),LENGTH
C
      IF (LENGTH.LE.0) RETURN
C
      DO 100 I=1,LENGTH
  100    INTI(I) = 0
C
      RETURN
      END
C  /* Deck izero8 */
      SUBROUTINE IZERO8(INTI8,LENGTH8)
C...................................................................
C   modified routine IZERO - S. Knecht - 19-Feb-2008
C
C   Subroutine IZERO sets an integer*8 array of length *LENGTH8*
C   to zero.
C...................................................................
      INTEGER*8 INTI8(*),LENGTH8
C
      IF (LENGTH8.LE.0) RETURN
C
      DO 100 I=1,LENGTH8
  100    INTI8(I) = 0
C
      RETURN
      END
C  /* Deck imino8 */
      SUBROUTINE IMINO8(INTI8,LENGTH8)
C...................................................................
C   modified routine IZERO - S. Knecht - 19-Feb-2008
C
C   Subroutine IMIN08 sets an integer*8 array of length *LENGTH8*
C   to minus one.
C...................................................................
      INTEGER*8 INTI8(*),LENGTH8, IMIN1
C
      IF (LENGTH8.LE.0) RETURN
      IMIN1 = - 1
C
      DO 100 I=1,LENGTH8
  100    INTI8(I) = IMIN1
C
      RETURN
      END


C  /* Deck jaco */
      SUBROUTINE JACO (F,V,NB,NMAX,NROWV,BIG,JBIG)
C
C Revisions:
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
#include "implicit.h"
      DIMENSION F(*),V(*)
      DIMENSION BIG(*) ,JBIG(*)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, ROOT2 = 0.707106781186548D0)
      PARAMETER(DP5 = 0.5D0, D1P5 = 1.5D0, D1P375 = 1.375D0,
     *          D3P875 = 3.875D0, DP25 = 0.25D0)
#include "thrzer.h"
      DATA C1,C2,C3,C4,C5,C6/THRZER,THRZER,1.D-20,1.D-14,1.D-9,1.D-5/
Cwas: DATA C1,C2,C3,C4,C5,C6/1.D-12,1.D-12,1.D-20,1.D-14,1.D-9,1.D-5/
      IF (NB.LE.1 .OR. NMAX.LE.0) RETURN
Cwas: IF (NB.EQ.1) RETURN !900727-hjaaj
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
         JAA=(IA-1)*NROWV
         JBB=(IB-1)*NROWV
#include "ivdep.h"
         DO 400 I=1,NROWV
            T=V(JBB+I)*SSX
            V(JBB+I)=-V(JAA+I)*SSX + V(JBB+I)*CX
  400       V(JAA+I)= V(JAA+I)*CX + T
      GO TO 410
  420 RETURN
      END
C  /* Deck norm */
      SUBROUTINE NORM(S,VC,N,M,W,THNORM,IRETUR)
C
C revised 14-May-1985 hjaaj (call MPAPV instead of CNTRC)
C
C     COMPUTES SCHMIDT-ORTHONORMALIZED SET OF VECTORS
C         CALLING SEQUENCE PARAMETERS ARE AS FOLLOWS
C              S    METRIC MATRIX STORED AS LOWER TRIANGLE (R*8)
C              VC   LOCATION OF ORIGINAL NON-ORTHONORMAL VECTORS (R*8)
C                   FINAL ORTHONORMALIZED VECTORS REPLACE ORIGINAL SET
C              N    DIMENSION OF BASIS SET (I*4)
C              M    NUMBER OF VECTORS TO BE ORTHONORMALIZED (I*4)
C              W    TEMPORARY WORKING AREA OF 2*N WORDS (R*8)
C         RETURNS
C              NORMAL RETURN ORTHONORMALIZED SET OBTAINED
C              RETURN 1       INITIAL VECTORS AT VC LINEARLY
C                             DEPENDENT WITHIN THRES HOLD (THNORM)
C         AUXILLIARY ENTRY
C
#include "implicit.h"
      DIMENSION S(*), VC(*), W(*)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
      IRETUR=0
C
C     N = 1 special case
C
      IF (N .EQ. 1) THEN
         IF (VC(1)*VC(1) .LT. THNORM) THEN
            IRETUR=1
         ELSE
            VC(1) = D1/SQRT(S(1))
         END IF
         RETURN
      END IF
C
C     BEGIN OUTERMOST LOOP OVER TRIAL VECTOR SET
C
      NP1 = N+1
      TNRML = THNORM
      IVCI = 1
      DO 20 I=1,M
         CALL MPAPV(N,S,VC(IVCI),W)
C        CALL CNTRC(S,VC(IVCI),W,N)
         T = D0
         TNORM= D0
         IVCJ = 1
         IVCW = N
C
C     BEGIN COEFFICIENTS AND NORMALIZATION LOOP
C
         DO 5 J=1,I
            TNORM = TNORM - T*T
            T = DDOT(N,VC(IVCJ),1,W(1),1)
            IVCJ = IVCJ + N
            IVCW = IVCW + 1
    5       W(IVCW) = -T
         IF (TNORM + T .LT. TNRML) THEN
            IRETUR = 2
            RETURN
         END IF
         TNORM = D1/SQRT(TNORM + T)
         W(IVCW) = D1
C
C     BEGIN REPLACE VC(I) OUTER LOOP
C
         DO 15 K=1,N
            IVCW = NP1
C
C     BEGIN REPLACE COMPONENT OF VC(I) LOOP
C
            T = D0
            DO 10 J=K,IVCI,N
               T = VC(J)*W(IVCW) + T
   10          IVCW = IVCW + 1
            VC(IVCI) = T*TNORM
   15       IVCI = IVCI + 1
   20 CONTINUE
      RETURN
      END

C  /* Deck mpapv */
      SUBROUTINE MPAPV(N,AP,VIN,VOUT)
C
C  28-Apr-1985 hjaaj
C
C  VOUT = AP * VIN
C
C  where VOUT and VIN are N vectors
C  and AP is a packed, symmetric N by N matrix.
C
CMI ... copied from pdpack/arhpack.F due to linking problems 
CMI    with DIRAC gp/ routines 
CMI    in ../../gp/qpack.o ../../gp/gptrygve.o  ../../gp/gphjj.o
C
#include "implicit.h"
      DIMENSION AP(*),VIN(*),VOUT(*)
      PARAMETER (D0 = 0.0D0)
C
      KI = 0
      DO 300 I = 1,N
         VINI = VIN(I)
         SUM = D0
         DO 100 J = 1,I-1
            SUM = SUM + AP(KI+J) * VIN(J)
            VOUT(J) = VOUT(J) + AP(KI+J) * VINI
  100    CONTINUE
         KI = KI + I
         VOUT(I) = SUM + AP(KI) * VINI
  300 CONTINUE
C
      RETURN
      END

C  /* Deck defopn */
      SUBROUTINE DEFOPN(LU)
C
C 29-Sep-1986 Hans Joergen Aa. Jensen
C
C DEFAULT OPEN :
C
C     OPEN(LU,STATUS='UNKNOWN',FORM='UNFORMATTED') with
C     default file name, following CWBOPN by C.W.Bauschlicher.
C
C     This is especially useful for CRAY machines which do not
C     allow for default names in the OPEN statement (i.e.
C     FILE='name' must be specified if you use OPEN).
C
      REWIND (LU)
      READ   (LU,ERR=1,END=1) A
    1 REWIND (LU)
      RETURN
      END
C  /* Deck mollab */
      SUBROUTINE MOLLAB(A,LU,LUERR)
C
C  16-Jun-1986 hjaaj
C  (as SEARCH but CHARACTER*8 instead of REAL*8 labels)
C
C  Purpose:
C     Search for MOLECULE labels on file LU
C
      CHARACTER*8 A, B(4), C
      DATA C/'********'/
      IRDERR = 0
    1 READ (LU,END=3,ERR=6,IOSTAT=IOSVAL) B
      IRDERR = 0
      IF (B(1).NE.C) GO TO 1
      IF (B(4).NE.A) GO TO 1
      IF (LUERR.LT.0) LUERR = 0
      RETURN
C
    3 IF (LUERR.LT.0) THEN
#if defined (VAR_MFDS)
!C 880916-hjaaj -- attempt to fix an IBM problem
!C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
!C     of FTxxF001), backspace makes LU ready for append.
!C 940510-hjaaj: same change for Cray's multifile datasets
         BACKSPACE LU
#endif
         LUERR = -1
         RETURN
      ELSE
         WRITE(LUERR,4)A,LU
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR (MOLLAB) MOLECULE label not found on file')
      END IF
    4 FORMAT(/' *** ERROR (MOLLAB), MOLECULE label ',A8,
     *        ' not found on unit',I4)
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      IF (LUERR.LT.0) THEN
         LUERR = -2
         RETURN
      ELSE
         WRITE (LUERR,7) LU,A,IOSVAL
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR (MOLLAB) error reading file')
      END IF
    7 FORMAT(/' *** ERROR (MOLLAB), error reading unit',I4,
     *       /T22,'when searching for label ',A8,
     *       /T22,'IOSTAT value :',I7)
      END
!C  /* Deck fndlab */
      LOGICAL FUNCTION FNDLAB(A,LU)
!C
!C 26-May-1985 hjaaj -- logical function version of SEARCH
!C 16-Jun-1986 hjaaj -- changed to CHARACTER*8 from REAL*8
!C
      CHARACTER*8 A, B(4), C
      DATA C/'********'/
      IRDERR = 0
    1 READ(LU,END=3,ERR=6)B
      IRDERR = 0
      IF (B(1).NE.C) GO TO 1
      IF (B(4).NE.A) GO TO 1
      FNDLAB = .TRUE.
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      GO TO 8
    3 CONTINUE
#if defined (VAR_MFDS)
!C 880916-hjaaj -- attempt to fix an IBM problem
!C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
!C     of FTxxF001), backspace makes LU ready for append.
!C 940510-hjaaj: same change for Cray's multifile datasets
      BACKSPACE LU
#endif
C
    8 FNDLAB = .FALSE.
C
   10 RETURN
      END
!C  /* Deck fndlab12 */
      LOGICAL FUNCTION FNDLAB12(A,LU)
!C
!C 26-May-1985 hjaaj -- logical function version of SEARCH
!C 16-Jun-1986 hjaaj -- changed to CHARACTER*8 from REAL*8
!C 02-May-2010 sk    -- copy/paste from fndlab + changed to CHARACTER*12 
!C                      from CHARACTER*8
!C
      CHARACTER*12 B(4), C
      CHARACTER*12 A
      DATA C/'************'/
      IRDERR = 0
    1 READ(LU,END=3,ERR=6)B
      IRDERR = 0
      IF (B(1).NE.C) GO TO 1
      IF (B(4).NE.A) GO TO 1
      FNDLAB12 = .TRUE.
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      GO TO 8
    3 CONTINUE
#if defined (VAR_MFDS)
!C 880916-hjaaj -- attempt to fix an IBM problem
!C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
!C     of FTxxF001), backspace makes LU ready for append.
!C 940510-hjaaj: same change for Cray's multifile datasets
      BACKSPACE LU
#endif
C
    8 FNDLAB12 = .FALSE.
C
   10 RETURN
      END
C  /* Deck mollb2 */
      SUBROUTINE MOLLB2(SRCLBL,RTNLBL,LU,LUERR)
C
C  28-Jun-1986 hjaaj
C  (as MOLLAB, but returns two middle labels in RTNLBL(2))
C
C  Purpose:
C     Search for MOLECULE labels on file LU
C
      CHARACTER*8 SRCLBL, RTNLBL(2), B(4), STAR8
      PARAMETER (STAR8 = '********')
C
      IRDERR = 0
    1 READ (LU,END=3,ERR=6,IOSTAT=IOSVAL) B
      IRDERR = 0
      IF (B(1).NE.STAR8)  GO TO 1
      IF (B(4).NE.SRCLBL) GO TO 1
C
         RTNLBL(1) = B(2)
         RTNLBL(2) = B(3)
         IF (LUERR.LT.0) LUERR = 0
         RETURN
C
    3 IF (LUERR.LT.0) THEN
#if defined (VAR_MFDS)
!C 880916-hjaaj -- attempt to fix an IBM problem
!C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
!C     of FTxxF001), backspace makes LU ready for append.
!C 940510-hjaaj: same change for Cray's multifile datasets
         BACKSPACE LU
#endif
         LUERR = -1
         RETURN
      ELSE
         WRITE (LUERR,4) SRCLBL,LU
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR (MOLLB2) MOLECULE label not found on file')
      END IF
    4 FORMAT(/' *** ERROR (MOLLB2), MOLECULE label ',A8,
     *        ' not found on unit',I4)
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      IF (LUERR.LT.0) THEN
         LUERR = -2
         RETURN
      ELSE
         WRITE (LUERR,7) LU,SRCLBL,IOSVAL
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR (MOLLB2) error reading file')
      END IF
    7 FORMAT(/' *** ERROR (MOLLB2), error reading unit',I4,
     *       /T22,'when searching for label ',A8,
     *       /T22,'IOSTAT value :',I7)
      END
C  /* Deck fndlb2 */
      LOGICAL FUNCTION FNDLB2(SRCLBL,RTNLBL,LU)
C
C  5-Aug-1986 hjaaj
C  (as FNDLAB, but returns two middle labels in RTNLBL(2))
C
      CHARACTER*8 SRCLBL, RTNLBL(2), B(4), STAR8
      PARAMETER (STAR8 = '********')
      IRDERR = 0
    1 READ (LU,END=3,ERR=6) B
      IRDERR = 0
      IF (B(1).NE.STAR8)  GO TO 1
      IF (B(4).NE.SRCLBL) GO TO 1
      FNDLB2    = .TRUE.
      RTNLBL(1) = B(2)
      RTNLBL(2) = B(3)
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      GO TO 8
    3 CONTINUE
#if defined (VAR_MFDS)
!C 880916-hjaaj -- attempt to fix an IBM problem
!C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
!C     of FTxxF001), backspace makes LU ready for append.
!C 940510-hjaaj: same change for Cray's multifile datasets
      BACKSPACE LU
#endif
    8 FNDLB2 = .FALSE.
C
   10 RETURN
      END
C  /* Deck nxtlab */
      LOGICAL FUNCTION NXTLAB(SRCLBL, RTNLBL, LU)
C
C  3-Nov-1986 hjaaj
C  (find and return next MOLECULE label on LU,
C   NXTLAB false if no label found)
C
      CHARACTER*8 SRCLBL, RTNLBL(2), B(4), STAR8
      PARAMETER ( STAR8 = '********' )
      IRDERR = 0
    1 READ(LU,END=3,ERR=6) B
      IRDERR = 0
      IF (B(1) .NE. STAR8) GO TO 1
      NXTLAB = .TRUE.
      SRCLBL = B(4)
      RTNLBL(1) = B(2)
      RTNLBL(2) = B(3)
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      GO TO 8
    3 CONTINUE
#if defined (VAR_MFDS)
!C 880916-hjaaj -- attempt to fix an IBM problem
!C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
!C     of FTxxF001), backspace makes LU ready for append.
!C 940510-hjaaj: same change for Cray's multifile datasets
      BACKSPACE LU
#endif
    8 NXTLAB = .FALSE.
C
   10 RETURN
      END
C  /* Deck dmplab */
      SUBROUTINE DMPLAB(LU,LUPRI)
C
C 27-Mar-1987 hjaaj -- dump remaining labels on file LU
C
      CHARACTER*8 B(4), C
      PARAMETER ( C = '********' )
C
      WRITE (LUPRI, '(//A,I5)') ' >>> DUMP OF LABELS ON UNIT',LU
      IRDERR = 0
      IREC = 0
    1 READ (LU,END=3,ERR=6,IOSTAT=IOSVAL) B
         IRDERR = 0
         IREC = IREC + 1
         IF (B(1).EQ.C) THEN
            WRITE (LUPRI, '(A,I6,4(2X,A8))') ' Rec. no.',IREC,B
         END IF
      GO TO 1
C
    6 CONTINUE
         IRDERR = IRDERR + 1
         IREC = IREC + 1
         WRITE (LUPRI, '(/A,I6,A,I7)')
     &      ' >>> ERROR reading rec. no.',IREC,'; IOSTAT value',IOSVAL
         IF (IRDERR .LT. 5) GO TO 1
         WRITE (LUPRI, '(/A)')
     &      ' >>> ERROR exit from DMPLAB: 5 consecutive read errors <<<'
    3 CONTINUE
      REWIND (LU)
C
      WRITE (LUPRI, '(/I10,A)') IREC,' records read from file.'
      RETURN
      END
C  /* Deck newlab */
      SUBROUTINE NEWLAB(LABEL,LU,LUERR)
C
C  29-Sep-1988 Hans Joergen Aa. Jensen
C
C  Write new MOLECULE-type label to LU
C
      CHARACTER*8 LABEL, LTIME, LDATE, LSTARS
      DATA LSTARS /'********'/
      character*100 msg
      CALL GETDAT(LDATE,LTIME)

#if !defined VAR_IFORT
!miro: ifort 10.1 does not know the "iomsg" parameter ....
      WRITE (LU,ERR=1000,IOSTAT=IOSVAL,iomsg=msg) 
     &              LSTARS,LDATE,LTIME,LABEL
#else
      WRITE (LU,ERR=1000,IOSTAT=IOSVAL) LSTARS,LDATE,LTIME,LABEL
#endif
      RETURN
C
 1000 IF (LUERR .LT. 0) THEN
         LUERR = -2
         RETURN
      ELSE
         WRITE (LUERR,'(/3A,I5/A,I7)')
     &      ' NEWLAB: error writing label "',LABEL,'" to unit',LU,
     &      '         IOSTAT value',IOSVAL
         WRITE(LUERR,*) 'LSTARS=',LSTARS
         WRITE(LUERR,*) 'LDATE=',LDATE,' LTIME=',LTIME
#if !defined VAR_IFORT
         WRITE (LUERR,'(A,A)') 'message:',msg
#endif
         CALL QTRACE(LUERR)
         CALL QUIT('NEWLAB: output error writing label')
      END IF
      END
C  /* Deck newlb2 */
      SUBROUTINE NEWLB2(LABEL,RTNLBL,LU,LUERR)
C
C  29-Sep-1988 Hans Joergen Aa. Jensen
C
C  Write new MOLECULE-type label to LU
C
      CHARACTER*8 LABEL, RTNLBL(2), LSTARS
      DATA LSTARS /'********'/
      WRITE (LU,ERR=1000,IOSTAT=IOSVAL) LSTARS,RTNLBL,LABEL
      RETURN
C
 1000 IF (LUERR .LT. 0) THEN
         LUERR = -2
         RETURN
      ELSE
         WRITE (LUERR,'(/3A,I5/A,I7)')
     &      ' NEWLAB2: error writing label "',LABEL,'" to unit',LU,
     &      '         IOSTAT value',IOSVAL
         CALL QTRACE(LUERR)
         CALL QUIT('NEWLB2: output error writing label')
      END IF
      END
C  /* Deck newlab12 */
      SUBROUTINE NEWLAB12(LABEL,LU,LUERR)
C
C  29-Sep-1988 Hans Joergen Aa. Jensen
C
C  Write new MOLECULE-type label to LU
C
      CHARACTER*12 LABEL, LTIME, LDATE, LSTARS
      DATA LSTARS /'************'/
      CALL GETDAT(LDATE,LTIME)
      WRITE (LU,ERR=1000,IOSTAT=IOSVAL) LSTARS,LDATE,LTIME,LABEL
      RETURN
C
 1000 IF (LUERR .LT. 0) THEN
         LUERR = -2
         RETURN
      ELSE
         WRITE (LUERR,'(/3A,I5/A,I7)')
     &      ' NEWLAB12: error writing label "',LABEL,'" to unit',LU,
     &      '           IOSTAT value',IOSVAL
         CALL QTRACE(LUERR)
         CALL QUIT('NEWLAB12: output error writing label')
      END IF
      END
C  /* Deck second */
#if !defined (SYS_CRAY) && !defined (VAR_SECSEC)
      FUNCTION SECOND ()
#if defined (SYS_AIX) || defined (SYS_ALLIANT) || defined (SYS_CONVEX) \
      || defined (SYS_HPUX) || defined (SYS_IRIX) || defined (SYS_DEC) \
      || defined (SYS_T3D) || defined (SYS_SUN) || defined (SYS_LINUX) \
      || defined (SYS_DEC) || defined (SYS_SX) || defined (SYS_DARWIN)
#include "implicit.h"
      real*4 etime,tarray(2)
      SECOND = etime(tarray)
#endif
      RETURN
      END
#endif
C  /* Deck sotmat */
      SUBROUTINE SOTMAT(NMO,UMO,IFAIL)
C
C  16-Feb-1986 Hans Jorgen Aa. Jensen
C
C  Purpose:
C
C    Construct the orbital transformation matrix
C    for transforming a CAS CI vector using a sequence
C    of single orbital transformations as described by
C    Per-Ake Malmquist.
C
C    The matrix is
C
C      C   +  C   = (1 - L) + U(inv);
C       L      U
C
C    where L and U constitute the LU decomposition of the
C    orthogonal orbital transformation matrix UMO.
C
#include "implicit.h"
      DIMENSION UMO(NMO,NMO)
#include "priunit.h"
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( THRESH = 1.D-4 )
C
      CALL QENTER('SOTMAT')
C
C     STEP 1: The LU decomposition
C     ============================
C
      DO 220 K = 1,NMO
         X = UMO(K,K)
         IF (ABS(X) .LE. THRESH) GO TO 960
         X = D1 / X
         DO 120 I = K+1,NMO
            UMO(I,K) = UMO(I,K) * X
  120    CONTINUE
         DO 200 J = (K+1),NMO
            Y = UMO(K,J)
            DO 180 I = K+1,NMO
               UMO(I,J) = UMO(I,J) - UMO(I,K)*Y
  180       CONTINUE
  200    CONTINUE
  220 CONTINUE
C
C     STEP 2: U inverse
C     =================
C
      DET = D1
      DO 300 K = 1,NMO
         DET = DET * UMO(K,K)
  300 CONTINUE
#if defined (VAR_SOTMATTEST)
C
C
      WRITE (LUPRI,'(//A,1P,D10.2)') ' SOTMAT: UMO determinant =',DET
      WRITE (LUPRI,'(//A)') ' SOTMAT: LU matrix'
      CALL OUTPUT (UMO,1,NMO,1,NMO,NMO,NMO,1,6)
#endif
      IF (ABS(DET) .LE. THRESH) GO TO 980
C
      DO 400 K = 1,NMO
         UIDIAG = D1 / UMO(K,K)
         UMO(K,K) = UIDIAG
         DO 380 J = 1,(K-1)
            SUM = D0
            DO 360 I = J,(K-1)
               SUM = SUM - UMO(J,I)*UMO(I,K)
  360       CONTINUE
            UMO(J,K) = UIDIAG * SUM
  380    CONTINUE
  400 CONTINUE
C
C     STEP 3: construct 1 - L
C     =======================
C
      DO 600 K = 1,NMO
         DO 580 J = (K+1),NMO
            UMO(J,K) = -UMO(J,K)
  580    CONTINUE
  600 CONTINUE
C
C     NORMAL AND ERROR RETURNS
C     ========================
C
      IFAIL = 0
      GO TO 9999
C
  960 CONTINUE
      WRITE (LUPRI,'(///A,1P,D10.2/)')
     1   ' >>> SOTMAT: DIAGONAL ELEMENT TOO SMALL:',X
      IFAIL = 1
      GO TO 9999
C
  980 CONTINUE
      WRITE (LUPRI,'(///2A,1P,D10.2/)')
     1   ' >>> SOTMAT: U MATRIX TOO CLOSE TO SINGULARITY,',
     2   ' DETERMINANT =',DET
      IFAIL = 2
C
 9999 CALL QEXIT('SOTMAT')
      RETURN
C
C     END OF SOTMAT.
C
      END
C  /* Deck nofdia */
      INTEGER FUNCTION NOFDIA(N,NDIM,AMAT,THREQL)
C
C  2-Jul-1992 hjaaj
C  This function returns the number of off-diagonal elements
C  with absolute value greater than THREQL.
C
#include "implicit.h"
      DIMENSION AMAT(NDIM,NDIM)
C
      CALL QENTER('NOFDIA')
      NELEM = 0
      DO 220 K = 1,N
         DO 210 I = 1,N
            IF (I .NE. K .AND. ABS(AMAT(I,K)) .GT. THREQL) THEN
               NELEM = NELEM + 1
            END IF
  210    CONTINUE
  220 CONTINUE
      NOFDIA = NELEM
      CALL QEXIT('NOFDIA')
      RETURN
      END
C  /* Deck fndmin */
      SUBROUTINE FNDMIN(NELMNT,IPLACE,VEC,NVEC,WRK,LWRK)
C
C 23-Nov-2000 hjaaj (FNDMIN = CIFNMN from 12-Aug-1990 hjaaj)
C (Find minimum elemnts)
C
C Purpose:
C   Return in IPLACE addresses on lowest NELMNT elements in VEC.
C
#include "implicit.h"
      DIMENSION VEC(NVEC), IPLACE(NELMNT), WRK(LWRK)
C
      IF (LWRK .LT. NELMNT) THEN
         CALL QUIT('FNDMIN: Insufficient memory (LWRK .lt. NELMNT)')
      END IF
      IF (NELMNT .GT. NVEC) THEN
         CALL QUIT('FNDMIN ERROR: NELMNT .gt. NVEC')
      END IF
C
C     Sort first NELMNT elements of VEC
C
      DO 120 I = 1,NELMNT
         VECI = VEC(I)
         DO 115 J = 1,(I-1)
         IF (WRK(J) .GT. VECI) THEN
            DO 112 K = I,(J+1),-1
               WRK(K)    = WRK(K-1)
               IPLACE(K) = IPLACE(K-1)
  112       CONTINUE
            IPLACE(J) = I
            WRK(J)    = VECI
         GO TO 120
         END IF
  115    CONTINUE
         IPLACE(I) = I
         WRK(I)    = VECI
  120 CONTINUE
C
C     Find lowest elements by insertion sort
C
      XHGH = WRK(NELMNT)
      DO 140 I = NELMNT+1,NVEC
      IF (VEC(I).GE.XHGH) GO TO 140
         DO 130 J = 1,NELMNT
         IF (VEC(I) .LT. WRK(J)) THEN
            DO 135 K = NELMNT,(J+1),-1
               WRK(K) = WRK(K-1)
               IPLACE(K) = IPLACE(K-1)
  135       CONTINUE
            IPLACE(J) = I
            WRK(J) = VEC(I)
            XHGH   = WRK(NELMNT)
            GO TO 140
         END IF
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
C  /* Deck rewmot */
      SUBROUTINE REWSPL(LU)
C
C     Short interface routine for rewinding a file that may have been 
C     split, to ensure that we search for labels on the first of the split
C     files. This routines preserves the UNIT number.
C
C     K.Ruud, April 19 (2000)
C
#if !defined (VAR_SPLITFILES)
      REWIND (LU)
#else
#include "implicit.h"
#include "dummy.h"
      CHARACTER FNNAME*80
C
      INQUIRE (UNIT=LU,NAME=FNNAME)
      LN = 1
 14   CONTINUE
      IF (FNNAME(LN:LN) .EQ. '-') THEN
         LN = LN - 1
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         CALL GPOPEN(LU,FNNAME(1:LN),'OLD','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
         INQUIRE (UNIT=LU,NAME=FNNAME)
         GOTO 15
      ELSE IF (FNNAME(LN:LN) .EQ. ' ') THEN
         GOTO 15
      END IF
      LN = LN + 1
      GOTO 14
 15   CONTINUE
      REWIND (LU)
#endif
      RETURN
      END

C  /* Deck gpinq */
      SUBROUTINE GPINQ(FILENM,TASK,VALUE)
C
C     Generalized routine for checking the existence of a file.
C     Needed for parallel I/O calculations, where the name of the file
C     may have been altered with a processor-specific appendage.
C     The routine is only needed when searching on a filename
C
C     K.Ruud, San Diego Aug 2000
C
#include "implicit.h"
      CHARACTER*(*) FILENM, TASK
      CHARACTER FILTMP*80
      INTEGER FILELN
      LOGICAL VALUE
#include "chrnos.h"
! maxorb for infpar.h
!include "maxorb.h"
!include "infpar.h"
C
      FILELN = LEN(FILENM)
!     IF (MYTID .GT. 0) THEN
!        FILTMP(1:FILELN+4) = FILENM(1:FILELN)//
!    &                        '#'//CHRNOS(MYTID/100)//
!    &                        CHRNOS(MOD(MYTID,100)/10)//
!    &                        CHRNOS(MOD(MOD(MYTID,100),10))
!        FILELN = FILELN + 4
!     ELSE
         FILTMP(1:FILELN) = FILENM(1:FILELN)
!     END IF
      IF (TASK(1:5) .EQ. 'EXIST') THEN
         INQUIRE(FILE=FILTMP(1:FILELN),EXIST=VALUE)
      ELSE IF (TASK(1:5) .EQ. 'OPENE') THEN
         INQUIRE(FILE=FILTMP(1:FILELN),OPENED=VALUE)
      ELSE
         CALL QUIT('Unknown inquire status requested')
      END IF
      RETURN
      END
      INTEGER FUNCTION LU_NOT_USED(LU_start)
      ! find a file logical unit number not in use
      ! hjaaj Jan 2017
      implicit none
      INTEGER  LU_start, LU, ios
      LOGICAL  is_open
!
      do LU = LU_start,99
         inquire (unit=LU, opened=is_open, iostat=ios)
         if (ios.ne.0) cycle
         if (.not.is_open) then
            LU_NOT_USED = LU
            RETURN
         end if
      end do
!
      CALL QUIT('LU_NOT_USED: no free i/o unit numbers')
      END
! --- end of gp/gphjj.F ---
