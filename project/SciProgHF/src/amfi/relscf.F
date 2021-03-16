C
C Authorship of this program should be dicussed before
C distributing it within the DIRAC package (MI,jan 2005);
C see email from Timo (Atomares SCF).
C
C NOTE: The previous  value of NDIMPQ was 500000. 
C  Replaced by a parameter in scfarr.h, NDIMPQPAR.
C This increase is needed for larger atoms (Rn for example).
C
CMI/June 08 ... array allocations are from DIRAC (WORK-style), no longer
C              big static allocations.
C
C

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADDMA(N,S,OVE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
      DIMENSION S(N),OVE(N)
      DO I=1,N
        OVE(I)=OVE(I)+S(I)
      END DO
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE AMFIALARM
     *(P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C According to  the entering NALARM print-out various convergence
C statuses (nonconvergence, divergence etc...) 
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "relscf_dim.h"
#include "amfi_if.h"
C     DIMENSION P(2),Q(2),S(2),U(2),T(2),DT(2),DOS(2),PCAP(2),
      DIMENSION P(*),Q(*),S(*),U(*),T(*),DT(*),DOS(*),PCAP(*),
     *QCAP(*),FC(*),FO(*),CM1SCF(MXBF,1),CM2SCF(MXBF,1),
     *CSAV1(MXBF,1),CSAV2(MXBF,1),GMIN(*),SMIN(*),THRE(*),EPS(2),
     *EC(2),CM1(2),Z1(2),Z2(2),DELZ(2),
CMI  *C(MXBF,1),CTRAN(2),PQN(*),ZETA(*),ZETINC(*),
     *C(MXBF,*),CTRAN(*),PQN(*),ZETA(*),ZETINC(*),
     *ITRAN(2),NROW(*),NBVAR(*),IORCN(*)
C
      NALAR1=NALARM-1

      IF (IPR_RELSCF.GE.5) THEN
        WRITE(LUPRI,*) 'ALARM:  NALARM=',NALARM
      ENDIF

      GO TO (20,35,50,65,75,85,95),NALAR1

   20 WRITE(LUPRI,25)
   25 FORMAT(/,5X,'DIAGONALIZATION DIVERGENCE')
      WRITE (LUPRI,30) THRE(N0),N0,DGATH
   30 FORMAT(6X, 'THRESHOLD',E9.2,1X,'FOR VECTOR',I2,1X,'GREATER THAN DG
     1ATH', E9.2)
      CALL FLSHFO(LUPRI)
      CALL QUIT('relscf error exit - DIAGONALIZATION DIVERGENCE')
      GO TO 105


   35 WRITE(LUPRI,40)
   40 FORMAT(/,5X,'SELF CONSISTENCE DIVERGENCE')
      WRITE(LUPRI,45) THRSCF,SCFAT
   45 FORMAT(6X,'SCF THRESHOLD',E9.2,1X,'GREATER THAN SCFAT', E9.2)
      CALL FLSHFO(LUPRI)
      CALL QUIT('relscf error exit - SELF CONSISTENCE DIVERGENCE')
      GO TO 105

CMI .... too many iterations ...
   50 WRITE(LUPRI,55)
   55 FORMAT(/,5X,'TOO MANY SCF ITERATIONS!')
      WRITE(LUPRI,60) NITSCF,MXXTRP
   60 FORMAT(6X,'ITERATIONS NUMBER=',I2,1X,'EQUAL MXXTRP=',I2)
      WRITE(LUPRI,'(2X,A)')
     &    'Change relscf atom occupation, or the basis set !'
CMI  ... rather error exit !
      CALL FLSHFO(LUPRI)
      CALL QUIT('relscf error exit - TOO MANY SCF ITERATIONS! !')
      GO TO 105

   65 WRITE(LUPRI,70) RED,RMNED
   70 FORMAT(/,6X,'RELATIVE ENERGY DIFFERENCE',E9.2,1X,
     *'SMALLER THAN RMNED', E9.2)
      CALL FLSHFO(LUPRI)
      CALL QUIT('relscf error exit on RELATIVE ENERGY DIFFERENCE')
      GO TO 105

   75 WRITE(LUPRI,80) RED,RMXED
   80 FORMAT(6X,'RELATIVE ENERGY DIFFERENCE',E9.2,1X,
     *'GREATER THAN RMXED',E9.2)
      CALL FLSHFO(LUPRI)
      CALL QUIT('relscf error exit on RELATIVE ENERGY DIFFERENCE 2')
      GO TO 105

   85 CONTINUE
      WRITE (LUPRI,90) ZETA(NBVAR1),NBVAR1,EXPMN
   90 FORMAT(/,6X,'EXPONENT ZETA',F9.5,1X,'FOR ORBITAL',I2,1X,
     &  'SMALLER THAN EXPMN', F9.5)
      CALL FLSHFO(LUPRI)
      CALL QUIT('relscf error exit on EXPONENT ZETA')
      GO TO 105

   95 CONTINUE
      WRITE (LUPRI,100) NBVAR1,IIII,ZETDIF,DMNXP
  100 FORMAT(/,6X,'EXPONENT REDUNDANCY FOR ORBITALS',I2,1X,'AND',
     *I2,2X,'ZEDIF',F9.5,1X,'SMALLER THAN DMNXP',F9.5)
      CALL FLSHFO(LUPRI)
      CALL QUIT('relscf error exit on EXPONENT REDUNDANCY FOR ORBITALS')

  105 RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   
      SUBROUTINE AT3(WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Internal Atomic SCF Code producing scalar relativistic atomic
C   orbitals for the AMFI mean-field summation.
C
C   Called from INT_AT_SCF (from AMFI). It is using the DIRAC memory manager throughout.
C
C   CHANGE FOLLOWING TWO VALUES TO INCREASE DIMENSIONS OF PROGRAM:
C
C     MXBAS IS MAXIMAL NUMBER OF BASIS FUNCTIONS
C     MXOCC IS MAXIMAL SUM OF ALL OCCUPATION NUMBERS (OPEN AND CLOSED)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "scfarr.h"
#include "relscf_v10.h"
      DIMENSION WORK(LWORK)
#include "amfi_if.h"
#include "memint.h"

      MXBAS=MXB
      MXOCC=MXO

      IF (IPR_RELSCF.GE.2) THEN
        WRITE (LUPRI,991) MXBAS,MXOCC
991     FORMAT(/' MAXIMAL NUMBER OF BASIS FUNCTIONS   =',I10/,
     *          ' MAXIMAL NUMBER OF OCCUPIED ORBITALS =',I10/)
      ENDIF
C
C     10000 CALCULATIONS IN ONE CHUNK SHOULD BE ENOUGH ...
C
      CALL MEMGET('REAL',KC     ,LDIMC+2  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCTRAN ,LDIMC  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPQN   ,MXB    ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KZETA  ,MXB    ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KZETINC,MXB    ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KITRAN ,LDIMC  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNROW  ,MXB    ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNBVAR ,MXB    ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIORCN ,MXB    ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCREAD ,MXB*MXO,WORK,KFREE,LFREE)
      CALL INPUT_AT34(WORK(KC),WORK(KCTRAN),WORK(KPQN),WORK(KZETA),
     &           WORK(KZETINC),WORK(KITRAN),WORK(KNROW),WORK(KNBVAR),
     &           WORK(KIORCN),WORK(KCREAD),WORK(KFREE),LFREE)
      CALL MEMREL('AT3',WORK,KWORK,KWORK,KFREE,LFREE)

      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE AT34R(N,ISIZE,CHARGE,
     *SMAT,V,H,EV2,BU,P,G,EIG,SINV,REVT,
     *AUX,OVE,EW,E,AA,RR,TT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     AT34R 1.0 - RELATIVISTIV ROUTINES FOR AT34 - 14.5.87
C     BERND HESS
C
C     INPUT: SMAT    OVERLAP MATRIX
C            V     POTENTIAL
C            H     RELATIVISTIC KINETIC ENERGY
C            EV2   PVP INTEGRALS
C
C
C Array dimensions fixes:
C--------------------------
C  MULT(ISIZE) -> MULT(*), MULT part of REL, which is from MEMGET
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      DIMENSION V(ISIZE),SMAT(ISIZE),P(ISIZE),G(ISIZE),
     *          H(ISIZE),BU(ISIZE),EV2(ISIZE),
     *          EIG(N,N),SINV(N,N),REVT(N,N),AUX(N,N),OVE(N,N),
     *          EW(N),E(N),AA(N),RR(N),TT(N)

      VELIT = 137.036D0
      TOL   = 1.0D-14
      CON   = VELIT*VELIT
      PREA  = 1.0D0/CON
      CON2  = PREA+PREA
C
C     SCHMIDT-ORTHOGONALIZE OVERLAP MATRIX
C
      CALL SOG(N,SMAT,SINV,P,OVE,EW)
      CALL FILLMA(N,SMAT,OVE)
C
C ** SINV CONTAINS TRANSFORMATION TO ORTHOGONAL AO-BASIS
C ** OVE  CONTAINS OVERLAP MATRIX IN FULL
C*******************************************************************
C     RELATIVISTIC OPERATORS
C*******************************************************************
C-------------------------------------------------------------------
C     MATRIX REPRESENTATION CALCULATED FROM DIAGONALIZED SQRA MATRIX
C-------------------------------------------------------------------
C.      DO 8 I=1,ISIZE
C.  8   P(I)=H(I)
C.      CALL DIAG(P,N,EIG,EW,SINV,AUX,0)
C.      DO 347 I=1,N
C.      TT(I)=0.5D0*( ((EW(I)+CON)**2-CON*CON)/CON)
C.  347 E(I)=EW(I)+CON
C--------------------------------------------------------------------
C     MATRIX REPRESENTATION CALCULATED FROM NONRELATIVISTIC T MATRIX
C--------------------------------------------------------------------
      CALL DIAGX(H,N,EIG,EW,SINV,AUX,0)
      CALL DCOPY(N,EW,1,TT,1)
!DIR$ NOVECTOR
      DO  I=1,N
C
C     IF T SUFFICIENTLY SMALL, USE SERIES EXPANSION TO AVOID CANCELLATION
C
        RATIO=TT(I)/VELIT
        IF (RATIO.GT.0.02D0) THEN
          EW(I)=CON*(DSQRT(1.D0+CON2*TT(I))-1.D0)
        ELSE
          TV1= TT(I)
          TV2= -TV1*TT(I)*PREA/2.D0
          TV3= -TV2*TT(I)*PREA
          TV4= -TV3*TT(I)*PREA*1.25D0
          EW(I) = TV1+TV2+TV3+TV4
        ENDIF
        E(I)=EW(I)+CON
      END DO
C---------------------------------------------------------------------
C     CALCULATE REVERSE TRANSFORMATION
C---------------------------------------------------------------------
C
C     CALCULATE TRANSFORMATION MATRICES
C
      DO I=1,N
        DO J=1,N
          AUX(I,J)=0.D0
          DO K=I,N
            AUX(I,J)=AUX(I,J)+SINV(I,K)*EIG(K,J)
          END DO
        END DO
      END DO
      DO I=1,N
        DO J=1,N
          REVT(I,J)=0.D0
          DO K=1,N
            REVT(I,J)=REVT(I,J)+OVE(I,K)*AUX(K,J)
          END DO
        END DO
      END DO
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          H(IJ)=0.D0
          DO  K=1,N
            H(IJ)=H(IJ)+REVT(I,K)*REVT(J,K)*EW(K)
          END DO
        END DO
      END DO
      DO I=1,N
        AA(I)=DSQRT((CON+E(I)) / (2.D0*E(I)))
        RR(I)=DSQRT(CON)/(CON+E(I))
      END DO
C
C     BEYOND THIS POINT, OVE IS USED AS SCRATCH ARRAY
C
C
C    TRANSFORM V TO T-BASIS
C
      CALL TRSM(V,SINV,G,N,AUX,OVE)
      CALL TRSM(G,EIG,BU,N,AUX,OVE)
C
C    MULTIPLY
C
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          P(IJ)=-BU(IJ)*CHARGE
          BU(IJ)= P(IJ)*AA(I)*AA(J)
        END DO
      END DO

      CALL TRSMT(BU,REVT,V,N,AUX,OVE)
C
C     PVP INTEGRALS
C
      CALL TRSM(EV2,SINV,G,N,AUX,OVE)
      CALL TRSM(G,EIG,BU,N,AUX,OVE)
C
C    MULTIPLY
C
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          G(IJ)=-BU(IJ)*CHARGE
          BU(IJ)= G(IJ)*AA(I)*RR(I)*AA(J)*RR(J)
        END DO
      END DO

      CALL TRSMT(BU,REVT,EV2,N,AUX,OVE)
      CALL ADDMA(ISIZE,EV2,V)
C
C     CALCULATE EVEN2 OPERATOR
C
      CALL EVEN2(N,P,G,E,AA,RR,TT,EIG,AUX,OVE)
C
C    TRANSFORM BACK
C
      CALL TRSMT(G,REVT,EV2,N,AUX,OVE)
      CR=1.0D0/CHARGE
      FAC_CRX = -1.0D0*CR
      CALL ADDMA(ISIZE,EV2,V)
      CALL DSCAL(ISIZE,FAC_CRX,V,1)

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE AT3DIM(NCOVAR,C,CTRAN,PQN,ZETA,ZETINC,ITRAN,
     &                  NROW,NBVAR,IORCN,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Contains many memory allocations (MEMGET) for the HFMAIN subroutine.
C
C   Called from INPUT_AT34.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "scfarr.h"
#include "relscf_ijpair.h"
#include "relscf_dim.h"
#include "relscf_v10.h"
C
Cvec  VERSION INCLUDING IDK ARRAY
C
#include "amfi_if.h"
c     logical relscf_verbose
      DIMENSION WORK(LWORK)
CMI   DIMENSION C(LDIMC),CTRAN(LDIMC),PQN(MXB),ZETA(MXB),ZETINC(MXB),
      DIMENSION C(*),CTRAN(LDIMC),PQN(MXB),ZETA(MXB),ZETINC(MXB),
     &          ITRAN(LDIMC),NROW(MXB),NBVAR(MXB),IORCN(MXB)
#include "memint.h"
C
C
C     THE SEQUENCE OF THE ARRAYS PASSED TO VARIOUS SUBROUTINES ARE
C     P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,CSAV1,CSAV2,
C     GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
C     C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC
C
      IDK(1)=0
      DO I=2,MXBF
        IDK(I)=IDK(I-1)+I-1
      END DO
      I10=1

CMI         ... P variable
      CALL MEMGET('REAL', K1,LQM  ,WORK,KFREE,LFREE)
CMI   CALL MEMGET('REAL', K1,LQM+1,WORK,KFREE,LFREE)

CMI         ... Q variable 
C     CALL MEMGET('REAL', N2,LQM  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL', N2,LQM*2,WORK,KFREE,LFREE)

CMI         ... S(LPQ) variable 
      CALL MEMGET('REAL', N3,LPQ  ,WORK,KFREE,LFREE)
CMI   CALL MEMGET('REAL', N3,LPQ+1,WORK,KFREE,LFREE)

CMI    ... U(LPQ) variable
      CALL MEMGET('REAL', N4,LPQ  ,WORK,KFREE,LFREE)
C     CALL MEMGET('REAL', N4,LPQ+3,WORK,KFREE,LFREE)


CMI    ... T(LPQ) variable
      CALL MEMGET('REAL', N5,LPQ  ,WORK,KFREE,LFREE)

CMI            DT(LPQ,2) 
      CALL MEMGET('REAL', N6,LPQ*2,WORK,KFREE,LFREE)

      CALL MEMGET('REAL', N7,LPQ*2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL', N8,LPQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL', N9,LPQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N10,LPQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N11,LPQ  ,WORK,KFREE,LFREE)
      MX=MXBF*NOR
      CALL MEMGET('REAL',N12,MX   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N13,MX   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N14,MX   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N15,MX   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N16,MX   ,WORK,KFREE,LFREE)
      MM=MX
      IF(MX.LT.LPQ*2) MM=LPQ*2
      CALL MEMGET('REAL',N17,MM   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N18,NOR  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N19,NOR  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N20,MXBF ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N21,MXBF ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N22,MXBF ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N23,MXBF ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N24,MXBF ,WORK,KFREE,LFREE)
C
C     ADDITIONAL ARRAYS FOR RELATIVISTIC CALCULATION
C
      MXLPQ=MX100*(MX100+1)/2
      NREL=7*MXLPQ+5*MX100*MX100+5*MX100
      IF (NREL.LT.4*LPQ) NREL=4*LPQ

      CALL MEMGET('REAL',N25,NREL   ,WORK,KFREE,LFREE)

      CALL MEMGET('REAL',N26,MX100  ,WORK,KFREE,LFREE)

!           ...  S100
      CALL MEMGET('REAL',N27,MXLPQ  ,WORK,KFREE,LFREE)

      CALL MEMGET('REAL',N28,MXLPQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N29,MXLPQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',N30,MX100*4,WORK,KFREE,LFREE)

C//SKIPCc1
      IF (IPR_RELSCF.GE.2) THEN 
        WRITE (LUPRI,'(/,I10,A,/)') KFREE,
     &  '(KFREE) DOUBLE WORDS ARE REQUIRED (COMMING IN WORK ARRAY)'
     & //'TO RUN THIS CASE'
        CALL FLSHFO(LUPRI)
      ENDIF

      IF (IPR_RELSCF.GE.5) THEN
      write(LUPRI,*) 'size(WORK(:K1))=',size(WORK(:K1))
      write(LUPRI,*) 'size(WORK(:N2))=',size(WORK(:N2))
      write(LUPRI,*) 'size(WORK(:N3))=',size(WORK(:N3))
      write(LUPRI,*) 'size(WORK(:N4))=',size(WORK(:N4))
      write(LUPRI,*) 'size(WORK(:N5))=',size(WORK(:N5))
      ENDIF

!      IF(NCOVAR.NE.0) GO TO 11

      CALL HFMAIN(WORK(K1),WORK(N2),WORK(N3),WORK(N4),WORK(N5),
     &            WORK(N6),WORK(N7),WORK(N8),WORK(N9),WORK(N10),
     &            WORK(N11),WORK(N12),WORK(N13),WORK(N14),WORK(N15),
     &            WORK(N16),WORK(N17),WORK(N18),WORK(N19),WORK(N20),
     &            WORK(N21),WORK(N22),WORK(N23),WORK(N24),WORK(N25),
     &            WORK(N26),WORK(N27),WORK(N28),WORK(N29),WORK(N30),
     &            C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN,
     &            WORK(KFREE),LFREE)

      IF (IPR_RELSCF.GE.5) THEN
      write(LUPRI,*) 'after size(WORK(:K1))=',size(WORK(:K1))
      write(LUPRI,*) 'after size(WORK(:N2))=',size(WORK(:N2))
      write(LUPRI,*) 'after size(WORK(:N3))=',size(WORK(:N3))
      write(LUPRI,*) 'after size(WORK(:N4))=',size(WORK(:N4))
      write(LUPRI,*) 'after size(WORK(:N5))=',size(WORK(:N5))
      ENDIF

      CALL MEMREL('AT34DIM',WORK,KWORK,KWORK,KFREE,LFREE)

      RETURN
      END

      SUBROUTINE BESSKA(AVAL,XVAL,KA,KA1)
#include "implicit.h"
#include "priunit.h"
      REAL*8 AVAL,XVAL,EPS,KA,KA1
C
C     SUBROUTINE TO CALCULATE MODIFIED BESSEL FUNCTION OF THE
C     THIRD KIND. FORTRAN SUBPROGRAM ADAPTED FROM ALGOL PROGRAM
C     IN: N.M.TEMME, J.COMP.PHYS 19 (1975) 324
C
C     MODIFICATIONS: RESULT IS DEXP(X) TIMES BESSEL FUNCTION
C                    INSTEAD OF BESSEL FUNCTION ITSELF.
C     RELATIVE ACCURACY FIXED TO 5*1.D-14
C
      REAL*8 A1,B,C,D,E,F,G,H,P,PI,Q,S,A,X,RCPG,SINHX
      INTEGER N,NA
      LOGICAL REC,REV
      A=AVAL
      X=XVAL
      EPS=5.D-14
      PI=4.0D0*DATAN(1.0D0)
      REV=A.LT.-5D0
      IF (REV) A=-A-1.D0
      REC=A.GE.0.5D0
      IF (.NOT. REC) GOTO 1
      NA=A+0.5D0
      A=A-DFLOAT(NA)
1     CONTINUE
      IF (A.NE.-.5D0) GOTO 2
      F=DSQRT(PI/X*.5D0)
      G=F
      GOTO 3
2     IF (X.GE.1.0D0) GOTO 4
      B=X*.5D0
      D=-DLOG(B)
      E=A*D
      C=A*PI
      IF (DABS(C).GE.1.D-15) GOTO 5
      C=1.D0
      GOTO 6
5     C=C/DSIN(C)
6     CONTINUE
      IF (DABS(E).GE.1.D-15) GOTO 7
      S=1.D0
      GOTO 8
7     S=SINHX(E)/E
8     CONTINUE
      E=DEXP(E)
      A1=(E+1.0D0/E)*.5D0
      G=RCPG(A,P,Q)*E
      KA=C*(P*A1+Q*S*D)
      F=KA
      E=A*A
      P=0.5D0*G*C
      Q=0.5D0/G
      C=1.D0
      D=B*B
      KA1=P
      N=1
      GOTO 10
9     N=N+1
10    F=(F*DFLOAT(N)+P+Q)/(DFLOAT(N*N)-E)
      C=C*D/DFLOAT(N)
      P=P/(DFLOAT(N)-A)
      Q=Q/(DFLOAT(N)+A)
      G=C*(P-DFLOAT(N)*F)
      H=C*F
      KA=KA+H
      KA1=KA1+G
      IF (H/KA+DABS(G)/KA.GT.EPS) GOTO 9
      F=KA*DEXP(X)
      G=KA1*DEXP(X)/B
      GOTO 3
4     C=0.25D0-A*A
      G=1.D0
      F=0.D0
      E=X*DCOS(A*PI)/PI/EPS
      N=1
      GOTO 13
12    N=N+1
13    H=(2.D0*(DFLOAT(N)+X)*G
     &   -(DFLOAT(N)-1.D0+C/DFLOAT(N))*F)/(DFLOAT(N)+1)
      F=G
      G=H
      IF (H*DFLOAT(N).LT.E) GOTO 12
      P=F/G
      Q=P
      B=X+X
      E=B-2.D0
15    P=(DFLOAT(N)-1.D0+C/DFLOAT(N)) / (E+(DFLOAT(N)+1.D0)*(2.D0-P))
      Q=P*(Q+1.D0)
      N=N-1
      IF (N.GT.0) GOTO 15
      F=DSQRT(PI/B)/(1.D0+Q)
      G=F*(A+X+0.5D0-P)/X
3     CONTINUE
      IF (.NOT.REC) GOTO 16
      X=2.D0/X
      N=1
17    H=F+(A+DFLOAT(N))*X*G
      F=G
      G=H
      N=N+1
      IF (N.LE.NA) GOTO 17
16    IF (.NOT. REV) GOTO 18
      KA1=F
      KA=G
      GOTO 19
18    KA=F
      KA1=G
19    RETURN
      END

      SUBROUTINE CHECK(ICORE,K,TEXT)
#include "implicit.h"
#include "priunit.h"
#include "amfi_if.h"
c     logical relscf_verbose
      CHARACTER*8 TEXT
      DATA MSK /0/, KSK /5000000/
C
C     GENUEGEND KERNSPEICHER VORHANDEN?
C
c     relscf_verbose = .false.
C
      IF (ICORE.EQ.0) GOTO 4

      L=IABS(ICORE)-K
      IF (KSK.GT.L) KSK=L
      IF (L .GE. 0) GOTO 1
      L=-L
      L=(L+1023)/1024

      WRITE (LUPRI,100) TEXT,K,ICORE,L
100   FORMAT(' ***'/,' *** NEED MORE CORE AT ',A8/,
     *       ' *** NEED ',I10,' IBM+ WORDS'/,
     *       ' *** HAVE ',I10,' IBM+ WORDS'/,
     *       ' *** NEED ',I10,' K WORDS MORE'/,
     *       ' ***')
      CALL FLSHFO(LUPRI)
      CALL QUIT('relscf/CHECK: error 99....')
C
C     PRINTOUT,ENDE ODER RETURN?
C
1     IC=IABS(ICORE)
      IF (MSK.LT.IC) MSK=IC
      IF (ICORE.GT.0) GOTO 3
      L=MSK-KSK
      IC=(L+1023)/1024

CMI   if (relscf_verbose) WRITE (LUPRI,101) L,IC
      IF (IPR_RELSCF.GE.2) THEN
          WRITE (LUPRI,101) L,IC
      ENDIF
101   FORMAT(' DYNAMIC CORE USED SO FAR ',I10,' WORDS (',I10,
     *'K)')
3     RETURN
C
C     ENDE DER RECHNUNG
C
4     IC=(MSK+1023)/1024
      KSK=MSK-KSK
      L=(KSK+1023)/1024
CMI   if (relscf_verbose) WRITE (LUPRI,102) MSK,IC,KSK,L
      if (IPR_RELSCF.GE.2)  THEN
          WRITE (LUPRI,102) MSK,IC,KSK,L
102   FORMAT(' DYNAMIC STORAGE: ',I8,' WORDS (',I5,'K ) - USED ',
     *I8,' WORDS (',I5,'K)'//,
     *' E N D   O F   C A L C U L A T I O N')
      ENDIF

CMI   STOP
      CALL QUIT('CHECK: error output, ICORE=0')
      END

      SUBROUTINE CPLAB (A,B,L,M,N,IA,IB,C,IC,IER)
#include "implicit.h"
#include "priunit.h"
      DIMENSION   A(IA,M),B(IB,N),C(IC,N)
      IF (IA .GE. L .AND. IB .GE. M .AND. IC .GE. L) GO TO 5
      IER=129
      GO TO 9000
    5 IER = 0
      DO I=1,L
           DO J=1,N
                TEMP=0.0D0
                DO K=1,M
                     TEMP=A(I,K)*B(K,J)+TEMP
                END DO
                C(I,J)=C(I,J)+TEMP
           END DO
      END DO
      GO TO 9005
 9000 CONTINUE
 9005 RETURN
      END

      DOUBLE PRECISION FUNCTION DCOF(E,LX,IX)
#include "implicit.h"
#include "relscf_crelop.h"
C
C     FUNCTION DCOF CALCULATES COEFFICIENT D NEEDED FOR FT
C     DCOF(E,LX,IX)= ((-E)**IX * FAK(LX))/(FAK(LX-2*IX)*FAK(IX))
C     V 1.0 - 12.3.86 - BERND HESS
C
      I=LX-2*IX
      D=1.D0
      IF (IX.EQ.0) GOTO 2
      DO K=1,IX
        D=-D*E
      END DO
2     CONTINUE
      D=D*FAK(LX+1)/(FAK(I+1)*FAK(IX+1))
      DCOF=D
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENSI(
     *    P,Q,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C
C     ONE CONFIGURATION
C
CMI - removed S from param.list
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "relscf_dim.h"
      DIMENSION P(*),Q(*),U(*),T(*),DT(*),DOS(*),PCAP(*),
     *QCAP(*),FC(*),FO(*),CM1SCF(MXBF,*),CM2SCF(MXBF,*),
     *CSAV1(MXBF,*),CSAV2(MXBF,*),GMIN(*),SMIN(*),THRE(*),EPS(*),
     *EC(*),CM1(*),Z1(*),Z2(*),DELZ(*),
     *C(MXBF,*),CTRAN(*),PQN(*),ZETA(*),ZETINC(*),
     *ITRAN(*),NROW(*),NBVAR(*),IORCN(*)
C
C     COMPUTE THE DENSITY MATRICES
C
      OCCUCS(1)= 2.D0 ! s-orbitals
      OCCUCS(2)= 6.D0 ! p-orbitals
      OCCUCS(3)=10.D0 ! d-orbitals
      OCCUCS(4)=14.D0 ! f-orbitals
C
      K=1
      NSTEP=0
      J1=0
      DO I=1,NSYM
        IS=0
        IK=0
        NBAS1 = NBAS(I) ! # basis functions
        NOSHIC= NOSH(I) ! # open shells
        NCSHIC= NCSH(I) ! # closed shell
        NOCSH=NOSHIC+NCSHIC
CMI    ... this is second index for C(MXBF,1) !
        J1=J1+NOCSH
        IF (NOSHIC.GT.0) IS=1 ! open shell
        IF (NCSHIC.GT.0) IK=1 ! closed shell
        DO M=1,NBAS1
          DO N=1,M
            TERM1=OCCUP(I)*C(M,J1)*C(N,J1) ! open shells
            TERM2=0.D0
            DO J=1,NCSHIC ! closed shells
              J2=J+NSTEP
              TERM2=TERM2+C(M,J2)*C(N,J2)
            END DO
            TERM2=TERM2*OCCUCS(I)
            IF (M.EQ.N) GO TO 35
            TERM1=2.D0*TERM1
            TERM2=2.D0*TERM2
   35       DOSA=TERM1*IS
            DT(K)=TERM2*IK+DOSA
            DOS(K)=DOSA
            K=K+1
          END DO
        END DO
        NSTEP=J1
      END DO
      RETURN
      END
C
C    CALCULATE INTEGRAL OVER DERIVATIVE OF THE FUNCTIONS
C
      DOUBLE PRECISION FUNCTION DER(IDER,IS1,IS2,AL,BE)
#include "implicit.h"
#include "relscf_crelop.h"
#include "amfi_if.h"
c     logical relscf_verbose
      INTEGER IS1(3),IS2(3),I1(2,3),I2(2,3)
      REAL*8 F(2),G(2)

c     relscf_verbose = .false.

      DO I=1,3
        DO J=1,2
          I1(J,I)=IS1(I)
          I2(J,I)=IS2(I)
        END DO
      END DO
      I1(1,IDER)=I1(1,IDER)+1
      I1(2,IDER)=I1(2,IDER)-1
      I2(1,IDER)=I2(1,IDER)+1
      I2(2,IDER)=I2(2,IDER)-1
      L1=IS1(IDER)+1
      GOTO (10,11,12,13),L1
101   if (IPR_RELSCF.GE.2) WRITE (LUPRI,100) IDER,IS1,IS2,AL,BE
100   FORMAT(' ILLEGAL ANGULAR MOMENTUM (DER)'/,
     &       ' IDER,IS1,IS2,AL,BE PRINTED'/,1X,7I5,3X,2D20.8)
CMI   STOP
      CALL QUIT('relscf error output !')
10    F(1)=-2.D0*AL
      J1=1
      GOTO 19
C
11    F(2)=1.D0
      F(1)=-2.D0*AL
      J1=2
      GOTO 19
C
12    F(2)=2.D0
      F(1)=-2.D0*AL
      J1=2
      GOTO 19
C
13    F(2)=3.D0
      F(1)=-2.D0*AL
      J1=2
C
19    L2=IS2(IDER)+1
      GOTO (20,21,22,23),L2
      GOTO 101
C
20    G(1)=-2.D0*BE
      J2=1
      GOTO 29
C
21    G(2)=1.D0
      G(1)=-2.D0*BE
      J2=2
      GOTO 29
C
22    G(2)=2.D0
      G(1)=-2.D0*BE
      J2=2
      GOTO 29
C
23    G(2)=3.D0
      G(1)=-2.D0*BE
      J2=2
C
29    SUM=0.D0
      DO I=1,J1
        DO J=1,J2
          II=I1(I,1)+I2(J,1)
          JJ=I1(I,2)+I2(J,2)
          KK=I1(I,3)+I2(J,3)
          ANG=THETA(II+JJ,KK)*PHI(JJ,II)
          IF (ANG.EQ.0.D0) GOTO 31
          EX=-DFLOAT(II+JJ+KK+2)*0.5D0
          SUM=SUM+F(I)*G(J)*0.5D0*ANG*GA(II+JJ+KK+2)*
     \            (AL+BE)**EX
31        CONTINUE
        END DO
      END DO
      DER=SUM
      RETURN
      END

      SUBROUTINE DIAGX(A,N,EIG,EW,SINV,AUX,IC)
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(*),AUX(N,N),SINV(N,N),EIG(N,N),EW(N)
      TOL=1.D-15
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUX(I,J)=A(IJ)
          AUX(J,I)=A(IJ)
        END DO
      END DO
      DO K=1,N
        DO J=1,N
          EIG(K,J)=0.D0
          DO L=1,J
            EIG(K,J)=EIG(K,J)+AUX(K,L)*SINV(L,J)
          END DO
        END DO
      END DO
      DO I=1,N
        DO J=1,I
          AUX(I,J)=0.D0
          DO K=1,I
            AUX(I,J)=AUX(I,J)+SINV(K,I)*EIG(K,J)
          END DO
          AUX(J,I)=AUX(I,J)
        END DO
      END DO
      CALL JACOB(AUX,EIG,EW,N,TOL,IC)
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EIGEA (F,
     *P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ONE CONFIGURATION
C     SOLVE FC=ESC
C     NBAS1=DIMENSION OF THE MATRIX F AND VECTOR C
C     N1=NBAS1*(NBAS1+1)/2
C     THE F AND S MATRICES ELEMENTS ARE THE ELEMENTS OF A SUPERVECTOR F
C     OR S,STARTING AT NSTEP+1 FOR THE F AND S SUPERVECTOR, AND OF
C     LENGTH N1
C     NO DEFINE THE COLUMN CORRESPONDING TO THE VECTOR C IN A MATRIX C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "amfi_if.h"
#include "relscf_dim.h"
#include "relscf_v10.h"
CMI   DIMENSION P(2),Q(2),S(2),U(2),T(2),DT(2),DOS(2),PCAP(2),
      DIMENSION P(*),Q(*),S(*),U(*),T(*),DT(*),DOS(2),PCAP(2),
     *QCAP(2),FC(2),FO(2),CM1SCF(MXBF,1),CM2SCF(MXBF,1),
     *CSAV1(MXBF,1),CSAV2(MXBF,1),                THRE(*),EPS(2),
     *EC(2),CM1(2),Z1(2),Z2(2),DELZ(2),GMIN(MXBF,1),SMIN(MXBF,1),
CMI  *C(MXBF,1),CTRAN(2),PQN(*),ZETA(*),ZETINC(*),
     *C(MXBF,*),CTRAN(2),PQN(*),ZETA(*),ZETINC(*),
     *ITRAN(2),NROW(*),NBVAR(*),IORCN(*)
      DIMENSION F(1)
C
C
C     FIRST NORMALIZE THE TRIAL VECTOR
    8 NSTART=1
      GO TO 81
C
    3 NITDIA=1
    1 M=0
      IF (C(1,N0)) 2,6,2
    6 IF (NITDIA-1) 7,7,125
    7 C(1,N0)=0.1D0
      GO TO 8
C     COMPUTE EPS=CT*F*C
    2 EPS(N0)=0.D0
      DO I=1,NBAS1
        DO J=1,I
          K=(I*(I-1)/2)+J
          KP=K+NSTEP
          TERM=C(I,N0)*F(KP)*C(J,N0)
          IF (I.NE.J) TERM=2.D0*TERM
          EPS(N0)=EPS(N0)+TERM
        END DO
      END DO
      IF (M) 96,11,96
C
C     COMPUTE SMIN AND GMIN
   11 DO I=1,NBAS1
        SMIN(I,N0)=0.D0
        GMIN(I,N0)=0.D0
        DO J=1,NBAS1
          IF (J-I) 15,15,20
   15     K=(I*(I-1)/2)+J
          GO TO 25
   20     K=(J*(J-1)/2)+I
   25     KP=K+NSTEP

          if (KP.GT.LPQ) call quit('KP > LPQ !')

          SMIN(I,N0)=SMIN(I,N0)+S(KP)*C(J,N0)
          GMIN(I,N0)=GMIN(I,N0)+(F(KP)-EPS(N0)*S(KP))*C(J,N0)
        END DO  
      END DO
C
C     COMPUTE DOS                                                      0
      DO I=1,NBAS1
        DO J=1,I
          K=(I*(I-1)/2)+J
          KP=K+NSTEP
           DOS(K)=F(KP)-EPS(N0)*S(KP)-SMIN(I,N0)*GMIN(J,N0)-
     &            GMIN(I,N0)*SMIN(J,N0)
        END DO
C
C     MOVE GMIN IN THE FIRST COLUMN OF DOS                             0
        K=(I*(I-1)/2)+1
        DOS(K)=GMIN(I,N0)
      END DO
C
C     COMPUTE T AND SUBSTITUTE IN DOS                                  0
      NBAS11=NBAS1-1
      DO I=1,NBAS11
        IP=NBAS1-I
        K1=IP+1
        DO J=1,IP
          JP=K1-J
          KP=(IP*(IP-1)/2)+JP
          TERM=0.D0
          DO K=K1,NBAS1
            K2=(K-1)*K/2
            KP1=K2+IP
            KP3=(K*(K+1))/2
            KP2=K2+JP
            TERM=TERM+ DOS(KP1)* DOS(KP2)/ DOS(KP3)
          END DO
           DOS(KP)= DOS(KP)-TERM
        END DO
      END DO
C
C     COMPUTE EC=DC
      EC(1)=0.D0
      DO I=2,NBAS1
        I1=I-1
        I2=I*I1/2
        K1=I2+1
        K2=I2+I
        TERM=0.D0
        IF (I1-2) 60,54,54
   54   DO J=2,I1
          K=I2+J
          TERM=TERM+ DOS(K)*EC(J)
        END DO
   60   EC(I)=-( DOS(K1)+TERM)/ DOS(K2)
      END DO
C
C     COMPUTE PROJ=CT*S*DC
      PROJ=0.D0
      DO I=1,NBAS1
        DO J=1,NBAS1
          IF (J-I) 65,65,70
   65     KP=NSTEP+(I*(I-1)/2)+J
          GO TO 75
   70     KP=NSTEP+(J*(J-1)/2)+I
          if (KP.GT.LPQ) call quit('KP > LPQ ! 1')
   75     PROJ=PROJ+C(I,N0)*S(KP)*EC(J)
        END DO
      END DO
C
C     COMPUTE THE NEW VECTOR,NORMALIZE IT AND SAVE THE OLD FOR TEST
      PROJ=1.D0-PROJ
      DO I=1,NBAS1
        CM1(I)=C(I,N0)
        C(I,N0)=PROJ*C(I,N0)+EC(I)
      END DO
   81 ENORM=0.D0
      DO I=1,NBAS1
        DO J=1,I
          KP=NSTEP+(I*(I-1)/2)+J
          if (KP.GT.LPQ) call quit('KP > LPQ ! 2')
          TERM=C(I,N0)*S(KP)*C(J,N0)
          IF (I-J) 85,90,85
   85     TERM=2.D0*TERM
   90     ENORM=ENORM+TERM
        END DO
      END DO
      ENORM=1.D0/DSQRT(ENORM)
      DO I=1,NBAS1
        C(I,N0)=ENORM*C(I,N0)
      END DO
      IF (NSTART) 98,98,97
   97 NSTART=0
      GO TO 3
C
C     COMPUTE THE NEW EPS AND SAVE THE OLD
   98 EPSM1=EPS(N0)
      M=1
      GO TO 2
C
C     TEST CONVERGENCE
   96 SUP=0.D0
      DO I=1,NBAS1
        TERM=ABS(C(I,N0)-CM1(I))
        IF (TERM-SUP) 105,105,100
  100   SUP=TERM
  105   CONTINUE  
      END DO
      IF (SUP-THRE(N0)) 130,130,110
  110 IF (NITDIA-NDIAG) 115,120,120
  115 NITDIA=NITDIA+1
      GO TO 1
  120 THRE(N0)=2.D0*THRE(N0)
      IF (THRE(N0)-DGATH) 121,121,125
  121 NITDIA=1
      GO TO 1
  125 NALARM=2
      GO TO 140
  130 IF (THRE(N0)-.1D-07) 140,140,135
  135 THRE(N0)=THRE(N0)/2.D0
  140 CONTINUE 

      IF (IPR_RELSCF.GE.7) THEN
        WRITE(LUPRI,*) 'EIGEA after sizeof(S(:1)   =',size(S(:1))
        WRITE(LUPRI,*) 'EIGEA after sizeof(S(:LPQ))=',size(S(:LPQ))
      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ENERGY_RELSCF
     *(P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     DT  - total density matrix
C     DOS - open-shell density matrix
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "scfarr.h"
#include "relscf_dim.h"
#include "amfi_if.h"
c     logical relscf_verbose
      DIMENSION P(*),Q(*),S(*),U(*),T(*),DT(*),DOS(*),PCAP(*),
     *QCAP(*),FC(*),FO(*),CM1SCF(MXBF,*),CM2SCF(MXBF,*),
     *CSAV1(MXBF,*),CSAV2(MXBF,*),GMIN(*),SMIN(*),THRE(*),EPS(*),
     *EC(*),CM1(*),Z1(*),Z2(*),DELZ(*),
     *           CTRAN(*),PQN(*),ZETA(*),ZETINC(*),
     *ITRAN(*),NROW(*),NBVAR(*),IORCN(*)
      COMMON /CLOCKX/ CPU,CLK
C
c     relscf_verbose = .false.
      POT  = 0.D0
      POTN = 0.D0
      CIN  = 0.D0
      IF(NDIMPQ.NE.NDIMPQPAR) REWIND IDATA
      IPQ=NDIMPQ
      IF(NDIMPQ.EQ.NDIMPQPAR) IPQ=0
C
      K=0
      ENERG=DT(1)
      K=K+1
      IPQ=IPQ+1
      POT=POT+0.5D0*(ENERG*P(IPQ)*DT(1)-DOS(1)*Q(IPQ)*DOS(1))
      POTN=POTN+U(1)*ENERG
      CIN=CIN+T(1)*ENERG
      IF(N1T.EQ.1) GOTO 21
C
      DO I=2,N1T
        ENERG=DT(I)
        DO J=1,I-1
          K=K+1
          IPQ=IPQ+1
          TERM=ENERG*P(IPQ)*DT(J)-DOS(I)*Q(IPQ)*DOS(J)
          POT=POT+TERM
        END DO
        K=K+1
        IPQ=IPQ+1
        POT=POT+0.5D0* ( ENERG*P(IPQ)*DT(I)-DOS(I)*Q(IPQ)*DOS(I) )
        POTN=POTN+U(I)*ENERG
        CIN=CIN+T(I)*ENERG
      END DO
21    CONTINUE
      POT=POT-ZN*POTN
      ENERG=CIN+POT
**
      If (CIN.eq.0.d0.and.POT.eq.0.d0) then ! MI: This is a problematic case.. 0/0 may occur !!
        VIR=0.d0
      Else
        VIR=POT/CIN
      EndIf
**
      IF(NDIMPQ.NE.NDIMPQPAR) REWIND IDATA
C
      IF (NUMVAR.EQ.0) RETURN
      CPU=SECONDX()

      IF (IPR_RELSCF.GE.2) THEN
       WRITE (LUPRI,65) NBVAR1,ZETA(NBVAR1),ENERG,CIN,POT,VIR,NITSCF,CPU
   65  FORMAT(I5,5(1X,D16.10),I5,F14.2)
      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EVEN2(N,V,G,E,A,R,TT,AUXF,AUXG,AUXH)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EVEN2 - BERND HESS - V 1.0 - 5.2.86
C     CALCULATE EVEN2 OPERATORS
C
C
C     N       DIMENSION OF MATRICES
C     V       POTENTIAL MATRIX
C     G       MATRIX OF PVP OPERATOR. WILL CONTAIN EVEN2 OPERATORS
C             ON OUTPUT
C     E       RELATIVISTIC ENERGY (DIAGONAL)
C     A       A-FACTORS (DIAGONAL)
C     R       R-FACTORS (DIAGONAL)
C     TT      NONREL. KINETIC ENERGY (DIAGONAL)
C     AUXF,AUXG,AUXH  SCRATCH ARAYS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      DIMENSION V(*),G(*),E(N),R(N),A(N),TT(N),AUXF(N,N),
     *AUXG(N,N),AUXH(N,N)
      M=N
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXH(I,J)=0.D0
          AUXH(J,I)=0.D0
          V(IJ)=V(IJ)/(E(I)+E(J))
          G(IJ)=G(IJ)/(E(I)+E(J))
          AUXF(I,J)=A(I)*R(I)*G(IJ)*A(J)*A(J)
          AUXF(J,I)=A(J)*R(J)*G(IJ)*A(I)*A(I)
          AUXG(I,J)=R(I)*V(IJ)*A(J)
          AUXG(J,I)=R(J)*V(IJ)*A(I)
        END DO
      END DO
C
C     ARQA ARQA
C
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
      IF (IE.NE.0) STOP 129
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXG(I,J)=-0.5D0/TT(I)*G(IJ)*A(J)*R(J)
          AUXG(J,I)=-0.5D0/TT(J)*G(IJ)*A(I)*R(I)
        END DO
      END DO
C
C     ARQA AQRA
C
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXF(I,J)=A(I)*V(IJ)*A(J)*A(J)*R(J)
          AUXF(J,I)=A(J)*V(IJ)*A(I)*A(I)*R(I)
          AUXG(I,J)=-2.D0*TT(I)*R(I)*V(IJ)*A(J)
          AUXG(J,I)=-2.D0*TT(J)*R(J)*V(IJ)*A(I)
        END DO
      END DO
C
C     AQRA ARQA
C
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXG(I,J)=G(IJ)*A(J)*R(J)
          AUXG(J,I)=G(IJ)*A(I)*R(I)
        END DO
      END DO
C
C     AQRA AQRA
C
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
C
C     1/2 EW*W + 1/2 W*WE
C
      DO I=1,N
        DO J=1,N
          AUXH(I,J)=0.5*( AUXH(I,J)*E(I) + AUXH(I,J)*E(J) )
        END DO
      END DO
C
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXF(I,J)=A(I)*R(I)*G(IJ)*A(J)*E(J)*A(J)
          AUXF(J,I)=A(J)*R(J)*G(IJ)*A(I)*E(I)*A(I)
          AUXG(I,J)=R(I)*V(IJ)*A(J)
          AUXG(J,I)=R(J)*V(IJ)*A(I)
        END DO
      END DO
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXG(I,J)=-0.5D0/TT(I)*G(IJ)*A(J)*R(J)
          AUXG(J,I)=-0.5D0/TT(J)*G(IJ)*A(I)*R(I)
        END DO
      END DO
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXF(I,J)=A(I)*V(IJ)*R(J)*A(J)*E(J)*A(J)
          AUXF(J,I)=A(J)*V(IJ)*R(I)*A(I)*E(I)*A(I)
          AUXG(I,J)=-2.D0*TT(I)*R(I)*V(IJ)*A(J)
          AUXG(J,I)=-2.D0*TT(J)*R(J)*V(IJ)*A(I)
        END DO
      END DO
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXG(I,J)=G(IJ)*A(J)*R(J)
          AUXG(J,I)=G(IJ)*A(I)*R(I)
        END DO
      END DO
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
C
C     SYMMETRISIEREN
C
408   IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          G(IJ)=-0.5D0*(AUXH(I,J)+AUXH(J,I))
        END DO
      END DO
      RETURN
      END

C *** PVP V 1.0 - 19.1.84 - BERND HESS
C *** PVP V 1.1 -  5.6.86 - BERND HESS
      DOUBLE PRECISION FUNCTION EXTC(LAMBDA,AL,BE,L1,M1,N1,L2,M2,N2)
#include "implicit.h"
#include "relscf_crelop.h"
      INTEGER IS1(3),IS2(3)
C
C    CALCULATE ANGULAR AND RADIAL PART
C
      II=L1+L2
      JJ=M1+M2
      KK=N1+N2
      IMAX=II+JJ+KK+3
      IF (IMAX.LE.20) GOTO 2

C
C    ERROR BRANCH: ANGULAR MOMENTUM  > MAXIMUM GIVEN BY ARRAY GAM
C
1001  WRITE (6,1002) L1,M1,N1,L2,M2,N2,LAMBDA
1002  FORMAT(' ILLEGAL ANGULAR MOMENTUM (PVP)'/,
     *       ' L1,M1,N1,L2,M2,N2,LAMBDA PRINTED'/,1X,7I5)
      STOP 1002
      CALL QUIT('error output !')
C
C    COMPUTE INTEGRAL OVER DERIVATIVE OF THE FUNCTIONS
C
2     IS1(1)=L1
      IS1(2)=M1
      IS1(3)=N1
      IS2(1)=L2
      IS2(2)=M2
      IS2(3)=N2

      SUM=DER(1,IS1,IS2,AL,BE)+DER(2,IS1,IS2,AL,BE)+
     \    DER(3,IS1,IS2,AL,BE)
C
C     NORMALIZATION
C
      II=L1+L1
      JJ=M1+M1
      KK=N1+N1
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
      EX=-0.5D0*DFLOAT(II+JJ+KK+3)
      OV1=0.5D0*ANG*GA(II+JJ+KK+3)*((AL+AL)**EX)

      II=L2+L2
      JJ=M2+M2
      KK=N2+N2
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
      EX=-0.5D0*DFLOAT(II+JJ+KK+3)
      OV2=0.5D0*ANG*GA(II+JJ+KK+3)*((BE+BE)**EX)
 
      EXTC=SUM/DSQRT(OV1*OV2)

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
       SUBROUTINE EXTRAP(IAV,SCR,
     * P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ONE CONFIGURATION
C     SCF ITERATIONS AND EXTRAPOLATIONS
C     MODIFIED 17.7.87 - BERND HESS - AVERAGING/ORTHOGONALIZATION
C
C
C   C(MXBF,1)  -->   C(MXBF,*)
C   CM1SCF(MXBF,1) --> CM1SCF(MXBF,*)
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "relscf_dim.h"
#include "relscf_v10.h"
      DIMENSION SCR(*),
C    *          P(2),Q(2),S(2),U(2),T(2),DT(2),DOS(2),PCAP(2),
     *          P(*),Q(*),S(*),U(*),T(*),DT(*),DOS(*),PCAP(*),
CMI  *QCAP(2),FC(2),FO(2),CM1SCF(MXBF,1),CM2SCF(MXBF,1),
     *QCAP(*),FC(*),FO(*),CM1SCF(MXBF,*),CM2SCF(MXBF,*),
     *CSAV1(MXBF,1),CSAV2(MXBF,1),GMIN(2),SMIN(2),THRE(*),EPS(2),
     *EC(2),CM1(2),Z1(2),Z2(2),DELZ(2),
CMI  *C(MXBF,1),CTRAN(2),PQN(*),ZETA(*),ZETINC(*),
     *C(MXBF,*),CTRAN(2),PQN(*),ZETA(*),ZETINC(*),
     *ITRAN(2),NROW(*),NBVAR(*),IORCN(*)
C
      NBIAS=BIAS
      THREM=0.D0
      DO I=1,NSHT
        IF (THRE(I)-THREM) 10,10,5
    5   THREM=THRE(I)
   10   CONTINUE
      END DO

      IF (NSVD.EQ.1) THREM=1.D-8
      THRSCF=THREM*(2.D0**NBIAS)
      IF (THRSCF-SCFAT) 20,15,15
   15 NALARM=3
      GO TO 130
   20 SUP=0.D0
      NSTEP=0
      DO I=1,NSYM
        NBAS1=NBAS(I)
        NSH=NCSH(I)+NOSH(I)
        DO J=1,NSH
          J1=J+NSTEP
          DO M=1,NBAS1
CMI          ... complains of J1 in CM1SCF !
            TERM=ABS(C(M,J1)-CM1SCF(M,J1))
            IF (TERM-SUP) 45,45,40
   40       SUP=TERM
   45       CONTINUE 
          END DO
        END DO
        NSTEP=NSTEP+NSH
      END DO

      IF (SUP-THRSCF) 50,50,55
   50 NCONV=1
      GO TO 130
   55 IF (NITSC2-NXTRP) 65,60,60
   60 BIAS=BIAS+1.D0
      NITSC2=0
   65 IF (NITSCF-MXXTRP) 75,70,70
   70 NALARM=4
      GO TO 130
   75 NITSCF=NITSCF+1
      NITSC2=NITSC2+1
C
C     MO VECTOR EXTRAPOLATION
C
      IF (NITSC1-4) 80,85,85
   80 NITSC1=NITSC1+1
      GO TO 130
85    CONTINUE
C
C     START EXTRAPOLATION
C
      NSTEP1=0
      NSTEP2=0
      DO I=1,NSYM
        NBAS1=NBAS(I)
        NSH=NCSH(I)+NOSH(I)
        DO J=1,NSH
          J1=J+NSTEP1
          DO M=1,NBAS1
            SUP=C(M,J1)
            ENORM=CM1SCF(M,J1)
            C(M,J1)=SUP+ENORM
            CM2SCF(M,J1)=SUP-ENORM
          END DO
C
C     RENORMALIZE
C
          ENORM=0.D0
          SUP=0.D0
          DO M=1,NBAS1
            DO K=1,M
              KP=NSTEP2+(M*(M-1)/2)+K

              if (KP.GT.LPQ) call quit('KP > LPQ !')

              TERP=C(M,J1)*S(KP)*C(K,J1)
              TERM=CM2SCF(M,J1)*S(KP)*CM2SCF(K,J1)
              IF (M.EQ.K) GOTO 116
              TERP=2.D0*TERP
              TERM=2.D0*TERM
  116         ENORM=ENORM+TERP
              SUP=SUP+TERM
            END DO
          END DO
          IF (SUP.GE.ENORM) THEN

              ENORM=SUP
          DO K=1,NBAS1
            C(K,J1)=CM2SCF(K,J1)
          END DO
          ENDIF
          ENORM=1.D0/DSQRT(ENORM)
          DO M=1,NBAS1
            C(M,J1)=ENORM*C(M,J1)
          END DO
C
C     ORTHOGONALIZE
C
          if (1+NSTEP2.GT.LPQ) call quit('1+NSTEP2 > LPQ !')

          CALL SCHMID(C(1,NSTEP1+1),J,NBAS1,MXBF,S(1+NSTEP2),SCR)

        END DO
        NSTEP2=NSTEP2+N1(I)
        NSTEP1=NSTEP1+NSH
      END DO
      NITSC1=1
  130 RETURN
      END

C***** NAME FILLMA
      SUBROUTINE FILLMA(N,S,OVE)
#include "implicit.h"
#include "priunit.h"
      DIMENSION S(*),OVE(N,N)
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          OVE(I,J)=S(IJ)
          OVE(J,I)=S(IJ)
        END DO
      END DO
      RETURN
      END

      DOUBLE PRECISION FUNCTION GAM(M)
#include "implicit.h"
#include "relscf_crelop.h"
      
      IF (MOD(M,2).EQ.0) GOTO 10
      MA=(M+1)/2
      G=1.D0
      IF (MA.EQ.1) GOTO 11
      DO I=2,MA
        G=G*DFLOAT(I-1)
      END DO
      GAM=G
      RETURN
10    MA=M
      G=SQPI
      IF (MA.EQ.0) GOTO 11
      DO I=1,MA,2
        G=G*0.5D0*DFLOAT(I)
      END DO
11    GAM=G
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GAUSS(N,X,W)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C Array dim.fixes:
C------------------
C
C  X(2) -> X(*)
C  W(2) -> W(*)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
CMI   DIMENSION X(2), W(2)
      DIMENSION X(*), W(*)
      DIMENSION GPT(166),GWT(166),P3(2),W3(2)
      EQUIVALENCE (P3(1),GPT(165)),(W3(1),GWT(165))
      INTEGER NS(16),LOCS(16)
      REAL*8 P2(1),P4(2),P6(3),P8(4),P10(5),P12(6),P14(7),P16(8),P20(10)
     1      ,P24(12),P28(14),P32(16),P40(20),P48(24),P64(32)
      REAL*8 W2(1),W4(2),W6(3),W8(4),W10(5),W12(6),W14(7),W16(8),W20(10)
     1      ,W24(12),W28(14),W32(16),W40(20),W48(24),W64(32)
      EQUIVALENCE (P2(1),GPT(1)),(P4(1),GPT(2)),(P6(1),GPT(4)),
     1    (P8(1),GPT(7)),(P10(1),GPT(11)),(P12(1),GPT(16)),
     2    (P14(1),GPT(22)),(P16(1),GPT(29)),(P20(1),GPT(37)),
     3    (P24(1),GPT(47)),(P28(1),GPT(59)),(P32(1),GPT(73)),
     4    (P40(1),GPT(89)),(P48(1) ,GPT(109)),(P64(1),GPT(133))
      EQUIVALENCE (W2(1),GWT(1)),(W4(1),GWT(2)),(W6(1),GWT(4)),
     1    (W8(1),GWT(7)),(W10(1),GWT(11)),(W12(1),GWT(16)),
     2    (W14(1),GWT(22)),(W16(1),GWT(29)),(W20(1),GWT(37)),
     3    (W24(1),GWT(47)),(W28(1),GWT(59)),(W32(1),GWT(73)),
     4    (W40(1),GWT(89)),(W48(1) ,GWT(109)),(W64(1),GWT(133))
      DATA P 2/0.577350269189626D 00/
      DATA W 2/0.999999999999999D 00/
      DATA P 3/0.774596669241483D 00,0.000000000000000D 00/
      DATA W 3/0.555555555555556D 00,0.888888888888889D 00/
      DATA P 4/0.861136311594052D 00,0.339981043584856D 00/
      DATA W 4/0.347854845137454D 00,0.652145154862546D 00/
      DATA P 6/0.932469514203152D 00,0.661209386466264D 00,
     1         0.238619186083197D 00/
      DATA W 6/0.171324492379170D 00,0.360761573048139D 00,
     1         0.467913934572691D 00/
      DATA P 8/0.960289856497536D 00,0.796666477413627D 00,
     1         0.525532409916329D 00,0.183434642495650D 00/
      DATA W 8/0.101228536290376D 00,0.222381034453374D 00,
     1         0.313706645877887D 00,0.362683783378362D 00/
      DATA P10/0.973906528517172D 00,0.865063366688984D 00,
     1         0.679409568299024D 00,0.433395394129247D 00,
     2         0.148874338981631D 00/
      DATA W10/0.066671344308688D 00,0.149451349150580D 00,
     1         0.219086362515982D 00,0.269266719309996D 00,
     2         0.295524224714753D 00/
      DATA P12/0.981560634246719D 00,0.904117256370475D 00,
     1         0.769902674194305D 00,0.587317954286617D 00,
     2         0.367831498998180D 00,0.125233408511469D 00/
      DATA W12/0.047175336386512D 00,0.106939325995318D 00,
     1         0.160078328543346D 00,0.203167426723066D 00,
     2         0.233492536538355D 00,0.249147045813403D 00/
      DATA P14/0.986283808696812D 00,0.928434883663573D 00,
     1         0.827201315069765D 00,0.687292904811685D 00,
     2         0.515248636358154D 00,0.319112368927890D 00,
     3         0.108054948707344D 00/
      DATA W14/0.035119460331752D 00,0.080158087159760D 00,
     1         0.121518570687903D 00,0.157203167158193D 00,
     2         0.185538397477938D 00,0.205198463721296D 00,
     3         0.215263853463158D 00/
      DATA P16/0.989400934991650D 00,0.944575023073232D 00,
     1         0.865631202387832D 00,0.755404408355003D 00,
     2         0.617876244402644D 00,0.458016777657227D 00,
     3         0.281603550779259D 00,0.095012509837637D 00/
      DATA W16/0.027152459411754D 00,0.062253523938648D 00,
     1         0.095158511682493D 00,0.124628971255534D 00,
     2         0.149595988816577D 00,0.169156519395002D 00,
     3         0.182603415044923D 00,0.189450610455068D 00/
      DATA P20/0.993128599185095D 00,0.963971927277914D 00,
     1         0.912234428251326D 00,0.839116971822219D 00,
     2         0.746331906460151D 00,0.636053680726515D 00,
     3         0.510867001950827D 00,0.373706088715419D 00,
     4         0.227785851141645D 00,0.076526521133497D 00/
      DATA W20/0.017614007139152D 00,0.040601429800387D 00,
     1         0.062672048334109D 00,0.083276741576705D 00,
     2         0.101930119817240D 00,0.118194531961518D 00,
     3         0.131688638449177D 00,0.142096109318382D 00,
     4         0.149172986472604D 00,0.152753387130726D 00/
      DATA P24/0.995187219997021D 00,0.974728555971309D 00,
     1         0.938274552002733D 00,0.886415527004401D 00,
     2         0.820001985973903D 00,0.740124191578554D 00,
     3         0.648093651936975D 00,0.545421471388839D 00,
     4         0.433793507626045D 00,0.315042679696163D 00,
     5         0.191118867473616D 00,0.064056892862606D 00/
      DATA W24/0.012341229799987D 00,0.028531388628934D 00,
     1         0.044277438817420D 00,0.059298584915437D 00,
     2         0.073346481411080D 00,0.086190161531953D 00,
     3         0.097618652104114D 00,0.107444270115966D 00,
     4         0.115505668053726D 00,0.121670472927803D 00,
     5         0.125837456346828D 00,0.127938195346752D 00/
      DATA P28/0.996442497573954D 00,0.981303165370873D 00,
     1         0.954259280628938D 00,0.915633026392132D 00,
     2         0.865892522574395D 00,0.805641370917179D 00,
     3         0.735610878013632D 00,0.656651094038865D 00,
     4         0.569720471811402D 00,0.475874224955118D 00,
     5         0.376251516089079D 00,0.272061627635178D 00,
     6         0.164569282133381D 00,0.055079289884034D 00/
      DATA W28/0.009124282593094D 00,0.021132112592771D 00,
     1         0.032901427782304D 00,0.044272934759004D 00,
     2         0.055107345675717D 00,0.065272923966999D 00,
     3         0.074646214234569D 00,0.083113417228901D 00,
     4         0.090571744393033D 00,0.096930657997930D 00,
     5         0.102112967578061D 00,0.106055765922846D 00,
     6         0.108711192258294D 00,0.110047013016475D 00/
      DATA P32/0.997263861849481D 00,0.985611511545268D 00,
     1         0.964762255587506D 00,0.934906075937740D 00,
     2         0.896321155766052D 00,0.849367613732570D 00,
     3         0.794483795967942D 00,0.732182118740290D 00,
     4         0.663044266930215D 00,0.587715757240762D 00,
     5         0.506899908932229D 00,0.421351276130635D 00,
     6         0.331868602282128D 00,0.239287362252137D 00,
     7         0.144471961582796D 00,0.048307665687738D 00/
      DATA W32/0.007018610009470D 00,0.016274394730906D 00,
     1         0.025392065309262D 00,0.034273862913021D 00,
     2         0.042835898022227D 00,0.050998059262376D 00,
     3         0.058684093478535D 00,0.065822222776362D 00,
     4         0.072345794108848D 00,0.078193895787070D 00,
     5         0.083311924226947D 00,0.087652093004404D 00,
     6         0.091173878695764D 00,0.093844399080804D 00,
     7         0.095638720079275D 00,0.096540088514728D 00/
      DATA P40/0.998237709710559D 00,0.990726238699457D 00,
     1         0.977259949983774D 00,0.957916819213792D 00,
     2         0.932812808278676D 00,0.902098806968874D 00,
     3         0.865959503212259D 00,0.824612230833312D 00,
     4         0.778305651426519D 00,0.727318255189927D 00,
     5         0.671956684614179D 00,0.612553889667980D 00,
     6         0.549467125095128D 00,0.483075801686179D 00,
     7         0.413779204371605D 00,0.341994090825758D 00,
     8         0.268152185007254D 00,0.192697580701371D 00,
     9         0.116084070675255D 00,0.038772417506051D 00/
      DATA W40/0.004521277098533D 00,0.010498284531153D 00,
     1         0.016421058381908D 00,0.022245849194167D 00,
     2         0.027937006980023D 00,0.033460195282548D 00,
     3         0.038782167974472D 00,0.043870908185673D 00,
     4         0.048695807635072D 00,0.053227846983937D 00,
     5         0.057439769099391D 00,0.061306242492929D 00,
     6         0.064804013456601D 00,0.067912045815234D 00,
     7         0.070611647391287D 00,0.072886582395804D 00,
     8         0.074723169057968D 00,0.076110361900626D 00,
     9         0.077039818164248D 00,0.077505947978425D 00/
      DATA P48/0.998771007252426D 00,0.993530172266351D 00,
     1         0.984124583722827D 00,0.970591592546247D 00,
     2         0.952987703160431D 00,0.931386690706554D 00,
     3         0.905879136715570D 00,0.876572020274248D 00,
     4         0.843588261624393D 00,0.807066204029443D 00,
     5         0.767159032515740D 00,0.724034130923815D 00,
     6         0.677872379632664D 00,0.628867396776514D 00,
     7         0.577224726083973D 00,0.523160974722233D 00,
     8         0.466902904750958D 00,0.408686481990717D 00,
     9         0.348755886292161D 00,0.287362487355455D 00,
     A         0.224763790394689D 00,0.161222356068892D 00,
     B         0.097004699209463D 00,0.032380170962869D 00/
      DATA W48/0.003153346052306D 00,0.007327553901276D 00,
     1         0.011477234579234D 00,0.015579315722944D 00,
     2         0.019616160457355D 00,0.023570760839324D 00,
     3         0.027426509708357D 00,0.031167227832798D 00,
     4         0.034777222564770D 00,0.038241351065831D 00,
     5         0.041545082943465D 00,0.044674560856694D 00,
     6         0.047616658492490D 00,0.050359035553854D 00,
     7         0.052890189485194D 00,0.055199503699984D 00,
     8         0.057277292100403D 00,0.059114839698396D 00,
     9         0.060704439165894D 00,0.062039423159893D 00,
     A         0.063114192286254D 00,0.063924238584648D 00,
     B         0.064466164435950D 00,0.064737696812684D 00/
      DATA P64/0.999305041735772D 00,0.996340116771955D 00,
     1         0.991013371476744D 00,0.983336253884626D 00,
     2         0.973326827789911D 00,0.961008799652054D 00,
     3         0.946411374858403D 00,0.929569172131939D 00,
     4         0.910522137078503D 00,0.889315445995114D 00,
     5         0.865999398154093D 00,0.840629296252580D 00,
     6         0.813265315122797D 00,0.783972358943341D 00,
     7         0.752819907260532D 00,0.719881850171611D 00,
     8         0.685236313054233D 00,0.648965471254657D 00,
     9         0.611155355172393D 00,0.571895646202634D 00,
     A         0.531279464019894D 00,0.489403145707053D 00,
     B         0.446366017253464D 00,0.402270157963992D 00,
     C         0.357220158337668D 00,0.311322871990211D 00,
     D         0.264687162208767D 00,0.217423643740007D 00,
     E         0.169644420423993D 00,0.121462819296120D 00,
     F         0.072993121787799D 00,0.024350292663424D 00/
      DATA W64/0.001783280721696D 00,0.004147033260562D 00,
     1         0.006504457968978D 00,0.008846759826364D 00,
     2         0.011168139460131D 00,0.013463047896719D 00,
     3         0.015726030476025D 00,0.017951715775697D 00,
     4         0.020134823153530D 00,0.022270173808383D 00,
     5         0.024352702568711D 00,0.026377469715055D 00,
     6         0.028339672614259D 00,0.030234657072402D 00,
     7         0.032057928354851D 00,0.033805161837142D 00,
     8         0.035472213256882D 00,0.037055128540240D 00,
     9         0.038550153178616D 00,0.039953741132720D 00,
     A         0.041262563242623D 00,0.042473515123653D 00,
     B         0.043583724529323D 00,0.044590558163756D 00,
     C         0.045491627927418D 00,0.046284796581314D 00,
     D         0.046968182816210D 00,0.047540165714830D 00,
     E         0.047999388596458D 00,0.048344762234803D 00,
     F         0.048575467441503D 00,0.048690957009140D 00/
      DATA      NS /2,4,6,8,10,12,14,16,20,24,28,32,40,48,64,3/
      DATA      LOCS /1,2,4,7,11,16,22,29,37,47,59,73,89,109,133,165/
    1 FORMAT ('0****** ERROR MESSAGE   GAUSS/',I4)
      DO I = 1, 16
        L1 = LOCS(I)
        IF (N.EQ.NS(I)) GO TO 200
      END DO

      WRITE (6,1) N

      CALL QUIT('relscf/GAUSS: error exit')

  200 L2 = (N + 1) / 2 - 1 + L1
      N1 = 1
      N2 = N
      DO L = L1, L2
CMI ....   for array 'x', upper bound of dimension 1 exceeded 
        X(N2) = GPT(L)
        X(N1) = - GPT(L)
        W(N1) = GWT(L)
CMI  ...  for array 'w', upper bound of dimension 1 exceeded
        W(N2) = GWT(L)
        N1 = N1 + 1
        N2 = N2 - 1
      END DO
      RETURN
      END

C     EXPONENTIATIONS REMOVED
      FUNCTION GINY(NY,NAB,NCD,ZAB,ZCD)
#include "implicit.h"
#include "relscf_dim.h"
C
C     CALCULATES THE INTEGRAL I(NY) NEEDED FOR THE SUPERMATRICES P AND Q
C     WHEN A GAUSSIAN BASIS SET IS USED
C
      TABCD=ZAB/ZCD
      NABNY=NAB-NY-1
      NCDNY=NCD+NY
Cvec  VNAB=FACTO(NABNY)/ZAB**(0.5D0*(NABNY+1))
      IPOW=(NABNY+1)/2
      JPOW=NABNY+1-2*IPOW
      RES=1.D0
      IF (JPOW.EQ.1) RES=DSQRT(ZAB)
      GOTO (110,111,112,113,114,115,116,117,118,119,120),IPOW
      IF (IPOW.EQ.0) GOTO 199
      WRITE (6,198) IPOW,NY,NAB,NCD
198   FORMAT(' GINY198 - POWER OVERFLOW ',4I10)
      STOP 'GINY'
120    RES=RES*ZAB
119    RES=RES*ZAB
118    RES=RES*ZAB
117    RES=RES*ZAB
116    RES=RES*ZAB
115    RES=RES*ZAB
114    RES=RES*ZAB
113    RES=RES*ZAB
112    RES=RES*ZAB
111    RES=RES*ZAB
110    RES=RES*ZAB
199   VNAB=FACTO(NABNY)/RES

      ICDPOW=(NCDNY+1)/2
      JCDPOW=NCDNY+1-2*ICDPOW
      RES=1.D0
      IF (JCDPOW.EQ.1) RES=DSQRT(ZCD)
      GOTO (210,211,212,213,214,215,216,217,218,219,220),ICDPOW
      IF (ICDPOW.EQ.0) GOTO 299
      WRITE (6,198) ICDPOW,NY,NAB,NCD
298   FORMAT(' GINY298 - POWER OVERFLOW ',4I10)
      STOP 'GINY'
220    RES=RES*ZCD
219    RES=RES*ZCD
218    RES=RES*ZCD
217    RES=RES*ZCD
216    RES=RES*ZCD
215    RES=RES*ZCD
214    RES=RES*ZCD
213    RES=RES*ZCD
212    RES=RES*ZCD
211    RES=RES*ZCD
210    RES=RES*ZCD
299   VNCD=FACTO(NCDNY)/RES
      TERM1=0.D0
      DO I=1,NABNY,2
        POW=(TABCD/(1.D0+TABCD))
        RES=1.D0
        IPOW=(I-1)/2
        GOTO (10,11,12,13,14,15,16,17,18,19,20),IPOW
        IF (IPOW.EQ.0) GOTO 99
        WRITE (6,98) IPOW,NY,NAB,NCD
98      FORMAT(' GINY - POWER OVERFLOW ',4I10)
        STOP 'GINY'
20      RES=RES*POW
19      RES=RES*POW
18      RES=RES*POW
17      RES=RES*POW
16      RES=RES*POW
15      RES=RES*POW
14      RES=RES*POW
13      RES=RES*POW
12      RES=RES*POW
11      RES=RES*POW
10      RES=RES*POW
99      TERM1=TERM1+FACTO(I+NCDNY-1)/(FACTO(I)*FACTO(NCDNY))*RES
      END DO

      POW=1.D0+TABCD
      RES=1.D0
      IF (JCDPOW.EQ.1) RES=DSQRT(POW)
      GOTO (310,311,312,313,314,315,316,317,318,319,320),ICDPOW
      IF (ICDPOW.EQ.0) GOTO 399
320    RES=RES*POW
319    RES=RES*POW
318    RES=RES*POW
317    RES=RES*POW
316    RES=RES*POW
315    RES=RES*POW
314    RES=RES*POW
313    RES=RES*POW
312    RES=RES*POW
311    RES=RES*POW
310    RES=RES*POW
399   CABCD=TERM1/RES
      GINY=VNAB*VNCD*CABCD
      RETURN
      END

C     V 2.0 - VECTORIZED - BERND HESS - MUST BE FORTRAN77
CICI
CICI  LOOP-STRUCTURE CHANGED TO AVAOID COMPILER-DEPENDENT PERFORMANCE
CICI
      SUBROUTINE HAMIL
     *(P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C
C       ONE CONFIGURATION
C
C MI:
C  C(MXBF,1) -> C(MXBF,*)
C  SMIN(MXBF,1) ->  SMIN(MXBF,*)
C  GMIN(MXBF,1) ->  GMIN(MXBF,*)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "scfarr.h"
#include "relscf_dim.h"
      DIMENSION P(*),Q(*),S(*),U(*),T(*),DT(*),DOS(*),PCAP(*),
     *QCAP(*),FC(*),FO(*),CM1SCF(MXBF,1),CM2SCF(MXBF,1),
     *CSAV1(MXBF,1),CSAV2(MXBF,1),                THRE(*),EPS(2),
     *EC(*),CM1(*),Z1(*),Z2(2),DELZ(2),GMIN(MXBF,*),SMIN(MXBF,*),
     *C(MXBF,*),CTRAN(2),PQN(*),ZETA(*),ZETINC(*),
     *ITRAN(2),NROW(*),NBVAR(*),IORCN(*)
C
C     COMPUTE IN THE FOLLOWING ORDER THE MATRICES QCAP,PCAP,SMIN,QMIN,
C     FC,FO
C     QCAP MEAN THE Q MATRIX OF ROOTHAAN
C     PCAP MEAN THE P MATRIX OF ROOTHAAN
C
C     COMPUTE QCAP AND PCAP
      N1T=0
      DO I=1,NSYM
        NBAS1=NBAS(I)
        NSTEP1=NBAS1*(NBAS1+1)/2
        N1T=N1T+NSTEP1
        N1(I)=NSTEP1
      END DO
      IF (NDIMPQ.NE.NDIMPQPAR) THEN
        write(lupri,*) 'HAMIL: NDIMPQ=',NDIMPQ
          CALL QUIT('HAMIL: NDIMPQ.NE.NDIMPQPAR!')
      ENDIF

      IF(NDIMPQ.EQ.NDIMPQPAR) IPQ=0
      DO I=1,N1T
        PCAP(I)=0.D0
        QCAP(I)=0.D0
      END DO

      TEMP1=PCAP(1)
      TEMP2=QCAP(1)
      IPQ=IPQ+1
      PCAP(1)=TEMP1+P(IPQ)*DT(1)
      QCAP(1)=TEMP2+Q(IPQ)*DOS(1)
      IF(N1T.EQ.1) GOTO 36
C
      DO I=2,N1T
        NTO=I-1
        TEMP1=PCAP(I)
        TEMP2=QCAP(I)
        DO J=1,NTO
          IPQ=IPQ+1
          TEMP1=TEMP1+P(IPQ)*DT(J)
          TEMP2=TEMP2+Q(IPQ)*DOS(J)
          PCAP(J)=PCAP(J)+P(IPQ)*DT(I)
          QCAP(J)=QCAP(J)+Q(IPQ)*DOS(I)
        END DO
        IPQ=IPQ+1
        PCAP(I)=TEMP1+P(IPQ)*DT(I)
        QCAP(I)=TEMP2+Q(IPQ)*DOS(I)
      END DO

36    IF (IPQ.GT.NDIMPQ) THEN
         write(LUPRI,*)'relscf/HAMIL: IPQ,NDIMPQ:',IPQ,NDIMPQ
         CALL QUIT('relscf/HAMIL: NDIMPQ < IPQ !')
      ENDIF
C
C      COMPUTES SMIN AND GMIN
      NSTEP1=0
      NSTEP2=0
      DO I=1,NSYM
        NSH=NCSH(I)+NOSH(I)
        NBAS1=NBAS(I)
        DO J=1,NSH
          J1=NSTEP2+J

              DO M=1,NBAS1-1
            FACT1=0.D0
            FACT2=0.D0

                  DO N=1,M
              K=(M*(M-1)/2)+N+NSTEP1
              FACT1=FACT1+S(K)*C(N,J1)
              FACT2=FACT2+QCAP(K)*C(N,J1)
            END DO

                  IHLP=M+1
            DO N=IHLP,NBAS1
              K=(N*(N-1)/2)+M+NSTEP1
              FACT1=FACT1+S(K)*C(N,J1)
              FACT2=FACT2+QCAP(K)*C(N,J1)
            END DO
            SMIN(M,J1)=FACT1

                  GMIN(M,J1)=FACT2
          END DO
CICI
CICI  CASE M=NBAS1
CICI
          FACT1=0.D0
          FACT2=0.D0
          DO N=1,NBAS1
            K=(M*(M-1)/2)+N+NSTEP1
            FACT1=FACT1+S(K)*C(N,J1)
            FACT2=FACT2+QCAP(K)*C(N,J1)
          END DO
          SMIN(NBAS1,J1)=FACT1
          GMIN(NBAS1,J1)=FACT2
        END DO
C
        NSTEP2=NSTEP2+NSH
        NSTEP1=NSTEP1+N1(I)
      END DO
C
C     COMPUTE FC AND FO
C
      K=1
      NSTEP=0
      DO I=1,NSYM
        NOSH1=NOSH(I)
        NSH=NCSH(I)+NOSH1
        FACT2=OCCUP(I)
        FACT1=OCCUCS(I)-FACT2
cfact
        if(fact1.eq.0.0) fact1=fact2
        NBAS1=NBAS(I)
        DO M=1,NBAS1
          DO N=1,M
            TERM1=0.D0
            TERM2=0.D0
CICI
CICI  I SEE NO EASY WAY OF COMPILER-INDEPENDENT VECTORIZATION
CICI
            DO J=1,NSH
              J1=J+NSTEP
              IF (J-NSH) 90,80,90
   80         IF (NOSH1) 85,90,85
   85         TERM2=SMIN(M,J1)*GMIN(N,J1)+GMIN(M,J1)*SMIN(N,J1)
              GO TO 95
   90         TERM1=TERM1+SMIN(M,J1)*GMIN(N,J1)+GMIN(M,J1)*SMIN(N,J1)
   95         CONTINUE
            END DO

                  HPPC=PCAP(K)+T(K)-ZN*U(K)
            FC(K)=HPPC+TERM2*FACT2/FACT1
            FO(K)=HPPC-QCAP(K)+TERM1*OCCUCS(I)/FACT1
            K=K+1
          END DO
        END DO
        NSTEP=NSTEP+NSH
      END DO

      IF(NDIMPQ.NE.NDIMPQPAR) REWIND IDATA
      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HFMAIN
     *    (P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,REL,Z100,
     *    S100,T100,U100,ZADD,C,ITRAN,CTRAN,NROW,PQN,ZETA,
     &    NBVAR,ZETINC,IORCN,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     MAIN PROGRAM ONE CONFIGURATION LCAO ATOMIC SCF CALCULATION
C
C  Called from AT3DIM where most parameters are allocated through
C  MEMGET.
C
C
C   C(MXBF,1) -> C(MXBF,*)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "relscf_dim.h"
#include "relscf_v10.h"
#include "scfarr.h"

CMI   DIMENSION P(LQM),Q(LQM),S(LPQ),U(LPQ),T(LPQ),DT(LPQ,2),
CMI   DIMENSION P(LQM),Q(LQM),S(LPQ),U(*),T(LPQ),DT(LPQ,*),
CMI   DIMENSION P(LQM),Q(LQM),S(*),U(*),T(LPQ),DT(LPQ,*),
      DIMENSION P(*),Q(*),S(*),U(*),T(LPQ),DT(LPQ,*),
CMI   DIMENSION P(LQM),Q(LQM),S(LPQ),U(0:*),T(LPQ),DT(LPQ,*),
C    &          DOS(LPQ,2),PCAP(LPQ),QCAP(LPQ),FC(LPQ),FO(LPQ),
     &          DOS(LPQ,*),PCAP(LPQ),QCAP(LPQ),FC(LPQ),FO(LPQ),
     &          CM1SCF(MXBF,NOR),CM2SCF(MXBF,NOR),CSAV1(MXBF,NOR),
     &          CSAV2(MXBF,NOR),GMIN(*),SMIN(*),THRE(NOR),EPS(NOR),
     &          EC(MXBF),CM1(MXBF),Z1(MXBF),Z2(MXBF),DELZ(MXBF),
     &          REL(*),Z100(MX100),S100(*),T100(*),U100(*),
CMI  &          ZADD(MX100,4),C(MXBF,1),CTRAN(LDIMC),PQN(MXB),ZETA(MXB),
     &          ZADD(MX100,4),C(MXBF,*),CTRAN(LDIMC),PQN(MXB),ZETA(MXB),
     &          ZETINC(MXB),ITRAN(*),NROW(MXB),NBVAR(MXB),IORCN(MXB)
      DIMENSION WORK(LWORK)
      COMMON /CLOCKX/ CPU,CLK

      DIMENSION ZSAV(4),ESAV(4),EOS(3)
      external AMFIALARM
   
      LOGICAL ZEROELECTRONS

#include "amfi_if.h"
#include "memint.h"
      CALL RELOP
      ncall=1
      IOS=0
      IPR=1
      IDEBUG=0
      JPARM=1
      KPARM=IPARM(JPARM)
      FACTO(1)=1.D0
      FACTO(2)=1.D0
      FATT(1)=1.D0
      FATT(2)=1.D0
      DO I=3,16
        IM1=I-1
        FACTO(I)=IM1*FACTO(I-2)
        FATT(I)=IM1*FATT(IM1)
      END DO
      DO I=17,21
        IM1=I-1
        FATT(I)=IM1*FATT(IM1)
      END DO
      IF(NVAR.NE.0) WRITE (6,1000)
 1000 FORMAT('1INTERMEDIATE OPTIMIZATION DATA',//)
      NUMVAR=0
      NALARM=0
      MVAR=0
      IEX=0
      NRUN=1
C=========================================================================
C
C             ITERATION STARTS - CALCULATION OF INTEGRALS
C
C=========================================================================
   15 CONTINUE 
      IF (IPR_RELSCF.GE.0) THEN
        WRITE(LUPRI,*) ' ### SCF ITERATIONS           ###'
        IF (KPARM.EQ.0) WRITE (LUPRI,*)
     &          ' ### NON-RELATIVISTIC APPROX. ###'
        IF (KPARM.EQ.1) WRITE (LUPRI,*)
     &           ' ### SQR APPROX.              ###'
        IF (KPARM.EQ.2) WRITE (LUPRI,*) 
     &           ' ### EV APPROX.               ###'
      ENDIF

      NDIM=10000
      CALL MEMGET('REAL',KUN100,NDIM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTN100,NDIM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHCORR,NDIM,WORK,KFREE,LFREE)
C.....One-electron integrals
      CALL OEISG(KPARM,
     * P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,REL,Z100,
     *    S100,T100,U100,ZADD,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN,
     &    NDIM,WORK(KUN100),WORK(KTN100),WORK(KHCORR))
      CALL MEMREL('HFMAIN-after QEISG',WORK,KWORK,KWORK,KFREE,LFREE)
      IF (ncall.eq.1) THEN
C.....Two-electron integrals
        CALL TEIG(P,Q,NSYM,NBAS,NOSH,ZETA,AJMN,U,OCCUP,ndimpq)
      ENDIF
      ncall=2
 1002 CONTINUE
C
C     SET ITERATION PARAMETERS TO ZERO
C
      CALL DCOPY(NSHT,0.1D-05,0,THRE,1)
C      DO I=1,NSHT
C        THRE(I)=0.1D-05
C      END DO
      NITSCF=1
      NITSC1=1
      NITSC2=1
      NCONV=0
C
C     CONSTRUCT DENSITY MATRIX FROM VECTOR
C
20    CONTINUE
      IF (IPR_RELSCF.GE.5)  THEN
        WRITE(LUPRI,'(/,4x,a,/)') 'Calling DENSI'
        WRITE(LUPRI,*) 'sizeof(S(:1)   =',size(S(:1))
        WRITE(LUPRI,*) 'sizeof(S(:LPQ))=',size(S(:LPQ))
      ENDIF
      CALL DENSI(P,Q,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
C
C     CALCULATE ENERGY, IF APPLICABLE
C
      IF (IPR_RELSCF.GE.5)  THEN
        WRITE(LUPRI,'(/,4x,a,/)') 'Calling ENERGY_RELSCF'
      ENDIF

      CALL ENERGY_RELSCF(
     *    P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)

      IF (IPR_RELSCF.GE.0) then
        write(LUPRI,"(2x,i3,a,1x,f26.12)") nitscf,
     & '. iteration,  total energy:', energ
      ENDIF
      IF (IPR_RELSCF.GE.15) THEN
C......Print intermediate iteration results
       WRITE(LUPRI,'(/,2x,a)') '...extra printout in advance:'
       CALL PREN(KPARM,
     *    P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
      ENDIF

      IF (NCONV.EQ.0) GOTO 30
C
C     NOW WE ARE CONVERGED. ANY ADDITIONAL OPERATORS?
C
      IOS=0
      IF (NRUN.EQ.1.AND.IPR.GT.0) THEN
        IF (JPARM.LT.3 .AND. KPARM.NE.IPARM(JPARM+1)) THEN

          if (IPR_RELSCF.GE.0) then
            IF (KPARM.EQ.0) WRITE (LUPRI,*)
     &          ' ### NON-RELATIVISTIC APPROX. ###'
            IF (KPARM.EQ.1) WRITE (LUPRI,*)
     &           ' ### SQR APPROX.              ###'
            IF (KPARM.EQ.2) WRITE (LUPRI,*) 
     &           ' ### EV APPROX.               ###'
          endif

          CALL PREN(KPARM,
     *    P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
        ELSE

          if (IPR_RELSCF.GE.0) then
            IF (KPARM.EQ.0) WRITE (LUPRI,*)
     &      ' ### NON-RELATIVISTIC RESULT  ###'
            IF (KPARM.EQ.1) WRITE (6,*) 
     &      ' ### SQR OPERATOR RESULT      ###'
            IF (KPARM.EQ.2) WRITE (6,*) 
     &      ' ### EV  OPERATOR RESULT      ###'
          endif

          CALL PREN(KPARM,
     *    P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
          IPR=0
        ENDIF
      ENDIF

      IF (JPARM.EQ.3) GOTO 25
      IF (KPARM.EQ.IPARM(JPARM+1)) GOTO 25
      JPARM=JPARM+1
      KPARM=IPARM(JPARM)
      BIAS=5.0D0

       CALL ENERGY_RELSCF
     *    (P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)

      IF (IPR_RELSCF.GE.2) THEN
         write(LUPRI,*) ' total energy', energ
      ENDIF
      GOTO 15

   25 CONTINUE
      CALL ENERGY_RELSCF
     *    (P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *      ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)

      IF(IEX) 85,85,102

      IF (IPR_RELSCF.GE.2) THEN
        WRITE(LUPRI,'(/,2x,a,/)') 'going to CALL HAMIL'
      ENDIF

   30 CONTINUE
      CALL HAMIL
     *(P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
C
C     SAVE THE RESULTS OF THE PRECEDING SCF RUNS
      NSTEP2=0
      DO I=1,NSYM
        NBAS1=NBAS(I)
        NSH=NCSH(I)+NOSH(I)
        DO J=1,NSH
          J1=J+NSTEP2
          DO M=1,NBAS1
            CM2SCF(M,J1)=CM1SCF(M,J1)
CMI         ... index J1 in C out of boundaries !
            CM1SCF(M,J1)=C(M,J1)
          END DO
        END DO
        NSTEP2=NSTEP2+NSH
      END DO

      IF (NSVD.EQ.1) GO TO 2

      NSTEP =0
      NSTEP2=0
      DO I=1,NSYM
        NBAS1=NBAS(I)
            NOSHIC=NOSH(I)
        NSH=NCSH(I)+NOSHIC
        DO J=1,NSH
          N0=J+NSTEP2
          IF (J.GE.NSH) GO TO 60
   55     CALL EIGEA(FC,
     *    P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *        CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *        C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)

           if (IPR_RELSCF.GE.7) write (LUPRI,*) 'after eigea 55'

           GO TO 70
   60     IF (NOSHIC.LE.0) GO TO 55
   65     CALL EIGEA(FO,
     *    P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *        CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *        C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
   70     IF (NALARM.GT.0) GO TO 80
        END DO
        NSTEP =NSTEP +N1(I)
        NSTEP2=NSTEP2+NSH
      ENDDO ! of I=1,NSYM
      GO TO 3
C
C****   JACOBI DIAGONALIZATION
C
    2  CONTINUE
      IF (IPR_RELSCF.GE.5)  THEN
        write(LUPRI,'(/,2x,a)') 'call JDAAA'
      ENDIF

       CALL JDAAA
     *    (P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)

3     CONTINUE
C
C     AVERAGING / EXTRAPOLATION
C
      IF (IPR_RELSCF.GE.5)  THEN
        write(LUPRI,'(/,2x,a)') 'call EXTRAP'
      ENDIF

      CALL EXTRAP(IAV,REL,
     * P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     * CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     * C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)

      IF (NALARM.LE.0) GO TO 20
C
C     IF SELF CONSISTENCE DIVERGENCE OR TOO MANY ITERATIONS,
C     USE DENSITY OF LAST BUT ONE ITERATION AND PROCEED WITH NEXT
C     EXPONENT
C
      IF (NRUN.GT.2 .AND. (NALARM.EQ.3 .OR. NALARM.EQ.4)) THEN

      IF (NALARM.EQ.3)
     *WRITE (LUPRI,*) ' SELF CONSISTENCE DIVERGENCE - NITSCF=',NITSCF
      IF (NALARM.EQ.4)
     *WRITE (LUPRI,*) ' TOO MANY ITERATIONS - NITSCF=',NITSCF
      ZETA(NBVAR1)=ZSAV(2)
      WRITE (LUPRI,*) ZSAV(2),ESAV(2)
      IDEBUG=1
C
C     RESTORE DENSITY
C
      NSTEP2=0
      DO I=1,NSYM
        NBAS1=NBAS(I)
        NSH=NCSH(I)+NOSH(I)
        DO J=1,NSH
          J1=J+NSTEP2
          DO M=1,NBAS1
            C(M,J1)=CSAV2(M,J1)
            CM1SCF(M,J1)=CSAV2(M,J1)
          END DO
        END DO
        NSTEP2=NSTEP2+NSH
      END DO
      BIAS=5.
      NRUN=1
      NALARM=0
      GOTO 15

      ENDIF

   80  CONTINUE
      IF (IPR_RELSCF.GE.5) THEN
        WRITE(LUPRI,'(/,2x,a,/)') 'Going to CALL ALARM'
      ENDIF

        CALL AMFIALARM
     *(P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
      GO TO 251

   85 IF (NRUN.GT.1) GO TO 130
C
C     OPEN SUBROUTINE EXPONENT VARIATION
C
90    CONTINUE
      IF(MXVAR-MVAR-1) 93,92,91
   91 IF (NVAR.GT.NUMVAR) GO TO 95
   94 MVAR=MVAR+1
      NUMVAR=0
      GO TO 90
   92 IF (NUMVAR.NE.0) GO TO 91
   96 IF (NEXTRA.LE.0) GO TO 91
   97 CALL ZEXTRA
     *(P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
      GO TO 91
   93 IF (NEXTRA.LE.0) GO TO 251
      IEX=1
  102 CALL ZEXTRA
     *(P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
      IF (IEX.GT.0) GO TO 15
   99 IF (NEXTRA.LE.0) GO TO 251
  101 MVAR=0
      GO TO 90
   95 NSAFE=0
      BIAS=5.D0
      NUMVAR=NUMVAR+1
      NBVAR1=NBVAR(NUMVAR)
      SIGN=+1.D0
      ESAV(1)=ENERG
      NRUN1=NRUN
      NRUN=2
      A=ZETINC(NUMVAR)
      TERM =ZETA(NBVAR1)
      ZETA(NBVAR1)=TERM+A
      ZSAV(1)=TERM
  109 NSTEP2=0
      DO I=1,NSYM
        NBAS1=NBAS(I)
        NSH=NCSH(I)+NOSH(I)
        DO J=1,NSH
          J1=J+NSTEP2
          DO M=1,NBAS1
            CSAV2(M,J1)=CSAV1(M,J1)
            CSAV1(M,J1)=C(M,J1)
            IF (NRUN1.LE.1) GO TO 125
  122       TERM=(CSAV2(M,J1)-CSAV1(M,J1))/(ZSAV(2)-ZSAV(1))
            C(M,J1)=CSAV1(M,J1)+(ZETA(NBVAR1)-ZSAV(1))*TERM
  125       CONTINUE
          END DO
        END DO
        NSTEP2=NSTEP2+NSH
      END DO
      GO TO 200
C
  130 CONTINUE
      RED=ABS((ENERG-ESAV(1))/ESAV(1))
      IF (RED.GE.RMNED) GO TO 140

  135 if (IPR_RELSCF.GE.2) WRITE(LUPRI,136) RED,RMNED

  136 FORMAT(6X,'RELATIVE ENERGY DIFFERENCE',E12.2,1X,
     *       'SMALLER THAN RMNED',E12.2)
      NRUN=1
      GO TO 90
  140 IF (RED.LE.RMXED) GO TO 150
  145 NALARM=6
      GO TO 80
  150 IF (NRUN.GT.2) GO TO 170
  155 IF (ENERG.LE.ESAV(1)) GO TO 165
  160 IF (NSAFE.GT.0) GO TO 175
C
  162 SIGN=-1.D0
      A=-2.D0*ZETINC(NUMVAR)
      TERM =ZETA(NBVAR1)
      ZETA(NBVAR1)=TERM+A
      ZSAV(2)=TERM
      ESAV(2)=ENERG
      NSAFE=1
      GO TO 200
C
  165 ZSAV(2)=ZSAV(1)
      ESAV(2)=ESAV(1)
      ESAV(1)=ENERG
      NRUN=3
      NRUN1=2
      A=SIGN*ZETINC(NUMVAR)
      TERM =ZETA(NBVAR1)
      ZETA(NBVAR1)=TERM+A
      ZSAV(1)=TERM
      GO TO 109
C
  170 IF (ENERG.LT.ESAV(1)) GO TO 165
  175 ZSAV(3)=ZSAV(2)
      ZSAV(2)=ZSAV(1)
      TERM =ZETA(NBVAR1)
      GO TO 13
   13 ZSAV(1)= TERM
      DZ3=ZSAV(3)-TERM
      DZ2=ZSAV(2)-TERM
      ESAV(3)=ESAV(2)
      ESAV(2)=ESAV(1)
      ESAV(1)=ENERG
      DE3=ESAV(3)-ENERG
      DE2=ESAV(2)-ENERG
      NRUN=1

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,180)
  180 FORMAT(6X,'OPTIMIZED ONE DIMENSION')

      ZUMER=DE2*(DZ3**2)-DE3*(DZ2**2)
      ZUMER=ZUMER/(2.D0*(DE2*DZ3-DE3*DZ2))+TERM
      ZETA(NBVAR1)=ZUMER
      GO TO 109
C
  200 IF (ZETA(NBVAR1).LT.EXPMN) GO TO 205
      IF (NBVAR1.GT.NSAVF(1)) GO TO 214
      N11=1
      N12=NSAVF(1)
      GO TO 213
  214 N12=NSAVF(1)+NSAVF(2)
      IF (NBVAR1.GT.N12) GO TO 216
      N11=NSAVF(1)+1
      GO TO 213
  216 N11=N12+1
      N12=N12+NSAVF(3)
  213 DO 225 I=N11,N12
      IF (I-NBVAR1) 217,225,217
  217 IF (PQN(I)-PQN(NBVAR1)) 225,215,225
  215 ZETDIF=ABS(ZETA(I)-ZETA(NBVAR1))
      IF (ZETDIF-DMNXP) 220,225,225
  220 NALARM=8
      IIII=I
      GO TO 80
  225 CONTINUE
      GO TO 15
  205 NALARM=7
      GO TO 80

251   CONTINUE
      CALL OUTPUT_RELSCF
     *(P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN,
     *    WORK(KFREE),LFREE)


      IF (IPR_RELSCF.GE.5) THEN
        WRITE(LUPRI,*) 'after sizeof(S(:1)   =',size(S(:1))
        WRITE(LUPRI,*) 'after sizeof(S(:LPQ))=',size(S(:LPQ))
      ENDIF
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INPAR(NR,N ,A,B,C,X,W,NP,XP,WP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Array dim. fixes
C---------------------
C X(2) -> X(*)
C W(2) -> W(*)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
CMI   DIMENSION N (2),A(2),B(2),C(2),X(2),W(2),NP(2),XP(2),WP(2),
      DIMENSION N (*),A(*),B(*),C(*),X(*),W(*),NP(*),XP(*),WP(*),
     1          X1(150),W1(150),X2(12),W2(12)
      DATA DMAX/1.D60/
      IF (NR.EQ.1) GO TO 40
      IF (NP(1).NE.0) GO TO 20
      DO I = 2, NR
        IF (NP(I).EQ.0) GO TO 10
        IFLAG = 1
        GO TO 7000
   10   CONTINUE
      END DO
      GO TO 40
   20 DO 30 I = 2, NR
      IF (NP(I).NE.0) GO TO 30
      IFLAG = 1
      GO TO 7000
   30 CONTINUE
   40 W1(1) = 0.D0
      X1(1) = -1.D0
      NX = 1
      NXP = 1
      DO I = 1, NR
        N1 = N(I)
        N2 = NP(I)
        IF (C(I).GE.B(I).OR.C(I).LE.A(I)) GO TO 6900
        CALL GAUSS(N1,X1(2),W1(2))
        X1(N1+2) = 1.D0
        IF (N2.NE.0) CALL GAUSS(N2,X2,W2)
        IF (A(I).LT.B(I)) GO TO 100
        IFLAG = 2
        GO TO 7000
  100   IF (B(I).NE.DMAX) GO TO 3000
        IF (A(I).NE.-DMAX) GO TO 2000
        IFLAG = 3
        GO TO 7000
 2000   T1 = (C(I) - A(I)) + (C(I) - A(I))
        DO J = 1, N1
          X(NX) = C(I) + T1 * X1(J+1) / (1.D0 - X1(J+1))
          W(NX) = T1 * W1(J+1) / ((1.D0 - X1(J+1)) * (1.D0 - X1(J+1)))
          NX = NX + 1
        END DO
        IF (N2.EQ.0) GO TO 1000
        X(NX) = B(I)
        W(NX) = 0.D0
        NX = NX + 1
        DO J = 1, N1
          C1 = (X1(J+1) + X1(J)) * .5D0
          C2 = (X1(J+1) - X1(J)) * .5D0
          T2 = C2 * T1
          DO K = 1, N2
            U = C1 + C2 * X2(K)
            XP(NXP) = C(I) + T1 * U / (1.D0 - U)
            WP(NXP) = T2 * W2(K) / ((1.D0 - U) * (1.D0 - U))
            NXP = NXP + 1
          END DO
        END DO
        GO TO 1000
 3000   IF (A(I).NE.-DMAX) GO TO 4000
        T1 = (B(I) - C(I)) + (B(I) - C(I))
        DO J = 1, N1
          X(NX) = C(I) + T1 * X1(J+1) / (1.D0 + X1(J+1))
          W(NX) = T1 * W1(J+1) / ((1.D0 + X1(J+1)) * (1.D0 + X1(J+1)))
          NX = NX + 1
        END DO
        IF (N2.EQ.0) GO TO 1000
        X(NX) = B(I)
        W(NX) = 0.D0
        NX = NX + 1
        DO J = 1, N1
          C1 = (X1(J+1) + X1(J)) * .5D0
          C2 = (X1(J+1) - X1(J)) * .5D0
          T2 = C2 * T1
          DO K = 1, N2
            U = C1 + C2 * X2(K)
            XP(NXP) = C(I) + T1 * U / (1.D0 + U)
            WP(NXP) = T2 * W2(K) / ((1.D0 + U) * (1.D0 + U))
            NXP = NXP + 1
          END DO
        END DO
        GO TO 1000
 4000   BETA = (B(I) + A(I) - (C(I) + C(I))) / (B(I) - A(I))
        IF (ABS(BETA).LE..1D0) GO TO 5000
        T1 = (B(I) - A(I)) * (1.D0 - BETA * BETA) * .5D0
        D2 = T1 / BETA
        D1 = ((A(I) + B(I)) + (A(I) - B(I)) / BETA) * .5D0
        DO J = 1, N1
CMI      ....   for array 'x', upper bound of dimension 1 exceeded 
          X(NX) = D1 + D2 / (1.D0 - BETA * X1(J+1))
CMI      ....   for array 'w', upper bound of dimension 1 exceeded
          W(NX) = T1 * W1(J+1) / ((1.D0 - BETA * X1(J+1))
     1                              * (1.D0 - BETA * X1(J+1)))
          NX = NX + 1
        END DO
        IF (N2.EQ.0) GO TO 1000
        X(NX) = B(I)
        W(NX) = 0.D0
        NX = NX + 1
        DO J = 1, N1
          C1 = (X1(J+1) + X1(J)) * .5D0
          C2 = (X1(J+1) - X1(J)) * .5D0
          T2 = C2 * T1
          DO K = 1, N2
            U = C1 + C2 * X2(K)
            XP(NXP) = D1 + D2 / (1.D0 - BETA * U)
            WP(NXP) = T2 * W2(K) / ((1.D0 - BETA * U) * 
     &                (1.D0 - BETA * U))
            NXP = NXP + 1
          END DO
        END DO
        GO TO 1000
 5000   D1 = (B(I) + A(I)) * .5D0
        D2 = (B(I) - A(I)) * .5D0
        DO J = 1, N1
          X(NX) = D1 + D2 * X1(J+1)
          W(NX) = D2 * W1(J+1)
          NX = NX + 1
        END DO
        IF (N2.EQ.0) GO TO 1000
        X(NX) = B(I)
        W(NX) = 0.D0
        NX = NX + 1
        DO J = 1, N1
          C1 = (X1(J+1) + X1(J)) * .5D0
          C2 = (X1(J+1) - X1(J)) * .5D0
          T2 = C2 * D2
          DO K = 1, N2
            U = C1 + C2 * X2(K)
            XP(NXP) = D1 + D2 * U
            WP(NXP) = T2 * W2(K)
            NXP = NXP + 1
          END DO
        END DO
 1000   CONTINUE
      END DO
      RETURN
 6900 IFLAG = 5
 7000 WRITE (6,1) IFLAG
    1 FORMAT('0******  ERROR MESSAGE   INPAR/',I3)
      CALL CHECK(0,0,0)
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INPUT_AT34(C,CTRAN,PQN,ZETA,ZETINC,ITRAN,NROW,NBVAR,
     &                 IORCN,CREAD,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Set up various parameters for relscf.
C
C   Called from AT3
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "scfarr.h"
#include "relscf_dim.h"
#include "relscf_v10.h"
C***********************************************************************
C     ARRAYS SIZES:
C     C: NOR*MXBF*NSYM AND  CTRAN,ITRAN: MXBAS*MXOCC     *
C***********************************************************************
CMI   DIMENSION C(LDIMC),CTRAN(LDIMC),PQN(MXB),ZETA(MXB),ZETINC(MXB),
CMI   DIMENSION C(0:LDIMC),CTRAN(LDIMC),PQN(MXB),ZETA(MXB),ZETINC(MXB),
CMI   DIMENSION C(*),CTRAN(LDIMC),PQN(MXB),ZETA(MXB),ZETINC(MXB),
      DIMENSION C(0:*),CTRAN(LDIMC),PQN(MXB),ZETA(MXB),ZETINC(MXB),
C     DIMENSION C(*),CTRAN(LDIMC),PQN(MXB),ZETA(MXB),ZETINC(MXB),
     *          ITRAN(LDIMC),NROW(MXB),NBVAR(MXB),IORCN(MXB),
     *          CREAD(MXB*MXO),WORK(LWORK)
C
      COMMON /CHAR/ NAME(12)
      CHARACTER*4 NAME,TITLE
C
      COMMON/WIGXXX/A1,A2,EXP1,EXP2,UP(MXB),DOWN(MXB),IUP(4),NGAUSS
C
      DIMENSION NOOC(5),VCC(24),TITLE(12)
      CHARACTER*8 OPER(3)
      EQUIVALENCE (AJMN(1),VCC(1)),(TITLE(1),NAME(1))
#include "amfi_if.h"
c     logical relscf_verbose
C
C     SET DEFAULT VALUES
C
      CHARACTER*4 NBLANK,BSTATE,STATE
      LOGICAL EX, RELSCF_RST, FNDLAB
      DATA OPER /'NON-REL ','SQR     ','EV      '/

CMI   DATA NBLANK /'    '/,BSTATE/'  '/,IBLANK /4H    /
      DATA NBLANK /'    '/,BSTATE/'  '/
CMI  ... gfortran gives warning, I don't know how to fix this ...
CMI   DATA IBLANK /4H    /

      DATA NCOUNT /0/

CMI   write(LUPRI,*) 'INPUT_AT34: NBLANK=',NBLANK
CMI   write(LUPRI,*) 'INPUT_AT34: IBLANK=',IBLANK
C
c     relscf_verbose = .false.
C
#include "memint.h"
      NDIMC=MXOCC*MXOCC
      NDIMZ1=MXBAS
      NDIMR=MXBAS
      NDIMT=MXOCC*MXBAS
      IF(NCOUNT.GT.0) GO TO 3
      DO I=1,4
        ISIZE(I)=0
      END DO
      NDIAG=10
      NXTRP=10

CMI  ... maximum number of iteration ...
CMI   MXXTRP=80
      MXXTRP=MXXTRP_INP

      DGATH=1.D-4
      SCFAT=1.D-3
      EXPMN=1.1D-4
      DMNXP=1.1D-3
      RMNED=1.1D-8
      RMXED=0.99D-1
      IPARM(1)=0
      IPARM(2)=2
      IPARM(3)=2
C
C     EXTREME EXPONENT FOR S,P,D,F
C
      TRLO(1)=1.D-6
      TRLO(2)=1.D-6
      TRLO(3)=1.D-6
      TRLO(4)=1.D-6
      TRHI(1)=1.D+8
      TRHI(2)=1.D+6
      TRHI(3)=1.D+5
      TRHI(4)=1.D+4

      MX100=100
      NVAR=0
      MXVAR=0
      NEXTRA=0
      NORBIT=0
      CON1=0.028816D0
      CON2=1.1D0
      NGAUSS=64
      A1=0.04412D0
      A2=0.0446D0
      EXP1=0.66666666666667D0
      EXP2=0.16666666666667D0
      NBLOCK=2
      NSIZE(1)=99
      NSIZE(2)=99
      RADINC(1)=0.1D0
      RADINC(2)=0.2D0
      NPRINT=0
      IPRINT=0
      NPVEFA=0
      NPOEFA=0
      NPENFA=0
      NPBAFI=0
      NPCOFI=0
      IDATA=2
      NDIMPQ=NDIMPQPAR
    3 DO 10 I=1,12
   10 NAME(I)=NBLANK

CMI ... fill/initialize with zeroes
      DO I=1,4
        IUP(I)=0
        OCCUP(I)=0.0D0
        NBAS(I)=0
        NOOC(I)=0
        NOSH(I)=0
        NCSH(I)=0
      END DO
      STATE=BSTATE
      DO I=1,24
        AJMN(I)=0.0D0
      END DO

CMI   write(LUPRI,*) 'IBLANK=',IBLANK
      DO I=1,NDIMR
CMI     IORCN(I)=IBLANK
        IORCN(I)=0
CMI     write(LUPRI,*) 'I,IORCN(I)',I,IORCN(I)
      END DO
C
      DO I=1,NDIMC
        C(I)=0.0D0
      END DO
      DO I=1,NDIMZ1
        ZETINC(I)=0.0D0
        NBVAR(I)=I
      END DO
      DO I=1,NDIMR
        NROW(I)=1
      END DO
      DO I=1,NDIMT
        ITRAN(I)=0
        CTRAN(I)=1.0D0
      END DO
      LPQ=0
      MXCC=1
      BIAS=5.0D0
      IGH=0
      GJ0=0.D0
      GK0=0.D0
      ZETRAT=0.05D0
      NEXPND=0
      IVEC=0
C
2080  FORMAT(12A4)

! ... reading from a very simple input file !
      LUTOSCF = 9
      OPEN(UNIT=LUTOSCF,FILE='TOSCF',STATUS='OLD',FORM='FORMATTED')
      READ(LUTOSCF,*,END=2002) charge
      READ(LUTOSCF,*) (NBAS(I),I=1,4)
      READ(LUTOSCF,*) (NCSH(I),I=1,4)
      READ(LUTOSCF,*) (NOOC(I),I=1,4)
      READ(LUTOSCF,*) IREL

      IF (IPR_RELSCF.GE.2) THEN
        WRITE(LUPRI,'(/,2X,A)')
     &  'INPUT_AT34: read data from the TOSCF file'
CSK     charge is relative atomic mass
        WRITE(LUPRI,'(2X,A,F12.6)') 'charge=',charge
        WRITE(LUPRI,'(2X,A,4I3)')   'NBAS:',(NBAS(I),I=1,4) 
        WRITE(LUPRI,'(2X,A,4I3)')   'NCSH:',(NCSH(I),I=1,4) 
        WRITE(LUPRI,'(2X,A,4I3)')   'NOOC:',(NOOC(I),I=1,4) 
        WRITE(LUPRI,'(2X,A,I3/)')   'parameter IREL=',IREL
      ENDIF

      if (IREL.eq.0) then
        Iparm(1)=0
        Iparm(2)=0
        Iparm(3)=0
      else
        Iparm(1)=0
        Iparm(2)=2
        Iparm(3)=2
      endif

CMI  ... reading s,p,d,f exponents from the AMFI prepared file...      
      nexp=NBAS(1)+NBAS(2)+NBAS(3)+NBAS(4)
CMI   read(5,*)(ZETA(I),I=1,nexp)
      READ(LUTOSCF,*)(ZETA(I),I=1,nexp)
      CLOSE(LUTOSCF,STATUS="KEEP")

      IF (IPR_RELSCF.GE.2) THEN
        WRITE(LUPRI,'(/,2X,A)') 'READ EXPONENTS:'
        WRITE(LUPRI,*) (ZETA(I),I=1,nexp)
        WRITE(LUPRI,*)
      ENDIF

      NVAR=0
      IVEC=1
      NPVEFA=1
      BIAS=7.0
      icount=0
      do jrun=1,4
      do irun=1,nbas(jrun)
      icount=icount+1
      PQN(icount)=float(jrun)
      enddo
      enddo
CBS   simplify INPUT
      IF (IPARM(2).EQ.0) THEN
      DO I=1,4
        TRHI(I)=0.D0
        TRLO(I)=0.D0
      END DO
      ENDIF
      IF (IGH.GT.1) GOTO 461
      IF (IGH.EQ.1) GOTO 460
      GJ0=1.D0
      GK0=1.D0
      GOTO 461
460   GK0=1.D0-GJ0
      GOTO 461
461   CONTINUE
      GJ00=GJ0
      GK00=GK0
C
      NCOUNT=NCOUNT+1

C ...  set up NSYM value
 777  NSYM=0
      DO I=1,4
CMI     IF(NCSH(I)+NOOC(I).EQ.0) GO TO 2
        IF(NBAS(I).EQ.0) GO TO 2
        NSYM=NSYM+1
        IF(NOOC(I).EQ.0) GO TO 2
        NOSH(I)=1
        if (occup(i).eq.0.d0) OCCUP(I)=NOOC(I)
 2      CONTINUE
      END DO

      NSVD=1
      NCOVAR=0
      ZN=CHARGE
C
      IF (IPR_RELSCF.GE.1) THEN

      WRITE(LUPRI,80) (NAME(I),I=1,12)
   80 FORMAT ('1',4X,'ATOMIC HARTREE-FOCK PROGRAM '/5X,
     1'MODIFIED BY LAURA GIANOLIO IST. DONEGANI'/5X,
     2'MONTEDISON SPA NOVARA (ITALY). MARCH 1,1978'/5X,
     4'DYNAMIC STORAGE BUILT IN BY BERND HESS'/5X,
     5'UNIVERSITY OF BONN - AUGUST 29,1986'/5X,
     4'IBM+ AND CRAY VERSION BY BERND HESS'/5X,
     5'UNIVERSITY OF BONN - OCTOBER 26,1986'/,5X,
     5'T.SAUE,M.ILIAS - 2007,2008         '/,5X,
     6//5X,12A4)
     
      WRITE (LUPRI,13)
  13  FORMAT ('0',5X,19HTYPE GAUSS ORBITALS)
      IF(NSVD.EQ.0)WRITE(LUPRI,504)
      IF(NSVD.EQ.1)WRITE(LUPRI,505)
  504 FORMAT('0',5X,'SINGLE VECTOR DIAGONALIZATION')
  505 FORMAT ('0',5X,'JACOBI DIAGONALIZATION')

  121 FORMAT('0',5X,'CUT OFF PARAMETERS =',F10.6,3X,F10.6)

      WRITE(LUPRI,85) ZN
   85 FORMAT('0',5X,'CHARGE =',F10.6)

      WRITE(LUPRI,90) BIAS,DGATH
   90 FORMAT('0',5X,'DIAGONALIZATION SCF THRESHOLD BIAS  =', E10.3,8X,
     1'DIAGONALIZATION DIVERGENCE =',7X,E10.3)
      
      WRITE(LUPRI,95) NXTRP,SCFAT
   95 FORMAT(6X,'NUMBER OF SCF ITERATIONS            =',I3,16X,
     *'SCF DIVERGENCE             =',7X,E10.3)

      WRITE(LUPRI,100) NDIAG,EXPMN
  100 FORMAT(6X,'NUMBER OF DIAGONALIZATION ITERATIONS=',I3,16X,
     *'MINIMUM EXPONENT           =',7X,E10.3)

      WRITE(LUPRI,105) MXXTRP,DMNXP
  105 FORMAT(6X,'MAXIMUM NUMBER OF ITERATIONS        =',I3,16X,
     *'MINIMUM EXPONENT DIFFERENCE=',7X,E10.3)

      WRITE(LUPRI,106) RMNED,RMXED
  106 FORMAT(6X,'RELATIVE MINIMUM ENERGY DIFFERENCE  =', E10.3,
     *8X,'RELATIVE MAXIMUM ENERGY DIFFERENCE=', E10.3)
      ENDIF

      IF (IPR_RELSCF.GE.0) THEN

        WRITE(LUPRI,115) (NBAS(I),I=1,NSYM)
  115   FORMAT (/,6X,'SYMMETRY SPECIES',12X,'S',7X,'P',7X,'D',7X,'F'/6X,
     &          'NUMBER OF BASIS FUNCTIONS:',4(I3,5X))

        WRITE(LUPRI,120) (NCSH(I),I=1,NSYM)
  120   FORMAT(6X,'NUMBER OF CLOSED SHELLS  :',10(I3,5X))

        WRITE(LUPRI,125 ) (NOOC(I) ,I=1,NSYM)
  125   FORMAT(6X,'OPEN SHELL OCCUPATION    :',10(I3,5X),/)

      ENDIF

C --- verify number of basis functions wrt shell occupation
      DO I=1,NSYM
       IF ((NCSH(I)+NOSH(I)).GT.NBAS(I)) THEN
         WRITE(LUPRI,'(/,2x,a,i2,a,i2,a,i2)')  
     &   'symmetry:',I,' basis fu:',NBAS(I),
     &   ' occupied:',NCSH(I)+NOOC(I)
         CALL QUIT('RELSCF: insufficient number of basis functions !')
       ENDIF
      ENDDO

      IF (IPR_RELSCF.GE.2) THEN
       WRITE(LUPRI,135)
  135  FORMAT (/,'0',5X,'VECTOR COUPLING COEFFICIENTS K')

       WRITE(LUPRI,140) (AJMN(I),I=1,24)
  140 FORMAT(5D21.8)
      ENDIF 

      DO i=1,24
        ajmn(i)=0.5*ajmn(i)
      END DO
      NBAST=0
      NOSHT=0
      NCSHT=0
      DO I=1,NSYM
        LPQ=LPQ+NBAS(I)*(NBAS(I)+1)/2
        NSAVF(I)=NBAS(I)
        NOSHT=NOSHT+NOSH(I)
        NCSHT=NCSHT+NCSH(I)
        NBAST=NBAST+NBAS(I)
      END DO
      NSHT=NCSHT+NOSHT

      DO I=1,NBAST
        IF (PQN(I).LE.0.) GO TO 370
        IF (ZETA(I).LT.EXPMN) GO TO 370
      END DO

      I=-1
      IF ((BIAS.LT.5.).OR.(BIAS.GT.10.)) GO TO 370
      I=I-1
      IF ((NDIAG.LT.5).OR.(NDIAG.GT.25)) GO TO 370
      I=I-1
      IF ((NXTRP.LT.5).OR.(NXTRP.GT.25)) GO TO 370
      I=I-1
      IF ((DGATH.LT.1.D-5).OR.(DGATH.GT.1.D-2)) GO TO 370
      I=I-1
      IF ((SCFAT.LT.1.D-4).OR.(SCFAT.GT.0.1)) GO TO 370
      I=I-1
      IF ((EXPMN.LT.1.D-4).OR.(EXPMN.GT.1.)) GO TO 370
      I=I-1
      IF ((DMNXP.LT.1.D-3).OR.(DMNXP.GT.1.)) GO TO 370
      I=I-1
      IF ((RMNED.LT.1.D-8).OR.(RMNED.GT.1.D-2)) GO TO 370
      I=I-1
      IF ((RMXED.LT.1.D-5).OR.(RMXED.GT.0.1)) GO TO 370
C
      if (IPR_RELSCF.GE.2) WRITE (LUPRI,145)
  145 FORMAT ('0',5X,'BASIS FUNCTIONS (PRINCIPAL QUANTUM NUMBER, ORBITAL
     1 EXPONENTS)'//15X,'S',28X,'P',28X,'D',28X,'F'/ )

  160 FORMAT(4(5X,F4.1,5X,F15.6))
C
      NBAS1=NBAS(1)
      IF (NSYM.GT.1) GO TO 155
       if (IPR_RELSCF.GE.2) then
          DO   I=1,NBAS1
             WRITE(LUPRI,160)  PQN(I),ZETA(I)
          END DO
       endif
       GO TO 215
C
  155 NBAS2=NBAS(2)
      IF (NSYM.GT.2) GO TO 156
      NBA3=1

  197 DO 185 I=NBA3,NBAS2
       J=NBAS1+I

  185 IF (IPR_RELSCF.GE.2) WRITE(LUPRI,160)
     &            PQN(I),ZETA(I),PQN(J),ZETA(J)

  196 IF (NBAS1.EQ.NBAS2) GO TO 215
       NBA2=NBAS2+1
       if (IPR_RELSCF.GE.2) then
          DO   I=NBA2,NBAS1
             WRITE(LUPRI,160) PQN(I),ZETA(I)
          END DO
       endif
       GO TO 215
C
  156 NBAS3=NBAS(3)
      IF(NSYM.GT.3) GO TO 157
      NBA4=1
  207 DO 195  I=NBA4,NBAS3
       J=NBAS1+I
       K=NBAS1+NBAS2+I

  195 if (IPR_RELSCF.GE.2) 
     &   WRITE(LUPRI,160) PQN(I),ZETA(I),PQN(J),ZETA(J),PQN(K),ZETA(K)

  206 IF (NBAS3.EQ.NBAS2) GO TO 196
      NBA3=NBAS3+1
      GO TO 197
C
  157 NBAS4=NBAS(4)
      if (IPR_RELSCF.GE.2) then
         DO I=1,NBAS4
            J=NBAS1+I
            K=NBAS1+NBAS2+I
            L=NBAS1+NBAS2+NBAS3+I
            WRITE(LUPRI,160) PQN(I),ZETA(I),PQN(J),
     &           ZETA(J),PQN(K),ZETA(K),PQN(L), ZETA(L)
         END DO
      endif
      IF(NBAS4.EQ.NBAS3) GO TO 206
      NBA4=NBAS4+1
      GO TO 207
C
  215 IF (NVAR.EQ.0) GO TO 1003
      IF(NVAR.LT.0) NVAR=NBAST
      IF(MXVAR.EQ.0) MXVAR=3
      IF(NEXTRA.EQ.0) NEXTRA=1
      DO I=1,NVAR
        J=NBVAR(I)
        IF(ZETINC(I).EQ.0.0D0) ZETINC(I)=ZETA(J)*ZETRAT
      END DO

 1003 NBASM=0
      DO I=1,NSYM
        NBAS1=NBAS(I)
        IF(NBASM.LT.NBAS1) NBASM=NBAS1
      END DO
      MXBF=NBASM
C
C     S. Knecht - Jan 2009: implementation of a working restart option; 
C     Atomic SCF calculations on neutral atoms are sometimes hard to 
C     converge. However, the ionic system may not. In this case, we want
C     to do the ion first and reuse the optimized coefficients for the
C     neutral atom.
C
C     restart from old coefficients for this atom?
C
C     initialize coefficients vector
      CALL DZERO(C,LDIMC+2)
CSK   IPR_RELSCF_SAVE = IPR_RELSCF
CSK   IPR_RELSCF   = 2
      LURELSCFCOEF = 7
C
C     search label RELSCFLABX (coefficients) on file RELSCF_COEF
      INQUIRE(FILE='RELSCF_COEF',EXIST=EX)
      RELSCF_RST = .FALSE.
      IF( EX )THEN
        CALL OPNFIL(LURELSCFCOEF,'RELSCF_COEF','OLD','INat34')
        REWIND LURELSCFCOEF
        IF (FNDLAB(RELSCFLABX,LURELSCFCOEF)) RELSCF_RST = .TRUE.
        CLOSE(LURELSCFCOEF,STATUS='KEEP')
      END IF
      IF(IPR_RELSCF.GE.2) WRITE(LUPRI,*) ' INPUT VECTORS'
      if(RELSCF_RST.and.IPR_RELSCF.GE.2) 
     &    write(LUPRI,*) 'starting orbitals read from file'
      IF( RELSCF_RST ) THEN 
        CALL OPNFIL(LURELSCFCOEF,'RELSCF_COEF','OLD','INat34')
        REWIND LURELSCFCOEF
        CALL IZERO(IRELSCF_INFO,MXRELSCF_INFO)
C       scratch memory
        CALL MEMGET('REAL',KRELCOEFFX,LDIMC+2,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KRELCOEFFX),LDIMC+2)
C
C       read info array
        CALL IREAKRMC(LURELSCFCOEF,RELSCFLABI,IRELSCF_INFO,
     &                MXRELSCF_INFO)
        IPLACEVLEN = 0
CSK     WRITE(LUPRI,*) ' IRELSCF_INFO(1)',IRELSCF_INFO(1)
        IPLACEVLEN = (IRELSCF_INFO(1) * 2) + 2
CSK     WRITE(LUPRI,*) ' IPLACEVLEN = ',IPLACEVLEN
        ISIZEVLEN  = IRELSCF_INFO(IPLACEVLEN)
C
C       read coefficients from file RELSCF_COEF to temporary array
        CALL REAAMFI(LURELSCFCOEF,RELSCFLABX,WORK(KRELCOEFFX),
     &               ISIZEVLEN)
CSK     WRITE(LUPRI,*) ' coeff matrix from file: len = ',ISIZEVLEN
CSK     CALL WRTMATMN(WORK(KRELCOEFFX),1,ISIZEVLEN,1,ISIZEVLEN,LUPRI)
C
C       FIXME: end of loop should be LMAX! - SK Jan 2009
        LBASL_MAX = 0
        DO L =1,4
          LBASL_MAX = MAX(LBASL_MAX,NBAS(L))
        END DO
        CALL MEMGET('REAL',KRELCLMAX,LBASL_MAX,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KRELCLMAX),LBASL_MAX)
        CLOSE(UNIT=LURELSCFCOEF,STATUS='KEEP')
      END IF
C
C     FIXME: end of loop should be LMAX! - SK Jan 2009
      IBEG  = 0
      IEND  = 0
      IENDL = 0
      DO L=1,4
        IF(NBAS(L).EQ.0) GO TO 665
        IF (IPR_RELSCF.GE.2) WRITE(LUPRI,*) ' SYMMETRY',L
        IF (IPR_RELSCF.GE.2) WRITE(LUPRI,*) ' NBAS(L)',NBAS(L)
        IF (IPR_RELSCF.GE.2) WRITE(LUPRI,*) ' MXBF',MXBF
        NTT=0
        IF(NOOC(L).GT.0) NTT=1
        NTTC=NCSH(L)+NTT
CSK     WRITE(LUPRI,*) ' NTTC, NCSH(L), NTT',NTTC, NCSH(L), NTT
        DO I=1,NTTC
          if (IPR_RELSCF.GE.2) WRITE(LUPRI,*) ' VECTOR',I
          IBEG=IEND+1
          IEND=IBEG+NBAS(L)-1
CSK       IF (IPR_RELSCF.GE.2) WRITE(LUPRI,*) 'IBEG, IEND',IBEG, IEND
          if (RELSCF_RST) THEN
            CALL DZERO(WORK(KRELCLMAX),LBASL_MAX)
            TEST_BASLNRM = 0.0D0
            do II=1,NBAS(L)
CSK           WRITE(LUPRI,*) ' index IENDL + I-1+(II-1)*NTTC',
CSK  &                               IENDL + I-1+(II-1)*NTTC
              C(II+IBEG-1)=WORK(KRELCOEFFX+ IENDL + I-1+(II-1)*NTTC)
              WORK(KRELCLMAX+II-1) = 
     &        WORK(KRELCOEFFX+IENDL+ I-1+(II-1)*NTTC)
            enddo
C
C           check whether this vector is a null vector, i.e. starting
C           from coefficients obtained for the ionized atom where this 
C           shell was empty.
            TEST_BASLNRM = DDOT(NBAS(L),WORK(KRELCLMAX),1,
     &                          WORK(KRELCLMAX),1)
            IF( TEST_BASLNRM .eq. 0.0D0 ) C(IBEG+I*NBAS(L)/NTTC-1)=1d0
          else 
            C(IBEG+I*NBAS(L)/NTTC-1)=1d0
          end if
          if (IPR_RELSCF.GE.2) WRITE(LUPRI,666) (C(II),II=IBEG,IEND)
Cbuggy? - SK - Jan 2009  IEND=IEND+MXBF-NBAS(L)
        END DO
 665    CONTINUE
        IENDL = IEND
      END DO
 666  FORMAT(4F18.14)
CSK   IPR_RELSCF = IPR_RELSCF_SAVE
      if (RELSCF_RST) CALL MEMREL('INat34',WORK,KWORK,KWORK,KFREE,LFREE)
C
 1005 IF (NVAR.EQ.0) GO TO 209
      if (IPR_RELSCF.GE.2)
     &   WRITE (6,216) (NBVAR(I),ZETA(NBVAR(I)),ZETINC(I),I=1,NVAR)
  216 FORMAT ('0',5X,'EXPONENT VARIATION DATA (EXPONENT INDICES,
     *EXPONENTS, INCREMENTS)',
     1//(6X,I2,3X,E12.6,3X,E12.6))
  209 DO 6 I=1,NDIMC
      IF(C(I).NE.0.0D0) GO TO 9
    6 CONTINUE
      GO TO 214
    9 if (IPR_RELSCF.GE.2) WRITE(LUPRI,210)
  210 FORMAT ('0',5X,'INPUT TRIAL VECTORS'/)
      KMIN=0
      KMAX=NBASM*(NSHT-1)+NBAS(NSYM)
      DO I=1,NBASM
        KMIN=KMIN+1
        if (IPR_RELSCF.GE.2) WRITE (6,240) (C(K),K=KMIN,KMAX,NBASM)
      END DO
  240 FORMAT(10(2X,F9.6))
  214 IF (NORBIT.EQ.0) GO TO 2000
      if (IPR_RELSCF.GE.2) WRITE (6,250)
  250 FORMAT ('0',5X,'NUMERICAL WAVE FUNCTIONS (POINT, INCREMENTS)'/)
      if (IPR_RELSCF.GE.2) WRITE (6,260) (NSIZE(I),RADINC(I),I=1,NBLOCK)
  260 FORMAT (' ',6X,I5,4X,F9.6)
C
 2000 NOR=NSHT
      NBF=NBAST
      LQM=LPQ*(LPQ+1)/2
      IF(NDIMPQ.NE.NDIMPQPAR) LQM=NDIMPQ
      J=0
       J=NGAUSS*NBLOCK
      GO TO 33
CMI ... ftnchek complains that "No path to this statement"
      IF(NORBIT.NE.1) GOTO 995
      J=0
      DO I=1,NBLOCK
        J=J+NSIZE(I)
      END DO
CMI  * end
   33 I=J*(NSHT+1)
      IF(I.GT.LQM) LQM=I
      IF(NDIMPQ.LT.I) NDIMPQ=I
C
C     WRITE OUT PARAMETERS OF CURRENT RUN
C

995   if (IPR_RELSCF.GE.2) THEN
         WRITE (LUPRI,996) MXBF,NBF,NOR,LPQ,LQM,GJ0,GK0
      endif
996   FORMAT(/' MXBF (LARGEST NO OF FCT. PER SYMM.) =',I10/,
     *        ' NBF  (NUMBER OF BASIS FUNCTIONS)    =',I10/,
     *        ' NOR  (NUMBER OF OCCUPIED ORBITALS)  =',I10/,
     *        ' LPQ  (MATRIX SIZE 1)                =',I10/,
     *        ' LQM  (MATRIX SIZE 2)                =',I10/,
     *        ' GJ0                                 =',F20.10/,
     *        ' GK0                                 =',F20.10/)

      IF (IPR_RELSCF.GE.2) WRITE (LUPRI,997) MX100,TRLO,TRHI
997   FORMAT(/' MX100                               =',I10/,
     *        ' LOW  EXPONENT THRESHOLDS            =',4(1X,D14.4)/,
     *        ' HIGH EXPONENT THRESHOLDS            =',4(1X,D14.4))
      DO I=2,3
        IF (IPARM(I).EQ.IPARM(I-1)) GOTO 993
      END DO
      I=4
993   I=I-1

      if (IPR_RELSCF.GE.2) 
     &   WRITE (LUPRI,992) (IPARM(J),OPER(IPARM(J)+1),J=1,I)
992   FORMAT(' FOLLOWING OPERATOR(S) WERE SPECIFIED BY IPARM:'/,
     *3(3X,I3,'=',A8))
      IF (NBF.GT.MXBAS) GOTO 999
      IF (NOR.GT.MXOCC) GOTO 999

      IF (NOR*MXBF*NSYM.GT.LDIMC) THEN
       WRITE (LUPRI,970) LDIMC
970    FORMAT(' #### DIMENSION OVERFLOW: ARRAY C #####'/,
     *' NOR*MXBF*NSYM > ',I6/)
        CALL FLSHFO(LUPRI)
       CALL QUIT( 'DIMENSION OVERFLOW FOR ARRAY C !!!')
      ENDIF
C
C     INPUTS CALLS AT3DIM, WHICH PASSES ARRAY ADDRESSES TO
C     OTHER SUBROUINES
C
      CALL AT3DIM(NCOVAR,C,CTRAN,PQN,ZETA,ZETINC,ITRAN,
     &            NROW,NBVAR,IORCN,WORK,LWORK)
      RETURN
C
C     ERROR MESSAGE (DIMENSION OVERFLOW)
C
999   WRITE (6,998) NBF,MXBAS,NOR,MXOCC
998   FORMAT(//' *** PROGRAM DIMENSION OVERFLOW ***'/,
     *' NBF      ',I10,' (MAX IS ',I10,' )'/,
     *' NOR      ',I10,' (MAX IS ',I10,' )'//,
     *' INCREASE DIMENSIONS IN FORTRAN SOURCE ***')
      RETURN
C
C     END OF DATA ON INPUT FILE
C
 2002 STOP
C
C     INPUT ERRORS
C
 370  WRITE (6,375) I
 375  FORMAT ('0',30X,'ERROR DATA',I10)
      STOP 'INPUTERR'
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE JACOB(A,B,EIG,N,RANGE,IC)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     IC=1 EIGNEVALUES AND VECTORS REARRANGED
C       =0 LEFT AS THEY ARE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(N,N),B(N,N),EIG(N)
      ENUI=0.D0
      U1=N
      DO I=1,N
        B(I,I)=1.D0
        EIG(I)=A(I,I)
        DO J=1,I
          IF (I.NE.J) THEN
            B(I,J)=0.D0
            B(J,I)=0.D0
            ENUI=ENUI+A(I,J)*A(I,J)
          END IF
        END DO
      END DO
      IF(ENUI)200,200,10
   10 ENUI=DSQRT(2.D0*ENUI)
      ENUF=ENUI*RANGE/U1
      IND=0
      THR=ENUI
   15 THR=THR/U1
   20 L=1
   25 M=L+1
   30 IF(DABS(A(M,L))-THR)90,35,35
   35 IND=1
      X=.5D0*(EIG(L)-EIG(M))
      Y=-A(M,L)/DSQRT(A(M,L)*A(M,L)+X*X)
      IF(X)40,45,45
   40 Y=-Y
   45 IF(Y.GT.1.D0)Y=1.D0
      IF(Y.LT.-1.D0)Y=-1.D0
      XY=1.D0-Y*Y
      SINT=Y/DSQRT(2.D0*(1.D0+DSQRT(XY)))
      SINT2=SINT*SINT
      COST2=1.D0-SINT2
      COST=DSQRT(COST2)
      SINCS=SINT*COST
      DO I=1,N
        IF(I-M)50,80,55
   50   IM=M
        MM=I
        GO TO 60
   55   IM=I
        MM=M
   60   IF(I-L)65,80,70
   65   IL=L
        LL=I
        GO TO 75
   70   IL=I
        LL=L
   75   X=A(IL,LL)*COST-A(IM,MM)*SINT
        A(IM,MM)=A(IL,LL)*SINT+A(IM,MM)*COST
        A(IL,LL)=X
   80   X=B(I,L)*COST-B(I,M)*SINT
        B(I,M)=B(I,L)*SINT+B(I,M)*COST
        B(I,L)=X
      END DO
      X=2.D0*A(M,L)*SINCS
      Y=EIG(L)*COST2+EIG(M)*SINT2-X
      X=EIG(L)*SINT2+EIG(M)*COST2+X
      A(M,L)=(EIG(L)-EIG(M))*SINCS+A(M,L)*(COST2-SINT2)
      EIG(L)=Y
      EIG(M)=X
   90 IF(M-N)95,100,95
   95 M=M+1
      GO TO 30
  100 IF(L-M+1)105,110,105
  105 L=L+1
      GO TO 25
  110 IF(IND-1)120,115,120
  115 IND=0
      GO TO 20
  120 IF(THR-ENUF)200,200,15
  200 IF(IC.EQ.0)GO TO 230
      DO I=1,N
        DO J=I,N
          IF(EIG(I)-EIG(J))225,225,210
  210     X=EIG(I)
          EIG(I)=EIG(J)
          EIG(J)=X
          DO K=1,N
            Y=B(K,I)
            B(K,I)=B(K,J)
            B(K,J)=Y
          END DO
  225     CONTINUE
        END DO
      END DO
  230 CONTINUE
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE JACVAO(A,N,EIVU,EIVR,THRSH)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     V 2.0 - VECTORIZED VERSION - BERND HESS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(N,N),EIVU(N),EIVR(N,N)

      IF(N-1) 2,2,1
    2 EIVR(1,1)=1.0D0
      EIVU(1)=A(1,1)
      RETURN
    1 DO J=1,N
        DO I=1,N
          EIVR(I,J)=0.0D0
        END DO
        EIVR(J,J)=1.0D0
      END DO
C        FIND THE ABSOLUTELY LARGEST ELEMENT OF A
  102 ATOP=0.0D0
      DO I=1,N
C
C      THE FOLLOWING SEQUENCE SHOULD EVENTUALLY BE REPLACED BY
C      A BUILT IN ROUTINE ...
        DO J=I,N
CMI       IF ((ATOP-DABS(A(I,J))) < 0.0d0) THEN
          IF ((ATOP-DABS(A(I,J))) .LT.  0.0d0) THEN
            ATOP=DABS(A(I,J))
          END IF
        END DO
      END DO
      CONTINUE

      IF(ATOP)109,109,113
  109 RETURN
C         CALCULATE THE STOPPING CRITERION -- DSTOP
  113 AVGF= DFLOAT(N*(N-1))*.55D0
      D=0.0D0
      DO JJ=2,N
        DO II=2,JJ
          S=A(II-1,JJ)/ATOP
          D=S*S+D
        END DO
      END DO
      DSTOP=(1.D-09)*D
C        CALCULATE THE THRESHOLD, THRSH
C
      THRSH=DSQRT(D/AVGF)*ATOP
C
C         START A SWEEP
C
  115 IFLAG=0
      DO JCOL=2,N
         JCOL1=JCOL-1
         DO IROW=1,JCOL1
            AIJ=A(IROW,JCOL)
C
C     COMPARE THE OFF-DIAGONAL ELEMENT WITH THRSH
C
            IF (DABS(AIJ)-THRSH) 130,130,117
 117        AII=A(IROW,IROW)
            AJJ=A(JCOL,JCOL)
            S=AJJ-AII
C
C     CHECK TO SEE IF THE CHOSEN ROTATION IS LESS THAN THE ROUNDING E
C     IF SO , THEN DO NOT ROTATE.
C
            IF (DABS(AIJ)-1.D-09*DABS(S)) 130,130,118
 118        IFLAG=1
C
C     IF THE ROTATION IS VERY CLOSE TO 45 DEGREES, SET SIN AND COS
C     TO 1/(ROOT 2).
C
            IF (1.0D-10*DABS(AIJ)-DABS(S)) 116,119,119
 119        S=.707106781186548D0
            C=S
            GO TO 120
C
C     CALCULATION OF SIN AND COS FOR ROTATION THAT IS NOT VERY CLOSE
C     TO 45 DEGREES
C
 116        T=AIJ/S
            S=0.25D0/DSQRT(0.25D0+T*T)
C
C     COS=C,  SIN=S
C
            C=DSQRT(0.5D0+S)
            S=2.0D0*T*S/C
C
C     CALCULATION OF THE NEW ELEMENTS OF MATRIX A
C
 120        DO I=1,IROW
               T=A(I,IROW)
               U=A(I,JCOL)
               A(I,IROW)=C*T-S*U
               A(I,JCOL)=S*T+C*U
            END DO
            I2=IROW+2
            IF(I2.GT.JCOL) GO TO 123

            I2M1=I2-1
            JCM1=JCOL-1
C     ICI
C     ICI  ATTENTION
C     ICI
            DO I=I2M1,JCM1
               T=A(I,JCOL)
               U=A(IROW,I)
               A(I,JCOL)=S*U+C*T
               A(IROW,I)=C*U-S*T
            END DO
 123        A(JCOL,JCOL)=S*AIJ+C*AJJ
            A(IROW,IROW)=C*A(IROW,IROW)-S*(C*AIJ-S*AJJ)
            DO J=JCOL,N
               T=A(IROW,J)
               U=A(JCOL,J)
               A(IROW,J)=C*T-S*U
               A(JCOL,J)=S*T+C*U
            END DO
C
C     ROTATION COMPLETED.
C
 131        DO I=1,N
               T=EIVR(I,IROW)
               EIVR(I,IROW)=C*T-EIVR(I,JCOL)*S
               EIVR(I,JCOL)=S*T+EIVR(I,JCOL)*C
            END DO
C
C     CALCULATE THE NEW NORM D AND COMPARE WITH DSTOP
C
            S=AIJ/ATOP
            D=D-S*S
            IF(D.GE.DSTOP) GO TO 129
C
C     RECALCULATE DSTOP AND THRSH TO DISCARD ROUNDING ERRORS
C
            D=0.D0
            DO JJ=2,N
               DO II=2,JJ
                  S=A(II-1,JJ)/ATOP
                  D=S*S+D
               END DO
            END DO
            DSTOP=(1.0D-09)*D
 129        THRSH=DSQRT(D/AVGF)*ATOP
  130       CONTINUE
         END DO
      END DO
      IF(IFLAG.GT.0) GO TO 115
C
C      PLACE EIGENVALUES IN EIVU
      DO J=1,N
        EIVU(J)=A(J,J)
      END DO

      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      SUBROUTINE JDAAA
     *(P,Q,S,U,T,DT,V,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,W,THRE,EPS,EPSI,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cvec  VERSION INCLUDING IDK ARRAY
C
C     THIS SUBROUTINE  DIAGONALIZES THE HARTREE-FOCK MATRIX  FC AND FO
C         SYMMETRY BY SYMMETRY,NBAS(N) IS SMALLER OR EQUAL TO 15
C
C
C  C(MXBF,1)  ->    C(MXBF,*)
C   U(2)  -> U(*)
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "scfarr.h"
#include "relscf_ijpair.h"
#include "relscf_dim.h"
#include "relscf_v10.h"
CMI   DIMENSION P(2),Q(2),S(2),U(2),T(2),DT(2),V(2),PCAP(2),
      DIMENSION P(*),Q(*),S(*),U(*),T(*),DT(*),V(*),PCAP(*),
CMI  *QCAP(*),FC(*),FO(*),CM1SCF(MXBF,1),CM2SCF(MXBF,1),
     *QCAP(*),FC(*),FO(*),CM1SCF(MXBF,*),CM2SCF(MXBF,*),
     *CSAV1(MXBF,1),CSAV2(MXBF,1),GMIN(*),W(*),THRE(*),EPS(*),
     *EPSI(*),CM1(*),Z1(*),Z2(*),DELZ(*),
CMI  *C(MXBF,1),CTRAN(2),PQN(*),ZETA(*),ZETINC(*),
     *C(MXBF,*),CTRAN(2),PQN(*),ZETA(*),ZETINC(*),
     *ITRAN(2),NROW(*),NBVAR(*),IORCN(*)

      M1=0
      M=0
      NSTEP1=1
      DO N=1,NSYM
        NSTEP2=1
        NBAS1=NBAS(N)
        DO I=1,NBAS1
          NSTEP=NBAS1*(I-1)
          MM=M
          DO J=1,I
            IJ=NSTEP+J
            MM=MM+1

C           write(LUPRI,*) 'IJ=',IJ, 'LPQ=',LPQ
            if (MM.GT.LPQ) call quit('MM > LPQ !')
            DT(IJ)=S(MM)

          END DO
          JI=I
          DO J=1,I
            M=M+1

            if (M.GT.LPQ) call quit('M > LPQ !')
            DT(JI)=S(M)

            JI=JI+NBAS1
          END DO
        END DO
C
C     NEW (EIGRS) DIAGONALIZATION ROUTINE
C
        CALL JACVAO (DT,NBAS1,EPSI,V,THRSH)

Cvec  INTRODUCE NEW VECTOR LOOP (INTERMEDIATE STORAGE OF TERM)
        eigval_min = 1.0d0
        DO I=1,NBAS1
          eigval_min = min(eigval_min,EPSI(I))
        END DO
        if(eigval_min.le.0.0d0) call quit(' *** linear dependency 
     &     detected you may want check your basis set input!***')
        DO I=1,NBAS1
          DT(I)=1.D0/DSQRT(EPSI(I))
        END DO
        DO I=1,NBAS1
          NSTEP=NBAS1*(I-1)
          DO J=1,NBAS1
            JI=NSTEP+J
            V(JI)=V(JI)*DT(I)
          END DO
        END DO
C
C**** V CONTAINS  VT/SQRT(S0) (S0=EIGENVALUE MATRIX)
C
        NOSHIC=NOSH(N)
        NCSHIC=NCSH(N)
        IF (NOSHIC.EQ.0) GO TO 20
        IF (NCSHIC.EQ.0) GO TO 20
        IFIN=2
        GO TO 21
   20   IFIN=1
C
C**** FORM DT=V*F*VT   (F IS EQUAL TO FC OR FO)
C
   21   DO ISH=1,IFIN
           DO I=1,NBAS1
              NSTEP=NBAS1*(I-1)
              DO J=1,I
                 JI=NSTEP+J
                 IK=NBAS1*(J-1)
                 IJ=IK+I
                 TERM=0.D0
                 DO K=1,NBAS1
                    IF (IFIN.EQ.1) GOTO 924
                    IF (ISH.LT.IFIN) GOTO 928
 929                IHLP=K-1
CICI
CICI  ATTENTION
CICI
                    DO L=1,IHLP
                       KL=IDK(K)+L+M1
                       LJ=IK+L
                       KI=NSTEP+K
                       TERM=TERM+V(KI)*V(LJ)*FO(KL)
                    END DO
                    DO L=K,NBAS1
                       KL=IDK(L)+K+M1
                       LJ=IK+L
                       KI=NSTEP+K
                       TERM=TERM+V(KI)*V(LJ)*FO(KL)
                    END DO
                    ILOOP=1
                    ITOS=0
                    GOTO 9
 924                IF (NOSHIC.GT.0) GOTO 929
 928                IHLP=K-1
CICI
CICI  ATTENTION
CICI
                    DO L=1,IHLP
                       KL=IDK(K)+L+M1
                       LJ=IK+L
                       KI=NSTEP+K
                       TERM=TERM+V(KI)*V(LJ)*FC(KL)
                    END DO
                    DO L=K,NBAS1
                       KL=IDK(L)+K+M1
                       LJ=IK+L
                       KI=NSTEP+K
                       TERM=TERM+V(KI)*V(LJ)*FC(KL)
                    END DO
                    ILOOP=NCSHIC
                    ITOS=1
 9                  CONTINUE
                 END DO
                 DT(IJ)=TERM
                 DT(JI)=TERM
              END DO
           END DO
C     
C**** DIAGONALIZES DT AND TRANSFORMATION FOR VECTORS DT=VT*U
C     
           CALL JACVAO(DT,NBAS1,EPSI,W,THRSH)

           DO I=1,NBAS1
              DO J=1,NBAS1
                 TERM=0.D0
                 NSTEP=NBAS1*(J-1)
                 JI=NSTEP+I
                 DO K=1,NBAS1
                    IK=NBAS1*(K-1)+I
                    KJ=NSTEP+K
                    TERM=TERM+V(IK)*W(KJ)
                 END DO
                 DT(JI)=TERM
              END DO
C     
           END DO
C     
C     
C**** ORDER EIGENVECTORS AND EIGENVALUES WITHIN SYMMETRY
           DO I=1,ILOOP
              IND=1
 29           CONTINUE
              iMIN=1
              DO J=2,NBAS1
                 IF (EPSI(iMIN).GT.EPSI(J)) iMIN=J
              END DO
              IF (ITOS.EQ.1) GO TO 28
              IF (IND.EQ.NSTEP2) GO TO 28
              EPSI(iMIN)=0.1D+39
              IND=IND+1
              GO TO 29
 28           EPS(NSTEP1)=EPSI(iMIN)
              EPSI(iMIN)=0.1D+39
              DO J=1,NBAS1
                 JJ=NBAS1*(iMIN-1)+J
                 C(J,NSTEP1)=DT(JJ)
              END DO
              NSTEP1=NSTEP1+1
              NSTEP2=NSTEP2+1
           END DO
        END DO
        M1=N1(N)+M1
      END DO
      RETURN
      END

      SUBROUTINE KBR(BETA,X,RESULT)
#include "implicit.h"
#include "priunit.h"
C
C     PROGRAM TO CALCULATE CONTINUED FRACTION
C     1/1+ A1*X/1+ A2*X/1+ ...
C
C     V 1.0 - 12.5.86 - BERND HESS
C
C     A(2N+1) = N+BETA
C     A(2N)   = N-1/2
C
C#    WRITE (6,900) BETA,X
C#900 FORMAT(' KBR - BETA,X ',2D16.8)
      THR=1.D-12
      THR1=1.D-8
      SUM=1.D0
      I=100
1     OLD=SUM
      SUM=1.D0
      J=I
3     SUM=X*(DFLOAT(J)-0.5D0)/SUM+1.D0
      SUM=X*(DFLOAT(J)+BETA )/SUM+1.D0
      J=J-1
      IF (J.GT.0) GOTO 3
      SUM=1.D0/(1.D0-0.5D0*X/SUM)
      IF (I.GT.2000) WRITE (6,991) I,BETA,X,SUM
C#    WRITE (6,991) I,BETA,X,SUM
  991 FORMAT(' KBR - I,BETA,X,SUM ',I5,F10.3,D20.10,D30.20)
      I=I+100
      IF (I.LT.300 .OR. DABS(OLD-SUM).GT.THR.AND. I.LT.3000) GOTO 1
C#    WRITE (6,991) I,SUM
      RESULT=SUM
      IF (I.LT.3000) RETURN
      DEL=OLD-SUM
      WRITE (6,999) DEL,THR
999   FORMAT(' CONTINUED FRACTION DEL=',D20.10,' LARGER THAN THR=',
     *D20.10)
      IF (DEL.GT.THR1) STOP 999
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   
CMI   PROGRAM MAIN
      SUBROUTINE INT_AT_SCF(WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C  DELC@
C     CRAY CHANGES ARE MARKED BY c
C     RELATIVISTIC VERSION
C     VERSION 1.2 - 14.2.87 - EIGRS
C     IBM+ VERSION - 24.10.86 - BERND HESS
C     VECTORIZATION CHANGES MARKED BY Cvec
C     THIS VERSION WORKS FOR F FUNCTIONS
C>@V  VERSION
C     FOLLOWING SUBROUTINES HAVE BEEN OPTIMIZED:
C
C     HAMIL     VECTORIZED
C     GINY      EXPONENTIATIONS REPLACED BY MULTIPLICATIONS
C     JDAAA     VECTORIZED, ARRAY IDK INTRODUCED
C     TEIPSG    EXPONENTIATIONS REPLACED BY MULTIPLICATIONS
C     JACVAO    VECTORIZED
C     AT3DIM    ARRAY IDK INITIALIZED
C
C     V 1.0 - 14.8.86 - VARIABLE DIMENSIONS - BERND HESS
C     DYNAMIC ALLOCATION OF STORAGE ADAPTED FROM G.C.LIE'S VERSION
C     UP(99) DOWN(99) NOT DYNAMIC - SHOULD BE OF DIM. NBF
C     INPUT ARRAYS IN S.R. INPUT ALSO STATIC
C     THIS VERSION SUPPORTS 99 BASIS FUNCTIONS AND 20 OCCUPIED
C     ORBITALS (CLOSED PLUS OPEN)
C     TO INCREASE DIMENSIONS CHANGE STATEMENTS ENCLOSED BY
C     C##DIM ... C##
C
C------------------------------------------------------------------------------
C  Instead of beeing separated program called from AMFI routine getAOs2.
C
C  Saue, Ilias : connected to DIRAC array WORK (a lot of memory saving)
C
CMI   PARAMETER (MAINST=1500000)
C     PARAMETER (MAINST=2400000)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "amfi_if.h"
c     logical relscf_verbose

      REAL*8 GJ0,GK0
      REAL*8 WORK
      DIMENSION WORK(LWORK)

!     if (relscf_verbose) WRITE (LUPRI,100)
C     if (IPR_RELSCF.GE.2) WRITE (LUPRI,100)
C100   FORMAT(' *** PROGRAM AT34 - ALLIANT - @V ***'/)
      IF (IPR_RELSCF.GE.0) THEN
        CALL FLSHFO(LUPRI)
        CALL HEADER('*** PROGRAM AT34 - ALLIANT - @V ***',-1)
      ENDIF

C     CALL FLSHFO(LUPRI)

      CALL AT3(WORK,LWORK)

      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE OEISG(IPARAM,
     * P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,REL,Z100,
     *    S100,T100,U100,ZADD,C,ITRAN,CTRAN,NROW,PQN,
     &    ZETA,NBVAR,ZETINC,IORCN,NDIM,UN100,TN100,HCORR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C  REL(2) -> REL(*) because REL is from MEMGET
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "amfi_if.h"
#include "relscf_dim.h"
#include "relscf_v10.h"
C
      DIMENSION P(*),Q(*),S(*),U(*),T(*),
     &          DT(*),DOS(*),PCAP(*),QCAP(*),FC(*),FO(*),
     &          CM1SCF(MXBF,1),CM2SCF(MXBF,1),
     &          CSAV1(MXBF,1),CSAV2(MXBF,1),GMIN(*),SMIN(*),
     &          THRE(*),EPS(*),EC(*),CM1(*),Z1(*),Z2(*),DELZ(*),
     &          REL(*),Z100(*),S100(*),T100(*),U100(*),ZADD(MX100,4),
     &          C(MXBF,*),CTRAN(*),PQN(*),ZETA(*),ZETINC(*),
     &          ITRAN(*),NROW(*),NBVAR(*),IORCN(*)
C
C     ONE CONFIGURATION
C     ONE ELECTRON INTEGRAL PROGRAM
C     COMPUTES THE MATRICES S,U,T WITH SLATER (NFLAG1=0) OR GAUSSIAN
C      (NFLAG1=1) ORBITALS AS  BASIS
C
      dimension un100(NDIM), tn100(NDIM), hcorr(NDIM)
      NSTEP1=0
      NSTEP=0
      rewind 11

      DO L=1,NSYM
        NBAS1=NBAS(L)
        I1=NSTEP1
C
C     LARGER BASIS SET
C
        IF (TRHI(L).NE.0.D0) THEN
C     
C     DETERMINE LARGEST AND SMALLEST EXPONENT
C     
           I1=NSTEP1
           ZSM=ZETA(I1+1)
           ZLA=ZETA(I1+1)
           NP=PQN(I1+1)
           DO I=1,NBAS1
              I1=I1+1
              NQ=PQN(I1)
              IF (NQ.NE.NP) THEN
                 WRITE (LUPRI,*)
     &           '##### DIFFERENT PQN NUMBERS WITHIN ONE SHELL ####'
                 STOP 'PQN ERROR'
              ENDIF
              IF (ZETA(I1).GT.ZLA) ZLA=ZETA(I1)
              IF (ZETA(I1).LT.ZSM) ZSM=ZETA(I1)
           END DO
           ILA=MX100-NBAS1
           ISM=ILA/2
           ILA=MX100-NBAS1-ISM
           I1=0
           IF(NEXPND.EQ.1) THEN
              DO I=1,ISM
                 ZSM=ZSM/4.D0
                 IF (ZSM.LT.TRLO(L)) GOTO 22
                 I1=I1+1
                 ZADD(I1,L)=ZSM
 22              CONTINUE
              END DO
              DO I=1,ILA
                 ZLA=ZLA*4.D0
                 IF (ZLA.GT.TRHI(L)) GOTO 23
                 I1=I1+1
                 ZADD(I1,L)=ZLA
 23              CONTINUE
              END DO
              TRHI(L)=0.D0
              WRITE (6,250) I1,L,(ZADD(I,L),I=1,I1)
 250          FORMAT(I5,' ADDITIONAL FUNCTIONS GENERATED FOR SYMMETRY',
     *               I3/,(1X,8D16.4))
           ENDIF
           ISIZE(L)=I1
        ENDIF
C
C     BASIS SET IS NOW DETERMINED
C
        DO I=1,NBAS1
          Z100(I)=ZETA(I+NSTEP1)
        END DO
        ILA=NBAS1
        DO I=1,ISIZE(L)
          ILA=ILA+1
          Z100(ILA)=ZADD(I,L)
        END DO
        IJ=0
        NP=PQN(NSTEP1+1)
        NQ=NP
        DO I=1,ILA
          ZP=Z100(I)
          WP=2.D0*(NP-L)/ZP
          DO J=1,I
            IJ=IJ+1
            K=NSTEP+(I*(I-1)/2)+J
            ZQ=Z100(J)
            ZPQ=0.5D0*(ZP+ZQ)
            NPQ=NP+NQ+1
            WQ=2.D0*(NQ-L)/ZQ
            NPQ1=NPQ-1
            NPQ2=NPQ-2
C************************************************************
C*                                                          *
C*     CALCULATE INTEGRALS FOR GAUSSIAN TYPE FUNCTIONS      *
C*                                                          *
C************************************************************
    1       VPQ=FACTO(NPQ1)/ZPQ**(0.5D0*NPQ)
            VP=FACTO(2*NP)/ZP**(NP+0.5D0)
            VQ=FACTO(2*NQ)/ZQ**(NQ+0.5D0)
            IF (NPQ1.GT.2) GO TO 3
            VPQM2=1.D0
            GO TO 4
    3       VPQM2=FACTO(NPQ2-1)/ZPQ**(0.5D0*NPQ2)
    4       VPQ1=FACTO(NPQ2)*1.595769121605731D0/ZPQ**(0.5D0*NPQ1)
            VPQP2=FACTO(NPQ+1)/ZPQ**(0.5D0*(NPQ+2))
            TERM2=VPQP2-VPQ*(WP+WQ)+WP*WQ*VPQM2
            TERM1=1.D0/DSQRT(VP*VQ)

            S100(IJ)  = TERM1*VPQ
            U100(IJ)  = TERM1*VPQ1
            UN100(IJ) = TERM1*VPQ1
            TN100(IJ) = 0.5D0*ZP*ZQ*TERM1*TERM2
            IF (IPARAM.GT.0) THEN
               I1=np-1
               I2=0
               I3=0
               J1=nq-1
               J2=0
               J3=0
C*
C*     PVP FOR XY-FUNCTIONS
C*
               if (np.eq.3) then
                  REL(IJ)=EXTC(L,ZP,ZQ,1,1,0,1,1,0)
                  GOTO 900
               endif
C*
C*     PVP FOR XYZ-FUNCTIONS
C*
               if (np.eq.4) then
                  REL(IJ)=EXTC(L,ZP,ZQ,1,1,1,1,1,1)
                  GOTO 900
               endif
C*
C*     PVP FOR X**(NP-1)-FUNCTIONS
C*
               REL(IJ)  = EXTC(L,ZP,ZQ,I1,I2,I3,J1,J2,J3) ! TROND: ok
 900           CONTINUE
               T100(IJ) = 0.5D0*ZP*ZQ*TERM1*TERM2 ! TROND: ok
            ELSE
               T100(IJ)=0.5D0*ZP*ZQ*TERM1*TERM2
            ENDIF
          END DO
        END DO

        IF (IPARAM.EQ.1) THEN
C
C     SQR OPERATOR - SAVE NONREL POTENTIAL RIGHTAWAY
C
           K=NBAS1*(NBAS1+1)/2
           DO I=1,K
              U(NSTEP+I)=U100(I)
              write(LUPRI,*) 'U index=',NSTEP+I
           END DO
        ENDIF

        IF (IPARAM.GT.0) THEN
C
C     CONVERT TO RELATIVISTIC INTEGRALS
C
          NBSZ=ILA*(ILA+1)/2
          NBSQ=ILA*ILA
          KEV2  = 1
          KBU   = KEV2  + NBSZ
          KP    = KBU   + NBSZ
          KG    = KP    + NBSZ
          KEIG  = KG    + NBSZ
          KSINV = KEIG  + NBSQ
          KREVT = KSINV + NBSQ
          KAUX  = KREVT + NBSQ
          KOVE  = KAUX  + NBSQ
          KEW   = KOVE  + NBSQ
          KE    = KEW   + ILA
          KAA   = KE    + ILA
          KRR   = KAA   + ILA
          KTT   = KRR   + ILA
          KLAST = KTT   + ILA
          CALL AT34R(ILA,NBSZ,ZN,S100,U100,T100,REL(KEV2),
     *               REL(KBU),REL(KP),REL(KG),REL(KEIG),REL(KSINV),
     *               REL(KREVT),REL(KAUX),REL(KOVE),REL(KEW),
     *               REL(KE),REL(KAA),REL(KRR),REL(KTT))
        ENDIF

        K=NBAS1*(NBAS1+1)/2
C
C     TRANSFER RELATIVISTIC KINETIC ENERGY AND POTENTIAL INTEGRALS
C
        IF (IPARAM.EQ.1) THEN
          DO I=1,K
            T(NSTEP+I)=T100(I)
            S(NSTEP+I)=S100(I)
          END DO
        ELSE
          DO I=1,K
            T(NSTEP+I)=T100(I)
            S(NSTEP+I)=S100(I)
            U(NSTEP+I)=U100(I)
          END DO
        ENDIF
        NSTEP1=NSTEP1+NBAS1
        ij=0
        DO i=1,nbas1
          DO j=1,i
            ij=ij+1
            hcorr(ij)=(t100(ij)-zn*u100(ij)) - 
     &                (tn100(ij) - zn*un100(ij))
          END DO
        END DO
        dum=0.0D0
        NSTEP=NSTEP+K
      END DO ! DO L=1,NSYM
       kbeg=1
 27    format(d24.18)
 31    format(5i5)

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE OUTPUT_RELSCF
     *(P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C MI: array dim. fixes:
C----------------------------
C  C(MXBF,1) -> C(MXBF,*)
C  Q(2)      -> Q(*)
C  P(2)      -> P(*)
C  U(2)       --> U(*)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "scfarr.h"
#include "relscf_dim.h"
#include "relscf_v10.h"
#include "dummy.h"

CMI   DIMENSION P(2),Q(2),S(2),U(2),T(2),DT(2),DOS(2),PCAP(2),
      DIMENSION P(*),Q(*),S(*),U(*),T(*),DT(*),DOS(*),PCAP(*),
     *QCAP(*),FC(*),FO(*),CM1SCF(MXBF,1),CM2SCF(MXBF,1),
     *CSAV1(MXBF,1),CSAV2(MXBF,1),GMIN(*),SMIN(*),THRE(*),EPS(*),
     *EC(*),CM1(*),Z1(*),Z2(*),DELZ(*),
CMI  *C(MXBF,1),CTRAN(2),PQN(*),ZETA(*),ZETINC(*),
     *C(MXBF,*),CTRAN(*),PQN(*),ZETA(*),ZETINC(*),
     *ITRAN(*),NROW(*),NBVAR(*),IORCN(*)
      DIMENSION WORK(LWORK)

      COMMON/WIGXXX/ A1,A2,EXP1,EXP2,UP(MXB),DOWN(MXB),IUP(4),NGAUSS
C
      COMMON /CHAR/ NAME(12)
      CHARACTER*4 NAME
C
      DIMENSION NOOC(5)
      DIMENSION NTYSM(4),ATCO(300)
C     PUNCH BASIS SET AND VECTORS FOR INPUT TO IBMOL-H
      LOGICAL LLIE(4), EX
#include "amfi_if.h"
c     logical relscf_verbose
CMI   DATA LLIE/1HS,1HP,1HD,1HF/, NTYSM/1,3,6,10/
      DATA LLIE/.false.,.false.,.false.,.false./, NTYSM/1,3,6,10/
#include "memint.h"
      IDCTR=MXCC*NEWBAS
C
C     write out vectors to file RELSCF_COEF
CSK   NSYM is equal to highest angular momentum orbital occupied, i.e.
C     if f^2 --> NSYM = 4
      CALL MEMGET('REAL',KRELCOEFFX,LDIMC+2,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KRELCOEFFX),LDIMC+2)
      CALL IZERO(IRELSCF_INFO,MXRELSCF_INFO)

      K          = 0
      LINF_X     = 1
      LRELSCFVEC = 1
      IRELSCF_INFO(LINF_X) = NSYM
      LINF_X = LINF_X + 1
      DO L=1,NSYM
        LNBF=NBAS(L)
        NTSHL=NCSH(L)+NOSH(L)
C        
        IRELSCF_INFO(LINF_X) = LNBF
        LINF_X = LINF_X + 1
        IRELSCF_INFO(LINF_X) = NTSHL
        LINF_X = LINF_X + 1
C
        DO J=1,LNBF
          IF (IPR_RELSCF.GE.2)WRITE (LUPRI,5037)(C(J,K1),K1=K+1,K+NTSHL)
          DO K1=K+1,K+NTSHL
            WORK(KRELCOEFFX+LRELSCFVEC-1) = C(J,K1)
            LRELSCFVEC = LRELSCFVEC + 1
          END DO
        END DO
        K=K+NTSHL
      END DO
 5037 FORMAT(4F18.14)
C
C     check if old 'RELSCF_COEF' exists.
      LURELSCFCOEF = 7
      INQUIRE ( FILE = 'RELSCF_COEF', EXIST = EX )
      IF( EX ) THEN
        CALL OPNFIL(LURELSCFCOEF,'RELSCF_COEF','OLD','OUTREL')
      ELSE
        CALL OPNFIL(LURELSCFCOEF,'RELSCF_COEF','UNKNOWN','OUTREL')
        CALL NEWLAB('SODLABEL',LURELSCFCOEF,LUPRI)
      END IF
      REWIND LURELSCFCOEF
      
      ISIZEV = LRELSCFVEC - 1
      IRELSCF_INFO(LINF_X) = ISIZEV
C     WRITE(LUPRI,*) ' array IRELSCF_INFO'
C     CALL IWRTMAMN(IRELSCF_INFO,1,MXRELSCF_INFO,1,MXRELSCF_INFO,LUPRI)
C     WRITE(LUPRI,*) ' ISIZEV is',ISIZEV
      CALL IWRTAMFI(LURELSCFCOEF,RELSCFLABI,IRELSCF_INFO,MXRELSCF_INFO)
C
C     write coefficients to file RELSCF_COEF
C     WRITE(LUPRI,*) ' coeff matrix to file'
C     CALL WRTMATMN(WORK(KRELCOEFFX),1,ISIZEV,1,ISIZEV,LUPRI)
      CALL WRTAMFI(LURELSCFCOEF,RELSCFLABX,WORK(KRELCOEFFX),ISIZEV)
      CLOSE(UNIT=LURELSCFCOEF,STATUS='KEEP')
C     ... release memmory
      CALL MEMREL('OUTPUT_RELSCF',WORK,KWORK,KWORK,KFREE,LFREE)

      IF(NPOEFA.EQ.0) GO TO 101
C     PUNCH ZETA OR CONTRACTION COEFFICIENTS
  502 FORMAT(2X,'ZETA=',7X,4(E13.6,','))
  503 FORMAT(2X,'C=',2X,4(F17.6,','))
  512 FORMAT(2X,'CTRAN=',2X,6(F10.6,','))


  101 CONTINUE

CMI    return ! ...good return !!!

      DO I=1,NSYM
        NBCON(I)=NSAVF(I)
      END DO
      IF(NPBAFI.EQ.0) GO TO 166
      JJ=1
      K=0
      L=0
      DO I=1,NSYM
        K=K+NBCON(I)*NTYSM(I)
        JMAX=NBCON(I)
        DO J=1,JMAX
          L=L+NROW(JJ)*NTYSM(I)
          JJ=JJ+1
        END DO
      END DO
      I=ZN+0.01D0

  165 CONTINUE

CMI    return ! ...bad return

      JJ=1
      KMIN=0
      DO I=1,NSYM
        JMAX=NBCON(I)
        DO J=1,JMAX
          KMAX=KMIN+NROW(JJ)
          KMIN=KMIN+1
          DO K=KMIN,KMAX
            L=ITRAN(K)
            L=K
CMI   WRITE (1,112) LLIE(I),ZETA(L),CTRAN(K)
  112       FORMAT(A1,F17.6,4X,F9.6)
          END DO
          KMIN=KMIN+MXCC-1
          JJ=JJ+1
        END DO
      END DO
C
  166 CONTINUE

CMI     return ! good return !!!

      IF(NPCOFI.EQ.0) GO TO 121
      NGF=0
      DO I=1,NSYM
        NGF=NGF+NBCON(I)*NTYSM(I)
      END DO

CMI    return ! bad return
C
      JM=0
      JMAX=0
      KK=0
      DO I=1,NSYM
        KMAX=NCSH(I)+NOSH(I)
        JI=JMAX+1
        JM=JMAX+NBCON(I)
        NTYI=NTYSM(I)
        DO K=1,KMAX
          KK=KK+1
          DO L=1,NTYI
            DO M=1,NGF
              ATCO(M)=0.0D0
            END DO
            J1=JI+(L-1)*NBCON(I)
            JMAX=JM+(L-1)*NBCON(I)
            JJ=0
            DO J=J1,JMAX
              JJ=JJ+1
              ATCO(J)=C(JJ,KK)
            END DO
          END DO
        END DO
      END DO
   22 FORMAT(11X,6(F10.6,','))
  122 FORMAT(2X,'ATCOEF=',2X,6(F10.6,','))


CMI    return ! bad return !
C
  121 IF (NALARM.EQ.0) GO TO 4
      IF (NPRINT.NE.0) GO TO 5

      IF (IPR_RELSCF.GE.2) WRITE (LUPRI,15) (S(I),I=1,N1T)
   15 FORMAT('1'//,5X,'S MATRIX'//(6X,10F12.6))

      GO TO 106

    4 IF (NPRINT.EQ.0) GO TO 106

    5 if (IPR_RELSCF.GE.2) WRITE (LUPRI,10) (S(I),I=1,N1T)
   10 FORMAT ('1'//,5X,'MATRICES AND SUPERMATRICES OF INTEGRALS BETWEEN
     1BASIS FUNCTIONS'//,6X,'S MATRIX'//(6X,10F12.6))
   25 FORMAT(6X,10F12.6)

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,30) (U(I),I=1,N1T)
   30 FORMAT ('0',5X,'U MATRIX'/(6X,10F12.6))

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,35)    (T(I),I=1,N1T)
   35 FORMAT ('0',5X,'T MATRIX'/(6X,10F12.6))
C
      IBLK=1
      IPQ=N1T*(N1T+1)/2
      IBM=IPQ/NDIMPQ+1
      IF(IBM.EQ.1) GO TO 33
      REWIND IDATA
      IPQ=NDIMPQ
      IBLK=0
      DO I=1,N1T
        DO K=1,I
          IF(IPQ-NDIMPQ) 31,32,32
   32     READ (IDATA) (P(LIE),LIE=1,NDIMPQ),(Q(LIE),LIE=1,NDIMPQ)
          IBLK=IBLK+1
          IPQ=0
          IF(IBLK.EQ.IBM)GOTO 31

          if (IPR_RELSCF.GE.2) THEN
            WRITE (LUPRI,40) IBLK,(P(LIE),LIE=1,NDIMPQ)
   40     FORMAT('1',5X,'P SUPERMATRIX BLOCK ',I4//(6X,10F12.6))
          endif
          if (IPR_RELSCF.GE.2) then
            WRITE (LUPRI,45) IBLK,(Q(LIE),LIE=1,NDIMPQ)
   45     FORMAT('1',5X,'Q SUPERMATRIX BLOCK ',I4//(6X,10F12.6))
          endif
 31       IPQ=IPQ+1
        END DO
      END DO
   33 if (IPR_RELSCF.GE.2) WRITE(LUPRI,40)IBLK,(P(I),I=1,IPQ)
      if (IPR_RELSCF.GE.2) WRITE(LUPRI,45)IBLK,(Q(I),I=1,IPQ)

      IF(NDIMPQ.NE.NDIMPQPAR) REWIND IDATA
C
      if (IPR_RELSCF.GE.2) WRITE(LUPRI,60) (S(I),I=1,N1T)
   60 FORMAT ('1'//,5X,'FINAL MATRICES'//5X,'S MATRIX'/(6X,10F12.6))
      DO I=1,N1T
        U(I)=-ZN*U(I)+T(I)
      END DO

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,65) (U(I),I=1,N1T)
   65 FORMAT ('0',5X,'H MATRIX'/(6X,10F12.6))
      if (IPR_RELSCF.GE.2) WRITE(LUPRI,70) (PCAP(I),I=1,N1T)
   70 FORMAT ('0',5X,'P MATRIX'/(6X,10F12.6))
      if (IPR_RELSCF.GE.2) WRITE(LUPRI,75) (QCAP(I),I=1,N1T)
   75 FORMAT ('0',5X,'Q MATRIX'/(6X,10F12.6))
      if (IPR_RELSCF.GE.2) WRITE(LUPRI,80)    (FC(I),I=1,N1T)
   80 FORMAT ('0',5X,'CLOSED SHELL F MATRIX'/(6X,10F12.6))
      if (IPR_RELSCF.GE.2) WRITE(LUPRI,90)  (FO(I),I=1,N1T)
   90 FORMAT ('0',5X,'OPEN SHELL F MATRIX'/(6X,10F12.6))
      if (IPR_RELSCF.GE.2) WRITE(LUPRI,95)    (DT(I),I=1,N1T)
   95 FORMAT ('0',5X,'TOTAL DENSITY MATRIX'/(6X,10F12.6))
      if (IPR_RELSCF.GE.2) WRITE(LUPRI,105)   (DOS(I),I=1,N1T)
  105 FORMAT ('0',5X,'OPEN SHELL DENSITY MATRIX'/(6X,10F12.6))
C
  106 if (IPR_RELSCF.GE.2) WRITE(LUPRI,110) (NAME(I),I=1,12)
  110 FORMAT ('1',//,5X,'FINAL RESULTS'//5X,12A4)


      !  return ! good return

      if (IPR_RELSCF.GE.2) WRITE (LUPRI,13)
   13 FORMAT ('0',5X,19HTYPE GAUSS ORBITALS)
   14 IF (NSVD.EQ.1) GO TO 500

      if (IPR_RELSCF.GE.2) WRITE (6,504)
  504 FORMAT('0',5X,'SINGLE VECTOR DIAGONALIZATION')
      GO TO 506

  500 if (IPR_RELSCF.GE.2) WRITE (LUPRI,505)
  505 FORMAT ('0',5X,'JACOBI DIAGONALIZATION')


  123 FORMAT('0',5X,'CUT OFF PARAMETERS =',F10.6,3X,F10.6)
  425 FORMAT('0',5X,'WIGNER PARAMETERS   A1 =',F10.6,3X,'EXP1 =',F10.6,
     *       /,26X,'A2 =',F10.6,3X,'EXP2 =',F10.6)

  506 if (IPR_RELSCF.GE.2) WRITE (LUPRI,120) ZN
  120 FORMAT ('0',5X,'CHARGE =',F10.6,//,6X,
     1'SYMMETRY SPECIES',11X,'S',5X,'P',5X,'D',5X,'F')

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,130) (NBAS(I),I=1,NSYM)
  130 FORMAT(6X,'NUMBER OF BASIS FUNCTIONS=',I3,3(4X,I3))

!     write(lupri,*) 'NSYM=',NSYM
!     write(lupri,*) 'NBAS(I):',(NBAS(I),I=1,NSYM)

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,135) (NCSH(I),I=1,NSYM)
  135 FORMAT(6X,'NUMBER OF CLOSED SHELLS  =',I3,3(4X,I3))

!      return ! ...good return !

      DO I=1,NSYM
        NOOC(I)=OCCUP(I)+0.01D0
      END DO

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,140) (NOOC (I),I=1,NSYM)
  140 FORMAT(6X,'OPEN SHELL OCCUPATION    =',10(1X,I1,4X),/)

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,150)
  150 FORMAT ('0',5X,'VECTOR COUPLING COEFFICIENTS K'/)
      fact=2.

!      return ! ...good return !

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,151) (fact*AJMN(I),I=1,24)
  151 FORMAT(5D21.8)

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,152)
  152 FORMAT ('0',5X,'BASIS FUNCTIONS (PRINCIPAL QUANTUM NUMBER, ORBITAL
     1 EXPONENTS)'//,15X,'S',28X,'P',28X,'D',28X,'F'/)
  160 FORMAT(4(5X,F4.1,5X,F15.6))
C
      NBAS1=NSAVF(1)
      IF (NSYM.GT.1) GO TO 155
       DO   I=1,NBAS1
        if (IPR_RELSCF.GE.2)  WRITE(LUPRI,160)  PQN(I),ZETA(I)
       END DO
       GO TO 215
C
  155 NBAS2=NSAVF(2)
      IF (NSYM.GT.2) GO TO 156
      NBA3=1
  197 DO 185 I=NBA3,NBAS2
       J=NBAS1+I

  185 continue 
      if (IPR_RELSCF.GE.2) THEN 
        WRITE(LUPRI,160) PQN(I),ZETA(I),PQN(J),ZETA(J)
      endif

  196 IF (NBAS1.EQ.NBAS2) GO TO 215
       NBA2=NBAS2+1
       DO   I=NBA2,NBAS1
        if (IPR_RELSCF.GE.2) WRITE(LUPRI,160) PQN(I),ZETA(I)
       END DO
       GO TO 215
C
  156 NBAS3=NSAVF(3)
      IF(NSYM.GT.3) GO TO 157
      NBA4=1
  207 DO 195  I=NBA4,NBAS3
       J=NBAS1+I
       K=NBAS1+NBAS2+I

  195 if (IPR_RELSCF.GE.2) 
     &    WRITE(LUPRI,160) PQN(I),ZETA(I),PQN(J),ZETA(J),PQN(K),ZETA(K)

  206 IF (NBAS3.EQ.NBAS2) GO TO 196
      NBA3=NBAS3+1
      GO TO 197
C
  157 NBAS4=NSAVF(4)
      DO I=1,NBAS4
        J=NBAS1+I
        K=NBAS1+NBAS2+I
        L=NBAS1+NBAS2+NBAS3+I

       if (IPR_RELSCF.GE.2) 
     &  WRITE(LUPRI,160) PQN(I),ZETA(I),PQN(J),ZETA(J),PQN(K),ZETA(K),
     &             PQN(L),ZETA(L)

      END DO

      IF(NBAS4.EQ.NBAS3) GO TO 206
      NBA4=NBAS4+1
      GO TO 207

  215 if (IPR_RELSCF.GE.2) WRITE(LUPRI,220) NITSCF,THRSCF

  220 FORMAT('0',5X,'FINAL SCF RESULTS OBTAINED AT COMPUTATION NO',I2,
     *       5X,'SCF THRESHOLD =',E11.4)

      if (IPR_RELSCF.GE.2) WRITE (LUPRI,1225)

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,230) ENERG,POT,CIN,VIR

!      return ! ...good return !

CMI   IF(NPENFA.EQ.1) WRITE(1,1230) ENERG,VIR
 1225 FORMAT('0',5X,'TOTAL HF ENERGY',4X,'POTENTIAL ENERGY',4X,
     *'KINETIC ENERGY',4X,'VIRIAL THEOREM')
 1226 FORMAT('0',5X,'TOTAL HF ENERGY',4X,'POTENTIAL ENERGY',4X,
     *'KINETIC ENERGY',4X,'VIRIAL THEOREM',4X,'COR. ENERGY')
 1227 FORMAT('0',5X,'TOTAL CHF ENERGY',3X,'POTENTIAL ENERGY',4X,
     *'KINETIC ENERGY',4X,'VIRIAL THEOREM')
  230 FORMAT(4X,5(D19.10))
 1229 FORMAT(2X,'P.E.=',D16.10,',',2X,'K.E.=',D16.10,',')
 1230 FORMAT(2X,'ENERG=',D16.10,',',2X,'VIR=',D16.10,',')
 1231 FORMAT(2X,'ENERG=',D16.10,',',2X,'VIR=',D16.10,',',2X,'ENCOR=',
     *D16.10,',')
 1232 FORMAT(2X,'EPS=',5X,5(F12.5,','))
 1233 FORMAT(2X,'WIG=',F12.6,',',2X,'WIGNER=',F12.6,',')

!      return ! ...good return !

      GO TO 231

 1002 if (IPR_RELSCF.GE.2) WRITE(LUPRI,1226)

      if (IPR_RELSCF.GE.2) WRITE(LUPRI,230) ENERG,POT,CIN,VIR,ENCOR

CMI   IF(NPENFA.EQ.1) WRITE(1,1231) ENERG,VIR,ENCOR
  231 if (IPR_RELSCF.GE.2) WRITE(LUPRI,255)

!      return ! ...good return !

  255 FORMAT('0',5X,'ORBITAL ENERGIES, EIGENVECTORS, AND DIAGONALIZATION
     * THRESHOLD')

!      return ! ...good return !

CMI   IF(NPENFA.EQ.1) WRITE(1,1229) POT,CIN
      IMAX=0
      DO LIE=1,NSHT,8
        IMIN=IMAX+1
        IMAX=IMAX+8
        IF(IMAX.GT.NSHT) IMAX=NSHT

        IF (IPR_RELSCF.GE.2) THEN
          WRITE(LUPRI,265) (EPS(I),I=IMIN,IMAX)
        ENDIF
  265   FORMAT('0',5X,F11.5,7(4X,F10.5))

        if (IPR_RELSCF.GE.2) THEN
           WRITE (LUPRI,20)
        ENDIF
   20   FORMAT (' ')

        DO I=1,NBASM
          IF (IPR_RELSCF.GE.2) THEN
            WRITE(LUPRI,290) (C(I,J),J=IMIN,IMAX)
          ENDIF
        END DO
  290   FORMAT(3X,8F14.6)


        IF(NSVD.EQ.0) then 
           IF (IPR_RELSCF.GE.2) THEN
              WRITE(LUPRI,295) (THRE(J),J=IMIN,IMAX)
           ENDIF
        endif
  295   FORMAT('0',5X,E11.4,7(4X,E10.4))
      END DO

!      return ! ...good return !

      IF(NSVD.EQ.1) then 
         IF (IPR_RELSCF.GE.2) THEN
            WRITE(LUPRI,295) THRSH
         ENDIF
      endif

CMI   IF(NPENFA.EQ.1) WRITE(1,1232) (EPS(I),I=1,NSHT)
C
C     SET UP R FOR ORBITAL OR WIGNER
C
      IF(NORBIT.EQ.1) then 
         IF (IPR_RELSCF.GE.2) THEN
           WRITE (LUPRI,310)
         ENDIF
      endif
  310 FORMAT ('1'//,5X,'TABULATION OF ORBITAL WAVE FUNCTIONS'//9X,'R'/)

!      return ! ...GOOD return !

      II=1
      JJ=NBLOCK*NGAUSS+1
      RF=0.0D0
      DO I=1,NBLOCK
        RI=RF
        RF=RF+NSIZE(I)*RADINC(I)
        RM=(RI+RF)*0.5D0-RADINC(I)*5.0D0

CMI      ...   for array 'q', upper bound of dimension 1 exceeded ....
        LQM1=LQM*2
        if (II.GT.LQM1.OR.JJ.GT.LQM1) then 
         write(LUPRI,*) 'LQM=',LQM
         write(LUPRI,*) 'LQM1=',LQM1
         write(LUPRI,*) 'II=',II
         write(LUPRI,*) 'JJ=',JJ
         call quit('before INPAR due to Q! increase LQM1 for MEMGET!')
        endif

        CALL INPAR(1,NGAUSS,RI,RF,RM,Q(II),Q(JJ),0,DUMMY,DUMMY)
        II=II+NGAUSS
        JJ=JJ+NGAUSS
      END DO

!      return ! ...BAD return !

      NUMBER=II-1
      II=1
      E=0.0D0
      DO I=1,NSYM
        IF(NCSH(I).EQ.0) GO TO 405
        NSH=NCSH(I)
        DO J=1,NSH
          UP(II)=2*(I-1)+1
          DOWN(II)=UP(II)
          E=E+UP(II)+DOWN(II)
          II=II+1
        END DO
  405   IF(NOSH(I).EQ.0) GO TO 406
        DOWN(II)=(NOOC(I) -IUP(I))/2.
        UP(II)=NOOC(I) -DOWN(II)
        E=E+UP(II)+DOWN(II)
        II=II+1
 406    CONTINUE
      END DO

!      return ! ...BAD return !

      GO TO 420


  410 IF(NORBIT.NE.1) RETURN
      II=2
      Q(1)=0.0D0
      NUMBER=0
      DO I=1,NBLOCK
        NUMBER=NUMBER+NSIZE(I)
        ITER=NSIZE(I)
        IF(I.EQ.1) ITER=NSIZE(I)-1
        DO J=1,ITER
          Q(II)=Q(II-1)+RADINC(I)
          II=II+1
        END DO
      END DO
C
C     COMPUTE THE ORBITAL WAVE FUNCTIONS AND PRINT
C
  420 CONTINUE

!      return ! ...bad return !

      K=1
      DO ITER=1,NUMBER
        R=Q(ITER)
        P(K)=R
        Q(ITER)=Q(ITER)**2
        K=K+1
        NSTEP1=0
        NSTEP2=0
        DO I=1,NSYM
          NBAS1=NBAS(I)
          NSH=NCSH(I)+NOSH(I)
          DO J=1,NSH
            J1=J+NSTEP1
 
          if (J1.GT.LPQ) call quit('for S exceeded range !')

            S(J1)=0.D0
            DO L=1,NBAS1
              L1=L+NSTEP2
              NQP=PQN(L1)
              NQP1=2*NQP+1
  322         ENORM=0.797884560802865D0*(2.D0**NQP1)*ZETA(L1)**
     &              (0.5D0*NQP1)/FACTO(NQP1-1)
              ENORM=DSQRT(ENORM)
              TERM=ENORM*R**NQP* EXP(-ZETA(L1)*R*R)
              S(J1)=S(J1)+C(L,J1)*TERM
            END DO
          END DO
          NSTEP1=NSTEP1+NSH
          NSTEP2=NSTEP2+NBAS1
        END DO
        DO I=1,NSHT
CMI  ....  for array 'p', upper bound of dimension 1 exceeded
          P(K)=S(I)**2
          K=K+1
        END DO
        IF(NORBIT.EQ.1) then 
           IF (IPR_RELSCF.GE.2) THEN
              WRITE (LUPRI,330) R,(S(I),I=1,NSHT)
           ENDIF
        endif
  330   FORMAT(2X,9(4X,F9.6))
      END DO
C
C     COMPUTE THE ORBITAL DENSITIES
C
      IF(NORBIT.EQ.0) GO TO 421

      if (IPR_RELSCF.GE.2) THEN
         WRITE(LUPRI,360)
  360    FORMAT ('1'//,5X,'TABULATION OF ORBITAL DENSITIES'//9X,'R'/)
      ENDIF

      K=1
      DO ITER=1,NUMBER
        TOT=0.D0
        K1=K+1
        DO I=1,NSHT
          TOT=TOT+P(K1)
          K1=K1+1
        END DO
        K1=K1-1
        IF (IPR_RELSCF.GE.2) THEN
           WRITE(LUPRI,375) (P(I),I=K,K1),TOT
  375      FORMAT(3X,10(3X,F9.6))
        ENDIF

        K=K+NSHT+1
      END DO
C
C     WIGNER-TYPE ESTIMATE OF CORRELATION
C
  421 UPR=0.0D0
      DOWNR=0.0D0
      WIGNER=0.0D0
      II=2
      JJ=NGAUSS*NBLOCK+1
      DO I=1,NUMBER
        UPL=0.0D0
        DWNL=0.0D0
        DO J=1,NSHT
          UPL=UPL+P(II)*UP(J)
          DWNL=DWNL+P(II)*DOWN(J)
          II=II+1
        END DO
        UPR=UPR+UPL*Q(JJ)
        DOWNR=DOWNR+DWNL*Q(JJ)
        TEMP=6.332573978D-3*UPL*DWNL/Q(I)**2
        WIGNER=WIGNER+TEMP**EXP1/(A2+TEMP**EXP2)*Q(I)*Q(JJ)
        II=II+1
        JJ=JJ+1
      END DO
      SUM=UPR+DOWNR
      CNST=A1*(2.0D0-1.0D0/( EXP((E -7.0D0)*0.5D0)+1.0D0))
      WIG   =WIGNER   *12.56637061D0
      WIGNER=WIG*CNST
      if (IPR_RELSCF.GE.2) THEN
         WRITE(LUPRI,424) UPR,DOWNR,SUM,WIG,WIGNER
         WRITE(LUPRI,*) 
      ENDIF
  424 FORMAT(' SPIN UP   DENSITY =',F12.6,
     */,' SPIN DOWN DENSITY =',F12.6,
     */,' TOTAL     DENSITY =',F12.6,
     */,' WIGNER CORRECTION =',F12.6,5X,F12.6)
C
CMI   IF(NPENFA.EQ.1) WRITE(1,1233) WIG,WIGNER

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      DOUBLE PRECISION FUNCTION PHI(M,N)
#include "implicit.h"
#include "relscf_crelop.h"
      IF (MOD(N,2).EQ.1.OR.MOD(M,2).EQ.1) GOTO 10
      PHI=2.D0*GA(M+1)*GA(N+1)/GA(M+N+2)
      RETURN
10    PHI=0.D0
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PREN(IPARAM,
     * P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PRINT INTERMEDIATE ITERATION RESULTS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "relscf_dim.h"
CMI   DIMENSION P(2),Q(2),S(2),U(2),T(2),DT(2),DOS(2),PCAP(2),
      DIMENSION P(*),Q(*),S(*),U(*),T(*),DT(*),DOS(*),PCAP(*),
     *QCAP(*),FC(*),FO(*),CM1SCF(MXBF,*),CM2SCF(MXBF,*),
CMI  *CSAV1(MXBF,1),CSAV2(MXBF,1),GMIN(*),SMIN(*),THRE(*),EPS(*),
     *CSAV1(MXBF,*),CSAV2(MXBF,*),GMIN(*),SMIN(*),THRE(*),EPS(*),
     *EC(*),CM1(*),Z1(*),Z2(*),DELZ(*),
CMI  *C(MXBF,1),CTRAN(*),PQN(*),ZETA(*),ZETINC(*),
     *C(MXBF,*),CTRAN(*),PQN(*),ZETA(*),ZETINC(*),
     *ITRAN(*),NROW(*),NBVAR(*),IORCN(*)
      COMMON /CLOCKX/ CPU,CLK
#include "amfi_if.h"
c     logical relscf_verbose
c     relscf_verbose = .false.
C
CMI   IF (NITSCF.EQ.1.and.relscf_verbose) THEN
CMI   IF (NITSCF.EQ.1.and.IPR_RELSCF.GE.2) THEN
      IF (NITSCF.EQ.1.and.IPR_RELSCF.GE.0) THEN
         IF (IPARAM.EQ.0) THEN
            WRITE (LUPRI,200)
200         FORMAT(' --- NON-RELATIVISTIC RESULT ---')
         ENDIF
         IF (IPARAM.EQ.1) THEN
            WRITE (LUPRI,201)
201         FORMAT(' --- SQR MATRIX RESULT ---')
         ENDIF
         IF (IPARAM.EQ.3) THEN
            WRITE (LUPRI,203)
203         FORMAT(' --- EXTERNAL FIELD NO-PAIR RESULT ---')
         ENDIF
         WRITE (LUPRI,100)
  100    FORMAT(8X,'IT',18X,'TOTAL',14X,'POTENTIAL',16X,'KINETIC',
     &          7X,'VIRIAL RELATION')
         RETURN
      ENDIF

      POT=0.D0
      POTN=0.D0
      CIN=0.D0
      IPQ=0
C
      DO I=1,N1T
        ENERG=DT(I)
        DO J=1,I
          K=(I*(I-1)/2)+J
          IPQ=IPQ+1
          TERM=ENERG*P(IPQ)*DT(J)-DOS(I)*Q(IPQ)*DOS(J)
          IF (I.EQ.J) TERM=0.5D0*TERM
          POT=POT+TERM
        END DO
        POTN=POTN+U(I)*ENERG
        CIN=CIN+T(I)*ENERG
      ENDDO

      POT=POT-ZN*POTN
      ENERG=CIN+POT
      VIR=POT/CIN
C
CMI   if (relscf_verbose) WRITE (6,65) NITSCF,ENERG,POT,CIN,VIR
CMI   if (IPR_RELSCF.GE.2) THEN
      if (IPR_RELSCF.GE.0) THEN
CMI      WRITE (LUPRI,65) NITSCF,ENERG,POT,CIN,VIR
         WRITE (LUPRI,'(5X,I5,4(3X,D20.10))') NITSCF,ENERG,POT,CIN,VIR
CMI65 FORMAT(5X,I5,4(3X,D20.10))
      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      DOUBLE PRECISION FUNCTION
     *PVP(CHARGE,AL,BE,L1,M1,N1,L2,M2,N2)
C *** PVP V 1.0 - 19.1.84 - BERND HESS
C *** PVP V 1.1 -  5.6.86 - BERND HESS
#include "implicit.h"
#include "relscf_crelop.h"
      INTEGER IS1(3),IS2(3)
C
C    CALCULATE ANGULAR AND RADIAL PART
C
      II=L1+L2
      JJ=M1+M2
      KK=N1+N2
CMI   ... ftnchek complains that LAMBDA  never set
      LAMBDA = 0
      IMAX=II+JJ+KK+LAMBDA+3
      IF (IMAX.LE.20) GOTO 2
C
C    ERROR BRANCH: ANGULAR MOMENTUM  > MAXIMUM GIVEN BY ARRAY GAM
C
1001  WRITE (6,1002) L1,M1,N1,L2,M2,N2,LAMBDA
1002  FORMAT(' ILLEGAL ANGULAR MOMENTUM (PVP)'/,
     *       ' L1,M1,N1,L2,M2,N2,LAMBDA PRINTED'/,1X,7I5)
      STOP 1002
C
C    COMPUTE INTEGRAL OVER DERIVATIVE OF THE FUNCTIONS
C
2     IS1(1)=L1
      IS1(2)=M1
      IS1(3)=N1
      IS2(1)=L2
      IS2(2)=M2
      IS2(3)=N2
      SUM=DER(1,IS1,IS2,AL,BE)+DER(2,IS1,IS2,AL,BE)+
     \    DER(3,IS1,IS2,AL,BE)
      SUM=SUM*CHARGE
      PVP=SUM
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      DOUBLE PRECISION FUNCTION RCPG(XVAL,ODD,EVEN)
      REAL*8 XVAL,ODD,EVEN
C
C     PROCEDURE TO CALCULATE RECIPROCAL OF GAMMA FUNCTION
C
      REAL*8 X,ALFA,BETA,X2,B(12)
      INTEGER I
      DATA B /-.28387 65422 76024D0, -.07685 28408 44786D0,
     *         .00170 63050 71096D0,  .00127 19271 36655D0,
     *         .00007 63095 97586D0, -.00000 49717 36704D0,
     *        -.00000 08659 20800D0, -.00000 00331 26120D0,
     *         .00000 00017 45136D0,  .00000 00002 42310D0,
     *         .00000 00000 09161D0, -.00000 00000 00170D0/
      X=XVAL
      X2=X*X*8.D0
      ALFA=-.00000 00000 00001D0
      BETA=0.D0
      I=12
1     BETA=-(ALFA*2.D0+BETA)
      ALFA=-BETA*X2-ALFA+B(I)
      I=I-2
      IF (I.GE.2) GOTO 1
      EVEN=(0.5D0*BETA+ALFA)*X2-ALFA+0.92187 02936 50453D0
      ALFA=-0.00000 00000 00034D0
      BETA=0.D0
      I=11
2     BETA=-(ALFA*2.D0+BETA)
      ALFA=-BETA*X2-ALFA+B(I)
      I=I-2
      IF (I.GE.1) GOTO 2
      ODD=(ALFA+BETA)*2.D0
      RCPG=ODD*X+EVEN
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RELOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      RELOP - PROGRAM TO CALCULATE RELATIVISTIC ONE-ELECTRON
C      OPERATORS IN A BASIS SET OF GAUSSIAN FUNCTIONS
C      V 1.0 - 12.3.86 - BERND HESS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "relscf_crelop.h"
C
C     SUBROUTINE RELOP INITIALIZES THE COMMON BLOCK USED BY
C     THE RELOP PACKAGE
C     V 1.0 - 12.3.86 - BERND HESS
C
#include "amfi_if.h"
c     logical relscf_verbose
c     relscf_verbose = .false.

CMI   if (relscf_verbose) WRITE (6,100)
      if (IPR_RELSCF.GE.2)  THEN
         WRITE (LUPRI,100)
100      FORMAT(/,
     X' ****** RELATIVISTIC OPERATORS V 1.0 - BERND HESS ******'
     X//)
      ENDIF

      PI=4.D0*DATAN(1.D0)
      ZWP=2.D0*PI
      ZWPH32=ZWP**1.5D0
      ZWPH12=DSQRT(ZWP)
      SQPI=DSQRT(PI)

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
CMI   ... rough value of the speed of light !!!!
      VELIT=137.036D0
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      PREA=1.D0/(VELIT*VELIT)
      CSQ=VELIT*VELIT
      FAK(1)=1.D0
      DO I=2,26
        FAK(I)=FAK(I-1)*DFLOAT(I-1)
      END DO
C
C     BINOMIALKOEFFIZIENTEN
C
      IMAX=20
      BCO(1)=1.D0
      IBIAS=1
      JBIAS=1
      K=IMAX-1
      DO I=1,K
        ADD=0.D0
        DO J=1,I
          JBIAS=JBIAS+1
          BCO(JBIAS)=ADD+BCO(IBIAS)
          ADD=BCO(IBIAS)
          IBIAS=IBIAS+1
        END DO
        JBIAS=JBIAS+1
        BCO(JBIAS)=1.D0
      END DO
C
      DO N=1,20
        GA(N)=GAM(N-1)
      ENDDO

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      SUBROUTINE SCHMID(C,N,NBAS1,MXBF,S,SCR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     V 1.0 - 17.7.87 - BERND HESS
C
C     SCHMIDT-ORTHOGONALIZATION OF N-TH VECTOR IN C TO THE N-1 VECTORS
C     WITH LOWER INDEX. THE N-TH VECTOR IS ASSUMED TO BE NORMALIZED
C
C     C         VECTOR MATRIX
C     N         NUMBER OF LINEAR INDEPENDENT VECTORS. THE N-TH VECTOR
C               WILL BE ORTHOGONALIZED TO THE REST
C     NBAS1     VECTOR DIMENSION
C     MXBF      MAXIMAL NUMBER OF BASIS SETS IN A SYMMETRY TYPE. ROW
C               DIMENSION OF C ARRAY
C     S         OVERLAP MATRIX FOR AO BASIS (LOWER TRIANGULAR)
C     SCR       SCRATCH ARRAY OF DIMENSION N-1.
C               ON OUTPUT, THIS ARRAY HOLDS OVERLAPS OF THE N-TH
C               VECTOR (BEFORE ORTHOGONALIZATION) WITH THE N-1 VECTORS
C               OF LOWER INDEX
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "relscf_v10.h"
C     DIMENSION C(MXBF,1),S(*),SCR(*)
      DIMENSION C(MXBF,*),S(*),SCR(*)
      IF (N.LE.1) RETURN
      NM1=N-1
C
C     CALCLUATE OVERLAP
C
      SUM=0.D0
      DO I=1,NM1
        SU=0.D0
        IG=0
        DO J=1,NBAS1
          DO K=1,J
            IG=IG+1

            if (IG.GT.LPQ) call quit('IG > LPQ !')

            SU=SU+C(J,N)*S(IG)*C(K,I)
            IF (J.EQ.K) GOTO 3
            SU=SU+C(K,N)*S(IG)*C(J,I)
 3          CONTINUE
          END DO
        END DO
        SCR(I)=SU
        SUM=SUM+SU*SU
      END DO
      SUM=1.D0/DSQRT(1.D0-SUM)
      DO I=1,NBAS1
        SU=C(I,N)
        DO J=1,NM1
          SU=SU-SCR(J)*C(I,J)
        END DO
        C(I,N)=SU*SUM
      END DO
      RETURN
      END

      DOUBLE PRECISION FUNCTION SECONDX()
      SECONDX=0.
      RETURN
      END
      DOUBLE PRECISION FUNCTION SINHX(XVAL)
      REAL*8 XVAL
      REAL*8 AX,Y,S,X
      X=XVAL
      AX=DABS(X)
      IF (AX.GE.0.3D0) GOTO 1
      IF (AX.GE.0.1D0) GOTO 2
      Y=X*X
      GOTO 3
2     Y=X*X/9.D0
3     CONTINUE
      X=(((1.D0/5040D0*Y+1.D0/120.D0)*Y+1.D0/6.D0)*Y+1.D0)*X
      IF (AX.GE.0.1D0) GOTO 5
      SINHX=X
      RETURN
5     SINHX=X*(1.D0+4.D0*X*X/27.D0)
      RETURN
1     AX=DEXP(AX)
      SINHX=DSIGN( (0.5D0*(AX-1.D0/AX)),X )
      RETURN
      END

      FUNCTION SINY(NY,NAB,NCD,ZAB,ZCD)
#include "implicit.h"
#include "relscf_dim.h"
C
C     CALCULATES THE INTEGRAL I(NY) NEEDED FOR THE SUPERMATRICES P AND Q
C     WHEN A SLATER BASIS SET IS USED
C     ******************************************************************
      TABCD=ZAB/ZCD
      NABNY=NAB-NY-1
      NCDNY=NCD+NY
      VNAB=FATT(NABNY+1)/ZAB**(NABNY+1)
      VNCD=FATT(NCDNY+1)/ZCD**(NCDNY+1)
      TERM1=0.D0
      NALF=NABNY+1
      DO I=1,NALF
        TERM1=TERM1+FATT(NABNY+NCDNY+2)*TABCD**(I-1)/(FATT(I)*
     &              FATT(NABNY+NCDNY+3-I))
      END DO
    1 CONTINUE
      CABCD=TERM1*(TABCD+1.D0)**(-NABNY-NCDNY-1)
      SINY=VNAB*VNCD*CABCD
      RETURN
      END

      SUBROUTINE SOG(N,SS,SINV,P,G,A1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE TO CALCULATE TRANSFORMATION TO SCHMIDT-
C     ORTHOGONALIZED BASIS.
C     N              DIMENSION OF MATRICES. ISIZE=N*(N+1)/2
C     SS(ISIZE)      ORIGINAL OVERLAP MATRIX (LOWER TRIANGULAR)
C                    WILL NOT BE DESTROYED
C     P (ISIZE)      OUTPUT TRANSFORMATION MATRIX
C     G (ISIZE)      SCRATCH
C     A1(N)          SCRATCH
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
CMI   DIMENSION SS(2),P(2),G(2),A1(2),SINV(N,N)
      DIMENSION SS(*),P(*),G(*),A1(*),SINV(N,N)
      JL=0
      IQ=0
      DO J=1,N
        IL=JL
        JQ=IQ
        S1KK=SS(IQ+J)
        G(IL+J)=1.D0
        IF(J.EQ.1)GO TO 341
        J1=J-1
        JL=0
        DO K=1,J1
          LG=JQ
          ETOT=0.D0
          DO L=1,K
            LG=LG+1
            JL=JL+1
            ETOT=ETOT+SS(LG)*G(JL)
          END DO
          S1KK=S1KK-ETOT*ETOT

              A1(K)=ETOT
        END DO
        IF=1
        JL=IL
        DO K=1,J1
          SUM=0.D0
          JL=JL+1
          IF=IF+K-1
          IH=IF
          DO L=K,J1
            IH=IH+L-1
            SUM=SUM+A1(L)*G(IH)
          END DO
          G(JL)=-SUM
        END DO
  341   S1KK=1.D0/DSQRT(S1KK)
        JL=IL
        DO K=1,J
          JL=JL+1
          IQ=IQ+1
          G(JL)=G(JL)*S1KK
          P(IQ)=G(JL)
        END DO
      END DO
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          SINV(I,J)=0.D0
          SINV(J,I)=P(IJ)
        END DO
      END DO
      RETURN
      END
C
C      SQROPY - PROGRAM TO CALCULATE RELATIVISTIC ONE-ELECTRON
C      OPERATORS IN A BASIS SET OF GAUSSIAN FUNCTIONS
C      V 1.0 - 12.3.86 - BERND HESS
C      V 1.1 - 21.5.87 - BERND HESS - NORMALIZED FUNCTIONS
C      SUBROUTINE "RELOP" MUST BE CALLED TO INITITALIZE COMMON
C      BLOCK BEFORE FIRST EXECUTION OF THIS ROUTINE
C      INTEGRALS ARE CALCULATED BETWEEN   NORMALIZED FUNCTIONS

      DOUBLE PRECISION FUNCTION
     *SQROPY(AA,AB,LA,MA,NA,LB,MB,NB)
#include "implicit.h"
#include "priunit.h"
#include "relscf_crelop.h"
C
C     SUBROUTINE SQROPY CALCULATES MATRIXELEMENT BETWEEN FUNCTION
C     XA**LA * YA**MA * ZA**NA * EXP(-AA*RA**2) AND
C     XB**LB * YB**MB * ZB**NB * EXP(-AB*RB**2)
C     WITH XA=X ,YA=Y ZA=Z
C     V 1.1 - 1.6.86 - BERND HESS
C
C
C     WE ALSO NEED OVERLAP. CALCULATE IT IN POSITION SPACE
C
      S=0.D0
      IMAX=LA+MA+NA+LB+MB+NB+3
      IF (IMAX.GT.20) GOTO 1100
      II=LA+LB
      JJ=MA+MB
      KK=NA+NB
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
C
C     AND USE ANGULAR PART TO DETERMINE IF IT VANISHES BY SYMMETRY
C
      IF (ANG.EQ.0.D0) THEN
       SQROPY=0.D0
       RETURN
      ENDIF
      EX=-0.5D0*DFLOAT(IMAX)
      OVL=0.5D0*ANG*GA(IMAX)*((AA+AB)**EX)
C
C     CALCULATE SQR INTEGRAL IN MOMENTUM SPACE
C
      EA=1.D0/(4.D0*AA)
      EB=1.D0/(4.D0*AB)
      DELTA=EA+EB
      LLA=LA+MA+NA
      LLB=LB+MB+NB
      C=((2.D0*AA)**LLA * (2.D0*AB)**LLB)
      C=C*(4.D0*AA*AB)**1.5D0
      C=1.D0/C
      LAX=LA/2+1
      MAX=MA/2+1
      NAX=NA/2+1
      LBX=LB/2+1
      MBX=MB/2+1
      NBX=NB/2+1
C
C     COMPUTE VALUE OF BESSEL FUNCTION
C
      ARG=0.5D0*DELTA*CSQ
      CALL BESSKA(0.D0,ARG,X0,X1)
C
C     COMPUTE ARGUMENT OF CONTINUED FRACTION
C
      ARG=1.D0/(CSQ*DELTA)
C
C     AND RUN THE TEDIOUS SUM OVER SIMPLE FOURIER INTEGRALS ...
C
      DO IA=1,LAX
        D1=DCOF(AA,LA,IA-1)
        DO JA=1,MAX
          D2=D1*DCOF(AA,MA,JA-1)
          DO KA=1,NAX
            D3=D2*DCOF(AA,NA,KA-1)
            DO IB=1,LBX
              D4=D3*DCOF(AB,LB,IB-1)
              DO JB=1,MBX
                D5=D4*DCOF(AB,MB,JB-1)
                DO KB=1,NBX
                  D6=D5*DCOF(AB,NB,KB-1)
                  II=LA-2*IA+LB-2*IB+4
                  JJ=MA-2*JA+MB-2*JB+4
                  KK=NA-2*KA+NB-2*KB+4
C
C     ANGULAR PART OF FOURIER INTEGRAL
C
                  IMAX=II+JJ+KK+3
                  ANG=THETA(II+JJ,KK)*PHI(JJ,II)
                  IF (ANG.EQ.0.D0) GOTO 6
C
C     RADIAL PART
C
                  U=X1*D6
                  IF (IMAX.EQ.3) GOTO 60
                  JMAX=IMAX-1
                  DO I=3,JMAX,2
                    BET=0.5D0*DFLOAT(I)
                    CALL KBR(BET,ARG,R)
                    U=U*R*BET/DELTA
                  END DO
60                CONTINUE
                  U=0.25D0*VELIT*U/DELTA
                  S=S+U*ANG
 6                CONTINUE
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
      IF (MOD(LLA-LLB,4).EQ.2) S=-S
      S=S*C
      S=CSQ*(S-OVL)
1000  CONTINUE
C
C     NORMALIZATION
C
      II=LA+LA
      JJ=MA+MA
      KK=NA+NA
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
      EX=-0.5D0*DFLOAT(II+JJ+KK+3)
      OV1=0.5D0*ANG*GA(II+JJ+KK+3)*((AA+AA)**EX)
C#    WRITE (6,*) ' SRQ  OV1',LA,MA,NA,AA,ANG,DSQRT(1/OV1)
      II=LB+LB
      JJ=MB+MB
      KK=NB+NB
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
      EX=-0.5D0*DFLOAT(II+JJ+KK+3)
      OV2=0.5D0*ANG*GA(II+JJ+KK+3)*((AB+AB)**EX)
C#    WRITE (6,*) ' SRQ  OV2',LB,MB,NB,AB,ANG,DSQRT(1/OV2)
      SQROPY=S/DSQRT(OV1*OV2)
      RETURN
C
C     ANGULAR MOMENTUM TOO LARGE
C
1100  WRITE (6,1101) LA,MA,NA,LB,MB,NB,AA,AB
1101  FORMAT(' *** ANGULA MOMENTUM TOO LARGE ***'/,
     *6I3,3X,2D10.10)
      STOP 1100
      END

      DOUBLE PRECISION FUNCTION THETA(M,N)
#include "implicit.h"
#include "relscf_crelop.h"
C
C     INTEGRATION OVER THETA. INCLUDES A FACTOR SIN(TH)
C     FOR THE VOLUME ELEMENT
C
      IF (MOD(N,2).EQ.1) GOTO 10
      THETA=GA(M+2)*GA(N+1)/GA(M+N+3)
      RETURN
10    THETA=0.D0
      RETURN
      END

      SUBROUTINE TRSM(A,B,C,N,H,W)
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(*),B(N,N),C(*),H(N,N),W(N,N)
C
C     TRANSFORM SYMMETRIC MATRIX A BY UNITARY TRANSFORMATION
C     IN B. RESULT IS IN C
C
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          C(IJ)=0.D0
          W(I,J)=A(IJ)
          W(J,I)=A(IJ)
          H(I,J)=0.D0
          H(J,I)=0.D0
        END DO
      END DO
      DO I=1,N
        DO L=1,N
          DO K=1,N
            H(I,L)=B(K,I)*W(K,L)+H(I,L)
          END DO
        END DO
      END DO
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          DO L=1,N
            C(IJ)=H(I,L)*B(L,J)+C(IJ)
          END DO
        END DO
      END DO
      RETURN
      END

      SUBROUTINE TRSMT(A,B,C,N,H,W)
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(*),B(N,N),C(*),H(N,N),W(N,N)
C
C     B*A*BT
C
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          C(IJ)=0.D0
          W(I,J)=A(IJ)
          W(J,I)=A(IJ)
          H(I,J)=0.D0
          H(J,I)=0.D0
        END DO
      END DO
      DO I=1,N
        DO L=1,N
          DO K=1,N
            H(I,L)=B(I,K)*W(K,L)+H(I,L)
          END DO
        END DO
      END DO
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          DO L=1,N
            C(IJ)=H(I,L)*B(J,L)+C(IJ)
          END DO
        END DO
      END DO
      RETURN
      END
      subroutine cpu
      return
      end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZEXTRA
     *(P,Q,S,U,T,DT,DOS,PCAP,QCAP,FC,FO,CM1SCF,CM2SCF,
     *    CSAV1,CSAV2,GMIN,SMIN,THRE,EPS,EC,CM1,Z1,Z2,DELZ,
     *    C,ITRAN,CTRAN,NROW,PQN,ZETA,NBVAR,ZETINC,IORCN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "relscf_dim.h"
C     DIMENSION P(2),Q(2),S(2),U(2),T(2),DT(2),DOS(2),PCAP(2),
      DIMENSION P(*),Q(*),S(*),U(*),T(*),DT(*),DOS(*),PCAP(*),
CMI  *QCAP(2),FC(2),FO(2),CM1SCF(MXBF,1),CM2SCF(MXBF,1),
     *QCAP(*),FC(*),FO(*),CM1SCF(MXBF,1),CM2SCF(MXBF,1),
CMI  *CSAV1(MXBF,1),CSAV2(MXBF,1),GMIN(2),SMIN(2),THRE(*),EPS(2),
     *CSAV1(MXBF,1),CSAV2(MXBF,1),GMIN(*),SMIN(*),THRE(*),EPS(*),
CMI  *EC(2),CM1(2),Z1(2),Z2(2),DELZ(2),
     *EC(*),CM1(*),Z1(*),Z2(*),DELZ(*),
CMI  *C(MXBF,1),CTRAN(2),PQN(*),ZETA(*),ZETINC(*),
     *C(MXBF,*),CTRAN(*),PQN(*),ZETA(*),ZETINC(*),
CMI  *ITRAN(2),NROW(*),NBVAR(*),IORCN(*)
     *ITRAN(*),NROW(*),NBVAR(*),IORCN(*)
      COMMON /CLOCKX/ CPU,CLK
C
C     PROGRAM FOR EXTRAPOLATION OF ORBITAL EXPONENTS
C     ******************************************************************
C
      IF(MXVAR-MVAR) 2,2,1
    1 SEN1=ENERG
      DO I=1,NVAR
        NBVAR1=NBVAR(I)
        Z1(I)=ZETA(NBVAR1)
      END DO
      RETURN
    2 IF(IEX-2) 3,5,8
    3 SUM=0
      DO I=1,NVAR
        B=Z1(I)
        NBVAR1=NBVAR(I)
        A=ZETA(NBVAR1)-B
        DELZ(I)=A
        SUM=SUM+(A/B)**2
      END DO
      SUM=DSQRT(SUM)
      IF(SUM-1.0D-06) 18,18,4
   18 WRITE(6,19) SUM
   19 FORMAT('0',5X,'RELATIVE CHANGE IN ORBITAL EXPONENTS',D8.2,
     *'LESS THAN 1.0D-06')
   27 NEXTRA=1
      GO TO 8
    4 SEN2=ENERG
      TEST=(SEN2-SEN1)/SEN1
      IF(TEST-0.1D-09) 24,24,25
   24 WRITE(6,26) TEST
   26 FORMAT('0',5X,'RELATIVE ENERGY DIFFERENCE',D8.2,
     *'LESS THAN 0.1D1-09')
      GO TO 27
   25 DO 13 I=1,NVAR
      NBVAR1=NBVAR(I)
      Z2(I)=ZETA(NBVAR1)
   13 ZETA(NBVAR1)=Z2(I)+DELZ(I)
      IEX=2
      WRITE(6,30) (Z2(I),I=1,NVAR)
   30 FORMAT('0',6D17.6)
      CPU=SECONDX()
      WRITE(6,31) SEN2,VIR,CPU
   31 FORMAT('0',5X,'T.E.=',D17.10,5X,'V.T.=',D17.10,
     *5X,'CPU=',F14.2)
      RETURN
    5 IF(ENERG-SEN2) 6,7,7
    6 SEN1=SEN2
      GO TO 4
    7 UBER=ENERG-SEN1
      UNTER=2*(ENERG+SEN1-2*SEN2)
      FA=UBER/UNTER
      IEX=3
      DO I=1,NVAR
        NBVAR1=NBVAR(I)
        A=Z2(I)
        B=DELZ(I)
        Z1(I)=A+B
        ZETA(NBVAR1)=A-FA*B
      END DO
      WRITE(6,30) (Z1(I),I=1,NVAR)
      WRITE(6,31) ENERG,VIR
      GO TO 10
    8 WRITE(6,20) NVAR
   20 FORMAT('0',5X,'ORBITAL EXPONENTS OPTIMIZED IN',I3,' DIMENSIONS')
      WRITE(6,21)
   21 FORMAT('0',5X,'ORBITAL EXPONENTS')
      WRITE(6,22) (ZETA(I),I=1,NBAST)
   22 FORMAT(1X,6D17.6)
      CPU=SECONDX()
      WRITE(6,31) ENERG,VIR,CPU
      IEX=0
      NEXTRA=NEXTRA-1
   10 RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE TEIG(P,Q,NSYM,NBAS,NOSH,ZETA,ajmn,U,OCCUP,ndimpq)
#include "implicit.h"
#include "priunit.h"
      LOGICAL OPEN,KLNEMN
C.......................................................................
C
C     TWO-ELECTRON INTEGRAL ROUTINE FOR S,P,D, AND F FUNCTIONS.
C.......................................................................
      LOGICAL RELFLG

      dimension p(*),q(*),occup(*)
      DIMENSION NBAS(*),NOSH(*),ZETA(*),AJMN(*),U(*),NCCUP(5)
      DIMENSION TWOPOW(14), PIFAC(14)

C.......................................................................
C
C     ANGULAR FACTORS FOR EXCHANGE INTEGRALS. OBTAINED AS SUM OF SQUARES
C     OF SLATER COEFFICIENTS C(KAPPA;L1,M1;L2,M2) DIVIDED BY
C     2*(2*L1+1)*(2*L2+1) FOR ANY KAPPA, L1, AND L2.
C.......................................................................
      DATA SS0,SP1,PP0,PP2,SD2,PD1,PD3,DD0,DD2,DD4,
     X     SF3,PF2,PF4,DF1,DF3,DF5,FF0,FF2,FF4,FF6
     X/ .50000000000D+00, .16666666667D+00, .16666666667D+00,
     X  .66666666667D-01, .10000000000D+00, .66666666667D-01,
     X  .42857142857D-01, .10000000000D+00, .28571428571D-01,
     X  .28571428571D-01, .71428571429D-01, .42857142857D-01,
     X  .31746031746D-01, .42857142857D-01, .19047619048D-01,
     X  .21645021645D-01, .71428571429D-01, .19047619048D-01,
     X  .12987012987D-01, .16650016650D-01/
      F0POL(A,B) = 3*(16*A**6 + 104*A**5*B + 286*A**4*B**2 +
     X             429*(A*B)**3 + 286*A**2*B**4 + 104*A*B**5 + 16*B**6)
      F2POL(A,B) = 8*A**4 + 52*A**3*B + 143*(A*B)**2 + 52*A*B**3
     X           + 8*B**4
      DF1POL(A,B) = 8*A**4 + 44*A**3*B + 99*(A*B)**2 + 44*A*B**3
     X            + 8*B**4
      DF3POL(A,B) = 2*A**2 + 11*A*B + 2*B**2
C.......................................................................
C     TWO-ELECTRON INTEGRAL ROUTINE FOR LCGO ATOM SCF.
C     RESTRICTED TO PRINCIPAL QUANTUM NUMBERS 1,2 AND 3 FOR RESPECTIVELY
C     S,P, AND D ORBITALS.
C.......................................................................

      DO i=1,5
            nccup(i)=occup(i) + .01
      END DO
      CALL SETABJ(AJMN,NSYM,NOSH,NCCUP)
      PI=3.1415926536
      ON6 = 1.0/6.0D00
      ON15 = 1.0/15.D00
      ON70 = 1.0/70.D00
      ON35 = 2*ON70
      PIFAC(1) = DSQRT(PI)
      TWOPOW(1) = 1
      DO I = 2,14
         TWOPOW(I) = TWOPOW(I-1)*0.5D00
         PIFAC(I) = PIFAC(I-1)*(2*I-1)
      END DO
C.......................................................................
C     
C     THIS PART SETS UP THE COEFFICIENTS LAMBDA,P,Q AND MU,R,S.
C.......................................................................
      IPQ=0
      J=0
      NSTEP1=0
      KMX=0
      FACTKL = 1
      KL = 0
      DO I=1,NSYM
         PRFAC1 = TWOPOW(I+1)*FACTKL
         KIN=KMX+1
         KMX=KIN+NBAS(I)-1
         DO K=KIN,KMX
            ZP = ZETA(K)
            DO L=KIN,K
               KL = KL + 1
               ZQ = ZETA(L)
               ZPQ = ZP + ZQ
               PRFAC2 = PRFAC1*U(KL)
               XFAC1 = PRFAC2*ZPQ**I*PIFAC(I)
               NSTEP2=0
               MMX=0
               FACTMN = 1.0
               MN = 0
               DO IM=1,I
                  OPEN = (NOSH(I).NE.0.AND.NOSH(IM).NE.0)
                  PRFAC3 = PRFAC2*PIFAC(IM)*TWOPOW(IM)*FACTMN
                  XFAC2 = XFAC1*FACTMN*TWOPOW(IM+1)
                  MIN=MMX+1
                  MMX=MIN-1+NBAS(IM)
                  MMXP = MMX
                  IF(IM.EQ.I) MMXP = K
                  DO M=MIN,MMXP
                     ZR = ZETA(M)
                     ZPR = ZP+ZR
                     ZQR = ZQ+ZR
                     NMX = L
                     IF(M.LT.K)NMX = M
                     DO N=MIN,NMX
                        IPQ=IPQ+1
                        MN = MN + 1
                        KLNEMN = (KL.NE.MN)
                        J=J+1
C.......................................................................
C     
C     J IS THE NUMBER LABEL OF THE MATRIX ELEMENTS TO BE CALCULATED
C     I=LAMBDA+1,K=P,L=Q,IM=MU+1,M=R,N=S
C.......................................................................
                        ZS = ZETA(N)
                        ZRS = ZR + ZS
                        ZQS = ZQ + ZS
                        ZPS = ZP + ZS
                        ZPQRS = ZPQ+ZRS
                        ZPQRS2 = 2*ZPQRS**2
                        ZPRZQS = ZPR*ZQS
                        ZPSZQR = ZPS*ZQR
                        XTERM = (1.0D00/DSQRT(ZPQRS))**
     &                       (2*(I+IM)-3)
                        PRFAC4 = PRFAC3*U(MN)*XTERM
                        XFAC3 = XFAC2*U(MN)*XTERM
                        XFAC11 = XFAC3*(ZRS/ZPRZQS)**IM
                        XFAC21 = XFAC3*(ZRS/ZPSZQR)**IM
                        XFSUM = XFAC11 + XFAC21
                        NTEST  = I*(I-1)/2 + IM
                        GO TO (110,120,130,140,150,160,
     &                       400,410,420,430,170),NTEST
C.......................................................................
C     
C     I=1,IM=1,(SS)-LOOP. X0=J0(SS),Y0=K0(SS)
C.......................................................................
 110                    X0 = PRFAC4
                        Y0 = XFSUM
                        PJ=X0-Y0*SS0
                        QJ=-AJMN(1)*Y0
                        GO TO 440
C.......................................................................
C     
C     I=2,IM=1,(SP)-LOOP. X0=J0(SP),Y1=K1(SP)
C.......................................................................
 120                    X0 = PRFAC4*(3*ZPQ+2*ZRS)
                        Y1 = XFSUM
                        PJ=X0-Y1*SP1
                        QJ=-AJMN(2)*Y1
                        GO TO 440
C.......................................................................
C     
C     I=2,IM=2,(PP)-LOOP. X0=J0(PP),Y0=K0(PP),Y2=K2(PP)
C.......................................................................
 130                    X0 = PRFAC4*(ZPQRS2+ZPQ*ZRS)
                        Y0 = XFSUM*ZPQRS2
                        XFSUM = XFAC11*ZPR*ZQS+
     &                       XFAC21*ZPS*ZQR
                        Y0 = Y0 + XFSUM
                        Y2 = XFSUM*5
                        PJ=X0-Y0*PP0-Y2*PP2
                        QJ=-(AJMN(3)*Y0+AJMN(4)*Y2)
                        GO TO 440
C.......................................................................
C     
C     I=3,IM=1,(SD)-LOOP. X0=J0(SD),Y2=K2(SD)
C.......................................................................
 140                    X0=PRFAC4*
     &                       (15*ZPQ**2+20*ZPQ*ZRS+8*ZRS**2)
                        Y2 = XFSUM
                        PJ=X0-Y2*SD2
                        QJ=-AJMN(5)*Y2
                        GO TO 440
C.......................................................................
C     
C     I=3,IM=2,(PD)-LOOP. X0=J0(PD),Y1=K1(PD),Y3=K3(PD),X2=J2(PD)
C.......................................................................
 150                    X0 = PRFAC4*(10*ZPQ**3+35*ZPQ**2*ZRS
     X                       +28*ZPQ*ZRS**2+8*ZRS**3)
                        Y1 = XFSUM*ZPQRS2
                        XFSUM = XFAC11*ZPR*ZQS+XFAC21*ZPS*ZQR
                        Y1 = Y1 + 3*XFSUM
                        Y3 = 7*XFSUM
                        PJ=X0-Y1*PD1-Y3*PD3
                        QJ=0.
                        IF (.NOT.OPEN) GO TO 440
                        X2 = PRFAC4*5*ZPQ*ZRS*(7*ZPQ+2*ZRS)
                        QJ=AJMN(21)*X2-(AJMN(6)*Y1+AJMN(7)*Y3)
                        GO TO 440
C.......................................................................
C     
C     I=3,IM=3,(DD)-LOOP. X0=J0(DD),Y0=K0(DD),Y2=K2(DD),Y4=K4(DD)
C.......................................................................
 160                    ZPRZQS = ZPR*ZQS
                        ZPSZQR = ZPS*ZQR
                        X0 = PRFAC4*((ZPQRS2+ZPQ*ZRS)*ZPQRS2*2
     X                       +7*(ZPQ*ZRS)**2)
                        Y01= XFAC11*((ZPQRS2+ZPRZQS)*ZPQRS2*2
     X                       +7*ZPRZQS**2)
                        Y02= XFAC21*((ZPQRS2+ZPSZQR)*ZPQRS2*2
     X                       +7*ZPSZQR**2)
                        Y0 = Y01+Y02
                        XFAC11 = XFAC11*7*ZPRZQS
                        XFAC21 = XFAC21*7*ZPSZQR
                        Y21 = XFAC11*(ZPQRS2 + 5*ZPRZQS)
                        Y22 = XFAC21*(ZPQRS2 + 5*ZPSZQR)
                        Y2 = Y21+Y22
                        Y4 = (XFAC11*ZPRZQS + XFAC21*ZPSZQR)*9
                        PJ=X0-Y0*DD0-(Y2+Y4)*DD2
                        QJ=-AJMN(8)*Y0-AJMN(9)*Y2-AJMN(10)*Y4
                        GO TO 440
C.......................................................................
C     
C     I=4,IM=1,(SF)-LOOP. X0=J0(SF),Y3=K3(SF)
C.......................................................................
 400                    X0 = PRFAC4*(105*ZPQ**3 + 210*ZPQ**2*ZRS +
     X                       168*ZPQ*ZRS**2 + 48*ZRS**3)
                        Y3 = XFSUM
                        PJ = X0 - Y3*SF3
                        QJ = -AJMN(11)*Y3
                        GO TO 440
C.......................................................................
C     
C     I=4,IM=2,(PF)-LOOP. X0=J0(PF),Y2=K2(PF),Y4=K4(PF),X2=J2(PF)
C.......................................................................
 410                    X0 =  PRFAC4*(70*ZPQ**4 + 315*ZPQ**3*ZRS +
     X                       378*(ZPQ*ZRS)**2 + 216*ZPQ*ZRS**3 +
     X                       48*ZRS**4)
C     Y2 = XFSUM*ZPQRS2
                        Y2 = XFAC11*(2*ZPR**2+9*ZPR*ZQS+2*ZQS**2)
     X                       + XFAC21*(2*ZPS**2+9*ZPS*ZQR+2*ZQR**2)
                        XFSUM = XFAC11*ZPR*ZQS + XFAC21*ZPS*ZQR
C     Y2 = Y2 + 5*XFSUM
                        Y4 = 9*XFSUM
                        PJ = X0 - Y2*PF2 - Y4*PF4
                        QJ = 0.0
                        IF (.NOT.OPEN) GO TO 440
                        X2 = PRFAC4*5*ZPQ*ZRS*(63*ZPQ**2 + 26*ZPQ*ZRS
     X                       + 8*ZRS**2)
                        QJ = AJMN(22)*X2 - AJMN(12)*Y2 - AJMN(13)*Y4
                        GO TO 440
C.......................................................................
C     
C     I=4,IM=3,(DF)-LOOP. X0=J0(DF),Y1=K1(DF),Y3=K3(DF),Y5=K5(DF)
C     X2=J2(DF),X4=J4(DF).
C.......................................................................
 420                    X0 = PRFAC4*(56*ZPQ**5 + 308*ZPQ**4*ZRS +
     X                       693*ZPQ**3*ZRS**2 + 594*ZPQ**2*ZRS**3 +
     X                       264*ZPQ*ZRS**4 + 48*ZRS**5)
                        Y1=XFAC11*DF1POL(ZPR,ZQS)+XFAC21*DF1POL(ZPS,ZQR)
                        XFAC11 = 9*ZPR*ZQS*XFAC11
                        XFAC21 = 9*ZPS*ZQR*XFAC21
                        Y3=XFAC11*DF3POL(ZPR,ZQS)+XFAC21*DF3POL(ZPS,ZQR)
                        Y5 = 11*(XFAC11*ZPR*ZQS+XFAC21*ZPS*ZQR)
                        PJ = X0 - Y1*DF1- Y3*DF3 - Y5*DF5
                        QJ = 0.0
                        IF(.NOT.OPEN) GO TO 440
                        PRFAC4 = PRFAC4*7*ZPQ*ZRS
                        X2 = PRFAC4*(18*ZPQ**3 + 99*ZPQ**2*ZRS +
     X                       44*ZPQ*ZRS**2 + 8*ZRS**3)
                        X4 = 9*PRFAC4*(11*ZPQ + 2*ZRS)*ZPQ*ZRS
                        QJ = X2*AJMN(23) + X4*AJMN(24) - Y1*AJMN(14)
     X                       - Y3*AJMN(15) - Y5*AJMN(16)
                        GO TO 440
C.......................................................................
C     
C     I=4,IM=4,(FF)-LOOP. X0=J0(FF),Y0=K0(FF),Y2=K2(FF),Y4=K4(FF)
C     Y6=K6(FF)
C.......................................................................
 430                    X0 = PRFAC4*F0POL(ZPQ,ZRS)
                        Y0 = XFAC11*F0POL(ZPR,ZQS)+XFAC21*F0POL(ZPS,ZQR)
                        XFAC11 = XFAC11*ZPR*ZQS*9
                        XFAC21 = XFAC21*ZPS*ZQR*9
                        Y2 = XFAC11*F2POL(ZPR,ZQS)+XFAC21*F2POL(ZPS,ZQR)
                        XFAC11 = XFAC11*ZPR*ZQS*11
                        XFAC21 = XFAC21*ZPS*ZQR*11
                        Y4 = XFAC11*(2*ZPR**2 + 13*ZPR*ZQS + 2*ZQS**2)
     X                       + XFAC21*(2*ZPS**2 + 13*ZPS*ZQR + 2*ZQR**2)
                        Y6 = 13*(XFAC11*ZPR*ZQS + XFAC21*ZPS*ZQR)
                        PJ = X0-Y0*FF0-Y2*FF2-Y4*FF4-Y6*FF6
                        QJ = -Y0*AJMN(17)-Y2*AJMN(18)-Y4*AJMN(19)
     X                       -Y6*AJMN(20)

                        if(ipq.gt.ndimpq) then
                           WRITE(LUPRI,*)
     &                          'relscf/TEIG: ipq, ndimpq:',ipq, ndimpq
                           CALL FLSHFO(LUPRI)
                           CALL QUIT(
     & 'relscf/TEIG: ipq too large ! enlarge ndimpq!')
                        endif

 440                    CONTINUE
                        P(IPQ)=PJ
                        Q(IPQ)=QJ
 170                    CONTINUE
                     END DO
                  END DO
                  FACTMN = FACTMN/IM
               END DO
            END DO
         END DO
         FACTKL = FACTKL/I
      END DO
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SETABJ(AJMN,NSYM,NOSH,NCCUP)
#include "implicit.h"
#include "priunit.h"
C.......................................................................
C
C     SET COUPLING COEFFICIENTS FOR OPEN SHELLS. FOR A TABULATION
C     OF THESE SEE
C                   G. MALLI AND J. P. OLIVE
C                   TECHN. REP. 1962-3, PART TWO,
C                   LAB. OF MOL. STRUCTURE AND SPECTRA,
C                   DEPT. OF PHYSICS, UNIV. OF CHICAGO.
C     (VERY HARD TO GET).
C.......................................................................
CMI   DIMENSION NOSH(2),AJMN(2),NCCUP(2),IANUM(90),IADEN(90)
CMI   DIMENSION NOSH(2),AJMN(24),NCCUP(2),IANUM(90),IADEN(90)
      DIMENSION NOSH(*),AJMN(*),NCCUP(*),IANUM(90),IADEN(90)

      DATA IANUM /-1,-5,2,-2,-1,-1,-2,-1,-1,-1,2,-9,2,2,-4,-26,9,-7,-58,
     X            -34,-3,-1,-1,-1,-2,-2,-2,-1,-1,-3,-174,-34,-1,-13,9,
     X            -1,2,2,-13,4,2,100,-6,-23,-53,3125,-11,-256,-722,
     X            17800,-5,-17,-74,-1325,-9,-344,-1598,-2272,-4,-34,-17,
     X            -850,-1,-4,-2,-100,-3,-17,-17,-425,-5,-344,-1598,
     X            -56800,-2,-17,-296,-53,-3,-256,-722,17800,-1,-23,-53,
     X            3125,-1,4,2,100/
      DATA IADEN /1,3,15,3,15,3,15,6,60,15,375,5,35,35,5,245,245,15,735,
     X            2205,10,20,20,5,35,35,15,45,45,35,12005,12005,20,1960,
     X            3920,45,2835,2835,7,105,77,3003,7,315,2541,99099,21,
     X            2835,22869,891891,14,252,2541,396396,35,7875,63525,
     X            99099,21,945,693,27027,7,105,77,3003,28,840,1232,
     X            24024,63,25515,205821,8027019,35,1575,63525,99099,
     X            77,38115,307461,11990979,42,11340,91476,3567564,91,
     X            17745,13013,507507/

#include "amfi_if.h"
c     logical relscf_verbose
c     relscf_verbose = .false.

      IEXCHG = 20
      DO I = 1,24  ! MI: AJMN has range 1-24 !
        AJMN(I) = 0.0
      END DO
      ICOUNT = 0
      JCOUNT = IEXCHG
      NOFF = 0
      DO I = 1,NSYM
        DO IM = 1,I
          IF(NOSH(I).EQ.0.OR.NOSH(IM).EQ.0) GO TO 30
             KOUNT = 0
             KSTART = NOFF + (NCCUP(I)-1)*I + 1

                 DO K = KSTART,KSTART+IM-1
               KOUNT = KOUNT + 1
               AJMN(ICOUNT+KOUNT) = DFLOAT(IANUM(K))/IADEN(K)

                   END DO
          IF (I.NE.IM.AND.IM.GT.1) THEN
                if (IPR_RELSCF.GE.2)
     &            WRITE(LUPRI,*) ' IN SETABJ I IM', I,IM
          ENDIF
   30     ICOUNT = ICOUNT + IM
          IF(I.NE.IM) JCOUNT = JCOUNT + IM - 1
        END DO
        NOFF = NOFF + (4*I-3)*I
      END DO
      DO I = 1,IEXCHG
        AJMN(I) = 0.5*AJMN(I)
      END DO
      IF (IPR_RELSCF.GE.2) THEN
        WRITE(LUPRI,135)
  135   FORMAT ('0',5X,'AUTOMATIC VECTOR COUPLING COEFFICIENTS K'/ )
         WRITE(LUPRI,200) (2*AJMN(I),I=1,24)
  200    FORMAT(4F18.14)
      ENDIF

  140 FORMAT(5D21.8)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LABAMFI(LUNIT,LABEL)
C***********************************************************************
C
C     Write label LABEL to LUNIT;
C       If label LABEL exists:
C         Seek to label SODLABEL and write from there.
C       Else:
C         Seek to label EOFLABEL, backspace and write from there
C       End If
C
C     Input:
C        LABEL   - the label
C
C     Output:
C
C     Written by S. Knecht - Jan 2009
C     inspired by LABKRMC written by J. Thyssen
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER LABEL*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C        label exists.
C        seek to label SODLABEL and write from there.
         REWIND(LUNIT)
         IF (.NOT. FNDLAB('SODLABEL',LUNIT)) THEN
            CALL QUIT('LABAMFI: wrong structure in file - no SODLABEL')
         END IF
      ELSE
C
C       label does not exist
C       seek to label EOFLABEL and write from there
         REWIND(LUNIT)
         IF (FNDLAB('EOFLABEL',LUNIT))
     &      BACKSPACE LUNIT
C
      END IF
C
C     write label
      CALL NEWLAB(LABEL,LUNIT,LUPRI)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WRTAMFI(LUNIT,LABEL,VECTOR,ISIZE)
C***********************************************************************
C
C     Write VECTOR to LUNIT under label LABEL.
C
C     Input:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Output:
C
C     written by S. Knecht - Jan 2009
C     inspired by WRTKRMC by J. Thyssen
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8
C
C     Write label at end or - if already existing - after SODLABEL
      CALL LABAMFI(LUNIT,LABEL)
C
C     Write vector
C     ... buggy on 32-bit systems 
C     CALL WRITSQ(LUNIT,ISIZE,VECTOR)
      WRITE(LUNIT) VECTOR
      CALL NEWLAB('EOFLABEL',LUNIT,LUPRI)
C
      REWIND LUNIT
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE IWRTAMFI(LUNIT,LABEL,IVECTOR,ISIZE)
C***********************************************************************
C
C     Write VECTOR to LUNIT under label LABEL.
C
C     Input:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Output:
C
C     written by S. Knecht - Jan 2009
C     inspired by IWRTKRMC by J. Thyssen
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION IVECTOR(ISIZE)
      CHARACTER LABEL*8
C
C     Write label at end or - if already existing - after SODLABEL
      CALL LABKRMC(LUNIT,LABEL)
C
C     Write vector
      CALL WRITI(LUNIT,ISIZE,IVECTOR)
      CALL NEWLAB('EOFLABEL',LUNIT,LUPRI)
C
      REWIND LUNIT
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE REAAMFI(LUNIT,LABEL,VECTOR,ISIZE)
C***********************************************************************
C
C     Read VECTOR from LUNIT under label LABEL.
C
C     Input:
C
C     Output:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     written by S. Knecht - Jan 2009
C     inspired by REAKRMC by J. Thyssen
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C        WRITE(LUPRI,*) ' start reading from file for LABEL ',LABEL
C        WRITE(LUPRI,*) ' LUNIT is',LUNIT
C        WRITE(LUPRI,*) ' size of vector is',ISIZE
C        Label exists.
         CALL READT(LUNIT,ISIZE,VECTOR)
C
      ELSE
         WRITE(LUPRI,'(//A/,3A,I3)')
     &        '*** ERROR in REAAMFI ***',
     &        'Label <',LABEL,'> does not exist on unit ',LUNIT
         CALL QUIT('*** ERROR in REAAMFI ***')
      END IF
C
      END
