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

C FILE    : dircnv.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dhfcnv */
      SUBROUTINE DHFCNV(FMO,FAO,DMAT,EVEC,FBUF,TMAT,BMAT,
     &                  ICHKCNV,TINV,WORK,LWORK)
C*****************************************************************************
C     Driver routine for checking convergence etc.
C
C     Written by T.Saue October 1995
C     Last revision 971028-jth: Lot of changes to make one open shell work.
C
C     Miro ILIAS, febr.2007, Tel Aviv - fix for 2c->4c restart
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
C Used from COMMON blocks
C   DCBBAS: N2BBASX
C   DCBORB: N2ORBT
C   DGROUP: NZ
C
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbpsi.h"
#include "dcbham.h"
#include "dcbgen.h"
      LOGICAL OLDFMO, DAMPON, roundv
      DIMENSION FMO(*), EVEC(*),FBUF(*),BMAT(*),TMAT(*),
     &          DMAT(*),TINV(*),WORK(LWORK), FAO(*)
      PARAMETER (N_OLD = 6)
      DIMENSION ERGOLD(N_OLD), EVCOLD(N_OLD)
      SAVE      ERGOLD, EVCOLD
      DATA      ERGOLD/N_OLD*1.0D20/, EVCOLD/N_OLD*1.0D20/

      CALL QENTER('DHFCNV')
C
C     Calculate error vector for this iteration
C     =========================================
C             FMO   - Fock matrix in MO-basis
C             FAO   - Fock matrix in AO-basis (FDao and FVao)
C             DMAT  - Density matrix in AO-basis (DCao and DVao)
C             TMAT  - MO-transformation matrix
C             TINV  - Right index transformed overlap matrix;
C             correponds to the inverse of the transformation matrix
C             FBUF  - work area .ge. NFMAT*N2BBASXQ
C
      IF (ICHKCNV .NE. 0)
     &   CALL DIISER(EVEC,FMO,FAO,DMAT,TINV,FBUF,TMAT,WORK,LWORK)
C
C
C     Read MO Fock matrix from previous cycle into FBUF...
C     Calculate largest change in AO Fock matrix for this iteration
C     =================================================
C
      OLDFMO = ICHKCNV .GT. 0
      IF(OLDFMO) THEN
! PSB+MI change (Aug2014):
        numfmat = max0 (5,mxdiis)
        IREC = MOD(NITER-2,numfmat-1) + 1
        CALL READAC(LUFOCK,N2TMOTQ,FBUF,IREC)
C
        IF(IPRSCF.GE.15) THEN
          DO I = 1,NFSYM
             WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
             CALL HEADER('DHFCNV: old Fock Matrix in FBUF:',-1)
             CALL PRQMAT(FBUF(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                   NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
          END DO
        ENDIF

        CALL FCKCMP(FMO,FBUF,WORK,LWORK)
      ELSE
        ERGOLD(1:N_OLD) = 1.0D20
        EVCOLD(1:N_OLD) = 1.0D20
C       ... reinitialize ERGOLD and EVCOLD /hjaaj
      ENDIF

C     ...and write current MO Fock matrix to file
C     ==============================================
! PSB+MI change/ Aug2014:
        numfmat = max0 (5,mxdiis)
        IREC = MOD(NITER-1,numfmat-1) + 1
        CALL WRTDAC(LUFOCK,N2TMOTQ,FMO,IREC)
c     DEBUG PRINT AFTER write Fock Matrix to file (PBS+MI)
            IF(IPRSCF.GE.15) THEN
               DO I = 1,NFSYM
                  WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
                  CALL HEADER('DHFCNV: MO Fock matrix to file:',-1)
                  CALL PRQMAT(FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                     NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
               END DO
            END IF
C
      ERGVAL = ERGBUF - DHFERG
C
      IF (ICHKCNV .EQ. 0) GOTO 9999
C
C     Check convergence
C     =================
C
!
!     Data for anti-oscillation test:
!     both gradient and energy must be lower compared to the
!     minimum gradient and energy until N_OLD iterations ago.
!
      ERGDIF = DHFERG - ERGOLD(N_OLD)
      EVCDIF = EVCVAL - EVCOLD(N_OLD)
      ERGOLD(N_OLD) = MIN(ERGOLD(N_OLD), ERGOLD(N_OLD-1))
      EVCOLD(N_OLD) = MIN(EVCOLD(N_OLD), EVCOLD(N_OLD-1))
      DO I = N_OLD-1,2,-1
         ERGOLD(I) = ERGOLD(I-1)
         EVCOLD(I) = EVCOLD(I-1)
      END DO
      ERGOLD(1) = DHFERG
      EVCOLD(1) = EVCVAL
!
      IF    (ERGCNV .AND. OLDFMO) THEN
        CONVRG = ERGVAL
      ELSEIF(FCKCNV .AND. OLDFMO) THEN
        CONVRG = FCKVAL
      ELSEIF(EVCCNV) THEN
        CONVRG = EVCVAL
      ELSE
C       no check on convergence possible this time
C       (must be called from PREDHF)
        IF (.NOT.OLDFMO) GOTO 9999
        CALL QUIT('DHFCNV ERROR: no DHF convergence criterium defined!')
      ENDIF
c     ... check if One-center app. should be turned off     /jkp
      IF(ONECAP.AND.(ABS(CONVRG).LT.ONECNV)) THEN
         ONECAP  = .FALSE.
         DOLVC   = .FALSE.
         ONECOFF = .TRUE.
c -->
c automatic conventional evaluation of one-center
c LS and SS integrals for ONECAP models disabled.
c see dirrdn.F
c         IF (.NOT. DIRSET) THEN
c            ILLDIR = 1
c            ISLDIR = 1
c            ISSDIR = 1
c            IDFLAG = ILLDIR+2*ISLDIR+4*ISSDIR
c         ENDIF
c <--
      END IF
c
C
C     Check convergence : we have two tresholds, the first one
C     is the desired treshold, the second one is the allowed
C     upper limit for proceeding with post-HF calculations.
C
CMI ... when continuing the BSS-SCF after the DC-SCF,
CMI      we have to reset the DHFCONV/DHFEXIT values before !
      DHFCONV(1) = ABS(CONVRG).LT.SCFCNV(1)
      DHFCONV(2) = ABS(CONVRG).LT.SCFCNV(2)
      DHFEXIT    = DHFCONV(1) .OR. (KITER.GT.MAXITR)

C
C     If wave function converged or DHFEXIT is true, then
C     return to avoid an extra diagonalization of the Fock matrix.
C
C     IF  (not old FMO matrix) exit,
C     because then neither DIIS nor DAMP possible
C     (OLDFMO can be false when called from PREDHF)
C
      IF ( DHFEXIT .OR. .NOT. OLDFMO) GOTO 9999
C
C
C     Convergence acceleration
C     ========================
C
      CACC = '        '
C     I (jth) want DIIS from iteration no. 2 and no damping.
C     hjj has done non-relativistic tests of DIIS, and it was
C     never found that DIIS slowed things down; in the worst
C     cases there was the same number of iterations.
C
      DIISON = DODIIS.AND.ABS(EVCVAL).LE.DIISTH
C
C     If change in which integral classes included (LL, LS, SS, GT),
C     then we must restart DIIS
C     Likewise if One-center app. is turned off in this iteration     /jkp
C
      IF (INTFLG.NE.INTBUF.OR.ONECOFF.OR..NOT.DIISON) THEN
        IF(ITDIIS.GT.0) THEN
          IF (.NOT.(START2C.AND..NOT.DO4C2C)) THEN
            CLOSE(LUDIIS,STATUS='DELETE')
          ENDIF
         ENDIF
        ITDIIS = 0
      END IF
C
C     DIIS --- Direct Inversion of Iterative Subspace
C     ===============================================
C
      roundv = dabs(ergval/dhferg) .lt. 5.0d-10
c     PSB revision of hjaaj rounding test
c     define roundv = .true. if ergval small and bypass turn off DIIS
c          roundv is .true.
C     ... avoid erroneous predicted energy increase (ERGVAL negative)
C         because of round-off errors /hjaaj
      IF(DIISON) THEN
        ITDIIS = ITDIIS + 1
        ERGDIF = ERGDIF + 1.D-10*DHFERG
!       ... avoid erroneous predicted energy increase (ERGDIF negative)
!           because of round-off errors /hjaaj
        IF (iprscf .GE. 2) THEN
          WRITE(LUPRI,'(/A)') 'DIIS anti-oscillation test data dump'
          WRITE(LUPRI,'(A,1P,10D10.2)')
     &      '- ERGDIF, ERGOLD(:) ',ERGDIF, ERGOLD(1:N_OLD)
          WRITE(LUPRI,'(A,1P,10D10.2)')
     &      '- EVCDIF, EVCOLD(:) ',EVCDIF, EVCOLD(1:N_OLD)
        END IF
        IF (KITER .GE. 8 .AND. ITDIIS .GE. 4 .AND.
     &     EVCDIF*ERGDIF .LT. D0 .and. .not. roundv 
     &     .and. .not. noswit ) THEN
c     DEBUG
        WRITE (6,9005) ITDIIS, EVCDIF*ERGVAL, ROUNDV, NOSWIT
 9005   format ('DEBUG in DHFCNV: ITDIIS=',i3,5x,'grad=',d10.3,5x,
     &         'round,noswit=',2l3/)
c     END DEBUG
C          ... often wild oscillations in the beginning,
C              which may settle down reasonably.
           WRITE(LUPRI,'(/A/A)')
     &        ' DIIS aborted because of the last two iterations',
     &        ' the lowest energy has the largest gradient and'//
     &        ' DIIS minimizes gradient !!!'
           CLOSE(LUDIIS,STATUS='DELETE')
           ITDIIS = 0
           IBSTAT = -1
        ELSE
           CALL DIISDR(BMAT,EVEC,FMO,FBUF,IBSTAT,WORK,LWORK)
        END IF
        IF(IBSTAT.NE.0) THEN
          DIISON = .FALSE.
C
C         If DODAMP Fock matrix from previous cycle has to be reread from
C         previous cycle after being overwritten in DIISDR
C
          IF(DODAMP .AND. IBSTAT.EQ.1) THEN
! NUMFMAT added to make usage consistent with other READAC & WRTDAC
! usage (PSB+MI, Aug2014)
            numfmat = max0 (5,mxdiis)
            IREC = MOD(NITER-2,numfmat-1) + 1
            CALL READAC(LUFOCK,N2TMOTQ,FBUF,IREC)
! DEBUG PRINT AFTER READ Fock Matrix from file (PBS+MI)
            IF(IPRSCF.GE.15) THEN
               DO I = 1,NFSYM
                  WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
                  CALL HEADER('DHFCNV: MO Fock matrix from file:',-1)
                  CALL PRQMAT(FBUF(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                     NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
               END DO
            END IF
          ENDIF
        ELSE
          CACC(1:4) = 'DIIS'
          WRITE(CACC(5:8),'(I4)') NELMBM - 1
C         NELMBM is current dim of B matrix
C         which is number of density/Fock matrics
C         plus one (for the constraint)
          IF(IPRSCF.GE.5) THEN
             DO I = 1,NFSYM
                WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
                CALL HEADER('DHFCNV: MO Fock matrix after DIIS:',-1)
                CALL PRQMAT(FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                   NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
             END DO
          END IF
        ENDIF
      ENDIF
C
C     Damping of total Fock matrix
C     ============================
C
      DAMPON = DODAMP.AND.(.NOT.DIISON).AND.INTFLG.EQ.INTBUF
C.....Do not turn on damping with atomic Huckel before DIIS gets activated
      IF(ATHUCK.AND.NITER.EQ.3) DAMPON=.FALSE.
      IF(DAMPON) THEN
        CACC(1:4) = 'DAMP'
        WRITE(CACC(5:7),'(I3)') INT(DAMPFC*100.0D0)
        CACC(8:8) = '%'
        DNEW = D1 - DAMPFC
        CALL DSCAL(N2TMOTQ,DNEW  ,FMO,1)
        CALL DAXPY(N2TMOTQ,DAMPFC,FBUF,1,FMO,1)
        IF(IPRSCF.GE.5) THEN
          DO I = 1,NFSYM
             WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
             CALL HEADER('DHFCNV: MO Fock matrix after DAMP:',-1)
             CALL PRQMAT(FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                   NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
          END DO
        ENDIF

!  save Damped FMO as New matrix on LUFOCK (PBS+MI/Aug2014)
        numfmat = max0 (5,mxdiis)
        IREC = MOD(NITER-1,numfmat-1) + 1
        CALL WRTDAC(LUFOCK,N2TMOTQ,FMO,IREC)
c
      ENDIF
 9999 CONTINUE
      CALL FLSHFO(LUPRI)
      CALL QEXIT('DHFCNV')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diiser */
      SUBROUTINE DIISER(EVEC,FMO,FAO,DMAT,TINV,EBUF,TMAT,WORK,LWORK)
C****************************************************************************
C
C  PURPOSE: Calculate the gradient (error) vector.
C
C  On input: DMAT, FAO, TMAT, TINV
C
C  On output: EVCVAL calculated (dcbdhf.h)
C
C   Last revision:   971029-jth
C                    Nov 2005 - MI
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2 = 2.0D0, DM1 = -1.00D00,
     &          D4 = 4.00D00, DP5 = 0.50D00)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbdhf.h"
C
C Used from COMMON blocks:
C    DCBBAS: N2BBASX,NTBAS(0)
C    DCBORB: NFORB,
C    DGROUP: NZ,NFSYM
C
C     Note that EBUF has dimension NFMAT*N2BBASXQ
C
      DIMENSION FAO(*),DMAT(*),EVEC(*),EBUF(*),TMAT(*),
     &          TINV(*),WORK(LWORK), FMO(*)
      CHARACTER*1 CTMP
C
#include "memint.h"
C
C
      CALL QENTER('DIISER')
C
C     DIIS in AO-basis
C     ================
C
      IF(DIISAO) THEN
C
C        The gradient is
C           g = DCao * FDao + f_o DVao * (FDao + (a_o-1) QVao) + H.C.
C             = (DCao+f_o DVao) FDao + f_o (a_o-1) DVao QVao + H.C.
C
C
         IF (AOC) THEN
C
C           Construct total density matrix in "DCao"
C
            DO IOPEN = 1,NOPEN
               CALL DAXPY(N2BBASXQ,DF(IOPEN),DMAT(1+N2BBASXQ*IOPEN),1,
     &                                       DMAT,1)
            END DO
         END IF
C
         CALL DZERO(EBUF,NFMAT*N2BBASXQ)
         CALL DZERO(EVEC,N2TMOTQ)
C
         DO 10 I = 1,NFSYM
            IF (NOCC(I).EQ.0) GOTO 10
            J2BASXI = I2BASX(I,I) + 1
C
C              (DCao+f_o DVao) FDao
C
            CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),D1,
     &                 'N','N',IPQTOQ(1,0),DMAT(J2BASXI),
     &                         NTBAS(0),NTBAS(0),NZ,
     &                 'N','N',IPQTOQ(1,0),FAO(J2BASXI),
     &                         NTBAS(0),NTBAS(0),NZ,
     &                 D0,IPQTOQ(1,0),EBUF(J2BASXI),
     &                         NTBAS(0),NTBAS(0),NZ)
            IF ( IPRSCF .GE. 10 ) THEN

             CALL HEADER(
     &       'DIISER: DMAT,total density matrix DTOTao ',-1)
             WRITE(LUPRI,'(4X,A,I1,A,I1)')
     &       'fermion symmetry: ', I,'/',NFSYM
             CALL PRQMAT(DMAT(J2BASXI),NFBAS(I,0),NFBAS(I,0),
     &                 NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)

             CALL HEADER('DIISER: Entering FDao matrix ',-1)
             WRITE(LUPRI,'(2X,A,I1,A,I1)')
     &       'fermion symmetry: ', I,'/',NFSYM
             CALL PRQMAT(FAO(J2BASXI),NFBAS(I,0),NFBAS(I,0),
     &                 NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)

CMI          ... the result of multiplication ...
             WRITE(CTMP,'(I1)') I
             CALL HEADER('DIISER: DTOTao . FDao,  corep '//CTMP,-1)
             CALL PRQMAT(EBUF(J2BASXI),NFBAS(I,0),NFBAS(I,0),
     &                   NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
            END IF
            IF (AOC) THEN
               DO 15 IOPEN = 1,NOPEN
                  IF (NACSH(I,IOPEN).EQ.0) GOTO 15
C
C                   (a_o-1) f_o DVao QVao
C                 = (a_o-1) f_o DVao f FVao
C                 where FVao is the Fock matrix from TWOINT
C
                  DSCALE = DF(IOPEN)*(DA(IOPEN)-D1)*DF(IOPEN)
                  CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),DSCALE,
     &                 'N','N',IPQTOQ(1,0),
     &                         DMAT(J2BASXI+N2BBASXQ*IOPEN),
     &                         NTBAS(0),NTBAS(0),NZ,
     &                 'N','N',IPQTOQ(1,0),
     &                         FAO(J2BASXI+N2BBASXQ*IOPEN),
     &                         NTBAS(0),NTBAS(0),NZ,
     &                D1,IPQTOQ(1,0),EBUF(J2BASXI),NTBAS(0),NTBAS(0),NZ)
                  IF ( IPRSCF .GE. 10 ) THEN
                     WRITE(CTMP,'(I1)') I
                     CALL HEADER('DIISER: (DCao+f_o DVao) FDao + '//
     &                           'f_o (a_o-1) DVao QVao, corep '//
     &                           CTMP,-1)
                     CALL PRQMAT(EBUF(J2BASXI), NFBAS(I,0),NFBAS(I,0),
     &                           NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
                  END IF
 15            CONTINUE
            END IF
C
C           Transform to MO-basis.: Ct Sao A C
C           use FAO as buffer
C
            IF(SUB_BL) THEN
              CALL TSUBBL(I,EBUF,EVEC,TINV,TMAT,
     &                      WORK(KFREE),LFREE,IPRSCF)
            ELSE
              CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                NTMO(I),NTMO(I),
     &                EBUF(J2BASXI),NTBAS(0),NTBAS(0),
     &                            NZ,IPQTOQ(1,0),
     &                EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &                TINV(I2TMT(I)+1),NFBAS(I,0),NTMO(I),
     &                            NZT,IPQTOQ(1,0),
     &                TMAT(I2TMT(I)+1),NFBAS(I,0),NTMO(I),
     &                            NZT,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,IPRSCF)

                IF (IPRSCF.GE.10) THEN
                 CALL HEADER('DIISER: Entering TINV',-1)
                 WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &           '*** Fermion ircop ',I,'/',NFSYM
                 CALL PRQMAT(TINV(I2TMT(I)+1),NFBAS(I,0),
     &           NTMO(I),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),LUPRI)

                 CALL HEADER('DIISER: Entering TMAT',-1)
                 WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &           '*** Fermion ircop ',I,'/',NFSYM
                 CALL PRQMAT(TMAT(I2TMT(I)+1),NFBAS(I,0),
     &           NTMO(I),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),LUPRI)
                ENDIF

            ENDIF
C
            IF (ZORA) CALL ZRDIIS (I,EVEC(I2TMOT(I)+1),WORK,KFREE,LFREE)
C
            IF (IPRSCF.GE.5) THEN
               CALL HEADER('DIISER DIISMO: FD*D matrix',-1)
               WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
               CALL PRQMAT(EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                     NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
C
C           Calculate ErrVec - ErrVec(Hermitian conj.)
C
            CALL DIFMMH(EVEC(I2TMOT(I)+1),NTMO(I),NZ,NTMO(I),NTMO(I))
C
            IF ( IPRSCF .GE. 10 ) THEN
               WRITE(CTMP,'(I1)') I
               CALL HEADER('DIISER: Error vector in MO-basis, corep '
     &                     //CTMP,-1)
               CALL PRQMAT(EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                  NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
C
 10      CONTINUE
C
      ELSE
C
C     DIIS in MO-basis
C     ================
C
C        Construct density matrix in MO-basis
C        ====================================
C        Read coefficients in orthonormal basis
         REWIND LUCMOS
         CALL READT(LUCMOS,N2TMOTQ,EVEC)
C        Observe that we use NISHMF which is the number of inactive orbitals
C        minus those that are frozen
         DO 20 I = 1,NFSYM
            IF(NOCCMF(I).EQ.0) GOTO 20
            IF(IPRSCF.GE.3) THEN
C             Transform Fock matrix to current orthonormal MO-basis
              CALL QTRANS('AOMO','S',D0,NTMO(I),NTMO(I),
     &                NTMO(I),NTMO(I),
     &                FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &                EBUF(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &                EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &                EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,IPRSCF)
              KOFF = I2TMOT(I) + NPSHMF(I)*NTMO(I) + 1
              CALL HEADER('DIISER: '//
     &         'Fock matrix (1:nesh,1:nocc) in current MO-basis',-1)
              WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
              CALL PRQMAT(EBUF(KOFF),NESHMF(I),NOCCMF(I),
     &                    NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
            ENDIF
C...........closed shell
            ISTART=NPSHMF(I)+1
            CALL DENST1(EBUF(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,
     &                  D1,D0,
     &                  EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                  ISTART,NISHMF(I),NTMO(I))
C...........open shells
            ISTART = ISTART+NISHMF(I)
            DO IOPEN = 1,NOPEN
              CALL DENST1(EBUF(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,
     &                    DF(IOPEN),D1,
     &                    EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                    ISTART,NACSHMF(I,IOPEN),NTMO(I))
              ISTART = ISTART + NACSHMF(I,IOPEN)
            ENDDO
            IF(IPRSCF.GE.5) THEN
               CALL HEADER('DIISER: Density matrix in MO-basis',-1)
               WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &          '*** Fermion ircop ',I,'/',NFSYM
               CALL PRQMAT(EBUF(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                   NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
            ENDIF
 20      CONTINUE
C
C        Form error vector
C        =================
C
         CALL DZERO(EVEC,N2TMOTQ)
         DO 30 I = 1,NFSYM
           IF(NOCCMF(I).EQ.0) GOTO 30
           IF(SUB_BL) THEN
             IBO = 0
             DO ISUB = 1,N_SUB_BL(I)
             NBO = NTMO_SUB(ISUB,I,0)
             IF(NBO.GT.0) THEN
               IMAT = I2TMOT(I) + (NTMO(I)+1)*IBO + 1
               CALL QGEMM(NBO,NBO,NBO,D1,
     &         'N','N',IPQTOQ(1,0),FMO(IMAT),NTMO(I),NTMO(I),1,
     &         'N','N',IPQTOQ(1,0),EBUF(IMAT),NTMO(I),NTMO(I),1,
     &         D0,IPQTOQ(1,0),EVEC(IMAT),NTMO(I),NTMO(I),1)
               IBO = IBO + NTMO_SUB(ISUB,I,0)
             ENDIF
             ENDDO
           ELSE
             CALL QGEMM(NTMO(I),NTMO(I),NTMO(I),D1,
     &        'N','N',IPQTOQ(1,0), FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,
     &        'N','N',IPQTOQ(1,0),EBUF(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,
     &        D0,IPQTOQ(1,0),EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ)
           ENDIF
           IF(IPRSCF.GE.5) THEN
             CALL HEADER('DIISER DIISMO: FD*D matrix',-1)
             WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
             CALL PRQMAT(EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                   NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
           ENDIF
C
           IF (ZORA) CALL ZRDIIS (I,EVEC(I2TMOT(I)+1),WORK,KFREE,LFREE)
C
           CALL DIFMMH(EVEC(I2TMOT(I)+1),NTMO(I),NZ,NTMO(I),NTMO(I))
           IF(IPRSCF.GE.5) THEN
             CALL HEADER('DIISER DIISMO: Error matrix',-1)
             WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
             CALL PRQMAT(EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                   NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
           ENDIF
   30    CONTINUE
C
      END IF
      IF(IPRSCF.GE.2) THEN
C     Transform error vector to current orthonormal MO-basis
        REWIND LUCMOS
        CALL READT(LUCMOS,N2TMOTQ,EBUF)
        DO 40 I = 1,NFSYM
          IF(NOCCMF(I).EQ.0) GOTO 40
          NDIM = NTMO(I)*NTMO(I)*NZ
          CALL MEMGET('REAL',KBUF,NDIM,WORK,KFREE,LFREE)
          CALL QTRANS('AOMO','S',D0,NTMO(I),NTMO(I),
     &                NTMO(I),NTMO(I),
     &                EVEC(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &                WORK(KBUF),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &                EBUF(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &                EBUF(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,IPRSCF)
          CALL HEADER
     &       ('DIISER DIISMO: Error matrix in current MO-basis',-1)
          WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
          KOFF = KBUF + NPSHMF(I)*NTMO(I)
          IF(IPRSCF.GE.3) THEN
            CALL PRQMAT(WORK(KOFF),NTMO(I),NOCCMF(I),
     &                  NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
          NDIM = NTMO(I)*NTMO(I)
          DO J1 = 1,NOCCMF(I)
            GMAX = D0
            J2 = 0
            KOFF = KOFF + NPSH(I)
            DO K = 1,NESHMF(I)
              GTMP=WORK(KOFF)*WORK(KOFF)
              DO IZ = 2,NZ
                KK = KOFF + NDIM
                GTMP = GTMP + WORK(KK)*WORK(KK)
              ENDDO
              IF(GTMP.GT.GMAX) THEN
                J2 = K
                GMAX = GTMP
              ENDIF
              KOFF = KOFF + 1
            ENDDO
            WRITE(LUPRI,'(A,I5,A,I5,E12.4)')
     &       'Occ. orb. ',J1, ' largest gradient partner: ',J2,
     &        SQRT(GMAX)
          ENDDO
          CALL MEMREL('DIISER.MOerr',WORK,KWORK,KBUF,KFREE,LFREE)
 40     CONTINUE
      ENDIF
C
C     Find EVCVAL = total norm of gradient in EVEC
C     ============================================
C
C        The factor:
C
C        Missing 2 in density matrices    : 2
C        from g,g* norm to g norm         : 1/sqrt(2)
C        --------------------------------------------
C        total factor to get g norm         sqrt(2)
C
      EVCVAL = SQRT(D2)*DNRM2(N2TMOTQ,EVEC,1)
C
      IF (IPRSCF.GE.5) THEN
        CALL HEADER('DIISER: Resulting value of EVCVAL ',-1)
        WRITE(LUPRI,'(4X,A,D15.6)')
     &  '>>> EVCVAL = SQRT(D2)*DNRM2(N2TMOTQ,EVEC,1) = ',EVCVAL
      ENDIF
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('DIISER')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fckcmp */
      SUBROUTINE FCKCMP(FNEW,FOLD,WORK,LWORK)
C*****************************************************************************
C
C     Find largest change in Fock matrix between two iterations
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from COMMON blocks:
C   DCBORB: N2ORBTQ
C   DCIDHF: NFMAT
C   DGROUP: NZ
C
#include "dcbdhf.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"
      DIMENSION FNEW(*),FOLD(*),WORK(LWORK)
C
      CALL QENTER('FCKCMP')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('REAL',KDIFF,N2TMOTQ,WORK,KFREE,LFREE)
C
      FCKVAL =  DIFFMX(FNEW,FOLD,WORK(KDIFF),N2TMOTQ)
C
C     Memory deallocation
      CALL MEMREL('FCKCMP.DIFFMX',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT ('FCKCMP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diffmx */
      FUNCTION DIFFMX(AVEC,BVEC,DIFF,NELM)
C*****************************************************************************
C
C     Find maximum differnce between two vectors
C
C     Written by T.Saue Aug 25 1995
C     Last revision Aug 25 1995 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION AVEC(NELM),BVEC(NELM),DIFF(NELM)
C
      DO I = 1,NELM
        DIFF(I) = AVEC(I) - BVEC(I)
      ENDDO
      IMAX = IDAMAX(NELM,DIFF,1)
      DIFFMX = DIFF(IMAX)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diisdr */
      SUBROUTINE DIISDR(BMAT,EVEC,FOCK,EBUF,IBSTAT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for DIIS - Direct Inversion of Iterative Subspaces
C
C     Written by Jon Laerdahl and T. Saue 1995
C     Last revision : Aug 28 1995 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from COMMON blocks:
C   dcbdhf: MXDIIS
C
#include "dcbdhf.h"
      DIMENSION BMAT(*),EVEC(*),FOCK(*),EBUF(*),WORK(LWORK)
C
      CALL QENTER('DIISDR')
#include "memint.h"
      CALL MEMGET('REAL',KCVEC,MXDIIS,WORK,KFREE,LFREE)
      CALL DIISD1(BMAT,EVEC,FOCK,EBUF,WORK(KCVEC),IBSTAT,
     &            WORK,KFREE,LFREE)
      CALL MEMREL('DIISDR.DIISD1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('DIISDR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diisd1 */
      SUBROUTINE DIISD1(BMAT,EVEC,FOCK,EBUF,CVEC,IBSTAT,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Driver for DIIS - Direct Inversion of Iterative Subspaces
C
C     Written by T.Saue Aug 28 1995
C     Last revision Aug 28 1995
C                   Miro ILIAS, febr.2007, Tel Aviv - fix for 2c->4c restart
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0,D0 = 0.0D0)
C
C Used from COMMON blocks:
C   dcbdhf: LUEVEC,IPRSCF,ITDIIS,MXDIIS
C
#include "dcbham.h"
#include "dcbdhf.h"
      DIMENSION BMAT(MXDIIS*MXDIIS),EVEC(*),FOCK(*),EBUF(*),CVEC(*),
     &          WORK(*)
      logical exists

      IPRINT = IPRSCF
C
C     Write current error vector to direct access unit
C     ================================================
C
      IREC = MOD(ITDIIS-1,MXDIIS-1) + 1
      CALL WRTDAC(LUEVEC,N2TMOTQ,EVEC,IREC)
C
C     First DIIS iteration; initialize
C     ================================
C
      IF(ITDIIS.EQ.1) THEN
C       BMAT is triangular packed.
C       Set first column in BMAT; used for constraint:
        NELMBM = 1
        BMAT(1) = D0
        II = 1
        DO I = 1,MXDIIS-1
          II = II + I
          BMAT(II) = DM1
        ENDDO
      ENDIF
      IF(ITDIIS.EQ.2) CALL OPNFIL(LUDIIS,'DFDIIS','UNKNOWN','DIISDR')
C
C     Solve DIIS equation for current iteration
C     ========================================
C
      CALL EQDIIS(BMAT,CVEC,EVEC,EBUF,IBSTAT,IPRINT,
     &            WORK(KFREE),LFREE)
C
C     Construct total Fock matrix using coefficients of CVEC
C
      IF(IBSTAT.EQ.0) THEN
        CALL FMDIIS(FOCK,EBUF,CVEC,IPRINT)
        REWIND LUDIIS
        WRITE (LUDIIS) ITDIIS,MXDIIS,NELMBM,BMAT
      ELSEIF(IBSTAT.EQ.1) THEN
         IF (.NOT.(START2C.AND..NOT.DO4C2C)) THEN
           inquire(LUDIIS,exist=exists)
           if(DO4C2C)then
!            stefan: if we delete the file here, the program crashes
!            (because the file does not exist on disc) when
!            we restart the 2c-SCF after we have done a 4c-SCF +
!            transformation of the 4c-Fock operator.
!            @ miro: did this ever work after your change 2007?
!            i tried it also with older trunk versions and it did not...
             if(exists) CLOSE(LUDIIS)
           else
             if(exists) CLOSE(LUDIIS,STATUS='DELETE')
           end if
         ENDIF
      ENDIF
C
      IF(IPRINT.GE.1) THEN
        CALL HEADER('Output from DIISD1:',-1)
        WRITE(LUPRI,'(A,I10)') '* DIIS iteration no.:',ITDIIS
        WRITE(LUPRI,'(A,I10)') '* B matrix dimension:',NELMBM
        WRITE(LUPRI,'(A,1P,D10.3)') '* B matrix condition:',BMCOND
        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,'(A)') '* Current B matrix:'
          CALL PRMUTC(BMAT,NELMBM,1,LUPRI)
        ENDIF
        IF(IBSTAT.EQ.0) THEN
          CFSUM = D0
          DO I = 2,NELMBM
            CFSUM = CFSUM + CVEC(I)
          ENDDO
          WRITE(LUPRI,'(/A,1P,D9.3/)') '* DIIS coefficients. Sum:',CFSUM
          IOFF = NITER - NELMBM
          WRITE(LUPRI,'(5X,A6,4X,1P,D10.3)')
     &       'lambda',CVEC(1)
          WRITE(LUPRI,'(5X,I5,5X,1P,D10.3)')
     &       (IOFF+I,CVEC(I),I=2,NELMBM)
        ENDIF
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck eqdiis */
      SUBROUTINE EQDIIS(BMAT,CVEC,EVEC,EBUF,IBSTAT,IPRINT,
     &                  WORK,LWORK)
C*****************************************************************************
C
C     Solve DIIS equations
C
C     Written by Jon Laerdahl & Trond Saue
C     Last revision : Aug 28 1995 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C Used from COMMON blocks:
C   dcbdhf: MXDIIS
#include "dcbdhf.h"
      DIMENSION BMAT(*),CVEC(*),EVEC(*),EBUF(*),WORK(LWORK)
C
      CALL QENTER('EQDIIS')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('REAL',KBFAC,(MXDIIS*(MXDIIS+1))/2,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KKPVT,MXDIIS,WORK,KFREE,LFREE)
      CALL EQDII1(BMAT,WORK(KBFAC),EVEC,EBUF,CVEC,IBSTAT,
     &            WORK(KKPVT),IPRINT,WORK,KFREE,LFREE)
C
C     Memory deallocation
      CALL MEMREL('EQDIIS.EQDII1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('EQDIIS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck eqdii1 */
      SUBROUTINE EQDII1(BMAT,BFAC,EVEC,EBUF,CVEC,IBSTAT,KPVT,
     &                  IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Solve DIIS equations
C
C     Written by Jon Laerdahl & Trond Saue
C     Last revision : Aug 28 1995 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0, D0 = 0.00D00)
      DIMENSION BMAT(*),BFAC(*),EVEC(*),EBUF(*),CVEC(*),
     &          KPVT(*),WORK(*)
C Used from COMMON blocks
C  dcbdhf: NELMBM
C
#include "dcbdhf.h"
C
C     Construct and factorize B - matrix
C     ==================================
C
      CALL BMDIIS(BMAT,BFAC,EVEC,EBUF,CVEC,KPVT,
     &            IBSTAT,IPRINT)
C
C     Solve DIIS equation
C     ===================
C
      IF(IBSTAT.EQ.0) THEN
         CALL DZERO(CVEC,NELMBM)
         CVEC(1) = DM1
         CALL DSPSL(BFAC,NELMBM,KPVT,CVEC)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bmdiis */
      SUBROUTINE BMDIIS(BMAT,BFAC,EVEC,EBUF,CVEC,KPVT,IBSTAT,
     &                  IPRINT)
C*****************************************************************************
C
C     Calculate and factor DIIS B matrix.
C     If the matrix is ill-conditioned
C     Elements are of form B_ij = <e_i|e_j> where e_i,e_j are error vectors.
C     The B matrix is symmetric and is therefore upper triangular
C     columnwise packed.
C
C     Written by Jon Laerdahl & T.Saue
C     Last revision: Aug 27 1995 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
C Used from COMMON blocks:
C  dcbdhf: LUEVEC,NELMBM,ITDIIS,MXDIIS
C  DGROUP: NZ
C  DCBBAS: N2VCTQ
C
#include "dcbdhf.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION BMAT(*),BFAC(*),EVEC(*),EBUF(*),CVEC(*),KPVT(*)
      CALL QENTER('BMDIIS')
C
      NELMBM  = NELMBM + 1
C
C     If necessary, make room for more elements in B matrix
C
      IF(NELMBM.GT.MXDIIS) THEN
        NELMBM = NELMBM - 1
        CALL BSHIFT(BMAT,NELMBM)
      ENDIF
C
C     Add new column to B matrix; IOFF is offset to
C     first element of new column
C
      IOFF  = (NELMBM*(NELMBM-1))/2 + 1
      NPREV = NELMBM - 2
      JREC  = ITDIIS - NPREV - 2
      DO I = 1, NPREV
        IREC = MOD(JREC+I,MXDIIS-1) + 1
        CALL READAC(LUEVEC,N2TMOTQ,EBUF,IREC)
        BMAT(IOFF+I) = BIJELM(EBUF,EVEC)
      ENDDO
      BMAT(IOFF+NPREV+1) = BIJELM(EVEC,EVEC)
C
C     Factor and check condition of B matrix
C     ======================================
C
   10 CONTINUE
      IF(NELMBM.LE.2) THEN
        IBSTAT = 1
        GOTO 999
      ELSE
        NXBM = (NELMBM*(NELMBM+1))/2
        CALL DCOPY(NXBM,BMAT,1,BFAC,1)
        CALL DSPCO(BFAC,NELMBM,KPVT,BMCOND,CVEC)
        IF(IPRINT.GE.5) THEN
          WRITE(LUPRI,'(A,I5)') 'ITDIIS: ',ITDIIS,
     & 'NELMBM: ',NELMBM
          WRITE(LUPRI,'(A,1P,E9.3)') 'BMCOND : ',BMCOND
        ENDIF
C
C       Ill-conditioned matrix; delete oldest elements
C
        IF (D1+BMCOND.EQ.D1) THEN
          IF(IPRINT.GE.2) THEN
            WRITE(LUPRI,'(A,1P,E9.3)') '***  INFO  *** B-cond:',BMCOND
            WRITE(LUPRI,'(/A/)') 'B matrix:'
            CALL PRMUTC(BMAT,NELMBM,1,LUPRI)
          ENDIF
          CALL BSHIFT(BMAT,NELMBM)
          NELMBM = NELMBM - 1
          GOTO 10
        ENDIF
        IBSTAT = 0
      ENDIF
C
 999  CONTINUE
      CALL QEXIT ('BMDIIS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bijelm */
      FUNCTION BIJELM(EVEC1,EVEC2)
C*****************************************************************************
C
C     Calculate an element of the DIIS B matrix
C
C     Written by T.Saue Aug 28 1995
C     Last revision: Oct 18 1995 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D2=2.0D0)
C
C Used from COMMON blocks:
C   dcbdhf.h: N2TMOTQ
C   dcborb.h: NASHT, NAELEC
C
#include "dcbdhf.h"
#include "dcborb.h"

      DIMENSION EVEC1(*),EVEC2(*)
      BIJELM = D2*DDOT(N2TMOTQ,EVEC1,1,EVEC2,1)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bshift */
      SUBROUTINE BSHIFT(BMAT,NELMBM)
C*****************************************************************************
C
C     Delete row 2 and column 2 of B matrix and shift all the
C     other elements accordingly diagonally upwards, that is
C       B(i,j) --> B(i-1,j-1)
C     The B matrix is upper triangularly columnwise packed, so
C     the transformation is
C       B([ij]) --> B([ij] - j)
C     The formula for a general shift
C       B(i,j) --> B(i-n,j-n)
C     is
C       B([ij]) --> B([ij]-N), N = [nj - n(n-1)/2]
C
C     Written by T.Saue Aug 28 1995
C     Last revision Aug 28 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION BMAT(*)
C
C     J runs over columns, I over rows;
C     ID gives the supervectorindex of diagonal element of previous
C     column
C
      ID = 3
      DO J = 3,NELMBM
        DO I = 3,J
          IJ = ID + I
          BMAT(IJ-J) = BMAT(IJ)
        ENDDO
        ID = ID + J
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fmdiis */
      SUBROUTINE FMDIIS(FOCK,FBUF,CVEC,IPRINT)
C*****************************************************************************
C
C     Construct total Fock matrix using coefficients obtained
C     by solving DIIS equations
C
C     Written by T.Saue Aug 28 1995
C     Last revision : Aug 28 1995 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from COMMON blocks:
C  dcbdhf: LUFOCK,NITER,NELMBM
C  DCIDHF: NFMAT
C  DGROUP: NZ
C
#include "dcbdhf.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"
      DIMENSION FOCK(*),FBUF(*),CVEC(*)
C
      CALL DSCAL(N2TMOTQ,CVEC(NELMBM),FOCK,1)
      IOFF  = 1
      NPREV = NELMBM - 2
      JREC  = NITER - NPREV - 2
      DO I = 1,NPREV
        IREC = MOD(JREC+I,MXDIIS-1) + 1
        CALL READAC(LUFOCK,N2TMOTQ,FBUF,IREC)
        CALL DAXPY(N2TMOTQ,CVEC(IOFF+I),FBUF,1,FOCK,1)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck zrdiis */
      SUBROUTINE ZRDIIS(I,EVEC,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Modify DIIS matrix to account for ZORA
C     When calculating the error vector we need to account for the
C     zero SS-metric
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dgroup.h"
C
      PARAMETER(D0=0.0D0)
      DIMENSION EVEC(*),WORK(*)
C
      NBRP = 4/NZ
      NL = NESH(I)
      NS = NPSH(I)
C
C     We get a chunk of memory, fill it with zeroes and insert it in
C     the full matrix so that the SS block becomes zero.
C
      CALL MEMGET('REAL',KESS,NZ*NS*NS,WORK,KFREE,LFREE)
      CALL DZERO (WORK(KESS),NZ*NS*NS)
      CALL EXTRSB (2,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &             0,0,2,2,1,NZ,
     &             WORK(KESS),NS,NS,NZ,
     &             EVEC,NTMO(I),NTMO(I),NZ)
      CALL MEMREL('ZRDIIS',WORK,KESS,KESS,KFREE,LFREE)
C
C     Same procedure for LS block
C
      CALL MEMGET('REAL',KELS,NZ*NL*NS,WORK,KFREE,LFREE)
      CALL DZERO (WORK(KELS),NZ*NL*NS)
      CALL EXTRSB (2,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &             0,0,1,2,1,NZ,
     &             WORK(KELS),NS,NS,NZ,
     &             EVEC,NTMO(I),NTMO(I),NZ)
      CALL MEMREL('ZRDIIS',WORK,KELS,KELS,KFREE,LFREE)
C
      RETURN
      END
! -- end of dirac/dircnv.F --
