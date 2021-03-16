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

!
!  -- File : huckel.F --
!
      SUBROUTINE HUCDRV(FMAT,NBAS_LS,NELECT,NOCCUP,WORK,LWORK,IPRINT)
C
C     Driver routine for solving the extended Huckel problem
C     panor+hjj 2010
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "huckel.h"
C
      LOGICAL   AUTO_OCC
      DIMENSION FMAT(*),NOCCUP(8),WORK(LWORK)
C
      CALL QENTER('HUCDRV')
#ifndef MOD_HUCKEL
      call quit('HUCKEL not active in this version')
#else
C
C     This routine does memory allocation for HUCDRV_1
C
      NSYM = MAXREP + 1          ! number boson irreps
      NBASHUC = NCOS(0,0)
      NHUCAO(1:8) = NCOS(1:8,0)
      IPRHUC = MAX(IPRINT,IPRHUC)
C
      IF (NBAS_LS.NE.NCOS(0,1)+NCOS(0,2) .OR. NBAS_LS.NE.NTBAS(0)) THEN
         WRITE(lupri,*) 'ERROR: NBAS_LS,NBAS_L,NBAS_S,NBAS_L+NBAS_S=',
     &      NBAS_LS,NCOS(0,1),NCOS(0,2),NCOS(0,1)+NCOS(0,2),NTBAS(0)
         CALL QUIT('Inconsistency in number of basis functions')
      END IF
C
      AUTO_OCC=.TRUE.
      DO I=1,NSYM
         AUTO_OCC = AUTO_OCC .AND. NOCCUP(I).LE.0
      END DO
C
      IF (IPRHUC.GE.0) THEN
         WRITE (LUPRI,'(/A)') ' Huckel start guess based on:'
         IF (EWMO) THEN
            WRITE(LUPRI,'(2A)')
     &         '   EWMO - Energy Weighted Maximum Overlap method,',
     &         ' see Linderberg and Ohrn (Wiley, 1973)'
         ELSE
            WRITE(LUPRI,'(2A/A,F10.3)')
     &         '   EHT - Extended Huckel method, see R. Hoffman',
     &         ' (J. Chem. Phys. 39, 1397, 1963)',
     &         '   Factor for resonance terms :',HUCCNT
         END IF
         WRITE(LUPRI,'(A,I5)') ' Number electrons       :', NELECT
         WRITE(LUPRI,'(A,I5)') ' Number boson irreps    :', NSYM
         WRITE(LUPRI,'(A,I5)') ' Number basis functions :', NBAS_LS
         WRITE(LUPRI,'(A,I5)') ' Number Huckel functions:', NBASHUC
         WRITE(LUPRI,'(A,8I5)')
     &        ' Number Huckel functions in symmetries:',
     &        (NHUCAO(I),I=1,NSYM)
         IF (AUTO_OCC) THEN
            WRITE(LUPRI,'(A)')
     &           ' Auto occupation requested in Huckel routine'
         ELSE
            WRITE(LUPRI,'(A,8I5)')
     &           ' Occupation specified:',
     &        (NOCCUP(I),I=1,NSYM)
         END IF
      END IF
C
      CALL FLSHFO(LUPRI)
C
      KFREE = 1
      LFREE = LWORK
      N2BASHUC = NBASHUC*NBASHUC
      N2BAS_LS  = NBAS_LS*NBAS_LS
      CALL MEMGET('REAL',KHUCMAT,N2BASHUC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHUCOVL,N2BASHUC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVEC,N2BASHUC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVAL,NBASHUC, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KF2MAT, N2BAS_LS,WORK,KFREE,LFREE)
      CALL HUCDRV_1(FMAT,WORK(KF2MAT),NELECT,NSYM,NOCCUP,AUTO_OCC,
     &            NBAS_LS,NBASHUC,WORK(KHUCMAT),WORK(KHUCOVL),
     &            WORK(KEIGVEC),WORK(KEIGVAL),WORK(KFREE),LFREE)
C
      CALL MEMREL('HUCDRV',WORK,1,1,KFREE,LFREE)
#endif 
      CALL QEXIT('HUCDRV')
      RETURN
      END
C ========================================================================
      SUBROUTINE HUCDRV_1(FMAT,F2MAT,NELECT,NSYM,NOCCUP,AUTO_OCC,
     &     NBAS_LS,NBASHUC,HUCMAT,HUCOVL,EIGVEC,EIGVAL,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
C
      LOGICAL   AUTO_OCC
      DIMENSION F2MAT(*),NOCCUP(8),WORK(LWORK)
      DIMENSION HUCMAT(*),HUCOVL(*),EIGVEC(*),EIGVAL(*)
C
C     We obtain a Huckel start guess in 3 steps:
C
C     1. Construct the Huckel matrix H and diagonalize it according to
C                          H*C = S*E*C
C     2. Given the eigenvectors C we obtain a density matrix D
C     3. From D we obtain the 2-electron Fock matrix in the basis set defined in
C        the user input
C
      KFREE = 1
      LFREE = LWORK
C
C     Generate Huckel Hamiltonian and overlap matrices by use of subroutine
C     PR1INT_1
C
      CALL HUCGETMAT(NSYM,NBASHUC,HUCMAT,HUCOVL,WORK(KFREE),LFREE)
C
C     Solve the Huckel eigenvalue equation and, if not user specified,
C     determine orbital occupation numbers to be stored in
C     NOCCUP(1:NSYM)
C
      CALL HUCEIG(NELECT,NSYM,NOCCUP,AUTO_OCC,NBASHUC,HUCMAT,HUCOVL,
     &     EIGVEC,EIGVAL,WORK(KFREE),LFREE)
C
C     Construct the Huckel density matrix (reuse HUCMAT for storage of
C     density matrix)
C
      CALL HUCDEN(HUCMAT,NBASHUC,NSYM,NOCCUP,EIGVEC)
C
C     Determine the Fock matrix in the user-defined basis set
C
      CALL HUCFCK(FMAT,F2MAT,HUCMAT,NBAS_LS,NBASHUC,WORK(KFREE),LFREE)
C
      RETURN
      END
C ========================================================================
      SUBROUTINE HUCGETMAT(NSYM,NBASHUC,HUCMAT,HUCOVL,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "huckel.h"
C
      DIMENSION HUCMAT(NBASHUC,NBASHUC),HUCOVL(NBASHUC,NBASHUC),
     &          WORK(LWORK)

      CHARACTER LABINT(2)*8, WORD*7
      INTEGER   INTREP(2)
      LOGICAL   TOFILE, DOINT(-MXBSETS:MXBSETS,2), PROPRI
      PARAMETER (TOFILE = .FALSE.)

      DOINT(:,:) = .FALSE.
      DOINT(0,1) = .TRUE.
      DOINT(0,2) = .TRUE.

      KFREE = 1
      LFREE = LWORK
      N2BASHUC = NBASHUC*NBASHUC
C
C     Memory allocation
C
      CALL MEMGET('REAL',KHUCMATS,2*N2BASHUC,WORK,KFREE,LFREE)
C
C     Generate Huckel and overlap matrices
C
      WORD   = 'HUCKEL '
      PROPRI = IPRHUC .GE. 5
      IPRINT = IPRHUC - 6
      NCOMP  = -2
      CALL PR1INT_1(WORK,KFREE,LFREE,INTREP,DUMMY,LABINT,
     &            WORD,0,0,.TRUE.,PROPRI,IPRINT,
     &            WORK(KHUCMATS),NCOMP,TOFILE,'SQUARE',
     &            DOINT,EXPVAL,.FALSE.,DUMMY)

      CALL DCOPY(N2BASHUC,WORK(KHUCMATS),1,HUCOVL,1)

      IF (EWMO) THEN
C     EWMO : generate H = SWS, requires C1 boson symmetry
         KW = KHUCMATS + N2BASHUC
         CALL DGEMM('N','N',NBASHUC,NBASHUC,NBASHUC,1.0D0,
     &              WORK(KW),NBASHUC, HUCOVL,NBASHUC,0.0D0,
     &              WORK(KHUCMATS),NBASHUC)
         CALL DGEMM('N','N',NBASHUC,NBASHUC,NBASHUC,1.0D0,
     &              HUCOVL,NBASHUC, WORK(KHUCMATS),NBASHUC,0.0D0,
     &              HUCMAT,NBASHUC)

      ELSE ! EHT

         CALL DCOPY(N2BASHUC,WORK(KHUCMATS+N2BASHUC),1,HUCMAT,1)

      END IF
C
C     Print
C
      IF (IPRHUC.GE.3) THEN
         IF (EWMO) THEN
            WRITE (LUPRI,'(/A)') ' EWMO AO Hamiltonian matrix W'//
     &           ' corresponding to separate atoms:'
            CALL OUTPUT(WORK(KW),1,NBASHUC,1,NBASHUC,NBASHUC,NBASHUC,
     &           -1,LUPRI)
         END IF
         WRITE (LUPRI,'(/A)') ' Huckel Hamiltonian matrix H:'
         CALL OUTPUT(HUCMAT,1,NBASHUC,1,NBASHUC,NBASHUC,NBASHUC,
     &        -1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Huckel overlap matrix S:'
         CALL OUTPUT(HUCOVL,1,NBASHUC,1,NBASHUC,NBASHUC,NBASHUC,
     &        -1,LUPRI)
      END IF
C
C     Memory deallocation
C
      CALL MEMREL('HUCGETMAT',WORK,1,1,KFREE,LFREE)
C
      RETURN
      END
C ========================================================================
      SUBROUTINE HUCEIG(NELECT,NSYM,NOCCUP,AUTO_OCC,
     &     NBASHUC,HUCMAT,HUCOVL,EIGVEC,EIGVAL,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "pgroup.h"
#include "huckel.h"
C
      LOGICAL   AUTO_OCC
      INTEGER   NOCCUP(8), NOCC_HUC(8)
      DIMENSION HUCMAT(NBASHUC,NBASHUC),
     &     HUCOVL(NBASHUC,NBASHUC),
     &     EIGVEC(NBASHUC,NBASHUC),EIGVAL(NBASHUC),
     &     WORK(LWORK)
C
      KFREE = 1
      LFREE = LWORK
C
      CALL MEMGET('REAL',KH,NBASHUC*NBASHUC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KS,NBASHUC*NBASHUC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KC,NBASHUC*NBASHUC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCR,2*NBASHUC,WORK,KFREE,LFREE)
C
      CALL DZERO(EIGVEC,NBASHUC*NBASHUC)
C
      IPOS=1
      DO ISYM=1,NSYM
C
         IF (NHUCAO(ISYM).GT.0) THEN
C
            CALL DCOPYRECT(NBASHUC,NBASHUC,HUCMAT,IPOS,IPOS,
     &           IPOS+NHUCAO(ISYM)-1,IPOS+NHUCAO(ISYM)-1,
     &           NHUCAO(ISYM),NHUCAO(ISYM),
     &           WORK(KH),1,1,NHUCAO(ISYM),NHUCAO(ISYM))
            CALL DCOPYRECT(NBASHUC,NBASHUC,HUCOVL,IPOS,IPOS,
     &           IPOS+NHUCAO(ISYM)-1,IPOS+NHUCAO(ISYM)-1,
     &           NHUCAO(ISYM),NHUCAO(ISYM),
     &           WORK(KS),1,1,NHUCAO(ISYM),NHUCAO(ISYM))
C
            CALL RSG(NHUCAO(ISYM),NHUCAO(ISYM),WORK(KH),WORK(KS),
     &           EIGVAL(IPOS),1,WORK(KC),
     &           WORK(KSCR),WORK(KSCR+NBASHUC),IERR)
            IF (IERR.NE.0) THEN
               CALL QUIT('HUCEIG: eigenvalue equation not converged')
            END IF
C
C           RSG returns eigenvalues in ascending order.
C
            CALL DCOPYRECT(NHUCAO(ISYM),NHUCAO(ISYM),
     &           WORK(KC),1,1,NHUCAO(ISYM),NHUCAO(ISYM),
     &           NBASHUC,NBASHUC,EIGVEC,
     &           IPOS,IPOS,IPOS+NHUCAO(ISYM)-1,IPOS+NHUCAO(ISYM)-1)
C
            IPOS = IPOS + NHUCAO(ISYM)
C
         END IF
C
      END DO
C
C     Find closed-shell orbital occupation predicted by HUCKEL
C
C     HJAaJ: TODO - use fractional occupations if E(HOMO) close to E(LUMO) !!!
C
      KEIG=KH     ! reuse temporary storage for ordered eigenvalues
      CALL DCOPY(NBASHUC,EIGVAL,1,WORK(KEIG),1)
      CALL ORDER(DUMMY,WORK(KEIG),NBASHUC,-1)
      E_HOMO=WORK(KEIG+NELECT/2-1)
      NOCC_HUC(1:8) = 0
      ICTR=0
      DO ISYM=1,NSYM
         DO I=1,NHUCAO(ISYM)
            ICTR=ICTR+1
            IF (EIGVAL(ICTR).LE.E_HOMO) THEN
               NOCC_HUC(ISYM) = NOCC_HUC(ISYM)+1
            END IF
         END DO
      END DO

      WRITE(LUPRI,'(/A,8I4)') ' Orbital occupation by Huckel: ',
     &     (NOCC_HUC(I),I=1,NSYM)
      IF (AUTO_OCC) THEN
         NOCCUP(1:8) = NOCC_HUC(1:8)
      ELSE
         WRITE(LUPRI,'(/A,8I4)') ' However, using user-specified:',
     &     (NOCCUP(I),I=1,NSYM)
      END IF
C
C     Print Huckel eigenvalues and eigenvectors
C
      IF (IPRHUC.GE.0) THEN
         I_END = 0
         DO ISYM = 1,NSYM
            I_ST  = I_END + 1
            I_END = I_END + NHUCAO(ISYM)
            WRITE (LUPRI,'(/A,I3,3A)')
     &      ' Huckel eigenvalues symmetry',ISYM,' (',REP(ISYM-1),')'
            WRITE (LUPRI,'(5X,5F15.6)') (EIGVAL(I), I=I_ST, I_END)
         END DO
      END IF
C
      IF (IPRHUC.GE.2) THEN
         WRITE (LUPRI,'(/A)') ' Huckel eigenvectors:'
         CALL OUTPUT(EIGVEC,1,NBASHUC,1,NBASHUC,
     &        NBASHUC,NBASHUC,-1,LUPRI)
      END IF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
C ========================================================================
      SUBROUTINE HUCDEN(DENHUC,NBASHUC,NSYM,NOCCUP,EIGVEC)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "huckel.h"
C
      DIMENSION DENHUC(NBASHUC,NBASHUC),EIGVEC(NBASHUC,NBASHUC),
     &     NOCCUP(8)
C
      CALL DZERO(DENHUC,NBASHUC*NBASHUC)
C
      IOFF=0
      DO ISYM=1,NSYM
         DO I=1,NHUCAO(ISYM)
         DO J=1,NHUCAO(ISYM)
            DO IOCC=1,NOCCUP(ISYM)
               DENHUC(I+IOFF,J+IOFF) = DENHUC(I+IOFF,J+IOFF) +
     &           2.0D0*EIGVEC(I+IOFF,IOCC+IOFF)*EIGVEC(J+IOFF,IOCC+IOFF)
            END DO
         END DO
         END DO
         IOFF=IOFF+NHUCAO(ISYM)
      END DO
C
C     Print Huckel density matrix
C
      IF (IPRHUC.GE.4) THEN
         WRITE (LUPRI,'(/A)') ' Huckel density matrix:'
         CALL OUTPUT(DENHUC,1,NBASHUC,1,NBASHUC,
     &        NBASHUC,NBASHUC,-1,LUPRI)
      END IF
C
      RETURN
      END
C ========================================================================
      SUBROUTINE HUCFCK(FMAT,F2MAT,DENHUC,NBAS_LS,NBASHUC,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION F2MAT(*),DENHUC(*),WORK(LWORK)
      DIMENSION DNTSKP(8)
#include "cbihr2.h"
#include "mxcent.h"
#include "nuclei.h"
#include "maxorb.h"
#include "huckel.h"
! dcbham.h : NOSMLV and BSS
#include "dcbham.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      CALL QENTER('HUCFCK')
C
      KFREE=1
      LFREE=LWORK
C
      N2BAS_LS = NBAS_LS*NBAS_LS
      F2MAT(1:N2BAS_LS) = 0.0D0
C
      NDMAT = 1
      IREPDM = 0 ! totally symmetric
!     IF (HFXFAC .EQ. D0) THEN
!        IFCTYP = 11 ! symmetric and Coulomb
!     ELSE
!        IFCTYP = 13 ! symmetric and Coulomb+exchange
!     END IF
      IRNTYP = 9 ! Direct calculation of Fock matrices in SO basis
      MAXDIF = 0 ! Undifferentiated integrals
      JATOM  = 0 ! All integrals
      SCRFCK_HUC = -1.0D-7

C      iprhuc = 7
      IF (NBAS_LS .NE. NBASIS_F) THEN
         WRITE (LUPRI,*) 'NBAS_LS .ne. NBASIS_F :',NBAS_LS,NBASIS_F
         CALL QUIT('NBAS_LS .ne. NBASIS_F')
      END IF

      IF (IPRHUC .GE. 7) THEN
         WRITE(LUPRI,*) 'HUCFCK: NBASIS_D, NBASIS_F ',NBASIS_D,NBASIS_F
         WRITE(LUPRI,'(/A)')
     &  ' 2-electron Fock matrix from Huckel density before TWOINT call'
         CALL OUTPUT(F2MAT,1,NBASIS_F,1,NBASIS_F,
     &                    NBASIS_F,NBASIS_F,-1,LUPRI)
      END IF

C
C     Set DMAT dimension in her2*.F routines in nuclei.h
C
      NBASIS_D_save = NBASIS_D
      NBASIS_D = NBASHUC
C
      IPRTWO_TMP = IPRTWO
C
      DO J2TYP = 5,6
      ! 5: (HH|**) integrals (Coulomb); 6: (H*|H*) integrals (exchange)
         IF (NOSMLV.OR.BSS.OR.LEVYLE) THEN
            WRITE(LUPRI,'(/A)') 'calling TWOINT with only LL'
            I2TYP = J2TYP + 100
         ELSE
            I2TYP = J2TYP
         END IF

         IF (J2TYP .EQ. 5) THEN
            IFCTYP = 19 ! symmetric(10) and Huckel(8) and Coulomb(1)
         ELSE
            IFCTYP = 18 ! symmetric(10) and Huckel(8) and special exchange code(0)
         END IF
         DNTSKP(:) = 0.0D0
C         write (lupri,*) 'DENHUC, I2TYP =',I2TYP
C         N2 = NBASIS_D**2
C         CALL PRMGN(N2,DENHUC,1,8,LUPRI)

         CALL GETTIM(CPU1,WALL1)
         CALL TWOINT(WORK,LWORK,F2MAT,DENHUC,NDMAT,IREPDM,IFCTYP,
     &          DUM,IDUM,IDUM,IRNTYP,MAXDIF,JATOM,.TRUE.,
     &          .TRUE.,.FALSE.,TKTIME,IPRTWO_TMP,IPRNTA,IPRNTB,IPRNTC,
     &          IPRNTD,RTNTWO,IDUM,I2TYP,ICEDIF,SCRFCK_HUC,
     &          GABRAO,DMRAO,DMRSO,DNTSKP,.TRUE.,.false.,IDUM,DUM)
         CALL GETTIM(CPU2,WALL2)
         CPU  = CPU2  - CPU1
         WALL = WALL2 - WALL1
         IF (J2TYP .EQ. 5) THEN
            CALL SCRSTA('SOfock:HT Coulomb  CPU ',DNTSKP,CPU)
            CALL SCRSTA('SOfock:HT Coulomb  WALL',DNTSKP,WALL)
         ELSE
            CALL SCRSTA('SOfock:HT Exchange CPU ',DNTSKP,CPU)
            CALL SCRSTA('SOfock:HT Exchange WALL',DNTSKP,WALL)
         END IF
!DBG        IF (IPRHUC .GE. 7) THEN
            IF (J2TYP .EQ. 5) THEN
               WRITE(LUPRI,'(/A)')
     &         ' 2-electron 2J matrix from Huckel density'
            ELSE
               WRITE(LUPRI,'(/A)')
     &         ' 2-electron (2J-K) matrix from Huckel density'
            END IF
            N2 = NBASIS_F**2
            CALL PRMGN(N2,F2MAT,1,8,LUPRI)
!DBG
         IF (IPRHUC .GE. 7) THEN ! moved to here for DBG
!DBG
            CALL OUTPUT(F2MAT,1,NBASIS_F,1,NBASIS_F,
     &                  NBASIS_F,NBASIS_F,-1,LUPRI)
         END IF
         CALL FLSHFO(LUPRI)
      END DO
C
C
      NBASIS_D = NBASIS_D_save
C
C     Reindex to sorted basis  (based on code in TWOFCK)
C
      CALL BUTOBS(F2MAT,1,WORK(KFREE),LFREE) ! NZ set to 1 here, only real part
C
      IF (NZ .LT. 4) THEN
         IREP = 0
         IZ = 1
         IQ = IPQTOQ(IZ,IREP)
          CALL Q2BPHASE('F', IQ,1,F2MAT)
      END IF
C     Symmetrize Fock matrix
C        F2MAT(I,J) = F2MAT(I,J) + F2MAT(J,I)
C        F2MAT(J,I) = F2MAT(I,J)
      CALL FULMAT('S',NBAS_LS,NBAS_LS,F2MAT)
C
C     Check for NAN elements in the Fock matrix
C
      N_NAN_F2MAT = N_NAN(N2BAS_LS,F2MAT)
      IF (N_NAN_F2MAT .GT. 0) IPRHUC = 6
      IF (IPRHUC.GE.6) THEN
         WRITE(LUPRI,'(/A)')
     &      ' 2-electron Fock matrix from Huckel density'
         IF (N_NAN_F2MAT .GT. 0) WRITE (LUPRI,'(/A,I5)')
     &   ' FATAL ERROR: number of NANs in matrix is > 0 :',N_NAN_F2MAT
         CALL OUTPUT(F2MAT,1,NBAS_LS,1,NBAS_LS,
     &               NBAS_LS,NBAS_LS,-1,LUPRI)
         IF (N_NAN_F2MAT .GT. 0) CALL QUIT('Hukel F2MAT has NANs')
         WRITE(LUPRI,'(/A,I2)') ' 1-electron Fock matrix; NZ =',NZ
         CALL PRQMAT(FMAT,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &               NZ,IPQTOQ(1,IREP),LUPRI)
      END IF
C
C     Add Huckel (2J-K)-Matrix to one-electron Fock matrix and return
C
      CALL DAXPY(N2BAS_LS,1.0D0,F2MAT,1,FMAT,1)
      IF (IPRHUC.GE.5) THEN
         WRITE(LUPRI,'(/A)') ' Total Fock matrix from Huckel density'
         CALL PRQMAT(FMAT,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &               NZ,IPQTOQ(1,IREP),LUPRI)
      END IF
C
      CALL QEXIT('HUCFCK')
      RETURN
      END
C ========================================================================
      SUBROUTINE HUCFUN(IZ,IQECP,L,ISHELL,ORBENERGY)
C
C panor 2010:
C     Orbital energies are taken from the relativistic DK3 Hartree-Fock
C     calculations by Nakajima and Hirao in J. Chem. Phys. 116, 8270 (2002)
C
C     Input:
C     IZ     = atomic number
C     IQECP  = number of core electrons described with ECP
C     L      = l_quantum_number + 1
C     ISHELL = orbital number for this L
C
C     Output:
C     ORBENERGY = orbital energy
C
#include <implicit.h>
#include <priunit.h>
C
      INTEGER IOFF(4)
      SAVE    IOFF
      DATA (IOFF(I),I=1,4)/0,7,12,16/
C
C     Ionization potentials for shells:
C     1s,2s,3s,4s,5s,6s,7s, 2p,3p,4p,5p,6p, 3d,4d,5d,6d, 4f,5f
      REAL*8 IP(18,103)
      SAVE   IP
      DATA ((IP(I,IZ),I=1,18),IZ=1,54) /
C  1 HYDROGEN
     &   0.5000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C  2 HELIUM
     &   9.1763300D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C  3 LITHIUM
     &   2.4778426D+00,  1.9633570D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C  4 BERYLLIUM
     &   4.7331763D+00,  3.0930520D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C  5 BORON
     &   7.6969897D+00,  4.9488260D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   3.0972950D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C  6 CARBON
     &   1.1329485D+01,  7.0607700D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   4.3314680D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C  7 NITROGEN
     &   1.5637118D+01,  9.4628230D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   5.6728190D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C  8 OXYGEN
     &   2.0683346D+01,  1.2461301D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   6.3147880D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C  9 FLUORINE
     &   2.6407408D+01,  1.5756466D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   7.2939390D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 10 NEON
     &   3.2811505D+01,  1.9353969D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   8.4960060D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 11 SODIUM
     &   4.0537168D+01,  2.8048982D+00,  1.8231850D-01,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   1.5169539D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0555560D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 12 MAGNESIUM
     &   4.9116404D+01,  3.7795227D+00,  2.5335200D-01,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   2.2804481D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0555560D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 13 ALUMINUM
     &   5.8620210D+01,  4.9283020D+00,  3.9419310D-01,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   3.2162907D+00,  2.0945310D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0555560D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 14 SILICON
     &   6.8976259D+01,  6.1821733D+00,  5.4139970D-01,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   4.2541924D+00,  2.9650420D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0555560D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 15 PHOSPHORUS
     &   8.0189737D+01,  7.5471324D+00,  6.9911100D-01,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   5.3995246D+00,  3.9100870D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0555560D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 16 SULFUR
     &   9.2294186D+01,  9.0535425D+00,  8.8378980D-01,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   6.6818383D+00,  4.3658190D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0555560D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 17 CHLORINE
     &   1.0525968D+02,  1.0673425D+01,  1.0791723D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   8.0728916D+00,  5.0553980D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0555560D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 18 ARGON
     &   1.1908910D+02,  1.2408737D+01,  1.2862131D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   9.5741405D+00,  5.9017130D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0555560D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 19 POTASSIUM
     &   1.3413515D+02,  1.4601125D+01,  1.7608585D+00,  1.4757160D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   1.1524240D+01,  9.5308690D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0555560D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 20 CALCIUM
     &   1.5011169D+02,  1.6963414D+01,  2.2615803D+00,  1.9625650D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   1.3637123D+01,  1.3388102D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   0.0555560D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 21 SCANDIUM
     &   1.6682374D+02,  1.9260995D+01,  2.5910602D+00,  2.1134960D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   1.5684377D+01,  1.5747561D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   3.3569200D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 22 TITANIUM
     &   1.8440142D+02,  2.1649433D+01,  2.9050657D+00,  2.2236360D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   1.7816808D+01,  1.7970860D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   4.3153800D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 23 VANADIUM
     &   2.0286900D+02,  2.4155276D+01,  3.2240602D+00,  2.3251490D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   2.0059761D+01,  2.0233451D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   4.9930180D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 24 CHROMIUM
     &   2.2203224D+02,  2.6559395D+01,  3.3410693D+00,  2.2612280D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   2.2197327D+01,  2.0621049D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   3.6666000D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 25 MANGANESE
     &   2.4248736D+02,  2.9526034D+01,  3.8805572D+00,  2.5056910D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   2.4881102D+01,  2.4892543D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   6.2561170D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 26 IRON
     &   2.6368389D+02,  3.2435920D+01,  4.2478864D+00,  2.6144240D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   2.7502545D+01,  2.7557483D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   6.3239820D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 27 COBALT
     &   2.8578069D+02,  3.5464669D+01,  4.6194346D+00,  2.7127770D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   3.0233036D+01,  3.0243525D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   6.5948160D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 28 NICKEL
     &   3.0879078D+02,  3.8623182D+01,  5.0019485D+00,  2.8073490D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   3.3082595D+01,  3.3010643D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   6.8934950D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 29 COPPER
     &   3.3248621D+02,  4.1656238D+01,  5.1539554D+00,  2.4450350D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   3.5799825D+01,  3.3612480D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   4.7939400D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 30 ZINC
     &   3.5756214D+02,  4.5328875D+01,  5.7973842D+00,  2.9833480D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   3.9135111D+01,  3.8757117D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   7.6102210D-01,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 31 GALLIUM
     &   3.8370754D+02,  4.9283976D+01,  6.5784162D+00,  4.3364480D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   4.2739785D+01,  4.5227988D+00,  2.0682330D-01,  0.0000000D+00,
     &   0.0000000D+00,
     &   1.1645916D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 32 GERMANIUM
     &   4.1083469D+02,  5.3431907D+01,  7.4032342D+00,  5.6619250D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   4.6522609D+01,  5.2075863D+00,  2.8588520D-01,  0.0000000D+00,
     &   0.0000000D+00,
     &   1.5983100D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 33 ARSENIC
     &   4.3895166D+02,  5.7775659D+01,  8.2742088D+00,  7.0282540D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   5.0485752D+01,  5.9330830D+00,  3.6805650D-01,  0.0000000D+00,
     &   0.0000000D+00,
     &   2.0670723D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 34 SELENIUM
     &   4.6808785D+02,  6.2339217D+01,  9.2137342D+00,  8.5939150D-01,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   5.4652370D+01,  6.7217583D+00,  4.0148190D-01,  0.0000000D+00,
     &   0.0000000D+00,
     &   2.5943392D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 35 BROMINE
     &   4.9822106D+02,  6.7096433D+01,  1.0195478D+01,  1.0205021D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   5.8995383D+01,  7.5476795D+00,  4.5584830D-01,  0.0000000D+00,
     &   0.0000000D+00,
     &   3.1542530D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 36 KRYPTON
     &   5.2935756D+02,  7.2048596D+01,  1.1220030D+01,  1.1872343D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   6.3515018D+01,  8.4114233D+00,  5.2303070D-01,  0.0000000D+00,
     &   0.0000000D+00,
     &   3.7474281D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 37 RUBIDIUM
     &   5.6177668D+02,  7.7467270D+01,  1.2554999D+01,  1.5655773D+00,
     &   1.3964420D-01,  0.0000000D+00,  0.0000000D+00,
     &   6.8481394D+01,  9.5783613D+00,  8.0778770D-01,  0.0000000D+00,
     &   0.0000000D+00,
     &   4.6401961D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C 38 STRONTIUM
     &   5.9523714D+02,  8.3107116D+01,  1.3953636D+01,  1.9476221D+00,
     &   1.7957030D-01,  0.0000000D+00,  0.0000000D+00,
     &   7.3648183D+01,  1.0802518D+01,  1.0943405D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   5.5862651D+00,  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C YTTRIUM
     &   6.2964887D+02,  8.8863962D+01,  1.5310553D+01,  2.2373520D+00,
     &   1.9908540D-01,  0.0000000D+00,  0.0000000D+00,
     &   7.8911807D+01,  1.1980148D+01,  1.3031848D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   6.4836320D+00,  2.3024790D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C ZIRCONIUM
     &   6.6507074D+02,  9.4796219D+01,  1.6685109D+01,  2.5044024D+00,
     &   2.1108970D-01,  0.0000000D+00,  0.0000000D+00,
     &   8.4328586D+01,  1.3169642D+01,  1.4937698D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   7.3896027D+00,  3.1455910D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C NIOBIUM
     &   7.0141304D+02,  1.0080375D+02,  1.7975420D+01,  2.6508427D+00,
     &   2.2733980D-01,  0.0000000D+00,  0.0000000D+00,
     &   8.9798138D+01,  1.4269836D+01,  1.5776316D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   8.2037553D+00,  2.8628790D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C MOLYBDENUM
     &   7.3889498D+02,  1.0711324D+02,  1.9408290D+01,  2.8962706D+00,
     &   2.3605980D-01,  0.0000000D+00,  0.0000000D+00,
     &   9.5543408D+01,  1.5505012D+01,  1.7494945D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   9.1473609D+00,  3.4115660D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C TECHNETIUM
     &   7.7758698D+02,  1.1379010D+02,  2.1047227D+01,  3.2970117D+00,
     &   2.3731090D-01,  0.0000000D+00,  0.0000000D+00,
     &   1.0162811D+02,  1.6938278D+01,  2.0622225D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   1.0282901D+01,  5.1314520D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C RUTHENIUM
     &   8.1708385D+02,  1.2040327D+02,  2.2457360D+01,  3.4378118D+00,
     &   2.3644000D-01,  0.0000000D+00,  0.0000000D+00,
     &   1.0762367D+02,  1.8137922D+01,  2.1405352D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   1.1182652D+01,  3.9224760D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C RHODIUM
     &   8.5777952D+02,  1.2736048D+02,  2.4048209D+01,  3.7114698D+00,
     &   2.3577850D-01,  0.0000000D+00,  0.0000000D+00,
     &   1.1393226D+02,  1.9509505D+01,  2.3380927D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   1.2247709D+01,  4.2735100D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C PALLADIUM
     &   8.9940787D+02,  1.3437670D+02,  2.5530582D+01,  3.8393104D+00,
     &   0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &   1.2026980D+02,  2.0765700D+01,  2.4012471D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   1.3193121D+01,  3.2705970D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C SILVER
     &   9.4241973D+02,  1.4191551D+02,  2.7370977D+01,  4.2707059D+00,
     &   2.3433670D-01,  0.0000000D+00,  0.0000000D+00,
     &   1.2709447D+02,  2.2369381D+01,  2.7417866D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   1.4477433D+01,  5.0960310D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C CADMIUM
     &   9.8657558D+02,  1.4972261D+02,  2.9308660D+01,  4.7440428D+00,
     &   2.8095130D-01,  0.0000000D+00,  0.0000000D+00,
     &   1.3415216D+02,  2.4060856D+01,  3.1169689D+00,  0.0000000D+00,
     &   0.0000000D+00,
     &   1.5842760D+01,  7.2022290D-01,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C INDIUM
     &   1.0319374D+03,  1.5784810D+02,  3.1391156D+01,  5.3058460D+00,
     &   3.9581680D-01,  0.0000000D+00,  0.0000000D+00,
     &   1.4149184D+02,  2.5888659D+01,  3.5780845D+00,  1.9430990D-01,
     &   0.0000000D+00,
     &   1.7337610D+01,  1.0107736D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C TIN
     &   1.0784371D+03,  1.6622060D+02,  3.3544893D+01,  5.8808385D+00,
     &   5.0645700D-01,  0.0000000D+00,  0.0000000D+00,
     &   1.4903999D+02,  2.7778130D+01,  4.0478971D+00,  2.6265080D-01,
     &   0.0000000D+00,
     &   1.8886852D+01,  1.3062919D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C ANTIMONY
     &   1.1260928D+03,  1.7484856D+02,  3.5775907D+01,  6.4749753D+00,
     &   6.1912520D-01,  0.0000000D+00,  0.0000000D+00,
     &   1.5680282D+02,  2.9734678D+01,  4.5326548D+00,  3.3277960D-01,
     &   0.0000000D+00,
     &   2.0495703D+01,  1.6138421D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C TELLURIUM
     &   1.1749309D+03,  1.8375374D+02,  3.8103799D+01,  7.1070697D+00,
     &   7.4639650D-01,  0.0000000D+00,  0.0000000D+00,
     &   1.6480050D+02,  3.1777418D+01,  5.0511518D+00,  3.5811400D-01,
     &   0.0000000D+00,
     &   2.2183122D+01,  1.9524601D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C IODINE
     &   1.2249486D+03,  1.9291907D+02,  4.0508674D+01,  7.7573211D+00,
     &   8.7615960D-01,  0.0000000D+00,  0.0000000D+00,
     &   1.7301287D+02,  3.3885781D+01,  5.5837714D+00,  4.0175090D-01,
     &   0.0000000D+00,
     &   2.3928189D+01,  2.3026219D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C XENON
     &   1.2761501D+03,  2.0234749D+02,  4.2991851D+01,  8.4269386D+00,
     &   1.0094342D+00,  0.0000000D+00,  0.0000000D+00,
     &   1.8144155D+02,  3.6060708D+01,  6.1317507D+00,  4.5617540D-01,
     &   0.0000000D+00,
     &   2.5731540D+01,  2.6654846D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00/
      DATA ((IP(I,IZ),I=1,18),IZ=55,103) /
C CESIUM
     &   1.3287682D+03,  2.1225672D+02,  4.5766901D+01,  9.3270229D+00,
     &   1.3085719D+00,  1.2786920D-01,  0.0000000D+00,
     &   1.9030111D+02,  3.8514373D+01,  6.9044163D+00,  6.8059450D-01,
     &   0.0000000D+00,
     &   2.7805532D+01,  3.2503720D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C BARIUM
     &   1.3826060D+03,  2.2244555D+02,  4.8630473D+01,  1.0253612D+01,
     &   1.6024384D+00,  1.6181880D-01,  0.0000000D+00,
     &   1.9938822D+02,  4.1042957D+01,  7.6982688D+00,  8.9849370D-01,
     &   0.0000000D+00,
     &   2.9945855D+01,  3.8535899D+00,  0.0000000D+00,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C LANTHANUM
     &   1.4376050D+03,  2.3283802D+02,  5.1503656D+01,  1.1128138D+01,
     &   1.8187271D+00,  1.7773870D-01,  0.0000000D+00,
     &   2.0862480D+02,  4.3567712D+01,  8.4363578D+00,  1.0518532D+00,
     &   0.0000000D+00,
     &   3.2073952D+01,  4.3994186D+00,  2.3430360D-01,  0.0000000D+00,
     &   0.0000000D+00,  0.0000000D+00,
C CERIUM
     &   1.4935571D+03,  2.4318973D+02,  5.4091779D+01,  1.1687720D+01,
     &   1.8789082D+00,  1.7900260D-01,  0.0000000D+00,
     &   2.1776807D+02,  4.5795114D+01,  8.8682024D+00,  1.0802555D+00,
     &   0.0000000D+00,
     &   3.3902492D+01,  4.6653173D+00,  2.6366290D-01,  0.0000000D+00,
     &   5.4973260D-01,  0.0000000D+00,
C PRASEODYMIUM
     &   1.5504790D+03,  2.5350754D+02,  5.6399175D+01,  1.1939763D+01,
     &   1.7930811D+00,  1.6981450D-01,  0.0000000D+00,
     &   2.2682205D+02,  4.7728795D+01,  9.0009123D+00,  9.9209910D-01,
     &   0.0000000D+00,
     &   3.5434542D+01,  4.6564865D+00,  0.0000000D+00,  0.0000000D+00,
     &   3.9570020D-01,  0.0000000D+00,
C NEODIUM
     &   1.6088973D+03,  2.6436389D+02,  5.9082413D+01,  1.2494972D+01,
     &   1.8488979D+00,  1.7210210D-01,  0.0000000D+00,
     &   2.3634802D+02,  5.0020776D+01,  9.4214663D+00,  1.0168170D+00,
     &   0.0000000D+00,
     &   3.7309024D+01,  4.9081547D+00,  0.0000000D+00,  0.0000000D+00,
     &   4.3660430D-01,  0.0000000D+00,
C PROMETHIUM
     &   1.6686293D+03,  2.7554574D+02,  6.1889780D+01,  1.3108314D+01,
     &   1.9210033D+00,  1.7499250D-01,  0.0000000D+00,
     &   2.4613327D+02,  5.2419949D+01,  9.8931350D+00,  1.0526833D+00,
     &   0.0000000D+00,
     &   3.9279570D+01,  5.2034363D+00,  0.0000000D+00,  0.0000000D+00,
     &   3.2093410D-01,  0.0000000D+00,
C SAMARIUM
     &   1.7295280D+03,  2.8688118D+02,  6.4628022D+01,  1.3624422D+01,
     &   1.9583371D+00,  1.7641300D-01,  0.0000000D+00,
     &   2.5600527D+02,  5.4733567D+01,  1.0268476D+01,  1.0631092D+00,
     &   0.0000000D+00,
     &   4.1155830D+01,  5.4096165D+00,  0.0000000D+00,  0.0000000D+00,
     &   4.9398570D-01,  0.0000000D+00,
C EUROPIUM
     &   1.7917731D+03,  2.9854607D+02,  6.7486597D+01,  1.4194129D+01,
     &   2.0105233D+00,  1.7851890D-01,  0.0000000D+00,
     &   2.6613304D+02,  5.7148818D+01,  1.0690192D+01,  1.0837231D+00,
     &   0.0000000D+00,
     &   4.3121782D+01,  5.6547841D+00,  0.0000000D+00,  0.0000000D+00,
     &   5.3084210D-01,  0.0000000D+00,
C GADOLINIUM
     &   1.8556002D+03,  3.1079456D+02,  7.0750551D+01,  1.5095701D+01,
     &   2.2289281D+00,  1.9471160D-01,  0.0000000D+00,
     &   2.7676563D+02,  5.9948756D+01,  1.1427918D+01,  1.2369285D+00,
     &   0.0000000D+00,
     &   4.5455126D+01,  6.1902134D+00,  2.4301620D-01,  0.0000000D+00,
     &   8.3464320D-01,  0.0000000D+00,
C TERBIUM
     &   1.9202677D+03,  3.2279297D+02,  7.3473389D+01,  1.5420084D+01,
     &   2.1329554D+00,  1.8313890D-01,  0.0000000D+00,
     &   2.8707442D+02,  6.2190214D+01,  1.1601357D+01,  1.1357321D+00,
     &   0.0000000D+00,
     &   4.7228215D+01,  6.1974946D+00,  0.0000000D+00,  0.0000000D+00,
     &   5.0546380D-01,  0.0000000D+00,
C DYSPROSIUM
     &   1.9865163D+03,  3.3534713D+02,  7.6562273D+01,  1.6042323D+01,
     &   2.1922127D+00,  1.8531010D-01,  0.0000000D+00,
     &   2.9785329D+02,  6.4775419D+01,  1.2058084D+01,  1.1593811D+00,
     &   0.0000000D+00,
     &   4.9327461D+01,  6.4652273D+00,  0.0000000D+00,  0.0000000D+00,
     &   5.0661730D-01,  0.0000000D+00,
C HOLMIUM
     &   2.0541384D+03,  3.4820485D+02,  7.9726180D+01,  1.6678746D+01,
     &   2.2527268D+00,  1.8755100D-01,  0.0000000D+00,
     &   3.0884763D+02,  6.7413375D+01,  1.2522797D+01,  1.1832946D+00,
     &   0.0000000D+00,
     &   5.1466342D+01,  6.7368886D+00,  0.0000000D+00,  0.0000000D+00,
     &   5.0529690D-01,  0.0000000D+00,
C ERBIUM
     &   2.1232078D+03,  3.6142570D+02,  8.3015838D+01,  1.7383451D+01,
     &   2.0300293D+00,  1.4558500D-01,  0.0000000D+00,
     &   3.2011477D+02,  7.0155518D+01,  1.3055072D+01,  1.2345478D+00,
     &   0.0000000D+00,
     &   5.3696569D+01,  7.0746198D+00,  0.0000000D+00,  0.0000000D+00,
     &   5.6416410D-01,  0.0000000D+00,
C THULIUM
     &   2.1938800D+03,  3.7514923D+02,  8.6563671D+01,  1.8279018D+01,
     &   1.8983011D+00,  1.7221560D-01,  0.0000000D+00,
     &   3.3179170D+02,  7.3133749D+01,  1.3774870D+01,  1.4489966D+00,
     &   0.0000000D+00,
     &   5.6150278D+01,  7.5964299D+00,  0.0000000D+00,  0.0000000D+00,
     &   8.0827660D-01,  0.0000000D+00,
C YTTERBIUM
     &   2.2657707D+03,  3.8899821D+02,  9.0005070D+01,  1.9006636D+01,
     &   1.8286991D+00,  1.6968650D-01,  0.0000000D+00,
     &   3.4349172D+02,  7.5977647D+01,  1.4312973D+01,  1.5129302D+00,
     &   0.0000000D+00,
     &   5.8453589D+01,  7.9302902D+00,  0.0000000D+00,  0.0000000D+00,
     &   8.7179470D-01,  0.0000000D+00,
C LUTECIUM
     &   2.3389472D+03,  4.0305891D+02,  9.3471108D+01,  1.9681457D+01,
     &   2.6959009D+00,  2.1838510D-01,  0.0000000D+00,
     &   3.5528697D+02,  7.8810386D+01,  1.4763357D+01,  1.4341766D+00,
     &   0.0000000D+00,
     &   6.0721397D+01,  8.1455680D+00,  1.8507650D-01,  0.0000000D+00,
     &   8.1663660D-01,  0.0000000D+00,
C HAFNIUM
     &   2.4139886D+03,  4.1783256D+02,  9.7384351D+01,  2.0732918D+01,
     &   2.9407691D+00,  2.3318990D-01,  0.0000000D+00,
     &   3.6769204D+02,  8.2068926D+01,  1.5596907D+01,  1.6008557D+00,
     &   0.0000000D+00,
     &   6.3402262D+01,  8.7427570D+00,  2.6065940D-01,  0.0000000D+00,
     &   1.1340000D+00,  0.0000000D+00,
C TANTALUM
     &   2.4905572D+03,  4.3297604D+02,  1.0141651D+02,  2.1826819D+01,
     &   3.1879975D+00,  2.4590640D-01,  0.0000000D+00,
     &   3.8034983D+02,  8.5416549D+01,  1.6462491D+01,  1.7678155D+00,
     &   0.0000000D+00,
     &   6.6155310D+01,  9.3642911D+00,  3.1829300D-01,  0.0000000D+00,
     &   1.4712775D+00,  0.0000000D+00,
C TUNGSTEN
     &   2.5686609D+03,  4.4849275D+02,  1.0556669D+02,  2.2961039D+01,
     &   3.4383440D+00,  2.5687330D-01,  0.0000000D+00,
     &   3.9326038D+02,  8.8851445D+01,  1.7357816D+01,  1.9357852D+00,
     &   0.0000000D+00,
     &   6.8978282D+01,  1.0008188D+01,  3.7026300D-01,  0.0000000D+00,
     &   1.8273602D+00,  0.0000000D+00,
C RHENIUM
     &   2.6483259D+03,  4.6438722D+02,  1.0983307D+02,  2.4132247D+01,
     &   3.6908919D+00,  2.6636950D-01,  0.0000000D+00,
     &   4.0642189D+02,  9.2370163D+01,  1.8279217D+01,  2.1041488D+00,
     &   0.0000000D+00,
     &   7.1867173D+01,  1.0671041D+01,  4.3130270D-01,  0.0000000D+00,
     &   2.1994323D+00,  0.0000000D+00,
C OSMIUM
     &   2.7295922D+03,  4.8068777D+02,  1.1423812D+02,  2.5361177D+01,
     &   3.9652004D+00,  2.7978360D-01,  0.0000000D+00,
     &   4.1985769D+02,  9.5993870D+01,  1.9247017D+01,  2.2897415D+00,
     &   0.0000000D+00,
     &   7.4842626D+01,  1.1372972D+01,  4.4836040D-01,  0.0000000D+00,
     &   2.6080211D+00,  0.0000000D+00,
C IRIDIUM
     &   2.8124821D+03,  4.9738035D+02,  1.1875994D+02,  2.6624190D+01,
     &   4.2406936D+00,  2.9145190D-01,  0.0000000D+00,
     &   4.3354452D+02,  9.9698370D+01,  2.0237065D+01,  2.4746100D+00,
     &   0.0000000D+00,
     &   7.7879862D+01,  1.2090085D+01,  4.8102470D-01,  0.0000000D+00,
     &   3.0292765D+00,  0.0000000D+00,
C PLATINUM
     &   2.8968928D+03,  5.1435176D+02,  1.2328073D+02,  2.7801754D+01,
     &   4.4021979D+00,  2.8818200D-01,  0.0000000D+00,
     &   4.4736591D+02,  1.0336558D+02,  2.1130921D+01,  2.5593436D+00,
     &   0.0000000D+00,
     &   8.0860603D+01,  1.2704922D+01,  4.1116500D-01,  0.0000000D+00,
     &   3.3454303D+00,  0.0000000D+00,
C GOLD
     &   2.9830810D+03,  5.3184518D+02,  1.2803665D+02,  2.9129500D+01,
     &   4.6788713D+00,  2.8671500D-01,  0.0000000D+00,
     &   4.6155601D+02,  1.0722820D+02,  2.2160972D+01,  2.7415456D+00,
     &   0.0000000D+00,
     &   8.4016132D+01,  1.3447829D+01,  4.5096990D-01,  0.0000000D+00,
     &   3.7876371D+00,  0.0000000D+00,
C MERCURY
     &   3.0710951D+03,  5.4990148D+02,  1.3306282D+02,  3.0640624D+01,
     &   5.0998297D+00,  3.2631650D-01,  0.0000000D+00,
     &   4.7615042D+02,  1.1131980D+02,  2.3359577D+01,  3.0456928D+00,
     &   0.0000000D+00,
     &   8.7379301D+01,  1.4350551D+01,  6.0350640D-01,  0.0000000D+00,
     &   4.3877293D+00,  0.0000000D+00,
C TALLIUM
     &   3.1609455D+03,  5.6848464D+02,  1.3831353D+02,  3.2286057D+01,
     &   5.6167056D+00,  4.4737130D-01,  0.0000000D+00,
     &   4.9110358D+02,  1.1559339D+02,  2.4679411D+01,  3.4425873D+00,
     &   1.8539970D-01,
     &   9.0902923D+01,  1.5366669D+01,  8.4097790D-01,  0.0000000D+00,
     &   5.0987354D+00,  0.0000000D+00,
C LEAD
     &   3.2525651D+03,  5.8751414D+02,  1.4370144D+02,  3.3976603D+01,
     &   6.1394866D+00,  5.6178620D-01,  0.0000000D+00,
     &   5.0632797D+02,  1.1995945D+02,  2.6029913D+01,  3.8384250D+00,
     &   2.5071080D-01,
     &   9.4496677D+01,  1.6405124D+01,  1.0732793D+00,  0.0000000D+00,
     &   5.8296970D+00,  0.0000000D+00,
C BISMUTH
     &   3.3460060D+03,  6.0700856D+02,  1.4923596D+02,  3.5719311D+01,
     &   6.6751582D+00,  6.7702800D-01,  0.0000000D+00,
     &   5.2183189D+02,  1.2442476D+02,  2.7417274D+01,  4.2406204D+00,
     &   3.1665500D-01,
     &   9.8166653D+01,  1.7471967D+01,  1.3087916D+00,  0.0000000D+00,
     &   6.5866512D+00,  0.0000000D+00,
C POLONIUM
     &   3.4413190D+03,  6.2699674D+02,  1.5493803D+02,  3.7532951D+01,
     &   7.2418010D+00,  8.0617200D-01,  0.0000000D+00,
     &   5.3763590D+02,  1.2900804D+02,  2.8859517D+01,  4.6672274D+00,
     &   3.3887450D-01,
     &   1.0193091D+02,  1.8585115D+01,  1.5657050D+00,  0.0000000D+00,
     &   7.3874532D+00,  0.0000000D+00,
C ASTATINE
     &   3.5385080D+03,  6.4746800D+02,  1.6079005D+02,  3.9398226D+01,
     &   7.8204747D+00,  9.3665380D-01,  0.0000000D+00,
     &   5.5372207D+02,  1.3368980D+02,  3.0336753D+01,  5.0994513D+00,
     &   3.7814390D-01,
     &   1.0576917D+02,  1.9724409D+01,  1.8254497D+00,  0.0000000D+00,
     &   8.2119125D+00,  0.0000000D+00,
C RADON
     &   3.6375341D+03,  6.6842323D+02,  1.6679368D+02,  4.1316674D+01,
     &   8.4129980D+00,  1.0700253D+00,  0.0000000D+00,
     &   5.7009451D+02,  1.3847225D+02,  3.1850551D+01,  5.5391982D+00,
     &   4.2741270D-01,
     &   1.0968286D+02,  2.0891199D+01,  2.0898461D+00,  0.0000000D+00,
     &   9.0613211D+00,  0.0000000D+00,
C FRANCIUM
     &   3.7387859D+03,  6.9008558D+02,  1.7314247D+02,  4.3474672D+01,
     &   9.2027413D+00,  1.3508737D+00,  1.3150520D-01,
     &   5.8693938D+02,  1.4353955D+02,  3.3584096D+01,  6.1680589D+00,
     &   6.2198460D-01,
     &   1.1385559D+02,  2.2268780D+01,  2.5398767D+00,  0.0000000D+00,
     &   1.0118860D+01,  0.0000000D+00,
C RADIUM
     &   3.8420271D+03,  7.1227148D+02,  1.7965477D+02,  4.5690626D+01,
     &   1.0008405D+01,  1.6244793D+00,  1.6423600D-01,
     &   6.0407812D+02,  1.4871123D+02,  3.5356290D+01,  6.8051944D+00,
     &   8.0726600D-01,
     &   1.1810595D+02,  2.3675601D+01,  2.9946726D+00,  0.0000000D+00,
     &   1.1202970D+01,  0.0000000D+00,
C ACTINIUM
     &   3.9472983D+03,  7.3495985D+02,  1.8629498D+02,  4.7925428D+01,
     &   1.0789982D+01,  1.8511503D+00,  1.8934730D-01,
     &   6.2147287D+02,  1.5394773D+02,  3.7127273D+01,  7.4116498D+00,
     &   9.5464930D-01,
     &   1.2239367D+02,  2.5071942D+01,  3.4155144D+00,  1.7586300D-01,
     &   1.2273720D+01,  0.0000000D+00,
C THORIUM
     &   4.0546073D+03,  7.5817801D+02,  1.9308746D+02,  5.0202678D+01,
     &   1.1571720D+01,  2.0614888D+00,  2.0631020D-01,
     &   6.3914888D+02,  1.5927251D+02,  3.8920020D+01,  8.0108387D+00,
     &   1.0882132D+00,
     &   1.2674125D+02,  2.6480321D+01,  3.8261220D+00,  2.3466100D-01,
     &   1.3353605D+01,  0.0000000D+00,
C PROTOACTINIUM
     &   4.1637953D+03,  7.8162491D+02,  1.9967918D+02,  5.2164855D+01,
     &   1.2001287D+01,  2.0424633D+00,  1.9372240D-01,
     &   6.5677606D+02,  1.6432594D+02,  4.0376568D+01,  8.2644276D+00,
     &   1.0412276D+00,
     &   1.3079105D+02,  2.7543351D+01,  3.9260720D+00,  2.2437520D-01,
     &   1.4088693D+01,  3.1895090D-01,
C URANIUM
     &   4.2751994D+03,  8.0576608D+02,  2.0658050D+02,  5.4323934D+01,
     &   1.2588770D+01,  2.1227648D+00,  1.9614210D-01,
     &   6.7483187D+02,  1.6962025D+02,  4.2007441D+01,  8.6644234D+00,
     &   1.0727001D+00,
     &   1.3505037D+02,  2.8770014D+01,  4.1567116D+00,  2.2863220D-01,
     &   1.4983208D+01,  3.9154130D-01,
C NEPTUNIUM
     &   4.3889338D+03,  8.3050975D+02,  2.1365155D+02,  5.6531057D+01,
     &   1.3183130D+01,  2.2025054D+00,  1.9884120D-01,
     &   6.9317732D+02,  1.7500421D+02,  4.3661057D+01,  9.0640350D+00,
     &   1.1029400D+00,
     &   1.3936817D+02,  3.0008746D+01,  4.3855711D+00,  2.3345040D-01,
     &   1.5886814D+01,  4.4463110D-01,
C PLUTONIUM
     &   4.5047033D+03,  8.5566444D+02,  2.2069877D+02,  5.8594595D+01,
     &   1.3595013D+01,  2.1562143D+00,  1.8590140D-01,
     &   7.1163873D+02,  1.8028850D+02,  4.5147420D+01,  9.2806309D+00,
     &   1.0338403D+00,
     &   1.4355549D+02,  3.1069618D+01,  4.4466969D+00,  0.0000000D+00,
     &   1.6610945D+01,  3.3858430D-01,
C AMERICIUM
     &   4.6230288D+03,  8.8161682D+02,  2.2810202D+02,  6.0887520D+01,
     &   1.4192833D+01,  2.2271542D+00,  1.8824950D-01,
     &   7.3056085D+02,  1.8584216D+02,  4.6835193D+01,  9.6690329D+00,
     &   1.0558336D+00,
     &   1.4797770D+02,  3.2320189D+01,  4.6615839D+00,  0.0000000D+00,
     &   1.7520029D+01,  3.8780340D-01,
C CURIUM
     &   4.7438982D+03,  9.0840491D+02,  2.3590040D+02,  6.3448459D+01,
     &   1.5013287D+01,  2.4415990D+00,  2.0971520D-01,
     &   7.4998456D+02,  1.9170471D+02,  4.8762777D+01,  1.0264431D+01,
     &   1.1886957D+00,
     &   1.5267311D+02,  3.3798329D+01,  5.0631826D+00,  2.1081330D-01,
     &   1.8651807D+01,  5.9387340D-01,
C BERKELIUM
     &   4.8671472D+03,  9.3570097D+02,  2.4370925D+02,  6.5883991D+01,
     &   1.5663572D+01,  2.5336977D+00,  2.1558730D-01,
     &   7.6954371D+02,  1.9748232D+02,  5.0535274D+01,  1.0686472D+01,
     &   1.2241876D+00,
     &   1.5724953D+02,  3.5109613D+01,  5.3043580D+00,  1.9036170D-01,
     &   1.9614604D+01,  6.0361210D-01,
C CALIFORNIUM
     &   4.9926504D+03,  9.6345617D+02,  2.5148380D+02,  6.8150306D+01,
     &   1.6101464D+01,  2.4674766D+00,  1.9636390D-01,
     &   7.8919818D+02,  2.0313209D+02,  5.2109368D+01,  1.0893826D+01,
     &   1.1335560D+00,
     &   1.6166294D+02,  3.6210536D+01,  5.3446302D+00,  0.0000000D+00,
     &   2.0364688D+01,  4.3140630D-01,
C EINSTEINIUM
     &   5.1210032D+03,  9.9212743D+02,  2.5966985D+02,  7.0686112D+01,
     &   1.6759872D+01,  2.5469395D+00,  1.9896550D-01,
     &   8.0935839D+02,  2.0908783D+02,  5.3918966D+01,  1.1304820D+01,
     &   1.1563260D+00,
     &   1.6634342D+02,  3.7533641D+01,  5.5693442D+00,  0.0000000D+00,
     &   2.1331864D+01,  4.5108880D-01,
C FERMIUM
     &   5.2519378D+03,  1.0215163D+03,  2.6806093D+02,  7.3284689D+01,
     &   1.7433290D+01,  2.6279385D+00,  2.0163950D-01,
     &   8.2983525D+02,  2.1514526D+02,  5.5758386D+01,  1.1720388D+01,
     &   1.1788286D+00,
     &   1.7108678D+02,  3.8873268D+01,  5.7951552D+00,  0.0000000D+00,
     &   2.2311652D+01,  4.6927540D-01,
C MENDELEVIUM
     &   5.3857182D+03,  1.0516798D+03,  2.7666864D+02,  7.5947020D+01,
     &   1.8119838D+01,  2.7087443D+00,  2.0425220D-01,
     &   8.5062772D+02,  2.2130155D+02,  5.7624376D+01,  1.2137510D+01,
     &   1.1994795D+00,
     &   1.7588898D+02,  4.0225772D+01,  6.0191440D+00,  0.0000000D+00,
     &   2.3300313D+01,  4.9025550D-01,
C NOBELIUM
     &   5.5223333D+03,  1.0826280D+03,  2.8549565D+02,  7.8673810D+01,
     &   1.8819634D+01,  2.7894043D+00,  2.0677370D-01,
     &   8.7173936D+02,  2.2755751D+02,  5.9516956D+01,  1.2556062D+01,
     &   1.2184736D+00,
     &   1.8074961D+02,  4.1590745D+01,  6.2410949D+00,  0.0000000D+00,
     &   2.4297369D+01,  5.1543850D-01,
C LAWRENCIUM
     &   5.6621151D+03,  1.1146484D+03,  2.9482186D+02,  8.1740758D+01,
     &   1.9804209D+01,  3.0688175D+00,  2.4442330D-01,
     &   8.9343257D+02,  2.3418636D+02,  6.1708942D+01,  1.3238803D+01,
     &   1.3918244D+00,
     &   1.8593924D+02,  4.3239915D+01,  6.7056161D+00,  1.4983630D-01,
     &   2.5573120D+01,  7.6717180D-01/
C
      IF ((IZ.LE.0) .OR. (IZ.GE.104)) THEN
         CALL QUIT('HUCFUN: Invalid element')
      END IF
C
      IPOS=IOFF(L)+ISHELL
C
      IF (IQECP.GT.0) THEN
C     Calculate offset for this l quantum number when IQECP core
C     electrons are included in the ECP
         CALL ECP_LCORE(IQECP,L,IECP)
         IPOS=IPOS+IECP
      END IF
      ORBENERGY = -IP(IPOS,IZ)
C
      RETURN
      END
C ========================================================================
      SUBROUTINE DCOPYRECT(NRA,NCA,A,KIRA,KICA,KFRA,KFCA,
     &                     NRB,NCB,B,KIRB,KICB,KFRB,KFCB)
C
C     Copy A(KIRA:KFRA, KICA:KFCA) to B(KIRB:KFRB, KICB:KFCB)
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION A(NRA,NCA),B(NRB,NCB)
C
      IF (NRA.LE.0 .OR. NCA.LE.0 .OR. NRB.LE.0 .OR. NCB.LE.0 .OR.
     &     KIRA.LE.0 .OR. KFRA.LT.KIRA .OR. KFRA.GT.NRA .OR.
     &     KICA.LE.0 .OR. KFCA.LT.KICA .OR. KFCA.GT.NCA .OR.
     &     KIRB.LE.0 .OR. KFRB.LT.KIRB .OR. KFRB.GT.NRB .OR.
     &     KICB.LE.0 .OR. KFCB.LT.KICB .OR. KFCB.GT.NCB .OR.
     &     (KFRA-KIRA+1)*(KFCA-KICA+1).NE.(KFRA-KIRA+1)*(KFCA-KICA+1)
     &     ) THEN
         WRITE(LUPRI,'(6(/A,2I5))')
     &        ' A full dimension :',NRA,NCA,
     &        ' A upper left pos :',KIRA,KICA,
     &        ' A lower right pos:',KFRA,KFCA,
     &        ' B full dimension :',NRB,NCB,
     &        ' B upper left pos :',KIRB,KICB,
     &        ' B lower right pos:',KFRB,KFCB
         CALL QUIT('DCOPYRECT: Illegal matrix specification')
      END IF
C
      IROWB=KIRB
      ICOLB=KICB
      DO ICOLA=KICA,KFCA
         DO IROWA=KIRA,KFRA
            B(IROWB,ICOLB)=A(IROWA,ICOLA)
            IROWB=IROWB+1
            IF (IROWB.GT.KFRB) THEN
               IROWB=KIRB
               ICOLB=ICOLB+1
            END IF
         END DO
      END DO
C
      RETURN
      END
C ==  end of huckel.F  ===================================================
