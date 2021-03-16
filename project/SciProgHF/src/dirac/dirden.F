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

C FILE    : dirden.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fckden */
      SUBROUTINE FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Get backtransformed (contravariant) density matrices.
C     Based on DENSTY
C
C     Written by J. Thyssen
C     Last revision: 1997/09/22 - jth
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
      DIMENSION DCAO(*),DVAO(*),CMO(*),DV(*)
      DIMENSION WORK(*)
      LOGICAL GETDC,GETDV
C
      CALL QENTER('FCKDEN')
      IF ( IPRINT .GE. 9 ) CALL HEADER('Output from FCKDEN',-1)
C
C     Get active density matrix
C
      IF (GETDV) THEN
C         CALL DENSDV(DVAO,CMO,DV,IPRINT,WORK,LWORK)
         CALL GENDEN(DVAO,CMO,2,IPRINT)
         CALL DSCAL(N2BBASXQ,0.50D00,DVAO,1)
      END IF
C
C     Get inactive density matrix
C   
      CALL DZERO(DCAO,N2BBASXQ)
      IF (GETDC) THEN
         CALL DENSTY(DCAO,CMO,IPRINT)
      END IF
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('FCKDEN')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck densdv */
      SUBROUTINE DENSDV(DVAO,CMO,DV,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Get backtransformed (contravariant) active density.
C     Based on DALTON/FCKDEN
C
C     Written by J. Thyssen
C     Last revision: 1997/09/22 - jth
C                    2002/09/26 - tf : minor cleanup
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      DIMENSION DVAO(*),CMO(*),DV(NASHT,NASHT,NZ)
      DIMENSION WORK(LWORK)
      PARAMETER ( D1 = 1.0 D00, D0 = 0.0 D00 )
C
      CALL QENTER('DENSDV')
#include "memint.h"
C
      IF ( IPRINT .GE. 10 ) CALL HEADER('Output from DENSDV',-1) 
C
C     Get active density matrix
C
      IF (NASHT .EQ. 0) THEN
         WRITE(LUPRI,*) 'WARNING!!! Active density matrix'//
     &                  ' requested but there are no active electrons.'
         WRITE(LUPRI,*) 'So DVAO = 0'
         GOTO 1000
      END IF
      CALL DZERO(DVAO,N2BBASXQ)
C 
C     Allocate memory for temp. DVMO matrix
C
      CALL MEMGET('REAL',KTEMP,N2ORBXQ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KTEMP),N2ORBXQ)
C
C
      CALL DVSCT(DV,WORK(KTEMP))
C
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('DENSDV "Expanded" DVMO',-1)
         CALL PRQMAT(WORK(KTEMP),NORBT,NORBT,NORBT,NORBT,
     &        NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      DO I = 1,NFSYM
         IF ( NASH(I) .GT. 0 ) THEN
            CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &               NORB(I),NORB(I),
     &               DVAO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),
     &                   NZ,IPQTOQ(1,0),
     &               WORK(KTEMP+I2ORBX(I,I)),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &               CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &               CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &               WORK(KFREE),LFREE,IPRINT)
         END IF
      END DO
C
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('DENSDV: Active density matrix',-1)
         DO I = 1,NFSYM
            WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
            CALL PRQMAT(DVAO(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
C
      CALL MEMREL('DENSDV',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
 1000 CONTINUE
      CALL QEXIT('DENSDV')
      CALL FLSHFO(LUPRI)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dgetdv */
      SUBROUTINE DGETDV(DV,NZ)
C*****************************************************************************
C
C     Get active density matrix 
C     This is a dummy routine; it just returns DV=[[1.0]]
C
C     Written by J. Thyssen
C     Last revision: 1997/10/22 - jth
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION DV(*)
C
      DV(1) = 0.50D00
      DO I = 2,NZ
         DV(I) = 0.00D00
      END DO
      RETURN
      END
C
C Dummy routine: writtes matrix to file for use by jth python program
C
      SUBROUTINE PYTHONMATRIX(A,NR,NC,NDIMR,NDIMC,FILENAME)
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(NDIMR,NDIMC)
      CHARACTER FILENAME*(*)
      WRITE(6,*) FILENAME
      OPEN(100, FILE = FILENAME, STATUS = 'NEW', FORM = 'FORMATTED' )
      WRITE(100,'(I10)') NR
      WRITE(100,'(I10)') NC
      DO J = 1,NC
         DO I = 1,NR
            WRITE(100,'(D22.12)') A(I,J)
         END DO
      END DO
      CLOSE (100)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck genden */
      SUBROUTINE GENDEN(DMAT,CMO,IFLAG,IPRINT)
C*****************************************************************************
C
C     Get (projection) projection density matrices (normed to unity)
C     IFLAG is a 6 bit pattern:
C       0 - 000000: No density matrices
C       1 - 000001: Inactive density matrix
C       2 - 000010: Active density matrix (see bit 6-?)
C       4 - 000100: Secondary electron density matrix
C       8 - 001000: Secondary positron density matrix
C      16 - 010000: C Ct (that is 1 in AO-basis)
C      32 - 100000: reference vacuum (for QED)
C
C     Important: Act. density matrix works only for one open shell.
C
C     Written by J. Thyssen
C     Last revision: 1997/11/23 - jth
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D1=1.0D00,DM1=-1.0D00,D0=0.0D00)
C
      LOGICAL LBIT
      DIMENSION DMAT(*), CMO(*)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
C
#include "ibtfun.h"
C
      CALL QENTER('GENDEN')
C
C
C     Initialize matrix
C
      IF (LBIT(IFLAG,6).AND.VACREF) THEN
        LUBVAC = 14
        CALL OPNFIL(LUBVAC,'DFBVAC','OLD','GENDEN')
        CALL READT(LUBVAC,N2BBASXQ,DMAT)
        CLOSE(LUBVAC,STATUS='KEEP')
        FAC    = DM1
      ELSE
        CALL DZERO(DMAT,N2BBASXQ)
        FAC    = D0
      ENDIF
C
C
      IF(IFLAG.EQ.0) GOTO 999
C
      ITEST = 0
      DO 10 I = 1,5
      IF (LBIT(IFLAG,I)) THEN
        DO 20 IFRP = 1,NFSYM
          IF (I.EQ.1) THEN
            ITEST  = NISH(IFRP)
            ISTART = NPSH(IFRP)+1
            ILOOP  = NISH(IFRP)
          ELSEIF (I.EQ.2) THEN
            ITEST  = NASH(IFRP)
            ISTART = NPSH(IFRP)+NISH(IFRP)+1
            ILOOP  = NASH(IFRP)
          ELSEIF (I.EQ.3) THEN
            ITEST  = NSSH(IFRP)
            ISTART = NPSH(IFRP)+NISH(IFRP)+NASH(IFRP)+1
            ILOOP  = NSSH(IFRP)
          ELSEIF (I.EQ.4) THEN
            ITEST  = NPSH(IFRP)
            ISTART = 1
            ILOOP  = NPSH(IFRP)
          ELSEIF (I.EQ.5) THEN
            ITEST  = NORB(IFRP)
            ISTART = 1
            ILOOP  = NORB(IFRP)
          ENDIF
          IF (ITEST.EQ.0) GOTO 20
CMI ... inserted control print-out of variables
          IF (IPRINT.GE.10) THEN
           WRITE(LUPRI,'(/,8X,A,/)')
     &     '* Output from GENDEN:  before calling DENST1'
           WRITE(LUPRI,*) 'IFRP,I2BASX(IFRP,IFRP):',
     &     IFRP,I2BASX(IFRP,IFRP)
           WRITE(LUPRI,*) 'FAC,NTBAS(0),NFBAS(IFRP,0),NORB(IFRP):',
     &      FAC,NTBAS(0),NFBAS(IFRP,0),NORB(IFRP)
           WRITE(LUPRI,*) 'ISTART,ILOOP:',ISTART,ILOOP
           WRITE(LUPRI,'(/,3X,A)')
     &     'GENDEN - entering molecular orbitals (CMO):'
           WRITE(LUPRI,'(/A,I1,A,I1)')
     &     '* Fermion irrep ',IFRP,'/',NFSYM
           CALL PRQMAT(CMO(1+ICMOQ(IFRP)),NFBAS(IFRP,0),NORB(IFRP),
     &            NFBAS(IFRP,0),NORB(IFRP),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

          CALL DENST1(DMAT(I2BASX(IFRP,IFRP)+1),NTBAS(0),NTBAS(0),NZ,
     &              D1,FAC,CMO(ICMOQ(IFRP)+1),NFBAS(IFRP,0),NORB(IFRP),
     &              ISTART,ILOOP,NFBAS(IFRP,0))
 20    CONTINUE
      IF (ITEST.NE.0) FAC= D1
      ENDIF
 10   CONTINUE
      IF (IPRINT.GE.10) THEN
        WRITE(LUPRI,'(/,8X,A,/)') '* Output from GENDEN:'
        DO I = 1,6
        IF (LBIT(IFLAG,I)) THEN
          IF (I.EQ.1) THEN
            WRITE(LUPRI,'(3X,A)') 
     &           '- inactive density matrix'
          ELSEIF (I.EQ.2) THEN
            WRITE(LUPRI,'(3X,A)') 
     &           '- active density matrix'
          ELSEIF (I.EQ.3) THEN
            WRITE(LUPRI,'(3X,A)') 
     &           '- Sec. (electron) density matrix'
          ELSEIF (I.EQ.4) THEN
            WRITE(LUPRI,'(3X,A)') 
     &           '- Sec. (positron) density matrix'
          ELSEIF (I.EQ.5) THEN
            WRITE(LUPRI,'(3X,A)') 
     &           '- C Ct (1 in AO-basis)'
          ELSEIF (I.EQ.6) THEN
            WRITE(LUPRI,'(3X,A)') '- reference vacuum'
          ENDIF
        ENDIF
        ENDDO
        DO IFRP = 1, NFSYM
          WRITE(LUPRI,'(/A,A3)')
     &      '* Fermion ircop ',FREP(IFRP)
          CALL PRQMAT(DMAT(I2BASX(IFRP,IFRP)+1),
     &              NFBAS(IFRP,0),NFBAS(IFRP,0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDDO
      END IF

C
C   
 999  CONTINUE
      CALL QEXIT('GENDEN')
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck actden */
      SUBROUTINE ACTDEN(DMAT,ALPHA,BETA,CMO,IOPEN,IPRINT)
C*****************************************************************************
C
C     Get active density matrix number IOPEN
C
C     Written by J. Thyssen
C     Last revision: 1997/12/27 - jth
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION DMAT(*), CMO(*)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
C
      CALL QENTER('ACTDEN')
C
      DO 20 ISYM = 1,NFSYM
         ITEST  = NACSH(ISYM,IOPEN)
         ISTART = NPSH(ISYM)+NISH(ISYM)+1
         ILOOP  = NACSH(ISYM,IOPEN)
         DO JOPEN = 1,IOPEN-1
            ISTART = ISTART + NACSH(ISYM,JOPEN)
         END DO
         IF (ITEST.EQ.0) GOTO 20
         CALL DENST1(DMAT(I2BASX(ISYM,ISYM)+1),NTBAS(0),NTBAS(0),NZ,
     &              ALPHA,BETA,
     &              CMO(ICMOQ(ISYM)+1),NFBAS(ISYM,0),NORB(ISYM),
     &              ISTART,ILOOP,NFBAS(ISYM,0))
         IF (IPRINT.GE.10) THEN
            WRITE(LUPRI,'(/3X,A,I2,A,I1,A,I1/)') 
     &         'ACTDEN: Active density matrix number ',IOPEN,
     &         ', corep ',ISYM,'/',NFSYM
            CALL PRQMAT(DMAT(I2BASX(ISYM,ISYM)+1),
     &                NFBAS(ISYM,0),NFBAS(ISYM,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
            CALL FLSHFO(LUPRI)
         END IF
 20   CONTINUE
C
      CALL QEXIT('ACTDEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Denmat */
      SUBROUTINE DENMAT(DMAT,CMO,IPRINT)
C***********************************************************************
C
C     Generate density matrix for SCF calculations both for closed-
C     and open-shell systems.
C
C     Entering CMO - MO coefficients, on output - DMAT, density matrix
C
C     Calls GENDEN, ACTDEN
C
C     Written by T. Saue Nov 23 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,DM1=-1.0D0)
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbdhf.h"
      DIMENSION DMAT(*),CMO(*)
C
      CALL QENTER('DENMAT')

      IF (IPRINT.GE.8) THEN
        CALL HEADER('Output from DENMAT',-1)
        WRITE(LUPRI,'(2X,A,I3)') 'NFMAT=',NFMAT
        WRITE(LUPRI,'(3X,A)')
     &  'Entering molecular orbitals (CMO):'
        DO IFRP = 1, NFSYM
          WRITE(LUPRI,'(/A,I1,A,I1)')
     &    '* Fermion irrep ',IFRP,'/',NFSYM
          CALL PRQMAT(CMO(1+ICMOQ(IFRP)),NFBAS(IFRP,0),NORB(IFRP),
     &            NFBAS(IFRP,0),NORB(IFRP),NZ,IPQTOQ(1,0),LUPRI)
        ENDDO
      ENDIF
      CALL DZERO(DMAT,N2BBASXQ*NFMAT)

      IDEN = 1
      IF(VACPOL.AND..NOT.MIXVAC) THEN
        IF(POSVAC) THEN
          IDEN = IDEN + 4 + 32
        ELSE
          IDEN = IDEN + 8 + 32
        ENDIF
      ENDIF
      CALL GENDEN(DMAT,CMO,IDEN,IPRINT)
CTROND...This only works for zero electrons !!!
      IF(POSVAC) CALL DSCAL(N2BBASXQ,DM1,DMAT,1)
      IF(MIXVAC) CALL GENMIXVAC(DMAT,CMO,IPRINT)
      IF(AOC) THEN
C.......Average of Configuration
        DO IOPEN = 1,NOPEN
          CALL ACTDEN(DMAT(1+N2BBASXQ*IOPEN),D1,D0,CMO,IOPEN,IPRINT)
        END DO
      ELSE
C.......Fractional Occupation
        DO IOPEN = 1,NOPEN
          CALL ACTDEN(DMAT,DF(IOPEN),D1,CMO,IOPEN,IPRINT)
        END DO
        IF (IPRINT.GE.8) THEN
          WRITE(LUPRI,"(2X,A)")
     &    'DENMAT: Total density matrix, open-shell fract. occupation'
          DO IFRP = 1, NFSYM
            WRITE(LUPRI,'(/A,A3)')
     &        '* Fermion ircop ',FREP(IFRP)
            CALL PRQMAT(DMAT(I2BASX(IFRP,IFRP)+1),
     &              NFBAS(IFRP,0),NFBAS(IFRP,0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
          ENDDO
        ENDIF
      ENDIF
C
      CALL QEXIT('DENMAT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck densty */
      SUBROUTINE DENSTY(DMAT,CMO,IPRINT)
C*****************************************************************************
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION DMAT(*),CMO(*)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
CTROND
      CALL DZERO(DMAT,N2BBASXQ)
      DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        CALL DENST1(DMAT(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,D1,D0,
     &              CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &              NPSH(I)+1,NISH(I),NFBAS(I,0))
   10 CONTINUE
C
C     Print section
C
      IF(IPRINT.GE.10) THEN
        CALL TITLER('Output from DENSTY','*',103)
        DO 20 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 20
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(DMAT(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   20   CONTINUE
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck denst1 */
      SUBROUTINE DENST1(DMAT,LRD,LCD,NZ,ALPHA,BETA,
     &                   CMO,LRC,LCC,JVEC,NVEC,NBAS)
C*****************************************************************************
C
C     Construct density matrix for fermion corep IFRP
C
C     DMAT:= alpha*CMO*CMO(H) +  beta*DMAT
C
C     Written by T.Saue Mar 3 1995
C
!------------------------------------------------------------------------
!  N.B. (sya, 2007.05.04)
!     This routine makes AO-basis density matrix assuming the
!     occupation is constant.  For closed-shell DHF, the occupation
!     is 1, which is a half of 2.  This halving is compensated in
!     subroutine PRPEX2.
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION DMAT(LRD,LCD,NZ),CMO(LRC,LCC,NZ),IQ(4)
      DATA IQ/1,2,3,4/
      IF(NVEC.LT.0) THEN
        WRITE(6,*) 'DENST1: Negative number of vectors: ', NVEC
        CALL QUIT('DENST1: Negative number of vectors !')
      ENDIF
      CALL QGEMM(NBAS,NBAS,NVEC,ALPHA,
     &     'N','N',IQ,CMO(1,JVEC,1),LRC,LCC,NZ,
     &     'H','N',IQ,CMO(1,JVEC,1),LRC,LCC,NZ,
     &          BETA,IQ,DMAT,LRD,LCD,NZ)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Density_attachment */
      SUBROUTINE ATTACHMENT_DENSITY(DAO,IOPT,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Create attachment/detachment densities from a difference density
C     The approach is inspired by section 4.4. of 
C     A. Dreuw and M. Head-Gordon:
C      "Single-Reference ab initio Methods for the Calculation of
C       Excited States of Large Molecules"
C       Chem. Rev. 105(2005) 4009-4037
C
C     The starting point is a (difference) density matrix in AO-basis (DAO).
C     It is transformed to MO-basis for numerical stability and 
C     then diagonalized. The diagonal matrix can be split into two parts:
C       1) The detachment density is generated from the eigenvectors 
C          corresponding to negative eigenvalues
C       2) The attachment density is generated from the eigenvectors
C          corresponding to positive eigenvalues
C
C     IOPT = 0: Simple backtransformation to restore the original matrix
C               (debug option)
C     IOPT < 0: Generate detachment density in AO-basis
C     IOPT > 0: Generate attachment density in AO-basis
C
C     The original references are:
C       Martin Head-Gordon, Ana M. Grana, David Maurice, and Christopher A. White:
C       "Analysis of Electronic Transitions as the Difference of Electron 
C        Attachment and Detachment Densities", 
C        J. Phys. Chem, 99 (1995) 14261
C
C       Ana M. Grana, Timothy J. Lee, and Martin Head-Gordon:
C       "Ab Initio Calculations of Singlet and Triplet Excited States of 
C        Chlorine Nitrate and Nitric Acid",
C        J. Phys. Chem, 99 (1995) 3493
C     Written by T. Saue

C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0 = 0.0D0,D2=2.0D0, D1 = 1.0D00)
C
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbdhf.h"
      DIMENSION DAO(*),WORK(*)
C
      KFRSAV = KFREE
      CALL MEMGET('REAL',KTMAT,N2TMT   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDVEC,N2BBASXQ,WORK,KFREE,LFREE)
C.....Read inverse MO-transformation
      CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','PREDHF')
      READ(LUTMAT)
      CALL READT(LUTMAT,N2TMT,WORK(KTMAT))
      CLOSE(LUTMAT)
C
      OCC = D0
      DO I = 1,NFSYM
        IF(IPRINT.GE.6) THEN
          WRITE(6,'(A)') 
     &     'ATTACHMENT_DENSITY: Density matrix in AO-basis'
          CALL PRQMAT(DAO(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                 NTBAS(I),NTBAS(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C       Transform density matrix to MO-basis
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &      NTMO(I),NTMO(I),
     &      DAO(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &      WORK(I2TMOT(I)+KDMAT),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &      WORK(KTMAT+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &      WORK(KTMAT+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &      WORK(KFREE),LFREE,0)
        IF(IPRINT.GE.6) THEN
          WRITE(6,'(A)') 
     &     'ATTACHMENT_DENSITY: Density matrix in MO-basis'
          CALL PRQMAT(WORK(I2TMOT(I)+KDMAT),NTMO(I),NTMO(I),
     &                 NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C       Diagonalize and extract the requested density
        CALL DENSEP(WORK(I2TMOT(I)+KDMAT),NTMO(I),IOPT,OCCBUF,
     &              IPRINT,WORK,KFREE,LFREE)
        OCC = OCC + OCCBUF
      ENDDO
C.....Read MO-transformation
      CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','PREDHF')
      CALL READT(LUTMAT,N2TMT,WORK(KTMAT))
      CLOSE(LUTMAT)
C.....and transform back to AO-basis
      DO I = 1,NFSYM
        CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &      NTMO(I),NTMO(I),
     &      DAO(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &      WORK(I2TMOT(I)+KDMAT),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &      WORK(KTMAT+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &      WORK(KTMAT+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &      WORK(KFREE),LFREE,0)
        IF(IPRINT.GE.5) THEN
          IF(IOPT.EQ.0) THEN
            WRITE(6,'(A)') 
     &      'ATTACHMENT_DENSITY: Density matrix restored to AO-basis'
          ELSEIF(IOPT.LT.0) THEN
            WRITE(6,'(A)') 
     &      'ATTACHMENT_DENSITY: Detachment density in AO-basis'
          ELSE
            WRITE(6,'(A)') 
     &      'ATTACHMENT_DENSITY: Attachment density in AO-basis'
          ENDIF
          CALL PRQMAT(DAO(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(I),NTBAS(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
      ENDDO
C
      WRITE(LUPRI,'(A,F16.9)') 
     &  '* ATTACHMENT_DENSITY: Electron number :',OCC
      CALL MEMREL('ATTACHMENT_DENSITY',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck densep */
      SUBROUTINE DENSEP(DMAT,NDIM,IOPT,OCC,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     The density matrix in MO-basis, DMAT, is diagonalized.
C     
C     IOPT = 0: Simple backtransformation to restore the original matrix
C               (debug option)
C     IOPT < 0: Generate detachment density in AO-basis
C     IOPT > 0: Generate attachment density in AO-basis
C     
C     Written by T. Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,DM1=-1.0D0)
C
#include "dgroup.h"
      DIMENSION DMAT(NDIM,NDIM,NZ),WORK(*)
C
      KFRSAV = KFREE
C
      N2DIM = NDIM*NDIM
      N2DIMQ = N2DIM*NZ
      CALL MEMGET('REAL',KEIG,NDIM ,WORK,KFREE,LFREE)            
      CALL MEMGET('REAL',KMAT,N2DIMQ,WORK,KFREE,LFREE)            
      CALL MEMGET('REAL',KVEC,N2DIMQ,WORK,KFREE,LFREE)            
      CALL MEMGET('REAL',KEMAT,N2DIM,WORK,KFREE,LFREE)            
C
      CALL QDIAG(NZ,NDIM,DMAT,NDIM,NDIM,WORK(KEIG),1,
     &           WORK(KVEC),NDIM,NDIM,WORK(KFREE),LFREE,IERR)
      IF(IPRINT.GE.2) THEN
        WRITE(6,'(A)') '* DENSEP:Eigenvalues...'
        CALL OUTPUT(WORK(KEIG),1,NDIM,1,1,NDIM,1,1,LUPRI)
      ENDIF
C
      CALL EIGSEP(WORK(KEIG),WORK(KEMAT),IOPT,OCC,NDIM)
C     Using selected eigenvectors, restore requested density in MO-basis
      CALL QTRANS('MOAO','S',D0,NDIM,NDIM,NDIM,NDIM,
     &      DMAT,NDIM,NDIM,NZ,IPQTOQ(1,0),
     &      WORK(KEMAT),NDIM,NDIM,1,IPQTOQ(1,0),
     &      WORK(KVEC),NDIM,NDIM,NZ,IPQTOQ(1,0),
     &      WORK(KVEC),NDIM,NDIM,NZ,IPQTOQ(1,0),
     &      WORK(KFREE),LFREE,0)
C
      IF(IPRINT.GE.6) THEN
        IF(IOPT.EQ.0) THEN
          WRITE(6,'(A)') 
     &    'DENSEP: Density matrix restored to MO-basis'
        ELSEIF(IOPT.LT.0) THEN
          WRITE(6,'(A)') 
     &    'DENSEP: Detachment density in MO-basis'
        ELSE
          WRITE(6,'(A)') 
     &    'DENSEP: Attachment density in MO-basis'
        ENDIF
        CALL PRQMAT(DMAT,NDIM,NDIM,NDIM,NDIM,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      CALL MEMREL('DENSEP',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Eigsep */
      SUBROUTINE EIGSEP(EIG,EMAT,IOPT,OCC,NDIM)
C***********************************************************************
C
C     IOPT = 0: Construct diagonal difference density matrix
C     IOPT < 0: Construct diagonal detachment density matrix
C     IOPT > 0: Construct diagonal attachment density matrix
C
C     This is not a memory-optimal way of doing things....
C
C     Written by T. Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION EIG(NDIM),EMAT(NDIM,NDIM)
C
      N2DIM = NDIM*NDIM
      CALL DZERO(EMAT,N2DIM)
C    
      OCC=D0
      IF(IOPT.EQ.0) THEN
C.......difference density (debug)
        DO I = 1,NDIM
          EMAT(I,I) = EIG(I)
          OCC = OCC + EMAT(I,I)
        ENDDO   
      ELSEIF(IOPT.LT.0) THEN
C.......detachment density
        DO I = 1,NDIM
          EMAT(I,I) = MIN(D0,EIG(I))
          OCC = OCC + EMAT(I,I)
        ENDDO   
      ELSE
C.......attachment density
        DO I = 1,NDIM
          EMAT(I,I) = MAX(D0,EIG(I))
          OCC = OCC + EMAT(I,I)
        ENDDO   
C
      ENDIF
      WRITE(LUPRI,'(A,F16.8)') 
     &   'EIGSEP: * Electron number: ',OCC
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck genaltvac */
      SUBROUTINE GENMIXVAC(DMAT,CMO,IPRINT)
C*****************************************************************************
C
C     Calculate alternative definition of density matrix of 
C     vacuum polarisation:
C
C     The standard definition is:
C
C     D(vacpol) = [D(-) - D(-)_ref]
C
C     Alternatively:
C
C     D(vacpol) = (1/2)[D(-) - D(+)]
C    
C
C     Note that the density matrix is NOT initialized to zero !
C
C     Written by  T. Saue Oct 9 2006
C     Last revision: 1997/11/23 - jth
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D1=1.0D00,DM1=-1.0D00,D0=0.0D00,DP5=0.5D0,DPM5=-0.5D0)
C
      DIMENSION DMAT(*), CMO(*)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
#include "ibtfun.h"
C
      CALL QENTER('GENMIXVAC')
C
C     Calculate negative energy part
C
      DO IFRP = 1,NFSYM
        ITEST  = NPSH(IFRP)
        ISTART = 1
        ILOOP  = NPSH(IFRP)
        IF(ITEST.GT.0) THEN
          CALL DENST1(DMAT(I2BASX(IFRP,IFRP)+1),NTBAS(0),NTBAS(0),NZ,
     &              DP5,D1,CMO(ICMOQ(IFRP)+1),NFBAS(IFRP,0),NORB(IFRP),
     &              ISTART,ILOOP,NFBAS(IFRP,0))
        ENDIF
        ITEST  = NESH(IFRP)
        ISTART = NPSH(IFRP)+1
        ILOOP  = NESH(IFRP)
        IF(ITEST.GT.0) THEN
          CALL DENST1(DMAT(I2BASX(IFRP,IFRP)+1),NTBAS(0),NTBAS(0),NZ,
     &              DPM5,D1,CMO(ICMOQ(IFRP)+1),NFBAS(IFRP,0),NORB(IFRP),
     &              ISTART,ILOOP,NFBAS(IFRP,0))
        ENDIF
      ENDDO
      IF (IPRINT.GE.10) THEN
        WRITE(LUPRI,'(A)') '* Output from GENMIXVAC:'
        DO IFRP = 1,NFSYM
          WRITE(LUPRI,'(/A,A3)')
     &      '* Fermion ircop ',FREP(IFRP)
          CALL PRQMAT(DMAT(I2BASX(IFRP,IFRP)+1),
     &              NFBAS(IFRP,0),NFBAS(IFRP,0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDDO
      END IF
C   
 999  CONTINUE
      CALL QEXIT('GENMIXVAC')
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck modens */
      SUBROUTINE MODENS(DMAT,OCC,IPRINT)
C***********************************************************************
C
C     Construct density matrix in MO-basis (first version)
C
C     Written by Trond Saue Aug 11 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION DMAT(NORBT,NORBT,NZ),OCC(NESHT)
      CALL DZERO(DMAT,N2ORBXQ)
      IOFF = 0
      JOFF = 0
      DO IFRP = 1,NFSYM
        IOFF = IOFF + NPSH(IFRP)
        IF(NISH(IFRP) .GT. 0) THEN
          DO J = 1,NESH(IFRP)
            DMAT(IOFF+J,IOFF+J,1) = OCC(JOFF+J)
          ENDDO
        ENDIF
        IOFF = IOFF + NESH(IFRP)
        JOFF = JOFF + NESH(IFRP)
      ENDDO
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(A)') '* MO density matrix:'
        CALL OUTPUT(DMAT,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
      ENDIF
      RETURN
      END
