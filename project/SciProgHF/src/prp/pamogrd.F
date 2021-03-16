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

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rorbgrd */
      SUBROUTINE RORBGRD(LORBTRI,DV,FC,FCX,FVX,FQX,GORB,JXOP,NXOP,
     &                   FPV,DFAC,IPQ,IPRINT)
C***********************************************************************
C
C     Calculate:
C     [ X_{pq}, H ]         (LORBTRI = .TRUE. ),
C     [ X_{pq}, \tilde{H} ] (LORBTRI = .TRUE. ), or
C     [ X_{pq}, \hat{H} ]   (LORBTRI = .FALSE.)
C
C     = DFAC(1) * F_{pq} + DFAC(2) * F_{qp}
C
C     Input:
C        DV      - one-electron active density matrix
C        FC      - the inactive Fock matrix for DV * FC term
C        FCX     - the inactive Fock matrix
C        FVX     - the active Fock matrix
C        FQX     - the Q Fock matrix
C        JXOP    - orbital rotations
C        FPV     - scratch memory used for F_{pv}
C        DFAC    - factors used in grad. expression
C        IPQ     - index array
C        LORBTRI - true for an orbital trial vector
C
C     Output:
C        GORB    - the gradient
C     
C     Written by J. Thyssen - Oct 22 1998
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
      PARAMETER (D2 = 2.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "maxorb.h"
#include "dcbidx.h"
C
      DIMENSION DV(NASHT,NASHT,NZ),  FC(NORBT,NORBT,NZ)
      DIMENSION FCX(NORBT,NORBT,NZ), FVX(NORBT,NORBT,NZ)
      DIMENSION FQX(NORBT,NASHT,NZ), FPV(NORBT,NASHT,NZ)
      DIMENSION GORB(NXOP,NZ), JXOP(2,*), DFAC(2), IPQ(*)
      LOGICAL   LORBTRI
C
      LOGICAL   NOFC, NOFV
C
      CALL QENTER('RORBGRD')
C
C
      NOFV = (NASHT .EQ. 0)
      NOFC = .NOT. LORBTRI
C
C
C     ***************************************************
C     *** Calculate F_{pv}                          *****
C     ***   = \sum_{u} D^V_{vu} F^C_{pu} + F^Q_{pv} *****
C     ***************************************************
C
C
      IF (.NOT.NOFV) CALL FGEN_PV(FC,FQX,DV,FPV,IPQ,IPRINT)
C
 9000 FORMAT(' Orbital trial vector: ',L1)
 9001 FORMAT(1X,A,' matrix'/)
 9002 FORMAT(' F_{',A,'} factor: ',F4.1)
 9003 FORMAT(' Fock matrix ',A,' not included in grad. calc.: ',L1)
C
      IF ( IPRINT .GE. 20 ) THEN
         CALL HEADER('Output from RORBGRD',-1)
         WRITE(LUPRI,9000) LORBTRI
         WRITE(LUPRI,9003) 'NOFC',NOFC
         WRITE(LUPRI,9003) 'NOFV',NOFV
         WRITE(LUPRI,9002) 'lk',DFAC(1)
         WRITE(LUPRI,9002) 'kl',DFAC(2)
C
         IF (.NOT.NOFC) THEN
            WRITE(LUPRI,9001) 'FCX'
            CALL PRQMAT(FCX,NORBT,NORBT,NORBT,NORBT,NZ,IPQ,LUPRI)
         END IF
C
         IF (.NOT.NOFV) THEN
            WRITE(LUPRI,9001) 'FVX'
            CALL PRQMAT(FVX,NORBT,NORBT,NORBT,NORBT,NZ,IPQ,LUPRI)
C        
            WRITE(LUPRI,9001) 'DV'
            CALL PRQMAT(DV,NASHT,NASHT,NASHT,NASHT,NZ,IPQ,LUPRI)
C
            WRITE(LUPRI,9001) 'FC'
            CALL PRQMAT(FC,NORBT,NORBT,NORBT,NORBT,NZ,IPQ,LUPRI)
C
            IF (NASHT .GT. 1) THEN
               WRITE(LUPRI,9001) 'FQX'
               CALL PRQMAT(FQX,NORBT,NASHT,NORBT,NASHT,NZ,IPQ,LUPRI)
            END IF
C
            WRITE(LUPRI,9001) 'F_{pv}'
            CALL PRQMAT(FPV,NORBT,NASHT,NORBT,NASHT,NZ,IPQ,LUPRI)
         END IF
C
      END IF
C     
C
C     ***************************************************
C     *** Calculate gradient elements *******************
C     ***************************************************
C
C
      CALL DZERO(GORB,NXOP * NZ)
C
      DO IG = 1,NXOP
         K = JXOP(1,IG)
         L = JXOP(2,IG)
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         ISYMK = IFSMO(K)
         ISYML = IFSMO(L)
C        
C        First index (K)
C        ---------------
C
C        g_{ip} = 2 ( F_{pi} - F_{ip}^* )
C        
C        g_{vp} = 2 ( F_{pv} - F_{vp}^* )
C        
         IF (ITYPK .EQ. JTINAC) THEN
C
C           First index is inactive
C           -----------------------
C
C           F_{pi} = 2 * ( FC_{pi} + FV_{pi} )
C
            IF (.NOT. NOFC) THEN
               DO IZ = 1,NZ 
                  GORB(IG,IZ) = DFAC(1) * D2 * FCX(L,K,IZ)
               END DO
            END IF
            IF (.NOT. NOFV) THEN
               DO IZ = 1,NZ 
                  GORB(IG,IZ) = GORB(IG,IZ) +
     &                 DFAC(1) * D2 * FVX(L,K,IZ)
               END DO
            END IF
         ELSE IF (ITYPK .EQ. JTACT) THEN
C
C           First index is active
C           ---------------------
C
C           add the precalculated F_{pv} (F_{LK})
C
            KU = IDXG2U(K)
            DO IZ = 1,NZ 
               GORB(IG,IZ) = DFAC(1) * FPV(L,KU,IZ)
            END DO
C        ELSE IF (ITYPK .EQ. JTSEC) THEN
C           
C           First index is secondary or positronic
C           --------------------------------------
C
C           Add 0
C
         END IF
C
C            
C        Second index (L)
C        ----------------
C
         IF (ITYPL .EQ. JTACT) THEN
C
C           Second index is active
C           ----------------------
C
C           Add F_{vp}^* (F_{KL})
C
            LU = IDXG2U(L)
            DO IZ = 1,NZ
               IP = IPQ(IZ)
               IF (IP .EQ. 1) THEN
                  GORB(IG,IZ) = GORB(IG,IZ) + DFAC(2) * FPV(K,LU,IZ)
               ELSE
                  GORB(IG,IZ) = GORB(IG,IZ) - DFAC(2) * FPV(K,LU,IZ)
               END IF
            END DO
C        ELSE
C
C           Second index is secondary
C           -------------------------
C           Add 0
C
         END IF

      END DO
C
C
  200 CALL QEXIT('RORBGRD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fgen_pv */
      SUBROUTINE FGEN_PV(FC,FQ,DV,FGEN,IPQ,IPRINT)
C***********************************************************************
C
C     Calculate (4.57:ii) of jth thesis
C               \sum_u F^C_{pu} D_{uv}^* + F^Q_{pv}
C               = (F^C * D^t + F^Q)_{pv}
C
C     Input:
C        FC      - inactive Fock matrix
C        DV      - active density matrix
C        FQ      - the Q Fock matrix
C
C     Output:
C        FGEN    - generalized Fock matrix with second index active
C
C     Written by J. Thyssen - Nov 16 1998
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (D1 = 1.00D00, D0 = 0.00D00 )
C
#include "dgroup.h"
#include "dcborb.h"
      DIMENSION FC(NORBT,NORBT,NZ),FQ(NORBT,NASHT,NZ),
     &          FGEN(NORBT,NASHT,NZ),IPQ(*),DV(NASHT,NASHT,NZ)
C
      CALL QENTER('FGEN_PV')
c      IF (IPRINT .GE. 5) CALL HEADER('Output from FGEN_PV',-1)
C
C     *****************************
C     *** F^Q_{pv} contribution ***
C     *****************************
C
C
      IF (NASHT .LE. 1) THEN
         FAC = D0
      ELSE
         FAC = D1
         CALL DCOPY(NORBT * NASHT * NZ, FQ, 1, FGEN, 1)
      END IF
C
C     ****************************
C     *** DV * FC contribution ***
C     ****************************
C
C
      DO I = 1,NFSYM
C
         IF (NASH(I) .GT. 0) THEN
            CALL QGEMM(NORBT,NASH(I),NASH(I),D1,
     &           'N','N',IPQTOQ(1,0),
     &           FC(1,1+IORB(I)+NPSH(I)+NISH(I),1),
     &           NORBT,NORBT,NZ,
     &           'T','N',IPQTOQ(1,0),DV(1+IASH(I),1+IASH(I),1),
     &           NASHT,NASHT,NZ,
     &           FAC,IPQTOQ(1,0),
     &           FGEN(1,1+IASH(I),1),NORBT,NASH(I),NZ)
         END IF
c         IF (IPRINT .GE. 5) THEN
c            WRITE(LUPRI,9000) I
c9000 FORMAT(/'Generalized Fock matrix irrep ',I1,
c    &        ' with second index active')
c            CALL PRQMAT(FGEN(1,1+IASH(I),1),NORBT,NASH(I),NORBT,NASHT,
c     &           NZ,IPQTOQ(1,0),LUPRI)
c         END IF
      END DO
C
      CALL QEXIT('FGEN_PV')
      RETURN
      END
