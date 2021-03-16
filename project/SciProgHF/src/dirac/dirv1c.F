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

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type1drv */
      SUBROUTINE TYPE1DRV(OP1INT,WORK,LWORK,DOINT,OMITVNUC,QMAT,IPRINT) 
C*****************************************************************************
C
C    TYPE1DRV : Driver for coulomb attraction integrals,
C               used to mimic the multicenter integrals
C               neglected in the one-center model 1 (Luuk's model)
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
      LOGICAL DOINT(2,2)
      CHARACTER*4 OMITVNUC
      DIMENSION OP1INT(NNBBASX),WORK(LWORK),QMAT(NUCDEP,NUCDEP,2)
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
C
      CALL QENTER('TYPE1DRV')
c
c     Start timer
c
      TIMSTR = SECOND()
#include "memint.h"
C
      CALL MEMGET('REAL',KVCMAT,  NNBBASX,WORK,KFREE,LFREE)
      LFCBA = NUCDEP**3
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXGEXP,  NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
      CALL DZERO(WORK(KVCMAT),NNBBASX)
c   
c     ... set charges and gaussian exponents.
c
      CALL DZERO(WORK(KCHRG),NUCIND)
      IOFF = 1
      DO I = 1,NUCIND
         WORK(KXGEXP-1+I) = GNUEXP(I)
         IF (ICTLV1C(2) .EQ. 2) THEN
c          ... Simulate terms neglected in Luuk's original model
c              (simulate the projection scheme)
            IF (IOFF .GT. NUCDEP)
     &          CALL QUIT('TYPE1DRV : IOFF .GT. NUCDEP')
c          ... Q    =   - q^L - q^S
            DO J = 1,NUCDEP
               WORK(KCHRG-1+I) = WORK(KCHRG-1+I) + QMAT(J,IOFF,1) 
     &                                           + QMAT(J,IOFF,2)
            END DO
            IOFF = IOFF + FMULT(ISTBNU(I))
         ELSE
            WORK(KCHRG-1+I) = -CHARGE(I)
         END IF
      ENDDO
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from TYPE1DRV',-1)
        WRITE(LUPRI,'(/1X,A)')
     &       '* Charges used for correction  :'
              WRITE(LUPRI,'(4X,A4,F12.8)') 
     &       (NAMN(I),WORK(KCHRG-1+I),I = 1,NUCIND)
      ENDIF
c
c     ... Make list of nuclei.
c
      CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &            WORK(KXGEXP),WORK(KCOORC), WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),OMITVNUC,IPRINT)
c
c     ... Calculate integrals
c
      CALL CATDR1(WORK(KVCMAT),DUMMY,NCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LWRK,IPRINT,
     &            .FALSE.,0,NNBBASX,DOINT)
C
C     Symmetry unpack integrals
C
      CALL SYMUPK(WORK(KVCMAT),OP1INT,1,NNBBASX)
C
C     Memory deallocation
C
      CALL MEMREL('TYPE1DRV.CATDR1',WORK,KWORK,KWORK,KFREE,LFREE)
c
c     Stop timer
c
      TIMEND = SECOND() 
      TIME   = TIMEND - TIMSTR
      WRITE (LUPRI,'()')
      CALL TIMTXT('>>> Time used in TYPE1DRV is',TIME,LUPRI)
      CALL QEXIT('TYPE1DRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type2drv */
      SUBROUTINE TYPE2DRV(OP1INT,WORK,LWORK,DOINT,OMITVNUC,
     &                     QMAT,IPRINT) 
C*****************************************************************************
C
C    TYPE2DRV : Driver for coulomb attraction integrals,
C               used to mimic the multicenter integrals
C               neglected in the one-center model 2.
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
      LOGICAL DOINT(2,2)
      CHARACTER*4 OMITVNUC
      DIMENSION OP1INT(NNBBASX),WORK(LWORK),QMAT(NUCDEP,NUCDEP,2)
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
C
      CALL QENTER('TYPE2DRV')
c
c     Start timer
c
      TIMSTR = SECOND()
c
      IF(DOINT(1,2) .OR. DOINT(2,1) ) THEN
         CALL QUIT ('TYPE2DRV called with LS or SL type integrals.')
      END IF
#include "memint.h"
C
      CALL MEMGET('REAL',KVCMAT,  NNBBASX,WORK,KFREE,LFREE)
      LFCBA = NUCDEP**3
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
      CALL DZERO(WORK(KVCMAT),NNBBASX)
c   
c     ... set charges. 
c         (We use the nuclear exponents)
c
      CALL DZERO(WORK(KCHRG),NUCIND)
      IF(DOINT(1,1) .OR. (DOINT(2,2).AND.OMITVNUC.EQ.'FTTT')) THEN
C            ... First case : Estimate (LL|SS) contributions
C                Last case  : Simulate one-center SS integrals since explicit
C                             evaluation has been omitted in this iteration.
         IF (ICTLV1C(2) .EQ. 1) THEN 
            CALL SCDENS_2(WORK(KCHRG),1,NUCIND,0)
         ELSE
            IOFF = 1
            DO I = 1,NUCIND
               IF (IOFF .GT. NUCDEP) 
     &             CALL QUIT('TYPE2DRV : IOFF .GT. NUCDEP')
               DO J = 1,NUCDEP
                  WORK(KCHRG-1+I) = WORK(KCHRG-1+I) + QMAT(J,IOFF,2)
               END DO
               IOFF = IOFF + FMULT(ISTBNU(I))
            ENDDO
         END IF
      ELSE IF(DOINT(2,2).AND.OMITVNUC.EQ.'TTFF') THEN
c        Estimate the error of this model <A | q_B^L + q_B^S - Z_B | C>
         IOFF = 1
         DO I = 1,NUCIND
            IF (IOFF .GT. NUCDEP) 
     &          CALL QUIT('TYPE2DRV : IOFF .GT. NUCDEP')
            DO J = 1,NUCDEP
               WORK(KCHRG-1+I) = WORK(KCHRG-1+I) 
     &                           + QMAT(J,IOFF,1) + QMAT(J,IOFF,2)
            END DO
            IOFF = IOFF + FMULT(ISTBNU(I))
            WORK(KCHRG-1+I) = WORK(KCHRG-1+I) + CHARGE(I)
         ENDDO
      ELSE 
C        ... Calculate terms to cancel SC nuclear attraction terms (previouly calculated)
C
         DO I = 1,NUCIND
C            ... Test adding penalty function :
C                <S_A | (-X/r_C) | S_B>
             IF (ICTLV1C(2).EQ.3.AND.OMITVNUC.EQ.'TTFF') THEN
                WORK(KCHRG-1+I) = -5.0D0
             ELSE
                WORK(KCHRG-1+I) = -CHARGE(I)
             ENDIF
         ENDDO
      END IF
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from TYPE2DRV',-1)
        IF (OMITVNUC(1:1).EQ.'F') THEN
           WRITE (LUPRI,'(/1X,A/)') 
     &    'INFO : One-center integrals estimated in this call'
        END IF
        WRITE(LUPRI,'(/1X,A)')
     &       '* Charges used for correction  :'
              WRITE(LUPRI,'(4X,A4,F12.8)') 
     &       (NAMN(I),WORK(KCHRG-1+I),I = 1,NUCIND)
      ENDIF
c
c     ... Make list of nuclei.
c
      CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &            GNUEXP,WORK(KCOORC), WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),OMITVNUC,IPRINT)
c
c    ... Calculate integrals
c
      CALL CATDR1(WORK(KVCMAT),DUMMY,NCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LWRK,IPRINT,
     &            .FALSE.,0,NNBBASX,DOINT)
C
C     Symmetry unpack integrals
C
      CALL SYMUPK(WORK(KVCMAT),OP1INT,1,NNBBASX)
C
C     Memory deallocation
C
      CALL MEMREL('TYPE2DRV.CATDR1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
c     Stop timer
c
      TIMEND = SECOND()
      TIME   = TIMEND - TIMSTR
      WRITE (LUPRI,'()')
      CALL TIMTXT('>>> Time used in TYPE2DRV is',TIME,LUPRI)
      CALL QEXIT('TYPE2DRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type3drv */
      SUBROUTINE TYPE3DRV(OP1INT,WORK,LWORK,DOINT,ADDQIC,
     &                    OMITVNUC,IPRINT,QMAT)
C*****************************************************************************
C
C    TYPE3DRV : Driver for coulomb attraction integrals,
C               used to mimic the multicenter integrals
C               neglected in the one-center model 3.
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
C      PARAMETER (D0 = 0.0D0, QTHRS = 1.0D-10)
      PARAMETER (D0 = 0.0D0, QTHRS = 1.0D-8)
      LOGICAL DOINT(2,2), ADDQIC(2)
      CHARACTER*4 OMITVNUC
      DIMENSION OP1INT(NNBBASX),WORK(LWORK),QMAT(NUCDEP,NUCDEP,2,*)
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
C
      CALL QENTER('TYPE3DRV')
c
c     Start timer
c
      TIMSTR = SECOND()
#include "memint.h"
C
      CALL MEMGET('REAL',KVCMAT,  NNBBASX,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KVCMAT),NNBBASX)
C
C     Initialization of charges and coordinates
C
      NNDEP = NUCDEP*(NUCDEP+1)/2
      IF(ICTLV1C(2).EQ. 2) THEN
C
C     ... Atomic centered mulliken charges
C
         MXCENTC = NUCDEP
         KCENTC  = NUCDEP
      ELSE 
C
C     ... Distributed charges
C
         CALL TYPE3MXC(MXCENTC,NNDEP,QMAT,ADDQIC,QTHRS)
         KCENTC  = 0
      END IF
C      IF (MXCENTC.LE.0) CALL QUIT('TYPE3DRV: MXCENTC.le.0' //
C     &   ' (Probably no charges are .GT. Threshold)')
      LFCBA = MXCENTC*NUCDEP**2
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  MXCENTC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*MXCENTC,WORK,KFREE,LFREE)
C     ... these are not used now, should be removed.
C         TODO: (write new AVECA1 etc. instead of AVENA1 etc.)
      CALL MEMGET('REAL',KSIGNC,3*KCENTC,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  KCENTC,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  KCENTC,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  KCENTC,WORK,KFREE,LFREE)
C
      IF(ICTLV1C(2).EQ. 2) THEN
C         ... Get charges and exponents
         CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KXGEXP,  NUCIND,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KCHRG),NUCIND)
         IOFF = 1
         DO I = 1,NUCIND
            IF (IOFF .GT. NUCDEP) 
     &          CALL QUIT('TYPE1DRV : IOFF .GT. NUCDEP')
            DO J = 1,NUCDEP
               IF (ADDQIC(1))
     &             WORK(KCHRG-1+I) = WORK(KCHRG-1+I) + QMAT(J,IOFF,1,1)
               IF (ADDQIC(2))
     &             WORK(KCHRG-1+I) = WORK(KCHRG-1+I) + QMAT(J,IOFF,2,1)
            END DO
            IOFF = IOFF + FMULT(ISTBNU(I))
            WORK(KXGEXP-1+I) = GNUEXP(I)
Ctest        ... finite size exponents for q_L (for Iodine,Carbon,Hydrogen)
C            IF(DOINT(2,2)) THEN
C               IF (NINT(CHARGE(I)) .EQ. 53)   WORK(KXGEXP-1+I) = 0.2D0
C               IF (NINT(CHARGE(I)) .EQ.  6)   WORK(KXGEXP-1+I) = 1.0D0
C               IF (NINT(CHARGE(I)) .EQ.  1)   WORK(KXGEXP-1+I) = 1.5D0
C            END IF
         ENDDO
         IF(IPRINT.GE.3) THEN
            CALL HEADER('Output from TYPE1DRV',-1)
            WRITE(LUPRI,'(/1X,A)')
     &         '* Charges used for correction  :'
                WRITE(LUPRI,'(4X,A4,F12.8)')
     &        (NAMN(I),WORK(KCHRG-1+I),I = 1,NUCIND)
         ENDIF
         CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &               WORK(KXGEXP),WORK(KCOORC), WORK(KSIGNC),
     &               WORK(KJSYMC),WORK(KJCENT),OMITVNUC,IPRINT)
         MXCENTC = NCENTC
      ELSE
         CALL MEMGET('REAL',KCOORA,3*NUCDEP,WORK,KFREE,LFREE)
         CALL GETACORD(WORK(KCOORA))
         CALL TYPE3DR1(MXCENTC,WORK(KCOORC),WORK(KCOORA),WORK(KGEXP),
     &                QMAT,ADDQIC,QTHRS,WORK(KFCBA),OMITVNUC,
     &                NNDEP,IPRINT)
      END IF
C
C    ... Calculate integrals
C
      CALL CATDR1(WORK(KVCMAT),DUMMY,MXCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LFREE,IPRINT,
     &            .FALSE.,0,NNBBASX,DOINT)
C
C     Symmetry unpack nuclear attraction integrals
C
      CALL SYMUPK(WORK(KVCMAT),OP1INT,1,NNBBASX)
C
C     Memory deallocation
C
      CALL MEMREL('TYPE3DRV.CATDR1',WORK,KWORK,KWORK,KFREE,LFREE)
C
c     Stop timer
c
      TIMEND = SECOND()
      TIME   = TIMEND - TIMSTR
      WRITE (LUPRI,'()')
      CALL TIMTXT('>>> Time used in TYPE3DRV is',TIME,LUPRI)
      CALL QEXIT('TYPE3DRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type3mxc */
      SUBROUTINE TYPE3MXC(MXCENTC,NNDEP,QMAT,ADDQIC,QTHRS)
C*****************************************************************************
C
C    TYPE3MXC : Find MXCENTC, the number of sites needed to
C               calculate one-center model 3
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
      PARAMETER (D0 = 0.0D0 ,  D2 = 2.0D0)
      LOGICAL ADDQIC(2)
      DIMENSION QMAT(NUCDEP,NUCDEP,2)
      CALL QENTER('TYPE3MXC')
C
      MXCENTC = 0
      IJ = 0
      DO J = 1, NUCDEP
         DO I = 1,J
            IJ = IJ + 1
            TMP = D0
            IF(ADDQIC(1)) TMP = QMAT(I,J,1) + TMP
            IF(ADDQIC(2)) TMP = QMAT(I,J,2) + TMP
            IF(I.NE.J) TMP = D2*TMP
            IF (ABS(TMP) .GT. QTHRS) MXCENTC = MXCENTC + 1
         ENDDO
      ENDDO
      IF (MXCENTC.EQ.0) THEN
         CALL HEADER('Output from TYPE3MXC',-1)
         WRITE(LUPRI,'(/3X,A,D10.3/)') 
     &   'No charges greater than threshold : ',QTHRS
      ENDIF
      CALL QEXIT('TYPE3MXC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type3dr1 */
      SUBROUTINE TYPE3DR1(MXCENTC,COORC,COORA,XGNUEXP,QMAT,
     &                   ADDQIC,QTHRS,FACCBA,OMITVNUC,NNDEP,IPRINT)
C*****************************************************************************
C
C    TYPE3DR1 : Driver for coulomb attraction integrals,
C               used to mimic the multicenter integrals
C               neglected in the one-center model 3
C               Centers with charges smaller than QTHRS
C               are neglected.
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dcbham.h"
      PARAMETER (D0 = 0.0D0 , D2 = 2.0D0 , DP5 = 0.5D0)
      DIMENSION FACCBA(MXCENTC,NUCDEP,NUCDEP),QMAT(NUCDEP,NUCDEP,2,*),
     &          XGNUEXP(MXCENTC),COORC(3,MXCENTC),COORA(3,NUCDEP)
      LOGICAL NO_CCC, NO_CAA,NO_CCA,NO_CBA, ADDQIC(2)
      CHARACTER*4 OMITVNUC
      CALL QENTER('TYPE3DR1')
C     <A | C | B> = CAB
      NO_CCC = (OMITVNUC(1:1) .EQ. 'T')
      NO_CAA = (OMITVNUC(2:2) .EQ. 'T')
      NO_CCA = (OMITVNUC(3:3) .EQ. 'T')
      NO_CBA = (OMITVNUC(4:4) .EQ. 'T')
      DO I = 1,MXCENTC
C     Initialize to point charge distributions
C     (change later if disired )
         XGNUEXP(I) = D0
      ENDDO
C
      IF(IPRINT.GE.2) THEN
         CALL HEADER('Output from TYPE3DR1',-1)
         WRITE(LUPRI,'(/A)') 
     &        '  ** Distribution of charges  **'
         IF(ICTLV1C(2).EQ.0) THEN
            WRITE(LUPRI,'(A/)') 
     &      '     (Charges placed where dipole moment is zero)'
         ELSE IF(ICTLV1C(2).EQ.3) THEN
            WRITE(LUPRI,'(A/)')
     &      '     (Charges placed at midpoints between nuclei)'
         ENDIF
         WRITE(LUPRI,'(/,4X,A,10X,A,17X,A,24X,A,24X,A,/)') 
     &        'NCENTC','-CHARGE(NCENTC)','X','Y','Z'
      END IF
C
      N2DEP = NUCDEP*NUCDEP
      QSKIP = D0
      NCENTC = 0
      IJ = 0
      DO J = 1, NUCDEP
C
C        off-diagonal terms
C
         DO I = 1,J-1
            IJ = IJ + 1
            QTMP = D0
            IF(ADDQIC(1)) QTMP = QTMP - QMAT(I,J,1,1)
            IF(ADDQIC(2)) QTMP = QTMP - QMAT(I,J,2,1)
            QTMP = D2*QTMP
            IF (ABS(QTMP) .GT. QTHRS) THEN
               NCENTC = NCENTC + 1
               IF (ICTLV1C(2).EQ.0) THEN
                  DO K = 1,3
                     CTMP = D0
                     IF(ADDQIC(1)) CTMP = CTMP - QMAT(I,J,1,1+K)
                     IF(ADDQIC(2)) CTMP = CTMP - QMAT(I,J,2,1+K)
                     COORC(K,NCENTC) = D2 * CTMP / QTMP
                  END DO
C                  IF (ADDQIC(2)) XGNUEXP(NCENTC) = 100.0D0
C                  IF (ADDQIC(1)) XGNUEXP(NCENTC) = 2.0D0
               ELSE IF (ICTLV1C(2).EQ.3) THEN
                  COORC(1,NCENTC) = DP5*(COORA(1,I) + COORA(1,J))
                  COORC(2,NCENTC) = DP5*(COORA(2,I) + COORA(2,J))
                  COORC(3,NCENTC) = DP5*(COORA(3,I) + COORA(3,J))
               END IF
               IF (IPRINT.GE.2) WRITE(LUPRI,'(4X,I3,5X,1P,4D25.12)')
     &         IJ,QTMP,COORC(1,NCENTC),COORC(2,NCENTC),COORC(3,NCENTC)
C
               CALL DCOPY(N2DEP,QTMP,0, FACCBA(NCENTC,1,1),MXCENTC)
C
C              Zero elements of FACCBA which should not be included
C              -- all off-diagonal included in this version
C
               IF (IPRINT.GT.5) THEN
                   WRITE(LUPRI,'(/4X,A)') 'FACCBA for this center:'
                   CALL OUTPUT(FACCBA(NCENTC,1,1),1,1,1,N2DEP,
     &                 MXCENTC,N2DEP,-1,LUPRI)
                   WRITE(LUPRI,'()')
               END IF
            ELSE
              IF (IPRINT.GE.2) WRITE(LUPRI,'(4X,I3,1P,D30.12,A)')
     &              IJ,QTMP, '         skipped!'
              QSKIP = QSKIP + ABS(QTMP)
            END IF
C
         END DO
C
C        diagonal terms
C
         IJ = IJ + 1
         QTMP = D0
         IF(ADDQIC(1)) QTMP = QTMP - QMAT(J,J,1,1)
         IF(ADDQIC(2)) QTMP = QTMP - QMAT(J,J,2,1)
         IF (ABS(QTMP) .GT. QTHRS) THEN
            NCENTC = NCENTC + 1
            IF (ICTLV1C(2).EQ.0) THEN
               DO K = 1,3
                  CTMP = D0
                  IF(ADDQIC(1)) CTMP = CTMP - QMAT(J,J,1,1+K)
                  IF(ADDQIC(2)) CTMP = CTMP - QMAT(J,J,2,1+K)
                  COORC(K,NCENTC) = CTMP / QTMP
               END DO
C               IF (ADDQIC(2)) XGNUEXP(NCENTC) = 100.0D0
C               IF (ADDQIC(1)) XGNUEXP(NCENTC) = 2.0D0
            ELSE IF (ICTLV1C(2).EQ.3) THEN
               COORC(1,NCENTC) = COORA(1,J)
               COORC(2,NCENTC) = COORA(2,J)
               COORC(3,NCENTC) = COORA(3,J)
            END IF
            IF (IPRINT.GE.2) WRITE(LUPRI,'(4X,I3,5X,4D25.12)')
     &          IJ,QTMP,COORC(1,NCENTC),COORC(2,NCENTC),COORC(3,NCENTC)
C
            CALL DCOPY(N2DEP,QTMP,0, FACCBA(NCENTC,1,1),MXCENTC)
C
C        Zero elements which should not be included
C
            IF (NO_CCC) FACCBA(NCENTC,J,J) = D0
            IF (NO_CAA) THEN
               DO I = 1, NUCDEP
                  IF(I.NE.J) FACCBA(NCENTC,I,I) = D0
               END DO
            END IF
            IF (NO_CCA) THEN
               DO I = 1, NUCDEP
                  IF (I.NE.J) THEN
                     FACCBA(NCENTC,J,I) = D0
                     FACCBA(NCENTC,I,J) = D0
                  END IF
               END DO
            END IF
            IF (NO_CBA) THEN
               DO K = 1, NUCDEP
                  IF (K.NE.J) THEN
                     DO I = 1, NUCDEP
                           IF (I.NE.J) FACCBA(NCENTC,I,K) = D0
                     END DO
                  END IF
               END DO
            END IF
            IF (IPRINT.GT.5) THEN
                WRITE(LUPRI,'(/4X,A)') 'FACCBA for this center:'
                CALL OUTPUT(FACCBA(NCENTC,1,1),1,1,1,N2DEP,
     &                      MXCENTC,N2DEP,-1,LUPRI)
                WRITE(LUPRI,'()')
            END IF
         ELSE
           IF (IPRINT.GE.2) WRITE(LUPRI,'(4X,I3,1P,D30.12,A)')
     &         IJ,QTMP, '         skipped!'
           QSKIP = QSKIP + ABS(QTMP)
         END IF
      END DO
C
      IF (IPRINT.GE.2) THEN
         WRITE(LUPRI,'(/A,I17)') 
     &               '     Number of charges skipped : ', NNDEP-NCENTC
         WRITE(LUPRI,'(/A,1P,D17.9)') 
     &               '     Amount of charge skipped  : ', QSKIP
         WRITE(LUPRI,'(A/)') 
     &               '     ( sum_i |Q_skip(i)| )'
      ENDIF
      IF (NCENTC.NE.MXCENTC) CALL QUIT('TYPE3DR1 : NCENTC .NE. MXCENTC')
      CALL QEXIT('TYPE3DR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type4drv */
      SUBROUTINE TYPE4DRV(OP1INT,WORK,LWORK,DOINT,OMITVNUC,
     &                     QMAT,IPRINT) 
C*****************************************************************************
C
C    TYPE4DRV : Driver for coulomb attraction integrals,
C               used to mimic the two-electron integrals
C               neglected in the one-center model 4
C
C               Written April 2002 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
      LOGICAL DOINT(2,2)
      CHARACTER*4 OMITVNUC
      DIMENSION OP1INT(NNBBASX),WORK(LWORK),QMAT(NUCDEP,NUCDEP,2)
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
C
      CALL QENTER('TYPE4DRV')
c
c     Start timer
c
      TIMSTR = SECOND()
c
      IF(DOINT(1,2) .OR. DOINT(2,1) ) THEN
         CALL QUIT ('TYPE4DRV called with LS or SL type integrals.')
      END IF
#include "memint.h"
C
      CALL MEMGET('REAL',KVCMAT,  NNBBASX,WORK,KFREE,LFREE)
      LFCBA = NUCDEP**3
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
      CALL DZERO(WORK(KVCMAT),NNBBASX)
c   
c     ... set charges. ( Q = q_s)
c         (We use the nuclear exponents)
c
      CALL DZERO(WORK(KCHRG),NUCIND)
      IF (ICTLV1C(2) .EQ. 1) THEN 
         CALL SCDENS_2(WORK(KCHRG),1,NUCIND,0)
      ELSE
         IOFF = 1
         DO I = 1,NUCIND
            IF (IOFF .GT. NUCDEP) 
     &          CALL QUIT('TYPE1DRV : IOFF .GT. NUCDEP')
            DO J = 1,NUCDEP
               WORK(KCHRG-1+I) = WORK(KCHRG-1+I) + QMAT(J,IOFF,2)
            END DO
            IOFF = IOFF + FMULT(ISTBNU(I))
         ENDDO
      END IF
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from TYPE4DRV',-1)
        WRITE(LUPRI,'(/1X,A)')
     &       '* Charges used for correction  :'
              WRITE(LUPRI,'(4X,A4,F12.8)') 
     &       (NAMN(I),WORK(KCHRG-1+I),I = 1,NUCIND)
      ENDIF
c
c     ... Make list of nuclei.
c
      CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &            GNUEXP,WORK(KCOORC), WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),OMITVNUC,IPRINT)
c
c    ... Calculate integrals
c
      CALL CATDR1(WORK(KVCMAT),DUMMY,NCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LWRK,IPRINT,
     &            .FALSE.,0,NNBBASX,DOINT)
C
C     Symmetry unpack integrals
C
      CALL SYMUPK(WORK(KVCMAT),OP1INT,1,NNBBASX)
C
C     Memory deallocation
C
      CALL MEMREL('TYPE4DRV.CATDR1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
c     Stop timer
c
      TIMEND = SECOND()
      TIME   = TIMEND - TIMSTR
      WRITE (LUPRI,'()')
      CALL TIMTXT('>>> Time used in TYPE4DRV is',TIME,LUPRI)
      CALL QEXIT('TYPE4DRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type5drv */
      SUBROUTINE TYPE5DRV(OP1INT,WORK,LWORK,DOINT,OMITVNUC,IPRINT) 
C*****************************************************************************
C
C    TYPE5DRV : Driver for coulomb attraction integrals,
C               used to mimic the multicenter integrals
C               neglected in the one-center model 5
C
C               Written June 2002 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
      PARAMETER ( D1 = 1.0D0 , DP5 = 0.5D0 )
      LOGICAL DOINT(2,2)
      CHARACTER*4 OMITVNUC
      DIMENSION OP1INT(NNBBASX),WORK(LWORK)
#include "dcbbas.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
C
      CALL QENTER('TYPE5DRV')
c
c     Start timer
c
      TIMSTR = SECOND()
#include "memint.h"
C
      CALL MEMGET('REAL',KVCMAT,  NNBBASX,WORK,KFREE,LFREE)
      LFCBA = NUCDEP**3
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXGEXP,  NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
      CALL DZERO(WORK(KVCMAT),NNBBASX)
c   
c     ---->  Set up calculation of integrals.
c
      CALL DZERO(WORK(KCHRG),NUCIND)
      DO I = 1,NUCIND
         WORK(KXGEXP-1+I) = GNUEXP(I)
         WORK(KCHRG-1+I)  = D1
      ENDDO
c
c     ... Make list of nuclei.
c
      CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &            WORK(KXGEXP),WORK(KCOORC), WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),OMITVNUC,IPRINT)
c
c    ... Calculate integrals
c
      CALL CATDR1(WORK(KVCMAT),DUMMY,NCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LWRK,IPRINT,
     &            .FALSE.,0,NNBBASX,DOINT)
c
c     <----
c
c     ----> Make X-matrix
c           
c           X_{AB} = \sum_{A} <S_A | S_B > D^{LS}_{AB}
c 
c     <----
c
c     ----> Combine elements of integralsmatrix and X-matrix
c           to produce the elements we wish to add to the 
c           2-elec. Fock-matrix :
c
c           (La Lb | Sc Sd) ~ 0.5 * <La | (1/rc) | Lb> <Sc | Sd>
c                           + 0.5 * <La | (1/rd) | Lb> <Sc | Sd>
c
C
C     Symmetry unpack integrals
C
      CALL SYMUPK(WORK(KVCMAT),OP1INT,1,NNBBASX)
C
C     Memory deallocation
C
      CALL MEMREL('TYPE5DRV.CATDR1',WORK,KWORK,KWORK,KFREE,LFREE)
c
      TIMEND = SECOND() 
      TIME   = TIMEND - TIMSTR
      WRITE (LUPRI,'()')
      CALL TIMTXT('>>> Time used in TYPE5DRV is',TIME,LUPRI)
      CALL QEXIT('TYPE5DRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck catcorr */
      SUBROUTINE CATCORR(F1MAT,F2MAT,DMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     CATCORR : Generate two-index (nuclear-attraction like)  integrals 
C     to mimic two-electron integrals for the one-center models.
C
C     Written by hjaaj & jkp   -   oct. 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
      PARAMETER(D1 = 1.0D0)
C
      LOGICAL DOINT(2,2),ADDQIC(2),LBIT,ICVTEST,DOLSQ
      CHARACTER*4 OMITVNUC
      DIMENSION WORK(LWORK),F1MAT(*),F2MAT(*),DMAT(*)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "cbihr1.h"
C
      CALL QENTER('CATCORR')
#include "memint.h"
C
C     Initialization
C     ==============
C
      ICVTEST    = .FALSE.
      DOINT(1,1) = .FALSE.
      DOINT(2,2) = .FALSE.
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
C
C     Memory allocation
C     =================
C
      IF(INTV1C.EQ.3 .AND. ICTLV1C(2) .EQ. 0) THEN
          NMPOLOP = 4
      ELSE IF (INTV1C.EQ.1) THEN
         IF(ICTLV1C(2) .EQ. 2) THEN
            NMPOLOP = 1
         ELSE
            NMPOLOP = 0
         ENDIF
      ELSE
         NMPOLOP = 1
      END IF
      CALL MEMGET('REAL',KQMAT,NUCDEP*NUCDEP*2*NMPOLOP,
     &             WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEMT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     ============================================
C     Which integral classes are needed in CATCORR
C     ============================================
C
      ILLINT=0
      ISLINT=0
      ISSINT=0
      IF (LBIT(INTDEF,1)) ILLINT = 1
      IF (LBIT(INTDEF,2)) ISLINT = 1
      IF (LBIT(INTDEF,3)) ISSINT = 1
      INTTEST=ILLINT + 2*ISLINT + 4*ISSINT
      ICVTEST=(INTTEST.NE.7).OR.(INTTEST.NE.INTFLG)
      IF(INTFLG.LT.1 .AND. .NOT.SMLV1C) 
     & CALL QUIT('  ** ERROR IN CATCORR : '//
     &     'Called wihout any integral-classes specified !')
C
C     ================================================
c     Get matrix of small- and large-component charges
C     ================================================
C
      DOLSQ=(INTV1C.EQ.1.AND.ICTLV1C(2).EQ.2).OR.
     &      (INTV1C.EQ.2.AND.ICTLV1C(2).NE.1).OR.
     &      (INTV1C.EQ.3)                    .OR.
     &      (INTV1C.EQ.4.AND.ICTLV1C(2).NE.1)
      KREL = KFREE
      IF(DOLSQ) THEN
        IF(NOPEN.GE.1) THEN
          CALL QUIT('open shell not implemented correctly for '//
     &   'this one-center model yet')
          CALL MEMGET('REAL',KBUF,N2BBASX,WORK,KFREE,LFREE)
          CALL DCOPY(N2BBASX,DMAT,1,WORK(KBUF),1)
          DO ISHELL = 1, NOPEN
            JDMAT = 1+N2BBASXQ*ISHELL
            CALL DAXPY(N2BBASX,DF(ISHELL),DMAT(JDMAT),1,
     &                 WORK(KBUF),1)
          ENDDO
          CALL LSQMAT(WORK(KQMAT),WORK(KBUF),NMPOLOP,WORK(KFREE),
     &                LFREE,IPRINT)
        ELSE
          CALL LSQMAT(WORK(KQMAT),DMAT,NMPOLOP,WORK(KFREE),LFREE,
     &                IPRINT)
        ENDIF
      ENDIF
c
C     Initialize OMITVNUC
c     (see CATNUC-routine (ababus/her1cat.F) for detail.)
c
      OMITVNUC    = 'FFFF'
c
c     ===============
c     Determine model
c     ===============
c
      IF (ONECAP) THEN
         IF (INTV1C .EQ. 1) THEN
c           ==============================================
c           Model 1 
c           Correct < X | V_C | Y > type integrals, X.ne.Y
c           ==============================================
            DOINT(2,2) = .TRUE.
            OMITVNUC    = 'TTFF'
            WRITE(LUPRI,'(//1X,A)')
     &      '<<<  ENTERING TYPE1DRV FOR SS-TYPE INTEGRALS  >>>'
            CALL TYPE1DRV(WORK(KONEIN),WORK(KFREE),LFREE,DOINT,
     &                  OMITVNUC,WORK(KQMAT),IPRINT)
            IF (.NOT.SMLV1C) THEN
c              ... Make full matrix - calculate energy
               CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
               CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
c           ... If smlv1c we dont have dmat at this point
               CALL EONECAP(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,2)
            ENDIF
         ELSE IF (INTV1C .EQ. 2) THEN
c           ==============================================
c           Model 2
c           Correct all but < C | V_C | C > type integrals
c           ==============================================
            DOINT(1,1)  = .TRUE.
            WRITE(LUPRI,'(//1X,A)') 
     & '<<<  ENTERING TYPE2DRV FOR APPROXIMATE LLSS-TYPE INTEGRALS  >>>'
c
c           > Check if one-center LS-integrals should be included 
c           > in this iteration
c
            IF(ICVTEST.AND.((INTFLG.NE.3).OR.(INTFLG.NE.7))) THEN
               OMITVNUC    = 'FFFF'
            ELSE
               OMITVNUC    = 'TFFF'
            ENDIF
            CALL TYPE2DRV(WORK(KONEIN),WORK(KFREE),LFREE,DOINT,
     &                  OMITVNUC,WORK(KQMAT),IPRINT)
c
c           ... calculate LL energy correction (Charge = q_s)
c
            CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL EONECAP(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,1)
c
c            correct <S|Vnuc|S> SS-type integrals
c
            DOINT(1,1)  = .FALSE.
            DOINT(2,2)  = .TRUE.
            WRITE(LUPRI,'(//1X,A)') 
     &      '<<<  ENTERING TYPE2DRV FOR <S|-Z|S> -TYPE INTEGRALS  >>>'
            CALL MEMGET('REAL',KONEIN2,NNBBASX,WORK,KFREE,LFREE)
c
c           > Check if one-center SS-integrals should be included 
c           > in this iteration
c
            IF(ICVTEST.AND.((INTFLG.NE.5).OR.(INTFLG.NE.7))) THEN
               OMITVNUC    = 'FTTT'
               CALL TYPE2DRV(WORK(KONEIN2),WORK(KFREE),LFREE,DOINT,
     &                       OMITVNUC,WORK(KQMAT),IPRINT)
               CALL DSPTSI(NTBAS(0),WORK(KONEIN2),WORK(KONEMT))
               CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
               CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,F2MAT,1)
               OMITVNUC    = 'TFFF'
            ELSE
               OMITVNUC    = 'TFFF'
            ENDIF
            CALL TYPE2DRV(WORK(KONEIN2),WORK(KFREE),LFREE,DOINT,
     &                  OMITVNUC,WORK(KQMAT),IPRINT)
c
c           ... calculate SS energy correction (Charge = -Z_nuc)
c
            CALL DSPTSI(NTBAS(0),WORK(KONEIN2),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL EONECAP(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,2)
c
            IF (ICTLV1C(2).EQ.3) THEN
c
c           Testing the addition of a penalty function (see type2drv)
c
               OMITVNUC = 'TTFF'
               CALL TYPE2DRV(WORK(KONEIN2),WORK(KFREE),LFREE,DOINT,
     &                       OMITVNUC,WORK(KQMAT),IPRINT)
               CALL DSPTSI(NTBAS(0),WORK(KONEIN2),WORK(KONEMT))
               CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
               CALL EONECAP(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,2)
            ENDIF
            CALL DAXPY(NNBBASX,D1,WORK(KONEIN2),1,WORK(KONEIN),1)
c
            IF (ICTLV1C(2).EQ.2) THEN
c
c           Estimate the error on Energy of this model <A | q_B^L + q_B^S - Z_B | C>
c
               DOINT(1,1)  = .FALSE.
               DOINT(2,2)  = .TRUE.
               OMITVNUC = 'TTFF'
               CALL TYPE2DRV(WORK(KONEIN2),WORK(KFREE),LFREE,DOINT,
     &                       OMITVNUC,WORK(KQMAT),IPRINT)
               CALL DSPTSI(NTBAS(0),WORK(KONEIN2),WORK(KONEMT))
               CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
               CALL EONECAP(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,3)
            ENDIF
C           ... Make full matrix of LL and SS corrections.
            CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL MEMREL('CATCORR.ONECAP2',WORK,KWORK,KONEIN2,KFREE,
     &                   LFREE)
         ELSE IF (INTV1C .EQ. 3) THEN
c           ==============================================
c           ONECAP type 3
c           Correct all but < C | V_C | C > type integrals
c           ==============================================
            OMITVNUC = 'TFFF'
            DOINT(1,1)  = .TRUE.
            DOINT(2,2)  = .FALSE.
            ADDQIC(1) = .FALSE.
            ADDQIC(2) = .TRUE.
            WRITE(LUPRI,'(//1X,A)') 
C
     &      '<<<  ENTERING TYPE3DRV FOR LL-TYPE INTEGRALS  >>>'
C
            CALL TYPE3DRV(WORK(KONEIN),WORK(KFREE),LFREE,DOINT,
     &                  ADDQIC,OMITVNUC,IPRINT,WORK(KQMAT))
c
c           ... calculate energy correction due to q_s
c
            CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL EONECAP(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,1)
c           
            OMITVNUC = 'TFFF'
            DOINT(1,1)  = .FALSE.
            DOINT(2,2)  = .TRUE.
            ADDQIC(1) = .TRUE.
            ADDQIC(2) = .TRUE.
            CALL MEMGET('REAL',KONEIN2,NNBBASX,WORK,KFREE,LFREE)
            WRITE(LUPRI,'(//1X,A)') 
C
     &      '<<<  ENTERING TYPE3DRV FOR SS-TYPE INTEGRALS  >>>'
C
            CALL TYPE3DRV(WORK(KONEIN2),WORK(KFREE),LFREE,DOINT,
     &                  ADDQIC,OMITVNUC,IPRINT,WORK(KQMAT))
c
c           ... calculate energy correction 
c
            CALL DSPTSI(NTBAS(0),WORK(KONEIN2),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL EONECAP(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,2)
c
c           ... Add LS and SS correction. (added to Fock-matrix later)
c
            CALL DAXPY(NNBBASX,D1,WORK(KONEIN2),1,WORK(KONEIN),1)
            CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL MEMREL('CATCORR.TYPE3',WORK,KWORK,KONEIN2,KFREE,
     &                   LFREE)
         ELSE IF (INTV1C .EQ. 4) THEN
c           ==============================================
c           ONECAP type 4
c           we approximate all LS and SS class integrals
c           with < A | V_C | B > type integrals
c           ==============================================
c
            DOINT(1,1)  = .TRUE.
            WRITE(LUPRI,'(//1X,A)') 
     &      '<<<  ENTERING TYPE4DRV FOR LL-TYPE INTEGRALS  >>>'
c
            OMITVNUC    = 'FFFF'
            CALL TYPE4DRV(WORK(KONEIN),WORK(KFREE),LFREE,DOINT,
     &                  OMITVNUC,WORK(KQMAT),IPRINT)
c
c           ... calculate LL energy correction (Charge = q_s)
c
            CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL EONECAP(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,1)
c
c           ... correct one-electron Fock-matrix
c
            CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,F1MAT,1)
            IF(IPRINT.GE.6) THEN
               CALL HEADER('CATCORR: ONECAP corrected F1MAT',-1)
               CALL PRQMAT(F1MAT,NTBAS(0),NTBAS(0),
     &                     NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
            END IF
            DOINT(1,1)  = .FALSE.
            DOINT(2,2)  = .TRUE.
            WRITE(LUPRI,'(//1X,A)') 
     &      '<<<  ENTERING TYPE4DRV FOR SS-TYPE INTEGRALS  >>>'
            OMITVNUC    = 'FFFF'
            CALL TYPE4DRV(WORK(KONEIN),WORK(KFREE),LFREE,DOINT,
     &                  OMITVNUC,WORK(KQMAT),IPRINT)
c
c           ... calculate SS energy correction (Charge = q_s)
c
            CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL EONECAP(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,2)
	     ELSE IF (INTV1C .EQ. 5) THEN
	        CALL QUIT(' ONECAP TYPE 5 NOT IMPLEMENTED YET!')
c           ... TO BE IMPLEMENTED
c           ONECAP type 5 :
c           Calculate one-center LS and SS explicitly.
c           Approximate multi-center LLSS contributions by :
c           (La Lb | Sc Sd) ~ 0.5 * <La | (1/rc) | Lb> <Sc | Sd>
c                           + 0.5 * <La | (1/rd) | Lb> <Sc | Sd>
c           Likewise for SSSS contributions.
            DOINT(1,1) = .TRUE.
			OMITVNUC = 'TFFF'
c            CALL TYPE5DRV(WORK(KONEIN),WORK(KFREE),LFREE,DOINT,
c     &                  OMITVNUC,IPRINT)
            DOINT(1,1) = .FALSE.
            DOINT(2,2) = .TRUE.
			OMITVNUC = 'TFFF'
c            CALL TYPE5DRV(WORK(KONEIN),WORK(KFREE),LFREE,DOINT,
c     &                  OMITVNUC,IPRINT)
         ELSE
	        CALL QUIT(' INVALID ONECAP MODEL SPECIFIED !')
         END IF
      ELSE
         CALL QUIT('CATCORR not implemented for .not. ONECAP')
      END IF
C     ====================
C     Print section
C     ====================
      IF(IPRINT.GE.6) THEN
         IF (INTV1C .EQ. 1) THEN
             CALL HEADER('CATCORR: ONECAP type 1 correction to '//
     &                    'nuclear attraction integrals',-1)
         ELSE IF (INTV1C .EQ. 2) THEN
             CALL HEADER('CATCORR: ONECAP type 2 correction to '//
     &                    'nuclear attraction integrals',-1)
         ELSE IF (INTV1C .EQ. 3) THEN
             CALL HEADER('CATCORR: ONECAP type 3 correction to '//
     &                    'nuclear attraction integrals',-1)
         END IF
         CALL PRQMAT(WORK(KONEMT),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     =============================
C     Add correction to Fock-matrix
C     =============================
C
      IF (INTV1C.EQ.3 .OR. INTV1C.EQ.4) THEN
c        ... Correct two-electron Fock-matrix
         CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,F2MAT,1)
         IF(IPRINT.GE.6) THEN
            CALL HEADER('CATCORR: ONECAP corrected F2MAT',-1)
            CALL PRQMAT(F2MAT,NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
         END IF
      ELSE IF (INTV1C.EQ.1 .OR. INTV1C.EQ.2) THEN
c        ... Correct one-electron Fock-matrix
         IF (SMLV1C) THEN
            CALL DAXPY(NNBBASX,D1,WORK(KONEIN),1,F1MAT,1)
            IF(IPRINT.GE.6) THEN
               CALL HEADER('CATCORR: F1MAT corrected for SMLV1C',-1)
               CALL OUTPAK(F1MAT,NTBAS(0),1,LUPRI)
            END IF
         ELSE
            CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,F1MAT,1)
            IF(IPRINT.GE.6) THEN
               CALL HEADER('CATCORR: ONECAP corrected F1MAT',-1)
               CALL PRQMAT(F1MAT,NTBAS(0),NTBAS(0),
     &                     NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
            END IF
         ENDIF
	  ELSE
	     CALL QUIT('Unknown .ONECAP type (INTV1C)')
      END IF
CC
CC     Memory deallocation
      CALL MEMREL('CATCORR',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('CATCORR')
      RETURN
      END
C*****************************************************************************
