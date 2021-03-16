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
C  /* Deck onesop */
      SUBROUTINE ONESOP(STHMAT,DENMAT,FOCMAT,FACINT,COORC,WORK,LWORK,
     &                  IPRINT,PROPTY,MAXDIF,IDENB0,CORBX0,CORBY0,
     &                  CORBZ0,DIFDIP,SECDER,NATOMC,TOLOG,TOLS,JSYMC,
     &                  JCENTC,NCENTC,SIGNC,GEXP,NNBASX,FCM,TLMD,
     &                  LOMITVNUC)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
C
#ifdef PRG_DIRAC
C     jth - length of STDER1 is now KCKTAB*6*2
      PARAMETER(NDER1=6,NDER2=6)
C     ... note that STDER2 is not implemented yet, and will require more than 6
#else
      PARAMETER(NDER1=3,NDER2=6)
#endif
C
      LOGICAL SECDER, DIFDIP, PROPTY, LOMITVNUC(2,5)
      DIMENSION WORK(LWORK)
      DIMENSION STHMAT(NNBASX,3), DENMAT(*), FOCMAT(*),
     &          FACINT(*), COORC(3,*), JSYMC(*),
     &          JCENTC(*), NCENTC(*), SIGNC(3,*),GEXP(*),
     &          FCM(*), TLMD(*)
#include "cbisol.h"
#include "onecom.h"
C
      IF (PROPTY) THEN
         NADER = 28*NATOMC*KCKTAB
      ELSE
         NADER = KCKTAB
      END IF
      KSDER0 = 1
      KSDER1 = KSDER0 +  2*KCKTAB
      KSDER2 = KSDER1 +  2*KCKTAB*NDER1
      KADER  = KSDER2 +  2*KCKTAB*NDER2
      KSINT0 = KADER  +    NADER
      KDINT1 = KSINT0 +    KCKTAB
      KDSHEL = KDINT1 +  9*KCKTAB
#ifdef PRG_DIRAC
C     jth -  make length of DSHELL = 4*KHKTAB (quaternion)
      KFSHEL = KDSHEL +  4*KHKTAB
#else 
      KFSHEL = KDSHEL +    KHKTAB
#endif 
      KLAST  = KFSHEL +    KHKTAB
      IF (SOLVNT) THEN
         IF (MAXDIF .EQ. 0) THEN
            KRLMNT = KLAST
            KRLMTB = KRLMNT + LMNTOT*KCKTAB
            KLAST  = KRLMTB
            KLMNO  = KLAST
         ELSE IF (MAXDIF .EQ. 1) THEN
            KRLMNT = KLAST
            KRLMTB = KRLMNT + 7*LMNTOT*KCKTAB
            KLAST  = KRLMTB
            KLMNO  = KLAST
         ELSE
            KRLMNT = KLAST
            KRLMTB = KRLMNT + 7*LMNTOT*KCKTAB
            KLMNO  = KRLMTB + 21*KCKTAB
            KLAST  = KLMNO  + (3*LMNTOT + 1)/IRAT
         END IF
      ELSE
         KRLMNT = KLAST
         KRLMTB = KLAST
         KLMNO  = KLAST
      END IF
      LWRK   = LWORK  - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ONESOP',' ',KLAST,LWORK)
C
      CALL ONESO1(STHMAT,WORK(KSDER0),WORK(KSDER1),WORK(KSDER2),
     &            WORK(KADER),WORK(KSINT0),WORK(KDINT1),
     &            DENMAT,FOCMAT,FACINT,COORC,WORK(KLAST),LWRK,IPRINT,
     &            PROPTY,MAXDIF,IDENB0,CORBX0,CORBY0,CORBZ0,DIFDIP,
     &            SECDER,NATOMC,TOLOG,TOLS,JSYMC,JCENTC,NCENTC,SIGNC,
     &            GEXP,NNBASX,WORK(KDSHEL),WORK(KFSHEL),WORK(KRLMNT),
     &            WORK(KRLMTB),FCM,TLMD,WORK(KLMNO),LOMITVNUC)
      RETURN
      END
C  /* Deck oneso1 */
      SUBROUTINE ONESO1(STHMAT,STDER0,STDER1,STDER2,ADER,SINT0,
     &                  DINT1,DENMAT,FOCMAT,FACINT,COORC,WORK,LWORK,
     &                  IPRINT,PROPTY,MAXDIF,IDENB0,CORBX0,CORBY0,
     &                  CORBZ0,DIFDIP,SECDER,NATOMC,TOLOG,TOLS,JSYMC,
     &                  JCENTC,NCENTC,SIGNC,GEXP,NNBASX,DSHELL,FSHELL,
     &                  RLMINT,RLMTAB,FCM,TLMD,LMNO,LOMITVNUC)
      use quaternion_algebra
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "iratdef.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
C
      LOGICAL SECDER, DIFDIP, PROPTY, FULMAT, ANTI, LOMITVNUC(2,5)
      DIMENSION WORK(LWORK)
#ifdef PRG_DIRAC
      PARAMETER(NDER1=6,NDER2=6)
C     ... note that STDER2 is not implemented yet, and will require more than 6
#else
      PARAMETER(NDER1=3,NDER2=6)
C     only DSHELL(*,1) used for Dalton
#endif
      DIMENSION STHMAT(NNBASX,3),STDER0(KCKTAB,2),
     &          STDER1(KCKTAB,NDER1,2),STDER2(KCKTAB,NDER2,2),ADER(*),
     &          DSHELL(KHKTAB,4),FSHELL(KHKTAB),
     &          DENMAT(*),FOCMAT(*),
     &          SINT0(KCKTAB), DINT1(KCKTAB,3,3),
     &          FACINT(NUCDEP,*), COORC(3,*), JSYMC(*),
     &          JCENTC(*), NCENTC(*), SIGNC(3,*),GEXP(*),
     &          RLMINT(*), RLMTAB(*), FCM(*), TLMD(*), LMNO(*)
#include "ccom.h"
#include "onecom.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbisol.h"
C
#ifdef PRG_DIRAC
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbgrd.h"
#else
#include "energy.h"
#endif 
#include "ibtfun.h"
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
#ifdef PRG_DIRAC
C     Element (I,J,IZ) in full quaternion matrix
      IQFULL(I,J,IZ) = MAX(I,J) + MIN(I-1,J-1) * NTBAS(0) +
     &     N2BBASX * (IZ-1) 
      IOLDQFULL(I,J,IZ) = I + (J-1) * NTBAS(0) + N2BBASX * (IZ-1) 
#endif 
      IF (IPRINT .GT. 4) CALL TITLER('Output from ONESO1','*',103)
C
C     *****************************************
C     ***** Loop over symmetry operations *****
C     *****************************************
C
      IDENB = IDENB0 - KHKTB
      DO 100 ISYMOP = 0, MAXOPR
      IF (IBTAND(ISYMOP,MAB) .NE. 0) THEN
         IF ( IBTAND ( ISYMOP, MULB ) .EQ. 0 ) THEN
            IDENB  = IDENB + KHKTB
         END IF
         GOTO 100
      ELSE
         IDENB  = IDENB + KHKTB
         ICENTB = NUCNUM(NCENTB,ISYMOP+1)
         ONECEN = ICENTA .EQ. ICENTB
         SIGNBX = PT(IBTAND(ISYMAX(1,1),ISYMOP))
         SIGNBY = PT(IBTAND(ISYMAX(2,1),ISYMOP))
         SIGNBZ = PT(IBTAND(ISYMAX(3,1),ISYMOP))
         CORBX  = SIGNBX*CORBX0
         CORBY  = SIGNBY*CORBY0
         CORBZ  = SIGNBZ*CORBZ0
         IF (IPRINT .GE. 05) WRITE (LUPRI, 1010) ISYMOP
         IF (IPRINT .GE. 10) THEN
            WRITE (LUPRI,'(A,I8)')    ' INTCLASS', INTCLASS
            WRITE (LUPRI,'(A,3F15.6)') ' CORBX/Y/Z ', CORBX,CORBY,CORBZ
         END IF
C
C *****
C   Omit selected terms depending on LOMITVNUC(INTTYP,1-5)
C   INTTYP = 1 : Specifies which terms to omit from the LL-part
C                of the nuclear-attraction integrals.
C   INTTYP = 2 : Specifies which terms to omit from the SS-part
C
C   LOMITVNUC(INTTYP,1) true : omit < C | V_C | C > 
C   LOMITVNUC(INTTYP,2) true : omit < A | V_C | A > 
C   LOMITVNUC(INTTYP,3) true : omit < A | V_C | C > and  < C | V_C | B > 
C   LOMITVNUC(INTTYP,4) true : omit < A | V_C | B > 
C   LOMITVNUC(INTTYP,5) true : any of the above true
C /hjaaj July 2001
C   added seperate treatment of SS- and LL-part     /jkp Sep 2001
C *****
         IFACINT = 1
         IF (INTCLASS .EQ. 0) THEN
            INTTYP = 1
         ELSE
            INTTYP = 2
         END IF
         IF (LOMITVNUC(INTTYP,5)) THEN
            IFACINT = 2
            IF (ONECEN) THEN
               IF (LOMITVNUC(INTTYP,1) .AND. LOMITVNUC(INTTYP,2)) THEN
                  CALL DZERO(FACINT(1,2),NATOMC)
               ELSE IF (LOMITVNUC(INTTYP,1)) THEN
                  CALL DCOPY(NATOMC,FACINT(1,1),1,FACINT(1,2),1)
C                 Omit < C | V_C | C >
                  DO I = 1,NATOMC
                     IF (NCENTC(I) .EQ. ICENTA)
     &                  FACINT(I,2) = D0
                  END DO
               ELSE IF (LOMITVNUC(INTTYP,2)) THEN
                  CALL DZERO(FACINT(1,2),NATOMC)
C                 Include < C | V_C | C >
                  DO I = 1,NATOMC
                     IF (NCENTC(I) .EQ. ICENTA)
     &                  FACINT(I,2) = FACINT(I,1)
                  END DO
               ELSE
                  IFACINT = 1
               END IF
            ELSE
C           not ONECEN
               IF (LOMITVNUC(INTTYP,3) .AND. LOMITVNUC(INTTYP,4)) THEN
                  CALL DZERO(FACINT(1,2),NATOMC)
               ELSE IF (LOMITVNUC(INTTYP,3)) THEN
                  CALL DCOPY(NATOMC,FACINT(1,1),1,FACINT(1,2),1)
C                 Omit < A | V_C | C > and < C | V_C | B >
                  DO I = 1,NATOMC
                     IF (NCENTC(I).EQ.ICENTA.OR.NCENTC(I).EQ.ICENTB)
     &                  FACINT(I,2) = D0
                  END DO
               ELSE IF (LOMITVNUC(INTTYP,4)) THEN
                  CALL DZERO(FACINT(1,2),NATOMC)
C                 Include < A | V_C | C > and < C | V_C | B >
                  DO I = 1,NATOMC
                     IF (NCENTC(I).EQ.ICENTA.OR.NCENTC(I).EQ.ICENTB)
     &                  FACINT(I,2) = FACINT(I,1)
                  END DO
               ELSE
                  IFACINT = 1
               END IF
            END IF
         ENDIF
C
C        **********************************************
C        ***** Calculation of Cartesian integrals *****
C        **********************************************
C
         CALL ONEPRM(STDER0,STDER1,STDER2,ADER,SINT0,DINT1,
     &               RLMINT,RLMTAB,FCM,WORK,LWORK,IPRINT,PROPTY,MAXDIF,
     &               NATOMC,TOLOG,TOLS,SECDER,DIFDIP,FACINT(1,IFACINT),
     &               COORC,GEXP,JCENTC,NCENTC)
C
C        *************************************************
C        ***** Transform to spherical harmonic basis *****
C        *************************************************
C
         IF (SPHRAB) THEN
            IF (PROPTY) THEN
               IF (ONECEN) THEN
                  CALL SPHRM1(ADER,ADER,10*NATOMC,WORK,LWORK,PROPTY,
     &                        IPRINT)
               ELSE
                  CALL SPHRM1(STDER0,STDER0,2,WORK,LWORK,PROPTY,IPRINT)
                  CALL SPHRM1(STDER1,STDER1,2*NDER1,
     &                        WORK,LWORK,PROPTY,IPRINT)
                  IF (SECDER) THEN
                     CALL SPHRM1(STDER2,STDER2,2*NDER2,
     &                           WORK,LWORK,PROPTY,IPRINT)
                  END IF
                  CALL SPHRM1(ADER,ADER,28*NATOMC,WORK,LWORK,PROPTY,
     &                        IPRINT)
                  IF (DIFDIP) THEN
                     CALL SPHRM1(DINT1,DINT1,9,WORK,LWORK,PROPTY,IPRINT)
                  END IF
               END IF
               CALL SPHRM1(SINT0,SINT0,1,WORK,LWORK,PROPTY,IPRINT)
               IF (SOLVNT) THEN
                  CALL SPHRM1(RLMINT,RLMINT,7*LMNTOT,WORK,LWORK,PROPTY,
     &                        IPRINT)
                  IF (SECDER) THEN
                     CALL SPHRM1(RLMTAB,RLMTAB,21,WORK,LWORK,PROPTY,
     &                           IPRINT)
                  END IF
               END IF
            ELSE
               CALL SPHRM1(STDER0,STDER0,2,WORK,LWORK,PROPTY,IPRINT)
               CALL SPHRM1(ADER,ADER,1,WORK,LWORK,PROPTY,IPRINT)
               IF (SOLVNT) THEN
                  CALL SPHRM1(RLMINT,RLMINT,LMNTOT,WORK,LWORK,PROPTY,
     &                        IPRINT)
               END IF
            END IF
         END IF
C
C        ******************************
C        ***** Expectation values *****
C        ******************************
C
         IF (PROPTY) THEN
C
C           One-electron Hamiltonian integrals and reorthonormalization
C           ===========================================================
C
C           Collect density and Fock elements
C
            MAXCMP = 0
            DO 200 IORBA = IDENA + 1, IDENA + KHKTA
               DO 250 IORBB = IDENB + 1, IDENB + KHKTB
                  MAXCMP = MAXCMP + 1
#ifdef PRG_DIRAC
Cjth              Due to the way FOCMAT and DENMAT are constructed in Dirac don't 
C                 do this...
C                 FAC = D1
Cjth - sib        IF (ONECEN .AND. LDIAG .AND. IORBB.NE.IORBA) FAC = DP5
                  IF (INTCLASS .EQ. 1) THEN
                     DO IZ = 2, 4
                        IDXAB = IQFULL(IORBA,IORBB,IZ)
                        DSHELL(MAXCMP,IZ) = DENMAT(IDXAB)
                     END DO
                  ELSE
                     IDXAB = IQFULL(IORBA,IORBB,1)
                     IDXAB = IOLDQFULL(IORBA,IORBB,1)
                     FSHELL(MAXCMP)   = FOCMAT(IDXAB)
                     DSHELL(MAXCMP,1) = DENMAT(IDXAB)
                  END IF
#else 
                  FAC = D1
                  IF (ONECEN .AND. LDIAG .AND. IORBB.NE.IORBA) FAC = DP5
                  IORBAB = ITRI(IORBA,IORBB)
                  DSHELL(MAXCMP) = FAC*DENMAT(IORBAB)
                  FSHELL(MAXCMP) = FAC*FOCMAT(IORBAB)
#endif 
  250          CONTINUE
  200       CONTINUE
C
#ifdef PRG_DIRAC
            IF (IPRINT .GT. 15) THEN
               CALL HEADER('Input density matrix',-1)
               CALL PRQMAT(DENMAT,NTBAS(0),NTBAS(0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
               CALL HEADER('Collected density matrix DSHELL',-1)
               IF (INTCLASS .EQ. 1) THEN
                  IZ_END = 4
               ELSE
                  IZ_END = 1
               END IF
               DO IZ = 1,IZ_END
                  WRITE(LUPRI,'(A,I3/)') ' *** IZ =',IZ
                  WRITE(LUPRI,'(10F12.6)') DSHELL(1:MAXCMP,IZ)
               END DO
            END IF
C
C           SL,LS : derivatives of the momentum
C           LL,SS : reorthonormalization and nuclear attraction
C           SS    : \beta-matrix derivatives
C 
            IF ( INTCLASS .EQ. 1 ) THEN
               IF ( .NOT. ONECEN ) THEN
C
C                 Derivatives of the momentum
C
                  CALL AVEMOM(STDER1,ISYMOP,MAXCMP,SECDER,DSHELL)
               END IF
            ELSE
C
C              GRD_NOSMLV: No small-small nuclear attraction 
C                          (eg. when using Levy-Leblond Hamiltonian)
C              GRD_NOSMLS: No small-small reorthonormalization term
C                          (eg. when using Levy-Leblond Hamiltonian)
C
               IF ( ONECEN ) THEN
C
C                 Derivatives of the nuclear attraction, one center case
C
                  IF (.NOT. (GRD_NOSMLV .AND. (INTCLASS .EQ. 2))) THEN
                     CALL AVENA1(ADER,NATOMC,SECDER,NCENTC,JCENTC,
     &                           MAXCMP,JSYMC,SIGNC,DSHELL)
                  END IF
C
               ELSE
C
C                 Derivatives of the nuclear attraction, two center case
C
                  IF (.NOT. (GRD_NOSMLV .AND. (INTCLASS .EQ. 2))) THEN
                     CALL AVENA2(ADER,NATOMC,ISYMOP,SECDER,NCENTC,
     &                           MAXCMP,JCENTC,JSYMC,SIGNC,DSHELL)
                  END IF
C
C                 Reorthonormalization
C
                  IF (.NOT. (GRD_NOSMLS .AND. (INTCLASS .EQ. 2))) THEN
                     CALL AVEFS(STDER1,ISYMOP,MAXCMP,SECDER,FSHELL)
                  END IF
C
C                 Derivatives of the \beta-matrix
C 
                  IF ( INTCLASS .EQ. 2 ) THEN
                     CALL AVEBET(STDER1,ISYMOP,MAXCMP,SECDER,DSHELL)
                  END IF
               END IF
            END IF
            IF (IPRINT .GT. 20) THEN
               IF ( DODERM  .AND. .NOT. ONECEN ) THEN
                  KCSTRA = 1
                  KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
                  KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
                  CALL HEADER('Momentum integral gradient',-1)
                  CALL PRIGRD(GRADKN(1,1),WORK(KCSTRA),WORK(KSCTRA))
               ELSE
                  KCSTRA = 1
                  KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
                  KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
                  CALL HEADER(
     &            'Nuclear attraction integral gradient - LL',-1)
                  CALL PRIGRD(GRADNU(1,1),WORK(KCSTRA),WORK(KSCTRA))
                  CALL HEADER(
     &            'Nuclear attraction integral gradient - SS',-1)
                  CALL PRIGRD(GRADNU(1,2),WORK(KCSTRA),WORK(KSCTRA))
                  IF (.NOT. ONECEN ) THEN
                     CALL HEADER(
     &               'Reorthonormalization gradient - LL',-1)
                     CALL PRIGRD(GRADRO(1,1),WORK(KCSTRA),WORK(KSCTRA))
                     CALL HEADER(
     &               'Reorthonormalization gradient - SS',-1)
                     CALL PRIGRD(GRADRO(1,2),WORK(KCSTRA),WORK(KSCTRA))
                     CALL HEADER('Beta matrix integral gradient',-1)
                     CALL PRIGRD(GRADKN(1,2),WORK(KCSTRA),WORK(KSCTRA))
                  END IF
               END IF
            END IF
#else 
            IF (ONECEN) THEN
C
C              Nuclear attraction
C
               CALL AVENA1(ADER,NATOMC,SECDER,NCENTC,JCENTC,MAXCMP,
     &                     JSYMC,SIGNC,DSHELL)
            ELSE
C
C              Kinetic energy and reorthonormalization
C
               CALL AVEKFS(STDER0,STDER1,STDER2,ISYMOP,MAXCMP,SECDER,
     &                     DSHELL,FSHELL)
C
C              Nuclear attraction
C
               CALL AVENA2(ADER,NATOMC,ISYMOP,SECDER,NCENTC,MAXCMP,
     &                     JCENTC,JSYMC,SIGNC,DSHELL)
            END IF
            IF (IPRINT .GT. 50) THEN
               CALL HEADER('Kinetic energy integral gradient',-1)
               CALL PRIGRD(GRADKE)
               CALL HEADER('Nuclear attraction integral gradient',-1)
               CALL PRIGRD(GRADNA)
               CALL HEADER('Reorthonormalization gradient',-1)
               CALL PRIGRD(GRADFS)
            END IF
#endif
C
C           Dipole gradient
C           ===============
C
            CALL AVEDIP(SINT0,DINT1,ISYMOP,DIFDIP,DSHELL,MAXCMP)
C
C           Solvent contributions
C           =====================
#if !defined (PRG_DIRAC)
            IF (SOLVNT) THEN
                CALL AVESOL(RLMINT,RLMTAB,FCM,MAXDIF,ISYMOP,
     &                      MAXCMP,DSHELL,TLMD,LMNO,IPRINT)
            END IF
#endif
         END IF
C
C        *******************************************
C        ***** Transform integrals to SO basis *****
C        *******************************************
C
         FULMAT = .TRUE.
         ANTI   = .FALSE.
C
C        Overlap integrals
C        =================
C
         CALL SYM1S(STDER0(1,1),STHMAT(1,1),ISYMOP,MULA,MULB,NHKTA,
     &              NHKTB,KHKTA,KHKTB,HKAB,LDIAG,FULMAT,DUMMY,IDUMMY,
     &              IPRINT)
C
C        Kinetic energy integrals
C        ========================
C
         CALL SYM1S(STDER0(1,2),STHMAT(1,3),ISYMOP,MULA,MULB,NHKTA,
     &              NHKTB,KHKTA,KHKTB,HKAB,LDIAG,FULMAT,DUMMY,IDUMMY,
     &              IPRINT)
C
C        Nuclear attraction integrals
C        ============================
C
         CALL SYM1S(ADER,STHMAT(1,2),ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &              KHKTA,KHKTB,HKAB,LDIAG,FULMAT,DUMMY,IDUMMY,
     &              IPRINT)
C
C        **************************************************
C        ***** Write differentiated integrals on file *****
C        **************************************************
C
         IF (IPRINT .GT. 5) THEN
            NADER = 28*NATOMC
            IF (ONECEN) NADER = 10*NATOMC
            CALL HEADER('First derivative overlap matrix',-1)
            CALL OUTPUT(STDER1(1,1,1),1,KHKTAB,1,3,KCKTAB,3,1,LUPRI)
            CALL HEADER('First derivative kinetic energy matrix',-1)
            CALL OUTPUT(STDER1(1,1,2),1,KHKTAB,1,3,KCKTAB,3,1,LUPRI)
            CALL HEADER('ADER',-1)
            CALL OUTPUT(ADER,1,KHKTAB,1,NADER,KCKTAB,NADER,1,LUPRI)
         END IF
         IF (PROPTY .AND. (SECDER .OR. DIFDIP)) THEN
            CALL WD1SYM(STDER1,ADER,RLMINT,FCM,WORK,LWORK,JSYMC,JCENTC,
     &                  ISYMOP,NATOMC,IPRINT)
         END IF
      END IF
  100 CONTINUE
      RETURN
 1010 FORMAT (//,2X,'***************************************',
     *         /,2X,'******** Symmetry operation ',I2,' ********',
     *         /,2X,'***************************************',/)
      END
c
C  /* Deck wd1sym */
      SUBROUTINE WD1SYM(STDER1,ADER,RLMINT,FCM,WORK,LWORK,JSYMC,JCENTC,
     &                  ISYMOP,NATOMC,IPRINT)
#include "implicit.h"
#include "priunit.h"
      DIMENSION STDER1(KCKTAB,3,2), ADER(*), JSYMC(*), JCENTC(*),
     &          RLMINT(*), FCM(*), WORK(LWORK)
#include "onecom.h"
C
#ifdef PRG_DIRAC
C     HJAaj Nov 2002: this routine is not programmed for Dirac (yet)
      CALL QUIT('WD1SYM called in Dirac, but not programmed for Dirac')
#endif
      KOMAT = 1
      KLAST  = KOMAT + 6*KCKTAB
      IF (KLAST .GT. LWORK) CALL STOPIT('WD1SYM',' ',KLAST,LWORK)
      CALL WD1SY1(STDER1,ADER,RLMINT,FCM,WORK(KOMAT),JSYMC,JCENTC,
     &            ISYMOP,NATOMC,IPRINT)
      RETURN
      END
C  /* Deck wd1sy1 */
      SUBROUTINE WD1SY1(STDER1,ADER,RLMINT,FCM,OMAT,JSYMC,JCENTC,ISYMOP,
     &                  NATOMC,IPRINT)
#include "implicit.h"
#include "priunit.h"
      DIMENSION STDER1(KCKTAB,3,2), ADER(KCKTAB,NATOMC,*), JSYMC(*),
     &          JCENTC(*), OMAT(KCKTAB,3,2), RLMINT(KCKTAB,LMNTOT,7),
     &          FCM(*)
#include "cbisol.h"
#include "onecom.h"
#include "ader.h"
C
C     Write differentiated integrals on file
C
C     Overlap and kinetic energy integrals
C     ====================================
C
      IF (.NOT.ONECEN) THEN
         CALL DCOPY(6*KCKTAB,STDER1,1,OMAT,1)
         CALL DRSYM1(OMAT(1,1,1),OMAT(1,1,2),NCENTA,NCENTB,ISYMOP,MULA,
     &               MULB,NHKTA,NHKTB,KHKTA,KHKTB,KHKTAB,KCKTAB,HKAB,
     &               LDIAG,IPRINT)
      END IF
C
C     Nuclear attraction integrals
C     ============================
C
      DO 200 IATOMC = 1, NATOMC
         DO 300 ICMPAB = 1, KHKTAB
            OMAT(ICMPAB,1,1) = ADER(ICMPAB,IATOMC,IA000X)
            OMAT(ICMPAB,2,1) = ADER(ICMPAB,IATOMC,IA000Y)
            OMAT(ICMPAB,3,1) = ADER(ICMPAB,IATOMC,IA000Z)
            IF (.NOT.ONECEN) THEN
               OMAT(ICMPAB,1,2) = ADER(ICMPAB,IATOMC,IA0X00)
               OMAT(ICMPAB,2,2) = ADER(ICMPAB,IATOMC,IA0Y00)
               OMAT(ICMPAB,3,2) = ADER(ICMPAB,IATOMC,IA0Z00)
            END IF
  300    CONTINUE
         ISYMC  = JSYMC(IATOMC)
         ICENTC = JCENTC(IATOMC)
         CALL DASYM1(OMAT(1,1,1),OMAT(1,1,2),ONECEN,NCENTA,NCENTB,
     &               ICENTC,ISYMOP,ISYMC,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &               KHKTB,KHKTAB,KCKTAB,HKAB,LDIAG,IPRINT)
  200 CONTINUE
C
C     Solvent contribution
C     ====================
C
      IF (SOLVNT) THEN
        DO 400 ICMPAB = 1, KHKTAB
          OMAT(ICMPAB,1,1)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,5),KCKTAB)
          OMAT(ICMPAB,2,1)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,6),KCKTAB)
          OMAT(ICMPAB,3,1)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,7),KCKTAB)
          OMAT(ICMPAB,1,2)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,2),KCKTAB)
          OMAT(ICMPAB,2,2)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,3),KCKTAB)
          OMAT(ICMPAB,3,2)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,4),KCKTAB)
  400   CONTINUE
        CALL DASYM1(OMAT(1,1,1),OMAT(1,1,2),ONECEN,NCENTA,NCENTB,
     &              NCNTCV,ISYMOP,0,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &              KHKTB,KHKTAB,KCKTAB,HKAB,LDIAG,IPRINT)
      END IF
      RETURN
      END
C  /* Deck dsym1 */
      SUBROUTINE DSYM1(DENMAT,FOCMAT,DSO,FSO,NBAST,IPRINT)
C
C     Take density matrix in symmetry orbital basis and generate
C     density matrix over distinct pairs of AOs
C
C                                          880418  PRT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION DSO(*), FSO(*), DENMAT(*), FOCMAT(*)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) CALL HEADER('Subroutine DSYM1',-1)
C
C     Loop over all irreps in molecule
C
      ISOFF = 0
      ISTR = 1
      NNBASX = NBAST*(NBAST + 1)/2
      CALL DZERO(DENMAT,NNBASX)
      CALL DZERO(FOCMAT,NNBASX)
      DO 100 IREP = 0, MAXREP
         NORBI = NAOS(IREP+1)
         IF (NORBI .EQ. 0) GOTO 110
         DO 200 I = ISTR,ISTR + NORBI - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            IOFF = KSTRT(IA)
            MULA = ISTBAO(IA)
            INDA = IOFF + NA
            DO 300 J = ISTR,I
               IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
               NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
               JOFF   = KSTRT(IB)
               NHKTB  = NHKT(IB)
               KHKTB  = KHKT(IB)
               MULB   = ISTBAO(IB)
               MAB    = IBTOR(MULA,MULB)
               KAB    = IBTAND(MULA,MULB)
               HKAB   = FMULT(KAB)
               ISOFF  = ISOFF + 1
               DSYMIJ = DSO(ISOFF)
               FSYMIJ = FSO(ISOFF)
               INDB   = JOFF + NB - KHKTB
               DO 400 ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MAB) .NE. 0) GOTO 400
                  INDB = INDB + KHKTB
C
C                 Weight and parity factor
C
                  FAC = HKAB*
     *                  PT(IBTAND(ISYMOP,IBTXOR(IREP,ISYMAO(NHKTB,NB))))
                  INDM = MAX(INDA,INDB)
                  IND  = (INDM*(INDM - 3))/2 + INDA + INDB
                  DENMAT(IND) = DENMAT(IND) + FAC*DSYMIJ
                  FOCMAT(IND) = FOCMAT(IND) + FAC*FSYMIJ
400            CONTINUE
300         CONTINUE
200      CONTINUE
110      CONTINUE
         ISTR = ISTR + NORBI
100   CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Total density matrix (sym. distinct AO basis)',-1)
         CALL OUTPAK(DENMAT,NBAST,1,LUPRI)
         CALL HEADER('Total Fock matrix (sym. distinct AO basis)',-1)
         CALL OUTPAK(FOCMAT,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck drsym1 */
      SUBROUTINE DRSYM1(SAODER,TAODER,ICENTA,ICENTB,ISYMOP,MULA,MULB,
     &                  NHKTA,NHKTB,KHKTA,KHKTB,KHKTAB,KCKTAB,HKAB,
     &                  LDIAG,IPRINT)
C
C     Arrange calculation of symmetry-adapted integral derivatives
C     from distinct AO integral derivatives
C                                                PRT & TUH  880428
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "csym1.h"
#include "ccom.h"
      DIMENSION SAODER(KCKTAB,*), TAODER(KCKTAB,*)
      LOGICAL LDIAG, FULMAT, ANTI
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine DRSYM1',-1)
         WRITE (LUPRI,'(A,2I5)') ' ICENTA/B ', ICENTA, ICENTB
         WRITE (LUPRI,'(A, I5)') ' ISYMOP   ', ISYMOP
         WRITE (LUPRI,'(A,2I5)') ' NHKTA/B  ', NHKTA, NHKTB
         WRITE (LUPRI,'(A,2I5)') ' KHKTA/B  ', KHKTA, KHKTB
         WRITE (LUPRI,'(A,2I5)') ' MULA/B   ', MULA, MULB
         WRITE (LUPRI,'(A,F12.6)') ' HKAB   ', HKAB
         WRITE (LUPRI,'(A, L5)') ' LDIAG    ', LDIAG
      END IF
      FULMAT = .FALSE.
      ANTI   = .FALSE.
      NMATS = 3*NUCDEP*(MAXREP+1)
      DO 10 ICL = 1,2
C
C        Determine factors to account for use of transl. invariance
C
         IF (ICL .EQ. 1) THEN
            JCENT = ICENTA
         ELSE
            JCENT = ICENTB
         ENDIF
         MULJ = ISTBNU(JCENT)
C
C     Run over Cartesian directions
C
      DO 20 JCDIR = 1,3
         ISYTYJ = ISYMAX(JCDIR,1)
C
C     Run over irrep's of the differentiation operator
C
      DO 30 IREPD = 0,MAXREP
         IF (IBTAND(MULJ,IBTXOR(IREPD,ISYTYJ)) .EQ. 0) THEN
            IF (ICL .EQ. 1) THEN
               FAC = HKAB
            ELSE
               FAC = - HKAB*PT(IBTAND(ISYTYJ,ISYMOP))
     *                     *PT(IBTAND(IREPD, ISYMOP))
            END IF
            IF (IPRINT .GT. 20) THEN
               WRITE (LUPRI, '(/A,3I5)') ' ICL, JCDIR, IREPD ',
     *                                     ICL, JCDIR, IREPD
               WRITE (LUPRI, '(A,2I5)') ' MULJ, ISYTYJ ', MULJ, ISYTYJ
               WRITE (LUPRI, '(A,F12.6)') ' FAC ', FAC
            END IF
            IMAT0 = (IPTCNT(3*(JCENT-1)+JCDIR,IREPD,1)-1)*(MAXREP+1)+1
            IF (IREPD .EQ. 0) THEN
C
C              Overlap matrix - totally symmetric perturbation
C
               CALL SYM1S(SAODER(1,JCDIR),DUMMY,ISYMOP,MULA,MULB,
     *                    NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,FULMAT,
     *                    THRS,IMAT0,IPRINT)
C
C              Kinetic energy  - totally symmetric perturbation
C
               IMAT0 = NMATS + IMAT0
               CALL SYM1S(TAODER(1,JCDIR),DUMMY,ISYMOP,MULA,MULB,
     *                    NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,FULMAT,
     *                    THRS,IMAT0,IPRINT)
            ELSE
C
C              Overlap matrix - non-symmetric perturbation
C
               CALL SYM1N(SAODER(1,JCDIR),DUMMY,IREPD,ISYMOP,MULA,
     *                    MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     *                    FULMAT,ANTI,THRS,IMAT0,IPRINT)
C
C              Kinetic energy  - non-symmetric perturbation
C
               IMAT0 = NMATS + IMAT0
               CALL SYM1N(TAODER(1,JCDIR),DUMMY,IREPD,ISYMOP,MULA,
     *                    MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     *                    FULMAT,ANTI,THRS,IMAT0,IPRINT)
            END IF
         END IF
30    CONTINUE
20    CONTINUE
10    CONTINUE
      RETURN
      END
C  /* Deck dasym1 */
      SUBROUTINE DASYM1(CAODER,AAODER,ONECEN,ICENTA,ICENTB,ICENTC,
     &                  ISYMOP,JSYMOP,MULA,MULB,NHKTA,NHKTB,
     &                  KHKTA,KHKTB,KHKTAB,KCKTAB,HKAB,LDIAG,
     &                  IPRINT)
C
C     Arrange calculation of symmetry-adapted integral derivatives
C     from distinct AO nuclear attraction integral derivatives
C                                                PRT & TUH  880502
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      LOGICAL LDIAG, ONECEN
      DIMENSION CAODER(KCKTAB,*), AAODER(KCKTAB,*)
      LOGICAL FULMAT, ANTI
#include "ccom.h"
#include "csym1.h"
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine DASYM1',-1)
         WRITE (LUPRI,'(A,2I5)') ' ICENTA/B ', ICENTA, ICENTB
         WRITE (LUPRI,'(A, I5)') ' ISYMOP   ', ISYMOP
         WRITE (LUPRI,'(A, I5)') ' ICENTC   ', ICENTC
         WRITE (LUPRI,'(A, I5)') ' JSYMOP   ', JSYMOP
         WRITE (LUPRI,'(A,2I5)') ' NHKTA/B  ', NHKTA, NHKTB
         WRITE (LUPRI,'(A,2I5)') ' KHKTA/B  ', KHKTA, KHKTB
         WRITE (LUPRI,'(A,2I5)') ' MULA/B   ', MULA, MULB
         WRITE (LUPRI,'(A,F12.6)') ' HKAB   ', HKAB
         WRITE (LUPRI,'(A, L5)') ' LDIAG    ', LDIAG
         WRITE (LUPRI,'(A, L5)') ' ONECEN   ', ONECEN
      END IF
      FULMAT = .FALSE.
      ANTI   = .FALSE.
      NMATS = 3*NUCDEP*(MAXREP+1)
      IF (ONECEN) THEN
         ICLMX = 2
      ELSE
         ICLMX = 3
      ENDIF
      DO 10 ICL = 1,ICLMX
C
C        Determine factors to account for use of transl. invariance
C
         IF (ICL .EQ. 1) THEN
            JCENT = ICENTC
         ELSE IF (ICL .EQ. 2) THEN
            JCENT = ICENTA
         ELSE
            JCENT = ICENTB
         ENDIF
         MULJ = ISTBNU(JCENT)
C
C     Run over Cartesian directions
C
      DO 20 JCDIR = 1,3
         ISYTYJ = ISYMAX(JCDIR,1)
         IF (ICL .EQ. 3) THEN
            DO 50 I = 1, KHKTAB
               CAODER(I,JCDIR) = CAODER(I,JCDIR) + AAODER(I,JCDIR)
50          CONTINUE
         ENDIF
C
C     Run over irreps of the differentiation operator
C
      DO 30 IREPD = 0,MAXREP
         IF (IBTAND(MULJ,IBTXOR(IREPD,ISYTYJ)) .EQ. 0) THEN
            IF (ONECEN) THEN
               IF (ICL .EQ. 1) THEN
                  FAC = HKAB*PT(IBTAND(ISYTYJ,JSYMOP))
     *                      *PT(IBTAND(IREPD ,JSYMOP))
               ELSE
                  FAC = - HKAB
               ENDIF
            ELSE
               IF (ICL .EQ. 1) THEN
                  FAC = HKAB*PT(IBTAND(ISYTYJ,JSYMOP))
     *                      *PT(IBTAND(IREPD ,JSYMOP))
               ELSE IF (ICL .EQ. 2) THEN
                  FAC = HKAB
               ELSE
                  FAC = - HKAB*PT(IBTAND(ISYTYJ,ISYMOP))
     *                        *PT(IBTAND(IREPD, ISYMOP))
               END IF
            END IF
            IF (IPRINT .GT. 20) THEN
               WRITE (LUPRI, '(/A,3I5)') ' ICL, JCDIR, IREPD ',
     *                                     ICL, JCDIR, IREPD
               WRITE (LUPRI, '(A,2I5)') ' MULJ, ISYTYJ ', MULJ, ISYTYJ
               WRITE (LUPRI, '(A,F12.6)') ' FAC ', FAC
            END IF
            IMAT0 = NMATS
     *            + (IPTCNT(3*(JCENT-1)+JCDIR,IREPD,1)-1)*(MAXREP+1)+1
            IF (IREPD .EQ. 0) THEN
               IF (ICL .EQ. 2 .AND. .NOT.ONECEN) THEN
                  CALL SYM1S(AAODER(1,JCDIR),DUMMY,ISYMOP,MULA,MULB,
     *                       NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,FULMAT,
     *                       THRS,IMAT0,IPRINT)
               ELSE
                  CALL SYM1S(CAODER(1,JCDIR),DUMMY,ISYMOP,MULA,MULB,
     *                       NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,FULMAT,
     *                       THRS,IMAT0,IPRINT)
               END IF
            ELSE
               IF (ICL .EQ. 2 .AND. .NOT.ONECEN) THEN
                  CALL SYM1N(AAODER(1,JCDIR),DUMMY,IREPD,ISYMOP,MULA,
     *                       MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     *                       FULMAT,ANTI,THRS,IMAT0,IPRINT)
               ELSE
                  CALL SYM1N(CAODER(1,JCDIR),DUMMY,IREPD,ISYMOP,MULA,
     *                       MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     *                       FULMAT,ANTI,THRS,IMAT0,IPRINT)
               END IF
            END IF
         END IF
30    CONTINUE
20    CONTINUE
10    CONTINUE
      RETURN
      END
C  /* Deck sym1s */
      SUBROUTINE SYM1S(AO,SO,KB,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     *                 HKAB,LDIAG,FULMAT,THRESH,IMAT0,IPRINT)
C
C     Take block of distinct AO integral (derivatives) and
C     generate symmetrized contributions to SO integral
C     (derivatives) for the totally symmetric case
C                                          880407  PRT
C     Modified tuh 880819
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
      PARAMETER (LUTEMP = 48, IREPO = 0)
      LOGICAL LDIAG, FULMAT
      DIMENSION AO(*), SO(*)
#include "csym1.h"
#include "symmet.h"
#include "symind.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine SYM1S',-1)
         WRITE (LUPRI, '(A, I5,A,L1)') ' KB ',KB,', FULMAT = ',FULMAT
         WRITE (LUPRI, '(A,2I5,A,F20.10)')
     &   ' NHKTA/B ',NHKTA, NHKTB,',  HKAB = ',HKAB
      END IF
C
C     Loop over all irreps in molecule
C
      DO 100 IREP = 0, MAXREP
         IF (FULMAT) THEN
            INDOFF = NPARSU(IREP + 1)
         ELSE
            IMAT = IMAT0 + IREP
         END IF
C
C        Loop over AOs which are of symmetry IREP in stabilizer MULA
C
         DO 200 NA = 1, KHKTA
         IF (IBTAND(MULA,IBTXOR(IREP,ISYMAO(NHKTA,NA))).EQ.0) THEN
            NAT = KHKTB*(NA - 1)
            IF (LDIAG) THEN
               KHKTBB = NA
            ELSE
               KHKTBB = KHKTB
            END IF
            DO 300 NB = 1,KHKTBB
            IF (IBTAND(MULB,IBTXOR(IREP,ISYMAO(NHKTB,NB))).EQ.0) THEN
C
C              Weight and parity factor
C
               FAC = HKAB*PT(IBTAND(KB,IBTXOR(IREP,ISYMAO(NHKTB,NB))))
C
C              Locate SO integrals to which AO's contribute
C
               INDA = INDFA(IREP + 1,NA)
               INDB = INDFB(IREP + 1,NB)
               INDM = MAX(INDA,INDB)
               RINT = FAC*AO(NAT+NB)
               IF (FULMAT) THEN
                  IND  = INDOFF + (INDM*(INDM - 3))/2 + INDA + INDB
                  SO(IND) = SO(IND) + RINT
                  IF (IPRINT.GT.10) WRITE(LUPRI,*)
     &               'INDA/B,IND,RINT ',INDA,INDB,IND,RINT,SO(IND)
               ELSE
                  IF (ABS(RINT) .GT. THRESH) THEN
                     IND  = (INDM*(INDM - 3))/2 + INDA + INDB
                     INDMAX = MAX(IND,INDMAX)
                     LABEL = IND*2**16 + IMAT
                     IF (IPRINT .GT. 20) THEN
                       WRITE (LUPRI,'(A,F12.6,2I3,2I2,I4,I2,I5)')
     *                     'SYM1S - NA/B,IREPA/B,IND,IREPO,IMAT',
     *                     RINT, NA, NB, IREP, IREP, IND, IREPO,
     *                     IMAT
                     END IF
                     IF (LENGTH .EQ. 600) THEN
                        WRITE (LUTEMP) BUF, IBUF, LENGTH
                        IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I4,A)')
     *                     ' Buffer of length',LENGTH,
     *                     ' has been written in SYM1S.'
                        LENGTH = 0
                     ENDIF
                     LENGTH = LENGTH + 1
                     BUF (LENGTH) = RINT
                     IBUF(LENGTH) = LABEL
                  ENDIF
               ENDIF
            END IF
300         CONTINUE
         END IF
200      CONTINUE
100   CONTINUE
      RETURN
      END
C  /* Deck sym1n */
      SUBROUTINE SYM1N(AO,SO,IREPO,KB,MULA,MULB,NHKTA,NHKTB,
     *                 KHKTA,KHKTB,HKAB,LDIAG,FULMAT,ANTI,THRESH,
     *                 IMAT0,IPRINT)
C
C     Take block of distinct AO integral (derivatives) and
C     generate symmetrized contributions to SO integral
C     (derivatives) over non-symmetric operators
C                                          880408  PRT
C     Modified tuh 880819
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      LOGICAL LDIAG, FULMAT, ANTI
      DIMENSION AO(*), SO(*)
      PARAMETER (LUTEMP = 48)
#include "csym1.h"
#include "symmet.h"
#include "symind.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine SYM1N',-1)
         WRITE (LUPRI, '(A,2I5,A,L1)') ' IREPO, KB ', IREPO, KB,
     &      ', FULMAT = ',FULMAT
         WRITE (LUPRI, '(A,2I5,A,F20.10)')
     &   ' NHKTA/B ',NHKTA, NHKTB,',  HKAB = ',HKAB
         WRITE(lupri,*) 'MULA/B, KHKTA/B',MULA,MULB,KHKTA,KHKTB
         WRITE(lupri,*) 'LDIAG,ANTI,THRESH,IMAT0,IPRINT',
     &                   LDIAG,ANTI,THRESH,IMAT0,IPRINT
      END IF
C
C     Loop over irreps for first basis function - those for second
C     are obtained from operator symmetry IREPO
C
      DO 100 IREPA = 0, MAXREP
         IREPB = IBTXOR(IREPO,IREPA)
         IF (FULMAT) INDOFF = NPARNU(IREPO+1,MAX(IREPA,IREPB)+1)
         IF (ANTI .AND. (IREPA .LT. IREPB)) THEN
            FAB = - HKAB
         ELSE
            FAB = HKAB
         END IF
C
C        Loop over AOs which are of symmetry IREPA in stabilizer MULA
C
         DO 200 NA = 1, KHKTA
         IF (IBTAND(MULA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))).EQ.0) THEN
            NAT = KHKTB*(NA - 1)
            IF (LDIAG) THEN
               KHKTBB = NA
            ELSE
               KHKTBB = KHKTB
            ENDIF
            DO 300 NB = 1,KHKTBB
            IF (NA.EQ.NB .AND. LDIAG .AND. IREPA .LT. IREPB) GOTO 300
            IF (IBTAND(MULB,IBTXOR(IREPB,ISYMAO(NHKTB,NB))).EQ.0) THEN
C
C              Weight and parity factor
C
               FAC = FAB*PT(IBTAND(KB,IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
C
C              Locate SO integrals to which AOs contribute
C
               INDA = INDFA(IREPA + 1,NA)
               INDB = INDFB(IREPB + 1,NB)
               RINT = FAC*AO(NAT+NB)
               IF (FULMAT) THEN
                  IF (IREPA .GE. IREPB) THEN
                     IND  = INDOFF + NAOS(IREPB+1)*(INDA-1) + INDB
                  ELSE
                     IND  = INDOFF + NAOS(IREPA+1)*(INDB-1) + INDA
                  ENDIF
                  SO(IND) = SO(IND) + RINT
                  IF (IPRINT.GT.10) WRITE(LUPRI,*)
     &               'INDA/B,IND,RINT ',INDA,INDB,IND,RINT,SO(IND)
               ELSE
                  IF (ABS(RINT) .GT. THRESH) THEN
                     IF (IREPA .GE. IREPB) THEN
                        IND  = NAOS(IREPB+1)*(INDA-1) + INDB
                     ELSE
                        IND  = NAOS(IREPA+1)*(INDB-1) + INDA
                     ENDIF
                     INDMAX = MAX(IND,INDMAX)
                     LABEL  = IND*2**16 + IMAT0 + MAX(IREPA,IREPB)
                     IF (IPRINT .GT. 20) THEN
                       WRITE (LUPRI,'(A,F12.6,2I3,2I2,I4,I2,I5)')
     *                     'SYM1N - NA/B,IREPA/B,IND,IREPO,IMAT',
     *                     RINT, NA, NB, IREPA, IREPB, IND, IREPO,
     *                     IMAT0 + MAX(IREPA,IREPB)
                     END IF
                     IF (LENGTH .EQ. 600) THEN
                        WRITE (LUTEMP) BUF, IBUF, LENGTH
                        IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I4,A)')
     *                     ' Buffer of length',LENGTH,
     *                     ' has been written in SYM1N.'
                        LENGTH = 0
                     ENDIF
                     LENGTH = LENGTH + 1
                     BUF(LENGTH)  = RINT
                     IBUF(LENGTH) = LABEL
                  ENDIF
               ENDIF
            END IF
300         CONTINUE
         END IF
200      CONTINUE
100   CONTINUE
      RETURN
      END
C  /* Deck sphrm1 */
      SUBROUTINE SPHRM1(CI,SPI,NTYPE,WORK,LWORK,PROPTY,IPRINT)
#include "implicit.h"
      LOGICAL PROPTY
      DIMENSION CI(KCKTAB,*), SPI(KCKTAB,*), WORK(LWORK)
#include "onecom.h"
      KTMP  = 1
      KTMQ  = KTMP + KCKTA*KCKTB
      KLAST = KTMQ + KCKTA*KHKTB
      IF (KLAST .GT. LWORK) CALL STOPIT('SPHRM1',' ',KLAST,LWORK)
      DO 100 I = 1, NTYPE
         CALL SPHRMX(CI(1,I),SPI(1,I),WORK(KTMP),WORK(KTMQ),PROPTY,
     &               IPRINT)
  100 CONTINUE
      RETURN
      END
C  /* Deck sphrmx */
      SUBROUTINE SPHRMX(CI,SPI,TMPINT,HALF,PROPTY,IPRINT)
C
C     Transform a block of Cartesian integrals to spherical harmonics
C
C                                          920511  PRT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION CI(*), SPI(*), HALF(*), TMPINT(*)
      LOGICAL PROPTY
#include "onecom.h"
#include "sphtrm.h"
#include "symmet.h"
#include "ibtfun.h"
C
      IF (IPRINT .GT. 10) THEN
         WRITE(LUPRI,'(A,2I6)') ' KCKT? ', KCKTA, KCKTB
         WRITE(LUPRI,'(A,2L6)') ' SPHR? ', SPHRA, SPHRB
         WRITE(LUPRI,'(A,4I6)') ' KHKT? ', KHKTA, KHKTB
         CALL HEADER('Cartesian integrals in SPHRMX',-1)
         IJ = 0
         DO 100 I = 1, KCKTA
            DO 110 J = 1,KCKTB
               IJ = IJ + 1
               WRITE(LUPRI,'(2I4,D13.6)') I,J, CI(IJ)
 110        CONTINUE
 100     CONTINUE
      END IF
C
      IF (SPHRA .OR. SPHRB) THEN
C
C        Collect integrals to be transformed
C        ===================================
C
         ICOFF = 0
         DO 200 ICOMPA = 1,KCKTA
            DO 220 ICOMPB = 1,KCKTB
               ICOFF = ICOFF + 1
               TMPINT((ICOMPA-1)*KCKTB + ICOMPB) = CI(ICOFF)
  220       CONTINUE
  200    CONTINUE
C
C        Transform second index (B) if required
C        ======================================
C
         IF (SPHRB) THEN
            CALL MXM(CSP(ISPADR(NHKTB)),KHKTB,TMPINT,KCKTB,HALF,KCKTA)
         ELSE
            CALL DCOPY(KCKTA*KCKTB,TMPINT,1,HALF,1)
         END IF
C
C        Transpose half transformed integrals
C        ====================================
C
         DO 300 ICOMPA = 1,KCKTA
            DO 310 ICOMPB = 1,KHKTB
               TMPINT((ICOMPB-1)*KCKTA + ICOMPA)
     &               = HALF((ICOMPA-1)*KHKTB + ICOMPB)
  310       CONTINUE
  300    CONTINUE
C
C        Transform first index (A) if required
C        =====================================
C
         IF (SPHRA) THEN
            CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,TMPINT,KCKTA,HALF,KHKTB)
         ELSE
            CALL DCOPY(KCKTA*KHKTB,TMPINT,1,HALF,1)
         END IF
C
C        Collect transformed integrals
C        =============================
C
         ISOFF = 0
         DO 400 ICOMPA = 1,KHKTA
            DO 410 ICOMPB = 1,KHKTB
               ISOFF = ISOFF + 1
               SPI(ISOFF) = HALF((ICOMPB-1)*KHKTA + ICOMPA)
  410        CONTINUE
  400     CONTINUE
      END IF
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Spherical integrals in SPHRMX',-1)
         IJ = 0
         DO 500 I = 1, KHKTA
            DO 510 J = 1,KHKTB
               IJ = IJ + 1
               WRITE(LUPRI,'(A,(2I4,D13.6))') '++',I,J,SPI(IJ)
 510        CONTINUE
 500     CONTINUE
      END IF
      RETURN
      END


