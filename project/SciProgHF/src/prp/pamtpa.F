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

      SUBROUTINE PRPTPA(WORK,LWRK)
C***********************************************************************
C
C     Driver routine for computing two-photon absorption.
C
C     Written by johhe 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
C
      PARAMETER ( LUCMO = 22 )
      DIMENSION RES_LR(3,MAXLQR,MAXLQR,MAXFQR)
      DIMENSION RES_TPA(3,NBSYM,MAXLQR,MAXLQR,MAXFQR)
      DIMENSION WORK(*)
C
      CALL QENTER('PRPTPA')
C
      KFREE = 1
      LFREE = LWRK
      WRITE(LUPRI,'(///2A/)') 
     &     '*** ENTERING TWO-PHOTON ',
     &     'ABSORPTION CALCULATION'
C
C     Get the MO coefficients
C     =======================
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL OPNFIL(LUCMO,'DFCOEF','OLD','PRPTPA')
      IF(SPINFR) THEN
        CALL MEMGET('INTE',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
        CALL REACMO(LUCMO,'DFCOEF',WORK(KCMO),DUM,WORK(KIBRP),
     &              TOTERG,11)
      ELSE
        KIBRP = KFREE
        CALL REACMO(LUCMO,'DFCOEF',WORK(KCMO),DUM,IDUM,TOTERG,3)
      ENDIF
      CLOSE(LUCMO,STATUS='KEEP')
C
C     Calculate allowed excitations
C     =============================
C
      WRITE(LUPRI,'(/A)') ' '
      CALL PRSYMB(LUPRI,'=',70,1)
      WRITE(LUPRI,'(2A)') 
     &     ' >>>>>>>>         E X C I T A T I O N  ',
     &     'E N E R G I E S         <<<<<<<<'
      CALL PRSYMB(LUPRI,'=',70,1)
C
      CALL TPAEXC(WORK(KFREE),LFREE)
C
C     Set up the linear response equations to be solved
C     =================================================
C
      CALL TPA_SETUP(WORK(KFREE),LFREE)
C
      WRITE(LUPRI,'(/A)') ' '
      CALL PRSYMB(LUPRI,'=',70,1)
      WRITE(LUPRI,'(2A)') 
     &     ' >>>>>>>>              L I N E A R  R E S P O N S E',
     &     '            <<<<<<<<'
      CALL PRSYMB(LUPRI,'=',70,1)
C
      CALL QRVEC(WORK(KFREE),LFREE)
C
C     Calculate linear response function values
C     =========================================
C
      CALL MEMGET('REAL',KVECAEE,MZYEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECAEP,MZYEP*NZ,WORK,KFREE,LFREE)
      CALL LRCALC(RES_LR,WORK(KVECAEE),WORK(KVECAEP),WORK(KCMO),
     &            WORK(KIBRP),WORK(KFREE),LFREE)
C
C     Calculate quadratic response function values
C     ============================================
C
      WRITE(LUPRI,'(/A)') ' '
      CALL PRSYMB(LUPRI,'=',70,1)
      WRITE(LUPRI,'(2A)') 
     &     ' >>>>>>>>       T W O - P H O T O N   ',
     &     'A B S O R P T I O N      <<<<<<<<'
      CALL PRSYMB(LUPRI,'=',70,1)
C
      CALL MEMGET('REAL',KVECBEE,MZYEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECBEP,MZYEP*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECEEE,MZYEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECEEP,MZYEP*NZ,WORK,KFREE,LFREE)
      CALL QRTPA(RES_TPA,WORK(KVECAEE),WORK(KVECBEE),WORK(KVECEEE),
     &           WORK(KVECAEP),WORK(KVECBEP),WORK(KVECEEP),
     &           WORK(KCMO),WORK(KIBRP),WORK(KFREE),LFREE)
C
C     Printout of the result
C     ======================
C
      CALL TPA_PRINT(RES_LR,RES_TPA,WORK(KFREE),LFREE)
C
      CALL QEXIT('PRPTPA')
      RETURN
      END
C
      SUBROUTINE TPAEXC(WORK,LWRK)
C***********************************************************************
C
C     Calculation of excitation energies.
C
C     Written by johhe 2002  (modified by E.T., aug 2004)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dcbxrs.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
C
      DIMENSION WORK(*)
      LOGICAL DOEXC
C
#include "ibtfun.h"
      MULD2H(I,J) = IBTXOR(I-1,J-1) + 1
C
      CALL QENTER('TPAEXC')
C
      KFREE = 1
      LFREE = LWRK
C
C     Loop over excitation symmetries
C     ===============================
C
      NFC = NBSYM/NFSYM
      JTIMOP = 1
      DO IFERE = 1,NFSYM
      DO IS = 1,NFC
         JOPSY  = IFERE
         JSYMOP = JFSYM(IS,JOPSY)
C
C     Check for non-trivial excitations (boson irreps)
C     ================================================
C
         IF (TPACNV(JSYMOP).NE.0) THEN
            DOEXC = .FALSE.
            DO IAOP = 1, NAQROP
            DO IBOP = 1, NBQROP
               IF (MULD2H(IPRPSYM(LAQROP(IAOP)),
     &              MULD2H(IPRPSYM(LBQROP(IBOP)),JSYMOP)).EQ.1)
     &              DOEXC = .TRUE.
            END DO
            END DO
            IF (.NOT.DOEXC) THEN
               WRITE(LUPRI,*) 'WARNING: TPA cannot be '//
     &              'calculated in symmetry group no. ', JSYMOP
               TPACNV(JSYMOP) = 0
               TPASIM(JSYMOP) = 0
               TPASTV(JSYMOP) = 0
            END IF
         END IF
      END DO
      END DO
C
C     Symmetry check done; now compute the excitation vectors
      CALL CALCEXCVEC(WORK(KFREE),LFREE,TPACNV,TPASIM,TPASTV)
C
      CALL QEXIT('TPAEXC')
      RETURN
      END
C
      SUBROUTINE CALCEXCVEC(WORK,LWRK,NEXCSYM,NEXCSIM,NEXCSTV)
C***********************************************************************
C
C     Calculation of excitation vectors and energies.
C
C     Written by johhe 2002 (modified by E.T., aug 2004)
C     Input:  NEXCSYM - number of excited states (in each symmetry)
C             NEXCSIM - number of simultaneous solutions
C             NEXCSTV - number of start vectors
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dcbxrs.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
C
      DIMENSION RCNV(MAXFQR,NBSYM)
      DIMENSION WORK(*)
      INTEGER NEXCSYM(8),NEXCSIM(8),NEXCSTV(8)
C
      CALL QENTER('CALCEXCVEC')
C
      KFREE = 1
      LFREE = LWRK
C
C     Setup of needed variables
C     =========================
C
      CALL EXCSET(WORK,KFREE,LFREE)
C
C     ************************************************
C     *** E X C I T A T I O N   F U N C T I O N S  ***
C     ************************************************
C
C     Memory allocation
C     =================
C     Type of trial vectors, C-, E-, P-type trial vectors to IBTYP and
C     eigenvectors of reduced system respectively.
C
      CALL MEMGET('INTE',KIBTYP,NREDM*2     ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBCVC,NREDM       ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEVC,NREDM       ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBPVC,NREDM       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVECR,NREDM*MAXFQR,WORK,KFREE,LFREE)
C
C     Loop over excitation symmetries
C     ===============================
C
      NFC = NBSYM/NFSYM
      JTIMOP = 1
      DO IFERE = 1,NFSYM
      DO IS = 1,NFC
         JOPSY  = IFERE
         JSYMOP = JFSYM(IS,JOPSY)
         WRITE(LUPRI,'(///A,A3/A,I3,A)')
     &        ' >>>  SOLVING PP EXCITATION ENERGIES - Symmetry: ',
     &        REP(JSYMOP-1),
     &        '      for the',NEXCSYM(JSYMOP),
     &        ' lowest excitation energies.'
C
C     Check for excitations
C     =====================
C
         IF (NEXCSYM(JSYMOP).GT.0) THEN
            CALL XRSPAR(XQR_INDSTR,XQR_SKIPEE,XQR_SKIPEP,
     &                  IPRXQR)
            IF (NZCONF.EQ.0) THEN
               NEXMAX = NZXOPE
            ELSE
               NEXMAX = NZCONF
            END IF
            IF (NEXCSYM(JSYMOP) .GT. NEXMAX) THEN
               WRITE(LUPRI,'(/A/A,I5)')
     &              ' WARNING: no. of excitation energies exceeds '//
     &              'no. of var. The number is decreased to',NEXMAX
               NEXCSYM(JSYMOP) = NEXMAX
            END IF
            NEXSIM = MIN(NEXCSIM(JSYMOP), NEXMAX)
            NFREQ  = NEXSIM
            NEXSTV = MIN(NEXCSTV(JSYMOP), NEXMAX)
C
C     *******************************************************
C     *** S O L V E   E I G E N V A L U E   E Q. S Y S T. ***
C     *******************************************************
C
            INDPRP = IEXCOFF + JSYMOP
            NFREQ  = NEXCSYM(JSYMOP)
            NCRED  = 0
            NERED  = 0
            NPRED  = 0
            NZRED  = NCRED + NERED + NPRED
            RSPLAB = 'PP EXCITATIONS  '
            PRPNAM(INDPRP) = 'EXCLABL'
C     
C     Solve for eigenvalues
C     =====================
C     
            CALL XRSCTL(DUMMY,DUMMY,DUMMY,
     &           WORK(KIBTYP),WORK(KIBCVC),WORK(KIBEVC),
     &           WORK(KIBPVC),RCNV(1,JSYMOP),
     &           EXCIT(1,JSYMOP),WORK(KEVECR),
     &           WORK,KFREE,LFREE)
C     
C     Initializing indices (record number)  before calling LRVEC
C     ==========================================================
C     
            DO IFR=1,NFREQ
               CALL INDQR(INDPRP,EXCIT(IFR,JSYMOP),IDUM)
            END DO
C     
C     Writing excitation vectors to file
C     ==================================
C     
            CALL LRVEC(INDPRP,WORK(KIBTYP),WORK(KIBCVC),
     &           WORK(KIBEVC),WORK(KIBPVC),EXCIT(1,JSYMOP),
     &           WORK(KEVECR),WORK(KFREE),LFREE)
         END IF
      END DO
      END DO
C
      CALL QEXIT('CALCEXCVEC')
      RETURN
      END
C
      SUBROUTINE EXCSET(WORK,KFREE,LFREE)
C***********************************************************************
C
C     Some settings used in CALCEXCVEC.
C
C     Written by johhe 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxlr.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbxqr.h"
#include "dcbxrs.h"
#include "dcbgen.h"
C
      LOGICAL EXST
      DIMENSION WORK(*)
C
      CALL QENTER('EXCSET')
C
C     Initiate variables in dcbxqr.h
C     ==============================
C
      NQROP = 0
      CALL DZERO(EXCIT,MAXFQR*8) 
C
C     Set orbital strings
C     ===================
C
      CALL OSTRING(XQR_INDSTR,XQR_SKIPEE,XQR_SKIPEP,IPRXQR,
     &             WORK,KFREE,LFREE)
C
C     Initiate variables in dcbxrs.h
C     ==============================
C
      CALL SETRSP
      LSVCFG(1) = .TRUE.
      LSVCFG(2) = .TRUE.
      CNVINT(1) = CNVXQR(1)
      CNVINT(2) = CNVXQR(2)
      ITRINT(1) = ITRIQR(1)
      ITRINT(2) = ITRIQR(2)
      LINEQ  = .FALSE.
      TKNORM = XQRNRM
      DIAGHE = XQRDIH
      IPRXRS = IPRXQR
      THCXRS = THCQR
      RESFAC = RESXQR
      MAXITR = ITRXQR
      MAXSIM = -1
      NREDM  = MAXQRM
      N2REDM = MAXQRM*MAXQRM
      LOFFTY = 0
      INTDEF = INTXQR
      STERNH = XSTERN
      STATIC = .FALSE.
C
C     Get the maximum length of the response vectors
C     ==============================================
C
      MZYEE  = 0
      MZYEP  = 0
      JTIMOP = 1
      NFC    = NBSYM/NFSYM
      DO IFERE = 1,NFSYM
      DO IS = 1,NFC
         JOPSY = IFERE
         CALL XRSPAR(XQR_INDSTR,XQR_SKIPEE,XQR_SKIPEP,
     &               IPRXQR)
         MZYEE = MAX(MZYEE,2*NZXOPE)
         MZYEP = MAX(MZYEP,2*NZXOPP)
      END DO
      END DO
C
C     Initiate file QRINFO
C     ====================
C
      INQUIRE(FILE='QRINFO',EXIST=EXST) 
      IF (EXST) THEN 
         WRITE(LUPRI,'(A,2(/A,A),/A,/)') ' *** WARNING ***',
     &        ' Response vector found in file QRINFO, and the',
     &        ' program will continue',
     &        ' assuming you have supplied the corresponding',
     &        ' direct access files',' QRVEC.EE and QRVEC.PP.'
      ELSE
         OPEN(UNIT=LUQRINFO,FILE='QRINFO',STATUS='NEW')
         WRITE(LUQRINFO,'(A16)') ' **End of File**'
      END IF
C
      CALL QEXIT('EXCSET')
      RETURN
      END
C
      SUBROUTINE TPA_SETUP(WORK,LWRK)
C***********************************************************************
C
C     Determine which linear response equations to be solved. Similar to 
C     QR_SETUP, but differences in loop structure and considerations of
C     different symmetries.
C
C     Written by johhe 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxpr.h"
#include "dcbxlr.h"
#include "dcbxqr.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbxrs.h"
C
      DIMENSION WORK(*)
      PARAMETER ( D2 = 2.0D0 )
C
#include "ibtfun.h"
      MULD2H(I,J) = IBTXOR(I-1,J-1) + 1
C
      CALL QENTER('TPA_SETUP')
C
      KFREE = 1
      LFREE = LWRK
C
C     Loop over excitation symmetries
C     ===============================
C
      NQRHYP = 0
      NFC = NBSYM/NFSYM
      DO IFERE = 1,NFSYM
      DO IS = 1,NFC
         JOPSY  = IFERE
         JSYMEX = JFSYM(IS,JOPSY)
C
C     Loop over frequencies, A- and B-operators
C     =========================================
C     Due to symmetry between S_{ab} and S_{ba} it is only necessary to
C     six out of nine elements in the tensor S.
C
         DO IFR = 1,TPACNV(JSYMEX)
         DO IAOP = 1,NAQROP
         DO IBOP = IAOP,NBQROP
            IAPT = LAQROP(IAOP)
            IBPT = LBQROP(IBOP)
            ISYMA = IPRPSYM(IAPT)
            ISYMB = IPRPSYM(IBPT)
C
C     Check if allowed due to symmetry and loop over frequencies
C     ==========================================================
C
            IF (MULD2H(ISYMA,MULD2H(ISYMB,JSYMEX)).EQ.1) THEN
               NQRHYP = NQRHYP +1
               IF (NQRHYP.LE.MAXQR) THEN
                  FREQ = EXCIT(IFR,JSYMEX)/D2
                  CALL INDQR(IAPT,FREQ,INDEX)
                  CALL INDQR(IBPT,FREQ,INDEX)
                  LQRHYP(NQRHYP,1) = IAPT
                  LQRHYP(NQRHYP,2) = IBPT
                  QRFRHYP(NQRHYP,1) = -FREQ
                  QRFRHYP(NQRHYP,2) = -FREQ
               ELSE
                  NQRHYP = MAXQR
                  WRITE(LUPRI,'(//A,/A,I5)')
     & '@ WARNING!!! TPA not evaluated since the total number',
     & '             exceeds MAXQR=',MAXQR
               END IF
            END IF
         END DO
         END DO
         END DO
      END DO
      END DO
C
      CALL QEXIT('TPA_SETUP')
      RETURN
      END
C
      SUBROUTINE QRTPA(RESULT,VECAEE,VECBEE,VECEEE,VECAEP,VECBEP,VECEEP,
     &     CMO,IBEIG,WORK,LWRK)
C***********************************************************************
C
C     Calculation of elements in two-photon transition tensor.
C
C     Written by johhe 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dcbgen.h"
#include "dgroup.h"
C
      PARAMETER ( D0=0.0D0 )
      CHARACTER ALAB*8,BLAB*8,ELAB*8
      DIMENSION VECAEE(MZYEE,NZ),VECBEE(MZYEE,NZ),VECEEE(MZYEE,NZ),
     &          VECAEP(MZYEP,NZ),VECBEP(MZYEP,NZ),VECEEP(MZYEP,NZ)
      DIMENSION CMO(*),IBEIG(*),WORK(*)
      DIMENSION RESULT(3,NBSYM,MAXLQR,MAXLQR,MAXFQR)
C
#include "ibtfun.h"
      MULD2H(I,J) = IBTXOR(I-1,J-1) + 1
C
      CALL QENTER('QRTPA')
C
      KFREE = 1
      LFREE = LWRK
C
C     Open the property file
C     ======================
C
      OPEN(LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
C     
      CALL DZERO(RESULT,3*NBSYM*MAXLQR*MAXLQR*MAXFQR)
C
      IQRF = 0
      NFC  = NBSYM/NFSYM
C
C     Loop over excitation symmetries
C     ===============================
C
      DO IFERE = 1,NFSYM
      DO IS = 1,NFC
         JSYMEX = JFSYM(IS,IFERE)
         INDPRP = IEXCOFF + JSYMEX
C
C     Loop over frequencies, A- and B-operators
C     =========================================
C     Due to symmetry between S_{ab} and S_{ba} only necessary to
C     six out of nine elements in the tensor S.
C         
         DO IFR = 1,TPACNV(JSYMEX)
         DO IAOP = 1,NAQROP
         DO IBOP = IAOP,NBQROP
            ISYMA = IPRPSYM(LAQROP(IAOP))
            ISYMB = IPRPSYM(LBQROP(IBOP))
            IF (MULD2H(ISYMA,MULD2H(ISYMB,JSYMEX)).EQ.1) THEN
               IQRF = IQRF +1
               IF (IQRF.LE.MAXQR) THEN
C
C     Initialize A- and B-operators and excitations
C     =============================================
C
                  CALL HYPINIT(IQRF,1,AFREQ,IAPT,ISYMA,IFERA,ITIMA,
     &                         ALAB,NZYAEE,NZYAEP,VECAEE,VECAEP,
     &                         WORK(KFREE),LFREE)
                  CALL HYPINIT(IQRF,2,BFREQ,IBPT,ISYMB,IFERB,ITIMB,
     &                         BLAB,NZYBEE,NZYBEP,VECBEE,VECBEP,
     &                         WORK(KFREE),LFREE)
                  CALL EXCINIT(INDPRP,IFR,EFREQ,JSYMEX,IFERE,ITIME,
     &                         ELAB,NZYEEE,NZYEEP,VECEEE,VECEEP,
     &                         WORK(KFREE),LFREE)
C
                  WRITE(LUPRI,'(/A,2(I2,A),3(/A,A10,I4,F10.6))')
     &                 ' Two-photon absorption calculation no. ',
     &                 IQRF,' out of ',NQRHYP,' in a.u.',
     &                 ' A operator, boson symmetry, frequency: ',
     &                 ALAB,ISYMA,AFREQ,
     &                 ' B operator, boson symmetry, frequency: ',
     &                 BLAB,ISYMB,BFREQ,
     &                 ' Excitation: boson symmetry, frequency: ',
     &                 ELAB,JSYMEX,EFREQ
                  CALL FLSHFO(LUPRI)
C
C     Calculation and summation of the different terms in the element
C     ===============================================================
C
                  CALL A2DRV(IAPT,VECBEE,VECBEP,VECEEE,VECEEP,
     &                       ISYMA,ISYMB,JSYMEX,ITIMA,ITIMB,ITIME,
     &                       IFERA,IFERB,IFERE,
     &                       NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYEEE,NZYEEP,
     &                       RESULT(1,JSYMEX,IAOP,IBOP,IFR),CMO,
     &                       IBEIG,WORK(KFREE),LFREE)
C     
                  CALL A2DRV(IAPT,VECEEE,VECEEP,VECBEE,VECBEP,
     &                       ISYMA,JSYMEX,ISYMB,ITIMA,ITIME,ITIMB,
     &                       IFERA,IFERE,IFERB,
     &                       NZYAEE,NZYAEP,NZYEEE,NZYEEP,NZYBEE,NZYBEP,
     &                       RESULT(1,JSYMEX,IAOP,IBOP,IFR),CMO,
     &                       IBEIG,WORK(KFREE),LFREE)
C     
                  CALL X2DRV(VECAEE,VECAEP,IBPT,VECEEE,VECEEP,
     &                       ISYMA,ISYMB,JSYMEX,ITIMA,ITIMB,ITIME,
     &                       IFERA,IFERB,IFERE,
     &                       NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYEEE,NZYEEP,
     &                       RESULT(1,JSYMEX,IAOP,IBOP,IFR),CMO,
     &                       IBEIG,WORK(KFREE),LFREE)
C     
                  CALL T3DRV(VECAEE,VECAEP,VECBEE,VECBEP,
     &                       VECEEE,VECEEP,
     &                       ISYMA,ISYMB,JSYMEX,ITIMA,ITIMB,ITIME,
     &                       IFERA,IFERB,IFERE,
     &                       NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYEEE,NZYEEP,
     &                       BFREQ,EFREQ,
     &                       RESULT(1,JSYMEX,IAOP,IBOP,IFR),CMO,
     &                       WORK(KFREE),LFREE)
C
                  WRITE(LUPRI,'(3(/A,F20.12))')
     &                 ' Value of electronic part:',
     &                 RESULT(2,JSYMEX,IAOP,IBOP,IFR),
     &                 ' Value of positronic part:',
     &                 RESULT(3,JSYMEX,IAOP,IBOP,IFR),
     &                 ' Value of total response :',
     &                 RESULT(1,JSYMEX,IAOP,IBOP,IFR)
               END IF
            END IF
         END DO
         END DO
         END DO
      END DO
      END DO
C
C     Close the property file
C     =======================
C
      CLOSE(LU1INT,STATUS='KEEP')
C
      CALL QEXIT('QRTPA')
      RETURN
      END
C
      SUBROUTINE TPA_PRINT(RES_LR,RES_TPA,WORK,LWRK)
C***********************************************************************
C
C     Printing of results aquired during TPA calculations.
C
C     Written by johhe 2002
C
C***********************************************************************
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "pgroup.h"
C
      PARAMETER ( LUCMO = 22, D0 = 0.0D0, D2R3 = (2.0D0/3.0D0),
     &            D2 = 2.0D0, D4 = 4.0D0 )
      CHARACTER PFMT*6,MXFORM*6,ALAB*8,BLAB*8
      DIMENSION RES_LR(3,MAXLQR,MAXLQR,MAXFQR)
      DIMENSION RES_TPA(3,NBSYM,MAXLQR,MAXLQR,MAXFQR)
      DIMENSION WORK(*)
C
      INTEGER EXCOPFROMSYM
C
#include "ibtfun.h"
      MULD2H(I,J) = IBTXOR(I-1,J-1) + 1
C
      CALL QENTER('TPA_PRINT')
C
      KFREE = 1
      LFREE = LWRK
C
      CALL TITLER('RESULTS OF TWO-PHOTON ABSORBTION CALCULATION',
     &            '*',116)
      CALL HEADER('Linear response properties',-1)
C
C     ***********************************
C     *** Print excitation properties ***
C     ***********************************
C
      CALL HEADER('Excitation properties',0)
C
      NFC = NBSYM/NFSYM
      JTIMOP = 1
C
C     Loop over excitation symmetries
C     ===============================
C
      DO IFERE = 1,NFSYM
      DO IS = 1,NFC
         JOPSY  = IFERE
         JSYMEX = JFSYM(IS,JOPSY)
C
         CALL PRSYMB(LUPRI,'-',80,0)
         WRITE(LUPRI,'(/A,I2,A,A3/)')
     &        ' *** Excited state boson symmetry',
     &        JSYMEX,' : ',REP(JSYMEX-1)
         WRITE(LUPRI,'(/2X,A,I2)')
     &        'No. of excitations in the symmetry: ',
     &        TPACNV(JSYMEX)
C
C     Loop over frequencies
C     =====================
C
         DO IFR = 1,TPACNV(JSYMEX)
            FREQ = EXCIT(IFR,JSYMEX)
C
            CALL PRSYMB(LUPRI,'.',80,0)
            WRITE(LUPRI,'(/3X,A,I2,A,F14.6,A)')
     &           '* Excitation energy no. ',IFR,' : ',
     &           FREQ, ' a.u.'
            WRITE(LUPRI,'(32X,F14.6,A)')
     &           FREQ*XTEV, ' eV'
C
C     Loop over A- and B-operators
C     ============================
C
            IBOP = EXCOPFROMSYM(JSYMEX)
C
            OSCILL = D0
            DO IAOP=IBOP,NQROP
               IAPT  = LQROP(IAOP,1)
               ALAB  = PRPNAM(IAPT)
               ISYMA = IPRPSYM(IAPT)
               IBPT  = LQROP(IBOP,1)
               BLAB  = PRPNAM(IBPT)
C
C     Check that A-operator is not an excitation operator
C     ===================================================
C
               IF (IAPT.GE.IEXCOFF) GOTO 110
C
C     Check if it is a transition moment and print if allowed
C     =======================================================
C
                  TRNMOM = RES_LR(1,IAOP,IBOP,IFR)
                  FREQ = EXCIT(IFR,JSYMEX)
C
                  IF ((TPACNV(JSYMEX).NE.0).AND.
     &                (MULD2H(ISYMA,JSYMEX).EQ.1)) THEN
C
                     PFMT = MXFORM(TRNMOM,12)
                     WRITE(LUPRI,'(/5X,A,A10,'//PFMT//',A)')
     &                     ' Operator, transition moment: ',
     &                     ALAB,TRNMOM,' a.u.'
C
C     Calculate and oscillator strength
C     =================================
C
                     OSCILL = OSCILL +
     &                        D2R3*TRNMOM*TRNMOM*FREQ
                  END IF
 110           CONTINUE
            END DO
C
C     Print oscillator strength
C     =========================
C
            IF (OSCILL.GE.1.0D-10) THEN
               PFMT = MXFORM(OSCILL,12)
               WRITE(LUPRI,'(/5X,A,'//PFMT//')')
     &              ' Oscillator strength: ',OSCILL
            ELSE
               WRITE(LUPRI,'(/5X,A)')
     &              'Forbidden transition'
            END IF
         END DO
      END DO
      END DO
      CALL PRSYMB(LUPRI,'=',80,0)
C
C     *****************************
C     *** Print linear response ***
C     *****************************
C     
      CALL HEADER('Linear response',0)
C
      WRITE(LUPRI,'(A,/A)')
     &     ' Boson symmetry    A-operator    B-operator    '//
     &     'Frequency [a.u]    Polarizability',
     &     ' ----------------------------------------------'//
     &     '---------------------------------'
C
C     Loop over A- and B-operators
C     ============================
C     Check that neither A- nor B-operator is excitation operator and
C     check if symmetry allowed.
C
      DO IBOP=1,NQROP
         IBPT = LQROP(IBOP,1)
         IF (IBPT.GE.IEXCOFF) GOTO 220
         DO IAOP=IBOP,NQROP
            IAPT  = LQROP(IAOP,1)
            ALAB  = PRPNAM(IAPT)(1:8)
            ISYMA = IPRPSYM(IAPT)
            BLAB  = PRPNAM(IBPT)(1:8)
            ISYMB = IPRPSYM(IBPT)
            IF ((IAPT.GE.IEXCOFF).OR.
     &          (MULD2H(ISYMA,ISYMB).NE.1)) GOTO 210
C
C     Loop over frequencies and print
C     ===============================
C
            DO IFR=1,LQROP(IBOP,2)
               FREQ = QRFREQ(IFR,IBOP)
               RSP  = RES_LR(1,IAOP,IBOP,IFR)
C
C     Polarizability is the negative value of the response.
C
               WRITE(LUPRI,
     &               '(7X,A3,10X,A8,6X,A8,6X,F14.6,4X,F14.6)')
     &               REP(ISYMA-1),ALAB,BLAB,FREQ,-RSP
            END DO
 210        CONTINUE
         END DO
 220     CONTINUE
      END DO
      CALL PRSYMB(LUPRI,'*',80,0)
C
C     *************************
C     *** Print TPA results ***
C     *************************
C
      CALL HEADER('Quadratic response properties',-1)
C
      IQRF = 0
      NFC  = NBSYM/NFSYM
C
C     Loop over excitation symmetries
C     ===============================
C
      DO IFERE = 1,NFSYM
      DO IS = 1,NFC
         JOPSY  = IFERE
         JSYMEX = JFSYM(IS,JOPSY)
         INDPRP = IEXCOFF + JSYMEX
C
         CALL PRSYMB(LUPRI,'-',80,0)
         WRITE(LUPRI,'(/A,I2,A,A3/)')
     &        ' *** Excited state boson symmetry',
     &        JSYMEX,' : ',REP(JSYMEX-1)
         WRITE(LUPRI,'(/2X,A,I2)')
     &        'No. of excitations in the symmetry: ',
     &        TPACNV(JSYMEX)
C
C     Loop over frequencies
C     =====================
C
         DO IFR = 1,TPACNV(JSYMEX)
            TPDF = 0
            TPDG = 0
            FREQ = EXCIT(IFR,JSYMEX)
C
            CALL PRSYMB(LUPRI,'.',80,0)
            WRITE(LUPRI,'(/3X,A,I2,A,F14.6,A)')
     &           '* Excitation energy no. ',IFR,' : ',
     &           FREQ, ' a.u.'
            WRITE(LUPRI,'(32X,F14.6,A)')
     &           FREQ*XTEV, ' eV'
C
            CALL HEADER('Two-photon transition tensor elements',0)
            WRITE(LUPRI,'(5X,A,/5X,A)')
     &            'A-operator    B-operator    S_{ab} [a.u]',
     &            '------------------------------------------'
C
C     Loop over A- and B-operators
C     ============================
C     Due to symmetry between S_{ab} and S_{ba} only necessary to
C     six out of nine elements in the tensor S.
C
            DO IAOP = 1,NAQROP
            DO IBOP = IAOP,NBQROP
               IAPT  = LAQROP(IAOP)
               ISYMA = IPRPSYM(IAPT)
               ALAB  = PRPNAM(IAPT)(1:8)
               IBPT  = LAQROP(IBOP)
               ISYMB = IPRPSYM(IBPT)
               BLAB  = PRPNAM(IBPT)(1:8)
C
C     Extracting tensor elements and calculating Df and Dg
C     ====================================================
C
               SAA = RES_TPA(1,JSYMEX,IAOP,IAOP,IFR)
               SAB = RES_TPA(1,JSYMEX,IAOP,IBOP,IFR)
               SBB = RES_TPA(1,JSYMEX,IBOP,IBOP,IFR)
               IF (IAOP.EQ.IBOP) THEN
                  TPDF = TPDF+SAA*SBB
                  TPDG = TPDG+SAB*SAB
               ELSE
                  TPDF = TPDF+2*SAA*SBB
                  TPDG = TPDG+2*SAB*SAB
               END IF
C
C     Check if allowed due to symmetry and print
C     ==========================================
C
               IF (MULD2H(ISYMA,MULD2H(ISYMB,JSYMEX)).EQ.1) THEN
                  IQRF = IQRF +1
                  IF (IQRF.LE.MAXQR)
     &                 WRITE(LUPRI,'(6X,A8,6X,A8,5X,F14.6)')
     &                       ALAB,BLAB,SAB
               END IF
            END DO
            END DO
            CALL HEADER('Two-photon transition probabilities',0)
C
C     Calculate and print linear and circular polarization transition moments
C     =======================================================================
C
            DTPLP = (TPDF+2*TPDG)/1.5D1
            DTPCP = (-TPDF+3*TPDG)/1.5D1
C
            PFMT = MXFORM(TPDF,12)
            WRITE(LUPRI,'(5X,A5,'//PFMT//',A)')
     &            'Df  =',TPDF,
     &            ' a.u.   (delta_f)'
            PFMT = MXFORM(TPDF,12)
            WRITE(LUPRI,'(5X,A5,'//PFMT//',A)')
     &            'Dg  =',TPDG,
     &            ' a.u.   (delta_g)'
            PFMT = MXFORM(DTPLP,12)
            WRITE(LUPRI,'(5X,A5,'//PFMT//',A)')
     &            'Dlp =',DTPLP,
     &            ' a.u.   (Linearly polarized light)'
            PFMT = MXFORM(DTPCP,12)
            WRITE(LUPRI,'(5X,A5,'//PFMT//',A)')
     &            'Dcp =',DTPCP,
     &            ' a.u.   (Circularly polarized light)'
            PFMT = MXFORM(DTPCP/DTPLP,12)
            WRITE(LUPRI,'(5X,A5,'//PFMT//',A)')
     &            'R   =',DTPCP/DTPLP,
     &            '        (Polarization ratio, R=Dcp/Dlp)'
         END DO
      END DO
      END DO
      CALL PRSYMB(LUPRI,'*',80,0)
C
      CALL QEXIT('TPA_PRINT')
      RETURN
      END
C
      SUBROUTINE DEF_TPA(WORK,LWRK,DOTPA)
C***********************************************************************
C
C     Input handeling for .TWO-PHOTON
C
C     Written by johhe 2002
!
!  Miro: July 2018, no more need for WORK(:) stuff
!
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dgroup.h"
      PARAMETER ( D1 = 1.0D0 )
      CHARACTER PNAME*16
      DIMENSION WORK(LWRK),PNAME(3)
      LOGICAL DOTPA
C
      CALL QENTER('DEF_TPA')
C
      LFREE = LWRK
      KFREE = 1
C
C     Initialize /XCIBQR/
C     ===================
C
      NAQROP = 0
      NBQROP = 0
      NCQROP = 0
      CALL IZERO(TPACNV,8)
      CALL IZERO(TPASIM,8)
      CALL IZERO(TPASTV,8)
C
      READ(LUCMD,*,END=999,ERR=999) (TPACNV(I),I=1,NBSYM)
      DO I=1,NBSYM
         IF (TPACNV(I).NE.0) DOTPA = .TRUE.
         IF (TPACNV(I).GT.MAXFQR) THEN
            WRITE(LUPRI,*)
     &           'WARNING: No. of excitations exceeds '//
     &           'maximum no. allowed,', MAXFQR
            CALL QUIT('.TWO-PHOTON: Too many excitations')
         END IF
      END DO
      IF (.NOT. DOTPA) THEN
         WRITE(LUPRI,*)
     &        'WARNING: No excitations given to '//
     &        '.TWO-PHOTON.'
         CALL QUIT('.TWO-PHOTON: No excitations')
      ELSE
         CALL ICOPY(8,TPACNV,1,TPASIM,1)
         CALL ICOPY(8,TPACNV,1,TPASTV,1)
      END IF
      PNAME(1) = 'XDIPLEN'
      PNAME(2) = 'YDIPLEN'
      PNAME(3) = 'ZDIPLEN'
      DO I=1,3
         CALL XPRIND(PNAME(I),1,1,D1,PNAME(I),INDXPR,ISYXPR,
     &               ITRXPR,IPRXQR)
         CALL OP1IND('NAQROP',IND1OP,LAQROP,NAQROP,INDXPR,MAXLQR)
         CALL OP1IND('NBQROP',IND1OP,LBQROP,NBQROP,INDXPR,MAXLQR)
      END DO
C
      CALL QEXIT('DEF_TPA')
      RETURN
C
 999  CONTINUE
      CALL QUIT('Input error: '//
     &          'Incorrect input following .TWO-PHOTON.')
C
      END
C
      SUBROUTINE EXCINIT(INDEX,NFRQ,FREQ,JSYM,IFER,ITIM,LAB,
     &                   NZYEE,NZYEP,VECEE,VECEP,WORK,LWRK)
C***********************************************************************
C
C     Initialize variables and read excitation vectors from file
C     for given TPA calculation.
C
C     Written by johhe 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dgroup.h"
#include "dcbxrs.h"
C
      PARAMETER ( D0=0.0D0 )
      CHARACTER LAB*8
      DIMENSION VECEE(MZYEE,NZ),VECEP(MZYEP,NZ)
      DIMENSION WORK(*)
C
      CALL QENTER('EXCINIT')
C
      KFREE = 1
      LFREE = LWRK
C
C     Open file with response vectors
C     ===============================
C
      IF (.NOT. XQR_SKIPEE) 
     &     OPEN(LUQREE,FILE='QRVEC.EE',FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*MZYEE*NZ,STATUS='OLD')
      IF (.NOT. XQR_SKIPEP) 
     &     OPEN(LUQREP,FILE='QRVEC.EP',FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*MZYEP*NZ,STATUS='OLD')
C
C     Set variables
C     =============
C
      FREQ = EXCIT(NFRQ,JSYM)
      ITIM = 1
      LAB  = 'EXCLABL '
C
      JOPSY  = IFER
      JSYMOP = JSYM
      JTIMOP = ITIM
      CALL XRSPAR(XQR_INDSTR,XQR_SKIPEE,XQR_SKIPEP,
     &            IPRXQR)
      NZYEE = 2*NZXOPE
      NZYEP = 2*NZXOPP
      CALL INDQR(INDEX,FREQ,IREC)
      IF (.NOT. XQR_SKIPEE) 
     &     CALL READAC(LUQREE,MZYEE*NZ,VECEE,IREC)
      IF (.NOT. XQR_SKIPEP) 
     &     CALL READAC(LUQREP,MZYEP*NZ,VECEP,IREC)
C
C     Print excitation vectors
C     ========================
C
      JSYMOP = JSYM
C
      IF (IPRXQR.GE.2) THEN
         CALL QMATNORM('Response (e-e) vector in EXCINIT',
     &        VECEE,NZYEE,1,NZ)
         CALL QMATNORM('Response (e-p) vector in EXCINIT',
     &        VECEP,NZYEP,1,NZ)
      END IF
      IF (IPRXQR.GE.10) THEN
         IF (.NOT. XQR_SKIPEE) THEN
            CALL HEADER('VECEE '//LAB//' in EXCINIT',-1)
            CALL PRBVEC(LUPRI,VECEE,1,NZYEE)
         ENDIF
         IF (.NOT. XQR_SKIPEP) THEN
            CALL HEADER('VECEP '//LAB//' in EXCINIT',-1)
            CALL PRBVEC(LUPRI,VECEP,1,NZYEP)
         ENDIF
      END IF
C
C     Close files
C     ===========
C
      IF (.NOT. XQR_SKIPEE) CLOSE(LUQREE,STATUS='KEEP')
      IF (.NOT. XQR_SKIPEP) CLOSE(LUQREP,STATUS='KEEP')
C
      CALL QEXIT('EXCINIT')
C
      RETURN
      END
