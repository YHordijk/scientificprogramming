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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Hermem */
      SUBROUTINE MEMHER(DIRFCK,IPRINT)
C***********************************************************************
C
C     Return memory needed for the specified batch indices
C
C     Memory expression of type
C       MWTOT = MEMCON + MEMAB + MEMCD + MABCD
C
C     Written by T.Saue Apr 6 1997
CTROND check for contraction and spherical
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
C
#include "twosta.h"
#include "ccom.h"
#include "symmet.h"
#include "blocks.h"
#include "twocom.h"
C
C
      LOGICAL DIRFCK,SOFOCK,DPATH1,DPATH2
C for now we assume ITYPE = 0
      MAXDER = 0
      MULE   = 0
      NOPREP = 0
      SOFOCK = .NOT.DIRFCK
      NINTYP = 1
      NHCINT = 1
      DPATH1 = .TRUE.
      DPATH2 = .FALSE.
      IDERIV = MAXDER
C
C     Initialization
C
      KFREE  = 1
      MWTOT  = 0
      MEMCON = 0
      MEMAB  = 0
      MEMCD  = 0
      MABCD  = 0
C
C     General information
C     ===================
C
      IF (MAXDER .EQ. 0) THEN
         NCFTYP = 1
      ELSE IF (MAXDER .EQ. 1) THEN
         NCFTYP = 3
      ELSE IF (MAXDER .EQ. 2) THEN
         NCFTYP = 6
      END IF
      JTOP  = 4*(NHTYP - 1) + MAXDER
      JTOP3 = (JTOP + 1)**3
      NRTOP = (JTOP + 1)*(JTOP + 2)*(JTOP + 3)/6
C
C
C     Information on index A
C     ======================
C
      MULA   = ISTBSH(ISHELA)
      NHKTA  = NHKTSH(ISHELA)
      KHKTA  = KHKTSH(ISHELA)
      NUCA   = NUCOSH(ISHELA)
      NORBA  = NORBSH(ISHELA)
      SPHRA  = SPHRSH(ISHELA)
      JMAXA  = NHKTA - 1 + MAXDER
C
C     Information on index B
C
      MULB   = ISTBSH(ISHELB)
      NHKTB  = NHKTSH(ISHELB)
      KHKTB  = KHKTSH(ISHELB)
      NUCB   = NUCOSH(ISHELB)
      NORBB  = NORBSH(ISHELB)
      SPHRB  = SPHRSH(ISHELB)
      JMAXB  = NHKTB - 1 + MAXDER
C
C     Information on index C
C
      MULC   = ISTBSH(ISHELC)
      NHKTC  = NHKTSH(ISHELC)
      KHKTC  = KHKTSH(ISHELC)
      NUCC   = NUCOSH(ISHELC)
      NORBC  = NORBSH(ISHELC)
      SPHRC  = SPHRSH(ISHELC)
      JMAXC  = NHKTC - 1 + MAXDER
C
C     Information on index D
C
      MULD   = ISTBSH(ISHELD)
      NHKTD  = NHKTSH(ISHELD)
      KHKTD  = KHKTSH(ISHELD)
      NUCD   = NUCOSH(ISHELD)
      NORBD  = NORBSH(ISHELD)
      SPHRD  = SPHRSH(ISHELD)
      JMAXD  = NHKTD - 1 + MAXDER
C
C     Information on electron 1
C     =========================
C
      NODCAB = NODSYM(MAXOPR,MULA,MULB)
      GENAB  = .NOT.(SEGMSH(ISHELA) .AND. SEGMSH(ISHELB))
      IGENAB = 1
      IF (.NOT.GENAB) IGENAB = 2
      DIAGAB = ISHELA .EQ. ISHELB
      SPHRAB = SPHRA .OR. SPHRB
      DIACAB = DIAGAB .AND. .NOT.SPHRAB
      MAXAB  = NHKTA + NHKTB - 2
      TCONAB = DIAGAB .AND. MAXAB .EQ. 0
      KHKTAB = KHKTA*KHKTB
      KCKTAB = KCKTA*KCKTB
      IF (DIAGAB) KHKTAB = KHKTA*(KHKTA + 1)/2
      IF (DIACAB) KCKTAB = KCKTA*(KCKTA + 1)/2
      IF(TCONAB) THEN
        NORBAB = NORBA*(NORBA+1)/2
        MXUCAB = NUCA*(NUCA+1)/2
        NUCAB  = MXUCAB
      ELSE
        NORBAB = NORBA*NORBB
        MXUCAB = NUCA*NUCB
        NUCAB  = MXUCAB
      ENDIF
C     In uncontracted basis NSETA is the number of exponents
C     in the block(IGENAB = 2)
      NSETA  = NSETSH(ISHELA,IGENAB)
      NSETB  = NSETSH(ISHELB,IGENAB)
C
C     Information on electron 2
C     =========================
C
      NODCCD = NODSYM(MAXOPR,MULC,MULD)
      GENCD  = .NOT.(SEGMSH(ISHELC) .AND. SEGMSH(ISHELD))
      IGENCD = 1
      IF (.NOT.GENCD) IGENCD = 2
      DIAGCD = ISHELC .EQ. ISHELD
      SPHRCD = SPHRC .OR. SPHRD
      DIACCD = DIAGCD .AND. .NOT.SPHRCD
      MAXCD  = NHKTC + NHKTD - 2
      TCONCD = DIAGCD .AND. MAXCD .EQ. 0
      KHKTCD = KHKTC*KHKTD
      KCKTCD = KCKTC*KCKTD
      IF (DIAGCD) KHKTCD = KHKTC*(KHKTC + 1)/2
      IF (DIACCD) KCKTCD = KCKTC*(KCKTC + 1)/2
      IF(TCONCD) THEN
        NORBCD = NORBC*(NORBC+1)/2
        MXUCCD = NUCC*(NUCC+1)/2
        NUCCD  = MXUCCD
      ELSE
        NORBCD = NORBC*NORBD
        MXUCCD = NUCC*NUCD
        NUCCD  = MXUCCD
      ENDIF
C     In uncontracted basis NSETA is the number of exponents
C     in the block(IGENAB = 2)
      NSETC  = NSETSH(ISHELC,IGENCD)
      NSETD  = NSETSH(ISHELD,IGENCD)
C
C     Information on both pairs
C     =========================
C
      JMAX0  = MAXAB + MAXCD
      JMAX   = JMAX0 + MAXDER
      NRTUV  = (JMAX + 1)**3
      NTUV   = (JMAX + 1)*(JMAX + 2)*(JMAX + 3)/6
      KCKMAX = MAX(KCKTA,KCKTB,KCKTC,KCKTD)
      NOABCD = NORBAB*NORBCD
      NCCINT = NOABCD*KHKTAB*KHKTCD
      NUABCD = NUCAB*NUCCD
      LCOFAB = MXUCAB*(JMAXA+JMAXB+1)*(JMAXA+1)*(JMAXB+1)*3*NCFTYP
      LCOFCD = MXUCCD*(JMAXA+JMAXB+1)*(JMAXA+1)*(JMAXB+1)*3*NCFTYP
      IF (DIRFCK) THEN
         NINTS  = 0
         NINTMX = 0
      ELSE
         CALL NINTSO(MULE,.FALSE.,.FALSE.,.TRUE.,SOFOCK,.FALSE.,.FALSE.,
     &               IPRINT)
      END IF
      NSOINT = NOABCD*NINTS
      LAOINT = NCCINT*NINTYP
      LRJ000 = (JMAX + 1)*NUABCD
      IF (JMAX .EQ. 0) THEN
         LHRWRK = 0
      ELSE
         LHRWRK = NTUV*NUABCD
      ENDIF
C
C     Allocations in TWOCAL
C
      KFRSAV  = KFREE
      CALL MEMGT2('INTE',KINDHR,JTOP3  ,KFREE)
      CALL MEMGT2('INTE',KINDSQ,NRTOP  ,KFREE)
      CALL MEMGT2('INTE',KIODHR,8*NRTOP,KFREE)
      MWHRND = KFREE - KFRSAV
      MEMCON = MEMCON + MWHRND
C
C     Allocations in TWOODS
C
      LUUAB = NUCA*NUCB*NODCAB
      LUUCD = NUCC*NUCD*NODCCD
      IF (GENAB) THEN
         LCONTA = 2*NORBA*NUCA*NODCAB
         LCONTB = 2*NORBB*NUCB*NODCAB
         LPNTAB = 2*LUUAB
         LREDAB = NORBA*NORBB
         LNCSAB = 0
      ELSE
         LCONTA = 0
         LCONTB = 0
         LPNTAB = 0
         LREDAB = 0
         LNCSAB = NORBA*NORBB*NODCAB
      END IF
      IF (GENCD) THEN
         LCONTC = 2*NORBC*NUCC*NODCCD
         LCONTD = 2*NORBD*NUCD*NODCCD
         LPNTCD = 2*LUUCD
         LREDCD = NORBC*NORBD
         LNCSCD = 0
      ELSE
         LCONTC = 0
         LCONTD = 0
         LPNTCD = 0
         LREDCD = 0
         LNCSCD = 0
         LNCSCD = NORBC*NORBD*NODCCD
      END IF
      CALL MEMGT2('INTE',KNPCOA,2*NSETA*(NODCAB + 1),KFREE)
      CALL MEMGT2('INTE',KNPCOB,2*NSETB*(NODCAB + 1),KFREE)
      CALL MEMGT2('INTE',KNPCOC,2*NSETC*(NODCCD + 1),KFREE)
      CALL MEMGT2('INTE',KNPCOD,2*NSETD*(NODCCD + 1),KFREE)
      CALL MEMGT2('INTE',KJSTRA,NSETA               ,KFREE)
      CALL MEMGT2('INTE',KJSTRB,NSETB               ,KFREE)
      CALL MEMGT2('INTE',KJSTRC,NSETC               ,KFREE)
      CALL MEMGT2('INTE',KJSTRD,NSETD               ,KFREE)
C
      CALL MEMGT2('REAL',KCORAB,9*LUUAB  ,KFREE)
      CALL MEMGT2('REAL',KCORCD,9*LUUCD  ,KFREE)
      CALL MEMGT2('REAL',KEXPAB,3*LUUAB  ,KFREE)
      CALL MEMGT2('REAL',KEXPCD,3*LUUCD  ,KFREE)
      CALL MEMGT2('REAL',KFACAB,  LUUAB  ,KFREE)
      CALL MEMGT2('REAL',KFACCD,  LUUCD  ,KFREE)
      CALL MEMGT2('INTE',KLMNVL,20*KCKMAX,KFREE)
C
      CALL MEMGT2('REAL',KCONTA,LCONTA,KFREE)
      CALL MEMGT2('REAL',KCONTB,LCONTB,KFREE)
      CALL MEMGT2('INTE',KPNTAB,LPNTAB,KFREE)
      CALL MEMGT2('INTE',KREDAB,LREDAB,KFREE)
      CALL MEMGT2('INTE',KNCSAB,LNCSAB,KFREE)
C
      CALL MEMGT2('REAL',KCONTC,LCONTC,KFREE)
      CALL MEMGT2('REAL',KCONTD,LCONTD,KFREE)
      CALL MEMGT2('INTE',KPNTCD,LPNTCD,KFREE)
      CALL MEMGT2('INTE',KREDCD,LREDCD,KFREE)
      CALL MEMGT2('INTE',KNCSCD,LNCSCD,KFREE)
C
      CALL MEMGT2('INTE',KINDAB,2*NORBA*NORBB,KFREE)
      CALL MEMGT2('INTE',KINDCD,2*NORBC*NORBD,KFREE)
C
C     Allocations in TWOOD1
C
      CALL MEMGT2('INTE',KCORBA,NORBA,KFREE)
      CALL MEMGT2('INTE',KCORBB,NORBB,KFREE)
      CALL MEMGT2('INTE',KCORBC,NORBC,KFREE)
      CALL MEMGT2('INTE',KCORBD,NORBD,KFREE)
C
      CALL MEMGT2('REAL',KCOFAB,LCOFAB,KFREE)
      CALL MEMGT2('REAL',KCOFCD,LCOFCD,KFREE)
C
      CALL MEMGT2('REAL',KSOINT,NSOINT         ,KFREE)
      CALL MEMGT2('INTE',KPNTAO,  NINTMX*NOPREP,KFREE)
      CALL MEMGT2('INTE',KPNTOP,3*NINTMX*NOPREP,KFREE)
      CALL MEMGT2('INTE',KPNTNO,4*NINTMX*NOPREP,KFREE)
      CALL MEMGT2('INTE',KPNTRP,3*NINTMX*NOPREP,KFREE)
      CALL MEMGT2('INTE',KPNTLG,3*NINTMX*NOPREP,KFREE)
C      IF(SOFOCK.OR.DIRFCK) THEN
C        LDNSBF = 2*NDMAT
C      ELSE
C        LDNSBF = 0
C      ENDIF
C      CALL MEMGT2('REAL',KDNSBF,LDNSBF,KFREE)
      MWFCAB = LCOFAB
      MWFCCD = LCOFCD
      MWPSO  = 0
      MWSOIN = NSOINT
C
C     Allocations in SYMLOP
C
      CALL MEMGT2('REAL',KAOINT,LAOINT,KFREE)
      MWAOIN = LAOINT
C
C     Allocations in INTDER
C
      LHRINT = NUABCD*NTUV
      CALL MEMGT2('REAL',KHRINT,LHRINT,KFREE)
      MWHRIN = LHRINT
C
C     ***************************************
C     ********** Hermite Integrals **********
C     ***************************************
C     --> Allocations in HR2DRV
C
      KHH = KFREE
        MWHRIN = 0
        CALL MEMGT2('REAL',KRJ000,LRJ000,KFREE)
        CALL MEMGT2('REAL',KPQX  ,NUABCD,KFREE)
        CALL MEMGT2('REAL',KPQY  ,NUABCD,KFREE)
        CALL MEMGT2('REAL',KPQZ  ,NUABCD,KFREE)
        MWRJ00 = LRJ000
C
C       Allocations in R000
C
        KFRSAV = KFREE
          CALL MEMGT2('REAL',KWVALU,NUABCD,KFREE)
          CALL MEMGT2('REAL',KALPHA,NUABCD,KFREE)
          CALL MEMGT2('REAL',KALPHJ,NUABCD,KFREE)
          CALL MEMGT2('INTE',KINDAD,NUABCD,KFREE)
          MWTOT = MAX(MWTOT,KFREE)
        KFREE = KFRSAV
C
C       Allocations in HERI
C
        KFRSAV = KFREE
          CALL MEMGT2('REAL',KHRWRK,LHRWRK,KFREE)
          MWHRSQ = KFREE - KFRSAV
          MWTOT = MAX(MWTOT,KFREE)
        KFREE = KFRSAV
      KFREE = KHH
C
C     *****************************************
C     ********** Cartesian Integrals **********
C     *****************************************
C..assuming DPATH1  
C     For contracted integrals define NUCR1, NORR1 etc. !!!!
C  hjaaj: NPCO[ABCD] not defined yet !!!
C         use instead NUC[ABCD] and NORB[ABCD] / July 2000
        IF (GENAB) THEN
c           NUCR1  = ISUM(NSETA,NPCOA(1,1),1)
c           NUCR2  = ISUM(NSETB,NPCOB(1,1),1)
c           NORR1  = ISUM(NSETA,NPCOA(1,2),1)
c           NORR2  = ISUM(NSETB,NPCOB(1,2),1)
            NUCR1  = NUCA
            NUCR2  = NUCB
            NORR1  = NORBA
            NORR2  = NORBB
        END IF
        IF (GENCD) THEN
c           NUCR3  = ISUM(NSETC,NPCOC(1,1),1)
c           NUCR4  = ISUM(NSETD,NPCOD(1,1),1)
c           NORR3  = ISUM(NSETC,NPCOC(1,2),1)
c           NORR4  = ISUM(NSETD,NPCOD(1,2),1)
            NUCR3  = NUCC
            NUCR4  = NUCD
            NORR3  = NORBC
            NORR4  = NORBD
        END IF
      KCC = KFREE    
        NCCPP = NORBAB*NUCCD
        NTUVCD = (MAXCD + 1)*(MAXCD + 2)*(MAXCD + 3)/6
        LHCINT = NORBAB*NUCCD*NTUVCD*KHKTAB*NHCINT
        CALL MEMGT2('INTE',KODD12,2*KCKTAB,KFREE)
        CALL MEMGT2('INTE',KODD34,2*KCKTCD,KFREE)
        CALL MEMGT2('REAL',KHCINT,LHCINT  ,KFREE)
        MWHCIN = LHCINT
C
C       Allocations in C1DRIV
C         DC101  = .TRUE.
C
        KFRSAV = KFREE
          IF (GENAB) THEN
            LSCR1 = NUCR1*NORR2*NUCCD
            LSCR2 = NUCR1*NUCR2*NUCCD
          ELSE
            LSCR1 = 0
            LSCR2 = 0
          END IF
          IF (KHKTAB .EQ. 1) THEN
            LETUV  = 0
            LHCPRM = 0
          ELSE
            LETUV  = NUABCD
            LHCPRM = NUABCD*NTUVCD
          END IF
          LCCONT = 0
          LCSINT = 0
          IF (SPHRAB)            LCCONT = NCCPP*NTUVCD*KCKTAB
          IF (SPHRA .AND. SPHRB) LCSINT = NCCPP*NTUVCD*KHKTAB
          LNGTH1 = LHCPRM + LETUV
          LNGTH2 = LHCPRM + LSCR1 + LSCR2
          LNGTH3 = LCSINT
          LBUF   = MAX(LNGTH1,LNGTH2,LNGTH3)
          CALL MEMGT2('REAL',KCCONT,LCCONT,KFREE)
          CALL MEMGT2('REAL',KBUF,LBUF,KFREE)
          MWC1DR = KFREE - KFRSAV
          MWTOT = MAX(MWTOT,KFREE)
        KFREE = KFRSAV
C
C       Allocations in C2DRIV
C
C         DC2H = DC2H1  = .TRUE.
C
        KFRSAV = KFREE
          NHCMAX = KHKTAB
          LNHCC  =  2*KCKTCD
          CALL MEMGT2('INTE',KNHCC,LNHCC,KFREE)
          MWC2DR = KFREE - KFRSAV
          KBF2 = KFREE
          LSSINT = NOABCD*NHCMAX*KHKTCD
          LCCPRM = NCCPP*NHCMAX
          IF (KHKTCD .EQ. 1) THEN
            LETUV = 0
          ELSE
            LETUV  = NCCPP
          END IF
          IF (SPHRC .AND. SPHRD) THEN
            LCSINT = NOABCD*NHCMAX*KHKTD
          ELSE
            LCSINT = 0
          END IF
          IF (SPHRCD) THEN
            LCCONT = NOABCD*NHCMAX*KCKTCD
          ELSE
            LCCONT = 0
          END IF
          IF (GENCD) THEN
            LSCR1= NUCR3*NUCR4*NHCMAX*NORBAB
            LSCR2= NORR3*NUCR4*NHCMAX*NORBAB
          ELSE
            LSCR1 = 0
            LSCR2 = 0
          END IF
C
          LNGTH1 = LCCPRM + LETUV
          LNGTH2 = LCCPRM + LSCR1 + LSCR2
          LNGTH3 = LCSINT
          LBUF   = MAX(LNGTH1,LNGTH2,LNGTH3)
          CALL MEMGT2('REAL',KSSINT,LSSINT,KFREE)
          CALL MEMGT2('REAL',KCCONT,LCCONT,KFREE)
          CALL MEMGT2('REAL',KBUF  ,LBUF  ,KFREE)
          MWC2HI = KFREE - KBF2
          MWTOT = MAX(MWTOT,KFREE)
        KFREE = KFRSAV
C
C       Print section
C
      IF(IPRINT.GE.1) THEN
         CALL HEADER('Maximum work space allocations in MEMHER',1)
         WRITE (LUPRI,'(1X,A,4I4,15(/1X,A,I12))')
     &      ' Shell combination :                 ', ISHELA,ISHELB,
     &                                               ISHELC,ISHELD,
     &      ' Total allocation:                   ', MWTOT,
     &      ' Expansion coefficients (electron 1):', MWFCAB,
     &      ' Expansion coefficients (electron 2):', MWFCCD,
     &      ' SO integrals:                       ', MWSOIN,
     &      ' AO integrals:                       ', MWAOIN,
     &      ' Hermite integrals:                  ', MWHRIN,
     &      ' Hermite integral pointer:           ', MWHRND,
     &      ' Hermite-Cartesian integrals:        ', MWHCIN,
     &      ' Gamma functions:                    ', MWRJ00,
     &      ' Work space for Hermite integrals:   ', MWHRSQ,
     &      ' Allocations in C1DRIV:              ', MWC1DR,
     &      ' Allocations in C2DRIV:              ', MWC2DR,
     &      ' Allocations in C2HINT:              ', MWC2HI
      ENDIF

C
      RETURN
      END
C  /* Deck memgt2 */
      SUBROUTINE MEMGT2(TYPE,KBASE,LENGTH,KFREE)
C
C Copyright 9-Jan-1988/27-Jul-1993 Hans Joergen Aa. Jensen
C
C     Memory allocation
C
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) TYPE
C
#include "iratdef.h"
C
      PARAMETER ( LENID = 1/IRAT + 1)
C
      CALL QENTER('MEMGT2')
      IF (KFREE .EQ. 1) THEN
         LENIDS = 2*LENID
      ELSE
         LENIDS = LENID
      END IF
      IF (TYPE(1:4) .EQ. 'REAL') THEN
         LREAL = LENGTH
      ELSE IF (TYPE(1:4) .EQ. 'INTE') THEN
         LREAL = (LENGTH-1)/IRAT + 1
      ELSE IF (TYPE(1:4) .EQ. 'LOGI') THEN
         LREAL = (LENGTH-1)/LRAT + 1
      ELSE IF (TYPE(1:4) .EQ. 'INT2') THEN
         LREAL = (LENGTH-1)/IRAT2 + 1
      ELSE
         WRITE (LUPRI,'(/2A)')
     *      ' MEMGT2 ERROR, illegal type :',TYPE
         CALL QUIT('MEMGT2 ERROR, illegal type')
      END IF
C
      IF (KFREE .EQ. 1) THEN
         KFREE = KFREE + LENID
      ENDIF
      KBASE  = KFREE
      KFREE  = KBASE + LREAL + LENID
      CALL QEXIT('MEMGT2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Memlop */
      SUBROUTINE MEMLOP(I2TYP,DIRFCK)
C***********************************************************************
C
C  Print out integral types and count the use of the memmory.
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "mxcent.h"
C
#include "blocks.h"
#include "twosta.h"
#include "twocom.h"
#include "nuclei.h"
      CHARACTER SPDCAR*1,COMP(2)*1
      LOGICAL DIRFCK
      DIMENSION IMSHL(4)
C
      IPRINT = 1
      COMP(1) = 'L'
      COMP(2) = 'S'
      IF(I2TYP.EQ.0) THEN
        IASTRT = 1
        IBSTRT = 1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
CMI    ... LL-integrals
        IASTRT = 1
        IBSTRT = 1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = NLRGBL
        IBSMAX = NLRGBL
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
CMI   ...  SL-integrals
        IASTRT = NLRGBL+1
        IBSTRT = NLRGBL+1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.3) THEN
CMI    ...  SS-integrals
        IASTRT = NLRGBL+1
        IBSTRT = NLRGBL+1
        ICSTRT = NLRGBL+1
        IDSTRT = NLRGBL+1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSE
        WRITE(LUPRI,'(A,I5)') 'MEMLOP: Unknown I2TYP =' ,I2TYP
        CALL QUIT('Unknown I2TYP !!!')
      ENDIF
C
      MWMAX = 0
      DO ISHELA = IASTRT,IASMAX
        NCENTA = NCNTSH(ISHELA)
        DO ISHELB = IBSTRT,ISHELA
          NCENTB = NCNTSH(ISHELB)
          ICMAX = ISHELA
          IF(I2TYP.EQ.2) ICMAX = NLRGBL
          DO ISHELC = ICSTRT, ICMAX
            NCENTC = NCNTSH(ISHELC)
            IDMAX = ISHELC
            IF (ISHELA.EQ.ISHELC) IDMAX = ISHELB
            DO ISHELD = IDSTRT,IDMAX
              NCENTD = NCNTSH(ISHELD)
              CALL MEMHER(DIRFCK,0)
              IF(MWTOT.GT.MWMAX) THEN
                IMSHL(1) = ISHELA
                IMSHL(2) = ISHELB
                IMSHL(3) = ISHELC
                IMSHL(4) = ISHELD
                MWMAX    = MWTOT
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     Print section
C
      ISHELA = IMSHL(1)
      ISHELB = IMSHL(2)
      ISHELC = IMSHL(3)
      ISHELD = IMSHL(4)
      NCENTA = NCNTSH(ISHELA)   
      NCENTB = NCNTSH(ISHELB)   
      NCENTC = NCNTSH(ISHELC)   
      NCENTD = NCNTSH(ISHELD)   
      ICA    = LCLASH(ISHELA)
      ICB    = LCLASH(ISHELB)
      ICC    = LCLASH(ISHELC)
      ICD    = LCLASH(ISHELD)
      WRITE(LUPRI,'(4(I4,A1,A4,1X,3A1),I20)')
     &        ISHELA,'(',NAMN(NCENTA),COMP(ICA),
     &        SPDCAR(NHKTSH(ISHELA)-1),')',
     &        ISHELB,'(',NAMN(NCENTB),COMP(ICB),
     &        SPDCAR(NHKTSH(ISHELB)-1),')',
     &        ISHELC,'(',NAMN(NCENTC),COMP(ICC),
     &        SPDCAR(NHKTSH(ISHELC)-1),')',
     &        ISHELD,'(',NAMN(NCENTD),COMP(ICD),
     &        SPDCAR(NHKTSH(ISHELD)-1),')',
     &        MWMAX
C
      MXHERMEM(I2TYP) = MWMAX
C
      RETURN
      END

