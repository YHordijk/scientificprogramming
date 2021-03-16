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
C  /* Deck un2out */
      SUBROUTINE UN2OUT(SO,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,THRESH,
     &                  NINDAB,NINDCD,IPRINT)
C
C     Write out blocks of symmetry integrals, eliminating duplicates
C
C                                          880412   PRT
C
C     Some low-brain work has been done by TUH
C
C     Rewritten to allow for triangular looping and to eliminate
C     all duplicates by testing  880601 TUH
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
      PARAMETER (LBUF = 600)
      LOGICAL DCMPAB, DCMPCD, DCMPAC, DRALTB, DRCLTD, FIRST, LAST,
     &        DRABAB, DCABAB, IAEQIC, IALTIC, IPNTLG(3,*), NOTEST
      DIMENSION SO(*), BUF(LBUF), IBUF(LBUF,2),
     &          IPNTNO(4,*), IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2)
#include "nuclei.h"
#include "twocom.h"
#include "symmet.h"
      SAVE BUF, IBUF, ICOUNT
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine UN2OUT',-1)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(2X,A,4I5)') 'NHKT? ', NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI,'(2X,A,4I5)') 'MUL?  ', MULA,  MULB,  MULC,  MULD
         WRITE (LUPRI,'(2X,A,4I5)') 'NORB? ', NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI,'(2X,A,4I5)') 'NSTR? ', NSTRA, NSTRB, NSTRC, NSTRD
         WRITE (LUPRI,'(2X,A,2I5)') 'NORBCD', NORBCD
         WRITE (LUPRI,'(2X,A,2I5)') 'NOABCD', NOABCD
         WRITE (LUPRI,'(2X,A,2L5)') 'DIAGAB/CD', DIAGAB, DIAGCD
         WRITE (LUPRI,'(2X,A,2L5)') 'TCONAB/CD', TCONAB, TCONCD
         WRITE (LUPRI,'(2X,A,2L5)') 'SHAEQB/CD', SHAEQB, SHCEQD
         WRITE (LUPRI,'(2X,A, L5)') 'SHABAB', SHABAB
      END IF
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (INT_STAR8)
      NIBUF = 1
      IF (NBASIS .LE. 255) THEN
         NBITS = 8
      ELSE
         NBITS = 16
      END IF
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
C     *******************************************************
C     ***** Initialization when subroutine first called *****
C     *******************************************************
C
      IF (FIRST) CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,-1,NBITS,0,
     &                       IPRINT)
C
      ISOFF  = 0
      NBUFCL = 0
      NSTART = ICOUNT
      NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
      DO 100 I = 1, NINTS
         NSTRNA = IPNTNO(1,I)
         NSTRNB = IPNTNO(2,I)
         NSTRNC = IPNTNO(3,I)
         NSTRND = IPNTNO(4,I)
         IREPA  = IPNTRP(1,I)
         IREPB  = IPNTRP(2,I)
         IREPC  = IPNTRP(3,I)
         IREPD  = IEOR(IEOR(IREPA,IREPB),IREPC)
         IF (NOTEST) THEN
            IF (NIBUF .EQ. 1) THEN
               INT = 0
               DO 200 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
                  DO 210 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT = SO(ISOFF+INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1)=MAX(INDAB,INDCD)*IBIT2+INDAB+INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
                          CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                       END IF
                     END IF
  210             CONTINUE
  200          CONTINUE
            ELSE
               INT = 0
               DO 205 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
                  DO 215 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT = SO(ISOFF+INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1) = MAX(INDAB,INDCD)
                       IBUF(ICOUNT,2) = MIN(INDAB,INDCD)
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
                          CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                       END IF
                     END IF
  215             CONTINUE
  205          CONTINUE
            END IF
         ELSE
            DCMPAB = IPNTLG(1,I)
            DCMPCD = IPNTLG(2,I)
            DCABAB = IPNTLG(3,I)
            DRALTB = IREPA .LT. IREPB
            DRCLTD = IREPC .LT. IREPD
            DRABAB = DCABAB .AND. IREPA.EQ.IREPC .AND. IREPB.EQ.IREPD
            INT = 0
            DO 300 IAB = 1, NORBAB
               IA = KHKTA*(NINDAB(IAB,1) - 1)
               IB = KHKTB*(NINDAB(IAB,2) - 1)
               IF (DCMPAB) THEN
                  IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                     INT = INT + NORBCD
                     GO TO 300
                  END IF
               END IF
               INDA = IPTSYM(NSTRNA + IA,IREPA)
               INDB = IPTSYM(NSTRNB + IB,IREPB)
               INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
               DO 310 ICD = 1,NORBCD
                  IC = KHKTC*(NINDCD(ICD,1) - 1)
                  ID = KHKTD*(NINDCD(ICD,2) - 1)
                  INT = INT + 1
                  IF (DCMPCD ) THEN
                     IF (ID.GT.IC) GO TO 310
                     IF (DRCLTD .AND. ID.EQ.IC) GO TO 310
                  END IF
                  IF (DRABAB) THEN
                     IF (IA.LT.IC.OR.(IA.EQ.IC.AND.IB.LT.ID)) GOTO 310
                  END IF
                  SOINT = SO(ISOFF+INT)
                  IF (ABS(SOINT) .GT. THRESH) THEN
                     IF (NIBUF .EQ. 1) THEN
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1)=MAX(INDAB,INDCD)*IBIT2+INDAB+INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                          NBUFCL = NBUFCL + 1
                       END IF
                     ELSE
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1) = MAX(INDAB,INDCD)
                       IBUF(ICOUNT,2) = MIN(INDAB,INDCD)
                       IF (ICOUNT.EQ.LBUF) THEN
                          CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                          NBUFCL = NBUFCL + 1
                       END IF
                     END IF
                  END IF
  310          CONTINUE
  300       CONTINUE
         END IF
         ISOFF = ISOFF + NOABCD
  100 CONTINUE
      NGINT = LBUF*NBUFCL + ICOUNT - NSTART
      CALL DELSTA(0,NGINT)
C
C     *************************************
C     ***** Last call to empty buffer *****
C     *************************************
C
      IF (LAST) CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,1,NBITS,0,
     &                      IPRINT)
      RETURN
      END
C  /* Deck un2wrt */
      SUBROUTINE UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,ITYPE,NBITS,INDA,
     &                  IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "ibtpar.h"
      DIMENSION BUF(LBUF), IBUF(LBUF,NIBUF)
#include "twosta.h"
#include "hrunit.h"
#include "nuclei.h"
      SAVE NBUF
C
      IF (ITYPE .EQ. -1) THEN
         REWIND LUINTA
         CALL NEWLAB('BASTWOEL',LUINTA,LUPRI)
         ICOUNT = 0
         NBUF = 0
      ELSE IF (ITYPE .EQ. 0) THEN
         IF (INDA .NE. 0) WRITE (LUINTA) INDA
         WRITE (LUINTA) BUF,IBUF,ICOUNT
         NBUF = NBUF + 1
C
         IF (IPRINT .GE. 6) THEN
            WRITE (LUPRI,'(2X,A,I5,A/)')
     &         'Integral buffer #',NBUF,' has been written.'
            IBIT1 = 2**NBITS - 1
            DO 100 INT = 1, ICOUNT
               IF (NIBUF .EQ. 1) THEN
                  I = IAND(ISHFT(IBUF(INT,1),-3*NBITS),IBIT1)
                  J = IAND(ISHFT(IBUF(INT,1),-2*NBITS),IBIT1)
                  K = IAND(ISHFT(IBUF(INT,1),  -NBITS),IBIT1)
                  L = IAND(      IBUF(INT,1),          IBIT1)
               ELSE
                  I = IAND(ISHFT(IBUF(INT,1),-NBITS),IBIT1)
                  J = IAND(      IBUF(INT,1),        IBIT1)
                  K = IAND(ISHFT(IBUF(INT,2),-NBITS),IBIT1)
                  L = IAND(      IBUF(INT,2),        IBIT1)
               END IF
               WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8)')
     &                      ' ## ', I, J, K, L, BUF(INT)
  100       CONTINUE
         END IF
         ICOUNT = 0
      ELSE
         IF (INDA .NE. 0) THEN
            WRITE (LUINTA) -INDA
            WRITE (LUINTA) BUF,IBUF,ICOUNT
         ELSE
            IF (ICOUNT .GT. 0) WRITE (LUINTA) BUF,IBUF,ICOUNT
            WRITE (LUINTA) BUF,IBUF,-1
            REWIND LUINTA
            CALL GPCLOSE(LUINTA,'KEEP')
         END IF
C
         IF (IPRINT .GE. 6) THEN
            IF (ICOUNT .GT. 0) THEN
               WRITE (LUPRI,'(2X,A,I5,A/)')
     &            'Integral buffer #',NBUF+1,' has been written.'
               IBIT1 = 2**NBITS - 1
               DO 200 INT = 1, ICOUNT
                  IF (NIBUF .EQ. 1) THEN
                     I = IAND(ISHFT(IBUF(INT,1),-3*NBITS),IBIT1)
                     J = IAND(ISHFT(IBUF(INT,1),-2*NBITS),IBIT1)
                     K = IAND(ISHFT(IBUF(INT,1),  -NBITS),IBIT1)
                     L = IAND(      IBUF(INT,1),          IBIT1)
                  ELSE
                     I = IAND(ISHFT(IBUF(INT,1),-NBITS),IBIT1)
                     J = IAND(      IBUF(INT,1),        IBIT1)
                     K = IAND(ISHFT(IBUF(INT,2),-NBITS),IBIT1)
                     L = IAND(      IBUF(INT,2),        IBIT1)
                  END IF
                  WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8)')
     &                         ' ## ', I, J, K, L, BUF(INT)
  200          CONTINUE
            END IF
         END IF
C
C        Statistics
C
         IF (INDA .EQ. 0) THEN
            N2WRIT = LBUF*NBUF + ICOUNT
            IF (ICOUNT.GT.0 .AND. INDA.NE.0) THEN
               NBUF = NBUF + 2
            ELSE
               NBUF = NBUF + 1
            END IF
            IF (IRAT .EQ. 1) LWORD = 8
            IF (IRAT .EQ. 2) LWORD = 4
            KBYTES = NBUF*LWORD*(LBUF*IRAT + NIBUF*LBUF + 1)
            KBYTES = (KBYTES + 999)/1000
            NALL   = (NBASIS*(NBASIS + 1))/2
            NALL   = (NALL*(NALL + 1))/2
            PERCNT = dble(100*N2WRIT) / dble(NALL)
            WRITE (LUPRI,'(/1X,A,I10,A,F4.1,A,/1X,A,I10//)')
     &         'Number of two-electron integrals written:',N2WRIT,
     &         ' (',PERCNT,'%)',
     &         'Kilobytes written:                       ',KBYTES
         END IF
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fckout */
      SUBROUTINE FCKOUT(FMAT,DMAT,NDMAT,SO,IPNTNO,IPNTRP,IPNTLG,
     &                  NINDAB,NINDCD,IFCTYP,DINTSKP,IREPDM,
     &                  DMRSO,DNSBUF,HFXFAC,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Direct calculation of Fock matrices in SO-basis
C
C     First routine written by P.Taylor with 
C     low-brain contributions by TUH
C
C     Major surgery by T.Saue October 1996
C       - major speed-up
C       - screening
C       - extended to include all Fock matrices
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
      LOGICAL IPNTLG(3,*)
      DIMENSION SO(*),IPNTNO(4,*), IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          FMAT(*), DMAT(*),DMRSO(*),DINTSKP(*),
     &          WORK(LWORK),DNSBUF(*),IFCTYP(NDMAT),
     &          IREPDM(NDMAT)
#include "twocom.h"
#include "symmet.h"
      CALL QENTER('FCKOUT')
#include "memint.h"
C
      CALL MEMGET('INTE',KIBUF,4*NOABCD,WORK,KFREE,LFREE)
C
      CALL FCKOU1(FMAT,DMAT,NDMAT,SO,WORK(KIBUF),
     &            IPNTNO,IPNTRP,IPNTLG,NINDAB,NINDCD,
     &            DMRSO,DNSBUF,DINTSKP,HFXFAC,IFCTYP,IREPDM,IPRINT)
      CALL MEMREL('FCKOUT',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('FCKOUT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fckou1 */
      SUBROUTINE FCKOU1(FMAT,DMAT,NDMAT,SO,IBUF,IPNTNO,IPNTRP,
     &                  IPNTLG,NINDAB,NINDCD,DMRSO,DNSBUF,
     &                  DINTSKP,HFXFAC,IFCTYP,IREPDM,IPRINT)
C*****************************************************************************
C
C     Direct calculation of Fock matrices in SO-basis
C
C     Revised by T.Saue August 1996
C     Write out blocks of symmetry integrals, eliminating duplicates
C
C                                          880412   PRT
C
C     Some low-brain work has been done by TUH
C
C     Rewritten to allow for triangular looping and to eliminate
C     all duplicates by testing  880601 TUH
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP5 = 0.5D0,DMP25 = -0.25D00,DMP5=-0.5D0,D1 = 1.0D0,
     &           D0 = 0.0D0)
C
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
      INTEGER A,B,C,D
      LOGICAL DCMPAB, DCMPCD, DCMPAC, DRALTB, DRCLTD, 
     &        DRABAB, DCABAB, IAEQIC, IALTIC, IPNTLG(3,*), NOTEST,
     &        LBIT,DOE,DOC,TST2IN
      DIMENSION SO(NOABCD,NINTS), IBUF(4,NOABCD), 
     &          IPNTNO(4,*), IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          FMAT(NBASIS,NBASIS,*), DMAT(NBASIS,NBASIS,*),
     &          DNSBUF(2,NDMAT),DMRSO(MAXSHL,MAXSHL,NDMAT),
     &          IFCTYP(NDMAT),DINTSKP(2,4),IREPDM(NDMAT)
#include "twosta.h"
#include "twocom.h"
#include "symmet.h"
#include "nuclei.h"
#include "blocks.h"
C
      IF (IPRINT .GE. 6) CALL HEADER('Subroutine FCKOU1',-1)
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI,'(2X,A,4I5)') 'NHKT? ', NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI,'(2X,A,4I5)') 'MUL? ', MULA,  MULB,  MULC,  MULD
         WRITE (LUPRI,'(2X,A,4I5)') 'NORB?', NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI,'(2X,A,4I5)') 'NSTR?', NSTRA, NSTRB, NSTRC, NSTRD
         WRITE (LUPRI,'(2X,A,2I5)') 'NORBCD', NORBCD
         WRITE (LUPRI,'(2X,A,2I5)') 'NOABCD', NOABCD
         WRITE (LUPRI,'(2X,A,2L5)') 'DIAGAB/CD', DIAGAB, DIAGCD
         WRITE (LUPRI,'(2X,A,2L5)') 'TCONAB/CD', TCONAB, TCONCD
         WRITE (LUPRI,'(2X,A,2L5)') 'SHAEQB/CD', SHAEQB, SHCEQD
         WRITE (LUPRI,'(2X,A, L5)') 'SHABAB', SHABAB
      END IF
C
C     Initialize screening
C     ====================
C
      DOC = DOSCRN.AND.LBIT(ICEFLG,1)
      DOE = DOSCRN.AND.LBIT(ICEFLG,2)
      DNSMAX   = D0
      IF(DOSCRN) THEN
        DO I = 1,NDMAT 
C Largest Coulomb contribution
          DNSBUF(1,I) = MAX(DMRSO(ISHELA,ISHELB,I),
     &                      DMRSO(ISHELB,ISHELA,I),
     &                      DMRSO(ISHELC,ISHELD,I),
     &                      DMRSO(ISHELD,ISHELC,I))
C Largest exchange contribution
          DNSBUF(2,I) = MAX(DMRSO(ISHELC,ISHELA,I),
     &                      DMRSO(ISHELA,ISHELC,I),
     &                      DMRSO(ISHELC,ISHELB,I),
     &                      DMRSO(ISHELB,ISHELC,I),
     &                      DMRSO(ISHELD,ISHELA,I),
     &                      DMRSO(ISHELA,ISHELD,I),
     &                      DMRSO(ISHELB,ISHELD,I),
     &                      DMRSO(ISHELD,ISHELB,I))
          DNSMAX = MAX(DNSMAX,DNSBUF(1,I),DNSBUF(2,I))
        ENDDO
        NSOINT = NINTS*NOABCD
        IF(DNSMAX.EQ.D0) THEN
          DINTSKP(1,2) = DINTSKP(1,2) + NSOINT
          DINTSKP(2,2) = DINTSKP(2,2) + NSOINT
          RETURN
        ENDIF
        FCKTOL = MAX((SCRTHR/DNSMAX),1.00D-15)
      ELSE
        FCKTOL = 1.00D-15
      ENDIF
      NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
      DO I = 1, NINTS
        NBUF   = 0
        NTOT   = 0
        NSTRNA = IPNTNO(1,I)
        NSTRNB = IPNTNO(2,I)
        NSTRNC = IPNTNO(3,I)
        NSTRND = IPNTNO(4,I)
        IREPA  = IPNTRP(1,I)
        IREPB  = IPNTRP(2,I)
        IREPAB = IEOR(IREPA,IREPB)
        IREPC  = IPNTRP(3,I)
        IREPD  = IEOR(IREPAB,IREPC)
        IREPAC = IEOR(IREPA,IREPC)
        IREPAD = IEOR(IREPA,IREPD)
        IF(NOTEST) THEN
          INT = 0
          SOMAX = D0
          DO IAB = 1, NORBAB
            IA = KHKTA*(NINDAB(IAB,1) - 1)
            IB = KHKTB*(NINDAB(IAB,2) - 1)
            INDA = IPTSYM(NSTRNA + IA,IREPA)
            INDB = IPTSYM(NSTRNB + IB,IREPB)
            DO ICD = 1, NORBCD
              INT = INT + 1
              SOINT = SO(INT,I)
              NTOT = NTOT + 1
              IF (ABS(SOINT) .GT. FCKTOL) THEN
                IC = KHKTC*(NINDCD(ICD,1) - 1)
                ID = KHKTD*(NINDCD(ICD,2) - 1)
                INDC = IPTSYM(NSTRNC + IC,IREPC)
                INDD = IPTSYM(NSTRND + ID,IREPD)
                NBUF = NBUF + 1
                SOMAX = MAX(SOMAX,ABS(SOINT))
                SO(NBUF,I)   = SOINT
                IBUF(1,NBUF) = INDA
                IBUF(2,NBUF) = INDB
                IBUF(3,NBUF) = INDC
                IBUF(4,NBUF) = INDD
              ENDIF
            ENDDO
          ENDDO
        ELSE
C***      DCMPAB = SHAEQB.AND.ICOMPA.EQ.ICOMPB
          DCMPAB = IPNTLG(1,I)
C***      DCMPCD = SHCEQD.AND.ICOMPC.EQ.ICOMPD
          DCMPCD = IPNTLG(2,I)
C***      DCABAB = .NOT.SQ12EL.AND.SHABAB.AND.ICOMPA.EQ.ICOMPC.AND.ICOMPB.EQ.ICOMPD
          DCABAB = IPNTLG(3,I)
          DRALTB = IREPA .LT. IREPB
          DRCLTD = IREPC .LT. IREPD
          DRABAB = DCABAB .AND. IREPA.EQ.IREPC .AND. IREPB.EQ.IREPD
          INT   = 0
          SOMAX = D0
          DO 300 IAB = 1, NORBAB
            FAB = D1
            IA  = KHKTA*(NINDAB(IAB,1) - 1)
            IB  = KHKTB*(NINDAB(IAB,2) - 1)
            INDA = IPTSYM(NSTRNA + IA,IREPA)
            INDB = IPTSYM(NSTRNB + IB,IREPB)
            IF(DCMPAB) THEN
              IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                INT = INT + NORBCD
                GO TO 300
              END IF
              IF(INDB.EQ.INDA) FAB = DP5*FAB
            ENDIF
            IF(DRABAB) THEN
              DO 310 ICD = 1,NORBCD
                INT   = INT + 1
                FABCD = FAB
                IC    = KHKTC*(NINDCD(ICD,1) - 1)
                ID    = KHKTD*(NINDCD(ICD,2) - 1)
                INDC  = IPTSYM(NSTRNC + IC,IREPC)
                INDD  = IPTSYM(NSTRND + ID,IREPD)
                IF (IA.LT.IC.OR.(IA.EQ.IC.AND.IB.LT.ID)) GOTO 310
                IF (INDA.EQ.INDC.AND.INDB.EQ.INDD) FABCD = DP5*FABCD
                IF(DCMPCD) THEN
                  IF (ID.GT.IC) GO TO 310
                  IF (DRCLTD .AND. ID.EQ.IC) GO TO 310
                  IF(INDD.EQ.INDC) FABCD = DP5*FABCD
                ENDIF
                NTOT  = NTOT + 1
                SOINT = SO(INT,I)
                IF (ABS(SOINT) .GT. FCKTOL) THEN
                  NBUF         = NBUF + 1
                  SOMAX        = MAX(SOMAX,ABS(SOINT))
                  SO(NBUF,I)   = FABCD*SOINT
                  IBUF(1,NBUF) = INDA
                  IBUF(2,NBUF) = INDB
                  IBUF(3,NBUF) = INDC
                  IBUF(4,NBUF) = INDD
                END IF
  310         CONTINUE
            ELSE
              DO 320 ICD = 1,NORBCD
                INT   = INT + 1
                FABCD = FAB
                IC    = KHKTC*(NINDCD(ICD,1) - 1)
                ID    = KHKTD*(NINDCD(ICD,2) - 1)
                INDC  = IPTSYM(NSTRNC + IC,IREPC)
                INDD  = IPTSYM(NSTRND + ID,IREPD)
                IF(DCMPCD) THEN
                  IF (ID.GT.IC) GO TO 320
                  IF (DRCLTD .AND. ID.EQ.IC) GO TO 320
                  IF(INDD.EQ.INDC) FABCD = DP5*FABCD
                ENDIF
                NTOT = NTOT + 1
                SOINT = SO(INT,I)
                IF (ABS(SOINT) .GT. FCKTOL) THEN
                  NBUF         = NBUF + 1
                  SO(NBUF,I)   = FABCD*SOINT
                  SOMAX        = MAX(SOMAX,ABS(SOINT))
                  IBUF(1,NBUF) = INDA
                  IBUF(2,NBUF) = INDB
                  IBUF(3,NBUF) = INDC
                  IBUF(4,NBUF) = INDD
                END IF
  320         CONTINUE
            ENDIF
  300     CONTINUE
        ENDIF
        NDEL = NTOT - NBUF
        DINTSKP(1,2) = DINTSKP(1,2) + NTOT
        DINTSKP(2,2) = DINTSKP(2,2) + NDEL
        DINTSKP(1,3) = DINTSKP(1,3) + NCM*NBUF
        DINTSKP(1,4) = DINTSKP(1,4) + NEM*NBUF
C
C       Print integrals
C       ===============
C
        IF(IPRINT.GE.6) THEN
          WRITE(LUPRI,'(A,I5)') '* FCKOU1: Integrals in NINTS =',I
          IF(NBUF.GT.0) WRITE(LUPRI,'(4I5,3X,F12.6)')
     &     (IBUF(1,J),IBUF(2,J),IBUF(3,J),IBUF(4,J),SO(J,I),J=1,NBUF)
        ENDIF
C
C       **********************************
C       ***** Contract Fock matrices *****
C       **********************************
C         IFCTYP = XY
C           X indicates symmetry about diagonal
C             X = 0 No symmetry
C             X = 1 Symmetric
C             X = 2 Anti-symmetric
C             X = 3 Some blocks symmetric and some antisymmetric
C           Y indicates contributions
C             Y = 1 Direct
C             Y = 2 Exchange
C             Y = 3 Direct + Exchange
C           G indicates whether we should build with the density (Coulomb)
C           or the current density matrix (Gaunt)
C             G = 0 Coulomb
C             G = 1 Gaunt
C
        DO 400 K = 1,NDMAT
          IF(IFCTYP(K) .EQ. 0) GOTO 400
C         Unpack the bitcode
          IX = IFCTYP(K)/10
          IY = MOD(IFCTYP(K),10)
          IG = IY/4
C         We can build Coulomb and Gaunt type matrices by specifying the
C         appropriate multiplication factors. 
          IF (IG.EQ.0) THEN
C............Coulomb interaction
             IC = MOD(IY,2)
             IE = IY - IC
             IF (IX.EQ.2) IC = 0
C            ... no direct term for antisymmetric density matrix
C
             IF(IREPDM(K).GE.0 .AND. IREPAB.NE.IREPDM(K)) IC = 0
C            ... IREPDM(K) is either the symmetry of the density matrix
C             or - in Dirac with NZ .lt. 4 and packed matrices (SOFOCK) -
C             the symmetry of the only symmetry block of the density matrix
C             which gives a Coulomb contribution (for IX .eq. 3 this will
C             be the symmetric component of the packed DMAT, for IX .eq. 0
C             in quadratic response all blocks are unsymmetric but only
C             the block of symmetry IREPDM(K) shall give a Coulomb contribution).
C             / hjaaj July 2004
C
          ELSE
C............Gaunt interaction
             IY = IY-4
             IC = MOD(IY,2)
             IE = IY - IC
             IF (IX.EQ.1) IC = 0
C            ... no direct term for symmetric density matrix
          ENDIF
C
C         Exchange may be switched off in DFT calculations
C
          IF (HFXFAC.EQ.D0) IE = 0
C
C         Screening on direct contributions
          IF(DOC.AND.IC.NE.0) THEN
            FCM = DNSBUF(1,K)*SOMAX
            IF(FCM.LT.SCRTHR) THEN
              DINTSKP(2,3) = DINTSKP(2,3) + NBUF
              IC = 0
            ENDIF
          ENDIF
C         Screening on exchange contributions
          IF(DOE.AND.IE.NE.0) THEN
            FEM = DNSBUF(2,K)*SOMAX
            IF(FEM.LT.SCRTHR) THEN
              DINTSKP(2,4) = DINTSKP(2,4) + NBUF
              IE = 0
            END IF
          ENDIF
          IY = IE + IC
          IF (IY.EQ.0) GOTO 400
          IF (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
            DIRFAC = D1
            EXCFAC = DMP25*HFXFAC
C
C           Symmetric singlet Fock matrix
C           =============================
C           F(i,j) = (1/4) * (FMAT(i,j) + FMAT(j,i))
C
            IF(IY.EQ.3) THEN
              DO INT = 1,NBUF
                A = IBUF(1,INT)
                B = IBUF(2,INT)
                C = IBUF(3,INT)
                D = IBUF(4,INT)
                DINT = DIRFAC*SO(INT,I)
                EINT = EXCFAC*SO(INT,I)
C...............Direct contributions
                FMAT(A,B,K) = FMAT(A,B,K) + DINT*DMAT(C,D,K)
                FMAT(C,D,K) = FMAT(C,D,K) + DINT*DMAT(A,B,K)
C...............Exchange (LL or SS for Gaunt)
                FMAT(C,A,K) = FMAT(C,A,K) + EINT*DMAT(D,B,K)
                FMAT(D,B,K) = FMAT(D,B,K) + EINT*DMAT(C,A,K)
C...............Exchange (LS or SL for Gaunt)
                FMAT(D,A,K) = FMAT(D,A,K) + EINT*DMAT(C,B,K)
                FMAT(C,B,K) = FMAT(C,B,K) + EINT*DMAT(D,A,K)
              ENDDO
C
C           Antisymmetric singlet Fock matrix OR
C           symmetric triplet Fock matrix OR
C           antisymmetric triplet Fock matrix
C           =========================================
C           F(i,j) = (1/4) (FMAT(i,j) +/- FMAT(j,i))
C
            ELSEIF(IY.EQ.2) THEN
              DO INT = 1,NBUF
                A = IBUF(1,INT)
                B = IBUF(2,INT)
                C = IBUF(3,INT)
                D = IBUF(4,INT)
                EINT = EXCFAC*SO(INT,I)
C...............Exchange (LL or SS for Gaunt)
                FMAT(C,A,K) = FMAT(C,A,K) + EINT*DMAT(D,B,K)
                FMAT(D,B,K) = FMAT(D,B,K) + EINT*DMAT(C,A,K)
C...............Exchange (LS or SL for Gaunt)
                FMAT(D,A,K) = FMAT(D,A,K) + EINT*DMAT(C,B,K)
                FMAT(C,B,K) = FMAT(C,B,K) + EINT*DMAT(D,A,K)
              ENDDO
C
C           Coulomb contributions only
C           ==========================
C
            ELSEIF(IY.EQ.1) THEN
              DO INT = 1,NBUF
                A = IBUF(1,INT)
                B = IBUF(2,INT)
                C = IBUF(3,INT)
                D = IBUF(4,INT)
C...............Direct contributions
                DINT = DIRFAC*SO(INT,I)
                FMAT(A,B,K) = FMAT(A,B,K) + DINT*DMAT(C,D,K)
                FMAT(C,D,K) = FMAT(C,D,K) + DINT*DMAT(A,B,K)
              ENDDO
C
            ELSE
              WRITE (LUPRI,'(/A,2(/A,I10))')
     &          'FCKOU1 ERROR, specified IFCTYP not implemented yet',
     &          '              specified IFCTYP was',IFCTYP(K),
     &          '              for F,D matrix no.  ',K
              CALL QUIT(
     &          'ERROR in FCKOU1: specified IFCTYP not implemented.')
            ENDIF
          ELSEIF(IX.EQ.0) THEN
            DIRFAC = D1
            EXCFAC = DMP5*HFXFAC
C
C           General singlet case - no permutational symmetry
C           ================================================
C           F(i,j) = (1/8) * FMAT(i,j)
C
            IF(IY.EQ.3) THEN
              DO INT = 1,NBUF
                A = IBUF(1,INT)
                B = IBUF(2,INT)
                C = IBUF(3,INT)
                D = IBUF(4,INT)
                DINT = DIRFAC*SO(INT,I)
C...............Direct contributions
                GCD  = DINT*(DMAT(C,D,K) + DMAT(D,C,K))
                FMAT(A,B,K) = FMAT(A,B,K) + GCD
                FMAT(B,A,K) = FMAT(B,A,K) + GCD
                GAB  = DINT*(DMAT(A,B,K) + DMAT(B,A,K))
                FMAT(C,D,K) = FMAT(C,D,K) + GAB
                FMAT(D,C,K) = FMAT(D,C,K) + GAB
                EINT = EXCFAC*SO(INT,I)
C...............Exchange (LL or SS for Gaunt)
                FMAT(C,A,K) = FMAT(C,A,K) + EINT*DMAT(D,B,K)
                FMAT(A,C,K) = FMAT(A,C,K) + EINT*DMAT(B,D,K)
                FMAT(D,B,K) = FMAT(D,B,K) + EINT*DMAT(C,A,K)
                FMAT(B,D,K) = FMAT(B,D,K) + EINT*DMAT(A,C,K)
C...............Exchange (LS or SL for Gaunt)
                FMAT(D,A,K) = FMAT(D,A,K) + EINT*DMAT(C,B,K)
                FMAT(A,D,K) = FMAT(A,D,K) + EINT*DMAT(B,C,K)
                FMAT(C,B,K) = FMAT(C,B,K) + EINT*DMAT(D,A,K)
                FMAT(B,C,K) = FMAT(B,C,K) + EINT*DMAT(A,D,K)
              ENDDO
C
C           General triplet case - no permutational symmetry
C           ================================================
C           F(i,j) = (1/8) * FMAT(i,j)
C
            ELSEIF(IY.EQ.2) THEN
              DO INT = 1,NBUF
                A = IBUF(1,INT)
                B = IBUF(2,INT)
                C = IBUF(3,INT)
                D = IBUF(4,INT)
                EINT = EXCFAC*SO(INT,I)
C...............Exchange (LL or SS for Gaunt)
                FMAT(C,A,K) = FMAT(C,A,K) + EINT*DMAT(D,B,K)
                FMAT(A,C,K) = FMAT(A,C,K) + EINT*DMAT(B,D,K)
                FMAT(D,B,K) = FMAT(D,B,K) + EINT*DMAT(C,A,K)
                FMAT(B,D,K) = FMAT(B,D,K) + EINT*DMAT(A,C,K)
C...............Exchange (LS or SL for Gaunt)
                FMAT(D,A,K) = FMAT(D,A,K) + EINT*DMAT(C,B,K)
                FMAT(A,D,K) = FMAT(A,D,K) + EINT*DMAT(B,C,K)
                FMAT(C,B,K) = FMAT(C,B,K) + EINT*DMAT(D,A,K)
                FMAT(B,C,K) = FMAT(B,C,K) + EINT*DMAT(A,D,K)
              ENDDO
C
C           General Coulomb case - no permutational symmetry
C           ================================================
C           F(i,j) = (1/8) * FMAT(i,j)
C
            ELSEIF(IY.EQ.1) THEN
              DO INT = 1,NBUF
                A = IBUF(1,INT)
                B = IBUF(2,INT)
                C = IBUF(3,INT)
                D = IBUF(4,INT)
                DINT = DIRFAC*SO(INT,I)
C...............Direct contributions
                GCD  = DINT*(DMAT(C,D,K) + DMAT(D,C,K))
                FMAT(A,B,K) = FMAT(A,B,K) + GCD
                FMAT(B,A,K) = FMAT(B,A,K) + GCD
                GAB  = DINT*(DMAT(A,B,K) + DMAT(B,A,K))
                FMAT(C,D,K) = FMAT(C,D,K) + GAB
                FMAT(D,C,K) = FMAT(D,C,K) + GAB
              ENDDO
            ELSE
              WRITE (LUPRI,'(/A,2(/A,I10))')
     &          'FCKOU1 ERROR, specified IFCTYP not implemented yet',
     &          '              specified IFCTYP was',IFCTYP(K),
     &          '              for F,D matrix no.  ',K
              CALL QUIT(
     &          'ERROR in FCKOU1: specified IFCTYP not implemented.')
            ENDIF
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'FCKOU1 ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(K),
     &        '              for F,D matrix no.  ',K
            CALL QUIT(
     &        'ERROR in FCKOU1: specified IFCTYP not implemented.')
          ENDIF
 400    CONTINUE
      ENDDO
C
      RETURN
      END
C  /* Deck dr2out */
      SUBROUTINE DR2OUT(SO,IATOM,MULE,FIRST,LAST,SPNORB,SQ12EL,THRESH,
     &                  IPNTNO,IPNTRP,IPNTLG,NINDAB,NINDCD,IPRINT)
C
C     tuh Sep 92
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
      LOGICAL NOTEST, DCMPAB, DCMPCD, DRALTB, DRCLTD, FIRST, LAST,
     &        DRABAB, DCABAB, SPNORB, SQ12EL, IPNTLG(3,NINTMX,*)
      DIMENSION IPNTNO(4,NINTMX,*), IPNTRP(3,NINTMX,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2)
      DIMENSION SO(*), BUF(600), IBUF(600)
#include "twocom.h"
#include "dorps.h"
#include "symmet.h"
#include "doxyz.h"
      SAVE BUF, IBUF, ICOUNT, NBUFCL
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine DR2OUT',-1)
C
      IF (FIRST) CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,-1,
     &                       SPNORB,IPRINT)
C
      ISOFF = 0
      DO 100 ICOOR = 1, 3
      IF (DOXYZ(ICOOR)) THEN
         DO 200 IREPE = 0,MAXREP
            IF (SPNORB) THEN
             IRPXYZ=IEOR(ISYMAX(1,1),IEOR(ISYMAX(2,1),ISYMAX(3,1)))
             IF (IREPE .NE. IEOR(ISYMAX(ICOOR,1),IRPXYZ)) GOTO 200
            ELSE
             IF(.NOT.DOREPS(IREPE)) GO TO 200
             IF(IAND(MULE,IEOR(IREPE,ISYMAX(ICOOR,1))).NE.0)GOTO 200
            END IF
            CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,2,SPNORB,
     &                  IPRINT)
            NBUFCL = 0
            IREPX = IPTREP(IREPE,2)
            NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
            DO 300 I = 1, NINTSR(IREPX)
               NSTRNA = IPNTNO(1,I,IREPX)
               NSTRNB = IPNTNO(2,I,IREPX)
               NSTRNC = IPNTNO(3,I,IREPX)
               NSTRND = IPNTNO(4,I,IREPX)
               IREPA  = IPNTRP(1,I,IREPX)
               IREPB  = IPNTRP(2,I,IREPX)
               IREPC  = IPNTRP(3,I,IREPX)
               IREPD  = IEOR(IEOR(IEOR(IREPA,IREPB),IREPC),IREPE)
               IF (NOTEST) THEN
                  INT = 0
                  DO 400 IAB = 1, NORBAB
                     IA = KHKTA*(NINDAB(IAB,1) - 1)
                     IB = KHKTB*(NINDAB(IAB,2) - 1)
                     INDA = IPTSYM(NSTRNA + IA,IREPA)
                     INDB = IPTSYM(NSTRNB + IB,IREPB)
                     INDAB = MAX(INDA,INDB)*255 + INDA + INDB
                     DO 410 ICD = 1, NORBCD
                        INT = INT + 1
                        SOINT = SO(ISOFF+INT)
                        IF (ABS(SOINT) .GT. THRESH) THEN
                           IC = KHKTC*(NINDCD(ICD,1) - 1)
                           ID = KHKTD*(NINDCD(ICD,2) - 1)
                           INDC = IPTSYM(NSTRNC + IC,IREPC)
                           INDD = IPTSYM(NSTRND + ID,IREPD)
                           IF (SQ12EL) THEN
                              IF (INDC .GT. INDD) THEN
                                 INDCD  = INDC*256 + INDD
                              ELSE IF (SPNORB) THEN
                                 SOINT = - SOINT
                                 INDCD  = INDD*256 + INDC
                              ELSE
                                 INDCD  = INDD*256 + INDC
                              END IF
                              IABCD = INDAB*65536 + INDCD
                           ELSE
                              INDCD  = MAX(INDC,INDD)*255 + INDC + INDD
                              IABCD = MAX(INDAB,INDCD)*65535+INDAB+INDCD
                           END IF
                           ICOUNT = ICOUNT + 1
                           BUF (ICOUNT) = SOINT
                           IBUF(ICOUNT) = IABCD
                           IF (ICOUNT.EQ.600) THEN
                              NBUFCL = NBUFCL + 1
                              CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                    IREPE,0,SPNORB,IPRINT)
                           END IF
                        END IF
  410                CONTINUE
  400             CONTINUE
               ELSE
                  DCMPAB = IPNTLG(1,I,IREPX)
                  DCMPCD = IPNTLG(2,I,IREPX)
                  DCABAB = IPNTLG(3,I,IREPX)
                  DRALTB = IREPA .LT. IREPB
                  DRCLTD = IREPC .LT. IREPD
                  DRABAB = DCABAB.AND.IREPA.EQ.IREPC.AND.IREPB.EQ.IREPD
                  INT = 0
                  DO 500 IAB = 1, NORBAB
                     IA = KHKTA*(NINDAB(IAB,1) - 1)
                     IB = KHKTB*(NINDAB(IAB,2) - 1)
                     IF (DCMPAB) THEN
                        IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                           INT = INT + NORBCD
                           GO TO 500
                        END IF
                     END IF
                     INDA = IPTSYM(NSTRNA + IA,IREPA)
                     INDB = IPTSYM(NSTRNB + IB,IREPB)
                     INDAB = MAX(INDA,INDB)*255 + INDA + INDB
                     DO 510 ICD = 1,NORBCD
                        IC = KHKTC*(NINDCD(ICD,1) - 1)
                        ID = KHKTD*(NINDCD(ICD,2) - 1)
                        INT = INT + 1
                        IF (DCMPCD ) THEN
                           IF (ID.GT.IC) GO TO 510
                           IF (DRCLTD .AND. ID.EQ.IC) GO TO 510
                        END IF
                        IF (DRABAB) THEN
                           IF (IA.LT.IC.OR.(IA.EQ.IC.AND.IB.LT.ID))
     &                        GOTO 510
                        END IF
                        SOINT = SO(ISOFF+INT)
                        IF (ABS(SOINT) .GT. THRESH) THEN
                           INDC = IPTSYM(NSTRNC + IC,IREPC)
                           INDD = IPTSYM(NSTRND + ID,IREPD)
                           IF (SQ12EL) THEN
                              IF (INDC .GT. INDD) THEN
                                 INDCD  = INDC*256 + INDD
                              ELSE IF (SPNORB) THEN
                                 SOINT = - SOINT
                                 INDCD  = INDD*256 + INDC
                              ELSE
                                 INDCD  = INDD*256 + INDC
                              END IF
                              IABCD = INDAB*65536 + INDCD
                           ELSE
                              INDCD  = MAX(INDC,INDD)*255 + INDC + INDD
                              IABCD = MAX(INDAB,INDCD)*65535+INDAB+INDCD
                           END IF
                           ICOUNT = ICOUNT + 1
                           BUF (ICOUNT) = SOINT
                           IBUF(ICOUNT) = IABCD
                           IF (ICOUNT.EQ.600) THEN
                              NBUFCL = NBUFCL + 1
                              CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                    IREPE,0,SPNORB,IPRINT)
                           END IF
                        END IF
  510                CONTINUE
  500             CONTINUE
               END IF
               ISOFF = ISOFF + NOABCD
  300       CONTINUE
  200    CONTINUE
      END IF
  100 CONTINUE
      IF (LAST) CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,1,SPNORB,
     &                      IPRINT)
      RETURN
      END
C  /* Deck dr2wrt */
      SUBROUTINE DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,ITYPE,SPNORB,
     &                  IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
      PARAMETER (TEN14 = 1.1D14, LU2DER = 50)
#include "hrunit.h"
      LOGICAL SPNORB
      DIMENSION BUF(600), IBUF(600)
#include "nodint.h"
#include "symmet.h"
#include "chrxyz.h"
      SAVE NBUF, NFLAGS, IRPOLD, ICROLD, ICTOLD, NBFOLD
C
C     Initialization
C     ==============
C
      IF (ITYPE .EQ. -1) THEN
         IF (SPNORB) THEN
            OPEN (LU2DER,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *            FILE='AO2SOINT')
            REWIND LU2DER
            CALL NEWLAB('AO2SOINT',LU2DER,LUPRI)
         ELSE IF (IATOM .EQ. 0) THEN
            OPEN (LU2DER,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *            FILE='AO2MGINT')
            REWIND LU2DER
            CALL NEWLAB('AO2MGINT',LU2DER,LUPRI)
         ELSE
            REWIND LU2DER
         END IF
         ICOUNT =  0
         ICTOLD =  0
         NBUF   =  0
         NBFOLD =  0
         NFLAGS =  0
         IRPOLD = -1
         ICROLD = -1
         DO 100 IREP = 0, MAXREP
            NDSINT(1,IREP) = 0
            NDSINT(2,IREP) = 0
            NDSINT(3,IREP) = 0
  100    CONTINUE
C
C     Write new buffer
C     ================
C
      ELSE IF (ITYPE .EQ. 0) THEN
         WRITE (LU2DER) BUF, IBUF, ICOUNT
         IF (IPRINT .GE. 6) CALL DR2PRI(BUF,IBUF,ICOUNT,NBUF)
         NBUF   = NBUF + 1
         ICOUNT = 0
C
C     Insert tag for new coordinate or symmetry
C     =========================================
C
      ELSE IF (ITYPE .EQ. 2) THEN
         IF (IREPE.NE.IRPOLD .OR. ICOOR.NE.ICROLD) THEN
            INTS = (NBUF - NBFOLD)*600 + ICOUNT - ICTOLD
            IF (INTS .GT. 0) NDSINT(ICROLD,IRPOLD) =
     &                       NDSINT(ICROLD,IRPOLD) + INTS
            ICOUNT = ICOUNT + 1
            BUF(ICOUNT)  = TEN14
            IBUF(ICOUNT) = IATOM*2**24+IREPE*2**16+ICOOR*2**8
            ICTOLD = ICOUNT
            NBFOLD = NBUF
            IRPOLD = IREPE
            ICROLD = ICOOR
            NFLAGS = NFLAGS + 1
            IF (ICOUNT .EQ. 600) THEN
               WRITE (LU2DER) BUF,IBUF,ICOUNT
               IF (IPRINT.GE.6) CALL DR2PRI(BUF,IBUF,ICOUNT,NBUF)
               NBUF   = NBUF + 1
               ICOUNT = 0
            END IF
         END IF
C
C     Empty final buffer and do statistics
C     ====================================
C
      ELSE
         INTS = (NBUF - NBFOLD)*600 + ICOUNT - ICTOLD
         IF (INTS .GT. 0) NDSINT(ICROLD,IRPOLD) =
     &                    NDSINT(ICROLD,IRPOLD) + INTS
         IF (ICOUNT .GT. 0) THEN
            WRITE (LU2DER) BUF, IBUF, ICOUNT
            IF (IPRINT .GE. 6) CALL DR2PRI(BUF,IBUF,ICOUNT,NBUF)
         END IF
         WRITE (LU2DER) BUF, IBUF, -1
         END FILE LU2DER
         NTOTAL = 600*NBUF + ICOUNT - NFLAGS
         NOINTS = NTOTAL
       IF (IPRINT .GE. 2) THEN
         IF (SPNORB) THEN
            CALL AROUND('Number of written 2-el. spin-orbit integrals')
         ELSE
            CALL AROUND('Number of written derivative integrals')
         END IF
         CALL HEADER('Symmetry           x         y         z',1)
         DO 500 IREP = 0, MAXREP
            IF (NDSINT(1,IREP)+NDSINT(2,IREP)+NDSINT(3,IREP).GT.0) THEN
               WRITE (LUPRI,'(2X,I5,5X,3I10)')
     &               IREP + 1, (NDSINT(I,IREP),I=1,3)
            END IF
 500     CONTINUE
         IF (SPNORB) THEN
            WRITE (LUPRI,'(/1X,I10,A,/,1X,I10,A//)')
     &         NTOTAL,' spin-orbit two-electron integrals and ',
     &         NFLAGS,' flags have been written on disk.'
         ELSE
            WRITE (LUPRI,'(/1X,I10,A,/,1X,I10,A//)')
     &         NTOTAL,' differentiated two-electron integrals and ',
     &         NFLAGS,' flags have been written on disk.'
         END IF
       END IF
       RETURN
      END IF
C
      IF (SPNORB) CLOSE(LU2DER,STATUS='KEEP')
      RETURN
      END
C  /* Deck dr2pri */
      SUBROUTINE DR2PRI(BUF,IBUF,ICOUNT,NBUF)
#include "implicit.h"
#include "priunit.h"
#include "ibtpar.h"
      DIMENSION BUF(600), IBUF(600)
#include "chrxyz.h"
      SAVE IATOM, ICOOR, IREPE
      WRITE(LUPRI,'(2X,A,I5,A/)')
     &   'Integral buffer #',NBUF + 1,' has been written.'
      DO 100 INT = 1, ICOUNT
         I = IAND(ISHFT(IBUF(INT),-24),IBT08)
         J = IAND(ISHFT(IBUF(INT),-16),IBT08)
         K = IAND(ISHFT(IBUF(INT),-8), IBT08)
         L = IAND(      IBUF(INT),     IBT08)
         IF (L .EQ. 0) THEN
            ICOOR = IAND(ISHFT(IBUF(INT), -8),IBT08)
            IREPE = IAND(ISHFT(IBUF(INT),-16),IBT08)
            IATOM = IAND(ISHFT(IBUF(INT),-24),IBT08)
         ELSE
            WRITE (LUPRI,'(10X,A,2X,4I4,5X,I2,A,I2,5X,1P,D16.8)')
     &         ' ## ',I,J,K,L,IATOM,CHRXYZ(ICOOR),IREPE,BUF(INT)
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck mg2out */
      SUBROUTINE MG2OUT(SO,FIRST,LAST,THRESH,IPNTNO,IPNTRP,IPNTLG,
     &                  NINDAB,NINDCD,IPRINT)
C
C     tuh Sep 92
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
      LOGICAL NOTEST, DCMPAB, DCMPCD, DRALTB, DRCLTD, FIRST, LAST,
     &        DRABAB, DCABAB, SPNORB, IPNTLG(3,NINTMX,*)
      DIMENSION IPNTNO(4,NINTMX,*), IPNTRP(3,NINTMX,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2)
      DIMENSION SO(*), BUF(600), IBUF(600)
#include "twocom.h"
#include "dorps.h"
#include "symmet.h"
#include "doxyz.h"
      SAVE BUF, IBUF, ICOUNT, NBUFCL
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine MG2OUT',-1)
C
      IATOM = 0
      SPNORB = .FALSE.
      IF (FIRST) CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,-1,
     &                       SPNORB,IPRINT)
C
      NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
      NBUFCL = 0
      ISOFF0  = 0
      DO 100 ICOOR = 1, 3
         IREPE = ISYMAX(ICOOR,2)
         CALL DR2WRT(BUF,IBUF,ICOUNT,0,ICOOR,IREPE,2,.FALSE.,IPRINT)
         IREPX = IPTREP(IREPE,2)
         NINTX = NINTSR(IREPX)
         ISOFF1 = ISOFF0
         ISOFF2 = ISOFF0 + NOABCD*NINTX
         DO 200 I = 1, NINTX
            NSTRNA = IPNTNO(1,I,IREPX)
            NSTRNB = IPNTNO(2,I,IREPX)
            NSTRNC = IPNTNO(3,I,IREPX)
            NSTRND = IPNTNO(4,I,IREPX)
            IREPA  = IPNTRP(1,I,IREPX)
            IREPB  = IPNTRP(2,I,IREPX)
            IREPC  = IPNTRP(3,I,IREPX)
            IREPD  = IEOR(IEOR(IEOR(IREPA,IREPB),IREPC),IREPE)
            IF (NOTEST) THEN
               INT = 0
               DO 400 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  DO 410 ICD = 1, NORBCD
                     IC = KHKTC*(NINDCD(ICD,1) - 1)
                     ID = KHKTD*(NINDCD(ICD,2) - 1)
                     INDC = IPTSYM(NSTRNC + IC,IREPC)
                     INDD = IPTSYM(NSTRND + ID,IREPD)
                     INT = INT + 1
                     SOINT = SO(ISOFF1 + INT)
C                    print *, '##n+',soint,inda,indb,indc,indd,icoor
                     IF (ABS(SOINT) .GT. THRESH) THEN
                        INDAB = INDA*256 + INDB
                        INDCD = INDC*256 + INDD
                        ICOUNT = ICOUNT + 1
                        BUF (ICOUNT) = SOINT
                        IBUF(ICOUNT) = INDAB*65536 + INDCD
                        IF (ICOUNT.EQ.600) THEN
                           NBUFCL = NBUFCL + 1
                           CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                 IREPE,0,SPNORB,IPRINT)
                        END IF
                     END IF
                     SOINT = SO(ISOFF2 + INT)
C                    print *, '##n-',soint,inda,indb,indd,indc,icoor
                     IF (ABS(SOINT) .GT. THRESH) THEN
                        INDAB = INDA*256 + INDB
                        INDCD = INDD*256 + INDC
                        ICOUNT = ICOUNT + 1
                        BUF (ICOUNT) = SOINT
                        IBUF(ICOUNT) = INDAB*65536 + INDCD
                        IF (ICOUNT.EQ.600) THEN
                           NBUFCL = NBUFCL + 1
                           CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                 IREPE,0,SPNORB,IPRINT)
                        END IF
                     END IF
  410             CONTINUE
  400          CONTINUE
            ELSE
               DCMPAB = IPNTLG(1,I,IREPX)
               DCMPCD = IPNTLG(2,I,IREPX)
               DCABAB = IPNTLG(3,I,IREPX)
               DRALTB = IREPA .LT. IREPB
               DRCLTD = IREPC .LT. IREPD
               DRABAB = DCABAB.AND.IREPA.EQ.IREPC.AND.IREPB.EQ.IREPD
               INT = 0
               DO 500 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  IF (DCMPAB) THEN
                     IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                        INT = INT + NORBCD
                        GO TO 500
                     END IF
                  END IF
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = MAX(INDA,INDB)*256 + INDA + INDB
                  DO 510 ICD = 1,NORBCD
                     IC   = KHKTC*(NINDCD(ICD,1) - 1)
                     ID   = KHKTD*(NINDCD(ICD,2) - 1)
                     INDC = IPTSYM(NSTRNC + IC,IREPC)
                     INDD = IPTSYM(NSTRND + ID,IREPD)
                     INT = INT + 1
                     IF (DCMPCD ) THEN
                        IF (ID.GT.IC) GO TO 510
                        IF (DRCLTD .AND. ID.EQ.IC) GO TO 510
                     END IF
                     IF (DRABAB) THEN
                        IF (IA.LT.IC.OR.(IA.EQ.IC.AND.IB.LT.ID))
     &                     GOTO 510
                     END IF
                     SOINT = SO(ISOFF1 + INT)
C                    print *, '##t+',soint,inda,indb,indc,indd,icoor
                     IF (ABS(SOINT) .GT. THRESH) THEN
                        INDAB = INDA*256 + INDB
                        INDCD = INDC*256 + INDD
                        ICOUNT = ICOUNT + 1
                        BUF (ICOUNT) = SOINT
                        IBUF(ICOUNT) = INDAB*65536 + INDCD
                        IF (ICOUNT.EQ.600) THEN
                           NBUFCL = NBUFCL + 1
                           CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                 IREPE,0,SPNORB,IPRINT)
                        END IF
                     END IF
                     SOINT = SO(ISOFF2 + INT)
C                    print *, '##t-',soint,inda,indb,indd,indc,icoor
                     IF (ABS(SOINT) .GT. THRESH) THEN
                     IF (INDA .EQ. INDB) GO TO 510
                     IF (INDC .EQ. INDD) GO TO 510
                        INDAB = INDA*256 + INDB
                        INDDC = INDD*256 + INDC
                        ICOUNT = ICOUNT + 1
                        BUF (ICOUNT) = SOINT
                        IBUF(ICOUNT) = INDAB*65536 + INDDC
                        IF (ICOUNT.EQ.600) THEN
                           NBUFCL = NBUFCL + 1
                           CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                 IREPE,0,SPNORB,IPRINT)
                        END IF
                     END IF
  510             CONTINUE
  500          CONTINUE
            END IF
            ISOFF1 = ISOFF1 + NOABCD
            ISOFF2 = ISOFF2 + NOABCD
  200    CONTINUE
         ISOFF0 = ISOFF0 + 2*NINTX*NOABCD
  100 CONTINUE
      IF (LAST) CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,1,SPNORB,
     &                      IPRINT)
      RETURN
      END
C  /* Deck ds2out */
      SUBROUTINE DS2OUT(SO,WRKBUF,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,
     &                  THRESH,NINDAB,NINDCD,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      LOGICAL FIRST, LAST, IPNTLG(*)
      DIMENSION SO(*), WRKBUF(*), IPNTNO(*), IPNTRP(*), NINDAB(*),
     &          NINDCD(*)
#include "disbuf.h"
C
C---------------------------------
C     Call sort and write routine.
C---------------------------------
C
      CALL DS2OU1(SO,WRKBUF(KDSBF),WRKBUF(KDSIBF),WRKBUF(KDSNCT),
     &            WRKBUF(KDSORB),WRKBUF(KORBDS),IPNTNO,IPNTRP,IPNTLG,
     &            FIRST,LAST,THRESH,NINDAB,NINDCD,LDSBUF,NDIST,IPRINT)
C
      RETURN
      END
C  /* Deck ds2ou1 */
      SUBROUTINE DS2OU1(SO,BUF,IBUF,NCOUNT,IDSORB,IORBDS,IPNTNO,IPNTRP,
     &                  IPNTLG,FIRST,LAST,THRESH,NINDAB,NINDCD,LBUF,
     &                  NDIST,IPRINT)
C
C     Write out blocks of symmetry integrals, eliminating duplicates
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
      LOGICAL DCMPCD, DRCLTD, FIRST, LAST, IPNTLG(3,*), NOTEST,
     &        DOINDX
      DIMENSION SO(*), BUF(LBUF,NDIST), IBUF(LBUF,2,NDIST),
     &          IPNTNO(4,*), IPNTRP(3,*), NCOUNT(NDIST),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          IDSORB(NDIST), IORBDS(NBASIS)
#include "nuclei.h"
#include "twocom.h"
#include "symmet.h"
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine DS2OUT',-1)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(2X,A,4I5)') 'NHKT? ', NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI,'(2X,A,4I5)') 'MUL?  ', MULA,  MULB,  MULC,  MULD
         WRITE (LUPRI,'(2X,A,4I5)') 'NORB? ', NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI,'(2X,A,4I5)') 'NSTR? ', NSTRA, NSTRB, NSTRC, NSTRD
         WRITE (LUPRI,'(2X,A,2I5)') 'NORBCD', NORBCD
         WRITE (LUPRI,'(2X,A,2I5)') 'NOABCD', NOABCD
         WRITE (LUPRI,'(2X,A,2L5)') 'DIAGAB/CD', DIAGAB, DIAGCD
         WRITE (LUPRI,'(2X,A,2L5)') 'TCONAB/CD', TCONAB, TCONCD
         WRITE (LUPRI,'(2X,A,2L5)') 'SHAEQB/CD', SHAEQB, SHCEQD
         WRITE (LUPRI,'(2X,A, L5)') 'SHABAB', SHABAB
         WRITE (LUPRI,'(2X,A, I5)') 'NDIST ', NDIST
      END IF
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (INT_STAR8)
      NIBUF = 1
      IF (NBASIS .LE. 255) THEN
         NBITS = 8
      ELSE
         NBITS = 16
      END IF
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
C     *******************************************************
C     ***** Initialization when subroutine first called *****
C     *******************************************************
C
      IF (FIRST) THEN
         CALL IZERO(NCOUNT,NDIST)
         CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,-1,NBITS,0,IPRINT)
         DOINDX = .TRUE.
         CALL AINDEX(ISHELA,NAINTS,IDSORB,DOINDX,IPRINT)
         DO 50 IDIST = 1, NDIST
             IORBDS(IDSORB(IDIST)) = IDIST
   50    CONTINUE
         NWRIT = 0
      END IF
C
      ISOFF  = 0
      NBUFCL = 0
      NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
      DO 100 I = 1, NINTS
         NSTRNA = IPNTNO(1,I)
         NSTRNB = IPNTNO(2,I)
         NSTRNC = IPNTNO(3,I)
         NSTRND = IPNTNO(4,I)
         IREPA  = IPNTRP(1,I)
         IREPB  = IPNTRP(2,I)
         IREPC  = IPNTRP(3,I)
         IREPD  = IEOR(IEOR(IREPA,IREPB),IREPC)
         IF (NOTEST) THEN
            IF (NIBUF .EQ. 1) THEN
               INT = 0
               DO 200 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDAB = INDA*(IBIT1 + 1) + IPTSYM(NSTRNB + IB,IREPB)
                  IDIST = IORBDS(INDA)
                  DO 210 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT = SO(ISOFF+INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       BUF (ICOUNT,IDIST) = SOINT
                       IBUF(ICOUNT,1,IDIST) = INDAB*(IBIT2 + 1) + INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
                          CALL UN2WRT(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,
     &                                NIBUF,ICOUNT,0,NBITS,INDA,IPRINT)
                          NCOUNT(IDIST) = 0
                       END IF
                     END IF
  210             CONTINUE
  200          CONTINUE
            ELSE
               INT = 0
               DO 205 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = INDA*(IBIT1 + 1) + INDB
                  IDIST = IORBDS(INDA)
                  DO 215 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT = SO(ISOFF+INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       BUF (ICOUNT,IDIST) = SOINT
                       IBUF(ICOUNT,1,IDIST) = INDAB
                       IBUF(ICOUNT,2,IDIST) = INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
                          CALL UN2WRT(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,
     &                                NIBUF,ICOUNT,0,NBITS,INDA,IPRINT)
                          NCOUNT(IDIST) = 0
                       END IF
                     END IF
  215             CONTINUE
  205          CONTINUE
            END IF
         ELSE
            DCMPCD = IPNTLG(2,I)
            DRCLTD = IREPC .LT. IREPD
            INT = 0
            DO 300 IAB = 1, NORBAB
               IA = KHKTA*(NINDAB(IAB,1) - 1)
               IB = KHKTB*(NINDAB(IAB,2) - 1)
               INDA = IPTSYM(NSTRNA + IA,IREPA)
               INDB = IPTSYM(NSTRNB + IB,IREPB)
               INDAB = INDA*(IBIT1 + 1) + INDB
               IDIST = IORBDS(INDA)
               DO 310 ICD = 1,NORBCD
                  IC = KHKTC*(NINDCD(ICD,1) - 1)
                  ID = KHKTD*(NINDCD(ICD,2) - 1)
                  INT = INT + 1
                  IF (DCMPCD ) THEN
                     IF (ID.GT.IC) GO TO 310
                     IF (DRCLTD .AND. ID.EQ.IC) GO TO 310
                  END IF
                  SOINT = SO(ISOFF+INT)
                  IF (ABS(SOINT) .GT. THRESH) THEN
                     IF (NIBUF .EQ. 1) THEN
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       BUF (ICOUNT,IDIST) = SOINT
                       IBUF(ICOUNT,1,IDIST) = INDAB*(IBIT2 + 1) + INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          CALL UN2WRT(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,
     &                                NIBUF,ICOUNT,0,NBITS,INDA,IPRINT)
                          NBUFCL = NBUFCL + 1
                          NCOUNT(IDIST) = 0
                       END IF
                     ELSE
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       BUF (ICOUNT,IDIST) = SOINT
                       IBUF(ICOUNT,1,IDIST) = INDAB
                       IBUF(ICOUNT,2,IDIST) = INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          CALL UN2WRT(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,
     &                                NIBUF,ICOUNT,0,NBITS,INDA,IPRINT)
                          NBUFCL = NBUFCL + 1
                          NCOUNT(IDIST) = 0
                       END IF
                     END IF
                  END IF
  310          CONTINUE
  300       CONTINUE
         END IF
         ISOFF = ISOFF + NOABCD
  100 CONTINUE
      NWRIT = NWRIT + LBUF*NBUFCL
C
C     *************************************
C     ***** Last call to empty buffer *****
C     *************************************
C
      IF (LAST) THEN
         DO 400 IDIST = 1, NDIST
            NWRIT = NWRIT + NCOUNT(IDIST)
            CALL UN2WRT(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,NIBUF,
     &                  NCOUNT(IDIST),1,NBITS,IDSORB(IDIST),IPRINT)
  400    CONTINUE
         NALL = NDIST*(NBASIS + 1)*(NBASIS**2)/2
         PERCNT = dble(100*NWRIT) / dble(NALL)
         WRITE (LUPRI,'(/1X,A,I10,A,F4.1,A)')
     &         'Number of two-electron integrals written:',NWRIT,
     &         ' (',PERCNT,'%)'
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dt2out */
      SUBROUTINE DT2OUT(SO,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,THRESH,
     &                  NINDAB,NINDCD,IPRINT,INDX,INDXAB,GMAT)
C*****************************************************************************
C
C     Write out distribution blocks of symmetry integrals
C     Luuk Visscher, may 1997, a rainy day
C
C     Two options : we write directly to gmat if npass = 1
C                   we write to a direct access file if npass > 1
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "twosta.h"
C
#include "nuclei.h"
#include "twocom.h"
#include "blocks.h"
#include "symmet.h"
#include "comdis.h"
      LOGICAL FIRST, LAST, IPNTLG(3,*)
      LOGICAL DRALTB,DCMPAB,NOTEST
      DIMENSION SO(*),GMAT(*),
     &          IPNTNO(4,*), IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          INDX(3,*),INDXAB(NINSHA,NINSHB,5)
C
      IF (NPASS.EQ.1) THEN
         CALL DT2UT1(SO,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,THRESH,
     &               NINDAB,NINDCD,IPRINT,INDX,INDXAB,GMAT)
      ELSE
C
C        In this case GMAT functions as a structure in which the 
C        buffer arrays and information are stored. The pointers are
C        kept in the common block comdis.
C
         CALL DT2UT2(SO,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,THRESH,
     &               NINDAB,NINDCD,IPRINT,INDX,INDXAB,
     &               GMAT(KLGREC),GMAT(KLGBUF),GMAT(KRGBUF),
     &               GMAT(KIGBUF))
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dt2ut1 */
      SUBROUTINE DT2UT1(SO,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,THRESH,
     &                  NINDAB,NINDCD,IPRINT,INDX,INDXAB,GMAT)
C*****************************************************************************
C
C     Write out distribution blocks of symmetry integrals to matrix gmat
C     Luuk Visscher, october 1996, a rainy day
C
C     Modified FCKOUT - we keep the cd array square, ab triangular.
C     No need to screen !
C
C     SO     - SO-integrals (pq|rs) with pq fixed
C     GMAT   - sorted batch of SO-integrals
C     INDXAB - pointers from given p and q (relative in shell)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "comdis.h"
#include "twosta.h"
      LOGICAL FIRST, LAST, IPNTLG(3,*)
      LOGICAL DRALTB,DCMPAB,NOTEST
      DIMENSION SO(*),GMAT(*),
     &          IPNTNO(4,*), IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          INDX(3,*),INDXAB(NINSHA,NINSHB,5)
#include "nuclei.h"
#include "twocom.h"
#include "blocks.h"
#include "symmet.h"
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine DT2OUT',-1)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(2X,A,4I5)') 'NHKT? ', NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI,'(2X,A,4I5)') 'MUL?  ', MULA,  MULB,  MULC,  MULD
         WRITE (LUPRI,'(2X,A,4I5)') 'NORB? ', NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI,'(2X,A,4I5)') 'NSTR? ', NSTRA, NSTRB, NSTRC, NSTRD
         WRITE (LUPRI,'(2X,A,2I5)') 'NORBCD', NORBCD
         WRITE (LUPRI,'(2X,A,2I5)') 'NOABCD', NOABCD
      END IF
C
      NALL = 0
C
C     **********************************
C     ***** Fill integrals in GMAT *****
C     **********************************
C
      NOTEST = .NOT.SHAEQB
      ICK = LCLASH(ISHELC)
      ICL = LCLASH(ISHELD)
C     Check whether we have Gaunt or Coulomb integrals
      IF (ICK.EQ.ICL) THEN
         ICKL = ICK
      ELSE
         ICKL = 3
      ENDIF
      ISOFF  = 0
      IF (NOTEST) THEN
         DO I = 1, NINTS
           NSTRNA = IPNTNO(1,I)
           NSTRNB = IPNTNO(2,I)
           NSTRNC = IPNTNO(3,I)
           NSTRND = IPNTNO(4,I)
           IREPA  = IPNTRP(1,I)
           IREPB  = IPNTRP(2,I)
           IREPC  = IPNTRP(3,I)
           IREPD  = IEOR(IEOR(IREPA,IREPB),IREPC)
           INT = 0
           DO IAB = 1, NORBAB
             IA = KHKTA*(NINDAB(IAB,1) - 1)
             IB = KHKTB*(NINDAB(IAB,2) - 1)
             INDA = IPTSYM(NSTRNA + IA,IREPA)
             INDB = IPTSYM(NSTRNB + IB,IREPB)
             INDAR = INDX(1,INDA)
             INDBR = INDX(1,INDB)
             IOFFD = INDXAB(INDAR,INDBR,2)
             DO ICD = 1, NORBCD
               INT   = INT + 1
               SOINT = SO(ISOFF+INT)
              IF (ABS(SOINT) .GT. THRESH) THEN
               IC    = KHKTC*(NINDCD(ICD,1) - 1)
               ID    = KHKTD*(NINDCD(ICD,2) - 1)
               INDC  = IPTSYM(NSTRNC + IC,IREPC)
               INDD  = IPTSYM(NSTRND + ID,IREPD)
               KRP   = INDX(2,INDC)
               LRP   = INDX(2,INDD)
               INDCR = INDC - ICOS(KRP+1,ICK)
               INDDR = INDD - ICOS(LRP+1,ICL)
               IOFFKL = ISPCK(KRP,LRP,ICKL)
     &                + (INDDR-1)*NCOS(KRP+1,ICK)+INDCR
               IOFFLK = ISPCK(LRP,KRP,ICKL)
     &                + (INDCR-1)*NCOS(LRP+1,ICL)+INDDR
               GMAT(IOFFD+IOFFKL) = SOINT
               IF (ICKL.NE.3) GMAT(IOFFD+IOFFLK) = SOINT
               NALL = NALL + 1
               IF (IOFFKL.NE.IOFFLK.AND.ICKL.NE.3)
     &         NALL = NALL + 1
              END IF
             ENDDO
           ENDDO
           ISOFF = ISOFF + NOABCD
         ENDDO
      ELSE
         DO 100 I = 1, NINTS
           NSTRNA = IPNTNO(1,I)
           NSTRNB = IPNTNO(2,I)
           NSTRNC = IPNTNO(3,I)
           NSTRND = IPNTNO(4,I)
           IREPA  = IPNTRP(1,I)
           IREPB  = IPNTRP(2,I)
           IREPC  = IPNTRP(3,I)
           IREPD  = IEOR(IEOR(IREPA,IREPB),IREPC)
           DCMPAB = IPNTLG(1,I)
           DRALTB = IREPA .LT. IREPB
           INT = 0
           DO 200 IAB = 1, NORBAB
             IA = KHKTA*(NINDAB(IAB,1) - 1)
             IB = KHKTB*(NINDAB(IAB,2) - 1)
             IF (DCMPAB) THEN
                IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                   INT = INT + NORBCD
                   GO TO 200
                END IF
             END IF
             INDA = IPTSYM(NSTRNA + IA,IREPA)
             INDB = IPTSYM(NSTRNB + IB,IREPB)
             INDAR = INDX(1,INDA)
             INDBR = INDX(1,INDB)
             IOFFD = INDXAB(INDAR,INDBR,2)
             DO 300 ICD = 1, NORBCD
               INT   = INT + 1
               SOINT = SO(ISOFF+INT)
              IF (ABS(SOINT) .GT. THRESH) THEN
               IC    = KHKTC*(NINDCD(ICD,1) - 1)
               ID    = KHKTD*(NINDCD(ICD,2) - 1)
               INDC  = IPTSYM(NSTRNC + IC,IREPC)
               INDD  = IPTSYM(NSTRND + ID,IREPD)
               KRP   = INDX(2,INDC)
               LRP   = INDX(2,INDD)
C
               INDCR = INDC - ICOS(KRP+1,ICK)
               INDDR = INDD - ICOS(LRP+1,ICL)
               IOFFKL = ISPCK(KRP,LRP,ICKL)
     &                + (INDDR-1)*NCOS(KRP+1,ICK)+INDCR
               IOFFLK = ISPCK(LRP,KRP,ICKL)
     &                + (INDCR-1)*NCOS(LRP+1,ICL)+INDDR
               GMAT(IOFFD+IOFFKL) = SOINT
               IF (ICKL.NE.3) GMAT(IOFFD+IOFFLK) = SOINT
               NALL = NALL + 1
               IF (IOFFKL.NE.IOFFLK.AND.ICKL.NE.3)
     &         NALL = NALL + 1
              END IF
  300        CONTINUE
  200      CONTINUE
           ISOFF = ISOFF + NOABCD
  100  CONTINUE
      ENDIF
C
      N2WRIT = N2WRIT + NALL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dt2ut2 */
      SUBROUTINE DT2UT2(SO,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,THRESH,
     &                  NINDAB,NINDCD,IPRINT,INDX,INDXAB,
     &                  LGREC,LGBUF,RGBUF,IGBUF)
C*****************************************************************************
C
C     Write out distribution blocks of symmetry integrals to 
C     direct acces file
C     Luuk Visscher, may 1997, another rainy day
C
C     Modified FCKOUT - we keep the cd array square, ab triangular.
C     No need to screen !
C
C     SO     - SO-integrals (pq|rs) with pq fixed
C     GMAT   - sorted batch of SO-integrals
C     INDXAB - pointers from given p and q (relative in shell)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "comdis.h"
#include "twosta.h"
      LOGICAL FIRST, LAST, IPNTLG(3,*)
      LOGICAL DRALTB,DCMPAB,NOTEST
      DIMENSION SO(*),
     &          IPNTNO(4,*), IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          INDX(3,*),INDXAB(NINSHA,NINSHB,5)
      DIMENSION LGREC(NPASS),LGBUF(NPASS)
      DIMENSION RGBUF(NGBFSZ,NPASS),IGBUF(NGBFSZ,NPASS)
#include "nuclei.h"
#include "twocom.h"
#include "blocks.h"
#include "symmet.h"
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine DT2OUT',-1)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(2X,A,4I5)') 'NHKT? ', NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI,'(2X,A,4I5)') 'MUL?  ', MULA,  MULB,  MULC,  MULD
         WRITE (LUPRI,'(2X,A,4I5)') 'NORB? ', NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI,'(2X,A,4I5)') 'NSTR? ', NSTRA, NSTRB, NSTRC, NSTRD
         WRITE (LUPRI,'(2X,A,2I5)') 'NORBCD', NORBCD
         WRITE (LUPRI,'(2X,A,2I5)') 'NOABCD', NOABCD
      END IF
C
      NALL = 0
C
C     **********************************
C     ***** Fill integrals in GMAT *****
C     **********************************
C
      NOTEST = .NOT.SHAEQB
      ISOFF  = 0
      ICK = LCLASH(ISHELC)
      ICL = LCLASH(ISHELD)
C     Check whether we have Gaunt or Coulomb integrals
      IF (ICK.EQ.ICL) THEN
         ICKL = ICK
      ELSE
         ICKL = 3
      ENDIF
      IF (NOTEST) THEN
         DO I = 1, NINTS
           NSTRNA = IPNTNO(1,I)
           NSTRNB = IPNTNO(2,I)
           NSTRNC = IPNTNO(3,I)
           NSTRND = IPNTNO(4,I)
           IREPA  = IPNTRP(1,I)
           IREPB  = IPNTRP(2,I)
           IREPC  = IPNTRP(3,I)
           IREPD  = IEOR(IEOR(IREPA,IREPB),IREPC)
           INT = 0
           DO IAB = 1, NORBAB
             IA = KHKTA*(NINDAB(IAB,1) - 1)
             IB = KHKTB*(NINDAB(IAB,2) - 1)
             INDA = IPTSYM(NSTRNA + IA,IREPA)
             INDB = IPTSYM(NSTRNB + IB,IREPB)
             INDAR = INDX(1,INDA)
             INDBR = INDX(1,INDB)
             IOFFD = INDXAB(INDAR,INDBR,2)
             IPASS = INDXAB(INDAR,INDBR,5)
             KBUF  = LGBUF(IPASS)
             DO ICD = 1, NORBCD
               INT   = INT + 1
               SOINT = SO(ISOFF+INT)
              IF (ABS(SOINT) .GT. THRESH) THEN
               IC    = KHKTC*(NINDCD(ICD,1) - 1)
               ID    = KHKTD*(NINDCD(ICD,2) - 1)
               INDC  = IPTSYM(NSTRNC + IC,IREPC)
               INDD  = IPTSYM(NSTRND + ID,IREPD)
               KRP   = INDX(2,INDC)
               LRP   = INDX(2,INDD)
               INDCR = INDC - ICOS(KRP+1,ICK)
               INDDR = INDD - ICOS(LRP+1,ICL)
               IOFFKL = ISPCK(KRP,LRP,ICKL)
     &                + (INDDR-1)*NCOS(KRP+1,ICK)+INDCR
               IOFFLK = ISPCK(LRP,KRP,ICKL)
     &                + (INDCR-1)*NCOS(LRP+1,ICL)+INDDR
C
               KBUF = KBUF + 1
               RGBUF(KBUF,IPASS) = SOINT
               IGBUF(KBUF,IPASS) = IOFFD+IOFFKL
               IF (KBUF.GE.NGBFSZ) THEN
                  IRECG = IRECG + 1
                  WRITE (LGFIL,REC=IRECG) LGREC(IPASS),KBUF,
     &                  (RGBUF(JBUF,IPASS),JBUF=1,KBUF),
     &                  (IGBUF(JBUF,IPASS),JBUF=1,KBUF)
                  LGREC(IPASS) = IRECG
                  LGBUF(IPASS) = 0
                  KBUF = 0
               ENDIF
               NALL = NALL + 1
C
               IF (IOFFKL.NE.IOFFLK.AND.ICKL.NE.3) THEN
                  KBUF = KBUF + 1
                  RGBUF(KBUF,IPASS) = SOINT
                  IGBUF(KBUF,IPASS) = IOFFD+IOFFLK
                  IF (KBUF.GE.NGBFSZ) THEN
                     IRECG = IRECG + 1
                     WRITE (LGFIL,REC=IRECG) LGREC(IPASS),KBUF,
     &                     (RGBUF(JBUF,IPASS),JBUF=1,KBUF),
     &                     (IGBUF(JBUF,IPASS),JBUF=1,KBUF)
                     LGREC(IPASS) = IRECG
                     LGBUF(IPASS) = 0
                     KBUF = 0
                  ENDIF
                  NALL = NALL + 1
               ENDIF
              END IF
             ENDDO
             LGBUF(IPASS) = KBUF
           ENDDO
           ISOFF = ISOFF + NOABCD
         ENDDO
      ELSE
         DO 100 I = 1, NINTS
           NSTRNA = IPNTNO(1,I)
           NSTRNB = IPNTNO(2,I)
           NSTRNC = IPNTNO(3,I)
           NSTRND = IPNTNO(4,I)
           IREPA  = IPNTRP(1,I)
           IREPB  = IPNTRP(2,I)
           IREPC  = IPNTRP(3,I)
           IREPD  = IEOR(IEOR(IREPA,IREPB),IREPC)
           DCMPAB = IPNTLG(1,I)
           DRALTB = IREPA .LT. IREPB
           INT = 0
           DO 200 IAB = 1, NORBAB
             IA = KHKTA*(NINDAB(IAB,1) - 1)
             IB = KHKTB*(NINDAB(IAB,2) - 1)
             IF (DCMPAB) THEN
                IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                   INT = INT + NORBCD
                   GO TO 200
                END IF
             END IF
             INDA = IPTSYM(NSTRNA + IA,IREPA)
             INDB = IPTSYM(NSTRNB + IB,IREPB)
             INDAR = INDX(1,INDA)
             INDBR = INDX(1,INDB)
             IOFFD = INDXAB(INDAR,INDBR,2)
             IPASS = INDXAB(INDAR,INDBR,5)
             KBUF  = LGBUF(IPASS)
             DO 300 ICD = 1, NORBCD
               INT   = INT + 1
               SOINT = SO(ISOFF+INT)
              IF (ABS(SOINT) .GT. THRESH) THEN
               IC    = KHKTC*(NINDCD(ICD,1) - 1)
               ID    = KHKTD*(NINDCD(ICD,2) - 1)
               INDC  = IPTSYM(NSTRNC + IC,IREPC)
               INDD  = IPTSYM(NSTRND + ID,IREPD)
               KRP   = INDX(2,INDC)
               LRP   = INDX(2,INDD)
C
               INDCR = INDC - ICOS(KRP+1,ICK)
               INDDR = INDD - ICOS(LRP+1,ICL)
               IOFFKL = ISPCK(KRP,LRP,ICKL)
     &                + (INDDR-1)*NCOS(KRP+1,ICK)+INDCR
               IOFFLK = ISPCK(LRP,KRP,ICKL)
     &                + (INDCR-1)*NCOS(LRP+1,ICL)+INDDR
C
               KBUF = KBUF + 1
               RGBUF(KBUF,IPASS) = SOINT
               IGBUF(KBUF,IPASS) = IOFFD+IOFFKL
               IF (KBUF.GE.NGBFSZ) THEN
                  IRECG = IRECG + 1
                  WRITE (LGFIL,REC=IRECG) LGREC(IPASS),KBUF,
     &                  (RGBUF(JBUF,IPASS),JBUF=1,KBUF),
     &                  (IGBUF(JBUF,IPASS),JBUF=1,KBUF)
                  LGREC(IPASS) = IRECG
                  LGBUF(IPASS) = 0
                  KBUF = 0
               ENDIF
               NALL = NALL + 1
C
               IF (IOFFKL.NE.IOFFLK.AND.ICKL.NE.3) THEN
                  KBUF = KBUF + 1
                  RGBUF(KBUF,IPASS) = SOINT
                  IGBUF(KBUF,IPASS) = IOFFD+IOFFLK
                  IF (KBUF.GE.NGBFSZ) THEN
                     IRECG = IRECG + 1
                     WRITE (LGFIL,REC=IRECG) LGREC(IPASS),KBUF,
     &                     (RGBUF(JBUF,IPASS),JBUF=1,KBUF),
     &                     (IGBUF(JBUF,IPASS),JBUF=1,KBUF)
                     LGREC(IPASS) = IRECG
                     LGBUF(IPASS) = 0
                     KBUF = 0
                  ENDIF
                  NALL = NALL + 1
               ENDIF
              END IF
C
  300        CONTINUE
             LGBUF(IPASS) = KBUF
  200      CONTINUE
           ISOFF = ISOFF + NOABCD
  100  CONTINUE
      ENDIF
C
      N2WRIT = N2WRIT + NALL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Tstint */
      LOGICAL FUNCTION TST2IN(I,J,K,L,INTCLS)
C***********************************************************************
C
C     Check for integral type
C
C      01 - XX: Class D   XY: Class D
C      02 - XX: Class Ca
C      03 - XX: Class Cb  XY: Class B
C      04 - XX: Class Bb2
C      05 - XX: Class Ba1
C      06 - XX: Class Ba2
C      07 - XX: Class Bb1 XY: Class A
C      08 - XX: Class A
C      09 -               XY: Class B
C      10 -               XY: Class C
C
C     Written by T.Saue June 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      GO TO ( 1,2,3,4,5,6,7,8,9,10),INTCLS
      STOP 'Error in TST2IN'
 1    CONTINUE
C    *Group D : all indices different(sorted in XX1SRT):
C       Class  1: I > K > J > L
C       Class  2: I > K > L > J
C       Class  3: I > J > K > L
        TST2IN = (I.NE.J).AND.(I.NE.K).AND.(J.NE.K).AND.
     &         (J.NE.L).AND.(K.NE.L)
        GOTO 999
 2    CONTINUE
C    *Group C: two indices equal
C       Class a1: I = K > J > L
C       Class a2: I > J = K > L
C       Class a3: I > K > J = L
        TST2IN = ((I.EQ.K).AND.(K.GT.J).AND.(K.GT.L).AND.(J.GT.L))
     &       .OR.((I.GT.J).AND.(J.EQ.K).AND.(I.GT.L).AND.(K.GT.L))
     &       .OR.((I.GT.K).AND.(I.GT.J).AND.(K.GT.J).AND.(J.EQ.L))
        GOTO 999
 3    CONTINUE
C    *Group C: two indices equal
C       Class b1: I > K = L > J
C       Class b2: I > J > K = L
C       Class b3: I = J > K > L
        TST2IN = ((I.GT.K).AND.(I.GT.J).AND.(K.EQ.L).AND.(L.GT.J))
     &       .OR.((I.GT.J).AND.(I.GT.K).AND.(J.GT.K).AND.(K.EQ.L))
     &       .OR.((I.EQ.J).AND.(I.GT.K).AND.(I.GT.L).AND.(K.GT.L))
        GOTO 999
 4    CONTINUE
C    *Group B: Three indices equal
C       Class b2: I = K > J = L
        TST2IN = (I.EQ.K).AND.(J.EQ.L).AND.(K.GT.J)
        GOTO 999
 5    CONTINUE
C    *Group B: Three indices equal
C       Class a1: I > J = K = L
        TST2IN = (I.GT.J).AND.(J.EQ.K).AND.(J.EQ.L)
        GOTO 999
 6    CONTINUE
C    *Group B: Three indices equal
C       Class a2: I = J = K > L
        TST2IN = (I.EQ.J).AND.(I.EQ.K).AND.(K.GT.L)
        GOTO 999
 7    CONTINUE
C    *Group B: Three indices equal
C       Class b1: I = J > K = L
        TST2IN = (I.EQ.J).AND.(J.GT.K).AND.(K.EQ.L)
        GOTO 999
 8    CONTINUE
C    *Group A: All indices equal
C       Class  1: I = J = K = L
        TST2IN = (I.EQ.J).AND.(I.EQ.K).AND.(I.EQ.L)
        GOTO 999
 9    CONTINUE
C    * XY: Class B: I > J,  K = L
        TST2IN = (I.GT.J).AND.(K.EQ.L)
        GOTO 999
 10   CONTINUE
C    * XY: Class C: I = J,  K > L
        TST2IN = (I.EQ.J).AND.(K.GT.L)
        GOTO 999
 999  CONTINUE
      RETURN
      END
