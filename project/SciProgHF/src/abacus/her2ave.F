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
C  /* Deck intexp */
      SUBROUTINE INTEXP(AOINT,DMAT,NDMAT,PSO,PSA,FT,FV,NINDAB,NINDCD,
     &                  NCCINT,NINTYP,WORK,LWORK,ISYMR,ISYMS,ISYMT,
     &                  ICORBA,ICORBB,ICORBC,ICORBD,THRESH,SYMFAC,
     &                  IPRINT,NOPMAT,NODV,EXPECT,LONDON,SUSCEP,DDFOCK,
     &                  DINTSKP,GENCNT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      LOGICAL NOPMAT, NODV, LONDON, SUSCEP, DDFOCK, EXPECT, GENCNT
      DIMENSION PSO(*), PSA(*), DMAT(*), AOINT(NCCINT,NINTYP),
     &          ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBC),
     &          WORK(LWORK), NINDAB(*), NINDCD(*), FT(*), FV(*), 
     &          DINTSKP(*)
#include "twocom.h"
#include "twosta.h"
C
C     Allocations
C
Cjth - something weird happens if KPAA etc. are undefined!
Cjth - they are now defined always.
      KDERIV = 1
      KPAO   = KDERIV + NINTYP
      IF (SUSCEP) THEN
         KPAA   = KPAO + NORBA*NORBB*NORBC*NORBD
         KLAST  = KPAA + NORBA*NORBB*NORBC*NORBD
      ELSE
         KPAA   = KPAO
         KLAST  = KPAO + NORBA*NORBB*NORBC*NORBD
      END IF
      IF (LONDON .AND. DDFOCK) THEN
         KPINT = KLAST
         KQINT = KPINT + 3*NCCINT
Cajt six components if second derivative integrals
         IF (SUSCEP) KQINT = KPINT + 6*NCCINT
         KLAST = KQINT + 3*NCCINT
Cajt six components if second derivative integrals
         IF (SUSCEP) KLAST = KQINT + 6*NCCINT
         IF (KLAST .GT. LWORK) CALL STOPIT('INTEXP',' ',KLAST,LWORK)
         LWRK   = LWORK  - KLAST + 1
         MWINTE = MAX(MWINTE,KLAST)
         LWTOT  = LWTOT + KLAST
         MWTOT  = MAX(MWTOT,LWTOT)
         CALL INTFCL(AOINT,DMAT,NDMAT,PSO,PSA,FT,FV,WORK(KDERIV),
     &            WORK(KPAO),WORK(KPAA),WORK(KPINT),
     &            WORK(KQINT),WORK(KLAST),LWRK,
     &            NCCINT,NINTYP,ISYMR,ISYMS,ISYMT,
     &            ICORBA,ICORBB,ICORBC,ICORBD,THRESH,
     &            SYMFAC,IPRINT,NOPMAT,NODV,NUCABQ,NUCCDQ,
     &            NINDAB,NINDCD,EXPECT,LONDON,SUSCEP,DDFOCK,
     &            DINTSKP)
CMI       ... I try this ...
         LWTOT  = LWTOT - KLAST
      ELSE
         KPINT = KLAST
         KQINT = KLAST
      END IF
      IF (KLAST .GT. LWORK) CALL STOPIT('INTEXP',' ',KLAST,LWORK)
      LWRK   = LWORK  - KLAST + 1
      MWINTE = MAX(MWINTE,KLAST)
      LWTOT  = LWTOT + KLAST
      MWTOT  = MAX(MWTOT,LWTOT)
      CALL INTEX1(AOINT,DMAT,NDMAT,PSO,PSA,FT,FV,WORK(KDERIV),
     &            WORK(KPAO),WORK(KPAA),WORK(KPINT),
     &            WORK(KQINT),WORK(KLAST),LWRK,
     &            NCCINT,NINTYP,ISYMR,ISYMS,ISYMT,
     &            ICORBA,ICORBB,ICORBC,ICORBD,THRESH,
     &            SYMFAC,IPRINT,NOPMAT,NODV,NUCABQ,NUCCDQ,
     &            NINDAB,NINDCD,EXPECT,LONDON,SUSCEP,DDFOCK,
     &            DINTSKP,GENCNT)
      LWTOT  = LWTOT - KLAST
      RETURN
      END
C  /* Deck intex1 */
      SUBROUTINE INTEX1(AOINT,DMAT,NDMAT,PSO,PSA,FT,FV,DERIV,PAO,PAA,
     &                  PINT,QINT,WORK,LWORK,NCCINT,NINTYP,
     &                  ISYMR,ISYMS,ISYMT,ICORBA,ICORBB,ICORBC,ICORBD,
     &                  THRESH,SYMFAC,IPRINT,NOPMAT,NODV,NUCABQ,NUCCDQ,
     &                  NINDAB,NINDCD,EXPECT,LONDON,SUSCEP,
     &                  DDFOCK,DINTSKP,GENCNT)
C
C     Calculates expectation values and Fock matrices
C     of differentiated integrals
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dcbham.h"
      PARAMETER (DP25 = 0.25 D00, DP5 = 0.5 D00, 
     &           D1 = 1.0 D00, D2 = 2.0 D00, ZERADD = 1.D-15,
     &           D0 = 0.0 D00, DP125 = 0.125D00,
     &           D4 = 4.0 D00, D8 = 8.0 D00)
      INTEGER A, B, C, D, X, Y, Z, Y2, Z2
      LOGICAL NODER, DCMPAB, DCMPCD, NOPMAT, NODV, DV, SUSCEP, DDFOCK,
     &        EXPECT, LONDON, DOQUAT, GENCNT
Cjth - bugfix ICORBD(NORBC)  -->  ICORBD(NORBD)
C      DIMENSION ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBC),
      DIMENSION ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBD),
     &          DERIV(NINTYP), PSO(*), PSA(*),
     &          PAO(NORBA,NORBB,NORBC,NORBD),
     &          PAA(NORBA,NORBB,NORBC,NORBD),
     &          AOINT(NCCINT,NINTYP),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          DMAT(NBASIS,NBASIS,NDMAT), WORK(LWORK),
     &          FT(NBASIS,NBASIS,*), FV(NBASIS,NBASIS,*), DINT(12),
     &          NCART(12), PINT(NCCINT,*), QINT(NCCINT,*),
     &          DIFAB(3), DIFCD(3), DINTSKP(*)
      DIMENSION DINTMAX(3)
#include "twocom.h"
#include "symmet.h"
#include "nuclei.h"
#include "expcom.h"
#include "dirprt.h"
#include "dcbgrd.h"
#include "doxyz.h"
Cjth
C
C     statement functions:
      XAND(I) = PT(IAND(ISYMAX(1,1),I))
      YAND(I) = PT(IAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IAND(ISYMAX(3,1),I))
      NEXT(I) = MOD(I,3) + 1
CMI  define calculation of SO contributions
      DOQUAT = .NOT.SPINFR
C
      CALL DZERO(DINTMAX,3) 
      DV = .NOT. NODV
      IF (IPRINT .GT. 9) THEN
         CALL HEADER('Output from INTEXP',-1)
         WRITE (LUPRI, '(A,4L5)') ' NODV, NOPMAT, SUSCEP, DDFOCK',
     &                              NODV, NOPMAT, SUSCEP, DDFOCK
         WRITE (LUPRI, '(A,3I5)') ' ISYMR/S/T  ', ISYMR,ISYMS,ISYMT
         WRITE (LUPRI, '(A,4I5)') ' NORB ', NORBA,NORBB,NORBC,NORBD
         WRITE (LUPRI, '(A,4I5)') ' ICENT ',ICENT1,ICENT2,ICENT3,ICENT4
         WRITE (LUPRI, '(A,2L5)') ' DIAGAB/CD ', DIAGAB,DIAGCD
         WRITE (LUPRI, '(A,3L5)') ' SHAEQB, SHCEQD, SHABAB ',
     &                              SHAEQB, SHCEQD, SHABAB
         WRITE (LUPRI, '(A,4I5)') ' NHKTA', NHKTA,NHKTB,NHKTC,NHKTD
         WRITE (LUPRI, '(A,4I5)') ' KHKTA', KHKTA,KHKTB,KHKTC,KHKTD
         WRITE (LUPRI, '(A,F12.6)') ' THRESH ', THRESH
         WRITE (LUPRI, '(A,F12.6)') ' SYMFAC ', SYMFAC
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals A '
         WRITE (LUPRI, '(20I5)') (ICORBA(I),I=1, NORBA)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals B '
         WRITE (LUPRI, '(20I5)') (ICORBB(I),I=1, NORBB)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals C '
         WRITE (LUPRI, '(20I5)') (ICORBC(I),I=1, NORBC)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals D '
         WRITE (LUPRI, '(20I5)') (ICORBD(I),I=1, NORBD)
      END IF
C
      IF (DDFOCK .AND. .NOT.LONDON) THEN
         IF (EXPECT) THEN
            IF (TWOCEN) THEN
               NCENTS = 1
               NCART(1) = 3*ICENT1 - 2
               NCART(2) = 3*ICENT1 - 1
               NCART(3) = 3*ICENT1
               NCART(4) = 3*ICENT2 - 2
               NCART(5) = 3*ICENT2 - 1
               NCART(6) = 3*ICENT2
            ELSE IF (THRCEN) THEN
               NCENTS = 2
               NCART(1) = 3*ICENT1 - 2
               NCART(2) = 3*ICENT1 - 1
               NCART(3) = 3*ICENT1
               NCART(4) = 3*ICENT2 - 2
               NCART(5) = 3*ICENT2 - 1
               NCART(6) = 3*ICENT2
               NCART(7) = 3*ICENT3 - 2
               NCART(8) = 3*ICENT3 - 1
               NCART(9) = 3*ICENT3
            ELSE
               NCENTS = 3
               NCART(1)  = 3*ICENT1 - 2
               NCART(2)  = 3*ICENT1 - 1
               NCART(3)  = 3*ICENT1
               NCART(4)  = 3*ICENT2 - 2
               NCART(5)  = 3*ICENT2 - 1
               NCART(6)  = 3*ICENT2
               NCART(7)  = 3*ICENT3 - 2
               NCART(8)  = 3*ICENT3 - 1
               NCART(9)  = 3*ICENT3
               NCART(10) = 3*ICENT4 - 2
               NCART(11) = 3*ICENT4 - 1
               NCART(12) = 3*ICENT4
            END IF
            NCARTS = 3*NCENTS
            NCARTZ = NCARTS + 3
         ELSE
            NCARTS = 3*NOATMS
            NCARTZ = NCARTS
            DO 5 I = 1, NCARTZ
               ICNT = (I + 2)/3
               IATM = ICNTDR(ICNT) - NUCNUM(NCNTDR(ICNT),1) + 1
               NCART(I) = 3*IATM + MOD(I-1,3) - 2
   5        CONTINUE
         END IF
      END IF
      ISYMTS   = IEOR(ISYMT,ISYMS)
      IF (DDFOCK .AND. LONDON) THEN
         DIFAB(1) =             CORAX0 - XAND(ISYMR )*CORBX0
         DIFAB(2) =             CORAY0 - YAND(ISYMR )*CORBY0
         DIFAB(3) =             CORAZ0 - ZAND(ISYMR )*CORBZ0
         DIFCD(1) = XAND(ISYMT)*CORCX0 - XAND(ISYMTS)*CORDX0
         DIFCD(2) = YAND(ISYMT)*CORCY0 - YAND(ISYMTS)*CORDY0
         DIFCD(3) = ZAND(ISYMT)*CORCZ0 - ZAND(ISYMTS)*CORDZ0
         DO 10 N = 1, 3
            Y   = NEXT(N)
            Z   = NEXT(Y)
            ABY = DIFAB(Y)
            ABZ = DIFAB(Z)
            CDY = DIFCD(Y)
            CDZ = DIFCD(Z)
            Y2  = Y + 3
            Z2  = Z + 3
            DO 20 I = 1, NCCINT
               AOAB = ABY*AOINT(I,Z ) - ABZ*AOINT(I,Y )
               AOCD = CDY*AOINT(I,Z2) - CDZ*AOINT(I,Y2)
               PINT(I,N) = AOAB + AOCD
               QINT(I,N) = AOAB - AOCD
  20        CONTINUE
  10     CONTINUE
      END IF
C
      SFAC = SYMFAC
      IF (.NOT.SHABAB) SFAC = D2*SFAC
      IF (.NOT.SHAEQB) SFAC = D2*SFAC
      IF (.NOT.SHCEQD) SFAC = D2*SFAC
C
      THRSH = MAX(THRESH,ZERADD)
C
      IF(DOSCRN .AND. .NOT.DDFOCK) THEN
C     ... DNSMAX fits EXPECT, not DDFOCK
         EXPTOL = MAX((SCRTHR/DNSMAX),1.00D-15)
      ELSE
         EXPTOL = 1.00D-15
      ENDIF
C
      NODER = .TRUE.
      CALL DZERO(DERIV,NINTYP)
C     initialize PAOVAL to zero for output in NOPMAT case
      PAOVAL = D0
C
C     ***** Loop over shell components *****
C
      IAOFF = 1
      DO 100 ICOMPA = 1, KHKTA
         KHKTBB = KHKTB
         IF (DIAGAB) KHKTBB = ICOMPA
         DO 200 ICOMPB = 1, KHKTBB
            DCMPAB = DIAGAB .AND. ICOMPA .EQ. ICOMPB
            FACAB  = D1
            IF (DIAGAB .AND. ICOMPA .NE. ICOMPB) FACAB = D2*FACAB
            DO 300 ICOMPC = 1, KHKTC
               KHKTDD = KHKTD
               IF (DIAGCD) KHKTDD = ICOMPC
               DO 400 ICOMPD = 1, KHKTDD
C
C                 Step 2 screening on gradient, hessian, London, and
C                 susceptibilities.
C                 Always do Step 2 screening 
C
                  AOMAX = D0
                  DO I = 1,NINTYP
                     IAOMAX = IDAMAX(NOABCD,AOINT(IAOFF,I),1)-1
                     AOMAX = MAX(AOMAX,ABS(AOINT(IAOFF+IAOMAX,I)))
                  END DO
                  IF (AOMAX .LT. EXPTOL) THEN
C
C                    If geometrical dist. then sum up statistics
C
                     IF (.NOT. LONDON .AND. EXPECT) THEN
                        DO IXYZ = 1,3
                           IF (DOXYZ(IXYZ))
     &                        CALL ADDDINT(DINTSKP,NOABCD,ICEN,3,IXYZ)
                        END DO
                     END IF
                     GO TO 510
                  END IF
C
C
                  DCMPCD = DIAGCD .AND. ICOMPC .EQ. ICOMPD
                  FACCD  = D1
                  IF (DIAGCD .AND. ICOMPC .NE. ICOMPD) FACCD = D2*FACCD
                  FCABCD = FACAB*FACCD*SFAC
C
C                 Transform P-matrix block from SO basis to AO basis
C
                  IF (.NOT.NOPMAT) THEN
                     CALL PBLOCK(PSO,PAO,ICOMPA,ICOMPB,ICOMPC,ICOMPD,
     &                           NHKTA,NHKTB,NHKTC,NHKTD,
     &                           KHKTA,KHKTB,KHKTC,MULA,MULB,MULC,MULD,
     &                           NORBA,NORBB,NORBC,NORBD,ISYMR,ISYMS,
     &                           ISYMT)
                  END IF
                  IF (.NOT.NOPMAT .AND. SUSCEP) THEN
                     CALL PBLOCK(PSA,PAA,ICOMPA,ICOMPB,ICOMPC,ICOMPD,
     &                           NHKTA,NHKTB,NHKTC,NHKTD,
     &                           KHKTA,KHKTB,KHKTC,MULA,MULB,MULC,MULD,
     &                           NORBA,NORBB,NORBC,NORBD,ISYMR,ISYMS,
     &                           ISYMT)
                  END IF
                  IF (IPRINT .GT. 10) THEN
                     WRITE (LUPRI,'(A,I10)') ' IAOFF ', IAOFF
                     WRITE (LUPRI,'(A,4I5)') ' ICOMP ',
     *                                       ICOMPA,ICOMPB,ICOMPC,ICOMPD
                     WRITE (LUPRI,'(A,2L5)') ' DCMPAB/CD ',DCMPAB,DCMPCD
                  END IF
C
                  INT = IAOFF
C
C
C                 ***** Loop over contracted functions *****
C
                  DO 500 IORBAB = 1, NORBAB
                     IORBA = NINDAB(IORBAB,1)
                     IORBB = NINDAB(IORBAB,2)
                     A = ICORBA(IORBA) + ICOMPA
                     B = ICORBB(IORBB) + ICOMPB
#ifndef PRG_DIRAC
                     DTAB = DMAT(A,B,1)
                     IF (DV) DVAB = DMAT(A,B,2)
#else
                     DTAB1 = DMAT(A,B,1)
C                    DTAB2 = DMAT(A,B,2)
C                    DTAB3 = DMAT(A,B,3)
C                    DTAB4 = DMAT(A,B,4)
                     IF (DV) THEN
                        DVAB1 = DMAT(A,B,5)
C                       DVAB2 = DMAT(A,B,6)
C                       DVAB3 = DMAT(A,B,7)
C                       DVAB4 = DMAT(A,B,8)
                     END IF
#endif
                                                      FAB = D1
                     IF (TCONAB .AND. IORBA.NE.IORBB) FAB = D2
                     DO 600 IORBCD = 1, NORBCD
                        IORBC = NINDCD(IORBCD,1)
                        IORBD = NINDCD(IORBCD,2)
                        C = ICORBC(IORBC) + ICOMPC
                        D = ICORBD(IORBD) + ICOMPD
#ifndef PRG_DIRAC
                        DTAC  = DMAT(A,C,1)
                        DTBC  = DMAT(B,C,1)
                        DTAD  = DMAT(A,D,1)
                        DTCD  = DMAT(C,D,1)
                        DTBD  = DMAT(B,D,1)
                        IF (DV) THEN
                           DVAC = DMAT(A,C,2)
                           DVBC = DMAT(B,C,2)
                           DVAD = DMAT(A,D,2)
                           DVBD = DMAT(B,D,2)
                           DVCD = DMAT(C,D,2)
                        END IF
#endif
                                                       FCD = D1
                        IF (TCONCD.AND.IORBC.NE.IORBD) FCD = D2
C
C                       ***********************************
C                       ***** Geometrical distortions *****
C                       ***********************************
C
                        IF (.NOT.LONDON) THEN
C
C                          Expectation values
C                          ==================
C
                           IF (EXPECT) THEN
#ifdef PRG_DIRAC
C
C                             For DMAT(1) include both coulombic and 
C                             exchange contributions
C                             For DMAT(2) to DMAT(4) include only
C                             exchange contributions
C
C                             For the sake of optimization only
C                             DMAT(A,B) is put in a variable.
C                             the others are not.
C
C ...try save comp. time by specifying spin-free/full relativity case
C
                              PVAL = DTAB1*DMAT(C,D,1)
                              XVAL = DMAT(C,A,1)*DMAT(D,B,1)
     &                             + DMAT(D,A,1)*DMAT(C,B,1)
                              IF (DOQUAT) THEN
                                 XVAL = XVAL
     &                                + DMAT(C,A,2)*DMAT(D,B,2)
     &                                + DMAT(D,A,2)*DMAT(C,B,2)
     &                                + DMAT(C,A,3)*DMAT(D,B,3)
     &                                + DMAT(D,A,3)*DMAT(C,B,3)
     &                                + DMAT(C,A,4)*DMAT(D,B,4)
     &                                + DMAT(D,A,4)*DMAT(C,B,4)
                              END IF
                              IF (DV) THEN
C                             ... subtract the act-act part, treated with PAO
                                 PVAL = PVAL - DVAB1*DMAT(C,D,5)
                                 XVAL = XVAL
     &                                - DMAT(C,A,5)*DMAT(D,B,5)
     &                                - DMAT(D,A,5)*DMAT(C,B,5)
                                IF (DOQUAT) THEN
                                 XVAL = XVAL
     &                                - DMAT(C,A,6)*DMAT(D,B,6)
     &                                - DMAT(D,A,6)*DMAT(C,B,6)
     &                                - DMAT(C,A,7)*DMAT(D,B,7)
     &                                - DMAT(D,A,7)*DMAT(C,B,7)
     &                                - DMAT(C,A,8)*DMAT(D,B,8)
     &                                - DMAT(D,A,8)*DMAT(C,B,8)
                                END IF
                              END IF

                              if (gencnt) then

!                               radovan: for the moment noddy,
!                               later i will make it more efficient

!                               there are two matrices now
!                               and we don't assume that the matrices are symmetric
!                               first matrix is on 1-4, second 5-8
!                               without SOC, operate only on 1 and 5

                                pval = (dmat(a, b, 1) + dmat(b, a, 1))
     &                                *(dmat(c, d, 5) + dmat(d, c, 5))
     &                               + (dmat(a, b, 5) + dmat(b, a, 5))
     &                                *(dmat(c, d, 1) + dmat(d, c, 1))

                                xval = dmat(d, b, 1)*dmat(a, c, 5)
     &                               + dmat(b, d, 1)*dmat(c, a, 5)
     &                               + dmat(c, b, 1)*dmat(a, d, 5)
     &                               + dmat(b, c, 1)*dmat(d, a, 5)
     &                               + dmat(d, a, 1)*dmat(b, c, 5)
     &                               + dmat(a, d, 1)*dmat(c, b, 5)
     &                               + dmat(c, a, 1)*dmat(b, d, 5)
     &                               + dmat(a, c, 1)*dmat(d, b, 5)

!                               SOC situation, only exchange part
                                if (doquat) then
                                xval = xval
!                                      2-6
     &                               - dmat(d, b, 2)*dmat(a, c, 6)
     &                               - dmat(b, d, 2)*dmat(c, a, 6)
     &                               - dmat(c, b, 2)*dmat(a, d, 6)
     &                               - dmat(b, c, 2)*dmat(d, a, 6)
     &                               - dmat(d, a, 2)*dmat(b, c, 6)
     &                               - dmat(a, d, 2)*dmat(c, b, 6)
     &                               - dmat(c, a, 2)*dmat(b, d, 6)
     &                               - dmat(a, c, 2)*dmat(d, b, 6)
!                                      3-7
     &                               - dmat(d, b, 3)*dmat(a, c, 7)
     &                               - dmat(b, d, 3)*dmat(c, a, 7)
     &                               - dmat(c, b, 3)*dmat(a, d, 7)
     &                               - dmat(b, c, 3)*dmat(d, a, 7)
     &                               - dmat(d, a, 3)*dmat(b, c, 7)
     &                               - dmat(a, d, 3)*dmat(c, b, 7)
     &                               - dmat(c, a, 3)*dmat(b, d, 7)
     &                               - dmat(a, c, 3)*dmat(d, b, 7)
!                                      4-8
     &                               - dmat(d, b, 4)*dmat(a, c, 8)
     &                               - dmat(b, d, 4)*dmat(c, a, 8)
     &                               - dmat(c, b, 4)*dmat(a, d, 8)
     &                               - dmat(b, c, 4)*dmat(d, a, 8)
     &                               - dmat(d, a, 4)*dmat(b, c, 8)
     &                               - dmat(a, d, 4)*dmat(c, b, 8)
     &                               - dmat(c, a, 4)*dmat(b, d, 8)
     &                               - dmat(a, c, 4)*dmat(d, b, 8)
                                end if

                                pval = dp125*pval
                                xval =  dp25*xval

                              end if !if (gencnt)

!TODO:                        lots of work above can be skipped if hfxfac is zero
!                             at the moment the zero is painfully calculated
                              PVAL = PVAL - DP25*XVAL*hfxfac
                              IF (.NOT. NOPMAT) THEN
                                 PAOVAL = PAO(IORBA,IORBB,IORBC,IORBD)
                                 PVAL = PVAL + PAOVAL
                              END IF
#endif /* ifdef PRG_DIRAC */
                              PVAL = DP5*FAB*FCD*FCABCD*PVAL
                              IF (ABS(PVAL) .GT. THRSH) THEN
                                 NODER = .FALSE.
                                 DO 700 I = 1, NINTYP
                                    DERIV(I)=DERIV(I)+PVAL*AOINT(INT,I)
  700                            CONTINUE
                              END IF
                           END IF
C
C                          Fock matrices
C                          =============
C
                           IF (DDFOCK) THEN
#ifdef PRG_DIRAC
                            CALL QUIT('DDFOCK not implemented in Dirac')
#else
                              CALL DCOPY(NCARTS,AOINT(INT,1),NCCINT,
     &                                   DINT,1)
                              IF (EXPECT) THEN
                                 DINT(NCARTS+1)=-DSUM(NCENTS,DINT(1),3)
                                 DINT(NCARTS+2)=-DSUM(NCENTS,DINT(2),3)
                                 DINT(NCARTS+3)=-DSUM(NCENTS,DINT(3),3)
                              END IF
C
C                         950210-hjaaj: new SKLFC1 symmetrizes with
C                         F(i,j) = (1/4) * (FMAT(i,j) + FMAT(j,i))
C
C old code/950210-hjaaj       DFAC = DP5*FAB*FCD*FCABCD
                              DFAC = FAB*FCD*FCABCD
                              EFAC = - DP25*DFAC
C
C                             Total Fock matrix
C
                              DCD  = DFAC*DTCD
                              DAB  = DFAC*DTAB
                              DBD  = EFAC*DTBD
                              DBC  = EFAC*DTBC
                              DAD  = EFAC*DTAD
                              DAC  = EFAC*DTAC
                              DO 800 I = 1, NCARTZ
                                 N = NCART(I)
                                 FT(A,B,N) = FT(A,B,N) + DINT(I)*DCD
                                 FT(C,D,N) = FT(C,D,N) + DINT(I)*DAB
                                 FT(A,C,N) = FT(A,C,N) + DINT(I)*DBD
                                 FT(A,D,N) = FT(A,D,N) + DINT(I)*DBC
                                 FT(B,C,N) = FT(B,C,N) + DINT(I)*DAD
                                 FT(B,D,N) = FT(B,D,N) + DINT(I)*DAC
  800                         CONTINUE
C
C                             Active Fock matrix
C
                              IF (DV) THEN
                                 DCD  = DFAC*DVCD
                                 DAB  = DFAC*DVAB
                                 DBD  = EFAC*DVBD
                                 DBC  = EFAC*DVBC
                                 DAD  = EFAC*DVAD
                                 DAC  = EFAC*DVAC
                                 DO 810 I = 1, NCARTZ
                                    N = NCART(I)
                                    FV(A,B,N) = FV(A,B,N) + DINT(I)*DCD
                                    FV(C,D,N) = FV(C,D,N) + DINT(I)*DAB
                                    FV(A,C,N) = FV(A,C,N) + DINT(I)*DBD
                                    FV(A,D,N) = FV(A,D,N) + DINT(I)*DBC
                                    FV(B,C,N) = FV(B,C,N) + DINT(I)*DAD
                                    FV(B,D,N) = FV(B,D,N) + DINT(I)*DAC
  810                            CONTINUE
                              END IF
                           END IF
                           IF (IPRINT .GT. 30) THEN
                              WRITE(LUPRI,'(//A,4I5,D12.4)')
     &                              ' IORBA... PAOVAL',
     &                              IORBA,IORBB,IORBC,IORBD,PAOVAL
                              WRITE(LUPRI,'(//A,4I5)')'A,B,C,D',A,B,C,D
                              WRITE(LUPRI,'(//A,6F12.6)') 'DTAB,...',
     &                              DTAB,DTAC,DTAD,DTBC,DTBD,DTCD
                              WRITE(LUPRI,'(//A,6F12.6)') 'DVAB,...',
     &                              DVAB,DVAC,DVAD,DVBC,DVBD,DVCD
                              WRITE (LUPRI, '(A,3F12.6)')
     &                              ' Fab, Fcd, FCabab',FAB, FCD, FCABCD
                              WRITE (LUPRI, '(A,F12.6)')' PVAL ',PVAL
                              IF (DDFOCK) THEN
                                 WRITE (LUPRI,'(A,12F6.3)') ' DINT ',
     &                               (DINT(I),I=1,NCARTZ)
                              END IF
#endif
                           END IF
C
C                       **************************
C                       ***** Magnetic field *****
C                       **************************
C
CMI/HJ (TEC) implementation of two-electron part of magnet. susceptibilities
                        ELSE
#ifdef PRG_DIRAC

!radovan: added gencnt stuff plus partial deobfuscation of the code

                          if (suscep .or. gencnt) then

                            if (.not. nodv) 
     &                      call quit('dv suscep not implemented')

                            if (suscep .and. .not. gencnt) then

                              ttcadb = dmat(c, a, 1)*dmat(d, b, 1)
                              ttdacb = dmat(d, a, 1)*dmat(c, b, 1)
                              if (doquat) then
                                ttcadb = ttcadb
     &                                 + dmat(c, a, 2)*dmat(d, b, 2)
     &                                 + dmat(c, a, 3)*dmat(d, b, 3)
     &                                 + dmat(c, a, 4)*dmat(d, b, 4)
                                ttdacb = ttdacb 
     &                                 + dmat(d, a, 2)*dmat(c, b, 2)
     &                                 + dmat(d, a, 3)*dmat(c, b, 3)
     &                                 + dmat(d, a, 4)*dmat(c, b, 4)
                              end if

                              pval1 = dmat(a, b, 1)*dmat(c, d, 1) 
     &                              - dp25*hfxfac*(ttcadb + ttdacb)
                              pval2 = dp25*hfxfac*(ttcadb - ttdacb)

                            end if

                            if (gencnt .and. .not. suscep) then

!                             coulomb part
                              pval1_c = (dmat(a, b, 1) - dmat(b, a, 1))
     &                                 *(dmat(c, d, 5) + dmat(d, c, 5))
     &                                + (dmat(a, b, 5) - dmat(b, a, 5))
     &                                 *(dmat(c, d, 1) + dmat(d, c, 1))

!                             exchange part
                              pval1_x = dmat(d, b, 1)*dmat(a, c, 5)
     &                                - dmat(b, d, 1)*dmat(c, a, 5)
     &                                + dmat(c, b, 1)*dmat(a, d, 5)
     &                                - dmat(b, c, 1)*dmat(d, a, 5)
     &                                - dmat(d, a, 1)*dmat(b, c, 5)
     &                                + dmat(a, d, 1)*dmat(c, b, 5)
     &                                - dmat(c, a, 1)*dmat(b, d, 5)
     &                                + dmat(a, c, 1)*dmat(d, b, 5)

!                             SOC situation, only exchange part
                              if (doquat) then
                              pval1_x = pval1_x
!                                       2-6
     &                                + dmat(d, b, 2)*dmat(a, c, 6)
     &                                - dmat(b, d, 2)*dmat(c, a, 6)
     &                                + dmat(c, b, 2)*dmat(a, d, 6)
     &                                - dmat(b, c, 2)*dmat(d, a, 6)
     &                                - dmat(d, a, 2)*dmat(b, c, 6)
     &                                + dmat(a, d, 2)*dmat(c, b, 6)
     &                                - dmat(c, a, 2)*dmat(b, d, 6)
     &                                + dmat(a, c, 2)*dmat(d, b, 6)
!                                       3-7
     &                                + dmat(d, b, 3)*dmat(a, c, 7)
     &                                - dmat(b, d, 3)*dmat(c, a, 7)
     &                                + dmat(c, b, 3)*dmat(a, d, 7)
     &                                - dmat(b, c, 3)*dmat(d, a, 7)
     &                                - dmat(d, a, 3)*dmat(b, c, 7)
     &                                + dmat(a, d, 3)*dmat(c, b, 7)
     &                                - dmat(c, a, 3)*dmat(b, d, 7)
     &                                + dmat(a, c, 3)*dmat(d, b, 7)
!                                       4-8
     &                                + dmat(d, b, 4)*dmat(a, c, 8)
     &                                - dmat(b, d, 4)*dmat(c, a, 8)
     &                                + dmat(c, b, 4)*dmat(a, d, 8)
     &                                - dmat(b, c, 4)*dmat(d, a, 8)
     &                                - dmat(d, a, 4)*dmat(b, c, 8)
     &                                + dmat(a, d, 4)*dmat(c, b, 8)
     &                                - dmat(c, a, 4)*dmat(b, d, 8)
     &                                + dmat(a, c, 4)*dmat(d, b, 8)
                              end if

!                             coulomb part (from pval1_c permute signs between dmats)
                              pval2_c = (dmat(a, b, 1) + dmat(b, a, 1))
     &                                 *(dmat(c, d, 5) - dmat(d, c, 5))
     &                                + (dmat(a, b, 5) + dmat(b, a, 5))
     &                                 *(dmat(c, d, 1) - dmat(d, c, 1))

!                             exchange part (from pval1_x permute some signs: 1, 2, 7, 8)
                              pval2_x =
     &                                - dmat(d, b, 1)*dmat(a, c, 5)
     &                                + dmat(b, d, 1)*dmat(c, a, 5)
     &                                + dmat(c, b, 1)*dmat(a, d, 5)
     &                                - dmat(b, c, 1)*dmat(d, a, 5)
     &                                - dmat(d, a, 1)*dmat(b, c, 5)
     &                                + dmat(a, d, 1)*dmat(c, b, 5)
     &                                + dmat(c, a, 1)*dmat(b, d, 5)
     &                                - dmat(a, c, 1)*dmat(d, b, 5)

!                             SOC situation, only exchange part
                              if (doquat) then
                              pval2_x = pval2_x
!                                       2-6
     &                                - dmat(d, b, 2)*dmat(a, c, 6)
     &                                + dmat(b, d, 2)*dmat(c, a, 6)
     &                                + dmat(c, b, 2)*dmat(a, d, 6)
     &                                - dmat(b, c, 2)*dmat(d, a, 6)
     &                                - dmat(d, a, 2)*dmat(b, c, 6)
     &                                + dmat(a, d, 2)*dmat(c, b, 6)
     &                                + dmat(c, a, 2)*dmat(b, d, 6)
     &                                - dmat(a, c, 2)*dmat(d, b, 6)
!                                       3-7
     &                                - dmat(d, b, 3)*dmat(a, c, 7)
     &                                + dmat(b, d, 3)*dmat(c, a, 7)
     &                                + dmat(c, b, 3)*dmat(a, d, 7)
     &                                - dmat(b, c, 3)*dmat(d, a, 7)
     &                                - dmat(d, a, 3)*dmat(b, c, 7)
     &                                + dmat(a, d, 3)*dmat(c, b, 7)
     &                                + dmat(c, a, 3)*dmat(b, d, 7)
     &                                - dmat(a, c, 3)*dmat(d, b, 7)
!                                       4-8
     &                                - dmat(d, b, 4)*dmat(a, c, 8)
     &                                + dmat(b, d, 4)*dmat(c, a, 8)
     &                                + dmat(c, b, 4)*dmat(a, d, 8)
     &                                - dmat(b, c, 4)*dmat(d, a, 8)
     &                                - dmat(d, a, 4)*dmat(b, c, 8)
     &                                + dmat(a, d, 4)*dmat(c, b, 8)
     &                                + dmat(c, a, 4)*dmat(b, d, 8)
     &                                - dmat(a, c, 4)*dmat(d, b, 8)
                              end if

                              pval1 = pval1_c - dp5*pval1_x
                              pval2 = pval2_c - dp5*pval2_x

                            end if

                            if (suscep .and. gencnt) then

!                             coulomb part
                              pval1_c = dmat(a, b, 1)*dmat(c, d, 5)
     &                                + dmat(b, a, 1)*dmat(d, c, 5)
     &                                + dmat(c, d, 1)*dmat(b, a, 5)
     &                                + dmat(d, c, 1)*dmat(a, b, 5)

!                             exchange part
                              pval1_x = dmat(a, c, 1)*dmat(b, d, 5)
     &                                + dmat(c, a, 1)*dmat(d, b, 5)
     &                                + dmat(b, d, 1)*dmat(a, c, 5)
     &                                + dmat(d, b, 1)*dmat(c, a, 5)

!                             SOC situation, only exchange part
                              if (doquat) then
                              pval1_x = pval1_x
!                                       2-6
     &                                + dmat(a, c, 2)*dmat(b, d, 6)
     &                                + dmat(c, a, 2)*dmat(d, b, 6)
     &                                + dmat(b, d, 2)*dmat(a, c, 6)
     &                                + dmat(d, b, 2)*dmat(c, a, 6)
!                                       3-7
     &                                + dmat(a, c, 3)*dmat(b, d, 7)
     &                                + dmat(c, a, 3)*dmat(d, b, 7)
     &                                + dmat(b, d, 3)*dmat(a, c, 7)
     &                                + dmat(d, b, 3)*dmat(c, a, 7)
!                                       4-8
     &                                + dmat(a, c, 4)*dmat(b, d, 8)
     &                                + dmat(c, a, 4)*dmat(d, b, 8)
     &                                + dmat(b, d, 4)*dmat(a, c, 8)
     &                                + dmat(d, b, 4)*dmat(c, a, 8)
                              end if

!gosia: hfxfac here untested
                              pval1 = hfxfac*pval1_x
                              pval1 = d2*pval1_c - hfxfac*pval1_x
                              pval2 = pval1

                            end if

                            if (.not. nopmat .and. .not. gencnt) then
                              pval1 = pval1 
     &                              + pao(iorba, iorbb, iorbc, iorbd)
                              pval2 = pval2 
     &                              + paa(iorba, iorbb, iorbc, iorbd)
                            end if

                            factor = dp125*fab*fcd*fcabcd
                            pval1  = factor*pval1
                            pval2  = factor*pval2

                            if (abs(pval1) .gt. thrsh) then
                              noder = .false.
                              do i = 7, 18
                                deriv(i) = deriv(i) 
     &                                   + pval1*aoint(int, i)
                              end do
                            end if
                            if (abs(pval2) .gt. thrsh) then
                              noder = .false.
                              do i = 19, 27
                                deriv(i) = deriv(i) 
     &                                   + pval2*aoint(int, i)
                              end do
                            end if

                          end if !if (suscep) then
#else
                           FACTOR = DP125*FAB*FCD*FCABCD
                           IF (SUSCEP) THEN
                              IF (NOPMAT) THEN
                                 PAOVL1 = D0
                                 PAOVL2 = D0
                              ELSE
                                 PAOVL1 = PAO(IORBA,IORBB,IORBC,IORBD)
                                 PAOVL2 = PAA(IORBA,IORBB,IORBC,IORBD)
                              ENDIF
                              IF (NODV) THEN
                                 PVAL1  = FACTOR*(DTAB*DTCD
     &                                  - DP25*(DTAC*DTBD + DTAD*DTBC)
     &                                  + PAOVL1)
                                 PVAL2  = FACTOR*(
     &                                    DP25*(DTAC*DTBD - DTAD*DTBC)
     &                                  + PAOVL2)
                              ELSE
                                 PVAL1  = FACTOR*
     &                                    (DTAB*DTCD - DVAB*DVCD
     &                                  - DP25*(DTAC*DTBD - DVAC*DVBD
     &                                        + DTAD*DTBC - DVAD*DVBC)
     &                                  + PAOVL1)
                                 PVAL2  = FACTOR*(
     &                                    DP25*(DTAC*DTBD - DVAC*DVBD
     &                                        - DTAD*DTBC + DVAD*DVBC)
     &                                  + PAOVL2)
                              END IF
                              IF (ABS(PVAL1) .GT. THRSH) THEN
                                 NODER = .FALSE.
                                 DO 900 I = 7, 18
                                    DERIV(I)=DERIV(I)+PVAL1*AOINT(INT,I)
  900                            CONTINUE
                              END IF
                              IF (ABS(PVAL2) .GT. THRSH) THEN
                                 NODER = .FALSE.
                                 DO 910 I = 19, 27
                                    DERIV(I)=DERIV(I)+PVAL2*AOINT(INT,I)
  910                            CONTINUE
                              END IF
                           END IF
                           IF (DDFOCK) THEN
                              DFAC = DP125*FAB*FCD*FCABCD
                              EFAC = DP5*DFAC
C
C                             Total Fock matrix
C
                              DAB  = DFAC*DTAB
                              DCD  = DFAC*DTCD
                              DAC  = EFAC*DTAC
                              DAD  = EFAC*DTAD
                              DBC  = EFAC*DTBC
                              DBD  = EFAC*DTBD
                              DO 920 N = 1, 3
                                 ABCD = PINT(INT,N)
                                 ABDC = QINT(INT,N)
                                 FT(A,B,N) = FT(A,B,N) + (ABCD+ABDC)*DCD
                                 FT(C,D,N) = FT(C,D,N) + (ABCD-ABDC)*DAB
                                 FT(A,C,N) = FT(A,C,N) - ABDC*DBD
                                 FT(A,D,N) = FT(A,D,N) - ABCD*DBC
                                 FT(B,C,N) = FT(B,C,N) + ABCD*DAD
                                 FT(B,D,N) = FT(B,D,N) + ABDC*DAC
  920                         CONTINUE
C
C                             Active Fock matrix
C
                              IF (DV) THEN
                                 DAB  = DFAC*DVAB
                                 DCD  = DFAC*DVCD
                                 DAC  = EFAC*DVAC
                                 DAD  = EFAC*DVAD
                                 DBC  = EFAC*DVBC
                                 DBD  = EFAC*DVBD
                                 DO 930 N = 1, 3
                                    ABCD = PINT(INT,N)
                                    ABDC = QINT(INT,N)
                                    FV(A,B,N) = FV(A,B,N)
     &                                        + (ABCD+ABDC)*DCD
                                    FV(C,D,N) = FV(C,D,N)
     &                                        + (ABCD-ABDC)*DAB
                                    FV(A,C,N) = FV(A,C,N) - ABDC*DBD
                                    FV(A,D,N) = FV(A,D,N) - ABCD*DBC
                                    FV(B,C,N) = FV(B,C,N) + ABCD*DAD
                                    FV(B,D,N) = FV(B,D,N) + ABDC*DAC
  930                            CONTINUE
                              END IF
                           END IF
                           IF (IPRINT .GT. 30) THEN
                              WRITE(LUPRI,'(//A,4I5,2D12.4)')
     &                             ' IORBA... PAOVAL',
     &                             IORBA,IORBB,IORBC,IORBD,PAOVL1,PAOVL2
                              WRITE(LUPRI,'(//A,4I5)')'A,B,C,D',A,B,C,D
                              WRITE(LUPRI,'(//A,6F12.6)') 'DTAB,...',
     &                              DTAB,DTAC,DTAD,DTBC,DTBD,DTCD
                              WRITE(LUPRI,'(//A,6F12.6)') 'DVAB,...',
     &                              DVAB,DVAC,DVAD,DVBC,DVBD,DVCD
                              WRITE (LUPRI, '(A,3F12.6)')
     &                              ' Fab, Fcd, FCabab',FAB, FCD, FCABCD
                              WRITE (LUPRI, '(A,2F12.6)')' PVAL ',
     &                               PVAL1,PVAL2
                              IF (DDFOCK) THEN
                                 WRITE (LUPRI,'(A,12F6.3)') ' DINT ',
     &                               (DINT(I),I=1,NCARTZ)
                              END IF
                           END IF
#endif  /* ifdef DIRAC */
                        END IF
                        INT = INT + 1
  600                CONTINUE
  500             CONTINUE
  510             IAOFF = IAOFF + NOABCD
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      IF (.NOT.NODER) THEN
         IF (.NOT.LONDON .AND. EXPECT) THEN
            IF (GENCNT) THEN
               CALL DROUTG(DERIV,FT,NINTYP,IPRINT)
            ELSE
               CALL DEROUT(DERIV,NINTYP,IPRINT)
            END IF
         END IF
         IF (LONDON .AND. GENCNT) THEN
            IF (.NOT.SUSCEP) THEN
               CALL LNDGC1(DERIV,FT,ISYMR,ISYMT,ISYMTS,IPRINT)
            ELSE IF (SUSCEP) THEN
               CALL LNDGC2(DERIV,FT,ISYMR,ISYMT,ISYMTS,IPRINT)
            END IF
         ELSE IF (LONDON .AND. SUSCEP) THEN
            CALL LNDOUT(DERIV,ISYMR,ISYMT,ISYMTS,IPRINT)
         END IF
      END IF
      RETURN
      END
C  /* Deck derout */
      SUBROUTINE DEROUT(DERIV,NINTYP,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER AX, AY, AZ, BX, BY, BZ, CX, CY, CZ, DX, DY, DZ
      DIMENSION DERIV(NINTYP)
C
#include "nuclei.h"
#include "expcom.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
C
      IF (IPRINT .GE.10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1110) TWOCEN, THRCEN
         WRITE (LUPRI, 1120) ICENT1, ICENT2, ICENT3, ICENT4
      END IF
C
C     ***************************
C     ***** Two-Center Case *****
C     ***************************
C
      IF (TWOCEN) THEN
         AX = 3*ICENT1 - 2
         AY = 3*ICENT1 - 1
         AZ = 3*ICENT1
         BX = 3*ICENT2 - 2
         BY = 3*ICENT2 - 1
         BZ = 3*ICENT2
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
C
C           A electron-repulsion gradient elements:
C
            GRADEE(AX) = GRADEE(AX) + DAX
            GRADEE(AY) = GRADEE(AY) + DAY
            GRADEE(AZ) = GRADEE(AZ) + DAZ
C
C           B electron-repulsion gradient elements:
C
            GRADEE(BX) = GRADEE(BX) - DAX
            GRADEE(BY) = GRADEE(BY) - DAY
            GRADEE(BZ) = GRADEE(BZ) - DAZ
         END IF
C
C        ***** Hessian *****
C
#ifndef PRG_DIRAC
         IF (DERTWO) THEN
C
            AAXX = DERIV(4)
            AAXY = DERIV(5)
            AAXZ = DERIV(6)
            AAYY = DERIV(7)
            AAYZ = DERIV(8)
            AAZZ = DERIV(9)
C
C           A-A electron-repulsion Hessian elements:
C
            HESSEE(AX,AX) = HESSEE(AX,AX) + AAXX
            HESSEE(AX,AY) = HESSEE(AX,AY) + AAXY
            HESSEE(AX,AZ) = HESSEE(AX,AZ) + AAXZ
            HESSEE(AY,AY) = HESSEE(AY,AY) + AAYY
            HESSEE(AY,AZ) = HESSEE(AY,AZ) + AAYZ
            HESSEE(AZ,AZ) = HESSEE(AZ,AZ) + AAZZ
C
C           A-B electron-repulsion Hessian elements:
C
            HESSEE(AX,BX) = HESSEE(AX,BX) - AAXX
            HESSEE(AX,BY) = HESSEE(AX,BY) - AAXY
            HESSEE(AX,BZ) = HESSEE(AX,BZ) - AAXZ
            HESSEE(AY,BX) = HESSEE(AY,BX) - AAXY
            HESSEE(AY,BY) = HESSEE(AY,BY) - AAYY
            HESSEE(AY,BZ) = HESSEE(AY,BZ) - AAYZ
            HESSEE(AZ,BX) = HESSEE(AZ,BX) - AAXZ
            HESSEE(AZ,BY) = HESSEE(AZ,BY) - AAYZ
            HESSEE(AZ,BZ) = HESSEE(AZ,BZ) - AAZZ
C
C           B-B electron-repulsion Hessian elements:
C
            HESSEE(BX,BX) = HESSEE(BX,BX) + AAXX
            HESSEE(BX,BY) = HESSEE(BX,BY) + AAXY
            HESSEE(BX,BZ) = HESSEE(BX,BZ) + AAXZ
            HESSEE(BY,BY) = HESSEE(BY,BY) + AAYY
            HESSEE(BY,BZ) = HESSEE(BY,BZ) + AAYZ
            HESSEE(BZ,BZ) = HESSEE(BZ,BZ) + AAZZ
         END IF
#endif
C
C     *****************************
C     ***** Three-Center Case *****
C     *****************************
C
      ELSE IF (THRCEN) THEN
         AX = 3*ICENT1 - 2
         AY = 3*ICENT1 - 1
         AZ = 3*ICENT1
         BX = 3*ICENT2 - 2
         BY = 3*ICENT2 - 1
         BZ = 3*ICENT2
         CX = 3*ICENT3 - 2
         CY = 3*ICENT3 - 1
         CZ = 3*ICENT3
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
            DBX = DERIV(4)
            DBY = DERIV(5)
            DBZ = DERIV(6)
C
C           A electron-repulsion gradient elements:
C
            GRADEE(AX) = GRADEE(AX) + DAX
            GRADEE(AY) = GRADEE(AY) + DAY
            GRADEE(AZ) = GRADEE(AZ) + DAZ
C
C           B electron-repulsion gradient elements:
C
            GRADEE(BX) = GRADEE(BX) + DBX
            GRADEE(BY) = GRADEE(BY) + DBY
            GRADEE(BZ) = GRADEE(BZ) + DBZ
C
C           C electron-repulsion gradient elements:
C
            GRADEE(CX) = GRADEE(CX) - (DAX + DBX)
            GRADEE(CY) = GRADEE(CY) - (DAY + DBY)
            GRADEE(CZ) = GRADEE(CZ) - (DAZ + DBZ)
         END IF
C
C        ***** Hessian *****
C
#ifndef PRG_DIRAC
         IF (DERTWO) THEN
C
            AAXX = DERIV(7)
            AAXY = DERIV(8)
            AAXZ = DERIV(9)
            AAYY = DERIV(10)
            AAYZ = DERIV(11)
            AAZZ = DERIV(12)
C
            BBXX = DERIV(13)
            BBXY = DERIV(14)
            BBXZ = DERIV(15)
            BBYY = DERIV(16)
            BBYZ = DERIV(17)
            BBZZ = DERIV(18)
C
            ABXX = DERIV(19)
            ABXY = DERIV(20)
            ABXZ = DERIV(21)
            ABYX = DERIV(22)
            ABYY = DERIV(23)
            ABYZ = DERIV(24)
            ABZX = DERIV(25)
            ABZY = DERIV(26)
            ABZZ = DERIV(27)
C
C           A-A electron-repulsion Hessian elements:
C
            HESSEE(AX,AX) = HESSEE(AX,AX) + AAXX
            HESSEE(AX,AY) = HESSEE(AX,AY) + AAXY
            HESSEE(AX,AZ) = HESSEE(AX,AZ) + AAXZ
            HESSEE(AY,AY) = HESSEE(AY,AY) + AAYY
            HESSEE(AY,AZ) = HESSEE(AY,AZ) + AAYZ
            HESSEE(AZ,AZ) = HESSEE(AZ,AZ) + AAZZ
C
C           A-B electron-repulsion Hessian elements:
C
            HESSEE(AX,BX) = HESSEE(AX,BX) + ABXX
            HESSEE(AX,BY) = HESSEE(AX,BY) + ABXY
            HESSEE(AX,BZ) = HESSEE(AX,BZ) + ABXZ
            HESSEE(AY,BX) = HESSEE(AY,BX) + ABYX
            HESSEE(AY,BY) = HESSEE(AY,BY) + ABYY
            HESSEE(AY,BZ) = HESSEE(AY,BZ) + ABYZ
            HESSEE(AZ,BX) = HESSEE(AZ,BX) + ABZX
            HESSEE(AZ,BY) = HESSEE(AZ,BY) + ABZY
            HESSEE(AZ,BZ) = HESSEE(AZ,BZ) + ABZZ
C
C           A-C electron-repulsion Hessian elements:
C
            HESSEE(AX,CX) = HESSEE(AX,CX) - (AAXX + ABXX)
            HESSEE(AX,CY) = HESSEE(AX,CY) - (AAXY + ABXY)
            HESSEE(AX,CZ) = HESSEE(AX,CZ) - (AAXZ + ABXZ)
            HESSEE(AY,CX) = HESSEE(AY,CX) - (AAXY + ABYX)
            HESSEE(AY,CY) = HESSEE(AY,CY) - (AAYY + ABYY)
            HESSEE(AY,CZ) = HESSEE(AY,CZ) - (AAYZ + ABYZ)
            HESSEE(AZ,CX) = HESSEE(AZ,CX) - (AAXZ + ABZX)
            HESSEE(AZ,CY) = HESSEE(AZ,CY) - (AAYZ + ABZY)
            HESSEE(AZ,CZ) = HESSEE(AZ,CZ) - (AAZZ + ABZZ)
C
C           B-B electron-repulsion Hessian elements:
C
            HESSEE(BX,BX) = HESSEE(BX,BX) + BBXX
            HESSEE(BX,BY) = HESSEE(BX,BY) + BBXY
            HESSEE(BX,BZ) = HESSEE(BX,BZ) + BBXZ
            HESSEE(BY,BY) = HESSEE(BY,BY) + BBYY
            HESSEE(BY,BZ) = HESSEE(BY,BZ) + BBYZ
            HESSEE(BZ,BZ) = HESSEE(BZ,BZ) + BBZZ
C
C           B-C electron-repulsion Hessian elements:
C
            HESSEE(BX,CX) = HESSEE(BX,CX) - (BBXX + ABXX)
            HESSEE(BX,CY) = HESSEE(BX,CY) - (BBXY + ABYX)
            HESSEE(BX,CZ) = HESSEE(BX,CZ) - (BBXZ + ABZX)
            HESSEE(BY,CX) = HESSEE(BY,CX) - (BBXY + ABXY)
            HESSEE(BY,CY) = HESSEE(BY,CY) - (BBYY + ABYY)
            HESSEE(BY,CZ) = HESSEE(BY,CZ) - (BBYZ + ABZY)
            HESSEE(BZ,CX) = HESSEE(BZ,CX) - (BBXZ + ABXZ)
            HESSEE(BZ,CY) = HESSEE(BZ,CY) - (BBYZ + ABYZ)
            HESSEE(BZ,CZ) = HESSEE(BZ,CZ) - (BBZZ + ABZZ)
C
C           C-C electron-repulsion Hessian elements:
C
            HESSEE(CX,CX) = HESSEE(CX,CX) + (AAXX + ABXX + ABXX + BBXX)
            HESSEE(CX,CY) = HESSEE(CX,CY) + (AAXY + ABXY + ABYX + BBXY)
            HESSEE(CX,CZ) = HESSEE(CX,CZ) + (AAXZ + ABXZ + ABZX + BBXZ)
            HESSEE(CY,CY) = HESSEE(CY,CY) + (AAYY + ABYY + ABYY + BBYY)
            HESSEE(CY,CZ) = HESSEE(CY,CZ) + (AAYZ + ABYZ + ABZY + BBYZ)
            HESSEE(CZ,CZ) = HESSEE(CZ,CZ) + (AAZZ + ABZZ + ABZZ + BBZZ)
         END IF
#endif
C
C     ****************************
C     ***** Four-Center Case *****
C     ****************************
C
      ELSE
         AX = 3*ICENT1 - 2
         AY = 3*ICENT1 - 1
         AZ = 3*ICENT1
         BX = 3*ICENT2 - 2
         BY = 3*ICENT2 - 1
         BZ = 3*ICENT2
         CX = 3*ICENT3 - 2
         CY = 3*ICENT3 - 1
         CZ = 3*ICENT3
         DX = 3*ICENT4 - 2
         DY = 3*ICENT4 - 1
         DZ = 3*ICENT4
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
            DBX = DERIV(4)
            DBY = DERIV(5)
            DBZ = DERIV(6)
            DCX = DERIV(7)
            DCY = DERIV(8)
            DCZ = DERIV(9)
C
C           A electron-repulsion gradient elements:
C
            GRADEE(AX) = GRADEE(AX) + DAX
            GRADEE(AY) = GRADEE(AY) + DAY
            GRADEE(AZ) = GRADEE(AZ) + DAZ
C
C           B electron-repulsion gradient elements:
C
            GRADEE(BX) = GRADEE(BX) + DBX
            GRADEE(BY) = GRADEE(BY) + DBY
            GRADEE(BZ) = GRADEE(BZ) + DBZ
C
C           C electron-repulsion gradient elements:
C
            GRADEE(CX) = GRADEE(CX) + DCX
            GRADEE(CY) = GRADEE(CY) + DCY
            GRADEE(CZ) = GRADEE(CZ) + DCZ
C
C           D electron-repulsion gradient elements:
C
            GRADEE(DX) = GRADEE(DX) - (DAX + DBX + DCX)
            GRADEE(DY) = GRADEE(DY) - (DAY + DBY + DCY)
            GRADEE(DZ) = GRADEE(DZ) - (DAZ + DBZ + DCZ)
         END IF
C
C        ***** Hessian *****
C
#ifndef PRG_DIRAC
         IF (DERTWO) THEN
C
            AAXX = DERIV(10)
            AAXY = DERIV(11)
            AAXZ = DERIV(12)
            AAYY = DERIV(13)
            AAYZ = DERIV(14)
            AAZZ = DERIV(15)
C
            BBXX = DERIV(16)
            BBXY = DERIV(17)
            BBXZ = DERIV(18)
            BBYY = DERIV(19)
            BBYZ = DERIV(20)
            BBZZ = DERIV(21)
C
            CCXX = DERIV(22)
            CCXY = DERIV(23)
            CCXZ = DERIV(24)
            CCYY = DERIV(25)
            CCYZ = DERIV(26)
            CCZZ = DERIV(27)
C
            ABXX = DERIV(28)
            ABXY = DERIV(29)
            ABXZ = DERIV(30)
            ABYX = DERIV(31)
            ABYY = DERIV(32)
            ABYZ = DERIV(33)
            ABZX = DERIV(34)
            ABZY = DERIV(35)
            ABZZ = DERIV(36)
C
            ACXX = DERIV(37)
            ACXY = DERIV(38)
            ACXZ = DERIV(39)
            ACYX = DERIV(40)
            ACYY = DERIV(41)
            ACYZ = DERIV(42)
            ACZX = DERIV(43)
            ACZY = DERIV(44)
            ACZZ = DERIV(45)
C
            BCXX = DERIV(46)
            BCXY = DERIV(47)
            BCXZ = DERIV(48)
            BCYX = DERIV(49)
            BCYY = DERIV(50)
            BCYZ = DERIV(51)
            BCZX = DERIV(52)
            BCZY = DERIV(53)
            BCZZ = DERIV(54)
C
C           A-A electron-repulsion Hessian elements:
C
            HESSEE(AX,AX) = HESSEE(AX,AX) + AAXX
            HESSEE(AX,AY) = HESSEE(AX,AY) + AAXY
            HESSEE(AX,AZ) = HESSEE(AX,AZ) + AAXZ
            HESSEE(AY,AY) = HESSEE(AY,AY) + AAYY
            HESSEE(AY,AZ) = HESSEE(AY,AZ) + AAYZ
            HESSEE(AZ,AZ) = HESSEE(AZ,AZ) + AAZZ
C
C           A-B electron-repulsion Hessian elements:
C
            HESSEE(AX,BX) = HESSEE(AX,BX) + ABXX
            HESSEE(AX,BY) = HESSEE(AX,BY) + ABXY
            HESSEE(AX,BZ) = HESSEE(AX,BZ) + ABXZ
            HESSEE(AY,BX) = HESSEE(AY,BX) + ABYX
            HESSEE(AY,BY) = HESSEE(AY,BY) + ABYY
            HESSEE(AY,BZ) = HESSEE(AY,BZ) + ABYZ
            HESSEE(AZ,BX) = HESSEE(AZ,BX) + ABZX
            HESSEE(AZ,BY) = HESSEE(AZ,BY) + ABZY
            HESSEE(AZ,BZ) = HESSEE(AZ,BZ) + ABZZ
C
C           A-C electron-repulsion Hessian elements:
C
            HESSEE(AX,CX) = HESSEE(AX,CX) + ACXX
            HESSEE(AX,CY) = HESSEE(AX,CY) + ACXY
            HESSEE(AX,CZ) = HESSEE(AX,CZ) + ACXZ
            HESSEE(AY,CX) = HESSEE(AY,CX) + ACYX
            HESSEE(AY,CY) = HESSEE(AY,CY) + ACYY
            HESSEE(AY,CZ) = HESSEE(AY,CZ) + ACYZ
            HESSEE(AZ,CX) = HESSEE(AZ,CX) + ACZX
            HESSEE(AZ,CY) = HESSEE(AZ,CY) + ACZY
            HESSEE(AZ,CZ) = HESSEE(AZ,CZ) + ACZZ
C
C           A-D electron-repulsion Hessian elements:
C
            HESSEE(AX,DX) = HESSEE(AX,DX) - (AAXX + ABXX + ACXX)
            HESSEE(AX,DY) = HESSEE(AX,DY) - (AAXY + ABXY + ACXY)
            HESSEE(AX,DZ) = HESSEE(AX,DZ) - (AAXZ + ABXZ + ACXZ)
            HESSEE(AY,DX) = HESSEE(AY,DX) - (AAXY + ABYX + ACYX)
            HESSEE(AY,DY) = HESSEE(AY,DY) - (AAYY + ABYY + ACYY)
            HESSEE(AY,DZ) = HESSEE(AY,DZ) - (AAYZ + ABYZ + ACYZ)
            HESSEE(AZ,DX) = HESSEE(AZ,DX) - (AAXZ + ABZX + ACZX)
            HESSEE(AZ,DY) = HESSEE(AZ,DY) - (AAYZ + ABZY + ACZY)
            HESSEE(AZ,DZ) = HESSEE(AZ,DZ) - (AAZZ + ABZZ + ACZZ)
C
C           B-B electron-repulsion Hessian elements:
C
            HESSEE(BX,BX) = HESSEE(BX,BX) + BBXX
            HESSEE(BX,BY) = HESSEE(BX,BY) + BBXY
            HESSEE(BX,BZ) = HESSEE(BX,BZ) + BBXZ
            HESSEE(BY,BY) = HESSEE(BY,BY) + BBYY
            HESSEE(BY,BZ) = HESSEE(BY,BZ) + BBYZ
            HESSEE(BZ,BZ) = HESSEE(BZ,BZ) + BBZZ
C
C           B-C electron-repulsion Hessian elements:
C
            HESSEE(BX,CX) = HESSEE(BX,CX) + BCXX
            HESSEE(BX,CY) = HESSEE(BX,CY) + BCXY
            HESSEE(BX,CZ) = HESSEE(BX,CZ) + BCXZ
            HESSEE(BY,CX) = HESSEE(BY,CX) + BCYX
            HESSEE(BY,CY) = HESSEE(BY,CY) + BCYY
            HESSEE(BY,CZ) = HESSEE(BY,CZ) + BCYZ
            HESSEE(BZ,CX) = HESSEE(BZ,CX) + BCZX
            HESSEE(BZ,CY) = HESSEE(BZ,CY) + BCZY
            HESSEE(BZ,CZ) = HESSEE(BZ,CZ) + BCZZ
C
C           B-D electron-repulsion Hessian elements:
C
            HESSEE(BX,DX) = HESSEE(BX,DX) - (ABXX + BBXX + BCXX)
            HESSEE(BX,DY) = HESSEE(BX,DY) - (ABYX + BBXY + BCXY)
            HESSEE(BX,DZ) = HESSEE(BX,DZ) - (ABZX + BBXZ + BCXZ)
            HESSEE(BY,DX) = HESSEE(BY,DX) - (ABXY + BBXY + BCYX)
            HESSEE(BY,DY) = HESSEE(BY,DY) - (ABYY + BBYY + BCYY)
            HESSEE(BY,DZ) = HESSEE(BY,DZ) - (ABZY + BBYZ + BCYZ)
            HESSEE(BZ,DX) = HESSEE(BZ,DX) - (ABXZ + BBXZ + BCZX)
            HESSEE(BZ,DY) = HESSEE(BZ,DY) - (ABYZ + BBYZ + BCZY)
            HESSEE(BZ,DZ) = HESSEE(BZ,DZ) - (ABZZ + BBZZ + BCZZ)
C
C           C-C electron-repulsion Hessian elements:
C
            HESSEE(CX,CX) = HESSEE(CX,CX) + CCXX
            HESSEE(CX,CY) = HESSEE(CX,CY) + CCXY
            HESSEE(CX,CZ) = HESSEE(CX,CZ) + CCXZ
            HESSEE(CY,CY) = HESSEE(CY,CY) + CCYY
            HESSEE(CY,CZ) = HESSEE(CY,CZ) + CCYZ
            HESSEE(CZ,CZ) = HESSEE(CZ,CZ) + CCZZ
C
C           C-D electron-repulsion Hessian elements:
C
            HESSEE(CX,DX) = HESSEE(CX,DX) - (ACXX + BCXX + CCXX)
            HESSEE(CX,DY) = HESSEE(CX,DY) - (ACYX + BCYX + CCXY)
            HESSEE(CX,DZ) = HESSEE(CX,DZ) - (ACZX + BCZX + CCXZ)
            HESSEE(CY,DX) = HESSEE(CY,DX) - (ACXY + BCXY + CCXY)
            HESSEE(CY,DY) = HESSEE(CY,DY) - (ACYY + BCYY + CCYY)
            HESSEE(CY,DZ) = HESSEE(CY,DZ) - (ACZY + BCZY + CCYZ)
            HESSEE(CZ,DX) = HESSEE(CZ,DX) - (ACXZ + BCXZ + CCXZ)
            HESSEE(CZ,DY) = HESSEE(CZ,DY) - (ACYZ + BCYZ + CCYZ)
            HESSEE(CZ,DZ) = HESSEE(CZ,DZ) - (ACZZ + BCZZ + CCZZ)
C
C           D-D electron-repulsion Hessian elements:
C
            HESSEE(DX,DX) = HESSEE(DX,DX) + (AAXX + ABXX + ACXX
     &                + ABXX + BBXX + BCXX + ACXX + BCXX + CCXX)
            HESSEE(DX,DY) = HESSEE(DX,DY) + (AAXY + ABXY + ACXY
     &                + ABYX + BBXY + BCXY + ACYX + BCYX + CCXY)
            HESSEE(DX,DZ) = HESSEE(DX,DZ) + (AAXZ + ABXZ + ACXZ
     &                + ABZX + BBXZ + BCXZ + ACZX + BCZX + CCXZ)
            HESSEE(DY,DY) = HESSEE(DY,DY) + (AAYY + ABYY + ACYY
     &                + ABYY + BBYY + BCYY + ACYY + BCYY + CCYY)
            HESSEE(DY,DZ) = HESSEE(DY,DZ) + (AAYZ + ABYZ + ACYZ
     &                + ABZY + BBYZ + BCYZ + ACZY + BCZY + CCYZ)
            HESSEE(DZ,DZ) = HESSEE(DZ,DZ) + (AAZZ + ABZZ + ACZZ
     &                + ABZZ + BBZZ + BCZZ + ACZZ + BCZZ + CCZZ)
         END IF
#endif
      END IF
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1110) TWOCEN, THRCEN
         WRITE (LUPRI, 1120) ICENT1, ICENT2, ICENT3, ICENT4
         WRITE (LUPRI,'(A,3I5)') ' AX/Y/Z ',AX,AY,AZ
         WRITE (LUPRI,'(A,3I5)') ' BX/Y/Z ',BX,BY,BZ
         WRITE (LUPRI,'(A,3I5)') ' CX/Y/Z ',CX,CY,CZ
         WRITE (LUPRI,'(A,3I5)') ' DX/Y/Z ',DX,DY,DZ
         WRITE (LUPRI, 4000) (DERIV(I),I = 1, NINTYP)
         NCDEP3 = 3*NUCDEP
         IF (IPRINT .GE. 10) THEN
            IF (DERONE) THEN
               WRITE(LUPRI, 4010)
               WRITE(LUPRI, 4020) (GRADEE(I),I=1,NCDEP3)
            END IF
#ifndef PRG_DIRAC
            IF (DERTWO) THEN
               WRITE(LUPRI, 4030)
               DO 800 I = 1,NCDEP3
                  WRITE (LUPRI, 4020)
     *                     (HESSEE(I,J) + HESSEE(J,I), J = 1,I)
  800          CONTINUE
            END IF
#endif
         END IF
      END IF
      RETURN
 1000 FORMAT (//,1X,' <<<<<<<<<< Subroutine DEROUT >>>>>>>>',/)
 1110 FORMAT (2X,'TWOCEN,...:  ',2L5)
 1120 FORMAT (2X,'ICENT1/2/3/4:',4I7)
 4000 FORMAT (2X,'DERIV        ',4E15.6,/(15X,4E15.6))
 4010 FORMAT (//,1X,' Two-electron integral gradient ',/)
 4020 FORMAT (6F12.6)
 4030 FORMAT (//,1X,' Two-electron integral Hessian ',/)
      END

!     radovan: copy-paste-adaptation of DEROUT
!              equivalent to DALTON's DROUTG
      subroutine droutg(deriv, expval, nintyp, iprint)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER AX, AY, AZ, BX, BY, BZ, CX, CY, CZ, DX, DY, DZ
      DIMENSION DERIV(NINTYP), expval(*)
C
#include "nuclei.h"
#include "expcom.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
C
      IF (IPRINT .GE.10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1110) TWOCEN, THRCEN
         WRITE (LUPRI, 1120) ICENT1, ICENT2, ICENT3, ICENT4
      END IF
C
C     ***************************
C     ***** Two-Center Case *****
C     ***************************
C
      IF (TWOCEN) THEN
         AX = 3*ICENT1 - 2
         AY = 3*ICENT1 - 1
         AZ = 3*ICENT1
         BX = 3*ICENT2 - 2
         BY = 3*ICENT2 - 1
         BZ = 3*ICENT2
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
C
C           A electron-repulsion gradient elements:
C
            EXPVAL(AX) = EXPVAL(AX) + DAX
            EXPVAL(AY) = EXPVAL(AY) + DAY
            EXPVAL(AZ) = EXPVAL(AZ) + DAZ
C
C           B electron-repulsion gradient elements:
C
            EXPVAL(BX) = EXPVAL(BX) - DAX
            EXPVAL(BY) = EXPVAL(BY) - DAY
            EXPVAL(BZ) = EXPVAL(BZ) - DAZ
         END IF
C
C
C     *****************************
C     ***** Three-Center Case *****
C     *****************************
C
      ELSE IF (THRCEN) THEN
         AX = 3*ICENT1 - 2
         AY = 3*ICENT1 - 1
         AZ = 3*ICENT1
         BX = 3*ICENT2 - 2
         BY = 3*ICENT2 - 1
         BZ = 3*ICENT2
         CX = 3*ICENT3 - 2
         CY = 3*ICENT3 - 1
         CZ = 3*ICENT3
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
            DBX = DERIV(4)
            DBY = DERIV(5)
            DBZ = DERIV(6)
C
C           A electron-repulsion gradient elements:
C
            EXPVAL(AX) = EXPVAL(AX) + DAX
            EXPVAL(AY) = EXPVAL(AY) + DAY
            EXPVAL(AZ) = EXPVAL(AZ) + DAZ
C
C           B electron-repulsion gradient elements:
C
            EXPVAL(BX) = EXPVAL(BX) + DBX
            EXPVAL(BY) = EXPVAL(BY) + DBY
            EXPVAL(BZ) = EXPVAL(BZ) + DBZ
C
C           C electron-repulsion gradient elements:
C
            EXPVAL(CX) = EXPVAL(CX) - (DAX + DBX)
            EXPVAL(CY) = EXPVAL(CY) - (DAY + DBY)
            EXPVAL(CZ) = EXPVAL(CZ) - (DAZ + DBZ)
         END IF
C
C
C     ****************************
C     ***** Four-Center Case *****
C     ****************************
C
      ELSE
         AX = 3*ICENT1 - 2
         AY = 3*ICENT1 - 1
         AZ = 3*ICENT1
         BX = 3*ICENT2 - 2
         BY = 3*ICENT2 - 1
         BZ = 3*ICENT2
         CX = 3*ICENT3 - 2
         CY = 3*ICENT3 - 1
         CZ = 3*ICENT3
         DX = 3*ICENT4 - 2
         DY = 3*ICENT4 - 1
         DZ = 3*ICENT4
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
            DBX = DERIV(4)
            DBY = DERIV(5)
            DBZ = DERIV(6)
            DCX = DERIV(7)
            DCY = DERIV(8)
            DCZ = DERIV(9)
C
C           A electron-repulsion gradient elements:
C
            EXPVAL(AX) = EXPVAL(AX) + DAX
            EXPVAL(AY) = EXPVAL(AY) + DAY
            EXPVAL(AZ) = EXPVAL(AZ) + DAZ
C
C           B electron-repulsion gradient elements:
C
            EXPVAL(BX) = EXPVAL(BX) + DBX
            EXPVAL(BY) = EXPVAL(BY) + DBY
            EXPVAL(BZ) = EXPVAL(BZ) + DBZ
C
C           C electron-repulsion gradient elements:
C
            EXPVAL(CX) = EXPVAL(CX) + DCX
            EXPVAL(CY) = EXPVAL(CY) + DCY
            EXPVAL(CZ) = EXPVAL(CZ) + DCZ
C
C           D electron-repulsion gradient elements:
C
            EXPVAL(DX) = EXPVAL(DX) - (DAX + DBX + DCX)
            EXPVAL(DY) = EXPVAL(DY) - (DAY + DBY + DCY)
            EXPVAL(DZ) = EXPVAL(DZ) - (DAZ + DBZ + DCZ)
         END IF
C
      END IF
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1110) TWOCEN, THRCEN
         WRITE (LUPRI, 1120) ICENT1, ICENT2, ICENT3, ICENT4
         WRITE (LUPRI,'(A,3I5)') ' AX/Y/Z ',AX,AY,AZ
         WRITE (LUPRI,'(A,3I5)') ' BX/Y/Z ',BX,BY,BZ
         WRITE (LUPRI,'(A,3I5)') ' CX/Y/Z ',CX,CY,CZ
         WRITE (LUPRI,'(A,3I5)') ' DX/Y/Z ',DX,DY,DZ
         WRITE (LUPRI, 4000) (DERIV(I),I = 1, NINTYP)
         NCDEP3 = 3*NUCDEP
         IF (IPRINT .GE. 10) THEN
            IF (DERONE) THEN
               WRITE(LUPRI, 4010)
               WRITE(LUPRI, 4020) (EXPVAL(I),I=1,NCDEP3)
            END IF
         END IF
      END IF
      RETURN
 1000 FORMAT (//,1X,' <<<<<<<<<< Subroutine DEROUT >>>>>>>>',/)
 1110 FORMAT (2X,'TWOCEN,...:  ',2L5)
 1120 FORMAT (2X,'ICENT1/2/3/4:',4I7)
 4000 FORMAT (2X,'DERIV        ',4E15.6,/(15X,4E15.6))
 4010 FORMAT (//,1X,' Two-electron integral gradient ',/)
 4020 FORMAT (6F12.6)
      END
C  /* Deck lndout */
      SUBROUTINE LNDOUT(DERIV,ISYMR,ISYMT,ISYMTS,IPRINT)
C
C     tuh Nov 24 92
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      INTEGER X1, Y1, Z1, X2, Y2, Z2, X1S, X2S
      DIMENSION DERIV(*), G11(3,3), G22(3,3), G12(3,3)
      DIMENSION DIFAB(3), DIFCD(3)
#include "twocom.h"
#include "symmet.h"
#include "suscpt.h"

      XAND(I) = PT(IAND(ISYMAX(1,1),I))
      YAND(I) = PT(IAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IAND(ISYMAX(3,1),I))
      NEXT(I) = MOD(I,3) + 1
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine LNDOUT',-1)
         WRITE (LUPRI,'(A,3I5)') ' ISYMR, ISYMT, ISYMTS ',
     &                             ISYMR, ISYMT, ISYMTS
      END IF
C
      DIFAB(1) =             CORAX0 - XAND(ISYMR )*CORBX0
      DIFAB(2) =             CORAY0 - YAND(ISYMR )*CORBY0
      DIFAB(3) =             CORAZ0 - ZAND(ISYMR )*CORBZ0
      DIFCD(1) = XAND(ISYMT)*CORCX0 - XAND(ISYMTS)*CORDX0
      DIFCD(2) = YAND(ISYMT)*CORCY0 - YAND(ISYMTS)*CORDY0
      DIFCD(3) = ZAND(ISYMT)*CORCZ0 - ZAND(ISYMTS)*CORDZ0
C
      G11(1,1) = DERIV (7)
      G11(1,2) = DERIV (8)
      G11(1,3) = DERIV (9)
      G11(2,1) = DERIV (8)
      G11(2,2) = DERIV(10)
      G11(2,3) = DERIV(11)
      G11(3,1) = DERIV (9)
      G11(3,2) = DERIV(11)
      G11(3,3) = DERIV(12)
C
      G22(1,1) = DERIV(13)
      G22(1,2) = DERIV(14)
      G22(1,3) = DERIV(15)
      G22(2,1) = DERIV(14)
      G22(2,2) = DERIV(16)
      G22(2,3) = DERIV(17)
      G22(3,1) = DERIV(15)
      G22(3,2) = DERIV(17)
      G22(3,3) = DERIV(18)
C
      CALL DCOPY(9,DERIV(19),1,G12,1)
C
      DO 100 X1 = 1, 3
         Y1   = NEXT (X1)
         Z1   = NEXT (Y1)
         ABY1 = DIFAB(Y1)
         ABZ1 = DIFAB(Z1)
         CDY1 = DIFCD(Y1)
         CDZ1 = DIFCD(Z1)
         X1S  = IPTAX(X1,2)
         DO 200 X2 = X1, 3
         IF (ISYMAX(X1,2).EQ.ISYMAX(X2,2)) THEN
            Y2   = NEXT (X2)
            Z2   = NEXT (Y2)
            ABY2 = DIFAB(Y2)
            ABZ2 = DIFAB(Z2)
            CDY2 = DIFCD(Y2)
            CDZ2 = DIFCD(Z2)
            X2S  = IPTAX(X2,2)
            SUS2EL(X1S,X2S) = SUS2EL(X1S,X2S)
     &                      - ABY1*ABY2*G11(Z1,Z2)
     &                      + ABY1*ABZ2*G11(Z1,Y2)
     &                      + ABZ1*ABY2*G11(Y1,Z2)
     &                      - ABZ1*ABZ2*G11(Y1,Y2)
C
     &                      - CDY1*CDY2*G22(Z1,Z2)
     &                      + CDY1*CDZ2*G22(Z1,Y2)
     &                      + CDZ1*CDY2*G22(Y1,Z2)
     &                      - CDZ1*CDZ2*G22(Y1,Y2)
C
     &                      - ABY1*CDY2*G12(Z2,Z1)
     &                      + ABY1*CDZ2*G12(Y2,Z1)
     &                      + ABZ1*CDY2*G12(Z2,Y1)
     &                      - ABZ1*CDZ2*G12(Y2,Y1)
C
     &                      - ABY2*CDY1*G12(Z1,Z2)
     &                      + ABY2*CDZ1*G12(Y1,Z2)
     &                      + ABZ2*CDY1*G12(Z1,Y2)
     &                      - ABZ2*CDZ1*G12(Y1,Y2)
         END IF
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 4000) (DERIV(I),I = 1, 27)
         WRITE(LUPRI,'(/A)') ' Unfinished two-electron susceptibilities'
         CALL OUTPUT(SUS2EL,1,3,1,3,3,3,1,LUPRI)
         IF (IPRINT .GT. 10) THEN
            WRITE (LUPRI,'(/A)') ' G11 '
            CALL OUTPUT(G11,1,3,1,3,3,3,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' G22 '
            CALL OUTPUT(G22,1,3,1,3,3,3,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' G12 '
            CALL OUTPUT(G12,1,3,1,3,3,3,1,LUPRI)
         END IF
      END IF
      RETURN
 1000 FORMAT (//,1X,' <<<<<<<<<< Subroutine LNDOUT >>>>>>>>',/)
 4000 FORMAT (2X,'DERIV        ',4E15.6,/(15X,4E15.6))
      END
C  /* Deck lndgc1 */
      SUBROUTINE LNDGC1(DERIV,EXPVAL,ISYMR,ISYMT,ISYMTS,IPRINT)
C
C     K.Ruud, May 2007, based on LNDOUT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      INTEGER X1, Y1, Z1, X2, Y2, Z2, X1S, X2S
      DIMENSION DERIV(*), G1(3), G2(3)
      DIMENSION DIFAB(3), DIFCD(3), EXPVAL(3)
#include "twocom.h"
#include "symmet.h"
      XAND(I) = PT(IAND(ISYMAX(1,1),I))
      YAND(I) = PT(IAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IAND(ISYMAX(3,1),I))
      NEXT(I) = MOD(I,3) + 1
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine LNDGC1',-1)
         WRITE (LUPRI,'(A,3I5)') ' ISYMR, ISYMT, ISYMTS ',
     &                             ISYMR, ISYMT, ISYMTS
      END IF
C
      DIFAB(1) =             CORAX0 - XAND(ISYMR )*CORBX0
      DIFAB(2) =             CORAY0 - YAND(ISYMR )*CORBY0
      DIFAB(3) =             CORAZ0 - ZAND(ISYMR )*CORBZ0
      DIFCD(1) = XAND(ISYMT)*CORCX0 - XAND(ISYMTS)*CORDX0
      DIFCD(2) = YAND(ISYMT)*CORCY0 - YAND(ISYMTS)*CORDY0
      DIFCD(3) = ZAND(ISYMT)*CORCZ0 - ZAND(ISYMTS)*CORDZ0
C
      G1(1) = DERIV(1)
      G1(2) = DERIV(2)
      G1(3) = DERIV(3)
C
      G2(1) = DERIV(4)
      G2(2) = DERIV(5)
      G2(3) = DERIV(6)
C
      DO 100 X1 = 1, 3
         Y1   = NEXT (X1)
         Z1   = NEXT (Y1)
         ABY1 = DIFAB(Y1)
         ABZ1 = DIFAB(Z1)
         CDY1 = DIFCD(Y1)
         CDZ1 = DIFCD(Z1)
         X1S  = IPTAX(X1,2)
         EXPVAL(X1S) = EXPVAL(X1S)
     &               - ABY1*G1(Z1)
     &               + ABZ1*G1(Y1)
     &               - CDY1*G2(Z1)
     &               + CDZ1*G2(Y1)
  100 CONTINUE
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 4000) (DERIV(I),I = 1, 27)
         WRITE(LUPRI,'(/A)') ' Unfinished contration with B1 integrals'
         CALL OUTPUT(EXPVAL,1,3,1,1,3,1,1,LUPRI)
         IF (IPRINT .GT. 10) THEN
            WRITE (LUPRI,'(/A)') ' G1 '
            CALL OUTPUT(G1,1,3,1,1,3,1,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' G2 '
            CALL OUTPUT(G2,1,3,1,1,3,1,1,LUPRI)
         END IF
      END IF
      RETURN
 1000 FORMAT (//,1X,' <<<<<<<<<< Subroutine LNDGC1 >>>>>>>>',/)
 4000 FORMAT (2X,'DERIV        ',4E15.6,/(15X,4E15.6))
      END
C  /* Deck lndgc2 */
      SUBROUTINE LNDGC2(DERIV,EXPVAL,ISYMR,ISYMT,ISYMTS,IPRINT)
C
C     K.Ruud, May 2007, based on LNDOUT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      INTEGER X1, Y1, Z1, X2, Y2, Z2, X1S, X2S
      DIMENSION DERIV(*), G11(3,3), G22(3,3), G12(3,3)
      DIMENSION DIFAB(3), DIFCD(3), EXPVAL(3,3)
#include "twocom.h"
#include "symmet.h"
      XAND(I) = PT(IAND(ISYMAX(1,1),I))
      YAND(I) = PT(IAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IAND(ISYMAX(3,1),I))
      NEXT(I) = MOD(I,3) + 1
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine LNDGC2',-1)
         WRITE (LUPRI,'(A,3I5)') ' ISYMR, ISYMT, ISYMTS ',
     &                             ISYMR, ISYMT, ISYMTS
      END IF
C
      DIFAB(1) =             CORAX0 - XAND(ISYMR )*CORBX0
      DIFAB(2) =             CORAY0 - YAND(ISYMR )*CORBY0
      DIFAB(3) =             CORAZ0 - ZAND(ISYMR )*CORBZ0
      DIFCD(1) = XAND(ISYMT)*CORCX0 - XAND(ISYMTS)*CORDX0
      DIFCD(2) = YAND(ISYMT)*CORCY0 - YAND(ISYMTS)*CORDY0
      DIFCD(3) = ZAND(ISYMT)*CORCZ0 - ZAND(ISYMTS)*CORDZ0
C
      G11(1,1) = DERIV (7)
      G11(1,2) = DERIV (8)
      G11(1,3) = DERIV (9)
      G11(2,1) = DERIV (8)
      G11(2,2) = DERIV(10)
      G11(2,3) = DERIV(11)
      G11(3,1) = DERIV (9)
      G11(3,2) = DERIV(11)
      G11(3,3) = DERIV(12)
C
      G22(1,1) = DERIV(13)
      G22(1,2) = DERIV(14)
      G22(1,3) = DERIV(15)
      G22(2,1) = DERIV(14)
      G22(2,2) = DERIV(16)
      G22(2,3) = DERIV(17)
      G22(3,1) = DERIV(15)
      G22(3,2) = DERIV(17)
      G22(3,3) = DERIV(18)
C
      CALL DCOPY(9,DERIV(19),1,G12,1)
C
      DO 100 X1 = 1, 3
         Y1   = NEXT (X1)
         Z1   = NEXT (Y1)
         ABY1 = DIFAB(Y1)
         ABZ1 = DIFAB(Z1)
         CDY1 = DIFCD(Y1)
         CDZ1 = DIFCD(Z1)
         X1S  = IPTAX(X1,2)
         DO 200 X2 = X1, 3
         IF (ISYMAX(X1,2).EQ.ISYMAX(X2,2)) THEN
            Y2   = NEXT (X2)
            Z2   = NEXT (Y2)
            ABY2 = DIFAB(Y2)
            ABZ2 = DIFAB(Z2)
            CDY2 = DIFCD(Y2)
            CDZ2 = DIFCD(Z2)
            X2S  = IPTAX(X2,2)
            EXPVAL(X1S,X2S) = EXPVAL(X1S,X2S)
     &                      - ABY1*ABY2*G11(Z1,Z2)
     &                      + ABY1*ABZ2*G11(Z1,Y2)
     &                      + ABZ1*ABY2*G11(Y1,Z2)
     &                      - ABZ1*ABZ2*G11(Y1,Y2)
C
     &                      - CDY1*CDY2*G22(Z1,Z2)
     &                      + CDY1*CDZ2*G22(Z1,Y2)
     &                      + CDZ1*CDY2*G22(Y1,Z2)
     &                      - CDZ1*CDZ2*G22(Y1,Y2)
C
     &                      - ABY1*CDY2*G12(Z2,Z1)
     &                      + ABY1*CDZ2*G12(Y2,Z1)
     &                      + ABZ1*CDY2*G12(Z2,Y1)
     &                      - ABZ1*CDZ2*G12(Y2,Y1)
C
     &                      - ABY2*CDY1*G12(Z1,Z2)
     &                      + ABY2*CDZ1*G12(Y1,Z2)
     &                      + ABZ2*CDY1*G12(Z1,Y2)
     &                      - ABZ2*CDZ1*G12(Y1,Y2)
         END IF
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 4000) (DERIV(I),I = 1, 27)
         WRITE(LUPRI,'(/A)') ' Unfinished two-electron susceptibilities'
         CALL OUTPUT(EXPVAL,1,3,1,3,3,3,1,LUPRI)
         IF (IPRINT .GT. 10) THEN
            WRITE (LUPRI,'(/A)') ' G11 '
            CALL OUTPUT(G11,1,3,1,3,3,3,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' G22 '
            CALL OUTPUT(G22,1,3,1,3,3,3,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' G12 '
            CALL OUTPUT(G12,1,3,1,3,3,3,1,LUPRI)
         END IF
      END IF
      RETURN
 1000 FORMAT (//,1X,' <<<<<<<<<< Subroutine LNDGC2 >>>>>>>>',/)
 4000 FORMAT (2X,'DERIV        ',4E15.6,/(15X,4E15.6))
      END
C  /* Deck pblock */
      SUBROUTINE PBLOCK(PSO,PAO,ICOMPA,ICOMPB,ICOMPC,ICOMPD,
     &                  NHKTA,NHKTB,NHKTC,NHKTD,
     &                  KHKTA,KHKTB,KHKTC,MULA,MULB,MULC,MULD,
     &                  NORBA,NORBB,NORBC,NORBD,ISYMR,ISYMS,ISYMT)
#include "implicit.h"
C
C     *********************************************************
C     ***  Transform a block of the P-matrix from SO basis  ***
C     ***  to AO basis.  This is done for fixed component   ***
C     ***  indices ICOMPx and fixed symmetry operations.    ***
C     *********************************************************
C
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "symmet.h"
      DIMENSION PSO(*), PAO(*)

      IAR(I,J,K,L) = KHKTA*(KHKTB*(KHKTC*(L-1)+K-1)+J-1)+I
      IRR(I,J,K,L) = MULTA*(MULTB*(MULTC*L+K)+J)+I+1
      ISYMTS = IEOR(ISYMT,ISYMS)
      MULTA = MULT(MULA)
      MULTB = MULT(MULB)
      MULTC = MULT(MULC)
      MULTD = MULT(MULD)
      NCABCD = NORBA*NORBB*NORBC*NORBD
      MULTOT = IRR(MULTA-1,MULTB-1,MULTC-1,MULTD-1)
      ITYNA = ISYMAO(NHKTA,ICOMPA)
      ITYNB = ISYMAO(NHKTB,ICOMPB)
      ITYNC = ISYMAO(NHKTC,ICOMPC)
      ITYND = ISYMAO(NHKTD,ICOMPD)
      CALL DZERO(PAO,NCABCD)
C
C     Loop over irreps
C
      IRCNTA = -1
      DO 100 IREPA = 0, MAXREP
      IF(IAND(MULA,IEOR(IREPA,ITYNA)) .EQ. 0) THEN
         IRCNTA = IRCNTA + 1
         IRCNTB = -1
         DO 200 IREPB = 0, MAXREP
         IF(IAND(MULB,IEOR(IREPB,ITYNB)) .EQ. 0) THEN
            SIGNB = PT(IAND(ISYMR, IEOR(IREPB,ITYNB)))
            IRCNTB = IRCNTB + 1
            IRCNTC = -1
            DO 300 IREPC = 0, MAXREP
            IF(IAND(MULC,IEOR(IREPC,ITYNC)) .EQ. 0) THEN
               SIGNBC = SIGNB*PT(IAND(ISYMT, IEOR(IREPC,ITYNC)))
               IRPABC = IEOR(IREPA,IEOR(IREPB,IREPC))
               IRCNTC = IRCNTC + 1
               IRCNTD = -1
               DO 400 IREPD = 0,MAXREP
               IF (IAND(MULD,IEOR(IREPD,ITYND)) .EQ. 0) THEN
                  IRCNTD = IRCNTD + 1
                  IF (IEOR(IREPD,IRPABC) .EQ. 0) THEN
                     SIGN =SIGNBC*PT(IAND(ISYMTS,IEOR(IREPD,ITYND)))
                     IOFF = NCABCD*(MULTOT
     &                      *(IAR(ICOMPA,ICOMPB,ICOMPC,ICOMPD) - 1)
     &                      + IRR(IRCNTA,IRCNTB,IRCNTC,IRCNTD) - 1)
                     DO 500 I = 1,NCABCD
                         PAO(I) = PAO(I) + SIGN*PSO(IOFF+I)
  500                CONTINUE
                  END IF
               END IF
  400          CONTINUE
            END IF
  300       CONTINUE
         END IF
  200    CONTINUE
      END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck pprim */
      SUBROUTINE PPRIM(PAO,PUNC,WORK,LWORK,NORBA,NORBB,NORBC,NORBD,
     &                 NUCA,NUCB,NUCC,NUCD,TPRIAB,TPRICD)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      LOGICAL TPRIAB,TPRICD
      DIMENSION PAO(*), PUNC(*), WORK(LWORK)
#include "twosta.h"
C
C     Allocations
C
      IF (TPRIAB) THEN
         LPHALF = (NUCA*(NUCA + 1)/2)*NORBC*NORBD
      ELSE
         LPHALF = NUCA*NUCB*NORBC*NORBD
      END IF
      KPHALF = 1
      KQTR   = KPHALF + LPHALF
      KHALF  = KQTR   + MAX(NUCA*NORBB,NUCC*NORBD)
      KLAST  = KHALF  + MAX(NUCA*NUCB,NUCC*NUCD)
      IF (KLAST .GT. LWORK) CALL STOPIT('PPRIM',' ',KLAST,LWORK)
      MWPPRI = MAX(MWPPRI,KLAST)
      LWTOT  = LWTOT + KLAST
      MWTOT  = MAX(MWTOT,LWTOT)
      CALL PPRIM1(PAO,PUNC,WORK(KHALF),WORK(KQTR),WORK(KPHALF),NORBA,
     &            NORBB,NORBC,NORBD,NUCA,NUCB,NUCC,NUCD,TPRIAB,TPRICD)
      LWTOT  = LWTOT - KLAST
      RETURN
      END
C  /* Deck pprim1 */
      SUBROUTINE PPRIM1(PAO,PUNC,QTR,HALF,PHALF,NORBA,NORBB,NORBC,NORBD,
     &                  NUCA,NUCB,NUCC,NUCD,TPRIAB,TPRICD)
#include "implicit.h"
C
C     *********************************************************
C     ***  Transform a block of the density matrix elements ***
C     ***  from contracted basis functions to primitive     ***
C     ***  basis functions.  The array PAO contains the     ***
C     ***  desired combination of one- and two-electron     ***
C     ***  density matrix elements.                         ***
C     *********************************************************
C
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      LOGICAL TPRIAB,TPRICD
#include "symmet.h"
      DIMENSION PAO(*), QTR(*), HALF(*), PUNC(*), PHALF(NORBC*NORBD,*)
C
C     Loop over second pair of indices and transform over the
C     first pair of indices
C
      NUCABQ = NUCA*NUCB
      IF (TPRIAB) NUCABQ = (NUCA*(NUCA+1))/2
      NUCCDQ = NUCC*NUCD
      IF (TPRICD) NUCCDQ = (NUCC*(NUCC+1))/2
      KLX = 0
      DO 10 KL = 1,NORBC*NORBD
         CALL MXM(CONTT1,NUCA,PAO(KLX+1),NORBA,QTR,NORBB)
         CALL MXM(QTR,NUCA,CONT2,NORBB,HALF,NUCB)
C
C     Scatter these values so that the outer index is IJ
C
         IF (.NOT. TPRIAB) THEN
            IJ = 0
            DO 20 I = 1,NUCA
               DO 23 J = 1,NUCB
                  IJ = IJ + 1
                  JI = NUCA*(J-1) + I
                  PHALF(KL,IJ) = HALF(JI)
23             CONTINUE
20          CONTINUE
         ELSE
            IJ = 0
            DO 25 I = 1,NUCA
               DO 27 J = 1,I
                  IJ = IJ + 1
                  INJ = NUCA*(J-1) + I
                  PHALF(KL,IJ) = HALF(INJ)
                  IF (I .NE. J) THEN
                     INJ = NUCA*(I-1) + J
                     PHALF(KL,IJ) = PHALF(KL,IJ) + HALF(INJ)
                  END IF
27             CONTINUE
25          CONTINUE
         END IF
         KLX = KLX + NORBA*NORBB
10    CONTINUE
C
C     Now loop over IJ pairs, transforming each block of KL
C
      IJ = 0
      DO 40 I = 1,NUCA
         JJ = NUCB
         IF (TPRIAB) JJ = I
         DO 50 J = 1,JJ
            IJ = IJ + 1
            CALL MXM(CONTT3,NUCC,PHALF(1,IJ),NORBC,QTR,NORBD)
            CALL MXM(QTR,NUCC,CONT4,NORBD,HALF,NUCD)
C
C     Scatter the final values into PUNC
C
            KL = 0
            DO 60 K = 1,NUCC
               IF (TPRICD) THEN
                  LL = K
               ELSE
                  LL = NUCD
               END IF
               DO 70 L = 1,LL
                  KL = KL + 1
                  KNL = NUCC*(L-1) + K
                  IPADR = NUCCDQ*(IJ-1)+KL
                  PUNC(IPADR) = HALF(KNL)
                  IF (TPRICD.AND.K.NE.L) THEN
                     KNL = NUCC*(K-1) + L
                     PUNC(IPADR) = PUNC(IPADR) + HALF(KNL)
                  END IF
70             CONTINUE
60          CONTINUE
50       CONTINUE
40    CONTINUE
      RETURN
      END

C  /* Deck intfcl */
      SUBROUTINE INTFCL(AOINT,DMAT,NDMAT,PSO,PSA,FT,FV,DERIV,PAO,PAA,
     &     PINT,QINT,WORK,LWORK,NCCINT,NINTYP,
     &     ISYMR,ISYMS,ISYMT,ICORBA,ICORBB,ICORBC,ICORBD,
     &     THRESH,SYMFAC,IPRINT,NOPMAT,NODV,NUCABQ,NUCCDQ,
     &     NINDAB,NINDCD,EXPECT,LONDON,SUSCEP,
     &     DDFOCK,DINTSKP)
C     
C     Calculates expectation values and Fock matrices
C     of differentiated integrals
C
C  MI: It's TEC's routine for implementation of London orbitals.
C   Odense, 2002,2003
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (DP25 = 0.25 D00, DP5 = 0.5 D00, 
     &     D1 = 1.0 D00, D2 = 2.0 D00, ZERADD = 1.D-15,
     &     D0 = 0.0 D00, DP125 = 0.125D00,
     &     D4 = 4.0 D00, D8 = 8.0 D00)
      INTEGER A, B, C, D, X, Y, Z, Y2, Z2
      LOGICAL NODER, DCMPAB, DCMPCD, NOPMAT, NODV, DV, SUSCEP, DDFOCK,
     &     EXPECT, LONDON
C     jth - bugfix ICORBD(NORBC)  -->  ICORBD(NORBD)
C     DIMENSION ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBC),
      DIMENSION ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBD),
     &     DERIV(NINTYP), PSO(*), PSA(*),
     &     PAO(NORBA,NORBB,NORBC,NORBD),
     &     PAA(NORBA,NORBB,NORBC,NORBD),
     &     AOINT(NCCINT,NINTYP),
     &     NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &     DMAT(NBASIS,NBASIS,NDMAT), WORK(LWORK),
     &     FT(NBASIS,NBASIS,NDMAT,*), FV(NBASIS,NBASIS,NDMAT,*),
     &     DINT(12), NCART(12), PINT(NCCINT,*), QINT(NCCINT,*),
     &     DIFAB(3), DIFCD(3), DINTSKP(*)
#include "twocom.h"
#include "symmet.h"
#include "nuclei.h"
#include "expcom.h"
#include "dirprt.h"
#include "dcbgrd.h"
#include "doxyz.h"
#include "dcbham.h"
      logical :: allow=.false.
C
C     jth
C     
      XAND(I) = PT(IAND(ISYMAX(1,1),I))
      YAND(I) = PT(IAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IAND(ISYMAX(3,1),I))
      NEXT(I) = MOD(I,3) + 1


      IF (IPRINT.GE.4) THEN
       WRITE(LUPRI,'(/8X,A)') 'Output from INTFCL'
       WRITE(LUPRI,'(10X,A,L2)') 'LONDON:',LONDON
       WRITE(LUPRI,'(10X,A,L2)') 'SUSCEP:',SUSCEP
       WRITE(LUPRI,'(10X,A,L2)') 'DDFOCK:',DDFOCK
       WRITE(LUPRI,'(10X,A,L2)') 'EXPECT:',EXPECT
       WRITE(LUPRI,'(10X,A,L2)') 'NOPMAT:',NOPMAT
       WRITE(LUPRI,'(10X,A,L2)') 'NODV:',NODV
      ENDIF
C     
      ISYMTS   = IEOR(ISYMT,ISYMS)
      DIFAB(1) =             CORAX0 - XAND(ISYMR )*CORBX0
      DIFAB(2) =             CORAY0 - YAND(ISYMR )*CORBY0
      DIFAB(3) =             CORAZ0 - ZAND(ISYMR )*CORBZ0
      DIFCD(1) = XAND(ISYMT)*CORCX0 - XAND(ISYMTS)*CORDX0
      DIFCD(2) = YAND(ISYMT)*CORCY0 - YAND(ISYMTS)*CORDY0
      DIFCD(3) = ZAND(ISYMT)*CORCZ0 - ZAND(ISYMTS)*CORDZ0


      IF (.NOT.SUSCEP) THEN
         DO 10 N = 1, 3
            Y   = NEXT(N)
            Z   = NEXT(Y)
            ABY = DIFAB(Y)
            ABZ = DIFAB(Z)
            CDY = DIFCD(Y)
            CDZ = DIFCD(Z)
            Y2  = Y + 3
            Z2  = Z + 3
            DO 20 I = 1, NCCINT
               AOAB = ABY*AOINT(I,Z ) - ABZ*AOINT(I,Y )
               AOCD = CDY*AOINT(I,Z2) - CDZ*AOINT(I,Y2)
               PINT(I,N) = AOAB + AOCD
               QINT(I,N) = AOAB - AOCD
  20        CONTINUE
  10     CONTINUE
         ELSE
         DO J = 1, 6
            K = (7+2*J)/5; M = J-K*(K-1)/2
            K = NEXT(K); L = NEXT(K)
            M = NEXT(M); N = NEXT(M)
            KK = 1+(8-MIN(K,M))*(MIN(K,M)-1)/2+ABS(K-M)
            KL = 1+(8-MIN(K,N))*(MIN(K,N)-1)/2+ABS(K-N)
            LK = 1+(8-MIN(L,M))*(MIN(L,M)-1)/2+ABS(L-M)
            LL = 1+(8-MIN(L,N))*(MIN(L,N)-1)/2+ABS(L-N)
            KM = 3*(K-1) + M; MK = 3*(M-1) + K
            KN = 3*(K-1) + N; NK = 3*(N-1) + K
            LM = 3*(L-1) + M; ML = 3*(M-1) + L
            LN = 3*(L-1) + N; NL = 3*(N-1) + L
            ABK = DIFAB(K); ABL = DIFAB(L)
            ABM = DIFAB(M); ABN = DIFAB(N)
            CDK = DIFCD(K); CDL = DIFCD(L)
            CDM = DIFCD(M); CDN = DIFCD(N)
            DO I = 1, NCCINT
               PPQQ = ABL*ABN*AOINT(I, 6+KK) - ABL*ABM*AOINT(I, 6+KL)
     &              - ABK*ABN*AOINT(I, 6+LK) + ABK*ABM*AOINT(I, 6+LL)
     &              + CDL*CDN*AOINT(I,12+KK) - CDL*CDM*AOINT(I,12+KL)
     &              - CDK*CDN*AOINT(I,12+LK) + CDK*CDM*AOINT(I,12+LL)
               PQPQ = ABL*CDN*AOINT(I,18+KM) - ABL*CDM*AOINT(I,18+KN)
     &              - ABK*CDN*AOINT(I,18+LM) + ABK*CDM*AOINT(I,18+LN)
     &              + ABN*CDL*AOINT(I,18+MK) - ABN*CDK*AOINT(I,18+ML)
     &              - ABM*CDL*AOINT(I,18+NK) + ABM*CDK*AOINT(I,18+NL)
               PINT(I,J) = -(PPQQ + PQPQ)/4
               QINT(I,J) = -(PPQQ - PQPQ)/4
            END DO
         END DO
      END IF

      SFAC = SYMFAC
      IF (.NOT.SHABAB) SFAC = D2*SFAC
      IF (.NOT.SHAEQB) SFAC = D2*SFAC
      IF (.NOT.SHCEQD) SFAC = D2*SFAC
C     
      THRSH = MAX(THRESH,ZERADD)
C     
c     IF(DOSCRN .AND. .NOT.DDFOCK) THEN
C     ... DNSMAX fits EXPECT, not DDFOCK
c     EXPTOL = MAX((FCKTHR/DNSMAX),1.00D-15)
c     ELSE
      EXPTOL = 1.00D-15
c     ENDIF
C     
C     
C     ***** Loop over shell components *****
C     
      IAOFF = 1
      DO 100 ICOMPA = 1, KHKTA
         KHKTBB = KHKTB
         IF (DIAGAB) KHKTBB = ICOMPA
         DO 200 ICOMPB = 1, KHKTBB
            DCMPAB = DIAGAB .AND. ICOMPA .EQ. ICOMPB
            FACAB  = D1
            IF (DIAGAB .AND. ICOMPA .NE. ICOMPB) FACAB = D2
            DO 300 ICOMPC = 1, KHKTC
               KHKTDD = KHKTD
               IF (DIAGCD) KHKTDD = ICOMPC
               DO 400 ICOMPD = 1, KHKTDD
C     
C     Step 2 screening on gradient, hessian, London, and
C     susceptibilities.
C     Always do Step 2 screening 
C     
                  AOMAX = D0
                  DO I = 1,NINTYP
                     IAOMAX = IDAMAX(NOABCD,AOINT(IAOFF,I),1)-1
                     AOMAX = MAX(AOMAX,ABS(AOINT(IAOFF+IAOMAX,I)))
                  END DO
                  IF (AOMAX .LT. EXPTOL) GO TO 510

C     
C     
                  DCMPCD = DIAGCD .AND. ICOMPC .EQ. ICOMPD
                  FACCD  = D1
                  IF (DIAGCD .AND. ICOMPC .NE. ICOMPD) FACCD = D2
                  FCABCD = FACAB*FACCD*SFAC
C     
                  INT = IAOFF
C     
C     
C     ***** Loop over contracted functions *****
C     
                  DO 500 IORBAB = 1, NORBAB
                     IORBA = NINDAB(IORBAB,1)
                     IORBB = NINDAB(IORBAB,2)
                     A = ICORBA(IORBA) + ICOMPA
                     B = ICORBB(IORBB) + ICOMPB
                     FAB = D1
                     IF (TCONAB .AND. IORBA.NE.IORBB) FAB = D2
                     DO 600 IORBCD = 1, NORBCD
                        IORBC = NINDCD(IORBCD,1)
                        IORBD = NINDCD(IORBCD,2)
                        C = ICORBC(IORBC) + ICOMPC
                        D = ICORBD(IORBD) + ICOMPD
                        FCD = D1
                        IF (TCONCD.AND.IORBC.NE.IORBD) FCD = D2
C     
C     *******************************************
C     ***** Two-electron Fock London matrix *****
C     *******************************************
C     
CMI/Dec 6 2004 Added comments to be clear how 
C             relativistic terms were derived 
C             based on DALTON integrals formalism
C
CMI TODO: Carefully analyze these formulas ! 
C

                        dfac = dp125*fab*fcd*fcabcd
                        efac = hfxfac*dp5*dfac

!radovan                cancelling a factor 2 somewhere in the code
                        dfac = 0.5d0*dfac
                        efac = 0.5d0*efac

!                                coulomb
!                                =======
            
                        dab = dfac*dmat(a, b, 1)
                        dcd = dfac*dmat(c, d, 1)
                                               
                        dba = dfac*dmat(b, a, 1)
                        ddc = dfac*dmat(d, c, 1)

                        if (suscep) then
                           do i = 1, 6

                              abcd = pint(int, i)
                              abdc = qint(int, i)
            
             ft(a, b, 1, i) = ft(a, b, 1, i) + abcd*ddc + abdc*dcd
             ft(c, d, 1, i) = ft(c, d, 1, i) + abcd*dba + abdc*dab

             ft(b, a, 1, i) = ft(b, a, 1, i) + abcd*dcd + abdc*ddc
             ft(d, c, 1, i) = ft(d, c, 1, i) + abcd*dab + abdc*dba
                           end do
                        else
                           do i = 1, 3
            
                              abcd = pint(int, i)
                              abdc = qint(int, i)
            
             ft(b, a, 1, i) = ft(b, a, 1, i) - abcd*dcd - abdc*ddc
             ft(d, c, 1, i) = ft(d, c, 1, i) - abcd*dab + abdc*dba
                                                              
             ft(a, b, 1, i) = ft(a, b, 1, i) + abcd*ddc + abdc*dcd
             ft(c, d, 1, i) = ft(c, d, 1, i) + abcd*dba - abdc*dab
            
                           end do
                        end if
            
!                                exchange
!                                ========

                        if (dabs(efac) > tiny(0.0d0)) then
                        !if (dabs(efac) > tiny(0.0d0).and.allow) then
                          do imat = 1, ndmat

                           dac = efac*dmat(a, c, imat)
                           dad = efac*dmat(a, d, imat)
                           dbc = efac*dmat(b, c, imat)
                           dbd = efac*dmat(b, d, imat)
                                                  
                           dca = efac*dmat(c, a, imat)
                           dda = efac*dmat(d, a, imat)
                           dcb = efac*dmat(c, b, imat)
                           ddb = efac*dmat(d, b, imat)

                           if (suscep) then
                              do i = 1, 6

                                 abcd = pint(int, i)
                                 abdc = qint(int, i)
                                
                        ft(a, c, imat, i) = ft(a, c, imat, i) - abdc*dbd
                        ft(a, d, imat, i) = ft(a, d, imat, i) - abcd*dbc
                        ft(b, c, imat, i) = ft(b, c, imat, i) - abcd*dad
                        ft(b, d, imat, i) = ft(b, d, imat, i) - abdc*dac
                                                                        
                        ft(c, a, imat, i) = ft(c, a, imat, i) - abdc*ddb
                        ft(d, a, imat, i) = ft(d, a, imat, i) - abcd*dcb
                        ft(c, b, imat, i) = ft(c, b, imat, i) - abcd*dda
                        ft(d, b, imat, i) = ft(d, b, imat, i) - abdc*dca
            
                              end do
                           else
                              do i = 1, 3
            
                                 abcd = pint(int, i)
                                 abdc = qint(int, i)
            
                        ft(c, a, imat, i) = ft(c, a, imat, i) + abdc*ddb
                        ft(d, a, imat, i) = ft(d, a, imat, i) + abcd*dcb
                        ft(c, b, imat, i) = ft(c, b, imat, i) - abcd*dda
                        ft(d, b, imat, i) = ft(d, b, imat, i) - abdc*dca
                    
                        ft(a, c, imat, i) = ft(a, c, imat, i) - abdc*dbd
                        ft(a, d, imat, i) = ft(a, d, imat, i) - abcd*dbc
                        ft(b, c, imat, i) = ft(b, c, imat, i) + abcd*dad
                        ft(b, d, imat, i) = ft(b, d, imat, i) + abdc*dac
            
                              end do
                           end if

                        end do
                        end if



CMI  ... control print-out
                        IF (IPRINT .GT. 30) THEN
                          WRITE(LUPRI,'(//A,4I5)')
     &                             ' IORBA...IORBD ',
     &                           IORBA,IORBB,IORBC,IORBD
                          WRITE(LUPRI,'(A,4I5,A,I5)')
     &                     'indexes A,B,C,D',A,B,C,D,
     &                     ' INT=',INT   
                          WRITE(LUPRI,'(A)')
     &                     'actual integrals and Fock matrix elements:'
                          DO IXYZ=1,3
                              ABCD = PINT(INT,IXYZ)
                              ABDC = QINT(INT,IXYZ)
                           WRITE(LUPRI,'(A,2I3,A,2I3,A,I1,A,D12.6,$)')
     &                      'int (',A,B,'|',C,D,')_',IXYZ,'=',ABCD
                           WRITE(LUPRI,'(A,2I3,A,2I3,A,I1,A,D12.6)')
     &                      '  int (',A,B,'|',D,C,')_',IXYZ,'=',ABDC
                           DO IMAT=1,NDMAT
                        WRITE(LUPRI,'(A,I3,A,I3,A,I1,A,I1,A,D12.6,$)')
     &                     ' FT(',B,',',A,',',IMAT,',',IXYZ,')=',
     &                      FT(B,A,IMAT,IXYZ)
                        WRITE(LUPRI,'(A,I3,A,I3,A,I1,A,I1,A,D12.6,$)')
     &                     '  FT(',D,',',C,',',IMAT,',',IXYZ,')=',
     &                      FT(D,C,IMAT,IXYZ)
                        WRITE(LUPRI,'(A,I3,A,I3,A,I1,A,I1,A,D12.6)')
     &                     '  FT(',C,',',A,',',IMAT,',',IXYZ,')=',
     &                      FT(C,A,IMAT,IXYZ)
                        WRITE(LUPRI,'(A,I3,A,I3,A,I1,A,I1,A,D12.6,$)')
     &                     ' FT(',D,',',A,',',IMAT,',',IXYZ,')=',
     &                      FT(D,A,IMAT,IXYZ)
                        WRITE(LUPRI,'(A,I3,A,I3,A,I1,A,I1,A,D12.6,$)')
     &                     '  FT(',C,',',B,',',IMAT,',',IXYZ,')=',
     &                      FT(C,B,IMAT,IXYZ)
                        WRITE(LUPRI,'(A,I3,A,I3,A,I1,A,I1,A,D12.6)')
     &                     '  FT(',D,',',B,',',IMAT,',',IXYZ,')=',
     &                      FT(D,B,IMAT,IXYZ)
                           ENDDO
                          ENDDO
                        END IF
                        INT = INT + 1
 600                 CONTINUE
 500              CONTINUE
 510              IAOFF = IAOFF + NOABCD
 400           CONTINUE
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE

      RETURN
      END


C  /* Deck polpri */
      SUBROUTINE POLPRI(AMAT,SPC,ITYPE1)
CMI ... stolen from dalton(abacus/abadip.F)
C      for making the dirac output of magnetizabilities 
C      the same as is in the older brother dalton

#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C#include <debye.h>
      PARAMETER (D1 = 1.0 D00)
      LOGICAL ALLREP, INERT, LARGE
      DIMENSION AMAT(3,3), IND(3)
      CHARACTER*(*) SPC
      CHARACTER LAB(3)*2, CHRABC(3)*1
#include "inirep.h"
#include "abainf.h"
#include "inftap.h"
#include "dorps.h"
#include "symmet.h"
#include "chrxyz.h"
C
      DATA CHRABC /'A','B','C'/
C
      INERT = INDEX(SPC,'PRI') .NE. 0
      LARGE = ITYPE1 .EQ. -2
      ITYPE = ABS(ITYPE1)
C
C     *****************************************
C     ***** Units (atomic or Angstrom**3) *****
C     *****************************************
C
      IF (INDEX(SPC,'EXP') .NE. 0) THEN
         FAC = XTANG*XTANG*XTANG
      ELSE IF (INDEX(SPC,'CGS') .NE. 0) THEN
         FAC = XTANG*DEBYE
      ELSE
         FAC = D1
      END IF
C
C     ***************************************
C     ***** Pick up components to print *****
C     ***************************************
C
      IF (.NOT.INERT) THEN
         NCOMP = 0
         DO 100 ICOOR = 1, 3
            IF (DOSYM(ISYMAX(ICOOR,ITYPE) + 1)
     &          .OR. ITYPE .EQ. 2) THEN
               NCOMP = NCOMP + 1
               IND(NCOMP) = IPTAX(ICOOR,ITYPE)
               IF (ITYPE .EQ. 1) THEN
                  LAB(NCOMP) = 'E'//CHRXYZ(-ICOOR)
               ELSE
                  LAB(NCOMP) = 'B'//CHRXYZ(-ICOOR)
               END IF
            END IF
  100    CONTINUE
      ELSE
         NCOMP = 0
         DO 200 ICOOR = 1, 3
            ALLREP = .TRUE.
            DO 210 IREPS = 1, NREPPI(ICOOR)
               IF (.NOT.DOREPS(IREPPI(ICOOR,IREPS))) ALLREP = .FALSE.
  210       CONTINUE
            IF (ALLREP) THEN
               NCOMP = NCOMP + 1
               IND(NCOMP) = ICOOR
               LAB(NCOMP) = 'E'//CHRABC(ICOOR)
            END IF
  200    CONTINUE
      END IF
C
C     *****************
C     ***** Print *****
C     *****************
C
      IF (NCOMP .GT. 0) THEN
         IF (INERT) WRITE (LUPRI,'(15X,A/)')
     &      '(Along principal axes of moments of inertia)'
         IF (LARGE) THEN
            WRITE (LUPRI,'(15X,3(18X,A2))') (LAB(I),I=1,NCOMP)
         ELSE
            WRITE (LUPRI,'(15X,3(10X,A2))') (LAB(I),I=1,NCOMP)
         END IF
         WRITE (LUPRI, '(/)')
         DO 300 I = 1, NCOMP
            IF (LARGE) THEN
               WRITE (LUPRI, '(12X,A2,3X,3F20.12)') LAB(I),
     &                  (FAC*AMAT(IND(I),IND(J)),J=1,NCOMP)
            ELSE
               WRITE (LUPRI, '(12X,A2,3X,3F12.6)') LAB(I),
     &                  (FAC*AMAT(IND(I),IND(J)),J=1,NCOMP)
            END IF
  300    CONTINUE
csonia 04/10/95
CMI ... LUCME not defined in Dirac!
CMI      IF (LUCME.GT.0) THEN
CMI         DO I = 1, NCOMP
CMI            WRITE (LUCME, '(1X,A2,1X,3D23.15)') LAB(I),
CMI  &                  (FAC*AMAT(IND(I),IND(J)),J=1,NCOMP)
CMI         END DO
CMI      END IF
csonia 04/10/95
      ELSE
         WRITE (LUPRI,'(2X,A)') ' Polarizabilities not calculated - '//
     &                       ' appropriate symmetries not requested.'
      END IF
      WRITE (LUPRI, '(/)')
      RETURN
      END
