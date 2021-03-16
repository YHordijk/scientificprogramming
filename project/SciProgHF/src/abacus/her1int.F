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

      SUBROUTINE PR1DRV(SOINT,NELMNT,WORK,LWORK,NPQUAD,LABINT,
     &                  INTTYP,INTREP,NOPTYP,NBAST,ANTI,IORDER,DOATOM,
     &                  INTADR,TRIANG,NATOM,SQUARE,IPRINT,DOINT)
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL ANTI, DOATOM(*), TRIANG, SQUARE, DOINT(2,2)
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION SOINT(NELMNT,NOPTYP), INTREP(NOPTYP),
     &          WORK(LWORK), INTADR(*)
C
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(A,/)')
     &     ' <<<<<<<<<< OUTPUT FROM PR1DRV >>>>>>>>>> '
      npoint = 1

      if      (inttyp == 12) then ! DSO
        npoint = npquad
        nshint = 9*natom*natom
      else if (inttyp == 101) then
!       g1n
        nshint = 3*noptyp
      else if (inttyp == 52) then
!       dplgra
        nshint = 18
      else if (inttyp == 53) then
!       quagra
        nshint = 36
      else
        nshint = noptyp
      end if

      KWEIGH = 1
      KABSCI = KWEIGH + NPOINT
      KLAST  = KABSCI + NPOINT
      IF (KLAST .GT. LWORK) CALL STOPIT('PR1DRV',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL PR1DR1(SOINT,WORK(KLAST),LWRK,NPOINT,LABINT,
     &            INTTYP,INTREP,NOPTYP,NBAST,NELMNT,ANTI,IORDER,DOATOM,
     &            WORK(KWEIGH),WORK(KABSCI),TRIANG,NATOM,INTADR,NSHINT,
     &            SQUARE,DOINT,IPRINT)
      RETURN
      END
C  /* Deck pr1dr1 */
      SUBROUTINE PR1DR1(SOINT,WORK,LWORK,NPOINT,LABINT,INTTYP,
     &                  INTREP,NOPTYP,NBAST,NELMNT,ANTI,IORDER,DOATOM,
     &                  WEIGHT,ABSCIS,TRIANG,NATOM,INTADR,NSHINT,SQUARE,
     &                  DOINT,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "pi.h"
#include "maxaqn.h"
      PARAMETER (D1 = 1.00 D00, D4 = 4.D0,D8 = 8.D0,DMTRA = .75D0)
C
      LOGICAL ANTI, DOATOM(*), HER, TRIANG, SQUARE, DONUC1, DOMOM1,
     &        DOINT(2,2)
C
      CHARACTER LABINT(NOPTYP)*8
C
      DIMENSION SOINT(NELMNT,NOPTYP), WORK(LWORK), INTREP(NOPTYP),
     &          WEIGHT(NPOINT), ABSCIS(NPOINT), INTADR(*)
C
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "onecom.h"
#include "orgcom.h"
#include "lmns.h"
#include "ccom.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "symind.h"
#include "huckel.h"
C
#include "ibtfun.h"
C
      TOLS   = THRS*THRS
      TOLOG  = - LOG(TOLS)
C
C   JMAXD  max j in (d/dx)**j
C   JMAXM : max j in ( x  )**j
C   DONUC1: true  - do < (d/dx)**j A | ... | B>
C           false - do <           A | ... | (d/dx)**j B>
C   DOMOM1: true  - do < ( x  )**j A | ... | B>
C           false - do <           A | ... | ( x  )**j B>
C   (comment added Nov. 2002 / HJAaJ)
C
      IF (INTTYP .EQ. 1 .OR. INTTYP .EQ. 45) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 2) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 3) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 4) THEN
         JMAXD  = 0
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 5) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 6) THEN
         JMAXD  = 0
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 7) THEN
         JMAXD  = 0
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 8) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 9) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.10) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.11) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.12) THEN ! DSO
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.13) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.14) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.15) THEN
         JMAXD  = 0
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.16) THEN
         JMAXD  = 0
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.17) THEN
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.18) THEN
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.19) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.20) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.21) THEN
         JMAXD  = 2
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.22) THEN
         JMAXD = 0
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.23) THEN
         JMAXD = 1
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.24) THEN
         JMAXD = 2
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.25) THEN
         JMAXD = 2
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.26) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.27) THEN
         JMAXD = 1
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.28) THEN
         JMAXD = 1
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.29 .OR. INTTYP.EQ.58) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.30) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.31) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.32) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.33) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.34) THEN
         JMAXD = 1
         JMAXM = 1
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.35) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1 = .FALSE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.36) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.37) THEN
         JMAXD = 0
         JMAXM = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.38) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.39) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.40) THEN
         JMAXD = 4
         JMAXM = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.41) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.42) THEN
         JMAXD = 0
         JMAXM = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.43) THEN
         JMAXD = 0
         JMAXM = 3
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.44) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
C     type 45 already defined
      ELSE IF (INTTYP .EQ.46) THEN
         JMAXD  = 0
         JMAXM  = 1 + IORDER
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.47) THEN
         JMAXD  = 0
         JMAXM  = 2 + IORDER
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 48) THEN
         JMAXD  = 0
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.49) THEN
         JMAXD  = 0
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.50) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSEIF (INTTYP .EQ.51) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSEIF (INTTYP .EQ.52) THEN
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSEIF (INTTYP .EQ.53) THEN
         JMAXD  = 1
         JMAXM  = 2
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ. 54) THEN
         JMAXD  = 0
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 55) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.56) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.57) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
C     type 58 already defined (NUCFIEL)
      ELSE IF (INTTYP .EQ. 62) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1 = .FALSE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.65) THEN
         JMAXD = 0
         JMAXM = 3
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.73) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.

      else if (inttyp == 100) then
!        g1o
         jmaxd  = 1
         jmaxm  = 0
         donuc1 = .true.
         domom1 = .false.
      else if (inttyp == 101) then
!        g1n
         jmaxd  = 1
         jmaxm  = 0
         donuc1 = .true.
         domom1 = .false.
      else if (inttyp == 102) then
!        g1b
         jmaxd  = 1
         jmaxm  = 0
         donuc1 = .true.
         domom1 = .false.
      else if (inttyp == 103) then
!        g1kx
         jmaxd  = 2
         jmaxm  = 0
         donuc1 = .true.
         domom1 = .true.
      else if (inttyp == 104) then
!        g1ky
         jmaxd  = 2
         jmaxm  = 0
         donuc1 = .true.
         domom1 = .true.
      else if (inttyp == 105) then
!        g1kz
         jmaxd  = 2
         jmaxm  = 0
         donuc1 = .true.
         domom1 = .true.
      else if (inttyp == 106) then
!        cxikr
         jmaxd  = 0
         jmaxm  = 0
         donuc1 = .true.
         domom1 = .true.
C
C hjaaj: INTTYP .gt. 200 are special relativistic operators for Dirac ...
C
      ELSE IF (INTTYP .EQ.201) THEN ! RM1H3
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 202) THEN ! RM1RN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 203) THEN ! RM1H2
         !JMAXD  = 0
         JMAXD  = 2 ! mi fix of out-of-bounds later on in ODC(..,..,..,here,..)
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 204) THEN ! RNST
         JMAXD = 0
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 205) THEN ! RDSUSNL - the same as DSUSNL !
         JMAXD = 0
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 206) THEN ! RDSUSLL
         JMAXD = 0
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .FALSE. ! false: take JMAXM moments to the right (nucleus B)
      ELSE IF (INTTYP .EQ. 207) THEN ! RM2H3
         JMAXD = 1
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE. ! Moment to left - does not commute with d/dq
      ELSE IF (INTTYP .EQ. 208) THEN ! RM2H2
!CMI/march '03 ... obviously should be JMAXD=0, but is gives unzero XY component !
         JMAXD = 2 ! Test...
         JMAXM = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.304) THEN
         JMAXD = 0
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.305) THEN
         JMAXD = 0
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.306) THEN
         JMAXD = 0
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.

      ELSE
         CALL QUIT('PR1DR1: Unknown INTTYP')
      END IF
C
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(A,I5)') ' NBAST  ', NBAST
         WRITE (LUPRI,'(A,I5)') ' NELMNT ', NELMNT
         WRITE (LUPRI,'(A,I5)') ' INTTYP ', INTTYP
         WRITE (LUPRI,'(A,I5)') ' NOPTYP ', NOPTYP
         WRITE (LUPRI,'(A,L5)') ' ANTI   ', ANTI
         WRITE (LUPRI,'(A,L5)') ' SQUARE ', SQUARE
         WRITE (LUPRI,'(A,(10I5))') ' INTREP ', (INTREP(I),I=1,NOPTYP)
         WRITE (LUPRI,'(A,(10(1X,A)))')' LABINT ',(LABINT(I),I=1,NOPTYP)
         WRITE (LUPRI,'(A,3F20.10)') ' ORIGIN ', (ORIGIN(I),I=1,3)
         WRITE (LUPRI,'(A,3F20.10)') ' GAUGE  ', (GAGORG(I),I=1,3)
         WRITE (LUPRI,'(A,3F20.10)') ' DIPOLE ', (DIPORG(I),I=1,3)
         WRITE (LUPRI,'(A,3F20.10)') ' CAVITY ', (CAVORG(I),I=1,3)
         WRITE (LUPRI,'(A,3F20.10)') ' CM     ', (CMXYZ(I),I=1,3)
         WRITE (LUPRI,'(A,I5)') ' JMAXD  ',JMAXD
         WRITE (LUPRI,'(A,I5)') ' JMAXM  ',JMAXM
      END IF
C
      CALL DZERO(SOINT,NELMNT*NOPTYP)
C
C     ****************************************************************
C     *** Prepare for Gaussian quadrature (diamagnetic spin-orbit) ***
C     ****************************************************************
C
      IF (INTTYP .EQ. 12) THEN ! DSO
C        HER = NPOINT .EQ. 5 .OR. NPOINT .EQ. 6 .OR.
C    &         NPOINT .EQ. 8 .OR. NPOINT .EQ. 10
         HER = .FALSE.
         IF (HER) THEN
            CALL GAUHER(ABSCIS,WEIGHT,NPOINT)
            DO 100 IPOINT = 1, NPOINT
               EXPD = DUMMY
               ABSCIS(IPOINT) = ABSCIS(IPOINT)**2
               WEIGHT(IPOINT) = D4*WEIGHT(IPOINT)*EXPD/SQRT(PI)
  100       CONTINUE
         ELSE
            CALL GAULEG(-D1,D1,ABSCIS,WEIGHT,NPOINT)
            DO 110 IPOINT = 1, NPOINT
               D1PA   = D1 + ABSCIS(IPOINT)
               D1MA   = D1 - ABSCIS(IPOINT)
               ABSCIS(IPOINT) = (DMTRA*D1PA/D1MA)**2
               WEIGHT(IPOINT) = D8*(DMTRA**3)*WEIGHT(IPOINT)*(D1PA**2)
     &                                           /(SQRT(PI)*(D1MA**4))
 110        CONTINUE
         END IF
      END IF
C
C     ************************************************************
C     ***** Triangular loop over symmetry independent shells *****
C     ************************************************************
C
      IBCNT = 1
      IF (INTTYP .EQ. 45) IBCNT = 2
      IORBA = 0
      CALL IZERO(ISOFRA, 8)
      IF (DOHUCKEL .AND. INTTYP .EQ. 1) THEN
         KMAXT = NSMLSH + NLRGSH
      ELSE
         KMAXT = KMAX
      END IF
      DO 200 ISHELA = 1,KMAXT
         NHKTA = NHKT(ISHELA)
         KHKTA = KHKT(ISHELA)
         KCKTA = KCKT(ISHELA)
         ICA   = LCLASS(ISHELA)
         SPHRA = SPHR(ISHELA)
         CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
         NCENTA = NCENT(ISHELA)
         ICENTA = NUCNUM(NCENTA,1)
         MULA   = ISTBAO(ISHELA)
         MULTA  = MULT(MULA)
         NUCA   = NUCO(ISHELA)
         NUMCFA = NUMCF(ISHELA)
         JSTA   = JSTRT(ISHELA)
         CORAX  = CENT(ISHELA,1,1)
         CORAY  = CENT(ISHELA,2,1)
         CORAZ  = CENT(ISHELA,3,1)
C
C        Compute symmetry integral pointers for contributions
C        from this block.
C
         IF (.NOT.SQUARE) THEN
            DO 210 I = 1, 8
               ISOFRB(I) = 0
               DO 220 J = 1, MXAQN
                 INDFA(I,J) = -10 000 000
 220           CONTINUE
 210        CONTINUE
            DO 230 NA = 1, KHKTA
               DO 240 IREP = 0, MAXREP
               IF (IBTAND(MULA,IBTXOR(IREP,ISYMAO(NHKTA,NA))).EQ.0) THEN
                  ISOFRA(IREP+1)    = ISOFRA(IREP+1) + 1
                  INDFA (IREP+1,NA) = ISOFRA(IREP+1)
               END IF
 240           CONTINUE
 230        CONTINUE
            IF (IPRINT .GT. 20) THEN
               WRITE(LUPRI,'(A,I4)')' IA offsets for shell ',ISHELA
               DO 250 NA = 1,KHKTA
                  WRITE(LUPRI,'(8(1X,I5))') (INDFA(I,NA),I=1,MAXREP+1)
 250           CONTINUE
            END IF
         END IF
         IORBB = 0
         ISHMXB = ISHELA
         IF (SQUARE) ISHMXB = KMAXT
         DO 300 ISHELB = 1, ISHMXB
            NHKTB = NHKT(ISHELB)
            KHKTB = KHKT(ISHELB)
            KCKTB = KCKT(ISHELB)
            ICB   = LCLASS(ISHELB)
            SPHRB = SPHR(ISHELB)
            CALL LMNVAL(NHKTB,KCKTB,LVALUB,MVALUB,NVALUB)
            NCENTB = NCENT(ISHELB)
            MULB   = ISTBAO(ISHELB)
            MULTB  = MULT(MULB)
            NUCB   = NUCO(ISHELB)
            NUMCFB = NUMCF(ISHELB)
            JSTB   = JSTRT(ISHELB)
            CORBX0 = CENT(ISHELB,1,IBCNT)
            CORBY0 = CENT(ISHELB,2,IBCNT)
            CORBZ0 = CENT(ISHELB,3,IBCNT)
            KHKTAB = KHKTA*KHKTB
            KCKTAB = KCKTA*KCKTB
            MAB    = IBTOR(MULA,MULB)
            KAB    = IBTAND(MULA,MULB)
            HKAB   = FMULT(KAB)
C
            SPHRAB = SPHRA .OR. SPHRB
C
C           Compute symmetry integral pointers.
C
            IF (.NOT.SQUARE) THEN
               DO 310 I = 1, 8
                  DO 320 J = 1, MXAQN
                    INDFB(I,J) = -10 000 000
 320              CONTINUE
 310           CONTINUE
               DO 330 NB = 1, KHKTB
               DO 340 IREP = 0, MAXREP
               IF (IBTAND(MULB,IBTXOR(IREP,ISYMAO(NHKTB,NB))).EQ.0) THEN
                  ISOFRB(IREP+1)    = ISOFRB(IREP+1) + 1
                  INDFB (IREP+1,NB) = ISOFRB(IREP+1)
               END IF
 340           CONTINUE
 330           CONTINUE
               IF (IPRINT .GT. 20) THEN
                  WRITE(LUPRI,'(A,I4)')' IB offsets for shell ',ISHELB
                  DO 350 NB = 1, KHKTB
                     WRITE(LUPRI,'(8(1X,I5))')(INDFB(I,NB),I=1,MAXREP+1)
 350              CONTINUE
               END IF
            END IF
            IF(.NOT.DOINT(ICA,ICB)) GOTO 300
            IF (IPRINT .GT. 05) WRITE (LUPRI, 1000) ISHELA, ISHELB
            IF (IPRINT .GE. 10) THEN
                WRITE (LUPRI,'(A,2I10)') ' NHKT   ', NHKTA, NHKTB
                WRITE (LUPRI,'(A,2I10)') ' KHKT   ', KHKTA, KHKTB
                WRITE (LUPRI,'(A,2I10)') ' KCKT   ', KCKTA, KCKTB
                WRITE (LUPRI,'(A,2I10)') ' NCENT  ', NCENTA, NCENTB
                WRITE (LUPRI,'(A,2I10)') ' ISTBAO ', MULA, MULB
                WRITE (LUPRI,'(A,2I10)') ' MULT   ', MULTA, MULTB
                WRITE (LUPRI,'(A,2I10)') ' NUC    ', NUCA, NUCB
                WRITE (LUPRI,'(A,2I10)') ' NUMCF  ', NUMCFA, NUMCFB
                WRITE (LUPRI,'(A,2I10)') ' JST    ', JSTA, JSTB
                WRITE (LUPRI,'(A,2F12.6)') ' CORAX    ', CORAX, CORBX0
                WRITE (LUPRI,'(A,2F12.6)') ' CORAY    ', CORAY, CORBY0
                WRITE (LUPRI,'(A,2F12.6)') ' CORAZ    ', CORAZ, CORBZ0
            END IF
C
C           ******************************
C           ***** Symmetry integrals *****
C           ******************************
C
            CALL PR1SOP(SOINT,WORK,LWORK,NPOINT,LABINT,
     &                  INTTYP,INTREP,NOPTYP,NBAST,NELMNT,ANTI,IORDER,
     &                  DOATOM,WEIGHT,ABSCIS,TRIANG,NATOM,INTADR,NSHINT,
     &                  SQUARE,IPRINT,JMAXD,JMAXM,ISHELA,ISHELB,IORBA,
     &                  IORBB,KAB,CORBX0,CORBY0,CORBZ0,TOLOG,TOLS,
     &                  DONUC1,DOMOM1,MULTA,MULTB)
C
  300    IORBB = IORBB + KHKTB
         IORBA = IORBA + KHKTA
  200 CONTINUE
C
C     ****************************
C     ***** Unpack integrals *****
C     ****************************
C
      IF (.NOT.SQUARE) THEN
         IF (NELMNT.GT.LWORK)CALL STOPIT('PR1DR1','SYMUPK',NELMNT,LWORK)
         DO 800 I = 1, NOPTYP
            CALL DCOPY(NELMNT,SOINT(1,I),1,WORK,1)
            CALL SYMUPK(WORK,SOINT(1,I),INTREP(I) + 1,NELMNT)
  800    CONTINUE
      END IF
      RETURN
 1000 FORMAT (//,2X,'***************************************',
     *         /,2X,'********** ISHELA/B =',I3,',',I3,' **********',
     *         /,2X,'***************************************',/)
      END
C  /* Deck pr1sop */
      SUBROUTINE PR1SOP(SOINT,WORK,LWORK,NPOINT,LABINT,
     &                  INTTYP,INTREP,NOPTYP,NBAST,NELMNT,ANTI,IORDER,
     &                  DOATOM,WEIGHT,ABSCIS,TRIANG,NATOM,INTADR,NSHINT,
     &                  SQUARE,IPRINT,JMAXD,JMAXM,ISHELA,ISHELB,IORBA,
     &                  IORBB,KAB,CORBX0,CORBY0,CORBZ0,TOLOG,TOLS,
     &                  DONUC1,DOMOM1,MULTA,MULTB)

#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
      LOGICAL ANTI, DOATOM(NUCIND), TRIANG, SQUARE, DONUC1, DOMOM1
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION SOINT(NELMNT,NOPTYP), WORK(LWORK), INTREP(NOPTYP),
     &          WEIGHT(NPOINT), ABSCIS(NPOINT), INTADR(*)
#include "onecom.h"
#include "nuclei.h"
C
      KSHINT = 1
      KSYMC  = KSHINT + KCKTAB*NSHINT
      KCENTC = KSYMC  + (NUCDEP + 1)/IRAT
      KNCNTC = KCENTC + (NUCDEP + 1)/IRAT
      KFACTR = KNCNTC + (NUCDEP + 1)/IRAT
      KCORCX = KFACTR + NUCDEP
      KCORCY = KCORCX + NUCDEP
      KCORCZ = KCORCY + NUCDEP
      KGEXP  = KCORCZ + NUCDEP
      KLAST  = KGEXP  + NUCDEP
      IF (KLAST .GT. LWORK) CALL STOPIT('PR1SOP',' ',KLAST,LWORK)
      LWRK = LWORK - KLAST + 1
      CALL PR1SO1(SOINT,WORK(KSHINT),WORK(KLAST),LWRK,
     &            NPOINT,LABINT,INTTYP,INTREP,NOPTYP,NBAST,NELMNT,ANTI,
     &            IORDER,DOATOM,WEIGHT,ABSCIS,TRIANG,NATOM,INTADR,
     &            NSHINT,SQUARE,IPRINT,WORK(KFACTR),JMAXD,JMAXM,ISHELA,
     &            ISHELB,IORBA,IORBB,KAB,CORBX0,CORBY0,CORBZ0,
     &            WORK(KCORCX),WORK(KCORCY),WORK(KCORCZ),WORK(KCENTC),
     &     WORK(KSYMC),WORK(KGEXP),TOLOG,TOLS,DONUC1,DOMOM1,MULTA,MULTB,
     &            WORK(KNCNTC))
      RETURN
      END
C  /* Deck pr1so1 */
      SUBROUTINE PR1SO1(SOINT,SHLINT,WORK,LWORK,NPOINT,
     &                  LABINT,INTTYP,INTREP,NOPTYP,NBAST,NELMNT,ANTI,
     &                  IORDER,DOATOM,WEIGHT,ABSCIS,TRIANG,NATOM,INTADR,
     &                  NSHINT,SQUARE,IPRINT,FACINT,JMAXD,JMAXM,ISHELA,
     &                  ISHELB,IORBA,IORBB,KAB,CORBX0,CORBY0,CORBZ0,
     &                  CORCX,CORCY,CORCZ,JCENTC,JSYMC,GEXP,TOLOG,TOLS,
     &                  DONUC1,DOMOM1,MULTA,MULTB,ncentc)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
      LOGICAL FULMAT, ANTI, DOATOM(NUCIND), TRIANG, SQUARE, DONUC1,
     &        MULCHA, DOMOM1, HUCORB
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION SOINT(NELMNT,NOPTYP), SHLINT(KCKTAB,NSHINT),
     &          WORK(LWORK), INTREP(NOPTYP),
     &          WEIGHT(NPOINT), ABSCIS(NPOINT), FACINT(*),
     &          CORCX(*), CORCY(*), CORCZ(*),
     &          INTADR(*), JSYMC(*), JCENTC(*), GEXP(*)
      integer :: ncentc(*)
C
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "onecom.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "huckel.h"
C
#include "ibtfun.h"
      XYZ(I,J) = PT(IBTAND(ISYMAX(I,1),J))
C
C     Initialization for Coulomb integrals
C
      NATOMC = 0
      IF ((INTTYP .EQ. 5) .OR. (INTTYP .EQ. 10)
     &                    .OR. (INTTYP .EQ. 11)
     &                    .OR. (INTTYP .EQ. 12) ! DSO
     &                    .OR. (INTTYP .EQ. 13)
     &                    .OR. (INTTYP .EQ. 19)
     &                    .OR. (INTTYP .EQ. 20)
     &                    .OR. (INTTYP .EQ. 24)
     &                    .OR. (INTTYP .EQ. 25)
     &                    .OR. (INTTYP .EQ. 26)
     &                    .OR. (INTTYP .EQ. 27)
     &                    .OR. (INTTYP .EQ. 28)
     &                    .OR. (INTTYP .EQ. 29) .OR. (INTTYP .EQ. 58)
     &                    .OR. (INTTYP .EQ. 30)
     &                    .OR. (INTTYP .EQ. 31)
     &                    .OR. (INTTYP .EQ. 35)
     &                    .OR. (INTTYP .EQ. 38)
     &                    .OR. (INTTYP .EQ. 48)
     &                    .OR. (INTTYP .EQ. 51)
     &                    .OR. (INTTYP .EQ. 203) ! RM1H2
     &                    .OR. (INTTYP .EQ. 204)
     &                    .OR. (INTTYP .EQ. 208) ! RM2H2
     &                    .OR. (INTTYP .EQ. 54)
     &                    .OR. (INTTYP .EQ. 55)
     &                    .OR. (INTTYP .EQ. 56)
     &                    .OR. (INTTYP .EQ. 57)
     &                    .or. (inttyp == 101)
     &                    .OR. (INTTYP .EQ. 62)) THEN
         JMAX = NHKTA + NHKTB - 2
         IF (INTTYP .EQ.  5) JMAX = JMAX + 1
         IF (INTTYP .EQ. 10) JMAX = JMAX + 1
         IF (INTTYP .EQ. 11) JMAX = JMAX + 2
         IF (INTTYP .EQ. 12) JMAX = JMAX + 2 ! DSO
         IF (INTTYP .EQ. 13) JMAX = JMAX + 2
         IF (INTTYP .EQ. 19) JMAX = JMAX + 1
         IF (INTTYP .EQ. 20) JMAX = JMAX + 1
         IF (INTTYP .EQ. 24) JMAX = JMAX + 2
         IF (INTTYP .EQ. 25) JMAX = JMAX + 2
         IF (INTTYP .EQ. 26) JMAX = JMAX + 2
         IF (INTTYP .EQ. 27) JMAX = JMAX + 2
         IF (INTTYP .EQ. 28) JMAX = JMAX + 2
         IF (INTTYP .EQ. 29) JMAX = JMAX + 1
         IF (INTTYP .EQ. 30) JMAX = JMAX + 2
         IF (INTTYP .EQ. 31) JMAX = JMAX + 2
         IF (INTTYP .EQ. 38) JMAX = JMAX + 2
         IF (INTTYP .EQ. 48) JMAX = JMAX + 2
         IF (INTTYP .EQ. 51) JMAX = JMAX + 1
         IF (INTTYP .EQ. 54) JMAX = JMAX + 1
         IF (INTTYP .EQ. 55) JMAX = JMAX + 2
         IF (INTTYP .EQ. 56) JMAX = JMAX + 4
         IF (INTTYP .EQ. 57) JMAX = JMAX + 4
         IF (INTTYP .EQ. 58) JMAX = JMAX + 1
         if (inttyp == 101) jmax = jmax + 1
         IF (INTTYP .EQ. 203) JMAX = JMAX + 1 ! RM1H2
         IF (INTTYP .EQ. 204) JMAX = JMAX + 2 ! RNST
         IF (INTTYP .EQ. 208) JMAX = JMAX + 2 ! RM2H2
         MULCHA = (INTTYP .EQ. 5) .OR. (INTTYP .EQ. 19)
     &                            .OR. (INTTYP .EQ. 20)
     &                            .OR. (INTTYP .EQ. 24)
     &                            .OR. (INTTYP .EQ. 25)
     &                            .OR. (INTTYP .EQ. 54)
     &                            .OR. (INTTYP .EQ. 58)
     &                            .or. (inttyp == 101)
     &                            .OR. (INTTYP .EQ. 203) ! RM1H2
     &                            .OR. (INTTYP .EQ. 208) ! RM2H2
         ISTEPU = JMAX + 1
         ISTEPV = ISTEPU*ISTEPU
         NAHGTF = ISTEPU*ISTEPV
         IF (INTTYP .NE. 12) THEN ! not DSO
            NATOMC = 0
            DO 300 IATOMC = 1, NUCIND
            IF (DOATOM(IATOMC)) THEN
               MULC   = ISTBNU(IATOMC)
               MABC   = IBTOR(MULC,KAB)
               CORCX0 =  CORD(1,IATOMC)
               CORCY0 =  CORD(2,IATOMC)
               CORCZ0 =  CORD(3,IATOMC)
               FACTOR = - FMULT(IBTAND(MULC,KAB))
               IF (MULCHA) FACTOR = FACTOR*CHARGE(IATOMC)/HKAB
               DO 310 ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MABC) .EQ. 0) THEN
                     NATOMC = NATOMC + 1
                     JSYMC(NATOMC)  = ISYMOP
                     JCENTC(NATOMC) = IATOMC
                     CORCX(NATOMC)  = XYZ(1,ISYMOP)*CORCX0
                     CORCY(NATOMC)  = XYZ(2,ISYMOP)*CORCY0
                     CORCZ(NATOMC)  = XYZ(3,ISYMOP)*CORCZ0
                     GEXP(NATOMC)   = GNUEXP(IATOMC)
                     FACINT(NATOMC) = FACTOR
                     NCENTC(NATOMC) = NUCNUM(IATOMC,ISYMOP+1)
                  END IF
  310          CONTINUE
            END IF
  300       CONTINUE
         END IF
      END IF
C
C     *****************************************
C     ***** Loop over symmetry operations *****
C     *****************************************
C

      DO 100 ISYMOP = 0, MAXOPR
      IF(IBTAND(ISYMOP,MAB) .EQ. 0) THEN
         LDIAG  = (ISHELA .EQ. ISHELB) .AND. (.NOT.SQUARE)
         ONECEN = .FALSE.
         CORBX  = XYZ(1,ISYMOP)*CORBX0
         CORBY  = XYZ(2,ISYMOP)*CORBY0
         CORBZ  = XYZ(3,ISYMOP)*CORBZ0
         ICENTB = NUCNUM(NCENTB,ISYMOP+1)
         DIFABX = CORAX - CORBX
         DIFABY = CORAY - CORBY
         DIFABZ = CORAZ - CORBZ
         DISTAB = DIFABX*DIFABX + DIFABY*DIFABY + DIFABZ*DIFABZ
C
C     Scale overlap integrals for the Huckel matrix part
C
         HUCFAC = 1.0D0
         IF (LCLASS(ISHELA) .EQ. 2 .AND. LCLASS(ISHELB) .EQ. 2 .AND.
     &       INTTYP .EQ. 1 .AND. DOHUCKEL) THEN
            IHUCA = IHUCPT(ISHELA)
            IF (ISHELA .EQ. ISHELB .AND. DISTAB .LT. 1.0D-3) THEN
               HUCFAC = HUCEXC(IHUCA)
            ELSE
               IHUCB = IHUCPT(ISHELB)
               HUCFAC = 0.50D0*HUCCNT*(HUCEXC(IHUCA) + HUCEXC(IHUCB))
            END IF
            SCAL   = DFLOAT(MULTA*MULTB)
            HUCFAC = HUCFAC/SQRT(SCAL)
         END IF
         IF (IPRINT .GT. 05) THEN
            WRITE (LUPRI, 1000) ISYMOP
            IF (IPRINT .GE. 10) THEN
               WRITE (LUPRI,'(A,1F12.6)') ' CORBX    ', CORBX
               WRITE (LUPRI,'(A,1F12.6)') ' CORBY    ', CORBY
               WRITE (LUPRI,'(A,1F12.6)') ' CORBZ    ', CORBZ
            END IF
         END IF
         CALL DZERO(SHLINT,KCKTAB*NSHINT)
C
C        *****************************************
C        ***** Calculate Cartesian integrals *****
C        *****************************************
C
         CALL PR1PRM(SHLINT,WORK,LWORK,LABINT,INTTYP,
     &               NOPTYP,NBAST,IORDER,DOATOM,IPRINT,TOLS,TOLOG,
     &               NATOMC,DISTAB,FACINT,CORCX,CORCY,CORCZ,JMAXD,
     &               JMAXM,DIFABX,DIFABY,DIFABZ,NPOINT,WEIGHT,ABSCIS,
     &               GEXP,TRIANG,KAB,DONUC1,DOMOM1,HUCFAC,ncentc)
C
C        ****************************************
C        ***** Transform to spherical basis *****
C        ****************************************
C
         IF (SPHRAB) THEN
            CALL SPHRM1(SHLINT,SHLINT,NSHINT,WORK,LWORK,.FALSE.,IPRINT)
         END IF
C
C        *************************************************
C        ***** Transform integrals to symmetry basis *****
C        *************************************************
C
         IF (INTTYP .EQ. 10) THEN
            CALL PSOTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 11 .OR. INTTYP .EQ. 13) THEN
            CALL SDTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                 NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 12) THEN ! DSO
            CALL DSOTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,DOATOM,
     &                  KAB,TRIANG,NATOM,INTADR,NSHINT,IPRINT)
         ELSE IF ((INTTYP .EQ. 14) .OR. (INTTYP .EQ. 44)) THEN
            CALL HDOTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,IORBA,IORBB,
     &                  NBAST,INTTYP,IPRINT)
         ELSE IF ((INTTYP .EQ. 26) .OR. (INTTYP .EQ. 27) .OR.
     &            (INTTYP .EQ. 28) .OR. (INTTYP .EQ. 38) .OR.
     &            (INTTYP .EQ. 48) .OR. (INTTYP .EQ. 55) .OR.
     &            (INTTYP .EQ. 204) ) THEN
            CALL NSTTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,IORBA,IORBB,NBAST,INTTYP,
     &                  IPRINT,INTADR)
         ELSE IF (INTTYP .EQ. 29) THEN
            CALL EFTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 30 .OR. INTTYP .EQ. 31) THEN
            CALL EFGTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 34) THEN
            CALL HDBTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,IORBA,IORBB,
     &                  NBAST,INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 35) THEN
            CALL NPETRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 62) THEN
            CALL PVCTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)

         ELSE IF (INTTYP .EQ. 52) THEN
            CALL DPGTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 53) THEN
            CALL QUGTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)

         ELSE IF (INTTYP .EQ. 56) THEN
            CALL EFTTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 57) THEN
            CALL EFTTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
         else if (inttyp == 100) then
!          g1o
           call s1htra(shlint, soint, isymop, 
     &                 nelmnt, noptyp, anti, iprint)
         else if (inttyp == 101) then
!          g1n
           call g1htra(shlint,soint,facint,jsymc,jcentc,ncentc,natomc,
     &                 isymop,nelmnt,noptyp,anti,intadr,iprint)
         else if (inttyp == 102) then
!          g1b
           call s1htra(shlint, soint, isymop, 
     &                 nelmnt, noptyp, anti, iprint)
         else if (inttyp == 103) then
!          g1kx
           call s1htra(shlint, soint, isymop, 
     &                 nelmnt, noptyp, anti, iprint)
         else if (inttyp == 104) then
!          g1ky
           call s1htra(shlint, soint, isymop, 
     &                 nelmnt, noptyp, anti, iprint)
         else if (inttyp == 105) then
!          g1kz
           call s1htra(shlint, soint, isymop, 
     &                 nelmnt, noptyp, anti, iprint)

         ELSE
            FULMAT = .TRUE.
            DO 200 I = 1, NOPTYP
               IF (SQUARE) THEN
                  CALL SYMSQR(SHLINT(1,I),SOINT(1,I),INTREP(I),ISYMOP,
     &                        IORBA,IORBB,-HKAB,NBAST,IPRINT)
               ELSE IF (INTREP(I) .EQ. 0) THEN
                  CALL SYM1S(SHLINT(1,I),SOINT(1,I),ISYMOP,MULA,MULB,
     &                       NHKTA,NHKTB,KHKTA,KHKTB,-HKAB,LDIAG,
     &                       FULMAT,DUM,IDUM,IPRINT)
               ELSE
                  CALL SYM1N(SHLINT(1,I),SOINT(1,I),INTREP(I),ISYMOP,
     &                       MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,-HKAB,
     &                       LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
               END IF
  200       CONTINUE
         END IF
      END IF
  100 CONTINUE
      RETURN
 1000 FORMAT (//,2X,'***************************************',
     &         /,2X,'******** Symmetry operation ',I2,' ********',
     &         /,2X,'***************************************',/)
      END
C  /* Deck pr1prm */
      SUBROUTINE PR1PRM(SHLINT,WORK,LWORK,LABINT,INTTYP,
     &                  NOPTYP,NBAST,IORDER,DOATOM,IPRINT,TOLS,TOLOG,
     &                  NATOMC,DISTAB,FACINT,CORCX,CORCY,CORCZ,JMAXD,
     &                  JMAXM,DIFABX,DIFABY,DIFABZ,NPOINT,WEIGHT,ABSCIS,
     &                  GEXP,TRIANG,KAB,DONUC1,DOMOM1,HUCFAC,ncentc)
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL DOATOM(*), TRIANG, DONUC1, DOMOM1, HUCORB
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION WORK(LWORK), FACINT(NATOMC), CORCX(NATOMC),
     &          CORCY(NATOMC), CORCZ(NATOMC), SHLINT(KCKTAB,*), GEXP(*)
      integer :: ncentc(*)
#include "onecom.h"
C
      JMAXA  = NHKTA - 1
      IF (INTTYP .EQ. 42) JMAXA = JMAXA + 1
      IF (INTTYP .EQ. 43) JMAXA = JMAXA + 2
      IF (INTTYP .EQ. 46) JMAXA = JMAXA + 1
      IF (INTTYP .EQ. 47) JMAXA = JMAXA + 2
      IF (INTTYP .EQ. 65) JMAXA = JMAXA + 2
      JMAXB  = NHKTB - 1

CMI     for RDSUSLL, wee need to multiply with x_N, y_N, z_N ...
      IF (INTTYP .EQ. 206) JMAXB = JMAXB + 1

      JMAXT  = JMAXA + JMAXB + JMAXD + JMAXM
C
      KODC   = 1
      KAHGTF = KODC+ 3*(JMAXA+1)*(JMAXB+1)*(JMAXT+1)*(JMAXD+1)*(JMAXM+1)
      KLAST  = KAHGTF + NAHGTF*(NATOMC + 1)
      IF (KLAST .GT. LWORK) CALL STOPIT('PR1PRM',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      IF (INTTYP .EQ. 62) THEN
C     Integrals for Parity Violation - chirality
         CALL PR1PRV(WORK(KODC),JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &               WORK(KAHGTF),SHLINT,WORK(KLAST),LWRK,
     &               LABINT,INTTYP,NOPTYP,NBAST,IORDER,DOATOM,IPRINT,
     &               TOLS,TOLOG,NATOMC,DISTAB,FACINT,CORCX,CORCY,CORCZ,
     &               DIFABX,DIFABY,DIFABZ,DONUC1,DOMOM1,HUCFAC)
      ELSE IF (INTTYP .NE. 12) THEN ! not DSO
         CALL PR1PRA(WORK(KODC),JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &               WORK(KAHGTF),SHLINT,WORK(KLAST),LWRK,
     &               LABINT,INTTYP,NOPTYP,NBAST,IORDER,DOATOM,IPRINT,
     &               TOLS,TOLOG,NATOMC,DISTAB,FACINT,CORCX,CORCY,CORCZ,
     &               DIFABX,DIFABY,DIFABZ,GEXP,DONUC1,DOMOM1,HUCFAC,
     &               ncentc)
      ELSE
C     DSO - Diamagnetic spin-orbit integrals
         CALL PR1PRD(WORK(KODC),JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &               WORK(KAHGTF),SHLINT,WORK(KLAST),LWRK,NPOINT,LABINT,
     &               INTTYP,NOPTYP,NBAST,DOATOM,WEIGHT,ABSCIS,TRIANG,
     &               IPRINT,TOLS,TOLOG,DISTAB,FACINT,CORCX,CORCY,
     &               CORCZ,DIFABX,DIFABY,DIFABZ,KAB,DONUC1,DOMOM1)
      END IF
      RETURN
      END
C  /* Deck pr1pra */
      SUBROUTINE PR1PRA(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  WORK,LWORK,LABINT,INTTYP,NOPTYP,
     &                  NBAST,IORDER,DOATOM,IPRINT,TOLS,TOLOG,NATOMC,
     &                  DISTAB,FACINT,CORCX,CORCY,CORCZ,DIFABX,DIFABY,
     &                  DIFABZ,GEXP,DONUC1,DOMOM1,HUCFAC,
     &                  ncentc)
C MI/HJAaJ - added new integrals for NMR & magnetizabilities with LAO
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "mxcent.h"
#include "pi.h"
      PARAMETER (D1 = 1.00 D00, D3 = 3.00 D00, D3INV = D1/D3,
     &           D2INV = 0.50D0, D1P5 = 1.5D0, GNUTHR = 1.D-4,
     &           ZERODIST = 1.0D-12)
C
      LOGICAL DOATOM(*), DIFODC, DONUC1, DOMOM1, ADDNAI, HUCORB, ONCN
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION WORK(LWORK), TEMP(3),
     &          FACINT(NATOMC), CORCX(NATOMC), CORCY(NATOMC),
     &          CORCZ(NATOMC), SHLINT(KCKTAB,*), AHGTF(*), GEXP(*),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
      integer :: ncentc(*)
C
#include "onecom.h"
#include "orgcom.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "nuclei.h"
#include "primit.h"
C     DIFODC true if we have JMAXD .gt. 0 (take derivatives)
      DIFODC = (INTTYP.EQ. 3) .OR. (INTTYP.EQ. 5) .OR. (INTTYP.EQ.10)
     &    .OR. (INTTYP.EQ.14) .OR. (INTTYP.EQ.17) .OR. (INTTYP.EQ.18)
     &    .OR. (INTTYP.EQ.19) .OR. (INTTYP.EQ.20) .OR. (INTTYP.EQ.21)
     &    .OR. (INTTYP.EQ.23) .OR. (INTTYP.EQ.24) .OR. (INTTYP.EQ.25)
     &    .OR. (INTTYP.EQ.27) .OR. (INTTYP.EQ.28) .OR. (INTTYP.EQ.34)
     &    .OR. (INTTYP.EQ.40) .OR. (INTTYP.EQ.44) .OR. (INTTYP.EQ.50)
     &    .OR. (INTTYP.EQ.51)
     &    .or. (inttyp ==  52)
     &    .or. (inttyp ==  53)
     &    .or. (inttyp == 100)
     &    .or. (inttyp == 101)
     &    .or. (inttyp == 102)
     &    .or. (inttyp == 103)
     &    .or. (inttyp == 104)
     &    .or. (inttyp == 105)
     &    .OR. (INTTYP.EQ.201) .OR. (INTTYP.EQ.207)

C     ADD Nuclear Attraction Integrals, true if we need V_ne operator.
      ADDNAI = (INTTYP.EQ. 5) .OR. (INTTYP.EQ.19) .OR. (INTTYP.EQ.20)
     &     .OR. (INTTYP.EQ.24) .OR. (INTTYP.EQ.25) .OR. (INTTYP.EQ.51)
     &     .OR. (INTTYP.EQ.54) .OR. (INTTYP.EQ.58)
     &     .OR. (INTTYP.EQ.203) .OR. (INTTYP.EQ.208)
C
      DO 100 IPRIMA = 1,NUCA
         JPRIMA = JSTA + IPRIMA
         CONTA = PRICCF(JPRIMA,NUMCFA)
         EXPA = PRIEXP(JPRIMA)
      DO 200 IPRIMB = 1,NUCB
         JPRIMB = JSTB + IPRIMB
         CONTB = PRICCF(JPRIMB,NUMCFB)
         EXPB = PRIEXP(JPRIMB)
         EXPP = EXPA + EXPB
         EXPPI = D1/EXPP
         EXPABQ = EXPA*EXPB*DISTAB*EXPPI
         IF (EXPABQ.GT.TOLOG) GO TO 200
         SAAB = CONTA*CONTB*EXP(-EXPABQ)
         ASAAB = ABS(SAAB)
         IF (ASAAB.LT.TOLS) GO TO 200
         SAAB13 = SIGN(ASAAB**D3INV,SAAB)
C
C        Calculate coordinates of product Gaussian P
C
         EXPAPI = EXPA*EXPPI
         EXPBPI = EXPB*EXPPI
         CORPX  = EXPAPI*CORAX + EXPBPI*CORBX
         CORPY  = EXPAPI*CORAY + EXPBPI*CORBY
         CORPZ  = EXPAPI*CORAZ + EXPBPI*CORBZ
C
C        *********************************************
C        ***** Overlap Distribution Coefficients *****
C        *********************************************
C
C        Expansion coefficients
C
         DIFPAX = CORPX - CORAX
         DIFPAY = CORPY - CORAY
         DIFPAZ = CORPZ - CORAZ
         DIFPBX = CORPX - CORBX
         DIFPBY = CORPY - CORBY
         DIFPBZ = CORPZ - CORBZ
C
C     Choose origin for expansions coefficients:
C     1) A certain nuclei
C     2) Dipole origin
C     3) Gauge origin
C     4) Molecular center of mass
C     5) Operator center ("molecular" origin)
C
         IF ((INTTYP .EQ. 17) .OR. (INTTYP .EQ. 20) .OR.
     &       (INTTYP .EQ. 22) .OR. (INTTYP .EQ. 26).OR. 
     &       (INTTYP .EQ. 202).OR. (INTTYP .EQ.205) ) THEN
            TEMP(1) = CORBX ! MI:shift the gauge to the right orbital
            TEMP(2) = CORBY
            TEMP(3) = CORBZ
         ELSE IF ((INTTYP .EQ. 2) .OR. (INTTYP .EQ. 4) .OR.
     &            (INTTYP .EQ. 6) .OR. (INTTYP .EQ. 7) .OR.
     &            (INTTYP .EQ. 8) .OR. (INTTYP .EQ. 29) .OR.
     &            (INTTYP .EQ. 30) .OR. (INTTYP .EQ. 31) .OR.
     &            (INTTYP .EQ. 53) .OR. (INTTYP .EQ. 56) .OR.
     &            (INTTYP .EQ. 57) .OR. (INTTYP .EQ. 58)) THEN
            TEMP(1) = DIPORG(1)
            TEMP(2) = DIPORG(2)
            TEMP(3) = DIPORG(3)
         ELSE IF ((INTTYP .EQ. 37) .OR. (INTTYP .EQ. 38)) THEN
            TEMP(1) = GAGORG(1)
            TEMP(2) = GAGORG(2)
            TEMP(3) = GAGORG(3)
         ELSE IF (INTTYP .EQ. 18) THEN
            TEMP(1) = CMXYZ(1)
            TEMP(2) = CMXYZ(2)
            TEMP(3) = CMXYZ(3)
         ELSE IF ((INTTYP .EQ. 46) .OR. (INTTYP .EQ. 47)) THEN
            TEMP(1) = CAVORG(1)
            TEMP(2) = CAVORG(2)
            TEMP(3) = CAVORG(3)
         ELSE IF (INTTYP .EQ. 48) THEN
C
C     We want to place the origin at the nucleus investigated
C
            DO 37 KR = 1, NUCIND
               IF (DOATOM(KR)) THEN
                  TEMP(1) = CORD(1,KR)
                  TEMP(2) = CORD(2,KR)
                  TEMP(3) = CORD(3,KR)
               END IF
 37         CONTINUE 
         ELSE
            TEMP(1) = ORIGIN(1)
            TEMP(2) = ORIGIN(2)
            TEMP(3) = ORIGIN(3)
         END IF
C
C        * NEXT IF-STATEMENT CHANGED TO INCLUDE DARWIN TERM*
C
         IF (INTTYP .NE. 9 .AND. INTTYP .NE. 41) THEN
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
         END IF
C
C        **********************************************
C        ***** Calculation of Hermitian integrals *****
C        **********************************************
C
C        Overlap integral
C
         SHGTF = SQRT(PI*EXPPI)
C
C        Nuclear attraction integrals for spin-orbit and london contribution
C        to angular momentum.
C
         IF ((INTTYP .EQ. 5) .OR. (INTTYP .EQ. 10)
     &                       .OR. (INTTYP .EQ. 11)
     &                       .OR. (INTTYP .EQ. 13)
     &                       .OR. (INTTYP .EQ. 19)
     &                       .OR. (INTTYP .EQ. 20)
     &                       .OR. (INTTYP .EQ. 24)
     &                       .OR. (INTTYP .EQ. 25)
     &                       .OR. (INTTYP .EQ. 26)
     &                       .OR. (INTTYP .EQ. 27)
     &                       .OR. (INTTYP .EQ. 28)
     &                       .OR. (INTTYP .EQ. 29) .OR. (INTTYP .EQ. 58)
     &                       .OR. (INTTYP .EQ. 30)
     &                       .OR. (INTTYP .EQ. 31)
     &                       .OR. (INTTYP .EQ. 35)
     &                       .OR. (INTTYP .EQ. 38)
     &                       .OR. (INTTYP .EQ. 48)
     &                       .OR. (INTTYP .EQ. 51)
     &                       .OR. (INTTYP .EQ. 54)
     &                       .OR. (INTTYP .EQ. 55)
     &                       .OR. (INTTYP .EQ. 56)
     &                       .OR. (INTTYP .EQ. 57)
     &                       .or. (inttyp == 101)
     &                       .OR. (INTTYP .EQ. 203) ! RM1H2
     &                       .OR. (INTTYP .EQ. 204) ! RNST
     &                       .OR. (INTTYP .EQ. 208) ! RM2H2 
     &                                             ) THEN
               CALL DZERO(AHGTF,(NATOMC + 1)*NAHGTF)
            IF (ADDNAI) THEN
               CALL DZERO(AHGTF,NAHGTF)
               IADR = 1 + NAHGTF
            ELSE
               IADR = 1
            END IF
            FACTOR = D1
            DO 300 IATOMC = 1, NATOMC
               if (addnai 
     &             .or. (inttyp == 101)) then
                 factor = facint(iatomc)
                 if (dabs(factor) < 1.0d-40) then
                   call dzero(ahgtf(iadr), nahgtf)
                   iadr = iadr + nahgtf
                   go to 300
                 end if
               end if

               DIFCPX = CORCX(IATOMC) - CORPX
               DIFCPY = CORCY(IATOMC) - CORPY
               DIFCPZ = CORCZ(IATOMC) - CORPZ
               if ((inttyp == 54 .or.
     &              inttyp == 58 .or.
!101: first geo derivatives of nuclear attraction integrals
     &              inttyp == 101 .or.
     &              inttyp == 203) .and.
     &             (abs(gexp(iatomc)) > gnuthr)) then
C                 58: NUCFIEL
                  EXPFAC = GEXP(IATOMC)
                  EXPFAC = EXPFAC/(EXPFAC + EXPP)
                  EXPPGN = EXPP*EXPFAC
                  FACTOR = FACTOR*(EXPFAC**D1P5)
                  CALL HERNAI(AHGTF,JMAX,EXPPGN,DIFCPX,DIFCPY,
     *                 DIFCPZ,FACTOR,IADR,ISTEPU,ISTEPV,NAHGTF,
     *                 IPRINT)
               ELSE
                  CALL HERNAI(AHGTF,JMAX,EXPP,DIFCPX,DIFCPY,DIFCPZ
     &                 ,FACTOR,IADR,ISTEPU,ISTEPV,NAHGTF,IPRINT)
               END IF
               
               IF (ADDNAI) THEN
                  CALL DAXPY(NAHGTF,D1,AHGTF(IADR),1,AHGTF(1),1)
               ELSE
                  IADR = IADR + NAHGTF
               END IF
  300       CONTINUE
         END IF
C
C        **********************************************
C        ***** Calculation of Cartesian integrals *****
C        **********************************************
C
C        Overlap integrals
C        -----------------
C
         IF (INTTYP .EQ. 1 .OR. INTTYP .EQ. 45) THEN
            CALL OVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,HUCFAC,
     &                  SHGTF,SHLINT)
C
C        Dipole length integrals
C        -----------------------
C
         ELSE IF (INTTYP .EQ. 2) THEN
            CALL DPLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIPORG)
C
C        Dipole velocity and half-derivative overlap integrals
C        -----------------------------------------------------
C
         ELSE IF ((INTTYP .EQ. 3) .OR. (INTTYP .EQ. 14) .OR.
     &            (INTTYP .EQ. 44)) THEN
            CALL DPVINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        Quadrupole integrals
C        --------------------
C
         ELSE IF (INTTYP .EQ. 4) THEN
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Spin-orbit integrals
C        --------------------
C
         ELSE IF (INTTYP .EQ. 5) THEN
            CALL SSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT)
C
C        Second moments integrals
C        ------------------------
C
         ELSE IF (INTTYP .EQ. 6) THEN
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Traceless theta quadrupole moments integrals
C        --------------------------------------------
C
         ELSE IF (INTTYP .EQ. 7) THEN
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Multipole moment integrals
C        --------------------------
C
         ELSE IF (INTTYP .EQ. 8) THEN
            CALL MOMINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,DIPORG,IORDER,WORK,
     &                  LWORK)
C
C        Fermi contact integrals
C        -----------------------
C
         ELSE IF (INTTYP .EQ. 9) THEN
            CALL FRMINT(SHLINT,NOPTYP,EXPP,CORPX,CORPY,CORPZ,SAAB,
     &                  DOATOM)
C
C        Paramagnetic spin-orbit integrals
C        ---------------------------------
C
         ELSE IF (INTTYP .EQ. 10) THEN
            CALL PSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C        Spin-dipole integrals
C        ---------------------
C
         ELSE IF (INTTYP .EQ. 11 .OR. INTTYP .EQ. 13) THEN
            CALL SDINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                 NOPTYP,NATOMC,INTTYP,SAAB,EXPP,CORCX,CORCY,CORCZ,
     &                 CORPX,CORPY,CORPZ)
C
C        Magnetic derivatives of overlap matrix
C        --------------------------------------
C
         ELSE IF (INTTYP .EQ. 15) THEN
            CALL SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C        Second order magnetic derivatives of overlap matrix
C        ---------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 16) THEN
            CALL SM2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C        Angular momentum around the nuclei
C        ----------------------------------
C
         ELSE IF (INTTYP .EQ. 17) THEN
            CALL AM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  D1)
C
C        Angular momentum around the molecular center of mass
C        ----------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 18) THEN
            CALL AM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  D1)
C
C        London orbital contribution to magnetic moment
C        ----------------------------------------------
C
         ELSE IF (INTTYP .EQ. 19) THEN
            CALL MM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,.TRUE.)
C
C        One-electron contribution to magnetic moment
C        --------------------------------------------
C
         ELSE IF (INTTYP .EQ. 20) THEN
            CALL AM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  D2INV)
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  ORIGIN,INTTYP)
            CALL MM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,.TRUE.)
C
C        Kinetic energy integrals
C        ------------------------
C
         ELSE IF (INTTYP .EQ. 21) THEN
            CALL KININT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        Diamagnetic susceptiblity
C        -------------------------
C
         ELSE IF (INTTYP .EQ. 22) THEN
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Angular London orbital contribution to magnetic susceptibility
C        --------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 23) THEN
            CALL DSUSAN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT,CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ)
C
C        London orbital contribution to magnetic susceptibility
C        ------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 24) THEN
            CALL DSUSLH(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,.TRUE.)
C
C        Diamagnetic susceptibility with all terms included ! - DIASUS
C      -------------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 25) THEN
            CALL DSUSAN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT,CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ)
            CALL DSUSLH(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,.TRUE.)
            TEMP(1) = CORBX
            TEMP(2) = CORBY
            TEMP(3) = CORBZ
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Nuclear shielding integrals
C        ---------------------------
C
         ELSE IF (INTTYP .EQ. 26) THEN
            CALL NSNLIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        London orbital contribution to nuclear shielding tensor
C        -------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 27) THEN
            CALL NSLOIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,DIFABX,DIFABY,DIFABZ)
C
C        Nuclear shielding tensor integrals
C        ----------------------------------
C
         ELSE IF (INTTYP .EQ. 28) THEN
            CALL NSLOIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,DIFABX,DIFABY,DIFABZ)
            TEMP(1) = CORBX
            TEMP(2) = CORBY
            TEMP(3) = CORBZ
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
            CALL NSNLIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        Electric field from the individual nuclei (type 29)
C        or from all nuclei (type 58)
C        ---------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 29 .OR. INTTYP .EQ. 58) THEN
            CALL EF1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,INTTYP)
C
C        Electric field gradient from the individual nuclei
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 30 .OR. INTTYP .EQ. 31) THEN
            CALL EFGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        Bra-differentiated overlap matrix with respect to magnetic field
C        ----------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 32) THEN
            DIFAOX = CORAX - GAGORG(1)
            DIFAOY = CORAY - GAGORG(2)
            DIFAOZ = CORAZ - GAGORG(3)
            CALL SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFAOX,DIFAOY,DIFAOZ)
C
C        Ket-differentiated overlap matrix with respect to magnetic field
C        ----------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 33) THEN
            DIFOBX = GAGORG(1) - CORBX
            DIFOBY = GAGORG(2) - CORBY
            DIFOBZ = GAGORG(3) - CORBZ
            CALL SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFOBX,DIFOBY,DIFOBZ)
C
C        Ket-differentiated HDO-integrals with respect to magnetic field
C        ---------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 34) THEN
            DIFOBX = GAGORG(1) - CORBX
            DIFOBY = GAGORG(2) - CORBY
            DIFOBZ = GAGORG(3) - CORBZ
            CALL HDBINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFOBX,DIFOBY,DIFOBZ)
C
C        Potential energy of interaction of electrons with individual nuclei
C        -------------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 35) THEN
            CALL NPEINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        Half B-differentiated overlap matrix
C        ------------------------------------
C
         ELSE IF (INTTYP .EQ. 36) THEN
            DIFOBX = CORBX - GAGORG(1)
            DIFOBY = CORBY - GAGORG(2)
            DIFOBZ = CORBZ - GAGORG(3)
            CALL SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFOBX,DIFOBY,DIFOBZ)
            DIFAOX = CORAX - GAGORG(1)
            DIFAOY = CORAY - GAGORG(2)
            DIFAOZ = CORAZ - GAGORG(3)
            CALL SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFAOX,DIFAOY,DIFAOZ)
C
C        Diamagnetic susceptiblity with common gauge origin
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 37) THEN
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Nuclear shielding integrals with common gauge origin
C        ----------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 38 .OR. INTTYP .EQ. 48) THEN
            CALL NSNLIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        Cosine and sine integrals
C        -------------------------
C
         ELSE IF (INTTYP .EQ. 39) THEN
            CALL GOSINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,ORIGIN,EXPPI)
C
C        Mass velocity integrals
C        -----------------------
C
         ELSE IF (INTTYP .EQ. 40) THEN
            CALL MVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        Darwin term integrals
C        ---------------------
C
         ELSE IF (INTTYP .EQ. 41) THEN
            CALL DWNINT(SHLINT,EXPP,CORPX,CORPY,CORPZ,SAAB)
C
C        Electric field contribution to magnetic moment
C        ----------------------------------------------
C
         ELSE IF (INTTYP .EQ. 42) THEN
            CALL CM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C        Electric field contribution to diamagnetic magnetizability
C        ----------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 43) THEN
            CALL CM2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C        First magnetic derivative of solvent integrals
C        ----------------------------------------------
C
         ELSE IF (INTTYP .EQ. 46) THEN
            CALL SL1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,WORK,LWORK)
C
C        Second magnetic derivative of solvent integrals
C        -----------------------------------------------
C
         ELSE IF (INTTYP .EQ. 47) THEN
            CALL SL2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,WORK,LWORK)
C
C        First electric derivative of overlap integrals. Type A
C        ------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 49) THEN
            CALL SE1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        First electric derivative of overlap integrals. Type B
C        ------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 50) THEN
            CALL SE1INB(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     $           SHLINT, EXPA, EXPB)
C
C        First electric deriv. of 1-electron Ham. integrals
C        --------------------------------------------------
C
         ELSEIF (INTTYP .EQ. 51) THEN
            CALL H1EINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     $           SHLINT,AHGTF,NATOMC)
C
         ELSE IF (INTTYP .EQ. 52) THEN
            CALL DPGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  SHGTF,SHLINT(1,1))
            DONUC1 = .FALSE.
            DOMOM1 = .TRUE.
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
            DONUC1 = .TRUE.
            DOMOM1 = .FALSE.
            CALL DPGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  SHGTF,SHLINT(1,10))

!radovan:   this was here but was unreachable !?
!           CALL RM1H3(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
!    &                  DIFABX,DIFABY,DIFABZ)
C
C
         ELSE IF (INTTYP .EQ. 53) THEN

            CALL QUGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  SHGTF,SHLINT(1,1))
            DONUC1 = .FALSE.
            DOMOM1 = .TRUE.
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
            DONUC1 = .TRUE.
            DOMOM1 = .FALSE.
            CALL QUGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  SHGTF,SHLINT(1,19))

!radovan:   this was here but was unreachable !?
!           CALL DPLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
!    &                  CORPX,CORPY,CORPZ,TEMP)
C
C        Potential energy integrals
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 54) THEN
            CALL RM1H2(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC)
C
C        diamagnetic contrib
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 55) THEN
            CALL RNSLOIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,DIFABX,DIFABY,DIFABZ)
C
C        Traceless electric field third derivatives
C        ----------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 56) THEN
            ONCN =   (CORAX .EQ. CORBX)
     &         .AND. (CORAY .EQ. CORBY)
     &         .AND. (CORAZ .EQ. CORBZ)
            CALL EFTINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,.TRUE.,ONCN,
     &                  CORCX,CORCY,CORCZ)
C
C        Non-traceless electric field third derivatives
C        ----------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 57) THEN
            ONCN =   (CORAX .EQ. CORBX)
     &         .AND. (CORAY .EQ. CORBY)
     &         .AND. (CORAZ .EQ. CORBZ)
            CALL EFTINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,.FALSE.,ONCN,
     &                  CORCX,CORCY,CORCZ)



C        Electric field contribution to magnetic moment
C        ----------------------------------------------
C
         ELSE IF (INTTYP .EQ. 65) THEN
            CALL QDBINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)


!     g1o: first geo derivatives of overlap integrals
!     ===============================================

         else if (inttyp == 100) then
           call dpvint(odc, 
     &                 jmaxa, 
     &                 jmaxb, 
     &                 jmaxt, 
     &                 jmaxd, 
     &                 jmaxm, 
     &                 shgtf, 
     &                 shlint)


!     g1n: first geo derivatives of nuclear attraction integrals
!     ==========================================================

         else if (inttyp == 101) then
           call g1n_integrals(jmaxa, 
     &                        jmaxb, 
     &                        jmaxt, 
     &                        jmaxd, 
     &                        jmaxm, 
     &                        odc, 
     &                        natomc, 
     &                        ncentc,
     &                        ahgtf, 
     &                        noptyp, 
     &                        shlint)


!     g1b: first geo derivatives of beta integrals
!     ============================================

         else if (inttyp == 102) then
           call dpvint(odc, 
     &                 jmaxa, 
     &                 jmaxb, 
     &                 jmaxt, 
     &                 jmaxd, 
     &                 jmaxm, 
     &                 shgtf, 
     &                 shlint)


!     g1kx: first geo derivatives of kinetic energy (x) integrals
!     ===========================================================

         else if (inttyp == 103) then
           if (icenta /= icentb) then
             call g1k_integrals(1, 0, 0,
     &                          jmaxa,
     &                          jmaxb,
     &                          jmaxt,
     &                          jmaxd,
     &                          jmaxm,
     &                          odc,
     &                          shgtf,
     &                          shlint)
           end if


!     g1ky: first geo derivatives of kinetic energy (y) integrals
!     ===========================================================

         else if (inttyp == 104) then
           if (icenta /= icentb) then
             call g1k_integrals(0, 1, 0,
     &                          jmaxa,
     &                          jmaxb,
     &                          jmaxt,
     &                          jmaxd,
     &                          jmaxm,
     &                          odc,
     &                          shgtf,
     &                          shlint)
           end if


!     g1kz: first geo derivatives of kinetic energy (z) integrals
!     ===========================================================

         else if (inttyp == 105) then
           if (icenta /= icentb) then
             call g1k_integrals(0, 0, 1,
     &                          jmaxa,
     &                          jmaxb,
     &                          jmaxt,
     &                          jmaxd,
     &                          jmaxm,
     &                          odc,
     &                          shgtf,
     &                          shlint)
           end if

!     CEXPIKR: expikr integrals for interaction with electromagnetic field
!     radiation
!     ===========================================================

         ELSE IF (INTTYP == 106) THEN
           CALL CEXPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &          CORPX,CORPY,CORPZ,EXPPI,NOPTYP,WORK,LWORK)
C
C        rp integrals, RM1H3  i c Q_MN < r (alpha.p) >
C        --------------------------------------------------
CMI      Thomas Enevoldsen (for 4-comp London orbitals)
C
         ELSE IF (INTTYP .EQ. 201) THEN
           CALL RM1H3(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                DIFABX,DIFABY,DIFABZ)

C
C     Dipole length integrals around the nuclei to which the atomic
C     orbitals are attached, RM1RN ... < r_n x alpha >
C     -------------------------------------------------------------
CMI        Thomas Enevoldsen (for 4-comp London orbitals)
C
         ELSE IF (INTTYP .EQ. 202) THEN
           CALL DPLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                 CORPX,CORPY,CORPZ,TEMP)
C
C        Dipole_x_Potential energy integrals, RM1H2
C        ... London orbital contribution to the potential V
C        --------------------------------------------------
CMI        Thomas Enevoldsen (for 4-comp London orbitals)
C
         ELSE IF (INTTYP .EQ. 203) THEN
           CALL MM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &         AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,.FALSE.)
C
C        relativ. London contrib. to the expect. value of NMR shielding, RNST (NSLO)
C      ---------------------------------------------------------------------------------
CMI           integrals i/c Q_MN <r (r_K x alp/r_K^3) >
         ELSE IF (INTTYP .EQ. 204) THEN
           CALL RNSLOIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                SHLINT,NOPTYP,NATOMC,DIFABX,DIFABY,DIFABZ)
C
C      London Relativistic Orbital contribution to Diamagnetic
C    Susceptibility tensor  'RDSUSNL' (analogous to nonrel. 'DSUSNL')
C  --------------------------------------------------------------------------
CMI      - use the same routine as INTTYP 22 (DSUSNL)

         ELSE IF (INTTYP .EQ. 205) THEN
           CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C       London Relativistic Orbital contribution to Diamagnetic
C       susceptibilty tensor  'RDSUSLL' (analogous to nonrel. 'DSUSLL')
C     -------------------------------------------------------------------
C      integrals i c Qmn <r (rn x alpha)>
C      HJAaJ & MI

         ELSE IF (INTTYP .EQ. 206) THEN
          IF (DISTAB .GT. ZERODIST)
     &     CALL RDSUSAN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &             SHLINT,CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ)
C
C       2nd derivative of relativ.kinetic term, RM2H3
C      --------------------------------------------------------------
C         HJAaJ & MI
         ELSE IF (INTTYP .EQ. 207) THEN
           CALL RM2H3(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                 DIFABX,DIFABY,DIFABZ)
C
C        QrrQ_x_Potential energy integrals, RM2H2
C        ... London orbital contribution to the potential V
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 208) THEN
           CALL DSUSLH(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                 AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,.FALSE.)



C        Second order magnetic derivatives of overlap matrix
C        ---------------------------------------------------
Cdj bra part
         ELSE IF (INTTYP .EQ. 304) THEN
            DIFAOX = CORAX - GAGORG(1)
            DIFAOY = CORAY - GAGORG(2)
            DIFAOZ = CORAZ - GAGORG(3)
            CALL SM2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFAOX,DIFAOY,DIFAOZ)
Cdj ket part
         ELSE IF (INTTYP .EQ. 305) THEN
            DIFOBX = GAGORG(1) - CORBX
            DIFOBY = GAGORG(2) - CORBY
            DIFOBZ = GAGORG(3) - CORBZ
            CALL SM2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFOBX,DIFOBY,DIFOBZ)
Cdj mixed bra-ket part
         ELSE IF (INTTYP .EQ. 306) THEN
            DIFAOX = CORAX - GAGORG(1)
            DIFAOY = CORAY - GAGORG(2)
            DIFAOZ = CORAZ - GAGORG(3)
            DIFOBX = GAGORG(1) - CORBX
            DIFOBY = GAGORG(2) - CORBY
            DIFOBZ = GAGORG(3) - CORBZ
            CALL SM2BIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFAOX,DIFAOY,DIFAOZ,DIFOBX,DIFOBY,DIFOBZ)
         END IF
  200    CONTINUE
  100 CONTINUE
C
C     Print
C
      IF (IPRINT .GT. 30) THEN
         DO 400 I = 1, NOPTYP
            CALL AROUND('SHLINT for '//LABINT(I)//' in PR1DR1')
            CALL OUTPUT(SHLINT(1,I),1,KCKTA,1,KCKTB,KCKTA,KCKTB,1,LUPRI)
  400    CONTINUE
      END IF
      RETURN
      END
C  /* Deck pr1prd */
      SUBROUTINE PR1PRD(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  WORK,LWORK,NPOINT,LABINT,INTTYP,NOPTYP,NBAST,
     &                  DOATOM,WEIGHT,ABSCIS,TRIANG,IPRINT,TOLS,TOLOG,
     &                  DISTAB,FACINT,CORCX,CORCY,CORCZ,DIFABX,
     &                  DIFABY,DIFABZ,KAB,DONUC1,DOMOM1)
!
!     Special subroutine for DSO - diamagnetic spin-orbit integrals
!
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (D1 = 1.00 D00, D3 = 3.00 D00, D3INV = D1/D3)
C
      LOGICAL DOATOM(NUCIND), TRIANG, SAMECD, DONUC1, DOMOM1
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION WORK(LWORK), WEIGHT(NPOINT), ABSCIS(NPOINT),
     &          FACINT(*), CORCX(*), CORCY(*),
     &          CORCZ(*), SHLINT(KCKTAB,*), AHGTF(*),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,3)
C
#include "onecom.h"
#include "orgcom.h"
#include "nuclei.h"
#include "primit.h"
#include "symmet.h"
C
#include "ibtfun.h"
      XYZ(I,J) = PT(IBTAND(ISYMAX(I,1),J))
C
      DO 100 IPRIMA = 1,NUCA
         JPRIMA = JSTA + IPRIMA
         CONTA  = PRICCF(JPRIMA,NUMCFA)
         EXPA   = PRIEXP(JPRIMA)
      DO 100 IPRIMB = 1,NUCB
         JPRIMB = JSTB + IPRIMB
         CONTB  = PRICCF(JPRIMB,NUMCFB)
         EXPB   = PRIEXP(JPRIMB)
C
C        Loop over fourth center
C        -----------------------
C
         ISTDSO = 1
         DO 200 IATOMD = 1, NUCIND
         IF (DOATOM(IATOMD)) THEN
            MULD = ISTBNU(IATOMD)
            MABD = IBTOR(MULD,KAB)
            KABD = IBTAND(MULD,KAB)
            DO 300 ISYMD = 0, MAXOPR
            IF (IBTAND(ISYMD,MABD).EQ.0) THEN
               CORDX = XYZ(1,ISYMD)*CORD(1,IATOMD)
               CORDY = XYZ(2,ISYMD)*CORD(2,IATOMD)
               CORDZ = XYZ(3,ISYMD)*CORD(3,IATOMD)
               DISTAD = (CORAX - CORDX)**2 + (CORAY - CORDY)**2
     &                                     + (CORAZ - CORDZ)**2
               DISTBD = (CORBX - CORDX)**2 + (CORBY - CORDY)**2
     &                                     + (CORBZ - CORDZ)**2
C
C              Prepare for loop over atom C
C              ----------------------------
C
               NATOMC = 0
               MXATMC = NUCIND
               IF (TRIANG) MXATMC = IATOMD
               DO 400 IATOMC = 1, MXATMC
               IF (DOATOM(IATOMC)) THEN
                  SAMECD = IATOMC .EQ. IATOMD
                  MULC   = ISTBNU(IATOMC)
                  MABCD  = IBTOR(MULC,KABD)
                  FACTOR = - FMULT(IBTAND(MULC,KABD))
                  DO 500 ISYMC = 0, MAXOPR
                  IF (SAMECD .AND. ISYMC.EQ.ISYMD) GO TO 500
                  IF (IBTAND(ISYMC,MABCD) .EQ. 0) THEN
                     NATOMC = NATOMC + 1
                     CORCX(NATOMC) = XYZ(1,ISYMC)*CORD(1,IATOMC)
                     CORCY(NATOMC) = XYZ(2,ISYMC)*CORD(2,IATOMC)
                     CORCZ(NATOMC) = XYZ(3,ISYMC)*CORD(3,IATOMC)
                     FACINT(NATOMC)= FACTOR
                  END IF
  500             CONTINUE
               END IF
  400          CONTINUE
               IF (NATOMC .EQ. 0) GO TO 300
C
C              Loop over quadrature points
C              ---------------------------
C
               DO 600 IPOINT = 1, NPOINT
                  EXPD   = ABSCIS(IPOINT)
                  EXPP   = EXPA + EXPB + EXPD
                  EXPPI  = D1/EXPP
                  EXPABQ = EXPPI*(EXPA*EXPB*DISTAB
     &                   + EXPA*EXPD*DISTAD + EXPB*EXPD*DISTBD)
               IF (EXPABQ.GT.TOLOG) GO TO 600
                  SAAB   = WEIGHT(IPOINT)*CONTA*CONTB
     &                                   *EXP(-EXPABQ)
                  ASAAB  = ABS(SAAB)
               IF (ASAAB.LT.TOLS) GO TO 600
                  SAAB13 = SIGN(ASAAB**D3INV,SAAB)
                  EXPAPI = EXPA*EXPPI
                  EXPBPI = EXPB*EXPPI
                  EXPCPI = EXPD*EXPPI
                  CORPX =EXPAPI*CORAX+EXPBPI*CORBX+EXPCPI*CORDX
                  CORPY =EXPAPI*CORAY+EXPBPI*CORBY+EXPCPI*CORDY
                  CORPZ =EXPAPI*CORAZ+EXPBPI*CORBZ+EXPCPI*CORDZ
C
C                 *********************************************
C                 ***** Overlap Distribution Coefficients *****
C                 *********************************************
C
                  CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,.FALSE.,
     &                        .FALSE.,.FALSE.,EXPA,EXPB,IPRINT,SAAB13,
     &                        EXPPI,WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,
     &                        DOMOM1,ORIGIN,INTTYP)
C
                  DIFPDX = CORPX - CORDX
                  DIFPDY = CORPY - CORDY
                  DIFPDZ = CORPZ - CORDZ
                  CALL DSOODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,EXPPI,
     &                        DIFPDX,DIFPDY,DIFPDZ,IPRINT)
C
C                 ********************************************
C                 **** Calculation of Hermitian integrals ****
C                 ********************************************
C
                  DO 700 I = 1, NATOMC
                     FACTOR = FACINT(I)
                     DIFCPX = CORCX(I) - CORPX
                     DIFCPY = CORCY(I) - CORPY
                     DIFCPZ = CORCZ(I) - CORPZ
                     IADR = 1 + (I - 1)*NAHGTF
                     CALL HERNAI(AHGTF,JMAX,EXPP,DIFCPX,DIFCPY,DIFCPZ,
     &                           FACTOR,IADR,ISTEPU,ISTEPV,NAHGTF,
     &                           IPRINT)
  700             CONTINUE
C
C                 ********************************************
C                 **** Calculation of Cartesian integrals ****
C                 ********************************************
C
                  CALL DSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                        SHLINT(1,ISTDSO),NOPTYP,NATOMC)
  600          CONTINUE
               ISTDSO = ISTDSO + 9*NATOMC
            END IF
  300       CONTINUE
         END IF
  200    CONTINUE
C
C     End loop over primitives
C
  100 CONTINUE
C
C     Print
C
      IF (IPRINT .GT. 30) THEN
         DO 900 I = 1, NOPTYP
            CALL AROUND('SHLINT for '//LABINT(I)//' in PR1DR1')
            CALL OUTPUT(SHLINT(1,I),1,KCKTA,1,KCKTB,KCKTA,KCKTB,1,LUPRI)
  900    CONTINUE
      END IF
      RETURN
      END
C  /* Deck ssoint */
      SUBROUTINE SSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT)
C
C     Spatial one-electron spin-orbit integrals
C
C     tuh 23 Nov 1989
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          AHGTF(*), SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               EE = EU*EV
               FE = FU*EV
               EF = EU*FV
               DO 400 IT = 0,MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  FEE = FT*EE
                  EFE = ET*FE
                  EEF = ET*EF
                  AH0T = AHGTF(IADRAU + IT + 1)
                  AH0U = AHGTF(IADRAU + IT + ISTEPU)
                  AH0V = AHGTF(IADRAU + IT + ISTEPV)
                  SHLINT(INT,1) = SHLINT(INT,1) + EFE*AH0V - EEF*AH0U
                  SHLINT(INT,2) = SHLINT(INT,2) + EEF*AH0T - FEE*AH0V
                  SHLINT(INT,3) = SHLINT(INT,3) + FEE*AH0U - EFE*AH0T
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
         INT = INT + 1
  100 CONTINUE
      RETURN
      END
C  /* Deck ovlint */
      SUBROUTINE OVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,HUCFAC,
     &                  SHGTF,SHLINT)
C
C     tuh 1989
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(*)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
C
         SHLINT(INT) = SHLINT(INT) - SX0*SY0*SZ0*HUCFAC
  100 CONTINUE
      RETURN
      END

      subroutine g1n_integrals(jmaxa,
     &                         jmaxb,
     &                         jmaxt,
     &                         jmaxd,
     &                         jmaxm,
     &                         odc,
     &                         natomc,
     &                         ncentc,
     &                         ahgtf,
     &                         noptyp,
     &                         shlint)

      implicit none

#include "mxcent.h"
#include "maxaqn.h"
#include "onecom.h"
#include "lmns.h"

!     --------------------------------------------------------------------------
      integer, intent(in)    :: jmaxa
      integer, intent(in)    :: jmaxb
      integer, intent(in)    :: jmaxt
      integer, intent(in)    :: jmaxd
      integer, intent(in)    :: jmaxm
      real(8) :: odc(0:jmaxa, 0:jmaxb, 0:jmaxt, 0:jmaxd, 0:jmaxm, 3)
      integer, intent(in)    :: natomc
      integer, intent(in)    :: ncentc(*)
      real(8), intent(in)    :: ahgtf(*)
      integer, intent(in)    :: noptyp
      real(8), intent(inout) :: shlint(kcktab, noptyp, 3)
!     --------------------------------------------------------------------------
      integer :: i, ioff, iatom
      integer :: ix,  iy,  iz
      integer :: la,  ma,  na
      integer :: lb,  mb,  nb
      integer :: ia,  ib
      integer :: iax, iay, iaz, ia0
      integer :: iu,  iv
      integer :: xa,  ya,  za
      integer :: xb,  yb,  zb
      integer :: xc,  yc,  zc
      real(8) :: ax,  ay,  az
      real(8) :: cx,  cy,  cz
      real(8) :: fx,  fy,  fz
      real(8) :: gx,  gy,  gz
      real(8) :: ahx, ahy, ahz, ah0
      real(8) :: ff, gf, fg, gg, fff, gff, fgf, ffg
      logical :: ab_same_center
!     --------------------------------------------------------------------------

      ab_same_center = (icenta == icentb)

      xa = 3*(icenta - 1) + 1
      ya = 3*(icenta - 1) + 2
      za = 3*(icenta - 1) + 3

      xb = 3*(icentb - 1) + 1
      yb = 3*(icentb - 1) + 2
      zb = 3*(icentb - 1) + 3

      i = 0
      do ia = 1, kckta

        la = lvalua(ia)
        ma = mvalua(ia)
        na = nvalua(ia)

        do ib = 1, kcktb
          i = i + 1

          lb = lvalub(ib)
          mb = mvalub(ib)
          nb = nvalub(ib)

          iv = 1
          do iz = 0, na + nb + 1
            iu = iv

            fz = odc(na, nb, iz, 0, 0, 3)
            if (.not. ab_same_center) then
              gz = odc(na, nb, iz, 1, 0, 3)
            end if

            do iy = 0, ma + mb + 1

              fy = odc(ma, mb, iy, 0, 0, 2)
              ff = fy*fz
              if (.not. ab_same_center) then
                gy = odc(ma, mb, iy, 1, 0, 2)
                gf = gy*fz
                fg = fy*gz
                gg = gy*gz
              end if

              do ix = 0, la + lb + 1

                fx = odc(la, lb, ix, 0, 0, 1)
                fff = fx*ff
                if (.not. ab_same_center) then
                  gx = odc(la, lb, ix, 1, 0, 1)
                  gff = gx*ff
                  fgf = fx*gf
                  ffg = fx*fg
                end if

                ia0 = iu  + ix
                iax = ia0 + 1
                iay = ia0 + istepu
                iaz = ia0 + istepv

                ioff = 0
                do iatom = 1, natomc
 
                  if (.not. ab_same_center) then

                    ah0 = ahgtf(ia0 + ioff)

                    ax =  gff*ah0
                    ay =  fgf*ah0
                    az =  ffg*ah0

                    shlint(i, xa, 1) = shlint(i, xa, 1) + ax
                    shlint(i, ya, 1) = shlint(i, ya, 1) + ay
                    shlint(i, za, 1) = shlint(i, za, 1) + az
                    shlint(i, xb, 2) = shlint(i, xb, 2) - ax
                    shlint(i, yb, 2) = shlint(i, yb, 2) - ay
                    shlint(i, zb, 2) = shlint(i, zb, 2) - az

                  end if

                  ahx = ahgtf(iax + ioff)
                  ahy = ahgtf(iay + ioff)
                  ahz = ahgtf(iaz + ioff)

                  cx = -fff*ahx
                  cy = -fff*ahy
                  cz = -fff*ahz

                  xc = 3*(ncentc(iatom) - 1) + 1
                  yc = 3*(ncentc(iatom) - 1) + 2
                  zc = 3*(ncentc(iatom) - 1) + 3
 
                  shlint(i, xb, 2) = shlint(i, xb, 2) - cx
                  shlint(i, yb, 2) = shlint(i, yb, 2) - cy
                  shlint(i, zb, 2) = shlint(i, zb, 2) - cz
                  shlint(i, xc, 3) = shlint(i, xc, 3) + cx
                  shlint(i, yc, 3) = shlint(i, yc, 3) + cy
                  shlint(i, zc, 3) = shlint(i, zc, 3) + cz
      
                  ioff = ioff + nahgtf
                end do
              end do
              iu = iu + istepu
            end do
            iv = iv + istepv
          end do
        end do
      end do

      end subroutine


C  /* Deck dplint */
      SUBROUTINE DPLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,ORIGIN)
C
C     tuh 1989
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3), ORIGIN(3)
#include "onecom.h"
#include "lmns.h"
      POX = CORPX - ORIGIN(1)
      POY = CORPY - ORIGIN(2)
      POZ = CORPZ - ORIGIN(3)
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = POX*SX0
         DY0 = POY*SY0
         DZ0 = POZ*SZ0
         IF (LVALA+LVALB.GT.0) DX0 = DX0+ SHGTF*ODC(LVALA,LVALB,1,0,0,1)
         IF (MVALA+MVALB.GT.0) DY0 = DY0+ SHGTF*ODC(MVALA,MVALB,1,0,0,2)
         IF (NVALA+NVALB.GT.0) DZ0 = DZ0+ SHGTF*ODC(NVALA,NVALB,1,0,0,3)
C
         INT = INT + 1
         SHLINT(INT,1) = SHLINT(INT,1) - DX0*SY0*SZ0
         SHLINT(INT,2) = SHLINT(INT,2) - SX0*DY0*SZ0
         SHLINT(INT,3) = SHLINT(INT,3) - SX0*SY0*DZ0
  100 CONTINUE
      RETURN
      END
C  /* Deck kinint */
      SUBROUTINE KININT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, Nov. 1991. Based on TUH's subroutine CINT0, but modified
C     for use in HERPRO.
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (DP5 = 0.5 D00)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C    **********************************************
C    ***** CALCULATE KINETIC ENERGY INTEGRALS *****
C    **********************************************
C
         X0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         Y0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         Z0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         X2 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         Y2 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         Z2 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         INT = INT + 1
         SHLINT(INT) = SHLINT(INT)
     &               + DP5*(X2*Y0*Z0 + X0*Y2*Z0 + X0*Y0*Z2)
 100  CONTINUE
      RETURN
      END
C  /* Deck mvlint */
      SUBROUTINE MVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     Sheela Kirpekar. jan. 93 (based on KININT)
C
C     K.Ruud, Nov. 1991. Based on TUH's subroutine CINT0, but modified
C     for use in HERPRO.
C
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB)
#include "onecom.h"
#include "lmns.h"
      DOUBLE PRECISION ALPMVL,ALPMV2

      ALPMVL = ALPHA2*0.125D0
      ALPMV2 = ALPHA2*0.25D0

      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C    **********************************************
C    ***** CALCULATE MASS VELOCITY INTEGRALS ******
C    **********************************************
C
         X0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         Y0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         Z0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         X2 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         Y2 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         Z2 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         X4 = SHGTF*ODC(LVALA,LVALB,0,4,0,1)
         Y4 = SHGTF*ODC(MVALA,MVALB,0,4,0,2)
         Z4 = SHGTF*ODC(NVALA,NVALB,0,4,0,3)
         INT = INT + 1
         SHLINT(INT) = SHLINT(INT)
     &               + ALPMVL*(X4*Y0*Z0 + X0*Y4*Z0 + X0*Y0*Z4)
     &               + ALPMV2*(X2*Y2*Z0 + X2*Y0*Z2 + X0*Y2*Z2)
 100  CONTINUE
      RETURN
      END
C  /* Deck hdbint */
      SUBROUTINE HDBINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFOBX,DIFOBY,DIFOBZ)
C
C     K.Ruud, Aug. 1992
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (DP5 = 0.5D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,9)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         INT = INT + 1
C
         SHLINT(INT,1) = SHLINT(INT,1) - DP5*(DIFOBY*SY0*SZ1
     &                                 -      DIFOBZ*SY1*SZ0)*DX0
         SHLINT(INT,2) = SHLINT(INT,2) - DP5*(DIFOBZ*DX1*SZ0
     &                                 -      DIFOBX*DX0*SZ1)*SY0
         SHLINT(INT,3) = SHLINT(INT,3) - DP5*(DIFOBX*DX0*SY1
     &                                 -      DIFOBY*DX1*SY0)*SZ0
         SHLINT(INT,4) = SHLINT(INT,4) - DP5*(DIFOBY*DY0*SZ1
     &                                 -      DIFOBZ*DY1*SZ0)*SX0
         SHLINT(INT,5) = SHLINT(INT,5) - DP5*(DIFOBZ*SX1*SZ0
     &                                 -      DIFOBX*SX0*SZ1)*DY0
         SHLINT(INT,6) = SHLINT(INT,6) - DP5*(DIFOBX*SX0*DY1
     &                                 -      DIFOBY*SX1*DY0)*SZ0
         SHLINT(INT,7) = SHLINT(INT,7) - DP5*(DIFOBY*SY0*DZ1
     &                                 -      DIFOBZ*SY1*DZ0)*SX0
         SHLINT(INT,8) = SHLINT(INT,8) - DP5*(DIFOBZ*SX1*DZ0
     &                                 -      DIFOBX*SX0*DZ1)*SY0
         SHLINT(INT,9) = SHLINT(INT,9) - DP5*(DIFOBX*SX0*SY1
     &                                 -      DIFOBY*SX1*SY0)*DZ0
 100  CONTINUE
      RETURN
      END
C  /* Deck sm1int */
      SUBROUTINE SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C     Kenneth Ruud's first subroutine. Modified Dec. 1991, KR
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D2I = 0.5D00)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
C
         DPLX = DX0*SY0*SZ0
         DPLY = SX0*DY0*SZ0
         DPLZ = SX0*SY0*DZ0
C
         SHLINT(INT,1) = SHLINT(INT,1) - D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
         SHLINT(INT,2) = SHLINT(INT,2) - D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
         SHLINT(INT,3) = SHLINT(INT,3) - D2I*(DIFABX*DPLY - DIFABY*DPLX)
 100  CONTINUE
      RETURN
      END
C  /* Deck sm2int */
      SUBROUTINE SM2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C     K. Ruud, Oct. 1991
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4INV = 0.25D0)
      DIMENSION SHLINT(KCKTAB,6)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - (D2*DIFABZ*DIFABY*DY1*DZ1*SX0 -
     &                   DIFABZ*DIFABZ*DY2*SX0*SZ0 -
     &                   DIFABY*DIFABY*DZ2*SX0*SY0)*D4INV
         SHLINT(INT,2) = SHLINT(INT,2) - (DIFABZ*DIFABZ*DX1*DY1*SZ0 -
     &                   DIFABY*DIFABZ*DZ1*DX1*SY0 -
     &                   DIFABZ*DIFABX*DY1*DZ1*SX0 +
     &                   DIFABY*DIFABX*DZ2*SX0*SY0)*D4INV
         SHLINT(INT,3) = SHLINT(INT,3) - (DIFABY*DIFABY*DX1*DZ1*SY0 -
     &                   DIFABX*DIFABY*DY1*DZ1*SX0 +
     &                   DIFABZ*DIFABX*DY2*SX0*SZ0 -
     &                   DIFABY*DIFABZ*DY1*DX1*SZ0)*D4INV
         SHLINT(INT,4) = SHLINT(INT,4) - (D2*DIFABZ*DIFABX*DX1*DZ1*SY0 -
     &                   DIFABZ*DIFABZ*DX2*SY0*SZ0 -
     &                   DIFABX*DIFABX*DZ2*SX0*SY0)*D4INV
         SHLINT(INT,5) = SHLINT(INT,5) - (DIFABZ*DIFABY*DX2*SY0*SZ0 -
     &                   DIFABZ*DIFABX*DX1*DY1*SZ0 -
     &                   DIFABY*DIFABX*DX1*DZ1*SY0 +
     &                   DIFABX*DIFABX*DY1*DZ1*SX0)*D4INV
         SHLINT(INT,6) = SHLINT(INT,6) - (D2*DIFABX*DIFABY*DX1*DY1*SZ0 -
     &                   DIFABY*DIFABY*DX2*SY0*SZ0 -
     &                   DIFABX*DIFABX*DY2*SX0*SZ0)*D4INV
 100  CONTINUE
      RETURN
      END

      SUBROUTINE SM2BIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFAX,DIFAY,DIFAZ,DIFBX,DIFBY,DIFBZ)
C
C     dj oct 08, based on SM2INT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4INV = 0.25D0)
      DIMENSION SHLINT(KCKTAB,9)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - (DIFAZ*DIFBY*SX0*DY1*DZ1 +
     &                   DIFAY*DIFBZ*SX0*DY1*DZ1 -
     &                   DIFAZ*DIFBZ*SX0*DY2*SZ0 -
     &                   DIFAY*DIFBY*SX0*SY0*DZ2)*D4INV
         SHLINT(INT,2) = SHLINT(INT,2) - (DIFAZ*DIFBZ*DX1*DY1*SZ0 -
     &                   DIFAY*DIFBZ*DX1*SY0*DZ1 -
     &                   DIFAZ*DIFBX*SX0*DY1*DZ1 +
     &                   DIFAY*DIFBX*SX0*SY0*DZ2)*D4INV
         SHLINT(INT,3) = SHLINT(INT,3) - (DIFAY*DIFBY*DX1*SY0*DZ1 -
     &                   DIFAY*DIFBX*SX0*DY1*DZ1 +
     &                   DIFAZ*DIFBX*DY2*SX0*SZ0 -
     &                   DIFAZ*DIFBY*DX1*DY1*SZ0)*D4INV
         SHLINT(INT,4) = SHLINT(INT,4) - (DIFAZ*DIFBZ*DX1*DY1*SZ0 -
     &                   DIFAZ*DIFBY*DX1*SY0*DZ1 -
     &                   DIFAX*DIFBZ*SX0*DY1*DZ1 +
     &                   DIFAX*DIFBY*SX0*SY0*DZ2)*D4INV
         SHLINT(INT,5) = SHLINT(INT,5) - (DIFAZ*DIFBX*DX1*SY0*DZ1 +
     &                   DIFAX*DIFBZ*DX1*SY0*DZ1 -
     &                   DIFAZ*DIFBZ*DX2*SY0*SZ0 -
     &                   DIFAX*DIFBX*SX0*SY0*DZ2)*D4INV
         SHLINT(INT,6) = SHLINT(INT,6) - (DIFAZ*DIFBY*DX2*SY0*SZ0 -
     &                   DIFAZ*DIFBX*DX1*DY1*SZ0 -
     &                   DIFAX*DIFBY*DX1*DZ1*SY0 +
     &                   DIFAX*DIFBX*DY1*DZ1*SX0)*D4INV
         SHLINT(INT,7) = SHLINT(INT,7) - (DIFAY*DIFBY*DX1*SY0*DZ1 -
     &                   DIFAX*DIFBY*SX0*DY1*DZ1 +
     &                   DIFAX*DIFBZ*SX0*DY2*SZ0 -
     &                   DIFAY*DIFBZ*DX1*DY1*SZ0)*D4INV
         SHLINT(INT,8) = SHLINT(INT,8) - (DIFAY*DIFBZ*DX2*SY0*SZ0 -
     &                   DIFAX*DIFBZ*DX1*DY1*SZ0 -
     &                   DIFAY*DIFBX*DX1*SY0*DZ1 +
     &                   DIFAX*DIFBX*SX0*DY1*DZ1)*D4INV
         SHLINT(INT,9) = SHLINT(INT,9) - (DIFAX*DIFBY*DX1*DY1*SZ0 +
     &                   DIFAY*DIFBX*DX1*DY1*SZ0 -
     &                   DIFAY*DIFBY*DX2*SY0*SZ0 -
     &                   DIFAX*DIFBX*SX0*DY2*SZ0)*D4INV
 100  CONTINUE
      RETURN
      END



C  /* Deck cm1int */
      SUBROUTINE CM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C     K. Ruud, Aug.-93.
C     Rewritten in order to separate dipole origin, aug.-94 KR
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "efield.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D2I = 1.0D0/D2)
      DIMENSION SHLINT(KCKTAB,9)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include "onecom.h"
#include "orgcom.h"
#include "lmns.h"
      INT = 0
      ADX = CORAX - DIPORG(1)
      ADY = CORAY - DIPORG(2)
      ADZ = CORAZ - DIPORG(3)
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
C
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
         DX1 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,0,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,0,1))
         DY1 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,0,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,0,2))
         DZ1 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,0,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,0,3))
         DX2 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,1,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,1,1))
         DY2 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,1,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,1,2))
         DZ2 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,1,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,1,3))

         SHLINT(INT,1) = SHLINT(INT,1) 
     &                 - D2I*(DIFABY*DX1*SY0*SZ1 - DIFABZ*DX1*SY1*SZ0)
         SHLINT(INT,2) = SHLINT(INT,2) 
     &                 - D2I*(DIFABZ*DX2*SY0*SZ0 - DIFABX*DX1*SY0*SZ1)
         SHLINT(INT,3) = SHLINT(INT,3) 
     &                 - D2I*(DIFABX*DX1*SY1*SZ0 - DIFABY*DX2*SY0*SZ0)

         SHLINT(INT,4) = SHLINT(INT,4) 
     &                 - D2I*(DIFABY*SX0*DY1*SZ1 - DIFABZ*SX0*DY2*SZ0)
         SHLINT(INT,5) = SHLINT(INT,5) 
     &                 - D2I*(DIFABZ*SX1*DY1*SZ0 - DIFABX*SX0*DY1*SZ1)
         SHLINT(INT,6) = SHLINT(INT,6) 
     &                 - D2I*(DIFABX*SX0*DY2*SZ0 - DIFABY*SX1*DY1*SZ0)

         SHLINT(INT,7) = SHLINT(INT,7) 
     &                 - D2I*(DIFABY*SX0*SY0*DZ2 - DIFABZ*SX0*SY1*DZ1)
         SHLINT(INT,8) = SHLINT(INT,8) 
     &                 - D2I*(DIFABZ*SX1*SY0*DZ1 - DIFABX*SX0*SY0*DZ2)
         SHLINT(INT,9) = SHLINT(INT,9) 
     &                 - D2I*(DIFABX*SX0*SY1*DZ1 - DIFABY*SX1*SY0*DZ1)
 100  CONTINUE
      RETURN
      END
C  /* Deck cm2int */
      SUBROUTINE CM2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C     K. Ruud, Aug.-93
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "efield.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4INV = 0.25D0)
      DIMENSION SHLINT(KCKTAB,18)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include "onecom.h"
#include "orgcom.h"
#include "lmns.h"
      INT = 0
      ADX = CORAX - DIPORG(1)
      ADY = CORAY - DIPORG(2)
      ADZ = CORAZ - DIPORG(3)
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
C
         DX1 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,0,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,0,1))
         DY1 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,0,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,0,2))
         DZ1 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,0,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,0,3))
         DX2 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,1,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,1,1))
         DY2 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,1,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,1,2))
         DZ2 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,1,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,1,3))
         DX3 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,2,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,2,1))
         DY3 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,2,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,2,2))
         DZ3 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,2,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,2,3))
C
         DPLXX = DX3*SY0*SZ0
         DPLXY = DX2*SY1*SZ0
         DPLXZ = DX2*SY0*SZ1
         DPLYY = DX1*SY2*SZ0
         DPLYZ = DX1*SY1*SZ1
         DPLZZ = DX1*SY0*SZ2

         ioff = 0

         SHLINT(INT,ioff+1) = SHLINT(INT,ioff+1) 
     &                - (D2*DIFABZ*DIFABY*DPLYZ -
     &                   DIFABZ*DIFABZ*DPLYY -
     &                   DIFABY*DIFABY*DPLZZ)*D4INV
         SHLINT(INT,ioff+2) = SHLINT(INT,ioff+2) 
     &                - (DIFABZ*DIFABZ*DPLXY -
     &                   DIFABY*DIFABZ*DPLXZ -
     &                   DIFABZ*DIFABX*DPLYZ +
     &                   DIFABY*DIFABX*DPLZZ)*D4INV
         SHLINT(INT,ioff+3) = SHLINT(INT,ioff+3) 
     &                - (DIFABY*DIFABY*DPLXZ -
     &                   DIFABX*DIFABY*DPLYZ +
     &                   DIFABZ*DIFABX*DPLYY -
     &                   DIFABY*DIFABZ*DPLXY)*D4INV
         SHLINT(INT,ioff+4) = SHLINT(INT,ioff+4) 
     &                - (D2*DIFABZ*DIFABX*DPLXZ -
     &                   DIFABZ*DIFABZ*DPLXX -
     &                   DIFABX*DIFABX*DPLZZ)*D4INV
         SHLINT(INT,ioff+5) = SHLINT(INT,ioff+5) 
     &                - (DIFABZ*DIFABY*DPLXX -
     &                   DIFABZ*DIFABX*DPLXY -
     &                   DIFABY*DIFABX*DPLXZ +
     &                   DIFABX*DIFABX*DPLYZ)*D4INV
         SHLINT(INT,ioff+6) = SHLINT(INT,ioff+6) 
     &                - (D2*DIFABX*DIFABY*DPLXY -
     &                   DIFABY*DIFABY*DPLXX -
     &                   DIFABX*DIFABX*DPLYY)*D4INV

         DPLXX = SX2*DY1*SZ0
         DPLXY = SX1*DY2*SZ0
         DPLXZ = SX1*DY1*SZ1
         DPLYY = SX0*DY3*SZ0
         DPLYZ = SX0*DY2*SZ1
         DPLZZ = SX0*DY1*SZ2

         ioff = 6

         SHLINT(INT,ioff+1) = SHLINT(INT,ioff+1) 
     &                - (D2*DIFABZ*DIFABY*DPLYZ -
     &                   DIFABZ*DIFABZ*DPLYY -
     &                   DIFABY*DIFABY*DPLZZ)*D4INV
         SHLINT(INT,ioff+2) = SHLINT(INT,ioff+2) 
     &                - (DIFABZ*DIFABZ*DPLXY -
     &                   DIFABY*DIFABZ*DPLXZ -
     &                   DIFABZ*DIFABX*DPLYZ +
     &                   DIFABY*DIFABX*DPLZZ)*D4INV
         SHLINT(INT,ioff+3) = SHLINT(INT,ioff+3) 
     &                - (DIFABY*DIFABY*DPLXZ -
     &                   DIFABX*DIFABY*DPLYZ +
     &                   DIFABZ*DIFABX*DPLYY -
     &                   DIFABY*DIFABZ*DPLXY)*D4INV
         SHLINT(INT,ioff+4) = SHLINT(INT,ioff+4) 
     &                - (D2*DIFABZ*DIFABX*DPLXZ -
     &                   DIFABZ*DIFABZ*DPLXX -
     &                   DIFABX*DIFABX*DPLZZ)*D4INV
         SHLINT(INT,ioff+5) = SHLINT(INT,ioff+5) 
     &                - (DIFABZ*DIFABY*DPLXX -
     &                   DIFABZ*DIFABX*DPLXY -
     &                   DIFABY*DIFABX*DPLXZ +
     &                   DIFABX*DIFABX*DPLYZ)*D4INV
         SHLINT(INT,ioff+6) = SHLINT(INT,ioff+6) 
     &                - (D2*DIFABX*DIFABY*DPLXY -
     &                   DIFABY*DIFABY*DPLXX -
     &                   DIFABX*DIFABX*DPLYY)*D4INV

         DPLXX = SX2*SY0*DZ1
         DPLXY = SX1*SY1*DZ1
         DPLXZ = SX1*SY0*DZ2
         DPLYY = SX0*SY2*DZ1
         DPLYZ = SX0*SY1*DZ2
         DPLZZ = SX0*SY0*DZ3

         ioff = 12

         SHLINT(INT,ioff+1) = SHLINT(INT,ioff+1) 
     &                - (D2*DIFABZ*DIFABY*DPLYZ -
     &                   DIFABZ*DIFABZ*DPLYY -
     &                   DIFABY*DIFABY*DPLZZ)*D4INV
         SHLINT(INT,ioff+2) = SHLINT(INT,ioff+2) 
     &                - (DIFABZ*DIFABZ*DPLXY -
     &                   DIFABY*DIFABZ*DPLXZ -
     &                   DIFABZ*DIFABX*DPLYZ +
     &                   DIFABY*DIFABX*DPLZZ)*D4INV
         SHLINT(INT,ioff+3) = SHLINT(INT,ioff+3) 
     &                - (DIFABY*DIFABY*DPLXZ -
     &                   DIFABX*DIFABY*DPLYZ +
     &                   DIFABZ*DIFABX*DPLYY -
     &                   DIFABY*DIFABZ*DPLXY)*D4INV
         SHLINT(INT,ioff+4) = SHLINT(INT,ioff+4) 
     &                - (D2*DIFABZ*DIFABX*DPLXZ -
     &                   DIFABZ*DIFABZ*DPLXX -
     &                   DIFABX*DIFABX*DPLZZ)*D4INV
         SHLINT(INT,ioff+5) = SHLINT(INT,ioff+5) 
     &                - (DIFABZ*DIFABY*DPLXX -
     &                   DIFABZ*DIFABX*DPLXY -
     &                   DIFABY*DIFABX*DPLXZ +
     &                   DIFABX*DIFABX*DPLYZ)*D4INV
         SHLINT(INT,ioff+6) = SHLINT(INT,ioff+6) 
     &                - (D2*DIFABX*DIFABY*DPLXY -
     &                   DIFABY*DIFABY*DPLXX -
     &                   DIFABX*DIFABX*DPLYY)*D4INV

 100  CONTINUE
      RETURN
      END


      SUBROUTINE QDBINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C     K. Ruud, February.-02.
C     Based on CM1INT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "efield.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D2I = 1.0D0/D2, D1P5 = 1.5D0,
     &           D1 = 1.0D0)
      DIMENSION SHLINT(KCKTAB,18)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
!radovan: dummy just to make it compile
      character(16) :: field3
#include "onecom.h"
#include "orgcom.h"
#include "lmns.h"
      INT = 0
      ADX = CORAX - DIPORG(1)
      ADY = CORAY - DIPORG(2)
      ADZ = CORAZ - DIPORG(3)
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
C
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX1 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,0,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,0,1))
         DY1 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,0,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,0,2))
         DZ1 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,0,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,0,3))
         DX2 = SHGTF*(ODC(LVALA + 2,LVALB,0,0,0,1) +
     &        D2*ADX* ODC(LVALA + 1,LVALB,0,0,0,1) +
     &        ADX*ADX*ODC(LVALA    ,LVALB,0,0,0,1))
         DY2 = SHGTF*(ODC(MVALA + 2,MVALB,0,0,0,2) +
     &        D2*ADY* ODC(MVALA + 1,MVALB,0,0,0,2) +
     &        ADY*ADY*ODC(MVALA    ,MVALB,0,0,0,2))
         DZ2 = SHGTF*(ODC(NVALA + 2,NVALB,0,0,0,3) +
     &        D2*ADZ* ODC(NVALA + 1,NVALB,0,0,0,3) +
     &        ADZ*ADZ*ODC(NVALA    ,NVALB,0,0,0,3))
         SDX2 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,1,1) +
     &            ADX* ODC(LVALA    ,LVALB,0,0,1,1))
         SDY2 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,1,2) +
     &            ADY* ODC(MVALA    ,MVALB,0,0,1,2))
         SDZ2 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,1,3) +
     &            ADZ* ODC(NVALA    ,NVALB,0,0,1,3))
         DX3 = SHGTF*(ODC(LVALA + 2,LVALB,0,0,1,1) +
     &        D2*ADX* ODC(LVALA + 1,LVALB,0,0,1,1) +
     &        ADX*ADX*ODC(LVALA    ,LVALB,0,0,1,1))
         DY3 = SHGTF*(ODC(MVALA + 2,MVALB,0,0,1,2) +
     &        D2*ADY* ODC(MVALA + 1,MVALB,0,0,1,2) +
     &        ADY*ADY*ODC(MVALA    ,MVALB,0,0,1,2))
         DZ3 = SHGTF*(ODC(NVALA + 2,NVALB,0,0,1,3) +
     &        D2*ADZ* ODC(NVALA + 1,NVALB,0,0,1,3) +
     &        ADZ*ADZ*ODC(NVALA    ,NVALB,0,0,1,3))
C
         FAC = 1.5d0

         XR2 = D2I*(DX3*SY0*SZ0 + SX1*DY2*SZ0 + SX1*SY0*DZ2)
         YR2 = D2I*(DX2*SY1*SZ0 + SX0*DY3*SZ0 + SX0*SY1*DZ2)
         ZR2 = D2I*(DX2*SY0*SZ1 + SX0*DY2*SZ1 + SX0*SY0*DZ3)

!        xx
         ioff = 0

         DPLX = - XR2 + FAC*DX3*SY0*SZ0
         DPLY = - YR2 + FAC*DX2*SY1*SZ0
         DPLZ = - ZR2 + FAC*DX2*SY0*SZ1

         SHLINT(INT,ioff+1) = SHLINT(INT,ioff+1) 
     &                 - D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
         SHLINT(INT,ioff+2) = SHLINT(INT,ioff+2) 
     &                 - D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
         SHLINT(INT,ioff+3) = SHLINT(INT,ioff+3) 
     &                 - D2I*(DIFABX*DPLY - DIFABY*DPLX)

!        yy
         ioff = ioff + 3

         DPLX = - XR2 + FAC*SX1*DY2*SZ0
         DPLY = - YR2 + FAC*SX0*DY3*SZ0
         DPLZ = - ZR2 + FAC*SX0*DY2*SZ1

         SHLINT(INT,ioff+1) = SHLINT(INT,ioff+1) 
     &                 - D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
         SHLINT(INT,ioff+2) = SHLINT(INT,ioff+2) 
     &                 - D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
         SHLINT(INT,ioff+3) = SHLINT(INT,ioff+3) 
     &                 - D2I*(DIFABX*DPLY - DIFABY*DPLX)

!        zz
         ioff = ioff + 3

         DPLX = - XR2 + FAC*SX1*SY0*DZ2
         DPLY = - YR2 + FAC*SX0*SY1*DZ2
         DPLZ = - ZR2 + FAC*SX0*SY0*DZ3

         SHLINT(INT,ioff+1) = SHLINT(INT,ioff+1) 
     &                 - D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
         SHLINT(INT,ioff+2) = SHLINT(INT,ioff+2) 
     &                 - D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
         SHLINT(INT,ioff+3) = SHLINT(INT,ioff+3) 
     &                 - D2I*(DIFABX*DPLY - DIFABY*DPLX)

!        xy
         ioff = ioff + 3

         DPLX = FAC*SDX2*DY1*SZ0
         DPLY = FAC*DX1*SDY2*SZ0
         DPLZ = FAC*DX1*DY1*SZ1

         SHLINT(INT,ioff+1) = SHLINT(INT,ioff+1) 
     &                 - D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
         SHLINT(INT,ioff+2) = SHLINT(INT,ioff+2) 
     &                 - D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
         SHLINT(INT,ioff+3) = SHLINT(INT,ioff+3) 
     &                 - D2I*(DIFABX*DPLY - DIFABY*DPLX)

!        xz
         ioff = ioff + 3

         DPLX = FAC*SDX2*SY0*SZ1
         DPLY = FAC*DX1*SY1*DZ1
         DPLZ = FAC*DX1*SY0*SDZ2

         SHLINT(INT,ioff+1) = SHLINT(INT,ioff+1) 
     &                 - D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
         SHLINT(INT,ioff+2) = SHLINT(INT,ioff+2) 
     &                 - D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
         SHLINT(INT,ioff+3) = SHLINT(INT,ioff+3) 
     &                 - D2I*(DIFABX*DPLY - DIFABY*DPLX)

!        yz
         ioff = ioff + 3

         DPLX = FAC*SX1*DY1*DZ1
         DPLY = FAC*SX0*SDY2*DZ1
         DPLZ = FAC*SX0*SY1*SDZ2

         SHLINT(INT,ioff+1) = SHLINT(INT,ioff+1) 
     &                 - D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
         SHLINT(INT,ioff+2) = SHLINT(INT,ioff+2) 
     &                 - D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
         SHLINT(INT,ioff+3) = SHLINT(INT,ioff+3) 
     &                 - D2I*(DIFABX*DPLY - DIFABY*DPLX)

 100  CONTINUE
      RETURN
      END



C  /* Deck am1int */
      SUBROUTINE AM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  FAC)
C
C     K.Ruud, Oct 1991, Modified Dec. 1991, KR
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - (SY1*DZ0 - DY0*SZ1)*SX0*FAC
         SHLINT(INT,2) = SHLINT(INT,2) - (SZ1*DX0 - SX1*DZ0)*SY0*FAC
         SHLINT(INT,3) = SHLINT(INT,3) - (SX1*DY0 - DX0*SY1)*SZ0*FAC
 100  CONTINUE
      RETURN
      END
C  /* Deck mm1int */
      SUBROUTINE MM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,DOKINE)
C
C     K.Ruud, Feb. 1992
C
CMI/HJ March 2003 Added parameter DOKINE to decide whether calculate the only
C       the potential contribution (DOKINE=.false.)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D4INV = 0.250D0, D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3), AHGTF(*)
      LOGICAL DOKINE
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1  = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1  = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1  = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX21 = SHGTF*ODC(LVALA,LVALB,0,2,1,1)
         SY21 = SHGTF*ODC(MVALA,MVALB,0,2,1,2)
         SZ21 = SHGTF*ODC(NVALA,NVALB,0,2,1,3)
         DX2  = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DY2  = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DZ2  = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
C
C       Kinetic energy contribution
CMI/HJ  .... calculate it only when DOKINE=.true.
        IF (DOKINE) THEN

         SHLINT(INT,1) = SHLINT(INT,1) + D4INV*( DIFABY*(DX2*SY0*SZ1 +
     &                                     SX0*DY2*SZ1 + SX0*SY0*SZ21)
     &                                   - DIFABZ*(DX2*SY1*SZ0 +
     &                                     SX0*SY21*SZ0 + SX0*SY1*DZ2))
         SHLINT(INT,2) = SHLINT(INT,2) + D4INV*( DIFABZ*(SX21*SY0*SZ0 +
     &                                     SX1*DY2*SZ0 + SX1*SY0*DZ2)
     &                                   - DIFABX*(DX2*SY0*SZ1 +
     &                                     SX0*DY2*SZ1 + SX0*SY0*SZ21))
         SHLINT(INT,3) = SHLINT(INT,3) + D4INV*( DIFABX*(DX2*SY1*SZ0 +
     &                                     SX0*SY21*SZ0 + SX0*SY1*DZ2)
     &                                   - DIFABY*(SX21*SY0*SZ0 +
     &                                     SX1*DY2*SZ0 + SX1*SY0*DZ2))
        ENDIF              
C
C     Nuclear attraction contribution
CMI/HJ  ... calculate it always

         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV + 1
            IADRAU = IADRAV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            IUMAX = MAXU + 1
            IF (IV .GT. MAXV) IUMAX = MAXU
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               ITMAX = MAXT + 1
               IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) ITMAX = MAXT
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  FX = FT*EU*EV
                  FY = ET*FU*EV
                  FZ = ET*EU*FV
                  ATUV = -D2INV*AHGTF(IADRAU + IT)
                  SHLINT(INT,1)=SHLINT(INT,1)+ATUV*(DIFABY*FZ-DIFABZ*FY)
                  SHLINT(INT,2)=SHLINT(INT,2)+ATUV*(DIFABZ*FX-DIFABX*FZ)
                  SHLINT(INT,3)=SHLINT(INT,3)+ATUV*(DIFABX*FY-DIFABY*FX)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck dsusan */
      SUBROUTINE DSUSAN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT,CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ)
C
C     K.Ruud, Feb. 1992
C
C MI note,Oct.2012: symmetrization of 3x3 tensor implemented inside as stated in 
C the article
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D2INV = 0.50D0, D4INV = 0.25D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6), ORIGIN(3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      PBX = CORBX - CORPX
      PBY = CORBY - CORPY
      PBZ = CORBZ - CORPZ
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1  = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1  = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1  = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SXT0 = SHGTF*ODC(LVALA,LVALB,1,0,0,1)
         SYT0 = SHGTF*ODC(MVALA,MVALB,1,0,0,2)
         SZT0 = SHGTF*ODC(NVALA,NVALB,1,0,0,3)
         SXT1 = SHGTF*ODC(LVALA,LVALB,1,0,1,1)
         SYT1 = SHGTF*ODC(MVALA,MVALB,1,0,1,2)
         SZT1 = SHGTF*ODC(NVALA,NVALB,1,0,1,3)
         SX11 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         SY11 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         SZ11 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         DX0  = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0  = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0  = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
C
         SX1N = SXT0 - PBX*SX0
         SY1N = SYT0 - PBY*SY0
         SZ1N = SZT0 - PBZ*SZ0
         SX2N = SXT1 - PBX*SX1
         SY2N = SYT1 - PBY*SY1
         SZ2N = SZT1 - PBZ*SZ1
C
         SHLINT(INT,1) = SHLINT(INT,1) + (DIFABY*SX0*
     &                                   (SZ11*SY1N - DY0*SZ2N )
     &                                 -  DIFABZ*SX0*
     &                                   (SY2N*DZ0  - SY11*SZ1N))*D2INV
         SHLINT(INT,2) = SHLINT(INT,2) + (DIFABY*SY0*
     &                                   (SZ2N*DX0  - SZ11*SX1N)
     &                                 - (DX0*SY1*SZ1N + SX1*DY0*SZ1N
     &                                 - (SX1*SY1N + SX1N*SY1)*DZ0)*
     &                                    DIFABZ
     &                                 -  DIFABX*SX0*
     &                                   (SZ11*SY1N - SZ2N*DY0))*D4INV
         SHLINT(INT,3) = SHLINT(INT,3) + (DIFABX*SX0*
     &                                   (SY2N*DZ0  - SY11*SZ1N )
     &                                 + (DY0*(SX1N*SZ1 + SX1*SZ1N)
     &                                 -  DX0*SY1N*SZ1 - SX1*SY1N*DZ0)*
     &                                    DIFABY
     &                                 -  DIFABZ*SZ0*
     &                                   (SY11*SX1N - DX0*SY2N))*D4INV
         SHLINT(INT,4) = SHLINT(INT,4) + (DIFABZ*SY0*
     &                                   (SZ1N*SX11  - SX2N*DZ0 )
     &                                 -  DIFABX*SY0*
     &                                   (SZ2N*DX0  - SX1N*SZ11))*D2INV
         SHLINT(INT,5) = SHLINT(INT,5) + (DIFABZ*SZ0*
     &                                   (SX2N*DY0  - SX11*SY1N)
     &                                 + (DX0*(SY1*SZ1N + SY1N*SZ1)
     &                                 -  SX1N*SY1*DZ0 - SX1N*DY0*SZ1)*
     &                                    DIFABX
     &                                 -  DIFABY*SY0*
     &                                   (SX11*SZ1N - SX2N*DZ0))*D4INV
         SHLINT(INT,6) = SHLINT(INT,6) + (DIFABX*SZ0*
     &                                   (SY11*SX1N - SY2N*DX0)
     &                                 -  DIFABY*SZ0*
     &                                   (SX2N*DY0  - SX11*SY1N))*D2INV
 100  CONTINUE
      RETURN
      END
C  /* Deck dsuslh */
      SUBROUTINE DSUSLH(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,DOEKIN)
C
C     K.Ruud, Feb. 1992
C
C     H.J.Aa.Jensen (&MI) Mar 2003: new flag DOEKIN, if false omit
C     kinetic energy terms
C     (used by Dirac for DSUSLH because Vnuc is same in 4-comp, but kin.en. is not)

#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D4INV = 0.250D0, D8INV = 0.1250D0, D2 = 2.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6), AHGTF(*)
      LOGICAL   DOEKIN
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
C  .... when wished, calculate the kinetic energy contribution
      IF (DOEKIN) THEN
         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1  = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1  = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1  = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2  = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2  = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2  = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
         SX21 = SHGTF*ODC(LVALA,LVALB,0,2,1,1)
         SY21 = SHGTF*ODC(MVALA,MVALB,0,2,1,2)
         SZ21 = SHGTF*ODC(NVALA,NVALB,0,2,1,3)
         SX22 = SHGTF*ODC(LVALA,LVALB,0,2,2,1)
         SY22 = SHGTF*ODC(MVALA,MVALB,0,2,2,2)
         SZ22 = SHGTF*ODC(NVALA,NVALB,0,2,2,3)
         DX2  = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DY2  = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DZ2  = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         XX   = SX22*SY0*SZ0 + SX2*DY2*SZ0 + SX2*SY0*DZ2
         XY   = SX21*SY1*SZ0 + SX1*SY21*SZ0 + SX1*SY1*DZ2
         XZ   = SX21*SY0*SZ1 + SX1*DY2*SZ1 + SX1*SY0*SZ21
         YY   = DX2*SY2*SZ0 + SX0*SY22*SZ0 + SX0*SY2*DZ2
         YZ   = DX2*SY1*SZ1 + SX0*SY21*SZ1 + SX0*SY1*SZ21
         ZZ   = DX2*SY0*SZ2 + SX0*DY2*SZ2 + SX0*SY0*SZ22
C
C       Kinetic energy contribution
C
         SHLINT(INT,1) = SHLINT(INT,1) + (D2*DIFABY*DIFABZ*YZ
     &                                 -  DIFABY*DIFABY*ZZ
     &                                 -  DIFABZ*DIFABZ*YY)*D8INV
         SHLINT(INT,2) = SHLINT(INT,2) + (DIFABZ*DIFABZ*XY
     &                                 +  DIFABY*DIFABX*ZZ
     &                                 -  DIFABZ*DIFABX*YZ
     &                                 -  DIFABY*DIFABZ*XZ)*D8INV
         SHLINT(INT,3) = SHLINT(INT,3) + (DIFABX*DIFABZ*YY
     &                                 +  DIFABY*DIFABY*XZ
     &                                 -  DIFABY*DIFABZ*XY
     &                                 -  DIFABX*DIFABY*YZ)*D8INV
         SHLINT(INT,4) = SHLINT(INT,4) + (D2*DIFABX*DIFABZ*XZ
     &                                 -  DIFABX*DIFABX*ZZ
     &                                 -  DIFABZ*DIFABZ*XX)*D8INV
         SHLINT(INT,5) = SHLINT(INT,5) + (DIFABX*DIFABX*YZ
     &                                 +  DIFABY*DIFABZ*XX
     &                                 -  DIFABX*DIFABY*XZ
     &                                 -  DIFABX*DIFABZ*XY)*D8INV
         SHLINT(INT,6) = SHLINT(INT,6) + (DIFABX*DIFABY*XY*D2
     &                                 -  DIFABY*DIFABY*XX
     &                                 -  DIFABX*DIFABX*YY)*D8INV

      ENDIF
C
C     Nuclear attraction contribution
C         ..... calculate it always

         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV + 2
            IADRAU = IADRAV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            GV = ODC(NVALA,NVALB,IV,0,2,3)
            IUMAX = MAXU + 2
            IF (IV .GT. MAXV + 1) THEN
               IUMAX = MAXU
            ELSE IF (IV .GT. MAXV) THEN
               IUMAX = MAXU + 1
            END IF
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               GU = ODC(MVALA,MVALB,IU,0,2,2)
               ITMAX = MAXT + 2
               IF (((IU .GT. MAXU + 1) .OR. (IV .GT. MAXV + 1))
     &               .OR. ((IU .EQ. MAXU + 1) .AND. (IV .EQ. MAXV + 1)))
     &               THEN
                  ITMAX = MAXT
               ELSE IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) THEN
                  ITMAX = ITMAX + 1
               END IF
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  GT = ODC(LVALA,LVALB,IT,0,2,1)
                  FXX = GT*EU*EV
                  FXY = FT*FU*EV
                  FXZ = FT*EU*FV
                  FYY = ET*GU*EV
                  FYZ = ET*FU*FV
                  FZZ = ET*EU*GV
                  ATUV = -D4INV*AHGTF(IADRAU + IT)
                  SHLINT(INT,1)=SHLINT(INT,1)+ATUV*(DIFABY*DIFABZ*FYZ*D2
     &                                            - DIFABY*DIFABY*FZZ
     &                                            - DIFABZ*DIFABZ*FYY)
                  SHLINT(INT,2)=SHLINT(INT,2)+ATUV*(DIFABZ*DIFABZ*FXY
     &                                            + DIFABY*DIFABX*FZZ
     &                                            - DIFABZ*DIFABX*FYZ
     &                                            - DIFABY*DIFABZ*FXZ)
                  SHLINT(INT,3)=SHLINT(INT,3)+ATUV*(DIFABX*DIFABZ*FYY
     &                                            + DIFABY*DIFABY*FXZ
     &                                            - DIFABY*DIFABZ*FXY
     &                                            - DIFABX*DIFABY*FYZ)
                  SHLINT(INT,4)=SHLINT(INT,4)+ATUV*(DIFABX*DIFABZ*FXZ*D2
     &                                            - DIFABX*DIFABX*FZZ
     &                                            - DIFABZ*DIFABZ*FXX)
                  SHLINT(INT,5)=SHLINT(INT,5)+ATUV*(DIFABX*DIFABX*FYZ
     &                                            + DIFABY*DIFABZ*FXX
     &                                            - DIFABX*DIFABY*FXZ
     &                                            - DIFABX*DIFABZ*FXY)
                  SHLINT(INT,6)=SHLINT(INT,6)+ATUV*(DIFABX*DIFABY*FXY*D2
     &                                            - DIFABY*DIFABY*FXX
     &                                            - DIFABX*DIFABX*FYY)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck nsnlin */
      SUBROUTINE NSNLIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C     Nuclear shielding tensor integrals without London orbital contribution
C
C     K. Ruud, March 1992
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  FEE = FT*EU*EV
                  EFE = ET*FU*EV
                  EEF = ET*EU*FV
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IXX = 9*(IATOM - 1) + 1
                     IYX = 9*(IATOM - 1) + 2
                     IZX = 9*(IATOM - 1) + 3
                     IXY = 9*(IATOM - 1) + 4
                     IYY = 9*(IATOM - 1) + 5
                     IZY = 9*(IATOM - 1) + 6
                     IXZ = 9*(IATOM - 1) + 7
                     IYZ = 9*(IATOM - 1) + 8
                     IZZ = 9*(IATOM - 1) + 9
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
C
                     SHLINT(INT,IXX) = SHLINT(INT,IXX) - D2INV*
     &                                 (EFE*AH0U + EEF*AH0V)
                     SHLINT(INT,IYX) = SHLINT(INT,IYX) + D2INV*FEE*AH0U
                     SHLINT(INT,IZX) = SHLINT(INT,IZX) + D2INV*FEE*AH0V
                     SHLINT(INT,IXY) = SHLINT(INT,IXY) + D2INV*EFE*AH0T
                     SHLINT(INT,IYY) = SHLINT(INT,IYY) - D2INV*
     *                                 (FEE*AH0T + EEF*AH0V)
                     SHLINT(INT,IZY) = SHLINT(INT,IZY) + D2INV*EFE*AH0V
                     SHLINT(INT,IXZ) = SHLINT(INT,IXZ) + D2INV*EEF*AH0T
                     SHLINT(INT,IYZ) = SHLINT(INT,IYZ) + D2INV*EEF*AH0U
                     SHLINT(INT,IZZ) = SHLINT(INT,IZZ) - D2INV*
     &                                 (FEE*AH0T + EFE*AH0U)
                     IOFF = IOFF + NAHGTF
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck nsloin */
      SUBROUTINE NSLOIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,DIFABX,DIFABY,DIFABZ)
C
C     London orbital contribution to nuclear shielding integrals
C
C     K. Ruud, March 1992
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 2
         MAXU = MVALA + MVALB + 2
         MAXV = NVALA + NVALB + 2
         IADRAV = 1
         DO 200 IV = 0, MAXV
            DV = ODC(NVALA,NVALB,IV,1,1,3)
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            GV = ODC(NVALA,NVALB,IV,0,1,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               DU = ODC(MVALA,MVALB,IU,1,1,2)
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               GU = ODC(MVALA,MVALB,IU,0,1,2)
               DO 400 IT = 0, MAXT
                  DT = ODC(LVALA,LVALB,IT,1,1,1)
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  GT = ODC(LVALA,LVALB,IT,0,1,1)
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IXX = 9*(IATOM - 1) + 1
                     IYX = 9*(IATOM - 1) + 2
                     IZX = 9*(IATOM - 1) + 3
                     IXY = 9*(IATOM - 1) + 4
                     IYY = 9*(IATOM - 1) + 5
                     IZY = 9*(IATOM - 1) + 6
                     IXZ = 9*(IATOM - 1) + 7
                     IYZ = 9*(IATOM - 1) + 8
                     IZZ = 9*(IATOM - 1) + 9
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
                     SHLINT(INT,IXX) = SHLINT(INT,IXX) - (DIFABY*ET*
     &                                 (FU*GV*AH0V - EU*DV*AH0U)
     &                               + (GU*FV*AH0U - DU*EV*AH0V)*
     &                                 DIFABZ*ET)*D2INV
                     SHLINT(INT,IYX) = SHLINT(INT,IYX) - (DIFABZ*GT*
     &                                 (FU*EV*AH0V - EU*FV*AH0U)
     &                               + (EU*DV*AH0U - FU*GV*AH0V)*
     &                                 DIFABX*ET)*D2INV
                     SHLINT(INT,IZX) = SHLINT(INT,IZX) - (DIFABX*ET*
     &                                 (DU*EV*AH0V - GU*FV*AH0U)
     &                               + (EU*FV*AH0U - FU*EV*AH0V)*
     &                                 DIFABY*GT)*D2INV
                     SHLINT(INT,IXY) = SHLINT(INT,IXY) + (DIFABY*EU*
     &                                 (FT*GV*AH0V - ET*DV*AH0T)
     &                               + (ET*FV*AH0T - FT*EV*AH0V)*
     &                                 DIFABZ*GU)*D2INV
                     SHLINT(INT,IYY) = SHLINT(INT,IYY) + (DIFABZ*EU*
     &                                 (DT*EV*AH0V - GT*FV*AH0T)
     &                               + (ET*DV*AH0T - FT*GV*AH0V)*
     &                                 DIFABX*EU)*D2INV
                     SHLINT(INT,IZY) = SHLINT(INT,IZY) + (DIFABX*GU*
     &                                 (FT*EV*AH0V - ET*FV*AH0T)
     &                               + (GT*FV*AH0T - DT*EV*AH0V)*
     &                                 DIFABY*EU)*D2INV
                     SHLINT(INT,IXZ) = SHLINT(INT,IXZ) - (DIFABY*GV*
     &                                 (FT*EU*AH0U - ET*FU*AH0T)
     &                               + (ET*DU*AH0T - FT*GU*AH0U)*
     &                                 DIFABZ*EV)*D2INV
                     SHLINT(INT,IYZ) = SHLINT(INT,IYZ) - (DIFABZ*EV*
     &                                 (DT*EU*AH0U - GT*FU*AH0T)
     &                               + (ET*FU*AH0T - FT*EU*AH0U)*
     &                                 DIFABX*GV)*D2INV
                     SHLINT(INT,IZZ) = SHLINT(INT,IZZ) - (DIFABX*EV*
     &                                 (FT*GU*AH0U - ET*DU*AH0T)
     &                               + (GT*FU*AH0T - DT*EU*AH0U)*
     &                                 DIFABY*EV)*D2INV
                     IOFF = IOFF + NAHGTF
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck npeint */
      SUBROUTINE NPEINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     Potential energy of interaction of electrons with individual nuclei
C
C     K.Ruud, July 1992
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &              SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
C
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  EEE = ET*EU*EV
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     SHLINT(INT,IATOM) = SHLINT(INT,IATOM)
     &                                 + EEE*AHGTF(IOFF + IADRAU + IT)
                     IOFF = IOFF + NAHGTF
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck efgint */
      SUBROUTINE EFGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     Electric field gradient integrals
C
C     K.Ruud, June 1992
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxmom.h"
#include "maxorb.h"
      PARAMETER (D1 = 1.0D0, D3 = 3.0D0, D3INV = D1/D3)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &              SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "nuclei.h"
#include "symmet.h"
#include "onecom.h"
#include "lmns.h"
#include "ibtfun.h"
C
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  EEE = ET*EU*EV
                  IOFF = 0
                  INTTYP = 0
                  DO 500 IATOM = 1, NATOMC
                     IADR0  = IOFF  + IADRAU + IT
                     IADRTT = IADR0 + 2
                     IADRTU = IADR0 + 1 + ISTEPU
                     IADRTV = IADR0 + 1 + ISTEPV
                     IADRUU = IADR0 + 2*ISTEPU
                     IADRUV = IADR0 + ISTEPU + ISTEPV
                     IADRVV = IADR0 + 2*ISTEPV
#ifndef NON_TRACELESS_EFG
                     SHLINT(INT,INTTYP+1) = SHLINT(INT,INTTYP+1)
     &                                  + EEE*D3INV*(2*AHGTF(IADRTT)
     &                                  - AHGTF(IADRUU) - AHGTF(IADRVV))
                     SHLINT(INT,INTTYP+2) = SHLINT(INT,INTTYP+2)
     &                                  + EEE*AHGTF(IADRTU)
                     SHLINT(INT,INTTYP+3) = SHLINT(INT,INTTYP+3)
     &                                  + EEE*AHGTF(IADRTV)
                     SHLINT(INT,INTTYP+4) = SHLINT(INT,INTTYP+4)
     &                                  + EEE*D3INV*(2*AHGTF(IADRUU)
     &                                  - AHGTF(IADRTT) - AHGTF(IADRVV))
                     SHLINT(INT,INTTYP+5) = SHLINT(INT,INTTYP+5)
     &                                  + EEE*AHGTF(IADRUV)
                     SHLINT(INT,INTTYP+6) = SHLINT(INT,INTTYP+6)
     &                                  + EEE*D3INV*(2*AHGTF(IADRVV)
     &                                  - AHGTF(IADRTT) - AHGTF(IADRUU))
#else
                     SHLINT(INT,INTTYP+1) = SHLINT(INT,INTTYP+1)
     &                                  + EEE*AHGTF(IADRTT)
                     SHLINT(INT,INTTYP+2) = SHLINT(INT,INTTYP+2)
     &                                  + EEE*AHGTF(IADRTU)
                     SHLINT(INT,INTTYP+3) = SHLINT(INT,INTTYP+3)
     &                                  + EEE*AHGTF(IADRTV)
                     SHLINT(INT,INTTYP+4) = SHLINT(INT,INTTYP+4)
     &                                  + EEE*AHGTF(IADRUU)
                     SHLINT(INT,INTTYP+5) = SHLINT(INT,INTTYP+5)
     &                                  + EEE*AHGTF(IADRUV)
                     SHLINT(INT,INTTYP+6) = SHLINT(INT,INTTYP+6)
     &                                  + EEE*AHGTF(IADRVV)
#endif
                     IOFF = IOFF + NAHGTF
                     INTTYP = INTTYP + 6
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck eftint */
      SUBROUTINE EFTINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC,TRLESS,ONCN,
     &                  CORCX,CORCY,CORCZ)
C
C     Electric field third derivatives.
C
C     Calculate the traceless (TRLESS = true) or 
C     non-traceless (TRLESS = false).
C
C     ONCN: one-center (coor a = coor b)
C     Don't calculate < s_A | d^4/dx^4 1/r_A | s_A >
C     
C
C     In Dirac the routine PRP_EFT calculates the 'traceless' tensor
C     using the Laplace equation. But for numerical reasons
C     it is best (I guess) to calculate the traceless one
C     directly.
C
C     Adapted from EFGINT.
C
C     J. Thyssen, September 30, 1999
C
C
#include "implicit.h"
      PARAMETER (D0  =  0.0D00, D1 = 1.0D00, D3 = 3.0D00,
     &           D4  =  4.0D00, D7 = 7.0D00,D27=27.0D00,
     &           D35 = 35.0D00, DTHR = 1.0D-10)
      PARAMETER (D8=8.0D0,D24=24.0D0,D6=6.0D0)

#include "maxaqn.h"
#include "mxcent.h"
#include "maxmom.h"
#include "maxorb.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &              SHLINT(KCKTAB,NOPTYP), AHGTF(*),
     &              CORCX(*),CORCY(*),CORCZ(*)
      DIMENSION II(0:2)
      LOGICAL   TRLESS,ONCN
#include "nuclei.h"
#include "symmet.h"
#include "onecom.h"
#include "lmns.h"
#include "ibtfun.h"
C
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  EEE = ET*EU*EV
                  IOFF = 0
                  INTTYP = 0
                  DO 500 IATOM = 1, NATOMC
C
C
C
                     IADR0 = IOFF + IADRAU + IT
                     VXXXX = AHGTF(IADR0 + 4)
                     VXXXY = AHGTF(IADR0 + 3 +     ISTEPU)
                     VXXXZ = AHGTF(IADR0 + 3              +     ISTEPV)
                     VXXYY = AHGTF(IADR0 + 2 + 2 * ISTEPU)
                     VXXYZ = AHGTF(IADR0 + 2 +     ISTEPU +     ISTEPV)
                     VXXZZ = AHGTF(IADR0 + 2              + 2 * ISTEPV)
                     VXYYY = AHGTF(IADR0 + 1 + 3 * ISTEPU)
                     VXYYZ = AHGTF(IADR0 + 1 + 2 * ISTEPU +     ISTEPV)
                     VXYZZ = AHGTF(IADR0 + 1 +     ISTEPU + 2 * ISTEPV)
                     VXZZZ = AHGTF(IADR0 + 1              + 3 * ISTEPV)
                     VYYYY = AHGTF(IADR0     + 4 * ISTEPU)
                     VYYYZ = AHGTF(IADR0     + 3 * ISTEPU + 1 * ISTEPV)
                     VYYZZ = AHGTF(IADR0     + 2 * ISTEPU + 2 * ISTEPV)
                     VYZZZ = AHGTF(IADR0     + 1 * ISTEPU + 3 * ISTEPV)
                     VZZZZ = AHGTF(IADR0                  + 4 * ISTEPV)
C
                     IF (ONCN) THEN
                        IF ( (CORCX(IATOM) .EQ. CORAX) 
     &                     .AND. (CORCY(IATOM) .EQ. CORAY) 
     &                     .AND. (CORCZ(IATOM) .EQ. CORAZ) 
     &                     .AND. (MAXU+MAXT+MAXV .LE. 3)) THEN
#ifdef UNDEF
                           dbla = 
     &                  (  D8 * VZZZZ 
     &                  - D24 * VXXZZ
     &                  - D24 * VYYZZ
     &                  +  D6 * VXXYY
     &                  +  D3 * VXXXX
     &                  +  D3 * VYYYY ) / D35
c                          if (dbla .ge. 1.0d-12) then
                              write(7,*) 'jth maxtuv ',lvala+mvala+nvala,
     &                            lvalb+mvalb+nvalb
                              write(7,*) 'jth tuv    ',it,iu,iv
                              write(7,*) 'jth cont   ',dbla
                              write(7,*) 'jth odc    ',eee
                              write(7,*) 'jth c * odc',eee*dbla
                              write(7,*) 'jth indi1  ',
     &                           vzzzz,vxxzz,vyyzz
                              write(7,*) 'jth indi2  ',
     &                           vxxyy,vxxxx,vyyyy
c                          end if
#endif
                           GOTO 600
                        END IF
                     END IF
                     IF (TRLESS) THEN
#ifndef OLD_CODE
                        SHLINT(INT,INTTYP+ 1) = SHLINT(INT,INTTYP+ 1) 
     &                     + EEE * (
     &                        D8 * VXXXX 
     &                     - D24 * VXXYY
     &                     - D24 * VXXZZ
     &                     +  D3 * VYYYY
     &                     +  D6 * VYYZZ
     &                     +  D3 * VZZZZ ) / D35
                        SHLINT(INT,INTTYP+ 2) = SHLINT(INT,INTTYP+ 2) 
     &                     + EEE * (
     &                        D4 * VXXXY 
     &                     -  D3 * VXYYY
     &                     -  D3 * VXYZZ ) / D7
                        SHLINT(INT,INTTYP+ 3) = SHLINT(INT,INTTYP+ 3) 
     &                     + EEE * (
     &                        D4 * VXXXZ 
     &                     -  D3 * VXYYZ
     &                     -  D3 * VXZZZ ) / D7
                        SHLINT(INT,INTTYP+ 4) = SHLINT(INT,INTTYP+ 4) 
     &                     + EEE * (
     &                     -  D4 * VXXXX 
     &                     + D27 * VXXYY
     &                     -  D3 * VXXZZ
     &                     -  D4 * VYYYY
     &                     -  D3 * VYYZZ
     &                     +       VZZZZ ) / D35
                        SHLINT(INT,INTTYP+ 5) = SHLINT(INT,INTTYP+ 5) 
     &                     + EEE * (
     &                        D6 * VXXYZ 
     &                     -       VYYYZ
     &                     -       VYZZZ ) / D7
                        SHLINT(INT,INTTYP+ 6) = SHLINT(INT,INTTYP+ 6) 
     &                     + EEE * (
     &                     -  D4 * VXXXX 
     &                     -  D3 * VXXYY
     &                     + D27 * VXXZZ
     &                     +       VYYYY
     &                     -  D3 * VYYZZ
     &                     -  D4 * VZZZZ ) / D35
                        SHLINT(INT,INTTYP+ 7) = SHLINT(INT,INTTYP+ 7) 
     &                     + EEE * (
     &                     -  D3 * VXXXY 
     &                     +  D4 * VXYYY
     &                     -  D3 * VXYZZ ) / D7
                        SHLINT(INT,INTTYP+ 8) = SHLINT(INT,INTTYP+ 8) 
     &                     + EEE * (
     &                        D6 * VXYYZ 
     &                     -       VXXXZ
     &                     -       VXZZZ ) / D7
                        SHLINT(INT,INTTYP+ 9) = SHLINT(INT,INTTYP+ 9) 
     &                     + EEE * (
     &                        D6 * VXYZZ 
     &                     -       VXXXY
     &                     -       VXYYY ) / D7
                        SHLINT(INT,INTTYP+10) = SHLINT(INT,INTTYP+10) 
     &                     + EEE * (
     &                        D4 * VXZZZ 
     &                     -  D3 * VXXXZ
     &                     -  D3 * VXYYZ ) / D7
                        SHLINT(INT,INTTYP+11) = SHLINT(INT,INTTYP+11) 
     &                     + EEE * (
     &                        D8 * VYYYY 
     &                     - D24 * VXXYY
     &                     - D24 * VYYZZ
     &                     +  D6 * VXXZZ
     &                     +  D3 * VXXXX
     &                     +  D3 * VZZZZ ) / D35
                        SHLINT(INT,INTTYP+12) = SHLINT(INT,INTTYP+12) 
     &                     + EEE * (
     &                        D4 * VYYYZ 
     &                     -  D3 * VXXYZ
     &                     -  D3 * VYZZZ ) / D7
                        SHLINT(INT,INTTYP+13) = SHLINT(INT,INTTYP+13) 
     &                     + EEE * (
     &                     -  D4 * VYYYY 
     &                     -  D3 * VXXYY
     &                     + D27 * VYYZZ
     &                     +       VXXXX
     &                     -  D3 * VXXZZ
     &                     -  D4 * VZZZZ ) / D35
                        SHLINT(INT,INTTYP+14) = SHLINT(INT,INTTYP+14) 
     &                     + EEE * (
     &                        D4 * VYZZZ 
     &                     -  D3 * VXXYZ
     &                     -  D3 * VYYYZ ) / D7
                        SHLINT(INT,INTTYP+15) = SHLINT(INT,INTTYP+15) 
     &                     + EEE * (
     &                        D8 * VZZZZ 
     &                     - D24 * VXXZZ
     &                     - D24 * VYYZZ
     &                     +  D6 * VXXYY
     &                     +  D3 * VXXXX
     &                     +  D3 * VYYYY ) / D35
#else
                        VXXAA = VXXXX + VXXYY + VXXZZ
                        VXYAA = VXXXY + VXYYY + VXYZZ
                        VXZAA = VXXXZ + VXYYZ + VXZZZ
                        VYYAA = VXXYY + VYYYY + VYYZZ
                        VYZAA = VXXYZ + VYYYZ + VYZZZ
                        VZZAA = VXXZZ + VYYZZ + VZZZZ
                        SHLINT(INT,INTTYP+ 1) = SHLINT(INT,INTTYP+ 1) 
     &                     + EEE * (
     &                       VXXXX 
     &                     - D27/D35 * VXXAA
     &                     +  D3/D35 * VYYAA
     &                     +  D3/D35 * VZZAA )
                        SHLINT(INT,INTTYP+ 2) = SHLINT(INT,INTTYP+ 2) 
     &                     + EEE * (
     &                       VXXXY 
     &                     -  D3/D7  * VXYAA )
                        SHLINT(INT,INTTYP+ 3) = SHLINT(INT,INTTYP+ 3) 
     &                     + EEE * (
     &                       VXXXZ 
     &                     -  D3/D7  * VXZAA )
                        SHLINT(INT,INTTYP+ 4) = SHLINT(INT,INTTYP+ 4) 
     &                     + EEE * (
     &                       VXXYY 
     &                     -  D4/D35 * VXXAA
     &                     -  D4/D35 * VYYAA
     &                     +  D1/D35 * VZZAA )
                        SHLINT(INT,INTTYP+ 5) = SHLINT(INT,INTTYP+ 5) 
     &                     + EEE * (
     &                       VXXYZ 
     &                     -  D1/D7  * VYZAA )
                        SHLINT(INT,INTTYP+ 6) = SHLINT(INT,INTTYP+ 6) 
     &                     + EEE * (
     &                       VXXZZ 
     &                     -  D4/D35 * VXXAA
     &                     +  D1/D35 * VYYAA
     &                     -  D4/D35 * VZZAA )
                        SHLINT(INT,INTTYP+ 7) = SHLINT(INT,INTTYP+ 7) 
     &                     + EEE * (
     &                       VXYYY 
     &                     -  D3/D7  * VXYAA )
                        SHLINT(INT,INTTYP+ 8) = SHLINT(INT,INTTYP+ 8) 
     &                     + EEE * (
     &                       VXYYZ 
     &                     -  D1/D7  * VXZAA )
                        SHLINT(INT,INTTYP+ 9) = SHLINT(INT,INTTYP+ 9) 
     &                     + EEE * (
     &                       VXYZZ 
     &                     -  D1/D7  * VXYAA )
                        SHLINT(INT,INTTYP+10) = SHLINT(INT,INTTYP+10) 
     &                     + EEE * (
     &                       VXZZZ 
     &                     -  D3/D7  * VXZAA )
                        SHLINT(INT,INTTYP+11) = SHLINT(INT,INTTYP+11) 
     &                     + EEE * (
     &                       VYYYY 
     &                     +  D3/D35 * VXXAA
     &                     - D27/D35 * VYYAA
     &                     +  D3/D35 * VZZAA )
                        SHLINT(INT,INTTYP+12) = SHLINT(INT,INTTYP+12) 
     &                     + EEE * (
     &                       VYYYZ 
     &                     -  D3/D7  * VYZAA )
                        SHLINT(INT,INTTYP+13) = SHLINT(INT,INTTYP+13) 
     &                     + EEE * (
     &                       VYYZZ 
     &                     +  D1/D35 * VXXAA
     &                     -  D4/D35 * VYYAA
     &                     -  D4/D35 * VZZAA )
                        SHLINT(INT,INTTYP+14) = SHLINT(INT,INTTYP+14) 
     &                     + EEE * (
     &                       VYZZZ 
     &                     -  D3/D7  * VYZAA )
                        SHLINT(INT,INTTYP+15) = SHLINT(INT,INTTYP+15) 
     &                     + EEE * (
     &                       VZZZZ
     &                     +  D3/D35 * VXXAA
     &                     +  D3/D35 * VYYAA
     &                     - D27/D35 * VZZAA )
#endif
                     ELSE
                        SHLINT(INT,INTTYP+ 1) =
     &                     SHLINT(INT,INTTYP+ 1) + EEE * VXXXX
                        SHLINT(INT,INTTYP+ 2) =
     &                     SHLINT(INT,INTTYP+ 2) + EEE * VXXXY
                        SHLINT(INT,INTTYP+ 3) =
     &                     SHLINT(INT,INTTYP+ 3) + EEE * VXXXZ
                        SHLINT(INT,INTTYP+ 4) =
     &                     SHLINT(INT,INTTYP+ 4) + EEE * VXXYY
                        SHLINT(INT,INTTYP+ 5) =
     &                     SHLINT(INT,INTTYP+ 5) + EEE * VXXYZ
                        SHLINT(INT,INTTYP+ 6) =
     &                     SHLINT(INT,INTTYP+ 6) + EEE * VXXZZ
                        SHLINT(INT,INTTYP+ 7) =
     &                     SHLINT(INT,INTTYP+ 7) + EEE * VXYYY
                        SHLINT(INT,INTTYP+ 8) =
     &                     SHLINT(INT,INTTYP+ 8) + EEE * VXYYZ
                        SHLINT(INT,INTTYP+ 9) =
     &                     SHLINT(INT,INTTYP+ 9) + EEE * VXYZZ
                        SHLINT(INT,INTTYP+10) =
     &                     SHLINT(INT,INTTYP+10) + EEE * VXZZZ
                        SHLINT(INT,INTTYP+11) =
     &                     SHLINT(INT,INTTYP+11) + EEE * VYYYY
                        SHLINT(INT,INTTYP+12) =
     &                     SHLINT(INT,INTTYP+12) + EEE * VYYYZ
                        SHLINT(INT,INTTYP+13) =
     &                     SHLINT(INT,INTTYP+13) + EEE * VYYZZ
                        SHLINT(INT,INTTYP+14) =
     &                     SHLINT(INT,INTTYP+14) + EEE * VYZZZ
                        SHLINT(INT,INTTYP+15) =
     &                     SHLINT(INT,INTTYP+15) + EEE * VZZZZ
                     END IF
 600                 CONTINUE
                     IOFF = IOFF + NAHGTF
                     INTTYP = INTTYP + 15
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
c     CALL HEADER('EFTINT - integrals',-1)
c     CALL OUTPUT(SHLINT,1,KCKTAB,1,NOPTYP,KCKTAB,NOPTYP,1,6)
      RETURN
      END
C  /* Deck ef1int */
      SUBROUTINE EF1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC,INTTYP)
C
C     Nuclear electric field integrals
C
C     K.Ruud, March 1992
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  EEE = ET*EU*EV
                  IF (INTTYP .EQ. 58) THEN
C                 ... total field from all nuclei
                     AH0T = AHGTF(IADRAU + IT + 1)
                     AH0U = AHGTF(IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IADRAU + IT + ISTEPV)
                     SHLINT(INT,1) = SHLINT(INT,1) - EEE*AH0T
                     SHLINT(INT,2) = SHLINT(INT,2) - EEE*AH0U
                     SHLINT(INT,3) = SHLINT(INT,3) - EEE*AH0V
                  ELSE
C                 ... (INTTYP 29): field from individual nuclei
                    IOFF = 0
                    DO 500 IATOM = 1, NATOMC
                     IX = 3*(IATOM - 1) + 1
                     IY = 3*(IATOM - 1) + 2
                     IZ = 3*(IATOM - 1) + 3
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
                     SHLINT(INT,IX) = SHLINT(INT,IX) - EEE*AH0T
                     SHLINT(INT,IY) = SHLINT(INT,IY) - EEE*AH0U
                     SHLINT(INT,IZ) = SHLINT(INT,IZ) - EEE*AH0V
                     IOFF = IOFF + NAHGTF
 500                CONTINUE
                  END IF
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck qdpint */
      SUBROUTINE QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C     tuh 1989
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D1P5 = 1.5D0, D4 = 4.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
C
         INT = INT + 1
         IF (INTTYP .EQ. 6) THEN
            SHLINT(INT,1) = SHLINT(INT,1) - DX2*SY0*SZ0
            SHLINT(INT,2) = SHLINT(INT,2) - DX1*DY1*SZ0
            SHLINT(INT,3) = SHLINT(INT,3) - DX1*SY0*DZ1
            SHLINT(INT,4) = SHLINT(INT,4) - SX0*DY2*SZ0
            SHLINT(INT,5) = SHLINT(INT,5) - SX0*DY1*DZ1
            SHLINT(INT,6) = SHLINT(INT,6) - SX0*SY0*DZ2
         ELSE IF (INTTYP .EQ. 7) THEN
            X2 = DX2*SY0*SZ0
            Y2 = SX0*DY2*SZ0
            Z2 = SX0*SY0*DZ2
            R2 = (X2 + Y2 + Z2)/D2
            SHLINT(INT,1) = SHLINT(INT,1) - D1P5*X2 + R2
            SHLINT(INT,2) = SHLINT(INT,2) - D1P5*DX1*DY1*SZ0
            SHLINT(INT,3) = SHLINT(INT,3) - D1P5*DX1*SY0*DZ1
            SHLINT(INT,4) = SHLINT(INT,4) - D1P5*Y2 + R2
            SHLINT(INT,5) = SHLINT(INT,5) - D1P5*SX0*DY1*DZ1
            SHLINT(INT,6) = SHLINT(INT,6) - D1P5*Z2 + R2
         ELSE
            X2 = DX2*SY0*SZ0
            Y2 = SX0*DY2*SZ0
            Z2 = SX0*SY0*DZ2
            R2 = (X2 + Y2 + Z2)/D4
            SHLINT(INT,1) = SHLINT(INT,1) - R2 + X2/D4
            SHLINT(INT,2) = SHLINT(INT,2) + DX1*DY1*SZ0/D4
            SHLINT(INT,3) = SHLINT(INT,3) + DX1*SY0*DZ1/D4
            SHLINT(INT,4) = SHLINT(INT,4) - R2 + Y2/D4
            SHLINT(INT,5) = SHLINT(INT,5) + SX0*DY1*DZ1/D4
            SHLINT(INT,6) = SHLINT(INT,6) - R2 + Z2/D4
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck momint */
      SUBROUTINE MOMINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,ORIGIN,IORDER,WORK,
     &                  LWORK)
C
C     tuh 1989
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      DIMENSION WORK(LWORK)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          ORIGIN(3), SHLINT(KCKTAB,(IORDER+1)*(IORDER+2)/2)
#include "onecom.h"
C
      KHMULT = 1
      KDX    = KHMULT + 3*((IORDER + 1)**2)
      KDY    = KDX    + IORDER + 1
      KDZ    = KDY    + IORDER + 1
      KLO    = KDZ    + IORDER + 1
      KMO    = KLO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KNO    = KMO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KLAST  = KNO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('MOMINT',' ',KLAST,LWORK)
      CALL MOMIN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,CORPX,
     &            CORPY,CORPZ,EXPPI,ORIGIN,IORDER,WORK(KHMULT),
     &            WORK(KDX),WORK(KDY),WORK(KDZ),WORK(KLO),WORK(KMO),
     &            WORK(KNO))
      RETURN
      END
C  /* Deck momin1 */
      SUBROUTINE MOMIN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,ORIGIN,IORDER,HMULT,DX,
     &                  DY,DZ,LO,MO,NO)
C
C     tuh 1989
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          ORIGIN(3), SHLINT(KCKTAB,(IORDER+1)*(IORDER+2)/2),
     &          HMULT(3,0:IORDER,0:IORDER),
     &          DX(0:IORDER), DY(0:IORDER), DZ(0:IORDER),
     &          LO((IORDER + 1)*(IORDER + 2)/2),
     &          MO((IORDER + 1)*(IORDER + 2)/2),
     &          NO((IORDER + 1)*(IORDER + 2)/2)
#include "onecom.h"
#include "lmns.h"
C
C
C     Hermitian integrals
C     -------------------
C
      HMULT(1,0,0) = SHGTF
      HMULT(2,0,0) = SHGTF
      HMULT(3,0,0) = SHGTF
      DO 100 IO = 1, IORDER
         DO 200 IT = 0, IO
            HX = D0
            HY = D0
            HZ = D0
            IF (IT .GT. 0) THEN
               HX = HX + dble(IT)*HMULT(1,IO-1,IT-1)
               HY = HY + dble(IT)*HMULT(2,IO-1,IT-1)
               HZ = HZ + dble(IT)*HMULT(3,IO-1,IT-1)
            END IF
            IF (IT .LE. IO-1) THEN
               HX = HX + (CORPX - ORIGIN(1))*HMULT(1,IO-1,IT)
               HY = HY + (CORPY - ORIGIN(2))*HMULT(2,IO-1,IT)
               HZ = HZ + (CORPZ - ORIGIN(3))*HMULT(3,IO-1,IT)
            END IF
            IF (IT .LE. IO-2) THEN
               HX = HX + DP5*EXPPI*HMULT(1,IO-1,IT+1)
               HY = HY + DP5*EXPPI*HMULT(2,IO-1,IT+1)
               HZ = HZ + DP5*EXPPI*HMULT(3,IO-1,IT+1)
            END IF
            HMULT(1,IO,IT) = HX
            HMULT(2,IO,IT) = HY
            HMULT(3,IO,IT) = HZ
  200    CONTINUE
  100 CONTINUE
C
C     Cartesian integrals
C
      INT = 0
      DO 300 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 300 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C        One-dimensional integrals
C
         DO 400 IO = 0, IORDER
            DX(IO) = D0
            DY(IO) = D0
            DZ(IO) = D0
            DO 500 IT = 0, MIN(LVALA+LVALB,IO)
               DX(IO) = DX(IO)+ ODC(LVALA,LVALB,IT,0,0,1)*HMULT(1,IO,IT)
  500       CONTINUE
            DO 510 IU = 0, MIN(MVALA+MVALB,IO)
               DY(IO) = DY(IO)+ ODC(MVALA,MVALB,IU,0,0,2)*HMULT(2,IO,IU)
  510       CONTINUE
            DO 520 IV = 0, MIN(NVALA+NVALB,IO)
               DZ(IO) = DZ(IO)+ ODC(NVALA,NVALB,IV,0,0,3)*HMULT(3,IO,IV)
  520       CONTINUE
  400    CONTINUE
C
C        Three-dimensional integrals
C
         INT = INT + 1
         CALL LMNVAL(IORDER+1,(IORDER + 1)*(IORDER + 2)/2,LO,MO,NO)
         DO 600 I = 1, (IORDER + 1)*(IORDER + 2)/2
            SHLINT(INT,I) = SHLINT(INT,I)-DX(LO(I))*DY(MO(I))*DZ(NO(I))
  600    CONTINUE
  300 CONTINUE
      RETURN
      END
C  /* Deck dpvint */
      SUBROUTINE DPVINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     tuh Jan. 17 1990
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         INT = INT + 1
         SHLINT(INT,1) = SHLINT(INT,1) - DX0*SY0*SZ0
         SHLINT(INT,2) = SHLINT(INT,2) - SX0*DY0*SZ0
         SHLINT(INT,3) = SHLINT(INT,3) - SX0*SY0*DZ0
  100 CONTINUE
      RETURN
      END

      subroutine g1k_integrals(ix,
     &                         iy,
     &                         iz,
     &                         jmaxa,
     &                         jmaxb,
     &                         jmaxt,
     &                         jmaxd,
     &                         jmaxm,
     &                         odc,
     &                         shgtf,
     &                         shlint)

      implicit none

#include "maxaqn.h"
#include "onecom.h"
#include "lmns.h"

!     --------------------------------------------------------------------------
      integer, intent(in)    :: ix
      integer, intent(in)    :: iy
      integer, intent(in)    :: iz
      integer, intent(in)    :: jmaxa
      integer, intent(in)    :: jmaxb
      integer, intent(in)    :: jmaxt
      integer, intent(in)    :: jmaxd
      integer, intent(in)    :: jmaxm
      real(8), intent(in)    ::
     &               odc(0:jmaxa, 0:jmaxb, 0:jmaxt, 0:jmaxd, 0:jmaxm, 3)
      real(8), intent(in)    :: shgtf
      real(8), intent(inout) :: shlint(kcktab, 3)
!     --------------------------------------------------------------------------
      integer                :: i, ia, ib, id, ixyz
      integer                :: lmn(2, 3)
      real(8)                :: derv(3, 0:2)
      real(8)                :: tx, ty, tz
!     --------------------------------------------------------------------------

      i = 0
      do ia = 1,kckta

        lmn(1, 1) = lvalua(ia)
        lmn(1, 2) = mvalua(ia)
        lmn(1, 3) = nvalua(ia)

        do ib = 1, kcktb
          i = i + 1

          lmn(2, 1) = lvalub(ib)
          lmn(2, 2) = mvalub(ib)
          lmn(2, 3) = nvalub(ib)

          do ixyz = 1, 3
            do id = 0, 2
              derv(ixyz, id) = 
     &                   odc(lmn(1, ixyz), lmn(2, ixyz), 0, id, 0, ixyz)
            end do
          end do
          derv = shgtf*derv

          tx = derv(1, ix+1)*derv(2, iy+0)*derv(3, iz+0)
          ty = derv(1, ix+0)*derv(2, iy+1)*derv(3, iz+0)
          tz = derv(1, ix+0)*derv(2, iy+0)*derv(3, iz+1)

          shlint(i, 1) = shlint(i, 1) - tx
          shlint(i, 2) = shlint(i, 2) - ty
          shlint(i, 3) = shlint(i, 3) - tz

        end do
      end do

      end subroutine

C  /* Deck frmint */
      SUBROUTINE FRMINT(SHLINT,NOPTYP,EXPP,XP,YP,ZP,SAAB,DOATOM)
C
C     tuh Feb 91
C
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "gfac.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
      PARAMETER (D0 = 0.0D0, D4 = 4.0D0, D3 = 3.0D0)
      PARAMETER (DFAC = D4*GFAC*PI/D3)
      LOGICAL DOATOM(NUCIND)
      DIMENSION SHLINT(KCKTAB,NOPTYP)
#include "onecom.h"
#include "lmns.h"
#include "symmet.h"
#include "nuclei.h"
C
#include "ibtfun.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         INTTYP = 0
         INT = INT + 1
         DO 200 IREP = 0, MAXREP
            DO 300 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
            IF (IBTAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
               INTTYP = INTTYP + 1
               FRM = D0
               DO 400 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
                  XC = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
                  YC = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
                  ZC = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
                  RCP2 = (XC-XP)**2 + (YC-YP)**2 + (ZC-ZP)**2
                  FINT = PT(IBTAND(IREP,ISYMOP))*EXP(-EXPP*RCP2)
                  IF (LVALA .GT. 0) FINT = FINT*((XC - CORAX)**LVALA)
                  IF (MVALA .GT. 0) FINT = FINT*((YC - CORAY)**MVALA)
                  IF (NVALA .GT. 0) FINT = FINT*((ZC - CORAZ)**NVALA)
                  IF (LVALB .GT. 0) FINT = FINT*((XC - CORBX)**LVALB)
                  IF (MVALB .GT. 0) FINT = FINT*((YC - CORBY)**MVALB)
                  IF (NVALB .GT. 0) FINT = FINT*((ZC - CORBZ)**NVALB)
                  FRM = FRM + FINT
               END IF
  400          CONTINUE
               SHLINT(INT,INTTYP) = SHLINT(INT,INTTYP) - SAAB*FRM*DFAC
            END IF
            END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck dwnint */
      SUBROUTINE DWNINT(SHLINT,EXPP,XP,YP,ZP,SAAB)
C
C     Sheela Kirpekar, Jan. 93 (based on FRMINT tuh feb 91)
C
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
      PARAMETER (D0 = 0.0D0)
      DOUBLE PRECISION DWNFAC
      DIMENSION SHLINT(KCKTAB)
#include "onecom.h"
#include "lmns.h"
#include "symmet.h"
#include "nuclei.h"
C
#include "ibtfun.h"
#include "pi.h"
C
      DWNFAC = 0.5D0*PI*ALPHA2
      DFAC = SAAB*DWNFAC
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         INT = INT + 1
         DO 300 IATOM = 1, NUCIND
            FRM = D0
            DO 400 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               XC = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
               YC = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
               ZC = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
               RCP2 = (XC-XP)**2 + (YC-YP)**2 + (ZC-ZP)**2
               FINT = EXP(-EXPP*RCP2)
               IF (LVALA .GT. 0) FINT = FINT*((XC - CORAX)**LVALA)
               IF (MVALA .GT. 0) FINT = FINT*((YC - CORAY)**MVALA)
               IF (NVALA .GT. 0) FINT = FINT*((ZC - CORAZ)**NVALA)
               IF (LVALB .GT. 0) FINT = FINT*((XC - CORBX)**LVALB)
               IF (MVALB .GT. 0) FINT = FINT*((YC - CORBY)**MVALB)
               IF (NVALB .GT. 0) FINT = FINT*((ZC - CORBZ)**NVALB)
               FRM = FRM + FINT
            END IF
  400       CONTINUE
            SHLINT(INT) = SHLINT(INT) - FRM*DFAC*CHARGE(IATOM)
  300    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck psoint */
      SUBROUTINE PSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     Paramagnetic spin-orbit integrals
C
C     tuh 8 Feb 1991
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               EE = EU*EV
               FE = FU*EV
               EF = EU*FV
               DO 400 IT = 0,MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  FEE = FT*EE
                  EFE = ET*FE
                  EEF = ET*EF
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IX = 3*(IATOM - 1) + 1
                     IY = 3*(IATOM - 1) + 2
                     IZ = 3*(IATOM - 1) + 3
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
                     SHLINT(INT,IX) = SHLINT(INT,IX)+EFE*AH0V-EEF*AH0U
                     SHLINT(INT,IY) = SHLINT(INT,IY)+EEF*AH0T-FEE*AH0V
                     SHLINT(INT,IZ) = SHLINT(INT,IZ)+FEE*AH0U-EFE*AH0T
                     IOFF = IOFF + NAHGTF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
         INT = INT + 1
  100 CONTINUE
      RETURN
      END
C  /* Deck nsttra */
      SUBROUTINE NSTTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,IORBA,IORBB,NBAST,INTTYP,
     &                  IPRINT,INTADR)
C  MI entry for LAO project
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      DO 100 IREPC = 0, MAXREP
         DO 200 IATOMC = 1, NATOMC
            ICENTC = JCENTC(IATOMC)
            ISYMC  = JSYMC(IATOMC)
            FACTOR = FACINT(IATOMC)
            DO 300 LCOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENTC - 1) + LCOOR,IREPC,2)
               IF (ISCOOR .GT. 0) THEN
                  DO 400 ICOOR = 1, 3
CMI                .... flag for RNST type of integrals
                     IF (INTTYP.EQ.55.OR.INTTYP.EQ.204) THEN
                        IREPCD = IBTXOR(IREPC,ISYMAX(ICOOR,1))
                        ISYMCA = IBTXOR(ISYMAX(ICOOR,1),ISYMAX(LCOOR,2))
                     ELSE
                        IREPCD = IBTXOR(IREPC,ISYMAX(ICOOR,2))
                        ISYMCA = IBTXOR(ISYMAX(ICOOR,2),ISYMAX(LCOOR,2))
                     END IF
                     FACSYM = - FACTOR*PT(IBTAND(ISYMCA,ISYMC))
     &                                *PT(IBTAND(IREPCD,ISYMC))
                     ITYP = 9*(IATOMC - 1) + 3*(LCOOR - 1) + ICOOR
                     IADR = INTADR(3*(ISCOOR - 1) + ICOOR)
CMI ... the RNST integral (204) is truly SQUARE !
                     IF (INTTYP .EQ. 26 .OR. INTTYP .EQ. 27
     &                     .OR. INTTYP .EQ. 204 ) THEN
                        CALL SYMSQR(SHLINT(1,ITYP),SOINT(1,IADR),IREPCD,
     &                              ISYMOP,IORBA,IORBB,FACSYM,NBAST,
     &                              IPRINT)
                     ELSE IF (IREPCD .EQ. 0) THEN
                        CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),ISYMOP,
     &                             MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                             FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                     ELSE
                        CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPCD,
     &                             ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                             KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                             IDUM,IPRINT)
                     END IF
 400              CONTINUE
               END IF
 300        CONTINUE
 200     CONTINUE
 100   CONTINUE
       RETURN
       END
C  /* Deck psotra */
      SUBROUTINE PSOTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      DO 100 IREPC = 0, MAXREP
         DO 200 IATOMC = 1, NATOMC
            ICENTC = JCENTC(IATOMC)
            ISYMC  = JSYMC (IATOMC)
            FACTOR = FACINT(IATOMC)
            DO 300 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENTC - 1) + ICOOR,IREPC,2)
               IF (ISCOOR .GT. 0) THEN
                  ISYMCR = ISYMAX(ICOOR,2)
                  FACSYM = - FACTOR*PT(IBTAND(ISYMCR,ISYMC))
     &                             *PT(IBTAND(IREPC ,ISYMC))
                  ITYP = 3*(IATOMC - 1) + ICOOR
                  IADR = INTADR(ISCOOR)
                  IF (IREPC .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPC,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck hdbtra */
      SUBROUTINE HDBTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,IORBA,IORBB,
     &                  NBAST,INTADR,IPRINT)
#include "implicit.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER B, X
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      ICOMP = 0
      DO 100 X = 1, 3
         DO 200 B = 1, 3
            ICOMP = ICOMP + 1
            DO 300 IREPO = 0, MAXREP
               ISCOOR = IPTCNT(3*(NCENTA - 1) + X,IREPO,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(3*(ISCOOR - 1) + B)
                  IREPCL = IBTXOR(IREPO,ISYMAX(B,2))
                  FACSYM = HKAB*PT(IBTAND(IREPCL,ISYMOP))
                  CALL SYMSQR(SHLINT(1,ICOMP),SOINT(1,ITYP),IREPCL,
     &                        ISYMOP,IORBA,IORBB,FACSYM,NBAST,IPRINT)
               END IF
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      RETURN
      END


C  /* Deck hdotra */
      SUBROUTINE HDOTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,IORBA,IORBB,
     &                  NBAST,INTTYP,IPRINT)
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      DO 100 IREPO = 0, MAXREP
         DO 200 ICOOR = 1, 3
            IF (INTTYP .EQ. 14) THEN
               ISCOOR = IPTCNT(3*(NCENTA - 1) + ICOOR,IREPO,1)
               FAC = HKAB
            ELSE
               ISCOOR = IPTCNT(3*(NCENTB - 1) + ICOOR,IREPO,1)
               FAC = HKAB*PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
     *                   *PT(IBTAND(IREPO, ISYMOP))
            END IF
            IF (ISCOOR .GT. 0) THEN
               CALL SYMSQR(SHLINT(1,ICOOR),SOINT(1,ISCOOR),IREPO,ISYMOP,
     &                     IORBA,IORBB,FAC,NBAST,IPRINT)
            END IF
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck symsqr */
      SUBROUTINE SYMSQR(AO,SO,IREPO,ISYMOP,IORBA,IORBB,FACTOR,NBAST,
     &                  IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION AO(KHKTAB), SO(NBAST,NBAST)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      DO 100 IREPA = 0, MAXREP
         IREPB = IBTXOR(IREPO,IREPA)
         DO 200 NA = 1, KHKTA
         IF (IBTAND(MULA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))).EQ.0) THEN
            NAT = KHKTB*(NA - 1)
            IA  = IPTSYM(IORBA + NA,IREPA)
            DO 300 NB = 1,KHKTB
            IF (IBTAND(MULB,IBTXOR(IREPB,ISYMAO(NHKTB,NB))).EQ.0) THEN
               IB  = IPTSYM(IORBB + NB,IREPB)
               FAC = PT(IBTAND(ISYMOP,IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
               SO(IA,IB) = SO(IA,IB) + FAC*FACTOR*AO(NAT+NB)
            END IF
300         CONTINUE
         END IF
200      CONTINUE
100   CONTINUE
      RETURN
      END
C  /* Deck sdint */
      SUBROUTINE SDINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                 NOPTYP,NATOMC,INTTYP,SAAB,EXPP,CORCX,CORCY,CORCZ,
     &                 CORPX,CORPY,CORPZ)
C
C     tuh Feb 09 91 (g_e factor 930726-hjaaj)
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "pi.h"
#include "gfac.h"
      PARAMETER (D2 = 2.0D0, D3 = 3.0D0)
      PARAMETER (DFAC11 = D2*GFAC*PI/D3, DFAC13 = D2*GFAC*PI)
      PARAMETER (DFACSD = GFAC / D2)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),AHGTF(*),
     &          CORCX(*), CORCY(*), CORCZ(*), SHLINT(KCKTAB,NOPTYP)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         INT = INT + 1
C
C        Subtract Fermi contact contribution
C
         IF (INTTYP .EQ. 11) THEN
            FACTOR = DFAC11
         ELSE
            FACTOR = DFAC13
         END IF
         DO 150 IATOM = 1,NATOMC
            I0 = 6*(IATOM - 1)
            XC = CORCX(IATOM)
            YC = CORCY(IATOM)
            ZC = CORCZ(IATOM)
            DIS = (XC-CORPX)**2 + (YC-CORPY)**2 + (ZC-CORPZ)**2
            DIR = SAAB*FACTOR*EXP(-EXPP*DIS)
            IF (LVALA.GT.0) DIR = DIR*((XC - CORAX)**LVALA)
            IF (MVALA.GT.0) DIR = DIR*((YC - CORAY)**MVALA)
            IF (NVALA.GT.0) DIR = DIR*((ZC - CORAZ)**NVALA)
            IF (LVALB.GT.0) DIR = DIR*((XC - CORBX)**LVALB)
            IF (MVALB.GT.0) DIR = DIR*((YC - CORBY)**MVALB)
            IF (NVALB.GT.0) DIR = DIR*((ZC - CORBZ)**NVALB)
            SHLINT(INT,I0+1) = SHLINT(INT,I0+1) + DIR
            SHLINT(INT,I0+3) = SHLINT(INT,I0+3) + DIR
            SHLINT(INT,I0+6) = SHLINT(INT,I0+6) + DIR
 150     CONTINUE
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = DFACSD*ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EE = ODC(MVALA,MVALB,IU,0,0,2)*EV
               DO 400 IT = 0,MAXT
                  EEE = ODC(LVALA,LVALB,IT,0,0,1)*EE
                  IADR00 = IADRAU + IT
                  IADR0T = IADR00 + 1
                  IADR0U = IADR00 + ISTEPU
                  IADR0V = IADR00 + ISTEPV
                  IADRTT = IADR0T + 1
                  IADRTU = IADR0T + ISTEPU
                  IADRTV = IADR0T + ISTEPV
                  IADRUU = IADR0U + ISTEPU
                  IADRUV = IADR0U + ISTEPV
                  IADRVV = IADR0V + ISTEPV
                  IADD = - NAHGTF
C
C                 ***** LOOP OVER NUCLEI *****
C
                  DO 500 IATOM = 1,NATOMC
                     I0 = 6*(IATOM - 1)
                     IADD = IADD + NAHGTF
                     SHLINT(INT,I0+1) = SHLINT(INT,I0+1)
     &                                + EEE*AHGTF(IADRTT + IADD)
                     SHLINT(INT,I0+2) = SHLINT(INT,I0+2)
     &                                + EEE*AHGTF(IADRTU + IADD)
                     SHLINT(INT,I0+3) = SHLINT(INT,I0+3)
     &                                + EEE*AHGTF(IADRUU + IADD)
                     SHLINT(INT,I0+4) = SHLINT(INT,I0+4)
     &                                + EEE*AHGTF(IADRTV + IADD)
                     SHLINT(INT,I0+5) = SHLINT(INT,I0+5)
     &                                + EEE*AHGTF(IADRUV + IADD)
                     SHLINT(INT,I0+6) = SHLINT(INT,I0+6)
     &                                + EEE*AHGTF(IADRVV + IADD)
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck sdtra */
      SUBROUTINE SDTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                 NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          INTADR(*), FACINT(NATOMC), JSYMC(NATOMC),
     &          JCENTC(NATOMC)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      DO 100 IREPC = 0, MAXREP
         DO 200 IATOMC = 1, NATOMC
            ICENTC =   JCENTC(IATOMC)
            ISYMC  =   JSYMC (IATOMC)
            FACTOR = - FACINT(IATOMC)
            DO 300 ICOOR1 = 1, 3
               ISCOR1 = IPTCNT(3*(ICENTC - 1) + ICOOR1,IREPC,2)
               IF (ISCOR1 .GT. 0) THEN
                  DO 400 ICOOR2 = 1, 3
                     ISCOOR = 3*(ISCOR1 - 1) + ICOOR2
                     IREPCD = IBTXOR(IREPC,ISYMAX(ICOOR2,2))
                     ISYMCR = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
                     FACSYM = FACTOR*PT(IBTAND(ISYMCR,ISYMC))
     &                              *PT(IBTAND(IREPCD,ISYMC))
                     IADR = INTADR(ISCOOR)
                     MXCR = MAX(ICOOR1,ICOOR2)
                     MNCR = MIN(ICOOR1,ICOOR2)
                     ITYP = 6*(IATOMC - 1) + MXCR*(MXCR - 1)/2 + MNCR
                     IF (IREPCD .EQ. 0) THEN
                        CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),
     &                             ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                             KHKTB,FACSYM,LDIAG,FULMAT,DUM,IDUM,
     &                             IPRINT)
                     ELSE
                        CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPCD,
     &                             ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                             KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                             IDUM,IPRINT)
                     END IF
  400             CONTINUE
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck pvctra */
      SUBROUTINE PVCTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,IPVCTBL,IPRINT)
C
C     J.Thyssen - Oct 7 1999
C
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          IPVCTBL(NUCIND,0:MAXREP)
#include "onecom.h"
#include "symmet.h"
#include "nuclei.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ITYP = 0
      DO 100 IATOMC = 1, NATOMC
         ICENTC = JCENTC(IATOMC)
         ISYMC  = JSYMC (IATOMC)
         FACTOR = FACINT(IATOMC)
         ITYP = ITYP + 1
         DO 200 IREPC = 0, MAXREP
            IF (IBTAND(ISTBNU(ICENTC),IREPC) .EQ. 0) THEN
               FACSYM = FACTOR*PT(IBTAND(IREPC,ISYMC))
               IADR = IPVCTBL(ICENTC,IREPC)
               IF (IREPC .EQ. 0) THEN
                  CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),ISYMOP,
     &                       MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                       FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
               ELSE
                  CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPC,
     &                       ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                       FACSYM,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck npetra */
      SUBROUTINE NPETRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,NPETBL,IPRINT)
C
C     K.Ruud, July 92
C
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          NPETBL(NUCIND,0:MAXREP)
#include "onecom.h"
#include "symmet.h"
#include "nuclei.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ITYP = 0
      DO 100 IATOMC = 1, NATOMC
         ICENTC = JCENTC(IATOMC)
         ISYMC  = JSYMC (IATOMC)
         FACTOR = FACINT(IATOMC)
         DO 200 IREPC = 0, MAXREP
            IF (IBTAND(ISTBNU(ICENTC),IREPC) .EQ. 0) THEN
               ITYP = ITYP + 1
               FACSYM = FACTOR*PT(IBTAND(IREPC,ISYMC))
               IADR = NPETBL(ICENTC,IREPC)
               IF (IREPC .EQ. 0) THEN
                  CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),ISYMOP,
     &                       MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                       FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
               ELSE
                  CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPC,
     &                       ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                       FACSYM,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck eftra */
      SUBROUTINE EFTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
C
C  KR, March 92
C     Symmetry properly treated october-95 (sic!)
C
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          INTADR(*)
#include "onecom.h"
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      DO 100 IREPC = 0, MAXREP
         DO 200 IATOMC = 1, NATOMC
            ICENTC = JCENTC(IATOMC)
            ISYMC  = JSYMC (IATOMC)
            FACTOR = FACINT(IATOMC)
            DO 300 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENTC - 1) + ICOOR,IREPC,1)
               IF (ISCOOR .GT. 0) THEN
                  ISYMCR = ISYMAX(ICOOR,1)
                  FACSYM = - FACTOR*PT(IBTAND(ISYMCR,ISYMC))
     &                             *PT(IBTAND(IREPC ,ISYMC))
                  ITYP = 3*(IATOMC - 1) +ICOOR
                  IADR = INTADR(ISCOOR)
                  IF (IREPC .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPC,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck efgtra */
      SUBROUTINE EFGTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,IFGTBL,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          IFGTBL(NUCIND,6,0:MAXREP)
#include "onecom.h"
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
C
      FULMAT = .TRUE.
      ITYP = 0
      DO 100 IATOMC = 1, NATOMC
         ICENTC = JCENTC(IATOMC)
         ISYMC  = JSYMC(IATOMC)
         FACTOR = FACINT(IATOMC)
         IJ = 0
         DO 200 ICOOR1 = 1, 3
         DO 200 ICOOR2 = ICOOR1, 3
            IJ = IJ + 1
            ITYP = ITYP + 1
            ISYMXY = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            DO 300 IREPC = 0, MAXREP
              IF (IBTAND(ISTBNU(ICENTC),IBTXOR(IREPC,ISYMXY)).EQ.0) THEN
                  FACSYM = FACTOR*PT(IBTAND(ISYMXY,ISYMC))
     &                           *PT(IBTAND(IREPC ,ISYMC))
                  IADR = IFGTBL(ICENTC,IJ,IREPC)
                  IF (IREPC .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPC,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
 300        CONTINUE
 200     CONTINUE
 100   CONTINUE
       RETURN
       END
C  /* Deck efttra */
      SUBROUTINE EFTTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,IFGTBL,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          IFGTBL(NUCIND,15,0:MAXREP)
#include "onecom.h"
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
C
      FULMAT = .TRUE.
      ITYP = 0
      DO 100 IATOMC = 1, NATOMC
         ICENTC = JCENTC(IATOMC)
         ISYMC  = JSYMC(IATOMC)
         FACTOR = FACINT(IATOMC)
         IJ = 0
         DO 200 ICOOR1 = 1, 3
         DO 200 ICOOR2 = ICOOR1, 3
         DO 200 ICOOR3 = ICOOR2, 3
         DO 200 ICOOR4 = ICOOR3, 3
            IJ = IJ + 1
            ITYP = ITYP + 1
            ISYMXY = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            ISYMXY = IBTXOR(ISYMXY          ,ISYMAX(ICOOR3,1))
            ISYMXY = IBTXOR(ISYMXY          ,ISYMAX(ICOOR4,1))
            DO 300 IREPC = 0, MAXREP
              IF (IBTAND(ISTBNU(ICENTC),IBTXOR(IREPC,ISYMXY)).EQ.0) THEN
                  FACSYM = FACTOR*PT(IBTAND(ISYMXY,ISYMC))
     &                           *PT(IBTAND(IREPC ,ISYMC))
                  IADR = IFGTBL(ICENTC,IJ,IREPC)
                  IF (IREPC .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPC,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
 300        CONTINUE
 200     CONTINUE
 100   CONTINUE
       RETURN
       END
C  /* Deck dsotra */
      SUBROUTINE DSOTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,DOATOM,
     &                  KAB,TRIANG,NATOM,INTADR,NSHINT,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI, DOATOM(NUCIND), SAMECD, TRIANG, TRICR
      DIMENSION SHLINT(KCKTAB,NSHINT), SOINT(NELMNT,NOPTYP),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "nuclei.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
C
      DO 100 IREPO = 0, MAXREP
        ISTAO = 0
C
C       Atom D
C
        DO 110 IATOMD = 1, NUCIND
        IF (DOATOM(IATOMD)) THEN
          ISTABD = ISTBNU(IATOMD)
          DO 120 ISYMD = 0, MAXOPR
          IF (IBTAND(ISYMD,IBTOR(ISTABD,KAB)) .EQ. 0) THEN
            KABD = IBTAND(ISTABD,KAB)
C
C           Atom C
C
C           (First determine NATOMC)
C
            NATOMC = 0
            MXATMC = NUCIND
            IF (TRIANG) MXATMC = IATOMD
            DO 130 IATOMC = 1, MXATMC
            IF (DOATOM(IATOMC)) THEN
               SAMECD = IATOMC .EQ. IATOMD
               MABCD  = IBTOR(ISTBNU(IATOMC),KABD)
               DO 135 ISYMC = 0, MAXOPR
                  IF (SAMECD .AND. ISYMC.EQ.ISYMD) GO TO 135
                  IF (IBTAND(ISYMC,MABCD) .EQ. 0) NATOMC = NATOMC + 1
  135          CONTINUE
            END IF
  130       CONTINUE
            IF (NATOMC .GT. 0) THEN
              ITYP0 = ISTAO
              MXATMC = NUCIND
              IF (TRIANG) MXATMC = IATOMD
              DO 200 IATOMC = 1, MXATMC
              IF (DOATOM(IATOMC)) THEN
                SAMECD = IATOMC .EQ. IATOMD
                TRICR = TRIANG .AND. SAMECD
                DO 210 ISYMC = 0, MAXOPR
                IF (SAMECD .AND. ISYMC.EQ.ISYMD) GO TO 210
                IF (IBTAND(ISYMC,IBTOR(ISTBNU(IATOMC),KABD)).EQ.0) THEN
C
C                 Cartesian directions
C
                  DO 300 ICOORD = 1, 3
                    ISCORD = IPTCNT(3*(IATOMD - 1) + ICOORD,IREPO,2)
                    IF (ISCORD .GT. 0) THEN
                      FACD = - PT(IBTAND(ISYMAX(ICOORD,2),ISYMD))
     &                        *PT(IBTAND(IREPO ,ISYMD))
                      MXCRC = 3
                      IF (TRICR) MXCRC = ICOORD
                      DO 310 ICOORC = 1, MXCRC
                        ISCORC = IPTCNT(3*(IATOMC - 1) + ICOORC,IREPO,2)
                        IF (ISCORC .GT. 0) THEN
                          IF (TRIANG) THEN
                             MXCOR = MAX(ISCORC,ISCORD)
                             MNCOR = MIN(ISCORC,ISCORD)
                             ISCOOR = MXCOR*(MXCOR - 1)/2 + MNCOR
                          ELSE
                             ISCOOR = 3*NUCDEP*(ISCORD - 1) + ISCORC
                          END IF
                          IADR = INTADR(ISCOOR)
                          FACSYM=FACD*PT(IBTAND(ISYMAX(ICOORC,2),ISYMC))
     &                               *PT(IBTAND(IREPO ,ISYMC))
                          ITYP = ITYP0 + 3*(ICOORC - 1) + ICOORD
                          CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),
     &                               ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &                               KHKTA,KHKTB,FACSYM,LDIAG,FULMAT,
     &                               DUM,IDUM,IPRINT)
                        END IF
  310                 CONTINUE
                    END IF
  300             CONTINUE
                  ITYP0 = ITYP0 + 9
                END IF
  210           CONTINUE
              END IF
  200         CONTINUE
              ISTAO = ISTAO + 9*NATOMC
            END IF
          END IF
  120     CONTINUE
        END IF
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck dsoodc */
      SUBROUTINE DSOODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,EXPPI,DIFPDX,
     &                  DIFPDY,DIFPDZ,IPRINT)
C
C     TUH Feb 19 1991
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.D0, D2 = 2.D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
      EXPPIH = EXPPI/D2
      DO 100 IA = 0, JMAXA
         DO 200 IB = 0, JMAXB
            DO 300 I = 0, IA + IB + 1
               FX = D0
               FY = D0
               FZ = D0
               IF (I .LE. IA + IB) THEN
                  FX = FX + DIFPDX*ODC(IA,IB,I,0,0,1)
                  FY = FY + DIFPDY*ODC(IA,IB,I,0,0,2)
                  FZ = FZ + DIFPDZ*ODC(IA,IB,I,0,0,3)
               END IF
               IF (I .GT. 0) THEN
                  FX = FX + EXPPIH*ODC(IA,IB,I - 1,0,0,1)
                  FY = FY + EXPPIH*ODC(IA,IB,I - 1,0,0,2)
                  FZ = FZ + EXPPIH*ODC(IA,IB,I - 1,0,0,3)
               END IF
               IF (I .LT. IA + IB) THEN
                  FX = FX + dble(I + 1)*ODC(IA,IB,I + 1,0,0,1)
                  FY = FY + dble(I + 1)*ODC(IA,IB,I + 1,0,0,2)
                  FZ = FZ + dble(I + 1)*ODC(IA,IB,I + 1,0,0,3)
               END IF
               ODC(IA,IB,I,1,0,1) = FX
               ODC(IA,IB,I,1,0,2) = FY
               ODC(IA,IB,I,1,0,3) = FZ
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 20) THEN
         JMAXAB = (JMAXA + 1)*(JMAXB + 1)
         ISTP   = JMAXA + JMAXB + 2
         CALL TITLER('Output from DSOODC','*',103)
         CALL AROUND('x component')
         CALL OUTPUT(ODC(0,0,0,1,0,1),1,JMAXAB,1,ISTP,JMAXAB,ISTP,1,
     &               LUPRI)
         CALL AROUND('y component')
         CALL OUTPUT(ODC(0,0,0,1,0,2),1,JMAXAB,1,ISTP,JMAXAB,ISTP,1,
     &               LUPRI)
         CALL AROUND('z component')
         CALL OUTPUT(ODC(0,0,0,1,0,3),1,JMAXAB,1,ISTP,JMAXAB,ISTP,1,
     &               LUPRI)
      END IF
      RETURN
      END
C  /* Deck dsoint */
      SUBROUTINE DSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     tuh & ov Feb 09 91
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),AHGTF(*),
     &          SHLINT(KCKTAB,*)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         INT = INT + 1
         IADRAV = 1
         DO 200 IV = 0, MAXV + 1
            IADRAU = IADRAV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            IUMAX = MAXU + 1
            IF (IV .GT. MAXV) IUMAX = MAXU
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               ITMAX = MAXT + 1
               IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) ITMAX = MAXT
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  FX = FT*EU*EV
                  FY = ET*FU*EV
                  FZ = ET*EU*FV
                  IADR0 = IADRAU + IT
                  IADRT = IADR0  + 1
                  IADRU = IADR0  + ISTEPU
                  IADRV = IADR0  + ISTEPV
C
C                 ***** LOOP OVER NUCLEI *****
C
                  DO 500 IATOM = 1, NATOMC
                     I0 = 9*(IATOM - 1)
                     SHLINT(INT,I0+1) = SHLINT(INT,I0+1)+FY*AHGTF(IADRU)
     &                                                  +FZ*AHGTF(IADRV)
                     SHLINT(INT,I0+2) = SHLINT(INT,I0+2)-FX*AHGTF(IADRU)
                     SHLINT(INT,I0+3) = SHLINT(INT,I0+3)-FX*AHGTF(IADRV)
                     SHLINT(INT,I0+4) = SHLINT(INT,I0+4)-FY*AHGTF(IADRT)
                     SHLINT(INT,I0+5) = SHLINT(INT,I0+5)+FX*AHGTF(IADRT)
     &                                                  +FZ*AHGTF(IADRV)
                     SHLINT(INT,I0+6) = SHLINT(INT,I0+6)-FY*AHGTF(IADRV)
                     SHLINT(INT,I0+7) = SHLINT(INT,I0+7)-FZ*AHGTF(IADRT)
                     SHLINT(INT,I0+8) = SHLINT(INT,I0+8)-FZ*AHGTF(IADRU)
                     SHLINT(INT,I0+9) = SHLINT(INT,I0+9)+FX*AHGTF(IADRT)
     &                                                  +FY*AHGTF(IADRU)
C                    SHLINT(INT,I0+1) = SHLINT(INT,I0+1)-FX*AHGTF(IADRT)
C                    SHLINT(INT,I0+5) = SHLINT(INT,I0+5)-FY*AHGTF(IADRU)
C                    SHLINT(INT,I0+9) = SHLINT(INT,I0+9)-FZ*AHGTF(IADRV)
                     IADRT = IADRT + NAHGTF
                     IADRU = IADRU + NAHGTF
                     IADRV = IADRV + NAHGTF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck symupk */
      SUBROUTINE SYMUPK(AINT,WORK,ISYMO,NELMNT)
C
C  290689 Henrik Koch
C
C  Purpose:
C          Symmetry unpack the matrix AINT into WORK
C
C
#include "implicit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
      DIMENSION AINT(NELMNT),WORK(NELMNT),IOFFBL(8)
C
#include "symmet.h"
C
#include "ibtfun.h"
C
C     Initialize work array.
C
      CALL DZERO(WORK,NELMNT)
C
C     Setup offset array to the different subblocks.
C
      IOFF = 0
      DO 50 I = 1,MAXREP + 1
         IOFFBL(I) = IOFF
         IOFF = IOFF + NAOS(I)
   50 CONTINUE
C
      IF (ISYMO .EQ. 1) THEN
         IOFF2 = 0
         DO 100 ISYMA = 1,MAXREP + 1
            NDIMA = NAOS(ISYMA)
            IOFF1 = IOFFBL(ISYMA)
            DO 110 J = 1,NDIMA
               DO 120 K = J,NDIMA
                  KT  = IOFF1 + K
                  JT  = IOFF1 + J
                  KJT = KT*(KT-1)/2 + JT
                  WORK(KJT) = AINT(IOFF2 + K*(K-1)/2 + J)
  120          CONTINUE
  110       CONTINUE
            IOFF2 = IOFF2 + NDIMA*(NDIMA + 1)/2
  100    CONTINUE
      ELSE
C tsaue 971112 - Bugfix here for SP2...
C ISYMA = IBTXOR(ISYMO - 1,ISYMB - 1) + 1 fungerer ikke...
         IREPO = ISYMO - 1     
         DO 200 IREPB = 0, MAXREP
            ISYMB = IREPB + 1
            ISYMA = IBTXOR(IREPO,IREPB) + 1
            IF (ISYMA .GT. ISYMB) THEN
               NDIMA = NAOS(ISYMA)
               NDIMB = NAOS(ISYMB)
               IF (NDIMA.GT.0 .AND. NDIMB.GT.0) THEN
                  IOFF = NPARNU(ISYMO,ISYMA)
                  DO 210 I = 1,NDIMB
                     DO 220 J = 1,NDIMA
                        JT  = IOFFBL(ISYMA) + J
                        JIT = JT*(JT - 1)/2 + IOFFBL(ISYMB) + I
                        WORK(JIT) = AINT(IOFF + I + (J - 1)*NDIMB)
  220                CONTINUE
  210             CONTINUE
               ENDIF
            ENDIF
  200    CONTINUE
      ENDIF
      RETURN
      END
C  /* Deck gosint */
      SUBROUTINE GOSINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,ORIGIN,EXPPI)
C
C     tuh 1993
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP25 = 0.25D0, D1 = 1.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6), ORIGIN(3)
      INTEGER T, U, V
#include "onecom.h"
#include "lmns.h"
      DKX = ORIGIN(1)
      DKY = ORIGIN(2)
      DKZ = ORIGIN(3)
      FACX = SHGTF/EXP(DP25*EXPPI*DKX**2)
      FACY = SHGTF/EXP(DP25*EXPPI*DKY**2)
      FACZ = SHGTF/EXP(DP25*EXPPI*DKZ**2)
      COSX = COS(CORPX*DKX)*FACX
      COSY = COS(CORPY*DKY)*FACY
      COSZ = COS(CORPZ*DKZ)*FACZ
      SINX = SIN(CORPX*DKX)*FACX
      SINY = SIN(CORPY*DKY)*FACY
      SINZ = SIN(CORPZ*DKZ)*FACZ
C
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         CX    = D0
         SSX   = D0
         SGNC  = D1
         SGNS  = D1
         FACKX = D1
         DO 200 T = 0, LVALA + LVALB
            IF (MOD(T,2) .EQ. 0) THEN
               CX  = CX + SGNC*FACKX*COSX*ODC(LVALA,LVALB,T,0,0,1)
               SSX = SSX + SGNS*FACKX*SINX*ODC(LVALA,LVALB,T,0,0,1)
               SGNC = - SGNC
            ELSE
               CX =   CX + SGNC*FACKX*SINX*ODC(LVALA,LVALB,T,0,0,1)
               SSX = SSX + SGNS*FACKX*COSX*ODC(LVALA,LVALB,T,0,0,1)
               SGNS = - SGNS
            END IF
            FACKX = FACKX*DKX
  200    CONTINUE
C
         CY    = D0
         SY    = D0
         SGNC  = D1
         SGNS  = D1
         FACKY = D1
         DO 300 U = 0, MVALA + MVALB
            IF (MOD(U,2) .EQ. 0) THEN
               CY = CY + SGNC*FACKY*COSY*ODC(MVALA,MVALB,U,0,0,1)
               SY = SY + SGNS*FACKY*SINY*ODC(MVALA,MVALB,U,0,0,1)
               SGNC = - SGNC
            ELSE
               CY = CY + SGNC*FACKY*SINY*ODC(MVALA,MVALB,U,0,0,1)
               SY = SY + SGNS*FACKY*COSY*ODC(MVALA,MVALB,U,0,0,1)
               SGNS = - SGNS
            END IF
            FACKY = FACKY*DKY
  300    CONTINUE
C
         CZ    = D0
         SZ    = D0
         SGNC  = D1
         SGNS  = D1
         FACKZ = D1
         DO 400 V = 0, NVALA + NVALB
            IF (MOD(V,2) .EQ. 0) THEN
               CZ = CZ + SGNC*FACKZ*COSZ*ODC(NVALA,NVALB,V,0,0,1)
               SZ = SZ + SGNS*FACKZ*SINZ*ODC(NVALA,NVALB,V,0,0,1)
               SGNC = - SGNC
            ELSE
               CZ = CZ + SGNC*FACKZ*SINZ*ODC(NVALA,NVALB,V,0,0,1)
               SZ = SZ + SGNS*FACKZ*COSZ*ODC(NVALA,NVALB,V,0,0,1)
               SGNS = - SGNS
            END IF
            FACKZ = FACKZ*DKZ
  400    CONTINUE
C
         OX = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         OY = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         OZ = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
C
         INT = INT + 1
         SHLINT(INT,1) = SHLINT(INT,1) -  CX*OY*OZ
         SHLINT(INT,2) = SHLINT(INT,2) -  OX*CY*OZ
         SHLINT(INT,3) = SHLINT(INT,3) -  OX*OY*CZ
         SHLINT(INT,4) = SHLINT(INT,4) - SSX*OY*OZ
         SHLINT(INT,5) = SHLINT(INT,5) -  OX*SY*OZ
         SHLINT(INT,6) = SHLINT(INT,6) -  OX*OY*SZ
  100 CONTINUE
      RETURN
      END
C/* Deck cexpint */
C***********************************************************************      
      SUBROUTINE CEXPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &     CORPX,CORPY,CORPZ,EXPPI,NINTS,WORK,LWORK)
C***********************************************************************      
C
C     nhl 2016 -- calculates integrals over exp(ikr)
C                 based on GOSINT routine.
C     References are:
C     Nanna Holmgaard List, Joanna Kauczor, Trond Saue,
C     Hans Jørgen Aagaard Jensen and Patrick Norman, J. Chem. Phys. 142 (2015) 244111
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "onecom.h"             !contains KCKTAB
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &     SHLINT(KCKTAB,NINTS),WORK(LWORK)
#include "memint.h"
      CALL MEMGET2('REAL','TMPARR',KTMP,8*KCKTAB,WORK,KFREE,LFREE)
      CALL CEXPIN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &     CORPX,CORPY,CORPZ,EXPPI,NINTS,WORK(KTMP))
      CALL MEMREL('CEXPINT',WORK,KWORK,KTMP,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      SUBROUTINE CEXPIN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &     CORPX,CORPY,CORPZ,EXPPI,NINTS,TMPARR)
C***********************************************************************      
C
C     nhl 2016 -- calculates integrals over exp(ikr)
C                 based on GOSINT routine.
C     References are:
C     Nanna Holmgaard List, Joanna Kauczor, Trond Saue,
C     Hans Jørgen Aagaard Jensen and Patrick Norman, J. Chem. Phys. 142 (2015) 244111
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP25 = 0.25D0, D1 = 1.0D0, NIREP = 8)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &     SHLINT(KCKTAB,NINTS),TMPARR(KCKTAB,8)
      INTEGER T, U, V, NINTS
#include "onecom.h" !contains KCKTAB
#include "lmns.h"
#include "symmet.h"
#include "inforb.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbwav.h"      
C
#include "ibtfun.h"
C
      CALL DZERO(TMPARR,KCKTAB*8)
C     
      DKX  = WAVEVEC(1)
      DKY  = WAVEVEC(2)
      DKZ  = WAVEVEC(3)
      FACX = SHGTF/EXP(DP25*EXPPI*DKX**2)
      FACY = SHGTF/EXP(DP25*EXPPI*DKY**2)
      FACZ = SHGTF/EXP(DP25*EXPPI*DKZ**2)
      COSX = COS(CORPX*DKX)*FACX
      COSY = COS(CORPY*DKY)*FACY
      COSZ = COS(CORPZ*DKZ)*FACZ
      SINX = SIN(CORPX*DKX)*FACX
      SINY = SIN(CORPY*DKY)*FACY
      SINZ = SIN(CORPZ*DKZ)*FACZ
C
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA) !power of x factor in basis a
         MVALA = MVALUA(ICOMPA) !power of y factor in basis a
         NVALA = NVALUA(ICOMPA) !power of z factor in basis a
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB) !power of x factor in basis b
         MVALB = MVALUB(ICOMPB) !power of y factor in basis b
         NVALB = NVALUB(ICOMPB) !power of z factor in basis b
C
         CX    = D0
         CXDX  = D0
         SX    = D0
         SGNC  = D1
         SGNS  = D1
         FACKX = D1
         DO 200 T = 0, LVALA + LVALB
            IF (MOD(T,2) .EQ. 0) THEN
               CX = CX + SGNC*FACKX*COSX*ODC(LVALA,LVALB,T,0,0,1)
               SX = SX + SGNS*FACKX*SINX*ODC(LVALA,LVALB,T,0,0,1)
               SGNC = - SGNC
            ELSE
               CX = CX + SGNC*FACKX*SINX*ODC(LVALA,LVALB,T,0,0,1)
               SX = SX + SGNS*FACKX*COSX*ODC(LVALA,LVALB,T,0,0,1)
               SGNS = - SGNS
            END IF
            FACKX = FACKX*DKX
  200    CONTINUE
         CY    = D0
         CYDY  = D0
         SY    = D0
         SGNC  = D1
         SGNS  = D1
         FACKY = D1
         DO 300 U = 0, MVALA + MVALB
            IF (MOD(U,2) .EQ. 0) THEN
               CY = CY + SGNC*FACKY*COSY*ODC(MVALA,MVALB,U,0,0,2)
               SY = SY + SGNS*FACKY*SINY*ODC(MVALA,MVALB,U,0,0,2)
               SGNC = - SGNC
            ELSE
               CY = CY + SGNC*FACKY*SINY*ODC(MVALA,MVALB,U,0,0,2)
               SY = SY + SGNS*FACKY*COSY*ODC(MVALA,MVALB,U,0,0,2)
               SGNS = - SGNS
            END IF
            FACKY = FACKY*DKY
  300    CONTINUE
C
         CZ    = D0
         CZDZ  = D0
         SZ    = D0
         SGNC  = D1
         SGNS  = D1
         FACKZ = D1
         DO 400 V = 0, NVALA + NVALB
            IF (MOD(V,2) .EQ. 0) THEN
               CZ = CZ + SGNC*FACKZ*COSZ*ODC(NVALA,NVALB,V,0,0,3)
               SZ = SZ + SGNS*FACKZ*SINZ*ODC(NVALA,NVALB,V,0,0,3)
               SGNC = - SGNC
            ELSE
               CZ = CZ + SGNC*FACKZ*SINZ*ODC(NVALA,NVALB,V,0,0,3)
               SZ = SZ + SGNS*FACKZ*COSZ*ODC(NVALA,NVALB,V,0,0,3)
               SGNS = - SGNS
            END IF
            FACKZ = FACKZ*DKZ
  400    CONTINUE
C
        INT = INT + 1
CC     Real part: eq. 44
        TMPARR(INT,1) = TMPARR(INT,1) + CX*CY*CZ ! Sym 0
        TMPARR(INT,2) = TMPARR(INT,2) - CX*SY*SZ ! Sym yz
        TMPARR(INT,3) = TMPARR(INT,3) - SX*CY*SZ ! Sym xz
        TMPARR(INT,4) = TMPARR(INT,4) - SX*SY*CZ ! Sym xy
CC     Imaginary part: eq. 45
        TMPARR(INT,5) = TMPARR(INT,5) + CX*CY*SZ ! Sym z
        TMPARR(INT,6) = TMPARR(INT,6) + CX*SY*CZ ! Sym y
        TMPARR(INT,7) = TMPARR(INT,7) + SX*CY*CZ ! Sym x
        TMPARR(INT,8) = TMPARR(INT,8) - SX*SY*SZ ! Sym xyz
 100  CONTINUE
      IOFF = NINTS/2
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
C     accumulating real part on symmetries
      II = IBTOSF(0,1)
      SHLINT(:,II) = SHLINT(:,II) + TMPARR(:,1) ! Sym 0
      II = IBTOSF(ISYMAX(1,2),1)
      SHLINT(:,II) = SHLINT(:,II) + TMPARR(:,2)  ! Sym yz
      II = IBTOSF(ISYMAX(2,2),1)
      SHLINT(:,II) = SHLINT(:,II) + TMPARR(:,3)  ! Sym xz
      II = IBTOSF(ISYMAX(3,2),1)
      SHLINT(:,II) = SHLINT(:,II) + TMPARR(:,4)  ! Sym xy
C     accumulating imaginary part on symmetries
      II = IBTOSF(ISYMAX(3,1),2) + IOFF
      SHLINT(:,II) = SHLINT(:,II) + TMPARR(:,5)  ! Sym z
      II = IBTOSF(ISYMAX(2,1),2) + IOFF
      SHLINT(:,II) = SHLINT(:,II) + TMPARR(:,6)  ! Sym y
      II = IBTOSF(ISYMAX(1,1),2) + IOFF
      SHLINT(:,II) = SHLINT(:,II) + TMPARR(:,7)  ! Sym x
      II = IBTOSF(IXYZ,2) + IOFF
      SHLINT(:,II) = SHLINT(:,II) + TMPARR(:,8)  ! Sym xyz
      RETURN
      END

C  /* Deck sl1int */
      SUBROUTINE SL1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,WORK,LWORK)
C
C     K.Ruud, 1994
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      DIMENSION WORK(LWORK)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3*(IORDER+1)*(IORDER+2)/2)
#include "onecom.h"
C
      KDX    = 1
      KDY    = KDX    + IORDER + 1
      KDZ    = KDY    + IORDER + 1
      KDX1   = KDZ    + IORDER + 1
      KDY1   = KDX1   + IORDER + 1
      KDZ1   = KDY1   + IORDER + 1
      KLO    = KDZ1   + IORDER + 1
      KMO    = KLO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KNO    = KMO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KLAST  = KNO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('SL1INT',' ',KLAST,LWORK)
      CALL SL1IN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,CORPX,
     &            CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,IORDER,
     &            WORK(KDX),WORK(KDY),WORK(KDZ),WORK(KDX1),
     &            WORK(KDY1),WORK(KDZ1),WORK(KLO),WORK(KMO),WORK(KNO))
      RETURN
      END
C  /* Deck sl1in1 */
      SUBROUTINE SL1IN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,DX,DY,DZ,DX1,DY1,DZ1,LO,
     &                  MO,NO)
C
C     K.Ruud-94
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3*(IORDER+1)*(IORDER+2)/2),
     &          DX(0:IORDER), DY(0:IORDER), DZ(0:IORDER),
     &          DX1(0:IORDER), DY1(0:IORDER), DZ1(0:IORDER),
     &          LO((IORDER + 1)*(IORDER + 2)/2),
     &          MO((IORDER + 1)*(IORDER + 2)/2),
     &          NO((IORDER + 1)*(IORDER + 2)/2)
#include "onecom.h"
#include "lmns.h"
#include "orgcom.h"
C
C
      INT = 0
      ADX = CORAX - ORIGIN(1)
      ADY = CORAY - ORIGIN(2)
      ADZ = CORAZ - ORIGIN(3)
C
C     Cartesian integrals
C
      DO 300 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 300 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C        One-dimensional integrals
C
         DO 400 IO = 0, IORDER
            DX(IO)  = SHGTF*ODC(LVALA,LVALB,0,0,IO,1)
            DY(IO)  = SHGTF*ODC(MVALA,MVALB,0,0,IO,2)
            DZ(IO)  = SHGTF*ODC(NVALA,NVALB,0,0,IO,3)
C
            DX1(IO) = SHGTF*(ODC(LVALA + 1,LVALB,0,0,IO,1) +
     &                       ODC(LVALA    ,LVALB,0,0,IO,1)*ADX)
            DY1(IO) = SHGTF*(ODC(MVALA + 1,MVALB,0,0,IO,2) +
     &                       ODC(MVALA    ,MVALB,0,0,IO,2)*ADY)
            DZ1(IO) = SHGTF*(ODC(NVALA + 1,NVALB,0,0,IO,3) +
     &                       ODC(NVALA    ,NVALB,0,0,IO,3)*ADZ)
 400     CONTINUE
C
C        Three-dimensional integrals
C
         INT = INT + 1
         CALL LMNVAL(IORDER+1,(IORDER + 1)*(IORDER + 2)/2,LO,MO,NO)
         DO 600 I = 1, (IORDER + 1)*(IORDER + 2)/2
            IX = 3*(I - 1) + 1
            IY = 3*(I - 1) + 2
            IZ = 3*(I - 1) + 3
            SHLINT(INT,IX) = SHLINT(INT,IX) - DP5*DX(LO(I))
     &                     *(DIFABY*DY(MO(I))*DZ1(NO(I)) -
     &                       DIFABZ*DZ(NO(I))*DY1(MO(I)))
            SHLINT(INT,IY) = SHLINT(INT,IY) - DP5*DY(MO(I))
     &                     *(DIFABZ*DZ(NO(I))*DX1(LO(I)) -
     &                       DIFABX*DX(LO(I))*DZ1(NO(I)))
            SHLINT(INT,IZ) = SHLINT(INT,IZ) - DP5*DZ(NO(I))
     &                     *(DIFABX*DX(LO(I))*DY1(MO(I)) -
     &                       DIFABY*DY(MO(I))*DX1(LO(I)))
  600    CONTINUE
  300 CONTINUE
      RETURN
      END
C  /* Deck sl2int */
      SUBROUTINE SL2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,WORK,LWORK)
C
C     K.Ruud, 1994
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      DIMENSION WORK(LWORK)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6*(IORDER+1)*(IORDER+2)/2)
#include "onecom.h"
C
      KDX    = 1
      KDY    = KDX    + IORDER + 1
      KDZ    = KDY    + IORDER + 1
      KDX1   = KDZ    + IORDER + 1
      KDY1   = KDX1   + IORDER + 1
      KDZ1   = KDY1   + IORDER + 1
      KDX2   = KDZ1   + IORDER + 1
      KDY2   = KDX2   + IORDER + 1
      KDZ2   = KDY2   + IORDER + 1
      KLO    = KDZ2   + IORDER + 1
      KMO    = KLO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KNO    = KMO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KLAST  = KNO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('SL2INT',' ',KLAST,LWORK)
      CALL SL2IN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,CORPX,
     &            CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,IORDER,
     &            WORK(KDX),WORK(KDY),WORK(KDZ),WORK(KDX1),
     &            WORK(KDY1),WORK(KDZ1),WORK(KDX2),WORK(KDY2),
     &            WORK(KDZ2),WORK(KLO),WORK(KMO),WORK(KNO))
      RETURN
      END
C  /* Deck sl2in1 */
      SUBROUTINE SL2IN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,DX,DY,DZ,DX1,DY1,DZ1,DX2,
     &                  DY2,DZ2,LO,MO,NO)
C
C     K.Ruud-94
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP25 = 0.25D0, D2=2.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6*(IORDER+1)*(IORDER+2)/2),
     &          DX(0:IORDER), DY(0:IORDER), DZ(0:IORDER),
     &          DX1(0:IORDER), DY1(0:IORDER), DZ1(0:IORDER),
     &          DX2(0:IORDER), DY2(0:IORDER), DZ2(0:IORDER),
     &          LO((IORDER + 1)*(IORDER + 2)/2),
     &          MO((IORDER + 1)*(IORDER + 2)/2),
     &          NO((IORDER + 1)*(IORDER + 2)/2)
#include "onecom.h"
#include "lmns.h"
#include "orgcom.h"
C
C
      INT = 0
      ADX = CORAX - ORIGIN(1)
      ADY = CORAY - ORIGIN(2)
      ADZ = CORAZ - ORIGIN(3)
C
C     Cartesian integrals
C
      DO 300 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 300 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C        One-dimensional integrals
C
         DO 400 IO = 0, IORDER
            DX(IO)  = SHGTF*ODC(LVALA,LVALB,0,0,IO,1)
            DY(IO)  = SHGTF*ODC(MVALA,MVALB,0,0,IO,2)
            DZ(IO)  = SHGTF*ODC(NVALA,NVALB,0,0,IO,3)
            DX1(IO) = SHGTF*(ODC(LVALA + 1,LVALB,0,0,IO,1) +
     &                   ADX*ODC(LVALA    ,LVALB,0,0,IO,1))
            DX2(IO) = SHGTF*(ODC(LVALA + 2,LVALB,0,0,IO,1) +
     &                D2*ADX*ODC(LVALA + 1,LVALB,0,0,IO,1) +
     &               ADX*ADX*ODC(LVALA    ,LVALB,0,0,IO,1))
            DY1(IO) = SHGTF*(ODC(MVALA + 1,MVALB,0,0,IO,2) +
     &                   ADY*ODC(MVALA    ,MVALB,0,0,IO,2))
            DY2(IO) = SHGTF*(ODC(MVALA + 2,MVALB,0,0,IO,2) +
     &                D2*ADY*ODC(MVALA + 1,MVALB,0,0,IO,2) +
     &               ADY*ADY*ODC(MVALA    ,MVALB,0,0,IO,2))
            DZ1(IO) = SHGTF*(ODC(NVALA + 1,NVALB,0,0,IO,3) +
     &                   ADZ*ODC(NVALA    ,NVALB,0,0,IO,3))
            DZ2(IO) = SHGTF*(ODC(NVALA + 2,NVALB,0,0,IO,3) +
     &                D2*ADZ*ODC(NVALA + 1,NVALB,0,0,IO,3) +
     &               ADZ*ADZ*ODC(NVALA    ,NVALB,0,0,IO,3))
  400    CONTINUE
C
C        Three-dimensional integrals
C
         INT = INT + 1
         CALL LMNVAL(IORDER+1,(IORDER + 1)*(IORDER + 2)/2,LO,MO,NO)
         DO 600 I = 1, (IORDER + 1)*(IORDER + 2)/2
            IXX = 6*(I - 1) + 1
            IXY = 6*(I - 1) + 2
            IXZ = 6*(I - 1) + 3
            IYY = 6*(I - 1) + 4
            IYZ = 6*(I - 1) + 5
            IZZ = 6*(I - 1) + 6
            SX0 = DX(LO(I))
            SY0 = DY(MO(I))
            SZ0 = DZ(NO(I))
            SX1 = DX1(LO(I))
            SY1 = DY1(MO(I))
            SZ1 = DZ1(NO(I))
            SX2 = DX2(LO(I))
            SY2 = DY2(MO(I))
            SZ2 = DZ2(NO(I))
            SHLINT(INT,IXX) = SHLINT(INT,IXX)
     &                      - (D2*DIFABZ*DIFABY*SY1*SZ1*SX0
     &                      - DIFABZ*DIFABZ*SY2*SX0*SZ0
     &                      - DIFABY*DIFABY*SZ2*SX0*SY0)*DP25
            SHLINT(INT,IXY) = SHLINT(INT,IXY)
     &                      - (DIFABZ*DIFABZ*SX1*SY1*SZ0
     &                      - DIFABY*DIFABZ*SZ1*SX1*SY0
     &                      - DIFABZ*DIFABX*SY1*SZ1*SX0
     &                      + DIFABY*DIFABX*SZ2*SX0*SY0)*DP25
            SHLINT(INT,IXZ) = SHLINT(INT,IXZ)
     &                      - (DIFABY*DIFABY*SX1*SZ1*SY0
     &                      - DIFABX*DIFABY*SY1*SZ1*SX0
     &                      + DIFABZ*DIFABX*SY2*SX0*SZ0
     &                      - DIFABY*DIFABZ*SY1*SX1*SZ0)*DP25
            SHLINT(INT,IYY) = SHLINT(INT,IYY)
     &                      - (D2*DIFABZ*DIFABX*SX1*SZ1*SY0
     &                      - DIFABZ*DIFABZ*SX2*SY0*SZ0
     &                      - DIFABX*DIFABX*SZ2*SX0*SY0)*DP25
            SHLINT(INT,IYZ) = SHLINT(INT,IYZ)
     &                      - (DIFABZ*DIFABY*SX2*SY0*SZ0
     &                      - DIFABZ*DIFABX*SX1*SY1*SZ0
     &                      - DIFABY*DIFABX*SX1*SZ1*SY0
     &                      + DIFABX*DIFABX*SY1*SZ1*SX0)*DP25
            SHLINT(INT,IZZ) = SHLINT(INT,IZZ)
     &                      - (D2*DIFABX*DIFABY*SX1*SY1*SZ0
     &                      - DIFABY*DIFABY*SX2*SY0*SZ0
     &                      - DIFABX*DIFABX*SY2*SX0*SZ0)*DP25
  600    CONTINUE
  300 CONTINUE
      RETURN
      END
C
C  /* Deck se1int */
      SUBROUTINE SE1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        H. Heiberg, -95
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
C
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
C
C$$$
         print *, 'shgtf',shgtf, 'int',int
         print *, 'odc(',LVALA,LVALB,0,0,0,1,')=',
     $        ODC(LVALA,LVALB,0,0,0,1)
         print *, 'odc(',MVALA,MVALB,0,0,0,2,')=',
     $        ODC(MVALA,MVALB,0,0,0,2)
         print *, 'odc(',NVALA,NVALB,0,0,0,3,')=',
     $        ODC(NVALA,NVALB,0,0,0,3)
         print *, 'odc(',LVALA,LVALB,0,0,1,1,')=',
     $        ODC(LVALA,LVALB,0,0,1,1)
         print *, 'odc(',MVALA,MVALB,0,0,1,2,')=',
     $        ODC(MVALA,MVALB,0,0,1,2)
         print *, 'odc(',NVALA,NVALB,0,0,1,3,')=',
     $        ODC(NVALA,NVALB,0,0,1,3)
C$$$
         DPLX = DX0*SY0*SZ0
         DPLY = SX0*DY0*SZ0
         DPLZ = SX0*SY0*DZ0
C$$$
         print *,'dplx=',dplx
         print *,'dply=',dply
         print *,'dplz=',dplz
C
         SHLINT(INT,1) = SHLINT(INT,1) - DPLX
         SHLINT(INT,2) = SHLINT(INT,2) - DPLY
         SHLINT(INT,3) = SHLINT(INT,3) - DPLZ
 100  CONTINUE
      RETURN
      END
C
C  /* Deck se1inb */
      SUBROUTINE SE1INB(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     $           SHLINT, EXPA, EXPB)
C
C        H. Heiberg, -95
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      PARAMETER (D4I = 0.25D0)
C
      EXPAIQ = D4I / (EXPA * EXPA)
      EXPBIQ = D4I / (EXPB * EXPB)
      DIFEXP = EXPAIQ - EXPBIQ
C
      INT = 0
C
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = DIFEXP*SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = DIFEXP*SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = DIFEXP*SHGTF*ODC(NVALA,NVALB,0,1,0,3)
C
         DPLX = DX0*SY0*SZ0
         DPLY = SX0*DY0*SZ0
         DPLZ = SX0*SY0*DZ0
C
         SHLINT(INT,1) = SHLINT(INT,1) - DPLX
         SHLINT(INT,2) = SHLINT(INT,2) - DPLY
         SHLINT(INT,3) = SHLINT(INT,3) - DPLZ
 100  CONTINUE
      RETURN
      END
C
C  /* Deck h1eint */
      SUBROUTINE H1EINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     $           SHLINT,AHGTF,NATOMC)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D4INV = 0.250D0, D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 0
C$$$
      print *,'shgtf=',shgtf
C$$$
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1  = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1  = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1  = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX2  = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DY2  = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DZ2  = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         DX21 = SHGTF*ODC(LVALA,LVALB,0,2,1,1)
         DY21 = SHGTF*ODC(MVALA,MVALB,0,2,1,2)
         DZ21 = SHGTF*ODC(NVALA,NVALB,0,2,1,3)
C
C       Kinetic energy contribution
C
         SHLINT(INT,1) = SHLINT(INT,1) - D2INV*( DX21*SY0*SZ0
     &                                         + SX1 *DY2*SZ0
     &                                         + SX1 *SY0*DZ2)
         SHLINT(INT,2) = SHLINT(INT,2) - D2INV*( DX2*SY1 *SZ0
     $                                         + SX0*DY21*SZ0
     $                                         + SX0*SY1 *DZ2)
         SHLINT(INT,3) = SHLINT(INT,3) - D2INV*( DX2*SY0*SZ1
     $                                         + SX0*DY2*SZ1
     $                                         + SX0*SY0*DZ21)
C
C       Nuclear attraction
C
c$$$         MAXT = LVALA + LVALB
c$$$         MAXU = MVALA + MVALB
c$$$         MAXV = NVALA + NVALB
c$$$         IADRAV = 1
c$$$         DO 200 IV = 0, MAXV + 1
c$$$            IADRAU = IADRAV
c$$$            EV = ODC(NVALA,NVALB,IV,0,0,3)
c$$$            FV = ODC(NVALA,NVALB,IV,0,1,3)
c$$$            IUMAX = MAXU + 1
c$$$            IF (IV .GT. MAXV) IUMAX = MAXU
c$$$            DO 300 IU = 0, IUMAX
c$$$               EU = ODC(MVALA,MVALB,IU,0,0,2)
c$$$               FU = ODC(MVALA,MVALB,IU,0,1,2)
c$$$               ITMAX = MAXT + 1
c$$$               IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) ITMAX = MAXT
c$$$               DO 400 IT = 0, ITMAX
c$$$                  ET = ODC(LVALA,LVALB,IT,0,0,1)
c$$$                  FT = ODC(LVALA,LVALB,IT,0,1,1)
c$$$                  FX = FT*EU*EV
c$$$                  FY = ET*FU*EV
c$$$                  FZ = ET*EU*FV
c$$$                  ATUV = -D2INV*AHGTF(IADRAU + IT)
c$$$                  SHLINT(INT,1)=SHLINT(INT,1)+ATUV*FX
c$$$                  SHLINT(INT,2)=SHLINT(INT,2)+ATUV*FY
c$$$                  SHLINT(INT,3)=SHLINT(INT,3)+ATUV*FZ
c$$$ 400           CONTINUE
c$$$               IADRAU = IADRAU + ISTEPU
c$$$ 300        CONTINUE
c$$$            IADRAV = IADRAV + ISTEPV
c$$$ 200     CONTINUE
 100  CONTINUE
      RETURN
      END
      SUBROUTINE RNSLOIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,DIFABX,DIFABY,DIFABZ)
C
C     Relativistic London orbital contribution to nuclear shielding integrals
C     ...the relativistic analogue to the NSLO nonrelativ.integrals..
C
C     T. Enevoldsen 10. Oct. 1997
C     MI febr.2006 - added multyplying by 1/2 to make integrals compatible with Dalton NSLO
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  FEE = FT*EU*EV
                  EFE = ET*FU*EV
                  EEF = ET*EU*FV
                  AH1T = (DIFABY*EEF-DIFABZ*EFE)
                  AH1U = (DIFABZ*FEE-DIFABX*EEF)
                  AH1V = (DIFABX*EFE-DIFABY*FEE)
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IXX = 9*(IATOM - 1) + 1
                     IXY = 9*(IATOM - 1) + 2
                     IXZ = 9*(IATOM - 1) + 3
                     IYX = 9*(IATOM - 1) + 4
                     IYY = 9*(IATOM - 1) + 5
                     IYZ = 9*(IATOM - 1) + 6
                     IZX = 9*(IATOM - 1) + 7
                     IZY = 9*(IATOM - 1) + 8
                     IZZ = 9*(IATOM - 1) + 9
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
C
                     SHLINT(INT,IXX) = SHLINT(INT,IXX) +
     &                    AH0T*AH1T*D2INV
                     SHLINT(INT,IYX) = SHLINT(INT,IYX) +
     &                    AH0T*AH1U*D2INV
                     SHLINT(INT,IZX) = SHLINT(INT,IZX) +
     &                    AH0T*AH1V*D2INV
                     SHLINT(INT,IXY) = SHLINT(INT,IXY) +
     &                    AH0U*AH1T*D2INV
                     SHLINT(INT,IYY) = SHLINT(INT,IYY) +
     &                    AH0U*AH1U*D2INV
                     SHLINT(INT,IZY) = SHLINT(INT,IZY) +
     &                    AH0U*AH1V*D2INV
                     SHLINT(INT,IXZ) = SHLINT(INT,IXZ) +
     &                    AH0V*AH1T*D2INV
                     SHLINT(INT,IYZ) = SHLINT(INT,IYZ) +
     &                    AH0V*AH1U*D2INV
                     SHLINT(INT,IZZ) = SHLINT(INT,IZZ) +
     &                    AH0V*AH1V*D2INV

                     IOFF = IOFF + NAHGTF
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END


      SUBROUTINE RM1H3(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C     Written by Thomas Enevoldsen 20. nov 1997
C     Integrals for Qmn<r(\alpha \cdot p)>
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,9)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - D2INV*
     &        (DIFABY*SZ1*SY0 - DIFABZ*SY1*SZ0)*DX0
         SHLINT(INT,2) = SHLINT(INT,2) - D2INV*
     &        (DIFABY*SZ1*DY0 - DIFABZ*DY1*SZ0)*SX0
         SHLINT(INT,3) = SHLINT(INT,3) - D2INV*
     &        (DIFABY*DZ1*SY0 - DIFABZ*SY1*DZ0)*SX0
         SHLINT(INT,4) = SHLINT(INT,4) - D2INV*
     &        (DIFABZ*DX1*SZ0 - DIFABX*SZ1*DX0)*SY0
         SHLINT(INT,5) = SHLINT(INT,5) - D2INV*
     &        (DIFABZ*SX1*SZ0 - DIFABX*SZ1*SX0)*DY0
         SHLINT(INT,6) = SHLINT(INT,6) - D2INV*
     &        (DIFABZ*SX1*DZ0 - DIFABX*DZ1*SX0)*SY0
         SHLINT(INT,7) = SHLINT(INT,7) - D2INV*
     &        (DIFABX*SY1*DX0 - DIFABY*DX1*SY0)*SZ0
         SHLINT(INT,8) = SHLINT(INT,8) - D2INV*
     &        (DIFABX*DY1*SX0 - DIFABY*SX1*DY0)*SZ0
         SHLINT(INT,9) = SHLINT(INT,9) - D2INV*
     &        (DIFABX*SY1*SX0 - DIFABY*SX1*SY0)*DZ0
C
 100  CONTINUE
      RETURN
      END
C  /* Deck RM1H2 */
      SUBROUTINE RM1H2(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC)
C
C     T. Enevoldsen
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D4INV = 0.250D0, D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
C
C     Nuclear attraction contribution
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV + 1
            IADRAU = IADRAV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            IUMAX = MAXU + 1
            IF (IV .GT. MAXV) IUMAX = MAXU
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               ITMAX = MAXT + 1
               IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) ITMAX = MAXT
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  FX = FT*EU*EV
                  FY = ET*FU*EV
                  FZ = ET*EU*FV
                  ATUV = -D2INV*AHGTF(IADRAU + IT)
                  SHLINT(INT,1)=SHLINT(INT,1)+ATUV*(DIFABY*FZ-DIFABZ*FY)
                  SHLINT(INT,2)=SHLINT(INT,2)+ATUV*(DIFABZ*FX-DIFABX*FZ)
                  SHLINT(INT,3)=SHLINT(INT,3)+ATUV*(DIFABX*FY-DIFABY*FX)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck pr1prv */
      SUBROUTINE PR1PRV(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  WORK,LWORK,LABINT,INTTYP,NOPTYP,
     &                  NBAST,IORDER,DOATOM,IPRINT,TOLS,TOLOG,NATOMC,
     &                  DISTAB,FACINT,CORCX,CORCY,CORCZ,DIFABX,DIFABY,
     &                  DIFABZ,DONUC1,DOMOM1,HUCFAC)
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "pi.h"
      PARAMETER (D1 = 1.00D00, D3 = 3.00D00, D3INV = D1/D3,
     &           D2INV = 0.50D0)
C
      PARAMETER (SMTHR = 1.00D-20)
      LOGICAL DOATOM(*), DIFODC, DONUC1, DOMOM1, ADDNAI, HUCFAC
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION WORK(LWORK), TEMP(3),
     &          FACINT(NATOMC), CORCX(NATOMC), CORCY(NATOMC),
     &          CORCZ(NATOMC), SHLINT(KCKTAB,*), AHGTF(*),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
C
#include "onecom.h"
#include "nuclei.h"
#include "orgcom.h"
#include "primit.h"
#include "symmet.h"
C
#include "ibtfun.h"
C
C Does only work for Tensor/Chiral effect!
C
      DIFODC = .FALSE.
C
C
C
C Exponent and coordinates of Spherical Nuclear distribution
      NUCPN1 = 0
      DO 50 IATOM = 1, NUCIND
      IF (DOATOM(IATOM)) THEN
      DO 60 IREP = 0, MAXREP
         IF (IBTAND(IREP,ISTBNU(IATOM)).NE.0) GOTO 60
      NUCPN1 = NUCPN1 + 1
      EXPN = GNUEXP(IATOM) 
      CORNX = CORCX(NUCPN1)
      CORNY = CORCY(NUCPN1)
      CORNZ = CORCZ(NUCPN1)
      CONTN = DSQRT((EXPN / PI)**3)
C
      DIFANX = CORAX - CORNX
      DIFANY = CORAY - CORNY
      DIFANZ = CORAZ - CORNZ
      DISTAN = DIFANX*DIFANX + DIFANY*DIFANY + DIFANZ*DIFANZ
      DIFBNX = CORBX - CORNX
      DIFBNY = CORBY - CORNY
      DIFBNZ = CORBZ - CORNZ
      DISTBN = DIFBNX*DIFBNX + DIFBNY*DIFBNY + DIFBNZ*DIFBNZ
C
      DO 100 IPRIMA = 1,NUCA
         JPRIMA = JSTA + IPRIMA
         CONTA = PRICCF(JPRIMA,NUMCFA)
         EXPA = PRIEXP(JPRIMA)
      DO 100 IPRIMB = 1,NUCB
         JPRIMB = JSTB + IPRIMB
         CONTB = PRICCF(JPRIMB,NUMCFB)
         EXPB = PRIEXP(JPRIMB)
         EXPP = EXPA + EXPB + EXPN
C
C     Tensor interaction
C
         EXPPI = D1/EXPP
         EXPABQ = EXPA*EXPB*DISTAB+EXPA*EXPN*DISTAN+EXPB*EXPN*DISTBN
         EXPABQ = EXPABQ*EXPPI
         IF (EXPABQ.GT.TOLOG) GO TO 200
         SAAB = CONTA*CONTB*CONTN*DEXP(-EXPABQ)
         ASAAB = DABS(SAAB)
         IF (ASAAB.LT.TOLS) GO TO 200
         SAAB13 = SIGN(ASAAB**D3INV,SAAB)
C
C        Calculate coordinates of product Gaussian P
C
         EXPAPI = EXPA*EXPPI
         EXPBPI = EXPB*EXPPI
         EXPNPI = EXPN*EXPPI
         CORPX  = EXPAPI*CORAX + EXPBPI*CORBX + EXPNPI*CORNX
         CORPY  = EXPAPI*CORAY + EXPBPI*CORBY + EXPNPI*CORNY
         CORPZ  = EXPAPI*CORAZ + EXPBPI*CORBZ + EXPNPI*CORNZ
C
C Check for same center (To avoid numerical noise, e.g. nonzero
C (p_x | \rho | s) on same center since CORPX - CORAX has noise
C of apprix 10^(-16)
C
         IF(DISTAB.LT.SMTHR) THEN
            IF(DISTAN.LT.SMTHR) THEN
               IF(DISTBN.LT.SMTHR) THEN
                  CORAX = CORNX
                  CORAY = CORNY
                  CORAZ = CORNZ
                  CORBX = CORNX
                  CORBY = CORNY
                  CORBZ = CORNZ
                  CORPX = CORNX
                  CORPY = CORNY
                  CORPZ = CORNZ
               ENDIF
            ENDIF
         ENDIF
C     
C        *********************************************
C        ***** Overlap Distribution Coefficients *****
C        *********************************************
C
         TEMP(1) = ORIGIN(1)
         TEMP(2) = ORIGIN(2)
         TEMP(3) = ORIGIN(3)
C
         CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
C        **********************************************
C        ***** Calculation of Hermitian integrals *****
C        **********************************************
C
C        Overlap integral
C
         SHGTF = DSQRT(PI*EXPPI)
C
C        **********************************************
C        ***** Calculation of Cartesian integrals *****
C        **********************************************
C
C        PV integrals
C        ------------
C
         CALL OVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &        D1, SHGTF,SHLINT(1,NUCPN1))
C
 200     CONTINUE
 100  CONTINUE

  60  CONTINUE
      ENDIF
  50  CONTINUE
C
C     Print
C
      IF (IPRINT .GT. 30) THEN
         DO 400 I = 1, NOPTYP
            CALL AROUND('SHLINT for '//LABINT(I)//' in PR1PRV')
            CALL OUTPUT(SHLINT(1,I),1,KCKTA,1,KCKTB,KCKTA,KCKTB,1,LUPRI)
  400    CONTINUE
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C    /* deck RM2H3  */
      SUBROUTINE RM2H3(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                 DIFABX,DIFABY,DIFABZ)
C*********************************************************************
C  Relativistic Magnetic 2nd derivative for Dirac Hamiltonian 
C             -- kinetic energy term
C
C  Integrals for 1/4 * <A_m | (Qmn r)(r Qmn)~ (c \alpha \cdot p) |B_n>
C
C  Written March 2003 by Miroslav Ilias and Hans Joergen Aa. Jensen.
C
C Note MI,Oct.2010: These integrals already contain derived factor of 1/4
C
C*********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D4INV = 0.25D0, D2 = 2.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,18)
#include "onecom.h"
#include "lmns.h"
C
      SHGTF3  = D4INV*SHGTF**3
      DIFABXX = DIFABX*DIFABX*SHGTF3
      DIFABXY = DIFABX*DIFABY*SHGTF3
      DIFABXZ = DIFABX*DIFABZ*SHGTF3
      DIFABYY = DIFABY*DIFABY*SHGTF3
      DIFABYZ = DIFABY*DIFABZ*SHGTF3
      DIFABZZ = DIFABZ*DIFABZ*SHGTF3
C
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = ODC(LVALA,LVALB,0,0,0,1)
         SY0 = ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = ODC(NVALA,NVALB,0,0,0,3)

         SX1 = ODC(LVALA,LVALB,0,0,1,1)
         SY1 = ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = ODC(NVALA,NVALB,0,0,1,3)

         SX2 = ODC(LVALA,LVALB,0,0,2,1)
         SY2 = ODC(MVALA,MVALB,0,0,2,2)
         SZ2 = ODC(NVALA,NVALB,0,0,2,3)

         DX0 = ODC(LVALA,LVALB,0,1,0,1)
         DY0 = ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = ODC(NVALA,NVALB,0,1,0,3)

         DX1 = ODC(LVALA,LVALB,0,1,1,1)
         DY1 = ODC(MVALA,MVALB,0,1,1,2)
         DZ1 = ODC(NVALA,NVALB,0,1,1,3)

         DX2 = ODC(LVALA,LVALB,0,1,2,1)
         DY2 = ODC(MVALA,MVALB,0,1,2,2)
         DZ2 = ODC(NVALA,NVALB,0,1,2,3)

C========================================================================
C
C       Qx = -Zmn.y + Ymn.z;  Qx~ = -Qx
C       Qy =  Zmn.x - Xmn.z;  Qy~ = -Qy
C       Qz = -Ymn.x + Xmn.y;  Qz~ = -Qz
C
C========================================================================
C    Qx*Qx~*grad(x,y,z)=-(-Zmn.y + Ymn.z)*(-Zmn.y + Ymn.z)*grad(x,y,z):
C
C           Qx.Qx~ = -(Zmn^2.y^2  - 2Ymn.Zmn.y.z + Ymn^2.z^2)
C========================================================================
         SHLINT(INT,1) = SHLINT(INT,1) + 
     &    (DIFABZZ*SY2*SZ0 - D2*DIFABYZ*SY1*SZ1 + DIFABYY*SZ2*SY0)*DX0
         SHLINT(INT,2) = SHLINT(INT,2) + 
     &    (DIFABZZ*DY2*SZ0 - D2*DIFABYZ*DY1*SZ1 + DIFABYY*SZ2*DY0)*SX0
         SHLINT(INT,3) = SHLINT(INT,3) + 
     &    (DIFABZZ*SY2*DZ0 - D2*DIFABYZ*SY1*DZ1 + DIFABYY*DZ2*SY0)*SX0
C========================================================================
C                      Qx*Qy*grad(x,y,z):
C
C     Qx.Qy~ = -(- Zmn^2.x.y  +  Xmn.Zmn.y.z + Ymn.Zmn.x.z - Xmn.Ymn.z^2)
C========================================================================
         SHLINT(INT,4) = SHLINT(INT,4)
     &         - DIFABZZ*DX1*SY1*SZ0 + DIFABXZ*DX0*SY1*SZ1
     &         + DIFABYZ*DX1*SY0*SZ1 - DIFABXY*DX0*SY0*SZ2
         SHLINT(INT,5) = SHLINT(INT,5)
     &         - DIFABZZ*SX1*DY1*SZ0 + DIFABXZ*SX0*DY1*SZ1
     &         + DIFABYZ*SX1*DY0*SZ1 - DIFABXY*SX0*DY0*SZ2
         SHLINT(INT,6) = SHLINT(INT,6)
     &         - DIFABZZ*SX1*SY1*DZ0 + DIFABXZ*SX0*SY1*DZ1
     &         + DIFABYZ*SX1*SY0*DZ1 - DIFABXY*SX0*SY0*DZ2
C========================================================================
C                     Qx*Qz*grad(x,y,z):
C
C     Qx.Qz~ = -(Ymn.Zmn.x.y - Xmn.Zmn.y^2 - Ymn^2.x.z + Xmn.Ymn.y.z)
C========================================================================
         SHLINT(INT,7) = SHLINT(INT,7)
     &         + DIFABYZ*DX1*SY1*SZ0 - DIFABXZ*DX0*SY2*SZ0
     &         - DIFABYY*DX1*SY0*SZ1 + DIFABXY*DX0*SY1*SZ1
         SHLINT(INT,8) = SHLINT(INT,8)
     &         + DIFABYZ*SX1*DY1*SZ0 - DIFABXZ*SX0*DY2*SZ0
     &         - DIFABYY*SX1*DY0*SZ1 + DIFABXY*SX0*DY1*SZ1
         SHLINT(INT,9) = SHLINT(INT,9)
     &         + DIFABYZ*SX1*SY1*DZ0 - DIFABXZ*SX0*SY2*DZ0
     &         - DIFABYY*SX1*SY0*DZ1 + DIFABXY*SX0*SY1*DZ1
C========================================================================
C                     Qy*Qy~*grad(x,y,z):
C
C     Qy  =  Zmn.x - Xmn.z
C
C     Qy.Qy~ = - ( Zmn^2.x^2  - 2.Xmn.Zmn.x.z + Xmn^2.z^2 )
C
C========================================================================
         SHLINT(INT,10) = SHLINT(INT,10)
     &         + DIFABZZ*DX2*SY0*SZ0 - D2*DIFABXZ*DX1*SY0*SZ1
     &         + DIFABXX*DX0*SY0*SZ2
         SHLINT(INT,11) = SHLINT(INT,11)
     &         + DIFABZZ*SX2*DY0*SZ0 - D2*DIFABXZ*SX1*DY0*SZ1
     &         + DIFABXX*SX0*DY0*SZ2
         SHLINT(INT,12) = SHLINT(INT,12)
     &         + DIFABZZ*SX2*SY0*DZ0 - D2*DIFABXZ*SX1*SY0*DZ1
     &         + DIFABXX*SX0*SY0*DZ2
C========================================================================
C                  Qy*Qz~*grad(x,y,z):
C
C   Qy.Qz~ = -(-Ymn.Zmn.x^2 + Zmn.Xmn.x.y + Xmn.Ymn.x.z - Xmn^2.y.z)
C========================================================================
         SHLINT(INT,13) = SHLINT(INT,13)
     &       - DIFABYZ*DX2*SY0*SZ0 + DIFABXZ*DX1*SY1*SZ0
     &       + DIFABXY*DX1*SY0*SZ1 - DIFABXX*DX0*SY1*SZ1
         SHLINT(INT,14) = SHLINT(INT,14)
     &       - DIFABYZ*SX2*DY0*SZ0 + DIFABXZ*SX1*DY1*SZ0
     &       + DIFABXY*SX1*DY0*SZ1 - DIFABXX*SX0*DY1*SZ1
         SHLINT(INT,15) = SHLINT(INT,15)
     &       - DIFABYZ*SX2*SY0*DZ0 + DIFABXZ*SX1*SY1*DZ0
     &       + DIFABXY*SX1*SY0*DZ1 - DIFABXX*SX0*SY1*DZ1
C========================================================================
C                 Qz*Qz~*grad(x,y,z):
C
C     Qz = -Ymn.x + Xmn.y
C
C    Qz.Qz~ =  -( Ymn^2.x^2 - 2.Xmn.Ymn.x.y + Xmn^2.y^2 )
C========================================================================
         SHLINT(INT,16) = SHLINT(INT,16) + 
     &     (DIFABYY*DX2*SY0 - D2*DIFABXY*DX1*SY1
     &          +  DIFABXX*DX0*SY2) * SZ0
         SHLINT(INT,17) = SHLINT(INT,17) +
     &    (DIFABYY*SX2*DY0 - D2*DIFABXY*SX1*DY1
     &          + DIFABXX*SX0*DY2) * SZ0
         SHLINT(INT,18) = SHLINT(INT,18) +
     &    (DIFABYY*SX2*SY0 - D2*DIFABXY*SX1*SY1
     &          + DIFABXX*SX0*SY2) * DZ0

 100  CONTINUE
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rdsusan */
      SUBROUTINE RDSUSAN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                SHLINT,CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ)
C*********************************************************************
C
C  RDSUSLL:
C  Integrals for i c Qmn <mi_M| r (r_N x alpha) |ni_N>
C
C  Written March 2003 by Miroslav Ilias and Hans Joergen Aa. Jensen,
C  after routine DSUSAN by Kenneth Ruud.
C
C MI,note,Oct.2012: All terms are multiplied with factor 1/2
C
C*********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D2INV = 0.50D0, D4INV = 0.25D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,18)
#include "onecom.h"
#include "lmns.h"

C     write(lupri,'(A,5I5)') 'RDSUSAN JMAXA/B/T/D/M',
C    *  JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
C     write(lupri,'(A,3F20.10)') 'RDSUSAN DIFAB',DIFABX,DIFABY,DIFABZ
C     write(lupri,'(A,5I5)') 'RDSUSAN KCKTA/B',KCKTA,KCKTB

      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX00 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)   ! 1
         SY00 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)   ! 1
         SZ00 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)   ! 1

         SX10 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)   ! x
         SY10 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)   ! y
         SZ10 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)   ! z

         SX01 = SHGTF*ODC(LVALA,LVALB+1,0,0,0,1) ! xn
         SY01 = SHGTF*ODC(MVALA,MVALB+1,0,0,0,2) ! yn
         SZ01 = SHGTF*ODC(NVALA,NVALB+1,0,0,0,3) ! zn

         SX11 = SHGTF*ODC(LVALA,LVALB+1,0,0,1,1) ! x.xn
         SY11 = SHGTF*ODC(MVALA,MVALB+1,0,0,1,2) ! y.yn
         SZ11 = SHGTF*ODC(NVALA,NVALB+1,0,0,1,3) ! z.zn
          
      
C     write(lupri,'(A,3I5)') 'RDSUSAN LMN A',LVALA,MVALA,NVALA
C     write(lupri,'(A,3I5)') 'RDSUSAN LMN B',LVALB,MVALB,NVALB
C     write(lupri,'(A,1P,3D15.5)') 'RDSUSAN S?00',SX00,SY00,SZ00
C     write(lupri,'(A,1P,3D15.5)') 'RDSUSAN S?10',SX10,SY10,SZ10
C     write(lupri,'(A,1P,3D15.5)') 'RDSUSAN S?01',SX01,SY01,SZ01
C     write(lupri,'(A,1P,3D15.5)') 'RDSUSAN S?11',SX11,SY11,SZ11
C=================================================================
C      1. Qx.yn.   alpha_z           XXRDSULY
C = (-Zmn.y.yn + Ymn.z.yn)
C=================================================================
         SHLINT(INT,1) = SHLINT(INT,1) - 
     & (-DIFABZ*SX00*SY11*SZ00+DIFABY*SX00*SY01*SZ10)*D2INV
C=================================================================
C      2. Qx.zn.   alpha_y           XXRDSULZ
C =   (-Zmn.y.zn + Ymn.z.zn)
C=================================================================
         SHLINT(INT,2) = SHLINT(INT,2) - 
     & (-DIFABZ*SX00*SY10*SZ01+DIFABY*SX00*SY00*SZ11)*D2INV
C=================================================================
C      3. Qx.zn. alpha_x             XYRDSULZ
C =     (-Zmn.y.zn + Ymn.z.zn)
C=================================================================
         SHLINT(INT,3) = SHLINT(INT,3) - 
     & (-DIFABZ*SX00*SY10*SZ01+DIFABY*SX00*SY00*SZ11)*D2INV
C=================================================================
C      4. Qx.xn. alpha_z             XYRDSULX
C =     (-Zmn.y.xn + Ymn.z.xn)
C=================================================================
         SHLINT(INT,4) = SHLINT(INT,4) - 
     & (-DIFABZ*SX01*SY10*SZ00+DIFABY*SX01*SY00*SZ10)*D2INV
C=================================================================
C      5. Qx.xn.  alpha_y            XZRDSULX
C =    (-Zmn.y.xn + Ymn.z.xn)
C=================================================================
         SHLINT(INT,5) = SHLINT(INT,5) - 
     & (-DIFABZ*SX01*SY10*SZ00+DIFABY*SX01*SY00*SZ10)*D2INV
C=================================================================
C      6. Qx.yn. alpha_x             XZRDSULY
C =     (-Zmn.y.yn + Ymn.z.yn)
C=================================================================
         SHLINT(INT,6) = SHLINT(INT,6) - 
     & (-DIFABZ*SX00*SY11*SZ00+DIFABY*SX00*SY01*SZ10)*D2INV
C=================================================================
C      7. Qy.yn.   alpha_z           YXRDSULY
C =   (Zmn.x.yn - Xmn.z.yn)
C=================================================================
         SHLINT(INT,7) = SHLINT(INT,7) - 
     & (DIFABZ*SX10*SY01*SZ00-DIFABX*SX00*SY01*SZ10)*D2INV
C=================================================================
C      8. Qy.zn.   alpha_y             YXRDSULZ
C =    (Zmn.x.zn - Xmn.z.zn)
C=================================================================
         SHLINT(INT,8) = SHLINT(INT,8) - 
     & (DIFABZ*SX10*SY00*SZ01-DIFABX*SX00*SY00*SZ11)*D2INV
C=================================================================
C      9. Qy.zn.   alpha_x              YYRDSULZ
C =      (Zmn.x.zn - Xmn.z.zn)
C=================================================================
         SHLINT(INT,9) = SHLINT(INT,9) - 
     & (DIFABZ*SX10*SY00*SZ01-DIFABX*SX00*SY00*SZ11)*D2INV
C=================================================================
C      10.  Qy.xn.  alpha_z                 YYRDSULX 
C =        (Zmn.x.xn - Xmn.z.xn)
C=================================================================
         SHLINT(INT,10) = SHLINT(INT,10) - 
     & (DIFABZ*SX11*SY00*SZ00-DIFABX*SX01*SY00*SZ10)*D2INV
C=================================================================
C      11. Qy.xn.  alpha_y                  YZRDSULX
C =       (Zmn.x.xn - Xmn.z.xn)
C=================================================================
         SHLINT(INT,11) = SHLINT(INT,11) - 
     & (DIFABZ*SX11*SY00*SZ00-DIFABX*SX01*SY00*SZ10)*D2INV
C=================================================================
C      12. Qy.yn.   alpha_x                 YZRDSULY
C =        (Zmn.x.yn - Xmn.z.yn)
C=================================================================
         SHLINT(INT,12) = SHLINT(INT,12) - 
     & (DIFABZ*SX10*SY01*SZ00-DIFABX*SX00*SY01*SZ10)*D2INV
C=================================================================
C      13.  Qz.yn.  alpha_z                 ZXRDSULY
C  =       (-Ymn.x.yn + Xmn.y.yn)
C=================================================================
         SHLINT(INT,13) = SHLINT(INT,13) - 
     & (-DIFABY*SX10*SY01*SZ00+DIFABX*SX00*SY11*SZ00)*D2INV
C=================================================================
C      14.  Qz.zn.   alpha_y                 ZXRDSULZ
C  =       (-Ymn.x.zn + Xmn.y.zn)
C=================================================================
         SHLINT(INT,14) = SHLINT(INT,14) - 
     & (-DIFABY*SX10*SY00*SZ01+DIFABX*SX00*SY10*SZ01)*D2INV
C=================================================================
C      15. Qz.zn.   alpha_x      ZYRDSULZ
C  =       (-Ymn.x.zn + Xmn.y.zn)
C=================================================================
         SHLINT(INT,15) = SHLINT(INT,15) - 
     & (-DIFABY*SX10*SY00*SZ01+DIFABX*SX00*SY10*SZ01)*D2INV
C=================================================================
C      16. Qz.xn.   alpha_z       ZYRDSULX
C  =       (-Ymn.x.xn + Xmn.y.xn)
C=================================================================
         SHLINT(INT,16) = SHLINT(INT,16) - 
     & (-DIFABY*SX11*SY00*SZ00+DIFABX*SX01*SY10*SZ00)*D2INV
C=================================================================
C      17. Qz.xn.  alpha_y        ZZRDSULX
C =        (-Ymn.x.xn + Xmn.y.xn)
C=================================================================
         SHLINT(INT,17) = SHLINT(INT,17) - 
     & (-DIFABY*SX11*SY00*SZ00+DIFABX*SX01*SY10*SZ00)*D2INV
C=================================================================
C      18.  Qz.yn.   alpha_x      ZZRDSULY
C =         (-Ymn.x.yn +  Xmn.y.yn)
C=================================================================
         SHLINT(INT,18) = SHLINT(INT,18) - 
     & (-DIFABY*SX10*SY01*SZ00 + DIFABX*SX00*SY11*SZ00)*D2INV

 100  CONTINUE
      RETURN
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dsuslh */
      SUBROUTINE RM2H2(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC)
C***************************************************************************************
C     K.Ruud, Feb. 1992
C     H.J.Aa.Jensen (&MI) Mar 2003: new flag DOEKIN, if false omit kinetic energy terms
C     (used by Dirac for DSUSLH because Vnuc is same in 4-comp, but kin.en. is not)
C
C Note,MI,Oct.2012:  integrals -Z_A/r_iA already contain factor of +1/4, resulting in -1/4
C
C***************************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D4INV = 0.250D0, D8INV = 0.1250D0, D2 = 2.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
C     Nuclear attraction contribution
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV + 2
            IADRAU = IADRAV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            GV = ODC(NVALA,NVALB,IV,0,2,3)
            IUMAX = MAXU + 2
            IF (IV .GT. MAXV + 1) THEN
               IUMAX = MAXU
            ELSE IF (IV .GT. MAXV) THEN
               IUMAX = MAXU + 1
            END IF
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               GU = ODC(MVALA,MVALB,IU,0,2,2)
               ITMAX = MAXT + 2
               IF (((IU .GT. MAXU + 1) .OR. (IV .GT. MAXV + 1))
     &               .OR. ((IU .EQ. MAXU + 1) .AND. (IV .EQ. MAXV + 1)))
     &               THEN
                  ITMAX = MAXT
               ELSE IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) THEN
                  ITMAX = ITMAX + 1
               END IF
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  GT = ODC(LVALA,LVALB,IT,0,2,1)
                  FXX = GT*EU*EV
                  FXY = FT*FU*EV
                  FXZ = FT*EU*FV
                  FYY = ET*GU*EV
                  FYZ = ET*FU*FV
                  FZZ = ET*EU*GV
                  ATUV = -D4INV*AHGTF(IADRAU + IT)
                  SHLINT(INT,1)=SHLINT(INT,1)+ATUV*(DIFABY*DIFABZ*FYZ*D2
     &                                            - DIFABY*DIFABY*FZZ
     &                                            - DIFABZ*DIFABZ*FYY)
                  SHLINT(INT,2)=SHLINT(INT,2)+ATUV*(DIFABZ*DIFABZ*FXY
     &                                            + DIFABY*DIFABX*FZZ
     &                                            - DIFABZ*DIFABX*FYZ
     &                                            - DIFABY*DIFABZ*FXZ)
                  SHLINT(INT,3)=SHLINT(INT,3)+ATUV*(DIFABX*DIFABZ*FYY
     &                                            + DIFABY*DIFABY*FXZ
     &                                            - DIFABY*DIFABZ*FXY
     &                                            - DIFABX*DIFABY*FYZ)
                  SHLINT(INT,4)=SHLINT(INT,4)+ATUV*(DIFABX*DIFABZ*FXZ*D2
     &                                            - DIFABX*DIFABX*FZZ
     &                                            - DIFABZ*DIFABZ*FXX)
                  SHLINT(INT,5)=SHLINT(INT,5)+ATUV*(DIFABX*DIFABX*FYZ
     &                                            + DIFABY*DIFABZ*FXX
     &                                            - DIFABX*DIFABY*FXZ
     &                                            - DIFABX*DIFABZ*FXY)
                  SHLINT(INT,6)=SHLINT(INT,6)+ATUV*(DIFABX*DIFABY*FXY*D2
     &                                            - DIFABY*DIFABY*FXX
     &                                            - DIFABX*DIFABX*FYY)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
 100  CONTINUE
      RETURN
      END


      SUBROUTINE G1HTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NCENTC,NATOMC,
     &                  ISYMOP,NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
C
C     tuh
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP,3), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          NCENTC(NATOMC),INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
C
      FULMAT = .TRUE.
C
C     Loops over ireps and Cartesian coordinats
C
      DO 100 IREP = 0, MAXREP
      DO 100 ICOOR = 1, 3
C
C        Center A
C        ========
C 
         ICOORA = 3*(ICENTA - 1) + ICOOR
         ISCORA = IPTCNT(3*(NCENTA - 1) + ICOOR,IREP,1)
         IF (ISCORA .GT. 0) THEN
            FC = HKAB
            IF (IREP .EQ. 0) THEN
               CALL SYM1S(SHLINT(1,ICOORA,1),SOINT(1,ISCORA),
     &                    ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,DUM,IDUM,IPRINT)
            ELSE
               CALL SYM1N(SHLINT(1,ICOORA,1),SOINT(1,ISCORA),
     &                    IREP,ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
            END IF
         END IF
C
C        Center B
C        ========
C
         ICOORB = 3*(ICENTB - 1) + ICOOR
         ISCORB = IPTCNT(3*(NCENTB - 1) + ICOOR,IREP,1)
         IF (ISCORB .GT. 0) THEN
            FC = HKAB*PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
     &               *PT(IBTAND(IREP,ISYMOP))
            IF (IREP .EQ. 0) THEN
               CALL SYM1S(SHLINT(1,ICOORB,2),SOINT(1,ISCORB),
     &                    ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,DUM,IDUM,IPRINT)
            ELSE
               CALL SYM1N(SHLINT(1,ICOORB,2),SOINT(1,ISCORB),
     &                    IREP,ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
            END IF
         END IF
C
C        Center C
C        ========
C
         DO 200 IATOM = 1, NATOMC
            ICOORC = 3*(NCENTC(IATOM) - 1) + ICOOR
            ISCORC = IPTCNT(3*(JCENTC(IATOM) - 1) + ICOOR,IREP,1)
            IF (ISCORC .GT. 0) THEN
               FC = HKAB*PT(IBTAND(ISYMAX(ICOOR,1),JSYMC(IATOM)))
     &                  *PT(IBTAND(IREP,JSYMC(IATOM)))
               IF (IREP .EQ. 0) THEN
                  CALL SYM1S(SHLINT(1,ICOORC,3),SOINT(1,ISCORC),
     &                       ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                       FC,LDIAG,FULMAT,DUM,IDUM,IPRINT)
               ELSE
                  CALL SYM1N(SHLINT(1,ICOORC,3),SOINT(1,ISCORC),
     &                       IREP,ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                       KHKTB,FC,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
               END IF
            END IF
  200    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END

      SUBROUTINE S1HTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,IPRINT)
C
C     tuh
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,3), SOINT(NELMNT,NOPTYP)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
C
      FULMAT = .TRUE.
C
C     Loops over ireps and Cartesian coordinats
C
      DO 100 IREP = 0, MAXREP
      DO 100 ICOOR = 1, 3
C
C        Center A
C        ========
C 
         ISCORA = IPTCNT(3*(NCENTA - 1) + ICOOR,IREP,1)
         IF (ISCORA .GT. 0) THEN
            FC = HKAB
            IF (IREP .EQ. 0) THEN
               CALL SYM1S(SHLINT(1,ICOOR),SOINT(1,ISCORA),
     &                    ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,DUM,IDUM,IPRINT)
            ELSE
               CALL SYM1N(SHLINT(1,ICOOR),SOINT(1,ISCORA),
     &                    IREP,ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
            END IF
         END IF
C
C        Center B
C        ========
C
         ISCORB = IPTCNT(3*(NCENTB - 1) + ICOOR,IREP,1)
         IF (ISCORB .GT. 0) THEN
            FC = -HKAB*PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
     &               *PT(IBTAND(IREP,ISYMOP))
            IF (IREP .EQ. 0) THEN
               CALL SYM1S(SHLINT(1,ICOOR),SOINT(1,ISCORB),
     &                    ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,DUM,IDUM,IPRINT)
            ELSE
               CALL SYM1N(SHLINT(1,ICOOR),SOINT(1,ISCORB),
     &                    IREP,ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
            END IF
         END IF
C
  100 CONTINUE
C
      RETURN
      END

C  /* Deck dpgint */
      SUBROUTINE DPGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, July-97
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,9)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         INT = INT + 1
C
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
C     
         SHLINT(INT,1) = SHLINT(INT,1) - DX1*SY0*SZ0
         SHLINT(INT,2) = SHLINT(INT,2) - DX0*SY1*SZ0
         SHLINT(INT,3) = SHLINT(INT,3) - DX0*SY0*SZ1
         SHLINT(INT,4) = SHLINT(INT,4) - SX1*DY0*SZ0
         SHLINT(INT,5) = SHLINT(INT,5) - SX0*DY1*SZ0
         SHLINT(INT,6) = SHLINT(INT,6) - SX0*DY0*SZ1
         SHLINT(INT,7) = SHLINT(INT,7) - SX1*SY0*DZ0
         SHLINT(INT,8) = SHLINT(INT,8) - SX0*SY1*DZ0
         SHLINT(INT,9) = SHLINT(INT,9) - SX0*SY0*DZ1
 100  CONTINUE
      RETURN
      END
C  /* Deck qugint */
      SUBROUTINE QUGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, October/November-97
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,18)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,1,2,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,1,2,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,1,2,3)
         INT = INT + 1
C     
         SHLINT(INT,1) = SHLINT(INT,1) - DX2*SY0*SZ0
         SHLINT(INT,2) = SHLINT(INT,2) - DX1*SY1*SZ0
         SHLINT(INT,3) = SHLINT(INT,3) - DX1*SY0*SZ1
         SHLINT(INT,4) = SHLINT(INT,4) - DX0*SY2*SZ0
         SHLINT(INT,5) = SHLINT(INT,5) - DX0*SY1*SZ1
         SHLINT(INT,6) = SHLINT(INT,6) - DX0*SY0*SZ2
         SHLINT(INT,7) = SHLINT(INT,7) - SX2*DY0*SZ0
         SHLINT(INT,8) = SHLINT(INT,8) - SX1*DY1*SZ0
         SHLINT(INT,9) = SHLINT(INT,9) - SX1*DY0*SZ1
         SHLINT(INT,10) = SHLINT(INT,10) - SX0*DY2*SZ0
         SHLINT(INT,11) = SHLINT(INT,11) - SX0*DY1*SZ1
         SHLINT(INT,12) = SHLINT(INT,12) - SX0*DY0*SZ2
         SHLINT(INT,13) = SHLINT(INT,13) - SX2*SY0*DZ0
         SHLINT(INT,14) = SHLINT(INT,14) - SX1*SY1*DZ0
         SHLINT(INT,15) = SHLINT(INT,15) - SX1*SY0*DZ1
         SHLINT(INT,16) = SHLINT(INT,16) - SX0*SY2*DZ0
         SHLINT(INT,17) = SHLINT(INT,17) - SX0*SY1*DZ1
         SHLINT(INT,18) = SHLINT(INT,18) - SX0*SY0*DZ2
 100  CONTINUE
      RETURN
      END

C  /* Deck dpgtra */
      SUBROUTINE DPGTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER B, X
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,18), SOINT(NELMNT,NOPTYP),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ICOMP = 0

      ! loop over cartesian displacements of the atoms
      DO 100 X = 1, 3
         JCOMP = 0
 
         ! loop over components of the dipole operator
         DO 200 B = 1, 3
            ICOMP = ICOMP + 1
            JCOMP = JCOMP + 1
            IREPC  = ISYMAX(B,1) ! irrep of dipole component

            ! loop over possible irreps of symmetry coordinates
            DO IREP = 0, MAXREP
               ISCOOR = IPTCNT(3*(NCENTA - 1) + X,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(3*(ISCOOR - 1) + JCOMP)
                  IREPO = IBTXOR(IREP,IREPC) ! irrep of operator
ckr                  FACSYM = HKAB*PT(IBTAND(IREPO,ISYMOP))
C                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
C     &                         *PT(IBTAND(IREP,ISYMOP))
ckr                  FACSYM = HKAB
                  FACSYM = HKAB
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP),SOINT(1,ITYP),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP),SOINT(1,ITYP),IREPO,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
               ISCOOR = IPTCNT(3*(NCENTB - 1) + X,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(3*(ISCOOR - 1) + JCOMP)
                  IREPO = IBTXOR(IREP,IREPC) ! irrep of operator
ckr                  FACSYM = HKAB*PT(IBTAND(IREPO,ISYMOP))
C                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
C    &                          *PT(IBTAND(IREP,ISYMOP))
                   FACSYM = HKAB*PT(IBTAND(ISYMAX(X,1),ISYMOP))
     &                          *PT(IBTAND(IREP,ISYMOP))
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP + 9),SOINT(1,ITYP),
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,DUM,IDUM,
     &                          IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP + 9),SOINT(1,ITYP),
     &                          IREPO,ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &                          KHKTA,KHKTB,FACSYM,LDIAG,FULMAT,ANTI,
     &                          DUM,IDUM,IPRINT)
                  END IF
               END IF
            END DO
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck qugtra */
      SUBROUTINE QUGTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER B, X, C
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,36), SOINT(NELMNT,NOPTYP),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ICOMP = 0
      DO 100 X = 1, 3
         JCOMP = 0
         DO 200 B = 1, 3
         DO 200 C = B, 3
            ICOMP = ICOMP + 1
            JCOMP = JCOMP + 1
            IREPC  = IBTXOR(ISYMAX(B,1),ISYMAX(C,1))
            DO IREPO = 0, MAXREP
               ISCOOR = IPTCNT(3*(NCENTA - 1) + X,IREPO,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(6*(ISCOOR - 1) + JCOMP)
                  IREPCL = IBTXOR(IREPO,IREPC)
ckr                  FACSYM = HKAB*PT(IBTAND(IREPCL,ISYMOP))
                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
     &                         *PT(IBTAND(IREPO,ISYMOP))
ckr                  FACSYM = HKAB
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP),SOINT(1,ITYP),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP),SOINT(1,ITYP),IREPO,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
               ISCOOR = IPTCNT(3*(NCENTB - 1) + X,IREPO,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(6*(ISCOOR - 1) + JCOMP)
                  IREPCL = IBTXOR(IREPO,IREPC)
ckr                  FACSYM = HKAB*PT(IBTAND(IREPCL,ISYMOP))
                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
     &                         *PT(IBTAND(IREPO,ISYMOP))
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP + 18),SOINT(1,ITYP),
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,DUM,IDUM,
     &                          IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP + 18),SOINT(1,ITYP),
     &                          IREPO,ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &                          KHKTA,KHKTB,FACSYM,LDIAG,FULMAT,ANTI,
     &                          DUM,IDUM,IPRINT)
                  END IF
               END IF
            END DO
 200     CONTINUE
 100  CONTINUE
      RETURN
      END


