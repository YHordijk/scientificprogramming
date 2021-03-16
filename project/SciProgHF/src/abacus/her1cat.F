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
C  /* Deck catdrv */
      SUBROUTINE CATDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,
     &                  NODV)
C*****************************************************************************
C
C   CATDRV: Coulomb attraction integrals
C           (based on ONEDRV)
C
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "iratdef.h"
C
      LOGICAL DIFINT, NODC, NODV, PROPTY, DOINT(2,2), TOFILE
      CHARACTER*4 OMITVNUC
      DIMENSION WORK(LWORK)
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "dorps.h"
C
#include "ibtfun.h"
C
      DOINT(1,1) = .TRUE.
      DOINT(2,1) = .TRUE.
      DOINT(1,2) = .TRUE.
      DOINT(2,2) = .TRUE.
      OMITVNUC   = 'FFFF'
      TOFILE     = .TRUE.
C
C     ***** Number of basis functions *****
C
      NBAST  = 0
      NNBAST = 0
      DO 100 KB = 0, MAXREP
         NBASI = 0
         DO 200 ISHELL = 1, KMAX
            IF (IBTAND(KB,ISTBAO(ISHELL)).EQ.0) THEN
               NBASI = NBASI + KHKT(ISHELL)
            END IF
  200    CONTINUE
         NBAST  = NBAST  + NBASI
         NNBAST = NNBAST + NBASI*(NBASI + 1)/2
  100 CONTINUE
      NNBASX = NBAST*(NBAST + 1)/2
C
      MXCENTC = NUCDEP
C
      KVCMAT = 1
      KDENMA = KVCMAT +   NNBASX
      KFACCBA= KDENMA +   NNBASX
      KCOORC = KFACCBA+   NUCDEP*NUCDEP*MXCENTC
C     allocate for nuclear charge AND an electronic core charge
C     for modified nuclear attraction (e.g. for the small charge
C     in relativistic calculations). /Mar2001 hjaaj
      KSIGNC = KCOORC + 3*NUCDEP
      KJSYMC = KSIGNC + 3*NUCDEP
      KJCENT = KJSYMC +  (NUCDEP + 1)/IRAT
      KGEXP  = KJCENT +  (NUCDEP + 1)/IRAT
      KLAST  = KGEXP  +  NUCDEP
      LWRK   = LWORK  - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('CATDRV',' ',KLAST,LWORK)
#if !defined(PRG_DIRAC)
C
C     **************************************************************
C     ***** Set up total density and Fock matrices in AO basis *****
C     **************************************************************
C
      IF (PROPTY) THEN
         KDSO  = KLAST
         KFSO  = KDSO  + NNBAST
         KFOCMA= KFSO  + NNBAST
         KLAST2= KFOCMA+ NNBASX
         IF(KLAST2.GT.LWORK) CALL STOPIT('CATDR1','DSOFSO',KLAST2,LWORK)
         CALL DSOFSO(WORK(KDSO),WORK(KFSO),WORK(KLAST),LWRK,IPRINT,NODC,
     &               NODV)
         CALL DSYM1(WORK(KDENMA),WORK(KFOCMA),WORK(KDSO),WORK(KFSO),
     &              NBAST,IPRINT)
      END IF
#endif
      CALL CATNUC(MXCENTC,NCENTC,CHARGE,GNUEXP, WORK(KFACCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),OMITVNUC,IPRINT)
      CALL DZERO(WORK(KVCMAT),NNBASX)
      CALL CATDR1(WORK(KVCMAT),WORK(KDENMA),NCENTC,
     &            WORK(KFACCBA),WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KLAST),
     &            LWRK,IPRINT,PROPTY,MAXDIF,
     &            NNBASX,DOINT)
C
C       Write undifferentiated one-electron integrals
C
      IF (PROPTY) THEN
         IF (MAXDIF .GE. 2) THEN
           IF (IPRINT .GE. 5) CALL SHDPRI(WORK,LWORK)
           IF (DOREPS(0)) CALL SHDCHK(WORK,LWORK,NODC,NODV,IPRINT)
         END IF
      ELSE
         IF (TOFILE) THEN
            CALL WRTUND(WORK(KVCMAT),NBAST,NNBASX,IPRINT)
         ENDIF
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck catnuc */
      SUBROUTINE CATNUC(MXCENTC,NCENTC,XCHARGE,FACCBA,GEXP,XGEXP,
     &                  COORC,SIGNC,JSYMC,JCENTC,OMITVNUC,IPRINT)
C*****************************************************************************
C
C     CATNUC: Make list of atomic centers for calculation
C             of nuclear attraction integrals with CATDR1.
C
C     Oct 2001 HJAaJ
C
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0)
C
      DIMENSION XCHARGE(NUCIND), COORC(3,MXCENTC), SIGNC(3,MXCENTC),
     &          JSYMC(MXCENTC), JCENTC(MXCENTC),GEXP(MXCENTC),
     &          XGEXP(NUCIND),FACCBA(MXCENTC,NUCDEP,NUCDEP)
      CHARACTER*4 OMITVNUC
#include "nuclei.h"
#include "symmet.h"
C
#include "ibtfun.h"
      IF (IPRINT .GE. 5) CALL TITLER('Output from CATNUC','*',103)
C
C
C     ************************************************************
C     ***** Set up list of centers for attraction integrals   ****
C     ************************************************************
C
         NCENTC = 0
         DO 120 IATOMC = 1,NUCIND
            MULC   = ISTBNU(IATOMC)
            CORCX0 = CORD(1,IATOMC)
            CORCY0 = CORD(2,IATOMC)
            CORCZ0 = CORD(3,IATOMC)
            CHARG1 = XCHARGE(IATOMC)
            DO 130 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,MULC) .EQ. 0) THEN
                  NCENTC = NCENTC + 1
                  IF (NCENTC.GT.MXCENTC) THEN
                     CALL QUIT('CATNUC: NCENTC.gt.MXCENTC')
                  END IF
                  JSYMC(NCENTC)   = ISYMOP
                  JCENTC(NCENTC)  = IATOMC
                  SIGNC(1,NCENTC) = PT(IBTAND(ISYMAX(1,1),ISYMOP))
                  SIGNC(2,NCENTC) = PT(IBTAND(ISYMAX(2,1),ISYMOP))
                  SIGNC(3,NCENTC) = PT(IBTAND(ISYMAX(3,1),ISYMOP))
                  COORC(1,NCENTC) = SIGNC(1,NCENTC)*CORCX0
                  COORC(2,NCENTC) = SIGNC(2,NCENTC)*CORCY0
                  COORC(3,NCENTC) = SIGNC(3,NCENTC)*CORCZ0
                  GEXP(NCENTC)    = XGEXP(IATOMC)
                  CALL DCOPY(NUCDEP*NUCDEP,-CHARG1,0,
     &                       FACCBA(NCENTC,1,1),MXCENTC)
               END IF
  130       CONTINUE
  120    CONTINUE
C
C *****
C   For example, to omit < C | V_C | C > 
C    FACCBA(i1,i,i) = D0 when i1 corresponds to center i
C
C   OMITVNUC(1:1) = 'T' : omit < C | V_C | C > 
C   OMITVNUC(2:2) = 'T' : omit < A | V_C | A > 
C   OMITVNUC(3:3) = 'T' : omit < A | V_C | C > and  < C | V_C | B > 
C   OMITVNUC(4:4) = 'T' : omit < A | V_C | B >
C *****
C
C     1) intialize to include everything
C
      IF (NCENTC.NE.MXCENTC) CALL QUIT('CATNUC: NCENTC.ne.MXCENTC')
C
C     2) set factor to 0 for terms to be omitted
C
      IF (OMITVNUC(4:4) .EQ. 'T') THEN
C        OMITVNUC(4:4) = 'T' : omit < A | V_C | B >
         DO K = 1, NUCDEP
            DO J = 1, NUCDEP
            IF (J.NE.K) THEN
               TMP1 = FACCBA(J,J,K)
               TMP2 = FACCBA(K,J,K)
               DO I = 1,NCENTC
                  FACCBA(I,J,K) = D0
               END DO
               FACCBA(J,J,K) = TMP1
               FACCBA(K,J,K) = TMP2
            END IF
            END DO
         END DO
      END IF
      IF (OMITVNUC(3:3) .EQ. 'T') THEN
C        OMITVNUC(3:3) = 'T' : omit < A | V_C | C > and  < C | V_C | B > 
         DO J = 1, NUCDEP
            TMP1 = FACCBA(J,J,J)
            DO I = 1,NCENTC
               FACCBA(I,I,J) = D0
               FACCBA(I,J,I) = D0
            END DO
            FACCBA(J,J,J) = TMP1
         END DO
      END IF
      IF (OMITVNUC(2:2) .EQ. 'T') THEN
C        OMITVNUC(2:2) = 'T' : omit < A | V_C | A > 
         DO J = 1, NUCDEP
            TMP1 = FACCBA(J,J,J)
            DO I = 1,NCENTC
               FACCBA(I,J,J) = D0
            END DO
            FACCBA(J,J,J) = TMP1
         END DO
      END IF
      IF (OMITVNUC(1:1) .EQ. 'T') THEN
C        OMITVNUC(1:1) = 'T' : omit < C | V_C | C > 
         DO I = 1,NCENTC
            FACCBA(I,I,I) = D0
         END DO
      END IF
C      IF(IPRINT .GE. 5) THEN
c        DO I = 1,NCENTC
chj sep 02: commented because PRQMAT does not print
c           the desired matrix !!!
c           WRITE(LUPRI,'(/A)') 
c    &      '*****************************************'
c           WRITE(LUPRI,'(A,I3,A)') 
c    &      '*********** NCENTC = ',I,' ****************'
c           WRITE(LUPRI,'(A/)') 
c    &      '*****************************************'
c           CALL PRQMAT(FACCBA(I,1,1),NUCDEP,NUCDEP,
c    &                  NUCDEP,NUCDEP,1,1,LUPRI) 
c        ENDDO
C      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck catdr1 */
      SUBROUTINE CATDR1(VCMAT,DENMAT,NCENTC,FACCBA,COORC,SIGNC,
     &                  JSYMC,JCENTC,GEXP,WORK,LWORK,IPRINT,
     &                  PROPTY,MAXDIF,NNBASX,DOINT)
C*****************************************************************************
C
C     CATDR1: Driver for Coulomb attraction integrals
C             (based on ONEDR1 by TUH)
C
C     < A | V_cat | B > = Sum(C = 1, NCENTC)) < A | V_C | B >
C
C     where V_C = Q_C / r_C (or spherical Gaussian charge distribution)
C
C     The integrals are ADDED to VCMAT.
C
C     Sep 2001 JKP
C
C     This program calculates Coulomb attraction
C     matrix elements and the first and second derivatives of these
C     elements using the McMurchie/Davidson scheme.  See L. E. McMurchie
C     & E. R. Davidson, J. Comp. Phys. 26 (1978) 218, and also V. R.
C     Saunders in "Methods in Computational Molecular Physics", Reidel
C     1983.
C
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
      PARAMETER (LUTEMP = 48)
      PARAMETER (D0 = 0.00D00)
C
      LOGICAL FRSDER, SECDER, PROPTY, DOINT(2,2)
      DIMENSION DENMAT(*), VCMAT(NNBASX), FACCBA(*),
     &          WORK(LWORK), COORC(3,NCENTC),
     &          SIGNC(3,NCENTC), JSYMC(NCENTC),
     &          JCENTC(NCENTC), GEXP(NCENTC)
#include "onecom.h"
#include "lmns.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
#include "dipole.h"
#include "ccom.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "symind.h"
#include "csym1.h"
C
#include "ibtfun.h"
      IF (IPRINT .GE. 5) CALL TITLER('Output from CATDR1','*',103)
C
      TIMHER = D0
      TIMINT = D0
      FRSDER = MAXDIF .GE. 1
      SECDER = MAXDIF .EQ. 2
      TOLS   = THRS*THRS
      TOLOG  = - LOG(TOLS)
      IF (PROPTY .AND. SECDER) THEN
         REWIND LUTEMP
         INDMAX = 0
         LENGTH = 0
      END IF
C
C
C     Gradient and Hessian elements equal to zero
C
      IF (PROPTY) THEN
         ENERNA = D0
#ifdef PRG_DIRAC
C
         CALL DZERO(GRADNU(1,0),3*NUCDEP)
         CALL DZERO(GRADNU(1,1),3*NUCDEP)
         CALL DZERO(GRADNU(1,2),3*NUCDEP)
C
#else 
         CALL DZERO(GRADNA,3*NUCDEP)
         IF (SECDER) HESSNA(:,:) = 0.0D0
#endif
      END IF
C
      KLAST = 1
      LWRK  = LWORK - KLAST + 1
C
C     ************************************************************
C     ***** Triangular loop over symmetry independent shells *****
C     ************************************************************
C
      CALL IZERO(ISOFRA, 8)
      IDENA = 0
      DO 100 ISHELA = 1,KMAX
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
         IDENB0 = 0
C
C        Compute symmetry integral pointers for contributions
C        from this block.  Note that at present this assumes all
C        components from a shell are included.
C
         DO 600 I = 1, 8
            ISOFRB(I) = 0
            DO 610 J = 1, MXAQN
              INDFA(I,J) = -10 000 000
610         CONTINUE
600      CONTINUE
         DO 620 NA = 1, KHKTA
            DO 630 IREP = 0, MAXREP
            IF (IBTAND(MULA,IBTXOR(IREP,ISYMAO(NHKTA,NA))).EQ.0) THEN
               ISOFRA(IREP+1)    = ISOFRA(IREP+1) + 1
               INDFA (IREP+1,NA) = ISOFRA(IREP+1)
            END IF
630         CONTINUE
620      CONTINUE
         IF (IPRINT .GT. 20) THEN
            WRITE(LUPRI,'(A,I4)')' IA address offsets for shell ',ISHELA
            DO 640 NA = 1,KHKTA
               WRITE(LUPRI,'(8(1X,I5))') (INDFA(I,NA), I = 1,MAXREP+1)
640         CONTINUE
         END IF
      DO 110 ISHELB = 1,ISHELA
         LDIAG = ISHELA .EQ. ISHELB
         NHKTB = NHKT(ISHELB)
         KHKTB = KHKT(ISHELB)
         KCKTB = KCKT(ISHELB)
         ICB   = LCLASS(ISHELB)
         SPHRB = SPHR(ISHELB)
         CALL LMNVAL(NHKTB,KCKTB,LVALUB,MVALUB,NVALUB)
         NCENTB = NCENT(ISHELB)
         NHKTAB = NHKTA + NHKTB
         MULB   = ISTBAO(ISHELB)
         MULTB  = MULT(MULB)
         NUCB   = NUCO(ISHELB)
         NUMCFB = NUMCF(ISHELB)
         JSTB   = JSTRT(ISHELB)
         CORBX0 = CENT(ISHELB,1,1)
         CORBY0 = CENT(ISHELB,2,1)
         CORBZ0 = CENT(ISHELB,3,1)
         KHKTAB = KHKTA*KHKTB
         KCKTAB = KCKTA*KCKTB
         MAB    = IBTOR(MULA,MULB)
         KAB    = IBTAND(MULA,MULB)
         HKAB   = FMULT(KAB)
C
         SPHRAB = SPHRA .OR. SPHRB
C
C        Compute symmetry integral pointers for contributions
C        from this block.  Note that at present this assumes all
C        components from a shell are included
C
         DO 700 I = 1, 8
            DO 710 J = 1, MXAQN
              INDFB(I,J) = -10 000 000
710         CONTINUE
700      CONTINUE
         DO 720 NB = 1, KHKTB
            DO 730 IREP = 0, MAXREP
            IF (IBTAND(MULB,IBTXOR(IREP,ISYMAO(NHKTB,NB))).EQ.0) THEN
               ISOFRB(IREP+1)    = ISOFRB(IREP+1) + 1
               INDFB (IREP+1,NB) = ISOFRB(IREP+1)
            END IF
730         CONTINUE
720      CONTINUE
         IF (IPRINT .GT. 20) THEN
            WRITE(LUPRI,'(A,I4)')' IB address offsets for shell ',ISHELB
            DO 740 NB = 1, KHKTB
               WRITE(LUPRI,'(8(1X,I5))') (INDFB(I,NB), I = 1,MAXREP+1)
740         CONTINUE
         ENDIF
         IF(.NOT.DOINT(ICA,ICB)) GOTO 110
#ifdef PRG_DIRAC
C        jth - CLASS (which integral class: LL,LS,SL,SS)
C                                           0  1  1  2
         IF (ICA .EQ. 1) THEN
            IF (ICB .EQ. 1) THEN
               INTCLASS = 0
            ELSE
               INTCLASS = 1
            END IF
         ELSE 
            IF (ICB .EQ. 1) THEN
               INTCLASS = 1
            ELSE
               INTCLASS = 2
            END IF
         END IF           
         IF (ICA.NE.ICB) THEN
            CALL QUIT('CATDR1 called with LS or SL shells')
         ENDIF
#endif
         IF (IPRINT .GE. 05) WRITE (LUPRI, 1000) ISHELA, ISHELB
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
C        Initialization for nuclear attraction integrals
C
         JMAX = NHKTAB - 2
         IF (PROPTY) JMAX = JMAX + MAXDIF
         ISTEPU = JMAX + 1
         ISTEPV = ISTEPU*ISTEPU
         NAHGTF = ISTEPU*ISTEPV
C
         CALL CATSOP(VCMAT,DENMAT,FACCBA,COORC,WORK(KLAST),LWRK,
     &               IPRINT,PROPTY,MAXDIF,IDENB0,CORBX0,CORBY0,CORBZ0,
     &               SECDER,NCENTC,TOLOG,TOLS,JSYMC,JCENTC,
     &               SIGNC,GEXP)
  110    IDENB0 = IDENB0 + KHKTB*MULTB
         IDENA = IDENA + KHKTA*MULTA
  100 CONTINUE
C
C     ***** End loop over symmetry independent orbitals *****
C
C     ****************************************************
C     ***** Write final buffers and dipole integrals *****
C     ***** ELSE write undifferentiated integrals    *****
C     ****************************************************
C
      IF (PROPTY .AND. SECDER) THEN
C
C        Write final buffer of first derivative integrals
C
         IF (LENGTH .GT. 0) WRITE (LUTEMP) BUF, IBUF, LENGTH
         WRITE (LUTEMP) BUF, IBUF, -LENGTH
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I4,A)')
     *       ' Last buffer of length',LENGTH,
     *       ' has been written on LUTEMP in CATDR1.'
C
C     *******************************
C     ***** Symmetrize Hessians *****
C     *******************************
C
#ifndef PRG_DIRAC
         DO 500 I = 2,3*NUCDEP
            DO 510 J = 1,I-1
               HESSNA(I,J) = HESSNA(I,J) + HESSNA(J,I)
               HESSNA(J,I) = HESSNA(I,J)
  510       CONTINUE
  500    CONTINUE
#endif
      END IF
C
C     *************************
C     ***** Print Section *****
C     *************************
C
#ifdef PRG_DIRAC
      IF (PROPTY .AND. (IPRINT .GE. 5)) THEN
         KCSTRA = 1
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
         CALL HEADER('Nuclear attraction integral gradient',-1)
         CALL PRIGRD(GRADNU(1,1),WORK(KCSTRA),WORK(KSCTRA))
         CALL PRIGRD(GRADNU(1,2),WORK(KCSTRA),WORK(KSCTRA))
#else 
      IF (PROPTY .AND. (IPRINT .GT. 0)) THEN
         CALL HEADER('Nuclear attraction integral gradient',-1)
         CALL PRIGRD(GRADNA)
#endif
         IF (SECDER) THEN
            CALL HEADER('Nuclear attraction integral Hessian',-1)
            CALL PRIHES(HESSNA,'CENTERS')
         END IF
      END IF
C
C     ********************************************************
C     ***** Sort gradient elements on direct access unit *****
C     ********************************************************
C
      IF (PROPTY .AND. SECDER) THEN
         NMATS = 3*NUCDEP*(MAXREP+1)
         CALL SORONE(WORK,LWORK,NMATS,INDMAX,IPRINT)
         CLOSE(LUTEMP, STATUS='DELETE')
      END IF
      RETURN
 1000 FORMAT (//,2X,'***************************************',
     *         /,2X,'********** ISHELA/B =',I3,',',I3,' **********',
     *         /,2X,'***************************************',/)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck catsop */
      SUBROUTINE CATSOP(VCMAT,DENMAT,FACCBA,COORC,WORK,LWORK,
     &                  IPRINT,PROPTY,MAXDIF,IDENB0,CORBX0,CORBY0,
     &                  CORBZ0,SECDER,NCENTC,TOLOG,TOLS,JSYMC,
     &                  JCENTC,SIGNC,GEXP)
C*****************************************************************************
C
C CATSOP: Coulomb attraction integrals, symmetry adaption
C         (based on ONESOP/ONESO1)
C
C output: VCMAT with Coulomb attraction integrals
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL SECDER, PROPTY
      DIMENSION WORK(LWORK)
      DIMENSION VCMAT(*), DENMAT(*), FACCBA(*),
     &          COORC(3,*), JSYMC(*), JCENTC(*), SIGNC(3,*),GEXP(*)
#include "onecom.h"
C
      IF (PROPTY) THEN
         NADER = 28*NCENTC*KCKTAB
      ELSE
         NADER = KCKTAB
      END IF
      KADER  = 1
      KDSHEL = KADER  +   NADER
      KLAST  = KDSHEL +   KHKTAB
      LWRK   = LWORK  - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('CATSOP',' ',KLAST,LWORK)
C
      CALL CATSO1(VCMAT, WORK(KADER),FACCBA,
     &            DENMAT,COORC,WORK(KLAST),LWRK,IPRINT,
     &            PROPTY,MAXDIF,IDENB0,CORBX0,CORBY0,CORBZ0,
     &            SECDER,NCENTC,TOLOG,TOLS,JSYMC,JCENTC,SIGNC,
     &            GEXP,WORK(KDSHEL))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck catso1 */
      SUBROUTINE CATSO1(VCMAT,ADER,FACCBA, DENMAT,COORC,WORK,LWORK,
     &                  IPRINT,PROPTY,MAXDIF,IDENB0,CORBX0,CORBY0,
     &                  CORBZ0,SECDER,NCENTC,TOLOG,TOLS,JSYMC,
     &                  JCENTC,SIGNC,GEXP,DSHELL)
C*****************************************************************************
C
C CATSOP: Coulomb attraction integrals, symmetry adaption
C         (based on ONESOP/ONESO1)
C
C output: VCMAT with Coulomb attraction integrals
C
C*****************************************************************************

      use quaternion_algebra
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
C
      LOGICAL SECDER, PROPTY, FULMAT
      DIMENSION WORK(LWORK)
      DIMENSION VCMAT(*), ADER(*), FACCBA(NCENTC,NUCDEP,NUCDEP),
     &          DSHELL(KHKTAB), DENMAT(*),
     &          COORC(3,*), JSYMC(*), JCENTC(*), SIGNC(3,*),GEXP(*)
#include "ccom.h"
#include "onecom.h"
#include "nuclei.h"
#include "symmet.h"
C
#ifdef PRG_DIRAC
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgrd.h"
#else
#include "energy.h"
#endif 
#include "ibtfun.h"
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
#ifdef PRG_DIRAC
C     Element (I,J,IZ) in full quaternion matrix
c     IQFULL(I,J,IZ) = MAX(I,J) + MIN(I-1,J-1) * NTBAS(0) +
c    &     N2BBASX * (IZ-1) 
      IOLDQFULL(I,J,IZ) = I + (J-1) * NTBAS(0) + N2BBASX * (IZ-1) 
#endif 
      IF (IPRINT .GT. 4) CALL TITLER('Output from CATSO1','*',103)
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
            WRITE (LUPRI,'(A,1I8)') ' INTCLASS', INTCLASS
            WRITE (LUPRI,'(A,1F12.6)') ' CORBX    ', CORBX
            WRITE (LUPRI,'(A,1F12.6)') ' CORBY    ', CORBY
            WRITE (LUPRI,'(A,1F12.6)') ' CORBZ    ', CORBZ
         END IF
C
C        **********************************************
C        ***** Calculation of Cartesian integrals *****
C        **********************************************
C
         CALL CATPRM(ADER, WORK,LWORK,IPRINT,PROPTY,MAXDIF,
     &               NCENTC,TOLOG,TOLS,SECDER,FACCBA(1,ICENTB,ICENTA),
     &               COORC,GEXP)
C
C        *************************************************
C        ***** Transform to spherical harmonic basis *****
C        *************************************************
C
         IF (SPHRAB) THEN
            IF (PROPTY) THEN
               IF (ONECEN) THEN
                  CALL SPHRM1(ADER,ADER,10*NCENTC,WORK,LWORK,PROPTY,
     &                        IPRINT)
               ELSE
                  CALL SPHRM1(ADER,ADER,28*NCENTC,WORK,LWORK,PROPTY,
     &                        IPRINT)
               END IF
            ELSE
               CALL SPHRM1(ADER,ADER,1,WORK,LWORK,PROPTY,IPRINT)
            END IF
         END IF
C
C        ******************************
C        ***** Expectation values *****
C        ******************************
C
         IF (PROPTY) THEN
C
C           Coulomb attraction integrals
C           ===========================================================
C
C           Collect density and Fock elements
C
            MAXCMP = 0
            DO 200 IORBA = IDENA + 1, IDENA + KHKTA
               DO 250 IORBB = IDENB + 1, IDENB + KHKTB
                  MAXCMP = MAXCMP + 1
                  FAC = D1
#ifdef PRG_DIRAC
Cjth              Due to the way DENMAT are constructed in Dirac don't 
C                 do this...
Cjth - sib        IF (ONECEN .AND. LDIAG .AND. IORBB.NE.IORBA) FAC = DP5
c                 IORBAB = IQFULL(IORBA,IORBB,1)
                  IORBAB = IOLDQFULL(IORBA,IORBB,1)
C                 ... only real part (IZ=1) needed for Coulomb attraction
#else 
                  IF (ONECEN .AND. LDIAG .AND. IORBB.NE.IORBA) FAC = DP5
                  IORBAB = ITRI(IORBA,IORBB)
#endif 
                  DSHELL(MAXCMP) = FAC*DENMAT(IORBAB)
  250          CONTINUE
  200       CONTINUE
C
            CALL QUIT('CAT1SOP: NCENTC not defined for AVENA1/2')
C           ... NCENTC renamed here to KCENTC(NUCDEP), KCENTC(I) = I
#ifdef PRG_DIRAC
C           SL,LS : nothing
C           LL,SS : Coulomb attraction
C 
            IF ( INTCLASS .EQ. 1 ) THEN
C
C                 SL, LS: nothing here in CATSO1
C
            ELSE
C
C              GRD_NOSMLV: No small-small nuclear attraction 
C                          (eg. when using Levy-Leblond Hamiltonian)
C
               IF ( ONECEN ) THEN
C
C                 Derivatives of the nuclear attraction, one center case
C
                  IF (.NOT. (GRD_NOSMLV .AND. (INTCLASS .EQ. 2))) THEN
                     CALL AVENA1(ADER,NCENTC,SECDER,KCENTC,JCENTC,
     &                           MAXCMP,JSYMC,SIGNC,DSHELL)
                  END IF
C
               ELSE
C
C                 Derivatives of the nuclear attraction, two center case
C
                  IF (.NOT. (GRD_NOSMLV .AND. (INTCLASS .EQ. 2))) THEN
                     CALL AVENA2(ADER,NCENTC,ISYMOP,SECDER,KCENTC,
     &                           MAXCMP,JCENTC,JSYMC,SIGNC,DSHELL)
                  END IF
               END IF
            END IF
            IF (IPRINT .GT. 20) THEN
               KCSTRA = 1
               KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
               KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
               CALL HEADER('Coulomb attraction integral gradient',-1)
               CALL PRIGRD(GRADNU(1,1),WORK(KCSTRA),WORK(KSCTRA))
               CALL PRIGRD(GRADNU(1,2),WORK(KCSTRA),WORK(KSCTRA))
            END IF
#else 
            IF (ONECEN) THEN
C
C              Nuclear attraction
C
               CALL AVENA1(ADER,NCENTC,SECDER,KCENTC,JCENTC,MAXCMP,
     &                     JSYMC,SIGNC,DSHELL)
            ELSE
C
C              Nuclear attraction
C
               CALL AVENA2(ADER,NCENTC,ISYMOP,SECDER,KCENTC,MAXCMP,
     &                     JCENTC,JSYMC,SIGNC,DSHELL)
            END IF
            IF (IPRINT .GT. 50) THEN
               CALL HEADER('Coulomb attraction integral gradient',-1)
               CALL PRIGRD(GRADNA)
            END IF
#endif
         END IF
C
C        *******************************************
C        ***** Transform integrals to SO basis *****
C        *******************************************
C
         FULMAT = .TRUE.
C
C        Coulomb attraction integrals
C        ============================
C
         CALL SYM1S(ADER,VCMAT,ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &              KHKTA,KHKTB,HKAB,LDIAG,FULMAT,DUMMY,IDUMMY,
     &              IPRINT)
C
C        **************************************************
C        ***** Write differentiated integrals on file *****
C        **************************************************
C
         IF (PROPTY .AND. SECDER) THEN
C            CALL WD1SYM(STDER1,ADER,RLMINT,FCM,WORK,LWORK,JSYMC,JCENTC,
C     &                  ISYMOP,NCENTC,IPRINT)
             CALL QUIT('CATSO1 : CALL WD1SYM NOT IMPLEMENTED ')
C            jkp sep 2001: need a flag in WD1SYM to only
C            calculate nuclear attraction type contributions.
C
         END IF
      END IF
  100 CONTINUE
      RETURN
 1010 FORMAT (//,2X,'***************************************',
     *         /,2X,'******** Symmetry operation ',I2,' ********',
     *         /,2X,'***************************************',/)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck catprm */
      SUBROUTINE CATPRM(ADER,WORK,LWORK,IPRINT,PROPTY,MAXDIF,NCENTC,
     &                  TOLOG,TOLS,SECDER,FACINT,COORC,GEXP)
C*****************************************************************************
C
C CATPRM: Coulomb attraction cartesian integrals,
C         based on ONEPRM/ONEPR1 in her1car.F
C
C Sep 2001 - hjaaj+jkp
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL PROPTY, SECDER
      DIMENSION WORK(LWORK)
      DIMENSION ADER(*), FACINT(*), COORC(3,*), GEXP(*)
#include "onecom.h"
C
      JMAXD = 2
      IF (PROPTY .AND. .NOT.ONECEN) JMAXD = 4
C
      JMAXA = NHKTA - 1
      JMAXB = NHKTB - 1
      JMAXT = JMAXA + JMAXB + JMAXD
      JMAXM = 0
C
      KAHGTF = 1
      KODC   = KAHGTF + (NCENTC + 1)*NAHGTF
      KLAST  = KODC   + 3*(JMAXA+1)*(JMAXB+1)*(JMAXT+1)*(JMAXD+1)
      LWRK   = LWORK  - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('CATPRM',' ',KLAST,LWORK)
      CALL CATPR1(ADER,
     &            WORK(KAHGTF),WORK(KODC),JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &            IPRINT,PROPTY,MAXDIF,NCENTC,TOLOG,TOLS,SECDER,
     &            FACINT,COORC,GEXP,NHKTAP,
     &            WORK(KLAST),LWRK)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck catpr1 */
      SUBROUTINE CATPR1(ADER,AHGTF,ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  IPRINT,PROPTY,MAXDIF,NCENTC,TOLOG,TOLS,SECDER,
     &                  FACINT,COORC,GEXP,NHKTAP,WORK,LWORK)
C*****************************************************************************
C
C CATPRM: Coulomb attraction cartesian integrals,
C         based on ONEPRM/ONEPR1 in her1car.F
C
C Sep 2001 - hjaaj+jkp
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "aovec.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.00 D00, D3 = 3.00 D00, 
     &           THIRD = D1/D3, GNUTHR = 1.D-4,D1P5 = 1.5D0)
C
      LOGICAL PROPTY, KINODC, SECDER, DOST, DONA
      PARAMETER (DOST = .FALSE., DONA = .TRUE.)
C     ... do not Sder and Tder, do "nucatt' der
      DIMENSION WORK(LWORK), ADER(*), AHGTF(*), 
     &          FACINT(*), COORC(3,*), GEXP(*)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
C
#include "onecom.h"
#include "ader.h"
#include "primit.h"
C
      IF (IPRINT .GE. 5) CALL TITLER('Output from CATPR1','*',103)
C
      DIFABX = CORAX - CORBX
      DIFABY = CORAY - CORBY
      DIFABZ = CORAZ - CORBZ
      DISTAB = DIFABX*DIFABX + DIFABY*DIFABY + DIFABZ*DIFABZ
C
      NATOT  = 0
      IF (PROPTY) THEN
         IF (ONECEN) THEN
            IA0000 =  1
            IA000X =  2
            IA000Y =  3
            IA000Z =  4
            IA00XX =  5
            IA00XY =  6
            IA00XZ =  7
            IA00YY =  8
            IA00YZ =  9
            IA00ZZ = 10
            NATOT  = 10*NCENTC*KCKTAB
            CALL DZERO(ADER,NATOT)
         ELSE
C
            IA0000 =  1
            IA0X00 =  2
            IA0Y00 =  3
            IA0Z00 =  4
            IAXX00 =  5
            IAXY00 =  6
            IAXZ00 =  7
            IAYY00 =  8
            IAYZ00 =  9
            IAZZ00 = 10
            IA000X = 11
            IA000Y = 12
            IA000Z = 13
            IA00XX = 14
            IA00XY = 15
            IA00XZ = 16
            IA00YY = 17
            IA00YZ = 18
            IA00ZZ = 19
            IA0X0X = 20
            IA0X0Y = 21
            IA0X0Z = 22
            IA0Y0X = 23
            IA0Y0Y = 24
            IA0Y0Z = 25
            IA0Z0X = 26
            IA0Z0Y = 27
            IA0Z0Z = 28
            NATOT  = 28*NCENTC*KCKTAB
            CALL DZERO(ADER,NATOT)
         END IF
      ELSE
         CALL DZERO(ADER,KCKTAB)
      END IF
C
C     ****************************************
C     ***** Loop over primitive orbitals *****
C     ****************************************
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
C
C           Calculate and test square root factor
C
            EXPABQ = EXPA*EXPB*DISTAB*EXPPI
         IF (EXPABQ.GT.TOLOG) GO TO 200
            SAAB = CONTA*CONTB*EXP(-EXPABQ)
            ASAAB = ABS(SAAB)
         IF (ASAAB.LT.TOLS) GO TO 200
            SAAB13 = SIGN(ASAAB**THIRD,SAAB)
C
C           Calculate coordinates of product Gaussian P
C
            EXPAPI = EXPA*EXPPI
            EXPBPI = EXPB*EXPPI
            CORPX  = EXPAPI*CORAX + EXPBPI*CORBX
            CORPY  = EXPAPI*CORAY + EXPBPI*CORBY
            CORPZ  = EXPAPI*CORAZ + EXPBPI*CORBZ
C
C           *********************************************
C           ***** Overlap Distribution Coefficients *****
C           *********************************************
C
C           Expansion coefficients for undifferentiated orbitals
C
            KINODC = .FALSE.
            IDUMMY = 0
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,PROPTY,KINODC,
     &                  ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,WORK,LWORK,
     &                  CORPX,CORPY,CORPZ,.TRUE.,.FALSE.,ORIGIN,IDUMMY)
C
C           **********************************************
C           ***** Calculation of Hermitian integrals *****
C           **********************************************
C
C           Overlap integral
C
            SHGTF = SQRT(PI*EXPPI)
C
C           Nuclear attraction integrals
C
               IADR = 1
               DO 300 ICENTC = 1,NCENTC
                  FACTOR = FACINT(ICENTC)
                  IF ( FACTOR .EQ. D0 ) THEN
                     CALL DZERO(AHGTF(IADR),NAHGTF)
                     IADR = IADR + NAHGTF
                     GO TO 300
                  END IF
                     DIFCPX = COORC(1,ICENTC) - CORPX
                     DIFCPY = COORC(2,ICENTC) - CORPY
                     DIFCPZ = COORC(3,ICENTC) - CORPZ
                     IF(ABS(GEXP(ICENTC)).GT.GNUTHR) THEN
                       EXPFAC = GEXP(ICENTC)
                       EXPFAC = EXPFAC/(EXPFAC + EXPP)
                       EXPPGN = EXPP*EXPFAC
                       FACTOR = FACTOR*(EXPFAC**D1P5)
                       CALL HERNAI(AHGTF,JMAX,EXPPGN,DIFCPX,DIFCPY,
     *                   DIFCPZ,FACTOR,IADR,ISTEPU,ISTEPV,NAHGTF,
     *                   IPRINT)
                     ELSE
                       CALL HERNAI(AHGTF,JMAX,EXPP,DIFCPX,DIFCPY,
     *                   DIFCPZ,FACTOR,IADR,ISTEPU,ISTEPV,NAHGTF,
     *                   IPRINT)
                     ENDIF
                     IADR = IADR + NAHGTF
C                  ENDIF
  300          CONTINUE
C
C           **********************************************
C           ***** Calculation of Cartesian integrals *****
C           **********************************************
C
            IF (PROPTY) THEN
               IF (ONECEN) THEN
                CALL QUIT('CAT1PRM: NCENTC not defined for CINT1')
C           ... NCENTC renamed here to KCENTC(NUCDEP), KCENTC(I) = I
                  CALL CINT1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,ADER,
     &                       AHGTF,KCENTC,NCENTC,SECDER)
               ELSE
                  CALL CINT2(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,STDER0,
     &                       STDER1,STDER2,ADER,SHGTF,AHGTF,NCENTC,
     &                       DOST,DONA,SECDER)
               END IF
            ELSE
               CALL CINT0(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,STDER0,
     &                    SHGTF,ADER,AHGTF,NCENTC,IPRINT,DOST,DONA)
            END IF
  200    CONTINUE
  100 CONTINUE
C
C     Print
C
      IF (IPRINT .GE. 5) THEN
chj sep 01: maybe print ADER here ?
c        IF (.NOT. ONECEN) THEN
c           WRITE (LUPRI,'(/A/)') ' Undifferentiated overlap integrals '
c           WRITE (LUPRI,'(1P,6D12.4)') (STDER0(I,1),I=1,KCKTAB)
c        END IF
      END IF
      RETURN
      END
