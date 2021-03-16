#define GASCIP_DEBUG -1
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

! FILE    : krmcgascip.F
! main authors: H.J.Aa. Jensen & J. Thyssen

! define task symbols for CALL DIRAC_PARCTL( task )
#include "dirac_partask.h"
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GASCIP_SETUP(UCIBOS,ICIREP,IPRINT)
C***********************************************************************
C
C     Setup and initialize variables in dcbgascip.h
C
C     Input : UCIBOS, ICIREP, IPRINT
C
C     Output: in common blocks
C
C
C     Written by J. Thyssen - Jan 2 2001
C     Last revision : hjaaj Aug 2014
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"

      LOGICAL UCIBOS
      INTEGER ICIREP, IPRINT
C
C Used from include files:
C  dcbopt.h : IOPT_MINMK2, ...
C  dcborb.h : NAELEC, NASHT, NASH(), IASH()
C
#include "dcbidx.h"
#include "dcbbos.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbgascip.h"
#include "pgroup.h"
C
C
      CALL QENTER('GASCIP_SETUP')
C
C     Setup symmetry mult. table
C     --------------------------
C
C     GMULTSF/GMULTA returns:
C       NREP     : #(boson irreps) = #(fermion irreps);
C       REPNA    : strings identifying the irreps
C       MULTB    : mult. table
C
      IF (UCIBOS) THEN
C        Use boson symmetry (for spin free and Levy-Leblond)
C        or use approximate boson symmetry
         CALL GMULTSF(NREP,REPNA,MULTB(1,1))
C        hjaaj aug 2003: NOTE that this info is NOT used in this
C        version of GASCIP !!!!!!!!!!!!!
C        (MULTB is not used as it includes spin-symmetry -- see GMULTSF --
C        and we treat the spin Ms separately as in non-rel codes)
      ELSE
         CALL GMULTA(NZ,NFSYM,NREP,REPNA,MULTB(1,1))
      ENDIF
C     We extend MULTB(1:2*NREP,1:2*NREP) with MULTB(1:2*NREP,0) for use in
C     GASCIP_SYMOK when we add the first spinor to a new
C     string (i.e. symmetry zero in second index means vacuum)
      DO I = 1, 2 * NREP
         MULTB(I,0) = I
      END DO
C
C     Setup symmetry
C     --------------
C
C     IRP  : ICIREP is user input for CI symmetry
C     IRRP : symmetry of wave function according to:
C            1 .. NREP           : fermion irreps
C            NREP + 1 .. 2 * NREP: boson irreps
C
C             fermion               boson
C     C1   :  A                     a
C     Ci   :  AG, AU                ag, au
C
C     C2   :  1E, 2E                a, b
C     Cs   :  1E, 2E                a, b
C     C2h  :  1Eg, 2Eg, 1Eu, 2Eu    ag, bg, au, bu
C
C     C2v  :  1E, 2E                a, b
C     D2h  :  1Eg, 2Eg, 1Eu, 2Eu    ag, bg, au, bu
C
C
      IRP = ICIREP
      IF (IRP.LE.0) THEN
C        ... code for include all symmetries
         IRRP = IRP
      ELSE IF ( MOD ( NAELEC, 2 ) .EQ. 0 ) THEN
C        ...boson symmetry
         IRRP = NREP + IRP
      ELSE
C        ...fermion symmetry
         IRRP = IRP
      END IF
C
C     Setup orbital symmetries
C     ------------------------
C
C     IRPAMO is the symmetry of each spinor.
C     IREPA counts fermion irrep of double group (see GMULTA)
C     (use boson irreps without spin for SPINFR .OR. LEVYLE,
C      which implies UCIBOS true, or use approx. boson symmetry
C      if UCIBOS true for relativistic calculation)
C
      II    = 0
      IREPA = 0
      DO I = 1, NFSYM
         IREPA = IREPA + 1
         DO ISTR = 1, NASH(I)
            II = II + 1
            IF (UCIBOS) THEN
               IG = IDXU2G(IASH(I)+ISTR)
               IRPAMO(II) = IBOSYM(IG)
            ELSE
               IRPAMO(II) = IREPA
            END IF
         ENDDO
         IF (NZ.NE.4.OR.SPINFR.OR.LEVYLE) IREPA = IREPA + 1
C        for NZ.eq.4: unbarred and barred in same fermion irrep !
C        (see GMULTA for def. of fermion irreps)
         DO ISTR = 1, NASH(I)
            II = II + 1
            IF (UCIBOS) THEN
               IG = IDXU2G(IASH(I)+ISTR)
               IRPAMO(II) = IBOSYM(IG)
C              to use spin symmetry here with MULTB from GMULTSF
C              we should have had NBSYM + IBOSYM(IG) here,
C              but we take care of the spin separately in GASCIP,
C              as in non-relativistic CI programs. /hjaaj aug 03
            ELSE
               IRPAMO(II) = IREPA
            END IF
         ENDDO
      END DO
C
C     Setup pointer from active orbital to a- and b-strings
C
      II = 0
      DO J = 1, NFSYM
         DO I = 1, NASH(J)
            II = II + 1
            IPTA2O(I + IASH(J) ) = II
            IPTO2A(II) = I + IASH(J)
            IPTO2B(II) = -1
         END DO
         DO I = 1, NASH(J)
            II = II + 1
            IPTB2O(I + IASH(J) ) = II
            IPTO2A(II) = -1
            IPTO2B(II) = I + IASH(J)
         END DO
      END DO
C
C     Transfer MK2 information from dcbopt.h to dcbgascip.h
C
      MINMK2 = IOPT_MINMK2
      MAXMK2 = IOPT_MAXMK2
      MK2DEL = IOPT_MK2DEL
      MK2REF = IOPT_MK2REF
C
C
C     Print section
C     -------------
C
      IF (IPRINT .GE. 2 ) THEN
         WRITE(LUPRI,'(/A,3I3/A)')
     &      ' (GASCIP_SETUP) NREP, IRP, IRRP = ', NREP, IRP, IRRP,
     &      ' -> IRRP.lt.0 means include all symmetries'
         IF (IRRP.GT.0) WRITE(LUPRI,'(2A)')
     &      ' (GASCIP_SETUP) REPNA(IRRP)     = ', REPNA(IRRP)
         WRITE(LUPRI,'(/A)')
     &      ' Symmetry of active orbitals (unbarred, barred):'
         DO I = 1,NASHT
            IA = IPTA2O(I)
            IB = IPTB2O(I)
            IF (.NOT. UCIBOS) THEN
               WRITE(LUPRI,'(I15,2(5X,A))')
     &         I, REPNA(IRPAMO(IA)), REPNA(IRPAMO(IB))
            ELSE
               WRITE(LUPRI,'(I15,2(5X,A))')
     &         I, REP(IRPAMO(IA)), REP(IRPAMO(IB))
            END IF
         END DO
      END IF
#if GASCIP_DEBUG > 0
      write(lupri,*) '*** DEBUG OUTPUT GASCIP_SETUP:'
      write(lupri,*) 'MULTB:'
      do i = 1, 2*nrep
         write(lupri,'(i5,a,65i3)') i,' :',multb(i,0:2*nrep)
      end do
      write(lupri,*) 'I, IPTA2O, IPTB2O:'
      do i = 1, nasht
         write(lupri,*) i, ipta2o(i), iptb2o(i)
      end do
      write(lupri,*) 'I, irpamo(i), irpamo(nasht+i)'
      do i = 1, nasht
         write(lupri,*) i, irpamo(i), irpamo(nasht+i)
      end do
#endif
C
      CALL QEXIT('GASCIP_SETUP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_gendet */
      SUBROUTINE GASCIP_GENDET(NDET,IDET,ISTR,MAXDET,MAXSTR,
     &                         UCIBOS,ICIREP,IPRINT_in)
C***********************************************************************
C
C     Generate determinants.
C
C     Output:
C        NDET           : number of determinants
C XXX    IDET(2,1:NDET) : string pointers for each determinant
C
C     Input :
C        MAXDET : max value for NDET
C        MAXSTR : max value for NSTR
C        UCIBOS : use boson symmetry for spinors
C        ICIREP : spatial symmetry (if .lt. 0, include all symmetries)
C        IRPINT_in : print level
C
C      Work space
C XXX    ISTR(1:MAXSTR) : array for strings (both a- and b-strings)
C
C     Written by J. Thyssen and H. J. Aa. Jensen - Jan 2001
C     Last revision : Aug. 2015 hjaaj
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from include files:
C  dcborb.h : NAELEC,NASHT, NGSH(), ?
C
#include "dcbgascip.h"
#include "dcborb.h"
#include "pgroup.h"
#include "dgroup.h"
C
      DIMENSION NA(MAXASH_GASCIP), NB(MAXASH_GASCIP)
      INTEGER*8 IDET(2,*), IASTR, IBSTR
      INTEGER*8 ISTR(*)
      LOGICAL GASCIP_GASOK, GASCIP_SYMOK, UCIBOS, MK_symmetric
C
      INTEGER*8 MSKGAS(MXGAS), L
      INTEGER IBTAB(0:255)

      iprint = max(GASCIP_DEBUG,iprint_in)
C
C     Use GASCIP:
C
      CALL GASCIP_SETUP(UCIBOS,ICIREP,IPRINT)
C
      IF ( IPRINT .GE. 1 ) THEN
         WRITE(LUPRI,'(/A)')
     &       ' (GASCIP_GENDET) Output from generation of determinants:'
         IF (UCIBOS) WRITE(LUPRI,'(10X,A)')
     &    'Using boson spinor symmetry (exact or approximate) in the CI'
         IF (IRRP .LE. 0) THEN
            IF (UCIBOS) THEN
               WRITE(LUPRI,'(10X,A,I4,A,I3)')
     &        'Generating determinants of all spatial symmetries',
     &        ' and 2*MS interval:',MINMK2,' to',MAXMK2
            ELSE
               WRITE(LUPRI,'(10X,A,I4,A,I3)')
     &        'Generating determinants of all symmetries',
     &        ' and 2*MK interval:',MINMK2,' to',MAXMK2
            END IF
         ELSE IF (UCIBOS) THEN
            WRITE(LUPRI,'(10X,3A,I4,A,I3)')
     &        'Generating dets. with spatial symmetry ', REP(IRP-1),
     &        ' and 2*MS interval:',MINMK2,' to',MAXMK2
         ELSE
            WRITE(LUPRI,'(10X,3A,I4,A,I3)')
     &        'Generating dets. with symmetry ', REPNA(IRRP),
     &        ' and 2*MK interval:',MINMK2,' to',MAXMK2
         END IF
         WRITE(LUPRI,'(10X,A,I10)')
     &        'Max no. of determinants',MAXDET
      END IF
C
      IF (NASHT .GT. MAXASH_GASCIP) THEN
         WRITE(LUPRI,'(//A,I3,A,I3)')
     &      ' *** ERROR in GASCIP_GENDET *** Number of active'//
     &      ' orbitals (',NASHT,') must be <',MAXASH_GASCIP+1
         CALL QUIT('*** ERROR in GASCIP_GENDET ***')
      END IF
      CALL FLSHFO(LUPRI)
C
C     Setup bit table
C     ---------------
C
      CALL IBTABINI(IBTAB)
C
C     Setup GAS stuff
C     ---------------
C
      IF (NGAS_DC .GT. 0) THEN
C
         DO J = 1, NGAS_DC
            MSKGAS(J) = 0
         END DO
C
         L = 1
         DO I = 1, NFSYM
            DO J = 1, NGAS_DC
               DO K = 1, NGSH(I,J)
                  MSKGAS(J) = IOR(MSKGAS(J),L)
                  L = L * 2
               END DO
            END DO
         END DO
C
#if GASCIP_DEBUG > 1
         DO I = 1, NGAS_DC
            write(lupri,'(I5,A,B64.64)')
     &      I,' gas bit mask            : ',mskgas(i)
         END DO
#endif
C
         DO I = 2, NGAS_DC
            MSKGAS(I) = MSKGAS(I-1) + MSKGAS(I)
         END DO
C
#if GASCIP_DEBUG > 1
         DO I = 1, NGAS_DC
            write(lupri,'(I5,A,B64.64)')
     &      I,' gas bit mask accumulated: ',mskgas(i)
         END DO
#endif
C
      ELSE
         CALL QUIT(
     &   'GASCIP module called, but number of GAS spaces is zero!')
      END IF

      IF (IPRINT .GE. 0) THEN
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(A)')
     &   '   Specification of GAS spaces in GASCIP calculation'
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(A,2I4)')
     &     ' * Inactive orbitals     : ',(NISH(I),I=1,NFSYM)
         WRITE(LUPRI,'(A,2I4)')
     &     ' * Active orbitals       : ',(NASH(I),I=1,NFSYM)
         WRITE(LUPRI,'(A,I4)')
     &     ' * Active electrons      : ',NAELEC
         WRITE(LUPRI,'(A,I3,A)')
     &     ' * GAS space setup for the ',NGAS_DC,' GAS space(s) :'
         DO I = 1, NGAS_DC
            WRITE(LUPRI,'(A,I3,A,2I4)')
     &           '   - GAS space ',I,'       : ',
     &           (NGSH(J,I),J=1,NFSYM)
            WRITE(LUPRI,'(A,I3,A1,I3,A)')
     &      '    (constraints: min/max active electrons after space : ',
     &           NGASSP(1,I),'/',NGASSP(2,I),')'
         END DO
         IF (NGASSP(2,NGAS_DC) .NE. NAELEC) THEN
            CALL QUIT('*** ERROR in GASCIP_GENDET *** ' //
     &      'Inconsistency between GAS specification and NAELEC')
         END IF
      END IF
C
C     Find min and max number of A & B electrons
C     ------------------------------------------
C     NB! NAELEC is number of active electrons, not number of A electrons
C
      MAXAEL = (MAXMK2 + NAELEC) / 2
      MAXAEL = MIN(MAXAEL,NASHT)
C
      MAXBEL = (NAELEC - MINMK2) / 2
      MAXBEL = MIN(MAXBEL,NASHT)
C
      MINAEL = NAELEC - MAXBEL
      MINBEL = NAELEC - MAXAEL
C
      MINOCC = MIN(MINAEL,MINBEL)
      MAXOCC = MAX(MAXAEL,MAXBEL)

      IF ( MOD ( NAELEC, 2 ) .EQ. 1  .AND. IRRP .gt. 0) THEN
         ! For symmetry and odd number of electrons we need to check
         ! all MK. For e.g. MAXMK2 = 3, we will get GASCIP_SYMOK true
         ! for MK = 3/2 and MK = -1/2, we will thus not get all if
         ! we only check positive MK2 values !!
         MK_symmetric = .false.
      ELSE
         MK_symmetric = (MAXMK2 + MINMK2 .EQ. 0)
      END IF

C
C
C     Generate determinants:
C     -----------------------
C
#if GASCIP_DEBUG > 1
      write(lupri,*) 'minmk2, maxmk2',minmk2,maxmk2
      write(lupri,*) 'NAELEC, NASHT ',NAELEC, NASHT
      write(lupri,*) 'MK_symmetric  ',MK_symmetric
#endif
      IEL_STR(:) = -123456789
      NEL_STR(:) = -123456789
      NEL_STR_T = 0
      DO IEL = MINOCC,MAXOCC
         IEL_STR(IEL) = NEL_STR_T
         J1 = NEL_STR_T + 1
         MAXSTR_I = MAXSTR - NEL_STR_T
         CALL GASCIP_GENSTR(IEL,ISTR(J1),NEL_STR(IEL),MSKGAS,MAXSTR_I)
#if GASCIP_DEBUG > 1
         write(lupri,*) 'GASCIP_GENSTR: IEL, NEL_STR',IEL,NEL_STR(IEL)
#endif
         NEL_STR_T = NEL_STR_T + NEL_STR(IEL)
      END DO

      NASTR_T = 0
      NBSTR_T = 0
      NDET    = 0
      IMK2(:,:) = -1234567890
      DO 10 IAEL = MINAEL, MAXAEL
         IBEL = NAELEC - IAEL
         MK2 = IAEL - IBEL
         IF (MK2 .GT. 0 .AND. MK_symmetric) GO TO 11
         IASTR_OFF = IEL_STR(IAEL)
         IBSTR_OFF = IEL_STR(IBEL)
         IMK2(MK2,1) = IASTR_OFF
         IMK2(MK2,2) = IBSTR_OFF
         IMK2(MK2,3) = NDET
C
         NASTR = NEL_STR(IAEL)
         NBSTR = NEL_STR(IBEL)

#if GASCIP_DEBUG > 1
         write(lupri,'(/A,3I10)') 'A-, B-electrons, MK2',IAEL,IBEL,MK2
         write(lupri,'(A,3I20)') 'NASTR,NBSTR,NDET',NASTR,NBSTR,NDET
#endif
         IF (NASTR .EQ. 0) GOTO 10
         IF (NBSTR .EQ. 0) GOTO 10
C
C HJ TODO : the next can easily be optimized by not redoing the IA
C           for each IB in GASCIP_GASOK and GASCIP_SYMOK /aug 03
C
         DO IA = 1, NASTR
            JASTR = IASTR_OFF + IA
            DO IB = 1, NBSTR
               JBSTR = IBSTR_OFF + IB
               IF (GASCIP_GASOK(ISTR(JASTR),ISTR(JBSTR),MSKGAS,IBTAB))
     &            THEN
                  IF (GASCIP_SYMOK(ISTR(JASTR),ISTR(JBSTR),UCIBOS)) THEN
                     NDET = NDET + 1
                     IF (NDET .GE. MAXDET) THEN
                        WRITE(LUPRI,'(/A,I12/A)')
     &                     ' *** ERROR in GASCIP_GENDET *** '//
     &                       'Maximum number of determinants reached!',
     &                       MAXDET,
     &                     'Allocate more memory!'
                        CALL QUIT('*** ERROR in GASCIP_GENDET ***')
                     END IF
                     IDET(1,NDET) = ISTR(JASTR)
                     IDET(2,NDET) = ISTR(JBSTR)
                  END IF
               END IF
            END DO
         END DO
#if GASCIP_DEBUG > 1
         write(lupri,'(A,I60)') 'after check NDET',NDET
         call flshfo(lupri)
#endif
         NMK2(MK2,1) = NASTR
         NMK2(MK2,2) = NBSTR
         NMK2(MK2,3) = NDET - IMK2(MK2,3)
         NASTR_T = NASTR_T + NASTR
         NBSTR_T = NBSTR_T + NBSTR
 10   CONTINUE

 11   CONTINUE
      IF (MK_symmetric) THEN
      ! generate MK2>0 determinants
         DO MK2 = 1, MAXMK2
            ! swap what is A-str and B-str
            IMK2(MK2,1) = IMK2(-MK2,2)
            IMK2(MK2,2) = IMK2(-MK2,1)
            NMK2(MK2,1) = NMK2(-MK2,2)
            NMK2(MK2,2) = NMK2(-MK2,1)
            NASTR_T     = NASTR_T + NMK2(MK2,1)
            NBSTR_T     = NBSTR_T + NMK2(MK2,2)

            NDET_MK2    = NMK2(-MK2,3)

            IMK2(MK2,3) = NDET
            NMK2(MK2,3) = NDET_MK2

            J = IMK2(-MK2,3)
            DO I = 1, NDET_MK2
               IDET(1,NDET+I) = IDET(2,J+I)
               IDET(2,NDET+I) = IDET(1,J+I)
            END DO
            NDET        = NDET + NDET_MK2

         END DO
      END IF
C
C
C     Write statistics about determinants
C
      IF (IPRINT .GE. 0) THEN
         WRITE(LUPRI,9100)
         DO MK2 = MINMK2, MAXMK2
            IF (NMK2(MK2,3) .GT. 0) THEN
               IF (MOD(MK2,2) .EQ. 0) THEN
                  WRITE(LUPRI,9101) MK2/2,NMK2(MK2,3),
     &                              NMK2(MK2,1), NMK2(MK2,2)
               ELSE
                  WRITE(LUPRI,9102) MK2,  NMK2(MK2,3),
     &                              NMK2(MK2,1), NMK2(MK2,2)
               END IF
            END IF
         END DO
         WRITE(LUPRI,9103) NDET, NASTR_T, NBSTR_T
         CALL PRSYMB(LUPRI,'=',75,0)
         CALL FLSHFO(LUPRI)
      END IF
 9100 FORMAT(//'   Distribution of determinants wrt. MK'
     &        /6X,'MK',T28,'#dets',T37,'#A-str',T47,'#B-str'
     &        /3X,56('-'))
 9101 FORMAT(I8,T23,3I10)
 9102 FORMAT(I6,'/2',T23,3I10)
 9103 FORMAT(3X,56('-') / 5X,'all',T23,3I10)
C
C
      IF (IPRINT .GE. 20 .OR. (IPRINT .GE. 3 .AND. NDET .LE. 50) ) THEN
         WRITE(LUPRI,9001)
         DO MK2 = MINMK2, MAXMK2
            IF (NMK2(MK2,3) .GT. 0) THEN
               IF (MOD(MK2,2) .EQ. 0) THEN
                  WRITE(LUPRI,9004) NMK2(MK2,3),MK2/2
               ELSE
                  WRITE(LUPRI,9005) NMK2(MK2,3),MK2
               END IF
            END IF
C
            DO I = 1, NMK2(MK2,3)
               IASTR = IDET(1,I+IMK2(MK2,3))
               DO J=1,NASHT
                  IF ( BTEST(IASTR,J-1) ) THEN
                     NA(J)=1
                  ELSE
                     NA(J)=0
                  END IF
               END DO
               IBSTR = IDET(2,I+IMK2(MK2,3))
               DO J=1,NASHT
                  IF ( BTEST(IBSTR,J-1) ) THEN
                     NB(J)=1
                  ELSE
                     NB(J)=0
                  END IF
               END DO
               WRITE(LUPRI,9002) I+IMK2(MK2,3),(NA(J),J=1,NASHT)
               WRITE(LUPRI,9003)               (NB(J),J=1,NASHT)
            END DO
         END DO
         WRITE(LUPRI,9006)
         CALL FLSHFO(LUPRI)
      END IF
 9001 FORMAT(//1X,79('=')/' Determinant',T15,'Strings'/1X,79('='))
 9002 FORMAT(I10,T15,100I1)
 9003 FORMAT(    T15,100I1)
 9004 FORMAT(I10,' determinants with MK =',I3,' :'/1X,39('-'))
 9005 FORMAT(I10,' determinants with MK =',I3,'/2 :'/1X,41('-'))
 9006 FORMAT(1X,79('='))
C
C
      RETURN
 1010 FORMAT(I6,I10,32I2)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_gasok */
      LOGICAL FUNCTION GASCIP_GASOK(IA,IB,MSKGAS,IBTAB)
C***********************************************************************
C
C     Check that determinant is allowed for GASSCF calculations.
C
C     Input :
C       IA,IB  : a- and b-string of determinant to be tested
C       MSKGAS : mask for counting occupation in each GAS space
C       IBTAB  : table of # of 1-bits for table lookup in FNBITS
C
C     Output:
C       function value true or false
C
C     Written by J. Thyssen - Feb 11 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from include files:
C  dcborb.h : NGASSP(),?

#include "dcbgascip.h"
#include "dcborb.h"
C
      INTEGER*8 MSKGAS(*), IA, IB, IGA, IGB
      DIMENSION IBTAB(*)
C
C
      GASCIP_GASOK = .TRUE.
      IF (NGAS_DC .EQ. 0) RETURN
C
      DO IGAS = 1, NGAS_DC
         IGA = IAND(MSKGAS(IGAS),IA)
         IGB = IAND(MSKGAS(IGAS),IB)
         CALL FNBITS(NA,IGA,IBTAB)
         CALL FNBITS(NB,IGB,IBTAB)
#if GASCIP_DEBUG > 10
         write(lupri,*) 'GASCIP_GASOK: fnbits igas,na,nb',igas,na,nb
         WRITE (LUPRI,'(10X,B64/10X,B64)') IGA,IGB
#endif
         NG = NA + NB
         IF ( NG .LT. NGASSP(1,IGAS) .OR. NG .GT. NGASSP(2,IGAS)) THEN
            GASCIP_GASOK = .FALSE.
            RETURN
         END IF
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_krlist */
      SUBROUTINE GASCIP_KRLIST(NDET,IDET,KRLIST,IPRINT)
C***********************************************************************
C
C     For each determinant find its Kramers partner.
C
C     Input :
C        NDET   : number of determinants
C        IDET   : determinants
C        IPRINT : print level
C
C     Output:
C        KRLIST(i) : Kramers partner determinant to determinant no. i
C
C     Written by J. Thyssen - Jan 23 2001
C     Revised Sep 2003 hjaaj (was with old, obsolete IDET format)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      INTEGER*8 IDET(2,*), IDA, IDB
      DIMENSION KRLIST(*)
C
      CALL IZERO(KRLIST,NDET)
C
      IERROR = 0
      DO 10 I = 1, NDET
         IF (KRLIST(I) .NE. 0) GOTO 10
         IDA = IDET(1,I)
         IDB = IDET(2,I)
         DO 20 J = NDET, I, -1
            IF (KRLIST(J) .NE. 0) GOTO 20
            IF (IDA .EQ. IDET(2,J) .AND. IDB .EQ. IDET(1,J)) THEN
               KRLIST(I) = J
               KRLIST(J) = I
               GOTO 30
            END IF
 20      CONTINUE
         WRITE(LUPRI,*) 'ERROR: No Kramers partner found for det.',I
         WRITE(LUPRI,'(A,B64)') '   a-string: ',IDET(1,I),
     &                          '   b-string: ',IDET(2,I)
         IERROR = IERROR + 1
 30      CONTINUE
 10   CONTINUE
C
      IF (IERROR .GT. 0)
     &   CALL QUIT('ERROR: not all det.s have a Kramers partner')
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('Output from GASCIP_KRLIST',-1)
         WRITE(LUPRI,'(T10,A,T30,A)') 'Determinant','Kramers partner'
         DO I = 1, NDET
            WRITE(LUPRI,'(T10,I9,T30,I9)') I,KRLIST(I)
         END DO
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ngascip_nbits */
      INTEGER FUNCTION NGASCIP_NBITS(JDET,IBTAB)
C***********************************************************************
C
C     Find number of bits in integer JDET.
C
C     Written by J. Thyssen - Jan 2 2001
C     Last revision : S. Knecht - Nov 17 2007
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgascip.h"
#include "dcborb.h"
C
      INTEGER   IBTAB(0:255)
C
C     17-11-2007-SK: INTEGER*1 is an illegal value on 64-Bit SX machines
C     Is this solution still working?
#if defined (SYS_SX64)
!     INTEGER*4 ID1(4)
      call qenter('NGASCIP_NBITS')
      call quit('ERROR: INTEGER*1 not defined on 64-bit SX machines')
#else
      INTEGER*8 JDET, ID8
      INTEGER*1 ID1(8)
      EQUIVALENCE (ID1(1),ID8)
#endif /* SYS_SX64 */
C
C     Anyone who can find with a faster routine?
C
      ID8 = JDET
C     write(lupri,*) 'jdet',jdet,' and id8',i84
C     write(lupri,*) 'id1',id1
      NGASCIP_NBITS
     &   = IBTAB(ID1(1)) + IBTAB(ID1(2)) + IBTAB(ID1(3)) + IBTAB(ID1(4))
     &   + IBTAB(ID1(5)) + IBTAB(ID1(6)) + IBTAB(ID1(7)) + IBTAB(ID1(8))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_symok */
      LOGICAL FUNCTION GASCIP_SYMOK(IA,IB,UCIBOS)
C***********************************************************************
C
C     Find symmetry of determinant (IA,IB).
C
C     Written by J. Thyssen - Jan 2 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      INTEGER*8 IA, IB, ID
      LOGICAL   UCIBOS
C
C dcbgascip.h : IRRP, *
C    dcborb.h : NASHT
#include "dcbgascip.h"
#include "dcborb.h"
C
      IF (IRP .LE. 0) THEN
C        ... include all symmetries
         GASCIP_SYMOK = .TRUE.
      ELSE IF (.NOT. UCIBOS) THEN
C
C        The normal symmetry scheme.
C
         JRRP = 0
         ID = IA
         DO I=1,NASHT
            IF (MOD(ID,2_8).EQ.1) THEN
               JRRP = MULTB(IRPAMO(IPTA2O(I)),JRRP)
            ENDIF
            ID = ID/2
         END DO
         ID = IB
         DO I=1,NASHT
            IF (MOD(ID,2_8).EQ.1) THEN
               JRRP = MULTB(IRPAMO(IPTB2O(I)),JRRP)
            ENDIF
            ID = ID/2
         END DO
         GASCIP_SYMOK = (JRRP .EQ. IRRP)
      ELSE
C
C        Boson symmetry (spin free or Levy-Leblond).
C        (do not use IRRP and MULTB which includes spin-symmetry -- see GMULTSF --
C         as we treat the spin Ms separately as in non-rel codes /hjaaj aug 03)
C
         JRRP = 0
         ID = IA
         DO I = 1, NASHT
            IF ( MOD(ID,2_8) .EQ. 1) THEN
               JRRP = IEOR(IRPAMO(IPTA2O(I)),JRRP)
            END IF
            ID = ID / 2
         END DO
         ID = IB
         DO I = 1, NASHT
            IF ( MOD(ID,2_8) .EQ. 1) THEN
               JRRP = IEOR(IRPAMO(IPTB2O(I)),JRRP)
            END IF
            ID = ID / 2
         END DO
         GASCIP_SYMOK = (JRRP .EQ. IRP-1)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_startvec */
      SUBROUTINE GASCIP_STARTVEC(IDET,CMO,CREF,ECI,NDET,I_STATE,
     &                           THR_PCI,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate CI start vectors from GASCIP by diagonalizing the
C     full CI matrix.
C
C     Input : CMO - MO coefficents
C
C     Output: CREF - CI vector
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      INTEGER*8 IA, IB, ID
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dgroup.h"
C
      INTEGER*8 IDET(2,*)
      DIMENSION WORK(*), CREF(NDET,NZ_in_CI), CMO(*)
C
      CALL QENTER('GASCIP_STARTVEC')
C
      KFRSAV = KFREE
      N2DET = NDET*NDET
C
      CALL MEMGET2('REAL','EIGVAL',KEIGVL,NDET,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EIGVEC',KEIGVC,N2DET*NZ_in_CI,
     &   WORK,KFREE,LFREE)
C
C
C     Do full CI
C
      CALL GASCIP_FULLCI(NDET,IDET,CMO,WORK(KEIGVL),WORK(KEIGVC),
     &                  ECORE,IPRINT,WORK,KFREE,LFREE)
C
C     Get CI vector:
C     --------------
C
C     The CI vector is stored in columns; the lowest eigenstate is in
c     the first column.
C
      IF (I_STATE .LE. 0 .OR. I_STATE .GT. NDET) THEN
         WRITE (LUPRI,*) 'ERROR: I_STATE .LE. 0 .OR. I_STATE .GT. NDET'
         WRITE (LUPRI,*) 'ERROR: I_STATE , NDET =',I_STATE,NDET
         CALL QUIT('Non-valid reference state specification')
      END IF
      IOFF = KEIGVC - 1 + (I_STATE-1) * NDET
      IF ( NZ .EQ. 1 ) THEN
         DO I = 1,NDET
            CREF(I,1) = WORK(IOFF + I)
         END DO
      ELSE
         DO I = 1,NDET
            CREF(I,1) = WORK(IOFF + I )
            CREF(I,2) = WORK((IOFF + N2DET) + I )
         END DO
      END IF
C
C DEBUG DEBUG DEBUG
c      call dzero(cref,ndetq)
c      if (nz .eq. 1) then
c      cref(2,1) = 1.0d0
c      cref(2,2) = 1.0d0/sqrt(2.0d0)
c      else
c         CREF(2,2) = 1.0D00
c      endif
c      DO I = 2, NDET
c         CREF(I,1) = 0.0D00
c      END DO
c      IF ( NDET .GE. 2) THEN
c         CREF(1,1) = 0.8
c         CREF(2,1) = 0.6
c      END IF
c      CREF(1,1) = 0.8D00
c      IF ( NDET .GE. 2) CREF(2,1) = 0.6D00
c      DO I = 3, NDET
c         CREF(I,1) = 0.0D00
c      END DO
C
C
C     Get CI energy
C
      ECI = WORK(KEIGVL + I_STATE - 1)
C
      WRITE(LUPRI,'(/A,F20.8,2(/,A,F20.8))')
     &     ' (GASCIP_STARTVEC) Core energy    =',ECORE,
     &     '                   Active energy  =',ECI-ECORE,
     &     '                   Total energy   =',ECI
C
C     Check norm of CI vector
C
      CALL RCINORM(CREF,IPRINT)
C
      IF (IPRINT .GE. 5 ) THEN
         WRITE(LUPRI,'(A)') ' (GASCIP_STARTVEC) CI start vector'
         CALL RPRCI(CREF,NDET,NZ,'GASCIP', THR_PCI,LUPRI)
      END IF
C
      CALL MEMREL('GASCIP_STARTVEC',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GASCIP_STARTVEC')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_fullci */
      SUBROUTINE GASCIP_FULLCI(NDET,IDET,CMO,EIGVAL,EIGVEC,ECORE,
     &                        IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate CI start vectors from GASCIP by diagonalizing the
C     full CI matrix.
C
C     Input : IDET   - determinant information
C             CMO    - MO coefficents
C
C     Output: EIGVAL - CI eigenvalues
C             EIGVEC - CI eigenvectors
C             ECORE  - core energy
C
C     Extracted from GASCIP_STARTVEC July 2002 /hjaaj.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION WORK(*), EIGVAL(*), EIGVEC(*), CMO(*)
      INTEGER*8 IDET(2,*)
C
      CALL QENTER('GASCIP_FULLCI')
      CALL GETTIM(CPU0,WALL0)
C
      KFRSAV = KFREE
      N2DET = NDET * NDET
      N2DETQ = N2DET * NZ_in_CI
C
C     Calculate core Hamiltonian:
C     ---------------------------
C
C     The core Hamiltonian is the active-active part of FC,
C     a.k.a. FCAC.
C
!     write(lupri,*) 'NASHT, NZ_in_CI,N2ASHXQ,N2ORBXQ',
!    &     NASHT, NZ_in_CI,N2ASHXQ,N2ORBXQ
      CALL MEMGET2('REAL','FCACM',KFCACM,(2*NASHT)*(2*NASHT)*NZ_in_CI,
     &     WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','FCAC',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','FC',  KFC,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Calculate FCmo (TOFILE false: do not write to LUKRM3) :
C
      CALL rGETFC(.FALSE.,CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
C     Get active-active part of FCmo:
C
      CALL RGETAC(WORK(KFC),WORK(KFCAC),IPRINT)
C
C     Transform quaternion FCACmo to molfdir type FCAC.
C
      CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPRINT)
C
      CALL MEMREL('GASCIP_FULLCI after rGETFC',WORK,1,KFCAC,KFREE,LFREE)
C
C     Read two-electron integrals:
C     ----------------------------
C
      CALL MEMGET2('REAL','MUUUU',KMUUUU,(2*NASHT)**4*NZ_in_CI,
     &     WORK,KFREE,LFREE)
C
      LH2AC = NASHT*NASHT*NNASHX*NZ*3
      CALL MEMGET2('REAL','H2AC',KH2AC,LH2AC,WORK,KFREE,LFREE)
C
C     read integrals from 4IND*: out --> work(kh2ac)
C     ==============================================

      call memget2('INTE','IBEIG',kibeig,norbt,work,kfree,lfree)

      if(spinfr.or.levyle)then
        call ireakrmc(lukrmc,'IBEIG   ',work(kibeig),norbt)
      else if(linear)then
        call ireakrmc(lukrmc,'MJVEC   ',work(kibeig),norbt)
      else
        call izero(work(kibeig),norbt)
      end if

      call rgeth2(dummy,work(kh2ac),dummy,work(kibeig),.false.,
     &            .true.,.true.,work(kfree),lfree)
C
C     Transform integrals to Molfdir format:
C
C     remove potential "garbage" from imaginary part...
      IF (NZ_in_CI .gt. 1)
     &   CALL DZERO(WORK(KMUUUU+(2*NASHT)**4),(2*NASHT)**4)

      if(naelec > 1)then
        CALL DNZ32M(WORK(KH2AC),WORK(KMUUUU),IPRINT)
      else
        CALL DZERO(WORK(KMUUUU),(2*NASHT)**4)
      end if

      CALL MEMREL('GASCIP_FULLCI.H2AC',WORK,1,KH2AC,KFREE,LFREE)
C
C     Allocate memory for matrices:
C
      CALL MEMGET2('REAL','HMAT',KHR,N2DETQ,WORK,KFREE,LFREE)
      KHI = KHR + N2DET
C
C     Construct Hamiltonian:
C     ----------------------
C
      CALL GETTIM(CPU1,WALL1)
      CALL GASCIP_MAKEH(.TRUE.,1,NDET,IDET,WORK(KHR),WORK(KHI),
     &                  WORK(KFCACM),WORK(KMUUUU))
      CALL GETTIM(CPU2,WALL2)
      WRITE (LUPRI,'(/A,2F20.2)')
     &  "CPU and WALL times for making CI matrix:",CPU2-CPU1,WALL2-WALL1
C
C
C     Diagonalize Hamiltonian:
C     ------------------------
C
      IF ( IPRINT .GE. 20 ) THEN
         WRITE(LUPRI,'(A)') ' (GASCIP_FULLCI) CI matrix - real part:'
         CALL OUTPUT(WORK(KHR),1,NDET,1,NDET,NDET,NDET,-1,LUPRI)
        IF (NZ_in_CI .GT. 1) THEN
         WRITE(LUPRI,'(A)') ' (GASCIP_FULLCI) CI matrix - imag part:'
         CALL OUTPUT(WORK(KHI),1,NDET,1,NDET,NDET,NDET,-1,LUPRI)
        END IF
      END IF
      CALL FLSHFO(LUPRI)
C
C
      CALL GETTIM(CPU1,WALL1)
      IF (NZ_in_CI .EQ. 1) THEN
         IJOB   = 1
         IORDER = 1
         IPACK  = 0
         CALL RSJACO(NDET,NDET,NDET,WORK(KHR),EIGVAL,
     &        IJOB,IORDER,IPACK,EIGVEC)
      ELSE
         MATZ = 1
         CALL QDIAG(2,NDET,WORK(KHR),NDET,NDET,
     &        EIGVAL,MATZ,EIGVEC,NDET,NDET,
     &        WORK(KFREE),LFREE,IERR)
         IF (IERR .NE. 0) THEN
            WRITE(LUPRI,'(/2A,I4)')
     &           ' *** ERROR in GASCIP_FULLCI ***: ',
     &           'QDIAG failed with error code ',IERR
            CALL QUIT('*** ERROR in GASCIP_FULLCI *** QDIAG failed')
         END IF
      END IF
      CALL GETTIM(CPU2,WALL2)
      WRITE (LUPRI,'(/A,2F20.2)')
     &   "CPU and WALL times for diagonalization:",CPU2-CPU1,WALL2-WALL1
C
C     Add core energy
C
      DO I = 1,NDET
         EIGVAL(I) = EIGVAL(I) + ECORE
      END DO
C
      IF ( IPRINT .GE. 3 ) THEN
         IF ( IPRINT .ge. 5 .OR. nDET.le.100) THEN
            WRITE(LUPRI,'(/A)') ' (GASCIP_FULLCI) CI eigenvalues:'
            n_eig = nDET
         ELSE
            n_eig = nDET
            WRITE(LUPRI,'(/A,I7,A,I7)') ' (GASCIP_FULLCI) First',
     &         n_eig,' CI eigenvalues out of',nDET
         END IF
         WRITE(LUPRI,'(5(I10,A,F13.6))') (I, ' :', EIGVAL(I), I=1,n_eig)
         EIGVAL_AVERAGE = DSUM(NDET,EIGVAL,1) / DBLE(NDET)
         WRITE(LUPRI,'(/A,F20.8)') ' CI average energy:',EIGVAL_AVERAGE
       IF ( IPRINT .GE. 20 .OR. NDET .LE. 12) THEN
         n_eig = min(nDET, IPRINT)
         WRITE(LUPRI,'(/A)')
     &   ' (GASCIP_FULLCI) first CI eigenvectors (real part)'
         CALL OUTPUT(EIGVEC,1,NDET,1,n_eig,
     &        NDET,NDET,-1,LUPRI)
         IF (NZ_in_CI .GT. 1) THEN
            WRITE(LUPRI,'(/A)')
     &      ' (GASCIP_FULLCI) first CI eigenvectors (imag part)'
            CALL OUTPUT(EIGVEC(1+N2DET),1,NDET,1,n_eig,
     &           NDET,NDET,-1,LUPRI)
         END IF
       END IF
      END IF
C
      CALL MEMREL('GASCIP_FULLCI',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL GETTIM(CPU2,WALL2)
      WRITE (LUPRI,'(/A,2F20.2)')
     &   "CPU and WALL times for full CI:",CPU2-CPU0,WALL2-WALL0
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GASCIP_FULLCI')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gascip_ireps */
      SUBROUTINE GASCIP_IREPS (NZ,NFSYM,NSYMRPA,REPA)
C***********************************************************************
C
C     Rip-off from GMULTA from MOLTRA./jth-03072000
C
C     Hard-wired multiplication table for the abelian subgroups of D2h
C     Defines representation names as well.
C
C     Luuk Visscher, Mar 18 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER*4 REPA(8)
C
      IF (NFSYM.EQ.1.AND.NZ.EQ.4) THEN
         NSYMRPA = 1
         REPA(1) = '   A'
         REPA(2) = '   a'
      ELSEIF (NFSYM.EQ.1.AND.NZ.LE.2) THEN
         NSYMRPA = 2
         REPA(1) = '  1E'
         REPA(2) = '  2E'
         REPA(3) = '   a'
         REPA(4) = '   b'
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.4) THEN
         NSYMRPA = 2
         REPA(1) = '  AG'
         REPA(2) = '  AU'
         REPA(3) = '  ag'
         REPA(4) = '  au'
      ELSEIF (NFSYM.EQ.2.AND.NZ.LE.2) THEN
         NSYMRPA = 4
         REPA(1) = ' 1Eg'
         REPA(2) = ' 2Eg'
         REPA(3) = ' 1Eu'
         REPA(4) = ' 2Eu'
         REPA(5) = '  ag'
         REPA(6) = '  bg'
         REPA(7) = '  au'
         REPA(8) = '  bu'
      ELSE
         CALL QUIT('Multiplication table is not defined')
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_rmakdm */
      SUBROUTINE GASCIP_RMAKDM(CALCDV,CALCPV,NDET,IDET,CL,CR,DV,PV,
     &                         IPRINT)
C***********************************************************************
C
C     Construct one- and two-electron active density matrix
C
C     D_{pq} = < CL | X^+_{pq} | CR >
C
C     P_{pqrs} = < CL | x^+_{pq,rs} | CR >
C
C     Input:
C        CL, CR  - CI vectors
C
C     Output:
C        DV      - one electron active density matrix
C        PV      - two electron active density matrix
C
C     DV has dimension (2*NASHT)**2 * MIN(NZ,2).
C     PV has dimension (2*NASHT)**4 * MIN(NZ,2). Triangular packing
C     would be better, but whatever...
C     The code is now very straightforward.
C
C     Written by J. Thyssen - Oct 22 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "thrzer.h"
C
#include "dcbgascip.h"
#include "dgroup.h"
#include "dcborb.h"
C
      DIMENSION CL(NDET,NZ_in_CI),CR(NDET,NZ_in_CI)
      DIMENSION DV(2*NASHT,2*NASHT,2)
      DIMENSION PV(2*NASHT,2*NASHT,2*NASHT,2*NASHT,2)
      INTEGER*8 IDET(2,*)
C
      LOGICAL CALCDV, CALCPV
C
C
      CALL QENTER('GASCIP_RMAKDM')
C
      IF (.NOT. CALCDV) THEN
         WRITE(LUPRI,'(/A)')
     $        ' *** ERROR in GASCIP_RMAKDM *** ' //
     $        'CALCDV false is not implemented.'
         CALL QUIT('*** ERROR in GASCIP_RMAKDM ***')
      END IF
C
C     Calculate density matrices
C     --------------------------
C
      CALL GASCIP_MAKETDM(CALCPV,NDET,IDET,CL,CR,DV,PV,IPRINT)
C
C     Print section
C     -------------
C
      IF ( IPRINT .GE. 10 ) THEN
C
C        Print DV
C
         IF (CALCDV) THEN
            WRITE(LUPRI,'(/A)')
     &           ' (GASCIP_RMAKDM) real(DV) in Molfdir basis'
            CALL OUTPUT(DV(1,1,1),1,2*NASHT,1,2*NASHT,
     &           2*NASHT,2*NASHT,-1,LUPRI)
            IF ( NZ .GE. 2 ) THEN
               WRITE(LUPRI,'(/A)')
     &              ' (GASCIP_RMAKDM) imag(DV) in Molfdir basis'
               CALL OUTPUT(DV(1,1,2),1,2*NASHT,1,2*NASHT,
     &              2*NASHT,2*NASHT,-1,LUPRI)
            END IF
         END IF
C
C        Print PV
C
         IF (CALCPV) THEN
            THROUT = 1.0D-10
            WRITE(LUPRI,'(/A,1P,D10.2)')
     &           ' (GASCIP_RMAKDM) PV in Molfdir basis,'//
     &           ' print threshold:',THROUT
            WRITE(LUPRI,'(4X,A,T20,A,10X,A)')
     $           'Index','Real(PV)','Imag(PV)'
            DO IY = 1,2*NASHT
               DO IX = 1,2*NASHT
                  DO IV = 1,2*NASHT
                     DO IU = 1,2*NASHT
                        IF ( NZ .EQ. 1 ) THEN
                           IF(ABS(PV(IU,IV,IX,IY,1)) .GT. THROUT) THEN
                              WRITE(LUPRI,'(4I3,F20.10)')
     $                             IU,IV,IX,IY,PV(IU,IV,IX,IY,1)
                           END IF
                        ELSE
                           IF((ABS(PV(IU,IV,IX,IY,1))
     &                       + ABS(PV(IU,IV,IX,IY,2))) .GT. THROUT) THEN
                              WRITE(LUPRI,'(4I3,2F20.10)')
     $                             IU,IV,IX,IY,
     $                             PV(IU,IV,IX,IY,1),PV(IU,IV,IX,IY,2)
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END IF
C
      CALL QEXIT('GASCIP_RMAKDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_maketdn */
      SUBROUTINE GASCIP_MAKETDM(DOPV,NDET,IDET,CL,CR,DV,PV,IPRINT)
C***********************************************************************
C
C     Calculate density matrices: <CL|Pq|CR>
C
C     Input:
C        CL, CR - CI vectors
C        DOPV - calculate 2-electron density matrix
C
C     Output:
C        DV, PV - density matrices
C
C     Written by H. J. Aa. Jensen & J. Thyssen - Dec 30 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL   DOPV
      DIMENSION CL(*), CR(*), DV(*), PV(*)
      INTEGER*8 IDET(2,*)
C
      CALL QENTER('GASCIP_MAKETDM')
C
      IHRM = 1
      CALL GASCIP_CI(
     &     .FALSE.,DUMMY,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,CL,CR,DV,PV,
     &     DOPV,NDET,IDET,IHRM,IPRINT)
C     CALL GASCIP_CI(
C    &     DODIA,DIAG,
C    &     DOSIG,SVEC,BVEC,HCORE,VMUUUU,
C    &     DODEN,CL,CR,DV,PV,
C    &     DO2EL,NDET,IDET,IPRINT)
C
      CALL QEXIT('GASCIP_MAKETDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_denmat */
      SUBROUTINE GASCIP_DENMAT(DOPV,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
     &     DMRIJ,DMIIJ,DMRJI,DMIJI,DV,PV)
C***********************************************************************
C
C     Calculate contributions determinant I and J.
C
C     Input:
C        DOPV   - calculate PV
C        DOIMAG - calculate imaginary contributions
C        IT     - the number of different spinors in determinants I and J
C        IOCCI, JOCCI - the strings for determinants I and J
C        K,L,M,N - IT = 0: not used
C                  IT = 2: K and M is the orbitals where I and J differ
C                  IT = 4: K,L,M,N  ditto.
C        DMRIJ, DMRJI - real part of cl^*_{i} cr_{j} and cl^*_{j} cr_{i}
C        DMIIJ, DMIJI - imag part of ditto.
C
C     Output:
C        DV, PV: density matrices
C
C     Written by H. J. Aa. Jensen & J. Thyssen - Dec 30 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from include files:
C  dcborb.h: NAELEC, NASHT, ?

#include "dcborb.h"
#include "dcbgascip.h"
C
      LOGICAL   DOPV, DOIMAG
      DIMENSION IOCCI(2*MAXASH_GASCIP), JOCCI(2*MAXASH_GASCIP)
      DIMENSION DV(2*NASHT,2*NASHT,*)
      DIMENSION PV(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      IF (IT .EQ. 0) THEN
C
C     ======================
C     DETERMINANTS ARE EQUAL
C     ======================
C
         DO 30 II = 1, NAELEC
            I = IOCCI(II)
            DV(I,I,1) = DV(I,I,1) + DMRIJ
            IF (DOIMAG) DV(I,I,2) = DV(I,I,2) + DMIIJ
            IF (DOPV) THEN
               DO 20 JJ = II + 1, NAELEC
                  J = IOCCI(JJ)
C---------------------
C                    H = H + (I, I, J, J)
C                    H = H - (I, J, J, I)
C---------------------
                  PV(J,J,I,I,1) = PV(J,J,I,I,1) + DMRIJ
                  PV(I,J,J,I,1) = PV(I,J,J,I,1) - DMRIJ
C---------------------
                  IF (DOIMAG) THEN
                     PV(J,J,I,I,2) = PV(J,J,I,I,2) + DMIIJ
                     PV(I,J,J,I,2) = PV(I,J,J,I,2) - DMIIJ
                  END IF
 20            CONTINUE
            END IF
 30      CONTINUE
C
      ELSE IF (IT.EQ.2) THEN
C
C        ===================================
C        DETERMINANTS DIFFER BY ONE FUNCTION K -> M
C        ===================================
C
         IPHASE = 0
         DO IELEC = 1,NAELEC
            IF (IOCCI(IELEC) .EQ. K) IPHASE = IPHASE + IELEC
            IF (JOCCI(IELEC) .EQ. M) IPHASE = IPHASE + IELEC
         END DO
         IPHASE = (-1)**IPHASE
C
         DMRIJ1 = IPHASE * DMRIJ
         DMRJI1 = IPHASE * DMRJI
         DV(K,M,1) = DV(K,M,1) + DMRIJ1
         DV(M,K,1) = DV(M,K,1) + DMRJI1
         IF (DOIMAG) THEN
            DMIIJ1 = IPHASE * DMIIJ
            DMIJI1 = IPHASE * DMIJI
            DV(K,M,2) = DV(K,M,2) + DMIIJ1
            DV(M,K,2) = DV(M,K,2) + DMIJI1
         END IF
C
         IF (DOPV) THEN
            DO IELEC = 1, NAELEC
               I = IOCCI(IELEC)
               IF (I .NE. K) THEN
C---------------------
C                 H = H + (I, I, K, M)
C                 H = H - (I, M, K, I)
C---------------------
                  PV(I,I,K,M,1) = PV(I,I,K,M,1) + DMRIJ1
                  PV(I,M,K,I,1) = PV(I,M,K,I,1) - DMRIJ1
                  PV(I,I,M,K,1) = PV(I,I,M,K,1) + DMRJI1
                  PV(I,K,M,I,1) = PV(I,K,M,I,1) - DMRJI1
C---------------------
                  IF (DOIMAG) THEN
                     PV(I,I,K,M,2) = PV(I,I,K,M,2) + DMIIJ1
                     PV(I,M,K,I,2) = PV(I,M,K,I,2) - DMIIJ1
                     PV(I,I,M,K,2) = PV(I,I,M,K,2) + DMIJI1
                     PV(I,K,M,I,2) = PV(I,K,M,I,2) - DMIJI1
                  END IF
               END IF
            END DO
         END IF
C
      ELSE IF(IT.EQ.4 .AND. DOPV) THEN
C
C        ====================================
C        DETERMINANTS DIFFER BY TWO FUNCTIONS  K,L -> M,N
C        ====================================
C
         IPHASE = 0
         DO IELEC = 1,NAELEC
            IF (IOCCI(IELEC) .EQ. K) IPHASE = IPHASE + IELEC
            IF (IOCCI(IELEC) .EQ. L) IPHASE = IPHASE + IELEC
            IF (JOCCI(IELEC) .EQ. M) IPHASE = IPHASE + IELEC
            IF (JOCCI(IELEC) .EQ. N) IPHASE = IPHASE + IELEC
         END DO
         IPHASE = (-1)**IPHASE
C
         DMRIJ1 = IPHASE * DMRIJ
         DMRJI1 = IPHASE * DMRJI
         IF (DOIMAG) THEN
            DMIIJ1 = IPHASE * DMIIJ
            DMIJI1 = IPHASE * DMIJI
         END IF
C
C---------------------
C        H = H + (K, M, L, N)
C---------------------
         PV(K,M,L,N,1) = PV(K,M,L,N,1) + DMRIJ1
         PV(K,N,L,M,1) = PV(K,N,L,M,1) - DMRIJ1
         PV(M,K,N,L,1) = PV(M,K,N,L,1) + DMRJI1
         PV(M,L,N,K,1) = PV(M,L,N,K,1) - DMRJI1
C---------------------
C        H = H - (K, N, L, M)
C---------------------
         IF (DOIMAG) THEN
            PV(K,N,L,M,2) = PV(K,N,L,M,2) + DMIIJ1
            PV(K,M,L,N,2) = PV(K,M,L,N,2) - DMIIJ1
            PV(M,K,N,L,2) = PV(M,K,N,L,2) + DMIJI1
            PV(M,L,N,K,2) = PV(M,L,N,K,2) - DMIJI1
         END IF
      ELSE
C
C     =============================================
C     DETERMINANTS DIFFER BY MORE THEN TWO FUNCTIONS
C     ==============================================
C
         CALL QUIT('GASCIP_DENMAT called with IT .gt. 4 !')
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_makeh */
      SUBROUTINE GASCIP_MAKEH(DO2EL,IHRM,NDET,IDET,HR,HI,HCORE,VMUUUU)
C***********************************************************************
C
C     Construct CI Hamiltonian explicitly.
C
C     Input:
C        DO2EL  - true: both one and two electron integrals
C                 false: only one electron integrals
C        IHRM   - .gt. 0: Hermitian operator
C                 .lt. 0: Anti-hermitian operator
C                 .eq. 0: non-Hermitian operator
C                         (or debug to calculate full CI Hamiltonian)
C        HCORE  - one electron integrals
C        VMUUUU - two electron integrals.
C
C     Output:
C        HR - real part of CI Hamiltonian
C        HI - imag part of CI Hamiltonian
C
C     Written by H. J. Aa. Jensen & J. Thyssen - Dec 30 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "thrzer.h"
C
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION HR(NDET,NDET),HI(NDET,NDET)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
      DIMENSION HCORE(*)
      INTEGER*8 IDET(2,*), IDA, IDB, JDA, JDB, IDAOLD, JDAOLD
      INTEGER*8 IDAV, IDBV, IAI, IBI, IAJ, IBJ
      LOGICAL   DO2EL,DOIMAG
      DIMENSION IBTAB(0:255)
C
      DIMENSION IOCCI(2*MAXASH_GASCIP), JOCCI(2*MAXASH_GASCIP)
      LOGICAL   HFULL
C
      CALL IBTABINI(IBTAB)
      HFULL = IHRM .EQ. 0
      DOIMAG  = NZ.GT.1
      CALL DZERO(HR,NDET*NDET)
      IF (DOIMAG) CALL DZERO(HI,NDET*NDET)
C
      IF (DO2EL) THEN
C        ... one and two electron terms
         ITMAX = 4
      ELSE
C        ... only one electron terms
         ITMAX = 2
      END IF
#if GASCIP_DEBUG > 30
      write (lupri,'(/A)')
     &   '  GASCIP_MAKEH: Test output of IDET(2,NDET):)'
      DO I = 1, NDET
         WRITE (LUPRI,'(/I5,5X,B64/10X,B64)') I, IDET(1,I), IDET(2,I)
      END DO
#endif
C
      IDAOLD = -1

      DO I = 1, NDET
C
C        Use that outer loop in gascip_gendet is over A-strings,
C        inner loop is over B-strings
C
         IDA = IDET(1,I)
         IF (IDA .NE. IDAOLD) THEN
            IDAOLD = IDA
            IELECA = 0
            DO K = 1, NASHT
               IF ( BTEST(IDA,K-1) ) THEN
                  IELECA = IELECA + 1
                  IOCCI(IELECA) = IPTA2O(K)
               END IF
            END DO
         END IF
C
         IDB = IDET(2,I)
         IELEC = IELECA
         DO K = 1, NASHT
            IF ( BTEST(IDB,K-1) ) THEN
               IELEC = IELEC + 1
               IOCCI(IELEC) = IPTB2O(K)
            END IF
         END DO
         IF (IELEC .NE. NAELEC) THEN
            WRITE(LUPRI,*) 'GASCIP_MAKEH: IELEC.ne.NAELEC',IELEC,NAELEC
            WRITE(LUPRI,'(A,B64)') 'IDA ',IDA
            WRITE(LUPRI,'(A,B64)') 'IDB ',IDB
            CALL QUIT('GASCIP_MAKEH: IELEC.ne.NAELEC')
         END IF
C
         IF (HFULL) THEN
            JEND = NDET
         ELSE
chj aug 2002: code generalized to non-hermitian operators
C
C        Diagonal elements!
C
            CALL GASCIP_HMAT(DO2EL,DOIMAG,0,IOCCI,IOCCI,0,0,0,0,
     &           HR(I,I),HI(I,I),HCORE,VMUUUU)
            JEND = I - 1
         END IF
         JDAOLD = -1
         DO J = 1, JEND
C
C           Use that outer loop in gascip_gendet is over A-strings,
C           inner loop is over B-strings
C
            JDA = IDET(1,J)
            IF (JDA .NE. JDAOLD) THEN
               JDAOLD = JDA
               IDAV = IEOR(IDA, JDA)
               CALL FNBITS(NDAV,IDAV,IBTAB)
C              Can we get a contribution from <I| H |J> ?
               IF (NDAV .LE. ITMAX) THEN
                  JELECA = 0
                  IAI = IAND(IDAV, IDA) ! a-orbs not in jdet
                  KA = 0
                  IAJ = IAND(IDAV, JDA) ! a-orbs not in idet
                  MA = 0
                  DO II = 1, NASHT
                     IF ( BTEST(JDA, II-1) ) THEN
                        JELECA = JELECA + 1
                        JOCCI(JELECA) = IPTA2O(II)
                     END IF
                     IF ( BTEST(IAI, II-1) ) THEN
                        L  = KA
                        KA = IPTA2O(II)
                     END IF
                     IF ( BTEST(IAJ, II-1) ) THEN
                        N  = MA
                        MA = IPTA2O(II)
                     END IF
                  END DO
               END IF
            END IF
C
C           Any reason to check the B-string difference?
C
            IF (NDAV .LE. ITMAX) THEN
C
C           Yes, do it.
C
               JDB = IDET(2,J)
               IDBV = IEOR(IDB, JDB)
               CALL FNBITS(NDBV,IDBV,IBTAB)
C
               IT = NDAV + NDBV
C
C           Can this give a non-zero contribution ?
C
            IF ( IT .LE. ITMAX ) THEN
C
C           Yes, calculate the contribution.
C
               IELEC = JELECA
               DO K = 1, NASHT
                  IF ( BTEST(JDB, K-1) ) THEN
                     IELEC = IELEC + 1
                     JOCCI(IELEC) = IPTB2O(K)
                  END IF
               END DO
C
               IBI = IAND(IDBV, IDB)
               IBJ = IAND(IDBV, JDB)
               K = KA
               M = MA
               DO II = 1, NASHT
                  IF ( BTEST(IBI, II-1) ) THEN
                     L = K
                     K = IPTB2O(II)
                  END IF
                  IF ( BTEST(IBJ, II-1) ) THEN
                     N = M
                     M = IPTB2O(II)
                  END IF
               END DO
C
#if GASCIP_DEBUG > 30
                  write (lupri,*)' dets i, j',i,j
                  write (lupri,*)' do2el, doimag, it',
     &                             do2el, doimag, it
                  write (lupri,*) 'k,l,m,n',k,l,m,n
                  write (lupri,*) 'iocci',iocci(1:2*naelec)
                  write (lupri,*) 'jocci',jocci(1:2*naelec)
                  write (lupri,*) 'b HR(i,j), HI(i,j)',HR(i,j),HI(i,j)
#endif
               CALL GASCIP_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
     &              HR(I,J),HI(I,J),HCORE,VMUUUU)
#if GASCIP_DEBUG > 30
               if (i .eq. j) then
                  write (lupri,*) 'a HR(i,j), HI(i,j)',HR(i,j),HI(i,j)
               end if
#endif
               IF( IHRM .GT. 0 ) THEN
                  HR (J,I) =  HR (I, J)
                  IF (DOIMAG) HI (J,I) = -HI (I, J)
               ELSE IF( IHRM .LT. 0 ) THEN
                  HR (J,I) = -HR (I, J)
                  IF (DOIMAG) HI (J,I) =  HI (I, J)
               END IF
            END IF
C           ... for IF ( IT .LE. ITMAX ) THEN
            END IF
C           ... for IF (NDAV .LE. ITMAX) THEN
         ENDDO
      ENDDO

C
c     write(LUPRI,*) 'real part of HMAT in GASCIP_MAKEH ='
c     call output(hr,1,NDET,1,NDET,NDET,NDET,-1,lupri)
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_makepmat */
      SUBROUTINE GASCIP_MAKE_P1MAT(CIVEC,N_CIVEC,NDET,IDET,P1MAT_CI,
     &                             N_P1MAT_CI,IHRM_P1,P1MAT_MO,
     &                             WORK,LWORK, IPRINT)
C***********************************************************************
C
C     Construct the property matrix between CI states from
C     P_ij = c_i^+ * sigma(P1)_j
C
C     INPUT:
C            CIVEC    :     CI eigenvectors
C            N_CIVEC  :     Number of CI eigenvectors
C            IDET     :     Determinant strings
C            P1MAT_MO :     One electron property integrals
C            IHRM_P1  :     Hermiticity of P1 when NZ_in_CI .eq. 1
C
C     OUTPUT:
C            P1MAT_CI :     Property matrix between the first
C                           N_P1MAT_CI CI states in CIVEC
C
C     Written by H. J. Aa. Jensen, M. S. Vad and
C                M. N. Pedersen - July 2011
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbesr.h"
#include "dcbgen.h"
C
!     output
      DIMENSION P1MAT_CI(N_P1MAT_CI,N_P1MAT_CI,*)
!     input
      DIMENSION CIVEC(NDET,N_CIVEC,*)
      DIMENSION P1MAT_MO(2*NASHT,2*NASHT,*), WORK(LWORK)
      INTEGER*8 IDET(2,*)
      LOGICAL   PAR_GASCIP

C***********************************************************************
C
C     Memory initialization
C     =====================
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
C


       CALL MEMGET2('REAL','SVECPRP',KSVECPRP ,NZ_in_CI*NDET*N_CIVEC,
     &      WORK,KFREE,LFREE)
       CALL MEMGET2('REAL','BVECPRP',KBVECPRP ,NZ_in_CI*NDET,
     &      WORK,KFREE,LFREE)
       CALL MEMGET2('REAL','MUUUU',  KPMUUUU  ,0,
     &      WORK,KFREE,LFREE) ! no 2-electron part
       KBVECPRP_I = KBVECPRP + NDET
#if GASCIP_DEBUG > 10
       WRITE(LUPRI,*) 'NZ_in_CI, IHRM_P1 check, NDET'
     &               , NZ_in_CI, IHRM_P1,NDET

       WRITE(LUPRI,*) 'GASCIP_MAKE_P1MAT: Real PART of P1MAT_MO'
       CALL OUTPUT(P1MAT_MO,1,2*NASHT,1,2*NASHT,
     &             2*NASHT,2*NASHT,-1,LUPRI)
       IF (NZ_in_CI .EQ. 2) THEN
          WRITE(LUPRI,*) 'GASCIP_MAKE_P1MAT: Imag PART of P1MAT_MO'
          CALL OUTPUT(P1MAT_MO(1,1,2),1,2*NASHT,1,2*NASHT,
     &                2*NASHT,2*NASHT,-1,LUPRI)
       END IF


       WRITE(LUPRI,*)
     &   'GASCIP_MAKE_P1MAT: Real part of the N_CIVEC CI vectors'
       CALL OUTPUT(CIVEC,1,NDET,1,N_CIVEC,
     &             NDET,N_CIVEC,-1,LUPRI)
       IF (NZ_in_CI .EQ. 2) THEN
          WRITE(LUPRI,*)
     &      'GASCIP_MAKE_P1MAT: Imag part of the N_CIVEC CI vectors'
          CALL OUTPUT(CIVEC(1,1,2),1,NDET,1,N_CIVEC,
     &                NDET,N_CIVEC,-1,LUPRI)
       END IF
#endif

      PAR_GASCIP = PARCAL .AND. NDET .GT. 10000
      THRSCR = 1.0D-10
      IF (PAR_GASCIP) THEN
         IPRINT_slaves = IPRINT - 1
         CALL GASCIP_WAKE_NODES()
         CALL GASCIP_PAR_CISIG_SETUP(
     &            P1MAT_MO,WORK(KPMUUUU),
     &            .FALSE.,NDET,IDET,THRSCR,IPRINT_slaves)
      END IF

C      Make sigma vectors from property matrix and CI vectors
       DO I = 1,N_CIVEC
          CALL DCOPY(NDET,CIVEC(1,I,1),1,WORK(KBVECPRP),1)
          IF (NZ_in_CI .EQ. 2) THEN
             CALL DCOPY(NDET,CIVEC(1,I,2),1,WORK(KBVECPRP_I),1)
          END IF
          IF (PAR_GASCIP) THEN
             CALL GASCIP_PAR_CISIG(1,NDET,
     &            WORK(KSVECPRP+(I-1)*NDET*NZ_in_CI),
     &            WORK(KBVECPRP),IPRINT_slaves)
         ELSE IF (NDET .GT. 1000) THEN
C hj 16-jun-2004 new sigma code with screening
            CALL GASCIP_CISIG_SCREEN(1,NDET,
     &           WORK(KSVECPRP+(I-1)*NDET*NZ_in_CI),
     &           WORK(KBVECPRP),
     &           P1MAT_MO,WORK(KPMUUUU),
     &           .FALSE.,NDET,IDET,THRSCR,IPRINT)
C      SUBROUTINE GASCIP_CISIG_SCREEN(IS_START, IS_END,
C     &           SVEC,BVEC,HCORE,VMUUUU,
C     &           DO2EL,NDET,IDET,THRSCR,IPRINT)
         ELSE
            CALL GASCIP_CI(.FALSE.,DUMMY,
     &           .TRUE.,WORK(KSVECPRP+(I-1)*NDET*NZ_in_CI),
     &           WORK(KBVECPRP),P1MAT_MO,WORK(KPMUUUU),
     &           .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &           .FALSE.,NDET,IDET,IHRM_P1,IPRINT)
         END IF
       END DO
#if GASCIP_DEBUG > 20
       WRITE(LUPRI,'(/A)')
     &   ' GASCIP_MAKE_P1MAT: Sigma vectors = P * CI vectors'
       IF (NZ_in_CI.GT.1)
     & WRITE(LUPRI,'(A)')' Order: Real s_1, Imag s_1, Real s_2, etc.'
       CALL OUTPUT(WORK(KSVECPRP),1,NDET,1,NZ_in_CI*N_CIVEC,
     &             NDET,NZ_in_CI*N_CIVEC,-1,LUPRI)
#endif
C
C      Construct N_CIVEC x N_CIVEC property matrix
C                P_ij = C_i^+*sigma_j
      DO I=1,N_P1MAT_CI

C        Calculate Diagonal elements here!

             P1MAT_CI(I,I,1) =
     &               DDOT(NDET,WORK(KSVECPRP+(I-1)*NDET*NZ_in_CI)
     &                   ,1,CIVEC(1,I,1),1)
          IF (NZ_in_CI .EQ. 2) THEN
             P1MAT_CI(I,I,1) =  P1MAT_CI(I,I,1)
     &            +  DDOT(NDET,WORK(KSVECPRP+(I-1)*NDET*NZ_in_CI+NDET)
     &                   ,1,CIVEC(1,I,2),1)
             P1MAT_CI(I,I,2) =
     &            -  DDOT(NDET,WORK(KSVECPRP+(I-1)*NDET*NZ_in_CI)
     &                   ,1,CIVEC(1,I,2),1)
     &            +  DDOT(NDET,WORK(KSVECPRP+(I-1)*NDET*NZ_in_CI+NDET)
     &                   ,1,CIVEC(1,I,1),1)
          END IF


         JEND = I-1
C        Calculate off diagonal elements
         DO J=1,JEND
             P1MAT_CI(I,J,1) =
     &               DDOT(NDET,WORK(KSVECPRP+(J-1)*NDET*NZ_in_CI)
     &                   ,1,CIVEC(1,I,1),1)
          IF (NZ_in_CI .EQ. 2) THEN
             P1MAT_CI(I,J,1) =  P1MAT_CI(I,J,1)
     &            +  DDOT(NDET,WORK(KSVECPRP+(J-1)*NDET*NZ_in_CI+NDET)
     &                   ,1,CIVEC(1,I,2),1)
             P1MAT_CI(I,J,2) =
     &            -  DDOT(NDET,WORK(KSVECPRP+(J-1)*NDET*NZ_in_CI)
     &                   ,1,CIVEC(1,I,2),1)
     &            +  DDOT(NDET,WORK(KSVECPRP+(J-1)*NDET*NZ_in_CI+NDET)
     &                   ,1,CIVEC(1,I,1),1)
          END IF
               IF( IHRM_P1 .GT. 0 ) THEN
                  P1MAT_CI(J,I,1) =  P1MAT_CI(I, J,1)
                  IF (NZ_in_CI .EQ. 2) THEN
                  P1MAT_CI(J,I,2) = -P1MAT_CI(I, J,2)
                  END IF
               ELSE IF( IHRM_P1 .LT. 0 ) THEN
C                  P1MAT_CI(I,J,1) = -P1MAT_CI(I,J,1) ! apparently we have minus imag(p1)
                  P1MAT_CI(J,I,1) = -P1MAT_CI(I,J,1)
                  IF (NZ_in_CI .EQ. 2) THEN
                  P1MAT_CI(J,I,2) =  P1MAT_CI(I, J,2)
                  END IF
               END IF

         END DO
      END DO

      IF (IPRINT .GT. 1) THEN
         WRITE(LUPRI,'(/A)')
     &   ' (GASCIP_MAKE_P1MAT) CI property matrix - real part'
         CALL OUTPUT(P1MAT_CI,1,N_P1MAT_CI,1,N_P1MAT_CI,
     &               N_P1MAT_CI,N_P1MAT_CI,-1,LUPRI)
       IF (NZ_in_CI .GT. 1) THEN
         WRITE(LUPRI,'(/A)')
     &   ' (GASCIP_MAKE_P1MAT) CI property matrix - imag part:'
         CALL OUTPUT(P1MAT_CI(1,1,2),1,N_P1MAT_CI,1,N_P1MAT_CI,
     &               N_P1MAT_CI,N_P1MAT_CI,-1,LUPRI)
       END IF
      END IF

      IF (PAR_GASCIP) THEN
         CALL GASCIP_RELEASE_NODES()
      END IF

C     MEMREL for sigma vectors
      CALL MEMREL('GASCIP_MAKE_P1MAT',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_hmat */
      SUBROUTINE GASCIP_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
     &     HR,HI,HCORE,VMUUUU)
C***********************************************************************
C
C     Calculate Hamiltonian matrix element
C
C     Written by H. J. Aa. Jensen & J. Thyssen - Dec 30 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgascip.h"
#include "dcborb.h"
C
      LOGICAL   DO2EL,DOIMAG
      DIMENSION IOCCI(2*MAXASH_GASCIP), JOCCI(2*MAXASH_GASCIP)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
      DIMENSION HCORE(2*NASHT,2*NASHT,*)
C
C
C
      IF (IT .EQ. 0) THEN
C
C     ======================
C     DETERMINANTS ARE EQUAL
C     ======================
C
         HR = 0.0D0
         IF (DOIMAG) HI = 0.0D0
         DO 30 II = 1, NAELEC
            I = IOCCI(II)
            HR = HR + HCORE(I,I,1)
            IF (DOIMAG) HI = HI + HCORE(I,I,2)
           IF (DO2EL) THEN
            DO 20 JJ = II + 1, NAELEC
               J = IOCCI(JJ)
C---------------------
C              H = H + (I, I, J, J) - (I, J, J, I)
C---------------------
               HR = HR + VMUUUU(J,J,I,I,1) - VMUUUU(I,J,J,I,1)
               IF (DOIMAG)
     &         HI = HI + VMUUUU(J,J,I,I,2) - VMUUUU(I,J,J,I,2)
 20         CONTINUE
           END IF
 30      CONTINUE
C
      ELSE IF(IT.EQ.2) THEN
C
C     ===================================
C     DETERMINANTS DIFFER BY ONE FUNCTION K -> M
C     ===================================
C
         HR = HCORE(K,M,1)
         IF (DOIMAG) HI = HCORE(K,M,2)
C
        IF (DO2EL) THEN
         DO IELEC = 1, NAELEC
            I = IOCCI(IELEC)
            IF (I .NE. K) THEN
C---------------------
C              H = H + (I, I, K, M)
C              H = H - (I, M, K, I)
C---------------------
               HR = HR + VMUUUU(I,I,K,M,1) - VMUUUU(I,M,K,I,1)
C---------------------
               IF (DOIMAG)
     &              HI = HI + VMUUUU(I,I,K,M,2) - VMUUUU(I,M,K,I,2)
            END IF
         END DO
        END IF
C
         H_TEST = ABS(HR)
         IF (DOIMAG) H_TEST = H_TEST + ABS(HI)
         IF (H_TEST .GT. 1.0D-12) THEN
            IPHASE = 0
            DO IELEC = 1,NAELEC
               IF (IOCCI(IELEC) .EQ. K) IPHASE = IPHASE + IELEC
               IF (JOCCI(IELEC) .EQ. M) IPHASE = IPHASE + IELEC
            END DO
            ! IPHASE = (-1)**IPHASE
            IF (IAND(IPHASE,1) .EQ. 1) THEN ! IPHASE odd or even?
               HR = -HR
               IF (DOIMAG) HI = -HI
            END IF
         END IF
C
      ELSE IF(IT.EQ.4 .AND. DO2EL) THEN
C
C     ====================================
C     DETERMINANTS DIFFER BY TWO FUNCTIONS  K,L -> M,N
C     ====================================
C
C
C---------------------
C        H = H + (K, M, L, N) - (K, N, L, M)
C---------------------
         HR = (VMUUUU(K,M,L,N,1) - VMUUUU(K,N,L,M,1))
         H_TEST = ABS(HR)
         IF (DOIMAG) THEN
            HI = (VMUUUU(K,M,L,N,2) - VMUUUU(K,N,L,M,2))
            H_TEST = H_TEST + ABS(HI)
         END IF
         IF (H_TEST .GT. 1.0D-12) THEN
            IPHASE = 0
            DO IELEC = 1,NAELEC
               IF (IOCCI(IELEC) .EQ. K) IPHASE = IPHASE + IELEC
               IF (IOCCI(IELEC) .EQ. L) IPHASE = IPHASE + IELEC
               IF (JOCCI(IELEC) .EQ. M) IPHASE = IPHASE + IELEC
               IF (JOCCI(IELEC) .EQ. N) IPHASE = IPHASE + IELEC
            END DO
            ! IPHASE = (-1)**IPHASE
            IF (IAND(IPHASE,1) .EQ. 1) THEN ! IPHASE odd or even?
               HR = -HR
               IF (DOIMAG) HI = -HI
            END IF
         END IF

      ELSE
C
C     =============================================
C     DETERMINANTS DIFFER BY MORE THEN TWO FUNCTIONS
C     ==============================================
C
         CALL QUIT('GASCIP_HMAT called with IT .gt. 4 !')
      END IF
 9999 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_svc */
      SUBROUTINE GASCIP_SVC(SAVDIA,DIAG,SVEC,NDET,IDET,BVEC,
     &                      HCORE,VMUUUU,IPRINT)
C***********************************************************************
C
C     Calculate CI sigma vectors: SVEC(I) = \sum_J H_{IJ} B_{J}
C
C     Input:
C        HCORE  - core Fock matrix
C        VMUUUU - two electron integrals.
C        BVEC   - trial vector
C        IDET   - determinant strings
C        SAVDIA - save diagonal of Hamiltonian (H_{II})
C
C     Output:
C        SVEC   - sigma vector
C        DIAG   - diagonal of Hamiltonian (if SAVDIA)
C
C     Written by J. Thyssen - Jan 3 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      LOGICAL   SAVDIA
      DIMENSION DIAG(*), SVEC(*), BVEC(*), HCORE(*), VMUUUU(*)
      INTEGER*8 IDET(2,*)
C
      CALL QENTER('GASCIP_SVC')
C
      IF (NDET .GT. 1000) THEN
C hj 16-jun-2004 new sigma code with screening
         THRSCR = 1.0D-10
         IHRM = 1
         CALL GASCIP_CISIG_SCREEN(1,NDET,
     &      SVEC,BVEC,HCORE,VMUUUU,.TRUE.,
     &      NDET,IDET,THRSCR,IPRINT)
         IF (SAVDIA) CALL GASCIP_CI(SAVDIA,DIAG,
     &     .FALSE.,SVEC,BVEC,HCORE,VMUUUU,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,NDET,IDET,IHRM,IPRINT)
      ELSE
         IHRM = 1
         CALL GASCIP_CI(SAVDIA,DIAG,
     &     .TRUE.,SVEC,BVEC,HCORE,VMUUUU,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,NDET,IDET,IHRM,IPRINT)
      END IF
C
C     CALL GASCIP_CI(
C    &     DODIA,DIAG,
C    &     DOSIG,SVEC,BVEC,HCORE,VMUUUU,
C    &     DODEN,CL,CR,DV,PV,
C    &     DO2EL,NDET,IDET,IHRM,IPRINT)
C
      CALL QEXIT('GASCIP_SVC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_hdiag */
      SUBROUTINE GASCIP_HDIAG(DIAG,NDET,IDET,HCORE,VMUUUU,IPRINT)
C***********************************************************************
C
C     Calculate diagonal of CI Hamiltonian.
C
C     Input:
C        HCORE  - core Fock matrix
C        VMUUUU - two electron integrals.
C        IDET   - determinant strings
C
C     Output:
C        DIAG   - diagonal of Hamiltonian (if SAVDIA)
C
C     Written by J. Thyssen - Jan 3 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION DIAG(*), IDET(*), HCORE(*), VMUUUU(*)
C
      CALL QENTER('GASCIP_HDIAG')
C
      IHRM = 1
      CALL GASCIP_CI(.TRUE.,DIAG,
     &     .FALSE.,DUMMY,DUMMY,HCORE,VMUUUU,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,NDET,IDET,IHRM,IPRINT)
C
C     CALL GASCIP_CI(
C    &     DODIA,DIAG,
C    &     DOSIG,SVEC,BVEC,HCORE,VMUUUU,
C    &     DODEN,CL,CR,DV,PV,
C    &     DO2EL,NDET,IDET,IHRM,IPRINT)
C
      CALL QEXIT('GASCIP_HDIAG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_rcistd */
      SUBROUTINE GASCIP_RCISTD(N_ROOTS_in,THR_CONV,MAX_CI_IT,
     &                         ECI_OUT,CVECS,NDET,IDET,CMO,
     &                         THR_CVEC,I_STATE,ICONV_TOT,
     &                         THR_PCI,IPRINT,WORK,LWORK)
C***********************************************************************
C
C     Do up to MAX_CI_IT CI iterations to converge the N_ROOTS lowest states
C     to a gradient threshold of THR_CONV .
C
C     Input :
C        N_ROOTS_in - abs value = number of roots to converge;
C                  if negative do not increase N_CSIM because of (near) degeneracy
C        THR_CONV  - CI gradient convergence threshold
C        MAX_CI_IT - max CI iterations
C        NDET    - number of determinants
C        IDET    - determinant strings
C        CMO     - MO coefficients
C        I_STATE - if > 0: which state converge to (not used p.t.)
C        IPRINT  - print level
C
C     Output:
C        ECI_OUT - CI energies
C        CVECS   - CI vectors
C        ICONV_TOT =1 all converged, <0 code for why not converged
C
C
C     Written by J. Thyssen - Jan 3 2001
C     Last revision :  Hans Jorgen A. Jensen - Oct 2010/parallel Aug 2014
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "thrzer.h"
C
C Used from include files:
C  dcbgen.h : PARCAL, ...
C  dcborb.h : NASHT,NCMOTQ,...
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
C
      INTEGER*8 IDET(2,*)
      DIMENSION ECI_OUT(*), CVECS(NDET*NZ_in_CI,*), CMO(*), WORK(*)
C
C     Local variales:
C
      PARAMETER ( MAX_ROOTS = 41, N_CSIM_EXTRA_START = 7)
      real*8    DNRES(MAX_ROOTS), ECIOLD(MAX_ROOTS), ECIA(MAX_ROOTS)
      INTEGER   J_HDIAG(MAX_ROOTS + N_CSIM_EXTRA_START)
      INTEGER   ICONV(MAX_ROOTS), LUGASCIP, IOSVAL
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
      real*8    ECI(MAX_ROOTS)
      LOGICAL   PAR_GASCIP
      LOGICAL   is_open,file_exists
C
#include "memint.h"
C
      CALL QENTER('GASCIP_RCISTD')

      N_ROOTS = abs(N_ROOTS_in)

      CALL GETTIM(CPU1,WALL1)

C
C     Set general parameters (which later maybe should be set by input?)
C
!d      write(lupri,*) 'THR_PCI, THR_CONV', THR_PCI,THR_CONV
      NDETQ  = NZ_in_CI*NDET
      MAX_HRED = MAX(100,N_ROOTS*MAX_CI_IT)
      MAX_HRED = MIN(2000, MAX_HRED)
C
      IF (N_ROOTS .LT. 1 .OR. N_ROOTS .GE. MAX_ROOTS) THEN
C        MAX_ROOTS must be greater than  N_ROOTS to be sure
C        to include all members of quasidegenerate sets.
         WRITE(LUPRI,*) 'N_ROOTS .lt. 1 .or. N_ROOTS .ge. MAX_ROOTS;'//
     &      ' NROOTS, MAX_ROOTS =',N_ROOTS,MAX_ROOTS
         CALL QUIT('N_ROOTS .lt. 1 .or N_ROOTS .gt. MAX_ROOTS')
      END IF
      IF (N_ROOTS .GT. NDET) THEN
         WRITE(LUPRI,*) 'N_ROOTS > number of determinants',
     &      N_ROOTS,NDET
         CALL QUIT('N_ROOTS .gt. number of determinants')
      END IF

C
C     ************************************
C     *** Calculate necesary integrals ***
C     ************************************
C
C     Calculate core Hamiltonian:
C     ---------------------------
C
C     The core Hamiltonian is the active-active part of FC,
C     a.k.a. FCAC.
C
      CALL MEMGET2('REAL','FCACM',KFCACM,(2*NASHT)*(2*NASHT)*NZ_in_CI,
     &     WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','FCAC',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','FC',  KFC,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Calculate FCmo (TOFILE false: do not write to LUKRM3) :
C
      CALL rGETFC(.FALSE.,CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
C     Get active-active part of FCmo:
C
      CALL RGETAC(WORK(KFC),WORK(KFCAC),IPRINT)
      IF (IPRINT.GE.6) THEN
         WRITE(LUPRI,'(A)') '(GASCIP_RCISTD) Active Fock matrix:'
         CALL OUTPUT(WORK(KFCAC),1,NASHT,1,NASHT,NASHT,NASHT,-1,LUPRI)
      END IF
C
C     Transform quaternion FCACmo to molfdir type FCAC.
C
      CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPRINT)
C
      CALL MEMREL('RCISTD after rGETFC etc.',WORK,1,KFCAC,KFREE,LFREE)
C
C     Read two-electron integrals:
C     ----------------------------
C
      CALL MEMGET2('REAL','MUUUU',KMUUUU,(2*NASHT)**4 * 2,
     &     WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','H2AC',KH2AC,NASHT*NASHT*NNASHX*NZ*3,
     &     WORK,KFREE,LFREE)
C
C     read integrals from 4IND*: out --> work(kh2ac)
C     ==============================================

      call memget2('INTE','IBEIG',kibeig,norbt,work,kfree,lfree)

      call izero(work(kibeig),norbt)
      if(spinfr.or.levyle)then
        call ireakrmc(lukrmc,'IBEIG   ',work(kibeig),norbt)
      else if(linear)then
        call ireakrmc(lukrmc,'MJVEC   ',work(kibeig),norbt)
      end if

      call rgeth2(dummy,work(kh2ac),dummy,work(kibeig),.false.,
     &            .true.,.true.,work(kfree),lfree)
C
C     Make sure imag(VMUUUU) = 0 for NZ=1
C     (DNZ32M only DZERO'es real part for NZ = 1)
C
      IF (NZ_in_CI .gt. 1)
     &  CALL DZERO(WORK(KMUUUU+(2*NASHT)**4),(2*NASHT)**4)
C
C     Transform integrals to Molfdir format:
C
      if(naelec > 1)then
        CALL DNZ32M(WORK(KH2AC),WORK(KMUUUU),IPRINT)
      else
        CALL DZERO(WORK(KMUUUU),(2*NASHT)**4)
      end if

      CALL MEMREL('GASCIP_RCISTD.H2AC',WORK,KWORK,KH2AC,KFREE,LFREE)

!     Allocate memory for reduced space and for b- and sigma-vectors

      LHRED = MAX_HRED*MAX_HRED*NZ_in_CI
      CALL MEMGET2('REAL','HRED',KHRED,LHRED,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','ERED',KEIGVL,MAX_HRED,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','VRED',KEIGVC,LHRED,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KHRED),LHRED)
      KHREDI  = KHRED + MAX_HRED*MAX_HRED
      KEIGVCI = KEIGVC + MAX_HRED*MAX_HRED
C
C     Assume we can have all sigma and trial vectors in memory.
C     This could get a little expensive for large CI's.
C
      CALL MEMGET2('REAL','BCVEC',KBCVEC,MAX_HRED * NDETQ,
     &     WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','SCVEC',KSCVEC,MAX_HRED * NDETQ,
     &     WORK,KFREE,LFREE)
C
C
C     ********************************
C     *** Find start trial vectors ***
C     ********************************
C
C
      CALL MEMGET2('REAL','HDIAG',KHDIAG,NDET,WORK,KFREE,LFREE)
C
C     Calculate diagonal of CI matrix:
C     --------------------------------
C
c     CALL GASCIP_HDIAG(WORK(KHDIAG),NDET,IDET,WORK(KFCACM),WORK(KMUUUU))
      IHRM = 1
      CALL GASCIP_CI(.TRUE.,WORK(KHDIAG),
     &     .FALSE.,DUMMY,DUMMY,WORK(KFCACM),WORK(KMUUUU),
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,NDET,IDET,IHRM,IPRINT)
C     CALL GASCIP_CI(DODIA,DIAG,
C          DOSIG,SVEC,BVEC,HCORE,VMUUUU,
C          DODEN,CL,CR,DV,PV,
C          DO2EL,NDET,IDET,IHRM,IPRINT)
C
      IF (IPRINT .GE. 10) THEN
         WRITE(LUPRI,'(/A)')
     &      ' (GASCIP_RCISTD) Diagonal of CI Hamiltonian'
         CALL OUTPUT(WORK(KHDIAG),1,1,1,NDET,1,NDET,-1,LUPRI)
      END IF

      ! Find out if this is a restart
      inquire(file='GASCIP.RST',opened=is_open,exist=file_exists,
     &              number=LUGASCIP)
      IF (file_exists) THEN
         IF (.NOT. is_open) THEN
            LUGASCIP = LU_NOT_USED(51)
            open(LUGASCIP,file='GASCIP.RST',
     &         status='OLD',form='UNFORMATTED')
         END IF
         REWIND (LUGASCIP)
         READ (LUGASCIP,IOSTAT=IOSVAL) NCMOTQ_r, NDET_r, NDETQ_r, N_CVEC
         IF (IOSVAL .NE. 0) THEN
            N_CVEC = 0
            REWIND (LUGASCIP)
         ELSE
         IF (N_CVEC .GT. 0 .AND. NCMOTQ .EQ. NCMOTQ_r .AND.
     &      NDET .EQ. NDET_r .AND. NDETQ_r .EQ. NDETQ) THEN
         ! GASCIP.RST seems to be OK for restart
            CALL READT(LUGASCIP,NDETQ*N_CVEC,WORK(KBCVEC))
            READ (LUGASCIP) DNRES(1:N_CVEC)
            N_CSIM = N_CVEC
            WRITE(LUPRI,'(/I5,A/,(1P,10D10.2))') N_CVEC,
     &      ' restart CI vectors read from GASCIP.RST; residual norms:',
     &      DNRES(1:N_CVEC)

            IF (N_CSIM .LT. N_ROOTS) GO TO 20 ! need to converge more roots

            ICONV_TOT = 1
            IF (I_STATE .GT. 0) THEN ! only converge this root
              IF ( DNRES(I_STATE) .GT. THR_CONV ) ICONV_TOT = 0 ! the requested state is not converged
            ELSE
              DO I_ROOT = 1, N_ROOTS
                IF ( DNRES(I_ROOT) .GT. THR_CONV ) ICONV_TOT = 0  ! not all converged
              END DO ! I_ROOT
            END IF

            IF (ICONV_TOT .EQ. 1) THEN ! requested roots are already converged
               CALL DCOPY(NDETQ*N_ROOTS,WORK(KBCVEC),1,CVECS,1)
               READ (LUGASCIP) ECIA(1:N_ROOTS)
               ITMIC = -1
               GO TO 210
            ELSE
               GO TO 50
            END IF
         ELSE
            N_CVEC = 0
         END IF
         END IF ! IOSVAL .eq. 0
      ELSE
         LUGASCIP = LU_NOT_USED(51)
         open(LUGASCIP,file='GASCIP.RST',
     &      status='NEW',form='UNFORMATTED')
         N_CVEC = 0
      END IF

      ! if here: no or not enough restart CI vectors
C
C     If N_ROOTS .gt. 0, then
C     find N_ROOTS + N_CSIM_EXTRA_START lowest elements in HDIAG:
C
C     (We start with extra diagonal elements to be fairly sure
C      that we cover the N_ROOTS lowest symmetries in the start.
C      Say that diagonal element no. N_ROOTS + 1 is the first of
C      symmetry 2, but after convergence symmetry 2 is eigenvalue
C      no. N_ROOTS - 2. This eigenvalue and eigenvector will be
C      missed if we don't start with some extra diagonal elements.
C      /hjaaj Oct 2010 )
C
   20 CONTINUE
      IF (N_ROOTS_in .GT. 0) THEN
         N_CSIM = MIN(NDET, N_ROOTS+N_CSIM_EXTRA_START, MAX_ROOTS)
         J_CSIM = MIN(NDET, N_CSIM+N_CSIM_EXTRA_START)
         ! always OK because DIMENSION J_HDIAG(MAX_ROOTS + N_CSIM_EXTRA_START)
      ELSE
         N_CSIM = N_ROOTS
         J_CSIM = N_ROOTS
      END IF
      CALL FNDMIN(J_CSIM,J_HDIAG,WORK(KHDIAG),NDET,WORK(KFREE),LFREE)
C
C     Make sure that we don't break symmetry by only having all
C     components of a (nearly) degenerate set - criterium for
C     degeneracy: energy difference less than 0.01 bohr)

      !> hardwire for debugging:
!     n_csim = 4
C
      if (N_CSIM .GT. N_ROOTS .AND. N_CSIM .LT. NDET) then
        HDIAG_N_CSIMp1 = WORK(KHDIAG-1 + J_HDIAG(N_CSIM+1))
   40   HDIAG_N_CSIM   = WORK(KHDIAG-1 + J_HDIAG(N_CSIM))
        call flshfo(lupri)
        IF (ABS(HDIAG_N_CSIMp1 - HDIAG_N_CSIM) .LT. 0.002D0) THEN
           N_CSIM = N_CSIM - 1
           HDIAG_N_CSIMp1 = HDIAG_N_CSIM
           GO TO 40
        END IF
      end if
      IF (N_CSIM .LT. N_ROOTS) THEN
          WRITE(LUPRI,*)
     &      'ERROR: N_CSIM .lt. N_ROOTS; increase MAX_ROOTS',
     &      N_CSIM,N_ROOTS
         CALL QUIT('N_CSIM .lt. N_ROOTS; increase MAX_ROOTS')
!        it will be necessary to increase MAX_ROOTS to be sure
!        to include all members of a quasi-degenerate set.
      END IF
C
C     Set first trial vectors equal to: BCVEC(I,J) = DELTA(I, J_HDIAG(J))
C
      WRITE(LUPRI,'(/1X,I0,A)') N_CSIM-N_CVEC,
     &   ' start vectors based on lowest diagonal elements of H_CI'
      JBCVEC = KBCVEC - 1
      DO I = 1, N_CSIM
         II = J_HDIAG(I)
         ECIA(I) = WORK(KHDIAG-1+II)
         ECI(I)  = ECIA(I) + ECORE
         IF (I .LE. N_CVEC) THEN
            WRITE (LUPRI,'(I10,F20.6)') II, ECI(I)
         ELSE
            WRITE (LUPRI,'(I10,F20.6,A)') II, ECI(I), ' (used)'
            CALL DZERO(WORK(JBCVEC+1),NDETQ)
            WORK(JBCVEC + II) = 1.0D0
         END IF
         JBCVEC = JBCVEC + NDETQ
      END DO

!     N_CSIM start vectors available, start micro iterations

   50 CONTINUE

      ITMIC = 0
      ICONV_TOT = 0
      N_CRED = 0
C     GASCIP screening, a calculation in C2 symmetry showed that
C     THRSCR = THR_CONV**2 was not good enough for the complex group /hjaaj July 04
      THRSCR = MIN(1.0D-8, THR_CONV**2 * 0.01D0)

      PAR_GASCIP = PARCAL .AND. NDET .GT. 10000
      IF (PAR_GASCIP) THEN
         IPRINT_slaves = IPRINT - 1
         CALL GASCIP_WAKE_NODES()
         CALL GASCIP_PAR_CISIG_SETUP(
     &            WORK(KFCACM),WORK(KMUUUU),
     &            .TRUE.,NDET,IDET,THRSCR,IPRINT_slaves)
      END IF

      !> for debug
!     iprint_save = iprint
!     iprint      = 5
C
C     **********************************
C     *** Start micro iteration loop ***
C     **********************************
C
      CALL GETTIM(CPU2,WALL2)
C
 100  CONTINUE  ! DO WHILE (not converged and no limits reached)
         N_CRED_OLD = N_CRED
         N_CRED = N_CRED + N_CSIM
         ITMIC = ITMIC + 1
         CALL GETTIM(CPU3,WALL3)
         CPU_st  = CPU3
         WALL_st = WALL3
C
         WRITE(LUPRI,'(//A,I0,A,I0)')
     &        ' (GASCIP_RCISTD) CI microiteration no. ',ITMIC,
     &        ' - number of trial vectors: ',N_CSIM
         IF (N_CSIM .LE. 0) THEN
            WRITE (LUPRI,'(/A,I0)')
     &         ' FATAL ERROR: no new trial vectors ',N_CSIM
            CALL QUIT(' FATAL ERROR: no new trial vectors')
         END IF
         CALL FLSHFO(LUPRI)
C
C        Calculate sigma vector(s)
C        -------------------------
C
         IF (PAR_GASCIP) THEN
            CALL GASCIP_PAR_CISIG(N_CSIM,NDET,
     &            WORK(KSCVEC+N_CRED_OLD*NDETQ),
     &            WORK(KBCVEC+N_CRED_OLD*NDETQ),IPRINT_slaves)
         END IF

         IHRM = 1
         DO I = 1, N_CSIM

            IF (.NOT. PAR_GASCIP) THEN ! if PAR_GASCIP then SCVEC is calculated above
            I_TOT = N_CRED_OLD + I
            IF (NDET .GT. 1000) THEN
C hj 16-jun-2004 new sigma code with screening
               CALL GASCIP_CISIG_SCREEN(1,NDET,
     &            WORK(KSCVEC+(I_TOT-1)*NDETQ),
     &            WORK(KBCVEC+(I_TOT-1)*NDETQ),
     &            WORK(KFCACM),WORK(KMUUUU),
     &            .TRUE.,NDET,IDET,THRSCR,IPRINT)
            ELSE
               CALL GASCIP_CI(.FALSE.,DUMMY,
     &            .TRUE.,WORK(KSCVEC+(I_TOT-1)*NDETQ),
     &            WORK(KBCVEC+(I_TOT-1)*NDETQ),
     &            WORK(KFCACM),WORK(KMUUUU),
     &            .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &            .TRUE.,NDET,IDET,IHRM,IPRINT)
C              CALL GASCIP_CI(DODIA,DIAG,
C                 DOSIG,SVEC,BVEC,
C                 HCORE,VMUUUU,
C                 DODEN,CL,CR,DV,PV,
C                 DO2EL,NDET,IDET,IHRM,IPRINT)
C
            END IF
            END IF ! (.NOT. PAR_GASCIP)
c
            IF (IPRINT .GE. 10) THEN
C
C              Trial vector
C
               WRITE(LUPRI,'(/A,I3,A,I3)')
     &              ' (GASCIP_RCISTD) Trial vector no. ',I,' of ',N_CSIM
               CALL GASCIP_ANACI(WORK(KBCVEC+(I_TOT-1)*NDETQ),
     &                     NDET,IDET,THR_PCI,WORK(KFREE),LFREE)
               IF (IPRINT .GE. 20)
     &            CALL OUTPUT(WORK(KBCVEC+(I_TOT-1)*NDETQ),
     &              1,NDET,1,NZ_in_CI,NDET,NZ_in_CI,1,LUPRI)
C
C              Sigma vector
C
               WRITE(LUPRI,'(/A,I3,A,I3)')
     &              ' (GASCIP_RCISTD) Sigma vector no. ',I,' of ',N_CSIM
               CALL GASCIP_ANACI(WORK(KSCVEC+(I_TOT-1)*NDETQ),
     &                     NDET,IDET,THR_PCI,WORK(KFREE),LFREE)
               IF (IPRINT .GE. 20)
     &            CALL OUTPUT(WORK(KSCVEC+(I_TOT-1)*NDETQ),
     &              1,NDET,1,NZ_in_CI,NDET,NZ_in_CI,1,LUPRI)
            END IF
         END DO
      IF (IPRINT .GE. 3) THEN
         CALL GETTIM(CPU_end,WALL_end)
         CPUTID  = SECTID(CPU_end-CPU_st)
         WALLTID = SECTID(WALL_end-WALL_st)
         WRITE(LUPRI,'(/5A)')
     &        ' (GASCIP_RCISTD) CPU (Wall) time for sigma vectors :',
     &        CPUTID,' (',WALLTID,')'
         CPU_st = CPU_end
         WALL_st = WALL_end
      END IF
C
C        Extend reduced Hamiltonian:
C        ---------------------------
C
C        HRED(I,J) = BCVEC(I)^+ * SCVEC(J)
C
C        (Old + New) sigma vectors * New trial vectors
C

         IF (IPRINT .GE. 20) THEN
             write (lupri,*) 'All N_CRED trial vectors',N_CRED
             CALL OUTPUT(WORK(KBCVEC),
     &                   1,NDETQ,1,N_CRED,NDETQ,N_CRED,-1,LUPRI)
             write (lupri,*) 'All N_CRED sigma vectors',N_CRED
             CALL OUTPUT(WORK(KSCVEC),
     &                   1,NDETQ,1,N_CRED,NDETQ,N_CRED,-1,LUPRI)

         END IF
         DO J = 1, N_CRED   ! All sigma vectors
            DO I = N_CRED_OLD+1, N_CRED   ! New trial vectors
               WORK(KHRED+(J-1)*MAX_HRED+(I-1)) =
     &              DDOT(NDETQ,WORK(KBCVEC+(I-1)*NDETQ),1,
     &                         WORK(KSCVEC+(J-1)*NDETQ),1)
               IF (NZ_in_CI .GT. 1) THEN
                  WORK(KHREDI+(J-1)*MAX_HRED+(I-1)) =
     &              DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ),1,
     &                        WORK(KSCVEC+(J-1)*NDETQ+NDET),1)
     &            - DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ+NDET),1,
     &                        WORK(KSCVEC+(J-1)*NDETQ),1)
               END IF
            END DO
         END DO
C
C        New sigma vectors * Old trial vectors
C
         DO J = N_CRED_OLD+1, N_CRED  ! New sigma vectors
            DO I = 1, N_CRED_OLD   ! Old trial vectors
               WORK(KHRED+(J-1)*MAX_HRED+(I-1)) =
     &              DDOT(NDETQ,WORK(KBCVEC+(I-1)*NDETQ),1,
     &                         WORK(KSCVEC+(J-1)*NDETQ),1)
               IF (NZ_in_CI .GT. 1) THEN
                  WORK(KHREDI+(J-1)*MAX_HRED+(I-1)) =
     &              DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ),1,
     &                        WORK(KSCVEC+(J-1)*NDETQ+NDET),1)
     &            - DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ+NDET),1,
     &                        WORK(KSCVEC+(J-1)*NDETQ),1)
               END IF
            END DO
         END DO
C
C        Check for inconsistencies:
C
chj      1.0D-10 can give problems with 1.0D-10 screening in GASCIP
C        THRDAS = 1.0D-10
         THRDAS = MAX(1.0D-09,THRSCR*100.0D0)
         DASR = 0.0D0
         DASI = 0.0D0
         DO I = 1, N_CRED
            DO J = 1, I - 1
               DASR = DASR +
     &            ABS( WORK(KHRED + (J-1) + (I-1)*MAX_HRED) -
     &                 WORK(KHRED + (I-1) + (J-1)*MAX_HRED) )
               IF (NZ_in_CI .GT. 1) DASI = DASI +
     &            ABS( WORK(KHREDI + (J-1) + (I-1)*MAX_HRED) +
     &                 WORK(KHREDI + (I-1) + (J-1)*MAX_HRED) )
            END DO
         END DO
         DAS = DASR + DASI

         IF (IPRINT .GE. 5 .OR. DAS .GT. THRDAS) THEN

            WRITE(LUPRI,'(/A)')
     &           ' (GASCIP_RCISTD) Reduced Hamiltonian, real part'
            CALL OUTPUT(WORK(KHRED),1,N_CRED,1,N_CRED,
     &                  MAX_HRED,MAX_HRED,-1,LUPRI)
           IF (NZ_in_CI .GT. 1) THEN
               WRITE(LUPRI,'(/A)')
     &           ' (GASCIP_RCISTD) Reduced Hamiltonian, imaginary part'
               CALL OUTPUT(WORK(KHREDI),1,N_CRED,1,N_CRED,
     &                  MAX_HRED,MAX_HRED,-1,LUPRI)
           END IF
         END IF
         IF (IPRINT .GE. 5 .OR. DAS .GT. 0.01D0*THRDAS) THEN
            WRITE(LUPRI,'(/A,1P,2D10.3)') ' (GASCIP_RCISTD) '//
     &         'Anti-hermiticity of reduced Hamiltonian (real,imag)',
     &         DASR, DASI
            IF (DAS .GT. THRDAS) THEN
               WRITE(LUPRI,'(/A,1P,D10.3)')
     &          ' FATAL ERROR, this is greater than threshold',THRDAS
               DO I = 1, N_CRED
                  DO J = 1, I
                     DASR = WORK(KHRED + (J-1) + (I-1)*MAX_HRED)
     &                    - WORK(KHRED + (I-1) + (J-1)*MAX_HRED)
                     WORK(KHRED + (J-1) + (I-1)*MAX_HRED) = DASR*0.5D0
                     WORK(KHRED + (I-1) + (J-1)*MAX_HRED) = DASR*0.5D0
                  IF (NZ_in_CI .GT. 1) THEN
                     DASI = WORK(KHREDI + (J-1) + (I-1)*MAX_HRED)
     &                    + WORK(KHREDI + (I-1) + (J-1)*MAX_HRED)
                     WORK(KHREDI + (J-1) + (I-1)*MAX_HRED) = DASI*0.5D0
                     WORK(KHREDI + (I-1) + (J-1)*MAX_HRED) = DASI*0.5D0
                  END IF
                  END DO
               END DO
               WRITE(LUPRI,'(/A)') ' (GASCIP_RCISTD) '//
     &            'Antihermicity of reduced Hamiltonian, real part'
               CALL OUTPUT(WORK(KHRED),1,N_CRED,1,N_CRED,
     &                     MAX_HRED,MAX_HRED,-1,LUPRI)
               IF (NZ_in_CI .GT. 1) THEN
                  WRITE(LUPRI,'(/A)') ' (GASCIP_RCISTD) '//
     &            'Antihermicity of reduced Hamiltonian, imaginary part'
                  CALL OUTPUT(WORK(KHREDI),1,N_CRED,1,N_CRED,
     &                        MAX_HRED,MAX_HRED,-1,LUPRI)
               END IF

               CALL QUIT('Anti-hermiticity of red. Ham. > threshold')
            END IF
         END IF
C
C        Diagonalize reduced Hamiltonian
C        -------------------------------
C
         IF (NZ_in_CI.EQ.1) THEN
            IF (IPRINT.GT.5) write(lupri,*) 'calling RSJACO'
            IJOB   = 1
            IORDER = 1
            IPACK  = 0
            CALL RSJACO(MAX_HRED,N_CRED,N_CRED,WORK(KHRED),WORK(KEIGVL),
     &        IJOB,IORDER,IPACK,WORK(KEIGVC))
         ELSE
            IF (IPRINT.GT.5) write(lupri,*) 'calling QDIAG'
            MATZ = 1
            CALL MEMGET2('REAL','HREDCP',KHREDCP,LHRED,WORK,KFREE,LFREE)
            CALL DCOPY(LHRED,WORK(KHRED),1,WORK(KHREDCP),1)
            CALL QDIAG(2,N_CRED,WORK(KHREDCP),MAX_HRED,MAX_HRED,
     &        WORK(KEIGVL),MATZ,WORK(KEIGVC),MAX_HRED,MAX_HRED,
     &        WORK(KFREE),LFREE,IERR)
            IF (IERR .NE. 0) THEN
               WRITE(LUPRI,'(/2A,I4)')
     &           ' *** ERROR in GASCIP_RCISTD ***: ',
     &           'QDIAG failed with error code ',IERR
               CALL QUIT('*** ERROR in GASCIP_RCISTD *** QDIAG failed')
            END IF
            CALL MEMREL('GASCIP_RCISTD.QDIAG',
     &         WORK,KWORK,KHREDCP,KFREE,LFREE)
         END IF
C
         IF (IPRINT .GE. 5) THEN
            WRITE(LUPRI,'(/A)')
     &         ' (GASCIP_RCISTD) Eigenvalues of reduced Hamiltonian'
            CALL OUTPUT(WORK(KEIGVL),1,1,1,N_CRED,1,MAX_HRED,-1,LUPRI)
         END IF
         IF (IPRINT .GE. 6) THEN
            WRITE(LUPRI,'(/A)') ' (GASCIP_RCISTD)'//
     &      ' Eigenvectors of reduced Hamiltonian, real part'
            CALL OUTPUT(WORK(KEIGVC),1,N_CRED,1,N_CRED,
     &                  MAX_HRED,MAX_HRED,-1,LUPRI)
           IF (NZ_in_CI .GT. 1) THEN
            WRITE(LUPRI,'(/A)') ' (GASCIP_RCISTD)'//
     &      ' Eigenvectors of reduced Hamiltonian, imaginary part'
            CALL OUTPUT(WORK(KEIGVCI),1,N_CRED,1,N_CRED,
     &                  MAX_HRED,MAX_HRED,-1,LUPRI)
           END IF
         END IF
         ECIOLD(1:N_CSIM) = ECI(1:N_CSIM)
         ECIA(1:N_CSIM)   = WORK(KEIGVL:KEIGVL-1 + N_CSIM)
         ECI(1:N_CSIM)    = ECORE + ECIA(1:N_CSIM)
         IF (IPRINT .GE. 10) THEN
            WRITE(LUPRI,'(/A,F20.8)')
     &        ' (GASCIP_RCISTD) Core energy        =',ECORE
            WRITE(LUPRI,'(A,4F20.8,/,(37X,4F20.8))')
     &        '                 Active energy      =',ECIA(1:N_CSIM)
         END IF
         WRITE(LUPRI,'(A,4F20.8/,(37X,4F20.8))')
     &        ' (GASCIP_RCISTD) Total energy       =',ECI(1:N_CSIM)
         WRITE(LUPRI,'(A,4F20.8/,(37X,4F20.8))')
     &        '                 Lowering of energy =',
     &        (ECIOLD(I)-ECI(I),I=1,N_CSIM)

         IF (IPRINT .GE. 3) THEN
            CALL GETTIM(CPU_end,WALL_end)
            CPUTID  = SECTID(CPU_end-CPU_st)
            WALLTID = SECTID(WALL_end-WALL_st)
            WRITE(LUPRI,'(/5A)')
     &     ' (GASCIP_RCISTD) CPU (Wall) time for reduced Hamiltonian :',
     &      CPUTID,' (',WALLTID,')'
            CPU_st = CPU_end
            WALL_st = WALL_end
         END IF
C
C        Construct residuals
C        -------------------
C
C        r = (H-E)x  (store in CVECS)
C          = Hx - Ex
C          = SUM_J C(J) SCVEC(J) - ECIA * SUM_J C(J) BCVEC(J)
C

         ! 1) construct current best solution vector SUM_J C(J,I_ROOT) BCVEC(J)
         !    (we cannot use QGEMM because of packing of BCVEC and CVECS)
         ! 1a) the real part of eigen vectors C(:,:)
         CALL DGEMM('N','N',NDETQ,N_ROOTS,N_CRED,1.0D0,
     &        WORK(KBCVEC),NDETQ,          ! real + imag. part
     &        WORK(KEIGVC),MAX_HRED,       ! real part
     &        0.0D0,CVECS,NDETQ)           ! store in real + imag. part
         ! 1b) the imaginary part of eigen vectors C(:,:)
         IF (NZ_in_CI .gt. 1) THEN
           CALL DGEMM('N','N',NDET,N_ROOTS,N_CRED,1.0D0,
     &        WORK(KBCVEC),NDETQ,          ! real part
     &        WORK(KEIGVCI),MAX_HRED,      ! imaginary part
     &        1.0D0,CVECS(1+NDET,1),NDETQ) ! add to imaginary part
           CALL DGEMM('N','N',NDET,N_ROOTS,N_CRED,-1.0D0,  ! i*i = -1
     &        WORK(KBCVEC+NDET),NDETQ,     ! imaginary part
     &        WORK(KEIGVCI),MAX_HRED,      ! imaginary part
     &        1.0D0,CVECS,NDETQ)           ! subtract from real part
         END IF

         ! 2) save them for restart
           ! hjaaj: ought to also save the corresponding SVECS ...
           ! will save one CI iteration at restart ..

         REWIND (LUGASCIP)
         WRITE (LUGASCIP) NCMOTQ, NDET, NDETQ, N_ROOTS
         CALL WRITT(LUGASCIP,N_ROOTS*NDETQ,CVECS)

         ! 3) multiply with -ECIA(i_root) (note: real eigenvalues) and save in BCVECS
         JBCVEC = KBCVEC + N_CRED*NDETQ - 1
         DO I_ROOT = 1, N_ROOTS
           WORK(JBCVEC+1:JBCVEC+NDETQ) =
     &        -ECIA(I_ROOT)*CVECS(1:NDETQ,I_ROOT)
           JBCVEC = JBCVEC + NDETQ
         END DO
         ! 1b) the imaginary part of eigen vectors C(:,:)

         ! 4) add SUM_J C(J,I_ROOT) SCVEC(J) to finish residual
         ! 4a) the real part of eigen vectors C(:,:)
         JBCVEC = KBCVEC + N_CRED*NDETQ
         CALL DGEMM('N','N',NDETQ,N_ROOTS,N_CRED,1.0D0,
     &        WORK(KSCVEC),NDETQ,          ! real + imag. part
     &        WORK(KEIGVC),MAX_HRED,       ! real part
     &        1.0D0,WORK(JBCVEC),NDETQ)    ! add to real + imag. part
         ! 4b) the imaginary part of eigen vectors C(:,:)
         IF (NZ_in_CI .gt. 1) THEN
           CALL DGEMM('N','N',NDET,N_ROOTS,N_CRED,1.0D0,
     &        WORK(KSCVEC),NDETQ,          ! real part
     &        WORK(KEIGVCI),MAX_HRED,      ! imaginary part
     &        1.0D0,WORK(JBCVEC+NDET),NDETQ) ! add to imaginary part
           CALL DGEMM('N','N',NDET,N_ROOTS,N_CRED,-1.0D0,  ! i*i = -1
     &        WORK(KSCVEC+NDET),NDETQ,     ! imaginary part
     &        WORK(KEIGVCI),MAX_HRED,      ! imaginary part
     &        1.0D0,WORK(JBCVEC),NDETQ)    ! subtract from real part
         END IF

         ! 5) calculate norms of residuals and check if converged
         ICONV_TOT = 1
         JBCVEC = KBCVEC + N_CRED*NDETQ
         DO I_ROOT = 1, N_ROOTS
           DNRES(I_ROOT) = DNORM2(NDETQ,WORK(JBCVEC),1)
           IF ( DNRES(I_ROOT) .GT. THR_CONV ) ICONV_TOT = 0  ! not all converged
           JBCVEC = JBCVEC + NDETQ
         END DO ! I_ROOT
         IF (I_STATE .GT. 0) THEN ! only converge this root
           IF ( DNRES(I_STATE) .LE. THR_CONV ) ICONV_TOT = 1 ! the requested state is converged
         END IF

         WRITE (LUGASCIP) DNRES(1:N_ROOTS)
         WRITE (LUGASCIP) ECIA(1:N_ROOTS)
         REWIND (LUGASCIP)

         ! 6) print
         WRITE(LUPRI,'(/A,1P,D8.2,A/,(I10,D20.5))')
     &     ' (GASCIP_RCISTD) Norm of CI residuals (thr = ',THR_CONV,')',
     &     (I_ROOT,DNRES(I_ROOT),I_ROOT=1,N_ROOTS)
         IF (IPRINT .GE. 10) THEN
            WRITE(LUPRI,'(/A,I4,A)')
     &       ' (GASCIP_RCISTD) CI residual vectors for',N_ROOTS,' roots'
            IF (NZ_in_CI.gt.1)
     &         WRITE(LUPRI,'(/A)') '(real 1, imag 1, real 2, ...)'
            JBCVEC = KBCVEC + N_CRED*NDETQ
            CALL OUTPUT(WORK(JBCVEC),1,NDET,1,NZ_in_CI*N_ROOTS,
     &         NDET,NZ_in_CI*N_ROOTS,-1,LUPRI)
         END IF

      IF (IPRINT .GE. 3) THEN
         CALL GETTIM(CPU_end,WALL_end)
         CPUTID  = SECTID(CPU_end-CPU_st)
         WALLTID = SECTID(WALL_end-WALL_st)
         WRITE(LUPRI,'(/5A)')
     &     ' (GASCIP_RCISTD) CPU (Wall) time for residuals :',
     &     CPUTID,' (',WALLTID,')'
         CPU_st = CPU_end
         WALL_st = WALL_end
      END IF
C
         IF (ICONV_TOT .EQ. 1) THEN
            WRITE(LUPRI,'(/A)')
     &         ' *** All requested CI roots are converged. ***'
            GOTO 200
         END IF
C
         CALL FLSHFO(LUPRI)
C
         IF ( ITMIC .EQ. MAX_CI_IT ) THEN
            ICONV_TOT = -2
            GOTO 200
         END IF
         IF (N_CRED .GE. MAX_HRED) THEN
            ICONV_TOT = -3
            GOTO 200
         END IF
C
C        Construct N_ROOTS new trial vectors
C        -----------------------------------
C
         N_ROOTS_NOT_CONVERGED = 0   !! counter for new trial vectors (only add non-converged roots!)
         JBCVEC = KBCVEC + N_CRED*NDETQ - 1
         DO I_ROOT = 1,N_ROOTS
            IF ( DNRES(I_ROOT) .LE. THR_CONV ) CYCLE ! this one is converged
            N_ROOTS_NOT_CONVERGED = N_ROOTS_NOT_CONVERGED + 1
            ECIA_I = ECIA(I_ROOT)
            DO I = 1, NDET
               HD = WORK(KHDIAG+I-1) - ECIA_I
               IF (ABS(HD) .LE. THRZER) HD = 1.0D0
               WORK(JBCVEC      + I) = - WORK(JBCVEC      + I) / HD
               IF (NZ_in_CI .EQ. 2)
     &         WORK(JBCVEC+NDET + I) = - WORK(JBCVEC+NDET + I) / HD
            END DO
            IF (IPRINT .GE. 10) THEN
              write (lupri,*) 'trial vector',N_ROOTS_NOT_CONVERGED,
     &                        ' for root',i_root,' before renorm'
              call output(work(jbcvec+1),1,ndet,1,NZ_in_CI,
     &                    ndet,NZ_in_CI,1,lupri)
            END IF
            CALL RCINORM(WORK(JBCVEC+1),IPRINT)
            IF (IPRINT .GE. 10) THEN
              write (lupri,*) 'trial vector',N_ROOTS_NOT_CONVERGED,
     &                        ' for root',i_root,' after renorm'
              call output(work(jbcvec+1),1,ndet,1,NZ_in_CI,
     &                 ndet,NZ_in_CI,1,lupri)
            END IF
            JBCVEC = JBCVEC + NDETQ
         END DO  ! I_ROOT = 1, N_ROOTS

      IF (IPRINT .GE. 3) THEN
         CALL GETTIM(CPU_end,WALL_end)
         CPUTID  = SECTID(CPU_end-CPU_st)
         WALLTID = SECTID(WALL_end-WALL_st)
         WRITE(LUPRI,'(/5A)')
     &     ' (GASCIP_RCISTD) CPU (Wall) time for new trial vectors :',
     &     CPUTID,' (',WALLTID,')'
         CPU_st = CPU_end
         WALL_st = WALL_end
      END IF
C
C        Orthogonalize trial vector(s)
C        and remove linear dependent vectors
C        -----------------------------------
C
C        WRITE(LUPRI,*) 'N_CRED, N_ROOTS_NOT_CONVERGED'
C     &                 ,N_CRED, N_ROOTS_NOT_CONVERGED
         N_CSIM = 0
         DO J = N_CRED + 1, N_CRED + N_ROOTS_NOT_CONVERGED
C           loop over new trial vectors
C      write (lupri,*) 'new vector no. J',J
C
C           Project out all previous trial vectors
C
            IROUND = 0
 400        CONTINUE
C
            OVLP_imag = 0.0D0
            DO I = 1, J - 1   ! loop over previous ("old") trial vectors
C
C              Calculate <B(old)|B(new)> / <B(old)|B(old)>
C              = <B(old)|B(new)>  (as trial vectors are normed to 1)
C
               OVLP_real = DDOT(NDETQ,WORK(KBCVEC+(I-1)*NDETQ),1,
     &                                WORK(KBCVEC+(J-1)*NDETQ),1)
               IF (NZ_in_CI .gt. 1) THEN
                  OVLP_imag =
     &                DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ),1,        ! +i <B(old,real)|B(new,imag)>
     &                          WORK(KBCVEC+(J-1)*NDETQ+NDET),1)
     &              - DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ+NDET),1,   ! -i <B(old,imag)|B(new,real)>
     &                          WORK(KBCVEC+(J-1)*NDETQ),1)
               END IF
C       write (lupri,*) 'old vector no. I, OVLP',I,OVLP_real,OVLP_imag
C
C              Calculate B(new) = B(new) - OVLP * B(old)
C
               IF (ABS(OVLP_real) .GT. THRZER) THEN
                  CALL DAXPY(NDETQ,-OVLP_real,
     &                       WORK(KBCVEC+(I-1)*NDETQ),1,
     &                       WORK(KBCVEC+(J-1)*NDETQ),1)
               END IF
               IF (ABS(OVLP_imag) .GT. THRZER) THEN
                  CALL DAXPY(NDET,-OVLP_imag,
     &                       WORK(KBCVEC+(I-1)*NDETQ),1,
     &                       WORK(KBCVEC+(J-1)*NDETQ+NDET),1)
                  CALL DAXPY(NDET,+OVLP_imag,
     &                       WORK(KBCVEC+(I-1)*NDETQ+NDET),1,
     &                       WORK(KBCVEC+(J-1)*NDETQ),1)
               END IF
            END DO  ! I = 1, J - 1 (previous trial vectors)
C
            DNAFT = DNORM2(NDETQ,WORK(KBCVEC+(J-1)*NDETQ),1)
C           write(lupri,*) 'DN after', DNAFT
            IF (ABS(DNAFT) .LE. THRZER) THEN
               WRITE (LUPRI,'(/A,I4,A)') 'INFO: New CI trial vector no.'
     &            ,J-N_CRED,' removed because of linear dependency'
               CYCLE ! do not include this vector because of it is linear dependent!
            END IF

            DNAFT = 1.0D0 / DNAFT
            CALL DSCAL(NDETQ,DNAFT,WORK(KBCVEC+(J-1)*NDETQ),1)

            IF (IPRINT .GE. 3) write(lupri,'(/A,1P,D10.2)')
     &      'Threshold for removing element in trial vectors', THR_CVEC
            IF (THR_CVEC .GT. 0.0D0 .AND. IROUND .EQ. 0) THEN
C              screen away insignificant elements here ...

               I_ZERO = 0
               I_REMOVED = 0
               J_OFF = KBCVEC-1 + (J-1)*NDETQ
               DO I= 1,NDETQ
                  IF (ABS(WORK(J_OFF+I)) .LT. THR_CVEC) THEN
                     IF (WORK(J_OFF+I) .EQ. 0.0D0) THEN
                        I_ZERO = I_ZERO + 1
                     ELSE
                        WORK(J_OFF+I) = 0.0D0
                        I_REMOVED = I_REMOVED + 1
                     END IF
                  END IF
               END DO
               IF (I_REMOVED .GT. 0) THEN
                  IF (IPRINT .GE. 3)
     &            WRITE (LUPRI,'(A,I8,A,I4,A,2I10)')
     &            'Removed',I_REMOVED,
     &            ' elements in new CI trial vector no.',J-N_CRED,
     &            '; zeroes, left :', I_ZERO,NDETQ-I_REMOVED-I_ZERO
                  IROUND = IROUND + 1 ! need a 2. round because some elements has been removed
               END IF

            ELSE
               IF (DNAFT .LE. 1.0D-4) THEN
                  IROUND = IROUND + 1 ! need a 2. round for numerical precision
               ELSE
                  IROUND = 0
               END IF
            END IF
            IF (IPRINT .GE. 5) WRITE(lupri,*)
     &         'IROUND', IROUND, ' for vector', J-N_CRED
C
C
            IF (IROUND .GT. 0) THEN
                  IF (IROUND .EQ. 1) GO TO 400
                  CALL QUIT('*** ERROR in GASCIP_RCISTD ***' //
     &              'IROUND .gt. 1')
            END IF
C
            N_CSIM = N_CSIM + 1
            J_NEW  = N_CRED + N_CSIM
            IF (J_NEW .LT. J)
     &         CALL DCOPY(NDETQ,WORK(KBCVEC+(J    -1)*NDETQ),1,
     &                          WORK(KBCVEC+(J_NEW-1)*NDETQ),1)
         END DO
C      write (lupri,*) 'All trial vectors test output after Gram-Schmidt'
C                 call output(work(kbcvec),1,ndetq,1,J_NEW,
C     &                    ndetq,J_NEW,-1,lupri)
      IF (IPRINT .GE. 3) THEN
         CALL GETTIM(CPU_end,WALL_end)
         CPUTID  = SECTID(CPU_end-CPU_st)
         WALLTID = SECTID(WALL_end-WALL_st)
         WRITE(LUPRI,'(/5A)')
     &     ' (GASCIP_RCISTD) CPU (Wall) time for orthonormalization :',
     &     CPUTID,' (',WALLTID,')'
         CPU_st = CPU_end
         WALL_st = WALL_end
      END IF

C
C        Print timings
C
         IF (IPRINT .GE. 2) THEN
            CALL GETTIM(CPU4,WALL4)
            CPUTID  = SECTID(CPU4-CPU3)
            WALLTID = SECTID(WALL4-WALL3)
            WRITE(LUPRI,'(/A,I4,5A)')
     &         ' (GASCIP_RCISTD) CPU (Wall) time for iteration ',
     &         ITMIC,': ',CPUTID,' (',WALLTID,')'
         END IF
C
      GOTO 100  ! iterations loop
C
C ==========================================================================
C
  200 CONTINUE
C
C     Print timings
C
      IF (IPRINT .GE. 2) THEN
         CALL GETTIM(CPU4,WALL4)
         CPUTID  = SECTID(CPU4-CPU3)
         WALLTID = SECTID(WALL4-WALL3)
         WRITE(LUPRI,'(/A,I3,5A)')
     &        ' (GASCIP_RCISTD) CPU (Wall) time for iteration ',
     &        ITMIC,': ',CPUTID,' (',WALLTID,')'
         CPUTID  = SECTID(CPU4-CPU2)
         WALLTID = SECTID(WALL4-WALL2)
         WRITE(LUPRI,'(/6A)')
     &        ' (GASCIP_RCISTD) Total CPU (Wall) '//
     &        'time for all microiterations: ',
     &        CPUTID,' (',WALLTID,')'
      END IF
C
  210 IF (ICONV_TOT .GT. 0) THEN
         WRITE(LUPRI,'(//A)')
     &        ' (GASCIP_RCISTD) Micro iterations converged.'
      ELSE
         IF (N_CSIM .EQ. 0) THEN
            WRITE(LUPRI,'(//A)') ' (GASCIP_RCISTD) '//
     &         'WARNING: linear dependency among all new trial vectors.'
            CALL QUIT('linear dependency among all new trial vectors.')
         ELSE IF (ICONV_TOT .EQ. -2) THEN
            WRITE(LUPRI,'(//A,I4,A)') ' (GASCIP_RCISTD) '//
     &           'WARNING: maximum number of micro iterations,',
     &           MAX_CI_IT, ', is reached, CI aborted.'
         ELSE IF (ICONV_TOT .EQ. -3) THEN
            WRITE(LUPRI,'(//A,I4,A)') ' (GASCIP_RCISTD) '//
     &           'WARNING: maximum reduced space dimension of',
     &           MAX_HRED, ' reached, CI aborted.'
         ELSE
            WRITE(LUPRI,'(//A,I4/A)') ' (GASCIP_RCISTD) '//
     &        'INFO: micro iterations not converged.',
     &        ' ICONV_TOT value :',ICONV_TOT
            ! CALL QUIT('micro iterations aborted while not converged')
         END IF
      END IF
C
      MAX_CI_IT = ITMIC
      ECI(1:N_ROOTS) = ECORE + ECIA(1:N_ROOTS)
      WRITE(LUPRI,'(/A,F20.8)')
     &     ' (GASCIP_RCISTD) Core energy        = ',ECORE
      WRITE(LUPRI,'(A,4F20.8/,(37X,4F20.8))')
     &     '                 Active energies    = ',ECIA(1:N_ROOTS)
      WRITE(LUPRI,'(A,4F20.8/,(37X,4F20.8))')
     &     ' (GASCIP_RCISTD) Final CI energies  = ',ECI(1:N_ROOTS)

C     Final solution vectors are available in CVECS

#if GASCIP_DEBUG > 20
       IF (NZ_in_CI .eq. 2) THEN
         WRITE(LUPRI,*)
     &   '(GASCIP_RCISTD): (Real, imag) part of the N_ROOTS CI vectors'
       ELSE
         WRITE(LUPRI,*)
     &   '(GASCIP_RCISTD): Real part of the N_ROOTS CI vectors'
       END IF
       CALL OUTPUT(CVECS,1,NDET,1,N_ROOTS*NZ_in_CI,
     &             NDET,N_ROOTS*NZ_in_CI,-1,LUPRI)
#endif
C
C     Check norm of CI vectors
C
      DO I_ROOT = 1, N_ROOTS

         WRITE(LUPRI,'(//A/A,I4,A,I3/A)')
     &   ' ============================================',
     &   ' (GASCIP_RCISTD) CI solution vector',I_ROOT,' of ',N_ROOTS,
     &   ' ============================================'

         CALL RCINORM(CVECS(1,I_ROOT),IPRINT)

         CALL GASCIP_ANACI(CVECS(1,I_ROOT),NDET,IDET,THR_PCI,
     &                     WORK(KFREE),LFREE)

         IF ( IPRINT .GE. 15 ) CALL RPRCI(CVECS(1,I_ROOT),
     &      NDET,NZ_in_CI,'GASCIP', THR_PCI,LUPRI)

      END DO
C
      IF (IPRINT .GE. 0) THEN
         CALL GETTIM(CPU2,WALL2)
         CPUTID  = SECTID(CPU2-CPU1)
         WALLTID = SECTID(WALL2-WALL1)
         WRITE(LUPRI,9000) CPUTID,WALLTID
      END IF
 9000 FORMAT(/' (GASCIP_RCISTD) Total CPU (Wall) time for routine: ',
     &     A,' (',A,')')
C
      CALL MEMREL('GASCIP_RCISTD',WORK,KWORK,KWORK,KFREE,LFREE)
C
!
!     return the final eigenvalues to the incoming array (needed in the calling routine)
      call dcopy(N_ROOTS,eci,1,eci_out,1)

      close(LUGASCIP, status='KEEP')
      LUGASCIP = -1

      CALL FLSHFO(LUPRI)

      IF (PAR_GASCIP) THEN
         CALL GASCIP_RELEASE_NODES()
      END IF
      CALL QEXIT('GASCIP_RCISTD')

      !> debug...
!     iprint      = iprint_save
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck idet2str */
      SUBROUTINE IDET2STR(IA,IB,JASTR,JBSTR,IAELEC,IBELEC)
C***********************************************************************
C
C     Transform bit packed IA and IB into alpha- and beta string.
C
C     Input :
C        JDET - determinant
C
C     Output:
C        JASTR,  JBSTR  - alpha and beta strings
C        IAELEC, IBELEC - alpha and beta electrons.
C
C     Written by J. Thyssen - Jan 5 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
C
      INTEGER*8 IA, IB
      DIMENSION JASTR(*), JBSTR(*)
C
      LOGICAL LBIT
C
C
      IAELEC = 0
      IBELEC = 0
C
      DO J = 1, NASHT
         IF (LBIT(IA,J)) THEN
            IAELEC = IAELEC + 1
            JASTR(IAELEC) = IPTA2O(J)
         END IF
         IF (LBIT(IB,J)) THEN
            IBELEC = IBELEC + 1
            JBSTR(IBELEC) = IPTB2O(J)
         END IF
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_anaci */
      SUBROUTINE GASCIP_ANACI(CREF,NDET,IDET,THR_PCI_in,WORK,LWORK)
C***********************************************************************
C
C     Print nice analysis of CI vector.
C
C     Input :
C        CREF - CI ref. vector
C        IDET - determintants
C
C     Output:
C
C     Written by J. Thyssen - Jan 5 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from include files:
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgascip.h"
C
      DIMENSION CREF(NDET,*), WORK(*)
      INTEGER*8 IDET(2,*)
C
      INTEGER JASTR(2*MAXASH_GASCIP), JBSTR(2*MAXASH_GASCIP)
C
      CALL QENTER('GASCIP_ANACI')
C
      THR_PCI_local = max(THR_PCI_in,1.0D-15) * 1.0000001D0
C
      DMX = 1.0D10 ! for print of sigma vectors
      DMN = 0.1D0
 10   CONTINUE
         WRITE(LUPRI,'(/A,1P,D7.1,A,D7.1,A/)')
     &     ' Print of CI coefficients with weight in range ',
     &     DMN,' to ',DMX,' :'
         K = 0
         DO I = 1, NDET
            DI = CREF(I,1) * CREF(I,1)
            IF (NZ_in_CI.eq.2)   DI = DI + CREF(I,2) * CREF(I,2)
            IF (DI .GT. DMN .AND. DI .LE. DMX) THEN
               K = K + 1
               IF (NZ_in_CI .EQ. 1) THEN
!                 WRITE(LUPRI,'(3X,A,I8,A,1P,D20.10,A,0P,F10.8)')
                  WRITE(LUPRI,'(3X,A,I8,A,F20.10,A,F10.8)')
     &                 'Determinant no. ',I,' coeff. ',
     &                 CREF(I,1),'; weight ',DI
               ELSE
!                 WRITE(LUPRI,
!    &                 '(3X,A,I8,A,1P,D20.10,A,D20.10,A,0P,F10.8)')
                  WRITE(LUPRI,
     &                 '(3X,A,I8,2(A,F20.10),A,F10.8)')
     &                 'Determinant no. ',I,' coeff. ',
     &                 CREF(I,1),' + ',CREF(I,2),' i; weight ',DI
               END IF
               CALL IDET2STR(IDET(1,I),IDET(2,I),JASTR,JBSTR,IA,IB)
               WRITE(LUPRI,'(3X,A,40I4)')
     &            '  P-string     : ',(JASTR(J),J=1,IA)
               WRITE(LUPRI,'(3X,A,40I4)')
     &            '  P(bar)-string: ',(JBSTR(J)-NASHT,J=1,IB)
            END IF
         END DO
         IF (K .EQ. 0) WRITE(LUPRI,'(3X,A/)')
     &      'No CI coefficients in this range.'
         DMX = DMN
         DMN = DMN * 0.1D0
      IF (DMX .GE. THR_PCI_local) GO TO 10
      ! make sure we do not miss a desired interval because of round-off
C
      WRITE(LUPRI,9100)
      N_NONZERO_TOT = 0
      DO MK2 = MINMK2, MAXMK2
C        ... MK2 is NAEL - NBEL
         IF (NMK2(MK2,3) .EQ. 0) CYCLE
         I_STRT = IMK2(MK2,3) + 1
         I_END  = IMK2(MK2,3) + NMK2(MK2,3)
C
         N_NONZERO = 0
         DO I = I_STRT, I_END
            IF (CREF(I,1) .NE. 0.0D0 .OR. CREF(I,NZ_in_CI) .NE. 0.0D0)
     &         N_NONZERO = N_NONZERO + 1
         END DO ! I
         N_NONZERO_TOT = N_NONZERO_TOT + N_NONZERO
C
         IF (MOD(MK2,2) .EQ. 0) THEN
            WRITE(LUPRI,9101) MK2/2,NMK2(MK2,3),N_NONZERO
         ELSE
            WRITE(LUPRI,9102) MK2,NMK2(MK2,3),N_NONZERO
         END IF
      END DO ! MK2
      WRITE(LUPRI,9103) NDET,N_NONZERO_TOT
 9100 FORMAT(//' Distribution of determinants wrt. MK'/1X,48('-')
     &        /6X,'MK',T23,'#dets',T34,'#non-zero'/1X,48('-'))
 9101 FORMAT(I8,T18,I10,T33,I10)
 9102 FORMAT(I6,'/2',T18,I10,T33,I10)
 9103 FORMAT(1X,48('-') / 5X,'all',T18,I10,T33,I10/1X,48('-'))
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('GASCIP_ANACI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_chkkr */
      SUBROUTINE GASCIP_CHKKR(BVEC,WORK,LWORK)
C***********************************************************************
C
C     Check that BVEC is Kramers restricted.
C
C     Input :
C       BVEC: CI vector
C
C     Output:
C
C     Written by J. Thyssen - Jan 23 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
C
      DIMENSION BVEC(*), WORK(*)
C
#include "memint.h"
C
      call quit('NDET not defined in GASCIP_CHKKR')
      IF ( MOD(NAELEC,2) .NE. 0) RETURN
C
      CALL MEMGET2('INTE','RLIST',KRLIST,NDET,WORK,KFREE,LFREE)
      CALL IREAKRMC(LUKRMC,'KRLIST  ',WORK(KRLIST),NDET)
      CALL GASCIP_CHKKR1(BVEC,NDET,WORK(KRLIST))
      CALL MEMREL('GASCIP_CHKKR',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_chkkr1 */
      SUBROUTINE GASCIP_CHKKR1(BVEC,NDET,KRLIST)
C***********************************************************************
C
C     Check that BVEC is Kramers restricted.
C
C     Input :
C       BVEC: CI vector
C       KRLIST: list of determinants and their Kramers partner.
C
C     Output:
C
C     Written by J. Thyssen - Jan 23 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
C
      DIMENSION BVEC(NDET,NZ_in_CI), KRLIST(*)
C
      DO I = 1, NDET
         IK = KRLIST(I)
         IF (IK .NE. I) THEN
            BVEC_DIFF = ABS(BVEC(I,1) - BVEC(IK,1))
ChjAaJ oct 08: first absolute check, if not OK, then relative check
C              (if not both, we might divide by zero!)
            IF ( BVEC_DIFF .GT. 1.0D-13 ) THEN
            IF ( BVEC_DIFF/ABS(BVEC(I,1)) .GT. 1.0D-13 ) THEN
               write(LUPRI,'(A,2I4,3D25.10)')
     &            ' GASCIP_CHKKR1: real diff',
     &            I,IK,BVEC(I,1),BVEC(IK,1),BVEC_DIFF
            END IF
            END IF
            IF (NZ_in_CI .eq. 2) THEN
               BVEC_DIFF = ABS(BVEC(I,2) + BVEC(IK,2))
ChjAaJ oct 08: first absolute check, if not OK, then relative check
C              (if not both, we might divide by zero!)
               IF ( BVEC_DIFF .GT. 1.0D-13 ) THEN
               IF ( BVEC_DIFF/ABS(BVEC(I,1)) .GT. 1.0D-13 ) THEN
                    write(LUPRI,'(A,2I4,3D25.10)')
     &                 ' GASCIP_CHKKR1: imag diff',
     &                 I,IK,BVEC(I,2),BVEC(IK,2),BVEC_DIFF
               END IF
               END IF
            END IF
         END IF
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_genstr */
      SUBROUTINE GASCIP_GENSTR(IEL,JSTR,NSTR,MSKGAS,MAXSTR)
C***********************************************************************
C
C     Generate all strings with IEL electrons and GAS constraints
C     fullfilled.
C
C     Input :
C       IEL: number of electrons
C       MSKGAS: bit mask for GAS spaces
C
C     Output:
C       JSTR: strings
C       NSTR: #strings
C
C     Written by J. Thyssen - Jan 23 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"

#include "dcborb.h"
!from dcborb.h NAELEC, NGAS_DC

#include "dcbgascip.h"
C
      INTEGER*8 JSTR(*), MSKGAS(*)
C
      !//INTEGER*4 nGASSpc(2,NGAS_DC), nGASTOT(NGAS_DC)
      INTEGER   nGASSpc(2,NGAS_DC), nGASTOT(NGAS_DC)
C
C     Calculate accumulated minimum and maximum number of
C     A electrons in each gas space. Total IEL A electrons.
C     NGASSP(1:2,i) give accumulated minimum and maximum number of
C     all electrons in gas space i.
C
#if GASCIP_DEBUG > 5
      write (lupri,*) 'gascip_genstr: iel, naelec,ngas_dc',
     &    iel,naelec,ngas_dc
#endif
      NGTOT = 0
      DO I = 1, NGAS_DC
         NGTOT  = NGTOT + NGSHT(I)
         nGASTOT(I) = NGTOT
C
         MAX_A_EL = MIN(       IEL, NGASSP(2,I), NGTOT)
         MAX_B_EL = MIN(NAELEC-IEL, NGASSP(2,I), NGTOT)
         ! max A/B electrons is min of
         !   total number of A/B electrons and
         !   max total acccumulated number of electrons for this GAS space and
         !   total number of available A/B orbitals
         MIN_A_EL = MAX(NGASSP(1,I) - MAX_B_EL,0)
         ! min A electrons is max of
         !   min total accumulated number of electrons minus max B electrons and
         !   zero

#if GASCIP_DEBUG > 5
         write (lupri,*) 'GAS space, constraints, NGTOT',
     &      I,NGASSP(1,I),NGASSP(2,I),NGTOT
         write (lupri,*) 'min_a_el,max_a_el,max_b_el',
     &                    min_a_el,max_a_el,max_b_el
#endif

         IF ( MAX_A_EL .LT. MIN_A_EL ) THEN
            ! IEL A electrons is not compatible with GAS constraint
#if GASCIP_DEBUG > 5
            WRITE (LUPRI,'(/A/A,3I10)')
     &       "This is not compatible with GAS constraints!",
     &       "gas space, MAX_A_EL, MIN_A_EL",I,MAX_A_EL,MIN_A_EL
#endif
            NSTR = 0
            GO TO 9000
         END IF
         nGASSpc(1,I) = MIN_A_EL
         nGASSpc(2,I) = MAX_A_EL
#if GASCIP_DEBUG > 5
         write(lupri,*)
     &       'gascip_genstr: gas, ngsht, nGASsp(1:2), nGASSpc(1:2) : ',
     &      I,ngsht(I),ngassp(1,i),ngassp(2,i),nGASSpc(1,I),nGASSpc(2,I)
#endif
      END DO
C
      CALL cgendet(IEL,NASHT,NSTR,JSTR,
     &             NGAS_DC,nGASTOT,MSKGAS,nGASSpc,MAXSTR)
#if GASCIP_DEBUG > 5
      write (lupri,*) 'After cgendet ..., NSTR ',NSTR
#endif
C
 9000 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_ci */
      SUBROUTINE GASCIP_CI(
     &     DODIA,DIAG,
     &     DOSIG,SVEC,BVEC,HCORE,VMUUUU,
     &     DODEN,CL,CR,DV,PV,
     &     DO2EL,NDET,IDET,IHRM,IPRINT)
C***********************************************************************
C
C     Calculate:
C     (1) diagonal of CI Hamiltonian (DODIA)
C     (2) CI sigma vector (DOSIG)
C     (3) Density matrices (DODEN)
C
C     Input:
C        DODIA, DOSIG, DODEN - control flags
C
C        BVEC - trial vector
C        HCORE  - core Fock matrix
C        VMUUUU - two electron integrals.
C
C        CL, CR - left and right vector for density matrices
C
C        DO2EL - include 2-electron part
C                -- DODIA and DOSIG: include VMUUUU
C                -- DODEN: calculate PV
C        IDET  - determinant strings
C
C     Output:
C        DIAG   - diagonal of Hamiltonian (if SAVDIA)
C
C        SVEC   - sigma vector
C
C        DV, PV - density matrices
C
C     Written by J. Thyssen and H.J.Aa. Jensen - Jan-Mar 2001
C     Last revision : July 2002/Sep 2003 hjaaj
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "thrzer.h"
C
C Used from include files:
C
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION DIAG(*)
C
      DIMENSION SVEC(NDET,*), BVEC(NDET,*)
      DIMENSION VMUUUU(*), HCORE(*)
C
      DIMENSION CL(NDET,*), CR(NDET,*)
      DIMENSION DV(2*NASHT,2*NASHT,*)
      DIMENSION PV(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      LOGICAL   DODIA, DOSIG, DODEN, DO2EL, DO_DENCHK
C
      INTEGER*8 IDET(2,*), IDA, IDB, JDA, JDB, IDAOLD, JDAOLD
      INTEGER*8 IDAV, IDBV, IAI, IBI, IAJ, IBJ
C
      LOGICAL   DOIMAG, DOIMAG1
C
      DIMENSION IOCCI(2*MAXASH_GASCIP), JOCCI(2*MAXASH_GASCIP)
      DIMENSION IBTAB(0:255)
      INTEGER*8 FNCALLS, CNCALLS
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
C
      CALL GETTIM(CPU1,WALL1)

      IF (IPRINT .GT. 1) WRITE (LUPRI,'(//A/A,3L10/A,L10,2I10)')
     &   ' --- Output from GASCIP_CI ---',
     &   ' DODIA, DOSIG, DODEN :',DODIA,DOSIG,DODEN,
     &   ' DO2EL, IHRM,  IPRINT:',DO2EL,IHRM,IPRINT
C
C     Initialization
C     --------------
C
      IF (IHRM .eq. 0) THEN
         CALL QUIT('IRHM = 0, NOT IMPLENTED YET IN GASCIP_CI')
      END IF

      NDETQ = NDET*NZ_in_CI
C
C     Bit table:
C
      CALL IBTABINI(IBTAB)
C
C     Zero counters
C
      CALL GETNBITSCALLS(1,FNCALLS,CNCALLS)
C
C     Zero vectors:
C
      IF (DODIA) CALL DZERO(DIAG,NDET)
      IF (DOSIG) CALL DZERO(SVEC,NDETQ)
      IF (DODEN) THEN
         CALL DZERO(DV,(2*NASHT)**2 * NZ_in_CI)
         IF (DO2EL) CALL DZERO(PV,(2*NASHT)**4 * NZ_in_CI)
         IF (IHRM .le. 0)
     &      CALL QUIT('GASCIP_CI: IHRM.le.0 not implemented for DODEN')
      END IF

#if GASCIP_DEBUG > 20
       IF (DOSIG) WRITE(LUPRI,*)
     &   'GASCIP_CI: This CI vector'
       IF (DOSIG) CALL OUTPUT(BVEC,1,NDET,1,NZ_in_CI,
     &             NDET,NZ_in_CI,1,LUPRI)
      write(lupri,*) 'DO2el, IHRM are ' , do2el,ihrm
#endif
C
      IF (DO2EL) THEN
C        max 4 differences between 2 determinants
C        for 2-electron terms
         ITMAX = 4
      ELSE
C        max 2 differences between 2 determinants
C        when only 1-electron terms
         ITMAX = 2
      END IF
C
      DOIMAG  = NZ.GT.1
      DOIMAG1 = DOIMAG
C
C     screening for density matrices:
C
      DO_DENCHK = DODEN .AND. .NOT.DOSIG .AND. .NOT.DODIA
      THRSCR = 1.0D-10
      NDSKIP = 0
C
      MK2I_STEP = 2 ! remeber that MK = MK2 / 2, step of 1 for MK is thus step of 2 for MK2
C
      DO 10 MK2I = MINMK2, MAXMK2, MK2I_STEP
C        ... MK2 is NAEL - NBEL
         IF (NMK2(MK2I,3) .EQ. 0) GOTO 10
C        write(lupri,*) '(GASCIP_CI) MK2I = ',MK2I
         ISTRT = IMK2(MK2I,3) + 1
         IEND  = IMK2(MK2I,3) + NMK2(MK2I,3)
C
         IF (NZ .LE. 2 .AND. ITMAX .EQ. 4) THEN
            MK2J_STEP = 4 ! only Delta(MK) = 2, i.e. Delta(MK2) = 4, couple for non-quaternion point groups
         ELSE
            MK2J_STEP = 2 ! this is a quaternion point group (NZ = 4), Delta(MK2) = 2 will also couple
         END IF
C
C        set MK2J limits for the off-diagonal blocks with MK2I .gt. MK2J
         MK2J_STRT = MK2I-ITMAX
         IF (MK2J_STRT .LT. MINMK2) MK2J_STRT = MK2J_STRT + MK2J_STEP
         IF (MK2J_STRT .LT. MINMK2) MK2J_STRT = MK2J_STRT + MK2J_STEP
         MK2J_END  = MK2I - MK2J_STEP ! because MK2J .lt. MK2I
C
C
C           Loop over determinants
C           (use that outer loop in gascip_gendet is over A-strings,
C           inner loop is over B-strings)
C
            IDAOLD = -1
            DO 100 I = ISTRT, IEND
               IF (DO_DENCHK) THEN
                  DENCHK = ABS(CL(I,1)) + ABS(CR(I,1))
                  IF (DOIMAG)
     &            DENCHK = DENCHK + ABS(CL(I,2)) + ABS(CR(I,2))
                  IF (DENCHK .LT. THRSCR) THEN
                     NDSKIP = NDSKIP + 1
                     GO TO 100
                  END IF
               ENDIF
C
C
               IDA = IDET(1,I)
               IF (IDA .NE. IDAOLD) THEN
C                 ... outer loop in gascip_gendet is over A-strings
                  IAELEC = 0
                  DO K = 1, NASHT
                     IF ( BTEST(IDA, K - 1) ) THEN
                        IAELEC = IAELEC + 1
                        IOCCI(IAELEC) = IPTA2O(K)
                     END IF
                  END DO
                  IDAOLD = IDA
               END IF
C
               IELEC = IAELEC
               IDB = IDET(2,I)
               DO K = 1, NASHT
                  IF (BTEST(IDB, K - 1) ) THEN
                     IELEC = IELEC + 1
                     IOCCI(IELEC) = IPTB2O(K)
                  END IF
               END DO
               IF (IELEC .NE. NAELEC)
     &            CALL QUIT('GASCIP_CI: IELEC.ne.NAELEC')
C
C        *************************
C        *** Off-diagonal loop *** (MK2J .lt. MK2I)
C        *************************
C
C
            IF (.NOT. (DODEN .OR. DOSIG)) GOTO 21
C           ...only off-diagonal blocks for dens.mat or sigma-vec.
C
            DO 20 MK2J = MK2J_STRT, MK2J_END, MK2J_STEP
               IF (NMK2(MK2J,3) .EQ. 0) GOTO 20
C              write(lupri,*) '(GASCIP_CI) MK2I MK2J = ',MK2I,MK2J
C              call flshfo(lupri)
C
               JSTRT = IMK2(MK2J,3) + 1
               JEND  = IMK2(MK2J,3) + NMK2(MK2J,3)
               MXNDAV = ITMAX + (MK2J - MK2I)/2
C              ... we know that NDBV .ge. (MK2I-MK2J)/2
C
               JDAOLD = -1
               DO 200 J = JSTRT, JEND
C
                  JDA = IDET(1,J)
                  IF (JDA .NE. JDAOLD) THEN
C                    ... outer loop in gascip_gendet is over A-strings
                     IDAV = IEOR(IDA, JDA)
                     CALL FNBITS(NDAV,IDAV,IBTAB)
                     IF (NDAV .GT. MXNDAV) GOTO 200
C
                     JAELEC = 0
                     IAI = IAND(IDAV, IDA)
                     KA  = 0
                     IAJ = IAND(IDAV, JDA)
                     MA  = 0
                     DO II = 1, NASHT
                        IF ( BTEST(JDA, II-1) ) THEN
                           JAELEC = JAELEC + 1
                           JOCCI(JAELEC) = IPTA2O(II)
                        END IF
                        IF ( BTEST(IAI, II-1) ) THEN
                           L  = KA
                           KA = IPTA2O(II)
                        END IF
                        IF ( BTEST(IAJ, II-1) ) THEN
                           N  = MA
                           MA = IPTA2O(II)
                        END IF
                     END DO
                     JDAOLD = JDA
                  END IF
                  IF (NDAV .GT. MXNDAV) GOTO 200
C                 ... see note above about MXNDAV
C
                  JDB  = IDET(2,J)
                  IDBV = IEOR(IDB, JDB)
                  CALL FNBITS(NDBV,IDBV,IBTAB)
C
                  IT = NDAV + NDBV
                  IF (IT .GT. ITMAX) GOTO 200
C                 ... IT = number of different spinors
C                          in det(I) and det(J)
                  IF (DO_DENCHK) THEN
                     DENCHK = ABS(CL(J,1)) + ABS(CR(J,1))
                     IF (DOIMAG)
     &               DENCHK = DENCHK + ABS(CL(J,2)) + ABS(CR(J,2))
                     IF (DENCHK .LT. THRSCR) GO TO 200
                  ENDIF
C
                  IELEC = JAELEC
                  IBI = IAND(IDBV, IDB)
                  K = KA
                  IBJ = IAND(IDBV, JDB)
                  M = MA
                  DO II = 1, NASHT
                     IF ( BTEST(JDB, II-1) ) THEN
                        IELEC = IELEC + 1
                        JOCCI(IELEC) = IPTB2O(II)
                     END IF
C
                     IF ( BTEST(IBI, II-1) ) THEN
                        L = K
                        K = IPTB2O(II)
                     END IF
                     IF ( BTEST(IBJ, II-1) ) THEN
                        N = M
                        M = IPTB2O(II)
                     END IF
                  END DO
C
C                 Hamiltonian matrix element
C                 --------------------------
C
                  IF (DOSIG) THEN
C
                     CALL GASCIP_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,
     &                    K,L,M,N,HR,HI,HCORE,VMUUUU)
C
C                    S_I += H_{IJ} B_{J}
C
                     SVEC(I,1) = SVEC(I,1) + HR * BVEC(J,1)
                     IF (DOIMAG) THEN
                        SVEC(I,1) = SVEC(I,1) - HI * BVEC(J,2)
                        SVEC(I,2) = SVEC(I,2)
     &                       + HR * BVEC(J,2) + HI * BVEC(J,1)
                     END IF
C
C                    S_J += H_{IJ}^* B_{I}
C
                     IF (IHRM .lt. 0) THEN
                        HR = -HR
                        HI = -HI
                     END IF
                     SVEC(J,1) = SVEC(J,1) + HR * BVEC(I,1)
                     IF (DOIMAG) THEN
                        SVEC(J,1) = SVEC(J,1) + HI * BVEC(I,2)
                        SVEC(J,2) = SVEC(J,2)
     &                       + HR * BVEC(I,2) - HI * BVEC(I,1)
                     END IF
                  END IF
C
C                 Density matrix elements
C                 -----------------------
C
                  IF (DODEN) THEN
C
                     DMRIJ = CL(I,1)*CR(J,1)
                     DMRJI = CL(J,1)*CR(I,1)
                     IF (DOIMAG) THEN
                        DMRIJ = DMRIJ + CL(I,2)*CR(J,2)
                        DMIIJ = CL(I,1)*CR(J,2) - CL(I,2)*CR(J,1)
                        DMRJI = DMRJI + CL(J,2)*CR(I,2)
                        DMIJI = CL(J,1)*CR(I,2) - CL(J,2)*CR(I,1)
                        DOIMAG1 = (ABS(DMIIJ)+ABS(DMIJI)) .GT. THRZER
                     END IF
                     CALL GASCIP_DENMAT(DO2EL,DOIMAG1,IT,IOCCI,JOCCI,
     &                    K,L,M,N,
     &                    DMRIJ,DMIIJ,DMRJI,DMIJI,DV,PV)
C
                  END IF
 200           CONTINUE
 20         CONTINUE
 21         CONTINUE
C
C
C
C        ***********************
C        *** Diagonal blocks *** (MK2J .eq. MK2I)
C        ***********************
C
C           Diagonal element (Det_I .eq. DET_J)
C
C           Calculate Hamiltonian matrix element
C           ------------------------------------
C
            IF (DODIA .OR. DOSIG) THEN
C
               CALL GASCIP_HMAT(DO2EL,.FALSE.,0,IOCCI,IOCCI,0,0,0,0,
     &              HR,DUMMY,HCORE,VMUUUU)
C
C              Sigma-vector elements
C
C              SVEC_I += H_{II} B_{I}
C
               IF (DOSIG) THEN
C
                  SVEC(I,1) = SVEC(I,1) + HR * BVEC(I,1)
                  IF (DOIMAG) SVEC(I,2) = SVEC(I,2) + HR * BVEC(I,2)
C
               END IF
C
               IF (DODIA) DIAG(I) = HR
C
            END IF
C
C           Calculate density matrix elements
C           ---------------------------------
C
            IF (DODEN) THEN
C
               IF (.NOT. DOIMAG) THEN
                  DMR = CL(I,1)*CR(I,1)
               ELSE
                  DMR = CL(I,1)*CR(I,1) + CL(I,2)*CR(I,2)
                  DMI = CL(I,1)*CR(I,2) - CL(I,2)*CR(I,1)
                  DOIMAG1 = ABS(DMI) .GT. THRZER
               END IF
               CALL GASCIP_DENMAT(DO2EL,DOIMAG1,0,IOCCI,IOCCI,0,0,0,0,
     &              DMR,DMI,DMR,DMI,DV,PV)
C
            END IF
C
C           Off-diagonal elements
C           ---------------------
C
            IF (.NOT. (DOSIG .OR. DODEN)) GOTO 100
            JDAOLD = -1
            DO 400 J = ISTRT, I - 1
C
               JDA = IDET(1,J)
               IF (JDA .NE. JDAOLD) THEN
C                 ... outer loop in gascip_gendet is over A-strings
                  IDAV = IEOR(IDA, JDA)
                  CALL FNBITS(NDAV,IDAV,IBTAB)
                  IF (NDAV .GT. ITMAX) GOTO 400
C
                  JAELEC = 0
                  IAI = IAND(IDAV, IDA)
                  KA  = 0
                  IAJ = IAND(IDAV, JDA)
                  MA  = 0
                  DO II = 1, NASHT
                     IF ( BTEST(JDA, II-1) ) THEN
                        JAELEC = JAELEC + 1
                        JOCCI(JAELEC) = IPTA2O(II)
                     END IF
                     IF ( BTEST(IAI, II-1) ) THEN
                        L  = KA
                        KA = IPTA2O(II)
                     END IF
                     IF ( BTEST(IAJ, II-1) ) THEN
                        N  = MA
                        MA = IPTA2O(II)
                     END IF
                  END DO
                  JDAOLD = JDA
               END IF
               IF (NDAV .GT. ITMAX) GOTO 400
C
               IF (DO_DENCHK) THEN
                  DENCHK = ABS(CL(J,1)) + ABS(CR(J,1))
                  IF (DOIMAG)
     &               DENCHK = DENCHK + ABS(CL(J,2)) + ABS(CR(J,2))
                  IF (DENCHK .LT. THRSCR) GO TO 400
               ENDIF
C
               JDB = IDET(2,J)
               IF (NDAV .EQ. ITMAX) THEN
                  IF (JDB .NE. IDB) GO TO 400
                  IT = NDAV
                  IBI = 0
                  IBJ = 0
                  K = KA
                  M = MA
                  DO IELEC = JAELEC + 1, NAELEC
                     JOCCI(IELEC) = IOCCI(IELEC)
                  END DO
               ELSE
                  IDBV = IEOR(IDB, JDB)
                  CALL FNBITS(NDBV,IDBV,IBTAB)
                  IT = NDAV + NDBV
                  IF (IT .GT. ITMAX) GOTO 400
C
                  IELEC = JAELEC
                  IBI = IAND(IDBV, IDB)
                  IBJ = IAND(IDBV, JDB)
                  K = KA
                  M = MA
                  DO II = 1, NASHT
                     IF ( BTEST(JDB, II-1) ) THEN
                        IELEC = IELEC + 1
                        JOCCI(IELEC) = IPTB2O(II)
                     END IF
C
                     IF ( BTEST(IBI, II-1) ) THEN
                        L = K
                        K = IPTB2O(II)
                     END IF
                     IF ( BTEST(IBJ, II-1) ) THEN
                        N = M
                        M = IPTB2O(II)
                     END IF
                  END DO
               END IF
C
C              Hamiltonian matrix element
C              --------------------------
C
               IF (DOSIG) THEN
C
                  CALL GASCIP_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
     &                 HR,HI,HCORE,VMUUUU)
C
C                 S_I += H_{IJ} B_{J}
C
                  SVEC(I,1) = SVEC(I,1) + HR * BVEC(J,1)
                  IF (DOIMAG) THEN
                     SVEC(I,1) = SVEC(I,1) - HI * BVEC(J,2)
                     SVEC(I,2) = SVEC(I,2)
     &                    + HR * BVEC(J,2) + HI * BVEC(J,1)
                  END IF
C
C                 S_J += H_{IJ}^* B_{I}
C
                  IF (IHRM .lt.0) THEN
                     HR = - HR
                     HI = - HI
                  END IF
                  SVEC(J,1) = SVEC(J,1) + HR * BVEC(I,1)
                  IF (DOIMAG) THEN
                     SVEC(J,1) = SVEC(J,1) + HI * BVEC(I,2)
                     SVEC(J,2) = SVEC(J,2)
     &                    + HR * BVEC(I,2) - HI * BVEC(I,1)
                  END IF
               END IF
C
C              Density matrix elements
C              -----------------------
C
               IF (DODEN) THEN
C
                  DMRIJ = CL(I,1)*CR(J,1)
                  DMRJI = CL(J,1)*CR(I,1)
                  IF (DOIMAG) THEN
                     DMRIJ = DMRIJ + CL(I,2)*CR(J,2)
                     DMIIJ = CL(I,1)*CR(J,2) - CL(I,2)*CR(J,1)
                     DMRJI = DMRJI + CL(J,2)*CR(I,2)
                     DMIJI = CL(J,1)*CR(I,2) - CL(J,2)*CR(I,1)
                     DOIMAG1 = (ABS(DMIIJ)+ABS(DMIJI)) .GT. THRZER
                  END IF
                  CALL GASCIP_DENMAT(DO2EL,DOIMAG1,IT,IOCCI,JOCCI,
     &                 K,L,M,N,DMRIJ,DMIIJ,DMRJI,DMIJI,DV,PV)
C
               END IF
 400        CONTINUE
 100     CONTINUE
C
 10   CONTINUE
C
      IF (IPRINT .GT. 1) THEN
        CALL GETTIM(CPU2,WALL2)
        IF ((WALL2-WALL1)*IPRINT .GT. 5.D0) THEN
C        ... do not print all the small numbers
          CPUTID  = SECTID(CPU2-CPU1)
          WALLTID = SECTID(WALL2-WALL1)
          WRITE(LUPRI,'(/5A)') ' CPU (Wall) time in GASCIP_CI kernel: ',
     &      CPUTID,' (',WALLTID,')'
          IF (NDSKIP .GT. 0) WRITE(LUPRI,'(A,I10,2(A,I10))')
     &      ' Number of configuration elements skipped:',NDSKIP,
     &      ' out of',NDET,'; remaining ',NDET-NDSKIP
          IF (IPRINT .GT. 2) THEN
            CALL GETNBITSCALLS(1,CNCALLS,FNCALLS)
            write(lupri,*)'#calls to cnbits and fnbits:',cncalls,fncalls
          END IF
        END IF
      END IF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_cisig_screen */
      SUBROUTINE GASCIP_CISIG_SCREEN(IS_START, IS_END,
     &                               SVEC,BVEC,HCORE,VMUUUU,
     &                               DO2EL,NDET,IDET,THRSCR,IPRINT)
C***********************************************************************
C
C     Calculate:
C         CI sigma vector in SVEC with screening
C
C     Input:
C
C        IS_START, IS_END - see below
C        BVEC   - trial vector
C        HCORE  - core Fock matrix
C        VMUUUU - two electron integrals.
C
C        DO2EL  - include 2-electron part : VMUUUU
C        IDET   - determinant strings
C        THRSCR - skip abs(BVEC(i)) .le. THRSCR
C
C     Output:
C
C        SVEC(IS:START,IS_END,NZ_in_CI)   - sigma vector
C
C
C     Written by J. Thyssen and H.J.Aa. Jensen - Jan-Mar 2001
C     Last revision : July 2002/Sep 2003/Aug 2014 hjaaj
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from include files:
C  infpar.h : MYTID
C
#include "dcbgascip.h"
#include "dcborb.h"
#include "infpar.h"
#include "dgroup.h"
C
      DIMENSION SVEC(NDET,*), BVEC(NDET,*)
      DIMENSION VMUUUU(*), HCORE(*)
C
C
      LOGICAL   DO2EL
C
      INTEGER*8 IDET(2,*), IDA, IDB, JDA, JDB, IDAOLD, JDAOLD
      INTEGER*8 IDAV, IDBV, IAI, IBI, IAJ, IBJ
C
      LOGICAL   DOIMAG, NEW_HMAT
C
      DIMENSION IOCCI(2*MAXASH_GASCIP), JOCCI(2*MAXASH_GASCIP)
      DIMENSION IBTAB(0:255)
      INTEGER*8 FNCALLS, CNCALLS
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
C
      CALL GETTIM(CPU1,WALL1)

      IF (IS_START .LT. 1 .OR. IS_END .GT. NDET) THEN
         WRITE (LUPRI,'(//A,I6/A,L10,1P,D10.2,3I12)')
     &   ' --- Interval out of bounds in GASCIP_CISIG_SCREEN --- node',
     &   MYTID,
     &   ' DO2EL, THRSCR, NDET, interval:',
     &   DO2EL,THRSCR,NDET,IS_START,IS_END
         CALL QUIT('IS_START or IS_END is '//
     &      'out of bounds in GASCIP_CISIG_SCREEN')
      END IF

      IF (IPRINT .GT. 1 .AND. MYTID .LE. 1)
     &   WRITE (LUPRI,'(//A,I6/A,L10,1P,D10.2,2I12)')
     &   ' --- Output from GASCIP_CISIG_SCREEN --- node',MYTID,
     &   ' DO2EL, THRSCR, interval:',DO2EL,THRSCR,IS_START,IS_END
C
C     Initialization
C     --------------
C
      NDETQ = NDET*NZ_in_CI
C
C     Bit table:
C
      CALL IBTABINI(IBTAB)
C
C     Zero counters
C
      CALL GETNBITSCALLS(1,FNCALLS,CNCALLS)
C
C     Zero vectors:
C
      CALL DZERO(SVEC,NDETQ)
C
      IF (DO2EL) THEN
C        max 4 differences between 2 determinants
C        for 2-electron terms
         ITMAX = 4
      ELSE
C        max 2 differences between 2 determinants
C        when only 1-electron terms
         ITMAX = 2
      END IF
C
      DOIMAG  = NZ.GT.1
C
      MK2I_STEP = 2 ! remeber that MK = MK2 / 2, step of 1 for MK is thus step of 2 for MK2
C
      NBSKIP = 0
      DO 10 MK2I = MINMK2, MAXMK2, MK2I_STEP
C        ... MK2 is NAEL - NBEL
         IF (NMK2(MK2I,3) .EQ. 0) GOTO 10
         ISTRT = IMK2(MK2I,3) + 1
         IEND  = IMK2(MK2I,3) + NMK2(MK2I,3)
C
         IF (ITMAX .EQ. 2) THEN
            MK2J_STEP = 2   ! 1-electron operators are often non-totally symmetric
                            ! thus they may couple with Delta(MK2) = 2
            MK2J_STRT = MAX(MINMK2,MK2I-ITMAX)
            MK2J_END  = MIN(MAXMK2,MK2I+ITMAX)
         ELSE ! we assume the 2-electron operators are totally symmetric
              ! (for non-totally symmetric operators we would need MK2J_STEP = 2 TODO FIXME)
            IF (NZ .LE. 2) THEN
               MK2J_STEP = 4 ! only Delta(MK) = 2, i.e. Delta(MK2) = 4, couple for non-quaternion point groups
            ELSE
               MK2J_STEP = 2 ! this is a quaternion point group (NZ = 4), Delta(MK2) = 2 will also couple
            END IF
            MK2J_STRT = MK2I-ITMAX
            IF (MK2J_STRT .LT. MINMK2) MK2J_STRT = MK2J_STRT + MK2J_STEP
            IF (MK2J_STRT .LT. MINMK2) MK2J_STRT = MK2J_STRT + MK2J_STEP
            MK2J_END  = MIN(MAXMK2,MK2I+ITMAX)
         END IF
C
C        Loop over determinants
C        (use that outer loop in gascip_gendet is over A-strings,
C        inner loop is over B-strings)
C
         IDAOLD = -1
!d       write(lupri,*) '(GASCIP_CISIG_SCREEN), MK2I, ISTRT,IEND',
!d   &                   MK2I,ISTRT,IEND
!d       write(lupri,*) '(GASCIP_CISIG_SCREEN), MK2J_STRT/END/STEP',
!d   &                   MK2J_STRT,MK2J_END,MK2J_STEP

         DO 100 I = ISTRT, IEND
C
C        Screening:
C
            IF (ABS(BVEC(I,1)) .GT. THRSCR) GOTO 12
            IF (DOIMAG) THEN
               IF (ABS(BVEC(I,2)) .GT. THRSCR) GOTO 12
            END IF
            NBSKIP = NBSKIP + 1
            GOTO 100
C
C           ---------------
C
   12       IDA = IDET(1,I)
            IF (IDA .NE. IDAOLD) THEN
C              ... outer loop in gascip_gendet is over A-strings
               IAELEC = 0
               DO K = 1, NASHT
                  IF ( BTEST(IDA, K - 1) ) THEN
                     IAELEC = IAELEC + 1
                     IOCCI(IAELEC) = IPTA2O(K)
                  END IF
               END DO
               IDAOLD = IDA
            END IF
C
            IELEC = IAELEC
            IDB = IDET(2,I)
            DO K = 1, NASHT
               IF ( BTEST(IDB, K - 1) ) THEN
                  IELEC = IELEC + 1
                  IOCCI(IELEC) = IPTB2O(K)
               END IF
            END DO
            IF (IELEC .NE. NAELEC)
     &         CALL QUIT('GASCIP_CISIG_SCREEN: IELEC.ne.NAELEC')
C
C
C        *************************
C        *************************
C
C
            DO 20 MK2J = MK2J_STRT, MK2J_END, MK2J_STEP
               IF (NMK2(MK2J,3) .EQ. 0) GOTO 20
C
               JSTRT = IMK2(MK2J,3) + 1
               JEND  = IMK2(MK2J,3) + NMK2(MK2J,3)
               MXNDAV = ITMAX - ABS(MK2J - MK2I)/2
C              ... we know that NDBV .ge. ABS(MK2I-MK2J)/2
C
!d          write(lupri,*) '(GASCIP_CISIG_SCREEN) MK2I MK2J JSTRT JEND='
!d   &         ,MK2I,MK2J,JSTRT,JEND
!d             call flshfo(lupri)
               JDAOLD = -1
               JSTRT = MAX(JSTRT, IS_START)
               JEND  = MIN(JEND,  IS_END)
               DO 200 J = JSTRT, JEND
C
                  JDA = IDET(1,J)
                  IF (JDA .NE. JDAOLD) THEN
                     JDAOLD = JDA
C                    ... outer loop in gascip_gendet is over A-strings
                     IDAV = IEOR(IDA, JDA)
                     CALL FNBITS(NDAV,IDAV,IBTAB)
                     IF (NDAV .GT. MXNDAV) GOTO 200
                     NEW_HMAT = .TRUE. ! new HR, HI from GASCIP_HMAT
C
                     JAELEC = 0
                     IAI =  IAND(IDAV, IDA)
                     KA  = 0
                     IAJ =  IAND(IDAV, JDA)
                     MA  = 0
                     DO II = 1, NASHT
                        IF ( BTEST(JDA, II-1) ) THEN
                           JAELEC = JAELEC + 1
                           JOCCI(JAELEC) = IPTA2O(II)
                        END IF
                        IF ( BTEST(IAI, II-1) ) THEN
                           L  = KA
                           KA = IPTA2O(II)
                        END IF
                        IF ( BTEST(IAJ, II-1) ) THEN
                           N  = MA
                           MA = IPTA2O(II)
                        END IF
                     END DO
                  END IF
                  IF (NDAV .GT. MXNDAV) GOTO 200
C                 ... see note above about MXNDAV
C
                  JDB  = IDET(2,J)
                  IF (NDAV .EQ. ITMAX) THEN
                     IF (JDB .NE. IDB) GO TO 200
                     IT = NDAV
                     IBI = 0
                     IBJ = 0
                     K = KA
                     M = MA
                     DO IELEC = JAELEC + 1, NAELEC
                        JOCCI(IELEC) = IOCCI(IELEC)
                     END DO
                  ELSE
                     IDBV = IEOR(IDB, JDB)
                     CALL FNBITS(NDBV,IDBV,IBTAB)
C
                     IT = NDAV + NDBV
                     IF (IT .GT. ITMAX) GOTO 200
C                    ... IT = number of different spinors
C                          in det(I) and det(J)
C
                     NEW_HMAT = .TRUE.
                     IELEC = JAELEC
                     IBI =  IAND(IDBV, IDB)
                     IBJ =  IAND(IDBV, JDB)
                     K = KA
                     M = MA
                     DO II = 1, NASHT
                        IF ( BTEST(JDB, II-1) ) THEN
                           IELEC = IELEC + 1
                           JOCCI(IELEC) = IPTB2O(II)
                        END IF
C
                        IF ( BTEST(IBI, II-1) ) THEN
                           L = K
                           K = IPTB2O(II)
                        END IF
                        IF ( BTEST(IBJ, II-1) ) THEN
                           N = M
                           M = IPTB2O(II)
                        END IF
                     END DO
                  END IF
C
C                 New Hamiltonian matrix element?
C                 -------------------------------
C
                  IF (NEW_HMAT) THEN
                     CALL GASCIP_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,
     &                    K,L,M,N,HR,HI,HCORE,VMUUUU)
                     NEW_HMAT = .FALSE.
                  END IF
C
C                 S_J += H_{IJ}^* B_{I}
C
                  SVEC(J,1) = SVEC(J,1) + HR * BVEC(I,1)
                  IF (DOIMAG) THEN
                     SVEC(J,1) = SVEC(J,1) + HI * BVEC(I,2)
                     SVEC(J,2) = SVEC(J,2)
     &                    + HR * BVEC(I,2) - HI * BVEC(I,1)
                  END IF
C              end loop 200 J and end loop 20 MK2J
 200           CONTINUE
 20         CONTINUE
C        end loop 100 I and end loop 10 MK2I
 100     CONTINUE
 10   CONTINUE
C     ------------------------------------------------
      IF (IPRINT .GT. 0) THEN
         CALL GETTIM(CPU2,WALL2)
       IF ((WALL2-WALL1)*IPRINT .GT. 5.D0) THEN
C        ... do not print all the small numbers
         CPUTID  = SECTID(CPU2-CPU1)
         WALLTID = SECTID(WALL2-WALL1)
         WRITE(LUPRI,'(/5A/A,I10,A,I10,A,I10)')
     &   ' CPU (Wall) time in GASCIP_CISIG_SCREEN kernel: ',CPUTID,' (',
     &   WALLTID,')',' Number of BVEC elements skipped:',NBSKIP,
     &   ' out of',NDET,'; remaining ',NDET-NBSKIP
         IF (NUMNOD .GT. 0) THEN
            WRITE(LUPRI,'(A,I5,A,2I15)')
     &      ' Node',MYTID,'; start and end in SVEC:',IS_START,IS_END
         END IF
        IF (IPRINT .GT. 2) THEN
         CALL GETNBITSCALLS(1,CNCALLS,FNCALLS)
         write(lupri,*) '#calls to cnbits and fnbits:',cncalls,fncalls
        END IF
       END IF
      END IF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck gascip_inp */
      SUBROUTINE GASCIP_INP(WORD,RESET,WORK,LWORK)
!***********************************************************************
!
!     Input section for GASCIP-module
!
!     Written by: S. Yamamoto + tsaue, Feb 06, 2007
!
!     Copied from KRMCINP and modified.
!     Common block data will be modified in this routine.
!     The data is stored in the common block in the cossya.h file.
!     Called by: DIRCTL->PAMINP->PSIINP->GASCIP_INP
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
!
      PARAMETER (NTABLE = 31)
!
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK),IMKSET(3)
      DIMENSION MULTB_TMP(64,64)
!
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcborb.h"
#include "dcbgascip.h"
!
!     ...The variable "SET" must have save attribute.
      SAVE SET
      DATA TABLE /'.PRINT ','.GASSH ','.GASSPC','.MK2REF','.MK2DEL',
     &            '.INACTI','.CASSCF','.CI SYM','.RELSYM','.MK2LIM',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX'/
      DATA SET/.FALSE./
!
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
C
C     Initialize DC_GASCIP_ORB
C     ================
C
      IPRGASCIP = 0
      IGASCIP_REP = -1
      DO IFRP = 1,NFSYM
        IGASCIP_NISH(IFRP) = 0
        IGASCIP_NASH(IFRP) = 0
      ENDDO
      IGASCIP_STATE = 1
      GASCIP_UCIBOS = .FALSE.
      CALL IZERO(IGASCIP_NGSHT,MAXASH_GASCIP)
C
C     Initialize DCBORB
C     =================
C     NOTE that we use the parts of it not used by HF !
C     We thus assume no MCSCF in the same run...
C
      IGASCIP_NGAS = 0
      DO I = 1,2
         IGASCIP_NISH(I) = 0
         IGASCIP_NASH(I) = 0
         DO J = 1, MAXASH_GASCIP
            IGASCIP_NGSH(I,J) = 0
         END DO
      END DO
      IGASCIP_NAELEC = -1
C
C     Process input
C     =============
C
      CALL IZERO(IMKSET,3)
      NEWDEF = (WORD .EQ. '*GASCIP')
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24,25,26,27,28,29,30,
     &                     31), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized for *GASCIP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword under *GASCIP.')
!
    1          CONTINUE
!&&&  PRINT : print level in GASCIP module
                  READ(LUCMD,*) IPRGASCIP
               GO TO 100
!
    2          CONTINUE
!&&&  GASSH : GAS setup
                  IF (IGASCIP_NGAS.GT.0) THEN
                     CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                 'You cannot specify both .GASSH and .CASSCF')
                  END IF
                  READ(LUCMD,*) IGASCIP_NGAS
                  IF (IGASCIP_NGAS .LT. 1 .OR.
     &                IGASCIP_NGAS .GT. MXGAS) THEN
                     CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                         'illegal value of IGASCIP_NGAS')
                  END IF
                  DO I = 1, IGASCIP_NGAS
                     READ(LUCMD,*) (IGASCIP_NGSH(J,I),J=1,NFSYM)
                  END DO
               GO TO 100
!
    3          CONTINUE
C&&&  GASSPC: GAS space constraints
                  IF (IGASCIP_NGAS .EQ. 0) THEN
                     WRITE(LUPRI,'(A)') ' *** ERROR in GASCIP_INP ***'
     &                  //' .GASSH must be specified before .GASSP!'
                     CALL QUIT('*** ERROR in GASCIP_INP ***')
                  END IF
                  DO I = 1, IGASCIP_NGAS
                     READ(LUCMD,*) (IGASCIP_NGASSP(J,I),J=1,2)
                  END DO
                  IF (IGASCIP_NAELEC.LT.0) IGASCIP_NAELEC =
     &               IGASCIP_NGASSP(2,IGASCIP_NGAS)
               GO TO 100
!
    4          CONTINUE
C&&&  MK2REFERENCE: 2 * M_K reference
                  READ(LUCMD,*) IGASCIP_MK2REF
                  IMKSET(2) = 1
               GO TO 100
    5          CONTINUE
C&&&  MK2DELTA: 2 * DELTA M_K
                  READ(LUCMD,*) IGASCIP_MK2DEL
                  IMKSET(3) = 1
               GO TO 100
    6          CONTINUE
C&&&  INACTIVE ORBITALS
               READ(LUCMD,*) (IGASCIP_NISH(I),I=1,NFSYM)
               GO TO 100
    7          CONTINUE
C&&&  CASSCF -- this defines a CAS space
               IF (IGASCIP_NGAS .GT. 0) THEN
                  CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                 'You cannot specify both .GASSH and .CASSCF')
               END IF
               READ(LUCMD,*) IGASCIP_NAELEC
               READ(LUCMD,*) (IGASCIP_NASH(I),I=1,NFSYM)
               IGASCIP_NGAS = 1
               IGASCIP_NGASSP(1,IGASCIP_NGAS) = IGASCIP_NAELEC
               IGASCIP_NGASSP(2,IGASCIP_NGAS) = IGASCIP_NAELEC
    8          CONTINUE
C&&&  SYMMETRY: symmetry of the wave function
               IF (IGASCIP_REP .GT. 0) THEN
                  CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                 'You can only specify symmetry once')
               END IF
               READ(LUCMD,*) IGASCIP_REP
               GO TO 100
!
    9          CONTINUE
C&&&  RELSYM: symmetry of the wave function
!.s/sya,2007.02.06
               IF (IGASCIP_REP .GT. 0) THEN
                  CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                 'You can only specify symmetry once')
               END IF
               READ(LUCMD,*) IGASCIP_REP
!.q
               GO TO 100
!
   10          CONTINUE
C&&&  MK2LIM: min/max value of 2 * MK
               READ(LUCMD,*) IGASCIP_MINMK2,IGASCIP_MAXMK2
               IMKSET(1) = 1
               GO TO 100
!
   11          CONTINUE
   12          CONTINUE
   13          CONTINUE
   14          CONTINUE
   15          CONTINUE
   16          CONTINUE
   17          CONTINUE
   18          CONTINUE
   19          CONTINUE
   20          CONTINUE
   21          CONTINUE
   22          CONTINUE
   23          CONTINUE
   24          CONTINUE
   25          CONTINUE
   26          CONTINUE
   27          CONTINUE
   28          CONTINUE
   29          CONTINUE
   30          CONTINUE
   31          CONTINUE
!&&&  XXXXXX: invalid data
                  CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                      'invalid data XXXXXX')
               GO TO 100
!
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in GASCIP_INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in GASCIP_INP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process input
C     =============
C
      IF (IGASCIP_NGAS .GT. 0) THEN
         IGASCIP_NASHT = 0
         DO J = 1, NFSYM
            IGASCIP_NASH(J) = 0
            DO I = 1, IGASCIP_NGAS
               IGASCIP_NASH(J) = IGASCIP_NASH(J) + IGASCIP_NGSH(J,I)
            END DO
            IGASCIP_NASHT = IGASCIP_NASHT + IGASCIP_NASH(J)
         END DO
      END IF
      IF (IGASCIP_NASHT.GT.MAXASH_GASCIP) THEN
         WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &      ' *** ERROR in GASCIP_INP *** Too many active orbitals: ',
     &      IGASCIP_NASHT, ' -- max is ',MAXASH_GASCIP
         CALL QUIT(
     &      '*** ERROR in GASCIP_INP *** Too many active orbitals!')
      ENDIF
CTROND: set UCIBOS !!!!
C
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' Set-up for *GASCIP calculation'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I4)') ' * General print level   : ',IPRGASCIP
      WRITE(LUPRI,'(A,2I4)')
     &     ' * Inactive orbitals     : ',(IGASCIP_NISH(I),I=1,NFSYM)
      WRITE(LUPRI,'(A,2I4)')
     &     ' * Active orbitals       : ',(IGASCIP_NASH(I),I=1,NFSYM)
      WRITE(LUPRI,'(A,I4)')
     &     ' * Active electrons      : ',IGASCIP_NAELEC
      IF (IGASCIP_NGAS .GT. 0) THEN
         WRITE(LUPRI,'(A,I3,A)')
     &        ' * GAS space setup for ',IGASCIP_NGAS,' GAS space(s) : '
         DO I = 1, IGASCIP_NGAS
            WRITE(LUPRI,'(A,I3,A,2I4)')
     &           '   - GAS space ',I,'       : ',
     &           (IGASCIP_NGSH(J,I),J=1,NFSYM)
            WRITE(LUPRI,'(A,I3,A1,I3,A)')
     &      '    (constraints: min/max active electrons after space : ',
     &           IGASCIP_NGASSP(1,I),'/',IGASCIP_NGASSP(2,I),')'
         END DO
         IF (IGASCIP_NGASSP(2,IGASCIP_NGAS) .NE. IGASCIP_NAELEC)
     &        CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &        'Inconsistency between .GASSP and .CASSCF')
      END IF
C
C     Symmetry of the wave function:
C     ------------------------------
C
C
C     MK constraints in the CI expansion.
C
      IF (IMKSET(1) .EQ. 1) THEN
C
C        MINMK2/MAXMK2 explicitly set
C
      ELSE
         IF (IMKSET(2) .EQ. 1 .AND. IMKSET(3) .EQ. 1) THEN
C
C           MK2REF/MK2DEL set, transform to MINMK2/MAXMK2
C
            IGASCIP_MINMK2 = IGASCIP_MK2REF - IGASCIP_MK2DEL
            IGASCIP_MAXMK2 = IGASCIP_MK2REF + IGASCIP_MK2DEL
C
         ELSE
C
C           Assign default values for MINMK2/MAXMK2
C
            IGASCIP_MAXMK2 = 2*MIN(IGASCIP_NAELEC,IGASCIP_NASHT)
     &                     - IGASCIP_NAELEC
            IGASCIP_MINMK2 = IGASCIP_NAELEC
     &                     - 2*MIN(IGASCIP_NAELEC,IGASCIP_NASHT)
            IGASCIP_MK2REF = 0
            IGASCIP_MK2DEL = IGASCIP_MAXMK2
         END IF
      END IF
C
C     Boson symmetry in the CI expansion
C
      IF(IGASCIP_REP.GT.0) THEN
        IF ((SPINFR.or.levyle).AND.GASCIP_UCIBOS) THEN
          NREP = NBSYM
          WRITE(LUPRI,'(A)')
     &    ' * Boson spatial spinor symmetry used in CI expansion.'
          WRITE(LUPRI,'(2A)')
     &    ' * Spatial symmetry (boson) of wave function : ',
     &    REP( IGASCIP_REP - 1 )
          WRITE(LUPRI,'(A,I4,A,I3)')
     &    ' * Allowed interval of 2 * MS :',
     &    IGASCIP_MINMK2, ' to ', IGASCIP_MAXMK2
        ELSE IF (.NOT.(SPINFR.or.levyle).AND.GASCIP_UCIBOS) THEN
          NREP = NBSYM
          WRITE(LUPRI,'(A)') ' * Approximate boson spatial spinor'//
     &    ' symmetry (as non-rel) used in CI expansion.'
          WRITE(LUPRI,'(2A)')
     &    ' * Approximate spatial symmetry (boson) of wave function : ',
     &    REP ( IGASCIP_REP - 1 )
          WRITE(LUPRI,'(A,I4,A,I3)')
     &    ' * Allowed interval of approximate 2 * MS :',
     &    IGASCIP_MINMK2, ' to ', IGASCIP_MAXMK2
        ELSE
          IF (SPINFR.or.levyle) THEN
            CALL GMULTSF(NREP,REPNA,MULTB_TMP)
          ELSE
            CALL GMULTA(NZ,NFSYM,NREP,REPNA,MULTB_TMP)
          ENDIF
          IF ( MOD ( IGASCIP_NAELEC, 2 ) .EQ. 0 ) THEN
            IRRP = NREP + IGASCIP_REP
            WRITE(LUPRI,'(2A)')
     &        ' * Symmetry (boson) of wave function: ',
     &        REPNA(IRRP)
          ELSE
            WRITE(LUPRI,'(2A)')
     &        ' * Symmetry (fermion) of wave function: ',
     &        REPNA(IGASCIP_REP)
          END IF
          WRITE(LUPRI,'(A,I4,A,I3)')
     &     ' * Allowed interval of 2 * MK :',
     &     IGASCIP_MINMK2, ' to ', IGASCIP_MAXMK2
        END IF
      ELSE
CTROND: imporive on this,...
        WRITE(LUPRI,'(A)')
     &    '* All symmetries used'
      ENDIF
C
      IF (IGASCIP_REP .GT. NREP) THEN
            WRITE(LUPRI,'(/A//A,2I3)')
     &         ' INPUT ERROR, requested .SYMMETRY is non-existent!',
     &         ' symmetry, max value of symmetry = ',
     &           IGASCIP_REP, NREP
            CALL QUIT('INPUT ERROR for *GASCIP  .SYMMETRY')
      END IF
C
C
      CALL GASCIP_COPY
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_copy */
      SUBROUTINE GASCIP_COPY
C***********************************************************************
C
C     Transfer information from common block in dcbgascip.h to
C     common block in dcbkrmc.h.
C
C     Written by: S. Yamamoto - Fri Feb  9 15:26:44 CET 2007
C     Last revision:
C     Called by: GASCIP_INP
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxgas.h"
C
#include "dgroup.h"
#include "dcbgascip.h"
#include "dcbkrmc.h"
C
C
      return ! we try ...
      CALL QENTER('GASCIP_COPY')
      CALL QUIT('GASCIP_COPY called, see comments subroutine')
C     hjaaj: NO, NO, NO! NKMRC variables are only for
C            MCSCF optimization and nothing else.
C            The correct implementation would be to
C            add a TGASCIPORB routine in dirac.F which
C            transfers the IGASCIP* information to dcborb.
C            Also in dirac.F, if DOGASCIP true then
C            TGASCIPORB should be called (same procedure as
C            TKRCIORB for DOKRCI etc.) followed by call of
C            GASCIP.  /hjaaj 14-Jul-2011
C
C
      NKRMCGAS = IGASCIP_NGAS
      NKRMC_MK2REF = IGASCIP_MK2REF
      NKRMC_MK2DEL = IGASCIP_MK2DEL
!radovan: IGASCIP_MINMK2 is not initialized
      NKRMC_MINMK2 = IGASCIP_MINMK2
!radovan: IGASCIP_MAXMK2 is not initialized
      NKRMC_MAXMK2 = IGASCIP_MAXMK2
ccc   CALL IZERO(NKRMC_GSHT,MXGAS)
      DO I = 1,NFSYM
         NKRMCISH(I) = IGASCIP_NISH(I)
         NKRMCASH(I) = IGASCIP_NASH(I)
         DO J = 1, NKRMCGAS
            NKRMCGSH(I,J) = IGASCIP_NGSH(I,J)
ccc         NGSHT(J) = NGSHT(J) + NGSH(I,J)
         END DO
ccc      NOCC(I) = IGASCIP_NISH(I) + IGASCIP_NASH(I)
      END DO
      NKRMCAELEC = IGASCIP_NAELEC
C     min and max electrons in each gas space:
      DO J = 1, NKRMCGAS
         NKRMCGSP(1,J) = IGASCIP_NGASSP(1,J)
         NKRMCGSP(2,J) = IGASCIP_NGASSP(2,J)
      END DO
C
!radovan: IPRGASCIP is not initialized
      IPRKRMC = IPRGASCIP
      IKRMC_STATE = 1
      IKRMC_SYMMETRY = IGASCIP_REP
C
      KRMC_UCIBOS = GASCIP_UCIBOS
      KRMC_CIPROGRAM = 'GASCIP'
C
      CALL QEXIT('GASCIP_COPY')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     MPI parallel master routines for GASCIP
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GASCIP_WAKE_NODES()
C
C     Written Aug. 2014 by Hans Joergen Aa. Jensen
C
C     Start GASCIP_NODE driver for GASCIP on slaves
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#ifdef VAR_MPI

      use interface_to_mpi

#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "infpar.h"
C
C send info to slaves with MPI_BCAST from these common blocks:
C
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"

      IF (.NOT.PARCAL) RETURN

C
C     Send ITASK = 11 for parallel GASCIP
C     ==================================
C
#if GASCIP_DEBUG > 0
      write (lupri,*)
     &  'GASCIP_WAKE_NODES: Calling DIRAC_PARCTL( GASCIP_PAR )'
      call flshfo(lupri)
#endif
      CALL DIRAC_PARCTL( GASCIP_PAR )
C
C     Send essential common block information
C     ==============================================
C
C     dcbgascip.h
C     --------
C
      NUMELM = ICOMMSIZE(I1_DCI_GASCIP,I2_DCI_GASCIP)
      call interface_MPI_BCAST(I1_DCI_GASCIP,NUMELM,MPARID,
     &                         global_communicator)
C
C     dcborb.h
C     --------
C
      NUMELM = ICOMMSIZE(I1_DCBORB,I2_DCBORB)
      call interface_MPI_BCAST(I1_DCBORB,NUMELM,MPARID,
     &                         global_communicator)
C
C     dgroup.h
C     --------
C
      NUMELM = ICOMMSIZE(I1_DGROUP,I2_DGROUP)
      call interface_MPI_BCAST(I1_DGROUP,NUMELM,MPARID,
     &                         global_communicator)
      call interface_mpi_bcast_l0(LINEAR,1,MPARID,
     &                         global_communicator)

#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GASCIP_RELEASE_NODES()
C
C     Written Aug. 2014 by Hans Joergen Aa. Jensen
C
C     master: release slaves from GASCIP tasks
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#ifdef VAR_MPI
      use interface_to_mpi

#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "infpar.h"
C
      IF (PARCAL) THEN
         itask = -1
         call interface_MPI_BCAST(itask,1,MPARID,
     &                         global_communicator)
         CALL DIRAC_PARCTL( RELEASE_NODES )
      END IF
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GASCIP_PAR_CISIG_SETUP(
     &               HCORE,VMUUUU,
     &               DO2EL,NDET,IDET,THRSCR,IPRINT)
C
C     Written Aug. 2014 by Hans Joergen Aa. Jensen
C
C     master: start parallel CI sigma calculation and send
C             CI information and one- and two-electron operator matrices
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#ifdef VAR_MPI
      use interface_to_mpi

#include "implicit.h"
#include "priunit.h"
C
#include "infpar.h"
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
C
      REAL*8    HCORE(*), VMUUUU(*)
      LOGICAL   DO2EL
      INTEGER*8 IDET(*)

!     Start parallel GASCIP CI sigma calculation on slaves

      itask = 1
      call interface_MPI_BCAST(itask,1,MPARID,
     &                         global_communicator)

!     Send CI information and operator matrices to slaves

      call interface_mpi_bcast_l0(DO2EL,1,MPARID,
     &                         global_communicator)

      call interface_MPI_BCAST(NDET,1,MPARID,
     &                         global_communicator)
      call interface_MPI_BCAST(THRSCR,1,MPARID,
     &                         global_communicator)
      call interface_MPI_BCAST(IPRINT,1,MPARID,
     &                         global_communicator)

      len_HCORE = (2*NASHT)**2 * NZ_in_CI
      call interface_MPI_BCAST_r1_work_f77( HCORE,len_HCORE,MPARID,
     &                         global_communicator)
      IF (DO2EL) THEN
         len_VMUUUU = (2*NASHT)**4 * NZ_in_CI
         call interface_MPI_BCAST_r1_work_f77( VMUUUU,len_VMUUUU,MPARID,
     &                         global_communicator)
      END IF
      len_IDET = 2*NDET
      call interface_MPI_BCAST_r1_work_f77( IDET,len_IDET,MPARID,
     &                         global_communicator)
#if GASCIP_DEBUG > 0
      write(lupri,*) 'Master, finished GASCIP_PAR_CISIG_SETUP'
      call flshfo(LUPRI)
#endif
#endif

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GASCIP_PAR_CISIG(N_CSIM,NDET,SCVEC, BCVEC, IPRINT)
C
C     Written Aug. 2014 by Hans Joergen Aa. Jensen
C
C     master: 1) send bc-vectors to all slaves
C             2) receive sc_vector pieces from slaves
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#ifdef VAR_MPI
      use interface_to_mpi

#include "implicit.h"

      INTEGER ISTAT(df_MPI_STATUS_SIZE)

#include "priunit.h"
C
#include "infpar.h"
#include "dgroup.h"
C
      REAL*8    SCVEC(NDET,NZ_in_CI,N_CSIM), BCVEC(*)
      INTEGER   SLAVE_S_INTERVAL(3)

      IF (N_CSIM .LE. 0) GO TO 9000

      IF (IPRINT .GE. 4) THEN
         WRITE(LUPRI,'(/A)')
     &   ' Timing output from GASCIP_PAR_CISIG. '//
     &   'Only receive sections using more than 30 seconds are printed.'
      END IF

C     master: 1) send bc-vectors to all slaves
      CALL GETTIM(CPU1,WALL1)

      call interface_MPI_BCAST(
     &     N_CSIM,1,MPARID,global_communicator)

      len_BCVEC = N_CSIM * NDET * NZ_in_CI
      call interface_MPI_BCAST_r1_work_f77(
     &     BCVEC,len_BCVEC,MPARID,global_communicator)

!     master: 2) receive sc_vector pieces from slaves
      CALL GETTIM(CPU2,WALL2)
      IF (IPRINT .GE. 4) THEN
         CALL TIMDAT
         WRITE(LUPRI,'(A,2F10.2)')
     &   'Master: CPU and WALL time for broadcasting CI trial vectors',
     &   CPU2-CPU1, WALL2-WALL1
      END IF

!     Select number of tasks
!     (1: currently fixed to 5 times number of slaves, hope this
!      will give OK load balancing)
!!    (2: could still be better, now set length to 2000 instead,
!!     independent of numer of nodes)
!!    N_TASKS = 5 * NUMNOD
!!    LEN_SCVEC_PER_TASK = (NDET - 1) / N_TASKS + 1
!!    IF (LEN_SCVEC_PER_TASK .LT. 2000) THEN
         LEN_SCVEC_PER_TASK = 2000
         N_TASKS = (NDET - 1) / LEN_SCVEC_PER_TASK + 1
         N_TASKS = MAX(N_TASKS, NUMNOD)
         LEN_SCVEC_PER_TASK = (NDET - 1) / N_TASKS + 1
!!    END IF

      call interface_MPI_BCAST(
     &     N_TASKS,1,MPARID,global_communicator)

      NEXT_TASK = NUMNOD + 1 ! nodes start automatically with their first task

      DO I_TASK = 1, N_TASKS
         CAll interface_mpi_RECV(SLAVE_S_INTERVAL,3,
     &                 df_mpi_any_source,40,global_communicator,ISTAT)
         J_START = SLAVE_S_INTERVAL(1)
         NWHO = ISTAT(df_MPI_SOURCE)
        IF (IPRINT .GE. 4) THEN
         CALL GETTIM(CPU3,WALL3)
        IF ((WALL3-WALL2) .GT. 30.D0) THEN
         CALL TIMDAT
         WRITE(LUPRI,'(A,I6,A,I6,2F10.2)')
     &   'Master task',I_TASK,
     &   ': CPU and WALL time to start of receive from node',
     &   NWHO,CPU3-CPU2, WALL3-WALL2
        END IF
        END IF
         DO I = 1, N_CSIM
            CAll interface_mpi_RECV(
     &                 SCVEC(J_START,1,I),SLAVE_S_INTERVAL(3),
     &                 NWHO,40,global_communicator)
         IF (NZ_in_CI .eq. 2) THEN
            CAll interface_mpi_RECV(
     &                 SCVEC(J_START,2,I),SLAVE_S_INTERVAL(3),
     &                 NWHO,40,global_communicator)
         END IF
         END DO
        IF (IPRINT .GE. 4) THEN
         CALL GETTIM(CPU2,WALL2)
        IF ((WALL2-WALL3) .GT. 30.D0) THEN
         CALL TIMDAT
         WRITE(LUPRI,'(A,I6,A,I6,2F10.2)')
     &   'Master task',I_TASK,
     &   ': CPU and WALL time to retrieve sigma vectors from node',
     &   NWHO,CPU2-CPU3, WALL2-WALL3
        END IF
        END IF
!        send next task to this node; nodes will automatically exit task loop when
!        NEXT_TASK is greater than N_TASKS
         CAll interface_mpi_SEND(NEXT_TASK,1,NWHO,
     &      40,global_communicator)
         NEXT_TASK = NEXT_TASK + 1
      END DO

      IF (IPRINT .GE. 4) THEN
         CALL GETTIM(CPU2,WALL2)
         CALL TIMDAT
         WRITE(LUPRI,'(A,2F10.2)')
     &   'Master: Total CPU and WALL time in GASCIP_PAR_CISIG',
     &   CPU2-CPU1, WALL2-WALL1
      END IF

 9000 CONTINUE
#endif
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     MPI parallel slave routines for GASCIP
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GASCIP_NODE()
C
C     Written Aug. 2014 by Hans Joergen Aa. Jensen
C
C     Slave node driver for GASCIP
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#ifdef VAR_MPI
      use memory_allocator
      use interface_to_mpi

#include "implicit.h"
#include "priunit.h"
C
#include "infpar.h"
C
C common blocks to receive info from MPI_BCAST:
C
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
C
      real(8), allocatable :: WORK(:)
C
      CALL QENTER('GASCIP_NODE')
C
C     Allocate WORK
C     =============
C
      call legacy_lwork_get(LWORK)
      call alloc(WORK,LWORK,id='WORK in GASCIP_NODE')
#include "memint.h"
      KFRSAV = KFREE
#if GASCIP_DEBUG > 0
      write(lupri,*) 'GASCIP_NODE started, slave ',MYTID
      call flshfo(lupri)
#endif
C
C     Receive essential common block information
C     ==============================================
C
C     dcbgascip.h
C     --------
C
      NUMELM = ICOMMSIZE(I1_DCI_GASCIP,I2_DCI_GASCIP)
      call interface_MPI_BCAST(I1_DCI_GASCIP,NUMELM,MPARID,
     &                         global_communicator)
C
C     dcborb.h
C     --------
C
      NUMELM = ICOMMSIZE(I1_DCBORB,I2_DCBORB)
      call interface_MPI_BCAST(I1_DCBORB,NUMELM,MPARID,
     &                         global_communicator)
C
C     dgroup.h
C     --------
C
      NUMELM = ICOMMSIZE(I1_DGROUP,I2_DGROUP)
      call interface_MPI_BCAST(I1_DGROUP,NUMELM,MPARID,
     &                         global_communicator)
      call interface_mpi_bcast_l0(LINEAR,1,MPARID,
     &                         global_communicator)
C
C     Receive information about what GASCIP task to do
C     ================================================
C
      call interface_MPI_BCAST(itask,1,MPARID,
     &                         global_communicator)
#if GASCIP_DEBUG > 0
      write(lupri,*) MYTID,' in GASCIP_NODE: received ITASK',ITASK
      call flshfo(lupri)
#endif
      IF (itask .lt. 0) THEN
!        exit from GASCIP_NODE
         CONTINUE
      ELSE IF (itask .eq. 1) then ! GASCIP_PAR_CISG
         CALL GASCIP_NODE_CISIG(WORK,KFREE,LFREE)
      ELSE
         CALL QUIT('Unknown task code in GASCIP_NODE')
      END IF

      call dealloc(WORK)
      CALL QEXIT('GASCIP_NODE')
#endif  /* VAR_MPI */
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GASCIP_NODE_CISIG(WORK,KFREE,LFREE)
C
C     Written Aug. 2014 by Hans Joergen Aa. Jensen
C
C     Slave node driver for GASCIP CI sigma vectos
C     Called from GASCIP_NODE
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#if defined (VAR_MPI)

      use interface_to_mpi

#include "implicit.h"
#include "priunit.h"
C
#include "infpar.h"
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
C
      REAL*8    WORK(*)
      LOGICAL   DO2EL
      INTEGER   MY_S_INTERVAL(3)

      KFRSAV = KFREE

!     Receive CI information from master

      call interface_mpi_bcast_l0(DO2EL,1,MPARID,
     &                         global_communicator)

      call interface_MPI_BCAST(NDET,1,MPARID,
     &                         global_communicator)
      call interface_MPI_BCAST(THRSCR,1,MPARID,
     &                         global_communicator)
      call interface_MPI_BCAST(IPRINT,1,MPARID,
     &                         global_communicator)
      IF (MYTID .GT. 1) THEN
         IPRINT = IPRINT / 10
      END IF

      len_HCORE = (2*NASHT)**2 * NZ_in_CI
      IF (DO2EL) THEN
         len_VMUUUU = (2*NASHT)**4 * NZ_in_CI
      ELSE
         len_VMUUUU = 0
      END IF
      len_IDET = 2*NDET
      CALL MEMGET2('REAL','HCORE', KHCORE, len_HCORE,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','VMUUUU',KVMUUUU,len_VMUUUU,WORK,KFREE,LFREE)
      CALL MEMGET2('INT8','IDET',  KIDET,  len_IDET,WORK,KFREE,LFREE)

#if GASCIP_DEBUG > 0
      write(lupri,*) MYTID,' GASCIP_NODE_CISIG DO2EL,NDET,THRSCR,IPRINT'
     & ,DO2EL, NDET,THRSCR,IPRINT,NASHT
      write(lupri,*)MYTID,' lengths',len_HCORE,len_VMUUUU,len_IDET
      call flshfo(lupri)
#endif

!     Receive CI information vector IDET and operator matrices to from master

      call interface_MPI_BCAST_r1_work_f77(
     &         WORK(KHCORE),len_HCORE,MPARID,
     &         global_communicator)
      IF (DO2EL) THEN
         call interface_MPI_BCAST_r1_work_f77(
     &         WORK(KVMUUUU),len_VMUUUU,MPARID,
     &         global_communicator)
      END IF
      call interface_MPI_BCAST_r1_work_f77(
     &         WORK(KIDET),len_IDET,MPARID, ! using that integer*8 and real*8 both are 8 bytes
     &         global_communicator)

!     Get N_CSIM, number of simultaneous trial vectors

  100 CONTINUE ! coming back to here for more batches

      call interface_MPI_BCAST(N_CSIM,1,MPARID,global_communicator)
#if GASCIP_DEBUG > 0
      write(lupri,*) MYTID,' GASCIP_NODE_CISIG: received N_CSIM',N_CSIM
      call flshfo(lupri)
#endif
      IF (N_CSIM .lt. 0) GO TO 9000
      IF (N_CSIM .eq. 0) GO TO 100

      NDETQ = NDET * NZ_in_CI
      len_BCVEC = N_CSIM * NDETQ
      CALL MEMGET2('REAL','BCVEC',KBCVEC,len_BCVEC,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','SCVEC',KSCVEC,len_BCVEC,WORK,KFREE,LFREE)

!     get trial vectors BCVEC

      call interface_MPI_BCAST_r1_work_f77(
     &     WORK(KBCVEC),len_BCVEC,MPARID,global_communicator)

!     do the work

      call interface_MPI_BCAST(N_TASKS,1,MPARID,global_communicator)
      LEN_SCVEC_PER_TASK = (NDET - 1) / N_TASKS + 1

      MY_TASK = MYTID ! first task is equal to my node number

  200 CONTINUE ! coming back here for more batches
      MY_S_INTERVAL(1) = 1 + (MY_TASK-1)*LEN_SCVEC_PER_TASK
      MY_S_INTERVAL(2) = MIN(MY_TASK*LEN_SCVEC_PER_TASK, NDET)
      MY_S_INTERVAL(3) = MY_S_INTERVAL(2) - MY_S_INTERVAL(1) + 1
      IF (IPRINT .GT. 0) THEN
         WRITE (LUPRI,'(A,3I6)')
     &   ' Node, MY_TASK, N_TASKS',MYTID,MY_TASK,N_TASKS
      END IF
      DO I = 1, N_CSIM
         CALL GASCIP_CISIG_SCREEN(MY_S_INTERVAL(1), MY_S_INTERVAL(2),
     &         WORK(KSCVEC+(I-1)*NDETQ),
     &         WORK(KBCVEC+(I-1)*NDETQ),
     &         WORK(KHCORE),WORK(KVMUUUU),
     &         DO2EL,NDET,WORK(KIDET),THRSCR,IPRINT)
      END DO

!     send the sigma vector components I have calculated

      CAll interface_mpi_SEND(MY_S_INTERVAL,3,
     &                 MPARID,40,global_communicator)
      DO I = 1, N_CSIM
         J_START = KSCVEC-1+(I-1)*NDETQ+MY_S_INTERVAL(1)
         CAll interface_mpi_SEND(WORK(J_START),MY_S_INTERVAL(3),
     &                 MPARID,40,global_communicator)
         IF (NZ_in_CI .eq. 2) THEN
            J_START = J_START + NDET
            CAll interface_mpi_SEND(WORK(J_START),MY_S_INTERVAL(3),
     &                 MPARID,40,global_communicator)
         END IF
      END DO
      CAll interface_mpi_RECV(MY_TASK,1,MPARID,
     &      40,global_communicator)

      IF (MY_TASK .LE. N_TASKS) GO TO 200

      CALL MEMREL('GASCIP_NODE_CISIG.1',WORK,1,KBCVEC,KFREE,LFREE)
      GO TO 100

 9000 CONTINUE
      CALL MEMREL('GASCIP_NODE_CISIG.2',WORK,1,KFRSAV,KFREE,LFREE)
#endif   /* VAR_MPI */
      RETURN
      END
! -- end of krmcgascip.F --
