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
      SUBROUTINE PAMGOS(CI,MAXCORE,NELECT,NGAS,IGAS,
     &                  MINE,MAXE,NGASO)
C     ===============================================================
C     SMALL CI PROGRAM TO PERFORM COMPLETE OPEN SHELL CI CALCULATIONS
C     METHOD DESCRIBED IN JCP 96, 2910, (1992).
C     FIRST VERSION WRITTEN NOVEMBER 1990 BY O. VISSER
C     THEORETICAL CHEMISTRY DEPARTMENT,
C     UNIVERSITY OF GRONINGEN (RUG), GRONINGEN,THE NETHERLANDS
C     Adaption for GAS by Tommi Matila, Apr 27 2000
C     ==========================================================
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "goscom.h"
C
      DIMENSION CI(MAXCORE)
      DIMENSION MINE (NGAS),MAXE (NGAS),IGAS(NGAS)
      DIMENSION NGASO(*)
C
      CHARACTER*8 DATEX*10,TIMEX
C
      EXTERNAL INIT
      CALL DAYTIME(DATEX,TIMEX)
      WRITE(6,1000) DATEX,TIMEX
      CALL CPUUSED(SEC)
C
      CALL IBTABINI1()
C     Read namelist input
      CALL READGOS(NELECT)
C
C     Set pointers to allocate blocks in memory
C
      NRB4    = NSPACT**4
      KINDTWR = 1
      KINDTWI = KINDTWR + NRB4
      IF (BREIT) THEN
         KINDGAR = KINDTWI + NRB4  
         KINDGAI = KINDGAR + NRB4  
         KTOP    = KINDGAI + NRB4
      ELSE
         KTOP    =  KINDTWI + NRB4
         KINDGAR = KTOP  
         KINDGAI = KTOP
      ENDIF
C
      NDETMX = 1
C.....Create masks fo GAS spaces
      CALL GASMASK(NGAS,IGAS,NGASO)
      WRITE(6,'(A)') 'Number of determinants:'
C.....loop over fermion irreps I Abelian subgroup
      IERR = 0
      DO IRP=1, NREP
C.......odd number of electrons: fermion irreps
        IRRP = IRP
C.......even number of electrons: boson irreps
        IF (MOD(NELACT,2).EQ.0) IRRP = IRP + NREP
C.......Generate list of determinants
        CALL GENDET(IRP,IRRP,NGAS,IGAS,MINE,MAXE,IERR)
        WRITE(6,'(3X,A,A4,A,I10)') 'Symmetry ',REPNA(IRRP),'  : ',NDET
        NDETMX = MAX0(NDETMX,NDET)
      ENDDO
      IF(IERR.NE.0) CALL Q2JOB(3,'P.GENDET','    N2',NDETMX)
      KHR  = KTOP
      KHI  = KHR + NDETMX*NDETMX
      KTOP = KHI + NDETMX*NDETMX
      MAXINT = (MAXCORE - KTOP - 2) / 2
C
C     Check if this still fits in the available memory
C
      IF (MAXINT.LE.0)  CALL Q2JOB (3,'P.PAMGOS','MAXINT',NRB4)
C
      KTWOINR = KTOP
      KTWOINI = KTWOINR + MAXINT
      KTOP    = KTWOINI + MAXINT 
C
      CALL RDTWO(CI(KINDTWR),CI(KINDTWI),CI(KINDGAI),CI(KINDGAR),
     &           CI(KTWOINR),CI(KTWOINI),MAXINT)
      DO 10 IRP=1,NREP
        IRRP = IRP
        IF (MOD(NELACT,2).EQ.0) IRRP = IRP + NREP
        CALL CPUUSED(SEC0)
        CALL GENDET(IRP,IRRP,NGAS,IGAS,MINE,MAXE,IERR)
        IF (NDET.EQ.0) GOTO 10
        CALL CPUUSED(SEC1)
        CALL MAKEH(CI(KHR),CI(KHI),CI(KINDTWR),CI(KINDTWI),
     &       CI(KINDGAI),CI(KINDGAI),CI(KTWOINR),CI(KTWOINI),MAXINT)
        CALL CPUUSED(SEC2)
C
        IF (IPRNT .GE. 11) THEN
          CALL HEADER('Real part of total H',-1)
          CALL OUTPUT(CI(KHR),1,NDET,1,NDET,NDET,NDET,1,6)
          CALL HEADER('Imag part of total H',-1)
          CALL OUTPUT(CI(KHI),1,NDET,1,NDET,NDET,NDET,1,6)
        END IF
C
        CALL DIAGH(IRP,CI(KHR),CI(KHI))
C
        IF (IPRNT .GE. 11) THEN
          CALL HEADER('Real part of diagonalized H',-1)
          CALL OUTPUT(CI(KHR),1,NDET,1,NDET,NDET,NDET,1,6)
          CALL HEADER('Imag part of diagonalized H',-1)
          CALL OUTPUT(CI(KHI),1,NDET,1,NDET,NDET,NDET,1,6)
        END IF
C
        CALL CPUUSED(SEC3)
        CALL WRVEC(IRP,IRRP,CI(KHR),CI(KHI))
        CALL CPUUSED(SEC4)
        WRITE (6,1010) SEC1-SEC0,SEC2-SEC1,SEC3-SEC2,SEC4-SEC3
        CALL PRTOUT(IRP,IRRP,CI(KHR),CI(KHI))
 10     CONTINUE
      INTREP = MOD(NELACT,2)
      CALL PRTOUT(0,INTREP,CI(KHR),CI(KHI))
      CALL CPUUSED(SEC5)
      WRITE (6,*)    '====================='
      WRITE (6,1020) SEC5-SEC
 1000 FORMAT(///1X,'GOSCIP (Version 2.0)'//
     &' Today is :',T15,A10/' The time is :',T17,A8//)
 1010 FORMAT(/' CPU time for this representation '//
     &' Generation of determinants :',T30,F12.4/
     &' Building the CI matrix :',T30,F12.4/
     &' Diagonalization :',T30,F12.4/
     &' Writing CI vectors :',T30,F12.4)
 1020 FORMAT (//'Total CPU time :',T30,F12.4//8X,
     &          '(NORMAL END OF PROGRAM)'///)
C
C VERSION 0.1: SEEMS TO BE WORKING FOR P-SHELLS OF ATOM
C              POPULATIONS ARE GIVEN ONLY IF $POPAN IS PRESENT
C                                                   OV, 19-11-1990
C VERSION 0.2: COMPRESSED LIST OF EIGENVALUES GIVEN AT THE END OF THE OUTPUT
C                                                   OV, 17-12-1990
C VERSION 0.3: TAKE CARE OF BREIT INTERACTION IN CALCULATION OF ENERGY
C                                                   LV, 21-12-1990
C VERSION 0.4: ADAPTED TO  CRAY
C                                                   LV, 23-1-1991
C VERSION 0.5: GIVE POPULATION ANALYSIS IN DIFFERENT FORM IF DESIRED
C                                                   OV, 7-2-1991
C VERSION 0.6: BUG FIX IN MAKONE (OPEN-CLOSED AND CLOSED-OPEN INTEGRALS)
C                                                   LV, 17-7-1991
C
C VERSION 0.7: READ ONE ELECTRON INTEGRALS AND OTHER INFORMATION FROM
C              MRCONEE                              LV, 10-9-1991
C
C VERSION 0.8: BUGFIX WRITING CI VECTORS (CASE OF 0 DETS IN REPRESENTATION)
C                                                   LV, 25-11-1991
C VERSION 1.0: WORKS WITH ABELIAN SYMMETRY
C                                                   LV, 17-02-1992
C VERSION 1.1: MINOR BUGFIXES LV 2-3-1992
C VERSION 1.2: OPTION TO SPECIFY NUMBER OF ELECTRONS IN IREF (LV 8-9-1992)
C VERSION 1.3: REPLACE ' AND " BY 1 AND 2 IN NAMES OF CIVECFILES (LV SEPT.'92)
C VERSION 1.4: ADAPTED TO HP WORKSTATIONS   (BDJ MAY 1993)
C VERSION 1.5: SX2 CALLS REMOVED. PROGRAM ONLY WORKS FOR HP, CRAY, CONVEX AND
C              CYBER  (BDJ JULY '93)
C VERSION 1.6: Cleaned up, Ported to IBM RS6000
C              LV (April 8, 1994)
C VERSION 1.7: Uses time-reversal unique canonical list of integrals
C              LV (June 28, 1994)
C VERSION 1.8: Now also real groups can be handled all real
C VERSION 2.0: Reads also from FCIDUMP files (LV, June 22, 2020)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sy_pamgos */
C    (LV, 2020) Note: this is a terrible example of code duplication, should be refactored
      SUBROUTINE SY_PAMGOS(CI,MAXCORE,NELECT,NGAS,IGAS,
     &                     MAXE,MINE,NGASO)
C***********************************************************************
C     ===============================================================
C     SMALL CI PROGRAM TO PERFORM COMPLETE OPEN SHELL CI CALCULATIONS
C     METHOD DESCRIBED IN JCP 96, 2910, (1992).
C     FIRST VERSION WRITTEN NOVEMBER 1990 BY O. VISSER
C     THEORETICAL CHEMISTRY DEPARTMENT,
C     UNIVERSITY OF GRONINGEN (RUG), GRONINGEN,THE NETHERLANDS
C     Adaption for GAS by Tommi Matila, Apr 27 2000
C     ==========================================================
!
!     Last-update: S. Yamamoto - 2007.06.18
!                  Copied and modified from PAMGOS.
!     Called by :  SY_GOSRES
!----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "goscom.h"
C
      DIMENSION CI(MAXCORE)
      DIMENSION MINE (NGAS),MAXE (NGAS),IGAS(NGAS)
      DIMENSION NGASO(*)
C
      CHARACTER*8 DATEX*10,TIMEX
C
!.s/sya,2007.01.30
#include "cosbit.h"
!.q
!
      EXTERNAL INIT
!
!.s/sya,2007.01.30
#include "cossya.h"
#include "ibtfun.h"
!.q
!----------------------------------------------------------------------
      CALL DAYTIME(DATEX,TIMEX)
      WRITE(6,1000) DATEX,TIMEX
      CALL CPUUSED(SEC)
C
C     Read namelist input
      CALL READGOS(NELECT)
!
!.s/sya,2007.02.02
      IF( KCOSINP ) THEN
         CALL RESETE(NGAS,NGASO,MAXE,MINE)
      END IF
!.q
!
!.s/sya,2007.01.30
      CALL IBTABINI1()
!.q
C
C     Set pointers to allocate blocks in memory
C
      NRB4    = NSPACT**4
      KINDTWR = 1
      KINDTWI = KINDTWR + NRB4
      IF (BREIT) THEN
         KINDGAR = KINDTWI + NRB4  
         KINDGAI = KINDGAR + NRB4  
         KTOP    = KINDGAI + NRB4
      ELSE
         KTOP    = KINDTWI + NRB4
         KINDGAR = KTOP  
         KINDGAI = KTOP
      ENDIF
C
      NDETMX = 1
C.....Create masks fo GAS spaces
      CALL GASMASK(NGAS,IGAS,NGASO)
      WRITE(6,'(A)') 'Number of determinants:'
C.....loop over fermion irreps I Abelian subgroup
      IERR = 0
      DO IRP=1, NREP
C.......odd number of electrons: fermion irreps
        IRRP = IRP
C.......even number of electrons: boson irreps
        IF (MOD(NELACT,2).EQ.0) IRRP = IRP + NREP
C.......Generate list of determinants
        CALL GENDET(IRP,IRRP,NGAS,IGAS,MINE,MAXE,IERR)
        WRITE(6,'(3X,A,A4,A,I10)') 'Symmetry ',REPNA(IRRP),'  : ',NDET
        NDETMX = MAX0(NDETMX,NDET)
      ENDDO
      IF(IERR.NE.0) CALL Q2JOB(3,'P.GENDET','    N2',NDETMX)
C
      KHR  = KTOP
      KHI  = KHR + NDETMX*NDETMX
      KTOP = KHI + NDETMX*NDETMX
      MAXINT = (MAXCORE - KTOP - 2) / 2
C
C     Check if this still fits in the available memory
C
      IF (MAXINT.LE.0)  CALL Q2JOB (3,'SY_PAMGOS','MAXINT',NRB4)
C
      KTWOINR = KTOP
      KTWOINI = KTWOINR + MAXINT
      KTOP    = KTWOINI + MAXINT 
C
      CALL RDTWO(CI(KINDTWR),CI(KINDTWI),CI(KINDGAI),CI(KINDGAR),
     &           CI(KTWOINR),CI(KTWOINI),MAXINT)
      DO 10 IRP=1,NREP
        IRRP = IRP
        IF (MOD(NELACT,2).EQ.0) IRRP = IRP + NREP
        CALL CPUUSED(SEC0)
        CALL GENDET(IRP,IRRP,NGAS,IGAS,MINE,MAXE,IERR)
        IF (NDET.EQ.0) GOTO 10
        CALL CPUUSED(SEC1)
        CALL MAKEH(CI(KHR),CI(KHI),CI(KINDTWR),CI(KINDTWI),
     &       CI(KINDGAI),CI(KINDGAI),CI(KTWOINR),CI(KTWOINI),MAXINT)
	CALL CPUUSED(SEC2)
C
        IF (IPRNT .GE. 11) THEN
          CALL HEADER('Real part of total H',-1)
          CALL OUTPUT(CI(KHR),1,NDET,1,NDET,NDET,NDET,1,6)
          CALL HEADER('Imag part of total H',-1)
          CALL OUTPUT(CI(KHI),1,NDET,1,NDET,NDET,NDET,1,6)
        END IF
C
        CALL DIAGH(IRP,CI(KHR),CI(KHI))
C
        IF (IPRNT .GE. 11) THEN
          CALL HEADER('Real part of diagonalized H',-1)
          CALL OUTPUT(CI(KHR),1,NDET,1,NDET,NDET,NDET,1,6)
          CALL HEADER('Imag part of diagonalized H',-1)
          CALL OUTPUT(CI(KHI),1,NDET,1,NDET,NDET,NDET,1,6)
        END IF
C
	CALL CPUUSED(SEC3)
!
!.s/sya,2007.02.05/
!       ...Writing CI-vector in a direct-access file.
!          Since the PRTOUT routine destroys CI(KHR) data,
!          we have to save the CI vectors in advance.
!.q
        CALL WRVEC(IRP,IRRP,CI(KHR),CI(KHI))
	CALL CPUUSED(SEC4)
        WRITE (6,1010) SEC1-SEC0,SEC2-SEC1,SEC3-SEC2,SEC4-SEC3
        CALL PRTOUT(IRP,IRRP,CI(KHR),CI(KHI))
 10     CONTINUE
      INTREP = MOD(NELACT,2)
      CALL PRTOUT(0,INTREP,CI(KHR),CI(KHI))
!
!.s/sya,2007.02.05/DV/
!
!       ...Calculation of 1-el transition density matries (including
!          usual density matrix)
!          Now, integrals are not necessary.  We abandon them.
!
      IF( SYA_KTRDMAT ) THEN
         KTOP = KTWOINR
         IPRINT = SYA_PRINT
         CALL SY_HI()
         CALL SY_TRDM(IPRINT,CI(KTOP),KTOP,MAXCORE)
         CALL SY_TMOM(IPRINT,CI(1),MAXCORE)
      END IF
!.q
      CALL CPUUSED(SEC5)
      WRITE (6,1020) SEC5-SEC
!
 1000 FORMAT(///1X,'GOSCIP (Version 1.8)'//
     &' Today is :',T15,A10/' The time is :',T17,A8//)
 1010 FORMAT(/' CPU time for this representation '//
     &' Generation of determinants :',T30,F12.4/
     &' Building the CI matrix :',T30,F12.4/
     &' Diagonalization :',T30,F12.4/
     &' Writing CI vectors :',T30,F12.4)
 1020 FORMAT (//'Total CPU time :',T30,F12.4//8X,
     &          '(NORMAL END OF PROGRAM)'///)
C
C VERSION 0.1: SEEMS TO BE WORKING FOR P-SHELLS OF ATOM
C              POPULATIONS ARE GIVEN ONLY IF $POPAN IS PRESENT
C                                                   OV, 19-11-1990
C VERSION 0.2: COMPRESSED LIST OF EIGENVALUES GIVEN AT THE END OF THE OUTPUT
C                                                   OV, 17-12-1990
C VERSION 0.3: TAKE CARE OF BREIT INTERACTION IN CALCULATION OF ENERGY
C                                                   LV, 21-12-1990
C VERSION 0.4: ADAPTED TO  CRAY
C                                                   LV, 23-1-1991
C VERSION 0.5: GIVE POPULATION ANALYSIS IN DIFFERENT FORM IF DESIRED
C                                                   OV, 7-2-1991
C VERSION 0.6: BUG FIX IN MAKONE (OPEN-CLOSED AND CLOSED-OPEN INTEGRALS)
C                                                   LV, 17-7-1991
C
C VERSION 0.7: READ ONE ELECTRON INTEGRALS AND OTHER INFORMATION FROM
C              MRCONEE                              LV, 10-9-1991
C
C VERSION 0.8: BUGFIX WRITING CI VECTORS (CASE OF 0 DETS IN REPRESENTATION)
C                                                   LV, 25-11-1991
C VERSION 1.0: WORKS WITH ABELIAN SYMMETRY
C                                                   LV, 17-02-1992
C VERSION 1.1: MINOR BUGFIXES LV 2-3-1992
C VERSION 1.2: OPTION TO SPECIFY NUMBER OF ELECTRONS IN IREF (LV 8-9-1992)
C VERSION 1.3: REPLACE ' AND " BY 1 AND 2 IN NAMES OF CIVECFILES (LV SEPT.'92)
C VERSION 1.4: ADAPTED TO HP WORKSTATIONS   (BDJ MAY 1993)
C VERSION 1.5: SX2 CALLS REMOVED. PROGRAM ONLY WORKS FOR HP, CRAY, CONVEX AND
C              CYBER  (BDJ JULY '93)
C VERSION 1.6: Cleaned up, Ported to IBM RS6000
C              LV (April 8, 1994)
C VERSION 1.7: Uses time-reversal unique canonical list of integrals
C              LV (June 28, 1994)
C VERSION 1.8: Now also real groups can be handled all real
! VERSION 1.9: Stand-alone COS_CI 
      RETURN
      END
C
C     ===========================
      SUBROUTINE DIAGH(IRP,HR,HI)
C     ===========================
      use memory_allocator
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "goscom.h"
      DIMENSION HR(NDET,NDET),HI(NDET,NDET)
      real*8, allocatable :: ee(:,:),e1(:),e2(:)
      real*8, allocatable :: e(:)
C
      call alloc(e, n2, id="e")
      call alloc(ee,n2*n2,2, id="ee")
      call alloc(e1,n2, id="e1")
      call alloc(e2,2*n2, id="e2")
      CALL CRDIAG(NDET,NDET,HR,HI,D(1,IRP),EE(1,1),EE(1,2),E,E1,E2,IERR)
      call dealloc(e2, id="e1")
      call dealloc(e1, id="e2")
      call dealloc(ee, id="ee")
      call dealloc(e, id="e")
      RETURN
      END

C     ====================================================
      SUBROUTINE GENDET(IRP,IRRP,NGAS,IGAS,MINE,MAXE,IERR)
C     ====================================================
!----------------------------------------------------------------------
!     Last-update: S. Yamamoto - 2007.01.30
!     Last-update: S. Yamamoto - 2007.09.27
!----------------------------------------------------------------------

      IMPLICIT REAL*8 (A-H, O-Z)
C
CTROND: IRP is not used !!!
CTROND: Convert ibtfun etc....
#include "goscom.h"
C
      DIMENSION IGAS (NGAS),MINE (NGAS),MAXE (NGAS)
      INTEGER N(N2)
!
!.s/sya,2007.01.30
#include "cosbit.h"
!.q
C
#include "ibtfun.h"
!----------------------------------------------------------------------
      IF (IPRNT.GE.10) WRITE(6,1000) REPNA(IRRP)
C
C     Loop over ALL possible determinants (ineffective!!)
C
      NDET = 0
      DO I = 1, 2**NSPACT
C
C       check that total number of electrons and symmetry is correct
C
!.s/sya,2007.01.30
!#      IF( (NBITS(I).EQ.NELACT).AND.(NBITS(IBTAND(I,IREF)).EQ.IREFE)
!#   +     .AND.(IRPDET(I).EQ.IRRP))  THEN
        CALL NCOSCI_NBITS(ITMP1,I)
        IF( ITMP1 .NE. NELACT ) CYCLE
cc      ITMP1 = IBTAND(I,IREF)
cc      CALL NCOSCI_NBITS(ITMP2,ITMP1)
cc      IF( ITMP2 .NE. IREFE ) CYCLE
        IF( IRPDET(I).NE. IRRP ) CYCLE
!.q
C
C       check that electron occupation in each GAS is correct
C
        DO JGAS = 1,NGAS
!.s/sya,2007.01.30
!#        IELGAS = NBITS(IBTAND(I,IGAS(JGAS)))
          ITMP1 = IBTAND(I,IGAS(JGAS))
          CALL NCOSCI_NBITS(IELGAS,ITMP1)
!.q
          IF( (IELGAS.GT.MAXE(JGAS)) .OR. (IELGAS.LT.MINE(JGAS)))
     &      GOTO 20
        END DO
        NDET = NDET + 1
        IF (NDET .GT. N2) THEN
          IERR = 1
          GOTO 20
        ENDIF
        IDET (NDET) = I
        IF(IPRNT.GE.10) CALL DETPRI(I,NDET)
!.s/sya,2007.01.30
!#    END IF
!.q
!.s/sya,2007.09.27
   20   CONTINUE
!.q
      END DO
      RETURN
 1000 FORMAT(/' List of determinants with symmetry ',A4,':'//)
      END
C
C     ===========================================================
      SUBROUTINE HMAT (ID,JD,HMR,HMI,INDTWR,INDTWI,INDGAR,INDGAI,
     &                 TWOINR,TWOINI,MAXINT)

C
CMiro Sept2016: consider c1,c2,c3 scaling parameters according  for DFT-COSCI
C  according to Roemelt et al, The Journal of Chemical Physics 138, 204101 (2013)
C
C     ===========================================================
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION TWOINR(0:MAXINT),TWOINI(0:MAXINT)
      DIMENSION INDTWR(NSPACT,NSPACT,NSPACT,NSPACT)
      DIMENSION INDTWI(NSPACT,NSPACT,NSPACT,NSPACT)
      DIMENSION INDGAR(NSPACT,NSPACT,NSPACT,NSPACT)
      DIMENSION INDGAI(NSPACT,NSPACT,NSPACT,NSPACT)
#include "goscom.h"
!
!.s/sya,2007.01.30
#include "cosbit.h"
!.q
C
#include "ibtfun.h"
C
      HMR = 0.0D0
      HMI = 0.0D0
C
      IDV = IBTXOR (ID, JD)
!.s/sya,2007.01.30
!#    IT = NBITS(IDV)
      CALL NCOSCI_NBITS(IT,IDV)
!.q
C
      IF (IT .EQ. 0) THEN
C
C     ======================
C     DETERMINANTS ARE EQUAL
C     ======================
C
      DO 30 I = 1, NSPACT
      IF (IBTAND (ID, IBTSHL (1, I - 1)) .NE. 0) THEN
        HMR = HMR + ONER (I, I)
        DO 20 J = I + 1, NSPACT
          IF (IBTAND (ID, IBTSHL (1, J - 1)) .NE. 0) THEN
C---------------------
C H = H + (I, I, J, J)
C---------------------
            INDR = INDTWR(I,I,J,J)
            INDI = INDTWI(I,I,J,J)
            SR = SIGN(1,INDR)
            SI = SIGN(1,INDI)
            INDR = ABS(INDR)
            INDI = ABS(INDI)
!Miro: c1 to Coulomb(Breit) term for DFT-COSCI
            HMR = HMR + SR * TWOINR (INDR)
            HMI = HMI + SI * TWOINI (INDI)
            IF (BREIT) THEN
               INDR = INDGAR(I,I,J,J)
               INDI = INDGAI(I,I,J,J)
               SR = SIGN(1,INDR)
               SI = SIGN(1,INDI)
               INDR = ABS(INDR)
               INDI = ABS(INDI)
               HMR = HMR + SR * TWOINR (INDR)
               HMI = HMI + SI * TWOINI (INDI)
            ENDIF
C---------------------
C H = H - (I, J, J, I)
C---------------------
            INDR = INDTWR(I,J,J,I)
            INDI = INDTWI(I,J,J,I)
            SR = SIGN(1,INDR)
            SI = SIGN(1,INDI)
            INDR = ABS(INDR)
            INDI = ABS(INDI)
!Miro: c2 to exchange(Breit) term for DFT-COSCI
            HMR = HMR - SR * TWOINR (INDR)
            HMI = HMI - SI * TWOINI (INDI)
            IF (BREIT) THEN
               INDR = INDGAR(I,J,J,I)
               INDI = INDGAI(I,J,J,I)
               SR = SIGN(1,INDR)
               SI = SIGN(1,INDI)
               INDR = ABS(INDR)
               INDI = ABS(INDI)
               HMR = HMR - SR * TWOINR (INDR)
               HMI = HMI - SI * TWOINI (INDI)
            ENDIF
          END IF
 20       CONTINUE
        END IF
 30     CONTINUE
      RETURN
      END IF
C
      IF(IT.EQ.2) THEN
C
C     ===================================
C     DETERMINANTS DIFFER BY ONE FUNCTION
C     ===================================
C
      IA = IBTAND(IDV, ID)
      IB = IBTAND(IDV, JD)
      DO 35 I = 1, NSPACT
        IF (IBTAND(IA, IBTSHL(1, I - 1)) .NE. 0) K = I
        IF (IBTAND(IB, IBTSHL(1, I - 1)) .NE. 0) L = I
35      CONTINUE
C
!Miro: c3 to damp all off-diagonal terms (one-el, COulomb, exchange) for DFT-COSCI
      HMR = HMR + ONER (K, L)
      HMI = HMI + ONEI (K, L)
C
      DO 40 I = 1, NSPACT
        IF ((I .EQ. K) .OR. (I .EQ. L)) GO TO 40
        IF (IBTAND (ID, IBTSHL (1, I - 1)) .NE. 0) THEN
C---------------------
C H = H + (I, I, K, L)
C---------------------
            INDR = INDTWR(I,I,K,L)
            INDI = INDTWI(I,I,K,L)
            SR = SIGN(1,INDR)
            SI = SIGN(1,INDI)
            INDR = ABS(INDR)
            INDI = ABS(INDI)
!Miro: c3 to damp all off-diagonal terms (one-el, Coulomb, exchange) for DFT-COSCI
            HMR = HMR + SR * TWOINR (INDR)
            HMI = HMI + SI * TWOINI (INDI)
            IF (BREIT) THEN
               INDR = INDGAR(I,I,K,L)
               INDI = INDGAI(I,I,K,L)
               SR = SIGN(1,INDR)
               SI = SIGN(1,INDI)
               INDR = ABS(INDR)
               INDI = ABS(INDI)
!Miro: c3 to damp all off-diagonal terms (one-el, Coulomb, exchange) for DFT-COSCI
               HMR = HMR + SR * TWOINR (INDR)
               HMI = HMI + SI * TWOINI (INDI)
            ENDIF
C---------------------
C H = H - (I, L, K, I)
C---------------------
            INDR = INDTWR(I,L,K,I)
            INDI = INDTWI(I,L,K,I)
            SR = SIGN(1,INDR)
            SI = SIGN(1,INDI)
            INDR = ABS(INDR)
            INDI = ABS(INDI)
!Miro: c3 to damp all off-diagonal terms (one-el, Coulomb, exchange) for DFT-COSCI
            HMR = HMR - SR * TWOINR (INDR)
            HMI = HMI - SI * TWOINI (INDI)
            IF (BREIT) THEN
               INDR = INDGAR(I,L,K,I)
               INDI = INDGAI(I,L,K,I)
               SR = SIGN(1,INDR)
               SI = SIGN(1,INDI)
               INDR = ABS(INDR)
               INDI = ABS(INDI)
!Miro: c3 to damp all off-diagonal terms (one-el, Coulomb, exchange) for DFT-COSCI
               HMR = HMR - SR * TWOINR (INDR)
               HMI = HMI - SI * TWOINI (INDI)
            ENDIF
        END IF
 40     CONTINUE
      PH = PHASE (K, ID) * PHASE (L, JD)
      HMR = PH * HMR
      HMI = PH * HMI
      RETURN
      END IF
C
      IF(IT.EQ.4) THEN
C
C     ====================================
C     DETERMINANTS DIFFER BY TWO FUNCTIONS
C     ====================================
C
      IA = IBTAND(IDV, ID)
      IB = IBTAND(IDV, JD)
      L = 0
      N = 0
      DO 50 I = 1, NSPACT
        IF (IBTAND(IA, IBTSHL(1, I - 1)) .NE. 0) THEN
          K = L
          L = I
        END IF
        IF (IBTAND(IB, IBTSHL(1, I - 1)) .NE. 0) THEN
          M = N
          N = I
        END IF
 50     CONTINUE
C
C---------------------
C H = H + (K, M, L, N)
C---------------------
      INDR = INDTWR(K,M,L,N)
      INDI = INDTWI(K,M,L,N)
      SR = SIGN(1,INDR)
      SI = SIGN(1,INDI)
      INDR = ABS(INDR)
      INDI = ABS(INDI)
!Miro: c3 to damp all off-diagonal terms (one-el, Coulomb, exchange) for DFT-COSCI
      HMR = HMR + SR * TWOINR (INDR)
      HMI = HMI + SI * TWOINI (INDI)
      IF (BREIT) THEN
         INDR = INDGAR(K,M,L,N)
         INDI = INDGAI(K,M,L,N)
         SR = SIGN(1,INDR)
         SI = SIGN(1,INDI)
         INDR = ABS(INDR)
         INDI = ABS(INDI)
!Miro: c3 to damp all off-diagonal terms (one-el, Coulomb, exchange) for DFT-COSCI
         HMR = HMR + SR * TWOINR (INDR)
         HMI = HMI + SI * TWOINI (INDI)
      ENDIF
C---------------------
C H = H - (K, N, L, M)
C---------------------
      INDR = INDTWR(K,N,L,M)
      INDI = INDTWI(K,N,L,M)
      SR = SIGN(1,INDR)
      SI = SIGN(1,INDI)
      INDR = ABS(INDR)
      INDI = ABS(INDI)
!Miro: c3 to damp all off-diagonal terms (one-el, Coulomb, exchange) for DFT-COSCI
      HMR = HMR - SR * TWOINR (INDR)
      HMI = HMI - SI * TWOINI (INDI)
      IF (BREIT) THEN
         INDR = INDGAR(K,N,L,M)
         INDI = INDGAI(K,N,L,M)
         SR = SIGN(1,INDR)
         SI = SIGN(1,INDI)
         INDR = ABS(INDR)
         INDI = ABS(INDI)
!Miro: c3 to damp all off-diagonal terms (one-el, Coulomb, exchange) for DFT-COSCI
         HMR = HMR - SR * TWOINR (INDR)
         HMI = HMI - SI * TWOINI (INDI)
      ENDIF
      PH = PHASE (K, ID) * PHASE (L, ID) *
     +     PHASE (M, JD) * PHASE (N, JD)
      HMR = HMR * PH
      HMI = HMI * PH
      RETURN
      END IF
C
C     =============================================
C     DETERMINANTS DIFFER BY MORE THEN TWO FUNCTIONS
C     ==============================================
C
      HMR = 0.0D0
      HMI = 0.0D0
      RETURN
      END
C
C     ===============
      BLOCK DATA INIT
C     ===============
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "goscom.h"
C
      DATA REPNA / N5A * '    ' /
      DATA MULTB / N5B * 0 /
      DATA ONER, ONEI / N3 * 0.0D0, N3 * 0.0D0 /
      DATA IDET / N2 * 0 /
      END
C
C     =====================
      FUNCTION IRPDET(JDET)
C     =====================
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "goscom.h"
C
      JRRP = 0
      ID=JDET
      DO 10 I=1,NSPACT
        IF (MOD(ID,2).EQ.1) THEN
          JRRP = MULTB(IRPAMO(I),JRRP)
        ENDIF
        ID=ID/2
 10     CONTINUE
      IRPDET=JRRP
      RETURN
      END
C
C     ================
      SUBROUTINE MAKEH(HR,HI,INDTWR,INDTWI,INDGAR,INDGAI,TWOINR,
     &                 TWOINI,MAXINT)
C     ================
      IMPLICIT REAL*8 (A-H, O-Z)
C
      DIMENSION HR(NDET,NDET),HI(NDET,NDET)
      DIMENSION TWOINR(0:MAXINT),TWOINI(0:MAXINT)
      DIMENSION INDTWR(NSPACT,NSPACT,NSPACT,NSPACT)
      DIMENSION INDTWI(NSPACT,NSPACT,NSPACT,NSPACT)
      DIMENSION INDGAR(NSPACT,NSPACT,NSPACT,NSPACT)
      DIMENSION INDGAI(NSPACT,NSPACT,NSPACT,NSPACT)
#include "goscom.h"
C
      DO I = 1, NDET
        DO J = 1, I
          CALL HMAT (IDET(I),IDET(J),HR(I,J),
     &               HI(I,J),INDTWR,INDTWI,
     &               INDGAR,INDGAI,TWOINR,TWOINI,MAXINT)
          HR (J,I) =  HR (I, J)
          HI (J,I) = -HI (I, J)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C     =============================
      INTEGER FUNCTION NBITS (JDET)
C     =============================
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "goscom.h"
C
      ID = JDET
      NB = 0
      DO 10 I = 1, NSPACT
        NB = NB + MOD (ID, 2)
 10     ID = ID / 2
      NBITS = NB
      RETURN
      END
C
C     =====================
      FUNCTION PHASE(I, ID)
C     =====================
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "goscom.h"
C
#include "ibtfun.h"
C
      P = 1.0D0
      DO 10 K = I + 1, NSPACT
        IF (IBTAND (ID, IBTSHL (1, K - 1)) .NE. 0) P = - P
 10     CONTINUE
      PHASE = P
      RETURN
      END
C
C     =================================
      SUBROUTINE PRTOUT(IRP,IRRP,HR,HI)
C     =================================
      use memory_allocator
      IMPLICIT REAL*8 (A-H, O-Z)
C
      CHARACTER FMT*2
      DIMENSION HR(NDET,NDET),HI(NDET,NDET)
      integer, allocatable :: IRPS(:),INDX(:),IND(:)
      real(8), allocatable :: DS(:)
C
#include "goscom.h"
C
      CHARACTER*1 N(N1)
C
#include "ibtfun.h"
 100  IF (IRP.EQ.0) GOTO 110
      NDETI(IRP)=NDET
C
C     ------------------------------------------------
C     PRINT POPULATION ANALYSIS FOR REPRESENTATION IRP
C     ------------------------------------------------
C
      DO 70 I=1,NDET
        IF (D(I,IRP)-D(1,IRP).LT.SELPOP) MX=I
 70   CONTINUE
      MX = MIN(MX,MAXPOP)
C
      WRITE(6,1006) REPNA(IRRP),MX
 1006 FORMAT(/' Population analysis for representation ',A4
     +       /'            The first',I6,' vectors are analyzed'/
     &       /T8,'energy',T22,'det#',T28,'determinant',
     &        T56,'norm',T66,'Re',T76,'Im'
     &       /T8,'------',T22,'----',T28,'-----------',
     &        T56,'----',T66,'--',T76,'--')
!.q
C
      allocate(DS(NSPACT))
      WRITE(FMT,'(I2)') N1
      DO 30 I=1,MX
        DS = 0.0D0
        WRITE(6,1000) D(I,IRP)+ECORE
        DO 40 J=1,NDET
          HNORM = HR(J,I)*HR(J,I) + HI(J,I)*HI(J,I)
          DO 50 JJ=1,NSPACT
            N(JJ)='0'
            IF (IBTAND(IDET(J),IBTSHL(1,JJ-1)).NE.0) THEN
              N(JJ)='1'
              DS(JJ)=DS(JJ)+HNORM
            ENDIF
 50       CONTINUE
          IF (HNORM.GT.THRESH) THEN
            DO 60 JJ=NSPACT+1,N1
              N(JJ)=' '
 60         CONTINUE
            WRITE(6,'(17X,I8,2X,'//FMT//'(A1),3F10.4)') 
     &           J,(N(JJ),JJ=1,N1),HNORM,HR(J,I),HI(J,I)
          ENDIF
 40     CONTINUE
        SUMOCC=0.0D0
        WRITE(6,'(A)') '    Orb.  Orbital energy       Occupation'
        DO JJ = 1,NSPACT
          WRITE(6,'(3X,I5,G20.10,5X,F8.4)') JJ,EPS(JJ),DS(JJ)
          SUMOCC=SUMOCC+DS(JJ)
        ENDDO
        WRITE(6,'(7X,A,F8.4)') 'Sum of occupations :',SUMOCC
        WRITE(6,*) 
 30   CONTINUE
      deallocate(DS)
      return
C
C     ------------------------------------
C     PRINT COMPRESSED LIST OF EIGENVALUES
C     ------------------------------------
C
 110  CONTINUE
      allocate(IRPS(N2*N5))
      allocate(INDX(N2*N5))
      allocate(IND(N5))
      allocate(DS(N2*N5))
C
      IRPS = 0
      INDX = 0
      IND  = 0
      ND = 0
      DO I = 1, NREP
         CALL DCOPY (NDETI(I),D(1,I),1,DS(ND+1),1)
         ND = ND + NDETI(I)
      ENDDO
C
      IF (IRRP.EQ.0) THEN
C        Even number of electrons, boson irreps
         CALL PRT_EV (6,NREP,NDETI,REPNA(NREP+1),ECORE,0.0D0,1,
     &                DS,N2*N5,INDX,IRPS,IND,xxdummy,-1)
      ELSE
C        Odd number of electrons, fermion irreps
         CALL PRT_EV (6,NREP,NDETI,REPNA,ECORE,0.0D0,1,
     &                DS,N2*N5,INDX,IRPS,IND,xxdummy,-1)
      ENDIF
C
 1000 FORMAT(F17.8)
      deallocate(IRPS)
      deallocate(INDX)
      deallocate(IND)
      deallocate(DS)
      RETURN
      END
C
C     ===========================================================
      SUBROUTINE RDTWO(INDTWR,INDTWI,INDGAR,INDGAI,TWOINR,TWOINI,
     &                 MAXINT)
C     ==========================================================
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "goscom.h"
C
      PARAMETER(D0 = 0.0D0)
      CHARACTER*10 DATEX,TIMEX*8
      DIMENSION INDTWR(NSPACT,NSPACT,NSPACT,NSPACT)
      DIMENSION INDTWI(NSPACT,NSPACT,NSPACT,NSPACT)
      DIMENSION INDGAR(NSPACT,NSPACT,NSPACT,NSPACT)
      DIMENSION INDGAI(NSPACT,NSPACT,NSPACT,NSPACT)
      DIMENSION TWOINR(0:MAXINT),TWOINI(0:MAXINT)
      DIMENSION KR(-N1:N1)
      COMPLEX*16 RPHASE(N1),ONE
      INTEGER INDK(N3),INDL(N3)
C
      TWOINR(0) = D0
      TWOINI(0) = D0

      INDTWR = 0
      INDTWI = 0
      IF (BREIT) INDGAR = 0
      IF (BREIT) INDGAI = 0

      if (in_fcidump /= 0) then
         call read_fcidump (2,n1,norb,nspact,irpamo,
     &                      ECORE,EPS,ONER,ONEI,
     &                      INDTWR,INDTWI,TWOINR,TWOINI,maxint)
         return
      end if
C
      MDINT = 1
      NUNIQ = 0
      OPEN(MDINT,FILE='MDCINT',FORM='UNFORMATTED')
      Read (mdint,err=10000,end=10000) datex,timex,nkr,
     &     (KR(I),KR(-I),I=1,NKR) 
C--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nkrcoul = nkr
  10  CONTINUE
      IF (REALAR) THEN
         read (mdint,ERR=10010,END=10010) ikr,jkr,nz,
     &                (indk(inz),indl(inz),inz=1,nz),
     &                (twoinr(inz),inz=nuniq+1,nuniq+nz)
         do inz = nuniq+1,nuniq+nz
            twoini(inz) = D0
         enddo
      ELSE
         read (mdint,ERR=10010,END=10010) ikr,jkr,nz,
     &                (indk(inz),indl(inz),inz=1,nz),
     &                (twoinr(inz),twoini(inz),inz=nuniq+1,nuniq+nz)
      ENDIF
      If (nuniq+nz.GT.MAXINT) Call Q2Job (3,'P.RDTWO','MAXINT',nuniq+nz)
      If (ikr.EQ.0) GoTo 20 ! End of File
C-->  Distribute this set and generate time-reversed integrals.
C-->  K|ia> = |ib> ; K|ib> = - |ia>
      I = kr(ikr)
      itr = kr(-ikr)
      j = kr(jkr)
      jtr = kr(-jkr)
      SignIJ = SIGN(1,ikr) * SIGN(1,jkr)
      do inz = 1, nz
         nuniq = nuniq + 1
         kkr = indk(inz)
         k = kr(kkr)
         ktr = kr(-kkr)
         lkr = indl(inz)
         l = kr(lkr)
         ltr = kr(-lkr)
         SignKL = SIGN(1,kkr) * SIGN(1,lkr)
C--> Original integral plus time-reversed partners
         INDTWR(I,J,K,L) = NUNIQ
         INDTWR(JTR,ITR,K,L) = NUNIQ * SignIJ
         INDTWR(I,J,LTR,KTR) = NUNIQ * SignKL
         INDTWR(JTR,ITR,LTR,KTR) = NUNIQ * SignIJ * SignKL
         INDTWI(I,J,K,L) = NUNIQ
         INDTWI(JTR,ITR,K,L) = NUNIQ * SignIJ
         INDTWI(I,J,LTR,KTR) = NUNIQ * SignKL
         INDTWI(JTR,ITR,LTR,KTR) = NUNIQ * SignIJ * SignKL
C--> Complex conjugate plus time-reversed partners
         INDTWR(J,I,L,K) = NUNIQ
         INDTWR(ITR,JTR,L,K) = NUNIQ * SignIJ
         INDTWR(J,I,KTR,LTR) = NUNIQ * SignKL
         INDTWR(ITR,JTR,KTR,LTR) = NUNIQ * SignIJ * SignKL
         INDTWI(J,I,L,K) = - NUNIQ
         INDTWI(ITR,JTR,L,K) = - NUNIQ * SignIJ
         INDTWI(J,I,KTR,LTR) = - NUNIQ * SignKL
         INDTWI(ITR,JTR,KTR,LTR) = - NUNIQ * SignIJ * SignKL
C--> Particle interchanged plus time-reversed partners
         INDTWR(K,L,I,J) = NUNIQ
         INDTWR(LTR,KTR,I,J) = NUNIQ * SignKL
         INDTWR(K,L,JTR,ITR) = NUNIQ * SignIJ
         INDTWR(LTR,KTR,JTR,ITR) = NUNIQ * SignIJ * SignKL
         INDTWI(K,L,I,J) = NUNIQ
         INDTWI(LTR,KTR,I,J) = NUNIQ * SignKL
         INDTWI(K,L,JTR,ITR) = NUNIQ * SignIJ
         INDTWI(LTR,KTR,JTR,ITR) = NUNIQ * SignIJ * SignKL
C--> Particle interchanged and complex conjugated plus time-reversed partners
         INDTWR(L,K,J,I) = NUNIQ
         INDTWR(KTR,LTR,J,I) = NUNIQ * SignKL
         INDTWR(L,K,ITR,JTR) = NUNIQ * SignIJ
         INDTWR(KTR,LTR,ITR,JTR) = NUNIQ * SignIJ * SignKL
         INDTWI(L,K,J,I) = - NUNIQ
         INDTWI(KTR,LTR,J,I) = - NUNIQ * SignKL
         INDTWI(L,K,ITR,JTR) = - NUNIQ * SignIJ
         INDTWI(KTR,LTR,ITR,JTR) = - NUNIQ * SignIJ * SignKL
      Enddo
      Goto 10
 20   Continue
      CLOSE(1)
      Write (6,1000) "Coulomb",datex,timex,nuniq
      IF (IPRNT .GT. 20) THEN
        call header('two electron integrals',-1)
        do l = 1,NSPACT
        do k = 1,NSPACT
        do j = 1,NSPACT
        do I = 1,NSPACT
           ii = indtwr(I,j,k,l)
           is = sign(1,ii)
           ii = abs(ii)
           ji = indtwi(I,j,k,l)
           js = sign(1,ji)
           ji = abs(ji)
           if (twoinr(ii) .ne. 0.0d0 .or. twoini(ji) .ne. 0.0d0)
     &    write(6,'(1P,1X,4I3,4X,2(D14.7,1X))') I,j,k,l,
     &        is*twoinr(ii),js*twoini(ji)
         end do
         end do
         end do
         end do
      END IF
C
      If (.NOT.BREIT) Return
      NUNIQC = NUNIQ
      OPEN(MDINT,FILE='MDBINT',FORM='UNFORMATTED')
      Read (mdint,err=10001,end=10001) datex,timex,nkr,
     &     (KR(I),KR(-I),I=1,NKR) 
      If (NKR.NE.NKRCOUL) Then
         Write (6,*) 'Basis set size on Coulomb file :',nkrcoul
         Write (6,*) 'Basis set size on Gaunt file :  ',nkr
         Stop 'Gaunt and Coulomb 2-el files incompatible'
      Endif
  30  CONTINUE
      IF (REALAR) THEN
         read (mdint,ERR=10010,END=10010) ikr,jkr,nz,
     &                (indk(inz),indl(inz),inz=1,nz),
     &                (twoinr(inz),inz=nuniq+1,nuniq+nz)
      ELSE
         read (mdint,ERR=10010,END=10010) ikr,jkr,nz,
     &                (indk(inz),indl(inz),inz=1,nz),
     &                (twoinr(inz),twoini(inz),inz=nuniq+1,nuniq+nz)
      ENDIF
      If (nuniq+nz.GT.N4A) Call Q2Job (3,'P.RDTWO','   N4A',nuniq+nz)
      If (ikr.EQ.0) GoTo 40 ! End of File
C-->  Distribute this set and generate time-reversed integrals.
C-->  K|ia> = |ib> ; K|ib> = - |ia>
      I = kr(ikr)
      itr = kr(-ikr)
      j = kr(jkr)
      jtr = kr(-jkr)
      SignIJ = - SIGN(1,ikr) * SIGN(1,jkr)
      do inz = 1, nz
         nuniq = nuniq + 1
         kkr = indk(inz)
         k = kr(kkr)
         ktr = kr(-kkr)
         lkr = indl(inz)
         l = kr(lkr)
         ltr = kr(-lkr)
         SignKL = - SIGN(1,kkr) * SIGN(1,lkr)
C--> Original integral plus time-reversed partners
         INDGAR(I,J,K,L) = NUNIQ
         INDGAR(JTR,ITR,K,L) = NUNIQ * SignIJ
         INDGAR(I,J,LTR,KTR) = NUNIQ * SignKL
         INDGAR(JTR,ITR,LTR,KTR) = NUNIQ * SignIJ * SignKL
         INDGAI(I,J,K,L) = NUNIQ
         INDGAI(JTR,ITR,K,L) = NUNIQ * SignIJ
         INDGAI(I,J,LTR,KTR) = NUNIQ * SignKL
         INDGAI(JTR,ITR,LTR,KTR) = NUNIQ * SignIJ * SignKL
C--> Complex conjugate plus time-reversed partners
         INDGAR(J,I,L,K) = NUNIQ
         INDGAR(ITR,JTR,L,K) = NUNIQ * SignIJ
         INDGAR(J,I,KTR,LTR) = NUNIQ * SignKL
         INDGAR(ITR,JTR,KTR,LTR) = NUNIQ * SignIJ * SignKL
         INDGAI(J,I,L,K) = - NUNIQ
         INDGAI(ITR,JTR,L,K) = - NUNIQ * SignIJ
         INDGAI(J,I,KTR,LTR) = - NUNIQ * SignKL
         INDGAI(ITR,JTR,KTR,LTR) = - NUNIQ * SignIJ * SignKL
C--> Particle interchanged plus time-reversed partners
         INDGAR(K,L,I,J) = NUNIQ
         INDGAR(LTR,KTR,I,J) = NUNIQ * SignKL
         INDGAR(K,L,JTR,ITR) = NUNIQ * SignIJ
         INDGAR(LTR,KTR,JTR,ITR) = NUNIQ * SignIJ * SignKL
         INDGAI(K,L,I,J) = NUNIQ
         INDGAI(LTR,KTR,I,J) = NUNIQ * SignKL
         INDGAI(K,L,JTR,ITR) = NUNIQ * SignIJ
         INDGAI(LTR,KTR,JTR,ITR) = NUNIQ * SignIJ * SignKL
C--> Particle interchanged and complex conjugated plus time-reversed partners
         INDGAR(L,K,J,I) = NUNIQ
         INDGAR(KTR,LTR,J,I) = NUNIQ * SignKL
         INDGAR(L,K,ITR,JTR) = NUNIQ * SignIJ
         INDGAR(KTR,LTR,ITR,JTR) = NUNIQ * SignIJ * SignKL
         INDGAI(L,K,J,I) = - NUNIQ
         INDGAI(KTR,LTR,J,I) = - NUNIQ * SignKL
         INDGAI(L,K,ITR,JTR) = - NUNIQ * SignIJ
         INDGAI(KTR,LTR,ITR,JTR) = - NUNIQ * SignIJ * SignKL
      Enddo
      Goto 30
 40   Continue
      CLOSE(1)
      Write (6,1000) "Gaunt",datex,timex,nuniq-nuniqc
      RETURN
 1000 Format (/1X,A7,' integral file was generated at ',A10,1x,A8,
     &/' Read',I8,' unique integrals')
10000 STOP 'ERROR READING HEADER OF MDCINT'
10001 STOP 'ERROR READING HEADER OF MDBINT'
10010 STOP 'ERROR READING INTEGRALS FROM MDCINT'
10011 STOP 'ERROR READING INTEGRALS FROM MDBINT'
      END
C
C     ==========================
      SUBROUTINE READGOS(NELECT)
C     ==========================
!----------------------------------------------------------------------
!     Last-update: S. Yamamoto - 2007/06/11 chukyo-u
!     Last-update: S. Yamamoto - 2007.09.14 chukyo-u
!----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "goscom.h"
!
!.s/sya,2007.01.30
#include "cosbit.h"
!.q
!
!.s/sya,2007.06.09
#include "cossya.h"
!.q
C
      CHARACTER*14 REPN(16)
C
      LOGICAL SPINFR
      INTEGER REFDET,IRPMO(N1),FCI
      REAL*8 SELPOP
      DIMENSION MINE(16),MAXE(16),IGAS(16)
      DIMENSION NGASO(256)
!.s/sya,2007.09.14
      EXTERNAL ISUM
      DATA IONE/1/
!.q
C
      NAMELIST/GOSCIP/NELACT,NGAS,MINE,MAXE,NGASO,REFDET,IREFE,IPRNT,FCI
      NAMELIST/POPAN/THRESH,DEGEN,SELPOP,MAXPOP
!----------------------------------------------------------------------
      MRCONEE=75
C
      REWIND(5)
      NELACT=0
      IPRNT=0
      REFDET=0
      IREFE=0
      MINE=0
      MAXE=0
      IGAS=0
      NGASO=0
      FCI=0
      READ(5,GOSCIP,END=200)
  200 WRITE(6,GOSCIP)
      NELACT = MAX(NELACT,NELECT)
      NELECT = NELACT
      IREF=REFDET
      IF( IREFE .LE. 0 ) THEN
        CALL NCOSCI_NBITS(IREFE,IREF)
      END IF
C
      REWIND(5)
      THRESH=1.0D-3
      DEGEN=1.0D-10
      SELPOP=1.0D2
      MAXPOP=100
      REWIND(5)
      READ(5,POPAN,END=100)
  100 WRITE(6,POPAN)
C
      IF (FCI == 0) THEN
         IN_FCIDUMP = 0
         OPEN(MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED')
         READ(MRCONEE) NSPACT,BREIT,ECORE,NFSYM,NZ,SPINFR
         REALAR = NZ .EQ. 1 .OR. SPINFR
         IF( NSPACT .GT. N1 ) CALL Q2JOB(3,'P.READIN','    N1',NSPACT)
         READ(MRCONEE) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
         READ(MRCONEE) NREP,(REPNA(IRP),IRP=1,2*NREP)
         READ(MRCONEE) ((MULTB(I,J),I=1,2*NREP),J=1,2*NREP)
         READ(MRCONEE) (IRPMO(IMO),IRPAMO(IMO),EPS(IMO),IMO=1,NSPACT)
         READ(MRCONEE)
     &    ((ONER(IMO,JMO),ONEI(IMO,JMO),IMO=1,NSPACT),JMO=1,NSPACT)
         CLOSE(MRCONEE)
      ELSE
         in_fcidump = 1
         nspact = 1 ! First read header to check symmetry and get norb
         call read_fcidump (0,n1,norb,nspact,irpamo,
     &                      ECORE,EPS,ONER,ONEI,
     &                      INDTWR,INDTWI,TWOINR,TWOINI,0)
         nrep = 1 ! interface works for C1 (NOSYM) only
         ! Could read symmetry info from FCITABLE, not implemented to
         ! keep code simple
         nfsym = 1
         repn  = ' A         '
         repna(1) = '   A'
         repna(2) = '   a'
         irpmo = irpamo
         multb(1,1) = 2
         multb(2,1) = 1
         multb(1,2) = 1
         multb(2,2) = 2
         nspact = norb ! Read again to get ecore, eps, oner and onei
         call read_fcidump (1,n1,norb,nspact,irpamo,
     &                      ECORE,EPS,ONER,ONEI,
     &                      INDTWR,INDTWI,TWOINR,TWOINI,0)
      ENDIF
!
!     Check NSPACT
!     ------------
!
      IACMO = ISUM(NFSYM,SYA_ACTIV,IONE)
      IACSP = IACMO * 2
      IF( IACSP .NE. NSPACT ) THEN
        IF( NSPACT .GT. 0 ) THEN
          IACSP = NSPACT
          IACMO = IACSP/2
        ELSE IF( IACSP .GT. 0 ) THEN
          NSPACT = IACSP
          WRITE(6,"(
     &      /1X,'* WARNING: NSPACT (number of active spinors)'
     &      /1X,'has been replaced by the cosci input data: ',I2)") 
     &      NSPACT
        ELSE
          WRITE(6,"(
     &      /1X,'# EROOR: Inconsistency between NSPACT and IACSP',
     &      /1X,2I5)") NSPACT,IACSP
          STOP
        END IF
      END IF
C
      DO 250 I = 1, 2 * NREP
        MULTB(I,0) = I
  250 CONTINUE
C
      IF (FCI == 0) THEN
         WRITE(6,1200)
      ELSE
         WRITE(6,1201)
      ENDIF

      WRITE(6,1210) NELACT,IREF
!
!     Check NELACT
!     ------------
!
      IACEL = ISUM(SYA_NOPEN,SYA_IELC,IONE)
      IF( IACEL .EQ. 0 ) THEN
        DO I = 1, SYA_NOPEN
          IACEL = MAX(IACEL,SYA_MAXE(I))
        END DO
        WRITE(6,"(/1X,'IACEL (number of active electrons) ',
     &            /1X,'has been calculated from MAXE: ',I2)") IACEL
      END IF
!
      IF( NELACT .NE. IACEL ) THEN
        IF( NELACT .GT. 0 ) THEN
          IACEL = NELACT 
        ELSE IF( IACEL .GT. 0 ) THEN
          NELACT = IACEL
          WRITE(6,"(/1X,'NELACT (number of active electrons) ',
     &              /1X,'has been replaced ',
     &                  'by the cosci input data: ',I2)") NELACT
        ELSE
          WRITE(6,"(
     &      /1X,'# EROOR: Inconsistency between NELACT and IACEL',
     &      /1X,2I5)") NELACT,IACEL
          STOP
        END IF
      END IF
!.q
      WRITE(6,1220) BREIT,ECORE,NSPACT
      DO 300 I=1,NSPACT
        WRITE(6,1230) I,REPN(IRPMO(I)),EPS(I)
 300  CONTINUE
C
      IF (IPRNT .GE. 11) THEN
        CALL HEADER('Real part of core H',-1)
        CALL OUTPUT(ONER,1,NSPACT,1,NSPACT,N1,N1,1,6)
        CALL HEADER('Imag part of core H',-1)
        CALL OUTPUT(ONEI,1,NSPACT,1,NSPACT,N1,N1,1,6)
      END IF
C
      RETURN
 1200 FORMAT(//' All information read from MRCONEE (and input)')
 1201 FORMAT(//' All information read from FCIDUMP (and input)')
 1210 FORMAT(//' Number of electrons:',I5/
     +         ' Decimal representation of frozen orbitals:',I4)
 1220 FORMAT(//' Breit interaction included :',T50,L1
     +        /' Core energy :',T30,G20.10
     +        /' Number of active orbitals :',T47,I4
     +       //' Orbital',T10,'Representation',T32,'Orbital energy')
 1230 FORMAT(I6,T10,A14,T30,G20.10)
      END
C
C     ==============================================
      SUBROUTINE SORTEW (DS,NDETT,IRPS)
C     ==============================================
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "goscom.h"
C
      REAL*8 DS(N2*N5)
      INTEGER NDETT,IND(N5),IRPS(N2*N5)
C
      DO IRP = 1, NREP
         IF (NDETI(IRP).LT.1) D(1,IRP)=1.D12
      ENDDO
C
      INDEX = 1
      IRPMIN = 1
      DO 1 IRP = 2, NREP
    1   IF (D(1,IRP).LT.D(1,IRPMIN)) IRPMIN = IRP
C
      DO 2 IRP = 1, NREP
    2   IND(IRP) = 1
C
      NDETT = 0
      DO 3 IRP = 1, NREP
        NDETT = NDETT + NDETI(IRP)
    3   CONTINUE
C
 10   DS(INDEX)=D(IND(IRPMIN),IRPMIN)
      IRPS(INDEX)=IRPMIN
      INDEX=INDEX+1
      IF (INDEX.GT.NDETT) GOTO 30
      IND(IRPMIN)=IND(IRPMIN)+1
      IF (IND(IRPMIN).GT.NDETI(IRPMIN)) D(IND(IRPMIN),IRPMIN)=1.D12
      DO 20 IRP = 1, NREP
      IF (D(IND(IRPMIN),IRPMIN).GT.D(IND(IRP),IRP)) THEN
        IRPMIN=IRP
      ENDIF
 20   CONTINUE
      GOTO 10
 30   CONTINUE
C
      RETURN
      END
C
C     =================================
      SUBROUTINE WRVEC (IRP,IRRP,HR,HI)
C     =================================
!
!----------------------------------------------------------------------
!     Writing CI vector into a direct-access file.
!     Last-update: S. Yamamoto - 2007.08.29, chukyo-u
!     Called by : PAMGOS
!                 SY_PAMGOS
!----------------------------------------------------------------------
!
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     WRITE TRIAL VECTORS FOR DIRRCI TO FILE MRCTRIV
C     RECORD 1 : RECORD LENGTH,NUMBER OF VECTORS, ENERGIES
C     RECORD 2 : BIT REPRESENTATION DETERMINANTS
C     RECORD 3 : FIRST VECTOR
C     RECORD 4 : SECOND VECTOR, ETC.
C
#include "goscom.h"
      DIMENSION HR(NDET,NDET),HI(NDET,NDET)
C
      CHARACTER*11 VECFIL
      DIMENSION DV (N2)
      DATA K8BYTE/8/
!.s/sya,2007.03.03
!#    NAMELIST/CIVEC/NVEC
!.q
!
!.s/sya,2007.08.29
      DATA NRIVEC/2/
!SK   DATA K4BYTE/4/, K1WORD/1/, K4WORD/4/
      DATA K1WORD/1/, K4WORD/4/
      INTEGER KXBYTE
!.q
      DATA LUCIVC/10/
!----------------------------------------------------------------------
C
      IF (NDET.EQ.0) RETURN
      DEGEN = 1.D-10
!.s/sya,2007.08.29
!#    NVEC = 5
!#    REWIND(5)
!#    READ (5,CIVEC,END=1,ERR=1)
!#  1 CONTINUE
!#    IF (NVEC.GT.NDET) NVEC = NDET
      NVEC = NDET
      IF (NVEC .GT. NDETMX) STOP 999
#if defined INT_STAR8
      KXBYTE = 8
#else
      KXBYTE = 4
#endif
!.q
C
C     DETERMINE RECORD LENGTH
C
!.s/sya,2007.08.29
!#    LEN = MAX0(NDET,2)
!#    LENREC = 2 * LEN * 8
      LENW   = NDETMX
      LENR0  = LENW * NRIVEC * K8BYTE 
      LENR2  = (K1WORD + NDETMX) * KXBYTE
      LENR1  = K4WORD * KXBYTE + NDETMX * K8BYTE
      LENREC = MAX(LENR0,LENR2,LENR1)
!.q
      VECFIL = 'MDTRIV_'//REPNA(IRRP)
      IF (VECFIL(8:8).EQ.'''') VECFIL(8:8) = '1'
      IF (VECFIL(8:8).EQ.'"') VECFIL(8:8) = '2'
!.s/sya,2007.06.07
      IF( IPRNT .GE. 10 ) THEN
        WRITE(6,"(1X,'>>>> Writing CI-vector ',A,' RECL=',I6)")
     &    VECFIL,LENREC
        WRITE(6,*) 'IRRP=',IRRP
      END IF
!.q
      OPEN(LUCIVC,FILE=VECFIL,ACCESS='DIRECT',RECL=LENREC)
      WRITE(LUCIVC,REC=2) NDET,(IDET(I),I=1,NDET)
      IF( IPRNT .GE. 20 ) THEN
        DO I = 1, NDET
          WRITE(6,"(I5,2X,B32)") I, IDET(I)
        END DO 
      END IF
!
      IVEC = 0
      EPREV = D(1,IRP) - DEGEN - 1.D0
      DO 10 I = 1, NDET
!.s/sya,2007.02.27,del
!#    IF (D(I,IRP)-EPREV.GT.DEGEN) THEN
!.q
         IVEC = IVEC + 1
         DV (IVEC) = D(I,IRP) + ECORE
         WRITE(LUCIVC,REC=IVEC+2) (HR(J,I),HI(J,I),J=1,NDET)
         IF( IPRNT .GE. 20 ) THEN
           WRITE(6,"(' IVEC,I=',2I5)") IVEC,I
           WRITE(6,"(5F15.5)") (HR(J,I),J=1,NDET)
         END IF
         EPREV = D(I,IRP)
!.s/sya,2007.02.27,del
!#    ENDIF
!#    IF (IVEC.EQ.NVEC) GOTO 100
!.q
   10 CONTINUE
  100 CONTINUE
!.s/sya,2007.02.27,del
!#    WRITE (10,REC=1) LENREC,NVEC,(DV(IVEC),IVEC=1,NVEC)
      WRITE(LUCIVC,REC=1) LENREC,NVEC,(DV(IVEC),IVEC=1,NVEC),NDETMX,NDET
!.q
      CLOSE(LUCIVC)
      WRITE(6,1000) REPNA(IRRP),NVEC,VECFIL,(DV(IVEC),IVEC=1,NVEC)
      RETURN
!.s/sya,2007.09.29/
!#1000 FORMAT (/' Representation ',A4,'; ',I2,
!#    &' vectors written to file ',A12//' Eigenvalues',(T20,F20.12))
 1000 FORMAT (/' Representation ',A4,';',I5,
     &' vectors written to file ',A12//' Eigenvalues',(T22,F20.12))
!.q
      END
C     ===================================
      SUBROUTINE GASMASK(NGAS,IGAS,NGASO)
C     ===================================
c     Create 'masks', 
C     that is a bitstring with ones for each orbital present in the GAS
C
C     Written by Tommi Matila, Apr 27 2000
C     Polish Trond Saue
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "goscom.h"
C
      DIMENSION IGAS(NGAS),NGASO(NREP,NGAS)
      DIMENSION INDX(N1)
C
      CALL IZERO(IGAS,NGAS)
C.....Construct ordering index
C.....This one is not particularly elegant, but does the job...
      II = 0
      DO ISYM = 1,NREP
        DO I = 1,NSPACT
          IF(IRPAMO(I).EQ.ISYM) THEN
            II = II + 1
            INDX(II) = I
          ENDIF
        ENDDO
      ENDDO
C
      II = 0
C.....Loop over fermion irreps of Abelian subgroup
      DO ISYM = 1,NREP
C.......Loop over GAS spaces
        DO JGAS = 1,NGAS
C.........Loop over orbitals in each irrep/GAS
          DO IEL = 1,NGASO(ISYM,JGAS)
            II = II + 1
            IORB = INDX(II)
            IF(IRPAMO(IORB).EQ.ISYM) THEN
              IGAS(JGAS) = IGAS(JGAS) + 2**(IORB-1)
            ELSE
              CALL QUIT('GASMASK: Error in NGASO')
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck detpri*/
      SUBROUTINE DETPRI(IIDET,IIREF)
C***********************************************************************
C
C    Print a determinant on bitstring form
C    A cut and paste job by Trond Saue
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H, O-Z)
#include "goscom.h"
      INTEGER N(N2)
#include "ibtfun.h"
      DO I=1,NSPACT
        N(I)=0
        IF (IBTAND (IIDET,IBTSHL (1,I-1)).NE.0)  N(I)=1
      ENDDO
      WRITE(6,1010) IIREF,(N(I),I=1,NSPACT)
      RETURN
 1010 FORMAT(I5,32I2)
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck ncosci_nbits */
      SUBROUTINE NCOSCI_NBITS(ICOUNT,JDET)
!***********************************************************************
!
!     Count number if bits in integer JDET.
!
!     Written by J. Thyssen - Jan 2 2001
!
!     Copied from NGOSCI_NBITS
!     "IBTAB1" should be set by subroutine IBTABINI1 beforehand.
!     Modified by sya, 2007.01.31
!     Last revision : sya, Fri Mar  2 16:06:00 JST 2007
!
!***********************************************************************
#include "cosbit.h"
      INTEGER IDD(4)
      INTRINSIC IBITS
!----------------------------------------------------------------------
      IDD(1) = IBITS(JDET, 0,8)
      IDD(2) = IBITS(JDET, 8,8)
      IDD(3) = IBITS(JDET,16,8)
      IDD(4) = IBITS(JDET,24,8)
!
      ICOUNT = IBTAB1(IDD(1)) + IBTAB1(IDD(2))
     &       + IBTAB1(IDD(3)) + IBTAB1(IDD(4))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ibtabini1 */
      SUBROUTINE IBTABINI1()
C***********************************************************************
C
C     Initialize table with IBTAB(I) = NBITS(I)
C
C     Input : none
C
C     Output:
C       IBTAB
C
C     Written by J. Thyssen - Feb 27 2001
C     Last revision :
!
!     Copied from IBTABINI.
!     Modified by sya, 2007.01.31.
C
C***********************************************************************
#include "implicit.h"
!
!.s/sya,2007.01.31
#include "cosbit.h"
!.q
C
#include "ibtfun.h"
C
!.s/sya,2007.01.31
!#    DIMENSION IBTAB(0:255)
!.q
!----------------------------------------------------------------------
C
      DO I = 0, 255
         ID = I
         NB = 0
         DO J = 1, 8
            NB = NB + IBTAND(1,ID)
            ID = ID / 2
         END DO
!.s/sya,2007.01.31
!#       IBTAB(I) = NB
         IBTAB1(I) = NB
!.q
      END DO
C
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck resete */
      SUBROUTINE RESETE(NGAS,NGASO,MAXE,MINE)
!***********************************************************************
!
!     Setting MAXE and MINE from COSCI input data.
!     Written by: S. Yamamoto - Feb 02, 2007, ulp
!     Called by:  SY_PAMGOS
!
!***********************************************************************
#include "implicit.h"
#include "cossya.h"
!
!     ...for NREP
#include "goscom.h"
!
      INTEGER NGAS
      INTEGER NGASO(NREP,NGAS)
      INTEGER MAXE(NGAS)
      INTEGER MINE(NGAS)
!----------------------------------------------------------------------
      IF( NGAS .NE. SYA_NOPEN ) THEN
        WRITE(6,*) 'NGAS .NE. SYA_NOPEN',NGAS,SYA_NOPEN
        NGAS = SYA_NOPEN
        WRITE(6,"(/1X,'#WARNING: NGAS has been replaced ',
     &            /1X,'by the cosci input data: ',I2)") 
     &                 NGAS
!
        IF (NREP .NE. SYA_NREP) THEN
           WRITE(6,*) 'NREP from CIGOSC .ne. .NREP from *COSCI input',
     &     NREP, SYA_NREP
           CALL QUIT('NREP from CIGOSC .ne. .NREP from *COSCI input')
        END IF
        DO I = 1, SYA_NOPEN
          DO J = 1, NREP
            NGASO(J,I) = SYA_GASO(J,I)
          END DO
        END DO
        WRITE(6,"(/1X,'#WARNING: NGASO has been replaced ',
     &            /1X,'by the cosci input data: ')") 
      END IF
!
!     ...In this case, the same number to be set for MAXE and MINE ?
!
      DO I = 1, SYA_NOPEN
         MAXE(I) = SYA_MAXE(I)
         MINE(I) = SYA_MAXE(I)
      END DO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck sy_hi */
      SUBROUTINE SY_HI()
!***********************************************************************
!     Starting calculation of transtion density matrix (TRDM) and 
!     transition dipole moments (TMOM).
!     Printing welcome message.
!
!     Written by:  S. Yamamoto - Feb 05, 2007, ulp
!     Last-update: S. Yamamoto - 2007.09.16, chukyo-u
!     Called by:   SY_PAMGOS
!***********************************************************************
#include "implicit.h"
!     ...for LUPRI
#include "priunit.h"
!----------------------------------------------------------------------
 1000 FORMAT(
     &  1X,'Transition density matrix and transition dipole moments',
     &      ' for COSCI.',
     &  /1X,'(version 1.00)',
     &  /1X,'Programmed by Shigeyoshi YAMAMOTO (2007.09.16)' )
!----------------------------------------------------------------------
      WRITE(LUPRI,'(/A)') REPEAT('=',70)
      WRITE(LUPRI,1000)
      WRITE(LUPRI,'(A/)') REPEAT('=',70)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck sy_trdm */
      SUBROUTINE SY_TRDM(IPRINT,WW,KTOP,MAXCORE)
!***********************************************************************
!     Calculation of transition density matrix (including usual density
!     matrix).
!
!     Written by:  S. Yamamoto - Feb 05, 2007, ulp
!     Last-update: S. Yamamoto - June 4, 2007, chukyo-u
!     Called by:   SY_PAMGOS
!***********************************************************************
#include "implicit.h"
#include "goscom.h"
#include "cossya.h"
!
!     ...for NZ
#include "dgroup.h"
!
!     ...for LUPRI
#include "priunit.h"
      DIMENSION WW(*)
!
      INTEGER KSAVE,KTRDM,LTRDM,MAXINT
!
      PARAMETER (NMPART = 20)
      CHARACTER*4 CC(NMPART)
      INTEGER LL(NMPART),KK(NMPART)
      DATA NRIVEC/2/
!----------------------------------------------------------------------
      IF( IPRINT .GE. 1 ) THEN
        WRITE(LUPRI,'(/A)') REPEAT('=',70)
        WRITE(LUPRI,'(A)')
     &   ' SY_TRDM: transition density matrix calculation '
        WRITE(LUPRI,'(A/)') REPEAT('=',70)
      END IF
!
      KSAVE = KTOP
      KFREE = 1
      LFREE = MAXCORE - KTOP 
!
!     ...setting common MSTATE in cossya.h
!
      MSTATE = 0
      DO I = 1, NREP
        MSTATE = MAX(MSTATE,SYA_NTRDM(I))
      END DO
      LDET2  = NDETMX * NRIVEC * MSTATE
!.s/sya,2007.0911
!#    LTRDM  = NSPACT * NSPACT * NZ * MSTATE * MSTATE
      LTRDM  = NSPACT * NSPACT * NRIVEC * MSTATE * MSTATE
!.q
!
!     >>>> CALL MEMGET by starting from KSAVE=1
!
      LL( 1) = LDET2
      CC( 1) = 'REAL'
      LL( 2) = LDET2
      CC( 2) = 'REAL'
      LL( 3) = NDETMX
      CC( 3) = 'REAL'
      LL( 4) = NDETMX
      CC( 4) = 'REAL'
      LL( 5) = LTRDM
      CC( 5) = 'REAL'
      CALL MEMPAT('CALTRDM',IPRINT,LUPRI,NMPART,1,5,
     &  CC,KK,LL,WW,KFREE,LFREE)
!
      CALL CALTRDM(IPRINT,NRIVEC,
     &  WW(KK( 1)), WW(KK( 2)), WW(KK( 3)), WW(KK( 4)), WW(KK( 5))
     &  )
!       CL,         CR,         DET_I,      DET_J,      TRDM
!
      KREL = KK(1)
      KFIRST = KREL
      CALL MEMREL('CALTRDM',WW,KFIRST,KREL,KFREE,LFREE)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck caltrdm */
      SUBROUTINE CALTRDM(IPRINT,NRIVEC,
     &  CL,         CR,         DET_I,      DET_J,      TRDM
     &  )
!***********************************************************************
!     Calculation of transition density matrix.
!
!     Written by:  S. Yamamoto - Feb 06, 2007, ulp
!     Last-update: S. Yamamoto - 2007.05.14, chukyo-u
!     Called by:   SY_TRDM
!----------------------------------------------------------------------
!     We keep TRDM on the memory only for a combination of two irreps,
!     not for all combinations of irreps.  After generating TRDM, we
!     write it on a file.  Then we advance to another combination of
!     irreps.
!***********************************************************************
#include "implicit.h"
#include "goscom.h"
!
!     ...for MSTATE
#include "cossya.h"
#include "cosbit.h"
#include "dgroup.h"
#include "priunit.h"
      DIMENSION CL(NDETMX,NRIVEC,MSTATE)
      DIMENSION CR(NDETMX,NRIVEC,MSTATE)
      INTEGER DET_I(NDETMX),DET_J(NDETMX)
!.s/sya,2007.0911
!#    DIMENSION TRDM(NSPACT,NSPACT,NZ,MSTATE,MSTATE)
      DIMENSION TRDM(NSPACT,NSPACT,NRIVEC,MSTATE,MSTATE)
!.q
!
      LOGICAL SPINFR
      LOGICAL DOIMAG
      DATA LUTRDM/11/  ! logical unit number for writing TRDM
!----------------------------------------------------------------------
      IRVEC = 0
      JRVEC = 0
      DOIMAG = NZ .GT. 1
      CALL IBTABINI1()
!.s/sya,2007.0911
!#    LX = NSPACT*NSPACT*NZ*MSTATE*MSTATE
      LX = NSPACT*NSPACT*NRIVEC*MSTATE*MSTATE
!.q
!
      CALL OPTRDM(IPRINT,LUTRDM,NREP,NSPACT,NREC)
!
      DO IRRP_I = 1, NREP
        IREP_I = IRRP_I
        IF( MOD(NELACT,2) .EQ. 0 ) IREP_I = IRRP_I + NREP
        NVEC_I = SYA_NTRDM(IRRP_I)
        IF( NVEC_I .EQ. 0 ) THEN
          WRITE(LUPRI,"(/'<bra|      irrep',I3,A,
     &                   '     out of scope.')") 
     &          IRRP_I,REPNA(IREP_I)
          CYCLE
        END IF
        CALL RDCIVEC(NRIVEC,MSTATE,IRRP_I,IREP_I,NVEC_I,CL,DET_I,
     &               NZCONF_I,KST,IRVEC,IPRINT)
        IF( KST .NE. 0 ) THEN
          WRITE(LUPRI,"(/'<bra|      irrep',I3,A,
     &                   '     skipped.',I4)") 
     &          IREP_I,REPNA(IREP_I),KST
          CYCLE
        END IF
!
        DO IRRP_J = 1, IRRP_I
          IREP_J = IRRP_J
          IF( MOD(NELACT,2) .EQ. 0 ) IREP_J = IRRP_J + NREP
          IF( IRRP_I .EQ. IRRP_J ) THEN
            NVEC_J = NVEC_I
          ELSE
            NVEC_J = SYA_NTRDM(IRRP_J)
          END IF
!
          IF( NVEC_J .EQ. 0 ) THEN
            WRITE(LUPRI,"(/'<bra|ket>  irrep',I3,A,I3,A,
     &                     '  out of scope.')") 
     &            IRRP_I,REPNA(IREP_I),IRRP_J,REPNA(IREP_J)
            CYCLE
          END IF
!
          CALL RDCIVEC(NRIVEC,MSTATE,IRRP_J,IREP_J,NVEC_J,CR,DET_J,
     &                 NZCONF_J,KST,JRVEC,IPRINT)
!
          IF( KST .NE. 0 ) THEN
            WRITE(LUPRI,"(/'<bra|ket>  irrep',I3,A,I3,A,
     &                     '  skipped.',I4)") 
     &        IRRP_I,REPNA(IREP_I),IRRP_J,REPNA(IREP_J),KST
            CYCLE
          END IF
!
          CALL DZERO(TRDM,LX)
!             
!         Loop over determinants "ID"
!         ===========================
!
          DO ICONF_I = 1, NZCONF_I
            ID = DET_I(ICONF_I)
            CALL COSCI_GENTRDM(NRIVEC,ID,DET_J,DOIMAG,ICONF_I,CL,CR,
     &           IRRP_I,IRRP_J,NVEC_I,NVEC_J,NZCONF_J,TRDM)
          END DO
!
          IF( IPRINT .GE. 2 ) THEN
            CALL OUTTRDM(NRIVEC,IRRP_I,IRRP_J,IREP_I,IREP_J,
     &                   NVEC_I,NVEC_J,TRDM)
          END IF
!.s/sya,2007.0911
!#        CALL WRTRDM(LUTRDM, NSPACT, MSTATE, TRDM, 
!#   &                IRRP_I, IRRP_J, NVEC_I, NVEC_J, NREC)
          CALL WRTRDM(NRIVEC,LUTRDM, NSPACT, MSTATE, TRDM, 
     &                IRRP_I, IRRP_J, NVEC_I, NVEC_J, NREC)
!.q
        END DO
      END DO
!
      CALL CSTRDM(LUTRDM,NREC)
!
      IF( IRVEC .EQ. 0 .AND. JRVEC .EQ. 0 ) THEN
        WRITE(LUPRI,"(/1X,'No CI-vectors read')")
        STOP
      END iF
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck rdcivec */
      SUBROUTINE RDCIVEC(NRIVEC,MSTATE,IRRP,IREP,NVEC,C,IXDET,
     &  NZCONF,KST,IREAD,IPRINT)
!**********************************************************************
!     Read CI vector from a direct-access file.
!
!     Written by:  S. Yamamoto - 2007.02.06, ulp
!     Last-update: S. Yamamoto - 2007.03.03, chukyo-univ.
!     Last-update: S. Yamamoto - 2007.09.07, chukyo-univ.
!     Called by:   CLTRDM
!----------------------------------------------------------------------
!     NRIVEC is 2.
!**********************************************************************
#include "implicit.h"
#include "goscom.h"
!
!     ...for NZ
#include "dgroup.h"
      DIMENSION C(NDETMX,NRIVEC,MSTATE)
      DIMENSION IXDET(NDETMX)
      CHARACTER*11 VECFIL
      DATA K8BYTE/8/
      DATA K4BYTE/4/, K1WORD/1/, K4WORD/4/ ! Miro: necessary to fix runtimecheck
      INTEGER KXBYTE
      DATA LUCIVC/10/
!----------------------------------------------------------------------
      KST = 1
#if defined INT_STAR8
      KXBYTE = 8
#else
      KXBYTE = 4
#endif
!.s/sya,2007.08.29
!#    LEN = MAX0(NDET,2)
!#    LENREC = LEN * 2 * 8
      LENW   = NDETMX
      LENR0  = LENW * NRIVEC * K8BYTE 
      LENR2  = (K1WORD + NDETMX) * KXBYTE
      LENR1  = K4WORD * KXBYTE + NDETMX * K8BYTE
      LENREC = MAX(LENR0,LENR2,LENR1)
!.q
      VECFIL = 'MDTRIV_'//REPNA(IREP)
      IF( VECFIL(8:8).EQ.'''' ) VECFIL(8:8) = '1'
      IF( VECFIL(8:8).EQ.'"' ) VECFIL(8:8) = '2'
      IF( IPRINT .GE. 10 ) THEN
        WRITE(6,"(1X,'>>>> Reading CI-vector ',A,' RECL=',I6)")
     &    VECFIL,LENREC
      END IF
      OPEN(LUCIVC,FILE=VECFIL,ACCESS='DIRECT',RECL=LENREC,
     &     STATUS='OLD',ERR=900,IOSTAT=IST)
      IF (IST.NE.0) THEN
        CALL QUIT('Error in opening file with lu=LUCIVC !')
      ENDIF
! miro: with the g95/ILP64/ownlibs/opt there is a reading hangup ! Please fix it !
      READ(LUCIVC,REC=2,ERR=901,IOSTAT=IST) 
     &  NZCONF,(IXDET(I),I=1,NZCONF)
      IF (IST.NE.0) THEN
! miro: TODO - fix ERR GOTO 901 
        CALL QUIT('Error in reading file with lu=LUCIVC !')
      ENDIF
      IF( IPRINT .GE. 20 ) THEN
        DO I = 1, NZCONF
          WRITE(6,"(I5,2X,B32)") I, IXDET(I)
        END DO
      END IF
!
      DO IVEC = 1, NVEC
        KST = IVEC
        READ(LUCIVC,REC=IVEC+2,ERR=902,IOSTAT=IST) 
     &       (( C(I,J,IVEC),J=1,NRIVEC),I=1,NZCONF)
        IF( IPRINT .GE. 20 ) THEN
          WRITE(6,"(' IVEC,NVEC=',2I5)") IVEC,NVEC
          DO I = 1, NZCONF
            WRITE(6,"(1X,2F20.10)") ( C(I,J,IVEC),J=1,NRIVEC)
          END DO
        END IF
        IREAD = IREAD + 1
      END DO
!
      CLOSE(LUCIVC)
      KST = 0
      RETURN
!
  900 CONTINUE
      KST = 900
      RETURN
  901 CONTINUE
      KST = 901
      RETURN
  902 CONTINUE
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck cosci_gentrdm */
      SUBROUTINE COSCI_GENTRDM(NRIVEC,ID,DET_J,DOIMAG,ICONF_I,CL,CR,
     &           IRRP_I,IRRP_J,NVEC_I,NVEC_J,NZCONF_J,TRDM)
!***********************************************************************
!     Generating 1-el transition density matrix.
!
!     Written by:  S. Yamamoto - 2007.02.06, ulp
!     Last-update: S. Yamamoto - 2007.06.10, chukyo-u
!     Called by:   CALTRDM
!---------------------------------------------------------------------- 
!     Determinant !DET_I> is now an integer variable named "ID".
!***********************************************************************
#include "implicit.h"
#include "goscom.h"
#include "cossya.h"
#include "cosbit.h"
#include "ibtfun.h"
!
!     ...for NZ
#include "dgroup.h"
      DIMENSION CL(NDETMX,NRIVEC,MSTATE),CR(NDETMX,NRIVEC,MSTATE)
!.s/sya,2007.0911
!#    DIMENSION TRDM(NSPACT,NSPACT,NZ,MSTATE,MSTATE)
      DIMENSION TRDM(NSPACT,NSPACT,NRIVEC,MSTATE,MSTATE)
!.q
      LOGICAL DOIMAG
      INTEGER DET_J(NDETMX)
      DATA ITMAX/2/  ! only 1-el density matrix.
!----------------------------------------------------------------------
!             
!     Loop over determinants "JD"
!     ---------------------------
!
      DO ICONF_J = 1, NZCONF_J
!
        JD = DET_J(ICONF_J)
        IDV = IBTXOR(ID,JD)
        CALL NCOSCI_NBITS(IT,IDV)  ! number of different bits
        IF( IT .GT. ITMAX ) CYCLE
        IF( IT .EQ. 0 ) THEN
!      
!         Same determinants 
!         -----------------
!
          CALL COSCI_DMAT0(NRIVEC,ID,
     &         ICONF_I,CL,CR,NVEC_I,NVEC_J,
     &         IRRP_I,IRRP_J,
     &         DOIMAG,MSTATE,TRDM)
!
        ELSE
!      
!         Determinants differ by one function
!         -----------------------------------
!
          CALL COSCI_DMAT2(NRIVEC,ID,JD,IDV,
     &         ICONF_I,ICONF_J,CL,CR,NVEC_I,NVEC_J,
     &         IRRP_I,IRRP_J,
     &         DOIMAG,MSTATE,TRDM)
        END IF
      END DO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck cosci_dmat0 */
      SUBROUTINE COSCI_DMAT0(NRIVEC,ID,
     &           ICONF_I,CL,CR,NVEC_I,NVEC_J,
     &           IRRP_I,IRRP_J,
     &           DOIMAG,MSTATE,TRDM)
!***********************************************************************
!     Calculate contributions determinant ID and JD (ID = JD).
!
!     Output:
!        TRDM    : density matrices
!
!     Written by:  S. Yamamoto - 2007.02.23, ulp
!     Last-update: S. Yamamoto - 2007.03.03, chukyo-u
!     Called by:   COSCI_GENTRDM
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "goscom.h"
!
!     ...for NZ
#include "dgroup.h"
      DIMENSION CL(NDETMX,NRIVEC,MSTATE),CR(NDETMX,NRIVEC,MSTATE)
!.s/sya,2007.0911
!#    DIMENSION TRDM(NSPACT,NSPACT,NZ,MSTATE,MSTATE)
      DIMENSION TRDM(NSPACT,NSPACT,NRIVEC,MSTATE,MSTATE)
!.q
      LOGICAL DOIMAG
!
#include "thrzer.h"
#include "ibtfun.h"
!----------------------------------------------------------------------
      DO I = 1, NSPACT
        IF( IBTAND(ID, IBTSHL (1, I - 1)) .EQ. 0) CYCLE
!
        DO IVEC_I = 1, NVEC_I
!
          JVEND = NVEC_J
          IF( IRRP_I .EQ. IRRP_J ) JVEND = IVEC_I
!
          DO IVEC_J = 1, JVEND
            IF (.NOT. DOIMAG) THEN
              DMR = CL(ICONF_I,1,IVEC_I) * CR(ICONF_I,1,IVEC_J)
            ELSE
              DMR = CL(ICONF_I,1,IVEC_I) * CR(ICONF_I,1,IVEC_J)
     &            + CL(ICONF_I,2,IVEC_I) * CR(ICONF_I,2,IVEC_J)
              DMI = CL(ICONF_I,1,IVEC_I) * CR(ICONF_I,2,IVEC_J) 
     &            - CL(ICONF_I,2,IVEC_I) * CR(ICONF_I,1,IVEC_J)
              IF( ABS(DMI) .GT. THRZER ) THEN
                TRDM(I,I,2,IVEC_I,IVEC_J) = 
     &            TRDM(I,I,2,IVEC_I,IVEC_J) + DMI
              END IF
            END IF
!
            TRDM(I,I,1,IVEC_I,IVEC_J) = 
     &        TRDM(I,I,1,IVEC_I,IVEC_J) + DMR
          END DO
        END DO
      END DO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck cosci_dmat2 */
      SUBROUTINE COSCI_DMAT2(NRIVEC,ID,JD,IDV,
     &           ICONF_I,ICONF_J,CL,CR,NVEC_I,NVEC_J,
     &           IRRP_I,IRRP_J,
     &           DOIMAG,MSTATE,TRDM)
!***********************************************************************
!     Calculate contributions determinant ID and JD (ID != JD).
!     Determinants differ by one function K -> M
!
!     Output:
!        TRDM    : density matrices
!
!     Written by:  S. Yamamoto - Fri Feb 23 13:21:14 CST 2007, ulp
!     Last-update: S. Yamamoto - 2007.06.10, chukyo-u
!     Called by:   COSCI_GENTRDM
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "goscom.h"
!
!     ...for NZ
#include "dgroup.h"
      DIMENSION CL(NDETMX,NRIVEC,MSTATE),CR(NDETMX,NRIVEC,MSTATE)
!.s/sya,2007.0911
!#    DIMENSION TRDM(NSPACT,NSPACT,NZ,MSTATE,MSTATE)
      DIMENSION TRDM(NSPACT,NSPACT,NRIVEC,MSTATE,MSTATE)
!.q
      LOGICAL DOIMAG
!
#include "thrzer.h"
#include "ibtfun.h"
!----------------------------------------------------------------------
      IA = IBTAND(IDV, ID)
      IB = IBTAND(IDV, JD)
      DO I = 1, NSPACT
        IF( IBTAND(IA, IBTSHL(1, I - 1)) .NE. 0) K = I
        IF( IBTAND(IB, IBTSHL(1, I - 1)) .NE. 0) M = I
      END DO
      PH = PHASE(K, ID) * PHASE(M, JD)
!
      DO IVEC_I = 1, NVEC_I
!
        JVEND = NVEC_J
        IF( IRRP_I .EQ. IRRP_J ) JVEND = IVEC_I
!
        DO IVEC_J = 1, JVEND
          DMRIJ = CL(ICONF_I,1,IVEC_I) * CR(ICONF_J,1,IVEC_J)
          IF( DOIMAG ) THEN
            DMRIJ = DMRIJ 
     &            + CL(ICONF_I,2,IVEC_I) * CR(ICONF_J,2,IVEC_J)
            DMIIJ = CL(ICONF_I,1,IVEC_I) * CR(ICONF_J,2,IVEC_J) 
     &            - CL(ICONF_I,2,IVEC_I) * CR(ICONF_J,1,IVEC_J)
            IF( ABS(DMIIJ) .GT. THRZER ) THEN
              TRDM(K,M,2,IVEC_I,IVEC_J) = TRDM(K,M,2,IVEC_I,IVEC_J)
     &                                  + PH * DMIIJ
            END IF
          END IF
!
          TRDM(K,M,1,IVEC_I,IVEC_J) = TRDM(K,M,1,IVEC_I,IVEC_J)
     &                              + PH * DMRIJ
        END DO
      END DO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck outtrdm */
!.s/sya,2007.0911
!#    SUBROUTINE OUTTRDM(IRRP_I,IRRP_J,IREP_I,IREP_J,
!#   &  NVEC_I,NVEC_J,TRDM)
      SUBROUTINE OUTTRDM(NRIVEC,IRRP_I,IRRP_J,IREP_I,IREP_J,
     &  NVEC_I,NVEC_J,TRDM)
!.q
!***********************************************************************
!     Printing 1-el transition density matrix.
!
!     Written by: S. Yamamoto - Sat Feb 24 15:40:31 CST 2007, ulp
!     Called by:  CALTRDM
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "goscom.h"
!
!     ...for MSTATE
#include "cossya.h"
#include "dgroup.h"
!.s/sya,2007.0911
!#    DIMENSION TRDM(NSPACT,NSPACT,NZ,MSTATE,MSTATE)
      DIMENSION TRDM(NSPACT,NSPACT,NRIVEC,MSTATE,MSTATE)
!.q
!
      DATA KREAL/1/, KIMAG/2/
!----------------------------------------------------------------------
      DO IVEC_I = 1, NVEC_I
!
        JVEND = NVEC_J
        IF( IRRP_I .EQ. IRRP_J ) JVEND = IVEC_I
!
        DO IVEC_J = 1, JVEND
          WRITE(LUPRI,'(/A,A,I3,1X,A,A,I4, A,I3,1X,A,A,I4,A)')
     &      'TRDM',
     &      ' <irrep',IRRP_I,REPNA(IREP_I),' vec#',IVEC_I,
     &      ' | irrep',IRRP_J,REPNA(IREP_J),' vec#',IVEC_J,'>'
!
          SP = 0.0D0
          DO K = 1, NSPACT
            SP = SP + TRDM(K,K,KREAL,IVEC_I,IVEC_J)
          END DO
          WRITE(LUPRI,"(/1X,'real trace(TRDM)',G20.8)") SP
!
          WRITE(LUPRI,"(/1X,'(SY_TRDM)  real(TRDM)')")
          CALL OUTPUT(TRDM(1,1,KREAL,IVEC_I,IVEC_J),
     &      1,NSPACT,1,NSPACT,
     &      NSPACT,NSPACT,1,LUPRI)
!
          IF( NZ .GE. KIMAG ) THEN
            SP = 0.0D0
            DO K = 1, NSPACT
              SP = SP + TRDM(K,K,KIMAG,IVEC_I,IVEC_J)
            END DO
            WRITE(LUPRI,"(/1X,'imag trace(TRDM)',G20.8)") SP
!
            WRITE(LUPRI,"(/1X,'(SY_TRDM) imag(TRDM)')")
            CALL OUTPUT(TRDM(1,1,KIMAG,IVEC_I,IVEC_J),
     &        1,NSPACT,1,NSPACT,
     &        NSPACT,NSPACT,1,LUPRI)
          END IF
        END DO
      END DO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck diagtrdm */
!.s/sya,2007.0911
!#    SUBROUTINE DIAGTRDM(IRRP_I,IRRP_J,NVEC_I,TRDM)
      SUBROUTINE DIAGTRDM(NRIVEC,IRRP_I,IRRP_J,NVEC_I,TRDM)
!.q
!***********************************************************************
!     Diagonalize 1-el density matrix and obtain occupation numbers.
!
!     Written by: S. Yamamoto - Thu Mar  1 13:05:09 CST 2007, ulp
!     Called by:  CALTRDM
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "goscom.h"
#include "cossya.h"
#include "dgroup.h"
!
!.s/sya,2007.0911
!#    DIMENSION TRDM(NSPACT,NSPACT,NZ,MSTATE,MSTATE)
      DIMENSION TRDM(NSPACT,NSPACT,NRIVEC,MSTATE,MSTATE)
!.q
      DIMENSION WAR(NSPACT,NSPACT),WAI(NSPACT,NSPACT)
      DIMENSION WW(NSPACT),WZR(NSPACT,NSPACT),WZI(NSPACT,NSPACT)
      DIMENSION WFV1(NSPACT,NSPACT)
      DIMENSION WFM1(NSPACT,NSPACT,2),WFV2(NSPACT,NSPACT)
!----------------------------------------------------------------------
      IF( IRRP_I .NE. IRRP_J ) RETURN
!.s/sya,2007.0911
      LX = NSPACT * NSPACT
      CALL DZERO(WAR,LX)
      CALL DZERO(WAI,LX)
!.q
!
      DO IVEC_I = 1, NVEC_I
        DO I = 1, NSPACT
          DO J = 1, NSPACT
            WAR(J,I) = TRDM(J,I,1,IVEC_I,IVEC_I)
            WAI(J,I) = TRDM(J,I,2,IVEC_I,IVEC_I)
          END DO
        END DO
        CALL CRDIAG(NSPACT,NSPACT,WAR,WAI,WW,
     &              WZR,WZI,WFV1,WFV2,WFM1,IERR)
        WRITE(*,"(/A,I3,A,I4)") 'occupation number  irrep',IRRP_I,
     &                          ' vec#',IVEC_I
        IF( IERR .NE. 0 ) THEN
          WRITE(*,*) IERR
          CALL QUIT('DIAGTRDM: Error in CRDIAG')
        END IF
        WRITE(*,"(5F14.6)") (WW(I),I=1,NSPACT)
      END DO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck optrdm */
      SUBROUTINE OPTRDM(IPRINT,LUTRDM,NREP,NSPACT,NREC)
!***********************************************************************
!     Open a direct-access file for saving transition density matrix
!     and write a dummy record number (1).
!
!     Written by:  S. Yamamoto - 2007.04.11, chukyo-u
!     Last-update: S. Yamamoto - Wed Apr 11 13:42:11 JST 2007
!     Last-update: S. Yamamoto - 2007.09.17, chukyo-u
!     Called by:   CALTRDM
!***********************************************************************
#include "implicit.h"
!
!     ...for MSTATE
#include "cossya.h"
!
!     ...for NZ
#include "dgroup.h"
!     CHARACTER*6 TRDFIL/'TRDMSQ'/
      CHARACTER*6 TRDFIL          
      DATA KTDINF/5/
      DATA KRCINF/4/
      DATA K8BYTE/8/
      DATA KREAL/1/, KIMAG/2/
!----------------------------------------------------------------------
      IF( NZ .GT. KREAL ) THEN
        KRANDI = KIMAG
      ELSE
        KRANDI = KREAL
      END IF
      LENW  = KRCINF + NSPACT * NSPACT * KRANDI
      LENREC = LENW * K8BYTE
      TRDFIL = 'TRDMSQ'
!     WRITE(6,*) 'OPTRDM: TRDFIL',TRDFIL
      OPEN(LUTRDM,FILE=TRDFIL,ACCESS='DIRECT',RECL=LENREC,
     &     STATUS='REPLACE',ERR=900,IOSTAT=IST)
!
!     ...estimate the file size by assuming integer as 8-byte
!
      IF( IPRINT .GE. 1 ) THEN
        NFSIZE = KTDINF * K8BYTE 
     &         + (NREP*(NREP+1))/2 * MSTATE * MSTATE * LENW * K8BYTE
        BYTEM = REAL(NFSIZE) / 1.D6
        WRITE(6,"(/,1X,A,F15.3)") 
     &        'Estimated file size (MB) for TRDM:',BYTEM 
      END IF
!
      NREC = 1
      WRITE(LUTRDM,REC=1,ERR=901,IOSTAT=IST) NREP,NSPACT,NZ,MSTATE,NREC
      RETURN
!
  900 CONTINUE
      WRITE(6,*) "OPTRDM: open error",IST
      STOP
!
  901 CONTINUE
      WRITE(6,*) "OPTRDM: writing error",IST
      STOP
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck cstrdm */
      SUBROUTINE CSTRDM(LUTRDM,NREC)
!***********************************************************************
!     Write NREC (total number of record) in the first record and 
!     close the file.
!
!     Written by:  S. Yamamoto - 2007.04.11, chukyo-u
!     Last-update: S. Yamamoto - Wed Apr 11 13:42:11 JST 2007, chukyo-u
!     Called by:   CALTRDM
!***********************************************************************
#include "implicit.h"
!----------------------------------------------------------------------
      READ(LUTRDM,REC=1)  NREP,NSPACT,KNZ,MSTATE,IREC
      WRITE(LUTRDM,REC=1) NREP,NSPACT,KNZ,MSTATE,NREC
      CLOSE(LUTRDM)
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck wrtrdm */
!.s/sya,2007.0911
!#    SUBROUTINE WRTRDM(LUTRDM, NSPACT, MSTATE, TRDM, 
!#   &                  IRRP_I, IRRP_J, NVEC_I, NVEC_J, NREC)
      SUBROUTINE WRTRDM(NRIVEC,LUTRDM, NSPACT, MSTATE, TRDM, 
     &                  IRRP_I, IRRP_J, NVEC_I, NVEC_J, NREC)
!.q
!***********************************************************************
!     Write transition density matrix in a direct-access file.
!
!     Written by:  S. Yamamoto - 2007.04.11, chukyo-u
!     Last-update: S. Yamamoto - May 17, 2007
!     Called by:   CALTRDM
!***********************************************************************
#include "implicit.h"
!
!     ...for NZ
#include "dgroup.h"
!.s/sya,2007.0911
!#    DIMENSION TRDM(NSPACT,NSPACT,NZ,MSTATE,MSTATE)
      DIMENSION TRDM(NSPACT,NSPACT,NRIVEC,MSTATE,MSTATE)
!.q
      DATA KREAL/1/, KIMAG/2/
!----------------------------------------------------------------------
!
!     ...write real part
!
      DO IVEC_I = 1, NVEC_I
!
        JVEND = NVEC_J
        IF( IRRP_I .EQ. IRRP_J ) JVEND = IVEC_I
!
        DO IVEC_J = 1, JVEND
          NREC = NREC + 1
          WRITE(LUTRDM,REC=NREC)
     &      IRRP_I, IRRP_J, IVEC_I, IVEC_J,
     &      ((TRDM(I,J,KREAL,IVEC_I,IVEC_J),I=1,NSPACT),J=1,NSPACT)
!
!     ...write imaginary part if exits
!
          IF( NZ .GT. KREAL ) THEN
            NREC = NREC + 1
            WRITE(LUTRDM,REC=NREC)
     &      IRRP_I, IRRP_J, IVEC_I, IVEC_J,
     &      ((TRDM(I,J,KIMAG,IVEC_I,IVEC_J),I=1,NSPACT),J=1,NSPACT)
          END IF
        END DO
      END DO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck op2trdm */
      SUBROUTINE OP2TRDM(LUTRDM,NREP,NSPACT,NZ,
     &           NREC)
!***********************************************************************
!     Open existing direct-access file to read transition density 
!     matrix.
!
!     Called by: CALTMOM1
!     Written by: S. Yamamoto - 2007.04.27, chukyo-u
!***********************************************************************
#include "implicit.h"
!
!     ...for MSTATE
#include "cossya.h"
!     CHARACTER*6 TRDFIL/'TRDMSQ'/
      CHARACTER*6 TRDFIL
      DATA KTDINF/5/
      DATA KRCINF/4/
      DATA K8BYTE/8/
      DATA KREAL/1/, KIMAG/2/
!----------------------------------------------------------------------
!.s/sya,2007.0911
!#    LENW  = KRCINF + NSPACT * NSPACT * NZ
      IF( NZ .GT. KREAL ) THEN
        KRANDI = KIMAG
      ELSE
        KRANDI = KREAL
      END IF
      LENW  = KRCINF + NSPACT * NSPACT * KRANDI
!.q
      LENREC = LENW * K8BYTE
      TRDFIL = 'TRDMSQ'
!     WRITE(6,*) 'OP2TRDM: TRDFIL',TRDFIL
      OPEN(LUTRDM,FILE=TRDFIL,ACCESS='DIRECT',RECL=LENREC,
     &     STATUS='OLD',ERR=900,IOSTAT=IST)
!
      READ(LUTRDM,REC=1,ERR=901,IOSTAT=IST) KREP,KORB,KNZ,KSTATE,NREC
!
      KOK = 0
      IF( KREP   .NE. NREP   ) KOK = KOK + 1
      IF( KORB   .NE. NSPACT ) KOK = KOK + 10
      IF( KNZ    .NE. NZ     ) KOK = KOK + 100
      IF( KSTATE .NE. MSTATE ) KOK = KOK + 1000
      IF( KOK .NE. 0 ) THEN
        WRITE(*,*) 'inconsistent data in Subroutine OP2TRDM', KOK
        WRITE(*,*) NREP,NSPACT,NZ,MSTATE
        WRITE(*,*) KREP,KORB,KNZ,KSTATE,NREC
        STOP
      END IF
!
      RETURN
!
  900 CONTINUE
      WRITE(6,*) "OP2TRDM: open error",IST
      STOP
!
  901 CONTINUE
      WRITE(6,*) "OP2TRDM: writing error",IST
      STOP
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck sy_tmom */
      SUBROUTINE SY_TMOM(IPRINT,WW,MAXCORE)
!***********************************************************************
!     Calculation of transition dipole moment (including usual dipole 
!     moment).
!
!     Written by:  S. Yamamoto - 2007.04.13, chukyo-u
!     Last-update: S. Yamamoto - 2007.06.08, chukyo-u
!     Last-update: S. Yamamoto - 2007.08.29, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.13, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.23, chukyo-u
!     Called by :  SY_PAMGOS
!-----------------------------------------------------------------------
!  Rough sketch of the procedure  (sya, 2007.09.13)
!
!     [Stage 0]
!     ...Making input data containing ".TRAPRP" for "MOLTRA".
!     ...Calling traprp and writing Kramers-pair based property 
!        molecular integrals (an array PRCMO) in a sequential file 
!        (LUMLF1).
!        The integrals are indexed by spinor number including the 
!        inactive (core) part. 
!
!     ---------
!     [Stage 1]
!     ...Reading PRCMO (Kramers-pair based dipole length integral 
!        matrix including core part) and save the x, y, and z 
!        components in DIPSP.
!
!     Loop over TRDM records
!
!         [Stage 2] 
!           ...Reading active-part Kramers-pair based transition 
!              density matrix TRDM from the direct-access file.
!
!         [Stage 3]
!           ...Obtaining TODM by inserting the core part into TRDM.
!              The core part is diagonal unities for the case of
!              permanent dipole moments, and zero's for the the case of 
!              transition dipole moments.
!
!         [Stage 4]
!           ...Obtaining expectation values of permanent and transition
!              dipole length operators by calculating the trace of 
!              an inner product (TODM . DIPSP).
!              This is a complex arithmetic.
!
!         [Stage 5]
!           ...Calculating total dipole moments and their norm.
!
!     End of LOOP 
!-----------------------------------------------------------------------
!  N.B.  Control constants and variables (sya, 2007.05.06)
!     NRIVEC   : 2 (real and imaginary components)
!     N3DIM    : 3 (x, y, and z components)
!     NACTMO   : number of active MOs
!     NSPACT   : number of active spinors = NACTMO*2
!     NOCCMO   : number of (inactive + active) MOs
!     NOCCSP   : number of (inactive + active) spinors = NOCCMO*2
!     NBAS     : total number of basis functions 
!                (large and small components)
!                = NTBAS(0)
!     NNBBASX  : NBAS * (NBAS+1)/2
!     N2BBASX  : NBAS * NBAS
!     N2BBASXQ : NBAS * NBAS * NZ
!     NORBT    : total number of MOs
!     N2ORBX   : NORBT * NORBT
!     N2ORBXQ  : N2ORBX * NZ
!
!     ICMOQ(ifrp) : offset for symmetry-blockd MO coefficient matrix
!                   CMO(NBAS,NORBT,NZ).
!                   Defined in subroutine SETDC2 (in dirac/dirset.F)
!              
!                NCMOT  = 0
!                DO IFRP = 1,NFSYM
!                  ICMOQ(IFRP) = NCMOT * NZ
!                  NCMO(IFRP)  = NFBAS(IFRP,0) * NORB(IFRP)
!                  NCMOT       = NCMOT + NCMO(IFRP)
!                END DO
!----------------------------------------------------------------------
!     Arrays (sya, 2007.09.14)
!       Name   Type Dimension/Meaning
!       TMOM   r8   (NRIVEC,N3DIM)
!                   Total transition dipole moments.
!                   Electronic and nuclear contributions to total 
!                   dipole moment.
!       TDLEN  r8   (NRIVEC,N3DIM)
!                   Expectation values of dipole length operators. 
!       TDIPE  r8   (NRIVEC,N3DIM)
!                   Electronic contribution to total dipole moments.
!       DIPSP  r8   (NOCCSP,NOCCSP,NRIVEC,N3DIM)
!                   Kramers-pair based dipole moment integral matrix
!                   of inactive (core) + active parts.
!                   Generated by TRAPRP routine.
!                   Generally complex matrix.
!       TODM   r8   (NOCCSP,NOCCSP,NRIVEC)
!                   Transition density matrix of core+active part.
!                   Generally complex matrix.
!       IA2CA  i4   (NSPACT)
!                   Index mapping table which converts active spinor
!                   index to inactive+active spinor index.
!       IC2CA  i4   (NCORS)
!                   Index mapping table which converts inactive
!                   (core) spinor index to inactive+active spinor 
!                   index.
!       TRDM   r8   (NSPACT,NSPACT,NRIVEC)
!                   Transition density matrix of active part only.
!                   Generated by SY_TRDM routine.
!                   Generally complex matrix.
!       IOCCMO i4   (NFSYM)
!                   Total number of inactive+active MOs belonging to 
!                   IFRP (fermion irrep).
!       IOFFSP i4   (NFSYM)
!                   Offset for inactive+active spinor numbering.
!                   (starting spinor number) -1 for ifrp.
!                   Kramers' pair-based spinors are indexed by fermion
!                   irrep firstr.
!
!***********************************************************************
#include "implicit.h"
#include "goscom.h"
#include "cossya.h"
#include "dgroup.h"
!
!     ...for NORB
c#include "dcborb.h" 
!
!     ...for NTBAS(0)
c#include "dcbbas.h"
!
!     ...for LUPRI
#include "priunit.h"
      DIMENSION WW(MAXCORE)
      PARAMETER (NMPART = 20)
      CHARACTER*4 CC(NMPART)
      INTEGER LL(NMPART),KK(NMPART)
!
      EXTERNAL ISUM
!!    INTRINSIC REPEAT
!
      DATA NRIVEC/2/  ! real and imaginary
      DATA N3DIM/3/   ! x, y, and z
      DATA IONE/1/
!----------------------------------------------------------------------
      LWORK = MAXCORE
#include "memint.h"
!
      IF( IPRINT .GE. 1 ) THEN
        WRITE(LUPRI,'(/A)') REPEAT('=',70)
        WRITE(LUPRI,'(A)')
     &   ' SY_TMOM: transition dipole moment calculation'
        WRITE(LUPRI,'(A/)') REPEAT('=',70)
      END IF
!
!     CALL QENTER('SY_TMOM')
!
!     ...set size control data
!
      NCORMO = ISUM(NFSYM,SYA_INACT,IONE)
      NCORSP = NCORMO * 2
!
      IACTMO = NSPACT / 2
      NACTMO = ISUM(NFSYM,SYA_ACTIV,IONE)
      IF( NACTMO .NE. IACTMO ) THEN
        NSPACT = NACTMO * 2
        IACTMO = NSPACT / 2
        WRITE(LUPRI,"(/1X,'## WARNING ##',
     &        /1X,'Inconsistency in number of active orbitals',
     &        3I5)") NACTMO,IACTMO,NSPACT
      END IF
      NOCCMO = NCORMO + NACTMO
      NOCCSP = NOCCMO * 2
!
      IF( IPRINT .GE. 1 ) THEN
        WRITE(LUPRI,"(1X,'Control constants and variables ',
     &                   'in subroutine SY_TMOM.')")
        WRITE(LUPRI,"(3X,'NFSYM :   ',I5)") NFSYM
        WRITE(LUPRI,"(3X,'NREP  :   ',I5)") NREP
        WRITE(LUPRI,"(3X,'NCORMO:   ',I5)") NCORMO
        WRITE(LUPRI,"(3X,'NCORSP:   ',I5)") NCORSP
        WRITE(LUPRI,"(3X,'NACTMO:   ',I5)") NACTMO
        WRITE(LUPRI,"(3X,'NSPACT:   ',I5)") NSPACT
        WRITE(LUPRI,"(3X,'NOCCMO:   ',I5)") NOCCMO
        WRITE(LUPRI,"(3X,'NOCCSP:   ',I5)") NOCCSP
        WRITE(LUPRI,"(3X,'SYA_INACT:',2I5)") (SYA_INACT(I),I=1,NFSYM)
        WRITE(LUPRI,"(3X,'SYA_ACTIV:',2I5)") (SYA_ACTIV(I),I=1,NFSYM)
        WRITE(LUPRI,"(1X)") 
      END IF
!
!     =================================================
!     Stage 1 (reading PMCMO and keeping them in DIPSP)
!     =================================================
!
!     ...TMOM(NRIVEC,N3DIM) 
      LL( 1) = NRIVEC * N3DIM
      CC( 1) = 'REAL'
!
!     ...TDLEN(NRIVEC,N3DIM) 
      LL( 2) = NRIVEC * N3DIM
      CC( 2) = 'REAL'
!
!     ...TDIPE(NRIVEC,N3DIM)
      LL( 3) = NRIVEC * N3DIM
      CC( 3) = 'REAL'
!
!     ...DIPSP(NOCCSP,NOCCSP,NRIVEC,N3DIM) 
      LL( 4) = NOCCSP * NOCCSP * NRIVEC * N3DIM 
      CC( 4) = 'REAL'
!
!     ...TODM(NOCCSP,NOCCSP,NRIVEC)
      LL( 5) = NOCCSP * NOCCSP * NRIVEC
      CC( 5) = 'REAL'
!
!     ...IA2CA(NSPACT)
      LL( 6) = NSPACT 
      CC( 6) = 'INTE'
!
!     ...IC2CA(NSPACT)
      LL( 7) = NCORSP
      CC( 7) = 'INTE'
!
!     ...TRDM(NSPACT,NSPACT,NRIVEC)
      LL( 8) = NSPACT * NSPACT * NRIVEC
      CC( 8) = 'REAL'
!
!     ...IOCCMO(NFSYM)
      LL( 9) = NFSYM
      CC( 9) = 'INTE'
!
!     ...IOFFSP(NFSYM)
      LL(10) = NFSYM
      CC(10) = 'INTE'
!
      CALL MEMPAT('SY_TMOM',IPRINT,LUPRI,NMPART,1,10,
     &  CC,KK,LL,WW,KFREE,LFREE)
!
      CALL CLTMOM1(
     &  KFREE,      LFREE,      IPRINT,     NRIVEC,     N3DIM,
     &              NOCCMO,     NOCCSP,     
     &                                      WW(KK( 4)), 
     &                                      WW(KK( 9)), WW(KK(10)), 
     &  WW(KK(11))
     &  )
!    &  TMOM,       TDLEN,      TDIPE,      DIPSP,      TODM,
!    &  IA2CA,      IC2CA,      TRDM,       IOCCMO,     IOFFSP,      
!    &  WORK
!
!
!     ======================
!     Stage 2 (reading TRDM)
!     ======================
!
      CALL CLTMOM2(
     &  KFREE,      LFREE,      IPRINT,     NRIVEC,     N3DIM,
     &  NACTMO,     NOCCMO,     NOCCSP,     NCORMO,     NCORSP,
     &  WW(KK( 1)), WW(KK( 2)), WW(KK( 3)), WW(KK( 4)), WW(KK( 5)),
     &  WW(KK( 6)), WW(KK( 7)), WW(KK( 8)),
     &  WW(KK(11))
     &  )
!    &  TMOM,       TDLEN,      TDIPE,      DIPSP,      TODM,
!    &  IA2CA,      IC2CA,      TRDM,                           
!    &  WORK
!
      KREL = KK(1)
      KFIRST = KREL
      CALL MEMREL('SY_TMOM',WW,KFIRST,KREL,KFREE,LFREE)
!     CALL QEXIT('SY_TMOM')
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck cltmom1 */
      SUBROUTINE CLTMOM1(
     &  KFREE,      LFREE,      IPRINT,     NRIVEC,     N3DIM,
     &              NOCCMO,     NOCCSP,     
     &                                      DIPSP, 
     &                                      IOCCMO,     IOFFSP,     
     &  WORK
     &  )
!----------------------------------------------------------------------
!     Reading PRCMO (Kramers-pair based dipole length integral 
!     matrix including core part) and save the x, y, and z.
!
!     Written by:  S. Yamamoto - 2007.09.13, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.13, chukyo-u
!     Called by :  SY_TMOM
!----------------------------------------------------------------------
!  N.B.  LUMLF1 (sya, 2007.09.13)
!     LUMLF1 = 75  (defined in dcbtra.h in moltra)
!----------------------------------------------------------------------
#include "implicit.h"
!
!     ...for NSPACT
#include "goscom.h"
!
!     ...for NFSYM
#include "dgroup.h"
#include "priunit.h"
!
      DIMENSION TMOM(NRIVEC,N3DIM)
      DIMENSION TDLEN(NRIVEC,N3DIM)
      DIMENSION TDIPE(NRIVEC,N3DIM)
      DIMENSION DIPSP(NOCCSP,NOCCSP,NRIVEC,N3DIM)
      DIMENSION TODM(NOCCSP,NOCCSP,NRIVEC)
      DIMENSION IA2CA(NSPACT)
      DIMENSION TRDM(NSPACT,NSPACT,NRIVEC)
      DIMENSION IOCCMO(NFSYM)
      DIMENSION IOFFSP(NFSYM)
      DIMENSION WORK(*)
!
      CHARACTER*32 ACHAR
      DATA LUMLF1/75/
      DATA KREAL/1/, KIMAG/2/
      DATA THR/1.0D-10/
!----------------------------------------------------------------------
      IF( IPRINT .GE. 5 ) WRITE(LUPRI,"(' >>>> CLTMOM1')")
      L1 = NOCCSP * NOCCSP * NRIVEC * N3DIM
      CALL DZERO(DIPSP,L1)
!
      OPEN(LUMLF1,FILE='MDPROP',STATUS='OLD',
     &     FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
!
      DO IXYZ = 1, N3DIM
        READ(LUMLF1,ERR=901) ACHAR
        CALL CKACHR(LUPRI,ACHAR(25:32),IXYZ)
        READ(LUMLF1,ERR=902) 
     &    ((DIPSP(I,J,KREAL,IXYZ),DIPSP(I,J,KIMAG,IXYZ),
     &      I=1,NOCCSP),J=1,NOCCSP)
      END DO
!
      READ(LUMLF1,ERR=901) ACHAR
      IF( ACHAR(25:32) .NE. 'EOFLABEL' ) THEN
        WRITE(LUPRI,"(' ## NOT EOF in CLTMOM1')")
        STOP
      END IF
      CLOSE(LUMLF1,STATUS='KEEP')
!
!     ...check symmetry of DIPSP
!
      CALL CKDIPSP(NOCCMO,NOCCSP,NRIVEC,N3DIM,LUPRI,DIPSP,
     &             IOCCMO,IOFFSP)
!
!     ...dump DIPSP
!
      IF( IPRINT .GE. 10 ) THEN
        WRITE(LUPRI,
     &    "(/1X,'Kramers-pair based dipole length integral, DIPSP')")
        DO IXYZ = 1, N3DIM
          WRITE(LUPRI,"(1X,'IXYZ:',I2)") IXYZ
          DO I = 1, NOCCSP
            DO J = 1, NOCCSP
              WRITE(LUPRI,"(2I5,F18.10,2X,F18.10)")
     &          I,J,(DIPSP(I,J,K,IXYZ),K=1,NRIVEC)
            END DO
          END DO
          WRITE(LUPRI,*)
        END DO
      END IF
!
      IF( IPRINT .GE. 5 ) WRITE(LUPRI,"(' <<<< CLTMOM1')")
      RETURN
!
  901 CONTINUE
      WRITE(LUPRI,"(/1X,'ERROR in reading LUMLF1',
     &              /1X,'Invalid ACHAR',I5)") IXYZ
      STOP 
  902 CONTINUE
      WRITE(LUPRI,"(/1X,'ERROR in reading LUMLF1',
     &              /1X,'Invalid integrals',I5)") IXYZ
      STOP 
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck ckdipsp */
      SUBROUTINE CKDIPSP(NOCCMO,NOCCSP,NRIVEC,N3DIM,LUPRI,DIPSP,
     &                   IOCCMO,IOFFSP)
!----------------------------------------------------------------------
!     Checking time-reversal symmetry of the DIPSP matrix.
!
!     Written by:  S. Yamamoto - 2007.09.16, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.16, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.22, chukyo-u
!     Called by :  CLTMOM1
!----------------------------------------------------------------------
!  N.B.  Structure of DIPSP (sya, 2007.09.16)
!     Kramers-pair based matrix of dipole length operators.
!     In the case of NFSYM is 1.
!
!               p-unbar          p-bar (Kramers' counterpart of p)
!               A-string         B-string
!              *----------------+----------------*
!     A-string | Op,q           | Op,q-bar       | 
!              | A=(a + bi)     | B=(c + di)     | 
!              +----------------+----------------+
!     B-string | Op-bar,q       | Op-bar,q-bar   | 
!              | -B*=(-c + di)  | A*=(a - bi)    |  
!              *----------------+----------------*
!
!     Op-bar,q     = -(Op,q-bar)*
!     Op-bar,q-bar = (Op,q)*
!     -----------------------------------------------------------------
!     
!     In the case of NFSYM (total number of fermion irreps) is 2.
!
!              IFSYM=1       IFSYM=2
!             <-----------> <----------->
!             p      p-bar  p'     p'-bar
!            *------+------+------+------*
!     p      |  A1  |  B1  |  A2  |  B2  |
!            +---- -+------+------+------+
!     p-bar  | -B1* |  A1* | -B2* |  A2* |
!            +---- -+------+------+------+
!     p'     |  A3  |  B3  |  A4  |  B4  |
!            +------+------+------+------+
!     p'-bar | -B3* |  A3* | -B4* |  A4* |
!            *------+------+------+------*
!
!----------------------------------------------------------------------
!   N.B.  (sya, 2007.09.23)
!     IOCCMO(ifrp) : Number of inactive+active MOs belonging to fermion
!                    irrep (ifrp).
!     IOFFSP(ifrp) : Offset for ifrp.  
!                    (starting spinor number for ifrp) - 1.
!                    The last spinor number for the (ifrp-1)-th 
!                    fermion irrep.
!
!     Example
!       IFSYM=1                        IFSYM=2        
!       -----------------------------  ------------------------
!       p-unbar        p-bar           p'-unbar     p'-bar
!       core  active   core  active    core active  core active
!       01 02 03 04 05 06 07 08 09 10  11   12 13   14   15 16  (CA)
!              |  |  |        |  |  |        |  |         |  | 
!              |  |  |        |  |  |        |  |         |  | 
!              1  2  3        4  5  6        7  8         9 10  (A)
!
!       IOCCMO(1)=5
!       IOCCMO(2)=3
!       IOFFSP(1)=0
!       IOFFSP(2)=10
!----------------------------------------------------------------------
#include "implicit.h"
!
!     ...for NFSYM
#include "dgroup.h"
!
!     ...for SYA_INACT,SYA_ACTIV
#include "cossya.h"
      DIMENSION DIPSP(NOCCSP,NOCCSP,NRIVEC,N3DIM)
      DIMENSION IOCCMO(NFSYM)
      DIMENSION IOFFSP(NFSYM)
      DATA KREAL/1/, KIMAG/2/
      DATA THR/1.0D-10/
      DATA NKPAIR/2/
!----------------------------------------------------------------------
 1000 FORMAT(1X,'Op,q vs Op-bar,q-bar (real):     ',5I5,2F20.10)
 1010 FORMAT(1X,'Op,q vs Op-bar,q-bar (imaginary):',5I5,2F20.10)
 2000 FORMAT(1X,'Op,q-bar vs Op-bar,q (real):     ',5I5,2F20.10)
 2010 FORMAT(1X,'Op,q-bar vs Op-bar,q (imaginary):',5I5,2F20.10)
 3000 FORMAT(1X,'Op,q vs Oq,p (Hermiteness):      ',4I5,2F20.10)
!----------------------------------------------------------------------
      DO IFRP = 1, NFSYM
        IOCCMO(IFRP) = SYA_INACT(IFRP) + SYA_ACTIV(IFRP)
      END DO
!
      IOFFSP(1) = 0
      DO IFRP = 2, NFSYM
        IOFFSP(IFRP) = IOFFSP(IFRP-1) + IOCCMO(IFRP-1) * NKPAIR
      END DO
!
      KERR = 0
!
!     ...Checking time-reversal symmetry
!
      DO IXYZ = 1, N3DIM
        DO IFRP = 1, NFSYM
          DO IPM = 1, IOCCMO(IFRP)
            IPU = IPM + IOFFSP(IFRP)
            IPB = IPU + IOCCMO(IFRP)
            DO JFRP = 1, NFSYM
              DO IQM = 1, IOCCMO(JFRP)
                IQU = IQM + IOFFSP(JFRP)
                IQB = IQU + IOCCMO(JFRP)
!
!               ...Op,q vs Op-bar,q-bar
                DIF = DIPSP(IPU,IQU,KREAL,IXYZ)
     &              - DIPSP(IPB,IQB,KREAL,IXYZ)
                IF( ABS(DIF) .GT. THR ) THEN
                  KERR = KERR + 1
                  WRITE(LUPRI,1000) IXYZ,IPU,IQU,IPB,IQB,
     &              DIPSP(IPU,IQU,KREAL,IXYZ),
     &              DIPSP(IQU,IPU,KREAL,IXYZ)
                END IF
!
                DIF = DIPSP(IPU,IQU,KIMAG,IXYZ)
     &              + DIPSP(IPB,IQB,KIMAG,IXYZ)
                IF( ABS(DIF) .GT. THR ) THEN
                  KERR = KERR + 1
                  WRITE(LUPRI,1010) IXYZ,IPU,IQU,IPB,IQB,
     &              DIPSP(IPU,IQU,KIMAG,IXYZ),
     &              DIPSP(IQU,IPU,KIMAG,IXYZ)
                END IF
!
!               ...Op,q-bar vs Op-bar,q
                DIF = DIPSP(IPU,IQB,KREAL,IXYZ)
     &              + DIPSP(IPB,IQU,KREAL,IXYZ)
                IF( ABS(DIF) .GT. THR ) THEN
                  KERR = KERR + 1
                   WRITE(LUPRI,2000) IXYZ,IPU,IQB,IPB,IQU,
     &               DIPSP(IPU,IQB,KREAL,IXYZ),
     &               DIPSP(IPB,IQU,KREAL,IXYZ)
                END IF
!
                DIF = DIPSP(IPU,IQB,KIMAG,IXYZ)
     &              - DIPSP(IPB,IQU,KIMAG,IXYZ)
                IF( ABS(DIF) .GT. THR ) THEN
                  KERR = KERR + 1
                  WRITE(LUPRI,2010) IXYZ,IPU,IQB,IPB,IQU,
     &              DIPSP(IPU,IQB,KIMAG,IXYZ),
     &              DIPSP(IPB,IQU,KIMAG,IXYZ)
                END IF
!
              END DO 
            END DO 
          END DO 
        END DO
      END DO
!
!     ...Checking symmetry for Op-unbar,q-unbar
!        within A1 and A4 blocks
!
      DO IXYZ = 1, N3DIM
        DO IFRP = 1, NFSYM
          DO IPM = 1, IOCCMO(IFRP)
            IPU = IPM + IOFFSP(IFRP)
            DO IQM = 1, IPM
              IQU = IQM + IOFFSP(IFRP)
              DO K = 1, NRIVEC
!
                IF( K .EQ. KREAL ) THEN
                  DIF = DIPSP(IPU,IQU,K,IXYZ) 
     &                - DIPSP(IQU,IPU,K,IXYZ)
                ELSE
                  DIF = DIPSP(IPU,IQU,K,IXYZ)
     &                + DIPSP(IQU,IPU,K,IXYZ)
                END IF
!
                IF( ABS(DIF) .GT. THR ) THEN
                  KERR = KERR + 1
                  WRITE(LUPRI,3000) IXYZ,IPU,IQU,K,
     &              DIPSP(IPU,IQU,K,IXYZ),
     &              DIPSP(IQU,IPU,K,IXYZ)
                END IF
!
              END DO
            END DO
          END DO
        END DO
      END DO
!
      IF( KERR .GT. 0 ) THEN
        WRITE(LUPRI,"('## ERROR!  Symmetry of DIPSP is broken.')")
        STOP
      END IF
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck ckachr */
      SUBROUTINE CKACHR(LUPRI,A8CHAR,IXYZ)
!----------------------------------------------------------------------
!     Check charaters representing the kind of property integrals.
!
!     Written by:  S. Yamamoto - 2007.09.13, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.13, chukyo-u
!     Called by :  CLTMOM1
!----------------------------------------------------------------------
#include "implicit.h"
      CHARACTER*8 A8CHAR
!----------------------------------------------------------------------
      IF( IXYZ .EQ. 1 ) THEN
        IF( A8CHAR .NE. 'XDIPLEN') THEN
          WRITE(LUPRI,"(/1X,'Invalid property file',
     &                  /1X,'XDIPLEN not included',A8)") A8CHAR
          STOP
        END IF
!
      ELSE IF( IXYZ .EQ. 2 ) THEN
        IF( A8CHAR .NE. 'YDIPLEN') THEN
          WRITE(LUPRI,"(/1X,'Invalid property file',
     &                  /1X,'YDIPLEN not included',A8)") A8CHAR
          STOP
        END IF
!
      ELSE IF( IXYZ .EQ. 3 ) THEN
        IF( A8CHAR .NE. 'ZDIPLEN') THEN
          WRITE(LUPRI,"(/1X,'Invalid property file',
     &                  /1X,'ZDIPLEN not included',A8)") A8CHAR
          STOP
        END IF
      END IF
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck mkia2ca */
      SUBROUTINE MKIA2CA(NCORSP,IA2CA,IC2CA)
!----------------------------------------------------------------------
!     Making an index mapping tables IA2CA and IC2CA.
!     IA2CA converts active spinor index to core+active spinor index.
!     IC2CA converts core spinor index to core+active spinor index.
!
!     Written by:  S. Yamamoto - 2007.09.14, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.14, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.23, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.24, chukyo-u
!     Called by :  SY_TMOM
!----------------------------------------------------------------------
!  N.B.  Strucutre of IA2CA (sya, 2007.09.23)
!     NFSYM=2            ! number of fermion irreps
!     SYA_INACT=(/2,1/)  ! inactive (core) shell (MO)
!     SYA_ACTIV=(/3,2/)  ! active shell (MO)
!     NSPACT=10          ! total number of active spinors
!
!     IFSYM=1                        IFSYM=2        
!     -----------------------------  ------------------------
!     p-unbar        p-bar           p'-unbar     p'-bar
!     core  active   core  active    core active  core active
!     01 02 03 04 05 06 07 08 09 10  11   12 13   14   15 16  (CA)
!            |  |  |        |  |  |        |  |         |  | 
!            |  |  |        |  |  |        |  |         |  | 
!            1  2  3        4  5  6        7  8         9 10  (A)
!                                                             IA2CA
!      ^ ^            ^ ^             ^            ^
!      | |            | |             |            |
!      1 2            3 4             5            6          (C)
!                                                             IC2CA
!----------------------------------------------------------------------
#include "implicit.h"
!
!     ...for NFSYM
#include "dgroup.h" 
!
!     ...for NSPACT
#include "goscom.h"
!
!     ...for SYA_INACT,SYA_ACTIV
#include "cossya.h" 
!
      integer IA2CA(NSPACT)
      integer IC2CA(NCORSP)
!
      INTEGER IKPAIR, IFSYM, IC_MO, IA_MO, ICA_SP, IA_SP
!
      DATA NKPAIR/2/
!----------------------------------------------------------------------
      CALL IZERO(IA2CA,NSPACT)
      CALL IZERO(IC2CA,NCORSP)
!
      ICA_SP = 0
      IA_SP = 0
      IC_SP = 0
!
      DO IFSYM = 1, NFSYM
        DO IKPAIR = 1, NKPAIR
          DO IC_MO = 1, SYA_INACT(IFSYM)
            ICA_SP = ICA_SP + 1
            IC_SP = IC_SP + 1
            IC2CA(IC_SP) = ICA_SP
          END DO
          DO IA_MO = 1, SYA_ACTIV(IFSYM)
            ICA_SP = ICA_SP + 1
            IA_SP = IA_SP + 1
            IA2CA(IA_SP) = ICA_SP
          END DO
        END DO
      END DO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck cltmom2 */
      SUBROUTINE CLTMOM2(
     &  KFREE,      LFREE,      IPRINT,     NRIVEC,     N3DIM,
     &  NACTMO,     NOCCMO,     NOCCSP,     NCORMO,     NCORSP,
     &  TMOM,       TDLEN,      TDIPE,      DIPSP,      TODM,
     &  IA2CA,      IC2CA,      TRDM,       WORK
     &  )
!----------------------------------------------------------------------
!     Reading active-part Kramers-pair based transition density matrix 
!     TRDM from the direct-access file.
!
!     Written by:  S. Yamamoto - 2007.09.13, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.13, chukyo-u
!     Called by :  SY_TMOM
!----------------------------------------------------------------------
#include "implicit.h"
!
!     ...for NZ
#include "dgroup.h"  
!
!     ...for MSTATE
#include "cossya.h"
!
!     ...for NSPACT
#include "goscom.h"
#include "priunit.h"
!
      DIMENSION TMOM(NRIVEC,N3DIM)
      DIMENSION TDLEN(NRIVEC,N3DIM)
      DIMENSION TDIPE(NRIVEC,N3DIM)
      DIMENSION DIPSP(NOCCSP,NOCCSP,NRIVEC,N3DIM)
      DIMENSION TODM(NOCCSP,NOCCSP,NRIVEC)
      DIMENSION IA2CA(NSPACT)
      DIMENSION IC2CA(NCORSP)
      DIMENSION TRDM(NSPACT,NSPACT,NRIVEC)
      DIMENSION WORK(*)
!
      LOGICAL DOIMAG
      LOGICAL IEQJ
!
      DATA LUTRDM/11/  ! logical unit number for reading TRDM
      DATA KREAL/1/,KIMAG/2/
!----------------------------------------------------------------------
      IF( IPRINT .GE. 5 ) WRITE(LUPRI,"(' >>>> CLTMOM2')")
!
!     Making IA2CA and IC2CA mapping tables
!     =====================================
!
      CALL MKIA2CA(NCORSP,IA2CA,IC2CA)
!
      IF( NZ .GT. KREAL ) THEN
        KRANDI = KIMAG
        DOIMAG = .TRUE.
      ELSE
        KRANDI = KREAL
        DOIMAG = .FALSE.
      END IF
!
      L1 = NSPACT * NSPACT * NRIVEC
!
!     WRITE(LUPRI,'(1X)')
!
!     ...Open TRDM file
!
      CALL OP2TRDM(LUTRDM,NREP,NSPACT,NZ,NREC)
!
!     ===================== 
!     Loop for TRDM records
!     ===================== 
!
      IREC = 1
      DO 
        CALL DZERO(TRDM,L1)
!
!       ...Read (spinor-base) TRDM for a set of states
!
        IREC = IREC + 1
        IF( IREC .GT. NREC ) EXIT
        READ(LUTRDM,REC=IREC)
     &    IRRP_I, IRRP_J, IVEC_I, IVEC_J,
     &    ((TRDM(I,J,KREAL),I=1,NSPACT),J=1,NSPACT)
        IF( DOIMAG ) THEN
          IREC = IREC + 1
          READ(LUTRDM,REC=IREC)
     &      IRRP_I, IRRP_J, IVEC_I, IVEC_J,
     &      ((TRDM(I,J,KIMAG),I=1,NSPACT),J=1,NSPACT)
        END IF
!
        IF( IRRP_I .EQ. IRRP_J .AND. IVEC_I .EQ. IVEC_J ) THEN
          IEQJ = .TRUE.
        ELSE
          IEQJ = .FALSE.
        END IF
!
!       ...dump TRDM
!
        IF( IPRINT .GE. 10 ) THEN
          WRITE(*,"(/1X,'Dump TRDM real part')")
          CALL OUTPUT(TRDM(1,1,KREAL),1,NSPACT,1,NSPACT,
     &      NSPACT,NSPACT,1,LUPRI)
          IF( NZ. GE. KIMAG ) THEN
            WRITE(*,"(/1X,'Dump TRDM imaginary part')")
            CALL OUTPUT(TRDM(1,1,KIMAG),1,NSPACT,1,NSPACT,
     &        NSPACT,NSPACT,1,LUPRI)
          END IF
        END IF
!
!
!       Making TODM
!       ===========
!
        CALL CLTMOM3(
     &    IPRINT, NOCCSP, NRIVEC, NCORMO, NCORSP,
     &    KRANDI, IEQJ,   
     &    TODM,   IA2CA,  TRDM,   IC2CA
     &    )
!
!       Calculating trace of [TODM . DIPSP]
!       ===================================
!
        CALL CLTMOM4(
     &    KFREE,      LFREE,      IPRINT,     NRIVEC,     N3DIM,
     &    NOCCMO,     NOCCSP,     IEQJ,
     &    TDLEN,      DIPSP,      TODM,
     &    WORK
     &  )
!
!       Calculating total dipole moments
!       ================================
!
        CALL CLTMOM5(
     &    IPRINT, NRIVEC, N3DIM, 
     &    IRRP_I, IRRP_J, IVEC_I, IVEC_J, IEQJ,
     &    TMOM,   TDLEN,  TDIPE  
     &    )
!
      END DO
!
      IF( IPRINT .GE. 5 ) WRITE(LUPRI,"(' <<<< CLTMOM2')")
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck cltmom3 */
      SUBROUTINE CLTMOM3(
     &    IPRINT, NOCCSP, NRIVEC, NCORMO,  NCORSP,
     &    KRANDI, IEQJ,   
     &    TODM,   IA2CA,  TRDM,   IC2CA
     &    )
!----------------------------------------------------------------------
!     Obtaining TODM by inserting the inactive (core) part into TRDM.
!     The inactive part is diagonal unities for the case of permanent
!     dipole moments, and zero's for the the case of 
!     transition dipole moments.
!
!     Written by:  S. Yamamoto - 2007.09.14, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.14, chukyo-u
!     Called by:   CLTMOM2
!----------------------------------------------------------------------
!  Structure of TRDM (sya, 2007.09.13)
! 
!              p-unbar  p-bar
!            (A-string) (B-string)
!               active  active
!              *-------+-------*
!     p-unbar  | ##### | $$$$$ |
!       active | ##### | $$$$$ |
!              | ##### | $$$$$ |
!              *-------+-------*
!     p-bar    | %%%%% | @@@@@ |
!       active | %%%%% | @@@@@ |
!              | %%%%% | @@@@@ |
!              *-------+-------*
!
!  Structure of TODM (sya, 2007.09.13)
!     In the case of diagonal (IEQJ), x is set to 1.
!     In the case of off-diagonal (NOT.IEQJ), x is set to 0.
!
!               p-unbar       p-bar
!               (A-string)    (B-string)
!               core  active  core active
!              *-----+-------+-----+-------*
!     p-unbar  | x   | 0     | 0   | 0     |
!       core   |   x |       |     |       |
!              +-----+-------+-----+-------+
!              | 0   | ##### | 0   | $$$$$ |
!       active |     | ##### |     | $$$$$ |
!              |     | ##### |     | $$$$$ |
!              +-----+-------+-----+-------+
!     p-bar    | 0   | 0     | x   | 0     |
!       core   |     |       |   x |       |
!              +-----+-------+-----+-------+
!              | 0   | %%%%% | 0   | @@@@@ |
!       active |     | %%%%% |     | @@@@@ |
!              |     | %%%%% |     | @@@@@ |
!              *-----+-------+-----+-------*
!
!----------------------------------------------------------------------
#include "implicit.h"
#include "goscom.h"
#include "priunit.h"
      DIMENSION IA2CA(NSPACT)
      DIMENSION IC2CA(NCORSP)
      DIMENSION TODM(NOCCSP,NOCCSP,NRIVEC)
      DIMENSION TRDM(NSPACT,NSPACT,NRIVEC)
!
      LOGICAL IEQJ
!
      DATA KREAL/1/, KIMAG/2/
      DATA D1/1.0D0/
!----------------------------------------------------------------------
      IF( IPRINT .GE. 6 ) WRITE(LUPRI,"(' >>>> CLTMOM3')")
      NOCCMO = NOCCSP / 2
      L1 = NOCCSP * NOCCSP * NRIVEC
      CALL DZERO(TODM,L1)
!
      DO K = 1, KRANDI
        DO JA_SP = 1, NSPACT
          JCA_SP = IA2CA(JA_SP)
          DO IA_SP = 1, NSPACT
            ICA_SP = IA2CA(IA_SP)
            TODM(ICA_SP,JCA_SP,K) = TRDM(IA_SP,JA_SP,K)
          END DO
        END DO
      END DO
!
      IF( IEQJ ) THEN
        DO IC_MO = 1, NCORSP
          IC_SP = IC2CA(IC_MO)
          TODM(IC_SP,IC_SP,KREAL) = D1
        END DO
      END IF
!
!     ...dump TODM
!
      IF( IPRINT .GE. 10 ) THEN
        WRITE(*,"(/1X,'Dump TODM real part')")
        CALL OUTPUT(TODM(1,1,KREAL),1,NOCCSP,1,NOCCSP,
     &    NOCCSP,NOCCSP,1,LUPRI)
        IF( KRANDI .EQ. KIMAG ) THEN
          WRITE(*,"(/1X,'Dump TODM imaginary part')")
          CALL OUTPUT(TODM(1,1,KIMAG),1,NOCCSP,1,NOCCSP,
     &      NOCCSP,NOCCSP,1,LUPRI)
        END IF
      END IF
!
      IF( IPRINT .GE. 6 ) WRITE(LUPRI,"(' <<<< CLTMOM3')")
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck cltmom4 */
      SUBROUTINE CLTMOM4(
     &  KFREE,      LFREE,      IPRINT,     NRIVEC,     N3DIM,
     &  NOCCMO,     NOCCSP,     IEQJ,
     &  TDLEN,      DIPSP,      TODM,
     &  WORK
     &  )
!***********************************************************************
!     Calculating the trace of (TODM * DIPSP).
!     <x> = DDOT [ TODM . DIPSP ]
!     This is a complex arithmetic.
!
!     Written by:  S. Yamamoto - 2007.05.12, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.15, chukyo-u
!     Called by:   CLTMOM2
!----------------------------------------------------------------------
!  N.B. (sya, 2007.09.15)
!     (a + bi) * (c + di) = ac - bd + (ad + bc)i
!***********************************************************************
#include "implicit.h"
!
!     ...for NZ
#include "dgroup.h"  
!
!     ...for LUPRI
#include "priunit.h"  
!
      LOGICAL   IEQJ
      DIMENSION TDLEN(NRIVEC,N3DIM)
      DIMENSION DIPSP(NOCCSP,NOCCSP,NRIVEC,N3DIM)
      DIMENSION TODM(NOCCSP,NOCCSP,NRIVEC)
      DIMENSION WORK(*)
!
      EXTERNAL DDOT
!
      DATA KREAL/1/,KIMAG/2/
      DATA IONE/1/
      DATA THR/1.0D-10/
!----------------------------------------------------------------------
      IF( IPRINT .GE. 6 ) WRITE(LUPRI,"(' >>>> CLTMOM4')")
      L1 = NOCCSP * NOCCSP
!
!
!       ...Loop for x, y, and z
!
      DO IXYZ = 1, N3DIM
! 
        VALRR = DDOT(L1,
     &    TODM(1,1,KREAL),IONE,
     &    DIPSP(1,1,KREAL,IXYZ),IONE)
! 
        VALRI = DDOT(L1,
     &    TODM(1,1,KREAL),IONE,
     &    DIPSP(1,1,KIMAG,IXYZ),IONE)
! 
        VALIR = 0.0D0
        VALII = 0.0D0
! 
        IF( NZ .GT. KREAL ) THEN
! 
          VALIR = DDOT(L1,
     &      TODM(1,1,KIMAG),IONE,
     &      DIPSP(1,1,KREAL,IXYZ),IONE)
! 
          VALII = DDOT(L1,
     &      TODM(1,1,KIMAG),IONE,
     &      DIPSP(1,1,KIMAG,IXYZ),IONE)
! 
        END IF
!
        TDLEN(KREAL,IXYZ) = VALRR - VALII
        TDLEN(KIMAG,IXYZ) = VALRI + VALIR
!
      END DO
!
!     ...Check realness of dipole length operator matrix.
!
      IF( IEQJ ) THEN
        KERR = 0
        DO IXYZ = 1, N3DIM
          IF( ABS(TDLEN(KIMAG,IXYZ)) .GT. THR ) THEN
            KERR = KERR + 1
            WRITE(LUPRI,"(/1X,'## ERROR: not real diagonal ',
     &        'dipole length operator'/1X,I5,F20.10)")
     &        IXYZ,TDLEN(KIMAG,IXYZ)
          END IF
        END DO
        IF( KERR .GT. 0 ) STOP
      END IF
!
      IF( IPRINT .GE. 6 ) WRITE(LUPRI,"(' <<<< CLTMOM4')")
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck cltmom5 */
      SUBROUTINE CLTMOM5(
     &  IPRINT, NRIVEC, N3DIM, 
     &  IRRP_I, IRRP_J, IVEC_I, IVEC_J, IEQJ,
     &  TMOM,   TDLEN,  TDIPE  
     &  )
!-----------------------------------------------------------------------
!     Calculation of the total dipole moments from the components.
!
!     Written by:  S. Yamamoto - 2007.06.04, chukyo-u
!     Last-update: S. Yamamoto - 2007.06.19, chukyo-u
!     Called by:   CLTMOM2
!-----------------------------------------------------------------------
#include "implicit.h"
#include "codata.h"
!
!     ...for MXCOOR
#include "mxcent.h"
!     ...for DIPMN
#include "dipole.h"
!
!     ...for LUPRI
#include "priunit.h"  
!
!
!             ...declare variable(s) in argument list
!
      DIMENSION TMOM(NRIVEC,N3DIM)
      DIMENSION TDLEN(NRIVEC,N3DIM)
      DIMENSION TDIPE(NRIVEC,N3DIM)
      LOGICAL IEQJ
!
!             ...declare constant(s) in code
!
      DATA IONE/1/
      DATA DM1/-1.0D0/
      DATA KREAL/1/
!
      EXTERNAL DDOT
!-----------------------------------------------------------------------
!
      IF( IPRINT .GE. 6 ) WRITE(LUPRI,"(' >>>> CLTMOM5')")
      L1 = NRIVEC * N3DIM
!
!     ...Electronic contribution for dipole moments.
!        Dipole length expectation values multiplied by -e.
!
      CALL DCOPY(L1,TDLEN,IONE,TDIPE,IONE)
      CALL DSCAL(L1,DM1,TDIPE,IONE)
!
!     ...Adding nuclear contribution for permanent dipole moments.
!        The nuclear part does not contribute to transition dipole 
!        moments.
!
      IF( IEQJ ) THEN
        CALL DIPNUC(0,.FALSE.)
        DO IXYZ = 1, N3DIM
          TMOM(KREAL,IXYZ) = TDIPE(KREAL,IXYZ) + DIPMN(IXYZ)
        END DO
      ELSE
        CALL DCOPY(L1,TDIPE,IONE,TMOM,IONE)
      END IF
!
!     ...Calculate the norm in atomic units.
      TDPAU = DDOT(L1,TMOM,IONE,TMOM,IONE)
      TDPAU = SQRT(TDPAU)
!
!#    TDPDB = TDPAU * CVEL * ECHARGE * BOHR
      TDPDB = TDPAU * DEBYE
!
!
!     Printing TMOM
!     =============
!
      CALL PRTMOM(
     &  IPRINT, NRIVEC, N3DIM, 
     &  IRRP_I, IRRP_J, IVEC_I, IVEC_J, IEQJ,
     &  TMOM,   TDLEN,  TDIPE,  TDPAU,  TDPDB )
!
      IF( IPRINT .GE. 6 ) WRITE(LUPRI,"(' <<<< CLTMOM5')")
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck prtmom */
      SUBROUTINE PRTMOM(
     &  IPRINT, NRIVEC, N3DIM, 
     &  IRRP_I, IRRP_J, IVEC_I, IVEC_J, IEQJ,
     &  TMOM,   TDLEN,  TDIPE,  TDPAU,  TDPDB 
     &  )
!----------------------------------------------------------------------
!     Printing TMOM (transition dipole moments).
!
!     Written by:  S. Yamamoto - 2007.05.12, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.15, chukyo-u
!     Called by:   CLTMOM2
!----------------------------------------------------------------------
#include "implicit.h"
#include "goscom.h"
#include "priunit.h"
      DIMENSION TMOM(NRIVEC,N3DIM)
      DIMENSION TDLEN(NRIVEC,N3DIM)
      DIMENSION TDIPE(NRIVEC,N3DIM)
      LOGICAL IEQJ
!
      DATA KREAL/1/, KIMAG/2/
!----------------------------------------------------------------------
 1000 FORMAT( 1X,3F18.10,' (real)',
     &       /1X,3F18.10,' (imaginary)')
!----------------------------------------------------------------------
      IREP_I = IRRP_I
      IREP_J = IRRP_J
      IF( MOD(NELACT,2) .EQ. 0 ) THEN
        IREP_I = IRRP_I + NREP
        IREP_J = IRRP_J + NREP
      END IF
!
      WRITE(LUPRI,'(/1X,A)') REPEAT('-',70)
      WRITE(LUPRI,'(1X,A,A,I3,1X,A,A,I4, A,I3,1X,A,A,I4,A)')
     &  'TMOM',
     &  ' <irrep',IRRP_I,REPNA(IREP_I),' vec#',IVEC_I,
     &  ' |x,y,z| irrep',IRRP_J,REPNA(IREP_J),' vec#',IVEC_J,'>'
!
!     ...Dipole length expectation values
!
      IF( IPRINT .GE. 5 ) THEN
        WRITE(LUPRI,"(1X,' Dipole length operators')")
        WRITE(LUPRI,1000)
     &    (TDLEN(KREAL,IXYZ),IXYZ=1,N3DIM),
     &    (TDLEN(KIMAG,IXYZ),IXYZ=1,N3DIM)
      END IF
!
!     ...Electronic contributions to total dipole moments.
!        (-1) multiplied to dipole length operators.
!        Printed only for permanent dipole moments.
!        Not printed for transition dipole moments.
!
      IF( IEQJ ) THEN
        WRITE(LUPRI,"(1X,
     &    ' Electronic contributions to dipole moments')")
        WRITE(LUPRI,1000)
     &    (TDIPE(KREAL,IXYZ),IXYZ=1,N3DIM),
     &    (TDIPE(KIMAG,IXYZ),IXYZ=1,N3DIM)
      END IF
!
!     ...Total dipole moments (nuclear contributions added)
!
      IF( IEQJ ) THEN
        WRITE(LUPRI,"(1X,
     &    ' Total permanent dipole moments ',
     &    '(nuclear contributions added)')")
      ELSE
        WRITE(LUPRI,"(1X,' Transition dipole moments')")
      END IF
      WRITE(LUPRI,1000)
     &  (TMOM(KREAL,IXYZ),IXYZ=1,N3DIM),
     &  (TMOM(KIMAG,IXYZ),IXYZ=1,N3DIM)
!
!     ...Norm in atomic units and in Debye unit.
!
      WRITE(LUPRI,"(1X,' Norm in atomic units and in Debye'
     &             /1X,2F18.10)")
     &  TDPAU,TDPDB
      WRITE(LUPRI,'(1X,A)') REPEAT('-',70)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck mempat */
      SUBROUTINE MEMPAT(
     &  MESG,IPRINT,IUNIT,NMPART,M1,M2,
     &  CC,KK,LL,WW,KFREE,LFREE)
!----------------------------------------------------------------------
!     Memory partitioning.
!
!     Written by: S. Yamamoto - May 28, 2007, chukyo-u
!----------------------------------------------------------------------
#include "implicit.h"
      CHARACTER*(*) MESG
      CHARACTER*4 CC(NMPART)
      DIMENSION KK(NMPART)
      DIMENSION LL(NMPART)
      DIMENSION WW(*)
!----------------------------------------------------------------------
      IF( M1 .LT. M2 ) THEN
        DO I = M1, M2
          CALL MEMGET(CC(I),KK(I),LL(I),WW,KFREE,LFREE)
        END DO 
        KK(M2+1) = KFREE
      END IF 
!
      IF( IPRINT .GE. 20 ) THEN
        WRITE(IUNIT,'(/1X,"Memory partitioning for calling ",A)') MESG
        WRITE(IUNIT,'(5I10)') (KK(I),I=1,(M2+1))
      END IF
      RETURN
      END
C
C     ===========================================================
      SUBROUTINE read_fcidump (mode,n1,norb,n,irpamo,
     & ECORE,EPS,ONER,ONEI,
     & INDTWR,INDTWI,TWOINR,TWOINI,MAXINT)
C     ==========================================================
      IMPLICIT none
C
C     Routine to read data from an FCIDUMP formatted file
C     Mode 0: header and symmetry check
C     Mode 1: one-electron integrals
C     Mode 2: two-electron integrals
C
      real(8), parameter   :: D0 = 0.0D0
      integer, parameter   :: FCIDUMP = 75 ! file unit
      integer, intent(in)  :: n1,n,mode
      integer, intent(out) :: norb,irpamo(N1)
      real(8), intent(out) :: ecore,eps(N1)
      integer, intent(in)  :: maxint
      integer, intent(out) :: INDTWR(N,N,N,N)
      integer, intent(out) :: INDTWI(N,N,N,N)
      real(8), intent(out) :: ONER(N1,N1),ONEI(N1,N1)
      real(8), intent(out) :: TWOINR(0:MAXINT),TWOINI(0:MAXINT)

      real*8 integral(2)
      integer nelec,isym,orbsym(32),iuhf,n_uniq,label(4)
      namelist/FCI/norb,nelec,orbsym,isym,iuhf

      orbsym = 0
      open(FCIDUMP,file='FCIDUMP',form='FORMATTED')
      read(FCIDUMP,FCI)
      if (MODE == 0) then
         write(6,FCI)
         if (norb > n1) then
            print*, " Error: GOSCIP can not handle",norb, " spinors"
            call quit ('dimensioning error')
         end if
         if (any (orbsym(1:norb) /= 1)) then
            print*, " Error: this interface does not handle symmetry"
            call quit ('incomplete symmetry information')
         end if
         irpamo(1:norb) = orbsym(1:norb)
      endif

      n_uniq = 0
      do
        read (FCIDUMP,*) integral,label
        if (all(label == 0)) then
           if (mode == 1) ecore = integral(1)
        elseif (all(label(2:4) == 0)) then
           if (mode == 1) eps(label(1)) = integral(1)
        elseif (all(label(3:4) == 0)) then
           if (mode == 1) oner(label(1),label(2)) = integral(1)
           if (mode == 1) onei(label(1),label(2)) = integral(2)
        else
           n_uniq = n_uniq + 1
           if (mode == 2) then
               TWOINR(n_uniq) = integral(1)
               TWOINI(n_uniq) = integral(2)
               indtwr(label(1),label(2),label(3),label(4)) = n_uniq
               indtwi(label(1),label(2),label(3),label(4)) = n_uniq
               indtwr(label(3),label(4),label(1),label(2)) = n_uniq
               indtwi(label(3),label(4),label(1),label(2)) = n_uniq
               indtwr(label(2),label(1),label(4),label(3)) = n_uniq
               indtwi(label(2),label(1),label(4),label(3)) = -n_uniq
               indtwr(label(4),label(3),label(2),label(1)) = n_uniq
               indtwi(label(4),label(3),label(2),label(1)) = -n_uniq
           end if
        end if
        if (all(label == 0)) exit
      end do

      close(FCIDUMP,status="KEEP")

      end
