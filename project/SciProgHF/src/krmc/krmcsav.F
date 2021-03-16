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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmcsav */
      SUBROUTINE KRMCSAV(IOPT,CMO,CREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Save info on LUKRMC
C
C     Input:
C        IOPT    - option (current only 0 supported)
C
C     Written by J. Thyssen - Nov 26 1998
C     Last revision :
C                S. Knecht - Feb 2010: print of mj-vector to file 
C                                      KRMCSCF added.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbopt.h"
#include "dcbgen.h"
C
      DIMENSION CMO(*), CREF(*)
      DIMENSION WORK(*)
      LOGICAL   FNDLAB, FND, RTRACI
C
      CHARACTER CNOLBL*6
C
      CALL QENTER('KRMCSAV')
      KFRSAV = KFREE
C
C     Write file KRMCSCF:
C
C     (1) MJVEC .or. IBEIG: orbital mj-vector or boson irreps orbital-vector
C     (2) OLDORB  : old MO orbitals
C     (3) RESTART : info needed for restart
C     (4) LREDUCED: info needed for backsteps
C     (5) CREF    : new CI vector
C     (6) NEWORB  : new MO orbitals
C
      if(.not.spinfr_krmc)then
        IF( OPT_CHCKJZ .AND. LINEAR )THEN
          CALL MEMGET('INTE',KMJVEC,NORBT,WORK,KFREE,LFREE)
          CALL IZERO(WORK(KMJVEC),NORBT)
          CALL IREAKRMC(LUKRMC,'MJVEC   ',WORK(KMJVEC),NORBT)
        END IF
      else
        CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
        CALL IZERO(WORK(KIBEIG),NORBT)
        CALL IREAKRMC(LUKRMC,'IBEIG   ',WORK(KIBEIG),NORBT)
      end if
C
C     Clean LUKRMC
C     ------------
C
      REWIND(LUKRMC)
      FND = FNDLAB('SODLABEL',LUKRMC)
      BACKSPACE LUKRMC
      CALL NEWLAB('SODLABEL',LUKRMC,LUPRI)
      REWIND(LUKRMC)

      if(.not.spinfr_krmc)then
C
C       (1).1 MJVEC: write orbital mj-vector
C       -------------------------------------------
C
        IF( OPT_CHCKJZ .AND. LINEAR )THEN
          CALL IWRTKRMC(LUKRMC,'MJVEC   ',WORK(KMJVEC),NORBT)
          CALL MEMREL('KRMCSAV.1',WORK,1,KFRSAV,KFREE,LFREE)
        END IF

      else
C
C       (1).2 IBEIG: write boson irreps orbital-vector
C       -------------------------------------------
C
        CALL IWRTKRMC(LUKRMC,'IBEIG   ',WORK(KIBEIG),NORBT)
        CALL MEMREL('KRMCSAV.1',WORK,1,KFRSAV,KFREE,LFREE)
      end if
C
C     (2) OLDORB & OLDCREF: write old MO orbitals
C     -------------------------------------------
C
      CALL WRTKRMC(LUKRMC,'OLDORB  ',CMO,NCMOTQ)
      CALL WRTKRMC(LUKRMC,'OLDCREF ',CREF,NZCONFQ)
C
C     (2) RESTART: write restart info
C     -------------------------------
C
      CALL LABKRMC(LUKRMC,'RESTART ')
      WRITE(LUKRMC)
     &     EMCOLD, DEPRED, BETA, GAMMA, RTRUST, ITMAC,
     &     EMY, EACTIV, EMCOLD, EMCSCF
C
C     (3) LREDUCED: info nededed for backsteps
C     ----------------------------------------
C
      REWIND(LUKRM4)
      READ(LUKRM4) NREDM
C
      CALL MEMGET('REAL',KREDE,NREDM*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDS,NREDM*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDGP,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVECR,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBTYP,2*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEVC,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBPVC,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBCVC,NREDM,WORK,KFREE,LFREE)
C
C
      REWIND(LUKRM4)
      CALL REALRED(LUKRM4,NZRED,NREDM,WORK(KREDE),WORK(KREDGP),
     &     WORK(KEVECR),WORK(KIBTYP),
     &     WORK(KIBEVC),WORK(KIBPVC),WORK(KIBCVC),NPRED)
C
      CALL LABKRMC(LUKRMC,'LREDUCED')
      CALL WRTLRED(LUKRMC,NZRED,NREDM,WORK(KREDE),WORK(KREDGP),
     &     WORK(KEVECR),WORK(KIBTYP),
     &     WORK(KIBEVC),WORK(KIBPVC),WORK(KIBCVC),NPRED)
C
      CALL LABKRMC(LUKRMC,'EOFLABEL')
C
C     (4) - (5) CREF and NEWORB
C     -------------------------
C
C     Rotate to get new CREF and CMO
C
      EVALR = D0
      CALL GTNORB(CMO,CREF,WORK(KIBTYP),
     &     WORK(KIBCVC),WORK(KIBEVC),WORK(KIBPVC),
     &     EVALR,WORK(KEVECR),WORK,KFREE,LFREE)
C
C     Rotate to Fock type and/or natural orbitals.
C
      CALL MEMGET('REAL',KEIG  ,NORBT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KEIG),NORBT)
      IF (IKRMCCNO .GE. -1) THEN
CSK   IF (.FALSE.) THEN
C     TODO Feb 2004 hjaaj: KRMCCNO not working correctly yet ...
C
C          July 2008 sk+hjaaj : fix for FOCKON only
C                               still TODO: implement RTRACI + 
C                                           fix TRTOCO/TRTONO branch
C
         IF (NATONL) THEN
            CNOLBL = 'ONLYNO'
         ELSE IF (FOCKON) THEN
            CNOLBL = 'ONLYFD'
         ELSE
            CNOLBL = 'FD+NO'
         END IF
C
         IF (IKRMCCNO .EQ. -1 .OR. NZCONF .LE. 1) THEN
            RTRACI = .FALSE.
         ELSE
            RTRACI = .TRUE.
         END IF
C
         WRITE (LUPRI,'(/2A)') ' - calling KRCCNO with keyword ',CNOLBL
         CALL KRMCCNO(CNOLBL,RTRACI,CREF,CMO,IPROPT,WORK(KEIG),
     &                WORK,KFREE,LFREE)
C
      END IF
C
      CALL WRTKRMC(LUKRMC,'CREF    ',CREF,NZCONFQ)
      CALL WRTKRMC(LUKRMC,'NEWORB  ',CMO,NCMOTQ)
C
      IF (NASHT .EQ. 0) THEN
C     ... also save on DFCOEF when SCF
C         (needed for e.g. RPA or restarts without changing
C          input from .HF to .MCSCF /hjaaj Feb 2004)
         CALL OPNFIL(LUCOEF,'DFCOEF','UNKNOWN','KRMCSAV')
         CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
         CALL IZERO(WORK(KIBEIG),NORBT)
         CALL WRICMO(LUCOEF,CMO,WORK(KEIG),WORK(KIBEIG),EMCSCF)
         CLOSE(LUCOEF,STATUS='KEEP')
      END IF
C
      REWIND(LUKRMC)
C
      CALL MEMREL('KRMCSAV.2',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('KRMCSAV')
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck save_moltra_inp */
      SUBROUTINE SAVE_MOLTRA_INP(SAVTRA,SAVLOG,SAVREAL,SAVINT,IPRINT)
C***********************************************************************
C
C    routine to save common block for MOLTRA
C
C     Written by JVP 2020
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "../moltra/dcbtra.h"
C
      CHARACTER SAVTRA(4,2,2)*(*)
      LOGICAL SAVLOG(*)
      REAL*8  SAVREAL(*)
      INTEGER SAVINT(*)
C
C     /DCBTC/
C
      DO I = 1, 4
         DO J = 1, 2
            SAVTRA(I,J,1)= TRA4_INDSTR(I,J)
         END DO
      END DO
      DO I = 1,2
         DO J = 1,2
            SAVTRA(I,J,2) = TRA2_INDSTR(I,J)
         END DO
         SAVTRA(3,I,2) = TRA_INDSTR(I)
         SAVTRA(4,I,2) = TRA_CORSTR(I)
      END DO
C
C     /DCBTRI/
C
      SAVINT(1) = ITRA_INTFLG 
      SAVINT(2) = ITRA_INTFL2 
      SAVINT(3) = ITRA_INTFL4
      SAVINT(4) = IPRTRA
      SAVINT(5) = ISTRAT
C
C     /DCBTRL/
C
      SAVLOG(1) = TRA_ANTIS
      SAVLOG(2) = TRA_MSOUT
      SAVLOG(3) = NO2IND
      SAVLOG(4) = NO4IND
      SAVLOG(5) = NOPAIR
      SAVLOG(6) = PRPTRA
      SAVLOG(7) = RCORBS
      SAVLOG(8) = FKRMC
      SAVLOG(9) = NOMDCINT
      SAVLOG(10) = MDCSCAT
C
C     /DCBTRR/
C
      SAVREAL(1) = THROUT
      SAVREAL(2) = SCRTRA
C
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reset_moltra_inp */
      SUBROUTINE RESET_MOLTRA_INP(SAVTRA,SAVLOG,SAVREAL,SAVINT,IPRINT)
C***********************************************************************
C
C    routine to reset common block for MOLTRA
C
C     Written by JVP 2020
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "../moltra/dcbtra.h"
C
      CHARACTER SAVTRA(4,2,2)*(*)
      LOGICAL SAVLOG(*)
      REAL*8  SAVREAL(*)
      INTEGER SAVINT(*)
C
C     /DCBTC/
C
      DO I = 1, 4
         DO J = 1, 2
            TRA4_INDSTR(I,J) = SAVTRA(I,J,1)
         END DO
      END DO
      DO I = 1,2
         DO J = 1,2
            TRA2_INDSTR(I,J) = SAVTRA(I,J,2)
         END DO
         TRA_INDSTR(I) = SAVTRA(3,I,2)
         TRA_CORSTR(I) = SAVTRA(4,I,2)
      END DO
C
C     /DCBTRI/
C
      ITRA_INTFLG = SAVINT(1)
      ITRA_INTFL2 = SAVINT(2)
      ITRA_INTFL4 = SAVINT(3)
      IPRTRA = SAVINT(4)
      ISTRAT = SAVINT(5)
C
C     /DCBTRL/
C
      TRA_ANTIS = SAVLOG(1)
      TRA_MSOUT = SAVLOG(2)
      NO2IND = SAVLOG(3)
      NO4IND = SAVLOG(4)
      NOPAIR = SAVLOG(5)
      PRPTRA = SAVLOG(6)
      RCORBS = SAVLOG(7)
      FKRMC  = SAVLOG(8)
      NOMDCINT = SAVLOG(9)
      MDCSCAT = SAVLOG(10)
C
C     /DCBTRR/
C
      THROUT = SAVREAL(1)
      SCRTRA = SAVREAL(2)
C
C
      RETURN
      END
