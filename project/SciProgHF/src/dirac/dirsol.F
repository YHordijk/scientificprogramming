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
C/* Deck solfck */
      SUBROUTINE SOLFCK(FOCK,DMAT,NFMAT,ESOLVE,ESOLVN,
     &                  WORK,LWORK,IPRINT)
C***********************************************************************
C
C     This routine will generate necessary solvent integrals and add
C     corresponding contributions to the Fock matrix
C
C     Written May 27 1998 by T.Saue and K.Mikkelsen
C
C***********************************************************************      
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
C
#include "dcbbas.h"
#include "dcbham.h"
C
      DIMENSION FOCK(*),DMAT(*),WORK(*)
#include "memint.h"
      CALL QENTER('SOLFCK')
C
      MAXYZ = (LSOLMX + 1)*(LSOLMX + 2)/2
      MAXLM  = 2*LSOLMX + 1
C
C     Get nuclear contribution to solvent polarization
C
      LMTOT  = (LSOLMX+1) ** 2
      LMNTOT = (LSOLMX+1)*(LSOLMX+2)*(LSOLMX+3) / 6
      CALL MEMGET('REAL',KRNUC,LMTOT ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTLMN,LMNTOT,WORK,KFREE,LFREE)
C     - Calculate TLMN = RNUC in Cartesian basis
      CALL SOLNUC(.FALSE.,0,.FALSE.,WORK(KTLMN),DUMMY,DUMMY,
     &            WORK,KFREE,LFREE,IPRINT)
C     - Transform TLMN to RNUC(l,m) in spherical basis
      CALL TLMTRA(1,WORK(KTLMN),WORK(KRNUC),WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('SOLFCK.nuc',WORK,KWORK,KTLMN,KFREE,LFREE)
C
      CALL MEMGET('REAL',KLMBUF,NNBBASX*MAXYZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLMTRI,NNBBASX      ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLMMAT,N2BBASX      ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCSPH ,MAXLM*MAXYZ  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIFUN ,N2BBASX      ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLVAL ,MAXYZ        ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KMVAL ,MAXYZ        ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNVAL ,MAXYZ        ,WORK,KFREE,LFREE)
      CALL SOLFC1(FOCK,DMAT,NFMAT,ESOLVE,ESOLVN,WORK(KLMTRI),
     &            WORK(KLMBUF),WORK(KLMMAT),
     &            WORK(KRNUC),WORK(KCSPH),WORK(KIFUN),
     &            WORK(KLVAL),WORK(KMVAL),WORK(KNVAL),
     &            MAXYZ,MAXLM,WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('SOLFCK',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('SOLFCK')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck solfc1 */
      SUBROUTINE SOLFC1(FOCK,DMAT,NFMAT,ESOLVE,ESOLVN,
     &            RLMTRI,RLMBUF,RLMMAT,RLMNUC,
     &            CSPH,IFUN,LVALUE,MVALUE,NVALUE,
     &            MAXYZ,MAXLM,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     This routine will generate necessary solvent integrals and add
C     corresponding contributions to the Fock matrix
C
C     Written May 27 1998 by T.Saue and K.Mikkelsen
C
C***********************************************************************      
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0,D2 = 2.0D0)
      PARAMETER(FLFAC = -0.5D0,DM2 = -2.0D0)
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "symmet.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "orgcom.h"
#include "dgroup.h"
C
      LOGICAL DOINT(2,2),FIRST
      DIMENSION FOCK(N2BBASX,NZ,NFMAT),DMAT(N2BBASX,NZ,NFMAT),
     &          RLMTRI(*),RLMBUF(NNBBASX,*),
     &          LVALUE(*),MVALUE(*),NVALUE(*),
     &          RLMMAT(*),RLMNUC(*),CSPH(MAXYZ,MAXLM),IFUN(*),WORK(*)
      DIMENSION DIPORG_save(3), IDOINT(2,2)
      CHARACTER SECTID*12,INTTID*12
#include "ibtfun.h"
C
C     Initialization
C     ==============
C     Note that Fock-matrix is not zeroed....
C
C
C     MORDER eq 0 : M_l order 0, +1, -1, +2, -2, ..., +LVAL, -LVAL
C     MORDER ne 0 : M_l order -LVAL, ..., -1, 0, 1, ..., +LVAL
C     MINTEG eq 1 : the spherical components are normalized such
C                   that all coefficients are integers.
C
      MORDER = 1
      MINTEG = 0
      DOINT(1,1) = .TRUE.
      DOINT(2,2) = .TRUE.
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
      IDOINT(1,1) = 1
      IDOINT(2,2) = 1
      IDOINT(1,2) = 0
      IDOINT(2,1) = 0
      ESOLVE = D0
      ESOLVN = D0
      CALL DCOPY(3,DIPORG,1,DIPORG_save,1)
      CALL DCOPY(3,CAVORG,1,DIPORG,1)
C     
C     Reindex real part of density matrix to unsorted basis
C
      DO I = 1,NFMAT
        CALL BSTOBU(DMAT(1,1,I),1,WORK,LWORK)
        CALL BSTOBU(FOCK(1,1,I),1,WORK,LWORK)
      ENDDO
C
C     Loop over orders in spherical expansion
C     =======================================
C
      IOFF = 0
      DO 100 L = 0,LSOLMX
        LL   = L + 1
        NXYZ = (L + 1)*(L + 2)/2
        NLM  = 2*L + 1
C
C       2a. Cartesian functions
C       -----------------------
C
        KXYZ = 0
        CALL LMNVAL(LL,NXYZ,LVALUE,MVALUE,NVALUE)
        DO ICMP = 1,NXYZ
          LVAR = MOD(LVALUE(ICMP),2)*ISYMAX(1,1)
          MVAR = MOD(MVALUE(ICMP),2)*ISYMAX(2,1)
          NVAR = MOD(NVALUE(ICMP),2)*ISYMAX(3,1)
          IREP = IBTXOR(LVAR,IBTXOR(MVAR,NVAR))
          IF(IREP.EQ.0) THEN
            KXYZ       = KXYZ + 1
            IFUN(KXYZ) = ICMP
          ENDIF
        ENDDO
        IF(KXYZ.EQ.0) GOTO 110
        ESOLLE = D0
        ESOLLN = D0
C
C       F_l factor
C
        RFL = FLFAC*RSOLAV**(-(2*L+1))
     &       *(L+D1)*(EPSOL-D1)/(L+EPSOL*(L+D1))
C
C       *********************************
C       ***  Electronic contribution  ***
C       *********************************
C
C
        CALL GETTIM(CPU1,WALL1)
C
        CALL RGETINT('CARMOM ',RLMBUF,DOINT,NXYZ,L,IDUMMY,
     &            'TRIANG',.FALSE.,IPRINT,WORK,LWORK)
C
        CALL GETTIM(CPU2,WALL2)
        CPU    = CPU2 - CPU1
        WALL   = WALL2 - WALL1
        INTTID = SECTID(CPU)
        WRITE(LUPRI,'(A,I2,A,A12)') 
     &     'SOLFCK: L = ',L,' - CPU : ',INTTID
        WRITE(LUPRI,'(A,I5)')
     &     ' Totally symmetric Cartesians:',KXYZ
C
C       Get spherical transformation matrix
C       ===================================
C
        CALL SPHCOM(L,CSPH,MAXLM,MAXYZ,MORDER,MINTEG,
     &                WORK,LWORK,IPRINT)
C
C       Loop over totally symmetric RLM contributions
C       =============================================
C
        DO ILM = 1,NLM
          M = ILM - LL
          IF(IREPLM(L,M).EQ.0) THEN
C
C             Spherical transformation
C             ------------------------
C
            FIRST = .TRUE.
            DO I = 1,KXYZ
              IXYZ = IFUN(I)
              FAC = CSPH(IXYZ,ILM)
              IF(FAC.NE.D0) THEN
                IF(FIRST) THEN
                  CALL DCOPY (NNBBASX,RLMBUF(1,IXYZ),1,RLMTRI,1)  
                  IF(FAC.NE.D1) THEN
                    CALL DSCAL(NNBBASX,FAC,RLMTRI,1)
                  ENDIF
                  FIRST = .FALSE.
                ELSE
                  CALL DAXPY(NNBBASX,FAC,RLMBUF(1,IXYZ),1,RLMTRI,1)
                ENDIF
              ENDIF
            ENDDO
C
C            Expand matrices
C            ---------------
C
            CALL QNPHASE(1,0,RLMTRI,IDOINT)
            CALL DSPTSI(NTBAS(0),RLMTRI,RLMMAT)
C
C            Scale with expectation value ++
C            -------------------------------
C
            DO I = 1,NFMAT
              RLMELC = DM2*DDOT(N2BBASX,DMAT(1,1,I),1,RLMMAT,1) 
              RLMTOT = RLMELC + RLMNUC(IOFF+ILM)
              FAC = DM2*RFL*RLMTOT
              ESOLLE = ESOLLE + RFL*RLMTOT*RLMELC
              ESOLLN = ESOLLN + RFL*RLMTOT*RLMNUC(IOFF+ILM)
              IF(IPRINT.GE.3) THEN
                 WRITE(LUPRI,'(A,I2,A1,I2,A1)') 
     &           '** Solvent contributions for (l,m) = (',
     &           L,',',M,')'
                WRITE(LUPRI,'(3X,A,F20.12)')
     &           '* F_l factor                     : ',RFL,
     &           '* Expectation value              : ',RLMTOT,
     &           '    - nuclear contribution       : ',
     &                RLMNUC(IOFF+ILM),
     &           '    - electronic contribution    : ',
     &                RLMELC,
     &           '* Solvent energy                 : ',
     &              RFL*RLMTOT*RLMTOT,
     &           '    - nuclear contribution       : ',
     &              RFL*RLMTOT*RLMNUC(IOFF+ILM),
     &           '    - electronic contribution    : ',
     &              RFL*RLMTOT*RLMELC
                IF(IPRINT.GE.6) THEN
                WRITE(LUPRI,'(/A/)') 
     &           '* Electronic solvent integrals:'
                CALL OUTPAK(RLMTRI,NTBAS(0),NTBAS(0),LUPRI)
                ENDIF
              ENDIF
C
C              Add to Fock matrix
C              ------------------
C
              CALL DAXPY(N2BBASX,FAC,RLMMAT,1,FOCK(1,1,I),1)
            ENDDO
          ENDIF
        ENDDO
        IOFF = IOFF + NLM
        ESOLVE = ESOLVE + ESOLLE
        ESOLVN = ESOLVN + ESOLLN
        IF(IPRINT.GE.3) THEN
           WRITE(LUPRI,'(A,I2)') 
     &     '*** Total Solvent contributions for L = ',L
           WRITE(LUPRI,'(3X,A,F20.12)')
     &           '* Solvent energy                 : ',ESOLLE+ESOLLN,
     &           '    - nuclear contribution       : ',ESOLLN,
     &           '    - electronic contribution    : ',ESOLLE
        ENDIF
  110   CONTINUE
        IOFF = IOFF + NLM
  100 CONTINUE
      CALL DCOPY(3,DIPORG_save,1,DIPORG,1)
C
C     Reindex to sorted basis
C     -----------------------
C
      DO I = 1,NFMAT
        CALL BUTOBS(DMAT(1,1,I),1,WORK,LWORK)
        CALL BUTOBS(FOCK(1,1,I),1,WORK,LWORK)
      ENDDO
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(3X,A,F20.12)')
     &    '* Total solvent energy           : ',
     &       ESOLVE+ESOLVN,
     &    '    - nuclear contribution       : ',
     &       ESOLVN,
     &    '    - electronic contribution    : ',
     &       ESOLVE
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck solvin */
      SUBROUTINE SOLVIN(WORD,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxmom.h"
#include "mxcent.h"
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 8)
      LOGICAL SET,CAVUSR,ADDUSR,RADUSR
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     &          MXFORM*6,FMT*6
      DIMENSION WORK(*)
#include "dcbham.h"
#include "orgcom.h"
#include "cbisol.h"
#include "cbiher.h"
#include "nuclei.h"
C
      SAVE SET

      DATA TABLE /'.CAVITY','.DIELEC','.MAX L ','.PRINT ',
     *            '.INERSI','.INERSF','.CAVORG','.CAVADD'/
      DATA SET/.FALSE./
      IF (SET) THEN
         RETURN
      END IF
      SET = .TRUE.
#include "memint.h"
C
C     Initialize
C
      RSOLAV = DUMMY
      EPSOL  = DUMMY
      EPSTAT = DUMMY
      LSOLMX = 10
      IPRSOL = 0
      INERSI = .FALSE.
      INERSF = .FALSE.
C
      CAVUSR = .FALSE.
      RADUSR = .FALSE.
      ADDUSR = .FALSE.
C
C     ***** INPUT FORMAT STATEMENTS *****
C
      INPEPS = 0
      WORD1 = WORD
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',6)
               GO TO 100
            END IF
            WRITE (6,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
               GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 9999
         ELSE
            WRITE (6,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',6)
         CALL QUIT(' ILLEGAL KEYWORD IN SOLINP ')
C
C *** Option 1 >CAVITY<  Cavity radius
  101 CONTINUE
        READ(LUCMD,*) RSOLAV
        RADUSR = .TRUE.
      GO TO 100
C
C *** Option 2 >DIELEC<  Dielectric constant
  102 CONTINUE
      READ(LUCMD,*) EPSOL
      EPSTAT = EPSOL
      INPEPS = INPEPS + 1
      GO TO 100
C
C *** Option 3 >MAX L <  Maximum L quantum number
  103 CONTINUE
      READ(LUCMD,*) LSOLMX
      NLMSOL = (LSOLMX+1)*(LSOLMX+1)
      GO TO 100
C
C *** Option 4 >PRINT <  Print level
  104 CONTINUE
      READ(LUCMD,*) IPRSOL
      GO TO 100
C
C *** Option 5 >INERSI<  Initial state defining inertial pol.
  105 CONTINUE
      INERSI = .TRUE.
      INPEPS = INPEPS + 1
      READ(LUCMD,*) EPSOL
      EPSTAT = EPSOL
      GO TO 100
C
C *** Option 6 >INERSF<  Final state with inertial pol. from
C                        previous calc. with ".INERSI"
  106 CONTINUE
      INPEPS = INPEPS + 1
      INERSF = .TRUE.
      READ(LUCMD,*) EPSTAT, EPSOL
      GO TO 100
C
C *** Option 7 >CAVORG<  Specify origo of solvent cavity
C                        
  107 CONTINUE
      READ (LUCMD,*) (CAVORG(I),I = 1, 3)
      CAVUSR = .TRUE.
      GO TO 100
C
  108 CONTINUE
C
C *** Option 8 >CAVADD< Add radius to distance from cavity origin to 
C     the most distant center
      READ (LUCMD,*) CAVADD    
      ADDUSR = .TRUE.
      GO TO 100
C
 9999 CONTINUE
      IF (INPEPS .GT. 1) THEN
         WRITE (LUCMD,'(/1X,A/A)') WORD1//' input error',
     *   ' Only one of ".DIELEC", ".INERSI", and ".INERSF"'//
     *   ' may be specified.'
         CALL QUIT('SOLVIN: Very ambigiuous...')
      END IF
      IF(.NOT.ADDUSR.AND.RSOLAV.EQ.DUMMY) THEN
        WRITE(LUCMD,'(A)')
     &   'Cavity radius NOT specified !'
        CALL QUIT('SOLVIN: No cavity radius.')
      ENDIF
      IF(EPSOL.EQ.DUMMY) THEN
        WRITE(LUCMD,'(A)')
     &   'No dielectric constant specified !'
        CALL QUIT('SOLVIN: No dielectric constant.')
      ENDIF
      IF(INERSF.AND.EPSTAT.EQ.DUMMY) THEN
        WRITE(LUCMD,'(A)')
     &   'No static dielectric constant specified !'
        CALL QUIT('SOLVIN: No static dielectric constant.')
      ENDIF
      IF((LSOLMX+1).GT.MXQNM) THEN
        WRITE(LUCMD,'(A,I5)')
     &   'Chosen LSOLMX value larger than max ',MXQNM
        CALL QUIT('SOLVIN: LSOLMX too large.')
      ENDIF
      IF(ADDUSR.AND.RADUSR) THEN
        WRITE(LUCMD,'(A,I5)')
     &   'Both CAVITY and CAVADD specified ..'
        CALL QUIT('SOLVIN: CAVADD and CAVITY.....')
      ENDIF
C
C     Process section
C
      IF (.NOT.CAVUSR) THEN
         CALL MEMGET('REAL',KGEOM,3*NATOMS,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KMASS,  NATOMS,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KNAT ,  NATOMS,WORK,KFREE,LFREE)
         CALL CMMASS(WORK(KGEOM),WORK(KMASS),WORK(KNAT),ISOTOP,IPRHAM)
         CALL DCOPY(3,CMXYZ,1,CAVORG,1)
         CALL MEMREL('SOLVIN',WORK,KWORK,KWORK,KFREE,LFREE)
      END IF
      IF(ADDUSR) THEN
        RSOLAV = DISCEN(CAVORG) + CAVADD        
      ENDIF
C
C     Print section
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)') 'SOLVIN: Solvent model: '
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A,I3)') '* Print level: ',IPRSOL
      WRITE(LUPRI,'(1X,A)') 
     &     '* Molecular system enclosed in spherical cavity'
      FMT = MXFORM(RSOLAV,10)
      WRITE(LUPRI,'(1X,A,'//FMT//',A)')'   of radius ',RSOLAV,' a.u.'
      WRITE(LUPRI,'(1X,A)') ' in dielectric medium'
      IF(CAVUSR) THEN
        WRITE(LUPRI,'(/,A,3F12.6)') ' Cavity center:',
     &                               (CAVORG(I), I = 1, 3)
      ELSE
        WRITE(LUPRI,'(/,A,3F12.6)') ' Cavity center (center of mass):',
     &                               (CAVORG(I), I = 1, 3)
      ENDIF
      WRITE(LUPRI,'(1X,A/1X,A,I3)') 
     &  '* The polarization potential is represented by ',
     &  '  a solid harmonic expansion to order L = ',LSOLMX
      FMT = MXFORM(EPSOL,10)
      WRITE(LUPRI,'(1X,A,'//FMT//')')
     &     ' Dielectric constant: ',EPSOL
      IF(INERSI) THEN
        WRITE(LUPRI,'(A)') 
     &  '* Initial state calculation.',
     &  '  Inertial polarization vector written to file'
      ENDIF
      IF(INERSF) THEN
        WRITE(LUPRI,'(A)') 
     &  '* Final state calculation.',
     &  '  Use of polarization vectors'
      ENDIF
C
C     Transfer information to CBISOL in HERMIT
C
C        transfer information to CBISOL
C        NCNTCV is defined in READIN
        SOLVNT = .TRUE.
        LCAVMX = LSOLMX
        LMTOT  = (LCAVMX+1) ** 2
        LMNTOT = (LCAVMX+1)*(LCAVMX+2)*(LCAVMX+3) / 6
        RCAV(1)   = RSOLAV
        RCAV(2)   = RSOLAV
        RCAV(3)   = RSOLAV
        EPDIEL    = EPSOL
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck discen */
      FUNCTION DISCEN(ORIGIN)
C***********************************************************************
C
C     Determine the most distant center with respect to a
C     given origin
C
C     Written by T.Saue Jun 19 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "symmet.h"
      DIMENSION ORIGIN(3)
C
#include "ibtfun.h"
C
C     First find the most distant center
C
      DISCEN = D0
      DO ICENTA = 1, NUCIND
         DO IA = 0, MAXOPR
         IF (IBTAND(IA,ISTBNU(ICENTA)) .EQ. 0) THEN
           CXA = PT(IBTAND(ISYMAX(1,1),IA))*CORD(1,ICENTA) - ORIGIN(1)
           CYA = PT(IBTAND(ISYMAX(2,1),IA))*CORD(2,ICENTA) - ORIGIN(2)
           CZA = PT(IBTAND(ISYMAX(3,1),IA))*CORD(3,ICENTA) - ORIGIN(3)
           RADBUF = CXA*CXA+CYA*CYA+CZA*CZA
           DISCEN = MAX(DISCEN,RADBUF)
         ENDIF
         ENDDO
      ENDDO
      DISCEN = SQRT(DISCEN)
C
      RETURN
      END

      
