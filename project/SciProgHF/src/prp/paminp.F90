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

! FILE    : paminp.F
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck prpinp */
      SUBROUTINE PRPINP(WRDSRC,WORK,LWORK)
!***********************************************************************
!
! <<< Property Input for DIRAC >>>
!
!
!    Called from:  PAMINP - first reading of the input
!                  MINPRP - called when optimize is on 
!                  (properties at final geometry)
!
!       Written by Trond Saue and Hans Joergen Aa.Jensen
!       Last revision: June 14 1996 - tsaue
!
!***********************************************************************
         use memory_allocator
#ifdef MOD_LAO_REARRANGED
         use london_helper
#endif
         use dirac_cfg

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 10, NTABLE = 41)
#include "maxorb.h"
#include "mxcent.h"
#include "dcbham.h"
#include "dcbgen.h"
!
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,      &
     &          WORD1*7, WRDSRC*7
      DIMENSION WORK(LWORK)
      LOGICAL CHGO
      double precision, allocatable :: GEOM(:,:),MASS(:)
      integer, allocatable :: NAT(:),ISOT(:)
!
#include "orgcom.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dgroup.h"
#include "nuclei.h"
!
      allocate(GEOM(NUCDEP,3))
      allocate(MASS(NATOMS))
      allocate(NAT(NATOMS))
      allocate(ISOT(NATOMS))
!
      DATA TABDIR /'*END OF','*EXPECT','*LINEAR','*EXCITA',             &
     &             '*MOLGRD','*NMR   ','*ESR   ','*QUADRA',             &
     &             '*STEX  ','*OPENRS'/
!                     1         2         3         4
      DATA TABLE  /'.DIPOLE','.QUADRU','.POLARI','.MAGNET',             &
!                     5         6         7         8
     &             '.NMR   ','.SPIN-S','.SHIELD','.PRINT ',             &
!                     9        10        11        12
     &             '.ABUNDA','.EFG   ','.NQCC  ','.DSO   ',             &
!                    13        14        15        16
     &             '.MOLGRD','.NSTDIA','.ESR   ','.PVC   ',             &
!                    17        18        19        20
     &             '.EFT   ','.EFTNTL','.PVCSHI','.EFN   ',             &
!                    21        22        23        24
     &             '.WAVE F','.RHONUC','.FIRST ','.TWO-PH',             &
!                    25        26        27        28
     &             '.VERDET','.ROTG  ','.MSTDIA','.PVCSPN',             &
!                    29        30        31        32
     &             '.NOPCTR','.EXCDIP','.MAGOUT','.OPTROT',             &
!                    33        34        35        36
     &             '.C6    ','.SPIN-R','.BXLAO ','.BYLAO ',             &
!                    37        38        39        40
     &             '.BZLAO ','.EFFDEN','.EFFDE2','.SHIEL2',             &
!                    41        42        43        44
     &             '.RDCCDE'/

!
      CALL QENTER('PRPINP')
#include "memint.h"
!
!     Initialize /CBIPRP/
!     ===================
!
      IPRPRP       = 0
      DOEXP        = .FALSE.
      DOESR        = .FALSE.
      DOXLR        = .FALSE.
      DOXPP        = .FALSE.
      DOXQR        = .FALSE.
      DOTPA        = .FALSE.
      DOVER        = .FALSE.
      DOEXCPRP     = .FALSE.
      EXCPRPTYP    = 0
      DIPOLE       = .FALSE.
      QUADRU       = .FALSE.
      EPOLAR       = .FALSE.
      DOVDW        = .FALSE.
      MSUSCP       = .FALSE.
      IMAGOUT      = 0
      SHIELD       = .FALSE.
      SPINRO       = .FALSE.
      PVC_SHIELD   = .FALSE.
      PVC_SPINSPIN = .FALSE.
      SPNSPN       = .FALSE.
      EFG          = .FALSE.
      EFN          = .FALSE.
      EFT          = .FALSE.
      EFTNTL       = .FALSE.
      OPTROT       = .FALSE.
      NQCC         = .FALSE.
      DSO          = .FALSE.
      NSTDIA       = .FALSE.
      MSUSCDIA     = .FALSE.
      ROTG         = .FALSE.
      LONDON       = .FALSE.
      PVC          = .FALSE.
      RHONUC       = .FALSE.
      ABUND        = D1
      NDSOINT      = 0
      MGRAD        = .FALSE.
      ESRGTENS     = .FALSE.
      ESR_HFCC     = .FALSE.
      NPRP_WF      = 1
      PRP_WF(1)    = 'DHF '
      BXLAO        = .FALSE.
      BYLAO        = .FALSE.
      BZLAO        = .FALSE.
      EFFDEN       = .FALSE.
      EFFDE2       = .FALSE.

!MI Activate picture change transformations
!   of four-component property matrices 
!   for the two-component (BSS) infinite order.
! TODO with TS/HJJ: BSS spin-free picture change transformation
!     set compulsory default values for X2C module.
      NOPCT = .FALSE.
      IF (X2C .or.                                                      &
     &   (BSS.AND.(IBSS.GE.99.OR.IBSS.EQ.009.OR.IBSS.GE.109))) THEN
        DOPCT = .TRUE.
!       ... by default DO the picture change transformation
      ELSE
        DOPCT = .FALSE.
      ENDIF

!
!     Read menu file
!     ==============
!
      REWIND (LUCMD,IOSTAT=IOS)
!     ... IOSTAT to avoid program abort on some systems
!         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      CALL UPCASE(WORD)
      IF (WORD .EQ. WRDSRC) THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
!        Tell user which property input is missing
!        (**PROPE or **PRP F or ?)
         WRITE(LUPRI,'(/3A/A)')                                         &
     &   ' * INFORMATION: No property input "',WRDSRC,                  &
     &   '" found in input file.'
         CALL QEXIT('PRPINP')
         RETURN
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no property input found')
  930 CONTINUE
      WORD1 = WORD
!
!     Process input for COMMON  /CBIPRP/
!     ==================================
!
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,                  &
     &                109,110,111,112,113,114,115,116,                  &
     &                117,118,119,120,121,122,123,124,                  &
     &                125,126,127,128,129,130,131,132,                  &
     &                133,134,135,136,137,138,139,140,                  &
     &                141), I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')                                    &
     &         ' Keyword ',WORD,' not recognized in PRPINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in PRPINP.')
  101    CONTINUE
!&&&& DIPOLE - dipole moment
            DIPOLE   = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  102    CONTINUE
!&&&& QUADRUpole moment
            DOEXP   = .TRUE.
            QUADRU  = .TRUE.
            GO TO 100
  103    CONTINUE
!&&&& POLAR - electronic dipole polarizabilities
            DOXLR  = .TRUE.
            EPOLAR = .TRUE.
            GO TO 100
  104    CONTINUE
!&&&& MAGNET - magnetic susceptibilities
            DOXLR  = .TRUE.
            MSUSCP = .TRUE.
            GO TO 100
  105    CONTINUE
!&&&& NMR parameters:
            DOXLR  = .TRUE.
            SHIELD = .TRUE.
            SPNSPN = .TRUE.
            GO TO 100
  106    CONTINUE
!&&&& SPIN-S: NMR indirect spin-spin coupling constants:
            DOXLR  = .TRUE.
            SPNSPN = .TRUE.
            GO TO 100
  107    CONTINUE
!&&&& SHIELD: NMR shielding constants:
            DOXLR = .TRUE.
            SHIELD = .TRUE.
            GO TO 100
  108    CONTINUE
!&&&& PRINT: print level
            READ(LUCMD,*) IPRPRP
            GO TO 100
  109    CONTINUE
!&&&& ABUND: minimum abundance in % of atomic isotopes to be included
            READ (LUCMD,*) ABUND
            GO TO 100
  110    CONTINUE
!&&&& EFG - Electric field gradients
            EFG      = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  111    CONTINUE
!&&&& NQCC - Nuclear Quadrupole Coupling Constants
            NQCC     = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  112    CONTINUE
!&&&& DSO - nonrel. diamagnetic spin-orbit contribution to indirect spin-spin couplings
            DSO   = .TRUE.
            DOEXP = .TRUE.
            GO TO 100
  113    CONTINUE
!&&&& MGRAD - molecular gradient
            MGRAD = .TRUE.
            GO TO 100
  114    CONTINUE
!&&&& NSTDIA - nonrel. diamagnetic contribution to nuclear shielding tensor
            NSTDIA   = .TRUE.
            DOEXP = .TRUE.
            GO TO 100
  115    CONTINUE
!&&&& ESR parameters
            DOESR = .TRUE.
            ESRGTENS = .TRUE.
            ESR_HFCC = .TRUE.
            GO TO 100
  116    CONTINUE
!&&&& PVC Parity violation - chirality
            PVC   = .TRUE.
            DOEXP = .TRUE.
            GO TO 100
  117    CONTINUE
!&&&& EFT - traceless electric field third derivatives
            EFT      = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  118    CONTINUE
!&&&& EFTNTL - non-traceless electric field third derivatives
            EFTNTL   = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  119    CONTINUE
!&&&& PV contribution to NMR shielding
            DOXLR      = .TRUE.
            PVC_SHIELD = .TRUE.
            IF(NBSYM .GT. 1) CALL QUIT('.PVCSHI only for C1 symmetry')
            GOTO 100
  120    CONTINUE
!&&&& EFN: Electric field at nuclei
            DOEXP = .TRUE.
            EFN   = .TRUE.
            GOTO 100
  121    CONTINUE
!&&&& WAVE FUNCTIONS - define which wave functions we should use
            IF (NPRP_WF .GT. 0) THEN
               WRITE(LUPRI,'(//2A/A/,(10X,A4))')                        &
     &         ' WARNING, list of wave functions for properties '//     &
     &         'reset under ',WRDSRC,' Old list which is discarded:',   &
     &         (PRP_WF(I), I=1, NPRP_WF)
            END IF
            READ(LUCMD,*) NPRP_WF
            IF (NPRP_WF .LT. 1 .OR. NPRP_WF .GT. MXPRP_WF) THEN
               WRITE(LUPRI,'(/3A,I10)') ' *** ERROR for ',WRDSRC,       &
     &            ' *** Illegal value for .WAVE FU : ',NPRP_WF
               CALL QUIT('*** ERROR in PRPINP ***')
            END IF
            DO I = 1, NPRP_WF
               READ(LUCMD,'(1X,A4)') PRP_WF(I)
            END DO
            GOTO 100
 122     CONTINUE
!&&&& RHONUC - electronic density at nuclei
            RHONUC   = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
 123     CONTINUE
!&&&& FIRST HYPERPOLARIZABILITY - compute complete beta tensor
            DOXQR = .TRUE.
            CALL DEF_1STHYP(WORK,LWORK)
            GO TO 100
 124     CONTINUE
!&&&& TWO-PHOTON ABSORPTION - compute all two-photon absorption matrix elements
            CALL DEF_TPA(WORK,LWORK,DOTPA)
            GO TO 100
 125     CONTINUE
!&&&& .VERDET - Calculates the Verdet constant with the QR module
            DOVER = .TRUE.
            DOXQR = .TRUE.
            CALL DEF_VERDET(WORK,LWORK)
            GO TO 100
 126     CONTINUE
!&&&& ROTG - rotational g-tensor according to Aucar et al
            ROTG  = .TRUE.
            DOXLR = .TRUE.
            GO TO 100
  127    CONTINUE 
!&&&& MSTDIA - diamagnetic contribution to the magnetizabilities
            MSUSCDIA = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  128    CONTINUE 
!&&&& PV contribution to spin-spin coupling
            DOXLR        = .TRUE.
            PVC_SPINSPIN = .TRUE.
            IF(NBSYM .GT. 1) CALL QUIT('.PVCSPN only for C1 symmetry')
            GO TO 100
  129    CONTINUE 
!&&&& NOPCT - even if DOPCT is .true. (because of default settings) 
!             do NOT activate picture change transformation...
            NOPCT = .TRUE.
            GO TO 100
  130    CONTINUE 
!&&&& .EXCDIP - Calculates excited state expectation values of dipole op.
            CALL EXCPRPINP(1,DOEXCPRP,EXCPRPTYP)
            CALL EXCPRPINP(4,DOEXCPRP,EXCPRPTYP)
            GO TO 100
  131    CONTINUE 
!&&&& .MAGOUT - specify the hard-output for the magnetizability
            READ(LUCMD,*) IMAGOUT
            GO TO 100
  132    CONTINUE 
!&&&& OPTROT - do optical rotation
            DOXLR = .TRUE.
            OPTROT = .TRUE.
            GO TO 100
  133    CONTINUE 
!&&&& .C6 - evaluate C6 dispersion coefficients
            DOXLR  = .TRUE.
            EPOLAR = .TRUE.
            DOVDW  = .TRUE.
            GO TO 100
  134    CONTINUE 
!&&&& SPIN-R - spin-rotation constant according to Aucar et al
            DOEXP = .TRUE.
            DOXLR = .TRUE.
            SPINRO = .TRUE.
            GO TO 100
  135    CONTINUE 
!&&&& BXLAO 
            bxlao  = .true.
            doxlr  = .true.
            london = .true.
            GO TO 100
  136    CONTINUE
!&&&& BYLAO 
            bylao  = .true.
            doxlr  = .true.
            london = .true.
            GO TO 100
  137    CONTINUE
!&&&& BZLAO 
            bzlao  = .true.
            doxlr  = .true.
            london = .true.
            GO TO 100
  138    CONTINUE
!&&&& EFFDEN; calculate effective electronic densities at the nuclei
            DOEXP  = .TRUE.
            EFFDEN = .TRUE. 
            GO TO 100
  139    CONTINUE
!&&&& EFFDE2; calculate effective electronic densities at the nuclei:old code
            DOEXP  = .TRUE.
            EFFDE2 = .TRUE. 
         GO TO 100
  140    CONTINUE
#ifdef MOD_LAO_REARRANGED
!&&&& SHIEL2, calculate shielding+london, but with rearranged terms
            DOXLR  = .TRUE.
            LONDON = .TRUE.
            call set_london_keywords('shiel2')
!gosia: todo now we have two logical variables: SHIEL2 and shielding_rearrange to indicate the same thing
! i will have to fix that...
#endif
            GO TO 100
  141    CONTINUE
! quick hack to read CC density matrix without the need to run CC as wavefunction
            call set_prop_flags_ccsd()
            GO TO 100
!
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 180
      ELSE
         WRITE (LUPRI,'(/3A/)') ' ERROR: Prompt "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in PRPINP.')
      END IF
  180 CONTINUE
!
!     Print section
!     =============
!
      CALL TITLER('Property module','*',129)
      WRITE(LUPRI,'(A,I4)') ' * Print level:',IPRPRP
      WRITE(LUPRI,'(2A)')   ' * Input label: ',WRDSRC
      WRITE(LUPRI,'(A)')                                                &
     &     ' * Properties calculated for the following wave functions:'
      DO I = 1, NPRP_WF
         WRITE(LUPRI,'(I6,2A)') I,': ',PRP_WF(I)
      END DO
!
!=========================================================================
!MI   For magnetic susceptibilities, NMR shielding constants and rotational
!MI   g-tensors, operator centers might be changed
      WRITE(LUPRI,'(1X,A)')                                             &
     &'These initial settings of center and origins'//                  &
     &' might be changed later:'
      WRITE (LUPRI,'(1X,A,3F18.10)')                                    &
     &     '* Operator center (a.u.):', (ORIGIN(I),I=1,3),              &
     &     '* Gauge origin    (a.u.):', (GAGORG(I),I=1,3),              &
     &     '* Dipole origin   (a.u.):', (DIPORG(I),I=1,3)
      GAGORGX = GAGORG(1)
      GAGORGY = GAGORG(2)
      GAGORGZ = GAGORG(3)

      IF (X2C .OR. BSS) THEN
      IF (DOPCT) THEN
        IF (.NOT.NOPCT) THEN
          WRITE(LUPRI,'(A)')                                            &
     & ' * Perform 4c->2c picture change transformation of'//           &
     & ' the four-component property operators'
        ELSE
          WRITE(LUPRI,'(A)')                                            &
     & ' * Do NOT perform 4c->2c picture change transformation of'//    &
     & ' the four-comp. property operators - take only LL block(s)'
        ENDIF
      ELSE
        WRITE (LUPRI,'(A)')                                             &
     & ' * Do NOT perform 4c->2c picture change transformation of'//    &
     & ' the four-component property operators'
      ENDIF
      ENDIF
!
!     Process input for various program sections
!     ==========================================
!
!     Now read input
!
  200 CONTINUE
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4,5,6,7,8,9,10), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/3A/)') ' Directory "',WORD,'" is not valid'
         CALL PRTAB(NDIR,TABDIR,WORD1//' directory keywords',LUPRI)
         CALL QUIT('Illegal directory keyword in PRPINP.')
      ELSE
         WRITE (LUPRI,'(/3A/)')                                         &
     &      ' Prompt in keyword "',WORD,'" illegal or out of order'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in PRPINP, error in prompt.')
      END IF
    2 CONTINUE
         CALL EXPINP(WORD,.FALSE.,WORK,LWORK)
         GO TO 200
    3 CONTINUE
         CALL XLRINP(WORD,.FALSE.,WORK,LWORK)
         GO TO 200
    4 CONTINUE
         CALL XPPINP(WORD,.FALSE.,WORK,LWORK)
         GO TO 200
    5 CONTINUE
         CALL GRDINP(WORD,.FALSE.)
         GO TO 200

!     *NMR
    6 CONTINUE
!MI ... changed for LAO-based magnetic properties
         CALL NMRINP(WORD,.FALSE.,WORK,LWORK)
         GO TO 200

    7 CONTINUE
#ifdef MOD_ESR
         CALL ESRINP(WORD,.FALSE.,WORK,LWORK)
#else
         CALL QUIT('*ESR    not included in this version')
#endif
         GO TO 200
    8 CONTINUE
         CALL XQRINP(WORD,.FALSE.)
         GO TO 200
!     *STEX
    9 CONTINUE
         CALL STEX_INPUT(WORD,.FALSE.,WORK,LWORK)
         GO TO 200

!     *OPENRSP    
   10 continue
         call quit('*OpenRSP not included in this version')
      go to 200

    1 CONTINUE
      IF(DOXLR) CALL XLRINP(WORD,.TRUE.,WORK,LWORK)
!MI ... read again the input
      CALL NMRINP(WORD,.TRUE.,WORK,LWORK)

!gosia: i have to reset these variables once again...
!or always remember to put *NMR .LONDON in the input...
      if (bxlao .or. bylao .or. bzlao) then
        london = .true.
        doxlr  = .true.
      end if

      CALL GRDINP(WORD,.TRUE.)
      IF(DOXPP) CALL XPPINP(WORD,.TRUE.,WORK,LWORK)
      IF(DOXQR) CALL XQRINP(WORD,.TRUE.)
#ifdef MOD_ESR
      IF(DOESR) CALL ESRINP(WORD,.TRUE.,WORK,LWORK)
#endif

!===========================================================================
! Define special operators for the list of properties programmed in DIRAC
!===========================================================================
      CALL PRPDEF(WORK,LWORK)

! .... check the change of the gauge origin ....
      THRS1=1.0D-7
      CHGO= DABS(GAGORGX-GAGORG(1)).GE.THRS1.OR.                        &
     &      DABS(GAGORGY-GAGORG(2)).GE.THRS1.OR.                        &
     &      DABS(GAGORGZ-GAGORG(3)).GE.THRS1
      IF (CHGO) THEN
         WRITE (LUPRI,'(A/A,3F18.10)')                                  &
     &     '   Gauge origin has been changed:',                         &
     &     ' * Gauge origin   :', (GAGORG(I),I=1,3)
      ENDIF

!     consistency checks
      if (twocomp .and. .not. (bss .or. x2c)) then
         if (doxlr .or. doxqr) then
            call quit('ERROR: response with .NONREL not possible')
         end if
      end if
! do not do shielding with open-shells
      if (shield) then
       if (NASHT.GT.0) then
        call quit('ERROR: NMR shielding with open shells not possible')
       end if
      end if

!gosia: probably here is a good place to add Miro's warning, 
! but i need to think of correct flags (for DC and SPINFREE Hamiltonians)
! so commented out for now
!      if (shield .and. .not. urkbal) then
!        write(lupri, *) 'warning: unsufficient flexibility of small '//
!     &  'component set in RKB scheme, you need very large basis sets,'//
!     &  '  or use simple magn.balance instead...'
!      end if

      CALL FLSHFO(LUPRI)
      CALL QEXIT('PRPINP')

      deallocate(ISOT)
      deallocate(NAT)
      deallocate(MASS)
      deallocate(GEOM)

      RETURN
!
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck expinp */
      SUBROUTINE EXPINP(WORD,RESET,WORK,LWORK)
!***********************************************************************
!
!     Input section for expectation values
!
!     Written by T.Saue - May 1996
!     Last revision: May 26  - 1996
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 8, DM1 = -1.0D00, D0 = 0.0D00)
!
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbxpr.h"
!
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
!
      SAVE SET
      DATA TABLE /'.PRINT ','.OPERAT','.ORBANA','.PROJEC',              &
     &            '.PRPCAN','.EPOLE ','.MPOLE ','.XXXXXX'/
      DATA SET/.FALSE./
!
      NEWDEF = (WORD .EQ. '*EXPECT')
      IF (SET) THEN
         IF (NEWDEF)                                                    &
     &      CALL QUIT('Only one "*EXPECT" input section allowed')
!        hjaaj: repeated input sections give infinite loop ...
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
!
!     Initialize /DCBEXP/
!     ===================
!
      IPREXP = 0
      NEXPP  = 0
      ORBANA = .FALSE.
      PRJANA = .FALSE.
      PRPCAN = .FALSE.
!
!     Process input for DCBEXP
!     =========================
!
      INPERR = 0
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,              &
     &            '" not recognized in EXPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in EXPINP.')
    1          CONTINUE
!&&&& PRINT:  Print level
                  READ(LUCMD,*) IPREXP
               GO TO 100
    2          CONTINUE
!&&&& OPERATOR: Define operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,   &
     &                        IPREXP)
                  CALL OP1IND('EXPINP',INDEXP,LEXPP,NEXPP,INDXPR,       &
     &                        MAXEXP)
                  DOEXP = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    3          CONTINUE
!&&&& ORBANA: Analyze expectation value in terms of individual orbitals
               ORBANA = .TRUE.
               GO TO 100
    4          CONTINUE
!&&&& PROJEC: Analyze expectation value in terms of fragment orbitals
               PRJANA = .TRUE.
               GO TO 100
    5          CONTINUE
!&&&& PRPCAN: Generate canonical occupied property orbitals
               PRPCAN = .TRUE.
               GO TO 100
    6          CONTINUE
!&&&& EPOLE : Electric multipoles of order L
               READ(LUCMD,*) IORDER
               NDEG=(IORDER+2)*(IORDER+1)/2
               IF((NEXPP+NDEG).GT.MAXEXP) THEN
                 WRITE(LUPRI,'(A,A,I3)') 'EXPINP: Pointer array LEXPP ',&
     &           'out of bounds for electric multipole of order ',      &
     &           IORDER
                 WRITE(LUPRI,'(A,I5)')                                  &
     &           'Increase MAXEXP to ',(NEXPP+NDEG)
               ENDIF
               CALL DEF_EPOLE(IORDER,LEXPP(NEXPP+1),IPREXP)
               DOEXP=.TRUE.
               NEXPP = NEXPP + NDEG
               GO TO 100
    7          CONTINUE
!&&&& MPOLE : Magnetic multipoles of order L (zero by time reversal symmetry)
               READ(LUCMD,*) IORDER
               NDEG=3*(IORDER+1)*IORDER/2
               IF((NEXPP+NDEG).GT.MAXEXP) THEN
                 WRITE(LUPRI,'(A,A,I3)') 'EXPINP: Pointer array LEXPP ',&
     &           'out of bounds for magnetic multipole of order ',      &
     &           IORDER
                 WRITE(LUPRI,'(A,I5)')                                  &
     &           'Increase MAXEXP to ',(NEXPP+NDEG)
               ENDIF
               CALL DEF_MPOLE(IORDER,LEXPP(NEXPP+1),IPREXP)
               DOEXP=.TRUE.
               NEXPP = NEXPP + NDEG
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,               &
     &            '" not recognized in EXPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in EXPINP.')
            END IF
      END IF
  300 CONTINUE
!
!     Print section
!     =============
!
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' EXPINP: Expectation values'
      IF(ORBANA.OR.PRJANA) THEN
        WRITE(LUPRI,'(A)')                                              &
     &   '* Expectation values analyzed in terms of '
        IF(ORBANA) WRITE(LUPRI,'(3X,A)') '- individual orbitals'
        IF(PRJANA) WRITE(LUPRI,'(3X,A)') '- fragment orbitals'
      ENDIF
      IF(PRPCAN) THEN
        WRITE(LUPRI,'(A)')                                              &
     &     '* Generating canonical orbitals for each property.'
      ENDIF
      CALL PRSYMB(LUPRI,'=',75,0)
      IF(ICHANG.GT.0) THEN
        WRITE(LUPRI,'(A)')                                              &
     &  ' * The following expectation values will be calculated:'
      ENDIF
      IF(NEXPP.GT.0) THEN
        DO I = 1,NEXPP
          INDXPR = LEXPP(I)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      ENDIF  
      IF(DIPOLE) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Dipole moment'
      ENDIF
      IF(QUADRU) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Quadrupole moment'
      ENDIF
      IF(PVC) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Parity violation - chirality'
      ENDIF
      IF(EFN) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Electric field at nuclei'
      ENDIF
      IF(EFG) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Electric field gradient at nuclei'
      ENDIF
      IF(EFT) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')                                           &
     &     '- traceless electric field third derivattives'
      ENDIF
      IF(EFTNTL) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')                                           &
     &     '- non-traceless electric field third derivattives'
      ENDIF
      IF(NQCC) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Nuclear Quadrupole Coupling Constant'
      ENDIF
!
      CALL PRSYMB(LUPRI,'-',75,0)
  999 CONTINUE
      IF (INPERR.GT.0) CALL QUIT('Input error in *EXPECT')
      RETURN
!
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck xprinp */
      SUBROUTINE XPRINP(LUCMDX,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,        &
     &                  IPRINT)
!***********************************************************************
!
! PURPOSE:
!  Assign a new property operator based on the input
!  and the /XCBPRP/ list.
!
!  An operator in 4-component theory is a linear combination
!  of products of a 4x4 matrix and some scalar operator P.
!  The 4x4 matrices are the alpha matrices and the identity.
!  Generic one-electron operators are defined in subroutine MOPDEF
!  XPRINP recognizes the following different kinds of operators:
!
!    1. DIAGONAL
!       P             * scalar operator
!
!    2. XALPHA
!       [alpha_x]P    * x-component of alpha times scalar operator
!
!    3. YALPHA
!       [alpha_y]P    * y-component of alpha times scalar operator
!
!    4. ZALPHA
!       [alpha_z]P    * z-component of alpha times scalar operator
!
!    5. XAVECTOR
!       [alpha x P]_x * vector product of alpha and vector operator,
!                       x-component    
!    6. YAVECTOR
!       [alpha x P]_y * vector product of alpha and vector operator,
!                       y-component    
!    7. ZAVECTOR
!       [alpha x P]_z * vector product of alpha and vector operator,
!                       z-component    
!    8. ALPHADOT
!       A.P           * dot-product of alpha and vector operator
!
!    9. GAMMA5        
!       gamma5 P      * The gamma5 matrix times a scalar operator
!                       gamma5 = | 0 I | 
!                                | I 0 |
!   10. XSIGMA
!       [sigma_x]p    * x-component of sigma times scalar operator
!
!   11. YSIGMA
!       [sigma_y]p    * y-component of sigma times scalar operator
!
!   12. ZSIGMA 
!       [sigma_z]p    * z-component of sigma times scalar operator
!
!   13. XBETASIG
!       [betasig_x]p  * x-component of beta sigma times scalar operator
!
!   14. YBETASIG
!       [betasig_y]p  * Y-component of beta sigma times scalar operator
!
!   15. ZBETASIG
!       [betasig_z]p  * z-component of beta sigma times scalar operator
!
!   16. XiBETAAL
!       i[betaalp_x]p * x-component of beta alpha times scalar operator
! 
!   17. YiBETAAL
!       i[betaalp_y]p * y-component of beta alpha times scalar operator
!
!   18. ZiBETAAL
!       i[betaalp_z]p * z-component of beta alpha times scalar operator
!
!   19. BETA   
!       [beta]P * beta times scalar operator
!
!   20. SIGMADOT 
!       iS.P * dot-product of sigma matrix and vector operator
!
!   21. iBETAGAM
!       i[beta*gamma5] P           * i beta*gamma5 times scalar operator
!
!
!  Operator types 5,6 and 7 form the components of a vectorproduct 
!  of the alpha matrices and some scalar operator.
!
!  The scalar operators may be supplied with factors, e.g. field strengths.
!
!  The pointer IKW controls what is to be read
!     0   Property name
!     1   Keywords
!     2   Property labels (to be found on file)
!     3   Class label
!     4   Separate factor for each component
!     5   Common factor for all components
!
!
!  Copyright (c) 27-Jan-1995 Hans Joergen Aa. Jensen
!  Last revision: May 14 1996 - tsaue
!                 Sept. 2012   M Ilias 
!
!
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
      CHARACTER*7 WORD
!
      PARAMETER (NTABKW = 28)
      CHARACTER*8 TABKW(NTABKW)
      CHARACTER   LINE*80, KEYWD*8, TOKEN16*16, BLANK8*8, BLANK16*16
      DIMENSION PFAC(3)
      CHARACTER   PLABEL(3)*8, PNAME*16
      LOGICAL     CLABEL
!
      PARAMETER  (D1 = 1.0D0)
      PARAMETER (BLANK8='        ', BLANK16='                ' )
      SAVE TABKW
      DATA TABKW/'CLASS   ', 'FACTORS ', 'COMFACTO', 'CMULT   ',        &
     &           'FOCKMAT ', 'DIAGONAL', 'XALPHA  ', 'YALPHA  ',        &
     &           'ZALPHA  ', 'XAVECTOR', 'YAVECTOR', 'ZAVECTOR',        &
     &           'ALPHADOT', 'GAMMA5  ', 'XSIGMA  ', 'YSIGMA  ',        &
     &           'ZSIGMA  ', 'XBETASIG', 'YBETASIG', 'ZBETASIG',        &
     &           'XiBETAAL', 'YiBETAAL', 'ZiBETAAL', 'BETA    ',        &
     &           'SIGMADOT', 'iBETAGAM', 'XXXXXXXX', 'XXXXXXXX'/
!
      CALL QENTER('XPRINP')
!
!     Initialization
!     ==============
!
      INDXPR = -1
      ISYXPR = -1
      ITRXPR = -1
      PFAC(1) = D1
      PFAC(2) = D1
      PFAC(3) = D1
      CLABEL  = .FALSE.
      !mi - position of the Fock operator in the TABKW list, or "border" index
      I_FOCKMAT=5
!
      IKW = 0
!
!     *****************
!     *** Read line ***
!     *****************
!
  100 READ (LUCMDX,'(A)',END=8000) LINE
      IF (IPRINT .GT. 20) THEN
         WRITE (LUPRI,'(A)') '***** Input line read in XPRINP:'
         WRITE (LUPRI,'(A)') LINE
      END IF
!
!     Line to be skipped
!     ==================
      IF (LINE(1:1) .EQ. '!' .OR. LINE(1:1) .EQ. '#') GO TO 100
!
!     Line is keyword from other section: analyze what you already have
!     =================================================================
!
      IF (LINE(1:1) .EQ. '*' .OR. LINE(1:1) .EQ. '.') THEN
         WORD = LINE(1:7)
         BACKSPACE LUCMDX
         GO TO 500
      END IF
!
!     Line must start with a blank
!     ============================
!
      IF (LINE(1:1) .NE. ' ') THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A/A/A)')                                       &
     &      ' XPRINP: INPUT ERROR following '//WORD,                    &
     &      ' First character is not blank; line in error:', LINE
      END IF
!
!     Line accepted; proceed
!     ======================
!
      ITKN = 2
!
!     XPRTKN returns initial adress(ITKN) and length (LTKN) of next token
!     on line. If ITKN -1, then no more tokens on LINE.
!
  200 CALL XPRTKN(LINE,ITKN,LTKN)
!
!     If no more tokens, read new line
!     ================================
!
      IF (LTKN .LE. 0) GO TO 100
!
!     ************************
!     *** Specify operator ***
!     ************************
!
!     Read property name
!     ==================
!
      IF (IKW .EQ. 0) THEN ! * Property name
!
         PNAME     = BLANK16
         LTKNR     = MIN(16,LTKN)
         READ(LINE(ITKN:80),'(A)') PNAME(1:LTKNR)
!        default diagonal operator name = label
         PLABEL(1) = PNAME(1:8)
         IPTYP     = 1
         NPCOMP    = 1
         NFAC      = 0
         IKW       = 1
      ELSE IF (IKW .EQ. 1) THEN  ! * Property operator keywords
         LKW = MIN(8,LTKN)
         READ(LINE(ITKN:80),'(A)') KEYWD(1:LKW)
         IF (LKW .LT. 3) THEN
            INPERR = INPERR + 1
             WRITE (LUPRI,'(/A/A)')                                     &
     &       ' XPRINP: INPUT ERROR following '//WORD,                   &
     &       ' Keyword "'//KEYWD(1:LKW)//                               &
     &       '" too short, min 3 char required.'
         END IF
!
!        Look for keyword from table
!        ===========================
!
         DO I = 1,NTABKW
         IF (KEYWD(1:LKW) .EQ. TABKW(I)(1:LKW)) THEN
!
!          If keyword defines operator type, then initialize
!          =================================================
!
           IF ( I .GT. I_FOCKMAT ) THEN
!             (the first fOUR entries in TABKW are not operator types)
!MI      ---- assignement of the IPTYP is needed in other routines -----
              IPTYP = I - I_FOCKMAT
              NLABEL= 0
              NFAC  = 0
              IKW = 2
           END IF
!
!          If keyword Fock defines matrix read, then initialize
!          ====================================================
!          (the matrix is assumed totally symmetric 
!           under point group and time reversal symmetry)
!
           IF (I .EQ. I_FOCKMAT) THEN
              IPTYP = 0
              NLABEL= 0
              NFAC  = 0
              IKW = 2
              ISYXPR = 1
              ITRXPR = 1
           END IF
 
           GO TO (4, &   ! CLASS
     &            5, &   ! FACTORS
     &            6, &   ! COMFACTOR
     &            7, &   ! CMULT
     &            1, &   ! FOCKMAT
     &            1, &   ! DIAGONAL
     &            1, &   ! XALPHA
     &            1, &   ! YALPHA
     &            1, &   ! ZALPHA
     &            2, &   ! XAVECTOR
     &            2, &   ! YAVECTOR
     &            2, &   ! ZAVECTOR
     &            3, &   ! ALPHADOT
     &            1, &   ! GAMMA5
     &            1, &   ! XSIGMA
     &            1, &   ! YSIGMA
     &            1, &   ! ZSIGMA
     &            1, &   ! XBETASIG
     &            1, &   ! YBETASIG
     &            1, &   ! ZBETASIG
     &            1, &   ! XiBETAAL
     &            1, &   ! YiBETAAL
     &            1, &   ! ZiBETAAL
     &            1, &   ! BETA
     &            3, &   ! SIGMADOT
     &            1), I ! iBETAGAM
         END IF
         END DO
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A/A)')                                         &
     &     ' XPRINP: INPUT ERROR following '//WORD,                     &
     &     ' Keyword "'//KEYWD(1:LKW)//'" not recognized.'
         GO TO 300
!&&&&    One component: DIAGONAL, [XYZ]ALPHA, GAMMA5, [XYZ]SIGMA,
!&&&&                   [XYZ]BETASIG, [XYZ]iBETAAL, BETA, FOCKMAT, iBETAGAMMA5
    1    CONTINUE
           NPCOMP = 1
         GO TO 300
!&&&&    Two components: [XYZ]AVECTOR
    2    CONTINUE
           NPCOMP = 2
!          note that the minus in the second component
!          is introduced later
         GO TO 300
!&&&&    Three components: ALPHADOT, SIGMADOT
    3    CONTINUE
           NPCOMP = 3
         GO TO 300
!&&&&    Class keyword, CLASS, TABKW(1)
    4    CONTINUE
           CLABEL = .TRUE.
           IKW = 3
         GO TO 300
!&&&&    Read a multiplicative factor for each component, FACTORS, TABKW(2)
    5    CONTINUE
           IKW  = 4
         GO TO 300
!&&&&    Read a common factor for each component, COMFACTOR, TABKW(3)
    6    CONTINUE
           IKW  = 5
         GO TO 300
!&&&     Multiply with current value of the speed of light, CMULT
    7    CONTINUE
           PFAC(1) = PFAC(1)*CVAL
           PFAC(2) = PFAC(2)*CVAL
           PFAC(3) = PFAC(3)*CVAL
         GO TO 300
!&&&     
    8    CONTINUE
         GO TO 300
         !GO TO 100
!
!     Read property labels for each component
!     =======================================
!
      ELSE IF (IKW .EQ. 2) THEN ! * Property labels (to be found on file)
         NLABEL = NLABEL + 1
         PLABEL(NLABEL) = BLANK8
         LTKNR = MIN(8,LTKN)
         READ(LINE(ITKN:80),'(A)') PLABEL(NLABEL)(1:LTKNR)
         IF (NLABEL .GE. NPCOMP) IKW = 1
!
!     Read class keyword
!     ==================
!
      ELSE IF (IKW .EQ. 3) THEN ! * Class label
!HJMAERKE: implement interpretation of class keyword
!miro - missing field, HJJ has an idea...
!
!     FACTORS, read factor for each component
!     ==================================================
!
      ELSE IF (IKW .EQ. 4) THEN ! * Separate factor for each component
         NFAC = NFAC + 1
         JTKN = ITKN - 1 + LTKN
         I    = 16 - LTKN + 1
         TOKEN16 = BLANK16
         TOKEN16(I:16) = LINE(ITKN:JTKN)
         READ(TOKEN16,'(F16.0)') TMP
         PFAC(NFAC) = TMP*PFAC(NFAC)
         IF (NFAC .GE. NPCOMP) IKW = 1
!
!     COMFACTOR, read a common factor for all components
!     =======================================
      ELSE IF (IKW .EQ. 5) THEN ! * Common factor for all components
         JTKN = ITKN - 1 + LTKN
         I    = 16 - LTKN + 1
         TOKEN16 = BLANK16
         TOKEN16(I:16) = LINE(ITKN:JTKN)
         READ(TOKEN16,'(F16.0)') COMFAC
         PFAC(1) = COMFAC*PFAC(1)
         PFAC(2) = COMFAC*PFAC(2)
         PFAC(3) = COMFAC*PFAC(3)
         IKW = 1

!
!     Unknown IKW
!     ===========
!
      ELSE
!        IKW .gt. 5
         INPERR = INPERR + 1
         WRITE (LUPRI,*)                                                &
     &    'Programming error in XPRINP, unknown IKW=',IKW
      END IF
!
  300 ITKN = ITKN + LTKN
      GO TO 200
!
!**********************************************************
!***** Look up in list of operators already defined. ******
!**********************************************************
!
!     Find INDXPR, ISYXPR, ITRXPR
!
  500 CONTINUE
      IF (CLABEL) THEN
!HJMAERKE: implement interpretation of class keyword
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A/A)')                                         &
     &      ' XPRINP: INPUT ERROR following '//WORD,                    &
     &      ' Keyword "'//KEYWD(1:LKW)//'" not implemented yet.'
      ELSE
      ! miro-output and XPRIND
         IF (IPRINT.GE.5) THEN
           WRITE(LUPRI,'(/,6x,a)') '*** Output from XPRINP ***'
           WRITE(LUPRI,'(2x,a,a16)') 'PNAME=',PNAME
           WRITE(LUPRI,'(2x,a,i2,a,i2)')                                &
     &     'IPTYP=',IPTYP,' NPCOMP=',NPCOMP
           WRITE(LUPRI,'(2x,a,3a8)')                                    &
     &     'PLABEL(1-NPCOMP)=',(PLABEL(I),I=1,NPCOMP)
           WRITE(LUPRI,'(2x,a,3f10.6)')                                 &
     &     'PFAC(1-NPCOMP)=',(PFAC(I),I=1,NPCOMP)
           WRITE(LUPRI,'(2x,a,i2,a,i2,a,i2)')                           &
     &     'initial values: INDXPR=',INDXPR,                            &
     &    ' ISYXPR=',ISYXPR,' ITRXPR=',ITRXPR
         ENDIF
!MI      ... assign remaining parameters in dedicated subroutine
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,                    &
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)

         IF (IPRINT.GE.8) THEN
           WRITE(LUPRI,'(/,6x,a)')                                      &
     &      '*** Output from XPRINP after XPRIND ***'
           WRITE(LUPRI,'(2x,a,i2,a,i2,a,i2)')                           &
     &     'INDXPR=',INDXPR,' ISYXPR=',ISYXPR,' ITRXPR=',ITRXPR
           WRITE(LUPRI,'(2x,a,i2,a,i2)')                                &
     &     'IPTYP=',IPTYP,' NPCOMP=',NPCOMP
           WRITE(LUPRI,*)
         ENDIF
      END IF
      GO TO 9999
!
 8000 INPERR = INPERR + 1
      WORD = ' INPEOF'
      WRITE (LUPRI,'(/A/A)')                                            &
     &' XPRINP: INPUT ERROR following '//WORD,                          &
     &' End of file encountered; input must finish with "**something".'
!
 9999 CALL QEXIT ('XPRINP')
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck xprtkn */
      SUBROUTINE XPRTKN(LINE,ITKN,LTKN)
!
!  Copyright (c) 27-Jan-1995 Hans Joergen Aa. Jensen
!
!  Find initial address (ITKN) and length (ITKN) of next token
!  If no more tokens return ITKN = -1
!
      PARAMETER (MAXLL = 80)
      CHARACTER*(MAXLL) LINE
      CHARACTER DELIM*2, DELIM2*1
      PARAMETER (DELIM = '''"')
!
      CALL QENTER('XPRTKN')
!
!     Find first character in next token
!
      DO I = ITKN,MAXLL
         IF (LINE(I:I) .NE. ' ') THEN
            IF (INDEX(DELIM,LINE(I:I)) .GT. 0) THEN
               DELIM2 = LINE(I:I)
               IND1 = I + 1
            ELSE
               DELIM2 = ' '
               IND1 = I
            END IF
            GO TO 100
         END IF
      END DO
      ITKN = -1
      LTKN = 0
      GO TO 9999
!
!     Find last character of this token
!
  100 CONTINUE
      INDL = MAXLL
      DO I = IND1,MAXLL
         IF (LINE(I:I) .EQ. DELIM2) THEN
            INDL = I - 1
            LINE(I:I) = ' '
            GO TO 200
         END IF
      END DO
  200 ITKN = IND1
      LTKN = (INDL + 1 - IND1)
!
 9999 CALL QEXIT('XPRTKN')
      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck xprind */
      SUBROUTINE XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,                 &
     &           INDXPR,ISYXPR,ITRXPR,IPRINT)
!*********************************************************************** 
!
!  Given a 4-component operator this routine returns an index
!  INDXPR to the list of operators in COMMON block XCBPRP and
!  returns its symmetry under time reversal(ITRXPR) and its
!  spatial symmetry(ISYXPR).
!
!  This routine is called from many places 
!
!  Copyright (c) 27-Jan-1995 Hans Joergen Aa. Jensen
!  Revised May 14 1996 - tsaue
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER PNAME*16, PLABEL(NPCOMP)*8
      DIMENSION PFAC(NPCOMP)
!
! Used from common blocks:
!  XCBPRP: everything
!  XCBPRL: PRPLBL()
!  XCBPRI: LUPRI,NINFO
!
#include "dcbxpr.h"
#include "dcbprl.h"
!
!     CALL QENTER('XPRIND')
      INDXPR = -1

      IF (IPRINT.GE.15) THEN
        WRITE(LUPRI,'(/,4x,a)') '***** Output from XPRIND *****'
        WRITE(LUPRI,'(2x,a,a16)') 'Entering PNAME=',PNAME
        WRITE(LUPRI,'(2x,a,i2)') 'Entering NPCOMP=',NPCOMP
        WRITE(LUPRI,'(2x,a,6a8)')                                       &
     &   'Entering PLABEL=',(PLABEL(I),I=1,NPCOMP)
        WRITE(LUPRI,'(2x,a,6f12.6)')                                    &
     &    'Entering PFAC=',(PFAC(I),I=1,NPCOMP)
        WRITE(LUPRI,'(2x,a,i2,a,i3)')                                   &
     &   'Entering IPTYP=',IPTYP,' INDXPR=',INDXPR
      ENDIF

!================================================================
!       Search for this property operator in existing list
!================================================================
!TF   print*,'Number of property operators : ',NPRPS
      DO 200 I = 1,NPRPS
!        Check operator type
         IF (IPRPTYP(I) .NE. IPTYP) GO TO 200
!        Check factors and property labels for each component
         DO 100 J = 1,NPCOMP
            IF (PFAC(J) .NE. FACPRP(J,I)) GO TO 200
            K = IPRPLBL(J,I)
            IF (PRPLBL(K) .NE. PLABEL(J)) GO TO 200
  100    CONTINUE
!
         IF (PNAME .NE. PRPNAM(I)) THEN
!        Operator exists, but under another name; rename
            WRITE (LUPRI,'(3(/A))')                                     &
     &      ' XPRIND INFO: property operator named '//PNAME,            &
     &      '   is already existing under the name '//PRPNAM(I),        &
     &      '   previous operator renamed to       '//PNAME
            PRPNAM(I) = PNAME
         END IF
         INDXPR = I
         GO TO 1000
  200 CONTINUE

!========================================================================
!     This property operator has not been defined previously.
!     Append the new property operator to /XCBPRP/ list.
!=========================================================================
!TF   print*,' Calling routine XPRADD, new operator.'
      CALL XPRADD(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,INDXPR,IPRINT)
!TF   print*,' INDXPR after XPRADD: ',INDXPR
!
 1000 CONTINUE
      IF(IPTYP.EQ.0) THEN
!       ... i.e. FOCKMAT 
        IPRPSYM(INDXPR) = ISYXPR
        IPRPTIM(INDXPR) = ITRXPR
      ELSE
        ISYXPR = IPRPSYM(INDXPR)
        ITRXPR = IPRPTIM(INDXPR)
      ENDIF
!TF   print*,' INDXPR,ISYXPR,ITRXPR :',INDXPR,ISYXPR,ITRXPR
!TF   print*,' List IPRPTIM : '
!TF   do IL = 1,22,1
!TF     print*,' IPRPTIM(',IL,') = ',IPRPTIM(INDXPR)
!TF   end do
!     CALL QEXIT('XPRIND')
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck xpradd */
      SUBROUTINE XPRADD(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,INDXPR,IPRINT)
!***********************************************************************
!
!  Add property to list of properties in COMMON block XCBPRP
!  Copyright (c) 6-Feb-1995 Hans Joergen Aa. Jensen
!
!  Entering: PNAME
!            IPTYP
!            NPCOMP
!
!  Routine is called only from XPRIND.
!
!  MI (&HAaaJ), 2003,2004: Extended for LAO integral types
!               as well as other integrals
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER PNAME*16, PLABEL(NPCOMP)*8
      DIMENSION PFAC(NPCOMP)
!
! Used from common blocks:
!  XCBPRP: everything
!  XCBPRI: LUPRI,NINFO
!
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbham.h"
!MI ... needs CVAL
#include "dcbgen.h"

      CHARACTER*4 PCOMB, PCOMB_SAVE
!     ... the four characters in PCOMB define inclusion of the
!         LL, SL, LS, and SS blocks of one-electron operators.
!
      CALL QENTER('XPRADD')
!
!     This property operator has not been defined previously.
!     Append the new property operator to /XCBPRP/ list.

      IF (IPRINT .GE. 6) THEN
         WRITE (LUPRI,'(/,4X,A/2A/A,I5/A,I5/A,1P,3(D20.10))')           &
     &   '*** XPRADD called for specifying PCOMB ***',                  &
     &   '  PNAME   =',PNAME,                                           &
     &   '  IPTYP   =',IPTYP,                                           &
     &   '  NPCOMP  =',NPCOMP,                                          &
     &   '  PFAC    =',(PFAC(I),I=1,NPCOMP)
         WRITE (LUPRI,'(A,3(A8,2X))')                                   &
     &   '  PLABEL  =',(PLABEL(I),I=1,NPCOMP)
      END IF
!MI  List of exceptions:
!MI/Dec.2005 I guess this routine could be reprogrammed as exceptions are before and after cycle.
      IF ((PLABEL(1)(1:7).EQ.'BETAMAT').OR.(PNAME(6:10).EQ.'RM1H1')     &
     &    .OR.(PNAME(7:11).EQ.'RM2H1')) THEN
!MI (&TEC): RM1H1/RM2H1 is first/second London derivative of BETAMAT
         PCOMB = '000+'
!MI/HJAaJ
! ... NEWS for LAO-ccRDSUSLL in Levy-Leblond - cancel the SL(2,1) block
!      to get compatibility with DALTON, plus doubled factor
!  nonrelativistic adaptation has to be set up explicitly 
      ELSE IF (PNAME(7:13).EQ.'RDSUSLL'.AND.RDSUSLLMOD) THEN
! ... fix the LS block for the LevyLe option
         PCOMB = '00+0'
        IF (IPRINT.GE.6) WRITE(LUPRI,'(/A)')                            &
     & 'XPRADD: Only the LS block kept for the RDSUSLL operator!'
      ELSE IF ((PNAME(1:5).EQ.'d|S>/'.OR.PNAME(1:3).EQ.'dS/').AND.BSS)  &
     &     THEN
        PCOMB='+000'
        IF (IPRINT.GE.6) WRITE(LUPRI,'(/A)')                            &
     &  'XPRADD: Only the LL block for connection'//                    &
     &  ' matrixes in two-component mode!'
      ELSE
         IF (IPTYP.EQ.0) THEN
!          ... for FOCKMAT
           PCOMB = 'xxxx'
           GOTO 200           
         ELSEIF (IPTYP.LT.0 .OR. IPTYP.GT.21) THEN
           WRITE(LUPRI,*)  'IPTYP out of the range, IPTYP=',IPTYP
           CALL QUIT('XPRADD: Illegal IPTYP (not 1-22)!')
         ENDIF
         GO TO (1,2,2,2,2,2,2,2,2,                                      &
     &          1,1,1,3,3,3,4,4,4,3,1,4),IPTYP
!
!
!        Operator types: DIAGONAL, SIGMA
!
    1    PCOMB = '+00+'
!JTH DEBUG!!
#ifdef UNDEF
         write(6,*) 'debug debug: operator ',PNAME
         IF (PNAME(1:2) .EQ. 'EF') THEN
            PCOMB = '+000'
            write(6,*) '           : setting comp to',pcomb
         END IF
#endif
         GO TO 100
!
!        Operator types: ALPHA, GAMMA5
!
    2    PCOMB = '0++0'
         GO TO 100
!
!        Operator types: BETA, BETA SIGMA
!
    3    PCOMB = '+00-'
         IF (LEVYLE) PCOMB = '+000'
         GO TO 100
!
!        Operator types: i BETA ALPHA, i BETA GAMMA5
!
    4    PCOMB = '0+-0'
         GO TO 100
!
  100    CONTINUE

!MI  ... for the Levy-Leblond mode cancel the SS block!
!MI ... check first if PLABEL(1) is not operator from
!MI ... pilot Dirac bare nucleus (from subroutine OP1INI)
        IF (LEVYLE.AND.PLABEL(1).NE.'OVERLAP '                          &
     &       .AND.PLABEL(1).NE.'MOLFIELD'.AND.                          &
     &       PLABEL(1).NE.'BETAMAT '                                    &
     &       .and.PLABEL(1).NE.'FDEVEMB ' ) PCOMB(4:4) = '0'
!LV ... same thing for ZORA!
        IF ((ZORA.AND..NOT.ZORA4).AND.PLABEL(1).NE.'OVERLAP '           &
     &       .AND.PLABEL(1).NE.'MOLFIELD'.AND.                          &
     &       PLABEL(1).NE.'BETAMAT '                                    &
     &       .and.PLABEL(1).NE.'FDEVEMB ' ) PCOMB(4:4) = '0'
      END IF

!MI ... inform about the assigned PCOMB
      IF (IPRINT .GE. 5) THEN
         WRITE(LUPRI,'(1X,A,A4)')                                       &
     &   'for this operator assigned PCOMB:',PCOMB
      END IF
!     ...  assign proper factor for EDM  ... accordingto E.D.Commins
      IF (PLABEL(1).EQ.'EDM') THEN
        PFAC(1)=+4.0D0 * CVAL * PFAC(1)
        WRITE(LUPRI,'(2X,A,D12.6)')                                     &
     &  '4c EDM2 operator: strength multiplied by +4*CVAL,'             &
     & //' PFAC(1)=', PFAC(1)
      ENDIF

! miro 7 June 2012 exception: for PLABEL(1:5)="e_EDM" assign PCOMB(1:1)='0'
      IF (PNAME(1:5).EQ.'e_EDM') THEN
         !todo: hard assign 000+ structure
         IF (PCOMB.NE.'+00+') THEN 
           PCOMB = '+00+' ! assig as for SIGMA type of operator
           WRITE(LUPRI,*) 'reassigned PCOMB=',PCOMB
         ENDIF
         PCOMB_SAVE = PCOMB
         PCOMB(1:1)='0'
         WRITE(LUPRI,'(2X,A,A8,A,A4,A,A4)') "Operator ",PNAME,          &
     &   'changed 4-component matrix structure: ',                      &
     &   PCOMB_SAVE," --> ",PCOMB
         ! adjust the operator strength
         DO I=1,NPCOMP
           PFAC_TEMP = PFAC(I)
           PFAC(I)=2.0d0*PFAC_TEMP
           write(LUPRI,'(3X,I1,A,D7.2,A,D7.2)')                         &
     &     I,'.factor of ',PFAC_TEMP,                                   &
     &     ' multiplied by 2, giving ', PFAC(I)
         ENDDO
      ENDIF

!
!     Define property
!
 200  CONTINUE
      NPRPS  = NPRPS + 1
      IF(NPRPS.GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)')                                           &
     &   ' ** ERROR in XPRADD **  Redimension MAXPRPS to ',NPRPS
        CALL QUIT('XPRADD: Too may properties !')
      ENDIF
!MI   ... assign key variables ...
      INDXPR = NPRPS
!     Property name
      PRPNAM(INDXPR)  = PNAME
!     Property type
      IPRPTYP(INDXPR) = IPTYP
!     Pointers to property labels
      IF(IPTYP.EQ.0) THEN
        DO 300 I = 1,NPRPLBL
!         Check property label
          IF (PLABEL(1).NE.PRPLBL(I)) GOTO 300
          IF (PCOMB.NE.PDOINT(I))  GOTO 300
          INDX = I
        GO TO 400
  300   CONTINUE
        NPRPLBL = NPRPLBL + 1
        INDX    = NPRPLBL
  400   CONTINUE      
        FACPRP(1,INDXPR)  = PFAC(1)
        IPRPLBL(1,INDXPR) = INDX
        PRPLBL(INDX) = PLABEL(1)
        PDOINT(INDX) = PCOMB
      ELSE
!MI     ... (IPTYP.NE.0) operator has more components
        DO 600 J = 1,NPCOMP
          FACPRP(J,INDXPR)  = PFAC(J)
          INDX = INDXPL(PLABEL(J),PCOMB,IPRINT)
          IF (INDX .LE. 0) THEN
            INDXPR = -1
            WRITE (LUPRI,'(3A)')                                        &
     &      'PLABEL, PCOMB not found in /XCBPRL/ : ', PLABEL(J), PCOMB
          END IF
          IPRPLBL(J,NPRPS) = INDX
  600   CONTINUE
        IF (INDXPR .EQ. -1) CALL QUIT('XPRADD: PLABEL, PCOMB not found')
!       Symmetry under time reversal and point group symmetry
        CALL PRPSYM(INDXPR,IPRINT)
      ENDIF
!
      CALL QEXIT('XPRADD')
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck indxpl */
      INTEGER FUNCTION INDXPL(PLABEL,PCOMB,IPRINT)
!***********************************************************************
!
!     9-Mar-1995 hjaaj
!     Complete revision : May 14 1996 - tsaue
!
!     Return value is pointer to list of property labels in 
!     COMMON block PRPLBL. The routine will furthermore update
!     list of integrals to calculate found in COMMON block
!     PRPCLB
!
!     If INDXPL .le. 0 then no match was found.
!
!     Last revision: May 15 1996 - tsaue
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER PLABEL*8, PCOMB*4,CLASLB*7
!
! Used from common blocks:
!   XCBPRL : PRPFILE
!
#include "dcbprl.h"
#include "dcbcls.h"
!
      CALL QENTER('INDXPL')
!
      INDXPL = -1
      DO 100 I = 1,NPRPLBL
!       Check property label
        IF (PLABEL.NE.PRPLBL(I)) GOTO 100
        IF (PCOMB.NE.PDOINT(I))  GOTO 100
!
!       Label and block specification already defined; return
!
        INDXPL = I
        GO TO 9999
  100 CONTINUE
!
!     Class label not yet defined; first check for class label
!
      IF(PLABEL.EQ.'MOLFIELD') THEN
        CLASLB = 'MOLFLD '
        IORDER  = 0
      ELSEIF(PLABEL(1:7).EQ.'FDEVEMB') THEN
        CLASLB = 'FDEVEMB'
        IORDER   = 0
      ELSEIF(PLABEL(1:7).EQ.'BETAMAT') THEN
        CLASLB = 'BETAMAT'
        IORDER   = 0
      ELSEIF(PLABEL(1:3).EQ.'EDM') THEN
        CLASLB = 'EDM    '
        IORDER   = 0
      ELSEIF(PLABEL(1:2).EQ.'ED') THEN
        CLASLB = 'EFFDEN '
        IORDER   = 0
      ELSE
        CALL GETCLS(PLABEL,CLASLB,IORDER,IPRINT)
        IF(CLASLB.EQ.'Unknown') THEN
          WRITE(LUPRI,'(A,A8)')                                         &
     &    '***** INDXPL ERROR ***** No class label found for ',         &
     &    PLABEL
          CALL QUIT('INDXPL ERROR: No class label found !')
        ENDIF
      ENDIF
      NPRPCLS = NPRPCLS + 1
      IF (NPRPCLS .GT. MAXCLS) THEN
         WRITE(LUPRI,'(/A/A,I5,A/A,1X,A)')                              &
     &      'INDXPL ERROR: MAXCLS in dcbcls.h too small',               &
     &      'Limit of',MAXCLS,' exceeded for operator:',                &
     &      PLABEL, PCOMB
         CALL QUIT('INDPXL: MAXCLS in dcbcls.h too small')
      END IF
      CLSINT(NPRPCLS) = CLASLB
      CLSCMB(NPRPCLS) = PCOMB
      CLSCAL(NPRPCLS) = .TRUE.
      IORDCL(NPRPCLS) = IORDER
      CALL XPLDEF(PLABEL,CLASLB,IORDER,PCOMB,INDX,IPRINT)
      INDXPL = INDX
 9999 CONTINUE
      CALL QEXIT('INDXPL')
      RETURN
!
 1000 CONTINUE
 1010 CONTINUE
!radovan: rest of the routine is unreachable
      WRITE (LUPRI,'(/A/2A/2A)')                                        &
     &'INDXPL ERROR: illegal PCOMB or PDOINT',                          &
     &'    PCOMB    = ',PCOMB,                                          &
     &'    PDOINT   = ',PDOINT(I)
      STOP 'INDXPL ERROR: illegal PCOMB or PDOINT(I)'
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck getcls */
      SUBROUTINE GETCLS(LABINT,CLASLB,IORDER,IPRINT)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!     Given a file label for a given integral, this subroutine
!     returns the label for the general class of integrals
!     the integral belongs to. For moment integrals the order is returned
!     as well.
!
!     Jan-18-1995 - tsaue
!
!  2003, MI,HJAaJ extended list for integrals over LAO
!        for NMR shielding and magnetic susceptibilities
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      CHARACTER LABINT*8,CLASLB*7

      CALL QENTER('GETCLS')

!     Overlap integral
!
      IF(LABINT.EQ.'OVERLAP') THEN
        CLASLB = 'OVERLAP'
!
!     Dipole length/velocity integrals
!
      ELSE IF(LABINT(2:4).EQ.'DIP') THEN
        IF(LABINT(5:7).EQ.'VEL') THEN
          CLASLB = 'DIPVEL '
        ELSE
          CLASLB = 'DIPLEN '
        ENDIF
!
!     Quadrupole moments integrals
!
      ELSE IF(LABINT(3:8).EQ.'QUADRU') THEN
        CLASLB = 'QUADRUP'
!
!     Spin-orbit integrals
!
      ELSE IF(LABINT(3:8).EQ.'SPNORB') THEN
        CLASLB = 'SPNORB '
!          
!     SOPP integrals for ECP
!          
      ELSE IF(LABINT(3:6) .EQ.'SOPP') THEN
        CLASLB = 'SOPPINT'
!
!     Second moments integrals
!
      ELSE IF(LABINT(3:8).EQ.'SECMOM') THEN
        CLASLB = 'SECMOM '
!
!     Traceless theta quadrupole integrals
!
      ELSE IF(LABINT(3:7).EQ.'THETA') THEN
        CLASLB = 'THETA  '

!requested by Malaya Kumar Nayak <mk.nayak72@gmail.com>
      ELSE IF(LABINT(2:8).EQ.'NUCFIEL') THEN
        CLASLB = 'NUCFIEL'

!
!     Cartesian/spherical/electronic solvent moments integrals
!
      ELSE IF(LABINT(2:2).EQ.'M') THEN
        IF(LABINT(1:1).EQ.'C') THEN
          CLASLB = 'CARMOM '    
          READ(LABINT(3:4),'(I2)') NX
          READ(LABINT(5:6),'(I2)') NY
          READ(LABINT(7:8),'(I2)') NZ
          IORDER = NX + NY + NZ
        ELSE
          CLASLB = 'SPHMOM '
          READ(LABINT(3:4),'(I2)') IORDER
        ENDIF
!
!     One-electron Fermi contact integrals
!
      ELSE IF(LABINT(1:2).EQ.'FC') THEN
        CLASLB = 'FERMI C'
!
!     Paramagnetic spin-orbit integrals
!
      ELSE IF(LABINT(1:3).EQ.'PSO') THEN
        CLASLB = 'PSO    '
!     
!     Spin-dipole integrals
!
      ELSE IF(LABINT(1:2).EQ.'SD') THEN
!
!     Spin-dipole + Fermi contact integrals
!
        IF(LABINT(3:3).EQ.'C') THEN
          CLASLB = 'SDFC   '
        ELSE
!     
!     Spin-dipole integrals
!
          CLASLB = 'SPIN-DI'
        ENDIF
!
!     Diamagnetic spin-orbit integrals
!
      ELSE IF(LABINT(1:2).EQ.'DS') THEN
        CLASLB = 'DSO    '
!
!     Half-derivative overlap integrals
!
      ELSE IF(LABINT(1:3).EQ.'HDO') THEN
        CLASLB = 'HDO   '
!
!     Contribution from overlap matrix to magnetic properties
!
      ELSE IF(LABINT(1:5).EQ.'dS/dB') THEN
        IF(LABINT(6:6).EQ.'2') THEN
!
!       Contribution from overlap matrix to magnetic properties
!
          CLASLB = 'S2MAG  '
        ELSE
!
!       Second order contribution from overlap matrix to 
!       magnetic properties
!
          CLASLB = 'S1MAG  '
        ENDIF
!
!     Angular momentum around the nuclei
!
      ELSE IF(LABINT(2:7).EQ.'ANGLON') THEN
        CLASLB = 'ANGLON '
!
!     Electronic angular momentum around origin
!
      ELSE IF(LABINT(2:7).EQ.'ANGMOM') THEN
        CLASLB = 'ANGMOM '
!
!     London orbital contribution to angular momentum
!
      ELSE IF(LABINT(2:7).EQ.'LONMOM') THEN
        CLASLB = 'LONMOM '
!
!     One-electron contributions to magnetic moment
!
      ELSE IF(LABINT(1:5).EQ.'dh/dB') THEN
        CLASLB = 'MAGMOM '
!
!     Electronic kinetic energy
!
      ELSE IF(LABINT(1:7).EQ.'KINENER') THEN   
        CLASLB = 'KINENER'
!
!     Diamagnetic susceptibility
!
      ELSE IF(LABINT(4:6).EQ.'SUS') THEN
        IF    (LABINT(7:8).EQ.'NL') THEN
!
!         Diamagnetic susceptibility without London contribution
!
          CLASLB = 'DSUSNOL'
        ELSEIF(LABINT(7:8).EQ.'LL') THEN
!
!         Angular London orbital contribution to diamagnetic
!         susceptibility
!
          CLASLB = 'DSUSLAN'
        ELSE
!
!       Angular London orbital contribution to diamagnetic 
!       susceptibility
!
          CLASLB = 'DSUSLH'
        ENDIF
!
!     Angular London orbital contribution to diamagnetic susceptibility
!
      ELSE IF(LABINT(3:8).EQ.'dh/dB2') THEN
        CLASLB = 'DIASUS '
!
!     Nuclear shielding integrals without London orbital contribution
!
      ELSE IF(LABINT(4:7).EQ.'NSNL') THEN
        CLASLB = 'NUCSNLO'
!
!     London orbital contribution to nuclear shielding tensor integrals
!
      ELSE IF(LABINT(4:7).EQ.'NSLO') THEN
        CLASLB = 'NUCSLO '
!
!     Nuclear shielding tensor integrals
!MI ... use this rather than 'NST'
      ELSE IF(LABINT(4:7).EQ.' NST') THEN
        CLASLB = 'NUCSHI '
!
!     Electric field at the individual nuclei
!
      ELSE IF(LABINT(1:3).EQ.'NEF') THEN
        CLASLB = 'NEFIELD'
!
!     traceless electric field third derivatives
!
      ELSE IF(LABINT(5:5).EQ.'V') THEN
        CLASLB = 'EFT    '
!
!     Parity violation - chirality operator
!
      ELSE IF(LABINT(1:3).EQ.'PVC') THEN
        CLASLB = 'PVC    '
!
!     Electric field gradient at the individual nuclei, spherical
!
      ELSE IF(LABINT(3:5).EQ.'FGS') THEN
        CLASLB = 'ELFGRDS'
!
!     Bra-differentiation of overlap matrix with respect to magnetic 
!     field
!
      ELSEIF(LABINT(1:7).EQ.'d<S|/dB') THEN
        CLASLB = 'S1MAGL'
!
!     Ket-differentiation of overlap matrix with respect to 
!     magnetic field
!
      ELSE IF(LABINT(1:7).EQ.'d|S>/dB') THEN
        CLASLB = 'S1MAGR'
!
!     Ket-differentation of HDO-integrals with respect to magnetic field
!
      ELSE IF(LABINT(4:6).EQ.'HDB') THEN
        CLASLB = 'HDOBR  '
!
!     Potential energy at the nuclei
!
      ELSE IF(LABINT(1:5).EQ.'POT.E') THEN
        CLASLB = 'NUCPOT '
!
!     Half B-differentiated overlap matrix
!
      ELSE IF(LABINT(2:5).EQ.'HBDO') THEN
        CLASLB = 'HBDO   '
!
!     Half-derivative overlap integrals not to be antisymmetrized
!
      ELSE IF(LABINT(1:5).EQ.'SQHDO') THEN
        CLASLB = 'SQHDO  '
!
!     Diamagnetic susceptibility with common gauge origin
!
      ELSE IF(LABINT(3:8).EQ.'SUSCGO') THEN
        CLASLB = 'DSUSCGO'
!
!     Nuclear shielding integrals with common gauge origin
!
      ELSE IF(LABINT(4:7).EQ.'NSCO') THEN
        CLASLB='NSTCGO '
!
!     Cosine and sine integrals
!
      ELSE IF(LABINT(2:7).EQ.'COSINE'.OR.LABINT(2:5).EQ.'SINE') THEN
        CLASLB = 'EXPIKR '
!
!     Mass velocity integrals
!
      ELSE IF(LABINT.EQ.'MASSVELO') THEN
        CLASLB = 'MASSVEL'
!
!     Darwin type integrals
!
      ELSE IF(LABINT(1:6).EQ.'DARWIN') THEN
        CLASLB = 'DARWIN '
!
!     First order magnetic field derivatives of electric field
!
      ELSE IF(LABINT(3:5).EQ.'CM1') THEN
        CLASLB = 'CM1    '
!
!     Second order magnetic field derivatives of electric field
!
      ELSE IF(LABINT(3:5).EQ.'CM2') THEN
        CLASLB = 'CM2    '
!
!     Half-derivative overlap integrals not to be anti-symmetrized
!
      ELSE IF(LABINT(1:6).EQ.'SQHDOR') THEN
        CLASLB = 'SQHDOR '
!
!
!
      ELSE IF(LABINT.EQ.'SQOVLAP') THEN
        CLASLB = 'SQOVLAP'
!
!     non-traceless electric field third derivatives
!
      ELSE IF(LABINT(5:5).EQ.'v') THEN
        CLASLB = 'EFTNTL '
!
!     Electric field gradient at the individual nuclei, cartesian
!
      ELSE IF(LABINT(3:5).EQ.'EFG') THEN
        CLASLB = 'ELFGRDC'
      ELSE IF(LABINT(2:6).EQ.'RM1H2') THEN
        CLASLB = 'RM1H2'
      ELSE IF(LABINT(3:7).EQ.'RM2H2') THEN
        CLASLB = 'RM2H2'
      ELSE IF(LABINT(2:6).EQ.'RM1H3') THEN
        CLASLB = 'RM1H3'
      ELSE IF(LABINT(3:7).EQ.'RM2H3') THEN
        CLASLB = 'RM2H3'
      ELSE IF(LABINT(3:7).EQ.'RDSUL') THEN
        CLASLB = 'RDSUSLL'
      ELSE IF(LABINT(2:6).EQ.'RM1RN') THEN
        CLASLB = 'RM1RN'
!MI  London orbital contribution to exp.value of nuclear shielding tensor integrals
      ELSE IF(LABINT(4:7).EQ.'RNST') THEN
        CLASLB = 'RM1N1H '
!MI  CAP integrals/labels to be used in the MOLTRA/PRPTRA input
      ELSE IF(LABINT(1:6).EQ.'CAP_RE') THEN
        CLASLB = 'CAP_RE  '
      ELSE IF(LABINT(1:6).EQ.'CAP_IM') THEN
        CLASLB = 'CAP_IM  '
      ELSE IF(LABINT(1:6).EQ.'CAP_OV') THEN
        CLASLB = 'CAP_OVL '
      ELSE IF(LABINT(2:7).EQ.'CAPD1R') THEN
        CLASLB = 'CAPD1R' 
      ELSE IF(LABINT(2:7).EQ.'CAPD1I') THEN
        CLASLB = 'CAPD1I' 
      ELSE IF(LABINT(2:7).EQ.'CAPDVE') THEN
        CLASLB = 'CAPDVE' 
      ELSE IF(LABINT(1:6).EQ.'CAPD_V') THEN
        CLASLB = 'CAPD_VD' 
!MI  EDM2 operator, see Refs. [Nataraj2007]_, [Mukherjee2009]_
      ELSE IF(LABINT(1:3).EQ.'EDM') THEN
        CLASLB = 'KINENER'
!    original  e_EDM operator of Flambaum
      ELSE IF(LABINT(2:6).EQ.'E_EDM') THEN
        CLASLB = 'E_EDM'
!    complex exponential operator
      ELSE IF(LABINT(1:4).EQ.'CEXP') THEN
        CLASLB = 'CXIKR'
!
!     Unknown LABINT
!   
      ELSE
       CLASLB = 'Unknown'
      ENDIF

!  Control print out
      IF (IPRINT.GE.5) THEN
        WRITE(LUPRI,'(/,2A)')                                           &
     &  'GETCLS called for LABINT=',LABINT
        WRITE(LUPRI,'(7X,2A)')                                          &
     &  ' assigned  CLASLB=',CLASLB
      ENDIF
!
      CALL QEXIT('GETCLS')
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck xpldef */
      SUBROUTINE XPLDEF(PLABEL,CLASLB,IORDER,PCOMB,INDX,IPRINT)
!***********************************************************************
!
!     Define property label: symmetries etc..
!
!     Called from: integ.function INDXPL
!
!     Written by T.Saue May 14 1996
!     Last revision: May 14 1996 - tsaue
!
!***********************************************************************
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
!
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
#include "dcbprl.h"
#include "cbiher.h"
      CHARACTER PLABEL*8,CLASLB*7,PCOMB*4
      LOGICAL ANTI, TRASPH, SOLVNT, SQUARE, LABERR
      integer, allocatable :: INTREP(:)
      integer, allocatable :: INTADR(:)
      logical, allocatable :: DOATOM(:)
      character (len=8), allocatable :: LABINT(:)
!
      CALL QENTER('XPLDEF')
!
      MAXTYP = 15*MXCENT
      IF ((CLASLB.EQ.'ELFGRDC').OR.(CLASLB.EQ.'ELFGRDS')) THEN
         NINTAD = 9*NUCIND*(MAXREP+1)
      ELSE
         NINTAD = MAXTYP
      END IF
!.....Memory allocation
      call alloc(INTREP,MAXTYP,id='INTREP in XPLDEF')
      call alloc(INTADR,NINTAD,id='INTADR in XPLDEF')
      allocate(DOATOM(NUCIND))
      allocate(LABINT(MAXTYP))
!
!     Get information on this class of integrals
!
      IF    (CLASLB.EQ.'MOLFLD ') THEN
        NOPTYP    = 1
        LABINT(1) = 'MOLFIELD'
        INTREP(1) = 0
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        IORDER    = 0
        DO I = 1,NUCIND
          DOATOM(I) = .TRUE.
        ENDDO
      ELSEIF(CLASLB.EQ.'FDEVEMB') THEN
        NOPTYP    = 1
        LABINT(1) = 'FDEVEMB '
        INTREP(1) = 0
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        IORDER    = 0
        DO I = 1,NUCIND
          DOATOM(I) = .TRUE.
        ENDDO
      ELSEIF(CLASLB.EQ.'BETAMAT') THEN
        NOPTYP    = 1
        LABINT(1) = 'BETAMAT '
        INTREP(1) = 0
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        IORDER    = 0
        DO I = 1,NUCIND
          DOATOM(I) = .TRUE.
        ENDDO
      ELSEIF(CLASLB .EQ. 'SPHMOM ') THEN
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        CALL SMTYP(IORDER,NOPTYP,INTREP,LABINT)
      ELSEIF(CLASLB .EQ. 'EFFDEN ') THEN
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        IORDER    = 0
        CALL SETATM(DOATOM,NATOM,-99)
        CALL EFFDE2TYP(NOPTYP,INTREP,LABINT,DOATOM,NATOM) 
      ELSE
        CALL LSET(NUCIND,.TRUE.,DOATOM)
        IF (PLABEL.EQ.'EDM') THEN
          PLABEL='KINENERG'
          WRITE(LUPRI,*)                                                &
     &     'EDM: PLABEL SET to KINENERG, factor +4*c*field_strength'
        ENDIF
       ! ... call abacus routine...
        CALL PR1DIR(CLASLB,INTTYP,NOPTYP,INTREP,                        &
     &            ANTI,SQUARE,INTADR,LABINT,TRIANG,TRASPH,              &
     &            SOLVNT,IORDER,DOATOM,NATOM,NBAST,NELMNT,IPRINT)
        IF (IPRINT.GE.5) THEN
           WRITE(LUPRI,'(/,2X,A)')                                      &
     &      '*** XPLDE1 after PR1DIR calling ***'
           WRITE(LUPRI,'(2X,A,I3)') 'INTTYP=',INTTYP
           WRITE(LUPRI,*) 'CLASLB=',CLASLB
           WRITE(LUPRI,*) 'PLABEL=',PLABEL
           WRITE(LUPRI,*) 'ANTI=',ANTI
           WRITE(LUPRI,*) 'SQUARE=',SQUARE
           WRITE(LUPRI,*) 'TRIANG=',TRIANG
           WRITE(LUPRI,*) 'NOPTYP=',NOPTYP
           DO ICMP = 1,NOPTYP
             WRITE(LUPRI,*) 'LABINT(1-NOPTYP):',ICMP,LABINT(ICMP)
           ENDDO
        ENDIF
      ENDIF
!
!     Loop over components of operator class
!     ======================================   
!
      LABERR = .TRUE.
      DO ICMP = 1, NOPTYP
        NPRPLBL = NPRPLBL + 1
        IF(NPRPLBL.GT.MAXPRPLBL) THEN
          CALL QUIT('XPLDEF:Redimension MAXPRPLBL !')
        ENDIF
        IF(PLABEL.EQ.LABINT(ICMP)) THEN
        !mi: the only place to reset LABERR, more printout needed for "else" branch
          LABERR = .FALSE.
          INDX   = NPRPLBL
        ELSE
          IF (IPRINT.GE.7) THEN
            WRITE(LUPRI,*) 'PLABEL.NE.LABINT(ICMP), ICMP=',ICMP
            WRITE(LUPRI,*)                                              &
     &      'PLABEL=',PLABEL,' LABINT(ICMP)=',LABINT(ICMP)
          ENDIF
        ENDIF
        PRPLBL (NPRPLBL) = LABINT(ICMP)
        PDOINT (NPRPLBL) = PCOMB
        IPRLREP(NPRPLBL) = INTREP(ICMP)
!       .....  
        IF(SQUARE) THEN
!       Non - symmetric operator
          IPRLTYP(NPRPLBL) = 0
        ELSEIF(ANTI) THEN
!       Anti-symmetric operator
          IPRLTYP(NPRPLBL) = -1
        ELSE
!       Symmetric operator
          IPRLTYP(NPRPLBL) =  1
        ENDIF
      ENDDO
!
!     Property label is not among the listed labels
!
      IF(LABERR) THEN
        WRITE(LUPRI,'(/,A,A8,A)') '*** XPLDE1 ERROR: Property label ',  &
     &        PLABEL,' not found ! *****'
        WRITE(LUPRI,'(A,A7)') 'Class label      : ',CLASLB
        WRITE(LUPRI,'(A,A8)') '  Property labels: ',LABINT(1)
        WRITE(LUPRI,'(19X,A8)') (LABINT(I),I=2,NOPTYP)
        CALL QUIT('XPLDE1:  error - operator label not found')
      ENDIF
!.....Deallocate memory
      call dealloc(INTREP)
      call dealloc(INTADR)
      deallocate(DOATOM)
      deallocate(LABINT)
      CALL QEXIT('XPLDEF')      
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck prpsym */
      SUBROUTINE PRPSYM(INDXPR,IPRINT)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This subroutine finds the symmetry of the unique operator no.
!     INDXPR in COMMON /DCBPRP/.
!
!     The routine assumes Hermitian operators and will therefore
!     consider symmetric operators purely real and antisymmetric operators
!     purely imaginary. Non-symmetric operators are handled by a special routine
!     which terminates in error if the exception has not been added.
!   
!     Written by T.Saue, January 1995
!     Last revision: tsaue - May 14 1996
!
!    MI,HJAaJ(TEC),2003 - added changes for LAO-based  integrals
!    MKN, MI, 2017 - added exception for the EDM operator
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
!
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
!
#include "dcbprl.h"
#include "dcbxpr.h"
#include "symmet.h"
#include "dgroup.h"
#include "pgroup.h"
!
      LOGICAL DIFSYM,DIFTIM, SKIP_JM4
      DIMENSION ICREP(3),ICTIM(3)
#include "ibtfun.h"

      CALL QENTER('PRPSYM')
!
!     Investigate property label no. INDXPR
!
!CTF
!      do IL = 1,22,1
!        print*,'JM4(1,',IL,') = ',JM4(1,IL)
!      end do
!CTF
!CTF
!      do IL = 0,7,1
!        print*,'JM4TIM(',IL,') = ',JM4TIM(IL)
!      end do
!CTF
      ITYP     = IPRPTYP(INDXPR)
      NLBL     = MCMP(ITYP)
      DIFSYM   = .FALSE.
      DIFTIM   = .FALSE.
      IAIND    = JM4(1,ITYP)
      IBIND    = IPRPLBL(1,INDXPR)
!
!     Point group symmetry
!     ====================
!
!     Matrix operator, JM4REP(JM4(1,ITYP))
      IAREP    = JM4REP(IAIND)
!     Scalar operator, IPRLREP (IPRPLBL(1,INDXPR))
      IBREP    = IPRLREP(IBIND)
!     Total point group symmetry of the operator, first component..
      ICREP(1) = IBTXOR(IAREP,IBREP)
!
!     Time reversal symmetry
!     ======================
!
!     Matrix operator :
      IATIM = JM4TIM(IAIND)
!TF
      IF (ITYP.EQ.21) THEN
        IATIM=-1 !
         WRITE(LUPRI,'(2X,A)')                                          &
     &    'i*BETA*GAMMA5 operator, setting IATIM=-1'
!    &    'FYI: i*BETA*GAMMA5 operator' 
!TF
      ENDIF
!     Scalar operator :
!     Here we use the following indirect argument:
!     The operator should be Hermitian, so if the scalar operator is
!     symmetric about the diagonal, then it is real and time symmetric.
!     If the scalar operator is anti-symmetric, then it is imaginary
!     and time anti-symmetric
!     However, if the operator has no symmetry about the diagonal 
!     we can not use this argument and has to add exceptions
!
!MKN/MI Exception for i*BETA*GAMMA5 operator: Malaya K. Nayak (09-07-15)
      IF (ITYP.EQ.21) THEN
         IATIM = -1
         WRITE(LUPRI,'(2X,A)')                                          &
     &   'exception: i*BETA*GAMMA5 operator, setting IATIM = -1'
      END IF
!MKN/MI Because the operator i*BETA*GAMMA5 is a T-odd by its definion

      IBTIM = IPRLTYP(IBIND)
      IF(IBTIM.EQ.0) IBTIM = ISCALAROP_TSYM(PRPLBL(IBIND))
!     Total time reversal symmetry of the 4-comp.operator
!TF   print*,'IATIM,IBTIM = ',IATIM,IBTIM
      ICTIM(1) = IATIM*IBTIM
      IF (IPRINT.GE.4) THEN
        CALL HEADER('Output from PRPSYM',-1)
        WRITE(LUPRI,*) '* Operator type: ',ITYP
        WRITE(LUPRI,*) '* First matrix operator ',IAIND,IAREP,IATIM
        WRITE(LUPRI,*) '* First scalar operator ',IBIND,IBREP,IBTIM
      ENDIF
!
!     Check against other components of operator
!     ==========================================
      DO 20 ILBL = 2,NLBL
        IAIND       = JM4(ILBL,ITYP)
        IBIND       = IPRPLBL(ILBL,INDXPR)
        IAREP       = JM4REP(IAIND)
        IBREP       = IPRLREP(IBIND)
        ICREP(ILBL) = IBTXOR(IAREP,IBREP)
        DIFSYM      = DIFSYM.AND.(ICREP(ILBL).NE.ICREP(1))
        IATIM       = JM4TIM(IAIND)
        IBTIM       = IPRLTYP(IBIND)
        IF(IBTIM.EQ.0) IBTIM = ISCALAROP_TSYM(PRPLBL(IBIND))
        ICTIM(ILBL) = IATIM*IBTIM
        DIFTIM      = DIFTIM.AND.(ICTIM(ILBL).NE.ICTIM(1))
        IF (IPRINT.GE.15) THEN
          WRITE(LUPRI,'(2X,4(A,I2))')                                   &
     &  'IAIND=JM4(ILBL,ITYP)',IAIND,' IBIND=IPRLBL(ILBL,INDXPR)',IBIND,&
     &  ' IATIM=',IATIM,' IBTIM=',IBTIM
        WRITE(LUPRI,'(2X,1(A,I2))')                                     &
     &  'ICTIM(ILBL)=',ICTIM(ILBL)
        ENDIF

   20 CONTINUE
!     error message
      IF(DIFSYM.OR.DIFTIM) THEN
!
!         The components of the operator have not the same symmetry; abort
!
        CALL HEADER('PRPSYM: Abort !!!',-1)
        WRITE(LUPRI,'(A,I5,A)') 'Components of operator no.',INDXPR,    &
     &    ' have not the same symmetry :'
        WRITE(LUPRI,'(A,I5)') 'Operator type:',ITYP
        WRITE(LUPRI,'(32X,A12,23X,A4)')                                 &
     &     'boson irrep','trev'
        DO 30 ILBL = 1,NLBL
          IAREP       = JM4REP(JM4(ILBL,ITYP))
          IBREP       = IPRLREP(IPRPLBL(ILBL,INDXPR))
          WRITE(LUPRI,'(3X,A8,3(3X,A3,A1,I1,A1,3X,A4),I5)')             &
     &      PRPLBL(IPRPLBL(ILBL,INDXPR)),                               &
     &      REP(IAREP),'(',(IAREP+1),')','X',                           &
     &      REP(IBREP),'(',(IBREP+1),')','=',                           &
     &      REP(ICREP(ILBL)),'(',(ICREP(ILBL)+1),')','  |  ',           &
     &      ICTIM(ILBL)
   30   CONTINUE
        WRITE(LUPRI,*) 'FYI: matrix oper, IAREP=JM4REP(:)=',JM4REP
        WRITE(LUPRI,*) '     scalar oper,IBREP=IPRLREP(:)=',IPRLREP
        CALL QUIT('PRPSYM: Components of different symmetry !!!')
      ENDIF

!MI   ... assign point group symmetry - only one !
      IPRPSYM(INDXPR) = ICREP(1)+1
!MI   ... assign time rev.symmetry
!TF   if (INDXPR.eq.7) print*,' ICTIM(1) = ',ICTIM(1)
      IPRPTIM(INDXPR) = ICTIM(1)

!MI ... control printout is always worth 
      IF (IPRINT.GE.15) THEN
        WRITE(LUPRI,'(2X,A,I3,A,I2,A,I2)')                              &
     &  'Operator index, INDXPR =',INDXPR,                              &
     &  ' and operator type (1-',NOPTYP,                                &
     &  '), IPRPTYP(INDXPR) =',IPRPTYP(INDXPR)

        WRITE(LUPRI,'(2X,A,I2)')                                        &
     & 'symmetry of scalar operator under matrix transpose(-1,0,1),'//  &
     & ' IPRLTYP() =',IPRLTYP(INDXPR)
        WRITE(LUPRI,'(2X,A,I1)')                                        &
     &   'Number of operator components, NLBL (1-3) =',NLBL   
        WRITE(LUPRI,'(2X,3(A,I2))')                                     &
     &  'IAIND=',IAIND,' IATIM=',IATIM,' IBTIM=',IBTIM
        WRITE(LUPRI,'(20X,A,14X,A)')                                    &
     &     '--- boson irrep ---','--- time rever symm ---'
        DO ILBL = 1, NLBL
          IAREP       = JM4REP(JM4(ILBL,ITYP))
          IBREP       = IPRLREP(IPRPLBL(ILBL,INDXPR))
          WRITE(LUPRI,'(3X,A8,1X,A4,1X,3(3X,A3,A1,I1,A1,3X,A4),I5)')    &
     &      PRPLBL(IPRPLBL(ILBL,INDXPR)),PDOINT(IPRPLBL(ILBL,INDXPR)),  &
     &      REP(IAREP),'(',(IAREP+1),')','X',                           &
     &      REP(IBREP),'(',(IBREP+1),')','=',                           &
     &      REP(ICREP(ILBL)),'(',(ICREP(ILBL)+1),')','  |  ',           &
     &      ICTIM(ILBL)
        ENDDO
        WRITE(LUPRI,'(2X,A,I2,A,I2/)')                                  &
     & '* assigned IPRPSYM() =',IPRPSYM(INDXPR),                        &
     & ' and time rev symm, IPRPTIM() =',IPRPTIM(INDXPR)
      ENDIF
!
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Output from PRPSYM',-1)
        CALL WRIXPR(0,INDXPR)
      ENDIF
!
      CALL QEXIT('PRPSYM')
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck wrixpr */
      SUBROUTINE WRIXPR(IND,INDXPR)
!***********************************************************************
!
!     Write out data on a operator specified by pointer INDXPR
!
!     Written by T.Saue - May 1996
!     Last revision: May 18 1996 - tsaue
!                    M.Ilias - Febr.2007 - new operator type 20
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0)
!
#include "dcbxpr.h"
#include "dcbprl.h"
#include "pgroup.h"
#include "dgroup.h"
      CHARACTER MXFORM*6,FMT*6,OPTYPE*12,TTYP(-1:1)*2
      SAVE TTYP
      DATA TTYP /'T-','T0','T+'/
!
!     Operator name and symmetry
!     ==========================
!
      ISYM = IPRPSYM(INDXPR)
      ITIM = IPRPTIM(INDXPR)
      IF(IND.EQ.0) THEN
        WRITE(LUPRI,'(6X,A16,5X,A3,5X,A2)')                             &
     &    PRPNAM(INDXPR),REP(ISYM-1),TTYP(ITIM)
      ELSE
        WRITE(LUPRI,'(2X,I3,2X,A16,5X,A3,5X,A2)')                       &
     &   IND,PRPNAM(INDXPR),REP(ISYM-1),TTYP(ITIM)
      ENDIF
      CALL PRSYMB(LUPRI,'.',75,0)
!
!     Operator type
!     =============
!
      ITYP = IPRPTYP(INDXPR)
      IF(ITYP.EQ.0) THEN
        WRITE(LUPRI,'(6X,A,A8)') 'Operator read from matrix ',          &
     &        PRPLBL(IPRPLBL(1,INDXPR))
        FMT = MXFORM(ABS(FACPRP(1,INDXPR)),16)
        WRITE(LUPRI,'(6X,A,4X,'//FMT//')')                              &
     &      'Factor     : ',FACPRP(1,INDXPR)
      GOTO 30
      ENDIF
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,                                &
     &       13,14,15,16,17,18,19,20,21), ITYP
        WRITE(LUPRI,'(6X,A,I5)')                                        &
     &  'Program error: unknown operator type:',ITYP
        CALL QUIT('WRIXPR: unknown operator type')
    1 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'DIAGONAL : scalar operator'
      GOTO 200
    2 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'XALPHA   : [alpha_x]P'
      GOTO 200
    3 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'YALPHA   : [alpha_y]P'
      GOTO 200
    4 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'ZALPHA   : [alpha_z]P'
      GOTO 200
    5 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'XAVECTOR : [alpha x P]_x'
      GOTO 200
    6 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'YAVECTOR : [alpha x P]_y'
      GOTO 200
    7 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'ZAVECTOR : [alpha x P]_z'
      GOTO 200
    8 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'ALPHADOT : [alpha.P]'
      GOTO 200
    9 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'GAMMA5   : [gamma_5]P'
      GOTO 200
   10 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'XSIGMA   : [Sigma_x]P'
      GOTO 200
   11 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'YSIGMA   : [Sigma_y]P'
      GOTO 200
   12 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'ZSIGMA   : [Sigma_z]P'
      GOTO 200
   13 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'XBETASIG : [beta Sigma_x]P'
      GOTO 200
   14 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'YBETASIG : [beta Sigma_y]P'
      GOTO 200
   15 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'ZBETASIG : [beta Sigma_z]P'
      GOTO 200
   16 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'XiBETAAL : i[beta alpha_x]P'
      GOTO 200
   17 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'YiBETAAL : i[beta alpha_y]P'
      GOTO 200
   18 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'ZiBETAAL : i[beta alpha_z]P'
      GOTO 200
   19 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'BETA     : scalar operator'
      GOTO 200
   20 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'SIGMADOT : [sigma.P]'
      GOTO 200
   21 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',                       &
     &   'BETAGAMMA5 : i[beta gamma5]P'
      GOTO 200
!
  200 CONTINUE
!
!     Labels and factors
!     ==================
!
      NLBL   = MCMP(ITYP)
      DO I = 1,NLBL
        FMT = MXFORM(ABS(FACPRP(I,INDXPR)),16)
        KLBL = IPRPLBL(I,INDXPR)
        IF    (IPRLTYP(KLBL).EQ. 1) THEN
          OPTYPE = ' (real)     '
        ELSEIF(IPRLTYP(KLBL).EQ.-1) THEN
          OPTYPE = ' (imaginary)'
        ELSE
          OPTYPE = ' (complex)  '
        ENDIF
        IF(I.EQ.1) THEN
          WRITE(LUPRI,'(6X,A,A8,1X,A4,4X,'//FMT//',A12)')               &
     &      'Labels and factors     : ',                                &
     &      PRPLBL(IPRPLBL(I,INDXPR)),PDOINT(IPRPLBL(I,INDXPR)),        &
     &      FACPRP(I,INDXPR),OPTYPE
        ELSE
          WRITE(LUPRI,'(31X,A8,1X,A4,4X,'//FMT//',A12)')                &
     &      PRPLBL(IPRPLBL(I,INDXPR)),PDOINT(IPRPLBL(I,INDXPR)),        &
     &      FACPRP(I,INDXPR),OPTYPE
        ENDIF
      ENDDO
 30   CONTINUE
      CALL PRSYMB(LUPRI,'.',75,0)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck op1ind */
      SUBROUTINE OP1IND(OP1KEY,IND1OP,IPR1OP,N1OPER,INDXPR,MAX1OP)
!***********************************************************************
!
!     INDXPR points to an operator in XCBPRP. 
!     Search through pointer array IPR1OP for elements pointing 
!     to this operator. If not found, extend the list by one element.
!     The program is stopped if the array goes out of the bound MAX1OP.
!     On output IND1OP will be the index of the operator in the 
!     pointer array IPR1OP, N1OPER the current number of operators
!     in the list.
!
!     Written by T.Saue June 28 1996
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER OP1KEY*(*) ! MI: instead of previous OP1KEY*6; fix for gfortran 4.9
      DIMENSION IPR1OP(MAX1OP)
!
!     Search for INDXPR
!     =================
!
      DO J = 1,N1OPER
        IF(INDXPR.EQ.IPR1OP(J)) THEN
          IND1OP = J
          RETURN
        ENDIF
      ENDDO
!
!     Not found --> Add to list
!     =========================
!
      N1OPER = N1OPER + 1
      IND1OP = N1OPER
      IF(N1OPER.GT.MAX1OP) THEN
        WRITE(LUPRI,'(A6,A)') OP1KEY,' out of bounds in OP1IND'
        WRITE(LUPRI,'(A,I5)') 'Current value of MAX1OP: ',MAX1OP
        CALL QUIT('OP1IND: Out of bounds !')
      ELSE
        IPR1OP(IND1OP) = INDXPR
      ENDIF
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck op1srt */
      SUBROUTINE OP1SRT(OP1KEY,IPOP1U,IPOP1S,NOP1T,NOP1,JOP1,           &
     &           WORK,LWORK)
!***********************************************************************
!
!     Sort a list of operators in symmetries
!
!        IPOP1U - input list of operators
!        IPOP1S - output list of operators sorted on symmetries        
!        NOP1T  - total number of operators
!        NOP1   - number of operators in ISYM
!        JOP1   - offset to operators in ISYM        
!      
!     Written by T.Saue June 28 1996
!     Last revision - June 26 1996 - tsaue
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
!
      CHARACTER OP1KEY*(*) ! MI fix for gfortran4.9, instead of  OP1KEY*6
      DIMENSION IPOP1U(*),IPOP1S(*),NOP1(8,2),JOP1(8,2),WORK(LWORK)
!
      CALL QENTER(OP1KEY)
#include "memint.h"
      NBUF = NOP1T*8
      CALL MEMGET('INTE',KIBUF,NBUF,WORK,KFREE,LFREE)
      CALL OP1SR1(IPOP1U,IPOP1S,NOP1T,NOP1,JOP1,WORK(KIBUF))
      CALL QEXIT (OP1KEY)      
!   
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck op1sr1 */
      SUBROUTINE OP1SR1(IPOP1U,IPOP1S,NOP1T,NOP1,JOP1,IBUF)
!***********************************************************************
!
!     Sort a list of operators in symmetries
!
!     Written by T.Saue June 28 1996
!     Last revision - June 26 1996 - tsaue
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
!
#include "dcbxpr.h"
#include "dgroup.h"
      DIMENSION IPOP1U(NOP1T),IPOP1S(NOP1T),NOP1(8),                    &
     &          JOP1(8),IBUF(NOP1T,8)
!
!     Sort in buffer
!
      DO J = 1,NOP1T
        INDXPR = IPOP1U(J)
        ISYM   = IPRPSYM(INDXPR)
        NOP1(ISYM) = NOP1(ISYM) + 1
        IBUF(NOP1(ISYM),ISYM) = J
      ENDDO
!   
!     Pack results in IPOP1U
!
      NFC  = NBSYM/NFSYM
      IOFF = 0
      DO IFRP = 1,NFSYM
        DO IS = 1,NFC
          ISYM = JFSYM(IS,IFRP)
          NBUF = NOP1(ISYM)
          JOP1(ISYM) = IOFF
          IF(NBUF.GT.0) THEN
            DO J = 1,NBUF
              IPOP1S(IOFF+J) = IBUF(J,ISYM)
            ENDDO            
            IOFF = IOFF + NBUF
          ENDIF
        ENDDO
      ENDDO
!
      RETURN
      END
!  /* Deck INIGRD */
      SUBROUTINE INIGRD(IPRGRD_inp)
!***********************************************************************
!
!     Initialization for molecular gradient
!
!     Written by J. Thyssen - 1997/09/01
!     Last revision: jth - 1997/09/01
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
!
#include "maxorb.h"
#include "mxcent.h"
!
#include "dcbgrd.h"
#include "dcbgen.h"
#include "dcbham.h"
!
#include "ibtfun.h"
!
!     Initialize /CBIGRD/
!     ===================
!
      IPRGRD = IPRGRD_inp
!
!     RELGRD is used in DHFGRD to control integrals
!
      RELGRD = .FALSE.
      DONGRD = .FALSE.
      IGRD_INTFLG = INTGEN
      IF (LEVYLE) IGRD_INTFLG = IBTAND ( 1, IGRD_INTFLG )
      DOTRCK    = .FALSE.
      SCRGRD    = -1.0D0      
!
      RETURN
      END
!***********************************************************************
!  /* Deck grdinp */
      SUBROUTINE GRDINP(WORD,RESET)
!***********************************************************************
!
!     Input section for molecular gradient
!
!     Written by J. Thyssen - 1997/07/09
!     Last revision: jth - 1997/07/09
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0=0.0D0)
#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 5)
!
#include "dcbham.h"
#include "dcbgrd.h"
#include "dcbprp.h"
!
      LOGICAL SET, NEWDEF, RESET, LBIT
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,             &
     &          PNAME*16, PLABEL(3)*8, STR*3
      DIMENSION PFAC(3),II(3)
!
      SAVE SET
      DATA TABLE /'.PRINT ','.INTFLG','.TRICK ','.SCREEN','.NUMGRA'/
      DATA SET/.FALSE./
!
      NEWDEF = (WORD .EQ. '*MOLGRD')
      IF (SET) THEN
         IF (NEWDEF)                                                    &
     &      CALL QUIT('Only one "*MOLGRD" input section allowed')
!        hjaaj: repeated input sections give infinite loop ...
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
!
      CALL INIGRD(IPRPRP)
!
!
!     Process input for DCBEXP
!     ========================
!
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,              &
     &            '" not recognized in GRDINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in GRDINP.')
    1          CONTINUE
!&&&& PRINT:  Print level
                  READ(LUCMD,*) IPRGRD
               GO TO 100
    2          CONTINUE
!&&&& INTFLG: Which two-electron gradient to calculate
                  IF (GAUNT) THEN 
                     READ(LUCMD,*) ILL,ISL,ISS,IGT
                  ELSE
                     READ(LUCMD,*) ILL,ISL,ISS
                  END IF
                  READ(LUCMD,*) ILL, ILS, ISS
                  IGRD_INTFLG = ILL + 2 * ILS + 4 * ISS + 8 * IGT
               GO TO 100
    3          CONTINUE
!&&&& TRICK : Skip LS- and SS- integrals if contribution negligible
                  DOTRCK = .TRUE.
               GO TO 100
    4          CONTINUE
!&&&& SCREEN: Threshold for screening on the gradient
                  READ (LUCMD,*) SCRGRD
               GO TO 100
    5          CONTINUE
!&&&& NUMGRA: Use numerical gradient
                  DONGRD = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,               &
     &            '" not recognized in GRDINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in GRDINP.')
            END IF
      END IF
  300 CONTINUE
!
!     Print section
!     =============
      IF (MGRAD) THEN
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(A)')                                             &
     &        ' GRDINP: Calculation of molecular gradient.'
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(1X,A,I3)')                                       &
     &        '* Print level                       : ', IPRGRD
         IF ( DONGRD ) THEN
            WRITE(LUPRI,'(1X,A)')                                       &
     &         '* Gradient evaluated with numerical differentiation.'
         ELSE
            if (bss) then
              call quit('analytical molecular '//                       &
     &                  'gradient not implemented for BSS/X2C/IOTC')
            end if
            WRITE(LUPRI,'(1X,A)')                                       &
     &         '* Gradient evaluated analytically.'
         END IF
         WRITE(LUPRI,'(1X,A)')                                          &
     &        '* Contributions to the molecular gradient from'
         IF ( LBIT( IGRD_INTFLG, 1 ) )                                  &
     &         WRITE(LUPRI,'(3X,A)') '- LL-integrals.'
         IF ( LBIT( IGRD_INTFLG, 2 ) ) THEN
            IF ( DOTRCK ) THEN
               WRITE(LUPRI,'(3X,A)')                                    &
     &            '- LS-integrals (skipped if estimated to be small).'
            ELSE
               WRITE(LUPRI,'(3X,A)') '- LS-integrals.'
            END IF
         END IF
         IF ( LBIT( IGRD_INTFLG, 3 ) ) THEN
            IF ( DOTRCK ) THEN
               WRITE(LUPRI,'(3X,A)')                                    &
     &            '- SS-integrals (skipped if estimated to be small).'
            ELSE
               WRITE(LUPRI,'(3X,A)') '- SS-integrals.'
            END IF
         END IF
         IF ( LBIT( IGRD_INTFLG, 4 ) ) THEN
            IF ( DOTRCK ) THEN
               WRITE(LUPRI,'(3X,A)')                                    &
     &            '- GT-integrals (skipped if estimated to be small).'
            ELSE
               WRITE(LUPRI,'(3X,A)')                                    &
     &         '- GT-integrals (WARNING: gradient not tested yet)'
            END IF
         END IF
         IF(SCRGRD.GT.D0) THEN
            WRITE(LUPRI,'(A,1P,E8.2)')                                  &
     &         ' * Screening threshold in gradient calculation ',SCRGRD
         ELSE
            WRITE(LUPRI,'(A)')                                          &
     &         ' * No screening in the calculation of '//               &
     &         'two-electron integrals.'
         END IF
         CALL PRSYMB(LUPRI,'-',75,0)
      END IF
  999 CONTINUE
      RETURN
!
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck minprp */
      SUBROUTINE MINPRP()
!***********************************************************************
!     Calculate specified properties at final, converged geometry
!     (called when .OPTIMIZE specified if geometry converges)
!
      use memory_allocator
#include "implicit.h"
#include "dcbxpr.h"
#include "priunit.h"
#include "dcbgen.h"
!
      real(8), allocatable :: WORK(:)
      LOGICAL FOPEN
!
      CALL QENTER('MINPRP')


!
      NPRPS = 0
!
!
!     Open LU1INT if not opened
!
      IF (DOPRP) THEN
         INQUIRE(FILE='AOPROPER',OPENED=FOPEN)
         IF (.NOT. FOPEN) THEN
            OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',           &
     &         FILE='AOPROPER')
         END IF
!
         call legacy_lwork_get(LWORK)
         call alloc(WORK,LWORK,id='WORK in MINPRP')
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL PRPINP('**PRP F',WORK,LWORK)
         CLOSE(LUCMD)
!
         CALL ONEGEN(WORK,LWORK)
         call dealloc(WORK)
!
!        Close AOPROPER if it was closed on entry to
!        MINPRP.
!
         IF (.NOT. FOPEN) CLOSE(LU1INT,STATUS = 'KEEP')
!
!
!***********************************************************************
!*****  R E S P O N S E    M O D U L E  ********************************
!***********************************************************************
!
         CALL PAMPRP()
!
      END IF
!
      CALL QEXIT('MINPRP')
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      INTEGER FUNCTION ISCALAROP_TSYM(PRPLBL)
      CHARACTER*8 PRPLBL
      INTEGER ITIM
      IF(PRPLBL(4:7).EQ.'RNST'.OR.                                      &
     &   PRPLBL(3:7).EQ.'RM2H3'.OR.                                     &
     &   PRPLBL(3:7).EQ.'RDSUL'.OR.                                     &
     &   PRPLBL(3:7).EQ.'RDSUS'.OR.                                     &
     &   PRPLBL(1:7).EQ.'d|S>/dB') THEN
         ITIM = -1
      ELSEIF(PRPLBL(2:6).EQ.'RM1H3'.OR.                                 &
     &   PRPLBL(4:7).EQ.'NSNL'.OR.                                      &
     &   PRPLBL(2:6).EQ.'1SOPP'.OR.                                     &
     &   PRPLBL(3:8).EQ.'DSUSNL'.OR.                                    &
     &   PRPLBL(2:6).EQ.'RM1RN') THEN
         ITIM = +1
      ELSE
        WRITE(6,*) 'ISCALAROP_TSYM ERROR :',PRPLBL
        CALL QUIT('ISCALAROP_TSYM: Can not assign ITIM')
      ENDIF
      ISCALAROP_TSYM=ITIM
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

