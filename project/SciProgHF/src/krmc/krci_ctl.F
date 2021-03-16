!      Copyright (c) 2014 by the authors of DIRAC.
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

! FILE : krmc/krci_ctl.F
!
! Specific routines for .KRCI calculations


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRCI_INP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     CI optimization control input section.
C
C     Written by S. Knecht and Hans Joergen Aa. Jensen
C
C                    'old statements':
C
C                    Added CI keyword CIROOT. Timo Fleig    05.01.2005
C                    Added CI keyword ANALYZ. Timo Fleig
C                    Added CI keyword PARINT. Stefan Knecht 29.06.2007
C                    Added CI keyword MEMFAC. Stefan Knecht 04.07.2007
C                    Added CI keyword RDFOCK. Stefan Knecht October 2007
C                    Added CI keyword WRTFCK. Stefan Knecht
C                    Added CI/KRMC keyword SHMEMO. Stefan Knecht 20.11.2007
C                    Added CI/KRMC keyword IJKLRO. Stefan Knecht
C
C     Major revision by S. Knecht: Aug 2008
C
C                    All CI-related keywords have been moved to this new
C                    input section *KRCICALC.
C                    Purpose: disentanglement of MCSCF/CI optimization
C                    input and start-up --> "modularize"
C
C                    Added CI keyword EEDM   Timo Fleig     30.04.2012
C                    Added CI keyword MHYP   Malaya Nayak   28.08.2013
C                    Added CI keyword ENSPS  Malaya Nayak   27.02.2014
C                    Added CI keyword TERSE  Timo Fleig     29.03.2016
C                    Added CI keyword NMQM   Malaya Nayak   28.09.2018
C
C***********************************************************************
      use symmetry_setup_krci
      use krci_cfg
      use os_utils
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxash.h"
      PARAMETER (NTABLE = 52)
      PARAMETER (RTOL = 1.0D-15,D1 = 1.0D0,D0=0.0D0,D2=2.0D00)
C
      LOGICAL SET, NEWDEF, RESET, LBIT
      CHARACTER*4 REPNA(64)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, TMPLAB*8
      CHARACTER CTEMP*72, LINE*80, TEXT*20, REPA(8)*4
      DIMENSION WORK(LWORK)
      INTEGER IKWSET(NTABLE)
C
#include "maxorb.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbopt.h"
#include "dcbprp.h"
#include "dcbxpr.h"
#include "nuclei.h"
#include "krciprop.h"
#include "dcborb.h"
#include "dcbkrci.h"
#include "dcbpsi.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "pgroup.h"
C
      SAVE SET
      DATA TABLE /'.MAX CI','.OPERAT','.RSTRCI','.ANALYZ','.CIROOT',
     &            '.MXCIVE','.KTRLVL','.PARINT','.MEMFAC','.NOOCCN',
     &            '.RDFOCK','.WRTFCK','.SHMEMO','.IJKLRO','.IJKLSP',
     &            '.PRINT ','.FCIDUM','.NATOCC','.XPSINT','.INTFLG',
     &            '.THRPCI','.MXVCBL','.MIN MK','.MAX MK','.GASSH ',
     &            '.GASSPC','.MK2REF','.MK2DEL','.CI PRO','.INACTI',
     &            '.CASSCF','.TRDM  ','.OMEGAQ','.THRECI','.CHECKP',
     &            '.USE PB','.GENFOC','.TRA NO','.LOWSRT','.SEL ST',
     &            '.DIPMOM','.MCANA ','.GAS SH','.SVRONO','.SVRONX',
     &            '.EEDM  ','.MHYP  ','.ENSPS ','.NMQM  ','.XXXXXX',
     &            '.XXXXXX','.TERSE '/ !MKN
      DATA SET/.FALSE./
      character (len= 4)              :: symmetry_str
      character (len= 4)              :: symmetry_str1
      character (len= 4)              :: symmetry_str2
      character (len= 4)              :: tmp_symmetry_str
      character (len= 4), allocatable :: nkrci_ciroots_symstr(:)
      character (len= 1)              :: last_char
      character (len= 2)              :: ferm_lab
      logical                         :: gerade_mcci_wf
      integer                         :: tmp_mij_val
      integer                         :: fermion_sym
      integer                         :: ascii_lchar
      integer                         :: double_group
      integer, allocatable            :: multb_tmp(:,:)
      character(len=400)              :: input_line
      integer                         :: ios, islash
C
#include "ibtfun.h"
#include "memint.h"
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
      CALL IZERO(IKWSET,NTABLE)
C
C     Initialize /CBLKRCI/
C     ===================
C
      KRCI_UCIBOS       = .FALSE.
      CRDFO_MAT_KRCI    = .FALSE.
      CWRTFO_MAT_KRCI   = .FALSE.
      CSHMEMO_KRCI      = .FALSE.
      CINT_REORD_KRCI   = .FALSE.
      CINT_SPLIT_KRCI   = .FALSE.
      NATOLCR_KRCI      = .FALSE.
      CHCKPT_KRCI       = .FALSE.
      CUSE_PBLFL_KRCI   = .FALSE.
      GENFOCK_KRCI      = .FALSE.
      TRANATO_KRCI      = .FALSE.
      LOWSORT_KRCI      = .FALSE.
      STATE_SELECT_KRCI = .FALSE.
      LOWSORT_KRCI      = .FALSE.
      save_reordered_nos_krci = .false.
      fcidump_krci      = .false.
      XPSINTERFACE_KRCI = .FALSE.
C
C     Initialize /CBIKRCI/
C     ===================
C
      DO I = 1,MAX_NKRCI_MAX_SYM
         NKRCI_CIROOTS(I) = 0
         NKRCI_SELOM(I)   = 0
      END DO
      NKRCI_MAX_SYM   = -1
      ITERSEOUT = 0

      IPRKRCI   = IPRGEN
      ILLINT    = IBTAND(INTGEN,1)
      ISLINT    = IBTAND(INTGEN/2,1)
      ISSINT    = IBTAND(INTGEN/4,1)
      IKRCI_INTBUF    = 0
      IKRCI_INTFLG    = 0
      IKRCI_ITRINT(1) = 0
      IKRCI_ITRINT(2) = 0

      NKRCIGAS        = -9999
      NKRCI_MK2REF    = -9999
      NKRCI_MK2DEL    = -9999
      NKRCI_MINMK2    = -9999
      NKRCI_MAXMK2    = -9999
C
C     Initialize /CBRKRCI/
C     ===================
C
      DKRCI_CNVINT(1) = DUMMY
      DKRCI_CNVINT(2) = DUMMY
      DKRCI_THRPCI = 1.0D-3
CTF   DKRCI_THRECI = 1.0d-8
CTF   DKRCI_THRECI = 0.3d-10   ! Seems to converge P,T-odd properties
C                                sufficiently well. More testing
C                                required.
      DKRCI_THRECI = 0.8d-12
C
C     Initialize /DCOKRCI/
C     ===================
C
      DO I = 1,2
         NKRCIISH(I) = 0
         NKRCIASH(I) = 0
         NKRCISSH(I) = 0
         NKRCIPSH(I) = 0
         NKRCIFRO(I) = 0
         IKRCI_SVRONO(i)  = 0
         DO J = 1, MXGAS
            NKRCIGSH(I,J) = 0
            NKRCIGSP(I,J) = 0
         END DO
      END DO
      NKRCIAELEC = -1
C
C     Initialize /DCCKRCI/
C     ===================
C
      KRCI_CIPROGRAM = 'GASCIP'
C
C     Initialize /DCICKRCI/
C     ===================
C
      MAXCIT_KRCI   =    4
      MXCIV_KRCI    =    0
      IANACI_KRCI   =    0
      KTRLVL_KRCI   =    3
      IRESTRK_KRCI  =    0
      IKRCI_MEMFAC  =    9
#ifdef VAR_MPI
      IDOPARIO_KRCI = get_environment_integer('GLBSCR',IDOPARIO_KRCI)
#endif
C     ... in "krciprop.h"
      NPROP_KRCI    =    0
      DOPROP_KRCI   = .FALSE.
      DOSYMOPPRP    = .FALSE.
      DOJZEXP       = .FALSE.
      DOOSCILLST    = .FALSE.
      DODIPMOM_KRCI = .FALSE.
      DOMCANA_KRCI  = .FALSE.
      RUNXPROP      = .FALSE.
      DOGENP_KRCI   = .FALSE.
      DOEEDM_KRCI   = .FALSE.
      DOMHYP_KRCI   = .FALSE.
      DOENSPS_KRCI  = .FALSE.
      DONMQM_KRCI   = .FALSE. !MKN
C     ... in dcbopt.h
      ISHMEM_TYPE =  - 1
C
C     some scratch
      IKRCI_SYMSEL      = 0
      NOMEGASEL_KRCI    = 0
      NOMEGASEL_statenr = 1


      if(IDOPARIO_KRCI .ne. 0 .or. IDOPARIO_KRCI.ne. 1)then
!       fall back to default: 0
        IDOPARIO_KRCI = 0
      end if

!     krci_cfg
      krci_cfg_max_vec_block = -1

C
C     Process input for *KRCICALC
C     ===========================
C
      NEWDEF = (WORD .EQ. '*KRCICA' .OR. WORD .EQ. '*KRCI  ')
      ICHANG = 0
      IF (NEWDEF) THEN
C
!       initialize symmetry information irreps required for input
!       processing. note, we assume the # double group irreps here to be 128. 
!  FIXME: make it general in Dirac for higher linear symmetry double groups. 
!         "global Dirac post-HF symmetry module".
C
        if(linear)then
          double_group = 11
          if(nfsym .eq. 1) double_group = 10
          call symmetry_setup_init(double_group,max_nkrci_max_sym)
          allocate(nkrci_ciroots_symstr(4*max_nkrci_max_sym))
        end if

        WORD1 = WORD
  100   CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24,25,26,27,28,29,30,
     &                     31,32,33,34,35,36,37,38,39,40,
     &                     41,42,43,44,45,46,47,48,49,50,
     &                     51,52)
     &                     , I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in *KRCICALC'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in *KRCICALC')
    1          CONTINUE
C&&& .MAX CI ITER...: maximum number of initial CI iterations.
                  READ(LUCMD,*) MAXCIT_KRCI
               GO TO 100
    2          CONTINUE
C&&& .OPERATOR: define one-electron operator to use in
C               CI property calculation
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        IPRKRCI)
                  CALL OP1IND('KRCIPR',LUCIINDTRP,LPROP_KRCI,
     &                        NPROP_KRCI,INDXPR,MXPROP_KRCI)
                  DOPROP_KRCI = .TRUE.
                  DOGENP_KRCI = .TRUE.
               GO TO 100
    3          CONTINUE
C&&& .RSTRCI: Restart CI from vector on file (1) or not (0)
                  READ(LUCMD,*) IRESTRK_KRCI
               GO TO 100
    4          CONTINUE
C&&& .ANALYZ: Analyze LUCIAREL CI vector(s)
                  IANACI_KRCI = 1
               GO TO 100
    5          CONTINUE
C&&& .CIROOTS: Symmetry and number of CI roots to be optimized in this
C              symmetry
                  if(.not.linear) then
                    READ(LUCMD,*) IKRCI_SYMMETRY,NCIROOT_KRCI
                    IF(IKRCI_SYMMETRY .LE. MAX_NKRCI_MAX_SYM .AND.
     &                 IKRCI_SYMMETRY .GT. 0) THEN
                      NKRCI_CIROOTS(IKRCI_SYMMETRY) = NCIROOT_KRCI
                      NKRCI_MAX_SYM = MAX(NKRCI_MAX_SYM, IKRCI_SYMMETRY)
                    ELSE
                      WRITE(LUPRI,'(//A/I10,A,I3)')
     &                '*** ERROR in KRCI_INP .CIROOTS *** ' //
     &                'Requested symmetry outside hardwired limits:',
     &                IKRCI_SYMMETRY, ' outside 1 to',MAX_NKRCI_MAX_SYM
                      CALL QUIT('*** ERROR in KRCI_INP .CIROOTS *** ' //
     &                'Requested symmetry outside hardwired limits.')
                    END IF
                  else
!                   read mij value and # roots in linear symmetry input format
                    read(lucmd,*) symmetry_str, NCIROOT_KRCI

!                   manipulate mij-value input string for further processing
!                   print *,'symmetry_str is ==> ',symmetry_str
                    tmp_symmetry_str = trim(symmetry_str)
                    length_of_string = len_trim(tmp_symmetry_str)
!                   print *,'length of symmetry_str is ==> ',
!    &                       length_of_string
                    last_char = tmp_symmetry_str(length_of_string:
     &                                           length_of_string)
!                   print *,'last char of symmetry_str is ==> ',
!    &                       last_char

!                   determine symmetry irrep of CI wave function
!                   --------------------------------------------

!                   a. check for inversion symmetry
                    ascii_lchar = iachar(last_char)
                    if(nfsym .gt. 1)then
                      if(ascii_lchar .eq. iachar('g'))then
                        gerade_mcci_wf = .true.
                      else if(ascii_lchar .eq. iachar('u'))then
                        gerade_mcci_wf = .false.
                      else
                        write(lupri,'(/a,/a)')
     & '  *** input error: Dinfh contains the inversion operation'//
     & ' as group element.',
     & '  please specify either g or u in addition to the mj value. ***'
                        call quit('  *** input error: Dinfh contains
     & the inversion operation as group element. ***')
                      end if
                    else
                      gerade_mcci_wf   = .true.
                      length_of_string = length_of_string + 1
                      if(ascii_lchar .eq. iachar('g')
     &                   .or. ascii_lchar .eq. iachar('u'))then
                        call quit(' *** input error: Cinfv does not
     & have inversion as group element. ***')
                      end if
                    end if

!                   b. get mij value in integer format
                    read(tmp_symmetry_str(1:length_of_string-1),
     &              fmt='(i6)') tmp_mij_val
                    fermion_sym = 1
                    if(.not.gerade_mcci_wf) fermion_sym = 2
!                   print *,'fermion_sym, tmp_mij_val are ==> ',
!    &                       fermion_sym, tmp_mij_val
!                   call quit('==> bla bla ***')

!                   retrieve active dbg irrep (routine is inside the module symmetry_setup_krci)
                    call convert_mj_ferm_2_dbg_irrep_linsym_krci(
     &                   mj2rep,max_nkrci_max_sym,tmp_mij_val,
     &                   fermion_sym,is_current_dbg_irrep)

!                   c. store # roots for the active irrep
                    nkrci_ciroots(is_current_dbg_irrep) = nciroot_krci
                    nkrci_max_sym = max(nkrci_max_sym,
     &                                  is_current_dbg_irrep)
                    nkrci_ciroots_symstr(is_current_dbg_irrep) =
     &              tmp_symmetry_str(1:length_of_string)
                  end if
               GO TO 100
    6          CONTINUE
C&&& .MXCIVE: Number of CI vectors allowed for subspace
                  READ(LUCMD,*) MXCIV_KRCI
               GO TO 100
    7          CONTINUE
C&&& .KTRLVL: default integral transformation level.
                 READ(LUCMD,*) KTRLVL_KRCI
               IF (KTRLVL_KRCI .LE. -2 .OR. KTRLVL_KRCI .GT. 5) THEN
                  WRITE(LUPRI,'(//2A,I10)')
     &                '*** ERROR in KRCI_INP *** ',
     &                'Illegal transformation level: ',KTRLVL_KRCI
                  CALL QUIT('*** ERROR in KRCI_INP *** ' //
     &                 'Illegal transformation level')
               END IF
               GO TO 100
    8          CONTINUE
C&&& .PARINT: integral handling for slaves in LUCIAREL calculations
               IDOPARIO_SAVE = IDOPARIO_KRCI
               READ(LUCMD,*) IDOPARIO_KRCI
               IDOPARIO_KRCI = ABS(IDOPARIO_KRCI)
               IF( IDOPARIO_KRCI .ne. IDOPARIO_SAVE )THEN
C                 we do not want to broadcast integrals if not
C                 necessary! - SK- 20-03-2008
                  IDOPARIO_KRCI = IDOPARIO_SAVE
               END IF
#if !defined (VAR_MPI)
               IF( IDOPARIO_KRCI .ne. 0 ) THEN
                  WRITE(LUPRI,'(//3A,I6)')
     &                  '*** ERROR in KRCI_INP *** ',
     &                  ' No slaves in sequential run, no integral
     &                    handling needed: ',
     &                   IDOPARIO_KRCI
                  IDOPARIO_KRCI = 0
               END IF
#else
               IF( IDOPARIO_KRCI .gt. 1 ) THEN
                  WRITE(LUPRI,'(//3A,I6)')
     &                  '*** ERROR in KRCI_INP *** ',
     &                  'Illegal integral handling option for
     &                   slaves in parallel run: ',
     &                   IDOPARIO_KRCI
                  IDOPARIO_KRCI = 0
                END IF
#endif
               GO TO 100
    9          CONTINUE
C&&& .MEMFAC: multiplier for subtracted scratch memory in
C             luciarel/ciinfo_r.F - Z_BLKFO_REL
               READ(LUCMD,*) IKRCI_MEMFAC
               IF( IKRCI_MEMFAC .LE. 0 ) THEN
                  WRITE(LUPRI,'(//3A,I6)')
     &                  '*** ERROR in KRCI_INP *** ',
     &                  ' Illegal multiplier for ',
     &                  ' subtracting scratch memory: ',
     &                   IKRCI_MEMFAC
                  CALL QUIT('*** ERROR in KRCI_INP *** ' //
     &                 'Illegal memory multiplier')
               END IF
               GO TO 100
   10          CONTINUE
C&&& .NOOCCN: calculate natural orbitals occupation numbers
                  NATOLCR_KRCI = .TRUE.
               GO TO 100
   11          CONTINUE
C&&& .RDFOCK: read fock matrix from file useful for LUCIAREL restarts
                  CRDFO_MAT_KRCI = .TRUE.
               GO TO 100
   12          CONTINUE
C&&& .WRTFCK: write fock matrix in CI start guess module
                  CWRTFO_MAT_KRCI = .TRUE.
               GO TO 100
   13          CONTINUE
C&&& .SHMEMO: run parallel CI/(KRMC) in MPI 'shared memory' mode
                  CSHMEMO_KRCI = .TRUE.
                  READ(LUCMD,*) ISHMEM_TYPE
#if !defined (VAR_MPI2) && !defined (VAR_SHMEM)
               IF( CSHMEMO_KRCI ) THEN
                  WRITE(LUPRI,'(//4A)')
     &               '*** ERROR in KRCI_INP *** ',
     &               ' shared memory calculations in sequential',
     &               ' mode and/or without enabling the shared-memory',
     &               ' parallelization are not permitted '
                  WRITE(LUPRI,'(/A)')
     &               ' re-configure with --sharedmem_par and'//
     &               ' a parallel Dirac.'
                  CSHMEMO_KRCI = .FALSE.
                  WRITE(LUPRI,'(/A)') '  User input will be ignored.'
               END IF
#else
               IF( ISHMEM_TYPE .gt. 6 ) THEN
                  WRITE(LUPRI,'(//A)')
     &               ' *** ERROR in KRCI_INP *** '
                  WRITE(LUPRI,'(//A)')
     &               ' unknown shared memory run-level! '
                  WRITE(LUPRI,'(A)')
     &               ' LEVEL 0: no shared memory in use'
                  WRITE(LUPRI,'(A)')
     &               ' LEVEL 1: sharing node-wise (ij|kl)'
                  WRITE(LUPRI,'(A)')
     &               ' LEVEL 2: sharing node-wise |C> '
                  WRITE(LUPRI,'(A)')
     &               ' LEVEL 3: sharing node-wise (ij|kl) + |C>'
                  WRITE(LUPRI,'(A)')
     &               ' LEVEL 4: sharing global (ij|kl) + node-wise |C>'
                  WRITE(LUPRI,'(A)')
     &               ' LEVEL 5: sharing global |C> + node-wise (ij|kl)'
                  WRITE(LUPRI,'(A)')
     &               ' LEVEL 6: sharing global (ij|kl) + |C> '
                  WRITE(LUPRI,'(//A,I3)')
     &               ' you asked for LEVEL ',ISHMEM_TYPE
                  CALL QUIT('*** ERROR in KRCI_INP *** ' //
     &                  ' unknown shared memory run-level')
               END IF
#endif
               GO TO 100
   14          CONTINUE
C&&& .IJKLRO: do GAS-scheme specific integral reordering resp. use
C             reordered integrals
                  CINT_REORD_KRCI = .TRUE.
C             ... this requires the use of MPI-2 I/O handling
#if !defined (VAR_MPI2)
                  WRITE(LUPRI,*) '  *** ERROR for *KRCICALC *** '//
     &                           ' Keyword .IJKLRO only valid in'//
     &                           ' parallel calculations'
                  WRITE(LUPRI,'(/A)') '  User input will be ignored.'
                  CINT_REORD_KRCI = .FALSE.
#endif
               GO TO 100
   15          CONTINUE
C&&& .IJKLSP: integral TT block allocation based on TTSS block
C             distribution scheme
                  CINT_SPLIT_KRCI = .TRUE.
               GO TO 100
   16          CONTINUE
C&&& PRINT : General print level in KRCI-module
                  READ(LUCMD,*) IPRKRCI
               GO TO 100
   17          CONTINUE
C&&& .FCIDUMP
                  fcidump_krci = .true.
               GO TO 100
   18          CONTINUE
C&&& .NATOCC: synonym for .NOOCCN
               GO TO 10
   19          CONTINUE
C&&& .XPSINTERFACE: interface information for Paul Bagus XPS system
                  XPSINTERFACE_KRCI = .TRUE.
               GO TO 100
   20          CONTINUE
C&&& INTFLG: Specify what two-itegrals should be included in this run
                  READ(LUCMD,*) ILLINT,ISLINT,ISSINT
               GO TO 100
   21          CONTINUE
C&&& THRPCI: threshold for printing CI vector
                  READ(LUCMD,*) KRCI_THRPCI
               GO TO 100
   22          CONTINUE
C&&& max vector block length to allocate
                  READ(LUCMD,*) krci_cfg_max_vec_block
               GO TO 100
   23          CONTINUE
C&&& MIN MK2: minimum value of 2 * MK
               READ(LUCMD,*) NKRCI_MINMK2
               IKWSET(23) = 1
               GO TO 100
   24          CONTINUE
C&&& MAX MK2: maximum value of 2 * MK
               READ(LUCMD,*) NKRCI_MAXMK2
               IKWSET(24) = 1
               GO TO 100
   25          CONTINUE
C&&& GASSH: GAS setup
                  IF (NKRCIGAS .GT. 0) THEN
                     CALL QUIT('*** ERROR in *KRCICALC *** ' //
     &                 'You cannot specify both .GASSH and .CASSCF')
                  END IF
                  READ(LUCMD,*) NKRCIGAS
                  IF (NKRCIGAS .LT. 1 .OR. NKRCIGAS .GT. MXGAS) THEN
                     CALL QUIT('*** ERROR in *KRCICALC *** ' //
     &                         'illegal value of NKRCIGAS')
                  END IF
                  DO I = 1, NKRCIGAS
                     READ(LUCMD,*) (NKRCIGSH(J,I),J=1,NFSYM)
                  END DO
               GO TO 100
   26          CONTINUE
C&&& GASSPC: GAS space constraints
                  IF (NKRCIGAS .EQ. 0) THEN
                     WRITE(LUPRI,'(A)') ' *** ERROR in *KRCICALC ***'
     &                  //' .GASSH must be specified before .GASSP!'
                     CALL QUIT('*** ERROR in *KRCICALC ***')
                  END IF
                  DO I = 1, NKRCIGAS
                     READ(LUCMD,*) (NKRCIGSP(J,I),J=1,2)
                  END DO
                  IF (NKRCIAELEC.LT.0) NKRCIAELEC = NKRCIGSP(2,NKRCIGAS)
               GO TO 100
   27          CONTINUE
C&&& MK2REFERENCE: 2 * M_K reference
                  READ(LUCMD,*) NKRCI_MK2REF
                  IKWSET(27) = 1
               GO TO 100
   28          CONTINUE
C&&& MK2DELTA: 2 * DELTA M_K
                  READ(LUCMD,*) NKRCI_MK2DEL
                  IKWSET(28) = 1
               GO TO 100
   29          CONTINUE
C&&&  CI PROGRAM: which CI program to use
                  READ(LUCMD,'(A)') KRCI_CIPROGRAM
               GO TO 100
   30          CONTINUE
C&&&  INACTIVE ORBITALS
                  READ(LUCMD,*) (NKRCIISH(I),I=1,NFSYM)
               GO TO 100
   31          CONTINUE
C&&&  CASSCF -- this defines a CAS space
               IF (NKRCIGAS .GT. 0) THEN
                  CALL QUIT('*** ERROR in *KRCICALC *** ' //
     &                 'You cannot specify both .GASSH and .CASSCF')
               END IF
                  READ(LUCMD,*) NKRCIAELEC
                  READ(LUCMD,*) (NKRCIASH(I),I=1,NFSYM)
                  NKRCIGAS = 1
                  NKRCIGSP(1,NKRCIGAS) = NKRCIAELEC
                  NKRCIGSP(2,NKRCIGAS) = NKRCIAELEC
               GO TO 100
   32          CONTINUE
C&&& .TRDM calculation of transition dipole moments and oscillator strengths
                  DOPROP_KRCI = .TRUE.
                  DOOSCILLST  = .TRUE.
               GO TO 100
   33          CONTINUE
C&&& .OMEGA quantum number: calculate <CI state|j_z|CI state> expectation value
C                          using j_z = l_z + s_z
                  DOPROP_KRCI  = .TRUE.
                  DOJZEXP = .TRUE.
               GO TO 100
C&&& .THRECI:  CI energy convergence threshold
   34          CONTINUE
               read(lucmd,*) DKRCI_THRECI
               GO TO 100
   35          CONTINUE
C&&& .CHECKP:  write intermediate check point file that is readable by
C              the MASTER - useful in parallel calculations
               CHCKPT_KRCI = .TRUE.
               GO TO 100
   36          CONTINUE
C&&& .USE PB:  use block file KRCI_BLOCKDIST.x (x=sym) (ascii format)
C              to write/get connection matrix and resulting block
C              distribution
               CUSE_PBLFL_KRCI = .TRUE.
               GO TO 100
   37          CONTINUE
C&&& .GENFOC:  transform positronic/inactive/virtual spinors to
C              Fock-type spinors
               GENFOCK_KRCI = .TRUE.
               GO TO 100
   38          CONTINUE
C&&& .TRA NO:  transform to natural CI spinors
               TRANATO_KRCI = .TRUE.
               GO TO 100
   39          CONTINUE
C&&& .LOWSRT:  use memory-saving path in the integral sorting step
C              can only be used in combination with .IJKLRO
               LOWSORT_KRCI = .TRUE.
               GO TO 100
   40          CONTINUE
C&&& .SEL ST:  select CI state(s) by <j_z> expectation value
C              the coefficients will be stored on the vector file
C              KRCI_CVECS."$sym_tag".select
C              this option may be used for "state-selection" in a
C              subsequent MCSCF...
               STATE_SELECT_KRCI = .TRUE.
               DOPROP_KRCI       = .TRUE.
               DOJZEXP           = .TRUE.
               if(.not.linear) then
                 READ(LUCMD,*) IKRCI_SYMSEL,NOMEGASEL_KRCI,
     &                         NOMEGASEL_statenr
               else
!                read mij value, root symmetry and # roots in linear symmetry input format
                 read(lucmd,*) symmetry_str1, symmetry_str2, 
     &                         NOMEGASEL_statenr

                 do i = 1,2
                 if(i == 1) symmetry_str = symmetry_str1
                 if(i == 2) symmetry_str = symmetry_str2
!                manipulate mij-value input string for further processing
!                print *,'symmetry_str is ==> ',symmetry_str
                 tmp_symmetry_str = trim(symmetry_str)
                 length_of_string = len_trim(tmp_symmetry_str)
!                print *,'length of symmetry_str is ==> ',
!    &                    length_of_string
                 last_char = tmp_symmetry_str(length_of_string:
     &                                        length_of_string)
!                print *,'last char of symmetry_str is ==> ',
!    &                    last_char

!                determine symmetry irrep of CI wave function
!                --------------------------------------------

!                a. check for inversion symmetry
                 ascii_lchar = iachar(last_char)
                 if(nfsym .gt. 1)then
                   if(ascii_lchar .eq. iachar('g'))then
                     gerade_mcci_wf = .true.
                   else if(ascii_lchar .eq. iachar('u'))then
                     gerade_mcci_wf = .false.
                   else
                     write(lupri,'(/a,/a)')
     & '  *** input error: Dinfh contains the inversion operation'//
     & ' as group element.',
     & '  please specify either g or u in addition to the mj value. ***'
                        call quit('  *** input error: Dinfh contains
     & the inversion operation as group element. ***')
                    end if
                  else
                    gerade_mcci_wf   = .true.
                    length_of_string = length_of_string + 1
                    if(ascii_lchar .eq. iachar('g')
     &                 .or. ascii_lchar .eq. iachar('u'))then
                      call quit(' *** input error: Cinfv does not
     & have inversion as group element. ***')
                    end if
                  end if

!                 b. get mij value in integer format
                  read(tmp_symmetry_str(1:length_of_string-1),
     &            fmt='(i6)') tmp_mij_val
                  fermion_sym = 1
                  if(.not.gerade_mcci_wf) fermion_sym = 2

!                 retrieve active dbg irrep (routine is inside the module symmetry_setup_krci)
                  call convert_mj_ferm_2_dbg_irrep_linsym_krci(
     &                 mj2rep,max_nkrci_max_sym,tmp_mij_val,
     &                 fermion_sym,is_current_dbg_irrep)

!                 c. store the active irrep
                  if(i == 1) IKRCI_SYMSEL   = is_current_dbg_irrep
                  if(i == 2) NOMEGASEL_KRCI = tmp_mij_val
                  end do
               end if
               IF (IKRCI_SYMSEL .LE. MAX_NKRCI_MAX_SYM .AND.
     &            IKRCI_SYMSEL .GT. 0) THEN
                  NKRCI_SELOM(IKRCI_SYMSEL) = NOMEGASEL_KRCI
               ELSE
                  WRITE(LUPRI,'(//A/I10,A,I3)')
     &            '*** ERROR in KRCI_INP .SEL ST  *** ' //
     &            'Requested symmetry outside hardwired limits:',
     &            IKRCI_SYMSEL, ' outside of',MAX_NKRCI_MAX_SYM
                  CALL QUIT('*** ERROR in KRCI_INP .SEL ST  *** ' //
     &            'Requested symmetry outside hardwired limits.')
               END IF
               GO TO 100
   41          CONTINUE
C&&& .DIPMOM:  calculate permanent dipole moment of electronic ground
C              and excited states.
               doprop_krci   = .true.
               dodipmom_krci = .true.
               GO TO 100
   42          CONTINUE
!&&& .MCANA :  tell luciarel that this is a analyze run where the CI vector is to be read from the MCSCF file
               domcana_krci  = .true.
               GO TO 100
   43          CONTINUE

C&&& .GAS SH:  the new input format in accordance with Dalton LUCITA (and soon Dirac LUCITA)
!              and Diracs opens-shell AOC Hartree-Fock format
               IF(NKRCIGAS .GT. 0)THEN
                  CALL QUIT('*** ERROR in *KRCICALC *** ' //
     &            'You cannot specify both GAS SH and '//
     &            '(.GASSH/GASSPC or .CASSCF)')
               END IF

!              step 1: read the # of GA spaces
               READ(LUCMD,*) NKRCIGAS
               IF(NKRCIGAS .LT. 1 .OR. NKRCIGAS .GT. MXGAS) THEN
                 CALL QUIT('*** ERROR in *KRCICALC *** ' //
     &                     'illegal value of NKRCIGAS')
               END IF

!              process the min max occupation / per orbital occupation in each GAS shell
               do i = 1, nkrcigas

                 read(lucmd,'(a)') input_line
                 call upcase(input_line)
                 islash = index(input_line,'/')

!                check for correct input format:
!                min e- max e- / #orbs_fermion_irrep1 #orbs_fermion_irrep2
                 if(islash <= 1)then

                   write(lupri,'(/a,i2,a,i2/a/2a)')
     &             'ERROR for *KRCICALC .GAS SHELL shell no.',
     &             i,'/',NKRCIGAS,
     &             '- the defining input line does not contain a "/":',
     &             '- the bad line : ',input_line

                   write(lupri,'(a)')
     &             '- Or the specification of the number of GAS-'//
     &             'shells might be wrong! Check the input.'
                   call quit('Input error for .GAS SHELL under'//
     &                       ' *KRCICALC')
                 end if

!                step 2: min max occupation
                 read(input_line(1:islash-1),*,iostat=ios)
     &           nkrcigsp(1,i), nkrcigsp(2,i)

                 if(ios /= 0)then
                   write(lupri,'(/a,i2,a,i2/a/2a)')
     &             'ERROR for *KRCICALC .GAS SHELL shell no.',
     &             i,'/',NKRCIGAS,
     &             '- the input line does not contain correct'//
     &             ' min max electrons','- the bad line : ',input_line
                   call quit('Input error for .GAS SHELL under '//
     &                       '*KRCICALC')
                 end if
!                step 3: GAS shell orbital occupation
                 read(input_line(islash+1:),*,iostat=ios)
     &           (NKRCIGSH(j,i), j=1,nfsym)

                 if(ios /= 0)then
                   write(lupri,'(/a,i2,a,i2/a,i2,a/2a)')
     &             'ERROR for *KRCICALC .GAS SHELL shell no.',
     &             i,'/',NKRCIGAS,
     &             '- the input line does not contain',nfsym,
     &             ' occupations','- the bad line : ',input_line
                   call quit('Input error for .GAS SHELL under '//
     &                       '*KRCICALC')
                 end if
               END DO
!              step 4: set # of active electrons ==> max occupation in last GAS space
               nkrciaelec = nkrcigsp(2,nkrcigas)

               GO TO 100
   44          CONTINUE
C&&& .SVRONO:  save reordered NOs rather than unsorted ones (the latter being the default)
               save_reordered_nos_krci = .true.
               GO TO 100
   45          CONTINUE
C&&& .SVRONX  -- save reordered NOs rather than the original NOs + active orbital offset
               save_reordered_nos_krci = .true.
               READ(LUCMD,*) (IKRCI_SVRONO(I),I=1,NFSYM)
               GO TO 100
   46          CONTINUE
C&&&  EEDM  :  calculate electron EDM effective electric field
               doprop_krci = .true.
               doeedm_krci = .true.
               GO TO 100
   47          CONTINUE
C&&&  MHYP  : calculate magnetic hyperfine structure constants
C           : alpha X electric field matrix elements
               doprop_krci = .true.
               domhyp_krci = .true.
               READ(LUCMD,*) (VKRCI_NUCSPIN(IATOM),IATOM=1,NUCIND,1)
               READ(LUCMD,*) (VKRCI_NUCMAGMOM(IATOM),IATOM=1,NUCIND,1)
               GO TO 100
CMKN The mHYP is newly implemented following eEDM
   48          CONTINUE
C&&&  ENSPS : calculate electron-Nucleus scalar-pseudoscalar P,T-odd
C           : interection constant i x gamma^0 gamma^5 \rho_N(r_e)
               doprop_krci  = .true.
               doensps_krci = .true.
               GO TO 100
CMKN The eNSPS is newly implemented following eEDM
   49          CONTINUE
C&&&  NMQM : calculate nuclear Magnetic Quadruple Moment interaction constant
C          : (\alpha X r)_z.(r_z) (i.e. \alpha X Electric Field Gradient)
               doprop_krci = .true.
               donmqm_krci = .true.
               GO TO 100
CMKN The nMQM is newly implemented following mHYP
   50          CONTINUE
               GO TO 100
   51          CONTINUE
               GO TO 100
   52          CONTINUE
C&&& TERSE  : Old LUCIAREL output for iterations and final info
               ITERSEOUT = 1
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in *KRCICALC.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in *KRCICALC')
            END IF
      END IF
  300 CONTINUE
C
      IF (.NOT. NEWDEF) GOTO 400
C
!     Warnings for options not implemented (yet)
!     ------------------------------------------
      if(lvnew)then
        write(lupri,'(/,a,/,a/)')
     & ' *** WARNING *** - (SS|SS)-integral approximation using LVNEW'//
     & ' is not implemented for KR-CI/KR-MCSCF/KR-CC.',
     & ' *** WARNING *** - Please use the LVCORR approximation.'
      end if
C
C     Consistency checks
C     ------------------

      ! symmetry: LUCIAREL does not work with C1/Ci symmetry
      if(KRCI_CIPROGRAM == 'LUCIAREL')then
        if(nfsym == 1 .and. nz == 4)then
          call quit('LUCIAREL does not work with C1/Ci symmetry.'//
     &              ' Either you use GASCIP or wait for a patch to'//
     &              ' Dirac14. For further information contact '//
     &              ' dirac-users@googlegroups.com .')
        end if 
      end if

C
C     ... integral flags
      IF(DOLVC) ISSINT = 0
      IF (LEVYLE) THEN
        ISLINT = 0
        ISSINT = 0
      ENDIF
      IF (ZORA.AND..NOT.ZORA4) THEN
        ISSINT = 0
      ENDIF
      IKRCI_INTDEF = ILLINT + 2*ISLINT + 4*ISSINT
C
C     ... calculating omega values only in linear symmetry
      if(doprop_krci.and.dojzexp)then
        if(.not.linear) then
          dojzexp = .false.
        write(lupri,'(/a)') '    *** WARNING: calculation of omega'//
     &    ' (<jz> expectation values) states turned off.***'
        write(lupri,'( a)') '    *** please specify linear symmetry'//
     &    ' in your basis set input otherwise j_z is not properly'//
     &    ' defined.'
        end if
      end if
C
C
      IF (NKRCIGAS .GT. 0) THEN
         NKRCIASHT = 0
         DO J = 1, NFSYM
            NKRCIASH(J) = 0
            DO I = 1, NKRCIGAS
               NKRCIASH(J) = NKRCIASH(J) + NKRCIGSH(J,I)
            END DO
            NKRCIASHT = NKRCIASHT + NKRCIASH(J)
         END DO
      END IF
C
C
C     Check that the number of orbitals is less than MAXASH
C     -----------------------------------------------------
C
C
      NTEST = 0
      DO J = 1, NFSYM
         NTEST = NTEST + NKRCIASH(J)
      END DO
      IF (NTEST .GT. MAXASH) THEN
        WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &       ' *** ERROR in *KRCICALC *** Too many active orbitals: ',
     &       NTEST, ' -- max is ',MAXASH
        CALL QUIT('*** ERROR in *KRCICALC : Too many active orbitals.')
      END IF
C
C     Use boson symmetry for CI expansion for the spin-free and
C     the Levy-Leblond Hamiltonian, if KRCI_UCIBOS not set in input.
C
      IF (IKWSET(22).LE.1 .AND. (SPINFR.or.levyle)) THEN
         KRCI_UCIBOS = .TRUE.
      END IF
c     krci_ucibos = .false.
C
C     Default: one CI root in symmetry 1
C
      IF (NKRCI_MAX_SYM .EQ. -1) THEN
         NKRCI_MAX_SYM    = 1
         NKRCI_CIROOTS(1) = 1
      END IF
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &   ' *KRCICALC: General set-up for KR-CI calculation:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
      WRITE(LUPRI,'(A,2I4)')
     &     ' * Inactive orbitals     : ',(NKRCIISH(I),I=1,NFSYM)
      WRITE(LUPRI,'(A,2I4)')
     &     ' * Active orbitals       : ',(NKRCIASH(I),I=1,NFSYM)
      WRITE(LUPRI,'(A,I4)')
     &     ' * Active electrons      : ',NKRCIAELEC
      IF (NKRCIGAS .GT. 0) THEN
         WRITE(LUPRI,'(A,I3,A)')
     &      ' * GAS space setup for ',NKRCIGAS,' GAS space(s) : '
         N_GAS_ERR = 0
         DO I = 1, NKRCIGAS
            WRITE(LUPRI,'(A,I3,A,2I4)')
     &           '   - GAS space ',I,'       : ',
     &           (NKRCIGSH(J,I),J=1,NFSYM)
            WRITE(LUPRI,'(A,I3,A1,I3,A)')
     &      '    (constraints: min/max active electrons after space : ',
     &           NKRCIGSP(1,I),'/',NKRCIGSP(2,I),')'
C
C           Consistency checks:
C
            IF (NKRCIGSP(1,I) .GT. NKRCIGSP(2,I))
     &         N_GAS_ERR = N_GAS_ERR + 1
            IF (I .GT. 1) THEN
               IF (NKRCIGSP(1,I-1) .GT. NKRCIGSP(1,I))
     &            N_GAS_ERR = N_GAS_ERR + 1
               IF (NKRCIGSP(2,I-1) .GT. NKRCIGSP(2,I))
     &            N_GAS_ERR = N_GAS_ERR + 1
            END IF
         END DO
         IF (N_GAS_ERR .GT. 0) THEN
            WRITE(LUPRI,'(//A)')
     &         'INPUT ERROR: The GAS specifications are not consistent'
            CALL QUIT('*** ERROR in *KRCICALC *** ' //
     &        'Inconsistency in GAS specification')
         END IF
      END IF
      WRITE(LUPRI,'(2A)')
     &     ' * CI program used       : ',KRCI_CIPROGRAM
C
C     Symmetry of the wave function:
C     ------------------------------
C
C
C     MK constraints in the CI expansion.
C
      IF (IKWSET(23) .EQ. 1 .AND. IKWSET(24) .EQ. 1) THEN
C
C        MINMK2/MAXMK2 explicitly set
C
         NKRCI_MK2REF = 0
         NKRCI_MK2DEL = NKRCI_MAXMK2
      ELSE
        IF(IKWSET(27) .EQ. 0 .AND. IKWSET(28) .EQ. 0)THEN
!         Assign default values for MINMK2/MAXMK2
          NKRCI_MAXMK2 = 2*MIN(NKRCIAELEC,NKRCIASHT) - NKRCIAELEC
          NKRCI_MINMK2 = NKRCIAELEC-2*MIN(NKRCIAELEC,NKRCIASHT)
          NKRCI_MK2REF = 0
          NKRCI_MK2DEL = NKRCI_MAXMK2
        ELSE
!         assign MINMK2/MAXMK2 defaults
          NKRCI_MINMK2 = NKRCI_MK2REF - NKRCI_MK2DEL
          NKRCI_MAXMK2 = NKRCI_MK2REF + NKRCI_MK2DEL
        END IF
      END IF
C
C     boson or fermion symmetries in the CI expansion
C
      WRITE(LUPRI,'(A)')
     & ' * Convergence of CI function(s) in the following '//
     &     'symmetries:'
      DO 202 IKRCI_SYMMETRY = 1, NKRCI_MAX_SYM
         IF( NKRCI_CIROOTS(IKRCI_SYMMETRY) .le. 0 ) GOTO 202
C     ... next symmetry
      IF ((SPINFR.or.levyle).AND.KRCI_UCIBOS) THEN
         NREP = NBSYM
         WRITE(LUPRI,'(A)')
     &     ' * Boson spatial spinor symmetry used in CI expansion.'
         WRITE(LUPRI,'(2A)')
     &     ' * Spatial symmetry (boson) of wave function : ',
     &     REP ( IKRCI_SYMMETRY - 1 )
         WRITE(LUPRI,'(A,I4,A,I3)')
     &     ' * Allowed interval of 2 * MS :',
     &     NKRCI_MINMK2, ' to ', NKRCI_MAXMK2
      ELSE IF (.NOT.(SPINFR.or.levyle).AND.KRCI_UCIBOS) THEN
         NREP = NBSYM
         WRITE(LUPRI,'(A)') ' * Approximate boson spatial spinor'//
     &     ' symmetry (as non-rel) used in CI expansion.'
         WRITE(LUPRI,'(2A)')
     &    ' * Approximate spatial symmetry (boson) of wave function : ',
     &     REP ( IKRCI_SYMMETRY - 1 )
         WRITE(LUPRI,'(A,I4,A,I3)')
     &     ' * Allowed interval of approximate 2 * MS :',
     &     NKRCI_MINMK2, ' to ', NKRCI_MAXMK2
      ELSE
         if(linear)then
!          store on common block in krciprop.h
           xrepeig(ikrci_symmetry) =
     &     nkrci_ciroots_symstr(ikrci_symmetry)
           if(mod(nkrciaelec, 2) .eq. 0)then
             ferm_lab = '  '
           else
             ferm_lab = '/2'
           end if
           WRITE(LUPRI,'(a,i3,a,a4,a2)')
     &        '    ** ',nkrci_ciroots(ikrci_symmetry),' eigenstate(s)'//
     &        ' for MJ-value (doubled): ',
     &          adjustr(xrepeig(ikrci_symmetry)),ferm_lab
         else
           allocate(multb_tmp(64,64))
           IF (SPINFR.or.levyle) THEN
              CALL GMULTSF(NREP,REPNA,MULTB_TMP)
           ELSE
              CALL GMULTA(NZ,NFSYM,NREP,REPNA,MULTB_TMP)
           ENDIF
           deallocate(multb_tmp)
           IF(MOD(NKRCIAELEC, 2) .EQ. 0)THEN
             IRRP = NREP + IKRCI_SYMMETRY
             WRITE(LUPRI,'(A,I3,2A)')
     &       '    ** ',NKRCI_CIROOTS(IKRCI_SYMMETRY),' eigenstate(s)'//
     &       ' in symmetry (boson):',REPNA(IRRP)
C            store on common block in krciprop.h
             XREPEIG(IKRCI_SYMMETRY) = REPNA(IRRP)
             IF (NKRCI_MAX_SYM .GT. NREP) THEN
                WRITE(LUPRI,'(/A//A,2I5)')
     &          ' INPUT ERROR, requested symmetry (boson) in .CIROOTS'//
     &          ' is non-existent!',
     &          ' Symmetry, max value of symmetry = ',
     &            NKRCI_MAX_SYM, NREP
                CALL QUIT('INPUT ERROR for *KRCICALC .CIROOTS')
             END IF
C
           ELSE
              WRITE(LUPRI,'(A,I3,2A)')
     &        '    ** ',NKRCI_CIROOTS(IKRCI_SYMMETRY),' eigenstate(s)'//
     &        ' in symmetry (fermion):',REPNA(IKRCI_SYMMETRY)
C               store on common block in krciprop.h
                XREPEIG(IKRCI_SYMMETRY) = REPNA(IKRCI_SYMMETRY)
C
C           FIXME: check proper abort criterion
           END IF
         end if
         WRITE(LUPRI,'(A,I4,A,I3)')
     &     '    -- Allowed interval of 2 * MK :',
     &     NKRCI_MINMK2,' to ',NKRCI_MAXMK2
      END IF
 202  CONTINUE
C
      IF(KRCI_CIPROGRAM .eq. 'LUCIAREL') THEN
         WRITE(LUPRI,'(A/A)')
     &' * Using symmetry nomenclature for LUCIAREL. ',
     &'      Boson and Fermion irreps of complex (sub)groups  '
      END IF
C
      WRITE(LUPRI,'(/A)')
     +    ' * Contributions from 2-electron integrals to Fock matrix:'
      IF(LBIT(IKRCI_INTDEF,1)) WRITE(LUPRI,'(A)') '   LL-integrals.'
      IF(LBIT(IKRCI_INTDEF,2)) WRITE(LUPRI,'(A)') '   SL-integrals.'
      IF(LBIT(IKRCI_INTDEF,3)) WRITE(LUPRI,'(A)') '   SS-integrals.'
      IF(LBIT(IKRCI_INTDEF,4)) WRITE(LUPRI,'(A)') '   GT-integrals.'
      WRITE(LUPRI,'(A,I4)') ' * General print level   : ',IPRKRCI
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &   ' Control parameters for KR-CI optimization'
      CALL PRSYMB(LUPRI,'=',75,0)
C
      WRITE(LUPRI,'(A,1P,D10.2)')
     &     ' * Energy convergence threshold: ', DKRCI_THRECI
      IF (MAXCIT_KRCI .GE. 0) THEN
         WRITE(LUPRI,'(A,I3)')
     &     ' * Maximum number of CI iterations for each symmetry: ',
     &     MAXCIT_KRCI
      ELSE
         WRITE(LUPRI,'(A)') ' * Full CI '
      END IF
C
      IF (IRESTRK_KRCI.EQ.1) THEN
         WRITE(LUPRI,'(A)')
     &     ' * Restarting CI from C vector on LUCIAREL file '
      END IF
C
      IF (IANACI.EQ.1) THEN
         WRITE(LUPRI,'(A)')
     &     ' * Analyzing CI vectors in terms of coefficients '
      END IF
C
      IF (MXCIV_KRCI .gt.0 ) WRITE(LUPRI,'(A,I3)')
     &     ' * Maximum subspace dimension set to ',MXCIV_KRCI
C
      IF( IKRCI_MEMFAC .ne. 9 ) WRITE (LUPRI, '(A,I3)' )
     &     ' * Multiplier for memory guess: ',IKRCI_MEMFAC

      IF( IDOPARIO_KRCI .eq. 1 )THEN
        WRITE(LUPRI,'(A)')
     &     ' * Integrals also for slaves assumed to be in place'
        WRITE(LUPRI,'(A)')
     &     '   No integral broadcast from MASTER to slaves.'
      ELSE
#if defined (VAR_MPI)
        WRITE(LUPRI,'(A)')
     &     ' * Integrals on slave nodes provided by the MASTER'
#endif
      END IF
C     no shared memory mode for non MPI-2 executables
#if !defined (VAR_MPI2)
      IF( CSHMEMO_KRCI ) CSHMEMO_KRCI = .FALSE.
#endif
      IF( CSHMEMO_KRCI )THEN
        WRITE(LUPRI,'(A)')
     &     ' * Running in parallel MPI-2 "shared memory" mode'
        IF( ISHMEM_TYPE .eq. 0) THEN
          WRITE(LUPRI,'(A)')
     &     ' --> chosen level 0: disabling shared memory usage'
        ELSE IF( ISHMEM_TYPE .eq. 1) THEN
          WRITE(LUPRI,'(A)')
     &     ' --> chosen level 1: sharing node-wise (ij|kl) '
        ELSE IF(  ISHMEM_TYPE .eq. 2 ) THEN
          WRITE(LUPRI,'(A)')
     &     ' --> chosen level 2: sharing node-wise |C> vector'
        ELSE IF(  ISHMEM_TYPE .eq. 3 ) THEN
          WRITE(LUPRI,'(A)')
     &     ' --> chosen level 3: sharing node-wise (ij|kl) + |C> vector'
        ELSE IF(  ISHMEM_TYPE .eq. 4 ) THEN
          WRITE(LUPRI,'(A)')
     &     ' --> chosen level 4: sharing global (ij|kl) + node-wise |C>'
        ELSE IF(  ISHMEM_TYPE .eq. 5 ) THEN
          WRITE(LUPRI,'(A)')
     &     ' --> chosen level 5: sharing global |C> + node-wise (ij|kl)'
        ELSE IF(  ISHMEM_TYPE .eq. 6 ) THEN
          WRITE(LUPRI,'(A)')
     &     ' --> chosen level 6: sharing global (ij|kl) + |C> vector  '
        END IF
      END IF
      IF( CINT_REORD_KRCI )THEN
        WRITE(LUPRI,'(A)')
     &     ' * GAS-scheme specific integral file in use (.IJKLRO)'
      END IF
      IF( CINT_SPLIT_KRCI )THEN
        WRITE(LUPRI,'(A)')
     &     ' * Revised integral allocation scheme (.IJKLSPLIT)'
      END IF
C     ... only natural orbitals occupation numbers
      IF( NATOLCR_KRCI )THEN
        WRITE(LUPRI,'(A)')
     &  ' * Calculation of natural orbital occupation numbers (.NOOCCN)'
      END IF
!     additional check because we might have turned off the <j_z>
!     calculation. see comment above for more information.
      if(dooscillst.or.dodipmom_krci.or.dojzexp) doprop_krci = .true.

      IF( DOPROP_KRCI ) THEN
C       ... obtain operator symmetry in LUCIARELs double group
C       subgroup notation
        DOSYMOPPRP = .TRUE.
        WRITE(LUPRI,'(A/)')
     &     ' * Property calculation for the following'//
     &     ' one-electron operators: '
        if(dooscillst.or.dodipmom_krci)call def_oscill(0)
        if(dojzexp)                    call def_omegaquant(0)
        if(doeedm_krci)                call def_eedm(iprkrci,work,lwork)
        if(domhyp_krci)                call def_mhyp(iprkrci,work,lwork)
        if(doensps_krci)              call def_ensps(iprkrci,work,lwork)
        if(donmqm_krci)                call def_nmqm(iprkrci,work,lwork) !MKN
        DO I = 1, NPROP_KRCI
C         ... initialize operator symmetry array
          ISYMOPPRP_KRCI(I) = 0
C
          INDXPR = LPROP_KRCI(I)
          CALL WRIXPR(I,INDXPR)
        ENDDO
C
      ENDIF
      IF(DOOSCILLST) THEN
        WRITE(LUPRI,'(/A)') ' * Calculation of oscillator strengths'//
     &                     ' for (allowed) electronic transitions'
      END IF
      if(dodipmom_krci)then
        write(lupri,'(/a)') ' * Calculation of permanent dipole'//
     &                     'moments of electronic states'
      end if
      if(doeedm_krci)then
        write(lupri,'(/a)') ' * Calculation of electron EDM'//
     &                     ' electric-field matrix elements'
      end if
      if(domhyp_krci)then
        write(lupri,'(/a)') ' * Calculation of Hyperfine Constants'//
     &                     ' Hyperfine Constants  matrix elements'
      end if
      if(doensps_krci)then
        write(lupri,'(/a)') ' * Calculation of electron-Nucleus'//
     &                     ' Scalar-Pseudoscalar matrix elements'
      end if
      if(donmqm_krci)then
        write(lupri,'(/a)') ' * Calculation of the Nuclear MQM'//
     &                     ' P,T-odd interaction matrix elements' !MKN
      end if
      IF(DOJZEXP) THEN
        WRITE(LUPRI,'(/A)') ' * Determine Omega quantum number for'//
     &                      ' electronic eigenstates'
      END IF
      IF( CHCKPT_KRCI ) THEN
        WRITE(LUPRI,'(/A)') ' * Write check point file in parallel'//
     &                      ' CI calculation'
      END IF
      IF( GENFOCK_KRCI ) THEN
        WRITE(LUPRI,'(/A)') ' * Transform positronic/inactive/virtual'//
     &                      ' spinors to Fock-type spinors'
      END IF
      IF( (.NOT. CINT_REORD_KRCI) .and. LOWSORT_KRCI ) THEN
        WRITE(LUPRI,'(/A)') ' *** Warning *** low-memory sorting'//
     &                      ' scheme only available in connection'//
     &                      ' with .IJKLRO'
        LOWSORT_KRCI = .FALSE.
      END IF
      IF( LOWSORT_KRCI ) THEN
        WRITE(LUPRI,'(/A)') ' * Use memory-saving integral reordering'//
     &                      ' scheme (.LOWSRT)'
      END IF
      if(state_select_krci)then
        write(lupri,'(/a)') ' * Select CI state(s) for subsequent'//
     &                      ' MCSCF optimizations based on the'//
     &                      ' respective Omega value(s).'
        do 303 ikrci_symmetry = 1, max_nkrci_max_sym
          if(nkrci_selom(ikrci_symmetry) .gt. 0)then
            if(nkrci_ciroots(ikrci_symmetry) .le. 0)then
              call quit(' *** error in KRCI_INP: ' //
     &                  ' CI root selection for inactive symmetry. ***')
            end if
            if(nkrci_ciroots(ikrci_symmetry) .lt.
     &         NOMEGASEL_statenr)then
              call quit(' *** error in KRCI_INP: '//
     &                  ' CI root selection for state > # of roots in'//
     &                  ' this symmetry. ***')
            end if
            write(lupri,'(/a,i4,a,a5,a)')
     &      ' *** symmetry:',ikrci_symmetry,
     &      ' Omega value = ',xrepeig(ikrci_symmetry),
     &      '(/2 for fermionic systems)'
          end if
  303   continue
      end if
C
      IF ( XPSINTERFACE_KRCI ) THEN
         WRITE (LUPRI,'(/A)')
     &   " * Writing interface file to Paul Bagus' XPS programs"
      END IF
C
  400 if(linear) then
        if(allocated(nkrci_ciroots_symstr))
     &    deallocate(nkrci_ciroots_symstr)
      end if

      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRCI_CALC(WFTYP,WFCONV)
C***********************************************************************
C
C     Driver routine to KRCI_CALC_1
C     Get CI wave function.
C     Either LUCIAREL or GASCIP or ...
C
C     Based on PSIOPT.
C
C     Written by S. Knecht and H. J. Aa. Jensen - Aug 2008
C
C     Last revision :
C
C***********************************************************************
      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
      PARAMETER(D0=0.0D0)
C
      CHARACTER WFTYP*(*), KRCLCTYP*8
      LOGICAL   WFCONV
C
#include "maxorb.h"
#include "dcbidx.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbopt.h"
#include "dcbgen.h"
#include "dcbkrci.h"
#include "krciprop.h"
#include "dcbham.h"
#include "frame.h"
!frame.h needed for CORRLV for LVCORR
C
      LOGICAL EX, FILEOPEN
      INTEGER RENAME
C
      real(8), allocatable :: WORK(:)

      CALL QENTER('KRCI_CALC')

      call legacy_lwork_get(LWORK)
#include "memint.h"
C
      call alloc(WORK,LWORK,id='WORK in KRCI_CALC')
C
C
C     Determine the wave function type
C     --------------------------------
C
      IF ( WFTYP .EQ. 'SCF  ' .or. WFTYP .EQ. 'MCSCF') THEN
C
C        that should not happen...
C
         WRITE(LUPRI,'(/A)')
     $   '*** ERROR in KRCI_CALC *** wrong wave function type '
         CALL QUIT(
     &   '*** ERROR in KRCI_CALC: wave function type not allowed ***')
C
      ELSE IF ( WFTYP .EQ. 'KR-CI ') THEN
        MCTYPE = JGAS
      END IF
C
C     Say Hi!
C     -------
C
      CALL KRCIHI
C
C     ***************************************************
C     *** Open files needed *****************************
C     ***************************************************
C     *** Set logical units ***
C
      CALL SETKRCI(LUPRI)
      CALL SETDCBIDX()
#if defined (VAR_MPI2)
C
C     ***************************************************
C     *** Initialize common block information         ***
C     ***************************************************
C
!     1. mpi
      CALL MPIXCALC_TYPESZ()
#endif
!     2. lvcorr correction
      if(dolvc.and.((iopt_intdef.eq.3).or.(onecap.and.intv1c.eq.2)))then
        kfree_save = kfree
        call memget('REAL',klvc,n2bbasxq*(nopen_mc+1),work,kfree,lfree)
        call lvcorr(work(klvc),work(kfree),lfree,ipropt)
        call memrel('lvcorr-ci',work,1,kfree_save,kfree,lfree)
!       set lvcorr value (otherwise it was initialized to 0.0d0)
        E_lvcorr_mc = corrlv
      else
        E_lvcorr_mc = 0.0d0
      end if

!
!     set convergence threshold for GASCIP CI calculations
      opt_thrgrd = DKRCI_THRGCI

      IF ( XPSINTERFACE_KRCI ) THEN
      ! write headers to XPS interface files
         CALL QUIT('Sorry, .XPSINTERFACE is not implemented yet')
      END IF
C
C     ***************************************************
C     *** Loop over symmetries for CI roots *************
C     ***************************************************
C
      CALL MEMGET('REAL',   KCMO  ,NCMOTQ  ,WORK,KFREE,LFREE)
      KFREE_100_LOOP = KFREE

      DO 100 IOPT_SYMMETRY = 1, NKRCI_MAX_SYM
        NCIROOT = NKRCI_CIROOTS(IOPT_SYMMETRY)
        IF (NCIROOT .LE. 0) GO TO 100
        WRITE (LUPRI,'(//80A1/A,I5/A,I8/80A1/)')
     &  ('&',I=1,80),
     &  '&&& KRCI calculation for symmetry no.',IOPT_SYMMETRY,
     &  '&&& Number of CI roots for this symmetry',NCIROOT,
     &  ('&',I=1,80)

C
C       ***************************************************
C       *** Initialize configurational and orbital      ***
C       *** parameters.                                 ***
C       ***************************************************
C
        CALL RSETWOP(WFTYP,WORK,KFREE,LFREE)
        CALL FLSHFO(LUPRI)
        IF (NZCONF .LE. 0) THEN
           WRITE(LUPRI,'(/A)')
     &      'INFO: CI skipped because no determinants in this symmetry'
           GO TO 100
        ELSE IF (NZCONF .LT. NCIROOT) THEN
           WRITE(LUPRI,'(/A,I0,A/A,I0)')
     &        'INFO: Only ',NZCONF,' determinants in this symmetry',
     %        'INFO: Number of CI roots therefore reduced to ',NZCONF
           NCIROOT = NZCONF
        END IF
C
C       ***************************************************
C       *** Memory allocation *****************************
C       ***************************************************
C
        IF (OPT_CIPROGRAM .EQ. 'LUCIAREL') THEN
           CALL MEMGET2('REAL','CVECS',KCVECS,0, WORK,KFREE,LFREE)
        ELSE
           CALL MEMGET2('REAL','CVECS',KCVECS,NCIROOT*NZCONFQ,
     &        WORK,KFREE,LFREE)
        END IF
!       disable all actions related to a spinfree MCSCF run
        spinfr_krmc = .false.
        call memget('INTE',kibeig,    0,work,kfree,lfree)
C
C       ***************************************************
C       *** Call DRIVER for CI calculation ****************
C       ***************************************************
C
        CALL KRCI_CALC_1(WORK(KCMO),WORK(KCVECS),WORK(KIBEIG),
     &                   WFTYP,WFCONV,WORK,KFREE,LFREE)
C
        CALL MEMREL('KRCI_CALC_1',WORK,KWORK,KFREE_100_LOOP,
     &              KFREE,LFREE)
  100 CONTINUE
C     ^-- end loop IOPT_SYMMETRY over CI state symmetries
C
C     ***************************************************
C     *** KR-CI property calculation ********************
C     ***************************************************
C
      IF( DOPROP_KRCI ) THEN
C
C        check symmetry with non-zero CI roots which
C        will also have NZCONF != 0.
C        IOPT_SYMMETRY = 1  ! this fails for fermionic systems
C                             (ungerade number of electrons)
         DO ISYM_CHECK = 1, NKRCI_MAX_SYM
            NCIROOT_CHECK = NKRCI_CIROOTS(ISYM_CHECK)
            IF( NCIROOT_CHECK .ne. 0 ) IOPT_SYMMETRY = ISYM_CHECK
         END DO
C        set run flag
         RUNXPROP = .TRUE.
         CALL RSETWOP(WFTYP,WORK,KFREE,LFREE)
         CALL KRCI_PROP(WORK(KCMO),WORK(KFREE),LFREE)
      END IF
C
C     ***************************************************
C     *** Close files ***********************************
C     ***************************************************
C
      CLOSE(LUKRM3,STATUS='KEEP')
      IF (LUKRM5.GT.0) CLOSE(LUKRM5,STATUS='KEEP')
C
C     Save orbitals on file KRMCSCF if it does not exist.
C     Other DIRAC modules may expect to find it, e.g., MULLPOP.
C
      LUNI = -1
      INQUIRE(FILE='KRMCSCF',EXIST=EX,OPENED=FILEOPEN,NUMBER=LUNI)
      IF(EX.and.FILEOPEN)THEN
        CALL WRTKRMC(LUNI,'NEWORB  ',WORK(KCMO),NCMOTQ)
      ELSE IF( EX.and.(.not. FILEOPEN))THEN
        LUNI = 99
        CALL OPNFIL(LUNI,'KRMCSCF','OLD    ','KR CI ')
        CALL WRTKRMC(LUNI,'NEWORB  ',WORK(KCMO),NCMOTQ)
        CLOSE(LUNI,STATUS='KEEP')
      ELSE IF(.NOT. EX)THEN
        LUNI = 99
        CALL OPNFIL(LUNI,'KRMCSCF','UNKNOWN','KR CI ')
        CALL NEWLAB('*KRCI   ',LUNI,LUPRI)
        CALL WRTKRMC(LUNI,'NEWORB  ',WORK(KCMO),NCMOTQ)
        CLOSE(LUNI,STATUS='KEEP')
      END IF
C
C
C     ***************************************************
C     *** Memory deallocation ***************************
C     ***************************************************
C
      CALL MEMREL('KRCI_CALC',WORK,KWORK,KWORK,KFREE,LFREE)
      call dealloc(WORK)
C
      CALL QEXIT('KRCI_CALC')
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRCI_CALC_1(CMO,CVECS,IBEIG,WFTYP,WFCONV,WORK,
     &                       KFREE,LFREE)
C***********************************************************************
C
C     Get CI wave function.
C
C     Based on PSIOP1.
C
C     Written by S. Knecht and H. J. Aa. Jensen - Aug 2008
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "maxorb.h"
#include "dcbpsi.h"
#include "dcbkrmc_itinfo.h"
#include "dcbgen.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "frame.h"
#include "krciprop.h"
C
      LOGICAL   WFCONV, DO4INDEX, TOBE
      DIMENSION WORK(*), CMO(*), CVECS(*), IBEIG(*)
      CHARACTER WFTYP*(*), SECTID*12, CPUTID*12, WALLTID*12
C
      CALL QENTER('KRCI_CALC_1')
C
      CALL GETTIM(CPUTOT1,WALLTOT1)
      KFRSAV = KFREE
C
C
C     ***************************************************
C     *** Initialize some parameters ********************
C     ***************************************************
C
      DO I = 1,2
         IOPT_ISYMOP(I) = 1
         IOPT_IFCKOP(I) = 1
         IOPT_IHRMOP(I) = 1
      END DO
      JTRLVL = 0
      DO4INDEX = ( MCTYPE .GE. JMCMIN )
C     ... MCTYPE .lt. JMCMIN is SCF
      EMCOLD = D0
      WFCONV = .FALSE.
C
C     Integrals for KR-CI calculation.
C
      CALL DZERO(GNORM,5)
      ITMAC  = 1
      CALL INTCON(IOPT_INTFLG,IOPT_INTBUF,IOPT_INTDEF,GNORM(5),
     &            OPT_CNVINT,ITMAC,IOPT_ITRINT,OPT_INTTYP)
C
C
C     ***************************************************
C     *** Get CI vectors ********************************
C     ***************************************************
C
C
#ifdef UNDEF
      write(LUPRI,*) '(1) iopt_intflg=',iopt_intflg
      write(LUPRI,*) '(1) iopt_intbuf=',iopt_intflg
      write(LUPRI,*) '(1) iopt_intdef=',iopt_intflg
      write(LUPRI,*) '(1) gnorm(5)=',gnorm(5)
      write(LUPRI,*) '(1) opt_cnvint=',opt_cnvint
      write(LUPRI,*) '(1) itmac=',itmac
      write(LUPRI,*) '(1) iopt_itrint=',iopt_itrint
      write(LUPRI,*) '(1) opt_inttyp=',opt_inttyp
#endif
C
C     Read start orbitals and get CI vectors:
C
C     Do not write orbitals to file KRMCSCF in ROPTST.
C     KRMCSCF may not be opened.
      IWRT = 0

      CALL ROPTST(IOPTST,IWRT,CMO,CVECS,IBEIG,IPROPT,WORK,KFREE,LFREE)
C     CALL ROPTST(ISTART,IWRT,CMO,CREF ,IBEIG,IPRINT,WORK,KFREE,LFREE)
C
      IF (KTRLVL.eq.5) GOTO 3000
C
C
C     ***************************************************
C     *** KR-CI finished - no task left *****************
C     ***************************************************
C
C
 3000 CONTINUE
C
      WFCONV = .TRUE.
C
      CALL GETTIM(CPUTOT2,WALLTOT2)
      CPUTID  = SECTID(CPUTOT2-CPUTOT1)
      WALLTID = SECTID(WALLTOT2-WALLTOT1)
      WRITE(LUPRI,9300) CPUTID,WALLTID
C
 9300 FORMAT(/'>>> TOTAL CPU (WALL) TIME IN KR-CI: ',A,'(',A,')'///)
C
      CALL MEMREL('KRCI_CALC_1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('KRCI_CALC_1')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE KRCIHI
C***********************************************************************
C
C     Write welcome message.
C
C     Written by J. Thyssen - Jan 18 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "maxorb.h"
#include "dcbidx.h"
#include "dcborb.h"
#include "dcbopt.h"
C
      CALL TITLER('KR-CI calculation ','*',199)
C
      WRITE(LUPRI,9000)
 9000 FORMAT(
     &     3X,'This is output from DIRAC KR-CI',/,
     &     3X,'- a relativistic four-component CI',
     &     1X,'wave function program.',
     &     ///,
     &     3X,'General structure: ',/,
     &     3X,'  Stefan Knecht and Hans Joergen Aa. Jensen ',//,
     &     3X,'Integral transformation: ',/,
     &     3X,'  Luuk Visscher, Jon K. Laerdahl, and Trond Saue',//,
     &     3X,'GASCIP CI code: ',/,
     &     3X,'  Joern Thyssen and Hans Joergen Aa. Jensen',//,
     &     3X,'LUCIAREL CI code: ',/,
     &     3X,'  Timo Fleig and Jeppe Olsen',//,
     &     3X,'Parallel LUCIAREL CI code: ',/,
     &     3X,'  Stefan Knecht, Hans Joergen Aa. Jensen and Timo Fleig',
     &   //3X,'Linear symmetry implementation (CI and MCSCF):',/,
     &     3X,'  Stefan Knecht and Hans Joergen Aa. Jensen',
     &     //,
     &     1X,79('*'))
!
      WRITE(LUPRI,9110)
 9110 FORMAT(/' This module is published in:',//,
     &     4X,'GASCIP: J Thyssen, T Fleig, and H J Aa Jensen',/,
     &    20X,'   J Chem Phys 129, 034109 (2008), suppl. material.',/,
     &     4X,'DIRAC-LUCIAREL: T Fleig, J Olsen, and L Visscher',/,
     &    20X,'   J Chem Phys, 119,6 (2003) 2963',/,
     &     4X,'PARALLEL LUCIAREL:',/,
     &     4X,'                S. Knecht, H J Aa Jensen, and T Fleig',/,
     &    20X,'   J Chem Phys, 132,1 (2010) 014108', //,
     &     1X,79('*'),/)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SETKRCI(LUPRI)
C***********************************************************************
C
C     Assign logical units and open files.
C
C     Based on SETKRMC.
C
C     Written by S. Knecht - Aug 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcbpsi.h"
C
      LOGICAL EX
C
C     LUFCK1 definition needs to be done here
C     as SETDHF is not called in the KRCI run:
C
      LUFCK1 = 10
C
C     Units 50-59 are reserved for KRMC and KRCI.
C     Units 50-56 are opened here, 57-58 in krmctra.F.
C     (hjaaj March 2002).
C     Unit 59 is currently used for KRMCOLD, restarting option.
C     (tf August 2003)
C
C     LUKRMC ('KRMCSCF') not used in KRCI, but not touched here
C            because LUKRMC = 50 is set globally in dirac/dirset.F
C     LUKRM1 ('KRMC_KAPPA') not used in KRCI
C     LUKRM2 ('KRMC_CMO') not used in KRCI
C     LUKRM3 ('KRMC_FOCK') used in KRCI for FOCK matrix and property
C            matrices
C     LUKRM4 ('KRMC_LRED') not used in KRCI
C     LUITFO ('KRMC_ITINFO') not used in KRCI
C
      LUKRM1 = -51
      LUKRM2 = -52
      LUKRM3 =  53 ! note the invisible plus
      LUKRM4 = -54
      LUITFO = -56
C
C     LUKRM5 is used for reading KR-MCSCF orbitals in ROPTST.
C     First check if KRMCSCF available, meaning this calculation follows
C     a KR-MCSCF calculation. Then we will use these orbitals. If
C     KRMCSCF is not available, we check for 'KRMCOLD' with orbitals
C     from a previous calculation.
C
      LUKRM5 = 59
      INQUIRE ( FILE = 'KRMCSCF', EXIST = EX )
      IF( EX ) THEN
         CALL OPNFIL(LUKRM5,'KRMCSCF','OLD','STKRCI')
         REWIND LUKRM5
      ELSE
         INQUIRE ( FILE = 'KRMCOLD', EXIST = EX )
         IF (EX) THEN
            CALL OPNFIL(LUKRM5,'KRMCOLD','OLD','STKRCI')
            REWIND LUKRM5
         ELSE
            LUKRM5 = -59
         END IF
      END IF
C
C     Delete 2-el integral files from KRMCSCF or from MOLTRA, they
C     are probably not correct for KRCI. However, if neither KRMCSCF or
C     MOLTRA have been active in this Dirac run, then we assume we are
C     restarting KRCI and that any integrals are correct. (aug08, HJAaJ)
C
      IF (DOKRMC .OR. DOTRA) THEN
         CALL OPNFIL(LUKRM3,'4INDINFO','UNKNOWN','STKRCI')
         CLOSE(LUKRM3,STATUS='DELETE')
         CALL OPNFIL(LUKRM3,'KRMC_FOCK','UNKNOWN','STKRCI')
         CLOSE(LUKRM3,STATUS='DELETE')
      END IF
C
C     If old 'KRMC_FOCK' still exists, then we are restarting a KRCI.
C
      INQUIRE ( FILE = 'KRMC_FOCK', EXIST = EX )
      IF( EX ) THEN
        CALL OPNFIL(LUKRM3,'KRMC_FOCK','OLD','STKRCI')
      ELSE
        CALL OPNFIL(LUKRM3,'KRMC_FOCK','UNKNOWN','STKRCI')
        CALL NEWLAB('SODLABEL',LUKRM3,LUPRI)
      END IF
      REWIND LUKRM3
C
      RETURN
      END

      SUBROUTINE fcidump_driver(CMO,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Generate FCI dump file
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION WORK(*), CMO(*)
C
      CALL QENTER('fcidump_driver')
      CALL GETTIM(CPU0,WALL0)
C
      KFRSAV = KFREE
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
C
C     Transform quaternion FCACmo to molfdir type FCAC.
C
      CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPRINT)
C
      CALL MEMREL('fcidump_driver after rGETFC',WORK,1,KFCAC,
     &            KFREE,LFREE)
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

      CALL DNZ32M(WORK(KH2AC),WORK(KMUUUU),IPRINT)

      call quit('DMRG interface not available in this Dirac version')
C
      CALL MEMREL('fcidump_driver',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL GETTIM(CPU2,WALL2)
      WRITE (LUPRI,'(/A,2F20.2)')
     &   "CPU and WALL times for generating the FCIDUMP file:",
     &    CPU2-CPU0,WALL2-WALL0
      CALL FLSHFO(LUPRI)
      CALL QEXIT('fcidump_driver')
C
      END
C=======================================================================
C     end of file krci_ctl.F
C=======================================================================
