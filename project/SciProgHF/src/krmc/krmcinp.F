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

!***********************************************************************

      SUBROUTINE KRMCINP(WORD,RESET,WORK,LWORK)
!***********************************************************************
!
!     Input section for KRMSCF-module
!
!     Written by J. Thyssen - Oct 26 1998
!     modified by Hans Joergen Aa. Jensen, Timo Fleig and Stefan Knecht
!
!***********************************************************************
      use symmetry_setup_krci
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxash.h"
!
      PARAMETER (NTABLE = 36)
      PARAMETER (RTOL = 1.0D-15,D1 = 1.0D0,D0=0.0D0,D2=2.0D00)
!
      LOGICAL SET, NEWDEF, LBIT, RESET
      CHARACTER*4 REPNA(64)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER CTEMP*72, LINE*80, TEXT*20, REPA(8)*4
      DIMENSION WORK(LWORK)
!
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbkrmc.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "pgroup.h"
#include "krciprop.h"
C
      character (len= 4)              :: symmetry_str
      character (len= 4)              :: tmp_symmetry_str
      character (len= 4), allocatable :: mcscf_ciroots_symstr(:)
      character (len= 1)              :: last_char
      character (len= 2)              :: ferm_lab
      logical                         :: gerade_mcscf_wf
      integer                         :: tmp_mij_val
      integer                         :: fermion_sym
      integer                         :: ascii_lchar
      integer                         :: double_group
      integer, allocatable            :: multb_tmp(:,:)
      integer, parameter              :: max_nkrmc_max_sym = 128
      character(len=400)              :: input_line
      character(len=400)              :: test_string
      integer                         :: ios, islash, icomment

      SAVE SET
      DATA TABLE /'.THRESH','.PRINT ','.MAX MA','.MAX MI','.NR ALW',
     &            '.NEO AL','.STATE ','.DELETE','.MAX BA','.MAKE C',
     &            '.SKIPEE','.SKIPEP','.CNVINT','.ITRINT','.INTFLG',
     &            '.NO CI ','.THRPCI','.FROZEN','.MIN MK','.MAX MK',
     &            '.GASSH ','.GASSPC','.MK2REF','.MK2DEL','.FREEZE',
     &            '.CI PRO','.INACTI','.CASSCF','.GAS SH','.SYMMET',
     &            '.NO1pDE','.MEMFAC','.NOCKJZ','.SVRONO','.WITHEP',
     &            '.MVOFAC'/
      DATA SET/.FALSE./
C
C
#include "ibtfun.h"
C
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
C     Initialize /CBLKRMC/
C     ===================
C
      DO I = 1,NMCFLAG
         MCFLAG(I)            = .FALSE.
      END DO
      KRMC_SKIPEE             = .FALSE.
      KRMC_SKIPEP             = .TRUE.  ! 31-Jan-17, new default /hjaaj
      KRMC_NOCI               = .FALSE.
      KRMC_UCIBOS             = .FALSE.
      KRMC_CHCKJZ             = .FALSE.
      IF(LINEAR) KRMC_CHCKJZ  = .TRUE.
      KRMC_no1pdens           = .FALSE.
      krmc_save_reordered_nos = .false.
!     initialize logical flag for new full-ci input 
      KRMC_full_ci  = .false.
C
C     Initialize /CBIKRMC/
C     ===================
C
      IPRKRMC   = IPRGEN
      MKRMC_MXMICRO   = 50
      MKRMC_MXMACRO   = 25
      MKRMC_MAXBCK    =  5
      ILLINT          = IBTAND(INTGEN,1)
      ISLINT          = IBTAND(INTGEN/2,1)
      ISSINT          = IBTAND(INTGEN/4,1)
      IKRMC_INTBUF    =  0
      IKRMC_INTFLG    =  0
      IKRMC_ITRINT(1) =  0
      IKRMC_ITRINT(2) =  0
      IKRMC_STATE     =  1 ! ground state
      IKRMC_SYMMETRY  =  1 ! totally symmetric
      IKRMC_MEMFAC    =  9
      NKRMCGAS        =  0
      NKRMC_MK2REF    = -1
      NKRMC_MK2DEL    = -1
      NKRMC_MINMK2    = -1
      NKRMC_MAXMK2    = -1
      IKRMC_SVRONO(1) =  0
      IKRMC_SVRONO(2) =  0
C
C     Initialize /CBRKRMC/
C     ===================
C
      DKRMC_CNVINT(1) = DUMMY
      DKRMC_CNVINT(2) = DUMMY
      DKRMC_THRGRD    = 0.5D-5
      DKRMC_THRPCI    = 1.0D-3
      DKRMC_MVOFAC    = 1.0D0
C
C     Initialize /DCOKRMC/
C     ===================
C
      DO I = 1,2
         NKRMCISH(I) = 0
         NKRMCASH(I) = 0
         NKRMCSSH(I) = 0
         NKRMCPSH(I) = 0
         NKRMCFRO(I) = 0
         DO J = 1, MXGAS
            NKRMCGSH(I,J) = 0
         END DO
      END DO
      NKRMCAELEC = -1
C
C     Initialize /DCCKRMC/
C     ===================
C
      KRMC_CIPROGRAM = 'GASCIP'
      KRMC_FRZSTR(1) = ' '
      KRMC_FRZSTR(2) = ' '
      KRMC_DELSTR(1) = ' '
      KRMC_DELSTR(2) = ' '
C
C     Process input for KRMCS
C     ========================
C
      NEWDEF = (WORD .EQ. '*KRMCSC')
      ICHANG = 0
      IF(NEWDEF)THEN

!       initialize symmetry information irreps required for input
!       processing. note, we assume
!       the # double group irreps here to be 128. FIXME: make it general
!       in
!       Dirac for higher linear symmetry double groups. "global Dirac
!       post-HF symmetry module".
        if(linear)then
          double_group = 11
          if(nfsym .eq. 1) double_group = 10
          call symmetry_setup_init(double_group,max_nkrmc_max_sym)
          allocate(mcscf_ciroots_symstr(4*max_nkrmc_max_sym))
        end if

        WORD1 = WORD
  100   CONTINUE
           READ (LUCMD, '(A7)') WORD
           CALL UPCASE(WORD)
           PROMPT = WORD(1:1)
           IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
              GO TO 100
           ELSE IF (PROMPT .EQ. '.') THEN
              ICHANG = ICHANG + 1
              DO 200 I = 1, NTABLE
                 IF (TABLE(I) .EQ. WORD) THEN
                    GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                    11,12,13,14,15,16,17,18,19,20,
     &                    21,22,23,24,25,26,27,28,29,30,
     &                    31,32,33,34,35,36), I
                 END IF
  200         CONTINUE
              IF (WORD .EQ. '.OPTION') THEN
                CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                GO TO 100
              END IF
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *           '" not recognized in *KRMCSCF.'
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
              CALL QUIT('Illegal keyword in *KRMCSCF.')
    1         CONTINUE
C&&& .THRESHOLD : convergence threshold
                  READ(LUCMD,*) DKRMC_THRGRD
               GO TO 100
    2          CONTINUE
C&&& .PRINT : General print level in KRMC-module
                  READ(LUCMD,*) IPRKRMC
               GO TO 100
    3          CONTINUE
C&&& .MAX MACRO ITERATIONS: Maximum number of macro iterations
                    READ(LUCMD,*) MKRMC_MXMACRO
               GO TO 100
    4          CONTINUE
C&&& .MAX MICRO ITERATIONS: Maximum number of micro iterations
                  READ(LUCMD,*) MKRMC_MXMICRO
               GO TO 100
    5          CONTINUE
C&&& .NR ALWAYS: Always use NR-iterations
                  MCFLAG(1) = .TRUE.
               GO TO 100
    6          CONTINUE
C&&& .NEO ALWAYS: Always use NEO-iterations
                  MCFLAG(2) = .TRUE.
               GO TO 100
    7          CONTINUE
C&&& .STATE : which state to converge to
                  READ(LUCMD,*) IKRMC_STATE
               GO TO 100
    8          CONTINUE
C&&& .DELETE SECONDARY ORBITALS:
                  DO I = 1, NFSYM
                     READ(LUCMD,'(A72)') KRMC_DELSTR(I)
                  END DO
               GO TO 100
    9          CONTINUE
C&&& .MAX BACKSTEPS : max. number of backsteps
                  READ(LUCMD,*) MKRMC_MAXBCK
               GO TO 100
   10          CONTINUE
C&&& .MAKE COFFEE:
                  WRITE(LUPRI,'(1X,A)')
     &              '*** ERROR in *KRMCSCF *** ',
     &              'Don''t know how to make coffee yet!'
                  CALL QUIT('*** ERROR in *KRMCSCF *** ' //
     &              'Don''t know how to make coffee yet!')
               GO TO 100
   11          CONTINUE
C&&& .SKIPEE: skip e-e rotations in the gradient and Hessian
                  KRMC_SKIPEE = .TRUE.
               GO TO 100
   12          CONTINUE
C&&& .SKIPEP: skip e-p rotations in the gradient and Hessian
                  KRMC_SKIPEP = .TRUE.
               GO TO 100
   13          CONTINUE
C&&& .CNVINT: Convergence thresholds for adding SL- and SS-integrals
                  READ(LUCMD,*) DKRMC_CNVINT(1),DKRMC_CNVINT(2)
               GO TO 100
   14          CONTINUE
C&&& .ITRINT: Number of iterations before adding SL- and SS-integrals
                  READ(LUCMD,*) IKRMC_ITRINT(1),IKRMC_ITRINT(2)
               GO TO 100
   15          CONTINUE
C&&& .INTFLG: Specify what two-integrals should be included in this run
                  READ(LUCMD,*) ILLINT,ISLINT,ISSINT
               GO TO 100
   16          CONTINUE
C&&& .NO CI : debug flag: no CI rotations.
                  KRMC_NOCI = .TRUE.
               GO TO 100
   17          CONTINUE
C&&& .THRPCI: threshold for printing CI vector
                  READ(LUCMD,*) KRMC_THRPCI
               GO TO 100
   18          CONTINUE
C&&& .FROZEN: freeze lowest inactive orbitals
                  READ(LUCMD,*) NKRMCFRO(1:NFSYM)
               GO TO 100
   19          CONTINUE
C&&& .MIN MK2: minimum value of 2 * MK
                  READ(LUCMD,*) NKRMC_MINMK2
               GO TO 100
   20          CONTINUE
C&&& .MAX MK2: maximum value of 2 * MK
                  READ(LUCMD,*) NKRMC_MAXMK2
               GO TO 100
   21          CONTINUE
C&&& .GASSH : GAS setup
                  IF (NKRMCGAS .GT. 0) THEN
                    CALL QUIT('*** ERROR in *KRMCSCF *** You may ' //
     &              'only specify one of .GASSH , .CASSCF, and .GAS SH')
                  END IF
                  READ(LUCMD,*) NKRMCGAS
                  IF (NKRMCGAS .LT. 1 .OR. NKRMCGAS .GT. MXGAS) THEN
                     WRITE(LUPRI,'(//A,I0)')
     &               '*** ERROR for *KRMCSCF/.GASSH  input *** ' //
     &               'invalid number of GAS spaces: ',NKRMCGAS
                     CALL QUIT('*** ERROR for *KRMCSCF input *** ' //
     &                     'invalid number of GAS spaces for .GASSH ')
                  END IF
                  DO I = 1, NKRMCGAS
                     READ(LUCMD,*) (NKRMCGSH(J,I),J=1,NFSYM)
                  END DO
               GO TO 100
   22          CONTINUE
C&&& .GASSPC: GAS space constraints
                  IF (NKRMCGAS .EQ. 0) THEN
                     WRITE(LUPRI,'(//A)') ' *** ERROR in *KRMCSCF ***'
     &                  //' .GASSH must be specified before .GASSP'
                     CALL QUIT('*** ERROR in *KRMCSCF ***')
                  END IF
                  DO I = 1, NKRMCGAS
                     READ(LUCMD,*) (NKRMCGSP(J,I),J=1,2)
                  END DO
                  IF (NKRMCAELEC.LT.0) NKRMCAELEC = NKRMCGSP(2,NKRMCGAS)
               GO TO 100
   23          CONTINUE
C&&& .MK2REFERENCE: 2 * M_K reference
                  READ(LUCMD,*) NKRMC_MK2REF
               GO TO 100
   24          CONTINUE
C&&& .MK2DELTA: 2 * DELTA M_K
                  READ(LUCMD,*) NKRMC_MK2DEL
               GO TO 100
 25            CONTINUE
C&&& .FREEZE INACTIVE ORBITALS:
                  WRITE(LUPRI,'(A)') ' *** ERROR in *KRMCSCF ***'//
     &            ' FREEZE option is not working yet! '//
     &            ' You have to be patient.'
C                 SK - Aug 2008
                  CALL QUIT(
     &         '*** ERROR in *KRMCSCF: FREEZE option not working yet.')
                  DO I = 1, NFSYM
                     READ(LUCMD,'(A72)') KRMC_FRZSTR(I)
                  END DO
               GO TO 100
 26            CONTINUE
C&&&  .CI PROGRAM: which CI program to use
                  READ(LUCMD,'(A)') KRMC_CIPROGRAM
               GO TO 100
 27            CONTINUE
C&&&  .INACTIVE ORBITALS
                  READ(LUCMD,*) (NKRMCISH(I),I=1,NFSYM)
               GO TO 100
 28            CONTINUE
C&&&  .CASSCF -- this defines a CAS space
                  IF (NKRMCGAS .GT. 0) THEN
                  CALL QUIT('*** ERROR in *KRMCSCF *** You may ' //
     &              'only specify one of .GASSH , .CASSCF, and .GAS SH')
                  END IF
                  NKRMCGAS = 1
                  READ(LUCMD,*) NKRMCAELEC
                  READ(LUCMD,*) (NKRMCGSH(I,NKRMCGAS),I=1,NFSYM)
                  NKRMCGSP(1,NKRMCGAS) = NKRMCAELEC
                  NKRMCGSP(2,NKRMCGAS) = NKRMCAELEC 
               GO TO 100
!
 29            CONTINUE
C&&&  .GAS SH:  the new input format in accordance with Dalton LUCITA (and soon Dirac LUCITA) 
!              and Diracs opens-shell AOC Hartree-Fock format
                  IF(NKRMCGAS .GT. 0)THEN
                     CALL QUIT('*** ERROR in *KRMCSCF *** You may ' //
     &              'only specify one of .GASSH , .CASSCF, and .GAS SH')
                  END IF

!              step 1: read the # of GA spaces
                  READ(LUCMD,*) NKRMCGAS
                  IF(NKRMCGAS .LT. 1 .OR. NKRMCGAS .GT. MXGAS) THEN
                    WRITE(LUPRI,'(//A,I0)')
     &              '*** ERROR for *KRMCSCF/.GAS SH input *** ' //
     &              'invalid number of GAS spaces: ',NKRMCGAS
                    CALL QUIT('*** ERROR for *KRMCSCF input *** ' //
     &                     'invalid number of GAS spaces for .GAS SH ')
                  END IF

!              process the min max occupation / per orbital occupation in each GAS shell
                  i = 0
                  do

                    i = i + 1
                    if(i > nkrmcgas) exit 

                    read(lucmd,'(a)') input_line
                    call upcase(input_line)
                    icomment = index(input_line,'!')
                    if (icomment > 0) input_line(icomment:) = ' '
                    icomment = index(input_line,'#')
                    if (icomment > 0) input_line(icomment:) = ' '
                    islash = index(input_line,'/')

!                check for correct input format: 
!                min e- max e- / #orbs_fermion_irrep1 #orbs_fermion_irrep2
                    if(islash <= 1)then

                      write(lupri,'(/a,i2,a,i2/a/2a)') 
     &             'ERROR for *KRMCSCF .GAS SHELL shell no.',
     &             i,'/',NKRMCGAS,
     &             '- the defining input line does not contain a "/":',
     &             '- the bad line : ',input_line

                      write(lupri,'(a)') 
     &             '- Or the specification of the number of GAS-'//
     &             'shells might be wrong! Check the input.'
                      call quit('Input error for .GAS SHELL under'//
     &                       ' *KRMCSCF')
                    end if
                
!                step 2: min max occupation
                    read(input_line(1:islash-1),*,iostat=ios) 
     &              nkrmcgsp(1,i), nkrmcgsp(2,i)

                    if(ios /= 0)then
                      write(lupri,'(/a,i2,a,i2/a/2a)') 
     &             'ERROR for *KRMCSCF .GAS SHELL shell no.',
     &             i,'/',NKRMCGAS, 
     &             '- the input line does not contain correct'//
     &             ' min max electrons','- the bad line : ',input_line
                      call quit('Input error for .GAS SHELL under '//
     &                       '*KRMCSCF')
                    end if
!                step 3: GAS shell orbital occupation
                    read(input_line(islash+1:),*,iostat=ios) test_string
                    islash = index(input_line,'/')
                    if(index(test_string,"ALL") .gt. 0)then
                      do j=1,nfsym
                        NKRMCGSH(j,i) = -9999
                      end do
                      KRMC_full_ci = .true.
                      exit ! jump out of loop
                    else
                      read(input_line(islash+1:),*,iostat=ios) 
     &                (NKRMCGSH(j,i), j=1,nfsym)
                    end if

                    if(ios /= 0)then
                      write(lupri,'(/a,i2,a,i2/a,i2,a/2a)') 
     &             'ERROR for *KRMCSCF .GAS SHELL shell no.',
     &             i,'/',NKRMCGAS,
     &             '- the input line does not contain',nfsym,
     &             ' occupations','- the bad line : ',input_line
                      call quit('Input error for .GAS SHELL under '//
     &                       '*KRMCSCF')
                    end if
                  END DO
!              step 4: set # of active electrons ==> max occupation in last GAS space
                  nkrmcaelec = nkrmcgsp(2,nkrmcgas)
               GO TO 100
 30            CONTINUE
C&&& .SYMMETRY: symmetry of the wave function
                  if(.not.linear) then
                    READ(LUCMD,*) IKRMC_SYMMETRY
                    IF (IKRMC_SYMMETRY .LE. 0) THEN
                       write(lupri,'(//A)') 'INPUT ERROR '//
     &                 'for *KRMCSF .SYMMETRY :',IKRMC_SYMMETRY
                       CALL QUIT('*KRMCSCF .SYMMETRY input error')
                    END IF
                  else
!                   read mij value and # roots in linear symmetry input format
                    read(lucmd,*) symmetry_str

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

!                   determine symmetry irrep of the wave function
!                   ---------------------------------------------

!                   a. check for inversion symmetry
                    ascii_lchar = iachar(last_char)
                    if(nfsym .gt. 1)then
                      if(ascii_lchar .eq. iachar('g'))then
                        gerade_mcscf_wf = .true.
                      else if(ascii_lchar .eq. iachar('u'))then
                        gerade_mcscf_wf = .false.
                      else 
                        write(lupri,'(/a,/a)') 
     & '  *** input error: Dinfh contains the inversion operation'//
     & ' as group element.',
     & '  please specify either g or u in addition to the mj value. ***'
                        call quit('  *** input error: Dinfh contains
     & the inversion operation as group element. ***')
                      end if
                    else
                      gerade_mcscf_wf  = .true.
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
                    if(.not.gerade_mcscf_wf) fermion_sym = 2
!                   print *,'fermion_sym, tmp_mij_val are ==> ',
!    &                       fermion_sym, tmp_mij_val
!                   call quit('==> bla bla ***')

!                   retrieve active dbg irrep (routine is inside the module symmetry_setup_krci)
                    call convert_mj_ferm_2_dbg_irrep_linsym_krci(
     &                   mj2rep,max_nkrmc_max_sym,tmp_mij_val,
     &                   fermion_sym,is_current_dbg_irrep)

!                   c. store symmetry of wave function == the active irrep
                    ikrmc_symmetry                             = 
     &              is_current_dbg_irrep
                    mcscf_ciroots_symstr(is_current_dbg_irrep) =
     &              tmp_symmetry_str(1:length_of_string)
                  end if
               GO TO 100
 31            CONTINUE
C&&& .NO1pDENS
                  KRMC_no1pdens = .true.
               GO TO 100
 32            CONTINUE
C&&& .MEMFAC: multiplier for subtracted scratch memory in
C             luciarel/ciinfo_r.F - Z_BLKFO_REL
                  READ(LUCMD,*) IKRMC_MEMFAC
                  IF( IKRMC_MEMFAC .LE. 0 ) THEN
                     WRITE(LUPRI,'(//3A,I6)')
     &                  '*** ERROR in *KRMCSCF *** ',
     &                  ' Illegal multiplier for ',
     &                  ' subtracting scratch memory: ',
     &                   IKRMC_MEMFAC
                     CALL QUIT('*** ERROR in *KRMCSCF *** ' //
     &                 'Illegal memory multiplier')
                  END IF
               GO TO 100
 33            CONTINUE
C&&&  .NOCKJZ  -- do not refine orbital rotation vector (kappa matrix) wrt
C                 the m_j values of each orbital
                 IF( .NOT. LINEAR )THEN
                   CALL QUIT('*** ERROR in *KRMCSCF *** .NOCKJZ ' //
     &             ' option is only allowed for linear symmetry')
                 END IF
                 KRMC_CHCKJZ = .FALSE.
               GO TO 100
 34            CONTINUE
C&&&  .SVRONO  -- save reordered NOs rather than the original NOs
               krmc_save_reordered_nos = .true.
               GO TO 100
! hjaaj Jan 2017: what was .SVRONX for? It was not activated in input TABLE
!C&&&  .SVRONX  -- save reordered NOs rather than the original NOs + active orbital offset
!               krmc_save_reordered_nos = .true.
!               READ(LUCMD,*) (IKRMC_SVRONO(I),I=1,NFSYM)
 35            CONTINUE
C&&&  .WITHEP  -- include e-p rotations in MCSCF optimization
                  KRMC_SKIPEP = .FALSE.
               GO TO 100
 36            CONTINUE
C&&&  .MVOFAC  -- factor on FV for modified virtual and inactive orbitals
                  READ(LUCMD,*) DKRMC_MVOFAC
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in *KRMCSCF.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in *KRMCSCF.')
            END IF
      END IF
  300 CONTINUE
C
      IF (.NOT.DOKRMC) GO TO 9999
C
!     warnings for options not implemented (yet)
!     ------------------------------------------
      if(lvnew)then
        write(lupri,'(/,a,/,a/)')
     & ' *** WARNING *** - (SS|SS)-integral approximation using LVNEW'//
     & ' is not implemented for KR-CI/KR-MCSCF/KR-CC.',
     & ' *** WARNING *** - Please use the LVCORR approximation.'
      end if
C
      IF(DOLVC) ISSINT = 0
      IF (LEVYLE) THEN
        ISLINT = 0
        ISSINT = 0
      ENDIF
      IF (ZORA.AND..NOT.ZORA4) THEN
        ISSINT = 0
      ENDIF
      IKRMC_INTDEF = ILLINT + 2*ISLINT + 4*ISSINT
      IF(DKRMC_CNVINT(1).LT.DUMMY) IKRMC_ITRINT(1) = 1
      IF(DKRMC_CNVINT(2).LT.DUMMY) IKRMC_ITRINT(2) = 1
C
C     set the # of active orbitals and Check that it is less than MAXASH
C     ------------------------------------------------------------------
      if(.not.KRMC_full_ci)then
        IF (NKRMCGAS .GT. 0) THEN
           NKRMCASHT = 0
           DO J = 1, NFSYM
              NKRMCASH(J) = 0
              DO I = 1, NKRMCGAS
                 NKRMCASH(J) = NKRMCASH(J) + NKRMCGSH(J,I)
              END DO
              NKRMCASHT = NKRMCASHT + NKRMCASH(J)
           END DO
        END IF
        NTEST = 0
        DO J = 1, NFSYM
           NTEST = NTEST + NKRMCASH(J)
        END DO
        IF (NTEST .GT. MAXASH) THEN
          WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &         ' *** ERROR in *KRMCSCF *** Too many active orbitals: ',
     &         NTEST, ' -- max is ',MAXASH
          CALL QUIT('*** ERROR in *KRMCSCF: Too many active orbitals!')
        END IF
      else
       ! skip the test - we set this later in TKRMCORB
      end if

C
C     Use boson symmetry for CI expansion for the spin-free and
C     the Levy-Leblond Hamiltonian, if KRMC_UCIBOS not set in input.
C
      IF(SPINFR)THEN
         KRMC_UCIBOS = .TRUE.
      END IF
C
C     spinfree MCSCF
      if((spinfr.or.levyle).and.(KRMC_CIPROGRAM.ne.'GASCIP'))then
        write(lupri,'(/a)') ' *** error in krmcinp: spinfree or'//
     *  ' Levy-Leblond Hamiltonian needs GASCIP as CI program. ***'
        call quit('*** error in krmcinp: spinfree or Levy-Leblond 
     &            Hamiltonian needs GASCIP as CI program.***')
      endif
c     krmc_ucibos = .false.
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &   ' *KRMCSCF: Set-up for KR-MCSCF calculation:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
      WRITE(LUPRI,'(A,I6)') ' * General print level   : ',IPRKRMC
      IF (NKRMCFRO(1) .ne. 0 .OR. NKRMCFRO(NFSYM) .ne. 0) THEN
         WRITE(LUPRI,'(A,T31,2I6)') ' * Frozen inactive orbitals',
     &      (NKRMCFRO(I),I=1,NFSYM)
         WRITE(LUPRI,'(A,T31,2I6)') ' * Relaxed inactive orbitals',
     &      (NKRMCISH(I)-NKRMCFRO(I),I=1,NFSYM)
      ELSE
         WRITE(LUPRI,'(A,T31,2I6)')
     &     ' * Inactive orbitals',(NKRMCISH(I),I=1,NFSYM)
      END IF
      if(.not.KRMC_full_ci)then
        WRITE(LUPRI,'(A,T31,2I6)')
     &     ' * Active orbitals  ',(NKRMCASH(I),I=1,NFSYM)
      else
        WRITE(LUPRI,'(A,T31,A)')
     &     ' * Active orbitals  ',' all remaining (code -9999)'
      end if
      WRITE(LUPRI,'(A,T31,I6)')
     &     ' * Active electrons ',NKRMCAELEC
      IF (NKRMCGAS .GT. 0) THEN
         WRITE(LUPRI,'(A,I3,A)')
     &     ' * GAS space setup for ',NKRMCGAS,' GAS space(s) : '
         N_GAS_ERR = 0
         DO I = 1, NKRMCGAS
            WRITE(LUPRI,'(A,I3,T31,2I6)')
     &     '   - GAS space',I,(NKRMCGSH(J,I),J=1,NFSYM)
            WRITE(LUPRI,'(A,I3,A1,I3)')
     &     '     min/max active electrons after this space :',
     &           NKRMCGSP(1,I),'/',NKRMCGSP(2,I)
C
C           Consistency checks:
C
            IF (NKRMCGSP(1,I) .GT. NKRMCGSP(2,I))
     &         N_GAS_ERR = N_GAS_ERR + 1
            IF (I .GT. 1) THEN
               IF (NKRMCGSP(1,I-1) .GT. NKRMCGSP(1,I))
     &            N_GAS_ERR = N_GAS_ERR + 1
               IF (NKRMCGSP(2,I-1) .GT. NKRMCGSP(2,I))
     &            N_GAS_ERR = N_GAS_ERR + 1
            END IF
         END DO
         IF (NKRMCGSP(2,NKRMCGAS) .NE. NKRMCAELEC)
     &      N_GAS_ERR = N_GAS_ERR + 1
         IF (N_GAS_ERR .GT. 0) THEN
            WRITE(LUPRI,'(//A)')
     &         'INPUT ERROR: The GAS specifications are not consistent'
            CALL QUIT('*** ERROR in *KRMCCALC *** ' //
     &        'Inconsistency in GAS specification')
         END IF

      END IF
      WRITE(LUPRI,'(A,T31,I4)')
     &     ' * State number',IKRMC_STATE
      WRITE(LUPRI,'(A,T31,A)')
     &     ' * CI program used :',KRMC_CIPROGRAM
C
C     Symmetry of the wave function:
C     ------------------------------
C
C
C     MK constraints in the CI expansion.
C
      if(.not.KRMC_full_ci)then
        IF(NKRMC_MINMK2 .ge. 0 .AND. NKRMC_MAXMK2.ge. 0)THEN
C
C          MINMK2/MAXMK2 explicitly set
C
        ELSE
           IF(NKRMC_MK2REF .ge. 0 .AND. NKRMC_MK2DEL .ge. 0)THEN
C
C             MK2REF/MK2DEL set, transform to MINMK2/MAXMK2
C
              NKRMC_MINMK2 = NKRMC_MK2REF - NKRMC_MK2DEL
              NKRMC_MAXMK2 = NKRMC_MK2REF + NKRMC_MK2DEL
C
           ELSE
C
C             Assign default values for MINMK2/MAXMK2
C
              NKRMC_MAXMK2 = 2*MIN(NKRMCAELEC,NKRMCASHT) - NKRMCAELEC
              NKRMC_MINMK2 = NKRMCAELEC-2*MIN(NKRMCAELEC,NKRMCASHT)
              NKRMC_MK2REF = 0
              NKRMC_MK2DEL = NKRMC_MAXMK2
           END IF
        END IF
      else
        NKRMC_MINMK2 = -9999
        NKRMC_MAXMK2 = -9999
        NKRMC_MK2REF =     0
        NKRMC_MK2DEL = -9999
      end if
C
C     Boson symmetry in the CI expansion
C
      IF ((SPINFR).AND.KRMC_UCIBOS) THEN
         NREP = NBSYM
         WRITE(LUPRI,'(A/2A/A,I4,A,I3)')
     &     ' * Boson spatial spinor symmetry used in CI expansion.',
     &     ' * Spatial symmetry (boson) of wave function : ',
     &        REP ( IKRMC_SYMMETRY - 1 ),
     &     ' * Allowed interval of 2 * MS :',
     &        NKRMC_MINMK2, ' to ', NKRMC_MAXMK2
      ELSE IF (.NOT.(SPINFR).AND.KRMC_UCIBOS) THEN
         NREP = NBSYM
         WRITE(LUPRI,'(A/2A/A,I4,A,I3)')
     &     ' * Approximate boson spatial spinor'//
     &        ' symmetry (as non-rel) used in CI expansion.',
     &     ' * Approximate spatial symmetry (boson) of wave function : '
     &        ,REP ( IKRMC_SYMMETRY - 1 ),
     &     ' * Allowed interval of approximate 2 * MS :',
     &        NKRMC_MINMK2, ' to ', NKRMC_MAXMK2
      ELSE
         if(linear)then
!          store on common block in krciprop.h
           xrepeig(ikrmc_symmetry) =
     &     mcscf_ciroots_symstr(ikrmc_symmetry)
           if(mod(nkrmcaelec, 2) .eq. 0)then
             ferm_lab = '  '
           else
             ferm_lab = '/2'
           end if
           WRITE(LUPRI,'(a,i3,a,a4,a2)')
     &        '    ** ', 1,' eigenstate for MJ-value (doubled): ',
     &          adjustr(xrepeig(ikrmc_symmetry)),ferm_lab
         else
           allocate(multb_tmp(64,64))
           IF (SPINFR) THEN
              CALL GMULTSF(NREP,REPNA,MULTB_TMP)
           ELSE
              CALL GMULTA(NZ,NFSYM,NREP,REPNA,MULTB_TMP)
           ENDIF
           deallocate(multb_tmp)
           IF ( MOD ( NKRMCAELEC, 2 ) .EQ. 0 ) THEN
              IRRP = NREP + IKRMC_SYMMETRY
              WRITE(LUPRI,'(2A)')
     &          ' * Symmetry (boson) of wave function: ',
     &          REPNA(IRRP)
C               store on common block in krciprop.h
                XREPEIG(IKRMC_SYMMETRY) = REPNA(IRRP)
           ELSE
              IRRP = IKRMC_SYMMETRY
              WRITE(LUPRI,'(2A)')
     &          ' * Symmetry (fermion) of wave function: ',
     &          REPNA(IRRP)
C               store on common block in krciprop.h
                XREPEIG(IKRMC_SYMMETRY) = REPNA(IRRP)
           END IF
         end if
         if(.not.KRMC_full_ci)then
           WRITE(LUPRI,'(A,I4,A,I3)')
     &       ' * Allowed interval of 2 * MK :',
     &       NKRMC_MINMK2, ' to ', NKRMC_MAXMK2
         else
           WRITE(LUPRI,'(A)')
     &       ' * Allowed interval of 2 * MK : will be set later!'
         end if
      END IF
C
      IF (.not.linear.and.IKRMC_SYMMETRY .GT. NREP) THEN
         WRITE(LUPRI,'(/A//A,2I3)')
     &      ' INPUT ERROR, requested .SYMMETRY is not valid!',
     &      ' symmetry requested, max value of symmetry = ', 
     &        IKRMC_SYMMETRY, NREP
         IF (KRMC_UCIBOS) THEN
            WRITE(LUPRI,'(I0, A,28(2X,A))')
     &      NREP,' spin-free boson symmetries known:',REP(0:NREP-1)
         ELSE
            WRITE(LUPRI,'(I0, A,28(2X,A))')
     &      NREP,' fermion symmetries known:',REPNA(1:NREP)
            WRITE(LUPRI,'(I0, A,28(2X,A))')
     &      NREP,' boson   symmetries known:',REPNA(NREP+1:2*NREP)
         END IF
         CALL QUIT('INPUT ERROR for *KRMCSCF  .SYMMETRY')
      END IF
C
      IF (KRMC_CIPROGRAM.EQ.'GASCIP  ') THEN
         if(linear)then
           call quit('*** error in KRMCINP: MCSCF calculations in'//
     &               ' linear symmetry are permitted only with the CI'//
     &               ' program LUCIAREL.***')
         end if
      END IF
C
      IF (DKRMC_MVOFAC .NE. 1.0D0) THEN
         WRITE(LUPRI,'(A,1P,D10.2)')
     &   ' * MVOs based on FC + a*FV,  a =    ',DKRMC_MVOFAC
      END IF
C
      WRITE(LUPRI,'(A,I4)')
     &   ' * Maximum number of macro iterations  : ',MKRMC_MXMACRO
      WRITE(LUPRI,'(A,I4)')
     &   ' * Maximum number of micro iterations  : ',MKRMC_MXMICRO
      WRITE(LUPRI,'(A,1P,D10.2)')
     &   ' * Convergence threshold               : ',DKRMC_THRGRD
      IF (MCFLAG(1) .AND. MCFLAG(2)) MCFLAG(1) = .FALSE.
      IF (MCFLAG(1)) WRITE(LUPRI,'(A)')
     &   ' * Only Newton-Raphson iterations.'
      IF (MCFLAG(2)) WRITE(LUPRI,'(A)')
     &   ' * Only NEO iterations.'
      WRITE(LUPRI,'(A)')
     +    ' * Contributions from 2-electron integrals to Fock matrix:'
      IF(LBIT(IKRMC_INTDEF,1)) THEN
        WRITE(LUPRI,'(3X,A)') 'LL-integrals.'
      ENDIF
      IF(LBIT(IKRMC_INTDEF,2)) THEN
        IF(DKRMC_CNVINT(1).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,E8.1)')
     &    'SL-integrals  below KR-MCSCF convergence ',DKRMC_CNVINT(1)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     &    'SL-integrals from iteration ',IKRMC_ITRINT(1)
        ENDIF
      ENDIF
      IF(LBIT(IKRMC_INTDEF,3)) THEN
        IF(DKRMC_CNVINT(2).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,E8.1)')
     &    'SS-integrals below KR-MCSCF convergence ',DKRMC_CNVINT(2)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     &    'SS-integrals from iteration ',IKRMC_ITRINT(2)
        ENDIF
      ENDIF
!     if number of positronic shells is zero, skip e-p rotations. hjj+sk - aug 2010
      IF(x2c.or.bss.or.levyle.or.freepj.or.vextpj) THEN
         IF(.NOT. KRMC_SKIPEP) THEN
            WRITE(LUPRI,*) ' * INFO: .WITHEP ignored for '//
     &      ' X2C, BSS, LEVYLE, FREEPJ, and VEXTPJ'
            KRMC_SKIPEP = .TRUE.
         END IF
      ELSE IF(KRMC_SKIPEP) THEN
         WRITE(LUPRI,'(A)')
     &     ' * No-pair KR-MCSCF: No rotations between electrons and'//
     &     ' positrons in 2nd order optimization.'
      END IF
      IF(KRMC_SKIPEE) WRITE(LUPRI,'(A)')
     &     ' * NB!!! No e-e rotations in 2nd order optimization'
      IF(KRMC_CHCKJZ) WRITE(LUPRI,'(A)')
     &     ' * Checking kappa matrix for orbital rotations'//
     &     ' with respect to the mj-values of the orbitals.'

C
C     Frozen orbitals:
C     ----------------
C
      NVECT = 0
      DO I = 1, NFSYM
         NVEC = 0
         CALL NUMLST(KRMC_FRZSTR(I),IDUMMY,NFBAS(I,0),
     &        -NFBAS(I,2),NFBAS(I,1),I,NVEC)
         NVECT = NVECT + NVEC
      END DO
      IF (NVECT .GT. 0) THEN
         WRITE(LUPRI,'(A)') ' * Frozen orbitals: '
         WRITE(LUPRI,'(A,I1,2A)')
     &        ('   - symmetry ',I,' : ',KRMC_FRZSTR(I),I=1,NFSYM)
      END IF
C
C
C     Deleted orbitals:
C     -----------------
C
      NVECT = 0
      DO I = 1, NFSYM
         NVEC = 0
         CALL NUMLST(KRMC_DELSTR(I),IDUMMY,NFBAS(I,0),
     &        -NFBAS(I,2),NFBAS(I,1),I,NVEC)
         NVECT = NVECT + NVEC
      END DO
      IF (NVECT .GT. 0) THEN
         WRITE(LUPRI,'(A)') ' * Deleted orbitals: '
         WRITE(LUPRI,'(A,I1,2A)')
     &        ('   - symmetry ',I,' : ',KRMC_DELSTR(I),I=1,NFSYM)
      END IF
      IF( KRMC_no1pdens ) WRITE (LUPRI, '(A)' )
     &     ' * DO NOT save the 1-particle density matrix obtained'//
     &     ' from the converged MCSCF wave function'
      IF( IKRMC_MEMFAC .ne. 9 ) WRITE (LUPRI, '(A,I3)' )
     &     ' * Multiplier for memory guess: ',IKRMC_MEMFAC
C
 9999 if(linear) then
        if(allocated(mcscf_ciroots_symstr))
     &    deallocate(mcscf_ciroots_symstr)
      end if
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rstpinp */
      SUBROUTINE RSTPINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for STEP CONTROL-module
C
C     Written by J. Thyssen - Nov 11 1998
C     Last revision:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dummy.h"
      PARAMETER (NTABLE = 16)
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
#include "dcbgen.h"
#include "dcbopt.h"
C
      SAVE SET
      DATA TABLE /'.DAMPIN','.MAX DA','.MIN DA','.NO EXT','.GOOD R',
     &            '.MIN RA','.REJECT','.TRUST ','.TOLERA','.INCREM',
     &            '.MAX ST','.DECREM','.TIGHT ','.THQMIN','.THQLIN',
     &            '.THQKVA'/
      DATA SET/.FALSE./
C
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
C
C     Process input for STEP CONTROL
C     ==============================
C
      NEWDEF = (WORD .EQ. '*STEP C')
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
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in DHFINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in STEP CONTROL.')
    1          CONTINUE
C *** Option 1 >DAMPIN <  Initial value of damping (BETA)
                  READ(LUCMD,*) BETA
               GO TO 100
    2          CONTINUE
C *** Option 2 >MAX DA<  Maximum damping value
                  READ(LUCMD,*) BETMAX
               GO TO 100
    3          CONTINUE
C *** Option 3 >MIN DA<  Minimum damping value
                  READ(LUCMD,*) BETMIN
               GO TO 100
    4          CONTINUE
C *** Option 4 >NO TER<  Skip extra termination tests
                  FLAG(20) = .TRUE.
               GO TO 100
    5          CONTINUE
C *** Option 5 >GOOD R<  "Good" trust radius
                  READ(LUCMD,*) RATGOD
               GO TO 100
    6          CONTINUE
C *** Option 6 >MIN RA<  "Bad" trust radius
                  READ(LUCMD,*) RATBAD
               GO TO 100
    7          CONTINUE
C *** Option 7 >REJ RA<  Ratio predicted/actual step that rejects step
                  READ(LUCMD,*) RATREJ
               GO TO 100
    8          CONTINUE
C *** Option 8 >TRUST <  Initial trust radius
                  READ(LUCMD,*) RTRUST
               GO TO 100
    9          CONTINUE
C *** Option 9 >TOLERA<  Tolerable trust radius
                  READ(LUCMD,*) RTTOL
               GO TO 100
   10          CONTINUE
C *** Option 10 >INCREM<  Increment factor on trust radius
                  READ(LUCMD,*) STPINC
               GO TO 100
   11          CONTINUE
C *** Option 11 >MAX ST<  Maximum acceptable step length
                  READ(LUCMD,*) STPMAX
               GO TO 100
   12          CONTINUE
C *** Option 12 >DECREM<  Decrement factor on trust radius
                  READ(LUCMD,*) STPRED
               GO TO 100
   13          CONTINUE
C *** Option 13 >TIGHT <  Tight step contol for ground state also
                  FLAG(19) = .TRUE.
               GO TO 100
   14          CONTINUE
C *** Option 14 >THQMIN<
                  READ(LUCMD,*) THQMIN
               GO TO 100
   15          CONTINUE
C *** Option 15 >THQLIN<
                  READ(LUCMD,*) THQLIN
               GO TO 100
   16          CONTINUE
C *** Option 16 >THQKVA<
                  READ(LUCMD,*) THQKVA
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in STEP CONTROL.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in STEP CONTROL.')
            END IF
      END IF
  300 CONTINUE
C
C
C     Print section
C     =============
C
      IF (.NOT. NEWDEF) GOTO 400
C
C     print stuff
C
C
  400 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck roptinp */
      SUBROUTINE ROPTINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     MCSCF optimization control input section.
C
C     Written by J. Thyssen - Dec 29 2000
C                    Added CI keyword NOOCCN. Stefan Knecht October 2007
C                    New start guess code 1531. Read MP2 natural 
C                    orbitals from KRMCOLD as default. Stefan Knecht 07.07.2008
C                    Added CI/KRMC keyword GENFOC. Read in orbitals
C                    and vector to generate Fock-type orbitals for 
C                    positronic/inactive/secondary space. Stefan Knecht 17.07.2008
C                    
C
C     Major revision by S. Knecht: 08.08.2008
C
C                    All CI-related keywords have been moved to a new 
C                    input section *KRCICALC. 
C                    Purpose: disentanglement of MCSCF/CI optimization
C                    input.
C
C***********************************************************************
      use os_utils
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 25)
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
#include "dcbgen.h"
#include "dcbopt.h"
C
      SAVE SET
      DATA TABLE /'.MAX CI','.IOPTST','.TRACI ','.NO TRA','.NATONL',
     &            '.FOCKON','.NO TRS','.COMPRE','.COMP F','.RSTRCI',
     &            '.NO PFQ','.NO FQX','.ANALYZ','-------','.MXCIVE',
     &            '.KTRLVL','.PARINT','-------','.NOOCCN','.RDFOCK',
     &            '.WRTFCK','-------','-------','-------','.GENFOC'/
      DATA SET/.FALSE./
C
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
C     Initialization
C     --------------
C
      MAXCIT      =    4
      MXCIV       =    0
      IOPTST      = 1563
      KTRLVL      =    3
      IKRMCCNO    =  - 3
      IANACI      =    0
      IRESTRK     =    0
      IDOPARIO    = get_environment_integer('GLBSCR',IDOPARIO)
      OPT_NOPFQ   = .FALSE.
      OPT_NOFQX   = .FALSE.
      NATONL      = .FALSE.
      FOCKON      = .FALSE.
      COMPROT     = .FALSE.
      GENFOCK     = .FALSE.
      NATOLCR     = .FALSE.
      TRA_NATO    = .FALSE.
      CRDFO_MAT   = .FALSE.
      CWRTFO_MAT  = .FALSE.
      COMPFAC     = 0.01D00

      if(IDOPARIO .ne. 0 .or. IDOPARIO.ne. 1)then
!       fall back to default: 0
        IDOPARIO = 0
      end if
C
C
C     Process input for *OPTIMIZE
C     ===========================
C
      NEWDEF = (WORD .EQ. '*OPTIMI')
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
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
     &                     , I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized on deck *OPTIMIZE.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword on deck *OPTIMIZE.')
    1          CONTINUE
C&&& .MAX CI ITER...: maximum number of initial CI iterations.
                  READ(LUCMD,*) MAXCIT
               GO TO 100
    2          CONTINUE
C&&& .IOPTST: start guess flag
C             (FIXME: make this human readable)
                   READ(LUCMD,*) IOPTST
               GO TO 100
    3          CONTINUE
C&&& .TRACI: call rtraci for transform of CI vectors
                  IKRMCCNO = 1
               GO TO 100
    4          CONTINUE
C&&& .NO TRACI: don't call rtraci for transform of CI vectors
                  IKRMCCNO = -1
               GO TO 100
    5          CONTINUE
C&&& .NATONL: only natural orbitals (no canonical orbitals)
                  IKRMCCNO = 1
                  NATONL = .TRUE.
                  FOCKON = .FALSE.
               GO TO 100
    6          CONTINUE
C&&& .FOCKON: only Fock-type orbitals (no natural orbitals)
                  IKRMCCNO = 1
                  FOCKON = .TRUE.
                  NATONL = .FALSE.
               GO TO 100
    7          CONTINUE
C&&& .NO TRANS: no transformation of orbitals
                  IKRMCCNO = -2
               GO TO 100
    8          CONTINUE
C&&& .COMPROT: compress rotations
                  COMPROT = .TRUE.
               GO TO 100
    9          CONTINUE
C&&& .COMP FACTOR: threshold factor for removal of orbitals
                  READ(LUCMD,*) COMPFAC
               GO TO 100
   10          CONTINUE
C&&& .RSTRCI: Restart CI from vector on file (1) or not (0)
                  READ(LUCMD,*) IRESTRK
               GO TO 100
   11          CONTINUE
C&&& .NO PFQ: no positronic FQ
                  OPT_NOPFQ = .TRUE.
               GO TO 100
   12          CONTINUE
C&&& .NO FQXT: no FQX and FQT (approximate with DV * FV)
                  OPT_NOFQX = .TRUE.
               GO TO 100
   13          CONTINUE
C&&& .ANALYZ: Analyze LUCIAREL CI vector(s)
                  IANACI = 1
               GO TO 100
   14          CONTINUE
               GO TO 100
   15          CONTINUE
C&&& .MXCIVE: Number of CI vectors allowed for subspace
                  READ(LUCMD,*) MXCIV
               GO TO 100
   16          CONTINUE
C&&& .KTRLVL: default integral transformation level.
               READ(LUCMD,*) KTRLVL
               IF (KTRLVL .LE. -2 .OR. KTRLVL .GT. 5) THEN
                  WRITE(LUPRI,'(//2A,I10)')
     &                '*** ERROR in ROPTINP *** ',
     &                'Illegal transformation level: ',KTRLVL
                  CALL QUIT('*** ERROR in ROPTINP *** ' //
     &                 'Illegal transformation level')
               END IF
               GO TO 100
   17          CONTINUE
C&&& .PARINT: integral handling for slaves in LUCIAREL calculations
                  IDOPARIO_SAVE = IDOPARIO
                  READ(LUCMD,*) IDOPARIO
                  IDOPARIO = ABS(IDOPARIO)
                  IF( IDOPARIO .ne. IDOPARIO_SAVE )THEN
C                   we do not want to broadcast integrals if not
C                   necessary! - SK- 20-03-2008
                    IDOPARIO = IDOPARIO_SAVE
                  END IF
#if !defined (VAR_MPI)
                  IF( IDOPARIO .ne. 0 ) THEN
                    WRITE(LUPRI,'(//3A,I6)')
     &                    '*** ERROR in ROPTINP *** ',
     &                    ' No slaves in sequential run, no integral
     &                      handling needed: ',
     &                     IDOPARIO
                    IDOPARIO = 0
                  END IF
#else
                  IF( IDOPARIO .gt. 1 ) THEN
                   WRITE(LUPRI,'(//3A,I6)')
     &                    '*** ERROR in ROPTINP *** ',
     &                    'Illegal integral handling option for
     &                     slaves in parallel run: ',
     &                     IDOPARIO
                   IDOPARIO = 0
                  END IF
#endif
               GO TO 100
   18          CONTINUE
               GO TO 100
   19          CONTINUE
C&&& .NOOCCN: calculate natural orbitals occupation numbers (with LUCIAREL)
               NATOLCR      = .TRUE.
               GO TO 100
   20          CONTINUE
C&&& .RDFOCK: read fock matrix from file useful for LUCIAREL restarts
               CRDFO_MAT = .TRUE.
               GO TO 100
   21          CONTINUE
C&&& .WRTFCK: write fock matrix in CI start guess module
               CWRTFO_MAT = .TRUE.
               GO TO 100
   22          CONTINUE
               GO TO 100
   23          CONTINUE
               GO TO 100
   24          CONTINUE
               GO TO 100
   25          CONTINUE
C&&& .GENFOC: generate Fock-type orbitals as start orbitals for CI 
C             calculations (starting, e.g., from a previous (converged) 
C             MCSCF wavefunction.
               GENFOCK = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/3A/)') ' Prompt "',WORD,
     *            '" not recognized in *OPTIMIZE.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in *OPTIMIZE.')
            END IF
      END IF
  300 CONTINUE
C
      IF (.NOT. NEWDEF) GOTO 400
C
C     If OPT_NOFQX then we only need 1st order integrals transformation:
C
      IF (OPT_NOFQX) THEN
         KTRLVL = 1
      END IF
C
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &   ' ROPTINP: Set-up for MCSCF wave function optimization'
      CALL PRSYMB(LUPRI,'=',75,0)
C
      IF (MAXCIT .GE. 0) THEN
         WRITE(LUPRI,'(A,I3)')
     &     ' * Maximum number of initial CI iterations: ',MAXCIT
      ELSE
         WRITE(LUPRI,'(A)')
     &     ' * CI start guess is based on full CI '
      END IF
C
      IF (KTRLVL.EQ.5) THEN
         WRITE(LUPRI,'(A)')
     &     ' * Running only integral transformation in KRMC '
      END IF
C
      IF (IRESTRK.EQ.1) THEN
         WRITE(LUPRI,'(A)')
     &     ' * Restarting CI from C vector on LUCIAREL file '
      END IF
C
      IF (IANACI.EQ.1) THEN
         WRITE(LUPRI,'(A)')
     &     ' * Analyzing CI vectors in terms of coefficients '
      END IF
C
      IF (MXCIV.GT.0) WRITE(LUPRI,'(A,I3)')
     &     ' * Maximum subspace dimension set to ',MXCIV
C
      IF( IDOPARIO .eq. 1 )THEN
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
C
      WRITE(LUPRI,'(A,I4)')
     &     ' * Start guess code : ',IOPTST
C
      IF ( OPT_NOPFQ ) WRITE (LUPRI, '(A)' )
     &     ' * Positronic part of positronic FQX ignored.'
C
      IF ( OPT_NOFQX ) WRITE (LUPRI, '(A)' )
     &     ' * FQX and FQT approximate with DV * FV.'
C
C     Transformation of CI vectors and orbitals.
C     ------------------------------------------
C
      IF (IKRMCCNO .EQ. -1) THEN
         FOCKON = .TRUE.
      END IF
C
      IF ( IKRMCCNO .GE. -1 ) THEN
         IF ( IKRMCCNO .GE. 1 ) WRITE(LUPRI,'(A)')
     &        ' * Call RTRACI for transformation of CI vectors.'
         IF ( IKRMCCNO .EQ. 0 ) WRITE(LUPRI,'(A)')
     &        ' * Call RTRACI for transformation of CI vectors ' //
     &        '(if necesary).'
         IF ( IKRMCCNO .EQ.-1 ) WRITE(LUPRI,'(A)')
     &        ' * Do not call RTRACI for transformation of CI vectors.'
C
         IF ( FOCKON ) THEN
            WRITE(LUPRI,'(A)')
     &           '   - transform to Fock-type orbitals.'
         ELSE IF (NATONL) THEN
            WRITE(LUPRI,'(A)')
     &           '   - transform to natural orbitals.'
         ELSE
            WRITE(LUPRI,'(A)')
     &           '   - transform to Fock-type and natural orbitals.'
         END IF
      ELSE
         WRITE(LUPRI,'(A)') ' * No transformation of orbitals'//
     &                      ' to Fock-type or natural orbitals.'
      END IF
C     ... only natural orbitals occupation numbers
      IF( NATOLCR )THEN
        WRITE(LUPRI,'(A)') ' * Calculation of nat. orb. occ. numbers'
      END IF
C
      IF (COMPROT) WRITE(LUPRI,'(A/A,1P,D12.2,A)')
     &     ' * Compress orbital rotations during optimization',
     &     '   (threshold is ',COMPFAC,' * 1-norm(g)).'
C
      WRITE(LUPRI,'(A,I3)')
     &     ' * Default integral transformation level  :',KTRLVL
      
      IF( GENFOCK )THEN
        WRITE(LUPRI,'(A)') 
     &     ' * Generate Fock-type positronic/inactive/virtual '//
     &     '   starting orbitals'
      END IF
C
  400 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck optini */
      SUBROUTINE OPTINI()
C***********************************************************************
C
C     Initialize data in dcbopt.h
C     Written by J. Thyssen - Nov 11 1998
C     Last revision:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
C
C     Initialize parameter info  pointers to negative value (for "not defined")
C
      KZCONF = -1
      LZCONF = -1
      LZXOPE = -1
      LZXOPP = -1
C
C     Default values for STPINP
C
      BETA   = 1.0D0
      BETMIN = 0.2D0
      BETMAX = 1.D6
      RTTOL  = 1.02D0
      RATREJ = 0.25D0
      RATMIN = 0.40D0
      RATGOD = 0.80D0
      STPMAX = 0.70D0
      STPINC = 1.2D0
      STPRED = 0.67D0
      RTRUST = STPMAX
      THQMIN = 0.1D0
      THQLIN = 0.2D0
      THQKVA = -1.0D0

      RETURN
      END
