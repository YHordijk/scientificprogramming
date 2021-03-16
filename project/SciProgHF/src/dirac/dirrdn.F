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
C FILE    : dirac/dirrdn.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck paminp */
      SUBROUTINE PAMINP()
C***********************************************************************
C
C     Driver for input module of DIRAC
C     Written by T.Saue April 26 1996
C
C     Called from:  DIRCTL (1.time in DIRAC run)
C                   EXEDIR (during optimization)
C
C     Last revision: 260496 - tsaue
C               Febr.2006 - MI to modify for pure two-comp. mode
C
C***********************************************************************

      use dirac_cfg
      use memory_allocator
      use num_grid_cfg
      use interface_functional_read
      use xmlout
#ifdef HAS_PCMSOLVER      
      use pcmmod_cfg
      use pcm_write
#endif
#ifdef HAS_PELIB
      use pe_variables, only: peqm
      use read_xyzfile, only: pelib_sym_xyz
#endif
      use x2cmod_cfg
#ifdef MOD_AOOSOC
      use atomic_oo_order_so_correction_cfg
#endif

#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbimp.h"
#include "dcborb.h"
#include "cbirea.h"

      integer, parameter   :: lwork = 100000
      real(8), allocatable :: WORK(:)
      character(100)       :: line
      integer, parameter   :: file_unit = 137
      logical              :: inp_exists, mol_exists, xyz_exists
      logical              :: input_found

C
      CALL QENTER('PAMINP')
      call alloc(WORK,LWORK, id="PAMINP WORK array")
C
C     Initializations
C
      CALL PAMINI()

!     print contents of input file to output
      inquire(file = 'DIRAC.INP', exist = inp_exists)
      if (inp_exists) then
         open(file_unit,
     &        file   = 'DIRAC.INP',
     &        status = 'unknown',
     &        form   = 'formatted',
     &        access = 'sequential')
         rewind(file_unit)
      else
         call quit('*.inp file not found')
      end if
      write(lupri, '(//a)') 'Contents of the input file'
      write(lupri, '(a/)')  '--------------------------'
      do while (.true.)
         read(file_unit, '(a100)', end=1) line
         write(lupri, '(a)') line
      end do
 1    continue
      close(file_unit, status = 'keep')

!     print contents of molecule file to output
      inquire(file = 'MOLECULE.MOL', exist = mol_exists)
      if (mol_exists) then
         open(file_unit,
     &        file   = 'MOLECULE.MOL',
     &        status = 'unknown',
     &        form   = 'formatted',
     &        access = 'sequential')
         rewind(file_unit)
      else
         inquire(file = 'MOLECULE.XYZ', exist = xyz_exists)
         if (xyz_exists) then
            open(file_unit,
     &           file   = 'MOLECULE.XYZ',
     &           status = 'unknown',
     &           form   = 'formatted',
     &           access = 'sequential')
            rewind(file_unit)
         else
            call quit('*.mol/*.xyz file not found')
         end if
      end if
      write(lupri, '(//a)') 'Contents of the molecule file'
      write(lupri, '(a/)')  '-----------------------------'
      do while (.true.)
         read(file_unit, '(a100)', end=2) line
         write(lupri, '(a)') line
      end do
 2    continue
      close(file_unit, status = 'keep')

#ifdef DO_KEYWORD_SCAN

! hjaaj Mar 2012: discovered that the NUMBER_OF_KEYWORDS is not used
!    for anything currently, thus I disabled the code.

      OPEN(LUCMD,FILE = 'DIRAC.INP',IOSTAT=IOS)
      IF (IOS.NE.0) THEN
        WRITE(LUPRI,*) 'Error in opening DIRAC.INP!'
        CALL QUIT('Error in opening DIRAC.INP!')
      ENDIF
      CALL KEYWORD_SCAN (.TRUE.,.TRUE.,.TRUE.,.TRUE.,
     &                   .FALSE.,NUMBER_OF_KEYWORDS)
      CLOSE(LUCMD)
#endif
C
C     Specify job
C
      OPEN(LUCMD,FILE = 'DIRAC.INP',IOSTAT=IOS)
      IF (IOS.NE.0) THEN
        WRITE(LUPRI,*) 'Error in opening DIRAC.INP!'
        CALL QUIT('Error in opening DIRAC.INP!')
      ENDIF
      CALL JOBINP
      CLOSE(LUCMD)

C     Start the xml output and echo the input to the xml output
      if (doxml) then
         call xml_begin('program')
         call xml_begin('name')
         write(luxml,'(A)') 'dirac'
         call xml_end('name')
         call xml_begin('version')
         write(luxml,'(A)') '2019'
         call xml_end('version')
         call xml_end('program')
         call xml_begin('input')
         call xml_embed_formatted('DIRAC.INP','dirac-input')
         if (mol_exists) then
            call xml_embed_formatted('MOLECULE.MOL','dirac-molecule')
         else
            call xml_embed_formatted('MOLECULE.XYZ','dirac-xyz')
         end if
         call xml_end('input')
      endif

C     General input

      OPEN(LUCMD,FILE = 'DIRAC.INP')
      REWIND LUCMD
      CALL GENINP(WORK,LWORK)
      CLOSE(LUCMD)

C     Integral and basis input

      OPEN(LUCMD,FILE = 'DIRAC.INP')
      CALL HERINP(WORK,LWORK)
      CLOSE(LUCMD)

      DOJACO = DOJACO.AND.(NZ.EQ.1)

C     Set DCB common blocks - first time
C     ========================================
C     (will be done the second time in PAMSET)

      CALL SETDC1(IPREAD)

C     Make primitive labels and use label
C     information to list SO-basis
C     ===================================

      CALL GETLAB(IPREAD)

C     Define Hamiltonian

      OPEN(LUCMD,FILE='DIRAC.INP')
      CALL HAMINP(WORK,LWORK)
      CLOSE(LUCMD)

#ifdef HAS_PELIB
C     Initialize PE library
       IF (PEQM) THEN
         ! first check if symmetry is disabled (for xyz input)
         IF (PELIB_SYM_XYZ) CALL QUIT(".PElib/xyz input:" //
     &                                  " Activate .NOSYM!")
         CALL PELIB_IFC_INIT()
       END IF
#endif
!     this is the new input reader
!     contains:
!     *DFT
!     *VISUAL
!     **RELCC
!     *OPENRSP
!     *X2C
!     *PCM
!     *PCMSOL      
!     *AooSOC
!     *QCORR
!     ignores rest
      call read_menu_input('DIRAC.INP', lucmd, 'ALL', input_found)

      if (dirac_cfg_dft_calculation) then
        call consistency_after_dft_input()
        call report_after_dft_input()
!       fixme num grid report should be independent of DFT
        call report_num_grid()
      end if
!     consistency check for fragment X2C and the proper X2C module
      if(x2cmod_fragment_x2c.and..not.x2c)then
        write(lupri,'(/a,a)')
     &  '  *** PAMINP reports an error after readin *X2C input:',
     &  '  a fragment X2C calculation requires the X2C module'//
     &  ' via .X2C under **HAMILTONIAN.'
        call quit('*** fragment X2C calculation not possible. ***')
      end if
!     consistency check for atomic oo-order SOC contributions
#ifdef MOD_AOOSOC
      if(aoomod)then
        write(lupri,*) 'Information about the oo-order SOC module:'
        write(lupri,*) '* atomic oo-order SOC module active'
        x2c_2c_mmf_mos = .false.
        if(x2cmod_mmf)then
          write(lupri,*) '  - generate 2e-SRCs/SOCs to be added to'//
     &                   ' the X2C Hamiltonian (of a molecule)'
          x2c_add_amfi   = -1
          x2c_2c_mmf_mos = .true.
        else
          write(lupri,*) '  - add 2e-SOCs to the X2C Hamiltonian'
        end if
        if(.not.x2c)then
          call quit(' aooSOC does not make sense in NON-2c '//
     &      ' relativistic runs')
        end if
      end if
#endif
! RDR Polarizable Continuum Model calculation
#ifdef HAS_PCMSOLVER      
      if (dirac_cfg_pcm) then
         call report_after_pcm_input(lupri)
      end if
#endif
CMI   ... instead of calling OP1INI in HAMINP, call it after
      CALL OP1INI(WORK,LWORK)
C
C     Wave function input
C
C     Call PSIINP always to get wave function information for
C     other modules. This could be changed if wave function information
C     is saved on file.  TODO ?
C
      IF ( DOPSI .OR. .NOT. INPTES ) THEN
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL PSIINP(WORK,LWORK)
         CLOSE(LUCMD)
      END IF
C
C     Input of wave function analysis module
C HJAaJ May 2002:
C   more logical to move CALL PAMANA to PAMPSI, so it follows directly HF/DFT in output
C   even if e.g. MP2 or RELCCSD is called, because .ANALYZE is currently an analysis
C   of the HF/DFT wave function, not of the following output.
C   - in a way *ANALYZE ought to be changed to *SCF ANALYZE and
C     .ANALYZE moved to *SCF INPUT ??
C   - one could use the analysis for other wave functions with some changes ?
C
      IF(DOANA) THEN
        OPEN(LUCMD,FILE = 'DIRAC.INP')
        CALL ANWINP('**ANALY',WORK,LWORK)
        CLOSE(LUCMD)
      ENDIF
C
C     Input of property module
C
      IF(DOPRP) THEN
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL PRPINP('**PROPE',WORK,LWORK)
         CLOSE(LUCMD)
      ENDIF
C
C     Input of transformation module
C
      IF(DOTRA) THEN
        OPEN(LUCMD,FILE = 'DIRAC.INP')
        CALL TRAINP(WORK,LWORK)
        CLOSE(LUCMD)
      ENDIF

C
C     Consistency checks
C
      CALL CONINP

C
C     set number of blocks in case of import
C     complex spinors from ReSpect (TODO) or Turbomole
C
      NKRBLK = 1
      IF (DOIMPMOS) THEN
       IF ((FAKE2C .OR. NZC1.EQ.4) .AND.
     &     INDEX(QCPACK,"TURBOMOLE").GT.0) THEN
        NKRBLK = 2
       END IF
      END IF
C
      CALL TITLER('End of input processing','*',124)
C
      call dealloc(WORK)
      CALL QEXIT('PAMINP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck geninp */
      SUBROUTINE GENINP(WORK,LWORK)
C***********************************************************************
C
C <<< General Input for DIRAC >>>
C
C       This routine is based on the corresponding HERINP routine
C       in HERMIT
C
C       Written by Trond Saue, Tromsoe October 1993
C       Last revision: Nov 14 1994 - tsaue
C
C***********************************************************************
        use codata
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "consts.h"
      PARAMETER (NDIR = 8, NTABLE = 24)
C
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7, MXFORM*6, FMT*6
      DIMENSION WORK(LWORK)
      LOGICAL   LBIT, FIND_KEYWORD, NOGENINP, USER_STOL, SET_CVALUE
      REAL*8    CSCALE
C
#include "dcbgen.h"

! MDIRAC is needed for check of STOL(2)
#include "dcbham.h"

#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcb_rkbimp.h"
#include "dcbimp.h"
#include "dcborb.h"
#include "logging.h"
#include "infpar.h"
C
      DATA TABDIR /'*END OF','*PARALL','*XXXXXX','*XXXXXX',
     &             '*XXXXXX','*XXXXXX','*XXXXXX','*XXXXXX'/
      DATA TABLE  /'.DIRECT','.CVALUE','.NOSET ','.DOJACO',
     &             '.PRINT ','.PRJTHR','.PCMOUT','.LINDEP',
     &             '.SPHTRA','.ECP   ','.LOGMEM','.QJACO ',
     &             '.ACMOUT','.ACMOIN','.IMPMOS','.FAMILY',
     &             '.SKIP2E','.LOWJAC','.RKBIMP','.CODATA',
     &             '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
C
C     Initialise NOGENINP and USER_STOL
C     NOGENINP is set, if there is no **GENERAL
C
      NOGENINP   = .FALSE.
      USER_STOL  = .FALSE.
      TWOCOMP    = .FALSE.
      TWOCOMPBSS = .FALSE.
      SKIP2E     = .FALSE.
      LOWJACO    = .FALSE.
      LOGMEMGET  = .FALSE.
      DOQJACO    = .FALSE.
      DOIMPMOS   = .FALSE.
C
      SET_CVALUE = .FALSE.
C
      CODSET     = 'NOCODATA'
      CALL SET_CODATA_VALUES(CODSET)
C
      RKBIMP_PRJTHR = 2.0D-5
C
C     Look for 2c Hamiltonians, ECPCALC, MDIRAC and more
C
      CALL HAMSCAN
C
C     Read menu file
C     ==============
C     **** Find General input *****
C
      IF (.NOT.FIND_KEYWORD('**GENER')) THEN
         NOGENINP = .TRUE.
         GOTO 180
      ENDIF
      WORD1 = '**GENERAL'
C
C     Process input for COMMON  /DCBGEN/
C     ==================================
C
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112,113,114,115,116,117,118,119,120,
     &                121,122,123,124),I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in GENINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in GENINP.')
  101    CONTINUE
C&&&& DIRECT: Direct calculation
C
C aspg, 2006-05-15: forcing flag for using gaunt interaction in the
C                   direct/conventional to zero if all other classes
C                   are set to zero (direct scf)... this is not the
C                   most satisfactory solution but at the moment allows
C                   for the use of conventional scf without changing
C                   the inputs.. it should be made consistent with
C                   the INTFLAG definition.
C
                 READ(LUCMD,*) ILLDIR,ISLDIR,ISSDIR
                 IF (ILLDIR.EQ.0.AND.ISLDIR.EQ.0.AND.ISSDIR.EQ.0) THEN
                    IGTDIR = 0
                 ENDIF
                 IDFLAG = ILLDIR+2*ISLDIR+4*ISSDIR+8*IGTDIR
                 DIRSET = .TRUE.
                 IF(IDFLAG.GT.0) THEN
                   DIRCAL = .TRUE.
                 ELSE
                   DIRCAL = .FALSE.
                 ENDIF
            GO TO 100
  102    CONTINUE
C&&&& CVALUE: Set a different value of speed of light (in a.u.)
            READ(LUCMD,*) CVAL
            SET_CVALUE = .TRUE.
            CSCALE = CVAL
            GO TO 100
  103    CONTINUE
C&&&& NOSET: Skip the PAMSET module
            NOSET = .TRUE.
            GO TO 100
  104    CONTINUE
C&&& DOJACO: Do Jacobi diagonalization if real group
C            (to avoid mixing of symmetries)
               DOJACO = .TRUE.
            GO TO 100
  105    CONTINUE
C&&&& PRINT : General print level
            READ(LUCMD,*) IPRGEN
C           transfer IPRGEN to IPRUSR originally defined in Dalton
C           and used in some abacus routines. /hjaaj aug 2005
            IPRUSR = IPRGEN
            USRIPR = .TRUE.
            GO TO 100
  106    CONTINUE
C&&&& PRJTHR
            READ(LUCMD,*) RKBIMP_PRJTHR
            GO TO 100
  107    CONTINUE
C&&&& PCMOUT: Write coefficients to formatted DFPCMO
            DOPUT = .TRUE.
            GO TO 100
  108    CONTINUE
C&&& LINDEP: Thresholds for linear dependence (eigenvalues of overlap matrix)
            READ(LUCMD,*) STOL(1),STOL(2)
            USER_STOL = .TRUE.
            GO TO 100
  109    CONTINUE
C&&&& SPHTRA: Embed transformation to spherical harmonics in
C&&&&         MO-transformation
C     This is transformation to
C
C     1. Spherical harmonics for Large components
C     2. Scalar restricted kinetic balance type for Small components
C
            READ(LUCMD,*) ILSPH,ISSPH
            ISPHTR = ILSPH+2*ISSPH
            GO TO 100
  110     CONTINUE
            GO TO 100
  111     CONTINUE
          LOGMEMGET = .TRUE.
            GO TO 100
  112     CONTINUE
!           Jacobi diagonalization of quaternion matrixes
            DOQJACO=.TRUE.
            GO TO 100
  113    CONTINUE
C&&&& ACMOUT: Print coefficients without symmetry (the C1 group).
            DOACUT = .TRUE.
            GO TO 100
  114    CONTINUE
C&&&& ACMOIN: Import coefficients in C1 to actual symmetry
            DOACIN = .TRUE.
            GO TO 100
  115    CONTINUE
C&&&& IMPMOS
            DOIMPMOS = .TRUE.
            READ(LUCMD,*) NFILES
            IF (NFILES.GT.2 .OR. NFILES.LT.1)  THEN
             CALL QUIT('1 or 2 MO files to import!')
            END IF
            READ(LUCMD,*) QCPACK, (FILIMP(IDEG),IDEG=1,NFILES)
            GO TO 100
  116    CONTINUE
C&&&& FAMILY: Family basis (Not yet implemented)
            FAMILY = .TRUE.
            GO TO 100
  117    CONTINUE
C&&&& SKIP2E
            SKIP2E = .TRUE.
            GO TO 100
  118    CONTINUE
C&&&& LOWJAC
            LOWJACO = .TRUE.
            GO TO 100
  119    CONTINUE
C&&&& RKBIMP: import RKB coefficients
            DORKBIMP=.TRUE.
            GO TO 100
  120    CONTINUE
C&&&& CODATA: Choose between CODATA sets
            READ(LUCMD,*) CODSET
            CALL UPCASE(CODSET)
            CALL SET_CODATA_VALUES(CODSET)
            IF(SET_CVALUE) THEN
             CVAL = CSCALE
            ELSE
             CVAL = CVEL
            END IF
            GO TO 100
  121    CONTINUE
            GO TO 100
  122    CONTINUE
C&&&& XXXXXX
            GO TO 100
  123    CONTINUE
C&&&& XXXXXX
            GO TO 100
  124    CONTINUE
C&&&& XXXXXX
            GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 180
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in GENINP.')
      END IF
  180 CONTINUE
C
      IF (PARCAL) THEN
ctec we require that it is a direct calculation
         DIRCAL = .TRUE.
         ILLDIR = 1
         ISLDIR = 1
         ISSDIR = 1
         IGTDIR = 1
         IDFLAG = 15
      ENDIF

!     InteRest uses a different normalization of small component
!     from Hermit, we must therefore correct STOL(2) in order to
!     get the desired threshold for linear dependence in LOWDIN module.

      if(mdirac .and. .not.USER_STOL)then
!       STOL(2) = STOL(1) / (2.0D0 * CVAL ** 2 )
        STOL(2) = 1.0d-04 / (2.0D0 * CVAL ** 2 )
      end if

C
C     process IMPMOS input
C

      IF (DOIMPMOS) THEN

       QCPACK = ADJUSTL(QCPACK)
       FILIMP(1) = ADJUSTL(FILIMP(1))

       if (nfiles.eq.2) then
        FILIMP(2) = ADJUSTL(FILIMP(2))
       else
        FILIMP(2) = ' -- unknown file name -- '
       end if

       CALL UPCASE(QCPACK)

       SELECT CASE (QCPACK(1:10))
        CASE ('MOLCAS    ')
         WRITE (LUPRI,'(/,A)') 
     &    "import MO coefficients from"
         WRITE (LUPRI,'(/,A,/)') 
     &    "  MOLCAS file "//TRIM(FILIMP(1))
        CASE ('RESPECT   ')
         WRITE (LUPRI,'(/,A)') 
     &    "import MO coefficients from"
         WRITE (LUPRI,'(/,A,/)') 
     &    "  ReSpect file "//TRIM(FILIMP(1))
        CASE ('TURBOMOLE ')
         WRITE (LUPRI,'(/,A)') 
     &    "import MO coefficients from"
         WRITE (LUPRI,'(/,A,/)') 
     &    "  Turbomole file "// TRIM(FILIMP(1))
        CASE DEFAULT
         CALL QUIT("unknown QCPACK type! Available options: "// 
     &             "MOLCAS, ReSpect")
       END SELECT

      ELSE
       QCPACK    = 'DIRAC '
       FILIMP(:) = ' -- unknown file name -- '
      END IF

C
C     Print section
C     =============
C
      CALL TITLER('General DIRAC set-up','*',126)
C
      call print_codata_reference(CODSET)

      WRITE(LUPRI,'(1X,A,F18.7)')
     +   '* The speed of light : ',CVEL
      IF(CVAL.NE.CVEL) THEN
        FMT = MXFORM(CVAL,14)
        WRITE(LUPRI,'(3X,A,1P,D9.2,0P,A,'//FMT//')')
     +   '* The speed of light scaled by ',(CVAL/CVEL),
     +   ' to : ',CVAL
      ENDIF

      IF (TWOCOMP) THEN
        WRITE(LUPRI,'(A)')
     &  ' * Running in two-component mode'
        INTGEN = 1
        ISLDIR = 0
        ISSDIR = 0
        IGTDIR = 0
      ELSE
        WRITE(LUPRI,'(A)')
     &  ' * Running in four-component mode'
      ENDIF

      IF(PARCAL) WRITE(LUPRI,'(A,I0,A)')
     &   ' * Parallel run with ',NUMNOD,' slaves.'
      IF(IDFLAG.GT.0) THEN
        WRITE(LUPRI,'(1X,A)') '* Direct evaluation of '//
     +    'the following two-electron integrals:'
        IF(ILLDIR.EQ.1) WRITE(LUPRI,'(3X,A)') '- LL-integrals'
        IF(ISLDIR.EQ.1) WRITE(LUPRI,'(3X,A)') '- SL-integrals'
        IF(ISSDIR.EQ.1) WRITE(LUPRI,'(3X,A)') '- SS-integrals'
        IF(IGTDIR.EQ.1) WRITE(LUPRI,'(3X,A)') '- GT-integrals'
      ENDIF
C
      IF(LBIT(ISPHTR,1)) THEN
        WRITE(LUPRI,'(1X,A)')
     &   '* Spherical transformation embedded in MO-transformation',
     &   '  for large components'
        IF(ISPHTR.EQ.3) THEN
           WRITE(LUPRI,'(1X,A)')
     &          '* Transformation to scalar RKB basis embedded in',
     &          '  MO-transformation for small components'
        ENDIF
      ELSEIF(ISPHTR.EQ.2) THEN
         WRITE(LUPRI,'(/1X,A)')
     &    'RKB with pre-projection not meaningful with non-spherical'//
     &    ' L component basis !'
         CALL QUIT('GENINP: Pre-RKB with non-spherical L-comp!')
      ENDIF
      WRITE(LUPRI,'(1X,A,2(/3X,A,1P,D9.2))')
     +     '* Thresholds for linear dependence:',
     +     'Large components:  ',STOL(1),
     +     'Small components:  ',STOL(2)

      IF (DOJACO.AND.DOQJACO) DOJACO=.FALSE.

      IF(DOJACO) WRITE(LUPRI,'(1X,A)')
     &    '* Matrix diagonalization by the Jacobi method '//
     &    'for real groups.'

      IF(DOQJACO) WRITE(LUPRI,'(1X,A)')
     &    '* Matrix diagonalization by the quaternion Jacobi method '//
     &    'for real, complex and quaternion groups.'

      IF(DOPUT) WRITE(LUPRI,'(1X,A)')
     &     '* MO-coefficients written to formatted file DFPCMO'
      IF(DOACUT) WRITE(LUPRI,'(1X,A)')
     &     '* MO-coefficients dumped in C1 format to'//
     &     ' unformatted file DFACMO'
      IF(DOACIN) WRITE(LUPRI,'(1X,A)')
     &     '* MO-coefficients imported in C1 format from'//
     &     ' unformatted file DFACMO'
      IF(DORKBIMP) THEN
        WRITE(LUPRI,'(1X,A,A)')
     &    '* Read RKB coefficients ',
     &    'and extended with complementary UKB space'
        WRITE(LUPRI,'(3X,A,E9.2)')
     &    ' Projection threshold: ',RKBIMP_PRJTHR
      ENDIF
      WRITE(LUPRI,'(1X,A,I3)')
     &     '* General print level   : ',IPRGEN
      IF(NOSET) WRITE(LUPRI,'(1X,A)')
     &     '* Skipping the PAMSET module. Crossing my fingers ...'
C
C     Process input for various program sections
C     ==========================================
C
      IF (NOGENINP) GOTO 300
  200 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4,5,6,7,8), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in GENINP.')
      ELSE
         WRITE (LUPRI,'(/,4A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in GENINP, error in prompt.')
      END IF
    2 CONTINUE
#if defined (VAR_MPI)
      CALL PRLINP(WORD)
#else
      CALL QUIT('*PARALLEL in input, but this is a non_MPI version')
#endif
        GO TO 200
    3 CONTINUE
        GO TO 200
    4 CONTINUE
        GO TO 200
    5 CONTINUE
        GO TO 200
    6 CONTINUE
        GO TO 200
    7 CONTINUE
        GO TO 200
    8 CONTINUE
        GO TO 200
C
    1 CONTINUE
  300 CONTINUE
#if defined (VAR_MPI)
      CALL PRLINP(WORD)
#endif
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck jobinp */
      SUBROUTINE JOBINP
C***********************************************************************
C
C <<< Job input for DIRAC >>>
C
C     Written by Trond Saue May 1996
C     Last revision: 1997/09/15 - jth (OPTIMIZE added)
C
C***********************************************************************
      use xmlout
      use dirac_cfg
#include "implicit.h"
#include "priunit.h"
      PARAMETER (NTABLE = 14, NDIR = 4)
      LOGICAL FIND_KEYWORD
C
      CHARACTER WORD*7, PROMPT*1, TABLE(NTABLE)*7, TABDIR(NDIR)*7,
     &          WORD1*7
C
#include "dcbgen.h"
C
      DATA TABDIR /'*END OF','*OPTIMI','*MINIMI','xXXXXXX'/
C
      DATA TABLE  /'.TITLE ','.INPTES','.4INDEX','.ONLY I',
     &             '.WAVE F','.ANALYZ','.PROPER','xxxxxx ',
     &             '.NO4IND','.NOSFCR','.OPTIMI','.MINIMI',
     &             '.XMLOUT','.xxxxxx'/
c     .XMLOUT can be removed after adapting the tests, is now always active
C
      NOSFCRASH=.FALSE.
C
C     Read menu file
C     ==============
C
C     **** Locate **DIRAC input
      IF (.NOT.FIND_KEYWORD('**DIRAC'))
     &   CALL QUIT('No **DIRAC input found')
      WORD1 = '**DIRAC'
C
C     Process job input for COMMON  /CBIPAM/
C     ======================================
C
      NJOBS = 0
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112,113,114),I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in JOBINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in JOBINP.')
  101    CONTINUE
C&&&& TITLE: Title of run
            READ(LUCMD,'(A50)') TITLE
            GO TO 100
  102    CONTINUE
C&&&& INPTES: Input test
            INPTES = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  103    CONTINUE
C&&&& activate 4-INDEX transformation (on by default for CC and CI)
            DOTRA = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  104    CONTINUE
C&&&& ONLY INTEGRALS: Only run integrals
            DOHRM = .TRUE.
            GO TO 100
  105    CONTINUE
C&&&& WAVE FUNCTION: call wavefunction module
            DOPSI = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  106    CONTINUE
C&&&& ANALYZE : call analysis module
            DOANA = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  107    CONTINUE
C&&&& PROPERTIES:  Call property module
            DOPRP = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  108    CONTINUE
C&&&
            GO TO 100
  109    CONTINUE
C&&&& deactivate 4-INDEX transformation (on by default for CC and CI)
            NOTRA = .TRUE.
            GO TO 100
  110    CONTINUE
C&&& .NOSFCR: no Segm.fault crash when quit for creating profiling output,gmon.out
            NOSFCRASH=.TRUE.
            GO TO 100
 111     CONTINUE
C&&&& OPTIMIZE: Geometry optimization
            OPTIMI = .TRUE.
            DOPSI  = .TRUE. ! no geometry optimization without a wave function
            GOTO 100
 112     CONTINUE
C&&&& MINIMIZE: Geometry optimization (tsaue: retained for backwards compatibility)
            GO TO 111
 113     CONTINUE
C&&& XMLOUT  : Retained for backwards compatibility (lv: is now always on)
         GOTO 100
!       xxxxxx keyword
 114     continue
         goto 100

      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 180
      ELSE
         WRITE (LUPRI,'(/,3A,A,/)') ' Prompter "',PROMPT,
     &    '" illegal in input word:',WORD
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in JOBINP.')
      END IF
C
C     Print section
C     ==============
C
  180 CONTINUE
      I  = LNBLNK(TITLE)
      IN = (272-I)/2
      CALL TITLER(TITLE(1:I),'*',IN)
      IF(NJOBS.GT.0) WRITE(LUPRI,'(A)') ' Jobs in this run:'
      IF(DOHRM) THEN
        WRITE(LUPRI,'(3X,A)') '* Only run one-electron integrals'
        GOTO 666
      ENDIF
      IF(OPTIMI)  WRITE(LUPRI,'(3X,A)') '* Geometry optimization'
      IF(DOPSI)   WRITE(LUPRI,'(3X,A)') '* Wave function'
      IF(DOANA)   WRITE(LUPRI,'(3X,A)') '* Analysis'
      IF(DOPRP)   WRITE(LUPRI,'(3X,A)') '* Properties'
      IF(DOTRA)   WRITE(LUPRI,'(3X,A)') '* Transformation to '//
     &    'Molecular Spinor basis'
      IF(INPTES)  WRITE(LUPRI,'(3X,A)') '* Input test run only'

C     LV: This option is badly out of place, has nothing to do with jobs !
      IF (NOSFCRASH) THEN
         WRITE(LUPRI,'(3X,A)')
     &   '* If crash, do NOT end with segmentation fault.'//
     &   'Necessary for getting profiling file, gmon.out.'
!radovan: commented out the following message about this default situation
!         this is strange for users to see under jobs
!     ELSE
!        WRITE(LUPRI,'(3X,A)')
!    &   '* If crash, end with a forced segmentation fault.'//
!    &   ' Useful for debugging.'
      END IF
C
 666  CONTINUE
C
C     Process input for various program sections
C     ==========================================
C
  200 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/3A/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in JOBINP.')
      ELSE
         WRITE (LUPRI,'(/4A/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in JOBINP, error in prompt.')
      END IF
    2 CONTINUE
C       *OPTIMI
        CALL OPINPU(WORD)
        GO TO 200
    3 CONTINUE
Ctsaue: *MINIMI for backwards compatibility
        CALL OPINPU(WORD)
        GO TO 200
    4 CONTINUE
        GO TO 200
    1 CONTINUE
C     **<something> or *END OF
      CALL OPINPU(WORD)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SCFINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for SCF module
C
C     Written by T.Saue - November 1992
C     Last revision: 1997/09/15 - jth
C
C***********************************************************************
       use dirac_cfg
       use x2cmod_cfg
#ifdef MOD_AOOSOC
       use atomic_oo_order_so_correction_cfg
#endif

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (NTABLE = 65)
      PARAMETER (RTOL = 1.0D-15, D1 = 1.0D0, D0=0.0D0, D2=2.0D0,
     &           DTHRS=1.0D-13, DP5 = 0.50D0)
C
      LOGICAL SET, NEWDEF, LBIT, RESET, NOAUTOCC, INTFLG_CHANGE
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER LINE*80, TEXT*20, CTEMP*72, temp_line*80
      DIMENSION WORK(LWORK)
      integer, external :: word_count
      INTEGER   KPMAX(2)
      INTEGER, allocatable :: SEL_NORB(:,:), KP_NORB(:)
C
#include "nuclei.h"
#include "dcbgen.h"
#include "huckel.h"
#include "dcbpsi.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbatom.h"
#include "dcbhoc.h"
#include "maxaqn.h"
#include "ccom.h"
C
      SAVE SET
      DATA TABLE /'.PROJEC','.PRINT ','.MAXITR','.ERGCNV','.EVCCNV',
     &            '.FCKCNV','.NODSCF','.NOSWIT','.FIXDIF','.NODAMP',      ! 10
     &            '.NODIIS','.MXDIIS','.CNVINT','.ITRINT','.INTFLG',
     &            '.MOSTAR','.XXXXXX','.CLOSED','.EIGPRI','.OPEN S',      ! 20
     &            '.AUTOCC','.OWNBAS','.NOBNCR','.DIISTH','.DAMPFC',
     &            '.OVLSEL','.NODYNS','.BNCORR','.DIISMO','.DIISAO',      ! 30
     &            '.XXXXXX','.MAX MA','.MAX MI','.2NDOPT','.SKIPEE',
     &            '.SKIPEP','.MAX BA','.PRJTHR','.FROZEN','.BOSONS',      ! 40
     &            '.MJSELE','.AOC   ','.FOCC  ','.LSHIFT','.MAXIT2',
     &            '.ERGCN2','.EVCCN2','.FCKCN2','.MOFREE','.OLEVEL',      ! 50
     &            '.ATOMST','.SCFPOP','.OPENFA','.BNC_FO','.FOMOUT',
     &            '.AD HOC','.HUCPAR','.KPSELE','.SELECT','.XXXXXX',      ! 60
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'
     &                                                               /
      DATA SET/.FALSE./
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
C *** Initialize in dcbdhf.h ***
C
C 1.  Initialize /DCCDHF/
C     ===================
C
      DHF_INTTYP= '            '
      CACC      = '        '
C
C 2.  Initialize /DCLDHF/
C     ===================
C
      SCFPOP    = .FALSE.
      ERGCNV    = .FALSE.
      ERGCNV2   = .FALSE.
      EVCCNV    = .TRUE.
      EVCCNV2   = .TRUE.
      FCKCNV    = .FALSE.
      FCKCNV2   = .FALSE.
      DHFCONV(1)= .FALSE.
      DHFCONV(2)= .FALSE.
      DHFEXIT   = .FALSE.

C     initial molecular orbitals:
      TRIVEC    = .TRUE. ! default!
      TRIFCK    = .FALSE.
      BARNUC    = .TRUE.
      BNCRON    = .FALSE.
      FORCE_BNC = .false.
      BNSPON    = .TRUE.  ! use either this or BNCRON
      SIRIFC    = .FALSE.
      DOHUCKEL  = .false.  
      EWMO      = .FALSE. ! in huckel.h
      HUCPROJCMO  = .FALSE. ! in huckel.h, projection method

      DODSCF    = .FALSE.
      FIXDIF    = .FALSE.
      IF(IDFLAG.GT.0) DODSCF    = .TRUE.
      DOCCNV    = .TRUE.
      DODAMP    = .TRUE.
      DODIIS    = .TRUE.
      NOSWIT    = .FALSE.
      DIISON    = .FALSE.
      OVLSEL    = .FALSE.
      DOMOFREEZE = .FALSE.
      DYNSEL    = .FALSE.
      DOBOSSEL  = .FALSE.
      NSUBBL    = NBSYM
      IBOSSEL   = 0
C     ... hjaaj: default = # boson sym.s
C         must be changed with .MJSELE for LINEAR
C         this can not be tested if correct because
C         GLINSM in PAMSET has not been called yet. /19.Oct.04
      DIISAO    = .FALSE.
      DIISMO    = .TRUE.
      AUTOCC    = .FALSE.
      INIOCC    = .TRUE.
      PROJEC  = .FALSE.
      PROOWN  = .FALSE.
      NOQCDHF   = .TRUE.
C     ... default do QCDHF if DIIS does not converge
      DHF_SKIPEE = .FALSE.
      DHF_SKIPEP = .FALSE.
!     if number of positronic shells is zero skip e-p rotations. hjj+sk - aug 2010
      if(x2c.or.bss.or.levyle.or.freepj.or.vextpj)
     &  DHF_SKIPEP = .TRUE.
C     ... Average-of-Configuration (AOC) is default for Hartree-Fock,
C     ... whereas Fractional Occupation is default for DFT
      AOC = .TRUE.
      IF (dirac_cfg_dft_calculation) AOC = .FALSE.
!     TODO: fix oo-order 2e-SOC-SC corrections for AOC - seems to work for FOCC - sk march 2014
#ifdef MOD_AOOSOC
      if (aoomod .and. x2cmod_mmf) AOC = .FALSE.
#endif
      ATOMST = .FALSE.
      ATHUCK = .FALSE.
      HUCPAR = 1.5D0
      WRITE_FMO_MATRIX=.FALSE.
C
C 3.  Initialize /DCIDHF/
C     ===================
C
      MAXITR    = 50
CMI  ... reduce the number iterations for preliminary BSS-SCF ...
      MAXITR2   = 25
      IPRSCF    = IPRGEN
      MXDIIS    = 10
      NBELM     = 1
      ITDIIS    = 0
      ILLINT    = IBTAND(INTGEN,1)
      ISLINT    = IBTAND(INTGEN/2,1)
      ISSINT    = IBTAND(INTGEN/4,1)
      IGTINT    = IBTAND(INTGEN/8,1)
      INTBUF    = 0
      INTFLG    = 0
      INTFLG_CHANGE = .FALSE.
      ITRINT(1) = 1
      ITRINT(2) = 1
      IPREIG    = 1
      DHF_INTTYP= 'Bare nuclei '
      NEWOCC  = 0
      MXMACRO = 25
      MXMICRO = 25
      MAXBCK  = 5
C
C 4.  Initialize /DCRDHF/
C     ===================
C
      SCFCNV(1) = 1.0D-7
      SCFCNV(2) = 1.0D-6
CMI   ... criteria of convergence for prelim.BSS-SCF
      SCFCNV2(1) = 1.0D-5
      SCFCNV2(2) = 1.0D-4
      ERGVAL    = DUMMY
      EVCVAL    = D0
      FCKVAL    = D0
      DHFERG    = D0
      ELERGY    = D0
      E1PART    = D0
      E2PART    = D0
      ERGBUF    = D0
      CONVRG    = D0
      TDF2      = D0
      TDDG      = D0
      BMCOND    = D0
      DAMPFC    = 0.25D0
      DIISTH    = DUMMY
      CNVINT(1) = DUMMY
      CNVINT(2) = DUMMY
      PRJTHR    = 1.0D-10
      DOLEVEL   = .FALSE.
      DLSHIF    = D0
      OPEN_FAC  = -0.5D0 ! abs(OPEN_FAC) seems to a good compromise value for good convergence
      ! the negative sign signifies that this is the default value, i.e. not user input
      CALL DZERO(OLEV,MXOPEN)
      CALL DZERO(SCFTID,8)
      CALL IZERO(ITRSCF,8)
C
C 5.  Initialize /DCODHF/
C     ===================
C
      NELEC_DHF(1) = 0
      NELEC_DHF(2) = 0
      NISH_DHF(1)  = 0
      NISH_DHF(2)  = 0
      NASH_DHF(1)  = 0
      NASH_DHF(2)  = 0
      CALL IZERO(NISH_BOS,MAX_BOS_BL)
      CALL IZERO(NACSH_BOS,MAX_BOS_BL*MXOPEN)
      CALL IZERO(NACSH,2*MXOPEN)
      NOPEN    = 0
      NFRO_DHF(1)  = 0
      NFRO_DHF(2)  = 0
      CALL DZERO(DA,1+MXOPEN)
      CALL DZERO(DF,1+MXOPEN)
      CALL DZERO(DALPHA,1+MXOPEN)
C
C     Process input for *SCF
C     ========================


      newdef = (word == '*SCF   ')
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
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24,25,26,27,28,29,30,
     &                     31,32,33,34,35,36,37,38,39,40,
     &                     41,42,43,44,45,46,47,48,49,50,
     &                     51,52,53,54,55,56,57,58,59), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/4A/)') ' Keyword "',WORD,
     *            '" not recognized for ',WORD1
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in *SCF.')
    1          CONTINUE
C&&&& PROJECtion: Reduce transformation matrix by projection
               PROJEC = .TRUE.
               READ(LUCMD,*) NFRAG
               IF(NFRAG.GT.MFRAG) THEN
                 WRITE(LUPRI,'(A,I5)')
     &           '* Number of fragments specified:',NFRAG,
     &           '* Current maximum              :',MFRAG
                 CALL QUIT('*SCF: Too many .PROJEC fragments !')
               ENDIF
               DO J = 1,NFRAG
                 READ(LUCMD,'(A6)') PRJFIL(J)
                 READ(LUCMD,'(I6)') NPRJNUC(J)
                 READ(LUCMD,'(A72)') (VCPROJ(I,J),I=1,NFSYM)
               ENDDO
               GO TO 100
    2          CONTINUE
C&&& PRINT : General print level in DHF-module
                  READ(LUCMD,*) IPRSCF
               GO TO 100
    3          CONTINUE
C&&& MAXITR: Maximum number of SCF-iterations
                  READ(LUCMD,*) MAXITR
               GO TO 100
    4          CONTINUE
C&&& ERGCNV: Convergence on energy difference
                  READ(LUCMD,'(A80)') LINE
                  CALL UPCASE(LINE)
                  READ(LINE,*,END=401,ERR=401) SCFCNV(1),SCFCNV(2)
                  GOTO 402
  401             READ(LINE,*) SCFCNV(1)
                  SCFCNV(2) = SCFCNV(1)
  402             CONTINUE
                  EVCCNV = .FALSE.
                  ERGCNV = .TRUE.
                  FCKCNV = .FALSE.
               GO TO 100
    5          CONTINUE
C&&& EVCCNV: Convergence on approximate gradient (error vector)
                  READ(LUCMD,'(A80)') LINE
                  CALL UPCASE(LINE)
                  READ(LINE,*,END=501,ERR=501) SCFCNV(1),SCFCNV(2)
                  GOTO 502
  501             READ(LINE,*) SCFCNV(1)
                  SCFCNV(2) = SCFCNV(1)
  502             CONTINUE
                  EVCCNV = .TRUE.
                  ERGCNV = .FALSE.
                  FCKCNV = .FALSE.
               GO TO 100
    6          CONTINUE
C&&& FCKCNV: Convergence on total Fock matrix
                  READ(LUCMD,'(A80)') LINE
                  CALL UPCASE(LINE)
                  READ(LINE,*,END=601,ERR=601) SCFCNV(1),SCFCNV(2)
                  GOTO 602
  601             READ(LINE,*) SCFCNV(1)
                  SCFCNV(2) = SCFCNV(1)
  602             CONTINUE
                  EVCCNV = .FALSE.
                  ERGCNV = .FALSE.
                  FCKCNV = .TRUE.
               GO TO 100
    7          CONTINUE
C&&&& NODSCF: Do not employ differential matrices in SCF-iterations
                  DODSCF = .FALSE.
               GO TO 100
    8          CONTINUE
C&&&& NOSWIT: PBS+MI Never turn DIIS off
                  NOSWIT=.TRUE.
               GO TO 100
    9          CONTINUE
C&&&& FIXDIF: Use fixed value of DIFDEN
                  FIXDIF = .TRUE.
               GO TO 100
   10          CONTINUE
C&&& NODAMP: No damping of Fock matrix
                  DODAMP = .FALSE.
               GO TO 100
   11          CONTINUE
C&&& NODIIS: Do not perform DIIS iterations
                  DODIIS = .FALSE.
               GO TO 100
   12          CONTINUE
C&&& MXDIIS: Maximum size of B-matrix in DIIS
                  READ(LUCMD,*) MXDIIS
                  DODIIS = .TRUE.
               GO TO 100
   13          CONTINUE
C&&& CNVINT: Convergence thresholds for adding SL- and SS-integrals
                  READ(LUCMD,*) CNVINT(1),CNVINT(2)
               GO TO 100
   14          CONTINUE
C&&& ITRINT: Number of iterations before adding SL- and SS-integrals
                  READ(LUCMD,*) ITRINT(1),ITRINT(2)
               GO TO 100
   15          CONTINUE
C&&& INTFLG: Specify what two-itegrals should be included in this run
                  IF (IGTINT.EQ.1) THEN
!                    ... if gaunt=true
                     READ(LUCMD,*,IOSTAT=IOS)
     &               ILLINT,ISLINT,ISSINT,IGTINT
                     IF (IOS.NE.0) THEN
                        WRITE(LUPRI,'(/,2X,A)')
     &                  'Error in *SCF .INTFLG reading !'//
     &                  ' 4 parameters needed for Gaunt term !'
                        CALL FLSHFO(LUPRI)
                        CALL QUIT( '*SCF: Error in .INTFLG'//
     &                  'reading for GAUNT=true!')
                     ENDIF
                  ELSE
                     READ(LUCMD,*,IOSTAT=IOS) ILLINT,ISLINT,ISSINT
                     IF (IOS.NE.0) THEN
                       WRITE(LUPRI,'(/,2X,A)')
     &                 'Error in *SCF .INTFLG reading !'//
     &                 ' 3 integer parameters needed  !'
                       CALL FLSHFO(LUPRI)
                       CALL QUIT('*SCF: Error in .INTFLG reading')
                     ENDIF
                  END IF
                  INTFLG_CHANGE=.TRUE.
               GO TO 100
   16          CONTINUE
C&&& MOSTARt: How to obtain initial MOs [ default: 1) TRIVEC, 2) BNCORR ]
               DOHUCKEL = .FALSE.
               BARNUC = .FALSE.
               TRIVEC = .FALSE.
               TRIFCK = .FALSE.
               SIRIFC = .FALSE.
! hjaaj March 2010: discovered gfortran only reads first word in line
! with '*' format  (Ex: '   BARE NUC' became 'BARE       ' in LINE !)
!              READ (LUCMD,*) LINE
               READ (LUCMD,'(A)') LINE
               CALL UPCASE(LINE)
               IF (INDEX(LINE,'HUCKEL').ne.0 .OR.
     &             INDEX(LINE,'EHT').ne.0) THEN
                  DOHUCKEL = .TRUE.
                  EWMO     = .FALSE.
               ELSE IF (INDEX(LINE,'HUCPRO').ne.0 .OR.
     &                  INDEX(LINE,'EHTPRO').ne.0) THEN
                  DOHUCKEL   = .TRUE.
                  HUCPROJCMO = .TRUE.
                  EWMO       = .FALSE.
               ELSE IF (INDEX(LINE,'EWMO').ne.0) THEN
                  DOHUCKEL = .TRUE.
                  EWMO     = .TRUE.
                  IF (NBSYM.GT.1) THEN
                     EWMO=.FALSE.
                     WRITE(LUPRI,*) ' EWMO deactivated for NBSYM.gt.1'
                  END IF
               ELSE IF (INDEX(LINE,'EQMOPR').ne.0 .OR.
     &                  INDEX(LINE,'EWMOPR').ne.0) THEN
                  DOHUCKEL   = .TRUE.
                  HUCPROJCMO = .TRUE.
                  EWMO       = .TRUE.
                  IF (NBSYM.GT.1) THEN
                     EWMO=.FALSE.
                     WRITE(LUPRI,*) ' EWMO deactivated for NBSYM.gt.1'
                  END IF
               ELSE IF (INDEX(LINE,'BARNUC').ne.0 .OR.
     &                  INDEX(LINE,'BARE NUC').ne.0) THEN
                  BARNUC = .TRUE.
                  BNCRON = .FALSE.
                  BNSPON = .FALSE.
               ELSE IF (INDEX(LINE,'BNCORR').ne.0) THEN
                  BARNUC = .TRUE.
                  BNCRON = .TRUE.
                  BNSPON = .FALSE.
               ELSE IF (INDEX(LINE,'SCRPOT').ne.0) THEN
                  BARNUC = .TRUE.
                  BNCRON = .FALSE.
                  BNSPON = .TRUE.
               ELSE IF (INDEX(LINE,'TRIFCK').ne.0) THEN
C              TRIFCK: Start on two-electron Fock matrix read from file
                  TRIFCK = .TRUE.
               ELSE IF (INDEX(LINE,'TRIVEC').ne.0) THEN
C              TRIVEC: Read trial vectors from DFCOEF file
                  TRIVEC = .TRUE.
               ELSE IF (INDEX(LINE,'SIRIFC').ne.0) THEN
                  SIRIFC = .TRUE.
               ELSE
                  WRITE(LUPRI,'(/A/A)') ' FATAL ERROR: '//
     &            '.MOSTART option not recognized. Option line was:',
     &            LINE
                  CALL QUIT('.MOSTART option not recognized')
               END IF
               GO TO 100
   17          CONTINUE
               GO TO 100

   18          CONTINUE
!                 .CLOSED
!                 number of closed shell electrons

                  ! sanity check: stop the code if .CLOSED does not match NFSYM
                  read(lucmd, '(a80)', iostat=ios) temp_line
                  backspace lucmd
                  if (word_count(temp_line) /= nfsym) then
                     call quit('.CLOSED does not match NFSYM')
                  end if

                  READ (LUCMD,*,IOSTAT=IOS) (NELEC_DHF(I),I=1,NFSYM)
                  IF (IOS.NE.0) THEN
                    WRITE(LUPRI,*)
     & ' *SCF: Error in reading .CLOSED occupation for NFSYM=',NFSYM
                    WRITE(LUPRI,*) 'Values of NELEC_DHF(I),I=1,NFSYM:',
     &              (NELEC_DHF(I),I=1,NFSYM)
                   CALL QUIT(
     &             '*SCF: Error in reading .CLOSED occupation')
                  ENDIF
                  INIOCC = .FALSE.
               GO TO 100

   19          CONTINUE
C&&&& EIGPRI: Print control for electron and positron solutions
                  READ(LUCMD,*,IOSTAT=IOS) IL,IS
                  IF (IOS.NE.0) THEN
                    CALL QUIT('*SCF ERROR: .EIGPRI input error')
                  ENDIF
                  IPREIG = IL+2*IS
               GO TO 100
   20          CONTINUE
C&&& OPEN SHELLS: read in number of active shells and their occupation
C
C             a) Read in number of open shells
C
               READ(LUCMD,*,IOSTAT=IOS) NOPEN
               IF (IOS.NE.0) THEN
                 CALL QUIT('*SCF: .OPEN SHELL input error !')
               ENDIF
               IF (NOPEN.LE.0.OR.NOPEN.GT.MXOPEN) THEN
                  WRITE(LUPRI,'(/A,I2)')
     &            ' *SCF: .OPEN SHELLS must be between 1 and',MXOPEN
                  CALL QUIT('*SCF ERROR: Open shell error!')
               END IF
C
C              b) Read in number of active shells and their occ.
C              Format:
C              If inversion symmetry: A/K,L
C              If not               : A/K
C
               DO IOPEN = 1,NOPEN
                  READ(LUCMD,'(A80)') LINE
                  CALL UPCASE(LINE)
C
                  ISLASH = INDEX(LINE,'/')
                  IF (ISLASH .LE. 1) THEN
                     WRITE (LUPRI,'(/A,I2,A,I2/A,A)')
     &               'ERROR for *SCF .OPEN SHELL shell no.',
     &                IOPEN,' out of',NOPEN,
     &               '- the input line is:',LINE
                     IF (ISLASH .LE. 0) THEN
                        WRITE(LUPRI,'(A)')
     &                '- error: the input line does not contain a "/"'
                     ELSE
                       WRITE(LUPRI,'(A)')
     &                '- error: no number of electrons before the "/"'
                     END IF
                     CALL QUIT(
     &               'Input error for .OPEN SHELL under *SCF')
                  END IF
                  TEXT = LINE(1:ISLASH-1)
                  READ (TEXT,*) FELEC
C
                  ICOMMA= INDEX(LINE,',')

                  IF (ICOMMA .GT. 0) THEN
                     IF (ICOMMA .LE. ISLASH+1) THEN
                        CALL QUIT('ERROR: Misplaced /')
                     END IF
                     TEXT = LINE(ISLASH+1:ICOMMA-1)
                     READ (TEXT,'(I20)') IAC1
                     TEXT = LINE(ICOMMA+1:)
                     READ (TEXT,'(I20)') IAC2
                  ELSE
                     TEXT = LINE(ISLASH+1:)
                     READ (TEXT,'(I20)') IAC1
                     IAC2 = 0
                  END IF

                  IF (ICOMMA.EQ.0.AND.NFSYM.GE.2) THEN
                   WRITE(LUPRI,*)
     &            '*SCF: Error in open-shell input reading !'
                   WRITE(LUPRI,*) 'NOPEN:',NOPEN,' NFSYM:',NFSYM
                   WRITE(LUPRI,'(A,A)') 'input line:',LINE
                   WRITE(LUPRI,*)
     &  'for NFSYM=2 has to be of the form xx/xx,xx (see the manual)'
                   CALL QUIT(
     &             '**SCF: Error in open-shell input reading !')
                  ELSE IF (ICOMMA.GT.0.AND.NFSYM.EQ.1) THEN
                   WRITE(LUPRI,*)
     &              '**SCF: Error in open-shell input reading'
                   WRITE(LUPRI,*)
     &              'NFSYM is 1 (no inversion symmetry)'
                   WRITE(LUPRI,'(A,A)')
     &              'you have specified (N/G,U):',LINE
                    WRITE(LUPRI,*)
     &                'for NFSYM=1 use N/O instead of N/G,U'
                    CALL QUIT(
     &                '**SCF: Error in open-shell input reading !')
                  ENDIF

                  IF (NFSYM .EQ. 1) IAC2 = 0
                  FACTOT = real(IAC1+IAC2,8)
                  NERR = 0
                  IF (FELEC .LT. D0 .OR. FELEC .GT. FACTOT) THEN
                     NERR = NERR + 1
                     WRITE(LUPRI,'(A,F8.2,A,I3)')
     &                  '.OPEN SHELL error: Number of electrons (',
     &                  FELEC,') must be a positive number less than',
     &                  IAC1+IAC2
                  END IF
                  IF ( (MOD(IAC1,2) .NE. 0) .OR.
     &                 (MOD(IAC2,2) .NE. 0) .OR.
     &                 (IAC1 .LT. 0) .OR. (IAC2 .LT. 0) .OR.
     &                 (IAC1+IAC2 .EQ. 0)) THEN
                     NERR = NERR + 1
                     WRITE(LUPRI,'(A,2I3)') '.OPEN SHELL error: '//
     &                  'Number of spinors must be even'//
     &                  ' numbers and at least 2 in total:',IAC1,IAC2
                  END IF
                  IF (NERR .GT. 0) THEN
                     CALL QUIT('Error in .OPEN SHELL specifications')
                  END IF
                  NACSH(1,IOPEN) = IAC1/2
                  NACSH(2,IOPEN) = IAC2/2
                  DF(IOPEN) = FELEC/FACTOT
               END DO
               INIOCC = .FALSE.
               GO TO 100
   21          CONTINUE
C&&& AUTOCC: Automatic DHF occupation
                  AUTOCC = .TRUE.
               GO TO 100
   22          CONTINUE
C&&&& OWNBASis: Use local basis for projection
                  PROOWN = .TRUE.
               GO TO 100
   23          CONTINUE
C&&&& NOBNCR : Deactivate screening for use with bare nucleus start guess
                  BNCRON = .FALSE.
                  BNSPON = .FALSE.
               GO TO 100
   24          CONTINUE
C&&& DIISTH: Change default threshold for DIIS
                  READ(LUCMD,*) DIISTH
               GO TO 100
   25          CONTINUE
C&&& DAMPFC: Change default damping factor
                  READ(LUCMD,*) DAMPFC
               GO TO 100
   26          CONTINUE
C&&& OVLSEL: Activate overlap selection
C            (default is dynamic overlap selection, see below)
                  OVLSEL = .TRUE.
                  DYNSEL = .TRUE.
               GO TO 100
   27          CONTINUE
C&&& NODYNSEL: No Dynamic overlap selection
C            (new overlap selection matrix formed in each iteration)
                  DYNSEL = .FALSE.
               GO TO 100
   28          CONTINUE
C&&& BNCORR: Bare nucleus screening (based on Slaters rules)
                  BNCRON = .TRUE.
                  BNSPON = .FALSE.
               GO TO 100
   29          CONTINUE
C&&& DIISMO: Do DIIS in MO-basis
                  DIISAO = .FALSE.
                  DIISMO = .TRUE.
               GO TO 100
   30          CONTINUE
C&&& DIISAO: Do DIIS in AO-basis
                  DIISAO = .TRUE.
                  DIISMO = .FALSE.
               GO TO 100
   31          CONTINUE
               GO TO 100
   32          CONTINUE
C&&& MAX MACRO ITERATIONS: Maximum number of macro iterations for QC-DHF
                  READ(LUCMD,*) MXMACRO
               GO TO 100
   33          CONTINUE
C&&& MAX MICRO ITERATIONS: Maximum number of micro iterations for QC-DHF
                  READ(LUCMD,*) MXMICRO
               GO TO 100
   34          CONTINUE
C&&& 2NDOPT: activate second-order optimization if regular SCF does not converge
                  NOQCDHF = .FALSE.
               GO TO 100
   35          CONTINUE
C&&& SKIPEE: no e-e rotations in 2nd order optimization
                  DHF_SKIPEE = .TRUE.
               GO TO 100
   36          CONTINUE
C&&& SKIPEP: no e-p rotations in 2nd order optimization
                  DHF_SKIPEP = .TRUE.
               GO TO 100
 37            CONTINUE
C&&& MAX BACKSTEPS: max. number of backsteps
                  READ(LUCMD,*) MAXBCK
               GOTO 100
 38            CONTINUE
C&&& PRJTHR: Smallest norm accepted in projection step
                  READ(LUCMD,*) PRJTHR
               GOTO 100
 39            CONTINUE
C&&& FROZEN: Frozen orbitals
               DO I = 1,NFSYM
                 READ(LUCMD,'(A72)') VCFROZ(I)
                 CALL NUMLS3(VCFROZ(I),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),NFRO_DHF(I),NZERO)
                 NFRO_DHF(I)=NFRO_DHF(I)-NZERO
               ENDDO
               GOTO 100
 40            CONTINUE
C&&& BOSONS: Selection based on (approx.) boson symmetry
               DOBOSSEL = .TRUE.
               IBOSSEL  = 1
!     READ(LUCMD,*) (NISH_BOS(I),I=1,NSUBBL) !old version
               allocate(SEL_NORB(NSUBBL,0:NOPEN))
               READ(LUCMD,*,IOSTAT=IOS) NSUBBL                
               READ(LUCMD,*) (SEL_NORB(I,0),I=1,NSUBBL)
               DO I =1,NSUBBL
                 IF(mod(SEL_NORB(I,0),2).eq.1)THEN
                   WRITE(LUPRI,*) 'Error in reading .BOSONS !'
                   CALL QUIT('Number of electrons must be even!')                    
                 ELSE   
                   NISH_BOS(I) = SEL_NORB(I,0)/2               
                 ENDIF
               ENDDO                 
               DO IOPEN = 1,NOPEN
!                 READ(LUCMD,*) (NACSH_BOS(I,IOPEN),I=1,NSUBBL) !old version
                 READ(LUCMD,*)(SEL_NORB(I,IOPEN),I=1,NSUBBL)
                 DO I =1,NSUBBL
                   IF(mod(SEL_NORB(I,IOPEN),2).eq.1)THEN
                     WRITE(LUPRI,*) 'Error in reading .BOSONS !'
                     CALL QUIT('Number of electrons must be even!')   
                   ELSE   
                     NACSH_BOS(I,IOPEN) = SEL_NORB(I,IOPEN)/2
                   ENDIF
                 ENDDO                    
               END DO
               deallocate(SEL_NORB)
               GOTO 100
 41            CONTINUE
C&&& MJSELE: Selection based on Mj symmetry,
C            NSUBBL = # of subblocks for MJ sym.
               DOBOSSEL = .TRUE.
               IBOSSEL  = -1
               READ(LUCMD,*,IOSTAT=IOS) NSUBBL
               IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*) 'Error in reading NSUBBL for .MJSELE !'
                 CALL QUIT('Error in reading NSUBBL for .MJSELE !')
               ENDIF
               allocate(SEL_NORB(NSUBBL,0:NOPEN))                             
!               READ(LUCMD,*,IOSTAT=IOS) (NISH_BOS(I),I=1,NSUBBL)
               READ(LUCMD,*,IOSTAT=IOS) (SEL_NORB(I,0),I=1,NSUBBL)
               DO I =1,NSUBBL
                 IF(mod(SEL_NORB(I,0),2).eq.1)THEN
                   WRITE(LUPRI,*) 'Error in reading .MJSELE !'
                   CALL QUIT('Number of electrons must be even!')                    
                 ELSE   
                   NISH_BOS(I) = SEL_NORB(I,0)/2               
                 ENDIF
               ENDDO
               IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*)
     &           'Error in reading NISH_BOS(I) for .MJSELE !'
                 CALL QUIT('Error in reading NISH_BOS(I) for .MJSELE !')
               ENDIF
               DO IOPEN = 1,NOPEN
!                 READ(LUCMD,*,IOSTAT=IOS)
!     &           (NACSH_BOS(I,IOPEN),I=1,NSUBBL)
                 READ(LUCMD,*,IOSTAT=IOS)
     &                (SEL_NORB(I,IOPEN),I=1,NSUBBL)                  
                 DO I =1,NSUBBL
                   IF(mod(SEL_NORB(I,IOPEN),2).eq.1)THEN
                     WRITE(LUPRI,*) 'Error in reading .MJSELE !'
                     CALL QUIT('Number of electrons must be even!')   
                   ELSE   
                     NACSH_BOS(I,IOPEN) = SEL_NORB(I,IOPEN)/2
                   ENDIF
                 ENDDO                  
                 IF (IOS.NE.0) THEN
                   WRITE(LUPRI,*)
     &             'Error in reading NACSH_BOS(I) for .MJSELE !'
                   CALL QUIT(
     &             'Error in reading NACSH_BOS(I) for .MJSELE !')
                 ENDIF
               END DO
               deallocate(SEL_NORB)
               GO TO 100
 42            CONTINUE
C&&& AOC: Average-of-Configuration
                  AOC = .TRUE.
               GOTO 100
 43            CONTINUE
C&&& FOCC: Fractional Occupation
                  AOC = .FALSE.
               GOTO 100
 44            CONTINUE
C&&& LSHIFT: Level shift
                  READ(LUCMD,*) DLSHIF
                  DOLEVEL = .TRUE.
                  IF (OPEN_FAC .LT. 0) OPEN_FAC = -1.0D0
               GOTO 100
 45            CONTINUE
C&&& Specify the maximum number of preliminary BSS-SCF iterations
                  READ(LUCMD,*) MAXITR2
               GOTO 100
 46            CONTINUE
CMI&&& Specify the ERGCNV2
C&&& ERGCNV2: Convergence on energy difference in preliminary BSS-SCF
                  READ(LUCMD,'(A80)') LINE
                  CALL UPCASE(LINE)
                  READ(LINE,*,END=403,ERR=403) SCFCNV2(1),SCFCNV2(2)
                  GOTO 404
  403             READ(LINE,*) SCFCNV2(1)
                  SCFCNV2(2) = SCFCNV2(1)
  404             CONTINUE
                  EVCCNV2 = .FALSE.
                  ERGCNV2 = .TRUE.
                  FCKCNV2 = .FALSE.
               GOTO 100
 47            CONTINUE
CMI&& Specify the EVCCNV2
C&&& EVCCNV: Convergence on approximate gradient (error vector)
                  READ(LUCMD,'(A80)') LINE
                  CALL UPCASE(LINE)
                  READ(LINE,*,END=503,ERR=503) SCFCNV2(1),SCFCNV2(2)
                  GOTO 504
  503             READ(LINE,*) SCFCNV2(1)
                  SCFCNV2(2) = SCFCNV2(1)
  504             CONTINUE
                  EVCCNV2 = .TRUE.
                  ERGCNV2 = .FALSE.
                  FCKCNV2 = .FALSE.
               GOTO 100
 48            CONTINUE
CMI&& Specify the FCKCNV2
C&&& FCKCNV2: Convergence on total Fock matrix
                  READ(LUCMD,'(A80)') LINE
                  CALL UPCASE(LINE)
                  READ(LINE,*,END=603,ERR=603) SCFCNV2(1),SCFCNV2(2)
                  GOTO 604
  603             READ(LINE,*) SCFCNV2(1)
                  SCFCNV2(2) = SCFCNV2(1)
  604             CONTINUE
                  EVCCNV2 = .FALSE.
                  ERGCNV2 = .FALSE.
                  FCKCNV2 = .TRUE.
               GOTO 100
 49            CONTINUE
C&&&  .MOFREEZE
               DOMOFREEZE = .TRUE.
               WRITE (LUPRI,*) 'WARNING: .MOFREEZE is experimental!'
               DO I = 1,NFSYM
                  READ(LUCMD,'(A72)') CTEMP
                  NMOFREEZE(I) = -1
                  CALL NUMLST(CTEMP,IMOFREEZE(1+(I-1)*(NMOFREEZE(1))),
     &                 MXMOFREEZE,1,MXMOFREEZE,
     &                 I,NMOFREEZE(I))
               END DO
               WRITE (LUPRI,*)'Freezing orbitals in fermion symmetry 1'
               WRITE (LUPRI,*) (IMOFREEZE(I),I=1,NMOFREEZE(1))
               IF (NFSYM.GT.1) THEN
                WRITE (LUPRI,*)'Freezing orbitals in fermion symmetry 2'
                WRITE (LUPRI,*) (IMOFREEZE(I),I=NMOFREEZE(1)+1,
     &               NMOFREEZE(1)+NMOFREEZE(2))
               ENDIF
               GOTO 100
 50            CONTINUE
C&&& OLEVEL: Open-shell level shifts
                  IF(NOPEN.EQ.0) THEN
                    CALL QUIT(' .OLEVEL: No open shells defined !')
                  ELSE
                    READ(LUCMD,*) (OLEV(IOPEN),IOPEN=1,NOPEN)
                    DOLEVEL = .TRUE.
                    IF (OPEN_FAC .LT. 0) OPEN_FAC = -1.0D0
                  ENDIF
               GOTO 100
 51            CONTINUE
C&&& ATOMST: Atomic start
               ATOMST=.TRUE.
               BARNUC=.FALSE.
               TRIVEC=.FALSE.
               TRIFCK=.FALSE.
               DO IATOM = 1,NONTYP
                 READ(LUCMD,*) ATOMFIL(IATOM),NVECATOM(IATOM)
                 DO ISHELL = 1,NVECATOM(IATOM)
                   READ(LUCMD,'(A72)') VECATOM(ISHELL,IATOM)
                   READ(LUCMD,*) OCCATOM(ISHELL,IATOM)
                 ENDDO
               ENDDO
               GOTO 100
 52            CONTINUE
C&&& SCFPOP: population analysis in each SCF iteration
                  SCFPOP = .TRUE.
               GOTO 100
 53            CONTINUE
C&&& OPENFA: extra factor on active-active correction block to open shell Fock matrix
                  READ(LUCMD,*) OPEN_FAC
               GOTO 100
 54            CONTINUE
C&&& BNC_FO: force bare nucleus correction (usefull for charged systems)
C    PBS+MI/Aug2014
               FORCE_BNC = .true.
               GOTO 100
 55            CONTINUE
!&&& FOMOUT: print out the Fock MO matrices for diagonalization for each
!symmetry into own formatted file (MI/Aug2014)
               WRITE_FMO_MATRIX=.TRUE.
               GOTO 100
 56            CONTINUE
C&&& AD HOC: extended Huckel based on atomic fragments
               IF(NBSYM.GT.1) THEN
                 CALL QUIT('Extended Huckel start '//
     &          'does not work with symmetry yet.')
               ENDIF
               ATHUCK = .TRUE.
               BNCRON = .FALSE.
               BARNUC=.FALSE.
               TRIVEC=.FALSE.
               TRIFCK=.FALSE.
               DHF_INTTYP='Atom. Huckel'
               IF(NONTYP.GT.MAXHOC) THEN
                 WRITE(LUPRI,'(A,I5)')
     &           '* Number of atomic types specified  :',NONTYP,
     &           '* Current maximum no. of types      :',MAXHOC
                 CALL QUIT('*AD HOC input: Too many atom types !')
               ENDIF
               DO ITYP = 1,NONTYP
                 READ(LUCMD,'(A6)')  HOCFIL(ITYP)
                 READ(LUCMD,'(A72)') VECHOC(ITYP)
               ENDDO
               GOTO 100
 57            CONTINUE
C&&& HUCPAR: modify Huckel parameter
               READ(LUCMD,*) HUCPAR
               GOTO 100
 58            CONTINUE
C&&& KPSELE: Selection based on Kappa symmetry
C    NSUBBL = # of subblocks for KP sym.
C    copied from MJSELE
               DOBOSSEL = .TRUE.
               IBOSSEL  = -1
C initla setup               
               IF(DOLVC.EQV..TRUE. .OR. IBTAND(INTGEN/4,1).EQ.1)THEN !4c
                 NSUBBL = (NHTYP-2)*2+1
               ELSE !2c
                 NSUBBL = (NHTYP-1)*2+1 
               ENDIF   
               KPBUF = NSUBBL/2+1
               IF(mod(KPBUF, 2).eq.0)THEN ! Max angular momentum (l) is s,d,g.. 
                 KPMAX(1) = KPBUF-1
                 KPMAX(2) = KPBUF
               ELSE ! Max l is p,f,h..
                 KPMAX(1) = KPBUF
                 KPMAX(2) = KPBUF-1      
               ENDIF
               IF(NFSYM.EQ.1)KPMAX(1) = MAX(KPMAX(1),KPMAX(2))
C readling input             
               READ(LUCMD,*,IOSTAT=IOS) NKP
               IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*) 'Error in reading NKP for .KPSELE !'
                 CALL QUIT('Error in reading NKP for .KPSELE !')
               ENDIF
               allocate(KP_NORB(NKP))
               allocate(SEL_NORB(NKP,0:NOPEN))
               READ(LUCMD,*,IOSTAT=IOS) (KP_NORB(I),   I=1,NKP)
               READ(LUCMD,*,IOSTAT=IOS) (SEL_NORB(I,0),I=1,NKP)
               DO I =1,NKP
                 IF(mod(SEL_NORB(I,0),2).eq.1)THEN
                   WRITE(LUPRI,*) 'Error in reading .KPSELE !'
                   CALL QUIT
     &             ('Number of electrons (closed shell) must be even!')                    
                 ELSE
                   KP = KP_NORB(I)
                   IF(NFSYM.EQ.1)THEN
                     J = ABS(KP)*2 + (KP/ABS(KP)-1)/2 !(KP,J)=(1,2)(-2,3)(2,4)etc
                     NISH_BOS(J) = SEL_NORB(I,0)/2
                   ELSE
                     IF(KP.LT.0)THEN
                       IF(mod(-KP,2).EQ.1)THEN
                         NISH_BOS(ABS(KP)         ) = SEL_NORB(I,0)/2
                       ELSE
                         NISH_BOS(ABS(KP)+KPMAX(1)) = SEL_NORB(I,0)/2        
                       ENDIF   
                     ELSE
                       IF(mod(KP,2).EQ.1)THEN
                         NISH_BOS(ABS(KP)+KPMAX(1)) = SEL_NORB(I,0)/2        
                       ELSE
                         NISH_BOS(ABS(KP)         ) = SEL_NORB(I,0)/2
                       ENDIF   
                     ENDIF 
                   ENDIF   
                 ENDIF
               ENDDO  
               IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*)
     &           'Error in reading NISH_BOS(I) for .KPSELE !'
                 CALL QUIT('Error in reading NISH_BOS(I) for .KPSELE !')
               ENDIF
               DO IOPEN = 1,NOPEN
                 READ(LUCMD,*,IOSTAT=IOS)
     &                (SEL_NORB(I,IOPEN),I=1,NKP)
                 DO I =1,NKP
                   IF(mod(SEL_NORB(I,IOPEN),2).eq.1)THEN
                     WRITE(LUPRI,*) 'Error in reading .KPSELE !'
                   CALL QUIT
     &             ('Number of orbitals (open shell) must be even!')                    
                   ELSE
                     KP = KP_NORB(I)  
                     IF(NFSYM.EQ.1)THEN
                       J = ABS(KP)*2 + (KP/ABS(KP)-1)/2 
                       NACSH_BOS(J,IOPEN) = SEL_NORB(I,IOPEN)/2
                     ELSE
                       IF(KP.LT.0)THEN
                         IF(mod(-KP,2).EQ.1)THEN
                           NACSH_BOS(ABS(KP)         ,IOPEN)
     &                     = SEL_NORB(I,IOPEN)/2
                         ELSE
                           NACSH_BOS(ABS(KP)+KPMAX(1),IOPEN)
     &                     = SEL_NORB(I,IOPEN)/2        
                         ENDIF   
                       ELSE
                         IF(mod(KP,2).EQ.1)THEN
                           NACSH_BOS(ABS(KP)+KPMAX(1),IOPEN)
     &                     = SEL_NORB(I,IOPEN)/2        
                         ELSE
                           NACSH_BOS(ABS(KP)         ,IOPEN)
     &                     = SEL_NORB(I,IOPEN)/2
                         ENDIF   
                       ENDIF 
                     ENDIF
                   ENDIF
                 ENDDO  
                 IF (IOS.NE.0) THEN
                   WRITE(LUPRI,*)
     &             'Error in reading NACSH_BOS(I) for .KPSELE !'
                   CALL QUIT(
     &             'Error in reading NACSH_BOS(I) for .KPSELE !')
                 ENDIF
               END DO
               deallocate(KP_NORB)               
               deallocate(SEL_NORB)               
               GO TO 100
 59            CONTINUE
C&&& SELECT: Generalized Selection based on Kappa symmetry
C    NSUBBL = # of subblocks.
C    copied from MJSELE
               DOBOSSEL = .TRUE.
               IBOSSEL  = -1
               READ(LUCMD,*,IOSTAT=IOS) NSUBBL
               IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*) 'Write # of subblocks in first line!'
                 CALL QUIT('Error in reading NSUBBL for .SELECT !')
               ENDIF
               allocate(SEL_NORB(NSUBBL,0:NOPEN))
               READ(LUCMD,*,IOSTAT=IOS) (SEL_NORB(I,0),I=1,NSUBBL)
               DO I =1,NSUBBL
                 IF(mod(SEL_NORB(I,0),2).eq.1)THEN
                   WRITE(LUPRI,*) 'Error in reading .SELECT !'
                   CALL QUIT('Number of electrons must be even!')                    
                 ELSE   
                   NISH_BOS(I) = SEL_NORB(I,0)/2               
                 ENDIF
               ENDDO  
               IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*)
     &           'Error in reading NISH_BOS(I) for .SELECT !'
                 CALL QUIT('Error in reading NISH_BOS(I) for .SELECT !')
               ENDIF
               DO IOPEN = 1,NOPEN
                 READ(LUCMD,*,IOSTAT=IOS)
     &                (SEL_NORB(I,IOPEN),I=1,NSUBBL)
                 DO I =1,NSUBBL
                   IF(mod(SEL_NORB(I,IOPEN),2).eq.1)THEN
                     WRITE(LUPRI,*) 'Error in reading .SELECT !'
                     CALL QUIT('Number of electrons must be even!')   
                   ELSE   
                     NACSH_BOS(I,IOPEN) = SEL_NORB(I,IOPEN)/2
                   ENDIF
                 ENDDO  
                 IF (IOS.NE.0) THEN
                   WRITE(LUPRI,*)
     &             'Error in reading NACSH_BOS(I) for .SELECT !'
                   CALL QUIT(
     &             'Error in reading NACSH_BOS(I) for .SELECT !')
                 ENDIF
               END DO
               deallocate(SEL_NORB)               
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/4A/)') ' Prompt "',WORD,
     *            '" not recognized for ',WORD1
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt under *SCF.')
            END IF
      END IF
  300 CONTINUE

C
C     Determine occupation
C     ====================
C
C     If there only is one symmetry there is no need to change
C     occupation during SCF cycles.
C
      IF (NFSYM .EQ. 1) AUTOCC = .FALSE.
C
      IF (INIOCC) THEN
C
         CALL RMOLCHR(ICHRG)
         WRITE(LUPRI,'(/A,3(/A,I3))')
     &      ' * Initial Automatic occupation based on:',
     &      '   Total charge of atoms = ',ICHRG,
     &      '   Charge of molecule    = ',KCHARG,
     &      '   i.e. no. of electrons = ',ICHRG-KCHARG
C
         NASHT_DHF  = 0
         NAELEC_DHF = 0
         IF (NFSYM .EQ. 1) THEN
C
C           The easy case:
C              all electrons go into the only symmetry!
C
            INIOCC = .FALSE.
            NELEC_DHF(1) = ICHRG - KCHARG
            IF (MOD(NELEC_DHF(1),2).EQ.1) THEN
               NOPEN = 1
               NACSH(1,1) = 1
               DF(1) = DP5
               NELEC_DHF(1) = NELEC_DHF(1) - 1
               WRITE(LUPRI,'(/,6X,A)')
     &            'WARNING! INIOCC assumes 1 electron in 1 open shell!'
               NASHT_DHF  = 1
            END IF
            NELECT_DHF = NELEC_DHF(1)
C
         ELSE
C
C           The difficult case:
C              The occupation is determined in the DHF routines.
C
            DF(0) = D1
            DA(0) = D1
            DALPHA(0) = D0
            NELECT_DHF = ICHRG - KCHARG
            IF (MOD(NELECT_DHF,2).EQ.1)  THEN
               NOPEN = 1
C
               DF(1) = DP5
               DA(1) = D0
               DALPHA(1) = D2
C
               NELECT_DHF = NELECT_DHF - 1
               NASHT_DHF  = 1
               WRITE(LUPRI,'(/,6X,A)')
     &            'WARNING! INIOCC assumes 1 electron in 1 open shell!'
            END IF
C
C           Since guessing the occupation is difficult, we should allow the programme to change
C           occupation during the SCF.
C
            AUTOCC = .TRUE.
         END IF
         NAELEC_DHF = NASHT_DHF
      END IF
      IF (.NOT. INIOCC) THEN
         DF(0)     = D1
         DA(0)     = D1
         DALPHA(0) = D0
         DTEMP     = D0
         NASHMFT   = 0
         DO IFRP = 1,NFSYM
            DO IOPEN = 1,NOPEN
               NASH_DHF(IFRP) = NASH_DHF(IFRP) + NACSH(IFRP,IOPEN)
               NACSHMF(IFRP,IOPEN) = NACSH(IFRP,IOPEN)
               NASHMFT = NASHMFT + NACSHMF(IFRP,IOPEN)
               DTEMP          = DTEMP +
     &                          real(NACSH(IFRP,IOPEN),8)*D2*DF(IOPEN)
            END DO
            IF (NELEC_DHF(IFRP).LT.0) THEN
               write(lupri,*) '*SCF error sym',ifrp,
     &             ' : negative number of '//
     &             'closed shell electrons :',nelec_dhf(ifrp)
               CALL QUIT('*SCF error: negative number of '//
     &             'closed shell electrons!')
            END IF
            IF (MOD(NELEC_DHF(IFRP),2).NE.0) THEN
               write(lupri,*) '*SCF error sym',ifrp,
     &             ' : odd even number of '//
     &             'closed shell electrons :',nelec_dhf(ifrp)
               CALL QUIT('*SCF error: odd even number of '//
     &             'closed shell electrons!')
            END IF
            NISH_DHF(IFRP)  = NELEC_DHF(IFRP)/2
         END DO
         NAELEC_DHF = NINT(DTEMP)
         DO IOPEN = 1,NOPEN
C
C           We do allow fractional occupation of 1
C           This corresponds to a inactive electron
C           and gives DALPHA = 0
C
            IF (ABS(DF(IOPEN)-D1).LE.DTHRS) THEN
               WRITE(LUPRI,'(1X,A,I2,A,F6.4)')
     &            'WARNING: Open shell no. ',IOPEN,
     &            ' has frac. occ. of ',DF(IOPEN)
               WRITE(LUPRI,'(1X,2A)')
     &            'This corresponds to an inactive shell. So ALPHA and',
     &            ' A is reset.'
               DA(IOPEN)     = D1
               DF(IOPEN)     = D1
               DALPHA(IOPEN) = D0
            ELSE IF (ABS(DF(IOPEN)) .LE. DTHRS) THEN
               WRITE(LUPRI,'(1X,A,I2,A,F6.4)')
     &            'WARNING: Open shell no. ',IOPEN,
     &            ' has frac. occ. of ',DF(IOPEN)
               WRITE(LUPRI,'(1X,2A)')
     &            'This corresponds to an secondary shell. So ALPHA',
     &            ' and A is reset.'
               DA(IOPEN)     = D0
               DF(IOPEN)     = D0
               DALPHA(IOPEN) = D0
            ELSE
               DM        = 2*(NACSH(1,IOPEN)+NACSH(2,IOPEN))
               DN        = DM*DF(IOPEN)
               DA(IOPEN) = (DM*(DN-D1))/(DN*(DM-D1))
               IF (ABS(DA(IOPEN)).LE.DTHRS) DA(IOPEN) = d0
               DALPHA(IOPEN) = (D1-DA(IOPEN))/(D1-DF(IOPEN))
!              print *,'dm,dn,df(iopen),da(iopen),dalpha(iopen)',
!    &                  dm,dn,df(iopen),da(iopen),dalpha(iopen)
            END IF
         ENDDO
         NASHT_DHF  = NASH_DHF(1)  + NASH_DHF(2)
         NFROT_DHF  = NFRO_DHF(1)  + NFRO_DHF(2)
         NELECT_DHF = NELEC_DHF(1) + NELEC_DHF(2)
      END IF
      IF(NASHT_DHF.EQ.0) THEN
        AOC     = .FALSE.
      ENDIF
      IF(PROJEC) THEN
        WRITE(LUPRI,'(A/A,1P,E13.4)')
     &      ' * Variational space reduced by projection'//
     &          ' using a set of fragment orbitals.',
     &      '   - Smallest norm accepted:',PRJTHR
        IF(PROOWN) THEN
          WRITE(LUPRI,'(3X,A)')
     &     '- Fragment orbitals are obtained from local bases !',
     &     '  The indexing of fragment spinor sets is assumed ',
     &     '  to follow list of symmetry independent nuclei.'
          DO J = 1,NFRAG
            WRITE(LUPRI,'(3X,A,I3,A,A4)') '- Fragment spinor set ',J,
     &             ' --> ',NAMN(J)
            WRITE(LUPRI,'(5X,A,A6)') ' - read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(7X,A,A3)')
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(7X,A,A3,A,A72)')
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO J = 1,NFRAG
            WRITE(LUPRI,'(3X,A,I3)') '- Fragment spinor set ',J
            WRITE(LUPRI,'(3X,A,A6)') '  read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(5X,A,A3)')
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(5X,A,A3,A,A72)')
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C     If no DIIS set MXDIIS to zero
C     =============================
C
      IF(.NOT.DODIIS) THEN
        MXDIIS = 0
      ENDIF
C
      IF(DOLVC .AND.
     &   .NOT.(ONECAP.AND.(INTV1C.EQ.2.OR.INTV1C.EQ.3)) ) THEN
        ISSINT = 0
      ENDIF

C   .... discard calculation of LS/SS/GT integrals
C          for the LevyLe/BSS/X2C Hamiltonians
      IF (LEVYLE.OR.
     &    ((BSS.or.x2c).AND..NOT.DO4C2C.AND..NOT.DO2C4C)
     &    .OR.TWOCOMP) THEN
        ISLINT = 0
        ISSINT = 0
        IGTINT = 0
      ENDIF
      IF (ZORA.AND..NOT.ZORA4) THEN
        ISSINT = 0
      ENDIF
      INTDEF = ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT
      IF(CNVINT(1).LT.DUMMY) ITRINT(1) = 1
      IF(CNVINT(2).LT.DUMMY) ITRINT(2) = 1
C
C     tsaue -
C     Temporary fix: load HF occupation numbers into dcborb.h
C
      CALL TDHFORB('NO SETDC2')
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      IF(dirac_cfg_dft_calculation) THEN
        WRITE(LUPRI,'(A)')
     &   ' *SCF: Set-up for Kohn-Sham calculation:'
      ELSE
        WRITE(LUPRI,'(A)')
     &   ' *SCF: Set-up for Hartree-Fock calculation:'
      ENDIF
      CALL PRSYMB(LUPRI,'=',75,0)
c
c     We disable evaluation of LS, SS, and GT integrals in
c     One-center app. type 4 /jkp
c
      IF (ONECAP.AND.(INTV1C.EQ.4).AND.INTDEF.NE.1) THEN
         INTDEF = 1
         WRITE(LUPRI,'(A)')
     &     ' * INTFLG reset to `1 0 0 0` since ONECAP type 4 has '//
     &     'been specified'
      END IF
C
C     Print info about #electrons etc.
C
      WRITE(LUPRI,'(A,I1)') ' * Number of fermion irreps: ',NFSYM
C
      IF (.NOT. INIOCC) CALL PROCC()
!     mi - initialize, otherwise runtime check failure
      NFROT_DHF  = NFRO_DHF(1)  + NFRO_DHF(2)
      IF(NFROT_DHF.GT.0) THEN
        WRITE(LUPRI,'(A)') ' * Number of frozen shells:'
        WRITE(LUPRI,'(3X,A,I2,A,3X,I3)')
     &   ('Fermion irrep ',I,':',NFRO_DHF(I),I=1,NFSYM)
      ENDIF
      IF (NELECT_DHF+NAELEC_DHF.EQ.0.AND.(.NOT.VACPOL)) THEN
        WRITE(LUPRI,'(2X,A)') 'Warning: zero electrons in SCF !'
      END IF

      IF(.NOT.dirac_cfg_scf_calculation) GOTO 999
C
      IF (BNCRON.OR.BNSPON) THEN
C     hjaaj nov 2004: bare nucleus correction is
C     designed for approx. neutral molecules,
C     and it will give very poor start guess for e.g.
C     a U(+90) calculation.
        IF (FORCE_BNC) THEN
           WRITE (LUPRI,'(A)')
     &     ' - INFO: bare nucleus correction used independently of '//
     &  'the system charge'
        ELSE
          CALL RMOLCHR(ICHRG)
          KCHARG = ICHRG - NELECT_DHF - NAELEC_DHF
          IF (ABS(KCHARG) .GT. 10) THEN
            BNCRON = .FALSE.
            BNCRON = .FALSE.
            WRITE (LUPRI,'(A)')
     &      ' - INFO: bare nucleus correction disabled because'//
     &      ' abs(molecular charge) .gt. 10'
          END IF
        END IF
        IF (BNCRON) THEN
          DHF_INTTYP='Scr. nuclei'
         WRITE (LUPRI,'(A)')
     &     ' * Bare nucleus screening correction used for start guess'
        END IF
        IF (BNSPON) THEN
          DHF_INTTYP='Atom. scrpot'
          WRITE (LUPRI,'(A)')
     &      ' * Sum of atomic potentials used for start guess'
        END IF
      ENDIF
C
      IF (DOBOSSEL) THEN
C     .BOSONS/.MJSELE were programmed by hjaaj Oct 2004
C     .KPSELE was added by ayaki Oct 2018.         
         NERR = 0
         IF (MAX_BOS_BL .NE. MAX_SUB_BL) THEN
            NERR = NERR + 1
            WRITE (LUPRI,*)
     &      'PROGRAMMING ERROR: MAX_BOS_BL .ne. MAX_SUB_BL'
         END IF
         IF (AUTOCC) THEN
            NERR = NERR + 1
            WRITE(LUPRI,*) 'INPUT ERROR! .AUTOCC combined with '//
     &           '.BOSONS/.MJSELE/.KPSELE does not work'
         END IF
         IF (INIOCC) THEN
            NERR = NERR + 1
            WRITE(LUPRI,*)
     &           'INPUT ERROR! .BOSONS/.MJSELE/.KPSELE requires'//
     &           'specification of .CLOSED/.OPEN SH'
         END IF
         IF (OVLSEL) THEN
            NERR = NERR + 1
            WRITE(LUPRI,*)
     &           'INPUT ERROR! .BOSONS/.MJSELE/.KPSELE combined with'
     &           //' .OVLSEL does not work'
         END IF
         IF (IBOSSEL .GT. 0) THEN
            IF (SPINFR) THEN
               WRITE(LUPRI,'(A)')
     *      ' * Vector selection based on spin-free boson symmetry'
            ELSE
               WRITE(LUPRI,'(A)')
     *      ' * Vector selection based on approximate boson symmetry'
            END IF
         ELSE IF (IBOSSEL .LT. 0) THEN
            IF(ATOMIC)THEN
              WRITE(LUPRI,'(A)')
     *        ' * Vector selection based on Kappa symmetry'
            ELSE  
              WRITE(LUPRI,'(A)')
     *        ' * Vector selection based on m_j symmetry'
            ENDIF  
            IF (.NOT. LINEAR) THEN
               NERR = NERR + 1
               WRITE(LUPRI,'(/A/)')
     &         'INPUT ERROR! Linear symmetry was not identified.'
            END IF
         END IF
         WRITE(LUPRI,'(A/4X,16I4)')
     &      '   - inactive shell:',(NISH_BOS(J)*2,J=1,NSUBBL)
         I1 = NISH_DHF(1) + NISH_DHF(2)
         I2 = ISUM(NSUBBL,NISH_BOS,1)
         IF (I1 .NE. I2) THEN
            NERR = NERR + 1
            WRITE(LUPRI,'(2(/A,I5,A))')
     &      'INPUT ERROR! This shell should have',I1,' orbitals!',
     &      '  ...while is actually has  ',I2,' orbitals - fix it!'
         END IF
         IF (NOPEN.GT.0) THEN
            WRITE(LUPRI,'(A)') '   - active shell(s):'
            DO IOPEN = 1,NOPEN
               WRITE(LUPRI,'(4X,16I4)')(NACSH_BOS(J,IOPEN)*2,J=1,NSUBBL)
               I1 = NACSH(1,IOPEN) + NACSH(2,IOPEN)
               I2 = ISUM(NSUBBL,NACSH_BOS(1,IOPEN),1)
               IF (I1 .NE. I2) THEN
                  NERR = NERR + 1
                  WRITE(LUPRI,'(/A,I3,A/)')
     &            'INPUT ERROR! This shell should have',I1,' orbitals!'
               END IF
            END DO
         END IF
C   INPUT      : Only Kappa
C   Computation: Kappa and MJ
C   So, we have to devide each array.
         IF(ATOMIC)THEN
           IF(NFSYM.EQ.2)THEN
             INVMAX   = 1
             IOFM = KPMAX(1)*(KPMAX(1)+1)/2
             IOFK = KPMAX(1)
           ELSE              
             INVMAX   = 2
             IOFM = 0
             IOFK = 0
           ENDIF
C          
           DO I =  NFSYM, 1, -1
             DO IKAP = KPMAX(I), 1, -1
               IF(NFSYM.EQ.2)THEN                 
                 KPADD = IKAP*(IKAP-1)/2 + 1
               ELSE
                 KPADD = IKAP*(IKAP-1)   + 1 + IKAP
               ENDIF
               DO INV=1,INVMAX
                 IOFKV = -INV+1+IOFK
                 IF(mod(NISH_BOS(IKAP*INVMAX+IOFKV),IKAP).eq.0)THEN
                   KCBUF = NISH_BOS(IKAP*INVMAX+IOFKV)/IKAP
                 ELSE
                   WRITE(LUPRI,'(2A,I2,A,3I3,A)')
     &              'Number of electrons in the closed shell',
     &              ' ABS(Kappa) =',IKAP,
     &              " should be",IKAP*2,IKAP*4,IKAP*6,
     &              " etc."
                   CALL QUIT('Number of elec. of closed shell is wrong') 
                 ENDIF
                 DO IK = IKAP, 1, -1  
                   NISH_BOS(KPADD+IK-1+IOFM) = KCBUF
                 ENDDO
C
                 DO IOPEN = 1, NOPEN
                   IF(mod(NACSH_BOS(IKAP*INVMAX+IOFKV,IOPEN),IKAP)
     &              .eq.0)THEN
                     KOBUF = NACSH_BOS(IKAP*INVMAX+IOFKV,IOPEN)/IKAP
                   ELSE
                     WRITE(LUPRI,'(A,I2,A,I1,A,3I3,A)')
     &              'Number of orbitals in ABS(Kappa) =',IKAP,
     &               " of open-shell No. ",IOPEN,
     &               " should be",IKAP*2,IKAP*4,IKAP*6,
     &               " etc."
                     CALL QUIT
     &                ('Number of orbitals of open-shell is wrong')   
                   ENDIF                  
                   DO IK = IKAP, 1, -1  
                     NACSH_BOS(KPADD+IK-1+IOFM,IOPEN) = KOBUF
                   ENDDO  
                 ENDDO
                 KPADD = KPADD-IKAP
               ENDDO  !INV            
             ENDDO    !IKAP
             IOFK = 0
             IOFM = 0
           ENDDO      !NFSYM
         ENDIF        ! IF(ATOMIC)
C  TODO HJAAJ: for NFSYM.eq.2
C  TODO HJAAJ: check NISH_DHF(1/2) against gerade/ungerade NISH_BOS etc.
         IF (NERR .GT. 0) THEN
            CALL QUIT('Input error for .BOSONSelection or .MJSELEction')
         END IF
      END IF  !(DOBOSSEL)
C
      WRITE(LUPRI,'(A,I3)') ' * General print level   : ',IPRSCF
      WRITE(LUPRI,'(/A)') ' ***** INITIAL TRIAL SCF FUNCTION *****'
      IF     (TRIFCK) THEN
        WRITE(LUPRI,'(A)')
     +     ' * Trial two-electron Fock matrix read from file DF2FCK'
      ELSE IF(TRIVEC) THEN
        WRITE(LUPRI,'(A)')
     +   ' * Trial vectors read from file DFCOEF'
      ELSE IF(ATOMST) THEN
        WRITE(LUPRI,'(A,I5//4X,A)') ' * Atomic start:',NONTYP,
     &   'atom   cffile     occ    orbitals'
        ICENT = 1
        DO I = 1,NONTYP
          WRITE(LUPRI,'(5X,A4,3X,A6,3X,F6.2,3X,A72)')
     &          NAMN(ICENT),ATOMFIL(I),OCCATOM(1,I),VECATOM(1,I)
          DO J = 2,NVECATOM(I)
            WRITE(LUPRI,'(5X,16X,F6.2,3X,A72)')
     &            OCCATOM(J,I),VECATOM(J,I)
          ENDDO
          ICENT = ICENT + NONT(I)
        ENDDO
        WRITE(LUPRI,'()')
      ELSE IF(ATHUCK) THEN
           WRITE(LUPRI,'(A)')
     &   ' * Trial vectors generated by the extended Huckel method'//
     &   ' using pre-calculated atomic orbitals'
           WRITE(LUPRI,'(A,F7.3)') '   Huckel parameter: ',HUCPAR
      ELSE IF(SIRIFC) THEN
        WRITE(LUPRI,'(A)')
     +   ' * Trial vectors read from Dalton file SIRIFC'
      ELSE IF(BARNUC) THEN
        WRITE(LUPRI,'(A)')
     +   ' * Trial vectors generated by the bare '//
     +   'nucleus approximation'
      ELSE IF(DOHUCKEL) THEN
        IF (EWMO) THEN
           WRITE(LUPRI,'(A)')
     +   ' * Trial vectors generated by the Energy Weighted'//
     +   ' Maximum Overlap start guess (EWMO)'
        ELSE
           WRITE(LUPRI,'(A)')
     +   ' * Trial vectors generated by the extended Huckel '//
     +   'start guess (EHT)'
        ENDIF
      ELSE
        CALL QUIT('ERROR, initial trial SCF function not defined!')
      ENDIF
      IF (abs(OPEN_FAC) .NE. D1) THEN
         WRITE (LUPRI,'(A,F12.6)')
     &   ' * Scaling of active-active block correction to '//
     &   'open shell Fock operator', abs(OPEN_FAC)
         IF (OPEN_FAC .LT. 0.0D0) THEN
            WRITE (LUPRI,'(A)')
     &   '   to improve convergence (default value).'
         END IF
      END IF

      WRITE(LUPRI,'(/A)') ' ***** SCF CONVERGENCE CRITERIA *****'
      IF(ERGCNV) THEN
        WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &    ' * Convergence on total energy.',
     &    '   Desired convergence:',SCFCNV(1),
     &    '   Allowed convergence:',SCFCNV(2)
      ELSEIF(FCKCNV) THEN
        WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &    ' * Convergence on total Fock matrix.',
     &    '   Desired convergence:',SCFCNV(1),
     &    '   Allowed convergence:',SCFCNV(2)
      ELSEIF(EVCCNV) THEN
        WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &    ' * Convergence on norm of error vector (gradient).',
     &    '   Desired convergence:',SCFCNV(1),
     &    '   Allowed convergence:',SCFCNV(2)
      ENDIF
      WRITE(LUPRI,'(A)') ' ***** CONVERGENCE CONTROL *****'
      IF(DODSCF) THEN
        WRITE(LUPRI,'(A)') ' * Fock matrix constructed using'//
     &       ' differential density matrix'
        IF(FIXDIF) THEN
          WRITE(LUPRI,'(4X,A)')
     &       'with fixed parameter 1.0'
        ELSE
          WRITE(LUPRI,'(4X,A)')
     &       'with optimal parameter.'
        ENDIF
      ENDIF
      IF(DODIIS) THEN
         IF (DIISAO) THEN
           WRITE(LUPRI,'(A)') ' * DIIS (in AO basis)'
         ELSEIF(DIISMO) THEN
           WRITE(LUPRI,'(A)') ' * DIIS (in MO basis)'
         ENDIF
         WRITE(LUPRI,'(A,1P,D8.2)')
     &      ' * DIIS will be activated when convergence reaches : ',
     &      DIISTH
         WRITE(LUPRI,'(3X,A,I5)') '- Maximum size of B-matrix:',
     &      MXDIIS
         IF (NOSWIT) THEN
           WRITE(LUPRI,'(A)') ' * DIIS is never switched off'
         ENDIF
         IF(DODAMP) THEN
           WRITE(LUPRI,'(A/A,F5.3)')
     &     ' * Damping of Fock matrix when DIIS is not activated. ',
     &     '   Weight of old matrix    : ',DAMPFC
         END IF
      ELSE IF(DODAMP) THEN
         WRITE(LUPRI,'(A,F5.3)')
     &     ' * Damping of Fock matrix. Weight of old matrix    : ',
     &     DAMPFC
      END IF
      IF(DOLEVEL) THEN
        IF(DLSHIF.NE.D0) THEN
          WRITE(LUPRI,'(A,F8.3)')
     &     ' * Level shift for virtuals activated : ',DLSHIF
        ENDIF
        DO IOPEN = 1,NOPEN
          IF(OLEV(IOPEN).NE.D0) THEN
            WRITE(LUPRI,'(A,I2,A,F8.3)')
     &      ' * Level shift for open-shell ',
     &      IOPEN,' activated : ',OLEV(IOPEN)
          ENDIF
        ENDDO
      ENDIF
      WRITE(LUPRI,'(A,I4)')
     &   ' * Maximum number of SCF iterations  : ',MAXITR
C
C     Quadratic convergent DHF?
C
      IF (NOQCDHF) THEN
         WRITE(LUPRI,'(A)')
     &   ' * No quadratic convergent Hartree-Fock'
      ELSE
         WRITE(LUPRI,'(A,2(/A,I3))')
     &   ' * Quadratic convergent Hartree-Fock',
     &   '   - Maximum number of macro iterations:',MXMACRO,
     &   '   - Maximum number of micro iterations:',MXMICRO
      END IF
#ifdef HAS_PCMSOLVER
      if (.not. noqcdhf .and. dirac_cfg_pcm) then
         call quit('Quadratic convergent Hartree-Fock not'//
     &   'available with PCM')
      endif              
#endif
      IF (DYNSEL .AND. .NOT. OVLSEL) OVLSEL = .TRUE.
      IF (OVLSEL) THEN
         IF (AUTOCC)
     &      CALL QUIT('ERROR! .AUTOCC combined with overlap '//
     &           'selection does not work')
         IF(.NOT. DYNSEL) THEN
            WRITE(LUPRI,'(A)') ' * Vector selection based on overlap'
         ELSE
            WRITE(LUPRI,'(A/3X,A)')
     &         ' * Vector selection based on overlap ',
     &            '(dynamic update of overlap-selection matrix)'
         END IF
      END IF
      IF (AUTOCC .AND. NASHT_DHF .GT. 1) THEN
         WRITE(LUPRI,'(A)')
     &      ' WARNING! .AUTOCC and NASHT > 1 does not work',
     &      '  ==> turning .AUTOCC off'
         AUTOCC = .FALSE.
      END IF
      IF (AUTOCC .AND. NOPEN .GT. 1) THEN
         WRITE(LUPRI,'(/A)')
     &      ' WARNING! .AUTOCC and NOPEN > 1 does not work',
     &      '  ==> turning .AUTOCC off'
         AUTOCC = .FALSE.
      END IF
      IF (AUTOCC .AND. NFSYM .GT. 1) THEN
         WRITE(LUPRI,'(A)')
     &      ' * DHF occupation is allowed to change during SCF cycles.'
      END IF
      IF (AUTOCC .AND. AOC) THEN
         WRITE(LUPRI,'(A)')
     &      ' WARNING! .AUTOCC and average-of-configuration active.',
     &      '  ==> you may get what you asked for, but no warranty!!!'
      END IF
      WRITE(LUPRI,'(A)')
     +    ' * Contributions from 2-electron integrals to Fock matrix:'
      IF(LBIT(INTDEF,1)) THEN
        WRITE(LUPRI,'(3X,A)') 'LL-integrals.'
      ENDIF
      IF(LBIT(INTDEF,2)) THEN
        IF(CNVINT(1).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,1P,D8.1)')
     +    'SL-integrals below SCF convergence ',CNVINT(1)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    'SL-integrals from iteration ',ITRINT(1)
        ENDIF
      ENDIF
      IF(LBIT(INTDEF,3)) THEN
        IF(CNVINT(2).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,1P,D8.1)')
     +    'SS-integrals below SCF convergence ',CNVINT(2)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    'SS-integrals from iteration ',ITRINT(2)
        ENDIF
      ENDIF
      IF(LBIT(INTDEF,4)) THEN
        IF(CNVINT(2).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,1P,D8.1)')
     +    'Gaunt integrals  below SCF convergence ',CNVINT(2)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    'Gaunt integrals from iteration ',ITRINT(2)
        ENDIF
      ENDIF
      IF (INTFLG_CHANGE) THEN
        WRITE(LUPRI,"(4X,A)")
     &  "---> accepted user's setting through .INTFLG"
      ELSE
        WRITE(LUPRI,"(4X,A)")
     &  "---> this is default setting from Hamiltonian input"
      ENDIF
C
C     Rotations in QC DHF
C
      IF(DHF_SKIPEE) THEN
         WRITE(LUPRI,'(A)')
     &      ' * NB!!! No e-e rotations in 2nd order optimization'
      ENDIF
      IF(DHF_SKIPEP) THEN
         WRITE(LUPRI,'(A)')
     &      ' * NB!!! No e-p rotations in 2nd order optimization.'
      ENDIF
C
      IF ((ONECNV .NE. D0) .AND. (ONECNV .LT. SCFCNV(1))) THEN
         WRITE(LUPRI,'(2(/A,D12.6),/A)')
     &   ' ** ONECAP-WARNING : ONECNV = ',ONECNV,
     &   '                     SCFCNV(1) = ',SCFCNV(1),
     &   '                     ONECNV re-initialized to SCFCNV(1)'
         ONECNV = SCFCNV(1)
      END IF

CMI  ... print out info on restarting 4c DC-SCF from previous 2c DC-SCF
CMI =====================================================================
      IF (DO2C4C.AND.START2C) THEN
       WRITE(LUPRI,'(A)')
     &  ' ***** Restarting four-component DC-SCF'//
     &  ' from the previous two-component IOTC/BSS-SCF *****'
       WRITE(LUPRI,*) '...restart code, INI2C=',INI2C,
     &  ' This means restart from:'
       IF (INI2C.EQ.1) THEN
         write(lupri,*)
     &   ' Fock MO matrix (backtransformed FD2C_EE+H1_PP)'
       ELSE IF (INI2C.EQ.2) THEN
         write(lupri,*)
     & ' Fock MO matrix (backtransformed FD2C_EE)'
       ELSE IF (INI2C.EQ.3) THEN
         write(lupri,*)
     &  ' backtransformed (electronic) MOs - THE BEST CHOICE (default)'
       ELSE IF (INI2C.EQ.4) THEN
         write(lupri,*)
     &  ' untransformed 2-component (electronic) MOs'
       ELSE
        CALL QUIT('*SCF: wrong value of INI2C (not 1/2/3/4)!')
       ENDIF

       WRITE(LUPRI,'(A,I4)')
     & ' * Maximum number of preliminary BSS-SCF iterations  : ',
     &  MAXITR2

       WRITE(LUPRI,'(A)')
     & ' ***** CONVERGENCE CRITERIA ON PRELIMINARY BSS-SCF *****'
      IF(ERGCNV2) THEN
        WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &    ' * Convergence on total energy.',
     &    '   Desired convergence:',SCFCNV2(1),
     &    '   Allowed convergence:',SCFCNV2(2)
      ELSEIF(FCKCNV2) THEN
        WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &    ' * Convergence on total Fock matrix.',
     &    '   Desired convergence:',SCFCNV2(1),
     &    '   Allowed convergence:',SCFCNV2(2)
      ELSEIF(EVCCNV2) THEN
        WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &    ' * Convergence on norm of error vector (gradient).',
     &    '   Desired convergence:',SCFCNV2(1),
     &    '   Allowed convergence:',SCFCNV2(2)
      ENDIF
      ENDIF

!     ... descend to the two-component level
      IF (DO4C2C) THEN
!      reset the X2C/BSS module flag
       bss = .false.
       if(x2c)then
         x2c        = .false.
!        keep the x2cmod_x2c value in order to identify to which module
!        to switch in dirscf.F --> subroutine PSISCF
       end if
       WRITE(LUPRI,'(A)')
     &  ' * after the four-component SCF step switch to'//
     &  ' the two-component framwork (X2C/BSS)'
       IF (USEDF) THEN
         WRITE(LUPRI,'(2X,A)')
     &   '...using the Fock-Dirac matrix as the defining h1'//
     &   ' Hamiltonian matrix for the X2C transformation.'
       ELSE
         WRITE(LUPRI,'(2X,A)')
     &   '...using the Dirac bare nucleus, h_D(1), matrix as the'//
     &   ' defining h1 Hamiltonian matrix for the X2C transformation.'
       ENDIF
       IF (CONT2C) THEN
         WRITE(LUPRI,'(2X,A)')
     &   '...continue the two-component SCF procedure'//
     &   '.'
       ENDIF
      ENDIF

      IF (WRITE_FMO_MATRIX) THEN
        WRITE(LUPRI,'(2X,A)')
     &'* Writing out Fock MO matrix(es) into formatted file(s) in'//
     &' linear symmetry'
      ENDIF
C
      WRITE(LUPRI,'(A)') ' ***** OUTPUT CONTROL *****'
      IF (SCFPOP) WRITE(LUPRI,'(A)')
     +  ' * Mulliken population analysis each SCF iteration.'
      IF(IPREIG.EQ.0) WRITE(LUPRI,'(A)')
     +  ' * No eigenvalues written out.'
      IF(IPREIG.EQ.1) WRITE(LUPRI,'(A)')
     +  ' * Only electron eigenvalues written out.'
      IF(IPREIG.EQ.2) WRITE(LUPRI,'(A)')
     +  ' * Only positron eigenvalues written out.'
      IF(IPREIG.EQ.3) WRITE(LUPRI,'(A)')
     +  ' * Both electron and positron component eigenvalues '
     +  //'written out.'
C
  999 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck idnnuc */
      SUBROUTINE IDNNUC(IPRINT)
C***********************************************************************
C
C       Count number of nuclei with identical coordinates
C
C       T.Saue, March 24 1994, University of Odense
C       Last revision : Nov 14 1994 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "nuclei.h"
#include "symmet.h"
#include "ccom.h"
#include "shells.h"
C
      CALL QENTER('IDNNUC')
C
C     Identify nuclei with identical coordinates and give
C     pointer to the first in the series
C     ===================================================
C
      CALL IZERO(INCENT,NUCIND)
      CALL IZERO(INUNIQ,NUCIND)
      NCLOSE = 0
      NUCNET = 0
      DO 20 I = 1,NUCIND
      IF(INCENT(I).NE.0) GOTO 20
        INCENT(I) = I
        NUCNET = NUCNET + 1
        INUNIQ(NUCNET) = I
        DO 30 J = (I+1),NUCIND
        IF(INCENT(J).EQ.J) GOTO 30
          DIST = (CORD(1,I)-CORD(1,J))**2
     +          +(CORD(2,I)-CORD(2,J))**2
     +          +(CORD(3,I)-CORD(3,J))**2
          IF(DIST.LT.0.1) THEN
            INCENT(J) = I
            NCLOSE = NCLOSE + 1
          ENDIF
   30   CONTINUE
   20 CONTINUE
      IF(NUCNET.LT.NUCIND) THEN
        CALL HEADER('Output from IDNNUC',-1)
        WRITE(LUPRI,'(3X,A,I5/9X,A,I5,A)')
     +    '* The number of centers with unique coordinates: ',
     +    NUCNET,
     +    '(out of ',NUCIND,' centers in total)'
        WRITE(LUPRI,'(3X,A)') 'Pointer INCENT:'
        WRITE(LUPRI,'(6X,2I5)') (I,INCENT(I),I=1,NUCIND)
        WRITE(LUPRI,'(3X,A)') 'Pointer INUNIQ:'
        WRITE(LUPRI,'(6X,2I5)') (I,INUNIQ(I),I=1,NUCNET)
      ENDIF
C
      CALL QEXIT('IDNNUC')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fldinp */
      SUBROUTINE FLDINP(WORD,INPERR)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER WORD*7
C
      WRITE(LUPRI,'(A)') 'FLDINP: Sorry, not yet written !'
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiinp */
      SUBROUTINE PSIINP(WORK,LWORK)
C***********************************************************************
C
C <<< Wave function Input for DIRAC >>>
C
C
C       Written by Trond Saue, Tromsoe April 1996
C       Last revision: April 27 1996 - tsaue
C
C***********************************************************************
      use dirac_cfg
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 21, NTABLE = 26)
C
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7,CTEMP*90
      DIMENSION WORK(LWORK)
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbpsi.h"
#include "dcbham.h"
#include "dcbreo.h"
#include "dcbrot.h"
C
      DATA TABDIR /'*END OF','xxxxxxx','*MP2CAL','*LUCITA',
     &             '*MVOCAL','*RESOLV','*KRMCSC','*STEP C',
     &             '*OPTIMI','*MP2 NO','*KRCICA','*COSCI ',
     &             '*GASCIP','*HSFS  ','*KRCCCA','*ARDUCC',
     &             '*SCF   ','*KRCI  ','*QCORR ','*LAPLCE',
     &             '*EXACC '/ 

      DATA TABLE  /'xxxxxxx','.MP2   ','.LUCITA','.REORDE',
     &             '.POST S','.RELCCS','.DIRRCI','.MVO   ',
     &             '.RESOLV','.ORBROT','.LUCIAR','.KRMCSC',
     &             '.PHCOEF','.MP2 NO','.RELADC','.KR CI ',
     &             '.KR CC ','.COSCI ','.GASCIP','.HSFS  ',
     &             '.ARDUCC','.SCF   ','.KRCI  ','.POLPRP',
     &             '.EXACC ','.LAPLCE'/
C
      N_WF = 0
C
C     Initialize /CBLPSI/
C     ===================
C
      DOMP2   = .FALSE.
      DORES   = .FALSE.
      DOEXACC = .FALSE.
      DOCCM   = .FALSE.
      DOCIM   = .FALSE.
      DOLUCT  = .FALSE.
      DOKRMC  = .FALSE.
      DOLUCIAR= .FALSE.
      L1ORBM  = .FALSE.
      L2ORBM  = .FALSE.
      PHCOEF  = .FALSE.
      DOADC   = .FALSE.
      DOKRCI  = .FALSE.
      DOKRCC  = .FALSE.
      CBKRINI = .FALSE.
      DOCOSCI = .FALSE.
      DOGASCIP= .FALSE.
      DOHSFS  = .FALSE.
      DOPOLPRP= .FALSE.
      DOLAPLCE= .FALSE.
C
C *** Initialize /CBIREO/ in dcbreo.h
C     ===================
C
      LMOORD    = .FALSE.
      LMOORF    = .FALSE.
C
C     Read menu file
C     ==============
C     **** Find Wave function input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      CALL UPCASE(WORD)
      IF (WORD .EQ. '**WAVE ') THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         CALL QUIT(
     &   'End of file on DIRAC.INP, no **WAVE function input found')
  920 CONTINUE
         CALL QUIT(
     &   'Error reading DIRAC.INP, no **WAVE function input found')
  930 CONTINUE
      WORD1 = WORD
C
C     Process input for COMMON  /CBIPSI/
C     ==================================
C
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,
     &                109,110,111,112,113,114,115,116,
     &                117,118,119,120,121,122,123,124,
     &                125,126),I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
             CALL PRTAB(NDIR,TABDIR, WORD1//' directory keywords',LUPRI)
             CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
             GO TO 100
            END IF
            WRITE (LUPRI,'(/3A/)')
     &         ' Keyword "',WORD,'" not recognized in PSIINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in PSIINP.')

  101    continue
            go to 100

C&&&& MP2: Perform second-order Moller-Plesset calculation
  102    CONTINUE
            DOMP2 = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'MP2'
            GO TO 100
C&&&& LUCITA: Perform spinfree CI calculation
  103    CONTINUE
!           stefan - nov 2011: check for integer compatibility for lucita on x86_64 systems
            call test_lucita_wrk_space_offset(work,lwork)
            IF (.NOT.NOTRA) DOTRA = .TRUE.
            DOLUCT = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'LUCITA'
            GO TO 100
  104    CONTINUE
C&&& REORDER: Reorder MO orbitals (only on restart on coefficients)
C             First read the number of orbital reorderings
C             in each fermion corep.
C
            LMOORD = .TRUE.
            L1ORBM = .TRUE.
            IREORD(1) = 0
            IREORD(2) = 0
            DO I = 1,NFSYM
              READ(LUCMD,'(A90)') CTEMP
              IREORD(I) = -1
              CALL NUMLST(CTEMP,IMOORD(1,I),MXREORD,1,MXREORD,
     &                    I,IREORD(I))
            END DO
            GO TO 100
  105    CONTINUE
C&&& POST SCF REORDER:
C           Reorder finished MO orbitals
C           First read the number of orbital reorderings
C           in each fermion corep.
C
            L2ORBM = .TRUE.
            LMOORF = .TRUE.
            IREORF(1) = 0
            IREORF(2) = 0
            DO I = 1,NFSYM
              READ(LUCMD,'(A72)') CTEMP
              IREORF(I) = -1
              CALL NUMLST(CTEMP,IMOORF(1,I),MXREORD,1,MXREORD,
     &                    I,IREORF(I))
                  END DO
            GO TO 100
  106    CONTINUE
C&&&& RELCCSD: Perfrom CCSD, CCSD(T) or MP2 with relccsd
            IF (.NOT.NOTRA) DOTRA = .TRUE.
            DOCCM = .TRUE.
            dirac_cfg_relcc = .true.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'RELCCSD'
            GO TO 100
  107    CONTINUE
C&&&& DIRRCI: Perform CISD or COSCI
            IF (.NOT.NOTRA) DOTRA = .TRUE.
            DOCIM = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'DIRRCI'
            GO TO 100
  108    CONTINUE
C&&&& MVO: Calculate modified virtual orbitals
            DOMVO = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'MVO'
            GO TO 100
  109    CONTINUE
C&&&& RESOLVE: Resolve open-shell states
            DORES = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'RESOLVE'
            GO TO 100
  110    CONTINUE
C&&&& ORBROT: Rotate orbitals prior to SCF calculations
            ROTORB = .TRUE.
            L1ORBM = .TRUE.
            READ(LUCMD,*) ROTANG
            READ(LUCMD,'(A72)') (VECROT(I),I=1,NFSYM)
            GO TO 100
  111    CONTINUE
C&&&& LUCIAREL: Call LUCIAREL
            DOLUCIAR = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'LUCIAREL'
            GO TO 100
  112    CONTINUE
C&&&& KRMCSCF: KR-MCSCF
            DOKRMC = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'KR-MCSCF'
            GO TO 100
  113    CONTINUE
C&&&& PHCOEF: Phase adjustment of coefficients
            L2ORBM = .TRUE.
            PHCOEF = .TRUE.
            GO TO 100
  114    CONTINUE
C&&&& MP2 NO: Calculate MP2 natural orbitals and their occupations
            DOMP2NO = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'MP2 NO'
            GO TO 100
  115    CONTINUE
C&&&& RELADC: Calculate relativistic propagators with ADC
            IF (.NOT.NOTRA) DOTRA = .TRUE.
            DOADC   = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'RELADC'
            GO TO 100
C&&&& KR CI :  perform (large-scale) KR-CI calculation
  116    CONTINUE
            DOKRCI  = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'KR-CI'
            GO TO 100
C&&&& KR CC : perform (large-scale) KR-CC calculation
  117    CONTINUE
            DOKRCC  = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'KR-CC'
            GO TO 100
C&&& COSCI: Do a COSCI calculation
  118    CONTINUE
            DOCOSCI = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'COSCI'
            GO TO 100
C&&& GASCIP: Do a GASCIP calculation (new code)
  119    CONTINUE
            DOGASCIP = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'GASCIP'
            GO TO 100
C&&& HSFS : Do High Sect FSCC calculations (new code)
  120    CONTINUE
            DOHSFS = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'HS-FSCC'
            GO TO 100
C&&&& ARDUCC: Perform MRCC calculation
  121    CONTINUE
            IF (.NOT.NOTRA) DOTRA = .TRUE.
            DOARDUC = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'AR-DU-CCA'
            GO TO 100
C&&&& SCF
  122    CONTINUE
            N_WF = N_WF + 1
            dirac_cfg_scf_calculation = .true.
            if (.not. dirac_cfg_dft_calculation) then
               dirac_cfg_hf_calculation = .true.
               DIRAC_WF(N_WF) = 'HF'
            else
               DIRAC_WF(N_WF) = 'DFT'
            end if
            go to 100
C&&&& KRCI :  alternative spelling of KR CI
  123    CONTINUE
            GO TO 116
C&&&& POLPRP: Calculate polarization propagator for excitations
  124    CONTINUE
            IF (.NOT.NOTRA) DOTRA = .TRUE.
            DOPOLPRP   = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'POLPRP'
            GO TO 100

C&&&& EXACC: Perfrom coupled cluster calculations with exacc
  125    CONTINUE
            DOEXACC = .TRUE.
            dirac_cfg_exacc = .true.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'EXA_CC'
            GO TO 100

C&&&& LAPLCE: Laplace transformation test program
  126    CONTINUE
            DOLAPLCE = .TRUE.
            N_WF = N_WF + 1
            DIRAC_WF(N_WF) = 'LAPLCE'
            GO TO 100

      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 180
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in PSIINP.')
      END IF
  180 CONTINUE
C
C     Print section
C     =============
C
      CALL TITLER('Wave function module','*',126)
      IF (N_WF.GT.0) THEN
         WRITE(LUPRI,'(A)')
     &      ' Wave function types requested (in input order):'
         DO I_WF = 1, N_WF
            WRITE(LUPRI,'(5X,A)') DIRAC_WF(I_WF)
         END DO

         WRITE(LUPRI,'(/A)')
     &      ' Wave function jobs in execution order (expanded):'

      ENDIF

      if (dirac_cfg_hf_calculation)  WRITE(LUPRI,'(A)')
     &   ' * Hartree-Fock calculation'
      if (dirac_cfg_dft_calculation) WRITE(LUPRI,'(A)')
     &   ' * Density functional calculation (Kohn-Sham method)'

      IF(DORES)    WRITE(LUPRI,'(A)')
     &   ' * Followed by resolution of open-shell states'
      IF(DOMP2)    WRITE(LUPRI,'(A)') ' * MP2 calculation'
      IF(DOMVO)    WRITE(LUPRI,'(A)')
     &   ' * Generate modified virtual orbitals'
      IF(DOMP2NO)  WRITE(LUPRI,'(A)')
     &   ' * Generate MP2 natural orbitals'
      IF(DOKRMC)   WRITE(LUPRI,'(A)')
     &   ' * Kramers restricted MCSCF calculation'
      IF(DOLUCIAR) THEN
C hj aug03: TODO: idea was/is to maybe follow MCSCF by big CI
         WRITE(LUPRI,'(A)') ' * Run LUCIAREL CI code'
         WRITE(LUPRI,'(A)') '   ERROR: not implemented yet'
         CALL QUIT('.LUCIAREL as stand-alone not implemented yet')
      END IF
      IF(DOCCM)    WRITE(LUPRI,'(A)') ' * Run RELCCSD code'
#ifdef MOD_HSFS
      IF (DOHSFS)  WRITE(LUPRI,'(A)') ' * Run Higher Sectors FSCC code'
#endif
      IF(DOCOSCI)  WRITE(LUPRI,'(A)') ' * Run COSCI CI code'
      IF(DOGASCIP) WRITE(LUPRI,'(A)') ' * Run GASCIP CI code'
      IF(DOCIM)    WRITE(LUPRI,'(A)') ' * Run DIRRCI CI code'
      IF(DOLUCT)   WRITE(LUPRI,'(A)') ' * Run LUCITA CI code'
      IF(DOARDUC)  WRITE(LUPRI,'(A)') ' * Run ARDUCCA CC code'
      IF(DOADC)    WRITE(LUPRI,'(A)') ' * Run RELADC code'
      IF(DOPOLPRP) WRITE(LUPRI,'(A)') ' * Run POLPRP code'
      IF(DOKRCI)   WRITE(LUPRI,'(A)')
     &   ' * Kramers restricted CI calculation'
      IF(DOKRCC)   WRITE(LUPRI,'(A)')
     &   ' * Kramers restricted CC calculation'
      IF(DOLAPLCE) WRITE(LUPRI,'(A)') ' * Run Laplace test'
C
C     Orbital manipulations
C
      IF (ROTORB) THEN
         WRITE(LUPRI,'(A)')
     &      ' * Jacobi rotations between pairs of orbitals'
         DO I = 1,NFSYM
            WRITE(LUPRI,'(A2,A1,3X,A72)') FREP(I),':',VECROT(I)
         ENDDO
      ENDIF
      IF (LMOORD) THEN
         WRITE(LUPRI,'(A)') ' * Reorder start orbitals'
         WRITE(LUPRI,'(3X,A)') '- new ordering of start orbitals...'
         DO I = 1,NFSYM
            WRITE(LUPRI,'(4X,A,I1)') 'Fermion corep. ',I
            IF (IREORD(I) .GT. 0) THEN
               WRITE(LUPRI,'(5X,50(I3,A))')
     &            (IMOORD(J,I),',',J=1,IREORD(I))
            ELSE
               WRITE(LUPRI,'(5X,A)') 'No reordering'
            END IF
         END DO
      END IF
      IF (LMOORF) THEN
         WRITE(LUPRI,'(A)') ' * Reorder converged orbitals...'
         WRITE(LUPRI,'(3X,A)') '- new ordering of converged orbitals...'
         DO I = 1,NFSYM
            WRITE(LUPRI,'(4X,A,I1)') 'Fermion corep. ',I
            IF (IREORF(I) .GT. 0) THEN
               WRITE(LUPRI,'(5X,50(I3,A))')
     &            (IMOORF(J,I),',',J=1,IREORF(I))
            ELSE
               WRITE(LUPRI,'(5X,A)') 'No reordering'
            END IF
         END DO
      END IF
C
C     Process input for various program sections
C     ==========================================
C
  200 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
     &                20), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in PSIINP.')
      ELSE
         WRITE (LUPRI,'(/3A/)')
     *      ' ERROR: Prompt in "',WORD,'" illegal or out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in PSIINP, error in prompt.')
      END IF

    2 continue
      go to 200

    3 CONTINUE
        CALL MP2INP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    4 CONTINUE
        CALL LUCITA_INP(WORD,.FALSE.)
        GO TO 200
    5 CONTINUE
        CALL MVOINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    6 CONTINUE
        CALL RESINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    7 CONTINUE
        CALL KRMCINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    8 CONTINUE
        CALL RSTPINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
C     '*OPTIMI' - control KRMCSCF optimization
    9 CONTINUE
        CALL ROPTINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
C     '*MP2 NO' - MP2 Natural Orbitals
   10 CONTINUE
        CALL MP2NOINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
C     '*KRCICALC' - Kramers Restricted CI CALCulation
   11 CONTINUE
        CALL KRCI_INP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
C
!.s/sya,2007.01.31
   12 CONTINUE
!
!       Input for COSCI module
!       Read input data
!
        CALL COSCI_INP(WORD,.FALSE.,WORK,LWORK)
!.q
        GO TO 200
C
!.s/sya,2007.02.06
   13 CONTINUE
        CALL GASCIP_INP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
C
   14 CONTINUE
#if defined (MOD_HSFS) || defined (HSFS_DYN_ALLOC)
      IF (DOHSFS) THEN
        CALL HSFSINP(WORD)
      ELSE
        WRITE(LUPRI,'(2X,A)')
     &'*HSFS input found, but not processed (because DOHSFS=.false.).'
        GOTO 1
      ENDIF
#else
C       prevent dirac.x running into an infinite loop here
C       if "*HSFS" is in the input but MOD_HSFS is not defined !!!
        CALL QUIT('*HSFS specified, '//
     &     'but HSFS module is not included in this version')
#endif
        GO TO 200
C     '*KRCCCALC' - Kramers Restricted CC CALCulation
   15 CONTINUE
#ifdef MOD_KRCC
        CALL KRCC_INP(WORD,.FALSE.,WORK,LWORK)
#else
        CALL QUIT('KRCC not activated.')
#endif
        GO TO 200
C
   16 CONTINUE
        CALL ARDINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200

!*SCF
   17 continue
      call SCFINP(word, .false., work, lwork)
      go to 200
!*KRCI
   18 continue
      CALL KRCI_INP(WORD,.FALSE.,WORK,LWORK)
      go to 200
!*QCORR
   19 continue
      ! this is a cheat so that *QCORR section is parsed
      ! by new input reader
      call move_to_next_star(word, lucmd)
      go to 200
!*LAPLCE
   20 CONTINUE
      CALL LAPINP(WORD,.FALSE.)
      GO TO 200

    1 CONTINUE
C
C     For geometry optimizations to work the INPUT routines
C     must always be called.
C


!from an email exchange: (in case you have the same question)
!radovan:     why is the DIRAC input read in twice, like here [...]
!Hans Joergen:
!             because the input routine is used to initialize its associated common block(s).
!             If the user hasn't specified this input module,
!             it still needs to be called in order to initialize with default values.
!             The SET / RESET code in the beginning of SCFINP takes care of that the variables
!             are not re-initialized if the user has specified something in input.
!             This way of doing it goes back to Trygve in ABACUS, and was adopted by Trond.
!             I have usually programmed it differently, with a separate subroutine to initialize variables.
!             You can argue for both ways of doing it.

      CALL SCFINP(WORD,.TRUE.,WORK,LWORK)
      CALL MVOINP(WORD,.TRUE.,WORK,LWORK)
      CALL MP2NOINP(WORD,.TRUE.,WORK,LWORK)
      CALL RESINP(WORD,.TRUE.,WORK,LWORK)
      CALL MP2INP(WORD,.TRUE.,WORK,LWORK)
      CALL LUCITA_INP(WORD,.TRUE.)
!
!     stefan: commented out this call. at present arducca is using the
!     same common blocks for setting input data as lucita; this is
!     WRONG!!! and causes all lucita tests to crash. :( :( :(
!
!     arducca needs to have its own input common blocks and set the
!     mandatory lucita data via an interface routine when it is taking advantage of lucita.
!     (transfer of important common block variables --> see response module for a template. response is
!     possible for both HF and MCSCF).
!
!     CALL ARDINP(WORD,.TRUE.,WORK,LWORK)

      CALL KRMCINP(WORD,.TRUE.,WORK,LWORK)
      CALL RSTPINP(WORD,.TRUE.,WORK,LWORK)
      CALL ROPTINP(WORD,.TRUE.,WORK,LWORK)

      IF(DOKRCI)THEN
        CALL KRCI_INP(WORD,.TRUE.,WORK,LWORK)
      END IF

#ifdef MOD_KRCC
      CALL KRCC_INP(WORD,.TRUE.,WORK,LWORK)
#endif

C
C     Input of relccsd module
C
      IF(DOCCM) THEN
C        DIRAC.INP is already open.
C        OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL CCINPT(LUCMD,LUPRI)
         CLOSE(LUCMD)
      ENDIF
C
C     Input of dirrci module
C
      IF(DOCIM) THEN
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL CIMINP(WORK,LWORK)
         CLOSE(LUCMD)
      ENDIF
C
C     Input of reladc module
C
      IF(DOADC) THEN
         if (parcal) then
            write (*,*) ' '
            write (*,*) ' Warning! Warning! Warning! Warning! Warning!'
            write (*,*) ' '
            write (*,*) '  RELADC does not run in parallel. To avoid'
            write (*,*) '  inconsistencies, it is suggested that jobs'
            write (*,*) '  such as this are rerun setting the number'
            write (*,*) '  of MPI processes to one (e.g. in calling'
            write (*,*) '  pam with the argument -mpi 1)'
            write (*,*) ' '
            write (*,*) '  If the 4-index transformation was performed'
            write (*,*) '  separately using scheme 6, it should also'
            write (*,*) '  be rerun with -mpi 1 (quantities from a '
            write (*,*) '  transformation using scheme 4 can be reused,'
            write (*,*) '  please refer to the section of the online'
            write (*,*) '  manual describing how to run coupled-cluster'
            write (*,*) '  calculations in separate steps for more)'
            write (*,*) ' '
            write (*,*) ' Warning! Warning! Warning! Warning! Warning!'
            write (*,*) ' '

            print *, 'controlled stop: only sequential'
            call quit ('RELADC not parallelized')
         endif
!
! MP we now comment out the default call of adcinpt since we have the
! new input structure.
!
!        OPEN(LUCMD,FILE = 'DIRAC.INP')
!        CALL ADCINPT(LUCMD,LUPRI)
!        CLOSE(LUCMD)
      ENDIF
!
!     Input of COSCI module
!     Reset and just return back.
!
      IF( DOCOSCI ) THEN
         CALL COSCI_INP(WORD,.TRUE.,WORK,LWORK)
      END IF

      IF( DOGASCIP ) THEN
         CALL GASCIP_INP(WORD,.TRUE.,WORK,LWORK)
      END IF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCITA_INP(WORD,RESET)
C***********************************************************************
C
C       Input section for LUCITA module
C
C       Written by Timo Fleig, December 2001
C       (modified SCFINP routine)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "luctinp.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
      LOGICAL   SET,RESET,NEWDEF
      CHARACTER WORD*7,WORD1*7,PROMPT*1,TABLE(NTABLE)*7
C
C
      DATA TABLE /'.TITLE ','.INIWFC','.CITYPE','.NROOTS','.SYMMET',
     &            '.NACTEL','.MULTIP','.PRINTG','.PRINTL','.SZCALC',
     &            '.INACTI','.GASSHE','.GASSPC','.FROORB','.RAS1  ',
     &            '.RAS2  ','.RAS3  ','.DENSI ','.RSTRCI','.MXCIVE',
     &            '.ANALYZ','xxxxxxx','.MAXITR','.LBLKSZ','xxxxxxx',
     &            '.DISTRT','.CONVER','.TRUNCF','.MEMFAC','XXXXXXX'/
      DATA SET/.FALSE./
C
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
C     Initialize keyword status flags driving interface.
C     (Therefore local keywords do not have to be initialized!)
C
      call izero(imokw,ntable)
C
C
C     Process input for LUCITA
C     ========================
C
      NEWDEF = (WORD .EQ. '*LUCITA')
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
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24,25,26,27,28,29,30), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &            '" not recognized in LUCITA.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in LUCITA_INP.')
C
C&&& TITLE : One title line for CI specification
    1          CONTINUE
                  READ(LUCMD,'(A)') TITLUC
                  IMOKW(1) = 1
               GO TO 100
    2          CONTINUE
C&&& INIWFC: Initial DHF wave function (closed- or open-shell)
                  READ(LUCMD,'(A)') WAFFCD
                  IMOKW(2) = 1
               GO TO 100
    3          CONTINUE
C&&& CITYPE: Type of CI calculation (FCI, SD, SDTQ, RAS, GAS)
                  READ(LUCMD,'(A)') CALCTP
                  IMOKW(3) = 1
               GO TO 100
    4          CONTINUE
C&&& NROOTS: Number of eigenvalues
                  READ(LUCMD,*) NROOTD
                  IMOKW(4) = 1
               GO TO 100
    5          CONTINUE
C&&& SYMMET: State symmetry
                  READ(LUCMD,*) ISSYMD
                  IMOKW(5) = 1
               GO TO 100
    6          CONTINUE
C&&& NACTEL: Number of active electrons
                  READ(LUCMD,*) NACTED
                  IMOKW(6) = 1
               GO TO 100
    7          CONTINUE
C&&& MULTIP: State spin multiplicity
                  READ(LUCMD,*) IMULTD
                  IMOKW(7) = 1
               GO TO 100
    8          CONTINUE
C&&& PRINTG: Global LUCIA print flag
                  READ(LUCMD,*) IPRNGD
                  IMOKW(8) = 1
               GO TO 100
    9          CONTINUE
C&&& PRINTL: Local print flag DIRLUC
                  READ(LUCMD,*) IPRNLD
                  IMOKW(9) = 1
               GO TO 100
   10          CONTINUE
C&&& SZCALC: Rough size of calculation
                  READ(LUCMD,'(A)') SZCALD
                  IMOKW(10) = 1
               GO TO 100
   11          CONTINUE
C&&& INACTI: Inactive function space
                  READ(LUCMD,'(A)') CRDINA
                  IMOKW(11) = 1
               GO TO 100
   12          CONTINUE
C&&& GASSHE: Number of GAS and distribution
                  READ(LUCMD,*) INGASD
                  DO IR = 1,INGASD,1
                    READ(LUCMD,'(A)') CRDGAS(IR)
                  END DO
                  IMOKW(12) = 1
               GO TO 100
   13          CONTINUE
C&&& GASSPC: Cumulated max. and min. GAS occupations
                  READ(LUCMD,*) NSEQCD
                  DO IR = 1,INGASD,1
                    READ(LUCMD,'(A)') CRDGOC(IR)
                  END DO
                  IMOKW(13) = 1
               GO TO 100
   14          CONTINUE
C&&& FROORB: Frozen orbital function space
                  READ(LUCMD,'(A)') CRDFRO
                  write(6,'(6A,1A72)') 'CRDFRO',CRDFRO
                  IMOKW(14) = 1
               GO TO 100
   15          CONTINUE
C&&& RAS1  : RAS1 function space
                  READ(LUCMD,'(A)') CRDRS1
                  READ(LUCMD,*) MXHL1D
                  IMOKW(15) = 1
               GO TO 100
   16          CONTINUE
C&&& RAS2  : RAS2 function space
                  READ(LUCMD,'(A)') CRDRS2
                  IMOKW(16) = 1
               GO TO 100
   17          CONTINUE
C&&& RAS3  : RAS3 function space
                  READ(LUCMD,'(A)') CRDRS3
                  READ(LUCMD,*) MXEL3D
                  IMOKW(17) = 1
               GO TO 100
   18          CONTINUE
C&&& DENSI : Level of density matrix calculation (1-,1- and 2-particle)
                  READ(LUCMD,*) IDENSD
                  IMOKW(18) = 1
               GO TO 100
   19          CONTINUE
C&&& RSTRCI: Restart option from CI vector
                  READ(LUCMD,*) IRSTLT
                  IMOKW(19) = 1
               GO TO 100
   20          CONTINUE
C&&& MXCIVE: Maximum number of subspace CI vectors
                  READ(LUCMD,*) MXCIVE
                  IMOKW(20) = 1
               GO TO 100
   21          CONTINUE
C&&& ANALYZE 
                  IMOKW(25) = 1
               GO TO 100
   22          CONTINUE
C&&& xxxxxxxxx:
               GO TO 100
   23          CONTINUE
C&&& ICIMAXITER: Maximum number of CI Iterations
                  READ(LUCMD,*) ICIMAXITER
                  IMOKW(23) = 1
               GO TO 100
   24          CONTINUE
C&&& IMAXBLKSIZE: Maximum batchsize for c-/sigma-vector
                  READ(LUCMD,*) IMAXBLKSIZE
                  IMOKW(24) = 1
               GO TO 100
   25          CONTINUE
C&&& xxxxxxxxx:
               GO TO 100
   26          CONTINUE
C&&& I_USE_DIST_ROUTE: Determination of distribution routine to use
                  READ(LUCMD,*) I_USE_DIST_ROUTE
                  IMOKW(26) = 1
               GO TO 100
   27          CONTINUE
C&&& CONVER: set convergence threshold
                  READ(LUCMD,*) my_convergence
                  IMOKW(27) = 1
               GO TO 100
   28          CONTINUE
C&&& TRUNCFactor: truncation factor for residual vectors for next trial
C                 vector
                  READ(LUCMD,*) CTRUNC_FAC
                  IMOKW(28) = 1
               GO TO 100
   29          CONTINUE
C&&& MEMFAC     : memory multiplier
                  READ(LUCMD,*) IN_MEMFAC
                  IMOKW(29) = 1
               GO TO 100
   30          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in LUCITA_INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in LUCITA_INP.')
            END IF
      END IF
  300 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ARDINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C       Input section for ARDUCCA module
C
C       Written by Lasse Sorensen
C       modified LUCTINP routine by Timo Fleig
C
C       Last revision:
C         Timo Fleig,    June 17, 2007
C             - Activation of CC linear response
C         Timo Fleig, February 8, 2011
C             - Modifications for adaptation to trunk
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "luctinp.h"
C only SPINFR used on next block
#include "dcbham.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
      LOGICAL   SET,RESET,NEWDEF
      CHARACTER WORD*7,WORD1*7,PROMPT*1,TABLE(NTABLE)*7
C
      DIMENSION WORK(LWORK)
      DIMENSION NCALCSEQ(MXNGAS)
C
C
      DATA TABLE /'.TITLE ','.INIWFC','.CITYPE','.NROOTS','.SYMMET',
     &            '.NACTEL','.MULTIP','.PRINTG','.PRINTL','.SZCALC',
     &            '.INACTI','.GASSHE','.GASSPC','.FROORB','.RAS1  ',
     &            '.RAS2  ','.RAS3  ','.DENSI ','.RSTRCI','.SEQUEN',
     &            '.CMB_CC','.NEWCCV','.RES_CC','.CCLR  ','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
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
C     Initialize keyword status flags driving interface.
C     (Therefore local keywords do not have to be initialized!)
C
      DO IKW = 1,NTABLE
         IMOKW(IKW) = 0
      END DO
C
C     Process input for ARDUCCA
C     =========================
C
      NEWDEF = (WORD .EQ. '*ARDUCC')
      ICHANG = 0
      IF (NEWDEF) THEN
C     Check that this is a spinfree calculation!
      IF (SPINFR.EQV..FALSE.) THEN
        CALL ABEND2('SPINFREE MUST BE SPECIFIED FOR ARDUCCA CALC')
      END IF
C
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
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &            '" not recognized in LUCITA.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in LUCTINP.')
C
C&&& TITLE : One title line for CI specification
    1          CONTINUE
                  READ(LUCMD,'(A)') TITLUC
                  IMOKW(1) = 1
               GO TO 100
    2          CONTINUE
C&&& INIWFC: Initial DHF wave function (closed- or open-shell)
                  READ(LUCMD,'(A)') WAFFCD
                  IMOKW(2) = 1
               GO TO 100
    3          CONTINUE
C&&& CITYPE: Type of CI calculation (FCI, SD, SDTQ, RAS, GAS)
C&&& CITYPE: Will also now contain cctype calculations added later
                  READ(LUCMD,'(A)') CALCTP
                  IMOKW(3) = 1
               GO TO 100
    4          CONTINUE
C&&& NROOTS: Number of eigenvalues
                  READ(LUCMD,*) NROOTD
                  IMOKW(4) = 1
               GO TO 100
    5          CONTINUE
C&&& SYMMET: State symmetry
                  READ(LUCMD,*) ISSYMD
                  IMOKW(5) = 1
               GO TO 100
    6          CONTINUE
C&&& NACTEL: Number of active electrons
                  READ(LUCMD,*) NACTED
                  IMOKW(6) = 1
               GO TO 100
    7          CONTINUE
C&&& MULTIP: State spin multiplicity
                  READ(LUCMD,*) IMULTD
                  IMOKW(7) = 1
               GO TO 100
    8          CONTINUE
C&&& PRINTG: Global LUCIA print flag
                  READ(LUCMD,*) IPRNGD
                  IMOKW(8) = 1
               GO TO 100
    9          CONTINUE
C&&& PRINTL: Local print flag DIRLUC
                  READ(LUCMD,*) IPRNLD
                  IMOKW(9) = 1
               GO TO 100
   10          CONTINUE
C&&& SZCALC: Rough size of calculation
                  READ(LUCMD,'(A)') SZCALD
                  IMOKW(10) = 1
               GO TO 100
   11          CONTINUE
C&&& INACTI: Inactive function space
                  READ(LUCMD,'(A)') CRDINA
                  IMOKW(11) = 1
               GO TO 100
   12          CONTINUE
C&&& GASSHE: Number of GAS and distribution
                  READ(LUCMD,*) INGASD
                  DO IR = 1,INGASD,1
                    READ(LUCMD,'(A)') CRDGAS(IR)
                  END DO
                  IMOKW(12) = 1
               GO TO 100
   13          CONTINUE
C&&& GASSPC: Cumulated max. and min. GAS occupations
C&&& Changed so several calculations can be performed in one run
                  READ(LUCMD,*) NSEQCD
                  DO IR = 1,INGASD*NSEQCD,1
                    READ(LUCMD,'(A)') CRDGOC(IR)
                  END DO
                  IMOKW(13) = 1
               GO TO 100
   14          CONTINUE
C&&& FROORB: Frozen orbital function space
                  READ(LUCMD,'(A)') CRDFRO
                  write(6,'(6A,1A72)') 'CRDFRO',CRDFRO
                  IMOKW(14) = 1
               GO TO 100
   15          CONTINUE
C&&& RAS1  : RAS1 function space
                  READ(LUCMD,'(A)') CRDRS1
                  READ(LUCMD,*) MXHL1D
                  IMOKW(15) = 1
               GO TO 100
   16          CONTINUE
C&&& RAS2  : RAS2 function space
                  READ(LUCMD,'(A)') CRDRS2
                  IMOKW(16) = 1
               GO TO 100
   17          CONTINUE
C&&& RAS3  : RAS3 function space
                  READ(LUCMD,'(A)') CRDRS3
                  READ(LUCMD,*) MXEL3D
                  IMOKW(17) = 1
               GO TO 100
   18          CONTINUE
C&&& DENSI : Level of density matrix calculation (1-,1- and 2-particle)
                  READ(LUCMD,*) IDENSD
                  IMOKW(18) = 1
               GO TO 100
   19          CONTINUE
C&&& RSTRCI: Restart option from CI vector
                  READ(LUCMD,*) IRSTLT
                  IMOKW(19) = 1
               GO TO 100
   20          CONTINUE
C&&& SEQUEN: Needed if NSEQCD > 1 hence if more than one calculation is
C&&& to be performed. This will always be the case when a cc-calculation
C&&& is desired
C&&& Notice the string in calctyp will not be decifered until readin_cc
C&&& Meaning no preliminary checking of the contents of this.
                  IUSEQ(1:6) = 'SEQUEN'
                  IF(IMOKW(13).EQ.1) THEN
                    NSEQUEN = 0
                    NTOTCALC = 0
                    DO IR = 1,NSEQCD
                      READ(LUCMD,*) NUMCALC
                      NCALCSEQ(IR) = NUMCALC
                      IF(NUMCALC.GE.1) THEN
                        DO IR2 =1,NUMCALC
                          READ(LUCMD,'(A)') CALCTYP(IR2 + NTOTCALC)
                        ENDDO
                      ELSE
                        WRITE(6,*) 'The number of different'
                        WRITE(6,*) 'calculations in this gasspace'
                        WRITE(6,*) 'must be specified'
                        write(6,*) 'check if the number of'
                        write(6,*) 'calculations in this gasspace'
                        write(6,*) 'matches',IR
C Commented out to allow for CC based on CI! Lasse
C                        CALL ABEND2('QUITTING.')
                      ENDIF
                      NTOTCALC = NUMCALC + NTOTCALC
                      NSEQUEN = 1 + NSEQUEN
                    ENDDO
                    IMOKW(20) = 1
                  ELSE
                    WRITE(6,*) 'NUMBER OF GASSPACES MUST BE SPECIFIED'
                    WRITE(6,*) 'BEFORE THE NUMBER OF SEQUENCES CAN'
                    WRITE(6,*) 'HENCE MOVE GASSPC KEYWORDS SO IT IS'
                    WRITE(6,*) 'BEFORE THE SEQUEN KEYWORD'
                    CALL ABEND2('QUITTING.')
                  ENDIF
               GO TO 100
   21          CONTINUE
C&&& CMB_CC: Use combinations for CC expansion
                  UCMBCC(1:6) = 'CMB_CC'
                  IMOKW(21) = 1
               GO TO 100
   22          CONTINUE
C&&& NEWCCV: Use new CC vector function routine
                  UNEWCCV(1:6) = 'NEWCCV'
                  IMOKW(22) = 1
               GO TO 100
   23          CONTINUE
C&&& RES_CC: Restart CC calc from CC vector
                  URES_CC(1:6) = 'RES_CC'
                  IMOKW(23) = 1
               GO TO 100
   24          CONTINUE
C&&& CCLR  : CC linear response excitation energies
                  READ(LUCMD,'(A)') CRDCCLR
                  IMOKW(24) = 1
               GO TO 100
C&&& XXXXXX: Free space for additional flag
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ARDINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ARDINP.')
            END IF
      END IF
  300 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ANWINP */
      SUBROUTINE ANWINP(WRDSRC,WORK,LWORK)
C***********************************************************************
C
C <<< Analysis Input for DIRAC >>>
C
C
C       Written by Trond Saue, April 1996
C       Last revision: April 27 1996 - tsaue
C
C***********************************************************************

      use dirac_cfg

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 12,NTABLE = 12)
#include "maxorb.h"
C
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7,WRDSRC*7
      DIMENSION WORK(LWORK)
C
#include "dcblab.h"
#include "dcbana.h"
C
      DATA TABDIR /'*END OF','*MULPOP','*PROJEC','*PRIVEC',
     &             '*RHO1  ','*DENSIT','*WT1   ','*LOCALI',
     &             '*VISUAL','*XXXXXX','*XXXXXX','*XXXXXX'/
      DATA TABLE  /'.MULPOP','.PROJEC','.PRIVEC','.RHO1  ',
     &             '.DENSIT','.WT1   ','.LOCALI','.XXXXXX',
     &             '.VISUAL','.XXXXXX','.XXXXXX','.XXXXXX'/
C
      CALL QENTER('ANWINP')
#include "memint.h"
C
      NJOBS = 0
C
C     Initialize /CBLPRP/
C     ===================
C
      DOVEC   = .FALSE.
      DOPOP   = .FALSE.
      DOPRJ   = .FALSE.
      DO1RHO  = .FALSE.
      DO3RHO  = .FALSE.
      DO1WT   = .FALSE.
      DOLOC   = .FALSE.


C
C     Read menu file
C     ==============
C
C
C     **** Find Analysis input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      CALL UPCASE(WORD)
      IF (WORD .EQ. WRDSRC) THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         WRITE(LUPRI,'(/2A)')
     &      'INFO: .ANALYZE ignored because no analysis input ',WRDSRC
         CALL QEXIT('ANWINP')
         RETURN
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no analysis input found')
  930 CONTINUE
      WORD1 = WORD
C
C     Process input for COMMON  /CBIANA/
C     ==================================
C
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,111,112),I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in ANWINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in ANWINP.')
C&&&& MULPOP: Perform Mulliken population analysis
  101    CONTINUE
            DOPOP = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
C&&&& PROJECTION: Perform projection of current solution onto another solution
  102    CONTINUE
            DOPRJ = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  103    CONTINUE
C&&&& PRIVEC: Print vectors
            DOVEC = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  104    CONTINUE
C&&&& RHO1: Punch density along the lines between centers
            DO1RHO = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  105    CONTINUE
C&&&& DENSITY: Write density to unformatted file (Gaussian cube format)
            DO3RHO = .TRUE.
            GO TO 100
  106    CONTINUE
C&&&& RHO1: Punch Becke integration weights along the lines between atomic centers
            DO1WT = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  107    CONTINUE
C&&&& LOCALIZE : Molecular orbital localization
            DOLOC = .TRUE.
            GO TO 100
  108    CONTINUE
            GO TO 100
  109    CONTINUE
!           do nothing
!           keep .VISUAL for backwards compatibility
            GO TO 100
  110    CONTINUE
            GO TO 100
  111    CONTINUE
            GO TO 100
  112    CONTINUE
            GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 180
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in ANWINP.')
      END IF
  180 CONTINUE
C
C     Print section
C     =============
C
      CALL TITLER('Analysis module','*',129)
      IF(NJOBS.GT.0) WRITE(LUPRI,'(A)') ' Jobs in this run:'
      IF(DOVEC) WRITE(LUPRI,'(1X,A)') '* Write vectors'
      IF(DOPOP) WRITE(LUPRI,'(1X,A)')
     &    '* Mulliken population analysis'
      IF(DOPRJ) WRITE(LUPRI,'(1X,A)')
     &    '* Projection onto another solution'
      IF(DO3RHO) WRITE(LUPRI,'(1X,A)')
     &    '* Write out density for visualization (Gaussian cube format)'
      IF(DO1RHO) WRITE(LUPRI,'(1X,A)')
     &    '* Punch density along lines between nuclear centers'
      IF(DO1WT) WRITE(LUPRI,'(1X,A)')
     &    '* Punch Becke weights along lines between nuclear centers'
      IF(DOLOC) WRITE(LUPRI,'(1X,A)')
     S    '* Molecular orbital localization'
C
C
C     Process input for various program sections
C     ==========================================
C
  200 CONTINUE
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4,5,6,7,8,9,10,11,12), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in ANWINP.')
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in ANWINP, error in prompt.')
      END IF
    2 CONTINUE
        CALL POPINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    3 CONTINUE
        CALL PRJINP(WORD,.FALSE.)
        GO TO 200
    4 CONTINUE
        CALL VECINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    5 CONTINUE
        CALL RHO1IN(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    6 CONTINUE
        CALL DENSIN(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    7 CONTINUE
        CALL WT1IN(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    8 CONTINUE
        CALL LOCINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200

!     *VISUAL
    9 continue
!     this is here for backwards compatibility, *VISUAL used to be under **ANALYZE
      call move_to_next_star(word, lucmd)
      go to 200

   10 CONTINUE
        GO TO 200
   11 CONTINUE
        GO TO 200
   12 CONTINUE
        GO TO 200
C
    1 CONTINUE
      IF(DOPOP) CALL POPINP(WORD,.TRUE.,WORK,LWORK)
      IF(DOPRJ) CALL PRJINP(WORD,.TRUE.)
      IF(DOVEC) CALL VECINP(WORD,.TRUE.,WORK,LWORK)
      IF(DO1RHO) CALL RHO1IN(WORD,.TRUE.,WORK,LWORK)
      IF(DO3RHO) CALL DENSIN(WORD,.TRUE.,WORK,LWORK)
      IF(DO1WT) CALL WT1IN(WORD,.TRUE.,WORK,LWORK)
      IF(DOLOC) CALL LOCINP(WORD,.TRUE.,WORK,LWORK)
C
      CALL QEXIT('ANWINP')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck popinp */
      SUBROUTINE POPINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for Mulliken population analysis
C
C     Written by T.Saue - April 1996
C     Last revision: April 27 - 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbpop.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbmp2no.h"
C
      LOGICAL   SET, NEWDEF,ORBDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, LINE*72
      DIMENSION WORK(LWORK)
C
! HJAAJ TODO : new parameter to set ADD_ALL = .TRUE. /Jan 2010
      SAVE SET
      DATA TABLE /'.AOLAB ','.LABDEF','.VECPOP','.LABEL ',
     &            '.PRINT ','.INDSML','.NETPOP','.XXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
#include "memint.h"
C
      LABOPT = 0
C
C     Initialize /CBIPOP/
C     ===================
C
      IPRPOP   = 0
      ILABDF   = 2
      DONETP   = .FALSE.
      ADD_ALL  = .FALSE.
      ADDSML   = .TRUE.
      ORBDEF   = .TRUE.
CSK   ... only if MP2-NO orbitals are available - set in MP2NOGE1
      MP2NATPOP = .FALSE.
C
C     Process input from CBIPOP
C     =========================
C
      NEWDEF = (WORD .EQ. '*MULPOP')
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
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in POPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in POPINP.')
    1          CONTINUE
C&&&& AOLAB:  Base definition of labels on AO-labels:
                 ILABDF = 1
               GO TO 100
    2          CONTINUE
C&&&& LABDEF: Define groups of primitive labels
                 LABDEF = .TRUE.
                 READ(LUCMD,*) NPOPLAB
                 CALL MEMGET('INTE',KIBUF,MAXLAB,WORK,KFREE,LFREE)
                 CALL PROLAB(IPOPLAB,POPLAB,MAXLAB,NPOPLAB,LINE,LUCMD,
     &                       WORK(KIBUF),WORK(KFREE),LFREE)
                 CALL MEMREL('LABDEF',WORK,KWORK,KWORK,KFREE,LFREE)
               GO TO 100
    3          CONTINUE
C&&&& VECPOP: Number of spinors to analyze
                 READ(LUCMD,'(A72)') (VECPOP(I),I=1,NFSYM)
                 ORBDEF = .FALSE.
                 IPRPOP = MAX(IPRPOP,1)
               GO TO 100
    4          CONTINUE
C&&&& LABEL: pre-define labels
                 LABDEF = .TRUE.
                 READ(LUCMD,'(A7)') WORD
                 IF    (WORD.EQ.'ATOM   ') THEN
                   LABOPT = 1
                 ELSEIF(WORD.EQ.'SHELL  ') THEN
                   LABOPT = 2
                 ELSE
                   CALL QUIT('*POPINP: Unknown word after LABEL')
                 ENDIF
               GO TO 100
    5          CONTINUE
C&&&& IPRPOP: Print level
                  READ(LUCMD,*) IPRPOP
               GO TO 100
    6          CONTINUE
C&&&& INDSML: Do not collect labels of small components of given center/symmetry
                  ADDSML = .FALSE.
               GO TO 100
    7          CONTINUE
C&&&& NETPOP: Mulliken Net populations
                 DONETP = .TRUE.
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in POPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in POPINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C
      IF(ORBDEF) THEN
        DO I = 1,NFSYM
          VECPOP(I) = 'UNDEFINED'
          IF(NOCC(I).GT.0) THEN
            IF(NOCC(I).EQ.1) THEN
              WRITE(VECPOP(I),'(I4)') 1
            ELSE
              WRITE(VECPOP(I),'(I4,A2,I4)') 1,'..',NOCC(I)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      IF(LABOPT.GT.0) CALL PRELAB(LABOPT,ILABDF,IPOPLAB,POPLAB,NPOPLAB)
C
C     Print section
C     =============
C
      IF(.NOT.DOPOP) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'POPINP: Mulliken population analysis'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') '* Gross populations'
      IF(DONETP) WRITE(LUPRI,'(1X,A)')
     &     '* Net populations and overlap'
      IF    (LABOPT.EQ.1) THEN
        WRITE(LUPRI,'(1X,A)') '* Label defined for individual atoms'
      ELSEIF(LABOPT.EQ.2) THEN
        WRITE(LUPRI,'(1X,A)') 
     &  '* Label defined for individual orbital shells'
      ENDIF
      IF(ILABDF.EQ.1) THEN
        WRITE(LUPRI,'(1X,A)') '* Label definitions based on AO-labels'
      ELSE
        WRITE(LUPRI,'(1X,A)') '* Label definitions based on SO-labels'
      ENDIF
      IF(.NOT.ADDSML) WRITE(LUPRI,'(3X,A)')
     &      '(individual small component labels will be used)'
      WRITE(LUPRI,'(1X,A)') '* Number of spinors analyzed:'
      DO I = 1,NFSYM
        NVEC = 0
        IF (VECPOP(I) .eq. 'UNDEFINED') THEN
          WRITE(LUPRI,'(4X,A,A3)')
     &      '- All occupied orbitals in fermion ircop ',FREP(I)
        ELSE
          CALL  NUMLST(VECPOP(I),IDUMMY,NFBAS(I,0),
     &                -NFBAS(I,2),NFBAS(I,1),I,NVEC)
          IF(NVEC.EQ.0) THEN
            WRITE(LUPRI,'(4X,A,A3)')
     &      '- No orbitals in fermion ircop ',FREP(I)
          ELSE
            WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &      '- Orbitals in fermion ircop ',FREP(I),' :',VECPOP(I)
          ENDIF
        ENDIF
      ENDDO
      WRITE(LUPRI,'(1X,A,I5)') '* Print level:',IPRPOP
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prjinp */
      SUBROUTINE PRJINP(WORD,RESET)
C***********************************************************************
C
C     Input section for module for wavefunction projection
C
C     Written by T.Saue - April 1996
C     Last revision: Apr 26 - 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 12)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbprj.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "nuclei.h"
C
#include "dcbxpr.h"
      LOGICAL SET, NEWDEF,ORBDEF,RESET,SYMSTP
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      SAVE SET
      DATA TABLE /'.PRINT ','.PROTHR','.VECPRJ','.VECREF',
     &            '.WGPOP ','.OWNBAS','.PRJADJ','.POLREF',
     &            '.MAYER ','.POLNRM','.ATOMS ','.XXXXXX'/
      DATA SET/.FALSE./
C
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
C     Initialize /CBIPRJ/
C     ===================
C
      IPRPRJ   = 0
C     Default is to project all occupied vectors onto all fragment vectors
      ORBDEF = .TRUE.
      NREFS    = 0
      PROTHR   = 1.0D-3
      LWGPOP = .FALSE.
      OWNBAS = .FALSE.
      PATOMS = .FALSE.
      PRJADJ = .FALSE.
      POLREF = .FALSE.
      MAYER  = .FALSE.
      POLNRM = .FALSE.
C
C     Process input from CBIPRJ
C     =========================
C
      NEWDEF = (WORD .EQ. '*PROJEC')
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
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized under *PROJEC.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword under *PROJEC.')
    1          CONTINUE
C&&&& PRINT: Print level for projections
                 READ(LUCMD,*) IPRPRJ
               GO TO 100
    2          CONTINUE
C&&&& PROTHR: Threshold for printout of projection coefficients
                 READ(LUCMD,*) PROTHR
               GO TO 100
    3          CONTINUE
C&&&& VECPRJ: Spinors to project onto fragment vectors
                  READ(LUCMD,'(A72)') (VECPRJ(I),I=1,NFSYM)
                  ORBDEF = .FALSE.
               GO TO 100
    4          CONTINUE
C&&&& VECREF: Fragment vectors
                  READ(LUCMD,*) NREFS
                  IF((NREFS+1).GT.MAXREF) THEN
                    WRITE(LUPRI,'(A,I5)')
     &         '* Number of fragment specified:',NREFS,
     &         '* Current maximum             :',MAXREF
                    CALL QUIT('*PROJECT input: Too many fragments !')
                  ENDIF
                  DO J = 1,NREFS
                    READ(LUCMD,'(A6)', IOSTAT=IOS1) REFFIL(J)
                    IF (IOS1.NE.0) THEN ! miro: specific error message
                      call quit('error in REFFIL(J) reading!')
                    ENDIF
                    READ(LUCMD, '(I8)', IOSTAT=IOS2) NPROJNUC(J)
                    IF (IOS2.NE.0) THEN ! miro: specific error message
                      write(LUPRI,*) 'J,NREFS:',J,NREFS
                      call quit('error in NPROJNUC(J) reading!')
                    ENDIF
                    READ(LUCMD,'(A72)', IOSTAT=IOS3) 
     &              (VECREF(I,J),I=1,NFSYM)
                    IF (IOS3.NE.0) THEN ! miro: specific error message
                      call quit('error in VECREF(I,J) reading!')
                    ENDIF

                 ENDDO
                  REFFIL(NREFS+1) = 'Polari'
               GO TO 100
    5          CONTINUE
C&&& WGPOP: Split overlap densities according to weight of contributions
                 LWGPOP = .TRUE.
               GO TO 100
    6          CONTINUE
C&&& OWNBAS: Fragments (defined from NUCIND)  are calculated in their own basis sets
                 OWNBAS = .TRUE.
               GO TO 100
    7          CONTINUE
C&&& PRJADJ: Phase adjustment according to selected reference orbital
                  PRJADJ = .TRUE.
                  READ(LUCMD,*) (IREFADJ(I),I=1,NFSYM)
               GO TO 100
    8          CONTINUE
C&&& POLREF: Polarize reference orbitals
                  POLREF = .TRUE.
               GO TO 100
 9             CONTINUE
C&&& MAYER: Polarize reference orbitals according to Mayers scheme
                  MAYER  = .TRUE.
                  OWNBAS = .TRUE.
               GO TO 100
 10            CONTINUE
C&&& POLNRM: Normalize polarization contribution
                  POLNRM = .TRUE.
               GO TO 100
 11            CONTINUE
C&&& ATOMS : Do projection analysis using atomic fragments
                  PATOMS = .TRUE.
                  IREFS=1
                  ININD=1
                  DO ITYP = 1,NONTYP
                    IF(.NOT.NOORBT(ININD)) THEN
                      READ(LUCMD,'(A6)') REFFIL(IREFS)
                      READ(LUCMD,'(A72)') VECREF(1,IREFS)
                      NPROJNUC(IREFS)=1
                      DO J = 0,NONT(ITYP)-1
                        IREFS = IREFS + NUCDEG(ININD+J)
                      ENDDO
                    ENDIF
                    ININD = ININD + NONT(ITYP)
                  ENDDO
                  IF(IREFS.GT.MAXREF) THEN
                    WRITE(LUPRI,'(A,I5)')
     &         '* Number of atomic types specified  :',IREFS,
     &         '* Current maximum no. of fragments  :',MAXREF
                    CALL QUIT('*PROJECT input: Too many fragments !')
                  ENDIF
                  REFFIL(IREFS) = 'Polari'
                  NREFS=IREFS-1
                  GO TO 100
 12            CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized under *PROJEC.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt under *PROJEC.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C
      IF(ORBDEF) THEN
        DO I = 1,NFSYM
          VECPRJ(I) = ' '
          IF(NOCC(I).GT.0) THEN
            IF(NOCC(I).EQ.1) THEN
              WRITE(VECPRJ(I),'(I4)') 1
            ELSE
              WRITE(VECPRJ(I),'(I4,A2,I4)') 1,'..',NOCC(I)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C     Print section
C     =============
C
      IF(.NOT.DOPRJ) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') '*PROJECT: Projection onto another solution'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,1P,D9.2)')
     &       '* Threshold for printing projection coefficients:',PROTHR
      WRITE(LUPRI,'(1X,A)') '* Number of spinors projected:'
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECPRJ(I),IDUMMY,NFBAS(I,0),
     &              -NFBAS(I,2),NFBAS(I,1),I,NVEC)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A,A3)')
     &       '- No orbitals in fermion ircop ',FREP(I)
        ELSE
          WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &      '- Orbitals in fermion ircop ',FREP(I),' :',VECPRJ(I)
        ENDIF
      ENDDO
      IF(POLREF) THEN
          WRITE(LUPRI,'(1X,A)') '* Reference orbitals will be polarized'
      ENDIF
      IF(MAYER) THEN
        IF(OWNBAS) THEN
          CALL QUIT
     &     ('ANAINP: Mayer polarization only possible with SELOWN=T')
        ELSE
          WRITE(LUPRI,'(1X,A,A)')
     &      '* Reference orbitals will be polarized according to ',
     &      ' a modified Mayer scheme'
        ENDIF
      ENDIF
      IF(OWNBAS) THEN
        WRITE(LUPRI,'(1X,A)')
     &   '* Fragments are calculated in their individual bases !',
     &   '    The indexing of fragment spinor sets is assumed ',
     &   '    to follow list of symmetry independent nuclei.'
        DO J = 1,NREFS
          WRITE(LUPRI,'(1X,A,I3,A,A6)') '* Fragment spinor set ',J,
     &           ' --> ',NAMDEP(J)
          WRITE(LUPRI,'(1X,A,A6)') ' - read from file ',REFFIL(J)
          DO I = 1,NFSYM
            NVEC = 0
            CALL  NUMLST(VECREF(I,J),IDUMMY,NFBAS(I,0),
     &                -NFBAS(I,2),NFBAS(I,1),I,NVEC)
            IF(NVEC.EQ.0) THEN
              WRITE(LUPRI,'(4X,A,A3)')
     &           '- No orbitals in fermion ircop ',FREP(I)
            ELSE
              WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &        '- Orbitals in fermion ircop ',FREP(I),' :',VECREF(I,J)
            ENDIF
          ENDDO
        ENDDO
      ELSEIF(PATOMS) THEN
        WRITE(LUPRI,'(1X,A)')
     &   '* Projection analysis based on atomic fragments, ',
     &   '  calculated in their individual bases'
        IREFS=1
        ININD=1
        DO ITYP = 1,NONTYP
          IF(.NOT.NOORBT(ININD)) THEN
            WRITE(LUPRI,'(6X,A6,3X,A72)') REFFIL(IREFS),VECREF(1,IREFS)
            DO J = 0,NONT(ITYP)-1
              IREFS = IREFS + NUCDEG(ININD+J)
            ENDDO
          ENDIF 
          ININD = ININD + NONT(ITYP)
        ENDDO
      ELSE
        DO J = 1,NREFS
          WRITE(LUPRI,'(1X,A,I3)') '* Fragment spinor set ',J
          WRITE(LUPRI,'(1X,A,A6)') ' - read from file ',REFFIL(J)
          DO I = 1,NFSYM
            NVEC = 0
            CALL  NUMLST(VECREF(I,J),IDUMMY,NFBAS(I,0),
     &                -NFBAS(I,2),NFBAS(I,1),I,NVEC)
            IF(NVEC.EQ.0) THEN
              WRITE(LUPRI,'(4X,A,A3)')
     &           '- No orbitals in fermion ircop ',FREP(I)
            ELSE
              WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &        '- Orbitals in fermion ircop ',FREP(I),' :',VECREF(I,J)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      IF(POLNRM) THEN
        WRITE(6,'(1X,A)')
     &    '* Polarization contributions will be normalized.'
      ENDIF
      IF(LWGPOP) WRITE(LUPRI,'(A,A)') '* Overlap densities split ',
     &  'between contributions according to their weight.'
      WRITE(LUPRI,'(1X,A,I5)') '* Print level:',IPRPRJ
      IF(PRJADJ) THEN
         WRITE(LUPRI,'(A,A)') 'Molecular coefficients phase adjusted ',
     &      'according to reference orbitals :'
         WRITE(LUPRI,'(A,A3,A,I5)')
     &      ('Fermion ircop ',FREP(I),': ',IREFADJ(I),I=1,NFSYM)
      ENDIF
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vecinp */
      SUBROUTINE VECINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for vector print
C
C     Written by T.Saue - May 1996
C     Last revision: May 2 - 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbvec.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.PRICMP','.AOLAB ','.VECPRI','.PRINT ',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
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
C     Initialize /CBIVEC/
C     ===================
C
      IPRCMP    = 1
      DO I = 1,NFSYM
        VECPRI(I) = ' '
        IF(NOCC(I).EQ.1) THEN
          WRITE(VECPRI(I),'(I1)') 1
        ELSEIF(NOCC(I).GT.1) THEN
          WRITE(VECPRI(I),'(I1,A2,I4)') 1,'..',NOCC(I)
        ENDIF
      ENDDO
      ILABDF  = 2
      IPRVEC  = 0
C
C     Process input from CBIPOP
C     =========================
C
      NEWDEF = (WORD .EQ. '*PRIVEC')
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
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in VECINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in VECINP.')
    1          CONTINUE
C&&&& PRICMP: Print control for large and small component vectors
                  READ(LUCMD,*) IL,IS
                  IPRCMP = IL+2*IS
               GO TO 100
    2          CONTINUE
C&&&& AOLAB : Print coefficients in AO-basis
                  ILABDF = 1
               GO TO 100
    3          CONTINUE
C&&&& VECPRI: Number of spinors to print
                  READ(LUCMD,'(A72)') (VECPRI(I),I=1,NFSYM)
               GO TO 100
    4          CONTINUE
C&&&& PRINT: General print level
                  READ(LUCMD,*) IPRVEC
               GO TO 100
    5          CONTINUE
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in VECINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in VECINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT.DOVEC) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'VECINP: Vector print'
      CALL PRSYMB(LUPRI,'=',75,0)
      IF(ILABDF.EQ.1) THEN
        WRITE(LUPRI,'(A)') ' * Coefficients written in AO-basis'
      ELSE
        WRITE(LUPRI,'(A)') ' * Coefficients written in SO-basis'
      ENDIF
      WRITE(LUPRI,'(A)')   ' * Vector print:'
      NTOT = 0
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECPRI(I),IDUMMY,NFBAS(I,0),
     &              -NFBAS(I,2),NFBAS(I,1),I,NVEC)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A,A3)')
     &       '- No orbital string specified in fermion ircop ',FREP(I)
        ELSE
          NTOT = NTOT + NVEC
          WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &      '- Orbitals in fermion ircop ',FREP(I),' :',VECPRI(I)
        ENDIF
      ENDDO
Cluuk IF(NTOT.EQ.0) IPRCMP = 0
      IF(IPRCMP.EQ.1) WRITE(LUPRI,'(1X,A)')
     +  '* Only large component coefficients written out.'
      IF(IPRCMP.EQ.2) WRITE(LUPRI,'(1X,A)')
     +  '* Only small component coefficients written out.'
      IF(IPRCMP.EQ.3) WRITE(LUPRI,'(1X,A)')
     +  '* Both large and small component coefficients written out.'
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PAMINI()
C***********************************************************************
C
C     Initialize data in DIRAC
C     Written by T.Saue May 14 1996
C     Last revision: May 14 1996 - tsaue
C
C***********************************************************************
      use codata
#include "implicit.h"
#include "priunit.h"
#include "consts.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "siripc.h"
#include "huckel.h"
#include "dcbgen.h"
#include "ccom.h"
#include "cbirea.h"
#include "dcbxpr.h"
#include "dcbcls.h"
#include "dcbprl.h"
C
      call set_codata_values(CODSET)
C
C     Initialize temporary variables
C     ==============================
C
      ILSPH  = 1
      ISSPH  = 1
      ILLDIR = 1
      ISLDIR = 1
      ISSDIR = 1
      IGTDIR = 1
      ILLINT = 1
      ISLINT = 1
      ISSINT = 1
      IGTINT = 0
C
C     Initialize /DCBGEN/
C     ===================
C
      TITLE    = 'DIRAC: No title specified !!!'
      FAMILY   = .FALSE.
      INPTES   = .FALSE.
      VACUUM   = .FALSE.
      NUCVAC   = .FALSE.
      DOPSI    = .FALSE.
      DOANA    = .FALSE.
      DOPRP    = .FALSE.
      DOHRM    = .FALSE.
      DOJACO   = .FALSE.
      DOPUT    = .FALSE.
      DOACUT   = .FALSE.
      DOACIN   = .FALSE.
      DORKBIMP = .FALSE.
      DOTRA    = .FALSE.
      NOTRA    = .FALSE.
      OPTIMI   = .FALSE.
      NOSET    = .FALSE.
      OPTWLK   = .FALSE.
      OPTNEW   = .FALSE.
      NMWALK   = .FALSE.
      DOHUCKEL = .false. ! TODO : make DOHUCKEL default
!
      IPRGEN   = 0
      ISPHTR   = ILSPH+2*ISSPH
      IDFLAG   = ILLDIR+2*ISLDIR+4*ISSDIR+8*IGTDIR
      INTGEN   = ILLINT+2*ISLINT+4*ISSINT+8*IGTINT
      INTGEN_SAVE = INTGEN
      CVAL     = CVEL
      STOL(1)  = 1.0D-6
      STOL(2)  = 1.0D-8
      PANAS  = D0
      HRINPC = .FALSE.
      RDINPC = .FALSE.
      IF(IDFLAG.GT.0) THEN
        DIRCAL = .TRUE.
      ELSE
        DIRCAL = .FALSE.
      ENDIF
C --->
C
C     Test : disable automatic conversional evaluation of
C            one-center LS and SS integrals in ONECAP models.
C            These can also take up quite some space !!!   /jkp
C
      DIRSET = .TRUE.
C <---
      IPRUSR = 0
      USRIPR = .FALSE.
      SEGBAS = .TRUE.
      NEWBAS = .TRUE.
      NEWPRP = .TRUE.
      RELCAL = .TRUE.
C
C     Initialize CNTMAT from cbirea
C
      CNTMAT = .FALSE.
C
C     Initialize SIRIPC from DALTON
C
      NEWGEO = .TRUE.
      INPPRC = .FALSE.
C
C     Initialize /HUCKEL/from DALTON
C
      IPRHUC = 1
      HUCCNT = 1.75D0
      CALL IZERO(NHUCAO,8)
      CALL IZERO(IHUCPT,MXSHEL)
C
C     Initialize XCBCLS - list of integrals to calculate
C     ==================================================
      NPRPCLS   = 0
C
C     Initialize XCBPRL - list of integral labels
C     ===========================================
C
      NPRPLBL    = 0
C
C     Initialize XCBXPR - list of properties
C     ======================================
      NPRPS        = 0
C
      CALL OPTINI()
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HAMINP(WORK,LWORK)
C***********************************************************************
C
C     Input section for definition of the Hamiltonian
C
C     Called from: PAMINP
C
C     Written by T.Saue - May 1996
C
C***********************************************************************
      use dirac_cfg
      use fde_mod, only: fde_input_init,fde_dirac_input
      use dft_cfg
      use xc_derv
      use interface_functional_read
      use x2cmod_cfg
#ifdef HAS_PCMSOLVER      
      use pcmmod_cfg
#endif
#ifdef HAS_PELIB
      use pe_variables, only: peqm
#endif
#ifdef MOD_AOOSOC
      use atomic_oo_order_so_correction_cfg
#endif

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 76,NDIR = 12)
      PARAMETER (D0 = 0.00D00, D1 = 1.0D0,DM2 = -2.0D0,D2 = 2.00D00)
C
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "gencon.h"
#include "nuclei.h"
#include "orgcom.h"
C
      LOGICAL  SYMSTP,TRSSTP,NOINPUT,SETFUN,SUCCESS,IOTC,GAUNT_SAVE,
     &         TWOCOMP_SAVE, TWOCOMPBSS_SAVE, find_keyword,INTFLG_CHANGE
      CHARACTER PROMPT*1, WORD*7, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7, LINE*80, LINE100*100
      DIMENSION KSOP(0:7),WORK(LWORK)
      real(8) :: hfx_out, mu_out, beta_out
C
      SAVE SETFUN
      DATA TABDIR /'*END OF','*SOLVEN','*DFT   ','*CAP   ',
     &             '*FDE   ','*AMFI  ','*X2C   ','*ECP   ',
     &             '*PCM   ','*aooSOC','*PEQM  ','*PCMSOL'/

!                    1         2         3         4
      DATA TABLE /'.PRINT ','.ONESYS','.OPERAT','.GAUNT ',
!                    5         6         7         8
     &            '.LVCORR','.LVNEW ','.SOLVEN','.FREEPJ',
!                    9        10        11        12
     &            '.VEXTPJ','.URKBAL','.NOSMLV','.LEVY-L',
!                   13        14        15        16
     &            '.SPINFR','.ZORA  ','.SMLV1C','.ALDA  ',
!                   17        18        19        20
     &            '.DFT   ','.HFXFAC','.SCQSET','.ONECNV',
!                   21        22        23        24
     &            '.ONECAP','.SPINF2','.NOSPIN','.USE_DF',
!                   25        26        27        28
     &            '.INTFLG','.BSS   ','.HFXMU ','.YREQ1 ',
!                   29        30        31        32
     &            '.ONESTE','.BLOCKD','.DO2C4C','.DO4C2C',
!                   33        34        35        36
     &            '.CONT2C','.BEG_2C','.VACPOL','.X2COLD',
!                   37        38        39        40
     &            '.X2C   ','.MIXVAC','.CAP   ','.NOVREF',
!                   41        42        43        44
     &            '.POSVAC','.FDE   ','.IOTC4 ','.NONREL',
!                   45        46        47        48
     &            '.BSS4  ','.CMPEIG','.DKH2  ','.DKH1  ',
!                   49        50        51        52
     &            '.NOAMFI','.QDOTS ','.IOTC  ','.X2C4  ',
!                   53        54        55        56
     &            '.NOSFMU','.DFTAUT','.MMF   ','.MO4C2C',
!                   57        58        59        60
     &            '.PCM   ','.HFXATT','.XB    ','.XC    ',
!                   61        62        63        64
     &            '.ECP   ','.X2CMMF','.DOSSSS','.MDIRAC',
!                   65        66        67        68
     &            '.TRSMIX','.PEQM  ','.JZOUT ', '.FAKE2C',
!                   69        70        71        72
     &            '.GAUGEO','.GO ANG','.DIPORG', '.PHASEO',
!                   73        74        75        76
     &            '.KOUT  ','.XXXXXX','.XXXXXX', '.XXXXXX'/
      DATA SETFUN/.FALSE./
C
#include "ibtfun.h"
C
      CALL QENTER('HAMINP')

      GAUNT = .FALSE.
      GAUNT_SAVE = GAUNT
      INTFLG_CHANGE = .FALSE.
      IOTC = .FALSE.
      SYMSTP = .FALSE.
      TRSSTP = .FALSE.
      INPERR = 0
      ICHANG  = 0 ! mi-to satisfy ifort runtime check
      CALL TITLER('Hamiltonian defined','*',127)
C
C
C     Initialize /CBIHAM/
C     ===================
C
      IPRHAM  = 0
      DOLVC   = .TRUE.
      QED     = .FALSE.
      MOLMF   = .FALSE.
      X2CMMF  = .FALSE.
      ONESYS  = .FALSE.
      FREEPJ  = .FALSE.
      VEXTPJ  = .FALSE.
      URKBAL  = .FALSE.
      NOSMLV  = .FALSE.
      MDIRAC  = .FALSE.
      LEVYLE  = .FALSE.
      NONREL  = .FALSE.
      FAKE2C  = .FALSE.
      SPINFR  = .FALSE.
      NOSPIN  = .FALSE.
      SPINFR2 = .FALSE.
      LVNEW   = .FALSE.
      SOLVEN  = .FALSE.
      ZORA    = .FALSE.
      ZORA4   = .FALSE.
      ZORASC  = .FALSE.
      ONECAP  = .FALSE.
      SMLV1C  = .FALSE.
      ONECOFF = .FALSE.
      VACPOL  = .FALSE.
      MIXVAC  = .FALSE.
      POSVAC  = .FALSE.
      VACREF  = .TRUE.
      CAP     = .FALSE.
!     ECPCALC = .FALSE. - already set in HAMSCAN
      ORIGIN(1:3) = 0.0d0
      GAGORG(1:3) = 0.0d0
      GAGORG_SET = .FALSE.
      DIPORG(1:3) = 0.0d0
      INTV1C  = 0
      ONECNV  = D0
      N1OPER  = 0
      IPOVRLAP = -1
      IPBETAMT = -1
      IPMOLFLD = -1
      IPKINERG = -1
      IPANGMOM = -1
      IPSPNMOM = -1
      IPSPNORB = -1
      IPVEMB0  = -1
      SSMTRC  = D1
      BSS     = .FALSE.
      x2c     = .false.
!     default for X2C module: add spin-same orbit MFSSO2 terms
      isorder_amfi_x2c = 2
      IBSS    = -1
      IMFCH   = 0
CMI ... default for BSS - do always the preliminary free particle basis transformation
      NOPRTR  = .FALSE.
C      ... for the infinite order solve the R-equation 2
      YREQ1   = .FALSE.
C      ... default - do not perform the numerical block diagonalization !
      BLOCKD  = .FALSE.
C     ... trasfer 2c-> 4c and 4c->2c
      DO2C4C  = .FALSE.
      START2C = .FALSE.
      INI2C   = -1
      DO4C2C  = .FALSE.
C     .... flag for 4cMO -> 2cMO transform.
      TRMO4C2C = .FALSE.

      TWOCOMP_SAVE    = TWOCOMP
      TWOCOMPBSS_SAVE = TWOCOMPBSS

      USEDF   = .FALSE.
      CONT2C  = .FALSE.
      CMPEIG  = .FALSE.
CMI   ... own keywords for 1. or 2.order Douglas-Kroll-Hess with AMFI SSO
      DKH1    = .FALSE.
      DKH2    = .FALSE.
      I2CHAM  = -1
C     ... default L&S components... needs discussion about relation to DO2C4C/DO4C2C
      MC      = 2
CMI     ... default initialization - always include AMFI contribution
      NOAMFI = .FALSE.
CMI   ... default value - no artificial atom
      QDOTS = .FALSE.
Culfek: Do not scale gaunt integrals with exchange scale factor by default.
      LSCALE_DFT_GAUNT = .FALSE.
      NOSFMU = .FALSE.
Ctsaue: Exchange factor must not be reset during geometry optimizations.
      IF (.NOT. SETFUN) THEN
        HFXFAC  = D1
        HFXMU   = D0
        HFXATT  = D0
      ENDIF
C ... get integral flags in advance
      ILLINT    = IBTAND(INTGEN,1)
      ISLINT    = IBTAND(INTGEN/2,1)
      ISSINT    = IBTAND(INTGEN/4,1)
      IGTINT    = IBTAND(INTGEN/8,1)
!mi:  by default do not allow mixing T+ with T-
      TRS_MIX=.FALSE.
!jmho: deafult - no polarizable embedding
#ifdef HAS_PELIB
      PEQM = .FALSE.
#endif
!mi:  by default, do not save Jz_MO matrices
      WRITE_Jz_MATRIX = .FALSE.
      WRITE_K_MATRIX = .FALSE. 
      NOINPUT   = .FALSE.
C
C     Read menu file
C     ==============
C     **** Find Hamiltonian input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      CALL UPCASE(WORD)
      IF (WORD .EQ. '**HAMIL') THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         NOINPUT = .TRUE.
         GOTO 300
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no HAMIL input found')
  930 CONTINUE
      WORD1 = WORD(1:1)
C
C     Process input for COMMON  /CBIHAM/
C     ==================================
C
C     Keep track of number of changes and added finite field operators
      I_FINFIELD = 0
      ICHANG  = 0
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112,113,114,115,116,117,118,119,120,
     &                121,122,123,124,125,126,127,128,129,130,
     &                131,132,133,134,135,136,137,138,139,140,
     &                141,142,143,144,145,146,147,148,149,150,
     &                151,152,153,154,155,156,157,158,159,160,
     &                161,162,163,164,165,166,167,168,169,170,
     &                171,172,173,174,175,176),I
            END IF
  99    CONTINUE
        IF (WORD .EQ. '.OPTION') THEN
           CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
           CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
           GO TO 100
         END IF
         WRITE (LUPRI,'(/,3A,/)')
     *       ' Keyword ',WORD,' not recognized under **HAMILTONIAN.'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword under **HAMILTONIAN.')
  101    CONTINUE
C&&&& PRINT: Print level for Hamiltonian
           READ(LUCMD,*) IPRHAM
           GO TO 100
  102    CONTINUE
C&&&& ONESYS: Ignore two-electron operator
           ONESYS = .TRUE.
           ICHANG = ICHANG + 1
           GO TO 100
  103    CONTINUE
C&&&& OPERATOR: Read in additional  operator, the print level is IPRGEN from the section **GENERAL
           CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                 IPRGEN)
           N1OPER = N1OPER + 1
           IF(N1OPER.GT.MAX1OPER) THEN
             WRITE(LUPRI,'(A)')
     &       '***** **HAMILTONIAN ERROR *****',
     &       'Redimension MAX1OPER of /CBIHAM/ !'
             CALL QUIT('**HAMILTONIAN: beyond MAX1OPER !')
           ENDIF
           IPR1OP(N1OPER) = INDXPR
           ICHANG = ICHANG + 1
           I_FINFIELD = I_FINFIELD + 1
           GO TO 100
  104    CONTINUE
C&&&& GAUNT: Gaunt operator
            GAUNT = .TRUE.
            GO TO 100
  105     CONTINUE
C&&&& LVCORR: Modeling of interatomic SS-integral contribution by classical
C             repulsion of small component atomic charges
            DOLVC = .TRUE.
            GO TO 100
  106     CONTINUE
C&&&& LVNEW: Get Small component charge via Mulliken analysis
            DOLVC = .TRUE.
            LVNEW = .TRUE.
            GO TO 100
  107     CONTINUE
C&&& SOLVENT: Solvent effects modelled by spherical cavity in dielectric medium
            SOLVEN = .TRUE.
            GO TO 100
  108     CONTINUE
C&&& FREEPJ: Reduce variational space by projecting out free positronic solutions
             FREEPJ = .TRUE.
             ICHANG = ICHANG + 1
            GO TO 100
  109     CONTINUE
C&&& VEXTPJ: Reduce variational space by projecting out positronic solutions of
C&&&         external field
             VEXTPJ = .TRUE.
             ICHANG = ICHANG + 1
             GO TO 100
  110      CONTINUE
C&&&& URKBAL : No restriction on kinetic balance; Unphysical solutions not removed
             URKBAL = .TRUE.
             GO TO 100
  111      CONTINUE
C&&& NOSMLV: Turn off small nuclear attraction integrals
             NOSMLV = .TRUE.
             ICHANG = ICHANG + 1
             GO TO 100
  112      CONTINUE
C&&& LEVYLE : Use Levy-Leblond Hamiltonian
             LEVYLE = .TRUE.
             DOLVC  = .FALSE.
C            Levy-Leblond implies also spinfree, specify that
             SPINFR = .TRUE.
             NOSMLV = .TRUE.
             SSMTRC = D0
C            Check that we have not defined operators yet, as they need to be modified in Levy-Leblond calcs.
             IF (I_FINFIELD.GT.0)
     &          CALL QUIT ('give .LEVY-LEBLOND option before .OPERAT')
C            No modifications of Hamiltonian, only new definition
C            ICHANG = ICHANG + 1
             GO TO 100
  113      CONTINUE
C&&& SPINFREE : Zero out spin-orbit terms
             SPINFR = .TRUE.
             ICHANG = ICHANG + 1
             GOTO 100
  114      CONTINUE
C&&& ZORA   : Use ZORA approximation
             READ(LUCMD,*) IZORA4,IZORAS
             ZORA   = .TRUE.
             ZORA4  = (IZORA4.NE.0)
             ZORASC = (IZORAS.NE.0)
             DOLVC  = .FALSE.
             SSMTRC = D0
             ICHANG = ICHANG + 1
             GO TO 100
  115      CONTINUE
C&&& SMLV1C : Neglect potential for off-diagonal SS-blocks
C    jkp    : This keyword enables Luuk's one-center
C             projection scheme
             SMLV1C = .TRUE.
             ONECAP = .TRUE.
             DOLVC  = .FALSE.
             INTV1C = 1
             ICTLV1C(2) = 1
             ICHANG = ICHANG + 1
             GO TO 100

!          old ALDA keyword
  116      continue
           write(lupri, *) 'place .ALDA keyword after *DFT'
           call quit('place .ALDA keyword after *DFT')
           go to 100

!          .DFT
  117      CONTINUE
C&&& DFT    : Kohn-Sham (DFT) Hamiltonian
!radovan: line can be longer than 80
               read(lucmd, '(a100)') line100
               IF (.NOT. SETFUN) THEN
C              ... hjaaj: DFTINPUT accumulates def. of functional
C                  so if we call more than once (as in geom.opt.)
C                  we get n*functional for n'th call !!!! /3-Jun-2003 hjaaj
                  CALL DFTINPUT(line100, INPERR)
                  if (inperr /= 0) then
                     call quit('failure in DFTINPUT after .DFT')
                  end if
                  call set_hf_exchange_factor(hfxfac)
               END IF
               dirac_cfg_dft_calculation = .true.
               SETFUN = .TRUE.
            GO TO 100

  118     CONTINUE
C&&& HFXFAC: weight of pure exchange
             READ(LUCMD,*) HFXFAC
            GO TO 100
  119     CONTINUE
C&&& SCQSET: Reset small component charge for nuclei IZ (for onecap)
             READ(LUCMD,*) IZ, SCQVAL
             CALL SCQSET(IZ,SCQVAL)
            GO TO 100
  120     CONTINUE
C&&& ONECNV: Set threshold for turning off ONECAP model
            READ(LUCMD,*) ONECNV
            GO TO 100
  121     CONTINUE
C&&& ONECAP : Approximate potential for off-diagonal SS-blocks
C             if no model specified, goto 199, choose onecap 2
            READ(LUCMD,'(I10)',END=199,ERR=199) JNTV1C
            INTV1C = MOD(JNTV1C,10)
            ICTLV1C(2) = JNTV1C / 10
            GOTO 198
  199       INTV1C = 2
            ICTLV1C(2) = 0
            BACKSPACE LUCMD
  198       CONTINUE
            ONECAP  = .TRUE.
            DOLVC   = .FALSE.
            IF(INTV1C.EQ.2) DOLVC = .TRUE.
C
C   ... Whether to evaluate integrals in a direct or concentional
C       matter has not been specified in **GENEREL
C       We write LS and SS integrals to disk by default in the
C       ONECAP type 2 & 3 model.            /jkp
C
C       Switched off in this version [DIRSET = .TRUE.]
C       - even one-center integrals become quite numerous in some cases!!!
            IF (.NOT.DIRSET .AND. INTV1C.NE.2) THEN
               ILLDIR = 1
               ISLDIR = 0
               ISSDIR = 0
               IDFLAG = ILLDIR+2*ISLDIR+4*ISSDIR+8*IGTDIR
            END IF
            ICHANG = ICHANG + 1
            GO TO 100
  122     CONTINUE
C&&& SPINFR2 : debug: Zero out 2-electron spin-orbit terms
            SPINFR2 = .TRUE.
            ICHANG = ICHANG + 1
            GO TO 100
  123     CONTINUE
C&&& NOSPIN: Take out all spin-interactions
             SPINFR = .TRUE.
             NOSPIN = .TRUE.
             ICHANG = ICHANG + 1
            GO TO 100
  124     CONTINUE
C&&  USEDF: Do the block-diagonalization of the converged Fock-Dirac operator
C     (and after that continue in two-component mode with non-canonical orbitals,
C       or continue SCF proces to get cannonical orbitals.)
            USEDF = .TRUE.
            GO TO 100
  125     CONTINUE
C&&& INTFLG: Specify what two-integrals should be included in this run
            IF (IGTINT .EQ. 1) THEN
               READ(LUCMD,*,IOSTAT=IOS) ILLINT,ISLINT,ISSINT,IGTINT
               IF (IOS.NE.0) THEN
                CALL QUIT(
     &   '**HAMILTONIAN: Error in reading .INTFLG with Gaunt')
               ENDIF
            ELSE
               READ(LUCMD,*,IOSTAT=IOS) ILLINT,ISLINT,ISSINT
               IF (IOS.NE.0) THEN
                CALL QUIT(
     &   '**HAMILTONIAN: Error in reading INTFLG without Gaunt')
               ENDIF
               INTFLG_CHANGE = .TRUE.
            END IF
            GO TO 100
  126     CONTINUE
C&&&&&  BSS:  prepare two-component Barysz-Sadlej-Snijders
C       transformed relativistic Hamiltonian
             BSS        = .TRUE.
             TWOCOMPBSS = .TRUE.
             TWOCOMP    = .TRUE.
             DOLVC   = .FALSE.
C        ... read 3-digit number
             READ(LUCMD,*,IOSTAT=IOS) IBSS
             IF (IOS.NE.0) THEN
               WRITE(LUPRI,'(2X,A)')
     &         '**HAMILTONIAN: Error in reading .BSS/IBSS'
               WRITE(LUPRI,*) 'IBSS=',IBSS
               CALL QUIT(
     &         '**HAMILTONIAN: Error in reading of .BSS/IBSS')
             ENDIF
C
C            Check that we have not defined operators yet,
C            as they need to be modified within BSS transformations
             IF (I_FINFIELD.GT.0) THEN
               CALL QUIT(
     &        '**HAMILTONIAN:give .BSS option before .OPERAT !')
             ENDIF
            GO TO 100

C&&& HFXMU: range separation
  127     CONTINUE
          READ(LUCMD,*) HFXMU
          GO TO 100

!C&&&&&  AMFI_CHARGE: Get the total charge of the molecule
!C     for mean-field summations of each atom in the system
!C          ... read integer number
!             READ(LUCMD,*) IMFCH
            GO TO 100
  128     CONTINUE
C&&&&   Solve R-equation 1
            YREQ1 = .TRUE.
            GO TO 100
  129     CONTINUE
C&&&&  ONESTEP: Do NOT do preliminary free-particle transformation
            NOPRTR = .TRUE.
            GO TO 100
  130     CONTINUE
C&&& BLOCKD:   Do use block diagonalization
CMI      ... not active yet
#ifdef MOD_UNRELEASED
             BLOCKD = .TRUE.
#else
             call quit('block diagonalizing feature not in release!')
#endif
            GO TO 100
  131     CONTINUE
C&&& DO2C4C: jump to 4c level after few iterations in the 2-component quasirel.mode...
            DO2C4C = .TRUE.
            START2C = .TRUE.
            INI2C   = 3
            call quit('2c -> 4c SCF functionality is broken -'//
     &                ' please wait for a patch to Dirac14 or contact'//
     &                ' dirac-users@googlegroups.com for further'//
     &                ' information.')
            !TWOCOMP=.FALSE.
            !TWOCOMPBSS=.TRUE.
            GO TO 100
  132     CONTINUE
C&&& DO4C2C:  descend to the 2c level after 4c SCF ...
            DO4C2C = .TRUE.
            GO TO 100
  133     CONTINUE
C&&& CONT2C:  after 4c DC-SCF continue SCF with BSS using picture-change transformed Fock-Dirac operator !
            CONT2C = .TRUE.
C          ... read the specification of 2c 1el Hamiltonian (3,4,5)
            READ(LUCMD,*) I2CHAM
            GO TO 100
  134     CONTINUE
C&&& BEG_2C: specify INI2C when DO2C4C is on
            START2C = .TRUE.
            DO2C4C = .TRUE.
            READ(LUCMD,*) INI2C
            call quit('2c -> 4c SCF functionality is broken -'//
     &                ' please wait for a patch to Dirac14 or contact'//
     &                ' dirac-users@googlegroups.com for further'//
     &                ' information.')
            GO TO 100
  135     CONTINUE
#ifdef MOD_UNRELEASED
C&&&& VACPOL: vacuum polarization included
            VACPOL = .TRUE.
#endif
            GO TO 100
  136     CONTINUE
C&&&& X2COLD: flag for the one-step Exact TWO-component relativistic Hamiltonian
!                                   ^     ^^^ ^
!             obtained by Miro's old BSS/X2C module
            IOTC       = .TRUE.
            BSS        = .TRUE.
            NOPRTR     = .TRUE.
            IBSS       = 2999
            TWOCOMPBSS = .TRUE.
            TWOCOMP    = .TRUE.
            DOLVC      = .FALSE.
C           Check that we have not defined operators yet,
C           as they need to be modified within BSS transformations
            IF (I_FINFIELD.GT.0) THEN
              CALL QUIT('**HAMILTONIAN:'//
     &         ' specify the .X2COLD option before .OPERATOR !')
            ENDIF
            GO TO 100
  137     CONTINUE
C&&&  X2C: one-step Exact TWO-component relativistic Hamiltonian
!                   ^     ^^^ ^
             NOPRTR     = .true.
             x2c        = .true.
             x2cmod_x2c = .true.
             DOLVC      = .FALSE.
C            Check that we have not defined operators yet,
C            as they need to be picture-changed transformed after the 2c-transformation
             IF (I_FINFIELD.GT.0) THEN
               CALL QUIT(
     &        '**HAMILTONIAN: specify the .X2C option before .OPERAT !')
             ENDIF
             GO TO 100
  138     CONTINUE
#ifdef MOD_UNRELEASED
C&&& MIXVAC: Alternative definition of vacuum polarisation
             VACPOL = .TRUE.
             MIXVAC = .TRUE.
#endif
             GO TO 100
  139     CONTINUE
C&&& Activate the CAP calculation for the complex Fock space method
#ifdef MOD_CAP
             CAP  =  .TRUE.
#else
             CAP  =  .FALSE.
#endif
             GO TO 100
  140     CONTINUE
C&&& NOVREF: Do not generate reference vacuum
             VACREF = .FALSE.
             GO TO 100
  141     CONTINUE
#ifdef MOD_UNRELEASED
C&&& POSVAC: Use positronic vacuum polarisation
             VACPOL = .TRUE.
             POSVAC = .TRUE.
#endif
             GO TO 100

!         .FDE
!         frozen density embedding
  142     continue
          dirac_cfg_fde = .true.
          go to 100

  143     CONTINUE
C&&& IOTC4: do IOTC in 4c-mode
!radovan: note that .X2C4 points to here
!         so be careful when changing "143"
             BSS = .TRUE.
             NOPRTR = .TRUE.
             IBSS = 2999
             DOLVC      = .FALSE.
             GO TO 100
  144     CONTINUE
C&&& NONREL: do 2-component non-relativistic calculation
             TWOCOMP = .TRUE.
             NONREL = .TRUE.
             SPINFR = .TRUE.
             DOLVC   = .FALSE.
             GO TO 100
  145     CONTINUE
             BSS = .TRUE.
             DOLVC = .FALSE.
C        ... read 3-digit number
             READ(LUCMD,*) IBSS
C
C    Check that we have not defined operators yet,
C      as they need to be modified within BSS transformations
             IF (I_FINFIELD.GT.0) THEN
               CALL QUIT(
     &        '**HAMILTONIAN:give .BSS option before .OPERAT !')
             ENDIF
C&&& BSS4: do BSS in 4c-mode
             GO TO 100
  146     CONTINUE
C&&& CMPEIG: add comparison between eigenvalues of 2c and 4c Hamiltonian
             CMPEIG = .TRUE.
             GO TO 100
  147     CONTINUE
C&&& DKH2: complex flag for second order Douglas-Kroll-Hess with AMFI SSO
             DKH2 = .TRUE.
             BSS = .TRUE.
             NOPRTR = .FALSE.
             IBSS = 2022
             TWOCOMPBSS = .TRUE.
             TWOCOMP = .TRUE.
             DOLVC = .FALSE.
             GO TO 100
  148     CONTINUE
C&&& DKH1: first order Douglas-Kroll-Hess with AMFI SSO, set other
             DKH1 = .TRUE.
             BSS = .TRUE.
             NOPRTR = .FALSE.
             IBSS = 2011
             TWOCOMPBSS = .TRUE.
             TWOCOMP = .TRUE.
             DOLVC = .FALSE.
             GO TO 100
  149     CONTINUE
C&&& .NOAMFI: do not include AMFI contribution to the 2comp. Hamiltonian
CMI    when the AMFI-module is selected through the preprocessor command
             NOAMFI = .TRUE.
             GO TO 100
  150     CONTINUE
C&&& QDOTS: activate artificial atom (quantum dots)
#ifdef MOD_QDOTS
             QDOTS = .TRUE.
#endif
             GO TO 100

!&&& IOTC: complex flag for the one-step infinite-order two-component relativistic Hamiltonian
CMI  it's BSS_RKB+AMFI; later fix this options when interacting
CMI  with other flags...
  151     continue
             IOTC = .TRUE.
             BSS = .TRUE.
             NOPRTR = .TRUE.
             IBSS = 2999
             TWOCOMPBSS = .TRUE.
             TWOCOMP = .TRUE.
             DOLVC = .FALSE.
             GO TO 100

!         .X2C4 is alias for .IOTC4
  152     continue
          go to 143
  153     CONTINUE
C&&& NOSFMU: do not generate spinfree multiplication table in spinfree calculations
C    (see MRCON1 of the MOLTRA module)
          NOSFMU = .TRUE.
          GOTO 100

!         .DFTAUTO (dft with automatic differentiation)
  154     continue
          read(lucmd, '(a80)') line
          call parse_functional(line, xc_fun,
     &                          hfx_out, mu_out, beta_out, .true.)
          hfxfac = hfx_out
          hfxmu  = mu_out
          hfxatt = beta_out
          call set_xc_fun_alda(xc_fun_alda)
          call set_xc_fun_xalda(xc_fun_xalda)
          fun_is_automatic = .true.
          dirac_cfg_dft_calculation = .true.
          go to 100

  155     CONTINUE
C&&& Activate : molecular mean field approximation for correlated calculations
C               using LL integrals only that is ^4 DC(G)**.
          MOLMF = .TRUE.
          GOTO 100
  156     CONTINUE
CMI&& Activate 4c-2c MO transformation
          TRMO4C2C=.TRUE.
          GOTO 100
! .PCM
  157     CONTINUE
#ifdef HAS_PCMSOLVER
          dirac_cfg_pcm = .true.
#else
          call quit(
     &     '**HAMILTONIAN: PCMSolver not compiled, no PCM available!')
#endif
          GOTO 100
! .HFXATT
  158     CONTINUE
          READ(LUCMD,*) HFXATT
          GOTO 100
! .xB
  159     CONTINUE
          call quit('.xB not implemented in this version')
          GOTO 100
! .xC
  160     CONTINUE
          call quit('.xC not implemented in this version')
          GOTO 100

! .ECP
  161     CONTINUE
          ECPCALC = .TRUE.
          TWOCOMP = .TRUE.
          DOLVC   = .FALSE.
          GO TO 100

  162     CONTINUE
!         X2Cmmf: molecular mean-field approach within the X2C framework
          noprtr                = .true.
          x2c                   = .true.
          x2cmod_x2c            = .true.
          x2cmod_mmf            = .true.
          x2c_2c_mmf_mos        = .true.
          x2c_4c_fock_mtx_defh1 = .true.
          x2c_add_amfi          = -1
          X2CMMF                = .TRUE.
C         Check that we have not defined operators yet,
C         as they need to be picture-changed transformed after the 2c-transformation
          IF(I_FINFIELD .gt. 0)THEN
            CALL QUIT(
     &     '**HAMILTONIAN: specify the .X2Cmmf option before .OPERAT !')
          ENDIF
          GO TO 100
  163     CONTINUE
C&&&& DOSSSS: Including SS-integral contribution explicitly
            DOLVC = .FALSE.
            GO TO 100
  164     CONTINUE
C&&& MDIRAC: use modified Dirac Hamiltonian as implemented in INTEREST integral package
            CALL QUIT('.MDIRAC option not implemented => STOP!')
            GO TO 100
C..TODO: check for inconsistencies and add stops
  165     CONTINUE
!&&&& .TRSMIX : allow mixing T+ Hamiltonian with added T- perturbation operator
            TRS_MIX = .TRUE.
          GO TO 100
!&&&& .PEQM: polarizable embedding
  166     CONTINUE
#ifdef HAS_PELIB
           IF(NBSYM.GT.1) THEN
             CALL QUIT('Polarized embedding '//
     &       'does not work with symmetry yet.')
            ENDIF
            WRITE(LUPRI,*) 'PE library will be used.'
            PEQM = .TRUE.
#else
            WRITE(LUPRI,*) 'ERROR for .PEQM: PE library not enabled.'
            CALL QUIT('ERROR for .PEQM: PE library not enabled.')
#endif
            GO TO 100
  167     CONTINUE
!&&&& .JZOUT : save Jz_MO matrices into formatted files for analysis
            WRITE_Jz_MATRIX = .TRUE.
          GO TO 100
  168     CONTINUE
!&&&& .FAKE2C : do a 2C calculation without generating a 2C
!               1-electron Hamilonian (can be used for correlated
!                energy calculations)
            FAKE2C  = .TRUE.
            TWOCOMP = .TRUE.
            NONREL  = .TRUE.
            SPINFR  = .TRUE.
            DOLVC   = .FALSE.
          GO TO 100
  169     CONTINUE
!&&&& GAUGEO
          READ (LUCMD, *,IOSTAT=IOS) (GAGORG(I),I=1,3)
          GAGORG_SET = .TRUE.
          IF (IOS.NE.0) THEN
            CALL QUIT('Error in reading (GAGORG(I),I=1,3)!')
          ENDIF
          GO TO 100
  170     CONTINUE
!&&&& GO ANG
          GO TO 100
!         read user specified gauge origin in angstrom
          read(lucmd, *, iostat=ios) gagorg
          GAGORG_SET = .TRUE.
          if (ios /= 0) then
            call quit('error in reading .GO ANG')
          end if
!         convert from angstrom to bohr
          gagorg = gagorg/xtang
  171     CONTINUE
!&&&& DIPORG
          READ (LUCMD, *,IOSTAT=IOS) DIPORG(1:3)
          IF (IOS.NE.0) THEN
            CALL QUIT('Error in reading DIPORG(1:3)')
          ENDIF
          GO TO 100
  172     CONTINUE
!&&&& PHASEO
          READ (LUCMD, *,IOSTAT=IOS) (ORIGIN(I),I=1,3)
          IF (IOS.NE.0) THEN
            WRITE(LUPRI,*) 'Error in reading ORIGIN(1:3) !'
            WRITE(LUPRI,*) 'ORIGIN(1:3) :', (ORIGIN(I),I=1,3)
            CALL QUIT('Error in reading ORIGIN(1:3) !')
          ENDIF
  173     CONTINUE
!&&&& .KOUT : save Kappa_MO matrices into formatted files for analysis
            WRITE_K_MATRIX = .TRUE.
          GO TO 100
  174     CONTINUE
          GO TO 100
  175     CONTINUE
          GO TO 100
  176     CONTINUE
          GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
        GO TO 300
      ELSE
        WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *         '" not recognized for **HAMILTONIAN.'
        CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
        CALL QUIT('Illegal prompt for **HAMILTONIAN.')
      END IF
  300 CONTINUE

      WRITE (LUPRI,'(/A/,(A,3F20.15))')
     &    ' One-electron operator origins:',
     &    ' - General operator origin (a.u.)       :', ORIGIN(1:3),
     &    ' - Magnetic gauge origin (a.u.)         :', GAGORG(1:3),
     &    ' - Dipole (and multipole) origin (a.u.) :', DIPORG(1:3)

      if (dirac_cfg_fde) then
         call fde_input_init(lupri,'DIRAC.OUT')
      end if 
      GAUNT_SAVE = GAUNT

!     take care of the correct choice of the X2C module
!     -------------------------------------------------
!     1. X2C+spinfree
!     2. X2C+compare eigenvalues 4c-Hamiltonian <--> 2c-Hamiltonian
      if(x2c.and.spinfr)then
        isorder_amfi_x2c = -1
        x2c_add_amfi     = -1
!     X2C+compare eigenvalues 4c-Hamiltonian <--> 2c-Hamiltonian
      else if(x2c.and.cmpeig)then
        write(lupri,'(a/)')
     &  ' *** Warning: the .cmpeig option is not implemented in'//
     &  ' the new X2C module. please rerun with .X2COLD ***'
        cmpeig = .false.
      end if
C
C     Remove inconsistencies
C     =====================
C
      IF (NOSMLV .AND. ONECAP) THEN
        IF(INTV1C.EQ.1.AND.ICTLV1C(2).EQ.1) THEN
          WRITE(LUPRI,'(A)')
     &  ' INFO: .SMLV1C ignored because .LEVY-L or .NOSMLV specified'
        ELSE
          WRITE(LUPRI,'(A)')
     &  ' INFO: .ONECAP ignored because .LEVY-L or .NOSMLV specified'
        ENDIF
        ONECAP = .FALSE.
      END IF
      IF (SMLV1C.AND.ONECAP.AND.(INTV1C.NE.1).AND.(ICTLV1C(2).NE.1))THEN
C     ... User asked for both SMLV1C and a non-equivalent ONECAP model
C         The ONECAP model overrides SMLV1C !!!
          WRITE(LUPRI,'(A)') ' INFO: .SMLV1C ignored - '//
     &   'another ONECAP model was specified.'
      ENDIF
      IF (LEVYLE) THEN
         IF(FREEPJ) CALL QUIT('Conflicting options : LEVYLE & FREEPJ !')
         IF(VEXTPJ) CALL QUIT('Conflicting options : LEVYLE & VEXTPJ !')
      ENDIF
C
      IF( (FREEPJ.OR.VEXTPJ.OR.LEVYLE) .AND. URKBAL )
     $   CALL QUIT('FREEPJ, VEXTPJ, or LEVYLE require URKBAL = .FALSE.')
C
      IF(FREEPJ.AND.VEXTPJ) THEN
        WRITE(LUPRI,'(A)') '***** **HAMILTONIAN ERROR *****'
        WRITE(LUPRI,'(A)') 'You cannot specify both VEXTPJ and FREEPJ !'
        CALL QUIT('Conflicting options : VEXTPJ and FREEPJ !')
      ENDIF

      IF (DO2C4C.AND.DO4C2C) THEN
        WRITE(LUPRI,'(A)') '.DO2C4C and .DO4C2C not accepted !'
        CALL QUIT('Conflicting options : .DO2C4C and .DO4C2C !')
      ENDIF

!     x2c module: check for the keyword indicating the 4c-DF operator as defining h1
      if(find_keyword('.X2Cmmf'))then
        do4c2c = .true.
!       default for the x2c module: if 4c-fock --> 2c then use the 4c-fock operator as defining h1
        usedf  = .true.
      end if
      if(find_keyword('**HAMIL')) call move_to_next_star(word, lucmd)

CMI   ... take care of 2c-4c/4c-2c transitions !
      IF (DO2C4C.OR.DO4C2C) THEN
         TWOCOMP    = TWOCOMP_SAVE
         TWOCOMPBSS = TWOCOMPBSS_SAVE
      ENDIF

!     quit if user wants .LEVY-LEBLOND and .LVCORR
      if (levyle .and. dolvc) call quit
     &   ('.LEVY-LEBLOND and .LVCORR should not be used together')


CMI     Remove inconsistencies for the TWOCOMP option
      IF (TWOCOMP.AND.(LEVYLE.OR.
     &              URKBAL.OR.FREEPJ.OR.VEXTPJ)) THEN
        WRITE(LUPRI,'(A)') '***** **HAMILTONIAN ERROR *****'
        WRITE(LUPRI,'(A)')
     &  'You cannot specify both TWOCOMP and other otions '//
     &  ' (LEVYLE.OR.URKBAL.OR.FREEPJ.OR.VEXTPJ)!'
        CALL QUIT(
     &  'Conflicting options : TWOCOMP and others'//
     &  ' (see the code here) !')
      ENDIF
C
C     Remove inconsistencies for the BSS/X2C option
C     =====================================================
      IF ((DO2C4C.OR.DO4C2C.OR.USEDF).AND..NOT.(BSS.or.x2c)) THEN
          CALL QUIT('**HAMILTONIAN: '//
     &     'DO2C4C/DO4C2C/USEDF options require BSS keyword !')
      ENDIF

      IF (BSS .or. x2c) THEN
        IF (LEVYLE.OR.FREEPJ.OR.VEXTPJ) THEN
          CALL QUIT('**HAMILTONIAN: '//
     &     'Conflicting options LEVYL/FREEPJ/VEXTPJ & BSS/X2C !')
        ENDIF
        IF ((SMLV1C.OR.DOLVC.OR.NOSMLV.OR.ONECAP).and.
     &      .not.x2cmod_mmf)THEN
          CALL QUIT('**HAMILTONIAN: '//
     &    'Conflicting options SMLV1C/DOLVC/NOSMLV/ONECAP
     &    & BSS/X2C !')
        ENDIF
        IF (URKBAL.OR.SPINFR2) THEN
          CALL QUIT('**HAMILTONIAN: '//
     &     'Conflicting options URKBAL/SPINFR2 & BSS/X2C!')
        ENDIF
C ...   DIRCON(1) means contracted L component basis set !
        IF (ZORA) THEN
          CALL QUIT('**HAMILTONIAN: '//
     &     'Conflicting options ZORA & BSS/X2C !')
        ENDIF

CMI  ... check not combining .IOTC, .DKH1, .DKH2
        IF ((IOTC.AND.DKH1).OR.(IOTC.AND.DKH2).OR.(DKH2.AND.DKH1)) THEN
          WRITE(LUPRI,'(/3(A,L6))')
     &    '**HAMILTONIAN: conflicting options; IOTC=',IOTC,
     &    ' DKH1=',DKH1,' DKH2=',DKH2
          CALL QUIT('**HAMILTONIAN: Conflicting options IOTC & DKH1 '//
     &    'or DKH2 & DKH2 or IOTC & DKH1 !')
        ENDIF

CMI  ... todo: more precise specification needed
        IF (.NOT.(DO4C2C.AND.CONT2C)) CONT2C=.FALSE.

CMI ... set up the default value of I2CHAM, when the BSS-postDC-SCF is on
        IF (DO4C2C.AND.I2CHAM.EQ.-1.AND..NOT.CONT2C) THEN
          if(.not.x2c)then
            IF (USEDF) THEN
              I2CHAM = 3 ! use Fock-Dirac matrix for 2c transformation
            ELSE
              I2CHAM = 4 ! use Dirac bare nucleus for 2c transformation
            ENDIF
          end if
        ENDIF

        IF (CONT2C.and.BSS) THEN
         IF (I2CHAM.NE.3.AND.I2CHAM.NE.4.AND.
     &       I2CHAM.NE.5.AND.I2CHAM.NE.6) THEN
           write(lupri,*) '**HAMILTONIAN: wrong value of I2CHAM=',I2CHAM
          CALL QUIT('**HAMILTONIAN: wrong value of I2CHAM !')
         ENDIF
        ENDIF

       IF (DOPRP) THEN
        WRITE(LUPRI,'(2X,A)') '- BSS with properties !'
       ENDIF

       IF(DO2C4C)THEN
         IF(GAUNT)THEN ! check whether it has been included already in INTGEN_SAVE
           IF(IBTAND(INTGEN_SAVE/8,1) /= 1) INTGEN_SAVE = INTGEN_SAVE+8
         END IF
       END IF

CMI ... fix the .IOTC keyword interaction with .GAUNT,.ONESYS,.SPINFR
        IF (((IBSS.EQ.2999).or.x2c).AND.NOPRTR) THEN
         IF (SPINFR) THEN
           IBSS = 009
C          write(lupri,*) '.IOTC for .SPINFR'
         ENDIF
         IF (GAUNT) THEN
CMI        ... AMFI SSO+SOO terms
CMI        ... we have to switch off this term as we are in the 2c framework
CMI        ... instead we use GAUNT_SAVE
           if(.not.x2c)then
             IBSS = 3999
           else
!            X2C and AMFI - handling of GAUNT; the order can also be set through
!            the *AMFI input deck.
             isorder_amfi_x2c = 3
           end if
           if(.not.x2cmod_mmf)then
             GAUNT = .FALSE.
           end if
         ENDIF
CMI      IF (ONESYS) THEN
CMI        ... no AMFI contribution
CMI        IBSS = 999
CMI      ENDIF
        ENDIF


CMI   ... fix the DKH2 keyword interaction with .GAUNT,.ONESYS,.SPINFR
        IF (DKH2) THEN
          IF (SPINFR) THEN
            IBSS = 002
          ENDIF
          IF (GAUNT) THEN
CMI        ... AMFI SSO+SOO terms
CMI        ... we have to switch off this term as we are in the 2c framework
CMI        ... instead we use GAUNT_SAVE
           GAUNT = .FALSE.
           IBSS = 3022
          ENDIF
        ENDIF ! end of DKH2

CMI   ... fix the DKH1 keyword interaction with .GAUNT,.ONESYS,.SPINFR
        IF (DKH1) THEN
          IF (SPINFR) THEN
            IBSS = 001
          ENDIF
          IF (GAUNT) THEN
CMI        ... AMFI SSO+SOO terms
CMI        ... we have to switch off this term as we are in the 2c framework
CMI        ... instead we use GAUNT_SAVE
           GAUNT = .FALSE.
           IBSS = 3011
          ENDIF
        ENDIF ! end of DKH1

      ENDIF ! End of .BSS
C
C     If Gaunt is set we always want to include these integrals
C
      IF (GAUNT) THEN
        IGTINT = 1
      ENDIF
C
C     NZC1 is the NZ value in C1 symmetry
C
      IF (NONREL .AND. .NOT. FAKE2C) THEN
       NZC1 = 1
      ELSE
       NZC1 = 4
      END IF
C
C.....Check RKBIMP
      IF (DORKBIMP) THEN
        IF (.NOT.URKBAL) THEN
          DORKBIMP = .FALSE.
          WRITE(LUPRI,'(1X,A)')
     &      '*** WARNING! RKBIMPort ignored since no URKBAL'
        ENDIF
      ENDIF

C
C     Print section
C     =============
C
 400  CONTINUE
      WRITE(LUPRI,'(A,I5)') ' * Print level:',IPRHAM
C
      IF(LEVYLE) THEN
        WRITE(LUPRI,'(A)') ' * Levy-Leblond Hamiltonian'
      ELSE IF (BSS.or.x2c) THEN
        IF(NOPRTR) THEN
          if (bss) then
            WRITE(LUPRI,'(A)')
     & ' * One-step Infinite-Order Two-component Hamiltonian'
          else if(x2c)then
            WRITE(LUPRI,'(A)')
     & ' * Exact-Two-Component (X2C) Hamiltonian'
          end if
          WRITE(LUPRI,'(A)')
     & '   Reference: ',
     & '    M. Ilias and T. Saue:',
     & '    "Implementation of an infinite-order two-component '//
     & 'relativistic Hamiltonian ',
     & '    by a simple one-step transformation." ',
     & '    J. Chem. Phys., 126 (2007) 064102.'
          if(x2c)then
          WRITE(LUPRI,'(a)')
     & '   additional reference for the new X2C module:'
          WRITE(LUPRI,'(a)')
     & '    S. Knecht and T. Saue:'
          WRITE(LUPRI,'(a/)')
     & '    manuscript in preparation, Strasbourg 2010.'
          end if
        ELSE
C         .... two-step method of Barysz
          WRITE(LUPRI,'(A,I5)')
     & ' * Barysz-Sadlej-Snijders Hamiltonian. Type code:',IBSS
          WRITE(LUPRI,'(A)')
     & '   Reference: ',
     & '    H. J. Aa. Jensen and M. Ilias,',
     & '    "Two-component relativistic methods based on the '//
     & 'quaternion modified Dirac equation. I:',
     & '     from Douglas-Kroll-Hess second order method to the '//
     &   'infinite order two-component method.",',
     & '     J. Chem. Phys., in preparation.'
          IF (DKH2) THEN
           WRITE(LUPRI,'(1X,A)')
     &     '*  second-order Douglas-Kroll-Hess method, .DKH2'
          ENDIF
          IF (DKH1) THEN
           WRITE(LUPRI,'(1X,A)')
     &     '* first-order Douglas-Kroll-Hess, .DKH1, very poor method !'
          ENDIF
        ENDIF

        IF (DO4C2C) THEN
         IF (CONT2C) THEN
           WRITE(LUPRI,'(A,/A)')
     & ' * after the four-component SCF procedure'//
     & ' perform a X2C/BSS/IOTC transformation',
     & '   of the 4c Fock operator'//
     & ' and continue the SCF step at the two-component level.'
         ELSE
           WRITE(LUPRI,'(A,/A)')
     & ' * after the four-component SCF procedure'//
     & ' perform a X2C/BSS/IOTC transformation',
     & '   of the 4c Fock operator.'
         ENDIF
         IF (TRMO4C2C) THEN
           WRITE(LUPRI,'(1X,A)')
     &     '* transforming 4c MO to 2c form, no extra 2c SCF iteration.'
         END IF
        ENDIF

        IF (DO2C4C) THEN
         WRITE(LUPRI,'(1X,A,I1)')
     &   '* after (few) two-component SCF iterations'//
     &   ' switch to the four-component level ! Code INI2C=',INI2C
        ENDIF

        IF (TWOCOMP) THEN
          WRITE(LUPRI,'(A)') ' * Running in two-component mode'
        ELSE
          WRITE(LUPRI,'(A)') ' * Running in four-component mode'
        ENDIF

C       ... take care of AMFI contribution
        IF (BSS.AND.NOAMFI) THEN
         WRITE(LUPRI,'(1X,A)')  '  - no AMFI contribution to the '
     &   //'2-component transformed Hamiltonian !'
        ENDIF
        IF (IBSS.GT.999.AND..NOT.NOAMFI) THEN
          WRITE(LUPRI,'(A)')
     &     ' * Spin-orbit corrections to the instantaneous Coulomb ',
     &     '   interaction provided by the AMFI code. ',
     &     '   References: ',
     &     '     1) B. Schimmelpfennig, program AMFI, ',
     &     '        Stockholm University, Sweden',
     &     '     2) B. A. Hess, C. M. Marian, U. Wahlgren'//
     &     ' and O. Gropen',
     &     '        Chem. Phys. Lett. 251, 365 (1996)',
     &     '     3)  M Ilias, V Kello, L Visscher'//
     &     ' and B Schimmelpfennig',
     &     '        J. Chem. Phys. Vol. 115 (2001) p.9667.',
     &     '     Adding :',
     &     '       some AMFI contribution'//
     &     ' (see output from EXTR_BSS_INFO below)'
          IF (GAUNT_SAVE) THEN
            WRITE(LUPRI,'(A)')
     &     '     N.B: Gaunt key is on'
          ENDIF

CMI       ... take care of AMFI CSC test transformation, if specified through the entering IBSS number
          IF (IBSS.GT.999.AND.IBSS.LE.1999) THEN
           WRITE(LUPRI,'(1X,A)')
     &     '* only C->S->C transformation of 2comp.'//
     &     'Hamiltonian; AMFI integrals calculated, but NOT added !'
          ENDIF

        ENDIF

      ELSE IF (ECPCALC) THEN
        WRITE(LUPRI,'(A)')
     &   ' * Employing effective-core-potentials (ECP)'
      ELSE IF (NONREL) THEN
        WRITE(LUPRI,'(A)') ' * Non-relativistic Hamiltonian'
      ELSE IF(MDIRAC) THEN
        WRITE(LUPRI,'(A)') ' * Modified Dirac-Coulomb Hamiltonian'
        WRITE(LUPRI,'(A)') 
     &    '   Switching to the INTEREST integral package (M. Repisky)'
      ELSE IF (GAUNT) THEN
        WRITE(LUPRI,'(A)') ' * Dirac-Coulomb-Gaunt Hamiltonian'
      ELSE
        WRITE(LUPRI,'(A)') ' * Dirac-Coulomb Hamiltonian'
      ENDIF

      if (x2c.and.ecpcalc) then
         call quit('Incompatible options: .X2C* and .ECP')
      end if
      if (GAUNT.and.ecpcalc) then
         call quit('Incompatible options: .GAUNT and .ECP')
      end if

      IF (ICHANG.GT.0) WRITE(LUPRI,'(3X,A)')
     &    'with the following modifications:'

#ifdef MOD_UNRELEASED
      IF (VACPOL) THEN
        WRITE(LUPRI,'(3X,A)') 'Vacuum polarization added.'
        IF(MIXVAC) THEN
          WRITE(LUPRI,'(6X,A)')
     &      '- alternative definition of vacuum polariation.'
        ENDIF
        IF(POSVAC) THEN
          WRITE(LUPRI,'(6X,A)')
     &      '- using positronic vacuum.'
        ENDIF
        IF(.NOT.VACREF) THEN
          WRITE(LUPRI,'(6X,A)')
     &      '- no reference vacuum generated.'
        ENDIF
      ENDIF
#endif
#ifdef MOD_CAP
CMI ...  CAP processing
      IF (CAP) THEN
        WRITE(LUPRI,'(1X,A)')
     &  '* Activated calculation of the Complex'//
     &  ' Absorption Potential (CAP) integrals'
      ENDIF
#endif

#ifdef MOD_QDOTS
CMI    ... quantum dots
      IF (QDOTS) THEN
        WRITE(LUPRI,'(2X,A)')
     &  '* Activated calculations of artificial system (quantum dots).'
      ENDIF
#endif

      IF(SPINFR) THEN
         WRITE(LUPRI,'(3X,A)') '- Spin-orbit interactions neglected'
         IF(NOSFMU) THEN
           WRITE(LUPRI,'(6X,A,A)')
     &      'Spinfree multiplication tables (spin x space) NOT ',
     &      'generated at the correlated level.'
         ENDIF
      ELSE IF(SPINFR2) THEN
         WRITE(LUPRI,'(3X,A)') '- 2-el. spin-orbit terms zeroed out '//
     &      '(NB : ONLY in the DHF calculation !!)'
      END IF

      if (dirac_cfg_dft_calculation) then

        if (.not. fun_is_automatic) then
!         this will report chosen functional
!         if nothing is chosen it will report lda
!         and if .dftauto functional is chosen
!         it would report lda
!         so let's not call it in this case
          call dftreport()
        end if

         IF (IGTINT.EQ.1) THEN
            IF (LSCALE_DFT_GAUNT) THEN
               WRITE (LUPRI,*) 'WARNING: Gaunt integrals will be'//
     &              'scaled as HF exchange according to the selected'
               WRITE (LUPRI,*) ' DFT functional.'
            ELSE
               WRITE (LUPRI,*) 'WARNING: Gaunt integrals will be'//
     &              'included to 100% independently of the selected'
               WRITE (LUPRI,*) ' DFT functional.'
            ENDIF
         ENDIF

      ENDIF
C
C     Remove integral classes not needed for this Hamiltonians
C
      IF (LEVYLE.OR.
     &   ((BSS.or.x2c).AND..NOT.DO4C2C.AND..NOT.DO2C4C)
     &    .OR.TWOCOMP) THEN
        ISLINT = 0
        ISSINT = 0
        IGTINT = 0
        INTFLG_CHANGE = .FALSE.
      ENDIF
      IF (ZORA.AND..NOT.ZORA4) THEN
        ISSINT = 0
        IGTINT = 0
        INTFLG_CHANGE = .FALSE.
      ENDIF
C
C     Write info about the modelling of SS interactions
C
      IF(DOLVC) THEN
C
C        Calculate LV-correction
C        (SS integrals modelled by Coulombic correction)
C        ===============================================
C
         WRITE(LUPRI,'(A)')
     &   ' * SS integrals neglected:',
     &   '   Interatomic Coulombic SS-contributions modelled by',
     &   '   classical repulsion of small component atomic charges'
         IF (LVNEW) THEN
            WRITE(LUPRI,'(A)')
     &   '   found from Mulliken population analysis.'
         ELSE
            WRITE(LUPRI,'(A)')
     &   '   using tabulated charges.'
         END IF
         ! miro : print warning message, DOSSSS deactivates default
         IF (ISSINT.NE.0.AND.INTFLG_CHANGE) THEN
           WRITE(LUPRI,"(2X,A)")
     &   "WARNING:LVCORR default (deactivate with .DOSSSS) -"//
     &   " not calculating wanted SS-integrals !"
         ENDIF
C
         ISSINT = 0
         INTFLG_CHANGE = .FALSE.
      END IF

!miro: default integral flag passed to all modules (stored in dcbgen.h)
      INTGEN = 1*ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT

C
C     Write information about default integrals flags
C
      WRITE(LUPRI,9001) '* Default integral flags passed to all modules'
      WRITE(LUPRI,9002) '- LL-integrals: ',ILLINT
      WRITE(LUPRI,9002) '- LS-integrals: ',ISLINT
      WRITE(LUPRI,9002) '- SS-integrals: ',ISSINT
      WRITE(LUPRI,9002) '- GT-integrals: ',IGTINT
      IF (INTFLG_CHANGE) WRITE(LUPRI,"(3X,A)")
     & "--> chosen by user's setting of INTFLG"
 9001 FORMAT(1X,A)
 9002 FORMAT(3X,A,I5)

      IF(BSS.AND.CMPEIG) THEN
        WRITE(LUPRI,'(1X,A)')
     & '* Comparing eigenvalues between parent 4c and derived 2c'//
     & ' one-electron Hamiltonians.'
      ENDIF

#ifdef NOT_USED

! hjaaj 11-Jun-2020: DIRCON(:) is not used in this version of Dirac

C ============ Output for the basis set ===============
      WRITE(LUPRI,'(A)') ' * Basis set:'
C     ... first large component
      IF ( DIRCON(1) ) THEN
         WRITE(LUPRI,'(A)')'   - contracted large component basis set'
CMI ... contracted L comp. basis set in NOT for generated matrix elemensts
CMI     but for precursive four-component matrix elements !!!
         IF (BSS.or.x2c.OR.DO2C4C.OR.DO4C2C) THEN
             WRITE(LUPRI,'(2X,A)')
     & 'WARNING: for BSS transformation the contracted'//
     & ' L component basis is not implemented!'
CTROND        CALL QUIT('HAMINIP: contracted basis not implemented'//
CTROND     &   ' for BSS/DO2C4C/DO4C2C')
         ENDIF
      ELSE
         WRITE(LUPRI,'(A)')'   - uncontracted large component basis set'
      END IF

CMI   ... treat the small component basis set only if not beeing
CMI   in the pure two-component mode
      IF (.NOT.TWOCOMP) THEN

C     ... then small component only if not being in the two-component mode
      IF ( DIRCON(2) ) THEN
         WRITE(LUPRI,'(A)')'   - contracted small component basis set'
      ELSE
         WRITE(LUPRI,'(A)')'   - uncontracted small component basis set'
      END IF
      IF ( DIRCON(2) .AND. .NOT. URKBAL ) THEN
         WRITE(LUPRI,'(A)') '**** ERROR in **HAMILTONIAN ****'
         WRITE(LUPRI,'(A)') 'Contracted small component basis set '//
     &      'and restricted kinetic balance is not implemented.'
         CALL QUIT('*** ERROR in **HAMILTONIAN ****')
      END IF

      END IF
#endif

      IF (.NOT.TWOCOMP) THEN

      IF(URKBAL) THEN
        WRITE(LUPRI,1000) 'unrestricted'
      ELSE
        WRITE(LUPRI,1000) 'restricted'
      ENDIF
C
C     Projections on basis set
C
      IF(ISPHTR.EQ.3) THEN
         IF(.NOT.URKBAL) THEN
            WRITE(LUPRI,'(A)')
     &      ' * Default RKB projection:',
     &      '   1: Pre-projection in scalar basis',
     &      '   2: Removal of unphysical solutions'//
     &      ' (via diagonalization of free particle Hamiltonian)'
         ELSE
            WRITE(LUPRI,'(A)')
     &      ' * RKB pre-projection in scalar basis',
     &      '  without removing unphysical'//
     &      ' solutions for no external field',
     &      '  This gives extended kinetic balance and is'//
     &      '  only recommended for CONTRACTED LC basis sets'
         ENDIF
      ELSEIF(ISPHTR.EQ.1) THEN
         IF(.NOT.URKBAL) THEN
            WRITE(LUPRI,'(A)')
     &      ' * RKB projection:',
     &      '   1: No Pre-projection in scalar basis',
     &      '   2: Removal of unphysical solutions'//
     &      ' (via diagonalization of free particle Hamiltonian)',
     &      ' *** WARNING **** This procedure may give severe linear'
     &    //'  dependence problems for large basis sets'
         ELSE
            WRITE(LUPRI,'(A)')
     &      '  Spherical large component basis and Cartesian'//
     &           ' small components ',
     &      ' *** WARNING **** This procedure may give severe linear'
     &    //'  dependence problems for large basis sets'
         ENDIF
      ELSEIF(ISPHTR.EQ.0) THEN
         IF(.NOT.URKBAL) THEN
            WRITE(LUPRI,'(A)')
     &      ' * Non-standard RKB projection:',
     &      '   1: Cartesian large components',
     &      '   2: Removal of unphysical solutions'//
     &      ' (via diagonalization of free particle Hamiltonian)'
         ELSE
            WRITE(LUPRI,'(A)') ' * Unrestricted Cartesian basis'
         ENDIF
      ELSE
         WRITE(LUPRI,'(/A,I5)')
     &        ' ERROR: ISPHTR not legal:',ISPHTR
         CALL QUIT('**HAMILTONIAN: Invalid ISPHTR!')
      ENDIF

      END IF ! .not.TWOCOMP

C ==================================================================
C   ... Check that the selected One-center model is reasonable,
C       and write appropriate output.
C ==================================================================
      IF(.NOT.LEVYLE .AND. ONECAP) THEN
         IF (INTV1C .LE. 0 .OR. INTV1C .GT. 4) THEN
            CALL QUIT('  Invalid one-center model specified.')
         END IF
         IF(INTV1C.EQ.1) THEN
            WRITE(LUPRI,'(1X,A)')
     &      '* Only one-center SS parts of the e-e repulsion and n-e',
     &      '  attraction are calculated explicitly.'
            IF (ICTLV1C(2) .EQ. 1)
     &          WRITE(LUPRI,'(1X,A)')
     &       '  - Using a projection scheme for the multicenter parts.'
            IF (ICTLV1C(2) .EQ. 2)
     &          WRITE(LUPRI,'(1X,A)')
     &       '  - Testing the simulation of terms neglected in this',
     &       '    way (cheap alternative to Luuk`s projection scheme)'
         ELSE IF (INTV1C .EQ. 2) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '* One-center approximation (model 2):',
     &          '  - Multicenter SS and LS integrals are modelled by',
     &          '    classical interaction of atomic charges'
             IF(ICTLV1C(2) .EQ. 0) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Using Mulliken small component charges.'
             ELSE IF (ICTLV1C(2) .EQ. 1) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Using tabulated small component charges.'
             ELSE IF(ICTLV1C(2) .EQ. 2) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Using Mulliken small component charges.',
     &          '  - Estimating error associated with this model.'
             ELSE IF (ICTLV1C(2) .EQ. 3) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Testing HJJs penalty suggestion.'
             ELSE
                WRITE(LUPRI,'(1X,/A/)')
     &          '    Wrong one-center model specification !'
                CALL QUIT('    HERRDN : ERROR in readin. ')
             END IF
         ELSE IF(INTV1C.EQ.3) THEN
             WRITE(LUPRI,'(1X,A)')
     &       '* One-center approximation (model 3):',
     &    '  - Multicenter SS and LS modelled by Coulomb attraction',
     &    '    integrals using large and small component charges',
     &    '    based on the density matrix.'
             IF(ICTLV1C(2) .EQ. 0) THEN
                WRITE(LUPRI,'(1X,A)')
     &    '  - Charges are distributed on the nuclei and on points',
     &    '    between the nuclei where the contribution to the',
     &    '    interaction from the dipole moment is zero.'
             ELSE IF(ICTLV1C(2) .EQ. 1) THEN
                WRITE(LUPRI,'(1X,A)')
     &    '  - Option not available (no tabulated large component '//
     &         'charges).'
               CALL QUIT('    HERRDN : ERROR in readin. ')
             ELSE IF(ICTLV1C(2) .EQ. 2) THEN
                WRITE(LUPRI,'(1X,A)')
     &    '  - Using atomic centered Mulliken large and small'//
     &        ' component charges.'
             ELSE IF(ICTLV1C(2) .EQ. 3) THEN
                WRITE(LUPRI,'(1X,A)')
     &    '  - Charges are distributed on the nuclei and on the ',
     &    '    midpoints between the nuclei.'
             ELSE
                WRITE(LUPRI,'(1X,/A/)')
     &          '    Invalid one-center model specified !'
                CALL QUIT('    ONECAP : ERROR in readin. ')
             END IF
         ELSE IF(INTV1C.EQ.4) THEN
             WRITE(LUPRI,'(1X,A)')
     &       '* One-center approximation (model 4):',
     &    '  - All SS and LS integrals modelled by Coulomb attraction',
     &    '    integrals using small component atomic charges'
             IF(ICTLV1C(2) .EQ. 0) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Using Mulliken small component charges.'
             ELSE IF (ICTLV1C(2) .EQ. 1) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Using tabulated small component charges.'
             ELSE
                WRITE(LUPRI,'(1X,/A/)')
     &          '    Wrong ONECAP type 4 specification !'
                CALL QUIT('    HERRDN : ERROR in readin. ')
             END IF
         END IF
         IF (.NOT.DIRSET .AND. INTV1C.NE.2) THEN
             WRITE(LUPRI,'(1X,A)')
     &    '  - Default direct evaluation of all integrals changed',
     &    '    by ONECAP since nothing was specified in **GENEREL.',
     &    '    -->  LS and SS integrals are written to disk.'
         END IF
         IF(ONECNV .GT. D0) THEN
           WRITE(LUPRI,'(1X,A,/1X,A,1P,D15.6)')
     &  '  - One-center app. will be turned off when wavefunction is',
     &  '    converged beyond the specified threshold of : ONECNV = ',
     &       ONECNV
         END IF
C
      END IF

C
      IF(N1OPER.GT.0) THEN
        WRITE(LUPRI,'(3X,A)')
     &    '- The following finite-field one-electron '//
     &    'operators have been added:'
        DO I = 1,N1OPER
          INDXPR = IPR1OP(I)
          CALL PRSYMB(LUPRI,'-',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &     'Operator no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'-',75,0)
          CALL WRIXPR(I,INDXPR)

! Added operator must be totally symmetric within molecular point group
!
          ISYXPR = IPRPSYM(INDXPR)-1
          IF(ISYXPR.NE.0) THEN
             SYMSTP = .TRUE.
             WRITE(LUPRI,'(4X,A)') '***** **HAMILTONIAN ERROR *****'
             WRITE(LUPRI,'(4X,A,A)')
     &         'The operator is not totally symmetric under ',
     &         'the molecular point group'
               CALL REDSYM(ISYXPR,IGROUP,KSOP,NELM)
             WRITE(LUPRI,'(4X,A)') '* Possible actions:'
             WRITE(LUPRI,'(6X,A)')
     &       '1. Include operator as perturbation.'
             WRITE(LUPRI,'(6X,A,A3,A,A3,10(1X,A3))')
     &       '2. Reduce symmetry: ',GROUP,' --> ',
     &       GROUPS(IGROUP),' :{',
     &       (SYMOP(KSOP(J-1)),J=1,NELM),'}'
          ENDIF

! Operator must be symmetric under time reversal only for Kramers restricted
!  wave-functions !
          ITRXPR = IPRPTIM(INDXPR)
          IF(ITRXPR.NE.1) THEN
            WRITE(LUPRI,'(4X,A)') '***** **HAMILTONIAN WARNING *****'
            WRITE(LUPRI,'(4X,A,A8,A)') 'The operator ',PRPNAM(INDXPR),
     &      ' is not totally symmetric under time reversal !'

            WRITE(LUPRI,'(4X,A)') '* Possible actions:'
            WRITE(LUPRI,'(6X,A)')
     &     '1. Use it with Kramers unrestricted wave functions (RelCC).'
            WRITE(LUPRI,'(6X,A)')
     &     '2. Include operator as perturbation for response methods.'
            IF (.NOT.TRS_MIX) THEN
              WRITE(LUPRI,'(A)')
     &       'T+ and T- operators mixing not allowed !'
              TRSSTP=.TRUE.
            ELSE
              WRITE(LUPRI,'(/,2X,A,/)')
     &"*** Mixing of T+ and T- operators allowed at user's risk ! ***"
            ENDIF
          ENDIF
        ENDDO

      IF (WRITE_Jz_MATRIX) THEN
        IF (LINEAR) THEN
         WRITE(LUPRI,'(1X,A)')
     &  '* Writing out Jzz MO matrices for diagonalization'
        ELSE
         WRITE(LUPRI,'(1X,A)')
     &  '* No linear symmetry => no writing out Jzz MO matrices !'
        ENDIF
      ENDIF

CMI  ... treat the INPERR as in paminp.F - proper input is required
        IF (INPERR.GT.0) THEN
          WRITE(LUPRI,'(/,2x,a,i4)') 'INPERR=',INPERR
          CALL QUIT('Operator input error in ***HAMILTONIAN')
        ENDIF
      ENDIF
      IF(ONESYS) THEN
        WRITE(LUPRI,'(3X,A)') '- Two-electron part ignored '//
     &   '(one-electron system)'
      ENDIF
      IF (.NOT.LEVYLE) THEN
         IF(FREEPJ) THEN
           WRITE(LUPRI,'(3X,A)')
     &        '- Free positronic solutions projected out'
         ENDIF
         IF(VEXTPJ) THEN
           WRITE(LUPRI,'(3X,A)') '- External potential positronic '//
     &          'solutions projected out'
         ENDIF
         IF(ZORA.AND..NOT.ZORA4) WRITE(LUPRI,'(3X,A,1P,D8.2)')
     &      '- ZORA approximation using approximate density'
         IF(ZORA4) WRITE(LUPRI,'(3X,A,1P,D8.2)')
     &      '- ZORA approximation using the full density'
         IF(ZORASC) WRITE(LUPRI,'(3X,A,1P,D8.2)')
     &      '- Calculating and printing scaled ZORA orbital energies'
         IF(NOSMLV) WRITE(LUPRI,'(3X,A)')
     &      '- No nuclear attraction integrals for small components !'
      ENDIF
C
      IF(PROJEC) THEN
        WRITE(LUPRI,'(A)')
     &      ' * Variational space reduced by projection'//
     &      '  using a set of fragment orbitals.'
        WRITE(LUPRI,'(3X,A,D12.4)') '- Smallest norm accepted: ',PRJTHR
        IF(PROOWN) THEN
          WRITE(LUPRI,'(3X,A)')
     &     '- Fragment orbitals are obtained from local bases !',
     &     '  The indexing of fragment spinor sets is assumed ',
     &     '  to follow list of symmetry independent nuclei.'
          DO J = 1,NFRAG
            WRITE(LUPRI,'(3X,A,I3,A,A6)') '- Fragment spinor set ',J,
     &             ' --> ',NAMDEP(J)
            WRITE(LUPRI,'(5X,A,A6)') ' - read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(7X,A,A3)')
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(7X,A,A3,A,A72)')
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO J = 1,NFRAG
            WRITE(LUPRI,'(3X,A,I3)') '- Fragment spinor set ',J
            WRITE(LUPRI,'(3X,A,A6)') '  read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(5X,A,A3)')
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(5X,A,A3,A,A72)')
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      IF(PROJEC) THEN
        WRITE(LUPRI,'(1X,A,A)')
     &      '* Variational space reduced by projection',
     &      '  using a set of fragment orbitals.'
        WRITE(LUPRI,'(3X,A,D12.4)') '- Smallest norm accepted: ',PRJTHR
        IF(PROOWN) THEN
          WRITE(LUPRI,'(3X,A)')
     &     '- Fragment orbitals are obtained from local bases !',
     &     '  The indexing of fragment spinor sets is assumed ',
     &     '  to follow list of symmetry independent nuclei.'
          DO J = 1,NFRAG
            WRITE(LUPRI,'(3X,A,I3,A,A6)') '- Fragment spinor set ',J,
     &             ' --> ',NAMDEP(J)
            WRITE(LUPRI,'(5X,A,A6)') ' - read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(7X,A,A3)')
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(7X,A,A3,A,A72)')
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO J = 1,NFRAG
            WRITE(LUPRI,'(3X,A,I3)') '- Fragment spinor set ',J
            WRITE(LUPRI,'(5X,A,A6)') ' - read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(7X,A,A3)')
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(7X,A,A3,A,A72)')
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
#ifdef HAS_PELIB
      IF (PEQM) THEN
        WRITE(LUPRI,'(1X,A,A)')
     &  '* Environment is modeled using polarizable embedding scheme'//
     &  ' (PE library)'
      END IF
#endif
      IF (NOINPUT) GOTO 500
C
C     Process input for various program sections
C     ==========================================
C
  200 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4,5,6,7,8,9,10,11,12), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in **HAMILTONIAN.')
      ELSE
         WRITE (LUPRI,'(/,4A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in **HAMILTONIAN, error in prompt.')
      END IF
    2 CONTINUE
!       work can be easily removed by changing internal memget to alloc
        CALL SOLVIN(WORD,WORK,LWORK)
        GO TO 200


!     *DFT
    3 continue
!     this is needed until the old input scheme is replaced
      call move_to_next_star(word, lucmd)
      go to 200



    4 CONTINUE
#ifdef MOD_CAP
CMI    ... reading the CAP integrals characteristic
        CALL CAPINP(WORD,.FALSE.)
#else
CMI       dirac.x went into an infinite loop here
C       if *CAP in input but MOD_CAP is not defined !!!
        CALL QUIT(
     &  '*CAP specified, but CAP is not included in this version')
#endif
        GO TO 200
!     .FDE, frozen density embedding
    5 CONTINUE
      call fde_dirac_input(word, .true.)
      GO TO 200
    6 CONTINUE
CMI  ... reading input routine for AMFI+RELSCF
        IF((IBSS.GT.999.or.x2c).AND..NOT.NOAMFI.AND..NOT.SPINFR) THEN
          CALL AMFI_RELSCF_INP(WORD,.FALSE.)
        ELSE
          GO TO 1
        ENDIF
        GO TO 200

!       *X2C
    7 CONTINUE
!       this is needed until the old input scheme is replaced
        call move_to_next_star(word, lucmd)
        GO TO 200

!       *ECP
    8 CONTINUE
        CALL RECP_LNK_RDMENU(WORD)
        GO TO 200

!       *PCM
    9 continue
      ! this is a cheat so that *PCM section is parsed
      ! by new input reader
        call move_to_next_star(word, lucmd)
        go to 200
   10 continue
      ! this is a cheat so that *aooSOC section is parsed
      ! by new input reader
        call move_to_next_star(word, lucmd)
        go to 200
!       *PEQM
   11 CONTINUE
#ifdef HAS_PELIB
        CALL PELIB_IFC_INPUT_READER(WORD)
#else
        CALL QUIT('ERROR for *PEQM: PE library not enabled.')
#endif
        GO TO 200
!       *PCMSOL
   12 continue
      ! this is a cheat so that *PCMSOL section is parsed
      ! by new input reader
        call move_to_next_star(word, lucmd)
        go to 200
C
    1 CONTINUE

!       work can be easily removed by changing internal memget to alloc
      IF(SOLVEN) CALL SOLVIN(WORD,WORK,LWORK)

#ifdef MOD_CAP
      IF (CAP) CALL CAPINP(WORD,.TRUE.)
#endif

      IF((IBSS.GT.999.or.x2c).AND..NOT.NOAMFI.and..not.DO4C2C)THEN
         if(.not.spinfr) CALL AMFI_RELSCF_INP(WORD,.TRUE.)
      ENDIF
C
 500  CONTINUE
C
      IF (SYMSTP) CALL QUIT('**HAMILTONIAN: Non-symmetric operator '//
     &                     'under  molecular point group symmetry!')
      IF (TRSSTP) CALL QUIT('**HAMILTONIAN: Non-symmetric operator '//
     &          'under  time-reversal symmetry (T+ and T- mixing) !')

C
      CALL QEXIT('HAMINP')
      RETURN
C
 1000 FORMAT (//' Information about the ',A,' kinetic balance scheme:')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck op1ini */
      SUBROUTINE OP1INI(WORK,LWORK)
      use dirac_cfg
C***********************************************************************
C
C     Initialize one-electron integrals for standard
C     Dirac-Coulomb Hamiltonian and for the non-relativistic
C     Hamiltonian as well (MI)
C
C     Called from:  PAMINP
C
C     Written by T.Saue May 21 1996
C
C     Last revision : may 21 1996
C        MI/Febr 2006 Added nonrelativistic kinetic
C                     energy integrals (TWOCOMP=.true)
C
C***********************************************************************
      use fde_mod

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0,DM2 = -2.0D0,DHALF=0.5D0)
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbham.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3),WORK(LWORK)
      logical is_2c_rel_run
      type(fde_import) :: itmp

      CALL QENTER('OP1INI')
CMI ... common for both nonrelativistic and Dirac Hamiltonians
C
C    Overlap integrals
C    =================
C
      PNAME            = 'Overlap matrix  '
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = D1
      PLABEL(1)        = 'OVERLAP '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRHAM)
      IPOVRLAP        = INDXPR

C
C     Nuclear attraction integrals
C     ============================
C
      PNAME            = 'Nuc. attraction '
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = D1
      PLABEL(1)        = 'MOLFIELD'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRHAM)
      IPMOLFLD         = INDXPR

C
C     (static) embedding pontential
C     =============================
C
C
      if (dirac_cfg_fde) then
         call fde_get_import_info(itmp)
         if (itmp%im_vemb) then
            PNAME         = 'static emb. potential '
            IPTYP         = 1
            NPCOMP        = 1
            PFAC(1)       = D1
            PLABEL(1)     = 'FDEVEMB '
            CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRHAM)
            IPVEMB0       = INDXPR
         endif
      ENDIF


      IF(NONREL.OR.ECPCALC)THEN
CMI   nonrelativistic  kinetic energy integrals ...
      PNAME            = 'Nonrel. kinetic energy '
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = D1
      PLABEL(1)        = 'KINENERG'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRHAM)
CMI    ... utilize existing ipkinerg  in dcbham.h
      IPKINERG         = INDXPR

c
        IF (ECPCALC) THEN
c
c     Spin-orbit integrals for ECP centres ...
c     ========================================
c
        PNAME            = 'Spin-orbit Pseudopotentials'
        IPTYP            = 20 ! for the new definition
        NPCOMP           = 3
        PFAC(1)          = +1.0
        PFAC(2)          = +1.0
        PFAC(3)          = +1.0
        PLABEL(1)        = 'X1SOPP  '
        PLABEL(2)        = 'Y1SOPP  '
        PLABEL(3)        = 'Z1SOPP  '
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRHAM)
c      ... use IPSOPP defined in dcbham.h
        IPSOPP           = INDXPR
        END IF
      ELSE
CMI  ... pure Dirac Hamiltonian stuff ...
C
C     Beta matrix integrals
C     =====================
C
      PNAME            = 'Beta matrix    '
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = DM2*CVAL*CVAL
      PLABEL(1)        = 'BETAMAT '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRHAM)
      IPBETAMT         = INDXPR
C
C     Kinetic energy integrals
C     ========================
C
      PNAME            = 'Kinetic energy  '
      IPTYP            = 8
      NPCOMP           = 3
      PFAC(1)          = -CVAL
      PFAC(2)          = -CVAL
      PFAC(3)          = -CVAL
      PLABEL(1)        = 'XDIPVEL '
      PLABEL(2)        = 'YDIPVEL '
      PLABEL(3)        = 'ZDIPVEL '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRHAM)
      IPKINERG         = INDXPR

      ENDIF
CMI   ... this is common
      IF (LINEAR) THEN
C
C        Angular momentum integrals (used in case of atomic or linear
C        symmetry !)
C        ============================================================
C
         PNAME            = 'Orbital z-momentum'
         IPTYP            = 1
         NPCOMP           = 1
         PFAC(1)          = -D1
         PLABEL(1)        = 'ZANGMOM '
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                  INDXPR,ISYXPR,ITRXPR,IPRHAM)
         IPANGMOM         = INDXPR
C
         PNAME            = 'Spin z-momentum'
         IPTYP            = 12
         NPCOMP           = 1
         PFAC(1)          = DHALF
         PLABEL(1)        = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                  INDXPR,ISYXPR,ITRXPR,IPRHAM)
         IPSPNMOM         = INDXPR
      ENDIF
C
      IF (ATOMIC) THEN
C
C        Spin-orbit integrals (used in case of atomic symmetry !)
C        ============================================================
C
         PNAME            = 'Spin-orbit matrix'
         IPTYP            = 20
         NPCOMP           = 3
         PFAC(1)          = -D1
         PFAC(2)          = -D1
         PFAC(3)          = -D1
         PLABEL(1)        = 'XANGMOM '
         PLABEL(2)        = 'YANGMOM '
         PLABEL(3)        = 'ZANGMOM '
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                  INDXPR,ISYXPR,ITRXPR,IPRHAM)
         IPSPNORB         = INDXPR
      ENDIF
C
      CALL QEXIT('OP1INI')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck CIMINP */
      SUBROUTINE CIMINP(WORK,LWORK)
C***********************************************************************
      implicit none
      integer :: lwork
      real(8) :: work(lwork)
C
C     Dummy routine while the interface is being implemented
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TWOMEM */
      SUBROUTINE TWOMEM(WORK,LWORK)
C***********************************************************************
C
C     Check for memory for processing of two-electron integrals
C
C     Written by T.Saue June 12 1997
C                 MI/Febr.2006 - if TWOCOMP=.true., do not evaluate
C                 integrals over small comp. basis functions
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
C
#include "dcbgen.h"
#include "cbihr2.h"
C
      DIMENSION WORK(LWORK)
      LOGICAL LBIT,DIRFCK
      CALL QENTER('TWOMEM')
C
      CALL TITLER('Memory required by HERMIT','*',127)
      CALL PAOVEC(WORK,LWORK,0,IPRTWO)
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &   ' Maximum memory load for two-electron integral processing:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
C     always DIRFCK false to calculate memory for SO-integrals
C     which is what we normally use in Dirac (both SOFOCK for
C     Fock matrices and for MOLTRA). DIRFCK true will calculate
C     the correct for AOFOCK only. /hjaaj May 2001
C

CMI ... memmory evaluation is based on INTGEN in dcbgen.h

      DIRFCK = .FALSE.
      IF(LBIT(INTGEN,1)) THEN
        WRITE(LUPRI,'(A)') ' * LL-integrals:'
        CALL MEMLOP(1,DIRFCK)
      ENDIF

CMI ... evaluate memmory for SL,SS integrals only if not beeing in two-component mode
      IF (.NOT.TWOCOMP) THEN

        IF(LBIT(INTGEN,2)) THEN
          WRITE(LUPRI,'(A)') ' * SL-integrals:'
          CALL MEMLOP(2,DIRFCK)
        ENDIF
        IF(LBIT(INTGEN,3)) THEN
          WRITE(LUPRI,'(A)') ' * SS-integrals:'
          CALL MEMLOP(3,DIRFCK)
        ENDIF

      ENDIF

      CALL QEXIT('TWOMEM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MINANA */
      SUBROUTINE MINANA()
C***********************************************************************
C
C     Readin analysis input for the final geometry
C     Call analysis module
C
C     Written by J. Thyssen, Nov 13 1997
C
C***********************************************************************
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
C
      real(8), allocatable :: WORK(:)

C
      CALL QENTER('MINANA')
C
      IF(DOANA) THEN
         call legacy_lwork_get(LWORK)
         call alloc(WORK,LWORK,id='WORK in MINPRP')
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL ANWINP('**ANA F',WORK,LWORK)
         CLOSE(LUCMD)
         call dealloc(WORK)
C
C***********************************************************************
C*****  A N A L Y S I S    M O D U L E  ********************************
C***********************************************************************
C
         IF(DOANA) CALL PAMANA()
C
      ENDIF
C
      CALL QEXIT('MINANA')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prolab */
      SUBROUTINE PROLAB(ILAB,LAB,NDIM,NLAB,LINE,IUNIT,IBUF,WORK,LWORK)
C***********************************************************************
C
C     Process labels:
C       Read a string of new label name and list of primitive labels.
C       Form pointer ILAB relating primitive labels with new labels.
C
C     Written by T.Saue - Feb 3 1998
C
C***********************************************************************
#include"implicit.h"
#include"priunit.h"
C
      CHARACTER LINE*72,LAB(*)*12
      DIMENSION ILAB(*),IBUF(*),WORK(*)
C
      CALL IZERO(ILAB,NDIM)
      DO I = 1,NLAB
        NGRP = 1
        LINE = ' '
        READ(IUNIT,'(A72)') LINE
        LAB(I) = LINE(1:12)
        CALL NUMLST(LINE(12:72),IBUF,NDIM,1,NDIM,1,NGRP)
        DO J = 1,NGRP
          ILAB(IBUF(J)) = I
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Prelab */
      SUBROUTINE PRELAB(LABOPT,ILABDF,ILAB,LAB,NLAB)
C***********************************************************************
C    
C     Define labels according to pre-defined types
C       LABOPT = 1: individual atoms
C       LABOPT = 2: individual orbital shells
C
C***********************************************************************
#include"implicit.h"
#include"priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "ccom.h"
#include "dcblab.h"
      CHARACTER LAB(*)*12,SPDCAR*1
      DIMENSION ILAB(*)
      integer, allocatable              :: ibuf(:),ibuf2(:)
C Bit operations
#include "dcbibt.h"
      allocate(IBUF(NPLAB(0)))
      IF    (LABOPT.EQ.1) THEN
C.....Individual atoms
        ILABDF = 1
        NLAB   = 0
C.......consider only large component labels
        ICLS=1 ! Miro: runtimecheck fix
        DO ILAB1 = 1,NPLAB(1)
          IC     = IGET(IATTR(ILAB1,ILABDF))          
          ICENT  = JGET(IATTR(ILAB1,ILABDF))          
          IDEG   = KGET(IATTR(ILAB1,ILABDF))          
          KATT   = IPACK(ICLS,ICENT,IDEG,0)
C.........check if label is already defined
          DO JLAB = 1,NLAB
          IF(KATT.EQ.IBUF(JLAB)) THEN
            ILAB(ILAB1) = JLAB
            GOTO 10
          ENDIF
          ENDDO
C.........define label
          NLAB = NLAB + 1
          ILAB(ILAB1) = NLAB
          IBUF(NLAB) = KATT
          LAB(NLAB) = PLABEL(ILAB1,ILABDF)(1:8)//'    '
 10       CONTINUE
        ENDDO
      ELSEIF(LABOPT.EQ.2) THEN
C.....Individual orbital shells
C     First make pointer to orbital angular momentum
        NDIM   = MXQN*(MXQN+1)*(MXQN+2)/6
        allocate(IBUF2(NDIM))
        ITYP = 0
        DO LVAL = 1,NHTYP
          DO ICOMP = 1,KHK(LVAL)
            ITYP = ITYP + 1
            IBUF2(ITYP)=LVAL
          ENDDO
        ENDDO
        ILABDF = 1
        NLAB   = 0
C.......consider only large component labels
        ICLS=1 ! Miro: runtimecheck fix
        DO ILAB1 = 1,NPLAB(1)
          IC     = IGET(IATTR(ILAB1,ILABDF))          
          ICENT  = JGET(IATTR(ILAB1,ILABDF))          
          IDEG   = KGET(IATTR(ILAB1,ILABDF))          
          ITYP   = LGET(IATTR(ILAB1,ILABDF))          
          LVAL   = IBUF2(ITYP) - 1
          KATT   = IPACK(ICLS,ICENT,IDEG,LVAL)
C.........check if label is already defined
          DO JLAB = 1,NLAB
          IF(KATT.EQ.IBUF(JLAB)) THEN
            ILAB(ILAB1) = JLAB
            GOTO 20
          ENDIF
          ENDDO
C.........define label
          NLAB = NLAB + 1
          ILAB(ILAB1) = NLAB
          IBUF(NLAB) = KATT
          LAB(NLAB)  = PLABEL(ILAB1,ILABDF)(1:8)//
     &                 SPDCAR(LVAL)//'   '
 20       CONTINUE
        ENDDO
        deallocate(ibuf2)
      ENDIF
      deallocate(ibuf)
      RETURN
      END 
#ifdef NOT_USED
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck READREORDER*/
      SUBROUTINE READREORDER(STRING,IREORD,IMOORD,MXREORD)
C***********************************************************************
C
C     Read in new order of orbitals
C
C     Written by J. Thyssen, Sep 1 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION IMOORD(MXREORD)
      CHARACTER STRING*90

      string(90:90) = char(0)
!     ... call C-routine
        write(LUPRI,*) 'got string:',STRING
      call creadreorder(string,ireord,imoord,mxreord)
      RETURN
      END
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mvoinp */
      SUBROUTINE MVOINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for generation of modified virtual orbitals
C
C     Written by T.Saue - June 1998
C     Last revision: June 16 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbmvo.h"
C
      LOGICAL SET, NEWDEF, RESET,LBIT
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.VECMVO','.PRINT ','.INTFLG','.F2CONT',
     &            '.ADDREP','.IONMVO','.POTMVO','.WEIGHT'/
      DATA SET/.FALSE./
C
#include "ibtfun.h"
C
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
C     Initialize /CBIMVO/
C     ===================
C
      IPRMVO = 0
      ILLINT = IBTAND(INTGEN,1)
      ISLINT = IBTAND(INTGEN/2,1)
      ISSINT = IBTAND(INTGEN/4,1)
      IGTINT = IBTAND(INTGEN/8,1)
      IFCMVO = 1
      IFEMVO = 1
C     If MVOWEIGHT < 1 the Fock operator will be constructed from
C     a weighted sum of the inactive density (with weight MVOWEIGHT - 1),
C     and the density of the orbitals given by VECMVO.
      MVOWEIGHT = 1.0D0
      ADDREP = .FALSE.
      IONMVO = .TRUE.
      DO I = 1,NFSYM
        VECMVO(I) = ' '
        IF(NOCC(I).EQ.0) THEN
          VECMVO(I) = '1 '
        ELSEIF(NOCC(I).GT.0) THEN
          WRITE(VECMVO(I),'(A3,I4)') '1..',NOCC(I)
        ENDIF
      ENDDO
C
C     Process input from CBIMVO
C     =========================
C
      NEWDEF = (WORD .EQ. '*MVOCAL')
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
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in MVOINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in MVOINP.')
    1          CONTINUE
C&&&& VECMVO: Spinors defining potential
                  READ(LUCMD,'(A72)') (VECMVO(I),I=1,NFSYM)
               GO TO 100
    2          CONTINUE
C&&&& PRINT: General print level
                  READ(LUCMD,*) IPRMVO
               GO TO 100
    3          CONTINUE
C&&& INTFLG: Specify what two-itegrals should be included in this run
                  IF (IGTINT.EQ.1) THEN
                     READ(LUCMD,*) ILLINT,ISLINT,ISSINT,IGTINT
                  ELSE
                     READ(LUCMD,*) ILLINT,ISLINT,ISSINT
                  END IF
               GO TO 100
    4          CONTINUE
C&&& F2CONT: Specify Coulomb/exchange contributions
                  READ(LUCMD,*) IFCMVO,IFEMVO
               GO TO 100
    5          CONTINUE
C&&& ADDREP: Add remaining occupied orbitals as repulsion
                  ADDREP = .TRUE.
               GO TO 100
    6          CONTINUE
C&&& IONMVO: Generate MVOS by ionic Fock operator
                  IONMVO = .TRUE.
               GO TO 100
    7          CONTINUE
C&&& POTMVO: Generate MVOS from Hartree-Fock potential
                  IONMVO = .FALSE.
               GO TO 100
    8          CONTINUE
C&&& WEIGHT: Give weight of the VECMVO density, the rest is inactive density
               READ(LUCMD,*) MVOWEIGHT
               IF (MVOWEIGHT.GT.1.OR.MVOWEIGHT.LT.0) THEN
                  WRITE(LUPRI,*) 'WARNING: MVO weight not in [0,1],',
     &                 MVOWEIGHT
               ENDIF
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in MVOINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in MVOINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C
      INTMVO = ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT
      IF2MVO = IFCMVO + 2*IFEMVO
C
C     Print section
C     =============
C
      IF(.NOT.DOMVO) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' MVOGEN: Modified virtual orbitals'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I3)') ' * General print level:',IPRMVO
      IF(IONMVO) THEN
        WRITE(LUPRI,'(A)') ' * MVOs generated by ionic Fock operator'
      ELSE
        WRITE(LUPRI,'(A)') ' * MVOs generated by from DHF potential:'
      ENDIF
      WRITE(LUPRI,'(A)') ' * Vectors defining the DHF potential:'
      NTOT = 0
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECMVO(I),IDUMMY,NFBAS(I,0),
     &            -NFBAS(I,2),NFBAS(I,1),I,NVEC)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A,A3)')
     &       '- No orbitals in fermion ircop ',FREP(I)
        ELSE
          NTOT = NTOT + NVEC
          WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &    '- Orbitals in fermion ircop ',FREP(I),' :',VECMVO(I)
        ENDIF
      ENDDO
      IF(ADDREP.AND..NOT.IONMVO) THEN
        WRITE(LUPRI,'(A)')
     &     ' * Remaining occupied orbitals added as repulsion.'
      ENDIF
      IF(ILLINT.EQ.0) WRITE(LUPRI,'(A)') ' * No LL contributions'
      IF(ISLINT.EQ.0) WRITE(LUPRI,'(A)') ' * No SL contributions'
      IF(ISSINT.EQ.0) WRITE(LUPRI,'(A)') ' * No SS contributions'
      IF(IGTINT.EQ.0) WRITE(LUPRI,'(A)') ' * No GT contributions'
      IF(IF2MVO.NE.3) THEN
        WRITE(LUPRI,'(A)')' * Two-electron contributions restricted to:'
        IF(LBIT(IF2MVO,1)) WRITE(LUPRI,'(3X,A)')
     &     '+ direct contribution'
        IF(LBIT(IF2MVO,2)) WRITE(LUPRI,'(3X,A)')
     &     '+ exchange contribution'
      ENDIF
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck coninp */
      SUBROUTINE CONINP
C***********************************************************************
C
C <<< Check concistency of the specified input >>>
C
C     Written by Luuk Visscher may 1999
C
C***********************************************************************
         use dirac_cfg
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "mxgas.h"
#include "dcbkrmc.h"
#include "dcbpsi.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbana.h"
#include "cbihr2.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dgroup.h"
      LOGICAL SS_COR, SS_COR_NEEDED
      logical is_2c_rel_run

#include "ibtfun.h"
C
      CALL TITLER('Input consistency checks','*',127)

      is_2c_rel_run = .false.
      is_2c_rel_run = bss.or.x2c

      IF (DOPRP.AND.TWOCOMP.AND..NOT.is_2c_rel_run) THEN
       WRITE(LUPRI,'(1X,A/)')
     & 'NB: Properties & TWOCOMP require non-relativistic'//
     & ' form of operators !'
      ENDIF

!radovan: i think they are now rather ok tested
!     IF (SPINFR.AND.DOPRP)WRITE (LUPRI,1000)'Spinfree properties'

#ifndef MOD_UNRELEASED
      IF (DORES .AND.dirac_cfg_dft_calculation)
     &   CALL QUIT ('.RESOLVE does not work with DFT !')
#endif


C hjaaj April 2002: why not DORES .and. DOTRA ??? TODO
C     the input routines use different common blocks,
C     so I think running both should be OK. I have modified the test:
C     if DORES or DOMP2 integral transformation, then DOTRA necessary
C     for DOCCM or DOCIM as integrals on disk now cannot be correct.
Chj   IF (DORES.AND.DOTRA)
      IF ( (DORES.OR.DOMP2) .AND. ((DOCCM.OR.DOCIM) .AND. .NOT.DOTRA) )
     &CALL QUIT('4INDEX for CI/CC must be specified when RESOLV or MP2')
C
C     Correct input from the naive user who neglects the SS without applying a correction.
C     The experienced user may subtract the LVCOR later on.
C
C     Check for Hamiltonians were SS do not appear anyway..
      SS_COR_NEEDED =
     &       .NOT.(LEVYLE.OR.is_2c_rel_run.OR.(ZORA.AND..NOT.ZORA4))
C     .. or are computed exactly
      SS_COR_NEEDED = SS_COR_NEEDED.AND.(IBTAND(INTGEN/4,1).EQ.0)
C     Check for possible ways of correcting..
      SS_COR = DOLVC.OR.DOTSC
#ifdef blubb
      IF (SS_COR_NEEDED.AND..NOT.SS_COR.AND.
     &   .NOT.TWOCOMP.AND..NOT.LEVYLE.AND..NOT.is_2c_rel_run) THEN
         DOLVC = .TRUE.
         WRITE (LUPRI,'(A/A)')
     &     " * WARNING: User asked for neglect of SS-integrals "//
     &     "without invoking a correction.",
     &     "   Activated SCC method (.LVCORR), "//
     &     "see L. Visscher TCA 98 (1997) 68."
      ENDIF
#endif
C
C     Check for C1 symmetry if localization requested
C
      IF(DOLOC .AND. (NBSYM .GT. 1)) THEN
         WRITE (LUPRI,'(/A/A)')
     &   ' * WARNING : Orbital localization not'//
     &      ' yet implemented for other groups than C1 !!!',
     &   ' * WARNING : -> .LOCALIZE option is ignored'
         DOLOC = .FALSE.
      ENDIF
C
C     Checks related to Gaunt. Allow dft+"hf gaunt" for developers.
C
      IF (GAUNT) THEN

         IF (dirac_cfg_dft_calculation) WRITE (LUPRI,'(/A)')
     &   ' * INFO: Gaunt only available for Hartree-Fock, '//
     &   ' the "HF Gaunt" term will be added in DFT.'

         IF (DOTRA) WRITE (LUPRI,'(/A)')
     &   ' * WARNING: Gaunt only available for Hartree-Fock, '//
     &   ' and will be ignored in 2-electron integral transformation.'

         IF (SOFOCK) THEN
            WRITE (LUPRI,'(/A)')
     &      ' * INFO: Gaunt only implemented for AO Fock builder'//
     &      ' ---> switched to AOFOCK'
            SOFOCK = .FALSE.
         ENDIF

         IF (LONDON)
     &   CALL QUIT ('London orbitals are not implemented for GAUNT')

      END IF ! (GAUNT)

!     check whether both .URKBAL and .SPINFREE are set and exit if true
!     =================================================================

      if (spinfr .and. urkbal) then
        write(lupri,*) '.URKBAL and .SPINFREE cannot be used together'
        call quit ('.URKBAL and .SPINFREE cannot be used together')
      end if


      RETURN
 1000 FORMAT (/' **** WARNING **** ',A,' not tested.')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck procc */
      SUBROUTINE PROCC
C***********************************************************************
C
C     Print DHF occupation information
C
C     Written by J. Thyssen Aug 12 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (D2=2.0D00,D0=0.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
C
      IF (NAELEC_DHF .EQ. 0) THEN
         IF (NFSYM .EQ. 1) THEN
            WRITE(LUPRI,'(A,I5,A,I5,A)')
     &      ' * Closed shell SCF calculation with ',NELECT_DHF,
     &      ' electrons in',NISH_DHF(1),' orbitals.'
         ELSE
            WRITE(LUPRI,'(A,I5,A)')
     &      ' * Closed shell SCF calculation with ',NELECT_DHF,
     &      ' electrons in'
            WRITE(LUPRI,'(I8,A,I5,A)')
     &          NISH_DHF(1),' orbitals in Fermion irrep 1 and',
     &          NISH_DHF(2),' orbitals in Fermion irrep 2'
         END IF
      ELSE
         DTEMP = D0
         DO IFRP = 1,NFSYM
            DO IOPEN = 1,NOPEN
               DTEMP = DTEMP + NACSH(IFRP,IOPEN)*D2*DF(IOPEN)
            END DO
         END DO
         IF(AOC) THEN
           WRITE(LUPRI,'(A/)')
     &  ' * Open shell SCF calculation using Average-of-Configuration'
         ELSE
           WRITE(LUPRI,'(A/)')
     &  ' * Open shell SCF calculation using Fractional Occupation'
         ENDIF
         WRITE(LUPRI,'(A//40X,A,2(/24X,A))')
     &      ' * Shell specifications:',
     &      'Orbitals',
     &      '#electrons  irrep 1  irrep 2     f        a      alpha',
     &      '----------  -------  -------  -------  -------  -------'
         IF (NFSYM .EQ. 1) THEN
            WRITE(LUPRI,1001)
     &         NELECT_DHF,NISH_DHF(1),DF(0),DA(0),DALPHA(0)
         ELSE
            WRITE(LUPRI,1002)
     &         NELECT_DHF,NISH_DHF(1),NISH_DHF(2),DF(0),DA(0),DALPHA(0)
         END IF
         DO IOPEN = 1,NOPEN
            IF (NFSYM .EQ. 1) THEN
               WRITE(LUPRI,1011)
     &            IOPEN,
     &            DF(IOPEN)*D2*NACSH(1,IOPEN),
     &            NACSH(1,IOPEN),DF(IOPEN),DA(IOPEN),DALPHA(IOPEN)
            ELSE
               WRITE(LUPRI,1012)
     &            IOPEN,
     &            DF(IOPEN)*D2*(NACSH(1,IOPEN)+NACSH(2,IOPEN)),
     &            NACSH(1,IOPEN),NACSH(2,IOPEN),
     &            DF(IOPEN),DA(IOPEN),DALPHA(IOPEN)
            END IF
         END DO
         WRITE(LUPRI,1020)
         DTEMP      = NELECT_DHF+DTEMP
         NELECT_TOT = NINT(DTEMP)
         IF (NFSYM .EQ. 1) THEN
            WRITE(LUPRI,1031)
     &         DTEMP,NISH_DHF(1)+NASH_DHF(1)
         ELSE
            WRITE(LUPRI,1032)
     &         DTEMP,NISH_DHF(1)+NASH_DHF(1),
     &         NISH_DHF(2)+NASH_DHF(2)
         END IF
         WRITE(LUPRI,'(/,3X,A,/)')
     &      'f is the fraction occupation; a and alpha open shell '//
     &      'coupling coefficients.'
      END IF
      RETURN
 1001 FORMAT(3X,'Closed shell',I13,I14,6X,'N/A ',3(3X,F6.4))
 1002 FORMAT(3X,'Closed shell',I13,I14,I9,1X,3(3X,F6.4))
 1011 FORMAT(3X,'Open shell no.',I2,F12.2,I11,6X,'N/A ',3(3X,F6.4))
 1012 FORMAT(3X,'Open shell no.',I2,F12.2,I11,I9,1X,3(3X,F6.4))
 1020 FORMAT(3X,76('-'))
 1031 FORMAT(3X,'Total',F23.2,I11)
 1032 FORMAT(3X,'Total',F23.2,I11,I9)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rho1in */
      SUBROUTINE RHO1IN(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for 1-dimensional density plot
C
C     Written by T.Saue - May 9 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbrho1.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.MESH  ','.PRINT ','.DIFDEN','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
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
C     Initialize /DCBRHO1/
C     ===================
C
      DSTEP  = 0.01D0
      IPRHO1 = 0
      IOPTDD = 0
C
C     Process input from DCBRHO1
C     ==========================
C
      NEWDEF = (WORD .EQ. '*RHO1  ')
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
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in RHO1IN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in RHO1IN.')
    1          CONTINUE
C&&&& MESH; Step length for grid
                  READ(LUCMD,*) DSTEP
               GO TO 100
    2          CONTINUE
C&&&& PRINT: Print level
                  READ(LUCMD,*) IPRHO1
               GO TO 100
    3          CONTINUE
C&&&& DIFDEN: Attachment/detahcment density plots
                  READ(LUCMD,*) IOPTDD
                  IF(IOPTDD.EQ.0) THEN
                    WRITE(LUPRI,'(A)')
     &               '* RHO1IN: IOPTDD.EQ.0: will ignore this.'
                  ENDIF
               GO TO 100
    4          CONTINUE
               GO TO 100
    5          CONTINUE
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in RHO1IN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in RHO1IN.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT.DO1RHO) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' RHO1IN: 1-dimensional density plot:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I4)') ' * Print level: ',IPRHO1
      IF(IOPTDD.GT.0) THEN
        WRITE(LUPRI,'(A)') 'Constructing attachment density'
      ELSEIF(IOPTDD.LT.0) THEN
        WRITE(LUPRI,'(A)') 'Constructing detachment density'
      ENDIF
      WRITE(LUPRI,'(A,F7.4)') ' * Grid mesh (in Angstroms): ',DSTEP
      WRITE(LUPRI,'(A)')
     &    ' (Remember that lines can be defined by ghost centers)'
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck densin */
      SUBROUTINE DENSIN(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for 3-dimensional density plot
C
C     Written by T.Saue - Nov 5 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbrho.h"
#include "dgroup.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,
     &     CTEMP*90
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.NCUBE ','.PRINT ','.CUBADJ','.DOCUBE',
     &            '.ORBITA','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
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
C     Initialize /DCBRHO1/
C     ===================
C
      ILSCUB = 3
      IPRRHO = 0
      NCUBE(1) = 80
      NCUBE(2) = 80
      NCUBE(3) = 80
      CUBADJ(1) = 4.0D0
      CUBADJ(2) = 8.0D0
      NCUBORB(1) = 0
      NCUBORB(2) = 0
C
C     Process input from DCBRHO
C     =========================
C
      NEWDEF = (WORD .EQ. '*DENSIT')
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
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in DENSIN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in DENSIN.')
    1          CONTINUE
C&&&& NCUBE: Number of points along the sides of the cube
               READ(LUCMD,*) NCUBE(1),NCUBE(2),NCUBE(3)
               GO TO 100
    2          CONTINUE
C&&&& PRINT: Print level
                  READ(LUCMD,*) IPRRHO
               GO TO 100
    3          CONTINUE
C&&&& CUBADJ: Adjust size of cube to hopefully include 'most' density.
C&&&&         CUBADJ(1) : adjust startpoint of cube
C&&&&         CUBADJ(2) : adjust stepvectors of cube
               READ(LUCMD,*)  CUBADJ(1), CUBADJ(2)
               GO TO 100
    4          CONTINUE
C&&&  DOCUBE: Specify generation of cube-file for L- and S-comp.
               READ(LUCMD,*) ILC,ISC
               ILSCUB = ILC + 2*ISC
               GO TO 100
    5          CONTINUE
C&&& ORBITAL: Read orbitals to make cubes from
               DO I = 1,NFSYM
                  READ(LUCMD,'(A90)') CTEMP
                  NCUBORB(I) = -1
                  CALL NUMLST(CTEMP,ICUBORB(1,I),1000,1,1000,
     &                 I,NCUBORB(I))
               END DO
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in DENSIN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in DENSIN.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT.DO3RHO) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'DENSITY: 3-dimensional density plot:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I4)') '* Print level: ',IPRRHO
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wt1in */
      SUBROUTINE WT1IN(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for 1-dimensional plot of Becke weights
C
C     Written by T.Saue - May 9 2000
C
C     Modified for WT1IN by O. Fossgaard - May 2003
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbwt1.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.MESH  ','.PRINT ','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
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
C     Initialize /DCBWT1/
C     ===================
C
      DSTEP  = 0.01D0
      IPWT1 = 0
C
C     Process input from DCBBWT1
C     ==========================
C
      NEWDEF = (WORD .EQ. '*WT1   ')
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
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in WT1IN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in WT1IN.')
    1          CONTINUE
C&&&& MESH; Step length for grid
                  READ(LUCMD,*) DSTEP
               GO TO 100
    2          CONTINUE
C&&&& PRINT: Print level
                  READ(LUCMD,*) IPWT1
               GO TO 100
    3          CONTINUE
               GO TO 100
    4          CONTINUE
               GO TO 100
    5          CONTINUE
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in WT1IN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in WT1IN.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT.DO1WT) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' WT1IN: 1-dimensional plot of Becke weights:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I4)') ' * Print level: ',IPWT1
      WRITE(LUPRI,'(A,F7.4)') ' * Grid mesh (in Angstroms): ',DSTEP
      WRITE(LUPRI,'(A)')
     &    ' (Remember that lines can be defined by ghost centers)'
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck locinp */
      SUBROUTINE LOCINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for Molecular orbital localization
C
C     Written by S. Dubillard - January 2004
C     Last revision:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (NTABLE = 12)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbprj.h"
#include "dcbloc.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "nuclei.h"
C
      LOGICAL SET, NEWDEF,ORBDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,LINE*72
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.PRINT ','.PRJLOC','.HESLOC','.MAXITR',
     &            '.VECPRJ','.VECREF','.OWNBAS','.THFULL',
     &            '.THDIAG','.THGRAD','.CHECK ','.SELECT'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
#include "memint.h"
C
C     Initialize /CBILOC/
C     ===================
C
      PRJLOC = .FALSE.
      OWNBAS = .FALSE.
      ITRLOC = 100
      IPRLOC = 1
      HESLOC = 'FULL'
      LGFULL = .FALSE.
      LGDIAG = .FALSE.
      LGGRAD = .FALSE.
      LGCHCK = .FALSE.
      THFULL = 1.0D-13
      THDIAG = 1.0D-05
      THGRAD = 1.0D-07
      SELMOS = 'allocc'
C
C     Process input from CBILOC
C     =========================
C
      NEWDEF = (WORD .EQ. '*LOCALI')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100         CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12), I
                  END IF
 200                        CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &            '" not recognized under *LOCALI.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword under *LOCALI.')
 1                       CONTINUE
C   Print level
                  READ(LUCMD,*) IPRLOC
               GO TO 100
 2                       CONTINUE
C   Do localisation from projection analysis
                  PRJLOC = .TRUE.
               GO TO 100
 3                       CONTINUE
C   Model of the Hessian
                  READ(LUCMD,'(A4)') HESLOC
               GO TO 100
 4                       CONTINUE
C   Maximum number of iterations
                  READ(LUCMD,*) ITRLOC
               GO TO 100
 5                       CONTINUE
C   VECPRJ: Spinors to project onto fragment vectors
                  READ(LUCMD,'(A72)') (VECPRJ(I),I=1,NFSYM)
                  ORBDEF = .FALSE.
               GO TO 100
 6                       CONTINUE
C   VECREF: Fragment vectors
                  READ(LUCMD,*) NREFS
                  IF(NREFS.GT.MAXREF) THEN
                    WRITE(LUPRI,'(A,I5)')
     &         '* Number of fragments specified:',NREFS,
     &         '* Current maximum              :',MAXREF
                    CALL QUIT('*LOCALI: Too many fragments !')
                  ENDIF
                  DO J = 1,NREFS
                    READ(LUCMD,'(A6)') REFFIL(J)
                    READ(LUCMD,'(A72)') (VECREF(I,J),I=1,NFSYM)
                  ENDDO
               GO TO 100
 7                       CONTINUE
C   OWNBAS: Fragments (defined from NUCIND)  are calculated in their own basis sets
               OWNBAS = .TRUE.
               GO TO 100
 8                       CONTINUE
C   Set convergence criterion in scheme using full hessian
C       - change of the functional between iterations
               LGFULL = .TRUE.
               READ(LUCMD,*) THFULL
               GO TO 100
 9                       CONTINUE
C   Set convergence criterion in scheme using diagonal aproximation of the hessian
C       - change of the functional between iterations
               LGDIAG = .TRUE.
               READ(LUCMD,*) THDIAG
               GO TO 100
 10                      CONTINUE
C   Set convergence criterion - norm of the gradient
               LGGRAD = .TRUE.
               READ(LUCMD,*) THGRAD
               GO TO 100
 11                      CONTINUE
C   After last iteration using diagonal approximation calculate one step using full hessian
               LGCHCK = .TRUE.
               GO TO 100
 12                      CONTINUE
C   SELMOS: Select MOs for localization
               READ(LUCMD,'(A72)') SELMOS
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     &            '" not recognized under *LOCALI.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt under *LOCALI.')
            END IF
      END IF
 300   CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT. DOLOC) GOTO 999
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)')
     & '*LOCALI: Molecular orbital localization (Pipek-Mezey criterion)'
      CALL PRSYMB(LUPRI,'=',75,0)
C     -----------------------------------------------------------------
C     Localization method
C
      IF(PRJLOC)THEN
        WRITE(LUPRI,'(1X,A)')
     &    '* Molecular orbital localization from projection analysis'
      ELSE
        WRITE(LUPRI,'(1X,A)')
     &    '* Molecular orbital localization from Mulliken analysis'
      ENDIF
C     -----------------------------------------------------------------
C     Convergence criterion for different schemes
C
      IF(LGCHCK.AND.(HESLOC /= 'COMB'))THEN
        WRITE(LUPRI,'(1X,A,A)')
     &   ' WARNING: STOP --> you requested check of the diagonal',
     &   ' convergence (.CHECK)'
        WRITE(LUPRI,'(1X,A,A)') ' WARNING: Use this keyword only with'
        WRITE(LUPRI,'(1X,A,A)') ' WARNING: .HESLOC'
        WRITE(LUPRI,'(1X,A,A)') ' WARNING: COMBINATION'
        STOP
      ENDIF
C     Diagonal approximation
      IF(HESLOC == 'DIAG')THEN
       WRITE(LUPRI,'(1X,A)')
     &  '* The diagonal approximation of the hessian will be used'
       WRITE(LUPRI,'(1X,A)')
     &  '* Convergence criteria of the localization process :'
       IF( .NOT.LGDIAG .AND. .NOT.LGGRAD)THEN
        WRITE(LUPRI,'(1X,A,A)')'  - If the maximum is reached',
     &                         ' (= there are no positive eigenvalues)'
       ELSE
        IF(LGDIAG)THEN
         WRITE(LUPRI,'(1X,A,A)')'  - If the maximum is reached',
     &                          ' (= there are no positive eigenvalues)'
         WRITE(LUPRI,'(5X,A,D9.2)')
     &    'and the change of the functional is below threshold: ',THDIAG
         WRITE(LUPRI,'(5X,A)')
     &    '(functional = inverse mean delocalization)'
        ENDIF
        IF(LGGRAD)THEN
         WRITE(LUPRI,'(1X,A,A)')'  - If the maximum is reached',
     &                          ' (= there are no positive eigenvalues)'
         WRITE(LUPRI,'(5X,A,D9.2)')
     &    'and the norm of the gradient is below threshold: ',THGRAD
        ENDIF
       ENDIF
      ENDIF
C     Full hessian
      IF(HESLOC == 'FULL')THEN
       WRITE(LUPRI,'(1X,A)')
     &  '* Full hessian will be used'
       WRITE(LUPRI,'(1X,A)')
     &  '* Convergence criteria of the localization process :'
       WRITE(LUPRI,'(1X,A,A)')'  - If the maximum is reached',
     &                        ' (= there are no positive eigenvalues)'
       IF(LGFULL)THEN
        WRITE(LUPRI,'(3X,A,D9.2)')
     &   '- If the change of the functional is below threshold: ',THFULL
        WRITE(LUPRI,'(1X,A)')
     &   '    (functional = inverse mean delocalization)'
       ENDIF
       IF(LGGRAD)THEN
        WRITE(LUPRI,'(1X,A,D9.2)')
     &   '  - If the norm of the gradient is below threshold: ',THGRAD
       ENDIF
      ENDIF
C     Combination of the diagonal approximation and the full hessian
      IF(HESLOC == 'COMB')THEN
       WRITE(LUPRI,'(1X,A)')
     &  '* First the diagonal approximation of the hessian will be used'
       WRITE(LUPRI,'(1X,A)')
     &  '  and after reaching some convergence criterion we will'
       WRITE(LUPRI,'(1X,A)')
     &  '  switch on the construction of the full hessian'
       WRITE(LUPRI,'(1X,A)')
     &  '* Convergence criteria of the localization process :'
C
       WRITE(LUPRI,'(1X,A)')
     &  '  * First part - diagonal approximation of the hessian'
       IF( .NOT.LGDIAG .AND. .NOT.LGGRAD)THEN
        WRITE(LUPRI,'(3X,A,A)')'  - If the maximum is reached',
     &                         ' (= there are no positive eigenvalues)'
       ELSE
        IF(LGDIAG)THEN
         WRITE(LUPRI,'(3X,A,A)')'  - If the maximum is reached',
     &                          ' (= there are no positive eigenvalues)'
         WRITE(LUPRI,'(7X,A,D9.2)')
     &    'and the change of the functional is below threshold: ',THDIAG
         WRITE(LUPRI,'(7X,A)')
     &    '(functional = inverse mean delocalization)'
        ENDIF
        IF(LGGRAD)THEN
         WRITE(LUPRI,'(3X,A,A)')'  - If the maximum is reached',
     &                          ' (= there are no positive eigenvalues)'
         WRITE(LUPRI,'(7X,A,D9.2)')
     &    'and the norm of the gradient is below threshold: ',THGRAD
        ENDIF
       ENDIF
C
       WRITE(LUPRI,'(1X,A)')
     &  '  * Second part - full hessian'
       IF(LGCHCK)THEN
        WRITE(LUPRI,'(3X,A)')
     &   '  - Only one step will be performed'
       ELSE
        WRITE(LUPRI,'(3X,A,A)')'  - If the maximum is reached',
     &                         ' (= there are no positive eigenvalues)'
        IF(LGFULL)THEN
         WRITE(LUPRI,'(5X,A,D9.2)')
     &   '- If the change of the functional is below threshold: ',THFULL
         WRITE(LUPRI,'(3X,A)')
     &   '    (functional = inverse mean delocalization)'
        ENDIF
       ENDIF
      ENDIF
C     -----------------------------------------------------------------
C     Other keywords
C
      WRITE(LUPRI,'(1X,A,I6)') '* Maximum number of iterations :',ITRLOC
      WRITE(LUPRI,'(1X,A,I5)') '* Print level :',IPRLOC
      IF(LGCHCK)THEN
        WRITE(LUPRI,'(1X,A,A)')
     &   '* After the diagonal approximation will converge',
     &   ' calculate one step using full hessian'
      ENDIF
      IF(SELMOS == 'allocc')THEN
        WRITE(LUPRI,'(1X,A)') '* Only occupied MOs will be localized'
      ELSE
        WRITE(LUPRI,'(1X,A)') '* User defined MOs for localization:'
        WRITE(LUPRI,'(1X,A)') '    --> ',SELMOS
      ENDIF
C     -----------------------------------------------------------------
C     Fragments of spinors used in the projection analysis
C
      IF(PRJLOC)THEN
       IF(OWNBAS)THEN
        WRITE(LUPRI,'(1X,A)')
     &   '* Fragments are calculated in their individual bases !',
     &   '    The indexing of fragment spinor sets is assumed ',
     &   '    to follow list of symmetry independent nuclei.'
        DO J = 1,NREFS
         WRITE(LUPRI,'(1X,A,I3,A,A6)') '* Fragment spinor set ',J,
     &         ' --> ',NAMDEP(J)
         WRITE(LUPRI,'(1X,A,A6)') ' - read from file ',REFFIL(J)
         DO I = 1,NFSYM
          NVEC = 0
          CALL NUMLST(VECREF(I,J),IDUMMY,NFBAS(I,0),
     &         -NFBAS(I,2),NFBAS(I,1),I,NVEC)
          IF(NVEC.EQ.0) THEN
           WRITE(LUPRI,'(4X,A,A3)')
     &      '- No orbitals in fermion ircop ',FREP(I)
          ELSE
           WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &      '- Orbitals in fermion ircop ',FREP(I),' :',VECREF(I,J)
          ENDIF
         ENDDO
        ENDDO
       ELSE
        DO J = 1,NREFS
         WRITE(LUPRI,'(1X,A,I3)') '* Fragment spinor set ',J
         WRITE(LUPRI,'(1X,A,A6)') ' - read from file ',REFFIL(J)
         DO I = 1,NFSYM
          NVEC = 0
          CALL  NUMLST(VECREF(I,J),IDUMMY,NFBAS(I,0),
     &         -NFBAS(I,2),NFBAS(I,1),I,NVEC)
          IF(NVEC.EQ.0) THEN
           WRITE(LUPRI,'(4X,A,A3)')
     &      '- No orbitals in fermion ircop ',FREP(I)
          ELSE
           WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &      '- Orbitals in fermion ircop ',FREP(I),' :',VECREF(I,J)
          ENDIF
         ENDDO
        ENDDO
       ENDIF
      ENDIF
      WRITE(LUPRI,'(1X,A)')
     &     '* Localised MO are written in file DFCOEF !!!'
 999  CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HAMSCAN
C***********************************************************************
C
C     Search for keywords in the Hamiltonian input
C
C     Called from: GENINP
C
C     Written by T. Saue Aug 21 2006
C     Last modifications: M.Ilias, Jan.2008 for 2c4c jumps
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "cbirea.h"
      CHARACTER WORD*7,PROMPT*1
      LOGICAL ISDO2C4C,ISDO4C2C,TWOCOMP_SAVE,TWOCOMPBSS_SAVE,CNTMAT_SAVE
      LOGICAL find_keyword

      ECPCALC = .FALSE.

      ISDO2C4C=.FALSE.
      ISDO4C2C=.FALSE.

      TWOCOMP_SAVE    = TWOCOMP
      TWOCOMPBSS_SAVE = TWOCOMPBSS
      CNTMAT_SAVE     = CNTMAT
C
C     Read menu file
C     ==============
C     **** Find Hamiltonian input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      CALL UPCASE(WORD)
      IF (WORD .EQ. '**HAMIL') THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
      RETURN
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no HAMIL input found')
  930 CONTINUE
      READ (LUCMD,'(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 930
      ELSE IF(PROMPT.EQ.'.') THEN
C
C        Look for 2-component Hamiltonians
C
         IF (WORD.EQ.'.IOTC  '.OR.WORD.EQ.'.BSS   '.OR.
     &       WORD.EQ.'.X2C   '.OR.WORD.EQ.'.X2COLD'.OR.
     &       WORD.EQ.'.X2CMMF') THEN
C............active 2-component modus
           TWOCOMP=.TRUE.
C............variable telling that this is a 2-component relativistic calculation
           TWOCOMPBSS=.TRUE.
           CNTMAT=.TRUE.
         ENDIF
         IF (WORD.EQ.'.NONREL' .OR. WORD.EQ.'.FAKE2C') THEN
            TWOCOMP=.TRUE.
            NONREL =.TRUE.
         ENDIF
         IF (WORD.EQ.'.ECP   ') THEN
            ECPCALC=.TRUE.
            TWOCOMP=.TRUE.
         END IF
         IF (WORD.EQ.'.MDIRAC') MDIRAC =.TRUE.
         IF (WORD.EQ.'.LEVY-L') LEVYLE =.TRUE.
         IF (WORD.EQ.'.DO2C4C') THEN
           ISDO2C4C =.TRUE.
         ENDIF
         IF (WORD.EQ.'.DO4C2C') THEN
            ISDO4C2C =.TRUE.
         ENDIF
         GOTO 930
      ELSE IF (PROMPT .EQ. '*') THEN
         GOTO 940 ! instead of RETURN
      ELSE
        GOTO 930
      END IF

 940  CONTINUE

!     x2c module: check for the keyword indicating the 4c-DF operator as defining h1
      if(find_keyword('.X2Cmmf'))then
        isdo4c2c = .true.
      end if

      IF(ISDO2C4C.or.ISDO4C2C)THEN
        TWOCOMP    = TWOCOMP_SAVE
        TWOCOMPBSS = TWOCOMPBSS_SAVE
        CNTMAT     = CNTMAT_SAVE
      ENDIF
C
      END
#ifdef DO_KEYWORD_SCAN
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck keyword_scan */
      SUBROUTINE KEYWORD_SCAN (DOUBLE_STAR,SINGLE_STAR,DOT,NAME_LIST,
     &                         FILL_LIST,
     &                         NUMBER_OF_KEYWORDS)
C***********************************************************************
C
C     Search for keywords in the DIRAC input file. These are given
C     in the list KEYWORD_LIST if FILL_LIST is TRUE. If FILL_LIST
C     is FALSE only the total number NUMBER_OF_KEYWORDS is given back.
C
C     Input switches for the scan (may be toggled individually)
C           - DOUBLE_STAR : Include keys corresponding to sections
C           - SINGLE_STAR : Include keys corresponding to subsections
C           - DOT         : Include keys from within subsections
C           - NAME_LIST   : Include NAMELIST type input (RELCCSD)
C
C     Written by L. Visscher, october 2006
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL DOUBLE_STAR,SINGLE_STAR,DOT,NAME_LIST
      LOGICAL FILL_LIST
      CHARACTER*80 LINE
      CHARACTER*7 KEYWORD_LIST(100)
      INTEGER NUMBER_OF_KEYWORDS, LEN_LINE
C
C     Scan menu file
C     ==============
C
      NUMBER_OF_KEYWORDS = 0
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  100 READ (LUCMD,'(A)',END=120,ERR=110) LINE
      LEN_LINE = LNBLNK(LINE)
C
C     Search for regular DIRAC input
C
      IF (DOUBLE_STAR .AND. LINE(1:2) .EQ. '**') THEN
        NUMBER_OF_KEYWORDS = NUMBER_OF_KEYWORDS + 1
        IF (FILL_LIST) KEYWORD_LIST(NUMBER_OF_KEYWORDS)(1:7) = LINE(3:9)
      ELSE IF (SINGLE_STAR .AND. LINE(1:1) .EQ. '*') THEN
        NUMBER_OF_KEYWORDS = NUMBER_OF_KEYWORDS + 1
        IF (FILL_LIST) KEYWORD_LIST(NUMBER_OF_KEYWORDS)(1:7) = LINE(2:8)
      ELSE IF (DOT. AND. LINE(1:1) .EQ. '.') THEN
        NUMBER_OF_KEYWORDS = NUMBER_OF_KEYWORDS + 1
        IF (FILL_LIST) KEYWORD_LIST(NUMBER_OF_KEYWORDS)(1:7) = LINE(2:8)
      ELSE
      ENDIF
C
C     Search for namelist input (first remove initial blanks)
C
      CALL BLANKR (80,LINE)
      LEN_LINE = LNBLNK(LINE)

      IF (NAME_LIST .AND. (LINE(1:1).EQ.'&' .OR. LINE(1:1).EQ.'$')) THEN
       NUMBER_OF_KEYWORDS = NUMBER_OF_KEYWORDS + 1
       IF (FILL_LIST) KEYWORD_LIST(NUMBER_OF_KEYWORDS)(1:7) = LINE(2:8)
      END IF
C
      GO TO 100
  110 CONTINUE
      CALL QUIT('Error scanning for keywords')
  120 CONTINUE
C
      RETURN
      END
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck find_keyword */
      LOGICAL FUNCTION FIND_KEYWORD (KEYWORD)
C***********************************************************************
C
C     Position menu file at the specified keyword
C
C     Written by L. Visscher, october 2006
C     Revision Jan. 2016 hjaaj: Comparison is now case insensitive.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL KEY_FOUND
      CHARACTER LINE*80, KEYWORD*(*), KEYWORD_UPCASE*80
C
      KEY_FOUND = .FALSE.
      REWIND (LUCMD,IOSTAT=IOS)
C
      LEN_KEYWORD = MIN(LEN(KEYWORD),80)
      KEYWORD_UPCASE(1:LEN_KEYWORD) = KEYWORD(1:LEN_KEYWORD)
      CALL UPCASE(KEYWORD_UPCASE(1:LEN_KEYWORD))

  200 READ (LUCMD,'(A80)',END=220,ERR=210) LINE
      CALL UPCASE(LINE)     ! change any lower case characters to upper case
      CALL BLANKR(80,LINE)  ! move any text all the way to the left of LINE
C
      IF (LINE(1:LEN_KEYWORD) .EQ. KEYWORD_UPCASE(1:LEN_KEYWORD)) THEN
         KEY_FOUND = .TRUE.
         GO TO 220
      ENDIF
      GO TO 200
C
  210 CONTINUE
      CALL QUIT('Error scanning for keywords')
  220 CONTINUE
C
      FIND_KEYWORD = KEY_FOUND
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck cosci_inp */
      SUBROUTINE COSCI_INP(WORD,RESET,WORK,LWORK)
!***********************************************************************
!
!     Input section for COSCI-module
!
!     Written by: S. Yamamoto - Jan 31, 2007
!     Last-update: S. Yamamoto - 2007.06.08, chukyo-u
!     Last-update: S. Yamamoto - 2007.09.16, chukyo-u
!
!     Copied from KRMCINP and modified.
!     Common block data will be modified in this routine.
!     The data is stored in the common block in the cossya.h file.
!     Called by: DIRCTL->PAMINP->PSIINP->COSCI_INP
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
!
      PARAMETER (NTABLE = 10)
!     PARAMETER (RTOL = 1.0D-15,D1 = 1.0D0,D0=0.0D0,D2=2.0D00)
!
      LOGICAL SET, NEWDEF, LBIT, RESET
!     CHARACTER*4 REPNA(64)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
!     CHARACTER REPA(8)*4
      DIMENSION WORK(LWORK)
!     DIMENSION MULTB_TMP(64,64)
!
!     ...for NFSYM
#include "dgroup.h"
!
#include "cossya.h"
!
      INTEGER ISETKW(NTABLE)
!
!     ...The variable "SET" must have save attribute.
      SAVE SET
      DATA TABLE /'.PRINT ','.INACT ','.ACTIVE','.NREP  ','.NOPEN ',
     &            '.GASO  ','.IELC  ','.MINMAX','.TRDM  ','.XXXXXX'/
      DATA SET/.FALSE./
!
      DATA MAXACT/31/
      DATA LUMLF1/75/
!----------------------------------------------------------------------
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
!
!     ...Zero clear
!
      CALL IZERO(ISETKW,NTABLE)
!
!     Initialize /COSSYA0/
!     ====================
!
      KCOSINP = .FALSE.
      SYA_KTRDMAT = .FALSE.
!
!
!     Initialize /COSSYA2/
!     ====================
!
      SYA_INACT = 0
      SYA_NOPEN = 0
      SYA_GASO  = 0
      SYA_IELC  = 0
      SYA_MAXE  = 0
      SYA_MINE  = 0
      SYA_NTRDM = 0
!
!     Process input for COSSYA
!     ========================
!
      NEWDEF = (WORD .EQ. '*COSCI ')
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
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in COSCI_INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in COSCI_INP.')
!
    1          CONTINUE
!&&& .PRINT : print level in COSCI module
                  READ(LUCMD,*) SYA_PRINT
                  ISETKW(1) = 1
               GO TO 100
!
    2          CONTINUE
!&&& .INACT: number of inactive orbitals for each irrep
                  READ(LUCMD,*) (SYA_INACT(I),I=1,NFSYM)
                  ISETKW(2) = 1
               GO TO 100
!
    3          CONTINUE
!&&& .ACTIVE: number of active orbitals for each irrep
                  READ(LUCMD,*) (SYA_ACTIV(I),I=1,NFSYM)
                  ISETKW(3) = 1
               GO TO 100
!
    4          CONTINUE
!&&& .NREP: number of fermion+boson irreps.
!           NREP is used only within this routine.
                  READ(LUCMD,*) SYA_NREP
                  ISETKW(4) = 1
!.s/sya,2007.0920/
!#                NREP = SYA_NREP
!.q
               GO TO 100
!
    5          CONTINUE
!&&& .NOPEN: number of gas shell (counted by spinor)
                  READ(LUCMD,*) SYA_NOPEN
                  ISETKW(5) = 1
               GO TO 100
!
    6          CONTINUE
!&&& .GASO: GAS orbital set
                  IF( ISETKW(4) .EQ. 0 ) THEN
                     WRITE(LUPRI,'(A)') ' *** ERROR in COSCI_INP ***'
     &                  //' .NREP must be specified before .GASO !'
                     CALL QUIT('*** ERROR in COSCI_INP ***')
                  END IF
                  IF( ISETKW(5) .EQ. 0 ) THEN
                     WRITE(LUPRI,'(A)') ' *** ERROR in COSCI_INP ***'
     &                  //' .NOPEN must be specified before .GASO !'
                     CALL QUIT('*** ERROR in COSCI_INP ***')
                  END IF
                  DO I = 1, SYA_NOPEN
                     READ(LUCMD,*) (SYA_GASO(J,I),J=1,SYA_NREP)
                  END DO
                  ISETKW(6) = 1
               GO TO 100
!
    7          CONTINUE
!&&& .IELC: number of electros for each GAS
                  IF( ISETKW(5) .EQ.0 .OR. ISETKW(6) .EQ. 0 ) THEN
                     WRITE(LUPRI,'(A)') ' *** ERROR in COSCI_INP ***'
     &                  //' .NOPEN must be specified before .IELC !'
     &                  //' .GASO  must be specified before .IELC !'
                     CALL QUIT('*** ERROR in COSCI_INP ***')
                  END IF
                  DO I = 1, SYA_NOPEN
                     READ(LUCMD,*) SYA_IELC(I)
                  END DO
                  ISETKW(7) = 1
               GO TO 100
!
    8          CONTINUE
!&&& .MINMAX: allowed minimum and maxmum number of electrons in each GAS shell
                  IF( ISETKW(6) .EQ. 0 ) THEN
                     WRITE(LUPRI,'(A)') ' *** ERROR in COSCI_INP ***'
     &                  //' .GASO must be specified before .MAXE !'
                     CALL QUIT('*** ERROR in COSCI_INP ***')
                  END IF
                  IF( ISETKW(7) .EQ. 1 ) THEN
                     WRITE(LUPRI,'(A)') ' *** ERROR in COSCI_INP ***'
     &                  //' .IELC and .MAXE are exclusive each other !'
                     CALL QUIT('*** ERROR in COSCI_INP ***')
                  END IF
                  DO I = 1, SYA_NOPEN
                      READ(LUCMD,*) SYA_MINE(I),SYA_MAXE(I)
                  END DO
                  ISETKW(8) = 1
               GO TO 100
!
    9          CONTINUE
!&&& .TRDM: number of states for which TRDM is to be calculated.
                  IF( ISETKW(4) .EQ. 0 ) THEN
                     WRITE(LUPRI,'(A)') ' *** ERROR in COSCI_INP ***'
     &                  //' .NREP must be specified before .TRDM !'
                     CALL QUIT('*** ERROR in COSCI_INP ***')
                  END IF
                  READ(LUCMD,*) (SYA_NTRDM(I),I=1,SYA_NREP)
                  ISETKW(9) = 1
                  SYA_KTRDMAT = .TRUE.
               GO TO 100
!
   10          CONTINUE
!&&&  XXXXXX: invalid data
                  CALL QUIT('*** ERROR in COSCI_INP *** ' //
     &                      'invalid data XXXXXX')
               GO TO 100
!
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in COSCI_INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in COSCI_INP.')
            END IF
      END IF
  300 CONTINUE
!
!
!     Check that the number of active orbitals is less than N1
!     --------------------------------------------------------
!
      MACTVM = 0
      DO I = 1, SYA_NOPEN
         DO J = 1, NFSYM
            MACTVM = MACTVM + SYA_GASO(J,I)
         END DO
      END DO
      MACTVS = MACTVM * 2
      IF( MACTVM .GT. MAXACT ) THEN
        WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &    ' *** ERROR in COSCI_INP *** Too many active orbitals: ',
     &    MACTVM, ' -- max is ',MAXACT
        CALL QUIT(
     &     '*** ERROR in COSCI_INP *** Too many active orbitals!')
      END IF
!
!
!     Check that the number of active electrons is less than MACTVS
!     ---------------------------------------------------------------
!
      IF( ISETKW(7) .EQ. 1 ) THEN
        MAEL = 0
        DO I = 1, SYA_NOPEN
           MAEL = MAEL + SYA_IELC(I)
        END DO
        IF( MAEL .GT. MACTVS ) THEN
          WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &      ' *** ERROR in COSCI_INP *** Too many active electrons: ',
     &      MAEL, ' -- max is ',MACTVS
          CALL QUIT(
     &      '*** ERROR in COSCI_INP *** Too many active electrons!')
        END IF
!
!
!     Check MINE and MAXE
!     -------------------
!
        DO I = 1, SYA_NOPEN
          SYA_MINE(I) = SYA_IELC(I)
          SYA_MAXE(I) = SYA_IELC(I)
        END DO
      END IF
!
      IF( (ISETKW(7) .EQ. 1) .AND. (ISETKW(8) .EQ. 1) ) THEN
        WRITE(LUPRI,'(A)') ' *** ERROR in COSCI_INP ***'
     &     //' .IELC and .MINMAX are mutually exclusive !'
        CALL QUIT('*** ERROR in COSCI_INP *** ')
      END IF
!
!     If everything is OK, then set KCOSINP to .TRUE. .
!
      KCOSINP = .TRUE.
!
!     Print section
!     =============
!
      IF( .NOT.KCOSINP ) GOTO 9999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &   ' COSCI_INP: Set-up for stand-alone COSCI calculation:'
      CALL PRSYMB(LUPRI,'=',75,0)
!
      WRITE(LUPRI,'(A,I4)')
     &   ' * General print level  : ', SYA_PRINT
!
      WRITE(LUPRI,'(A,I4)')
     &   ' * # fermion irreps     : ', NFSYM
!
      WRITE(LUPRI,'(A,16I4)')
     &   ' * Inactive orbitals    : ',
     &   (SYA_INACT(I),I=1,NFSYM)
!
      WRITE(LUPRI,'(A,16I4)')
     &   ' * Active orbitals      : ',
     &   (SYA_ACTIV(I),I=1,NFSYM)
!
      WRITE(LUPRI,'(A,I4)')
     &   ' * # fermi+boson irreps : ', SYA_NREP
!
      WRITE(LUPRI,'(A,I4)')
     &   ' * Number of gas shells : ', SYA_NOPEN
!
      WRITE(LUPRI,'(A)')
     &   ' * GAS orbital set      : '
      DO I = 1, SYA_NOPEN
         WRITE(LUPRI,'(26X,32I4)') (SYA_GASO(J,I),J=1,SYA_NREP)
      END DO
!
      IF( ISETKW(7) .EQ. 1 ) THEN
      WRITE(LUPRI,'(A)')
     &   ' * IELC                 : '
      DO I = 1, SYA_NOPEN
         WRITE(LUPRI,'(26X,I4)') SYA_IELC(I)
      END DO
      END IF
!
      IF( ISETKW(8) .EQ. 1 ) THEN
      WRITE(LUPRI,'(A)')
     &   ' * MINE,MAXE            : '
      DO I = 1, SYA_NOPEN
         WRITE(LUPRI,'(26X,2I4)') SYA_MINE(I),SYA_MAXE(I)
      END DO
      END IF
!
      IF( ISETKW(9) .EQ. 1 ) THEN
      WRITE(LUPRI,'(A)')
     &   ' * Flag for TR-DMAT     : '
      WRITE(LUPRI,'(26X,L7)') SYA_KTRDMAT
      END IF
!
 9999 CONTINUE
!
      RETURN
      END

      subroutine set_lscale_dft_gaunt()
      implicit none
#include "dcbham.h"
      lscale_dft_gaunt = .true.
      end subroutine

      function word_count(s)

!       ----------------------------------------------------------------------------
        character(*), intent(in) :: s
        integer                  :: word_count
!       ----------------------------------------------------------------------------
        integer                  :: i
        logical                  :: is_blank
!       ----------------------------------------------------------------------------

        word_count = 0

        if (len(s) <= 0) return

        is_blank = .true.

        do i = 1, len(s)
          if (s(i:i) == ' ') then
            is_blank = .true.
          else if (is_blank) then
            word_count = word_count + 1
            is_blank = .false.
          end if
        end do

      end function
C======================================================================C
      SUBROUTINE LAPINP(WORD,RESET)
C***********************************************************************
C
C     Input section for Laplace module (read in old-fashioned Dirac way)
C
C     Written by B. Helmich-Paris (Sep 2015)
C
C***********************************************************************

      IMPLICIT NONE

C parameters
#include "priunit.h"
#include "lapdim.h"

C common blocks
#include "dcblap.h"

C constants:
      INTEGER NTABLE
      CHARACTER*7, CHRERR
      PARAMETER (NTABLE = 9)
      PARAMETER (CHRERR = 'LAPINP>')

C input: 
      CHARACTER*7 WORD
      LOGICAL   RESET
      INTENT(IN) RESET
      INTENT(INOUT) WORD

C local:
      INTEGER ILAP, MAXCOL, MINCOL, IOS, ITAB, ICHANG
      LOGICAL SET, NEWDEF, LWGHTS, LXPNTS
      CHARACTER PROMPT*1, TABLE(NTABLE)*7, WORD1*7
C
      SAVE SET
      DATA TABLE /'.PRINT ','.MXITER','.NUMPTS','.TOLRNG','.TOLPAR',
     &            '.TOLLAP','.STEPMX','.XPNTS ','.WGHTS '/
      DATA SET /.FALSE./

C
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
C     Initialize /DCBLAP/
C     ===================
C
      IPRLAP = 1
      NLAP   = 0
      MXITER = 100
      STEPMX = 0.3D0
      TOLRNG = 1.D-10
      TOLPAR = 1.D-15

      ! rather conservative accuracy threshold
      TOLLAP = 1.D-8

      LWGHTS = .FALSE.
      LXPNTS = .FALSE.

C
C     Process input for LAPLCE
C     ========================
C
      NEWDEF = (WORD .EQ. '*LAPLCE')
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
               DO 200 ITAB = 1, NTABLE
                  IF (TABLE(ITAB) .EQ. WORD) THEN
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9), ITAB
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in LAPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in LAPINP.')
    1          CONTINUE
C&&&& PRINT  - print level in Laplace module
                  READ(LUCMD,*) IPRLAP
               GO TO 100
    2          CONTINUE
C&&&& MXITER - maximum number of iterations for Newton type algorithms
                  READ(LUCMD,*) MXITER
               GO TO 100
    3          CONTINUE
C&&&& NLAP   - number of quadrature points
                  READ(LUCMD,*) NLAP
               GO TO 100
    4          CONTINUE
C&&&& TOLRNG - threshold for converging Newton type algorithms
                  READ(LUCMD,*) TOLRNG
               GO TO 100
    5          CONTINUE
C&&&& TOLPAR - threshold for converging Newton type algorithms
                  READ(LUCMD,*) TOLPAR
               GO TO 100
    6          CONTINUE
C&&&& TOLLAP - threshold for converging Newton type algorithms
                  READ(LUCMD,*) TOLLAP
               GO TO 100
    7          CONTINUE
C&&&& STEPMX - step length in line search algorithm
                  READ(LUCMD,*) STEPMX
               GO TO 100
    8          CONTINUE
C&&&& XPNTS  - Laplace exponents
               LXPNTS = (NLAP.GT.0)
               MAXCOL = 0
               DO WHILE (MAXCOL.LT.NLAP)
                MINCOL = MAXCOL+1
                MAXCOL = MIN(MAXCOL+5,NLAP)
                READ(LUCMD,*,IOSTAT=IOS) 
     &            (XPNTS(ILAP),ILAP=MINCOL,MAXCOL)
                IF (IOS.NE.0) THEN
                 CALL QUIT(CHRERR//'Error in reading Laplace exponents')
                ENDIF
               END DO 
               GO TO 100
    9          CONTINUE
C&&&& WGHTS  - Laplace weights
               LWGHTS = (NLAP.GT.0)
               MAXCOL = 0
               DO WHILE (MAXCOL.LT.NLAP)
                MINCOL = MAXCOL+1
                MAXCOL = MIN(MAXCOL+5,NLAP)
                READ(LUCMD,*,IOSTAT=IOS) 
     &            (WGHTS(ILAP),ILAP=MINCOL,MAXCOL)
                IF (IOS.NE.0) THEN
                 CALL QUIT(CHRERR//'Error in reading Laplace exponents')
                ENDIF
               END DO 
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     &            '" not recognized in LAPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in LAPINP.')
            END IF
      END IF
  300 CONTINUE

C
C     Post-processing section
C     =======================
C
      ! check if weights and exponents were provided in input file
      IF (LWGHTS.NEQV.LXPNTS) THEN
       CALL QUIT(CHRERR//'You need to specify both Laplace '//
     &                   'exponents and weights!')
      END IF

      ! if quadrature points not available in input read 
      ! pre-tabulated start values
      INILAP = (.NOT.(LWGHTS.AND.LXPNTS))

      ! if number of quadrature points not given or silly
      ! determine according to accuracy threshold
      FNDLAP = (NLAP.LT.1)

      IF (IPRLAP.GT.1) THEN
       WRITE(LUPRI,*) CHRERR,"print level:",IPRLAP
       WRITE(LUPRI,*) CHRERR,"maximum iteratations:",MXITER
       WRITE(LUPRI,*) CHRERR,"number of quadrature points:",NLAP
       WRITE(LUPRI,*) CHRERR,
     &  "threshold for converging Newton-Maehly algorithms:",TOLRNG 
       WRITE(LUPRI,*) CHRERR,
     &  "threshold for converging Newton-Raphson algorithms:",TOLPAR
       WRITE(LUPRI,*) CHRERR,"accuracy threshold:",TOLLAP 
       WRITE(LUPRI,*) CHRERR,
     &  "step length in line search algorithm:",STEPMX 
       IF (LWGHTS) THEN
        WRITE(LUPRI,*) CHRERR, "start weights were given:"
        DO ILAP = 1,NLAP
         WRITE(LUPRI,*) CHRERR,ILAP,WGHTS(ILAP)
        END DO
       END IF
       IF (LXPNTS) THEN
        WRITE(LUPRI,*) CHRERR, "start exponents were given:"
        DO ILAP = 1,NLAP
         WRITE(LUPRI,*) CHRERR,ILAP,XPNTS(ILAP)
        END DO
       END IF
      END IF

      RETURN
      END
C======================================================================C
!  -- end of dirac/dirrdn.F --
