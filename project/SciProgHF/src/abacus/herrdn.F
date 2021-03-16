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
#define HERRDN_DEBUG -1
C
      SUBROUTINE REAINP(WORD,RELCAL,TSTINP)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxsymm.h"
      PARAMETER ( NTABLE = 12, D1 = 1.0D0)
      LOGICAL NEWDEF, RELCAL, TSTINP
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
!
!     Parameter MAXPRD in cbirea.h is default for MAXPRI
#include "cbirea.h"
#include "ccom.h"
Caspg adding numder to set NOMOVE
#include "numder.h"
Caspg
C
      DATA TABLE /'.PRINT ', '.OLDNOR', '.MAXPRI', '.BIGVEC',
     *            '.MOLINP', '.SYMTHR', '.CM FUN', '.ZCMVAL',
     *            '.MOLPRI', '.UNCONT', '.NORTSD', '.CONTRA'/
C
      NEWDEF = (WORD .EQ. '*READIN' .OR. WORD .EQ. '*MOLBAS')
      IPREAD_ini = IPREAD
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)',IOSTAT=IOS) WORD
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
               WRITE (LUPRI,'(/4A/)') ' ERROR: Keyword "',WORD,
     *            '" not recognized for ',WORD1
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword for '//WORD1)
C                 '.PRINT '
    1          CONTINUE ! .PRINT
                  READ (LUCMD, *,IOSTAT=IOS) IPREAD
                  IF (IOS.NE.0) THEN
                     CALL QUIT('Error in reading .PRINT')
                  ENDIF
                  IF (IPREAD .EQ. IPREAD_ini) ICHANG = ICHANG - 1
               GO TO 100
    2          CONTINUE !
               GO TO 100
    3          CONTINUE ! .MAXPRI
                  READ (LUCMD, *,IOSTAT=IOS) MAXPRI
                  IF (IOS.NE.0) THEN
                    CALL QUIT('Error reading .MAXPRI')
                  ENDIF
               GO TO 100
    4          CONTINUE ! .BIGVEC
                  BIGVC = .TRUE.
               GO TO 100
    5          CONTINUE ! .MOLINP
                  LUMLCL = LUCMD
               GO TO 100
    6          CONTINUE ! .SYMTHR
                  READ (LUCMD, *,IOSTAT=IOS) TOL_SYMADD
                  IF (IOS.NE.0) THEN
                    CALL QUIT('Error in reading .SYMTHR')
                  ENDIF
               GO TO 100
    7          CONTINUE ! .CM FUN
                  READ (LUCMD,*,IOSTAT=IOS) LCMMAX, CMSTR, CMEND
                  IF (IOS.NE.0) THEN
                    CALL QUIT('Error reading .CM FUN')
                  ENDIF
                  NCMSTR = NINT(2*CMSTR + 0.01D0)
                  NCMEND = NINT(2*CMEND + 0.01D0)
               GO TO 100
    8          CONTINUE ! .ZCMVAL
                  READ (LUCMD,*,IOSTAT=IOS) ZCMVAL
                  IF (IOS.NE.0) THEN
                    CALL QUIT('Error reading .ZCMVAL')
                  ENDIF
               GO TO 100
    9          CONTINUE ! .MOLPRI
                  PRIBAS = .TRUE.
               GO TO 100
   10          CONTINUE ! .UNCONT
                  UNCONT = .TRUE.
               GO TO 100
   11          CONTINUE ! .NORTSD
C setting NOMOVE to true disables redefinition of coordinates for automatically
C determined symmetry
                  NOMOVE = .TRUE.
               GO TO 100
   12          CONTINUE ! .CONTRA
                  UNCONT = .FALSE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/4A/)') ' ERROR: Keyword "',WORD,
     *            '" not recognized for ',WORD1
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt for '//WORD1)
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for '//WORD1,0)
         IF (IPREAD .NE. IPREAD_ini) THEN
            WRITE (LUPRI,'(A,I5)')
     &         ' Print level in molecule setup (READIN):',IPREAD
         END IF
         IF (TOL_SYMADD .NE. 0.0D0) THEN
            WRITE (LUPRI,'(A,1P,G10.2)')
     &      ' Symmetry detection threshold:',TOL_SYMADD
         END IF
         IF (LUMLCL .EQ. LUCMD) THEN
            WRITE (LUPRI,'(/A)') ' MOLECULE input read from menu file'
         END IF
         IF (UNCONT) THEN
            WRITE (LUPRI,'(/A)')
     &         ' Uncontracted basis forced, irrespective of '//
     &         'basis input file.'
         END IF
         IF (MAXPRI .NE. MAXPRD) THEN
            WRITE (LUPRI,'(/A,I5)')
     *         ' Maximum number of primitives per integral block :',
     *         MAXPRI
         END IF
         IF (BIGVC) THEN
            WRITE (LUPRI,'(/A)')
     *         ' Primitives from different centers treated '//
     *         'simultaneously.'
            WRITE (LUPRI,'(/A)')
     &       ' Option .BIGVEC not implemented in present version.'
             CALL QUIT('Option .BIGVEC not allowed.')
         END IF
         IF (LCMMAX .GT. -1) THEN
            WRITE (LUPRI,'(/A,I5)')
     &           ' Center of mass basis functions used with l(max) =',
     &           LCMMAX
            WRITE (LUPRI,'(2(A,F4.1))')
     &           ' Sequence starts at n=',CMSTR,' and ends at n=', CMEND
            IF (ZCMVAL .NE. D1)
     &           WRITE (LUPRI,'(A,F10.3)') ' Net charge of ionic core',
     &           ZCMVAL
            WRITE(LUPRI,'(/A/A/A/A/A)')
     &         ' Proper reference for these basis functions:'
     &        ,' K. Kaufmann, W. Baumeister, and M. Jungen'
     &        ,' "Universal Gaussian basis sets for an optimum'
     &        ,' representation of Rydberg and continuum wavefunctions"'
     &        ,' J. Phys. B: At. Mol. Opt. Phys. 22 (1989) 2223-2240'

            IF (LCMMAX .GT. 4) THEN
               WRITE (LUPRI,'(/A)')
     &           ' FATAL ERROR: Center of mass basis functions is'//
     &           ' not defined for l(max) .gt. 4'
               CALL QUIT('Error in processing for '//WORD1)
            END IF
         END IF
         IF (PRIBAS) THEN
            WRITE (LUPRI,'(/A)')
     &           ' The molecular input will be dumped on DALTON.BAS'
         END IF
C        IF (GENCON) THEN
C           WRITE (LUPRI,'(/A)')
C    *         ' Routines for general contraction used.'
C        END IF
         WRITE (LUPRI,'()')
      END IF
      RETURN
      END
C  /* Deck reaini */
      SUBROUTINE REAINI(IPREAD_ini,RELCAL,TSTINP)
C
C     Define default THRS in CCOM
C
C
C     Initialize /CBIREA/
C
#include "implicit.h"
      LOGICAL RELCAL, TSTINP
#include "infpar.h"
#include "maxaqn.h"
#include "ccom.h"
#include "mxsymm.h"
#include "cbirea.h"
Caspg adding numder to set NOMOVE
#include "numder.h"
Caspg
      THRS = 1.0D-15
!     THRS   = 1.0D-12 ! default integral accuracy
C
      NOMOVE = .FALSE.
      IPREAD = MAX(HERRDN_DEBUG,IPREAD_ini)
      LUMLCL = 9
      BIGVC  = .FALSE.
      DIRAC  = RELCAL
      GENCON = .TRUE.
      INPTST = TSTINP
      IF (TSTINP) IPREAD = MAX(4,IPREAD)
      BASIS  = .FALSE.
      ATOMBA = .FALSE.
      PRIBAS = .FALSE.
      TOL_SYMADD = 0.0D0
      ZCMVAL = 1.0D0
      LCMMAX = -1
C
C     ***** Cartesian powers *****
C
      CALL CARPOW
C
C     In MPI runs this initialization needs to be done by the master
C     because this is the only provess to read the DIRAC.INP file.
C
C     TODO : Check that the other variables are indeed not modified
C     in DIRAC because the slaves will always get the default, unless
C     the variable is explictly communicated after reading !!
C
      IF (MYTID.EQ.MPARID) THEN
C        This is the master, initialize as usual
C        Slaves are already initialized before coming here.
         UNCONT = .FALSE.
         MAXPRI = MAXPRD
      ENDIF
      RETURN
      END
C  /* Deck readin */
      SUBROUTINE READIN(HERMIT)
C     Based on
C       READIN:Input processing routine for  M O L E C U L E.
C         Jan Almlof, Stockholm, Dec. 1971.
C         - adaption for ABACUS (T.Helgaker, University of Oslo)
C         - Occams razor by T.Saue, University of Oslo, March 10 - 1993
C           Major surgery includes:
C               - common block CCOM :
C                       - added KHK(MXQN),KCK(MXQN),NHKOFF(MXQN)
C                       - added GTOTYP(MXQN*(MXQN+1)*(MXQN+2)/6)
C               - common block NUCLEI:
C                       - added GAUNUC, GNUEXP(MXCENT)
C                       - added NAMN(MXCENT)
C               - NEW common block FRAME
C                       - added POTNUC, DPNUC(3)
C               - common block CBIREA:
C                       - IPRINT changed to IPREAD
C                       - added INPTST,FAMILY
C               - common block SHELLS:
C                       - added NLRGSH,NSMLSH,NLARGE,NSMALL,NORBS
C                       - added LCLASS(MXSHEL),NOTWOC(MXSHEL)
C               - common block MOLINP: deleted !!!
C               - common block MOLINC: deleted !!!
C               - common block HRUNIT: deleted !!!
C******************************************************************************
        USE READ_XYZFILE
        use memory_allocator
        use codata
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
C
#include "abainf.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbirea.h"
C
      CHARACTER(LEN=8) :: moltype
      LOGICAL          :: HERMIT,NEWFILE,SET,FEXIST
      CHARACTER(LEN=6) :: KEYWRD
      INTEGER          :: inpfile

C     JHS
      INTEGER :: i
      REAL :: time1,time2

#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "frame.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "chrsgn.h"
#include "aosotr.h"
#include "infpar.h"

      integer, parameter   :: LWORK = 15000000 ! do not use WORK_memsize because WORK may be allocated in calling routine
      real(8), allocatable :: WORK(:)          ! and we do not need that much memory in READIN
C
      call set_codata_values(CODSET)
C
C     A SET variable has been replaced by RDINPC. This allows external control
C     of READIN processing. The processing should be performed once every
C     iteration.
C
      IF (RDINPC) RETURN
      CALL QENTER('READIN')
      RDINPC = .TRUE.
C
!     call alloc(WORK,LWORK,id='WORK in READIN')
      allocate(WORK(LWORK))
C
      IF (IPREAD .GE. 4) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (IPREAD .GT. 0) THEN
         CALL TITLER('Output from MOLECULE input processing','*',118)
      END IF
      LUINFO = LUMLCL
C
      NEWFILE = LUINFO .NE. LUCMD
      IF (NEWFILE) THEN
          moltype='MOL'
          CALL GPINQ('MOLECULE.MOL','EXIST',FEXIST)
        IF (FEXIST) GO TO 100
          moltype='XYZ'
          CALL GPINQ('MOLECULE.XYZ','EXIST',FEXIST)
        IF (FEXIST) GO TO 100
          CALL QUIT('Unable to open a .MOL or .XYZ input file')
  100   CONTINUE
        LUINFO = -1
        CALL GPOPEN (LUINFO,'MOLECULE.'//moltype,'OLD',' ',
     &               'FORMATTED',IDUMMY,.FALSE.)
      ELSE
         moltype='MOL' ! .mol input is in the .inp file
      END IF

      REWIND(LUINFO)
C     JHS TMP reset the variables with the read .mol file
      NMLINE=0

      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(A6)',IOSTAT=IOS) KEYWRD
      IF (IOS.NE.0) THEN
         WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
         WRITE(LUPRI,*) MLINE(NMLINE)
         CALL QUIT('Error in reading...KEYWRD')
      ENDIF

C     In Dirac we do the readin the same way whatever KEYWRD is.
C
C     ********************************************
C     ***** Allocate memory **********************
C     ********************************************
C
      IF (DIRAC) THEN
C        Max number of basis sets : large, small, fit LL dens, fit SS dens,
C        Huckel
         KSETS = 5
         ISETHUCKEL=5
      ELSE
C        Max number of basis sets : regular, Huckel, fit
         IF (DOHUCKEL) THEN
            KSETS = 3
            ISETHUCKEL=2
         ELSE
            KSETS = 2
            ISETHUCKEL = -100
         END IF
      ENDIF
      KATOM  = MXATOM
      KANG   = 2*MXQN+1 ! max ang momentum needed is 2 * max(L) + 1
      KPRIM  = MAXPRI
      KBLOCK = MXSHEL
C
C     ********************************************
C     ***** Read input file and process data *****
C     ********************************************
C
C     JHS050707 .MOL files are treated in the old way
C               .XYZ files are treated using the Module READ_XYZFILE. In the long
C                 run this module should also be able to read .MOL files
C               These are the only two extensions at the moment.

      SELECT CASE(TRIM(moltype))
      CASE ('MOL')

        CALL READ_MOL(LUINFO,WORK,LWORK,
     &              KATOM,KANG,KSETS,KBLOCK,KPRIM, HERMIT)


      CASE ('XYZ')

        CALL READ_XYZ(LUINFO, LUCMD, LUPRI, WORK, LWORK,
     &                KANG,KSETS,KBLOCK,KPRIM,HERMIT)

      CASE DEFAULT
        WRITE(LUPRI, '(//A)') ' *** ERROR *** extension of MOLECULE'
     &     //' input file could not be detected.'
        CALL QUIT(' *** ERROR (READIN) unable to detect extension')
      END SELECT

C     Determine center of mass

      KGEOM = 1
      KMASS = KGEOM + 3*NUCDEP
      KNAT  = KMASS +   NUCDEP
      KNUMIS= KNAT  +   NUCDEP
      KLAST = KNUMIS+   NUCDEP
      IF (KLAST .GT. LWORK)
     &   CALL STOPIT('READIN','CMMASS',KLAST,LWORK)
      CALL CMMASS(WORK(KGEOM),WORK(KMASS),WORK(KNAT),WORK(KNUMIS),
     &            IPREAD)

      CALL GPCLOSE(LUINFO,'KEEP')
      IF (IPREAD .GT. 20) CALL REAPRI

!     call dealloc(WORK)
      deallocate(WORK)

      IF (IPREAD .GE. 4) CALL TIMER('READIN',TIMSTR,TIMEND)
      CALL QEXIT('READIN')
      RETURN
      END
C  /* Deck read_mol */
      SUBROUTINE READ_MOL(LUINFO,WORK,LWORK,
     &           KATOM,KANG,KSETS,KBLOCK,KPRIM,HERMIT)
C******************************************************************************
C
C     Reads and processes input file of molecular data
C
C       **** Temporary variables used in READIN ****
C
C       KATOM   - max number of atomic types
C       KANG    - max number of different angular momenta
C       KBLOCK  - max number of AO-blocks
C       KPRIM   - max number of primitives
C
C      (Large/small)
C       NONT            (KATOM)         - number of symmetry independent
C                                         centers for this atomic type
C       IQM             (KATOM)         - highest L-value
C       NBLCK           (KATOM)         - number of AO-blocks for a given
C                                         atomic type
C       JCO             (KANG,KATOM)    - number of AO-blocks for a given
C                                         atomic type and L-value
C       NUC             (KBLOCK)        - number of uncontracted functions
C                                         in a given AO-block
C       NRC             (KBLOCK)        - number of contracted functions
C                                         in a given AO-block
C       SEG             (KBLOCK)        - TRUE: segmented contraction
C       ALPHA           (KPRIM,KBLOCK)  - exponents
C       CPRIM           (KPRIM,KPRIM,KBLOCK)  - normalized contraction
C                                               coefficients
C       CPRIMU          (KPRIM,KPRIM,KBLOCK)  - contraction coefficients
C
C     LOCAL VARIABLES:
C       IPRIMA - runs over primitives for a given atomic type
C       ISHELL - runs over all shells
C       IPRIM  - runs over all primitives
C       NUCIND - runs over all symmetry independent centers
C******************************************************************************
      USE RECP_NTR
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "argoscom.h"
#include "argoscomch.h"
#include "consts.h"
      PARAMETER (THRMIN = 1.D-15)
#include "cbirea.h"
      REAL(KIND=8):: WORK(LWORK)

C     JHS The following local arrays are allocated at runtime to make memory
C         management easier.
      INTEGER,ALLOCATABLE :: JCO(:,:,:),
     &                       JCO2(:,:), 
     &                       JBLOCK(:),
     &                       NUC(:,:),
     &                       IQM(:,:),
     &                       NRC(:,:),
     &                       NBLCK(:,:),
     &                       ISGEN(:)
      LOGICAL,ALLOCATABLE :: SEG(:,:)
      CHARACTER(LEN=80),ALLOCATABLE :: BASREF(:,:,:)
      REAL(KIND=8),ALLOCATABLE :: ALPHA(:,:,:), 
     &                CPRIM(:,:,:,:),
     &                CPRIMU(:,:,:,:)

      LOGICAL :: HERMIT,DOOWN,ZMAT

C     JHS The following variables come from a common block and are thus kept in F77
      CHARACTER*1  KASYM(3,3),ID3,CRT
      CHARACTER*72 TTITLE(2)
      CHARACTER*11 SYMTXT
      CHARACTER*15 CLASS
      INTEGER   IFXYZ(3), blk
      LOGICAL   ANG, ADDSYM

#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#include "orgcom.h"
#include "infpar.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
      LOGICAL LINEAR, ATOMIC
#endif
#include "dgroup.h"
#include "memint.h"
C
      CALL QENTER('READ_MOL')
      ZMAT = .FALSE.
C
C***********************************************************************
C     Card 1: Keyword (INTGRL/BASIS) - read in READIN
C***********************************************************************
C***********************************************************************
C     Card 2-3: Title cards
C***********************************************************************
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(A72)',IOSTAT=IOS) TTITLE(1)
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(A72)',IOSTAT=IOS) TTITLE(2)
      IF (IPREAD .GT. 0) THEN
         CALL HEADER('Title Cards',1)
         WRITE (LUPRI,'(2X,A,/,2X,A)') TTITLE(1),TTITLE(2)
      END IF
C***********************************************************************
C     Card 4:
C       CRT     - flag for spherical harmonics or "your own scheme"
C       NONTYP  - number of atomic types
C       SYMTXT  - Symmetry operations read as a text to test if
C                 symmetry is to be added automatically
C       KCHARG  - Molecular charge
C       NSYMOP  - number of independent twofold symmetry operations of
C                 point group to be used
C       KASYM   - 3x3 character array specifying the basic symmetry operations;
C                 specify axis that are reversed under the operation
C       THRS    - integral threshold
C***********************************************************************
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(BN,A1,I4,I3,A11,A1,D10.2,6I5)',IOSTAT=IOS)
     &     CRT,NONTYP,KCHARG,SYMTXT,ID3,THRS
      IF (IOS.NE.0) THEN
       WRITE(LUPRI,*)
       WRITE(LUPRI,*) 'Error reading line no ',NMLINE,':'
       WRITE(LUPRI,*) MLINE(NMLINE)
       WRITE(LUPRI,*) 'Wanted format BN,A1,I4,I3,A11,A1,D10.2,6I5'
       CALL QUIT('Error in reading .mol ..BN,A1,I4,I3,A11,A1,D10.2,6I5')
      ENDIF
      IF(NONTYP.GT.MXATOM) GOTO 5000
C***********************************************************************
C Check integral threshold: THRS must be .gt. 0
C       ( log(thrs*thrs) taken in oneint) (851005-hjaaj)
C***********************************************************************
C
      IF (MLINE(NMLINE)(21:30) .EQ. '          ') THEN
C        nothing specified, default minimum value used
         THRS = THRMIN
      ELSE IF (THRS .LT. THRMIN) THEN
         IF (.NOT.SLAVE) THEN
            WRITE (LUPRI,'(/2X,A)') '*** WARNING from READIN ***'
            WRITE (LUPRI,'(2(/2X,A,1P,D12.2))')
     &      'Threshold for discarding integrals was',THRS,
     &      'Threshold is reset to minimum value   ',THRMIN
         END IF
         THRS = THRMIN
      ENDIF
C***********************************************************************
C Check if symmetry is to be added automatically
C       ADDSYM  - flag for automatic determination of symmetry
C***********************************************************************
      ADDSYM = (SYMTXT(1:2) .EQ. '  ')
C     Also automatic symmetry detection for Dinfh or Cinfv
      ADDSYM = (SYMTXT(1:1) .EQ. 'D') .OR. ADDSYM
      ADDSYM = (SYMTXT(1:1) .EQ. 'C') .OR. ADDSYM
      IF (.NOT. ADDSYM) THEN
         READ(SYMTXT,'(BN,I2,9A1)',IOSTAT=IOS) 
     &          NSYMOP,((KASYM(I,J),I=1,3),J=1,3)
         IF (IOS.NE.0) THEN
          WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
          WRITE(LUPRI,*) MLINE(NMLINE)
          CALL QUIT('Error in decoding SYMTXT')
         ENDIF
      ELSE
         NSYMOP = 0
      END IF
      NMLAU = NMLINE
      IF (ID3 .NE. ' ') THEN
        ANG = .TRUE.
        IF (.NOT. SLAVE)
     &       WRITE (LUPRI,'(/2X,A/10X,A,F11.8,A2/)')
     &          'Coordinates are entered in Angstroms'//
     &          ' and converted to atomic units.',
     &          '- Conversion factor : 1 bohr =',XTANG,' A'
      ELSE
         ANG = .FALSE.
      ENDIF
C***********************************************************************
C Process Cartesian to spherical or your own transformation
C     note: SPHINP will read additional input lines if DOOWN true
C***********************************************************************
      DOCART = CRT .EQ. 'C' .OR. CRT .EQ. 'c'
      DOSPHE = CRT .EQ. 'S' .OR. CRT .EQ. 's'
      DOOWN  = CRT .EQ. 'X'
      IF(DIRAC) THEN
         IF(TWOCOMP.AND.DOSPHE) THEN
C        In 2-component mode spherical Gaussians can be used,
C        and we switch off the embedded cartesian to spherical
C        gaussian transformation in LOWDIN()
            ISPHTR = 0
            IF (IPREAD .GT. 0) THEN
               WRITE(LUPRI,'(/A/A)')
     &      ' * Two-component calculations will be performed'//
     &      ' in spherical gaussians.',
     &      '   Embedded C->S transformation'//
     &      ' in the Lowdin will be skipped in this step.'
            ENDIF
         ELSE IF (.NOT. (DOCART .OR. DOOWN)) THEN
C        In 4-component mode Cartesian Gaussians are always used.  
            IF (IPREAD .GT. 0) THEN
               WRITE(LUPRI,'(/A/A)')
     &      'WARNING!!! Cartesian Gaussians must be specified',
     &      ' in molecule input file.',
     &      '           DIRAC proceeds with READIN assuming '//
     &      ' Cartesian Gaussians.'
            ENDIF
            DOCART = .TRUE.
            IF (TWOCOMPBSS.AND.DOSPHE) THEN
               ISPHTR = 3
               IF (IPREAD .GT. 0) THEN
                  WRITE(LUPRI,'(/A/A/A//A/A)')
     &      ' * Two-component calculations will be performed'//
     &      ' in spherical gaussians,',
     &      '   but the generation of the hamiltonian from  '//
     &      ' the Four-component will',
     &      '   be done in cartesian gaussians.',
     &      ' * The embedded C->S transformation in the Lowdin'//
     &      ' orthogonalization procedure', 
     &      '   will be performed for large and small components'// 
     &      ' in this step.'
               ENDIF
            ENDIF
        ENDIF
      ENDIF
      CALL SPHINP(LUINFO,WORK,LWORK,DOOWN,MXSPD)

      KATOM  = NONTYP
      KANG   = 2*MXQN+1 
      KPRIM  = MAXPRI
      KCMAT  = KPRIM*KPRIM
      ALLOCATE(IQM(KATOM,KSETS))
      ALLOCATE(NBLCK(KATOM,KSETS))
      ALLOCATE(JCO(KANG,KATOM,KSETS))
      ALLOCATE(NUC(KBLOCK,KSETS))
      ALLOCATE(NRC(KBLOCK,KSETS))
      ALLOCATE(ISGEN(KBLOCK))
      ISGEN(:) = 0
      ALLOCATE(SEG(KBLOCK,KSETS))
      ALLOCATE(BASREF(10,MXCENT,MXBSETS))
      ALLOCATE(ALPHA(KPRIM,KBLOCK,KSETS))
      ALLOCATE(CPRIM(KPRIM,KPRIM,KBLOCK,KSETS))
      ALLOCATE(CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS))

CMI  ... with gfortran the allocated space is 'spoiled' with debris
CMI      and must be cleaned before (in F90 syntax)
      JCO = 0
      NUC = 0
      NRC = 0
C
C     ************************************************
C     ***** Read orbital and geometry input data *****
C     ************************************************
C
C
      CALL BASINP(LUINFO,WORK,LWORK,
     &            IQM,NBLCK,JCO,BASREF,
     &            NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,
     &            KATOM,KANG,KSETS,KBLOCK,KPRIM,ANG,ZMAT)
C
C     *************************************************
C     ***** Add symmetry to molecule if requested *****
C     *****   Returns CLASS, the full symmetry    *****
C     *************************************************
C

      LINEAR = .FALSE.
      ATOMIC = .FALSE.
      CLASS = 'N/A'
      IF (ADDSYM) THEN
         IF (IPREAD .GE. 3) THEN
            CALL HEADER('Copy of input in .mol file before ADDSYM',0)
            WRITE (LUPRI,'(A)') (MLINE(I), I = 1, NMLINE)
            WRITE (LUPRI,'(/)')
         END IF
         CALL SYMADD(WORK,LWORK,NSYMOP,KATOM,KASYM,
     &               CLASS,TOL_SYMADD,IPREAD,.true.)
         IF (CLASS(3:4).EQ.'oo') LINEAR = .TRUE.
      END IF

C     Fix a problem with ADDSYM for some D(2) cases, also corrects for wrong user inputs (e.g. specifying X Y and XY as generators)
      CALL FIX_DUPLICATE_GENERATORS (NSYMOP,KASYM)
C
C     ***************************************
C     ***** Process symmetry input data *****
C     ***************************************
C
      CALL SYMINP(NSYMOP,KASYM,IFXYZ, CLASS)
C
C     ***************************************************
C     ***** Process orbital and geometry input data *****
C     ****************************************************
C
      CALL BASPRO(WORK,LWORK,NSYMOP,IQM,NBLCK,JCO,NUC,NRC,
     &            SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,KBLOCK,KPRIM,
     &     DOOWN,blk)
#ifdef BUILD_GEN1INT
!     initialize basis sets (large and small components) used in Gen1Int interface
      call gen1int_host_init(2, NONTYP, KATOM, NONT, IQM,
     &                       NBLCK, KANG, JCO, KBLOCK,
     &                       NUC, NRC, KPRIM, ALPHA, CPRIMU)
#endif

C
C     *************************************************
C     ***** Print orbital and geometry input data *****
C     *************************************************
C
      CALL BASOUT(WORK,LWORK,IQM,NBLCK,JCO,BASREF,
     &            NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,
     &            KBLOCK,KPRIM,HERMIT,blk)
C
      IF (ALLOCATED(ALPHA)) DEALLOCATE(ALPHA)
      IF (ALLOCATED(CPRIM)) DEALLOCATE(CPRIM)
      IF (ALLOCATED(CPRIMU)) DEALLOCATE(CPRIMU)

      IF (IPREAD .GE. 2) THEN
         CALL HEADER('Copy of .mol input',0)
         WRITE (LUPRI,'(A)') (MLINE(I), I = 1, NMLINE)
         WRITE (LUPRI,'(80A1/)') ('-',I = 1, 80)
      END IF

C
C     ****************************************
C     ***** Output on LUONEL *****************
C     ****************************************
C
      IF (HERMIT) THEN
         ALLOCATE(JCO2(KANG,KATOM))
         CALL WRONEL(TTITLE,IQM,IFXYZ,KATOM,JCO2,KANG)
         IF (ALLOCATED(JCO2)) DEALLOCATE(JCO2)
      END IF

C     ********************************
C     **** This is for ECP module ****
C     ********************************
      CALL RECP_LNK_CHECKMOL(DOCART)
      CALL RECP_LNK_IRREP(MAXREP)
      CALL RECP_LNK_RDGEO(SYMTXT,NSYMOP,KASYM,blk,
     &     IQM,JCO,KATOM,KANG)
C
C     JHS The check whether the arrays are still allocated is in principle redundant.
      IF (ALLOCATED(JCO)) DEALLOCATE(JCO)
      IF (ALLOCATED(NUC)) DEALLOCATE(NUC)
      IF (ALLOCATED(IQM)) DEALLOCATE(IQM)
      IF (ALLOCATED(NRC)) DEALLOCATE(NRC)
      IF (ALLOCATED(NBLCK)) DEALLOCATE(NBLCK)
      IF (ALLOCATED(ISGEN)) DEALLOCATE(ISGEN)
      IF (ALLOCATED(BASREF)) DEALLOCATE(BASREF)
      IF (ALLOCATED(SEG)) DEALLOCATE(SEG)

      CALL QEXIT('READ_MOL')
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     &   ' >>> READ_MOL error, no. of atomic types    ',NONTYP,
     &   '                     current maximum number ',MXATOM
      CALL QUIT('Too many atomic types')
      END ! SUBROUTINE READ_MOL
C  /* Deck baspar */
      SUBROUTINE BASPAR(NSYMOP)
C***********************************************************************
C
C       Set various basis parameters:
C               KHK(J) - number of spherical (cartesian) components for given J
C               KCK(J) - number of Cartesian components for given J
C                      - tabulate incomplete Gamma function
C                      - determine Cartesian powers
C               NHKOFF(J) - offset for components in list of l-functions
C               MAXOPR - maximum number of operations to loop over
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "cbirea.h"
C
#include "ccom.h"
#include "nuclei.h"
#include "symmet.h"
#include "aosotr.h"
C
C     ***** Tabulate Incomplete Gamma Function *****
C
      JMAX = MAX(4,4*(NHTYP - 1) + 2)
      CALL GAMTAB(JMAX)
C
C     ***** Cartesian powers *****
C
CTROND      CALL CARPOW
C
C     ***** NHKOFF() *****
C     - offset for components in list of l-functions
C
      IOFF = 0
      DO 100 I = 1, NHTYP
         NHKOFF(I) = IOFF
         IOFF = IOFF + KHK(I)
  100 CONTINUE
C
C     ***** MAXOPR *****
C     - maximum number of operations necessary to loop over
C
      MAXOPR = 1
      IF (NSYMOP .GT. 0) THEN
         MAXLO = 7
         DO 200 I = 1, NUCIND
            MAXLO = IAND(ISTBNU(I),MAXLO)
  200    CONTINUE
         II = 1
         DO 300 I = 1, NSYMOP
            II = 2*II
            IF(IAND(II,MAXLO).EQ.0) MAXOPR = II
  300    CONTINUE
      END IF
      MAXOPR = MAXOPR - 1
C
      END ! SUBROUTINE BASPAR
C  /* Deck basinp */
      SUBROUTINE BASINP(LUINFO,WORK,LWORK,IQM,NBLCK,
     &                  JCO,BASREF,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,
     &                  KATOM,KANG,KSETS,KBLOCK,KPRIM,ANG,ZMTTST)
C***********************************************************************
C
C       Read orbital and geometry input data
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (DSM = 1.0D-30)
#include "cbirea.h"
      DIMENSION WORK(LWORK)
      LOGICAL   SEG,ANG,NOORBTS,NOFITFS,ZMTTST,CNTBAS
      DIMENSION IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),
     &          NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &          SEG(KBLOCK,KSETS), ALPHA(KPRIM,KBLOCK,KSETS),
     &          CPRIM( KPRIM,KPRIM,KBLOCK,KSETS),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS), ISGEN(KBLOCK)
      CHARACTER SPDCAR*1,BSET*5,ZMATL*1,BASFIL*20,AUXFIL*20,BSKEYWORD*14
      CHARACTER*80 BASREF(10,KATOM,KSETS)
      DIMENSION IBLOCK(MXBSETS)
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbisol.h"
#include "huckel.h"
#include "aosotr.h"
#include "gencon.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "argoscom.h"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Additional scratch arrays needed for MULTIBASIS option
C
      PARAMETER (MAXBSFIL=5,KNG=MXAQN,KBLK=500,KPRM=30)
C
      INTEGER IQMT
      integer :: nblck1(maxbsfil)
      integer :: start_nuc,start_nrc,end_nuc,end_nrc
      DIMENSION JCOT(KNG,MAXBSFIL)
      DIMENSION IBLKAGNMOM(KBLK)
      DIMENSION NUCT(KBLK)
      DIMENSION NRCT(KBLK)
      DIMENSION SEGT(KBLK)
      DIMENSION ALPHAT(KPRM,KBLK)
      DIMENSION CPRIMT(KPRM,KPRM,KBLK)
      DIMENSION CPRIMUT(KPRM,KPRM,KBLK)
      CHARACTER*80 BASNAM
      CHARACTER*80 BASREFT(10,KATOM,MAXBSFIL)
      CHARACTER*20 BSFIL(MAXBSFIL)
      LOGICAL SEGT
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     MBSIP = basis-set identifier (WK/UniKA/04-11-2002).
      MBSIP = 0
C
      NUCIND = 0
      NHTYP  = 0
      NHTYPF = 0
      IBLOCK(:) = 1
      IQM(:,:)  = 0
      DIRCON(:) = .FALSE.
C
C     Run over atomic centers
C     =======================
C
      DO 100 I = 1, NONTYP
         NOORBTS = .TRUE. ! NO ORBiTalS
         NOFITFS = .TRUE. ! NO FIT FunctionS
         INUC = NUCIND
C
C***********************************************************************
C        Card 6: ATOMIC DATA
C         Q     - nuclear charge
C         NONT  - number of symmetry independent atoms
C         IQM   - maximum angular quantum number (s=1,p=2,d=3 etc.)
C         NHTYP - maximum angular quantum number for ALL orbitals
C         JCO   - number of AO-blocks for each l-value
C***********************************************************************
C
         CALL RDLINE(LUINFO)
C
C        Test if Z-matrix for first line (I.eq.1):
C

         IF (I .EQ. 1 .AND. .NOT. DIRAC) THEN
           READ (MLINE(NMLINE),'(A1)',IOSTAT=IOS) ZMATL
           IF (IOS.NE.0) THEN
             WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
             WRITE(LUPRI,*) MLINE(NMLINE)
             CALL QUIT('Error read READ..akasjs')
           ENDIF
           IF (ZMATL .NE. ' ') THEN
              ZMTTST = .TRUE.
              WRITE (LUPRI,'(/A)')
     &        'Atoms and their geometry is read in Z-matrix format.'
              IF (.NOT. BASIS) THEN
                 WRITE (LUPRI,'(/A)') '*** READ_MOL ERROR:'
                 CALL QUIT('Z-matrix only implemented with BASIS')
              END IF
           END IF
         END IF
C
C       Read/process data on centers/blocks
C       ===================================
C
        QEXP = D0
        IF(.NOT.ZMTTST) THEN
          IF (DIRAC) THEN
            DOHUCKEL = .FALSE.
            READ (MLINE(NMLINE),'(BN,1X,F9.0,I5,F20.5)',IOSTAT=IOS)
     &            Q,NONT(I),QEXP
            IF (IOS.NE.0) THEN
              WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
              WRITE(LUPRI,*) MLINE(NMLINE)
              CALL QUIT('Error read READ..akasjsxxkd')
            ENDIF
            IF (NONT(I) .LE. 0) THEN
              WRITE(LUPRI,*) 'Error in line',NMLINE,':'
              WRITE(LUPRI,*) MLINE(NMLINE)
              WRITE(LUPRI,*) 'Number of atoms became ',NONT(I)
              CALL QUIT('Error read READ..zz1')
            END IF
          ELSE IF (.NOT.ZMTTST) THEN
            NSETS = 1
            IF (DOHUCKEL) NSETS = 2
            BSET  = 'Basis'
            READ (MLINE(NMLINE),'(BN,1X,F9.0,25I5)',IOSTAT=IOS)
     &                           Q,NONT(I),
     &                           IQM(I,1), (JCO(J,I,1), J=1, IQM(I,1))
            IF (IOS.NE.0) THEN
              WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
              WRITE(LUPRI,*) MLINE(NMLINE)
              CALL QUIT('Error read READ..zzz')
            ENDIF
            IF (NONT(I) .LE. 0) THEN
              WRITE(LUPRI,*) 'Error in line',NMLINE,':'
              WRITE(LUPRI,*) MLINE(NMLINE)
              WRITE(LUPRI,*) 'Number of atoms became ',NONT(I)
              CALL QUIT('Error read READ..zz1')
            END IF
            IF (NINT(Q) .GT. 18) THEN
              NSETS = 1
              DOHUCKEL = .FALSE.
culfek: Not relevant for Dirac at this point.
c              WRITE (LUPRI,'(/A)') ' No Huckel start guess due to '//
c     &             'nuclear charge larger than 18'
            END IF
          ENDIF
        ENDIF
C
C       Read/process data on symmetry independent centers
C       =================================================
C
        IF (ZMTTST) THEN
           CALL ZMAT
           Q = CHARGE(NUCIND)
           ISOTOP(NUCIND) = 1
C          ... "Isotope=#" not implemented for ZMAT yet,
C              thus always use isotope no. 1 /hjaaj feb 2003
        ELSE
           IF (INPTST) THEN
              WRITE(LUPRI,'(/A,I5)') 'Atomic type: ',I
              WRITE(LUPRI,'(3X,A,F4.1)') 'Nuclear charge:  ',Q
              WRITE(LUPRI,'(3X,A,I5)') 'Number of symmetry '//
     &             'independent centres: ',NONT(I)
           ENDIF
           IATOMTYP(NUCIND+1:NUCIND+NONT(I)) = I
           CALL CNTINP(LUINFO,NONT(I),Q,QEXP,ANG)
C
C
C          Initialize reference strings for this atom type.
C
           BASREF(1,I,1) = ' Reference not found in input'
           DO ITYP = 1, 3
               DO IREF = 2, 10
                  BASREF(IREF,I,ITYP) = ' Not initialized'
               ENDDO
           ENDDO
C
C          Read the basis set. The input differs in DALTON/DIRAC due to
C          the optional reading of small component basis sets. The code 
C          switch starts here with first the processing of the DIRAC format.
C
           IF(DIRAC) THEN
C
C            Construct basis set for dummy centers if requested
C            ==================================================
C
             NQ = NINT(Q + 0.01D0)
             QB = Q
             IF ((LCMMAX .GT. -1) .AND. (NQ .EQ. 0)) THEN
               WRITE(LUPRI,'(/A)') 'BASINP: CM basis not implemented !'
               CALL QUIT('BASINP: CM basis not implemented in Dirac !')
C
C            Get basis set from basis set library
C            ===================================
C
             ELSE IF (BASIS) THEN
               WRITE(LUPRI,'(/A/A)')
     &           'BASINP: "BASIS" in first line not allowed in Dirac.',
     &           'Use "LARGE BASIS <basis>" for each atom type instead.'
               CALL QUIT(
     &           'BASINP: "BASIS" in first line not allowed in Dirac!')
C
             ELSE
C
C            Read basis set info from file MOLECULE.INP
C            ==========================================
C
C              Start with LARGE Component
C              ==========================
C
C              The following possibilities exist:
C
C              LARGE INTGRL    2    1    1  (IBSFLAG = 1)
C              LARGE    2    1    1 (IBSFLAG = 1)
C              LARGE EXPLICIT    2    1    1 (IBSFLAG = 1)

C              LARGE MOLFBAS H.bas (IBSFLAG = 2)
C              LARGE BASIS basisset (IBSFLAG = 3)
C              LARGE GEOM alpha beta N 2 1 1 (IBSFLAG = 4)
C              LARGE EVENTEMP alpha beta N 2 1 1 (IBSFLAG = 4)
c               (alias for LARGE GEOM)
C              LARGE WELLTEMP alpha beta gamma delta N 2 1 1 (IBSFLAG = 4)
C              LARGE FAMILY    N    2    1    1 (IBSFLAG = 4)
C              LARGE DUALFAMILY    N    2    1    1 (IBSFLAG = 4)

C              LARGE POINTCHARGE (IBSFLAG = 1)
C              LARGE NOBASIS     (IBSFLAG = 1) (alias for LARGE POINTCHARGE)
C
               CALL RDLINE(LUINFO)
               READ (MLINE(NMLINE),'(A5,1X,A14)',IOSTAT=IOS) 
     &         BSET,BSKEYWORD
C
               IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*) 'Error in reading BSET,BSKEYWORD:',
     &           BSET,BSKEYWORD
                 WRITE(LUPRI,*) 'Line',NMLINE,': ',MLINE(NMLINE)
                 CALL QUIT('Error in reading BSET,BSKEYWORD')
               ENDIF
               
               IF ( BSET .NE. 'LARGE') THEN
                  WRITE(LUPRI,'(A,I5,A,/,A)')
     $                 'Expected keyword LARGE, read line no. ',
     $                 NMLINE,':',MLINE(NMLINE)
                  GOTO 5100
               END IF
C
C              Initialize basis set flag
C
               IBSFLAG = 0
C
C              Check what kind of basis set we have
C
               IF ( BSKEYWORD(1:6) .EQ. 'POINTC'  .OR.
     &              BSKEYWORD(1:7) .EQ. 'NOBASIS') THEN
                  IBSFLAG = 1
                  IQM(I,1) = 0
                  IZATOM(I) = 0
               ELSE IF ( BSKEYWORD(1:6) .EQ. 'INTGRL' ) THEN
                  IBSFLAG = 1
                  READ (MLINE(NMLINE)(13:),*,IOSTAT=IOS)
     &                 IQM(I,1),(JCO(K,I,1),K=1,IQM(I,1))
                  IF (IOS.NE.0)THEN
                    WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                    WRITE(LUPRI,*) MLINE(NMLINE)
                    CALL QUIT('Error in reading INTGRL...')
                  ENDIF
               ELSE IF ( BSKEYWORD(1:8) .EQ. 'EXPLICIT') THEN
                  IBSFLAG = 1
                  READ (MLINE(NMLINE)(15:),*,IOSTAT=IOS)
     &                 IQM(I,1),(JCO(K,I,1),K=1,IQM(I,1))
                  IF (IOS.NE.0)THEN
                    WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                    WRITE(LUPRI,*) MLINE(NMLINE)
                    CALL QUIT('Error in reading EXPLICIT...')
                  ENDIF
               ELSE IF ( BSKEYWORD(1:7) .EQ. 'MOLFBAS') THEN
                  IBSFLAG = 2
                  READ (MLINE(NMLINE),'(14X,A20)',IOSTAT=IOS) BASFIL
                  IF (IOS.NE.0)THEN
                    WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                    WRITE(LUPRI,*) MLINE(NMLINE)
                    CALL QUIT('Error in reading MOLFBAS..')
                  ENDIF
               ELSE IF ( BSKEYWORD(1:5) .EQ. 'BASIS') THEN
                  IBSFLAG = 3
                  IF(MLINE(NMLINE)(13:14).EQ.'Q=') THEN
                    READ (MLINE(NMLINE)(15:),*,IOSTAT=IOS) QB,BASFIL
                    WRITE(LUPRI,*) '* INFO: Reading basis ',BASFIL,
     &                 ' using requested charge Q=',QB,
     &                 ' although this center has a charge of ',Q
                  ELSE
                    READ (MLINE(NMLINE),'(12X,A20)',IOSTAT=IOS) BASFIL
                  ENDIF
                  IF (IOS.NE.0)THEN
                    WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                    WRITE(LUPRI,*) MLINE(NMLINE)
                    CALL QUIT('Error in reading BASIS...')
                  ENDIF
               ELSE IF ( BSKEYWORD(1:8) .EQ. 'WELLTEMP' ) THEN
                  IBSFLAG = 4
                  READ (MLINE(NMLINE)(16:),*,ERR=9000)
     $                 (FAMPAR(K),K=1,4),NFAMEXP(1),
     $                 IQM(I,1), (JCO(K,I,1),K=1,IQM(I,1))
                  CALL GENFAMEXP(LUINFO,0)
               ELSE IF ( BSKEYWORD(1:8) .EQ. 'EVENTEMP' ) THEN
                  IBSFLAG = 4
                  READ (MLINE(NMLINE)(16:),*,ERR=9000)
     $                 (FAMPAR(K),K=1,2),NFAMEXP(1),
     $                 IQM(I,1), (JCO(K,I,1),K=1,IQM(I,1))
                  CALL GENFAMEXP(LUINFO,0)
               ELSE IF ( BSKEYWORD(1:4) .EQ. 'GEOM' ) THEN
                  IBSFLAG = 4
                  READ (MLINE(NMLINE)(11:),*,ERR=9000)
     $                 (FAMPAR(K),K=1,2),NFAMEXP(1),
     $                 IQM(I,1), (JCO(K,I,1),K=1,IQM(I,1))
                  FAMPAR(3) = D0
                  FAMPAR(4) = D0
                  CALL GENFAMEXP(LUINFO,0)
               ELSE IF ( BSKEYWORD(1:6) .EQ. 'FAMILY') THEN
                  IBSFLAG = 4
                  READ (MLINE(NMLINE),'(12X,12I5)',ERR=9000)
     $                 NFAMEXP(1), IQM(I,1), (JCO(K,I,1),K=1,IQM(I,1))
                  CALL GENFAMEXP(LUINFO,1)
               ELSE IF ( BSKEYWORD(1:10) .EQ. 'DUALFAMILY') THEN
                  IBSFLAG = 4
                  READ (MLINE(NMLINE),'(16X,12I5)',ERR=9000)
     $                 NFAMEXP(1), NFAMEXP(2),
     &                 IQM(I,1), (JCO(K,I,1),K=1,IQM(I,1))
                  CALL GENFAMEXP(LUINFO,2)
C
               ELSE IF ( BSKEYWORD(1:10) .EQ. 'MULTIBASIS') THEN
C
                  IBSFLAG = 5
                  READ (MLINE(NMLINE)(18:),*,ERR=9000)
     &                   NFILES,(BSFIL(K),K=1,NFILES)
                  IF ( NFILES .EQ. 0 ) THEN
                     WRITE(*,*) 
     &                 'PLEASE SPECIFY THE NUMBER OF BASIS SET FILES'
                     GOTO 9000
                  END IF
C
               ELSE
C                 Assume it is 'LARGE   2    1....'
                  IBSFLAG = 1
                  READ (MLINE(NMLINE)(6:),*,IOSTAT=IOS)
     &                 IQM(I,1),(JCO(K,I,1),K=1,IQM(I,1))
                  IF (IOS.NE.0) THEN
                    WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                    WRITE(LUPRI,*) MLINE(NMLINE)
                    CALL QUIT('Error in basis set reading')
                  ENDIF
               END IF

C
               IF ( IBSFLAG .EQ. 1 ) THEN
C
C                 Read explicitly given basis set.
C
                  KAOVEC = KBLOCK + 1 - IBLOCK(1)
                  CALL GTOINP(LUINFO,IQM(I,1),JCO(1,I,1),
     &                 NUC(IBLOCK(1),1),NRC(IBLOCK(1),1),
     &                 SEG(IBLOCK(1),1),
     &                 ALPHA(1,IBLOCK(1),1),CPRIM(1,1,IBLOCK(1),1),
     &                 CPRIMU(1,1,IBLOCK(1),1),ISGEN(IBLOCK(1)),
     &                 NBLOCK,KAOVEC,KPRIM)
                  NBLCK(I,1) = NBLOCK
                  BASREF(1,I,1) = 
     &            "Basis set typed explicitly in input file "
C
               ELSE IF ( IBSFLAG .EQ. 2 ) THEN
C
C                 Molfdir basis set
C
                  KBFPR = 1
                  KBFCO = KBFPR + KPRIM
                  KLAST = KBFCO + KPRIM*KPRIM
                  IF (KLAST.GT.LWORK)
     &               CALL STOPIT('BASINP','BFGINP',KLAST,LWORK)
                  DIRCON(1) = .TRUE.
                  CALL BFGINP(I,IBLOCK(1),1,WORK(KBFPR),WORK(KBFCO),
     &                IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &                NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,BASFIL,
     &                BASREF(1,1,1))
               ELSE IF ( IBSFLAG .EQ. 3 ) THEN
C
                  BASNAM = BASFIL
C
                  KAOVEC = KBLOCK + 1 - IBLOCK(1)
                  
                  CALL BASLIB(IQM(I,1),JCO(1,I,1),NUC(IBLOCK(1),1),
     &                 NRC(IBLOCK(1),1),SEG(IBLOCK(1),1),
     &                 ALPHA(1,IBLOCK(1),1),
     &                 CPRIM(1,1,IBLOCK(1),1),CPRIMU(1,1,IBLOCK(1),1),
     &                 NBLOCK,KAOVEC,KPRIM,QB,QB,DSM,
     &                 UNCONT,BASNAM,BASREF(1,I,1),IPREAD)
                  NBLCK(I,1) = NBLOCK
C
               ELSE IF ( IBSFLAG .EQ. 4 ) THEN
C
C                 Generate well-temperered/family basis set
C
                  KAOVEC = KBLOCK + 1 - IBLOCK(1)
                  CALL FAMBAS(LUINFO,IQM(I,1),JCO(1,I,1),
     &                 NUC(IBLOCK(1),1),NRC(IBLOCK(1),1),
     &                 SEG(IBLOCK(1),1),
     &                 ALPHA(1,IBLOCK(1),1),CPRIM(1,1,IBLOCK(1),1),
     &                 CPRIMU(1,1,IBLOCK(1),1),ISGEN(IBLOCK(1)),
     &                 NBLOCK,KAOVEC,KPRIM)
                  NBLCK(I,1) = NBLOCK
                  BASREF(1,I,1) = 
     &   "Well-tempered basis set typed explicitly in input file "
C
               ELSE IF ( IBSFLAG .EQ. 5 ) THEN
C
C reading of basis sets from a number of files, including when the files
C do not necessarily have the same maximum and/or minimum orbital 
C angular momenta (last modif. 11/09/2003 by gomes) 
C
C we will move through the basis set files using the dummy variable iblock1
C 
                  IBLOCK1 = 1
                  MAXIQM  = 0
C
C we make sure we don't have any junk on the temporary variables 
C
                  JCOT(1:KNG,NFILES)=0
                  NUCT(1:KBLK)=0
                  NRCT(1:KBLK)=0
                  NBLCK1(:)=0
C
                  BASREF(1,I,1) =
     &            "Basis set read from multiple files"
C
                  DO M=1,NFILES
                     BASNAM = BSFIL(M)
                     KAOVEC  = KBLOCK + 1 - IBLOCK1 
C
                     BASREF(3*(M-1)+2,I,1) = " " 
                     BASREF(3*(M-1)+3,I,1) = " ref. in file "//BSFIL(M) 
                     
C Note that the Variables with T have dimension KPRM instead of KPRIM!
                     CALL BASLIB(IQMT,JCOT(1,M),
     &                   NUCT(IBLOCK1),
     &                   NRCT(IBLOCK1),SEGT(IBLOCK1),
     &                   ALPHAT(1,IBLOCK1),
     &                   CPRIMT(1,1,IBLOCK1),
     &                   CPRIMUT(1,1,IBLOCK1),
     &                   NBLOCK,KAOVEC,KPRM,Q,Q,DSM,
     &                   UNCONT,BASNAM,BASREFT(1,I,M),IPREAD)

                     BASREF(3*(M-1)+4,I,1) = BASREFT(1,I,M)
C
C to help in the reordering, we create a table mapping the orbital
C angular momentum of a block to the number of the block in the
C temp storage 
C 
                     IDBLK = IBLOCK1
                     IDANG = 1
                     DO IRR=1,IQMT
                        DO IRS=1,JCOT(IRR,M)
                           IBLKAGNMOM(IDBLK) = IDANG
                           IDBLK = IDBLK + 1
                        END DO
                        IDANG = IDANG + 1
                     END DO
C
                     IF ( IQMT .GT. MAXIQM ) THEN
                        MAXIQM = IQMT
                     ENDIF
                     NBLCK1(M)=NBLOCK
                     IBLOCK1 = IBLOCK1 + NBLOCK
                  END DO
C
C after reading all files, we reorder the data wrt orbital angular 
C momentum. in this reordering variable kt is the "final" index of a given block
C
                  IQM(I,1) = MAXIQM
                  KT = IBLOCK(1)
                  DO IQ=1,MAXIQM
                    DO IR=1,IBLOCK1 - 1
C
C in filling up the permanent storage with data, care should be taken to skip 
C zero-sized blocks that arise when we have only parts of a basis
C set in one file, as when the polarization functions are in one file and
C the scf set is in another...
C
                      IF ( IBLKAGNMOM(IR) .EQ. IQ ) THEN
                        IF ( NUCT(IR) .GT. 0 ) THEN
                 
                          start_nuc=nuc(kt,1)+1
                          end_nuc=nuc(kt,1)+nuct(ir)
                          start_nrc=nrc(kt,1)+1
                          end_nrc=nrc(kt,1)+nrct(ir)
CMI           ... filling up !!! the JCO must be INITIALIZED before !!! 
              ALPHA(start_nuc:end_nuc,KT,1)
     &                         = ALPHAT(1:NUCT(IR),IR)
              CPRIM(start_nuc:end_nuc,start_nrc:end_nrc,KT,1)
     &                         = CPRIMT(1:NUCT(IR),1:NRCT(IR),IR)    
              CPRIMU(start_nuc:end_nuc,start_nrc:end_nrc,KT,1)
     &                         = CPRIMUT(1:NUCT(IR),1:NRCT(IR),IR)    

                          JCO(IQ,I,1)= JCO(IQ,I,1)+ 1
                          NUC(KT,1)  = NUC(KT,1) + NUCT(IR)
                          NRC(KT,1)  = NRC(KT,1) + NRCT(IR)
                          SEG(KT,1)  = SEGT(IR)
                          ISGEN(KT)  = 0
                          KT = KT + 1
                        ENDIF 
                      ENDIF
                    END DO
                  END DO
C
C finally, we update nblck(i,1) with the total number of non-zero sized
C blocks read, subtracting one only because on the loops above kt was 
C incremented after the last block was put in place ...
C
                 NBLCK(I,1) = KT - 1 
C
C ... and we are done
C
               ENDIF
C
C              Print out basis set
C
               NOORBTS    = NOORBTS.AND.IQM(I,1).EQ.0
               IF (IQM(I,1).GT.0) THEN
                  IF(INPTST) THEN
                     WRITE(LUPRI,'(3X,A)') 'Large component set:'
                     WRITE(LUPRI,'(6X,A,I5,A,12I5)')
     &                    'Max.ang.quantum no.:',(IQM(I,1)-1),
     &                    '  Blocks:',(JCO(K,I,1), K=1,IQM(I,1))
                  ENDIF
                  NHTYP  =   MAX(NHTYP,IQM(I,1))
                  IF(NHTYP.GT.MXQN) GOTO 5000
               ENDIF
C
C              Check if large component basis is contracted
C              ============================================
C
               CNTBAS = .FALSE.
               DO JBLCK = IBLOCK(1), NBLCK(I,1)+IBLOCK(1)-1
                  IF ((NUC(JBLCK,1).NE.NRC(JBLCK,1)).AND.
     &                (NRC(JBLCK,1).GT.0)) CNTBAS = .TRUE.
               END DO
               IF ( CNTBAS .AND.
     $              (IBSFLAG .EQ. 3 .OR. IBSFLAG .EQ. 1 )) THEN
                  DIRCON(1) = .TRUE.
               END IF

               IF (.NOT.TWOCOMP) THEN
C
C
C
C              Next is the SMALL Component
C              ==========================
C
C
C              The following four possibilities exist:
C
C              SMALL INTGRL     2    1    1  (IBSFLAG = 1)
C              SMALL    2    1    1 (IBSFLAG = 1)
C              SMALL EXPLICIT     2    1    1  (IBSFLAG = 1)
C              SMALL MOLFBAS H.bas (IBSFLAG = 2)
C              SMALL KINBAL (IBSFLAG = 4)
C
               CALL FLSHFO(LUPRI)
               CALL RDLINE(LUINFO)
               READ (MLINE(NMLINE),'(A5,1X,A8)',IOSTAT=IOS)
     &             BSET,BSKEYWORD
               IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                 WRITE(LUPRI,*) MLINE(NMLINE)
                 CALL QUIT('Error in reading...jsjuduu')
               ENDIF
C
C              Initialize basis set flag
C
               IBSFLAG = 0
C
C              Check whether the small component is specified in input
C
               IF ( BSET .NE. 'SMALL' ) THEN
C
C                    No 'SMALL' keyword given.
C                    Assume kinetic balance
C
C                    Discard the line we just read.
C
                     NMLINE = NMLINE - 1
                     IBSFLAG = 4
C
               ELSE
C
C
C                 Check what kind of basis set we have
C
                  IF ( BSKEYWORD(1:6) .EQ. 'INTGRL' ) THEN
                     IBSFLAG = 1
                     READ (MLINE(NMLINE),'(12X,12I5)',IOSTAT=IOS)
     &                     IQM(I,2),(JCO(K,I,2),K=1,IQM(I,2))
                     IF (IOS.NE.0) THEN
                       WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                       WRITE(LUPRI,*) MLINE(NMLINE)
                       CALL QUIT('Error in reading...ssdk123uu')
                     ENDIF
                  ELSE IF ( BSKEYWORD(1:8) .EQ. 'EXPLICIT') THEN
                     IBSFLAG = 1
                     READ (MLINE(NMLINE)(15:),*,IOSTAT=IOS)
     &                     IQM(I,2),(JCO(K,I,2),K=1,IQM(I,2))
                     IF (IOS.NE.0) THEN
                       WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                       WRITE(LUPRI,*) MLINE(NMLINE)
                       CALL QUIT('Error in reading...3323dd')
                     ENDIF
                  ELSE IF ( BSKEYWORD .EQ. 'MOLFBAS') THEN
                     IBSFLAG = 2
                     READ (MLINE(NMLINE),'(14X,A20)',IOSTAT=IOS) BASFIL
                     IF (IOS.NE.0) THEN
                       WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                       WRITE(LUPRI,*) MLINE(NMLINE)
                       CALL QUIT('Error in reading...asxdd')
                     ENDIF
                  ELSE IF ( BSKEYWORD(1:6) .EQ. 'KINBAL' ) THEN
                     IBSFLAG = 4
                  ELSE
C                    Assume it is 'SMALL   2    1....'
                     IBSFLAG = 1
                     READ (MLINE(NMLINE),'(BN,A5,12I5)',IOSTAT=IOS)
     &                    BSET,IQM(I,2),(JCO(K,I,2),K=1,IQM(I,2))
                     IF (IOS.NE.0) THEN
                       WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                       WRITE(LUPRI,*) MLINE(NMLINE)
                       CALL QUIT('Error in reading...askdjdcn')
                     ENDIF

                  END IF
 
               END IF
C
C              We know the type, call the appropriate reader
C
               IF ( IBSFLAG .EQ. 1 ) THEN
C
C                 Read explicitly given small component basis set.
C
                  KAOVEC = KBLOCK + 1 - IBLOCK(2)
                  CALL GTOINP(LUINFO,IQM(I,2),JCO(1,I,2),
     &                 NUC(IBLOCK(2),2),NRC(IBLOCK(2),2),
     &                 SEG(IBLOCK(2),2),
     &                 ALPHA(1,IBLOCK(2),2),CPRIM(1,1,IBLOCK(2),2),
     &                 CPRIMU(1,1,IBLOCK(2),2),ISGEN(IBLOCK(2)),
     &                 NBLOCK,KAOVEC,KPRIM)
                  NBLCK(I,2) = NBLOCK
                  BASREF(1,I,2) = 
     &            "Basis set typed explicitly in input file "
C
               ELSE IF ( IBSFLAG .EQ. 2 ) THEN
C
C                 Read MOLFDIR small component basis set.
C
                  KBFPR = 1
                  KBFCO = KBFPR + KPRIM
                  KLAST = KBFCO + KPRIM*KPRIM
                  IF (KLAST.GT.LWORK)
     &                 CALL STOPIT('BASINP','BFGINP',KLAST,LWORK)
                  CALL BFGINP(I,IBLOCK(2),2,WORK(KBFPR),WORK(KBFCO),
     &                IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &                NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,BASFIL,
     &                BASREF(1,1,2))
C
               ELSE IF ( IBSFLAG .EQ. 4 ) THEN
C
C                 When uncontracted we generate the small component basis
C                 using the kinetic balance relation
C                 ...or...
C                 Make uncontracted small component basis
C                 set from kinetic balance.
C                 =======================================================
C
                  CALL KINBAL(I,IBLOCK(1),IBLOCK(2),
     &                 IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &                 NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,INPTST,
     &                 CNTBAS)
                  BASREF(1,I,2) = ' Derived from large component'
C
               ENDIF
               NOORBTS    = NOORBTS.AND.IQM(I,2).EQ.0
C
C              Print small comp. basis set
C
               IF (IQM(I,2).GT.0) THEN
                 IF(INPTST) THEN
                   WRITE(LUPRI,'(3X,A)') 'Small component set:'
                   WRITE(LUPRI,'(6X,A,I5,A,12I5)')
     &                 'Max.ang.quantum no.:',(IQM(I,2)-1),
     &                 '  Blocks:',(JCO(K,I,2), K=1,IQM(I,2))
                 ENDIF
                 NHTYP     = MAX(NHTYP,IQM(I,2))
                 IF(NHTYP.GT.MXQN) GOTO 5000
               ENDIF
C
C              Check if small component basis is contracted
C              ============================================
C
               CNTBAS = .FALSE.
               DO JBLCK = IBLOCK(2), NBLCK(I,2)+IBLOCK(2)-1 
                  IF ((NUC(JBLCK,2).NE.NRC(JBLCK,2)).AND.
     &                (NRC(JBLCK,2).GT.0)) CNTBAS = .TRUE.
               END DO
               IF ( CNTBAS .AND. IBSFLAG .NE. 4 ) DIRCON(2) = .TRUE.
C
CMI        ... end of small component case ... TWOCOMP
             ENDIF

               IBLOCK(1)  = IBLOCK(1) + NBLCK(I,1)
               IBLOCK(2)  = IBLOCK(2) + NBLCK(I,2)
C              IBLOCK(3)  = IBLOCK(3) + NBLCK(I,3)
             ENDIF

           ELSE
C
C            Dalton (NOT DIRAC) case
C
             DO 300 J = 1,NSETS
C
C              Construct basis set for dummy centers if requested
C              ==================================================
C
               NQ = NINT(Q + 0.01D0)
               IF ((LCMMAX .GT. -1) .AND. (NQ .EQ. 0)) THEN
                  IQM(I,J) = LCMMAX + 1
                  KAOVEC = KBLOCK + 1 - IBLOCK(J)
                  CALL CMBAS(IQM(I,J),JCO(1,I,J),NUC(IBLOCK(J),J),
     &                 NRC(IBLOCK(J),J),SEG(IBLOCK(J),J),
     &                 ALPHA(1,IBLOCK(J),J),
     &                 CPRIM(1,1,IBLOCK(J),J),CPRIMU(1,1,IBLOCK(J),J),
     &                 NBLOCK,KAOVEC,KPRIM)
C
C              Get basis set from basis set library
C              ===================================
C
               ELSE IF ((BASIS .OR. (DOHUCKEL .AND. J .EQ. 2))) THEN
                  KAOVEC = KBLOCK + 1 - IBLOCK(2)
                  BASNAM = 'HUCKEL'
                  CALL BASLIB(IQM(I,J),JCO(1,I,J),NUC(IBLOCK(J),J),
     &                  NRC(IBLOCK(J),J),SEG(IBLOCK(J),J),
     &                  ALPHA(1,IBLOCK(J),J),
     &                  CPRIM(1,1,IBLOCK(J),J),CPRIMU(1,1,IBLOCK(J),J),
     &                  NBLOCK,KAOVEC,KPRIM,Q,Q,DSM,
     &                  UNCONT,BASNAM,BASREF(1,I,1),IPREAD)
C
               ELSE
C
C              Read basis set from file
C              ========================
C
                  KAOVEC = KBLOCK + 1 - IBLOCK(J)
                  CALL GTOINP(LUINFO,IQM(I,J),JCO(1,I,J),
     &                NUC(IBLOCK(J),J),NRC(IBLOCK(J),J),
     &                SEG(IBLOCK(J),J),
     &                ALPHA(1,IBLOCK(J),J),CPRIM(1,1,IBLOCK(J),J),
     &                CPRIMU(1,1,IBLOCK(J),J),IDUMMY,
     &                NBLOCK,KAOVEC,KPRIM)
               ENDIF
               NBLCK(I,J) = NBLOCK
               IBLOCK(J)  = IBLOCK(J) + NBLCK(I,J)
               NOORBTS    = NOORBTS.AND.IQM(I,J).EQ.0
               IF (IQM(I,J).GT.0) THEN
                 IF(INPTST) THEN
                    WRITE(LUPRI,'(3X,A5,A)') BSET,' set:'
                    WRITE(LUPRI,'(6X,A,I5,A,12I5)')
     &                   'Max.ang.quantum no.:',(IQM(I,J)-1),
     &                   '  Blocks:',(JCO(K,I,J), K=1,IQM(I,J))
                 ENDIF
C
C     Check angular momentum quantum number
C
                 NHTYP     = MAX(NHTYP,IQM(I,J))
                 IF(NHTYP.GT.MXQN) GOTO 5000
               ENDIF
 300         CONTINUE

           ENDIF ! End of switch DALTON/DIRAC
        END IF 

C       Point charge (i.e. NO ORBITALS) ?
        DO K = 1,NONT(I)
           NOORBT(K+INUC) = NOORBTS
        END DO
C
        DO ISET = 3, 4
C
C       Read basis sets for density fitting
C       ===================================
C
C       The following possibilities exist:
C
C       FTSET INTGRL     2    1    1  (IBSFLAG = 1)
C       FTSET     2    1    1 (IBSFLAG = 1)
C       FTSET DEFAULT (IBSFLAG = 2)
C       FTSET basisset (IBSFLAG = 3)
C       FTSET GEOM alpha beta N 2 1 1 (IBSFLAG = 4)
C       FTSET EVENTEMP alpha beta N 2 1 1 (IBSFLAG = 4)
c               (alias for FTSET GEOM)
C       FTSET WELLTEMP alpha beta gamma delta N 2 1 1 (IBSFLAG = 4)
C       FTSET FAMILY    N    2    1    1 (IBSFLAG = 4)
C       FTSET DUALFAMILY    N    2    1    1 (IBSFLAG = 4)
C 
C       We also distinguish between the keywords :
C       FTSET : fitset to be used for both large and small component densities
C       FTLRG : fitset to be used for only the large component density
C       FTSML : fitset to be used for only the small component density
C
        CALL RDLINE(LUINFO)
        READ (MLINE(NMLINE),'(A5,1X,A14)',IOSTAT=IOS) BSET,BSKEYWORD
        IF (IOS.NE.0) THEN
          WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
          WRITE(LUPRI,*) MLINE(NMLINE)
          CALL QUIT('Error in reading..csmxn')
        ENDIF
C
C       Initialize basis set flag
C
        IBSFLAG = 0
C
C       Check whether the fit function set is specified in input
C                 
        IF ( BSET.EQ.'FTLRG' .OR. BSET.EQ.'FTSML' ) THEN

           WRITE(LUPRI,'(//2A)') BSET,
     &     ' specified, but only FTSET is supported at the moment.'
           CALL QUIT('Unsupported option in MOLECULE input')

        ELSE IF ( BSET .NE. 'FTSET' ) THEN

C          No 'FTSET' keyword given.
C          Take standard set
C
C          Discard the line we just read.
C
           NMLINE = NMLINE - 1   
           IBSFLAG = 2     

        ELSE 
C
C          Check what kind of basis set we have
C
           IF ( BSKEYWORD(1:6) .EQ. 'INTGRL' ) THEN
              IBSFLAG = 1
              READ (MLINE(NMLINE),'(12X,12I5)',IOSTAT=IOS)
     &             IQM(I,ISET),(JCO(K,I,ISET),K=1,IQM(I,ISET))
              IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                 WRITE(LUPRI,*) MLINE(NMLINE)
                 CALL QUIT('Error in reading...ajdcn')
              ENDIF
           ELSE IF ( BSKEYWORD(1:8) .EQ. 'EXPLICIT') THEN
              IBSFLAG = 1
              READ (MLINE(NMLINE),'(14X,12I5)',IOSTAT=IOS)
     &             IQM(I,ISET),(JCO(K,I,ISET),K=1,IQM(I,ISET))
              IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                 WRITE(LUPRI,*) MLINE(NMLINE)
                 CALL QUIT('Error in reading...ajxxcn')
              ENDIF
           ELSE IF ( BSKEYWORD(1:7) .EQ. 'DEFAULT') THEN
              IBSFLAG = 2
           ELSE IF ( BSKEYWORD(1:5) .EQ. 'BASIS') THEN
              IBSFLAG = 3
              READ (MLINE(NMLINE),'(12X,A20)') BASFIL
           ELSE IF ( BSKEYWORD(1:8) .EQ. 'WELLTEMP' ) THEN
              IBSFLAG = 4
              READ (MLINE(NMLINE)(16:),*,ERR=9000)
     $             (FAMPAR(K),K=1,4),NFAMEXP(1),
     $             IQM(I,ISET), (JCO(K,I,ISET),K=1,IQM(I,ISET))
              CALL GENFAMEXP(LUINFO,0)
           ELSE IF ( BSKEYWORD(1:8) .EQ. 'EVENTEMP' ) THEN
              IBSFLAG = 4
              READ (MLINE(NMLINE)(16:),*,ERR=9000)
     $             (FAMPAR(K),K=1,2),NFAMEXP(1),
     $             IQM(I,ISET), (JCO(K,I,ISET),K=1,IQM(I,ISET))
              CALL GENFAMEXP(LUINFO,0)
           ELSE IF ( BSKEYWORD(1:4) .EQ. 'GEOM' ) THEN
              IBSFLAG = 4
              READ (MLINE(NMLINE)(11:),*,ERR=9000)
     $             (FAMPAR(K),K=1,2),NFAMEXP(1),
     $             IQM(I,ISET), (JCO(K,I,ISET),K=1,IQM(I,ISET))
              FAMPAR(3) = D0
              FAMPAR(4) = D0
              CALL GENFAMEXP(LUINFO,0)
           ELSE IF ( BSKEYWORD(1:6) .EQ. 'FAMILY') THEN
              IBSFLAG = 4
              READ (MLINE(NMLINE),'(11X,12I5)',ERR=9000)
     $        NFAMEXP(1), IQM(I,ISET), (JCO(K,I,ISET),K=1,IQM(I,ISET))
              CALL GENFAMEXP(LUINFO,1)
           ELSE IF ( BSKEYWORD(1:10) .EQ. 'DUALFAMILY') THEN
              IBSFLAG = 4
              READ (MLINE(NMLINE),'(15X,12I5)',ERR=9000)
     $             NFAMEXP(1), NFAMEXP(2),
     &             IQM(I,ISET), (JCO(K,I,ISET),K=1,IQM(I,ISET))
              CALL GENFAMEXP(LUINFO,2)
C
           ELSE
C             Assume it is 'FTSET   2    1....'
              IBSFLAG = 1
              READ (MLINE(NMLINE),'(5X,12I5)',IOSTAT=IOS)
     &             IQM(I,ISET),(JCO(K,I,ISET),K=1,IQM(I,ISET))
              IF (IOS.NE.0) THEN
                WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                WRITE(LUPRI,*) MLINE(NMLINE)
                CALL QUIT('Error in reading...kdjcn')
              ENDIF

           END IF
C
           IF ( IBSFLAG .EQ. 1 ) THEN
C
C             Read explicitly given basis set.
C
              KAOVEC = KBLOCK + 1 - IBLOCK(ISET)
              CALL GTOINP(LUINFO,IQM(I,ISET),JCO(1,I,ISET),
     &          NUC(IBLOCK(ISET),ISET),NRC(IBLOCK(ISET),ISET),
     &          SEG(IBLOCK(ISET),ISET),
     &          ALPHA(1,IBLOCK(ISET),ISET),CPRIM(1,1,IBLOCK(ISET),ISET),
     &          CPRIMU(1,1,IBLOCK(ISET),ISET),ISGEN(IBLOCK(ISET)),
     &          NBLOCK,KAOVEC,KPRIM)
              NBLCK(I,ISET) = NBLOCK
              BASREF(1,I,ISET) = 
     &        "Fit set typed explicitly in input file   "
C
           ELSE IF ( IBSFLAG .EQ. 2 ) THEN
              BASREF(1,I,ISET) = ' Default fit function set    '
C
           ELSE IF ( IBSFLAG .EQ. 3 ) THEN
C
              BASNAM = BASFIL
C
              KAOVEC = KBLOCK + 1 - IBLOCK(ISET)
              CALL BASLIB(IQM(I,ISET),JCO(1,I,ISET),
     &             NUC(IBLOCK(ISET),ISET),
     &             NRC(IBLOCK(ISET),ISET),SEG(IBLOCK(ISET),ISET),
     &             ALPHA(1,IBLOCK(ISET),ISET),
     &             CPRIM(1,1,IBLOCK(ISET),ISET),
     &             CPRIMU(1,1,IBLOCK(ISET),ISET),
     &             NBLOCK,KAOVEC,KPRIM,Q,Q,DSM,
     &             UNCONT,BASNAM,BASREF(1,I,ISET),IPREAD)
C
              NBLCK(I,ISET) = NBLOCK
           ELSE IF ( IBSFLAG .EQ. 4 ) THEN
C
C             Generate well-temperered/family basis set
C
              KAOVEC = KBLOCK + 1 - IBLOCK(ISET)
              CALL FAMBAS(LUINFO,IQM(I,ISET),JCO(1,I,ISET),
     &             NUC(IBLOCK(ISET),ISET),
     &             NRC(IBLOCK(ISET),ISET),SEG(IBLOCK(ISET),ISET),
     &             ALPHA(1,IBLOCK(ISET),ISET),
     &             CPRIM(1,1,IBLOCK(ISET),ISET),
     &             CPRIMU(1,1,IBLOCK(ISET),ISET),
     &             ISGEN(IBLOCK(ISET)),
     &             NBLOCK,KAOVEC,KPRIM)
              NBLCK(I,ISET) = NBLOCK
              BASREF(1,I,ISET) = 
     &   "Well-tempered fit set typed explicitly in input file "
C
           ENDIF
C
C          Print out fit set
C
           NOFITFS  = NOFITFS.AND.IQM(I,ISET).EQ.0
           IF (IQM(I,ISET).GT.0) THEN
              IF(INPTST) THEN
                 WRITE(LUPRI,'(3X,A)') 'Density fit set:'
                 WRITE(LUPRI,'(6X,A,I5,A,12I5)')
     &                'Max.ang.quantum no.:',(IQM(I,ISET)-1),
     &                '  Blocks:',(JCO(K,I,ISET), K=1,IQM(I,ISET))
              ENDIF
C aspg, 2006-07-03
C replacing NHTYPF for NHTYP while a complete rework of the reading
C of basis for density fitting is not done.
              NHTYPF    = MAX(NHTYPF,IQM(I,ISET))
C             Since it's the density the angular momentum may be twice
C             the maximum value for the orbitals : 2 * MXQN - 2
              IF(NHTYPF-1.GT.2*MXQN-2) GOTO 5001
C              NHTYP    = MAX(NHTYP,IQM(I,ISET))
C              IF(NHTYP -1.GT.2*MXQN-2) GOTO 5001
           ENDIF
C
C       (no need to check whether the set is contracted, so code differs
C        from analogous basis set read : this is not an error !)
C
        ENDIF
        IBLOCK(ISET)  = IBLOCK(ISET) + NBLCK(I,ISET)
C
C       End of loop over small/large fitsets.
        ENDDO
C
C
C       Read RECP parameter if exist
C
        CALL RECP_LNK_READCP(LUINFO,I,NONT(I),IPREAD)
C
C     End of loop over unique atom types
  100 CONTINUE !  DO 100 I = 1, NONTYP
C
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(6X,A,I3,3A/9X,2A/9X,2(A,I3),A)')
     &     '*  Input specifies highest orbital of atomic type ',
     &     I,' as "',SPDCAR(NHTYP-1),'".',
     &     ' Highest allowed orbital in this version: ',
     &     SPDCAR(MXQN-1),
     &     ' Increase MXQN from',MXQN,' to',NHTYP/2,' and recompile.'
        CALL QUIT('Too high angular specified in input.')
 5001 CONTINUE
        WRITE (LUPRI,'(6X,A,I3,3A/9X,2A/9X,2(A,I3),A)')
     &     '*  Input specifies highest fit function of atomic type ',
     &     I,' as "',SPDCAR(NHTYPF-1),'".',
     &     ' Highest allowed fit function in this version: ',
     &     SPDCAR(2*MXQN-2),
     &     ' Increase MXQN from',MXQN,' to',NHTYPF,' and recompile.'
        CALL QUIT('Too high angular specified in input.')
 5100 CONTINUE
        CALL QUIT('Error reading keyword LARGE.')
 5200 CONTINUE
        CALL QUIT('Error reading keyword SMALL or SBFIL.')
 9000 CONTINUE
      WRITE(LUPRI,'(A,/,A,I3,/,A)')
     $     '*** ERROR in READIN ***',
     &     'Error reading line ',NMLINE,
     &     MLINE(NMLINE)
         CALL QUIT('*** ERROR in READIN ***')
      END
C  /* Deck cntinp */
      SUBROUTINE CNTINP(LUINFO,NONTVC,Q,QEXP,ANG)
C*****************************************************************************
C
C     Read and process data about symmetry independent centers
C
C*****************************************************************************
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "consts.h"
      PARAMETER(CORMAX = 1.D5)
#if defined (SYS_HPUX)
      PARAMETER (LUERR = 7)
#else
      PARAMETER (LUERR = 0)
#endif
#include "cbirea.h"
C
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "symmet.h"
#include "dcbgen.h"
#include "argoscom.h"
C
      CHARACTER*(LEN_MLINE) TMP_LINE
      LOGICAL ANG
C
      NQ = NINT(Q)

      IF (GAUNUC .AND. NQ.GT.0) THEN
        IF(QEXP.EQ.D0) THEN
           CALL NUCSIZ(NQ,GEXP)
        ELSE
           GEXP = QEXP
        ENDIF
C
C       Print out warning if GEXP < 1.0D8
C
        WRITE(LUPRI,'(A,F8.3,A,1P,D20.10)')
     &       '  Nuclear Gaussian exponent for atom of charge',
     &       Q,' :',GEXP
        IF (GEXP .LE. 1.0D8) WRITE(LUPRI,'(A)')
     &       '  -- WARNING!!! Nuclear Gaussian exponent is small!'
      ELSE
         GEXP = D0
      END IF
      ECPNUC = NUCIND
      NUCSTR = NUCIND + 1
      IF (INPTST)
     &   WRITE(LUPRI,'(A)') '   Symmetry independent centers:'
      DO 100 N = 1,NONTVC
        NUCIND = NUCIND + 1
        IF (NUCIND .GT. MXCENT) GOTO 5000
C
C*****************************************************************************
C       Card 7:
C          NAMN - name of symmetry independent center
C          CORR - coordinates of symmetry independent center
C*****************************************************************************
C
        CALL RDLINE(LUINFO)
        READ (MLINE(NMLINE),'(A4)') NAMN(NUCIND)
C        951115-hjaaj: NAMN() must never be read with free format
C          (otherwise it will on e.g. AIX be filled with nulls:
C           e.g. ' C  ' would become 'C\NULL\NULL\NULL'; and on
C           e.g. IRIX free format read of char. var. is not defined).
#if !defined (SYS_CRAY) && !defined (VAR_NOFREE) && !defined (SYS_T3D)
        READ (MLINE(NMLINE)(5:80),*,ERR=101) (CORD(J,NUCIND),J = 1, 3)
        GO TO 104
#endif
  101   READ (MLINE(NMLINE)(1:80),'(BN,4X,3F20.0)',ERR=102)
     &        (CORD(J,NUCIND), J = 1,3)
        GO TO 104
  102   READ (MLINE(NMLINE)(1:80),'(BN,4X,3F10.0)',ERR=103)
     &        (CORD(J,NUCIND), J = 1,3)
        GO TO 104
  103   CONTINUE
            WRITE(LUPRI,'(/A,I5/A,I5,A)')
     &      ' ERROR: Unable to read Cartesian coordinates of atom no.',
     &      NUCIND,' from line',NMLINE,' in the MOLECULE input file:'
            WRITE(LUPRI,'(A)') MLINE(NMLINE)
         CALL QUIT('ERROR reading atomic coordinates in MOLECULE input')
  104   CONTINUE

        TMP_LINE = MLINE(NMLINE)
        CALL UPCASE(TMP_LINE)
        IPOS = INDEX(TMP_LINE,'ISOTOPE=')
        IF (IPOS .NE. 0) THEN
          IPOS = IPOS + 8
          READ (MLINE(NMLINE)(IPOS:),'(I3)') MASSNM
          ISOTOP(NUCIND) = MASS2ISOTOP(NQ,MASSNM)
        ELSE
          ISOTOP(NUCIND) = 1
        END IF
C
        IF(INPTST) THEN
           WRITE(LUPRI,'(6X,A4,3F20.15)') NAMN(NUCIND),
     &          (CORD(J,NUCIND), J = 1,3)
        ENDIF
        NCLINE(NUCIND) = NMLINE
        NAMEX(3*NUCIND)     = NAMN(NUCIND)//' z'
        NAMEX(3*NUCIND - 1) = NAMN(NUCIND)//' y'
        NAMEX(3*NUCIND - 2) = NAMN(NUCIND)//' x'
        DO 200 J = 1,3
           IF(ANG) CORD(J,NUCIND) = CORD(J,NUCIND)/XTANG
           IF(ABS(CORD(J,NUCIND)).GT.CORMAX) GOTO 5010
  200   CONTINUE
C
C*****************************************************************************
C       CHARGE  - charge of center
C       NOORBT  - TRUE: no orbitals on this center
C       GNUEXP  - exponent of Gaussian nuclear charge distribution
C*****************************************************************************
C
        CHARGE(NUCIND) = Q
        IZATOM(NUCIND) = NQ
C       ... if point charge this is reset outside to -1234567890
C           to tell this is a point charge /hjaaj Nov 2003
        GNUEXP(NUCIND) = GEXP
  100 CONTINUE   !  DO 100 N = 1,NONTVC
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A/A,I5)')
     &    ' Too many atomic centers: MXCENT exceed in CNTINP,',
     &    ' Current limit:',MXCENT
        CALL QUIT('*** ERROR *** MXCENT exceeded in CNTINP')
 5010 CONTINUE
        WRITE (LUPRI,'(A,1P,E12.5,A/A/A,E12.5)')
     &    ' Atomic coordinate ',CORD(J,NUCIND),
     &    ' too large in CNTINP.',
     &    ' Note: Program is unstable for large coordinate values.',
     &    ' Maximum coordinate value:',CORMAX
        CALL QUIT('*** ERROR: Atomic coordinate too large in CNTINP')
      END ! SUBROUTINE CNTINP
C  /* Deck baspro */
      SUBROUTINE BASPRO(WORK,LWORK,NSYMOP,IQM,NBLCK,JCO,NUC,
     &                  NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,
     &                  KBLOCK,KPRIM,DOOWN,blk)
C
C*****************************************************************************
C
C       Process orbital and geometry input data
C
C*****************************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "cbirea.h"
#include "dcbgen.h"
#include "argoscom.h"
#include "argoscomch.h"

      DIMENSION WORK(LWORK)
      LOGICAL SEG,DOOWN
      DIMENSION IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),
     &          NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &          SEG(KBLOCK,KSETS),
     &          ALPHA(KPRIM,KBLOCK,KSETS),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,KSETS),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS)
      integer blk
#include "ccom.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "primit.h"
#include "huckel.h"
#include "aosotr.h"
#include "gencon.h"
#include "fitshells.h"
#include "fitprimit.h"

! local arrays:
      REAL*8    CORD_DEP(3,4) ! (xyz,1:nsymop+1)

C
C*****************************************************************************
C       MULK    - bitstring of basic operations that stabilise a center
C*****************************************************************************
C
      II = 0
      DO 100 I = 1, NONTYP        ! loop over atomic types
         DO 110 J = 1, NONT(I)    ! loop over symmetry-independent centers
            II = II + 1
            II_DEP = 1
            CORD_DEP(1:3,II_DEP) = CORD(1:3,II)
            MULK = 0
            DO 140 L = 1, NSYMOP  ! loop over generators
               LL = 2**(L-1)
               II_DEP = II_DEP + 1
               DO 150 M = 1,3     ! loop over Cartesian coordinates
                  IF(IAND(LL,ISYMAX(M,1)) .NE. 0) THEN
C                 Basic operation LL moves coordinate M
                     CORD_DEP(M,II_DEP) = -CORD(M,II)
                  ELSE
                     CORD_DEP(M,II_DEP) =  CORD(M,II)
                  END IF
 150           CONTINUE
               DO K = 1, II_DEP-1
                  DIST = (CORD_DEP(1,K) - CORD_DEP(1,II_DEP))**2
     &                 + (CORD_DEP(2,K) - CORD_DEP(2,II_DEP))**2
     &                 + (CORD_DEP(3,K) - CORD_DEP(3,II_DEP))**2
                  IF (DIST .LT. 1.D-12) THEN ! not a new center
                     MULK = MULK + LL
                     EXIT
                  END IF
               END DO
 140        CONTINUE
            ISTBNU(II) = MULK
 110     CONTINUE
 100  CONTINUE      

C*****************************************************************************
C     Set various basis parameters
C*****************************************************************************
      CALL BASPAR(NSYMOP)

C*****************************************************************************
C     Process nuclear data
C*****************************************************************************
      CALL NUCPRO(WORK,LWORK)

C*****************************************************************************
C     Process orbital data for large components
C*****************************************************************************
      ISHELL = 0
      IPRIM  = 0
      IPRIMD = 0
      IORB   = 0
      IORBD  = 0

      LCOMP  = 1
      CALL ORBPRO(IQM(1,1),NBLCK(1,1),JCO(1,1,1),NUC(1,1),
     &            NRC(1,1),SEG(1,1),
     &            ALPHA(1,1,1),CPRIM(1,1,1,1),CPRIMU(1,1,1,1),
     &            KATOM,KANG,KBLOCK,KPRIM,ISHELL,
     &            IPRIM,IPRIMD,IORB,IORBD,LCOMP)
C      For relativistic RECP integrals
       CALL RECP_LNK_RDORB(IQM(1,1),JCO(1,1,1),NUC(1,1),NRC(1,1),
     &            ALPHA(1,1,1),CPRIM(1,1,1,1),CPRIMU(1,1,1,1),
     &            KATOM,KANG,KBLOCK,KPRIM)
      NLRGSH = ISHELL
      NLARGE = IORBD
      NORBL  = IORB
      NPLRG  = IPRIMD
      NPLSH  = IPRIM
CTROND      IF (DIRCON(1).AND.CNTMAT) THEN
      IF (.FALSE.) THEN
         KLPRIM = 1
         KDKWRK = KLPRIM + KPRIM
         LDKWRK = LWORK - KDKWRK
         IF (LDKWRK .LT. NPLRG*NLARGE)
     &     CALL STOPIT('READI1 ','MAKE_CNTMAT',NPLRG*NLARGE,LDKWRK)
         CALL MAKE_CNTMAT(IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &              NUC(1,1),NRC(1,1),SEG(1,1),ALPHA(1,1,1),
     &              CPRIM(1,1,1,1),WORK(KLPRIM),WORK(KDKWRK),
     &              KATOM,KANG,KBLOCK,KPRIM,NPLRG,NLARGE)
      END IF
C*****************************************************************************
C       Process orbital data for small components
C*****************************************************************************
      LCOMP = 2
      CALL ORBPRO(IQM(1,2),NBLCK(1,2),JCO(1,1,2),NUC(1,2),
     &            NRC(1,2),SEG(1,2),
     &            ALPHA(1,1,2),CPRIM(1,1,1,2),CPRIMU(1,1,1,2),
     &            KATOM,KANG,KBLOCK,KPRIM,
     &            ISHELL,IPRIM,IPRIMD,IORB,IORBD,LCOMP)

      IF (ISHELL.GT. MXSHEL) GOTO 5000
      IF (IORBD .GT. MXCORB) GOTO 5010
      IF (IPRIM .GT. MXPRIM) GOTO 5020

      IF (DOHUCKEL) THEN
         NSMLSH = ISHELL - NLRGSH
         NORBS  = IORB
         NSMALL = IORBD  - NLARGE
         NPSML  = IPRIMD - NPLRG
         NPSHEL = IPRIM  - NPLRG
C
         KMAX   = NLRGSH
         NBASIS = NLARGE
         NPBAS  = NPLRG
C
         NHUCSH = NSMLSH
         NHUCBA = NSMALL
         NHUCNP = NPSML
C
      ELSE
         KMAX   = ISHELL
         NSMLSH = KMAX - NLRGSH
         NORBS  = IORB
         NBASIS = IORBD
         IF (IPREAD .GT. 20) THEN
            WRITE(LUPRI,'(A)') 'PRIEXP:'
            WRITE(LUPRI,'(I5,F15.5)') (KK,PRIEXP(KK),KK=1,IPRIM)
         END IF
         NSMALL = NBASIS - NLARGE
         NPBAS  = IPRIMD
         NPSML   = IPRIMD - NPLRG
         NPSHEL  = IPRIM
         NPSSH   = IPRIM - NPLSH
      END IF

C*****************************************************************************
C       Process fit set data for first (LL or LL+SS) density fitset
C*****************************************************************************
      LCOMP_FIT  = 1
      ISHELL_FIT = 0
      IPRIM_FIT  = 0
      IPRIMD_FIT = 0
      IORB_FIT   = 0
      IORBD_FIT  = 0
      CALL FITPRO(IQM(1,3),NBLCK(1,3),JCO(1,1,3),NUC(1,3),
     &            NRC(1,3),SEG(1,3),
     &            ALPHA(1,1,3),CPRIM(1,1,1,3),CPRIMU(1,1,1,3),
     &            KATOM,KANG,KBLOCK,KPRIM,
     &            ISHELL_FIT,IPRIM_FIT,IPRIMD_FIT,
     &            IORB_FIT,IORBD_FIT,LCOMP_FIT)
      NLRGSH_FIT = ISHELL_FIT
      NLARGE_FIT = IORBD_FIT
      NORBL_FIT  = IORB_FIT
      NPLRG_FIT  = IPRIMD_FIT
      NPLSH_FIT  = IPRIM_FIT
C*****************************************************************************
C       Process fit set data for secondary (SS) density fitset
C*****************************************************************************
      LCOMP_FIT = 2
      CALL FITPRO(IQM(1,4),NBLCK(1,4),JCO(1,1,4),NUC(1,4),
     &            NRC(1,4),SEG(1,4),
     &            ALPHA(1,1,4),CPRIM(1,1,1,4),CPRIMU(1,1,1,4),
     &            KATOM,KANG,KBLOCK,KPRIM,
     &            ISHELL_FIT,IPRIM_FIT,IPRIMD_FIT,
     &            IORB_FIT,IORBD_FIT,LCOMP_FIT)
      IF (ISHELL_FIT.GT. MXSHEL_FIT) GOTO 5000
      IF (IORBD_FIT .GT. MXCORB_FIT) GOTO 5010
      IF (IPRIM_FIT .GT. MXPRIM_FIT) GOTO 5020
      KMAX_FIT   = ISHELL_FIT
      NSMLSH_FIT = KMAX_FIT - NLRGSH_FIT
      NBASIS_FIT = IORBD_FIT
      NSMALL_FIT = NBASIS_FIT - NLARGE_FIT
      NPBAS_FIT  = IPRIMD_FIT
      NPSML_FIT  = IPRIMD_FIT - NPLRG_FIT
      NPSHEL_FIT = IPRIM_FIT
      NPSSH_FIT  = IPRIM_FIT - NPLSH_FIT
C
C     ***************************************
C     ***** Orbital Symmetry Processing *****
C     ***************************************
C
      KKVAL = 1
      KMVAL = KKVAL + MXAQN
      KNVAL = KMVAL + MXAQN
      KIREP = KNVAL + MXAQN
      KLAST = KIREP + MXCORB
      IF (KLAST .GT. LWORK) GOTO 5030
      CALL SYMPRO(WORK(KKVAL),WORK(KMVAL),WORK(KNVAL),WORK(KIREP),
     &            DOOWN,blk)
C
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' BASPRO error, number of contracted shells      ',ISHELL,
     *  '               current maximum number (MXSHEL)  ',MXSHEL,
     *  ' Increase MXSHEL in DIRAC/include/maxorb.h and recompile'
        CALL QUIT('BASPRO: Too many contracted shells')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' BASPRO error, number of contracted basis functions',IORBD,
     *  '               current maximum number (MXCORB)     ',MXCORB,
     *  ' Increase MXCORB in DIRAC/include/maxorb.h and recompile'
        CALL QUIT('BASPRO: Too many contracted basis functions')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' BASPRO error, number of primitive basis shells   ',IPRIM,
     *  '               current maximum number (MXPRIM)    ',MXPRIM,
     *  ' Increase MXPRIM in .../include/maxorb.h and recompile'
        CALL QUIT('BASPRO: Too many primitive shells')
 5030 CONTINUE
        CALL STOPIT('BASPRO','SYMPRO',KLAST,LWORK)
 5070 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A,I6/A,I6)')
     *  ' BASPRO error, # of shells when adding Density-fitting'//
     &     ' orbitals',ISHELL,
     *  '                   current maximum number                 ',
     &     MXSHEL, ' IORBD=', IORBD, 'IPRIMD=', IPRIMD
        CALL QUIT('BASPRO: Too many shells')
 5080 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A,I6)')
     *  ' BASPRO error, # of contracted bf when adding Density-'//
     &     'fitting orbitals',IORBD,
     *  '                   current maximum number                 ',
     &     MXCORB, 'IPRIMD=', IPRIMD
        CALL QUIT('BASPRO: Too many contracted basis functions')
 5090 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' BASPRO error, # of primitive bf when adding Density-'//
     &     'fitting orbitals',IPRIMD,
     *  '                   current maximum number                 ',
     &     MXPRIM
        CALL QUIT('BASPRO: Too many primitive basis functions')
      END
C  /* Deck basout */
      SUBROUTINE BASOUT(WORK,LWORK,IQM,NBLCK,JCO,BASREF,
     &                  NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,
     &                  KBLOCK,KPRIM,HERMIT,blk)
C*****************************************************************************
C
C     Print orbital and geometry input data
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "cbirea.h"
#include "dcbgen.h"
      DIMENSION WORK(LWORK)
      LOGICAL   SEG,HERMIT,NRMPRI
      DIMENSION IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),
     &          NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &          SEG(KBLOCK,KSETS),
     &          ALPHA (KPRIM,KBLOCK,KSETS),
     &          CPRIM (KPRIM,KPRIM,KBLOCK,KSETS),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS)
      CHARACTER*(6*MXQN+3) TSTRING(MXBSETS)
      CHARACTER*80 BASREF(10,KATOM,KSETS)
      CHARACTER*31 SET_TYPE(MXBSETS)
      DIMENSION NPRIM(MXBSETS), NCONT(MXBSETS),
     &          NTPRIM(MXBSETS),NTBAS(MXBSETS)
      DIMENSION IBLOCK(MXBSETS)
      integer blk
#include "ccom.h"
#include "nuclei.h"
#include "shells.h"
#include "primit.h"
#include "symmet.h"
#include "aosotr.h"
C


      IF(DIRAC) THEN
         SET_TYPE(1) = 'large component basis functions'
         SET_TYPE(2) = 'small component basis functions'
         SET_TYPE(3) = 'primary density fit functions  '
         SET_TYPE(4) = 'secondary density fit functions'
         SET_TYPE(5) = 'extended Huckel basis functions'
C
      ELSE
         SET_TYPE(1) = 'basis functions'
         SET_TYPE(2) = 'extended Huckel basis functions'
         SET_TYPE(3) = 'auxiliary density fit basis functions'
      END IF
C
C     ******************************************
C     ***** Printing of basis information  *****
C     ******************************************
C
      IF (NBASIS .LE. 0) IPREAD = MAX(1,IPREAD)
      IF (IPREAD .GT. 0) THEN
         CALL HEADER('Atoms and basis sets',1)
         WRITE (LUPRI,'(A,I5 )') '  Number of atom types :',NONTYP
         WRITE (LUPRI,'(A,I5/)') '  Total number of atoms:',NUCDEP
      END IF
C
      IBLOCK(1:KSETS) = 1
      ICENT  = 0
      IF (IPREAD .GT. 0) THEN
         WRITE (LUPRI,'(2X,A)')
     &        'label    atoms   charge   prim    cont     basis   '
         WRITE (LUPRI,'(2X,70A1)') ('-',I=1,70)
      END IF
      NCHTOT = 0
      DO ISET = 1,KSETS
        NTPRIM(ISET) = 0
        NTBAS(ISET)  = 0
      ENDDO
      NBAS_ATOMTYPE(1:NONTYP,0:mxbsets) = 0
      DO 100 I = 1,NONTYP
C
         NUCLEI_IN_TYPE = 0
         DO 20 N = 1,NONT(I)
            ICENT  = ICENT + 1
            NCHARG = IZATOM(ICENT)
            NCHTOT = NCHTOT + NUCDEG(ICENT)*NCHARG
            NUCLEI_IN_TYPE = NUCLEI_IN_TYPE + NUCDEG(ICENT)
   20    CONTINUE
C
C        Convert information to a string that can be understood by the innocent user
C
         DO ISET = 1, KSETS
            TSTRING(ISET) = ' '
            IF (IQM(I,ISET).GT.0) THEN
               CALL BASTYP(IQM(I,ISET),JCO(1,I,ISET),
     &                     NRC(IBLOCK(ISET),ISET),
     &                     NUC(IBLOCK(ISET),ISET),
     &                     NPRIM(ISET),NCONT(ISET),TSTRING(ISET))
            ELSE
               TSTRING(ISET) = 'No '//SET_TYPE(ISET)//' attached'
               NPRIM(ISET) = 0
               NCONT(ISET) = 0
            END IF
            NTPRIM(ISET) = NTPRIM(ISET)+NUCLEI_IN_TYPE*NPRIM(ISET)
            NTBAS(ISET)  = NTBAS(ISET) +NUCLEI_IN_TYPE*NCONT(ISET)
            IBLOCK(ISET) = IBLOCK(ISET) + NBLCK(I,ISET)
            IF (ISET.EQ.1 .OR. ISET.EQ.2) THEN
               NBAS_ATOMTYPE(I,ISET) = NCONT(ISET)
               NBAS_ATOMTYPE(I,0) = NBAS_ATOMTYPE(I,0) + NCONT(ISET)
            END IF
         END DO

         IF (IPREAD .GT. 0) THEN
            IF(DIRAC) THEN
               WRITE (LUPRI,'(2X,A4,1X,4I8,6X,2A)')
     &            NAMN(ICENT),NUCLEI_IN_TYPE,NCHARG,
     &            NPRIM(1),NCONT(1),'L  - ',TSTRING(1)
               IF(.NOT.TWOCOMPBSS) THEN
                 WRITE(LUPRI,'(23X,2I8,6X,2A)')
     &               NPRIM(2),NCONT(2),'S  - ',TSTRING(2)
               ENDIF
               IF (NPRIM(3).GT.0)
     &            WRITE(LUPRI,'(23X,2I8,6X,2A)')
     &               NPRIM(3),NCONT(3),'FL - ',TSTRING(3)
               IF (NPRIM(4).GT.0)
     &            WRITE(LUPRI,'(23X,2I8,6X,2A)')
     &               NPRIM(4),NCONT(4),'FS - ',TSTRING(4)
            ELSE
               WRITE (LUPRI,'(2X,A4,1X,4I8,6X,A)') NAMN(ICENT),
     &              NUCLEI_IN_TYPE,NCHARG,NPRIM(1),NCONT(1),TSTRING(1)
               WRITE(LUPRI,'(23X,2I8,6X,2A)')
     &               NPRIM(3),NCONT(2),'Aux - ',TSTRING(2)
               WRITE(LUPRI,'(23X,2I8,6X,2A)')
     &               NPRIM(3),NCONT(3),'Fit - ',TSTRING(3)
            ENDIF
            flush(lupri)
         END IF
  100 CONTINUE ! I = 1, NONTYP
      IF (IPREAD .GT. 0) THEN
         IF(DIRAC) THEN
            CALL PRSYMB(LUPRI,'-',70,2)
            WRITE(LUPRI,'(23X,2I8,3X,A)')
     &            NTPRIM(1),NTBAS(1),'   L  - large components'
            IF(.NOT.TWOCOMP) THEN
              WRITE(LUPRI,'(23X,2I8,3X,A)')
     &            NTPRIM(2),NTBAS(2),'   S  - small components'
            ENDIF
            IF (KSETS.EQ.3) THEN
               WRITE(LUPRI,'(23X,2I8,3X,A)')
     &            NPLRG_FIT,NLARGE_FIT, '   FL - global components fit'
            ELSEIF (KSETS.EQ.4) THEN
               WRITE(LUPRI,'(23X,2I8,3X,A)')
     &            NPLRG_FIT,NLARGE_FIT, '   FL - large components fit'
               WRITE(LUPRI,'(23X,2I8,3X,A)')
     &            NPSML_FIT,NSMALL_FIT, '   FS - small components fit'
            ENDIF
C
         ENDIF
         CALL PRSYMB(LUPRI,'-',70,2)
C     
         NTPBAS  = NPBAS  + NPBAS_FIT 
         NTBASIS = NBASIS + NBASIS_FIT
         WRITE (LUPRI,'(2X,A,I7,3I8,6X,/)')
     &         'total:',NUCDEP,NCHTOT,NTPBAS,NTBASIS
C
         IF (NHTYP .GT. 2) THEN
            IF(DOCART) THEN
               WRITE (LUPRI,'(2X,A)') 'Cartesian basis used.'
            ELSE
               WRITE (LUPRI,'(2X,A)') 'Spherical harmonic basis used.'
            ENDIF
         END IF
         WRITE (LUPRI,'(1X,A,1P,D10.2)') 
     &    ' Threshold for integrals (to be written to file):',THRS
      END IF

      IF (NBASIS.LE.0) THEN
         CALL QUIT('Fatal error: no basis functions!')
      END IF

CMI   ... prepare the MNF.INP basis set input file for AMFI
      CALL MNFOUT(IQM(1,1),NBLCK(1,1),
     &           JCO(1,1,1),NUC(1,1),NRC(1,1),SEG(1,1),
     &           KATOM,KANG,KBLOCK,
     &           KPRIM,CPRIMU)

#ifdef MOD_CAP
CMI   ... prepare the basis set input files for the CAP program
      CALL CAPOUT(IPREAD,IQM,NBLCK,
     &           JCO,NUC,NRC,SEG,
     &           KATOM,KANG,KSETS,KBLOCK,
     &           KPRIM,CPRIMU,WORK,LWORK)
#endif
C
C     Print the references to the basis sets
C

      IF (IPREAD .GT. 0) CALL PRBASREF(LUPRI,NONTYP,BASREF)

C
C     ******************************************
C     ***** Printing of atomic coordinates *****
C     ******************************************
C
      CALL PRICAR(WORK,LWORK)
C
 1010 FORMAT(I14,' large AO-blocks of ',A1,' type')
 1020 FORMAT(I14,' small AO-blocks of ',A1,' type')
      IF (IPREAD .GE. 2) THEN
         KIPCON = 1
         KLAST  = KIPCON + KPRIM
         IF (KLAST.GT.LWORK) CALL STOPIT('BASOUT','ORBOUT',KLAST,LWORK)
C*****************************************************************************
C        Print information about large component basis
C*****************************************************************************
         IF(NLARGE.GT.0) THEN
            IF (DIRAC) CALL HEADER('Large Components Basis',-1)
            NRMPRI = .FALSE.
            CALL ORBOUT(1,IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &                  NUC(1,1),NRC(1,1),SEG(1,1),
     &                  WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                  ALPHA,CPRIMU,NRMPRI)
            IF (IPREAD .GT. 2) THEN
               NRMPRI = .TRUE.
               CALL ORBOUT(1,IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &                     NUC(1,1),NRC(1,1),SEG(1,1),
     &                     WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                     ALPHA(1,1,1),CPRIMU,NRMPRI)
            END IF
         END IF
C*****************************************************************************
C        Print information about small component basis
C*****************************************************************************
         IF (NBASIS.GT.NLARGE) THEN
            CALL HEADER('Small Components Basis',-1)
            NRMPRI = .FALSE.
            CALL ORBOUT(2,IQM(1,2),NBLCK(1,2),JCO(1,1,2),
     &                  NUC(1,2),NRC(1,2),SEG(1,2),
     &                  WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                  ALPHA(1,1,2),CPRIMU(1,1,1,2),NRMPRI)
            IF (IPREAD .GT. 2) THEN
               NRMPRI = .TRUE.
               CALL ORBOUT(2,IQM(1,2),NBLCK(1,2),JCO(1,1,2),
     &                     NUC(1,2),NRC(1,2),SEG(1,2),
     &                     WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                     ALPHA(1,1,2),CPRIMU(1,1,1,2),NRMPRI)
            END IF
         END IF
C 
C aspg, 2006-06-13
C adding printout for the fit functions. when things are stable this
C and the calls to orbout above should be put in a separate function.
C
         IF (NLARGE_FIT.GT.0) THEN
            CALL HEADER('Fit Functions, Global/Large Component',-1)
            NRMPRI = .FALSE.
            CALL ORBOUT(3,IQM(1,3),NBLCK(1,3),JCO(1,1,3),
     &               NUC(1,3),NRC(1,3),SEG(1,3),
     &               WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &               ALPHA(1,1,3),CPRIMU(1,1,1,3),NRMPRI)
            IF (IPREAD .GT. 2) THEN
               NRMPRI = .TRUE.
              CALL ORBOUT(3,IQM(1,3),NBLCK(1,3),JCO(1,1,3),
     &                  NUC(1,3),NRC(1,3),SEG(1,3),
     &                  WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                  ALPHA(1,1,3),CPRIMU(1,1,1,3),NRMPRI)
            END IF
            IF (NBASIS_FIT.GT.NLARGE_FIT) THEN 
               CALL HEADER('Fit Functions, Small Component',-1)
               NRMPRI = .FALSE.
              CALL ORBOUT(4,IQM(1,4),NBLCK(1,4),JCO(1,1,4),
     &                  NUC(1,4),NRC(1,4),SEG(1,4),
     &                  WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                  ALPHA(1,1,4),CPRIMU(1,1,1,4),NRMPRI)
               IF (IPREAD .GT. 2) THEN
                  NRMPRI = .TRUE.
              CALL ORBOUT(4,IQM(1,4),NBLCK(1,4),JCO(1,1,4),
     &                     NUC(1,4),NRC(1,4),SEG(1,4),
     &                     WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                     ALPHA(1,1,4),CPRIMU(1,1,1,4),NRMPRI)
               END IF
            ENDIF
         ENDIF
C
      END IF
C
C     Symmetry basis:
C
      IF (MAXREP.GT.0.AND.IPREAD.GT.2) CALL SYMOUT
C
C     Symmetrized magnetic moments
C
CTROND      IF (HERMIT) CALL MAGCOR(IPREAD)
      CALL MAGCOR(IPREAD)
C
C     ***************************
C     ***** Determine IAOAO *****
C     ***************************
C
      CALL AOTOAO(IAOAO,IPREAD)
C
C     **************************
C     **** Determine ICNTAO ****
C     **************************
C
      CALL CNTAO(IPREAD)
C
      RETURN
      END ! subroutine BASOUT
C  /* Deck bastyp */
      SUBROUTINE BASTYP(IQM,JCO,NRC,NUC,NPRIM,NCONT,TSTRING)
C*****************************************************************************
C
C     Generate string : Uncontracted/contracted basis functions
C
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
      CHARACTER TSTRING*(6*MXQN+3)
      CHARACTER SPDCAR*1 ! type declaration for external function
      DIMENSION JCO(IQM),NRC(*),NUC(*)
#include "ccom.h"
#include "gencon.h"
C*****************************************************************************
      IND = 1
      TSTRING(IND:IND) = '['
C*****************************************************************************
C * P R I M I T I V E S:
C*****************************************************************************
      JBLOCK = 0
      NPRIM = 0
      DO 100 I = 1,IQM
         IPRIM  = 0
         ISHELL = 0
         NCOMP = KHK(I)
C
         DO 200 J = 1,JCO(I)
            JBLOCK = JBLOCK + 1
            IPRIM  = IPRIM  + NUC(JBLOCK)
            NPRIM  = NPRIM  + NCOMP*NUC(JBLOCK)
  200    CONTINUE
C
         IF (IPRIM.GT.0) THEN
            IF (IPRIM.GE.10) THEN
               NDIG = INT(LOG10(FLOAT(IPRIM))+1)
               ITEN = 10**NDIG
               DO 300 K = NDIG,1,-1
                  ITEN = ITEN/10
                  IND = IND + 1
                  IDIG = ICHAR('0') + INT(IPRIM/ITEN)
                  TSTRING(IND:IND) = CHAR(IDIG)
                  IPRIM=MOD(IPRIM,ITEN)
  300          CONTINUE
            ELSE
               IND = IND + 1
               IDIG = ICHAR('0') + IPRIM
               TSTRING(IND:IND) = CHAR(IDIG)
            ENDIF
            IND = IND + 1
            TSTRING(IND:IND) = SPDCAR(I-1)
         END IF
  100 CONTINUE
C*****************************************************************************
      IND = IND + 1
      TSTRING(IND:IND) = '|'
C*****************************************************************************
C * S H E L L S:
C*****************************************************************************
      JBLOCK = 0
      NCONT  = 0
      DO 400 I = 1,IQM
         ISHELL = 0
         NCOMP = KHK(I)
         DO 500 J = 1,JCO(I)
            JBLOCK = JBLOCK + 1
            ISHELL = ISHELL + NRC(JBLOCK)
            NCONT  = NCONT  + NCOMP*NRC(JBLOCK)
  500    CONTINUE
         IF(ISHELL.GT.0) THEN
           IF (ISHELL.GE.10) THEN
              NDIG = INT(LOG10(FLOAT(ISHELL))+1)
              ITEN = 10**NDIG
              DO 600 K = NDIG,1,-1
                 ITEN = ITEN/10
                 IND = IND + 1
                 IDIG = ICHAR('0') + INT(ISHELL/ITEN)
                 TSTRING(IND:IND) = CHAR(IDIG)
                 ISHELL=MOD(ISHELL,ITEN)
  600         CONTINUE
           ELSE
              IND = IND + 1
              IDIG = ICHAR('0') + ISHELL
              TSTRING(IND:IND) = CHAR(IDIG)
           ENDIF
           IND = IND + 1
           TSTRING(IND:IND) = SPDCAR(I-1)
         ENDIF
  400 CONTINUE
C*****************************************************************************
      IND = IND + 1
      TSTRING(IND:IND) = ']'
      RETURN
      END
C  /* Deck gtoinp */
      SUBROUTINE GTOINP(LUINFO,IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,
     &                  CPRIMU,ISGEN,NBLOCK,KAOVEC,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (DSM=1.0D-30)
      CHARACTER*1 FRMT
      CHARACTER*8 FMT1
      CHARACTER*12 FMT2
C
C
      LOGICAL SEG, FREE, SEGIJ
      DIMENSION JCO(IQM),NUC(KAOVEC),NRC(KAOVEC),SEG(KAOVEC),
     &          ALPHA(KPRIM,KAOVEC),ISGEN(KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC),CPRIMU(KPRIM,KPRIM,KAOVEC)
! dcbham.h : ECPCALC
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "gencon.h"
C
      IF(DIRAC) THEN
C
C     DIRAC allows user to split blocks !
C
        NBLOCK = 0
        DO 100 I = 1, IQM
          DO 200 J = 1, JCO(I)
            NBLOCK = NBLOCK + 1
            IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
            CALL DZERO(ALPHA (1,  NBLOCK),KPRIM)
            CALL DZERO(CPRIMU(1,1,NBLOCK),KPRIM*KPRIM)
            CALL DZERO(CPRIM (1,1,NBLOCK),KPRIM*KPRIM)
C
C           Card 8:
C
C           NUC - number of uncontracted shells in AO-block
C           NRC - number of contracted shells in AO-block
C           ISGEN - bit packed information on generation of small
C                   components using the kinetic balance relation
C             1 - downwards: small L+1 from large L
C             2 - upwards  : small L-1 from large L
C             0 - if not specified same as 3 (both down and up)
C
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),'(BN,A1,I4,2I5)',IOSTAT=IOS)
     &           FRMT,NUCIJ,NRCIJ,ISGEN(NBLOCK)
            IF (IOS.NE.0) THEN
               WRITE(LUPRI,*) 'Error in reading .mol input, line',NMLINE
               WRITE(LUPRI,*) 'Read MLINE(NMLINE)=',MLINE(NMLINE)
               IF (LUERR .NE. LUPRI) THEN
               WRITE(LUERR,*) 'Error in reading .mol input, line',NMLINE
               WRITE(LUERR,*) 'Read MLINE(NMLINE)=',MLINE(NMLINE)
               END IF
               CALL QUIT('Error in reading BN,A,I4,2I5 ... ')
            ENDIF

            IF (INPTST) THEN
              WRITE(LUPRI,'(3X,A,I2,A,I4)') 'L= ',(I-1),'  Block: ',J
              WRITE(LUPRI,'(6X,A,I5)')
     $             'Primitives:  ',NUCIJ,
     &             'Shells    :  ',NRCIJ,
     $             'KinBal    :  ',ISGEN(NBLOCK)
            END IF
            SEGIJ = NRCIJ .GE. 0
            NUCIJ = ABS(NUCIJ)
            NRCIJ = ABS(NRCIJ)
            IF (NUCIJ.GT.KPRIM) GOTO 5010
            IF (NRCIJ.GT.KPRIM) GOTO 5020
C
C           Read in exponents and contraction coefficients
C
            CALL ACPORB(LUINFO,FRMT,NUCIJ,NRCIJ,ALPHA(1,NBLOCK),
     &                  CPRIMU(1,1,NBLOCK),KPRIM)
C
C           Identify segmented contractions
C
            CALL SEGORB(SEGIJ,NUCIJ,NRCIJ,CPRIMU(1,1,NBLOCK),KPRIM,DSM)
            SEG(NBLOCK) = SEGIJ
C
C           Reorder primitive orbitals
C
            if (.not. ECPCALC) then
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUCIJ,NRCIJ,
     &                  SEG(NBLOCK),KPRIM,DSM)
            endif
C
C           Normalize orbitals
C
            CALL NRMORB(I,NRCIJ,NUCIJ,ALPHA(1,NBLOCK),
     &                  CPRIM(1,1,NBLOCK),CPRIMU(1,1,NBLOCK),
     &                  KPRIM,NBLOCK)
C
            NUC(NBLOCK) = NUCIJ
            NRC(NBLOCK) = NRCIJ
  200      CONTINUE
  100   CONTINUE
      ELSE
      NBLOCK = 0
      DO 300 I = 1, IQM
         NBLOCK = NBLOCK + 1
         IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
         CALL DZERO(ALPHA (1,  NBLOCK),KPRIM)
         CALL DZERO(CPRIMU(1,1,NBLOCK),KPRIM*KPRIM)
         CALL DZERO(CPRIM (1,1,NBLOCK),KPRIM*KPRIM)
C
         NUCI = 0
         NRCI = 0
         SEGIJ  = .TRUE.
         DO 400 J = 1, JCO(I)
C
C           Card 8:
C
C           NUC - number of uncontracted shells in AO-block
C           NRC - number of contracted shells in AO-block
C
            CALL RDLINE(LUINFO)
            IF (J .EQ. 1) THEN
               READ (MLINE(NMLINE),'(BN,A,I4,1I5)',IOSTAT=IOS)
     &              FRMT,NUCIJ,NRCIJ
               IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                 WRITE(LUPRI,*) MLINE(NMLINE)
                 CALL QUIT('Error in reading...ddxm23jdcn')
               ENDIF
            ELSE
               READ (MLINE(NMLINE),'(BN,A,I4,1I5)',IOSTAT=IOS) 
     &              FRMT,NUCIJ,NRCIJ
               IF (IOS.NE.0) THEN
                 WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                 WRITE(LUPRI,*) MLINE(NMLINE)
                 CALL QUIT('Error in reading...d123s0dm3jdcn')
               ENDIF

            END IF
            IF (INPTST) THEN
               WRITE(LUPRI,'(3X,A,I2,A,I4)') 'L= ',(I-1),'  Block: ',J
               WRITE(LUPRI,'(6X,A,I5)')
     $              'Primitives:  ',NUCIJ,
     &              'Shells    :  ',NRCIJ
            END IF
            SEGIJ = SEGIJ .AND. (NRCIJ .GE. 0)
            NUCIJ = ABS(NUCIJ)
            NRCIJ = ABS(NRCIJ)
            IF (NUCI + NUCIJ.GT.KPRIM) GOTO 5010
            IF (NRCI + NRCIJ.GT.KPRIM) GOTO 5020
C
C           Read in exponents and contraction coefficients
C
            CALL ACPORB(LUINFO,FRMT,NUCIJ,NRCIJ,ALPHA(NUCI+1,NBLOCK),
     &                  CPRIMU(NUCI+1,NRCI+1,NBLOCK),KPRIM)
C
            NUCI = NUCI + NUCIJ
            NRCI = NRCI + NRCIJ
  400    CONTINUE
         NUC(NBLOCK) = NUCI
         NRC(NBLOCK) = NRCI
C
         JCO(I) = 1
C
C        Identify segmented contractions
C
         CALL SEGORB(SEGIJ,NUCI,NRCI,CPRIMU(1,1,NBLOCK),KPRIM,DSM)
         SEG(NBLOCK) = SEGIJ
C
C        Reorder primitive orbitals
C
         if (.not. ECPCALC) then
         CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUCI,NRCI,
     &               SEG(NBLOCK),KPRIM,DSM)
         endif
C
C        Normalize orbitals
C
         CALL NRMORB(I,NRCI,NUCI,ALPHA(1,NBLOCK),CPRIM(1,1,NBLOCK),
     &               CPRIMU(1,1,NBLOCK),KPRIM,NBLOCK)
C
  300 CONTINUE
      ENDIF
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' >>> GTOINP error, number of AO-blocks            ',NBLOCK,
     *  '                   number of blocks not used yet  ',KAOVEC,
     *  '                ** Increase WRKMEM and try again!'
        CALL QUIT('Too many AO-blocks')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error,'//
     *  'number of primitives per block      ',NUCIJ,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error, number of contracted functions      ',
     *  NRCIJ,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many contracted functions')
      END
C  /* Deck gtoinp */
      SUBROUTINE FAMBAS(LUINFO,IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,
     &                  CPRIMU,ISGEN,NBLOCK,KAOVEC,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (DSM=1.0D-30)
      CHARACTER*1 FRMT
      CHARACTER*8 FMT1
      CHARACTER*12 FMT2
C
C
      LOGICAL SEG, FREE, SEGIJ
      DIMENSION JCO(IQM),NUC(KAOVEC),NRC(KAOVEC),SEG(KAOVEC),
     &          ALPHA(KPRIM,KAOVEC),ISGEN(KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC),CPRIMU(KPRIM,KPRIM,KAOVEC)
! dcbham.h : ECPCALC
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "gencon.h"
C
      NBLOCK = 0
      DO 100 I = 1, IQM
         DO 200 J = 1, JCO(I)
            NBLOCK = NBLOCK + 1
            IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
            CALL DZERO(ALPHA (1,  NBLOCK),KPRIM)
            CALL DZERO(CPRIMU(1,1,NBLOCK),KPRIM*KPRIM)
            CALL DZERO(CPRIM (1,1,NBLOCK),KPRIM*KPRIM)
C
            ISGEN (NBLOCK) = 0
C
C           Get exponents
C
            IC = MOD(I + 1, 2) + 1
            CALL FAMCOPY(LUINFO,NUCIJ,NRCIJ,ALPHA(1,NBLOCK),
     &           CPRIMU(1,1,NBLOCK),KPRIM,IC)
C
            IF (INPTST) THEN
               WRITE(LUPRI,'(3X,A,I2,A,I4)') 'L= ',(I-1),'  Block: ',J
               WRITE(LUPRI,'(6X,A,I5)')
     $              'Primitives:  ',NUCIJ,
     &              'Shells    :  ',NRCIJ,
     $              'KinBal    :  ',ISGEN(NBLOCK)
            END IF
            SEGIJ = NRCIJ .GE. 0
            NUCIJ = ABS(NUCIJ)
            NRCIJ = ABS(NRCIJ)
            IF (NUCIJ.GT.KPRIM) GOTO 5010
            IF (NRCIJ.GT.KPRIM) GOTO 5020
C
C
C           Identify segmented contractions
C
            CALL SEGORB(SEGIJ,NUCIJ,NRCIJ,CPRIMU(1,1,NBLOCK),KPRIM,DSM)
            SEG(NBLOCK) = SEGIJ
C
C           Reorder primitive orbitals
C
            if (.not. ECPCALC) then
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUCIJ,NRCIJ,
     &           SEG(NBLOCK),KPRIM,DSM)
            endif
C
C           Normalize orbitals
C
            CALL NRMORB(I,NRCIJ,NUCIJ,ALPHA(1,NBLOCK),
     &           CPRIM(1,1,NBLOCK),
     &           CPRIMU(1,1,NBLOCK),KPRIM,NBLOCK)
C
            NUC(NBLOCK) = NUCIJ
            NRC(NBLOCK) = NRCIJ
  200      CONTINUE
  100   CONTINUE
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6,A)')
     *  ' >>> FAMBAS error, number of AO-blocks            ',NBLOCK,
     *  '                   number of blocks not used yet  ',KAOVEC,
     *  '                ** Increase WRKMEM and try again!'
        CALL QUIT('Too many AO-blocks')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> FAMBAS error,'//
     *  'number of primitives per block      ',NUCIJ,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> FAMBAS error, number of contracted functions      ',
     *  NRCIJ,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many contracted functions')
      END
C  /* Deck famcopy */
      SUBROUTINE FAMCOPY(LUINFO,NUCIJ,NRCIJ,ALPHA,CPRIMU,KPRIM,IC)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D1 = 1.0D0)
      DIMENSION ALPHA(KPRIM), CPRIMU(KPRIM,KPRIM)
#include "cbirea.h"
      DIMENSION ILIST(MXFAMEXP),TMP(MXFAMEXP+2)
#include "molinp.h"
C
C     Card 9:
C
C     ALPHA  - exponent of primitive
C     CPRIMU - unnormalized contraction coefficients
C
C
C     NUMLST needs an work array for dimension NFALPHA.
C
C     Read string in format as described in SUBR NUMLST
C     Eg.: 1..20
C
      CALL RDLINE(LUINFO)
      NUCIJ = 1
      CALL NUMLST(MLINE(NMLINE)(1:72),
     &        ILIST,NFAMEXP(IC),1,NFAMEXP(IC),
     &     1,NUCIJ)
      NRCIJ = NUCIJ
      DO 300 L = 1, NUCIJ
         ALPHA(L) = FAMEXP(ILIST(L), IC)
         CPRIMU(L,L) = D1
  300 CONTINUE
C
C     Test print
C
      IF (INPTST) THEN
         WRITE(LUPRI,'(A,I4,A,I4,A)')
     &   'Contraction matrix (',NUCIJ,'x',NRCIJ,') :' 
         DO 400 L = 1, NUCIJ
            WRITE (LUPRI,'(E16.10)') ALPHA(L)
            WRITE (LUPRI,'(6(2X,F10.8))') (CPRIMU(L,M),M=1,NRCIJ)
  400    CONTINUE
      ENDIF
C
      RETURN
      END
C  /* Deck cmbas */
      SUBROUTINE CMBAS(IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,
     &                  CPRIMU,NBLOCK,KAOVEC,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
      LOGICAL SEG, FREE, SEGIJ
      DIMENSION JCO(IQM),NUC(KAOVEC),NRC(KAOVEC),SEG(KAOVEC),
     &          ALPHA(KPRIM,KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC),CPRIMU(KPRIM,KPRIM,KAOVEC),
     &          BASA(0:4), BASB(0:4)
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
      DATA (BASA(I),I=0,4) /0.584342, 0.452615, 0.382362, 0.337027,
     &                      0.304679/
      DATA (BASB(I),I=0,4) /0.424483, 0.309805, 0.251333, 0.215013,
     &                      0.189944/
C
      NBLOCK = 0
      DO 100 I = 1, IQM
         NBLOCK = NBLOCK + 1
         IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
         CALL DZERO(ALPHA (1,  NBLOCK),KPRIM)
         CALL DZERO(CPRIMU(1,1,NBLOCK),KPRIM*KPRIM)
         CALL DZERO(CPRIM (1,1,NBLOCK),KPRIM*KPRIM)
C
         NUCI = 0
         NRCI = 0
         SEGIJ  = .TRUE.
         DO 110 J = NCMSTR, NCMEND
            NUCI = NUCI + 1
            NRCI = NRCI + 1
            IF (NUCI.GT.KPRIM) GOTO 5010
            IF (NRCI.GT.KPRIM) GOTO 5020
            ALPHA(NUCI,NBLOCK) = ((ZCMVAL/J)/
     &                           ((BASA(I-1)*J)/2+BASB(I-1)))**2
            CPRIMU(NUCI,NRCI,NBLOCK) = 1.0D0
  110    CONTINUE
         NUC(NBLOCK) = NUCI
         NRC(NBLOCK) = NRCI
C
         JCO(I) = 1
         SEG(NBLOCK) = SEGIJ
C
C        Normalize orbitals
C
         CALL NRMORB(I,NRCI,NUCI,ALPHA(1,NBLOCK),CPRIM(1,1,NBLOCK),
     &               CPRIMU(1,1,NBLOCK),KPRIM,NBLOCK)
C
  100 CONTINUE
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' >>> CMBAS  error, number of AO-blocks            ',NBLOCK,
     *  '                   number of blocks not used yet  ',KAOVEC,
     *  '                ** Increase WRKMEM and try again!'
        CALL QUIT('Too many AO-blocks')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> CMBAS  error,'//
     *  'number of primitives per block      ',NUCI,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> CMBAS  error, number of contracted functions      ',
     *  NRCI,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many contracted functions')
      END
#ifdef NOT_USED
C old version of GTOINP ??? -hjaaj July 2000
C  /* Deck orbinx */
      SUBROUTINE ORBINX(LUINFO,IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,
     &                  CPRIMU,NBLOCK,KAOVEC,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
! dcbham.h : ECPCALC
#include "dcbham.h"
      PARAMETER (DSM=1.0D-30)
      CHARACTER*1 FRMT
      CHARACTER*8 FMT1
      CHARACTER*12 FMT2
      LOGICAL SEG, FREE
      DIMENSION JCO(IQM),NUC(KAOVEC),NRC(KAOVEC),SEG(KAOVEC),
     &          ALPHA(KPRIM,KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC),CPRIMU(KPRIM,KPRIM,KAOVEC)
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
C
      NBLOCK = 0
      DO 100 I=1,IQM
         DO 200 J = 1, JCO(I)
            NBLOCK = NBLOCK + 1
            IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
C           Card 8:
C
C           NUC - number of uncontracted shells in AO-block
C           NRC - number of contracted shells in AO-block
C
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),'(BN,A,I4,I5)',IOSTAT=IOS)
     &           FRMT,NUC(NBLOCK),NRC(NBLOCK)
            IF (IOS.NE.0) THEN
              WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
              WRITE(LUPRI,*) MLINE(NMLINE)
              CALL QUIT('Error in reading...dd283cn')
            ENDIF

            IF (INPTST) THEN
               WRITE(LUPRI,'(3X,A,I2,A,I4)') 'L= ',(I-1),'  Block: ',J
               WRITE(LUPRI,'(6X,A,I5)') 'Primitives:  ',NUC(NBLOCK),
     &                                  'Shells    :  ',NRC(NBLOCK)
            END IF
            SEG(NBLOCK) = NRC(NBLOCK).GE.0
            NUC(NBLOCK) = ABS(NUC(NBLOCK))
            NRC(NBLOCK) = ABS(NRC(NBLOCK))
            NUCIJ = NUC(NBLOCK)
            NRCIJ = NRC(NBLOCK)
            IF (NUCIJ.GT.KPRIM) GOTO 5010
            IF (NRCIJ.GT.KPRIM) GOTO 5020
C
C           Read in exponents and contraction coefficients
C
            CALL ACPORB(LUINFO,FRMT,NUCIJ,NRCIJ,ALPHA(1,NBLOCK),
     &                  CPRIMU(1,1,NBLOCK),KPRIM)
C
C           Identify segmented contractions
C
            CALL SEGORB(SEG(NBLOCK),NUCIJ,NRCIJ,CPRIMU(1,1,NBLOCK),
     &                  KPRIM,DSM)
C
C           Reorder primitive orbitals
C
            if (.not. ECPCALC) then
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUCIJ,NRCIJ,
     &                  SEG(NBLOCK),KPRIM,DSM)
            endif
C
C           Normalize orbitals
C
            CALL NRMORB(I,NRCIJ,NUCIJ,ALPHA(1,NBLOCK),CPRIM(1,1,NBLOCK),
     &                  CPRIMU(1,1,NBLOCK),KPRIM,NBLOCK)
C
  200    CONTINUE
  100 CONTINUE
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' >>> GTOINP error, number of AO-blocks            ',NBLOCK,
     *  '                   number of blocks not used yet  ',KAOVEC,
     *  '                ** Increase WRKMEM and try again!'
        CALL QUIT('Too many AO-blocks')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error,'//
     *  'number of primitives per block      ',NUC(NBLOCK),
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error, number of contracted functions      ',
     *  NRC(NBLOCK),
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many contracted functions')
      END
#endif /* NOT_USED */
C  /* Deck acporb */
      SUBROUTINE ACPORB(LUINFO,FRMT,NUCIJ,NRCIJ,ALPHA,CPRIMU,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0)
      CHARACTER*1 FRMT
      CHARACTER*8 FMT1
      CHARACTER*12 FMT2
      CHARACTER*7 FMT3
      LOGICAL FREE
      DIMENSION ALPHA(KPRIM), CPRIMU(KPRIM,KPRIM)
#include "cbirea.h"
#include "molinp.h"
C
C     Card 9:
C
C     ALPHA  - exponent of primitive
C     CPRIMU - unnormalized contraction coefficients
C
      IF (FRMT .EQ. 'F' .OR. FRMT .EQ. 'f') THEN
         FREE = .TRUE.
      ELSE IF (FRMT .EQ. 'H' .OR. FRMT .EQ. 'h') THEN
         FREE = .FALSE.
         NCOL = 4
         FMT1 = '(4F20.8)'
         FMT2 = '(20X,3F20.8)'
         FMT3 = '(F20.8)'
      ELSE
         FREE = .FALSE.
         NCOL = 8
         FMT1 = '(8F10.4)'
         FMT2 = '(10X,7F10.4)'
         FMT3 = '(F10.4)'
      END IF
      IF(NRCIJ.NE.0) THEN
        DO 100 L = 1, NUCIJ
           CALL RDLINE(LUINFO)
           IF (FREE) THEN
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D)
              ISTART = 1
              CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,ALPHA(L),'REA')
              DO M = 1, NRCIJ
                 CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,
     &                       CPRIMU(L,M),'REA')
              END DO
#else
              READ (MLINE(NMLINE),*,IOSTAT=IOS)
     &         ALPHA(L), (CPRIMU(L,M),M=1,NRCIJ)
              IF (IOS.NE.0) THEN
                WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                WRITE(LUPRI,*) MLINE(NMLINE)
                CALL QUIT(
     &    'Error in reading ALPHA(L), (CPRIMU(L,M),M=1,NRCIJ)')
              ENDIF
#endif
           ELSE
              READ (MLINE(NMLINE),FMT1,IOSTAT=IOS) ALPHA(L),
     &           (CPRIMU(L,M), M = 1, MIN(NRCIJ,NCOL - 1))
              IF (IOS.NE.0) THEN
                WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                WRITE(LUPRI,*) MLINE(NMLINE)
                CALL QUIT(
     &         'Error in basis set reading ALPHA(L),'//
     &         '(CPRIMU(L,M), M = 1, MIN(NRCIJ,NCOL - 1))')
              ENDIF

              DO K = 2, (NRCIJ - 1)/(NCOL - 1) + 1
                 CALL RDLINE(LUINFO)
                 READ (MLINE(NMLINE),FMT2,IOSTAT=IOS) (CPRIMU(L,M),
     &              M = 1 + (NCOL-1)*(K-1), MIN(NRCIJ,(NCOL-1)*K))
                 IF (IOS.NE.0) THEN
                   WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                   WRITE(LUPRI,*) MLINE(NMLINE)
                   CALL QUIT('Error in basis set reading...') 
                 ENDIF
              END DO
           END IF
  100   CONTINUE
      ELSE
C
C     Read only exponents
C
        NRCIJ = NUCIJ
        DO 300 L = 1, NUCIJ
           CALL RDLINE(LUINFO)
           IF (FREE) THEN
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D)
              ISTART = 1
              CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,ALPHA(L),'REA')
#else
              READ (MLINE(NMLINE),*,IOSTAT=IOS) ALPHA(L)
              IF (IOS.NE.0) THEN
                WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
                WRITE(LUPRI,*) MLINE(NMLINE)
                CALL QUIT('Error in basis set reading of ALPHA(L)')
              ENDIF
#endif
           ELSE
              READ (MLINE(NMLINE),FMT3,IOSTAT=IOS) ALPHA(L)
              IF (IOS.NE.0) THEN
                WRITE(LUPRI,*) 'Error reading ALPHA(L), line',NMLINE,':'
                WRITE(LUPRI,*) MLINE(NMLINE)
                CALL QUIT('Error in basis set reading of ALPHA(L)')
              ENDIF
           END IF
        CPRIMU(L,L) = D1
  300   CONTINUE
      ENDIF
C
C     Test print
C
      IF (INPTST) THEN
         WRITE(LUPRI,'(A,I4,A,I4,A)')
     &   'Contraction matrix (',NUCIJ,'x',NRCIJ,') :' 
         DO 400 L = 1, NUCIJ
            WRITE (LUPRI,'(E16.10)') ALPHA(L)
            WRITE (LUPRI,'(6(2X,F10.8))') (CPRIMU(L,M),M=1,NRCIJ)
  400    CONTINUE
      END IF
      RETURN
      END
C  /* Deck segorb */
      SUBROUTINE SEGORB(SEGIJ,NUCIJ,NRCIJ,CPRIMU,KPRIM,DSM)
#include "implicit.h"
#include "priunit.h"
      LOGICAL SEGIJ
      DIMENSION CPRIMU(KPRIM,KPRIM)
C
      DO 100 L = 1, NUCIJ
         NONZER = 0
         DO 200 M = 1, NRCIJ
            IF (ABS(CPRIMU(L,M)).GT.DSM) NONZER = NONZER + 1
  200    CONTINUE
         SEGIJ = SEGIJ .AND. (NONZER .LE. 1)
  100 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nrmorb */
      SUBROUTINE NRMORB(I,NRCIJ,NUCIJ,ALPHA,CPRIM,CPRIMU,KPRIM,NBLOCK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Purpose: do the normalization of incoming atomic orbitals, GTO -> CGTP
C
C On entry:
C           I - of 1 - IQM
C           NRCIJ, NUCIJ - number of contracted/uncontracted AO
C           ALPHA - exponents  
C           CPRIMU - contraction coefficients read from the basis set input
C           KPRIM, NBLOCK - dimensions
C
C On output: 
C            CPRIM - HERMIT style normalized contraction coefficients
C
C Description added by Miro ILIAS, March 2007, Tel Aviv
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "cbirea.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0,
     &           DP25 = 0.25D0, DP5 = 0.5D0, DP75 = 0.75D0,
     &           THRMIN = 1.D-17)
      CHARACTER*1 SPDCAR

      DIMENSION ALPHA(KPRIM), CPRIM(KPRIM,KPRIM), CPRIMU(KPRIM,KPRIM)

      LOGICAL ZERONORM
C
      ZERONORM=.FALSE.
      PIPPI = (DP5/PI)**DP75
      DO 100 N = 1, NRCIJ
         SUM = D0
         DO 200 L = 1, NUCIJ
         DO 200 M = 1, NUCIJ
            T = D2*SQRT(ALPHA(L)*ALPHA(M))/(ALPHA(L)+ALPHA(M))
            SUM = SUM + CPRIMU(L,N)*CPRIMU(M,N)*(T**(I + DP5))
  200    CONTINUE
         IF (SQRT(SUM) .LT. THRMIN) THEN ! error branch - print output
           ZERONORM=.TRUE.
           IPREAD = 10
           GOTO 900  ! ensure the full control print-out
         ENDIF
         SUM=D1/SQRT(SUM)
         DO 300 L=1, NUCIJ
            CPRIM(L,N)=CPRIMU(L,N)*SUM*(D4*ALPHA(L))**(DP5*I+DP25)*PIPPI
  300    CONTINUE
  100 CONTINUE

C     ... control print-out
 900   CONTINUE
      IF (IPREAD.GE.10) THEN
       WRITE(LUPRI,'(/8X,A,I3/)') '***  Output from NRMORB, block #:',I
       WRITE(LUPRI,'(2X,A,I3,A,I3)')
     & 'uncontracted / contracted shells, NUCIJ/NRCIJ: ',NUCIJ,'/',NRCIJ
C      WRITE(LUPRI,'(2X,A)')
C    & 'CPRIMU - read contr.coeff. and CPRIM - normalized contr.coeff:'
       DO N = 1, NRCIJ
        WRITE(LUPRI,'(2X,A,I2)') ' contracted shell #',N
        DO L = 1,NUCIJ
         WRITE(LUPRI,'(2X,I3,A,F20.7,A,2I3,F12.7,A,F12.7)')
     &   L,'. exponent:',ALPHA(L),' contr.coeff:',L,N,CPRIMU(L,N),
     &   ' normalized contr.coeff.', CPRIM(L,N)
        ENDDO
       ENDDO
       IF (ZERONORM) GOTO 1000 ! error brach 
      ENDIF

      flush(lupri)
      RETURN
 1000 CONTINUE
        WRITE (LUPRI,'(/1X,A,2(I3,A),A1,A/A,2I4)')
     &    'NRMORB INPUT ERROR: CGTO no.',N,' for block',NBLOCK,
     &    ' of ',SPDCAR(I-1),' type has zero norm.',
     &    ' Contraction matrix: ',NUCIJ,NRCIJ
        WRITE(LUPRI,'(2X,A)')
     &  'This might indicate broken basis set input, please check it!'
        CALL QUIT('CGTO with zero norm. Check the basis set input !')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck orbpro */
      SUBROUTINE ORBPRO(IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           ALPHA,CPRIM,CPRIMU,KATOM,KANG,KBLOCK,KPRIM,
     &           ISHELL,IPRIM,IPRIMD,IORB,IORBD,LCOMP)
C******************************************************************************
C
C     Process basis data
C
C
CMI LCOMP is 1 (L) or 2 (S component)
C
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
C
#include "infpar.h"
      PARAMETER (D0 = 0.0D0)
C
#include "cbirea.h"
      LOGICAL SEG,SPHER,SEGAUX
      DIMENSION IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          SEG(KBLOCK),ALPHA(KPRIM,KBLOCK),
     &          CPRIM(KPRIM,KPRIM,KBLOCK),CPRIMU(KPRIM,KPRIM,KBLOCK)
#include "ccom.h"
#include "nuclei.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
      ICENT  = 0
      JBLOCK = 0
      DO 10 I = 1,NONTYP ! atomic type
         DO 20 N = 1,NONT(I) ! symmetry-independent center
            ICENT = ICENT + 1
            KBCH = JBLOCK
            NDEG = NUCDEG(ICENT)
            DO 30 J = 1,IQM(I) ! angular momentum value
               KKK = 0
               NCOMP = KHK(J)
               NCCMP = KCK(J)
               SPHER = SPH(J)
               DO 40 K = 1, JCO(J,I) ! blocks of given atom and L-value
                  KBCH = KBCH + 1
                  NCONT = NRC(KBCH)
                  IF(NCONT.GT.MXCONT) GOTO 5000
C                 ...  dimension PRICCF(MXPRIM,MXCONT)
                  DO 50 KK = 1, NCONT ! contracted shells
                     KKK = KKK + 1
                     ISHELL = ISHELL + 1
C*****************************************************************************
C     Data on the AO-block associated with a given shell(ISHELL):
C       NUCO    - number of uncontracted functions
C       NRCO    - number of contracted functions
C     Data on a given shell in an AO-block:
C       NCENT   - index of symmetry independent center
C       NUMCF   - index of shell in AO-block
C       NBCH    - index of block in AO-vector
C       ISTBAO  - stabiliser: basic sym. op. that do not move center
C       NHKT    - angular quantum number (s=1,p=2,d=3 etc.)
C       KHKT    - number of spherical (Cartesian) components
C       KCKT    - number of Cartesian components
C       SEGM    - segmented contraction
C       LCLASS  - component: large(1) or small(0)
C       CENT    - coordinates of center
C        Card 6: ATOMIC DATA
C         Q     - nuclear charge
C         NONT  - number of symmetry independent atoms
C         IQM   - maximum angular quantum number (s=1,p=2,d=3 etc.)
C         NHTYP - maximum angular quantum number for ALL orbitals
C         JCO   - number of AO-blocks for each l-value

C*****************************************************************************
                     NUCO  (ISHELL) = NUC(KBCH)
                     NRCO  (ISHELL) = NCONT
                     JSTRT (ISHELL) = IPRIM
                     NSTRT (ISHELL) = IORB
                     KSTRT (ISHELL) = IORBD
                     NCENT (ISHELL) = ICENT
                     NUMCF (ISHELL) = KK
                     NUMCFT(ISHELL) = KKK
                     NBCH  (ISHELL) = KBCH
                     SHARE (ISHELL) = .FALSE.
                     ISTBAO(ISHELL) = ISTBNU(ICENT)
                     NHKT  (ISHELL) = J
                     KHKT  (ISHELL) = NCOMP
                     SPHR  (ISHELL) = SPHER
                     KCKT  (ISHELL) = NCCMP
                     SEGM  (ISHELL) = SEG(KBCH)
                     LCLASS(ISHELL) = LCOMP
                     CENT(ISHELL,1,1) = CORD(1,ICENT)
                     CENT(ISHELL,2,1) = CORD(2,ICENT)
                     CENT(ISHELL,3,1) = CORD(3,ICENT)
C                    Basis-set identifier (WK/UniKa/04-11-2002).
C                    MBSID(ISHELL) = MULBSI(ICENT)
                     IORB  = IORB  + NCOMP
                     IORBD = IORBD + NCOMP*NDEG
                     
                     IF (.NOT. SEGM(ISHELL) .AND. LCOMP .EQ. 1)
     &                         SEGBAS = .FALSE.
                     IF (.NOT. SEGM(ISHELL) .AND. LCOMP .EQ. 3)
     &                         SEGAUX = .FALSE.
   50             CONTINUE
                  DO 60 L = 1, NUC(KBCH)
                     IPRIM = IPRIM + 1
C*****************************************************************************
C       PRIEXP  - exponent of primitive shell
C       PRICCF  - normalized contraction coefficient
C       PRICRX  - x-coordinate of center
C       PRICRY  - y-coordiante of center
C       PRICRZ  - z-coordinate of center
C*****************************************************************************
                     PRIEXP(IPRIM) = ALPHA(L,KBCH)
                     DO 70  M = 1, NCONT
                        PRICCF(IPRIM,M) = CPRIM(L,M,KBCH)
   70                CONTINUE
                     PRICRX(IPRIM) = CORD(1,ICENT)
                     PRICRY(IPRIM) = CORD(2,ICENT)
                     PRICRZ(IPRIM) = CORD(3,ICENT)
                     IPRIMD = IPRIMD + NCOMP*NDEG
   60             CONTINUE
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
         JBLOCK = JBLOCK + NBLCK(I)
   10 CONTINUE
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *   ' >>> ORBPRO error, '//
     *   'no. of contracted functions per AO-block',NCONT,
     *   '                   '//
     *   'current maximum number MXCONT           ',MXCONT
      CALL QUIT('Too many contracted functions per AO-block')
      END
C  /* Deck fitpro */
      SUBROUTINE FITPRO(IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           ALPHA,CPRIM,CPRIMU,KATOM,KANG,KBLOCK,KPRIM,
     &           ISHELL_FIT,IPRIM_FIT,IPRIMD_FIT,
     &           IORB_FIT,IORBD_FIT,LCOMP_FIT)
C******************************************************************************
C
C     Process basis data
C     written by Luuk Visscher, trivial modification of orbpro
C
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
C
#include "cbirea.h"
      LOGICAL SEG,SPHER
      DIMENSION IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          SEG(KBLOCK),ALPHA(KPRIM,KBLOCK),
     &          CPRIM(KPRIM,KPRIM,KBLOCK),CPRIMU(KPRIM,KPRIM,KBLOCK)
#include "ccom.h"
#include "nuclei.h"
#include "fitprimit.h"
#include "fitshells.h"
#include "symmet.h"
#include "aosotr.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
      ICENT  = 0
      JBLOCK = 0
      DO 10 I = 1,NONTYP
         DO 20 N = 1,NONT(I)
            ICENT = ICENT + 1
            KBCH = JBLOCK
            NDEG = NUCDEG(ICENT)
            DO 30 J = 1,IQM(I)
               KKK = 0
               NCOMP = KHK(J)
               NCCMP = KCK(J)
               SPHER = SPH(J)
               DO 40 K = 1, JCO(J,I)
                  KBCH = KBCH + 1
                  NCONT = NRC(KBCH)
                  IF(NCONT.GT.MXCONT) GOTO 5000
C                 ...  dimension PRICCF(MXPRIM,MXCONT)
                  DO 50 KK = 1, NCONT
                     KKK = KKK + 1
                     ISHELL_FIT = ISHELL_FIT + 1
C*****************************************************************************
C     Data on the AO-block associated with a given shell(ISHELL):
C       NUCO    - number of uncontracted functions
C       NRCO    - number of contracted functions
C     Data on a given shell in an AO-block:
C       NCENT   - index of symmetry independent center
C       NUMCF   - index of shell in AO-block
C       NBCH    - index of block in AO-vector
C       ISTBAO  - stabiliser: basic sym. op. that do not move center
C       NHKT    - angular quantum number (s=1,p=2,d=3 etc.)
C       KHKT    - number of spherical (Cartesian) components
C       KCKT    - number of Cartesian components
C       SEGM    - segmented contraction
C       LCLASS  - component: large(1) or small(0)
C       CENT    - coordinates of center
C*****************************************************************************
                     NUCO_FIT  (ISHELL_FIT) = NUC(KBCH)
                     NRCO_FIT  (ISHELL_FIT) = NCONT
                     JSTRT_FIT (ISHELL_FIT) = IPRIM_FIT
                     NSTRT_FIT (ISHELL_FIT) = IORB_FIT
                     KSTRT_FIT (ISHELL_FIT) = IORBD_FIT
                     NCENT_FIT (ISHELL_FIT) = ICENT
                     NUMCF_FIT (ISHELL_FIT) = KK
                     NUMCFT_FIT(ISHELL_FIT) = KKK
                     NBCH_FIT  (ISHELL_FIT) = KBCH
                     SHARE_FIT (ISHELL_FIT) = .FALSE.
                     ISTBAO_FIT(ISHELL_FIT) = ISTBNU(ICENT)
                     NHKT_FIT  (ISHELL_FIT) = J
                     KHKT_FIT  (ISHELL_FIT) = NCOMP
                     SPHR_FIT  (ISHELL_FIT) = SPHER
                     KCKT_FIT  (ISHELL_FIT) = NCCMP
                     SEGM_FIT  (ISHELL_FIT) = SEG(KBCH)
                     LCLASS_FIT(ISHELL_FIT) = LCOMP_FIT
                     CENT_FIT(ISHELL_FIT,1,1) = CORD(1,ICENT)
                     CENT_FIT(ISHELL_FIT,2,1) = CORD(2,ICENT)
                     CENT_FIT(ISHELL_FIT,3,1) = CORD(3,ICENT)
                     IORB_FIT  = IORB_FIT  + NCOMP
                     IORBD_FIT = IORBD_FIT + NCOMP*NDEG
   50             CONTINUE
                  DO 60 L = 1, NUC(KBCH)
                     IPRIM_FIT = IPRIM_FIT + 1
C*****************************************************************************
C       PRIEXP  - exponent of primitive shell
C       PRICCF  - normalized contraction coefficient
C       PRICRX  - x-coordinate of center
C       PRICRY  - y-coordiante of center
C       PRICRZ  - z-coordinate of center
C*****************************************************************************
                     PRIEXP_FIT(IPRIM_FIT) = ALPHA(L,KBCH)
                     DO 70  M = 1, NCONT
                        PRICCF_FIT(IPRIM_FIT,M) = CPRIM(L,M,KBCH)
   70                CONTINUE
                     PRICRX_FIT(IPRIM_FIT) = CORD(1,ICENT)
                     PRICRY_FIT(IPRIM_FIT) = CORD(2,ICENT)
                     PRICRZ_FIT(IPRIM_FIT) = CORD(3,ICENT)
                     IPRIMD_FIT = IPRIMD_FIT + NCOMP*NDEG
   60             CONTINUE
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
         JBLOCK = JBLOCK + NBLCK(I)
   10 CONTINUE
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *   ' >>> FITPRO error, '//
     *   'no. of contracted functions per AO-block',NCONT,
     *   '                   '//
     *   'current maximum number                  ',MXCONT
      CALL QUIT('Too many contracted functions per AO-block')
      END
C  /* Deck orbout */
      SUBROUTINE ORBOUT(ICLS,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           IPCON,KATOM,KANG,KBLOCK,KPRIM,ALPHA,CPRIMU,NRMPRI)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
C
#include "ccom.h"
#include "cbirea.h"
#include "nuclei.h"
#include "primit.h"
#include "fitprimit.h"
      CHARACTER    CHRSEG*10,DEG*2
      CHARACTER*60 FMTA, FMTB, FMTC, FMTD
      LOGICAL      SPHER, SEG, NRMPRI
!miro: out-of-bounds fixes: KBLOCK -> *
      DIMENSION IQM(KATOM),NBLCK(KATOM),
!    &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          JCO(KANG,KATOM),NUC(*),NRC(*),
!    &          CPRIMU(KPRIM,KPRIM,KBLOCK),
     &          CPRIMU(KPRIM,KPRIM,*),
     &          ALPHA(KPRIM,KBLOCK),
!    &          SEG(KBLOCK)
     &          SEG(*)
      DIMENSION IPCON(KPRIM)
C
C     P R I M I T I V E S
C     ===================
C
      IF (IPREAD .GT. 0) THEN
       IF (.NOT.NRMPRI) THEN
        CALL HEADER('Orbital exponents and contraction coefficients',1)
       ELSE
        CALL HEADER
     &  ('Orbital exponents and normalized contraction coefficients',1)
       END IF
      END IF
      IF(ICLS.EQ.1) THEN
        IPRIMD = 0
        IPRIM =  0
      ELSEIF(ICLS.EQ.2) THEN
        IPRIMD = NPLRG
        IPRIM  = NPLSH
      ELSEIF(ICLS.EQ.3) THEN
        IPRIMD = 0
        IPRIM  = 0
      ELSEIF(ICLS.EQ.4) THEN
        IPRIMD = NPLRG_FIT
        IPRIM  = NPLSH_FIT
      ENDIF
      ICENT  = 0
      JBLOCK = 0
      DO 100 I = 1, NONTYP
         DO 110 N = 1, NONT(I)
            ICENT = ICENT + 1
            NDEG  = NUCDEG(ICENT)
            KBCH  = JBLOCK
            DO 200 J = 1, IQM(I)
            DO 200 K = 1, JCO(J,I)
               KBCH = KBCH + 1
               NNUC  = NUC(KBCH)
               NNRC  = NRC(KBCH)
            IF (NNUC .EQ. 0) GO TO 200
               ITYP = NHKOFF(J)
               IPSTRT = IPRIM + 1
               IPRIM  =  IPRIM + NNUC
               IF (UNCONT) THEN
                  CHRSEG = 'uncontrac.'
                  NNRC = 0
               ELSE IF (SEG(KBCH)) THEN
                  CHRSEG = 'seg. cont.'
               ELSE
                  CHRSEG = 'gen. cont.'
               END IF
               IF (NNRC .LE. 5) THEN
                  FMTA='(/2X,A4,A2,1X,A4,I5,1P,D18.10,0P,2X,5F8.4)'
!                 FMTB='(/2X,A4,3X,   A4,I5,1P,D18.10,0P,2X,5F8.4)'
                  FMTC='( 3X,A10,        I5,1P,D18.10,0P,2X,5F8.4)'
                  FMTD='(13X,            I5,1P,D18.10,0P,2X,5F8.4)'
               ELSE
          FMTA='(/2X,A4,A2,1X,A4,I5,1P,D18.10,0P,2X,5F8.4/,(37X,5F8.4))'
!         FMTB='(/2X,A4,3X,   A4,I5,1P,D18.10,0P,2X,5F8.4/,(37X,5F8.4))'
          FMTC='( 3X,A10,        I5,1P,D18.10,0P,2X,5F8.4/,(37X,5F8.4))'
          FMTD='(13X,            I5,1P,D18.10,0P,2X,5F8.4/,(37X,5F8.4))'
               END IF
               DEG = '  '
               DO 400 L = 1, NDEG
                 IF (NDEG .GT. 1) WRITE(DEG,'(A1,I1)') '#',L
                 ITYP = NHKOFF(J)
                 DO 300 ICOMP = 1, KHK(J)
                     ITYP = ITYP + 1
                     IPRIMD = IPRIMD + 1
                     IF (IPREAD .GT. 0) THEN
                      IF (ICLS.LE.2) THEN
                       WRITE (LUPRI,FMTA)  NAMN(ICENT),DEG,
     &                                     GTOTYP(ITYP),
     &                                     IPRIMD,PRIEXP(IPSTRT),
     &                                     (CPRIMU(1,MM,KBCH),MM=1,NNRC)
                      ELSE
                       WRITE (LUPRI,FMTA)  NAMN(ICENT),DEG,
     &                                     GTOTYP(ITYP),
     &                                     IPRIMD,PRIEXP_FIT(IPSTRT),
     &                                     (CPRIMU(1,MM,KBCH),MM=1,NNRC)
                      ENDIF
                     ENDIF
                     IF (NNUC .GT. 1) THEN
                        IPRIMD = IPRIMD + 1
                        IF (IPREAD .GT. 0) THEN
                          IF (ICLS.LE.2) THEN
                           WRITE(LUPRI,FMTC) CHRSEG,
     &                                    IPRIMD,PRIEXP(IPSTRT+1),
     &                                    (CPRIMU(2,MM,KBCH),MM=1,NNRC)
                          ELSE
                           WRITE(LUPRI,FMTC) CHRSEG,
     &                                     IPRIMD,PRIEXP_FIT(IPSTRT+1),
     &                                     (CPRIMU(2,MM,KBCH),MM=1,NNRC)
                          ENDIF
                        END IF
                     END IF
                     DO 410 M = 3, NNUC
                        IPRIMD = IPRIMD + 1
                        IF (IPREAD .GT. 0) THEN
                          IF (ICLS.LE.2) THEN
                           WRITE (LUPRI,FMTD) IPRIMD,
     &                                      PRIEXP(IPSTRT-1+M),
     &                                     (CPRIMU(M,MM,KBCH),MM=1,NNRC)
                         ELSE
                           WRITE (LUPRI,FMTD) IPRIMD,
     &                                      PRIEXP_FIT(IPSTRT-1+M),
     &                                     (CPRIMU(M,MM,KBCH),MM=1,NNRC)
                          END IF
                       END IF
  410              CONTINUE
  300            CONTINUE
  400          CONTINUE
  200       CONTINUE
  110    CONTINUE
         JBLOCK = JBLOCK + NBLCK(I)
  100 CONTINUE
C
      IF (NRMPRI) RETURN
C
C       C O N T R A C T I O N
C       =====================
C
C       Looping is over
C                       - atomic type
C                         - symmetry independent center
C                           - shell
C                             - degeneracy of center
C                               - component
C       [NOTE that in DALTON the looping over degeneracy and components
C        is reversed !]

C
      CALL HEADER('Contracted Orbitals',1)

      IF (UNCONT) THEN
         WRITE(LUPRI,'(A)') ' All orbitals are uncontracted.'
         RETURN
      END IF
      IF(ICLS.EQ.1) THEN
        IPRIMD  = 0
        IORBD   = 0
        IPRIM   = 0
      ELSE
        IPRIMD  = NPLRG
        IORBD   = NLARGE
        IPRIM   = NPLSH
      ENDIF
      ICENT   = 0
      JBLOCK  = 0
      DO 500 I = 1,NONTYP
         DO 510 N = 1,NONT(I)
            ICENT = ICENT + 1
            NDEG  = NUCDEG(ICENT)
            KBCH  = JBLOCK
            DO 600 J = 1,IQM(I)
               NCOMP = KHK(J)
               DO 610 K = 1,JCO(J,I)
                  KBCH = KBCH + 1
                  NNUC = NUC(KBCH)
                  NNRC  = NRC(KBCH)
                  IPSTRT = IPRIM + 1
                  IPRIM  = IPRIM + NNUC
                  DO 700 L = 1,NNRC
                     JPRIM = 0
                     JPRIMD = IPRIMD
                     DO 750 M = IPSTRT,IPRIM
                        JPRIMD = JPRIMD + 1
                        IF(PRICCF(M,L).NE.0.00D0) THEN
                           JPRIM = JPRIM + 1
                           IPCON(JPRIM) = JPRIMD
                        ENDIF
  750                CONTINUE
                     IOFF = 0
                     DEG = '  '
                     DO 810 LL = 1,NDEG
                       IF(NDEG.GT.1) WRITE(DEG,'(A1,I1)') '#',LL
                       ITYP = NHKOFF(J)
                       DO 800 ICOMP = 1, NCOMP
                         ITYP = ITYP + 1
                         IORBD = IORBD + 1
                         IF (IPREAD .GT. 0) THEN
                            WRITE(LUPRI,1040)
     &                      IORBD,NAMN(ICENT),DEG,GTOTYP(ITYP),
     &                      ((IPCON(M)+IOFF),M=1,JPRIM)
                         END IF
                        IOFF = IOFF + NNUC
  800                  CONTINUE
  810                CONTINUE
  700             CONTINUE
                  IPRIMD = IPRIMD + NNUC*NDEG*NCOMP
  610          CONTINUE
  600       CONTINUE
  510    CONTINUE
         JBLOCK = JBLOCK + NBLCK(I)
  500 CONTINUE
      IF (IPREAD .GT. 0) THEN
         WRITE (LUPRI,'(/)')
      END IF
 1040 FORMAT(I5,2X,A4,A2,2X,A4,24I4)
 1050 FORMAT(I4,2X,A4,4X,A4,12I6)
C
 4242 CONTINUE
      RETURN
      END
C  /* Deck nucpro */
      SUBROUTINE NUCPRO(WORK,LWORK)
C******************************************************************************
C
C     Process molecular data
C
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
C
#include "cbirea.h"
      DIMENSION WORK(LWORK)
#include "ccom.h"
#include "nuclei.h"
#include "frame.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#include "dgroup.h"
C     
C
C
C     Statement function to "convert" 0 - 7 to a printable
C     "binary" form
C
      IFAKBN(I) = 25*IAND(I,4) + 5*IAND(I,2) + IAND(I,1)

      CALL QENTER('NUCPRO')
C
C     *************************************************************
C     ***** Calculate Distances ,Nuclear Potential Energy and *****
C     ***** Nuclear Contributions to Dipole Moment            *****
C     * LOOPS:
C     * 10      - loop over first symmetry independent centre
C     * 20+IF   - loop over dependent centers
C     * 30      - loop over second symmetry independent centre
C     * 40+IF   - loop over dependent centers
C     ************************************************************
C
      DPNUC(1) = D0
      DPNUC(2) = D0
      DPNUC(3) = D0
      NUCDEP    = 0
      NUCPRE(1) = 0
      NFLOAT    = 0

      DO 100 N = 1, NUCIND
         IVARB = ISTBNU(N)
         NDEG  = MULT(IVARB) !degeneracy of center
         IF (IPREAD .GE. 5)
     &      WRITE (LUPRI,1000) NAMN(N),(CORD(I,N),I=1,3),
     &                        IFAKBN(ISTBNU(N)), NDEG - 1
         NUCDEG(N) = NDEG
         NUCPRE(N + 1) = NUCPRE(N) + NDEG
         IF (NINT(CHARGE(N)) .EQ. 0 .AND. NAMN(N) .NE. 'cav')
     &        NFLOAT = NFLOAT + 1
         DPNUC(1) = DPNUC(1) + CORD(1,N)*CHARGE(N)*FMULT(IVARB)
         DPNUC(2) = DPNUC(2) + CORD(2,N)*CHARGE(N)*FMULT(IVARB)
         DPNUC(3) = DPNUC(3) + CORD(3,N)*CHARGE(N)*FMULT(IVARB)
         II = 0
         DO 200 LA = 0,MAXREP
           IF (IAND(IVARB,LA) .EQ. 0) THEN ! filter out generator products containing generators that stabilise the center
              II = II + 1
              NUCDEP = NUCDEP + 1
              NUCNUM(N,LA+1) = NUCDEP
              IF (NUCDEP .GT. MXCENT) GO TO 200
C             ... otherwise we may get segmentation fault from NAMDEP,
C             NAMDPX, and NDEGNM assignments. We exit after 200 loop
C             when total NUCDEP is known. /hjaaj Jan 2010
              IF (NDEG .EQ. 1 ) THEN
                 NAMDEP(NUCDEP)     = NAMEX(3*N)(1:4)//'  '
                 NAMDPX(3*NUCDEP-2) = NAMEX(3*N-2)
                 NAMDPX(3*NUCDEP-1) = NAMEX(3*N-1)
                 NAMDPX(3*NUCDEP  ) = NAMEX(3*N  )
                 NDEGNM(NUCDEP)     = 1
              ELSE
                 ICHARD = ICHAR('0') + II
                 NAMDEP(NUCDEP)     = NAMEX(3*N)(1:4)//' '
     &               //CHAR(ICHARD)
                 NAMDPX(3*NUCDEP-2) = NAMEX(3*N-2)(1:3)
     &               //CHAR(ICHARD)//NAMEX(3*N-2)(5:6)
                 NAMDPX(3*NUCDEP-1) = NAMEX(3*N-1)(1:3)
     &               //CHAR(ICHARD)//NAMEX(3*N-1)(5:6)
                 NAMDPX(3*NUCDEP  ) = NAMEX(3*N  )(1:3)
     &               //CHAR(ICHARD)//NAMEX(3*N  )(5:6)
                 NDEGNM(NUCDEP)     = II
              END IF
           ELSE
              NUCNUM(N,LA+1) = 0
           END IF
  200    CONTINUE
  100 CONTINUE
      NATOMS = NUCDEP - NFLOAT
      IF (LINEAR .EQV. .TRUE.) THEN
        IF((NATOMS .Eq. 1).OR.((NATOMS.EQ.0).AND.(NFLOAT.EQ.1))) 
     &    ATOMIC = .TRUE.
      ENDIF
      IF (NUCDEP .GT. MXCENT) GOTO 5010
C     
C     Nuclear repulsion energy
C     ========================
C
      POTNUC =  GETPOT(CHARGE,IPREAD)
C
C     **********************************************************
C     ***** Calculate symmetry-adapted nuclear coordinates *****
C     **********************************************************
C
      CALL SYMNCO(LWORK,WORK)
C
      CALL QEXIT('NUCPRO')
      RETURN
C
C     Error messages:
C
 5010 CONTINUE
      WRITE (LUPRI,'(///A/2(A,I5/)/A,I5,A)')
     &    ' Allowed number of atoms exceeded.',
     &    ' Number of atoms in input:',NUCDEP,
     &    ' Number of atoms allowed: ',MXCENT,
     &    ' Increase MXCENT to at least',NUCDEP,' and recompile DIrac.'
      CALL QUIT('NUCPRO: Too many atoms in READIN.')
C
 1000 FORMAT(2X,A4,5X,3F10.5,/,
     &       16X, 'Stabilizer ',I3,',   with ',I2,
     &       ' symmetry equivalent atoms',/)
      END
C  /* Deck sympro */
      SUBROUTINE SYMPRO(LVALUE,MVALUE,NVALUE,IRREP,DOOWN,blk)
C*****************************************************************************
C
C     SYMPRO generates the symmetry orbitals (SOs) from the list of
C     symmetry independent atomic orbitals(AOs).
C     A given SO of irrep L is given as:
C
C         SO(L,a) = SUM{i} [CTRAN(NSORB,i)*G(i)*AO(NSORB)]
C
C     Here {G(i)} is a set of symmetry operators transforming AO(NSORB)
C     into all corresponding orbitals centered on symmetry related sites
C     in the molecule. The loop over i is only allowed to encounter
C     elements that do not contain any basic operations belonging to the
C     stabiliser of the centre; thus any G uniquely defines transformation
C     to a given centre.
C
C     By allowing all symmetry operations to work on the list of
C     symmetry independent AOs a list of symmetry dependent AOs
C     is generated (counted by NAORB). Redundancies are removed by
C     limiting symmetry operations to the set G(i) and the list thus
C     reduced to the list of non-trivial symmetry dependent AOs
C     (counted by NSORB).
C
C******************************************************************************
C
C     This subroutine determines the following parameters:
C
C     ISYMAO  -  indicates how an atomic orbital behaves under the basic
C                symmetry operations. When a basic operation changes the
C                sign of the AO (centered at origo) the corresponding bit
C                gets the value 1.
C
C     CTRAN   -  transformation coefficients (+1 or -1)
C     IPIND   -  packed integer
C     IPTSYM  -  pointer from redundant list of symmetry dependent AOs
C                to list of non-trivial symmetry dependent AOs. The
C                pointer is zero if the AO does not contribute to the
C                irrep.
C     ITRAN   -  given an operation G(i): pointer from list of symmetry
C                dependent AOs to list of symmetry independent AOs
C     JPRX    -  AO label
C     JTRAN   -  index of operation G moving a given AO to another centre
C     MAMN    -  name of AO centre
C     NPARNU  -  offset pointer from non-symmetric operators for given
C                irrep
C     NPARSU  -  offset pointer for symmetry dependent AOs for given irrep
C     NSORB   -  number of symmetry dependent AOs
C     NSYM    -  number of irreps that the AOs contribute to
C*************************************************************************
C
C     The code originally resided inside READIN
C     tuh 120988
C
C     Polished and annotated - tsaue March 10 1993
C*************************************************************************
      USE RECP_NTR
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.0D0)
C
      LOGICAL DOOWN
      DIMENSION LVALUE(MXAQN), MVALUE(MXAQN), NVALUE(MXAQN),
     &          IRREP(MXCORB)
      integer soef, blk
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "pincom.h"
#include "aosotr.h"
#include "sphtrm.h"
#include "cbirea.h"
#include "huckel.h"
#include "argoscom.h"
#include "argoscomch.h"
C
C*****************************************************************************
C Determine ISYMAO:
C       ISYMAO is a bitstring indicating how atomic orbitals centered in
C       origo behaves under the basic symmetry operations
C*****************************************************************************
C
C     I. Cartesian functions
C     ----------------------
C
      IF (DOCART) THEN
         DO 100 I = 1, NHTYP
            CALL LMNVAL(I,I*(I + 1)/2,LVALUE,MVALUE,NVALUE)
            DO 100 J = 1, I*(I + 1)/2
                LVAR = MOD(LVALUE(J),2)*ISYMAX(1,1)
                MVAR = MOD(MVALUE(J),2)*ISYMAX(2,1)
                NVAR = MOD(NVALUE(J),2)*ISYMAX(3,1)
                ISYMAO(I,J) = IEOR(LVAR,IEOR(MVAR,NVAR))
  100    CONTINUE
C
C     II. Your own scheme
C     -------------------
C
      ELSE IF (DOOWN) THEN
         IJK = 0
         DO 200 I = 1, NHTYP
            DO 210 K = 1, KHK(I)
               ISYMAO(I,K) = -1
  210       CONTINUE
            CALL LMNVAL(I,I*(I + 1)/2,LVALUE,MVALUE,NVALUE)
            DO 220 J = 1, I*(I + 1)/2
                LVAR = MOD(LVALUE(J),2)*ISYMAX(1,1)
                MVAR = MOD(MVALUE(J),2)*ISYMAX(2,1)
                NVAR = MOD(NVALUE(J),2)*ISYMAX(3,1)
                DO 220 K = 1, KHK(I)
                   IJK = IJK + 1
                   IF (ABS(CSP(IJK)).GT.D0) THEN
                      IF (ISYMAO(I,K) .EQ. -1) THEN
                         ISYMAO(I,K) = IEOR(LVAR,IEOR(MVAR,NVAR))
                      ELSE
                         WRITE (LUPRI,'(/A/A/A,2I5)')
     &                   ' Incorrect Cartesian transformation input:',
     &                   '  Components belonging to different irreps '
     &                   //'have been combined.',
     &                   '  Angular quantum number and component:',I-1,K
                         CALL QUIT('Error in Cartesian->own bf input')
                      END IF
                   END IF
  220       CONTINUE
  200    CONTINUE
      ELSE
C
C     III. Spherical harmonics
C     ------------------------
C
         DO L = 0,NHTYP-1
           LL = L+1
           NLM = 2*L + 1
           DO ILM = 1,NLM
             M = MDEF(L,ILM)
             ISYMAO(LL,ILM) = IREPLM(L,M)
           ENDDO
         ENDDO
      END IF
C*****************************************************************************
C Determine transformation from AOs to SOs :
C   * DO-loops:
C       400     - run over irreducible representations
C       410     - run over shells
C       420     - run over components, thus (60+70):run over orbitals
C       430     - run over the symmetry operations that have no part in the
C               stabilizer of the center, that is run over a unique set
C               of symmetry operations transforming between symmetry
C               dependent centers
C   * Counting variables:
C       NSORB   - SO orbital number
C       NAORBD  - AO orbital number
C       NAORB   - AO orbital shell number
C       NAA     -
C       JKB     -
C       IRREPN  -
C*****************************************************************************
      IHNMAO    = 0
      NPARSU(1) = 0
      NSYM      = 0
      NPARLA    = 0
      CALL IZERO(NAOS,8)
      CALL IZERO(NCOS,16)
      CALL IZERO(IRREP,MXCORB)
C
C     Run over Irreducible Representations LAMBDA
C
      NSORB     = 0 ! counter for total SO number over all symmetries
      DO 400 LA = 0, MAXREP
         NAORBD = 0
         NAORB  = 0
         JBLOCK = 0
         DO 410 IA = 1, KMAX
            MULA   = ISTBAO(IA)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            NUMCFA = NUMCF(IA)
            IC     = LCLASS(IA)
            IF (NUMCFA .EQ. 1) JBLOCK = JBLOCK + 1
            IOFF = NAORBD
            DO 420 NA = 1, KHKTA
               NAORB = NAORB+1
               IVARB = IEOR(LA,ISYMAO(NHKTA,NA))
C
C              If orbital contributes to this representation:
C
c                (IVARB is a bitstring of basic operations with
C                1 in the positions where the basic operations
C                has a different parity for the irrep and the
C                AO(when centered in origo). These basic
C                operations with different parities must then
C                not be part of the stabilizer of the center
C                for the AO).
C
               IF (IAND(MULA,IVARB) .EQ. 0) THEN
                  NSORB = NSORB + 1
                  JKB = 0
C
C                 loop over symmetry dependent centers
C
                  DO 430 KB = 0, MAXREP
                  IF (IAND(KB,MULA) .EQ. 0) THEN
                     JKB              = JKB + 1
                     NAORBD           = NAORBD + 1
                     CTRAN(NSORB,JKB) = PT(IAND(KB,IVARB))
                     ITRAN(NSORB,JKB) = IOFF + (JKB-1)*KHKTA + NA
                  END IF
  430             CONTINUE
                  NAOS(LA+1)    = NAOS(LA+1) + 1
                  NCOS(LA+1,IC) = NCOS(LA+1,IC) + 1
                  IPTYP(NSORB)  = NHKOFF(NHKTA) + NA
                  IPCEN(NSORB)  = NCENT(IA)
                  ICLASS(NSORB) = IC
#ifndef PRG_DIRAC
                  IF (IC .EQ. 2 .AND. DOHUCKEL) THEN
                     IHNMAO = IHNMAO + 1
                     IHUCPT(IA) = IHNMAO
                     HUCCHR = CHARGE(NCENT(IA))
                     NHUCCO = NUMCF(IA)
                     CALL HUCFUN(HUCCHR,NHKTA,NHUCCO,HUCEXC(IHNMAO))
                  END IF
#endif
                  IRREPN        = IRREP(NAORB)
                  IRREP(NAORB)  = IRREPN + 1
C                 IPIND(NSORB)  = IA*2**8 + NA*2**4 + IRREPN
                  IPIND(NSORB)  = IA*2**16 + NA*2**8 + IRREPN
                  JTRAN(NSORB)  = JKB
                  INDGEN(NSORB) = IRREPN*2**24 + NA*2**16 + NUMCFA*2**8
     &                                                    + JBLOCK
                  IPTSYM(NAORB,LA) = NSORB
               ELSE
                  IPTSYM(NAORB,LA) = 0
                  NAORBD = NAORBD + MULT(MULA)
               END IF
  420       CONTINUE
  410    CONTINUE
         IF (LA .GT. 0) NPARSU(LA+1) = NPARSU(LA) + NPARLA
         NBASI  = NAOS(LA+1)
         NPARLA = NBASI*(NBASI + 1)/2
         IF (NBASI .GT. 0) NSYM = NSYM + 1
  400 CONTINUE
C
C     Compute symmetry offsets for non-symmetric operators
C     (Note case IREPO = 0 is thereby excluded)
C
      DO 500 IREPO = 1,MAXREP
         IBLK = 0
         DO 500 IREPA = 0,MAXREP
            IREPB = IEOR(IREPO,IREPA)
            IF (IREPA .GT. IREPB) THEN
               NPARNU(IREPO+1,IREPA+1) = IBLK
               IBLK = IBLK + NAOS(IREPA+1)*NAOS(IREPB+1)
            ENDIF
 500  CONTINUE
C
C     Compute symmetry offsets for basis functions
C
      IBAS = 0
      DO ISYM = 1,MAXREP+1
        DO IC = 1,2
          ICOS(ISYM,IC) = IBAS
          IBAS = IBAS + NCOS(ISYM,IC)
        ENDDO
      ENDDO
      DO JC = 1,2
        DO IC = 1,2
         IF (IPREAD.GT.5)write(lupri,*) 'I2COSX: ic1,ic2', ic,jc
          DO JSYM = 1,MAXREP+1
            DO ISYM = 1,MAXREP+1
              I2COSX(ISYM,JSYM,IC,JC) =
     &               NSORB*ICOS(JSYM,JC)+ICOS(ISYM,IC)
            ENDDO
            if (ipread.gt.5)
     &         write(lupri,*) (i2cosx(isym,jsym,ic,jc),isym=1,maxrep+1)
          ENDDO
        ENDDO
      ENDDO

C
C     For ECP integral
C
      CALL RECP_LNK_CTRAN(blk)
C
      RETURN
      END
C  /* Deck symout */
      SUBROUTINE SYMOUT
C*****************************************************************************
C
C     Print routine for symmetry
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "pgroup.h"
#include "pincom.h"
#include "cbirea.h"
#include "aosotr.h"
#include "chrsgn.h"
      CHARACTER COMP(2)*1
C
C*****************************************************************************
C       S Y M M E T R Y      O R B I T A L S
C*****************************************************************************
      IF (MAXREP.GT.0.AND.IPREAD.GT.0) THEN
         CALL HEADER('Symmetry Orbitals',1)
         WRITE (LUPRI,'(A,6X,8I5)')
     &     '  Number of orbitals in each symmetry: ',
     &     (NAOS(I),I=1,MAXREP+1)
      END IF
      IF(DIRAC) THEN
        COMP(1) = 'L'
        COMP(2) = 'S'
        IF (IPREAD .GT. 0) THEN
           WRITE (LUPRI,'(A,8I5)')
     &        '  Number of large orbitals in each symmetry: ',
     &        (NCOS(I,1),I=1,MAXREP+1)
           WRITE (LUPRI,'(A,8I5)')
     &        '  Number of small orbitals in each symmetry: ',
     &        (NCOS(I,2),I=1,MAXREP+1)
        END IF
      ELSE
        COMP(1) = ' '
        COMP(2) = ' '
      ENDIF
      I=0
      DO 10 LA=1,MAXREP+1
        NBI=NAOS(LA)
        IF (NBI.GT.0) THEN
           IF (MAXREP.GT.0.AND.IPREAD.GT.0)
     &          WRITE (LUPRI,'(//A,2X,A3,A1,I2,A1/)')
     &          '  Symmetry',REP(LA-1),'(',LA,')'
          DO 20 L=1,NBI
            I=I+1
            ICENT = IPCEN(I)
            J=NUCDEG(ICENT)
            IF (MAXREP.GT.0.AND.IPREAD.GT.0)
     &         WRITE (LUPRI,1060) I,NAMN(ICENT),COMP(ICLASS(I)),
     +            GTOTYP(IPTYP(I)), ITRAN(I,1),
     *            (CHRSGN(NINT(CTRAN(I,K))),ITRAN(I,K),K=2,J)
   20     CONTINUE
        ELSE
           IF (MAXREP.GT.0.AND.IPREAD.GT.0)
     &        WRITE (LUPRI,'(//2X,A,2X,A3,A1,I2,A1)')
     &        'No orbitals in symmetry',REP(LA-1),'(',LA,')'
        END IF
  10  CONTINUE
      IF (IPREAD .GT. 10) THEN
         WRITE(LUPRI,'(10X,A,/)') 'Symmetry pointer indices'
         WRITE(LUPRI,'(4X,8I4)')
     &      ((IPTSYM(I,J),J = 0,MAXREP),I=1,NORBS)
         IF (NFLOAT .GT. 0) WRITE (LUPRI,'(/10X,A,I4)')
     *       'Number of floating orbitals :',NFLOAT
      END IF
      WRITE (LUPRI,'(/1X,A,3(2X,A3,A1,I1,A1))')
     &        ' Symmetries of electric field:',
     &        (REP(ISYMAX(I,1)),'(',ISYMAX(I,1)+1,')',I=1,3)
      WRITE (LUPRI,'(/1X,A,3(2X,A3,A1,I1,A1))')
     &        ' Symmetries of magnetic field:',
     &        (REP(ISYMAX(I,2)),'(',ISYMAX(I,2)+1,')',I=1,3)
 1060 FORMAT(I5,3X,A6,1X,A1,3X,A4,5X,I5,7(2X,A,1X,I5))
      RETURN
      END
C  /* Deck symnco */
      SUBROUTINE SYMNCO(LWORK,WORK)
C     **********************************************************
C     ***** Calculate symmetry-adapted nuclear coordinates *****
C     **********************************************************
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "maxorb.h"
C
#include "chrxyz.h"
#include "cbirea.h"
#include "nuclei.h"
#include "symmet.h"
      DIMENSION WORK(LWORK)

      LENGTH = 48*NUCIND
      CALL IZERO(IPTCNT,LENGTH)
      IOFFT = 0
      IOFFR = 0
      IOFFTN = 0
      DO 100 IREP = 0,MAXREP ! loop over irreps
         ICENTA = 0
         ICENTB = 0
         DO 110 ICENT = 1,NUCIND ! loop over symmetry-independent centers
            MULC = ISTBNU(ICENT)
            DO 200 IDIRN = 1,3
               ICENTA = ICENTA + 1
               IF (IAND(MULC,IEOR(IREP,ISYMAX(IDIRN,1))).EQ.0)THEN
                  IOFFT = IOFFT + 1
                  IPTCNT(ICENTA,IREP,1) = IOFFT
               ENDIF
               IF (IAND(MULC,IEOR(IREP,ISYMAX(IDIRN,2))).EQ.0)THEN
                  IOFFR = IOFFR + 1
                  IPTCNT(ICENTA,IREP,2) = IOFFR
               ENDIF
  200       CONTINUE
            DO 300 IDIR1 = 1, 3
            DO 300 IDIR2 = IDIR1, 3
               ICENTB = ICENTB + 1
               ISYMIJ = IEOR(ISYMAX(IDIR1,1),ISYMAX(IDIR2,1))
               IF (IAND(MULC,IEOR(IREP,ISYMIJ)) .EQ. 0) THEN
                  IOFFTN = IOFFTN + 1
               ENDIF
  300       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IOFF = 0
      DO 400 IREP = 0, MAXREP
         DO 410 ICENT = 1, NUCIND
            IPTNUC(ICENT,IREP) = 0
  410    CONTINUE
         DO 420 ICENT = 1, NUCIND
            MULC = ISTBNU(ICENT)
            IF (IAND(MULC,IREP) .EQ. 0) THEN
               IOFF = IOFF + 1
               IPTNUC(ICENT,IREP) = IOFF
            END IF
  420    CONTINUE
  400 CONTINUE
C
      IF (IPREAD .GT. 10) THEN
         CALL HEADER('Symmetry-adapted nuclear coordinates',-1)
         WRITE(LUPRI,'(11X,8I4)') (IREP, IREP = 0,MAXREP)
         ICENTA = 0
         DO 500 ICENT = 1,NUCIND
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'X', (IPTCNT(ICENTA+1,IREP,1),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'Y', (IPTCNT(ICENTA+2,IREP,1),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'Z', (IPTCNT(ICENTA+3,IREP,1),IREP = 0,MAXREP)
            ICENTA = ICENTA + 3
  500    CONTINUE
         CALL HEADER('Symmetry-adapted nuclear magnetic moments',-1)
         WRITE(LUPRI,'(11X,8I4)') (IREP, IREP = 0,MAXREP)
         ICENTA = 0
         DO 510 ICENT = 1,NUCIND
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'X', (IPTCNT(ICENTA+1,IREP,2),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'Y', (IPTCNT(ICENTA+2,IREP,2),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'Z', (IPTCNT(ICENTA+3,IREP,2),IREP = 0,MAXREP)
            ICENTA = ICENTA + 3
  510    CONTINUE
         CALL HEADER('Symmetry-adapted nuclei',-1)
         WRITE (LUPRI,'(11X,8I4)') (IREP, IREP = 0, MAXREP)
         DO 520 ICENT = 1, NUCIND
            WRITE(LUPRI,'(4X,A6,2X,A1,8I6)')
     &      NAMN(ICENT),CHRXYZ(1),(IPTNUC(ICENT,IREP),IREP = 0,MAXREP)
  520    CONTINUE
      END IF
C
C     *****************************************************************
C     **** Transformation matrix Cartesian to symmetry coordinates ****
C     *****************************************************************
C
      NCOOR  = 3*NUCDEP
      KCSTRA = 1
      KSCTRA = KCSTRA + NCOOR*NCOOR
      KLAST  = KSCTRA + NCOOR*NCOOR
      IF (KLAST.GT.LWORK) CALL STOPIT('NUCPRO','TRACOR',KLAST,LWORK)

C     Initialize NCRREP(irep,1)
      CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),1,NCOOR,IPREAD)

C     Initialize NCRREP(irep,2)
      CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),2,NCOOR,IPREAD)

C     Initialize IPTCOR
      CALL TRACR(NCOOR)

      END
C  /* Deck pricar */
      SUBROUTINE PRICAR(WORK,LWORK)
C
C     This subroutine prints information about atomic coordinates
C
C     tuh 081188 - Bush elected
Clf       031104 - Bush Jr. re-elected
C
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "chrxyz.h"
#include "chrsgn.h"
C
      CHARACTER*4 NAME
      INTEGER NBASE(8), NSIGN(8)
      DIMENSION WORK(LWORK)
C
#include "nuclei.h"
#include "frame.h"
#include "symmet.h"
#include "pgroup.h"
#include "cbirea.h"
C
#include "chrnos.h"
C
#include "memint.h"
C
      NCOOR = 3*NUCDEP
      NWARN = 0
C
C     *********************************
C     ***** Cartesian Coordinates *****
C     *********************************
C
      IF (IPREAD .GT. 0) THEN
         CALL HEADER('Cartesian Coordinates (bohr)',1)
         WRITE (LUPRI,'(A,I3//)')'  Total number of coordinates:',NCOOR
         ICRX = 1
         ICRY = 2
         ICRZ = 3
         DO 100 ICENT = 1, NUCIND
            MULCNT = ISTBNU(ICENT)
            NAME   = NAMEX(3*ICENT)(1:4)
            IF (MULT(MULCNT) .EQ. 1) THEN
               WRITE (LUPRI,'(I4,3X,A,5X,A,3X,F15.10)')
     *               ICRX, NAME,CHRXYZ(-1), CORD(1,ICENT)
               WRITE (LUPRI,'(I4,12X,A,3X,F15.10)')
     *               ICRY, CHRXYZ(-2), CORD(2,ICENT)
               WRITE (LUPRI,'(I4,12X,A,3X,F15.10,/)')
     *               ICRZ, CHRXYZ(-3), CORD(3,ICENT)
               ICRX = ICRX + 3
               ICRY = ICRY + 3
               ICRZ = ICRZ + 3
            ELSE
               JATOM = 0
               DO 200 ISYMOP = 0, MAXOPR
                  IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
                     JATOM = JATOM + 1
                     WRITE (LUPRI,'(I4,3X,A,I2,3X,A,3X,F15.10)')
     *                     ICRX, NAME,JATOM,CHRXYZ(-1),
     *                     PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
                     WRITE (LUPRI,'(I4,12X,A,3X,F15.10)')
     *                     ICRY,CHRXYZ(-2),
     *                     PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
                     WRITE (LUPRI,'(I4,12X,A,3X,F15.10,/)')
     *                     ICRZ,CHRXYZ(-3),
     *                     PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
                     ICRX = ICRX + 3
                     ICRY = ICRY + 3
                     ICRZ = ICRZ + 3
                  END IF
  200          CONTINUE
            END IF
  100    CONTINUE
      END IF

!     print also coordinates in xyz format for easy copy-paste-ability
      if (ipread > 0) then
         call header('Cartesian coordinates in XYZ format (Angstrom)',
     &   1)
         write(lupri, '(i5)') nucdep
         write(lupri, *)
         do icent = 1, nucind
            mulcnt = istbnu(icent)
            name   = namex(3*icent)(1:4)
            if (mult(mulcnt) == 1) then
               write(lupri, '(a, 3f15.10)')
     &               name, 
     *               cord(1, icent)*xtang,
     *               cord(2, icent)*xtang,
     *               cord(3, icent)*xtang
            else
               do isymop = 0, maxopr
                  if (iand(isymop, mulcnt) == 0) then
                     write(lupri, '(a, 3f15.10)')
     &                  name, 
     *            pt(iand(isymax(1, 1), isymop))*cord(1, icent)*xtang,
     *            pt(iand(isymax(2, 1), isymop))*cord(2, icent)*xtang,
     *            pt(iand(isymax(3, 1), isymop))*cord(3, icent)*xtang
                  end if
               end do
            end if
         end do
      end if

C
C     ********************************
C     ***** Symmetry Coordinates *****
C     ********************************
C
      IF (LWORK .LT. 2*NCOOR*NCOOR) CALL STOPIT('PRICAR','TRACOR',
     &                                          2*NCOOR*NCOOR,LWORK)
      CALL MEMGET('REAL',KCS ,NCOOR*NCOOR,   WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSC ,NCOOR*NCOOR,   WORK,KFREE,LFREE)
      CALL TRACOR(WORK(KCS),WORK(KSC),1,NCOOR,IPREAD)
      IF (MAXREP .GT. 0) THEN
         IF (IPREAD .GT. 0) THEN
            CALL HEADER('Symmetry Coordinates',1)
            WRITE (LUPRI,'(A,8I5)')
     &            '  Number of coordinates in each symmetry: ',
     &            (NCRREP(I,1),I=0,MAXREP)
         END IF
         DO 300 ISYM = 0, MAXREP
         IF (NCRREP(ISYM,1) .GT. 0) THEN
            IF (IPREAD .GT. 0) THEN
                WRITE (LUPRI,'(//A,2X,A3,A1,I2,A1/)')
     &          '  Symmetry',REP(ISYM),'(',ISYM+1,')'
            END IF
            DO 400 IATOM = 1, NUCIND
               DO 500 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,ISYM,1)
                  IF (ISCOOR .GT. 0) THEN
                     NB = 0
                     DO 600 I = 1, NCOOR
                        IADR = KCS + (I - 1)*NCOOR + ISCOOR - 1
                        NSGN = NINT(WORK(IADR))
                        IF (NSGN .NE. 0) THEN
                           NB = NB + 1
                           NBASE(NB) = I
                           NSIGN(NB) = NSGN
                        END IF
  600                CONTINUE
                     IF (IPREAD .GT. 0) THEN
                     IF (NB .EQ. 1) THEN
                        WRITE (LUPRI,'(I5,3X,A,2X,A,I5)') ISCOOR,
     &                     NAMEX(3*IATOM)(1:4), CHRXYZ(-ICOOR), NBASE(1)
                     ELSE
                        WRITE (LUPRI,'(I5,3X,A,2X,A,3X,A,I3,'
     &                   //CHRNOS(NB-1)//'(2X,A,I5),A,I1)')
     &                   ISCOOR,NAMEX(3*IATOM)(1:4),CHRXYZ(-ICOOR),' [',
     &                   NBASE(1),(CHRSGN(NSIGN(I)),NBASE(I),I=2,NB),
     &                   ' ]/',NB
                     END IF
                     END IF
                  END IF
  500          CONTINUE
  400       CONTINUE
         END IF
  300    CONTINUE
      END IF
C
      CALL MEMREL('PRICAR',WORK,KWORK,KWORK,KFREE,LFREE)
      IF (NUCDEP .EQ. 1) THEN
         IF (IPREAD .GT. 0)
     &      WRITE (LUPRI,'(//A)') '  This is an atomic calculation.'
      ELSE
C
C     Print internuclear distances (if not an atom)
C
         IF (IPREAD .GT. 0) THEN
            CALL GEOANA(CORD,.TRUE.,.FALSE.,NBONDS,.FALSE.,WORK,LWORK)
            IF (NUCDEP .GT. 2 .AND. NBONDS .LE. 0) THEN
               NWARN = NWARN + 1
               WRITE (LUPRI,'(//A/A/)') '  WARNING:  No bonds - '//
     &       ' no atom pairs are within normal bonding distances'
     &       ,'  WARNING:  maybe coordinates were in Bohr, '//
     &        'but DIRAC were told they were in Angstrom ?'
            END IF
         END IF
C
         IF (IPREAD .GT. 0) WRITE(LUPRI,'(/A,F18.12,A)')
     &   '   Nuclear repulsion energy :',POTNUC,' Hartree'
      END IF
      RETURN
      END
C  /* Deck priord */
      SUBROUTINE PRIORD(ALPHA,CPRIMU,NPRI,NRCI,SEG,KPRIM,DSM)
C
C     Order primitive basis functions
C
#include "implicit.h"
#include "maxorb.h"
      LOGICAL SEG
      DIMENSION ALPHA(KPRIM), CPRIMU(KPRIM,KPRIM)
C
C     *****************************
C     *** Segmented contraction ***
C     *****************************
C
      IDONE = 0
      IF (SEG) THEN
         DO 100 ICONTR = 1, NRCI
            ISTART = IDONE + 1
            NLEFT  = NPRI - IDONE
C
C           Find first primitive
C           ====================
C
            IMXA = ISTART + IDAMAX(NLEFT,ALPHA(ISTART),1) - 1
            CALL DSWAP(1,ALPHA(ISTART),1,ALPHA(IMXA),1)
            CALL DSWAP(NRCI,CPRIMU(ISTART,1),KPRIM,CPRIMU(IMXA,1),
     &                 KPRIM)
C
C           Find corresponding contracted function
C           ======================================
C
            IMXC = IDAMAX(NRCI,CPRIMU(ISTART,1),KPRIM)
            CALL DSWAP(NPRI,CPRIMU(1,ICONTR),1,CPRIMU(1,IMXC),1)
C
C           Collect other primitives contributing to same contracted
C           ========================================================
C
            IPRI = 1
            DO 200 I = ISTART + 1, NPRI
               IF (ABS(CPRIMU(I,ICONTR)) .GT. DSM) THEN
                  CALL DSWAP(1,ALPHA(I),1,ALPHA(ISTART+IPRI),1)
                  CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,
     &                            CPRIMU(ISTART+IPRI,1),KPRIM)
                  IPRI = IPRI + 1
               END IF
  200       CONTINUE
C
C           Sort primitives
C           ===============
C
            IF (IPRI .GT. 2) THEN
               DO 300 I = ISTART + 1, ISTART + IPRI - 2
                  DO 400 J = I + 1, ISTART + IPRI - 1
                     IF (ALPHA(J) .GT. ALPHA(I)) THEN
                        CALL DSWAP(1,ALPHA(I),1,ALPHA(J),1)
                        CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,
     &                                  CPRIMU(J,1),KPRIM)
                     END IF
  400             CONTINUE
  300          CONTINUE
            END IF
C
            IDONE = IDONE + IPRI
  100    CONTINUE
C
C     ***************************
C     *** General contraction ***
C     ***************************
C
      ELSE
         DO 500 I = 1, NPRI - 1
            DO 600 J = I + 1, NPRI
            IF(ALPHA(J) .GT. ALPHA(I)) THEN
               CALL DSWAP(1,ALPHA(I),1,ALPHA(J),1)
               CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,CPRIMU(J,1),KPRIM)
           END IF
  600      CONTINUE
  500   CONTINUE
      END IF

!#define DEBUG_PRIMITIVES
#ifdef DEBUG_PRIMITIVES
      write(*, *) 'debug: primitives'
      do i = 1, npri
        write(*, '(e16.9)') alpha(i)
      end do
#endif /* ifdef DEBUG_PRIMITIVES */

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nucsiz */
      SUBROUTINE NUCSIZ(INUC,DELTA)
C.......................................................................
C
C     For a nucleus of atomic no. "INUC", calculate the exponent DELTA
C     for a Gaussian charge distribution approximation to the finite
C     Nucleus. Uses the Formula and Nuclear Masses from
C     REHE Newsletter No. 13 (14. June 1995), L. Visscher and K. Dyall.
C     For Charge.GT.109 Nuclear Mass = CHARGE*2.556.
C     Conversion Factor Bohr/fm 52917.7249.
C                                      T. Saue and J.K. Laerdahl 12.02.1996
C.......................................................................
#include "implicit.h"
#include "priunit.h"
#include "nucdata.h"
#include "codata.h"
      PARAMETER (DP3=1.0D0/3.0D0)
      PARAMETER (FAC=1.5D10,BF=0.836D0,CF=0.570D0)
      PARAMETER (FAC2=2.556D0)
C
      IF (INUC.LE.109.AND.INUC.GT.0) THEN
         A = NUCMASS(INUC)
      ELSE
         A = INUC*FAC2
      END IF
      A = A**DP3
      DELTA = (XTANG/(BF*A+CF))
      DELTA = FAC*DELTA*DELTA
      RETURN
      END
C  /* Deck sphlab */
      SUBROUTINE SPHLAB(IORDER,LABINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      CHARACTER LABINT(MXQN**2)*4
      CHARACTER SPDCAR*1, SIGNJ*1
#include "chrnos.h"
C
      IOFF = 0
      DO 100 I = 0, IORDER
         IF (I .EQ. 0) THEN
            LABINT(1) = '1s  '
            IOFF = IOFF + 1
         ELSE IF (I .EQ. 1) THEN
            LABINT(2) = '2px '
            LABINT(3) = '2py '
            LABINT(4) = '2pz '
            IOFF = IOFF + 3
         ELSE
            DO 200 J = -I, I
               IADR = IOFF + J + I + 1
               IF (J .LT. 0) THEN
                  SIGNJ = '-'
               ELSE IF (J .EQ. 0) THEN
                  SIGNJ = ' '
               ELSE
                  SIGNJ = '+'
               END IF
               Jabs = ABS(J)
               IF ((I+1) .LE. 9) THEN ! boundaries fix
                  LABINT(IADR) = CHRNOS(I+1)//SPDCAR(I)//
     &                        CHRNOS(Jabs)//SIGNJ
               ELSE IF (Jabs .GE. 10) THEN
                  Jabs = Jabs - 10
                  LABINT(IADR) = SPDCAR(I)//'1'//
     &                        CHRNOS(Jabs)//SIGNJ
               ELSE
                  LABINT(IADR) = SPDCAR(I)//'0'//
     &                        CHRNOS(Jabs)//SIGNJ
               ENDIF
  200       CONTINUE
            IOFF = IOFF + 2*I + 1
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck carlab */
      SUBROUTINE CARLAB(IORDER,LABINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION IX(MXAQN), IY(MXAQN), IZ(MXAQN)
      CHARACTER LABINT(MXQN*(MXQN+1)*(MXQN+2)/6)*4
      CHARACTER SPDCAR*1
#include "chrnos.h"
      J = 0
      DO 100 I = 0, IORDER
         NFUN = (I + 1)*(I + 2)/2
         IF (I .LE. 3) THEN
            J = J + NFUN
            CALL SETCLB(I,LABINT,MXQN*(MXQN+1)*(MXQN+2)/6)
         ELSE
            CALL LMNVAL(I+1,NFUN,IX,IY,IZ)
            DO 200 K = 1, NFUN
               J  = J + 1
               NX = IX(K)
               NY = IY(K)
               NZ = IZ(K)
               LABINT(J) = SPDCAR(I)//CHRNOS(NX)//CHRNOS(NY)//CHRNOS(NZ)
  200       CONTINUE
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck aotoao */
      SUBROUTINE AOTOAO(IAOAO,IPREAD)
C
C     This subroutine sets up pointer IAOAO which convertes between
C     two different orderings of AO's: from the ordering
C     of MOLECULE to an ordering in which the outer loop is over
C     atoms. (This orderings are identical when no symmetry is used.)
C
C     The purpose of this ordering is to make the AO's appear in the
C     same order as in the corresponding run with no symmetry.
C
C     tuh 120988
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
      DIMENSION IAOAO(MXCORB)
C
#include "symmet.h"
#include "nuclei.h"
#include "shells.h"
C
C
      IF (IPREAD .GE. 10) CALL TITLER('Output from AOTOAO','*',103)
C
C     Loop over all atoms
C
      IAO   = 0
      JATOM = 0
      DO 100 IATOM = 1, NUCIND
         ISTABA = ISTBNU(IATOM)
         DO 200 ISYMOP = 0, MAXREP
         IF (IAND(ISYMOP,ISTABA) .EQ. 0) THEN
            JATOM = JATOM + 1
C
C           Loop over all orbitals
C
            JAO = 0
            DO 300 ISHELL = 1, KMAX
               ICENT  = NCENT(ISHELL)
               ISTABO = ISTBAO(ISHELL)
               DO 400 ICMP = 1, KHKT(ISHELL)
                  DO 500 JSYMOP = 0, MAXREP
                  IF (IAND(JSYMOP,ISTABO) .EQ. 0) THEN
                     JAO = JAO + 1
                     JCENT = NUCNUM(ICENT,JSYMOP+1)
                     IF (JATOM .EQ. JCENT) THEN
                        IAO = IAO + 1
                        IAOAO(JAO) = IAO
                     END IF
                  END IF
  500             CONTINUE
  400          CONTINUE
  300       CONTINUE
C
C           End loop over orbitals
C
         END IF
  200    CONTINUE
  100 CONTINUE
C
C     End loop over atoms
C
      IF (IPREAD .GE. 10) THEN
         CALL HEADER('I - IAOAO',6)
         DO 600 I = 1, IAO
            WRITE (LUPRI,'(4X,I5,1X,I5)') I, IAOAO(I)
  600    CONTINUE
      END IF
      RETURN
      END
C/* Deck cntao */
      SUBROUTINE CNTAO(IPRINT)
C 
C     tuh
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "symmet.h"
#include "nuclei.h"
#include "shells.h"
C
C
      IORB = 0
      DO 100 LA = 0, MAXREP
      DO 200 IA = 1, KMAX
         DO 300 NA = 1, KHKT(IA)
            IVARB = IEOR(LA,ISYMAO(NHKT(IA),NA))
            IF (IAND(ISTBAO(IA),IVARB) .EQ. 0) THEN
               IORB = IORB + 1
               ICNTAO(IORB)  = NCENT(IA)
            END IF
  300    CONTINUE
  200 CONTINUE
  100 CONTINUE
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ICNTAO',-1)
         DO 400 I = 1, IORB
            WRITE (LUPRI,'(4X,I5,1X,I5)') I, ICNTAO(I)
  400    CONTINUE
      END IF
      RETURN
      END
C  /* Deck magcor */
      SUBROUTINE MAGCOR(IPREAD)
C
C     Printing of symmetrized nuclear magnetic moments
C
C     tuh March 1991
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "chrxyz.h"
#include "chrsgn.h"
C
      INTEGER NBASE(8), NSIGN(8), NRTREP(0:7)
C
#include "cbiher.h"
#include "nuclei.h"
#include "symmet.h"
#include "pgroup.h"
C
C
      IF (NMRISS .AND. MAXREP .GT. 0) THEN
         IF (IPREAD .GE. 0) THEN
            WRITE (LUPRI,'(/)')
            CALL HEADER('Symmetrized nuclear magnetic moments',1)
         END IF
         DO 100 IREP = 0, MAXREP
            NRTREP(IREP) = 0
            DO 110 ICENT = 1, NUCIND
               DO 120 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(ICENT - 1) + ICOOR,IREP,2)
                  IF (ISCOOR .GT. 0) NRTREP(IREP) = NRTREP(IREP) + 1
 120           CONTINUE
 110        CONTINUE
 100     CONTINUE
         IF (IPREAD .GE. 0) THEN
            WRITE (LUPRI,'(1X,A,8I3)')
     &           ' Number of components in each symmetry: ',
     &           (NRTREP(I),I=0,MAXREP)
         END IF
C
         DO 200 IREP = 0, MAXREP
         IF (NRTREP(IREP) .GT. 0) THEN
            IF (IPREAD .GE. 0) THEN
              WRITE (LUPRI,'(//2X,A,2X,A3,A1,I2,A1/)')
     &           'Symmetry',REP(IREP),'(',IREP,')'
            END IF
            DO 210 ICENT = 1, NUCIND
               DO 220 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(ICENT - 1) + ICOOR,IREP,2)
                  IF (ISCOOR .GT. 0) THEN
                     IVAR = IEOR(IREP,ISYMAX(ICOOR,2))
                     NB = 0
                     DO 230 ISYMOP = 0, MAXOPR
                     IF (IAND(ISYMOP,ISTBNU(ICENT)) .EQ. 0) THEN
                        NB = NB + 1
                        NSIGN(NB) = NINT(PT(IAND(ISYMOP,IVAR)))
                        NBASE(NB) = 3*(NUCNUM(ICENT,ISYMOP+1)-1) + ICOOR
                     END IF
 230                 CONTINUE
                     IF (IPREAD .GE. 0) THEN
                        WRITE (LUPRI,
     &                  '(2X,A,I3,3X,A,2X,A,3X,I2,7(2X,A,1X,I2))')
     &                  'I',ISCOOR, NAMEX(3*ICENT)(1:4), CHRXYZ(-ICOOR),
     &                  NBASE(1), (CHRSGN(NSIGN(I)),NBASE(I),I=2,NB)
                     END IF
                  END IF
  220          CONTINUE
  210       CONTINUE
         END IF
  200    CONTINUE
      END IF
      RETURN
      END
C  /* Deck tracr */
      SUBROUTINE TRACR(NCOOR)
C
C     Set up IPTCOR - points from symmetry coordinate to generating
C     Cartesian coordinate. Extracted from old TRACOR routine by
C     K.Ruud, Dec-96
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "symmet.h"
C
      DO 400 ITYPE = 1, 2
         IJ = 0
         DO 410 J = 0, MAXREP
            DO 420 I = 1, NCOOR
               IF (IPTCNT(I,J,ITYPE) .NE. 0) THEN
                  IJ = IJ + 1
                  IPTCOR(IJ,ITYPE) = I
               END IF
  420       CONTINUE
  410    CONTINUE
  400 CONTINUE
      RETURN
      END
C  /* Deck setclb */
      SUBROUTINE SETCLB(I,LABINT,NDIM)
#include "implicit.h"
      CHARACTER LABINT(NDIM)*4
      IF (I .EQ. 0) THEN
         LABINT(1) = 's   '
      ELSE IF (I .EQ. 1) THEN
         LABINT(2) = 'px  '
         LABINT(3) = 'py  '
         LABINT(4) = 'pz  '
      ELSE IF (I .EQ. 2) THEN
         LABINT(5)  = 'dxx '
         LABINT(6)  = 'dxy '
         LABINT(7)  = 'dxz '
         LABINT(8)  = 'dyy '
         LABINT(9)  = 'dyz '
         LABINT(10) = 'dzz '
      ELSE IF (I .EQ. 3) THEN
         LABINT(11)  = 'fxxx'
         LABINT(12)  = 'fxxy'
         LABINT(13)  = 'fxxz'
         LABINT(14)  = 'fxyy'
         LABINT(15)  = 'fxyz'
         LABINT(16)  = 'fxzz'
         LABINT(17)  = 'fyyy'
         LABINT(18)  = 'fyyz'
         LABINT(19)  = 'fyzz'
         LABINT(20)  = 'fzzz'
      END IF
      RETURN
      END
C  /* Deck spdcar */
      FUNCTION SPDCAR(IORDER)
#include "implicit.h"
      CHARACTER*1 SPDCAR
      IF (IORDER .EQ. 0) THEN
         SPDCAR = 's'
      ELSE IF (IORDER .EQ. 1) THEN
         SPDCAR = 'p'
      ELSE IF (IORDER .EQ. 2) THEN
         SPDCAR = 'd'
      ELSE IF (IORDER .LE. 6) THEN
         SPDCAR = CHAR(ICHAR('f') + IORDER - 3)
      ELSE
C
C        Note: j is not used - therefore this special case
C
         SPDCAR = CHAR(ICHAR('k') + IORDER - 7)
      END IF
      RETURN
      END
C  /* Deck sphinp */
C aspg, 2006-05-12
C updating to dalton's version
C
      SUBROUTINE SPHINP(LUINFO,WORK,LWORK,DOOWN,MXSPD)
C=======================================================================
C     Calculate spherical harmonic transformation coefficients
C     or, if DOOWN true, read user specified coefficients.
C
C     Modifications to DIRAC by Trond Saue Aug 18 2006
C
C     Miro,Sept.2017 - parallel job was setting LUINFO value to 0 for GNU compilers 
C     upon fixing I also made routine "implicit none" style
C
C=======================================================================
      implicit none
#include "priunit.h"
      real*8, parameter :: D1 = 1.0D0
#include "mxcent.h"
#include "maxaqn.h"
#include "sphtrm.h"
#include "molinp.h"
#include "ccom.h"
#include "cbirea.h"
      integer, intent(in) :: LUINFO
      CHARACTER*1 KEY
      LOGICAL DOOWN
      integer :: LWORK,i,j,IOFF,IOS,KCKI,KHKI,MXSPD
      real*8 :: WORK(LWORK)

C
      IF (DOOWN) THEN
         WRITE(LUPRI,'(/A)')
     &   '  Reading user specified transformation matrices for GTOs.'
      ELSE IF (IPREAD .GT. 4) THEN
         WRITE(LUPRI,'(/A)')
     &   '  Calculation of transformation matrices for spherical GTOs.'
      END IF
C
!Miro      integer MXQN is in include/maxaqn.h
      MXSPD = MXQN
C
      IOFF = 1
      DO I = 1, MXQN
        KCK(I)    = I*(I+1)/2
        ISPADR(I) = IOFF
        IOFF      = IOFF + KCK(I)**2
        IF (DOOWN) THEN
          SPH(I) = .TRUE.
        ELSE
          IF (DOCART) THEN
            KHK(I) = I*(I+1)/2
            SPH(I) = .FALSE.
          ELSE
            KHK(I) = 2*I - 1
            SPH(I) = I .GT. 2
          END IF
        END IF
      END DO
C
      SPHNRM = .TRUE.

C
C     Always make Cartesian transformation matrix
C     
      CALL DZERO(CSP,NCSP)
      IF (.NOT.DOOWN) THEN
        CALL SPHLAB(MXQN-1,GTOTYP)
        CALL GETSPH(WORK,LWORK,IPREAD)
      ELSE
        CALL RDLINE(LUINFO)
        READ (MLINE(NMLINE),'(BN,I5)',IOSTAT=IOS) MXSPD
        IF (IOS.NE.0) THEN
           WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
           WRITE(LUPRI,*) MLINE(NMLINE)
           CALL QUIT('Error in reading...1jsj38')
        ENDIF

        IF (MXSPD .GT. MXQN) THEN
          WRITE (LUPRI,'(/A/I10,A,I4//A//A/A)')
     &    ' ERROR: specified MXSPD for GTO transformation matrix',
     &    MXSPD,' is greater than allowed max value of',MXQN,
     &    ' Dump of last two lines read from .mol file:',
     &    MLINE(NMLINE-1),MLINE(NMLINE)
          CALL QUIT('ERROR in reading .mol file, see output')
        END IF
        IOFF = 1
        DO 300 I = 1, MXSPD
          CALL RDLINE(LUINFO)
          READ (MLINE(NMLINE),'(BN,4X,A1,I5)',IOSTAT=IOS)  KEY, KHK(I)
          IF (IOS.NE.0) THEN
             WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
             WRITE(LUPRI,*) MLINE(NMLINE)
             CALL QUIT('Error in reading...aj2sjjsn dcn')
          ENDIF
          IF (KHK(I) .GT. KCK(I)) THEN
            WRITE (LUPRI,'(/A,I5,A,I3/A,I5,A//A//A)')
     &        ' ERROR: you specified',KHK(I),
     &        ' components for GTO transformation matrix for L =',
     &         I-1,' and that is greater than',KCK(I),
     &         ', the number of Cartesian components.',
     &         ' Dump of last line read from .mol file:',MLINE(NMLINE)
               CALL QUIT('ERROR in reading .mol file, see output')
          END IF
          CALL RDCSP(LUINFO,CSP(ISPADR(I)),GTOTYP(IOFF),
     &               KHK(I),KCK(I))
          IOFF = IOFF + KHK(I)
  300   CONTINUE
        CALL HEADER('Cartesian transformation matrices',-1)
        WRITE (LUPRI,'(A/)') '  to user specified combinations'
        IOFF = 0
        DO 400 I = 1, MXSPD
          KHKI = KHK(I)
          KCKI = KCK(I)
          WRITE (LUPRI,'(A,I5//A/,(5X,15(A4,1X)))')
     &      '  Coefficients for angular quantum number ',I-1,
     &      '  to GTOs with labels:',(GTOTYP(IOFF+J),J=1,KHKI)
          IOFF = IOFF + KHKI
          CALL OUTPUT(CSP(ISPADR(I)),1,KHKI,1,KCKI,KHKI,KCKI,1,
     &                LUPRI)
 400    CONTINUE
        CALL FLSHFO(LUPRI)
      END IF
      IF (DOCART) CALL CARLAB(MXQN-1,GTOTYP)

      RETURN
      END
C  /* Deck rdcsp */
      SUBROUTINE RDCSP(LUINFO,CSP,GTOTYP,KHKI,KCKI)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "molinp.h"
      integer, intent(in) :: LUINFO
      CHARACTER*4 GTOTYP(*)
      DIMENSION CSP(KHKI,KCKI)
      DO 100 I = 1, KHKI
         CALL RDLINE(LUINFO)
         READ (MLINE(NMLINE),'(1X,A,15F5.3)',IOSTAT=IOS)
     &          GTOTYP(I), (CSP(I,J),J=1,KCKI)
         IF (IOS.NE.0) THEN
           WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
           WRITE(LUPRI,*) MLINE(NMLINE)
           CALL QUIT('Error in reading...sajksxn2773')
         ENDIF
  100 CONTINUE
      RETURN
      END
C  /* Deck rdline */
      SUBROUTINE RDLINE(LUINFO)
C
C     Purpose: Read a line from MOLECULE.INP
C
C     We keep track of which line we read last time (NMLN_LASTREAD).
C     This allows some routines to read a line, check if it want to use
C     that line, and then discard it by decrementing NMLINE.
C     Note that we only allow decrements of one /jth-20000621.
C     (This is used to check if SMALL basis keyword has been specified,
C     or if default - kinetic balance - should be used.)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "molinp.h"
#include "infpar.h"

      PARAMETER (len_NEXT_INPUT_LINE = 120)
      CHARACTER*(len_NEXT_INPUT_LINE) NEXT_INPUT_LINE
      INTEGER    NMLN_LASTREAD
      DATA       NMLN_LASTREAD /0/
      SAVE       NMLN_LASTREAD
      INTEGER :: NMLINE_STORED(1:2)=0
      integer, intent(in) :: LUINFO
C
C     If NMLINE is zero, reset NMLN_LASTREAD since this is a new call of readin,
C     e.g. during a geometry optimization.
C
      IF ( NMLINE .EQ. 0) NMLN_LASTREAD = 0
C
   10 CONTINUE
      NMLINE = NMLINE + 1

      NMLINE_STORED(1) = NMLINE_STORED(2)
      NMLINE_STORED(2) = NMLINE

      IF ( NMLINE .LT. NMLN_LASTREAD ) THEN
         WRITE(LUPRI,'(//A,2(/A,I5))')
     &        '*** ERROR ***: NMLINE < NMLN_LASTREAD in RDLINE',
     &        'NMLINE = ',NMLINE,
     &        'NMLN_LASTREAD = ',NMLN_LASTREAD
         CALL QUIT('*** ERROR reading .mol file')
      ENDIF

      IF (NMLINE.GT.NMLN_LASTREAD.AND
     &   .NMLINE_STORED(1).NE.NMLINE_STORED(2)) THEN

         READ (LUINFO,'(A)',END=12,IOSTAT=IOS) MLINE(NMLINE)

!miro print lengthy error message for debug purposes
         IF (IOS.NE.0) THEN
          WRITE(LUPRI,*) 
     &    'ERROR(subroutine RDLINE): Error in reading molecule.MOL. '
     &    //'NODE=',MYTID
          WRITE(LUPRI,*) 'IOS=',IOS 
          WRITE(LUPRI,*) 'NMLINE, line no:', NMLINE
          WRITE(LUPRI,*) 'NMLN_LASTREAD',NMLN_LASTREAD
          WRITE(LUPRI,*) 'read MLINE(NMLINE)=',MLINE(NMLINE)
          WRITE(LUPRI,*) 'LEN(MLINE(NMLINE))=',LEN(MLINE(NMLINE))
          WRITE(LUPRI,*) 'NMLINE_STORED(1),(2):',
     &    NMLINE_STORED(1),NMLINE_STORED(2)
          IF (NMLINE.GT.1) THEN
            WRITE(LUPRI,*) 'Previous line:',MLINE(NMLINE-1)
          ENDIF
         ENDIF
      ENDIF

      IF (NMLINE.GT.KMLINE) THEN
         WRITE (LUPRI,'(//A,/A,I5)')
     &      ' READIN ERROR: Too many lines in input -,',
     &      ' Increase dimension in molinp.h to at least',NMLINE
         CALL QUIT('*** ERROR *** Too many input lines in .mol file')
      END IF

      IF((MLINE(NMLINE)(1:1) .EQ. '#') .OR.
     $   (MLINE(NMLINE)(1:1) .EQ. '!') .OR.
     $   (MLINE(NMLINE)(1:1) .EQ. '$')) GOTO 10

      NMLN_LASTREAD = NMLINE

   12 CONTINUE
      RETURN
      END
C  /* Deck wronel */
      SUBROUTINE WRONEL(TITLE,IQM,IFXYZ,KATOM,JCO2,KANG)
C
C     Write interface records
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
#include "maxorb.h"
#include "maxaqn.h"
C
      CHARACTER*72 TITLE(2)
      CHARACTER*8 CDATE,CTIME
      CHARACTER*8 TABLE1(4), TABLE2(4), TABLE3(4), TABLE4(4), TABLE5(4)
      DIMENSION QPOL(6), QQ(3), IFXYZ(3), JFXYZ(3)
      DIMENSION CHRN(MXCENT), COOO(MXCENT,3), IPRXYZ(MXCORB),
     &          INAMN(MXCORB), IGTO(MXCORB), RTITLE(24)
      DIMENSION IQM(KATOM), JCO2(KANG,KATOM)
      DATA TABLE1 /'********','        ','        ','ISORDK  '/,
     &     TABLE2 /'********','        ','        ','SCFINP  '/,
     &     TABLE3 /'********','        ','        ','SYMINPUT'/,
     &     TABLE4 /'********','        ','        ','TESTDATA'/,
     &     TABLE5 /'********','        ','        ','SPHERICA'/
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#include "hrunit.h"
#include "nuclei.h"
#include "frame.h"
#include "sphtrm.h"
C
      CALL GETDAT(CDATE,CTIME)
      TABLE1(2) = CDATE
      TABLE2(2) = CDATE
      TABLE3(2) = CDATE
      TABLE4(2) = CDATE
      TABLE5(2) = CDATE
      TABLE1(3) = CTIME
      TABLE2(3) = CTIME
      TABLE3(3) = CTIME
      TABLE4(3) = CTIME
      TABLE5(3) = CTIME
C
      READ (TITLE(1),'(12A6)',IOSTAT=IOS) (RTITLE(I),I= 1,12)
      READ (TITLE(2),'(12A6)',IOSTAT=IOS) (RTITLE(I),I=13,24)
C
      ICENT = 0
      DO 100 N = 1, NUCIND
      DO 100 IREP = 0,MAXREP
         IF (IAND(ISTBNU(N),IREP) .EQ. 0) THEN
            ICENT = ICENT + 1
            CHRN(ICENT)   = CHARGE(N)
            COOO(ICENT,1) = PT(IAND(ISYMAX(1,1),IREP))*CORD(1,N)
            COOO(ICENT,2) = PT(IAND(ISYMAX(2,1),IREP))*CORD(2,N)
            COOO(ICENT,3) = PT(IAND(ISYMAX(3,1),IREP))*CORD(3,N)
      END IF
  100 CONTINUE
C
      I = 0
      DO 200 ISYMOP = 0, MAXREP
         NBI = NAOS(ISYMOP+1)
         DO 220 L = 1,NBI
            I = I + 1
            READ(  NAMN(IPCEN(I)),'(A4)',IOSTAT=IOS) INAMN(I)
            READ(GTOTYP(IPTYP(I)),'(A4)',IOSTAT=IOS) IGTO(I)
  220    CONTINUE
  200 CONTINUE
C
C     0 - unlabeled first record
C     ==========================
C
      REWIND LUONEL
      WRITE (LUONEL) RTITLE,MAXREP+1,(NAOS(I),I=1,MAXREP+1), POTNUC
#ifndef PRG_DIRAC
C
C     1 - ISORDK
C     ==========
C
      WRITE (LUONEL) TABLE1
      WRITE (LUONEL) DUMMY, DUMMY, DUMMY, DUMMY
      WRITE (LUONEL) CHRN, NUCDEP, COOO
C
C     2 - SCFINP
C     ==========
C
      WRITE (LUONEL) TABLE2
      WRITE (LUONEL) RTITLE,
     &               MAXREP+1, (NAOS(I), I = 1,MAXREP+1),
     &               POTNUC,
     &               KMAX, (NCENT(I), I = 1,KMAX),
     &               NBASIS, ( JTRAN(I),             I = 1,NBASIS),
     &                       ((ITRAN(I,J), J = 1,8), I = 1,NBASIS),
     &                       ((CTRAN(I,J), J = 1,8), I = 1,NBASIS),
     &               NBASIS, ( INAMN(I),             I = 1,NBASIS),
     &                       (IPRXYZ(I),             I = 1,NBASIS),
     &               DPNUC,
     &               NUCDEP, ((COOO(I,J),J=1,3),I=1,NUCDEP),
CTROND: JFXYZ, QPOL and QQ is not defined !!!
     &               IFXYZ, DUMMY, QPOL, QQ, JFXYZ
#endif
C
C     3 - SYMINPUT
C     ============
C
      WRITE (LUONEL) TABLE3
      WRITE (LUONEL) NBASIS,(INAMN(I),I=1,NBASIS),(IGTO(I),I=1,NBASIS),
     &               POTNUC,DUMMY,DUMMY,DUMMY
C
C     4 - TESTDATA
C     ============
C
C     Information for testing magnetic integrals
C
      WRITE (LUONEL) TABLE4
      WRITE (LUONEL) NONTYP, (NONT(I), I = 1, NONTYP),
     &               (IQM(I), (JCO2(I,J),J=1,IQM(I)), I = 1, NONTYP),
     &               (IDUMMY,I=1,5)
C
C     5 - SPHERICA
C     ============
C
C     Information for testing spherical transformation
C
      IDOCRT = 1
      IF (.NOT.DOCART) IDOCRT = 0
      WRITE (LUONEL) TABLE5
      WRITE (LUONEL) KMAX, IDOCRT, (NHKT(I),I = 1,KMAX), (IDUMMY,I=1,5)
      WRITE (LUONEL) MAXREP, (ISTBAO(I),I=1,KMAX), (IDUMMY,I=1,6)
      WRITE (LUONEL) MXQN, MXAQN, ((ISYMAO(I,J),I=1,MXQN),J=1,MXAQN)
      WRITE (LUONEL) (CSP(I),I=1,NCSP),(ISPADR(I),I=1,MXQN)
C
      RETURN
      END
C  /* Deck reapri */
      SUBROUTINE REAPRI
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
C
#include "abainf.h"
#include "cbirea.h"
C
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "chrsgn.h"
#include "aosotr.h"
C
      CALL HEADER('Test output from REAPRI',-1)
      WRITE (LUPRI,'(A,1P,D15.6)') ' THRS  ', THRS
      WRITE (LUPRI,'(A,I5)') ' NHTYP ', NHTYP
      WRITE (LUPRI,'(A,L5)') ' DOCART', DOCART
      print *, 'nucind,nucdep',NUCIND, NUCDEP
      print *, 'charge ',(CHARGE(i),i=1,nucind)
      print *, 'cord1',(cord(1,i),i=1,nucind)
      print *, 'cord2',(cord(2,i),i=1,nucind)
      print *, 'cord2',(cord(3,i),i=1,nucind)
      print *, 'noorbt ',(noorbt(i),i=1,nucind)
      print *, 'nucpre ',(nucpre(i),i=1,nucind)
      print *, 'nucnum ',((NUCNUM(i,j),i=1,nucind),j=1,8)
      print *, 'nucdeg ',(nucdeg(i),i=1,nucind)
      print *, 'istbnu ',(istbnu(i),i=1,nucind)
      print *, 'ntraco,itraco',NTRACO, (ITRACO(i),i=1,3)
      print *, 'NATOMS, NFLOAT, NBASIS, NPBAS',
     &          NATOMS, NFLOAT, NBASIS, NPBAS
      print *, 'namex ',(namex(i),i=1,3*nucind)
      print *, 'namdep',(namdep(i),i=1,nucdep)
      print *, 'namdpx',(namdpx(i),i=1,nucdep)
      call header('priexp',-1)
      call output(priexp,1,1,1,mxprim,1,mxprim,1,lupri)
      call header('priccf',-1)
      call output(priccf,1,mxprim,1,mxcont,mxprim,mxcont,1,lupri)
      call header('pricrx',-1)
      call output(pricrx,1,mxprim,1,3,mxprim,3,1,lupri)
      print *, ' kmax ', kmax
      call header('cent',-1)
      call output(cent,1,kmax,1,3,kmax,3,1,lupri)
      print *, ' nhkt ',(nhkt(i),i=1,kmax)
      print *, ' khkt ',(khkt(i),i=1,kmax)
      print *, ' kckt ',(kckt(i),i=1,kmax)
      print *, ' istbao ',(istbao(i),i=1,kmax)
      print *, ' nuco ',(nuco(i),i=1,kmax)
      print *, ' jstrt ',(jstrt(i),i=1,kmax)
      print *, ' nstrt ',(nstrt(i),i=1,kmax)
      print *, ' ncent ',(ncent(i),i=1,kmax)
      print *, ' share ',(share(i),i=1,kmax)
      print *, ' nrco ',(nrco(i),i=1,kmax)
      print *, ' numcf ',(numcf(i),i=1,kmax)
      print *, ' nbch ',(nbch(i),i=1,kmax)
      print *, ' kstrt ',(kstrt(i),i=1,kmax)
      print *, ' segm ',(segm(i),i=1,kmax)
      print *, ' iptshl ',(iptshl(i),i=1,kmax)
      print *, ' numcft ',(numcft(i),i=1,kmax)
      print *, ' sphr ',(sphr(i),i=1,kmax)
      print *, ' fmult ',(fmult(i),i=0,7)
      print *, ' pt ',(pt(i),i=0,7)
      print *, ' mult ',(mult(i),i=0,7)
      print *, ' maxrep, maxopr ',MAXREP, MAXOPR
      print *, ' isymax ',((isymax(i,j),i=1,3),j=1,2)
      print *, ' isymao ',((isymao(i,j),i=1,mxqn),j=1,mxaqn)
      print *, ' nparsu ',(nparsu(i),i=1,8)
      print *, ' naos ',(naos(i),i=1,8)
      print *, ' nparnu ',((nparnu(i,j),i=1,8),j=1,8)
      print *, ' iptsym ',((iptsym(i,j),i=1,norbs),j=0,7)
      print *, ' iptcnt ',(((iptcnt(i,j,k),i=1,3*nucdep),j=0,7),k=1,2)
      print *, ' ncrrep ',((ncrrep(i,j),i=0,7),j=1,2)
      print *, ' iptcor ',((iptcor(i,j),i=1,3*nucdep),j=1,2)
      print *, ' naxrep ',((naxrep(i,j),i=0,7),j=1,2)
      print *, ' iptax ',((iptax(i,j),i=1,3),j=1,2)
      print *, ' iptxyz ',(((iptxyz(i,j,k),i=1,3),j=0,7),k=1,2)
      print *, ' iptnuc ',((iptnuc(i,j),i=1,nucdep),j=0,7)
      call header('ctran',-1)
      call output(ctran,1,nbasis,1,8,mxcorb,8,1,lupri)
      print *, ' itran ',((itran(i,j),i=1,nbasis),j=1,8)
      print *, ' iaoao ',(iaoao(i),i=1,nbasis)
      return
      end
C  /* Deck zmat */
      SUBROUTINE ZMAT
C
C  Modified for incorporation in HERMIT/SIRIUS/ABACUS, K.Ruud, May-95
C  Based on program ZMAT from CCQC/ Nov. 1988 (H.J.Aa.Jensen)
C**********************************************************
C*** LAST UPDATED ON SEPTEMBER 11, 1986 BY R. REMINGTON ***
C* REASON: TO OBTAIN HIGHER ACCURACY BY INCREASING THE  ***
C* NUMBER OF DIGITS FOR PI AND CHANGE ATAN TO DATAN #143 **
C**********************************************************
C***LAST UPDATED ON FEBRUARY 01, 1985 BY YUKIO YAMAGUCHI***
C**********************************************************
C   THIS PROGRAM CALCULATES CARTESIAN COORDINATES OF A SYSTEM
C   SPECIFIED BY INTERNAL COORDINATES.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      LOGICAL FIRST
#include "nuclei.h"
#include "molinp.h"
      COMMON/CBZMAT/IZ(MXCENT,4),BL(MXCENT),ALP(MXCENT),BET(MXCENT)
      DATA FIRST /.TRUE./
      SAVE FIRST
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_NOFREE)
      DUMMY = 1.0D20
      IDUMMY = - 9 999 999
#endif
      NUCIND = NUCIND + 1
      IF (NUCIND .GT. MXCENT) GOTO 5000
      IZ(NUCIND,1) = 0
      IZ(NUCIND,2) = 0
      IZ(NUCIND,3) = 0
      IZ(NUCIND,4) = 0
      BL(NUCIND)=D0
      ALP(NUCIND)=D0
      BET(NUCIND)=D0
      READ (MLINE(NMLINE),'(A4)',IOSTAT=IOS) NAMN(NUCIND)
      IF (FIRST) THEN
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D)
         ISTART = 5
         CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT')
         CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA')
#else
         READ (MLINE(NMLINE)(5:80),*,IOSTAT=IOS) INUM, CHARGE(NUCIND)
         IF (IOS.NE.0) THEN
           WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
           WRITE(LUPRI,*) MLINE(NMLINE)
           CALL QUIT('Error in reading...12jsj')
         ENDIF
#endif
         FIRST = .FALSE.
         GOTO 204
      ELSE
         IF(NUCIND .EQ. 2) GO TO 201
         IF(NUCIND .EQ. 3) GO TO 202
         GO TO 203
      END IF
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D)
 201  ISTART = 5
      CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,1),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BL(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA')
      GO TO 204
 202  ISTART = 5
      CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,1),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BL(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,2),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,ALP(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA')
      GO TO 204
 203  ISTART = 5
      CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,1),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BL(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,2),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,ALP(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,3),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BET(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,4),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA')
#else
  201 READ (MLINE(NMLINE)(5:80),*,IOSTAT=IOS)
     &     INUM,IZ(NUCIND,1),BL(NUCIND),CHARGE(NUCIND)
      IF (IOS.NE.0) THEN
        WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
        WRITE(LUPRI,*) MLINE(NMLINE)
        CALL QUIT('Error in reading...2339jdcn')
      ENDIF
      GO TO 204
  202 READ (MLINE(NMLINE)(5:80),*,IOSTAT=IOS)
     &     INUM,IZ(NUCIND,1),BL(NUCIND),IZ(NUCIND,2),
     &     ALP(NUCIND),CHARGE(NUCIND)
      IF (IOS.NE.0) THEN
        WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
        WRITE(LUPRI,*) MLINE(NMLINE)
        CALL QUIT('Error in reading...2dcn')
      ENDIF
      GO TO 204
  203 READ (MLINE(NMLINE)(5:80),*,IOSTAT=IOS)
     &     INUM,IZ(NUCIND,1),BL(NUCIND),IZ(NUCIND,2),
     &     ALP(NUCIND),IZ(NUCIND,3),BET(NUCIND),
     &     IZ(NUCIND,4),CHARGE(NUCIND)
      IF (IOS.NE.0) THEN
        WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
        WRITE(LUPRI,*) MLINE(NMLINE)
        CALL QUIT('Error in reading...2xsydcn')
      ENDIF
#endif
 204  CONTINUE
      NMLINE = NMLINE + 1
      NCLINE(NUCIND) = NMLINE
      NONT(NUCIND) = 1
      NAMEX(3*NUCIND)     = NAMN(NUCIND)//' z'
      NAMEX(3*NUCIND - 1) = NAMN(NUCIND)//' y'
      NAMEX(3*NUCIND - 2) = NAMN(NUCIND)//' x'
C
      RETURN
C
C     Error message:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A/A,I5)')
     &    ' Too many atomic centers: MXCENT exceed in ZMAT,',
     &    ' Current limit:',MXCENT
        CALL QUIT('*** ERROR *** MXCENT exceeded in ZMAT')
      END
#ifndef PRG_DIRAC
C  /* Deck buildz */
      SUBROUTINE BUILDZ(IPRINT,NSYMOP,IQM,NBLCK,JCO,
     &                NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KBLOCK,
     &                KPRIM,DOOWN,UNCONT,BASNAM,BASREF)
C
C     Modified for HERMIT/SIRIUS/ABACUS -> K.Ruud, May-95
C
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "consts.h"
#include "pi.h"
      PARAMETER (DSM = 1.0D-30)
      PARAMETER (AR=PI/180.0D+00)
      PARAMETER (CORMAX = 1.D5)
      LOGICAL   UNCONT, NOORBTS, DOOWN
      COMMON/CBZMAT/IZ(MXCENT,4),BL(MXCENT),ALP(MXCENT),BET(MXCENT)
      DIMENSION A(MXCENT),B(MXCENT),CZ(MXCENT,3),D(MXCENT)
      DIMENSION U1(3),U2(3),U3(3),U4(3),VJ(3),VP(3),V3(3)
      DIMENSION IQM(KATOM,2),JCO(KANG,KATOM,2),
     &          NBLCK(KATOM,2),NUC(KBLOCK,2),NRC(KBLOCK,2),
     &          SEG(KBLOCK,2),
     &          ALPHA(KPRIM,KBLOCK,2),CPRIM(KPRIM,KPRIM,KBLOCK,2),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,2), CHNON(MXCENT)
      CHARACTER SPDCAR*1, NAVNK*4, CRT*1
      CHARACTER*80 BASREF(10,KATOM), BASNAM
#include "nuclei.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "ccom.h"
#include "molinp.h"
#include "cbisol.h"
    1 FORMAT(/'   Internal coordinates (Z matrix)'/)
    2 FORMAT(I5,62X,F15.5)
    3 FORMAT(I5,I5,F14.7,43X,F15.5)
    4 FORMAT(I5,I5,F14.7,I5,F14.7,24X,F15.5)
    5 FORMAT(I5,I5,F14.7,I5,F14.7,I5,F14.7,I5,F15.5)
C   6 FORMAT(//25X,' CARTESIAN COORDINATES in Angstrom',
C    &        /14X,'X',24X,'Y',24X,'Z'/)
    7 FORMAT(I5,3(F20.10,5X),F15.5)
    8 FORMAT(//25X,' Cartesian coordinates in bohr',
     &        /14X,'X',24X,'Y',24X,'Z'/)
    9 format(3f20.10)
C
      IF (SOLVNT) NUCIND = NUCIND - 1
C
      IF (IPRINT .GT. 2) THEN
         WRITE(LUPRI,1)
         WRITE(LUPRI,2) 1,CHARGE(1)
         IF (NUCIND .LE. 1) GOTO 201
         WRITE(LUPRI,3) 2,IZ(2,1),BL(2),CHARGE(2)
         IF(NUCIND.LE.2) GO TO 201
         WRITE(LUPRI,4) 3,IZ(3,1),BL(3),IZ(3,2),ALP(3),CHARGE(3)
         IF(NUCIND.LE.3) GO TO 201
         DO I=4,NUCIND
            WRITE(LUPRI,5) I,IZ(I,1),BL(I),IZ(I,2),ALP(I),
     &                     IZ(I,3),BET(I),IZ(I,4),CHARGE(I)
         END DO
 201  CONTINUE
      END IF
      DO J = 1,3
         DO I = 1,NUCIND
            CZ(I,J) = D0
         END DO
      END DO
      DO I = 1,NUCIND
         ALP(I) = ALP(I)*AR
         BET(I) = BET(I)*AR
      END DO
      CZ(2,3) = BL(2)
      IF (NUCIND .GE. 3) THEN
         CZ(3,1) = BL(3)*DSIN(ALP(3))
         IF ((IZ(3,1)-1) .EQ. 0) THEN
            CZ(3,3) = BL(3)*DCOS(ALP(3))
         ELSE
            CZ(3,3) = CZ(2,3)-BL(3)*DCOS(ALP(3))
         END IF
         DO I = 4,NUCIND
            IF((ABS(CZ(I-1,1))-1.0D-6) .LT. D0) THEN
               CZ(I,1) = BL(I)*DSIN(ALP(I))
               ITEMP = IZ(I,1)
               JTEMP = IZ(I,2)
               CZ(I,3) = CZ(ITEMP,3) - BL(I)*DCOS(ALP(I))*
     &                   DSIGN(D1,CZ(ITEMP,3) - CZ(JTEMP,3))
            ELSE
               GOTO 90
            END IF
         END DO
 90      CONTINUE
         K=I
         IF (K .LE. NUCIND) THEN
            DO 250 J = K,NUCIND
               CAJ = COS(ALP(J))
               SAJ = SIN(ALP(J))
               CBJ = COS(BET(J))
               SBJ = SIN(BET(J))
               IF(IZ(J,4) .EQ. 0) THEN
                  CALL V3VEC(U1,CZ,IZ(J,2),IZ(J,3))
                  CALL V3VEC(U2,CZ,IZ(J,1),IZ(J,2))
                  CALL V3PROD(VP,U1,U2)
                  R = SQRT(D1 - (U1(1)*U2(1) + U1(2)*U2(2)
     &                 + U1(3)*U2(3))**2)
                  DO I = 1,3
                     U3(I) = VP(I)/R
                  END DO
                  CALL V3PROD(U4,U3,U2)
                  DO I = 1,3
                     VJ(I) = BL(J)*(-U2(I)*CAJ + U4(I)*SAJ*CBJ
     &                    + U3(I)*SAJ*SBJ)
                     ITEMP = IZ(J,1)
                     CZ(J,I)=VJ(I)+CZ(ITEMP,I)
                  END DO
               ELSE
                  IF (IABS(IZ(J,4)) .EQ. 1) THEN
                     CALL V3VEC(U1,CZ,IZ(J,1),IZ(J,3))
                     CALL V3VEC(U2,CZ,IZ(J,2),IZ(J,1))
                     ZETA = -(U1(1)*U2(1) + U1(2)*U2(2) +U1(3)*U2(3))
                     A(J) = (-CBJ + ZETA*CAJ)/(D1 - ZETA*ZETA)
                     B(J) = (CAJ - ZETA*CBJ)/(D1 - ZETA*ZETA)
                     R=D0
                     GAMMA=PI/D2
                     IF (ZETA .LT. D0) THEN
                        R=PI
                        GAMMA = DATAN(SQRT(D1-ZETA*ZETA)/ZETA)+R
                     ELSE IF (ZETA .GT. D0) THEN
                        GAMMA = DATAN(SQRT(D1-ZETA*ZETA)/ZETA)+R
                     END IF
                     D(J) = D0
                     IF ((ABS(GAMMA+ALP(J)+BET(J)-D2*PI)-1.0D-6)
     &                    .GE. D0) THEN
                        D(J) = IZ(J,4)*(SQRT(D1+A(J)*CBJ-B(J)*CAJ))
     &                       /SQRT(D1-ZETA*ZETA)
                     END IF
                     CALL V3PROD(V3,U1,U2)
                     DO I = 1,3
                        U3(I) = A(J)*U1(I)+B(J)*U2(I)+D(J)*V3(I)
                        VJ(I)=BL(J)*U3(I)
                        ITEMP=IZ(J,1)
                        CZ(J,I)=VJ(I)+CZ(ITEMP,I)
                     END DO
                  ELSE
                     CALL V3VEC(U1,CZ,IZ(J,1),IZ(J,3))
                     CALL V3VEC(U2,CZ,IZ(J,2),IZ(J,1))
                     ZETA = -(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
                     CALL V3PROD(V3,U1,U2)
                     V3MAG = SQRT(V3(1)*V3(1)+V3(2)*V3(2)
     &                    + V3(3)*V3(3))
                     A(J) = V3MAG*CBJ/(D1-ZETA*ZETA)
                     B(J) = SQRT((D1-CAJ*CAJ-A(J)*CBJ*V3MAG)
     &                    /(D1-ZETA*ZETA))
                     IF((IZ(J,4)-2) .NE. 0) THEN
                        B(J)=-B(J)
                     END IF
                     D(J) = B(J)*ZETA+CAJ
                     DO I = 1,3
                        U3(I) = B(J)*U1(I)+D(J)*U2(I)+A(J)*V3(I)
                        VJ(I) = BL(J)*U3(I)
                        ITEMP = IZ(J,1)
                        CZ(J,I) = VJ(I)+CZ(ITEMP,I)
                     END DO
                  END IF
               END IF
 250        CONTINUE
         END IF
      END IF
C
C   ELIMINATE IMAGINARY ATOM(S)
C   Hmmm, what do we do here? Remove them or keep them as dummy atoms?
C   Check to see if there is a basis set associated with the center?
C   Is it at all possible to place a dummy atom on top of an already
C   existing center in internal coordinates? K.Ruud, May-95
C
C      NATOMS=0
C      DO 290 I=1,NUCIND
C         IF(ANZ(I) .NE. D0) THEN
C            NATOMS=NATOMS+1
C            CHRG(NATOMS) = DFLOAT( ANZ(I) )
C            DO 280 J=1,3
C               C(NATOMS,J)=CZ(I,J)
C 280        CONTINUE
C         END IF
C 290  CONTINUE
C
C
      IF (IPRINT .GT. 2) THEN
         WRITE(LUPRI,8)
         DO I=1,NUCIND
            WRITE(LUPRI,7) I,(CZ(I,J)/XTANG,J=1,3),CHARGE(I)
         END DO
      END IF
C
C     Move information to HERMIT common blocks
C
      DO I=1,NUCIND
         DO J = 1, 3
            CORD(J,I) = CZ(I,J)/XTANG
            IF (ABS(CORD(J,I)).GT. CORMAX) GOTO 5010
         END DO
      END DO
C
C     Update MOLECULE.INP. In order to work properly with the symmetry
C     detection routines we need to sort atoms with similar charge
C     together in a block.
C
      ITYP = 1
      NONT(ITYP) = 1
      CHNON(ITYP) = CHARGE(1)
      DO 600 I = 2, NUCIND
         DO 601 J = I-1, 1, -1
            IF (CHARGE(I) .EQ. CHARGE(J)) THEN
               NONTYP = NONTYP - 1
               DO L = 1, ITYP
                  IF (CHARGE(I) .EQ. CHNON(L)) NONT(L) = NONT(L) + 1
               END DO
               CHARGK = CHARGE(I)
               XCORD  = CORD(1,I)
               YCORD  = CORD(2,I)
               ZCORD  = CORD(3,I)
               NAVNK  = NAMN(I)
               DO K = I - 1, J + 1, -1
                  CHARGE(K + 1) = CHARGE(K)
                  CORD(1,K + 1) = CORD(1,K)
                  CORD(2,K + 1) = CORD(2,K)
                  CORD(3,K + 1) = CORD(3,K)
                  NAMN(K + 1)   = NAMN(K)
               END DO
               CHARGE(J + 1) = CHARGK
               CORD(1,J + 1) = XCORD
               CORD(2,J + 1) = YCORD
               CORD(3,J + 1) = ZCORD
               NAMN(J + 1)   = NAVNK
               GOTO 600
            END IF
 601     CONTINUE
         ITYP = ITYP + 1
         NONT(ITYP) = 1
         CHNON(ITYP) = CHARGE(I)
 600  CONTINUE
C     define atom number in IZATOM /hjaaj Mar 2004
C     TODO : define multiple basis sets and point charges
      DO I = 1, NUCIND
         IZATOM(I) = NINT(CHARGE(I))
      END DO
      IF (ITYP.NE.NONTYP) CALL QUIT('Problem ITYP.ne.NONTYP in BUILDZ')
      IF (DOCART) THEN
         CRT = 'C'
      ELSE IF (DOOWN) THEN
         CRT = 'X'
      ELSE
         CRT = ' '
      END IF
      WRITE(MLINE(5),'(A1,I4,I3,I2,10X,1P,D10.2)')
     &     CRT,NONTYP,KCHARG,NSYMOP,THRS
      NMLINE = 6
      NHTYP  = 0
      JBLOCK = 1
      IATOM = 1
      DO 605 I = 1, NONTYP
         NOORBTS = .TRUE.
         WRITE(MLINE(NMLINE),'(6X,F5.0,I4)') CHARGE(IATOM),NONT(I)
         NMLINE = NMLINE + 1
         DO 606 J = 1, NONT(I)
            NCLINE(IATOM) = NMLINE
            WRITE(MLINE(NMLINE)(1:4),'(A4)') NAMN(IATOM)
            WRITE(MLINE(NMLINE)(5:24),'(F20.16)')  CORD(1,IATOM)
            WRITE(MLINE(NMLINE)(25:44),'(F20.16)') CORD(2,IATOM)
            WRITE(MLINE(NMLINE)(45:64),'(F20.16)') CORD(3,IATOM)
            WRITE(MLINE(NMLINE)(65:80),'(A16)') '                '
            NMLINE = NMLINE + 1
            IATOM = IATOM + 1
 606     CONTINUE
C
C     Add basis set to newly order atomic set
C
         KAOVEC = KBLOCK + 1 - JBLOCK
         Q = CHARGE(IATOM-1)
         CALL BASLIB(IQM(I,1),JCO(1,I,1),NUC(JBLOCK,1),
     &               NRC(JBLOCK,1),SEG(JBLOCK,1),
     &               ALPHA(1,JBLOCK,1),
     &               CPRIM(1,1,JBLOCK,1),CPRIMU(1,1,JBLOCK,1),
     &               NBLOCK,KAOVEC,KPRIM,Q,Q,DSM,
     &               UNCONT,BASNAM,BASREF(1,I),IPREAD)
         NBLCK(I,1) = NBLOCK
         JBLOCK = JBLOCK + NBLCK(I,1)
         NOORBTS = NOORBTS .AND. IQM(I,1).EQ.0
         IF (IQM(I,1) .GT. 0) THEN
            NHTYP = MAX(NHTYP,IQM(I,1))
            IF (NHTYP .GT. MXQN) GOTO 5000
         END IF
         DO K = 1,NONT(I)
            NOORBT(K + IATOM - 2) = NOORBTS
         END DO
 605  CONTINUE
      NMLINE = NMLINE - 1
      IF (SOLVNT) NUCIND = NUCIND + 1
      RETURN
 5000 CONTINUE
         WRITE (LUPRI,'(6X,A,I3,3A/9X,2A/9X,2(A,I3),A)')
     &      '*  Input specifies highest orbital of atomic type ',
     &      I,' AS "',SPDCAR(NHTYP - 1),'".',
     &      ' Highest allowed orbital in this version: ',
     &      SPDCAR(MXQN - 1),
     &      ' Increase MXQN from',MXQN,' to',NHTYP,' and recompile.'
         CALL QUIT('Too high angular quantum no. specified in input.')
 5010 CONTINUE
        WRITE (LUPRI,'(A,1P,E12.6,A/A/A,E12.6)')
     &    ' Atomic coordinate ',CORD(J,NUCIND),
     &    ' too large in BUILDZ',
     &    ' Note: Program is unstable for large coordinate values.',
     &    ' Maximum coordinate value:',CORMAX
        CALL QUIT('*** ERROR: Atomic coordinate too large in BUILDZ')
      END
C  /* Deck V3vec */
      SUBROUTINE V3VEC(U,C,J,K)
#include "implicit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      DIMENSION C(MXCENT,3),R(3),U(3)
C
      R2=D0
      DO 101 I=1,3
         R(I)=C(J,I)-C(K,I)
         R2=R2+R(I)*R(I)
 101  CONTINUE
      R2=SQRT(R2)
      DO 102 I=1,3
         U(I)=R(I)/R2
  102 CONTINUE
      RETURN
      END
C  /* Deck v3prod */
      SUBROUTINE V3PROD(VP,X,Y)
#include "implicit.h"
      DIMENSION VP(3),X(3),Y(3)
C
      VP(1)=X(2)*Y(3)-X(3)*Y(2)
      VP(2)=X(3)*Y(1)-X(1)*Y(3)
      VP(3)=X(1)*Y(2)-X(2)*Y(1)
      RETURN
      END
#endif /* PRG_DIRAC */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Getpot */
      FUNCTION GETPOT(Q,IPRINT)
C***********************************************************************
C
C     Get nuclear potential energy
C     Surgery by T.Saue June 7 1997
C     Point charges are not included. /June2000,jth+hjaaj
C     Yes, but not if too close, and not pc-pc (pc=point charge). /May2010,hjaaj
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
      PARAMETER (DIST_IDENT = 0.01D0)
      ! if DIST .lt. DIST_IDENT then same center and no self-repulsion.
      ! Could be used for small component charge modeled as a point
      ! charge.
      ! Maybe also useful for getting basis set from an atom with
      ! another charge (e.g. for equivalent atom model in XPS). Example: 
      ! to use nitrogen basis on oxygen, first specify nitrogen and
      ! then a point charge of +1 with the same coordinates.
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "symmet.h"
#include "nuclei.h"
#include "infpar.h"
      DIMENSION Q(*)
C
      CALL QENTER('GETPOT')
      N_IDENT   = 0
      NTOTAL   = 0
      EPOT_QM_QM = D0 ! QM is charge with orbitals
      EPOT_QM_PC = D0 ! PC is point charge without orbitals
      EPOT_PC_PC = D0
      DO 10 N = 1, NUCIND
         DO 20 M = N, NUCIND
            Q_MN = Q(M)*Q(N)
         IF (Q_MN .EQ. D0) GO TO 20
            DO 30 KB = 0, MAXREP
            IF (IAND(KB,ISTBNU(M)) .NE. 0) GO TO 30
            IF (IAND(KB,ISTBNU(N)) .NE. 0) GO TO 30
            IF (M.EQ.N .AND. KB.EQ.0) GO TO 30 ! no self-repulsion
               DS = (CORD(1,N)-CORD(1,M)*PT(IAND(ISYMAX(1,1),KB)))**2
     &            + (CORD(2,N)-CORD(2,M)*PT(IAND(ISYMAX(2,1),KB)))**2
     &            + (CORD(3,N)-CORD(3,M)*PT(IAND(ISYMAX(3,1),KB)))**2
               DIST = SQRT(DS)
               IF (DIST.LT.DIST_IDENT) THEN
                  IF ( .NOT.(NOORBT(M).OR.NOORBT(N)) ) GOTO 5000
C                  ... exit if close and no floating and no point charge
                  N_IDENT = N_IDENT + 1
                  GO TO 30
               END IF
               HKAB  = Q_MN * FMULT(IAND(ISTBNU(M),ISTBNU(N))) / DIST
               IF (M .EQ. N) HKAB = DP5*HKAB
               IF (NOORBT(N) .AND. NOORBT(M)) THEN
                  EPOT_PC_PC = EPOT_PC_PC + HKAB
               ELSE IF (NOORBT(M) .OR. NOORBT(N)) THEN
                  EPOT_QM_PC = EPOT_QM_PC + HKAB
               ELSE
                  EPOT_QM_QM = EPOT_QM_QM + HKAB
               END IF
   30       CONTINUE ! DO 30 KB = 0, MAXREP
   20    CONTINUE ! DO 20 M = N, NUCIND
   10 CONTINUE ! DO 10 N = 1, NUCIND
C
#ifdef VAR_MPI
!     print only from one process
      if (mytid == 0) then
#endif
        IF (IPRINT .GE. 0) THEN


           if (      (dabs(epot_pc_pc) == 0.0d0)
     &         .and. (dabs(epot_qm_pc) == 0.0d0)
     &        ) then
              WRITE(LUPRI,'(4(/A,F20.12))')
     &        '  QM-QM nuclear repulsion energy :',EPOT_QM_QM
           else if ( (dabs(epot_qm_qm) == 0.0d0)
     &         .and. (dabs(epot_qm_pc) == 0.0d0)
     &         .and. (dabs(epot_pc_pc) == 0.0d0)
     &        ) then
!             don't print anything if there is
!             only one atom and no point charges
           else
              WRITE(LUPRI,'(4(/A,F20.12))')
     &        '  QM-QM nuclear repulsion energy :',EPOT_QM_QM,
     &        '  QM-pc interaction energy       :',EPOT_QM_PC,
     &        '  pc-pc interaction energy       :',EPOT_PC_PC,
     &        '  Total interaction energy       :',
     &              EPOT_QM_QM+EPOT_QM_PC+EPOT_PC_PC
              IF (N_IDENT .GT. 0) WRITE (LUPRI,*) N_IDENT,
     &    ' interactions neglected because between same physical center'
           end if
        END IF
#ifdef VAR_MPI
      end if
#endif
      GETPOT = EPOT_QM_QM + EPOT_QM_PC
      CALL QEXIT('GETPOT')
      RETURN
C
 5000 CONTINUE
      WRITE (LUPRI,'(//3(A,I5),A,D15.5/,2(/A,4F12.5))')
     &     'ERROR: Nucleus',N,' is too close to the',KB,
     &     'th transformation of nucleus no.',M,' DIST=',DIST,
     &     'Charge annd coordinates of first  center:',Q(N),CORD(1:3,N),
     &     'Charge annd coordinates of second center:',Q(M),CORD(1:3,M)
      CALL QUIT('GETPOT: Nuclei too close.')
      END
C/*  Deck frefrm */
      SUBROUTINE FREFRM(STRING,ISTART,IVAL,RVAL,TYPE)
C
C     Simulate free format input for Cray computers when reading from an
C     internal file, K.Ruud-Feb.97
C
!
!    IVAL - entering INTORB (no of contracted shells for given L)
!
! MI: Added debugging output (June,2009)
!
#include "implicit.h"
#include "priunit.h"
#include "cbirea.h"
C
      CHARACTER*80 STRING
      CHARACTER TYPE*3, FRMT*7
C
      CALL QENTER('FREFRM')
C
      IF (IPREAD.GE.10) THEN
        CALL HEADER('Output from FREFRM',-1)
        WRITE(LUPRI,'(A,A)') 'STRING=',STRING
        WRITE(LUPRI,*) 'entering ISTART=',ISTART
        WRITE(LUPRI,*) 'entering IVAL=',IVAL
        WRITE(LUPRI,*) 'TYPE=',TYPE
      ENDIF

      IF (TYPE .EQ. 'INT') THEN
         IVAL = 0
      ELSE
         RVAL = 0.0D0
      END IF
      DO 10 IPOS = ISTART, 80
         ! leave cycle when found first non-white character
         IF (STRING(IPOS:IPOS) .NE. ' ') THEN
            GO TO 20
         ENDIF
 10   CONTINUE
C
 20   CONTINUE
      IEND = INDEX(STRING(IPOS:),' ')
      ILEN = IEND - ISTART - 1 + IPOS
      IF (IPREAD.GE.10) THEN
        WRITE(LUPRI,*) 'IPOS=',IPOS
        WRITE(LUPRI,*) 'IEND=',IEND
        WRITE(LUPRI,*) 'ILEN=',ILEN
      ENDIF

      IF (TYPE .EQ. 'INT') THEN
         IF (ILEN .GT. 9) THEN
            WRITE (FRMT,'(2X,A2,I2,A1)') '(I',ILEN,')'
         ELSE
            WRITE (FRMT,'(3X,A2,I1,A1)') '(I',ILEN,')'
         END IF

!mi: crucial - filling IVAL !
         READ (STRING(ISTART:80),FRMT,IOSTAT=IOS) IVAL

         IF (IOS.NE.0.OR.IPREAD.GE.10) THEN
           IF(IOS.NE.0) WRITE(LUPRI,*) 'Error in reading of IVAL !'
           WRITE(LUPRI,*) 'IOS=',IOS
           WRITE(LUPRI,*) 'ISTART=',ISTART
           WRITE(LUPRI,*) 'ILEN=',ILEN
           WRITE(LUPRI,*) 'IVAL=',IVAL
           WRITE(LUPRI,'(a,a)') 'determined FRMT=',FRMT
           WRITE(LUPRI,'(a,a)') 'STRING(ISTART:80)=',STRING(ISTART:80)
           WRITE(LUPRI,'(a,a)') 'STRING=',STRING
           IF (IOS.NE.0) CALL QUIT('FREFRM: Error in reading of IVAL')
         ENDIF
      ELSE IF (TYPE .EQ. 'REA') THEN
         IF (ILEN .GT. 9) THEN
            WRITE (FRMT,'(A2,I2,A3)') '(F',ILEN,'.0)'
         ELSE
            WRITE (FRMT,'(A3,I1,A3)') ' (F',ILEN,'.0)'
         END IF
         READ (STRING(ISTART:80),FRMT,IOSTAT=IOS) RVAL
         IF (IOS.NE.0) THEN
           WRITE(LUPRI,*) 'Error reading string section : "',
     &          STRING(ISTART:80),'"'
           WRITE(LUPRI,*) 'Whole string is "',STRING,'"'
           CALL QUIT('Error in reading...RVAL')
         ENDIF
      ELSE
         WRITE (LUPRI,'(/A)') 'Illegal variable type in FREFRM : ',TYPE
         CALL QUIT('Illegal variable type in FREFRM')
      END IF

!      update ISTART !
      !IF (IPOS.LT.80) THEN ! mi: this does not help !!!
        ISTART = IEND + IPOS
      ! ENDIF

      IF (IPREAD.GE.10) THEN
        WRITE(LUPRI,*) 'updated ISTART=',ISTART
      ENDIF

      CALL QEXIT('FREFRM')
      RETURN
      END
C  /* Deck MASS2ISOTOP */
      INTEGER FUNCTION MASS2ISOTOP(ICHARG,MASSNM)
C
C     Function to switch from mass number to isotope number sorted
C     according to abundance, K.Ruud-02
C
#include "implicit.h"
#include "priunit.h"
C
      CALL QENTER('MASS2ISOTOP')
      IORD = 0
      DO I = 1, 5
         MASS_I = NINT( DISOTP(ICHARG,I,'MASS') )
         IF (MASS_I .EQ. MASSNM) IORD = I
      END DO
      IF (IORD .EQ. 0) THEN
         WRITE (LUPRI,'(/A,I4,A,I4)') 'ERROR: unknown isotope mass',
     &        MASSNM,' for atom with charge',ICHARG
         CALL QUIT('Unknown isotope mass for chosen atomic charge')
      ELSE
         MASS2ISOTOP = IORD
      END IF
      CALL QEXIT('MASS2ISOTOP')
      RETURN
      END
C  /* Deck genfamily */
      SUBROUTINE GENFAMEXP(LUINFO,IFLAG)
C***********************************************************************
C
C     Read in or calculate exponents for family basis sets
C
C     Written by J. Thyssen - Jun 4 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER ( D0 = 0.00 D00 , D1 = 1.00 D00 )
      CHARACTER KEYWRD*6
C
#include "cbirea.h"
#include "molinp.h"
C
C
C
      IF ( IFLAG .EQ. 0 ) THEN
C
C        Generate well-temperered exponents
C        Geometric progression is just a special case
C        with gamma = 0 and delta = 0.
C
C        See Huzinaga and Klobukowski, CPL 212, 260 (1993).
C
         IF (NFAMEXP(1) .LT. 1 .OR. NFAMEXP(1) .GT. MXFAMEXP) GOTO 5000
         IF (FAMPAR(1) .LE. D0) GOTO 5010
         IF (FAMPAR(2) .LE. D1) GOTO 5010
         IF (FAMPAR(3) .LT. D0) GOTO 5010
         IF (FAMPAR(4) .LT. D0) GOTO 5010
C
         IF ( INPTST ) THEN
            WRITE(LUPRI,'(A,/,4(A,F20.10,/),A,I3)')
     $           'Parameters for well-tempered progression: ',
     $           'alpha = ',FAMPAR(1),
     $           'beta  = ',FAMPAR(2),
     $           'gamma = ',FAMPAR(3),
     $           'delta = ',FAMPAR(4),
     $           'N     = ',NFAMEXP(1)
         END IF
C
C        Formula is:
C
C        zeta_N = alpha,
C
C        zeta_N-k+1 = zeta_N-k+2 beta ( 1 + gamma (k/N)^delta), k = 2,N
C
         FAMEXP(NFAMEXP(1), 1) = FAMPAR(1)
         DO K = 2,NFAMEXP(1)
            FAMEXP(NFAMEXP(1) - K + 1, 1) =
     $           FAMEXP(NFAMEXP(1) - K + 2, 1) * FAMPAR(2) *
     $           ( D1 + FAMPAR(3) *
     $           ( ( dble(K) / dble(NFAMEXP(1)) ) ** FAMPAR(4) ) )
         END DO
C
C        Copy into set 2.
C
         NFAMEXP(2) = NFAMEXP(1)
         DO K = 1, NFAMEXP(1)
            FAMEXP( K, 2 ) = FAMEXP( K, 1 )
         END DO
C
      ELSE IF ( IFLAG .EQ. 1 ) THEN
C
C        Read exponents from file
C
         DO I = 1,NFAMEXP(1)
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),*,IOSTAT=IOS) FAMEXP(I, 1)
            IF (IOS.NE.0) THEN
              WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
              WRITE(LUPRI,*) MLINE(NMLINE)
              CALL QUIT('Error in reading...askc233n')
            ENDIF
         END DO
C
C        Copy into set 2.
C
         NFAMEXP(2) = NFAMEXP(1)
         DO K = 1, NFAMEXP(1)
            FAMEXP( K, 2 ) = FAMEXP( K, 1 )
         END DO
C
      ELSE IF ( IFLAG .EQ. 2 ) THEN
C
C        Read dual set of exponents from file
C
         DO I = 1,NFAMEXP(1)
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),*,IOSTAT=IOS) FAMEXP(I, 1)
            IF (IOS.NE.0) THEN
              WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
              WRITE(LUPRI,*) MLINE(NMLINE)
              CALL QUIT('Error in reading...3daskc233n')
            ENDIF
         END DO
         DO I = 1,NFAMEXP(2)
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),*,IOSTAT=IOS) FAMEXP(I, 2)
            IF (IOS.NE.0) THEN
              WRITE(LUPRI,*) 'Error reading line',NMLINE,':'
              WRITE(LUPRI,*) MLINE(NMLINE)
              CALL QUIT('Error in reading...3cxzx33n')
            ENDIF
         END DO
C
      ELSE
         WRITE(LUPRI,'(/,A,/,A,I3)')
     $        '*** ERROR in GENFAMEXP ***:',
     &        'Unknown IFLAG: ',IFLAG
         CALL QUIT('Unknown IFLAG in GENFAMEXP')
      END IF
C
      IF (INPTST) THEN
         CALL HEADER('Exponents for family basis set(1)',-1)
         CALL OUTPUT(FAMEXP,1,NFAMEXP(1),1,1,NFAMEXP(1),1,1,LUPRI)
         IF ( IFLAG .EQ. 2 ) THEN
            CALL HEADER('Exponents for family basis set(2)',-1)
            CALL OUTPUT(FAMEXP(1,2),1,NFAMEXP(2),1,1,NFAMEXP(2),1,
     &           -1,LUPRI)
         END IF
      END IF
      RETURN
C
 5000 CONTINUE
      WRITE(LUPRI,'(/A/,2(A,I3/))')
     &     '*** ERROR in GENFAMEXP *** ',
     $     'Illegal value for number of exponents: ',NFAMEXP(1),
     &     'Value must be between 1 and ',MXFAMEXP
      CALL QUIT('Illegal value for number of exponents')
 5010 CONTINUE
      WRITE(LUPRI,'(2(A,/),4(A,F20.10,/))')
     &     '*** ERROR in GENFAMEXP *** ',
     $     'Illegal non-positive value for parameters in'//
     &     ' well-tempered progression',
     &     'Alpha = ',FAMPAR(1),
     &     'Beta  = ',FAMPAR(2),
     &     'Gamma = ',FAMPAR(3),
     &     'Delta = ',FAMPAR(4)
      CALL QUIT('Illegal value for parameter in ' //
     $     'well-tempered progression')
      END
C  /* Deck PRBASREF */
      SUBROUTINE  PRBASREF (LUPRI,NONTYP,BASREF)
C***********************************************************************
C
C     Print out the reference information for basis sets
C     Sets that are used for more than one atom type are printed
C     only once.
C
C     Written by L. Visscher -  May 14, 2003
C
C***********************************************************************
#include "mxcent.h"
      DIMENSION I_SET(MXCENT),I_EQUAL(MXCENT)
      LOGICAL FIRST(MXCENT), BASREFS_EQUAL
      CHARACTER*80 BASREF(10,*)
C
C     Check for common sets, the first set is always unique
C
      FIRST(1) = .TRUE.
      DO I = 2, NONTYP
         FIRST(I) = .TRUE.
         DO J = 1, I - 1
            BASREFS_EQUAL = .TRUE.
            DO K = 1, 10
               IF (BASREF(K,I).NE.BASREF(K,J)) THEN
                   BASREFS_EQUAL = .FALSE.
                   EXIT
               ENDIF
            ENDDO
!        If equal is still true, then we have a match and can exit the loop over J
            IF (BASREFS_EQUAL) THEN
               FIRST(I) = .FALSE.
               I_SET(I) = J
               EXIT
            ENDIF
         ENDDO
      ENDDO
C
      WRITE (LUPRI,1000)
      DO I = 1, NONTYP
C        We will print the information only once
         IF (FIRST(I)) THEN
C           Find out which other types also use the same basis set
            N = 1
            I_EQUAL(N) = I
            DO J = I+1, NONTYP
               IF (I_SET(J).EQ.I) THEN
                  N = N + 1
                  I_EQUAL(N) = J
               ENDIF
            ENDDO
C           print line with info about the equal sets
            WRITE (LUPRI,1001) (I_EQUAL(J),J=1,N)
C           print basis set information
            DO ILINE = 1, 10
               IF (BASREF(ILINE,I)(1:16).NE.' Not initialized')
     &         WRITE (LUPRI,1002) BASREF(ILINE,I)
            ENDDO
         ENDIF
      ENDDO
C
 1000 FORMAT (//"  References for the basis sets"/
     &          "  -----------------------------")
 1001 FORMAT (/"  Atom type",16I4,(/11X,16I4))
 1002 FORMAT (2X,A80)
      RETURN
      END

C  /* Deck mnfout */
      SUBROUTINE MNFOUT(IQM,NBLCK,JCO,NUC,NRC,SEG,
     &     KATOM,KANG,KBLOCK,KPRIM,CPRIMU)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Now also dumps a file MNF.INP for use for SO mean-field approximation
C     Written by K.Ruud 160297 
C
C     Added and modified by Miro Ilias for the AMFI part
C     Strasbourg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "ccom.h"
#include "cbirea.h"
#include "nuclei.h"
#include "primit.h"
#include "molinp.h"

#include "infpar.h"

      LOGICAL SEG
C     CHARACTER KASYM(3,3)*1, CRT*1
      DIMENSION IQM(KATOM),NBLCK(KATOM),
     &     JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &     CPRIMU(KPRIM,KPRIM,KBLOCK),SEG(KBLOCK)

C     .... only for master or single process !!! 
      IF (MYTID .EQ. MPARID) THEN
C
      LUTMP2 = -1
      CALL GPOPEN(LUTMP2,'MNF.INP',' ',' ','FORMATTED',IDUMMY,.FALSE.)
C    
      IPRIMD = 0
      IPRIM =  0
      ICENT  = 0
      JBLOCK = 0
      ISYCNT = 1

      write(lutmp2,*)
     & 'Formatted input for AMFI - decontracted basis sets !'
      write(lutmp2,*) NONTYP

      DO 100 I = 1, NONTYP
         DO N = 1, NONT(I)
            ISYCNT = ISYCNT + 1
         END DO
         MNCENT = ICENT + 1
C   .... prevent segmentation of L shells...
         DO 110 N = 1, NONT(I)
            ICENT = ICENT + 1
            NDEG  = NUCDEG(ICENT)
            KBCH  = JBLOCK
CMI         ... charge, highest angular momentum
!          skip AMFI contributions with l > 4 (h, i, ... functions)
!          stefan - nov 2011: simple fix/approximation for the time
!                             being...
           imax_l = min(IQM(I)-1,4)
           IF (N.EQ.1) THEN
            WRITE(LUTMP2,*) NINT(CHARGE(MNCENT)) 
            WRITE (LUTMP2,*) CHARGE(MNCENT),imax_l
!sk         WRITE (LUTMP2,*) CHARGE(MNCENT),IQM(I)-1
           ENDIF
            II = 0
            DO 200 J = 1, IQM(I)
!          IF (N.EQ.1) THEN
!            write(lupri,*)
!    &       '  >>>  momentum J=',J,' JCO(J,I)=',JCO(J,I)
!          ENDIF
CMI           .... segmentation of L blocks ...not desirable for AMFI...
               IF (JCO(J,I).GT.1) THEN
                KBCH_SAVE=KBCH
                NNUCT=0
                NNRCT=0
                DO KI=1, JCO(J,I)
                  KBCH = KBCH + 1
                  NNUCT = NNUCT + NUC(KBCH) 
                  NNRCT = NNRCT + NRC(KBCH) 
                ENDDO
                 KBCH = KBCH_SAVE
               ENDIF
               DO 210 K = 1, JCO(J,I)
                  KBCH = KBCH + 1
                  NNUC  = NUC(KBCH)
                  NNRC  = NRC(KBCH)
                  IF (NNUC .EQ. 0) GO TO 200
                IF (N.EQ.1) THEN
                 IF (JCO(J,I).EQ.1) THEN
                   if(j .lt. 6)then ! up to g functions
                     WRITE (LUTMP2,*) NNUC,NNRC
                   end if
!                   WRITE (lupri,*) 'K(1-JCO(J,I))=',K,
!    &                  '  -> NNUC,NNRC:',NNUC,NNRC
                 ELSE
                  IF (K.EQ.1) THEN
                    if(j .lt. 6)then ! up to g functions
                      WRITE (LUTMP2,*) NNUCT,NNRCT
!                     WRITE (lupri,*) 
!    &                 '  -> NNUCT,NNRCT:',NNUCT,NNRCT
                    end if
                  ENDIF
                 ENDIF
                ENDIF
                  ITYP = NHKOFF(J)
                  IPSTRT = IPRIM + 1
                  IPRIM =  IPRIM + NNUC
                  ITYP = ITYP + 1
                  IPRIMD = IPRIMD + 1
                  IF (N .EQ. 1) THEN
                     if(j .lt. 6)then ! up to g functions
                       DO 410 M = 1, NNUC
                         IPRIMD = IPRIMD + 1
                         WRITE (LUTMP2,'(F24.11)') PRIEXP(IPSTRT-1+M)
!                        WRITE (lupri,'(F20.10)') PRIEXP(IPSTRT-1+M)
 410                   CONTINUE
                     else
                       IPRIMD = IPRIMD + NNUC
#ifdef X2C_debug
                       WRITE (lupri,'(a,i4)') 'mnf input writing:'//
     &             ' skipped angular momentum functions > 5 (h,i,...)',
     &               NNUC
#endif
                     end if
                  END IF
 210           CONTINUE
 200        CONTINUE
 110     CONTINUE
         JBLOCK = JBLOCK + NBLCK(I)
 100  CONTINUE
C
      WRITE (LUTMP2,'(A)') 'END'
      CALL GPCLOSE(LUTMP2,'KEEP')

      CALL FLSHFO(LUPRI)

      ENDIF ! MYTID .eq. MPARID

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Getsph */
      SUBROUTINE GETSPH(WORK,LWORK,IPRINT)
C***********************************************************************
C!
C!     Get spherical transformation matrices up to LMAX
C!     Written by T.Saue May 1 1998 
C!     Rewrite by cut'n paste 2006 Trond Saue
C!
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER(D1 = 1.0D0)
#include "ccom.h"
#include "sphtrm.h"
#include "cbirea.h"
      DIMENSION WORK(LWORK)
#include "memint.h"
      CALL MEMGET('REAL',KCSPT,(MXQN*(MXQN + 1)/2)*(2*MXQN - 1),
     &             WORK,KFREE,LFREE)
      DO LVAL = 0, MXQN - 1
        NLM = 2*LVAL + 1
        NXYZ = (LVAL+1)*(LVAL+2)/2
        IF (LVAL .EQ. 0) THEN
          CSP(1)  = D1
        ELSE IF (LVAL .EQ. 1) THEN
          CSP (2) = D1
          CSP (6) = D1
          CSP(10) = D1
        ELSE
C         generate normalized s,p,d,f,g,h,i GTO's
          MINTEG = 2
          CALL SPHCOM(LVAL,WORK(KCSPT),NLM,NXYZ,1,MINTEG,
     &                WORK(KFREE),LFREE,IPREAD)
C.........transpose  CSPT(NXYZ,NLM) --> CSP(NLM,NXYZ)
C.........note, however, 
C.........that ISPADR gives offsets in blocks of NXYZ*NXYZ
C
          CALL MTRSP(NXYZ,NLM,WORK(KCSPT),NXYZ,
     &               CSP(ISPADR(LVAL+1)),NLM)
        END IF
      ENDDO
C
      IF (IPRINT .GT. 4) THEN
        CALL HEADER('Cartesian transformation matrices',-1)
        WRITE (LUPRI,'(A/)') '  to spherical harmonics'
        IOFF = 0
        DO I = 1, MXQN
          KHKI = 2*I - 1
          KCKI = KCK(I)
          WRITE (LUPRI,'(A,I5//A/,(5X,15(A4,1X)))')
     &      '  Coefficients for angular quantum number ',I-1,
     &      '  to GTOs with labels:',(GTOTYP(IOFF+J),J=1,KHKI)
          IOFF = IOFF + KHKI
          CALL OUTPUT(CSP(ISPADR(I)),1,KHKI,1,KCKI,KHKI,KCKI,1,
     &                LUPRI)
        ENDDO
        CALL FLSHFO(LUPRI)
      END IF
C
      CALL MEMREL('GETSPH',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
C     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Mdef */
      FUNCTION MDEF(L,ICMP)
C***********************************************************************
C
C     An integer function that returns the M-value for component ICMP
C     of angular momentum L of solid harmonics.
C
C     The main purpose of the function is to centralize the definition
C     of the Cartesian to spherical transformation.
C
C     Current definition is: -L,..,+L
C
C     Written by Trond Saue Aug 19 2006
C
C***********************************************************************
#include "priunit.h"
      IF(L.EQ.1) THEN
        GO TO (1,2,3), ICMP
 1      CONTINUE
          MDEF = 1
          RETURN
 2      CONTINUE
          MDEF =-1
          RETURN
 3      CONTINUE
          MDEF = 0
          RETURN
      ELSE
        MDEF = ICMP-L-1
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION KHK_SC(L,DELTA_L)
C***********************************************************************
C
C     An integer function that returns the degeneracy of the small
C     component for a given large component with angular momentum L
C
C     Note L is given here as L + 1 to start with one (instead of zero)
C     Isolated in a function to avoid using ccom.h directly
C
C     Written by Lucas Visscher, 26-6-2007
C
C***********************************************************************
#include "maxaqn.h"
#include "ccom.h"
      INTEGER DELTA_L
      IF (L+DELTA_L.GT.0) THEN
         KHK_SC = KHK(L+DELTA_L)
      ELSE
         KHK_SC = 0
      ENDIF
      RETURN
      END
C  /* Deck dkpro */
      SUBROUTINE MAKE_CNTMAT(IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           ALPHA,CPRIM,LPRIM,CONTC,KATOM,KANG,KBLOCK,KPRIM,
     &           IPRIMD,IORBD)
C******************************************************************************
C    
C     Set up contraction matrix to transform from decontracted to
C     contracted basis. Note that in decontracted basis the primitives
C     have been normalized, so this contribution has to be removed
C     from the normalization of the contracted basis.
C
C     This is DKCON, written by K. Ruud, 
C     and imported from DALTON (with polish) by T. Saue Jan 15 2008
C
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4 = 4.0D0, DP5 = 0.50D0,
     &           DP75 = 0.75D0)
C
#include "cbirea.h"
      LOGICAL SEG,SPHER
      DIMENSION IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          SEG(KBLOCK),ALPHA(KPRIM,KBLOCK),
     &          CPRIM(KPRIM,KPRIM,KBLOCK)
      DIMENSION CONTC(IPRIMD,IORBD)
      DIMENSION LPRIM(MAXPRD)
#include "ccom.h"
#include "nuclei.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#include "dcbgen.h"
C#include "gnrinf.h"
C
C
C     PIPPINV is the inverse of PIPPI (see routine NRMORB)
      PIPPINV = (D2*PI)**DP75
      ICENT  = 0
      JBLOCK = 0
      IPR    = 1
      IRB    = 1
      CALL DZERO(CONTC,IPRIMD*IORBD)
C.....loop over atomic types
      DO 10 I = 1,NONTYP
C........loop over symmetry independent atoms
         DO 20 N = 1,NONT(I)
            ICENT = ICENT + 1
            KBCH = JBLOCK
            NDEG = NUCDEG(ICENT)
C...........loop over angular quantum numbers
            DO 30 J = 1,IQM(I)
               KKK = 0
               NCOMP = KHK(J)
               NCCMP = KCK(J)
               SPHER = SPH(J)
C..............loop over AO-blocks for given l-value
               DO 40 K = 1, JCO(J,I)
                  KBCH = KBCH + 1
                  NCONT = NRC(KBCH)
C
C     Quick sort of primitive exponents if necessary
C
                  DO L1 = 1, NUC(KBCH)
                     LPRIM(L1) = L1
                  END DO
                  DO L1 = 1, NUC(KBCH) - 1
                     DO L2 = L1 + 1, NUC(KBCH)
                        IF(ALPHA(L2,KBCH) .GT. ALPHA(L1,KBCH)) THEN
                           IDUM = LPRIM(L2)
                           LPRIM(L2) = LPRIM(L1)
                           LPRIM(L1) = IDUM
                        END IF
                     END DO
                  END DO
C
C     Core for setting up contraction matrix
C
                  DO ICOMP = 1, NCOMP
                     IF (ICOMP .GT. 1) THEN
                        IPR = IPR - NDEG*(NCOMP*NUC(KBCH) - 1)
                        IRB = IRB - NDEG*NCONT*NCOMP + 1
                     END IF
                     DO LJ = 1, NUC(KBCH)
                        L = LPRIM(LJ)
                        IF (L .GT. 1) IRB = IRB - NCONT*NCOMP*NDEG
     &                       - (NDEG - 1)
                        PRFA  = PIPPINV*
     &                          ((D4*ALPHA(L,KBCH))**(-DP5*(J + DP5)))
                        DO IDEG = 1, NDEG
                           IF (IDEG .GT. 1) IRB =IRB
     &                          -NCONT*NCOMP*NDEG + 1
                           DO M = 1, NCONT
                              CONTC(IPR,IRB) = CPRIM(L,M,KBCH)*PRFA
                              IRB = IRB + NCOMP*NDEG
                           END DO
                           IPR = IPR + 1
                        END DO
                        IPR = IPR + NDEG*NCOMP - NDEG
                     END DO
                  END DO
                  IPR = IPR - NCOMP*NDEG + NDEG
                  IRB = IRB - NCOMP*NDEG + 1
C
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
         JBLOCK = JBLOCK + NBLCK(I)
 10   CONTINUE
      CALL OUTPUT(CONTC,1,IPRIMD,1,IORBD,IPRIMD,IORBD,1,LUPRI) 
C
C     Write contraction matrix to file
C
      LUCNMT = -9000
      CALL GPOPEN(LUCNMT,'CNTMAT','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      WRITE (LUCNMT,*) ((CONTC(I,J),I=1,IPRIMD),J=1,IORBD)
      CALL GPCLOSE(LUCNMT,'KEEP')
C
      RETURN
      END

      subroutine fix_duplicate_generators (nsymop,kasym)
      implicit none
      integer, intent(inout)    :: nsymop
      character, intent(inout)  :: kasym(3,3)*1
      integer :: i, j, k, igen12, igen(3)
#include "ibtfun.h"

C     We only need to check for 3 generators, 2 generators are always independent
      if (nsymop.lt.3) return

      do j = 1, nsymop
        igen(j) = 0
        do i = 1, 3
C         Convert the string into a bit representation (e.g. "XY"  becomes 3 = 011) 
          if (kasym(i,j).NE.' ') THEN
            k = ICHAR(kasym(i,j)) - ICHAR('W')
            igen(j) = igen(j) + 2**(k-1)
          end if
        end do
      end do

C     Check and correct for cases with redundant generators (generator 1 times generator 2 yields generator 3)
      igen12 = IBTXOR(igen(1),igen(2))
      if (igen12 .eq. igen(3)) then
          nsymop = nsymop - 1
      end if
      end

C --- end of herrdn.F ---
