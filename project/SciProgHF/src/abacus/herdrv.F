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

! file: abacus/herrdn.F

C  /* Deck herinp */
      SUBROUTINE HERINP(WORK,LWORK)
C
C --- General Input for HERMIT ---
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "efield.h"
C
      PARAMETER (NDIR = 6, NTABLE = 80)
      PARAMETER (D0 = 0.0D0)
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7
      DIMENSION WORK(LWORK)
C
#ifdef PRG_DIRAC
! dcbham.h : ECPCALC, NONREL
#include "dcbgen.h"
#include "dcbham.h"
#else
#include "gnrinf.h"
#endif
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
#include "nuclei.h"
#include "orgcom.h"
#include "huckel.h"
#include "ccom.h"
#include "hrunit.h"

#ifdef BUILD_GEN1INT
#include "gen1int.h"
#endif
C
      DATA TABDIR /'*END OF', '*READIN', '*ONEINT', '*TWOINT',
     &             '*SUPINT', '*XXXXXX'/
C     adds test suite of Gen1Int interface as .GENINT
C     Bin Gao, Nov. 9, 2011
      DATA TABLE  /'.PRINT ', '.INPTES', '.NOSUP ', '.SPIN-O',
     &             '.DIPLEN', '.NO HAM', '.SOTEST', '.DIPVEL',
     &             '.QUADRU', 'xXXXXXX', '.SECMOM', '.SUPONL',
     &             '.CARMOM', '.SPHMOM', '.FC    ', '.PSO   ',
     &             '.SD    ', '.DSO   ', '.POINTS', '.SELECT',
     &             '.QUASUM', '.SD+FC ', '.PROPRI', '.HDO   ',
     &             '.S1MAG ', '.S2MAG ', '.ANGMOM', '.ANGLON',
     &             '.LONMOM', '.MAGMOM', '.S1MAGT', '.MGMOMT',
     &             '.KINENE', '.S2MAGT', '.DSUSNL', '.DSUSLL',
     &             '.DSUSLH', '.DIASUS', '.DSUTST', '.NSTNOL',
     &             '.NSTLON', '.NST   ', '.NSNLTS', '.NSLTST',
     &             '.NELFLD', '.NSTTST', '.EFGCAR', '.EFGSPH',
     &             '.S1MAGL', '.S1MAGR', '.HDOBR ', '.S1MLT ',
     &             '.HDOBRT', '.S1MRT ', '.NUCPOT', '.NPOTST',
     &             '.MGMO2T', '.MGMTHR', '.HBDO  ', '.SUSCGO',
     &             '.NSTCGO', '.EXPIKR', '.MASSVE', '.DARWIN',
     &             '.CM-1  ', '.CM-2  ', '.SQHDOL', '.SQHDOR',
     &             '.NOTWO ', 'xXXXXXX', 'xXXXXXX', '.GFACDI',
     &             '.S1ELE ', '.S1ELB ', '.ONEELD', '.THETA ',
     &             '.NUCMOD', '.MAGCOR', '.GENINT', 'xXXXXXX'/
C
C 950602-vebjorn:
C     Added flag HRINPC to ensure that HERMIT input processing is done
C     only once.
C
      IPRDEF = IPRUSR + 1
      IF (HRINPC) GOTO 1000
      HRINPC = .TRUE.
C
C     Initialize /CBIHER/
C
      IF(RELCAL) THEN
        HAMILT = .FALSE.
        ONEPRP = .FALSE.
        DIPVEL = .FALSE.
        DIRAC  = .TRUE.
        SUPMAT = .FALSE.
      ELSE
        HAMILT = .TRUE.
        ONEPRP = .FALSE.
        DIPVEL = .FALSE.
        DIRAC  = .FALSE.
        SUPMAT = .NOT. DIRCAL
      ENDIF
      SKIP   = .FALSE.
      TSTINP = .FALSE.
      DIPLEN = .FALSE.
      QUADRU = .FALSE.
      SPNORB = .FALSE.
      SOTEST = .FALSE.
      NOTWO  = DIRCAL
      SECMOM = .FALSE.
      CARMOM = .FALSE.
      SPHMOM = .FALSE.
      OCTUPO = .FALSE.
      FERMI  = .FALSE.
      PSO    = .FALSE.
      SPIDIP = .FALSE.
      DSO    = .FALSE.
      SDFC   = .FALSE.
      PROPRI = .FALSE.
      HDO    = .FALSE.
      S1MAG  = .FALSE.
      S2MAG  = .FALSE.
      ANGMOM = .FALSE.
      ANGLON = .FALSE.
      LONMOM = .FALSE.
      MAGMOM = .FALSE.
      S1MAGT = .FALSE.
      MGMOMT = .FALSE.
      KINENE = .FALSE.
      S2MAGT = .FALSE.
      DSUSNL = .FALSE.
      DSUSLL = .FALSE.
      DSUSLH = .FALSE.
      DIASUS = .FALSE.
      DSUTST = .FALSE.
      NUCSNL = .FALSE.
      NUCSLO = .FALSE.
      NUCSHI = .FALSE.
      NSTTST = .FALSE.
      NSLTST = .FALSE.
      NELFLD = .FALSE.
      NSNLTS = .FALSE.
      EFGCAR = .FALSE.
      EFGSPH = .FALSE.
      S1MAGL = .FALSE.
      S1MAGR = .FALSE.
      HDOBR  = .FALSE.
      S1MLT  = .FALSE.
      S1MRT  = .FALSE.
      HDOBRT = .FALSE.
      NUCPOT = .FALSE.
      NPOTST = .FALSE.
      MGMO2T = .FALSE.
      HBDO   = .FALSE.
      SUSCGO = .FALSE.
      NSTCGO = .FALSE.
      MASSVL = .FALSE.
      DARWIN = .FALSE.
      CM1    = .FALSE.
      CM2    = .FALSE.
      SQHDOL = .FALSE.
      SQHDOR = .FALSE.
      S1ELE  = .FALSE.
      S1ELB  = .FALSE.
      ONEELD = .FALSE.
      GFACDI = .FALSE.
      THETA  = .FALSE.
      PRTHRS = 1.0D-10
      NPQUAD = 40
      ALLATM = .TRUE.
      TRIANG = .TRUE.
      EXPIKR = .FALSE.
      MGCOOR = .FALSE.
      EXPKR (1:3) = 0.0d0
      NPATOM = 0
      CALL IZERO(IPATOM,MXCENT)

      CALL TITLER('Output from HERMIT input processing','*',118)

C    Initialize nuclear model in NUCLEI
C      1 - Point nucleus
C      2 - Gaussian nucleus
C
      IF (ECPCALC) THEN
         GAUNUC = .FALSE.
      ELSE IF (RELCAL) THEN
        GAUNUC = .TRUE.
      ELSE
        GAUNUC = .FALSE.
      ENDIF
C
C     Initialize /CBISOL/ (10-Dec-92 th+hjaaj)
C     -- not used in Hermit; SOLVNT must be false
C        in order to skip solvent modules in ONEDRV
C        and cavity center in VIBMAS
C
      SOLVNT = .FALSE.
C
C     **** Find Hermit input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      CALL UPCASE(WORD)
      IF ((WORD .EQ. '**HERMI') .OR. (WORD .EQ. '*HERMIT') .OR.
     &    (WORD .EQ. '**INTEG')) THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         CALL REAINI(IPRDEF,RELCAL,TSTINP)
         WORD = '*END OF'
         GOTO 1
  920 CONTINUE
         CALL REAINI(IPRDEF,RELCAL,TSTINP)
         WORD = '*END OF'
         GOTO 1
  930 CONTINUE
      WORD1 = WORD
C
C     ***** Process input for COMMON  /CBIHER/  *****
C
      INUC = 0
  100 READ (LUCMD, '(A7)') WORD
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
     &                171,172,173,174,175,176,177,178,179,180), I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in HERINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in HERINP.')
  101    CONTINUE
            READ (LUCMD, '(I5)') IPRDEF
            GO TO 100
  102    CONTINUE
            TSTINP = .TRUE.
            GO TO 100
  103    CONTINUE
            SUPMAT = .FALSE.
            GO TO 100
  104    CONTINUE
            SPNORB = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  105    CONTINUE
            DIPLEN = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  106    CONTINUE
            HAMILT = .FALSE.
            SUPMAT = .FALSE.
            GO TO 100
  107    CONTINUE
            SOTEST = .TRUE.
            GO TO 100
  108    CONTINUE
            DIPVEL = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  109    CONTINUE
            QUADRU = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  110    CONTINUE
            GO TO 100
  111    CONTINUE
            SECMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  112    CONTINUE
            SUPMAT = .TRUE.
            HAMILT = .FALSE.
            NOTWO  = .TRUE.
            ONEPRP = .FALSE.
            GO TO 100
  113    CONTINUE
            CARMOM = .TRUE.
            READ (LUCMD,'(I5)') IORCAR
            ONEPRP = .TRUE.
            GO TO 100
  114    CONTINUE
            SPHMOM = .TRUE.
            READ (LUCMD,'(I5)') IORSPH
            ONEPRP = .TRUE.
            GO TO 100
  115    CONTINUE
            FERMI = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  116    CONTINUE
            PSO = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  117    CONTINUE
            SPIDIP = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  118    CONTINUE
            DSO = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  119    CONTINUE
            READ (LUCMD,'(I5)') NPQUAD
            GO TO 100
  120    CONTINUE ! .SELECT
            READ (LUCMD, *) NPATOM
            IF (NPATOM .GT. MXCENT) THEN
               WRITE (LUPRI,'(/A,/A,I6,/A,I6)')
     &             ' Too many atoms selected with .SELECT.',
     &             ' Number of atoms selected: ',NPATOM,
     &             ' Number of atoms allowed:  ',MXCENT
               CALL QUIT('Error in HERINP')
            END IF
            READ (LUCMD, *) (IPATOM(I),I=1,NPATOM)
            ALLATM = .FALSE.
            GO TO 100
  121    CONTINUE
            TRIANG = .FALSE.
            GO TO 100
  122    CONTINUE
            SDFC = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  123    CONTINUE
            PROPRI = .TRUE.
            GO TO 100
  124    CONTINUE
            HDO = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  125    CONTINUE
            S1MAG = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 126     CONTINUE
            S2MAG = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  127    CONTINUE
            ANGMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  128    CONTINUE
            ANGLON = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 129     CONTINUE
            LONMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 130     CONTINUE
            MAGMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 131     CONTINUE
            S1MAG  = .TRUE.
            S1MAGT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 132     CONTINUE
            MGMOMT = .TRUE.
            LONMOM = .TRUE.
            HAMILT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 133     CONTINUE
            KINENE = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 134     CONTINUE
            S2MAG  = .TRUE.
            S2MAGT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 135     CONTINUE
            DSUSNL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 136     CONTINUE
            DSUSLL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 137     CONTINUE
            DSUSLH = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 138     CONTINUE
            DIASUS = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 139     CONTINUE
            DSUTST = .TRUE.
            DSUSLL = .TRUE.
            ANGLON = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 140     CONTINUE
            NUCSNL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 141     CONTINUE
            NUCSLO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 142     CONTINUE
            NUCSHI = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 143     CONTINUE
            NSNLTS = .TRUE.
            NUCSLO = .TRUE.
            PSO    = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 144     CONTINUE
            NSLTST = .TRUE.
            NELFLD = .TRUE.
            NUCSNL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 145     CONTINUE
            NELFLD = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 146     CONTINUE
            NSTTST = .TRUE.
            NUCSLO = .TRUE.
            NUCSNL = .TRUE.
            NUCSHI = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 147     CONTINUE
            EFGCAR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 148     CONTINUE
            EFGSPH = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 149     CONTINUE
            S1MAGL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 150     CONTINUE
            S1MAGR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 151     CONTINUE
            HDOBR  = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 152     CONTINUE
            S1MLT  = .TRUE.
            S1MAGL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 153     CONTINUE
            HDOBR  = .TRUE.
            HDOBRT = .TRUE.
            DIPVEL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 154     CONTINUE
            S1MRT  = .TRUE.
            S1MAGR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 155     CONTINUE
            NUCPOT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 156     CONTINUE
            NUCPOT = .TRUE.
            HAMILT = .TRUE.
            NPOTST = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 157     CONTINUE
            MGMO2T = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 158     CONTINUE
            READ (LUCMD,*) PRTHRS
            GO TO 100
 159     CONTINUE
            HBDO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 160     CONTINUE
            SUSCGO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 161     CONTINUE
            NSTCGO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 162     CONTINUE
            EXPIKR = .TRUE.
            ONEPRP = .TRUE.
            READ (LUCMD,*) (EXPKR(I),I=1,3)
            GOTO 100
 163     CONTINUE
            MASSVL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 164     CONTINUE
            DARWIN = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 165     CONTINUE
            READ (LUCMD,'(A7)') FIELD1
            IF (.NOT. ((FIELD1 .EQ. 'X-FIELD')
     &          .OR. (FIELD1 .EQ. 'Y-FIELD')
     &          .OR. (FIELD1 .EQ. 'Z-FIELD'))) THEN
               WRITE (LUPRI,'(/,3A,/)') ' Field direction "',FIELD1,
     &               '" illegal'
               CALL QUIT('Illegal field directions for CM-1 integrals')
            END IF
            CM1    = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 166     CONTINUE
            READ (LUCMD,'(A7)') FIELD2
            IF (.NOT. ((FIELD2 .EQ. 'X-FIELD')
     &          .OR. (FIELD2 .EQ. 'Y-FIELD')
     &          .OR. (FIELD2 .EQ. 'Z-FIELD'))) THEN
               WRITE (LUPRI,'(/,3A,/)') ' Field direction "',FIELD2,
     &               '" illegal'
               CALL QUIT('Illegal field directions for CM-2 integrals')
            END IF
            CM2    = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 167     CONTINUE
            SQHDOL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 168     CONTINUE
            SQHDOR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 169     CONTINUE ! '.NOTWO '
            NOTWO = .TRUE.
            SUPMAT = .FALSE.
            GO TO 100
 170     CONTINUE
            GO TO 100
 171     CONTINUE
            GO TO 100
 172     CONTINUE ! '.GFACDI'
            GFACDI = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 173     CONTINUE
            S1ELE  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 174     CONTINUE
            S1ELB  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 175     CONTINUE
            ONEELD = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 176     CONTINUE
            THETA  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 177     CONTINUE
C&&& Nuclear model .NUCMOD
            READ (LUCMD, *,IOSTAT=IOS) INUC
            IF (IOS.NE.0) THEN
             CALL QUIT('Error in reading nuclear model INUC  !')
            ENDIF
            IF    (INUC.EQ.1) THEN
              GAUNUC = .FALSE.
            ELSEIF(INUC.EQ.2) THEN
              GAUNUC = .TRUE.
            ELSE
              WRITE(LUPRI,'(A,I5)')
     &         '*** ERROR *** Unknown nuclear model:',INUC
              CALL QUIT('*** ERROR *** HERINP: Unknown nuclear model')
            ENDIF
            GO TO 100
 178     CONTINUE
C&&& MAGCOR - print magnetic coordinates
            MGCOOR = .TRUE.
            GO TO 100
 179     CONTINUE
#ifdef BUILD_GEN1INT
            TEST_GEN1INT = .true.
#endif
            goto 100
 180     CONTINUE
            goto 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 199
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in '//WORD1)
      END IF
 199  CONTINUE
      NMRISS = FERMI.OR.PSO.OR.SPIDIP.OR.DSO.OR.SDFC.OR.MGCOOR
      WRITE (LUPRI,'(A,I5)') ' Default print level:    ',IPRDEF
      IF (TSTINP) WRITE (LUPRI,'(/A/)') ' Input test run only ...'
      IF (INUC .EQ. 0) THEN
         IF(GAUNUC) THEN
           WRITE(LUPRI,'(A)') ' Using default nuclear model: '//
     &    'Gaussian charge distribution.'
         ELSE
           WRITE(LUPRI,'(A)') ' Using default nuclear model: '//
     &     'Point charge.'
           IF (NONREL) WRITE(LUPRI,'(A)')
     &     ' (.NONREL calculations have point nucleus as default)'
           IF (ECPCALC) WRITE(LUPRI,'(A)')
     &     ' (.ECP calculations have point nucleus as default)'
         ENDIF
      ELSE
         IF(GAUNUC) THEN
           WRITE(LUPRI,'(A)') ' Nuclear model requested in input: '//
     &    'Gaussian charge distribution.'
         ELSE
           WRITE(LUPRI,'(A)') ' Nuclear model requested in input: '//
     &     'Point charge.'
         ENDIF
      END IF
      IF (SPNORB .AND. SOTEST) THEN
         WRITE (LUPRI,'(/A/)')
     *    ' Error in input: .SPIN-ORBIT and .SOTEST cannot '/
     *    /'both be specified.'
          CALL QUIT('Error in HERMIT input.')
      END IF
      IF (HAMILT) THEN
          IF (NOTWO) THEN
             WRITE (LUPRI,'(/A/)')
     &          ' Calculation of one-electron Hamiltonian integrals.'
          ELSE
             WRITE (LUPRI,'(/A/)')
     &    ' Calculation of one- and two-electron Hamiltonian integrals.'
          END IF
      ELSE
         IF (NOTWO) WRITE (LUPRI,'(/A/)')
     &          ' Two-electron integrals not calculated.'
      END IF

      IF (ONEPRP) THEN
         WRITE (LUPRI,'(/A/)')
     &  ' The following one-electron property integrals are calculated:'
                     WRITE (LUPRI,'(10X,A)') '- overlap integrals'
         IF (DIPLEN) WRITE (LUPRI,'(10X,A)') '- dipole length integrals'
         IF (DIPVEL) WRITE (LUPRI,'(10X,A)')
     &                         '- dipole velocity integrals'
         IF (QUADRU) WRITE (LUPRI,'(10X,A)')
     &                         '- quadrupole moment integrals'
         IF (THETA)  WRITE (LUPRI,'(10X,A)')
     &                         '- traceless quadrupole moment integrals'
         IF (SECMOM) WRITE (LUPRI,'(10X,A)')'- second moments integrals'
         IF (SPNORB) WRITE (LUPRI,'(10X,A)')
     &                                  '- spatial spin-orbit integrals'
         IF (OCTUPO) WRITE (LUPRI,'(10X,A)')
     &                                    '- octupole moment integrals'
         IF (CARMOM) THEN
            IF (IORCAR .GT. 0) THEN
               WRITE (LUPRI,'(10X,A,I2,A)')
     &          '- Cartesian multipole moment integrals of orders',
     &          ABS(IORCAR),' and lower'
            ELSE
               WRITE (LUPRI,'(10X,A,I2)')
     &          '- Cartesian multipole moment integrals of order',
     &           ABS(IORCAR)
            END IF
         END IF
         IF (SPHMOM) THEN
            IF (IORSPH .GT. 0) THEN
               WRITE (LUPRI,'(10X,A,I2,A)')
     &          '- Spherical multipole moment integrals of orders',
     &          ABS(IORSPH),' and lower'
            ELSE
               WRITE (LUPRI,'(10X,A,I2)')
     &          '- Spherical multipole moment integrals of order',
     &           ABS(IORSPH)
            END IF
         END IF
         IF (KINENE) WRITE (LUPRI,'(10X,A)')
     &      '- electronic kinetic energy'
         IF (MASSVL) WRITE (LUPRI,'(10X,A)')
     &       '- mass velocity integrals'
         IF (DARWIN) WRITE (LUPRI,'(10X,A)')
     &       '- 1-electron Darwin integrals'
         IF (NMRISS) THEN
            IF (FERMI) THEN
               WRITE (LUPRI,'(10X,A)')'- Fermi contact integrals'
               WRITE (LUPRI,'(10X,A)')
     &              '  (Dirac delta function integrals)'
            END IF
            IF (PSO) THEN
               WRITE (LUPRI,'(10X,A)')
     &              '- paramagnetic spin-orbit integrals'
               WRITE (LUPRI,'(10X,A)')
     &              '  (nuclear moment - electron orbit coupling)'
            END IF
            IF (SPIDIP) THEN
               WRITE (LUPRI,'(10X,A)')'- spin-dipole integrals'
               WRITE (LUPRI,'(10X,A)')
     &              '  (electron spin - nuclear moment coupling)'
            END IF
            IF (DSO) THEN
               WRITE (LUPRI,'(10X,A)')
     &              '- diamagnetic spin-orbit integrals'
               WRITE (LUPRI,'(10X,A)')
     &              '  (indirect nuclear dipole - dipole coupling)'
            END IF
            IF (SDFC) THEN
               WRITE (LUPRI,'(10X,A)')
     &             '- spin-dipole + Fermi contact integrals'
               WRITE (LUPRI,'(10X,A)')
     &             '  (electron spin - nuclear magnetic field coupling)'
            END IF
         END IF
         IF (HDO) WRITE (LUPRI,'(10X,A)')
     &                            '- half-derivative overlap integrals'
         IF (S1MAG) WRITE (LUPRI,'(10X,A)')
     &       '- first magnetic derivatives of overlap integrals'
         IF (S1MAGT) WRITE (LUPRI,'(10X,A)')
     &       '- test of first magnetic derivative of overlap integrals'
         IF (S2MAG) WRITE (LUPRI,'(10X,A)')
     &       '- second magnetic derivatives of overlap integrals'
         IF (S2MAGT) WRITE (LUPRI,'(10X,A)')
     &      '- test of second magnetic derivatives of overlap integrals'
         IF (ANGMOM) WRITE (LUPRI,'(10X,A)')
     &      '- electronic angular momentum around the origin'
         IF (ANGLON) WRITE (LUPRI,'(10X,A)')
     &      '- electronic angular momentum around the nuclei'
         IF (LONMOM) WRITE (LUPRI,'(10X,A)')
     &      '- London orbital contribution to angular momentum'
         IF (MAGMOM) WRITE (LUPRI,'(10X,A)')
     &      '- one-electron contribution to magnetic moment'
         IF (MGMOMT) WRITE (LUPRI,'(10X,A)')
     &      '- test of London contribution to angular momentum'
         IF (DSUSNL) WRITE (LUPRI,'(10X,A)')
     &   '- Magnetic susceptibility without London orbital contribution'
         IF (DSUSLL) WRITE (LUPRI,'(10X,A)')
     &      '- Angular London orbital contribution to magnetic susc.'
         IF (DSUSLH) WRITE (LUPRI,'(10X,A)')
     &      '- London orbital contribution to magnetic susceptibility'
         IF (DIASUS) WRITE (LUPRI,'(10X,A)')
     &      '- Magnetic susceptibility integrals'
         IF (DSUTST) WRITE (LUPRI,'(10X,A)')
     &     '- Test of London orbital contr. to magnetic susc. integrals'
         IF (NUCSNL) WRITE (LUPRI,'(10X,A)')
     &     '- Nuclear shieldings without London orbital contribution'
         IF (NUCSLO) WRITE (LUPRI,'(10X,A)')
     &     '- London orbital contribution to nuclear shieldings'
         IF (NUCSHI) WRITE (LUPRI,'(10X,A)')
     &     '- Nuclear shielding tensor integrals'
         IF (NSNLTS) WRITE (LUPRI,'(10X,A)')
     &     '- Test of London orbital contribution to nuclear shieldings'
         IF (NELFLD) WRITE (LUPRI,'(10X,A)')
     &     '- Electric field at the nucleus'
         IF (NSNLTS) WRITE(LUPRI,'(10X,A)')
     &     '- Test of non-London orbital contr. to nuclear shieldings'
         IF (NSTTST) WRITE (LUPRI,'(10X,A)')
     &     '- Test of nuclear shielding tensor integrals'
         IF (EFGCAR) WRITE (LUPRI,'(10X,A)')
     &            '- Cartesien electric field gradient integrals'
         IF (EFGSPH) WRITE (LUPRI,'(10X,A)')
     &            '- Spherical electric field gradient integrals'
         IF (S1MAGL) WRITE (LUPRI,'(10X,A)')
     &        '- Bra-differentiated overlap matrix with respect to B'
         IF (S1MAGR) WRITE (LUPRI,'(10X,A)')
     &        '- Ket-differentiated overlap matrix with respect to B'
         IF (HBDO) WRITE (LUPRI,'(10X,A)')
     &        '-Half B-differentiated overlap matrix'
         IF (HDOBR) WRITE (LUPRI,'(10X,A)')
     &        '- Ket-differentiated hdo-integrals with respect to B'
         IF (S1MLT) WRITE (LUPRI,'(10X,A)')
     &        '- Test of bra-diff. overlap matrix with respect to B'
         IF (S1MRT) WRITE (LUPRI,'(10X,A)')
     &        '- Test of ket-diff. overlap matrix with respect to B'
         IF (HDOBRT) WRITE (LUPRI,'(10X,A)')
     &        '- Test og ket-diff. hdo-integrals with respect to B'
         IF (SQHDOL) WRITE (LUPRI,'(10X,A)')
     &        '- Bra differentiated half-derivative overlap matrix'
         IF (SQHDOR) WRITE (LUPRI,'(10X,A)')
     &        '- Ket differentiated half-derivative overlap matrix'
         IF (NUCPOT) WRITE (LUPRI,'(10X,A)')
     &    '- Potenial energy of interaction of electrons with a nucleus'
         IF (NPOTST) WRITE (LUPRI,'(10X,A)')
     &            '- Test of nuclear potential energy'
         IF (MGMO2T) WRITE (LUPRI,'(10X,A)')
     &            '- Test of two-electron part of magnetic moment'
         IF (SUSCGO) WRITE (LUPRI,'(10X,A)')
     &      '- Diamagnetic magnetizability using common gauge origin'
         IF (NSTCGO) WRITE (LUPRI,'(10X,A)')
     &      '- Diamagnetic shielding tensor using common gauge origin'
         IF (GFACDI) WRITE (LUPRI,'(10X,A)')
     &      '- Diamagnetic London contribution to rotational g-factors'
         IF (EXPIKR) WRITE (LUPRI,'(10X,A)')
     &      '- Cosine and sine integals'
         IF (CM1) THEN
            WRITE (LUPRI,'(10X,A)')
     &         '- First order magnetic derivative of electric field'
            WRITE (LUPRI,'(12X,A,A1,A)')
     &         'Electric field applied in ',FIELD1(1:1),'-direction'
         END IF
         IF (CM2) THEN
            WRITE (LUPRI,'(10X,A)')
     &         '- Second order magnetic derivative of electric field'
            WRITE (LUPRI,'(12X,A,A1,A)')
     &         'Electric field applied in ',FIELD2(1:1),'-direction'
         END IF
         IF (PROPRI) WRITE (LUPRI,'(/A)')
     &      ' All one-electron property integrals are printed.'
         IF (S1ELE) WRITE (LUPRI,'(10X,A)')
     &       '- first electric derivatives of overlap integrals,'//
     &        'Type A'
         IF (S1ELB) WRITE (LUPRI,'(10X,A)')
     &       '- first electric derivatives of overlap integrals,'//
     &        'Type B'
         IF (ONEELD) THEN
            WRITE (LUPRI,'(10X,A)')
     &       '- first electric derivatives of one-electron'
            WRITE (LUPRI,'(12X,A)') 'Hamiltonian integrals'
         ENDIF
      END IF

      IF (EXPIKR) WRITE (LUPRI,'(/,A,3F20.15)')
     &    ' Wave numbers:', (EXPKR(I),I=1,3)
      IF (SOTEST) WRITE (LUPRI,'(/,A,/)')
     *    ' Test of spatial spin-orbit integrals.'
      IF (.NOT.HAMILT) WRITE (LUPRI,'(/,A,/)')
     *    ' Ordinary (field-free non-relativistic) Hamiltonian '/
     *     /'integrals not calculated.'
      IF (SUPMAT) WRITE (LUPRI,'(/,A,/,A,/)') ' In this run '//
     *     ' precalculated two-electron integrals are transformed',
     *     ' to P-supermatrix elements.'
      IF (MGMO2T) THEN
         WRITE (LUPRI,'(/A,D12.6)')
     &     ' Threshold for testing two-electron integrals:',PRTHRS
      END IF
      IF (NMRISS) THEN
         IF (ALLATM) THEN
            WRITE (LUPRI,'(/2A)')
     &         ' Integrals for all indirect spin-spin',
     &         ' coupling tensors are calculated.'
         ELSE
            WRITE (LUPRI,'(/2A,/)')
     &         ' Indirect spin-spin integrals involving the following',
     &         ' nuclei are calculated:'
            WRITE (LUPRI,'(10X,20I3)') (IPATOM(I),I = 1, NPATOM)
         END IF
         IF (DSO) THEN
            WRITE (LUPRI,'(/2A,I3)')
     &        ' Number of integration points for diamagnetic',
     &        ' spin-orbit integrals: ',NPQUAD
            IF (.NOT.TRIANG) WRITE (LUPRI,'(A)')
     &        ' Integrals for symmetry related coupling tensors'
     &           //' JAB and JBA calculated.'
         END IF
      END IF
      IF (.NOT. ALLATM) THEN
            WRITE (LUPRI,'(/A/2A,/)')
     &         ' For property integrals involving atomic centers',
     &         '     only integrals involving the following',
     &         ' nuclei are calculated:'
            WRITE (LUPRI,'(10X,20I3)') (IPATOM(I),I = 1, NPATOM)
      END IF
C
C     **** Process input for various program sections  *****
C
  200 CONTINUE
      CALL REAINI(IPRDEF,RELCAL,TSTINP)
 201  PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         IF (WORD(1:2) .EQ. '**') GO TO 1
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4,5,6), I
            END IF
  210    CONTINUE
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in '//WORD1)
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,
     *      '" illegal or out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in '//WORD1//', error in prompt.')
      END IF
    2   CALL REAINP(WORD,RELCAL,TSTINP)
      GO TO 201
    3   CALL HR1INP(WORD)
      GO TO 201
    4   CALL HR2INP(WORD)
      GO TO 201
#ifndef PRG_DIRAC
    5   CALL HRSINP(WORD)
#else
    5   CALL QUIT('*SUPINT not implemented in Dirac')
#endif /* ifndef PRG_DIRAC */
      GO TO 201
    6 CONTINUE
C
    1 CONTINUE
C
      CALL HR1INP(WORD)
      CALL HR2INP(WORD)
#ifndef PRG_DIRAC
      CALL HRSINP(WORD)
#endif /* ifndef PRG_DIRAC */
 1000 CONTINUE
C     We should check if LUONEL is opened and closed properly,
C     kr and vb, April-96
C
C      CLOSE (LUONEL)
C
C     Read geometry, orbital spec., etc. (MOLECULE or XYZ format).
C     Parameter .TRUE. in call indicates that LUONEL for Dalton must be written.
C
      CALL READIN(.FALSE.)
      CALL SETDCH
      IF (DORLM .AND. .NOT. CAVUSR) THEN
         KGEOM = 1
         KMASS = KGEOM + 3*NATOMS
         KNAT  = KMASS + NATOMS
         KNUMIS= KNAT  + (NATOMS + 1)/IRAT
         KLAST = KNUMIS+ (NATOMS + 1)/IRAT
         IF (KLAST .GT. LWORK) CALL STOPIT('HERDRV','CMMASS',KLAST,
     &                                     LWORK)
         CALL CMMASS(WORK(KGEOM),WORK(KMASS),WORK(KNAT),
     &      WORK(KNUMIS),IPRDEF)
         CAVORG(1:3) = CMXYZ(1:3)
         WRITE(LUPRI,'(/A,3F15.6)')
     &      ' Cavity center (at center of mass):',
     &      CAVORG(1:3)
      END IF
C
      CLOSE(UNIT=LUCMD,STATUS='KEEP')
C
      RETURN
      END

C  /* Deck hr1inp */
      SUBROUTINE HR1INP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include "orgcom.h"
#include "cbiher.h"
#include "cbihr1.h"
C
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.SOLVEN', '.ALLRLM', '.CAVORG',
     &            'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      SET = .TRUE.
C
C     Initialize /CBIHR1/
C
      RUNONE = .TRUE.
      IPRONE = IPRDEF
      DORLM  = .FALSE.
      ALLRLM = .FALSE.
      CAVUSR = .FALSE.
      CAVORG(1:3) = 0.0d0
C
      NEWDEF = WORD .EQ. '*ONEINT'
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
                     GO TO (1,2,3,4,5,6,7,8,9,10), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ONEINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ONEINP.')
    1          CONTINUE
                  RUNONE = .FALSE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRONE
                  IF (IPRONE .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  DORLM = .TRUE.
                  READ (LUCMD, '(I5)') LMAX
               GO TO 100
    4          CONTINUE
                  ALLRLM = .TRUE.
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) (CAVORG(I),I = 1, 3)
                  CAVUSR = .TRUE.
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
    9          CONTINUE
               GO TO 100
   10          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ONEINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ONEINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .EQ. 0) RETURN
      IF (NEWDEF) THEN
         CALL HEADER('Changes of defaults for ONEINP:',1)
         IF (.NOT.RUNONE) THEN
            WRITE (LUPRI,'(A)') ' No one-electron integrals calculated.'
         ELSE
            IF (IPRONE .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     &         ' Print level in ONEINT:',IPRONE
         END IF
         IF (DORLM) THEN
            WRITE (LUPRI,'(A/A,I2)')
     &         ' One-electron RLM integrals calculated.',
     &         ' Maximum L quantum number: ', LMAX
            IF (ALLRLM) THEN
               WRITE (LUPRI,'(A)') ' All symmetries saved on file.'
            ELSE
               WRITE (LUPRI,'(A)')
     &            ' Only totally symmetric integrals saved on file.'
            END IF
            IF (CAVUSR) WRITE(LUPRI,'(A,3F15.10)')
     &         ' User supplied cavity center',(CAVORG(I),I=1,3)
         END IF
      END IF
      RETURN
      END
C  /* Deck hr2inp */
      SUBROUTINE HR2INP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D1 = 1.0D0,D0 = 0.0D0)
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include "cbiher.h"
#include "cbihr2.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
! infpar.h : NUMNOD
#include "infpar.h"
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.PANAS ', '.RETURN', '.SOFOCK',
     &            '.TIME  ', '.ICEDIF', '.SCREEN', '.THRFAC', '.AOFOCK'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      SET = .TRUE.
C
C     Initialize /CBIHR2/
C
      RUNTWO = .NOT.NOTWO
      IPRTWO = IPRDEF
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      IPRSUM = 0
      RTNTWO = .FALSE.
      TKTIME = .FALSE.
      IF (RELCAL .AND. NUMNOD .LE. 24) THEN
         SOFOCK = .TRUE.
      ELSE
         SOFOCK = .FALSE.
      ENDIF
      ICDIFF = 1
      IEDIFF = 1
      USRSCR = .FALSE.
      SCRFCK = 1.0D-12
      THRFAC(1) = D1
      THRFAC(2) = D1
C
      NEWDEF = WORD .EQ. '*TWOINT'
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
                     GO TO (1,2,3,4,5,6,7,8,9,10), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in HR2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in HR2INP.')
    1          CONTINUE
                  RUNTWO = .FALSE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(5I5)') IPRTWO,
     &                     IPRNTA, IPRNTB, IPRNTC, IPRNTD
                  IPRSUM = IPRNTA + IPRNTB + IPRNTC + IPRNTD
                  IF (IPRTWO .EQ. IPRDEF .AND. IPRSUM .EQ. 0) THEN
                     ICHANG = ICHANG - 1
                  END IF
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*,ERR=35) PANAS
                  GOTO 36
 35               PANAS = 0.25D0
 36               CONTINUE 
C     
C     We cannot use new integral code for Panas correction
C
                  SEGBAS = .FALSE.
               GO TO 100
    4          CONTINUE
                  RTNTWO = .TRUE.
               GO TO 100
    5          CONTINUE
C&&&& SOFOCK - construction of Fock matrices in SO-basis
                  SOFOCK = .TRUE.
               GO TO 100
    6          CONTINUE
                  TKTIME = .TRUE.
               GO TO 100
    7          CONTINUE
C&&&& ICEDIF Separate screening of Coulomb and exchange contributions
C&&&& in direct SCF
                  READ (LUCMD,*) ICDIFF,IEDIFF
               GO TO 100
    8          CONTINUE
C&&&& SCREEN: Screening threshold in direct construction of Fock matrices
                  READ (LUCMD,*) SCRFCK
                  USRSCR = .TRUE.
               GO TO 100
    9          CONTINUE
C&&& THRFAC: Factors to multiply LL-integral threshold for SL- and SS - integrals
C&&& This option only used in DIRAC
                 READ(LUCMD,*) THRFAC(1),THRFAC(2)
               GO TO 100
   10          CONTINUE
C&&& AOFOCK: Direct construction of Fock-matrix in AO-basis
                 SOFOCK = .FALSE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in HR2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in HR2INP.')
            END IF
      END IF
  300 CONTINUE
      ICEDIF = ICDIFF + 2*IEDIFF
CTROND      IF (ICHANG .EQ. 0) RETURN
      IF (NEWDEF) THEN
         CALL HEADER('Set-up from HR2INP:',1)
         IF (.NOT.(RUNTWO.OR.DIRCAL)) THEN
            WRITE (LUPRI,'(A)') ' No two-electron integrals calculated.'
         ELSE
            WRITE (LUPRI,'(A,I5)') ' Print level in TWOINT:',IPRTWO
            IF (IPRSUM .GT. 0) THEN
                 WRITE (LUPRI,'(A,4I3)')
     &                ' Extra output for the following shells:',
     &                 IPRNTA, IPRNTB, IPRNTC, IPRNTD
                IF (RTNTWO) WRITE (LUPRI,'(A)')
     &               ' Program will exit TWOINT after these shells.'
            END IF
            IF (TKTIME) WRITE (LUPRI,'(/,2A)') ' Detailed timing for',
     &         ' integral calculation will be provided.'
            IF (PANAS .NE. 0.0D0) WRITE (LUPRI,'(/,A,F10.5)')
     &           ' Coulomb integrals screened with a factor of',PANAS
            IF (DIRCAL) THEN
              IF(SOFOCK) THEN
                WRITE(LUPRI,'(1X,A)')
     &            '* Direct calculation of Fock matrices in SO-basis.' 
              ELSE 
                WRITE(LUPRI,'(1X,A)')
     &            '* Direct calculation of Fock matrices in AO-basis.' 
              ENDIF
              IF(SCRFCK.GT.D0) THEN
                IF (USRSCR) THEN
                  WRITE(LUPRI,'(1X,A,A,1P,E8.2)') 
     &            '* User specified screening threshold in direct Fock',
     &            ' matrix construction: ',SCRFCK
                ELSE
                  WRITE(LUPRI,'(1X,A,A,1P,E8.2)') 
     &            '* Default screening threshold in direct Fock',
     &            ' matrix construction: ',SCRFCK
                END IF
                IF(ICDIFF.EQ.1) WRITE(LUPRI,'(1X,A)')
     &    '* Separate density screening of Coulomb integral batches'
                IF(IEDIFF.EQ.1) WRITE(LUPRI,'(1X,A)')
     &    '* Separate density screening of exchange integral batches'
                CALL FLSHFO(LUPRI)
              ELSE
                 WRITE(LUPRI,'(4X,A)') 
     &           '---> WARNING : Integral screening turned off !'
                 CALL FLSHFO(LUPRI)
              ENDIF
            ENDIF
            IF(RELCAL .AND. (THRFAC(1).NE.D1.OR.THRFAC(2).NE.D1)) THEN
              WRITE(LUPRI,'(1X,A,2(/3X,A,1P,D9.3))')
     +           '* Threshold factors for omitting integrals:',
     +           'SL-integrals: ',THRFAC(1),
     +           'SS-integrals: ',THRFAC(2)
            ENDIF
         END IF
      END IF
      RETURN
      END
#ifndef PRG_DIRAC
C  /* Deck hrsinp */
      SUBROUTINE HRSINP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include "cbiher.h"
#include "cbihrs.h"
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.NOSYMM', '.OLDSUP', '.THRESH',
     &            'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      SET = .TRUE.
C
C     Initialize /CBIHRS/
C
      RUNSUP = SUPMAT
      IPRSUP = IPRDEF
      NOSSUP = .FALSE.
      OLDSUP = .FALSE.
      THRSUP = -1.0D0
C
      NEWDEF = WORD .EQ. '*SUPINT'
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
                     GO TO (1,2,3,4,5,6,7,8,9,10), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in SUPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in SUPINP.')
    1          CONTINUE
                  RUNSUP = .FALSE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRSUP
                  IF (IPRSUP .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  NOSSUP = .TRUE.
               GO TO 100
    4          CONTINUE
                  OLDSUP = .TRUE.
               GO TO 100
    5          CONTINUE
                  READ (LUCMD, '(D12.6)') THRSUP
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
    9          CONTINUE
               GO TO 100
   10          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in SUPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in SUPINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .EQ. 0) RETURN
      IF (RUNSUP. AND. NEWDEF) THEN
         CALL HEADER('Changes of defaults for SUPINT:',1)
         IF (IPRSUP .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     &         ' Print level in SUPINT:',IPRSUP
         IF (THRSUP .NE. -1.0D0) WRITE (LUPRI,'(A,D12.2)')
     &         ' Threshold for supermatrix integrals:', THRSUP
         IF (OLDSUP) THEN
            WRITE (LUPRI,'(A)')
     &         ' Old format for supmatrix integral file.'
            IF (NOSSUP) THEN
               NOSSUP = .FALSE.
               WRITE (LUPRI,'(/A/)' )
     &           ' Note: NOSSUP has been set .FALSE. to conform with'
     &           //' .OLDSUP in input.'
            END IF
         END IF
         IF (NOSSUP) THEN
            WRITE (LUPRI,'(A)') ' No symmetry used in SUPINT.'
         ELSE
            WRITE (LUPRI,'(A)') ' Symmetry used in SUPINT.'
         END IF
      END IF
      RETURN
      END
#endif /* ifndef PRG_DIRAC */
C  /* Deck setdch */
      SUBROUTINE SETDCH
C
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "dorps.h"
C
      DO 100 IREP = 0, MAXREP
         DOREPS(IREP) = .TRUE.
  100 CONTINUE
      RETURN
      END
