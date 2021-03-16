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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mp2inp */
      SUBROUTINE MP2INP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for MP2-module
C
C     Written by Jon K. Laerdahl and T. Saue (Dec 1996)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
      PARAMETER(D0 = 0.0D0)
      PARAMETER (NTABLE = 20)
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*2 WS
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbmp2.h"
#include "dcbmpt.h"
#include "dcbham.h"
      DIMENSION WORK(LWORK),IR(2,2)
C
      SAVE SET
      DATA TABLE /'.PRINT ','.OCCUP ','.VIRTUA','.INTFLG','xMSOUT ',
     &            '.POSITR','.SCHEME','xTHROUT','.SPECIA','.VIRTHR',
     &            '.IJTSK ','.SCLMEM','.SCREEN','.ORGALG','.NOSORT',
     &            '.XXXXXX','.COMMIN','.TPRI34','.TPRI44','.xxxxxx'/
      DATA SET/.FALSE./
#include "ibtfun.h"
#include "memint.h"
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
C     Initialize /DCBTRA/
C     ===================
C
      ISTRMP2 = 2
      IPRMP2 = 0
      ILLINT = IBTAND(INTGEN,1)
      ILSINT = IBTAND(INTGEN/2,1)
      ISSINT = IBTAND(INTGEN/4,1)
      IGTINT = IBTAND(INTGEN/8,1)
      THROUT = 1.0D-14
C     Default orbital strings is all occupied and all virtuals
      DO I = 1,NFSYM
        MP2_INDSTR(1,I) = ' '
        MP2_INDSTR(2,I) = ' '
        IF(NOCC(I).GT.0) THEN
          WRITE(MP2_INDSTR(1,I),'(I4,A2,I4)') 1,'..',NOCC(I)
        ENDIF
        WRITE(MP2_INDSTR(2,I),'(A)') 'all'
      ENDDO
      MP2_INTFLG = ILLINT + 2*ILSINT + 4*ISSINT + 8*IGTINT
      IF (LEVYLE.OR.BSS.or.x2c) MP2_INTFLG = ILLINT
C     MP2_ANTIS  = .TRUE.
      MP2_ANTIS  = .FALSE.
C     note TODO: ANTIS is not implemented!
      MP2_MSOUT  = .FALSE.
C     note TODO: MSOUT is not implemented!
      MP2NOP = .TRUE.
      TESTLS = .FALSE.
      TESTSL = .FALSE.
      DMP2_VIRTHR = DUMMY
      IJTSK = 0
      MAXSCL = 0
Ctsaue      SCRMP2 = 1.D-14
Ctsaue: turning of default screening for now since further testing is needed
      SCRMP2 = -1.0D0
      SRTSHL = .TRUE.
      CMMNFO = .FALSE.
      TPRI34 = -1.D0
      TPRI44 = -1.D0
C
C     Activate/deactivate original algorithm
C
      MP2ORG = .FALSE.
C
C     Process input for MP2CAL
C     ========================
C
      NEWDEF = (WORD .EQ. '*MP2CAL')
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
     &                    11,12,13,14,15,16,17,18,19,20), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in MP2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in MP2INP.')
    1          CONTINUE
C&&&& PRINT - print level in MP2 module
                  READ(LUCMD,*) IPRMP2
               GO TO 100
    2          CONTINUE
C&&&& OCCUP - String of active occupied pairs in MP2 module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') MP2_INDSTR(1,I)
                  ENDDO
               GO TO 100
    3          CONTINUE
C&&&& VIRTUA -String of active virtual pairs in MP2 module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') MP2_INDSTR(2,I)
                  ENDDO
               GO TO 100
    4          CONTINUE
C&&&& INTFLG -Integral classes not included in 4-index transformation
                  IF (IGTINT .EQ. 1) THEN
                     READ(LUCMD,*,IOSTAT=IOS) 
     &                  ILLINT,ILSINT,ISSINT,IGTINT
                     IF (IOS.NE.0) THEN
                       CALL QUIT(
     &       'MP2INP: Error in reading ILLINT,ILSINT,ISSINT,IGTINT')
                     ENDIF
                  ELSE
                     READ(LUCMD,*,IOSTAT=IOS) ILLINT,ILSINT,ISSINT
                     IF (IOS.NE.0) THEN
                       CALL QUIT(
     &       'MP2INP: Error in reading ILLINT,ILSINT,ISSINT')
                     ENDIF
                  END IF
                  MP2_INTFLG = ILLINT + 2*ILSINT + 4*ISSINT + 8*IGTINT
               GO TO 100
    5          CONTINUE
C&&&& MSOUT :Write 4-index transformed integrals to file
C     aug 2001 hjaaj TODO: disabled this option because only false implemented
                  MP2_MSOUT = .TRUE.
               GO TO 100
    6          CONTINUE
C&&&& POSITR: Allow positronic orbitals
                  MP2NOP = .FALSE.
               GO TO 100
    7          CONTINUE
C&&&& SCHEME - Choose transformation scheme
                  READ(LUCMD,*) ISTRMP2
               GO TO 100
    8          CONTINUE
C&&&& THROUT - Treshold for output integrals
C     aug 2001 hjaaj TODO: disabled this option because not implemented
                  READ(LUCMD,*,IOSTAT=IOS) THROUT
                  IF (IOS.NE.0) THEN
                       CALL QUIT(
     &       'MP2INP: Error in reading THROUT !')
                  ENDIF
               GO TO 100
    9          CONTINUE
C&&&& SPECIAL: Mostly for testing
C
C   'LS' gives the (LL|SS) contribution to the MP2 energy
C   'SL' gives the (SS|LL) contribution to the MP2 energy
C
                  IF (LEVYLE.OR.BSS.or.x2c) CALL QUIT(
     &         '.SPECIAL testing in MP2 not allowed for Levy-Leblond')
                  READ(LUCMD,'(A2)') WS
                  IF (WS.EQ.'LS') THEN
                     TESTLS = .TRUE.
                  ELSEIF (WS.EQ.'SL') THEN
                     TESTSL = .TRUE.
                  ELSE
                     CALL QUIT('Unrecognized .SPECIAL option for MP2')
                  ENDIF
               GO TO 100
   10          CONTINUE
C&&&& VIRTHR: Only include Kramers pairs with energy below
C     threshold in Virtual space.
               READ(LUCMD,'(F20.5)',IOSTAT=IOS) DMP2_VIRTHR
               IF (IOS.NE.0) THEN
                 CALL QUIT(
     &         'MP2INP: Error in reading DMP2_VIRTHR !')
               ENDIF
               GO TO 100
   11          CONTINUE
C&&&& IJTSK: Max number of active I (and J) Kramers
C     pairs in each task.
               READ(LUCMD,*) IJTSK
               GO TO 100
   12          CONTINUE
C&&&& SCLMEM: Max internal memory available for scalar integrals
               READ(LUCMD,*) MAXSCL
               GO TO 100
   13          CONTINUE
C&&&& SCREEN: Screening threshold for 4-index transformation
               READ(LUCMD,*) SCRMP2
               GO TO 100
   14          CONTINUE
C&&&& ORGALG: Activate/deactivate original MP2 algorithm
               MP2ORG = .TRUE.
               GO TO 100
   15          CONTINUE
C&&&& NOSORT: Controls the sorting of shells in decreasing size
               SRTSHL = .FALSE.
               GO TO 100
   16          CONTINUE
               GO TO 100
   17          CONTINUE
C&&&& COMMINfo: Active printing of communication status (for parallel jobs)
               CMMNFO = .TRUE.
               GO TO 100
   18          CONTINUE
C&&&& TPRI34: Threshold for printing 3/4 transformed inetgrals
               READ(LUCMD,*) ITHR
               TPRI34 = 10.D0**(ITHR)
               GO TO 100
   19          CONTINUE
C&&&& TPRI44: Threshold for printing fully-transformed inetgrals
               READ(LUCMD,*) ITHR
               TPRI44 = 10.D0**(ITHR)
               GO TO 100
   20          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in MP2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in MP2INP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT.DOMP2) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   'MP2INP: Set-up for RMP2 calculation:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I3)') '* General print level   : ',IPRMP2
      IF(TESTLS) THEN
         ILLINT = 0
         ISSINT = 0
         ILSINT = 1
         IGTINT = 0
         MP2_INTFLG = ILLINT + 2*ILSINT + 4*ISSINT + 8*IGTINT
         WRITE(LUPRI,'(1X,A)')
     +     '* Testing (LL|SS) ONLY in MP2.'
      ENDIF
      IF(TESTSL) THEN
         ILLINT = 0
         ISSINT = 0
         ILSINT = 1
         IGTINT = 0
         MP2_INTFLG = ILLINT + 2*ILSINT + 4*ISSINT + 8*IGTINT
         WRITE(LUPRI,'(1X,A)')
     +     '* Testing (SS|LL) ONLY in MP2.'
      ENDIF
C
      IF (LEVYLE.OR.BSS.or.x2c) MP2_INTFLG = ILLINT
      IF(SCRMP2.GT.D0) THEN
        WRITE(LUPRI,'(1X,A,1P,E8.2)') '* Screening threshold : ',SCRMP2
      ELSE
        WRITE(LUPRI,'(1X,A)') '* No screening.'
      ENDIF
C
      IF (MP2ORG) THEN
         WRITE(LUPRI,'(1X,A)') '* Original MP2 algorithm.'
         IF ((.NOT. SRTSHL) .OR. CMMNFO .OR. (TPRI34 .GT. 0.D0)
     &        .OR. (TPRI44 .GT. 0.D0))
     &        WRITE(LUPRI,'(1X,A/A)')
     &    '  - Keywords .NOSORT, .COMMIN, .TPRI34 and .TPRI44',
     &    '    have no effect when running original MP2 algorithm.'
      ELSE
         WRITE(LUPRI,'(1X,A)')
     &        '* New MP2 algorithm (better parallelization).'
         IF (SRTSHL) THEN
            WRITE(LUPRI,'(1X,A)') '  - Shells will be sorted.'
         ELSE
            WRITE(LUPRI,'(1X,A)') '  - No sorting of shells.'
         END IF
         IF (CMMNFO) WRITE(LUPRI,'(1X,A)')
     &           '  - Detailed (MPI) communication info.'
         IF (TPRI34 .GT. 0.D0) WRITE(LUPRI,'(1X,A,G10.3)')
     &        '  - Print 3/4 transf. integrals larger than:   ', TPRI34
         IF (TPRI44 .GT. 0.D0) WRITE(LUPRI,'(1X,A,G10.3)')
     &        '  - Print fully-transf. integrals larger than: ', TPRI44
      END IF
C
      IF (ILLINT.NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* LL Integrals Not included.'
      IF (ILSINT.NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* LS Integrals Not included.'
      IF (ISSINT.NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* SS Integrals Not included.'
      IF (IJTSK.GT.0) WRITE(LUPRI,'(1X,A,I4)')
     +     '* Max number of Kramers pairs I (J) in task: ', IJTSK
      IF (MP2_MSOUT) THEN
         WRITE(LUPRI,'(1X,A)')
     +        '* 4-index transformed integrals written to file.'
         CALL QUIT('MSOUT not implemented in MP2 program')
      ENDIF
C     Active occupied and virtual space: only electrons allowed
      WRITE(LUPRI,'(1X,A)') '* Active spaces:'
      IF(MP2NOP) THEN
        WRITE(LUPRI,'(1X,A)') '* Electronic orbitals only.'
        IC   = 1
        DO I = 1,NFSYM
          IR(1,I) = 1
          IR(2,I) = NFBAS(I,1)
        ENDDO
      ELSE
        IC   = 0
        DO I = 1,NFSYM
          IR(1,I) = -NFBAS(I,2)
          IR(2,I) =  NFBAS(I,1)
        ENDDO
      ENDIF
      DO I = 1,NFSYM
        WRITE(LUPRI,'(4X,A,A3)') 'Fermion ircop:',FREP(I)
        NVEC = 0
        CALL  NUMLST(MP2_INDSTR(1,I),IDUMMY,NFBAS(I,IC),
     &               IR(1,I),IR(2,I),I,NVEC)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A)')
     &       '- Occupied space : No electrons'
        ELSE
          WRITE(LUPRI,'(4X,A,A72)')
     &       '- Occupied space : ', MP2_INDSTR(1,I)
        ENDIF
        IF(DMP2_VIRTHR.NE.DUMMY) MP2_INDSTR(2,I) = 'all '
        NVEC = 0
        CALL  NUMLST(MP2_INDSTR(2,I),IDUMMY,NFBAS(I,IC),
     &               IR(1,I),IR(2,I),I,NVEC)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A)')
     &       '- Virtual  space : No electrons'
        ELSE
          WRITE(LUPRI,'(4X,A,A72)')
     &       '- Virtual  space : ', MP2_INDSTR(2,I)
        ENDIF
      ENDDO
C
      IF(DMP2_VIRTHR.NE.DUMMY) WRITE(LUPRI,'(1X,A,F7.2,A)')
     &     '* Only virtuals with energy below ', DMP2_VIRTHR,
     &     ' a.u. included in active space.'
      IF(MAXSCL.NE.0) WRITE(LUPRI,'(2(1X,A),I10,F10.3,A)')
     &     '* Max memory available for scalar two electron',
     &     'integrals:',MAXSCL,MAXSCL*8./(1024.*1024.),' Mb'
C
      IF (ISTRMP2.NE.2)
     &     CALL QUIT('Can currently only do MP2 with .SCHEME 2')

C
 999  CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck trainp */
      SUBROUTINE TRAINP(WORK,LWORK)
C*****************************************************************************
C
C     Input section for module for 4-index transformation
C
C     Written by Jon K. Laerdahl, T. Saue and L. Visscher (Dec 1996)
C
C     Called from PAMINP/dirrdn.F
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
      PARAMETER(D0=0.0D0)
      PARAMETER (NDIR = 3, NTABLE = 29)
C
      LOGICAL SET, NEWDEF,LSPES,LSPES2,LSPES4
      CHARACTER PROMPT*1, WORD*7, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7
      CHARACTER*7 LINE
      DIMENSION WORK(LWORK),IR(2,2),NTVEC(4)
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbham.h"

C
      SAVE SET
      DATA TABDIR /'*END OF','*PRPTRA','*HSCC  '/
C
      DATA TABLE /'.PRINT ','.ACTIVE','.INTFLG','xANTIS ','.4INDEX',
     &            '.EXACOR','.POSITR','.CORE  ','.SCHEME','.THROUT',
     &            '.NO2IND','.NO4IND','.2INDEX','.INTFL2','.INTFL4',
     &            '.PRPTRA','.SCREEN','.RCORBS','.NOMDCI','.NOSCAT',
     &            '.PAR4BS','.SCATTE','.MDCINT','.HSCC  ','.HTSORT',
     &            '.ASCII ','.MOFILE','.CORE2 ','.CHOLES'/
      DATA SET/.FALSE./
C
#include "ibtfun.h"
#include "memint.h"
      IF (SET) RETURN
!Miro: nasty workaround to enable repeated input reading in the RelCC geometry optimization
      !SET = .TRUE. 
C
C     Local variables
C
      LSPES  = .FALSE.
      LSPES2 = .FALSE.
      LSPES4 = .FALSE.

!       ... no HSCC
      DOHSCC = .FALSE.

C
C     Initialize /DCBTRA/
C     ===================
C
      IPRTRA  = 0
Caspg, default in MOLTRA is scheme 6 unless it's a GOSCIP, DIRRCI or LUCITA run... 
C      if problems are encountered, should be reverted to scheme 4
Clv,   reason is that scheme 6 produces a list of integrals that is distributed over nodes and the
clv,   authors of the above-mentioned modules nor anyone else have taken the time to adjust the
clv,   integral processing routines accordingly (if anyone is inclined to do so: this is not rocket
clv    science, just have a look on how RELCCSD does this in ccints.F).
      if (DOLUCT.OR.DOCIM) then
         ISTRAT = 4
      else
         ISTRAT = 6
      endif

c
c     activate/deactivate fine sorting of 1HT data for scheme6
c
      is6sort           = 0
      do_finesort       = .false.
      set_auto_finesort = .false.
c

      ITRA_INTFLG = INTGEN
      IGTINT = IBTAND(INTGEN/8,1)
      ITRA_INTFL2 = -8
      ITRA_INTFL4 = -8
      IPAR4BS = -1
      NCORE2 = -1
      DO I = 1,NFSYM
        WRITE(TRA_CORSTR(I),'(A)') 'not specified'
        WRITE(TRA_CORSTR2(I),'(A)') 'not specified'
        WRITE(TRA_INDSTR(I),'(A)') 'energy -10.0 20.0 1.0'
        DO J = 1,2
          WRITE(TRA2_INDSTR(J,I),'(A)') 'energy -10.0 20.0 1.0'
        ENDDO
        DO J = 1,4
          WRITE(TRA4_INDSTR(J,I),'(A)') 'energy -10.0 20.0 1.0'
        ENDDO
      ENDDO
      TRA_ANTIS  = .FALSE.
      TRA_EXACORR = .FALSE.
      TH_CHOLESKY = -1.0D0
      TRA_ASCII = .FALSE.
      NO2IND = .FALSE.
      NO4IND = .FALSE.
      NOPAIR = .TRUE.
      IF (SPINFR) THEN
         PRPTRA = .FALSE. ! Phase factor scaling needed for spinfree is not yet debugged for properties (see in traprp)
      ELSE
         PRPTRA = .TRUE.  ! Always on to allow intensity calculation in RELCCSD
      END IF
      RCORBS = .FALSE.
      NOMDCINT = .FALSE. ! will be reset depending on where PAMTR1 is called
      NOMDCINT_TRAINP = .FALSE. ! for **MOLTRA input, used to set NOMDCINT for DOTRA call of PAMTRA
      MDCSCAT = .FALSE.
!.s/sya,2007.0907/ SK - 30-11-2009
      PRPSYA = .FALSE.
!.q
      THROUT = 1.0D-14
      SCRTRA = 1.0D-14
      MOFILE_TRAINP = 'UNKNOWN'
C
C     Process input for MOLTRA
C     ========================
C
C
C     Read menu file
C     ==============
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',ERR=910,END=920) WORD
      CALL UPCASE(WORD)
      NEWDEF = (WORD .EQ. '**MOLTR')
      IF (NEWDEF) THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
        CALL QUIT('Error reading LUCMD, no transformation input found')
  920 CONTINUE
        WRITE (LUPRI,'(/A)') 
     &  'No input for integrals transformation, using defaults'
        GOTO 400
  930 CONTINUE
      WORD1 = WORD
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
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
     &                      17,18,19,20,21,22,23,24,25,26,27,28,29), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in TRAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in TRAINP.')
    1          CONTINUE
C&&&& PRINT - print level in 4-index transformation module
                  READ(LUCMD,*) IPRTRA
               GO TO 100
    2          CONTINUE
C&&&& ACTIVE - Range of ACTIVE orbitals in 4-index transformation module
                  IF (LSPES) THEN
                     WRITE (LUPRI,*) ' Warning : encountered .ACTIVE '//
     &                               'block twice , using first block'
                     GO TO 100
                  ENDIF
                  LSPES = .TRUE.
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') TRA_INDSTR(I)
                  ENDDO
C                 Allow for one-line input also for groups with inversion symmetry
                  IF (NFSYM.EQ.2) THEN
                     PROMPT = TRA_INDSTR(2)
                     IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#' 
     &                   .OR. PROMPT .EQ. '*') THEN
                        TRA_INDSTR(2) = TRA_INDSTR(1)
                        BACKSPACE(LUCMD)
                     END IF
                  END IF
               GO TO 100
    3          CONTINUE
C&&&& INTFLG - Integral classes not included in transformation
                 IF (IGTINT .EQ. 1) THEN
                    READ(LUCMD,*,IOSTAT=IOS) 
     &                ILLINT,ILSINT,ISSINT,IGTINT
                    IF (IOS.NE.0) THEN
                     WRITE(LUPRI,*)
     &                'TRAINP: Error in reading .INTFLG data,'//
     &                ' ILLINT,ILSINT,ISSINT,IGTINT:'
                     WRITE(LUPRI,*) 'read values :',
     &               ILLINT,ILSINT,ISSINT,IGTINT
                     CALL FLSHFO(LUPRI)
                     CALL QUIT(
     &               'TRAINP: Error in reading 4 .INTFLG paramaters')
                    ENDIF
                 ELSE
                    READ(LUCMD,*,IOSTAT=IOS) ILLINT,ILSINT,ISSINT
                    IF (IOS.NE.0) THEN
                      WRITE(LUPRI,*)
     &                'TRAINP: Error in reading .INTFLG data,'//
     &                ' ILLINT,ILSINT,ISSINT:'
                     WRITE(LUPRI,*) 'read values :',
     &               ILLINT,ILSINT,ISSINT
                     CALL FLSHFO(LUPRI)
                     CALL QUIT(
     &               'TRAINP: Error in reading 3 .INTFLG paramaters')
                    ENDIF
                 END IF
                  ITRA_INTFLG = ILLINT + 2*ILSINT + 4*ISSINT + 8*IGTINT
               GO TO 100
    4          CONTINUE
C&&&& "ANTIS " - Anti-Symmetrize integrals
C     aug 2001 hjaaj TODO: disabled this option because not implemented
                  TRA_ANTIS = .TRUE.
               GO TO 100
    5          CONTINUE
                  CALL QUIT ('Option disabled')
C&&&& 4INDEX
C     Ranges of active orbitals in 4-index transformation module
C     specified for index 1 to 4 and Fermion irrep 1 and 2
C
C     FORMAT:
C             INDEX1:
C             1..4, 6..9, 5,73..111
C             1..5, 6..11, 73..11
C             INDEX2:
C             1..9
C             5..111
C             etc....
C
                  IF (LSPES4) THEN
                     WRITE (LUPRI,*) ' Warning : encountered .4INDEX '//
     &                               'block twice , using first block'
                     GO TO 100
                  ENDIF
                  LSPES4 = .TRUE.
                  DO J=1,4
                     READ(LUCMD,'(A7)') LINE
                     IF(LINE(1:5).NE.'INDEX') GOTO 1010
                     DO I=1,NFSYM
                        READ(LUCMD,'(A)') TRA4_INDSTR(J,I)
                     ENDDO
                  ENDDO
               GO TO 100
    6          CONTINUE
C&&&& Use exacorr transformation routine
                  TRA_EXACORR = .TRUE.
                  ISTRAT = 4 ! This corresponds to non-distributed integrals
               GO TO 100
    7          CONTINUE
C&&&& POSITRONS - allow positronic solutions
                  NOPAIR = .FALSE.
               GO TO 100
    8          CONTINUE
C&&&& "CORE  " - Range of CORE orbitals in 2-index transformation module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') TRA_CORSTR(I)
                  ENDDO
               GO TO 100
    9          CONTINUE
C&&&& SCHEME - choose method in 4-index transformation
                  READ(LUCMD,*,IOSTAT=IOS) ISTRAT
                  IF (IOS.NE.0) THEN
                   WRITE(LUPRI,*) 'TRAINP: read ISTRAT=',ISTRAT
                   CALL QUIT(
     &             'TRAINP: Error in reading **MOLTRA .SCHEME')
                  ENDIF
               GO TO 100
   10          CONTINUE
C&&&& THROUT - Treshold for output of integrals to MDCINT
                  READ(LUCMD,*,IOSTAT=IOS) THROUT
                  IF (IOS.NE.0) THEN
                   WRITE(LUPRI,*) 'TRAINP: read THROUT=',THROUT
                   CALL QUIT(
     &             'TRAINP: Error in reading **MOLTRA .THROUT')
                  ENDIF
               GO TO 100
   11          CONTINUE
C&&&& NO2IND - Skip 2-index transformation of effective Fock matrix
                  NO2IND = .TRUE.
               GO TO 100
   12          CONTINUE
C&&&& NO4IND - Skip 4-index transformation 2-electron integrals
                  NO4IND = .TRUE.
               GO TO 100
   13          CONTINUE
                  CALL QUIT ('Option disabled')
C&&&& 2INDEX
C     Ranges of active orbitals in 2-index transformation module
C     specified for index 1 to 2 and Fermion irrep 1 and 2
C
C     FORMAT:
C             INDEX1:
C             1..4, 6..9, 5,73..111
C             1..5, 6..11, 73..11
C             INDEX2:
C             1..9
C             5..111
C
                  IF (LSPES2) THEN
                     WRITE (LUPRI,*) ' Warning : encountered .2INDEX '//
     &                               'block twice , using first block'
                     GO TO 100
                  ENDIF
                  LSPES2 = .TRUE.
                  DO J=1,2
                     READ(LUCMD,'(A7)') LINE
                     IF(LINE(1:5).NE.'INDEX') GOTO 1010
                     DO I=1,NFSYM
                        READ(LUCMD,'(A)') TRA2_INDSTR(J,I)
                     ENDDO
                  ENDDO
               GO TO 100
   14          CONTINUE
C&&&& INTFL2 - Integral classes not included in 2-index transformation
                  IF (IGTINT .EQ. 1) THEN
                     READ(LUCMD,*,IOSTAT=IOS) 
     &               ILLINT,ILSINT,ISSINT,IGTINT
                     IF (IOS.NE.0) THEN
                       CALL QUIT(
     &      'TRAINP: Error in reading  ILLINT,ILSINT,ISSINT,IGTINT'//
     &      ' for INTFL2 !')
                     ENDIF
                  ELSE
                     READ(LUCMD,*,IOSTAT=IOS) 
     &               ILLINT,ILSINT,ISSINT
                     IF (IOS.NE.0) THEN
                       CALL QUIT(
     &     'TRAINP: Error in reading ILLINT,ILSINT,ISSINT for INTFL2!')
                     ENDIF
                  END IF
                  ITRA_INTFL2 = ILLINT + 2*ILSINT + 4*ISSINT + 8*IGTINT
               GO TO 100
   15          CONTINUE
C&&&& INTFL4 - Integral classes not included in 4-index transformation
                  IF (IGTINT .EQ. 1) THEN
                     READ(LUCMD,*,IOSTAT=IOS) 
     &               ILLINT,ILSINT,ISSINT,IGTINT
                     IF (IOS.NE.0) THEN
                       CALL QUIT(
     &      'TRAINP: Error in reading ILLINT,ILSINT,ISSINT,IGTINT'
     &      //' for INTFL4 !')
                     ENDIF
                  ELSE
                     READ(LUCMD,*,IOSTAT=IOS) 
     &               ILLINT,ILSINT,ISSINT
                     IF (IOS.NE.0) THEN
                       CALL QUIT(
     &      'TRAINP: Error in reading ILLINT,ILSINT,ISSINT'
     &       //' for INTFL4 !')
                     ENDIF
                  END IF
                  ITRA_INTFL4 = ILLINT + 2*ILSINT + 4*ISSINT + 8*IGTINT
               GO TO 100
   16          CONTINUE
C&&&& PRPTRA - Transform property integrals to MO basis
                  PRPTRA = .TRUE.
               GO TO 100
   17          CONTINUE
C&&&& SCREEN - Screening threshold
                  READ(LUCMD,*) SCRTRA
               GO TO 100
   18          CONTINUE
C&&&& RCORBS - Recanonize orbitals before transforming
                  RCORBS = .TRUE.
               GO TO 100
   19          CONTINUE
C&&&& NOMDCINT - Don't write MDCINT file, instead save scratch files
C              and write an additional 4INDINFO file which contains various
C              pointers'n'stuff. This options is used within the MCSCF
C              program to avoid the generation of the enormous MDCINT file.
                  NOMDCINT_TRAINP = .TRUE.
               GO TO 100
   20          CONTINUE
C&&&& NOSCATTER - Don't Scatter MDCINT/4IND to all nodes
                  MDCSCAT = .FALSE.
               GO TO 100
   21          CONTINUE
C&&&& PAR4BS  - size of batches for strategy 4 (parallel only)
                  READ(LUCMD,*) IPAR4BS
               GO TO 100
   22          CONTINUE
C&&&& SCATTER - Do Scatter MDCINT/4IND to all nodes
                  MDCSCAT = .TRUE.
               GO TO 100
   23          CONTINUE
C&&&& MDCINT - Force writing of the MDCINT file in MOLFDIR format
                  NOMDCINT_TRAINP = .FALSE.
               GO TO 100
   24          CONTINUE
!&&&& "HSCC  " - HSMRCC interface
                  DOHSCC = .TRUE.
               GO TO 100
   25          CONTINUE
!&&&& HTSORT - MOLTRA scheme 6 stuff
               READ(LUCMD,*) IS6SORT
#ifndef VAR_PFS
               if (is6sort.eq.0) then
                  do_finesort       = .false.
                  set_auto_finesort = .false.
               elseif (is6sort.eq.1) then
                  do_finesort       = .true.
                  set_auto_finesort = .false.
               endif
#else
               write (*,*) 'MPI-IO enabled! 1HT sort deactivated'
               do_finesort       = .false.
               set_auto_finesort = .false.
#endif
               GO TO 100
   26          CONTINUE
!&&&& "ASCII " - Write ASCII file with transformed MO integrals
                  TRA_ASCII = .TRUE.
               GO TO 100
   27          CONTINUE
!&&&& "MOFILE" - which file to read MO coefficients from
!      note: list directed read will ignore initial blanks, so MOFILE_TRAINP
!            will always begin with the first non-blank character in input line.
                  READ(LUCMD,*) MOFILE_TRAINP
               GO TO 100
   28          CONTINUE
C&&&& "CORE2 " - Range of open shell CORE orbitals in 2-index transformation module
                  WRITE (LUPRI,'(A)') '* use open-shell core:'
                  READ(LUCMD,*) NCORE2
                  WRITE (LUPRI,'(I4,A)') NCORE2,' el. in'
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') TRA_CORSTR2(I)
                     WRITE (LUPRI,'(A)') TRA_CORSTR2(I)
                  ENDDO
               GO TO 100
   29          CONTINUE
C&&&& ."CHOLES" - CHOLESKY threshold
                  READ(LUCMD,*) TH_CHOLESKY
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in TRAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in TRAINP.')
            END IF
      END IF
C
C     Process input for various program sections
C     ==========================================
C
  300 CONTINUE
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 300
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 310 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (400,312,313,314), I
            END IF
  310    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 400
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in TRAINP.')
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in TRAINP, error in prompt.')
      END IF
  312 CONTINUE
        CALL TRPINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 300
  313 CONTINUE
#ifdef MOD_HSCC
      IF (DOHSCC) THEN
        CALL HSCCINP(WORD,.false.)
      ELSE
        WRITE(LUPRI,'(2X,A)') 
     &  '*HSCC input found, but not processed (DOHSCC=.false.).'
        GOTO 400
      ENDIF
#else
C       prevent dirac.x running into an infinite loop here
C       if "*HSCC" is in the input but MOD_HSCC is not defined !!!
        CALL QUIT('*HSCC specified, '//
     &     'but HSCC module is not included in this version')
#endif
        GO TO 300
  314 CONTINUE
        GO TO 300
C
  400 CONTINUE

      CALL TRPINP(WORD,.TRUE.,WORK,LWORK)  ! always call this routine, sets transformation for dipole integrals
!Miro: set RESET to .true. to enable geometry optimization

C
C     Process defaults for 2- and 4-index transformations
C
      IF (LSPES.AND..NOT.LSPES2) THEN
        DO I = 1, NFSYM
           DO J = 1, 2
             TRA2_INDSTR(J,I) = TRA_INDSTR(I)
           ENDDO
        ENDDO
      ENDIF
C
      IF (LSPES.AND..NOT.LSPES4) THEN
        DO I = 1, NFSYM
           DO J = 1, 4
             TRA4_INDSTR(J,I) = TRA_INDSTR(I)
           ENDDO
        ENDDO
      ENDIF
C
      IF (ITRA_INTFL2.EQ.-8) THEN
         ITRA_INTFL2 = ITRA_INTFLG
      ENDIF
C
      IF (ITRA_INTFL4.EQ.-8) THEN
         ITRA_INTFL4 = ITRA_INTFLG
      ENDIF

      IF (LEVYLE.OR.BSS.or.x2c.OR.DO4C2C) THEN
         ITRA_INTFL2 = MOD(ITRA_INTFL2,2)
         ITRA_INTFL4 = MOD(ITRA_INTFL4,2)
      ENDIF
C
      IF(NOPAIR) THEN
        IC   = 1
        DO I = 1,NFSYM
          IR(1,I) = 1
          IR(2,I) = NFBAS(I,1)
        ENDDO
      ELSE
        IC   = 0
        DO I = 1,NFSYM
          IR(1,I) = -NFBAS(I,2)
          IR(2,I) =  NFBAS(I,1)
        ENDDO
      ENDIF
C
C     Print section
C     =============
C
      WRITE (LUPRI,'(//)')
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   'TRAINP: Set-up for index transformation'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(/1X,A,I3)') '* General print level   : ',IPRTRA
      IF(RCORBS) WRITE(LUPRI,'(1X,A)') '* Recanonization of orbitals.'
      IF(NOPAIR) WRITE(LUPRI,'(1X,A)') '* Electronic orbitals only.'
      IF(NO2IND) WRITE(LUPRI,'(1X,A)') '* Skip 2-index transformation.'
      IF(NO4IND) WRITE(LUPRI,'(1X,A)') '* Skip 4-index transformation.'
      IF(NOMDCINT_TRAINP) WRITE(LUPRI,'(1X,A)')
     &     '* MDCINT is not going to be written, 4IND1* is kept.'
      IF(MDCSCAT) WRITE(LUPRI,'(1X,2A)')
     &     '* MDCINT (and in principle the 4IND*-files) ',
     &     '  will be scattered to all nodes.'
      IF ( MOFILE_TRAINP .NE. 'UNKNOWN' ) THEN
         WRITE(LUPRI,'(2A)')
     &   ' * MO coeffients will be read from ',MOFILE_TRAINP
      END IF
      WRITE(LUPRI,'(1X,A)') '* Total active space.'
      NTVEC(1) = 0
      DO I = 1,NFSYM
        NVEC = 0
        WRITE(LUPRI,'(4X,A,A3)') 'Fermion ircop:',FREP(I)
        CALL  NUMLST(TRA_INDSTR(I),IDUMMY,NFBAS(I,IC),
     &             IR(1,I),IR(2,I),I,NVEC)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A)') ' No explicit orbitals specified '
        ELSE
          WRITE(LUPRI,'(4X,A72)') TRA_INDSTR(I)
        ENDIF
        NTVEC(1) = NTVEC(1) + NVEC
      ENDDO
cluuk IF (NTVEC(1).EQ.0) THEN
cluuk    CALL QUIT('TRAINP: Error defining active orbitals')
cluuk ENDIF
C
      WRITE(LUPRI,'(/1X,A)')
     &   '* Set-up for 2-index transformation'
      IF (MOD(ITRA_INTFL2,2).NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* LL Integrals not included in core Fock-matrix'
      IF (MOD(ITRA_INTFL2/2,2).NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* LS Integrals not included in core Fock-matrix'
      IF (MOD(ITRA_INTFL2/4,2).NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* SS Integrals not included in core Fock-matrix'
      IF (GAUNT.AND.MOD(ITRA_INTFL2/8,2).NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* Gaunt Integrals not included in core Fock-matrix'
C
      WRITE(LUPRI,'(1X,A)') '* Active spaces:'
      CALL IZERO (NTVEC,2)
      DO I = 1,NFSYM
        WRITE(LUPRI,'(4X,A,A3)') 'Fermion ircop:',FREP(I)
        DO J = 1,2
          NVEC = 0
          CALL  NUMLST(TRA2_INDSTR(J,I),IDUMMY,NFBAS(I,IC),
     &               IR(1,I),IR(2,I),I,NVEC)
          IF(NVEC.EQ.0) THEN
            WRITE(LUPRI,'(4X,A,I2)')
     &      ' No explicit orbitals specified for index ',J
          ELSE
            WRITE(LUPRI,'(4X,A,I2,A,A72)')
     &      '- Index ',J,':  ',TRA2_INDSTR(J,I)
          ENDIF
          NTVEC(J) = NTVEC(J) + NVEC
        ENDDO
      ENDDO
C
      WRITE(LUPRI,'(/1X,A)')
     &   '* Set-up for 4-index transformation'
      WRITE(LUPRI,'(1X,A,I3)') '* Following scheme      : ',ISTRAT
      GOTO (501,502,503,504,505,506) ISTRAT
  501 CONTINUE
        WRITE(LUPRI,'(A)')
     &    ' - all half-transformed integrals (ij|rs) kept in core'
        GOTO 500
  502 CONTINUE
        WRITE(LUPRI,'(A)')
     &    ' - accumulate (ij|ks) in core; calculate MP2 energy'
        GOTO 500
  503 CONTINUE
        WRITE(LUPRI,'(A)')
     &    ' - write half-transformed integrals (ij|rs) to disk'

        if (set_auto_finesort) then
           WRITE(LUPRI,'(A)')
     &    ' - sorting of intermediate 1HT integrals to be determined'//
     &    ' on the fly from the dimension of the problem'
        else
           if (do_finesort) then
              WRITE(LUPRI,'(A)')
     &    ' - sorting of intermediate 1HT integrals is enabled'
           else
              WRITE(LUPRI,'(A)')
     &    ' - sorting of intermediate 1HT integrals is disabled'
           endif
        endif

        GOTO 500
  504 CONTINUE
        WRITE(LUPRI,'(A)')
     &    ' - write (rs)-batches of half-transformed integrals ',
     &    ' (ij|rs) to disk; parallel scheme'
        GOTO 500
  505 CONTINUE
        WRITE(LUPRI,'(A)')
     &    ' - scheme 5 is meant for MP2 properties'
        CALL QUIT ('Specified wrong scheme for index transform')
        GOTO 500
  506 CONTINUE
        WRITE(LUPRI,'(A)')
     &    ' - write half-transformed integrals (ij|rs) to disk'

        if (set_auto_finesort) then
           WRITE(LUPRI,'(A)')
     &    ' - sorting of intermediate 1HT integrals to be determined'//
     &    ' on the fly from the dimension of the problem'
        else
           if (do_finesort) then
              WRITE(LUPRI,'(A)')
     &    ' - sorting of intermediate 1HT integrals is enabled'
           else
              WRITE(LUPRI,'(A)')
     &    ' - sorting of intermediate 1HT integrals is disabled'
           endif
        endif

        GOTO 500
  500 CONTINUE
      IF(SCRTRA.GT.D0) THEN
        WRITE(LUPRI,'(1X,A,1P,E8.2)') '* Screening threshold :',SCRTRA
      ELSE
        WRITE(LUPRI,'(1X,A)') '* No screening.'
      ENDIF
      WRITE(LUPRI,'(1X,A,1P,E8.2)') '* MO integral threshold :',THROUT
      IF (MOD(ITRA_INTFL4,2).NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* LL Integrals not transformed'
      IF (MOD(ITRA_INTFL4/2,2).NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* LS Integrals not transformed.'
      IF (MOD(ITRA_INTFL4/4,2).NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* SS Integrals not transformed.'
      IF (MOD(ITRA_INTFL4/8,2).NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* Gaunt Integrals not transformed.'
      IF (TRA_ANTIS) WRITE(LUPRI,'(1X,A)')
     +     '* Calculating Anti-symmetrized integrals.'
      WRITE(LUPRI,'(1X,A)') '* Active spaces:'
      CALL IZERO (NTVEC,4)
      DO I = 1,NFSYM
        WRITE(LUPRI,'(4X,A,A3)') 'Fermion ircop:',FREP(I)
        DO J = 1,4
          NVEC = 0
          CALL  NUMLST(TRA4_INDSTR(J,I),IDUMMY,NFBAS(I,IC),
     &               IR(1,I),IR(2,I),I,NVEC)
          IF(NVEC.EQ.0) THEN
            WRITE(LUPRI,'(4X,A,I2)')
     &      ' No explicit orbitals specified for index ',J
          ELSE
            WRITE(LUPRI,'(4X,A,I2,A,A72)')
     &      '- Index ',J,':  ',TRA4_INDSTR(J,I)
          ENDIF
          NTVEC(J) = NTVEC(J) + NVEC
        ENDDO
      ENDDO
C
C     Batch size of strategy 4.
C
      IF (PARCAL .AND. ISTRAT.EQ.4) THEN
        IF (IPAR4BS .GT. 0) THEN
          WRITE(LUPRI,'(1X,A,I4)')
     &    '* Size of batches for strategy 4 (parallel only): ', IPAR4BS
        ELSE
          WRITE(LUPRI,'(1X,A,I4)')
     &    '* Size of batches for strategy 4 determined dynamically.'
        END IF
      END IF
  999 CONTINUE
      RETURN
 1010 CONTINUE
      WRITE(LUPRI,'(A)') 'TRAINP: Error in reading .INDEX'
      WRITE(LUPRI,'(A,A)') 'Read: ',LINE
      CALL QUIT('TRAINP: Error in reading .INDEX')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck trpinp */
      SUBROUTINE TRPINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for operators that will be transformed to MO-basis
C
C     Written by L. Visscher - Dec 1996 (adapted EXPINP)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "dummy.h"

#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbprp.h"
#include "dcbtrp.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcbtra.h"
! Added temporarily to test for bss (where transform is not working)
#include "dcbham.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER PNAME*16
      DIMENSION PNAME(3)
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.PRINT ','.OPERAT','.MORANG','.NMQM  ', !MKN
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
      INPERR = 0
C
C
C     Initialize /CBIEXP/
C     ===================
C
CMKN
      DONMQM = .FALSE.
CMKN
      IPRTRP   = IPRTRA
!     Always transform the dipole length operator
      IF (.NOT. BSS) THEN
        PNAME(1) = 'XDIPLEN'
        PNAME(2) = 'YDIPLEN'
        PNAME(3) = 'ZDIPLEN'
        DO I=1,3
          CALL XPRIND(PNAME(I),1,1,D1,PNAME(I),INDXPR,ISYXPR,
     &               ITRXPR,IPRTRP)
          CALL OP1IND('TRPINP',INDTRP,LTRPP,NTRPP,INDXPR,MAXTRP)
        END DO
      END IF

C
C     Process input from CBITRP
C     =========================
C
      NEWDEF = (WORD .EQ. '*PRPTRA')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
  110       CONTINUE
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
     *            '" not recognized in TRPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in TRPINP.')
    1          CONTINUE
C&&&& PRINT:  Print level
                  READ(LUCMD,*) IPRTRP
               GO TO 100
    2          CONTINUE
C&&&& OPERATOR: Define operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        IPRTRP)
                  CALL OP1IND('TRPINP',INDTRP,LTRPP,NTRPP,INDXPR,
     &                        MAXTRP)
                  PRPTRA = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    3          CONTINUE
!.s/sya,2007.0907,add/ SK - 30-11-2009
!&&&& MORANG - Range of MOs in property integral transformation
                  PRPSYA = .TRUE.
                  DO I = 1, NFSYM
                    READ(LUCMD,'(A)') TRA2_INDPRP(I)
                  END DO
!.q
               GO TO 100
    4          CONTINUE
CMKN NMQM: Define operator Nuclear Magnetic Quadruple Moment
                  DONMQM = .TRUE.
                  PRPTRA = .TRUE.
                  ICHANG = ICHANG + 1
CMKN NMQM: Defination added by M. K. Nayak on July 19 2018
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
     *            '" not recognized in TRPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in TRPINP.')
            END IF
      END IF
  300 CONTINUE
CMKN
      IF(DONMQM) THEN
        CALL DEF_PRP_NMQM(LTRPP,NTRPP,MAXTRP,IPRTRP)
      ENDIF
CMKN
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'TRPINP: Property integral transformation'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I5)') '* Print level:',IPRTRP
      WRITE(LUPRI,'(1X,A)')
     &  '*The following operators will be transformed:'
      IF(NTRPP.GT.0) THEN
        DO I = 1,NTRPP
          INDXPR = LTRPP(I)
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
      IF(EFG) THEN
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Electric field gradient'
      ENDIF
      IF(NQCC) THEN
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Nuclear Quadrupole Coupling Constant'
      ENDIF
      IF(DONMQM) THEN
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Nuclear Magnetic Quadruple Moment'
      ENDIF
C
      CALL PRSYMB(LUPRI,'-',75,0)
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck resinp */
      SUBROUTINE RESINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for RESOLVE module (GOSCI)
C
C     Written by Jon K. Laerdahl and T. Saue (Dec 1996)
C     Last changes: Miro Ilias, Aug.2016 (for DFT-COSCI)
C
C*****************************************************************************
      use include_dcbtra_h, only : set_tra_ascii
      use dirac_cfg
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
      PARAMETER(D0 = 0.0D0)
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*2 WS
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbres.h"
#include "dcbham.h"
      DIMENSION WORK(LWORK),IR(2,2)
C
      SAVE SET
      DATA TABLE /'.PRINT ','.INTFLG','.SCREEN','.SCHEME','xTHROUT',
     &            '.ASCII ','.NOVXC ','.COREDN','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
#include "ibtfun.h"
#include "memint.h"
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
C     Initialize /DCLRES, DCBRES, DCIRES/
C     ===================================
C
      ISTRES = 4
      IPRRES = 0
      ILLINT = IBTAND(INTGEN,1)
      ILSINT = IBTAND(INTGEN/2,1)
      ISSINT = IBTAND(INTGEN/4,1)
      IGTINT = IBTAND(INTGEN/8,1)
      SCRRES = 1.D-14
      RESOUT = 1.D-14
      NOVXC  = .FALSE.
      COREDENS = .FALSE.
C
C     Process input for RESOLV
C     ========================
C
      NEWDEF = (WORD .EQ. '*RESOLV')
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
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8),I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in RESINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in RESINP.')
    1          CONTINUE
C&&&& PRINT - print level in resolve module
                  READ(LUCMD,*) IPRRES
               GO TO 100
    2          CONTINUE
C&&&& Integral classes included in 4-indextransformation
                  IF (IGTINT .EQ. 1) THEN
                     READ(LUCMD,*) ILLINT,ILSINT,ISSINT,IGTINT
                  ELSE
                     READ(LUCMD,*) ILLINT,ILSINT,ISSINT
                  END IF
               GO TO 100
    3          CONTINUE
C&&&& SCREEN: Screening threshold for 4-index transformation
                  READ(LUCMD,*) SCRRES
               GO TO 100
    4          CONTINUE
C&&&& SCHEME - Choose transformation scheme
                  READ(LUCMD,*) ISTRES
               GO TO 100
    5          CONTINUE
C&&&& THROUT - Treshold for output integrals
CTROND: disabled, because THROUT is not implemented 
                  READ(LUCMD,*) RESOUT
               GO TO 100
!Miro: .ASCII keyword for writing the "MO_integrals.txt" file
    6          CONTINUE
                  call set_tra_ascii(.true.)
               GO TO 100
!Miro: .NOVXC keyword for not adding DFT-Vxc terms into Fij (rough approx)
    7          CONTINUE
                  NOVXC = .TRUE.
               GO TO 100
!Miro: .COREDN keyword to have only core density for Vxc (approximation)
    8          CONTINUE
                  COREDENS= .TRUE.
               GO TO 100
! .XXXXX
    9          CONTINUE
               GO TO 100
! .XXXXX
   10          CONTINUE
               GO TO 100

            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in RESINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in RESINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Processing section
C     ------------------
C
CMI ... TODO:  should be adapted properly for 2c-4c/4c-2c stuff
      IF (LEVYLE.OR.BSS.or.x2c.OR.DO4C2C) THEN
        ILSINT = 0
        ISSINT = 0
        IGTINT = 0
      ENDIF
      INTRES = ILLINT + 2*ILSINT + 4*ISSINT + 8*IGTINT
C
C     Print section
C     =============
C
      IF(.NOT.DORES) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   'RESINP: Set-up for resolution of open-shell states:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I3)') '* General print level   : ',IPRRES
      WRITE(LUPRI,'(1X,A,I3)')
     &  '* 4-index transformation follows scheme      : ',ISTRES
      GOTO (501,502,503,504) ISTRES
  501 CONTINUE
        WRITE(LUPRI,'(A)')
     &    ' - all half-transformed integrals (ij|rs) kept in core'
        GOTO 500
  502 CONTINUE
        WRITE(LUPRI,'(A)')
     &    ' - accumulate (ij|ks) in core; calculate MP2 energy'
        GOTO 500
  503 CONTINUE
        WRITE(LUPRI,'(A)')
     &    ' - write half-transformed integrals (ij|rs) to disk'
        GOTO 500
  504 CONTINUE
        WRITE(LUPRI,'(A)')
     &    ' - write (rs)-batches of half-transformed integrals ',
     &    ' (ij|rs) to disk; parallel scheme'
  500 CONTINUE
      IF(SCRRES.GT.D0) THEN
        WRITE(LUPRI,'(1X,A,1P,E8.2)') '* Screening threshold : ',SCRRES
      ELSE
        WRITE(LUPRI,'(1X,A)') '* No screening.'
      ENDIF
      IF (ILLINT.NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* LL integrals not included.'

! .. always print out what integrals are included
         IF (ILSINT.NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* LS integrals not included.'
         IF (ISSINT.NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* SS integrals not included.'
         IF (IGTINT.NE.1) WRITE(LUPRI,'(1X,A)')
     +     '* Gaunt integrals not included.'

!Miro: printouts concerning DFT-COSCI
      IF (dirac_cfg_dft_calculation) THEN
         IF (NOVXC) THEN
            WRITE(LUPRI,"(2X,A)")
     &     "No Vxc term in transformed Fij elements for COSCI!"
            WRITE(LUPRI,"(2X,A)")
     &     "This is crude approximation for DFT-COSCI ..."
         ENDIF
         IF (COREDENS) THEN
            WRITE(LUPRI,"(2X,A)")
     &      "Only core-density for Vxc in Fij terms ... approximation "
         ENDIF
         IF (PARCAL) CALL QUIT('DFT-COSCI not parallelized')
      ENDIF
      IF ((NOVXC.OR.COREDENS).AND..NOT.dirac_cfg_dft_calculation) THEN
          WRITE(LUPRI,"(2X,A)")
     &    "COSCI-DFT works only with open-shell DFT !"
      ENDIF

 999  CONTINUE
      RETURN
      END
