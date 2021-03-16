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
C FILE    : gptrygve.F
C
C General Purpose utility routines (originally) introduced in Dalton by Trygve Helgaker
C or one of his students
C
! Dirac uses INT_STAR8 and Dalton used VAR_INT64;
! If one set then set the other:
#ifdef INT_STAR8
#define VAR_INT64
#else
#ifdef VAR_INT64
*define INT_STAR8
#endif
#endif
!
C  /* Deck header */
      SUBROUTINE HEADER(HEAD,IN)
      CHARACTER HEAD*(*)
#include "priunit.h"
      integer, intent(in) :: in
C
      LHEAD = LNBLNK(HEAD)
      IF (IN .GE. 0) THEN
         INDENT = IN + 1
      ELSE
         INDENT = MAX(1,(80 - LHEAD)/2 + 1)
      END IF
      WRITE (LUPRI, '(//,150A)') (' ',I=1,INDENT), HEAD(1:LHEAD)
      WRITE (LUPRI, '(   150A)') (' ',I=1,INDENT), ('-',I=1,LHEAD)
      WRITE (LUPRI, '()')
      CALL FLSHFO(LUPRI)
      RETURN
      END
      SUBROUTINE HEADER2(HEAD,IN,unitn)
      CHARACTER HEAD*(*)
      integer, intent(in) :: unitn, in
C
      LHEAD = LNBLNK(HEAD)
      IF (IN .GE. 0) THEN
         INDENT = IN + 1
      ELSE
         INDENT = MAX(1,(80 - LHEAD)/2 + 1)
      END IF
      WRITE (unitn, '(//,150A)') (' ',I=1,INDENT), HEAD(1:LHEAD)
      WRITE (unitn, '(   150A)') (' ',I=1,INDENT), ('-',I=1,LHEAD)
      WRITE (unitn, '()')
      RETURN
      END
C  /* Deck timer */
      SUBROUTINE TIMER(TEXT,TIMSTR,TIMEND)
      ! print CPU time used
      ! template:
      !   CALL TIMER('START',TIMSTR,TIMEND)
      !   ... work ...
      !   CALL TIMER('work', TIMSTR,TIMEND
      !   ... more work ...
      !   CALL TIMER('more work',TIMSTR,TIMEND)
      implicit none
#include "priunit.h"
      CHARACTER TEXT*(*)
      REAL*8    TIMSTR, TIMEND
      ! local variables:
      REAL*8    TIME, DUMMY
      INTEGER   LTEXT, LOUT, LNBLNK
      CHARACTER OUT_TEXT*50
C
      IF (TEXT(1:5) .EQ. 'START') THEN
         CALL GETTIM(TIMSTR,DUMMY)
      ELSE
         CALL GETTIM(TIMEND,DUMMY)
         TIME   = TIMEND - TIMSTR
         TIMSTR = TIMEND
         IF (TIME .GT. 0.1D0) THEN
C        Jun 06 hjaaj: do not print times less than 0.1 seconds
            LTEXT = LNBLNK(TEXT)
            LTEXT = MIN(30,LTEXT)
            OUT_TEXT = '>>> CPU time used in '//TEXT(1:LTEXT)//' is'
            LOUT = 21 + LTEXT + 3
            CALL TIMTXT(OUT_TEXT(1:LOUT),TIME,LUPRI)
            flush(LUPRI)
         END IF
      END IF
      RETURN
      END
      SUBROUTINE TIMER2(TEXT,TIMSTR,TIMEND)
      ! print CPU and WALL time used
      ! template:
      !   CALL TIMER('START',TIMSTR,TIMEND)
      !   ... work ...
      !   CALL TIMER('work', TIMSTR,TIMEND
      !   ... more work ...
      !   CALL TIMER('more work',TIMSTR,TIMEND)
      implicit none
#include "priunit.h"
      CHARACTER TEXT*(*)
      REAL*8    TIMSTR(2), TIMEND(2)
      ! local variables:
      REAL*8    CTIME, WTIME
      INTEGER   LTEXT, LOUT, LNBLNK
      CHARACTER OUT_TEXT*50
C
      IF (TEXT(1:5) .EQ. 'START') THEN
         CALL GETTIM(TIMSTR(1),TIMSTR(2))
      ELSE
         CALL GETTIM(TIMEND(1),TIMEND(2))
         CTIME     = TIMEND(1) - TIMSTR(1)
         TIMSTR(1) = TIMEND(1)
         WTIME     = TIMEND(2) - TIMSTR(2)
         TIMSTR(2) = TIMEND(2)
         IF (CTIME .GT. 0.1D0 .OR. WTIME .GT. 0.1D0) THEN
C        Jun 06 hjaaj: do not print times less than 0.1 seconds
            LTEXT = LNBLNK(TEXT)
            LTEXT = MIN(30,LTEXT)
            LOUT = 22 + LTEXT + 3
            OUT_TEXT = '>>> CPU  time used in '//TEXT(1:LTEXT)//' is'
            CALL TIMTXT(OUT_TEXT(1:LOUT),CTIME,LUPRI)
            OUT_TEXT = '>>> WALL time used in '//TEXT(1:LTEXT)//' is'
            CALL TIMTXT(OUT_TEXT(1:LOUT),WTIME,LUPRI)
            flush(LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck timdat */
      SUBROUTINE TIMDAT
#include "implicit.h"
#include "priunit.h"
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      INTEGER CDATE, CTIME
      CDATE = DATE()
      CTIME = CLOCK()
      WRITE(LUPRI,100) CDATE, CTIME
  100 FORMAT(/,'  Date and time (CRAY) : ',A8,2X,A8)
#endif
#if defined (SYS_CONVEX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Convex) : ',A)
#endif
#if defined (SYS_HAL)
      CHARACTER*(24) CFDATE
      CALL FDATE(CFDATE)
      WRITE (LUPRI,100) CFDATE
  100 FORMAT(/,' Date and time (HAL)    : ',A)
#endif
#if defined (SYS_AIX)
C     AIX XL FORTRAN version 2.3+
      CHARACTER*(24) fdate
      WRITE (LUPRI,100) fdate()
  100 FORMAT(/,' Date and time (IBM-AIX): ',A)
#endif
#if defined (SYS_HPUX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (HP-UX)  : ',A)
#endif
#if defined (SYS_DEC)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (DEC)  : ',A)
#endif
#if defined (SYS_IRIX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (IRIX)   : ',A)
#endif
#if defined (SYS_LINUX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Linux) : ',A)
#endif
#if defined (SYS_WINDOWS)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Windows) : ',A)
#endif
#if defined (SYS_DARWIN)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Darwin) : ',A)
#endif
#if defined (SYS_SX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Nec-SX) : ',A)
#endif
#if defined (SYS_SUN)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (SUN)   : ',A)
#endif
      RETURN
      END
C  /* Deck timpri */
      SUBROUTINE TIMPRI(TEXT,TIME,TIMALL)
#include "implicit.h"
#include "priunit.h"
      CHARACTER TEXT*6
      PARAMETER (HUN = 100.0D00)
      SAVE TREST
      DATA TREST /0.0D0/
      IF (TEXT .EQ. 'REST  ') THEN
         TIME = TREST
C 980826-hjaaj: initialize TREST for next abacus call
         TREST = 0.0D0
      END IF
      IF (TIMALL .EQ. 0.0D0) RETURN
      ITIME = NINT(HUN*TIME/TIMALL)
      IF (ITIME .GT. 0 .AND. INT(TIME) .GT. 0) THEN
         MINUTE = INT(TIME)/60
         IHOURS = MINUTE/60
         MINUTE = MINUTE - 60*IHOURS
         ISECND = NINT(TIME) - 3600*IHOURS - 60*MINUTE
         WRITE(LUPRI,100) TEXT, IHOURS, MINUTE, ISECND, ITIME
      ELSE
         TREST = TREST + TIME
      END IF
  100 FORMAT(1X,A6,'     ',I2.2,':',I2.2,':',I2.2,5X,I3,' %')
      RETURN
      END
C  /* Deck titler */
      SUBROUTINE TITLER(HEAD,A,IN)
#include "priunit.h"
      CHARACTER HEAD*(*), A*(*)
C
      LHEAD  = LNBLNK(HEAD)
      LENGTH = LHEAD
      IF (IN .EQ. 200) THEN
         LENGTH = LENGTH + 2
      ELSE IF (IN .GE. 100) THEN
         MARG = (80 - LENGTH)/2 - 1
         MARG = MIN(MARG, IN - 100)
         IF (MARG .GT. 0) MARG = MARG + 1
         LENGTH = LENGTH + 2*MARG
      END IF
      IF (IN .GE. 0 .AND. IN .LT. 100) THEN
         INDENT = IN + 1
      ELSE
         INDENT = MAX(1,(80 - LENGTH)/2 + 1)
      END IF
      IF (IN .EQ. 200) THEN
         WRITE (LUPRI, '(//150A)')
     *      (' ',I=1,INDENT),'.', ('-',I=1,LENGTH),'.'
         WRITE (LUPRI, '(150A)')
     &      (' ',I=1,INDENT),'| ', HEAD(1:LHEAD), ' |'
         WRITE (LUPRI, '(150A)')
     *      (' ',I=1,INDENT),'`', ('-',I=1,LENGTH),"'"
      ELSE IF (IN .EQ. 100) THEN
         WRITE (LUPRI, '(//150A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(150A)') (' ',I=1,INDENT), HEAD(1:LHEAD)
         WRITE (LUPRI, '(150A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE IF (IN .GT. 100) THEN
         WRITE (LUPRI, '(//150A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(150A)') (' ',I=1,INDENT),
     *      (A,I=1,MARG-1), ' ', HEAD(1:LHEAD), ' ', (A,I=1,MARG-1)
         WRITE (LUPRI, '(150A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE
         WRITE (LUPRI, '(//150A)') (' ',I=1,INDENT), HEAD(1:LHEAD)
         WRITE (LUPRI, '(150A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      END IF
      WRITE (LUPRI, '()')
      CALL FLSHFO(LUPRI)
      RETURN
      END
C  /* Deck around */
      SUBROUTINE AROUND(HEAD)
      CHARACTER HEAD*(*)
#include "priunit.h"
      LHEAD  = LNBLNK(HEAD)
      LNG    = LHEAD + 2
      IND = MAX(1,(80 - LNG)/2 + 1)
      WRITE (LUPRI,'(//150A)') (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
      WRITE (LUPRI,'(150A)')   (' ',I=1,IND), '! ', HEAD(1:LHEAD), ' !'
      WRITE (LUPRI,'(150A)')   (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
Cx    WRITE (LUPRI,'(//150A)') (' ',I=1,IND), '.', ('-',I=1,LNG), '.'
Cx    WRITE (LUPRI,'(150A)')   (' ',I=1,IND), '| ', HEAD(1:LHEAD), ' |'
Cx    WRITE (LUPRI,'(150A)')   (' ',I=1,IND), '`', ('-',I=1,LNG), ''''
      WRITE (LUPRI,'()')
      RETURN
      END
C  /* Deck stopit */
      SUBROUTINE STOPIT(SUB,PLACE,INT1,INT2)
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) SUB, PLACE
      WRITE (LUPRI,'(//3A)') ' Work space exceeded in subroutine ',
     *                         SUB,'.'
      IF (LNBLNK(PLACE) .GT. 0) THEN
         WRITE (LUPRI,'(/2A)') ' Location: ',PLACE
      END IF
      JNT1  = ABS(INT1)
      JNT2  = ABS(INT2)
      LWRKR = MAX(JNT1,JNT2)
      LWRKA = MIN(JNT1,JNT2)
      IF (INT1 .LT. 0 .OR. INT2 .LT. 0) THEN
         WRITE (LUPRI,'(/A,I10)  ') ' Space required  >',LWRKR
      ELSE
         WRITE (LUPRI,'(/A,I10)  ') ' Space required  =',LWRKR
      END IF
      WRITE (LUPRI,'( A,I10)') ' Space available =',LWRKA
      CALL QUIT('Work space exceeded.')
      END
C  /* Deck quit */
      SUBROUTINE QUIT(TEXT)
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
#include "priunit.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#endif
      CHARACTER TEXT*(*), STRING*70, BLA*1
C
#if defined (VAR_MPI)
      IF (MYTID .GT. 0) THEN
         WRITE (STRING,'(A,I4,A)') ' Slave node',MYTID,
     &   ' :  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---'
      ELSE
         WRITE (STRING,'(A)')
     &   '  Master node : --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---'
      END IF
#else
      WRITE (STRING,'(A)')
     &   '  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---'
#endif
      WRITE (LUPRI,'(/A)') STRING
      CALL TIMDAT
      WRITE (LUPRI,'(1X,A)') TEXT
#if defined (SYS_CRAY) || defined (SYS_CONVEX)\
|| defined (SYS_UNIX) || defined (SYS_AIX)\
|| defined (SYS_HPUX) || defined (SYS_IRIX)\
|| defined (SYS_NEC) || defined (SYS_DEC)\
|| defined (SYS_T3D) || defined (SYS_LINUX)\
|| defined (SYS_SUN) || defined (SYS_HPUX)\
|| defined (SYS_HAL) || defined (SYS_T90)\
|| defined (SYS_SX) || defined (SYS_DARWIN)
C     Write to stderr
      WRITE (0,'(/A/1X,A)') STRING, TEXT
      CALL FLSHFO(0)
#endif
      CALL FLSHFO(LUPRI)

#if defined (SYS_LINUX)
C aspg, replacing trace with os-provided traceback facilities. linux support for now
C 3 = pretty print traceback to stdout...
      CALL GET_TRACEBACK_INFO (3,BLA,1) 
#else
      CALL TRACE
#endif /* SYS_LINUX */

#if defined (VAR_MPI)
      ierr = -9999
#   if defined (PRG_DIRAC)
      IF (PARCAL) CALL interface_MPI_ABORT(global_communicator,IERR)
#   else
      CALL interface_MPI_ABORT(global_communicator,IERR)
#   endif
#endif

CMI ... .*F -> *.p dependency brings -DNO_SF_QUIT
#if defined (NO_SF_QUIT)
      CALL FLSHFO(LUPRI)
      STOP 1000
#else
CMI ... might be input defined restriction for some reasons
      IF (.NOT.NOSFCRASH) THEN
        CALL FLSHFO(LUPRI)
        CALL CRASH()
      ELSE
        CALL FLSHFO(LUPRI)
        STOP 1001
      ENDIF
#endif
      END

C  /* Deck opendx */
      SUBROUTINE OPENDX (LUDX,NAME,NELEM,STATUS,OLDDX)
C
C 15-Jun-1985 hjaaj
C
C Revisions :  9-Dec-1987 hjaaj (Alliant version)
C
C Purpose:
C   Open files for direct access through WRITDX and READDX routines.
C   The ....DX routines enables direct access, even when the number
C   of elements per record (the logical record length) is greater
C   than the maximum physical record length.
C   THIS IS MACHINE DEPENDENT 
C
C Input:
C  LUDX     file unit number
C  NELEM    number of integer words per logical record
C  STATUS   file status: 'OLD', 'NEW', or 'UNKNOWN'
C
C Output:
C  OLDDX    logical, true if old LUDX file was opened
C
C
      CHARACTER*(*) NAME, STATUS
      LOGICAL OLDDX, FEXIST
#include "priunit.h"
      LRECL  = 4*NELEM
C
      IF (STATUS .EQ. 'NEW') THEN
C        some systems abort if we try to open as NEW a file,
C        which already exist
         INQUIRE(FILE=NAME,EXIST=FEXIST)
         IF (FEXIST) THEN
            OPEN(UNIT=LUDX,STATUS='OLD',FILE=NAME)
            CLOSE(UNIT=LUDX,STATUS='DELETE')
         END IF
         GO TO 300
      END IF
C
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 9000
C
C     OPEN OLD FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',ERR=300,
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .TRUE.
      GO TO 600
C
  300 CONTINUE
      IF (STATUS .EQ. 'OLD') GO TO 9100
C
C     OPEN NEW FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .FALSE.
  600  CONTINUE
C
C error branches
C
 9000 CONTINUE
      WRITE (LUPRI,'(//2A/2A/A,I5)')
     *   ' --> ERROR (OPENDX) INVALID STATUS KEYWORD: ',STATUS,
     *   '                    FILE NAME   =',NAME,
     *   '                    FILE NUMBER =',LUDX
      CALL QTRACE(LUPRI)
      CALL QUIT('ERROR (OPENDX) INVALID STATUS KEYWORD')
C
 9100 CONTINUE
      WRITE (LUPRI,'(//A/2A/A,I5/A)')
     *   ' --> ERROR (OPENDX) OLD FILE NOT FOUND',
     *   '                    FILE NAME   =',NAME,
     *   '                    FILE NUMBER =',LUDX,
     *   ' --- or wrong record length on old file.'
      CALL QTRACE(LUPRI)
      CALL QUIT('ERROR (OPENDX) FILE NOT FOUND')
C
C end of OPENDX
C
      END
C  /* Deck finddx */
      LOGICAL FUNCTION FINDDX(LU,I,LEN,IVEC)
C
C 27-Jun-1985 Hans Jorgen Aa. Jensen
C
C For direct access find record,
C when LEN may be greater than maximum record length.
C
      INTEGER IVEC(LEN)
      READ (LU, REC=I, IOSTAT=IOS) IVEC
      IF (IOS .NE. 0) GO TO 900
      FINDDX = .TRUE.
      RETURN
C
  900 CONTINUE
      FINDDX = .FALSE.
      RETURN
      END
C  /* Deck readdx */
      SUBROUTINE READDX(LU,I,LEN,IVEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
C For direct access read
C
#include "implicit.h"
#include "priunit.h"
#include "2gbdef.h"
#include "dummy.h"
      LOGICAL OLDDX
      CHARACTER*80 FNNAME, FNNM2
      INTEGER IVEC(LEN)
#include "inftap.h"
#include "chrnos.h"
C
      IF (LEN .EQ. 0) RETURN
#if defined (VAR_SPLITFILES)
      INQUIRE(UNIT=LU,RECL=LRECL)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90) || defined (SYS_LINUX)
      LRECL = LRECL/8
#else
#if !defined (SYS_NEC) && !defined (SYS_IRIX)
      LRECL = LRECL/4
#endif
#endif
      IRECPF = I2GB/LRECL
C
C     Error branch for very large record lengths
C
      IF (IRECPF .EQ. 0) THEN
         WRITE (LUPRI,'(/A,I14,A)')
     &        ' ERROR: The requested record length of ', LRECL,
     &        ' is larger than 2Gb'
         WRITE (LUPRI,'(/A)') ' The program cannot write even a '//
     &        'single record to file!'
         CALL QUIT('Record length of file too long for this computer')
      END IF
C
      IF (I .GT. IRECPF) THEN
C
C     Ooops, this record is not in this file
C
         LFILNM = I/IRECPF - 1
         IF (LFILNM .GT. 9) THEN
            WRITE (LUPRI,'(/A)') ' DALTON needs to read from a file '//
     &           ' split more than 11 times.',
     &           ' This is currently not supported'
            CALL QUIT('Too many splittings of a file')
         END IF
         IRECT = MOD(I,IRECPF)
         IF (IRECT .EQ. 0) IRECT = IRECPF
         INQUIRE(UNIT=LU,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         FNNM2 = FNNAME(1:LN)//'-'//CHRNOS(LFILNM)
         CALL GPOPEN(LU,FNNM2(1:(LN+2)),'UNKNOWN','DIRECT',' ',LRECL,
     &               OLDDX)
      ELSE
         IRECT = I
      END IF
#else
      IRECT = I
#endif
      READ (LU, REC = IRECT) IVEC
#if defined (VAR_SPLITFILES)
      IF (I .GT. IRECPF) THEN
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         CALL GPOPEN(LU,FNNAME(1:LN),'UNKNOWN','DIRECT',' ',LRECL,OLDDX)
      END IF
#endif
      RETURN
      END
C  /* Deck writdx */
      SUBROUTINE WRITDX(LU,I,LEN,IVEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
#include "implicit.h"
#include "dummy.h"
#include "2gbdef.h"
#include "priunit.h"
      LOGICAL OLDDX
      CHARACTER*80 FNNAME, FNNM2
      INTEGER IVEC(LEN)
#include "inftap.h"
#include "chrnos.h"
C
#if defined (VAR_SPLITFILES)
      INQUIRE(UNIT=LU,RECL=LRECL)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90) || defined (SYS_LINUX)
      LRECL = LRECL/8
#else
#if !defined (SYS_NEC) && !defined (SYS_IRIX)
      LRECL = LRECL/4
#endif
#endif
      IRECPF = I2GB/LRECL
C
C     Error branch for very large record lengths
C
      IF (IRECPF .EQ. 0) THEN
         WRITE (LUPRI,'(/A,I14,A)')
     &        ' ERROR: The requested record length of ', LRECL,
     &        ' is larger than 2Gb'
         WRITE (LUPRI,'(/A)') ' The program cannot write even a '//
     &        'single record to file!'
         CALL QUIT('Record length of file too long for this computer')
      END IF
C
      IF (I .GT. IRECPF) THEN
C
C     Ooops, this record will not fit in this file
C
         LFILNM = I/IRECPF - 1
         IF (LFILNM .GT. 9) THEN
            WRITE (LUPRI,'(/A)') ' DALTON need to split a file more '//
     &           ' than 11 times.', ' This is currently not supported'
            CALL QUIT('Too many splittings of a file')
         END IF
         IRECT = MOD(I,IRECPF)
         IF (IRECT .EQ. 0) IRECT = IRECPF
         INQUIRE(UNIT=LU,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         FNNM2 = FNNAME(1:LN)//'-'//CHRNOS(LFILNM)
         LN2 = LN + 2
         CALL GPOPEN(LU,FNNM2(1:LN2),'UNKNOWN','DIRECT',' ',LRECL,
     &               OLDDX)
      ELSE
         IRECT = I
      END IF
#else
      IRECT = I
#endif
      WRITE (LU, REC = IRECT) IVEC
#if defined (VAR_SPLITFILES)
      IF (I .GT. IRECPF) THEN
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         CALL GPOPEN(LU,FNNAME(1:LN),'UNKNOWN','DIRECT',' ',LRECL,OLDDX)
      END IF
#endif
      RETURN
      END
C  /* Deck aolab4 */
      SUBROUTINE AOLAB4(IINDPK,NMAX,NIBUF,NBITS,IINDX4,N)
C
C     Written by Henrik Koch 22-Nov-1991 (as AOLABE)
C     Generalized 25-Oct-1993 hjaaj
C                 31-Jul-1996 hjaaj (new NIBUF parameter)
C
C     Unpack 4 integer indices packed in NIBUF *4 integer(s)
C     or one *8 integer.
C
#include "implicit.h"
      DIMENSION IINDPK(*), IINDX4(4,*)
#include "ibtdef.h"
C
      N = IINDPK(NIBUF*NMAX+1)
C
      IF (NIBUF .EQ. 1) THEN
        IF (NBITS .EQ. 8) THEN
          DO I = 1,N
            LABEL = IINDPK(I)
            IINDX4(1,I) = IBTAND(IBTSHR(LABEL,24),IBT08)
            IINDX4(2,I) = IBTAND(IBTSHR(LABEL,16),IBT08)
            IINDX4(3,I) = IBTAND(IBTSHR(LABEL, 8),IBT08)
            IINDX4(4,I) = IBTAND(       LABEL,    IBT08)
          END DO
        ELSE IF (NBITS .EQ. 16) THEN
#ifdef INT_STAR8
          DO I = 1,N
            LABEL = IINDPK(I)
            IINDX4(1,I) = IBTAND(IBTSHR(LABEL,48),IBT16)
            IINDX4(2,I) = IBTAND(IBTSHR(LABEL,32),IBT16)
            IINDX4(3,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
            IINDX4(4,I) = IBTAND(       LABEL,    IBT16)
          END DO
#else
          CALL QUIT('AOLAB4 error; NBITS.eq.16 for INTEGER*4')
#endif
        ELSE
          CALL QUIT('AOLAB4 error; NBITS neither 8 nor 16')
        END IF
      ELSE
        DO I = 1,N
          LABEL = IINDPK(I)
          IINDX4(1,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
          IINDX4(2,I) = IBTAND(       LABEL    ,IBT16)
          LABEL = IINDPK(NMAX+I)
          IINDX4(3,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
          IINDX4(4,I) = IBTAND(       LABEL    ,IBT16)
        END DO
      END IF
C
      RETURN
      END
C  /* Deck GPOPEN */
      SUBROUTINE GPOPEN(LUNIT,FILEIN,STATIN,ACCEIN,FORMIN,NELEM,OLDDX)
C
C 23-Feb-2000 K.Ruud
C
C Purpose:
C   General purpose routine for opening files in the Dalton program.
C   The routine will dynamically allocate unit numbers that will
C   become available again when the file is closed using GPCLOSE
C
C   Direct access files are provided with a UNIT-number, but the rest
C   of the file opening process is taken care of by a call to OPENDX.
C
C   These files are strongly machine dependent, although care has been
C   taken to avoid using unit numbers that are illegal or reserved on a
C   particular architecture. However, the use of GPOPEN and GPCLSE should
C   remove much of the machine dependence in the rest of the Dalton
C   program
C
C Input:
C  LUNIT    Suggested unit number (OPTIONAL, but is mandatory if an unnamed
C           file is reopened after having been closed with STATUS='KEEP')
C           Otherwise it is recommended to not assign this.
C  FILEIN   Suggested name for the file (OPTIONAL, but strongly recommended)
C  STATIN   Suggested status of the file (OPTIONAL and maybe not recommended)
C  ACCEIN   Access method. 'DIRECT' or 'SEQUENTIAL'. Default is sequential.
C  FORMIN   Formatted or unformatted file format. Default is 'UNFORMATTED'
C
C Direct access-specific input:
C  NELEM    number of integer words per logical record
C  NREC     number of logical records
C
C
C Output:
C  LUNIT    Assigned file unit number
C
C  Direct access-specific output (through 
C  LRDX     physical record length (in integers)
C  NBDX     number of physical records per logical record
C  OLDDX    logical, true if old LUDX file was opened
C
C
#include "implicit.h"
#include "priunit.h"
#include "chrnos.h"
      CHARACTER*80  FILENM, STATUF, ACCESF, FORMF
      CHARACTER*(*) FILEIN, STATIN, ACCEIN, FORMIN
      INTEGER       FILELN, STATLN, ACCELN, FORMLN
      LOGICAL OLDDX, FEXST, FIRST, access_is_direct
      COMMON /UNITAR/ IUNTAB(1:99)
      DATA FIRST /.TRUE./
      SAVE FIRST
C
C     Initialization
C     
      IF (FIRST) THEN
         CALL IZERO(IUNTAB,99)
         FIRST = .FALSE.
      END IF
C
C     We first deal with the unit number
C
      FILELN = LNBLNK(FILEIN)
      STATLN = LNBLNK(STATIN)
      ACCELN = LNBLNK(ACCEIN)
      FORMLN = LNBLNK(FORMIN)
      FILENM(1:FILELN) = FILEIN(1:FILELN)
      STATUF(1:STATLN) = STATIN(1:STATLN)
      ACCESF(1:ACCELN) = ACCEIN(1:ACCELN)
      FORMF (1:FORMLN) = FORMIN(1:FORMLN)
      IF ((LUNIT .LT. 1) .OR. (LUNIT .GT. 99)) THEN
C
C     Unit number left unassigned, we get to decide!
C
         IUN = 0
 10      CONTINUE
         IUN = IUN + 1
         IF ((IUN .EQ. 5) .OR. (IUN .EQ. 6))
     &        GOTO 10
         IF (IUNTAB(IUN) .NE. 0) GOTO 10
         LUNIT = IUN
         IF (LUNIT .GT. 99) GOTO 9001
      ELSE
C
C     The user has requested a specific unit number. We don't 
C     quite trust the user, so we check that
C     1) It is not unit 5, or 6
C     2) The file either has been closed with status='KEEP' or is
C        not currently in use
C
         IF ((LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9002
         IF (IUNTAB(LUNIT) .EQ. 1) GOTO 9003
      END IF
C
C     We've got a file number now, get rid of the direct access case
C

      access_is_direct = .false.
      if (acceln > 5) then
         if (accesf(1:6) == 'DIRECT') then
            access_is_direct = .true.
         end if
      end if

      IF (access_is_direct) THEN
         IF (FILELN .EQ. 0) THEN
            FILELN = 13
            FILENM(1:FILELN) = 'GPOPEN_UNIT'//CHRNOS(LUNIT/10)//
     &                    CHRNOS(MOD(LUNIT,10))
         END IF
         IF (STATLN .LE. 1) THEN
            STATLN = 7
            STATUF(1:7) = 'UNKNOWN'
         END IF
         CALL OPENDX(LUNIT,FILENM(1:FILELN),NELEM,STATUF(1:STATLN),
     &               OLDDX)
      ELSE
         IF (STATLN .LE. 1) THEN
            STATLN = 7
            STATUF(1:7) = 'UNKNOWN'
         END IF
         IF (STATUF(1:3) .EQ. 'OLD' .AND. IUNTAB(LUNIT) .EQ. 0) THEN
C
C     This better be a file with a name, and it better exist
C
            IF (FILELN .EQ. 0) GOTO 9005
            INQUIRE(FILE=FILENM(1:FILELN),EXIST=FEXST,IOSTAT=IOS,
     &              ERR=9006)
            IF (.NOT. FEXST) GOTO 9007
            IF (FORMLN .EQ. 0) THEN
               FORMLN = 9
               FORMF(1:FORMLN)='FORMATTED'
            END IF
            OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),STATUS='OLD',
     &           ACCESS='SEQUENTIAL',FORM=FORMF(1:FORMLN),ERR=9004,
     &           IOSTAT=IOS)
         ELSE
            IF (STATUF(1:3) .EQ. 'NEW' .AND. IUNTAB(LUNIT) .NE. 0)
     &           GOTO 9006
            IF (FILELN .EQ. 0) THEN
               FILELN = 13
               FILENM(1:FILELN) = 'GPOPEN_UNIT'//CHRNOS(LUNIT/10)//
     &                               CHRNOS(MOD(LUNIT,10))
            END IF
            IF (FORMLN .EQ. 0) THEN
               FORMLN = 11
               FORMF(1:11)='UNFORMATTED'
            END IF
            IF (STATUF(1:7) .EQ. 'SCRATCH') THEN
               OPEN(UNIT=LUNIT,STATUS='SCRATCH',
     &              ACCESS='SEQUENTIAL',FORM=FORMF(1:FORMLN),
     &              ERR=9004,IOSTAT=IOS)
            ELSE IF (ACCESF(1:8) == 'READONLY') THEN
               OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),
     &              STATUS=STATUF(1:STATLN),ACCESS='SEQUENTIAL',
     &              ACTION='READ',
     &              FORM=FORMF(1:FORMLN),ERR=9004,IOSTAT=IOS)
            ELSE
               OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),
     &              STATUS=STATUF(1:STATLN),ACCESS='SEQUENTIAL',
     &              FORM=FORMF(1:FORMLN),ERR=9004,IOSTAT=IOS)
            END IF
         END IF
      END IF
      IUNTAB(LUNIT) = 1
      RETURN
C     
C error branches
C
 9001 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPOPEN) NO MORE AVAILABLE FILENUMBERS! ',
     &     ' *** THIS CALCULATION EITHER NEEDS TOO MANY SIMULTANEOUS '//
     &     'FILES OR', ' *** SOMEBODY HAS FORGOT TO CLOSE FILES IN'//
     &     'THE SOURCE CODE'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      GOTO 9100
C
 9002 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPOPEN) TRYING TO OPEN AN ILLEGAL FILE NUMBER',
     &     ' *** SOMEBODY HAS TRIED TO OPEN UNITS 0, 3, 5, OR 6',
     &     ' *** THE PROGRAM DOES NOT ALLOW THE USE OF THESE RESERVED'//
     &     ' UNIT NUMBERS'
      WRITE (6,'(//A,I5)') ' Offending UNIT number: ',LUNIT
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      GOTO 9100
C
 9003 CONTINUE
      WRITE (6,'(//A/A)')
     &     ' *** ERROR (GPOPEN) TRYING TO OPEN A FILE ALREADY IN USE',
     &     ' *** SOMEBODY IS TRYING TO USE A FILENUMBER THAT IS '//
     &     'ALREADY IN USE'
      WRITE (6,'(//A,I5)') ' Offending UNIT number: ',LUNIT
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      GOTO 9100
C
 9004 CONTINUE
#ifdef VAR_GFORTRAN
C
      call perror(' GPOPEN calling perror:')
#endif
C
      WRITE (6,'(//A,I3/A,I7)')
     &     ' *** ERROR (GPOPEN) UPON TRYING TO OPEN FILE ON UNIT',
     &     LUNIT, ' *** ERROR CODE RETURNED ',IOS
      GOTO 9100
C
 9005 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPOPEN) TRYING TO OPEN A NON-EXISTING OLD FILE',
     &     ' *** A FILE HAS BEEN SPECIFIED TO BE OLD, BUT THE TABLE',
     &     ' *** ENTRY INDICATE THAT IT DOES NOT EXIST'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      GOTO 9100
C
 9006 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPOPEN) TRYING TO OPEN AN EXISTING NEW FILE',
     &     ' *** A FILE HAS BEEN SPECIFIED TO BE NEW, BUT THE TABLE',
     &     ' *** ENTRY INDICATE THAT IT ALREADY EXISTS'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      GOTO 9100
C
 9007 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPOPEN) TRYING TO OPEN A NON-EXISTING OLD FILE',
     &     ' *** A FILE HAS BEEN SPECIFIED TO BE OLD, BUT THE FILE',
     &     ' *** DOES NOT EXIST'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      GOTO 9100
 9100 CONTINUE
      write(6,'(/A/A)') 'Parameters specified to GPOPEN:',
     &                  '==============================='
      write(6,*) 'File unit ',LUNIT
      write(6,*) 'File name ',FILEIN,' -> ',FILENM(1:FILELN)
      write(6,*) 'Status    ',STATIN,' -> ',STATUF(1:STATLN)
      write(6,*) 'Access    ',ACCEIN,' -> ',ACCESF(1:ACCELN)
      write(6,*) 'Format    ',FORMIN,' -> ',FORMF(1:FORMLN)
      write(6,*) 'NELEM     ',NELEM
      write(6,*) 'OLDDX     ',OLDDX
      inquire(unit=lunit,name=filenm)
      write(6,*) 'File name from inquire for this file unit: ',FILENM
      CALL QTRACE(6)
      CALL QUIT('*** (GPOPEN) ERROR TRYING TO OPEN A FILE')
C
C end of GPOPEN
C
      END
C  /* Deck gpclse */
      SUBROUTINE GPCLOSE(LUNIT,DISP)
C
C 23-Feb-2000 K.Ruud
C
C Purpose:
C   General purpose routine for closing files in the Dalton program.
C   The routine will ensure that files that are closed and that will not
C   be reopened in a later part of the program will make their unit 
C   number available for reuse.
C
C   These files are strongly machine dependent, although care has been 
C   taken to avoid using unit numbers that are illegal or reserved on a
C   particular architecture. However, the use of GPOPEN and GPCLSE should
C   remove much of the machine dependence in the rest of the Dalton 
C   program
C
C Input:
C  LUNIT    Mandatory
C  STATUS   Indicate whether the file should be removed or kept
C           ('KEEP' or 'DELETE'). Note that the default is that the file
C           is to be deleted
C
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) DISP
      CHARACTER*16 FNNAME, FNNM2
      CHARACTER*3 IODIR
      LOGICAL FEXIST, OPND, OLDDX
      COMMON /UNITAR/ IUNTAB(1:99)
#include "chrnos.h"
C
C     We first deal with the unit number
C
      IF ((LUNIT .LT. 1) .OR. (LUNIT .GT. 99) .OR.
     &    (LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9001
C
C     Check that the file actually have been opened
C
      IF (IUNTAB(LUNIT) .EQ. 0) GOTO 9002
C
      INQUIRE(UNIT=LUNIT,EXIST=FEXIST,OPENED=OPND)
      IF (.NOT. FEXIST .OR. .NOT. OPND) WRITE (LUPRI,'(/A/A,I3)')
     &     ' WARNING: Tried to close a non-existent or already '//
     &     'closed file', ' Unit number was: ',LUNIT
      CALL FLSHFO(LUPRI)
      IF (DISP .EQ. 'KEEP') THEN
         CLOSE (LUNIT,STATUS='KEEP',ERR=9003)
         GOTO 30
      ELSE
#if defined (VAR_SPLITFILES)
         INQUIRE(UNIT=LUNIT,DIRECT=IODIR)
         INQUIRE(UNIT=LUNIT,NAME=FNNAME)
         IF (IODIR .EQ. 'YES') THEN
            INQUIRE(UNIT=LUNIT,RECL=LRECL)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90) || defined (SYS_LINUX)
            LRECL = LRECL/8
#endif
#if !defined (SYS_NEC) && !defined (SYS_IRIX)
            LRECL = LRECL/4
#endif
         END IF
         LN = 1
 15      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 15
         END IF
         LN = LN - 1
         IF (IODIR .NE. 'YES') THEN
            IF (FNNAME((LN-1):(LN-1)) .EQ. '-') THEN
               CLOSE(LUNIT,STATUS='KEEP')
               LN = LN - 2
               OPEN (UNIT=LUNIT,FILE=FNNAME(1:LN),STATUS='OLD')
            END IF
         END IF
         CLOSE (LUNIT,STATUS='DELETE',ERR=20)
         IUNIT = 0
 16      CONTINUE
         FNNM2 = FNNAME(1:LN)//'-'//CHRNOS(IUNIT)
         LN2 = LN + 2
         INQUIRE(FILE=FNNM2(1:LN2),EXIST=FEXIST)
         IF (FEXIST) THEN
            IF (IODIR .EQ. 'YES') THEN
               CALL OPENDX(LUNIT,FNNM2(1:LN2),LRECL,'OLD',OLDDX)
            ELSE
               OPEN(UNIT=LUNIT,FILE=FNNM2(1:LN2),STATUS='OLD')
            END IF
            CLOSE(LUNIT,STATUS='DELETE')
            IUNIT = IUNIT + 1
            GOTO 16
         END IF
#else
         CLOSE (LUNIT,STATUS='DELETE',ERR=20)
#endif
         GOTO 30
      END IF
C
C     We might get an error for trying to close a SCRATCH file with STATUS
C     We just close it
C     
 20   CLOSE (LUNIT)
C
C     We release that saved unit number by resetting it to 0
C
 30   CONTINUE
      IUNTAB(LUNIT) = 0
      LUNIT = 0
      RETURN
C     
C error branches
C
 9001 CONTINUE
      WRITE (6,'(//A/A/A)')
     &    ' *** ERROR (GPCLOSE) TRYING TO CLOSE AN ILLEGAL FILE NUMBER',
     &    ' *** SOMEBODY HAS TRIED TO CLOSE UNITS 0, 3, 5, OR 6',
     &    ' *** THE PROGRAM DOES NOT ALLOW THE USE OF THESE RESERVED'//
     &    ' UNIT NUMBERS'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPCLOSE) ILLEGAL FILE NUMBER REQUESTED')
C
 9002 CONTINUE
      WRITE (6,'(//A/A)')
     &     ' *** ERROR (GPCLOSE) TRYING TO CLOSE A FILE NOT IN USE',
     &     ' *** SOMEBODY IS TRYING TO USE A FILENUMBER THAT HAS '//
     &     'NOT BEEN USED'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPCLOSE) TRYING TO CLOSE A FILE NOT '//
     &          'IN USE')
C
 9003 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPCLOSE) TRYING TO KEEP A SCRATCH FILE',
     &     ' *** A FILE HAS BEEN INDICATED TO BE CLOSED AND KEPT,',
     &     ' *** BUT IT APPEARS THE FILE IS A SCRATCH FILE'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPCLOSE) TRYING TO CLOSE A SCRATCH FILE ')
C
C end of GPCLOSE
C
      END
#if ! defined (PRG_DIRAC)
C  /* Deck gpioin */
      SUBROUTINE GPIOIN
C
C 1-Mar-2000 K.Ruud
C
#include "priunit.h"
#include "inftap.h"
C
      LUCMD  = 0
      LUPRI  = 0
C
      LUCME  = 0
      LUMOL  = 0
      LUPROP = 0
      LUSOL  = 0
      LUINTA = 0
      LUONEL = 0
      LUSUPM = 0
      LUTLM  = 0
      LUDA1  = 0
      LUITMP = 0
      LU2DER = 0
      LUDASP = 0
      LUAAT  = 0
      LURDR  = 0
      LURDI  = 0
      LUGDR  = 0
      LUGDI  = 0
      LUGDT  = 0
      LURDT  = 0
      LUDFCK = 0
      LUSFDA = 0
      LUFDC  = 0
      LUWLK  = 0
      LUPAO  = 0
      LUPAS  = 0
      LUNR1  = 0
      LUNR3  = 0
      LUNR5  = 0
      LUCMD  = 0
      LUINTR = 0
      LUAORC = 0
      LUMOLDEN=0
C
C     /RSPTAP/
C
      LUAHSO = 0
      LUCRV1 = 0
      LUCRV2 = 0
      LUXYVE = 0
      LUCRVE = 0
      LUQRVE = 0
      LURSP3 = 0
      LURSP4 = 0
      LURSP5 = 0
      LUMHSO = 0
C
C     /SIRTAP/
C
      LUINTM = 0
      LUIT1  = 0
      LUINF = 0
      LUSIFC = 0
C
      FNSOL  = 'AOSOLINT'
      ABARDR = 'ABACUS.RD'
      ABARDI = 'ABACUS.RDI'
      ABAGDR = 'ABACUS.GD'
      ABAGDI = 'ABACUS.GDI'
      ABAGDT = 'ABACUS.GDT'
      ABARDT = 'ABACUS.RDT'
      ABADFK = 'ABACUS.DFK'
      ABASF  = 'ABACUS.SF'
      ABATLM = 'ABACUS.TLM'
      ABAWLK = 'ABACUS.WLK'
      ABANR1 = 'ABAENR.RST'
      ABANR3 = 'ABAENR.BVC'
      ABANR5 = 'ABAENR.SVC'
      FNINTM = 'MOTWOINT'
      FNSUPM = 'AOSUPINT'
      FNONEL = 'AOONEINT'
      FNSIFC = 'SIRIFC'
      LBSIFC = 'SIR IPH '
C
      RETURN
      END
C /* Deck daintf */
      SUBROUTINE DAINTF(LUNIT,NAME)
C
C     Interface for assigning a unit number to a unit before opening it in
C     DANAME. Needed to ensure compatibility with use of GPOPEN/GPCLOSE
C     K.Ruud, March 28 1999
C
#include "implicit.h"
      CHARACTER NAME*(*)
      COMMON /UNITAR/ IUNTAB(1:99)
C      
      IF (LUNIT .EQ. 0) THEN
         IUN = 0
 10      CONTINUE
         IUN = IUN + 1
         IF ((IUN .EQ. 5) .OR. (IUN .EQ. 6))
     &        GOTO 10
         IF (IUNTAB(IUN) .NE. 0) GOTO 10
         LUNIT = IUN
         IF (LUNIT .GT. 99) GOTO 9001
      END IF
      RETURN
C     
C error branch
C
 9001 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (DAINTF) NO MORE AVAILABLE FILENUMBERS! ',
     &     ' *** THIS CALCULATION EITHER NEEDS TOO MANY SIMULTANEOUS '//
     &     'FILES OR', ' *** SOMEBODY HAS FORGOT TO CLOSE FILES IN'//
     &     'THE SOURCE CODE'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (DAINTF) NO MORE FILE NUMBERS')
C
      END
#endif
C /* Deck UPCASE */
      SUBROUTINE UPCASE(WORD)
C
C     Routine that uppercases the word provided. Allows the input easily to
C     be made case insensitive. F90 routine supplied by Anthony Stone and 
C     77-ified by K.Ruud
C
C     Revision March 2010 hjaaj: remove any inline comment
C
#include "implicit.h"
      CHARACTER*(*) word
      CHARACTER*26  UPPER, LOWER
C
      UPPER = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      LOWER = "abcdefghijklmnopqrstuvwxyz"
C
C     Remove in-line comment /hjaaj March 2010
      K = INDEX(WORD,'!')
      IF (K .EQ. 0) K = INDEX(WORD,'#')
      LEN_WORD = LNBLNK(WORD)
      IF ( K .GT. 1 .AND. LEN_WORD .GT. 1) THEN
C        ... NB! do not remove the "!" or "#" if in column 1
C            because then the tests in the *INP routines do not work!
         WORD(K:LEN_WORD) = ' '
         LEN_WORD = LNBLNK(WORD)
      END IF
C
      DO I = 1, LEN_WORD
         K = INDEX(LOWER,WORD(I:I))
         IF (K .NE. 0) WORD(I:I) = UPPER(K:K)
         IF (ICHAR(WORD(I:I)) .EQ. 9) WORD(I:I) = ' ' ! change tab to blank
      END DO
C
      RETURN
      END
! -- end of gptrygve.F --
