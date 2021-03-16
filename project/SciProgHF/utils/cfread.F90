      PROGRAM CFREAD
!***********************************************************************
!
!     Get basis information from a coefficient file
!
!     Written by T. Saue Oct 7 2010
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (LUINP = 1)
      CHARACTER TEXT*74,FMT*6,MXFORM*6
      CHARACTER*8 B(4), STARS
      PARAMETER (STARS = '********')
      DIMENSION IDIM(3,2)      
      INTEGER   NARG,IARG,JARG
      CHARACTER*80 ARGC
      LOGICAL FNDLAB
      NARG = command_argument_count()
      IF (NARG .EQ. 0) THEN
         CALL get_command_argument(0,ARGC)
         JARG = LEN_TRIM(ARGC)
         WRITE(6,'(/3A/)') 'Usage: ',ARGC(1:JARG),' file1 file2 ...'
         STOP ' '
      END IF
      DO IARG = 1,NARG
        CALL get_command_argument(IARG,ARGC)
        JARG = LEN_TRIM(ARGC)
        WRITE (LUPRI, '(/2A/)')                                         &
      ' CFREAD; Information read from file : ',ARGC(1:JARG)
!
        OPEN(LUINP,FILE=ARGC,STATUS='UNKNOWN',FORM='UNFORMATTED')
        REWIND LUINP
        IREC = 0
        IERR = 0
   1    READ (LUINP,END=3,ERR=2) B
        IREC = IREC + 1
        IF (B(1) .NE. STARS) GO TO 1
        WRITE (LUPRI, '(I5,3X,4(2X,A8))')  IREC, B
        IF(B(4).EQ.'INFO    ') THEN
          READ (LUINP,END=10,ERR=20) TEXT,NSYM,NZBUF,                         &
          ((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG
          WRITE(LUPRI,'(/10X,2A)') ' - Heading :',TEXT
          FMT = MXFORM(TOTERG,20)
          WRITE(LUPRI,'(10X,A,'//FMT//')') ' - Total energy: ',TOTERG
          WRITE(LUPRI,'(10X,A,2I6)')                                          &
          ' - Fermion ircops          : ',(I,I=1,NSYM)
          WRITE(LUPRI,'(10X,A,2A )')                                          &
          '-----------------------------',('------',I=1,NSYM)
          WRITE(LUPRI,'(10X,A,2I6)')                                          &
          ' - Negative-energy orbitals: ',(IDIM(1,I),I=1,NSYM)
          WRITE(LUPRI,'(10X,A,2I6)')                                          &
          ' - Positive-energy orbitals: ',(IDIM(2,I),I=1,NSYM)
          WRITE(LUPRI,'(10X,A,2I6)')                                          &
               ' - AO-basis functions      : ',(IDIM(3,I),I=1,NSYM)
          WRITE(LUPRI,'()')
        ENDIF
        GO TO 1
!
   2    CONTINUE
        IREC = IREC + 1
        IERR = IERR + 1
        WRITE (LUPRI, '(/A,I5/)') ' ERROR READING RECORD NO.',IREC
        REWIND LUINP
        DO 102 I = 1,IREC
          READ (LUINP) J
  102    CONTINUE
        IF (IERR .LE. 2) GO TO 1
  202   CONTINUE
           READ (LUINP,END=3) J
           IREC = IREC + 1
        GO TO 202
!
  3  CONTINUE
         WRITE (LUPRI,'(/I10,A)') IREC,                                 &
       ' records read before EOF on file.'
        GOTO 30
 10     CONTINUE
        CALL QUIT('REACMO: END OF FILE reading TEXT')
 20     CONTINUE
        CALL QUIT('REACMO: ERROR reading TEXT')
 30     CONTINUE
        CLOSE(LUINP)
      ENDDO
      END

