!      Copyright (c) 2018 by the authors of DIRAC.
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

      PROGRAM RSPREAD
!
! Written by Trond Saue July 4 2001
! based on program LABREAD
!
#include "implicit.h"
#include "priunit.h"
      PARAMETER (LUINP = 1)
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2,LBUF*16
      INTEGER   NARG,IARG,JARG
      CHARACTER*80 ARGC
!
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
         WRITE (LUPRI, '(/2A/)')                                        &
       ' Solution vectors found on file : ',ARGC(1:JARG)
!
         OPEN(LUINP,FILE=ARGC,STATUS='UNKNOWN',FORM='UNFORMATTED')
         REWIND LUINP
         IREC = 0
         IERR = 0
         IVEC = 0
         LBUF = '                '
    1    CONTINUE
         READ (LUINP,END=3,ERR=2) LAB1
         IF(LAB1.EQ.'END_OF_THIS_FILE') GOTO 10
         BACKSPACE LUINP
         READ (LUINP,END=3,ERR=2) LAB1,LAB2,TYP1,TYP2,                  &
           FREQ1,FREQ2,JSYMOP1,JSYMOP2,                                 &
           JTIMOP1,JTIMOP2,RNORM,LEN,INTFLG,                            &
           ERGRSP,NBAS,NORB
         IVEC = IVEC + 1
         IREC = IREC + 1
         IF(LAB1.NE.LBUF) THEN
           WRITE(LUPRI,'(A,A16,2(A,I3))')                               &
       '** Solution vector : ',LAB1,                                    &
            ' Irrep: ',JSYMOP1,' Trev: ',JTIMOP1
           LBUF = LAB1
         ENDIF
         WRITE(LUPRI,'(I3,A,A2,A,F12.6,A,E10.2,A,I10)')                 &
      IVEC, ' Type : ',TYP1,' Freq.: ',FREQ1,                           &
      ' Rnorm: ',RNORM,' Length:', LEN
         READ(LUINP,END=3,ERR=2)
         READ(LUINP,END=3,ERR=2)
         IREC = IREC + 1
         GO TO 1
!
    2    CONTINUE
         IREC = IREC + 1
         IERR = IERR + 1
         WRITE (LUPRI, '(/A,I5/)') ' ERROR READING RECORD NO.',IREC
         GOTO 10
!
    3    CONTINUE
         WRITE (LUPRI,'(/I10,A)') IREC,                                 &
       ' records read before EOF on file.'
         GOTO 10
 10      CONTINUE
         CLOSE(LUINP)
      ENDDO
!
      END
