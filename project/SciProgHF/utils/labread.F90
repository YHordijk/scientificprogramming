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

!
      PROGRAM READLABL
!
! Written by Hans Jorgen Aa. Jensen, Nov. 1983
! CRAY version 10. Oct. 1986
! Alliant version 30. Dec. 1987
! Linux/IBM, +updates, Joern Thyssen, Nov. 1999
!
!#include "implicit.h"
#include "../include/implicit.h"
!TROND      DOUBLE PRECISION TIME, SECOND
      CHARACTER*8 B(4), STARS
      PARAMETER (STARS = '********')
      PARAMETER (LUINP = 1, LUPRI = 6)
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
       ' MOLECULE labels found on file :',ARGC(1:JARG)
!
         OPEN(LUINP,FILE=ARGC,STATUS='UNKNOWN',FORM='UNFORMATTED')
         REWIND LUINP
         IREC = 0
         IERR = 0
    1    READ (LUINP,END=3,ERR=2) B
            IREC = IREC + 1
         IF (B(1) .NE. STARS) GO TO 1
            WRITE (LUPRI, '(5X,I5,3X,4(2X,A8))')  IREC, B
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
  202    CONTINUE
            READ (LUINP,END=3) J
            IREC = IREC + 1
         GO TO 202
!
    3    CONTINUE
!      TIME = SECOND() - TIME
         WRITE (LUPRI,'(/I10,A)') IREC,                                 &
       ' records read before EOF on file.'
!        WRITE (LUPRI,'(F10.2,A)') TIME,' CPU seconds used.'
         CLOSE(LUINP)
      END DO
      STOP ' '
!
      END
