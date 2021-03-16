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

C WAIO : Word addressable I/O package, written by Luuk Visscher
C Initial version (adapted from RELCCSD) september 2002
C  
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WAIO_INIT (IW)
C
C     Initialize list of word-addressable files
C     This variable gives the number of records that has been written.
C
#include "waio.h"
C
#if defined (VAR_WAIO_C)
      WRITE (IW,1000) "      C",NWORDS_ON_RECORD*4/1024
#else
      WRITE (IW,1000) "FORTRAN",NWORDS_ON_RECORD*4/1024
#endif
      DO ITAPE = 1, MAX_WAIO_FILES-1
         IWRIT_WAIO(ITAPE) = -1
      ENDDO
C
C     Initialize timers and other statistics
C
      WAIO_RTIME = 0.0
      WAIO_WTIME = 0.0
      WAIO_RCALL = 0.0 
      WAIO_WCALL = 0.0
      WAIO_RNLEN = 0.0
      WAIO_WNLEN = 0.0
C
 1000 FORMAT (/" Initializing word-addressable I/O : the ",A7,
     &        "-interface is used with",I6," KB records")
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WAIO_PRST (IW)
C
C     Print statistics for the Word-addressable I/O
C
#include "waio.h"
      REAL*8 WAIO_RMB,WAIO_WMB,WAIO_RAMB,WAIO_WAMB
C
C     Compute average speed and totals in MB
C
      WAIO_WMB = WAIO_WNLEN / (256 * 1024) 
      WAIO_RMB = WAIO_RNLEN / (256 * 1024) 
      WAIO_WAMB = 0.0
      WAIO_RAMB = 0.0
      IF (WAIO_WTIME.NE.0.0) WAIO_WAMB = WAIO_WMB / WAIO_WTIME
      IF (WAIO_RTIME.NE.0.0) WAIO_RAMB = WAIO_RMB / WAIO_RTIME
C
      WRITE (IW,1000) WAIO_WCALL,WAIO_RCALL,WAIO_WMB,WAIO_RMB,
     &                WAIO_WTIME,WAIO_RTIME,WAIO_WAMB,WAIO_RAMB
C
 1000 FORMAT (/' Statistics for the word-addressable I/O'
     &/' Number of write calls',T40,F20.0
     &/' Number of read calls',T40,F20.0
     &/' Megabytes written',T40,F20.3
     &/' Megabytes read',T40,F20.3
     &/' Seconds spent in reads',T40,F20.3
     &/' Seconds spent in writes',T40,F20.3
     &/' average I/O speed for write (Mb/s)',T40,F20.3
     &/' average I/O speed for read  (Mb/s)',T40,F20.3
     & )
     
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WAIO_OPEN (ITAPE)
C
C     Opens a word-addressable file.
C
      implicit none
      integer itape
C
#include "waio.h"
#if defined (VAR_MPI)
#include "ccpar.inc"
#endif
C
      CHARACTER*6 IDENT
      INTEGER*8 SIZE
      INTEGER LENREC,NINTE,NLOGI,NREAL,IRECLEN
C
      IF (ITAPE.GT.MAX_WAIO_FILES-1.OR.ITAPE.LT.1) THEN
         WRITE(6,*) ' INVALID TAPE NUMBER GIVEN: ITAPE=' ,ITAPE
         CALL QUIT('WAIO_OPEN called with invalid unit number')
      ENDIF
C
C     Make a unique node identifier (use that char(0) marks the 
C     string end in C, and can thus be used to end the file name).
C
#if defined (VAR_MPI)
      IF (NMPROC .LT. 10) THEN  !  MPI ID has one digit
          WRITE(IDENT,'(A1,I1,A1,A3)') '.',MYPROC,CHAR(0),'   '
      ELSEIF (NMPROC .LT. 100) THEN  ! MPI ID has two digits
          WRITE(IDENT,'(A1,I2,A1,A2)') '.',MYPROC,CHAR(0),'  '
      ELSEIF (NMPROC .LT. 1000) THEN  ! MPI ID has three digits
          WRITE(IDENT,'(A1,I3,A1,A1)') '.',MYPROC,CHAR(0),' '
      ELSEIF (NMPROC .LT. 10000) THEN  ! MPI ID has four digits
          WRITE(IDENT,'(A1,I4,A1)') '.',MYPROC,CHAR(0)
      ENDIF
#else
      WRITE(IDENT,'(A1,I1,A1,A3)') '.',0,CHAR(0),'   '
#endif
C
      IF (ITAPE.LT.10) THEN
         WRITE (DN(ITAPE),'(A3,I1,A6)') 'ft0',ITAPE,IDENT
      ELSE
         WRITE (DN(ITAPE),'(A2,I2,A6)') 'ft',ITAPE,IDENT
      ENDIF
C
#if defined (VAR_WAIO_C)
      STATUS = 1
      SIZE = 0
      CALL OPENC(ITAPE,DN(ITAPE),SIZE,STATUS)
#else
      NINTE = NWORDS_ON_RECORD
      NREAL = 0
      NLOGI = 0
      LENREC = IRECLEN (NREAL,NINTE,NLOGI)
      OPEN (ITAPE,FILE=DN(ITAPE),ACCESS='DIRECT',RECL=LENREC)
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WAIO_CLOSE(ITAPE,JCODE)
C
C     JCODE = 3     CLOSE AND SAVE FILE
C     JCODE = 4     CLOSE AND DELETE FILE
C
      implicit none
      integer itape,jcode
C
#include "waio.h"
#if defined (VAR_MPI)
#include "ccpar.inc"
#endif
C
      INTEGER*8 STAT
      IF(JCODE.NE.3.AND.JCODE.NE.4) THEN
         WRITE(*,*) ' INVALID JCODE IN WAIO_CLOSE,  JCODE = ',JCODE
         WRITE(*,*) ' FILE ',ITAPE,'  CLOSED AND SAVED.'
         JCODE = 3
      END IF
C
      IF(ITAPE.EQ.6) CALL QUIT(' YOU CANNOT CLOSE A FILE ON UNIT 6')
C
#if defined (VAR_WAIO_C)
      CALL CLOSEC(ITAPE)
      IF (JCODE.EQ.4) THEN
         STAT = UNLINK(DN(ITAPE))
         IF(STAT.NE.0) THEN
           WRITE(6,*) '  Cannot delete file ',ITAPE,'  STAT = ',STAT
           WRITE(6,*) ' DN = ',DN(ITAPE)
         END IF
      ENDIF
#else
      IF (JCODE.EQ.3) THEN
          CLOSE(ITAPE,STATUS='KEEP')
      ELSEIF (JCODE.EQ.4) THEN
          CLOSE(ITAPE,STATUS='DELETE')
          IWRIT_WAIO(ITAPE) = -1
      ENDIF
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WAIO_SAVE(ITAPE)
C
C     Flush the word-addressable file.
C
      IMPLICIT NONE
      integer itape,jcode
C
      JCODE = 3
      CALL WAIO_CLOSE (ITAPE,JCODE)
      CALL WAIO_OPEN  (ITAPE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WAIO_INTIO (ITAPE,ARRAY,NLEN,OFFSET,DO_READ)
C
C     This routine reads or writes NLEN integer words into an array.
C     The file position is specified by the variable OFFSET that should
C     either be zero (start at begin of file) or positive.
C     This word-addressable I/O is done using either standard c byte 
C     addressable I/O or using FORTRAN direct access files.
C
      implicit none
C
#include "waio.h"
      integer itape,array,nlen,extra,init_waio,irec,irec_first
      integer irec_last,nwords4,offb
      INTEGER BUFF(NWORDS_ON_RECORD)
      INTEGER*8 OFFSET,OFF,NWORDS
      REAL*8 WALL1,WALL2,CPU1,CPU2,WALL
C
      DIMENSION ARRAY(NLEN)
      LOGICAL DO_READ
C
C     Check for trivial cases
C
cdbg  print*,"waio",itape,do_read,offset,nlen
      IF (NLEN.EQ.0) THEN
         RETURN
      ELSEIF (NLEN.LT.0) THEN
         CALL QUIT ("NLEN smaller than zero in WAIO_INTIO")
      ENDIF
C
      IF (OFFSET.LT.0) THEN
         CALL QUIT ("OFFSET smaller than zero in WAIO_INTIO")
      ENDIF
C
C     Initialize timer
C
      CALL GETTIM(CPU1,WALL1)
C
#if defined (VAR_WAIO_C)
      IF (DO_READ) THEN
         CALL RDABSF(ITAPE,ARRAY,NLEN,OFFSET)
      ELSE
         CALL WRABSF(ITAPE,ARRAY,NLEN,OFFSET)
      ENDIF
#else
      INIT_WAIO = 12345678
C
C     Since the amount of words and the offset is usually not
C     an integer multiple of the number of records we first read 
C     the first partially filled record, then all
C     full records, and finally the remaining 
C     partially filled record.
C
C     Variables : NWORDS - Number of words still to be read
C                 OFF    - Offset indicating where the I/O data 
C                          should be placed in the array
C                 IREC   - Current record on direct acces file 
C
      NWORDS = NLEN
      OFF = 1
      IREC = OFFSET / NWORDS_ON_RECORD + 1
C
      OFFB =  MOD(OFFSET,NWORDS_ON_RECORD8)
      IF (OFFB.NE.0) THEN
C
C        Read the first partial record  
C
         EXTRA = NWORDS_ON_RECORD-OFFB
         IF (NWORDS.LT.EXTRA) EXTRA=NWORDS
         OFFB = OFFB + 1
         IF (DO_READ) THEN
            READ (ITAPE,REC=IREC) BUFF
cdbg        print*,"first record",irec,offb,extra,'-->',off 
            CALL ICOPY (EXTRA,BUFF(OFFB),1,ARRAY(OFF),1)
         ELSE
            IF (IWRIT_WAIO(ITAPE).GE.IREC) THEN
C              This record has been written before, copy the
C              old information first before partially overwriting
C              with new information.
               READ (ITAPE,REC=IREC) BUFF
            ELSE
C              This record has NOT been written before, this is not
C              allowed in the current version of WAIO (all records
C              should first be written sequentially before accessing
C              them in a random fashion).
               CALL QUIT ('Trying to read from non-existing record')
            ENDIF
            CALL ICOPY (EXTRA,ARRAY(OFF),1,BUFF(OFFB),1)
cdbg        print*,"first record",irec,offb,extra,'-->',off 
            WRITE (ITAPE,REC=IREC) BUFF
         ENDIF
         IREC = IREC + 1
         NWORDS = NWORDS - EXTRA
         OFF = OFF + EXTRA
      END IF
C
C     Read all full records
C
      IREC_FIRST = IREC
      IREC_LAST = IREC + NWORDS / NWORDS_ON_RECORD - 1
C
      IF (DO_READ) THEN
         DO IREC = IREC_FIRST, IREC_LAST 
cdbg     print*,"middle record",irec,1,NWORDS_ON_RECORD,'-->',off 
            CALL READDI (ITAPE,IREC,NWORDS_ON_RECORD,ARRAY(OFF))
            OFF = OFF + NWORDS_ON_RECORD
            NWORDS = NWORDS - NWORDS_ON_RECORD
         ENDDO
      ELSE
         DO IREC = IREC_FIRST, IREC_LAST 
cdbg     print*,"middle record",irec,1,NWORDS_ON_RECORD,'-->',off 
            CALL WRITDI (ITAPE,IREC,NWORDS_ON_RECORD,ARRAY(OFF))
            OFF = OFF + NWORDS_ON_RECORD
            NWORDS = NWORDS - NWORDS_ON_RECORD
         ENDDO
      ENDIF
C
C     Read the last partial record  
C
      IF (NWORDS.NE.0) THEN
C        Copy NWORDS to an integer*4 variable before passing to the read routine !
         NWORDS4 = NWORDS
         IREC_LAST = IREC_LAST + 1
         IF (DO_READ) THEN
cdbg        print*,"last record",irec_last,1,nwords,'-->',off 
            CALL READDI (ITAPE,IREC_LAST,NWORDS4,ARRAY(OFF))
         ELSE
C           Be careful when writing : the file may already
C           have information on it.
C
            IF (IWRIT_WAIO(ITAPE).GE.IREC_LAST) THEN
C              This record has been written before, copy the
C              old information first before partially overwriting
C              with new information.
               READ (ITAPE,REC=IREC_LAST) BUFF
            ELSE
C              This record has NOT been written before, put an
C              arbitrary number in BUFF (may be useful for debugging)
               CALL ICOPY (NWORDS_ON_RECORD,INIT_WAIO,0,BUFF,1)
            ENDIF
            CALL ICOPY (NWORDS4,ARRAY(OFF),1,BUFF,1)
cdbg        print*,"last record",irec_last,1,nwords,'-->',off 
            WRITE (ITAPE,REC=IREC_LAST) BUFF
         ENDIF
      END IF
C     Update the information about the highest record that is written
      IWRIT_WAIO(ITAPE) = MAX0(IWRIT_WAIO(ITAPE),IREC_LAST)
#endif
C
C     Take statistics
C
      CALL GETTIM(CPU2,WALL2)
      WALL   = WALL2 - WALL1
      IF (DO_READ) THEN
         WAIO_RTIME = WAIO_RTIME + WALL
         WAIO_RCALL = WAIO_RCALL + 1
         WAIO_RNLEN = WAIO_RNLEN + NLEN
      ELSE
         WAIO_WTIME = WAIO_RTIME + WALL
         WAIO_WCALL = WAIO_RCALL + 1
         WAIO_WNLEN = WAIO_WNLEN + NLEN
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WAIO_REALIO (ITAPE,ARRAY,NLEN,OFFSET,DO_READ)
C
C     This routine reads or writes NLEN real words into an array.
C     The file position is specified by the variable OFFSET that should
C     either be zero (start at begin of file) or positive.
C
      implicit none
C
#include "waio.h"
#include "complex.inc"
      INTEGER*8 OFFSET,OFFSETI
C
      REAL*8 ARRAY(NLEN)
      LOGICAL DO_READ
      integer itape,nlen,nleni
C
C     Check for trivial cases
C
      IF (NLEN.EQ.0) THEN
         RETURN
      ELSEIF (NLEN.LT.0) THEN
         CALL QUIT ("NLEN smaller than zero in WAIO_REALIO")
      ENDIF
C
      OFFSETI = IRW * OFFSET
      NLENI = IRW * NLEN
      CALL WAIO_INTIO (ITAPE,ARRAY,NLENI,OFFSETI,DO_READ)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RREAD (ITAPE,ARRAY,NLEN,IREC)
C
C     Used for compatibility with the old code.
C
      implicit none
C
#include "waio.h"
      integer itape,array,nlen,irec
      DIMENSION ARRAY(NLEN)
      LOGICAL DO_READ
      INTEGER*8 OFFSET,TEMP1,TEMP2
C
      TEMP1 = IREC - 1
      TEMP2 = NWORDS_ON_RECORD
      OFFSET = TEMP1 * TEMP2
      DO_READ = .TRUE.
cdbg  print*,"in rread",itape,irec,nlen,offset
      CALL WAIO_INTIO (ITAPE,ARRAY,NLEN,OFFSET,DO_READ)
      RETURN
C
      ENTRY RWRIT (ITAPE,ARRAY,NLEN,IREC)
      TEMP1 = IREC - 1
      TEMP2 = NWORDS_ON_RECORD
      OFFSET = TEMP1 * TEMP2
      DO_READ = .FALSE.
cdbg  print*,"in rwrit",itape,irec,nlen,offset
      CALL WAIO_INTIO (ITAPE,ARRAY,NLEN,OFFSET,DO_READ)
      RETURN
C
      END

