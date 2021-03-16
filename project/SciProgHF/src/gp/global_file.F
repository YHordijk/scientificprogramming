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

C aspg: at this point (2007-05-08) i have not yet setup a way yo 
C       automatically enable VAR_PFS in the make system, so if
C       you need this you have to add "-DVAR_PFS" to the
C
C       DEPEND_FLAGS=
C       CPP_FLAGS=
C       CCP_FLAGS=
C
C       variables in Makefile.config, after it has been generated
C

#if defined (VAR_MPI2) 
C
C Allow for the use of a global file (or parallel file system) in
C moltra/tradr6.F and friends.
C This eliminates the message passing for the matrix transposition
C during the MOLTRA step completely. More importantly, it removes
C the synchronization wait time (most probably due to the limited
C number of buffers).
C
C The current implementation uses (partly asynchronous) MPI I/O
C (part of MPI-2), which is at least available on IBM and SGI Altix
C machines.
C
C ***IMPORTANT*** ***IMPORTANT*** ***IMPORTANT*** ***IMPORTANT***
C
C In order to use this modification you have to define VAR_PFS
C (furthermore, we assume that you already have VAR_MPI2 defined).
C
C ***IMPORTANT*** ***IMPORTANT*** ***IMPORTANT*** ***IMPORTANT***
C
C Note (WML TODO?): it is possible to make a more portable implementation
C using Fortran direct access and multiple files (possibly local) and do
C an intermediate more-or-less perfectly load-balanced message passing
C and re-ordering of local files. (Comparable with the current
C SELECT_LOCAL_FROM_GLOBAL routine.)
C
C Implemented as "enabling work" for the TRIPOD project
C (a DECI / DEISA project).
C Written by Walter.Lioen@sara.nl, 2007.
C
      SUBROUTINE GLOBAL_FILE_OPEN_WRITE
      use interface_to_mpi
      IMPLICIT NONE

#include "bucket_sort.h"
#include "infpar.h"

      INTEGER IERR
      INTEGER GF_INFO
      INTEGER I, J, K

      CALL interface_MPI_INFO_CREATE(GF_INFO)
#if defined (SYS_AIX)
C     ... IBMs GPFS specific file hint - SK Aug 2008
      CALL interface_MPI_INFO_SET(GF_INFO, "IBM_largeblock_io", "true")
#else
      CALL interface_MPI_INFO_SET(GF_INFO, "direct_write", "true")
      CALL interface_MPI_INFO_SET(GF_INFO, "romio_cb_write", "disable")
      CALL interface_MPI_INFO_SET(GF_INFO, "romio_ds_write", "disable")
#endif

      CALL interface_MPI_FILE_OPEN(
     +   global_communicator, GF_NAME,
     +   df_MPI_MODE_WRONLY + df_MPI_MODE_CREATE + df_MPI_MODE_EXCL,
     +   GF_INFO, GF_HANDLE
     +)
      CALL interface_MPI_INFO_FREE(GF_INFO)

      DO K = 0, NUMNOD
         DO J = 0, NUMNOD
            GF_NREC(J,K) = 0
            DO I = 1, N_BUFFER
               GF_BREC(I,J,K) = 0
            ENDDO
         ENDDO
      ENDDO

      WRITE_REQUEST = df_MPI_REQUEST_NULL

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_OPEN_READ
      use interface_to_mpi
      IMPLICIT NONE

#include "bucket_sort.h"
#include "infpar.h"

      INTEGER IERR
      INTEGER GF_INFO
      INTEGER I, NMPROC, NMBUFF

      NMPROC = NUMNOD + 1
      NMBUFF = N_BUFFER * NMPROC

      CALL interface_MPI_INFO_CREATE(GF_INFO)
#if defined (SYS_AIX)
C     ... IBMs GPFS specific file hint - SK Aug 2008
      CALL interface_MPI_INFO_SET(GF_INFO, "IBM_largeblock_io", "true")
#else
      CALL interface_MPI_INFO_SET(GF_INFO, "ind_rd_buffer_size", 
     &                            "16777216")
      CALL interface_MPI_INFO_SET(GF_INFO, "romio_cb_read", "disable")
      CALL interface_MPI_INFO_SET(GF_INFO, "romio_ds_read", "enable")
#endif

      CALL interface_MPI_FILE_OPEN(
     +   global_communicator, GF_NAME,
     +   df_MPI_MODE_RDONLY + df_MPI_MODE_DELETE_ON_CLOSE,
     +   GF_INFO, GF_HANDLE
     +)

      CALL interface_MPI_INFO_FREE(GF_INFO)

      DO I=0,NUMNOD
         CALL interface_MPI_BCAST(GF_NREC(0,I),NMPROC,I,
     +                  global_communicator)
         CALL interface_MPI_BCAST(GF_BREC(1,0,I),NMBUFF,I,
     +                  global_communicator)
      ENDDO

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_CLOSE_WRITE
      use interface_to_mpi
      IMPLICIT NONE

#include "bucket_sort.h"

      INTEGER istat(df_MPI_STATUS_SIZE)

      CALL interface_MPI_WAIT(WRITE_REQUEST,istat)
      CALL interface_MPI_FILE_CLOSE(GF_HANDLE)

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_CLOSE_READ
      use interface_to_mpi
      IMPLICIT NONE

#include "bucket_sort.h"

      CALL interface_MPI_FILE_CLOSE(GF_HANDLE)

      RETURN
      END

      INTEGER FUNCTION GLOBAL_FILE_RECORD(
     +   TARGET_NODE, SOURCE_NODE, LOCAL_RECORD)
      use interface_to_mpi
      IMPLICIT NONE
      INTEGER TARGET_NODE, SOURCE_NODE, LOCAL_RECORD

#include "bucket_sort.h"
#include "infpar.h"
      INTEGER IERR

      IF (TARGET_NODE .LT. 0 .OR. TARGET_NODE .GT. NUMNOD) THEN
         PRINT *, 'ERROR: GLOBAL_FILE_RECORD: TARGET_NODE = ',
     +      TARGET_NODE
         call interface_MPI_ABORT(global_communicator,1)
      ENDIF
      IF (SOURCE_NODE .LT. 1 .OR. SOURCE_NODE .GT. NUMNOD) THEN
         PRINT *, 'ERROR: GLOBAL_FILE_RECORD: SOURCE_NODE = ',
     +      SOURCE_NODE
         call interface_MPI_ABORT(global_communicator,1)
      ENDIF
      IF (LOCAL_RECORD .LT. 1 .OR. LOCAL_RECORD .GT. GF_MAXREC) THEN
C increase GF_MAXREC in bucket_sort.h
         PRINT *, 'ERROR: GLOBAL_FILE_RECORD: LOCAL_RECORD = ',
     +      LOCAL_RECORD
         call interface_MPI_ABORT(global_communicator,1)
      ENDIF

      GLOBAL_FILE_RECORD =
     +   (
     +      (
     +         (
     +            TARGET_NODE * NUMNOD
     +         ) + SOURCE_NODE - 1
     +      ) * GF_MAXREC
     +   ) + LOCAL_RECORD

      RETURN
      END

      FUNCTION GLOBAL_FILE_OFFSET(
     +   TARGET_NODE, SOURCE_NODE, LOCAL_RECORD)
      use interface_to_mpi
      IMPLICIT NONE
      INTEGER(KIND=df_MPI_OFFSET_KIND) GLOBAL_FILE_OFFSET
      INTEGER TARGET_NODE, SOURCE_NODE, LOCAL_RECORD

#include "bucket_sort.h"

      INTEGER GLOBAL_FILE_RECORD
      EXTERNAL GLOBAL_FILE_RECORD

C     Prevent INTEGER overflow by splitting the following computation
C     in two parts (the first part performs an implicit cast).
      GLOBAL_FILE_OFFSET =
     +   GLOBAL_FILE_RECORD(TARGET_NODE,SOURCE_NODE,LOCAL_RECORD) - 1
      GLOBAL_FILE_OFFSET = GLOBAL_FILE_OFFSET * GF_RECLEN

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_SEEK(TARGET_NODE,SOURCE_NODE,LOCAL_RECORD)
      use interface_to_mpi
      IMPLICIT NONE
      INTEGER TARGET_NODE, SOURCE_NODE, LOCAL_RECORD

#include "bucket_sort.h"
      INTEGER(KIND=df_MPI_OFFSET_KIND) RECORD_OFFSET
      INTEGER IERR

      INTEGER(KIND=df_MPI_OFFSET_KIND) GLOBAL_FILE_OFFSET
      EXTERNAL GLOBAL_FILE_OFFSET

      RECORD_OFFSET =
     +   GLOBAL_FILE_OFFSET(TARGET_NODE,SOURCE_NODE,LOCAL_RECORD)

      CALL interface_mpi_file_seek(GF_HANDLE,RECORD_OFFSET,
     +                             df_MPI_SEEK_SET)

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_WRITE(
     +   RGBUF,IGBUF,N,SOURCE_NODE,TARGET_NODE,IBUF)
      use interface_to_mpi
      IMPLICIT NONE
      INTEGER N, SOURCE_NODE, TARGET_NODE, IBUF
      INTEGER istat(df_MPI_STATUS_SIZE)

#include "bucket_sort.h"

      REAL*8  RGBUF(NGBFSZ)
      INTEGER IGBUF(2,NGBFSZ)
      INTEGER PREV_REC, THIS_REC, INCOUNT, POSITION, IERR
      INTEGER(KIND=df_MPI_OFFSET_KIND) RECORD_OFFSET

      INTEGER(KIND=df_MPI_OFFSET_KIND) GLOBAL_FILE_OFFSET
      EXTERNAL GLOBAL_FILE_OFFSET

      CALL interface_MPI_WAIT(WRITE_REQUEST,istat)

      PREV_REC = GF_BREC(IBUF,TARGET_NODE,SOURCE_NODE)
      THIS_REC = GF_NREC(TARGET_NODE,SOURCE_NODE) + 1

      POSITION = 0

      INCOUNT = 1
      call interface_mpi_pack_i0(PREV_REC,INCOUNT,
     +  GF_RECBUF,GF_RECLEN,POSITION,self_communicator)

      INCOUNT = 1
      call interface_mpi_pack_i0(N,INCOUNT,
     +  GF_RECBUF,GF_RECLEN,POSITION,self_communicator)

      INCOUNT = NGBFSZ
      call interface_mpi_pack_r(RGBUF,INCOUNT,
     +  GF_RECBUF,GF_RECLEN,POSITION,self_communicator)

      INCOUNT = 2*NGBFSZ
      call interface_mpi_pack_i2(IGBUF,INCOUNT,
     +  GF_RECBUF,GF_RECLEN,POSITION,self_communicator)

      RECORD_OFFSET =
     +  GLOBAL_FILE_OFFSET(TARGET_NODE,SOURCE_NODE,THIS_REC)

      IERR = 0

      call interface_MPI_FILE_IWRITE_AT_pck(
     +  GF_HANDLE,RECORD_OFFSET,GF_RECBUF,GF_RECLEN,WRITE_REQUEST)

!     check for successful writing. SK - Aug 2010
!     if it fails you need to change the general record offset 
!     for a given node; for more information consult the 
!     comments in bucket_sort.h.
      IF(IERR .ne. 0)THEN
        PRINT *, 'ERROR: GLOBAL_FILE_WRITE: writing failed at
     +  RECORD_OFFSET (mytid) =',RECORD_OFFSET,SOURCE_NODE
!       errorcode 112
        call interface_MPI_ABORT(global_communicator,112)
      END IF

      GF_BREC(IBUF,TARGET_NODE,SOURCE_NODE) = THIS_REC
      GF_NREC(TARGET_NODE,SOURCE_NODE)      = THIS_REC

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_GET_RECORD(
     +   TARGET_NODE,SOURCE_NODE,LOCAL_RECORD)
      use interface_to_mpi
      IMPLICIT NONE
      INTEGER TARGET_NODE, SOURCE_NODE, LOCAL_RECORD

#include "bucket_sort.h"
      INTEGER   STATUS(df_MPI_STATUS_SIZE)
      INTEGER NUMREAD
      INTEGER IERR
      INTEGER(KIND=df_MPI_OFFSET_KIND) RECORD_OFFSET

      INTEGER(KIND=df_MPI_OFFSET_KIND) GLOBAL_FILE_OFFSET
      EXTERNAL GLOBAL_FILE_OFFSET

      RECORD_OFFSET = GLOBAL_FILE_OFFSET(
     +  TARGET_NODE,SOURCE_NODE,LOCAL_RECORD)

      CALL interface_MPI_FILE_READ_AT(
     +  GF_HANDLE,RECORD_OFFSET,GF_RECBUF,GF_RECLEN,STATUS)
CSK   ... zero NUMREAD, otherwise it may fail
      NUMREAD = 0
      call interface_mpi_get_count(STATUS,df_mpi_packed,NUMREAD)
      IF (NUMREAD .NE. GF_RECLEN) THEN
         PRINT *, 'ERROR: GLOBAL_FILE_GET_RECORD: NUMREAD = ', NUMREAD
         call interface_MPI_ABORT(global_communicator,1)
      ENDIF

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_READ(A,NROW,ICOL,IBUF)
      use interface_to_mpi
      IMPLICIT NONE
      INTEGER NROW, ICOL,IBUF
      REAL*8 A(NROW)
C
C     Read back the sorted integrals and put them at the right position
C
#include "bucket_sort.h"
#include "infpar.h"
      INTEGER SOURCE_NODE
      REAL*8 RGBUF(NGBFSZ)
      INTEGER IGBUF(2,NGBFSZ)
      INTEGER I, THIS_REC, PREV_REC, N, POSITION, OUTCOUNT, IERR

      DO I=1,NROW
         A(I) = 0D0
      ENDDO

      DO SOURCE_NODE = 1, NUMNOD

         THIS_REC = GF_BREC(IBUF,MYTID,SOURCE_NODE)

         IF (THIS_REC .EQ. 0) GOTO 20

   10    CALL LOCAL_FILE_GET_RECORD(SOURCE_NODE,THIS_REC)

            POSITION = 0

            OUTCOUNT = 1
            call interface_mpi_unpack_i0(GF_RECBUF,GF_RECLEN,POSITION,
     +         PREV_REC,OUTCOUNT,global_communicator)

            OUTCOUNT = 1
            call interface_mpi_unpack_i0(GF_RECBUF,GF_RECLEN,POSITION,
     +         N,OUTCOUNT,global_communicator)

            OUTCOUNT = NGBFSZ
            call interface_mpi_unpack_r(GF_RECBUF,GF_RECLEN,POSITION,
     +         RGBUF,OUTCOUNT,global_communicator)

            OUTCOUNT = 2*NGBFSZ
            call interface_mpi_unpack_i2(GF_RECBUF,GF_RECLEN,POSITION,
     +         IGBUF,OUTCOUNT,global_communicator)

C           Pick out the integrals that we need.

            IF (N.GT.NGBFSZ) THEN
               PRINT *, 'ERROR: GLOBAL_FILE_READ: N, NGBFSZ =',
     +                  N, NGBFSZ
               call interface_MPI_ABORT(global_communicator,1)
            ENDIF
            DO I = 1, N
               IF (IGBUF(2,I).EQ.ICOL) THEN
                  IF (IGBUF(1,I) .LT. 1 .OR. IGBUF(1,I) .GT. NROW) THEN
                     PRINT *, 'ERROR: GLOBAL_FILE_READ: IGBUF(1,I) =',
     +                        IGBUF(1,I),'I =',I
                     call interface_MPI_ABORT(global_communicator,1)
                  ENDIF
                  A(IGBUF(1,I)) = RGBUF(I)
               ENDIF
            ENDDO

            THIS_REC = PREV_REC

         IF (THIS_REC .NE. 0) GOTO 10
   20    CONTINUE 

      ENDDO

      RETURN
      END

      SUBROUTINE LOCAL_FILE_OPEN_WRITE
      use interface_to_mpi
      IMPLICIT NONE

#include "bucket_sort.h"
#include "infpar.h"

      INTEGER IERR
      INTEGER GF_INFO

      WRITE (LF_NAME,'(a,i5.5)') 'local', MYTID

      CALL interface_MPI_INFO_CREATE(GF_INFO)
#if defined (SYS_AIX)
C     ... IBMs GPFS specific file hint - SK Aug 2008
      CALL interface_MPI_INFO_SET(GF_INFO, "IBM_largeblock_io", "true")
#else
      CALL interface_MPI_INFO_SET(GF_INFO, "direct_write", "true")
      CALL interface_MPI_INFO_SET(GF_INFO, "romio_cb_write", "disable")
      CALL interface_MPI_INFO_SET(GF_INFO, "romio_ds_write", "disable")
#endif

      CALL interface_MPI_FILE_OPEN(
     +   self_communicator, LF_NAME,
     +   df_MPI_MODE_WRONLY + df_MPI_MODE_CREATE,
     +   GF_INFO, LF_HANDLE
     +)

      CALL interface_MPI_INFO_FREE(GF_INFO)

      WRITE_REQUEST = df_MPI_REQUEST_NULL

      RETURN
      END

      SUBROUTINE LOCAL_FILE_CLOSE_WRITE
      use interface_to_mpi
      IMPLICIT NONE

#include "bucket_sort.h"

      INTEGER istat(df_MPI_STATUS_SIZE)

      CALL interface_MPI_WAIT(WRITE_REQUEST,istat)
      CALL interface_MPI_FILE_CLOSE(LF_HANDLE)

      RETURN
      END

      SUBROUTINE SELECT_LOCAL_FROM_GLOBAL
      use interface_to_mpi
      IMPLICIT NONE

#include "bucket_sort.h"
#include "infpar.h"
      CHARACTER LF_RECBUF(GF_RECLEN)
      INTEGER SOURCE_NODE,LOCAL_RECORD
      INTEGER I
      INTEGER ISTAT(df_MPI_STATUS_SIZE)
      INTEGER NUMREAD
      INTEGER IERR

      CALL GLOBAL_FILE_OPEN_READ
      CALL LOCAL_FILE_OPEN_WRITE

      DO SOURCE_NODE = 1, NUMNOD

         CALL GLOBAL_FILE_SEEK(MYTID,SOURCE_NODE,1)

         DO LOCAL_RECORD = 1, GF_NREC(MYTID,SOURCE_NODE)

            CALL interface_MPI_FILE_READ_pck(
     +        GF_HANDLE,GF_RECBUF,GF_RECLEN,ISTAT)
CSK         ... zero NUMREAD, otherwise it may fail
            NUMREAD = 0
            call interface_mpi_get_count(ISTAT,df_mpi_packed,NUMREAD)
            IF (NUMREAD .NE. GF_RECLEN) THEN
               PRINT *, 'ERROR: SELECT_LOCAL_FROM_GLOBAL: NUMREAD = ',
     +         NUMREAD
               call interface_MPI_ABORT(global_communicator,1)
            ENDIF

            CALL interface_MPI_WAIT(WRITE_REQUEST,istat)

            DO I = 1, GF_RECLEN
              LF_RECBUF(I) = GF_RECBUF(I)
            ENDDO

            call interface_mpi_file_iwrite(
     +           lf_handle,lf_recbuf,gf_reclen,write_request)

         ENDDO

      ENDDO

      CALL LOCAL_FILE_CLOSE_WRITE
      CALL GLOBAL_FILE_CLOSE_READ

      RETURN
      END

      SUBROUTINE LOCAL_FILE_OPEN_READ
      use interface_to_mpi
      IMPLICIT NONE

#include "bucket_sort.h"
#include "infpar.h"

      INTEGER SOURCE_NODE

      OPEN(UNIT=LF_UNIT,FILE=LF_NAME,
     +     STATUS='OLD',ACTION='READ',
     +     ACCESS='DIRECT',RECL=GF_RECLEN)

      LF_OFFSET(0) = 0
      LF_OFFSET(1) = 0
      DO SOURCE_NODE = 2, NUMNOD
         LF_OFFSET(SOURCE_NODE) = LF_OFFSET(SOURCE_NODE-1) +
     +      GF_NREC(MYTID,SOURCE_NODE-1)
      ENDDO

      RETURN
      END

      INTEGER FUNCTION LOCAL_FILE_RECORD(SOURCE_NODE, LOCAL_RECORD)
      use interface_to_mpi
      IMPLICIT NONE
      INTEGER SOURCE_NODE, LOCAL_RECORD

#include "bucket_sort.h"
#include "infpar.h"

      IF (SOURCE_NODE .LT. 1 .OR. SOURCE_NODE .GT. NUMNOD) THEN
         PRINT *, 'ERROR: LOCAL_FILE_RECORD: SOURCE_NODE = ',
     +      SOURCE_NODE
         CALL interface_MPI_ABORT(global_communicator,1)
      ENDIF
      IF (LOCAL_RECORD .LT. 1 .OR. LOCAL_RECORD .GT. GF_MAXREC) THEN
C increase GF_MAXREC in bucket_sort.h
         PRINT *, 'ERROR: LOCAL_FILE_RECORD: LOCAL_RECORD = ',
     +      LOCAL_RECORD
         CALL interface_MPI_ABORT(global_communicator,1)
      ENDIF

      LOCAL_FILE_RECORD = LF_OFFSET(SOURCE_NODE) + LOCAL_RECORD

      RETURN
      END

      SUBROUTINE LOCAL_FILE_GET_RECORD(SOURCE_NODE,LOCAL_RECORD)
      IMPLICIT NONE
      INTEGER SOURCE_NODE, LOCAL_RECORD

#include "bucket_sort.h"

      INTEGER IREC
      INTEGER LOCAL_FILE_RECORD
      EXTERNAL LOCAL_FILE_RECORD

      IREC = LOCAL_FILE_RECORD(SOURCE_NODE,LOCAL_RECORD)
      READ(UNIT=LF_UNIT,REC=IREC) GF_RECBUF

      RETURN
      END

      SUBROUTINE LOCAL_FILE_CLOSE_READ
      IMPLICIT NONE

#include "bucket_sort.h"

      CLOSE(UNIT=LF_UNIT,STATUS='DELETE')

      RETURN
      END
C
#else /* VAR_MPI2 */ 
C aspg: dummy subroutine to avoid compilation errors in case mpi is 
C       not enabled and this file somehow ends up in a makefile
      SUBROUTINE GLOBAL_FILE_STUB
      END
#endif /* VAR_MPI2 */
