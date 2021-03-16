      SUBROUTINE PRQMAT(QMAT,NROW,NCOL,LRQ,LCQ,NZ,IQP,IUNIT)
#include "implicit.h"
      DIMENSION QMAT(LRQ,LCQ,NZ),IQP(*)
      CHARACTER M(4)*6
      M(1) = 'A real'
      M(2) = 'A imag'
      M(3) = 'B real'
      M(4) = 'B imag'
      DO 10 IZ = 1,NZ
        IQ = IQP(IZ)
        WRITE(IUNIT,'(/A,A6,A,4I4/)')
     &  '*** ',M(IQ),' part ***; NROW,NCOL,LRQ,LCQ:',NROW,NCOL,LRQ,LCQ
        CALL OUTPUT(QMAT(1,1,IZ),1,NROW,1,NCOL,LRQ,LCQ,-1,IUNIT)
   10 CONTINUE
      END
!C  /* Deck readi */
      SUBROUTINE READI (IT,N,INTX)
C
C (30-Jan-1984) hjaaj
C
      DIMENSION INTX(N)
      READ (IT) INTX
      RETURN
      END
C  /* Deck readi2 */
      SUBROUTINE READI2 (IT,N,INTX)
C
C (30-Jan-1984) hjaaj
C
#if defined (VAR_STAR2)
      INTEGER*2 INTX(N)
#else
      INTEGER   INTX(N)
#endif
      READ (IT) INTX
      RETURN
      END
C  /* Deck readdi */
      SUBROUTINE READDI(IT,IU,N,IX)
      DIMENSION IX(N)
#if defined (VAR_FASTIO)
      ERROR - READDI is not programmed for VAR_FASTIO
#else
      READ(IT,REC=IU) IX
#endif
      RETURN
      END
C  /* Deck readt */
      SUBROUTINE READT (IT,N,X)
#include "implicit.h"
      CHARACTER*30 STRING
      DIMENSION X(N)
      READ (IT,END=10,ERR=20) X
      RETURN
 10   CONTINUE
      WRITE (STRING,'(A25,I5)') 'READT: END reading unit  ',IT
      CALL QUIT(STRING)
 20   CONTINUE
      WRITE (STRING,'(A25,I5)') 'READT: Error reading unit',IT
      CALL QUIT(STRING)
      END
!C  /* Deck writi */
      SUBROUTINE WRITI (IT,N,INTX)
!C
!C (30-Jan-1984) hjaaj
!C
      DIMENSION INTX(N)
      WRITE (IT) INTX
      RETURN
      END
C  /* Deck writi2 */
      SUBROUTINE WRITI2 (IT,N,INTX)
C
C (30-Jan-1984) hjaaj
C
#if defined (VAR_STAR2)
      INTEGER*2 INTX(N)
#else
      INTEGER   INTX(N)
#endif
      WRITE (IT) INTX
      RETURN
      END
C  /* Deck writdi */
      SUBROUTINE WRITDI(IT,IU,N,IX)
      DIMENSION IX(N)
#if !defined (VAR_FASTIO)
      WRITE(IT,REC=IU) IX
#endif
      RETURN
      END
C  /* Deck writt */
      SUBROUTINE WRITT(IT,N,X)
#include "implicit.h"
      DIMENSION X(N)
      WRITE (IT) X
      RETURN
C
      END
C  /* Deck open_dx */
      SUBROUTINE OPEN_DX (LU_DX,NAME,NELEM,STATUS,OLD_DX)
C
C Original: 15-Jun-1985 hjaaj
C
C Purpose:
C   Open files for direct access through WRIT_DX and READ_DX routines.
C   The ...._DX routines enables direct access, even when the number
C   of elements per record (the logical record length) is greater
C   than the maximum physical record length.
C   >>> THIS IS MACHINE DEPENDENT <<<
C
C Input:
C  LU_DX     file unit number
C  NAME      name of file to be opened
C  NELEM     number of integer words per logical record
C  STATUS    file status: 'OLD', 'NEW', or 'UNKNOWN'
C
C Output on common block:
C  LRECL_DX  physical record length (in integers)
C  NBUF_DX   number of physical records per logical record
C Output:
C  OLD_DX    logical, true if old LU_DX file was opened
C
C  Why these routines ???
C     We (Stefan Knecht and I) ran into problems with KRMCSCF
C     when record length in bytes became greater then 2Gb,
C     even though there is no 2Gb limit on the file size.
C     Apparently compiler programmers have used
C     32bit integers internally for the record length !
C     Therefore I dug up this old package for direct
C     access with logical record length possibly bigger
C     than the physical recored length.
C     It might also not be optimal to write 1Gb in one shot,
C     so we chose for now to define the maximum physical
C     record length to 8*10**7 bytes (for INT_STAR8) instead of
C     2*10**9 bytes (=2Gb).
C     -- Odense, May 18, 2011, Hans Joergen Aa. Jensen
C     
C
      CHARACTER*(*) NAME, STATUS
      LOGICAL OLD_DX
      PARAMETER (MAX_LRECL_DX = 10000000) ! 10 000 000
#include "priunit.h"
      integer NBUF_DX, LRECL_DX
      COMMON /CB_DX/ NBUF_DX(99), LRECL_DX(99)
C
      IF (LU_DX .GT. 99 .OR. LU_DX .LE. 0) THEN
         CALL QENTER('OPEN_DX')
         CALL QUIT('File unit number not in range 1:99')
      END IF
C
      NBUF_DX(LU_DX)   = (NELEM-1)/MAX_LRECL_DX + 1
      LRECL_DX(LU_DX)  = (NELEM-1)/NBUF_DX(LU_DX) + 1
      LRECL            = 8*LRECL_DX(LU_DX)
C
      IF (STATUS .EQ. 'NEW') GO TO 300
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 9000
C
C     OPEN OLD FILE
C
         OPEN(LU_DX,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',ERR=300,
     *        ACCESS='DIRECT',RECL=LRECL)
         OLD_DX = .TRUE.
      GO TO 600
C
  300 CONTINUE
      IF (STATUS .EQ. 'OLD') GO TO 9100
C
C     OPEN NEW FILE
C
         OPEN(LU_DX,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',RECL=LRECL)
         OLD_DX = .FALSE.
  600 CONTINUE

      RETURN
C
C error branches
C
 9000 CONTINUE
      WRITE (LUPRI,'(//A,A/A,I5/,A,A)')
     *   ' *** ERROR (OPEN_DX) INVALID STATUS KEYWORD: ',STATUS,
     *   '                    FILE NUMBER =',LU_DX,
     &   '                    FILE NAME = ',NAME
      CALL QTRACE(LUPRI)
      CALL QUIT('*** ERROR (OPEN_DX) INVALID STATUS KEYWORD')
C
 9100 CONTINUE
      WRITE (LUPRI,'(//A/A,I5/A)')
     *   ' *** ERROR (OPEN_DX) OLD FILE NOT FOUND',
     *   '                    FILE NUMBER =',LU_DX,
     *   ' --- or wrong record length on old file.'
      CALL QTRACE(LUPRI)
      CALL QUIT('*** ERROR (OPEN_DX) FILE NOT FOUND')
C
C end of OPEN_DX
C
      END
C  /* Deck FIND_DX */
      LOGICAL FUNCTION FIND_DX(LU_DX,record_number,LEN_VEC,VEC)
C
C 27-Jun-1985 Hans Jorgen Aa. Jensen
C
C For direct access find record,
C when len_vec may be greater than maximum record length.
C
C See comments in OPEN_DX for more explanations.
C
      integer NBUF_DX, LRECL_DX
      COMMON /CB_DX/ NBUF_DX(99), LRECL_DX(99)
      integer, intent(in)    :: LU_DX
      integer, intent(in)    :: record_number, len_vec
      real(8), intent(out)   :: VEC(LEN_VEC)

      IF (len_vec .LE. LRECL_DX(LU_DX)) THEN
         READ (LU_DX, REC=record_number, IOSTAT=IOS) VEC
         IF (IOS .NE. 0) GO TO 900
      ELSE
         NBUF_DX_TEST = (len_vec-1)/LRECL_DX(LU_DX) + 1
         IF (NBUF_DX_TEST .gt. NBUF_DX(LU_DX)) THEN
            CALL QUIT(' recalculated # of buffers > # of buffers'//
     &                ' calculated when opening the file.')
         END IF
         IREC = 1 + NBUF_DX(LU_DX)*(record_number-1)
         JADD = 0
         DO 100 IBUF = 1,NBUF_DX(LU_DX)-1
            READ (LU_DX, REC=IREC, IOSTAT=IOS) ! replace with READ_DI call
     &         (VEC(JADD+J), J = 1,LRECL_DX(LU_DX))
            IF (IOS .NE. 0) GO TO 900
            IREC = IREC + 1
            JADD = JADD + LRECL_DX(LU_DX)
  100    CONTINUE
         READ (LU_DX, REC=IREC, IOSTAT=IOS) (VEC(J), J = JADD+1,len_vec)
         IF (IOS .NE. 0) GO TO 900
      END IF
      FIND_DX = .TRUE.
      RETURN
C
  900 CONTINUE
      FIND_DX = .FALSE.
      RETURN
      END
C  /* Deck READ_DX */
      SUBROUTINE READ_DX(LU_DX,record_number,len_vec,VEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
C For direct access read
C when len_vec may be greater than maximum record length.
C
C See comments in OPEN_DX for more explanations.
C
      integer NBUF_DX, LRECL_DX
      COMMON /CB_DX/ NBUF_DX(99), LRECL_DX(99)
      integer, intent(in)    :: LU_DX
      integer, intent(in)    :: record_number, len_vec
      real(8), intent(out)   :: VEC(LEN_VEC)
      IF (len_vec .LE. LRECL_DX(LU_DX)) THEN
         READ (LU_DX, REC = record_number) VEC
      ELSE
         NBUF_DX_TEST = (len_vec-1)/LRECL_DX(LU_DX) + 1
         IF (NBUF_DX_TEST .gt. NBUF_DX(LU_DX)) THEN
            CALL QUIT(' recalculated # of buffers > # of buffers'//
     &                ' calculated when opening the file.')
         END IF
         IREC = 1 + NBUF_DX(LU_DX)*(record_number-1)
         JADD = 0
         DO 100 IBUF = 1,NBUF_DX(LU_DX)-1
            CALL READAC(LU_DX,LRECL_DX(LU_DX),VEC(JADD+1),IREC)
!           READ (LU_DX, REC = IREC) 
!    &           (VEC(JADD+J), J = 1,LRECL_DX(LU_DX))
            IREC = IREC + 1
            JADD = JADD + LRECL_DX(LU_DX)
  100    CONTINUE
         CALL READAC(LU_DX,len_vec-JADD,VEC(JADD+1),IREC)
!        READ (LU_DX, REC = IREC) (VEC(J), J = JADD+1,len_vec)
      END IF
      RETURN
      END
C  /* Deck WRIT_DX */
      SUBROUTINE WRIT_DX(LU_DX,record_number,len_vec,VEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
C For direct access write
C when len_vec may be greater than maximum record length.
C
C See comments in OPEN_DX for more explanations.
C
      integer NBUF_DX, LRECL_DX
      COMMON /CB_DX/ NBUF_DX(99), LRECL_DX(99)
      integer, intent(in)    :: LU_DX
      integer, intent(in)    :: record_number, len_vec
      real(8), intent(in)    :: VEC(LEN_VEC)

      IF (len_vec .LE. LRECL_DX(LU_DX)) THEN
         WRITE (LU_DX, REC = record_number) VEC
      ELSE
         NBUF_DX_TEST = (len_vec-1)/LRECL_DX(LU_DX) + 1
         IF (NBUF_DX_TEST .gt. NBUF_DX(LU_DX)) THEN
            CALL QUIT(' recalculated # of buffers > # of buffers'//
     &                ' calculated when opening the file.')
         END IF
         IREC = 1 + NBUF_DX(LU_DX)*(record_number-1)
         JADD = 0
         DO 100 IBUF = 1,NBUF_DX(LU_DX)-1
            CALL WRTDAC(LU_DX,LRECL_DX(LU_DX),VEC(JADD+1),IREC)
            IREC = IREC + 1
            JADD = JADD + LRECL_DX(LU_DX)
  100    CONTINUE
         CALL WRTDAC(LU_DX,len_vec-JADD,VEC(JADD+1),IREC)
!        WRITE (LU_DX, REC = IREC) (VEC(J), J = JADD+1,len_vec)
      END IF
      RETURN
      END
