C
C     Define common blocks and other data for bucket sort
C
C     Offset for the file unit
      INTEGER I_BF_UNIT
      PARAMETER (I_BF_UNIT=70)
C
C     The maximum number of files that can be opened.
      INTEGER MAX_BUFF_FILES
      PARAMETER (MAX_BUFF_FILES=5)
C
C     The length of records, measured in integer words.
C     Performance tip : Make this an multiple of 1024 MINUS 1 (to leave space for the 2 additional administrative ints)
C
      INTEGER NADMIN,NKWORDS,NKBUFFS,NGBFSZ
      PARAMETER (NADMIN=1)
!
!     !!! CAUTION!!! the current value for this parameter leads to 
!     crashes on various machines BRUTUS/ETH Zuerich, Horseshoe/SDU Odense, my laptop, ... 
!     (I could not test IBM machines since I do not have access anymore).
!     SK - Aug 2010.
!
!     example: node 1 will always try to read/write the first record at:
!     (( NKWORDS * 1024 *4)**2) * 4 * 64 > 4.38 * 10^12.
!
      PARAMETER (NKWORDS=32)
      PARAMETER (NKBUFFS=1)
      PARAMETER (NGBFSZ=NKWORDS*1024-NADMIN)
C
C     The maximum number of buckets per file.
      INTEGER N_BUFFER
      PARAMETER (N_BUFFER=NKBUFFS*1024)
C
C     The memory needed to store all buckets for one sort (need MAX_MEM
C     times one real and two integer numbers)
C
      INTEGER MAX_MEM
      PARAMETER (MAX_MEM=NGBFSZ*N_BUFFER)
C
C     The file name for each file. This name is made unique
C     by appending the ID of the node in parallel runs, in
C     in serial runs we append 0.
C
      COMMON/ BF_FLNAME/ BF_NAME(MAX_BUFF_FILES)
      CHARACTER*10 BF_NAME
C
C     Information about the sort and the files :
C
C     - The highest record written
C     - For each buffer
C        - the last record written
C        - the amount of unflushed data in the buffer 
C
      INTEGER I_BF_REC, LGREC, LGBUF
      COMMON /BF_RECS/ I_BF_REC(MAX_BUFF_FILES),
     &                 LGREC(N_BUFFER,MAX_BUFF_FILES),
     &                 LGBUF(N_BUFFER,MAX_BUFF_FILES)
C
C     Information necessary for distributed write and associated
C     communication
C
C     The maximum number of send buffers that is used.
      INTEGER MAX_SEND_BUFFERS
      PARAMETER (MAX_SEND_BUFFERS=5)
C     Buffers for sending integrals
      INTEGER ISGBUF, IREQ
      COMMON /BF_SENDI/ ISGBUF(4+4*NGBFSZ,MAX_SEND_BUFFERS),
     &                  IREQ(MAX_SEND_BUFFERS)
C
C     Statistics
      INTEGER*8 N_WRITTEN, N_DELETED
      COMMON /BF_SENDI8/ N_WRITTEN, N_DELETED

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c aspg, 2009-04-24
c
c variables used in the fine-sorting of the 1HT of scheme6
c
      INTEGER MAX_BUF_files, I_SORTED_UNIT
      PARAMETER (MAX_BUF_files=1)
      PARAMETER (I_SORTED_UNIT=400)

      integer KLpairsInBuffer, max_records_per_buffer

      integer node_for_buffer(N_BUFFER)
      integer records_in_buffer(N_BUFFER)
      integer irec_md_lastwritten(N_BUFFER)
      integer irec_md_global_lastwritten

      INTEGER i_md_unit
      parameter(i_md_unit=931)

      common/kl_buf_info/KLpairsInBuffer,max_records_per_buffer,
     &       records_in_buffer, node_for_buffer,
     &       irec_md_lastwritten, irec_md_global_lastwritten
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER GF_UNIT
      PARAMETER (GF_UNIT = 91)
      CHARACTER*(*) GF_NAME
      PARAMETER (GF_NAME = 'gf-matrix')

      INTEGER LF_UNIT
      PARAMETER (LF_UNIT = 92)
      CHARACTER*(10) LF_NAME

      INTEGER BLOCKSIZE
      PARAMETER (BLOCKSIZE = (NGBFSZ+1)*4)

      INTEGER MAX_MPI_SIZE
      PARAMETER (MAX_MPI_SIZE = 32)

      INTEGER GF_RECLEN, GF_MAXREC
#if defined (INT_STAR8)
      PARAMETER (GF_RECLEN = 8 * BLOCKSIZE)
#else
      PARAMETER (GF_RECLEN = 4 * BLOCKSIZE)
#endif
      PARAMETER (GF_MAXREC = 64 * BLOCKSIZE)

      INTEGER GF_BREC(N_BUFFER, 0:MAX_MPI_SIZE-1, 0:MAX_MPI_SIZE-1)
      INTEGER GF_NREC(0:MAX_MPI_SIZE-1, 0:MAX_MPI_SIZE-1)
      INTEGER LF_OFFSET(0:MAX_MPI_SIZE-1)
      INTEGER GF_HANDLE, LF_HANDLE
      INTEGER WRITE_REQUEST

      CHARACTER GF_RECBUF(GF_RECLEN)

      COMMON /GF_STATE/  GF_BREC, GF_NREC, LF_OFFSET,
     +   GF_HANDLE, LF_HANDLE, WRITE_REQUEST
      COMMON /GF_BUFFER/ GF_RECBUF
      COMMON /LF_NAME/ LF_NAME
