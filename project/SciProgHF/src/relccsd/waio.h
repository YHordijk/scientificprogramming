C
C     Define common blocks and other data for word-addressable I/O
C
C     The length of records, measured in integer words.
      INTEGER*8 NWORDS_ON_RECORD8
      INTEGER   NWORDS_ON_RECORD
      PARAMETER (NWORDS_ON_RECORD8=4096)
      PARAMETER (NWORDS_ON_RECORD=4096)
C
C     The maximum number fo files that can be opened.
      INTEGER MAX_WAIO_FILES
      PARAMETER (MAX_WAIO_FILES=100)
C
C     The file name for each file. This name is made unique 
C     by appending the ID of the node in parallel runs, in
C     in serial runs we append 0.
C
      COMMON/ FLNAME/ DN(MAX_WAIO_FILES)
      CHARACTER*10 DN
C
C     The position of the file : this is the highest record written
C     by the WAIO routines.
C
      COMMON /POINTR/ IWRIT_WAIO(MAX_WAIO_FILES)
      INTEGER IWRIT_WAIO
C
C     Some statistics to compute I/O performance
C
      REAL*8 WAIO_WCALL,WAIO_RCALL,WAIO_WTIME,WAIO_RTIME,
     &       WAIO_WNLEN,WAIO_RNLEN
      COMMON /WAIO_STAT/ WAIO_WCALL,WAIO_RCALL,
     & WAIO_WTIME,WAIO_RTIME,WAIO_WNLEN,WAIO_RNLEN


