      INTEGER PRINT_LEVEL
      INTEGER RECP_OUT      ! RECP output
      INTEGER RECP_DBG      ! For debug program
      INTEGER RECP_FNUM(10) ! I/O file number
      CHARACTER*10 RECP_FNAM(10) ! I/O file number
      COMMON  /RECP_VAR_PRINT/                                          &
     &        PRINT_LEVEL, RECP_OUT, RECP_DBG, RECP_FNUM, RECP_FNAM
