! FILE: relccsd/freeze.inc
!
      INTEGER   MAXFROZ, MAXIREP, NFROZEN,IFROZ
      PARAMETER (MAXFROZ=100, MAXIREP=32)
      LOGICAL   IFROZ_set_in_input
      COMMON /RELCCSD_FREEZE/ NFROZEN,IFROZ(MAXFROZ,MAXIREP),           &
     &          IFROZ_set_in_input
! -- end of relccsd/freeze.inc
