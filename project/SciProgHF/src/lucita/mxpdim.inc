!comdeck mxpdim.inc $Revision$ 
!
      INTEGER IMAXLBLKSZ,MXPWRD
      COMMON /CMXPWRD/ MXPWRD,IMAXLBLKSZ
!
! contains all PARAMETERS defining LUCIA
      INTEGER :: MXPIRR
      PARAMETER (MXPIRR = 8)
      INTEGER :: MXPOBS
      PARAMETER (MXPOBS = 35)
!TF  Globally included here: 
      INTEGER :: MXTSOB
      parameter (MXTSOB = 35)
!SK   parameter (MXTSOB = 50)
!
      INTEGER :: MXPR4T
      PARAMETER (MXPR4T = 10)
      INTEGER :: MXPR0T
      PARAMETER (MXPR0T = 1)
      INTEGER :: MXPORB
      PARAMETER (MXPORB = 510)
      INTEGER :: MXPXOT
      PARAMETER (MXPXOT = 9)
      INTEGER :: MXPXST
      PARAMETER (MXPXST = 100)
      INTEGER :: MXPSHL
      PARAMETER (MXPSHL = 100)
      INTEGER :: MXPL
      PARAMETER (MXPL = 20)
      INTEGER :: MXPXT
      PARAMETER (MXPXT = 25)
      INTEGER :: MXPICI
      PARAMETER (MXPICI = 10)
      INTEGER :: MXPSTT
      PARAMETER (MXPSTT = 1000)
      INTEGER :: MXPCSM
      PARAMETER (MXPCSM = 20)
      INTEGER :: MXPCTP
      PARAMETER (MXPCTP = 30)
      INTEGER :: MXCNSM
      PARAMETER (MXCNSM = 8)
      INTEGER :: MTYP
      PARAMETER (MTYP = 30)
!. Note : MXPNGAS = MXPR4T+6 !!
!. Required in order to handle GAS and RAS within /LUCINP/   
      INTEGER :: MXPNGAS
      PARAMETER (MXPNGAS = 16)
      INTEGER :: MXPNSMST
      PARAMETER(MXPNSMST = 8)
!. Largest allowed division of space for perturbation operator
      INTEGER :: MXPPTSPC
      PARAMETER (MXPPTSPC=20)
!. Largest cycle of spaces in multispace method
      INTEGER :: MXPCYC
      PARAMETER(MXPCYC = 100)
!
!  For MOLUC
      INTEGER :: MXTIT
      parameter (MXTIT   =          10)
      INTEGER :: MXSEQ
      parameter (MXSEQ   =           5)
      INTEGER :: MXITEM
      parameter (MXITEM  =          32)
!     parameter (MXROOT  =          25) ! defined in rasdim.inc
