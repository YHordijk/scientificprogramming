!*
      INTEGER IMAXLBLKSZ
#ifdef LUCIAREL_STANDALONE
      PARAMETER (MXPWRD = 8 000 000)
#else
      INTEGER MXPWRD
      COMMON /CMXPWRD/ MXPWRD, IMAXLBLKSZ
#endif
!
!     this file contains all required PARAMETERS for 
!     static dimensioning in LUCIAREL
!
!     max. number of roots per irrep
      INTEGER :: MXROOT
      parameter (MXROOT = 64  )

!     MXPSTT == max. number of symmetry type-types
!     decrease/increase when running smaller/larger systems
      INTEGER :: MXPSTT
!TF   PARAMETER (MXPSTT = 12500 )
      PARAMETER (MXPSTT = 50000 )

!     MXINT_TP == max. number of SO-integral types 
      INTEGER MXINT_TP
      PARAMETER (MXINT_TP = 10000)
!     MXINT_DIA == max. number of diagonal integrals  
!     decrease/increase when running smaller/larger systems
      INTEGER :: MXINT_DIA
      PARAMETER (MXINT_DIA = 300000)
      INTEGER :: MXPICI
      PARAMETER (MXPICI = 10  )
!
!     max. number of GASpaces
      INTEGER :: MXPNGAS
      PARAMETER (MXPNGAS = 8)
!     max. no. of double group irreps, limited to D32H: 128
      INTEGER :: MXNDGIRR
      parameter (MXNDGIRR = 128 )
!     max. no. of point group irreps, limited to D32H: MXNDGIRR/2
      INTEGER :: MXPIRR
      PARAMETER (MXPIRR = MXNDGIRR/2)
!     max. no. of orbital symmetries - to be consistent with the code (as it is): MXNDGIRR
      INTEGER :: MXPOBS
      PARAMETER (MXPOBS = MXNDGIRR  )
!     max. no. of spinors
      INTEGER :: MXPORB
      PARAMETER (MXPORB = 800)
!     whatever the original meaning of MXPCSM was, it is now used as equivalent to MXNDGIRR, 
!     therefore:  (discovered while implementing linear symmetry)
      INTEGER :: MXPCSM
      PARAMETER (MXPCSM = MXNDGIRR)
!
!     ALLOWD NUMBER OF MK SPACES
      INTEGER :: MXPNMK
      PARAMETER(MXPNMK=100)
!     preliminary for SIGDEN route (this seems like a random guess to me - FIXME?) 
!                S. Knecht March 2011
      INTEGER :: MXNOP
      INTEGER :: MXPLCCOP
      parameter (MXNOP    =  32)    ! Max. no. of exc. operators
      parameter (MXPLCCOP =  MXNOP)
