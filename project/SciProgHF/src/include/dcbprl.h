!
!     PROPERTY LABELS:
!
!       PRPLBL  - property label
!       IPRLREP  - boson irrep of scalar operator (0 - 7)
!       IPRLTYP - symmetry of scalar operator under matrix transpose (-1,0,1)
!
      INTEGER MAXPRPLBL
      PARAMETER (MAXPRPLBL = 3200)

      CHARACTER*8 PRPLBL(MAXPRPLBL)
      CHARACTER*4 PDOINT(MAXPRPLBL)
      INTEGER NPRPLBL, IPRLREP(MAXPRPLBL), IPRLTYP(MAXPRPLBL)
      COMMON /XCBPRL/ NPRPLBL, IPRLREP, IPRLTYP

      COMMON /XCBLBL/ PRPLBL, PDOINT
