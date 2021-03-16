!*** Information for generating two-electron distributions
!    with two indices fixes
!     Length of common block : keep up to date !!
      INTEGER NCOMDIS
      PARAMETER (NCOMDIS = 224)

      INTEGER NSPCK(0:7, 0:3), ISPCK(0:7, 0:7, 3)
      COMMON /COMDIS/ NSPCK, ISPCK

      INTEGER NPASS, NGBFSZ, LGFIL, IRECG, KLGREC, KLGBUF, KRGBUF,      &
     &        KIGBUF
      COMMON /DISBUF/ NPASS, NGBFSZ, LGFIL, IRECG, KLGREC, KLGBUF,      &
     &                KRGBUF, KIGBUF

      INTEGER ISHLA, ISHLB, ISETA, ISETB, NINSHA, NINSHB
      COMMON /DIS2AB/ ISHLA, ISHLB, ISETA, ISETB, NINSHA, NINSHB

