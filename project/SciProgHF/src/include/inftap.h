! include/inftap.h (from Dalton, needed for ABACUS routines)
      INTEGER LUINP, LUINTA, LUINTM, LUINTR, LUSUPM, LUONEL, LUIT1,     &
     &        LUTEMP, LUDA, LUINDX, LUPROP, LUINFO, LUSOL, LUSIFC,      &
     &        LBINTM
      COMMON /INFTAP/ LUINP, LUINTA, LUINTM, LUINTR, LUSUPM, LUONEL,    &
     &                LUIT1, LUTEMP, LUDA, LUINDX, LUPROP, LUINFO,      &
     &                LUSOL, LUSIFC, LBINTM

      INTEGER LUHALF
      EQUIVALENCE (LUTEMP, LUHALF)

      CHARACTER*8     FNONEL, FNSOL, FNSUPM, FNSIFC, LBSIFC
      COMMON /CHRTAP/ FNONEL, FNSOL, FNSUPM, FNSIFC, LBSIFC
! end of include/inftap.h
