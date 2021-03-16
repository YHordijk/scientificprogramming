!*. CSTATE

      REAL*8 ETHRES,CTHRES,EROOT
      
      INTEGER IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS,                     &
     &        NROOT,IROOT,IDC,PSSIGN,PLSIGN,INTSEL,IALIGN,              &
     &        NMS2VAL,MS2VAL,                                           &
     &        NBLK_MS2,NBAT_MS2,IBLK_MS2,                               &
     &        IBAT_MS2,NBLK_MS2_C, NMAX_SYM_CI

      COMMON/CSTATER/IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS,              &
     &              NROOT,IROOT(100),IDC,PSSIGN,PLSIGN,INTSEL,IALIGN,   &
     &              ETHRES,CTHRES,EROOT(MXROOT),NMS2VAL,MS2VAL(100),    &
     &              NBLK_MS2(100),NBAT_MS2(100),IBLK_MS2(100),          &
     &              IBAT_MS2(100),NBLK_MS2_C(100), NMAX_SYM_CI
      CHARACTER*1 SYMLABEL
