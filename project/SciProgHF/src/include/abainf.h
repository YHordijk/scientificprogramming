! FILE: include/abainf.h
! (Control information for ABACUS routines)
      INTEGER IPRDEF, NWNABA, NTCSY(8), NSCSY(8)
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR, TSTINP, VIB, RESTAR,       &
     &        DOWALK, GDALL, CCSD, H2MO, DOSYM(8), DOLRES, DOEXCI,      &
     &        SHIELD, SPNSPN, MAGSUS, VCD, NACME, AAT, NOLOND, FCKDDR,  &
     &        ECD, NODIFC, DODRCT, SUPMAT, PARALL, ROTG, SPINRO,        &
     &        MASSVE, DARWIN, ABALNR, VROA, NOCMC, EXPFCK, RAMAN,       &
     &        QUADRU, NQCC, HYPER, VERDET, MCD, GSDIP, GSQUAD, GSOCT,   &
     &        GSDIDI, GSDIQU, GSQUQU,                                   &
     &        lincpl, numhes
      COMMON /ABAINF/ IPRDEF, NWNABA, NTCSY, NSCSY,                     & ! integers
     &                MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP, VIB,      & ! logicals
     &                RESTAR, DOWALK, GDALL, CCSD, H2MO,                &
     &                DOSYM, DOLRES, DOEXCI, SHIELD,                    &
     &                SPNSPN, MAGSUS, VCD, NACME, AAT, NOLOND, FCKDDR,  &
     &                ECD, DODRCT, SUPMAT, PARALL, ROTG, SPINRO,        &
     &                MASSVE, DARWIN, ABALNR, VROA, NODIFC,             &
     &                NOCMC, EXPFCK, RAMAN, QUADRU, NQCC, HYPER, VERDET,&
     &                MCD, GSDIP, GSQUAD, GSOCT,                        &
     &                GSDIDI, GSDIQU, GSQUQU,                           &
     &                lincpl, numhes
! -- end of include/abainf.h --
