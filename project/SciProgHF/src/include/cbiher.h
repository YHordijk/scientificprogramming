! FILE: include/cbiher.h
! (Control information for HERMIT)
      REAL*8  EXPKR(3), THRESH, PRTHRS
      INTEGER IPRDEF, IORCAR, IORSPH, NPQUAD, NPATOM, IPATOM(MXCENT)
      LOGICAL SKIP,   TSTINP, HAMILT, DIPLEN, SPNORB, SOTEST, SUPMAT,   &
     &        DIPVEL, QUADRU, SECMOM, ONEPRP, CARMOM, SPHMOM, OCTUPO,   &
     &        FERMI,  PSO, SPIDIP, DSO, ALLATM, NMRISS, TRIANG, SDFC,   &
     &        PROPRI, HDO, S1MAG, S2MAG, ANGMOM, ANGLON, LONMOM, MAGMOM,&
     &        MGCOOR, S1MAGT, MGMOMT, KINENE, S2MAGT, DSUSNL, DSUSLL,   &
     &        DSUSLH, DIASUS, DSUTST, NUCSNL, NUCSLO, NUCSHI, NSNLTS,   &
     &        NSLTST, NELFLD, NSTTST, EFGCAR, EFGSPH, NUCPOT, S1MAGL,   &
     &        S1MAGR, HDOBR, S1MLT, S1MRT, HDOBRT, NPOTST, MGMO2T, HBDO,&
     &        SUSCGO, NSTCGO, EXPIKR, DIRAC, DARWIN, MASSVL, CM1, CM2,  &
     &        NOTWO, SQHDOL, SQHDOR, S1ELE, S1ELB, ONEELD, GFACDI,      &
     &        THETA, CXIKR
      COMMON /CBIHER/ EXPKR,  THRESH, PRTHRS,                           & ! real*8
     &                IPRDEF, IORCAR, IORSPH, NPQUAD, NPATOM, IPATOM,   & ! integers
     &                SKIP, TSTINP, HAMILT, DIPLEN, SPNORB,             & ! LOGICALS
     &                SOTEST, SUPMAT, DIPVEL, QUADRU, SECMOM, ONEPRP,   &
     &                CARMOM, SPHMOM, OCTUPO, FERMI, PSO, SPIDIP, DSO,  &
     &                ALLATM, NMRISS, TRIANG, SDFC,                     &
     &                PROPRI, HDO, S1MAG, S2MAG, ANGMOM, ANGLON, LONMOM,&
     &                MAGMOM, S1MAGT, MGMOMT, KINENE, S2MAGT, DSUSNL,   &
     &                DSUSLL, DSUSLH, DIASUS, DSUTST, NUCSNL, NUCSLO,   &
     &                NUCSHI, NSNLTS, NSLTST, NELFLD, NSTTST, EFGCAR,   &
     &                EFGSPH, NUCPOT, S1MAGL, S1MAGR, HDOBR, S1MLT,     &
     &                S1MRT, HDOBRT, NPOTST, MGMO2T, HBDO,              &
     &                SUSCGO, NSTCGO, EXPIKR, DIRAC, DARWIN, MASSVL,    &
     &                CM1, CM2, NOTWO, SQHDOL, SQHDOR, S1ELE, S1ELB,    &
     &                ONEELD, GFACDI, THETA, MGCOOR, CXIKR
! -- end of include/cbiher.h --
