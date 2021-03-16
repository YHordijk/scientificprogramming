      LOGICAL SOLVNT
      REAL*8 ORICAV(3), RCAV(3), EPDIEL
      INTEGER LCAVMX, LMTOT, LMNTOT, NCNTCV
      COMMON /CBISOL/ ORICAV, RCAV, EPDIEL, LCAVMX, LMTOT, LMNTOT,      &
     &                NCNTCV, SOLVNT

!     ... LCAVMX = max l in moment expansion for solvent cavity
!         LMTOT  = number of spherical components for R(l,m)
!         LMNTOT = number of cartesian components for RC(L,M,N)
