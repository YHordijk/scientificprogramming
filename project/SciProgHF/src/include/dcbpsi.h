!
! file: dcbpsi.h -- Dirac Common Block for PSI (selection of wave functions)
!
! *** Wave function module ***
!
!  DOMP2    - perform direct second-order Moeller-Plesset calculation
!  L1ORBM   - orbital pre-modifications (reordering, rotations)
!  L2ORBM   - orbital post-modifications (reordering, phase adjustment)
!  DORES    - resolution of open-shell manifold with GOSCIP
!  DOMVO    - make modified virtual orbitals
!  DOMP2NO  - make MP2 natural orbitals
!  DOEXACC  - perform coupled-cluster calculation with EXA_CC
!  DOCCM    - perform coupled-cluster calculation with RELCCSD
!  DOHSCC   - perform HSMRCC calculation with external code
!  DOHSFS   - perform higher sector FSCC calculations with external code
!  DOCOSCI  - perform CI calculation with GOSCIP
!  DOCIM    - perform CI calculation with DIRRCI/GOSCIP
!  DOLUCT   - perform CI calculation with LUCITA
!  DOARDUC  - perform CC calculation with ARDUCCA
!  DOADC    - compute ionization spectrum via Green's functions
!  DOKRCI   - perform CI calculation with LUCIAREL/...
!  DOKRCC   - perform coupled-cluster calculation with LUCIAREL
!  DOCOSCI  - do COSCI calculation
!  DOGASCIP - do GASCIP calculation
!  DOPOLPRP - do relativistic polarization propagator for excitations
!  DOLAPLCE - do Laplace transformation of orbital energy denominator

      LOGICAL                      DOMP2, DOMVO, DORES, DOCCM, DOCIM,   &
     &        DOHSCC, DOHSFS, DOLUCT, DOKRMC, DOLUCIAR, DOMP2NO,        &
     &        L1ORBM, L2ORBM, PHCOEF, DOCOSCI, DOADC, DOKRCI, DOKRCC,   &
     &        DOGASCIP, DOARDUC, DOPOLPRP, DOLAPLCE, DOEXACC

      COMMON /DCLPSI/                                                   &
     &                             DOMP2, DOMVO, DORES, DOCCM, DOCIM,   &
     &        DOHSCC, DOHSFS, DOLUCT, DOKRMC, DOLUCIAR, DOMP2NO,        &
     &        L1ORBM, L2ORBM, PHCOEF, DOCOSCI, DOADC, DOKRCI, DOKRCC,   &
     &        DOGASCIP, DOARDUC, DOPOLPRP, DOLAPLCE, DOEXACC

! --- end of dcbpsi.h ---
