!
!  DCBRES:  Parameters for resolution of open-shell manifold
!
! Updated by Miro Ilias, August 2016 for DFT-COSCI
!
     
      INTEGER IPRRES,ISTRES,INTRES
      COMMON/DCIRES/IPRRES,ISTRES,INTRES

!
! NOVXC - do not add Vxc term to transformed Fij - no DFT at COSCI level ! (rough approx.)
! COREDENS - only core-density for Vxc terms (poor approxim)
!
!
      LOGICAL NOVXC,COREDENS
      COMMON/DCLRES/NOVXC,COREDENS
    
      REAL*8 SCRRES,RESOUT
      COMMON/DCRRES/SCRRES,RESOUT


