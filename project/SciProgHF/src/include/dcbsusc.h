!-----------------------------------------------------------------------
!  Variables for calculating the magnetic susceptibility;
!
! ... when using London atomic orbitals (LAO)
!   SUSDIA  One-electron expectation values
!   SUS2EL  Two-electron expectation values
!   SUSFS   Highest-order reorthonormalization
!   SUSFSY  Lowest-order reorthonormalization
!   SUSREL  Relaxation
!
! ... without using LAO
!   SUSCOM
!
!
      REAL*8 SUSDIA(3,3), SUS2EL(3,3), SUSFS(3,3), SUSREL(3,3),
     &       SUSTOT(3,3), SUSFSY(3,3), SUSCOM(3,3)
      COMMON /SUSCPT/ SUSDIA, SUS2EL, SUSFS, SUSREL,
     &                SUSTOT, SUSFSY, SUSCOM
