! ******************************************
! ***            NMR module              ***
! ***
! ***   LUTWOL - file number for two-ele.contrib. of LAO
! ***   LUFCKL - file number  of LAO Fock matrix
! 
!  NOTWOL - do not include two-electron London contribution
!  NOOONEI - do not include reorthonormalization contributions
!  INTNMR - what kind of integrals for LAO 2-el.contribution
! 
! ******************************************
      INTEGER LUTWOL, LUFCKL
      PARAMETER (LUTWOL = 69,LUFCKL = 70)
      LOGICAL NOTWOL, NOONEI, EPREORTH
      COMMON/CBLNMR/NOTWOL, NOONEI, EPREORTH
      INTEGER INTNMR
      COMMON/CBINMR/INTNMR
