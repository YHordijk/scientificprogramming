!  FILE : krmcluci_inf.h
!  Information used for KRMC and LUCIAREL interface
!  ------------------------------------------------
!
!  1) CI energy
!  2) IRIQ for real/imaginary mismatch when running real groups
!     in LUCIAREL (preliminary implementation)
!
      REAL*8 ECI
      INTEGER IRIQ
      COMMON /CIVCINF/ ECI, IRIQ
