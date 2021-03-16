!
!     File: cbirea.h - common block for reading of MOLECULE.INP (herrdn.F)
!     MAXPRD = default value for MAXPRI

      INTEGER MAXPRD
      PARAMETER ( MAXPRD = 35 )

      REAL*8  ZCMVAL, TOL_SYMADD
      INTEGER IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND
      LOGICAL BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  PRIBAS, ATOMBA,   &
     &        UNCONT, CNTMAT, LCNTNUUM
      COMMON /CBIREA/ ZCMVAL, TOL_SYMADD,                               &
     &        IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND,           &
     &        BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  PRIBAS, ATOMBA,   &
     &        UNCONT, CNTMAT, LCNTNUUM

!     MAXFAMEXP = maximum number of exponents in family basis sets
      INTEGER MXFAMEXP
      PARAMETER (MXFAMEXP = 100)
      REAL*8 FAMEXP(MXFAMEXP, 2), FAMPAR(4)
      INTEGER NFAMEXP(2)
      COMMON /CBIFAM/ FAMEXP, FAMPAR, NFAMEXP
! --  end of cbirea.h --
