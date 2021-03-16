#if defined (VAR_STAR2)
      INTEGER*2 JWOP, KLWOP
#endif
C     MAXWOP = maximum number of rotations (dimension of JWOP)
C     MAXOCC = maximum number of occupied orbitals
#if defined (VAR_SIRBIG)
      PARAMETER ( MAXWOP = 30 000 , MAXOCC = 120 )
#else
      PARAMETER ( MAXWOP = 5 000 , MAXOCC = 120 )
#endif
      COMMON /INFVAR/ NCONF,NWOPT,NVAR,JWOPSY,NWOP(8),NWOPH,NVARH,
     *                JWOP(2,MAXWOP),KLWOP(MAXOCC,MXCORB)
