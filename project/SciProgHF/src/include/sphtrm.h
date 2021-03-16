! -- file sphtrm.h --
!  CSP : Cartesian->Spherical transformation matrices
!  CSP( ISPADR(L+1) ) : first element of transformation matrix for ang.mom. = L
!  -- CSP and ISPADR are set in SPHINP(herrdn.F)
!
      INTEGER NCSP
      PARAMETER (NCSP =                                                 &
     &           ((MXQN+2)*(MXQN+1)*(3*MXQN**2+6*MXQN+1)*MXQN)/60)

      REAL*8  CSP(NCSP)
      INTEGER ISPADR(MXQN)
      COMMON /SPHTRM/ CSP, ISPADR
! -- end of sphtrm.h --
