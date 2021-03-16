!
!     CLASS LABELS
!
!     CLSINT - keyword for a class of integrals
!     CLSCMB - indicates what blocks to compute : LL,SL,LS,SS
!     CLSCAL - indicates whether class is active, that is to be calculated
!
      INTEGER MAXCLS
! MI: increase from 15 to 45 
      PARAMETER (MAXCLS = 45)

      LOGICAL CLSCAL(MAXCLS)
      CHARACTER*7 CLSINT(MAXCLS)
      CHARACTER*4 CLSCMB(MAXCLS)
      INTEGER NPRPCLS, IORDCL(MAXCLS)
      COMMON /XCBCLS/ NPRPCLS, IORDCL, CLSCAL

      COMMON /CLSLBL/ CLSINT, CLSCMB
