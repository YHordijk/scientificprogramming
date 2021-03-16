!
! FILE    : gencon.h
!
!   - Parameter for maximum number of functions in a block:
!     MXCON_PARAM = NPRIM * l(l+1)/2 * NCONTR
!   - DIRCON true if contracted basis set is used
!
      INTEGER MXCON_PARAM
      PARAMETER (MXCON_PARAM = 300)

      LOGICAL DIRCON(4)
      COMMON /GENCON/ DIRCON
