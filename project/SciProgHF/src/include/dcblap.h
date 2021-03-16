!
! FILE    : ./include/dcblap.h
!
!  DCBLAP:  Parameters for Laplace module
!
!
      LOGICAL :: INILAP, FNDLAP
      INTEGER :: NLAP, MXITER, IPRLAP
      REAL(8) :: STEPMX, TOLRNG, TOLPAR, TOLLAP, XPNTS(MXLAP),          &
     &           WGHTS(MXLAP)
!
      COMMON /DCBLAP/ STEPMX, TOLRNG, TOLPAR, TOLLAP, XPNTS, WGHTS,     &
     &                NLAP, MXITER, IPRLAP,                             &
     &                INILAP, FNDLAP
