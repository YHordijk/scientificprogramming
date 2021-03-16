!
! FILE    : parkrmc.h
!
!                     PARALLEL INFORMATION
!
!             common blocks file for parallel KR-MCSCF
!     ======================================================
!
!
      PARAMETER (IKRMCPARF = 8)
!
!     integer block         
!     -------------
      COMMON/KRMCPARICB/NUMCPU,LMCPRI
!
!     double precision block
!     ----------------------
!     DOUBLE PRECISION  
!     COMMON/KRMCPARRCB/ 
!
!     logical block
!     -------------
      LOGICAL PARMCSCF
      COMMON/KRMCPARLCB/PARMCSCF
!
!     character block
!     ---------------
!     CHARACTER 
!     COMMON/KRMCPARCCB/ 
