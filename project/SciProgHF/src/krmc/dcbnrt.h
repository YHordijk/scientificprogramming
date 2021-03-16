C -- FILE : dcbnrt.h
!
! FILE : dcbnrt.h
!
! For freezing ("NO ROTATION") orbitals in MCSCF
!
!     Dependencies: maxorb.h
!
      LOGICAL LNOROT
      COMMON /DCINRT/ NOROT(MXCORB), NOROTC(MXCORB)
      COMMON /DCLNRT/ LNOROT
! -- end of dcbnrt.h --
