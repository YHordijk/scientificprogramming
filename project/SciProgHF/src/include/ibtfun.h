!
! FILE:  ibtfun.h
!
      INTEGER IBTAND,IBTOR,IBTSHL,IBTSHR,IBTXOR
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = IEOR(I,J)
