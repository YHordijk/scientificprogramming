!
! .. some description 
!
!   PRIEXP  - exponents of primitive shells
!   PRICCF  - normalized contraction coefficients
!   PRICRX  - x-coordinate of center
!   PRICRY  - y-coordiante of center
!   PRICRZ  - z-coordinate of center
!
!
      REAL*8 PRIEXP(MXPRIM), PRICCF(MXPRIM, MXCONT), PRICRX(MXPRIM),    &
     &       PRICRY(MXPRIM), PRICRZ(MXPRIM)

      INTEGER NPSHEL, NPLSH, NPSSH

      COMMON /PRIMIT/ PRIEXP, PRICCF, PRICRX, PRICRY, PRICRZ, NPSHEL,   &
     &                NPLSH, NPSSH
