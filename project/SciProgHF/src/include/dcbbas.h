!*** Dimensions for AO/SO (non-orthomormal) basis ***
!     Length of common block : keep up to date !!
      INTEGER NDCBBAS
      PARAMETER (NDCBBAS = 338)

      INTEGER I1_DCBBAS, NBBAS(0:7, 0:2), NFBAS(2, 0:2), NTBAS(0:2),    &
     &        NNBBAST, N2BBAST, NNBBASX, N2BBASX, N2BBASXQ, MXBBAS,     &
     &        MXFBAS, N2BAS(2), N2BASQ(2), N2BAST, N2BASTQ,             &
     &        IBBAS(0:7, 0:2), I2BBASX(0:7, 0:7, 2, 2), IBAS(2),        &
     &        I2BASX(2, 2), I2BAST(2), NZT, I2TMT(2), N2TMT, I2_DCBBAS
      COMMON /DCBBAS/ I1_DCBBAS, NBBAS, NFBAS, NTBAS, NNBBAST, N2BBAST, &
     &                NNBBASX, N2BBASX, N2BBASXQ, MXBBAS, MXFBAS, N2BAS,&
     &                N2BASQ, N2BAST, N2BASTQ, IBBAS, I2BBASX, IBAS,    &
     &                I2BASX, I2BAST, NZT, I2TMT, N2TMT, I2_DCBBAS

!    Do NOT add any variables after the end tag !!!

