      integer    MAXSUB
      PARAMETER (MAXSUB = 20)
      LOGICAL DIPOL, SKIP, DOVCD, HESFIL, HESPUN, ROACID
      real(8)         SCALE, COEF
      integer         NCARD,  NINTCM, IPRINT, MAXDIF, NISOTP,           &
     &                iatoms, isotp
      COMMON /CBIVIB/ SCALE(MXCOOR), COEF(MXCOOR), IATOMS(4,MXCOOR),    &
     &                NCARD,  NINTCM, DIPOL,  IPRINT, MAXDIF, NISOTP,   &
     &                ISOTP(MAXSUB,MXCENT), SKIP, DOVCD, HESFIL,        &
     &                HESPUN, ROACID
      CHARACTER KWORD*1, ITYPCD*4, ITYPCM*4
      COMMON /CBIVBC/ KWORD(MXCOOR), ITYPCD(MXCOOR), ITYPCM(MXCOOR)
