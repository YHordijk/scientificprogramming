!     Common block for LUCITA input information
!     read in dirrdn.F/PSIINP/LUCTINP
!     Information is transferred to LUCITA routine dirluct.F
!     for further processing.
!
      INTEGER NTABLE
      PARAMETER (NTABLE = 30)

      INTEGER MXNGAS
      PARAMETER (MXNGAS = 16)

      CHARACTER*72 WAFFCD, CALCTP, SZCALD, TITLUC, CRDINA, CRDCCLR,     &
     &             CRDGAS(MXNGAS), CRDGOC(MXNGAS), CRDFRO, CRDRS1,      &
     &             CRDRS2, CRDRS3, CALCTYP(MXNGAS*4),IUSEQ,UCMBCC,      &
     &             UNEWCCV, URES_CC

      COMMON /LUCTINFC/ CRDINA, CRDCCLR, CRDGAS, CRDGOC, CRDFRO,        &
     &                  CRDRS1, CRDRS2, CRDRS3, WAFFCD, CALCTP,         &
     &                  SZCALD, TITLUC,                                 &
     &                  CALCTYP,IUSEQ,UCMBCC,UNEWCCV,URES_CC

      INTEGER IMOKW, NROOTD, ISSYMD, NACTED, IMULTD, IPRNGD, IPRNLD,    &
     &        IDENSD, MXHL1D, MXEL3D, INGASD, NSEQCD, IRSTLT, MXCIVE,   &
     &        IPARMODEL, ICIMAXITER, IMAXBLKSIZE,                       &
     &        I_USE_DIST_ROUTE, IN_MEMFAC
      COMMON /LUCTINFI/IMOKW(NTABLE), NROOTD, ISSYMD, NACTED, IMULTD,   &
     &                 IPRNGD,IPRNLD,IDENSD,MXHL1D,MXEL3D,INGASD,NSEQCD,&
     &                 IRSTLT, MXCIVE, IPARMODEL, ICIMAXITER,           &
     &                 IMAXBLKSIZE, I_USE_DIST_ROUTE, IN_MEMFAC
      REAL*8 CTRUNC_FAC, my_convergence
	COMMON /LUCTINFR/ CTRUNC_FAC, my_convergence
