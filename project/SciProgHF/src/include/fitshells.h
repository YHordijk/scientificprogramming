!     Data structure similar to that of shells but now for the
!     set of fit functions.
      LOGICAL SHARE_FIT(MXSHEL_FIT), SEGM_FIT(MXSHEL_FIT)
      LOGICAL SPHR_FIT(MXSHEL_FIT)
      REAL*8 CENT_FIT(MXSHEL_FIT, 3, 2)
      INTEGER NHKT_FIT(MXSHEL_FIT), KHKT_FIT(MXSHEL_FIT),               &
     &        KCKT_FIT(MXSHEL_FIT), ISTBAO_FIT(MXSHEL_FIT),             &
     &        NUCO_FIT(MXSHEL_FIT), JSTRT_FIT(MXSHEL_FIT),              &
     &        NSTRT_FIT(MXSHEL_FIT), NCENT_FIT(MXSHEL_FIT),             &
     &        NRCO_FIT(MXSHEL_FIT), NUMCF_FIT(MXSHEL_FIT),              &
     &        NBCH_FIT(MXSHEL_FIT), KSTRT_FIT(MXCORB_FIT),              &
     &        LCLASS_FIT(MXSHEL_FIT), IPTSHL_FIT(MXSHEL_FIT),           &
     &        NUMCFT_FIT(MXSHEL_FIT), KMAX_FIT, NLRGSH_FIT, NSMLSH_FIT
      COMMON /FITSHELLS/ CENT_FIT, NHKT_FIT, KHKT_FIT, KCKT_FIT,        &
     &                   ISTBAO_FIT, NUCO_FIT, JSTRT_FIT, NSTRT_FIT,    &
     &                   NCENT_FIT, SHARE_FIT, NRCO_FIT, NUMCF_FIT,     &
     &                   NBCH_FIT, KSTRT_FIT, SEGM_FIT, LCLASS_FIT,     &
     &                   IPTSHL_FIT, NUMCFT_FIT, SPHR_FIT, KMAX_FIT,    &
     &                   NLRGSH_FIT, NSMLSH_FIT
