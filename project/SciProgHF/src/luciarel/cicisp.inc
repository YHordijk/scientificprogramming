! 
!     CI space information
!
      integer :: NICISP, MXSB, MXSOOB, MXNTTS, LCJRES, MXADKBLK
      integer :: ISYM_KRMC, LCJRES_SAVE
      integer :: ISMOST
      real(8) :: XISPSM

      COMMON/CICISPI/NICISP,ISMOST(MXNDGIRR,MXNDGIRR),MXSB,MXSOOB,      &
     &               MXNTTS,LCJRES,MXADKBLK,ISYM_KRMC,LCJRES_SAVE

      COMMON/CICISPR/XISPSM(MXNDGIRR,MXPICI)
