!comdeck strinp.inc $Revision$ 
!./STRINP/

      INTEGER NSTTYP,MNRS1,MXRS1,MNRS3,MXRS3,NELEC,IZORR,IAZTP,IBZTP,    &
     &        IARTP,IBRTP,NZSTTP,NRSTTP,ISTTP                

      COMMON/STRINP/NSTTYP,MNRS1(MXPSTT),MXRS1(MXPSTT),                 &
     &              MNRS3(MXPSTT),MXRS3(MXPSTT),NELEC(MXPSTT),          &
     &              IZORR(MXPSTT),IAZTP,IBZTP,IARTP(3,10),IBRTP(3,10),  &
     &              NZSTTP,NRSTTP,ISTTP(MXPSTT)
