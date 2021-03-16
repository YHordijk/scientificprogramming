!comdeck cgas.inc $Revision$

      INTEGER IDOGAS,NGAS,NGSSH,NGSOBT,NGSOB,IGSOCC,IGSINA,IGSDEL,      &
     &        NCISPC,NCMBSPC,LCMBSPC,ICMBSPC,NMXOCCLS,IPHGAS,IHPVGAS,   &
     &        IGSOCCX
 
      COMMON/CGAS/IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),                    &
     &            NGSOB(MXPOBS,MXPNGAS),                                &
     &            NGSOBT(MXPNGAS),IGSOCC(MXPNGAS,2),IGSINA,IGSDEL,      &
     &            IGSOCCX(MXPNGAS,2,MXPICI),NCISPC,                     &
     &            NCMBSPC, LCMBSPC(MXPICI),ICMBSPC(MXPSTT,MXPICI),      &
     &            NMXOCCLS,IPHGAS(MXPNGAS),IHPVGAS(MXPNGAS)
