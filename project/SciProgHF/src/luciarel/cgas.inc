      integer NACTEL, IDOGAS, NGAS, IGSINA, IGSDEL, NCISPC, NCMBSPC
      integer NMXOCCLS
      integer NGSSH, NGSSH2
      integer IGSOCC, IGSOCCX
      integer LCMBSPC, ICMBSPC
      integer MXNELGS
      COMMON/CGASR/NACTEL,IDOGAS,NGAS,IGSINA,IGSDEL,NCISPC,NCMBSPC,     &
     &             NMXOCCLS,                                            &
     &             NGSSH(MXNDGIRR,MXPNGAS),NGSSH2(MXNDGIRR,MXPNGAS),    &
     &             IGSOCC(MXPNGAS,2),IGSOCCX(MXPNGAS,2,MXPICI),         &
     &             LCMBSPC(MXPICI),ICMBSPC(MXPSTT,MXPICI),              &
     &             MXNELGS(MXPNGAS)
