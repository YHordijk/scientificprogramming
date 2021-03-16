!
! FILE: include/nuclei.h
!
      REAL*8  CHARGE(MXCENT), CORD(3, MXCENT), GNUEXP(MXCENT)
      INTEGER NUCPRE(MXCENT), NUCNUM(MXCENT, 8), NUCDEG(MXCENT),        &
     &    ISTBNU(MXCENT), NUCIND, NUCDEP, NTRACO, ITRACO(3), NATOMS,    &
     &    NFLOAT, NBASIS, NLARGE, NSMALL, NPBAS, NPLRG, NPSML,          &
     &    NUCNET, INCENT(MXCENT), INUNIQ(MXCENT), NDEGNM(MXCENT),       &
     &    ISOTOP(MXCENT), IZATOM(MXCENT),                               &
     &    NAMFA, IAMFA(MXCENT), IATOMTYP(MXCENT),                       &
     &    NONTYP, NONT(MXATOM), NBAS_ATOMTYPE(MXATOM,0:MXBSETS),        & ! <--- atom type info
     &    NBASIS_FIT, NLARGE_FIT, NSMALL_FIT, NPBAS_FIT, NPLRG_FIT,     &
     &    NPSML_FIT
      LOGICAL NOORBT(MXCENT), GAUNUC

      COMMON /NUCLEI/     CHARGE, CORD,   GNUEXP,                       &
     &    NUCPRE, NUCNUM, NUCDEG, ISTBNU, NUCIND, NUCDEP, NTRACO,       &
     &    ITRACO, NATOMS, NFLOAT, NBASIS, NLARGE, NSMALL,               &
     &    NPBAS,  NPLRG,  NPSML,  NUCNET, INCENT, INUNIQ,               &
     &    NDEGNM, ISOTOP, IZATOM, NAMFA,  IAMFA,  IATOMTYP,             &
     &    NONTYP, NONT        , NBAS_ATOMTYPE            ,              & ! <--- atom type info
     &    NBASIS_FIT, NLARGE_FIT, NSMALL_FIT, NPBAS_FIT, NPLRG_FIT,     &
     &    NPSML_FIT, NOORBT, GAUNUC

      CHARACTER*4 NAMN(MXCENT)
      CHARACTER*6 NAMEX(MXCOOR), NAMDEP(MXCENT), NAMDPX(MXCOOR)
      COMMON /NUCLEC/ NAMN, NAMEX, NAMDEP, NAMDPX
! --- end of nuclei.h ---
