!comdeck gasstr.inc $Revision$
      

      INTEGER MNGSOC,MXGSOC,NGPSTR,                                      &
     &              IBGPSTR,NELFGP,IGSFGP,                              &
     &              NSTFGP,MNELFGP,MXELFGP,                             &
     &              NELFTP,NSPGPFTP,IBSPGPFTP,                          &
     &              ISPGPFTP,NELFSPGP,                                  &
     &              NSTFSMSPGP,                                         &
     &              NGRP,NSTTP,MXNSTR,NTSPGP,MAX_STR_SPGP,              &
     &              MAX_STR_OC_BLK,MXSMCLS,MXSMCLSE,MXSMCLSE1,          &
     &              NHLFSPGP,MNHL,NSTFSMGP,                             &
     &              MINMAX_SM_GP


      COMMON/GASSTR/MNGSOC(MXPNGAS),MXGSOC(MXPNGAS),NGPSTR(MXPNGAS),    &
     &              IBGPSTR(MXPNGAS),NELFGP(MXPSTT),IGSFGP(MXPSTT),     &
     &              NSTFGP(MXPSTT),MNELFGP(MXPNGAS),MXELFGP(MXPNGAS),   &
     &              NELFTP(MXPSTT),NSPGPFTP(MXPSTT),IBSPGPFTP(MXPSTT),  &
     &              ISPGPFTP(MXPNGAS,MXPSTT),NELFSPGP(MXPNGAS,MXPSTT),  &
     &              NSTFSMSPGP(MXPNSMST,MXPSTT),                        &
     &              NGRP,NSTTP,MXNSTR,NTSPGP,MAX_STR_SPGP,              &
     &              MAX_STR_OC_BLK,MXSMCLS,MXSMCLSE,MXSMCLSE1,          &
     &              NHLFSPGP(MXPSTT),MNHL,NSTFSMGP(MXPNSMST,MXPSTT),    &
     &              MINMAX_SM_GP(2,MXPSTT)
