!* File: gasstr.inc
!*  MAX_NSYMDIST added August 2003


      INTEGER MNGSOC,MXGSOC,NGPSTR,                                     &
     &        IBGPSTR,NELFGP,IGSFGP,                                    &
     &        NSTFGP,MNELFGP,MXELFGP,                                   &
     &        NELFTP,NSPGPFTP,IBSPGPFTP,                                &
     &        ISPGPFTP,NELFSPGP,                                        &
     &        NSTFSMSPGP,                                               &
     &        NGRP,NSTTP,MXNSTR,NTSPGP,MAX_STR_SPGP,                    &
     &        MAX_STR_OC_BLK,NELFTP_REF,                                &
     &        ITPM1,ITPM2,                                              &
     &        IST_FOR_DT,NSTFSMSPGP2,                                   &
     &        NSTFSMGP,NSTFSMGP2,                                       &
     &        ISTFSMGP,ISTFSMGP2,                                       &
     &        MINMAX_SM_GP,MINMAX_SM_GP2,                               &
     &        MAX_NSYMDIST

      COMMON/GASSTRR/MNGSOC(MXPNGAS),MXGSOC(MXPNGAS),NGPSTR(MXPNGAS),   &
     &            IBGPSTR(MXPNGAS),NELFGP(MXPSTT),IGSFGP(MXPSTT),       &
     &            NSTFGP(MXPSTT),MNELFGP(MXPNGAS),MXELFGP(MXPNGAS),     &
     &            NELFTP(MXPSTT),NSPGPFTP(MXPSTT),IBSPGPFTP(MXPSTT),    &
     &            ISPGPFTP(MXPNGAS,MXPSTT),NELFSPGP(MXPNGAS,MXPSTT),    &
     &            NSTFSMSPGP(MXNDGIRR,MXPSTT),                          &
     &            NGRP,NSTTP,MXNSTR,NTSPGP,MAX_STR_SPGP,                &
     &            MAX_STR_OC_BLK,NELFTP_REF(MXPSTT),                    &
     &            ITPM1(2,MXPSTT),ITPM2(2,MXPSTT),                      &
     &            IST_FOR_DT(2,MXPSTT),NSTFSMSPGP2(MXNDGIRR,MXPSTT),    &
     &            NSTFSMGP(MXNDGIRR,MXPSTT),NSTFSMGP2(MXNDGIRR,MXPSTT), &
     &            ISTFSMGP(MXNDGIRR,MXPSTT),ISTFSMGP2(MXNDGIRR,MXPSTT), &
     &            MINMAX_SM_GP(2,MXPSTT),MINMAX_SM_GP2(2,MXPSTT),       &
     &            MAX_NSYMDIST
