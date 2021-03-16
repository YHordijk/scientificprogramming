!
! File: optinf.h
!
! Information for geometry optimization
! in abaopt.F, abaop2.F, and abarint.F
!
      REAL*8  TRSTRA, TRSTIN, TRSTDE, RTENBD, RTENGD, RTRJMN, RTRJMX,   &
     &        ENERGY, ERGOLD, ERGPRD, ERGPRO, STPNRM, STPNRO, GRADNM,   &
     &        THRERG, GRDTHR, THRSTP, THRSYM, EVLINI, DISPLA,           &
     &        STPDIA(8*MXCENT), STPSYM(8*MXCENT), STPINT(8*MXCENT),     &
     &        GRDDIA(8*MXCENT), EVAL(8*MXCENT), EVALOL(8*MXCENT),       &
     &        GRDINT(8*MXCENT), CRDIN1(8*MXCENT), CRDINT(8*MXCENT),     &
     &        CNDHES(0:7)
      INTEGER INDHES(0:7), INTCRD(8*MXCENT, 6), ICONF(0:5),             &
     &        ICNSTR(8*MXCENT), IAUXRD, ITOTRJ, KEPTIT, NSPMOD, NCNSTP, &
     &        INDTOT, ITRNMR, ITRMAX, MAXREJ, IPRINT, NCRTOT, NCART,    &
     &        NPROJ, NTMAT, IINTCR, IREDIC, ICRTCR, ICONDI, ITRBRK,     &
     &        NUMPRE, IPRE, NENECALCS
      LOGICAL GECONV, NOTRST, NOBRKS, BRKSYM, NWSYMM, DOSPE, DOPRE,     &
     &        FINPRE, VRML, VRBOND, VREIGV, VRCORD, VRVIBA, VRML_SYM,   &
     &        VISUAL, BOFILL, INITHS, HSFILE, BFGSR1, STEEPD, RANKON,   &
     &        PSB, DFP, BFGS, SCHLEG, NEWTON, QUADSD, KEEPHE, BAKER,    &
     &        REDINT, CARTCO, INRDHS, FSTORD, SNDORD, REJINI, GRDINI,   &
     &        MULTI, CHGRDT, CONOPT, MODHES, INMDHS, FINDRE, TRSTRG,    &
     &        RATFUN, GDIIS, DELINT, RSTARR, LNSRCH, SADDLE, REBILD,    &
     &        CMBMOD, HFPROP, CONFRM, NOAUX, NODIHE, LINDHD
      COMMON /OPTINF/ TRSTRA, TRSTIN, TRSTDE, RTENBD, RTENGD, RTRJMN,   & ! we start with double precision variables
     &                RTRJMX, ENERGY, ERGOLD, ERGPRD, ERGPRO, STPNRM,   &
     &                STPNRO, GRADNM, THRERG, GRDTHR, THRSTP, THRSYM,   &
     &                EVLINI, DISPLA, STPDIA, STPSYM, STPINT, GRDDIA,   &
     &                EVAL, EVALOL, GRDINT, CRDIN1, CRDINT, CNDHES,     &
     &                INDHES, INTCRD, ICONF, ICNSTR, IAUXRD, ITOTRJ,    & ! first line with integer variables
     &                KEPTIT, NSPMOD, NCNSTP, INDTOT, ITRNMR, ITRMAX,   &
     &                MAXREJ, IPRINT, NCRTOT, NCART, NPROJ, NTMAT,      &
     &                IINTCR, IREDIC, ICRTCR, ICONDI, ITRBRK, NUMPRE,   &
     &                IPRE, NENECALCS,                                  &
     &                GECONV, NOTRST, NOBRKS,                           & ! first line with logical variables
     &                BRKSYM, NWSYMM, DOSPE, DOPRE, FINPRE, VRML,       &
     &                VRBOND, VREIGV, VRCORD, VRVIBA, VRML_SYM, VISUAL, &
     &                INITHS, HSFILE, BFGSR1, STEEPD, RANKON, PSB, DFP, &
     &                BFGS, SCHLEG, NEWTON, QUADSD, KEEPHE, BAKER,      &
     &                REDINT, CARTCO, INRDHS, FSTORD, SNDORD, REJINI,   &
     &                GRDINI, MULTI, CHGRDT, CONOPT, MODHES, INMDHS,    &
     &                FINDRE, TRSTRG, RATFUN, GDIIS, DELINT, RSTARR,    &
     &                LNSRCH, SADDLE, REBILD, BOFILL, CMBMOD, HFPROP,   &
     &                CONFRM, NOAUX, NODIHE, LINDHD

      INTEGER MAXPRE
      PARAMETER (MAXPRE = 10)
      CHARACTER*80            PREBTX,         SPBSTX
      COMMON /OPTINF_C/       PREBTX(MAXPRE), SPBSTX
! --- end of optinf.h ---
