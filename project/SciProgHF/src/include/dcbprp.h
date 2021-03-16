! File: include/dcbprp.h
!
! Purpose: control information for property module (in .../src/prp/ )
!*****************************************************************************
! *** Property module ***
!  DOEXP - do expectation values
!  DOESR - do ESR properties
!  DOXLR - do linear response
!  DOXPP - do excitation energies 
!          (by the polarization propagator method)
!  DOXQR - do quadratic response
!  DOTPA - do two-photon absortion calculation
!  DOEXCPRP  - do excited state properties calculation
!  EXCPRPTYP - type of excited state property calculation
!  DOVER - do verdet calculations (via DOXQR)
!  MGRAD - calculate molecular gradient
!  DOPCT - do perform picture change transformation of four-component
!          property operators
!  NOPCT - do NOT perform the picture change transformation when DOPCT is on
!          - take only LL block of the four-component property operator
           
      INTEGER LUCORP
      PARAMETER (LUCORP = 91)

      LOGICAL DOEXP, DOXLR, DOXPP, DOXQR, DIPOLE, EPOLAR, MSUSCP,       &
     &        SHIELD, SPNSPN, EFG, NQCC, MGRAD, NSTDIA, PVC, DOESR,     &
     &        ESRGTENS, ESR_HFCC, EFT, EFTNTL, EFN, RHONUC,             &
     &        DOTPA,DOVER,DOEXCPRP,DOVDW,REGRAD,                        &
     &        ROTG, LONDON, SYMCON, NOORTH, DOPCT, MSUSCDIA, NOPCT,     & 
     &        PVC_SHIELD, PVC_SPINSPIN, DOSTEX, OPTROT,                 &
     &        BXLAO, BYLAO, BZLAO, EFFDEN, EFFDE2,SEPTEP,FULMAT,EXPPED, &
     &        DONMQM,SPINRO
      COMMON /CBLPRP/ DOXLR, DOXPP, DOEXP, DOXQR, DIPOLE, EPOLAR,DOVDW, &
     &                MSUSCP, SHIELD, SPNSPN, EFG, NQCC, MGRAD, NSTDIA, &
     &                PVC, DOESR, ESRGTENS, ESR_HFCC, EFT, EFTNTL,      &
     &                EFN, RHONUC, DOTPA, DOVER, DOEXCPRP,REGRAD,       &
     &                ROTG,LONDON,SYMCON,NOORTH,DOPCT,MSUSCDIA,NOPCT,   & 
     &                PVC_SHIELD,PVC_SPINSPIN,DOSTEX,OPTROT,            &
     &                BXLAO, BYLAO, BZLAO, EFFDEN, EFFDE2,SEPTEP,FULMAT,&
     &                EXPPED, DONMQM,SPINRO

      INTEGER IPRPRP, IPEPOLAR(3, 2), IPSPNSPN(2, MXCOOR),              &
     &        IPDSO(MXCOOR*(MXCOOR+1)/2), NDSOIN, IPSHIELD(MXCOOR+3),   &
     &        IPNSTDIA(MXCOOR*3), NNSTIN, IPESR(MXCOOR+3, 2), NPRP_WF,  &
     &        IPROTG(9),IPMSUSC(2,3),IMAGOUT, EXCPRPTYP,                &
     &        IPCON21(3),IPCON22(3),IPCON23(6),                         &
     &        IPLONDON(4,3),IPCON(3),IPRNST(MXCOOR*3),                  &
     &        IP_PVC_SHIELD(MXCOOR+3),                                  &
     &        IP_PVC_SPINSPIN(2,MXCOOR),                                &
     &        IP_BLAO(2,3),IPSPINRO(MXCOOR+6)
      COMMON /CBIPRP/ IPRPRP, IPEPOLAR, IPSPNSPN, IPDSO, NDSOIN,        &
     &                IPSHIELD, IPNSTDIA, NNSTIN, IPESR, NPRP_WF,       &
     &                         IPROTG,IPMSUSC,IMAGOUT, EXCPRPTYP,       &
     &                IPCON21,IPCON22,IPCON23,                          &
     &                IPLONDON,IPCON,IPRNST,                            &
     &                IP_PVC_SHIELD,                                    &
     &                IP_PVC_SPINSPIN,                                  &
     &                IP_BLAO,IPSPINRO

      REAL*8 ABUND
      COMMON /CBRPRP/ ABUND

      INTEGER MXPRP_WF
      PARAMETER (MXPRP_WF = 4)
      CHARACTER*4 PRP_WF(MXPRP_WF)
      COMMON /CBCPRP/ PRP_WF

#include "cbiher.h"
! end of dcbprp.h
