!
! FILE    : dcbkrmc.h
!
!*** Information for KR-MCSCF calculations ***

      INTEGER NMCFLAG
      PARAMETER (NMCFLAG = 20)

      LOGICAL MCFLAG(NMCFLAG), KRMC_SKIPEE, KRMC_SKIPEP
      LOGICAL KRMC_NOCI, KRMC_UCIBOS, KRMC_CHCKJZ
      LOGICAL KRMC_NO1pdens, krmc_save_reordered_nos, KRMC_full_ci
      COMMON /DCLKRMC/ MCFLAG, KRMC_SKIPEE, KRMC_SKIPEP, KRMC_NOCI,     &
     &                 KRMC_UCIBOS, KRMC_CHCKJZ,                        &
     &                 KRMC_NO1pdens, krmc_save_reordered_nos,          &
     &                 KRMC_full_ci

      INTEGER NKRMCISH(2), NKRMCASH(2), NKRMCSSH(2),                    &
     &        NKRMCGSH(2, MXGAS), NKRMCGSP(2, MXGAS), NKRMCESH(2),      &
     &        NKRMCPSH(2), NKRMCFRO(2), NKRMCAELEC, NKRMCGAS,           &
     &        NKRMC_MK2REF, NKRMC_MK2DEL, NKRMC_MINMK2, NKRMC_MAXMK2
      COMMON /DCOKRMC/ NKRMCISH, NKRMCASH, NKRMCSSH, NKRMCGSH, NKRMCGSP,&
     &                 NKRMCESH, NKRMCPSH, NKRMCFRO, NKRMCAELEC,        &
     &                 NKRMCGAS, NKRMC_MK2REF, NKRMC_MK2DEL,            &
     &                 NKRMC_MINMK2, NKRMC_MAXMK2
  
      INTEGER IKRMC_ITRINT(2), IKRMC_INTDEF, IKRMC_INTFLG, IKRMC_INTBUF,&
     &        MKRMC_MXMICRO, MKRMC_MXMACRO, MKRMC_MAXBCK, IKRMC_STATE,  &
     &        IKRMC_SYMMETRY, IPRKRMC, IKRMC_MEMFAC, IKRMC_SVRONO(2)
      COMMON /DCIKRMC/ IKRMC_ITRINT, IKRMC_INTDEF, IKRMC_INTFLG,        &
     &                 IKRMC_INTBUF, MKRMC_MXMICRO, MKRMC_MXMACRO,      &
     &                 MKRMC_MAXBCK, IKRMC_STATE, IKRMC_SYMMETRY,       &
     &                 IPRKRMC, IKRMC_MEMFAC, IKRMC_SVRONO

      REAL*8 DKRMC_CNVINT(2), DKRMC_THRGRD, DKRMC_THRPCI, DKRMC_MVOFAC
      COMMON /DCRKRMC/ DKRMC_CNVINT, DKRMC_THRGRD, DKRMC_THRPCI,        &
     &                 DKRMC_MVOFAC

      CHARACTER*8 KRMC_CIPROGRAM
      CHARACTER*72 KRMC_FRZSTR(2), KRMC_DELSTR(2)
      COMMON /DCCKRMC/ KRMC_CIPROGRAM, KRMC_FRZSTR, KRMC_DELSTR
! -- end of dcbkrmc.h --
