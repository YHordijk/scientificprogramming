!*** Information for 4-index transformation program **
!
! Adapted for include-into-Fortran90 form by Miro Ilias
!
!     Length of common block : keep up to date !!
      LOGICAL TRA_ANTIS,TRA_EXACORR,TRA_MSOUT
      LOGICAL NOPAIR,NO2IND,NO4IND,PRPTRA,RCORBS
      LOGICAL FKRMC, MDCSCAT, TRA_ASCII, PRPSYA
      LOGICAL NOMDCINT, NOMDCINT_TRAINP

!     aspg, 20090424, for fine-sort of 1HT
      logical do_finesort,set_auto_finesort
!
      CHARACTER*72 TRA_INDSTR,TRA2_INDSTR,TRA4_INDSTR,TRA_CORSTR,       &
     &             TRA2_INDPRP,TRA_CORSTR2
      CHARACTER*7  MOFILE_TRAINP

! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!     NDCBTRI, NDCBTRR, NDCBTRL are *4 lengths of the three
!     common blocks DCBTRI, DCBTRR, DCBTRL. They are used
!     to transfer the common blocks to slaves when parallel
!     with MPI.
!     NBNBNB! THEY MUST BE UPDATED IF YOU CHANGE THE COMMON BLOCKS!!!!
      INTEGER NDCBTRI,NDCBTRR,NDCBTRL
      PARAMETER (NDCBTRI=284,NDCBTRR=259,NDCBTRL=16)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      COMMON/DCBTC/TRA_INDSTR(2),TRA2_INDSTR(2,2),TRA4_INDSTR(4,2),     &
     &             TRA_CORSTR(2), TRA2_INDPRP(2), MOFILE_TRAINP,        &
     &             TRA_CORSTR2(2)

      INTEGER IPRTRA,ISTRAT,ITRA_INTFLG,ITRA_INTFL2,ITRA_INTFL4,        &
     &              NQQCLASS,ISAME,ICLARR,                              &
     &              NFPCK12,NFPCK34,NFPCK12T,NFPCK34T,                  &
     &              IFPCK12,IFPCK34,N2EFIL,NTKL,IPAR4BS,NCORE2
      COMMON/DCBTRI/IPRTRA,ISTRAT,ITRA_INTFLG,ITRA_INTFL2,ITRA_INTFL4,  &
     &              NQQCLASS,ISAME(4),ICLARR(4,4,16),                   &
     &              NFPCK12(2),NFPCK34(2),NFPCK12T,NFPCK34T,            &
     &              IFPCK12(2,2),IFPCK34(2,2),N2EFIL,NTKL,IPAR4BS,      &
     &              NCORE2

      COMMON/DCBTRL/TRA_ANTIS,TRA_EXACORR,TRA_MSOUT,                              &
     &              NOPAIR,NO2IND,NO4IND,PRPTRA,RCORBS,                 &
     &              FKRMC, MDCSCAT, TRA_ASCII, PRPSYA,                  &
     &              NOMDCINT,NOMDCINT_TRAINP,                           &
     &              do_finesort,set_auto_finesort


      REAL*8 THROUT,SCRTRA,RINVICLARR,TH_CHOLESKY
      COMMON/DCBTRR/THROUT,SCRTRA,RINVICLARR(4,4,16),TH_CHOLESKY

!     Unit for ASCII interface file
      INTEGER LUASCII
      PARAMETER (LUASCII=74)
!     Units for MOLFDIR interface files
      INTEGER LUMLF1,LUMLF2,LUMLF3
      PARAMETER (LUMLF1=75,LUMLF2=76,LUMLF3=76)
!     Units for scratch files (LUTRA1 is an offset (-1) for 2 files !)
      INTEGER LUTRA1,LUTRA2
      PARAMETER (LUTRA1=78,LUTRA2=81)
#include <comdis.h>

! -- end of dcbtra.h --
