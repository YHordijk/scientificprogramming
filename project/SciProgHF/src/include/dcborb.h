!
! FILE    : include/dcborb.h
!
! DCBORB  : Dirac Common Block - ORBital informaiton
!
!  NKRBLK : 1 for codes that store MOs as quaternions or
!           2                             complex spinors
!

#include "mxgas.h"

      INTEGER I1_DCBORB,                                                &
     &        NORB(2), N2ORB(2), NORBT, N2ORBT, N2ORBX,                 &
     &                 NBORB(4, 2, 0:2),NFORB(2, 0:2), NTORB(0:2),      &
     &                 N2ORBTQ,  N2ORBXQ,                               &
     &        NISH(2), N2ISH(2), NISHT, N2ISHT, N2ISHX,                 &
     &        NASH(2), N2ASH(2), NASHT, N2ASHT, N2ASHX,                 &
     &                 N2ASHTQ,  N2ASHXQ, NNASHX,                       &
     &        NOCC(2), N2OCC(2), NOCCT, N2OCCT, N2OCCX,                 &
     &        NSSH(2), N2SSH(2), NSSHT, N2SSHT, N2SSHX,                 &
     &        NESH(2), N2ESH(2), NESHT, N2ESHT, N2ESHX,                 &
     &        NPSH(2), N2PSH(2), NPSHT, N2PSHT, N2PSHX,                 &
     &        NFRO(2),           NFROT,                                 &
     &        IORB(2), I2ORBT(2), I2ORBX(2, 2),                         &
     &        IASH(2), I2ASHT(2), I2ASHX(2, 2),                         &
     &        IOCC(2), I2OCCT(2), I2OCCX(2, 2),                         &
     &        NCMO(2), NCMOQ(2), NCMOT, NCMOTQ, ICMO(2), ICMOQ(2),      &
     &        NAOCC(2), NAOCCT, NIOCC(2), NIOCCT,                       &
     &        NAVIR(2), NAVIRT, NIVIR(2), NIVIRT,                       &
     &        NELEC(2), NAELEC, NELECT, NKRBLK,                         &
     &        NGAS_DC, NGSH(2, MXGAS),NGSHT(MXGAS), NGASSP(2, MXGAS),   &  ! NGAS_DC to avoid conflict with NGAS in LUCI* modules  /SK-Aug 2008
     &        NPAS_DC, NPASSP(MXGAS, MXPAS),                            &
     &        I_DCBORB_SET, I2_DCBORB
!
!     IMPORTANT !!!
!     I1_DCBORB and I2_DCBORB are tags used to calculate the size of
!     this COMMON block; must remain at start/end 
!
      COMMON /DCBORB/  I1_DCBORB,                                       &
     &        NORB   , N2ORB   , NORBT, N2ORBT, N2ORBX,                 &
     &                 NBORB           ,NFORB        , NTORB     ,      &
     &                 N2ORBTQ,  N2ORBXQ,                               &
     &        NISH   , N2ISH   , NISHT, N2ISHT, N2ISHX,                 &
     &        NASH   , N2ASH   , NASHT, N2ASHT, N2ASHX,                 &
     &                 N2ASHTQ,  N2ASHXQ, NNASHX,                       &
     &        NOCC   , N2OCC   , NOCCT, N2OCCT, N2OCCX,                 &
     &        NSSH   , N2SSH   , NSSHT, N2SSHT, N2SSHX,                 &
     &        NESH   , N2ESH   , NESHT, N2ESHT, N2ESHX,                 &
     &        NPSH   , N2PSH   , NPSHT, N2PSHT, N2PSHX,                 &
     &        NFRO   ,           NFROT,                                 &
     &        IORB   , I2ORBT   , I2ORBX      ,                         &
     &        IASH   , I2ASHT   , I2ASHX      ,                         &
     &        IOCC   , I2OCCT   , I2OCCX      ,                         &
     &        NCMO   , NCMOQ   , NCMOT, NCMOTQ, ICMO   , ICMOQ   ,      &
     &        NAOCC   , NAOCCT, NIOCC   , NIOCCT,                       &
     &        NAVIR   , NAVIRT, NIVIR   , NIVIRT,                       &
     &        NELEC   , NAELEC, NELECT, NKRBLK,                         &
     &        NGAS_DC, NGSH          ,NGSHT       , NGASSP          ,   &
     &        NPAS_DC, NPASSP              ,                            &
     &        I_DCBORB_SET, I2_DCBORB         ! Do NOT place any variables after the end tag "I2_DCBORB" !!!

!     Information on orbital subblocks. Used at the moment only for linear symmetry
!     but can be generalized to eliminate couplings between orbital spaces.

      INTEGER MAX_SUB_BL
      PARAMETER (MAX_SUB_BL = 72)

      INTEGER N_SUB_BL(2), ID_SUB_BL(MAX_SUB_BL, 2),                    &
     &        NORB_SUB(MAX_SUB_BL,2,0:2),NTMO_SUB(MAX_SUB_BL,2,0:2)
      COMMON /DCBORB_SUB_BL/ N_SUB_BL, ID_SUB_BL, NORB_SUB, NTMO_SUB

! -- end of dcborb.h --
