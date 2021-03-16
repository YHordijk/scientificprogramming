!      Copyright (c) 2019 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org

!*******************************************************************
!
!     RECP : RECP integral routine
!
!     This program is based on argos integral routine
!     with the permission of R. M. Pitzer
!     and maintained by Yoon Sup Lee and Young Choon Park (KAIST)
!
!*******************************************************************
      SUBROUTINE RECP

      USE RECP_OUTPUT
      USE RECP_CHK
      USE RECP_INP_READ
      IMPLICIT NONE
#include "priunit.h"
#include "argoscom.h"
#include "inc_print.h"

!     Enter RECP
      CALL QENTER('RECP')
      CALL TITLER('RECP integral start','*',126)

!     Parameters from menu (RECP_LNK_RDMENU)
      PRINT_LEVEL = ECPINT_PRINT ! .PRINT
      RECP_DBG    = ECPINT_DEBUG ! .DEBUG
      RECP_OUT    = LUPRI        ! DIRAC output

      CALL OUTPUT_TITLE   ! Print title
      CALL RECP_CHK_MAIN  ! Check consistency
      CALL RECP_1E        ! Main part

!     Exit RECP
      CALL TITLER('RECP integral end','*',126)
      CALL QEXIT('RECP')
      END

! ------------------------------------------------

      SUBROUTINE RECP_1E 
!  read the user input and compute the 1-e integrals.
      USE RECP_CSO_TABLE
      USE RECP_INP
      USE RECP_IPT
      USE RECP_NTR
      USE RECP_ONE
      USE RECP_OUTPUT
      USE RECP_SOCFPD
      implicit logical(a-z)
#include "inc_mxvalue.h"
#include "inc_print.h" 
      integer
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      common /parmi/
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      integer
     & mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,inam,nnam,mdum
      common /ntgr/
     & mconsu,mru,mcru,msfu,msfru,ngcs,nu,   
     & lxyzir(3),inam(5),nnam,mdum(32)
      integer
     & ica(mcup,msup,mgup),icb(4,24,mccup),icxst(2,mnsfup),
     & iprst(mnsfup),lmnp1(mconsp),mcons(msfup),nc(msup),
     & ncon(mconsp),nf(msup),nfct(mccup),ngw(mccup),npair(2,mnsfup),
     & nprir(2,mstup,mnsfup),nrcr(mconsp),nt(msfup),ntl(msfup),
     & ntu(msfup)
      real*8
     & zet(mconup,mconsp), eta(mrcrup,mconup,mconsp),
     & x(mcup,msup), y(mcup,msup), z(mcup,msup)
!     # local:
      integer
     & lcr(msup),       lls(msup),      nkcrl(6,msup), nkcru(6,msup),
     & nklsl(4,msup),   nklsu(4,msup),  mcrs(msup),    nir(kaordp),
     & maords(mgcsup),  nd(mstup),     nso(mstup),
     & nsopr(mstup),    nblpr(mstup),   nopir(mstup),  ilxyz(3,mstup),
     & la(mrup,kaordp), ncr(mcrup),     nct(msfup),    mgcs(msfup),
     & lb(msfrup),      ms(msfrup),     mnl(msfrup),   mau(mnrup),
     & icxast(mnsfup),
     & idp(mstup,mstup,mstup),lmn1u,ncru,lproju,ndfac,lmax
      integer mgcsu,MAXCX,i,n,INT_AREP
      real*8 chg(msup), zcr(mcrup), ccr(mcrup), eps
      character*3 ityp(mstup), mtype(msup)

      CALL RECP_SETPARM
      mgcsu  = mgcsup

!     ================
!     input processing
!     ================
      CALL RECP_INP_MAIN(ccr,chg,eta,ica,idp,ityp,la,lb,lcr,
     &     lls,lmnp1,lmn1u,lproju,maords,mcons,
     &     mcrs,mgcs,mnl,ms,mtype,nblpr,nc,
     &     ncon,ncr,ncru,nct,nd,nf,nir,
     &     nkcrl,nkcru,nklsl,nklsu,nrcr,nso,nsopr,
     &     nt,ntl,ntu,x,y,z,zcr,zet,MAXCX)

      ALLOCATE(IPT_IL(nnbft))  !allocate space for il(*)

!     Compute symmetry coefficient products.
      CALL RECP_SOCFPD_MAIN(icxast,icxst,idp,iprst,la,
     &     lb,maords,mau,mcons,mgcs,nblpr,nc,nd,nf,nir,
     &     npair,nprir,nrcr,nsopr,nt,ntl,ntu,MAXCX)

!     Allocate variables
      eps = 1.0d-12
      CALL RECP_IPTA_CORTAB(ndfac,lmn1u,lproju,ncru,lmax)
      CALL RECP_CSO_TABLE0(eps,lmax,lmn1u,lproju,ndfac)

!     =========================
!     compute the ecp integrals
!     =========================
      INT_AREP = 0
      IF (nu.eq.4) INT_AREP = 1

      ALLOCATE (IPT_A(10))
      DO n = 4, nu
         nnam  = inam(n)
         CALL RECP_ONE_MAIN(IPT_A,ccr,chg,eta,ica,
     &        icb,icxast,icxst,ilxyz,iprst,lcr,lls,lmnp1,
     &        lproju,mcons,mcrs,nc,ncon,ncr,nf,nfct,ngw,nkcrl,
     &        nkcru,nklsl,nklsu,nopir,npair,nprir,nrcr,nt,ntl,
     &        ntu,x,y,z,zcr,zet,nso,INT_AREP)
      ENDDO

!     Deallocate variables
      DEALLOCATE (IPT_CX)
      CALL RECP_IPTD_CORTAB
      DEALLOCATE(IPT_IL)
      DEALLOCATE(IPT_LMNV)
      DEALLOCATE(IPT_A)
      IF (ALLOCATED(RECP_SET)) DEALLOCATE(RECP_SET)
      END


      SUBROUTINE RECP_SETPARM
      IMPLICIT NONE
#include "inc_mxvalue.h"
#include "inc_print.h" 
      integer
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      common /parmi/
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      integer
     & mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,inam,nnam,mdum
      common /ntgr/
     & mconsu,mru,mcru,msfu,msfru,ngcs,nu,
     & lxyzir(3),inam(5),nnam,mdum(32)
      integer i

      DATA RECP_FNUM / 71, 72, 73, 74, 75, 76, 77, 78, 79, 80 /
      DATA RECP_FNAM / 'RECP_XXXXX', 'RECP_XXXXX', 'RECP_XXXXX',
     &                 'RECP_H2___', 'RECP_XXXXX', 'RECP_XXXXX',
     &                 'RECP_XXXXX', 'RECP_XXXXX', 'RECP_XXXXX',
     &                 'RECP_XXXXX' /

!     parmi values
!     ------------
      mccu  = mccup
      mconu = mconup
      mcu   = mcup
      mrcru = mrcrup
      mstu  = mstup
      msu   = msup

!     ntgr values
!     -----------
      mconsu = mconsp
      mru    = mrup
      mcru   = mcrup
      msfu   = msfup
      msfru  = msfrup

      do i = 1, 3
        inam(i) = i
      enddo
      END SUBROUTINE RECP_SETPARM

