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

MODULE RECP_CSO
CONTAINS
SUBROUTINE CINTS(a,ccr,etai,etaj,g2,lcr,lambu, & 
           ltot1,mcrs,mproju,nc,ncr,nkcrl,nkcru,x,y,z,zcr,zet,tmp_ic0,tmp_ic1)
! lcru is the max l value + 1 for the potential.
! ncr contains the value of n for each term.
! zcr contains the value of alpha for each term.
! ccr contains the coefficient of each term.
  USE RECP_CSO_PS 
  USE RECP_CSO_SUB 
  USE RECP_IPT 
  USE RECP_OUTPUT 
  implicit real*8 (a-h,o-z)
#include "inc_print.h"
  logical esf, esfc, igueq1, jgueq1
  parameter (a1s2=0.5d0, a1=1.0d0, a4=4.0d0)
  common /parmr/ tol
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /one/ fnfct,rr,xij,xijm,yij,yijm,zij,zijm,ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic, &
               icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl, &
               jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
  common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
  common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa, &
                  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
  dimension a(*),ccr(*),etai(mrcru,*),etaj(mrcru,*),g2(*),lcr(*), &
            mcrs(*),nc(*),ncr(*),nkcrl(6,*),nkcru(6,*),x(mcu,*),  &
            y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
  integer i 
  integer tmp_ic0,tmp_ic1
  CALL OUTPUT_LOC('cint','E')

  CALL RECP_CSO_INIT1(iu,lit,mrcru,fctr1,etai,etaj,esfc,igueq1,jgueq1)
 
!     i primitive
!     -----------
      do 490 ig=1,igu

      CALL RECP_CSO_G1ZERO(nc1,tmp_ic0,igueq1)
      ai=zet(ig,icons)
      tai=ai+ai
 
!     j primitive
!     -----------
      if(esfc) jgu=ig
      do 390 jg=1,jgu

      CALL RECP_CSO_GOUTZERO(ij,jgueq1)
      aj=zet(jg,jcons)
      taj=aj+aj
      aa=ai+aj
      taa=aa+aa
      aaa=(ai-aj)/aa
      apr=ai*aj/aa
      aarr1=apr*rr
 
!     form core potential integrals
!     =============================
      xp=xij+aaa*xijm
      yp=yij+aaa*yijm
      zp=zij+aaa*zijm
      fctr2=fctr1
      if((esfc.and..not.igueq1).and.ig.eq.jg) fctr2=a1s2*fctr2
      do 350 ks=1,ns
        kcrs=mcrs(ks)
        if(kcrs.eq.0) goto 350
        lcru=lcr(kcrs)
        if(lcru.lt.0) goto 340
        do kc=1,nc(ks)
           xc=x(kc,ks)
           yc=y(kc,ks)
           zc=z(kc,ks)
           CALL RECP_CSO_SETCRDA(xka,yka,zka,xc,yc,zc,xi,yi,zi,ca,lit)
           CALL RECP_CSO_SETCRDB(xkb,ykb,zkb,xc,yc,zc,xj,yj,zj,cb,ljt)

!          type1 pseudo
!          ------------
           if(aarr1.le.tol) then
             xk=xp-xc
             yk=yp-yc
             zk=zp-zc
             call pseud1(ccr,IPT_CZ_CRDA,IPT_CZ_CRDB,IPT_CZ_GOUT, &
                  IPT_LMNV,ltot1,ncr,nkcrl,nkcru,zcr)
           endif

!          type2 pseudo
!          ------------
           aarr2=apr*(ca-cb)**2
           if (aarr2.le.tol.and.lcru.ne.0) then
              call pseud2(a,ccr,IPT_CZ_CRDA,IPT_CZ_CRDB, &
                   IPT_CZ_GOUT,lambu,ltot1,mproju,ncr,nkcrl,nkcru,zcr)
           endif
        enddo
  340   continue
  350 continue
 
!     j transformation
!     ----------------
      if(jgueq1) goto 400
 
      j1=0
      DO jrcr=1,jrcru
         DO i=1,ij
!           IF ((j1+i).GT.tmp_ic0) WRITE(RECP_OUT,*) 'G1 exceed(J)',(j1+1),'>',tmp_ic0
            IPT_CZ_G1(j1+i)=IPT_CZ_G1(j1+i)+IPT_CZ_GOUT(i)*etaj(jrcr,jg)
!           WRITE(RECP_OUT,*)'G1/GOUT',(j1+i),IPT_CZ_G1(j1+i), i,IPT_CZ_GOUT(i)
         ENDDO
         j1=j1+ij
      ENDDO
  390 continue
 
!     i transformation
!     ----------------
  400 if(igueq1) return
 
      if(esfc) i1=0
      ij2=0
      DO ircr=1,ircru
         j1=0
         IF (esfc) THEN
           CALL RECP_CSO_G2CAL1(ircr,ig,iu,i1,j1,ij2,tmp_ic0,tmp_ic1,mrcru,g2,etai)
           i1=i1+ij
         ELSE
           CALL RECP_CSO_G2CAL3(ircr,jrcru,ig,ij,ij2,j1,tmp_ic0,tmp_ic1,mrcru,g2,etai)
         ENDIF
      ENDDO
      
  490 continue
  CALL OUTPUT_LOC('RECP_CINT','X')
END SUBROUTINE CINTS


SUBROUTINE LSINTS(a,ccr,etai,etaj,g2,lls,lambu,ltot1,mcrs,mproju,nc,ncr,nklsl, &
           nklsu,x,y,z,zcr,zet,tmp_ic0,tmp_ic1)
! lcru is the max l value + 1 for the potential.
! ncr contains the value of n for each term.
! zcr contains the value of alpha for each term.
! ccr contains the coefficient of each term.
  USE RECP_CSO_PS
  USE RECP_CSO_SUB
  USE RECP_IPT
  USE RECP_OUTPUT
  implicit real*8 (a-h,o-z)
#include "inc_print.h"
  logical esf, esfc, igueq1, jgueq1
  parameter (a1s2=0.5d0, a1=1.0d0, a4=4.0d0)
  common /parmr/ tol
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /one/fnfct,rr,xij,xijm,yij,yijm,zij,zijm,ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,  &
              icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,  &
              jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
  common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
  common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa, &
                  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
  dimension a(*),ccr(*),etai(mrcru,*),etaj(mrcru,*),g2(*),lls(*),mcrs(*),nc(*),ncr(*), &
            nklsl(4,*),nklsu(4,*),x(mcu,*),y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
  integer i,tmp_ic0,tmp_ic1

  CALL OUTPUT_LOC('LSINTS','E') 
  CALL RECP_CSO_INIT1(iu,lit,mrcru,fctr1,etai,etaj,esfc,igueq1,jgueq1)
 
!     i primitive
!     -----------
      do 490 ig=1,igu

      CALL RECP_CSO_G1ZERO(nc1,tmp_ic0,igueq1)
      ai=zet(ig,icons)
      tai=ai+ai
 
!     j primitive
!     -----------
      if(esfc) jgu=ig
      do 390 jg=1,jgu

      CALL RECP_CSO_GOUTZERO(ij,jgueq1)
      aj=zet(jg,jcons)
      taj=aj+aj
      aa=ai+aj
      taa=aa+aa
      apr=ai*aj/aa
 
!     form spin-orbit potential integrals
!     ===================================
      fctr2=fctr1
      if((esfc.and..not.igueq1).and.ig.eq.jg) fctr2=a1s2*fctr2
      do 350 ks=1,ns
      kcrs=mcrs(ks)
      if (kcrs.eq.0) goto 350
      lcru=lls(kcrs)+1
      if (lcru.le.1) goto 340
      DO kc=1,nc(ks)
         xc=x(kc,ks)
         yc=y(kc,ks)
         zc=z(kc,ks)

         CALL RECP_CSO_SETCRDA(xka,yka,zka,xc,yc,zc,xi,yi,zi,ca,lit)
         CALL RECP_CSO_SETCRDB(xkb,ykb,zkb,xc,yc,zc,xj,yj,zj,cb,ljt)

         aarr2=apr*(ca-cb)**2
         IF (aarr2.LE.tol) THEN
            CALL PSEUD3(a,ccr,IPT_CZ_CRDA,IPT_CZ_CRDB,IPT_CZ_GOUT,lambu, &
                 ltot1,mproju,ncr,nklsl,nklsu,zcr)
         ENDIF
      ENDDO
  340 continue
  350 continue
 
!     j transformation
!     ----------------
      if (jgueq1) goto 400
 
      j1=0
      DO jrcr=1,jrcru
         DO i=1,ij
            IPT_CZ_G1(j1+i)=IPT_CZ_G1(j1+i)+IPT_CZ_GOUT(i)*etaj(jrcr,jg)
         ENDDO   
         j1=j1+ij
      ENDDO
  390 continue
  400 continue
 
!     i transformation
!     ----------------
      IF (igueq1) GOTO 500
 
      IF (esfc) i1=0
      ij2=0
      DO ircr=1,ircru
         j1=0
         IF (esfc) THEN
            CALL RECP_CSO_G2CAL2(ircr,ig,i1,j1,iu,tmp_ic0,tmp_ic1,ij2,mrcru,g2,etai)
            i1=i1+ij
         ELSE
            CALL RECP_CSO_G2CAL3(ircr,jrcru,ig,ij,ij2,j1,tmp_ic0,tmp_ic1,mrcru,g2,etai) 
         ENDIF
      ENDDO
  490 continue

  500 CONTINUE 
      CALL OUTPUT_LOC('LSINTS','X') 
END SUBROUTINE LSINTS
END MODULE RECP_CSO
