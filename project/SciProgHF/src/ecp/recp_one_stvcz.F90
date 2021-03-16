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

MODULE RECP_ONE_STVCZ
CONTAINS
SUBROUTINE RECP_ONE_STVCZ0(a,ccr,chg,eta,ilxyz,lcr,lls,lmnv,lproju,mcrs,nc,ncr, &
           nkcrl,nkcru,nklsl,nklsu,nopir,nprir,nt,x,y,z,zcr,zet)
  USE RECP_ONE_SUB
  USE RECP_IPT
  USE RECP_CSO
  USE RECP_OUTPUT
  implicit real*8 (a-h,o-z)
#include "inc_print.h"
  logical chsign,ec,esf,sfb,esfc,esfb,esfbc,igueq1,jgueq1
  parameter (a1s2=0.5d0)
  common /parmr/ tol
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  integer        mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,   inam,   nnam,mdum
  common /ntgr/  mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir(3),inam(5),nnam,mdum(32)
  common /one/ fnfct,rr,xij,xijm,yij,yijm,zij,zijm,ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic, &
               icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl, &
               jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
  common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
  dimension a(*),ccr(*),chg(*),eta(mrcru,mconu,*),ilxyz(3,*),lcr(*), &
            lls(*),lmnv(3,*),mcrs(*),nc(*),ncr(*),nkcrl(6,*),nkcru(6,*), &
            nklsl(4,*),nklsu(4,*),nopir(*),nprir(2,mstu,*),nt(*),x(mcu,*), &
            y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
  integer ic0,ic1,ltot1,mproju,lproju,lamau,lambu,i,IPQ_TEMP 
  integer tmp_ic0,tmp_ic1

  CALL OUTPUT_LOC('RECP_STVCZ','E')

  lijt = lit*ljt
  ec   = ic.eq.jc
  esfc = esf.and.ec
  icx2 = icxi2+((ic-1)*nc(js)+jc-1)*ntij2
  if(.not.esfc.and.ircru.ne.1.and.esf) then
    ibld=ibl2-ibl1
    icx3=icx2-(ic-jc)*(nc(is)-1)*ntij2
  endif
  if(esf) then
    IPQ_TEMP = ((ic-1)*nt(isf)+1)*(((ic-1)*nt(isf)+1)-1)/2
    icx1=icxi1+IPQ_TEMP*nblt1+(jc-1)*ntij1
  endif
 
  xi=x(ic,is)  ! ishell
  yi=y(ic,is)
  zi=z(ic,is)
 
  xj=x(jc,js)  ! jshell
  yj=y(jc,js)
  zj=z(jc,js)

  xij=(xi+xj)/2.D0
  yij=(yi+yj)/2.D0
  zij=(zi+zj)/2.D0

  xijm=(xi-xj)/2.D0
  yijm=(yi-yj)/2.D0
  zijm=(zi-zj)/2.D0

  rr = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj)
!
! prepare items for pairs of (i,j) functions
!
  ij=0
  do it=itl,itu
!   print *,'stvcz-it',it
    nx=ljt*lmnv(1,it)+1
    ny=ljt*lmnv(2,it)+1
    nz=ljt*lmnv(3,it)+1
    do jt=jtl,jtu
      ij=ij+1
    enddo
  enddo
  ij=nop*ij


  igueq1=igu.eq.1
  jgueq1=jgu.eq.1

! Allocate nc1/2 , ic1/2, tmp_ic0/1
! ---------------------------------
  CALL RECP_1ESUB_NC12(nc1,nc2,ircru,jrcru,ij,esfc)   !set nc1/nc2
  CALL RECP_1ESUB_IC12(ic0,ic1,nc1,nc2,igueq1,jgueq1) !set ic0/ic1

  tmp_ic0 = ic0
  tmp_ic1 = ic1

  IF (RECP_DBG.GE.10) THEN
     WRITE(RECP_OUT,'(3X,A)') '* Size of variables, G2/G1/GOUT'
     WRITE(RECP_OUT,'(3X,5(3X,A3))')'nc1','nc2','ic0','ic1','ij'
     WRITE(RECP_OUT,'(3X,5I6)') nc1,nc2,ic0,ic1,ij
  ENDIF

! Allocate variables
  CALL RECP_IPTA_STVCZ(ic0,ic1,ij,lit,ljt,ltot1,mproju,lproju,lamau,lambu)

! ------------------
! Calculate integral 
! ------------------
  if     (nnam.eq.1) then
    IF (RECP_DBG.GE.1) WRITE(6,'(X,A)') '* Overlap(S) integral calculation is skipped'
  elseif (nnam.eq.2) then
    IF (RECP_DBG.GE.1) WRITE(6,'(X,A)') '* Kinetic(T) integral calculation is skipped'
  elseif (nnam.eq.3) then
    IF (RECP_DBG.GE.1) WRITE(6,'(X,A)') '* Coulomb(V) integral calculation is skipped'
  elseif (nnam.eq.4) then
        call cints( a,ccr,eta(1,1,icons),eta(1,1,jcons),IPT_CZ_G2, &
                    lcr,lambu,ltot1,mcrs,mproju,nc,ncr,nkcrl, &
                    nkcru,x,y,z,zcr,zet,tmp_ic0,tmp_ic1 )
  elseif (nnam.eq.5) then
        call lsints( a,ccr,eta(1,1,icons),eta(1,1,jcons),IPT_CZ_G2, & 
                     lls,lambu,ltot1,mcrs,mproju,nc,ncr,nklsl, &
                     nklsu,x,y,z,zcr,zet,tmp_ic0,tmp_ic1 )
  endif

! --------------------------------------------

  ndx=0
  ndxa=0
  DO ircr=1,ircru     !#cont.coeff(i)
     IF (esfc) jrcru=ircr
     DO jrcr=1,jrcru  !#cont.coeff(j)
        esfb  =esf.and.ircr.eq.jrcr
        esfbc =esfb.and.ec
        chsign=nnam.ge.5.and.esf.and.ircr.lt.jrcr

!       Set index for IPT_ONE_H2
        CALL RECP_1ESUB_H2INDEX(icx,icx1,icx2,icx3,ircr,jrcr,ndx, &
         ndxa,ndxi,ibl2,ibld,ngti,ngtj,nop,isf,jsf,nt,esfb,esfc,esf)

!       Write to integral into IPT_ONE_H2
        DO it=itl,itu
           ndxj=ndxi
           IF (esfbc) jtu=it
           ndxi=ndxi+ngti
           DO jt=jtl,jtu
              index1=ndxj
              ndxj=ndxj+ngtj
              ndxb=ndxa
              IF (RECP_DBG.GE.1) WRITE(6,'(A,4(X,A,I5))') &
                 '* H2:','ircr',ircr,'jrcr',jrcr,'it',it,'jt',jt
              CALL RECP_1ESUB_H2WRITE(ic1,ic0,INDEX1,ilxyz,nst,ijsf, &
                   ij,ndxb,icx,nopir,mstu,nprir,chsign,esfb,fnfct)
           ENDDO
        ENDDO

        ndxa=ndxb
        ndx=ndx+ij
     ENDDO
  ENDDO

! --------------------------------------------

! Deallocate variables
  CALL RECP_IPTD_STVCZ
  CALL OUTPUT_LOC('RECP_STVCZ','X')

END SUBROUTINE RECP_ONE_STVCZ0
END MODULE RECP_ONE_STVCZ
