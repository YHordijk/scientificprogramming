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

MODULE RECP_CSO_PS
CONTAINS
SUBROUTINE pseud1(ccr,crda,crdb,gout,lmnv,ltot1,ncr,nkcrl,nkcru,zcr)
! compute type 1 core potential integrals
  USE RECP_CSO_ANG
  USE RECP_CSO_RAD
  USE RECP_FUNCTION1
  USE RECP_IPT
  USE RECP_OUTPUT
  implicit real*8 (a-h,o-z)
#include "inc_print.h"
  logical esf, esfc, igueq1, jgueq1
  parameter (a0=0.0d0)
  common /parmr/ tol
  common /one/fnfct,rr,xij,xijm,yij,yijm,zij,zijm,ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic, &
              icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl, &
              jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
  common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,aarr1,aarr2, &
                  xk,yk,zk,fctr2,kcrs,lcru
  dimension ccr(*),crda(lit,3),crdb(ljt,3),gout(*),lmnv(3,*),ncr(*),nkcrl(6,*),nkcru(6,*),zcr(*)
  integer i,j

! # should include kcr loop in recur1 and use recurrence relations on qsum
  rp2=xk*xk+yk*yk+zk*zk

  if(rp2.eq.0.0d0) then
     rp=0.0d0
     arp2=0.0d0
     alpt=0.0d0
     rk=0.0d0
     lamu=1
  else
     rp=sqrt(rp2)
     xk=xk/rp
     yk=yk/rp
     zk=zk/rp
     arp2=aa*rp2
     alpt=aa*arp2
     rk=taa*rp
     lamu=ltot1
  endif
 
! compute radial integrals and sum over potential terms
! =====================================================

  do i = 1, ltot1
!    do j = 1, lamu
     do j = 1, ltot1
        IPT_P1_QSUM(I,J) = 0.0d0
     enddo
  enddo

  kcrl=nkcrl(1,kcrs)
  kcru=nkcru(1,kcrs)
!     do 40 kcr=kcrl,kcru
!       alpha=aa+zcr(kcr)
!       # exponential factor from q functions included in dum
!       dum=aarr1+zcr(kcr)*arp2/alpha
!       if(dum.gt.tol) go to 40
!       prd=fctr2*ccr(kcr)*exp(-dum)
!       t=alpt/alpha
!       call recur1(alpha,a(ipt(11)),ncr(kcr),ltot1,q,t,rk)
!       do 30 lam=1,lamu
!         nhi=ltot1-mod(ltot1-lam,2)
!         do 20 n=lam,nhi,2
!           IPT_P1_QSUM(n,lam)=IPT_P1_QSUM(n,lam)+prd*IPT_P1_Q(n,lam)
!  20     continue
!  30   continue
!  40 continue

! type1 integral
  call RECP_CSO_RAD1(aa,aarr1,alpt,arp2,ccr,fctr2,kcrl,kcru,lamu,ltot1,ncr,rk,tol,zcr)

  ijt=0
  do it=itl,itu
!    print *,'pseud1-it',it 
     na1=lmnv(1,it)+1
     la1=lmnv(2,it)+1
     ma1=lmnv(3,it)+1
     do jt=jtl,jtu
!       print *,'pseud1-jt',jt 
        nb1=lmnv(1,jt)+1
        lb1=lmnv(2,jt)+1
        mb1=lmnv(3,jt)+1
 
!       compute angular integrals
!       -------------------------
        CALL RECP_FUNCTION_FACAB(IPT_SH_BINOM,na1,nb1,crda(1,1),crdb(1,1),IPT_P1_XAB)
        CALL RECP_FUNCTION_FACAB(IPT_SH_BINOM,la1,lb1,crda(1,2),crdb(1,2),IPT_P1_YAB)
        CALL RECP_FUNCTION_FACAB(IPT_SH_BINOM,ma1,mb1,crda(1,3),crdb(1,3),IPT_P1_ZAB)
        CALL RECP_CSO_ANG1(na1+nb1-1,la1+lb1-1,ma1+mb1-1,lamu,ltot1,xk,yk,zk)

!       combine angular and radial integrals
!       ------------------------------------
        s=a0
        ijt=ijt+1
        DO lam=1,lamu
           nhi=ltot1-mod(ltot1-lam,2)
           DO n=lam,nhi,2
              s=s+IPT_P1_ANG(n,lam)*IPT_P1_QSUM(n,lam)
           ENDDO
        ENDDO
        gout(ijt)=gout(ijt)+s
     ENDDO
  ENDDO
END SUBROUTINE pseud1


SUBROUTINE pseud2(a,ccr,crda,crdb,gout,lambu,ltot1,mproju,ncr,nkcrl,nkcru,zcr)
! compute type 2 core potential integrals
  USE RECP_CSO_ANG
  USE RECP_CSO_RAD
  USE RECP_IPT
  USE RECP_OUTPUT
  implicit real*8 (a-h,o-z)
#include "inc_print.h"
  logical esf, esfc, igueq1, jgueq1
  parameter (a0=0.0d0)
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /one/   fnfct,rr,xij,xijm,yij,yijm,zij,zijm,ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic, &
                 icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl, &
                 jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
  common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa, &
                  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
  dimension a(*),ccr(*),crda(lit,3),crdb(ljt,3),gout(*),ncr(*),nkcrl(6,*),nkcru(6,*),zcr(*)
  CALL OUTPUT_LOC('RECP_PSEUD2','E')
  if(ca.eq.a0) then
    rka=a0
    lmau=1
  else
    xka=-xka/ca
    yka=-yka/ca
    zka=-zka/ca
    rka=tai*ca
    lmau=lcru+(lit-1)
  endif
  if(cb.eq.a0) then
    rkb=a0
    lmbu=1
  else
    xkb=-xkb/cb
    ykb=-ykb/cb
    zkb=-zkb/cb
    rkb=taj*cb
    lmbu=lcru+(ljt-1)
  endif
  if((ca.eq.a0).and.(cb.eq.a0)) then
    lhi=min(lcru,lit,ljt)
    llo=mod((lit-1),2)+1
    if(llo.ne.mod((ljt-1),2)+1.or.llo.gt.lhi) return
    inc=2
  elseif(ca.eq.a0) then
    lhi=min(lcru,lit)
    llo=mod((lit-1),2)+1
    if(llo.gt.lhi) return
    inc=2
  elseif(cb.eq.a0) then
    lhi=min(lcru,ljt)
    llo=mod((ljt-1),2)+1
    if(llo.gt.lhi) return
    inc=2
  else
    lhi=lcru
    llo=1
    inc=1
  endif
      do 88 l=llo,lhi,inc
      mhi=l+l-1
      lmalo=max(l-(lit-1),1)
      lmahi=min(lmau,l+(lit-1))
      lmblo=max(l-(ljt-1),1)
      lmbhi=min(lmbu,l+(ljt-1))
 
!     compute radial integrals
!     ------------------------
      kcrl=nkcrl(l+1,kcrs)
      kcru=nkcru(l+1,kcrs)
      CALL RECP_CSO_RAD2(ccr,kcrl,kcru,l,lambu,lmahi,lmalo,lmbhi,lmblo,ltot1,ncr,rka,rkb,zcr )
 
!     compute angular integrals and combine with radial integrals
!     -----------------------------------------------------------
      ijt=0
      do 84 it=itl,itu
      CALL RECP_CSO_ANG2(IPT_P23_ANGA,crda,it,l,lit,lmalo,lmahi,IPT_LMNV,mproju,xka,yka,zka)
      do 80 jt=jtl,jtu
      ijt=ijt+1
      s=a0
      CALL RECP_CSO_ANG2(IPT_P23_ANGB,crdb,jt,l,ljt,lmblo,lmbhi,IPT_LMNV,mproju,xkb,ykb,zkb)
      do 76 lama=lmalo,lmahi
        ldifa1=abs(l-lama)+1
        nlmau=lit-mod(lit-ldifa1,2)
        do 72 lamb=lmblo,lmbhi
          ldifb=abs(l-lamb)
          nlmbu=(ljt-1)-mod((ljt-1)-ldifb,2)
          nlo=ldifa1+ldifb
          nhi=nlmau+nlmbu
          do 68 n=nlo,nhi,2
            nlmalo=max(ldifa1,n-nlmbu)
            nlmahi=min(nlmau,n-ldifb)
            angp=a0
            do 60 m=1,mhi
              do 56 nlma=nlmalo,nlmahi,2
                angp=angp+IPT_P23_ANGA(nlma,m,lama)*IPT_P23_ANGB((n+1)-nlma,m,lamb)
   56         continue
   60       continue
            s=s+angp*IPT_P23_QSUM(n,lamb,lama)
   68     continue
   72   continue
   76 continue
      gout(ijt)=gout(ijt)+s
   80 continue
   84 continue
   88 continue
      CALL OUTPUT_LOC('RECP_PSEUD2','X')
      return
END SUBROUTINE PSEUD2


SUBROUTINE PSEUD3(a,ccr,crda,crdb,gout,lambu,ltot1,mproju,ncr,nklsl,nklsu,zcr )
! computes spin-orbit potential integrals
  USE RECP_CSO_ANG
  USE RECP_CSO_RAD
  USE RECP_IPT
  USE RECP_OUTPUT
  implicit real*8 (a-h,o-z)
#include "inc_print.h"
  logical esf, esfc, igueq1, jgueq1
  parameter (a0=0.0d0)
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /one/   fnfct,rr,xij,xijm,yij,yijm,zij,zijm,ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic, &
                 icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl, &
                 jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
  common /callin/xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
  dimension a(*),ccr(*),crda(lit,3),crdb(ljt,3),gout(*),ncr(*),nklsl(4,*),nklsu(4,*),zcr(*)

  CALL OUTPUT_LOC('RECP_PSEUD3','E')
  if(ca.eq.a0) then
    rka=a0
    lmau=1
  else
    xka=-xka/ca
    yka=-yka/ca
    zka=-zka/ca
    rka=tai*ca
    lmau=lcru+(lit-1)
  endif
  if(cb.eq.a0) then
    rkb=a0
    lmbu=1
  else
    xkb=-xkb/cb
    ykb=-ykb/cb
    zkb=-zkb/cb
    rkb=taj*cb
    lmbu=lcru+(ljt-1)
  endif
  if((ca.eq.a0).and.(cb.eq.a0)) then
    lhi=min(lcru,lit,ljt)
    llo=mod(lit,2)+2
    if(llo.ne.(mod(ljt,2)+2).or.llo.gt.lhi) return
    inc=2
  elseif(ca.eq.a0) then
    lhi=min(lcru,lit)
    llo=mod(lit,2)+2
    if(llo.gt.lhi) return
    inc=2
  elseif(cb.eq.a0) then
    lhi=min(lcru,ljt)
    llo=mod(ljt,2)+2
    if(llo.gt.lhi) return
    inc=2
  else
    lhi=lcru
    llo=2
    inc=1
  endif
      do 88 l=llo,lhi,inc
      mhi=l+l-3
      lmalo=max(l-(lit-1),1)
      lmahi=min(lmau,l+(lit-1))
      lmblo=max(l-(ljt-1),1)
      lmbhi=min(lmbu,l+(ljt-1))
 
!     compute radial integrals
!     ------------------------
      kcrl=nklsl(l-1,kcrs)
      kcru=nklsu(l-1,kcrs)
      CALL RECP_CSO_RAD2( ccr,kcrl,kcru,l,lambu,lmahi,lmalo,lmbhi,lmblo,ltot1,ncr,rka,rkb,zcr )
 
!     compute angular integrals and combine with radial integrals
!     -----------------------------------------------------------
      ijt=0
      do 84 it=itl,itu
      CALL RECP_CSO_ANG2(IPT_P23_ANGA,crda,it,l,lit,lmalo,lmahi,IPT_LMNV,mproju,xka,yka,zka)
      do 80 jt=jtl,jtu
      CALL RECP_CSO_ANG2(IPT_P23_ANGB,crdb,jt,l,ljt,lmblo,lmbhi,IPT_LMNV,mproju,xkb,ykb,zkb)
      do 76 lama=lmalo,lmahi
      ldifa1=abs(l-lama)+1
      nlmau=lit-mod(lit-ldifa1,2)
      do 72 lamb=lmblo,lmbhi
      ldifb=abs(l-lamb)
      nlmbu=(ljt-1)-mod((ljt-1)-ldifb,2)
      nlo=ldifa1+ldifb
      nhi=nlmau+nlmbu
      do 68 n=nlo,nhi,2
      nlmalo=max(ldifa1,n-nlmbu)
      nlmahi=min(nlmau,n-ldifb)
      do 64 kt=1,3
        s=a0
        do 60 m=1,mhi
          angp=a0
          do 56 nlma=nlmalo,nlmahi,2
            angp=angp + IPT_P23_ANGA(nlma,      IPT_SH_MR(kt,m),lama)  &
                      * IPT_P23_ANGB((n+1)-nlma,IPT_SH_MC(kt,m),lamb)  &
                      - IPT_P23_ANGA(nlma,      IPT_SH_MC(kt,m),lama)  &
                      * IPT_P23_ANGB((n+1)-nlma,IPT_SH_MR(kt,m),lamb)
   56     continue
          s=s+angp*IPT_SH_FLMTX(kt,(l-2)**2+m)
   60   continue
        gout(ijt+kt)=gout(ijt+kt)+s*IPT_P23_QSUM(n,lamb,lama)
   64 continue
   68 continue
   72 continue
   76 continue
      ijt=ijt+3
   80 continue
   84 continue
   88 continue
      CALL OUTPUT_LOC('RECP_PSEUD3','X')
END SUBROUTINE PSEUD3


END MODULE RECP_CSO_PS
