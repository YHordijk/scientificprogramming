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

MODULE RECP_CSO_RAD
CONTAINS

SUBROUTINE RECP_CSO_RAD1(aa,aarr1,alpt,arp2,ccr,fctr2,  &
           kcrl,kcru,lamu,ltot1,ncr,rk,tol,zcr)
! compute type 1 radial integrals
  USE RECP_IPT 
  USE RECP_FUNCTION2
  USE RECP_OUTPUT 
  implicit real(8) (a-h,o-z)
#include "inc_print.h"
  dimension ccr(*), ncr(*), zcr(*)
! global variables
  REAL(8) :: aa
! local variables
  INTEGER :: AREP_INDEX, AREP_N
  REAL(8) :: alpha

  CALL OUTPUT_LOC('RECP_CSO_RAD1','E')
! ---------------------  
! AREP block(kcrl~kcru) 
! ---------------------  
  DO AREP_INDEX = kcrl,kcru
!    AREP : AREP_C * r^(AREP_N-2) * exp(-AREP_Z*r*r)
     AREP_N = ncr(AREP_INDEX)     
!    alpha_j= alpha_A + alpha_B + AREP_Z ?? 
     alpha  = aa+zcr(AREP_INDEX)  

!    # exponential factor from q functions included in dum
     dum=aarr1+zcr(AREP_INDEX)*arp2/alpha
     IF (dum.LE.tol) THEN
 
        prd=fctr2*ccr(AREP_INDEX)*exp(-dum)
 
        IF (rk.EQ.0.0d0) THEN
           t=0.0d0
           DO n=1,ltot1-mod(ltot1-1,2),2
              IPT_P1_QSUM(n,1) = IPT_P1_QSUM(n,1)  &
                               + prd * qcomp(alpha,IPT_SH_DFAC,AREP_N+n-1,0,t,rk)
           ENDDO
        ELSE
           t=alpt/alpha
           DO lam=1,lamu
              IPT_P1_QSUM(lam,lam) = IPT_P1_QSUM(lam,lam)  &
                                   + prd * qcomp(alpha,IPT_SH_DFAC,AREP_N+lam-1,lam-1,t,rk)
           ENDDO
        ENDIF
     ENDIF
  ENDDO
 
! ---------------------  
! QSUM
! ---------------------  
  if(rk.ne.0.0d0) then
    f2lam3=(lamu+lamu-3)
    do lam=lamu-2,1,-1
      do n=lam+2,lamu-mod(lamu-lam,2),2
!              sum(d*Q)   =       sum(d*Q)     + ((2l+3)/k)* sum(d*Q)
        IPT_P1_QSUM(n,lam)=IPT_P1_QSUM(n,lam+2)+(f2lam3/rk)*IPT_P1_QSUM(n-1,lam+1)
      enddo
      f2lam3=f2lam3-2.0d0
    enddo
  endif
  CALL OUTPUT_LOC('RECP_CSO_RAD1','X')
END SUBROUTINE RECP_CSO_RAD1


SUBROUTINE RECP_CSO_RAD2(ccr,kcrl,kcru,l,lambu,lmahi,lmalo,lmbhi,lmblo,ltot1,ncr,rka,rkb,zcr )
! compute type 2 radial integrals.
  USE RECP_IPT
  USE RECP_FUNCTION2
  implicit real*8 (a-h,o-z)
  logical esf, esfc, igueq1, jgueq1
  parameter (a0=0.0d0, eps1=1.0d-15, a1=1.0d0, a2=2.0d0,a4=4.0d0, a50=50.0d0)
  common /parmr/ tol
  common /one/fnfct,rr,xij,xijm,yij,yijm,zij,zijm,ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic, &
              icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl, &
              jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
  common /callin/xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
  dimension ccr(*),ncr(*),  zcr(*)

  do i = 1, ltot1
     do j = 1, lambu
        do k = 1, lmahi
           IPT_P23_QSUM(i,j,k) = 0.0d0
        enddo
     enddo
  enddo

!     # sum over potential terms
      do 68 kcr=kcrl,kcru
      npi=ncr(kcr)
      alpha=aa+zcr(kcr)
      rc=(rka+rkb)/(alpha+alpha)
      arc2=alpha*rc*rc
      dum=aarr2+zcr(kcr)*arc2/aa
      if(dum.gt.tol) go to 68
      prd=fctr2*ccr(kcr)*exp(-dum)
 
      if(rka.eq.a0.and.rkb.eq.a0) then
!       # rka=0 and rkb=0
        rk=a0
        t=a0
        IPT_P23_QSUM(ltot1,1,1) = IPT_P23_QSUM(ltot1,1,1)  &
                                + prd*qcomp(alpha,IPT_SH_DFAC,npi+ltot1-1,0,t,rk)
      elseif(rka.eq.a0) then
!       # rka=0 and rkb>0
        rk=rkb
        t=arc2
        do lamb=l,lmbhi
          IPT_P23_QSUM(lamb-l+lit,lamb,1) = IPT_P23_QSUM(lamb-l+lit,lamb,1)  & 
              + prd*qcomp(alpha,IPT_SH_DFAC,npi+lamb-l+lit-1,lamb-1,t,rk)
        enddo
      elseif(rkb.eq.a0) then
!       # rka>0 and rkb=0
        rk=rka
        t=arc2
        do lama=l,lmahi
          IPT_P23_QSUM(lama-l+ljt,1,lama) = IPT_P23_QSUM(lama-l+ljt,1,lama)  &
             + prd*qcomp(alpha,IPT_SH_DFAC,npi+lama-l+ljt-1,lama-1,t,rk)
        enddo
      elseif(npi.eq.2) then
!       # rka>0 and rkb>0; use bessel function formula.
!       # to be applicable for a set of integrals, must have nu.le.l and
!       # nu.eq.(integer), where nu=l+1-npi/2, so it is used here only
!       # for the npi=2 case.  It can't be used at all for npi=(odd) and
!       # only for partial sets for npi=0
        nu=l
        call qbess( alpha,l,lambu,lmahi,lmbhi,ltot1,nu,prd,rka,rkb,ljt,lit )
      elseif(arc2.ge.a50) then
!       # rka>0 and rkb>0; use pts and wts method
!       # estimate radial integrals and compare to threshold
        qlim=abs(prd)/(max(a1,(rc+rc)*rka)*max(a1,(rc+rc)*rkb))*  &
          sqrt(a4*(tai+tai)**lit*(taj+taj)**ljt*sqrt(tai*taj)/alpha)
        if(rc.lt.ca) then
          nlim=npi
          qlim=qlim*ca**(lit-1)
        else
          nlim=npi+(lit-1)
        endif
        if(rc.lt.cb) then
          qlim=qlim*cb**(ljt-1)
        else
          nlim=nlim+(ljt-1)
        endif
        if(qlim*rc**nlim.ge.eps1) then
          call ptwt(arc2,npi,l,lambu,ltot1,lmahi,lmbhi,alpha,rc,rka,rkb,prd)
        endif
      else
!       # rka>0 and rkb>0; use partially asymptotic method
        call qpasy(alpha,npi,l,lambu,lmahi,lmbhi,ltot1,rka,rkb,fctr2*ccr(kcr),dum+arc2)
      endif
   68 continue
 
      if(rka.eq.a0.and.rkb.ne.a0) then
!       # rka=0 and rkb>0
        f2lmb3=(2*lmbhi-3)
        do lamb=lmbhi-2,lmblo,-1
          nlo=abs(lamb-l+1)+lit+1
          nhi=ljt-mod((ljt-1)-abs(lamb-l),2)+lit-1
          do n=nlo,nhi,2
            IPT_P23_QSUM(n,lamb,1) = IPT_P23_QSUM(n,lamb+2,1)  &
                                   + (f2lmb3/rkb)*IPT_P23_QSUM(n-1,lamb+1,1)
          enddo
          f2lmb3=f2lmb3-a2
        enddo
      elseif(rka.ne.a0.and.rkb.eq.a0) then
!       # rka>0 and rkb=0
        f2lma3=(2*lmahi-3)
        do lama=lmahi-2,lmalo,-1
          nlo=abs(lama-l+1)+ljt+1
          nhi=lit-mod((lit-1)-abs(lama-l),2)+ljt-1
          do n=nlo,nhi,2
            IPT_P23_QSUM(n,1,lama) = IPT_P23_QSUM(n,1,lama+2)  &
                                   + (f2lma3/rka)*IPT_P23_QSUM(n-1,1,lama+1)
          enddo
          f2lma3=f2lma3-a2
        enddo
      elseif(rka.ne.a0.and.rkb.ne.a0) then
!       # rka>0 and rkb>0
        f2lma3=(lmahi+lmahi+1)
        do 96 lama=lmahi,lmalo,-1
          ldifa1=abs(l-lama)+1
          f2lmb3=(2*lmbhi+1)
          do 92 lamb=lmbhi,lmblo,-1
            ldifb=abs(l-lamb)
            nlo=ldifa1+ldifb
            nhi=(ltot1-mod(lit-ldifa1,2))-mod((ljt-1)-ldifb,2)
            do 88 n=nlo,nhi,2
              if(n-(lama+lamb).eq.(1-l-l)) go to 88
              if(lama.gt.(lmahi-2).or.n.le.(abs(l-lama-2)+ldifb)) then
!               # lamb recursion
                IPT_P23_QSUM(n,lamb,lama) = IPT_P23_QSUM(n,lamb+2,lama)  &
                                          + (f2lmb3/rkb)*IPT_P23_QSUM(n-1,lamb+1,lama)
              else
!               # lama recursion
                IPT_P23_QSUM(n,lamb,lama) = IPT_P23_QSUM(n,lamb,lama+2)  &
                                          + (f2lma3/rka)*IPT_P23_QSUM(n-1,lamb,lama+1)
              endif
   88       continue
            f2lmb3=f2lmb3-a2
   92     continue
          f2lma3=f2lma3-a2
   96   continue
      endif
      return
END SUBROUTINE RECP_CSO_RAD2

!-------------------------------------------------

SUBROUTINE QBESS( alpha,l,lambu,lmahi,lmbhi,ltot1,nu,prd,rka,rkb,ljt,lit )
! compute type 2 radial integrals, scaled by exp(-arc2)/sqrt(pi),
! using the bessel function formula for
! lama=max(l,nu) to lmahi, lamb=max(l,nu) to lmbhi, n=lama+lamb-l-l
  USE RECP_IPT
  implicit real*8 (a-h,o-z)
  integer  ljt,lit
  parameter (a0=0.0d0, a1=1.0d0, a2=2.0d0, a4=4.0d0)

! nu=l+1-npi/2
! # bessel function formula applies to all npi=2 cases, no npi=1
! # cases, and some npi=0 cases.
  num1=nu-1
  lmlo=max(l,nu)
  lmaphi=lmahi-num1
  lmbphi=lmbhi-num1

  fcta = rka/(alpha+alpha)
  fctb = rkb/(alpha+alpha)

! Set BPREF
  IPT_P23_BPREF(1)=fctb**num1
  DO lambp=2,lmbphi
     IPT_P23_BPREF(lambp)=fctb*IPT_P23_BPREF(lambp-1)
  ENDDO

! Set APWR 
  fct = rka*fcta
  IPT_P23_APWR(1)=1.0d0
  IF (ljt.NE.1) THEN
     IPT_P23_APWR(2)=fct
  ENDIF
  DO lambp=3,lmbphi
     IPT_P23_APWR(lambp)=fct*IPT_P23_APWR(lambp-1)
  ENDDO

! Set BPWR 
  fct = rkb*fctb
  IPT_P23_BPWR(1)=1.0d0
  IF (lit.NE.1) THEN
     IPT_P23_BPWR(2)=fct
  ENDIF
  DO lamap=3,lmaphi
     IPT_P23_BPWR(lamap)=fct*IPT_P23_BPWR(lamap-1)
  ENDDO

  lmihi=lmaphi+lmbphi+(nu-2)
  call ssibfn(lmihi-1,rka*fctb,IPT_P23_SSI)
  DO lami=1,lmihi
    IPT_P23_SSI(lami)=IPT_P23_SSI(lami)/IPT_SH_DFAC(lami+lami-1)
  ENDDO
  lmplo=lmlo-num1
  fctra=(alpha+alpha)**(nu-2)*fcta**(lmlo-1)*prd/sqrt(a4*alpha)*    &
   ((IPT_SH_DFAC(2*(2*lmplo+num1)-1)/IPT_SH_DFAC(2*(lmplo+num1)+1)) &
    *IPT_SH_DFAC(2*num1+1)) / IPT_SH_DFAC(2*(lmplo+num1)+1) 
  fctran=(2*(2*lmplo+num1)-1)
  fctrad=(2*(lmplo+num1)+1)

  DO lamap=lmplo,lmaphi
     fctru=a1
     fctrun=(nu+num1)
     fctrud=(2*(lamap+num1)+1)
     DO iu=1,lmbphi
        IPT_P23_ATERM1(iu)=fctru*IPT_P23_APWR(iu)
        fctru=fctru*fctrun/fctrud
        fctrun=fctrun+a2
        fctrud=fctrud+a2
     ENDDO
     DO it=1,lamap
        IPT_P23_BTERM1(it)=IPT_SH_BINOM((lamap*(lamap-1))/2+it) * IPT_P23_BPWR(it)
     ENDDO
     fctrb=fctra
     fctrbn=(2*(lamap+lmplo+num1)-1)
     fctrbd=(2*(lmplo+num1)+1)
     DO lambp=lmplo,lmbphi
        n=((2*(nu-l)-1)+lamap)+lambp
        DO iu=1,lambp
           IPT_P23_ATERM2(iu)=IPT_SH_BINOM((lambp*(lambp-1))/2+iu) * IPT_P23_ATERM1(iu)
        ENDDO
        sum=a0
        fctrt=a1
        fctrtn=(nu+num1)
        fctrtd=(2*(lambp+num1)+1)
        DO it=1,lamap
           DO iu=1,lambp
              sum=sum+IPT_P23_ATERM2(iu)*(fctrt*IPT_P23_BTERM1(it))*IPT_P23_SSI((it+(nu-2))+iu)
           ENDDO
           fctrt=fctrt*fctrtn/fctrtd
           fctrtn=fctrtn+a2
           fctrtd=fctrtd+a2
        ENDDO
        IPT_P23_QSUM(n,lambp+num1,lamap+num1) = IPT_P23_QSUM(n,lambp+num1,lamap+num1) &
                                              + fctrb * IPT_P23_BPREF(lambp) * sum
        fctrb=fctrb*fctrbn/fctrbd
        fctrbn=fctrbn+a2
        fctrbd=fctrbd+a2
     ENDDO
     fctra=fcta*fctra*fctran/fctrad
     fctran=fctran+a2
     fctrad=fctrad+a2
  ENDDO
END SUBROUTINE QBESS


subroutine ptwt(arc2,npi,l,lambu,ltot1,lmahi,lmbhi,alpha,rc,rka,rkb,prd)
! compute type 2 radial integrals, scaled by exp(-arc2)/sqrt(pi),
! using the points and weights method,
! for lama=l to lmahi, lamb=l to lmbhi, n=lama+lamb-l-l
  USE RECP_IPT
  implicit real*8 (a-h,o-z)
  parameter (a500=500.0d0, a50000=50000.0d0)

  DO I = 1, lambu
     DO J = 1, lmahi
        IPT_P23_Q2(I,J) = 0.0d0
     ENDDO
  ENDDO

  if(arc2.gt.a50000) then
    npt=5
    idif=0
  else if(arc2.gt.a500) then
    npt=10
    idif=5
  else
    npt=20
    idif=15
  endif
  sqalp=sqrt(alpha)
  prd=prd/sqalp

  do i=1,npt
     pt=rc+IPT_SH_HPT(i+idif)/sqalp
     call ssibfn(lmahi-1,rka*pt,IPT_P23_ABESS)
     call ssibfn(lmbhi-1,rkb*pt,IPT_P23_BBESS)
     if((npi+l+l-2).eq.0) then
        IPT_P23_PTPOW(1)=prd
     else
        IPT_P23_PTPOW(1)=prd*pt**(npi+l+l-2)
     endif
     do n=2,ltot1
        IPT_P23_PTPOW(n)=(pt*pt)*IPT_P23_PTPOW(n-1)
     enddo
     do lama=l,lmahi
        do lamb=l,lmbhi
           n=((1-l-l)+lama)+lamb
           IPT_P23_Q2(lamb,lama) =IPT_P23_Q2(lamb,lama)  &
              + (IPT_SH_HWT(i+idif)*IPT_P23_ABESS(lama)) &
              * IPT_P23_BBESS(lamb)*IPT_P23_PTPOW(n)
        enddo
     enddo
  enddo

  fctr=rkb**(l-1)
  do lamb=l,lmbhi
     IPT_P23_BBESS(lamb)=fctr/IPT_SH_DFAC(lamb+lamb+1)
     fctr=rkb*fctr
  enddo

  fctr=rka**(l-1)
  do lama=l,lmahi
     do lamb=l,lmbhi
        n=((1-l-l)+lama)+lamb
        IPT_P23_QSUM(n,lamb,lama) = IPT_P23_QSUM(n,lamb,lama)  &
           +(fctr/IPT_SH_DFAC(lama+lama+1)) *IPT_P23_BBESS(lamb) *IPT_P23_Q2(lamb,lama)
     enddo
     fctr=rka*fctr
  enddo
end subroutine ptwt


subroutine ssibfn(nmax,x,ssi)
! scaled spherical i bessel functions
  implicit real*8 (a-h,o-z)
  parameter (a0=0.0d0, a1=1.0d0, a3s2=1.5d0, a2=2.0d0, a3=3.0d0,a20=20.0d0)
  dimension ssi(*)
 
      x2=x*x
      xmin=(abs(3*nmax-1))
      if(x.gt.xmin) go to 5
      n=nmax
      f2np1=(n+n+1)
      f2kp3=f2np1
      pkm1=a0
      pk=a1
      qkm1=a1
      qk=a1
      aprod=a1
    1 f2kp1=f2kp3
      f2kp3=f2kp3+a2
      ak=x2/(f2kp1*f2kp3)
      akpkm2=ak*pkm1
      pkm1=pk
      pk=pkm1+akpkm2
      akqkm2=ak*qkm1
      qkm1=qk
      qk=qkm1+akqkm2
      aprod=ak*aprod
      if(((pk*qkm1)+aprod).ne.(pk*qkm1)) go to 1
      ssi(n+1)=pk/qk
    2 if(n.eq.0) go to 3
      n=n-1
      f2np3=f2np1
      f2np1=f2np1-a2
      ssi(n+1)=(f2np1*f2np3)/((f2np1*f2np3)+x2*ssi(n+2))
      go to 2
    3 ssi(1)=ssi(1)/(a1+x*ssi(1))
      do 4 n=1,nmax
        ssi(n+1)=ssi(n+1)*ssi(n)
    4 continue
      return
 
    5 if(x.ge.a20) then
        ex=a0
      else
        ex=exp(-(x+x))
      endif
      ssi(1)=(a1-ex)/(x+x)
      if(nmax.eq.0) return
      ssi(2)=a3s2*(a1+ex+(ex-a1)/x)/x2
      if(nmax.eq.1) return
      f2np1=a3
      do n=2,nmax
        f2nm1=f2np1
        f2np1=f2np1+a2
        ssi(n+1)=(ssi(n-1)-ssi(n))*(f2nm1*f2np1)/x2
      enddo
      return
end subroutine ssibfn


subroutine qpasy(alpha,npi,l,lambu,lmahi,lmbhi,ltot1,xka,xkb,prd,dum)
! compute type 2 radial integrals, scaled by exp(-arc2)/sqrt(pi),
! using the partially asymptotic method,
! for lama=l to lmahi, lamb=l to lmbhi, n=lama+lamb-l-l
  USE RECP_IPT
  USE RECP_FUNCTION2
  implicit real*8 (a-h,o-z)
  parameter (a0=0.0d0, accrcy=1.0d-13, a1s4=0.25d0, a1s2=0.5d0,a1=1.0d0)
 
  sqalpi=a1/sqrt(alpha)
  alf1=a1
  if(xka.gt.xkb) go to 42
 
! xka is smaller: set up parameters for qcomp using xkb
  xk=xkb*sqalpi
  t=a1s4*xk*xk
  prde=prd*exp(t-dum)*sqalpi**(npi+l)
  if(l.ge.2) prde=prde*xka**(l-1)
  tk=xka*xka/(alpha+alpha)
  do 30 lama=l,lmahi
  la=lama-1
  prefac=prde
  do 28 lamb=l,lmbhi
  lb=lamb-1
  n=((1-l-l)+lama)+lamb

! # run power series using xka, obtaining initial
! # q(n,l) values from qcomp, then recurring upwards
! # j=0 term in sum
  nprime=npi+n+la-1
  qold1=qcomp(alf1,IPT_SH_DFAC,nprime,lb,t,xk)/IPT_SH_DFAC(la+la+3)
  sum=qold1
  if(tk.eq.a0) go to 24

! # j=1 term in sum
  nprime=nprime+2
  qnew=qcomp(alf1,IPT_SH_DFAC,nprime,lb,t,xk)/IPT_SH_DFAC(la+la+3)
  f1=(la+la+3)
  qold2=(tk/f1)*qold1
  qold1=(tk/f1)*qnew
  sum=sum+qold1
  j=1

! # increment j for next term
   22 j=j+1
      nprime=nprime+2
      f1=(nprime+nprime-5)
      f2=((lb-nprime+4)*(lb+nprime-3))
      qnew=(t+a1s2*f1)*qold1+a1s4*f2*qold2
      f1=(j*(la+la+j+j+1))
      qold2=(tk/f1)*qold1
      qold1=(tk/f1)*qnew
      sum=sum+qold1
      if(qold1.gt.accrcy*sum) go to 22
   24 IPT_P23_QSUM(n,lamb,lama)=IPT_P23_QSUM(n,lamb,lama)+prefac*sum
      prefac=prefac*sqalpi
   28 continue
      prde=prde*(xka/alpha)
   30 continue
      return
 
!  xkb is smaller: set up parameters for qcomp using xka
   42 xk=xka*sqalpi
      t=a1s4*xk*xk
      prde=prd*exp(t-dum)*sqalpi**(npi+l)
      if(l.ge.2) prde=prde*xkb**(l-1)
      tk=xkb*xkb/(alpha+alpha)
      do 60 lama=l,lmahi
      la=lama-1
      prefac=prde
      do 58 lamb=l,lmbhi
      lb=lamb-1
      n=((1-l-l)+lama)+lamb
!     # run power series using xkb, obtaining initial
!     # q(n,l) values from qcomp, then recurring upwards
!     # j=0 term in sum
      nprime=npi+n+lb-1
      qold1=qcomp(alf1,IPT_SH_DFAC,nprime,la,t,xk)/IPT_SH_DFAC(lb+lb+3)
      sum=qold1
      if(tk.eq.a0) go to 54
!     # j=1 term in sum
      nprime=nprime+2
      qnew=qcomp(alf1,IPT_SH_DFAC,nprime,la,t,xk)/IPT_SH_DFAC(lb+lb+3)
      f1=(lb+lb+3)
      qold2=(tk/f1)*qold1
      qold1=(tk/f1)*qnew
      sum=sum+qold1
      j=1
!     # increment j for next term
   52 j=j+1
      nprime=nprime+2
      f1=(nprime+nprime-5)
      f2=((la-nprime+4)*(la+nprime-3))
      qnew=(t+a1s2*f1)*qold1+a1s4*f2*qold2
      f1=(j*(lb+lb+j+j+1))
      qold2=(tk/f1)*qold1
      qold1=(tk/f1)*qnew
      sum=sum+qold1
      if(qold1.gt.accrcy*sum) go to 52
   54 IPT_P23_QSUM(n,lamb,lama)=IPT_P23_QSUM(n,lamb,lama)+prefac*sum
      prefac=prefac*(xkb/alpha)
   58 continue
      prde=prde*sqalpi
   60 continue
      return
end subroutine qpasy

END MODULE RECP_CSO_RAD
