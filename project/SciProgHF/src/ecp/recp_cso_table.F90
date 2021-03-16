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

MODULE RECP_CSO_TABLE
CONTAINS
SUBROUTINE RECP_CSO_TABLE_HERMIT(nn,x,a,eps)
! calculates the zeros  x(i)  of the nn-th order hermite polynomial.  
! the largest zero will be stored in x(1).  also calculates the corresponding coefficients a(i)
! of the nn-th order gauss-hermite quadrature formula of degree 2*nn-1.  
! the factor of sqrt(pi) has been removed from the a(i).
! a. h. stroud & d. secrest, gaussian quadrature formulas, prentice-hall, 1966
  USE RECP_OUTPUT 
  implicit real*8 (a-h,o-z)
#include "inc_print.h"
  parameter (sixth = 1.0d0/6.0d0)
  dimension x(*), a(*)

  fn = (nn)
  n1 = nn - 1
  n2 = (nn+1)/2
  cc = 1.0d0
  s  = 0.0d0
  do i=1,n1
     s  = s + 0.5d0
     cc = s*cc
  enddo
  s  = (2.0d0*fn+1.0d0)**sixth
  do i=1,n2
     if( i.eq.1 ) then      ! largest zero
       xt = s**3 - 1.85575d0/s
     elseif( i.eq.2 ) then  ! second zero
       xt = xt - 1.14d0*fn**0.426d0/xt
     elseif( i.eq.3 ) then  ! third zero
       xt = 1.86d0*xt - 0.86d0*x(1)
     elseif( i.eq.4 ) then  ! fourth zero
       xt = 1.91d0*xt - 0.91d0*x(2)
     else                   ! all other zeros
       xt = 2.0d0*xt - x(i-2)
     endif
 
     call RECP_CSO_TABLE_HROOT(xt,nn,dpn,pn1,eps)
     x(i) = xt
     a(i) = cc/dpn/pn1
!    write (6,'(2i4,2d25.17)') nn, i, xt, a(i)
     ni = nn-i+1
     x(ni) = -xt
     a(ni) = a(i)
  enddo
  return
END SUBROUTINE RECP_CSO_TABLE_HERMIT


SUBROUTINE RECP_CSO_TABLE_HROOT(x,nn,dpn,pn1,eps)
! improves the approximate root  x
! in addition we also obtain
! dpn = derivative of h(n) at x
! pn1 = value of h(n-1) at x
  implicit real*8 (a-h,o-z)
! # iter = 5 sufficient for 8-byte accuracy up to nn = 7
  do iter=1,10
    call RECP_CSO_TABLE_HRECUR(p,dp,pn1,x,nn)
    d  = p/dp
    x  = x - d
    if(abs(d).le.eps) goto 16
  enddo
  
  16 CONTINUE
  dpn = dp
  return
END SUBROUTINE RECP_CSO_TABLE_HROOT


SUBROUTINE RECP_CSO_TABLE_HRECUR(pn,dpn,pn1,x,nn)
  implicit real*8 (a-h,o-z)
  p1 = 1.0d0
  p  = x
  dp1 = 0.0d0
  dp  = 1.0d0
  do j=2,nn
    fj = (j)
    fj2 = 0.5d0*(fj-1.0d0)
    q  = x*p - fj2*p1
    dq = x*dp + p - fj2*dp1
    p1 = p
    p  = q
    dp1 = dp
    dp  = dq
  enddo
  pn  = p
  dpn = dp
  pn1 = p1
  return
END SUBROUTINE RECP_CSO_TABLE_HRECUR


SUBROUTINE RECP_CSO_TABLE0(eps,lmax,lmn1u,lproju,ndfac)
! # tables for core potential and spin-orbit integrals.
  USE RECP_IPT
  USE RECP_FUNCTION1
  implicit real*8 (a-h,o-z)
#include "inc_print.h"
  parameter (a0=0.0d0, a1s2=0.5d0, a1=1.0d0, a2=2.0d0, a3=3.0d0)

! # compute gauss-hermite points and weights for c, z integrals.
  igh =1
  nn = 5
  do i = 1, 3
     CALL RECP_CSO_TABLE_HERMIT(nn,IPT_SH_HPT(igh),IPT_SH_HWT(igh),eps)
     igh = igh + nn
     nn = 2*nn
  enddo

! compute double factorials.
  CALL RECP_FN1_DFAC(ndfac)

! # compute binomial coefficients.
  inew=1
  IPT_SH_BINOM(1)=a1
  do j=1,lmn1u-1
    inew=inew+1
    IPT_SH_BINOM(inew)=a1
    do i=1,j-1
      inew=inew+1
      IPT_SH_BINOM(inew)=((j-i+1)*IPT_SH_BINOM(inew-1))/(i)
    enddo
    inew=inew+1
    IPT_SH_BINOM(inew)=a1
  enddo

! # compute tables by recursion for real spherical harmonics.  they
! # are indexed by l, m and sigma.  the sequence number of the
! # harmonic with quantum numbers l, m and sigma is given by
! #            l**2+2*m+1-sigma
! # IPT_SH_LMF(index) and IPT_SH_LML(index) hold the positions 
! # of the first and last terms of the harmonic in the 
! # IPT_SH_LMX, IPT_SH_LMY, IPT_SH_LMZ, and IPT_SH_ZLM arrays.
! # the harmonics with angular momentum l are generated from those
! # with angular momenta l-1 and l-2.
! # for m = 0,1,2,...,l-1, the recursion relation
! z*Z(l-1,m,s) = sqrt(((l-m)*(l+m))/((2*l-1)*(2*l+1)))*Z(l,m,s)+
!              sqrt(((l+m-1)*(l-m-1))/((2*l-3)*(2*l-1)))*Z(l-2,m,s)
! # is used.
! # for m = l, the recursion relation
! x*Z(l-1,l-1,s)+(-1)**(1-s)*y*Z(l-1,l-1,1-s) =
!              sqrt((2*l))/((2*l+1)))*Z(l,l,s)
! # is used.
! # l=0
  IPT_SH_LMF(1) = 1
  IPT_SH_LML(1) = 1
  IPT_SH_LMX(1) = 0
  IPT_SH_LMY(1) = 0
  IPT_SH_LMZ(1) = 0
  IPT_SH_ZLM(1) = a1
! # l=1
  IPT_SH_LMF(2) = 2
  IPT_SH_LML(2) = 2
  IPT_SH_LMX(2) = 0
  IPT_SH_LMY(2) = 0
  IPT_SH_LMZ(2) = 1
  IPT_SH_ZLM(2) = sqrt(a3)
  IPT_SH_LMF(3) = 3
  IPT_SH_LML(3) = 3
  IPT_SH_LMX(3) = 0
  IPT_SH_LMY(3) = 1
  IPT_SH_LMZ(3) = 0
  IPT_SH_ZLM(3) = IPT_SH_ZLM(2)
  IPT_SH_LMF(4) = 4
  IPT_SH_LML(4) = 4
  IPT_SH_LMX(4) = 1
  IPT_SH_LMY(4) = 0
  IPT_SH_LMZ(4) = 0
  IPT_SH_ZLM(4) = IPT_SH_ZLM(2)
  nterm=4
      do 270 lang=2,lmax
        do 240 mang=0,lang-1
          anum = ((2*lang-1)*(2*lang+1))
          aden = ((lang-mang)*(lang+mang))
          coef1 = sqrt(anum/aden)
          anum = ((lang+mang-1)*(lang-mang-1)*(2*lang+1))
          aden = (2*lang-3)*aden
          coef2 = sqrt(anum/aden)
          nsigma=min(1,mang)
          do 230 isigma=nsigma,0,-1
            indexh=lang**2+2*mang+1-isigma
            lone=lang-1
            ltwo=lang-2
            ione=lone**2+2*mang+1-isigma
            itwo=ltwo**2+2*mang+1-isigma
            IPT_SH_LMF(indexh)=IPT_SH_LML(indexh-1)+1
            IPT_SH_LML(indexh)=IPT_SH_LML(indexh-1)
            nxy=(mang-isigma+2)/2
            iu=IPT_SH_LMF(ione)+nxy-1
            do i=IPT_SH_LMF(ione),iu
               IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
               j=IPT_SH_LML(indexh)
               IPT_SH_LMX(j)=IPT_SH_LMX(i)
               IPT_SH_LMY(j)=IPT_SH_LMY(i)
               IPT_SH_LMZ(j)=IPT_SH_LMZ(i)+1
               IPT_SH_ZLM(j)=IPT_SH_ZLM(i)*coef1
               nterm=nterm+1
            enddo
            if(ltwo.ge.mang) then
              il=iu+1
              do i=il,IPT_SH_LML(ione)
                 IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
                 j=IPT_SH_LML(indexh)
                 k=IPT_SH_LMF(itwo)+i-il
                 IPT_SH_LMX(j)=IPT_SH_LMX(k)
                 IPT_SH_LMY(j)=IPT_SH_LMY(k)
                 IPT_SH_LMZ(j)=IPT_SH_LMZ(k)
                 IPT_SH_ZLM(j)=IPT_SH_ZLM(i)*coef1-IPT_SH_ZLM(k)*coef2
                 nterm=nterm+1
              enddo
              il=IPT_SH_LML(itwo)-nxy+1
              if(mod(lang-mang,2).eq.0) then
                do i=il,IPT_SH_LML(itwo)
                   IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
                   j=IPT_SH_LML(indexh)
                   IPT_SH_LMX(j)=IPT_SH_LMX(i)
                   IPT_SH_LMY(j)=IPT_SH_LMY(i)
                   IPT_SH_LMZ(j)=IPT_SH_LMZ(i)
                   IPT_SH_ZLM(j)=-IPT_SH_ZLM(i)*coef2
                   nterm=nterm+1
                enddo
              endif
            endif
  230     end do
  240   end do
        anum = (2*lang+1)
        aden = (2*lang)
        coef = sqrt(anum/aden)
        mang=lang
        isigma=1
        indexh=lang**2+2*mang+1-isigma
        IPT_SH_LMF(indexh)=IPT_SH_LML(indexh-1)+1
        IPT_SH_LML(indexh)=IPT_SH_LML(indexh-1)
!       # isig:  index of the harmonic (l-1),(m-1),sigma
!       # isigm: index of the harmonic (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        k=IPT_SH_LMF(isigm)
        do i=IPT_SH_LMF(isig),IPT_SH_LML(isig)
           IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
           j=IPT_SH_LML(indexh)
           IPT_SH_LMX(j)=IPT_SH_LMX(i)+1
           IPT_SH_LMY(j)=IPT_SH_LMY(i)
           IPT_SH_LMZ(j)=IPT_SH_LMZ(i)
           IPT_SH_ZLM(j)=(IPT_SH_ZLM(i)+IPT_SH_ZLM(k))*coef
           k=k+1
           nterm=nterm+1
        enddo
        if(mod(mang,2).eq.1) then
          IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
          j=IPT_SH_LML(indexh)
          IPT_SH_LMX(j)=IPT_SH_LMX(k)
          IPT_SH_LMY(j)=IPT_SH_LMY(k)+1
          IPT_SH_LMZ(j)=IPT_SH_LMZ(k)
          IPT_SH_ZLM(j)=IPT_SH_ZLM(k)*coef
          nterm=nterm+1
        endif
        isigma=0
        indexh=lang**2+2*mang+1-isigma
!       # isig:  index of the harmonic (l-1),(m-1),sigma
!       # isigm: index of the harmonuc (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        IPT_SH_LMF(indexh)=IPT_SH_LML(indexh-1)+1
        IPT_SH_LML(indexh)=IPT_SH_LMF(indexh)
        j=IPT_SH_LML(indexh)
        i=IPT_SH_LMF(isig)
        IPT_SH_LMX(j)=IPT_SH_LMX(i)+1
        IPT_SH_LMY(j)=IPT_SH_LMY(i)
        IPT_SH_LMZ(j)=IPT_SH_LMZ(i)
        IPT_SH_ZLM(j)=IPT_SH_ZLM(i)*coef
        nterm=nterm+1
        k=IPT_SH_LMF(isigm)
        do i=IPT_SH_LMF(isig)+1,IPT_SH_LML(isig)
          IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
          j=IPT_SH_LML(indexh)
          IPT_SH_LMX(j)=IPT_SH_LMX(i)+1
          IPT_SH_LMY(j)=IPT_SH_LMY(i)
          IPT_SH_LMZ(j)=IPT_SH_LMZ(i)
          IPT_SH_ZLM(j)=(IPT_SH_ZLM(i)-IPT_SH_ZLM(k))*coef
          k=k+1
          nterm=nterm+1
        enddo
        if(mod(mang,2).eq.0) then
          IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
          j=IPT_SH_LML(indexh)
          k=IPT_SH_LML(isigm)
          IPT_SH_LMX(j)=IPT_SH_LMX(k)
          IPT_SH_LMY(j)=IPT_SH_LMY(k)+1
          IPT_SH_LMZ(j)=IPT_SH_LMZ(k)
          IPT_SH_ZLM(j)=-IPT_SH_ZLM(k)*coef
          nterm=nterm+1
        endif
  270 end do

   ixy = 0
   iz = 0
   do lang=1,lproju
      do mang=0,lang-1
         nsigma=min(1,mang)
         ndelta=max(0,1-mang)
         anum = ((lang-mang)*(lang+mang+1))
         aden = (2*(2-ndelta))
         coef=sqrt(anum/aden)
         do isigma=nsigma,0,-1
            isign=2*isigma-1
            ixy = ixy+1
            IPT_SH_FLMTX(1,ixy) = (isign)*coef
            IPT_SH_FLMTX(2,ixy) = coef
            if(mang.ne.0) then
               iz=iz+1
               IPT_SH_FLMTX(3,iz) = -(mang*isigma)
            endif
         enddo
      enddo
      iz=iz+1
      IPT_SH_FLMTX(3,iz) = -(lang)
   enddo

! # column and row indices for angular momentum matrix elements.
  iadd = 1
  do i=1,2*lproju-1
     IPT_SH_MC(1,i) = i
     IPT_SH_MC(2,i) = i
     IPT_SH_MC(3,i) = i+1
     IPT_SH_MR(1,i) = i+iadd
     IPT_SH_MR(2,i) = i+2
     IPT_SH_MR(3,i) = i+2
     iadd = 4-iadd
  enddo
  return
END SUBROUTINE RECP_CSO_TABLE0

END MODULE RECP_CSO_TABLE
