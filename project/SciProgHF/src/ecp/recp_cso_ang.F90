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

MODULE RECP_CSO_ANG
CONTAINS
SUBROUTINE RECP_CSO_ANG1(nanb,lalb,mamb,lamu,ltot1,xk,yk,zk)
! compute type 1 angular integrals
  USE RECP_IPT 
  USE RECP_OUTPUT
  implicit real*8 (a-h,o-z)
#include "inc_print.h"
  parameter (a0=0.0d0, a1=1.0d0)
  integer i,j

  CALL OUTPUT_LOC('ANG1','E')
  do i = 1,ltot1 
!    do j = 1,lamu
     do j = 1,ltot1 
        IPT_P1_ANG(i,j) = 0.0d0
     enddo
  enddo

      do 96 n=1,nanb
      if(IPT_P1_XAB(n).eq.a0) go to 96
      do 94 l=1,lalb
      if(IPT_P1_YAB(l).eq.a0) go to 94
      do 92 m=1,mamb
      if(IPT_P1_ZAB(m).eq.a0) go to 92
      nlm=((n-2)+l)+m
      lamlo=mod(nlm-1,2)+1
      lamhi=min(nlm,lamu)
      if(lamlo.gt.lamhi) go to 92
      do 90 lam=lamlo,lamhi,2
      l2=lam+lam-1
      angt=a0
      loc=(lam-1)**2
      do 80 mu1=1,l2
      istart=IPT_SH_LMF(loc+mu1)
      if(mod(n,2).eq.mod(IPT_SH_LMX(istart),2).or.  &
         mod(l,2).eq.mod(IPT_SH_LMY(istart),2).or.  &
         mod(m,2).eq.mod(IPT_SH_LMZ(istart),2)) go to 80
      pre=a0
      aint=a0
      iend=IPT_SH_LML(loc+mu1)
      do 70 i=istart,iend
        indx=IPT_SH_LMX(i)
        indy=IPT_SH_LMY(i)
        indz=IPT_SH_LMZ(i)
        if(indx.eq.0) then
          xkp=a1
        else
          xkp=xk**indx
        endif
        if(indy.eq.0) then
          ykp=a1
        else
          ykp=yk**indy
        endif
        if(indz.eq.0) then
          zkp=a1
        else
          zkp=zk**indz
        endif
        pre=pre+IPT_SH_ZLM(i)*xkp*ykp*zkp
        aint = aint + IPT_SH_ZLM(i)        & 
                    * IPT_SH_DFAC(n+indx)  &
                    * IPT_SH_DFAC(l+indy)  &
                    * IPT_SH_DFAC(m+indz)  &
                    / IPT_SH_DFAC((n+indx)+(l+indy)+(m+indz))
   70 continue
      angt=angt+pre*aint
   80 continue
      IPT_P1_ANG(nlm,lam)=IPT_P1_ANG(nlm,lam)+((IPT_P1_XAB(n)*IPT_P1_YAB(l))*IPT_P1_ZAB(m))*angt
   90 continue
   92 continue
   94 continue
   96 continue
  CALL OUTPUT_LOC('ANG1','X')
  return
END SUBROUTINE RECP_CSO_ANG1


SUBROUTINE RECP_CSO_ANG2(ang,crda,it,l,lit,lmlo,lmhi,lmnv,mproju,xk,yk,zk)
! compute type 2 angular integrals
  USE RECP_IPT
  implicit real*8 (a-h,o-z)
  parameter (a0=0.0d0, a1=1.0d0, a2=2.0d0, a3=3.0d0)
  dimension ang(lit,mproju,*),crda(lit,3),lmnv(3,*)
  integer i,j,k
 
  do i = 1, lit
     do j = 1, mproju
        do k = 1, lmhi
           ang(i,j,k) = 0.0d0 
        enddo
     enddo
  enddo

!     print *,'ang2-it',it
      na1=lmnv(1,it)+1
      la1=lmnv(2,it)+1
      ma1=lmnv(3,it)+1
      naind=(na1*(na1-1))/2
      laind=(la1*(la1-1))/2
      maind=(ma1*(ma1-1))/2
      loc1=(l-1)**2
      mhi=l+l-1
      do 80 ia=1,na1
      pab1=IPT_SH_BINOM(naind+ia)*crda((na1+1)-ia,1)
      if(pab1.eq.a0) go to 80
      do 70 ib=1,la1
      pab2=pab1*IPT_SH_BINOM(laind+ib)*crda((la1+1)-ib,2)
      if(pab2.eq.a0) go to 70
      do 60 ic=1,ma1
      pab3=pab2*IPT_SH_BINOM(maind+ic)*crda((ma1+1)-ic,3)
      if(pab3.eq.a0) go to 60
      n=((ia-3)+ib)+ic
      lamlo=max(l-n,lmlo+mod(l+n+lmlo,2))
      lamhi=min(l+n,lmhi-mod(l+n+lmhi,2))
      if(lamlo.gt.lamhi) go to 60
      do 55 m=1,mhi
      mstart=IPT_SH_LMF(loc1+m)
      mend=IPT_SH_LML(loc1+m)
      do 50 lam=lamlo,lamhi,2
      l2=lam+lam-1
      angt=a0
      loc2=(lam-1)**2
      do 40 mu=1,l2
      istart=IPT_SH_LMF(loc2+mu)
      if(mod(ia+IPT_SH_LMX(mstart)+IPT_SH_LMX(istart),2).ne.1.or.  &
         mod(ib+IPT_SH_LMY(mstart)+IPT_SH_LMY(istart),2).ne.1.or.  &
         mod(ic+IPT_SH_LMZ(mstart)+IPT_SH_LMZ(istart),2).ne.1) goto 40
      pre=a0
      iend=IPT_SH_LML(loc2+mu)
      aint=a0
      DO i=istart,iend
         indx=IPT_SH_LMX(i)
         indy=IPT_SH_LMY(i)
         indz=IPT_SH_LMZ(i)
         if(indx.eq.0) then
           xkp=a1
         else
           xkp=xk**indx
         endif
         if(indy.eq.0) then
           ykp=a1
         else
           ykp=yk**indy
         endif
         if(indz.eq.0) then
           zkp=a1
         else
           zkp=zk**indz
         endif
         pre=pre+IPT_SH_ZLM(i)*xkp*ykp*zkp
         DO j=mstart,mend
            mndx=IPT_SH_LMX(j)
            mndy=IPT_SH_LMY(j)
            mndz=IPT_SH_LMZ(j)
            aint=aint + IPT_SH_ZLM(i) * IPT_SH_ZLM(j)  &
                      * IPT_SH_DFAC(ia+indx+mndx)      &  
                      * IPT_SH_DFAC(ib+indy+mndy)      &
                      * IPT_SH_DFAC(ic+indz+mndz)      &
                      / IPT_SH_DFAC(ia+indx+mndx+ib+indy+mndy+ic+indz+mndz)
         ENDDO
      ENDDO
      angt=angt+pre*aint
   40 continue
      ang(n+1,m,lam)=ang(n+1,m,lam)+angt*pab3
   50 continue
   55 continue
   60 continue
   70 continue
   80 continue
      return
END SUBROUTINE RECP_CSO_ANG2
END MODULE RECP_CSO_ANG
