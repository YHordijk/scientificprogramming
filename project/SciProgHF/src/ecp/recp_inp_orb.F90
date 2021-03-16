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

MODULE RECP_INP_ORBITAL
CONTAINS

SUBROUTINE RECP_INP_INTCHECK(nu,ns,ncrs,mcrs,lcr,lls,inam)
! determine the number of types of one-electron integrals
!   inam=4 : arep
!   inam=5 : sorep
  IMPLICIT NONE
  INTEGER :: I,nu,ns,ncrs,mcrs(*),lcr(*),lls(*),inam(5)

  nu = 3 !initialize nu for s(*),t(*),v(*)

  if ( ncrs .gt. 0 ) then
    do I = 1,ns
      if ( mcrs(I) .ne. 0 ) then
        if ( lcr(mcrs(I)) .gt. 0 ) then
          nu = 4
          inam(4) = 4
          go to 200
        endif
      endif
    enddo
    200 continue
    do I = 1, ns
      if ( mcrs(I) .ne. 0 ) then
        if ( lls(mcrs(I)) .gt. 0 ) then
          nu = nu + 1
          inam(nu) = 5
          go to 208
        endif
      endif
    enddo
    208 continue
  endif
END SUBROUTINE RECP_INP_INTCHECK

SUBROUTINE RECP_INP_ORBITAL_CCNORM(mrcru,mconu,ncons,nrcr,ncon,lmnp1,eta,zet)
! normalize the contraction coefficients
  IMPLICIT NONE
#include "inc_print.h"
! variables : local 
  INTEGER :: mrcru,mconu ! from common (parmi)
  INTEGER :: ncons    ! I  Number of basis-set blocks
  INTEGER :: nrcr(*)  ! J  Number of contraction coefficient
  INTEGER :: ncon(*)  ! K  Number of primitive
  INTEGER :: lmnp1(*) !
  REAL(8) :: eta(mrcru,mconu,*)
  REAL(8) :: zet(mconu,*)
! variables : local 
  INTEGER :: I, J, K, L
  REAL(8) :: prtint, temp

  DO I = 1, ncons
    DO J = 1, nrcr(I)
!     make eta=1 if ncon=1
      if (ncon(I).eq.1) then
         eta(1,1,I) = 1.0d0
         GOTO 312               
      endif

      prtint = eta(J,1,I)**2
      DO K = 2, ncon(I)
         DO L = 1, (K-1)
            temp = (zet(K,I) + zet(L,I)) * 0.5d0
            temp = zet(K,I)/temp * zet(L,I)/temp
            prtint = prtint + 2.0d0 * eta(J,K,I)*eta(J,L,I) * sqrt(temp**lmnp1(I) * sqrt(temp))
         ENDDO
         prtint = prtint + eta(J,K,I)**2
      ENDDO

      prtint = sqrt(prtint)
      DO K = 1, ncon(I)
         eta(J,K,I) = eta(J,K,I) / prtint
      ENDDO

      312 CONTINUE !yc-check
    ENDDO
  ENDDO

! print
! IF (RECP_DBG.GE.5) THEN
!    DO I = 1, ncons
!       print *,'From RECP_INP(CCNORM)'
!       DO J = 1, ncon(I)
!          print *, zet(J,I),(eta(K,J,I),K=1,nrcr(I))
!       ENDDO
!    ENDDO
! ENDIF
END SUBROUTINE RECP_INP_ORBITAL_CCNORM



SUBROUTINE RECP_INP_ORBITAL_OFFSET(nsopr,nso,nblpr,nst,FILE_NUM)
! compute contracted orbital and lower-triangle-packed matrix symmetry offsets.
  USE RECP_FUNCTION1
  IMPLICIT NONE
#include "inc_print.h"
! variable : global
  INTEGER  :: nsopr(*),nso(*),nblpr(*),nst,FILE_NUM(*)
! variable : local
  INTEGER  :: i 

  IF(RECP_DBG.GE.1) WRITE(FILE_NUM(1),'(/,A)') ' *** RECP_INP_ORBITAL_OFFSET ***'
  IF(RECP_DBG.GE.1) WRITE(FILE_NUM(1),'(6X,2A)') 'I',':  nsopr  nblpr'
  nsopr(1) = 0
  nblpr(1) = 0
  DO I = 2,nst
     nsopr(i) = nsopr(i-1) + nso(i-1)
     nblpr(i) = nblpr(i-1) + ( (nso(i-1)+1)*((nso(i-1)+1)-1)/2 )
     IF(RECP_DBG.GE.1) WRITE(FILE_NUM(1),'(I7,A,2I7)') I,':',nsopr(i),nblpr(i)
  ENDDO
END SUBROUTINE RECP_INP_ORBITAL_OFFSET


SUBROUTINE RECP_INP_ORBITAL_ANGSYM(icu,itl,itu,ICS,igcs,lai,ixyzir)
! # determine angular momentum symmetries.
  USE RECP_IPT
  USE RECP_FUNCTION2
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER :: icu,itl,itu,ICS,igcs,lai,ixyzir(3)
! local variables
  REAL*8  :: c_sum 
  INTEGER :: I1,I2,I3,ict,ic,it
  LOGICAL :: L1,L2,L3
  DO I1 = 1, 3
     I2 = mod(I1,  3) + 1
     I3 = mod(I1+1,3) + 1
     c_sum = 0.0d0
     ict = 0
     DO ic = 1, icu
        DO it = itl,itu
           ict = ict + 1
           L1 = mod(IPT_LMNV(I1,it),2).eq.1
           L2 = mod(IPT_LMNV(I2,it),2).eq.0
           L3 = mod(IPT_LMNV(I3,it),2).eq.0
           IF(L1.AND.L2.AND.L3) c_sum=c_sum+IPT_AO2SO(ict,ICS,igcs) 
        ENDDO
     ENDDO

     IF ( (c_sum.NE.0.0d0).AND.(lai.NE.ixyzir(I1)) ) THEN 
        IF ( ixyzir(I1).NE.0 ) THEN
           WRITE (RECP_OUT,'(3X,A)') 'axis choice is unsuitable for spin-orbit integrals'
           CALL bummer('RECP_INP_ORBANGSYM: ixyzir=',ixyzir(I1),2)
        ENDIF
        ixyzir(I1) = lai
     ENDIF
  ENDDO
END SUBROUTINE RECP_INP_ORBITAL_ANGSYM


SUBROUTINE RECP_INP_ORBITAL_MULTNORMCOEF(ncons,ncon,nrcr,mrcru,mconu,zet,eta,lmnp1)
! # multiply normalization constants into contraction coefficients.
  IMPLICIT NONE
#include "inc_print.h"
! global variables
  INTEGER ncons,ncon(*),nrcr(*),mrcru,mconu,lmnp1(*) 
  REAL*8  zet(mconu,*),eta(mrcru,mconu,*)
! local variables
  INTEGER I,J,K
  REAL*8  tmp
  
  DO I = 1,ncons
     DO K = 1,ncon(I)
        DO J = 1,nrcr(I)
           tmp = 4.0d0 * zet(K,I)
           eta(J,K,I) = SQRT( tmp**lmnp1(I) /2.0d0*SQRT(tmp/2.0d0) ) * eta(J,K,I)
        ENDDO   
!       print *,zet(K,I),(eta(J,K,I),J=1,nrcr(I))
     ENDDO   
  ENDDO   
END SUBROUTINE RECP_INP_ORBITAL_MULTNORMCOEF


SUBROUTINE RECP_INP_ORBITAL_SETLXYZIR(nu,inam,nst,ixyzir,lxyzir,mstu,idp,ityp)
! set lxyzir
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER nu,inam(5),nst,ixyzir(3),lxyzir(3),mstu,idp(mstu,mstu,*)
  CHARACTER*3 ityp(*)
! local variables
  INTEGER I,I2,I3,J 
  IF (inam(nu).EQ.5) THEN
     DO I = 1, 3
        I2 = mod(I,  3) + 1
        I3 = mod(I+1,3) + 1
        DO J = 1, nst
           IF ( idp(J,ixyzir(I2),ixyzir(I3)) .NE. 0 ) THEN
              lxyzir(I) = J
              GOTO 860
           ENDIF
        ENDDO
        860 CONTINUE
     ENDDO

     IF (PRINT_LEVEL.GE.5) THEN
        WRITE (RECP_OUT,'(/,3(3X,A,A3))') &
          'lx: ',ityp(lxyzir(1)), 'ly: ',ityp(lxyzir(2)), 'lz: ',ityp(lxyzir(3))
     ENDIF
  ENDIF
END SUBROUTINE RECP_INP_ORBITAL_SETLXYZIR


SUBROUTINE RECP_INP_ORBITAL_L2OPXV(lmn,itl,itu,lmnv,eval,leig,ics,igcs)
!  compute the matrix vector product, (L^2) * v, in an unnormalized cartesian basis, 
!  and determine if v(*) is an eigenvector of the total angular momentum operator.
!  input:
!  lmn = (l + m + n) where l, m, and n are the exponents of x, y, and z
!        respectively.  lmn=0 for cartesian "s" functions, 1 for
!        cartesian "p" functions, 2 for cartesian "d" functions, etc.
!  v(1:xyzdim) = input vector.  the elements are coefficients of the
!              cartesian basis functions (x**l * y**m * z**n).
!              where xyzdim = ((lmn + 1) * (lmn + 2))/2
!  lmnv(1:3,1:lmnvmx) = cartesian exponents of each basis function.
!  scr(1:lenscr) = scratch work vector.  This must be at least as large
!                  as the cartesian subspace dimension.
!                  lenscr >= xyzdim
!  output:
!  eval = vector expectation value = <v| L^2 |v> / <v|v>
!  leig = eigenvector return code.
!       = -4 for vector norm error.
!       = -3 for lmnv(*,*) inconsistency.
!       = -1 if v(*) is not an eigenvector of L^2
!       = principal quantum number if v(*) is an eigenvector.
!         0 for s-type vectors, 1 for p-type vectors, etc.
!         the eigenvalue is given by (leig*(leig+1)) to within numerical
!         precision.
 
  USE RECP_IPT
  USE RECP_FUNCTION1
  USE RECP_FUNCTION2
  implicit logical(a-z)
  integer itl, itu, lmn, leig, lmnv(3,*)
  real*8 eval
  integer ics,igcs 
! # local:
  integer l, m, n, i, j, it, jt, xyzdim, ll, mm, nn, p, q
  real*8 vi, rnorm, vnorm2
! # small is used to determine if a vector is close enough to be
! #       called an eigenvector.  this should be about 1000x the
! #       machine epsilon to avoid sqrt() precision problems.
  real*8    zero,     one,     two,     fourth,        small
  parameter(zero=0d0, one=1d0, two=2d0, fourth=0.25d0, small=1d-10)
 
! # the [l,m,n] basis function is assigned the internal index:
! #  (lx*(lx+1))/2 + n + 1   with lx=(lmn-l)=(m+n)
! # this convention is inconsistent with the rest of this program,
! # but it is used anyway since it is convenient and since the
! # computed results are local.
  integer llp1
  llp1(l) = (l * (l + 1)) / 2 + 1
 
  leig = -1
  eval = -one
 
! # check to make sure scr(*) is large enough.
  xyzdim = llp1( lmn) + lmn
  IF(ALLOCATED(IPT_L2OPXV_SCR)) DEALLOCATE(IPT_L2OPXV_SCR)
  ALLOCATE(IPT_L2OPXV_SCR(xyzdim))

! # initialize the scratch vector.
  CALL RECP_SETZERO_R1(IPT_L2OPXV_SCR,xyzdim)
 
! # compute the matrix-vector product and the vector norm.
! # matrix elements of lop(*,*) are computed from the operator
! # definition in cartesian space.
  vnorm2 = zero
  do it = itl, itu
     i = it - (itl - 1)
     l = lmnv(1,it)
     m = lmnv(2,it)
     n = lmnv(3,it)
!    vi = v(i)
     vi = IPT_AO2SO(i,ics+1,igcs)
     do jt = itl, it-1
        j = jt - (itl - 1)
!       vnorm2 = vnorm2 + (two * vi) * v(j) *
        vnorm2 = vnorm2 + (two * vi) * IPT_AO2SO(j,ics+1,igcs) &
                        * iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
!       print *,'v(',i,')',vi,'v(',j,')',v(j)
     enddo

     vnorm2 = vnorm2 + vi * vi * iodfac(l,l,m,m,n,n)
!    print *,'vnorm2',vnorm2,'vi',vi,'iodfac',iodfac(l,l,m,m,n,n)
     if ( (l+m+n) .ne. lmn ) then
!       # inconsistent exponents in the basis function.
        call bummer('l2opxv: (l+m+n-lmn)=', (l+m+n-lmn), 0 )
        leig = -3
        IF(ALLOCATED(IPT_L2OPXV_SCR)) DEALLOCATE(IPT_L2OPXV_SCR)
        return
     endif
 
!    # L^2 is sparse in this representation.  each vi contributes
!    # to, at most, only 7 elements in the matrix-vector product.

     if ( l .ge. 2 ) then
!       # include -l*(l-1) * ( [l-2,m,n+2] + [l-2,m+2,n] ) terms.
        ll = l * (l - 1)
        p = llp1( lmn - l + 2 ) + n
        IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (ll)
        p = p + 2
        IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (ll)
     endif
 
     if ( m .ge. 2 ) then
!       # include -m*(m-1) * ( [l,m-2,n+2] + [l+2,m-2,n] ) terms.
        mm = m * (m - 1)
        p = llp1( lmn - l ) + n + 2
        IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (mm)
        p = llp1( lmn - l - 2 ) + n
        IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (mm)
     endif
 
     if ( n .ge. 2 ) then
!       # include -n*(n-1) * ( [l,m+2,n-2] + [l+2,m,n-2] ) terms.
        nn = n * (n - 1)
        p = llp1( lmn - l ) + n - 2
        IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (nn)
        p = llp1( lmn - l - 2 ) + n - 2
        IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (nn)
     endif
 
!    # include the 2*(l*m+l*n+m*n+l+m+n)*[l,m,n] diagonal term.
     p = llp1( lmn - l ) + n
     IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) + vi * (2 * (l*(m+n) + m*n + l + m + n))
  enddo
 
  if ( vnorm2 .le. small ) then
     call bummer('l2mxv: small vector norm, lnm=', lmn, 0)
     leig = -4
     IF(ALLOCATED(IPT_L2OPXV_SCR)) DEALLOCATE(IPT_L2OPXV_SCR)
     return
  endif
 
! # compute the expectation value.
  eval  = zero
  do it = itl, itu
     i = it - (itl - 1)
     l = lmnv(1,it)
     m = lmnv(2,it)
     n = lmnv(3,it)
     p = llp1( lmn - l ) + n
     do jt = itl, itu
        j = jt - (itl - 1)
!       eval = eval + v(j) * IPT_L2OPXV_SCR(p) *
        eval = eval + IPT_AO2SO(j,ics+1,igcs) * IPT_L2OPXV_SCR(p) &
                    * iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
     enddo
  enddo
  eval = eval / vnorm2
 
! # compute the residual norm.
  rnorm = zero
 
  do it = itl, itu
     i = it - (itl - 1)
     l = lmnv(1,it)
     m = lmnv(2,it)
     n = lmnv(3,it)
     p = llp1( lmn - l ) + n
     do jt = itl, it-1
        j = jt - (itl - 1)
        q = llp1( lmn - lmnv(1,jt) ) + lmnv(3,jt)
!       rnorm = rnorm + two * ( IPT_L2OPXV_SCR(p) - eval * v(i) ) 
!                           * ( IPT_L2OPXV_SCR(q) - eval * v(j) ) 
!                           * iodfac
        rnorm = rnorm + two * ( IPT_L2OPXV_SCR(p) - eval * IPT_AO2SO(i,ics+1,igcs) ) &
                            * ( IPT_L2OPXV_SCR(q) - eval * IPT_AO2SO(j,ics+1,igcs) ) &
                            * iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
     enddo
!   rnorm = rnorm + (IPT_L2OPXV_SCR(p)-eval*v(i))**2 
    rnorm = rnorm + (IPT_L2OPXV_SCR(p)-eval*IPT_AO2SO(i,ics+1,igcs))**2 *iodfac(l,l,m,m,n,n)
  enddo
 
! # normalize w.r.t. |v|=1.
  rnorm = rnorm / vnorm2
  rnorm = sqrt( rnorm )

  if ( rnorm .gt. small ) then
!    # v(*) is not an eigenvector.
     leig = -1
  else
!    # v(*) is an eigenvector.
!    # determine leig such that eval = leig*(leig+1)
!    # the following assignment should truncate to the
!    # next smaller integer value...
     leig = ( sqrt( eval + fourth ) )
  endif
  IF(ALLOCATED(IPT_L2OPXV_SCR)) DEALLOCATE(IPT_L2OPXV_SCR)
  return
END SUBROUTINE RECP_INP_ORBITAL_L2OPXV


function iodfac(l1,l2,m1,m2,n1,n2)
! compute the product of odd number factorials which gives the
! one-center overlaps of cartesian gaussian aos.
  if(mod(l1+l2,2).eq.1 .or.  mod(m1+m2,2).eq.1 .or. mod(n1+n2,2).eq.1) then
    iodfac = 0
  else
    iprd = 1
    do i=3,l1+l2-1,2
      iprd = i*iprd
    enddo
    do i=3,m1+m2-1,2
      iprd = i*iprd
    enddo
    do i=3,n1+n2-1,2
      iprd = i*iprd
    enddo
    iodfac = iprd
  endif
  return
end function iodfac



! --------------- SOCFPD ---------------
SUBROUTINE RECP_SOCFPD_SETNPR(nst,npr,isfr,jsfr,isfrib,jsfrjb,iru,jru,iaords,jaords,ipair, &
           ijsf,mru,la,lb,mstu,idp,nsopr,mau,MAXICS,MAXJCS,nd,nprir,npair,ilfact,esfb)
  USE RECP_IPT
  USE RECP_FUNCTION2
  IMPLICIT NONE
#include "inc_mxvalue.h"
#include "inc_print.h"
! global variable
  INTEGER nst,npr,isfr,jsfr,isfrib,jsfrjb,iru,jru,iaords,jaords,ipair,ijsf
  INTEGER mru,la(mru,*),lb(*),mstu,idp(mstu,mstu,*),nsopr(*),mau(*),MAXICS(*),MAXJCS(*),nd(*)
  INTEGER nprir(2,mstu,*),npair(2,*),ilfact
  LOGICAL esfb
! local variable
  INTEGER ist,ipr,ics,jcs,ir,jr,la1,la2,npru,iesfb

 
  if(esfb) then
    iesfb = 1
    npru = (iru * (iru + 1)) / 2
  else
    iesfb = 2
    npru = iru * jru
  endif
 
  if (npru.gt.mnrup) call bummer('change mnrup (one place) to ',npru,2)

  npr = 0 
  DO ist = 1, nst
     ipr = npr
     ics = 0
     isfr = isfrib
     DO ir = 1, iru
        isfr = isfr + 1
        la1 = la(ir,iaords)
        IF(esfb) jru = ir
        jsfr = jsfrjb
        jcs = 0
        DO jr = 1, jru
           jsfr = jsfr + 1
           la2 = la(jr,jaords)
!          print *,'socfpd: idp',la1,la2,ist,'=',idp(la1,la2,ist)
           IF ( idp(la1,la2,ist).NE.0 ) THEN
              ipair = ipair + 1
              IPT_IL(ipair) = (nsopr(la1)+lb(isfr))*ilfact + (nsopr(la2)+lb(jsfr))
 
              npr = npr + 1
              mau(npr)  = nd(la1)
              MAXICS(npr) = ics
              MAXJCS(npr) = jcs
           ENDIF
           jcs = jcs + nd(la2)
        ENDDO
        ics = ics + nd(la1)
     ENDDO
     nprir(iesfb,ist,ijsf) = npr - ipr
  ENDDO
  npair(iesfb,ijsf) = npr
END SUBROUTINE RECP_SOCFPD_SETNPR


SUBROUTINE RECP_SOCFPD_SETICXU(RETURNVALUE,igcs,jgcs,ijsf,icx,icu,icxu,isf,jsf, &
           npr,MAXCX,icxa,ngcs,icxsv1,icxsv2,icxst,icxast,nt,jcu,esfb,SO_CALC)
  USE RECP_FUNCTION2
  IMPLICIT NONE
! global variable
  INTEGER RETURNVALUE,igcs,jgcs,ijsf,icx,icu,icxu,isf,jsf,npr,MAXCX,icxa,ngcs
  INTEGER icxsv1(ngcs),icxsv2(ngcs,ngcs)
  INTEGER icxst(2,*),icxast(*),nt(*),jcu
  LOGICAL esfb,SO_CALC
! local variable
  INTEGER isv1,isv2,SUM0

  RETURNVALUE = 0  !set initial return value  

  if (esfb) then
     isv1 = icxsv1(igcs)
     if (isv1.ne.0) then
        icxst(1,ijsf) = icxst(1,isv1)
        if(SO_CALC) icxast(ijsf) = icxast(isv1)
        RETURNVALUE = 1
        RETURN
     endif

     icxsv1(igcs) = ijsf
     icxst(1,ijsf) = icx
     SUM0 = (icu*nt(isf)+1)*((icu*nt(isf)+1)-1)/2
     icxu = icx+(npr*SUM0)
     if (SO_CALC) then
        MAXCX = MAXCX-(npr*SUM0)
        icxa = MAXCX 
        icxast(ijsf) = icxa
     endif
  else
     isv2 = icxsv2(jgcs,igcs)
     if (isv2.ne.0) then
        icxst(2,ijsf) = icxst(2,isv2)
        RETURNVALUE = 1
        RETURN
     endif
     icxsv2(jgcs,igcs) = ijsf
     icxst(2,ijsf) = icx
     SUM0 = icu*nt(isf)*jcu*nt(jsf)
     icxu = icx+npr*SUM0
  endif

  if (icxu.gt.MAXCX) call bummer('socfpd: (icxu-MAXCX)=',(icxu-MAXCX),2)
! print *,'SOCFPD-MAXCX :',MAXCX
END SUBROUTINE RECP_SOCFPD_SETICXU


SUBROUTINE RECP_SOCFPD_CALCX0(icu,jcu,itl,itu,jtl,jtu,npr,icx,icxa,MAXICS,MAXJCS, &
           mau,igcs,jgcs,FILE_NUM,SO_CALC,msfbct,esfb)
! # compute new cx block
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
! global variable
  INTEGER icu,jcu,itl,itu,jtl,jtu,npr,icx,icxa,MAXICS(*),MAXJCS(*),mau(*),igcs,jgcs,FILE_NUM(3)
  LOGICAL SO_CALC,msfbct,esfb
! local variable
  INTEGER ic,jc,ict,jct,it,jt,ictic,jctjc
  LOGICAL easfb,esfbc
  easfb = SO_CALC.and.esfb
  ict = 0
  DO ic = 1, icu
     ictic = ict
     IF(esfb) jcu = ic
     jct = 0
     DO jc = 1, jcu
        jctjc = jct
        esfbc = esfb .and. ic.eq.jc
        ict = ictic
        DO it = itl, itu
           ict = ict + 1
           IF(esfbc) jtu = it
           jct = jctjc
           DO jt = jtl, jtu
              jct = jct + 1
              msfbct = (.not.esfb) .or. (ic.eq.jc .and. it.eq.jt)

              CALL RECP_SOCFPD_CX0(npr,icx,icxa,MAXICS,MAXJCS,  &
                   mau,ict,jct,igcs,jgcs,FILE_NUM,easfb,msfbct)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE RECP_SOCFPD_CALCX0


SUBROUTINE RECP_SOCFPD_CX0(npr,icx,icxa,MAXICS,MAXJCS,mau,ict,jct,igcs,jgcs, &
           FILE_NUM,easfb,msfbct)
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
! global variable
  INTEGER npr,icx,icxa,MAXICS(*),MAXJCS(*),mau(*),ict,jct,igcs,jgcs,FILE_NUM(*)
  LOGICAL easfb,msfbct
! local variable
  INTEGER I,J,ICS,JCS
  REAL(8) CX1
  
  DO I = 1, npr
!    set icx/IPT_CX0 (icxa/IPT_CX0)
!    ------------------------------
     icx = icx + 1
     IPT_CX0(icx) = 0.0d0
     if (easfb) then
       icxa = icxa + 1
       IPT_CX0(icxa) = 0.0d0
     endif

!    set ICS/JCS
!    -----------
     ICS = MAXICS(I)
     JCS = MAXJCS(I)

!    write file for debug
!    --------------------
     IF (RECP_DBG.GE.1) THEN
        WRITE (FILE_NUM(2),'(A,4I10)') 'I/icx/ICS/JCS',I,icx,ICS,JCS
        IF (easfb) THEN
           WRITE (FILE_NUM(2),'(A,4I10)') 'I/icxa/ICS/JCS',I,icxa,ICS,JCS
        ENDIF
     ENDIF

!    get IPT_CX0
!    -----------
     DO J = 1, mau(I)
        ICS = ICS + 1
        JCS = JCS + 1

        CX1 = (IPT_AO2SO(ict,ICS,igcs) * IPT_AO2SO(jct,JCS,jgcs))
        IPT_CX0(icx) = IPT_CX0(icx) + CX1 
        IF (easfb) IPT_CX0(icxa) = IPT_CX0(icxa) + CX1

        IF (.NOT.msfbct) THEN
           CX1 = (IPT_AO2SO(jct,ICS,igcs) * IPT_AO2SO(ict,JCS,jgcs))
           IPT_CX0(icx) = IPT_CX0(icx) + CX1 
           IF (easfb) IPT_CX0(icxa) = IPT_CX0(icxa) - CX1
        ENDIF

!       write file for debug
        IF (RECP_DBG.GE.1) THEN
           WRITE (FILE_NUM(2),'(A,I7,A,F15.7)')'IPT_CX0(',icx,') = ',IPT_CX0(icx)
           IF (easfb) THEN
              WRITE (FILE_NUM(3),'(A,I7,A,F15.7)')'IPT_CX0(',icxa,') = ',IPT_CX0(icxa)
           ENDIF
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE RECP_SOCFPD_CX0


SUBROUTINE RECP_SOCFPD_SHIFTCX(MAXCX,MAXCX0,icx,icxu,ijsf,icxast,SO_CALC)
! move antisymmetric coefficient products down.

!       offsets: 0               MAXCX0  mblu
! before: a(*)-->................c(*)
!                0          MAXCX MAXCX0  mblu
! now:    a(*)-->cx(*)......x(*),c(*)
!                0    MAXCX2 MAXCX       (mblu+MAXCX)
! after:  a(*)-->cx(*),x(*),........
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
! global variable
  INTEGER MAXCX,MAXCX0,icx,icxu,ijsf,icxast(*)
  LOGICAL SO_CALC
! local variable
  INTEGER I,icxl
  IF (SO_CALC) THEN
     icxl  = icx + 1
     icxu  = MAXCX - icx
     MAXCX = MAXCX0 - icxu
     do icx = icxl, MAXCX
        IPT_CX0(icx) = IPT_CX0(icx+icxu)
     enddo

     do I = 1, ijsf
        icxast(I) = icxast(I) - icxu
     enddo
  ELSE
     MAXCX = icx
  ENDIF

  IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,I9)') ' * SOCFPD: MAXCX=',MAXCX
END SUBROUTINE RECP_SOCFPD_SHIFTCX


END MODULE RECP_INP_ORBITAL 
