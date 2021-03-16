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

MODULE RECP_SOCFPD

! local
  INTEGER, ALLOCATABLE :: icxsv1(:)
  INTEGER, ALLOCATABLE :: icxsv2(:,:)
CONTAINS

SUBROUTINE RECP_SOCFPD_MAIN(icxast,icxst,idp,iprst,la,   &
           lb,maords,mau,mcons,mgcs,nblpr,nc,nd,nf,nir,npair,nprir, &  
           nrcr,nsopr,nt,ntl,ntu,MAXCX)
! Compute symmetry orbital coefficient products.
! on return from socfpd(*), cx(1:MAXCX) was computed.
  USE RECP_FUNCTION1
  USE RECP_FUNCTION2
  USE RECP_INP_ORBITAL
  USE RECP_IPT
  USE RECP_OUTPUT
  USE RECP_WRITE
  implicit logical(a-z)
! il(*) packing factor. see also oneint() and wtint2().
  integer    ilfact
  parameter( ilfact=1024 )
#include "inc_mxvalue.h"
#include "inc_print.h"
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  integer        mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,   inam,   nnam,mdum
  common /ntgr/  mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir(3),inam(5),nnam,mdum(32)
  integer MAXCX

  integer icxast(*),icxst(2,*),idp(mstu,mstu,*),  &
          iprst(*),la(mru,*),lb(*),maords(*),mau(*),mcons(*),  &
          mgcs(*),nblpr(*),nc(*),nd(*),nf(*),nir(*),  &
          npair(2,*),nrcr(*),nsopr(*),nprir(2,mstu,*),nt(*),ntl(*),ntu(*)
! local variables
  integer MAXICS(mnrup),MAXJCS(mnrup)
  integer icx, ipair, ijsf, isf, isfr, igcs, jgcs, is, &
          icu, isf1, isfr1, jsf, jsfr, js, jsf1, jsfr1, jfu, &
          jcu, if, isfr2, iaords, iru, ircru, itl, itu, jf, jaords, &
          jtl, jsfr2, jtu, jru, ist, jrcru, ircr, isfrib, jrcr, &
          jsfrjb,iesfb,npr,icxu,icxa,ict,ic,jct,jc,it, &
          jt,I,J,K,FILE_NUM(3),RETURNVALUE
  INTEGER MAXCX0
  LOGICAL SO_CALC,esfb,msfbct
  CHARACTER(10) FILE_NAME(3)
  DATA FILE_NUM  / 64, 65, 66 /
  DATA FILE_NAME / 'RECP_SOC_0','RECP_SOC_1','RECP_SOC_2' /

  CALL OUTPUT_LOC('RECP_SOCFPD_MAIN','E')

! write file in debug mode
  IF (RECP_DBG.GE.1) THEN
     CALL RECP_WRITE_FILEOPEN(FILE_NUM(2),FILE_NAME(2),'REPLACE')
     CALL RECP_WRITE_FILEOPEN(FILE_NUM(3),FILE_NAME(3),'REPLACE')
  ENDIF

  CALL RECP_SETZERO_I3(nprir,2,mstup,mnsfup)

! # check the il(*) packing factor.
  if (nbft.ge.ilfact) call bummer('socfpd: larger ilfact required, nbft=',nbft,2)
 
  SO_CALC = inam(nu).eq.5
  MAXCX0= MAXCX   !original MAXCX offset
  icx   = 0
  ipair = 0
  ijsf = 0  !ijsf = one dimensionalization of is, js, if, jf

  ALLOCATE(icxsv1(ngcs))
  ALLOCATE(icxsv2(ngcs,ngcs))
  icxsv1=0
  icxsv2=0

  isf  = 0
  isfr = 0  !isfr = one dimensionalization of is, if, ir
  DO is = 1, ns      ! icenter
     if(nf(is).eq.0) goto 340
     isf1 = isf
     isfr1 = isfr

     jsf = 0
     jsfr = 0
     DO js = 1,is     ! jcenter
        jsf1 = jsf
        jsfr1 = jsfr

        if(nf(js).eq.0) goto 330

        isf = isf1
        isfr = isfr1
        DO if = 1,nf(is)   ! i-function
           isf = isf+1
           isfr2 = isfr

           jsf = jsf1
           jsfr = jsfr1

           jfu = nf(js)
           if(is.eq.js) jfu = if
           DO jf = 1,jfu ! j-function
              jsf = jsf+1
              jsfr2 = jsfr

              igcs = mgcs(isf)
              jgcs = mgcs(jsf)
              iaords = maords(igcs)
              jaords = maords(jgcs)

              iru = nir(iaords)
              ircru = nrcr(mcons(isf))
              jtl = ntl(jsf)
              ijsf = ijsf+1
              if (isf.eq.jsf) then
                 if(ircru.ge.2) then
                   jtu = ntu(jsf)
                   jru = nir(jaords)
                 else
                   npair(2,ijsf) = 0
                   icxst(2,ijsf) = 0
                   do ist = 1,nst
                     nprir(2,ist,ijsf) = 0
                   enddo
                 endif
              else
                 jrcru = nrcr(mcons(jsf))
                 jtu   = ntu(jsf)
                 jru   = nir(jaords)
              endif
              iprst(ijsf) = ipair
              isfr = isfr2
              DO ircr = 1,ircru
                 isfrib = isfr
                 jsfr = jsfr2
                 if(isf.eq.jsf) jrcru = ircr
                 DO jrcr = 1,jrcru
                    jsfrjb = jsfr
                    esfb = ((isf.eq.jsf).and.(ircr.eq.jrcr))

                    CALL RECP_SOCFPD_SETNPR(nst,npr,isfr,jsfr,isfrib, &
                         jsfrjb,iru,jru,iaords,jaords,ipair,ijsf,mru,la, &
                         lb,mstu,idp,nsopr,mau,MAXICS,MAXJCS,nd,nprir,npair,ilfact,esfb)
                  
                    IF (npr.NE.0) THEN 
                       icu = nc(is)
                       jcu = nc(js)
                       CALL RECP_SOCFPD_SETICXU(RETURNVALUE,igcs,jgcs,ijsf,icx, &
                       icu,icxu,isf,jsf,npr,MAXCX,icxa,ngcs, &
                       icxsv1,icxsv2,icxst,icxast,nt,jcu,esfb,SO_CALC)
                       IF (RETURNVALUE.EQ.0) THEN
!                         compute new cx block
                          itl = ntl(isf)
                          itu = ntu(isf)
                          CALL RECP_SOCFPD_CALCX0(icu,jcu,itl,itu,jtl,jtu,npr,icx, &
                               icxa,MAXICS,MAXJCS,mau,igcs,jgcs,FILE_NUM,SO_CALC,msfbct,esfb)
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO     
        330 CONTINUE
     ENDDO
     340 CONTINUE
  ENDDO 

  DEALLOCATE(icxsv1)
  DEALLOCATE(icxsv2)

  CALL RECP_SOCFPD_SHIFTCX(MAXCX,MAXCX0,icx,icxu,ijsf,icxast,SO_CALC) 

! Allocate variables cx(*)
  ALLOCATE (IPT_CX(MAXCX))
  DO I = 1, MAXCX
     IPT_CX(I) = IPT_CX0(I)
  ENDDO
! Deallocate variables cx_temp(*),c(*)
  DEALLOCATE(IPT_CX0)
  DEALLOCATE(IPT_AO2SO)

  IF (RECP_DBG.GE.1) THEN
     CLOSE(FILE_NUM(2))
     CLOSE(FILE_NUM(3))
  ENDIF
  CALL OUTPUT_LOC('RECP_SOCFPD_MAIN','X')
END SUBROUTINE RECP_SOCFPD_MAIN


END MODULE RECP_SOCFPD
