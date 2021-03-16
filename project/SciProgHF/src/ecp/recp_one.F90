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

MODULE RECP_ONE
CONTAINS
SUBROUTINE RECP_ONE_MAIN(a,ccr,chg,eta,ica,icb,icxast,icxst,ilxyz,iprst,lcr,lls,lmnp1,    &
           lproju,mcons,mcrs,nc,ncon,ncr,nf,nfct,ngw,nkcrl,nkcru,nklsl,nklsu,nopir,npair, &
           nprir,nrcr,nt,ntl,ntu,x,y,z,zcr,zet,n_bfn_sym,INT_AREP)
  USE RECP_FUNCTION1
  USE RECP_FUNCTION2
  USE RECP_IPT
  USE RECP_ONE_SUB
  USE RECP_ONE_STVCZ
  USE RECP_OUTPUT
  USE RECP_WRITE
  implicit logical(a-z)
#include "inc_mxvalue.h"
#include "inc_print.h"
  INTEGER INT_AREP
  real*8         tol
  common /parmr/ tol
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  integer        mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,   inam,   nnam,mdum
  common /ntgr/  mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir(3),inam(5),nnam,mdum(32)
  real*8  fnfct,rr,xij,xijm,yij,yijm,zij,zijm
  integer ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,icons,igu,ircru,is,isf,itl,itu,jc,jcons,  &
          jgu,jrcru,js,jsf,jtl,jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2
  logical esf,esfc,igueq1,jgueq1
  common /one/   fnfct,rr,xij,xijm,yij,yijm,zij,zijm,ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,icons, &
                 igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,jtu, &
                 lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
! # /bufout/ holds some output integral file parameters.
  integer         itypea,itypeb,ibuf,numout,nrec,ntape
  common /bufout/ itypea,itypeb,ibuf,numout,nrec,ntape
  integer    nipv
  parameter( nipv=2 )
 
  integer ica(mcu,msu,*),icb(4,24,*),icxast(*),icxst(2,*),ilxyz(3,*),iprst(*), &
          lcr(*),lls(*),lmnp1(*),lproju,mcons(*),mcrs(*),nc(*),ncon(*),ncr(*),nf(*),nfct(*), &
          ngw(*),nkcrl(6,*),nkcru(6,*),nklsl(4,*),nklsu(4,*),nopir(*),npair(2,*), &
          nprir(2,mstu,*),nrcr(*),nt(*),ntl(*),ntu(*)
  real*8  a(*),ccr(*),chg(*),x(mcu,*),y(mcu,*),z(mcu,*),eta(mrcru,mconu,*),zcr(*),zet(mconu,*)
! local variables
  integer ixyz,ist,ifu,icu,isfis,jsfjs,ncc,igwu,SYM_I,ICA_MIN,icc,igw,jfu,if,jf,iop, &
          ipr,itst,iprm,inx,jnx,nwt,ipr1,ipr2,iju,i,j,last,iapt(3),nblu,NUC_J
  integer n_bfn_sym(mstup)  ! (nbpsy or nso)
  integer N_INTEGRAL(3), MN_INTEGRAL
  logical  es
  CHARACTER(11) FILEFORM(2)
  DATA FILEFORM /'UNFORMATTED','FORMATTED  '/

  CALL OUTPUT_LOC('RECP_ONE_MAIN','E') 
! --------------------------------------
! # operator symmetry
! --------------------------------------
  CALL RECP_SETZERO_I1(nopir,nst)
  CALL RECP_SETZERO_I2(ilxyz,3,nst)

  CALL RECP_SETZERO_I1(N_INTEGRAL,3)
  MN_INTEGRAL = 0
 
      IF (nnam.EQ.5) THEN
!        # x, y, and z components of spin-orbit integrals are all
!        # computed together.
         DO ixyz = 1, 3
            ist                   = lxyzir(ixyz)
            nopir(ist)            = nopir(ist) + 1
            ilxyz(nopir(ist),ist) = ixyz
         ENDDO
         nop     = 3
         iapt(1) = 0
         iapt(2) = 0
         iapt(3) = 0
         CALL RECP_WRITE_FILEOPEN(64,'RECP_ISA1','REPLACE')
         CALL RECP_WRITE_FILEOPEN(65,'RECP_ISA2','REPLACE')
         CALL RECP_WRITE_FILEOPEN(66,'RECP_ISA3','REPLACE')
         CALL RECP_WRITE_FILEOPEN(67,'RECP_ISB1','REPLACE')
         CALL RECP_WRITE_FILEOPEN(68,'RECP_ISB2','REPLACE')
         CALL RECP_WRITE_FILEOPEN(69,'RECP_ISB3','REPLACE')
      ELSE
         nopir(1)   = 1
         ilxyz(1,1) = 1
         nop        = 1
         iapt(1)    = 0
         CALL RECP_WRITE_FILEOPEN(64,'RECP_ISA0','REPLACE')
         CALL RECP_WRITE_FILEOPEN(67,'RECP_ISB0','REPLACE')
      ENDIF
 
      ijsf = 0
      isf  = 0

      do 712 is = 1,ns
      ifu = nf(is)
      if(ifu.eq.0) go to 712
      icu=nc(is)
      isfis=isf
      jsf = 0
      do 708 js = 1,is
      jsfjs=jsf
      es = is.eq.js
      ic=icu
      ncc = 0

!     --------------------------------------
!     center combinations for is=js or is>js
!     --------------------------------------
      IF (es) THEN 
!        center combinations for is=js
         CALL RECP_1ESUB_ISEQJS(is,js,ic,ng,ncc,mcu,msu,ica,icb,nfct,ngw)
      ELSE
!        center combinations for is>js
         jfu = nf(js) !YCP FIX
         if(nf(js).eq.0) goto 708
         CALL RECP_1ESUB_ISGTJS(is,js,nc,ng,ic,mcu,msu,ica,icb,ncc,nfct)
      ENDIF
      IF (ncc.GT.mccu) call bummer('change mccup (one place) to ',ncc,2)

!     --------------------------------------
!     evaluate integrals
!     --------------------------------------
      IF (RECP_DBG.GE.2) WRITE(RECP_OUT,*) "evaluate integrals" 

      isf=isfis
      do 704 if = 1, ifu

      isf = isf + 1
      CALL RECP_1ESUB_IBLOCK(icons,ircru,lit,igu,isf,mcons,nrcr,lmnp1,ncon)
      itl=ntl(isf)
      itu=ntu(isf)

      jsf=jsfjs
      if(es) jfu=if
      do 700 jf = 1, jfu

      jsf = jsf + 1
      CALL RECP_1ESUB_JBLOCK(jcons,jrcru,ljt,jsf,mcons,nrcr,lmnp1)

      ijsf = ijsf + 1
      ibl1=0
      ibl2=0
      esf=isf.eq.jsf
      do ist=1,nst
        iop=nopir(ist)
        if(iop.ne.0) then
          if(esf) ibl1=ibl1+iop*nprir(1,ist,ijsf)
          ibl2=ibl2+iop*nprir(2,ist,ijsf)
        endif
      enddo
      iju=nt(isf)*nt(jsf)
      ntij2=iju*npair(2,ijsf)
      icxi2=icxst(2,ijsf)
!     print *,'isf/jsf',is,js,if,jf
      if(esf) then
        nblt1=npair(1,ijsf)
        ntij1=iju*nblt1
        nblu =ircru*ibl1 + ( ircru*(ircru-1)/2 )*ibl2
        if(nnam.eq.5) then
          icxi1=icxast(ijsf)
        else
          icxi1=icxst(1,ijsf)
        endif
      else
        nblu=ircru*jrcru*ibl2
      endif

!     -----------------------
      if(nblu.eq.0) go to 696

!     # allocate space.
      CALL RECP_IPTA_ONEINT(nblu)
      do i = 1, nblu
         IPT_ONE_H2(i) = 0.0d0
      enddo

      ipr=iprst(ijsf)
 
      jgu=ncon(jcons)
      jtl=ntl(jsf)
      jtu=ntu(jsf)

      if ((.not.es).or.esf) then
        do icc = 1, ncc
           ic = icb(1,1,icc)
           jc = icb(2,1,icc)
           fnfct = nfct(icc)
           CALL RECP_ONE_STVCZ0(a,ccr,chg,eta,ilxyz,lcr,lls,IPT_LMNV,lproju,mcrs,nc,ncr, &
                nkcrl,nkcru,nklsl,nklsu,nopir,nprir,nt,x,y,z,zcr,zet )
        enddo
      else
        do icc = 1, ncc
          igwu = ngw(icc)
          itst = 0
          do 520 iprm = 1, 2
            if (iprm.eq.1) then
              inx=1
              jnx=2
            else
              inx=2
              jnx=1
            endif
            jc = icb(jnx,1,icc)
            if (jc.eq.itst) goto 520
            ic = icb(inx,1,icc)
            nwt = 1
            do igw = 2, igwu
              if(icb(jnx,igw,icc).lt.jc) goto 520
              nwt = nwt + 1
            enddo
            fnfct = (nwt * nfct(icc))
            CALL RECP_ONE_STVCZ0( a,ccr,chg,eta,ilxyz,lcr,lls,IPT_LMNV,  &
                 lproju,mcrs,nc,ncr,nkcrl,nkcru,nklsl,nklsu,nopir,nprir,nt,x,y,z,zcr,zet )
            itst = jc
520       continue
        enddo
      endif
!     # put integrals and labels into the intermediate storage arrays.
      CALL RECP_WRITE_ISA(ircru,jrcru,nst,mstu,nprir,ijsf,  &
              nopir,ilxyz,iapt,mpru,ipr,N_INTEGRAL,IPT_ONE_H2,esf)
!     Deallocate IPT_ONE_H2
      CALL RECP_IPTD_ONEINT
696   continue
700   continue
704   continue
708   continue
712   continue

      IF (nnam.EQ.5) THEN
         IF (RECP_DBG.GE.1) WRITE(RECP_OUT,'(A)')'*End writing ISA(1/2/3),ISB(1/2/3)' 
         CLOSE (UNIT=64)
         CLOSE (UNIT=65)
         CLOSE (UNIT=66)
         CLOSE (UNIT=67)
         CLOSE (UNIT=68)
         CLOSE (UNIT=69)
      ELSE
         IF (RECP_DBG.GE.1) WRITE(RECP_OUT,'(A)')'*End writing ISA0,ISB0' 
         CLOSE (UNIT=64)
         CLOSE (UNIT=67)
      ENDIF
!     # write 1-e integral arrays.
      MN_INTEGRAL = N_INTEGRAL(1) + N_INTEGRAL(2) + N_INTEGRAL(3)
      CALL RECP_WRITE_MAIN(nop,iapt,ipr,ipr1,ipr2,n_bfn_sym,N_INTEGRAL,  &
           MN_INTEGRAL,INT_AREP,nnam,mpru,nbft,nnbft,nst)
      CALL OUTPUT_LOC('RECP_ONE_MAIN','X') 
END SUBROUTINE RECP_ONE_MAIN
END MODULE RECP_ONE
