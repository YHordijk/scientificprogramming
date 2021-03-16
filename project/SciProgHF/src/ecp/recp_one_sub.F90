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

MODULE RECP_ONE_SUB
CONTAINS

SUBROUTINE RECP_1ESUB_ISEQJS(is,js,ic,ng,ncc,mcu,msu,ica,icb,nfct,ngw)
! center combinations for is=js, Y.C.Park
! ica(NUC_I/J, JS/IS, SYM/SYMUNQ)
! icb(NUC_I/J, JS/IS, SYM/SYMUNQ)
!        ^       ^       ^
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "inc_print.h"
! vaiables : global
  INTEGER :: is             ! DO is->ns
  INTEGER :: js             ! DO js->is
  INTEGER :: ic             ! CENTER_I nc(is)
  INTEGER :: ng             ! SYM_N : symmetry_number from ngen?
  INTEGER :: ncc            ! number of symmetry unique center?
  INTEGER :: mcu, msu       ! from common(parmi)
  INTEGER :: ica(mcu,msu,*)
  INTEGER :: icb(4,24,*)
  INTEGER :: nfct(*)
  INTEGER :: ngw(*)
! vaiables : local
  INTEGER :: NUC_J          ! NUC_CENTER_J (jc)
  INTEGER :: SYM_I          ! SYM_I : symmetry_ith (ig)
  INTEGER :: ICA_MIN
  INTEGER :: igw, igwu
  INTEGER :: NUC_SYMUNQ_I   ! NUC_CENTER_SYM_UNIQUE_I (icc)
  INTEGER :: TEMP

      CALL OUTPUT_LOC('ONEINT_CENTER_ISEQJS','E')
      do 140 NUC_J = 1, ic
      igwu=0
      do 132 SYM_I = 1, ng
        if(ica(ic,is,SYM_I)   .ne.ic .and. &
           ica(NUC_J,is,SYM_I).ne.ic) goto 132
        ICA_MIN=min(ica(ic,is,SYM_I),ica(NUC_J,is,SYM_I))
 
        do 112 NUC_SYMUNQ_I=ncc,1,-1
          if(ICA_MIN-icb(2,1,NUC_SYMUNQ_I)) 112,136,116
  112   continue
 
  116   if(ICA_MIN.ne.NUC_J) go to 132
 
        do igw=igwu,1,-1
          if(ica(ic,is,SYM_I)   .eq.icb(1,igw,ncc+1).and. &
             ica(NUC_J,is,SYM_I).eq.icb(2,igw,ncc+1)) goto 132
        enddo
 
        igwu=igwu+1
        icb(1,igwu,ncc+1)=ica(ic,is,SYM_I)
        icb(2,igwu,ncc+1)=ica(NUC_J,is,SYM_I)
  132 continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc)=igwu
      go to 140
  136 nfct(NUC_SYMUNQ_I) = nfct(NUC_SYMUNQ_I) + ic
  140 continue

      do NUC_SYMUNQ_I = 1, ncc-1
        nfct(NUC_SYMUNQ_I) = nfct(NUC_SYMUNQ_I)/2
      enddo

      CALL OUTPUT_LOC('ONEINT_CENTER_ISEQJS','X')
END SUBROUTINE RECP_1ESUB_ISEQJS


SUBROUTINE RECP_1ESUB_ISGTJS(is,js,nc,ng,ic,mcu,msu,ica,icb,ncc,nfct)
! Center combinations for is>js
! Y. C. Park
!
! ica(CENTER_I/J, JS/IS, SYM/SYMUNQ)
! icb(CENTER_I/J, JS/IS, SYM/SYMUNQ)
!        ^          ^       ^
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "inc_print.h"
! vaiables : global
  INTEGER :: is
  INTEGER :: js
  INTEGER :: nc(*)
  INTEGER :: ng             ! SYM_N : symmetry_number from ngen?
  INTEGER :: ic
  INTEGER :: mcu, msu       ! from common(parmi)
  INTEGER :: ica(mcu,msu,*)
  INTEGER :: icb(4,24,*)
  INTEGER :: ncc            ! number of symmetry unique center?
  INTEGER :: nfct(*)
! vaiables : local
  INTEGER :: NUC_J         ! CENTER_J (jc)
  INTEGER :: SYM_I          ! SYM_I : symmetry_ith (ig)
  INTEGER :: NUC_SYMUNQ_I  ! CENT_SYM_UNIQUE_I (icc)

      CALL OUTPUT_LOC('ONEINT_CENTER_ISGTJS','E')

      do 240 NUC_J = 1, nc(js)
      do 232 SYM_I = 1, ng
         if(ica(ic,is,SYM_I) .eq. ic) then
           do 212 NUC_SYMUNQ_I=ncc,1,-1
              if(ica(NUC_J,js,SYM_I) - icb(2,1,NUC_SYMUNQ_I)) 212,236,232
  212      continue
         endif
  232 continue
      ncc = ncc + 1
      icb(1,1,ncc) = ic
      icb(2,1,ncc) = NUC_J
      nfct(ncc)    = ic
      goto 240
  236 nfct(NUC_SYMUNQ_I) = nfct(NUC_SYMUNQ_I) + ic
  240 continue

      CALL OUTPUT_LOC('ONEINT_CENTER_ISGTJS','X')
END SUBROUTINE RECP_1ESUB_ISGTJS


SUBROUTINE RECP_1ESUB_IBLOCK(icons,ircru,lit,igu,isf,mcons,nrcr,lmnp1,ncon)
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER icons,ircru,lit,igu,isf,mcons(*),nrcr(*),lmnp1(*),ncon(*)

  icons=mcons(isf)  !mcons : mcons' basis block
  ircru=nrcr(icons) !nrcr  : number of contraction coefficient
  lit=lmnp1(icons)  !lmnp1 : angular momentum value
  igu=ncon(icons)   !ncon  : number of primitive in a contraction
! print
  IF (RECP_DBG.GE.5) THEN
     WRITE(RECP_OUT,*) ' ------------------------------------------- '
     WRITE(RECP_OUT,*) ' * basis block (i)          : ',icons
     WRITE(RECP_OUT,*) ' * num. cont. coeff. (in i) : ',ircru
     WRITE(RECP_OUT,*) ' * num. ang momentum (in i) : ',lit
     WRITE(RECP_OUT,*) ' * num. prim.        (in i) : ',igu
     WRITE(RECP_OUT,*) ' ------------------------------------------- '
  ENDIF
END SUBROUTINE RECP_1ESUB_IBLOCK

!------------------------------------------------------------------------------

SUBROUTINE RECP_1ESUB_JBLOCK(jcons,jrcru,ljt,jsf,mcons,nrcr,lmnp1)
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER jcons,jrcru,ljt,jsf,mcons(*),nrcr(*),lmnp1(*)

  jcons=mcons(jsf)
  jrcru=nrcr(jcons)
  ljt=lmnp1(jcons)
! print
  IF (RECP_DBG.GE.5) THEN
     WRITE(RECP_OUT,*) ' ------------------------------------------- '
     WRITE(RECP_OUT,*) ' * basis block (j)          : ',jcons
     WRITE(RECP_OUT,*) ' * num. cont. coeff. (in j) : ',jrcru
     WRITE(RECP_OUT,*) ' * num. ang momentum (in j) : ',ljt
     WRITE(RECP_OUT,*) ' ------------------------------------------- '
  ENDIF
END SUBROUTINE RECP_1ESUB_JBLOCK


!------------------------------------------------------------------------------
!         RECP_STVCZ
!------------------------------------------------------------------------------


SUBROUTINE RECP_1ESUB_NC12(nc1,nc2,ircru,jrcru,ij,esfc)
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER nc1,nc2,ircru,jrcru,ij
  LOGICAL esfc 
  INTEGER TEMP !local

  nc1=jrcru*ij
  IF (esfc) THEN
    TEMP = (ircru+1)*ircru/2
    nc2  = TEMP*ij
  ELSE
    nc2  = ircru*nc1
  ENDIF
END SUBROUTINE RECP_1ESUB_NC12

!------------------------------------------------------------------------------

SUBROUTINE RECP_1ESUB_IC12(ic0,ic1,nc1,nc2,igueq1,jgueq1)
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER ic0,ic1,nc1,nc2
  LOGICAL igueq1,jgueq1

  IF(igueq1) THEN
    ic1=0
  ELSE
    ic1=nc2
  ENDIF

  IF(jgueq1) THEN
    ic0=0
  ELSE
    ic0=nc1
  ENDIF
END SUBROUTINE RECP_1ESUB_IC12

!------------------------------------------------------------------------------

SUBROUTINE RECP_1ESUB_H2INDEX(icx,icx1,icx2,icx3,ircr,jrcr,ndx,ndxa,ndxi, &
           ibl2,ibld,ngti,ngtj,nop,isf,jsf,nt,esfb,esfc,esf)
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER icx,icx1,icx2,icx3,ircr,jrcr,ndx,ndxa,ndxi,ibl2,ibld,ngti,ngtj,nop
  INTEGER isf,jsf,nt(*)
  LOGICAL esfb,esfc,esf
! local variable
  INTEGER IPQ_TEMP

! WRITE(RECP_OUT,*)'H2INDEX:icx',icx1,icx2,icx3
  IF (esfb) THEN
     icx=icx1
  ELSE
     IF ((esf).AND.(.NOT.esfc)) THEN
        IF (ircr.LE.jrcr) THEN
           IPQ_TEMP = (jrcr)*((jrcr)-1)/2
           ndxa=(IPQ_TEMP+(ircr-1))*ibl2-(jrcr-1)*ibld
           icx=icx3
           ngti=nop           !ngti/ngtj/ndxi
           ngtj=(nop*nt(isf))
           ndxi=ndx
           RETURN   !return
        ENDIF

        IF (jrcr.EQ.1) THEN
           IPQ_TEMP = (ircr)*((ircr)-1)/2
           ndxa=(IPQ_TEMP+(jrcr-1))*ibl2-(ircr-1)*ibld
        ENDIF
     ENDIF
     icx=icx2
  ENDIF

  ngti=(nop*nt(jsf))
  ngtj=nop
  ndxi=ndx
END SUBROUTINE RECP_1ESUB_H2INDEX


SUBROUTINE RECP_1ESUB_H2WRITE(ic1,ic0,INDEX1,ilxyz,nst,ijsf,ij,ndxb,icx,nopir, &
           mstu,nprir,chsign,esfb,fnfct)
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER ic1,ic0,INDEX1,ilxyz(3,*),nst,ijsf,ij,ndxb,icx
  INTEGER nopir(*),mstu,nprir(2,mstu,*)
  LOGICAL chsign,esfb
  REAL(8) fnfct
! local variable
  INTEGER ist,iop,npr,J 
  REAL(8) VAL
  CHARACTER(5) STR0

  DO ist=1,nst     !nst:symmetry
     IF (esfb) THEN
        npr=nprir(1,ist,ijsf)
     ELSE
        npr=nprir(2,ist,ijsf)
     ENDIF

!    IF ((nopir(ist).ne.0) .and. (npr.ne.0)) THEN
        DO iop=1,nopir(ist)
           IF (ic1.NE.0) THEN
              STR0= 'G2  :' 
              IF ((INDEX1+ilxyz(iop,ist)).GT.ic1) THEN
                 WRITE(6,*) 'IPT_CZ_G2 INDEX1 exceeded maximum'
                 WRITE(6,*) (INDEX1+ilxyz(iop,ist)),'>',ic1,'(ic1)'
                 CALL QUIT('RECP_1ESUB_H2WRITE: IPT_CZ_G2 error')
              ENDIF
              VAL = IPT_CZ_G2(INDEX1+ilxyz(iop,ist))
           ELSEIF (ic0.NE.0) THEN
              STR0= 'G1  :' 
              IF ((INDEX1+ilxyz(iop,ist)).GT.ic0) THEN
                 WRITE(6,*) 'IPT_CZ_G1 INDEX1 exceeded maximum'
                 WRITE(6,*) (INDEX1+ilxyz(iop,ist)),'>',ic0,'(ic0)'
                 CALL QUIT('RECP_1ESUB_H2WRITE: IPT_CZ_G1 error')
              ENDIF
              VAL = IPT_CZ_G1(INDEX1+ilxyz(iop,ist))
           ELSE
              STR0= 'GOUT:' 
              IF ((INDEX1+ilxyz(iop,ist)).GT.ij) THEN
                 WRITE(6,*) 'IPT_CZ_GOUT INDEX1 exceeded maximum'
                 WRITE(6,*) (INDEX1+ilxyz(iop,ist)),'>',ij,'(ij)'
                 CALL QUIT('RECP_1ESUB_H2WRITE: IPT_CZ_GOUT error')
              ENDIF
              VAL = IPT_CZ_GOUT(INDEX1+ilxyz(iop,ist))
           ENDIF
           if(chsign) VAL=-VAL

!          * print values
           IF (RECP_DBG.GE.1) WRITE(6,'(A,L1,F15.7,2(A,I2))') & 
              STR0,chsign,VAL,' iop=',iop,' ist=',ist
           VAL=fnfct*VAL

           DO J=1,npr
              IPT_ONE_H2(ndxb+J) = IPT_ONE_H2(ndxb+J) + IPT_CX(icx+J)*VAL
           ENDDO
           ndxb=ndxb+npr
        ENDDO
!    ENDIF
     icx=icx+npr
  ENDDO
END SUBROUTINE RECP_1ESUB_H2WRITE

END MODULE RECP_ONE_SUB
