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

MODULE RECP_CSO_SUB
CONTAINS

!------------------------------------------------------------------------------
SUBROUTINE RECP_CSO_INIT1(iu,lit,mrcru,fctr1,etai,etaj,esfc,igueq1,jgueq1)
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER iu,lit,mrcru
  REAL(8) fctr1,etai(mrcru,*),etaj(mrcru,*)
  LOGICAL esfc,igueq1,jgueq1
  INTEGER i
! set CRDA/CRDB
  DO i=1,3
     IPT_CZ_CRDA(1,i)=1.0D0
     IPT_CZ_CRDB(1,i)=1.0D0
  ENDDO

! set iu
  IF (esfc) iu=(lit+1)*lit/2

! set fctr1
  IF (igueq1) THEN
     fctr1=4.0D0*etai(1,1)
  ELSE
     fctr1=4.0D0
  ENDIF
  IF (jgueq1) THEN
     fctr1=fctr1*etaj(1,1)
  ENDIF
END SUBROUTINE RECP_CSO_INIT1

!------------------------------------------------------------------------------

SUBROUTINE RECP_CSO_G1ZERO(nc1,tmp_ic0,igueq1)
! Set G1(and GOUT) zero in i primitive
! Called by cints,lsints
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER nc1,tmp_ic0
  LOGICAL igueq1
  INTEGER i  !local variable

! size of variables  
! +--G1---+-----GOUT-------+
! i       tmp_ic0    nc1   (imp_ic0+ij)

  IF (.NOT.igueq1) THEN
     DO i = 1, tmp_ic0  
        IPT_CZ_G1(i) = 0.0D0
     ENDDO
     IF (nc1.GT.tmp_ic0) THEN
        DO i = 1,(nc1-tmp_ic0)
            IPT_CZ_GOUT(i) = 0.0D0
        ENDDO
     ENDIF
  ENDIF
END SUBROUTINE RECP_CSO_G1ZERO

!------------------------------------------------------------------------------

SUBROUTINE RECP_CSO_GOUTZERO(ij,jgueq1)
! Set GOUT zero in j primitive
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER ij
  LOGICAL jgueq1
  INTEGER i  !local variable

  IF (.NOT.jgueq1) THEN
     DO i = 1, ij
        IPT_CZ_GOUT(i) = 0.0d0
     ENDDO
  ENDIF
END SUBROUTINE RECP_CSO_GOUTZERO

!------------------------------------------------------------------------------

SUBROUTINE RECP_CSO_SETCRDA(xka,yka,zka,xc,yc,zc,xi,yi,zi,ca,lit)
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  REAL(8) xka,yka,zka,xc,yc,zc,xi,yi,zi,ca
  INTEGER lit
  INTEGER i,j  !local variables
  xka=xc-xi
  yka=yc-yi
  zka=zc-zi
  ca=SQRT(xka*xka+yka*yka+zka*zka)
  IF (lit.GE.2) THEN
     IPT_CZ_CRDA(2,1)=xka
     IPT_CZ_CRDA(2,2)=yka
     IPT_CZ_CRDA(2,3)=zka
     IF (lit.GE.3) THEN
        DO i=1,3
           DO j=3,lit
              IPT_CZ_CRDA(j,i)=IPT_CZ_CRDA(2,i)*IPT_CZ_CRDA(j-1,i)
           ENDDO
        ENDDO
     ENDIF
  ENDIF
END SUBROUTINE RECP_CSO_SETCRDA

!------------------------------------------------------------------------------

SUBROUTINE RECP_CSO_SETCRDB(xkb,ykb,zkb,xc,yc,zc,xj,yj,zj,cb,ljt)
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  REAL(8) xkb,ykb,zkb,xc,yc,zc,xj,yj,zj,cb
  INTEGER ljt
  INTEGER i,j  !local variables 
  xkb=xc-xj
  ykb=yc-yj
  zkb=zc-zj
  cb=SQRT(xkb*xkb+ykb*ykb+zkb*zkb)
  IF (ljt.GE.2) THEN
     IPT_CZ_CRDB(2,1)=xkb
     IPT_CZ_CRDB(2,2)=ykb
     IPT_CZ_CRDB(2,3)=zkb
     IF (ljt.GE.3) THEN
        DO i=1,3
           DO j=3,ljt
              IPT_CZ_CRDB(j,i)=IPT_CZ_CRDB(2,i)*IPT_CZ_CRDB(j-1,i)
           ENDDO
        ENDDO
     ENDIF
  ENDIF
END SUBROUTINE RECP_CSO_SETCRDB

!------------------------------------------------------------------------------

SUBROUTINE RECP_CSO_G2CAL1(ircr,ig,iu,i1,j1,ij2,tmp_ic0,tmp_ic1,mrcru,g2,etai)
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER ircr,ig,iu,i1,j1,ij2,tmp_ic0,tmp_ic1,mrcru
  REAL(8) g2(*),etai(mrcru,*)
! DIMENSION g2(*),etai(mrcru,*)
! local variables
  INTEGER jrcr,i,j,i1n
  REAL(8) G1MOD1,G1MOD2

  DO jrcr=1,ircr
     DO i=1,iu
        i1n=i1+i
        DO j=1,iu
!          * Check G1 index
           IF ((j1+j).GT.tmp_ic0) THEN
              IF(RECP_DBG.GE.5) WRITE(RECP_OUT,*)'G1 index exceed(i11)',(j1+j),'>',tmp_ic0
              G1MOD1 = IPT_CZ_GOUT((j1+j)-tmp_ic0)
           ELSE
              G1MOD1 = IPT_CZ_G1(j1+j)
           ENDIF

!          * Check G1 index
           IF (i1n.GT.tmp_ic0) THEN
              IF(RECP_DBG.GE.5) WRITE(RECP_OUT,*)'G1 index exceed(i12)',i1n,'>',tmp_ic0
              G1MOD2 = IPT_CZ_GOUT(i1n-tmp_ic0)
           ELSE
              G1MOD2 = IPT_CZ_G1(i1n)
           ENDIF

!          * Check G2 index
           IF ((ij2+j).GT.tmp_ic1) THEN
              IF(RECP_DBG.GE.5) WRITE(RECP_OUT,*)'G2 index exceed(i11)',(ij2+j),'>',tmp_ic1
!             Error termination
!             CALL RECP_EXIT
           ENDIF

!          * Get G2 
           g2(ij2+j) = g2(ij2+j) + G1MOD1*etai(ircr,ig) + G1MOD2*etai(jrcr,ig)
!          IF(RECP_DBG.GE.5) WRITE(RECP_OUT,*)'g1/g1',(j1+j),G1MOD1,(i1n),G1MOD2

           i1n=i1n+iu
        ENDDO
        ij2=ij2+iu
        j1=j1+iu
     ENDDO   
  ENDDO   
END SUBROUTINE RECP_CSO_G2CAL1

!------------------------------------------------------------------------------

SUBROUTINE RECP_CSO_G2CAL2(ircr,ig,i1,j1,iu,tmp_ic0,tmp_ic1,ij2,mrcru,g2,etai)
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER ircr,ig,i1,j1,iu,tmp_ic0,tmp_ic1,ij2,mrcru
  REAL(8) g2(*),etai(mrcru,*)
! local variables
  INTEGER jrcr,i1m,i1n,i,j,k
  REAL(8) G1MOD1,G1MOD2
  DO jrcr=1,ircr
     i1m=i1
     DO i=1,iu
        i1n=i1m
        DO j=1,iu
           DO k=1,3
!             * Check G1 index
              IF ((j1+k).GT.tmp_ic0) THEN
                 IF(RECP_DBG.GE.5) WRITE(RECP_OUT,*)'G1 index exceed(i11)',(j1+k),'>',tmp_ic0
                 G1MOD1 = IPT_CZ_GOUT((j1+k)-tmp_ic0)
              ELSE
                 G1MOD1 = IPT_CZ_G1(j1+k)
              ENDIF

!             * Check G1 index
              IF ((i1n+k).GT.tmp_ic0) THEN
                 IF(RECP_DBG.GE.5) WRITE(RECP_OUT,*)'G1 index exceed(i12)',(i1n+k),'>',tmp_ic0
                 G1MOD2 = IPT_CZ_GOUT((i1n+k)-tmp_ic0)
              ELSE
                 G1MOD2 = IPT_CZ_G1(i1n+k)
              ENDIF

!             * Check G2 index
              IF ((ij2+k).GT.tmp_ic1) THEN
                 IF(RECP_DBG.GE.5) WRITE(RECP_OUT,*)'G2 index exceed(i11)',(ij2+k),'>',tmp_ic1
!                Error termination
!                CALL RECP_EXIT
              ENDIF

!             * Get G2              
              g2(ij2+k)=g2(ij2+k) + G1MOD1*etai(ircr,ig) - G1MOD2*etai(jrcr,ig)
!             IF(RECP_DBG.GE.5) WRITE(RECP_OUT,*)'g1/g1',(j1+k),G1MOD1,(i1n+k),G1MOD2
           ENDDO
           ij2=ij2+3
           j1=j1+3
           i1n=i1n+(3*iu)
        ENDDO
        i1m=i1m+3
     ENDDO
  ENDDO
END SUBROUTINE RECP_CSO_G2CAL2

!------------------------------------------------------------------------------

SUBROUTINE RECP_CSO_G2CAL3(ircr,jrcru,ig,ij,ij2,j1,tmp_ic0,tmp_ic1,mrcru,g2,etai)
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER ircr,jrcru,ig,ij,ij2,j1,tmp_ic0,tmp_ic1,mrcru
  REAL(8) g2(*),etai(mrcru,*)
! DIMENSION g2(*),etai(mrcru,*)
! local variables
  INTEGER jrcr,i
  REAL(8) G1MOD
  DO jrcr=1,jrcru
     DO i=1,ij
!       * Check G1 index
        IF ((j1+i).GT.tmp_ic0) THEN
           IF(RECP_DBG.GE.5) WRITE(RECP_OUT,*)'G1 index exceed(i21)',(j1+i),'>',tmp_ic0
           G1MOD = IPT_CZ_GOUT((j1+i)-tmp_ic0)
        ELSE
           G1MOD = IPT_CZ_G1(j1+i)
        ENDIF

!       * Check G2 index
        IF ((ij2+i).GT.tmp_ic1) THEN
           IF(RECP_DBG.GE.5) WRITE(RECP_OUT,*)'G2 index exceed(i21)',(ij2+i),'>',tmp_ic1
        ENDIF

        g2(ij2+i)=g2(ij2+i)+G1MOD*etai(ircr,ig)
!       WRITE(RECP_OUT,*)'2* g2/g1',(ij2+i),g2(ij2+i),(j1+i),G1MOD
     ENDDO
     ij2=ij2+ij
     j1=j1+ij
  ENDDO
END SUBROUTINE RECP_CSO_G2CAL3

!------------------------------------------------------------------------------

END MODULE RECP_CSO_SUB
