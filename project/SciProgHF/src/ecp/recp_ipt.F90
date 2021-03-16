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

MODULE RECP_IPT
  REAL(8), ALLOCATABLE :: IPT_L2OPXV_SCR(:)
  INTEGER, ALLOCATABLE :: IPT_WTHEAD_MAP(:,:)

  REAL(8), ALLOCATABLE :: IPT_A(:)
  INTEGER, ALLOCATABLE :: IPT_LMNV(:,:)       ! 01:lmnv(1:3,1:lmnvmx)
  INTEGER, ALLOCATABLE :: IPT_IL(:)           ! 02:il(1:nnbft)
! in RECP_INP_READ_AO2SO, RECP_SOCFPD
  REAL(8), ALLOCATABLE :: IPT_AO2SO(:,:,:)    ! ipt.55 -> c(1:mctu,mcsu,mgcsu)
  REAL(8), ALLOCATABLE :: IPT_CX0(:)          !* temporary before 3:cx(1:mcxu) 
  REAL(8), ALLOCATABLE :: IPT_CX(:)           ! 03:cx(1:mcxu)
! main
  REAL(8), ALLOCATABLE :: IPT_ARRAY(:,:)      ! 04:array(1:mpru,1:narray)
  INTEGER, ALLOCATABLE :: IPT_IARRAY(:,:)     ! 05:iarray(1:iamax,1:narray) 
! skip 6,7,8 
! cortab (compute tables by recursion for real spherical harmonics)
  REAL(8), ALLOCATABLE :: IPT_SH_HPT(:)       ! 09:hpt(1:ihwsp)
  REAL(8), ALLOCATABLE :: IPT_SH_HWT(:)       ! 10:hwt(1:ihwsp)
  REAL(8), ALLOCATABLE :: IPT_SH_DFAC(:)      ! 11:dfac(1:ndfac)
  REAL(8), ALLOCATABLE :: IPT_SH_BINOM(:)     ! 12:binom(1:xyzdim)
  INTEGER, ALLOCATABLE :: IPT_SH_LMF(:)       ! 13:lmf(1:l1max2)
  INTEGER, ALLOCATABLE :: IPT_SH_LML(:)       ! 14:lml(1:l1max2)
  INTEGER, ALLOCATABLE :: IPT_SH_LMX(:)       ! 15:lmx(1:lmnpwr)
  INTEGER, ALLOCATABLE :: IPT_SH_LMY(:)       ! 16:lmy(1:lmnpwr)
  INTEGER, ALLOCATABLE :: IPT_SH_LMZ(:)       ! 17:lmz(1:lmnpwr)
  REAL(8), ALLOCATABLE :: IPT_SH_ZLM(:)       ! 18:zlm(1:lmnpwr)
  REAL(8), ALLOCATABLE :: IPT_SH_FLMTX(:,:)   ! 19:flmtx(1:3,1:lmnu2)
  INTEGER, ALLOCATABLE :: IPT_SH_MC(:,:)      ! 20:mc(1:3,1:l2m1)
  INTEGER, ALLOCATABLE :: IPT_SH_MR(:,:)      ! 21:mr(1:3,1:l2m1)
! recp_one_main
  REAL(8), ALLOCATABLE :: IPT_ONE_H2(:)       ! 22:h2(1:nblu)
! skip 23,24,25 : only s,t,v
! cints,lsints
  REAL(8), ALLOCATABLE :: IPT_CZ_G2(:)        ! 26:g2(1:ic1)
  REAL(8), ALLOCATABLE :: IPT_CZ_G1(:)        ! 27:g1(1:ic0)
  REAL(8), ALLOCATABLE :: IPT_CZ_GOUT(:)      ! 28:gout(1:ij)
  REAL(8), ALLOCATABLE :: IPT_CZ_CRDA(:,:)    ! 29:crda(1:lit,1:3)
  REAL(8), ALLOCATABLE :: IPT_CZ_CRDB(:,:)    ! 30:crdb(1:ljt,1:3)
! pseud1
  REAL(8), ALLOCATABLE :: IPT_P1_ANG(:,:)     ! 31:ang(1:ltot1,1:ltot1)
  REAL(8), ALLOCATABLE :: IPT_P1_Q(:,:)       ! 32:q(1:ltot1,1:ltot1)
  REAL(8), ALLOCATABLE :: IPT_P1_QSUM(:,:)    ! 33:qsum(1:ltot1,1:ltot1)
  REAL(8), ALLOCATABLE :: IPT_P1_XAB(:)       ! 34:xab(1:ltot1)
  REAL(8), ALLOCATABLE :: IPT_P1_YAB(:)       ! 35:yab(1:ltot1)
  REAL(8), ALLOCATABLE :: IPT_P1_ZAB(:)       ! 36:zab(1:ltot1)
! pseud2 and pseud3
  REAL(8), ALLOCATABLE :: IPT_P23_ANGA(:,:,:) ! 37:anga(1:lit,1:mproju,1:lamau)
  REAL(8), ALLOCATABLE :: IPT_P23_ANGB(:,:,:) ! 38:angb(1:ljt,1:mproju,1:lambu)
  REAL(8), ALLOCATABLE :: IPT_P23_QSUM(:,:,:) ! 39:qsum(1:ltot1,1:lambu,1:lamau)
  REAL(8), ALLOCATABLE :: IPT_P23_APWR(:)     ! 40:apwr(1:ljt)
  REAL(8), ALLOCATABLE :: IPT_P23_ATERM1(:)   ! 41:aterm1(1:ljt) 
  REAL(8), ALLOCATABLE :: IPT_P23_ATERM2(:)   ! 42:aterm2(1:ljt) 
  REAL(8), ALLOCATABLE :: IPT_P23_BPREF(:)    ! 43:bpref(1:ljt)
  REAL(8), ALLOCATABLE :: IPT_P23_BPWR(:)     ! 44:bpwr(1:lit) 
  REAL(8), ALLOCATABLE :: IPT_P23_BTERM1(:)   ! 45:bterm1(1:lit) 
  REAL(8), ALLOCATABLE :: IPT_P23_SSI(:)      ! 46:ssi(ltot1+lproju)
  REAL(8), ALLOCATABLE :: IPT_P23_ABESS(:)    ! 47:abess(1:lamau)
  REAL(8), ALLOCATABLE :: IPT_P23_BBESS(:)    ! 48:bbess(1:lambu)
  REAL(8), ALLOCATABLE :: IPT_P23_PTPOW(:)    ! 49:ptpow(1:ltot1)
  REAL(8), ALLOCATABLE :: IPT_P23_Q2(:,:)     ! 50:q2(1:lambu,1:lamau)
CONTAINS


SUBROUTINE RECP_IPTA_CORTAB(ndfac,lmn1u,lproju,ncru,lmax)
! Allocate variables
  IMPLICIT NONE
  INTEGER  ndfac,lmn1u,lproju,ncru,lmax
! local variable
  INTEGER  l1max2,lmnpwr,l2m1,sum_lmn1u

  ndfac = max(4*lmn1u+2*lproju-3, 6*lproju+3, 4*lmn1u-1,  &
              2*lmn1u+2*lproju+1, 4, ncru+4*lmn1u+2*lproju-1)
  lmax =  max(1, lmn1u-1 + max(lmn1u-1,lproju))
  l1max2 = (lmax+1)**2
  lmnpwr = (((lmax*(lmax+2)*(lmax+4))/3)*(lmax+3) + (lmax+2)**2*(lmax+4))/16
  l2m1 = 2*lproju-1
  sum_lmn1u = (lmn1u*(lmn1u+1))/2
  ALLOCATE (IPT_SH_HPT(35))            ! 09:hpt(1:ihwsp)
  ALLOCATE (IPT_SH_HWT(35))            ! 10:hwt(1:ihwsp)
  ALLOCATE (IPT_SH_DFAC(ndfac))        ! 11:dfac(1:ndfac)
  ALLOCATE (IPT_SH_BINOM(sum_lmn1u))   ! 12:binom(1:xyzdim)
  ALLOCATE (IPT_SH_LMF(l1max2))        ! 13:lmf(1:l1max2)
  ALLOCATE (IPT_SH_LML(l1max2))        ! 14:lml(1:l1max2)
  ALLOCATE (IPT_SH_LMX(lmnpwr))        ! 15:lmx(1:lmnpwr)
  ALLOCATE (IPT_SH_LMY(lmnpwr))        ! 16:lmy(1:lmnpwr)
  ALLOCATE (IPT_SH_LMZ(lmnpwr))        ! 17:lmz(1:lmnpwr)
  ALLOCATE (IPT_SH_ZLM(lmnpwr))        ! 18:zlm(1:lmnpwr)
  ALLOCATE (IPT_SH_FLMTX(3,lproju**2)) ! 19:flmtx(1:3,1:lmnu2)
  ALLOCATE (IPT_SH_MC(3,l2m1))         ! 20:mc(1:3,1:l2m1)
  ALLOCATE (IPT_SH_MR(3,l2m1))         ! 21:mr(1:3,1:l2m1)
END SUBROUTINE RECP_IPTA_CORTAB

SUBROUTINE RECP_IPTD_CORTAB
! Deallocate variables
  IMPLICIT NONE
  DEALLOCATE (IPT_SH_HPT)       ! 09:hpt(1:ihwsp)
  DEALLOCATE (IPT_SH_HWT)       ! 10:hwt(1:ihwsp)
  DEALLOCATE (IPT_SH_DFAC)      ! 11:dfac(1:ndfac)
  DEALLOCATE (IPT_SH_BINOM)     ! 12:binom(1:xyzdim)
  DEALLOCATE (IPT_SH_LMF)       ! 13:lmf(1:l1max2)
  DEALLOCATE (IPT_SH_LML)       ! 14:lml(1:l1max2)
  DEALLOCATE (IPT_SH_LMX)       ! 15:lmx(1:lmnpwr)
  DEALLOCATE (IPT_SH_LMY)       ! 16:lmy(1:lmnpwr)
  DEALLOCATE (IPT_SH_LMZ)       ! 17:lmz(1:lmnpwr)
  DEALLOCATE (IPT_SH_ZLM)       ! 18:zlm(1:lmnpwr)
  DEALLOCATE (IPT_SH_FLMTX)     ! 19:flmtx(1:3,1:lmnu2)
  DEALLOCATE (IPT_SH_MC)        ! 20:mc(1:3,1:l2m1)
  DEALLOCATE (IPT_SH_MR)        ! 21:mr(1:3,1:l2m1)
END SUBROUTINE RECP_IPTD_CORTAB


SUBROUTINE RECP_IPTA_ONEINT(nblu)
  INTEGER  nblu
  ALLOCATE (IPT_ONE_H2(nblu))       ! 22:h2(1:nblu)
END SUBROUTINE RECP_IPTA_ONEINT

SUBROUTINE RECP_IPTD_ONEINT
  DEALLOCATE (IPT_ONE_H2)           ! 22:h2(1:nblu)
END SUBROUTINE RECP_IPTD_ONEINT


SUBROUTINE RECP_IPTA_STVCZ(ic0,ic1,ij,lit,ljt,ltot1,mproju,lproju,lamau,lambu)
! Allocate variables
  USE RECP_OUTPUT
  IMPLICIT NONE  
  INTEGER  ic0,ic1,ij,lit,ljt,ltot1,mproju,lproju,lamau,lambu,IOS 
  INTEGER  I,J,K
  ltot1   = lit+ljt-1
  mproju  = 2*lproju+1
  lamau   = lit+lproju
  lambu   = ljt+lproju

  ALLOCATE (IPT_CZ_G2(ic1),STAT=IOS)                  ! stvcz,26:g2(1:ic1)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_CZ_G2')
  DO I = 1, ic1
     IPT_CZ_G2(I) = 0.0d0
  ENDDO

  ALLOCATE (IPT_CZ_G1(ic0),STAT=IOS)                  ! c/lsints,27:g1(1:ic0)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_CZ_G1')
  DO I = 1, ic0 
     IPT_CZ_G1(I) = 0.0d0
  ENDDO

  ALLOCATE (IPT_CZ_GOUT(ij),STAT=IOS)                 ! c/lsints,28:gout(1:ij)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_CZ_GOUT')
  DO I = 1, ij
     IPT_CZ_GOUT(I) = 0.0d0
  ENDDO

  ALLOCATE (IPT_CZ_CRDA(lit,3),STAT=IOS)              ! c/lsints,29:crda(1:lit,1:3)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_CZ_CRDA')
  DO I = 1, 3
     DO J = 1, lit
        IPT_CZ_CRDA(J,I) = 0.0d0
     ENDDO
  ENDDO

  ALLOCATE (IPT_CZ_CRDB(ljt,3),STAT=IOS)              ! c/lsints,30:crdb(1:ljt,1:3)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_CZ_CRDB')
  ALLOCATE (IPT_P1_ANG(ltot1,ltot1),STAT=IOS)         ! pseud1,31:ang(1:ltot1,1:ltot1)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P1_ANG')
! ALLOCATE (IPT_P1_Q(ltot1,ltot1),STAT=IOS)           ! pseud1,32:q(1:ltot1,1:ltot1)
! IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P1_Q')
  ALLOCATE (IPT_P1_QSUM(ltot1,ltot1),STAT=IOS)        ! pseud1,33:qsum(1:ltot1,1:ltot1)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P1_QSUM')
  ALLOCATE (IPT_P1_XAB(ltot1),STAT=IOS)               ! pseud1,34:xab(1:ltot1)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P1_XAB')
  ALLOCATE (IPT_P1_YAB(ltot1),STAT=IOS)               ! pseud1,35:yab(1:ltot1)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P1_YAB')
  ALLOCATE (IPT_P1_ZAB(ltot1),STAT=IOS)               ! pseud1,36:zab(1:ltot1)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P1_ZAB')
  ALLOCATE (IPT_P23_ANGA(lit,mproju,lamau),STAT=IOS)  ! pseud2/3,37:anga(1:lit,1:mproju,1:lamau)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P23_ANGA')
  ALLOCATE (IPT_P23_ANGB(ljt,mproju,lambu),STAT=IOS)  ! pseud2/3,38:angb(1:ljt,1:mproju,1:lambu)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P23_ANGB')
  ALLOCATE (IPT_P23_QSUM(ltot1,lambu,lamau),STAT=IOS) ! pseud2/3,39:qsum(1:ltot1,1:lambu,1:lamau)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P23_QSUM')
  ALLOCATE (IPT_P23_APWR(ljt),STAT=IOS)               ! qbess,40:apwr(1:ljt)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P23_APWR')
  ALLOCATE (IPT_P23_ATERM1(ljt),STAT=IOS)             ! qbess,41:aterm1(1:ljt) 
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P23_ATERM1')
  ALLOCATE (IPT_P23_ATERM2(ljt),STAT=IOS)             ! qbess,42:aterm2(1:ljt) 
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P23_ATERM2')
  ALLOCATE (IPT_P23_BPREF(ljt),STAT=IOS)              ! qbess,43:bpref(1:ljt)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_')
  ALLOCATE (IPT_P23_BPWR(lit),STAT=IOS)               ! qbess,44:bpwr(1:lit) 
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_')
  ALLOCATE (IPT_P23_BTERM1(lit),STAT=IOS)             ! qbess,45:bterm1(1:lit) 
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_')
  ALLOCATE (IPT_P23_SSI(ltot1+lproju),STAT=IOS)       ! qbess,46:ssi(ltot1+lproju)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_')
  ALLOCATE (IPT_P23_ABESS(lamau),STAT=IOS)            ! ptwt, 47:abess(1:lamau)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_')
  ALLOCATE (IPT_P23_BBESS(lambu),STAT=IOS)            ! ptwt, 48:bbess(1:lambu)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_')
  ALLOCATE (IPT_P23_PTPOW(ltot1),STAT=IOS)            ! ptwt, 49:ptpow(1:ltot1)
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P23_PTPOW')
  ALLOCATE (IPT_P23_Q2(lambu,lamau),STAT=IOS)         ! ptwt, 50:q2(1:lambu,1:lamau) 
  IF(IOS.NE.0) CALL RECP_EXIT_IPT('IPT_P23_Q2')

END SUBROUTINE RECP_IPTA_STVCZ

SUBROUTINE RECP_IPTD_STVCZ
  IMPLICIT NONE
  DEALLOCATE (IPT_CZ_G2)         ! 26:g2(1:ic1)
  DEALLOCATE (IPT_CZ_G1)         ! 27:g1(1:ic0)
  DEALLOCATE (IPT_CZ_GOUT)       ! 28:gout(1:ij)
  DEALLOCATE (IPT_CZ_CRDA)       ! 29:crda(1:lit,1:3)
  DEALLOCATE (IPT_CZ_CRDB)       ! 30:crdb(1:ljt,1:3)
  DEALLOCATE (IPT_P1_ANG)        ! pseud1,31:ang(1:ltot1,1:ltot1)
! DEALLOCATE (IPT_P1_Q)          ! pseud1,32:q(1:ltot1,1:ltot1)
  DEALLOCATE (IPT_P1_QSUM)       ! pseud1,33:qsum(1:ltot1,1:ltot1)
  DEALLOCATE (IPT_P1_XAB)        ! pseud1,34:xab(1:ltot1)
  DEALLOCATE (IPT_P1_YAB)        ! pseud1,35:yab(1:ltot1)
  DEALLOCATE (IPT_P1_ZAB)        ! pseud1,36:zab(1:ltot1)
  DEALLOCATE (IPT_P23_ANGA)      ! pseud2/3,37:anga(1:lit,1:mproju,1:lamau)
  DEALLOCATE (IPT_P23_ANGB)      ! pseud2/3,38:angb(1:ljt,1:mproju,1:lambu)
  DEALLOCATE (IPT_P23_QSUM)      ! pseud2/3,39:qsum(1:ltot1,1:lambu,1:lamau)
  DEALLOCATE (IPT_P23_APWR)      ! qbess,40:apwr(1:ljt)
  DEALLOCATE (IPT_P23_ATERM1)    ! qbess,41:aterm1(1:ljt) 
  DEALLOCATE (IPT_P23_ATERM2)    ! qbess,42:aterm2(1:ljt) 
  DEALLOCATE (IPT_P23_BPREF)     ! qbess,43:bpref(1:ljt)
  DEALLOCATE (IPT_P23_BPWR)      ! qbess,44:bpwr(1:lit) 
  DEALLOCATE (IPT_P23_BTERM1)    ! qbess,45:bterm1(1:lit) 
  DEALLOCATE (IPT_P23_SSI)       ! qbess,46:ssi(ltot1+lproju)
  DEALLOCATE (IPT_P23_ABESS)     ! ptwt, 47:abess(1:lamau)
  DEALLOCATE (IPT_P23_BBESS)     ! ptwt, 48:bbess(1:lambu)
  DEALLOCATE (IPT_P23_PTPOW)     ! ptwt, 49:ptpow(1:ltot1)
  DEALLOCATE (IPT_P23_Q2)        ! ptwt, 50:q2(1:lambu,1:lamau)
END SUBROUTINE RECP_IPTD_STVCZ

END MODULE RECP_IPT

!     # allocate space.
!     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(1:mcxu),
!     #          4:array(1:mpru,1:narray), 5:iarray(1:iamax,1:narray),
!     #          6:buffer(1:l1rec), 7:values(1:n1max),
!     #          8:labels(1:2,1:n1max), 9:hpt(1:ihwsp), 10:hwt(1:ihwsp),
!     #          11:dfac(1:ndfac), 12:binom(1:xyzdim), 13:lmf(1,l1max2),
!     #          14:lml(1:l1max2), 15:lmx(1:lmnpwr), 16:lmy(1:lmnpwr),
!     #          17:lmz(1:lmnpwr), 18:zlm(1:lmnpwr),
!     #          19:flmtx(1:3,1:lmnu2), 20:mc(1:3,1:l2m1),
!     #          21:mr(1:3,1:l2m1), 22:h2(1:nblu), 23:ijx(1:iiju),
!     #          24:ijy(1,iiju), 25:ijz(1:iiju), 26:g2(1:ic1),
!     #          27:g1(1:ic0), 28:gout(1:ij),
!     #        if s, t, or v integrals, then
!     #          29:xin(1:nxyzin), 30:yin(1:nxyzin), 31:zin(1:nxyzin),
!     #          32:uf(1:nroots), 33:wf(1,nroots), 34:ff(1:nroots+n1),
!     #          (cf and sf must be kept together for droot)
!     #          35:cf(1:n1,1:n1), 36:sf(1:n1,1:n1), 37:af(1:n1),
!     #          38:rt(1:n1), 39:r(1:nroots,1:nroots),
!     #          40:w(1:nroots,1:nroots)
!     #        elseif c or z integrals, then
!     #          29:crda(1:lit,1:3), 30:crdb(1:ljt,1:3),
!     #          (for pseud1)
!     #          31:ang(1:ltot1,1:ltot1), 32:q(1:ltot1,1:ltot1),
!     #          33:qsum(1:ltot1,1:ltot1), 34:xab(1:ltot1),
!     #          35:yab:(1:ltot1), 36:zab(1:ltot1),
!     #          (for pseud2 and pseud3)
!     #          37:anga(1:lit,1:mproju,1:lamau),
!     #          38:angb(1:ljt,1:mproju,1:lambu),
!     #          39:qsum(1:ltot1,1:lambu,1:lamau), 40:apwr(1:ljt),
!     #          41:aterm1(1:ljt), 42:aterm2(1:ljt), 43:bpref(1:ljt),
!     #          44:bpwr(1:lit), 45:bterm1(1:lit), 46:ssi(ltot1+lproju),
!     #          47:abess(1:lamau), 48:bbess(1:lambu),
!     #          49:ptpow(1:ltot1), 50:q2(1:lambu,1:lamau)
!     #        endif
