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

MODULE RECP_INP
CONTAINS
SUBROUTINE RECP_INP_MAIN(ccr,chg,eta,ica,idp,ityp,la,lb,lcr,  &
           lls,lmnp1,lmn1u,lproju,maords,mcons,mcrs,mgcs,mnl,ms,mtype,nblpr,nc, &
           ncon,ncr,ncru,nct,nd,nf,nir,nkcrl,nkcru,nklsl,nklsu,nrcr,nso,nsopr,  &
           nt,ntl,ntu,x,y,z,zcr,zet,CX_TEMP)
! read and verify the user input: symmetry info, basis sets,
! geometry, etc.
  USE RECP_INP_READ
  USE RECP_INP_ORBITAL 
  USE RECP_INP_PRINT
  USE RECP_IPT
  USE RECP_NTR
  USE RECP_OUTPUT
  USE RECP_FUNCTION1
  USE RECP_FUNCTION2
  USE RECP_WRITE
  implicit logical(a-z)
#include "inc_mxvalue.h"
#include "inc_print.h" 
! /bufout/ holds some output integral file parameters.
  integer         itypea,itypeb,ibuf,numout,nrec,ntape
  common /bufout/ itypea,itypeb,ibuf,numout,nrec,ntape
  real*8         tol
  common /parmr/ tol
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  integer        mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,   inam,   nnam,mdum
  common /ntgr/  mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir(3),inam(5),nnam,mdum(32)
  integer ierr
  integer ica(mcu,msu,*),idp(mstu,mstu,*),la(mru,*),lb(*),lcr(*),lls(*),  &
          lmnp1(*),lmn1u,lproju,maords(*),mcons(*),mcrs(*),mgcs(*),mnl(*),ms(*),   &
          nblpr(*),nc(*),ncon(*),ncr(*),ncru,nct(*),nd(*),nf(*),nir(*),nkcrl(6,*), &
          nkcru(6,*),nklsl(4,*),nklsu(4,*),nrcr(*),nso(*),nsopr(*),nt(*),ntl(*),ntu(*)
  real*8  ccr(*),chg(*),eta(mrcru,mconu,*),x(mcu,*),y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
  character*3 ityp(*), mtype(*)
! local variables
  integer I,J,K,itol,inrm,ncrs,ngen,ncons,naords,iaords,iru,igcs,icsu, &
          ic,icons,iconu,ircru,icon,l,isf,isfr,is,icu,ig,ilmnp1,itl,itu,jcts,it,jct,jt,js
  integer ixyzir(3), lblso(16),CX_TEMP
  real*8  t, eval
  character*1 ibl, ipc, isk, lblsh(21)
  INTEGER  FILE_NUM(2)
  CHARACTER(10) FILE_NAME(2)

  DATA FILE_NUM  / 64, 65 /
  DATA FILE_NAME /'RECP_INP_0','RECP_INP_P'/
  DATA ibl, isk /' ','0'/ 
  DATA lblsh /'s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u','v','w','x','y','z'/

  CALL OUTPUT_LOC('RECP_INP_MAIN','E')
! write file in debug mode
  IF (RECP_DBG.GE.1) THEN
     CALL RECP_WRITE_FILEOPEN(FILE_NUM(1),FILE_NAME(1),'REPLACE')
  ENDIF

! intitalize some input variables.
  itol   = 20
  tol    = 2.30258d0 * itol   ! ln10 *itol

! ----------------------------------------------------------------
  CALL OUTPUT_DBG_TITLE('INP_MAIN : start process input data')

! # read general information 
  CALL RECP_INP2_READGEN(ngen,ns,naords,ncons,ngcs,ncrs)

! # Number of distinct products of irrep 
  CALL RECP_INP2_READIRREP(idp,nd,ityp)

! # Read AO to SO transformation matrices 
  CALL RECP_INP2_READAO2SO(naords,ngcs,maords,la,nir,CX_TEMP)

! # read in exponents and contraction coefficients.
  CALL RECP_INP2_READORB(zet,eta,lmnp1,nrcr,ncon)

! read in expansions for core and spin-orbit potentials.
  CALL RECP_INP2_READCSO(lcr,lls,nkcrl,nkcru,nklsl,nklsu,ncr,ncrs,ncru,lproju,ns,zcr,ccr)

  ng = ngen + 1

! read in atomic labels, charges, coordinates
  CALL RECP_INP2_READGEO(ns,ng,mcu,msu,nf,nc,ica,X,Y,Z,chg,mtype,nst,nso,nbft, &
       mcons,mgcs,nt,ntl,ntu,la,lb,nrcr,nir,lmn1u,lmnp1,maords,ncrs,mcrs)

  DEALLOCATE ( RECP_SET ) 
  CALL OUTPUT_DBG_TITLE('INP_MAIN : end process input data')
! ----------------------------------------------------------------

! generate cartesian gaussian exponent array. (mklmnv)
  CALL RECP_INP_CARTEXP
 
! determine the number of types of one-electron integrals.
  CALL RECP_INP_INTCHECK(nu,ns,ncrs,mcrs,lcr,lls,inam) 
 
! normalization of symmetry orbitals. (deleted)
! if (inrm.ne.0) CALL RECP_INP_ORB_SONORM
 
! contracted orbital and lower-triangle-packed matrix symmetry offsets.
  CALL RECP_INP_ORBITAL_OFFSET(nsopr,nso,nblpr,nst,FILE_NUM)

! nbft  = total number of basis functions. 
! nnbft = total number of unique so function pairs.
! mpru  = the number of elements in a diagonal,
!         symmetry-blocked, lower-triangle-packed matrix.  this
!         is also an upper bound to the size of any
!         symmetry-blocked array in which only unique
!         off-diagonal elements are stored and which is
!         associated with an operator that transforms as an
!         irreducible representation.
  nnbft = (nbft * (nbft + 1)) / 2
  mpru  = nblpr(nst) + (nso(nst)+1)*((nso(nst)+1)-1)/2

! normalize the contraction coefficients.
  CALL RECP_INP_ORBITAL_CCNORM(mrcru,mconu,ncons,nrcr,ncon,lmnp1,eta,zet)

  CALL gcentr(ica,nc)
 
  IF ((ng+1).GT.mgup) call bummer('change mgup (two places) to ',(ng+1),2)
! IF (PRINT_LEVEL.GE.2) THEN
!    WRITE (RECP_OUT,'(A,A,I2,A)')' * primitive ao integrals neglected', &
!          ' if exponential factor below 10**(-',itol,')'
! ENDIF 

  CALL OUTPUT_DBG_TITLE('INP_MAIN : start print input info')
  CALL RECP_SETZERO_I1(ixyzir,3)

! ----------------------------------------------------------------
  isf = 0
  isfr = 0
  DO is = 1,ns
     icu = nc(is)

!    INFO1 : atom name/charge/geometry/center interchanges
     CALL RECP_INP_INFO1(is,nc,mcu,msu,ica,ng,ngen,chg,x,y,z,mtype)

     DO J = 1,nf(is)
        isf = isf+1
        itl = ntl(isf)
        itu = ntu(isf)
        igcs = mgcs(isf)
        iaords = maords(igcs)
        iru = nir(iaords)
        icons = mcons(isf)
        ircru = nrcr(icons)
        iconu = ncon(icons)
        icsu = 0
        DO K = 1,iru
           icsu = icsu + nd(la(K,iaords))
        ENDDO

!       INFO2 : basisset/symmetry orbital labels/symmetry orbitals
        CALL RECP_INP_INFO2(lmnp1,icons,ircru,iconu,zet,eta,icsu,igcs, &
             isfr,iru,la,nd,lb,ityp,icu,itl,itu,IPT_LMNV,iaords)
!       INFO3 : orbital 
        CALL RECP_INP_INFO3(mru,ircru,iru,isfr,iaords,la,nsopr,lb,ms,is,ixyzir, &
             icons,lmnp1,itl,itu,igcs,mnl,nd,nu,inam,icu,eval,ityp)
     ENDDO
     CALL RECP_INP_PRINTCSO(is,mcrs,lcr,nkcrl,nkcru,ncr,nklsl,nklsu,lls,zcr,ccr,mtype,lblsh)
  ENDDO 
  IF (PRINT_LEVEL.GE.2) WRITE(RECP_OUT,'(/,10X,60A,/)') ('=', I=1,60)
! ----------------------------------------------------------------

! multiply normalization constants into contraction coefficients.
  CALL RECP_INP_ORBITAL_MULTNORMCOEF(ncons,ncon,nrcr,mrcru,mconu,zet,eta,lmnp1)

! set lxyzir
  CALL RECP_INP_ORBITAL_SETLXYZIR(nu,inam,nst,ixyzir,lxyzir,mstu,idp,ityp)
 
! write orbital labels 
  CALL RECP_INP_PRINT_ORBLABEL(nst,ns,isfr,ityp,nso,mtype,ms,mnl)

  IF (RECP_DBG.GE.1) CLOSE(74) !close file (RECP_DBG_INP)
  CALL OUTPUT_LOC('RECP_INP_MAIN','X')
END SUBROUTINE RECP_INP_MAIN
END MODULE RECP_INP
