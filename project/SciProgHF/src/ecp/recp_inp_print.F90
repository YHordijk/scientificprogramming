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

MODULE RECP_INP_PRINT
CONTAINS

SUBROUTINE RECP_INP_PRINT_ORBLABEL(nsym,ns,nbft,ityp,nso,mtype,ms,mnl)
!  input:
!   nsym = number of symmetry blocks.
!   ns   = number of symmetry inequivalent atom labels.
!   nbft = total number of basis functions.
!   repnuc = nuclear repulsion energy.
!   ityp(1:nsym) = input symmetry labels.
!   nso(1:nsym) = number of basis functions in each symmetry block.
!   mtype(1:ns) = atom labels.
!   ms(1:nbft) = basis_function-to-center mapping vector.
!   mnl(1:nbft) = integer codes for basis function types.
!                 1:s, 2:p, 3:3s, 4:d, 5:4p, 6:f, 7:5s, 8:5d, 9:g
!   map(1:2,1:nbft) = temporary scratch array for writing the output map(*) vectors.
  USE RECP_IPT
  USE RECP_FUNCTION2
  implicit logical(a-z)
#include "inc_print.h"
  integer    nbfmx      
  parameter( nbfmx=1000 )  !nbfmx=500->1000
  integer    nsym,ns,nbft,nso(*),ms(1:nbft),mnl(1:nbft)
  character*3   ityp(1:nsym),mtype(1:ns)
! local variables
  integer    bmin,   bmax
  parameter( bmin=0, bmax=17 )
  integer      clip,i,j
  character*8  bfnlab(nbfmx)
  character*4  slabel(8)
  character*2  btype(bmin:bmax)
! # basis function types.
  data btype/' ?','_s','_p','3s','_d','4p','_f','5s','5d',  &
             '_g','6p','6f','_h','7s','7d','7g','_i',' ?' /
 
  ALLOCATE(IPT_WTHEAD_MAP(nbft,2))

! Check if basis function number exceeds maxium
  IF (nbft.GT.nbfmx) THEN
     WRITE(RECP_OUT,'(A)') ' * Error : In this version of ECP integral routine,'
     WRITE(RECP_OUT,'(3X,A,X,I4)') 'the number of basis function are limited to',nbfmx
     WRITE(RECP_OUT,'(3X,I4,X,A)') (nbft-nbfmx),'basis functions are exceed'
     CALL QUIT('Number of basis function are exceed.')
  ENDIF
 
! Set map(*,*) and create bfnlab(*) from the input arrays.
  DO i = 1, nbft
     IPT_WTHEAD_MAP(i,1) = ms(i)
     IPT_WTHEAD_MAP(i,2) = mnl(i)
     clip = min( max( bmin, mnl(i) ), bmax )
     WRITE(bfnlab(i),'(i3,a3,a2)') i, mtype(ms(i)), btype(clip)
     DO j = 4, 8
!       replace embedded spaces with a printing character to
!       avoid confusion when printing labels later.
        if (bfnlab(i)(j:j).eq.' ') bfnlab(i)(j:j) = '_'
     ENDDO   
  ENDDO    
 
! Create the output symmetry labels.
  DO i = 1, nsym
     slabel(i)(1:1) = ' '
     slabel(i)(2:4) = ityp(i)
  ENDDO

! =============
! PRINT SECTION 
! =============
  IF (PRINT_LEVEL.GE.2) THEN

!    basic infomation
!    ----------------
     WRITE(RECP_OUT,'(/,a,2x,i2)')  ' * symmetry operator No. : ',nsym
     WRITE(RECP_OUT,'(a,i4)')       ' * basis function No.    : ',nbft
     WRITE(RECP_OUT,'(a,8i5)')      ' * symmetry number       : ',(i,i=1,nsym)
     WRITE(RECP_OUT,'(a,8(1x,a4))') ' * symmetry label        : ',(slabel(i),i=1,nsym)
     WRITE(RECP_OUT,'(a,8i5)')      ' * basis function        : ',(nso(i),i=1,nsym)

!    orbital labels
!    --------------
     WRITE(RECP_OUT,'(/,a)') ' * output orbital labels (i:bfnlab(i))'
     WRITE(RECP_OUT,'(5(i4,a1,a8))') (i,':',bfnlab(i),i=1,nbft)

!    basis function to atomic center map
!    -----------------------------------
     WRITE(RECP_OUT,'(/,a)') ' * basis function to center (i:map(i))'
     WRITE(RECP_OUT,'(5(i4,a1,i4,4x))') (i,':',IPT_WTHEAD_MAP(i,1),i=1,nbft)

!    basis function to orbital type map
!    ----------------------------------
     WRITE(RECP_OUT,'(/,a)') ' * basis function to orbital type (i:map(i))'
     WRITE(RECP_OUT,'(5(i4,a1,i4,4x))') (i,':',IPT_WTHEAD_MAP(i,2),i=1,nbft)

!    symmetry orbitals
!    -----------------
     WRITE(RECP_OUT,'(/,I9,A,3X,4(3X,A3,A,I4),(/,31X,4(3X,A3,A,I4)))') &
          nbft,' symmetry orbitals,', (ityp(I),':',nso(I), I=1,nsym)
  ENDIF

  DEALLOCATE(IPT_WTHEAD_MAP)
END SUBROUTINE RECP_INP_PRINT_ORBLABEL


SUBROUTINE RECP_INP_INFO1(is,nc,mcu,msu,ica,ng,ngen,chg,x,y,z,mtype)
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER     is,nc(*),mcu,msu,ica(mcu,msu,*),ng,ngen
  REAL*8      chg(*),x(mcu,*),y(mcu,*),z(mcu,*)
  CHARACTER*3 mtype(*)
  INTEGER I,J  !local variables

  IF (PRINT_LEVEL.LE.1) RETURN

  WRITE(RECP_OUT,'(/,10X,60A)') ('=', I=1,60)
  WRITE(RECP_OUT,'(/36X,A3,A)') mtype(is),' atoms'
  WRITE(RECP_OUT,'(/30X,A,F7.2)') 'nuclear charge',chg(is)
  WRITE(RECP_OUT,'(/11x,A,12X,A,15X,A,15X,A)') 'center','x','y','z'
  WRITE(RECP_OUT,'(I14,7X,3F16.8)') (I, x(I,is), y(I,is), z(I,is), I=1,nc(is))

  IF (nc(is).GT.1) THEN
     WRITE (RECP_OUT,'(A,6X,A)') ' operator','center interchanges'
     WRITE (RECP_OUT,'(A,6X,24I3)') '  1(id.)',(ica(I,is,1), I=1,nc(is))

     DO J = 2,(ngen+1)
        WRITE (RECP_OUT,'(I3,A,5X,24I3)') J,'(gen.)',(ica(I,is,J), I=1,nc(is))
     ENDDO

     DO J = (ngen+2),ng
        WRITE(RECP_OUT,'(I3,11X,24I3)') J,(ica(I,is,J), I=1,nc(is))
     ENDDO
  ENDIF
END SUBROUTINE RECP_INP_INFO1


SUBROUTINE RECP_INP_INFO2(LMNP1,ICONS,IRCRU,ICONU,ZET,ETA,ICSU,IGCS,  &
                ISFR,IRU,LA,ND,LB,ITYP,ICU,ITL,ITU,LMNV,IAORDS)
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  integer        mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,inam,nnam,mdum
  common /ntgr/  mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir(3),inam(5),nnam,mdum(32)
  INTEGER     LMNP1(*), I
  REAL*8      ZET(MCONU,*), ETA(MRCRU,MCONU,*)
  INTEGER     ICONS, ICONU, ICON
  INTEGER     IRCRU, IRCRH, IRCR, IRCRL
  INTEGER     JSFR, ISFR, IR, IRU
  INTEGER     LAI, LA(MRU,*), IAORDS, IA, IAU, ND(*), LB(*)
  INTEGER     ICSR, ICOL, NCOL, LBLA(16), ICT, IC, IT
  INTEGER     ICU, ITL, ITU, LMNV(3,*), ICSU, IGCS
  CHARACTER*3 ITYP(*), LTYP(16)
  CHARACTER*1 LBLSH(21)
  DATA        LBLSH  /'s','p','d','f','g','h','i','k','l','m','n',  &
                      'o','q','r','t','u','v','w','x','y','z'/

  IF (PRINT_LEVEL.LE.1) RETURN

  WRITE(RECP_OUT,'(/,10X,60A)') ('-', I=1,60)

! orbital exponents  contraction coefficients for each orbital
! ------------------------------------------------------------
! WRITE(RECP_OUT,'(/,I21,A1,X,A)') lmnp1(icons), LBLSH(lmnp1(icons)),'orbitals'
  WRITE(RECP_OUT,'(/,36X,A1,X,A)') LBLSH(lmnp1(icons)),'orbitals'
  WRITE(RECP_OUT,'(/,10X,A,X,A)') 'orbital exponents',' contraction coefficients'
! use larger value than 4 for output lines longer than 80
  ircrh = min(4,ircru)
  DO icon = 1,iconu
     WRITE (RECP_OUT,'(10X,8(1pg16.7))') zet(icon,icons),(eta(ircr,icon,icons), ircr=1,ircrh)
  ENDDO   

! symmetry orbital labels
! -----------------------
  write (RECP_OUT,"(/,10X,'symmetry orbital labels')")
  jsfr = isfr
  do ir = 1,iru
     jsfr = jsfr + 1
     lai = la(ir,iaords)
     do ia = 1,nd(lai)
        write (RECP_OUT,'(10x,7(i5,a3,i1))')  & 
              (lb(jsfr+(ircr-1)*icsu),ityp(lai),ia, ircr=1,ircrh)
     enddo
  enddo

  724 CONTINUE
  if (ircrh.lt.ircru) then
     ircrl = ircrh + 1
     ircrh = min(ircrh+4,ircru)
     write (6,"(/20x,'contraction coefficients')")
     do icon = 1,iconu
        write (RECP_OUT,'(16x,7(1pg16.7))') (eta(ircr,icon,icons), ircr=ircrl,ircrh)
     enddo   
     write (RECP_OUT,"(/21x,'symmetry orbital labels')")
     jsfr = isfr
     do ir = 1,iru
        jsfr = jsfr + 1
        lai = la(ir,iaords)
        do ia = 1,nd(lai)
           write (RECP_OUT,'(10x,7(i12,a3,i1))')  &
                 (lb(jsfr+(ircr-1)*icsu),ityp(lai),ia, ircr=ircrl,ircrh)
        enddo   
     enddo   
     goto 724
  endif

! symmetry orbitals
! -----------------
  WRITE (RECP_OUT,"(/10x,'symmetry orbitals')")
  icsr = 0
  ncol = 0
  do ir = 1,iru
     lai = la(ir,iaords)
     iau = nd(lai)
     do ia = 1,iau
        ncol = ncol+1
        ltyp(ncol) = ityp(lai)
        lbla(ncol) = ia
        if ( ncol.ne.10 .and. (ir.ne.iru .or. ia.ne.iau) ) goto 756
        write(RECP_OUT,'(10X,A,16(3X,A3,I1))') 'ctr, ao',  &
             (ltyp(icol),lbla(icol), icol=1,ncol)
        ict = 0
        do ic = 1,icu
           do it = itl,itu
              ict = ict+1
              write(RECP_OUT,'(9X,i3,a2,3i1,16f7.3)') ic,', ',(lmnv(i,it), i=1,3),  &
                   (IPT_AO2SO(ict,icsr+icol,igcs), icol=1,ncol)
           enddo   
        enddo   
        icsr = icsr+ncol
        ncol = 0
        756 continue
     enddo
  enddo   
END SUBROUTINE RECP_INP_INFO2


SUBROUTINE RECP_INP_INFO3(mru,ircru,iru,isfr,iaords,la,nsopr,lb,ms,is,ixyzir, &
           icons,lmnp1,itl,itu,igcs,mnl,nd,nu,inam,icu,eval,ityp)
  USE RECP_INP_ORBITAL
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER  mru,ircru,iru,isfr,iaords,la(mru,*),nsopr(*),lb(*),ms(*),is,ixyzir(3)
  INTEGER  icons,lmnp1(*),itl,itu,igcs,mnl(*),nd(*),nu,inam(5),icu
  REAL(8)  eval
  CHARACTER ityp(*)
! local variables
  INTEGER  I,J,K,ics,la0,so_bfn,lmn,leig

  IF (PRINT_LEVEL.GE.2) WRITE(RECP_OUT,*)  !for the blank line

  DO I = 1,ircru
     ics = 0
     DO J = 1,iru
        isfr   = isfr+1
        la0    = la(J,iaords)
        so_bfn = nsopr(la0)+lb(isfr)
        ms(so_bfn) = is
 
!       # compute the L^2-operator expectation value and
!       # determine if it is an eigenvalue.
        lmn = lmnp1(icons) - 1
        CALL RECP_INP_ORBITAL_L2OPXV(lmn,itl,itu,IPT_LMNV,eval,leig,ics,igcs)

!       # not an eigenvalue. print a warning.
        IF (leig.LT.0) THEN
           IF (PRINT_LEVEL.GE.2) THEN
              WRITE(RECP_OUT,'(3x,a,i3,a3,a,f12.8)') &
              's.o. basis function',so_bfn,ityp(la0), &
              ' is not an L^2 eigenvector. : <v|L^2|v>/<v|v> =', eval
           ENDIF 
           leig = lmn  !assign the maximum value
        ENDIF

!       # assign the mnl(*) values to the so functions.
!       # 1:s, 2:p, 3:3s, 4:d, 5:4p, 6:f, 7:5s, 8:5d, 9:g,...
        mnl(so_bfn) = ( lmnp1(icons)**2 + 2 * leig + 5 ) / 4

        DO K = 1, nd(la0)
           ics = ics+1
           IF ( (inam(nu).ge.5).and.(leig.eq.1) ) THEN
!             # determine angular momentum symmetries.
              CALL RECP_INP_ORBITAL_ANGSYM(icu,itl,itu,ics,igcs,la0,ixyzir)
           ENDIF
        ENDDO
     ENDDO    
  ENDDO
END SUBROUTINE RECP_INP_INFO3

END MODULE RECP_INP_PRINT
