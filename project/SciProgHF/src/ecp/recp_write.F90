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

MODULE RECP_WRITE
CONTAINS 

SUBROUTINE RECP_WRITE_ISA(ircru,jrcru,nst,mstu,nprir,ijsf,nopir,ilxyz,iapt, &
                          mpru,ipr,N_INTEGRAL,h2,esf)

! # put integrals and labels into the intermediate storage arrays.
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
! global variables
  INTEGER :: ircru,jrcru,nst,mstu,nprir(2,mstu,*),ijsf,nopir(*),ilxyz(3,*),iapt(3)
  INTEGER :: mpru,ipr
  INTEGER :: N_INTEGRAL(3)     ! number of integral element of each type 
  REAL(8) :: h2(*)
  LOGICAL :: esf
! local variables
  INTEGER :: ist,iop,INDEX1,NPAIR,IPAIR,I,J,INT_TYPE 

  INDEX1 = 0
   
  do I = 1, ircru
    if(esf) jrcru = I
    do J = 1, jrcru
      do ist = 1, nst
!       setting NPAIR ?
        if ( esf .and. J.eq.I ) then
          NPAIR = nprir(1,ist,ijsf)
        else
          NPAIR = nprir(2,ist,ijsf)
        endif

        if ( (nopir(ist).ne.0) .and. (NPAIR .ne. 0) ) then
          do iop = 1, nopir(ist)
!           INT_TYPE : integral type: S,T,V = 1 / SO(x,y,z) = 1,2,3
            INT_TYPE = ilxyz(iop,ist)
!           print *,'ISA-DOLOOP',INDEX1,I,J,ist,iop,NPAIR
            do IPAIR = 1, NPAIR
              INDEX1 = INDEX1 + 1
              iapt(INT_TYPE) = iapt(INT_TYPE) + 1
              IF (RECP_DBG.GE.1) THEN
                 IF (INT_TYPE.EQ.1) WRITE(64,'(x,f30.16)') h2(INDEX1)
                 IF (INT_TYPE.EQ.2) WRITE(65,'(x,f30.16)') h2(INDEX1)
                 IF (INT_TYPE.EQ.3) WRITE(66,'(x,f30.16)') h2(INDEX1)
                 IF (INT_TYPE.EQ.1) WRITE(67,'(x,i10)') (ipr+IPAIR) 
                 IF (INT_TYPE.EQ.2) WRITE(68,'(x,i10)') (ipr+IPAIR)
                 IF (INT_TYPE.EQ.3) WRITE(69,'(x,i10)') (ipr+IPAIR)
              ELSE
                 IF (INT_TYPE.EQ.1) WRITE(64,*) h2(INDEX1)
                 IF (INT_TYPE.EQ.2) WRITE(65,*) h2(INDEX1)
                 IF (INT_TYPE.EQ.3) WRITE(66,*) h2(INDEX1)
                 IF (INT_TYPE.EQ.1) WRITE(67,*) (ipr+IPAIR) 
                 IF (INT_TYPE.EQ.2) WRITE(68,*) (ipr+IPAIR)
                 IF (INT_TYPE.EQ.3) WRITE(69,*) (ipr+IPAIR)
              ENDIF
              N_INTEGRAL(INT_TYPE) = N_INTEGRAL(INT_TYPE) + 1
            enddo
          enddo
        endif

        ipr = ipr + NPAIR

      enddo
    enddo
  enddo

END SUBROUTINE RECP_WRITE_ISA

! ------------------------------------------------

SUBROUTINE RECP_WRITE_MAIN(nop,iapt,ipr,ipr1,ipr2,n_bfn_sym, &
                           N_INTEGRAL,MN_INTEGRAL,INT_AREP,nnam,mpru,nbft,nnbft,nst)
  USE RECP_IPT
  USE RECP_FUNCTION2
  USE RECP_OUTPUT
! write 1-e integral arrays, Y.C. Park
  IMPLICIT NONE
! variables
#include "inc_mxvalue.h"
#include "inc_print.h"

! global variable
  INTEGER :: mpru, nbft, nnbft, nst   ! common block (parmi)
  INTEGER :: nnam                     ! common block (ntgr)
  INTEGER :: INT_AREP
  INTEGER :: N_INTEGRAL(3) ! number of integral element of each type
  INTEGER :: MN_INTEGRAL   ! maximum of N_INTEGRAL(i,i=1-3)
  INTEGER :: nop, iapt(3), ipr, ipr1, ipr2

! local variable
  REAL(8), ALLOCATABLE :: array_stvc(:)
  REAL(8), ALLOCATABLE :: array_so(:,:)
  REAL(8), ALLOCATABLE :: values_stvc(:)    ! to be fixed
  INTEGER, ALLOCATABLE :: labels_stvc(:,:)  ! from labels(nipv,*)

  INTEGER :: INT_TYPE     ! integral type 
  INTEGER :: ilfact = 1024
  INTEGER :: itypea, itypeb
  INTEGER :: lblop
  INTEGER :: n_bfn_sym(mstup)  ! (nbpsy or nso)
  INTEGER :: num, I,value_temp,DBG1
  CHARACTER(8) :: chrtyp
  CHARACTER(11) :: FILEFORM(2)
  DATA FILEFORM /'UNFORMATTED','FORMATTED  '/
  
  INTEGER :: itypan(7), itypbn(7)
! #             s    t    v  veff so:x so:y so:z
  data itypan / 0,   0,   0,   0,   2,   2,   2 /
  data itypbn / 0,   1,   2,   3,   0,   1,   2 /

  CALL OUTPUT_LOC('RECP_WRITE_MAIN','E')
  DBG1=1
  IF (RECP_DBG.GE.1) DBG1=2

  IF (nnam.EQ.5) THEN
     CALL RECP_WRITE_FILEOPEN(74,'RECP_ISA1','OLD')
     CALL RECP_WRITE_FILEOPEN(75,'RECP_ISA2','OLD')
     CALL RECP_WRITE_FILEOPEN(76,'RECP_ISA3','OLD')
     CALL RECP_WRITE_FILEOPEN(77,'RECP_ISB1','OLD')
     CALL RECP_WRITE_FILEOPEN(78,'RECP_ISB2','OLD')
     CALL RECP_WRITE_FILEOPEN(79,'RECP_ISB3','OLD')
     REWIND(74)
     REWIND(75)
     REWIND(76)
     REWIND(77)
     REWIND(78)
     REWIND(79)
  ELSE
     CALL RECP_WRITE_FILEOPEN(74,'RECP_ISA0','OLD')
     CALL RECP_WRITE_FILEOPEN(77,'RECP_ISB0','OLD')
     REWIND(74)
     REWIND(77)
  ENDIF

  ALLOCATE( array_stvc(MN_INTEGRAL) )
  ALLOCATE( array_so(nbft,nbft) )
  ALLOCATE( values_stvc(MN_INTEGRAL) )
  ALLOCATE( labels_stvc(2,MN_INTEGRAL) )

  IF (PRINT_LEVEL.GE.5) WRITE (RECP_OUT,'(A,I3)') ' nop :',nop
! # write out the components of the operator arrays.
  DO lblop = 1, nop

!    # check integral type
     itypea = itypan(nnam + lblop - 1)
     itypeb = itypbn(nnam + lblop - 1)
     call siftyp( itypea, itypeb, chrtyp )
     call RECP_WRITE_TYPE(itypea, itypeb, INT_TYPE, lblop,PRINT_LEVEL)

    DO I = 1, N_INTEGRAL(lblop) 
       IF (RECP_DBG.GE.1) THEN
          IF (lblop.EQ.1) READ(74,'(x,f30.16)') values_stvc(I)
          IF (lblop.EQ.2) READ(75,'(x,f30.16)') values_stvc(I)
          IF (lblop.EQ.3) READ(76,'(x,f30.16)') values_stvc(I)
          IF (lblop.EQ.1) READ(77,'(x,i10)') value_temp 
          IF (lblop.EQ.2) READ(78,'(x,i10)') value_temp 
          IF (lblop.EQ.3) READ(79,'(x,i10)') value_temp
       ELSE
          IF (lblop.EQ.1) READ(74,*) values_stvc(I)
          IF (lblop.EQ.2) READ(75,*) values_stvc(I)
          IF (lblop.EQ.3) READ(76,*) values_stvc(I)
          IF (lblop.EQ.1) READ(77,*) value_temp 
          IF (lblop.EQ.2) READ(78,*) value_temp 
          IF (lblop.EQ.3) READ(79,*) value_temp
       ENDIF
       labels_stvc(1,I) =     IPT_IL(value_temp)/ilfact
       labels_stvc(2,I) = mod(IPT_IL(value_temp),ilfact)
!      print *,'labels',lblop,labels_stvc(1,I),labels_stvc(2,I),I
    ENDDO
!   # write integral to file
    CALL RECP_WRITE_ARRAY( nst,INT_TYPE,n_bfn_sym,nbft,labels_stvc,mpru,array_stvc, &
                           array_so,values_stvc,N_INTEGRAL,MN_INTEGRAL,lblop)
    CALL RECP_WRITE_FILE(INT_TYPE,N_INTEGRAL,nbft,array_stvc,array_so,INT_AREP, &
                         MN_INTEGRAL,LBLOP)

    IF(PRINT_LEVEL.GE. 5)  &
       WRITE (RECP_OUT,'(4x,a,a)') chrtyp,' integrals were written in a file'
    IF(PRINT_LEVEL.GE.10) PRINT *,'   N_INTEGRAL(',LBLOP,')',N_INTEGRAL(LBLOP)

  ENDDO

  DEALLOCATE( array_stvc )
  DEALLOCATE( array_so )
  DEALLOCATE( values_stvc )
  DEALLOCATE( labels_stvc )

  IF (nnam.EQ.5) THEN
     IF (RECP_DBG.GE.1) THEN
        CLOSE (UNIT=74)
        CLOSE (UNIT=75)
        CLOSE (UNIT=76)
        CLOSE (UNIT=77)
        CLOSE (UNIT=78)
        CLOSE (UNIT=79)
     ELSE
        CLOSE (UNIT=74,STATUS='DELETE')
        CLOSE (UNIT=75,STATUS='DELETE')
        CLOSE (UNIT=76,STATUS='DELETE')
        CLOSE (UNIT=77,STATUS='DELETE')
        CLOSE (UNIT=78,STATUS='DELETE')
        CLOSE (UNIT=79,STATUS='DELETE')
     ENDIF
  ELSE
     IF (RECP_DBG.GE.1) THEN
        CLOSE (UNIT=74)
        CLOSE (UNIT=77)
     ELSE
        CLOSE (UNIT=74,STATUS='DELETE')
        CLOSE (UNIT=77,STATUS='DELETE')
     ENDIF
  ENDIF
  CALL OUTPUT_LOC('RECP_WRITE_MAIN','X')
END SUBROUTINE RECP_WRITE_MAIN

! ------------------------------------------------


SUBROUTINE RECP_WRITE_ARRAY(n_sym,INT_TYPE,n_bfn_sym,nbft,labels_stvc,mpru,array_stvc,array_so, & 
                             values_stvc,N_INTEGRAL,MN_INTEGRAL,lblop)
! From cnvrt_stvz, Y.C. Park
  USE RECP_FUNCTION1
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "inc_mxvalue.h"
#include "inc_print.h"
  integer :: n_bfn_sym(mstup) ! nbpsy or nso : number of basis functions per symmetry block
  integer :: INT_TYPE         ! integral type 
  integer :: n_sym            ! nsym or nst : number of symmetry blocks
  integer :: nbft
  integer :: mpru
  integer :: lblop
  integer :: N_INTEGRAL(3), MN_INTEGRAL
  real(8) :: array_stvc(:)   ! array_stvc(MN_INTEGRAL) 
  real(8) :: array_so(:,:)     ! array_so(nbft,nbft) 
  real(8) :: values_stvc(:)  ! values_stvc(MN_INTEGRAL)
  integer :: labels_stvc(2,*) ! labels(1:2,1:n1max) orbital label buffer. 
! variables : local
  integer :: i, j, ij
  integer :: n_tot        ! n_total
  integer :: n_n_tot      ! n! of n_total
  integer :: n_skp(8)     ! n_skip
  integer :: n_n_skp(8)   ! n! of n_skip
  integer :: i_lab, j_lab
  integer :: SYM1, SYM2, SYM_12   ! isym, jsym, ijsym
  integer, ALLOCATABLE :: sym_bfn(:) ! symb(1:nbft) symmetry index of each basis function 'n_bfn_sym' is needed
  integer :: nndxf
  nndxf(i) = (i*(i-1))/2     ! define nndxf

  CALL OUTPUT_LOC('RECP_WRITE_ARRAY','E')
  IF (PRINT_LEVEL .GE. 5) THEN 
     WRITE(RECP_OUT,*)'   INT_TYPE       =',INT_TYPE
     WRITE(RECP_OUT,*)'   N_INTEGRAL     =',N_INTEGRAL(LBLOP)
     WRITE(RECP_OUT,*)'   N_SYM          =',N_SYM
     WRITE(RECP_OUT,*)'   MN_INTEGRAL    =',MN_INTEGRAL
  ENDIF

! Check and allocate sym_bfn  
  n_tot  = 0
  DO SYM1 = 1, n_sym 
     n_tot = n_tot + n_bfn_sym(SYM1)               
  ENDDO
  ALLOCATE ( sym_bfn(n_tot) )  !nbfmx(100000)->n_tot
  CALL RECP_SETZERO_I1(sym_bfn,n_tot)

! make sym_bfn 
  n_tot  = 0
  n_n_tot = 0
  DO SYM1 = 1, n_sym
     n_skp(SYM1)   = n_tot  
     n_n_skp(SYM1) = n_n_tot 
     n_tot         = n_tot   + n_bfn_sym(SYM1)               
     n_n_tot       = n_n_tot + nndxf(n_bfn_sym(SYM1) + 1)  
!    print *, '1-n_bfn_sym :',SYM1,n_bfn_sym(SYM1)
     DO j = (n_skp(SYM1)+1), n_tot
        sym_bfn(j)  =SYM1
     ENDDO
  ENDDO

  IF (INT_TYPE.LE.4) THEN
! # symmetric matrix, totally symmetric operator.
!    #  initialize matrix element 
     array_stvc = 0

     do i = 1, N_INTEGRAL(1)
       i_lab  = max( labels_stvc(1,i), labels_stvc(2,i) )
       j_lab  = min( labels_stvc(1,i), labels_stvc(2,i) )
       SYM1  = sym_bfn(i_lab)
       SYM2  = sym_bfn(j_lab)
       SYM_12 = nndxf(SYM1) + SYM2
       if (SYM1.eq.SYM2) then
!         # only valid array(*) elements are referenced.
          ij  =  n_n_skp(SYM1) + nndxf(i_lab - n_skp(SYM1)) + j_lab - n_skp(SYM1)
!       | ij  :  skip for       +                             +       -              |
!       |        each symmetry                                                       | 
!         print *,'ij',ij,i_lab,j_lab,SYM1,SYM2
          array_stvc(ij) = values_stvc(i)
       endif
     enddo

  elseif ( INT_TYPE .ge. 5 ) then
! # antisymmetric matrix, nontotally symmetric operator.
! # initialize spin-orbit matrix element 
     array_so = 0
     do i = 1, N_INTEGRAL(LBLOP)
        i_lab  = labels_stvc(1,i)
        j_lab  = labels_stvc(2,i)
        if (i_lab .lt. j_lab) then
           i_lab = labels_stvc(2,i)
           j_lab = labels_stvc(1,i)
!          the sign change upon transposition.
           values_stvc(i) = -values_stvc(i)
        endif
        array_so(i_lab,j_lab) =  values_stvc(i)
        array_so(j_lab,i_lab) = -values_stvc(i)
!       print *,'array_so',i_lab,j_lab,array_so(i_lab,j_lab)
     enddo
  endif

  DEALLOCATE ( sym_bfn )

  CALL OUTPUT_LOC('RECP_WRITE_ARRAY','X')
END SUBROUTINE RECP_WRITE_ARRAY


! ------------------------------------------------

SUBROUTINE RECP_WRITE_FILE(INT_TYPE,N_INTEGRAL,nbft,array_stvc,array_so, &
                            INT_AREP,MN_INTEGRAL,LBLOP)
! write ro file
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER :: LBLOP
  integer :: i,j
  integer :: INT_TYPE, INT_AREP
  integer :: N_INTEGRAL(3), MN_INTEGRAL
  integer :: nbft
  real(8) :: array_stvc(:)  ! array_stvc(MN_INTEGRAL)
  real(8) :: array_so(:,:)    ! array_so(nbft,nbft)

  LOGICAL :: FILE_EXIST  = .FALSE.
  LOGICAL :: ARGOS_DEBUG = .FALSE.
  INTEGER :: FILE_NUM(7)
  CHARACTER(10) :: FILE_INT(7),FILE_DBG(7),FILE_NAME(7)
  DATA FILE_NUM / 64, 64, 64, 64, 65, 66, 67 /
  DATA FILE_INT /'RECP_INT_S','RECP_INT_T','RECP_INT_V', &
                 'RECP_INT_C','RECP_INT_X','RECP_INT_Y', &
                 'RECP_INT_Z'/
  DATA FILE_DBG /'RECP_DBG_S','RECP_DBG_T','RECP_DBG_V', &
                 'RECP_DBG_C','RECP_DBG_X','RECP_DBG_Y', &
                 'RECP_DBG_Z'/
  DATA FILE_NAME/'Overlap   ','Kinetic   ','Potential ', &
                 'RECP AVE  ','RECP SO(X)','RECP SO(Y)', &
                 'RECP SO(Z)'/

  CALL OUTPUT_LOC('RECP_WRITE_FILE','E')

! # Skip writing S,T,V integral
  IF(INT_TYPE .LE. 3) GOTO 100

! # Delete if the file exist (for geometry optimization)
  INQUIRE(FILE=FILE_INT(INT_TYPE),EXIST=FILE_EXIST)
  IF (FILE_EXIST) THEN
!    WRITE(6,5X,A,A,A) 'File ',FILE_INT(INT_TYPE),' will be updated'
     OPEN (UNIT=FILE_NUM(INT_TYPE),FILE=FILE_INT(INT_TYPE),FORM='UNFORMATTED',STATUS='OLD')
     CLOSE(UNIT=FILE_NUM(INT_TYPE),STATUS='DELETE')
  ENDIF

! # Write RECP integrals to file 
  OPEN(UNIT=FILE_NUM(INT_TYPE),FILE=FILE_INT(INT_TYPE),FORM='UNFORMATTED')
  IF (INT_TYPE.EQ.4) THEN
     WRITE(FILE_NUM(INT_TYPE)) (ARRAY_STVC(I),I=1,N_INTEGRAL(LBLOP)) !yc-check
     WRITE(RECP_OUT,'(A,A)')   ' * VCORE was generated : ', FILE_INT(INT_TYPE)
  ENDIF
  IF (INT_TYPE.GE.5) WRITE(FILE_NUM(INT_TYPE)) ((ARRAY_SO(J,I),J=1,NBFT),I=1,NBFT)
     IF (INT_TYPE.EQ.5) THEN
        WRITE(RECP_OUT,'(A,A)')' * SOX was generated   : ', FILE_INT(INT_TYPE)
     ELSEIF (INT_TYPE.EQ.6) THEN
        WRITE(RECP_OUT,'(A,A)')' * SOY was generated   : ', FILE_INT(INT_TYPE)
     ELSEIF (INT_TYPE.EQ.7) THEN
        WRITE(RECP_OUT,'(A,A)')' * SOZ was generated   : ', FILE_INT(INT_TYPE)
     ENDIF
  CLOSE(UNIT=FILE_NUM(INT_TYPE))

! --------------------------------------
! Write SOREP as 0 when AREP calculation
! --------------------------------------
  IF(INT_TYPE.EQ.4 .AND. INT_AREP.eq.1) THEN 
     ARRAY_SO = 0
     DO INT_TYPE = 5, 7  !INT_TYPE = 5~7
!       Delete if the file exist (for geometry optimization)
        INQUIRE(FILE=FILE_INT(INT_TYPE),EXIST=FILE_EXIST)
        IF(FILE_EXIST) THEN
           OPEN (UNIT=FILE_NUM(INT_TYPE),FILE=FILE_INT(INT_TYPE), FORM='UNFORMATTED',STATUS='OLD')
           CLOSE(UNIT=FILE_NUM(INT_TYPE),STATUS='DELETE')
        ENDIF
!       Make 'ARRAY_SO = 0'
        OPEN(UNIT=FILE_NUM(INT_TYPE),FILE=FILE_INT(INT_TYPE), FORM='UNFORMATTED')
        WRITE(FILE_NUM(INT_TYPE)) ((ARRAY_SO(J,I),J=1,NBFT),I=1,NBFT) 
        CLOSE(UNIT=FILE_NUM(INT_TYPE))
     ENDDO
     INT_TYPE = 4
  ENDIF

! ----------------------------------------
! Write RECP integrals to file (for debug)
! ----------------------------------------
  IF (RECP_DBG.GE.5) THEN
     OPEN (UNIT=FILE_NUM(INT_TYPE),FILE=FILE_DBG(INT_TYPE), FORM='FORMATTED')
     IF (INT_TYPE.EQ.4) THEN
        WRITE(FILE_NUM(INT_TYPE),'(a)') FILE_NAME(INT_TYPE) 
        DO I = 1,N_INTEGRAL(LBLOP)
           WRITE(FILE_NUM(INT_TYPE),'(i7,2x,f20.16)') I,ARRAY_STVC(I)
        ENDDO
     ELSEIF (INT_TYPE.GE.5) THEN
        WRITE(FILE_NUM(INT_TYPE),'(a)') FILE_NAME(INT_TYPE) 
        DO I = 1,NBFT
           DO J = 1,NBFT
              WRITE(FILE_NUM(INT_TYPE),'(i7,i7,2x,f20.16)') I,J,ARRAY_SO(J,I)
           ENDDO
        ENDDO
     ENDIF
     CLOSE(UNIT=FILE_NUM(INT_TYPE))
  ENDIF
 
  100  CONTINUE

  CALL OUTPUT_LOC('RECP_WRITE_FILE','X')

END SUBROUTINE RECP_WRITE_FILE


! ------------------------------------------------


SUBROUTINE RECP_WRITE_TYPE(itypea, itypeb, INT_TYPE, lblop,PRINT_LEVEL)
! Check integral type 
! INT_TYPE = 1-4 : symmetric matrix, totally symmetric operator 
! INT_TYPE =     : symmetric matrix, nontotally symmetric operator 
! INT_TYPE = 5-7 : antisymmetric matrix, nontotally symmetric operator 
  IMPLICIT NONE
  INTEGER :: PRINT_LEVEL
  INTEGER :: itypea,  itypeb, INT_TYPE, lblop

  if (itypea.eq.0) then
     if     (itypeb.eq.0) then
        INT_TYPE = 1        !   S
     elseif (itypeb.eq.1) then
        INT_TYPE = 2        !   T 
     elseif (itypeb.eq.2) then
        INT_TYPE = 3        !   V
     elseif (itypeb.eq.3) then
        INT_TYPE = 4        !   Core or Veff
     endif
  elseif (itypea.eq.1) then
     print *, 'This integral type is not supported for the now'
     stop
  elseif (itypea.eq.2) then
     if     (lblop.eq.1) then
        INT_TYPE = 5         !   SO(x)
     elseif (lblop.eq.2) then
        INT_TYPE = 6         !   SO(y)
     elseif (lblop.eq.3) then
        INT_TYPE = 7         !   SO(z)
     endif
  endif
END SUBROUTINE RECP_WRITE_TYPE


! ------------------------------------------------


SUBROUTINE RECP_WRITE_FILEOPEN(UNIT0,FILE0,STAT0)
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER :: UNIT0,IOS0
  CHARACTER(LEN=*) :: FILE0 
  CHARACTER(LEN=*) :: STAT0 
  IF (RECP_DBG.GE.1) THEN
     OPEN (UNIT=UNIT0,FILE=FILE0,IOSTAT=IOS0,STATUS=STAT0,FORM='FORMATTED')
     WRITE(RECP_OUT,'(A,A,A)') ' * File (',FILE0,') is opened.'  
  ELSE
     OPEN (UNIT=UNIT0,FILE=FILE0,IOSTAT=IOS0,STATUS=STAT0)
  ENDIF
! Error message
  IF (IOS0.NE.0) THEN
     WRITE(RECP_OUT,'(/,A,A)') ' * Error in opening ',FILE0
     CALL QUIT('Error in opening a file')
  ENDIF
END SUBROUTINE RECP_WRITE_FILEOPEN

END MODULE RECP_WRITE
