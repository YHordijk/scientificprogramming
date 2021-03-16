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

***********************************************************************

      SUBROUTINE COPNBLKD(LUIN,LUOUT,SEGMNT,NBLOCK,IREW,LBLK)
*
* If (LUOUT .gt. 0) Copy Nblocks from FILE LUIN to file LUOUT
*
* If (LUOUT .le. 0) Skip Nblocks on FILE LUIN
*
*
* Type of file LUOUT is inherited from LUIN
*
* If less than NBLOCK blocks are left in the given vector
* on LUIN a warning is printed,
* but normal return
*
*
* NO EOV mark written
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        REWIND  LUIN
        if (LUOUT .gt. 0) REWIND  LUOUT
      END IF

      DO IIBLK = 1, NBLOCK
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          READ(LUIN) LBL
          if (LUOUT .gt. 0 .AND. LBL .ge. 0) WRITE(LUOUT) LBL
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          IF (LUOUT .gt. 0) THEN
            CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
            IF(IAMPACK.EQ.0) THEN
              CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            CALL SKPRCD2(LBL,-1,LUIN)
          END IF
        END IF
        IF(LBL.EQ.-1) THEN
           WRITE(6,*) ' END OF VECTOR mark for BLOCK'
           write(6,*) 'IBLOCK not defined!!'
           call quit('Stop in COPNBLKD')
           GOTO 3333
        END IF
      END DO
*     ^End of loop over blocks
 3333 CONTINUE
*
      END
***********************************************************************

      SUBROUTINE COPVCD_REL(LUIN,LUOUT,SEGMNT,IREW,LBLK)
C
C COPY VECTOR ON FILE LUIN TO FILE LUOUT
C
C
C LBLK DEFINES STRUCTURE OF FILE
*
* Structure of output file is inherited by output file,
* if input file is packed, so is output file
*
*
C Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
      integer iscr(2)
      integer LUIN,LUOUT, myluout
C
      IF( IREW .NE. 0 ) THEN
        REWIND(LUIN)
        REWIND(LUOUT)
      END IF
 
      call izero(iscr,2)
C
C LOOP OVER BLOCKS
C
        myluout = LUOUT
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(myluout) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS( LBL,1,-1,myluout)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &         NO_ZEROING)
          IF(IAMPACK.NE.0) THEN
            WRITE(6,*) ' COPVCD_REL, IMZERO,IAMPACK,FILE = ', 
     &      IMZERO,IAMPACK,LUIN
          END IF
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              CALL TODSC(SEGMNT,LBL,KBLK,myluout)
            ELSE  
              CALL TODSCP(SEGMNT,LBL,KBLK,myluout)
            END IF
          ELSE
            ISCR(1) = 1
            ISCR(2) = 0
            CALL ITODS(ISCR,2,2,myluout)
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      END
***********************************************************************

      SUBROUTINE COPVCDC(LUIN,LUOUT,VEC1,IREW,IRI,LBLK)
*
* Copy complex vector from LUIN to LUOUT
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION VEC1(*)
*
      IF(IREW.NE.0) THEN
        CALL REWINE(LUIN,LBLK)
        CALL REWINE(LUOUT,LBLK)
      END IF
*
*. Real part
      CALL COPVCD_REL(LUIN,LUOUT,VEC1,0,LBLK)
*
      if (IRI.eq.2) then
*. Imaginary part
        CALL COPVCD_REL(LUIN,LUOUT,VEC1,0,LBLK)
      end if
*
      END
!**********************************************************************

      subroutine krci_cp_sel_root(xflabel,istate,iri,omega_val)
!**********************************************************************
      implicit none

      character (len=3),  intent(in)     :: xflabel
      integer,            intent(in)     :: istate
      integer,            intent(in)     :: iri
      integer,            intent(in)     :: omega_val
!----------------------------------------------------------------------
      real(8), allocatable               :: vecblock(:)
      integer                            :: in_fh
      integer                            :: out_fh
      integer                            :: block_length
      integer                            :: iszero
      integer                            :: ispacked
      integer                            :: iskip
      integer                            :: ipart
      integer                            :: loutf
      character (len=14)                 :: infile_name
      character (len=21)                 :: outfile_name
      logical                            :: is_opened
!**********************************************************************
!
!     initialize scratch variables
      in_fh  = 98
      out_fh = 99
      block_length = 0

      loutf = -1
!     initialize file names
      write(infile_name,'(a11,a3)')     'KRCI_CVECS.',xflabel
      select case(len_trim(xflabel))
        case (1)
          write(outfile_name,'(a11,a1,a7)') 
     &    'KRCI_CVECS.',xflabel,'.select'
          loutf = 19
        case (2)
          write(outfile_name,'(a11,a2,a7)') 
     &    'KRCI_CVECS.',xflabel,'.select'
          loutf = 20
        case (3)
          write(outfile_name,'(a11,a3,a7)') 
     &    'KRCI_CVECS.',xflabel,'.select'
          loutf = 21
        case default
          call quit('non-existing file label in krci_cp_sel_root')
      end select

!     step 1: open files
      inquire(file=infile_name, opened=is_opened)
      if(.not.is_opened)then
        open(in_fh,file=infile_name,status='old',form='unformatted',
     &       position='rewind')
      else
        inquire(file=infile_name, number=in_fh)
        rewind in_fh  
      end if
      open(out_fh,file=outfile_name,status='new',form='unformatted',
     &     position='rewind')

!     step 2: forward to the right position in the vector file
      do iskip = 1, ((istate-1)*iri)
 101    continue
          read(in_fh) block_length
          if(block_length.ge.0)then
            allocate(vecblock(block_length))
            call frmdsc(vecblock,block_length,-1,in_fh,iszero,ispacked)
            deallocate(vecblock)
          end if
        if(block_length.ge.0) goto 101
      end do

!     step 3: copy the selected vector
!     real part + imaginary part (iri == 2)     
      do ipart = 1, iri
 202    continue
        call ifrmds(block_length,1,-1,in_fh)
        call itods(block_length,1,-1,out_fh)
        if(block_length.ge.0)then
          allocate(vecblock(block_length))
          call frmdsc2(vecblock,block_length,-1,in_fh,iszero,ispacked,1)
          if(iszero.eq.0)then
            if(ispacked.eq.0)then
              call todsc(vecblock,block_length,-1,out_fh)
            else
              call todscp(vecblock,block_length,-1,out_fh)
            end if
          else
            call zerorc_rel(block_length,out_fh)
          end if
          deallocate(vecblock)
        end if
        if(block_length.ge.0) goto 202
      end do

!     step 4: close files
      if(.not.is_opened) close(in_fh)
      close(out_fh)

      print '(/2x,a,i3,a,i2,a,a21)', 
     &      ' *** selected CI vector #',istate,'(omega = ',omega_val,
     &      '[/2 for fermionic systems]) saved on file ',outfile_name
      end
***********************************************************************

      SUBROUTINE DISKUN_REL
      IMPLICIT REAL*8           (A-H,O-Z)
#include "ipoist8.inc"
*
* Assign logical unit numbers for LUCIA:
*
*==================================================
*  Complete revision motivated by MCSCF development
*==================================================
*
*   1 -  9 :  Retained for I/O
*  ---------
*        5 :  Standard IN
*        6 :  Standard OUT
*        7 :  Standard OUT (Slaves)
*
*  60 - 69 :  Important files retained (DIAG,C,HC etc.)
*  ---------
*       60 :    CI diagonal
*       61 :    CI reference vector
*       62 :    CI sigma vector
*
*  70 - 79 :  Scratch files
*  ---------
*
*  80 - 99 :  Files in use in different context 
*             (e.g. standalone, INS, MOLCAS etc.)
*  ---------
*
*  Timo Fleig    Mar 19, 2002
*
*  changes for parallel business: Stefan Knecht   Apr 24, 2007
*
*==============================================
*
*  Old setup:
*
* All files with some kind of input information  :  10 - 19
* Internal files (retained through job)          :  20 - 29
* Scratch files                                  :  30 - 40
* Files for machine information			 :  50 - 59
* Additional scratch files		       	 :  41,61,62
* Files for relativistic business (SO ints etc.) :  70 - 79
* All files containing final results             :  90 - 99
*
*
#include "clunit.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#endif
#include "mxpdim.inc"
#include "parluci.h"
#include "files_r.inc"
#include "cstate.inc"
#include "krciprop.h"
      INTEGER LUSCR_MPIL,LUDIA_MPIL,LUC_MPIL, LUHC_MPIL, LUSC123_MPIL
      INTEGER LUSC344_MPIL
      LOGICAL EX
*
* =========================
* Standard input and output
* =========================
*. Input file
      LUIN = 5
      LU5 = LUIN
*. Output file
      LUOUT = 6
* =================
* Input information
* =================
* Input file containing MO-AO transformation matrix
      LUMOIN = - 1
*. Input file for CI-vectors
*. restart from file 21 is assumed
*. Input , two electron integrals
      LU2INT = 81
      LU40 = LU2INT
*. Input , one electron integrals
      LU1SO = 82
      LU1INT = 83
      LU2 = LU1INT        ! for MOLCAS environment
      LU9 = LU1INT        ! for MOLECULE environment
*. Input , property one-electron integral files
      LUPRP  = 84
*. Sirius interface file
      LUSIR1 = 85
*. Input MOAO coefficients
      LU19 = 86
*. Input Spin-orbit integrals (X)
      LU69 = 87
*. Input Spin-orbit integrals (X)
      LU70 = 88
*. Input Spin-orbit integrals (X)
      LU71 = 89
* =================
* Internal files
* =================
*. CI diagonal
      LUDIA = 60
      LUDIA_MPIN ="CIDIA_REL"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUDIA_MPI,'(A9,A1,I1)') LUDIA_MPIN,'.',MYPROC
         LUDIA_MPIL=11
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUDIA_MPI,'(A9,A1,I2)') LUDIA_MPIN,'.',MYPROC
         LUDIA_MPIL=12
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUDIA_MPI,'(A9,A1,I3)') LUDIA_MPIN,'.',MYPROC
         LUDIA_MPIL=13
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUDIA_MPI=LUDIA_MPIN
      LUDIA_MPIL=9
#endif
*
      Open(Unit=LUDIA,File=NLUDIA_MPI(1:LUDIA_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
*. CI vector
      LUC =   61
      LUC_MPIN ="KRCI_CVECS"
      LUC_MPIL=14
      IF( MYPROC .eq. MASTER )THEN
        IF( IREFSM .le. NMAX_SYM_CI )THEN
          SYMFLABEL = XSYMFLAB(IREFSM)
C         WRITE(6,*) ' SYMFLABEL is now ',SYMFLABEL
        ELSE 
          WRITE(LUWRT,*) '  *** Error in DISKUN_REL *** assignment of'//
     &                   ' symmetry label to KRCI_CVECS file failed'   
          CALL QUIT("*** Error in DISKUN_REL *** assignment of 
     &              symmetry label to KRCI_CVECS file failed ")
        END IF
        IF( .not. RUNXPROP )THEN
           WRITE (NLUC_MPI,'(A10,A1,A3)') LUC_MPIN,'.',SYMFLABEL
           Open(Unit=LUC,File=NLUC_MPI(1:LUC_MPIL),Status='UNKNOWN',
     &          Form='UNFORMATTED')
        END IF
      END IF
*. Sigma vector file
      LUHC =  62
      LUHC_MPIN ="HCFILE_REL"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUHC_MPI,'(A10,A1,I1)') LUHC_MPIN,'.',MYPROC
         LUHC_MPIL=12
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUHC_MPI,'(A10,A1,I2)') LUHC_MPIN,'.',MYPROC
         LUHC_MPIL=13
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUHC_MPI,'(A10,A1,I3)') LUHC_MPIN,'.',MYPROC
         LUHC_MPIL=14
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUHC_MPI=LUHC_MPIN
      LUHC_MPIL=10
#endif
      Open(Unit=LUHC,File=NLUHC_MPI(1:LUHC_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
* =================
* Scratch files
* =================
      LUSC1 = 70
      LUSC1_MPIN ="LUSC1_REL"
      LUSC2 = 71
      LUSC2_MPIN ="LUSC2_REL"
      LUSC3 = 72
      LUSC3_MPIN ="LUSC3_REL"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC1_MPI,'(A9,A1,I1)') LUSC1_MPIN,'.',MYPROC
         WRITE (NLUSC2_MPI,'(A9,A1,I1)') LUSC2_MPIN,'.',MYPROC
         WRITE (NLUSC3_MPI,'(A9,A1,I1)') LUSC3_MPIN,'.',MYPROC
         LUSC123_MPIL=11
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC1_MPI,'(A9,A1,I2)') LUSC1_MPIN,'.',MYPROC
         WRITE (NLUSC2_MPI,'(A9,A1,I2)') LUSC2_MPIN,'.',MYPROC
         WRITE (NLUSC3_MPI,'(A9,A1,I2)') LUSC3_MPIN,'.',MYPROC
         LUSC123_MPIL=12
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC1_MPI,'(A9,A1,I3)') LUSC1_MPIN,'.',MYPROC
         WRITE (NLUSC2_MPI,'(A9,A1,I3)') LUSC2_MPIN,'.',MYPROC
         WRITE (NLUSC3_MPI,'(A9,A1,I3)') LUSC3_MPIN,'.',MYPROC
         LUSC123_MPIL=13
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUSC1_MPI=LUSC1_MPIN
      NLUSC2_MPI=LUSC2_MPIN
      NLUSC3_MPI=LUSC3_MPIN
      LUSC123_MPIL=9
#endif
      Open(Unit=LUSC1,File=NLUSC1_MPI(1:LUSC123_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC2,File=NLUSC2_MPI(1:LUSC123_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC3,File=NLUSC3_MPI(1:LUSC123_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
*. Scratch space for subspace handling
      LUSC34 = 73
      LUSC34_MPIN ="LUSC34_REL"
      LUSC35 = 74
      LUSC35_MPIN ="LUSC35_REL"
      LUSC36 = 75
      LUSC36_MPIN ="LUSC36_REL"
      LUSC37 = 76
      LUSC37_MPIN ="LUSC37_REL"
      LUSC38 = 77
      LUSC38_MPIN ="LUSC38_REL"
      LUSC39 = 78
      LUSC39_MPIN ="LUSC39_REL"
      LUSC41 = 79
      LUSC41_MPIN ="LUSC41_REL"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC34_MPI,'(A10,A1,I1)') LUSC34_MPIN,'.',MYPROC
         WRITE (NLUSC35_MPI,'(A10,A1,I1)') LUSC35_MPIN,'.',MYPROC
         WRITE (NLUSC36_MPI,'(A10,A1,I1)') LUSC36_MPIN,'.',MYPROC
         WRITE (NLUSC37_MPI,'(A10,A1,I1)') LUSC37_MPIN,'.',MYPROC
         WRITE (NLUSC38_MPI,'(A10,A1,I1)') LUSC38_MPIN,'.',MYPROC
         WRITE (NLUSC39_MPI,'(A10,A1,I1)') LUSC39_MPIN,'.',MYPROC
         WRITE (NLUSC41_MPI,'(A10,A1,I1)') LUSC41_MPIN,'.',MYPROC
         LUSC344_MPIL=12
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC34_MPI,'(A10,A1,I2)') LUSC34_MPIN,'.',MYPROC
         WRITE (NLUSC35_MPI,'(A10,A1,I2)') LUSC35_MPIN,'.',MYPROC
         WRITE (NLUSC36_MPI,'(A10,A1,I2)') LUSC36_MPIN,'.',MYPROC
         WRITE (NLUSC37_MPI,'(A10,A1,I2)') LUSC37_MPIN,'.',MYPROC
         WRITE (NLUSC38_MPI,'(A10,A1,I2)') LUSC38_MPIN,'.',MYPROC
         WRITE (NLUSC39_MPI,'(A10,A1,I2)') LUSC39_MPIN,'.',MYPROC
         WRITE (NLUSC41_MPI,'(A10,A1,I2)') LUSC41_MPIN,'.',MYPROC
         LUSC344_MPIL=13
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC34_MPI,'(A10,A1,I3)') LUSC34_MPIN,'.',MYPROC
         WRITE (NLUSC35_MPI,'(A10,A1,I3)') LUSC35_MPIN,'.',MYPROC
         WRITE (NLUSC36_MPI,'(A10,A1,I3)') LUSC36_MPIN,'.',MYPROC
         WRITE (NLUSC37_MPI,'(A10,A1,I3)') LUSC37_MPIN,'.',MYPROC
         WRITE (NLUSC38_MPI,'(A10,A1,I3)') LUSC38_MPIN,'.',MYPROC
         WRITE (NLUSC39_MPI,'(A10,A1,I3)') LUSC39_MPIN,'.',MYPROC
         WRITE (NLUSC41_MPI,'(A10,A1,I3)') LUSC41_MPIN,'.',MYPROC
         LUSC344_MPIL=14
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUSC34_MPI=LUSC34_MPIN
      NLUSC35_MPI=LUSC35_MPIN
      NLUSC36_MPI=LUSC36_MPIN
      NLUSC37_MPI=LUSC37_MPIN
      NLUSC38_MPI=LUSC38_MPIN
      NLUSC39_MPI=LUSC39_MPIN
      NLUSC41_MPI=LUSC41_MPIN
      LUSC344_MPIL=10
#endif
      Open(Unit=LUSC34,File=NLUSC34_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC35,File=NLUSC35_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC36,File=NLUSC36_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC37,File=NLUSC37_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC38,File=NLUSC38_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC39,File=NLUSC39_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC41,File=NLUSC41_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
*
*. Reference vectors for (transition) densities
*
      LUSC61 = 65
C     LUSC62 = 66
C  avoid unit 66. Used in DALTON
CSK      LUSC62 = 95    ! conflict might occur with LU18 in genh1.F
      LUSC62 = 67    ! might this work? SK April 24, 2007
      LUSC61_MPIN ="LUSC61_REL"
      LUSC62_MPIN ="LUSC62_REL"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC61_MPI,'(A10,A1,I1)') LUSC61_MPIN,'.',MYPROC
         WRITE (NLUSC62_MPI,'(A10,A1,I1)') LUSC62_MPIN,'.',MYPROC
         LUSC344_MPIL=12
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC61_MPI,'(A10,A1,I2)') LUSC61_MPIN,'.',MYPROC
         WRITE (NLUSC62_MPI,'(A10,A1,I2)') LUSC62_MPIN,'.',MYPROC
         LUSC344_MPIL=13
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC61_MPI,'(A10,A1,I3)') LUSC61_MPIN,'.',MYPROC
         WRITE (NLUSC62_MPI,'(A10,A1,I3)') LUSC62_MPIN,'.',MYPROC
         LUSC344_MPIL=14
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUSC61_MPI=LUSC61_MPIN
      NLUSC62_MPI=LUSC62_MPIN
      LUSC344_MPIL=10
#endif
      Open(Unit=LUSC61,File=NLUSC61_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN', Form='UNFORMATTED')
      Open(Unit=LUSC62,File=NLUSC62_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN', Form='UNFORMATTED')
*. Some more general purpose scratch
      LUSC92 = 92
      LUSC93 = 93
      LUSC92_MPIN ="LUSC92_REL"
      LUSC93_MPIN ="LUSC93_REL"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC92_MPI,'(A10,A1,I1)') LUSC92_MPIN,'.',MYPROC
         WRITE (NLUSC93_MPI,'(A10,A1,I1)') LUSC93_MPIN,'.',MYPROC
         LUSC929_MPIL=12
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC92_MPI,'(A10,A1,I2)') LUSC92_MPIN,'.',MYPROC
         WRITE (NLUSC93_MPI,'(A10,A1,I2)') LUSC93_MPIN,'.',MYPROC
         LUSC929_MPIL=13
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC92_MPI,'(A10,A1,I3)') LUSC92_MPIN,'.',MYPROC
         WRITE (NLUSC93_MPI,'(A10,A1,I3)') LUSC93_MPIN,'.',MYPROC
         LUSC929_MPIL=14
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUSC92_MPI=LUSC92_MPIN
      NLUSC93_MPI=LUSC93_MPIN
      LUSC929_MPIL=10
#endif
*
      LUSC94 = 94
      LUSC95 = 95
      LUSC94_MPIN ="INTEGRALS"
      LUSC95_MPIN ="INDICES"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC94_MPI,'(A9,A1,I1)') LUSC94_MPIN,'.',MYPROC
         WRITE (NLUSC95_MPI,'(A7,A1,I1)') LUSC95_MPIN,'.',MYPROC
         LUSC94_MPIL=11
         LUSC95_MPIL=9
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC94_MPI,'(A9,A1,I2)') LUSC94_MPIN,'.',MYPROC
         WRITE (NLUSC95_MPI,'(A7,A1,I2)') LUSC95_MPIN,'.',MYPROC
         LUSC94_MPIL=12
         LUSC95_MPIL=10
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC94_MPI,'(A9,A1,I3)') LUSC94_MPIN,'.',MYPROC
         WRITE (NLUSC95_MPI,'(A7,A1,I3)') LUSC95_MPIN,'.',MYPROC
         LUSC94_MPIL=13
         LUSC95_MPIL=11
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUSC94_MPI=LUSC94_MPIN
      NLUSC95_MPI=LUSC95_MPIN
      LUSC94_MPIL=9
      LUSC95_MPIL=7
#endif
*     open and close of LUSC92 - LUSC95 at relevant places ...
      LUSC96 = 96     ! File MRCONEE
      LUSC97 = 97     ! File MDCINT
      LUSC98 = 98
C
C     file for reordered GAS-specific integrals
C
      IIJKL_ROD    = 99
      NIIJKL_ROD   = "IJKL_REOD"
      IIJKL_ROD_LL = 9
C           
C     =================
C       Output files
C     =================
C
C     output file containing energies for (converged) roots 
C     in each symmetry - file handle stored on common block 
C     in file parluci.h
      LU_INFO  = 80
      IF( MYPROC .eq. MASTER )THEN
C
C       If old 'KRCI_CVECS.INFO' exists, then we are restarting a KRCI.
        INQUIRE ( FILE = 'KRCI_CVECS.INFO', EXIST = EX )
        IF( EX ) THEN
          CALL OPNFIL(LU_INFO,'KRCI_CVECS.INFO','OLD','DISKUN')
        ELSE
          CALL OPNFIL(LU_INFO,'KRCI_CVECS.INFO','UNKNOWN','DISKUN')
          CALL NEWLAB('SODLABEL',LU_INFO,LUOUT)
        END IF
        REWIND LU_INFO
      END IF
C
C     file containing the property integrals (KR-CI/KR-MCSCF: KRMC_FOCK)
      LUPROP_MAT = 53
C 
C     natural spinors in terms of double group adapted atomic spin-orbitals
      LUMOUTR = 90
      LUMOUTI = 91
      LUMOUTR_MPIN ="NATSPIR" 
      LUMOUTI_MPIN ="NATSPII" 
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUMOUTR_MPI,'(A7,A1,I1)') LUMOUTR_MPIN,'.',MYPROC
         WRITE (NLUMOUTI_MPI,'(A7,A1,I1)') LUMOUTI_MPIN,'.',MYPROC
         LUMOUTRI_MPIL=9
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUMOUTR_MPI,'(A7,A1,I2)') LUMOUTR_MPIN,'.',MYPROC
         WRITE (NLUMOUTI_MPI,'(A7,A1,I2)') LUMOUTI_MPIN,'.',MYPROC
         LUMOUTRI_MPIL=10
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUMOUTR_MPI,'(A7,A1,I3)') LUMOUTR_MPIN,'.',MYPROC
         WRITE (NLUMOUTI_MPI,'(A7,A1,I3)') LUMOUTI_MPIN,'.',MYPROC
         LUMOUTRI_MPIL=11
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUMOUTR_MPI = LUMOUTR_MPIN
      NLUMOUTI_MPI = LUMOUTI_MPIN
      LUMOUTRI_MPIL = 7
#endif
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE GET_LUC_FILE(LUC,IXSYM,EX)
**********************************************************************
C
C     open coefficient file (fh LUC) for symmetry irrep IXSYM
C
C     written by S. Knecht - Nov 2008
C
C     last revision:
C
C**********************************************************************
#include "implicit.h"
#include "mxpdim.inc"
#include "cstate.inc"
#include "krciprop.h"
      CHARACTER LUC_NAME*10, my_LUC_NAME*14, SYMFLABEL*3
      LOGICAL EX, open_already
      EX = .FALSE.
      LUC_NAME    = "KRCI_CVECS"
      LUC_NAMEL   = 14
      IF( IXSYM .le. NMAX_SYM_CI )THEN
         SYMFLABEL = XSYMFLAB(IXSYM)
C         WRITE(6,*) ' SYMFLABEL is now ',SYMFLABEL
      ELSE 
        WRITE(6,*) '  *** Error in GET_LUC_FILE *** assignment of'//
     &              ' symmetry label to KRCI_CVECS file failed'   
        CALL QUIT("*** Error in GET_LUC_FILE *** assignment of 
     &            symmetry label to GET_LUC_FILE file failed ")
      END IF
C     set file name and open file (if existing)
      WRITE (my_LUC_NAME,'(A10,A1,A3)') LUC_NAME,'.',SYMFLABEL
      INQUIRE (FILE=my_LUC_NAME,EXIST=EX,opened=open_already)
      if(open_already) print *,'*** attention: file already open***'
      IF( EX )THEN 
        OPEN(UNIT=LUC,FILE=my_LUC_NAME(1:LUC_NAMEL),STATUS='OLD',
     &       FORM='UNFORMATTED')
        REWIND LUC
      END IF
C
      END
