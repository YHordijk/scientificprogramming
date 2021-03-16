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

      SUBROUTINE GET_IAIB_DIRAC(IASP_KRCI_SYM,IBSP_KRCI_SYM,ELM_ABSMAX,
     &                          AMAT,NROW,NCOL,LRQ,LCQ,NA1,NA2,
     &                          JZ,IQP,NTEST)
C***********************************************************************
C
C     get symmetry indices of spinors which are connected by the 
C     largest absolute element ELM_ABSMAX in matrix AMAT(LROW,LCOL,JZ).
C
C     order in DIRAC (Molfdir format) matrix:
C                                 _
C           p spinors followed by p
C                                                 _     _
C           if NFSYM = 2 (inversion symmetry): p, p, q, q
C
C
C     Written by S. Knecht - Sep 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION AMAT(NROW,NCOL,JZ), IQP(*)
      logical :: indices_found

      indices_found = .false.
!
      DO 10 IZ = 1, JZ

#ifdef LUCI_DEBUG
        call header('property mat',-1)
        call prqmat(amat,NROW,NCOL,NROW,NCOL,iz,iqp,lupri)
#endif

        CALL GET_KRCI_SPINOR_IND(AMAT(1,1,IZ),1,NROW,1,NCOL,LRQ,LCQ,NA1,
     &                           NA2,1,ELM_ABSMAX,IASP_KRCI_SYM,
     &                           IBSP_KRCI_SYM,INDICES_FOUND,NTEST)
 10   CONTINUE
!
!     error exit if no indices pair is found
      IF(.NOT. INDICES_FOUND)THEN
        WRITE(LUPRI,'(/A)')'  *** ERROR in GET_IAIB_DIRAC ***'//
     &                     ' failed to find indices pair'//
     &                     ' of DIRAC spinors for largest element in'//
     &                     ' property matrix.'
        CALL QUIT('  *** ERROR in GET_IAIB_DIRAC *** failed to find 
     &             DIRAC spinor indices pair')

      END IF
C
      END
!***********************************************************************

      SUBROUTINE GET_KRCI_SPINOR_IND(AMATRX,ROWLOW,ROWHI,COLLOW,
     &                               COLHI,ROWDIM,COLDIM,NA1,NA2,NCTL,
     &                               ELM_ABSMAX,
     &                               IASP_KRCI_SYM,IBSP_KRCI_SYM,
     &                               INDICES_FOUND,NTEST)
!***********************************************************************
!
!     find spinor indices IASP_DC and IBSP_DC and determine the 
!     corresponding symmetry indices IASP_KRCI_SYM and IBSP_KRCI_SYM.
!
!     order in DIRAC (Molfdir format) matrix:
!                                 _
!           p spinors followed by p
!                                                 _     _
!           if NFSYM = 2 (inversion symmetry): p, p, q, q
!
!     Written by S. Knecht - Sep 2008
!
!     Last revision :
!
!***********************************************************************
      use symmetry_setup_krci
      use mospinor_info
      implicit none
#include "priunit.h"
!-------------------------------------------------------------------------------
      integer, intent(in)    :: rowlow, rowhi, collow, colhi 
      integer, intent(in)    :: rowdim, coldim
      integer, intent(in)    :: na1, na2
      integer, intent(in)    :: nctl
      real(8), intent(in)    :: amatrx(rowdim,coldim)
      real(8), intent(in)    :: elm_absmax
      integer, intent(inout) :: iasp_krci_sym, ibsp_krci_sym
      integer, intent(inout) :: ntest
      logical, intent(inout) :: indices_found
!-------------------------------------------------------------------------------
      integer, parameter     :: kcolp = 4, kcoln = 6
      integer                :: iasp_dc, ibsp_dc
      integer                :: irassign, icassign
      integer                :: kcol, last, begin, k, i
!-------------------------------------------------------------------------------
!#define LUCI_DEBUG
 
#ifdef LUCI_DEBUG
      ntest = 100 ! debug
#endif
      CALL QENTER('GET_KRCI_SPINOR_IND')
!
!     ... check if we already have what we want...
      IF (INDICES_FOUND) GOTO 999
 
      IASP_KRCI_SYM = 0
      IBSP_KRCI_SYM = 0
      IASP_DC       = 0
      IBSP_DC       = 0
 
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
 
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (AMATRX(K,I).NE.0.0d0) GO TO 5
    4       CONTINUE
         GO TO 1
    5       CONTINUE
            DO 7 I = BEGIN,LAST 
              IF( INDICES_FOUND ) GO TO 7
              IF( ABS(AMATRX(K,I)) .eq. ELM_ABSMAX )THEN
                 IASP_DC = K
                 IBSP_DC = I
                 INDICES_FOUND = .TRUE.
!                WRITE(LUPRI,*) ' indices K, I',K,I
!                WRITE(LUPRI,*) ' AMATRX(K,I)',AMATRX(K,I)
              END IF
    7       CONTINUE
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
 
!     row and column index determined - translate to LUCIAREL symmetry indices
      IF( .NOT. INDICES_FOUND ) GOTO 999
!
!     test print:    
      if(ntest .GT. 0)then
        write(lupri,'(/a,i5,a,i5)') '  SPINOR indices: index a',
     &                                 IASP_DC,' and SPINOR '//
     &                                 'index b',IBSP_DC
      END IF

      select case(dougrp)
        case(5,6,8,10) ! double group without inversion symmetry

          if(iasp_dc.gt.rowdim/2)then
            irassign = imosp_dirac_counter2(iasp_dc-rowdim/2)
          else
            irassign = imosp_dirac_counter1(iasp_dc)
          end if
          if(ibsp_dc.gt.coldim/2)then
            icassign = imosp_dirac_counter2(ibsp_dc-coldim/2)
          else
            icassign = imosp_dirac_counter1(ibsp_dc)
          end if

        case(4,7,11) ! double group with inversion symmetry

!         row index - ungerade or gerade
#ifdef LUCI_DEBUG
          print *, 'na1, na2, iasp_dc, ibsp_dc',
     &              na1, na2, iasp_dc, ibsp_dc
#endif
          if(iasp_dc.gt.2*na1)then ! ungerade
            if(iasp_dc.gt.(2*na1+na2))then ! barred
              irassign = imosp_dirac_counter2(iasp_dc-na1-na2)
            else ! unbarred
              irassign = imosp_dirac_counter1(iasp_dc-na1)
            end if
          else ! gerade
            if(iasp_dc.gt.na1)then ! barred
              irassign = imosp_dirac_counter2(iasp_dc-na1)
            else ! unbarred
              irassign = imosp_dirac_counter1(iasp_dc)
            end if
          end if
!         column index - ungerade or gerade
          if(ibsp_dc.gt.2*na1)then ! ungerade
            if(ibsp_dc.gt.(2*na1+na2))then ! barred
              icassign = imosp_dirac_counter2(ibsp_dc-na1-na2)
            else ! unbarred
              icassign = imosp_dirac_counter1(ibsp_dc-na1)
            end if
          else ! gerade
            if(ibsp_dc.gt.na1)then ! barred
              icassign = imosp_dirac_counter2(ibsp_dc-na1)
            else ! unbarred
              icassign = imosp_dirac_counter1(ibsp_dc)
            end if
          end if
        case default
          call quit('*** error in GET_KRCI_SPINOR_IND: unknown double'//
     &              ' group. ***')
      end select ! double group switch (o/w inversion symmetry)

      if(ntest .gt. 0)then
        print *, 'my irrep assigment row - col ==> ',irassign,icassign
      end if

!     assign luciarel symmetry irreps
!     ===============================
      IASP_KRCI_SYM = irassign
      IBSP_KRCI_SYM = icassign
!     ===============================

      if(ntest .gt. 0)then
        write(lupri,'(/a,/a)') 
     &  '  GET_KRCI_SPINOR_IND:        IASP_DC          IBSP_DC   ',
     &  '  (indices pair found)     _____________    _____________'
C                                                  
        write(lupri,'(i35,i17)') iasp_dc,ibsp_dc
        write(lupri,'(/a,/a)') 
     &  '                           IASP_KRCI_SYM    IBSP_KRCI_SYM',
     &  '  (irreps assigned)        _____________    _____________'
!                                                  
        write(lupri,'(i35,i17)') iasp_krci_sym,ibsp_krci_sym
      end if

 999  CALL QEXIT('GET_KRCI_SPINOR_IND')
#ifdef LUCI_DEBUG
      ntest = 000 ! debug
#endif
#undef LUCI_DEBUG
      end
!***********************************************************************

      SUBROUTINE GET_XPROP_REOD_KRCI(LUF1MAT,WRK,LWRK,CIRUN)
!***********************************************************************
!
!     property operator matrix in DIRAC (Molfdir format) order. 
!     Reorder to KR-CI format according to SIGDEN formalism and store 
!     on file KRMC_FOCK (file handle LUF1MAT).
!                                  
!     Written by S. Knecht - Sep 2008
!
!     Last revision :
!
!***********************************************************************
      use symmetry_setup_krci
      use memory_allocator
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
* boson symmetry info
#include "maxorb.h"
#include "../krmc/dcbbos.h"
#include "dcborb.h"
#include "dgroup.h"
* memory
#include "wrkspc.inc"
#include "ipoist8.inc"
C print flags
#include "cprnt.inc"
C required common blocks
#include "mxpdim.inc"
#include "crun.inc"
#include "cintfo.inc"
#include "ctcc.inc"
#include "glbbas.inc"
C information for property run
#include "krciprop.h"
#include "dcbxpr.h"
#if defined (VAR_MPI)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#endif
#include "parluci.h"
      CHARACTER CIRUN*6, FILELAB*8, FILELAB_REOD*8, RODLAB*2
      LOGICAL MCSCF, RUN_SHMEM_T
      DIMENSION WRK(LWRK)
      real*8, allocatable :: XXXM(:)
C
      CALL QENTER('GET_XPROP_REOD_KRCI')
C
      RODLAB = 'RO'
      MZ = MIN(NZ,2)
      IMEMCK = 0
      MCSCF = .true.
#if defined (VAR_MPI)
      MASTER = MPARID
      MYPROC = MYTID
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
C     output file handle in KR-CI module
      LUWRT = LUPRI

C     Compute offset from WRK in parameter list
C     to WORK in common in wrkspc.inc. WORK is used internally in LUCIAREL.
C     This offset is added to all Ksomething addresses in WORK(Ksomething)
C     so that WRK(1:LWRK) are the memory addresses really used.

C     Note that K_OFFSET is in ipoist8.inc and is always i*8
C     because the offset can be bigger than what can be stored in i*4.
C     ** All Ksomething used in WORK(Ksomething) must always be i*8 and declared in ipoist8.inc **
C     See also subroutine MEMMAR.

      CALL COMPOFF(WRK,WORK,K_OFFSET)
C     *k_offset = (long) (wrk-work+1);

      MXPWRD = LWRK
      KADD = MXPWRD

C     initialize memory
      call memmar(K_OFFSET,KADD,'INI   ',IDUMMY,'DUMMY')
C
C     Internal string information (stored in WORK, bases in /STRBAS/)
      CALL STRINF_GAS_REL(WORK(1),IPRSTR)
C
C     Internal subspaces
      CALL LCISPC_REL(LCSBLK,IPRCIX,CIRUN)
C
C     IDGSTYP determines the integral classes and corresponding operator
C     classes depending on the double group in use.
C
      IDGSTYP = -1
C
      IF ( PNTGRP .le. 4 ) THEN
        IDGSTYP = 1
      ELSE IF ( PNTGRP .eq. 6 .or. PNTGRP .eq. 5 ) THEN
        IDGSTYP = 3
      ELSE IF ( PNTGRP .ge. 7)THEN
        IDGSTYP = 2
      END IF
C
C     operator type
      IF (IDGSTYP .le. 2 ) THEN
         IHTYPE = 4
      ELSE IF (IDGSTYP.eq.3) THEN
         IHTYPE = 5
      ELSE
        WRITE(LUWRT,'(/A,I5)')'  *** ERROR in GET_XPROP_REOD_KRCI.'//
     &                        ' Unallowed value of IDGSTYP :',IDGSTYP
        CALL QUIT(' *** ERROR in GET_XPROP_REOD_KRCI: 
     &               Unallowed value of IDGSTYP.')
      END IF
C     transfer to common block in krciprop.h
      IHTYPE_X  = IHTYPE
C     set memory mark
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'MARK ',IDUM,'PROPRO')
C
C     allocate for property operator matrix
C
      CALL MEMMAR(KXPROPM,4*N2ASHX*MZ,'ADDL  ',2,'IXPRPM')
C     do reordering - this flag controls behavior in subroutine ITRCTL
      DOPROPREOD  = .TRUE.
C
C     ... loop over property operators
      DO 100 II = 1, NPROP_KRCI
C
        ISYM_T = ISYMOPPRP_KRCI(II)
C       ... nothing to do for "inactive" operators
        IF( ISYM_T .le. 0 ) GOTO 100
C
C       ... no two-electron terms
        NO_TTWO         = 1
        IPRHAM_CI       = 0
        INDXPR          = LPROP_KRCI(II)
        ISYM_XXX        = IPRPSYM(INDXPR)
C
C
        CALL SET_HOP_DBG(IHTYPE,NO_TTWO,ISYM_T,IPRHAM_CI)
C
C       read property matrix from file
        CALL DZERO(WORK(KXPROPM),4*N2ASHX*MZ)
        FILELAB = PRPNAM(INDXPR)(1:8)
        CALL KRCI_PRPFILE(LUF1MAT,FILELAB,WORK(KXPROPM),
     &                    4*N2ASHX*MZ,1)
C       ... reorder always at this point: RUN_SHMEM_T = .FALSE.
        RUN_SHMEM_T = .FALSE.
#if defined DEBUG_LUCIPROP
        WRITE(LUWRT,'(/A,A8/)') '  Reorder step for property matrix: ',
     &                             FILELAB
        CALL WRTMATMN(WORK(KXPROPM),1,4*N2ASHX*MZ,1,4*N2ASHX*MZ,LUWRT)
#endif
        CALL PICASSO(IPRTRA,MCSCF,NL2D,IUB2D,ISM2D,
     &               WORK(KXPROPM),DUMMY,DUMMY,DUMMY,IDUMMY,
     &               ibosym,CIRUN,RUN_SHMEM_T)
C
C       write sorted property matrix back to file
        WRITE(FILELAB_REOD,'(A6,A2)') PRPNAM(INDXPR)(1:6),RODLAB
#if defined DEBUG_LUCIPROP
!SK DEBUG
        WRITE(LUWRT,'(/A,A8/)') '  Final matrix on file with label: ',
     &                             FILELAB_REOD
        WRITE(LUWRT,'(A,I18/)') '  total number of integrals on KT_CC ',
     &                             NALLINT
        IF( FILELAB_REOD .eq. 'Y dipoRO')THEN
          WRITE(LUWRT,'(/A)') '  DEBUG READING DUMMY Y-OPERATOR MATRIX'
!radovan: presently alloc cannot handle integer(8) size with ILP32
!         call alloc(XXXM,NALLINT)
          allocate(XXXM(NALLINT))
          open(unit=19,file='y_int_test', 
     &         status='UNKNOWN',form='FORMATTED')
          rewind(19)
          READ(19,'(1D22.14)') (XXXM(I),I=1,NALLINT)
          close(19)
          LINT_WRITE = NALLINT
          WRITE(LUWRT,'(/A)') '  DEBUG WRITING DUMMY Y-OPERATOR MATRIX'
          CALL WRTMATMN(XXXM,1,LINT_WRITE,1,LINT_WRITE,LUWRT)
          CALL DCOPY(LINT_WRITE,XXXM,1,WORK(KT_CC),1)
          WRITE(LUWRT,'(/A/)') '  USING DEBUG DUMMY Y-OPERATOR MATRIX'
!radovan: presently alloc cannot handle integer(8) size with ILP32
!         call dealloc(XXXM)
          deallocate(XXXM)
        END IF
        IF( FILELAB_REOD .eq. 'X dipoRO')THEN
          WRITE(LUWRT,'(/A)') '  DEBUG READING DUMMY X-OPERATOR MATRIX'
!radovan: presently alloc cannot handle integer(8) size with ILP32
!         call alloc(XXXM,NALLINT)
          allocate(XXXM(NALLINT))
          open(unit=19,file='x_int_test', 
     &         status='UNKNOWN',form='FORMATTED')
          rewind(19)
          READ(19,'(1D22.14)') (XXXM(I),I=1,NALLINT)
          close(19)
          LINT_WRITE = NALLINT
          WRITE(LUWRT,'(/A)') '  DEBUG WRITING DUMMY X-OPERATOR MATRIX'
          CALL WRTMATMN(XXXM,1,LINT_WRITE,1,LINT_WRITE,LUWRT)
          CALL DCOPY(LINT_WRITE,XXXM,1,WORK(KT_CC),1)
          WRITE(LUWRT,'(/A/)') '  USING DEBUG DUMMY X-OPERATOR MATRIX'
!radovan: presently alloc cannot handle integer(8) size with ILP32
!         call dealloc(XXXM)
          deallocate(XXXM)
        END IF
!SK DEBUG
#endif
C
C       store number of non-zero integrals for this operator
        IN1ELPRP_KRCI(II) = 0
        LINT_WRITE        = NALLINT
        IN1ELPRP_KRCI(II) = LINT_WRITE
#if defined DEBUG_LUCIPROP
        CALL WRTMATMN(WORK(KT_CC),1,LINT_WRITE,1,LINT_WRITE,LUWRT)
#endif
        CALL KRCI_PRPFILE(LUF1MAT,FILELAB_REOD,WORK(KT_CC),
     &                    LINT_WRITE,0)
        XT_TEST = 0.0D0
        XT_TEST = DDOT(LINT_WRITE,WORK(KT_CC),1,WORK(KT_CC),1)
C
C       reset property operator entries for operator with zero entries 
C       in reordered integral matrix
        IF( XT_TEST .eq. 0.0D0 ) THEN
           ISYMOPPRP_KRCI(II) = 0
           IN1ELPRP_KRCI(II)  = 0
        END IF
 100  CONTINUE
C     reordering done.
      DOPROPREOD  = .FALSE.
C
C     release memory
C
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'FLUSM',IDUM,'PROPRO')
C
 999  CALL QEXIT('GET_XPROP_REOD_KRCI')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XPROP_KRCI(XPROPMAT,IHTYPE,IPROCLIST,IGROUPLIST,CIRUN,
     &                      MZ)
C***********************************************************************
C
C     calculate property matrix. 
C     Reordered property integrals (KR-CI SIGDEN formalism) are stored 
C     on file KRMC_FOCK (file handle LUPROP_MAT).
C                                  
C     Written by S. Knecht - Oct 2008
C
C     Last revision :
C
C***********************************************************************
      use krci_properties
      use symmetry_setup_krci
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
* boson symmetry info
#include "maxorb.h"
#include "../krmc/dcbbos.h"
#include "dcborb.h"
#include "dgroup.h"
* memory
#include "wrkspc.inc"
#include "ipoist8.inc"
C print flags
#include "cprnt.inc"
C required common blocks
#include "mxpdim.inc"
#include "crun.inc"
#include "cands.inc"
#include "cicisp.inc"
#include "cgas.inc"
#include "oper.inc"
#include "cintfo.inc"
#include "clunit.inc"
#include "ctcc.inc"
#include "glbbas.inc"
#include "cstate.inc"
C information for property run
#include "krciprop.h"
#include "dcbxpr.h"
#if defined (VAR_MPI)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
#endif
#include "parluci.h"
     
      DIMENSION XPROPMAT(NPROP_ROOTS_KRCI,NPROP_ROOTS_KRCI,
     &                   NPROP_KRCI,MZ)
      DIMENSION IPROCLIST(*), IGROUPLIST(*)
      CHARACTER CIRUN*6,FILELAB*8,FILELAB_REOD*8,RODLAB*2,x_NAME*800
      LOGICAL MCSCF, RUN_SHMEM_T, COPY_NEW_C, EX
      INTEGER*8 IL_KTCC8
      integer, intent(in) :: mz
C
      CALL QENTER('XPROP_KRCI')
      NTESTX = 10
C
C     set memory mark
      IDUM = 0
      CALL MEMMAR(KDUMMY,IDUM,'MARK  ',IDUM,'XPROP ')
C
      RODLAB = 'RO'
C     set required variables
      IMEMCK  = 0
      INTIMP  = 7
C     ... no two-electron excitations
      NO_TTWO   = 1
      MCSCF     = .TRUE.
      DOSIGPROP = .TRUE.
C     hardwired to 1
      ICSPC  = NCMBSPC
      ISSPC  = NCMBSPC
C     real or complex eigenvectors
      IRC = 2
      IF( NZ .eq. 1 ) IRC = 1
      IRC_SAVE = IRC
C     scratch
      NBATXS = 0
C
C     property calculation setup section
C     __________________________________
C
      IF(MYPROC .eq. MASTER ) REWIND LUPROP_MAT
C
C     determine max. dimension of T and VECx( x = 1-3)
C     ------------------------------------------------
      IL_KTCC          = 1
      L_COMBI_MAX      = 0
      L_COMBI_MAX_SCR  = 0
      L_COMBI_MAX_SCRX = 0
C
      DO I = 1, NPROP_KRCI
        IL_KTCC  = MAX(IL_KTCC,IN1ELPRP_KRCI(I))
      END DO
      IL_KTCC8 = IL_KTCC
C
      CALL IZERO(IXSYMLIST,MXPROPKRCI_SYM)
C
C     IXSYMLIST will now be initialized and later adapted in 
C     Z_BLKFO_XPROP such that, e.g. IXSYMLIST(1) == '4', which is 
C     therefore the first "active" symmetry irrep, gives the 
C     right offset for NPARBLOCK(*,'1').
C
      DO IPROP_X = 1, NPROP_KRCI
         IF( ISYMOPPRP_KRCI(IPROP_X) .gt. 0 ) THEN
           DO I = 1, NPROP_ROOTS_KRCI
              ILC_C = XISPSM(ISYMEIG_KRCI(I),ICSPC)
              JSSM   = 
     &        IDBGMULT(ISYMOPPRP_KRCI(IPROP_X),ISYMEIG_KRCI(I))
              ILC_S  = XISPSM(JSSM,ISSPC)
C             mark "active" C/sigma symmetry irreps
              IXSYMLIST(ISYMEIG_KRCI(I)) = ISYMEIG_KRCI(I)
              IXSYMLIST(JSSM)  = JSSM
              L_COMBI_MAX_SCRX = MAX(ILC_C,ILC_S)
              L_COMBI_MAX_SCR  = MAX(L_COMBI_MAX_SCR,L_COMBI_MAX_SCRX)
           END DO
         END IF
      END DO
C
C     store symmetry pairs - active pairs are marked by 1.
      CALL MEMMAR(KXSYMPAIRS,MXPROPKRCI_SYM**2,'ADDL  ',1,'KXPAIR')
      CALL IZERO(WORK(KXSYMPAIRS),MXPROPKRCI_SYM**2)
C
      CALL GET_SC_PAIRS(WORK(KXSYMPAIRS),MXPROPKRCI_SYM,NPROP_KRCI,
     &                  ISYMOPPRP_KRCI,NPROP_ROOTS_KRCI,ISYMEIG_KRCI)
C
C     largest number of dets - used for MPI file offset calculation
      L_COMBI_MAX = MAX(L_COMBI_MAX_SCR,L_COMBI_MAX)
C     initialize L_COMBI to avoid crash in block distribution 
      L_COMBI     = L_COMBI_MAX
C
C     T block allocation
      CALL MEMMAR_I8(KT_CC,IL_KTCC8,'ADDL  ',2,'KT_CC ')
      CALL DZERO8(WORK(KT_CC),IL_KTCC8)
C     number of blocks per symmetry irrep
      CALL MEMMAR(KXSYM_BLK,MXPROPKRCI_SYM,'ADDL  ',1,'KXSYMB')
      CALL IZERO(WORK(KXSYM_BLK),MXPROPKRCI_SYM)
C
C     set block (parallel: distribution list) arrays and vector dimensions 
      CALL Z_BLKFO_XPROP(IDC,NMS2VAL,ICSPC,KNODELIST,KBLOCKLIST,
     &                   KRCCTOS,KXSYMDISTR,WORK(KXSYM_BLK),
     &                   WORK(KXSYMPAIRS),IRC,IGROUPLIST)
C
      CALL GET_3BLKS_REL(KVEC1,KVEC2,LSCR2,CIRUN)
C
C     set file offset/view for MPI-I/O files
#if defined (VAR_MPI2)
C     set active blocks for C scratch file ILUC
      MY_ACT_BLK_ALL = NUM_BLOCKS
C     general offset
      MY_LUC_OFF = 0
C
      CALL CALC_OFF_MPI_FILE_XPRP(IGROUPLIST,WORK(KBLOCKLIST),
     &                            WORK(KNODELIST),NPROP_ROOTS_KRCI,
     &                            NPROP_ROOTS_KRCI)
C
C     length for allocation of file arrays
C     WRITE(LUWRT,*) ' NPROP_ROOTS_KRCI,MY_ACT_BLK2,NUM_BLOCKS,...',
C    &NPROP_ROOTS_KRCI,MY_ACT_BLK2,NUM_BLOCKS,NUM_BLOCKS2
      IALL_LU1 = NPROP_ROOTS_KRCI * MY_ACT_BLK2
      IALL_LU2 = NPROP_ROOTS_KRCI * MY_ACT_BLK2
      IALL_LUC = NUM_BLOCKS2
C     WRITE(LUWRT,*) ' RMEMCHK before list alloc',
C    &                 IALL_LU1,IALL_LU2,IALL_LUC
C     CALL RMEMCHK
C
C     allocate file arrays
      CALL MEMMAR(KILU1LIST,IALL_LU1,'ADDS  ',1,'LU1LST')
      CALL MEMMAR(KILU2LIST,IALL_LU2,'ADDS  ',1,'LU2LST')
      CALL MEMMAR(KILUCLIST,IALL_LUC,'ADDS  ',1,'LUCLST')
C     initialize ...
      CALL IZERO(WORK(KILU1LIST),IALL_LU1)
      CALL IZERO(WORK(KILU2LIST),IALL_LU2)
      CALL IZERO(WORK(KILUCLIST),IALL_LUC)
C
C     T block distribution - reordered arrays of property integrals
      CALL DIST_XPROP(WORK(KT_CC),LUPROP_MAT)
C     WRITE(LUWRT,*) ' RMEMCHK before DIST_XCVEC'
C     CALL RMEMCHK
C
C     distribute C vectors from MASTER to slaves for all active irreps
      CALL DIST_XCVEC(WORK(KVEC1),LUC,ILU1,WORK(KILU1LIST),
     &                WORK(KNODELIST),WORK(KBLOCKLIST),WORK(KXSYM_BLK),
     &                IXSYMLIST,IXSYMACT,IRC)
C     WRITE(LUWRT,*) ' RMEMCHK after DIST_XCVEC'
C     CALL RMEMCHK
C     scratch arrays
      CALL MEMMAR(KXR,NPROP_ROOTS_KRCI,'ADDL  ',2,'KXR   ')
      IF(MZ.eq.2) THEN 
        CALL MEMMAR(KXI,NPROP_ROOTS_KRCI,'ADDL  ',2,'KXI   ')
      ELSE
        CALL MEMMAR(KXI,               0,'ADDL  ',2,'KXI   ')
      END IF
      CALL DZERO(WORK(KXR),NPROP_ROOTS_KRCI)
      IF( MZ .eq. 2 ) CALL DZERO(WORK(KXI),NPROP_ROOTS_KRCI)
C
C     end of setup section
C     --------------------
#endif
C
C     computation section
C     ___________________
C
C     initialize start codes for ICSM_OLD and NVC_LUC (# of vec per irrep)
      ICSM_OLD  = -1
      NVC_LUC   = -3
#if defined (VAR_MPI2)
C     ... C batches of blocks 
      CALL MEMMAR(KPCLBXX,MXNTTS,  'ADDL  ',1,'CLBTX ')
      CALL MEMMAR(KPCLEXX,MXNTTS,  'ADDL  ',1,'CLEBTX')
      CALL MEMMAR(KPCI1XX,MXNTTS,  'ADDL  ',1,'CI1BTX')
      CALL MEMMAR(KPCIBXX,8*MXNTTS,'ADDL  ',1,'CIBTX ')
C     ... sigma batches of blocks
      CALL MEMMAR(KPCLBXXS,MXNTTS,  'ADDL  ',1,'CLBTXY')
      CALL MEMMAR(KPCLEXXS,MXNTTS,  'ADDL  ',1,'CLEBTY')
      CALL MEMMAR(KPCI1XXS,MXNTTS,  'ADDL  ',1,'CI1BTY')
      CALL MEMMAR(KPCIBXXS,8*MXNTTS,'ADDL  ',1,'CIBTXY')
#endif
      WRITE(LUWRT,'(/A/)') 
     & '  ********* start of property matrix calculation *********'
C
C     outer loop over all eigenstates
      DO 100 ISTATE_X = 1, NPROP_ROOTS_KRCI
C
         ICSM     = ISYMEIG_KRCI(ISTATE_X)
         L_COMBI  = XISPSM(ICSM,ICSPC)
         IF( ICSM .ne. ICSM_OLD )
     &     WRITE(LUWRT,'(/A,I4/)') 
     &     '  * looping now over all states of symmetry',ICSM
#if defined (VAR_MPI2)
C        erase T operator file offset
         IOFFSET_SCRATCH = 0
C        get block partitioning for new active irrep
         IF( ICSM .ne. ICSM_OLD )THEN
C          which "active" C irrep?
           IX_CT = IGET_ACTIVE_IRREP(IXSYMLIST,IXSYMACT,ICSM)
C
           CALL Z_BLKFO_XPRP_PAR(IDC,NMS2VAL,ICSPC,ICSM,WORK(KPCLBXX),
     &                           WORK(KPCLEXX),WORK(KPCI1XX),
     &                           WORK(KPCIBXX),NBATX,WORK(KNODELIST),
     &                           IX_CT)
         END IF
#else
         NVC_LUC = NVC_LUC + 1
         IF( ICSM .ne. ICSM_OLD )THEN
C          close C vector file for previous irrep (-2 == start code +1)
           IF( NVC_LUC .ne. -2 ) CLOSE(UNIT=LUC,STATUS='KEEP')
C          reset NVC_LUC and open new C vector file
           NVC_LUC = 1
           CALL GET_LUC_FILE(LUC,ICSM,EX)
         END IF
         REWIND LUPROP_MAT
#endif
C        new ICSM_OLD
         ICSM_OLD = ICSM
C        copy new C vector blocks from ILU1 to ILUC
         COPY_NEW_C = .TRUE.
C
C        loop over all property operators
         DO 200 IPROP_X = 1, NPROP_KRCI
C
            ISYM_T = ISYMOPPRP_KRCI(IPROP_X)
C           only zero elements in property operator matrix
            IF( ISYM_T .le. 0 ) GOTO 200
C
C           property pointer
            INDXPR      = LPROP_KRCI(IPROP_X)
C
C           reset operator specific common blocks
            CALL SET_HOP_DBG(IHTYPE,NO_TTWO,ISYM_T,IPRHAM_CI)
C
            LENGTH_BLK = IN1ELPRP_KRCI(IPROP_X)
C           erase old property integrals
            CALL DZERO8(WORK(KT_CC),IL_KTCC8)
C
C           read new property integrals from disk
#if defined (VAR_MPI2)
            call interface_mpi_FILE_READ_AT(ILPRP_X,IOFFSET_SCRATCH,
     &                            WORK(KT_CC),LENGTH_BLK,ISTAT)
            IOFFSET_SCRATCH = IOFFSET_SCRATCH + LENGTH_BLK
C
#else
            WRITE(FILELAB_REOD,'(A6,A2)') PRPNAM(INDXPR)(1:6),RODLAB
            CALL KRCI_PRPFILE(LUPROP_MAT,FILELAB_REOD,WORK(KT_CC),
     &                        LENGTH_BLK,1)
#endif
!#define MOD_DEBUG
#ifdef MOD_DEBUG
C
C           debug print
            WRITE(LUWRT,*) ' matrix of property:',PRPNAM(INDXPR)
            CALL WRTMATMN(WORK(KT_CC),1,LENGTH_BLK,1,LENGTH_BLK,LUWRT)
            WRITE(LUWRT,*) ' symmetry of property    ==>', ISYM_T
            WRITE(LUWRT,*) ' symmetry of sigma state ==>', 
     &                       IDBGMULT(ICSM,ISYM_T)
#ifdef VAR_MPI2
            WRITE(LUWRT,*) ' active sigma irrep      ==>', 
     &                       IGET_ACTIVE_IRREP(IXSYMLIST,IXSYMACT,ISSM)
#endif
#endif
            ISSM      = IDBGMULT(ICSM,ISYM_T)
            L_COMBI_S =   XISPSM(ISSM,ISSPC)
C           which "active" S irrep?
#ifdef VAR_MPI2
            IX_ST     = IGET_ACTIVE_IRREP(IXSYMLIST,IXSYMACT,ISSM)
#endif
C
C           compute new sigma vector = T_prop x C
C           -------------------------------------
            CALL XSIGDEN_CTRL(WORK(KVEC1),WORK(KVEC2),WORK(KT_CC),
     &                        COPY_NEW_C,NVC_LUC,MZ
#if defined (VAR_MPI2)
     &                        ,WORK(KILU1LIST),WORK(KILU2LIST),
     &                         WORK(KILUCLIST),WORK(KPCLBXX),
     &                         WORK(KPCLEXX),WORK(KPCI1XX),
     &                         WORK(KPCIBXX),NBATX,WORK(KBLOCKLIST),
     &                         WORK(KNODELIST),IGROUPLIST,IPROCLIST,
     &                         WORK(KRCCTOS),WORK(KXSYMPAIRS),ICSM,ISSM,
     &                         IX_CT,IX_ST,ISTATE_X-1
#endif
     &                        )

C
C           inner loop over all eigenstates
#if defined (VAR_MPI2)
C 
!#define MOD_DEBUG
#ifdef MOD_DEBUG
            WRITE(LUWRT,*) 'par - <S_ISSM|prop|C_ICSM> for C_#,IPROP_X',
     &                         ISSM,PRPNAM(INDXPR),ICSM,ISTATE_X,IPROP_X
#endif
           CALL Z_BLKFO_XPRP_PAR(IDC,NMS2VAL,ISSPC,ISSM,WORK(KPCLBXXS),
     &                           WORK(KPCLEXXS),WORK(KPCI1XXS),
     &                           WORK(KPCIBXXS),NBATXS,WORK(KNODELIST),
     &                           IX_ST)

            CALL INPRDD_REAL_CPLX_XPRP(ILU2,ILU1,WORK(KVEC1),WORK(KVEC2)
     &                                ,WORK(KXR),WORK(KXI),NBATXS,
     &                                 WORK(KPCLBXXS),WORK(KPCLEXXS),
     &                                 WORK(KPCI1XXS),WORK(KPCIBXXS),
     &                                 MY_LU2_OFF,MY_LU1_OFF,
     &                                 WORK(KILU2LIST),WORK(KILU1LIST),
     &                                 NPROP_ROOTS_KRCI,ISSM,
     &                                 ISYMEIG_KRCI,MZ,1,ISTATE_X)
C
C           debug printing
#ifdef MOD_DEBUG
C
            WRITE(LUWRT,*) ' WORK(KXR) '
            CALL WRTMATMN(WORK(KXR),1,NPROP_ROOTS_KRCI,1,
     &                    NPROP_ROOTS_KRCI,LUWRT)
            IF( MZ.eq.2)THEN
              WRITE(LUWRT,*) ' WORK(KXI) '
              CALL WRTMATMN(WORK(KXI),1,NPROP_ROOTS_KRCI,1,
     &                      NPROP_ROOTS_KRCI,LUWRT)
            END IF
#endif
C
C           lz is pure imaginary operator - scale with -1.0
C           FIXME: introduce proper scaling for each imaginary 
C                  operator
            DMONE = -1.0D0
            IF( PRPNAM(INDXPR) .eq. 'z-angular mom' )THEN
              CALL DSCAL(NPROP_ROOTS_KRCI,DMONE,WORK(KXR),1)
              IF(MZ.eq.2) CALL DSCAL(NPROP_ROOTS_KRCI,DMONE,WORK(KXI),1)
            END IF
C
            CALL DAXPY(NPROP_ROOTS_KRCI,1.0d0,WORK(KXR),1,
     &                 XPROPMAT(1,ISTATE_X,IPROP_X,1),1)
            IF( MZ .eq. 2 ) 
     &      CALL DAXPY(NPROP_ROOTS_KRCI,1.0d0,WORK(KXI),1,
     &                 XPROPMAT(1,ISTATE_X,IPROP_X,2),1)
C
C           DO 300 ISTATE_Y = 1, NPROP_ROOTS_KRCI
C             IF( ISSM .ne. ISYMEIG_KRCI(ISTATE_Y) ) GOTO 300
C             XPROPMAT(ISTATE_Y,ISTATE_X,IPROP_X,1) =
C    &                                              WORK(KXR+ISTATE_Y-1)
C             WRITE(LUWRT,*) 'ISTATE_Y,ISTATE_X, PRPNAM(INDXPR) and'//
C    &                       ' XPROPMAT(Y,X,PROP,1)',ISTATE_Y,ISTATE_X,
C    &                        PRPNAM(INDXPR),
C    &        XPROPMAT(ISTATE_Y,ISTATE_X,IPROP_X,1)
C             IF( MZ .eq. 2 ) XPROPMAT(ISTATE_Y,ISTATE_X,IPROP_X,2) =
C    &                                              WORK(KXI+ISTATE_Y-1)
C             WRITE(LUWRT,*) 'ISTATE_Y,ISTATE_X, PRPNAM(INDXPR) and'//
C    &                       ' XPROPMAT(Y,X,PROP,2)',ISTATE_Y,ISTATE_X,
C    &                        PRPNAM(INDXPR),
C    &        XPROPMAT(ISTATE_Y,ISTATE_X,IPROP_X,2)
!           DO 300 ISTATE_Y = 1, NPROP_ROOTS_KRCI
! stefan test: symmetrize
!             xpropmat(istate_x,istate_y,iprop_x,1) = 
!    &        xpropmat(istate_y,istate_x,iprop_x,1)   
!             if(mz == 2)
!    &        xpropmat(istate_x,istate_y,iprop_x,2) = 
!    &        xpropmat(istate_y,istate_x,iprop_x,2)   
!300        CONTINUE
C           ^ end of inner loop over all eigenstates
#else
!#define MOD_DEBUG
#ifdef MOD_DEBUG
            WRITE(LUWRT,*) 'seq - <S_ISSM|prop|C_ICSM> for C_#,IPROP_X',
     &                         ISSM,PRPNAM(INDXPR),ICSM,ISTATE_X,IPROP_X
#endif
!           close current c-vector file
            close(unit=luc,status="keep")
            is_correct_sym = 0
            
            do istate_y = 1, nprop_roots_krci
C
              call get_inner_product_sigmaprop_dot_c(work(kvec1),
     &                                               work(kvec2),
     &                                               xpropmat(
     &                                               istate_y,
     &                                               istate_x,
     &                                               iprop_x,1),
     &                                               xpropmat(
     &                                               istate_y,
     &                                               istate_x,
     &                                               iprop_x,mz),
     &                                               issm,            
     &                                               isymeig_krci(
     &                                               istate_y),
     &                                               luhc,           
     &                                               is_correct_sym,    
     &                                               mz,     
     &                                               xsymflab(
     &                                               isymeig_krci(
     &                                               istate_y)))

              if(prpnam(indxpr).eq.'z-angular mom')then
                            xpropmat(istate_y,istate_x,iprop_x,1) = 
     &                     -xpropmat(istate_y,istate_x,iprop_x,1)
                if(mz.eq.2) xpropmat(istate_y,istate_x,iprop_x,2) = 
     &                     -xpropmat(istate_y,istate_x,iprop_x,2)
              end if

C
            end do
C           ^ end of inner loop over all eigenstates

C           reopen current c-vector file
            call get_luc_file(luc,icsm,ex)
#endif
 200     CONTINUE
C        ^ end of loop over all properties
 100  CONTINUE
C     ^ end of outer loop over all eigenstates
C
      WRITE(LUWRT,'(/A/)') 
     & '  *********  end of property matrix calculation  *********'
C
#if defined (VAR_MPI2)
C
C     update of complete property matrix
C     ----------------------------------
      CALL MEMMAR(KXPROP_SCR,LXPRPKRCI,'ADDL  ',2,'XPRPSC')
      CALL DZERO(WORK(KXPROP_SCR),LXPRPKRCI)
C
      CAll redvec(XPROPMAT,WORK(KXPROP_SCR),LXPRPKRCI,2,
     &                op_MPI_SUM,global_communicator,MASTER)
      IF( MYPROC .eq. MASTER) 
     &    CALL DCOPY(LXPRPKRCI,WORK(KXPROP_SCR),1,XPROPMAT,1)
#endif

!     symmetrize
      DO I = 1, NPROP_ROOTS_KRCI
        DO j = 1, NPROP_ROOTS_KRCI
          DO k = 1, NPROP_KRCI
            DO l = 1, mz
              xpropmat(i,j,k,mz) = xpropmat(j,i,k,mz)
            end do
          end do
        end do
      end do
C
C     print section
C     _____________
C
!#define MOD_DEBUG
#ifdef MOD_DEBUG
      NTESTX = 10
#else
      NTESTX = 00
#endif
      IF( NTESTX .gt. 5 .and. myproc.eq.0)THEN
        WRITE(LUWRT,*) ' calculation of XPROPMAT finished'
        WRITE(LUWRT,*) ' content of XPROPMAT'
        CALL WRTMATMN(XPROPMAT,1,LXPRPKRCI,1,LXPRPKRCI,LUWRT)
      END IF
C     release memory
 9999 call MEMMAR(KDUMMY,IDUMMY,'FLUSM ',IDUM,'XPROP ')
C
#if defined (VAR_MPI2)
C
C     ... close scratch file with reordered property integrals (opened
C     in DIST_XPROP)
      call interface_mpi_FILE_CLOSE(ILPRP_X)
#endif
C
C     set back to default
      DOSIGPROP = .FALSE.
C      
      CALL QEXIT('XPROP_KRCI')
#undef MOD_DEBUG
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     end of file property_r.F
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
