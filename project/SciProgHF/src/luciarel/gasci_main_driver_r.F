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
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GASCI_REL(F1,F2,ISM,ISPC,IPRNT,EREF,
     &                     ECORE_KRMC,CIRUN,IPROCLIST,IGROUPLIST)
*
!     external modules
      use memory_allocator
      use interface_to_mpi
*
* CI optimization in GAS space number ISPC for symmetry ISM
*
* Jeppe Olsen, Winter of 1995
*              Winter of 1998 : Modified for relativistic calculations
*
!     internal modules
      use mospinor_info
      use symmetry_setup_krci

      IMPLICIT REAL*8(A-H,O-Z)
#if defined (VAR_MPI2)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "ipoist8.inc"
*
      character*6 CIRUN
      real(8), intent(inout) :: eref
*
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "clunit.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "ctcc.inc"
#include "cintfo.inc"
*
*. Common block for communicating with sigma
#include "cands.inc"
*  nat. orb. occ. num.
#include "noccn_inf.inc"
#include "integrals_off.inc"
#include "dgroup.h"
#include "maxorb.h"
#include "../krmc/dcbbos.h"
*     /CECORE/
      COMMON/CECORER/ECORE
      DIMENSION F1(*), F2(*), IPROCLIST(*), IGROUPLIST(*)
      INTEGER*8 L_MAX_DENS_BLK
      LOGICAL CPIC_PASS_D, DIAGONAL_ALLOC
C 
C     some local scratch
      real*8,    allocatable :: eci_start(:)
      real*8,    allocatable :: scratch_dav1(:)
      real*8,    allocatable :: scratch_dav2(:)
      integer,   allocatable :: iscratch_dav1(:)
      integer, parameter     :: pert_dummy = 1
C     memory pointer for TT-block distribution among CPU's
C     memory pointer for TT-block offsets
C     memory pointer for T coefficients
C     memory pointer for T diag block
C     memory pointer for C1 vector block
C     memory pointer for C2 vector block
#if defined (VAR_MPI2)
#include "krmc_shmem.h"
#include "mxdim_mpi2.h"
      INTEGER(KIND=df_MPI_OFFSET_KIND) ITEST_OFF
#if defined VAR_SHMEM
      INTEGER(KIND=df_MPI_OFFSET_KIND) IT_TTOL
      POINTER (MY_TTPL_PTR, IT_TTPL(IDUMMY_BUFF_SZ))
      POINTER (MY_TTOL_PTR, IT_TTOL(IDUMMY_BUFF_SZ))
      POINTER (MY_T_PTR   , T_BUFF(IDUMMY_BUFF_SZ))
      POINTER (MY_T_D_PTR , T_BUFF_D(IDUMMY_BUFF_SZ))
      POINTER (MY_XT_PTR  , XT_BUFF(IDUMMY_BUFF_SZ))
      POINTER (MY_VEC1_PTR, VEC1_BUFF(IDUMMY_BUFF_SZ))
      POINTER (MY_VEC2_PTR, VEC2_BUFF(IDUMMY_BUFF_SZ))
      POINTER (MY_DL_PTR  , IT_DLIST_TMP(IDUMMY_BUFF_SZ))
#else
      integer,   allocatable :: IT_TTPL(:)
      integer,   allocatable :: IT_DLIST_TMP(:)
      integer*8, allocatable :: IT_TTOL(:)
      real*8,    allocatable :: T_BUFF(:)
      real*8,    allocatable :: T_BUFF_D(:)
#endif /* VAR_SHMEM */
#endif /* VAR_MPI2 */
      integer(8)             :: NCALC_MAX_BLK, IGET_MY_T_LEN
      integer*8              :: keci_start

      character*12 SECTID ! external function
C
      CALL GETTIM(CPU0,WALL0)
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'GASCI ')
C
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
C     set some variables ...
      MIN2  = - 2
      NZERO =   0
      IIONE =   1
C     CPIC_PASS_D controls whether the T_DIAG part is 
C     stored on T_BUFF_D or not
C     DIAGONAL_ALLOC keeps track of allocation of T_BUFF_D
      CPIC_PASS_D    = .FALSE.
      DIAGONAL_ALLOC = .FALSE.
C
CSK   NTEST = 3
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ====================================='
        WRITE(LUWRT,*) ' Control has been transferred to GASCI'
        WRITE(LUWRT,*) ' ====================================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,'(A)') '  A few pertinent data : '
        WRITE(LUWRT,*)
        WRITE(LUWRT,'(A,I2)') '  CI space         ',ISPC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of GAS spaces included ',LCMBSPC(ISPC)
        WRITE(LUWRT,'(A,10I3)') ' GAS spaces included           ',
     &               (ICMBSPC(II,ISPC),II=1,LCMBSPC(ISPC))
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Occupation constraints : '
        WRITE(LUWRT,*) '========================= '
        WRITE(LUWRT,*)
        WRITE(LUWRT,*)
        DO JJGASSPC = 1, LCMBSPC(ISPC)
         JGASSPC = ICMBSPC(JJGASSPC,ISPC)
        WRITE(LUWRT,*)
     &  ' Gas space  Min acc. occupation Max acc. occupation '
        WRITE(LUWRT,*)
     &  ' ================================================== '
        DO IGAS = 1, NGAS
          WRITE(LUWRT,'(3X,I2,13X,I3,16X,I3)') IGAS,
     &     IGSOCCX(IGAS,1,JGASSPC),IGSOCCX(IGAS,2,JGASSPC)
        END DO
        END DO
C
       END IF
CSK    NTEST = 0
C
!     already done one routine up
!     ECORE = ECORE_KRMC
      NDET = XISPSM(ISM,ISPC)
      IF(NTEST > 25) 
     &  write(LUWRT,*) ' Number of determinants/combinations  ',NDET
C     transfer NDET to common block LUCIPARREL
      L_COMBI     = NDET
      L_COMBI_MAX = NDET
C
C     real or complex eigenvectors
      IF(ispinfree.EQ.1) THEN
C       real
        IRC = 1
      ELSE
        IRC = 2
C       well, generally complex, but
        if (NZ.eq.1) then
C         quaternion real
          IRC = 1
        end if
      END IF
C     control variable for general file structure
      LBLK = -1
C
C     transfer symmetry of C and sigma to CANDS
      ICSM  = ISM
      ISSM  = ISM
      ICSPC = ISPC
      ISSPC = ISPC
!     determinant code -- nvar == ndet
      NVAR  = NDET
C
C     transfer to LUCIPARREL
      IRC_SAVE = IRC
C
!     define whether it is a fresh CI (== 0) or a restart (== 1)
      istart_nroot =  1 * nroot
      IF(IRESTR.EQ.0) THEN
        INICI             = 0
        istart_nroot =  min(6 * nroot,MXCIV_CI-nroot,ndet)
!old SK istart_nroot =  nroot
      ELSE
        INICI             = -1
      END IF

      !> for debugging 
!     istart_nroot = 3

C
#if defined (VAR_MPI2)
      IF( IT_SHL .ge. 0 .or. SPLIT_IJKL )THEN
C
C       allocate T-block processor and offset list
        MY_TTPL_LEN = NSPOBEX_TP * IRC
        MY_TTOL_LEN = NMPROC
#if defined VAR_SHMEM
        MY_TTOL_LEN = NSPOBEX_TP * IRC
        CALL MPIXMEM_ALLOC( MY_TTPL_PTR, 2, MY_TTPL_LEN,
     &                      df_MPI_INFO_NULL, .FALSE.)
        CALL MPIXMEM_ALLOC( MY_TTOL_PTR, 3, MY_TTOL_LEN,
     &                      df_MPI_INFO_NULL, .FALSE.)
        CALL IMINO8(IT_TTOL,MY_TTOL_LEN)
#else
        allocate(IT_TTPL(MY_TTPL_LEN))
        allocate(IT_TTOL(MY_TTOL_LEN))
        IT_TTPL =  0
        IT_TTOL = -1
#endif
      END IF
#endif
C
C     partition CI vector - for all MS2 spaces into batches
      CALL Z_BLKFO_REL_1(IDC,NMS2VAL,ICSPC,ICSM,KPCIBT,KPCBLTP,
     &                   KNODELIST,KBLOCKLIST,KRCCTOS,
     &                   NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,
     &                   IBAT_MS2,1,0,IDUMMY,IRC,IT_TTPL,
     &                   IGROUPLIST)
C
#if defined (VAR_MPI2)
C
C     define length of T_DIAG block and store on common block 
C     variable MY_T_LEN_D
      LEN_T_D_BUFF = (3*(NTOOB**2)) + (2*NTOOB)
      MY_T_LEN_D   = LEN_T_D_BUFF
C     
C     INFO about T block distribution on IT_TTPL
C     ==> if (shared memory route active) then
C            allocate T buff, (read (ij|kl) in), open T window
C         end if
C     <==
C
      IF( IT_SHL .ge. 0 .or. SPLIT_IJKL )THEN
C
C       calculate length of T_BUFF
        MY_T_LEN     = IGET_MY_T_LEN(IT_TTPL,NSPOBEX_TP,
     &                               N1ELINT,IRC)
        LEN_T_BUFF   = MY_T_LEN
        IF( SPLIT_IJKL ) THEN
C         real part of T array
          LEN_T_BUFF = MY_T_LEN/IRC
        END IF
csk     write(luwrt,*) ' LEN_T_BUFF, MY_T_LEN, MY_T_LEN_D',
csk  &                   LEN_T_BUFF, MY_T_LEN, MY_T_LEN_D
csk     CALL IWRTMAMN(IT_TTPL,1,NSPOBEX_TP,1,NSPOBEX_TP,LUWRT)
C
C       FIXME: provide information if available
C
        IS_DUMMY = 0
C
C       allocate T_BUFF array (replaces use of WORK(KT_CC))
#if defined VAR_SHMEM
        CALL GET_MEM_T( IS_DUMMY )
        CALL GET_MEM_T_D( IS_DUMMY )
C       zero arrays
        CALL DZERO8(T_BUFF,MY_T_LEN)
        CALL DZERO8(T_BUFF_D,MY_T_LEN_D)
#else
        allocate(T_BUFF(MY_T_LEN))
        allocate(T_BUFF_D(MY_T_LEN_D))
        DIAGONAL_ALLOC = .TRUE.
        T_BUFF   =  0
        T_BUFF_D =  0
#endif
C
#ifdef INT_STAR8
        IF( SPLIT_IJKL ) THEN
         call interface_MPI_GATHER(MY_T_LEN,1,IT_TTOL,1,MASTER,
     &                             global_communicator)
        END IF
#else
        call quit('gasci main driver: the chosen integral split '//
     &            'routine requires 64-bit integer in dirac')
#endif
C
        WRITE(LUWRT,'(/2X,A)')
     & '=============================================================='
        WRITE(LUWRT,'(2X,A/)')
     & '==>  allocation of process unique T (integral) elements    <=='
        DO II = 1, NMPROC
        WRITE(LUWRT,'(2X,A,1X,I4,1X,A,I12,1X,A)')
     & 'process',II-1,'allocates for ',IT_TTOL(II),'elements          '
        END DO
      
        WRITE(LUWRT,'(2X,A/)')
     & '=============================================================='

        IF(SPLIT_IJKL) deallocate(IT_TTOL)
C
C       read in all (active) blocked integrals (T blocks)
C       ==> integrals required for diagonal H on T_BUFF_D
C
        IF (CIRUN /= 'ANALYZ' .and. CIRUN /= 'QCORR ' .and. 
     &      CIRUN /= 'REFVEC') THEN
C           FIXME: 11-02-08 SK: definition of IBU2D?
          call picasso(iprtra,.TRUE.,
     &                 NL2D,IBU2D,ISM2D,F1,F2,T_BUFF,
     &                 T_BUFF_D,IT_TTPL,ibosym,CIRUN,.FALSE.)
          CPIC_PASS_D = .TRUE.
        END IF
C
csk     CALL WRTMATMN8(T_BUFF,1,MY_T_LEN,1,MY_T_LEN,LUWRT)
csk     CALL WRTMATMN8(T_BUFF_D,1,MY_T_LEN_D,1,MY_T_LEN_D,LUWRT)
C
#if defined VAR_SHMEM
        IF( IT_SHL .ge. 0 )THEN
C
C         open memory window (collective call for communicator
C         group MY_NEWCOMM_SM)
C
          CALL GET_T_WIN( T_BUFF, MYNEW_COMM_SM, MY_T_LEN )
C
          CALL CALC_WIN_OFFSET(IT_TTPL,IT_TTOL,NSPOBEX_TP,IRC)
C
csk        WRITE(LUWRT,*) ' printing offset array IT_TTOL'
csk        CALL IWRTMAMN8(IT_TTOL,1,MY_TTOL_LEN,1,MY_TTOL_LEN,LUWRT)
        END IF
#endif
C
C       transfer IS_LENGTH_TT to common block LUCIPAR_I8
C       == max. T block length
C
        IS_LENGTH_TT = 0
        IS_LENGTH_TT = NCALC_MAX_BLK(IRC)
C
      END IF
C     ... read-in reordered integrals for the old-school way.
      IF( .NOT. CPIC_PASS_D .and. REORD_IJKL ) THEN
C
         IF( .NOT. DIAGONAL_ALLOC ) allocate(T_BUFF_D(MY_T_LEN_D))
         DIAGONAL_ALLOC = .TRUE. 
C
CSK      CALL DZERO8(T_BUFF,MY_T_LEN)
         IF (CIRUN /= 'ANALYZ' .and. CIRUN /= 'QCORR ' .and. 
     &      CIRUN /= 'REFVEC') THEN
           call picasso(iprtra, .TRUE. ,
     &                  NL2D,IBU2D,ISM2D,F1,F2,T_BUFF,
     &                  T_BUFF_D,IT_TTPL,ibosym,CIRUN,.FALSE.)
           CPIC_PASS_D = .TRUE.
         END IF
      END IF
#endif
C
      CALL Z_BLKFO_REL_2(IDC,NMS2VAL,ICSPC,ICSM,KPCIBT,KPCBLTP,
     &                   NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,
     &                   IBAT_MS2,1,0,IDUMMY,IRC)
C
#if defined (VAR_MPI2)
C
      IF ( CIRUN .eq. 'DIAG  ' .or. CIRUN .eq. 'CIINII' .or.
     &     CIRUN .eq. 'KR-CI ' .or.
     &     CIRUN .eq. 'SIGMA ' .or. CIRUN .eq. 'SIGMAD' .or. 
     &     CIRUN .eq. 'DENS1 ' .or. CIRUN .eq. 'DENS2 ' ) THEN
C
C       dummy file array allocation for 'remaining' CPUs
        IF( IAM_NOT_INV .eq. 0 ) THEN
          IALL_LU1 = 1
          IALL_LU2 = 1
          IALL_LU3 = 1
          IALL_LU4 = 1
          IALL_LU5 = 1
          IALL_LUC = NUM_BLOCKS2
        END IF
C
        CALL CALC_OFF_MPI_FILE('ILUC  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                         WORK(KNODELIST),1)
C
        CALL CALC_OFF_MPI_FILE('ILU1  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                         WORK(KNODELIST),MXCIV_CI+istart_nroot)
C
        CALL CALC_OFF_MPI_FILE('ILU2  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                         WORK(KNODELIST),MXCIV_CI+istart_nroot)
C
C       initialize general offset for ILUC
        MY_LUC_OFF = 0
C
C       length for allocation of file arrays
        IALL_LU1 = ( MXCIV_CI + istart_nroot ) * MY_ACT_BLK2
        IALL_LU2 = ( MXCIV_CI + istart_nroot ) * MY_ACT_BLK2
        IALL_LUC = NUM_BLOCKS2
C
C       allocate file arrays
        CALL MEMMAR(KILU1LIST,IALL_LU1,'ADDS  ',1,'LU1LST')
        CALL MEMMAR(KILU2LIST,IALL_LU2,'ADDS  ',1,'LU2LST')
        CALL MEMMAR(KILUCLIST,IALL_LUC,'ADDS  ',1,'LUCLST')
C       initialize ...
        CALL IZERO(WORK(KILU1LIST),IALL_LU1)
        CALL IZERO(WORK(KILU2LIST),IALL_LU2)
        CALL IZERO(WORK(KILUCLIST),IALL_LUC)
C
        IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'SIGMAD' .or. 
     &      CIRUN .eq. 'DIAG  ' .or. CIRUN .eq. 'KR-CI ' ) THEN 
C
          CALL CALC_OFF_MPI_FILE('LUDIA ',IGROUPLIST,WORK(KBLOCKLIST),
     &                           WORK(KNODELIST),1)
C
          IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' ) THEN
C
            CALL CALC_OFF_MPI_FILE('ILU3  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                             WORK(KNODELIST),NROOT)
C
            CALL CALC_OFF_MPI_FILE('ILU4  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                             WORK(KNODELIST),1)
C
            CALL CALC_OFF_MPI_FILE('ILU5  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                             WORK(KNODELIST),1)
C
C           some more file arays ... 
            IALL_LU3 = ( NROOT ) * MY_ACT_BLK2
            IALL_LU4 = ( 1 )     * MY_ACT_BLK2
            IALL_LU5 = ( 1 )     * MY_ACT_BLK2
C           allocate ...
            CALL MEMMAR(KILU3LIST,IALL_LU3,'ADDS  ',1,'LU3LST')
            CALL MEMMAR(KILU4LIST,IALL_LU4,'ADDS  ',1,'LU4LST')
            CALL MEMMAR(KILU5LIST,IALL_LU5,'ADDS  ',1,'LU5LST')
C           initialize ...
            CALL IZERO(WORK(KILU3LIST),IALL_LU3)
            CALL IZERO(WORK(KILU4LIST),IALL_LU4)
            CALL IZERO(WORK(KILU5LIST),IALL_LU5)
C
          END IF
C         ^ CIRUN == CIINII or KR-CI
        END IF            
C       ^ CIRUN == CIINII .or. SIGMAD .or. DIAG .or. KR-CI
      END IF
C     ^ CIRUN == DIAG .or. CIRUN == CIINII .or. CIRUN == SIGMA ...
#endif
C
C     allocate memory : 2 vector blocks 
C     store veclen for  1 vector/scratch block
      CALL GET_3BLKS_REL(KVEC1,KVEC2,LSCR2,CIRUN)
C
C     CI diagonal - if required
*+++++++++++++++++++++++++++++++++++++++++++++++++++
      if (CIRUN.eq.'DIAG  '.or.CIRUN.eq.'CIINII'.or.
     &    CIRUN.eq.'KR-CI '.or.
     &    CIRUN.eq.'SIGMA '.or.CIRUN.eq.'SIGMAD') then
*
        if (CIRUN.ne.'SIGMA ') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
*
      IF(IDIAG.EQ.2) THEN
        LUDIA = LUSC1
      END IF
#if defined (VAR_MPI2)
      IF( USE_EX_IDIA ) THEN
        WRITE(LUWRT,*) '  saving CPU/WALL-time in this MPI run   '
        WRITE(LUWRT,*) '  using existing diagonal file DIAPAR.1. '
        WRITE(LUWRT,*) '                                         '
        GOTO 55
      END IF
#endif
      IF(.NOT.(IDIAG.EQ.2.AND.IRESTR.EQ.1)) THEN
        ECOREP = ECORE
C
        IF(ICISTR.GE.2) CALL REWINE(LUDIA,-1)
        IRECO = 1
        CALL GASDIAT_REL(WORK(KVEC1),LUDIA,ECOREP,ICISTR,IHAM12,
     &               WORK(KPCBLTP),NBLOCK,WORK(KPCIBT),
     &               T_BUFF,T_BUFF_D,IDENSI,
     &               nirr_dg,ISPINFREE,
     &               ENVIRO,SIGDEN_ROUTE,IRECO,WORK(KNODELIST),
     &               CPIC_PASS_D)
C       end of real part
        CALL ITODS(-1,1,-1,LUDIA)
        IF (IRC.EQ.2) THEN
C         diagonal once again- for the imaginary part
          IRECO = 2
          CALL GASDIAT_REL(WORK(KVEC1),LUDIA,ECOREP,ICISTR,IHAM12,
     &                 WORK(KPCBLTP),NBLOCK,WORK(KPCIBT),
     &                 T_BUFF,T_BUFF_D,IDENSI,
     &                 nirr_dg,ISPINFREE,
     &                 ENVIRO,SIGDEN_ROUTE,IRECO,WORK(KNODELIST),
     &                 CPIC_PASS_D)
C       end of imaginary  part
        CALL ITODS(-1,1,-1,LUDIA)
        END IF
        IF(IPRCIX.GE.2) WRITE(LUWRT,*) ' Diagonal constructed  '
      ELSE
         WRITE(LUWRT,*) ' Diagonal not calculated '
      END IF
c
#if defined (VAR_MPI2)
C
 55   CONTINUE
C
      NPTEST_VAR = 00
      IF( NPTEST_VAR .ge. 10 ) THEN
      IRECO = IRC_SAVE
      ITEST_OFF = 0
      DO ILPRECO = 1, IRECO
        ITEST_OFF = MY_DIA_OFF + MY_VEC1_IOFF * ( ILPRECO - 1 )
        DO IIBLK = 1, NUM_BLOCKS
          ILEN = 0
          CALL GET_BLOCK_PROC(WORK(KNODELIST),IIBLK,IPROC)
          IF( MYPROC .eq. IPROC ) THEN
            CALL GET_BLOCK_LENGTH(WORK(KBLOCKLIST),IIBLK,ILEN)
            call interface_mpi_FILE_READ_AT(IDIA,ITEST_OFF,
     &                            WORK(KVEC1),ILEN,ISTAT)
            WRITE(LUWRT,*) ' Diagonal elements',ILEN
            CALL WRTMATMN(WORK(KVEC1),1,ILEN,1,ILEN,LUWRT)
            ITEST_OFF = ITEST_OFF + ILEN
          END IF
        END DO
      END DO
      
      ENDIF
      NPTEST_VAR = 00
#endif
c
*+++++++++++++++++++++++++++++++++++++++++++++++++++
        end if
        if (CIRUN.ne.'DIAG  ') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     transfer control to optimization routine
c
      MINST  = 1
      NPRDET = 0
C
C     restart from previous vector(s)
C
#if defined (VAR_MPI2)
C
      IF( INICI .lt. 0 )THEN
        if ( CIRUN.eq.'SIGMAD'.or.CIRUN.eq.'SIGMA ') GOTO 88 
          WRITE(LUWRT,'(/A,A)') '  MPI restart file will be'//
     &    ' prepared from C-vector file: KRCI_CVECS.',SYMFLABEL
C
         IF( MYPROC .eq. MASTER ) CALL REWINE(LUC,-1)
C
         CALL COP_REST_VEC_REL(WORK(KVEC1),LUC,ILU1,MY_LU1_OFF,
     &                         WORK(KILU1LIST),
     &                         WORK(KNODELIST),WORK(KBLOCKLIST),
     &                         NBLOCK,NROOT,IRC)
C
        WRITE(LUWRT,'(A/)') '  file conversion from standard I/O'//
     &                      ' to MPI-IO finished.'
C
      END IF
 88   CONTINUE
C
#endif
C
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
          if (CIRUN.eq.'SIGMAD'.or.CIRUN.eq.'SIGMA ') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
C
C         C-vector resides on MASTER in file LUC
C
C                        ============
C                        PARALLEL RUN
C                        ============
C
C         NODES recieve the appropriate C-blocks
C
C         MASTER recieves the corresponding SIGMA-blocks to 
C         write on LUHC
C
C         IF( SIGMAD )   MASTER recieves the corresponding diagonal
C                        elements to write on LUDIA
C
#if defined (VAR_MPI2)
C
            IF( MYPROC .eq. MASTER ) CALL REWINE(LUC,-1)
C
            CALL COP_REST_VEC_REL(WORK(KVEC1),LUC,ILU1,MY_LU1_OFF,
     &                            WORK(KILU1LIST),
     &                            WORK(KNODELIST),WORK(KBLOCKLIST),
     &                            NBLOCK,IIONE,IRC)
C
C            Arrays for partitioning of the CI vector
C            (*_MS2 arrays are in CSTATE )
             call z_blkfo_rel_par(IDC,NMS2VAL,ISPC,ICSM,
     &                            KLVLBT,KLVLEBT,KLVI1BT,KLVIBT,KDUM,
     &                            NBATCH_LD,NM_BLK_LD,NBLK_MS2_LD,
     &                            IBLK_MS2_LD,NBAT_MS2_LD,IBAT_MS2_LD,
     &                            NZERO,IIONE,WORK(KNODELIST))
C
C           set offset for sigma-file
            JVEC_SF = 0
C
C           erase old LUCLIST
            CALL IZERO(WORK(KILUCLIST),IALL_LUC)
C
            IF( IRC .eq. 1 ) THEN
C             real quaternion double group
              CALL COPVCD_PP_CC_B(ILU1,ILUC,WORK(KVEC1),NBATCH_LD,
     &                            WORK(KLVLBT),WORK(KLVLEBT),
     &                            WORK(KLVI1BT),WORK(KLVIBT),
     &                            MY_LU1_OFF,MY_LUC_OFF,
     &                            WORK(KILU1LIST),WORK(KILUCLIST),
     &                            WORK(KBLOCKLIST),NZERO)
            ELSE
C             complex quaternion double group
              CALL COPVCD_PP_CC_B_C(ILU1,ILUC,WORK(KVEC1),NBATCH_LD,
     &                              WORK(KLVLBT),WORK(KLVLEBT),
     &                              WORK(KLVI1BT),WORK(KLVIBT),
     &                              MY_LU1_OFF,MY_LUC_OFF,
     &                              WORK(KILU1LIST),WORK(KILUCLIST),
     &                              WORK(KBLOCKLIST),NZERO)
            END IF
C
#endif
C           start calculation: sigma = H x C
C
************************************************************************
            call sigden_ctrl(WORK(KVEC1),WORK(KVEC2),
#if defined (VAR_MPI2)
     &                       ILUC,ILU2,
#else
     &                       LUC,LUHC,
#endif
     &                       WORK(KT_CC),T_BUFF,1
#if defined (VAR_MPI2)
     &                      ,WORK(KILUCLIST),WORK(KILU2LIST),
     &                       WORK(KBLOCKLIST),WORK(KNODELIST),
     &                       IGROUPLIST,IPROCLIST,WORK(KRCCTOS),
     &                       IT_TTPL,IT_TTOL
     &                      ,WORK(KBLOCKLIST),WORK(KNODELIST)
#endif
     &                       )
************************************************************************
C
#if defined (VAR_MPI2)
C
            IF( IAM_NOT_INV .eq. 0 ) GOTO 405
C
            IF( MYPROC .eq. MASTER )  CALL REWINE(LUHC,-1)
C
            LBLK = -1
C
            CALL COPVCD_PAR_BDRIV_REL(ILU2,LUHC,WORK(KVEC1),
     &                                WORK(KNODELIST),NUM_BLOCKS,
     &                                WORK(KBLOCKLIST),
     &                                global_communicator,LBLK,IIONE,
     &                                WORK(KILU2LIST),MY_LU2_OFF,IRC)
C
C
            IF( CIRUN .eq. 'SIGMAD' ) THEN
C
C             save H diagonal on MASTER
C
              IF( MYPROC .eq. MASTER ) CALL REWINE(LUDIA,LBLK)
C             H is a hermitian matrix...
              CALL COPVCD_PAR_BDRIV_DIA_REL2(IDIA,LUDIA,WORK(KVEC1),
     &                                       WORK(KNODELIST),NUM_BLOCKS,
     &                                       WORK(KBLOCKLIST),
     &                                       global_communicator,
     &                                       LBLK,IIONE,MY_DIA_OFF)
            END IF
C           ^ CIRUN .eq. SIGMAD
 405    CONTINUE
C
#endif
          else
 
            if (NROOT.gt.MXROOT) then
              WRITE(LUWRT,'(A,I4)') ' Too many roots. MXROOT is ',MXROOT
              WRITE(LUWRT,'(A,I4)') ' Increase MXROOT to ',NROOT
              CALL QUIT(' *** ERROR in GASCI_REL: too many roots ***')
            end if
C
C           !> allocate some scratch space
            NFINDM = min(ndet,istart_nroot+3*nroot)

!           allocate(scratch_dav1(4*NFINDM+1),iscratch_dav1(5*NFINDM))
!           allocate(scratch_dav2(NROOT))
            call alloc(scratch_dav1, 4*NFINDM+1)
            call alloc(iscratch_dav1,5*NFINDM)
            call alloc(scratch_dav2, NROOT)
            call alloc(eci_start,istart_nroot,id='eci_start')
C
            CALL CIEIG5_REL(INICI,EROOT,WORK(KVEC1),WORK(KVEC2),
     &                      MINST,LUDIA,LUC,LUHC,LUSC1,LUSC2,LUSC3,
     &                      LUSC34,LUSC35,NVAR,NBLK,NROOT,MXCIV_CI,
     &                      MAXIT,MXROOT,LUCIVI,istart_nroot,IPRNT,
     &                      NPRDET,pert_dummy,MXP1,MXP2,MXQ,
     &                      scratch_dav1,iscratch_dav1,scratch_dav2,
     &                      WORK(KT_CC),T_BUFF,ECORE,ICISTR,LBLK,
     &                      IDIAG,THRES_E,NBATCH,WORK(KPCLBT),
     &                      WORK(KPCLEBT),WORK(KPCIBT),INIDEG,
     &                      E_THRE,C_THRE,E_CONV,C_CONV,ICLSSEL,NOCCLS,
     &                      IRC,IDCOMH,IHAM12,
     &                      WORK(KBLOCKLIST),WORK(KNODELIST),
     &                      IPROCLIST,IGROUPLIST,eci_start,ITERSEOUT
#if defined (VAR_MPI2)
     &                      ,WORK(KILU1LIST),WORK(KILU2LIST),
     &                      WORK(KILU3LIST),WORK(KILU4LIST),
     &                      WORK(KILU5LIST),WORK(KILUCLIST),
     &                      WORK(KRCCTOS),IT_TTPL,IT_TTOL
#endif
     &                      )
C
C           save energy of the lowest eigenstate
            EREF = EROOT(1)
C
C           release local scratch memory
            call dealloc(eci_start)
            call dealloc(scratch_dav2)
            call dealloc(iscratch_dav1)
            call dealloc(scratch_dav1)
!           deallocate(scratch_dav2,iscratch_dav1,scratch_dav1)
C
#if defined (VAR_MPI2)
C
C           delete big scratch files when doing large-scale CI run
C           (CIRUN .eq. KR-CI) before collecting all vectors
            CALL SETUNITS_PAR_CLOSE_REL(CIRUN,1)
C
            IF( IAM_NOT_INV .eq. 0 ) GOTO 505
C
            LBLK = -1
C           copy c-vectors from nodes and master back to the master
C           purpose: master does the ana-stuff
C
C           NOTE: only C-vectors will be saved.
C
            IF( MYPROC .eq. MASTER ) CALL REWINE(LUC,-1)
            CALL REWINE(LUSC41,-1)
            DO JROOT = 1, NROOT
              CALL COPVCD_PAR_BDRIV_REL(ILU1,LUSC41,WORK(KVEC1),
     &                                  WORK(KNODELIST),NUM_BLOCKS,
     &                                  WORK(KBLOCKLIST),
     &                                  global_communicator,LBLK,JROOT,
     &                                  WORK(KILU1LIST),MY_LU1_OFF,IRC)
            END DO
            CALL SETUNITS_PAR_CLOSE_REL(CIRUN,2)
            IF(MYPROC.EQ.MASTER) THEN
              CALL REWINE(LUSC41,-1)
              CALL REWINE(LUC,-1)
              DO JROOT = 1, NROOT
                CALL copvcdc(LUSC41,LUC,WORK(KVEC1),0,IRC,LBLK)
              ENDDO
              CALL REWINE(LUC,-1)
            END IF
            NPTESTVAR = 00
            IF( NPTESTVAR .ge. 10 )THEN
              IF(MYPROC.EQ.MASTER) THEN
                DO IVEC = 1, NROOT
                  WRITE(LUWRT,*) '  final new C vector',IVEC
                  CALL WRTVCDC(WORK(KVEC1),LUC,0,IRC,LBLK)
                END DO
              END IF
            END IF
            NPTESTVAR = 00
C
 505        CONTINUE
            IF( IAM_NOT_INV .eq. 0 ) 
     &      CALL SETUNITS_PAR_CLOSE_REL(CIRUN,2)
#endif
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
          end if
C         ^ CIRUN ne SIGMAD or SIGMA
        end if
C       ^ CIRUN ne DIAG
      else if (CIRUN.eq.'ANALYZ') then
C     analysis of optimized wave function
*+++++++++++++++++++++++++++++++++++++++++++++++++++
C
#if defined (VAR_MPI2)
C
C     only MASTER analyzes the optimized wave function ...
C
      IF( MYPROC .ne. MASTER ) GOTO 888
C
#endif
      LBLK = -1
      call REWINE(LUC,-1)
      NPTESTVAR = 00
      IF( NPTESTVAR .ge. 10 )THEN
        IF(MYPROC.EQ.MASTER) THEN
          DO IVEC = 1, NROOT
            WRITE(LUWRT,*) '  final new C vector before ana',IVEC
            CALL WRTVCDC(WORK(KVEC1),LUC,0,IRC,LBLK)
          END DO
        END IF
      END IF
      NPTESTVAR = 00
      call REWINE(LUC,-1)
      do JROOT=1,NROOT,1
        call REWINE(LUSC61,-1)
        call copvcdc(LUC,LUSC61,WORK(KVEC2),0,IRC,LBLK)
        call REWINE(LUSC61,-1)
        WRITE(LUWRT,'(/A,I3)')
     &  ' ***************************************************'
        WRITE(LUWRT,'(A,I3)')
     &  ' Analysis of Optimized Wave Function for ROOT = ',JROOT
        WRITE(LUWRT,'(A,I3)')
     &  ' ***************************************************'
C
        call gasana_rel(WORK(KVEC1),NBLOCK,WORK(KPCIBT),
     &                  WORK(KPCBLTP),LUSC61,LUSC37,
     &                  ICISTR,IRC)
      end do
      call REWINE(LUC,-1) 
      NPTESTVAR = 00
      IF( NPTESTVAR .ge. 10 )THEN
        IF(MYPROC.EQ.MASTER) THEN
          DO IVEC = 1, NROOT
            WRITE(LUWRT,*) '  final new C vector after ana',IVEC
            CALL WRTVCDC(WORK(KVEC1),LUC,0,IRC,LBLK)
          END DO
        END IF
      END IF
      NPTESTVAR = 00

 888  CONTINUE 
#if defined (VAR_MPI2)
C     let the slaves wait... it is only a short way to go so why not?
      call interface_mpi_BARRIER(global_communicator)
#endif

      else if (CIRUN.eq.'QCORR ') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     only MASTER computes the +Q correction ...
C
        if(myproc == master)then
          call rewine(LUC,-1)
          call q_correction(WORK(KVEC1),NBLOCK,WORK(KPCIBT),
     &                      WORK(KPCBLTP),LUC,LUSC37,
     &                      ICISTR,IRC,WORK,.true.)
        end if

#if defined (VAR_MPI2)
C     let the slaves wait... it is only a short way to go so why not?
      call interface_mpi_BARRIER(global_communicator)
#endif
C
      else if (CIRUN.eq.'REFVEC') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     only MASTER writes the reference vector data to file
C
        if(myproc == master)then
          call rewine(LUC,-1)
          call put_refvec_for_qcorr(WORK(KVEC1),NBLOCK,WORK(KPCIBT),
     &                              WORK(KPCBLTP),LUC,LUSC37,
     &                              ICISTR,IRC,WORK,.true.,ndet,eref)
        end if

#if defined (VAR_MPI2)
C     let the slaves wait... it is only a short way to go so why not?
      call interface_mpi_BARRIER(global_communicator)
#endif
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
      else if (CIRUN.eq.'DENS1 '.or.CIRUN.eq.'DENS2 ') then
C
C     analysis of density and occupation
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     real or complex eigenvectors
      IF(ispinfree.EQ.1) THEN
C       real
        IRC = 1
      ELSE
C       well, complex
        IRC = 2
        if (NZ.eq.1) then
C         quaternion real
          IRC = 1
        end if
      END IF
C
      LBLK = -1
C
      if (IDENSI.ne.0) then
        if (IDENSI.eq.1.or.IDENSI.eq.2) then
          write(LUWRT,*) 'IDENSI is ',IDENSI
          write(LUWRT,*) 
     &    'I am not allowed to calculate non-relativistic'
          write(LUWRT,*) ' density matrices YET.'
          write(LUWRT,*) ' Please ask the experts for instructions'
          Call Abend2( 'LUCIAREL quitting' )
        end if
#if defined (VAR_MPI2)
C
C       arrays for partitioning of the CI vector
C       (*_MS2 arrays are in CSTATE )
        call z_blkfo_rel_par(IDC,NMS2VAL,ISPC,ICSM,
     &                       KLVLBT,KLVLEBT,KLVI1BT,KLVIBT,KDUM,
     &                       NBATCH_LD,NM_BLK_LD,NBLK_MS2_LD,
     &                       IBLK_MS2_LD,NBAT_MS2_LD,IBAT_MS2_LD,
     &                       0,1,WORK(KNODELIST))
C
#endif
C       
C       calculation of densities for natural orbital occ. numbers
C       ---------------------------------------------------------
        IF( NOOCCN_DENS_RUN ) THEN
C
          IF( MYPROC .eq. MASTER ) THEN
            CALL REWINE(LUSC61,LBLK)
            CALL REWINE(LUSC62,LBLK)
            CALL SKPVCD(LUC,IRC * (IDENSI_STATE - 1),WORK(KVEC1),1,LBLK)
          END IF 
#if defined (VAR_MPI2)
          IVEC_TMP = 1
C         right hand-side vector
          CALL COP_REST_VEC_REL(WORK(KVEC1),LUC,ILU1,MY_LU1_OFF,
     &                          WORK(KILU1LIST),
     &                          WORK(KNODELIST),WORK(KBLOCKLIST),
     &                          NUM_BLOCKS,IVEC_TMP,IRC)
C          
          IF( MYPROC .eq. MASTER) THEN
            CALL SKPVCD(LUC,IRC * (IDENSI_STATE - 1),WORK(KVEC1),1,LBLK)
          END IF
C         left hand-side vector
          CALL COP_REST_VEC_REL(WORK(KVEC1),LUC,ILU2,MY_LU2_OFF,
     &                          WORK(KILU2LIST),
     &                          WORK(KNODELIST),WORK(KBLOCKLIST),
     &                          NUM_BLOCKS,IVEC_TMP,IRC)
C       files are in place now... say hello
        WRITE(LUWRT,*) ' '
        WRITE(LUWRT,'(12X,A)') ' ==================================== '
        WRITE(LUWRT,'(12X,A)') ' PARALLEL 1p DENSITY MAT. CALCULATION '
        WRITE(LUWRT,'(12X,A)') ' ==================================== '
        WRITE(LUWRT,*) ' '
C
C
C
#else
C         right hand-side vector to LUSC61
          CALL COPVCDC(LUC,LUSC61,WORK(KVEC1),0,IRC,LBLK)
C
          CALL SKPVCD(LUC,IRC * (IDENSI_STATE - 1),WORK(KVEC1),1,LBLK)
C         left hand-side vector on LUSC62
          CALL COPVCDC(LUC,LUSC62,WORK(KVEC1),0,IRC,LBLK)
#endif
C
          GOTO 895
C
        END IF ! natural orbital occupation numbers
C
        if(is_ana_mcscf)then
!         right hand-side vector to LUSC61
          CALL COPVCDC(LUC,LUSC61,WORK(KVEC1),1,IRC,LBLK)
!         left hand-side vector on LUSC62
          CALL COPVCDC(LUC,LUSC62,WORK(KVEC1),1,IRC,LBLK)          
        end if

C       prepare vector files for pure dens. matrix calc. run (parallel)
C       ---------------------------------------------------------------
#if defined (VAR_MPI2)
        IF( MYPROC .eq. MASTER ) THEN
          CALL REWINE(LUSC61,LBLK)
          CALL REWINE(LUSC62,LBLK)
        END IF
C       right hand-side vector on LUSC61
        CALL COP_REST_VEC_REL(WORK(KVEC1),LUSC61,ILU1,MY_LU1_OFF,
     &                        WORK(KILU1LIST),
     &                        WORK(KNODELIST),WORK(KBLOCKLIST),
     &                        NUM_BLOCKS,1,IRC)
C       left hand-side vector on LUSC62
        CALL COP_REST_VEC_REL(WORK(KVEC1),LUSC62,ILU2,MY_LU2_OFF,
     &                        WORK(KILU2LIST),
     &                        WORK(KNODELIST),WORK(KBLOCKLIST),
     &                        NUM_BLOCKS,1,IRC)
C
C       files are in place now... say hello
        WRITE(LUWRT,*) ' '
        WRITE(LUWRT,'(12X,A)') ' =================================== '
        WRITE(LUWRT,'(12X,A)') ' PARALLEL DENSITY MATRIX CALCULATION '
        WRITE(LUWRT,'(12X,A)') ' =================================== '
        WRITE(LUWRT,*) ' '
C
#endif
C
 895    CONTINUE
C
C       one-body density matrix
C       spin-dependent first-order density matrix
C
        JROOT = 1
        JROOT_SCR = 0
C
        IF( NOOCCN_DENS_RUN ) THEN
           JROOT_SCR = IDENSI_STATE 
        ELSE
           JROOT_SCR = 1
        END IF
C
        if (IDENSI.eq.3) then
          IDEN12 = 1
          write(LUWRT,*) ' '
          write(LUWRT,'(2X,A)') 
     &    ' Computing spin-dependent one-particle density'
          write(LUWRT,*) ' '
        else if (IDENSI.eq.4) then
          IDEN12 = 2
          write(LUWRT,*) ' '
          write(LUWRT,'(2X,A)') 
     &   ' Computing spin-dependent one- and two-particle densities.'
          write(LUWRT,*) ' '
        else
          write(LUWRT,*) 'IDENSI =  ',IDENSI
          Call Abend2( 'Problemo in gasci. Quitting.' )
        end if
C
C
!       IPRDEN_save = IPRDEN
!       IPRDEN      = 500
        if (IPRDEN .ge. 500 .and. .not. NOOCCN_DENS_RUN) then
          IF( MYPROC .eq. MASTER )THEN
            write(LUWRT,*) 'Current vector',JROOT_SCR,' on unit',LUSC61
            call wrtvcdc(WORK(KVEC1),LUSC61,1,IRC,LBLK)
            write(LUWRT,*) 'Current vector ',JROOT_SCR,' on unit',LUSC62
            call wrtvcdc(WORK(KVEC2),LUSC62,1,IRC,LBLK)
          END IF
        end if
!       IPRDEN = IPRDEN_save
C
        NTERMS = NALLINT
        if (CIRUN.eq.'DENS1 ') NTERMS = N1ELINT
C
C       erase T amplitudes from sigma calculation
C
        IF( IT_SHL .ge. 0 .or. SPLIT_IJKL )THEN
          call DZERO8(T_BUFF,MY_T_LEN)
        ELSE
          call DZERO8(WORK(KT_CC),NALLINT)
        END IF
C
!      WRITE(LUWRT,*) ' start: MY T (BUFF) array : ',NALLINT
!      IF( SPLIT_IJKL )THEN
!        CALL WRTMATMN8(T_BUFF,1,MY_T_LEN,1,MY_T_LEN,LUWRT)
!      ELSE
!        CALL WRTMATMN8(WORK(KT_CC),1,NALLINT,1,NALLINT,LUWRT)
!      END IF
C
#if defined (VAR_MPI2)
C       erase old LUCLIST
        CALL IZERO(WORK(KILUCLIST),IALL_LUC)
C       copy right-hand vector to scratch file
        IF( IRC .eq. 1 ) THEN
C         real quaternion double group
          CALL COPVCD_PP_CC_B(ILU1,ILUC,WORK(KVEC1),NBATCH_LD,
     &                        WORK(KLVLBT),
     &                        WORK(KLVLEBT),WORK(KLVI1BT),
     &                        WORK(KLVIBT),MY_LU1_OFF,MY_LUC_OFF,
     &                        WORK(KILU1LIST),WORK(KILUCLIST),
     &                        WORK(KBLOCKLIST), JROOT - 1 )
        ELSE
C         complex quaternion double group
          CALL COPVCD_PP_CC_B_C(ILU1,ILUC,WORK(KVEC1),NBATCH_LD,
     &                          WORK(KLVLBT),WORK(KLVLEBT),
     &                          WORK(KLVI1BT),WORK(KLVIBT),
     &                          MY_LU1_OFF,MY_LUC_OFF,
     &                          WORK(KILU1LIST),WORK(KILUCLIST),
     &                          WORK(KBLOCKLIST), JROOT - 1 )
        END IF
C
#endif
C
C       set offset for left-hand vector on ILU2 
C       (parallel calculation)
C
        JVEC_SF = JROOT - 1
C
C       compute density matrix elements
C       -------------------------------
C
*=======================================================================
        call sigden_ctrl(WORK(KVEC1),WORK(KVEC2),
#if defined (VAR_MPI2)
     &                   ILUC,ILU2,
#else
     &                   LUSC61,LUSC62,
#endif
     &                   WORK(KT_CC),T_BUFF,2
#if defined (VAR_MPI2)
     &                   ,WORK(KILUCLIST),WORK(KILU2LIST),
     &                   WORK(KBLOCKLIST),WORK(KNODELIST),
     &                   IGROUPLIST,IPROCLIST,WORK(KRCCTOS),
     &                   IT_TTPL,IT_TTOL
     &                   ,WORK(KBLOCKLIST),WORK(KNODELIST)
#endif
     &                   )
*=======================================================================
C
C         file INDICES.x
C
C         IF( .NOT. SHARED_M ) THEN
C           read indices from cpu-specific file INDICES.'MYPROC'
C           ( written in dist_ints_sigden --> interface_r.F )
C         ELSE
C           read indices from node-master specific file 
C           ( written in dist_ints_sigden --> interface_r.F )
C         ENDIF
C
C         SK-30-11-2007: since only LUCIAREL is parallel at present
C                        it is sufficient to allocate the 4-INDEX-array
C                        on the real (!) master
C                         
C         FIXME!!!       HAS TO BE CHANGED AS SOON AS THE COMPLETE
C                        KR-MCSCF IS RUNNING IN PARALLEL!
C
          IF(MYPROC .eq. MASTER)THEN
            WRITE(LUWRT,'(/a,i16/)')
     &      ' Allocating 4-index array of length',NTERMS*4
            call memmar_i8(KINDEX,NTERMS*4,'ADDL  ',1,'KINDEX')
            call IZERO8(WORK(KINDEX),NTERMS*4)
            call read_indices(WORK(KINDEX),NTERMS)
          END IF ! MASTER

          IPRDEN = 1
C
#if defined (VAR_MPI2)
C
csk       WRITE(LUWRT,*) ' MY index array '
csk       DO IJINDEX = 1, NTERMS 
csk       WRITE(LUWRT,'(A,I4)') 'index: ',WORK(KINDEX+IJINDEX-1)
csk       END DO
csk    WRITE(LUWRT,*) ' MY T (BUFF) array'
csk    IF( SPLIT_IJKL )THEN
csk      CALL WRTMATMN8(T_BUFF,1,MY_T_LEN,1,MY_T_LEN,LUWRT)
csk    ELSE
csk      CALL WRTMATMN8(WORK(KT_CC),1,NALLINT,1,NALLINT,LUWRT)
csk    END IF
          ISHARE_DENSM = 0
C
C         FIXME: control ISHARE_DENSM from the outside
C
C         case 1: only MASTER will get all density matrix elements
C             --> ISHARE_DENSM = 0
C         case 2: density matrix elements on all CPUs
C             --> ISHARE_DENSM = 1
C
          IF( IT_SHL .lt. 0 )THEN
C
            L_MAX_DENS_BLK = 0
            L_MAX_DENS_BLK = NCALC_MAX_BLK(IRC)
C
            IF(CIRUN .eq.'DENS1 ' .or. SPLIT_IJKL)THEN
              CALL MEMMAR_I8(KMAX_DENS_BLK,L_MAX_DENS_BLK,'ADDL  ',
     &                       1,'KDMBLK')
              CALL DZERO8(WORK(KMAX_DENS_BLK),L_MAX_DENS_BLK)
            END IF
C
            IF(CIRUN .eq. 'DENS1 ')THEN
              IF( SPLIT_IJKL ) THEN
                CALL COLLECT_DENS_MAT_SPLIT(WORK(KMAX_DENS_BLK),T_BUFF,
     &                                      IT_TTPL,IRC)
              ELSE 
              CALL COLLECT_DENS_MAT(WORK(KMAX_DENS_BLK),WORK(KT_CC),
     &                               IRC,ISHARE_DENSM)
              END IF
            ELSE
              IF( SPLIT_IJKL ) THEN
                CALL COLLECT_DENS_MAT_SPLIT(F2,T_BUFF,IT_TTPL,IRC)
              ELSE

                CALL COLLECT_DENS_MAT(F2,WORK(KT_CC),IRC,ISHARE_DENSM)
              END IF
            END IF ! 1- or 2-particle density matrix run

          ELSE ! MPI one-sided based T split
#if defined VAR_SHMEM
C
C           make sure each CPU gets here before starting T block 
C           collection
C
            call interface_mpi_BARRIER(global_communicator)
C
            IF( IT_SHL .eq. 0 )THEN
C
               MY_XT_BUFF_LEN = NCALC_MAX_BLK(IRC)
C
               CALL MPIXMEM_ALLOC(MY_XT_PTR,1,MY_XT_BUFF_LEN,
     &                            df_MPI_INFO_NULL,.FALSE.)
               CALL DZERO8(XT_BUFF,MY_XT_BUFF_LEN)
               CALL COLLECT_DENS_MAT_SM_LOC(T_BUFF,XT_BUFF,IT_TTPL,
     &                                      IT_TTOL,IRC,ISHARE_DENSM)
               CALL MPIXMEM_FREE(XT_BUFF)
C
            END IF
C
C           parallel MCSCF: provide WORK(KINDEX) to all CPU's
C
#endif /* VAR_SHMEM */
          END IF
C         ^ IT_SHL .lt. 0
C
#endif
C
!         IPRDEN = 10 ! debug
C
          IF (IPRDEN .ge. 10 ) THEN
            IF(MYPROC .eq. MASTER) THEN
              write(LUWRT,*) '======================================='
              write(LUWRT,'(A,I4)') ' Density matrix for root ',
     &                                JROOT_SCR
              write(LUWRT,*) '======================================='
C
              call wrtdens(WORK(KT_CC),WORK(KINDEX),WORK(KLIHIND),
     &                     WORK(KSIGNNHX),NTOOB2,WORK(KLABEXTP),
     &                     WORK(KLOP_REO),WORK(KSIOPREO),IRC)
            END IF
          END IF
C
#if defined (VAR_MPI2)
C
C         parallel MCSCF: use all CPU's to transfer the density 
C         matrices...
C
          IF( SPLIT_IJKL )THEN
C
C           collect information about density T blocks and their 
C           "contributors". Choose process with the lowest tag 
C           (MYPROC) as sender for the master (receiver). 
C           IT_TTPL --> gather --> sort --> write back on IT_TTPL
C           previous information on IT_TTPL is lost!
C
C           NOTE: --> parallel MCSCF: receiver may be all processes 
C           or a subgroup depending on the MCSCF integral distribution.  
C           That is, why we already use MPI_ALLGATHER!
C
            MY_DL_TMP_LEN = NSPOBEX_TP*NMPROC*IRC_SAVE
C
            allocate(IT_DLIST_TMP(MY_DL_TMP_LEN))
            IT_DLIST_TMP = 0
            call interface_mpi_ALLGATHER(IT_TTPL,NSPOBEX_TP*IRC_SAVE,
     &                         IT_DLIST_TMP,
     &                         NSPOBEX_TP*IRC_SAVE,global_communicator)
C
            CALL FIND_SENDER(IT_TTPL,IT_DLIST_TMP,NSPOBEX_TP)
C
csk         WRITE(LUWRT,*) ' AFTER SENDER matching IT_TTPL'
csk         CALL IWRTMAMN(IT_TTPL,1,NSPOBEX_TP*IRC,1,
csk  &                    NSPOBEX_TP*IRC,LUWRT)
C
C           release scratch memory
            deallocate(IT_DLIST_TMP)
          ELSE
            IF( MYPROC .ne. MASTER ) GOTO 919
          END IF
#endif
C
C         transfer density matrices to KRMC format
C         ----------------------------------------
C
C         local or global T distribution
          IF( IT_SHL .ge. 0 )THEN
#if defined (VAR_MPI2)
            IF(CIRUN .eq. 'DENS1 ')THEN
C
              call trnsfdens_no_sm(F1,T_BUFF,IT_TTPL,IT_TTOL,
     &                             WORK(KINDEX),
     &                             WORK(KLABEXTP),WORK(KLIHIND),
     &                             WORK(KSIGNNHX),WORK(KLOP_REO),
     &                             WORK(KSIOPREO),IRC,NTOOB2,IPRDEN)
C
            ELSE
C
              CALL DZERO(F2,LF2_ZERO)
C
C             FIXME: WRITE ME!!!
C
csk           call trnsfdens_sm(F1,F2,T_BUFF,IT_TTPL,IT_TTOL,
csk     &                       WORK(KINDEX),
csk     &                       WORK(KLABEXTP),WORK(KLIHIND),
csk     &                       WORK(KSIGNNHX),WORK(KLOP_REO),
csk     &                       WORK(KSIOPREO),IRC,NTOOB2,IPRDEN)
              call quit('GASCI_REL: missing a density transfer routine')
            END IF
#endif
C         T split / no T distribution
          ELSE
            IF(CIRUN .eq. 'DENS1 ')THEN
              IF( SPLIT_IJKL )THEN
#if defined (VAR_MPI2)
                call trnsfdens_no_split(F1,T_BUFF,WORK(KMAX_DENS_BLK),
     &                            IT_TTPL,WORK(KINDEX),WORK(KLABEXTP),
     &                            WORK(KLIHIND),
     &                            WORK(KSIGNNHX),WORK(KLOP_REO),
     &                            WORK(KSIOPREO),IRC,NTOOB2,IPRDEN)
#endif
              ELSE
                call trnsfdens_no(F1,WORK(KT_CC),WORK(KINDEX),
     &                            WORK(KLABEXTP),WORK(KLIHIND),
     &                            WORK(KSIGNNHX),WORK(KLOP_REO),
     &                            WORK(KSIOPREO),IRC,NTOOB2,IPRDEN)
              END IF
            ELSE ! ... two-particle (+ one-particle) density matrix
              CALL DZERO(F2,LF2_ZERO)
C
              IF( SPLIT_IJKL )THEN
#if defined (VAR_MPI2)
                call trnsfdens_split(F1,F2,T_BUFF,WORK(KMAX_DENS_BLK),
     &                        IT_TTPL,
     &                        WORK(KINDEX),WORK(KLABEXTP),WORK(KLIHIND),
     &                        WORK(KSIGNNHX),WORK(KLOP_REO),
     &                        WORK(KSIOPREO),IRC,NTOOB2,IPRDEN)
#endif
              ELSE
csk             IF( IDEN12 .eq. 1 ) THEN
csk               write(LUWRT,*) ' '
csk               write(LUWRT,*) '    this is my INDEX array        '
csk               write(LUWRT,*) ' '
csk               CALL IWRTMAMN8(WORK(KINDEX),1,NTERMS,1,NTERMS,LUWRT)
csk             END IF
                call trnsfdens(F1,F2,WORK(KT_CC),
     &                        WORK(KINDEX),WORK(KLABEXTP),WORK(KLIHIND),
     &                        WORK(KSIGNNHX),WORK(KLOP_REO),
     &                        WORK(KSIOPREO),IRC,NTOOB2,IDEN12,IPRDEN)
              END IF
            END IF ! 1-/2-particle density matrix switch
          END IF
C
#ifdef LUCI_DEBUG
          WRITE(LUWRT,*) ' MY updated densities'
          WRITE(LUWRT,*) '           F1        '
          CALL WRTMATMN(F1,1,N1ELINT,1,N1ELINT,LUWRT)
          IF(.NOT. CIRUN .eq. 'DENS1 ')THEN
            WRITE(LUWRT,*) '           F2        '
            CALL WRTMATMN8(F2,1,N2ELINT,1,N2ELINT,LUWRT)
          END IF
#endif
C
C
#if defined (VAR_MPI2)  
C
C        all slaves are waiting for the master who is transferring the
C        densities...
C
 919     CONTINUE
         call interface_mpi_BARRIER(global_communicator)
#endif        
C
      end if
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
      else
        WRITE(LUWRT,*) 'CIRUN has unallowed value :',CIRUN
        CALL QUIT(' *** ERROR in GASCI *** Unallowed CIRUN type.')
      end if ! CIRUN cases for MCSCF calls.
*+++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
#if defined (VAR_MPI2)
      IF( IT_SHL .ge. 0 ) THEN
#if defined VAR_SHMEM
C
C       de-allocate T array and destroy T memory window
        CALL MPIXMEM_FREE(IT_TTPL)
        CALL MPIXMEM_FREE(IT_TTOL)
        CALL REL_MEM_T_WIN(T_BUFF)
        IF ( CPIC_PASS_D ) MPIXMEM_FREE(T_BUFF_D)
#endif
      ELSE IF (SPLIT_IJKL ) THEN
        deallocate(IT_TTPL)
        deallocate(T_BUFF)
      END IF
#if !defined VAR_SHMEM
      IF( DIAGONAL_ALLOC ) deallocate(T_BUFF_D)
#endif
C
C     reset partitioning of CI vector to 'sequential' fashion
C     (*_MS2 arrays are in CSTATE )
      CALL Z_BLKFO_REL(IDC,NMS2VAL,
     &                 ISPC,ICSM,KPCLBT,KPCLEBT,KPCI1BT,KPCIBT,KPCBLTP,
     &                 NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,
     &                 IBAT_MS2,IDUMMY,.FALSE.,0,0,IDUMMY)
C
#endif
      if (CIRUN /= 'REFVEC') THEN
         CALL GETTIM(CPU1,WALL1)
         WALL1 = WALL1 - WALL0
         write(LUWRT,'(/5A/)') '  I have completed task ',CIRUN,
     &      ' in ',SECTID(WALL1),' (wall time)'
      END IF

C     release memory
      IDUM = 0
      call memmar(KDUM,  IDUM,  'FLUSM ',IDUM,'GASCI ')
C
      END
