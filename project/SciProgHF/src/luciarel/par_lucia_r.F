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

C
C
#if defined (VAR_MPI2)
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE BLOCK_DISTR_DRV_REL(NPARBLOCK,IBLOCKL_S,IBLOCKL_C,
     &                               RCCTOS,IT_TTPL,IGROUPLIST)
C
C     Driver routine for non-vanishing TTSS-block distribution among 
C     all processors.  
C
C     Written by  S. Knecht         - April 27 2007
C
C     T-block distribution for MPI based shared memory runs added.
C
C     S. Knecht                     - January 31 2008
C
C     new TTSS block distribution routine added - purpose:
C     fill up all node groups first with each block in a row 
C     --> try to get as even as possible afterwards.
C
C     S. Knecht                     - May 07 2008
C
C     optional reading from/writing to block file KRCI_BLOCKDIST.x 
C     added.
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE), IT_TTPL(*), IGROUPLIST(*)
      INTEGER   NPARBLOCK(*), IBLOCKL_S(*), IBLOCKL_C(*)
      INTEGER   RCCTOS(*)
#include "parluci.h"
#include "typesz_mpi2.h"
#include "wrkspc.inc"
#include "ipoist8.inc"
#include "krmc_shmem.h"
#include "mxdim_mpi2.h"
#include "ctcc.inc"
#include "cands.inc"
#include "mxpdim.inc"
#include "integrals_off.inc"
      POINTER (MY_TCTR_PTR,   IT_CTONE(IDUMMY_BUFF_SZ))
      POINTER (MY_TRED_PTR,   IT_CTALL(IDUMMY_BUFF_SZ))
      POINTER (MY_RCCTOS_PTR, IRCCTOS_2(IDUMMY_BUFF_SZ))
      POINTER (MY_IAPROC_PTR, I_APROC(IDUMMY_BUFF_SZ))
      POINTER (MY_BWEIGHT_PTR,I_BWEIGHT(IDUMMY_BUFF_SZ))
      POINTER (MY_BPROC_PTR,  I_BPROC(IDUMMY_BUFF_SZ))
      INTEGER(KIND=df_MPI_OFFSET_KIND)  I_BWEIGHT, I_BPROC
      INTEGER(KIND=df_MPI_ADDRESS_KIND) MY_TCTR_LEN
      INTEGER(KIND=df_MPI_ADDRESS_KIND) MY_TRED_LEN
      INTEGER(KIND=df_MPI_ADDRESS_KIND) MY_RCCTOS_LEN
      INTEGER(KIND=df_MPI_ADDRESS_KIND) MY_IAPROC_LEN
      CHARACTER BLOCKFLE*18
      LOGICAL EX_BFL, DIST_COMPLETE, SPLIT_IJKL_FILE, SHARED_M_FILE
      LOGICAL PRINT_BLOCK_DIST
      integer, allocatable :: GROUPCHECK(:)
C
      IMIN1 = -1
      IT_ACTIVE_R  = 0
      IT_ACTIVE_I  = 0
      IF( CPBLCK_FILE_LUCIX ) 
     &    WRITE(BLOCKFLE,'(A15,A3)') 'KRCI_BLOCKDIST.',SYMFLABEL
      PRINT_BLOCK_DIST = .TRUE.
C
C     set mark for local memory
      IDUM = 0
      CALL MEMMAR(KDUM,  IDUM,    'MARK  ',IDUM,'BLOCKD')
C
C     calculate # of 'active blocks'
      I_NZERO_LEN_SCR = 0
      I_NZERO_LEN_C   = 0
      I_NZERO_LEN_S   = 0
C
C     initialize RCCTOS - only real part considered
      CALL IZERO(RCCTOS,NUM_BLOCKS2)
      CALL I_DO_RCCTOS(IBLOCKL_S,RCCTOS,NUM_BLOCKS)
C
      DO IBLK = 1, NUM_BLOCKS
        IF( IBLOCKL_C(IBLK) .ne. 0 ) THEN 
          I_NZERO_LEN_SCR = I_NZERO_LEN_SCR + 1
        END IF
      END DO
C     transfer to common block /LUCIPARREL/
      I_NZERO_LEN_C = I_NZERO_LEN_SCR
      IF( CXPROPRUN ) THEN 
        I_NZERO_LEN_SCR = 0
        DO IBLK = 1, NUM_BLOCKS
          IF( IBLOCKL_S(IBLK) .ne. 0 ) THEN 
            I_NZERO_LEN_SCR = I_NZERO_LEN_SCR + 1
          END IF
        END DO
        CALL MEMMAR(KXRCCTOS,NUM_BLOCKS,'ADDL  ',1,'IXRCCT')
        CALL IZERO(WORK(KXRCCTOS),NUM_BLOCKS)
        CALL I_DO_RCCTOS(IBLOCKL_C,WORK(KXRCCTOS),NUM_BLOCKS)
      END IF
C     transfer to common block
      I_NZERO_LEN_S = I_NZERO_LEN_SCR
C
C     needed for T-split distribution
      ICOMM_SIZE_TMP = ICOMM_SIZE
      call interface_mpi_BCAST(ICOMM_SIZE_TMP,1,N_MASTER,
     &               MYNEW_COMM)
C
      CALL MEMMAR(KCWEIGHTF  ,I_NZERO_LEN_S        ,'ADDL  ',1,'ICWGHT')
      CALL MEMMAR(KACT_NUM   ,I_NZERO_LEN_S        ,'ADDL  ',1,'IACTNM')
      CALL MEMMAR(KAM_BLK_ACTC,NUM_BLOCKS          ,'ADDL  ',1,'IM_BLC')
      CALL MEMMAR(KAM_BLK_ACTS,NUM_BLOCKS          ,'ADDL  ',1,'IM_BLS')
      CALL MEMMAR(KICCTOS,I_NZERO_LEN_C*I_NZERO_LEN_S,
     &                                              'ADDL  ',1,'ICCTOS')
C
      CALL IZERO(WORK(KICCTOS)    ,I_NZERO_LEN_C*I_NZERO_LEN_S)
      CALL IZERO(WORK(KAM_BLK_ACTC),NUM_BLOCKS)
      CALL IZERO(WORK(KAM_BLK_ACTS),NUM_BLOCKS)
      CALL IZERO(WORK(KCWEIGHTF)  ,I_NZERO_LEN_S)
      CALL IZERO(WORK(KACT_NUM)  , I_NZERO_LEN_S)
      allocate(GROUPCHECK(NMPROC*ICOMM_SIZE_TMP))
      GROUPCHECK = 0
C
      CALL ISET_ARRAY_ACT_BLK(WORK(KAM_BLK_ACTS),IBLOCKL_S,NUM_BLOCKS)
      IF( CXPROPRUN ) THEN
        CALL ISET_ARRAY_ACT_BLK(WORK(KAM_BLK_ACTC),IBLOCKL_C,NUM_BLOCKS)
      ELSE
        CALL ICOPY(NUM_BLOCKS,WORK(KAM_BLK_ACTS),1,WORK(KAM_BLK_ACTC),1)
      END IF
C
C     check if we want to reuse (write) the block distribution file
C     KRCI_BLOCKDIST.x
      DIST_COMPLETE = .FALSE.
      IF( CPBLCK_FILE_LUCIX .and. (.NOT. CXPROPRUN ))THEN
        IF( MYPROC .eq. MASTER )THEN
          INQUIRE(FILE=BLOCKFLE,EXIST=EX_BFL)
          IF( EX_BFL ) 
     &       CALL GET_BLOCKFILE_DATA(BLOCKFLE,DIST_COMPLETE,NMPROC_FILE,
     &                               SPLIT_IJKL_FILE,SHARED_M_FILE,
     &                               WORK(KICCTOS),RCCTOS,NPARBLOCK,
     &                               WORK(KCWEIGHTF),IT_TTPL,
     &                               NSPOBEX_TP*IRC_SAVE)
        END IF
C       ^ end if MASTER
C    
C       slave update
        call interface_mpi_BCAST_l0(EX_BFL,1,MASTER,global_communicator)
        IF( EX_BFL )THEN
          CALL SYNC_NODES_BLOCKFILE_DATA(NMPROC_FILE,SPLIT_IJKL_FILE,
     &                                   SHARED_M_FILE,DIST_COMPLETE,
     &                                   WORK(KICCTOS),RCCTOS,NPARBLOCK,
     &                                   WORK(KCWEIGHTF),IT_TTPL,
     &                                   NSPOBEX_TP*IRC_SAVE)
          IF( .NOT. DIST_COMPLETE )THEN
C
C           incomplete information on file available, that is either 
C           - a block distribution where NMPROC(old) != NMPROC(act) 
C           - or an inconsistency between keywords as SPLIT_IJKL or 
C             SHARED_M (and the corresponding level)
            IF( NMPROC_FILE     .ne. NMPROC ) GOTO 100
            IF((SPLIT_IJKL_FILE .neqv. SPLIT_IJKL) .or.
     &        (SHARED_M_FILE .neqv. SHARED_M ) )GOTO 200
          END IF
        END IF
      END IF
C
C     complete information on file available, that is 
C     NMPROC(old) on file KRCI_BLOCKDIST.x matches NMPROC(act)
C     we have all required information
      IF(DIST_COMPLETE) GOTO 9999
C
 100  CHECK_TC = .FALSE.
C
C     1. check for sigma <==> C connections
      IF( CXPROPRUN )THEN
        CALL FIND_IMAT_SC_REL(WORK(KXRCCTOS),RCCTOS,WORK(KICCTOS),
     &                         WORK(KCWEIGHTF),WORK(KAM_BLK_ACTC),
     &                         WORK(KAM_BLK_ACTS),NUM_BLOCKS)
      ELSE
        CALL FIND_IMAT_SC_REL(RCCTOS,RCCTOS,WORK(KICCTOS),
     &                        WORK(KCWEIGHTF),WORK(KAM_BLK_ACTC),
     &                        WORK(KAM_BLK_ACTS),NUM_BLOCKS)
      END IF
C
      IF( SPLIT_IJKL )THEN
C
C       alternative block distribution
C       ------------------------------
C
C       gather on each node-master
        IF( MYNEW_ID .eq. N_MASTER )THEN
C
          call interface_mpi_ALLGATHER(IGROUPLIST,NMPROC,
     &                       GROUPCHECK,NMPROC,
     &                       ICOMM)
        END IF
C
        call interface_mpi_bcast_i1_work_f77(GROUPCHECK,
     &                 NMPROC*ICOMM_SIZE_TMP,
     &                 N_MASTER,MYNEW_COMM)
C
C       order group array
        CALL ORDER_ARRAY(GROUPCHECK,NMPROC*ICOMM_SIZE_TMP)
C
        MY_BWEIGHT_LEN = I_NZERO_LEN_S
        CALL MPIXMEM_ALLOC( MY_BWEIGHT_PTR, 3, MY_BWEIGHT_LEN,
     &                      df_MPI_INFO_NULL, .FALSE.)
        MY_BPROC_LEN   = NMPROC
        CALL MPIXMEM_ALLOC( MY_BPROC_PTR, 3, MY_BPROC_LEN,
     &                      df_MPI_INFO_NULL, .FALSE.)
C
        CALL IZERO8(I_BWEIGHT,MY_BWEIGHT_LEN)
        CALL IZERO8(I_BPROC  ,MY_BPROC_LEN)
        CALL BLOCK_DISTR_2(NUM_BLOCKS,IBLOCKL_S,NPARBLOCK,
     &                     WORK(KCWEIGHTF),
     &                     GROUPCHECK,I_BWEIGHT,I_BPROC,
     &                     WORK(KACT_NUM),ICOMM_SIZE_TMP,ISSM)
C       release memory
        CALL MPIXMEM_FREE(I_BWEIGHT)
        CALL MPIXMEM_FREE(I_BPROC)
C
        PRINT_BLOCK_DIST = .FALSE.
C
      ELSE ! split_ijkl
C
C       default block distribution - if required
        IF( COMPDISTL ) CALL BLOCK_DISTR_1(NUM_BLOCKS,IBLOCKL_S,
     &                                     NPARBLOCK,WORK(KCWEIGHTF),
     &                                     ISSM)
        PRINT_BLOCK_DIST = .FALSE.
C
      END IF
C
C     find all c-blocks connecting to a given sigma-block,
C     store information in RCCTOS
      CALL IZERO(RCCTOS,NUM_BLOCKS2)
      CALL FIND_CCTOS_REL(RCCTOS,NPARBLOCK,WORK(KICCTOS),IBLOCKL_S,
     &                    IBLOCKL_C,IRC_SAVE,NUM_BLOCKS)
C
 200  IF( SPLIT_IJKL )THEN
C
C       find required T-blocks
        MY_TCTR_LEN   = NSPOBEX_TP * IRC_SAVE
        CALL MPIXMEM_ALLOC( MY_TCTR_PTR ,   2, MY_TCTR_LEN, 
     &                      df_MPI_INFO_NULL, .FALSE.)
C
        CALL IZERO(IT_CTONE, NSPOBEX_TP * IRC_SAVE )
        CALL IZERO(IT_TTPL , NSPOBEX_TP * IRC_SAVE )
C       --> use connection vector RCCTOS and NPARBLOCK to produce
C           output on IT_CTONE
C
C       2. check for T-connections
        CALL FIND_IMAT_SC_REL(RCCTOS,RCCTOS,DUMMY,IT_CTONE,NPARBLOCK,
     &                          NPARBLOCK,-1)
C
C       save connections
        LEN_COPY = NSPOBEX_TP*IRC_SAVE
        CALL ICOPY(LEN_COPY,IT_CTONE,1,IT_TTPL,1)
C
csk     WRITE(LUWRT,*) ' IT_CTONE'
csk     CALL IWRTMAMN(IT_CTONE,1,NSPOBEX_TP*IRC_SAVE,1,
csk  &                NSPOBEX_TP*IRC_SAVE,LUWRT)
csk     WRITE(LUWRT,*) ' IT_TTPL'
csk     CALL IWRTMAMN(IT_TTPL,1,NSPOBEX_TP*IRC_SAVE,1,
csk  &                NSPOBEX_TP*IRC_SAVE,LUWRT)
C
C       release scratch memory
        CALL MPIXMEM_FREE(IT_CTONE)
      END IF
C
      IF( SHARED_M ) THEN
        IF( IT_SHL .ge. 0 )THEN
C
C         allocate memory via MPI calls for scratch matrices needed
C
          MY_TCTR_LEN   = NSPOBEX_TP * IRC_SAVE
          MY_TRED_LEN   = NSPOBEX_TP * NEWCOMM_PROC_SM * IRC_SAVE
          MY_RCCTOS_LEN = NUM_BLOCKS
          MY_IAPROC_LEN = NEWCOMM_PROC_SM
C
          CALL MPIXMEM_ALLOC( MY_TCTR_PTR ,   2, MY_TCTR_LEN, 
     &                        df_MPI_INFO_NULL, .FALSE.)
          CALL MPIXMEM_ALLOC( MY_TRED_PTR ,   2, MY_TRED_LEN, 
     &                        df_MPI_INFO_NULL, .FALSE.)
          CALL MPIXMEM_ALLOC( MY_RCCTOS_PTR , 2, MY_RCCTOS_LEN, 
     &                        df_MPI_INFO_NULL, .FALSE.)
          CALL MPIXMEM_ALLOC( MY_IAPROC_PTR , 2, MY_IAPROC_LEN, 
     &                        df_MPI_INFO_NULL, .FALSE.)
C
C         in case of MPI based "shared memory" runs, it is important
C         to know how to distribute all integrals among the CPUs
C         |
C         --> use connection vector RCCTOS and NPARBLOCK to produce
C             output on IT_CTONE / IT_CTALL <==> input for T_BLOCK_DISTR
C
          CALL IZERO(IT_CTONE, NSPOBEX_TP * IRC_SAVE )
          CALL IZERO(IT_CTALL, NSPOBEX_TP * IRC_SAVE * NEWCOMM_PROC_SM)
          CALL IZERO(IRCCTOS_2,NUM_BLOCKS)
          CALL ISETVC(I_APROC,IMIN1,NEWCOMM_PROC_SM)
C         initialize for a second time
          CALL I_DO_RCCTOS(IBLOCKL_C,IRCCTOS_2,NUM_BLOCKS)
          CHECK_TC = .TRUE.
C
C         3. check for T-connections
          CALL FIND_IMAT_SC_REL(IRCCTOS_2,IRCCTOS_2,DUMMY,IT_CTONE,
     &                          DUMMY,DUMMY,NUM_BLOCKS)
C
C
          call interface_mpi_ALLGATHER(IT_CTONE,NSPOBEX_TP*IRC_SAVE,
     &                       IT_CTALL,
     &                       NSPOBEX_TP*IRC_SAVE,
     &                       MYNEW_COMM_SM)
C
C         find active T-blocks...
          CALL IZERO(IT_CTONE,NSPOBEX_TP*IRC_SAVE)
C
          IT_ACTIVE_R = 0
          IT_ACTIVE_R = IPICK_T_ACTIVE(IT_CTALL,IT_CTONE,
     &                                 NEWCOMM_PROC_SM,NSPOBEX_TP,1)
          IT_ACTIVE_I = 0
          IF( IRC_SAVE .eq. 2 ) THEN
            IT_ACTIVE_I = IPICK_T_ACTIVE(IT_CTALL,IT_CTONE,
     &                                   NEWCOMM_PROC_SM,NSPOBEX_TP,2)
          END IF
C
csk       WRITE(LUWRT,*) ' IT_ACTIVE_R, IT_ACTIVE_I, NSPOBEX_TP',
csk  &                     IT_ACTIVE_R, IT_ACTIVE_I, NSPOBEX_TP
C
          CALL T_BLOCK_DISTR(IT_TTPL,IT_CTALL,IT_CTONE,I_APROC,
     &                       NSPOBEX_TP,IT_ACTIVE_R,IT_ACTIVE_I)
csk       WRITE(LUWRT,*) ' showing array IT_TTPL'
csk       CALL IWRTMAMN(IT_TTPL,1,NSPOBEX_TP,1,NSPOBEX_TP,LUWRT)
C
C         release scratch memory
          CALL MPIXMEM_FREE(IT_CTONE)
          CALL MPIXMEM_FREE(IT_CTALL)
          CALL MPIXMEM_FREE(IRCCTOS_2)
          CALL MPIXMEM_FREE(I_APROC)
        END IF
      END IF
C
C     write block distribution file if asked for
      IF( CPBLCK_FILE_LUCIX .and. (.NOT. CXPROPRUN ))THEN
        IF( MYPROC .eq. MASTER )THEN
          CALL PUT_BLOCKFILE_DATA(BLOCKFLE,WORK(KICCTOS),RCCTOS,
     &                            NPARBLOCK,WORK(KCWEIGHTF),
     &                            IT_TTPL,NSPOBEX_TP*IRC_SAVE)
        END IF
      END IF
C
C     print block distribution (if we restart from a previous run)
9999  IF( PRINT_BLOCK_DIST .and. MYPROC .eq. MASTER ) 
     &    CALL PRINT_BLOCK_DISTR(NUM_BLOCKS,IBLOCKL_S,NPARBLOCK,
     &                           WORK(KCWEIGHTF),ISSM,BLOCKFLE)
C
C     determine whether this process is "active" or not
      IAM_NOT_INV = 0
      IAM_NOT_INV = ISACTIVE_CPU(NPARBLOCK,NUM_BLOCKS,MYPROC)
      IF( IAM_NOT_INV .eq. 0) CALL ISETVC(NPARBLOCK,-2,NUM_BLOCKS)
C
C     flush local memory
      deallocate(GROUPCHECK)
      IDUM = 0
      CALL MEMMAR(KDUM ,IDUM,'FLUSM ',2,'BLOCKD')
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE BLOCK_DISTR_1(NDIM,IBLOCKL,NPARBLOCK,ICWEIGHTF,
     &                         JSYM_DISTBLK)
C
C     Block distribution routine for non-vanishing TTSS-block - initial 
C     version.
C
C     Written by  S. Knecht         - April 30 2007
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "parluci.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
C
C     -----
C     INPUT
C     -----
C
C     NDIM:            total number of TTSS-blocks
C     IBLOCKL(NDIM):   total length of each block
C     ICWEIGHTF(NDIM): total 'weight factors' for each block
C     JSYM_DISTBLK   : symmetry irrep for which blocks are to be
C                      distributed
C
C     ------
C     OUTPUT
C     ------
C
C     NPARBLOCK(NDIM): list of blocks containing the corresponding assigned
C     node
C
      INTEGER ICOUNTABLK, NPTEST
      DIMENSION ICWEIGHTF(*)
      DIMENSION NPARBLOCK(NDIM),IBLOCKL(NDIM)
C     -------
C     SCRATCH
C     -------
C     MXSIZE: actual maximum size for one of all blocks
C     MXNUMB: number of the current largest block
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
      INTEGER IMINUS2
      INTEGER*8 ITLTND,MXSIZE,MXSZTMP
      INTEGER*8 ITEMPL, ITEMPW1, ITEMPW2
      INTEGER*8 ILEN, ITOTBLCKL, IPRODLEN
      INTEGER*8 IIWEIGHTBLK, IILENGTH, IILENGTH_MAX
C
      NPTEST = MAX(NPTEST_VAR,0)
      IMINUS2 = 0
      IIWEIGHTBLK = 0
      IILENGTH = 0
      ITLTND = 0
      ICOUNTABLK = 0
      MXSIZE = 0
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 0
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      ITOTBLCKL = 0
      ILEN = 0
      IONE = 1
      I_AM_ACTIVE = 0
      IPRODLEN = 0
      IILENGTH_MAX = 0
*     print some information
      NPTEST = 1
C     end of initialization
      DO ICBL = 1, NDIM
C
        IILENGTH = IBLOCKL( ICBL )
C
        IF( IILENGTH .gt. 0 )THEN
C
          I_AM_ACTIVE = I_AM_ACTIVE + IONE
C
          IIWEIGHTBLK = ICWEIGHTF( I_AM_ACTIVE )
C
          IPRODLEN = IILENGTH * IIWEIGHTBLK
C
          ILEN = IPRODLEN
          IILENGTH_MAX = MAX(IILENGTH_MAX,IILENGTH)
C
          IF( ILEN .lt. 0 ) THEN
            WRITE(6,*)'Attention, minus block detected',ICBL
            WRITE(6,*)'ILEN = ',ILEN
            WRITE(6,*)'IBLOCKL(ICBL) = ',IBLOCKL(ICBL)
            WRITE(6,*)'ICWEIGHTF(ICBL) = ',ICWEIGHTF( I_AM_ACTIVE )
            call quit('*** error in block_distr_1: block type 
     & (active/inactive) not detected.***')
          END IF
          IF( ILEN .ne. 0 )THEN
             ICOUNTABLK = ICOUNTABLK + 1
             ITOTBLCKL = ITOTBLCKL + ILEN
             NPARBLOCK( ICBL ) = -1
          END IF
          IF(ILEN.GE.MXSIZE) MXSIZE = ILEN
        END IF
*       ^ IILENGTH > 0
      END DO
*
CSK      WRITE(LUWRT,*)'MXSIZE of MYPROC',MXSIZE,MYPROC
      ITEMPW1 = MXSIZE
*
      IF(NPTEST.GE.0) THEN
        WRITE(LUWRT,'(/7X,A)')
     &  '==================================================='
        WRITE(LUWRT,'(7X, A,I3)')
     &  ' parallel distribution setup for symmetry irrep ',JSYM_DISTBLK
        WRITE(LUWRT,'(7X,A/)')
     &  '==================================================='
         WRITE(LUWRT,'(A,I19)')
     &'  total number of processes to distribute on :',NMPROC
         WRITE(LUWRT,'(A,I19)')
     &'  total number of blocks                     :',NDIM
         WRITE(LUWRT,'(A,I19)') 
     &'  total number of active blocks              :',ICOUNTABLK
         WRITE(LUWRT,'(A,I19)')
     &'  size of largest TTSS block                 :',IILENGTH_MAX 
         WRITE(LUWRT,'(A,I19)')
     &'  overall weighted active block length       :',ITOTBLCKL
         WRITE(LUWRT,'(A,I19)')
     &'  Maximum weighted block size                :',MXSIZE
      END IF
*
      IMINNP = NMPROC
      NTEMPP = IMINNP
C     store on common block
      ICOUNTABLK_C = ICOUNTABLK
C
      IF( ICOUNTABLK .lt. NMPROC ) THEN
        IMINNP = MIN( ICOUNTABLK,IMINNP )
        if(IMINNP.lt.NMPROC)then
        write(luwrt,'(/a,/a,/a,i6,/a)') 
     &  '  **************** info from block_distr_1 *****************',
     &  '  number of active blocks lower than the number of processes.',
     &  '  please decrease the number of processes to:',IMINNP,
     &  '  or can you afford to let the other processes idle? ;)'
        write(luwrt,'(a,/a,/a)') 
     &  '  alternatively change the GAS specification (more GA spaces)',
     &  '  to increase the number of TTSS blocks.',
     &  '  **********************************************************'
        end if
        NTEMPP = IMINNP
        IF( MYPROC .ge. ICOUNTABLK ) THEN 
          IAM_NOT_INV = 0
          IMINUS2 = - 2
          CALL ISETVC(NPARBLOCK,IMINUS2,NDIM)
          GOTO 8888
        END IF
      END IF
      IAM_NOT_INV  = 1
C
C     starting the treausure quest for the ?optimal? block distribution
C
100   CONTINUE
C
      IRUN = IRUN + 1
      IF(IRUN.LE.ICOUNTABLK) THEN
C
        MXSZTMP = 0
C
        I_AM_ACTIVE = 0
C
        DO 3000 II = 1, NDIM
C
          IILENGTH = IBLOCKL(II)
C
          IF( IILENGTH .gt. 0 ) THEN
C
            I_AM_ACTIVE = I_AM_ACTIVE + IONE
C
            IIWEIGHTBLK = ICWEIGHTF( I_AM_ACTIVE )
C
            ITEMPL = IILENGTH * IIWEIGHTBLK
C
            ITEMPD = NPARBLOCK(II)
            ITEMPN = II

CSK          IF(MYPROC.EQ.0)
CSK     & WRITE(LUWRT,*) 'ITEMPL,ITEMPD,ITEMPN:',ITEMPL,ITEMPD,ITEMPN
C
CSK          IF(MYPROC.EQ.0) WRITE(LUWRT,*) 'MXSIZE:',MXSIZE
CSK          IF(MYPROC.EQ.0) WRITE(LUWRT,*) 'MXSZTMP:',MXSZTMP
            IF(ITEMPL.GT.0)THEN
              IF(ITEMPD.EQ.-1) THEN
                IF(ITEMPL.LE.MXSIZE)THEN
                  IF(ITEMPL.GT.MXSZTMP) THEN
                     MXSZTMP = ITEMPL
                     MXNUMB = ITEMPN
CSK                    WRITE(LUWRT,*) 'MXSZTMP (2):',MXSZTMP
                  END IF
                END IF
              END IF
            END IF
          END IF
C         ^ IILENGTH > 0
3000    CONTINUE
C
C
        DO 4000 IPR = 1, NTEMPP
C
          I_AM_ACTIVE = 0
C
          ITEMPW2 = 0
          DO IBLK = 1, NDIM
            IILENGTH = IBLOCKL(IBLK)
            IF( IILENGTH .ne. 0 ) THEN
              I_AM_ACTIVE = I_AM_ACTIVE + IONE
              IF( NPARBLOCK(IBLK) .eq. IPR-1 ) THEN
                IIWEIGHTBLK = ICWEIGHTF( I_AM_ACTIVE )
                ITEMPW2 = ITEMPW2 + ( IILENGTH * IIWEIGHTBLK )
              END IF
            END IF
          END DO

C         attention: IAMTPROC = 1 --> MASTER = 0 !
          IAMTPROC = IPR
C
          IF(IPR.EQ.1) ITEMPW1 = ITEMPW2
C
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC
          END IF
C
4000    CONTINUE
C      now we should have found a proc and a block --> put both together !
C
C      calculation of block MXNUMB by proc IAMGPROC - 1
       IF( NPTEST .gt. 10 ) THEN
         WRITE(LUWRT,*) 'calculation of block MXNUMB by proc',
     &                   MXNUMB,IAMGPROC-1
       ENDIF
C
       IAMTPROC = IAMGPROC-1
       MXSIZE = MXSZTMP
       NPARBLOCK(MXNUMB)     = IAMTPROC
C
      ELSE
        GOTO 101
      END IF
C     /\ IRUN !!!
C
      GOTO 100
C
101   CONTINUE
C
      I_AM_ACTIVE = 0
C
      IF(NPTEST.GE.0) THEN
        DO II = 1, NDIM
         IF( NPARBLOCK(II) .eq. -1 ) THEN
           WRITE(LUWRT,*) 'AAAHHHH, block',II,' is not distributed'
           WRITE(LUWRT,*) 'Since I do not know how to proceed, 
     &                     I will stop!'
             call quit('*** error in block_distr_1: active block
     & not distributed among the processes.***')
          END IF
        END DO
        WRITE(LUWRT,'(/20X,A)')'================================'
        WRITE(LUWRT,'(20X, A)')' Summation of even distribution '
        WRITE(LUWRT,'(20X,A/)')'================================'
        DO ISTI = 1, NTEMPP
          ITLTND = 0
          ITLBND = 0
          I_AM_ACTIVE = 0
          DO IBLK = 1, NDIM
            IILENGTH = IBLOCKL(IBLK)
            IF( IILENGTH .gt. 0 )THEN
              I_AM_ACTIVE = I_AM_ACTIVE + IONE
              IF( NPARBLOCK(IBLK) .eq. ISTI -1 ) THEN
C
                 IIWEIGHTBLK = ICWEIGHTF( I_AM_ACTIVE )
C
                 ITLTND = ITLTND + ( IILENGTH * IIWEIGHTBLK )
                 ITLBND = ITLBND + 1
              END IF
            END IF
          END DO
        WRITE(LUWRT,'(2X,A,1X,I4,1X,A,1X,I6,1X,A,1X,I18)')
     &'CPU',ISTI-1,' computes',ITLBND,'blocks with a total weight of',
     &   ITLTND
        END DO
      END IF
8888  CONTINUE
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE BLOCK_DISTR_2(NDIM,IBLOCKL,NPARBLOCK,ICWEIGHTF,
     &                         IALLGRPLST,I_BWEIGHT,I_BPROC,
     &                         IACT_NUM,ICOMM_SZ_TMP,JSYM_DISTBLK)
C
C     Block distribution routine for non-vanishing TTSS-block.
C
C     JSYM_DISTBLK   : symmetry irrep for which blocks are to be
C                      distributed
C
C     purpose: 
C
C     'total weight'/#procs = 'averaged weight'
C
C     start with process group 1 and assign as many blocks as possible
C     then continue with next group...
C
C     Finally, assign all blocks that are left
C
C     Written by  S. Knecht         - April 30 2007
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "parluci.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
C
C     -----
C     INPUT
C     -----
C
C     NDIM:                     total number of TTSS-blocks
C     IBLOCKL(NDIM):            total length of each block
C     ICWEIGHTF(active blocks): total 'weight factors' for each block
C     IACT_NUM(active blocks):  corresponding block number in total list
C     IALLGRPLST(NMPROC*ICOMM_SIZE_TMP) : 
C                               list of all node groups and their CPU tags
C
C     ------
C     OUTPUT
C     ------
C
C     NPARBLOCK(NDIM): list of blocks containing the corresponding assigned
C     node
C
      INTEGER ICOUNTABLK, NPTEST
      DIMENSION ICWEIGHTF(*), IACT_NUM(*)
      DIMENSION NPARBLOCK(NDIM), IBLOCKL(NDIM), IALLGRPLST(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) I_BWEIGHT(*), I_BPROC(*)
C     -------
C     SCRATCH
C     -------
      INTEGER IMINUS2
      INTEGER*8 ITLTND, MXSIZE, ICURR_PLOAD, ICURR_BLOAD
      INTEGER*8 ILEN, ITOTBLCKL, IPRODLEN, IMIN_CNT, IMIN_CNT_TMP
      INTEGER*8 IIWEIGHTBLK, IILENGTH, IILENGTH_MAX
C
      NPTEST       = MAX(NPTEST_VAR,0)
      IMINUS2      = 0
      IIWEIGHTBLK  = 0
      IILENGTH     = 0
      ITLTND       = 0
      ICOUNTABLK   = 0
      MXSIZE       = 0
      IRUN         = 0
      ITOTBLCKL    = 0
      ILEN         = 0
      IONE         = 1
      I_AM_ACTIVE  = 0
      IPRODLEN     = 0
      IILENGTH_MAX = 0
*     print some information
      NPTEST = 1
C     end of initialization
      DO ICBL = 1, NDIM
C
        IILENGTH = IBLOCKL(ICBL)
C
        IF( IILENGTH .gt. 0 )THEN
C
          I_AM_ACTIVE              = I_AM_ACTIVE + IONE
C
          IIWEIGHTBLK              = ICWEIGHTF(I_AM_ACTIVE)
          ILEN                     = IILENGTH * IIWEIGHTBLK
          I_BWEIGHT( I_AM_ACTIVE ) = ILEN
csk       WRITE(LUWRT,*) ' ILEN, I_AM_ACTIVE, I_BWEIGHT(...)',
csk  &                     ILEN, I_AM_ACTIVE, I_BWEIGHT(I_AM_ACTIVE)
          IILENGTH_MAX             = MAX(IILENGTH_MAX,IILENGTH)
          MXSIZE                   = MAX(MXSIZE,ILEN)
          NPARBLOCK( ICBL )        = -1
          IACT_NUM( I_AM_ACTIVE )  = ICBL
          ITOTBLCKL                = ITOTBLCKL + ILEN
        END IF
*       ^ IILENGTH > 0
      END DO
C
C     debug print      
C
#ifdef MOD_DEBUG
      WRITE(LUWRT,*) '  I_AM_ACTIVE is ',I_AM_ACTIVE
      WRITE(LUWRT,*) '  *********** BLOCK WEIGHTs ************'
      CALL IWRTMAMN8(I_BWEIGHT,1,I_AM_ACTIVE,1,I_AM_ACTIVE,LUWRT)
      WRITE(LUWRT,*) ' ICOMM_SZ_TMP, ICOMM_SZ_TMP*NMPROC', 
     &                 ICOMM_SZ_TMP, ICOMM_SZ_TMP*NMPROC
#endif
C
      ICOUNTABLK = I_AM_ACTIVE
*
      IMINNP = NMPROC
      NTEMPP = IMINNP
C     store on common block
      ICOUNTABLK_C = ICOUNTABLK
C
      IAM_NOT_INV = 0
      IF( ICOUNTABLK .lt. NMPROC ) THEN
        IMINNP = MIN( ICOUNTABLK,IMINNP )
        if(IMINNP.lt.NMPROC)then
        write(luwrt,'(/a,/a,/a,i6,/a)') 
     &  '  **************** info from block_distr_2 *****************',
     &  '  number of active blocks lower than the number of processes.',
     &  '  please decrease the number of processes to:',IMINNP,
     &  '  or can you afford to let the other processes idle? ;)'
        write(luwrt,'(a,/a,/a)') 
     &  '  alternatively change the GAS specification (more GA spaces)',
     &  '  to increase the number of TTSS blocks.',
     &  '  **********************************************************'
        end if
        NTEMPP = IMINNP
        DO JJ = 1, NTEMPP
          IF( MYPROC .eq. IALLGRPLST(JJ) ) IAM_NOT_INV = 1 
        END DO
        IF( IAM_NOT_INV .eq. 0 ) THEN 
          IMINUS2 = - 2
          CALL ISETVC(NPARBLOCK,IMINUS2,NDIM)
          GOTO 8888
        END IF
      END IF
      IAM_NOT_INV  = 1
C
      IPRODLEN = ITOTBLCKL/NTEMPP
C
      IF(NPTEST.GE.0) THEN
        WRITE(LUWRT,'(/7X,A)')
     &  '==================================================='
        WRITE(LUWRT,'(7X, A,I3)')
     &  ' parallel distribution setup for symmetry irrep ',JSYM_DISTBLK
        WRITE(LUWRT,'(7X,A/)')
     &  '==================================================='
         WRITE(LUWRT,'(A,I19)')
     &'  total number of processes to distribute on :',NMPROC
         WRITE(LUWRT,'(A,I19)')
     &'  total number of blocks                     :',NDIM
         WRITE(LUWRT,'(A,I19)') 
     &'  total number of active TTSS blocks         :',ICOUNTABLK
         WRITE(LUWRT,'(A,I19)')
     &'  size of largest TTSS block                 :',IILENGTH_MAX 
         WRITE(LUWRT,'(A,I19)')
     &'  overall weighted active block length       :',ITOTBLCKL
         WRITE(LUWRT,'(A,I19)')
     &'  Maximum weighted block size                :',MXSIZE
         WRITE(LUWRT,'(A,I19)')
     &'  aiming at an block weight per process of   :',IPRODLEN
      END IF
C
C     starting the treausure quest for the ?optimal? block distribution
C
C     start with MASTER (and its group of processes)
C
      ICURR_PR    = 0
      ICURR_PLOAD = 0
      ICURR_BLOAD = 0
      ICURR_PTAG  = 0
C
      DO 300 II = 1, ICOUNTABLK
         ICURR_BLOAD = I_BWEIGHT(II)
csk      WRITE(LUWRT,*) ' ICURR_PTAG,ICURR_PLOAD,ICURR_PR,ICURR_BLOAD',
csk  &                    ICURR_PTAG,ICURR_PLOAD,ICURR_PR,ICURR_BLOAD
         IF( ICURR_PTAG .lt. NTEMPP )THEN
           IF( ICURR_PLOAD + ICURR_BLOAD .le. IPRODLEN)THEN
               ICURR_PLOAD       = ICURR_PLOAD + ICURR_BLOAD
               NPARBLOCK(IACT_NUM(II))   = ICURR_PR
csk            WRITE(LUWRT,*) ' assigning block II to ICURR_PR',
csk  &                          II,ICURR_PR
           ELSE
C              write processor load to array
               I_BPROC(ICURR_PR+1) = ICURR_PLOAD
C              next processor
               ICURR_PTAG  = ICURR_PTAG + 1
               IF( ICURR_PTAG .lt. NTEMPP ) THEN
csk              WRITE(LUWRT,*) ' 3. ICURR_PTAG+1',ICURR_PTAG+1
                 ICURR_PR                  = IALLGRPLST(ICURR_PTAG+1)
                 ICURR_PLOAD               = ICURR_BLOAD
                 NPARBLOCK(IACT_NUM(II))   = ICURR_PR
csk              WRITE(LUWRT,*) ' 3. assigning block II to ICURR_PR',
csk  &                          II,ICURR_PR
               ELSE
C                end of all temporary available CPU's;
C                assign current block to CPU with lowest count
                 IMIN_CNT_TMP = 0
                 IMIN_CNT     = ICURR_PLOAD
                 DO JJ = 1, NTEMPP
                   IMIN_CNT_TMP = I_BPROC(JJ)
                   IF( IMIN_CNT_TMP .lt. IMIN_CNT ) ICURR_PR = JJ -1
                 END DO
C
                 I_BPROC(ICURR_PR+1) = I_BPROC(ICURR_PR+1)+ICURR_BLOAD
                 NPARBLOCK(IACT_NUM(II))   = ICURR_PR
                 ICURR_PLOAD = I_BPROC(ICURR_PR+1)
csk              WRITE(LUWRT,*) ' 2. assigning block II to ICURR_PR',
csk  &                          II,ICURR_PR
               END IF
           END IF
         ELSE
C
C          find CPU with lowest count - start with MASTER
C
           ICURR_PR     = 0
           IMIN_CNT_TMP = 0
           IMIN_CNT     = I_BPROC(ICURR_PR+1)
           DO JJ = 1, NTEMPP
              IMIN_CNT_TMP = I_BPROC(JJ)
              IF( IMIN_CNT_TMP .lt. IMIN_CNT ) ICURR_PR = JJ -1
           END DO
C
           I_BPROC(ICURR_PR+1)       = I_BPROC(ICURR_PR+1)+ICURR_BLOAD
           NPARBLOCK(IACT_NUM(II))   = ICURR_PR
csk        WRITE(LUWRT,*) ' 4. assigning block II to ICURR_PR',
csk  &                      II,ICURR_PR
C
         END IF
 300     CONTINUE
C
      IF(NPTEST.GE.0) THEN
        DO II = 1, NDIM
         IF( NPARBLOCK(II) .eq. -1 ) THEN
           WRITE(LUWRT,*) 'AAAHHHH, block',II,' is not distributed'
           WRITE(LUWRT,*) 'Since I do not know how to proceed, 
     &                     I will stop!'
             call quit('*** error in block_distr_2: active block
     & not distributed among the processes.***')
          END IF
        END DO
        WRITE(LUWRT,'(/20X,A)')'================================'
        WRITE(LUWRT,'(20X, A)')' Summation of even distribution '
        WRITE(LUWRT,'(20X,A/)')'================================'
        DO ISTI = 1, NTEMPP
          ITLTND = 0
          ITLBND = 0
          I_AM_ACTIVE = 0
          DO IBLK = 1, NDIM
            IILENGTH = IBLOCKL(IBLK)
            IF( IILENGTH .gt. 0 )THEN
              I_AM_ACTIVE = I_AM_ACTIVE + IONE
              IF( NPARBLOCK(IBLK) .eq. ISTI -1 ) THEN
C
                 IIWEIGHTBLK = ICWEIGHTF( I_AM_ACTIVE )
C
                 ITLTND = ITLTND + ( IILENGTH * IIWEIGHTBLK )
                 ITLBND = ITLBND + 1
              END IF
            END IF
          END DO
        WRITE(LUWRT,'(2X,A,1X,I4,1X,A,1X,I6,1X,A,1X,I18)')
     &'CPU',ISTI-1,' computes',ITLBND,'blocks with a total weight of',
     &   ITLTND
        END DO
      END IF
8888  CONTINUE
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PRINT_BLOCK_DISTR(NDIM,IBLOCKL,NPARBLOCK,ICWEIGHTF,
     &                             JSYM_DISTBLK,BLKFILE)
C
C     print block distribution as read from file
C     KRCI_BLOCKDIST."JSYM_DISTBLK"
C
C     Written by  S. Knecht         - Dec 2008
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "parluci.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE), IBLOCKL(*), NPARBLOCK(*)
      INTEGER   ICWEIGHTF(*)
      CHARACTER BLKFILE*16
      INTEGER*8 ITLTND, MXSIZE, ILEN, ITOTBLCKL, IPRODLEN
      INTEGER*8 IIWEIGHTBLK, IILENGTH, IILENGTH_MAX
C
      I_AM_ACTIVE  = 0
      IILENGTH_MAX = 0
      IPRODLEN     = 0
      ICOUNTABLK   = 0
      ITOTBLCKL    = 0
      ILEN         = 0
      MXSIZE       = 0
C
      WRITE(LUWRT,'(/7X,A)')
     &'=================================================='
      WRITE(LUWRT,'(7X,A,A16)')
     &'Block distribution read from file ',BLKFILE
      WRITE(LUWRT,'(7X,A)')
     &'=================================================='
C
      DO ICBL = 1, NDIM
        IILENGTH = IBLOCKL(ICBL)
        IF( IILENGTH .gt. 0 )THEN
C
          I_AM_ACTIVE  = I_AM_ACTIVE + 1
          IIWEIGHTBLK  = ICWEIGHTF(I_AM_ACTIVE)
          IPRODLEN     = IILENGTH * IIWEIGHTBLK
          ILEN         = IPRODLEN
          IILENGTH_MAX = MAX(IILENGTH_MAX,IILENGTH)
C
          IF( ILEN .lt. 0 ) THEN
            WRITE(LUWRT,'(/A)')' *** ERROR in PRINT_BLOCK_DISTR *** '
            WRITE(LUWRT,'(A,I10)')' minus block detected. # =',ICBL
            WRITE(LUWRT,*)' ILEN = ',ILEN
            WRITE(LUWRT,*)' IBLOCKL(ICBL) = ',IBLOCKL(ICBL)
            WRITE(LUWRT,*)' ICWEIGHTF(ICBL) = ',ICWEIGHTF( I_AM_ACTIVE )
            CALL QUIT(' *** ERROR in PRINT_BLOCK_DISTR *** Minus 
     & block detected!')
          END IF
          IF( ILEN .ne. 0 )THEN
             ICOUNTABLK = ICOUNTABLK + 1
             ITOTBLCKL = ITOTBLCKL + ILEN
          END IF
          IF(ILEN.GE.MXSIZE) MXSIZE = ILEN
        END IF
C       ^ IILENGTH > 0
      END DO
C
      WRITE(LUWRT,'(/7X,A)')
     &'==================================================='
      WRITE(LUWRT,'(7X, A,I3)')
     &' parallel distribution setup for symmetry irrep ',JSYM_DISTBLK
      WRITE(LUWRT,'(7X,A/)')
     &'==================================================='
       WRITE(LUWRT,'(A,I19)')
     &'  total number of processes to distribute on :',NMPROC
       WRITE(LUWRT,'(A,I19)')
     &'  total number of blocks                     :',NDIM
       WRITE(LUWRT,'(A,I19)') 
     &'  total number of active blocks              :',ICOUNTABLK
       WRITE(LUWRT,'(A,I19)')
     &'  size of largest TTSS block                 :',IILENGTH_MAX 
       WRITE(LUWRT,'(A,I19)')
     &'  overall weighted active block length       :',ITOTBLCKL
       WRITE(LUWRT,'(A,I19)')
     &'  Maximum weighted block size                :',MXSIZE
C
C     final distribution
      I_AM_ACTIVE  = 0
      DO II = 1, NDIM
       IF( NPARBLOCK(II) .eq. -1 ) THEN
         WRITE(LUWRT,*) 'AAAHHHH, block',II,' is not distributed'
         WRITE(LUWRT,*) 'Since I do not know how to proceed,'// 
     &                  ' I will stop!'
         WRITE(LUWRT,*) ' *** ERROR in PRINT_BLOCK_DISTR ***'
           CALL QUIT(' *** ERROR in PRINT_BLOCK_DISTR : Block not 
     & distributed ***')
       END IF
      END DO
      WRITE(LUWRT,'(/20X,A)')'================================'
      WRITE(LUWRT,'(20X, A)')' Summation of even distribution '
      WRITE(LUWRT,'(20X,A/)')'================================'
      DO ISTI = 1, NMPROC
        ITLTND = 0
        ITLBND = 0
        I_AM_ACTIVE = 0
        DO IBLK = 1, NDIM
          IILENGTH = IBLOCKL(IBLK)
          IF(IILENGTH.gt.0)THEN
            I_AM_ACTIVE = I_AM_ACTIVE + 1
            IF(NPARBLOCK(IBLK).eq.ISTI-1)THEN
              IIWEIGHTBLK = ICWEIGHTF(I_AM_ACTIVE)
              ITLTND      = ITLTND + (IILENGTH*IIWEIGHTBLK)
              ITLBND      = ITLBND + 1
            END IF
          END IF
        END DO
        WRITE(LUWRT,'(2X,A,1X,I4,1X,A,1X,I6,1X,A,1X,I18)')
     &  'CPU',ISTI-1,' computes',ITLBND,'blocks with a total weight of',
     &      ITLTND
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE T_BLOCK_DISTR(IT_TTPL,IT_CTALL,IT_CTONE,I_APROC,NDIM,
     &                         IT_ACTIVE_R,IT_ACTIVE_I)
C
C     Block distribution routine for T-coefficient blocks  
C
C     Written by  S. Knecht         - February 01 2008
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
#include "mxpdim.inc"
#include "integrals_off.inc"
#include "ctcc.inc"
C
C     -----
C     INPUT
C     -----
C
C     NDIM       : total number of T-blocks
C     IT_CTALL(*): weight factors for each block on each node
C     IT_CTONE(*): weight factors for each block on the calling node
C
C     ------
C     OUTPUT
C     ------
C
C     IT_TTPL(NDIM): list of blocks containing the corresponding 
C                    assigned CPU
C
      DIMENSION IT_TTPL(*), IT_CTALL(*), IT_CTONE(*), I_APROC(*)
C
      INTEGER ICOUNTABLK, NPTEST, ITMPFAC
C     -------
C     SCRATCH
C     -------
C     MXSIZE: actual maximum size for one of all blocks
C     MXNUMB: number of the current largest block
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN, IOFF_CTALL
      INTEGER IMINUS2, IXTYPE, IACT_PART, IJKEY, IJCOLOR, IJ_LIST_COMM
      INTEGER*8 ITLTND,MXSIZE,MXSZTMP
      INTEGER*8 ITEMPL, ITEMPW1, ITEMPW2
      INTEGER*8 ILEN, ITOTBLCKL, IPRODLEN
      INTEGER*8 IIWEIGHTBLK, IILENGTH, IILENGTH_MAX
      INTEGER*8 NINT_TOT, nelm_cc
C
C     initialize    
C    
      NPTEST       = MAX(NPTEST_VAR,0)
      IONE         =  1
      IMIN1        = -1
      NINT_TOT     =  0
      IMINUS2      =  0
      IIWEIGHTBLK  =  0
      IILENGTH     =  0
      ITLTND       =  0
      ICOUNTABLK   =  0
      MXSIZE       =  0
      MXSZTMP      =  0
      MXNUMB       =  0
      IRUN         =  0
      ITEMPW1      =  0
      ITEMPW2      =  0
      ITEMPN       =  0
      ITEMPL       =  0
      IAMGPROC     =  0
      IAMTPROC     =  0
      ITOTBLCKL    =  0
      ILEN         =  0
      I_AM_ACTIVE  =  0
      IPRODLEN     =  0
      IILENGTH_MAX =  0
      IXTYPE       =  0
      IOFF_CTALL   =  0
      IACT_PART    =  0
      ITMPFAC      =  0
      IMIN2        = -2
      NINT_TOT     =  N1ELINT + N2ELINT
      NINT_TOTAL   =  NINT_TOT
      ICOUNTABLK   =  IT_ACTIVE_R + IT_ACTIVE_I
C
      CALL ISETVC(IT_TTPL,IMIN2,NDIM*IRC_SAVE)
*     print some information
      NPTEST = 0
C    
C     end of initialization
C
      DO IMULT = 1, IRC_SAVE
        IF( IMULT .eq. 1 ) IACT_PART = IT_ACTIVE_R
        IF( IMULT .eq. 2 ) IACT_PART = IT_ACTIVE_I
        DO ICBL = 1, IACT_PART
C
csk           WRITE(LUWRT,*) 'checking for active at',
csk     &     ICBL + ( NDIM * ( IMULT - 1 ) )
C
          IXTYPE = IT_CTONE( ICBL + ( NDIM * ( IMULT - 1 ) ) )
C
          IILENGTH = NELM_CC(IXTYPE,NDIM,NINT_TOT)
C
          IF( IILENGTH .gt. 0 )THEN
            I_AM_ACTIVE  = I_AM_ACTIVE + IONE
csk            WRITE(LUWRT,*) ' setting IXTYPE active at ',IXTYPE,
csk     &                       IXTYPE + ( NDIM * ( IMULT - 1 ) )
            IT_TTPL( IXTYPE + ( NDIM * ( IMULT - 1 ) ) ) = - 1
            IILENGTH_MAX = MAX(IILENGTH_MAX,IILENGTH)
          ELSE
            WRITE(LUWRT,*)' zero-length active T-block detected:',IXTYPE
            WRITE(LUWRT,*)' reset to inactive'
            IT_CTONE( ICBL + ( NDIM * ( IMULT - 1 ) ) ) = 0
          END IF
C         ^ IILENGTH > 0
        END DO
      END DO
C     ^ real/imag part
C
      ITEMPW1 = IILENGTH_MAX
      MXSIZE  = IILENGTH_MAX
C     reset
      IF( ICOUNTABLK .gt. I_AM_ACTIVE ) ICOUNTABLK = I_AM_ACTIVE
C
C
      IF(NPTEST.GE.1) THEN
         WRITE(LUWRT,'(A,I19)')
     &'  total number of processes to distribute on (group) :',
     &   NEWCOMM_PROC_SM
         WRITE(LUWRT,'(A,I19)')
     &'  total number of T-blocks                           :',
     &   NDIM * IRC_SAVE
         WRITE(LUWRT,'(A,I19)')
     &'  total number of active T-blocks                    :',
     &   ICOUNTABLK
         WRITE(LUWRT,'(A,I19)')
     &'  size of largest T-block                            :',
     &   IILENGTH_MAX
      END IF
*
      IMINNP = NEWCOMM_PROC_SM
      NTEMPP = IMINNP
      IF( ICOUNTABLK .lt. IMINNP ) THEN
        IMINNP = MIN( ICOUNTABLK,IMINNP )
        write(luwrt,'(/a,i6)') '*** error in t_block_distr:'//
     &    ' number of active T blocks lower than the total number of'//
     &    ' processes. please decrease the number of processes to:',
     &    IMINNP
        call quit('*** error in t_block_distr: number of processes >
     & number of active T blocks. the distribution algorithm
     & will therefore fail (see output for more information).***')
        NTEMPP = IMINNP
        IF( MYPROC .ge. ICOUNTABLK ) THEN 
          IAM_NOT_INV_T = 0
          IMINUS2 = - 2
          GOTO 8888
        END IF
      END IF

      IAM_NOT_INV_T = 1
C
C     starting the treausure quest for the ?optimal? T-block distribution
C
      DO IMULT = 1, IRC_SAVE
        IF( IMULT .eq. 1 ) IACT_PART = IT_ACTIVE_R
        IF( IMULT .eq. 2 ) IACT_PART = IT_ACTIVE_I
        ITEMPW1 = IILENGTH_MAX
        MXSIZE  = IILENGTH_MAX
        IRUN = 0
100     CONTINUE
C
        IRUN = IRUN + 1
csk        WRITE(LUWRT,*) ' IRUN, IMULT, IACT_PART',IRUN, IMULT, IACT_PART
        IF( IRUN .le. IACT_PART ) THEN
C
          MXSZTMP = 0
          I_AM_ACTIVE = 0
C
C         find largest T block
C
          DO 3000 II = 1, IACT_PART
C
            IXTYPE = 0
            IXTYPE = IT_CTONE( II + ( NDIM * ( IMULT - 1 ) ) )
C
            IF( IXTYPE .eq. 0 ) GOTO 3000
            IILENGTH = NELM_CC(IXTYPE,NDIM,NINT_TOT)
C
            IF( IILENGTH .gt. 0 ) THEN
C
              I_AM_ACTIVE = I_AM_ACTIVE + IONE
C
              ITEMPL = IILENGTH
C
              ITEMPD = IT_TTPL( IXTYPE + ( NDIM * ( IMULT - 1 ) ) )
              ITEMPN = IXTYPE
csk              WRITE(LUWRT,*) ' ITEMPD,IXTYPE,IXTYPE+(NDIM*(IMULT-1) )',
csk     &                         ITEMPD,IXTYPE,IXTYPE+(NDIM*(IMULT-1) )

              IF(ITEMPL.GT.0)THEN
                IF(ITEMPD.EQ.-1) THEN
                  IF(ITEMPL.LE.MXSIZE)THEN
                    IF(ITEMPL.GT.MXSZTMP) THEN
                       MXSZTMP = ITEMPL
                       MXNUMB = ITEMPN
                    END IF
                  END IF
                END IF
              END IF
            END IF
C           ^ IILENGTH > 0
3000      CONTINUE
C
C         find all CPUs asking for this T block
C
          CALL ISETVC(I_APROC,IMIN1,NTEMPP)
C
          ITMPFAC_MX  = 0
          DO 3500 IPR = 1, NTEMPP
C
            IOFF_CTALL = NDIM*IRC_SAVE*(IPR-1)+(NDIM*(IMULT-1))+MXNUMB
csk         IOFF_CTALL = IPR  * MXNUMB + ( NDIM * ( IMULT - 1 ) )
csk            WRITE(LUWRT,*) 'checking for ITMPFAC_MX at IOFF_CTALL',
csk     &                      IOFF_CTALL
            ITMPFAC = IT_CTALL( IOFF_CTALL )
            IF( ITMPFAC .gt. 0 ) I_APROC( IPR ) = ITMPFAC
            ITMPFAC_MX = MAX(ITMPFAC_MX,ITMPFAC)
C
3500      CONTINUE
csk          WRITE(LUWRT,*) ' ITMPFAC_MX for block MXNUMB',
csk     &                 ITMPFAC_MX,MXNUMB, MXNUMB + (NDIM * ( IMULT - 1))
C
C         pick out the CPU with current lowest T load
C         FIXME: take ITMPFAC_MX into account?
C
          IPROC_RUN = 0
C
          DO 4000 IPR = 1, NTEMPP
C
            ITMPFAC = I_APROC( IPR )
            ITEMPW2 = 0
            IF( ITMPFAC .le. 0) GOTO 4000
C
            DO IBLK = 1, IACT_PART
C            
              IXTYPE = 0
              IXTYPE = IT_CTONE( IBLK + ( NDIM * ( IMULT - 1 ) ) )
C
              IF( IXTYPE .gt. 0 ) THEN
                IILENGTH = NELM_CC(IXTYPE,NDIM,NINT_TOT)
C
                
                IF( IT_TTPL(IXTYPE+(NDIM*(IMULT - 1))).eq.IPR-1) THEN
CSK?              IOFF_CTALL = IPR * IXTYPE + ( NDIM * ( IMULT - 1 ) )
CSK?              IIWEIGHTBLK = IT_CTALL( IOFF_CTALL )
CSK?              ITEMPW2  = ITEMPW2 + IILENGTH * IIWEIGHTBLK
                  ITEMPW2 = ITEMPW2 + IILENGTH
                END IF
              END IF
            END DO

C         attention: IAMTPROC = 1 --> MASTER = 0 !
          IAMTPROC = IPR
C
          IPROC_RUN = IPROC_RUN + 1
          IF(IPROC_RUN .eq. 1 ) ITEMPW1 = ITEMPW2
C
          IF( ITEMPW2 .le. ITEMPW1 ) THEN
             ITEMPW1  = ITEMPW2
             IAMGPROC = IAMTPROC
          END IF
C
4000      CONTINUE
C         now we should have found a proc and a block --> put both together !
C
C         calculation of block MXNUMB by proc IAMGPROC -1
csk          IF( NPTEST .gt. 10 ) THEN
csk            WRITE(LUWRT,*) 'calculation of block MXNUMB by proc',
csk     &               MXNUMB, MXNUMB + (NDIM * ( IMULT - 1)) ,IAMGPROC-1
csk          ENDIF
C
          IAMTPROC = IAMGPROC - 1
          MXSIZE   = MXSZTMP
          IT_TTPL( MXNUMB + (NDIM * ( IMULT - 1) ) )     = IAMTPROC
C
        ELSE
          GOTO 101
        END IF
C     ^ IRUN !!!
C
        GOTO 100
C
101     CONTINUE
      END DO
C     ^ real/imag part
      I_AM_ACTIVE = 0
C
      IF(NPTEST.GE.0) THEN
        DO IMULT = 1, IRC_SAVE
          DO II = 1, NDIM
           IF( IT_TTPL( II + (NDIM * (IMULT - 1 ) ) ) .eq. -1 ) THEN
             WRITE(LUWRT,*) 'AAAHHHH, block',II+(NDIM * (IMULT - 1 ) ),
     &' is not distributed'
             WRITE(LUWRT,*) 'Since I do not know how to proceed, 
     &                       I will stop!'
             call quit('*** error in t_block_distr: active T block
     & not distributed among the processes.***')
            END IF
          END DO
        END DO
        WRITE(LUWRT,*)' '
        WRITE(LUWRT,'(20X,A)')
     &  '++++++++++++++++++++++++++++++++'
        WRITE(LUWRT,'(20X,A)')
     &  '+    T element distribution    +'
        IF( IT_SHL .eq. 0 )THEN
        WRITE(LUWRT,'(20X,A)')
     &  '+       (on MASTER node)       +'
        ELSE
        WRITE(LUWRT,'(20X,A)')
     &  '+           (global)           +'
        END IF
        WRITE(LUWRT,'(20X,A)')
     &  '++++++++++++++++++++++++++++++++'
        WRITE(LUWRT,*)' '
        DO ISTI = 1, NTEMPP
          ITLTND = 0
          ITLBND = 0
          DO IMULT = 1, IRC_SAVE
            DO IBLK = 1, NDIM
              IF( IT_TTPL( IBLK+(NDIM * (IMULT - 1 ) ) ).eq.ISTI -1)THEN
                  IILENGTH = NELM_CC(IBLK,NDIM,NINT_TOT) 
                  ITLTND = ITLTND + IILENGTH
                  ITLBND = ITLBND + 1
              END IF
            END DO
          END DO
          WRITE(LUWRT,'(2X,A,1X,I4,1X,A,1X,I5,1X,A,1X,I12)')
     &'CPU',ISTI-1,' allocates for',ITLBND,
     &' T blocks with a total length of', ITLTND
        END DO
      END IF
8888  CONTINUE
      IF( NTEMPP .lt. NEWCOMM_PROC_SM )THEN
        IF( IAM_NOT_INV_T .eq. 0 .or. MYNEW_ID_SM .eq. N_MASTER_SM )THEN
           IJKEY   = MYNEW_ID_SM
           IJCOLOR = 5
        ELSE 
           IJKEY   = MYNEW_ID_SM
           IJCOLOR = 6
        END IF
        call interface_mpi_COMM_SPLIT(MYNEW_COMM_SM,IJCOLOR,IJKEY,
     &                      IJ_LIST_COMM)
        IF( IJCOLOR .eq. 6 ) GOTO 8899
        call interface_mpi_bcast_i1_work_f77(IT_TTPL,NDIM*IRC_SAVE,
     &                 N_MASTER_SM,IJ_LIST_COMM)
8899    CONTINUE
        call interface_mpi_COMM_FREE(IJ_LIST_COMM)
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE CALC_OFF_MPI_FILE(FILENAME,IGROUPLIST,IBLOCKL,
     &                             IBLOCKD,IVEC_IN)
C
      use interface_to_mpi
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
C
      DIMENSION IGROUPLIST(*), IBLOCKL(*), IBLOCKD(*)
      CHARACTER*6 FILENAME
      INTEGER(KIND=df_MPI_OFFSET_KIND) ND_VEC_IOFF, ND_VEC_IOFF2
      INTEGER(KIND=df_MPI_OFFSET_KIND) MULT1, MULT2, JVEC
      INTEGER(KIND=df_MPI_OFFSET_KIND) I_VEC_LEN_SCR
C
      I_VEC_LEN_SCR = 0
      JVEC  = IVEC_IN
      MULT1 = 0
      MULT2 = 0
C
      IF( NPTEST_VAR .ge. 10 ) 
     & WRITE(LUWRT,*) 'Calculating offset for file', FILENAME
C
C     count active blocks for each cpu, calculate individual offset
      ND_VEC_IOFF = 0
      ND_VEC_IOFF2 = 0
      NBLK_ACT = 0
C
C     scratch file for c-vector
      IF( FILENAME .eq. 'ILUC  ' ) THEN
        DO JBLK = 1, NUM_BLOCKS
          IF( IBLOCKL(JBLK) .gt. 0 ) THEN
            NBLK_ACT = NBLK_ACT + 1
          END IF
        END DO
C
        MY_ACT_BLK_ALL = NUM_BLOCKS
        GOTO 101
      END IF
C
      NBLK_ACT = 0
C
C     we have group files, so only cpus within a group are relevant
      MY_VEC1_IOFF = 0
      MY_VEC2_IOFF = 0
C
      DO IPROC = 1, NEWCOMM_PROC
C     
        JPROC = IGROUPLIST(IPROC)
C
        DO JBLK = 1, NUM_BLOCKS
          IF(IBLOCKD(JBLK) .eq. JPROC ) THEN
C?            WRITE(LUWRT,*) 'assigning a block to JPROC',JBLK,JPROC
            ND_VEC_IOFF = ND_VEC_IOFF + IBLOCKL(JBLK)
            IF( MYPROC .eq. JPROC ) THEN
              MY_VEC1_IOFF = MY_VEC1_IOFF + IBLOCKL(JBLK) 
              MY_VEC2_IOFF = MY_VEC2_IOFF + IBLOCKL(JBLK) 
              NBLK_ACT = NBLK_ACT + 1
            END IF
          END IF
        END DO
C
        IF( MYPROC .eq. JPROC ) THEN
          MULT1 = IRC_SAVE 
C
C         MY_VEC1_IOFF = length of one vector
          MY_VEC1_IOFF = MY_VEC1_IOFF
C
C         MY_VEC1_IOFF = length of one vector x complex multiplier ( x 1 = real!)
          I_VEC_LEN_SCR = MY_VEC2_IOFF
          MY_VEC2_IOFF = I_VEC_LEN_SCR * MULT1
C
          MY_ACT_BLK1 = NBLK_ACT
          MY_ACT_BLK2 = NBLK_ACT * MULT1
C
          IF( FILENAME .eq. 'LUDIA ' ) THEN
C
C           initialize MY_DIA_OFF
            MY_DIA_OFF = 0
            MY_DIA_OFF = ND_VEC_IOFF2 * MULT1
          ELSE IF( FILENAME .eq. 'ILU1  ') THEN
C
C           initialize MY_LU1_OFF
            MY_LU1_OFF = 0
            MY_LU1_OFF = ND_VEC_IOFF2 * JVEC * MULT1
          ELSE IF( FILENAME .eq. 'ILU2  ') THEN
C
C           initialize MY_LU2_OFF
            MY_LU2_OFF = 0
            MY_LU2_OFF = ND_VEC_IOFF2 * JVEC * MULT1
          ELSE IF( FILENAME .eq. 'ILU3  ') THEN
C
C           initialize MY_LU3_OFF
            MY_LU3_OFF = 0
            MY_LU3_OFF = ND_VEC_IOFF2 * JVEC * MULT1
          ELSE IF( FILENAME .eq. 'ILU4  ') THEN
C
C           initialize MY_LU4_OFF
            MY_LU4_OFF = 0
            MY_LU4_OFF = ND_VEC_IOFF2 * JVEC * MULT1
          ELSE IF( FILENAME .eq. 'ILU5  ') THEN
C
C           initialize MY_LU5_OFF
            MY_LU5_OFF = 0
            MY_LU5_OFF = ND_VEC_IOFF2 * JVEC * MULT1
          END IF
C         ^ filenames
        END IF
C
C       reset ND_VEC_IOFF2 and ND_VEC_IOFF values
        ND_VEC_IOFF2 = ND_VEC_IOFF
        NBLK_ACT = 0
C
      END DO
C
      IF( NPTEST_VAR .ge. 10 ) THEN
        IF( FILENAME .eq. 'LUDIA ' ) THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'LUDIA offset :', MY_DIA_OFF
        ELSE IF( FILENAME .eq. 'ILU1  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU1  offset :', MY_LU1_OFF
        ELSE IF( FILENAME .eq. 'ILU2  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU2  offset :', MY_LU2_OFF
        ELSE IF( FILENAME .eq. 'ILU3  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU3  offset :', MY_LU3_OFF
        ELSE IF( FILENAME .eq. 'ILU4  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU4  offset :', MY_LU4_OFF
        ELSE IF( FILENAME .eq. 'ILU5  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU5  offset :', MY_LU5_OFF
        END IF
        WRITE(LUWRT,*) '  MY_VEC1_IOFF :', MY_VEC1_IOFF
        WRITE(LUWRT,*) '  MY_VEC2_IOFF :', MY_VEC2_IOFF
        WRITE(LUWRT,*) '  MY_ACT_BLK1  :', MY_ACT_BLK1
        WRITE(LUWRT,*) '  MY_ACT_BLK2  :', MY_ACT_BLK2
      END IF
C
 101  CONTINUE     
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CALC_OFF_MPI_FILE_XPRP(IGROUPLIST,IBLOCKL_X,IBLOCKD_X,
     &                                  IVECIN_ILU1,IVECIN_ILU2)
C**********************************************************************
C
C     calculate file offsets and position pointers for MPI-I/O files 
C     used in KR-CI property run.
C
C     find maximum vector length for a process within a group which 
C     serves as MY_VEC1_IOFF for all processes in this group.
C
C     written by S. Knecht - Nov 2008
C
C     last revision:
C
C**********************************************************************
C
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
#include "parluci.h"
#include "krciprop.h"
C
      DIMENSION IGROUPLIST(*), IBLOCKL_X(NUM_BLOCKS,*)
      DIMENSION IBLOCKD_X(NUM_BLOCKS,*)
      INTEGER NBLK_MAX_GP
      integer(kind=df_MPI_OFFSET_KIND) ivec_max_gp_tmp2, ivec_max_gp_tmp
      integer                       nblk_max_gp_tmp2, nblk_max_gp_tmp
      INTEGER(KIND=df_MPI_OFFSET_KIND) MULT1, MULT2, JVEC1, JVEC2
      INTEGER(KIND=df_MPI_OFFSET_KIND) LVECLEN, IVEC_MAX_GP, IPROC_OFF
C
      IVEC_MAX_GP      = 0
      NBLK_MAX_GP      = 0
      ivec_max_gp_tmp2 = 0
      nblk_max_gp_tmp2 = 0
      ivec_max_gp_tmp  = 0
      nblk_max_gp_tmp  = 0
      JVEC1            = IVECIN_ILU1
      JVEC2            = IVECIN_ILU2
      MULT1            = IRC_SAVE 
      MULT2            = 0
C
      LVECLEN          = 0
      NBLK_ACT         = 0
      MY_VEC1_IOFF     = 0
      MY_VEC2_IOFF     = 0
C     WRITE(LUWRT,*) '  NEWCOMM_PROC:', NEWCOMM_PROC
C
C     determine IVEC_MAX_GP and NBLK_MAX_GP
C     --------------------------------
      DO IPROC = 1, NEWCOMM_PROC
C       next process in this group
        JPROC = IGROUPLIST(IPROC)
        DO IJSYM = 1, IXSYMACT
C         WRITE(LUWRT,*) '  IJSYM:', IJSYM
C         next active symmetry irrep IJSYM (no absolute symmetry label)
          LVECLEN  = 0
          NBLK_ACT = 0
          DO JBLK = 1, NUM_BLOCKS
!           WRITE(6,*) '  IBLOCKD_X(JBLK,IJSYM):',
!    &                        IBLOCKD_X(JBLK,IJSYM)
            IF(IBLOCKD_X(JBLK,IJSYM) .eq. JPROC ) THEN
              LVECLEN  = LVECLEN + IBLOCKL_X(JBLK,IJSYM)
              NBLK_ACT = NBLK_ACT + 1
            END IF
          END DO
C         max vec length + block count for this symmetry irrep 
C         and process JPROC
          IVEC_MAX_GP_tmp2 = MAX(IVEC_MAX_GP_tmp2,LVECLEN)
          NBLK_MAX_GP_tmp2 = MAX(NBLK_MAX_GP_tmp2,NBLK_ACT)
!         WRITE(6,*) '  IVEC_MAX_GP,NBLK_MAX_GP:',
!    &                      IVEC_MAX_GP,NBLK_MAX_GP
        END DO
      END DO

!     make sure that all processes in this group obtain the correct
!     offset
      call interface_mpi_allreduce(ivec_max_gp_tmp2,ivec_max_gp_tmp,
     &                             1,op_mpi_max,mynew_comm)
     
      call interface_mpi_allreduce(nblk_max_gp_tmp2,nblk_max_gp_tmp,
     &                             1,op_mpi_max,mynew_comm)
      
      ivec_max_gp = ivec_max_gp_tmp
      nblk_max_gp = nblk_max_gp_tmp

!     WRITE(6,*) '  IVEC_MAX_GP,NBLK_MAX_GP:',
!    &              IVEC_MAX_GP,NBLK_MAX_GP,myproc
C
C     set MY_VEC1_IOFF, MY_VEC2_IOFF, MY_ACT_BLK1, MY_ACT_BLK2  
C     --------------------------------------------------------

C     MY_VEC1_IOFF = length of one vector (here: max per group)
      MY_VEC1_IOFF = IVEC_MAX_GP
      MY_VEC2_IOFF = IVEC_MAX_GP * MULT1
C     MY_ACT_BLK1  = number of active blocks (here: max per group)
      MY_ACT_BLK1  = NBLK_MAX_GP
      MY_ACT_BLK2  = NBLK_MAX_GP * MULT1
C
C     set MY_LU1_OFF and MY_LU2_OFF
      DO IPROC = 1, NEWCOMM_PROC
C       next process in this group
        JPROC      = IGROUPLIST(IPROC)
        IPROC_OFF  = IPROC - 1
        IF( MYPROC .eq. JPROC )THEN
          MY_LU1_OFF = 0
          MY_LU2_OFF = 0
          MY_LU1_OFF = IVEC_MAX_GP * JVEC1 * MULT1 * IPROC_OFF
          MY_LU2_OFF = IVEC_MAX_GP * JVEC2 * MULT1 * IPROC_OFF
        END IF
      END DO
C
      IF( NPTEST_VAR .ge. 10 ) THEN
        WRITE(LUWRT,'(2X,A,1X,I18,i2)') 'ILU1  offset :', MY_LU1_OFF,
     &  myproc
        WRITE(LUWRT,'(2X,A,1X,I18,i2)') 'ILU2  offset :', MY_LU2_OFF,
     &  myproc
        WRITE(LUWRT,*) '  MY_VEC1_IOFF :', MY_VEC1_IOFF, myproc
        WRITE(LUWRT,*) '  MY_VEC2_IOFF :', MY_VEC2_IOFF, myproc
        WRITE(LUWRT,*) '  MY_ACT_BLK1  :', MY_ACT_BLK1, myproc
        WRITE(LUWRT,*) '  MY_ACT_BLK2  :', MY_ACT_BLK2, myproc
      END IF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DIST_XCVEC(VEC1,LUC,LUOUT,LU1LIST,NPARBLOCK_X,
     &                      IBLOCKL_X,IXBLK,IXSYMLIST,IXSYMACT,IZ)
C**********************************************************************
C
C     distribute C vector blocks from MASTER to slaves for a given 
C     symmetry irrep.
C
C     written by S. Knecht - Nov 2008
C
C     last revision:
C
C**********************************************************************
C
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      DIMENSION VEC1(*), IXBLK(*), IXSYMLIST(*)
      DIMENSION LU1LIST(*), NPARBLOCK_X(NUM_BLOCKS,*)
      DIMENSION IBLOCKL_X(NUM_BLOCKS,*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) ILU1_IOFFX, MULTX
      LOGICAL EX
#include "parluci.h"
#include "mxgas.h"
#include "dcbkrci.h"
      ILU1_IOFFX = 0
      MULTX      = 0
      NVEC_TOT   = 0
      NVEC       = 0
C
C     update slaves of #roots per irrep
      CALL interface_MPI_BCAST(NKRCI_CIROOTS,MAX_NKRCI_MAX_SYM,
     &                         MASTER,global_communicator)
C
c     loop over active irreps
      DO IXA = 1, IXSYMACT
        IXSYM = IXSYMLIST(IXA)
        NVEC  = NKRCI_CIROOTS(IXSYM)
        IF( MYPROC .eq. MASTER .and. NVEC .gt. 0)
     &  CALL GET_LUC_FILE(LUC,IXSYM,EX)
        MULTX    = NVEC_TOT
C       offset for ILU1 and LU1LIST
        ILU1_IOFFX     = MY_LU1_OFF + MULTX * MY_VEC2_IOFF
        INT_ILU1_IOFFX = 1 + NVEC_TOT * MY_ACT_BLK2
        NBLOCK         = IXBLK(IXSYM)
C       WRITE(LUWRT,*) ' copy #vector(s) for sym',NVEC,IXSYM
!       WRITE(6,*) ' ILU1_IOFFX and INT_ILU1_IOFFX (MYPROC)',
!    &                   ILU1_IOFFX,INT_ILU1_IOFFX,MYPROC
C       WRITE(LUWRT,*) ' LUOUTLIST NOW',MYPROC
C       CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,LUWRT)
C       distribute C vector(s)
        CALL COP_REST_VEC_REL(VEC1,LUC,LUOUT,ILU1_IOFFX,
     &                        LU1LIST(INT_ILU1_IOFFX),
     &                        NPARBLOCK_X(1,IXA),IBLOCKL_X(1,IXA),
     &                        NBLOCK,NVEC,IZ)
C       close C vector file (opened in GET_LUC_FILE)
        IF( MYPROC .eq. MASTER .and. EX ) CLOSE(LUC,STATUS='KEEP')
C       keep track of correct offset
        NVEC_TOT = NVEC_TOT + NVEC
      END DO
!     WRITE(6,*) ' MY LUOUTLIST AFTER ALL',MYPROC
!     CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,6)
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE COLLECT_DENS_MAT(DENS_SCR,DENS,ISF,ISHARE)
C
C     Written by  S. Knecht         - June 27 2007
C
C**********************************************************************
C
C     collect density matrix elements from all CPUs to:
C
C                                           ISHARE = 0: MASTER  
C                                           ISHARE = 1: MASTER + NODES 
C
C     density matrix DENS for current eigenvector is in SIGDEN ordering
C     output: updated DENS
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "parluci.h"
*
      dimension DENS(*), DENS_SCR(*)
      INTEGER ISHCOMM, ISHCOLOR, ISHKEY
      INTEGER*8 NINT_TOT, NINT_TP_I8, IRIOFF, nelm_cc
*
      ISHCOLOR = 0
      ISHKEY   = 0
*
      NINT_TOT = N1ELINT + N2ELINT
*     check for quaternion algebra
      if (ISF.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
*     loop over real and imag part of density matrix elements
*
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
        IF( NPTESTVAR .GE. 10 )THEN
        write(LUWRT,*) '  '
        write(LUWRT,*) '  '
        write(LUWRT,'(2X,A)') 
     &            'updating one- and two-particle densities '
        write(LUWRT,*) '  '
        write(LUWRT,'(2X,A,1X,I3)')' parallel update level ', ISHARE
        write(LUWRT,*) '  '
        if (IRI.eq.1) write(LUWRT,*) ' Real density elements     : '
        if (IRI.eq.2) write(LUWRT,*) ' Imaginary density elements: '
        write(LUWRT,*) ' -------------------------------'
        END IF
*       
*       loop over excitation types
*
        do IXTP = 1,NSPOBEX_TP,1
CSK          write(LUWRT,*) '  '
CSK          write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
CSK          write(LUWRT,*)        ' -------------------------'
CSK          write(LUWRT,*) '  '
             NINT_TP_I8 = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
             XNINT_TP_I8 = REAL(NINT_TP_I8)
             NINT_TP = IGIVE_I_B(XNINT_TP_I8)
*
csk          write(LUWRT,*) ' NINT_TP =',NINT_TP
*
          CALL DZERO(DENS_SCR,NINT_TP)
*
*
          IF( ISHARE .eq. 0 )THEN          
*
            XXX = 0.0D0
            XXX = DDOT(NINT_TP,DENS(IRIOFF+IRUN),1,DENS(IRIOFF+IRUN),1)
*
            IF( MYPROC .eq. MASTER )THEN
              ISHCOLOR = 7
            ELSE
              IF( XXX .eq. 0.0D0 ) ISHCOLOR = 6
              IF( XXX .ne. 0.0D0 ) ISHCOLOR = 7
            END IF
*           key value: MASTER will be MASTER again
            ISHKEY = MYPROC
*           new communicator for REDUCE on MASTER
            call interface_mpi_COMM_SPLIT(global_communicator,
     &                          ISHCOLOR,ISHKEY,ISHCOMM)
*
            IF( ISHCOLOR .eq. 6 ) GOTO 100
*
          END IF
*
*         start update 
          IF( ISHARE .eq. 0 )THEN
*           sum up partial density matrix elements on MASTER
            CAll redvec(DENS(IRIOFF+IRUN),DENS_SCR,NINT_TP,2,
     &                      op_MPI_SUM,ISHCOMM,MASTER)
*
*           copy back updated matrix elements
            CALL DCOPY(NINT_TP,DENS_SCR,1,DENS(IRIOFF+IRUN),1)
*
          ELSE IF( ISHARE .eq. 1 )THEN
*           sum up partial density matrix elements on ALL CPUs
            CAll redvec(DENS(IRIOFF+IRUN),DENS_SCR,NINT_TP,2,
     &                      op_MPI_SUM,global_communicator, - 1 )
*
*           copy back updated matrix elements
            CALL DCOPY(NINT_TP,DENS_SCR,1,DENS(IRIOFF+IRUN),1)
          END IF
*         ^ update either MASTER or all CPUs
  100     CONTINUE
*
          IF( ISHARE .eq. 0 ) call interface_mpi_COMM_FREE(ISHCOMM)
*
*         keep track of correct offset
          IRUN = IRUN + NINT_TP
        end do
*       ^ loop over excitation types
      end do
*     ^ loop over real and imag part of density matrix elements
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE COLLECT_DENS_MAT_SPLIT(DENS_SCR,DENS,IT_TTPL,ISF)
C
C     Written by  S. Knecht         - May 28 2008
C
C**********************************************************************
C
C     collect density matrix elements
C
C     density matrix DENS for current eigenvector is in SIGDEN ordering
C     output: updated DENS
C
C     Last revision:     S. Knecht       - May 2008
C
************************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "parluci.h"
#include "krmc_shmem.h"
*
      dimension DENS(*), DENS_SCR(*), IT_TTPL(*)
      INTEGER ISHCOMM, ISHCOLOR, ISHKEY
      INTEGER*8 NINT_TOT, ITB_BUFF, ITB_BUFF_ADD, nelm_cc
*
      ISHCOLOR     = 0
      ISHKEY       = 0
      ITB_BUFF_ADD = 0
      ID_GET       = - 1
*
      NINT_TOT = N1ELINT + N2ELINT
csk   WRITE(LUWRT,*) ' NINT_TOT, N1ELINT, N2ELINT',
csk  &                 NINT_TOT, N1ELINT, N2ELINT
*     check for quaternion algebra
      if (ISF.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
*     loop over real and imag part of density matrix elements
*
*
      do IRI = 1,IRILP,1
        ITB_BUFF = 1
        do JRI = 1,IRI-1,1
          ITB_BUFF = ITB_BUFF + LEN_T_BUFF
        end do
csk     write(LUWRT,*) '  ITB_BUFF, LEN_T_BUFF',ITB_BUFF, LEN_T_BUFF
        IF( NPTESTVAR .GE. 10 )THEN
        write(LUWRT,*) '  '
        write(LUWRT,*) '  '
        write(LUWRT,'(2X,A)') 
     &            'updating one- and two-particle densities '
        write(LUWRT,*) '  '
        if (IRI.eq.1) write(LUWRT,*) ' Real density elements     : '
        if (IRI.eq.2) write(LUWRT,*) ' Imaginary density elements: '
        write(LUWRT,*) ' -------------------------------'
        END IF
*       
*       loop over excitation types
*
        do IXTP = 1,NSPOBEX_TP,1
csk        write(LUWRT,*) '  '
csk        write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
csk        write(LUWRT,*)        ' -------------------------'
csk        write(LUWRT,*) '  '
           IXXTTP = (( IRI - 1 ) * NSPOBEX_TP) + IXTP
C          active T block? (complex case included)
           ID_GET = IT_TTPL( IXXTTP )
csk        WRITE(6,*) ' type IXXTTP and ID_GET',IXXTTP,ID_GET,MYPROC
           IF( ID_GET .le. 0 ) THEN
             ITB_BUFF_ADD = 0
             NINT_TP      = 0
             ISHCOLOR     = 6
           ELSE
C            length of T block to compute T_BUFF offset
             ITB_BUFF_ADD = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
             XNINT_TP_I8  = REAL(ITB_BUFF_ADD)
             NINT_TP      = IGIVE_I_B(XNINT_TP_I8)
             XXX = 0.0D0
             XXX = DDOT(NINT_TP,DENS(ITB_BUFF),1,DENS(ITB_BUFF),1)
             ISHCOLOR     = 7
             IF( XXX .eq. 0.0D0 ) THEN 
                ISHCOLOR        = 6
                NINT_TP         = 0
                IT_TTPL(IXXTTP) = - 1
             END IF
           END IF
*
           CALL DZERO(DENS_SCR,NINT_TP)
*
           ISHKEY = MYPROC
*          new communicator for REDUCE on MASTER
           call interface_mpi_COMM_SPLIT(global_communicator,
     &                         ISHCOLOR,ISHKEY,ISHCOMM)
*
*          "global" update  
*
           CAll redvec(DENS(ITB_BUFF),DENS_SCR,NINT_TP,2,
     &                     op_MPI_SUM,ISHCOMM, - 1 )
           CALL DCOPY(NINT_TP,DENS_SCR,1,DENS(ITB_BUFF),1)
*
  100      CONTINUE
           call interface_mpi_COMM_FREE(ISHCOMM)
*
*          keep track of correct DENS offset
           ITB_BUFF = ITB_BUFF + ITB_BUFF_ADD
        end do
*       ^ loop over excitation types
      end do
*     ^ loop over real and imag part of density matrix elements
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE COLLECT_DENS_MAT_SM_LOC(DENS,DENS_SCR,IT_TTPL,IT_TTOL,
     &                                   ISF,ISHARE)
C
C     Written by  S. Knecht         - February 26 2008
C
C**********************************************************************
C
C     collect density matrix elements from all CPUs to:
C
C                                           ISHARE = 0: MASTER  
C
C     IT_TTPL: T block list (group internal CPU)
C     IT_TTOL: T block window offset list (group internal)
C
C     density matrix DENS for current eigenvector is in SIGDEN ordering
C     output: updated DENS
C
C     Last revision:     S. Knecht       - February  2008
C
************************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "parluci.h"
*
      DIMENSION DENS(NALLINT), DENS_SCR(*), IT_TTPL(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IT_TTOL(*)
      INTEGER ISHCOMM, ISHCOLOR, ISHKEY
      INTEGER*8 NINT_TOT, NINT_TP_I8,ITB_I8, nelm_cc
*
      ISHCOLOR =   0
      ISHKEY   =   0
      ITB_I8   =   0
      ID_GET   = - 1
*
      NINT_TOT = N1ELINT + N2ELINT
*     check for quaternion algebra
      if (ISF.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
C
C     loop over real and imag part of density matrix elements
C
      do IRI = 1,IRILP,1
C
        IF( NPTESTVAR .GE. 10 )THEN
          write(LUWRT,*) '  '
          write(LUWRT,*) '  '
          write(LUWRT,'(2X,A)') 
     &            'updating one- and two-particle densities '
          write(LUWRT,*) '  '
          write(LUWRT,'(2X,A,1X,I3)')' parallel update level ', ISHARE
          write(LUWRT,*) '  '
          if (IRI.eq.1) write(LUWRT,*) ' Real density elements     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary density elements: '
          write(LUWRT,*) ' -------------------------------'
        END IF
C       
C       loop over excitation types
C
        do IXTP = 1,NSPOBEX_TP,1
C
           NINT_TP_I8  = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
           XNINT_TP_I8 = REAL(NINT_TP_I8)
           NINT_TP = IGIVE_I_B(XNINT_TP_I8)
C
C          determine T type for use in T list
           IXXTTP = IXTP
           IF( IRI .eq. 2 ) IXXTTP = ITTP + NSPOBEX_TP
C          process id and window offset (complex case included)
           ID_GET = IT_TTPL(IXXTTP)
C
           IF( ID_GET .ge. 0 )THEN
             IF( ISHARE .eq. 0 )THEN
C
               CALL DZERO(DENS_SCR,NINT_TP)
               IF( ID_GET .eq. MYNEW_ID_SM )THEN
                 ISHCOLOR = 7
               ELSE
                 ISHCOLOR = 6
               END IF
C              key value: MASTER will be MASTER again
               ISHKEY = MYPROC
C              new communicator for REDUCE task
               call interface_mpi_COMM_SPLIT(global_communicator,
     &                             ISHCOLOR,ISHKEY,ISHCOMM)
C
               IF( ISHCOLOR .eq. 6 ) GOTO 100
C              window offset
               ITB_I8 = IT_TTOL(IXXTTP)
C
C              start update 
C              ------------
C
C              1. sum up partial density matrix elements
C              ---------------------------------------------------
               CAll redvec(DENS(ITB_I8+1),DENS_SCR,NINT_TP,2,
     &                         op_MPI_SUM,ISHCOMM,-1)
C
C              2. copy back updated matrix elements
C              ------------------------------------
               CALL DCOPY(NINT_TP,DENS_SCR,1,DENS(ITB_I8+1),1)
  100          CONTINUE
C
               call interface_mpi_COMM_FREE(ISHCOMM)
             END IF
           END IF
C
        end do
C       ^ loop over excitation types
      end do
C     ^ loop over real and imag part of density matrix elements
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE DISTBLKND_REL(NDIM,NBLOCKL,NPARBLOCK)
*
*. Find distribution for (nonvanishing blocks) among the nodes
*. Non vanishing block has a nonzero blocklength
*
      use interface_to_mpi
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
*
      INTEGER*8 NPARBLKWT
      DIMENSION NBLOCKL(NDIM), NPARBLKWT(2,NMPROC)
*. Scratch
      INTEGER INONVAN,MAXLBL,LBL,IWEIGHT,LABEL,IBLOCKN,IMINNP
      INTEGER MXSZTMP,MXNUMB,ITEMPL,ITEMPD
      INTEGER*8 ITOTBLL 
*. ======
*. Output
*. ======
*. array NPARBLOCK that contains the final block distribution
      DIMENSION NPARBLOCK(NDIM)
*. 
      MXSIZE = 1000000000
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 1000000000
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      INONVAN = 0
      ITOTBLL = 0
      IMINNP = 0
      MAXLBL = 0 
      LABEL = 0
      MIN2 = -2
      NZERO = 0
      NTEST = 0
      CALL ISETVC(NPARBLOCK,MIN2,NDIM)
      DO 100 II = 1, NDIM
        LBL = NBLOCKL(II)
        IF(LBL.GT.0) THEN
          INONVAN = INONVAN + 1
          NPARBLOCK(II) = -1
          ITOTBLL = ITOTBLL + LBL
        ENDIF
        IF(LBL.GT.MAXLBL) THEN
          MAXLBL = LBL
          LABEL = II
        END IF
 100  CONTINUE
      DO II = 1,2
        DO JJ = 1, NMPROC
          NPARBLKWT(II,JJ) = 0
        END DO
      END DO
*
*
      NTEST = 10
      IF(NTEST.GT.0) THEN
         WRITE(LUWRT,*) '  '
         WRITE(LUWRT,*) '  '
         WRITE(LUWRT,*) '  total number of processes to distribute on:',
     &                 NMPROC
         WRITE(LUWRT,*) '  total number of active blocks :', INONVAN
         WRITE(LUWRT,*) '  overall active block length   :', ITOTBLL 
         WRITE(LUWRT,*) '  number of largest active block:', LABEL
         WRITE(LUWRT,*) '  size of largest active block  :', MAXLBL
         WRITE(LUWRT,*) '  '
      END IF
*
      IF(INONVAN.LT.NMPROC) THEN
        IMINNP = MIN(INONVAN,NMPROC)
        write(luwrt,'(/a,i6)') '*** error in distblknd_rel:'//
     &    ' number of active blocks lower than the total number of'//
     &    ' processes. please decrease the number of processes to:',
     &    IMINNP
        write(luwrt,'(/a)') ' alternatively change the GAS'//
     &    ' specification (more GA spaces) to increase the number of'//
     &    ' TTSS blocks.***'
        call quit('*** error in distblknd_rel:: number of processes >
     & number of active TTSS blocks. the distribution algorithm
     & will therefore fail (see output for more information).***')
        IF(MYPROC.GE.INONVAN) THEN
          IAM_NOT_INV = 0
          GOTO 1001
        END IF
      ELSE
        IMINNP = NMPROC
      END IF
      IAM_NOT_INV = 1
*
*. starting the treausure quest for the ?optimal? c-block distribution
*
 200  CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.INONVAN) THEN
*
        MXSZTMP = 0
*
        DO 300 II = 1, NDIM
*
          ITEMPL = NBLOCKL(II)
          ITEMPD = NPARBLOCK(II)
          ITEMPN = II
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.EQ.-1) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GE.MXSZTMP) THEN
                  MXSZTMP = ITEMPL
                  MXNUMB = II
                END IF
              END IF
            END IF
          END IF
*
300     CONTINUE
*
        DO 400 IPR = 1, IMINNP
*
          ITEMPW2 = NPARBLKWT(2,IPR)
* attention: IAMTPROC = 1 <--> MYPROC = 0 !
          IAMTPROC = IPR
*
          IF(IPR.EQ.1) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
            ITEMPW1 = ITEMPW2
            IAMGPROC = IAMTPROC
          END IF
*
 400    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*. calculation of block MXNUMB by proc IAMGPROC
        IF(NTEST.GT.0) THEN
          IF(MYPROC.EQ.MASTER) THEN
            WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,
     &      IAMGPROC-1
          END IF
        ENDIF
*. raising NPARBLKWT(1,proc) by 1
*. adding on NPARBLKWT(2,proc) the weight of the new block
        MXSIZE = MXSZTMP
        NPARBLKWT(1,IAMGPROC) = NPARBLKWT(1,IAMGPROC) + 1
        NPARBLKWT(2,IAMGPROC) = NPARBLKWT(2,IAMGPROC) + MXSZTMP 
        NPARBLOCK(MXNUMB) = IAMGPROC - 1
      ELSE
        GOTO 500
      END IF
*     /\ IRUN !!!
      GOTO 200
*
500   CONTINUE
      IF(NTEST.GE.0) THEN
          WRITE(LUWRT,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF 
     &    MYPROC:',MYPROC
          DO ISTI = 1, IMINNP
            WRITE(LUWRT,'(3X,A,I17,A,I17,A,1X,I17)')
     &      'process',ISTI -1,' calculates',NPARBLKWT(1,ISTI),
     &      ' blocks with a total length of',NPARBLKWT(2,ISTI)
          END DO
        END IF
1001  CONTINUE
*      
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FIND_CCTOS_REL(ISCALFAC,IBLOCKD,ICCTOS,IBLOCKL_S,
     &                          IBLOCKL_C,I_RUN_COMPLEX,NBLOCK)
C***********************************************************************
C
C     Find all c-blocks connecting to a given sigma-block for each cpu
C     using the connection matrix ICCTOS.
C     Each cpu stores the information in ISCALFAC.
C     Connection is marked by 1
C
C     written by S. Knecht - Jun 2007
C
C     Last revision : S. Knecht - Oct 2008 
C                   - adaption for property calculations
C
C***********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
C     ==========
C       INPUT
C     ==========
C                      CBLOCKS    , SBLOCKS
      DIMENSION ICCTOS(I_NZERO_LEN_C,I_NZERO_LEN_S), IBLOCKD( NBLOCK )
      DIMENSION IBLOCKL_S(NBLOCK), IBLOCKL_C(NBLOCK)
C
C     ===========
C       OUTPUT
C     ===========
C
C     NOTE: IF NZ == 2 (complex quaternion algebra needed):
C     DIMENSION ISCALFAC(2*NBLOCK) else ISCALFAC(NBLOCK)
      DIMENSION ISCALFAC(*)
      IONE          = 1
      IPUT_BLK      = 0
C
      DO II_RC = 1, I_RUN_COMPLEX
C
        I_AM_ACTIVE_S = 0
        I_AM_ACTIVE_C = 0
C
        DO 100 IBLK = 1, NBLOCK
C
          IF( IBLOCKL_S(IBLK) .eq. 0 ) GOTO 100
C
          I_AM_ACTIVE_S = I_AM_ACTIVE_S + IONE
C
          IF( IBLOCKD(IBLK) .ne. MYPROC ) GOTO 100
            I_AM_ACTIVE_C = 0
            DO 200 JBLK = 1,NBLOCK
C
              IF( IBLOCKL_C(JBLK) .eq. 0 ) GOTO 200 
              I_AM_ACTIVE_C = I_AM_ACTIVE_C + IONE
C
              IPUT_BLK = JBLK + (NBLOCK *(II_RC-1))
C
              IF( ICCTOS(I_AM_ACTIVE_C,I_AM_ACTIVE_S) .ne. 0 ) 
     &            ISCALFAC(IPUT_BLK) = IONE
C
 200        CONTINUE
C
 100      CONTINUE
      ENDDO
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_GROUP_OF_PROCS_REL(MYOWNNAME,MYNAMELENGTH,
     &                                   IPROCLIST,NAMELENGTHLIST)
C
      use interface_to_mpi
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
C
      DIMENSION IPROCLIST(NMPROC),NAMELENGTHLIST(NMPROC)
      CHARACTER*255 PNAMELIST(NMPROC),MYOWNNAME,MYSCRNAME
      CHARACTER*255 OLD_NAME
      INTEGER INUMBER, IFINI, NTEST,ICOUNT_GROUPS
C
      NTEST = 000
      INUMBER = 0
      IFINI = 0
      ICOUNT_GROUPS = 0
C
C     start with an all_gather of all name length
C
      call interface_mpi_ALLGATHER(MYNAMELENGTH,1,
     &                   NAMELENGTHLIST,1,
     &                   global_communicator)
C
C     ... all name length are stored now!
C     ... save own name in scratch array!
C
      PNAMELIST(MYPROC+1) = MYOWNNAME
C
C     ... we start to gather the names ... 
      DO IPROC = 1, NMPROC
C
         MYSCRNAME = MYOWNNAME
         call interface_mpi_BCAST(MYSCRNAME,NAMELENGTHLIST(IPROC),
     &                 IPROC-1,global_communicator)
         IF( MYPROC .NE. IPROC -1 ) THEN
           PNAMELIST(IPROC) = MYSCRNAME(1:NAMELENGTHLIST(IPROC))
         END IF
C
      END DO
C
C     test writing if needed
C
      IF( NTEST .GE. 100 ) THEN
       DO IPROC = 1, NMPROC
        MYSCRNAME = PNAMELIST(IPROC)(1:NAMELENGTHLIST(IPROC))
        WRITE(LUWRT,'(2X,A,1X,I4,1X,A7,1X,A)')'Processor ',IPROC -1,
     & ' alias ',MYSCRNAME(1:NAMELENGTHLIST(IPROC))
C    & MYSCRNAME
       END DO
      END IF
C
C     find all processors on the same deck and reorder (if necessary) 
C     to get the processors as close as possible, starting with the
C     master 
C
      INUMBER = MASTER + 1
 100  CONTINUE
      IFINI = 0
      ICOUNT_GROUPS = ICOUNT_GROUPS + 1
      OLD_NAME = PNAMELIST( INUMBER)(1:NAMELENGTHLIST( INUMBER))
C
      DO IPROC = 1, NMPROC
C
        MYSCRNAME(1:NAMELENGTHLIST(IPROC)) = 
     &  PNAMELIST(IPROC)(1:NAMELENGTHLIST(IPROC))
C
        IF(MYSCRNAME(1:NAMELENGTHLIST(IPROC)) .eq. 
     &     OLD_NAME(1:NAMELENGTHLIST(INUMBER)))
     &     IPROCLIST(IPROC) = INUMBER
C
      END DO
C
C     if finished: IFINI = 1, else 0
C     search for the next lowest cpu building the 'local group master'
C
      IFINI = 1
      DO IPRC = 1, NMPROC
        IF(IFINI .ne. 0 .and. IPROCLIST(IPRC) .eq. -1) THEN
          IFINI = 0
          INUMBER = IPRC
        END IF
      END DO
      IF(IFINI .eq. 0 ) GOTO 100
C
C     write the count of 'groups' formed by processes
C
      IF( NTEST .ge. 10) THEN
        WRITE(LUWRT,'(2X,I4,1X,A)') ICOUNT_GROUPS,'group(s) of 
     & processors found! '
        IF( NTEST .ge. 50) THEN
          WRITE(LUWRT,'(2X,A)') 'Complete listing:'
          WRITE(LUWRT,*) '  '
          CALL IWRTMAMN(IPROCLIST,1,NMPROC,1,NMPROC,LUWRT)
        END IF
        WRITE(LUWRT,*) '  '
      END IF
C
C     Transfer information to common block
C
      NFLGRPS_REL = ICOUNT_GROUPS
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_GROUP_OF_PROCS_REL_SM(MYOWNNAME,MYNAMELENGTH,
     &                                     IPROCLIST,NAMELENGTHLIST)
C
      use interface_to_mpi
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
C
      DIMENSION IPROCLIST(NMPROC),NAMELENGTHLIST(NMPROC)
      CHARACTER*255 PNAMELIST(NMPROC),MYOWNNAME,MYSCRNAME
      CHARACTER*255 OLD_NAME
      INTEGER INUMBER, IFINI, NTEST,ICOUNT_GROUPS
C
      NTEST =  0
csk      NTEST = 100
      INUMBER = 0
      IFINI = 0
      ICOUNT_GROUPS = 0
C
C
C     start with an all_gather of all name length
C
      call interface_mpi_ALLGATHER(MYNAMELENGTH,1,
     &                   NAMELENGTHLIST,1,
     &                   global_communicator)
C
C     ... all name length are stored now!
C     ... save own name in scratch array!
C
      PNAMELIST(MYPROC+1) = MYOWNNAME
C
C     ... we start to gather the names ... 
      DO IPROC = 1, NMPROC
C
         MYSCRNAME = MYOWNNAME
         call interface_mpi_BCAST(MYSCRNAME,NAMELENGTHLIST(IPROC),
     &                 IPROC-1,global_communicator)
         IF( MYPROC .NE. IPROC -1 ) THEN
           PNAMELIST(IPROC) = MYSCRNAME(1:NAMELENGTHLIST(IPROC))
         END IF
C
      END DO
C
C     test writing if needed
C
      IF( NTEST .GE. 100 ) THEN
       DO IPROC = 1, NMPROC
        MYSCRNAME = PNAMELIST(IPROC)(1:NAMELENGTHLIST(IPROC))
        WRITE(LUWRT,'(2X,A,1X,I4,1X,A7,1X,A)')'Processor ',IPROC -1,
     & ' alias ',MYSCRNAME(1:NAMELENGTHLIST(IPROC))
C    & MYSCRNAME
       END DO
      END IF
C
C     find all processors on the same deck and reorder (if necessary) 
C     to get the processors as close as possible, starting with the
C     master 
C
      INUMBER = MASTER + 1
 100  CONTINUE
      IFINI = 0
      ICOUNT_GROUPS = ICOUNT_GROUPS + 1
      OLD_NAME = PNAMELIST( INUMBER)(1:NAMELENGTHLIST( INUMBER))
C
      DO IPROC = 1, NMPROC
C
        MYSCRNAME(1:NAMELENGTHLIST(IPROC)) = 
     &  PNAMELIST(IPROC)(1:NAMELENGTHLIST(IPROC))
C
        IF(MYSCRNAME(1:NAMELENGTHLIST(IPROC)) .eq. 
     &     OLD_NAME(1:NAMELENGTHLIST(INUMBER)))
     &     IPROCLIST(IPROC) = INUMBER
C
      END DO
C
C     if finished: IFINI = 1, else 0
C     search for the next lowest cpu building the 'local group master'
C
      IFINI = 1
      DO IPRC = 1, NMPROC
        IF(IFINI .ne. 0 .and. IPROCLIST(IPRC) .eq. -1 ) THEN
          IFINI = 0
          INUMBER = IPRC
        END IF
      END DO
      IF(IFINI .eq. 0 ) GOTO 100
C
C     write the count of 'groups' formed by processes
C
      IF( NTEST .ge. 10) THEN
        WRITE(LUWRT,'(2X,I4,1X,A)') ICOUNT_GROUPS,'group(s) of 
     & processors found! '
        IF( NTEST .ge. 50) THEN
          WRITE(LUWRT,'(2X,A)') 'Complete listing:'
          WRITE(LUWRT,*) '  '
          CALL IWRTMAMN(IPROCLIST,1,NMPROC,1,NMPROC,LUWRT)
        END IF
        WRITE(LUWRT,*) '  '
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_SENDER(ID_LIST,ID_LIST_ALL,NDIM)
*
*  Find density block sender == lowest processor tag 
*
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION ID_LIST(*), ID_LIST_ALL(*)
C
      DO IMULT = 1, IRC_SAVE
         JMULT = (IMULT-1)*NDIM
        DO IXTP = 1, NDIM
          NLOW = 0
 100      CONTINUE
          NLOW = NLOW + 1
          IF(ID_LIST_ALL(((NLOW-1)*NDIM*IRC_SAVE)+IXTP+JMULT).gt.0)THEN
            ID_LIST(IXTP+JMULT) = NLOW - 1
            GOTO 200
          ELSE
            IF(NLOW .le. NMPROC-1) THEN 
               GOTO 100
            ELSE
               ID_LIST(IXTP+JMULT) = -1
               GOTO 200
            END IF
          END IF
 200      CONTINUE
        END DO
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GROUP_DESTRUCTOR_REL(JCOMM1,JCOMM2,JCOMM3,JCOMM4)
C
      use interface_to_mpi
      IMPLICIT REAL*8 (A-H,O-Z)
C
      call interface_mpi_COMM_FREE(JCOMM1)
      call interface_mpi_COMM_FREE(JCOMM3)
      call interface_mpi_COMM_FREE(JCOMM2)
      call interface_mpi_COMM_FREE(JCOMM4)
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      INTEGER*8 FUNCTION NCALC_MAX_BLK(ISF)
C
C     Written by  S. Knecht         - October 1 2007
C
C**********************************************************************
C
C     calculate max. density matrix elements block:
C
C     Last revision:     S. Knecht       - October 1 2007
C
************************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "parluci.h"
*
      INTEGER*8 NINT_TOT, NINT_TP, MAX_SCR_ELEMENT
      MAX_SCR_ELEMENT = 0
      NINT_TOT = N1ELINT + N2ELINT
*
*     check for quaternion algebra
*
      if (ISF.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
*     loop over real and imag part of density matrix elements
*
      do IRI = 1,IRILP,1
*       
*       loop over excitation types
*
        do IXTP = 1,NSPOBEX_TP,1
CSK          write(LUWRT,*) '  '
CSK          write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
CSK          write(LUWRT,*)        ' -------------------------'
CSK          write(LUWRT,*) '  '
          NINT_TP = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          MAX_SCR_ELEMENT = MAX(MAX_SCR_ELEMENT,NINT_TP)
*
csk          write(LUWRT,*) ' NINT_TP =',NINT_TP
csk          write(LUWRT,*) ' MAX_SCR_ELEMENT =',MAX_SCR_ELEMENT
*
        end do
*       ^ loop over excitation types
      end do
*     ^ loop over real and imag part of density matrix elements
      NCALC_MAX_BLK = MAX_SCR_ELEMENT
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SETUNITS_PAR_CLOSE_REL(CIRUN,JTYPE)
C
C     OUTPUT
C     ======
C
C     close files
C     files handles stored on common block LUCIPFILE are set to MPI_NULL
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      CHARACTER*6 CIRUN
C
C     ... close
C
      IF( JTYPE .gt. 0 )THEN
        IF( CIRUN .eq. 'KR-CI ')THEN
C
          IF( JTYPE .eq. 1 )THEN
            call interface_mpi_FILE_CLOSE(IDIA)
            call interface_mpi_FILE_CLOSE(ILU2)
            call interface_mpi_FILE_CLOSE(ILU3)
            call interface_mpi_FILE_CLOSE(ILU4)
            call interface_mpi_FILE_CLOSE(ILU5)
            call interface_mpi_FILE_CLOSE(ILU7)
            call interface_mpi_FILE_CLOSE(ILUC)
C
          ELSE IF( JTYPE .eq. 2 )THEN
            call interface_mpi_FILE_CLOSE(ILU1)
          END IF
        END IF
      ELSE
        call interface_mpi_FILE_CLOSE(IDIA)
        call interface_mpi_FILE_CLOSE(ILU1)
        call interface_mpi_FILE_CLOSE(ILU2)
        call interface_mpi_FILE_CLOSE(ILU3)
        call interface_mpi_FILE_CLOSE(ILU4)
        call interface_mpi_FILE_CLOSE(ILU5)
        call interface_mpi_FILE_CLOSE(ILU7)
        call interface_mpi_FILE_CLOSE(ILUC)
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SETUNITS_PAR_OPEN_REL(NFILE_ID)
C
C     OUTPUT
C     ======
C
C     open files ready for MPI-I/O
C     files handles stored on common block LUCIPFILE
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
#include "mxgas.h"
#include "dcbkrmc.h"
      CHARACTER*6 PLU1BL,PLU2BL,PLU3BL,PLU4BL,PLU7BL,PLU5BL, PLUCBL
      CHARACTER*10 PDIA, PLU1, PLU2, PLU3, PLU4, PLU5, PLU6, PLU7, PLUC
      CHARACTER*6 PDIABL 
      CHARACTER*4 FILE_INFO_GROUPSZ
      INTEGER PFILELENGTH
      INTEGER(KIND=df_MPI_OFFSET_KIND) IDISP, IDIA_SIZE
      IDISP = 0
      PFILELENGTH = 0
C
      PDIABL = 'DIAPAR'
      PLU1BL = 'LU1PAR'
      PLU2BL = 'LU2PAR'
      PLU3BL = 'LU3PAR'
      PLU4BL = 'LU4PAR'
      PLU5BL = 'LU5PAR'
      PLU7BL = 'LU7PAR'
      PLUCBL = 'LUCPAR'
       
C     set correct names
      IF( NFILE_ID .lt. 10 ) THEN    ! MPI ID has one digit
        WRITE (PDIA,'(A6,A1,I1)') PDIABL,'.',NFILE_ID
        WRITE (PLU1,'(A6,A1,I1)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I1)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I1)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I1)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I1)') PLU5BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I1)') PLU7BL,'.',NFILE_ID
        WRITE (PLUC,'(A6,A1,I1)') PLUCBL,'.',NFILE_ID
        PFILELENGTH = 8
      ELSE IF( NFILE_ID .lt. 100)THEN  ! MPI ID has two digits
        WRITE (PDIA,'(A6,A1,I2)') PDIABL,'.',NFILE_ID
        WRITE (PLU1,'(A6,A1,I2)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I2)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I2)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I2)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I2)') PLU5BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I2)') PLU7BL,'.',NFILE_ID
        WRITE (PLUC,'(A6,A1,I2)') PLUCBL,'.',NFILE_ID
        PFILELENGTH = 9
      ELSE                        ! MPI ID has three digits
        WRITE (PDIA,'(A6,A1,I3)') PDIABL,'.',NFILE_ID
        WRITE (PLU1,'(A6,A1,I3)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I3)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I3)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I3)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I3)') PLU5BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I3)') PLU7BL,'.',NFILE_ID
        WRITE (PLUC,'(A6,A1,I3)') PLUCBL,'.',NFILE_ID
        PFILELENGTH = 10
      END IF
C
C     ... open
C
C     if MCSCF --> delete on close! save diskspace AND
C     may be necessary --> DIAPAR.1 trick!
C
C     file info object - provide useful hints for the MPI implementation
C
      call interface_mpi_INFO_CREATE(FILE_INFO_OBJ)
C
C     ... number of CPUs sharing the following MPI-I/O files
      WRITE (FILE_INFO_GROUPSZ,'(I4)') NEWCOMM_PROC
      call interface_mpi_INFO_SET(FILE_INFO_OBJ,"nb_proc",
     &                            FILE_INFO_GROUPSZ)
C
#if defined (VAR_PFS)
C
C     special information on IBMs GPFS to enhance I/O performance
C
      call interface_mpi_INFO_SET(FILE_INFO_OBJ, "IBM_largeblock_io", 
     &                            "true")
#endif
C
      call interface_mpi_FILE_OPEN(MYNEW_COMM,PDIA(1:PFILELENGTH),
     &     df_MPI_MODE_CREATE + df_MPI_MODE_RDWR,FILE_INFO_OBJ,IDIA)
      call interface_mpi_FILE_OPEN(MYNEW_COMM,PLU1(1:PFILELENGTH),
     &     df_MPI_MODE_CREATE + df_MPI_MODE_RDWR 
     &   + df_MPI_MODE_DELETE_ON_CLOSE,
     &     FILE_INFO_OBJ,ILU1)
      call interface_mpi_FILE_OPEN(MYNEW_COMM,PLU2(1:PFILELENGTH),
     &     df_MPI_MODE_CREATE + df_MPI_MODE_RDWR
     &   + df_MPI_MODE_DELETE_ON_CLOSE,
     &     FILE_INFO_OBJ,ILU2)
      call interface_mpi_FILE_OPEN(MYNEW_COMM,PLU3(1:PFILELENGTH),
     &     df_MPI_MODE_CREATE + df_MPI_MODE_RDWR
     &   + df_MPI_MODE_DELETE_ON_CLOSE,
     &     FILE_INFO_OBJ,ILU3)
      call interface_mpi_FILE_OPEN(MYNEW_COMM,PLU4(1:PFILELENGTH),
     &     df_MPI_MODE_CREATE + df_MPI_MODE_RDWR
     &   + df_MPI_MODE_DELETE_ON_CLOSE,
     &     FILE_INFO_OBJ,ILU4)
      call interface_mpi_FILE_OPEN(MYNEW_COMM,PLU5(1:PFILELENGTH),
     &     df_MPI_MODE_CREATE + df_MPI_MODE_RDWR
     &   + df_MPI_MODE_DELETE_ON_CLOSE,
     &     FILE_INFO_OBJ,ILU5)
      call interface_mpi_FILE_OPEN(MYNEW_COMM,PLU7(1:PFILELENGTH),
     &     df_MPI_MODE_CREATE + df_MPI_MODE_RDWR
     &   + df_MPI_MODE_DELETE_ON_CLOSE,
     &     FILE_INFO_OBJ,ILU7)
      call interface_mpi_FILE_OPEN(MYNEW_COMM,PLUC(1:PFILELENGTH),
     &     df_MPI_MODE_CREATE + df_MPI_MODE_RDWR
     &   + df_MPI_MODE_DELETE_ON_CLOSE,
     &     FILE_INFO_OBJ,ILUC)
C
      call interface_mpi_INFO_FREE(FILE_INFO_OBJ)
C
C     ... set fileview
C
      call interface_mpi_FILE_SET_VIEW(IDIA,IDISP,
     &         df_MPI_REAL8,df_MPI_REAL8,
     &        "native",df_MPI_INFO_NULL)
      call interface_mpi_FILE_SET_VIEW(ILU1,IDISP,
     &         df_MPI_REAL8,df_MPI_REAL8,
     &        "native",df_MPI_INFO_NULL)
      call interface_mpi_FILE_SET_VIEW(ILU2,IDISP,
     &         df_MPI_REAL8,df_MPI_REAL8,
     &        "native",df_MPI_INFO_NULL)
      call interface_mpi_FILE_SET_VIEW(ILU3,IDISP,
     &         df_MPI_REAL8,df_MPI_REAL8,
     &        "native",df_MPI_INFO_NULL)
      call interface_mpi_FILE_SET_VIEW(ILU4,IDISP,
     &         df_MPI_REAL8,df_MPI_REAL8,
     &        "native",df_MPI_INFO_NULL)
      call interface_mpi_FILE_SET_VIEW(ILU5,IDISP,
     &         df_MPI_REAL8,df_MPI_REAL8,
     &        "native",df_MPI_INFO_NULL)
      call interface_mpi_FILE_SET_VIEW(ILU7,IDISP,
     &         df_MPI_REAL8,df_MPI_REAL8,
     &        "native",df_MPI_INFO_NULL)
      call interface_mpi_FILE_SET_VIEW(ILUC,IDISP,
     &         df_MPI_REAL8,df_MPI_REAL8,
     &        "native",df_MPI_INFO_NULL)
C
C     IF (NFLGRPS_REL == 1 .and. ICIONLY .eq. 1) THEN
C       check for existing DIAPAR.1 which could be reused...
C     END IF
C
      USE_EX_IDIA = .FALSE.
      END
***********************************************************************

      SUBROUTINE PART_CIV_PAR3(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,ISCLFAC,ISTART_BLK,
     &                         IEXCLBLK,IDEBUGPRNT)
C
C     Partition a CI vector into batches of blocks. The length of a
C     batch must be atmost MXLNG.
C     C-vector routine.
C
C     IF ICOMP. eq. 1: the complete CI vector is constructed in just one
C     batch.
C
C     OUTPUT
C     ======
C
C     NBATCH      : Number of batches
C     LBATCH(*)   : Number of blocks in a given batch
C     LEBATCH(*)  : Number of elements in a given batch ( packed ) !
C     I1BATCH(*)  : Number of first block in a given batch
C     IBATCH(8,*) : TTS blocks in Start of a given TTS block with respect to
C                   start
C     of batch --
C     IBATCH(1,*) : Alpha type
C     IBATCH(2,*) : Beta  type
C     IBATCH(3,*) : Sym of alpha
C     IBATCH(4,*) : Sym of beta
C     IBATCH(5,*) : Offset of block with respect to start of block in
C                   expanded form
C     IBATCH(6,*) : Offset of block with respect to start of block in
C                   packed form
C     IBATCH(7,*) : Length of block, expandend form
C     IBATCH(8,*) : Length of block, packed form
C    
C     original version : Jeppe Olsen     - August 1995
C     parallel adaption: S. Knecht       - May    2007 
C
C     Last revision:     S. Knecht       - May    2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
C     input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION ISCLFAC(*)
C     output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
C     scratch
      INTEGER LBLOCK_SAVE, LBLOCKP_SAVE 
C
      I_SAVE_IBLOCK = 0
      N_TRY = 0
C
      LBLOCK_SAVE  = 0
      LBLOCKP_SAVE = 0 
C
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ==================='
        WRITE(LUWRT,*) '    PART_CIV_PAR3   '
        WRITE(LUWRT,*) ' ==================='
        WRITE(LUWRT,*) ' IDC = ', IDC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' IOCOC Array '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        if (NTEST.ge.500) then
          WRITE(LUWRT,*) ' NSSOA array ( input ) '
          CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
          WRITE(LUWRT,*) ' NSSOB array ( input ) '
          CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
          write(LUWRT,*) ' IBLTP array: '
          call iwrtma(IBLTP,1,NSMST,1,NSMST)
        end if
      END IF
C
C     block zero
C
      ILOOPBLK = 0
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = 0
      IBLOCK = 0
      IFINI = 0
C     loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
      IF( NBATCH .gt. 1 ) ILOOPBLK = ILOOPBLK - 1 
C     loop over blocks in batch
 1000 CONTINUE
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF (IFINI.EQ.1) GOTO 2002
C     should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
      ILOOPBLK = ILOOPBLK + 1
C     can this block be included
      IBSM   = ISMOST(ISM)
      NSTA   = NSSOA(ISM,IA)
      NSTB   = NSSOB(IBSM,IB)
      LBLOCK = NSTA*NSTB
      LBLOCK_SAVE = LBLOCK
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
      LBLOCKP_SAVE = LBLOCKP
C     check if we need this block, if not: 
C     LBLOCKP = 0
C     LBLOCK  = 0
CSK      IF(IEXCLBLK .eq. 1) THEN
CSK       WRITE(6,*) 'ILOOPBLK,SCLFAC(ILOOPBLK),IBLOCK',ILOOPBLK,
CSK     &             SCLFAC(ILOOPBLK),IBLOCK
CSK      END IF
*
*     IEXCLBLK == 1: DEBUG ME!!!
*
      IF( IEXCLBLK .eq. 1 ) THEN
         IF( ISCLFAC(ILOOPBLK) .eq. 0 ) THEN
           LBLOCKP = 0
           LBLOCK  = 0
           GOTO 1000
         END IF
      ELSE
        IF( ISCLFAC(ISTART_BLK + IBLOCK + 1 ) .eq. 0 ) THEN
          LBLOCKP = 0
          LBLOCK  = 0
        END IF
      END IF
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK,myproc
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
C       keep length information on 7!
        IBATCH(7,IBLOCK) = LBLOCK_SAVE
        IBATCH(8,IBLOCK) = LBLOCKP
C       all blocks are included, but only active blocks have a length
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(LUWRT,*) ' Not enough scratch space to include a 
     &                   single block'
        WRITE(LUWRT,*) ' Since I cannot procede I will stop '
        WRITE(LUWRT,*) ' Insufficient buffer detected in PART_CIV_PAR3'
        write(LUWRT,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(LUWRT,*) ' Alter GAS space of raise buffer from ', MXLNG
        call quit( ' Insufficient buffer space in PART_CIV_PAR3. ' )
      ELSE
*
*     this batch is finished, goto next batch
*
*       is this the TT-1 block?
        IF( ISM .ne. 1 )THEN
*
*         go back to ISM == 1 since we loop over NSMST in GNSIDE_REL
*
          N_TRY = N_TRY + 1
*
          IF( N_TRY .gt. 1 .and. I_SAVE_IBLOCK .eq. IBLOCK ) THEN
            WRITE(LUWRT,*) ' Not enough scratch space to include a '//
     &                     ' single Block'
            WRITE(LUWRT,*) ' Since I cannot procede I will stop '
            WRITE(LUWRT,*) ' Insufficient buffer detected in '//
     &                     ' PART_CIV_PAR3'
            WRITE(LUWRT,*) '  LENGTH ',LENGTH
            WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from',MXLNG
            stop ' In PART_CIV_PAR3 because of N_TRY .ge. 1. '
          END IF
*
          I_SAVE_IBLOCK = IBLOCK
*         set batch back to first available TT-1 block
          call bck_tts(ITTSS_ORD,ISM,LBATCH,LEBATCH,IBATCH,IBLOCK,
     &                 NBLOCK,LENGTH,LENGTHP,NBATCH,LUWRT)
*
          GOTO 2000
*
        ELSE
          N_TRY = 0
          I_SAVE_IBLOCK = 0
          GOTO 2000
        END IF
*       ^ ISM check
      END IF
*     ^ length of batch < MXLNG?
 2002 CONTINUE
C
      IF( IDEBUGPRNT .ne. 0 ) THEN
        WRITE(LUWRT,*) 'Output from PART_CIV_PAR3'
        WRITE(LUWRT,*) '========================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of batches (MYPROC)', NBATCH, MYPROC
        DO JBATCH = 1, NBATCH
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Info on batch ', JBATCH
          WRITE(LUWRT,*) ' *********************** '
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '      Length of batch           ', 
     &                          LEBATCH(JBATCH)
          WRITE(LUWRT,*) '      Number of blocks included ', 
     &                          LBATCH(JBATCH)
          WRITE(LUWRT,*) '      TTSS and offsets and lengths of 
     &                          each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(LUWRT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
C
      END
***********************************************************************

      SUBROUTINE PART_CIV4_PAR(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         IOFFBLK,IOFFBTC,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,IBLOCKD,IDEBUGPRNT,
     &                         ILENB_OLD,INUMB_OLD,IRUN_THIS)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
* Stefan Knecht, May 2007
*
* further modified from PART_CIV3 for parallel purposes
*
* Compared to PART_CIV3 the changes are
*
* 1 : partitioning with respect to block <--> node distribution
*
* 2 : to be continued ...
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "parluci.h"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION IBLOCKD(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ====================='
        WRITE(LUWRT,*) '     PART_CIV4_PAR    '
        WRITE(LUWRT,*) ' ====================='
        WRITE(LUWRT,*) ' IDC = ', IDC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' NSSOA array ( input ) '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        WRITE(LUWRT,*) ' NSSOB array ( input ) '
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
*
*. block  zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      IFINI = 0
      IF( IRUN_THIS .eq. 1 ) THEN
        NBATCH = IOFFBTC - 1
CSK        WRITE(LUWRT,*) 'NBATCH is ',NBATCH
        IBLOCK = IOFFBLK - 1
      ELSE
        NBATCH = IOFFBTC 
C       be careful - is that correct???
        NBATCH = IOFFBTC
CSK        WRITE(LUWRT,*) 'NBATCH is now',NBATCH
C
        IBLOCK = IOFFBLK - 1
        LBATCH(  NBATCH ) = INUMB_OLD
        LEBATCH( NBATCH ) = ILENB_OLD
        LENGTH  = ILENB_OLD  
        LENGTHP = ILENB_OLD
        GOTO 1000
      END IF
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, ISM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(ISM)
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
C     set unpacked length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCK = 0
C
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C     set packed length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCKP = 0
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(LUWRT,*) ' Not enough scratch space to include a 
     &                  single Block'
        WRITE(LUWRT,*) ' Since I cannot procede I will stop '
        WRITE(LUWRT,*) ' Insufficient buffer detected in PART_CIV4_PAR'
        WRITE(LUWRT,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from ', MXLNG
        stop ' In PART_CIV4_PAR. '
      ELSE
*. This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
*
      ILENB_OLD = LEBATCH(NBATCH)
      INUMB_OLD = LBATCH(NBATCH)
*
      IF( IDEBUGPRNT .ne. 0 ) THEN
        WRITE(LUWRT,*) 'Output from PART_CIV4_PAR'
        WRITE(LUWRT,*) '========================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Info on batch ', JBATCH
          WRITE(LUWRT,*) ' *********************** '
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '      Number of blocks included ', 
     &                          LBATCH(JBATCH)
          WRITE(LUWRT,*) '      TTSS and offsets and lengths of 
     &                          each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(LUWRT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      END
***********************************************************************

      SUBROUTINE FIND_IMAT_SC_REL(ISCALFAC_C,ISCALFAC_S,ICCTOS,
     &                            ICWEIGHTF,IAM_BLK_ACT_C,
     &                            IAM_BLK_ACT_S,NDIM)
C**********************************************************************
C
C     Compute connection matrix ICCTOS between the sigma and c vector
C     Store contribution of one c-block to sigma-blocks in ICWEIGHTF.
C
C     NOTE: only active blocks are counted in the connection matrix!
C
C     Written by  S. Knecht         - June 20 2007 
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
#include "krciprop.h"
      DIMENSION ISCALFAC_C(*),ISCALFAC_S(*),ICCTOS(I_NZERO_LEN_C,*)
      DIMENSION ICWEIGHTF(*), IAM_BLK_ACT_S(*), IAM_BLK_ACT_C(*)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "ipoist8.inc"
#include "cstate.inc"
#include "cands.inc"
C
C
      IDUM = 0
      CALL MEMMAR(KXDUM,  IDUM,    'MARK  ',IDUM,'FINDSC')
C
C     arrays for partitioning of sigma
C     (*_MS2 arrays are in CSTATE )
      JCMBSPC = 1
C     we use I_SET_L2BLOCK as ICSM carrier to
C     z_blkfo_rel. note the minus sign!
      ICSM_SCR = 0
      IF( DOSIGPROP ) ICSM_SCR =  - ICSM
      CALL Z_BLKFO_REL(IDC,NMS2VAL,JCMBSPC,ISSM,
     &                 KLSLBT,KLSLEBT,KLSI1BT,KLSIBT,KDUM,
     &                 NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,
     &                 NBAT_MS2,IBAT_MS2,NBLK_MS2_C,DOSIGPROP,
     &                 ICSM_SCR,3,IDUMMY)
C     SUBROUTINE Z_BLKFO_REL(IDC,NMS2VAL,ISPC,ISM,KPCLBT,KPCLEBT,
C    &                       KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK,
C    &                       NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2,
C    &                       NBLK_MS2_C,CALC_MS2_C,I_SET_L2BLOCK,
C    &                       I_USE_PC,NPARBLOCK)

C
      ISIGDEN = 1
C
      CALL FIND_IMAT_DRV1(NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                    WORK(KLSI1BT),WORK(KLSIBT),ISIGDEN,
     &                    ISCALFAC_C,ISCALFAC_S,ICCTOS,ICWEIGHTF,NDIM,
     &                    IAM_BLK_ACT_C,IAM_BLK_ACT_S)
*
      CALL MEMMAR(KXDUM,IDUM,'FLUSM ',IDUM,'FINDSC')
*
      END
***********************************************************************

      SUBROUTINE FIND_IMAT_DRV1(NBATS,LBATS,LEBATS,I1BATS,IBATS,ISIGDEN,
     &                          ISCLFAC_C,ISCLFAC_S,ICCTOS,ICWEIGHTF,
     &                          NDIM_B,IAM_BLK_ACT_C,IAM_BLK_ACT_S)
*
      use symmetry_setup_krci
      use mospinor_info
      use interface_to_mpi
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
* =====
*.Input
* =====
#include "cands.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "cprnt.inc"
#include "ctcc.inc"
#include "clunit.inc"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE), ICWEIGHTF(*), ISCLFAC_C(*)
      INTEGER   IAM_BLK_ACT_S(*),IAM_BLK_ACT_C(*),ISCLFAC_S(*)
      INTEGER   ICCTOS(I_NZERO_LEN_C,I_NZERO_LEN_S) 
#include "parluci.h"
#include "krciprop.h"
#include "dgroup.h"
C     batches of sigma
      DIMENSION LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
C     scratch
      INTEGER IONE, INT_IOFF1, INT_IOFF2, INT_IOFF1_C, INT_IOFF2_C
      LOGICAL CHECK_TSPLIT
C
      IONE = 1
      INT_IOFF1 = 0
      INT_IOFF2 = 0
C
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'FNDSC1')
C
*     assumption: no matter if NZ ==1 or NZ == 2, 
*     the contribution of a block is CONTRIB * NZ
*     ==> IRILP = 1 even for NZ == 2
      IRILP = 1
      IF( CHECK_TC ) IRILP = IRC_SAVE
*
*     NOTE: IF( SPLIT_IJKL .and. NDIM_B .eq. - 1 ) T-connection run 
*
      CHECK_TSPLIT = .FALSE.
      IF ( SPLIT_IJKL .and. NDIM_B .eq. -1 ) THEN 
           CHECK_TSPLIT = .TRUE.
           IRILP = IRC_SAVE
      END IF
*
*     loop over vectors of sigma and C : MK2 space
      DO IRIS = 1,IRILP,1
*
        INT_IOFF1  = 0
        INT_IOFF2  = 0
*
        INT_IOFF1  = 1 +  I_NZERO_LEN_S * ( IRIS - 1 )
*
*       FIXME: find max allocation down in sigden_rel2
        IVAL_ALLOC_SCR = 0
*
*       count the sigma blocks
        NUM_SBLK = 0
*
        JBATABS = 0
        DO IMK2_S = 1, NMS2VAL
          MK2_S = MS2VAL(IMK2_S)
*
          IASTP = IST_FOR_DT(1,IMK2_S)
          IBSTP = IST_FOR_DT(2,IMK2_S)
          NOCTPA_S = NOCTYP(IASTP)
          NOCTPB_S = NOCTYP(IBSTP)
*         arrays giving allowed type combinations (of alpha and beta)
          call memmar(KSIOIO,NOCTPA_S*NOCTPB_S,'ADDL  ',2,'SIOIO ')
          CALL IAIBCM_REL(ISSPC,IASTP,IBSTP,WORK(KSIOIO))
*         arrays giving block type
          KSVST = 1
          call memmar(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
          CALL ZBLTP_REL(ISMOST(1,ISSM),NSMST,WORK(KSBLTP))
*
*         loop over batches of sigma for this MK2_S projection value
*
          do JBATS = 1,NBAT_MS2(IMK2_S)
*
            JBATABS = JBATABS + 1
*
*           check whether we need to call FIND_IMAT_DRV2
*
            KBATSEND = I1BATS(JBATABS) + LBATS(JBATABS)
*
            ICOMPUTE = 0
*
            LS = 0
*
            IF( CHECK_TSPLIT ) THEN
              DO ISBLK = I1BATS(JBATABS),KBATSEND-1,1
                IF( IAM_BLK_ACT_S(ISBLK) .eq. MYPROC ) ICOMPUTE = 1
              END DO
            ELSE 
              DO ISBLK = I1BATS(JBATABS),KBATSEND-1,1
!               print *,'icompute check: ISBLK',ISBLK
                IF( ISCLFAC_S(ISBLK) .eq. 1 ) ICOMPUTE = 1
!               print *,'icompute check: isblk, IS.. ==> ',
!    &                   ISBLK,ISCLFAC_S(ISBLK)
              END DO
            END IF
*
!           WRITE(LUWRT,*) ' ICOMPUTE is ', ICOMPUTE
*
            IF( ICOMPUTE .eq. 0 ) THEN
              GOTO 60
            ENDIF
*
            DO IRIC = 1,IRILP,1
*
              INT_IOFF1_C = 0
              INT_IOFF2_C = 0
*
              INT_IOFF1_C = 1 + I_NZERO_LEN_C * ( IRIC -1 )
              NUM_BLK = 0
*
              DO IMK2_C = 1, NMS2VAL
                MK2_C = MS2VAL(IMK2_C)
C
                MK2_DIFF = IABS(MK2_C - MK2_S)
                max_MK2_DIFF = 4
                IF( NTEST .ge. 50 )THEN
                  IF (MK2_DIFF .gt. max_MK2_DIFF) THEN
                    write(LUWRT,*) 'MK2 skip, MK2_DIFF =',MK2_DIFF
                  END IF
                END IF
*               ... for any two-electron operator.
*                   For one-electron operators max_MK2_DIFF = 1 or 0
*                   (or zero, if NZ .lt. 4 and operator is not j- or
*                   k-imaginary) /HJAAJ Aug 2008. TODO TODO to define
*                   max_MK2_DIFF
*
                IACTP = IST_FOR_DT(1,IMK2_C)
                IBCTP = IST_FOR_DT(2,IMK2_C)
                NBL_C = NBLK_MS2(IMK2_C)
C               C and sigma may have different symmetry in general
                IF( DOSIGPROP ) NBL_C = NBLK_MS2_C(IMK2_C)
*               set new offset
                INT_IOFF1_C = INT_IOFF1_C + INT_IOFF2_C
*
!               print *,'NUM_SBLK is ...',NUM_SBLK
                IF (MK2_DIFF .le. max_MK2_DIFF)
     &          CALL FIND_IMAT_DRV2(LBATS(JBATABS),
     &                              IBATS(1,I1BATS(JBATABS)),1,
     &                              1,IASTP,IBSTP,IRIS,IRIC,
     &                              INT_IOFF1_C,ISCLFAC_C,ISCLFAC_S,
     &                              ICCTOS,ICWEIGHTF,INT_IOFF1,
     &                              NUM_SBLK,IAM_BLK_ACT_C,
     &                              IAM_BLK_ACT_S,CHECK_TSPLIT)
C
C               keep track of correct offset
                NUM_BLK = NUM_BLK + NBL_C
                INT_IOFF2_C = NBL_C
*
              END DO
*             ^ End of loop over IMK2_C
            END DO
*           ^ End of loop over IRIC
*
 60         CONTINUE
*           count the sigma blocks
*
            NUM_SBLK = NUM_SBLK + LBATS(JBATABS)
*
          end do
*         ^ End loop over S batches of given MK2_S value
        END DO
*       ^ End of loop over IMK2_S
      END DO
*     ^ End of loop over IRIS
      call memmar(KDUM ,IDUM,'FLUSM ',2,'FNDSC1')
*
      END
***********************************************************************

      SUBROUTINE FIND_IMAT_DRV2(NBLOCK,IBLOCK,IBOFF,ISIGDEN,IASTP_SCR,
     &                         IBSTP_SCR,IRIS,IRIC,
     &                         INT_IOFF1_C,ISCLFAC_C,ISCLFAC_S,ICCTOS,
     &                         ICWEIGHTF,INT_IOFF1,NUM_SBLK,
     &                         IAM_BLK_ACT_C,IAM_BLK_ACT_S,CHECK_TSPLIT)
*
      use symmetry_setup_krci
      use mospinor_info
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
#include "cands.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "ctcc.inc"
      COMMON/CMXCJ/MXCJ
*
#include "cintfo.inc"
#include "parluci.h"
      LOGICAL CHECK_TSPLIT
*
!               print *,'NUM_SBLK - 2 - is ...',NUM_SBLK
      IASTP = IASTP_SCR
      IBSTP = IBSTP_SCR
*
      KDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'FNDSC2')
*
      L0BLOCK = MXSOOB_AS
      LSCR1 = MAX(L0BLOCK,L2BLOCK)
*. Offsets for alpha and beta supergroups
*
      NOCTPA_C = NOCTYP(IACTP)
      NOCTPB_C = NOCTYP(IBCTP)
*
      NAEL_C = NELEC(IACTP)
      NBEL_C = NELEC(IBCTP)
*
*. Arrays giving allowed type combinations
      call memmar(KCIOIO,NOCTPA_C*NOCTPB_C,'ADDL  ',1,'CIOIO ')
      CALL IAIBCM_REL(ICSPC,IACTP,IBCTP,WORK(KCIOIO))
*
*. Arrays for storing NEL consecutive annihilations/creations
*  (resolution matrices)
* Number of alpha/beta excitation types
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &                      1,NAOBEX_TP,NBOBEX_TP,IDUMMY,IDUMMY)
*. And the alpha-and beta-excitations
      LENA = 2*NGAS*NAOBEX_TP
      LENB = 2*NGAS*NBOBEX_TP
      call memmar(KLAOBEX,LENA,'ADDL  ',1,'IAOBEX')
      call memmar(KLBOBEX,LENB,'ADDL  ',1,'IAOBEX')
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &                      0,NAOBEX_TP,NBOBEX_TP,
     &                      WORK(KLAOBEX),WORK(KLBOBEX))
*
*. For alpha excitations, C and S vector (might be flipped around!)
      IOCTPA_C = IBSPGPFTP(IACTP)
      CALL LEN_GENOP_STR_MAP(NAOBEX_TP,WORK(KLAOBEX),NOCTPA_C,
     &                       NELFSPGP(1,IOCTPA_C),NOBPT,NGAS,
     &                       MAXLENAC)
      IOCTPA_S = IBSPGPFTP(IASTP)
      NOCTPA_S = NOCTYP(IASTP)
      CALL LEN_GENOP_STR_MAP(NAOBEX_TP,WORK(KLAOBEX),NOCTPA_S,
     &                       NELFSPGP(1,IOCTPA_S),NOBPT,NGAS,
     &                       MAXLENAS)
*
*. For beta excitations
      IOCTPB_C = IBSPGPFTP(IBCTP)
      CALL LEN_GENOP_STR_MAP(NBOBEX_TP,WORK(KLBOBEX),NOCTPB_C,
     &                       NELFSPGP(1,IOCTPB_C),NOBPT,NGAS,
     &                       MAXLENBC)
      IOCTPB_S = IBSPGPFTP(IBSTP)
      NOCTPB_S = NOCTYP(IBSTP)
      CALL LEN_GENOP_STR_MAP(NBOBEX_TP,WORK(KLBOBEX),NOCTPB_S,
     &                       NELFSPGP(1,IOCTPB_S),NOBPT,NGAS,
     &                       MAXLENBS)
*
      MAXLEN_I1 = MAX(MAXLENAC,MAXLENBC,MAXLENAS,MAXLENBS)
C?    WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1
*
      LSCR3 = MAXLEN_I1
C?    print*,'MAXLENAC,MAXLENBC,MAXLENAS,MAXLENBS ',
C?   &        MAXLENAC,MAXLENBC,MAXLENAS,MAXLENBS
*
*. Arrays for storing occupations of T-operators
      call memmar(KTOCC1,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'TOCC1 ')
      call memmar(KTOCC2,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'TOCC2 ')
      call memmar(KTOCC3,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'TOCC3 ')
      call memmar(KTOCC4,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'TOCC4 ')
*
*.TTS arrays for partitioning of vector 
      NTTS = MXNTTS
      call memmar(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      call memmar(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      call memmar(KLI1BT,2*NTTS,'ADDL  ',1,'I1BTC ')
      call memmar(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
*. Arrays giving block type
      call memmar(KCBLTP,NSMST,'ADDL  ',1,'CBLTP ')
      CALL ZBLTP_REL(ISMOST(1,ICSM),NSMST,WORK(KCBLTP))
*. Use all symmetryblocks of given TT
      ITTSS_ORD = 2
*
      CALL PART_CIV2_SPC(IDC,WORK(KCBLTP),WORK(KNSTSO(IACTP)),
     &     WORK(KNSTSO2(IBCTP)),NOCTPA_C,NOCTPB_C,NSMST,LSCR1,
     &     WORK(KCIOIO),ISMOST(1,ICSM),NCBATCH,WORK(KLLBT),
     &     WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),0,ITTSS_ORD)
 
!     allocation for string length handling of unbarred (alpha)/ barred
!     (beta) creation / annihilation strings
      call memmar(KNICA    ,MXNDGIRR              ,'ADDL  ',1,'KNICA ')
      call memmar(KNICB    ,MXNDGIRR              ,'ADDL  ',1,'KNICB ')
      call memmar(KNIAA    ,MXNDGIRR              ,'ADDL  ',1,'KNIAA ')
      call memmar(KNIAB    ,MXNDGIRR              ,'ADDL  ',1,'KNIAB ')
      call memmar(KNKA     ,MXNDGIRR              ,'ADDL  ',1,'KNKA  ')
      call memmar(KNKB     ,MXNDGIRR              ,'ADDL  ',1,'KNKB  ')
      
      call memmar(KIB_CA   ,MXNDGIRR**2 * MXNOP   ,'ADDL  ',1,'KIB_CA')
      call memmar(KIB_CB   ,MXNDGIRR**2 * MXNOP   ,'ADDL  ',1,'KIB_CB')
      call memmar(KIB_AA   ,MXNDGIRR**2 * MXNOP   ,'ADDL  ',1,'KIB_AA')
      call memmar(KIB_AB   ,MXNDGIRR**2 * MXNOP   ,'ADDL  ',1,'KIB_AB')
      call memmar(KIB_T    ,MXNDGIRR**2 * MXNDGIRR,'ADDL  ',1,'KIB_T ')
      call memmar(KNIKAINTM,MXPOBS      * MXNOP   ,'ADDL  ',1,'KNIKA ')
      call memmar(KNIKBINTM,MXPOBS      * MXNOP   ,'ADDL  ',1,'KNIKB ')
      call memmar(KNJKAINTM,MXPOBS      * MXNOP   ,'ADDL  ',1,'KNJKA ')
      call memmar(KNJKBINTM,MXPOBS      * MXNOP   ,'ADDL  ',1,'KNJKB ')

      CALL MEMMAR(KLISMDST,MAX_NSYMDIST,'ADDL  ',1,'ISMDST')
*
!               print *,'NUM_SBLK - 3 - is ...',NUM_SBLK
      CALL FIND_IMAT_DRV3(NBLOCK,IBLOCK(1,IBOFF),
     &     WORK(KNSTSO(IASTP)),WORK(KNSTSO2(IBSTP)),
     &     WORK(KNSTSO(IACTP)),WORK(KNSTSO2(IBCTP)),
     &     IASTP,IBSTP,IACTP,IBCTP,MXINKA,
     &     WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &     ISSM,ICSM,ISIGDEN,NCBATCH,IDIAG,IRIS,IRIC,
     &     NSPOBEX_TP,WORK(KLSOBEX),
     &     WORK(KTOCC1),WORK(KTOCC2),WORK(KTOCC3),WORK(KTOCC4),
     &     WORK(KNICA),WORK(KNICB),WORK(KNIAA),WORK(KNIAB),
     &     WORK(KNKA),WORK(KNKB),
     &     WORK(KIB_CA),WORK(KIB_CB),WORK(KIB_AA),
     &     WORK(KIB_AB),WORK(KIB_T),WORK(KNIKAINTM),
     &     WORK(KNIKBINTM),WORK(KNJKAINTM),WORK(KNJKBINTM),
     &     WORK(KLISMDST),INT_IOFF1_C,ISCLFAC_C,ISCLFAC_S,ICCTOS,
     &     ICWEIGHTF,INT_IOFF1,NUM_SBLK,IAM_BLK_ACT_C,IAM_BLK_ACT_S,
     &     CHECK_TSPLIT)
C
C     release local memory
      IDUM = 0
      call memmar(KDUM ,IDUM,'FLUSM ',2,'FNDSC2')
      END
***********************************************************************

      SUBROUTINE FIND_IMAT_DRV3(NSBLOCK,ISBLOCK,
     &                          NSSOA_S,NSSOB_S,NSSOA_C,NSSOB_C,
     &                          IASTP,IBSTP,IACTP,IBCTP,MAXK,
     &                          LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,
     &                          ISSM,ICSM,ISIGDEN,NCBATCH,IDIAG,IRIS,
     &                          IRIC,NSPOBEX_TP,ITSPOBEX_TP,
     &                          TOCC1,TOCC2,TOCC3,TOCC4,
     &                          NICA,NICB,NIAA,NIAB,NKA,NKB,
     &                          IB_CA,IB_CB,IB_AA,IB_AB,IB_T,
     &                          NIKAINTM,NIKBINTM,NJKAINTM,NJKBINTM,
     &                          ISMDST,
     &                          INT_IOFF1_C,ISCLFAC_C,ISCLFAC_S,ICCTOS,
     &                          ICWEIGHTF,INT_IOFF1,NUM_SBLK,
     &                          IAM_BLK_ACT_C,IAM_BLK_ACT_S,
     &                          CHECK_TSPLIT)
*
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included 
*
* NSSOA : Number of strings per type and symmetry for alpha strings
* NSSOB : Number of strings per type and symmetry for beta strings
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
      use symmetry_setup_krci
      use mospinor_info
      use interface_to_mpi
#include "implicit.h"
#include "ipoist8.inc"
*
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE), ICWEIGHTF(*),ISCLFAC_C(*)
      INTEGER   ICCTOS(I_NZERO_LEN_C,I_NZERO_LEN_S), ISCLFAC_S(*)
      INTEGER   IAM_BLK_ACT_S(*), IAM_BLK_ACT_C(*)
#include "parluci.h"
*. General input
#include "mxpdim.inc"
#include "gasstr.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "cstate.inc"
#include "strbas.inc"
*. Specific input
      INTEGER ISBLOCK(8,*)
      INTEGER LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*.General input
      INTEGER NSSOA_S(NSMST,*),NSSOB_S(NSMST,*)
      INTEGER NSSOA_C(NSMST,*),NSSOB_C(NSMST,*)
*.Scratch
      DIMENSION ISMDST(*)
*. T-coefficients
      DIMENSION ITSPOBEX_TP(4*NGAS,*)
      LOGICAL CHECK_TSPLIT
*
      ONEM = -1.0D0
      IZERO = 0
      IST_C_BLK_CT = 0
*
      IOCTPA = IBSPGPFTP(IASTP)
      IOCTPB = IBSPGPFTP(IBSTP)
      JOCTPA = IBSPGPFTP(IACTP)
      JOCTPB = IBSPGPFTP(IBCTP)
*
      NSTT_BLK = NSBLOCK/NSMST
* Loop over batches over C blocks      
*
      DO 20000 JCBATCH = 1, NCBATCH             
*
*. Read C blocks into core
*
        NJBLOCK = LCBLOCK(JCBATCH)
        NCTT_BLK = NJBLOCK/NSMST
        I1C = I1CBLOCK(JCBATCH)
CSK        WRITE(LUWRT,*) ' JCBATCH, NJBLOCK, NCTT_BLK, I1C',
CSK     & JCBATCH, NJBLOCK, NCTT_BLK, I1C
*
*. Loop over TT blocks of sigma and C in batches and 
*  obtain contributions 
        DO 9000 ICTT_BLK = 1, NCTT_BLK
*. first block of next TT block of C    
          ICBLK = I1C + (ICTT_BLK-1)*NSMST
          JATP  = ICBLOCK(1,ICBLK)
          JBTP  = ICBLOCK(2,ICBLK)
          ICOFF = ICBLOCK(5,ICBLK)
*
          IF(NTEST.GE.500) THEN
            WRITE(LUWRT,*) ' Next block of C, ICBLK,JATP,JBTP,ICOFF ',
     &      ICBLK,JATP,JBTP
          END IF
*
*         check if TT-blocks are 'active'! if not, 
*         skip all loops that will follow.
*
          ICOMPUTE_C = 0
          NCBLK = NSMST
          IC_OFF_BLK = INT_IOFF1_C
          IST_C_BLK_CT = IC_OFF_BLK
          DO JBLK = ICBLK, ICBLK + NCBLK - 1
csk         WRITE(LUWRT,*) ' Checking for c block',IC_OFF_BLK+JBLK-1
            IF( ISCLFAC_C(IC_OFF_BLK+JBLK-1) .gt. 0 ) THEN
              ICOMPUTE_C = 1
csk            WRITE(LUWRT,*) ' active c block found',IC_OFF_BLK+JBLK-1
            END IF
          END DO
*
          IF( ICOMPUTE_C .eq. 0 ) GOTO 8999
*
          IRATP = JATP + JOCTPA - 1
          IRBTP = JBTP + JOCTPB - 1
*
          DO 10000 ISTT_BLK = 1, NSTT_BLK
*
*. first block of next TT block of Sigma
            ISBLK = (ISTT_BLK-1)*NSMST + 1
*
*           check if TT-blocks need to be calculated! if not,
*           skip all loops that follow.
*
            ICOMPUTE_S = 0
            NSBLK      = NSMST
            IS_OFF_BLK = 0
            IS_OFF_BLK = NUM_SBLK + 1
!           WRITE(LUWRT,*) 'starting at IS_OFF_BLK',IS_OFF_BLK
            IF( CHECK_TSPLIT ) THEN
              DO JJBLK = ISBLK, ISBLK + NSBLK - 1
csk               WRITE(LUWRT,*) 'Checking for sigma block',
csk  &                            IS_OFF_BLK+JJBLK-1
                  IF( IAM_BLK_ACT_S( JJBLK + IS_OFF_BLK-1 ).eq. 
     &                MYPROC ) THEN
                    ICOMPUTE_S = 1
csk                 WRITE(LUWRT,*) ' active s block found',
csk  &                               IS_OFF_BLK+JJBLK-1
                  END IF
              END DO
            ELSE 
              DO JJBLK = ISBLK, ISBLK + NSBLK - 1
!                 WRITE(LUWRT,*) 'Checking for sigma block',
!    &                            IS_OFF_BLK+JJBLK-1
                  IF( ISCLFAC_S( JJBLK + IS_OFF_BLK-1 ) .gt. 0 ) THEN
                    ICOMPUTE_S = 1
!                   WRITE(LUWRT,*) ' active s block found',
!    &                               IS_OFF_BLK+JJBLK-1
                  END IF
              END DO
            END IF
            IF( ICOMPUTE_S .eq. 0 ) GOTO 9999
*
            IATP = ISBLOCK(1,ISBLK)
            IBTP = ISBLOCK(2,ISBLK)
            ISOFF = ISBLOCK(5,ISBLK)
CSK            WRITE(LUWRT,*) ' Next block of S, ISBLK,IATP,IBTP,ISOFF ',
CSK     &      ISBLK,IATP,IBTP,ISOFF
*
            ILATP = IATP + IOCTPA - 1
            ILBTP = IBTP + IOCTPB - 1
*. Connections ?
            DO 8000 ITTP = 1, NSPOBEX_TP
*
              ICA_OFF = 1
              ICB_OFF = 1 +  NGAS
              IAA_OFF = 1 +2*NGAS
              IAB_OFF = 1 +3*NGAS
*. Connections ?
               CALL GXFSTR(NELFSPGP(1,ILATP),NELFSPGP(1,ILBTP),
     &                     NELFSPGP(1,IRATP),NELFSPGP(1,IRBTP),
     &                     ITSPOBEX_TP(ICA_OFF,ITTP), 
     &                     ITSPOBEX_TP(IAA_OFF,ITTP),
     &                     ITSPOBEX_TP(ICB_OFF,ITTP), 
     &                     ITSPOBEX_TP(IAB_OFF,ITTP),
     &                     NGAS,ICON)
*
               IF(ICON .eq. 1)THEN
*
CSK                 WRITE(LUWRT,*) 'CONNECTION FOUND for ITTP', ICON, ITTP
C
C                ... ready for connection matrix determination
C
!                print *,' ITTP is ==> ',ITTP
!                call RMEMCHK_opt('next type           ')
                 CALL FIND_IMAT_CALC(ISIGDEN,IRIS,IRIC,
     &                ITSPOBEX_TP(ICA_OFF,ITTP),
     &                ITSPOBEX_TP(ICB_OFF,ITTP),
     &                ITSPOBEX_TP(IAA_OFF,ITTP),
     &                ITSPOBEX_TP(IAB_OFF,ITTP),
     &                NELFSPGP(1,ILATP),NELFSPGP(1,ILBTP),
     &                NELFSPGP(1,IRATP),NELFSPGP(1,IRBTP),
     &                NSSOA_S(1,IATP),NSSOB_S(1,IBTP),
     &                NSSOA_C(1,JATP),NSSOB_C(1,JBTP),
     &                ITTP,ISSM,ICSM,
     &                MAXK,TOCC1,TOCC2,TOCC3,TOCC4,
     &                NICA,NICB,NIAA,NIAB,NKA,NKB,
     &                IB_CA,IB_CB,IB_AA,IB_AB,IB_T,
     &                NIKAINTM,NIKBINTM,NJKAINTM,NJKBINTM,
     &                ISMDST,INT_IOFF1_C+ICBLK-1,IS_OFF_BLK+ISBLK-1,
     &                ISCLFAC_C,ISCLFAC_S,ICCTOS,ICWEIGHTF,
     &                IAM_BLK_ACT_C,IAM_BLK_ACT_S,CHECK_TSPLIT)
               END IF
*              ^ End if connection 
 8000       CONTINUE
*.          ^  End of loop over TT sigma blocks 
 9999       CONTINUE
*           ^ skip TT-blocks for sigma if not included in nodelist
10000     CONTINUE
*         ^ End of loop over blocks of sigma 
*
 8999     CONTINUE
*         ^ skip TT-c-blocks if not 'active'
 9000   CONTINUE
*.      ^ End of loop over TT C blocks in Batch
20000 CONTINUE
*.    ^End of loop over batches of C blocks
*
C         
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_IMAT_CALC(ISD,IRIS,IRIC,ICA,ICB,IAA,IAB,
     &                          IAOC,IBOC,JAOC,JBOC,
     &                          NIA,NIB,NJA,NJB,
     &                          IXTP,ISSM,ICSM,
     &                          MAXLB,ICA_STR,ICB_STR,IAA_STR,IAB_STR,
     &                          NICA,NICB,NIAA,NIAB,NKA,NKB,
     &                          IB_CA,IB_CB,IB_AA,IB_AB,IB_T,
     &                          NIKAINTM,NIKBINTM,NJKAINTM,NJKBINTM,
     &                          ISMDST,IOFF_BLK_C,IOFF_BLK_S,IACT_C,
     &                          IACT_S,ICCTOS,ICWEIGHTF,IAM_BLK_ACT_C,
     &                          IAM_BLK_ACT_S,CHECK_TSPLIT)
*
* IX_RESTRICT = 1: BLocks in X are restricted to those with IASM.GE.IBSM
* IX_RESTRICT =-1: BLocks in X are restricted to those with IASM.LT.IBSM
*
      use mospinor_info
      use symmetry_setup_krci
#include "implicit.h"
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "ctcc.inc"
#include "parluci.h"
*.Input
      DIMENSION NIA(*),NIB(*),NJA(*),NJB(*),IAOC(*),IBOC(*)
      DIMENSION ICA(*),ICB(*),IAA(*),IAB(*),JAOC(*),JBOC(*)
      DIMENSION IACT_C(*), IACT_S(*),IAM_BLK_ACT_C(*),IAM_BLK_ACT_S(*)
      DIMENSION ICCTOS(I_NZERO_LEN_C,*),ICWEIGHTF(*)
*.Local Scratch
      INTEGER KAOC(MXPNGAS),KBOC(MXPNGAS)
*. 
      INTEGER NICA(MXNDGIRR),NICB(MXNDGIRR) 
      INTEGER NIAA(MXNDGIRR),NIAB(MXNDGIRR)
      INTEGER  NKA(MXNDGIRR), NKB(MXNDGIRR)

      INTEGER ICA_EXP(MXNOP),ICB_EXP(MXNOP),
     &        IAA_EXP(MXNOP),IAB_EXP(MXNOP)
      INTEGER ICAGP(MXPNGAS),ICBGP(MXPNGAS)
      INTEGER IAAGP(MXPNGAS),IABGP(MXPNGAS)
      integer KAGP(MXPNGAS),KBGP(MXPNGAS)
      integer JAGP(MXPNGAS),JBGP(MXPNGAS)

      INTEGER IB_CA(MXNDGIRR,MXNDGIRR,MXNOP)
      INTEGER IB_CB(MXNDGIRR,MXNDGIRR,MXNOP)
      INTEGER IB_AA(MXNDGIRR,MXNDGIRR,MXNOP)
      INTEGER IB_AB(MXNDGIRR,MXNDGIRR,MXNOP)
      INTEGER IB_T(MXNDGIRR,MXNDGIRR,MXNDGIRR)
      INTEGER NIKAINTM(MXPOBS,MXNOP),NIKBINTM(MXPOBS,MXNOP)
      INTEGER NJKAINTM(MXPOBS,MXNOP),NJKBINTM(MXPOBS,MXNOP)

      INTEGER IB_C(MXNDGIRR),IB_S(MXNDGIRR)
*. Scratch through input 
      DIMENSION ISMDST(*)
*. ... Must hold excitations for all intermediate strings of given sym 
*. and all orbitals of given type
      INTEGER ICA_STR(*),ICB_STR(*),IAA_STR(*),IAB_STR(*)
      INTEGER IFOUND_T_BLOCK_ACT, IOFF_XTYPE
C
      COMMON/KKKDUMMY/LEN_C2,LEN_S2
      COMMON/CMXCJ/MXCJ
      LOGICAL CHECK_TSPLIT
C
      ISKIP              = 0
      IPLAC_S            = 0
      IPLAC_C            = 0
      IFOUND_T_BLOCK_ACT = 0
      IOFF_XTYPE         = 0
      IONE               = 1
C      
CSK   WRITE(LUWRT,*) ' IAM_BLK_ACT_C in FIND_IMAT_CALC'
CSK   CALL IWRTMAMN(IAM_BLK_ACT_C,1,NUM_BLOCKS,1,NUM_BLOCKS,LUWRT)
CSK   WRITE(LUWRT,*) ' IAM_BLK_ACT_S in FIND_IMAT_CALC'
CSK   CALL IWRTMAMN(IAM_BLK_ACT_S,1,NUM_BLOCKS,1,NUM_BLOCKS,LUWRT)
CSK   WRITE(LUWRT,*) ' IACT_C in FIND_IMAT_CALC'
CSK   CALL IWRTMAMN(IACT_C,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUWRT)
CSK   WRITE(LUWRT,*) ' IACT_S in FIND_IMAT_CALC'
CSK   CALL IWRTMAMN(IACT_S,1,NUM_BLOCKS,1,NUM_BLOCKS,LUWRT)
CSK   WRITE(LUWRT,*) ' IOFF_BLK_S in FIND_IMAT_CALC',IOFF_BLK_S
C
CSK   write(LUWRT,*) 'Coupling with type ',IXTP
C
CSK   WRITE(LUWRT,*) ' ======================== '
CSK   WRITE(LUWRT,*) '  FIND_IMAT_CALC speaking '
CSK   WRITE(LUWRT,*) ' ======================== '
CSK   write(LUWRT,*) ' Real / imaginary contribution :'
CSK   write(LUWRT,*) ' Left , Right ',IRIS,IRIC
CSK   WRITE(LUWRT,*)
CSK   WRITE(LUWRT,*) ' IAOC and IBOC in GENSIG'
CSK   CALL IWRTMAMN(IAOC,1,NGAS,1,NGAS,LUWRT)
CSK   CALL IWRTMAMN(IBOC,1,NGAS,1,NGAS,LUWRT)
CSK   WRITE(LUWRT,*)
CSK   WRITE(LUWRT,*) ' JAOC and JBOC in GENSIG'
CSK   CALL IWRTMAMN(JAOC,1,NGAS,1,NGAS,LUWRT)
CSK   CALL IWRTMAMN(JBOC,1,NGAS,1,NGAS,LUWRT)
CSK   WRITE(LUWRT,*)
CSK   WRITE(LUWRT,*) ' Type of operator '
CSK   CALL WRT_TP_GENOP(ICA,ICB,IAA,IAB,NGAS)
C
C     T operator symmetry
      IOPSM = IDBGMULT(INVELM(ICSM),ISSM)
      LI    = INVELM(ICSM)
!     WRITE(LUWRT,*) ' ISSM and ICSM and INVELM(ICSM)   ',ISSM,ICSM,LI
!     WRITE(LUWRT,*) ' operator symmetry in FIND_IMAT...',IOPSM
!     WRITE(LUWRT,*) ' IXTP                          ...',IXTP
C
C     types of strings in T in groupnotation
      CALL OCC_TO_GRP(ICA,ICAGP,1)
      CALL OCC_TO_GRP(ICB,ICBGP,1)
      CALL OCC_TO_GRP(IAA,IAAGP,1)
      CALL OCC_TO_GRP(IAB,IABGP,1)
C          OCC_TO_GRP(IOCC,IGRP,IWAY)
C
C     operators in T in expanded form
      CALL REF_OP(ICA,ICA_EXP,NCA_OP,NGAS,1,MXNOP)
      CALL REF_OP(ICB,ICB_EXP,NCB_OP,NGAS,1,MXNOP)
      CALL REF_OP(IAA,IAA_EXP,NAA_OP,NGAS,1,MXNOP)
      CALL REF_OP(IAB,IAB_EXP,NAB_OP,NGAS,1,MXNOP)
C          REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY)
C
C     type of Ka and Kb
      CALL CCEX_OCC_OCC(JAOC,KAOC,NGAS,IAA,IKA_ZERO)
      CALL CCEX_OCC_OCC(JBOC,KBOC,NGAS,IAB,IKB_ZERO)
!     print *, 'NCA_OP 2 ==> NCA_OP',NCA_OP
C
C     check that annihilation has not acted upon empty GAS 
      IF(IKA_ZERO.EQ.0.AND.IKB_ZERO.EQ.0) THEN
C
C      get group of K mapped occupations
       call occ_to_grp(KAOC,KAGP,1)
       call occ_to_grp(KBOC,KBGP,1)
!     print *, 'NCA_OP 3 ==> NCA_OP',NCA_OP
C
C      number of strings in this supergroup
       NEL_KA = IELSUM(KAOC,NGAS)
       NEL_KB = IELSUM(KBOC,NGAS)
!     print *, 'NCA_OP 3.5 ==> NCA_OP',NCA_OP
       do KSM = 1,NSMST,1
         call nst_spgrp_dbg(NGAS,KAGP,KSM,NSMST,NKA(KSM),NKADST,1)
!        call getstr2_totsm_spgp_rel(1,KAGP,NGAS,KSM,NEL_KA,
!    &                               NKA(KSM),NKADST)
!     print *, 'NCA_OP 3.x  ==> NCA_OP',NCA_OP, KSM, NEL_KA, NKA(KSM)
         call nst_spgrp_dbg(NGAS,KBGP,KSM,NSMST,NKB(KSM),NKBDST,2)
!        call getstr2_totsm_spgp_rel(2,KBGP,NGAS,KSM,NEL_KB,
!    &                               NKB(KSM),NKBDST)
!     print *, 'NCA_OP 3.xx ==> NCA_OP',NCA_OP, KSM, NEL_KB, NKB(KSM)
       end do
       NEL_CA = IELSUM(ICA,NGAS)
       NEL_CB = IELSUM(ICB,NGAS)
       NEL_AA = IELSUM(IAA,NGAS)
       NEL_AB = IELSUM(IAB,NGAS)
C       number of strings in excitation ops
        do ISYM = 1,NSMST,1

!         call getstr2_totsm_spgp_rel(1,ICAGP,NGAS,ISYM,NEL_CA,
!    &                                NICA(ISYM),NCADST)
!         call getstr2_totsm_spgp_rel(2,ICBGP,NGAS,ISYM,NEL_CB,
!    &                                NICB(ISYM),NCBDST)
!         call getstr2_totsm_spgp_rel(1,IAAGP,NGAS,ISYM,NEL_AA,
!    &                                NIAA(ISYM),NAADST)
!         call getstr2_totsm_spgp_rel(2,IABGP,NGAS,ISYM,NEL_AB,
!    &                                NIAB(ISYM),NABDST)
         call nst_spgrp_dbg(NGAS,ICAGP,ISYM,NSMST,NICA(ISYM),NCADST,1)
         call nst_spgrp_dbg(NGAS,ICBGP,ISYM,NSMST,NICB(ISYM),NCBDST,2)
         call nst_spgrp_dbg(NGAS,IAAGP,ISYM,NSMST,NIAA(ISYM),NAADST,1)
         call nst_spgrp_dbg(NGAS,IABGP,ISYM,NSMST,NIAB(ISYM),NABDST,2)

        end do
C
C      offset for symmetryblocks of T
       CALL Z_TCC_OFF_opt(IB_T,NICA,NICB,NIAA,NIAB,IOPSM)
C
C      loop over symmetry-blocks of C and sigma
C      ________________________________________
C
C      initialize offsets for sigma- and c-block loops
       I_BLK_S_SCR = 0
       I_BLK_C_SCR = 0
       I_BLK_S_SCR = IOFF_BLK_S
       I_BLK_C_SCR = IOFF_BLK_C
C
       DO IASM = 1, NSMST 
C
!       print *, ' I_BLK_S_SCR ==> ',I_BLK_S_SCR
!       print *, ' I_BLK_C_SCR ==> ',I_BLK_C_SCR
        IF( IACT_S( I_BLK_S_SCR ) .eq. 0 ) THEN
!         WRITE(LUWRT,*) ' I SKIP S_BLOCK',I_BLK_S_SCR,'with IACT_S',
!    &                     IACT_S( I_BLK_S_SCR )
          GOTO 3003
!       ELSE
!         WRITE(LUWRT,*) ' I CALC S_BLOCK',I_BLK_S_SCR,'with IACT_S',
!    &                     IACT_S( I_BLK_S_SCR )
        END IF
C
        IBSM = IDBGMULT(ISSM,INVELM(IASM))
C
C       reset offset for c-block loop
        I_BLK_C_SCR = 0
        I_BLK_C_SCR = IOFF_BLK_C
!       WRITE(LUWRT,*) 'starting c-block check at block',I_BLK_C_SCR
        DO JASM = 1, NSMST
C
         IF( IACT_C( I_BLK_C_SCR ) .eq. 0 ) GOTO 2002
C
         JBSM = IDBGMULT(ICSM,INVELM(JASM))
C
         DO 100 KASM = 1, NSMST
C         IAA_SM is symmetry of ANNIHILATOR string
          IAA_SM  = IDBGMULT(KASM,INVELM(JASM))
          ICA_SM  = IDBGMULT(IASM,INVELM(KASM))
C         to obtain the number of annnihilator strings of a given sym
C         we find the number of creator strings of conjugated
C         symmetry
          LIAA_SM = NIAA(INVELM(IAA_SM))
          LICA_SM = NICA(ICA_SM)
          LKA_SM  = NKA(KASM)
          if(LIAA_SM*LICA_SM*LKA_SM.le.0) goto 100
C
          DO KBSM = 1, NSMST
C
           IAB_SM = IDBGMULT(KBSM,INVELM(JBSM))
           ICB_SM = IDBGMULT(IBSM,INVELM(KBSM))
C
           LIAB_SM = NIAB(INVELM(IAB_SM))
           LICB_SM = NICB(ICB_SM)
C
           LKB_SM = NKB(KBSM)
!          print *, ' KBSM is ==> ',KBSM
!      call rmemchk_opt('inside find 3       ')
C
           IF(LIAB_SM*LICB_SM*LKB_SM.NE.0)THEN
C
csk             WRITE(LUWRT,*) '  GREAT, THIS IS AN ALLOWED CONNECTION'
             IF( CHECK_TC .or. CHECK_TSPLIT ) THEN 
csk               WRITE(LUWRT,*) ' this is T-block type ',IXTP
               IFOUND_T_BLOCK_ACT = 1
               GOTO 111
             END IF
             
CSK          WRITE(LUWRT,*) ' ICA_SM, ICB_SM, IAA_SM,IAB_SM =',
CSK  &                        ICA_SM, ICB_SM, IAA_SM,IAB_SM
CSK          WRITE(LUWRT,*) ' now i have to put in the connection'
             IPLAC_S = IAM_BLK_ACT_S( I_BLK_S_SCR )
             IPLAC_C = IAM_BLK_ACT_C( I_BLK_C_SCR )
             ICCTOS( IPLAC_C, IPLAC_S ) = 1 
             ICWEIGHTF( IPLAC_S ) = ICWEIGHTF( IPLAC_S ) + IONE
  111        CONTINUE
C
           END IF
C          ^ End if symmetry combinations have nonvanishing dimensions
         end do
  100    CONTINUE
C        ^ End of loop over KASM, KBSM
 2002  CONTINUE
C      ^ skip c-block not included in active-block-list; only in parallel
C        calculations
       I_BLK_C_SCR = I_BLK_C_SCR + 1
       END DO
C      ^ End of loop over JASM
 3003  CONTINUE
C      ^ skip s-block not included in node-list; only in parallel
C        calculations
       I_BLK_S_SCR = I_BLK_S_SCR + 1
       END DO
C      ^ End of loop over IASM
      END IF
C     ^ End if Ka and Kb are nontrivial strings
C
C
      IF( IFOUND_T_BLOCK_ACT .eq. 1 ) THEN
         IRIMULTIP = 1
         IF( (IRIS .eq. 1 .and. IRIC .eq. 2)
     &       .or. ( IRIS .eq. 2 .and. IRIC .eq. 1) ) THEN
           IRIMULTIP = 2
         END IF
        IOFF_XTYPE = NSPOBEX_TP * ( IRIMULTIP - 1 ) + IXTP
csk     WRITE(LUWRT,*) ' IXTP to IOFF_XTYPE, NSPOBEX_TP, IRIMULTIP',
csk  &                    IXTP, IOFF_XTYPE, NSPOBEX_TP, IRIMULTIP
        ICWEIGHTF( IOFF_XTYPE ) = ICWEIGHTF( IOFF_XTYPE ) + 1
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE FNDMND_PAR_REL(LU,LBLK,SEGMNT,NSUBMX,NSUB,ISCR,
     &                          SCR,ISCAT,SUBVAL,IBLOCKL,IBLOCKD,
     &                          NBLOCK,NTESTG)
C
C     Written by  S. Knecht         - April 30 2007 
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
*
* FIND NSUB LOWEST ELEMENTS OF VECTOR RESIDING ON FILE
* LU. ENSURE THAT NO DEGENERENCIES ARE SPLIT
*
*
* INPUT
*=======
* LU :    FILE WHERE VECTOR OF INTEREST IS RESIDING, REWOUND
* LBLK :  DEFINES FILE STRUCTURE ON FILE LU
* NSUBMX: LARGEST ALLOWED NUMBER OF SORTED ELEMENTS
*
* OUTPUT
*=======
* NSUB : ACTUAL NUMBER OF ELEMENTS OBTAINED. CAN BE SMALLER
*        THAN NSUBMX IF THE LAST ELEMENT BELONGS TO A DEGENERATE
*        SET
*ISCAT:  SCATTERING ARRAY, ISCAT(I) GIVES FULL ADRESS OF SORTED
*        ELEMENT I
*SUBVAL: VALUE OF SORTED ELEMENTS
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
#include "ipoist8.inc"
#include "wrkspc.inc"

      DIMENSION SEGMNT(*), ISCAT(*),SUBVAL(*),SCR(*),ISCR(*)
      DIMENSION IBLOCKL(NBLOCK), IBLOCKD(NBLOCK)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFF_SCR
*
*.LOOP OVER BLOCKS
*
C     write(6,*) ' FNDMND NSUBMX = ', NSUBMX
*
Cold     CALL REWINE(LU,-1)
      NTESTL = 0000
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 000
*
      IOFF_SCR = 0
      IOFF_SCR = IOFF_SCR + MY_DIA_OFF
      IBASE = 1
      LSUB = 0
*. Loop over blocks
      DO 1000 II = 1, NBLOCK
*
        IF( IBLOCKD(II) .eq. MYPROC )THEN
          LBL = IBLOCKL(II)
        ELSE
*.        useful to set all other blocks to 0?
          LBL = 0
        ENDIF
*
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Info about block ',II
          WRITE(6,*) ' Number of elements ',LBL
        END IF
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          IF( IBLOCKD(II) .eq. MYPROC )THEN
            call interface_mpi_FILE_READ_AT_r(IDIA,IOFF_SCR,SEGMNT,LBL,
     &                            ISTAT)
          ENDIF
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Elements read in '
            CALL WRTMATMN(SEGMNT,1,LBL,1,LBL,LUWRT)
          END IF
          IF(LBL .GE. 0 ) THEN
*. LOWEST ELEMENTS IN SEGMNT  ( ADD TO PREVIOUS LIST )
            MSUBMX = MIN(NSUBMX,LBL)
            IF(LBL.GE.1) THEN
              CALL SORLOW(SEGMNT,SCR(1+LSUB),ISCR(1+LSUB),LBL,
     &                    MSUBMX,MSUB,NTEST)
            ELSE
              MSUB = 0
            END IF
            DO 10 I = 1, MSUB
   10         ISCR(LSUB+I) = ISCR(LSUB+I) + IBASE - 1
* SORT COMBINED LIST
            MSUBMX = MIN(NSUBMX,LSUB+MSUB)
            IF(MSUBMX.GT.0) THEN
              CALL SORLOW(SCR,SUBVAL,ISCAT,LSUB+MSUB,MSUBMX,LLSUB,
     &                    NTEST)
            ELSE
              LLSUB = 0
            END IF
            LSUB = LLSUB
            DO 20 I = 1, LSUB
              ISCR(I+2*NSUBMX) = ISCR(ISCAT(I))
   20       CONTINUE
*
            CALL ICOPVE(ISCR(1+2*NSUBMX),ISCR(1),LSUB)
            CALL DCOPY(LSUB,SUBVAL,1,SCR,1)
*
            IF(NTEST .GE. 20 ) THEN
              WRITE(6,*) ' Lowest elements and their original place '
              WRITE(6,*) ' Number of elements obtained ', LSUB
              CALL WRTMATMN(SUBVAL,1,LSUB,1,LSUB,LUWRT)
              CALL IWRTMAMN(ISCR,1,LSUB,1,LSUB,LUWRT)
            END IF
          END IF
*
        END IF
        IOFF_SCR = IOFF_SCR + LBL
*. set to lbl to true value
        LBL = IBLOCKL(II)
        IBASE = IBASE + LBL
*
 1000 CONTINUE
*
      NSUB = LSUB
      CALL ICOPVE(ISCR,ISCAT,NSUB)
      NTEST = 00
      IF(NTEST .GE. 20) THEN
        WRITE(LUWRT,*) ' Lowest elements and their original place '
        WRITE(LUWRT,*) ' Number of elements obtained ', NSUB
        WRITE(LUWRT,*) ' WRITTEN BY PROC', MYPROC
        CALL WRTMATMN(SUBVAL,1,NSUB,1,NSUB,LUWRT)
        CALL IWRTMAMN(ISCAT,1,NSUB,1,NSUB,LUWRT)
      END IF
*
*     current dummy solution --> does this always work?
*
      IF( NSUB .lt. NSUBMX ) THEN
        IMAXVALUE = 1 000 000 000
        CALL SETVEC(SUBVAL(NSUB+1),1.0D99,NSUBMX-NSUB)
        CALL ISETVC(ISCAT(NSUB+1),IMAXVALUE,NSUBMX-NSUB)
      END IF
*
*
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'MARK  ',IDUM,'GATHER')
*
      CALL MEMMAR(KGATHERA,NMPROC*NSUBMX,'ADDL  ',2,'PARRA1')
      CALL MEMMAR(KGATHERB,NMPROC*NSUBMX,'ADDL  ',2,'PARRA2')
      CALL MEMMAR(KGATHERC,NMPROC*NSUBMX,'ADDL  ',1,'PARIA1')
      CALL MEMMAR(KGATHERD,NMPROC*NSUBMX,'ADDL  ',1,'PARIA2')
      CALL MEMMAR(KGATHERE,NSUBMX,'ADDL  ',1,'PARIA3')
*. We gather all lowest values from each node
*. and build up a combined list of those
      CALL GATHER_LOW_PAR_REL(NSUB,NSUBMX,SUBVAL,ISCAT,
     &                        WORK(KGATHERA),WORK(KGATHERB),
     &                        WORK(KGATHERC),WORK(KGATHERD),
     &                        WORK(KGATHERE),global_communicator,NTESTG)
*. update SCR1 and ISCR1
      CALL DCOPY(NSUBMX,SUBVAL,1,SCR,1)
      CALL ICOPVE(ISCAT,ISCR,NSUBMX)

      NTEST = 00
      IF(NTEST.GE.20)THEN
        WRITE(LUWRT,*)'after search: SUBVAL and ISCAT'
        CALL WRTMATMN(SUBVAL,1,NSUBMX,1,NSUBMX,LUWRT)
        CALL IWRTMAMN(ISCAT,1,NSUBMX,1,NSUBMX,LUWRT)
      END IF
      NTEST = 0
      IF(NSUB.NE.NSUBMX.AND.NTEST.GE.20)THEN
        WRITE(LUWRT,*)'Warning! NSUB is lower than NSUBMX'
        WRITE(LUWRT,*)'NSUB is set to be equal to NSUBMX'
        NSUB = NSUBMX
      END IF
*
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'FLUSM ',IDUM,'GATHER')
      WRITE(LUWRT,'(2X,A)') ' '
      WRITE(LUWRT,'(2X,A)') ' '
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GATHER_LOW_PAR_REL(NSUB,NSUBMX,SUBVAL,ISCAT,
     &                              RECVARRAY,TESTARR,
     &                              IRECVARRAY,ITESTARRAY,IRANKARR,
     &                              IGATCOMM,NTESTG)
C
C     gather all lowest elements in an array - sorted!!!
C
C     Written by  S. Knecht         - April 30 2007
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SUBVAL(*),ISCAT(*),RECVARRAY(*),TESTARR(*)
      DIMENSION IRECVARRAY(*),ITESTARRAY(*),IRANKARR(NSUBMX)
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)'SUBROUTINE GATHER_LOW_PAR entered'
      END IF
      ITESTPURPOS = 0
      IZERO = 0
      VMAX = 1.0D99
      VMIN = -1.0D99
      IF(NTEST.GE.100)
     & WRITE(LUWRT,*)'SET VMAX TO 1.0D99 and VMIN TO -1.0D99',VMAX,VMIN
      CALL SETVEC(TESTARR,VMAX,NSUBMX*NMPROC)
      CALL SETVEC(RECVARRAY,VMAX,NSUBMX*NMPROC)
      CALL ISETVC(IRECVARRAY,IZERO,NSUBMX*NMPROC)
      CALL ISETVC(ITESTARRAY,IZERO,NSUBMX*NMPROC)
      ITESTPURPOS = 0
      IF(MYPROC.NE.MASTER.AND.ITESTPURPOS.EQ.1)THEN
        CALL SETVEC(SUBVAL,VMAX,NSUBMX)
        CALL ISETVC(ISCAT,IZERO,NSUBMX)
      END IF
      ICOUNT = 1
      DO IL = 1, NSUBMX
        IRANKARR(IL) = ICOUNT
        ICOUNT = ICOUNT + 1
      ENDDO
*
*.    gather lowest elements from nodes on all nodes - we need no master
      call interface_mpi_ALLGATHER(SUBVAL,NSUBMX,RECVARRAY,NSUBMX,
     &                  IGATCOMM)
      call interface_mpi_ALLGATHER(ISCAT,NSUBMX,IRECVARRAY,NSUBMX,
     &                  IGATCOMM)
      call interface_mpi_ALLGATHER(IRANKARR,NSUBMX,ITESTARRAY,
     &                   NSUBMX,IGATCOMM)
*
      IF(MYPROC.EQ.MASTER.AND.NTEST.GE.10)THEN
        WRITE(LUWRT,*)'before search: RECVARRAY and IRECVARRAY'
        CALL WRTMATMN(RECVARRAY,1,NSUBMX*NMPROC,1,NSUBMX*NMPROC,LUWRT)
        CALL IWRTMAMN(IRECVARRAY,1,NSUBMX*NMPROC,1,NSUBMX*NMPROC,LUWRT)
        WRITE(LUWRT,*)'before search: ITESTARRAY'
        CALL IWRTMAMN(ITESTARRAY,1,NSUBMX*NMPROC,1,NSUBMX*NMPROC,LUWRT)
      END IF

*. we need the NSUBMX lowest elements from all over the world - sorted
      ISEARCH = 0
      VMINLAST = VMIN
1000  CONTINUE
*
      VMINTMP = VMAX
      IMINRANK = 0
      IMINPLACE = NSUBMX + 1

      DO 300 II = 1,NSUBMX*NMPROC
*
         TEMPMIN   = RECVARRAY(II)
         ITEMPMIN  = IRECVARRAY(II)
         ITEMPRANK = ITESTARRAY(II)
         ITEMPN    = II
*
         IF(ITEMPMIN.GT.0)THEN
           IF(TEMPMIN.LE.VMINTMP.AND.TEMPMIN.NE.VMAX)THEN
             IF(ITEMPRANK.EQ.1.AND.IMINRANK.EQ.1)THEN
               IF(TEMPMIN.LT.VMINTMP)THEN
*
                 IMINPLACE = ITEMPMIN
                 VMINTMP   = TEMPMIN
                 IMINNUMB  = ITEMPN
                 IMINRANK  = ITEMPRANK
*
               ELSE IF(TEMPMIN.EQ.VMINTMP)THEN
                 IF(ITEMPMIN.LT.IMINPLACE)THEN
*
                 IMINPLACE = ITEMPMIN
                 VMINTMP   = TEMPMIN
                 IMINNUMB  = ITEMPN
                 IMINRANK  = ITEMPRANK
*
                 END IF
               ELSE
                 GOTO 300
               END IF
             ELSE IF(ITEMPRANK.EQ.1.AND.IMINRANK.NE.1)THEN
*
               IMINPLACE = ITEMPMIN
               VMINTMP   = TEMPMIN
               IMINNUMB  = ITEMPN
               IMINRANK  = ITEMPRANK
*
             ELSE IF(ITEMPRANK.NE.1.AND.IMINRANK.EQ.1)THEN
               IF(TEMPMIN.LT.VMINTMP) THEN
*
                 IMINPLACE = ITEMPMIN
                 VMINTMP   = TEMPMIN
                 IMINNUMB  = ITEMPN
                 IMINRANK  = ITEMPRANK
*
               END IF
             ELSE IF(ITEMPRANK.NE.1.AND.IMINRANK.NE.1)THEN
               IF(ITEMPRANK.EQ.IMINRANK) THEN
                 IF(TEMPMIN.LT.VMINTMP) THEN
*
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*
                 ELSE IF(TEMPMIN.EQ.VMINTMP)THEN
                   IF(ITEMPMIN.GT.IMINPLACE)THEN
*
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*
                   END IF
                 ENDIF
               ELSE IF(ITEMPRANK.LT.IMINRANK) THEN
                 IF(TEMPMIN.LE.VMINTMP) THEN
*
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*
                 END IF
               ELSE IF(ITEMPRANK.GT.IMINRANK) THEN
                 IF(TEMPMIN.LT.VMINTMP) THEN
*
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*
                 END IF
               END IF
             END IF
           END IF
         END IF
*
*
 300  CONTINUE
*
*
      ISEARCH = ISEARCH + 1
*.    VMINTMP should be the lowest value w.r.t. ISEARCH
      VMINLAST = VMINTMP
*. test writing
      IF(NTEST.GE.10) THEN
        WRITE(LUWRT,*)'VMINLAST',VMINLAST
        WRITE(LUWRT,*)'PLACE FOR VMINLAST',IMINPLACE
        WRITE(LUWRT,*)'NUMBER',IMINNUMB
        WRITE(LUWRT,*)'IMINRANK',IMINRANK
      END IF
*. end of test writing
*. put the result to permanent storage in memory
      SUBVAL(ISEARCH) = VMINTMP
      ISCAT(ISEARCH) = IMINPLACE
      RECVARRAY(IMINNUMB) = VMAX
      IRECVARRAY(IMINNUMB) = 0
      ITESTARRAY(IMINNUMB) = -1

      IF(ISEARCH.LT.NSUBMX) GOTO 1000
1001  CONTINUE
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE H0M1TD_REL_PAR_CX(LUOUT,LUDIA,LUIN,SHIFT,VEC1,VEC2,
     &                             LUOUTLIST,LUINLIST,
     &                             NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                             MY_IOFF_LUOUT,MY_IOFF_LUDIA,
     &                             MY_IOFF_LUIN,INV,IRI_RUN,XNORM)
C
C     calculate inverted general preconditioner matrix times vector
C
C     original written by Jeppe Olsen - September 1993
C
C     adaption of sequential routine for parallel purposes  
C                      by S. Knecht   -  April 30 2007
C     extension for complex quaternion algebra ( NZ == IRI_RUN == 2 ) 
C                         S. Knecht   -  June  24 2007
C
C     disc version
C
C     vecout=  (H0 + shift )-1 * vecin
C
C      LUOUT       LUDIA        LUIN
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
*
* =====
* Input
* =====
*
* LUOUT : File to contain output vector
* LUDIA : File Containing diagonal of H0
* LUIN  : File Containing input vector
* SHIFT : constant ADDED to diagonal
*
* ======
* Output
* ======
*
* LUOUT : contains output vector, not rewinded
*
* =======
* Scratch
* =======
*
* VEC1,VEC2 : Must each be able to hold largest segment of vector
C
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUTLIST(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUDIA
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN, IOFFSET_IN_LUDIA
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
      XXXNORM = 0.0D0
      XNORM   = 0.0D0
C
      DO IRILP = 1, IRI_RUN
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
      IOFFSET_IN_LUDIA  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN    = 0
      IOFFSET_INT_LUOUT = 0
      NUM_ACTIVE_BATCH  = 0
C
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
C
        IOFFSET_IN_LUIN = MY_IOFF_LUIN + MY_VEC1_IOFF * ( IRILP - 1 )  +
     &                    IOFFSET_SCRATCH 
C
        IOFFSET_INT_IN  = 1 + MY_ACT_BLK1 * ( IRILP - 1 ) + NUM_BLK 
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK     &                   IOFFSET_IN_LUIN
C
        CALL RDVEC_BATCH_DRV4(LUIN,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                        LUINLIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC2 on LUIN'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
        IOFFSET_IN_LUDIA = MY_IOFF_LUDIA + IOFFSET_SCRATCH
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUDIA',
CSK     &                   IOFFSET_IN_LUDIA
C
C       read in batch ISBATCH from LUDIA to VEC1
C
        CALL RDVEC_BATCH_DRV5(LUDIA,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUDIA)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUDIA'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
C       calculate inverse diagonal on VEC1
C
        ILEN_BATCH =  0
        ILEN_BATCH =  LEBATCH(ISBATCH)
C
        IF( ILEN_BATCH .gt. 0 )THEN
            CALL DIAVC2(VEC1,VEC2,VEC1,SHIFT,LEBATCH(ISBATCH))
        END IF
C
C       calculate XXXNORM
C      
        XXXNORM = XXXNORM + DDOT(LEBATCH(ISBATCH),VEC1,1,VEC1,1)
C
C       write VEC1 to LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + IOFFSET_SCRATCH +
     &                      MY_VEC1_IOFF * ( IRILP - 1)  
C
        IOFFSET_INT_LUOUT = 1 + MY_ACT_BLK1 * ( IRILP - 1 ) + NUM_BLK
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches
      END DO
C     ^ loop over IRI_RUN 
C
      CAll redvec(XXXNORM,XNORM,1,2,op_MPI_SUM,global_communicator,-1)
C      
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE INPROD_B_PAR_CPLX(LUIN1,LUIN2,VEC1,VEC2,SUBSPH_R,
     &                             SUBSPH_I,NBATCH,
     &                             LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                             MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                             LUIN1LIST,LUIN2LIST,JOFF)
C
C     Written by  S. Knecht         - June 24 2007 on a sunny SUNDAY 
C                                                  in Odense
C
C**********************************************************************
C
C     calculating GENERAL dot product between two vectors on 
C     file LUIN1 resp. LUIN2
C
C     GENERAL means (in my order of vectors): 
C
C     real part == real * real + imag * imag
C     imag part == imag * real - real * imag 
C
C     NOTE: JOFF = IVEC
C
C     active blocks on the MPI-files are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH_R(*), SUBSPH_I(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN1,IOFFSET_IN_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_INT_IN1  = 0
      IOFFSET_INT_IN2  = 0
C
C
C     =====================
C     GENERAL INNER PRODUCT 
C     =====================
C
C     -------------
C      REAL * REAL
C     -------------
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C       position in file is at the end of vector JOFF - 1
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK  +
     &                     ( JOFF - 1 ) * MY_ACT_BLK2
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                       LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
        DO 100 IVEC = 1, JOFF
C
C          set new offset and zero read-in vector
C
           CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
           IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                        ( IVEC - 1 )  * MY_VEC2_IOFF
           IOFFSET_INT_IN2  = 1 + NUM_BLK  + 
     &                        ( IVEC - 1 ) * MY_ACT_BLK2
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                  IOFFSET_IN_LUIN2
C
C
C          read in batch ISBATCH from LUIN1 to VEC1
C
           CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                           LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C
C          IJ = JOFF*(JOFF-1)/2 + IVEC
C          SUBSPH_R(IJ) == VEC1 * VEC2
C
           IJ = JOFF*(JOFF-1)/2 + IVEC
C
           SUBSPH_R(IJ) = SUBSPH_R(IJ) + 
     &                  DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  100   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
C
      IOFFSET_SCRATCH = 0
      NUM_BLK = 0
      NUM_ACTIVE_BATCH = 0
C
C
C     -------------
C      REAL * IMAG
C     -------------
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C       position in file is at the end of vector JOFF - 1
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK  +
     &                     ( JOFF - 1 ) * MY_ACT_BLK2
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                       LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
        DO 200 IVEC = 1, JOFF
C
C          set new offset and zero read-in vector
C
           CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
           IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                        ( IVEC - 1 )  * MY_VEC2_IOFF    +
     &                                        MY_VEC1_IOFF 
           IOFFSET_INT_IN2  = 1 + NUM_BLK  + 
     &                        ( IVEC - 1 ) * MY_ACT_BLK2      +
     &                                       MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                  IOFFSET_IN_LUIN2
C
C
C          read in batch ISBATCH from LUIN1 to VEC1
C
           CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                           LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C
C          IJ = JOFF*(JOFF-1)/2 + IVEC
C          SUBSPH_I(IJ) == VEC1 * VEC2
C
           IJ = JOFF*(JOFF-1)/2 + IVEC
C
           SUBSPH_I(IJ) = SUBSPH_I(IJ) - 
     &                  DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  200   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
C
      IOFFSET_SCRATCH = 0
      NUM_BLK = 0
      NUM_ACTIVE_BATCH = 0
C
C
C     -------------
C      IMAG * REAL
C     -------------
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C       position in file is at the end of vector JOFF - 1
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF    +
     &                                     MY_VEC1_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK  +
     &                     ( JOFF - 1 ) * MY_ACT_BLK2      +
     &                                    MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                       LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
        DO 300 IVEC = 1, JOFF
C
C          set new offset and zero read-in vector
C
           CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
           IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                        ( IVEC - 1 )  * MY_VEC2_IOFF
CSK     &                                        MY_VEC1_IOFF 
           IOFFSET_INT_IN2  = 1 + NUM_BLK  + 
     &                        ( IVEC - 1 ) * MY_ACT_BLK2
CSK     &                                       MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                  IOFFSET_IN_LUIN2
C
C
C          read in batch ISBATCH from LUIN1 to VEC1
C
           CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                           LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C
C          IJ = JOFF*(JOFF-1)/2 + IVEC
C          SUBSPH_I(IJ) == VEC1 * VEC2
C
           IJ = JOFF*(JOFF-1)/2 + IVEC
C
           SUBSPH_I(IJ) = SUBSPH_I(IJ) +
     &                  DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  300   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
C
      IOFFSET_SCRATCH = 0
      NUM_BLK = 0
      NUM_ACTIVE_BATCH = 0
C
C
C     -------------
C      IMAG * IMAG
C     -------------
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C       position in file is at the end of vector JOFF - 1
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF    +
     &                                     MY_VEC1_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK  +
     &                     ( JOFF - 1 ) * MY_ACT_BLK2      +
     &                                    MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                       LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
        DO 400 IVEC = 1, JOFF
C
C          set new offset and zero read-in vector
C
           CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
           IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                        ( IVEC - 1 )  * MY_VEC2_IOFF    +
     &                                        MY_VEC1_IOFF 
           IOFFSET_INT_IN2  = 1 + NUM_BLK  + 
     &                        ( IVEC - 1 ) * MY_ACT_BLK2      +
     &                                       MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                  IOFFSET_IN_LUIN2
C
C
C          read in batch ISBATCH from LUIN1 to VEC1
C
           CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                           LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C
C          IJ = JOFF*(JOFF-1)/2 + IVEC
C          SUBSPH_I(IJ) == VEC1 * VEC2
C
           IJ = JOFF*(JOFF-1)/2 + IVEC
C
           SUBSPH_R(IJ) = SUBSPH_R(IJ) +
     &                  DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  400   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
C     finally ...
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE INPROD_B_PAR_CPLX2(LUIN1,LUIN2,VEC1,VEC2,SUBSPH_R,
     &                              SUBSPH_I,NBATCH,
     &                              LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                              MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                              LUIN1LIST,LUIN2LIST,NVEC,IADD)
C
C     Written by  S. Knecht         - June 25 2007 on a rainy MONDAY
C                                                  in Odense
C
C**********************************************************************
C
C     calculating GENERAL dot product between two vectors on 
C     file LUIN1 resp. LUIN2
C
C     GENERAL means (in my order of vectors): 
C
C     real part == real * real + imag * imag
C     imag part == real * imag - imag * real 
C
C     NOTE: NVEC = NVEC
C
C     active blocks on the MPI-files are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH_R(*), SUBSPH_I(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN1,IOFFSET_IN_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_INT_IN1  = 0
      IOFFSET_INT_IN2  = 0
C
C
C     =====================
C     GENERAL INNER PRODUCT 
C     =====================
C
      DO JOFF = 1, NVEC + IADD - 1
C
C     zero scratch offsets
      IOFFSET_SCRATCH  = 0
      NUM_BLK          = 0
      NUM_ACTIVE_BATCH = 0
C
      DO ISBATCH = 1, NBATCH
C
C     -------------
C      REAL * REAL
C     -------------
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C       position in file is at the end of vector JOFF - 1
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK  +
     &                     ( JOFF - 1 ) * MY_ACT_BLK2
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                       LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       set new offset and zero read-in vector
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
        IOFFSET_INT_IN2  = 1 + NUM_BLK
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                  IOFFSET_IN_LUIN2
C
C
C       read in batch ISBATCH from LUIN1 to VEC1
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       SUBSPH_R(IJ) == VEC1 * VEC2
C
CSK        IJ = ( 2 * ( JOFF - 1 ) + 1 )
        IJ = ( JOFF )
C
        SUBSPH_R(IJ) = SUBSPH_R(IJ) + 
     &                 DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  100   CONTINUE
C
C     -------------
C      REAL * IMAG
C     -------------
C
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     MY_VEC1_IOFF
C
        IOFFSET_INT_IN2  = 1 + NUM_BLK  + MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                       LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       SUBSPH_I(IJ) == VEC1 * VEC2
C
CSK        IJ = ( 2 * ( JOFF - 1 ) + 2 )
        IJ = ( JOFF )
C
        SUBSPH_I(IJ) = SUBSPH_I(IJ) +
     &                 DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  200   CONTINUE
C
C
C     -------------
C      IMAG * IMAG
C     -------------
C
C       set new offset and zero read-in vector
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF    +
     &                                     MY_VEC1_IOFF 
        IOFFSET_INT_IN1  = 1 + NUM_BLK  + 
     &                     ( JOFF - 1 ) * MY_ACT_BLK2      +
     &                                    MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
C
C       read in batch ISBATCH from LUIN1 to VEC2
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                        LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       SUBSPH_I(IJ) == VEC1 * VEC2
C
CSK        IJ = ( 2 * ( JOFF - 1 ) + 1 )
        IJ = ( JOFF )
C
        SUBSPH_R(IJ) = SUBSPH_R(IJ) + 
     &                 DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  300   CONTINUE
C
C     -------------
C      IMAG * REAL
C     -------------
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
        IOFFSET_INT_IN2  = 1 + NUM_BLK
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                  IOFFSET_IN_LUIN2
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                       LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       SUBSPH_I(IJ) == VEC1 * VEC2
C
CSK        IJ = ( 2 * ( JOFF - 1 ) + 2 )
        IJ = ( JOFF )
C
        SUBSPH_I(IJ) = SUBSPH_I(IJ) -
     &                 DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  400   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches ...
C
      END DO
C     ^ loop over NVEC + IADD - 1
C
C     finally ...
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE INPROD_B_PAR_CPLX3(LUIN1,LUIN2,VEC1,VEC2,SUBSPH_R,
     &                              SUBSPH_I,NBATCH,
     &                              LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                              MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                              LUIN1LIST,LUIN2LIST,NVEC,IVEC,
     &                              IMUSTRED,ISTRED)
C
C     Written by  S. Knecht         - June 26 2007 
C
C**********************************************************************
C
C     calculating GENERAL dot product between two vectors on 
C     file LUIN1 resp. LUIN2
C
C     GENERAL means (in my order of vectors): 
C
C     real part == real * real + imag * imag
C     imag part == real * imag - imag * real 
C
C     NOTE: NVEC = NVEC
C           IVEC = IVEC
C
C     active blocks on the MPI-files are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH_R(*), SUBSPH_I(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN1,IOFFSET_IN_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_INT_IN1  = 0
      IOFFSET_INT_IN2  = 0
C
C
C     =====================
C     GENERAL INNER PRODUCT 
C     =====================
C
      DO JOFF = 1, NVEC + IVEC
C
C
      IOFFSET_SCRATCH  = 0
      NUM_BLK          = 0
      NUM_ACTIVE_BATCH = 0
C
C
      DO ISBATCH = 1, NBATCH
C
C     -------------
C      REAL * REAL
C     -------------
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C       position in file is at the end of vector NVEC + IVEC - 1
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( NVEC + IVEC - 1 )  * MY_VEC2_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK  +
     &                     ( NVEC + IVEC - 1 ) * MY_ACT_BLK2
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                       LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       set new offset and zero read-in vector
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF
        IOFFSET_INT_IN2  = 1 + NUM_BLK   + 
     &                     ( JOFF - 1 )  * MY_ACT_BLK2
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                  IOFFSET_IN_LUIN2
C
C
C       read in batch ISBATCH from LUIN1 to VEC1
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       SUBSPH_R(IJ) == VEC1 * VEC2
C
        IJ = ( (IVEC + NVEC ) * ( IVEC + NVEC - 1 ) / 2 + JOFF )
C
CSK        XRR = DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
CSK        WRITE(LUWRT,*) 'my partial XRR for IJ',XRR,IJ
C
        SUBSPH_R(IJ) = SUBSPH_R(IJ) + 
     &                 DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  100   CONTINUE
C
C     -------------
C      REAL * IMAG
C     -------------
C
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF + MY_VEC1_IOFF
C
        IOFFSET_INT_IN2  = 1 + NUM_BLK   + 
     &                     ( JOFF - 1 )  * MY_ACT_BLK2  + MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                       LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       SUBSPH_I(IJ) == VEC1 * VEC2
C
CSK        XRI = DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
CSK        WRITE(LUWRT,*) 'my partial XRI for IJ',XRI,IJ
C
CSK        SUBSPH_I(IJ) = SUBSPH_I(IJ) + XRI
        SUBSPH_I(IJ) = SUBSPH_I(IJ) +
     &                 DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  200   CONTINUE
C
C
C     -------------
C      IMAG * IMAG
C     -------------
C
C       set new offset and zero read-in vector
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( NVEC + IVEC - 1 ) * MY_VEC2_IOFF +
     &                                           MY_VEC1_IOFF 
        IOFFSET_INT_IN1  = 1 + NUM_BLK  + 
     &                     ( NVEC + IVEC - 1 ) * MY_ACT_BLK2  +
     &                                           MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
C
C       read in batch ISBATCH from LUIN1 to VEC2
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                        LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       SUBSPH_I(IJ) == VEC1 * VEC2
C
CSK        XII = DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
CSK        WRITE(LUWRT,*) 'my partial XII for IJ',XII,IJ
C
CSK        SUBSPH_R(IJ) = SUBSPH_R(IJ) + XII
        SUBSPH_R(IJ) = SUBSPH_R(IJ) +
     &                 DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  300   CONTINUE
C
C     -------------
C      IMAG * REAL
C     -------------
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF
        IOFFSET_INT_IN2  = 1 + NUM_BLK   + 
     &                     ( JOFF - 1 )  * MY_ACT_BLK2
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                  IOFFSET_IN_LUIN2
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                       LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       SUBSPH_I(IJ) == VEC1 * VEC2
C
CSK        XIR = DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
CSK        WRITE(LUWRT,*) 'my partial XIR for IJ',XIR,IJ
C
CSK        SUBSPH_I(IJ) = SUBSPH_I(IJ) - XIR
        SUBSPH_I(IJ) = SUBSPH_I(IJ) -
     &                 DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  400   CONTINUE
C
C       keep track of memory offset and the 'reduction' counter
C
        IF( ISBATCH .eq. 1 ) THEN
           IMUSTRED = IMUSTRED + 1
           IF( IVEC .eq. 1 .AND. JOFF .eq. 1 ) ISTRED = IJ
        END IF
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches ...
C
      END DO
C     ^ loop over NVEC + IVEC
C
C     finally ...
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INPRDD_REAL_CPLX_XPRP(LUIN1,LUIN2,VEC1,VEC2,XR,XI,
     &                                 NBATCH,LBATCH,LEBATCH,I1BATCH,
     &                                 IBATCH,MY_IOFF_LUIN1,
     &                                 MY_IOFF_LUIN2,LUIN1LIST,
     &                                 LUIN2LIST,NVEC,ISSM,ISYM_CTRL,
     &                                 MZ,IXOFF,ISTATE_X)
C
C
C**********************************************************************
C
C     calculating GENERAL dot product between two vectors on 
C     file LUIN1 resp. LUIN2
C
C     GENERAL means (in my order of vectors): 
C
C     real part == real * real + imag * imag
C     imag part == real * imag - imag * real 
C
C     NOTE: NVEC = NVEC
C
C     active blocks on the MPI-files are flagged by a nonzero length
C
C     Written by  S. Knecht         - Oct 2008
C
C     Last revision:     S. Knecht  - Oct 2008
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), XR(*), XI(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*), ISYM_CTRL(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN1,IOFFSET_IN_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER(KIND=df_MPI_OFFSET_KIND) IXOFF_internal
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK          = 0
      IOFFSET_SCRATCH  = 0
      IOFFSET_IN_LUIN1 = 0
      IOFFSET_IN_LUIN2 = 0
      IOFFSET_INT_IN1  = 0
      IOFFSET_INT_IN2  = 0
!
      IXOFF_internal   =  IXOFF
C
C     zero final dot product vectors
C
      CALL DZERO(XR,NVEC)
      IF( MZ .eq. 2 ) CALL DZERO(XI,NVEC)
C     WRITE(LUWRT,*) ' XR at start'
C     CALL WRTMATMN(XR,1,NVEC,1,NVEC,LUWRT)
C     WRITE(LUWRT,*) ' XI at start'
C     CALL WRTMATMN(XI,1,NVEC,1,NVEC,LUWRT)
C
C     =====================
C     GENERAL INNER PRODUCT 
C     =====================
C
      DO 500 JOFF = 1, NVEC
C
        IF( ISYM_CTRL(JOFF) .ne. ISSM ) GOTO 500
C
!       WRITE(LUWRT,'(/A,I4,A,I4)') ' dot product for sigma vector',
!    &        ISTATE_X,' and C vector ',JOFF
C
        IOFFSET_SCRATCH  = 0
        NUM_BLK          = 0
        NUM_ACTIVE_BATCH = 0
C
      DO ISBATCH = 1, NBATCH
C
C     -------------
C      REAL * REAL
C     -------------
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C
C       set new offset
C       position in file is at the end of vector NVEC + IVEC - 1
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( IXOFF_internal - 1 )  * MY_VEC2_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK  +
     &                     ( IXOFF          - 1 )  * MY_ACT_BLK2
C
!       WRITE(LUWRT,*) 'This is my OFFSET for sigma ',IOFFSET_IN_LUIN1
!       WRITE(LUWRT,*) 'This is my block OFFSET     ',IOFFSET_INT_IN1
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                       LUIN1LIST,NUM_ACTIVE_BATCH)
C
!       WRITE(LUWRT,*) 'initial VEC2 (real) of sigma vec',ISTATE_X
!       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       set new offset and zero read-in vector
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF
        IOFFSET_INT_IN2  = 1 + NUM_BLK   + 
     &                     ( JOFF - 1 )  * MY_ACT_BLK2
C
!       WRITE(LUWRT,*) 'This is my OFFSET for cstate ',IOFFSET_IN_LUIN2
C
C       read in batch ISBATCH from LUIN1 to VEC1
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
!       WRITE(LUWRT,*) 'initial VEC1 (real) of cstate # =',JOFF
!       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
        XR(JOFF) = XR(JOFF) + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC1,1)
C
  100   CONTINUE
C
        IF( MZ .eq. 1) GOTO 444
C
C     -------------
C      REAL * IMAG
C     -------------
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF + MY_VEC1_IOFF
        IOFFSET_INT_IN2  = 1 + NUM_BLK   + 
     &                     ( JOFF - 1 )  * MY_ACT_BLK2  + MY_ACT_BLK1
C
C       WRITE(LUWRT,*) 'This is my OFFSET for JOFF',IOFFSET_IN_LUIN2
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                       LUIN2LIST,NUM_ACTIVE_BATCH)
C
C       WRITE(LUWRT,*) 'initial VEC1 (imag) of JOFF'
C       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
        XI(JOFF) = XI(JOFF) + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC1,1)
C
  200   CONTINUE
C
C     -------------
C      IMAG * IMAG
C     -------------
C
C       set new offset and zero read-in vector
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( IXOFF_internal  - 1 ) * MY_VEC2_IOFF +
     &                                               MY_VEC1_IOFF 
        IOFFSET_INT_IN1  = 1 + NUM_BLK  + 
     &                     ( IXOFF           - 1 ) * MY_ACT_BLK2  +
     &                                               MY_ACT_BLK1
C
C       WRITE(LUWRT,*) 'This is my OFFSET for ISTATE_X',IOFFSET_IN_LUIN1
C
C       read in batch ISBATCH from LUIN1 to VEC2
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                        LUIN1LIST,NUM_ACTIVE_BATCH)
C
C       WRITE(LUWRT,*) 'initial VEC2 (imag) of ISTATE_X'
C       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       SUBSPH_I(IJ) == VEC1 * VEC2
        XR(JOFF) = XR(JOFF) + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC1,1)
  300   CONTINUE
C
C     -------------
C      IMAG * REAL
C     -------------
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC2_IOFF
        IOFFSET_INT_IN2  = 1 + NUM_BLK   + 
     &                     ( JOFF - 1 )  * MY_ACT_BLK2
C
C       WRITE(LUWRT,*) 'This is my OFFSET for JOFF',IOFFSET_IN_LUIN2
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                       LUIN2LIST,NUM_ACTIVE_BATCH)
C
C       WRITE(LUWRT,*) 'initial VEC1 (real) of JOFF'
C       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       SUBSPH_I(IJ) == VEC1 * VEC2
        XI(JOFF) = XI(JOFF) - DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
  400   CONTINUE
C
  444   CONTINUE
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
      END DO
C     ^ loop over batches ...
C
!       WRITE(LUWRT,*) ' final XR(JOFF) --> JOFF,ISTATE_X',
!    &                   XR(JOFF),JOFF,ISTATE_X
!       IF( MZ .eq. 2) WRITE(LUWRT,*) ' final XI(JOFF)',XI(JOFF)
C
  500 CONTINUE
C     ^ loop over NVEC
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE INPROD_B_PAR_RL2(LUIN1,LUIN2,VEC1,VEC2,SUBSPH,NBATCH,
     &                            LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                            MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                            LUIN1LIST,LUIN2LIST,IVEC,NVEC,
     &                            IMUSTRED,ISTRED)
C
C     Written by  S. Knecht         - June 7 2007
C
C**********************************************************************
C
C     calculating dot product between two vectors on file LUIN1 resp.
C     LUIN2
C
C     NOTE: IVEC = IVEC
C           NVEC = NVEC
C
C     active blocks on the MPI-files are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN1,IOFFSET_IN_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_INT_IN1  = 0
      IOFFSET_INT_IN2  = 0
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C       position in file is at the end of vector JOFF - 1
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( NVEC + IVEC - 1 )  * MY_VEC1_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK  +
     &                     ( NVEC + IVEC - 1 ) * MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                  IOFFSET_IN_LUIN1
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                       LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
        DO 100 JVEC = 1, NVEC+IVEC
C
C          set new offset and zero read-in vector
C
           CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
           IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                        ( JVEC - 1 )  * MY_VEC1_IOFF
           IOFFSET_INT_IN2  = 1 + NUM_BLK  + 
     &                        ( JVEC - 1 ) * MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                  IOFFSET_IN_LUIN2
C
C
C          read in batch ISBATCH from LUIN1 to VEC1
C
           CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                           LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C
C          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
C          SUBSPH(IJ) == VEC1 * VEC2
C
           IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
C
           SUBSPH(IJ) = SUBSPH(IJ) + 
     &                  DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
C          keep track of memory offset and the 'reduction' counter
C
           IF( ISBATCH .eq. 1 ) THEN
             IMUSTRED = IMUSTRED + 1
             IF(IVEC.EQ.1.AND.JVEC.EQ.1) ISTRED = IJ
           END IF
C
C
  100   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE LUCI_TASK_NODES(CIRUN)
C
C     Introduce LUCIAREL action - Node part
C
C     Written by  S. Knecht         - April 13 2007 (Friday!)
C
C**********************************************************************
#include "implicit.h"
#include "parluci.h"
C
      character*6 CIRUN
*
      NTEST = 0
*
      if (NTEST.ge.1) then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' ////////////////////////////////// '
        write(MSLVOUT_REL,*) ' //                              // '
        write(MSLVOUT_REL,*) ' //                              // '
        write(MSLVOUT_REL,*) ' //  KR-MCSCF calling LUCIAREL   // '
        write(MSLVOUT_REL,*) ' //                              // '
        write(MSLVOUT_REL,*) ' //                              // '
        write(MSLVOUT_REL,*) ' ////////////////////////////////// '
      end if
*
      if (CIRUN.eq.'NDET  ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' NDET   flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) ' The number of CI dets will be computed.'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'CIINII') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' CIINII flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Determining configurational start 
     & guess'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'KR-CI ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' KR-CI  flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Running large-scale KR-CI calculation'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'RSTRMC') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' RSTRMC flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Determining info for restart of MCSCF.'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'SIGMA ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' SIGMA  flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Computing a sigma vector.'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'SIGMAD') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' SIGMAD flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Computing a sigma vector.'
        write(MSLVOUT_REL,*) '     Saving H diagonal also.'
      else if (CIRUN.eq.'DENS1 ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' DENS1  flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Computing only 1-particle density 
     & matrices.'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'DENS2 ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' DENS2  flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Computing 1- and 2-particle dens. 
     & matrices.'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'ANALYZ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' ANALYZ flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Analyzing the current CI vector(s).'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'DIAG  ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' DIAG   flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Computing CI diagonal.'
        call quit('NOT IMPLEMENTED YET!')
      end if
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GSTTBLD_PAR(CTT,IATP,IASM,IBTP,IBSM,
     &                       NSASO,NSBSO,PSSIGN,IDC,
     &                       PLSIGN,LUIN,SCR,NSMST,ISCALE,SCLFAC,
     &                       IOFFSET_IN_LUIN,ILEN)
C
C     Adapted by  S. Knecht         - June 12 2007
C
C**********************************************************************
C
C     obtain  determinant block (iatp iasm, ibtp ibsm )
C     from vector packed in combination format according to IDC
C
C.    If ISCALE = 1, the routine scales and returns the block
C     in determinant normalization, and SCLFAC = 1.0D0
C
C     If ISCALE = 0, the routine does not perform any overall
C     scaling, and a scale factor is returned in SCLFAC
C
C     IF ISCALE = 0, zero blocks are not set explicitly to zero,
C     instead  zero is returned in SCLFAC
C
C     ISCALE, SCLFAC added May 97
C
C     Simplified version working only for vectors on disc
C
C
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "ipoist8.inc"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION CTT(*), SCR(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN
      DIMENSION NSASO( NSMST, * ), NSBSO( NSMST, * )
C
      NTEST = 000
C
      IF(NTEST.GE.100) THEN
        write(LUWRT,*) ' GSTTBLD_PAR, IATP,IASM,IBTP,IBSM,ISCALE'
        write(LUWRT,*)            IATP,IASM,IBTP,IBSM,ISCALE
        WRITE(LUWRT,*) ' LUIN = ', LUIN
      END IF
* =================
* Read in from disc
* =================
      call interface_mpi_FILE_READ_AT_r(LUIN,IOFFSET_IN_LUIN,
     &                     SCR,ILEN,ISTAT)
*
      NAST = NSASO(IASM,IATP)
      NBST = NSBSO(IBSM,IBTP)
      IF(ILEN .ne. 0 ) THEN
        CALL SDCMRF(CTT,SCR,2,IATP,IBTP,IASM,IBSM,NAST,NBST,
     &       IDC,PSSIGN,PLSIGN,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
      ELSE
        SCLFAC = 0.0D0
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GSTTBLD_PAR2(CTT,LUIN,IOFFSET_IN_LUIN,ILEN)
C
C     Adapted by  S. Knecht         - July 09 2007
C
C**********************************************************************
C
C     GET BLOCK FROM DISC using MPI-I/O
C
C     Last revision:     S. Knecht       - July  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION CTT(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN
*
      NTEST = 000
*
*     read in from disc
      call interface_mpi_FILE_READ_AT_r(LUIN,IOFFSET_IN_LUIN,
     &                      CTT,ILEN,ISTAT)
*
*
      IF( NTEST .GE. 100 )THEN
        CALL WRTMATMN(CTT,1,ILEN,1,ILEN,LUWRT)
      END IF
*
      NTEST = 000
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE P1_B_PAR_RL_1(VEC1,VEC2,SUBSPH,EIGAPR,RNRM,EIGSHF,
     &                      EIG,TEST,E_CONV,RTCNV,CONVER,ITER,MAXIT,
     &                      IROOT,LUIN1LIST,LUIN2LIST,LUOUTLIST,
     &                      NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                      MY_IOFF_LUIN1,MY_IOFF_LUIN2,MY_IOFF_LUOUT,
     &                      SCRRED,NVEC,LUIN1,LUIN2,LUOUT,MAXVEC)
C
C     Written by  S. Knecht         - June 4 2007
C
C**********************************************************************
C
C     calculating dot product between two vectors on file LUIN1 resp.
C     LUIN2
C
C     NOTE: IROOT = IROOT
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*), LUOUTLIST(*)
      DIMENSION RNRM(MAXIT,*), EIG(MAXVEC,MAXIT)
      DIMENSION SCRRED(*)
      LOGICAL CONVER, RTCNV(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN1,IOFFSET_IN_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN1   = 0
      IOFFSET_INT_IN2   = 0
      IOFFSET_INT_LUOUT = 0
C
      RNORM     = 0.0D0
      REDSCRVAR = 0.0D0
      CALL DZERO(SCRRED,NVEC)
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 100 IVEC = 1, NVEC
C
          FACTOR = SUBSPH( ( IROOT-1 ) * NVEC + IVEC )
CSK          WRITE(LUWRT,*) ' THIS IS THE FIRST FACTOR from AVEC',FACTOR
C
C           set new offset
C
C           position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )  * MY_VEC1_IOFF
          IOFFSET_INT_IN1  = 1 + NUM_BLK   +
     &                     ( IVEC - 1 )  * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                   IOFFSET_IN_LUIN1
C
          IF( IVEC .eq. 1 ) THEN
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
C
          ELSE
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC1 on LUIN1'
CSK            CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C           VEC2 == VEC2 + VEC1 * FACTOR
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
          END IF
C
C         new scal-FACTOR
C
          FACTOR = - EIGAPR * FACTOR
CSK          WRITE(LUWRT,*) ' THIS IS THE 2nd FACTOR from AVEC',FACTOR
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )    * MY_VEC1_IOFF
          IOFFSET_INT_IN2  = 1 + NUM_BLK   + 
     &                     ( IVEC - 1 )    * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                   IOFFSET_IN_LUIN2
C
C
C         read in batch ISBATCH from LUIN2 to VEC1
C
          CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                          LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
C
C         VEC2 == VEC2 + VEC1 * scal-FACTOR 
C
          CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
  100   CONTINUE
C
C       calculate partial RNORM
C
        REDSCRVAR = REDSCRVAR + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC2,1)
C
C       write VEC2 to LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT = 1 + NUM_BLK
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUOUT',
CSK     &                  IOFFSET_LUOUT
CSK        WRITE(LUWRT,*) 'final VEC2 to write on LUOUT'
C
C       needed for unknown reason. compiler optimization problem?
C
CSK          WRITE(LUWRT,*) '  ( HX - EX ) '
CSK          CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                LUWRT)
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
C     communicate REDSCRVAR to get full RNORM
C
      CAll redvec(REDSCRVAR,SCRRED,1,2,op_MPI_SUM,
     &            global_communicator,-1)
      CALL DCOPY(1,SCRRED,1,REDSCRVAR,1)
CSK      WRITE(LUWRT,*) 'RNORM**2',REDSCRVAR
      RNORM = SQRT(REDSCRVAR)
C
      RNRM(ITER-1,IROOT) = RNORM
C
C     check for convergence
C
      IF(RNORM.LT. TEST .OR. (ITER.GT.2.AND. 
     &  ABS(EIG(IROOT,ITER-2)-EIG(IROOT,ITER-1)).LT.E_CONV)) THEN
C
        RTCNV(IROOT) = .TRUE.
      ELSE
C
        RTCNV(IROOT) = .FALSE.
        CONVER       = .FALSE.
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE P1_B_PAR_CX_1(VEC1,VEC2,SUBSPH,SUBSPH_2,EIGAPR,RNRM,
     &                        EIGSHF,RTCNV,CONVER,CONVLIM,ITER,NROOT,
     &                        MAXIT,IROOT,LUIN1LIST,LUIN2LIST,LUOUTLIST,
     &                        NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                        MY_IOFF_LUIN1,MY_IOFF_LUIN2,MY_IOFF_LUOUT,
     &                        SCRRED_R,SCRRED_I,NVEC,NVEC_D,
     &                        LUIN1,LUIN2,LUOUT)
C
C     Written by  S. Knecht         - June 25 2007
C
C**********************************************************************
C
C     calculating residue and HX - EIGAPR*X vector using LUIN1 and LUIN2
C
C     NOTE: IROOT  = IROOT
C           NVEC   = NVEC
C           NVEC_D = NVEC * 2
C
C     general routine for calculations in the complex quaternion algebra
C     regime
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*), SUBSPH_2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*), LUOUTLIST(*)
      DIMENSION RNRM(NROOT,MAXIT)
      DIMENSION SCRRED_R(*), SCRRED_I(*)
      LOGICAL CONVER, RTCNV(*)
C     LOGICAL CONVER, RTCNV(100)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN1,IOFFSET_IN_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN1   = 0
      IOFFSET_INT_IN2   = 0
      IOFFSET_INT_LUOUT = 0
C
      RNORM     = 0.0D0
      REDSCRVAR_R = 0.0D0
      REDSCRVAR_I = 0.0D0
      CALL DZERO(SCRRED_R,NVEC_D)
      CALL DZERO(SCRRED_I,NVEC_D)
      FACTOR   = 0.0D0
      FACTOR_N = 0.0D0
C
C
C     ================
C      COMPLEX VECTOR
C     ================
C
C     .................
C      REAL PART FIRST
C     .................
C
C
C     ===========================
C     Real part of HX - EIGAPR*X
C     ===========================
C
      DO IVEC = 1, NVEC
           SUBSPH_2(2*(IVEC-1)+1) =   SUBSPH((IROOT-1)*2*NVEC+IVEC)
           SUBSPH_2(2*(IVEC-1)+2) = - SUBSPH((IROOT-1)*2*NVEC+NVEC+IVEC)
      END DO
C
C
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 100 IVEC = 1, NVEC_D
C
          FACTOR = SUBSPH_2( IVEC )
CSK          WRITE(LUWRT,*) ' THIS IS THE FIRST FACTOR from WORK',FACTOR
C
C           set new offset
C
C           position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )  * MY_VEC1_IOFF
          IOFFSET_INT_IN1  = 1 + NUM_BLK   +
     &                     ( IVEC - 1 )  * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                   IOFFSET_IN_LUIN1
C
          IF( IVEC .eq. 1 ) THEN
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
C
          ELSE
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC1 on LUIN1'
CSK            CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C           VEC2 == VEC2 + VEC1 * FACTOR
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
          END IF
C
 100    CONTINUE
C
C       new scal-FACTOR
C
        FACTOR = - EIGAPR
CSK        WRITE(LUWRT,*) ' THIS IS THE 2nd FACTOR',FACTOR
C
C       scale VEC2 ... == -EIGAPR * X
C
        CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
C
C
        DO 200 IVEC = 1, NVEC_D
C
C         set new offset, new scaling factor and zero read-in vector
C
          FACTOR_N = SUBSPH_2( IVEC )
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )    * MY_VEC1_IOFF
          IOFFSET_INT_IN2  = 1 + NUM_BLK   + 
     &                     ( IVEC - 1 )    * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                   IOFFSET_IN_LUIN2
C
C
C         read in batch ISBATCH from LUIN2 to VEC1
C
          CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                          LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
C
C                  -EIGAPR * X + HX * FACTOR
C
C         VEC2  == VEC2 + VEC1 * scal-FACTOR 
C
          CALL DAXPY(LEBATCH(ISBATCH),FACTOR_N,VEC1,1,VEC2,1)
C
  200   CONTINUE
C
C       calculate partial RNORM ( real part )
C
        REDSCRVAR_R = REDSCRVAR_R + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC2,1)
C
C       write VEC2 to LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT = 1 + NUM_BLK
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUOUT',
CSK     &                  IOFFSET_LUOUT
CSK        WRITE(LUWRT,*) 'final VEC2 to write on LUOUT'
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches ( real part )
C
C     ..............
C        IMAG PART
C     ..............
C
      IOFFSET_SCRATCH  = 0
      NUM_BLK          = 0
      NUM_ACTIVE_BATCH = 0
C
C
C     ===========================
C     imag part of HX - EIGAPR*X
C     ===========================
C
      DO IVEC = 1, NVEC
        SUBSPH_2(2*(IVEC-1)+1) = SUBSPH((IROOT-1)*2*NVEC+NVEC+IVEC)
        SUBSPH_2(2*(IVEC-1)+2) = SUBSPH((IROOT-1)*2*NVEC     +IVEC)
      END DO
C
C
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 300 IVEC = 1, NVEC_D
C
          FACTOR = SUBSPH_2( IVEC )
CSK          WRITE(LUWRT,*) ' THIS IS THE FIRST FACTOR from WORK',FACTOR
C
C           set new offset
C
C           position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )  * MY_VEC1_IOFF
          IOFFSET_INT_IN1  = 1 + NUM_BLK   +
     &                     ( IVEC - 1 )  * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                   IOFFSET_IN_LUIN1
C
          IF( IVEC .eq. 1 ) THEN
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
C
          ELSE
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC1 on LUIN1'
CSK            CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C           VEC2 == VEC2 + VEC1 * FACTOR
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
          END IF
C
 300    CONTINUE
C
C       new scal-FACTOR
C
        FACTOR = - EIGAPR
CSK        WRITE(LUWRT,*) ' THIS IS THE 2nd FACTOR',FACTOR
C
C       scale VEC2 ... == -EIGAPR * X
C
        CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
C
C
        DO 400 IVEC = 1, NVEC_D
C
C         set new offset, new scaling factor and zero read-in vector
C
          FACTOR_N = SUBSPH_2( IVEC )
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )    * MY_VEC1_IOFF
          IOFFSET_INT_IN2  = 1 + NUM_BLK   + 
     &                     ( IVEC - 1 )    * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                   IOFFSET_IN_LUIN2
C
C
C         read in batch ISBATCH from LUIN2 to VEC1
C
          CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                          LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
C
C                  -EIGAPR * X + HX * FACTOR
C
C         VEC2  == VEC2 + VEC1 * scal-FACTOR 
C
          CALL DAXPY(LEBATCH(ISBATCH),FACTOR_N,VEC1,1,VEC2,1)
C
  400   CONTINUE
C
C       calculate partial RNORM ( imag part )
C
        REDSCRVAR_I = REDSCRVAR_I + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC2,1)
C
C       write VEC2 to LUOUT - imag part
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + MY_VEC1_IOFF + 
     &                      IOFFSET_SCRATCH
C
        IOFFSET_INT_LUOUT = 1 + NUM_BLK + MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUOUT',
CSK     &                  IOFFSET_LUOUT
CSK        WRITE(LUWRT,*) 'final VEC2 to write on LUOUT'
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches ( imag part )
C
C     communicate real and imag part of REDSCRVAR to get full RNORM
C
      CAll redvec(REDSCRVAR_R,SCRRED_R,1,2,op_MPI_SUM,
     &                global_communicator,-1)
C
      CAll redvec(REDSCRVAR_I,SCRRED_I,1,2,op_MPI_SUM,
     &                global_communicator,-1)
C
      CALL DCOPY(1,SCRRED_R,1,REDSCRVAR_R,1)
      CALL DCOPY(1,SCRRED_I,1,REDSCRVAR_I,1)
CSK      WRITE(LUWRT,*) 'RNORM**2',REDSCRVAR
      RNORM = SQRT( REDSCRVAR_R + REDSCRVAR_I )
C
      RNRM(IROOT,ITER) = RNORM
C
C     check for convergence
C
      IF( RNORM .lt. CONVLIM ) THEN
C
        RTCNV(IROOT) = .TRUE.
C
      ELSE
C
        RTCNV(IROOT) = .FALSE.
        CONVER       = .FALSE.
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE P1_B_PAR_RL_2(VEC1,VEC2,SUBSPH,SHIFT,IROOT,
     &                         LUINLIST,LUOUT1LIST,LUOUT2LIST,
     &                         LUOUT3LIST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         MY_IOFF_LUIN,MY_IOFF_LUOUT1,
     &                         MY_IOFF_LUOUT2,MY_IOFF_LUOUT3,
     &                         MY_IOFF_LUDIA,
     &                         NVEC,LUIN,
     &                         LUOUT1,LUOUT2,LUOUT3,LUDIA,INV)
C
C     Written by  S. Knecht         - June 4 2007
C
C**********************************************************************
C
C     calculating dot product between two vectors on file LUIN1 resp.
C     LUIN2
C
C     NOTE: IROOT = IROOT
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUT1LIST(*), LUOUT2LIST(*)
      DIMENSION LUOUT3LIST(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUOUT1
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUOUT2
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUOUT3
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUDIA
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN, IOFFSET_LUOUT1
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_LUOUT2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_LUOUT3
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUDIA
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK            = 0
      IOFFSET_SCRATCH    = 0
      IOFFSET_IN_LUIN    = 0
      IOFFSET_LUOUT1     = 0
      IOFFSET_LUOUT2     = 0
      IOFFSET_LUOUT3     = 0
      IOFFSET_IN_LUDIA   = 0
      IOFFSET_INT_IN     = 0
      IOFFSET_INT_LUOUT1 = 0
      IOFFSET_INT_LUOUT2 = 0
      IOFFSET_INT_LUOUT3 = 0
C
      REDSCRVAR = 0.0D0
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 100 IVEC = 1, NVEC
C
          FACTOR = SUBSPH( ( IROOT-1 ) * NVEC + IVEC )
C
C           set new offset
C
C           position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )   * MY_VEC1_IOFF
          IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                     ( IVEC - 1 )   * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK     &                   IOFFSET_IN_LUIN
C
          IF( IVEC .eq. 1 ) THEN
C
            CALL RDVEC_BATCH_DRV4(LUIN,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                           LUINLIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC2 on LUIN'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
CSK            WRITE(LUWRT,*) 'scaling factor for this vector',FACTOR
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
CSK            WRITE(LUWRT,*) 'this is VEC2 after first scaling in P1...2'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
          ELSE
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
            CALL RDVEC_BATCH_DRV4(LUIN,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                           LUINLIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC1 on LUIN'
CSK            CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C           VEC2 == VEC2 + VEC1 * FACTOR
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
CSK            WRITE(LUWRT,*) 'final VEC2 after DAXPY in P1...2'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
          END IF
C         ^ IVEC == 1 ?
C
  100   CONTINUE
C
CSK        WRITE(LUWRT,*) ' WE SURVIVED FIRST PART IN ..._2'
CSK        WRITE(LUWRT,*) ' SHIFT IS',SHIFT
C
C       calculate inverse diagonal on VEC1
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       new offset for file containing diagonal
C
        IOFFSET_IN_LUDIA = MY_IOFF_LUDIA + IOFFSET_SCRATCH
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUDIA',
CSK     &                   IOFFSET_IN_LUDIA
C
C       read in batch ISBATCH from LUDIA to VEC1
C
        CALL RDVEC_BATCH_DRV5(LUDIA,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUDIA)
C
CSK          WRITE(LUWRT,*) 'current X in ..._2'
CSK          CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
C         scratch writing, needed!!! reason???
C
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
        IF( LEBATCH(ISBATCH) .gt. 0 )THEN
          IF( INV .ne. 0 ) THEN
CSK            WRITE(LUWRT,*) 'callign DIAVC2 in ..._2'
            CALL DIAVC2(VEC1,VEC2,VEC1,SHIFT,LEBATCH(ISBATCH))
          ELSE
            CALL VVTOV(VEC1,VEC2,VEC1,LEBATCH(ISBATCH))
            CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,SHIFT,LEBATCH(ISBATCH))
          END IF
        END IF
CSK          WRITE(LUWRT,*) ' Lu4 after H0...'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
CSK          WRITE(LUWRT,*) ' Lu5 after H0...'
CSK          CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
C
C
C       write VEC1 to LUOUT1 and VEC2 to LUOUT2
C
        IOFFSET_LUOUT1     = MY_IOFF_LUOUT1  + IOFFSET_SCRATCH
        IOFFSET_LUOUT2     = MY_IOFF_LUOUT2  + IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT1 = 1 + NUM_BLK
        IOFFSET_INT_LUOUT2 = 1 + NUM_BLK
C
C       VEC1
        CALL WTVEC_BATCH_DRV4(LUOUT1,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT1,IOFFSET_INT_LUOUT1,
     &                       LUOUT1LIST,NUM_ACTIVE_BATCH)
C       VEC2
        CALL WTVEC_BATCH_DRV4(LUOUT2,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT2,IOFFSET_INT_LUOUT2,
     &                       LUOUT2LIST,NUM_ACTIVE_BATCH)
C
C
CSK        WRITE(LUWRT,*) 'THIS IS my partial REDSCRVAR',REDSCRVAR
CSK        WRITE(LUWRT,*) 'THIS IS LEBATCH(ISBATCH)',LEBATCH(ISBATCH)
C       calculate partial GAMMA
        IF( LEBATCH(ISBATCH) .gt. 0 ) THEN
          REDSCRVAR = REDSCRVAR + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC1,1)
        END IF
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
C     communicate REDSCRVAR to get full GAMMA
C
      CAll redvec(REDSCRVAR,GAMMA,1,2,op_MPI_SUM,global_communicator,-1)
C
CSK      WRITE(LUWRT,*) 'THIS IS GAMMA',GAMMA
C
C     continue with VNORM ...
C
C     reset scratch offsets
      NUM_BLK            = 0
      IOFFSET_SCRATCH    = 0
      REDSCRVAR = 0.0D0
CSK      WRITE(LUWRT,*) 'THIS IS REDSCRVAR',REDSCRVAR
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       read VEC1 from LUOUT2 and VEC2 from LUOUT1
C
        IOFFSET_LUOUT1     = MY_IOFF_LUOUT1  + IOFFSET_SCRATCH
        IOFFSET_LUOUT2     = MY_IOFF_LUOUT2  + IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT1 = 1 + NUM_BLK
        IOFFSET_INT_LUOUT2 = 1 + NUM_BLK
C
C       VEC1
        CALL RDVEC_BATCH_DRV4(LUOUT2,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_LUOUT2,IOFFSET_INT_LUOUT2,
     &                        LUOUT2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUOUT2 in P1..._2 again'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)
C       VEC2
        CALL RDVEC_BATCH_DRV4(LUOUT1,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_LUOUT1,IOFFSET_INT_LUOUT1,
     &                        LUOUT1LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) ' VEC2 before DAXPY call'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)
CSK        WRITE(LUWRT,*) ' VEC1 before DAXPY call'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)
C
C       VEC2 == VEC2 + VEC1 * FACTOR
C
        CALL DAXPY(LEBATCH(ISBATCH),-GAMMA,VEC1,1,VEC2,1)
C
CSK        WRITE(LUWRT,*) ' VEC2 after DAXPY call'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)
CSK        WRITE(LUWRT,*) ' VEC1 after DAXPY call'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)
C
C       calculate partial VNORM_Q
C
        REDSCRVAR = REDSCRVAR + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC2,1)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      VNORM_Q = 0.0D0
      VNORM   = 0.0D0
C
C     communicate REDSCRVAR to get full VNORM_Q
C
      CAll redvec(REDSCRVAR,VNORM_Q,1,2,op_MPI_SUM,
     &            global_communicator,-1)
C
C     is X an eigen vector for (H0 - 1 ) - 1 ???
C
      VNORM = SQRT(VNORM_Q)
C
CSK      WRITE(LUWRT,*) 'GAMMA ',GAMMA
CSK      WRITE(LUWRT,*) 'VNORM ',VNORM
C
      IF( VNORM .GT. 1.0D-7 ) THEN
        IOLSAC = 1
      ELSE
        IOLSAC = 0
      END IF
      IF(IOLSAC .EQ. 1 ) THEN
C
CSK        WRITE(LUWRT,*) ' Olsen correction active'
        DELTA = 0.0D0 
C
C       continue with DELTA ...
C
C       reset scratch offsets
        NUM_BLK            = 0
        IOFFSET_SCRATCH    = 0
        REDSCRVAR = 0.0D0
C
        DO ISBATCH = 1, NBATCH
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
          CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C         read VEC1 from LUOUT2 and VEC2 from LUOUT3
C
          IOFFSET_LUOUT3     = MY_IOFF_LUOUT3  + IOFFSET_SCRATCH
          IOFFSET_LUOUT2     = MY_IOFF_LUOUT2  + IOFFSET_SCRATCH
          IOFFSET_INT_LUOUT3 = 1 + NUM_BLK
          IOFFSET_INT_LUOUT2 = 1 + NUM_BLK
C
C         VEC1
          CALL RDVEC_BATCH_DRV4(LUOUT2,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT2,IOFFSET_INT_LUOUT2,
     &                          LUOUT2LIST,NUM_ACTIVE_BATCH)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUOUT2'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                   LUWRT)
C         VEC2
          CALL RDVEC_BATCH_DRV4(LUOUT3,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT3,IOFFSET_INT_LUOUT3,
     &                          LUOUT3LIST,NUM_ACTIVE_BATCH)
C
CSK         WRITE(LUWRT,*) 'initial VEC2 on LUOUT3'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
C         calculate partial DELTA
C
          REDSCRVAR = REDSCRVAR + DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
C         keep track of correct offset
          IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
          NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
        END DO
C
C
C       communicate REDSCRVAR to get full DELTA
C
        CAll redvec(REDSCRVAR,DELTA,1,2,op_MPI_SUM,
     &            global_communicator,-1)
C
C
CSK        WRITE(LUWRT,*) ' THIS IS DELTA'
C
        FACTOR = - DELTA / GAMMA
CSK        WRITE(LUWRT,*) 'FACTOR, DELTA, GAMMA',FACTOR, DELTA, GAMMA
C
C       reset scratch offsets
        NUM_BLK            = 0
        IOFFSET_SCRATCH    = 0
C
        DO ISBATCH = 1, NBATCH
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
          CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C         read VEC1 from LUOUT1 and VEC2 from LUOUT3
C
          IOFFSET_LUOUT3     = MY_IOFF_LUOUT3  + IOFFSET_SCRATCH
          IOFFSET_LUOUT1     = MY_IOFF_LUOUT1  + IOFFSET_SCRATCH
          IOFFSET_INT_LUOUT3 = 1 + NUM_BLK
          IOFFSET_INT_LUOUT1 = 1 + NUM_BLK
C
C         VEC1
          CALL RDVEC_BATCH_DRV4(LUOUT1,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT1,IOFFSET_INT_LUOUT1,
     &                          LUOUT1LIST,NUM_ACTIVE_BATCH)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUOUT1'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                   LUWRT)
C         VEC2
          CALL RDVEC_BATCH_DRV4(LUOUT3,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT3,IOFFSET_INT_LUOUT3,
     &                          LUOUT3LIST,NUM_ACTIVE_BATCH)
C
CSK         WRITE(LUWRT,*) 'initial VEC2 on LUOUT3'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
C         VEC2 == VEC2 + VEC1 * FACTOR
C
          CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
C         write VEC2 on LUOUT3
C      
          IOFFSET_LUOUT3     = MY_IOFF_LUOUT3  + IOFFSET_SCRATCH
          IOFFSET_INT_LUOUT3 = 1 + NUM_BLK
C
C         needed?
C
CSK          CALL ISETVC(LUOUT3LIST(IOFFSET_INT_LUOUT3),0,NUM_ACTIVE_BATCH)
C
CSK         WRITE(LUWRT,*) 'final VEC2 to write on LUOUT3'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
          CALL WTVEC_BATCH_DRV4(LUOUT3,VEC2,LBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         IOFFSET_LUOUT3,IOFFSET_INT_LUOUT3,
     &                         LUOUT3LIST,NUM_ACTIVE_BATCH)
C
C         keep track of correct offset
          IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
          NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
        END DO
C
      END IF
C     ^ IOLSAC ?
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE P1_B_PAR_RL_3(VEC1,VEC2,SUBSPH,
     &                         LUIN1LIST,LUIN2LIST,LUOUTLIST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                         MY_IOFF_LUOUT,SCRRED,NVEC,IADD,
     &                         LUIN1,LUIN2,LUOUT)
C
C     Written by  S. Knecht         - June 4 2007
C
C**********************************************************************
C
C     calculating dot product between two vectors on file LUIN1 resp.
C     LUIN2
C
C     NOTE: IROOT = IROOT
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*), LUOUTLIST(*)
      DIMENSION SCRRED(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN1,IOFFSET_IN_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
      THR_CVEC = 1.0d-04
!     THR_CVEC = -1.0d-03
!     write(luwrt,*) 'THR_CVEC ==>', THR_CVEC

      IROUND = 0
 10   CONTINUE 
      IROUND = IROUND + 1
      I_ZERO = 0
      I_REMOVED = 0
    
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN1   = 0
      IOFFSET_INT_IN2   = 0
      IOFFSET_INT_LUOUT = 0
C
      REDSCRVAR = 0.0D0
      CALL DZERO(SCRRED,NVEC+IADD)
      CALL DZERO(SUBSPH,NVEC+IADD)
CSK      WRITE(LUWRT,*) ' NVEC + IADD - 1',  NVEC + IADD - 1
CSK      WRITE(LUWRT,*) ' LUIN1,LUIN2,LUOUT', LUIN1,LUIN2,LUOUT
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C           set new offset
C
C       position in file is at the end of vector IVEC - 1
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
        IOFFSET_INT_IN2  = 1 + NUM_BLK 
C
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN2 in P1..._3 100',
CSK     &                   IOFFSET_IN_LUIN2
C
         CALL RDVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                         LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK         WRITE(LUWRT,*) 'initial VEC2 on LUIN2 in P1..._3 100'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)

C
        DO 100 IVEC = 1, NVEC + IADD - 1
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C           set new offset
C
C           position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )  * MY_VEC1_IOFF
          IOFFSET_INT_IN1  = 1 + NUM_BLK   +
     &                     ( IVEC - 1 )  * MY_ACT_BLK1
C
CSK          WRITE(LUWRT,*) 'This is my OFFSET for LUIN1 in P1..._3 100',
CSK     &                   IOFFSET_IN_LUIN1
CSK          WRITE(LUWRT,*) 'This is my INT_OFFSET for LUIN1 in 
CSK     &                    P1..._3 100',IOFFSET_INT_IN1
CSK          WRITE(LUWRT,*) 'THIS IS MY LU1LIST inside P1_B_PAR_RL_3 100'
CSK          CALL IWRTMAMN(LUIN1LIST,1,IALL_LU1,1,IALL_LU1,LUWRT)
C
          CALL RDVEC_BATCH_DRV4(LUIN1,VEC1,LBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                         LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUIN1 in P1..._3 100'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
          SCRRED(IVEC) = SCRRED(IVEC) + 
     &                   DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  100   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
C     communicate SCRRED to get full OVERLAP matrix
C
      CAll redvec(SCRRED,SUBSPH,NVEC+IADD-1,2,op_MPI_SUM,
     &                global_communicator,-1)
C
C
C
CSK      WRITE(LUWRT,*) ' THIS IS MY SUBSPH in P1..._3'
CSK      CALL WRTMATMN(SUBSPH,1,NVEC+IADD-1,1,NVEC+IADD-1,LUWRT)
C
C     zero scratch offsets
C
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
C
C
CSK      WRITE(LUWRT,*) ' THIS IS MY LUIN2LIST in P1..._3'
CSK      CALL IWRTMAMN(LUIN2LIST,1,IALL_LU3,1,IALL_LU3,LUWRT)
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C           set new offset
C
C       position in file is at the end of vector IVEC - 1
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
        IOFFSET_INT_IN2  = 1 + NUM_BLK
C
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                   IOFFSET_IN_LUIN2
C
         CALL RDVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                         LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK         WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)

C
        DO 200 IVEC = 1, NVEC + IADD - 1
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C           set new offset
C
C           position in file is at the end of vector IVEC - 1
C
            IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                       ( IVEC - 1 )  * MY_VEC1_IOFF
            IOFFSET_INT_IN1  = 1 + NUM_BLK   +
     &                       ( IVEC - 1 )  * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN1 in P1..._3 200',
CSK     &                   IOFFSET_IN_LUIN1
C
          CALL RDVEC_BATCH_DRV4(LUIN1,VEC1,LBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                         LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUIN1 in P1..._3 200'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
          FACTOR = SUBSPH(IVEC)
C
C         VEC2 == VEC2 + VEC1 * FACTOR * ( - 1 )
C
          CALL DAXPY(LEBATCH(ISBATCH), -FACTOR, VEC1, 1, VEC2, 1) 
C
  200   CONTINUE
C
        REDSCRVAR = REDSCRVAR 
     &            + DDOT( LEBATCH(ISBATCH), VEC2, 1, VEC2, 1)
C
CSK         WRITE(LUWRT,*) 'final VEC2 to write on LUIN2 == orhtogonalized'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
C       new offset for writing on LUIN2
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
        IOFFSET_INT_IN2  = 1 + NUM_BLK
C
        CALL WTVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                       LUIN2LIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      SCALEVEC = 0.0D0
C
C     communicate REDSCRVAR to get full scale factor
C       
      CAll redvec(REDSCRVAR,SCALEVEC,1,2,op_MPI_SUM,
     &            global_communicator,-1)
C
C     1.4 normalizing the new vector
C
C     zero scratch offsets
C
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
C
C
      FACTOR = 1.0D0 / SQRT( SCALEVEC )
CSK      WRITE(LUWRT,*) 'THIS IS MY SCALING FACTOR',FACTOR
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C
C       position in file is at the end of vector IVEC - 1
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
        IOFFSET_INT_IN2  = 1 + NUM_BLK
C
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                   IOFFSET_IN_LUIN2
C
         CALL RDVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                         LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK         WRITE(LUWRT,*) 'initial VEC2 on LUIN2'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)

C
         CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)

         
         IF (THR_CVEC .GT. 0.0D0 .AND. IROUND .EQ. 1) THEN
           DO I = 1,LEBATCH(ISBATCH) 
             IF(ABS(VEC2(I)) .LT. THR_CVEC) THEN
               IF(VEC2(I) .EQ. 0.0D0) THEN
                 I_ZERO = I_ZERO + 1
               ELSE
                  VEC2(I) = 0.0D0
                  I_REMOVED = I_REMOVED + 1
               END IF
             END IF
           END DO
         END IF
C
C        set new offset
C
C        position in file is at the end of vector NVEC + IADD - 1
C
         IOFFSET_LUOUT  = MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                  (  NVEC + IADD - 1 )  * MY_VEC1_IOFF
C
         IOFFSET_INT_LUOUT = 1 + NUM_BLK +
     &                     ( NVEC + IADD - 1 )  * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUOUT',
CSK     &                   IOFFSET_LUOUT
C
CSK         WRITE(LUWRT,*) 'absolute final new vec on VEC2 to LUOUT'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)
C
         IDEBUGPRNT = 0
C
         CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                         LUOUTLIST,NUM_ACTIVE_BATCH)
C
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
      ! need to do a MPI_ALLREDUCE to check whether any slave has removed an element in the trial vector
      IF(THR_CVEC .GT. 0.0D0 .AND. IROUND .EQ. 1) THEN
         I_REMOVED_total = 0
         I_ZERO_total    = 0
         CAll redvec(I_REMOVED,I_REMOVED_total,1,1,
     &                op_MPI_SUM,global_communicator,-1)
         CAll redvec(I_ZERO,I_ZERO_total,1,1,
     &                op_MPI_SUM,global_communicator,-1)
         IF(I_REMOVED_total .GT. 0) THEN
#ifdef LUCI_DEBUG
            WRITE (luwrt,'(/A,I12,A,I4,A,1P,D10.2,I14)')
     &       'info: Removed',I_REMOVED_total,
     &       ' elements in new CI trial vector no.',IADD,
     &       '; threshold & zeroes',THR_CVEC, I_ZERO_total
#endif
            GO TO 10
         END IF
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE P1_B_PAR_CX_2(VEC1,VEC2,SUBSPH,SCAL_FAC,
     &                         LUIN1LIST,LUIN2LIST,LUOUTLIST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                         MY_IOFF_LUOUT,NVEC,IADD,
     &                         LUIN1,LUIN2,LUOUT)
C
C     Written by  S. Knecht         - June 25 2007
C
C**********************************************************************
C
C     vector orthogonalization and normalization step in MICDVC_PAR ... 
C     
C
C     NOTE: NVEC = NVEC
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*), LUOUTLIST(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN1,IOFFSET_IN_LUIN2
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
      ONEM = - 1.0D0
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN1   = 0
      IOFFSET_INT_IN2   = 0
      IOFFSET_INT_LUOUT = 0
C
CSK      WRITE(LUWRT,*) ' NVEC + IADD - 1',  NVEC + IADD - 1
CSK      WRITE(LUWRT,*) ' LUIN1,LUIN2,LUOUT', LUIN1,LUIN2,LUOUT
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C           set new offset
C
        DO 100 IVEC = 1, ( 2 * (NVEC + IADD - 1) )
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IF( IVEC .eq. 1 ) THEN
C
C           set new offset
C           position in file is at the end of vector IVEC - 1
C
            IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                       ( IVEC - 1 )  * MY_VEC1_IOFF
            IOFFSET_INT_IN1  = 1 + NUM_BLK   +
     &                       ( IVEC - 1 )  * MY_ACT_BLK1
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUIN1 in P1..._3 100'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
            FACTOR = SUBSPH( IVEC )
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
C
          ELSE
C           set new offset
C           position in file is at the end of vector IVEC - 1
C
            IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                       ( IVEC - 1 )  * MY_VEC1_IOFF
            IOFFSET_INT_IN1  = 1 + NUM_BLK   +
     &                       ( IVEC - 1 )  * MY_ACT_BLK1
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK         WRITE(LUWRT,*) 'initial VEC1 on LUIN1 in P1..._3 100'
CSK         CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C
            FACTOR = SUBSPH( IVEC )
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
          END IF
C
  100   CONTINUE
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       =========
C       REAL PART
C       =========
C
C       set new offset
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH  
        IOFFSET_INT_IN2  = 1 + NUM_BLK
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK      WRITE(LUWRT,*) 'initial VEC1 on LUIN2 in P1..._3 100'
CSK      CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)
C
        CALL DSCAL(LEBATCH(ISBATCH),ONEM,VEC2,1) 
C
        CALL DAXPY(LEBATCH(ISBATCH),SCAL_FAC,VEC1,1,VEC2,1)
C      
C       write final vector to disk ...
C
C       position in file is at the end of vector NVEC + IADD - 1
C
        IOFFSET_LUOUT  = MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                 (  2 * (NVEC + IADD - 1) )  * MY_VEC1_IOFF
C
        IOFFSET_INT_LUOUT = 1 + NUM_BLK +
     &                    ( 2 * (NVEC + IADD - 1) )  * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUOUT real',
CSK     &                   IOFFSET_LUOUT
C
CSK         WRITE(LUWRT,*) 'absolute final new vec on VEC2 to LUOUT real'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)
C
        IDEBUGPRNT = 0
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                        LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ end of loop for real part
C
C
      BUF = 0.0D0
      DO JVEC = 1, NVEC+IADD-1
           BUF                   =   SUBSPH( 2*(JVEC-1) + 1 )
        SUBSPH( 2*(JVEC-1) + 1 ) = - SUBSPH( 2*(JVEC-1) + 2 )
        SUBSPH( 2*(JVEC-1) + 2 ) =      BUF
      END DO
C
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN1   = 0
      IOFFSET_INT_IN2   = 0
      IOFFSET_INT_LUOUT = 0
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C           set new offset
C
        DO 200 IVEC = 1, ( 2 * ( NVEC + IADD - 1) )
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IF( IVEC .eq. 1 ) THEN
C
C           set new offset
C           position in file is at the end of vector IVEC - 1
C
            IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                       ( IVEC - 1 )  * MY_VEC1_IOFF
            IOFFSET_INT_IN1  = 1 + NUM_BLK   +
     &                       ( IVEC - 1 )  * MY_ACT_BLK1
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUIN1 in P1..._3 100'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
            FACTOR = SUBSPH( IVEC )
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
C
          ELSE
C           set new offset
C           position in file is at the end of vector IVEC - 1
C
            IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                       ( IVEC - 1 )  * MY_VEC1_IOFF
            IOFFSET_INT_IN1  = 1 + NUM_BLK   +
     &                       ( IVEC - 1 )  * MY_ACT_BLK1
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK         WRITE(LUWRT,*) 'initial VEC1 on LUIN1 in P1..._3 100'
CSK         CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C
            FACTOR = SUBSPH( IVEC )
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
          END IF
C
  200   CONTINUE
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       =========
C       IMAG PART
C       =========
C
C       set new offset
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     MY_VEC1_IOFF    
        IOFFSET_INT_IN2  = 1 + NUM_BLK + MY_ACT_BLK1
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK      WRITE(LUWRT,*) 'initial VEC1 on LUIN2 in P1..._3 100'
CSK      CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)
C
        CALL DSCAL(LEBATCH(ISBATCH),ONEM,VEC2,1)
C
        CALL DAXPY(LEBATCH(ISBATCH),SCAL_FAC,VEC1,1,VEC2,1)
C      
C       write final vector to disk ...
C
C       position in file is at the end of vector NVEC + IADD - 1
C
        IOFFSET_LUOUT  = MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                 (  2 * (NVEC + IADD - 1) )  * MY_VEC1_IOFF +
     &                                               MY_VEC1_IOFF
C
        IOFFSET_INT_LUOUT = 1 + NUM_BLK             + MY_ACT_BLK1 +
     &                    ( 2 * (NVEC + IADD - 1) ) * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUOUT imag',
CSK     &                   IOFFSET_LUOUT
C
CSK         WRITE(LUWRT,*) 'absolute final new vec on VEC2 to LUOUT imag'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                 LUWRT)
C
        IDEBUGPRNT = 0
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                        LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ end of loop for imag part
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE P3_B_PAR_RL_1(VEC1,VEC2,SUBSPH,
     &                         LUINLIST,LUOUTLIST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         MY_IOFF_LUIN,
     &                         MY_IOFF_LUOUT,NVEC,IROOT,LUIN,LUOUT)
C
C     Written by  S. Knecht         - June 9 2007
C
C**********************************************************************
C
C     calculating scaled vecsum between two vectors on file LUIN resp.
C     LUIN; saving on LUOUT
C
C     NOTE: IROOT = IROOT
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUTLIST(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN    = 0
      IOFFSET_INT_LUOUT = 0
C
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 100 IVEC = 1, NVEC
C
          FACTOR = SUBSPH( IVEC )
C
C           set new offset
C
C           position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )   * MY_VEC1_IOFF
          IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                     ( IVEC - 1 )   * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK     &                   IOFFSET_IN_LUIN
C
          IF( IVEC .eq. 1 ) THEN
C
            CALL RDVEC_BATCH_DRV4(LUIN,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                           LUINLIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC2 on LUIN'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
CSK            WRITE(LUWRT,*) 'scaling factor for this vector',FACTOR
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
CSK            WRITE(LUWRT,*) 'this is VEC2 after first scaling in
CP1...2'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
          ELSE
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
            CALL RDVEC_BATCH_DRV4(LUIN,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                           LUINLIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC1 on LUIN'
CSK            CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C           VEC2 == VEC2 + VEC1 * FACTOR
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
CSK            WRITE(LUWRT,*) 'final VEC2 after DAXPY in P1...2'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
          END IF
C         ^ IVEC == 1 ?
C
  100     CONTINUE
C
CSK            WRITE(LUWRT,*) 'final VEC2 to write on position',IROOT -1
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C       new offset for writing on LUOUT
C
        IOFFSET_LUOUT      =  MY_IOFF_LUOUT + IOFFSET_SCRATCH + 
     &                        ( IROOT - 1 ) * MY_VEC1_IOFF
C
        IOFFSET_INT_LUOUT  = 1 + NUM_BLK    +
     &                        ( IROOT - 1 ) * MY_ACT_BLK1
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE P3_B_PAR_CPX(VEC1,VEC2,SUBSPH,SUBSPH_2,
     &                       LUINLIST,LUOUTLIST,
     &                       NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                       MY_IOFF_LUIN,MY_IOFF_LUOUT,NVEC,NVEC_D,
     &                       IROOT,LUIN,LUOUT)
C
C     Written by  S. Knecht         - June 26 2007
C
C**********************************************************************
C
C     calculating scaled vecsum between two vectors on file LUIN resp.
C     LUIN; saving on LUOUT
C
C     NOTE: IROOT  = IROOT
C            NVEC  = NVEC
C           NVEC_D = NVEC * 2
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*), SUBSPH_2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUTLIST(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=df_MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_IN_LUIN
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN    = 0
      IOFFSET_INT_LUOUT = 0
C
C     ==============
C     COMPLEX VECTOR
C     ==============
C
C     ---------
C     REAL PART
C     ---------
C
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 100 IVEC = 1, NVEC_D
C
          FACTOR = SUBSPH( IVEC )
C
C           set new offset
C
C           position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )   * MY_VEC1_IOFF
          IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                     ( IVEC - 1 )   * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK     &                   IOFFSET_IN_LUIN
C
          IF( IVEC .eq. 1 ) THEN
C
            CALL RDVEC_BATCH_DRV4(LUIN,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                           LUINLIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC2 on LUIN'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
CSK            WRITE(LUWRT,*) 'scaling factor for this vector',FACTOR
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
CSK            WRITE(LUWRT,*) 'this is VEC2 after first scaling in
CP1...2'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
          ELSE
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
            CALL RDVEC_BATCH_DRV4(LUIN,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                           LUINLIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC1 on LUIN'
CSK            CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C           VEC2 == VEC2 + VEC1 * FACTOR
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
CSK            WRITE(LUWRT,*) 'final VEC2 after DAXPY in P1...2'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
          END IF
C         ^ IVEC == 1 ?
C
  100     CONTINUE
C
CSK            WRITE(LUWRT,*) 'final VEC2 to write on position',IROOT -1
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C       new offset for writing on LUOUT
C
        IOFFSET_LUOUT      =  MY_IOFF_LUOUT + IOFFSET_SCRATCH + 
     &                        ( IROOT - 1 ) * MY_VEC2_IOFF
C
        IOFFSET_INT_LUOUT  = 1 + NUM_BLK    +
     &                        ( IROOT - 1 ) * MY_ACT_BLK2
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches - real part
C
C
C     zero scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      NUM_ACTIVE_BATCH  = 0
C
      DO JVEC = 1, NVEC
        SUBSPH( (JVEC-1)*2+1 ) = SUBSPH_2( (IROOT-1)* 2 *NVEC+NVEC+JVEC)
        SUBSPH( (JVEC-1)*2+2 ) = SUBSPH_2( (IROOT-1)* 2 *NVEC+JVEC)
      END DO
C
C     ---------
C     IMAG PART
C     ---------
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 200 IVEC = 1, NVEC_D
C
          FACTOR = SUBSPH( IVEC )
C
C           set new offset
C
C           position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )   * MY_VEC1_IOFF
          IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                     ( IVEC - 1 )   * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK     &                   IOFFSET_IN_LUIN
C
          IF( IVEC .eq. 1 ) THEN
C
            CALL RDVEC_BATCH_DRV4(LUIN,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                           LUINLIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC2 on LUIN'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
CSK            WRITE(LUWRT,*) 'scaling factor for this vector',FACTOR
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
CSK            WRITE(LUWRT,*) 'this is VEC2 after first scaling in
CP1...2'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
          ELSE
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
            CALL RDVEC_BATCH_DRV4(LUIN,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                           LUINLIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC1 on LUIN'
CSK            CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C           VEC2 == VEC2 + VEC1 * FACTOR
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
CSK            WRITE(LUWRT,*) 'final VEC2 after DAXPY in P1...2'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
          END IF
C         ^ IVEC == 1 ?
C
  200     CONTINUE
C
CSK            WRITE(LUWRT,*) 'final VEC2 to write on position',IROOT -1
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C       new offset for writing on LUOUT - imag part
C
        IOFFSET_LUOUT      =  MY_IOFF_LUOUT + IOFFSET_SCRATCH + 
     &                        ( IROOT - 1 ) * MY_VEC2_IOFF    + 
     &                                        MY_VEC1_IOFF
C
        IOFFSET_INT_LUOUT  = 1 + NUM_BLK    + MY_ACT_BLK1     +
     &                        ( IROOT - 1 ) * MY_ACT_BLK2
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches - imag part
C
      END
***********************************************************************

      SUBROUTINE SYNC_NODES_REL_P1(ECORE_KRMC,CIRUN,LUWRT_FILE)
C
C     Synchronize MASTER and NODES - general routine
C
C     Written by  S. Knecht         - April 13 2007 (Friday!)
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
#include "mxpdim.inc"
#include "cicisp.inc"
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "maxorb.h"
#include "maxash.h"
#include "dcbidx.h"
#include "krmcluci_inf.h"
      CHARACTER*6 CIRUN
!
!     introduce CIRUN and JKRRUNTYPE for NODES
      call interface_mpi_BCAST(CIRUN,6,MASTER,global_communicator)
      call interface_mpi_BCAST(JKRRUNTYPE,1,MASTER,
     &               global_communicator)
!
      IF( MYPROC .ne. MASTER )THEN
        CALL LUCI_TASK_NODES(CIRUN)
      END IF
 
!     communicate some CIRUN specific variables, arrays ...
      call interface_mpi_BCAST(NZ,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(NFSYM,1,MASTER,
     &               global_communicator)

!     ISYM_KRMC == IOPT_SYMMETRY
      call interface_mpi_BCAST(ISYM_KRMC,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(ECORE_KRMC,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(NASHT,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(NNASHX,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(IASH,NFSYM,MASTER,
     &               global_communicator)
!
!     all needed for proper read-in of integrals...
      call interface_mpi_BCAST(NBSYM,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(JBTOF,8*2,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(IPQTOQ,4*8,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(IQTOPQ,4*8,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(JTRLVL,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(nzxopp,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST_l0(OPT_NOPFQ,1,MASTER,
     &               global_communicator)
!     big integer arrays ...
      call interface_mpi_BCAST(IDXT2G,4*MXCORB,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(IDXU2G,MAXASH,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(IDXG2U,MXCORB,MASTER,
     &               global_communicator)
!
      call interface_mpi_BCAST(NIDX3,2,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(NIDX4,2,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(IIDX3,2,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(IIDX4,2,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(NORBT,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(N2ORBX,1,MASTER,
     &               global_communicator)

      END
***********************************************************************

      SUBROUTINE SYNC_NODES_REL_P2(F1,F2,CIRUN)
C
C     Synchronize MASTER and NODES - CIRUN specific routine
C
C     Written by  S. Knecht         - April 16 2007
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION F1(*), F2(*)
      CHARACTER*6 CIRUN
C
      MZ = MIN(NZ,2)
C
C     communicate arrays F1 and F2
C
C     FIXME: F2 --> distribute integrals but only those 
C                   that are really needed on a CPU
C
C     Fock matrix ...
C
      call interface_mpi_bcast_r1_work_f77(F1,(2*NASHT)*(2*NASHT)*MZ,
     &               MASTER,global_communicator)
C
C     debug printing
C
      IPRNT = 0
      IF( MYPROC .eq. 0 .and. IPRNT .gt. 5 ) THEN
C
        WRITE(LUWRT,*) ' I am printing the FOCK matrix', MYPROC
        WRITE(LUWRT,'(A)') 'Molfdir format FCAC:'
        WRITE(LUWRT,'(A)') 'Real part:'
        CALL WRTMATMN(F1,1,(2*NASHT)*(2*NASHT),1,
     &                (2*NASHT)*(2*NASHT),LUWRT)
      END IF
      IPRNT = 0
C
C     store some important information 
      LF2_ZERO = (2*NASHT)**4 * MZ
C
C     broadcast only for SIGMA and SIGMAD (current fix for large-scale 
C     CIINII runs). FIXME: large-scale KR-MCSCF - SK - july 10, 2007
      IF( CIRUN(1:5) .eq. 'SIGMA' )THEN
C       ... integrals
        call interface_mpi_bcast_r1_work_f77(F2,NASHT*NASHT*NNASHX*NZ*3,
     &                 MASTER,global_communicator)
      END IF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SYNC_NODES_REL_XPROP()
C**********************************************************************
C
C     Synchronize MASTER and NODES in KR-CI property calculation.
C
C     written by S. Knecht - Nov 2008
C
C     last revision:
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krciprop.h"
#include "dcbkrci.h"
#include "dcbxpr.h"
C
C     common block in dcbkrci.h
      call interface_mpi_BCAST(NKRCI_CIROOTS,MAX_NKRCI_MAX_SYM,
     &               MASTER,global_communicator)
C  
C     common block in krciprop.h - integer part
      call interface_mpi_BCAST(NPROP_KRCI,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(LPROP_KRCI,MXPROP_KRCI,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(NPROP_ROOTS_KRCI,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(ISYMOPPRP_KRCI,MXPROP_KRCI,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(IN1ELPRP_KRCI,MXPROP_KRCI,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(ISYMEIG_KRCI,MXPRPE_KRCI,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(IHTYPE_X,1,MASTER,
     &               global_communicator)
C  
C     common block in krciprop.h - logicals
      call interface_mpi_BCAST_l0(DOPROPREOD,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST_l0(DOPROPREOD,1,MASTER,
     &               global_communicator)
C
C     common block in dcbxpr.h - integer part
      call interface_mpi_BCAST(IPRPTIM,MAXPRPS,MASTER,
     &               global_communicator)
C
C     common block in dcbxpr.h - character array
      call interface_mpi_BCAST(PRPNAM,MAXPRPS*16,MASTER,
     &               global_communicator)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SYNC_NODES_BLOCKFILE_DATA(NMPROC_FILE,SPLIT_IJKL_FILE,
     &                                     SHARED_M_FILE,DIST_COMPLETE,
     &                                     ICCTOS,RCCTOS,NPARBLOCK,
     &                                     ICWEIGHTF,IT_TTPL,LEN_TTPL)

C**********************************************************************
C
C     Synchronize MASTER and NODES about block file information.
C
C     written by S. Knecht - Dec 2008
C
C     last revision:
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE), NPARBLOCK(*)
      INTEGER   RCCTOS(*)
#include "parluci.h"
      DIMENSION ICCTOS(I_NZERO_LEN_C,I_NZERO_LEN_S), IT_TTPL(LEN_TTPL)
      LOGICAL SHARED_M_FILE, SPLIT_IJKL_FILE, DIST_COMPLETE
C
      call interface_mpi_BCAST(NMPROC_FILE,1,
     &               MASTER,global_communicator)
C
      IF( NMPROC_FILE .ne. NMPROC)THEN
        call interface_mpi_bcast_i2_work_f77(ICCTOS,
     &                 I_NZERO_LEN_C*I_NZERO_LEN_S,
     &                 MASTER,global_communicator)
        call interface_mpi_bcast_i1_work_f77(ICWEIGHTF,
     &                 I_NZERO_LEN_S,MASTER,
     &                 global_communicator)
      END IF
C
C     reduced connection list and block distribution list
      call interface_mpi_bcast_i1_work_f77(RCCTOS,NUM_BLOCKS2,MASTER,
     &               global_communicator)
      call interface_mpi_bcast_i1_work_f77(NPARBLOCK,NUM_BLOCKS,MASTER,
     &               global_communicator)
C
C     T list?
      call interface_mpi_BCAST_l0(SHARED_M_FILE,1,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST_l0(SPLIT_IJKL_FILE,1,MASTER,
     &               global_communicator)
      IF( SHARED_M_FILE .or. SPLIT_IJKL_FILE )THEN
        call interface_mpi_BCAST(LEN_TTPL,1,MASTER,
     &                 global_communicator)
        call interface_mpi_BCAST(IT_TTPL,LEN_TTPL,MASTER,
     &                 global_communicator)
      END IF
      call interface_mpi_BCAST_l0(DIST_COMPLETE,1,MASTER,
     &               global_communicator)
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE READ_INT_TTYPES(T_BUFF,IT_TTPL,NDIM,N1_INT,IOFF_T)
C
C     read in integrals to T_BUFF according to T-type list IT_TTPL
C
C     NDIM    : number of T-types
C     IOFF_T  : file offset
C
C     Written by  S. Knecht         - 15 February 2008             
C
C     Last revision: 
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
#include "typesz_mpi2.h"
#include "ctcc.inc"
      DIMENSION T_BUFF(*), IT_TTPL(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFF_T, IOFF_SCR
#if defined (INT_STAR8)
      INTEGER(KIND=df_MPI_OFFSET_KIND)  IREAD_BYTE
#else
      INTEGER IREAD_BYTE
#endif
      INTEGER(KIND=df_MPI_OFFSET_KIND) IS_BYTE8, ISLEN_TBL, NINT_8, N1_8
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFF_T_B, INUM_BYTE, nelm_cc
C
      N1_8       = N1_INT
      NINT_8     = N1_8 + N2ELINT
      IOFF_SCR   = 0
      IREAD_BYTE = 0
      INUM_BYTE  = 0
      IS_BYTE8   = ISIZE_dp
      ISLEN_TBL  = 0
      IOFF_SCR   = IOFF_T 
      IOFF_T_B   = 1
csk   WRITE(LUWRT,*) ' IOFF_SCR, IS_BYTE8, NINT_8, NDIM, N1_INT',
csk  &                 IOFF_SCR, IS_BYTE8, NINT_8, NDIM, N1_INT
csk   WRITE(LUWRT,*) ' showing IT_TTPL with IRC_SAVE',IRC_SAVE
csk   CALL IWRTMAMN(IT_TTPL,1,NDIM*IRC_SAVE,1,NDIM*IRC_SAVE,LUWRT)
      IF( SPLIT_IJKL )THEN
        DO IMULT = 1, IRC_SAVE
          DO IJKL = 1, NDIM
C
            ISLEN_TBL = NELM_CC(IJKL,NDIM,NINT_8)
csk       WRITE(LUWRT,*) ' ISLEN_TBL',ISLEN_TBL
C
            ITYPE = IJKL + NDIM * ( IMULT - 1)
            IF( IT_TTPL( ITYPE ) .gt. 0 )THEN
#if defined (INT_STAR8)
              IREAD_BYTE = IS_BYTE8 * ISLEN_TBL
#else
              XLENTIT      = REAL(ISLEN_TBL)
              ISLEN_TBL_I4 = IGIVE_I_B(XLENTIT) 
              IREAD_BYTE   = ISIZE_dp * ISLEN_TBL_I4
#endif
csk         WRITE(LUWRT,*) ' reading at IOFF_SCR: # elements (byte):',
csk  &                                        IOFF_SCR,   IREAD_BYTE
csk         WRITE(LUWRT,*) ' memory offset: IOFF_T_B',IOFF_T_B
              call interface_mpi_FILE_READ_AT_br(IIJKL_ROD,IOFF_SCR,
     &                              T_BUFF(IOFF_T_B),
     &                              IREAD_BYTE,df_MPI_BYTE,ISTAT)
              IOFF_T_B = IOFF_T_B + ISLEN_TBL
            ENDIF
            INUM_BYTE  = IS_BYTE8 * ISLEN_TBL
            IOFF_SCR   = IOFF_SCR + INUM_BYTE
C
          END DO
        END DO
      ELSE
        DO IMULT = 1, IRC_SAVE
          DO IJKL = 1, NDIM
C
            ISLEN_TBL = NELM_CC(IJKL,NDIM,NINT_8)
csk       WRITE(LUWRT,*) ' ISLEN_TBL',ISLEN_TBL
C
            ITYPE = IJKL + NDIM * ( IMULT - 1)
            IF( IT_TTPL( ITYPE ) .eq. MYNEW_ID_SM )THEN
#if defined (INT_STAR8)
              IREAD_BYTE = IS_BYTE8 * ISLEN_TBL
#else
              XLENTIT      = REAL(ISLEN_TBL)
              ISLEN_TBL_I4 = IGIVE_I_B(XLENTIT) 
              IREAD_BYTE   = ISIZE_dp * ISLEN_TBL_I4
#endif
csk         WRITE(LUWRT,*) ' reading at IOFF_SCR: # elements (byte):',
csk  &                                        IOFF_SCR,   IREAD_BYTE
csk         WRITE(LUWRT,*) ' memory offset: IOFF_T_B',IOFF_T_B
              call interface_mpi_FILE_READ_AT_br(IIJKL_ROD,IOFF_SCR,
     &                              T_BUFF(IOFF_T_B),
     &                              IREAD_BYTE,df_MPI_BYTE,ISTAT)
              IOFF_T_B = IOFF_T_B + ISLEN_TBL
            ENDIF
            INUM_BYTE  = IS_BYTE8 * ISLEN_TBL
            IOFF_SCR   = IOFF_SCR + INUM_BYTE
C
          END DO
        END DO
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE TEST_PAR_CALC(PARCAL_SAVE,ISET)
C
C     Test if DIRAC is running in parallel
C
C     ISET = 0: save  PARCAL (dcbgen.h) in   value PARCAL_SAVE
C     ISET = 1: reset PARCAL (dcbgen.h) from value PARCAL_SAVE
C
C     Written by  S. Knecht         - June 16 2007  - on a saturday
C     afternoon ...
C
C**********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "dcbgen.h"
      LOGICAL PARCAL_SAVE
C
      IF( ISET .eq. 0 )THEN
        PARCAL_SAVE = PARCAL
        PARCAL = .FALSE.
      ELSE
        PARCAL = PARCAL_SAVE
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE UPDATE_LUC_LIST(ISCLFAC_GROUP,LUCLIST,RCCTOS,CB,
     &                           NPARBLOCK,IBLOCKL,IGROUPLIST,
     &                           IPROCLIST,IRILP,BLOCKTIME)
C
C     make an update of of grouplist for c-vector file based on 
C     different list gathered from global_communicator
C
C
C     Written by  S. Knecht         - July 04 2007 
C
C     OUTPUT: ISCLFAC_GROUP and updated file ILUC
C
C**********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
#include "wrkspc.inc"
#include "ipoist8.inc"
      DIMENSION ISCLFAC_GROUP(*), LUCLIST(*), CB(*)
      INTEGER   RCCTOS(*)
      DIMENSION NPARBLOCK(*), IBLOCKL(*)
      DIMENSION IGROUPLIST(*), IPROCLIST(*)
      CHARACTER*12 WALLTID3, SECTID
C
C     set mark for local memory
      IDUM = 0
      CALL MEMMAR(KDUM,  IDUM,    'MARK  ',IDUM,'UPLIST')
C
      CALL MEMMAR(KSCALLOC2,NUM_BLOCKS2,'ADDL  ',1,'ICLLC2')
      CALL MEMMAR(KSCALLOC3,NUM_BLOCKS2,'ADDL  ',1,'ICLLC3')
C
C     fill complete local iscalfac arrays with zero
      CALL IZERO(WORK(KSCALLOC2), NUM_BLOCKS2)
      CALL IZERO(WORK(KSCALLOC3), NUM_BLOCKS2)
      CALL IZERO(ISCLFAC_GROUP  , NUM_BLOCKS2)
C
      starttimer = interface_MPI_WTIME()
C
C     "mpi_allsum" local LUCLIST which then on all
C     nodes will contain the number of non-zero C-blocks in
C     the complete CI-vector
      CAll redvec(LUCLIST,WORK(KSCALLOC2),NUM_BLOCKS2,1,
     &                op_MPI_SUM,global_communicator,-1)
C
C     find all c-blocks connecting to all sigma-blocks on each cpu
      CALL ICOPY(NUM_BLOCKS2,RCCTOS,1,WORK(KSCALLOC3),1)

!     WRITE(LUWRT,*) '  subroutine UPDATE_LUC_LIST speaking'
!     WRITE(LUWRT,*) 'RCCTOS:'
!     CALL IWRTMAMN(RCCTOS,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUWRT)

C
C     case 1: number of CPUs in new group != to total number of CPUs
C     case 2: number of CPUs in new group == to total number of CPUs
C
      IF( NEWCOMM_PROC .ne. NMPROC ) THEN
C
        CAll redvec(WORK(KSCALLOC3),ISCLFAC_GROUP,NUM_BLOCKS2,1,
     &                  op_MPI_SUM,MYNEW_COMM,0)
C
C       all local node-masters call this routine!
C
        IF( MYNEW_ID .eq. 0 ) THEN
           CALL COPVCD_PAR_BDRIV5_REL(ILUC,ILUC,CB,NPARBLOCK,
     &                                WORK(KSCALLOC2),ISCLFAC_GROUP,
     &                                IBLOCKL,NUM_BLOCKS,ICOMM,
     &                                IGROUPLIST,IPROCLIST,IRILP)
C          CALL COPVCD_PAR_BDRIV5_REL(LUIN,LUOUT,SEGMNT,IBLOCKD,
C     &                               ISCALFAC,ISCALFAC_GROUP,
C     &                               IBLOCKL,NBLOCK,JCOMM,
C     &                               IGROUPLIST,IPROCLIST,IRILP)
C

        END IF
        call interface_mpi_bcast_i1_work_f77(ISCLFAC_GROUP,NUM_BLOCKS2,
     &                 0,MYNEW_COMM)
      ELSE
C
       CALL UPDATE_GEN_LIST(WORK(KSCALLOC3),WORK(KSCALLOC2),NUM_BLOCKS2)
C
C      to be consistent with output of case 1
       CALL IZERO(ISCLFAC_GROUP,NUM_BLOCKS2)
       CALL ICOPY(NUM_BLOCKS2,WORK(KSCALLOC3),1,ISCLFAC_GROUP,1)
C
      END IF
C     ^ NEWCOMM_PROC == NMPROC ?
C
!     WRITE(LUWRT,*) 'LUCLIST:'
!     CALL IWRTMAMN(LUCLIST,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUWRT)
!     WRITE(LUWRT,*) 'ISCLFAC_GROUP:'
!     CALL IWRTMAMN(ISCLFAC_GROUP,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUWRT)
C
C     final timing for block distribution
      blocktime = blocktime + interface_MPI_WTIME() - starttimer
C
C     flush local memory
C
      CALL MEMMAR(KDUM ,IDUM,'FLUSM ',2,'UPLIST')
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE Z_BLKFO_REL_PAR(IDC,NMS2VAL,ISPC,ISM,KPCLBT,KPCLEBT,
     &                           KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK,
     &                           NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2,
     &                           I_SET_L2BLOCK,I_USE_PC,NPARBLOCK)
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM.
*
* Output is given in the form of pointers to vectors in WORK
* and  NBLK_MS2,IBLK_MS2 which should be dimensioned outside
*
* KPCLBT : Length of each Batch ( in blocks)
* KPCLEBT : Length of each Batch ( in elements)
* KPCI1BT : Length of each block
* KPCIBT  : Info on each block
* KPCBLTP : Block type for each symmetry
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* NBLK_MS2 : Number of blocks with a given MS2 values
* IBLK_MS2 : First block with a given MS2 value
*
* Jeppe Olsen, Feb. 98
*
* For relativistic program
      use symmetry_setup_krci
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
* new
#include "parluci.h"
*
*.output
      INTEGER NBLK_MS2(*), IBLK_MS2(*)
      INTEGER NBAT_MS2(*), IBAT_MS2(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' Output from Z_BLKFO '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        WRITE(6,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
*. Pointers to output arrays
      NTTS = MXNTTS
      call memmar(KPCLBT ,  MXNTTS,'ADDL  ',1,'CLBTPA')
      call memmar(KPCLEBT,  MXNTTS,'ADDL  ',1,'CLEBPA')
      call memmar(KPCI1BT,  MXNTTS,'ADDL  ',1,'CI1BPA')
      call memmar(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBTPA')
      call memmar(KPCBLTP,  NSMST, 'ADDL  ',2,'CBLTPA')
*.    ^ These should be preserved after exit so put mark for
*       flushing after
*
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'Z_BLp1')
*
      IDUM = 0
*
*     determine allowed length of each batch
*   
      IF(ICISTR.EQ.1) THEN
*       assumes determinants
        L0BLOCK = XISPSM(ISM,ISPC)
      ELSE IF (ICISTR.EQ.2) THEN
        L0BLOCK = MXSB
      ELSE IF( ICISTR.EQ.3) THEN
        L0BLOCK = MXSOOB
      END IF
*
      IF( I_SET_L2BLOCK .eq. 1 ) THEN
*
*       L2BLOCK = max memory for c-vec and sigma-vec from
*       ( curr. total mem - 3 mio for string info etc.) / 4
*
        L2BLOCK = 0
*
        CALL MEMMAR(L2BLOCK,0,'SFREEM',2,'SEEFRM')
        L2BLOCK        = LMEMFREE_PTR
CSK        WRITE(LUWRT,'(1X,A,1X,I20)')
CSK     & '  Current available free memory in double words:',L2BLOCK
*
*       we want to keep three blocks in memory at the same time
*       CB,SB,VEC3(=C2). estimated scratch memory: 3 000 000 real*8
*       division by a factor of 4 = safety!
*
        L2BLOCK = (L2BLOCK - 5 000 000 )/4
*
CSK        WRITE(LUWRT,*) '  L0BLOCK,L2BLOCK,LCSBLK ',
CSK     &                    L0BLOCK,L2BLOCK,LCSBLK
*
        L2BLOCK = MIN(LCSBLK,L2BLOCK)
*
      ELSE IF (I_SET_L2BLOCK .eq. -1 ) THEN
        L2BLOCK = 100 000
      END IF
*    /\ I_SET_L2BLOCK
*
*     set LBLOCK value
*
      LBLOCK  = MAX(L0BLOCK,L2BLOCK)
*
CACTU      LBLOCK  = MAX(L0BLOCK,L2BLOCK)
CTEST      LBLOCK  = MIN(L0BLOCK,L2BLOCK)
*          old versions works without L0BLOCK and L2BLOCK
CVOLD     LBLOCK = MAX(LBLOCK,LCSBLK)
*
      IF(NTEST.GE.1000) THEN
         WRITE(LUWRT,*) ' LBLOCK = ', LBLOCK
      END IF
      IRUN_THIS = 0
      ILENB_OLD = 0
      INUMB_OLD = 0
*     loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
*
        IRUN_THIS = IRUN_THIS + 1
*
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
*
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
*
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*       info needed for generation of block info
        call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
*
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,WORK(KPCBLTP))
*
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK = NBLOCK + 1
Cold          IOFFBTC = NBATCH + 1
          IOFFBTC   = NBATCH
        END IF
*. Batches of C vector
        ITTSS_ORD = 2
*
*         special routine that uses only node-blocks to determine
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
CSK       replace LBLOCK by L0BLOCK for test purposes --> leads to 
CSK       creation of many many batches
*
          CALL PART_CIV4_PAR(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                     WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,LBLOCK,WORK(KLCIOIO),
     &                     ISMOST(1,ISM),
     &                     NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                     WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                     NPARBLOCK,IDEBUGPRNT,
     &                     ILENB_OLD,INUMB_OLD,IRUN_THIS)
*
*       number of BLOCKS hitherto
        NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &         + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
*
      END DO
*     ^ End of loop over MS2 values

      NTEST = 00
      IF(NTEST.GE.1) THEN
        WRITE(LUWRT,*) ' Total number of batches', NBATCH
        WRITE(LUWRT,*) ' Total number of blocks ', NBLOCK
        IF( NTEST .ge. 10 ) THEN
          CALL INFO_PRINT_BATCH_REL(WORK(KPCLBT), WORK(KPCLEBT),
     &                              WORK(KPCI1BT),WORK(KPCIBT),NBATCH,
     &                              LUWRT)
        END IF
      END IF
      NTEST = 0
*
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'Z_BLp1')
*
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE Z_BLKFO_XPRP_PAR(IDC,NMS2VAL,ISPC,ISM,LBATX,LEBATX,
     &                            I1BATX,IBATX,NBATX,NPARBLOCK_X,IXCT)
**********************************************************************
C
C     Construct information about batch and block structure of CI space
C     defined by ISPC and ISM.
C
C     LBATX : length of each batch (in blocks)
C     LEBATX: length of each batch (in elements)
C     I1BATX: starting block of each batch
C     IBATX : info on each block
C     NBATX : number of batches
C
C      Jeppe Olsen, Feb 1998
C
C      S. Knecht - Nov 2008: modified for CI property runs
C
C**********************************************************************
      use symmetry_setup_krci
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
#include "parluci.h"
      DIMENSION NPARBLOCK_X(NUM_BLOCKS,*)
C
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'MARK  ',IDUM,'XZBLKF')
C
C     get batches and blocks
      IRUN_THIS = 0
      ILENB_OLD = 0
      INUMB_OLD = 0
      NTTS      = MXNTTS
      CALL MEMMAR(KPCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
C
C     loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
*
        IRUN_THIS = IRUN_THIS + 1
        IATP      = IST_FOR_DT(1,IMS2)
        IBTP      = IST_FOR_DT(2,IMS2)
        IOCTPA    = IBSPGPFTP(IATP)
        IOCTPB    = IBSPGPFTP(IBTP)
        NOCTPA    = NOCTYP(IATP)
        NOCTPB    = NOCTYP(IBTP)
C       info needed for generation of block info
        call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
C
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,WORK(KPCBLTP))
C
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK   = NBLOCK + 1
          IOFFBTC   = NBATX
        END IF
C
        ITTSS_ORD = 2
C
C       special routine that uses only node-blocks to determine
C       batches of a CI-vector
        CALL PART_CIV4_PAR(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                     WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,NSMST,LBLOCK,WORK(KLCIOIO),
     &                     ISMOST(1,ISM),NBATX,
     &                     LBATX,LEBATX,I1BATX,IBATX,0,ITTSS_ORD,
     &                     NPARBLOCK_X(1,IXCT),0,ILENB_OLD,
     &                     INUMB_OLD,IRUN_THIS)
C
C       number of BLOCKS
        NBLOCK = IFRMR_REL(I1BATX,1,NBATX)+IFRMR_REL(LBATX,1,NBATX)-1
      END DO
C     ^ End of loop over MS2 values
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'FLUSM ',IDUM,'XZBLKF')
      END
#else
      SUBROUTINE PAR_LUCIAREL_DUMMY
      END
#endif
