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
      subroutine picasso(IPRTRA,MCSCF,NL2D,IUB2D,ISM2D,F1,F2,T_BUFF,
     &                   T_BUFF_D,IT_TTPL,IBOSYM,CIRUN,
     &                   RUN_WO_T_BUFF)
****************************************************************
***                     PICASSO                              ***
***                                                          ***
***   P rogram                                               ***
***     for                                                  ***
***   I ntelligent                                           ***
***   C omplete                                              ***
***   A ctive                                                ***
***   S pace                                                 ***
***   S pin                                                  ***
***   O rbit                                                 ***
***     calculations                                         ***
***         (at least what it was once planned to become...) ***
***                                                          ***
****************************************************************
***
*** Written by Timo Fleig, Bonn
***        and Jeppe Olsen, Lund
***
*** October 1996 - February 1998
!
!     general clean-up: all of what is left for picasso is being 
!                       an interface to integral_distribution_driver...
!     S. Knecht - December 2010
!
      implicit real*8 (A-H,O-Z)
#include "mxpdim.inc"
#include "crun.inc"
      logical mcscf, run_wo_t_buff
      character*6 CIRUN
 
      if(.not.run_wo_t_buff)then
        call integral_distribution_driver(ISPINFREE,NL2D,IUB2D,ISM2D,
     &                                    ISYM_T,F1,F2,T_BUFF,T_BUFF_D,
     &                                    IT_TTPL,IBOSYM,MCSCF,iprtra,
     &                                    CIRUN)
      end if
      end

***********************************************************************
      subroutine integral_distribution_driver(ISPINFREE,NL2D,IUB2D,
     &                                        ISM2D,ISYM_T,F1,F2,T_BUFF,
     &                                        T_BUFF_D,IT_TTPL,
     &                                        IBOSYM,MCSCF,iprtra,CIRUN)
*
*---------------------------------------------------
*  Written by Timo Fleig (including all subroutines)
*---------------------------------------------------
*
*  Implementation of MPI "shared" memory version 
*  (use of T_BUFF) by Stefan Knecht  
*
*
      use interface_to_mpi
      use symmetry_setup_krci
      use memory_allocator
#include "implicit.inc"
#include "ipoist8.inc"
#include "infpar.h"
#include "parluci.h"
#include "krciprop.h"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "clunit.inc"
#include "oper.inc"
#include "ctcc.inc"
#include "cgas.inc"
#include "integrals_off.inc"
*
#include "dgroup.h"
#include "dcborb.h"
*
      character*6 CIRUN
      dimension T_BUFF(*), T_BUFF_D(*), IT_TTPL(*)
      integer              :: ndim_tot
      integer, allocatable :: ibeig(:)
      integer, allocatable :: icheck_file(:)
      logical mcscf, ex, alloc_int, alloc_int_scr, fileopen, ex_mc
      logical int_file_is_master_file
      integer, parameter   :: integral_comm_limit = 250000000 !approx.  2GB limit
      integer              :: batch_ints_loops
      integer              :: remaining_batch_ints
      integer(8)           :: kloop_pointer
*
C
      INTIMP = 6
      if(MCSCF) INTIMP = 7
*
      CALL QENTER('idistdrv')
*
      NTESTL = 00
      NTEST = max(iprtra/100,NTESTL)
!     NTEST = 000
*
      IF( NTEST .ge. 5 ) THEN
      WRITE(LUWRT,*) ' ************************************************'
      WRITE(LUWRT,*) ' *                                              *'
      WRITE(LUWRT,*) ' *  Transformer man                             *'
      WRITE(LUWRT,*) ' *  You run the show                            *'
      WRITE(LUWRT,*) ' *  Remote control                              *'
      WRITE(LUWRT,*) ' *  Direct the action with the push of a button *'
      WRITE(LUWRT,*) ' *                                              *'
      WRITE(LUWRT,*) ' *                               Neil Young     *'
      WRITE(LUWRT,*) ' *                                              *'
      WRITE(LUWRT,*) ' ************************************************'
      END IF
*
!     set total number of spin-orbitals/Kramers-pairs
      ndim_tot = 0
      do i = 1, ngas
        ndim_tot = ndim_tot + ngsht(i)
      end do
*
*     1. allocate static memory
!
!     calculate # of diagonal integrals
      IRIMLT = MIN(NZ,2)
      NDIAMINT = (2*NDIM_TOT + 3*((2*NDIM_TOT)**2)) * IRIMLT/2
!
      IF( NDIAMINT .gt. MXINT_DIA ) THEN
        WRITE(LUWRT,*) 
     &  ' STOP in subroutine integral_distribution_driver !!!'
        WRITE(LUWRT,*) ' Number of diagonal integrals larger'
        WRITE(LUWRT,*) ' than MXINT_DIA: ', NDIAMINT, MXINT_DIA
        WRITE(LUWRT,*) ' increase MXINT_DIA (mxpdim.inc): to', 
     &                   NDIAMINT
        CALL QUIT(
     &  ' *** Error in integral_distribution_driver *** value'//
     &  ' of MXINT_DIA too low.')
      END IF
C
C     ... property run. This subroutine is called twice. 
C     1st time: do reordering - no slaves enter this routine
C     2nd time: all CPUs enter this routine - calculate 
C               necessary information and allocate required arrays. 
C               Reordered property matrix is communicated outside 
C               this routine.
      IF( CIRUN .eq. 'PROP1 ')THEN
         ALLOC_INT = .TRUE.
         IF( .NOT. DOPROPREOD ) ALLOC_INT = .FALSE.
         GOTO 100
      END IF
C     ... other run types
      ALLOC_INT = .TRUE.
C
C     allocate integral/density array (if necessary)
C
C     IT_SHL == - 2    :  if MPI-shared memory AND CIRUN == .IJKLRO
C     IT_SHL == - 1    :  if not MPI-shared memory
C     IT_SHL == 0 or 1 :  if MPI-shared memory
C
      IF( IT_SHL .eq. - 2 )THEN
C
C       what may happen if we do a little DIRTY trick?
C       simulating local disks but having a global disk!!!
C       will it crash??? SK - Feb08
C
        IF( NFLGRPS_REL .eq. 1 ) THEN
          IF( MYPROC .ne. MASTER ) ALLOC_INT = .FALSE.
        ELSE
          IF( MYNEW_ID .ne. N_MASTER ) ALLOC_INT = .FALSE.
        END IF
C
      ELSE IF( IT_SHL .eq. - 1 )THEN
C
        IF( CIRUN .eq. 'IJKLRO' )THEN
          IF( NFLGRPS_REL .eq. 1 ) THEN
            IF( MYPROC .ne. MASTER ) ALLOC_INT = .FALSE.
          ELSE
            IF( MYNEW_ID .ne. N_MASTER ) ALLOC_INT = .FALSE.
          END IF
        ELSE
          IF( SPLIT_IJKL ) ALLOC_INT = .FALSE.
        END IF
      ELSE
        ALLOC_INT = .FALSE.
      END IF
C
      EX = .FALSE.
      IF( CIRUN .eq. 'IJKLRO' )THEN
        IF( ALLOC_INT )THEN
          INQUIRE(FILE='IJKL_REOD',EXIST=EX)
C
          IF( EX ) THEN
            WRITE(LUWRT,*)' '
            WRITE(LUWRT,*)' using provided integral file '
            WRITE(LUWRT,*)' IJKL_REOD! '
            WRITE(LUWRT,*)' '
          END IF
          IXXX_CLR = 1
        ELSE  
          IXXX_CLR = 2
        END IF
CSK     WRITE(6,*) ' my EX result: ALLOC_INT',myproc,ex,ALLOC_INT
#if defined(VAR_MPI2)
        IXXX_KEY  = MYPROC
        IXCOMM_SZ = 0
        IXCOMM_RK = 0
        IXCOMM    = 0
        call interface_mpi_COMM_SPLIT(global_communicator,IXXX_CLR,
     &                       IXXX_KEY,IXCOMM)
        call interface_mpi_COMM_SIZE(IXCOMM,IXCOMM_SZ)
        call interface_mpi_COMM_RANK(IXCOMM,IXCOMM_RK)
        IF(ALLOC_INT)THEN

          call alloc(icheck_file,ixcomm_sz)
          call izero(icheck_file,ixcomm_sz)
          IEX_GAT = 0 
          IF(EX) IEX_GAT = 1
          call interface_MPI_GATHER(iex_gat,1,icheck_file,1,0,
     &                              ixcomm)
#ifdef LUCI_DEBUG
          NPTESTVAR = 0
          IF( NPTESTVAR .ge. 20 )THEN
            call interface_mpi_BARRIER(IXCOMM)
            DO IJK = 1, IXCOMM_SZ
              call interface_mpi_BARRIER(IXCOMM)
              IF( IXCOMM_RK .eq. IJK - 1 )THEN
                WRITE(6,*) ' my icheck_file array', 
     &                       MYPROC, MYNEW_ID, IXCOMM_RK
                DO JJJ = 1, IXCOMM_SZ
                   WRITE(6,*) ' icheck_file(',jjj,') = ',
     &                          icheck_file(jjj)
                END DO
              END IF
            END DO
          END IF
#endif
          if(IXCOMM_RK .eq. 0)then
            call check_consist(icheck_file,IXCOMM_SZ)
          end if
          call dealloc(icheck_file)
        END IF ! alloc_int
C       consistency is checked; update global_communicator
        call interface_mpi_BCAST_l0(EX,1,MASTER,
     &       global_communicator)
        call interface_mpi_COMM_FREE(IXCOMM)
#endif
      END IF
C     CIRUN == IJKLRO
      IF( EX ) THEN
        GOTO 1011
      END IF
C
 100  CONTINUE
C
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
      IF( CIRUN .eq. 'IJKLRO' ) THEN
        ALLOC_INT_SCR = ALLOC_INT
        ALLOC_INT     = .FALSE.
      END IF
#endif
      IF( CIRUN .eq. 'IJKLRO' )THEN
!       using low-memory sorting scheme (disk-I/O driven)
        IF(LOWSRT_IJKL) THEN 
          CALL MEMMAR(KT_CC,           0,'ADDS  ',2,'KT_CCX')
          CALL MEMMAR(KMDIASV,4*NDIAMINT,'ADDS  ',1,'MDIAXX')
          CALL IZERO(WORK(KMDIASV),4*NDIAMINT)
          CALL IZERO(IPOSDIA,NDIAMINT)
          GOTO 200
        END IF
      END IF
!
!     allocate array for diagonal integrals
      IF( ALLOC_INT )THEN
        call memmar(KMDIASV,4*NDIAMINT,'ADDS  ',1,'MDIASV')
        CALL IZERO(IPOSDIA,NDIAMINT)
        CALL IZERO(WORK(KMDIASV),4*NDIAMINT)
        IF(CIRUN .ne. 'PROP1 ' ) WRITE(LUWRT,'(A,I16,A/)') 
     &       '  Allocating for ',NDIAMINT,' diagonal integrals.'
      ELSE
        call memmar(KMDIASV, 0,'ADDS  ',1,'MDIASV')
        CALL IZERO(IPOSDIA,NDIAMINT)
      END IF
!
!     allocate full integral array
      IF( ALLOC_INT ) THEN
        IF(CIRUN .ne. 'PROP1 ' ) WRITE(LUWRT,'(/A,I16,A)') 
     &       '  Allocating for ',NALLINT,' integrals.'
        CALL MEMMAR_I8(KT_CC,NALLINT,'ADDS  ',2,'KT_CC ')
        CALL DZERO8(WORK(KT_CC),NALLINT)
      ELSE
        CALL MEMMAR_I8(KT_CC, 0,'ADDS  ',2,'KT_CC ')
      END IF
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
        IF( CIRUN .eq. 'IJKLRO' ) ALLOC_INT = ALLOC_INT_SCR
#endif
 200  CONTINUE
*
*
************************************************
* Set mark for flushing local memory
************************************************
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'ITRCTL')
*
      if (NTEST.ge.7) then
       write(LUWRT,*)
       write(LUWRT,*) '********************************************'
       write(LUWRT,*) '* Distributing MOLTRA transformed          *'
       write(LUWRT,*) '* dirac one- (and two-) electron integrals *'
       write(LUWRT,*) '* using input integral reordering          *'
       write(LUWRT,*) '********************************************'
      end if
*
      IRIOFF = (2*NDIM_TOT)**2
      NDIM_TOT2 = NDIM_TOT*2
C
C     distribute 1-/2- e- integrals in SIGDEN formalism
!     allocate local memory
      call memmar(KISTR_CA,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_CA')
      call memmar(KISTR_CB,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_CB')
      call memmar(KISTR_AA,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_AA')
      call memmar(KISTR_AB,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_AB')
C
C     orbital symmetry info
      call alloc(ibeig,norbt)
      call izero(ibeig,norbt)
#if defined (VAR_MPI)
      if(CIRUN.ne.'PROP1 ')then
        call interface_mpi_bcast_l0(linear,1,master,
     &                 global_communicator)
      end if
#endif
      if(linear)then
        if(myproc.eq.master)then
          luni = -1
          inquire(file='KRMCSCF',exist=ex_mc,opened=fileopen,
     &            number=luni)
          if(ex_mc.and.fileopen)then
            rewind(luni)
          else if(ex_mc.and.(.not.fileopen))then
            luni = 99
            call opnfil(luni,'KRMCSCF','OLD    ','itrctl')
            rewind(luni)
          end if
          call ireakrmc(luni,'MJVEC   ',ibeig,norbt)
          if(.not.fileopen) close(luni,status='keep')
        end if
#if defined (VAR_MPI)
        if(CIRUN.ne.'PROP1 ')then
          call interface_mpi_bcast(ibeig,norbt,master,
     &                   global_communicator)
        end if
#endif
      end if
        
 
C
      IF( REORD_IJKL ) THEN
        IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' ) THEN
#if defined (VAR_MPI2)
C
            WRITE(LUWRT,'(/A)') 
     &      '  ===================================='
            WRITE(LUWRT,'(A)') 
     &      '  Reading from integral file IJKL_REOD'
            WRITE(LUWRT,'(A/)') 
     &      '  ===================================='
#endif
C
          CALL READ_IJKL_GAS(WORK(KT_CC),T_BUFF,T_BUFF_D, 
     &                       IT_TTPL,WORK(KMDIASV),NDIAMINT,
     &                       NSPOBEX_TP,N1ELINT)
C
          call dealloc(ibeig)
          GOTO 1001
C
        ELSE IF( CIRUN .eq. 'IJKLRO' ) THEN
C
#if defined (VAR_MPI2)
C
          IF( ALLOC_INT )THEN
            IXXX_CLR = 1
          ELSE  
            IXXX_CLR = 2
          END IF
          IXXX_KEY = MYPROC
          call interface_mpi_COMM_SPLIT(global_communicator,IXXX_CLR,
     &                        IXXX_KEY,IXCOMM)
C
          IF( IXXX_CLR .eq. 2 ) GOTO 797
#endif
C
          WRITE(LUWRT,*)' '
          WRITE(LUWRT,*)' Resorting integrals...'
          WRITE(LUWRT,*)' '
C
          CALL WRITE_IJKL_GAS( WORK(KT_CC), T_BUFF,
     &                         INTIMP,ISPINFREE,ISYM_T,
     &                         WORK(KLSOBEX),WORK(KISTR_CA),
     &                         WORK(KISTR_CB),WORK(KISTR_AA),
     &                         WORK(KISTR_AB),NDIM_TOT2,
     &                         F1,WORK(KMDIASV),IBOSYM,
     &                         WORK(KLABEXTP),CIRUN,NTEST,
     &                         NDIAMINT,'DIRAC',ibeig,'SIGDEN')
C
 797      CONTINUE
C
#if defined (VAR_MPI2)
          call interface_mpi_COMM_FREE(IXCOMM)
#endif
C
          call dealloc(ibeig)
          GOTO 1001
        END IF ! CIRUN ???
      END IF ! REORD_IJKL
C
      IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' )THEN
C
        IF( IT_SHL .ge. 0 .or. SPLIT_IJKL )THEN
          WRITE(LUWRT,*) ' Running in MPI-shared memory mode     '
          WRITE(LUWRT,*) ' or making use of revised integral     '
          WRITE(LUWRT,*) ' allocation (.IJKLSP) requires         '
          WRITE(LUWRT,*) ' the integral file IJKL_REOD.          '
          WRITE(LUWRT,*) ' Programm will abort now.              '
          WRITE(LUWRT,*) ' Restart job with option .IJKLRO !     '
          CALL QUIT(
     &    '*** Error in integral_distribution_driver ***'//
     &    ' Option .IJKLRO required. ')
        END IF
C
        CALL MEMMAR_I8(KINTI12,LEN_ALL_INT,'ADDL  ',2,'XTIT12')
C
        if(iiomod_rel.eq.1 .or. myproc.eq.master)then

!         get dimensions for scratch matrices
          call rgeth2_dim(ndim_ij,ndim_kl,ndim_h12_tmp)

          call memmar(kdumx12,idum,'MARK  ',idum,'geth12')
!         allocate space for scratch matrices (+ dummy space)
          is_space_req = 3*ndim_h12_tmp + 4*ndim_ij + 2*ndim_kl + 100

          call memmar(khmat,is_space_req,'ADDL  ',2,'hmat12')
          lfree = is_space_req - 2

!         make sure we always read from the integral file(s) associated
!         with the master process
          int_file_is_master_file = iiomod_rel .eq. 1 
     &                              .or. myproc .eq. master

!         get integrals on work(kinti12)
          call rgeth2(dummy,work(kinti12),dummy,ibeig,.false.,
     &                .true.,int_file_is_master_file,work(khmat),lfree)

!         release scratch space
          call memmar(kdumx12,idum,'FLUSM ',idum,'geth12')
        end if ! iiomod_rel == 1 or myproc == master
#if defined (VAR_MPI2)
        IF( IIOMOD_REL .eq. 0 )THEN
*           ... co-workers get the integrals
          LEN_INT_BCAST = 0
#if defined (INT_STAR8)
          LEN_INT_BCAST = LEN_ALL_INT
#else
          XLENALL = REAL(LEN_ALL_INT)
          LEN_INT_BCAST = IGIVE_I_B(XLENALL)
#endif
!         check for safeguard bcast limit
          if(LEN_INT_BCAST > integral_comm_limit)then
            batch_ints_loops     = LEN_INT_BCAST / integral_comm_limit
            remaining_batch_ints = mod(LEN_INT_BCAST, 
     &                                 integral_comm_limit)
            kloop_pointer        = 0
            do i = 1, batch_ints_loops 
              call interface_mpi_BCAST_r1_work_f77(
     &                       WORK(KINTI12+kloop_pointer),
     &                       integral_comm_limit,
     &                       MASTER,global_communicator)
              kloop_pointer = kloop_pointer + integral_comm_limit
            end do
!           the rest
              call interface_mpi_BCAST_r1_work_f77(
     &                       WORK(KINTI12+kloop_pointer),
     &                       remaining_batch_ints,
     &                       MASTER,global_communicator)
          else
            call interface_mpi_BCAST_r1_work_f77(WORK(KINTI12),
     &                     LEN_INT_BCAST,
     &                     MASTER,global_communicator)
          end if
        END IF
#endif
C
C       integrals on WORK(KINTI12) / F1-F2 --> 
C       distribute to WORK(KT_CC) resp. T_BUFF!
        call dist_ints_sigden_opt(INTIMP,ISPINFREE,ISYM_T,
     &                        WORK(KLSOBEX),WORK(KT_CC),
     &                        T_BUFF,
     &                        WORK(KISTR_CA),WORK(KISTR_CB),
     &                        WORK(KISTR_AA),WORK(KISTR_AB),
     &                        NDIM_TOT2,F1,WORK(KINTI12),
     &                        WORK(KMDIASV),IBOSYM,
     &                        WORK(KLABEXTP),CIRUN,
     &                        NTEST)
      ELSE IF( CIRUN .eq. 'PROP1 ') THEN
        IF( DOPROPREOD )THEN
CSK     NTEST = 20
        call dist_ints_sigden_opt(INTIMP,ISPINFREE,ISYM_T,
     &                        WORK(KLSOBEX),WORK(KT_CC),
     &                        T_BUFF,
     &                        WORK(KISTR_CA),WORK(KISTR_CB),
     &                        WORK(KISTR_AA),WORK(KISTR_AB),
     &                        NDIM_TOT2,F1,DUMMY,
     &                        WORK(KMDIASV),IBOSYM,
     &                        WORK(KLABEXTP),CIRUN,
     &                        NTEST)
CSK     NTEST = 00
        END IF
      ELSE
        call dist_ints_sigden_opt(INTIMP,ISPINFREE,ISYM_T,
     &                        WORK(KLSOBEX),WORK(KT_CC),
     &                        T_BUFF,
     &                        WORK(KISTR_CA),WORK(KISTR_CB),
     &                        WORK(KISTR_AA),WORK(KISTR_AB),
     &                        NDIM_TOT2,F1,F2,
     &                        WORK(KMDIASV),IBOSYM,
     &                        WORK(KLABEXTP),CIRUN,
     &                        NTEST)
C
      END IF ! CIRUN eq CIINII .or. CIRUN .eq. 'KR-CI '?
      call dealloc(ibeig)
C
      if(CIRUN .ne. 'PROP1 ' ) then
        write(LUWRT,'(/A,I16)') ' Number of integrals/CC amplitudes: ',
     &                            NALLINT
        write(LUWRT,'(A/)') ' (including imaginary parts) '
      end if
C
1001  CONTINUE
C
C     flush local memory to mark set at beginning of do loop
      CALL MEMMAR(KDUM,IDUM,'FLUSM ',IDUM,'ITRCTL')
C
1011  CALL QEXIT('idistdrv')
      END
***********************************************************************

      subroutine dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,ISOBEX,T,
     &                            T_BUFF,ISTR_CA,ISTR_CB,ISTR_AA,
     &                            ISTR_AB,NASHT,F1,F2,MDIASV,IBOSYM,
     &                            LABEL,CIRUN,IPRNT)
***********************************************************************
* Driver routine for distribution of DIRAC integrals to
* double group GAS lists
*
*   Routine using excitation operator sequence for integral
*   storage definition (SIGDEN route)
*
*----------------------------------------------
*  Timo Fleig, Mar. 2001 (Ellen Tabea is here!)
*  Jeppe Olsen modifies it. (Apr. 2001)
*  In complete form, June 2001, Timo Fleig
*
*  Implemented real groups (DIRAC) treated like complex groups
*  with zero imaginary parts.
*  Timo Fleig, August 25, 2003
*
*  implemented correct symmetry handling of one-electron property 
*  operator elements wrt to time-reversal symmetry
*  Stefan Knecht - Nov 2008
*
*----------------------------------------------
*
*  Implementation of T_BUFF use for MPIs "shared" memory version
* 
      use mospinor_info
      use symmetry_setup_krci
      use interface_to_mpi
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "ctcc.inc"
#include "cgas.inc"
#include "glbbas.inc"
#include "clunit.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IDISP
#endif
#include "parluci.h"
#include "files_r.inc"
#include "integrals_off.inc"
*
#include "dgroup.h"
*
      character*6 CIRUN
      dimension ISOBEX(NGAS,4,NSPOBEX_TP),
     &          LABEL(4,NSPOBEX_TP)
      dimension I_ADD(4),INDEX(4),IND_SAVE(4),
     &          MDIASV(4,*)           
*
* Group notation of the strings:
      integer   IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS),
     &          IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
      dimension T(NALLINT), T_BUFF(*)
*
* Space for the construction of SUBSTRINGS
* dimensioned outside as integers MX_ST_TSOSO_BLK_MX
      dimension ISTR_CA(*),ISTR_CB(*),ISTR_AA(*),ISTR_AB(*)
      INTEGER*8 ICOUNT_I8
      LOGICAL NOINT_FILE
      integer :: is_zero
*
      NTESTL = 0000
      NTEST = max(NTESTL,IPRNT/1000)
!     NTEST = 1
 
!     constant
      is_zero = 0
!
!     control writing of file INDICES.x (x =0,... last slave tag)
      NOINT_FILE = .FALSE.
      if( lowsrt_ijkl .or. reord_ijkl ) NOINT_FILE = .TRUE.
C     NOINT_FILE = .TRUE.
      IF( NOINT_FILE )THEN
        write(LUWRT,'(/A)') ' NOTICE from DIST_INTS_SIGDEN: '
        write(LUWRT,'( A)') ' INDICES file is not written as'
        write(LUWRT,'(/A)') ' requested by ".IJKLRO".       '
        write(LUWRT,'( A)') ' The calculation of natural    '
        write(LUWRT,'( A)') ' orbital occupations numbers   '
        write(LUWRT,'( A)') ' is therefore not allowed.     '
      END IF
*
!     check for inversion symmetry of group:
!     --------------------------------------
!     - no inversion symmetry: invsym == 0
!     -    inversion symmetry: invsym == 1

      invsym = 0
      if(DOUGRP.eq.1.or.DOUGRP.eq.4.or.DOUGRP.eq.7.or.
     &   DOUGRP.eq.9.or.DOUGRP.eq.11) INVSYM = 1

      IRILP = 2
      if(ISPINFREE.eq.1.or.nz.eq.1) IRILP = 1
!
!     open ASCII file for checking integral sorting information
!     ntest = 1
      if (NTEST.ge.1) then
        if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
          LUINGR = LUSC94
          open(unit=LUINGR,file=NLUSC94_MPI(1:LUSC94_MPIL),
     &         status='UNKNOWN',form='FORMATTED')
          rewind(LUINGR)
        end if
      end if
!
!     open file for saving index sequence
      LUINSV = LUSC95
      open(unit=LUINSV,file=NLUSC95_MPI(1:LUSC95_MPIL),
     &     status='UNKNOWN',form='FORMATTED')
      rewind(LUINSV)
!
!     open file with pre-fetched integrals
      if( lowsrt_ijkl ) then
#if defined (VAR_MPI2)
        call interface_mpi_file_open(self_communicator,"IJKL_PRESORT",
     &                     df_mpi_mode_rdonly + 
     &                     df_mpi_mode_delete_on_close,
     &                     df_mpi_info_null,ijkl_pre)
        idisp = 0
        call interface_mpi_file_set_view(ijkl_pre,idisp,
     &                         df_mpi_real8,df_mpi_real8,
     &                         "native",df_mpi_info_null)
#endif
      end if
!
!     set counter to start value and initialize
      ICOUNT_I8 = 1
      NDIAINT   = 0
      N1DIAINT  = 0
      N2DIAINT  = 0
!
!     Loop over real/imaginary part
!
      do IRI=1,IRILP,1
        if (NTEST.ge.1) then
          if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
            if (IRI.eq.1.and.NTEST.ge.2) write(LUWRT,*) 'REAL LOOP'
            if (IRI.eq.1) write(LUINGR,'(A)') 'Real  integrals.'
            if (IRI.eq.2.and.NTEST.ge.2) write(LUWRT,*) 'IMAG LOOP'
            if (IRI.eq.2) write(LUINGR,'(A)') 'Imag. integrals.'
          end if
        end if
        ITPSV = 0
        if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
          if (NTEST.ge.1) 
     &      write(LUINGR,'(A)') 'INTEGRAL FILE'
        end if

 
        print *, 'loop over NSPOBEX_TP: ',NSPOBEX_TP
* Loop over excitation types
        do IXTP=1,NSPOBEX_TP,1
*
* Check for correct pointers to type blocks
          if (IRI.eq.1.and.ICOUNT_I8.ne.IBTSOSO_TP(IXTP)) then
            write(LUWRT,*) 'Pointers to types do not match numbers'
            write(LUWRT,*) 'of integrals generated.'
            write(LUWRT,*) 'LABEL array: ',(LABEL(I,IXTP),I=1,4,1)
            write(LUWRT,*) 'int.no.,offset: ',ICOUNT_I8,IBTSOSO_TP(IXTP)
            write(LUWRT,*) 'Current type: ',IXTP
            call abend2('Quitting.')
          end if
          if (NTEST.ge.1) then
            write(LUWRT,*) '======================'
            write(LUWRT,*) 'Processing type ',IXTP
            write(LUWRT,*) '======================'
            write(LUWRT,*) 'ISOBEX partition:'
            call iwrtma(ISOBEX(1,1,IXTP),NGAS,4,NGAS,4)
          end if
*
* write statement about next integral type (to file)
          if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
            if (NTEST.ge.1) 
     &          call integral_type(ISOBEX(1,1,IXTP),NGAS,ITPSV,LUINGR)
            if (NTEST.ge.1) 
     &        write(LUINGR,'(A,I5)') 'Type ',IXTP
          end if
*
*. Change from Occupation to group notation of strings
          call occ_to_grp(ISOBEX(1,1,IXTP),IGRP_CA,1)
          call occ_to_grp(ISOBEX(1,2,IXTP),IGRP_CB,1)
          call occ_to_grp(ISOBEX(1,3,IXTP),IGRP_AA,1)
          call occ_to_grp(ISOBEX(1,4,IXTP),IGRP_AB,1)
          if (NTEST.ge.5.and.IRI.eq.1) then
            write(LUWRT,*) 'Group information:'
            write(LUWRT,'(A,4I6)') 'IGRP_CA : ',(IGRP_CA(I),I=1,NGAS,1)
            write(LUWRT,'(A,4I6)') 'IGRP_CB : ',(IGRP_CB(I),I=1,NGAS,1)
            write(LUWRT,'(A,4I6)') 'IGRP_AA : ',(IGRP_AA(I),I=1,NGAS,1)
            write(LUWRT,'(A,4I6)') 'IGRP_AB : ',(IGRP_AB(I),I=1,NGAS,1)
          end if
*
* Corresponding number of electrons
          NEL_CA = ielsum(ISOBEX(1,1,IXTP),NGAS)
          NEL_CB = ielsum(ISOBEX(1,2,IXTP),NGAS)
          NEL_AA = ielsum(ISOBEX(1,3,IXTP),NGAS)
          NEL_AB = ielsum(ISOBEX(1,4,IXTP),NGAS)
          if (NTEST.ge.2.and.IRI.eq.1) then
            write(LUWRT,*) 'NEL_CA = ',NEL_CA
            write(LUWRT,*) 'NEL_CB = ',NEL_CB
            write(LUWRT,*) 'NEL_AA = ',NEL_AA
            write(LUWRT,*) 'NEL_AB = ',NEL_AB
          end if
*
* Symmetry loops
* (Annihilator symm. is assumed in terms of corresponding creators!)
*  Symmetry restrictions according to number of electrons
          NEL_C = NEL_CA + NEL_CB
          NEL_A = NEL_AA + NEL_AB
          call evenodd(IEONELC,NEL_C)
          if (IEONELC.eq.1) then
            ISMCSTA = NIRR_PN + 1
            ISMCEND = NIRR_DG
          else if (IEONELC.eq.2) then
            ISMCSTA = 1
            ISMCEND = NIRR_PN
          end if
*
          call evenodd(IEONELCA,NEL_CA)
          if (IEONELCA.eq.1) then
            ISMCASTA = NIRR_PN + 1
            ISMCAEND = NIRR_DG
          else if (IEONELCA.eq.2) then
            ISMCASTA = 1
            ISMCAEND = NIRR_PN
          end if
*
          if (NTEST.ge.5.and.IRI.eq.1) then
            write(LUWRT,*) 'C  symmetry loop: ',ISMCSTA,ISMCEND
            write(LUWRT,*) 'CA symmetry loop: ',ISMCASTA,ISMCAEND
          end if
*
* Symmetry loops
          do ISM_C=ISMCSTA,ISMCEND,1
            ISM_A = IDBGMULT(ISM_C,IADJSYM(ISYM_T))
            do ISM_CA=ISMCASTA,ISMCAEND,1
              ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
              do ISM_AA=1,NIRR_DG,1
                ISM_AB = IDBGMULT(ISM_A,INVELM(ISM_AA))
                if (NTEST.ge.2000.and.IRI.eq.1) then
                  write(LUWRT,*) 'Creator symm.      :',ISM_C
                  write(LUWRT,*) 'Annihilator symm.  :',ISM_A
                  write(LUWRT,*) 'Creator alpha symm.:',ISM_CA
                  write(LUWRT,*) 'Creator beta symm. :',ISM_CB
                  write(LUWRT,*) 'Annihi. alpha symm.:',ISM_AA
                  write(LUWRT,*) 'Annihi. beta symm. :',ISM_AB
                end if
*
*. Generate strings of these symmetries
*. obtain strings
                IUB = 1
                call getstr2_totsm_spgp_rel(IUB,IGRP_CA,NGAS,ISM_CA,
     &                                      NEL_CA,NSTR_CA,ISTR_CA)
                IUB = 2
                call getstr2_totsm_spgp_rel(IUB,IGRP_CB,NGAS,ISM_CB,
     &                                      NEL_CB,NSTR_CB,ISTR_CB)
                IUB = 1
                call getstr2_totsm_spgp_rel(IUB,IGRP_AA,NGAS,ISM_AA,
     &                                      NEL_AA,NSTR_AA,ISTR_AA)
                IUB = 2
                call getstr2_totsm_spgp_rel(IUB,IGRP_AB,NGAS,ISM_AB,
     &                                      NEL_AB,NSTR_AB,ISTR_AB)
*
* Loops over strings for the given symmetries and types
* Loop over T elements as matrix T(I_CA,I_CB,I_AA,I_AB)
*  Determine possible connections
                write(LUWRT,*) 'Creator alpha strings:',NSTR_CA
                write(LUWRT,*) 'Creator beta  strings:',NSTR_CB
                write(LUWRT,*) 'Annihi. alpha strings:',NSTR_AA
                write(LUWRT,*) 'Annihi. beta  strings:',NSTR_AB

                if ((NSTR_CA.ge.1.or.NSTR_CB.ge.1).and.
     &              (NSTR_AA.ge.1.or.NSTR_AB.ge.1)) then
*
                  do I_AB=1,NSTR_AB,1
                    do I_AA=1,NSTR_AA,1
                      do I_CB=1,NSTR_CB,1
                        do I_CA=1,NSTR_CA,1
*
* Notes: 
*  1) Standard notation in second quantization for 2-el. ints is:
*     a+I a+K_ aL aJ_ (IJ|KL) = - a+I a+K_ aJ_ aL (IJ|KL)
*     The code proceeds according to the right hand side. The
*     signs are located in WORK(KSIGNNHX) for each excitation
*     class. So first indices 2 and 3 are interchanged, then
*     3 and 4 (implicit in the code).
*  2) Integral classes with 2 alpha or 2 beta creators (or both)
*     are treated differently: (due to exchange character)
*     a+I a+K aL aJ [(IJ|KL) - (IL|KJ)]
*  3) In the full spin-dependent formalism, extra classes of
*     integrals are required. 
* 
*. Indices of STRING I_CA are stored in IOCC_CA and first element is 
*  IOCC_CA(1+NEL_CA*(I_CA-1)) 
*
                        select case(NEL_C)

                          case(1)
                          IXCH = 0
*
                            if (NEL_CA.eq.1.and.NEL_CB.eq.0) then
                              INDEX(1) 
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                              INDEX(2) = 0
                            else if (NEL_CA.eq.0.and.NEL_CB.eq.1) then
                              INDEX(1) = 0
                              INDEX(2) 
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of creators.'
                              write(LUWRT,*) 
     &                              'NEL_CA,NEL_CB = ',NEL_CA,NEL_CB
                              call abend2('Quitting.')
                            end if
*
                            if (NEL_AA.eq.1.and.NEL_AB.eq.0) then
                              INDEX(3) 
     &                              = ISTR_AA(1+NEL_AA*(I_AA-1))
                              INDEX(4) = 0
                            else if (NEL_AA.eq.0.and.NEL_AB.eq.1) then
                              INDEX(3) = 0
                              INDEX(4) 
     &                              = ISTR_AB(1+NEL_AB*(I_AB-1))
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of annihil.'
                              write(LUWRT,*) 
     &                              'NEL_AA,NEL_AB = ',NEL_AA,NEL_AB
                              call abend2('Quitting.')
                            end if
*
                          case(2)
*
                            if (NEL_CA.eq.2.and.NEL_CB.eq.0) then
                              IXCH = 0
                              if (NEL_AA.eq.2) then
*                               Exchange integrals
                                IXCH = 1
                                do IEL=1,NEL_CA,1
                                  INDEX(IEL) = 
     &                                ISTR_CA(IEL+NEL_CA*(I_CA-1))
                                end do
                              else if (NEL_AB.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(1)
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                                INDEX(2)
     &                              = ISTR_CA(2+NEL_CA*(I_CA-1))
                              else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
*                               Spin-orbit integral
                                IXCH = 0
                                do IEL=1,NEL_CA,1
                                  INDEX(IEL) =
     &                                ISTR_CA(IEL+NEL_CA*(I_CA-1))
                                end do
                              end if
                            else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
                              IXCH = 0
                              INDEX(1) 
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                              INDEX(3) 
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                            else if (NEL_CA.eq.0.and.NEL_CB.eq.2) then
                              IXCH = 0
                              if (NEL_AB.eq.2) then
*                               Exchange integrals
                                IXCH = 1
                                do IEL=1,NEL_CB,1
                                  INDEX(IEL) = 
     &                                ISTR_CB(IEL+NEL_CB*(I_CB-1))
                                end do
                              else if (NEL_AA.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(1)
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                                INDEX(2)
     &                              = ISTR_CB(2+NEL_CB*(I_CB-1))
                              else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
*                               Spin-orbit integral
                                INDEX(1)
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                                INDEX(3)
     &                              = ISTR_CB(2+NEL_CB*(I_CB-1))
                              end if
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of creators.'
                              write(LUWRT,*) 
     &                              'NEL_CA,NEL_CB = ',NEL_CA,NEL_CB
                              call abend2('Quitting.')
                            end if
*
                            if (NEL_AA.eq.2.and.NEL_AB.eq.0) then
                              if (NEL_CA.eq.2) then
                                IXCH = 1
                                do IEL=1,NEL_AA,1
                                 INDEX(IEL+2) = 
     &                                ISTR_AA(IEL+NEL_AA*(I_AA-1))
                                end do
                              else if (NEL_CB.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(3)
     &                              = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(4)
     &                              = ISTR_AA(2+NEL_AA*(I_AA-1))
                              else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spin-orbit integral
                                do IEL=1,NEL_AA,1
                                 INDEX(IEL*2) =
     &                                ISTR_AA(IEL+NEL_AA*(I_AA-1))
                                end do
                              end if
                            else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
                              IXCH = 0
                              if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spinfree integral
                                INDEX(2) 
     &                                = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(4) 
     &                                = ISTR_AB(1+NEL_AB*(I_AB-1))
                              else
*                               Spin-orbit integral
                                INDEX(4)
     &                                = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(2)
     &                                = ISTR_AB(1+NEL_AB*(I_AB-1))
                              end if
                            else if (NEL_AA.eq.0.and.NEL_AB.eq.2) then
                              if (NEL_CB.eq.2) then
                                IXCH = 1
                                do IEL=1,NEL_AB,1
                                 INDEX(IEL+2) = 
     &                                ISTR_AB(IEL+NEL_AB*(I_AB-1))
                                end do
                              else if (NEL_CA.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                INDEX(3)
     &                              = ISTR_AB(1+NEL_AB*(I_AB-1))
                                INDEX(4)
     &                              = ISTR_AB(2+NEL_AB*(I_AB-1))
                              else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spin-orbit integral
                                do IEL=1,NEL_AB,1
                                 INDEX(IEL*2) =
     &                                ISTR_AB(IEL+NEL_AB*(I_AB-1))
                                end do
                              end if
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of annihil.'
                              write(LUWRT,*) 
     &                              'NEL_AA,NEL_AB = ',NEL_AA,NEL_AB
                              call abend2('Quitting.')
                            end if
*
                          case default
                            write(LUWRT,'(/a,i3)')
     &                      ' Number of e- for creator strings out of'//
     &                      ' bounds 0 < NEL_C < 3. current NEL_C ==>',
     &                        NEL_C
                            call quit('*** error in dist_ints_sigden:'//
     &                      ' wrong number of e- for creator strings')
                          end select
*
*  Stored integrals (DIRAC) refer to symmetry-type ordering.
*  Internal ordering is type-symmetry, so rename:
                          IR_IND = 0
                          JR_IND = 0
                          KR_IND = 0
                          LR_IND = 0
!                         stefan: refering to imosp_luci2dirac1 / IREOTS
!                                 is not consistent for a Kramers unrestricted implementation.
!                                 instead the reordering should occur
!                                 above where we know whether its
!                                 unbarred/barred. march 2011.
                          if(INDEX(1).ne.is_zero)then
                            IR_IND = imosp_luci2dirac1(IREOTS(INDEX(1)))
                          end if
                          if(INDEX(2).ne.is_zero)then
                            JR_IND = imosp_luci2dirac1(IREOTS(INDEX(2)))
                          end if
                          if(INDEX(3).ne.is_zero)then
                            KR_IND = imosp_luci2dirac1(IREOTS(INDEX(3)))
                          end if
                          if(INDEX(4).ne.is_zero)then
                            LR_IND = imosp_luci2dirac1(IREOTS(INDEX(4)))
                          end if
                          INDEX(1) = IR_IND
                          INDEX(2) = JR_IND
                          INDEX(3) = KR_IND
                          INDEX(4) = LR_IND
*
* Index reordering and second integral in exchange case:
* (Reason: CACA string generated e.g. for two alpha excitations,
*          but this does not correspond to index ordering in integral!)
*  Apply only to 2-electron integrals:
*
                          YXINT = 0.0d0
                          YINT  = 0.0d0
                          if (IXCH.eq.1) then
* Index sequence reads:
*  (IJ|KL) (which is not the operator sequence for the desired int!)
* Reorder for Coulomb integral:
*  (IJ|KL) --> (IK|LJ)
* Now create exchange integral:
*  (IJ|KL) --> (IJ|LK)
* (Overall, indices L and J are interchanged and a "-" introduced,
*  due to I+ K+ L J = - I+ K+ J L
* Coulomb:
                            I_SV = INDEX(1)
                            L_SV = INDEX(2)
                            K_SV = INDEX(3)
                            J_SV = INDEX(4)
                            INDEX(1) = I_SV
                            INDEX(2) = K_SV
                            INDEX(3) = L_SV
                            INDEX(4) = J_SV
* Exchange:
                            I_ADD(1) = I_SV
                            I_ADD(2) = J_SV
                            I_ADD(3) = L_SV
                            I_ADD(4) = K_SV
                            do IRUN = 1,4,1
                              I_ADD(IRUN) = 
     &                              I_ADD(IRUN) * 
     &                              LABEL(IRUN,IXTP)
                            end do
                            if (IRI.eq.1.or.(IRI.eq.2.and.NZ.ne.1)) then
                              if (CIRUN.ne.'DENS1 '
     &                            .and.CIRUN.ne.'DENS2 ')
     &                          call getint_sigden(YXINT,I_ADD,IRI,
     &                                    nasht/2,NASHT,F1,F2,NTEST)
                            end if
                          end if
*
                          if (CIRUN.ne.'DENS1 '.and.
     &                        CIRUN.ne.'DENS2 ') then
*
*  Modify if integrals are over relativistic spinors:
                            if ((INTIMP.eq.6.or.INTIMP.eq.7).and.
     &                           ISPINFREE.eq.0) then
                              if (INDEX(1).ne.0.and.
     &                            INDEX(2).ne.0) then
                                do IRUN = 1,4,1
                                    INDEX(IRUN) = 
     &                                  INDEX(IRUN) * 
     &                                  LABEL(IRUN,IXTP)
                                end do
*
* Correct internal ordering of two-bar exchange type integrals:
* (This might be implemented in some more elegant way later on.)
                                if (INDEX(1).gt.0.and.
     &                              INDEX(2).lt.0.and.
     &                              INDEX(3).lt.0.and.
     &                              INDEX(4).gt.0) then
                                  J_SV = INDEX(2)
                                  L_SV = INDEX(4)
                                  INDEX(2) = -L_SV
                                  INDEX(4) = -J_SV
                                end if
                              end if
                            end if
*
* Check for diagonal integrals
                            if (IRI.eq.1) then
                              ICNT = 0
                              do IR = 1,4,1
                                if (INDEX(IR).ne.0) then
                                  ICNT = ICNT + 1
                                  IND_SAVE(ICNT) = INDEX(IR)
                                end if
                              end do
                              if (ICNT.eq.2) then
                                if (IND_SAVE(1).eq.IND_SAVE(2)) then
                                  if (INDEX(1).ne.is_zero) IUBDIA = 1
                                  if (INDEX(2).ne.is_zero) IUBDIA = 2
                                  if (IUBDIA.eq.1) then
                                    NDIAINT = NDIAINT + 1
                                    N1DIAINT = N1DIAINT + 1
                                    IPOSDIA(NDIAINT) = ICOUNT_I8
                                    MDIASV(1,NDIAINT) = IND_SAVE(1)
                                    MDIASV(2,NDIAINT) = is_zero
                                    MDIASV(3,NDIAINT) = IND_SAVE(2)
                                    MDIASV(4,NDIAINT) = is_zero
C           print*,'D INT ',NDIAINT,' is ',(MDIASV(JR,NDIAINT),JR=1,4),
C    &             ' at ',ICOUNT_I8
                                  else if (IUBDIA.eq.2) then
                                    NDIAINT = NDIAINT + 1
                                    N1DIAINT = N1DIAINT + 1
                                    IPOSDIA(NDIAINT) = ICOUNT_I8
                                    MDIASV(1,NDIAINT) = is_zero
                                    MDIASV(2,NDIAINT) = IND_SAVE(1)
                                    MDIASV(3,NDIAINT) = is_zero
                                    MDIASV(4,NDIAINT) = IND_SAVE(2)
C           print*,'D INT ',NDIAINT,' is ',(MDIASV(JR,NDIAINT),JR=1,4),
C    &             ' at ',ICOUNT_I8
                                  end if
                                end if
                              else if (ICNT.eq.4) then
                                if ((INDEX(1).eq.INDEX(4)).and.
     &                              (INDEX(2).eq.INDEX(3)).or.
     &                              (INDEX(1).eq.INDEX(2)).and.
     &                              (INDEX(3).eq.INDEX(4)).or.
     &                              (INDEX(1).eq.-INDEX(3)).and.
     &                              (INDEX(2).eq.-INDEX(4))) then
                                  NDIAINT = NDIAINT + 1
                                  N2DIAINT = N2DIAINT + 1
                                  IPOSDIA(NDIAINT) = ICOUNT_I8
                                  do IR = 1,4,1
                                    MDIASV(IR,NDIAINT) = INDEX(IR)
                                  end do
C           print*,'D INT ',NDIAINT,' is ',(MDIASV(JR,NDIAINT),JR=1,4),
C    &             ' at ',ICOUNT_I8
                                end if
                              else
                                write(6,*) ' Diagonal check failed in '
                                write(6,*) ' dist_ints_sigden.        '
                                write(6,*) '  ICNT = ',ICNT
                                stop 'Quitting.'
                              end if
                            end if
*
*                           Write INDEX sequence to line in file
*                           for use in diagonal routine (saving memory):
                            IF(.NOT. NOINT_FILE )
     &                         write(LUINSV,'(4I4)') (INDEX(JJ),JJ=1,4)
*
                            if (IRI.eq.1.or.(IRI.eq.2.and.NZ.ne.1)) then
                              call getint_sigden(YINT,INDEX(1),
     &                                  IRI,nasht/2,NASHT,F1,F2,NTEST)
                            end if
C
                            T(ICOUNT_I8) = YINT - YXINT
C FIXME !! Uncommented sign here!!
                            if (IXCH.eq.1) T(ICOUNT_I8) = - T(ICOUNT_I8)
C FIXME !! Uncommented sign here!!
                              if (NTEST.ge.4) then
                                 write(LUWRT,'(A,I4,A,4I3,A,1P,1E16.7)')
     &                              'INTEGRAL no. ',ICOUNT_I8,'(',
     &                              (INDEX(M),M=1,4),
     &                              ') Value: ',T(ICOUNT_I8)
                                 write(LUWRT,'(A,3F18.14)') 
     &                              'INTEGRALs: YINT,YXINT,YINT-YXINT',
     &                                          YINT,YXINT,YINT-YXINT
*                             write integral to file (for debugging):
                                if (IXCH.eq.0) then
                                  write(LUINGR,'(I4,4X,4I3,1F18.14)') 
     &                                  ICOUNT_I8,INDEX(1),
     &                                         INDEX(2),
     &                                         INDEX(3),
     &                                         INDEX(4),T(ICOUNT_I8)
                                else
                       write(LUINGR,'(I4,4X,4I3,1X,A,1X,4I3,1P,1E16.7)')
     &                                  ICOUNT_I8,INDEX(1),
     &                                         INDEX(2),
     &                                         INDEX(3),
     &                                         INDEX(4),'-',
     &                              I_ADD(1),I_ADD(2),I_ADD(3),I_ADD(4),
     &                                         T(ICOUNT_I8)
                                end if
                              end if
                          end if
                          ICOUNT_I8 = ICOUNT_I8 + 1
*
                        end do
*                       ^ End loop over alpha creation strings
                      end do
*                     ^ End loop over beta creation strings
                    end do
*                   ^ End loop over alpha annihilation strings
                  end do
*                 ^ End loop over beta annihilation strings
                end if
*               ^ End if connection condition
              end do
*             ^ End loop over alpha annihilator symmetry
            end do
*           ^ End loop over alpha creator symmetry
          end do
*         ^ End loop over creation part symmetry
        end do
*       ^ End loop over excitation types
      end do
*     ^ End loop over real/imaginary
      if (NTEST.ge.1.and.(CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 '))
     &   close(LUINGR)
      close(LUINSV)
      NALLINT = ICOUNT_I8 - 1
C
 300  CONTINUE
C
      if (NTEST.ge.1) then
        write(LUWRT,*) ' DIST_INTS_SIGDEN reporting   : '
        write(LUWRT,*) '  Total number of integrals   : ',NALLINT
        write(LUWRT,*) '  Number of diagonal integrals: ',NDIAINT
        write(LUWRT,*) '    Number of diagonal 1-el. integrals: ',
     &                      N1DIAINT
        write(LUWRT,*) '    Number of diagonal 2-el. integrals: ',
     &                      N2DIAINT
        write(LUWRT,*) '    sorted integrals: '
        CALL WRTMATMN8(T,1,NALLINT,1,NALLINT,LUWRT)
      end if
      NTEST = 0
!
!     close (and destroy) file with pre-fetched integrals
      if( lowsrt_ijkl ) then
#if defined (VAR_MPI2)
        call interface_mpi_file_close(ijkl_pre)
#endif
      end if
!     call quit('*** bla bla: stop after dist_ints_sigden... ***')
*
      end

***********************************************************************

      subroutine dist_ints_sigden_opt(INTIMP,ISPINFREE,ISYM_T,ISOBEX,T,
     &                                T_BUFF,ISTR_CA,ISTR_CB,ISTR_AA,
     &                                ISTR_AB,NASHT,F1,F2,MDIASV,IBOSYM,
     &                                LABEL,CIRUN,IPRNT)
*
      use mospinor_info
      use symmetry_setup_krci
      use interface_to_mpi
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "ctcc.inc"
#include "cgas.inc"
#include "glbbas.inc"
#include "clunit.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IDISP
#endif
#include "parluci.h"
#include "files_r.inc"
#include "integrals_off.inc"
*
#include "dgroup.h"
*
      character*6 CIRUN
      dimension ISOBEX(NGAS,4,NSPOBEX_TP),
     &          LABEL(4,NSPOBEX_TP)
      dimension I_ADD(4),INDEX(4),IND_SAVE(4),
     &          MDIASV(4,*)           
*
* Group notation of the strings:
      integer   IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS),
     &          IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
      dimension T(NALLINT), T_BUFF(*)
*
* Space for the construction of SUBSTRINGS
* dimensioned outside as integers MX_ST_TSOSO_BLK_MX
      dimension ISTR_CA(*),ISTR_CB(*),ISTR_AA(*),ISTR_AB(*)
      INTEGER*8 ICOUNT_I8
      LOGICAL NOINT_FILE
      integer :: is_zero
*
!#define LUCI_DEBUG
!     NTESTL = 5000
      NTESTL = 0000
      NTEST = max(NTESTL,IPRNT/1000)
      NTEST = 0
 
!     constant
      is_zero = 0
!
!     control writing of file INDICES.x (x =0,... last slave tag)
      NOINT_FILE = .FALSE.
      if( lowsrt_ijkl .or. reord_ijkl ) NOINT_FILE = .TRUE.
C     NOINT_FILE = .TRUE.
      IF( NOINT_FILE )THEN
        write(LUWRT,'(/A)') ' NOTICE from DIST_INTS_SIGDEN: '
        write(LUWRT,'( A)') ' INDICES file is not written as'
        write(LUWRT,'(/A)') ' requested by ".IJKLRO".       '
        write(LUWRT,'( A)') ' The calculation of natural    '
        write(LUWRT,'( A)') ' orbital occupations numbers   '
        write(LUWRT,'( A)') ' is therefore not allowed.     '
      END IF
*
!     check for inversion symmetry of group:
!     --------------------------------------
!     - no inversion symmetry: invsym == 0
!     -    inversion symmetry: invsym == 1

      invsym = 0
      if(DOUGRP.eq.1.or.DOUGRP.eq.4.or.DOUGRP.eq.7.or.
     &   DOUGRP.eq.9.or.DOUGRP.eq.11) INVSYM = 1

      IRILP = 2
      if(ISPINFREE.eq.1.or.nz.eq.1) IRILP = 1
!
!     open ASCII file for checking integral sorting information
      if (NTEST.ge.2) then
        if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
          LUINGR = LUSC94
          open(unit=LUINGR,file=NLUSC94_MPI(1:LUSC94_MPIL),
     &         status='UNKNOWN',form='FORMATTED')
          rewind(LUINGR)
        end if
      end if
!
!     open file for saving index sequence
      LUINSV = LUSC95
      open(unit=LUINSV,file=NLUSC95_MPI(1:LUSC95_MPIL),
     &     status='UNKNOWN',form='FORMATTED')
      rewind(LUINSV)
!
!     open file with pre-fetched integrals
      if( lowsrt_ijkl ) then
#if defined (VAR_MPI2)
        call interface_mpi_file_open(self_communicator,"IJKL_PRESORT",
     &                 df_mpi_mode_rdonly + df_mpi_mode_delete_on_close,
     &                 df_mpi_info_null,ijkl_pre)
        idisp = 0
        call interface_mpi_file_set_view(ijkl_pre,idisp,
     &                         df_mpi_real8,
     &                         df_mpi_real8,
     &                         "native",df_mpi_info_null)
#endif
      end if
!
!     set counter to start value and initialize
      ICOUNT_I8 = 1
      NDIAINT   = 0
      N1DIAINT  = 0
      N2DIAINT  = 0
!
!     Loop over real/imaginary part
!
      do IRI=1,IRILP,1
        if (NTEST.ge.2) then
          if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
            if (IRI.eq.1.and.NTEST.ge.2) write(LUWRT,*) 'REAL LOOP'
            if (IRI.eq.1) write(LUINGR,'(A)') 'Real  integrals.'
            if (IRI.eq.2.and.NTEST.ge.2) write(LUWRT,*) 'IMAG LOOP'
            if (IRI.eq.2) write(LUINGR,'(A)') 'Imag. integrals.'
          end if
        end if
        ITPSV = 0
        if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
          if (NTEST.ge.1)then
            write(LUINGR,'(A)') 'INTEGRAL FILE'
            call dump_integral_info(f1,f2,IRI,nasht/2,NASHT,nz,
     &                              irilp,ipqtoq(1,0),luingr)
          end if
        end if

 
* Loop over excitation types
        do IXTP=1,NSPOBEX_TP,1
*
          if (NTEST.ge.1) then
            write(LUWRT,*) '======================'
            write(LUWRT,*) 'Processing type ',IXTP
            write(LUWRT,*) '======================'
            write(LUWRT,*) 'ISOBEX partition:'
            call iwrtma(ISOBEX(1,1,IXTP),NGAS,4,NGAS,4)
          end if
* Check for correct pointers to type blocks
          if (IRI.eq.1.and.ICOUNT_I8.ne.IBTSOSO_TP(IXTP)) then
            write(LUWRT,*) 'Pointers to types do not match numbers'
            write(LUWRT,*) 'of integrals generated.'
            write(LUWRT,*) 'LABEL array: ',(LABEL(I,IXTP),I=1,4,1)
            write(LUWRT,'(a,i12,a,i12)') 
     &      ' current offset:',ICOUNT_I8,' <==> pre-calculated:',
     &       IBTSOSO_TP(IXTP)
            write(LUWRT,*) 'Current type: ',IXTP
            call quit('*** error in dist_ints_...: offset mismatch.***')
          end if
*
* write statement about next integral type (to file)
          if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
            if (NTEST.ge.2) 
     &          call integral_type(ISOBEX(1,1,IXTP),NGAS,ITPSV,LUINGR)
            if (NTEST.ge.2) 
     &        write(LUINGR,'(A,I5)') 'Type ',IXTP
          end if
*
*. Change from Occupation to group notation of strings
          call occ_to_grp(ISOBEX(1,1,IXTP),IGRP_CA,1)
          call occ_to_grp(ISOBEX(1,2,IXTP),IGRP_CB,1)
          call occ_to_grp(ISOBEX(1,3,IXTP),IGRP_AA,1)
          call occ_to_grp(ISOBEX(1,4,IXTP),IGRP_AB,1)
          if (NTEST.ge.5.and.IRI.eq.1) then
            write(LUWRT,*) 'Group information:'
            write(LUWRT,'(A,4I6)') 'IGRP_CA : ',(IGRP_CA(I),I=1,NGAS,1)
            write(LUWRT,'(A,4I6)') 'IGRP_CB : ',(IGRP_CB(I),I=1,NGAS,1)
            write(LUWRT,'(A,4I6)') 'IGRP_AA : ',(IGRP_AA(I),I=1,NGAS,1)
            write(LUWRT,'(A,4I6)') 'IGRP_AB : ',(IGRP_AB(I),I=1,NGAS,1)
          end if
*
* Corresponding number of electrons
          NEL_CA = ielsum(ISOBEX(1,1,IXTP),NGAS)
          NEL_CB = ielsum(ISOBEX(1,2,IXTP),NGAS)
          NEL_AA = ielsum(ISOBEX(1,3,IXTP),NGAS)
          NEL_AB = ielsum(ISOBEX(1,4,IXTP),NGAS)

!         write(LUWRT,*) 'NEL_CA = ',NEL_CA
!         write(LUWRT,*) 'NEL_CB = ',NEL_CB
!         write(LUWRT,*) 'NEL_AA = ',NEL_AA
!         write(LUWRT,*) 'NEL_AB = ',NEL_AB
*
* Symmetry loops
* (Annihilator symm. is assumed in terms of corresponding creators!)
*  Symmetry restrictions according to number of electrons
          NEL_C = NEL_CA + NEL_CB
          NEL_A = NEL_AA + NEL_AB

!         set loop limits according to an even/odd number of e:
          ISMCSTA  = 1
          ISMCASTA = 1
          ISMCEND  = NIRR_PN
          ISMCAEND = NIRR_PN

          if(mod(NEL_C,2) .ne. 0)then
            ISMCSTA  = 1 + NIRR_PN
            ISMCEND  =     NIRR_DG
          end if
          if(mod(NEL_CA,2) .ne. 0)then
            ISMCASTA = 1 + NIRR_PN
            ISMCAEND =     NIRR_DG
          end if

          if(NTEST.ge.5.and.IRI.eq.1)then
            write(LUWRT,*) 'C  symmetry loop: ',ISMCSTA,ISMCEND
            write(LUWRT,*) 'CA symmetry loop: ',ISMCASTA,ISMCAEND
          end if
*
* Symmetry loops
          do 101 ISM_CA = ISMCASTA, ISMCAEND

            IUB = 1
            call getstr2_totsm_spgp_rel(IUB,IGRP_CA,NGAS,ISM_CA,
     &                                  NEL_CA,NSTR_CA,ISTR_CA)
            if(NSTR_CA .le. 0) goto 101

            do 202 ISM_AA = 1, NIRR_DG

              IUB = 1
              call getstr2_totsm_spgp_rel(IUB,IGRP_AA,NGAS,ISM_AA,
     &                                    NEL_AA,NSTR_AA,ISTR_AA)
              if(NSTR_AA .le. 0) goto 202

              do 303 ISM_C = ISMCSTA, ISMCEND

!               annihilator symmetry
                ISM_A =  IDBGMULT(ISM_C,IADJSYM(ISYM_T))
!               creator     symmetry "beta"/barred e-
                ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
!               annihilator symmetry "beta"/barred e-
                ISM_AB = IDBGMULT(ISM_A,INVELM(ISM_AA))

*. Generate strings of these symmetries
*. obtain strings
                call getstr2_totsm_spgp_rel(2,IGRP_CB,NGAS,ISM_CB,
     &                                      NEL_CB,NSTR_CB,ISTR_CB)
                call getstr2_totsm_spgp_rel(2,IGRP_AB,NGAS,ISM_AB,
     &                                      NEL_AB,NSTR_AB,ISTR_AB)

!             print *, 'ISM_CB:',ISM_CB,' NSTR_CB ==> ',NSTR_CB
!             print *, 'ISM_AB:',ISM_AB,' NSTR_AB ==> ',NSTR_AB

!             print *, 'ISM_A :',ISM_A 
!             print *, 'ISM_C :',ISM_C
*
* Loops over strings for the given symmetries and types
* Loop over T elements as matrix T(I_CA,I_CB,I_AA,I_AB)
*  Determine possible connections
                if ((NSTR_CB.ge.1).and.(NSTR_AB.ge.1))then
*
#ifdef LUCI_DEBUG
                  if(NTEST.ge.2000.and.IRI.eq.1) then
                    write(LUWRT,*) 'Creator symm.      :',ISM_C
                    write(LUWRT,*) 'Annihilator symm.  :',ISM_A
                    write(LUWRT,*) 'Creator alpha symm.:',ISM_CA
                    write(LUWRT,*) 'Creator beta symm. :',ISM_CB
                    write(LUWRT,*) 'Annihi. alpha symm.:',ISM_AA
                    write(LUWRT,*) 'Annihi. beta symm. :',ISM_AB
                  end if
#endif
*
                  do I_AB=1,NSTR_AB,1
                    do I_AA=1,NSTR_AA,1
                      do I_CB=1,NSTR_CB,1
                        do I_CA=1,NSTR_CA,1
*
* Notes: 
*  1) Standard notation in second quantization for 2-el. ints is:
*     a+I a+K_ aL aJ_ (IJ|KL) = - a+I a+K_ aJ_ aL (IJ|KL)
*     The code proceeds according to the right hand side. The
*     signs are located in WORK(KSIGNNHX) for each excitation
*     class. So first indices 2 and 3 are interchanged, then
*     3 and 4 (implicit in the code).
*  2) Integral classes with 2 alpha or 2 beta creators (or both)
*     are treated differently: (due to exchange character)
*     a+I a+K aL aJ [(IJ|KL) - (IL|KJ)]
*  3) In the full spin-dependent formalism, extra classes of
*     integrals are required. 
* 
*. Indices of STRING I_CA are stored in ISTR_CA and first element is 
*  ISTR_CA(1+NEL_CA*(I_CA-1)) 
*
                        select case(NEL_C)

                          case(1)
                            IXCH = 0
*
                            if (NEL_CA.eq.1.and.NEL_CB.eq.0) then
                              INDEX(1) 
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                              INDEX(2) = 0
                            else if (NEL_CA.eq.0.and.NEL_CB.eq.1) then
                              INDEX(1) = 0
                              INDEX(2) 
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of creators.'
                              write(LUWRT,*) 
     &                              'NEL_CA,NEL_CB = ',NEL_CA,NEL_CB
                              call abend2('Quitting.')
                            end if
*
                            if (NEL_AA.eq.1.and.NEL_AB.eq.0) then
                              INDEX(3) 
     &                              = ISTR_AA(1+NEL_AA*(I_AA-1))
                              INDEX(4) = 0
                            else if (NEL_AA.eq.0.and.NEL_AB.eq.1) then
                              INDEX(3) = 0
                              INDEX(4) 
     &                              = ISTR_AB(1+NEL_AB*(I_AB-1))
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of annihil.'
                              write(LUWRT,*) 
     &                              'NEL_AA,NEL_AB = ',NEL_AA,NEL_AB
                              call abend2('Quitting.')
                            end if
*
                          case(2)
                            IXCH = 0
*
                            if (NEL_CA.eq.2.and.NEL_CB.eq.0) then
                              if (NEL_AA.eq.2) then
*                               Exchange integrals
                                IXCH = 1
                                do IEL=1,NEL_CA,1
                                  INDEX(IEL) = 
     &                                ISTR_CA(IEL+NEL_CA*(I_CA-1))
                                end do
                              else if (NEL_AB.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(1)
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                                INDEX(2)
     &                              = ISTR_CA(2+NEL_CA*(I_CA-1))
                              else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
*                               Spin-orbit integral
                                do IEL=1,NEL_CA,1
                                  INDEX(IEL) =
     &                                ISTR_CA(IEL+NEL_CA*(I_CA-1))
                                end do
                              end if
                            else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
                              INDEX(1) 
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                              INDEX(3) 
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                            else if (NEL_CA.eq.0.and.NEL_CB.eq.2) then
                              if (NEL_AB.eq.2) then
*                               Exchange integrals
                                IXCH = 1
                                do IEL=1,NEL_CB,1
                                  INDEX(IEL) = 
     &                                ISTR_CB(IEL+NEL_CB*(I_CB-1))
                                end do
                              else if (NEL_AA.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(1)
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                                INDEX(2)
     &                              = ISTR_CB(2+NEL_CB*(I_CB-1))
                              else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
*                               Spin-orbit integral
                                INDEX(1)
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                                INDEX(3)
     &                              = ISTR_CB(2+NEL_CB*(I_CB-1))
                              end if
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of creators.'
                              write(LUWRT,*) 
     &                              'NEL_CA,NEL_CB = ',NEL_CA,NEL_CB
                              call abend2('Quitting.')
                            end if
*
                            if (NEL_AA.eq.2.and.NEL_AB.eq.0) then
                              if (NEL_CA.eq.2) then
                                IXCH = 1
                                do IEL=1,NEL_AA,1
                                 INDEX(IEL+2) = 
     &                                ISTR_AA(IEL+NEL_AA*(I_AA-1))
                                end do
                              else if (NEL_CB.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(3)
     &                              = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(4)
     &                              = ISTR_AA(2+NEL_AA*(I_AA-1))
                              else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spin-orbit integral
                                do IEL=1,NEL_AA,1
                                 INDEX(IEL*2) =
     &                                ISTR_AA(IEL+NEL_AA*(I_AA-1))
                                end do
                              end if
                            else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
                              IXCH = 0
                              if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spinfree integral
                                INDEX(2) 
     &                                = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(4) 
     &                                = ISTR_AB(1+NEL_AB*(I_AB-1))
                              else
*                               Spin-orbit integral
                                INDEX(4)
     &                                = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(2)
     &                                = ISTR_AB(1+NEL_AB*(I_AB-1))
                              end if
                            else if (NEL_AA.eq.0.and.NEL_AB.eq.2) then
                              if (NEL_CB.eq.2) then
                                IXCH = 1
                                do IEL=1,NEL_AB,1
                                 INDEX(IEL+2) = 
     &                                ISTR_AB(IEL+NEL_AB*(I_AB-1))
                                end do
                              else if (NEL_CA.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                INDEX(3)
     &                              = ISTR_AB(1+NEL_AB*(I_AB-1))
                                INDEX(4)
     &                              = ISTR_AB(2+NEL_AB*(I_AB-1))
                              else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spin-orbit integral
                                do IEL=1,NEL_AB,1
                                 INDEX(IEL*2) =
     &                                ISTR_AB(IEL+NEL_AB*(I_AB-1))
                                end do
                              end if
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of annihil.'
                              write(LUWRT,*) 
     &                              'NEL_AA,NEL_AB = ',NEL_AA,NEL_AB
                              call abend2('Quitting.')
                            end if
*
                          case default
                            write(LUWRT,'(/a,i3)')
     &                      ' Number of e- for creator strings out of'//
     &                      ' bounds 0 < NEL_C < 3. current NEL_C ==>',
     &                        NEL_C
                            call quit('*** error in dist_ints_sigden:'//
     &                      ' wrong number of e- for creator strings')
                          end select
*
*  Stored integrals (DIRAC) refer to symmetry-type ordering.
*  Internal ordering is type-symmetry, so rename:
                          IR_IND = 0
                          JR_IND = 0
                          KR_IND = 0
                          LR_IND = 0
!                         stefan: refering to imosp_luci2dirac1 / IREOTS
!                                 is not consistent for a Kramers unrestricted implementation.
!                                 instead the reordering should occur
!                                 above where we know whether its
!                                 unbarred/barred. march 2011.
                          if(INDEX(1).ne.is_zero)then
                            IR_IND = imosp_luci2dirac1(IREOTS(INDEX(1)))
                          end if
                          if(INDEX(2).ne.is_zero)then
                            JR_IND = imosp_luci2dirac1(IREOTS(INDEX(2)))
                          end if
                          if(INDEX(3).ne.is_zero)then
                            KR_IND = imosp_luci2dirac1(IREOTS(INDEX(3)))
                          end if
                          if(INDEX(4).ne.is_zero)then
                            LR_IND = imosp_luci2dirac1(IREOTS(INDEX(4)))
                          end if
#ifdef LUCI_DEBUG
                          print '(a,4i4)',' ==> ts indices: i,j,k,l',
     &                    INDEX(1), INDEX(2), INDEX(3), INDEX(4)
                          print '(a,4i4)',' ==> st indices: i,j,k,l',
     &                    IR_IND, JR_IND, KR_IND, LR_IND
#endif
                          INDEX(1) = IR_IND
                          INDEX(2) = JR_IND
                          INDEX(3) = KR_IND
                          INDEX(4) = LR_IND
*
* Index reordering and second integral in exchange case:
* (Reason: CACA string generated e.g. for two alpha excitations,
*          but this does not correspond to index ordering in integral!)
*  Apply only to 2-electron integrals:
*
                          YXINT = 0.0d0
                          YINT  = 0.0d0
                          if (IXCH.eq.1) then
* Index sequence reads:
*  (IJ|KL) (which is not the operator sequence for the desired int!)
* Reorder for Coulomb integral:
*  (IJ|KL) --> (IK|LJ)
* Now create exchange integral:
*  (IJ|KL) --> (IJ|LK)
* (Overall, indices L and J are interchanged and a "-" introduced,
*  due to I+ K+ L J = - I+ K+ J L
* Coulomb:
                            I_SV = INDEX(1)
                            L_SV = INDEX(2)
                            K_SV = INDEX(3)
                            J_SV = INDEX(4)
                            INDEX(1) = I_SV
                            INDEX(2) = K_SV
                            INDEX(3) = L_SV
                            INDEX(4) = J_SV
* Exchange:
                            I_ADD(1) = I_SV
                            I_ADD(2) = J_SV
                            I_ADD(3) = L_SV
                            I_ADD(4) = K_SV
                            do IRUN = 1,4,1
                              I_ADD(IRUN) = 
     &                              I_ADD(IRUN) * 
     &                              LABEL(IRUN,IXTP)
                            end do
                            if (IRI.eq.1.or.(IRI.eq.2.and.NZ.ne.1)) then
                              if (CIRUN.ne.'DENS1 '
     &                            .and.CIRUN.ne.'DENS2 ')
     &                          call getint_sigden(YXINT,I_ADD,IRI,
     &                                    nasht/2,NASHT,F1,F2,NTEST)
                            end if
                          end if
*
                          if (CIRUN.ne.'DENS1 '.and.
     &                        CIRUN.ne.'DENS2 ') then
*
*  Modify if integrals are over relativistic spinors:
                            if ((INTIMP.eq.6.or.INTIMP.eq.7).and.
     &                           ISPINFREE.eq.0) then
                              if (INDEX(1).ne.0.and.
     &                            INDEX(2).ne.0) then
                                do IRUN = 1,4,1
                                    INDEX(IRUN) = 
     &                                  INDEX(IRUN) * 
     &                                  LABEL(IRUN,IXTP)
                                end do
*
* Correct internal ordering of two-bar exchange type integrals:
* (This might be implemented in some more elegant way later on.)
                                if (INDEX(1).gt.0.and.
     &                              INDEX(2).lt.0.and.
     &                              INDEX(3).lt.0.and.
     &                              INDEX(4).gt.0) then
                                  J_SV = INDEX(2)
                                  L_SV = INDEX(4)
                                  INDEX(2) = -L_SV
                                  INDEX(4) = -J_SV
                                end if
                              end if
                            end if
*
* Check for diagonal integrals
                            if (IRI.eq.1) then
                              ICNT = 0
                              do IR = 1,4,1
                                if (INDEX(IR).ne.0) then
                                  ICNT = ICNT + 1
                                  IND_SAVE(ICNT) = INDEX(IR)
                                end if
                              end do
                              if (ICNT.eq.2) then
                                if (IND_SAVE(1).eq.IND_SAVE(2)) then
                                  if (INDEX(1).ne.is_zero) IUBDIA = 1
                                  if (INDEX(2).ne.is_zero) IUBDIA = 2
                                  if (IUBDIA.eq.1) then
                                    NDIAINT = NDIAINT + 1
                                    N1DIAINT = N1DIAINT + 1
                                    IPOSDIA(NDIAINT) = ICOUNT_I8
                                    MDIASV(1,NDIAINT) = IND_SAVE(1)
                                    MDIASV(2,NDIAINT) = is_zero
                                    MDIASV(3,NDIAINT) = IND_SAVE(2)
                                    MDIASV(4,NDIAINT) = is_zero
C           print*,'D INT ',NDIAINT,' is ',(MDIASV(JR,NDIAINT),JR=1,4),
C    &             ' at ',ICOUNT_I8
                                  else if (IUBDIA.eq.2) then
                                    NDIAINT = NDIAINT + 1
                                    N1DIAINT = N1DIAINT + 1
                                    IPOSDIA(NDIAINT) = ICOUNT_I8
                                    MDIASV(1,NDIAINT) = is_zero
                                    MDIASV(2,NDIAINT) = IND_SAVE(1)
                                    MDIASV(3,NDIAINT) = is_zero
                                    MDIASV(4,NDIAINT) = IND_SAVE(2)
C           print*,'D INT ',NDIAINT,' is ',(MDIASV(JR,NDIAINT),JR=1,4),
C    &             ' at ',ICOUNT_I8
                                  end if
                                end if
                              else if (ICNT.eq.4) then
                                if ((INDEX(1).eq.INDEX(4)).and.
     &                              (INDEX(2).eq.INDEX(3)).or.
     &                              (INDEX(1).eq.INDEX(2)).and.
     &                              (INDEX(3).eq.INDEX(4)).or.
     &                              (INDEX(1).eq.-INDEX(3)).and.
     &                              (INDEX(2).eq.-INDEX(4))) then
                                  NDIAINT = NDIAINT + 1
                                  N2DIAINT = N2DIAINT + 1
                                  IPOSDIA(NDIAINT) = ICOUNT_I8
                                  do IR = 1,4,1
                                    MDIASV(IR,NDIAINT) = INDEX(IR)
                                  end do
C           print*,'D INT ',NDIAINT,' is ',(MDIASV(JR,NDIAINT),JR=1,4),
C    &             ' at ',ICOUNT_I8
                                end if
                              else
                                write(6,*) ' Diagonal check failed in '
                                write(6,*) ' dist_ints_sigden.        '
                                write(6,*) '  ICNT = ',ICNT
                                stop 'Quitting.'
                              end if
                            end if
*
*                           Write INDEX sequence to line in file
*                           for use in diagonal routine (saving memory):
!#ifdef BLUBB
                            IF(.NOT. NOINT_FILE )
     &                         write(LUINSV,'(4I4)') (INDEX(JJ),JJ=1,4)
!#endif
*
                            if (IRI.eq.1.or.(IRI.eq.2.and.NZ.ne.1)) then
                              call getint_sigden(YINT,INDEX(1),
     &                                  IRI,nasht/2,NASHT,F1,F2,NTEST)
                            end if
C
                            T(ICOUNT_I8) = YINT - YXINT
C FIXME !! Uncommented sign here!!
                            if (IXCH.eq.1) T(ICOUNT_I8) = - T(ICOUNT_I8)
C FIXME !! Uncommented sign here!!

#ifdef DMRG_DEBUG
!                           DMRG DEBUG: no 2e-ints...
                            if(icnt > 2) T(ICOUNT_I8) = 0.0d0
!                           DMRG DEBUG: no 2e-ints...
#endif

                              if (NTEST.ge.2) then
                                write(LUWRT,'(A,I4,A,4I3,A,1P,1E16.7)') 
     &                              'INTEGRAL no. ',ICOUNT_I8,'(',
     &                              (INDEX(M),M=1,4),
     &                              ') Value: ',T(ICOUNT_I8)
                                 write(LUWRT,'(A,3F18.14)') 
     &                              'INTEGRALs: YINT,YXINT,YINT-YXINT',
     &                                          YINT,YXINT,YINT-YXINT
*                               write integral to file (for debugging):
!                               write(LUINGR,'(2x,a,2i5)')
!    &                          'NSTR_CA, NSTR_CB :',NSTR_CA, NSTR_CB
!                               write(LUINGR,'(2x,a,2i5)')
!    &                          'NSTR_AA, NSTR_AB :',NSTR_AA, NSTR_AB
                                if (IXCH.eq.0) then
                                  write(LUINGR,'(I4,4X,4I3,1F18.14)') 
     &                                  ICOUNT_I8,INDEX(1),
     &                                         INDEX(2),
     &                                         INDEX(3),
     &                                         INDEX(4),T(ICOUNT_I8)
                                else
                       write(LUINGR,'(I4,4X,4I3,1X,A,1X,4I3,1P,1E16.7)')
     &                                  ICOUNT_I8,INDEX(1),
     &                                         INDEX(2),
     &                                         INDEX(3),
     &                                         INDEX(4),'-',
     &                              I_ADD(1),I_ADD(2),I_ADD(3),I_ADD(4),
     &                                         T(ICOUNT_I8)
                                end if
                              end if
                          end if
                          ICOUNT_I8 = ICOUNT_I8 + 1
!                         print *, 'new offset for type IXTP',
!    &                    ICOUNT_I8,IXTP
*
                        end do
*                       ^ End loop over alpha creation strings
                      end do
*                     ^ End loop over beta creation strings
                    end do
*                   ^ End loop over alpha annihilation strings
                  end do
*                 ^ End loop over beta annihilation strings
                end if
*               ^ End if connection condition
 303          continue
*             ^ End loop over alpha annihilator symmetry
 202        continue
*           ^ End loop over alpha creator symmetry
 101      continue
*         ^ End loop over creation part symmetry
        end do
*       ^ End loop over excitation types
      end do
*     ^ End loop over real/imaginary
      if (NTEST.ge.2.and.(CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 '))
     &   close(LUINGR)
      close(LUINSV)
      NALLINT = ICOUNT_I8 - 1
C
 300  CONTINUE
C
      if (NTEST.ge.2) then
        write(LUWRT,*) ' DIST_INTS_SIGDEN reporting   : '
        write(LUWRT,*) '  Total number of integrals   : ',NALLINT
        write(LUWRT,*) '  Number of diagonal integrals: ',NDIAINT
        write(LUWRT,*) '    Number of diagonal 1-el. integrals: ',
     &                      N1DIAINT
        write(LUWRT,*) '    Number of diagonal 2-el. integrals: ',
     &                      N2DIAINT
        write(LUWRT,*) '    sorted integrals: '
        CALL WRTMATMN8(T,1,NALLINT,1,NALLINT,LUWRT)
      end if
!
!     close (and destroy) file with pre-fetched integrals
      if( lowsrt_ijkl ) then
#if defined (VAR_MPI2)
        call interface_mpi_file_close(ijkl_pre)
#endif
      end if
*
#undef LUCI_DEBUG
      end
***********************************************************************
      subroutine dump_integral_info(h1,h2ac,iri,nashth,nasht,nz,
     &                              irilp,ipqtoq,luwrt)

!     dump all 1e/2e integrals on file INTEGRALS in nice format

      implicit none
!     input
      real(8), intent(in)    :: h1((2*nasht)**2 * irilp)
      real(8), intent(in)    :: h2ac(nashth,nashth,
     &                              (nashth*(nashth+1))/2,nz,*)
      integer, intent(in)    :: iri, irilp,nashth, nasht, nz, luwrt
      integer, intent(in)    :: ipqtoq(4,0:7)

!     scratch
      integer                :: i, j , ii, jj

      write(luwrt,'(/A)')' print of unsorted integrals from MOLTRA'
      write(luwrt,'(A/)')' ---------------------------------------'

      ! 1e-integrals
      write(luwrt,'(a/)') ' list of 1e-integrals'
      write(luwrt,'(a)' ) ' indices         elements       '
      write(luwrt,'(a)' ) '--------------------------------'
      do i = 1, nasht
        do j = 1, nasht
          ii = i
          jj = j
          if(i > nashth) ii = -(i-nashth)
          if(j > nashth) jj = -(j-nashth)
          write(luwrt,'(i4,i4,5x,1E19.12)') 
     &    jj,ii,h1(((iri-1)*((nasht)**2))+(j-1)*nasht + i)
        end do  
      end do 

      ! 2e-integrals
      write(luwrt,'(/A)')' list of 2e-integrals in (NZ,3) format'
      call prdnz3(h2ac(1,1,1,1,1),nashth,(nashth*(nashth+1))/2,nz,
     &            ipqtoq(1,0),luwrt)

      write(luwrt,'(/A)')' start of sorted integrals for KRCI'
      write(luwrt,'(A/)')' ----------------------------------'

      end 
***********************************************************************

      subroutine getind_nz3(YINT,INDEX,IRI,IDIM,NASHT,NNASHX,H2AC,IPRNT)
!
!     purpose: determine 2- particle integral and fetch integral.
*
*=========================
* Timo Fleig, Oct 20, 2004
*=========================
*
*
      use interface_to_mpi
#include "implicit.h"
#include "dgroup.h"
#if defined (VAR_MPI2)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFFSET
#endif
#include "parluci.h"
*
      dimension H2AC(NASHT/2,NASHT/2,NNASHX,NZ,*)
      dimension ICSAVE(4),INDEX(4),INDEX_LOC(4)
*
      do IL = 1,4
        INDEX_LOC(IL) = INDEX(IL)
      end do
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT/100)
*
      if (IRI.eq.1) then
        MZ = 1
      else if (IRI.eq.2) then
        MZ = NZ
      end if
*
      if (NTEST.ge.2) then
        write(6,'(/A,6X,4I3)') '  INDEX array: ',(INDEX(I),I=1,4,1)
      end if
      if (NTEST.ge.5) then
        write(6,'(/A)') ' GETIND_NZ3 : '
        write(6,*) '  NZ,MZ : ',NZ,MZ
        write(6,*) '  No. of active orbitals : ',NASHT/2
        write(6,*) '  Compound index dimension: ',NNASHX
      end if
*
      if (NTEST.ge.200) then
        if( .not. lowsrt_ijkl ) then
          write(6,'(/A)')'2-el. integrals in getind_nz3, (NZ,3) format:'
          call prdnz3(H2AC(1,1,1,1,1),NASHT/2,NNASHX,NZ,IPQTOQ(1,0),6)
        end if
      end if
*
      NCT = 0
      do IC=1,4,1
        if (INDEX_LOC(IC).ne.0) then
          NCT = NCT + 1
          ICSAVE(NCT) = IC
        end if
      end do
*
* Interface to (NZ,3) format for 2-electron integrals
*
* Use Kramers symmetry/particle symmetry to generate a stored  
* index sequence on (NZ,3) format:
*  NZ3TP = 1     (pq|rs)
*  NZ3TP = 2     (pq_|rs_)
*  NZ3TP = 3     (p_q|rs_)
      IRIFAC = 1
      if (INDEX_LOC(1).lt.0.and.INDEX_LOC(2).lt.0.and.
     &    INDEX_LOC(3).lt.0.and.INDEX_LOC(4).lt.0) then
        call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,-1,1)
        call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,-1,1)
        NZ3TP = 1
      else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).gt.0.and.
     &         INDEX_LOC(3).gt.0.and.INDEX_LOC(4).gt.0) then
        NZ3TP = 1
      else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).gt.0.and.
     &         INDEX_LOC(3).lt.0.and.INDEX_LOC(4).lt.0) then
        call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,-1,1)
        NZ3TP = 1
      else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).lt.0.and.
     &         INDEX_LOC(3).lt.0.and.INDEX_LOC(4).gt.0) then
        call index_exch(INDEX_LOC(1),INDEX_LOC(3),1,1,1)
        call index_exch(INDEX_LOC(2),INDEX_LOC(4),1,1,1)
        INDEX_LOC(1) = -INDEX_LOC(1)
        INDEX_LOC(4) = -INDEX_LOC(4)
        NZ3TP = 3
      else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).lt.0.and.
     &         INDEX_LOC(3).gt.0.and.INDEX_LOC(4).lt.0) then
        INDEX_LOC(2) = -INDEX_LOC(2)
        INDEX_LOC(4) = -INDEX_LOC(4)
        NZ3TP = 2
      else if (INDEX_LOC(1).lt.0.and.INDEX_LOC(2).gt.0.and.
     &         INDEX_LOC(3).lt.0.and.INDEX_LOC(4).gt.0) then
        if (IRI.eq.2) IRIFAC = -1
        call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,1,1)
        call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,1,1)
        INDEX_LOC(2) = -INDEX_LOC(2)
        INDEX_LOC(4) = -INDEX_LOC(4)
        NZ3TP = 2
      end if
      if (NTEST.ge.2) write(6,'(A,1X,4I3)') 
     &     '  After 3 list step:', (INDEX_LOC(I),I=1,4,1)
*
* Find smallest index and its position in Mulliken notation
      call findlow(INDEX_LOC,4,ILOW,IPOS)
*
* Transfer smallest to fourth index by using 
*   particle symmetry 
*   Kramers symmetry
*   complex conjugation symmetry (IRIFAC2 = -1 if imaginary part)
*   CAUTION HERE:
*    The sign has to be derived from the original index sequence!
*    Thus the below sign arrangements may look queer.
*    They are nevertheless correct.
*
      IRIFAC2 = 1
      if (IPOS.eq.1) then 
        call index_exch(INDEX_LOC(1),INDEX_LOC(3),1,1,1)
        call index_exch(INDEX_LOC(2),INDEX_LOC(4),1,1,1)
        if (IRI.eq.2.and.NZ3TP.ne.2) IRIFAC2 = -1
        call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,1,1)
        call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,1,1)
      else if (IPOS.eq.2) then
        if (IRI.eq.2.and.NZ3TP.eq.3) IRIFAC2 = -1
        call index_exch(INDEX_LOC(1),INDEX_LOC(3),1,1,1)
        call index_exch(INDEX_LOC(2),INDEX_LOC(4),1,1,1)
      else if (IPOS.eq.3) then
        if (IRI.eq.2.and.NZ3TP.eq.1) IRIFAC2 = -1
        call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,1,1)
        call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,1,1)
      end if
      if (NTEST.ge.2) write(6,'(A,1X,4I3)') 
     &   '  After posit. step:', (INDEX_LOC(I),I=1,4,1)
*
* Determine address in (NZ,3) format
      NASHTH = NASHT/2
*
* Determine compound index KL: NDIM2
      NDIM2 = 0
      do IC = 1,INDEX_LOC(3),1
        if (IC.eq.INDEX_LOC(3)) then
          JCEND = INDEX_LOC(4)
        else if (IC.ne.INDEX_LOC(3)) then
          JCEND = IC
        end if
        do JC = 1,JCEND,1
          NDIM2 = NDIM2 + 1
        end do
      end do
*
      IIN = (NZ3TP-1)        *MZ*NNASHX*NASHTH*NASHTH +
     &      (MZ-1)           *NNASHX*NASHTH*NASHTH +
     &      (NDIM2-1)        *NASHTH*NASHTH +
     &      (INDEX_LOC(2)-1) *NASHTH +
     &       INDEX_LOC(1)
*
C     YINT = H2AC(IIN) * IRIFAC * IRIFAC2
      if( lowsrt_ijkl ) then
#if defined (VAR_MPI2)
        yint4 = 0.0D0
!       file offset is 0!
        ioffset = iin - 1
        call interface_mpi_file_read_at(ijkl_pre,ioffset,yint4,1,
     &                        istat)
!       write(6,*) 'iin ,ioffset, yint4',iin ,ioffset, yint4
        yint = yint4 * irifac * irifac2
#endif
      else
        YINT = H2AC(INDEX_LOC(1),INDEX_LOC(2),NDIM2,MZ,NZ3TP) 
     &              * IRIFAC * IRIFAC2
      end if
!     ntest = 2
      if (NTEST.ge.2) write(6,'(A,5I3,2X,4I3,2X,1F18.10)') 
     &                      'MZ,NZ3TP,IRIFAC,IRIFAC2,IIN  integral',
     &                       MZ,NZ3TP,IRIFAC,IRIFAC2,IIN,
     &                       (INDEX_LOC(ID),ID=1,4,1),YINT
!     ntest = 0
      end
***********************************************************************

      subroutine getint_nts(yint,index,iri,idim,dmfc,iprnt)
C***********************************************************************
C
C     fetch integral to YINT for non-totally symmetric operators 
C     from Molfdir-type array DMFC.  
C
C     order in DIRAC (Molfdir format) matrix DMFC:
C                                 _
C           p spinors followed by p
C                                                 _     _
C           if NFSYM = 2 (inversion symmetry): p, p, q, q
C
C    *********************************************************
C    * NOTE: the ordering is GA space independent, e.g.      *
C    *                                                       *
C    * .GASSH   yields the very same DMFC matrix as   .GASSH *
C    * 3                                              2      *
C    * 1 1                                            1 1    *
C    * 0 2                                            1 2    *
C    * 1 0                                                   *
C    *********************************************************
C
C     The present routine exhibits a new and more general offset 
C     calculation for accessing the integral matrix DMFC compared to 
C     GETIND_NZ3. The old offset calculation in GETIND_NZ3 is *NOT* valid 
C     for general one-electron operators.
C     The new one looks a little more complicated though...
C
C     Written by S. Knecht - Dec 2009
C
C     Last revision :
C
C***********************************************************************
!     use memory_allocator
#include "implicit.h"
#include "dgroup.h"
#include "dcborb.h"
!
      DIMENSION DMFC(*), ICSAVE(4), INDEX(4)
      integer jc_fsym, ja_fsym, nasht2, n2ashx2
!
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRNT)
!
!     set dimensions for DMFC matrix. Note that we here access the
!                         _   _
!     matrix in spinor (p,p,q,q) and NOT quaternion-packed Kramers-paired 
!     basis.
!     E.g. IASH(2) refers in the usual Dirac notation to the 
!     offset in fermion symmetry 2 wrt to the Kramers-paired basis. 
!     This means that we must multiply the offset by 2 to get the 
!     correct spinor offset.
!
      nasht2  = NASHT  * 2
      n2ashx2 = nasht2 * nasht2
 
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
!     input print
      WRITE(6,*) ' GETIND_NTS : '
      WRITE(6,*) '  No. of active orbitals : ',NASHT
      WRITE(6,*) '  Compound index dimension: ',NNASHX
      WRITE(6,'(A,6X,4I3)') '  INDEX array: ',(INDEX(I),I=1,4,1)
      WRITE(6,*) '1-el. integral array in GETIND_NTS:'
!     CALL WRTMAT(DMFC,IDIM*2,IDIM*2,IDIM*2,IDIM*2)
      if(iri == 1)then
        write(6,'(a/)') ' real part of 1e-integrals'
      else
        write(6,'(a/)') ' imaginary part of 1e-integrals'
      end if
      write(6,*) ' indices         elements       '
      write(6,*) '--------------------------------'
      do i = 1, nasht * 2
        do j = 1, nasht * 2
          ii = i
          jj = j
          if(i > nasht) ii = -ii
          if(j > nasht) jj = -jj
          write(6,'(i4,i4,5x,1E19.12)') 
     &    jj,ii,DMFC(((IRI-1)*n2ashx2)+(j-1)*nasht*2 + i)
        end do  
      end do  
#endif
!#undef LUCI_DEBUG
!
      NCT = 0
      DO IC=1,4,1
        IF( INDEX(IC).ne. 0)THEN
          NCT = NCT + 1
          ICSAVE(NCT) = IC
        END IF
      END DO
!
!     NCT = 2: pick out one-electron integral identified by 
!     two indices IORI (creator index) and JORI (annihilator index)
      IORI = INDEX(ICSAVE(1))
      JORI = INDEX(ICSAVE(2))
!
!     get fermion symmetry of creator and annihilator 
      jc_fsym = 1
      ja_fsym = 1
      IF( INDEX(ICSAVE(1)) .gt. NASH(1)) jc_fsym = 2
      IF( INDEX(ICSAVE(2)) .gt. NASH(1)) ja_fsym = 2
!
!     compute input address on array DMFC:
!     compute a relative creator and annihilator spinor 
!     (not Kramers pair) index
!
!     1. creator index
      I_REL = IORI + ((ICSAVE(1)-1)*NASH(jc_fsym))+IASH(jc_fsym)
!
!     2. annihilator index
      J_REL = JORI + ((ICSAVE(2)-3)*NASH(ja_fsym))+IASH(ja_fsym)
!
!     combine I_REL and J_REL to compute the address IIN
!     --------------------------------------------------
!     first  line: offset wrt real/imaginary
!     second line: offset wrt annihilator index
!     third  line: offset wrt the creator index
      IIN =  ((IRI-1)  * n2ashx2) +
     &       ((J_REL-1)* nasht2)  +
     &         I_REL
!     fetch integral
      YINT = DMFC(IIN)
 
!#define LUCI_DEBUG
!     print section
#ifdef LUCI_DEBUG
      write(6,*) 'orig. index pair ',IORI,JORI,' address :',IIN
      write(6,'(A,I5)') ' I_REL = ',I_REL
      write(6,'(A,I5)') ' J_REL = ',J_REL
!     DO LLLL = 1, (n2ashx2 * NZ)
!       WRITE(6,*) ' element ',LLLL,'in matrix DMFC is = ',DMFC(LLLL)
!     END DO
      write(6,*) '1-electron integral.'
      write(6,'(3I4,2X,1E14.8)') I_REL,J_REL,IIN,YINT
#endif
#undef LUCI_DEBUG
      end
***********************************************************************

      subroutine getint_sigden(yint,index,iri,itotdim,nasht,f1,f2,iprnt)
!
!     purpose: fetch integral from integral list according to SIGDEN
!              excitation class formalism - driver routine
      implicit none
      integer, intent(in)    :: index(4)
      integer, intent(in)    :: iri, itotdim, nasht, iprnt
      real(8), intent(inout) :: yint, f1, f2
!----------------------------------------------------------------------
      integer                :: ic, nct, nnashx
!----------------------------------------------------------------------

      nct = 0
      do ic = 1,4
        if(INDEX(ic).ne.0) nct = nct + 1
      end do

      select case(nct)
        case (2)
!         1e- ints
          call getint_nts(yint,index,iri,itotdim,f1,iprnt)
        case (4) ! 2-e integral
!         2e- ints in DIRAC (NZ,3) format
          nnashx = ((nasht/2) * ((nasht/2)+1)) / 2
          call getind_nz3(yint,index,iri,itotdim,nasht,nnashx,f2,iprnt)
        case default
          print '(a,i4)', '*** error in getint_sigden: wrong number'//
     &          ' of indices, neither 2 or 4 but ==>',nct
          call quit('*** error in getint_sigden: wrong number of'//
     &              ' indices.***')
      end select

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_PRPOPSYM_KRCI(WORK,LWORK,NTEST)
C***********************************************************************
C
C     calculate operator symmetry of all property operators NPROP_KRCI 
C     in list LPROP_KRCI in LUCIAREL sub-double group notation.
C     
C     symmetry of property operator IVOP according to
C      _           _              _      
C     | (IVOP) =  | (IBSP_KRCI) x | (IASP_KRCI)*
C
C     final symmetry stored in ISYMOPPRP_KRCI (krciprop.h).
C     
C     symmetry setup in LUCIAREL
C     ==========================
C
C             boson               fermion
C             =====               =======
C
C     C2   :  a , b               1E,  2E
C     C2h  :  ag, bg, au, bu      1Eg, 2Eg, 1Eu, 2Eu  
C
C     Written by S. Knecht - Sep 2008
C
C     Last revision :
C
C***********************************************************************
      use symmetry_setup_krci
#include "implicit.h"
#include "priunit.h"
#include "consts.h"
C
#include "dcbgen.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "krciprop.h"
#include "krmcluci_inf.h"
#include "dcbkrci.h"
#include "mxpdim.inc"
      DIMENSION WORK(*)
      CHARACTER REPVOP_KRCI*3, FILELAB*8
      integer, allocatable :: nroots_temp(:)
      integer              :: sym_xy
C
#include "memint.h"

      CALL QENTER('GET_PRPOPSYM_KRCI')
      MZ = MIN(NZ,2)

!     consistency check for property runs (TDM/perm. dipmom) + linear symmetry
!
!     check for either all roots being the unbarred or barred Kramers partner of 
!     a Kramers pair - if mixed we might get funny results 
!     for the (x,y) (unbarred/barred-classified) operators.
!      
      if(linear.and.(dooscillst.or.dodipmom_krci))then
        allocate(nroots_temp(nirr_dg))
        call icopy(nirr_dg,nkrci_ciroots,1,nroots_temp,1)

        nr_unbarred = 0
        nr_barred   = 0
        if(mod(nkrciaelec,2).eq.0)then ! even # e-
          nr_unbarred = nr_unbarred + nroots_temp(1)
          nroots_temp(1) = 0
          if(nfsym.eq.2) 
     &    nr_unbarred = nr_unbarred + nroots_temp((nirr_dg/4)+1)
          nroots_temp((nirr_dg/4)+1) = 0
          do i = 2, nirr_dg, 2
            nr_unbarred = nr_unbarred + nroots_temp(i)
          end do
          do i = 3, nirr_dg, 2
            nr_barred = nr_barred + nroots_temp(i)
          end do
        else ! odd # e-
          nr_unbarred = nr_unbarred + nroots_temp((nirr_dg/2)+1)
          nroots_temp((nirr_dg/2)+1) = 0
          if(nfsym.eq.2) 
     &    nr_unbarred = nr_unbarred + nroots_temp((nirr_dg/4)*3+1)
          nroots_temp((nirr_dg/4)*3+1) = 0
          do i = 1, nirr_dg, 2
            nr_unbarred = nr_unbarred + nroots_temp(i)
          end do
          do i = 2, nirr_dg, 2
            nr_barred = nr_barred + nroots_temp(i)
          end do
        end if
        deallocate(nroots_temp)

        sym_xy = 2 ! unbarred x/y as default
        if(nr_unbarred.lt.nr_barred) sym_xy = 1

        if(nr_unbarred.ne.0.and.nr_barred.ne.0)then
          print *, ' *** warning: unbarred/barred roots have been'//
     &             ' mixed in a TDM/perm. dipmom run. do not'//
     &             ' necessarily trust the x/y components;'
          print *, ' this is because of symmetry instability with the'//
     &             ' x/y operators which due to numerical noise might'
          print *, ' be wrongly classified as *unbarred*/*barred* '//
     &             '(u/b) or vice versa instead of ub/ub or b/b.'
        end if
      end if
!
!     allocate for VOP matrix
      CALL MEMGET('REAL',KVOP,4*N2ASHX*MZ,WORK,KFREE,LFREE)
!
      IASP_DC       = 0
      IBSP_DC       = 0
      IASP_KRCI_SYM = 0
      IBSP_KRCI_SYM = 0
!
!     loop over list of operators  
      DO 20 I = 1, NPROP_KRCI
 
        CALL DZERO(WORK(KVOP),  4*N2ASHX*MZ)
!
!       read VOP matrix from file
        INDXPR  = LPROP_KRCI(I)
        FILELAB = PRPNAM(INDXPR)(1:8)

!       print *,'FILELAB is ==> ',FILELAB
!       print *,'INDXPR  is ==> ',INDXPR
 
        CALL KRCI_PRPFILE(LUKRM3,FILELAB,WORK(KVOP),
     &                    4*N2ASHX*MZ,1)
C   
C       find largest absolute element in VOP matrix
C
        INDEX_VOP   = IDAMAX(4*N2ASHX*MZ,WORK(KVOP),1)
        ABS_MAX_ELM = DABS(WORK(KVOP+INDEX_VOP-1))
C       ... if largest absolute element is zero, 
C           nothing to do for this operator
        IF(ABS_MAX_ELM .eq. 0.0d0) GOTO 20
C
        IF( NTEST .gt. 10 ) THEN
           WRITE(LUPRI,'(/A,A8)') 
     &     '  absolute largest element for operator ',FILELAB
           WRITE(LUPRI,*) ' found to be:', ABS_MAX_ELM
        END IF
C
C       which two spinors (in LUCIAREL symmetry notation) 
C       IASP_KRCI_SYM and IBSP_KRCI_SYM does this operator connect?
C
        CALL GET_IAIB_DIRAC(IASP_KRCI_SYM,IBSP_KRCI_SYM,ABS_MAX_ELM,
     &                      WORK(KVOP),2*NASHT,2*NASHT,2*NASHT,2*NASHT,
     &                      NASH(1),NASH(2),MZ,IPQTOQ(1,0),NTEST)
C            GET_IAIB_DIRAC(IASP_KRCI_SYM,IBSP_KRCI_SYM,ELM_ABSMAX,
C    &                      AMAT,NROW,NCOL,LRQ,LCQ,JZ,IQP,NTEST)
                  
!
!       compute symmetry of operator IVOP and store in ISYMOPPRP_KRCI
        ISYMOPPRP_KRCI(I)=IDBGMULT(INVELM(IASP_KRCI_SYM),IBSP_KRCI_SYM)

        if((PRPNAM(INDXPR)(1:8).eq.'X dipole'.or.
     &      PRPNAM(INDXPR)(1:8).eq.'Y dipole').and.linear)then

!         # of unbarred > #   barred roots - ensure 'unbarred' x/y
          if(sym_xy.eq.2)then
            if(mod(ISYMOPPRP_KRCI(I),2).ne.0)then 
              ISYMOPPRP_KRCI(I) = ISYMOPPRP_KRCI(I) - 1
            end if
!         # of   barred > # unbarred roots - ensure '  barred' x/y
          else
            if(mod(ISYMOPPRP_KRCI(I),2).eq.0)then 
              ISYMOPPRP_KRCI(I) = ISYMOPPRP_KRCI(I) + 1
            end if
          end if
        end if

        if(ntest.gt.0)then
          print *, ' final double group irrep of property operator ==>',
     &               ISYMOPPRP_KRCI(I)
        end if
 20   CONTINUE
C
C     print section
C
      WRITE (LUPRI,'(/A,/A)') 
     &   '  Determination of operator representation in KRCI module',
     &   '  _______________________________________________________'
     
      WRITE (LUPRI,'(/A/A)')
     &'  operator     boson rep in DIRAC     boson rep in KR-CI',
     &'  ________     __________________     __________________'
C            ZDIPLEN             B1u                     Au
      DO I = 1, NPROP_KRCI
        INDXPR       = LPROP_KRCI(I)
        ISYMVOP      = IPRPSYM(INDXPR)
        ISYMVOP_KRCI = ISYMOPPRP_KRCI(I)
        IF(DOUGRP .eq. 5.or. 
     &     DOUGRP .eq. 6.or.
     &     DOUGRP .eq. 7)THEN
            IF( ISYMVOP_KRCI .eq. 1 )THEN
               REPVOP_KRCI = 'A  '
            ELSE
               REPVOP_KRCI = 'B  '
            END IF
        ELSE IF(DOUGRP .eq. 4)THEN
            IF( ISYMVOP_KRCI .eq. 1 )THEN
               REPVOP_KRCI = 'Ag '
            ELSE IF( ISYMVOP_KRCI .eq. 2 )THEN
               REPVOP_KRCI = 'Bg '
            ELSE IF( ISYMVOP_KRCI .eq. 3 )THEN
               REPVOP_KRCI = 'Au '
            ELSE 
               REPVOP_KRCI = 'Bu '
            END IF
        END IF
        if(linear) REPVOP_KRCI = REP(ISYMVOP-1)
        WRITE(LUPRI,'(A11,A15,A21)') 
     &        PRPNAM(INDXPR)(1:8), REP(ISYMVOP-1), REPVOP_KRCI(1:3)
      END DO
C
C     release memory
      CALL MEMREL('KRCIGETPRP',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('GET_PRPOPSYM_KRCI')
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine integral_type(ISOB,NGAS,ITPSV,LUINGR)
*
#include "implicit.h"
#include "ipoist8.inc"
      character*11 CELEC2
      character*13 CELEC
      dimension ISOB(NGAS,4)
      ISUM = 0
      do J=1,4,1
        do I=1,NGAS,1
          ISUM = ISUM + ISOB(I,J)
        end do
      end do
      if (ISUM.eq.2) CELEC = 'One-electron '
      if (ISUM.eq.4) CELEC = 'Two-electron '
      NAC = 0
      NBC = 0
      NAA = 0
      NBA = 0
      do I=1,NGAS,1
        NAC = NAC + ISOB(I,1)
        NBC = NBC + ISOB(I,2)
        NAA = NAA + ISOB(I,3)
        NBA = NBA + ISOB(I,4)
      end do
      if (NAC.eq.NAA.and.NBC.eq.NBA) CELEC2 = 'spin-free  '
      if (NAC.ne.NAA.or.NBC.ne.NBA)  CELEC2 = 'spin-orbit '
      if (ISUM.ne.ITPSV) 
     &     write(LUINGR,'(A,A,A)') CELEC,CELEC2,'integrals.'
      ITPSV = ISUM
      end

***********************************************************************
* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*  1-particle: molfdir full spinor format
*  2-particle: molfdir full spinor format
*
*--------------------------------------------------------------------*
*   November 12, 2001, Timo Fleig
*                - updated for correct one-particle function ordering
*                  Timo Fleig, June 2004
*                - Revised FSM/IUB ordering section             
*                  Timo Fleig, November 2006
*                - revised for proper NO-occupation number runs
*                  Stefan Knecht, October 2007
*                - CINC_NCT0 controls NCT = 0 option
*                  one-particle density run may fail otherwise
*                  Stefan Knecht, July 2008
*--------------------------------------------------------------------*
*
      subroutine trnsfdens(F1,F2,DENS,INDEX,LABEL,IHINDX,SIGN,IOP_REO,
     &                     SIGN_OPREO,ISF,NASHTL,IDEN_TYPE,IPRNT)
*
*  Input:  DENS
*  Output: F1,F2
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
#include "dgroup.h"
#include "dcborb.h"
#include "parluci.h"
#include "noccn_inf.inc"
*
      dimension F1(*),F2(*)
      dimension DENS(NALLINT),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,NALLINT),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),IOFFKF(2,MXPNGAS,2),
     &          IREOKF(MXPORB)
      integer*8 NINT_TOT, NINT_TP_I8, IRIOFF, nelm_cc
      integer   NINT_TP
      LOGICAL   CINC_NCT0
*
      NTESTL = 000
*
      CINC_NCT0 = .TRUE.
      IF ( IDEN_TYPE .eq. 1 ) THEN
         CINC_NCT0 = .FALSE.
         NTESTL = 00
      END IF
*      
      NTEST = max(IPRNT,NTESTL)
*
      NINT_TOT = N1ELINT + N2ELINT
      XNINT_TP_I8 = 0.0D0
*
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHTL2 = NASHTL * 2
*
      if (NTEST.ge.2) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================================='
        write(LUWRT,*) ' Transferring densities to MOLFDIR format.'
        write(LUWRT,*) '=========================================='
      end if
*
*  Fermion symmetry/Kramers reordering (for output density indexing)
*  KRMC wants - IFSYM loop
*              - unbarred/barred loop
*               - GAS loop
* Determine offset array
      IADDKF = 1
      do IFSM = 1,NFSYM
        do IUBT = 1,2
          do KGAS = 1,NGAS_DC
            IOFFKF(IFSM,KGAS,IUBT) = IADDKF
            IADDKF = IADDKF + NGSH(IFSM,KGAS)
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        do IUBT = 1,2
          if (IUBT.eq.1) write(LUWRT,*) ' IOFFKF for unbarred.'
          if (IUBT.eq.2) write(LUWRT,*) ' IOFFKF for   barred.'
          call iwrtmamn(IOFFKF(1,1,IUBT),NFSYM,NGAS_DC,2,MXPNGAS,LUWRT)
        end do
      end if
*
* Determine reordering array
      II = 1
      do IUBT = 1,2
        do IFSM = 1,NFSYM
          do KGAS = 1,NGAS_DC
            do IRN = 1,NGSH(IFSM,KGAS)
              IREOKF(II) = IOFFKF(IFSM,KGAS,IUBT) + IRN - 1
              II = II + 1
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        write(LUWRT,*) ' IREOKF array'
        call iwrtmamn(IREOKF,1,II-1,1,MXPORB,LUWRT)
      end if
*
* Process densities
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
*
        if (NTEST.ge.4) then
          write(6,*)
          if (IRI.eq.1) write(LUWRT,*) ' Real part     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary part: '
          write(6,*) ' -------------------'
        end if
*
        do IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(LUWRT,*)
            write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
            write(LUWRT,*)        ' -------------------------'
            write(LUWRT,*)
            write(LUWRT,*) ' LABEL array (for integral indices!):  '
            write(LUWRT,*) ' Internally, LABEL is reord. for operators.'
            call iwrtmamn(LABEL(1,IXTP),1,4,1,4,LUWRT)
            write(LUWRT,*) 'Input LUCIAREL operator order, output dens.'
            call iwrtmamn(IOP_REO(1,IXTP),1,4,1,4,LUWRT)
          end if
*
          NINT_TP_I8  = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8 = REAL(NINT_TP_I8)
          NINT_TP_X   = IGIVE_I_B(XNINT_TP_I8)
          NINT_TP     = NINT_TP_X
csk       write(LUWRT,*) ' NINT_TP is for IRUN',NINT_TP,IRUN
          do INT = 1,NINT_TP,1
*
*           determine type of density element
*
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
csk         write(luwrt,*) 'NCT is',NCT,'for IRUN',IRUN
*
            if (NCT.eq.2) then
              I = IREOKF(INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHTL)
              J = IREOKF(INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHTL)
*
              IOUT =  (IRI-1)*NASHTL2*NASHTL2 +
     &                 (J-1) *NASHTL2 +
     &                  I
              F1(IOUT) = DENS(IRUN+IRIOFF)
              if (NTEST.ge.10) then
                write(LUWRT,'(2I3,1F18.14,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
              end if
*
            else if (NCT.eq.4) then
              if (SIGN_OPREO(IXTP).eq.ONE.or.
     &            SIGN_OPREO(IXTP).eq.ONEM) then
*
*
*  INDEX array refers to integral indexing.
*  The corresponding second quantization operators then result
*  from the index reordering as given in GET_HX_RELA. 
*  Retaining integral indexing for densities as they come out of
*  the program, and labelled by the INDEX array; this has been
*  generated correctly within DIST_INTS_SIGDEN for all excitation 
*  types. An overall sign for all 2-particle densities is accounted
*  for in GET_HX_RELA (here called SIGN_OPREO):
                IND_R(1) = IREOKF(INDEX_L(1)-(LABEL(1,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator I '
                  write(LUWRT,*) ' INDEX_L(1) ',INDEX_L(1)
                  write(LUWRT,*) ' LABEL(1,IXTP) ',LABEL(1,IXTP)
                end if
*
                IND_R(2) = IREOKF(INDEX_L(2)-(LABEL(2,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator L '
                  write(LUWRT,*) ' INDEX_L(2) ',INDEX_L(2)
                  write(LUWRT,*) ' LABEL(2,IXTP) ',LABEL(2,IXTP)
                end if
*
                IND_R(3) = IREOKF(INDEX_L(3)-(LABEL(3,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator K '
                  write(LUWRT,*) ' INDEX_L(3) ',INDEX_L(3)
                  write(LUWRT,*) ' LABEL(3,IXTP) ',LABEL(3,IXTP)
                end if
*
                IND_R(4) = IREOKF(INDEX_L(4)-(LABEL(4,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator J '
                  write(LUWRT,*) ' INDEX_L(4) ',INDEX_L(4)
                  write(LUWRT,*) ' LABEL(4,IXTP) ',LABEL(4,IXTP)
                end if
*
                IOUT = (IRI-1)*NASHTL2*NASHTL2*NASHTL2*NASHTL2 +
     &                  (IND_R(4)-1) *NASHTL2*NASHTL2*NASHTL2 +
     &                  (IND_R(3)-1) *NASHTL2*NASHTL2 +
     &                  (IND_R(2)-1) *NASHTL2 +
     &                   IND_R(1)
                F2(IOUT) = DENS(IRUN+IRIOFF) * SIGN(IXTP) 
     &                                       * SIGN_OPREO(IXTP)
                if (NTEST.ge.10) then
                  write(LUWRT,*)
                  write(LUWRT,'(4I3,1F18.14,A,I4)') 
     &                      IND_R(1),IND_R(2),IND_R(3),IND_R(4),
     &                      F2(IOUT),' to pos.',IOUT
                end if
*
* ... and copy to the entry for allowed index permutation
* (ab|cd) = -(ad|cb)  in terms of integral indexing
*
                I_R = IND_R(1)
                L_R = IND_R(2)
                K_R = IND_R(3)
                J_R = IND_R(4)
                SIGN_CP = -1 * ONE
                IOUT_CP = (IRI-1) *NASHTL2*NASHTL2*NASHTL2*NASHTL2 +
     &                    (L_R-1) *NASHTL2*NASHTL2*NASHTL2 +
     &                    (K_R-1) *NASHTL2*NASHTL2 +
     &                    (J_R-1) *NASHTL2 +
     &                     I_R
                F2(IOUT_CP) = F2(IOUT) * SIGN_CP
                if (NTEST.ge.10) then
                  write(LUWRT,'(4I3,1F18.14,A,I4)') 
     &                      I_R,J_R,K_R,L_R,F2(IOUT_CP),
     &                      ' to pos.',IOUT_CP
                end if
*
              end if
*             ^ Non-redundant density type (see get_hx_rela)
            else if( (NCT.eq. 0) .and. .NOT. CINC_NCT0 )then
csk               write(LUWRT,*) ' what to do with 4 zero inidces???'
csk               write(LUWRT,*) ' continuing ...'
            else 
              write(LUWRT,*) 'No QFT implemented YET in LUCIAREL.'
              write(LUWRT,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting in trnsfdens.')
            end if
            IRUN = IRUN + 1
          end do ! integrals in excitation type
        end do ! Loop over excitation types
      end do ! Real/Imaginary
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*
*  1-particle: molfdir full spinor format
*
*--------------------------------------------------------------------*
*   November 12, 2001, Timo Fleig
*                - updated for correct one-particle function ordering
*                  Timo Fleig, June 2004
*                - Revised FSM/IUB ordering section             
*                  Timo Fleig, November 2006
*                - revised for proper NO-occupation number runs
*                  Stefan Knecht, October 2007
*      
*                - only 1-particle ordering, 
*                  Stefan Knecht - October 2007
*
*--------------------------------------------------------------------*
*
      subroutine trnsfdens_no(F1,DENS,INDEX,LABEL,IHINDX,SIGN,IOP_REO,
     &                     SIGN_OPREO,ISF,NASHTL,IPRNT)
*
*  Input:  DENS
*  Output: F1
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
#include "dgroup.h"
#include "dcborb.h"
#include "parluci.h"
#include "noccn_inf.inc"
*
      dimension F1(*)
      dimension DENS(NALLINT),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,N1ELINT),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),IOFFKF(2,MXPNGAS,2),
     &          IREOKF(MXPORB)
      integer*8 NINT_TOT, NINT_TP_I8, IRIOFF, nelm_cc
*
!     NTESTL = 10 ! debug
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
      NINT_TOT   = N1ELINT + N2ELINT
*
C
      CALL DZERO(F1,MIN(NZ,2)*(2*NASHT)**2 )
CSK   WRITE(LUWRT,*) ' F1 before transfer'
CSK   CALL WRTMATMN(F1,1,MIN(NZ,2)*(2*NASHT)**2,1,
CSK  &                   MIN(NZ,2)*(2*NASHT)**2,LUWRT)
C
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHTL2 = NASHTL * 2
*
      if (NTEST.ge.2) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================================='
        write(LUWRT,*) ' Transferring densities to MOLFDIR format.'
        write(LUWRT,*) '=========================================='
      end if
*
*  Fermion symmetry/Kramers reordering (for output density indexing)
*  KRMC wants - IFSYM loop
*              - unbarred/barred loop
*               - GAS loop
* Determine offset array
      IADDKF = 1
      do IFSM = 1,NFSYM
        do IUBT = 1,2
          do KGAS = 1,NGAS_DC
            IOFFKF(IFSM,KGAS,IUBT) = IADDKF
            IADDKF = IADDKF + NGSH(IFSM,KGAS)
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        do IUBT = 1,2
          if (IUBT.eq.1) write(LUWRT,*) ' IOFFKF for unbarred.'
          if (IUBT.eq.2) write(LUWRT,*) ' IOFFKF for   barred.'
          call iwrtmamn(IOFFKF(1,1,IUBT),NFSYM,NGAS_DC,2,MXPNGAS,LUWRT)
        end do
      end if
*
* Determine reordering array
      II = 1
      do IUBT = 1,2
        do IFSM = 1,NFSYM
          do KGAS = 1,NGAS_DC
            do IRN = 1,NGSH(IFSM,KGAS)
              IREOKF(II) = IOFFKF(IFSM,KGAS,IUBT) + IRN - 1
              II = II + 1
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        write(LUWRT,*) ' IREOKF array'
        call iwrtmamn(IREOKF,1,II-1,1,MXPORB,LUWRT)
      end if
*
* Process densities
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
*
        if (NTEST.ge.4) then
          write(6,*)
          if (IRI.eq.1) write(LUWRT,*) ' Real part     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary part: '
          write(6,*) ' -------------------'
        end if
*
*
        do IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(LUWRT,*)
            write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
            write(LUWRT,*)        ' -------------------------'
            write(LUWRT,*)
            write(LUWRT,*) ' LABEL array (for integral indices!):  '
            write(LUWRT,*) ' Internally, LABEL is reord. for operators.'
            call iwrtmamn(LABEL(1,IXTP),1,4,1,4,LUWRT)
            write(LUWRT,*) 'Input LUCIAREL operator order, output dens.'
            call iwrtmamn(IOP_REO(1,IXTP),1,4,1,4,LUWRT)
          end if
*
          NINT_TP_I8 = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8 = REAL(NINT_TP_I8)
          NINT_TP_X = IGIVE_I_B(XNINT_TP_I8)
          NINT_TP = NINT_TP_X
csk          write(LUWRT,*) ' NINT_TP is for IRUN and iRI',NINT_TP,IRUN,IRI
*
          IF( IRI .eq. 1 .and. IRUN .eq. (N1ELINT + 1) ) THEN
*
*           real part finished
*
            GOTO 222
*           
          ELSE IF( IRI. eq. 2 .and. IRUN .eq. ( N1ELINT + 1) ) THEN
*
*           well, we have finished the 1-particle part 
*
            GOTO 333
*
          END IF
          do INT = 1,NINT_TP,1
*
* Determine type of density element
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
csk            write(luwrt,*) 'NCT is',NCT,'for IRUN',IRUN
*
            if (NCT.eq.2) then
*
              I = IREOKF(INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHTL)
              J = IREOKF(INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHTL)
*
              IOUT =  (IRI-1)*NASHTL2*NASHTL2 +
     &                 (J-1) *NASHTL2 +
     &                  I
              F1(IOUT) = DENS(IRUN+IRIOFF)
              if (NTEST.ge.10) then
                write(LUWRT,'(2I3,1F18.14,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
                write(LUWRT,*) ' original position IRUN+IRIOFF',
     &                           IRUN+IRIOFF
              end if
*
            else if( NCT .eq. 0 .or. NCT .eq. 4)then
csk               write(LUWRT,*) ' what to do with 4 zero inidces???'
csk               write(LUWRT,*) ' continuing ...'
csk               write(LUWRT,*) ' what to do with 4 nonzero inidces???'
csk               write(LUWRT,*) ' continuing ...'
            else 
              write(LUWRT,*) 'No QFT implemented YET in LUCIAREL.'
              write(LUWRT,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting in trnsfdens_no.')
            end if
            IRUN = IRUN + 1
          end do
*         ^ Loop over integrals in excitation type
        end do
*       ^ Loop over excitation types
 222    CONTINUE
      end do
*     ^ loop over real/imaginary
*
 333  return
      end
#if defined (VAR_MPI2)
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*
*  1-particle: molfdir full spinor format
*
*--------------------------------------------------------------------*
*   November 12, 2001, Timo Fleig
*                - updated for correct one-particle function ordering
*                  Timo Fleig, June 2004
*                - Revised FSM/IUB ordering section             
*                  Timo Fleig, November 2006
*                - revised for proper NO-occupation number runs
*                  Stefan Knecht, October 2007
*      
*                - only 1-particle ordering, 
*                  Stefan Knecht - October 2007
*
*                - shared memory mode: Stefan Knecht - February 2008
*
*--------------------------------------------------------------------*
*
      subroutine trnsfdens_no_sm(F1,DENS,IT_TTPL,IT_TTOL,INDEX,
     &                           LABEL,IHINDX,SIGN,IOP_REO,
     &                           SIGN_OPREO,ISF,NASHTL,IPRNT)
*
*  Input:  DENS
*  Output: F1
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
#include "krmc_shmem.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
#include "dgroup.h"
#include "dcborb.h"
#include "noccn_inf.inc"

*
      dimension F1(*), IT_TTPL(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IT_TTOL(*)
      dimension DENS(NALLINT),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,N1ELINT),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),IOFFKF(2,MXPNGAS,2),
     &          IREOKF(MXPORB)
      INTEGER*8 NINT_TOT, NINT_TP_I8, IRIOFF, ITB_I8
      INTEGER*8 IRUN_XTP, IOFF_GET , nelm_cc
      INTEGER ID_GET
*
      IRUN_XTP = 0
      IOFF_GET = 0
      ID_GET   = 0
      ITB_I8   = 0
      NTESTL   = 00
      NTEST = max(IPRNT,NTESTL)
*
      NINT_TOT   = N1ELINT + N2ELINT
*
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHTL2 = NASHTL * 2
*
      if (NTEST.ge.2) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================================='
        write(LUWRT,*) ' Transferring densities to MOLFDIR format.'
        write(LUWRT,*) '=========================================='
      end if
*
*  Fermion symmetry/Kramers reordering (for output density indexing)
*  KRMC wants - IFSYM loop
*              - unbarred/barred loop
*               - GAS loop
* Determine offset array
      IADDKF = 1
      do IFSM = 1,NFSYM
        do IUBT = 1,2
          do KGAS = 1,NGAS_DC
            IOFFKF(IFSM,KGAS,IUBT) = IADDKF
            IADDKF = IADDKF + NGSH(IFSM,KGAS)
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        do IUBT = 1,2
          if (IUBT.eq.1) write(LUWRT,*) ' IOFFKF for unbarred.'
          if (IUBT.eq.2) write(LUWRT,*) ' IOFFKF for   barred.'
          call iwrtmamn(IOFFKF(1,1,IUBT),NFSYM,NGAS_DC,2,MXPNGAS,LUWRT)
        end do
      end if
*
* Determine reordering array
      II = 1
      do IUBT = 1,2
        do IFSM = 1,NFSYM
          do KGAS = 1,NGAS_DC
            do IRN = 1,NGSH(IFSM,KGAS)
              IREOKF(II) = IOFFKF(IFSM,KGAS,IUBT) + IRN - 1
              II = II + 1
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        write(LUWRT,*) ' IREOKF array'
        call iwrtmamn(IREOKF,1,II-1,1,MXPORB,LUWRT)
      end if
*
* Process densities
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
*
        if (NTEST.ge.4) then
          write(6,*)
          if (IRI.eq.1) write(LUWRT,*) ' Real part     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary part: '
          write(6,*) ' -------------------'
        end if
*
*
        do IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(LUWRT,*)
            write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
            write(LUWRT,*)        ' -------------------------'
            write(LUWRT,*)
            write(LUWRT,*) ' LABEL array (for integral indices!):  '
            write(LUWRT,*) ' Internally, LABEL is reord. for operators.'
            call iwrtmamn(LABEL(1,IXTP),1,4,1,4,LUWRT)
            write(LUWRT,*) 'Input LUCIAREL operator order, output dens.'
            call iwrtmamn(IOP_REO(1,IXTP),1,4,1,4,LUWRT)
          end if
*
          NINT_TP_I8 = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8 = REAL(NINT_TP_I8)
          NINT_TP_X = IGIVE_I_B(XNINT_TP_I8)
          NINT_TP = NINT_TP_X
csk          write(LUWRT,*) ' NINT_TP is for IRUN and iRI',NINT_TP,IRUN,IRI
*
          IF( IRI .eq. 1 .and. IRUN .eq. (N1ELINT + 1) ) THEN
*
*           real part finished
*
            GOTO 222
*           
          ELSE IF( IRI. eq. 2 .and. IRUN .eq. ( N1ELINT + 1) ) THEN
*
*           well, we have finished the 1-particle part 
*
            GOTO 333
*
          END IF
          IRUN_XTP = 0
          do INT = 1,NINT_TP,1
*
* Determine type of density element
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
csk            write(luwrt,*) 'NCT is',NCT,'for IRUN',IRUN
*
            if (NCT.eq.2) then
*
              I = IREOKF(INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHTL)
              J = IREOKF(INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHTL)
*
              IOUT =  (IRI-1)*NASHTL2*NASHTL2 +
     &                 (J-1) *NASHTL2 +
     &                  I
csk              F1(IOUT) = DENS(IRUN+IRIOFF)
              IXXTTP = IXTP
              IF( IRI .eq. 2 ) THEN
                IXXTTP = ITTP + NSPOBEX_TP
              END IF
              IF( IRUN_XTP .eq. 0 )THEN
C               process id and window offset (complex case included)
                ID_GET = IT_TTPL( IXXTTP )
C               check for availability
                IF( ID_GET .lt. 0) THEN
                  IRUN = IRUN + NINT_TP
                  GOTO 155
                END IF
C               window offset
                ITB_I8 = IT_TTOL( IXXTTP )
              END IF
              IOFF_GET = ITB_I8 + IRUN_XTP
              YDENS = 0.0D0
              IF( MYNEW_ID_SM .eq. ID_GET )THEN
                YDENS = DENS(IOFF_GET+1)
              ELSE
                CALL MPIXGET(YDENS,1,2,ID_GET,
     &                       IOFF_GET,1,2,MY_T_WIN)
              END IF
              F1(IOUT) = YDENS
              if (NTEST.ge.10) then
                write(LUWRT,'(2I3,1F18.14,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
                write(LUWRT,*) ' original position IRUN+IRIOFF',
     &                           IRUN+IRIOFF, IOFF_GET
              end if
*
            else if( NCT .eq. 0 .or. NCT .eq. 4)then
csk               write(LUWRT,*) ' what to do with 4 zero inidces???'
csk               write(LUWRT,*) ' continuing ...'
csk               write(LUWRT,*) ' what to do with 4 nonzero inidces???'
csk               write(LUWRT,*) ' continuing ...'
            else 
              write(LUWRT,*) 'No QFT implemented YET in LUCIAREL.'
              write(LUWRT,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting in trnsfdens_no_sm.')
            end if
            IRUN_XTP = IRUN_XTP + 1
            IRUN = IRUN + 1
          end do
*         ^ Loop over integrals in excitation type
 155      CONTINUE
        end do
*       ^ Loop over excitation types
 222    CONTINUE
      end do
*     ^ Loop over Real/Imaginary
*
 333  CONTINUE
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*
*  1-particle: molfdir full spinor format
*
*--------------------------------------------------------------------*
*   November 12, 2001, Timo Fleig
*                - updated for correct one-particle function ordering
*                  Timo Fleig, June 2004
*                - Revised FSM/IUB ordering section             
*                  Timo Fleig, November 2006
*                - revised for proper NO-occupation number runs
*                  Stefan Knecht, October 2007
*      
*                - only 1-particle ordering, 
*                  Stefan Knecht - October 2007
*
*--------------------------------------------------------------------*
*
      subroutine trnsfdens_no_split(F1,DENS,DENS_SCR,IT_TTPL,INDEX,
     &                              LABEL,IHINDX,SIGN,
     &                              IOP_REO,SIGN_OPREO,ISF,NASHTL,
     &                              IPRNT)
*
*  Input:  DENS
*  Output: F1
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
#include "dgroup.h"
#include "dcborb.h"
#include "parluci.h"
#include "noccn_inf.inc"
#include "krmc_shmem.h"
*
      dimension F1(*), DENS_SCR(*), IT_TTPL(*)
      dimension DENS(*),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,N1ELINT),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),IOFFKF(2,MXPNGAS,2),
     &          IREOKF(MXPORB)
      integer*8 NINT_TOT, ITB_BUFF, ITB_BUFF_ADD, nelm_cc
*
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
      NINT_TOT     = N1ELINT + N2ELINT
      ITB_BUFF_ADD = 0
      ITB_BUFF     = 0
*
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHTL2 = NASHTL * 2
*
      if (NTEST.ge.2) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================================='
        write(LUWRT,*) ' Transferring densities to MOLFDIR format.'
        write(LUWRT,*) '=========================================='
      end if
*
*  Fermion symmetry/Kramers reordering (for output density indexing)
*  KRMC wants - IFSYM loop
*              - unbarred/barred loop
*               - GAS loop
* Determine offset array
      IADDKF = 1
      do IFSM = 1,NFSYM
        do IUBT = 1,2
          do KGAS = 1,NGAS_DC
            IOFFKF(IFSM,KGAS,IUBT) = IADDKF
            IADDKF = IADDKF + NGSH(IFSM,KGAS)
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        do IUBT = 1,2
          if (IUBT.eq.1) write(LUWRT,*) ' IOFFKF for unbarred.'
          if (IUBT.eq.2) write(LUWRT,*) ' IOFFKF for   barred.'
          call iwrtmamn(IOFFKF(1,1,IUBT),NFSYM,NGAS_DC,2,MXPNGAS,LUWRT)
        end do
      end if
*
* Determine reordering array
      II = 1
      do IUBT = 1,2
        do IFSM = 1,NFSYM
          do KGAS = 1,NGAS_DC
            do IRN = 1,NGSH(IFSM,KGAS)
              IREOKF(II) = IOFFKF(IFSM,KGAS,IUBT) + IRN - 1
              II = II + 1
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        write(LUWRT,*) ' IREOKF array'
        call iwrtmamn(IREOKF,1,II-1,1,MXPORB,LUWRT)
      end if
*
* Process densities
      do IRI = 1,IRILP,1
        IRUN = 1
        ITB_BUFF = 1
        do JRI = 1,IRI-1,1
          ITB_BUFF = ITB_BUFF + LEN_T_BUFF
        end do
csk     write(LUWRT,*) '  ITB_BUFF, LEN_T_BUFF',ITB_BUFF,LEN_T_BUFF

*
        if (NTEST.ge.4) then
          write(LUWRT,*)
          if (IRI.eq.1) write(LUWRT,*) ' Real part     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary part: '
          write(LUWRT,*) ' -------------------'
        end if
*
*
        do 111 IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(LUWRT,*)
            write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
            write(LUWRT,*)        ' -------------------------'
            write(LUWRT,*)
            write(LUWRT,*) ' LABEL array (for integral indices!):  '
            write(LUWRT,*) ' Internally, LABEL is reord. for operators.'
            call iwrtmamn(LABEL(1,IXTP),1,4,1,4,LUWRT)
            write(LUWRT,*) 'Input LUCIAREL operator order, output dens.'
            call iwrtmamn(IOP_REO(1,IXTP),1,4,1,4,LUWRT)
          end if
*
          IXXTTP = (( IRI - 1 ) * NSPOBEX_TP) + IXTP
C         active T block? (complex case included)
          ID_GET = IT_TTPL( IXXTTP )
csk       WRITE(6,*) ' type IXXTTP and ID_GET',IXXTTP,ID_GET,MYPROC
          ITB_BUFF_ADD = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8  = REAL(ITB_BUFF_ADD)
          NINT_TP      = IGIVE_I_B(XNINT_TP_I8)
          IF( ID_GET .lt. 0 ) THEN
            IRUN         = IRUN + NINT_TP
            ITB_BUFF_ADD = 0
          ELSE
C           receive or copy DENS block - at present: MASTER only
            IF( MYPROC .eq. MASTER )THEN
              CALL DZERO(DENS_SCR,NINT_TP)
              IF( ID_GET .eq. 0 )THEN
                CALL DCOPY(NINT_TP,DENS(ITB_BUFF),1,DENS_SCR,1)
csk             WRITE(LUWRT,*) ' density block at ITB_BUFF'
csk             CALL WRTMATMN(DENS_SCR,1,NINT_TP,1,NINT_TP,LUWRT)
              ELSE
                call interface_mpi_RECV_r1_work_f77(DENS_SCR,NINT_TP,
     &                        ID_GET,821,
     &                        global_communicator)
                ITB_BUFF_ADD = 0
              END IF
            ELSE
              IRUN         = IRUN + NINT_TP
              IF( MYPROC .eq. ID_GET )THEN
                call interface_mpi_SEND(DENS(ITB_BUFF),NINT_TP,
     &                        MASTER,821,global_communicator)
                ITB_BUFF     = ITB_BUFF + ITB_BUFF_ADD
              ELSE
                ITB_BUFF_ADD = 0
              END IF
            END IF
          END IF
*
*
          IF( IRI .eq. 1 .and. IRUN .eq. (N1ELINT + 1) ) THEN
*
*           real part finished
*
            GOTO 222
*           
          ELSE IF( IRI. eq. 2 .and. IRUN .eq. ( N1ELINT + 1) ) THEN
*
*           well, we have finished the 1-particle part 
*
            GOTO 333
*
          END IF
C
          IF( MYPROC .ne. MASTER ) THEN
              GOTO 111
          ELSE
              IF( ID_GET .lt. 0 ) GOTO 111
          END IF
C
          IOFF_SET = 1
C
          do INT = 1,NINT_TP,1
*
*           Determine type of density element
*
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
            if (NCT.eq.2) then
*
              I = IREOKF(INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHTL)
              J = IREOKF(INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHTL)
*
              IOUT =  (IRI-1)*NASHTL2*NASHTL2 +
     &                 (J-1) *NASHTL2 +
     &                  I
              F1(IOUT) = DENS_SCR(IOFF_SET)
              if (NTEST.ge.10) then
                write(LUWRT,'(2I3,1F18.14,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
                write(LUWRT,*) ' original position IOFF_SET',
     &                           IOFF_SET
              end if
*
            else if( NCT .eq. 0 .or. NCT .eq. 4)then
csk               write(LUWRT,*) ' what to do with 4 zero inidces???'
csk               write(LUWRT,*) ' continuing ...'
csk               write(LUWRT,*) ' what to do with 4 nonzero inidces???'
csk               write(LUWRT,*) ' continuing ...'
            else 
              write(LUWRT,*) 'No QFT implemented YET in LUCIAREL.'
              write(LUWRT,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting in trnsfdens_no.')
            end if
            IRUN     = IRUN     + 1
            IOFF_SET = IOFF_SET + 1
          end do
*         ^ Loop over integrals in excitation type
          ITB_BUFF     = ITB_BUFF + ITB_BUFF_ADD
 111    CONTINUE
*       ^ Loop over excitation types
 222    CONTINUE
      end do
*     ^ Loop over Real/Imaginary
*
 333  CONTINUE
*
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*  1-particle: molfdir full spinor format
*  2-particle: molfdir full spinor format
*
*--------------------------------------------------------------------*
*   November 12, 2001, Timo Fleig
*                - updated for correct one-particle function ordering
*                  Timo Fleig, June 2004
*                - Revised FSM/IUB ordering section             
*                  Timo Fleig, November 2006
*                - T split version                             
*                  Stefan Knecht, May 2008
*--------------------------------------------------------------------*
*
      subroutine trnsfdens_split(F1,F2,DENS,DENS_SCR,IT_TTPL,
     &                           INDEX,LABEL,IHINDX,SIGN,IOP_REO,
     &                           SIGN_OPREO,ISF,NASHTL,IPRNT)
*
*  Input:  DENS
*  Output: F1,F2
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
#include "dgroup.h"
#include "dcborb.h"
#include "parluci.h"
#include "noccn_inf.inc"
#include "krmc_shmem.h"
*
      dimension F1(*),F2(*), DENS_SCR(*), IT_TTPL(*)
      dimension DENS(*),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,NALLINT),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),IOFFKF(2,MXPNGAS,2),
     &          IREOKF(MXPORB)
      integer*8 NINT_TOT, ITB_BUFF, ITB_BUFF_ADD, nelm_cc
      integer   NINT_TP
*
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
      NINT_TOT     = N1ELINT + N2ELINT
      ITB_BUFF_ADD = 0
      ITB_BUFF     = 0
*
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHTL2 = NASHTL * 2
*
      if (NTEST.ge.2) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================================='
        write(LUWRT,*) ' Transferring densities to MOLFDIR format.'
        write(LUWRT,*) '=========================================='
      end if
*
*  Fermion symmetry/Kramers reordering (for output density indexing)
*  KRMC wants - IFSYM loop
*              - unbarred/barred loop
*               - GAS loop
* Determine offset array
      IADDKF = 1
      do IFSM = 1,NFSYM
        do IUBT = 1,2
          do KGAS = 1,NGAS_DC
            IOFFKF(IFSM,KGAS,IUBT) = IADDKF
            IADDKF = IADDKF + NGSH(IFSM,KGAS)
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        do IUBT = 1,2
          if (IUBT.eq.1) write(LUWRT,*) ' IOFFKF for unbarred.'
          if (IUBT.eq.2) write(LUWRT,*) ' IOFFKF for   barred.'
          call iwrtmamn(IOFFKF(1,1,IUBT),NFSYM,NGAS_DC,2,MXPNGAS,LUWRT)
        end do
      end if
*
* Determine reordering array
      II = 1
      do IUBT = 1,2
        do IFSM = 1,NFSYM
          do KGAS = 1,NGAS_DC
            do IRN = 1,NGSH(IFSM,KGAS)
              IREOKF(II) = IOFFKF(IFSM,KGAS,IUBT) + IRN - 1
              II = II + 1
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        write(LUWRT,*) ' IREOKF array'
        call iwrtmamn(IREOKF,1,II-1,1,MXPORB,LUWRT)
      end if
*
* Process densities
      do IRI = 1,IRILP,1
        IRUN     = 1
        ITB_BUFF = 1
        do JRI = 1,IRI-1,1
          ITB_BUFF = ITB_BUFF + LEN_T_BUFF
        end do
*
        if (NTEST.ge.4) then
          write(LUWRT,*)
          if (IRI.eq.1) write(LUWRT,*) ' Real part     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary part: '
          write(LUWRT,*) ' -------------------'
        end if
*
        do 111 IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(LUWRT,*)
            write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
            write(LUWRT,*)        ' -------------------------'
            write(LUWRT,*)
            write(LUWRT,*) ' LABEL array (for integral indices!):  '
            write(LUWRT,*) ' Internally, LABEL is reord. for operators.'
            call iwrtmamn(LABEL(1,IXTP),1,4,1,4,LUWRT)
            write(LUWRT,*) 'Input LUCIAREL operator order, output dens.'
            call iwrtmamn(IOP_REO(1,IXTP),1,4,1,4,LUWRT)
          end if
*
          IXXTTP = (( IRI - 1 ) * NSPOBEX_TP) + IXTP
C         active T block? (complex case included)
          ID_GET = IT_TTPL( IXXTTP )
csk       WRITE(6,*) ' type IXXTTP and ID_GET',IXXTTP,ID_GET,MYPROC
          ITB_BUFF_ADD = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8  = REAL(ITB_BUFF_ADD)
          NINT_TP      = IGIVE_I_B(XNINT_TP_I8)
          IF( ID_GET .lt. 0 ) THEN
            IRUN         = IRUN + NINT_TP
            ITB_BUFF_ADD = 0
            GOTO 111
          ELSE
C           receive or copy DENS block - at present: MASTER only
            IF( MYPROC .eq. MASTER )THEN
              CALL DZERO(DENS_SCR,NINT_TP)
              IF( ID_GET .eq. 0 )THEN
                CALL DCOPY(NINT_TP,DENS(ITB_BUFF),1,DENS_SCR,1)
              ELSE
                call interface_mpi_RECV_r1_work_f77(DENS_SCR,NINT_TP,
     &                        ID_GET,821,
     &                        global_communicator)
                ITB_BUFF_ADD = 0
              END IF
            ELSE
              IRUN         = IRUN + NINT_TP
              IF( MYPROC .eq. ID_GET )THEN
                call interface_mpi_SEND(DENS(ITB_BUFF),NINT_TP,
     &                        MASTER,821,global_communicator)
                ITB_BUFF     = ITB_BUFF + ITB_BUFF_ADD
              ELSE
                ITB_BUFF_ADD = 0
              END IF
            END IF
          END IF
C
C         MASTER enters loop over integrals
C
          IF( MYPROC .ne. MASTER ) GOTO 111
C
          IOFF_SET = 1
C
          do INT = 1,NINT_TP,1
*
*           Determine type of density element
*
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
*
            if (NCT.eq.2) then
              I = IREOKF(INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHTL)
              J = IREOKF(INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHTL)
*
              IOUT =  (IRI-1)*NASHTL2*NASHTL2 +
     &                 (J-1) *NASHTL2 +
     &                  I
              F1(IOUT) = DENS_SCR(IOFF_SET)
              if (NTEST.ge.10) then
                write(LUWRT,'(2I3,1F18.14,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
              end if
*
            else if (NCT.eq.4) then
              if (SIGN_OPREO(IXTP).eq.ONE.or.
     &            SIGN_OPREO(IXTP).eq.ONEM) then
*
*
*  INDEX array refers to integral indexing.
*  The corresponding second quantization operators then result
*  from the index reordering as given in GET_HX_RELA. 
*  Retaining integral indexing for densities as they come out of
*  the program, and labelled by the INDEX array; this has been
*  generated correctly within DIST_INTS_SIGDEN for all excitation 
*  types. An overall sign for all 2-particle densities is accounted
*  for in GET_HX_RELA (here called SIGN_OPREO):
                IND_R(1) = IREOKF(INDEX_L(1)-(LABEL(1,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator I '
                  write(LUWRT,*) ' INDEX_L(1) ',INDEX_L(1)
                  write(LUWRT,*) ' LABEL(1,IXTP) ',LABEL(1,IXTP)
                end if
*
                IND_R(2) = IREOKF(INDEX_L(2)-(LABEL(2,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator L '
                  write(LUWRT,*) ' INDEX_L(2) ',INDEX_L(2)
                  write(LUWRT,*) ' LABEL(2,IXTP) ',LABEL(2,IXTP)
                end if
*
                IND_R(3) = IREOKF(INDEX_L(3)-(LABEL(3,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator K '
                  write(LUWRT,*) ' INDEX_L(3) ',INDEX_L(3)
                  write(LUWRT,*) ' LABEL(3,IXTP) ',LABEL(3,IXTP)
                end if
*
                IND_R(4) = IREOKF(INDEX_L(4)-(LABEL(4,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator J '
                  write(LUWRT,*) ' INDEX_L(4) ',INDEX_L(4)
                  write(LUWRT,*) ' LABEL(4,IXTP) ',LABEL(4,IXTP)
                end if
*
                IOUT = (IRI-1)*NASHTL2*NASHTL2*NASHTL2*NASHTL2 +
     &                  (IND_R(4)-1) *NASHTL2*NASHTL2*NASHTL2 +
     &                  (IND_R(3)-1) *NASHTL2*NASHTL2 +
     &                  (IND_R(2)-1) *NASHTL2 +
     &                   IND_R(1)
                F2(IOUT) = DENS_SCR(IOFF_SET) * SIGN(IXTP) 
     &                                        * SIGN_OPREO(IXTP)
                if (NTEST.ge.10) then
                  write(LUWRT,*)
                  write(LUWRT,'(4I3,1F18.14,A,I4)') 
     &                      IND_R(1),IND_R(2),IND_R(3),IND_R(4),
     &                      F2(IOUT),' to pos.',IOUT
                end if
*
* ... and copy to the entry for allowed index permutation
* (ab|cd) = -(ad|cb)  in terms of integral indexing
*
                I_R = IND_R(1)
                L_R = IND_R(2)
                K_R = IND_R(3)
                J_R = IND_R(4)
                SIGN_CP = -1 * ONE
                IOUT_CP = (IRI-1) *NASHTL2*NASHTL2*NASHTL2*NASHTL2 +
     &                    (L_R-1) *NASHTL2*NASHTL2*NASHTL2 +
     &                    (K_R-1) *NASHTL2*NASHTL2 +
     &                    (J_R-1) *NASHTL2 +
     &                     I_R
                F2(IOUT_CP) = F2(IOUT) * SIGN_CP
                if (NTEST.ge.10) then
                  write(LUWRT,'(4I3,1F18.14,A,I4)') 
     &                      I_R,J_R,K_R,L_R,F2(IOUT_CP),
     &                      ' to pos.',IOUT_CP
                end if
*
              end if
*             ^ Non-redundant density type (see get_hx_rela)
            else 
              write(LUWRT,*) 'No QFT implemented YET in LUCIAREL.'
              write(LUWRT,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting in trnsfdens.')
            end if
            IRUN = IRUN + 1
            IOFF_SET = IOFF_SET + 1
          end do
*         ^ Loop over integrals in excitation type
        ITB_BUFF     = ITB_BUFF + ITB_BUFF_ADD
 111    CONTINUE
*       ^ Loop over excitation types
      end do
*     ^ Loop over Real/Imaginary
*
      return
      end
#endif
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Write density matrix for current eigenvector in 
* SIGDEN ordering.
*
*   October 18, 2001, Timo Fleig
*
*   Revised for operator ordering, Mar 20, 2002
*
      subroutine wrtdens(DENS,INDEX,IHINDX,SIGN,NASHT,LABEL,IOP_REO,
     &                   SIGN_REO,ISF)
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
#include "parluci.h"
*
      dimension DENS(NALLINT),SIGN(NSPOBEX_TP)
      dimension INDEX(4,NALLINT),IND_R(4),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),LABEL(4,NSPOBEX_TP),
     &          SIGN_REO(NSPOBEX_TP)
      integer*8 NINT_TOT, NINT_TP_I8, IRIOFF , nelm_cc
*
      NINT_TOT = N1ELINT + N2ELINT
      if (ISF.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
        write(LUWRT,*) '  '
        write(LUWRT,*) '  '
        write(LUWRT,'(2X,A)') 
     &            'Writing densities in reordered form (for MCSCF):'
        write(LUWRT,*)'  - a+ a b b+  '
        write(LUWRT,'(2X,A)') 
     &            '(means sign multiplied to all 2-particle densities)'
        write(LUWRT,*)'  (done in GET_HX_RELA)'
        write(LUWRT,*) '  '
        if (IRI.eq.1) write(LUWRT,*) ' Real density elements     : '
        if (IRI.eq.2) write(LUWRT,*) ' Imaginary density elements: '
        write(LUWRT,*) ' -------------------------------'
        do IXTP = 1,NSPOBEX_TP,1
          write(LUWRT,*) '  '
          write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
          write(LUWRT,*)        ' -------------------------'
          NINT_TP_I8 = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8 = REAL(NINT_TP_I8)
          NINT_TP_X = IGIVE_I_B(XNINT_TP_I8)
          NINT_TP = NINT_TP_X
csk          write(LUWRT,*) ' NINT_TP is for IRUN',NINT_TP,IRUN
          do INT = 1,NINT_TP,1
*
            if (INDEX(1,IRUN).ne.0.and.INDEX(2,IRUN).ne.0) then
*  Two-particle density element.
*
*  INDEX array refers to integral indexing.
              IUBFAC = 1
              if (INDEX(1,IRUN).lt.0) IUBFAC = -1
              IND_R(1) = INDEX(1,IRUN) - (LABEL(1,IXTP) - 1)/2 * NASHT
*
              IUBFAC = 1
              if (INDEX(2,IRUN).lt.0) IUBFAC = -1
              IND_R(2) = INDEX(2,IRUN) - (LABEL(2,IXTP) - 1)/2 * NASHT
*
              IUBFAC = 1
              if (INDEX(3,IRUN).lt.0) IUBFAC = -1
              IND_R(3) = INDEX(3,IRUN) - (LABEL(3,IXTP) - 1)/2 * NASHT
*
              IUBFAC = 1
              if (INDEX(4,IRUN).lt.0) IUBFAC = -1
              IND_R(4) = INDEX(4,IRUN) - (LABEL(4,IXTP) - 1)/2 * NASHT
            else
*  One-particle density element. No modification.
              write(LUWRT,*) ' One-particle density element for IRUN',
     &                         IRUN
              IND_R(1) = INDEX(1,IRUN)
              IND_R(2) = INDEX(2,IRUN)
              IND_R(3) = INDEX(3,IRUN)
              IND_R(4) = INDEX(4,IRUN)
            end if
            write(LUWRT,'(1X,4I3,2X,1F18.10)') 
     &           (IND_R(I),I=1,4,1),DENS(IRUN+IRIOFF) * SIGN(IXTP)
     &                                                * SIGN_REO(IXTP)
            IRUN = IRUN + 1
          end do
        end do
      end do
c     stop 'AFTER WRTDENS.'
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
      subroutine LCR_RD_MRCMDC(TINTR,TINTI,INDK,INDL,RKLR,RKLI)
*
#include "implicit.h"
#include "mxpdim.inc"
#include "clunit.inc"
*
      character*4 REPNA(14),REPN(14)
      character*7 FNAM1
      character*6 FNAM2,FNAM3
      character*8 TIMEX
      character*10 DATEX,DATEXB,TIMEXB(8)
      character*10 NEXTFIL
      logical BREIT,SPINFR,REALAR
*
      dimension MULTB(MXPIRR*2,MXPIRR*2)
      dimension IRPMO(MXPORB),IRPAMO(MXPORB),EPS(MXPORB)
      dimension TINTR(*),TINTI(*)
      dimension KR(-MXPORB:MXPORB)
      dimension INDK(*),INDL(*)
      dimension RKLR(*),RKLI(*)
*
      FNAM1 = 'MRCONEE'
      FNAM2 = 'MDCINT'
*
* Read one-electron integrals
*
      open (LUSC96,FILE=FNAM1,FORM='UNFORMATTED')
       print*,'File ',FNAM1,' opened.'
      rewind (LUSC96)
      read (LUSC96) NORB,BREIT,ECORE,NFSYM,NZ,SPINFR
       print*,'NORB,BREIT,ECORE,NFSYM: ',NORB,BREIT,ECORE,NFSYM
       print*,' (NORB is the number of spinors.)'
       print*,'NZ,SPINFR: ',NZ,SPINFR
      REALAR = NZ.eq.1.or.SPINFR
       print*,'REALAR ',REALAR
      read (LUSC96) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
       print*,'NSYMRP: ',NSYMRP
       print*,'REPN: ',(REPN(IRP),IRP=1,NSYMRP)
      read (LUSC96) NREP,(REPNA(IRP),IRP=1,2*NREP)
       print*,'NREP: ',NREP
       print*,'REPNA: ',(REPNA(IRP),IRP=1,2*NREP)
      read (LUSC96) ((MULTB(I,J),I=1,2*NREP),J=1,2*NREP)
       call iwrtma(MULTB(1,1),NREP*2,NREP*2,MXPIRR*2,MXPIRR*2)
      read (LUSC96) (IRPMO(IMO),IRPAMO(IMO),EPS(IMO),IMO=1,NORB)
c      print*,'IRPMO: ',(IRPMO(I),I=1,NORB)
c      print*,'IRPAMO: ',(IRPAMO(I),I=1,NORB)
c      print*,'EPS: ',(EPS(I),I=1,NORB)
      if (REALAR) then
         read (LUSC96) (TINTR(JI),TDUM,JI=1,NORB*NORB)
      else
         read (LUSC96) (TINTR(JI),TINTI(JI),JI=1,NORB*NORB)
      end if
       do IP=1,NORB**2
        print*,'TINTR(',IP,') = ',TINTR(IP)
       end do
      close (LUSC96)
*
* Read two-electron integrals
*
      open (LUSC97,file=FNAM2,form='UNFORMATTED')
       print*,'File ',FNAM2,' opened.'
      rewind (LUSC97)
      read (LUSC97,err=10000,end=10000) DATEX,TIMEX,NKR,
     &                                  (KR(I),KR(-I),I=1,NKR)
       print*,'DATEX,TIMEX,NKR ',DATEX,TIMEX,NKR
       print*,'KR(NKR)..KR(-NKR) ',(KR(I),I=-NKR,NKR)
      IF (2*NKR.NE.NORB)
     &   CALL QUIT('MRCONEE AND MDCINT ARE INCOMPATIBLE')
      IF (BREIT) THEN
        call quit ('BREIT interaction not yet implemented.')
c       open (MDBINT,file=FNAM(3),form='UNFORMATTED')
c       read (mdbint,err=10001,end=10001) datexb,timexb,nkrb
c       IF (NKR.NE.NKRB) CALL QUIT('MDCINT AND MDBINT ARE INCOMPATIBLE')
      END IF
      IF (REALAR) THEN
c        DO I = 1, M1
c           RKLI(I) = 0.0
c        ENDDO
      ENDIF

c     IF (REALAR) THEN
c        read (LUSC97,ERR=10010,END=10010) ikr,jkr,nz,
c    &                (indk(inz),indl(inz),inz=1,nz),
c    &                (rklr(inz),inz=1,nz)
c     ELSE
c        read (LUSC97,ERR=10010,END=10010) ikr,jkr,nz,
c    &                (indk(inz),indl(inz),inz=1,nz),
c    &                (rklr(inz),rkli(inz),inz=1,nz)
c     END IF
       print*,'IKR,JKR,NZ ',IKR,JKR,NZ
       print*,'INDK : ',(INDK(I),I=1,NZ)
       print*,'INDL : ',(INDL(I),I=1,NZ)
       print*,'RKLR : ',(RKLR(I),I=1,NZ)



      close (LUSC97)
       print*,'File ',FNAM2,' closed.'


      return
*
10000 CALL QUIT('ERROR READING HEADER OF MDCINT')
10001 CALL QUIT('ERROR READING HEADER OF MDBINT')
10010 CALL QUIT('ERROR READING INTEGRALS FROM MDCINT')
10011 CALL QUIT('ERROR READING INTEGRALS FROM MDBINT')
*
      end
