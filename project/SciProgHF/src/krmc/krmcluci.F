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

! define task symbols for CALL DIRAC_PARCTL( task )
#include "dirac_partask.h"
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_SETUP()
C***********************************************************************
C
C     Set up starting information for LUCIAREL
C      - assign disk units
C      - determinant types
C      - string types
C      - symmetry information
C      - integral dimensioning
C
C
C     Written by J. Thyssen and T. Fleig         - Jan 12 2001
C     Last revision :
C
C***********************************************************************
      use memory_allocator
      use mospinor_info
      use symmetry_setup_krci
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "../krmc/dcbbos.h"
#include "../luciarel/cprnt.inc"
#include "../luciarel/mxpdim.inc"
#include "../luciarel/cgas.inc"
#include "dcborb.h"
#include "dgroup.h"
#include "parluci.h"
!     scratch
      integer, allocatable  :: orb_sym_vec(:)
      integer               :: orb_sym_fh
      character (len=11)    :: orb_sym_fn
      character (len= 4)    :: mynum_str
C
!
      is_ana_mcscf = .false.

C     set "I/O node master" + "shared memory node master"
      N_MASTER    = 0
      N_MASTER_SM = 0
 
!     assign disk units
      CALL DISKUN_REL

!     setup double group multiplication table
      call symmetry_setup_init(dougrp,nirr_dg)

!     initialize mo-spinor information arrays
      call mospinor_info_init()
      
!     place spinors (as MOs) in appropriate boson irreps  for further processing (see next step)
      if(dougrp .ge. 10)then
        call num2str(myproc,mynum_str)
        call alloc(orb_sym_vec,norbt)
        write(orb_sym_fn,'(a7,a4)') 'orbsym.',mynum_str
        orb_sym_vec = 0
        orb_sym_fh = 99
        open(orb_sym_fh,file=orb_sym_fn,status='old',
     &       form='unformatted',action='readwrite',position='rewind')
        read(orb_sym_fh) (orb_sym_vec(i), i=1,norbt)
        close(orb_sym_fh)
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
        print *, 'content of orb_sym_vec is: ',myproc
        call iwrtma(orb_sym_vec,1,norbt,1,norbt)
#endif
!#undef LUCI_DEBUG
      end if

!     generalized routine for handling of all double groups including
!     linear symmetry (dougrp == 10 or 11)
!     write(lupri,*) 'bla - go in here???? '
#if defined VAR_IFORT && !defined INT_STAR8
!     write(lupri,*) 'bla - go in here, yes!!!! '
      if(myproc .eq. master)then
        luxxx = lupri
      else
        luxxx = luwrt
      end if
      write(luxxx,*) ' '
      call flshfo(luxxx)
#endif
      call match_spinor_2_boson_irrep(ngsob,ngsob2,ibosym,ngsh,
     &                                orb_sym_vec,iorb,npsh,nish,
     &                                mxndgirr,mxpngas,nfsym,ngas,
     &                                dougrp,nirr_dg,mj2rep,
     &                                imosp_dirac_counter1,
     &                                imosp_dirac_counter2,
     &                                imosp_dirac_mjub,
     &                                imosp_dirac_mjb)

      if(dougrp .ge. 10) call dealloc(orb_sym_vec)

!     setup shell, orbital and spinor index lists 
!     (in short: fill arrays stored in module mospinor_info)
      CALL ORBINF_REL(LUPRI,IPRORB)

!     types of determinants to be included
      CALL DETTYP(IPRSTR)

#if defined VAR_IFORT && !defined INT_STAR8
      write(luxxx,*) ' '
      call flshfo(luxxx)
#endif
!     number of string types
      CALL STRTYP_GAS_REL(IPRSTR)

!     number of integrals
      CALL INTDIM_REL(IPRTRA)
C
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_CLOSEDOWN()
!**********************************************************************
!     purpose: 
!
!       - close LUCIAREL files and keep them
!       - destroy arrays containing symmetry information
!
!     written by Stefan Knecht   April 24, 2007
!
!**********************************************************************
      use mospinor_info
      use symmetry_setup_krci
      IMPLICIT REAL*8           (A-H,O-Z)
#include "../luciarel/clunit.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#endif
#include "parluci.h"
#include "krciprop.h"
*
!     close files
      close(unit=ludia,status='KEEP')
      if(myproc.eq.master.and..not.runxprop)
     &   close(unit=luc,status='KEEP')
      if(myproc.eq.master) close(unit=lu_info,status='KEEP')
      close(unit=luhc,  status='KEEP')
      close(unit=lusc1, status='KEEP')
      close(unit=lusc2, status='KEEP')
      close(unit=lusc3, status='KEEP')
      close(unit=lusc34,status='KEEP')
      close(unit=lusc35,status='KEEP')
      close(unit=lusc36,status='KEEP')
      close(unit=lusc37,status='KEEP')
      close(unit=lusc38,status='KEEP')
      close(unit=lusc39,status='KEEP')
      close(unit=lusc41,status='KEEP')
      close(unit=lusc61,status='DELETE')
      close(unit=lusc62,status='DELETE')

!     deallocate mo-spinor information arrays
      call mospinor_info_delete()

!     deallocate symmetry arrays
      call symmetry_setup_delete()
     
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRMC_LUCI(F1,F2,WRK,LWRK,ECORE_KRMC,CIRUN)
C***********************************************************************
C
C
C     L U C I A R E L
C
C     parallel driver routine,
C     written by Stefan Knecht - April 2007
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
#if defined (VAR_MPI)
#include "infpar.h"
#endif
#include "parluci.h"
      CHARACTER*6 CIRUN
C     initial definitions
#if defined (VAR_MPI)
      MASTER = MPARID
      MYPROC = MYTID
C     Add the master node
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
C
C     Summon the co-workers, who are waiting in the general menu routine.
      CALL LUMSTI_REL
C
C     communicate necessary information with co-workers ...
#if defined (VAR_MPI2)
      CALL SYNC_NODES_REL_P1(ECORE_KRMC,CIRUN,LUPRI)
C
      IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'SIGMA ' .or.
     &    CIRUN .eq. 'SIGMAD' .or. CIRUN .eq. 'IJKLRO' .or. 
     &    CIRUN .eq. 'KR-CI ' ) THEN
C
        CALL SYNC_NODES_REL_P2(F1,F2,CIRUN)
      ELSE IF( CIRUN .eq. 'PROP1 ' )THEN
        CALL SYNC_NODES_REL_XPROP
      END IF
C
C     Transfer data to LUCIAREL common blocks (only for NODES) 
C     --> "1" means that communication is enabled.
      CALL TRKRLUCI(CIRUN,1)
#endif
C
C     Transfer correct output file handle to common block
      LUWRT = LUPRI
C
C     Enter the generic LUCIAREL routine 
C     (master becomes now somewhat like primus inter parens)
      CALL KRMC_LUCI2(F1,F2,WRK,LWRK,ECORE_KRMC,CIRUN)
C
C     Return the co-workers to the general menu routine.
      CALL LUMSTE_REL
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRMC_LUCI2(F1,F2,WRK,LWRK,ECORE_KRMC,CIRUN)
C***********************************************************************
*
*
* L U C I A R E L
*
*
* CI for program for :FCI
*                     RASCI
*                     GAS GAS GAS GAS GAS GAS
*
* Written by Jeppe Olsen , winter of 1991
*                          GAS version in action summer of 95
*
* Implementation of relativistic double group CI,
*                          Timo Fleig + Jeppe Olsen
*                          1995 - 1999
*
* parallel adaption, Stefan Knecht, April 2007 - ?
*                  
*
* Simplified and modified from luciarel.f
* Execute task given by task parameter CIRUN
*
C***********************************************************************
*
      use mc_energies
      use symmetry_setup_krci
      use os_utils
      use interface_to_mpi
      IMPLICIT REAL*8(A-H,O-Z)
*. Parameters for dimensioning
#include "../luciarel/mxpdim.inc"
*. Boson symmetry info
#include "maxorb.h"
#include "dcbbos.h"
#include "dgroup.h"
*.Memory
#include "../luciarel/wrkspc.inc"
#include "../luciarel/ipoist8.inc"
*.File numbers
#include "../luciarel/clunit.inc"
*.Print flags
#include "../luciarel/cprnt.inc"
#include "../luciarel/cstate.inc"
#include "../luciarel/crun.inc"
#include "../luciarel/cicisp.inc"
#include "../luciarel/oper.inc"
#include "../luciarel/cgas.inc"
#include "../luciarel/glbbas.inc"
#include "../luciarel/cintfo.inc"
#include "../luciarel/ctcc.inc"
C information for property run
#include "krciprop.h"
      COMMON/CECORER/ECORE
CTF
      character*6 CIRUN
      CHARACTER*255 MACHINENAME
      INTEGER NAMELENGTH, MZ_NUM
      logical MCSCF, RUN_SHMEM_T
CTF  global reference values for Z_BLKFO needed in MV7 !!
      common/glbref/ISPC,ISM,ITERGLB
      DIMENSION F1(*),F2(*),WRK(LWRK)
#include "infpar.h"
#include "parluci.h"
      integer global_scratch_disk
      real(8) :: eref = 0.0d0
!     integer(8) :: kbla_blubb
C
C*******************************************
      CALL QENTER('KRMC_LUCI2')
C*******************************************
C
!     transfer core energy to common block
      ECORE = ECORE_KRMC
#ifdef LUCI_DEBUG
      if(myproc .eq. master) call print_ci_info()
#endif

      MZ_NUM = 0
      MZ_NUM = MIN(NZ,2)
      MIN_ONE = - 1

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
      call memmar(K_OFFSET,KADD,'INI   ',IDUMMY,'DUMMY ')

!     write(luwrt,*) 'K_OFFSET ==> ',K_OFFSET
!     IDUM = 0
!     CALL MEMMAR(KDUM,IDUM,'MARK ',IDUM,'GROUPL')
!     call memmar(kbla_blubb,4,'ADDS  ',1,'BLABLU')
!     call izero(work(kbla_blubb),4)
C     ... no memchecks as default
      IMEMCK = 0
      NTEST = 0
      MCSCF = .true.
C
C     debugging print level for parallel runs -- control from outside?
C
      NPTEST_VAR = 0
C
C     Internal string information (stored in WORK, bases in /STRBAS/)
      CALL STRINF_GAS_REL(WORK(1),IPRSTR)
C
C     Internal subspaces
      CALL LCISPC_REL(LCSBLK,IPRCIX,CIRUN)
C
******************************************
      if (CIRUN.eq.'NDET  ') GO TO 9999
******************************************
      if (CIRUN.eq.'RSTRMC') then
        ISPC = 1
        CALL Z_BLKFO_REL(IDC,NMS2VAL,ISPC,IREFSM,KPCLBT,KPCLEBT,
     &                   KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK,
     &                   NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2,
     &                   IDUMMY,.FALSE.,2,0,IDUMMY)
        GO TO 9999
      end if
*******************************************
C     OK, we want to do something, not just calculate number of
C     determinants.
*******************************************
C
C     Prepare for parallelization:
C     setting up groups of processors (stored in /LUCIPARREL/)
C
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'MARK ',IDUM,'GROUPL')
C
      CALL MEMMAR(KPROCLIST,  NMPROC,'ADDS  ',1,'IHOSTL')
!     call RMEMCHK_opt('bla bla check 1     ')
      CALL MEMMAR(KGROUPLIST, NMPROC,'ADDS  ',1,'IHOST2')
      CALL MEMMAR(KPROCLENGTH,NMPROC,'ADDS  ',1,'INAMEL')
      CALL MEMMAR(KPROCLIST2, NMPROC,'ADDS  ',1,'IHOST3')
C
      CALL ISETVC(WORK(KPROCLIST),  MIN_ONE,NMPROC)
      CALL ISETVC(WORK(KGROUPLIST), MIN_ONE,NMPROC)
      CALL ISETVC(WORK(KPROCLENGTH),MIN_ONE,NMPROC)
      CALL ISETVC(WORK(KPROCLIST2), MIN_ONE,NMPROC)
C     initialize number of process groups sharing a c-vector file
      NFLGRPS_REL = 0
C
#if defined (VAR_MPI2)
C
C     hardwired to 1, otherwise it fails...
C
      IPARALLELIO = 1
C
C
      global_scratch_disk = get_environment_integer('GLBSCR',
     &global_scratch_disk)

      if(global_scratch_disk .ne. 0 .or. global_scratch_disk.ne. 1)then
!       fall back to default: 0
        global_scratch_disk = 0
      end if

      call interface_mpi_bcast(global_scratch_disk,1,
     &               MASTER,global_communicator)


      if(global_scratch_disk .eq. 0)then
        NAMELENGTH = 0
C
        call interface_mpi_GET_PROCESSOR_NAME(MACHINENAME,NAMELENGTH)

        CALL FIND_GROUP_OF_PROCS_REL(MACHINENAME,NAMELENGTH,
     &                               WORK(KPROCLIST),
     &                               WORK(KPROCLENGTH))
C
C       number of groups sharing a file ( = NFLGRPS ) is determined
C       in the subroutine above. If there is only one group
C       than we can simulate a shared filesystem.
        IF( NFLGRPS_REL .eq. 1 ) THEN 
          CALL ISETVC(WORK(KPROCLIST),1,NMPROC)
        END IF
      else
!       in case of a parallel filesystem, we should use it
        NFLGRPS_REL = 1
        CALL ISETVC(WORK(KPROCLIST),1,NMPROC)
      end if
!
C     no integral broadcast necessary
      IF( NFLGRPS_REL .eq. 1 ) IIOMOD_REL = 1
C
C     MPI "shared memory" mode
C
      IF( SHARED_M )THEN
C
        IF( LEVEL_SM .lt. 6 ) THEN
C
C          level 0 and 5: build node-groups for shared memory usage
C
           NAMELENGTH = 0
C
           call interface_mpi_GET_PROCESSOR_NAME(MACHINENAME,NAMELENGTH)
C
           CALL ISETVC(WORK(KPROCLENGTH),-1,NMPROC)
C
           CALL FIND_GROUP_OF_PROCS_REL_SM(MACHINENAME,NAMELENGTH,
     &                        WORK(KPROCLIST2),WORK(KPROCLENGTH))
C
        ELSE
C
C          level 6:   modelling global shared memory 
C                     NOTE: you need a fast inter-node connect,
C                         e.g. InfiniBand, Myrinet, etc ...
C
           CALL ISETVC(WORK(KPROCLIST2),1,NMPROC)
C         
        END IF
C
      ELSE
C
         CALL ISETVC(WORK(KPROCLIST2),1,NMPROC)
C
      END IF
C
C     set-up all useful MPI stuff: communication groups, MPI-files ...
C                                  all important variables are stored on 
C                                  common block /LUPARGROUP/
C
C     construct new communication groups
      CALL GROUP_CONSTRUCTOR_REL(WORK(KGROUPLIST),WORK(KPROCLIST),
     &                           WORK(KPROCLIST2))
C
      NPTESTVAR = 00
      IF( NPTESTVAR .ge. 20 ) THEN
        WRITE(LUWRT,*) 'this is my old and new rank ',myproc, MYNEW_ID
        WRITE(LUWRT,*) 'this is my "shared memory" rank ', MYNEW_ID_SM
        WRITE(LUWRT,*) 'this is my file grouplist:'
        CALL IWRTMAMN(WORK(KGROUPLIST),1,NMPROC,1,NMPROC,LUWRT)
        WRITE(LUWRT,*) 'this is my complete proclist:'
        CALL IWRTMAMN(WORK(KPROCLIST),1,NMPROC,1,NMPROC,LUWRT)
        WRITE(LUWRT,*) 'this is my shared memory proclist:'
        CALL IWRTMAMN(WORK(KPROCLIST2),1,NMPROC,1,NMPROC,LUWRT)
      END IF
      NPTESTVAR = 0
C
C     open the files for parallel file I/O
C
      CALL SETUNITS_PAR_OPEN_REL(MY_GROUPN)
C
#else
      CALL ISETVC(WORK(KPROCLIST),1,NMPROC)
      NFLGRPS_REL = 0
#endif
C     Finished with preparation for parallelization.

C     IDGSTYP determines the integral classes and corresponding operator
C     classes depending on the double group in use.
      IDGSTYP = -1
      if(PNTGRP.eq.6.or.PNTGRP.eq.5) then
        IDGSTYP = 3
      else if (PNTGRP.ge.7)then
        IDGSTYP = 2
      end if
C
C     operator type
C
      if (IDGSTYP.le.2) then
        IHTYPE = 4
      else if (IDGSTYP.eq.3) then
        IHTYPE = 5
      else
        WRITE(LUWRT,*) 'Unallowed value of IDGSTYP : ',IDGSTYP
        CALL QUIT('*** KRMC_LUCI2: unallowed value of IDGSTYP ***')
      end if
C
C     property calculation - we split here from the "usual CIRUN" route.
C     optimal (what we think - sk + hjaaj) loop structure:
C
C     loop iroot = 1, nroot
C        loop ipropT = 1, npropT
C                              
C           |sigma_iroot> = ipropT|C_iroot>
C            
C           loop jroot = 1, nroot
C              PROPMAT(jroot,iroot,ipropT,MZ) = <sigma_iroot|C_jroot>
C           end loop
C        end loop
C     end loop
C             
      CXPROPRUN = .FALSE.
C
      IF( CIRUN .eq. 'PROP1 ') THEN 
         CXPROPRUN = .TRUE.
         CALL XPROP_KRCI(F2,IHTYPE,WORK(KPROCLIST),WORK(KGROUPLIST),
     &                   CIRUN,mz_num)
         GOTO 8888
      END IF
C
C     loop over operator
C     __________________
C
      NLOOPOP_T = 1
C
      DO II = 1, NLOOPOP_T
C
C       set operator symmetry (Hamiltonian will be called 'T' 
C       throughout the KR-CI module)
C       default: totally symmetric operator
        ISYM_T      = 1
C
C       ... no two-electron excitations
        NO_TTWO = 0
        IF( CIRUN .eq. 'DENS1 ') NO_TTWO = 1
C
C
        CALL SET_HOP_DBG(IHTYPE,NO_TTWO,ISYM_T,IPRHAM_CI)
C
C       no integrals - return to main routine 
        IF( NOINT .eq. 1 ) GO TO 8888
C
C       loop over natural spinors (currently hardwired to 1)
        do INATIT=1,NATITER,1
C
          RUN_SHMEM_T = .FALSE.
C
          IF( IT_SHL .ge. 0 .or. SPLIT_IJKL ) THEN
C           ... read-in comes later in GASCI_REL
            RUN_SHMEM_T = .TRUE.
          END IF
          IF( REORD_IJKL ) RUN_SHMEM_T = .TRUE.
C         ... reorder always at this point: RUN_SHMEM_T = .FALSE.
          IF( CIRUN .eq. 'IJKLRO' ) RUN_SHMEM_T = .FALSE.
C
************************************************************
          IF (CIRUN /= 'ANALYZ' .and. CIRUN /= 'QCORR ' .and. 
     &        CIRUN /= 'REFVEC')THEN
CSK         IDUMMY ==> IT_TTPL
            call picasso(IPRTRA,MCSCF,
     &                   NL2D,IBU2D,ISM2D,F1,F2,DUMMY1,DUMMY2,
     &                   IDUMMY,ibosym,CIRUN,RUN_SHMEM_T)
          END IF
************************************************************
          IF ( CIRUN .eq. 'IJKLRO' ) GOTO 300
C
          ISTOSPC = 0
          IF(IRESTR.EQ.1) ISTOSPC = 1
C
C         Restarted Lambda calculations need special attention
C         since restart is realized in second calculation
C         First calculation is used to establish H0
C
          IRESTR_ORIG=IRESTR
          IF(IRESTR.EQ.1.AND.XLAMBDA.NE.1.0D0) THEN
C
            write(LUWRT,*) 'THIS ROUTE IS STILL FORBIDDEN !!'
            write(LUWRT,'(A16,I2,1F3.1)')
     &            'IRESTR,XLAMBDA =',IRESTR,XLAMBDA
            Call Abend2( 'Quitting.' )
C
          END IF
C         ^ End of special handling of restarted calc with
C           lambda modified op.
C
          ISKIPEI_INI = ISKIPEI
C         note: NCMBSPC is hardwired to 1
          DO JCMBSPC = 1, NCMBSPC
            ISKIPEI = 0
C
C           Perform GASCI/GASPT calculation for internal space ICISPC
C
            IF (NTEST.ge.1) THEN
              WRITE(LUWRT,'(/A)')
     &        ' =============================================='
              WRITE(LUWRT,'(A,I3)')
     &        ' Information on calculations in space ', JCMBSPC
              WRITE(LUWRT,'(A)')
     &        ' =============================================='
              WRITE(LUWRT,'(/A,I3)')
     &        ' Number of calculation in this CI space ',
     &          NSEQCI(JCMBSPC)
            END IF
C
            I_EXPAND = 1
            IF(XLAMBDA.NE.1.0D0 .AND.JCMBSPC.GT.1) THEN
              WRITE(LUWRT,*) ' =================================='
              WRITE(LUWRT,*) '   Modified operator will be used'
              WRITE(LUWRT,*) ' =================================='
              IF(JCMBSPC.EQ.2) THEN
C               Modify the operator
                WRITE(LUWRT,*) ' Operator will be modified '
                WRITE(LUWRT,*) ' This route is still forbidden! '
                Call Abend2( 'Quitting.' )
              END IF
            END IF
            IF(JCMBSPC.EQ.2.AND.IRESTR_ORIG.EQ.1
     &                     .AND.XLAMBDA.NE.1.0D0)  THEN
              WRITE(LUWRT,*) ' Restart vectors will be copied to LUC'
              WRITE(LUWRT,*)
     &              ' CI will restart with vectors from prev. calc'
              IF(ICISTR.EQ.1) THEN
                LBLK = XISPSM(IREFSM,2)
              ELSE
                LBLK = -1
              END IF
              call memmar(IDUM,IDUM,'MARK ',IDUM,'VCSAVE')
              IF(ICISTR.EQ.1) THEN
                LBLOCK = XISPSM(IREFSM,2)
              ELSE IF (ICISTR.EQ.2) THEN
                LBLOCK = MXSB
              ELSE IF (ICISTR.EQ.3) THEN
                LBLOCK = MXSOOB
              END IF
              call memmar(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
              IF( MYPROC .eq. MASTER ) CALL REWINE(LUC,-1)
                CALL REWINE(LUSC39,-1)
              DO JROOT = 1, NROOT
                CALL COPVCD_REL(LUSC39,LUC,WORK(KVEC1),0,LBLK)
              END DO
              call memmar(IDUM,IDUM,'FLUSM',IDUM,'VCSAVE')
C             No expansion should follow
              I_EXPAND = 0
C             But normal restart
              IRESTR = 1
            END IF
C
            DO JSEQ = 1,NSEQCI(JCMBSPC)
              IF((ISEQCI(JSEQ,JCMBSPC).GT.0
     &            .OR.ISEQCI(JSEQ,JCMBSPC).EQ.-5)
     &            .AND. JCMBSPC.NE.1.AND.ISTOSPC.NE.JCMBSPC
     &            .AND.I_EXPAND.EQ.1)                   THEN
C               Restart from previous spaces 
C               ( Assuming a progressing sequence :
C               spaces are just added, not subtracted )
                LUIN = LUC
                LUOUT = LUSC1
                IF(ICISTR.EQ.1) THEN
                  LBLK = XISPSM(IREFSM,JCMBSPC)
                ELSE
                  LBLK = -1
                END IF
                CALL EXPCIV_REL(IREFSM,ISTOSPC,LUIN,JCMBSPC,LUOUT,
     &                          LBLK,LUSC2,
     &                          NROOT,1,IDC,IPRDIA)
C               last space where vectors were stored
                ISTOSPC = JCMBSPC
                ISKIPEI = ISKIPEI_INI
C               Expanded vector will be used as initial vector in the
C               zero space calculation. Tell next CI to restart from
C               CI vectors
                IRESTR = 1
              END IF
              IF(ISEQCI(JSEQ,JCMBSPC).GT.0) THEN
C
C               Good old normal CI !!!!
C               do CI in space JCMBPSC
                MAXIT = ISEQCI(JSEQ,JCMBSPC)
C
                ISPC = JCMBSPC
                ISM  = IREFSM
                eref = emc_ref ! initialize...
C
C###############################################################
                CALL GASCI_REL(F1,F2,
     &                         ISM,ISPC,IPRDIA,
     &                         EREF,ECORE_KRMC,CIRUN,
     &                         WORK(KPROCLIST),WORK(KGROUPLIST))
C###############################################################
C
                if(cirun == 'CIINII' .or. cirun == 'KR-CI ') 
     &          eci_init = emc_ref ! transfer new value to the MCSCF world outside
C               last space where vectors were stored
                ISTOSPC = JCMBSPC
              END IF
            END DO
*
          END DO
*         ^ end loop over CI spaces
 300      CONTINUE
C
        END DO
C       ^ end loop over natural spinor iterations (default=1)
      END DO
C     ^ end loop over different T operators (default=1)
C
C     release memory
 8888 CALL MEMMAR(IDUM,IDUM,'FLUSM',IDUM,'GROUPL')
C
#if defined (VAR_MPI2)
C
C     free all MPI-2 groups
C
      IF(CIRUN .ne. 'KR-CI ') CALL SETUNITS_PAR_CLOSE_REL(CIRUN,0)
      CALL GROUP_DESTRUCTOR_REL(MYNEW_COMM,MYNEW_COMM_SM,ICOMM,
     &                          MYNEW_COMM_SM_C)
C
#endif
C
 9999 CALL QEXIT('KRMC_LUCI2')

C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_TASK(CIRUN)
C***********************************************************************
C
C     Introduce LUCIAREL action.
C
C     Written by  T. Fleig         - Nov 09 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      character*6 CIRUN
*
      NTEST = 1
*
      if (NTEST.gt.1) then
        write(LUPRI,*)
        write(LUPRI,*) ' ////////////////////////////////// '
        write(LUPRI,*) ' //'
      end if
      if (NTEST.gt.0) then
        write(LUPRI,*) ' //  LUCIAREL called'
        write(LUPRI,*) ' //  CIRUN = ',CIRUN
      end if
*
      if (NTEST.gt.0) then
      if (CIRUN.eq.'NDET  ') then
        write(LUPRI,*)
     &  ' //  The number of CI determinants will be computed'
      else if (CIRUN.eq.'CIINII') then
        write(LUPRI,*) ' //  Determining configurational start guess'
      else if (CIRUN.eq.'KR-CI ') then
        write(LUPRI,*) ' //  Running large-scale CI calculation'
      else if (CIRUN.eq.'RSTRMC') then
        write(LUPRI,*) ' //  Determining info for restart of MCSCF.'
      else if (CIRUN.eq.'SIGMA ') then
        write(LUPRI,*) ' //  Computing a sigma vector.'
      else if (CIRUN.eq.'SIGMAD') then
        write(LUPRI,*) ' //  Computing a sigma vector.'
        write(LUPRI,*) ' //      Saving H diagonal also.'
      else if (CIRUN.eq.'DENS1 ') then
        write(LUPRI,*) ' //  Computing only 1-particle density matrices'
      else if (CIRUN.eq.'DENS2 ') then
        write(LUPRI,*) ' //  Computing 1- and 2-particle dens. matrices'
      else if (CIRUN.eq.'ANALYZ') then
        write(LUPRI,*) ' //  Analyzing the current CI vector(s).'
      else if (CIRUN.eq.'IJKLRO') then
        write(LUPRI,*) ' //  Reordering all (ij|kl) according to GAS'
        write(LUPRI,*) ' //    scheme. Saving to file IJKL_REOD.'
      else if (CIRUN.eq.'DIAG  ') then
        write(LUPRI,*) ' //  Computing CI diagonal.'
      else if (CIRUN.eq.'PROP1 ') then
        write(LUPRI,*) ' //  Computing one-electron properties.'
      end if
      end if

      if (NTEST.gt.1) then
        write(LUPRI,*) ' //'
        write(LUPRI,*) ' ////////////////////////////////// '
      end if

      if (CIRUN.eq.'DIAG  ') then
        call quit('NOT IMPLEMENTED YET!')
      end if
*
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_NDET(NDET,ISYM,WORK,LWORK)
C***********************************************************************
C
C     Calculate number of determinants.
C
C     Written by J. Thyssen and T. Fleig         - Jan 12 2001
C     Last revision : 23.01.2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
C
#include "../luciarel/mxpdim.inc"
#include "../luciarel/cicisp.inc"
C
      CHARACTER*6 CIRUN
C
      DIMENSION WORK(*)
C
      CIRUN = 'NDET  '
      ISYM_KRMC = ISYM
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Compute the number of determinants for this setup and return.
      CALL KRMC_LUCI(DUMMY,DUMMY,WORK,LWORK,DUMMY,CIRUN)
C
C     transfer number of dets to common block variable
      NDET = INT(XISPSM(ISYM_KRMC,1))
C
C     close files and release symmetry arrays
      CALL LUCI_CLOSEDOWN()
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_RCIST(CREF,CMO,WORK,LWORK)
C***********************************************************************
C
C     Driver routine for LUCIAREL CI module.
C
C     Written by S. Knecht - Aug 2008
C
C     Last revision :
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbopt.h"
      DIMENSION CREF(*), CMO(*), WORK(*)
C
C     which module is calling?
C
      IF( JKRRUNTYPE .eq. 0 ) THEN
         WRITE(LUPRI,'(A,I3)') ' *** ERROR in LUCI_RCIST ***'//
     &   ' LUCI_RCIST not in use for 2nd order SCF!',JKRRUNTYPE
         CALL QUIT('*** ERROR in LUCI_RCIST ***')
      ELSE IF ( JKRRUNTYPE .eq. 1 ) THEN
C
C        get CI start guess for KR-MCSCF wave function optimization
C
         CALL LUCI_CIINII(CREF,CMO,WORK,LWORK)
      ELSE IF ( JKRRUNTYPE .eq. 2) THEN
C
C        run large-scale KR-CI calculation
C
         IF( MAXCIT .gt. 0 ) CALL LUCI_KRCI(CMO,WORK,LWORK)
C
      ELSE
         WRITE(LUPRI,'(A,I18)') ' *** ERROR in LUCI_RCIST ***'//
     &   ' wrong KR-run type in LUCI_RCIST !',JKRRUNTYPE
         CALL QUIT('*** ERROR in LUCI_RCIST ***')
      END IF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_KRCI(CMO,WORK,LWORK)
C***********************************************************************
C
C     Run large-scale CI calculations. Called from module KR-CI.
C
C     Based on LUCI_CIINII.
C
C     Written by S. Knecht - Aug 2008
C
C     Last revision :
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
      CHARACTER*6 CIRUN,RCISTFILE
      LOGICAL PARCAL_SAVE
C
      DIMENSION WORK(*),CMO(*)
C
#include "../luciarel/clunit.inc"
#include "parluci.h"
C
#include "memint.h"
C
      CIRUN = 'KR-CI '
      MZ = MIN(NZ,2)
C
      KFRSAV = KFREE
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
#if defined (VAR_MPI2)
C
C     test fock matrix work for slaves
C
      ISET_PARCAL = 0
      IF( ISET_PARCAL .eq. 1 )THEN
        PARCAL_SAVE = .TRUE.
        CALL TEST_PAR_CALC(PARCAL_SAVE,0)
      END IF
#endif
C
C     Flag for type of integral interface:
C     I_IITP = 1: Old 4IND1XXXX0 driven interface
C     I_IITP = 2: New MRCONEE/MDCINT driven interface
C
      N_FCACM = (2*NASHT)*(2*NASHT)
C
      I_IITP = 1
C
      if (I_IITP.eq.2) then
C
        print*,' Testing new integral interface (MRCONEE/MDCINT)'
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDK,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDL,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRKLR,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRKLI,NZ,WORK,KFREE,LFREE)
        print*,'Allocated for ',N_FCACM,' integrals, MZ=',MZ
        CALL LCR_RD_MRCMDC(WORK(KFCACM),WORK(KFCACM+N_FCACM),
     &                     WORK(KINDK),WORK(KINDL),
     &                     WORK(KRKLR),WORK(KRKLI))
        CALL MEMGET('REAL',KH2AC,NASHT*NASHT*NNASHX*NZ*3,
     &               WORK,KFREE,LFREE)
        print*,'Dimension of MOLFDIR array : ',(2*NASHT)**4 * 2
        print*,'Dimension of (NZ,3)  array : ',NASHT*NASHT*NNASHX*NZ*3

        stop 'FIXME: this interface is not complete'
        call memget('INTE',kibeig,norbt,work,kfree,lfree)
        call izero(work(kibeig),norbt)

        call rgeth2(dummy,work(kh2ac),dummy,work(kibeig),.false.,
     &              .true.,.true.,work(kfree),lfree)
C
      else if (I_IITP.eq.1) then
C
C       Calculate core Hamiltonian:
C       ---------------------------
C
C       The core Hamiltonian is the active-active part of FC,
C       a.k.a. FCAC.
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
        IF( CRDFO_MAT ) THEN
C
C         read FCmo from file LUKRM3
C
          CALL KRMC_GETFCK(WORK(KFC),IPROPT)
C
          ECORE = ECORE_LR
C
        ELSE
C
C         Calculate FCmo ( CWRTFO_MAT false: do not write to LUKRM3) :
C
CSK       WRITE(LUPRI,'(//A/)') 'CMO input '
CSK       CALL WRTMATMN(CMO,1,NCMOTQ,1,NCMOTQ,LUPRI)
          CALL rGETFC(CWRTFO_MAT,CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
        END IF

        IF (IPROPT.GE.1) THEN
          WRITE(LUPRI,'(A,1F18.10)') '   Core energy : ',ECORE
        END IF
C
C       Get active-active part of FCmo:
C
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(//A/)') 'complete Fock matrix:'
          CALL WRTMATMN(WORK(KFC),1,N2ORBXQ,1,N2ORBXQ,LUPRI)
        END IF
C
        CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
C
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(//A/)') 'Active Fock matrix:'
          CALL WRTMATMN(WORK(KFCAC),1,N2ASHXQ,1,N2ASHXQ,LUPRI)
        END IF
C
C       Transform quaternion FCACmo to molfdir type FCAC.
C
CSK     IPROPT = 30
        CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPROPT)
C
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(A,I6,I6)') 'Molfdir format FCAC, dimension:',
     &                        (2*NASHT)*(2*NASHT),NASHT
          WRITE(LUPRI,'(A)') 'Real part:'
          CALL WRTMATMN(WORK(KFCACM),1,N_FCACM,1,N_FCACM,LUPRI)
          IF (NZ.GT.1) THEN
            WRITE(LUPRI,'(A)') 'Imaginary part:'
            CALL WRTMATMN(WORK(KFCACM+N_FCACM),1,
     &                    N_FCACM,1,N_FCACM,LUPRI)
          END IF
        END IF
CSK     IPROPT = 0
C
C       WORK(KFCACM) now contains all 1-e integrals
C
        CALL MEMREL('after rGETFC',WORK,1,KFCAC,KFREE,LFREE)
C
C       read two-electron integrals in luciarel/transform_r.F 
C       (subroutine itrctl)
C
        LEN_ALL_INT = 0
        LEN_ALL_INT = NASHT*NASHT*NNASHX*NZ*3
        CALL MEMGET('REAL', KH2AC, 0, WORK, KFREE, LFREE)
      else
        CALL QUIT('LUCI_KRCI: Invalid flag I_IITP.')
      end if
#if defined (VAR_MPI2)
C
C     end of slave reading test
C
      IF( ISET_PARCAL .eq. 1 ) THEN
        CALL TEST_PAR_CALC(PARCAL_SAVE,1)
      END IF
#endif
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Call LUCIAREL
      CALL KRMC_LUCI(WORK(KFCACM),WORK(KH2AC),WORK(KFREE),LFREE,ECORE,
     &               CIRUN)
C
      CALL MEMREL('LUCI_KRCI',WORK,1,KFRSAV,KFREE,LFREE)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_CIINII(CREF,CMO,WORK,LWORK)
C***********************************************************************
C
C     Get CI start guess for KR-MCSCF optimization.
C
C     Written by J. Thyssen and T. Fleig                  - Jan 12 2001
C     Last revision :                           T. Fleig  - Nov 30 2001
C                      Read integral files      T. Fleig  - May    2007
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
      CHARACTER*6 CIRUN,RCISTFILE
      LOGICAL PARCAL_SAVE
C
      DIMENSION WORK(*),CMO(*),CREF(*)
C
#include "../luciarel/clunit.inc"
#include "parluci.h"
C
#include "memint.h"
C
      CIRUN = 'CIINII'
      MZ = MIN(NZ,2)
C
      KFRSAV = KFREE
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
#if defined (VAR_MPI2)
C
C     test fock matrix work for slaves
C
      ISET_PARCAL = 0
      IF( ISET_PARCAL .eq. 1 )THEN
        PARCAL_SAVE = .TRUE.
        CALL TEST_PAR_CALC(PARCAL_SAVE,0)
      END IF
#endif
C
C     Flag for type of integral interface:
C     I_IITP = 1: Old 4IND1XXXX0 driven interface
C     I_IITP = 2: New MRCONEE/MDCINT driven interface
C
      N_FCACM = (2*NASHT)*(2*NASHT)
C
      I_IITP = 1
C
      if (I_IITP.eq.2) then
C
        print*,' Testing new integral interface (MRCONEE/MDCINT)'
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDK,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDL,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRKLR,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRKLI,NZ,WORK,KFREE,LFREE)
        print*,'Allocated for ',N_FCACM,' integrals, MZ=',MZ
        CALL LCR_RD_MRCMDC(WORK(KFCACM),WORK(KFCACM+N_FCACM),
     &                     WORK(KINDK),WORK(KINDL),
     &                     WORK(KRKLR),WORK(KRKLI))
        CALL MEMGET('REAL',KH2AC,NASHT*NASHT*NNASHX*NZ*3,
     &               WORK,KFREE,LFREE)
        print*,'Dimension of MOLFDIR array : ',(2*NASHT)**4 * 2
        print*,'Dimension of (NZ,3)  array : ',NASHT*NASHT*NNASHX*NZ*3

        stop 'FIXME: this interface is not complete'
        call memget('INTE',kibeig,norbt,work,kfree,lfree)
        call izero(work(kibeig),norbt)

        call rgeth2(dummy,work(kh2ac),dummy,work(kibeig),.false.,
     &              .true.,.true.,work(kfree),lfree)
C
      else if (I_IITP.eq.1) then
C
C       Calculate core Hamiltonian:
C       ---------------------------
C
C       The core Hamiltonian is the active-active part of FC,
C       a.k.a. FCAC.
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
        IF( CRDFO_MAT ) THEN
C
C         read FCmo from file LUKRM3
C
          CALL KRMC_GETFCK(WORK(KFC),IPROPT)
C
          ECORE = ECORE_LR
C
        ELSE
C
C         Calculate FCmo ( CWRTFO_MAT false: do not write to LUKRM3) :
C
          CALL rGETFC(CWRTFO_MAT,CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
        END IF

        IF (IPROPT.GE.1) THEN
          WRITE(LUPRI,'(A,1F18.10)') '   Core energy : ',ECORE
        END IF
C
C       Get active-active part of FCmo:
C
        CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(A)') 'Active Fock matrix:'
          CALL WRTMATMN(WORK(KFCAC),1,N2ASHXQ,1,N2ASHXQ,LUPRI)
        END IF
C
C       Transform quaternion FCACmo to molfdir type FCAC.
C
        CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPROPT)
C
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(A,I6,I6)') 'Molfdir format FCAC, dimension:',
     &                        (2*NASHT)*(2*NASHT),NASHT
          WRITE(LUPRI,'(A)') 'Real part:'
          CALL WRTMATMN(WORK(KFCACM),1,N_FCACM,1,N_FCACM,LUPRI)
          IF (NZ.GT.1) THEN
            WRITE(LUPRI,'(A)') 'Imaginary part:'
            CALL WRTMATMN(WORK(KFCACM+N_FCACM),1,
     &                    N_FCACM,1,N_FCACM,LUPRI)
          END IF
        END IF
C
C       WORK(KFCACM) now contains all 1-e integrals
C
        CALL MEMREL('after rGETFC',WORK,1,KFCAC,KFREE,LFREE)
C
C       read two-electron integrals in luciarel/transform_r.F 
C       (subroutine itrctl)
C
        LEN_ALL_INT = 0
        LEN_ALL_INT = NASHT*NASHT*NNASHX*NZ*3
        CALL MEMGET('REAL', KH2AC, 0, WORK, KFREE, LFREE)
C
C       FIXME: write me
C              SPINFR=.TRUE. not yet implemented in RGETH2 (why?)
C
      else
        CALL QUIT('LUCI_CIINII: Invalid flag I_IITP.')
      end if
#if defined (VAR_MPI2)
C
C     end of slave reading test
C
      IF( ISET_PARCAL .eq. 1 ) THEN
        CALL TEST_PAR_CALC(PARCAL_SAVE,1)
      END IF
#endif
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Call information transfer routine
C FIXME: Direct handling hardwired for the moment. Later, this
C        should be made dynamical and controlled from the outside.
C        Possibly also integrals reside on file in LUCIAREL
      IO_RCIST_IN = 0
C
C     SK - 03-04-08: analyzing converged MCSCF vector with LUCIAREL
C                    set MOPT_MXMACRO to 0 and MAXCIT to 0
C
      IF( MOPT_MXMACRO .eq. 0 .and. MAXCIT .eq. 0 ) IO_RCIST_IN = 1
      IF (IO_RCIST_IN.EQ.1) THEN
CC      CALL KR_LUCI_IO(0,0,1,DUMMY,1,DUMMY,1,
CC   &                  WORK(KFCACM),N2ASHXQ,WORK(KH2AC),
CC   &                  (2*NASHT)**4*2,
CC   &                  2,IRIQ,
CC   &                  IDUMMY,IDUMMY,IDUMMY,IDUMMY,IPROPT)
        REWIND LUC
        CALL KR_LUCI_IO(1,0,0,CREF,NZCONF*MZ,
     &                  DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  2,MZ,
     &                  LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      END IF
C
C
C
C     Call LUCIAREL
      CALL KRMC_LUCI(WORK(KFCACM),WORK(KH2AC),WORK(KFREE),LFREE,ECORE,
     &               CIRUN)
C
C     Read CI vector from LUCIARELs file
      IF (MOPT_MXMACRO .ge. 0 .and. MAXCIT.GE.0.AND.NCIROOT.EQ.1) THEN
csk+hjaaj, 03JUL07: .MAX MACRO negative is special code for not writing
c   CREF to disk. Dirty fix for doing large scale CI in LUCIAREL, where
C   memory is too small for CREF.
        IO_RCIST_OUT = 1
* FIXME: Future use of file names instead of the primitive fortran
*        files. Activation flags (IO_RCIST_OUT) should be controlled
*        from the outside, not hardwired.
        IF (IO_RCIST_OUT.EQ.1) THEN
          CALL KR_LUCI_IO(1,0,0,CREF,NZCONF*MZ,DUMMY,1,
     &                    DUMMY,1,DUMMY,1,
     &                    1,MZ,
     &                    LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
        END IF
      END IF
C
      CALL MEMREL('LUCI_CIINII',WORK,1,KFRSAV,KFREE,LFREE)
C
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_IJKLRO(CMO,WORK,LWORK)
C***********************************************************************
C
C     Do integral resorting and write to file
C    
C     !!!!!!!!!!!!!!!!!!!!!!  GAS-scheme specific !!!!!!!!!!!!!!!!!!!!!!
C
C     Written by S. Knecht                                - Jan 09 2008
C     Last revision :                                                  
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CHARACTER*6 CIRUN
C
      DIMENSION WORK(*), CMO(*)
C
#include "parluci.h"
C
#include "memint.h"
C
      CIRUN = 'IJKLRO'
      MZ = MIN(NZ,2)
C
      KFRSAV = KFREE
C
      N_FCACM = (2*NASHT)*(2*NASHT)
C
C     Calculate core Hamiltonian:
C     ---------------------------
C
C     The core Hamiltonian is the active-active part of FC,
C     a.k.a. FCAC.
C
      CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
      IF( CRDFO_MAT ) THEN
C
C       read FCmo from file LUKRM3
C
        CALL KRMC_GETFCK(WORK(KFC),IPROPT)
C
        ECORE = ECORE_LR
C
      ELSE
C
C       Calculate FCmo ( CWRTFO_MAT false: do not write to LUKRM3) :
C
        CALL rGETFC(CWRTFO_MAT,CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
      END IF

      IF (IPROPT.GE.1) THEN
        WRITE(LUPRI,'(A,1F18.10)') 'Core energy : ',ECORE
      END IF
C
C     Get active-active part of FCmo:
C
      CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
      IF (IPROPT.GE.5) THEN
        WRITE(LUPRI,'(A)') 'Active Fock matrix:'
        DO I=1,N2ASHXQ,1
          WRITE(6,'(1I3,1F18.14)') I,WORK(KFCAC+I-1)
        END DO
      END IF
C
C     Transform quaternion FCACmo to molfdir type FCAC.
C
      CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPROPT)
C
      IF (IPROPT.GE.5) THEN
        WRITE(LUPRI,'(A)') 'Molfdir format FCAC:'
        WRITE(LUPRI,'(A)') 'Real part:'
        DO I=1,N_FCACM
          WRITE(6,'(1I3,1F18.14)') I,WORK(KFCACM+I-1)
        END DO
        IF (NZ.GT.1) THEN
          WRITE(LUPRI,'(A)') 'Imaginary part:'
          DO I=N_FCACM+1,MZ*N_FCACM
            WRITE(6,'(1I3,1F18.14)') I,WORK(KFCACM+I-1)
          END DO
        END IF
      END IF
C
      CALL MEMREL('after rGETFC',WORK,1,KFCAC,KFREE,LFREE)
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     read two-electron integrals in luciarel/transform_r.F 
C     (subroutine itrctl) 
C
      LEN_ALL_INT = 0
      LEN_ALL_INT = NASHT*NASHT*NNASHX*NZ*3
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Call LUCIAREL
      CALL KRMC_LUCI(WORK(KFCACM),DUMMY,WORK(KFREE),LFREE,DUMMY,
     &               CIRUN)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_SIGMA(SAVDIA,DIAG,CB,SB,FCACM,H2ACM,
     &                      WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate sigma vector: SB = H x CB, where H is the Hamiltonian
C     with integrals FCACM and H2ACM.
C     If (SAVDIA) also calculate H diagonal in DIAG
C
C     Written by J. Thyssen and T. Fleig         - Jan 12 2001
C     Last revision :       T. Fleig             - Mar 16 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
      LOGICAL   SAVDIA
      CHARACTER*6 CIRUN
      DIMENSION DIAG(*),CB(*),SB(*),FCACM(*),H2ACM(*)
      DIMENSION WORK(*)
C
#include "../luciarel/clunit.inc"
C
      KFRSAV = KFREE
C
      MZ = MIN(NZ,2)
C
      CIRUN = 'SIGMA '
      IF (SAVDIA) THEN
        CIRUN = 'SIGMAD'
      ELSE
        CIRUN = 'SIGMA '
      END IF
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Setup LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Write C vector to LUCIAREL restart file:
C           FIXME: Control flag setting (IO_SIGMA_IN) from outside.
C                  (If modification to storage mode is attempted.)
C
      IO_SIGMA_IN = 1
CSK      IPROPT = 10
      IF (IO_SIGMA_IN.eq.1) THEN
        REWIND LUC
        CALL KR_LUCI_IO(1,0,0,CB,NZCONF*MZ,
     &                  DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  2,MZ,
     &                  LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      END IF
CSK      IPROPT = 0
C
C     Call LUCIAREL
      CALL KRMC_LUCI(FCACM,H2ACM,WORK(KFREE),LFREE,DUMMY,CIRUN)
C
C     Read sigma vector on file to SB array.
      IO_SIGMA_OUT = 1
      IF (IO_SIGMA_OUT.eq.1) THEN
        CALL KR_LUCI_IO(1,0,0,SB,NZCONF*MZ,
     &                  DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  1,MZ,
     &                  LUHC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      END IF
C     ... only MASTER should do it
      CALL REWINE(LUC,-1)
      IO_C_VEC_OUT = 1
      IF (IO_C_VEC_OUT .eq. 1 ) THEN
        CALL KR_LUCI_IO(1,0,0,CB,NZCONF*MZ,
     &                  DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  1,MZ,
     &                  LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      END IF
C
      IF (SAVDIA) THEN
C       Read diagonal on file LUDIA to WORK(KDIAG) array.
C
C       only "real part" of diagonal needed!
C       loop over ILOOP_DIA instead of IRIQ! /SK - June 30 2007
        ILOOP_DIA = 1
        CALL KR_LUCI_IO(0,1,0,
     &                  DUMMY,1,DIAG,NZCONF,
     &                  DUMMY,1,DUMMY,1,
     &                  1,ILOOP_DIA,
     &                  IDUMMY,LUDIA,IDUMMY,IDUMMY,IPROPT)
      END IF
CSK      IPROPT = 0
C
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_DENS(CALCR1,CALCR2,RHO1,RHO2,CREFL,CREFR,
     &                     WORK,KFREE,LFREE)
C***********************************************************************
C
C     Read a bra vector CL and a ket vector CR and compute requested
C     density matrices RHO1 (and RHO2)
C
C     CIRUN = 'DENS1 ' : only one-particle density matrix
C     CIRUN = 'DENS2 ' : complete          density matrix
C
C     MOPT_MXMACRO < 0 : calculation of complete density matrix for 
C                        final root IDENSLR_STATE (dcbopt.h)
C
C     CREFL            : Left CI vector
C     CREFR            : Right CI vector
C       general also for transition densities
C
C     Written by T. Fleig                  - Jan 24 2001
C     Last revision : KR_LUCI_IO calls     - Dec 04 2001    T. Fleig
C                     CREF added           - Mar 17 2002    T. Fleig
C                     Only R1 possible     - Feb 01 2005    T. Fleig
C                     Only one common call 
C                     to KRMC_LUCI since
C                     RHO1,RHO2 are always
C                     allocated.           - Jun 29 2007    S. Knecht
C                     implementation of 
C                     multi-root
C                     natural orbital 
C                     occupation numbers   - Sep 28 2007    S. Knecht
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
      LOGICAL CALCR1,CALCR2, natolcr_temp
      CHARACTER*6 CIRUN
C
      DIMENSION WORK(*),RHO1(*),RHO2(*),CREFL(*),CREFR(*)
C
#include "../luciarel/clunit.inc"
#include "parluci.h"
C
      KFRSAV = KFREE
C
      MZ = MIN(NZ,2)
C
      if (CALCR2) then
        CIRUN = 'DENS2 '
        if (.not.CALCR1)
     &       call quit('Not implemented option in LUCI_DENS.')
      else if (CALCR1) then
        CIRUN = 'DENS1 '
      else
        call quit('Check CALCR setting in LUCI_DENS')
      end if

!     check for natural orbital occupation number flag
      natolcr_temp = natolcr 
      if(calcr2) natolcr = .false.
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Initialize density matrices to zero
      L1 = (2*NASHT)**2 * MZ
      CALL DZERO(RHO1,L1)
      if (CALCR2) then
        L2 = (2*NASHT)**4 * MZ
        CALL DZERO(RHO2,L2)
      end if
C
      LF2_ZERO = L2
C
#ifdef MCSCF_DEBUG
C     Write current CI vector(s) to file - debug version
      idothishere = 1
      REWIND LUC
      CALL KR_LUCI_IO(1,0,0,CREFL,NZCONF*MZ,
     &                DUMMY,1,
     &                DUMMY,1,DUMMY,1,
     &                2,MZ,
     &                LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      CALL KR_LUCI_IO(1,1,0,CREFL,NZCONF*MZ,
     &                CREFR,NZCONF*MZ,
     &                DUMMY,1,DUMMY,1,
     &                2,MZ,
     &                LUSC61,LUSC62,IDUMMY,IDUMMY,IPROPT)
#else
!     print *, 'decide branch: natolcr    ==> ',natolcr
!     print *, 'decide branch: cana_mcscf ==> ',cana_mcscf
!     print *, 'decide branch: calcr1     ==> ',calcr1
!     print *, 'decide branch: calcr2     ==> ',calcr2
      if((calcr1.and..not.natolcr.and..not.cana_mcscf) .or. calcr2)then ! 

        ! mcscf 1-/2-particle density matrix calculation

!       print *, 'decide branch: calcr2     ==> ',calcr2
!       xxx = 0.0d0
!       print *, 'NZCONF is ==> ',NZCONF
!       xxx = ddot(NZCONF*MZ,CREFL,1,CREFL,1)
!       print *, 'xxx lhs is ==> ',xxx
!       xxx = ddot(NZCONF*MZ,CREFR,1,CREFR,1)
!       print *, 'xxx hhs is ==> ',xxx
      
!       Write current CI vector(s) to file - choice 1: MCSCF run setup:
!       lhs --> LUSC61
!       rhs --> LUSC62
        CALL KR_LUCI_IO(1,1,0,CREFL,NZCONF*MZ,
     &                  CREFR,NZCONF*MZ,
     &                  DUMMY,1,DUMMY,1,
     &                  2,MZ,
     &                  LUSC61,LUSC62,IDUMMY,IDUMMY,IPROPT)
      end if ! mcscf run setup

      if(cana_mcscf)then ! ci analysis of mcscf vector (restart from the file KRMCOLD)
!       Write current CI vector(s) to file - choice 1: MCSCF run setup:
!       vector --> LUC 
!       internally we copy the 
!       lhs --> LUSC61
!       rhs --> LUSC62
        is_ana_mcscf = .true.
        print *, 'read mcscf vector from file KRMCOLD (unit 59)...'
        CALL MEMGET('REAL',KRSCR,NZCONF*MZ,WORK,KFREE,LFREE)
        CALL REAKRMC(59,'CREF    ',WORK(KRSCR),NZCONF*MZ)
        REWIND LUC
!       xxx = 0.0d0
!       xxx = ddot(NZCONF*MZ,WORK(KRSCR),1,WORK(KRSCR),1)
!       print *, 'xxx is ==> ',xxx
        CALL KR_LUCI_IO(1,0,0,WORK(KRSCR),NZCONF*MZ,
     &                  DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  2,MZ,
     &                  LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      
        CALL MEMREL('LUCI_DENS',WORK,1,KFRSAV,KFREE,LFREE)
      END IF
#endif
C
C     RHO1 and RHO2 are always allocated in the calling routines.
C     RHO2 is essential for communication in parallel runs.
C     / SK June 29 2007.
C
      CALL KRMC_LUCI(RHO1,RHO2,WORK(KFREE),LFREE,DUMMY,CIRUN)

!
!     possibly re-set the natural orbital occupation number flag to its entry value
      natolcr = natolcr_temp
C
C     close files ...
      CALL LUCI_CLOSEDOWN()
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_ANAL(CREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Analyze CI vectors
C
C     Written by T. Fleig                  - Jan 24, 2001
C     Last revision : T. Fleig
C                     activated            - Jan 05, 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
#include "../luciarel/clunit.inc"
C
      CHARACTER*6 CIRUN
C
      DIMENSION WORK(*)
C
      KFRSAV = KFREE
C
      MZ = MIN(NZ,2)
C
      CIRUN = 'ANALYZ'
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C      WRITE(6,*) 'NCIROOT and MOPT_MXMACRO',NCIROOT,MOPT_MXMACRO
C
      IF ( NCIROOT.EQ.1 .AND. MOPT_MXMACRO .gt. 0) THEN
        IO_ANAL_IN = 1
        IF (IO_ANAL_IN.EQ.1) THEN
          REWIND LUC
          CALL KR_LUCI_IO(1,0,0,CREF,NZCONF*MZ,
     &                          DUMMY,1,
     &                    DUMMY,1,DUMMY,1,
     &                    2,MZ,
     &                    LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
        END IF
      ELSE IF (NCIROOT.GT.1) THEN
        CONTINUE
      END IF
C
C     Call LUCIAREL for analysis of CI vectors.
      CALL KRMC_LUCI(DUMMY,DUMMY,WORK(KFREE),LFREE,DUMMY,CIRUN)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      END

      SUBROUTINE LUCI_QCORR(CREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate Davidson-type energy correction
C
C     Written by S. Knecht                 - Jan 12, 2014
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
#include "../luciarel/clunit.inc"
C
      CHARACTER*6 CIRUN
C
      DIMENSION WORK(*)
C
      KFRSAV = KFREE
C
      MZ = MIN(NZ,2)
C
      CIRUN = 'QCORR '
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
      IF ( NCIROOT.EQ.1 .AND. MOPT_MXMACRO .gt. 0) THEN
        IO_ANAL_IN = 1
        IF (IO_ANAL_IN.EQ.1) THEN
          REWIND LUC
          CALL KR_LUCI_IO(1,0,0,CREF,NZCONF*MZ,
     &                          DUMMY,1,
     &                    DUMMY,1,DUMMY,1,
     &                    2,MZ,
     &                    LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
        END IF
      ELSE IF (NCIROOT.GT.1) THEN
        CONTINUE
      END IF
C
C     Call LUCIAREL for analysis of CI vectors.
      CALL KRMC_LUCI(DUMMY,DUMMY,WORK(KFREE),LFREE,DUMMY,CIRUN)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      END

      SUBROUTINE LUCI_refvec(CREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     save reference vector in format for QCORR module
C
C     Written by S. Knecht                 - Jan 13, 2014
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
#include "../luciarel/clunit.inc"
C
      CHARACTER*6 CIRUN
C
      DIMENSION WORK(*)
C
      KFRSAV = KFREE
C
      MZ = MIN(NZ,2)
C
      CIRUN = 'REFVEC'
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
      IF(NCIROOT > 1) 
     &call quit('refvec writing only for MCSCF wave function')
      REWIND LUC
      CALL KR_LUCI_IO(1,0,0,CREF,NZCONF*MZ,
     &                DUMMY,1,
     &                DUMMY,1,DUMMY,1,
     &                2,MZ,
     &                LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
C
C     Call LUCIAREL for writing the reference vector in QCORR format
      CALL KRMC_LUCI(DUMMY,DUMMY,WORK(KFREE),LFREE,DUMMY,CIRUN)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_RSTRMC(WORK,LWORK)      
C***********************************************************************
C
C     Restart information for LUCIAREL in MCSCF runs.
C
C     Written by T. Fleig                 - Aug 28, 2003
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
C
#include "dcbopt.h"
C
      CHARACTER*6 CIRUN
      DIMENSION WORK(*)
C
      CIRUN = 'RSTRMC'
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Call LUCIAREL for setting up vector partitioning information.
      CALL KRMC_LUCI(DUMMY,DUMMY,WORK,LWORK,DUMMY,CIRUN)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRCIGETPRP(VPROP_KRCI,WORK,LWORK)
C***********************************************************************
C
C     CIRUN = 'PROP1 ' : compute one-electron properties using 
C                        property integrals residing on file KRMC_FOCK
C
C     Written by S. Knecht - Sep 2008
C
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dummy.h"
#include "dcbgen.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krciprop.h"
#include "krmcluci_inf.h"
      CHARACTER*6 CIRUN
      DIMENSION VPROP_KRCI(*), WORK(*)
C               VPROP_KRCI(NPROP_ROOTS_KRCI,NPROP_ROOTS_KRCI,NPROP_KRCI,MZ)
#include "parluci.h"
      integer :: print_level
C
#include "memint.h"
C
      CALL QENTER('KRCIGETPRP')
C
      print_level = 0
      MZ          = MIN(NZ,2)
C     output file handle
      LUWRT = LUPRI
C
C     Set CIRUN
      CIRUN = 'PROP1 '
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Determine operator symmetries in sub-double group used in LUCIAREL
      CALL GET_PRPOPSYM_KRCI(WORK(KFREE),LFREE,print_level)
C
C     Get reordered operator matrices - write back to KRMC_FOCK on label 
C     LABEL = PRPNAM(property operator)(1:6)RO (first 6 chars of
C     property name + RO)
      CALL GET_XPROP_REOD_KRCI(LUKRM3,WORK,LWORK,CIRUN)
C
C     calculate property matrix VPROP_KRCI - from now on we also run in 
C     parallel if it is asked for.
      CALL KRMC_LUCI(DUMMY,VPROP_KRCI,WORK,LWORK,DUMMY,CIRUN)
C
C     CALL MEMREL('KRCIGETPRP',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     close files ...
      CALL LUCI_CLOSEDOWN()
C
      CALL QEXIT('KRCIGETPRP')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KR_LUCI_IO(IVA,IVB,II,VA,IDIMA,VB,IDIMB,
     &                      X1,IDIM1,X2,IDIM2,
     &                      IOTYPE,IRILP,
     &                      LU1,LU2,LU3,LU4,IPRNT)
C***********************************************************************
C
C     Transfer vectors on file between LUCIAREL and KR-MCSCF
C
C     Written by T. Fleig                  - Nov 28 2001
C     Last revision : T. Fleig             - Mar    2002
C
C
C  Different cases:
C
C    1) Vector 1           IVA
C    2) Vector 2           IVB
C    3) Integral quantity  II
C     Quantity specifier:  0 no quantity
C                          1 quantity
C
C   I/O ?
C
C   Read file to array     IOTYPE = 1
C          for output   (KRMC <- LUCIAREL).
C
C   Write array to file    IOTYPE = 2
C          for input    (KRMC -> LUCIAREL)
C
C***********************************************************************
#include "implicit.h"
*
      dimension VA(IDIMA),VB(IDIMB),X1(IDIM1),X2(IDIM2)
*
      NTESTL = 00
CSK      IPRNT = 000
      NTEST = max(NTESTL,IPRNT/1000)
*
      if (NTEST.ge.5) then
*   Print section:
*===========================================
        write(6,*)
        write(6,*) '#########################'
        write(6,*) ' KR_LUCI_IO activity     '
        write(6,*) '#########################'
        if (IVA.eq.1.or.IVB.eq.1) then
          write(6,*) '  Vector task   '
        else if (II.eq.1) then
          write(6,*) '  Integral / density task   '
        end if
*
* A) Output from LUCIAREL (read from file)
        if (IOTYPE.eq.1) then
          if (NTEST.ge.1) write(6,*) '  Output from LUCIAREL '
*
* B) Input for LUCIAREL (write to file)
        else if (IOTYPE.eq.2) then
          if (NTEST.ge.1) write(6,*) '  Input for LUCIAREL '
        end if
*===========================================
      end if
*
*===============
*  Load vectors
*===============
      if (IOTYPE.eq.1) then
*  Load vector A (ket)
        if (IVA.eq.1) then
          if (NTEST.ge.10) write(6,*) '   *** Loading vector A ***'
          call rdvcdc(VA,LU1,-1,IRILP,IPRNT)
        end if
*
*  Load vector B (bra)
        if (IVB.eq.1) then
          if (NTEST.ge.10) write(6,*) '   *** Loading vector B ***'
          call rdvcdc(VB,LU2,-1,IRILP,IPRNT)
        end if
*
*  Load 1-particle quantity
*  Load 2-particle quantity
        if (II.eq.1) then
          call rdvcdc(X1,LU3,-1,IRILP,IPRNT)
          call rdvcdc(X2,LU4,-1,IRILP,IPRNT)
        end if
*
*===============
*  Save vectors
*===============
*  (Write as REAL PART, -1 separator, IMAG PART, -1 trailer)
      else if (IOTYPE.eq.2) then
*  Save vector A (ket)
        if (IVA.eq.1) then
          if (NTEST.ge.10) write(6,*) '   *** Saving vector A ***'
          call savevcdc(VA,LU1,IRILP,IPRNT)
        end if
*
*  Save vector B (bra)
        if (IVB.eq.1) then
          if (NTEST.ge.10) write(6,*) '   *** Saving vector B ***'
          call savevcdc(VB,LU2,IRILP,IPRNT)
        end if
*
*  Save integral quantity
        if (II.eq.1) then
          write(6,*) 'Saving integrals to file not implemented yet.'
          stop 'Quitting.'
        end if
*
      end if
*     ^ type of IO
*
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUMSTI_REL
C
C     Call the slaves if we run in parallel for KRMC-LUCIAREL
C
C     adaption of the corresp. RELCCSD routine written by Luuk Visscher
C     KRMC-LUCIAREL adaption: S. Knecht
C
C     Last revision: S. Knecht           - March  2007
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "dcbgen.h"
#include "parluci.h"
C
C     Send task code for parallel KRMC-LUCIAREL
C
      IF (PARCAL) CALL DIRAC_PARCTL( LUCIAREL_PAR )
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUMSTE_REL
C
C     Release the slaves if we run in parallel for KRMC-LUCIAREL
C
C     adaption of the corresp. RELCCSD routine written by Luuk Visscher
C     KRMC-LUCIAREL adaption: S. Knecht
C
C     Last revision: S. Knecht           - March  2007
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "dcbgen.h"
#include "infpar.h"
#include "parluci.h"
C
      IF (.NOT.PARCAL) RETURN
C
C     Let the slaves sleep (if relevant for this operating system)
C
      CALL DIRAC_PARCTL( RELEASE_NODES )
C
      END
#if defined (VAR_MPI2)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCIAREL_NODE()
C
C     KRMC-LUCIAREL driver routine for the slaves
C
C     adapted version of the slave RELCCSD routine 
C     written by Luuk Visscher.
C     KRMC-LUCIAREL adaption: S. Knecht
C
C     Last revision: S. Knecht           - Nov  2008
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use memory_allocator
C
#include "implicit.h"
#include "priunit.h"
#include "infpar.h"
#include "parluci.h"
      LOGICAL LUEND
      CHARACTER LUCITABASF_REL*12, LUCIFILN_REL*16
      CHARACTER CIRUN*6
C
#include "../luciarel/mxpdim.inc"
#include "../luciarel/cicisp.inc"
#include "../luciarel/noccn_inf.inc"
#include "krciprop.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"

      real(8), allocatable :: WORK(:)
      call legacy_lwork_get(LWORK)
      call alloc(WORK,LWORK,id='WORK in LUCIAREL_NODE')

#include "memint.h"
C
      KFRSAV = KFREE
C
C       START 
C     =========
C
C     arrange for the MPI stuff and correct node number
C     to the total number of running invocations.
      MASTER = MPARID
      MYPROC = MYTID
C     Add the master node, NUMNOD = number of slaves
      NMPROC = NUMNOD + 1
C
C     initialize core energy
      ECORE_KRMC = 0.0D0
C
C     create a node-unique filename as output file. Important on
C     shared file systems. Otherwise all the output gets mingled in one
C     file. You don't really want to do this.
      LUCITABASF_REL="LUCIAREL.OUT"
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (LUCIFILN_REL,'(A12,A1,I1)') LUCITABASF_REL,'.',MYPROC
         LUFIL=14
      ELSE IF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (LUCIFILN_REL,'(A12,A1,I2)') LUCITABASF_REL,'.',MYPROC
         LUFIL=15
      ELSE IF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (LUCIFILN_REL,'(A12,A1,I3)') LUCITABASF_REL,'.',MYPROC
         LUFIL=16
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND LUCIAREL_NODE MODULE")
      ENDIF
C
C     open the local input file and the node specific output file.
C     Every access to the local stdout handle then automatically writes
C     to the corresponding output file.
      OPEN(MSLVOUT_REL,FILE = LUCIFILN_REL(1:LUFIL))
C
C     transfer file handle to common block
      LUWRT      = MSLVOUT_REL
      lupri_save = lupri
      lupri      = LUWRT
C
C     synchronize nodes, part 1
      CALL SYNC_NODES_REL_P1(ECORE_KRMC,CIRUN,MSLVOUT_REL)
C
C     set property module run flag
      RUNXPROP = .FALSE.
C     MZ value
      MZ = MIN(NZ,2)
C
      IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'SIGMA ' .or. 
     &    CIRUN .eq. 'SIGMAD' .or. CIRUN .eq. 'KR-CI ' ) THEN
C
C       read-in of integrals by nodes... no broadcast in SYNC_NO...P2
        IF( (CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' ) 
     &      .and. IIOMOD_REL .eq. 1 )THEN
          PARCAL = .FALSE.
        END IF
C
C     Flag for type of integral interface:
C     I_IITP = 1: Old 4IND1XXXX0 driven interface
C     I_IITP = 2: New MRCONEE/MDCINT driven interface
C
      N_FCACM = (2*NASHT)*(2*NASHT)
      I_IITP = 1
C
      if (I_IITP.eq.2) then
C
        print*,' Testing new integral interface (MRCONEE/MDCINT)'
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDK,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDL,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRKLR,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRKLI,NZ,WORK,KFREE,LFREE)
        print*,'Allocated for ',N_FCACM,' integrals, MZ=',MZ
        CALL LCR_RD_MRCMDC(WORK(KFCACM),WORK(KFCACM+N_FCACM),
     &                     WORK(KINDK),WORK(KINDL),
     &                     WORK(KRKLR),WORK(KRKLI))
        CALL MEMGET('REAL',KH2AC,NASHT*NASHT*NNASHX*NZ*3,
     &               WORK,KFREE,LFREE)
        print*,'Dimension of MOLFDIR array : ',(2*NASHT)**4 * 2
        print*,'Dimension of (NZ,3)  array : ',NASHT*NASHT*NNASHX*NZ*3

!       get dimensions for scratch matrices
!       SK - Aug 2010: FIXME: this structure as it is cannot work 
!       if the integral files are not available to the co-workers... 
!       @ Timo: before continuing (if ever) we have to think about how
!       to handle this issue properly.

        stop 'FIXME: this interface is not complete'
        call memget('INTE',kibeig,norbt,work,kfree,lfree)
        call izero(work(kibeig),norbt)

        call rgeth2(dummy,work(kh2ac),dummy,work(kibeig),.false.,
     &              .true.,.true.,work(kfree),lfree)
C
      else if (I_IITP.eq.1) then
C
C       Calculate core Hamiltonian:
C       ---------------------------
C
C       The core Hamiltonian is the active-active part of FC,
C       a.k.a. FCAC.
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KFCACM),MZ*N_FCACM)
C
C       Read two-electron integrals:
C       ----------------------------
C
C       ... now inside LUCIAREL (subroutine itrctl)
        LEN_ALL_INT = NASHT*NASHT*NNASHX*NZ*3
C
        IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' )THEN
          CALL MEMGET('REAL',KH2AC,          0,WORK,KFREE,LFREE)
        ELSE
          CALL MEMGET('REAL',KH2AC,LEN_ALL_INT,WORK,KFREE,LFREE)
        END IF
C
      else
        CALL QUIT('LUCI_RCIST: Invalid flag I_IITP.')
      end if
C
        PARCAL = .TRUE.
C
C       synchronize nodes, part 2 (now CIRUN dependent!)
        CALL SYNC_NODES_REL_P2(WORK(KFCACM),WORK(KH2AC),CIRUN)
      ELSE IF( CIRUN .eq. 'DENS1' .or. CIRUN .eq. 'DENS2' ) THEN 
C
C       Initialize density matrices to zero; currently I want to test
C       this CIRUN, therefore we make no difference between DENS1 and
C       DENS2 - SK June 2007
C
        L1 = (2*NASHT)**2 * MZ
        CALL MEMGET('REAL',KRHO1,L1,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KRHO1),L1)
        L2       = (2*NASHT)**4 * MZ
        LF2_ZERO = (2*NASHT)**4 * MZ
C
      ELSE IF( CIRUN .eq. 'IJKLRO' ) THEN 
C
        N_FCACM     = (2*NASHT)*(2*NASHT)
        LEN_ALL_INT = NASHT*NASHT*NNASHX*NZ*3
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KH2AC,          0,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KFCACM),MZ*N_FCACM)
C
        CALL SYNC_NODES_REL_P2(WORK(KFCACM),WORK(KH2AC),CIRUN)
      ELSE IF( CIRUN .eq. 'PROP1 ') THEN 
        RUNXPROP = .TRUE.
        CALL SYNC_NODES_REL_XPROP
C
C       allocate for property matrix
        LXPRPKRCI = MZ * NPROP_KRCI * NPROP_ROOTS_KRCI**2
        CALL MEMGET('REAL',KXPRPKRCI,LXPRPKRCI,WORK,KFREE,LFREE)
C
C       ... initialize
        CALL DZERO(WORK(KXPRPKRCI),LXPRPKRCI)
      END IF
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,1)
C
C     attention to nat. orb. occ. num. calculation
      IF( CIRUN .eq. 'DENS1' .or. CIRUN .eq. 'DENS2' ) THEN
       IF( .NOT. NOOCCN_DENS_RUN ) THEN
          CALL MEMGET('REAL',KRHO2,L2,WORK,KFREE,LFREE)
          CALL DZERO(WORK(KRHO2),L2)
       ELSE
          CALL MEMGET('REAL',KRHO2, 0,WORK,KFREE,LFREE)
       END IF
      END IF
C
C     Setup routine for running LUCIAREL
      CALL LUCI_SETUP()
C
C     initialize common block
      CALL MPIXCALC_TYPESZ()
C
C     Run LUCIAREL
C
      IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'SIGMA ' .or. 
     &    CIRUN .eq. 'SIGMAD' .or. CIRUN .eq. 'KR-CI ' ) THEN
C
C       Call LUCIAREL for calculating sigma vector(s).
        CALL KRMC_LUCI2(WORK(KFCACM),WORK(KH2AC),WORK(KFREE),LFREE,
     &                  ECORE_KRMC,CIRUN)
        CALL MEMREL('LUCIAREL_NODE',WORK,1,KFRSAV,KFREE,LFREE)
C
      ELSE IF( CIRUN .eq. 'DENS1' .or. CIRUN .eq. 'DENS2' ) THEN 
C
C       Call LUCIAREL for calculating density matrix elements.
        CALL KRMC_LUCI2(WORK(KRHO1),WORK(KRHO2),WORK(KFREE),
     &                  LFREE,DUMMY,CIRUN)
        CALL MEMREL('LUCIAREL_NODE',WORK,1,KFRSAV,KFREE,LFREE)
C
      ELSE IF(CIRUN == 'ANALYZ' .or. CIRUN == 'QCORR ' .or. 
     &        CIRUN == 'REFVEC') THEN 
C
C       Call LUCIAREL for analyzing the wavefunction.
        CALL KRMC_LUCI2(DUMMY,DUMMY,WORK,LWORK,DUMMY,CIRUN)
C
      ELSE IF( CIRUN .eq. 'NDET  ' ) THEN 
C
C       Call LUCIAREL for determining number of dets.
        CALL KRMC_LUCI2(DUMMY,DUMMY,WORK,LWORK,DUMMY,CIRUN)
C
        NDET = INT(XISPSM(ISYM_KRMC,1))
C
      ELSE IF( CIRUN .eq. 'RSTRMC' ) THEN
C
C       Call LUCIAREL for setting up vector partitioning information.
        CALL KRMC_LUCI2(DUMMY,DUMMY,WORK,LWORK,DUMMY,CIRUN)
C
      ELSE IF( CIRUN .eq. 'IJKLRO' ) THEN 
C
C       Call LUCIAREL for integral reordering.
        CALL KRMC_LUCI2(WORK(KFCACM),WORK(KH2AC),WORK(KFREE),LFREE,
     &                  DUMMY,CIRUN)
        CALL MEMREL('LUCIAREL_NODE',WORK,1,KFRSAV,KFREE,LFREE)
      ELSE IF( CIRUN .eq. 'PROP1 ' ) THEN 
C
C       Call LUCIAREL for property calculation.
        CALL KRMC_LUCI2(DUMMY,WORK(KXPRPKRCI),WORK(KFREE),LFREE,
     &                  DUMMY,CIRUN)
        CALL MEMREL('LUCIAREL_NODE',WORK,1,KFRSAV,KFREE,LFREE)
C
      END IF
C
C     close files ...
      CALL LUCI_CLOSEDOWN()
C
C     LUCIAREL run is finished for the slaves ...
      CLOSE (MSLVOUT_REL,STATUS='KEEP')

      lupri = lupri_save
C
      call dealloc(WORK)
      END
#endif    /* ifdef VAR_MPI2 */
C --- end of krmcluci.F ---
