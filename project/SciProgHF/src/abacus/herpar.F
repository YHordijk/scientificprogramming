#define DEBUG_IPRPAR -1
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
#ifdef UNDEF
/* Comdeck notes */
961125 - tec RELCAL: FIRST replaced by COMMON parfir remember to set FIRST1
961003 - tsaue: DOINT is removed and replaced by I2TYP.
                Some polishing av code:
                  Note that NPOS enters PARDRV as argument
                  and that Fock matrices are not initialized at master !
950620-kr: Ported to MPI, based on previous efforts by Bjarne Herland at
           Para//ab
  --- in PARLOP is checked for DOINT(*,*,1) on the IJSHEL;
      this will call NULL proccesses on nodes if some DOINT values are
      false. It would be better to initialize with DOINT in INDEKS,
      but that requires a little more coding.
      SORT renamed to PARSRT (it is likely that library routines exist
      called SORT).
950428-kr+hjaaj;
 Merged with Odense version. Implemented DOINT in PARLOP.
#endif
C  /* Deck pardrv */
      SUBROUTINE HER_PARDRV(WORK,LWORK,FMAT,DMAT,NDMAT,IREPDM,IFCTYP,
     &                  ITYPE,MAXDIF,IATOM,NODV,NOPV,NOCONT,
     &                  TKTIME,RETUR,FIRST,NPOS,NTOTTK,I2TYP,
     &                  ICEDIF,SCREEN,GABRAO,DMRAO,DMRSO,DINTSKP)
C
C     Written by Henrik Koch and Paal Dahle Nov-1993 -> Okt-1994
C     Extended with I2TYP Oct 1 1996 by T.Saue
C     Extended with screening Oct 11 1996 by T.Saue
C
C
C     **************************************************************
C     *        Driver routine for the parallel calculation.        *
C     **************************************************************
C
C
C      A short survey of the different messagetags (MSGTAG) :
C
C      10 - Send initialization to nodes.
C      20 - Receive request from node.
C      30 - Send batch to node.
C      40 - Receive final results from nodes.
C      50 - Get timing results from nodes.
C
C           ITYPE  - Calc. type: 2 = expectation values.
C                                3 = direct SCF (AO-basis)
C                                4 = (rs|**)-distributions
C                               -5 = Integrals derivated with respect
C                                    to magnetic field.
C                                7 = direct CC (DALTON)
C                                8 = Fock matrices derivated with respect
C                                    to atomic coordinates.
C                                9 = direct SCF (SO-basis)
C
C           NTASK  - Number of tasks
C           IPRPAR - Print level during parallelization.
C           MAXDIF -
C           MAXREP -
C           IREPDM - Symmetries of density matrices
C           IFCTYP - Fock matrix types (see twoint).
C           IATOM  -
C           NODV   - Density matrix neglected in TWOEXP
C           NOPV   -           -   "   -
C           NOCONT -           -   "   -
C           RETUR  - Program will exit after spec. shells
C           TKTIME - Take time in TWOINT
C           TIMING - Take time for each integralbatch IJ
C           DOREPS -
C           DOSYM  -
C           DOCOOR -
C           NMLINE - Number of lines in MOLECULE.INP
C           MLINE  - MOLECULE.INP as internal file
C           NDMAT  - Number of density-matrices
C           DMAT   - Density-matrices.
C
C
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
C
C dirac_partask.h : symbols for DIRAC_PARCTL codes
C
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "molinp.h"
#include "twocom.h"
#if defined (PRG_DIRAC)
#include "dirac_partask.h"
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "inforb.h"
#include "infpar.h"
C
!     DIMENSION FMAT(*), DMAT(N2BASX,NDMAT), WORK(LWORK),
      DIMENSION FMAT(*), DMAT(*), WORK(LWORK),
     &          NPOS(NTOTTK), IREPDM(*), IFCTYP(*),
     &          GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*)
      LOGICAL FIRST,NODV,NOPV,NOCONT,RETUR,TKTIME
      real(8)  :: TIMSTR(2), TIMEND(2)
C
      CALL QENTER('HER_PARDRV')
C
#if !defined (VAR_MPI)
      CALL QUIT
     &('program error: parallel code called but MPI not activated')
#endif
C
      IPRPAR = max(DEBUG_IPRPAR,IPRPAR)
      IF (IPRPAR .GT. 0) THEN
         WRITE(LUPRI,'(/5X,A,I5/)')
     &    'HER_PARDRV: Number of MPI worker nodes:', NODES
      ENDIF
C
C     Take calculation time of each IJ pair?
C
      TIMING = .TRUE.
! initialization
      KFREE=1
C
C     Determine total number of tasks - MTOTTK.
      CALL HERPAR_DET_TASK (I2TYP,1,MTOTTK,NPOS,WORK(KFREE))
      IF(MTOTTK.GT.NTOTTK) THEN
        WRITE(LUPRI,'(A,I8)') 'HER_PARDRV: Inconsistent MTOTTK:',MTOTTK
        WRITE(LUPRI,'(3X,A,I8,A,I8)') ' - should be less or equal to',
     &   NTOTTK,' for I2TYP =',I2TYP
        CALL QUIT('HER_PARDRV: Inconsistent MTOTTK !!!')
      ENDIF
C
C     Allocate local variables
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
      CALL MEMGET2('INTE','NSTAT',KNSTAT,NODES ,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','INDEX',KINDEX,MTOTTK,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','TIMES',KTIMES,MTOTTK,WORK,KFREE,LFREE)
C
C     Test of work space allocation.
C
      LDINTSKP = 8
      IF (ITYPE .EQ. 2) THEN
         NUMGRA = MXCOOR
         LDINTSKP = 3*3*14
         IF (MAXDIF .GT. 1) THEN
            NUMHES = NUMGRA*NUMGRA
            NREQUI = NUMGRA + NUMHES
         ELSE
            NUMHES = 0
            NREQUI = NUMGRA
         END IF
      ELSE IF (ITYPE .EQ. 3) THEN
         NREQUI = NDMAT*N2BASX
      ELSE IF (ITYPE .EQ. 9) THEN
         NREQUI = NDMAT*N2BASX
      ELSE IF (ITYPE .EQ. -5) THEN
         IF (MAXDIF.EQ.1) THEN
           NREQUI = 3*NDMAT*N2BASX
         ELSEIF(MAXDIF.EQ.2) THEN
           NREQUI = 0
         ELSE
           CALL QUIT
     &     ('HER_PARDRV: higher magnetic derivatives not programmed !')
         ENDIF
      ELSE IF (ITYPE .EQ. 8) THEN
         NREQUI = 3*NNBASX
      END IF

      LDMAT = NDMAT*N2BASX
      IF (NREQUI .LE. LDMAT .AND. LFREE .LT. NREQUI) THEN
         ! use DMAT as scratch space in RVRES
         LREQUI = -1
         open(unit=78, file='DMAT_AO_SAVE', form='UNFORMATTED')
         write(78) DMAT(1:LDMAT)
         rewind (78)
      ELSE
         LREQUI = NREQUI
      END IF
C
      IF (LFREE .LT. LREQUI) THEN
         WRITE(LUPRI,'(/5X,A,2(/5X,A,I12),/5X,A,2I8)')
     &    'Insufficient work space for subroutine RVRES',
     &    'Free space     :',LFREE,
     &    'Required space :',LREQUI,
     &    'ITYPE, NDMAT   :',ITYPE,NDMAT
         CALL QUIT('Insufficient work space for RVRES')
      ENDIF
C
C--------------------------
C     Make the index array.
C--------------------------
C
      CALL HERPAR_INDEKS(FIRST,WORK(KINDEX),WORK(KTIMES),
     &   ITYPE,MTOTTK,NBATCH,NPOS,I2TYP)
C
C---------------------------------------
C     Send initialization data to nodes.
C---------------------------------------
C
      CALL TIMER2('START',TIMSTR,TIMEND)
      CALL SDINIT(DMAT,NDMAT,IREPDM,IFCTYP,ITYPE,MAXDIF,IATOM,NODV,
     &            NOPV,NOCONT,TKTIME,RETUR,I2TYP,ICEDIF,SCREEN,
     &            GABRAO,DMRAO,DMRSO)
C     and transfer screening information to twocom.h
C     (such that master will print correct scr.thr. in summary) /hjaaj apr 01
      SCRTHR = SCREEN
      IF (IPRPAR.GT.0) CALL TIMER2('SDINIT',TIMSTR,TIMEND)
C
C-----------------------
C     Start calculation.
C-----------------------
C
      CALL DO_HERPAR(WORK(KINDEX),MTOTTK,NBATCH,WORK(KNSTAT))
      IF (IPRPAR.GT.0) CALL TIMER2('DO_HERPAR',TIMSTR,TIMEND)
C.....Get results from nodes
      IF (LREQUI .EQ. -1) THEN
         ! use DMAT as scratch space in RVRES
         CALL RVRES(DMAT,FMAT,ITYPE,MTOTTK,NREQUI,
     &           NUMGRA,NUMHES,MAXDIF,WORK(KTIMES),DINTSKP,LDINTSKP)
         rewind (78)
         read(78) DMAT(1:LDMAT)
         CLOSE(78, status='DELETE')
      ELSE
         CALL RVRES(WORK(KFREE),FMAT,ITYPE,MTOTTK,NREQUI,
     &           NUMGRA,NUMHES,MAXDIF,WORK(KTIMES),DINTSKP,LDINTSKP)
      END IF
      IF (IPRPAR.GT.0) CALL TIMER2('RVRES',TIMSTR,TIMEND)
C
C-------------------------------------
C     Get timing results if requested.
C-------------------------------------
C
      IF (IPRPAR .GT. 0) THEN
         CALL TIMRES(WORK(KNSTAT))
         CALL TIMER2('TIMRES',TIMSTR,TIMEND)
      END IF
C
C--------------------------------------
C     Release nodes from this task
C     and let them return to main menu.
C--------------------------------------
C
      CALL DIRAC_PARCTL( EXIT_NODEMENU )
C
C-----------------------------------------
C     Sort the integral-calculation times.
C-----------------------------------------
C
      IF(TIMING) CALL PARSR2(MTOTTK,NPOS,WORK(KTIMES),
     &                       WORK(KINDEX),IPRPAR)
C
      IF (FIRST)  FIRST  = .FALSE.
      IF (NEWBAS) NEWBAS = .FALSE.
C

 90   CONTINUE
      CALL MEMREL('HER_PARDRV',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('HER_PARDRV')
      RETURN
      END
C  /* Deck herpar_indeks */
      SUBROUTINE HERPAR_INDEKS(FIRST,INDEX_IJ,TIMES,
     &   ITYPE,MTOTTK,NBATCH,NPOS,I2TYP)
C
C     Written by Paal Dahle May 1994.
C     Revised Oct 1996 by T.Saue
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "infpar.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "blocks.h"
#include "dcbfir.h"
C
      INTEGER INDEX_IJ(MTOTTK), NPOS(MTOTTK)
      REAL*8  TIMES(MTOTTK)
      LOGICAL FIRST,DOFIRST
C
C     If new run, make initial guess as to run times
C
      DOFIRST = .FALSE.
C
      IF (RELCAL) THEN
         IF (I2TYP.EQ.1.AND.FIRST1) THEN
            DOFIRST = .TRUE.
            FIRST1 = .FALSE.
         ELSEIF (I2TYP.EQ.2.AND.FIRST2) THEN
            DOFIRST = .TRUE.
            FIRST2 = .FALSE.
         ELSEIF (I2TYP.EQ.3.AND.FIRST3) THEN
            DOFIRST = .TRUE.
            FIRST3 = .FALSE.
         ELSEIF (I2TYP.EQ.4.AND.FIRST4) THEN
            DOFIRST = .TRUE.
            FIRST4 = .FALSE.
         ENDIF
      ENDIF
C
      IF (NEWBAS .OR. DOFIRST) THEN
         ! make a list of tasks in INDEX_IJ and a task weight for each in TIMES
         CALL HERPAR_DET_TASK(I2TYP,2,MTOTTK,INDEX_IJ,TIMES)
         ! make a sorted list of tasks in NPOS based on the weights in TIMES
         CALL PARSR2(MTOTTK,NPOS,TIMES,INDEX_IJ,IPRPAR)
      ENDIF
C
C     Distribute jobs on nodes
C
      ICOUNT = 0
      IBLOCK = NODES*NTASK
      DO I = MTOTTK, 1, -1
        MORE1  = (ICOUNT/IBLOCK)*IBLOCK
        MORE2  = MOD(ICOUNT,NODES)*NTASK
        MORE3  = ((ICOUNT-MORE1)/NODES) + 1
        IPLACE = MORE1+MORE2+MORE3
        INDEX_IJ(IPLACE) = NPOS(I)
        ICOUNT=ICOUNT+1
      ENDDO
C
C     Determine number of batches (NBATCH).
C
      NOTALL = (MTOTTK/IBLOCK)*IBLOCK
      NREST  = MTOTTK - NOTALL
      IF (NREST .GT. NODES) NREST = NODES
      NBATCH = NOTALL/NTASK + NREST
C
C     Print current settings.
C
      NEACH = NBATCH/NODES
      NREST = NBATCH - NEACH*NODES
C
      IF (FIRST .AND. IPRPAR .GT. 2) THEN
         WRITE(LUPRI,'(/A,5(/5X,A,I5),2X,A7,I3,A1/)')
     &        'Output from HERPAR_INDEKS:',
     &        'Number of nodes   :',NODES,
     &        'Number of tasks   :',MTOTTK,
     &        'Tasks pr. batch   :',NTASK,
     &        'Number of batches :',NBATCH,
     &        'Batches pr. node  :',NEACH,
     &        '(rest =',NREST,')'
      END IF
C
C     Print the INDEX_IJ array.
C
      IF (IPRPAR .GE. 10) THEN
         WRITE(LUPRI,'(/5X,A/)') 'The order of task-distribution :'
         DO 400 I=1, NBATCH*NTASK
            WRITE(LUPRI,'(5X,A3,I5,A15,I5)')
     &           'i =',I,'index_ij(i) =',INDEX_IJ(I)
 400      CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck sdinit */
      SUBROUTINE SDINIT(DMAT,NDMAT,IREPDM,IFCTYP,ITYPE,MAXDIF,IATOM,
     &                  NODV,NOPV,NOCONT,TKTIME,RETUR,I2TYP,
     &                  ICEDIF,SCREEN,GABRAO,DMRAO,DMRSO)
C
C     Written by Henrik Koch and Paal Dahle
C
      use interface_to_mpi
#include "implicit.h"
C
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#else
#include "gnrinf.h"
#endif
#include "molinp.h"
#include "aovec.h"
C
#include "hrunit.h"
#include "dorps.h"
#include "abainf.h"
#include "blocks.h"
#include "symmet.h"
#include "nuclei.h"
#include "infpar.h"
#include "parint.h"
C
#include "inforb.h"
#include "cbisol.h"
#include "cbirea.h"
#include "siripc.h"
C
      LOGICAL NODMAT,NOFMAT
      DIMENSION DMAT(*), IREPDM(NDMAT), IFCTYP(NDMAT),
     &          DMRAO(NSYMBL,NSYMBL,*),GABRAO(NSYMBL,NSYMBL,*),
     &          DMRSO(MAXSHL,MAXSHL,*), DSEND(2), ISEND(2)
      LOGICAL NODV,NOPV,NOCONT,TKTIME,RETUR, LSEND(6), PRIMUPD_NODES
C
C--------------------------------------
C     Calculate number of Fock matrices
C--------------------------------------
C
      IF (ITYPE .EQ. 2 .OR. ITYPE .EQ. 4) THEN
         NFMAT  = 0
      ELSE IF (ITYPE .EQ. 3 .OR. ITYPE .EQ. 9 .OR. ITYPE .EQ. 7) THEN
         NFMAT  = NDMAT
      ELSE IF (ITYPE .EQ. -5) THEN
         IF(MAXDIF.EQ.1) THEN
           NFMAT  = 3*NDMAT
         ELSEIF(MAXDIF.EQ.2) THEN
           NFMAT  = 6*NDMAT
         ENDIF
      ELSE IF (ITYPE .EQ. 8) THEN
         NFMAT  = 3*NDMAT*NUCDEG(IATOM)
      ELSE
         WRITE (LUPRI,'(/A,I5)')
     &' SDINIT error; specified ITYPE for TWOINT not defined. ITYPE =',
     &   ITYPE
         CALL QUIT('Specified ITYPE for TWOINT not defined in SDINIT.')
      END IF
C
C----------------------------
C     Set common-block PARINT
C----------------------------
C
      JCODE  = 0
C     ... not used any more (was used for PVM)
      JATOM  = IATOM
      JLUDAS = LUDASP
      JLUINT = LUINTA
      JLUONE = LUONEL
      JLUSOL = LUSOL
      JLUSUP = LUSUPM
      JMXDIF = MAXDIF
      JMXREP = MAXREP
      JNODES = NODES
      JPRINT = IPRPAR
      JTASK  = NTASK
      JTYPE  = ITYPE

      JNODV  = 0
      IF (NODV)   JNODV  = 1

      JNOPV  = 0
      IF (NOPV)   JNOPV  = 1

      JNOCNT = 0
      IF (NOCONT) JNOCNT = 1

      JRETUR = 0
      IF (RETUR)  JRETUR = 1

      JTKTIM = 0
      IF (TKTIME) JTKTIM = 1

      JTIMIN = 0
      IF (TIMING) JTIMIN = 1

      JSOLVN = 0
      IF (SOLVNT) JSOLVN = 1

      JRELCL = 0
      IF (RELCAL) JRELCL = 1

      JLEVYL = 0
      IF (LEVYLE) JLEVYL = 1

      JNONREL = 0
      IF (NONREL) JNONREL = 1

      JTWOCOMP = 0
      IF (TWOCOMP) JTWOCOMP = 1

      JECPCALC = 0
      IF (ECPCALC) JECPCALC = 1

      J2TYP  = I2TYP
      JCEDIF = ICEDIF
C
      CALL IZERO(JDOSYM,8)
      CALL IZERO(JDOREP,8)
      CALL IZERO(JDOCOR,MXCOOR)
C
      DO 100 I = 0, MAXREP
         IF (DOREPS(I))  JDOREP(I+1) = 1
         IF (DOSYM(I+1)) JDOSYM(I+1) = 1
 100  CONTINUE
C
      DO 110 J = 1, MXCENT
         DO 120 I = 1, 3
            IF (DOCOOR(I,J)) JDOCOR(I,J) = 1
 120     CONTINUE
 110  CONTINUE
C
#if defined (VAR_MPI)
      PRIMUPD_NODES = .FALSE.
      IF (NEWGEO) THEN
         CAll interface_mpi_BCAST(NMLINE,1,MPARID,
     &       global_communicator)

C........inforb.h
         NUMELM=ICOMMSIZE(I1_INFORB,I2_INFORB)
         CAll interface_mpi_BCAST(NUMELM,1,MPARID,
     &               global_communicator)
C........note that with icommsize COMMON blocks must be sent as integers...
         CAll interface_mpi_BCAST(I1_INFORB,NUMELM,MPARID,
     &                  global_communicator)
C
         CAll interface_mpi_bcast_l0(UNCONT,1,MPARID,
     &                  global_communicator)
         CAll interface_mpi_BCAST(MAXPRI,1,MPARID,
     &                  global_communicator)
         NEWGEO = .FALSE.
         PRIMUPD_NODES = .TRUE.
      ELSE
         NFINIS = 0
         CAll interface_mpi_BCAST(NFINIS,1,MPARID,
     &       global_communicator)
      END IF
C
C     Send Panas value and screening threshold

      DSEND(1) = PANAS
      DSEND(2) = SCREEN
      CAll interface_mpi_BCAST(DSEND,2,MPARID,
     &               global_communicator)
#if defined (PRG_DIRAC)
C     Send SMLV1C and ONECAP,INTV1C
      LSEND(1) = SMLV1C
      LSEND(2) = ONECAP
      CAll interface_mpi_bcast_l0(LSEND,2,MPARID,
     &               global_communicator)
      CAll interface_mpi_BCAST(INTV1C,1,MPARID,
     &               global_communicator)
#endif
C
C     Send parint.h
      NUMELM = ICOMMSIZE(I1_PARINT,I2_PARINT)
      CAll interface_mpi_BCAST(NUMELM,1,MPARID,
     &               global_communicator)
C........note that with icommsize COMMON blocks must be sent as integers...
      CAll interface_mpi_BCAST(I1_PARINT,NUMELM,MPARID,
     &               global_communicator)
C
C     ... otherwise xyz reader fails in parallel mode
C
      IF( PRIMUPD_NODES )THEN
#if defined (VAR_MPI)
        CALL GET_PRIMITF()
#endif
      END IF
      ISEND(1) = NDMAT
      ISEND(2) = NFMAT
      CAll interface_mpi_BCAST(ISEND,2,MPARID,global_communicator)
      IF(NDMAT.GT.0.AND.ITYPE.NE.4) THEN
        CAll interface_mpi_bcast_r1_work_f77(DMAT,N2BASX*NDMAT,MPARID,
     &               global_communicator)
      ENDIF
C
      IF(NFMAT.GT.0) THEN
        CAll interface_mpi_BCAST(IREPDM,NDMAT,MPARID,
     &               global_communicator)
      ENDIF
      IF(ITYPE.NE.4) THEN
        IF (ITYPE.EQ.2) THEN
          NFCTYP = NDMAT
        ELSE
          NFCTYP = NFMAT
        ENDIF
        CAll interface_mpi_BCAST(IFCTYP,NFCTYP,MPARID,
     &                 global_communicator)
      END IF
      CAll interface_mpi_BCAST(HFXFAC,1,MPARID,
     &               global_communicator)
      CAll interface_mpi_BCAST(HFXATT,1,MPARID,
     &               global_communicator)
      CAll interface_mpi_BCAST(HFXMU,1,MPARID,
     &               global_communicator)
      call interface_mpi_bcast_l0(spinfr, 1, mparid,
     &               global_communicator)



C
C     ITYPE = 3 : direct SCF (AO-basis)
C     ITYPE = 9 : direct SCF (SO-basis)
C     ITYPE = 2 : gradient
C     ITYPE = 4 : (rs|**)-distributions
C
      IF(SCREEN.GT.-1.0D0) THEN
        IF(ITYPE.EQ.3.OR.ITYPE.EQ.9.OR.ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
C         ... broadcast screening matrix
          N2RED = NSYMBL*NSYMBL
          CAll interface_mpi_BCAST(N2RED,1,MPARID,
     &                   global_communicator)
          CAll interface_mpi_bcast_r1_work_f77(DMRAO,N2RED*NDMAT,MPARID,
     &                   global_communicator)
          IF (ITYPE .EQ. 2) THEN
             CAll interface_mpi_bcast_r1_work_f77(GABRAO,4*N2RED,MPARID,
     &                      global_communicator)
          ELSE
             CAll interface_mpi_bcast_r1_work_f77(GABRAO,N2RED,MPARID,
     &                      global_communicator)
          END IF
        END IF
        IF(ITYPE.EQ.9) THEN
          N2RED = MAXSHL*MAXSHL
          CAll interface_mpi_BCAST(N2RED,1,MPARID,
     &                   global_communicator)
          CAll interface_mpi_bcast_r1_work_f77(DMRSO,N2RED*NDMAT,MPARID,
     &                   global_communicator)
        END IF
      ENDIF
#endif
      RETURN
      END
C  /* Deck do_herpar */
      SUBROUTINE DO_HERPAR(INDEX_IJ,MTOTTK,NBATCH,NSTAT)
C
C     Written by Henrik Koch and Paal Dahle Nov-1993
C
      use interface_to_mpi
#include "implicit.h"
#ifdef VAR_MPI
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#endif
#include "infpar.h"
#include "priunit.h"
#include "maxorb.h"
C
      DIMENSION INDEX_IJ(MTOTTK), NSTAT(NODES), ISHELL(2)
C

      CALL IZERO(NSTAT,NODES)
C
C------------------------------------------
C     Start loop over two outermost shells.
C------------------------------------------
C
      NTEST  = 1
      NUMTSK = 0
      IPLACE = 0
C
      DO 100, I = 1,NBATCH
C
C        Ready to receive request from any node.
C
#if defined (VAR_MPI)


         CAll interface_mpi_RECV(NODE,1,df_MPI_ANY_SOURCE,20,
     &                 global_communicator)
C
         CAll interface_mpi_SEND(NTEST,1,NODE,30,
     &                 global_communicator)

#endif
CSK         write(6,*)'send and recieve done',mytid
C
         DO 200 J = 1,NTASK
            IPLACE = IPLACE + 1
            NUMTSK = NUMTSK + 1
            ISHELL(1) = INDEX_IJ(IPLACE)
            ISHELL(2) = IPLACE
C
#if defined (VAR_MPI)
            CAll interface_mpi_SEND(ISHELL,2,NODE,30,
     &                    global_communicator)
#endif
C
C           Pack only one null-task.
C
            IF (INDEX_IJ(IPLACE) .EQ. 0) THEN
               IPLACE = IPLACE + NTASK - J
               NUMTSK = NUMTSK - 1
               GO TO 300
            END IF
C
            IF (IPRPAR .GT. 10) THEN
               WRITE(LUPRI,'(5X,I5,A10,I5,A8,I4)')
     &         NUMTSK, 'Send task', INDEX_IJ(IPLACE), 'to node', NODE
            END IF
  200    CONTINUE
C
  300    CONTINUE
C
C
C        Who-got-the-last-batch statistics
C
         NSTAT(NODE) = NSTAT(NODE) + 1
C
  100 CONTINUE
C
C------------------------------------------------
C     Instruct nodes to return the final results.
C------------------------------------------------
C
      NTEST = -1
      DO I = 1,NODES
C
C        Ready to receive request from any node.
C
#if defined (VAR_MPI)
         CAll interface_mpi_RECV(NODE,1,df_MPI_ANY_SOURCE,20,
     &                 global_communicator)
C
         CAll interface_mpi_SEND(NTEST,1,NODE,30,
     &                 global_communicator)
#endif
C
      END DO ! I
      RETURN
      END
C  /* Deck rvres */
      SUBROUTINE RVRES(WORK,FMAT,ITYPE,MTOTTK,NREQUI,
     &                 NUMGRA,NUMHES,MAXDIF,TIMES,DINTSKP,LDINTSKP)
C
C           ITYPE  - Calc. type: 2 = expectation values.
C                                3 = direct SCF (AO-basis)
C                                4 = (rs|**)-distributions
C                               -5 = Integrals derivated with respect
C                                    to magnetic field.
C                                7 = direct CC (DALTON)
C                                8 = Fock matrices derivated with respect
C                                    to atomic coordinates.
C                                9 = direct SCF (SO-basis)
C     Written by Henrik Koch and Paal Dahle July-1994
C
C   MI: Added TEC changes....
C
      use interface_to_mpi
#include "implicit.h"
#if defined (VAR_MPI)
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#endif
      PARAMETER (D1 = 1.0D0)
C
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
#include "suscpt.h"
#include "inforb.h"
#include "infpar.h"

C
      REAL*8    FMAT(*),WORK(NREQUI),TIMES(MTOTTK),DINTSKP(LDINTSKP),
     &          DNSKBF(3*3*14),SUS2ELBUF(9)

C
C     Clean matrices.
C
      IF (ITYPE .EQ. 2) THEN
         CALL DZERO(GRADEE,NUMGRA)
         IF (MAXDIF .GT. 1)
     &      CALL DZERO(HESSEE,NUMHES)
      ELSE IF (ITYPE .EQ. -5) THEN
         IF(MAXDIF.EQ.2) CALL DZERO(SUS2ELBUF,9)
      ELSE IF (ITYPE .EQ. 8) THEN
         CALL DZERO(FMAT,3*N2BASX)
      END IF
C
C---------------------------
C    use MPI_REDUCE for type 3 ( direct SCF)
C---------------------------
C
#if defined (VAR_MPI)
      IF (ITYPE .EQ. 3) THEN
         ! using _f77 version
         CALL dcopy(NREQUI,FMAT,1,WORK,1)
         CALL interface_mpi_reduce_r1_work_f77(
     &        WORK, FMAT, NREQUI,
     &        op_MPI_SUM, MPARID, global_communicator)
         call dcopy(ldintskp,DINTSKP,1,WORK,1)
         CALL interface_mpi_reduce_r1_work_f77(
     &        WORK, DINTSKP, LDINTSKP,
     &        op_MPI_SUM, MPARID, global_communicator)
         IF (TIMING) THEN
            call dzero(WORK,MTOTTK)
            CALL interface_mpi_reduce_r1_work_f77(
     &        WORK, TIMES, MTOTTK,
     &        op_MPI_SUM, MPARID, global_communicator)
         END IF
         GO TO 9000 ! RETURN
      END IF
C
C---------------------------
C     Start loop over nodes.
C---------------------------
C
      ICOUNT = 0
      DO 100 I = 1,NODES
C
         CALL interface_mpi_RECV(ITASK,1,df_mpi_any_source,40,
     &                 global_communicator,ISTAT)
C
         IF (ITASK .EQ. 0) GOTO 100 ! more nodes than tasks

         ICOUNT = ICOUNT + ITASK
C
         NWHO = ISTAT(df_MPI_SOURCE)

C
         IF(ITYPE.EQ.3.OR.ITYPE.EQ.9.OR.ITYPE.EQ.2) THEN
             CAll interface_mpi_RECV(DNSKBF,LDINTSKP,
     &           NWHO,40,global_communicator)
         ENDIF
C
Cjth & trond - this is a bit risky
C
         IF(.NOT.(ITYPE.EQ.-5.AND.MAXDIF.EQ.2)) THEN
           CALL interface_mpi_RECV(WORK,NREQUI,NWHO,40,
     &                 global_communicator)
         ENDIF
C
C        Add results to final matrices.
C
         IF (ITYPE .EQ. 2) THEN
            CALL DAXPY(NUMGRA,D1,WORK,1,GRADEE,1)
            DO J = 1,LDINTSKP
              DINTSKP(J) = DINTSKP(J) + DNSKBF(J)
            ENDDO
            IF (MAXDIF .GT. 1)
     &         CALL DAXPY(NUMHES,D1,WORK(NUMGRA+1),1,HESSEE,1)
C
         ELSE IF (ITYPE .EQ. 3.OR.ITYPE.EQ.9) THEN
            CALL DAXPY(NREQUI,D1,WORK,1,FMAT,1)
            DO J = 1,8
              DINTSKP(J) = DINTSKP(J) + DNSKBF(J)
            ENDDO
C
         ELSE IF (ITYPE .EQ. 8) THEN
            CALL DAXPY(NREQUI,D1,WORK,1,FMAT,1)
         ELSE IF (ITYPE .EQ. -5) THEN
            IF(MAXDIF.EQ.1) THEN
              CALL DAXPY(NREQUI,D1,WORK,1,FMAT,1)
            ELSEIF(MAXDIF.EQ.2) THEN
              CALL interface_mpi_RECV(SUS2ELBUF,9,NWHO,40,
     &                global_communicator)
              CALL DAXPY(9,D1,SUS2ELBUF,1,SUS2EL,1)
            ENDIF
         END IF
C
C        The time each index-pair calculation need is collected in TIMES
C
         IF (TIMING) THEN
            DO 300, K = 1, ITASK
C
               CAll interface_mpi_RECV(NUMBER,1,NWHO,40,
     &                       global_communicator)
               CAll interface_mpi_RECV(VALUE,1,NWHO,40,
     &                       global_communicator)
C
               TIMES(NUMBER) = VALUE
C
  300       CONTINUE
         END IF
  100 CONTINUE
#endif
C
      IF (ICOUNT .NE. MTOTTK) THEN
        WRITE(LUPRI,'(/5X,A)') 'Error in parallel calculation!'
        WRITE(LUPRI,'(5X,A,I5)') 'Number of tasks to calculate  : ',
     &       MTOTTK
        WRITE(LUPRI,'(5X,A,I5)') 'Number of tasks    calculated : ',
     &       ICOUNT
        CALL QUIT('Inconsistency in RVRES : MTOTTK .NE. ICOUNT')
      ENDIF
C
 9000 RETURN
      END
C  /* Deck parsrt */
      SUBROUTINE PARSRT(MTOTTK,NPOS,SORTED,TIMES)
C
C     Written by Paal Dahle December 1993.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "infpar.h"
C
      DIMENSION  TIMES(MTOTTK), NPOS(0:MTOTTK+1), SORTED(0:MTOTTK+1)
C
C------------------------------------------------------
C     The array TIMES containing CPU-times is sorted by
C     increasing CPUs and stored in the array SORTED.
C------------------------------------------------------
C
C     SORTED(0) must be less than any possible CPU-time,
C     and TIMMAX must be greater than any possible CPU-time
C
      ITEMS     = MTOTTK
      TIMMAX    = 99999D00
      SORTED(0) =  -1.0D00
      SORTED(1) = TIMES(1)
      NPOS  (0) = 0
      NPOS  (1) = 1
C
      DO 100,ITEM = 2,ITEMS
C
         LASTLO       = 0
         LASTHI       = ITEM
         NUMBER       = INT(ITEM/2)
         SORTED(ITEM) = TIMMAX
C
 200     CONTINUE
C
         IF (TIMES(ITEM) .EQ. SORTED(NUMBER) .OR.
     &      (TIMES(ITEM) .GT. SORTED(NUMBER) .AND.
     &       TIMES(ITEM) .LT. SORTED(NUMBER+1))) THEN
C
            DO 300, NCOUNT = ITEM,NUMBER,-1
               NPOS  (NCOUNT+1) = NPOS  (NCOUNT)
               SORTED(NCOUNT+1) = SORTED(NCOUNT)
 300        CONTINUE
C
            NPOS  (NUMBER+1) = ITEM
            SORTED(NUMBER+1) = TIMES(ITEM)
            GO TO 100
C
         ELSE IF (TIMES(ITEM) .LT. SORTED(NUMBER)) THEN
C
            LASTHI = NUMBER
            NUMBER = INT((LASTLO+LASTHI)/2)
            GO TO 200
C
         ELSE IF (TIMES(ITEM) .GT. SORTED(NUMBER)) THEN
C
            LASTLO = NUMBER
            NUMBER = INT((LASTLO+LASTHI)/2)
            GO TO 200
C
         END IF
 100  CONTINUE
C
      IF (IPRPAR .GE. 9) THEN
         DO 400 I = 1, ITEMS
            WRITE(LUPRI,'(15X,A10,I5,A9,F8.3)')
     &           'IJ-shell =', NPOS(I), '   time =', SORTED(I)
 400      CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck timres */
      SUBROUTINE TIMRES(NSTAT)
C
C     Written by Henrik Koch and Paal Dahle Nov-1993
      use interface_to_mpi
#include "implicit.h"
      PARAMETER ( D0 = 0.00D00 )
#if defined (VAR_MPI)
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#endif
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "infpar.h"
C
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#endif
C
      INTEGER   NSTAT(NODES)

      REAL*8    RBUF(4), TTOT, WTOT
      CHARACTER SECTID*12, INTTID*12

      REAL*8, SAVE        :: TSOFAR=0.0D0, WSOFAR=0.0D0
      REAL*8, ALLOCATABLE :: CPUS(:,:),WALLS(:,:)
C
      IF (NODES .LE. 0) RETURN

#if defined (VAR_MPI)

      allocate( CPUS(2,NODES), WALLS(2,NODES) )
#ifdef PRG_DIRAC
      IF (RELGRD) DMXCPU = D0
      IF (RELGRD) DMXWLL = D0
#endif
      TTOT = 0.0D0
      WTOT = 0.0D0
      DO I = 1,NODES
         CAll interface_mpi_RECV(RBUF,4,df_mpi_any_source,
     &                 50,global_communicator,istat)
         NODE = ISTAT(df_MPI_SOURCE)
         CPUS(1,NODE)  = RBUF(1)
         CPUS(2,NODE)  = RBUF(2)
         WALLS(1,NODE) = RBUF(3)
         WALLS(2,NODE) = RBUF(4)
C
         TTOT = TTOT + CPUS(1,NODE)
         WTOT = WTOT + WALLS(1,NODE) + WALLS(2,NODE)
#ifdef PRG_DIRAC
         IF ( RELGRD ) THEN
            DMXCPU = MAX(DMXCPU,CPUS(1,NODE))
            DMXWLL = MAX(DMXWLL,WALLS(2,NODE))
         END IF
#endif
C
      END DO ! I
      TSOFAR = TSOFAR + TTOT
      WSOFAR = WSOFAR + WTOT
C
      IF (IPRPAR .GT. 2) THEN
         CALL AROUND('Overall statistics for distribution of batches')
C
         IF (IPRPAR .GT. 4 .OR. NODES .le. 10) THEN
            NODES_print = NODES
         ELSE
            NODES_print = 10
            WRITE(LUPRI,'(11X,A/)') '(Only first 10 nodes printed).'
         END IF

         WRITE(LUPRI,'(11X,7(A,4X))')   'Node',
     &                                  'Batches',
     &                                  'CPU time',
     &                                  'Wall time',
     &                                  'Efficiency',
     &                                  'Wall time w/o MPI',
     &                                  'Efficiency w/o MPI'
         WRITE(LUPRI,'(11X,7(A,4X)/)')  '----',
     &                                  '-------',
     &                                  '--------',
     &                                  '---------',
     &                                  '----------',
     &                                  '-----------------',
     &                                  '------------------'
C
         DO NODE = 1, NODES_print
            WALLS_TOT = WALLS(1,NODE)+WALLS(2,NODE)
            IF (WALLS_TOT .NE. 0.0D0) THEN
               EFFI_T = (CPUS(1,NODE)/WALLS_TOT)*100.D0
            ELSE
               EFFI_T = 0.0D0
            END IF
            IF (WALLS(1,NODE) .NE. 0.0D0) THEN
               EFFI_1 = (CPUS(1,NODE)/WALLS(1,NODE))*100.D0
            ELSE
               EFFI_1 = 0.0D0
            END IF
            WRITE(LUPRI,'(I14,I10,F12.2,F13.2,F12.2,2F20.2)')
     &         NODE,NSTAT(NODE),CPUS(1,NODE),WALLS_TOT,EFFI_T,
     &         WALLS(1,NODE),EFFI_1
         END DO ! NODE
         WRITE(LUPRI,'(/11X,A)')
     &   'Note: CPU times are without CPU time used in MPI.'
      END IF
C
#ifdef PRG_DIRAC
      IF ( RELGRD ) THEN
         DMXWLT = DMXWLL
         DMXCPT = DMXCPU
      END IF
#endif
C
      INTTID = SECTID(TTOT)
      WRITE(LUPRI,'(/A,A12,A)')
     &     '>>>> Total CPU  time used in all Hermit nodes this call :',
     &     INTTID, ' (MPI wait excluded)'
      INTTID = SECTID(WTOT)
      WRITE(LUPRI,'(A,A12/A,F12.2)')
     &     '>>>> Total WALL time used in all Hermit nodes this call :',
     &     INTTID,
     &     '>>>> Efficiency this call :', (TTOT/WTOT) * 100.0D0

      INTTID = SECTID(TSOFAR)
      WRITE(LUPRI,'(/A,A12,A)')
     &     '>>>> Total CPU  time used in all Hermit nodes so far    :',
     &     INTTID, ' (MPI wait excluded)'
      INTTID = SECTID(WSOFAR)
      WRITE(LUPRI,'(A,A12/A,F12.2)')
     &     '>>>> Total WALL time used in all Hermit nodes so far    :',
     &     INTTID,
     &     '>>>> Efficiency so far    :', (TSOFAR/WSOFAR) * 100.0D0

C
      deallocate( CPUS, WALLS )
#endif    /* VAR_MPI */
      RETURN
      END
C  /* Deck hernod */
      SUBROUTINE HERNOD()
C
C     Written by Henrik Koch and Paal Dahle July-1994
C
C    *****************************************************************
C    *    This is the node program for the construction of fock-     *
C    *   matrices, derivated fock-matrices and expectation values.   *
C    *****************************************************************
C
CMI called from DIRNOD,DIRNOD2 (main/dirac.F)
C
CMI  TEC changes for LAO (2003,2006)

      use memory_allocator
      use interface_to_mpi
#include "implicit.h"
C
C
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "suscpt.h"
#include "inforb.h"
#include "infpar.h"
#include "molinp.h"
#include "blocks.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcbgrd.h"
#else
#include "gnrinf.h"
#include "energy.h"
#endif
C
      DIMENSION DINTSKP(14*3*3)
      LOGICAL   NODV,NOPV,NOCONT,RETUR,TKTIME,NEWGEO,FINISH
C
      real(8), allocatable :: WORK(:)
      call legacy_lwork_get(LWORK)
      call alloc(WORK,LWORK,id='WORK in HERNOD')
C
      CALL QENTER('HERNOD')
#include "memint.h"
C
      PARHER = .TRUE.
C
      CALL DZERO(DINTSKP,3*3*14)
C
 100  CONTINUE
C
C     Start timing.
C
      CALL GETTIM(CPU1,WALL1)
C
C
C     ***************************************************
C     ***** Receive initialization data from master *****
C     ***************************************************
C

      CALL RVINIT(WORK,KFREE,LFREE,NDMAT,NFMAT,ITYPE,
     &            IATOM,MAXDIF,MYTID,NODV,NOPV,NOCONT,RETUR,
     &            TKTIME,NEWGEO,FINISH,I2TYP,ICEDIF,SCREEN,
     &            KFMAT,KDMAT,KIFC,KIRD,KGABAO,KDMRAO,KDMRSO)
C
      IF (FINISH) GOTO 999
C
C     ***************************
C     ***** Allocate memory *****
C     ***************************
C
C

      IF (ITYPE .EQ. 2) THEN
         NUMGRA = MXCOOR
         IF (MAXDIF .GT. 1) THEN
            NUMHES = NUMGRA*NUMGRA
            NREQUI = NUMGRA + NUMHES
         ELSE
            NUMHES = 0
            NREQUI = NUMGRA
         END IF
         LFMAT  = 0
      ELSE IF (ITYPE .EQ.3.OR.ITYPE.EQ.9) THEN
         NREQUI = 0
         LFMAT  = NFMAT*N2BASX
      ELSE IF (ITYPE .EQ. -5) THEN
         LFMAT  = NFMAT*N2BASX
      ELSE IF (ITYPE .EQ. 8) THEN
         NREQUI = 3*NNBASX
         LFMAT  = 0
      ELSE
         WRITE (LUPRI,*) ' ERROR in HERNOD: unknown ITYPE:',ITYPE
         CALL QUIT('Unknown ITYPE in HERNOD')
      END IF
C
C     Determine total number of tasks - MTOTTK.
      CALL HERPAR_DET_TASK (I2TYP,1,MTOTTK,NPOS,WORK(KFREE))
C
      CALL MEMGET2('INTE','INDEX',KINDEX,MTOTTK,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','TIMES',KTIMES,MTOTTK,WORK,KFREE,LFREE)
C
C     ********************************
C     ***** Initialize workspace *****
C     ********************************
C
      CALL DZERO(WORK(KFMAT),LFMAT)
      IF (ITYPE .EQ. 2) THEN
         CALL DZERO(GRADEE,NUMGRA)
         IF (MAXDIF.GT.1) CALL DZERO(HESSEE,NUMHES)
      END IF
C
C     ***********************************
C     ***** Calculate the integrals *****
C     ***********************************
C
      GMAT   = 0
      JPRINT = 0
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      ISHLA  = 0
C
      CALL TWOINT(WORK(KFREE),LFREE,WORK(KFMAT),WORK(KDMAT),NDMAT,
     &     WORK(KIRD),WORK(KIFC),GMAT,IDUMMY,IDUMMY,ITYPE,
     &     MAXDIF,IATOM,NODV,NOPV,NOCONT,TKTIME,JPRINT,
     &     IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,ISHLA,I2TYP,
     &     ICEDIF,SCREEN,WORK(KGABAO),
     &     WORK(KDMRAO),WORK(KDMRSO),DINTSKP,RELCAL,.false.,
     &     WORK(KINDEX),WORK(KTIMES))
C
C     ****************************************
C     ***** Send final results to master *****
C     ****************************************
C

      CALL SDRES(WORK(KFMAT),LFMAT,NUMGRA,MAXDIF,
     &     NUMHES,NUMSUS,CPU1,WALL1,ITYPE,DINTSKP,
     &     MTOTTK,WORK(KINDEX),WORK(KTIMES))

C
C     We want the slave program to go until there is nothing
C     more to calculate, K.Ruud, May-95
C
 20   CONTINUE
      CALL MEMREL('HERNOD',WORK,1,KDMAT,KFREE,LFREE)
      GOTO 100
C
 999  CONTINUE
C
      call dealloc(WORK)
      PARHER = .FALSE.
      CALL QEXIT('HERNOD')
      END
C  /* Deck rvinit */
      SUBROUTINE RVINIT(WORK,KFREE,LFREE,NDMAT,NFMAT,ITYPE,
     &                  IATOM,MAXDIF,MYTI,NODV,NOPV,NOCONT,RETUR,
     &                  TKTIME,NEWGEO,FINISH,I2TYP,ICEDIF,SCREEN,
     &                  KFMAT,KDMAT,KIFC,KIRD,KGABAO,KDMRAO,KDMRSO)
C
C     Written by Henrik Koch and Paal Dahle
C     Revised Oct 4 1996 - tsaue
C
      use interface_to_mpi
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "hrunit.h"
#include "dorps.h"
#include "abainf.h"
#include "symmet.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#else
#include "gnrinf.h"
#endif
#include "molinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "parint.h"
C
#include "inforb.h"
#include "cbisol.h"
#include "cbirea.h"
C
      DIMENSION WORK(*), DSEND(2), ISEND(2)
      CHARACTER*6 NAVN
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME,NEWGEO,FINISH, LSEND(6)

C
C-----------------------------------------
C     Receive initializations from master.
C-----------------------------------------
C
#if defined (VAR_MPI)

      CAll interface_mpi_BCAST(NMLINE,1,MPARID,global_communicator)
      FINISH = .FALSE.
      IF (NMLINE .LT. 0) THEN
         FINISH = .TRUE.
         RETURN
      ELSE IF (NMLINE .GT. 0) THEN

C........inforb.h
         CAll interface_mpi_BCAST(NUMELM,1,MPARID,
     &                  global_communicator)
         CAll interface_mpi_BCAST(I1_INFORB,NUMELM,MPARID,
     &                  global_communicator)
C........
         CAll interface_mpi_bcast_l0(UNCONT,1,MPARID,
     &                  global_communicator)
         CAll interface_mpi_BCAST(MAXPRI,1,MPARID,
     &                  global_communicator)
         NEWGEO = .TRUE.
      ELSE
        NEWGEO = .FALSE.
      END IF


C     Receive PANAS and screening threshold
      CAll interface_mpi_BCAST(DSEND,2,MPARID,
     &               global_communicator)
      PANAS  = DSEND(1)
      SCREEN = DSEND(2)
#if defined (PRG_DIRAC)
C     Receive SMLV1C and ONECAP,INTV1C
      CAll interface_mpi_bcast_l0(LSEND,2,MPARID,
     &               global_communicator)
      SMLV1C = LSEND(1)
      ONECAP = LSEND(2)
      CAll interface_mpi_BCAST(INTV1C,1,MPARID,
     &               global_communicator)
#endif
C.....receive parint.h
      CAll interface_mpi_BCAST(NUMELM,1,MPARID,
     &               global_communicator)
      CAll interface_mpi_BCAST(I1_PARINT,NUMELM,MPARID,
     &               global_communicator)
#endif
C
C------------------------------------
C     Retain information from PARINT.
C------------------------------------
C
      IATOM  = JATOM
      ITYPE  = JTYPE
      IPRPAR = JPRINT
      LUDASP = JLUDAS
      LUINTA = JLUINT
      LUONEL = JLUONE
      LUSOL  = JLUSOL
      LUSUPM = JLUSUP
      MAXDIF = JMXDIF
      MAXREP = JMXREP
      NODES  = JNODES
      NTASK  = JTASK
      NOCONT = JNOCNT .EQ. 1
      NODV   = JNODV  .EQ. 1
      NOPV   = JNOPV  .EQ. 1
      RETUR  = JRETUR .EQ. 1
      TIMING = JTIMIN .EQ. 1
      TKTIME = JTKTIM .EQ. 1
      SOLVNT = JSOLVN .EQ. 1
      I2TYP  = J2TYP
      RELCAL = JRELCL .EQ. 1
      TWOCOMP = JTWOCOMP .EQ. 1
      LEVYLE  = JLEVYL .EQ. 1
      NONREL  = JNONREL .EQ. 1
      ECPCALC = JECPCALC .EQ. 1
      ICEDIF  = JCEDIF
C
      DO 200 I = 0, MAXREP
         DOREPS(I)  = JDOREP(I+1) .EQ. 1
         DOSYM(I+1) = JDOSYM(I+1) .EQ. 1
 200  CONTINUE
      DO 210 J = 1, MXCENT
         DO 220 I = 1, 3
            DOCOOR(I,J) = JDOCOR(I,J) .EQ. 1
 220     CONTINUE
 210  CONTINUE
C
C     Set hermit.
C
      IF (NEWGEO) RDINPC = .FALSE.
C
C     call SETHER to read basis set input and set hermit variables in common blocks for this slave
C
      CALL SETHER(0,NEWGEO,RELCAL,WORK(KFREE),LFREE)
C
      IF (NEWGEO .AND. SOLVNT) THEN
         NEWGEO = .TRUE.
         NUCIND = NUCIND + 1
         NUCDEP = NUCDEP + 1
         NATOMS = NATOMS + 1
         NCNTCV = NUCIND
         NAVN(1:4) = 'cav '
         NCLINE(NUCIND) = 0
         NAMEX(3*NUCIND)     = NAVN(1:4)//' z'
         NAMEX(3*NUCIND - 1) = NAVN(1:4)//' y'
         NAMEX(3*NUCIND - 2) = NAVN(1:4)//' x'
         NAMDEP(NUCDEP)     = NAMEX(3*NUCIND)(1:4)//'  '
         NAMDPX(3*NUCDEP-2) = NAMEX(3*NUCIND-2)
         NAMDPX(3*NUCDEP-1) = NAMEX(3*NUCIND-1)
         NAMDPX(3*NUCDEP  ) = NAMEX(3*NUCIND  )
         IF (NUCDEP .GT. MXCENT) THEN
            WRITE (LUPRI,'(//2A,/A,I5)')
     &         ' Too many atomic centers: MXCENT exceed in READIN for',
     &         ' solvent cavity,',' Current limit:',MXCENT
            CALL QUIT('*** ERROR *** MXCENT exceeded in READIN')
         END IF
         READ (NAVN,'(A4)') NAMN(NUCIND)
         CORD(1,NUCIND) = D0
         CORD(2,NUCIND) = D0
         CORD(3,NUCIND) = D0
         ISTBNU(NUCIND) = 7
         CHARGE(NUCIND) = D0
         CALL NUCPRO(WORK(KFREE),LFREE)
         NEWGEO = .FALSE.
      END IF
C
C     ... xyz input reader fails/gives wrong results from now on
C         in parallel if there is no update
C         for primitives. Does it change because of reordered geoemtry?
C         FIXME: check enabling reordered geoemtry instead of broadcast of
C         primit.h. - SK 14 Sep 08
      IF( NEWGEO )THEN
#if defined (VAR_MPI)
        CALL GET_PRIMITF()
#endif
      END IF
C
C     Set common-block BLOCKS
C
      CALL PAOVEC(WORK(KFREE),LFREE,0,0)
C

#if defined (VAR_MPI)
      CAll interface_mpi_BCAST(ISEND,2,MPARID,global_communicator)
      NDMAT = ISEND(1)
      NFMAT = ISEND(2)

      IF(NDMAT.GT.0.AND.ITYPE.NE.4) THEN
        LDMAT = NDMAT*N2BASX
        CALL MEMGET('REAL',KDMAT,LDMAT,WORK,KFREE,LFREE)
        CAll interface_mpi_BCAST(WORK(KDMAT),LDMAT,
     &                 MPARID,global_communicator)
      ELSE
        KDMAT = KFREE
      ENDIF

      LFMAT = NFMAT*N2BASX

      IF(NFMAT.GT.0) THEN
        CALL MEMGET('REAL',KFMAT,LFMAT,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIRD,NDMAT,WORK,KFREE,LFREE)
        call interface_mpi_bcast_i1_work_f77(WORK(KIRD),NDMAT,MPARID,
     &               global_communicator)
      ELSE
        KFMAT = KFREE
        KIRD  = KFREE
      ENDIF
      IF(ITYPE.EQ.4) THEN
        KIFC = KFREE
      ELSEIF (ITYPE.EQ.2) THEN
        CALL MEMGET('INTE',KIFC,NDMAT,WORK,KFREE,LFREE)
        call interface_mpi_bcast_i1_work_f77(WORK(KIFC),NDMAT,MPARID,
     &               global_communicator)
      ELSE
        CALL MEMGET('INTE',KIFC,NFMAT,WORK,KFREE,LFREE)
        call interface_mpi_bcast_i1_work_f77(WORK(KIFC),NFMAT,MPARID,
     &               global_communicator)
      END IF
      CAll interface_mpi_BCAST(HFXFAC,1,MPARID,
     &               global_communicator)
      CAll interface_mpi_BCAST(HFXATT,1,MPARID,
     &               global_communicator)
      CAll interface_mpi_BCAST(HFXMU,1,MPARID,
     &               global_communicator)

      call interface_mpi_bcast_l0(spinfr, 1, mparid,
     &               global_communicator)
      KDMRAO = KFREE
      KGABAO = KFREE
      KDMRSO = KFREE
      IF(SCREEN.GT.-1.0D0) THEN
        IF(ITYPE.EQ.3.OR.ITYPE.EQ.9.OR.ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
          CAll interface_mpi_BCAST(N2RED,1,MPARID,
     &                 global_communicator)
          LDMRAO = NDMAT*N2RED
          CALL MEMGET('REAL',KDMRAO,LDMRAO,WORK,KFREE,LFREE)
          CAll interface_mpi_BCAST(WORK(KDMRAO),LDMRAO,
     &                 MPARID,global_communicator)
          IF (ITYPE.EQ.2) THEN
             CALL MEMGET('REAL',KGABAO,4*N2RED,WORK,KFREE,LFREE)
             CAll interface_mpi_BCAST(WORK(KGABAO),4*N2RED,
     &                    MPARID,global_communicator)
          ELSE
             CALL MEMGET('REAL',KGABAO,N2RED,WORK,KFREE,LFREE)
             CAll interface_mpi_BCAST(WORK(KGABAO),N2RED,
     &                    MPARID,global_communicator)
          END IF
        ENDIF
        IF(ITYPE.EQ.9) THEN
          CAll interface_mpi_BCAST(N2RED,1,
     &               MPARID,global_communicator)
          LDMRSO = NDMAT*N2RED
          CALL MEMGET('REAL',KDMRSO,LDMRSO,WORK,KFREE,LFREE)
          CAll interface_mpi_BCAST(WORK(KDMRSO),LDMRSO,
     &               MPARID,global_communicator)
        ENDIF
      ENDIF
#endif
C
      RETURN
      END
#if defined(VAR_MPI)
      SUBROUTINE GET_PRIMITF()
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "primit.h"
#include "infpar.h"
C........primit.h
         CAll interface_mpi_BCAST(PRIEXP,MXPRIM,MPARID,
     &                  global_communicator)
         CAll interface_mpi_BCAST(PRICCF,MXPRIM*MXCONT,MPARID,
     &                  global_communicator)
         CAll interface_mpi_BCAST(PRICRX,MXPRIM,MPARID,
     &                  global_communicator)
         CAll interface_mpi_BCAST(PRICRY,MXPRIM,MPARID,
     &                  global_communicator)
         CAll interface_mpi_BCAST(PRICRZ,MXPRIM,MPARID,
     &                  global_communicator)
      END
#endif
C  /* Deck parlop */
      SUBROUTINE PARLOP(WORK,LWORK,FMAT,DMAT,NDMAT,GMAT,MAXDER,
     &     EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,SOFOCK,
     &     DISTRI,LONDON,
     &     SPNORB,PERTUR,IATOM,MULE,MULTE,NODV,NOPV,NOCONT,
     &     THRESH,JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &     SQ12EL,INDHER,INDHSQ,IODDHR,
     &     IREPDM,IFCTYP,ADISTR,I2TYP,ICEDIF,SCREEN,
     &     GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL,GENCNT,
     &     INDEX_IJ,TIMES)
C
C     Copied from TWOLOP and rewritten by Paal Dahle Nov.1994
C     Minor revisions Oct 1996 T.Saue
C
CMI changes from TEC for LAO
      use interface_to_mpi
#include "implicit.h"
#if defined (VAR_MPI)
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#endif
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "dummy.h"
      PARAMETER (LUPAS = 27, LUPAO = 29)
      PARAMETER(D0=0.0D0)
      LOGICAL PRINTA, PRINTB, PRINTC, PRINTD, NOPV, NODV, PERTUR,
     &        EXPECT, UNDIFF, DDFOCK, DIRFCK, DISTRI, NOCONT, SPNORB,
     &        RETUR, FIRST, SQ12EL, LONDON, SUSCEP, ADISTR, RELCAL,
     &        GENCNT, SOFOCK
      DIMENSION DMAT(*), FMAT(*), GMAT(*), INDHSQ(*), IODDHR(*),
     &          INDHER(*), WORK(LWORK), IREPDM(*), IFCTYP(*),
     &          GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*),
     &          IJSHEL(2),INDEX_IJ(*),TIMES(*)
C
#include "cbisol.h"
#include "twocom.h"
#include "nuclei.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#include "dcbham.h"
#else
#include "energy.h"
#endif
#include "taymol.h"
#include "taysol.h"
#include "suscpt.h"
#include "blocks.h"
#include "symmet.h"
#include "infpar.h"
C
#include "inforb.h"
C
      IF (JPRINT .GT. 5) CALL TITLER('Output from PARLOP','*',103)
C
      FIRST = .TRUE.
      DIRAC = RELCAL
      IF (EXPECT .AND. .NOT.NOPV) THEN
         REWIND LUPAO
      END IF
      IF (SUSCEP) THEN
         IF (.NOT.NOPV) REWIND LUPAO
         IF (.NOT.NOPV) REWIND LUPAS
         CALL DZERO(SUS2EL,9)
      END IF
C
      DOSCRN = .FALSE.
      IF(DIRFCK.OR.SOFOCK.OR.(EXPECT.AND.MAXDER.EQ.1))  THEN
         IF (DIRFCK.OR.SOFOCK) THEN
           CALL DZERO(DINTSKP,8)
         ELSE
           CALL DZERO(DINTSKP,3*14*3)
         END IF
        IF(SCREEN.GT.D0) THEN
          DOSCRN = .TRUE.
          SCRTHR = SCREEN
          ICEFLG = ICEDIF
          NCM = 0
          NEM = 0
          DO I = 1,NDMAT
            IX  = IFCTYP(I) / 10
            IY  = MOD(IFCTYP(I),10)
            IC  = MOD(IY,2)
            IE  = (IY - IC)/2
C           ... no direct Coulomb term for antisymmetric density matrix
            IF (I2TYP.LT.4.AND.IX.EQ.2) IC = 0
C           ... no direct Gaunt term for symmetric density matrix
            IF (I2TYP.EQ.4.AND.IX.EQ.1) IC = 0
            NCM = NCM + IC
            NEM = NEM + IE
          ENDDO
        ENDIF
      ENDIF
      IF(I2TYP.EQ.0) THEN
        ICSTRT = 1
        IDSTRT = 1
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
        ICSTRT = 1
        IDSTRT = 1
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
        ICSTRT = NLRGBL+1
        IDSTRT = NLRGBL+1
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.3) THEN
        ICSTRT = NLRGBL+1
        IDSTRT = NLRGBL+1
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.4) THEN
        ICSTRT = NLRGBL+1
        IDSTRT = 1
        ICSMAX = MAXSHL
        IDSMAX = NLRGBL
      ENDIF
C
C     ***********************************
C     ***** Send request for a task *****
C     ***********************************
C
      ICOUNT = 0
C
#if defined (VAR_MPI)
  100 CONTINUE
C
      CAll interface_mpi_SEND(MYTID,1,MPARID,20,
     &              global_communicator)
C
      CAll interface_mpi_RECV(NTEST,1,MPARID,30,
     &               global_communicator)
C
      IF (NTEST .EQ. -1) GOTO 300
C
      DO 200 I = 1, NTASK
C
         CAll interface_mpi_RECV(IJSHEL,2,MPARID,30,
     &                 global_communicator)
C
C        **************************************
C        ***** Check if it is a null task *****
C        **************************************
C
         IF (IJSHEL(1).EQ. 0) THEN
            ICOUNT = ICOUNT + I - 1
            GO TO 100
         END IF
C
C        The compressed index IJSHEL is split to indices ISHELA and ISHELB.
C
         CALL UNPKIJ(IJSHEL(1),ISHELA,ISHELB)
!Miro: runtimecheck problem - ISHELA gets wrong value, see below !
C
         IF (TIMING) CPUTK1 = SECOND()
C
C        *****************************
C        ***** First Shell Index *****
C        *****************************
C
! Fortran runtime error: Index '-2147483648' of dimension 1 of array 'lclash' below lower bound of 1
! Miro: this means ISHELA get wrong value, needs fix !!!
         ICA    = LCLASH(ISHELA)
         NHKTA  = NHKTSH(ISHELA)
         KHKTA  = KHKTSH(ISHELA)
         KCKTA  = KCKTSH(ISHELA)
         SPHRA  = SPHRSH(ISHELA)
         NCENTA = NCNTSH(ISHELA)
         MULA   = ISTBSH(ISHELA)
         MULTA  = MULT(MULA)
         NSTRA  = IORBSB(IORBSH(ISHELA,1))
         NUCA   = NUCOSH(ISHELA)
         NORBA  = NORBSH(ISHELA)
         IF (.NOT.BIGVEC) THEN
            CORAX0 = CENTSH(ISHELA,1)
            CORAY0 = CENTSH(ISHELA,2)
            CORAZ0 = CENTSH(ISHELA,3)
         END IF
         PRINTA = .TRUE.
         IF ((ISHELA .NE. IPRNTA).AND.(IPRNTA .NE. 0)) PRINTA = .FALSE.
C
C        ******************************
C        ***** Second Shell Index *****
C        ******************************
C
         ICB    = LCLASH(ISHELB)
         NHKTB  = NHKTSH(ISHELB)
         KHKTB  = KHKTSH(ISHELB)
         KCKTB  = KCKTSH(ISHELB)
         SPHRB  = SPHRSH(ISHELB)
         NCENTB = NCNTSH(ISHELB)
         MULB   = ISTBSH(ISHELB)
         MULTB  = MULT(MULB)
         NSTRB  = IORBSB(IORBSH(ISHELB,1))
         NUCB   = NUCOSH(ISHELB)
         NORBB  = NORBSH(ISHELB)
         IF (.NOT.BIGVEC) THEN
            CORBX0 = CENTSH(ISHELB,1)
            CORBY0 = CENTSH(ISHELB,2)
            CORBZ0 = CENTSH(ISHELB,3)
         END IF
         GENAB  = .NOT.(SEGMSH(ISHELA) .AND. SEGMSH(ISHELB))
         IGENAB = 1
         IF (.NOT.GENAB) IGENAB = 2
         NSETA  = NSETSH(ISHELA,IGENAB)
         NSETB  = NSETSH(ISHELB,IGENAB)
         PRINTB = PRINTA
         IF ((ISHELB.NE.IPRNTB).AND.(IPRNTB.NE.0)) PRINTB = .FALSE.
         IF (PRINTB) THEN
            IPRINT = JPRINT
         ELSE
            IPRINT = 0
         END IF
C
C        *****************************
C        ***** Third Shell Index *****
C        *****************************
C
         ICMAX = ISHELA
         IF (SPNORB.OR.I2TYP.EQ.2) ICMAX = MAXSHL
         DO 400 ISHELC = ICSTRT,ICMAX
C
            ICC    = LCLASH(ISHELC)
            NHKTC  = NHKTSH(ISHELC)
            KHKTC  = KHKTSH(ISHELC)
            KCKTC  = KCKTSH(ISHELC)
            SPHRC  = SPHRSH(ISHELC)
            NCENTC = NCNTSH(ISHELC)
            MULC   = ISTBSH(ISHELC)
            MULTC  = MULT(MULC)
            NSTRC  = IORBSB(IORBSH(ISHELC,1))
            NUCC   = NUCOSH(ISHELC)
            NORBC  = NORBSH(ISHELC)
            IF (.NOT.BIGVEC) THEN
               CORCX0 = CENTSH(ISHELC,1)
               CORCY0 = CENTSH(ISHELC,2)
               CORCZ0 = CENTSH(ISHELC,3)
            END IF
            PRINTC = PRINTB
            IF ((ISHELC.NE.IPRNTC).AND.(IPRNTC.NE.0)) PRINTC=.FALSE.
C
C           ******************************
C           ***** Fourth Shell Index *****
C           ******************************
C
            IDMAX = ISHELC
C           Note the order of these tests : can not be changed !
            IF (I2TYP.EQ.4) IDMAX = NLRGBL
            IF (.NOT.SPNORB.AND.(ISHELA.EQ.ISHELC)) IDMAX = ISHELB
            DO 500 ISHELD = IDSTRT,IDMAX
C
               ICD    = LCLASH(ISHELD)
               NHKTD  = NHKTSH(ISHELD)
               KHKTD  = KHKTSH(ISHELD)
               KCKTD  = KCKTSH(ISHELD)
               SPHRD  = SPHRSH(ISHELD)
               NCENTD = NCNTSH(ISHELD)
               MULD   = ISTBSH(ISHELD)
               MULTD  = MULT(MULD)
               NSTRD  = IORBSB(IORBSH(ISHELD,1))
               NUCD   = NUCOSH(ISHELD)
               NORBD  = NORBSH(ISHELD)
               IF (.NOT.BIGVEC) THEN
                  CORDX0 = CENTSH(ISHELD,1)
                  CORDY0 = CENTSH(ISHELD,2)
                  CORDZ0 = CENTSH(ISHELD,3)
               END IF
               GENCD = .NOT.(SEGMSH(ISHELC) .AND. SEGMSH(ISHELD))
               IGENCD = 1
               IF (.NOT.GENCD) IGENCD = 2
               NSETC = NSETSH(ISHELC,IGENCD)
               NSETD = NSETSH(ISHELD,IGENCD)
               PRINTD = PRINTC
               IF ((ISHELD .NE. IPRNTD).AND.(IPRNTD .NE. 0))
     &              PRINTD = .FALSE.
               IF (PRINTD) THEN
                  IPRINT = JPRINT
               ELSE
                  IPRINT = 0
               END IF
C
               SHAEQB = ISHELA .EQ. ISHELB
               SHCEQD = ISHELC .EQ. ISHELD
               SHABAB = (ISHELA.EQ.ISHELC) .AND. (ISHELB.EQ.ISHELD)
C
               IF (SMLV1C) THEN
                  IF (NCENTC.NE.NCENTD.AND.ICC.EQ.2.AND.ICD.EQ.2)
     &               GO TO 500
               ENDIF
C
C              *******************************
C              ***** Calculate integrals *****
C              *******************************
C
               CALL TWOODS(FMAT,DMAT,NDMAT,GMAT,WORK,LWORK,UNDIFF,
     &              PERTUR,LONDON,SPNORB,EXPECT,SUSCEP,DDFOCK,
     &              DIRFCK,SOFOCK,DISTRI,IATOM,MULE,MULTE,MAXDER,
     &              NOCONT,NODV,NOPV,THRESH,IPRINT,FIRST,
     &              SQ12EL,INDHSQ,IODDHR,INDHER,IFCTYP,
     &              ADISTR,DUMMY,ICEDIF,DINTSKP,
     &              GABRAO,DMRAO,DMRSO,IREPDM,IDUMMY,IDUMMY,GENCNT)
C
               IF (RETUR) THEN
                  IF (ISHELA .EQ. IPRNTA .AND.
     &                ISHELB .EQ. IPRNTB .AND.
     &                ISHELC .EQ. IPRNTC .AND.
     &                ISHELD .EQ. IPRNTD) RETURN
               END IF
 500        CONTINUE
 400     CONTINUE
C
         IF (TIMING) THEN
            CPUTK2 = SECOND()
            CPUTSK = CPUTK2 - CPUTK1
            IPLACE = ICOUNT + I
C
            INDEX_IJ(IPLACE) = IJSHEL(2)
            TIMES(IPLACE) = CPUTSK
         END IF
 200  CONTINUE
C
      ICOUNT = ICOUNT + NTASK
      GOTO 100
C
 300  CONTINUE
#endif
c                 print *,'dmat-start'
c                 print *, (dmat(ii),ii=1,10)
c                 print *,'dmat-slutt'
c                 print *, (dmat(ii),ii=n2basx-10,n2basx)ccc
c
c                 print *,'fmat-start'
c                 print *, (fmat(ii),ii=1,10)
c                 print *,'fmat-slutt'
c                 print *, (fmat(ii),ii=n2basx-10,n2basx)c
c
c                 call flshfo(6)
C
C     Symmetrize skeleton Fock matrices for itypes -5 and 8.
C     This reduces the size of data transferred to Master.
C     For itypes 2 and 3 this routine is called outside TWOINT,
C     to avoid doing the same work more than once.
C
      IF (DDFOCK.AND..NOT.LONDON) THEN
         CALL SKLFCK(FMAT,WORK,LWORK,JPRINT,DIRFCK,DDFOCK,EXPECT,PERTUR,
     &               NODV,MAXDER,LONDON,NDMAT,IREPDM,IFCTYP,IATOM)
      END IF
C
C     <<<<< Print Section - Gradient and Hessian Elements >>>>>
C
      IF (EXPECT) THEN
         IF (JPRINT .GT. 0) THEN
            CALL HEADER('Two-electron integral gradient',-1)
            CALL PRIGRD(GRADEE)
            CALL HEADER('Potential energy (NN + NE + EE) gradient',-1)
            CALL ZERGRD
            CALL ADDGRD(GRADNN)
            CALL ADDGRD(GRADNA)
            CALL ADDGRD(GRADEE)
            CALL PRIGRD(GRDMOL)
            CALL HEADER('Molecular gradient',-1)
            CALL ADDGRD(GRADFS)
            CALL ADDGRD(GRADKE)
            IF (SOLVNT) THEN
               CALL ADDGRD(GSOLTT)
               CALL ADDGRD(GSOLNN)
            END IF
            CALL PRIGRD(GRDMOL)
            NCDEP3 = 3*NUCDEP
            GRDNRM = DDOT(NCDEP3,GRDMOL,1,GRDMOL,1)
            GRDNRM = SQRT(GRDNRM)
            WRITE (LUPRI,'(/19X,A,1P,E10.2)')
     *         'Molecular gradient norm:', GRDNRM
            CALL ZERGRD
            IF (MAXDER.EQ.2) THEN
               CALL HEADER('Two-electron integral Hessian',-1)
               CALL PRIHES(HESSEE,'CENTERS')
               CALL HEADER('Potential energy (NN + NE + EE) Hessian',-1)
               CALL ZERHES
               CALL ADDHES(HESSNN)
               CALL ADDHES(HESSNA)
               CALL ADDHES(HESSEE)
               IF (SOLVNT) CALL ADDHES(HSOLT2)
               IF (SOLVNT) CALL ADDHES(HSOLNN)
               CALL PRIHES(HESMOL,'CENTERS')
               CALL ZERHES
            END IF
         END IF
      END IF
      IF (LONDON .AND. MAXDER.EQ.2) THEN
         SUS2EL(2,1) = SUS2EL(1,2)
         SUS2EL(3,1) = SUS2EL(1,3)
         SUS2EL(3,2) = SUS2EL(2,3)
         IF (JPRINT .GT. 1) THEN
            CALL HEADER('Two-electron integral susceptibilities',-1)
            CALL OUTPUT(SUS2EL,1,3,1,3,3,3,1,LUPRI)
         END IF
      END IF
C
C     Print Fock matrices
C
      IF (DIRFCK .AND. JPRINT.GT.3) THEN
         CALL HEADER('Fock matrix in PARLOP',-1)
         DO 600 I = 1, NDMAT
            ISTR = NBAST*NBAST*(I - 1) + 1
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  600    CONTINUE
      END IF
      RETURN
      END
C  /* Deck sdres */
      SUBROUTINE SDRES(FMAT,LFMAT,NUMGRA,MAXDIF,
     &                 NUMHES,NUMSUS,CPU1,WALL1,ITYPE,
     &                 DINTSKP,MTOTTK,INDEX_IJ,TIMES)
C
C           ITYPE  - Calc. type: 2 = expectation values.
C                                3 = direct SCF (AO-basis)
C                                4 = (rs|**)-distributions
C                               -5 = Integrals derivated with respect
C                                    to magnetic field.
C                                7 = direct CC (DALTON)
C                                8 = Fock matrices derivated with respect
C                                    to atomic coordinates.
C                                9 = direct SCF (SO-basis)
C
C     Written by Henrik Koch and Paal Dahle Nov-1993.
C
CMi TECs changes for LAO
C
      use interface_to_mpi
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
#include "suscpt.h"
#include "infpar.h"
#include "inforb.h"
C
      INTEGER   LFMAT, NUMGRA, MAXDIF, NUMHES, NUMSUS, ITYPE
      REAL*8    CPU1, WALL1
      REAL*8    FMAT(LFMAT), DINTSKP(*), TIMES(*)
      INTEGER   INDEX_IJ(*)

#if defined (VAR_MPI)
      REAL*8    RBUF(4)
      REAL*8, allocatable :: my_TIMES(:)

      CALL GETTIM(CPU2,WALL2)
C
C---------------------------
C    use MPI_REDUCE for type 3 ( direct SCF)
C---------------------------
C
      IF (ITYPE .EQ. 3) THEN
         ! using _f77 version
         CALL interface_mpi_reduce_r1_work_f77(
     &        FMAT, DUMMY, LFMAT,
     &        op_MPI_SUM, MPARID, global_communicator)
         CALL interface_mpi_reduce_r1_work_f77(
     &        DINTSKP, DUMMY, 8,
     &        op_MPI_SUM, MPARID, global_communicator)
         IF (TIMING) THEN ! used for better load balancing in next iteration
            allocate (my_TIMES(ICOUNT))
            my_TIMES(1:ICOUNT) = TIMES(1:ICOUNT)
            TIMES(1:MTOTTK) = 0.0D0
            DO I = 1, ICOUNT
               TIMES( INDEX_IJ(I) ) = my_TIMES(I)
            END DO 
            deallocate (my_TIMES)
            CALL interface_mpi_reduce_r1_work_f77(
     &        TIMES, DUMMY, MTOTTK,
     &        op_MPI_SUM, MPARID, global_communicator)
         END IF
         GO TO 6000 ! to send general timing info
      END IF

      CAll interface_mpi_SEND(ICOUNT,1,MPARID,40,
     &                 global_communicator)

      IF (ICOUNT .EQ. 0) RETURN ! this node has done nothing

      IF(ITYPE.EQ.3.OR.ITYPE.EQ.9) THEN
         CAll interface_mpi_send_r1_work_f77(DINTSKP,8,MPARID,40,
     &                 global_communicator)
         CAll interface_mpi_SEND(FMAT,LFMAT,MPARID,40,
     &                 global_communicator)
      ELSE IF (ITYPE.EQ.2) THEN
         CAll interface_mpi_send_r1_work_f77(DINTSKP,3*3*14,MPARID,40,
     &                 global_communicator)
         CAll interface_mpi_SEND(GRADEE,NUMGRA,MPARID,40,
     &                 global_communicator)
      ELSE IF(ITYPE.EQ.-5) THEN
         IF(MAXDIF.EQ.1) THEN
           Call interface_mpi_SEND(FMAT,LFMAT,MPARID,40,
     &          global_communicator)
         ELSEIF(MAXDIF.EQ.2) THEN
           Call interface_mpi_SEND(SUS2EL,9,MPARID,40,
     &          global_communicator)
        ENDIF
      ENDIF
C
C
      IF (TIMING) THEN

         DO I = 1, ICOUNT
            Call interface_mpi_SEND(INDEX_IJ(I),1,MPARID,40,
     &                    global_communicator)
            Call interface_mpi_SEND(TIMES(I),1,MPARID,40,
     &                    global_communicator)
         END DO

      END IF
C
C----------------------------------------------
C     Send general timing results if requested.
C----------------------------------------------
C
 6000 CONTINUE
      IF (IPRPAR .GT. 0) THEN
         CALL GETTIM(CPU3,WALL3)
         RBUF(1) = CPU2 - CPU1
         RBUF(2) = CPU3 - CPU2
         RBUF(3) = WALL2 - WALL1
         RBUF(4) = WALL3 - WALL2
         CAll interface_mpi_SEND(RBUF,4,MPARID,50,
     &                 global_communicator)
C
      END IF
#endif
C
      RETURN
      END
C  /* Deck prlinp */
      SUBROUTINE PRLINP(WORD)
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (NTABLE = 4)
#include "infpar.h"
      LOGICAL FIRST, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      SAVE FIRST
      DATA TABLE /'.NODES ', '.NTASK ', '.PRINT ', 'xXXXXXX'/
      DATA FIRST /.TRUE./
C
C-------------------------
C     Initialize /PRLINP/.
C-------------------------
C
      IF (.NOT. FIRST) THEN
         IF ( WORD .EQ. '*PARALL' ) THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      FIRST = .FALSE.
C
C     For infpar.h
C
! SLAVE is set and PARHER is initialized to .false. in dirac.F at beginning of execution
!     SLAVE  = .FALSE.
!     PARHER = .FALSE.
#if defined (VAR_MPI)
      NODES = NUMNOD
#else
      NODES  = 0
#endif
      NTASK  = 1
      IPRPAR = 0
      TOTWAL = 0.0D0
C
      NEWDEF = (WORD .EQ. '*PARALL')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in PRLINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in PRLINP')
    1          CONTINUE
                  READ(LUCMD,*) NODES
               GO TO 100
    2          CONTINUE
                  READ(LUCMD,*) NTASK
               GO TO 100
    3          CONTINUE
                  READ(LUCMD,*) IPRPAR
               GO TO 100
    4          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in PRLINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in PRLINP')
            END IF
      END IF
  300 CONTINUE
C
#if defined (VAR_MPI)
      IF (NODES .NE. NUMNOD) THEN
         WRITE (LUPRI,'(/A,I3,A,I3)') 'Number of nodes in input ',
     &        NODES,' .ne. nodes in MPI calculation ',NUMNOD
         CALL QUIT('Wrong number of nodes in input')
      END IF
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Parsr2 */
      SUBROUTINE PARSR2(NDIM,NPOS,TIMES,INDEX_IJ,IPRINT)
C*****************************************************************************
C
C     Sort IJ-indices according to increasing run time...
C     Written by T.Saue Oct 3 1996
C     Rewritten 3-Aug-2017 by hjaaj because old code sorted after decreasing run time
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(THRES = 1.0D-2)
      DIMENSION NPOS(NDIM),TIMES(NDIM),INDEX_IJ(NDIM)
      INTEGER, ALLOCATABLE :: IPLACE(:)


!     Sort IPLACE pointer array
!     such that TIMES( IPLACE(I) ) .le. TIMES( IPLACE(I+1) )
      allocate(IPLACE(NDIM))
      DO I = 1, NDIM
         IPLACE(I) = I
         TI = TIMES( I )
         DO J = 1,(I-1)
            TJ = TIMES( IPLACE(J) )
            IF (TJ .GT. TI) THEN
               DO K = I, (J+1), -1
                  IPLACE(K) = IPLACE(K-1)
               END DO ! K
               IPLACE(J) = I
               EXIT ! exit loop, found the right place to put it in
            END IF
         END DO ! J
      END DO ! I

      IF (IPRINT .GE. 6) THEN
         DO I = 1, NDIM
           WRITE(LUPRI,'(15X,A,I8,A,F8.3)')
     &           'IJ-shell =', INDEX_IJ(IPLACE(I)),
     &           '    time =', TIMES(IPLACE(I))
        ENDDO
      ENDIF

C
C     Transfer sorted order of INDEX_IJ to NPOS
C
      DO I = 1,NDIM
         NPOS(I) = INDEX_IJ(IPLACE(I))
      ENDDO
C
      deallocate(IPLACE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck HERPAR_det_task */
      SUBROUTINE HERPAR_DET_TASK (I2TYP,IMODE,MTOTTK,NPOS,TIMES)
C*****************************************************************************
C
C     Determine tasks that are to be distributed.
C     IMODE 1 : Determine only total number
C     IMODE 2 : Fill the array NPOS
C
C     Written by Luuk Visscher, august 2001
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "nuclei.h"
#ifdef PRG_DIRAC
#include "dcbham.h"
#endif
#include "blocks.h"
C
      INTEGER  NPOS(*)
      REAL*8   TIMES(*)
C
      IF(I2TYP.EQ.0) THEN
        IASTRT = 1
        IBSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
        IASTRT = 1
        IBSTRT = 1
        IASMAX = NLRGBL
        IBSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
        IASTRT = 1
        IBSTRT = 1
        IASMAX = NLRGBL
        IBSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.3) THEN
        IASTRT = NLRGBL+1
        IBSTRT = NLRGBL+1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.4) THEN
        IASTRT = NLRGBL+1
        IBSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = NLRGBL
      ELSE
        WRITE(LUPRI,'(A,I5)') 'DET_TOTTK: Unknown I2TYP =' ,I2TYP
        CALL QUIT('Unknown I2TYP !!!')
      ENDIF
C
C     *****************************
C     ***** First Shell Index *****
C     *****************************
C
      IND = 0
      DO 100 ISHELA = IASTRT,IASMAX
C
         ICA    = LCLASH(ISHELA)
         NCENTA = NCNTSH(ISHELA)
         NORBA  = NORBSH(ISHELA) * KHKTSH(ISHELA)
C
C        ******************************
C        ***** Second Shell Index *****
C        ******************************
C
         IBMAX = ISHELA
         IF (I2TYP.EQ.4) IBMAX = NLRGBL
         DO 200 ISHELB = IBSTRT,IBMAX
C
            ICB    = LCLASH(ISHELB)
            NCENTB = NCNTSH(ISHELB)
            NORBB  = NORBSH(ISHELB) * KHKTSH(ISHELB)
C
C           Skip two-center SS contributions if requested
C           This test is provisory, in SYMLOP we take out more integrals
C           in case of symmetry-degenerate nuclei. Testing at this place,
C           will, however, take out the work in ODCVEC and EXCOEF that
C           otherwise would cause a quartic scaling term.
C
            IF (SMLV1C) THEN
               IF (NCENTA.NE.NCENTB.AND.ICA.EQ.2.AND.ICB.EQ.2) GO TO 200
            ENDIF
C
            IJ  = (ISHELA*(ISHELA-1))/2 + ISHELB
            IF (IMODE.EQ.2) THEN
               NPOS(MTOTTK-IND) = IJ
               ! hjaaj Aug 2017: assume time proportional to number of (ij/ functions for this (IJ/ shell combination
               IF (ISHELB .EQ. ISHELA) THEN
                  TIMES(MTOTTK-IND) = (NORBA*(NORBA+1))/2
               ELSE
                  TIMES(MTOTTK-IND) = NORBA*NORBB
               END IF
            END IF
            IND = IND + 1
C
  200    CONTINUE
  100 CONTINUE
C
      IF (IMODE.EQ.1) MTOTTK = IND
C
      RETURN
      END
