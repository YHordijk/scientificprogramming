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

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PAMCCM ()
C
C---------------Description-----------------------------------------------
C
C    Double group symmetry adapted unrestricted coupled cluster code
C    This version is called from within DIRAC
C
C---------------Routines called-------------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C              minor modifications for use of RELCCSD-MP2 for 
C              MP2-NO generation module. S. Knecht - June 2008
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      use interface_to_mpi

#include "implicit.h"
#include "maxorb.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbpsi.h"
#include "dcbmp2no.h"
#include "infpar.h"
#include "ccpar.inc"
#include "files.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 SEC,SEC0
      LOGICAL CA
      CHARACTER*10 INTFCE
C
C---------------Executable code--------------------------------------
C
      CALL CPUUSED(SEC0)
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      IN = LUCMD
      IW = LUPRI
C aspg 20080903
C iface set to DIRAC6 as default, in line with having scheme 6 as default in moltra
c     INTFCE = 'DIRAC     '
      INTFCE = 'DIRAC6    '
      IF( IAOMOSCHEME .eq. 4 ) INTFCE = 'DIRAC     '
      MP2_DENRUN = .FALSE.
      IF( DOMP2NO )THEN
        IF( RUN_CCMOD ) MP2_DENRUN = .TRUE.
      END IF
      CA = NZ.GE.2.AND..NOT.(SPINFR.OR.LEVYLE)
CMI   ... if CAP is on, switch to complex algebra
      IF (CAP) CA = .TRUE.
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
C     Summon the slaves, who are waiting in the general menu routine.
C
      CALL CCMSTI
C
C     Tell them the correct arithmetic (NZ is not initialized on the slaves
C     if you do not run HF !)
C
#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         call interface_mpi_BCAST_l0(CA,1,MASTER,global_communicator)
         call interface_mpi_BCAST_l0(MP2_DENRUN,1,MASTER,
     &                  global_communicator)
      END IF
C
C     if MP2-NO is enabled, tell the slaves about the correct interface
C
      IF( MP2_DENRUN )THEN
        call interface_mpi_BCAST(IAOMOSCHEME,1,MASTER,
     &                 global_communicator)
        IF( IAOMOSCHEME .eq. 4 ) INTFCE = 'DIRAC     '
      END IF
#endif
C
C     Enter the generic CCMAIN routine (master becomes now primus inter parens)
C
      CALL CCMAIN(IN,IW,CA,INTFCE)
C
C     Return the slaves to the general menu routine.
C
      CALL CCMSTE
C
      CALL CPUUSED(SEC)
      WRITE (IW,1000) SEC-SEC0,SEC0,SEC
      WRITE (IW,'(//80A1)') ('#',I=1,80)
      CALL FLSHFO(IW)
      CLOSE(LUCMD)
C
 1000 FORMAT (//' CPU time (seconds) used in RELCCSD:',T50,F14.4/
     & ' CPU time (seconds) used before RELCCSD:',T50,F14.4/
     & ' CPU time (seconds) used in total sofar:',T50,F14.4//
     & '  --- Normal end of RELCCSD Run ---')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CCMSTI */
      SUBROUTINE CCMSTI
C
C     Written by Luuk Visscher, August 1997
C     Call the slaves if we run in parallel
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "maxorb.h"
#include "dcbgen.h"
C
C     Send task code for parallel RELCCSD
C
      IF (PARCAL) CALL DIRAC_PARCTL( RELCC_PAR ) 
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CCMSTE */
      SUBROUTINE CCMSTE
C
C     Written by Luuk Visscher, August 1997
C     Release the slaves if we run in parallel
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use interface_to_mpi
#include "implicit.h"
#include "infpar.h"
#ifdef VAR_MPI
      integer   istat(df_mpi_status_size)
#endif
C
      if(numnod.eq.0) return
C
C     Find the slaves and release them
C
      NTEST = -1
#if defined (VAR_MPI)
      DO I = 1,NUMNOD
         CALL interface_MPI_IRECV(NODE,1,df_MPI_ANY_SOURCE,20,
     &                            global_communicator,IREQ)
         CALL interface_MPI_WAIT(IREQ,ISTAT)
         CALL interface_MPI_SEND(NTEST,1,NODE,30,
     &                           global_communicator)
      ENDDO
#endif
C
C     Let the slaves sleep (if relevant for this operating system)
C
      CALL DIRAC_PARCTL( RELEASE_NODES )
      RETURN
      END
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CCNOD */
      SUBROUTINE CCNOD
C
C     Written by Luuk Visscher, August 1997
C     The slave RELCCSD routine (DIRAC interface)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use interface_to_mpi

#include "implicit.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dcbmp2no.h"
#include "infpar.h"
#if defined (VAR_MPI)
      integer   istat(df_mpi_status_size)
#endif
#include "priunit.h"
#include "dgroup.h"
#include "ccpar.inc"
#include "files.inc"
      CHARACTER*10 INTFCE, MP2NOF
      LOGICAL CA,CCEND
C
C     arrange for the MPI stuff and correct node number
C     to the total number of running invocations.
C
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

C aspg 20080903
C iface set to DIRAC6 as default, in line with having scheme 6 as default in moltra
c     INTFCE = 'DIRAC     '
      INTFCE = 'DIRAC6    '
      IF( IAOMOSCHEME .eq. 4 ) INTFCE = 'DIRAC     '
      MP2_DENRUN = .FALSE.
C
C     Get the correct arithmetic (NZ is not initialized on slaves if you do not run HF !)
C
#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         call interface_mpi_BCAST_l0(CA,1,MASTER,global_communicator)
         call interface_mpi_BCAST_l0(MP2_DENRUN,1,MASTER,
     &                  global_communicator)
      END IF
      IF( MP2_DENRUN )THEN
        call interface_mpi_BCAST(IAOMOSCHEME,1,MASTER,
     &                  global_communicator)
        IF( IAOMOSCHEME .eq. 4 ) INTFCE = 'DIRAC     '
      END IF
#endif
C
C     create a node-unique filename as output file. Important on
C     shared file systems. Otherwise all the output gets mingled in one
C     file
C
      CCBASF="RELCCSD.OUT"
      
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (CCFILN,'(A11,A1,I1)') CCBASF,'.',MYPROC
         LCCFIL    = 13
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (CCFILN,'(A11,A1,I2)') CCBASF,'.',MYPROC
         LCCFIL    = 14
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (CCFILN,'(A11,A1,I3)') CCBASF,'.',MYPROC
         LCCFIL    = 15
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND CCNOD.F MODULE")
      ENDIF
#else
      CCFILN    = CCBASF
      LCCFIL    = 11
#endif
C
C
C     open the local input file and the node specific output file.
C     The associated file handle MSLVOUT is transferred to CCMAIN.
C     Every access to the local stdout handle then automatically writes
C     to the corresponding output file.
C
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      OPEN(MSLVOUT,FILE = CCFILN(1:LCCFIL))
      CALL CCMAIN(LUCMD,MSLVOUT,CA,INTFCE)
      CLOSE (MSLVOUT,STATUS='KEEP')
      CLOSE (LUCMD,STATUS='KEEP')
C
C     The slave is finished, wait for the master
C
      CCEND = .FALSE.
#if defined (VAR_MPI)
      CALL interface_MPI_ISEND(MYTID,1,MPARID,20,
     &                         global_communicator,IREQ)
C
    1    CONTINUE
         CALL SLEEP(1)
         call interface_mpi_iprobe(MPARID,30,CCEND,global_communicator,
     &                             istat)
         IF (.NOT.CCEND) GOTO 1
C
      CALL interface_MPI_RECV(NTEST,1,MPARID,30,global_communicator)
#endif
C
      RETURN
      END
#endif    /* ifdef VAR_MPI */
