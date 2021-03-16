!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!      *****************************************************
!      ****     POLARIZATION  PROPAGATOR  MODULE        ****
!      ****                   for                       ****
!      ****           ELECTRONIC EXCITATIONS            ****
!      *****************************************************
!
!                       AUTHOR: M.Pernpointner
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! define task symbols for CALL DIRAC_PARCTL( task )
#include "dirac_partask.h"

      SUBROUTINE PAMPOLPRP ()
!
!---------------Description-----------------------------------------------
!
!    Polarization propagator code for the calculation of 
!    electronic excitations.
!    In parallel execution this is the part called by the dirac(...)
!    subroutine, the MASTER PART. 
!    It has to wake up the slaves with the correct ITASK
!    number in order to run the shared calculation.
!    Afterwards it gives itself to the shared calculation.
!
!---------------Common Blocks--------------------------------------
!
      use interface_to_mpi
#include "implicit.h"
#include "maxorb.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbpsi.h"
#include "dcbgen.h"
#include "dcbmp2no.h"
#include "../moltra/dcbtra.h"
#include "infpar.h"
#include "../relccsd/ccpar.inc"
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8 SEC,SEC0
      LOGICAL CA
      CHARACTER*10 INTFCE
#ifdef VAR_MPI
      integer istat(df_mpi_status_size)
#endif
!
!---------------Executable code--------------------------------------
!
      CALL CPUUSED(SEC0)
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      IN = LUCMD
      IW = LUPRI
!
!
!  determine condition when we have complex arithmetic !
!
      CA = NZ.GE.2.AND..NOT.(SPINFR.OR.LEVYLE)
!MI   ... if CAP is on, switch to complex algebra
      IF (CAP) CA = .TRUE.
#if defined (VAR_MPI)
      MASTER = MPARID
      MYPROC = MYTID
!     Add the master node
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
!
!  Do a "wake up" with the correct ITASK in order to guide the slaves
!  to the corresponding POLPRP section.
!
      IF (PARCAL) CALL DIRAC_PARCTL( POLPRP_PAR )

! iface for integral reading taken from common block dcbtra.h
      IF ( ISTRAT .eq. 4 ) THEN
          INTFCE = 'DIRAC     '
      ELSEIF ( ISTRAT .eq. 6 ) THEN
          INTFCE = 'DIRAC6    '
      ELSE
          CALL QUIT ('Unknown DIRAC MO integrals interface')
      END IF

#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         call interface_mpi_BCAST_l0(CA,1,MASTER,global_communicator)
         call interface_mpi_BCAST_l0(MP2_DENRUN,1,MASTER,
     &                  global_communicator)
        call interface_mpi_BCAST(ISTRAT,1,MASTER,
     &                 global_communicator)
      END IF
!
!     if MP2-NO is enabled, tell the slaves about the correct interface
!
      IF( MP2_DENRUN )THEN
        call interface_mpi_BCAST(IAOMOSCHEME,1,MASTER,
     &                 global_communicator)
        call interface_mpi_BCAST(INFO_STRING(1),62,MASTER,
     &                 global_communicator)
        call interface_mpi_BCAST(INFO_STRING(2),62,MASTER,
     &                 global_communicator)
        call interface_mpi_BCAST(INFO_STRING(3),62,MASTER,
     &                 global_communicator)
      END IF
#endif
      IF (ISTRAT /= IAOMOSCHEME) THEN
         print*,"WARNING: fixed incompatible schemes for MP2 and MOLTRA"
         print*,' Dirac MO integrals interface used: ',INTFCE
         IAOMOSCHEME = ISTRAT
      END IF
!
!     Enter the generic POLPRP_MAIN routine (master becomes now primus inter pares)
!
      CALL POLPRP_MAIN(IN,IW,CA,INTFCE)
!
!     Return the slaves to the general menu routine.
!
!
      if(numnod.eq.0) return
!
!     Find the slaves and release them
!
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
!
!     Let the slaves sleep (if relevant for this operating system)
!
      CALL DIRAC_PARCTL( RELEASE_NODES )

      CALL CPUUSED(SEC)
      WRITE (IW,1000) SEC-SEC0,SEC0,SEC
      WRITE (IW,'(//80A1)') ('#',I=1,80)
      CLOSE(LUCMD)
!
 1000 FORMAT (//' CPU time (seconds) used in POLPRP:',T50,F14.4/
     & ' CPU time (seconds) used before POLPRP:',T50,F14.4/
     & ' CPU time (seconds) used in total sofar:',T50,F14.4//
     & '  --- Normal end of POLPRP Run ---')
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE POLPRP_NODE
!
!  similar to CC_NOD due to identical requirements
!
      use interface_to_mpi

#include "implicit.h"
#include "maxorb.h"
#include "dcbham.h"
#include "../moltra/dcbtra.h"
#include "dcbmp2no.h"
#include "infpar.h"
#if defined (VAR_MPI)
      integer istat(df_mpi_status_size)
#endif
#include "priunit.h"
#include "dgroup.h"
#include "../relccsd/ccpar.inc"
#include "../relccsd/files.inc"
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

C
C     Get the correct arithmetic (NZ is not initialized on slaves if you do not run HF !)
C
#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         call interface_mpi_BCAST_l0(CA,1,MASTER,global_communicator)
         call interface_mpi_BCAST_l0(MP2_DENRUN,1,MASTER,
     &                  global_communicator)
        call interface_mpi_BCAST(ISTRAT,1,MASTER,
     &                 global_communicator)
      END IF
      IF( MP2_DENRUN )THEN
        call interface_mpi_BCAST(IAOMOSCHEME,1,MASTER,
     &                  global_communicator)
      END IF
#endif
! iface for integral reading is now communicated from master, set the string as well
      IF ( ISTRAT .eq. 4 ) THEN
          INTFCE = 'DIRAC     '
      ELSEIF ( ISTRAT .eq. 6 ) THEN
          INTFCE = 'DIRAC6    '
      ELSE
          CALL QUIT ('Unknown DIRAC MO integrals interface')
      END IF
C
C     create a node-unique filename as output file. Important on
C     shared file systems. Otherwise all the output gets mingled in one
C     file
C
      CCBASF="POL_PRP.OUT"
      
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
         CALL QUIT("NMPROC.GT.1000! EXTEND POLPRP_NOD Routine!")
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

      CALL POLPRP_MAIN(LUCMD,MSLVOUT,CA,INTFCE)

      CLOSE (MSLVOUT,STATUS='KEEP')
      CLOSE (LUCMD,STATUS='KEEP')
      write(*,*) 'Node',MYPROC,' has arrived after POLPRP_MAIN!'
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
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE POLPRP_MAIN (IN123,IW123,CA,INTFCE)

      use memory_allocator
      use interface_to_mpi
      use polprp_cfg

      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Driver for Polarization Propagator 
!
!---------------Calling variables--------------------------------------
!
      INTEGER IN123,IW123
      LOGICAL CA
      CHARACTER*10 INTFCE
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/files.inc"
#include "../relccsd/inpt.inc"
#include "../relccsd/symm.inc"
#include "../relccsd/eqns.inc"
#include "../relccsd/complex.inc"
#include "../relccsd/ccpar.inc"
#include "../relccsd/waio.h"
#include "dcbham.h"
!
!---------------Local variables--------------------------------------
!
      integer , parameter        :: MAXSP=1000
      INTEGER IREPSP(MAXSP*4),IREPSPI(MAXSP*MXREP*2)
      REAL*8 EPS(MAXSP),EPS1(MAXSP)
      REAL*8 CPU0,CPU1,CPU2,CPU3

      real(8), allocatable :: CC(:)

      integer :: ncore_elec_adc = 0  ! new argument for RDSYMI, needed in RELCCSD, not in ADC.

      integer           :: maxcor,mxcorr,mxcorr1,mword_limit
      real*8            :: total_mem_new = 0   !number of bytes, if number of megawords is given in --aw
      real*8            :: new_mword_limit  = 0

      integer     :: MAXCORE    ! contains --mw * 1000 * 1000 (number of words, decimal mega)
      integer     :: maxcore_info(2)
      character*6 :: maxcore_infc(2)

      logical     :: skip_ccseti
!
!---------------interface region --------------------------------------
!

      interface

        SUBROUTINE POLPRP_PAREXE(i1,i2,ra1,ia1,i3,r1)
          integer                               :: i1,i2,i3
          real*8, dimension(:)                  :: ra1
          integer,dimension(:)                  :: ia1
          real*8                                :: r1
        END SUBROUTINE

      end interface
!
!---------------executable code --------------------------------------
!
      maxcore_info = 0
      maxcore_infc = (/'CCSETI','POLPRP'/)

#if defined (VAR_MPI)
      IF(MYPROC.eq.MASTER) then
        skip_ccseti = polprp_skipccseti
        write(*,*) 'Master will communicate skip flag.'
      ENDIF
#else
      skip_ccseti = polprp_skipccseti
#endif

      call legacy_lwork_get(MAXCORE)

!
!---------------Executable code--------------------------------------

!
!     The INTERFACE/EQNS variable is defined in ccinpt.inc and is needed
!     by CCSETI. The default scheme is DIRAC6. This is also hardcoded
!     in the RELCCSD code.
!
      CARITH = CA
      IN=IN123
      IW=IW123
      INTERFACE = INTFCE
      EQNS='RELADC'
!
!    indicate next module run
!
      write(IW,1999)
 1999 format(//4X,90('*')/
     &   4X,21('*'),5X,'Start of Polarization Propagator Module',4X,
     &      21('*')/
     &   4X,90('*')//)
      WRITE(IW,*)
      WRITE(IW,*) ' The theory for this module is outlined in:'
      WRITE(IW,*)
      WRITE(IW,*) ' M. Pernpointner, J. Chem. Phys. 140, 084108 (2014).'
      WRITE(IW,*) ' doi: 10.1063/1.4865964'
      WRITE(IW,*)
      WRITE(IW,*) ' The parallelization and transition moments'
      WRITE(IW,*) ' are described in'
      WRITE(IW,*)
      WRITE(IW,*) ' M. Pernpointner, L. Visscher and A. B. Trofimov,'
      WRITE(IW,*) ' J. Chem. Theory Comput. 14, 1510 (2017).'
      WRITE(IW,*) ' doi: 10.1021/acs.jctc.7b01056'
      WRITE(IW,*)
      WRITE(IW,*)
      WRITE(IW,*)' ________________________________________________ '
      WRITE(IW,*)'|                                                |'
      WRITE(IW,*)'|       Integral preparation section             |'
      WRITE(IW,*)'|________________________________________________|'
      WRITE(IW,*)

      WRITE (IW,1011) MYPROC+1,NMPROC
 1011 FORMAT (/'---< Process ',I5,' of ',I5,'----<'/)
!
!     Initialize low-level I/O
!
      CALL WAIO_INIT (IW)
!
!     Initialize the timing.
!
      CALL CPUUSED(CPU0)
!
!     Write out memory information and
!     call new mem alloc routines
!
      WRITE (IW,*) 'Memory info:'
      WRITE (IW,*) 'Megawords (MB) given by --mw: ',
     &             MAXCORE/(1000*1000),MAXCORE / (128 * 1024)
      call allocator_get_available_mem(total_mem_new)
      WRITE (IW,*)
      WRITE (IW,*) 'Maximally allowed dynamical memory (MW,B):',
     &      total_mem_new/(8*1000*1000),total_mem_new
      maxcor = nint((total_mem_new/8.0d0),kind(8))
      mword_limit = maxcor / (1024*1024)
!
!     Read symmetry information from MRCONEE
!
      CALL RDSYMI (NSP,EPS1,IREPSP,ncore_elec_adc)
      IF (NSP.GT.MAXSP) THEN
         PRINT*,' INCREASE MAXSP IN POLPRP_MAIN TO',NSP
         CALL QUIT('NSP.GT.MAXSP')
      ENDIF
!
      IPRNT = 2   ! print orbital energies (variable in inpt.inc!)
!
!     Number of integer words for REAL*8 and 
!     number of REAL*8 words per variable (2 for complex groups, 
!     1 for real)
!
      IRW = INTOWP(1)
      IF (CARITH) THEN
         RCW = 2
      ELSE
         RCW = 1
      ENDIF
!
!     Set up symmetry tables and parallel distribution arrays
!
      CALL SYMTAB (NELEC,NFROZ,NINACT,NACT,NSP,EPS1,EPS,IREPSP,
     &             IREPSPI,.FALSE.,.TRUE.,NELEC_F1,NELEC_F2)

      write(iw,*) 'NSP  :',NSP
      write(iw,*) 'MAXSP:',maxsp
      write(iw,*) 'MXREP:',mxrep
      write(iw,*) 'Length of irepspi array:',MAXSP*MXREP*2,' integers'
!
!  determine number of electrons and print general information
!
      NELECT = 0
      DO I = 1, 16
         NELECT = NELECT + NELEC(I)
      ENDDO
!
      WRITE (IW,1004) NELECT,NSP-NELECT,CARITH
1004  FORMAT (//' Number of electrons :',T40,I5
     &        /' Number of virtual spinors :',T40,I5
     &        /' Complex arithmetic mode :',T44,L1/)

!
!  initialize USEOE=T default also for POLPRP
      USEOE   = .true.
!  because it relies on CCSETI and needs the orbital
!  energies from the original HF calculation
!  stefan: i am not sure whether you always want to set USEOE = .true.
!          but let us keep it.
!    
!  in principle this allows for the use of X2Cmmf
!  but the NORECMP parameter is NOT irrelevant - stefan march 2014
!
!
      if(X2CMMF)then
!       X2Cmmf hamiltonian in use - by default always set USEOE and NORECMP to .true.
!       everything else results in nonsense energies
        NORECMP = .true.
        useoe   = .true. 
        write(IW,'(/A,/A,/A,/A,/A/)')
     &  ' ************************************************************',
     &  ' *** molecular-mean-field X2C Hamiltonan (X2Cmmf) active  ***',
     &  ' *** use of SCF orbital energies and no re-computation of ***',
     &  ' *** Fock matrix elements set explicitly                  ***',
     &  ' ************************************************************'
      else
        if(USEOE) WRITE (IW,*) 'Orbital energies from HF are used'
      end if

      ICALC = 2
!
! we offer to skip ccseti completely. This is provided because for really
! large files this step can take unexpectedly large amount of time.
! If the user provides the ftxx.yy files on the shared scratch file
! system, it should work. Attention: the files are opened in CCSETI on the
! master and slave. Therefore we have to open them explicitly, if CCSETI
! is to be skipped.
!

! Input variables only set for master. In case of parallel run we have
! to communicate the flag.

#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         call interface_mpi_BCAST_l0(skip_ccseti,1,MASTER,
     &                               global_communicator)
      END IF
!     print*,'+++++ node',MYPROC,' has skip flag:',skip_ccseti
#endif
      
      if(.not.skip_ccseti) then
        write(iw,*) 'Entering CCSETI'
        CALL CCSETI (ICALC,maxcore_info(1),EPS,NSP,IREPSP,IREPSPI)
      else
        write(iw,*)
        write(iw,*) 'Skipping CCSETI, opening files explicitly.'
        write(iw,*)
        write(iw,*) '*** All the ftXX.YY files have to be on the'
        write(iw,*) '*** shared scratch file system. This option'
        write(iw,*) '*** will not work on strictly distributed'
        write(iw,*) '*** file systems.'
        write(iw,*)
        write(iw,*)
        CALL WAIO_OPEN(ITAPT+0)
        CALL WAIO_OPEN(ITAPT+1)
        CALL WAIO_OPEN(ITAPT+2)
        CALL WAIO_OPEN(ITAPT+3)
        CALL WAIO_OPEN(ITAPT+4)
        CALL WAIO_OPEN(ITAPT+5)
        CALL WAIO_OPEN(ITAPT+6)
      endif

!     WRITE(IW,*) 'Memory requirements of ',
!    &             maxcore_infc(1),maxcore_info(1)
!
!  call *parallel* polarization propagator code
!  
      CALL POLPRP_PAREXE(IW,ITAPADC,EPS,IREPSPI,NSP,total_mem_new)
!
!     return after closing low-level I/O files
!
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+0,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+1,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+2,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+3,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+4,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+5,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+6,JCODE)
!
!     Print I/O statistics
!
      CALL WAIO_PRST (IW)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE POLPRP_PAREXE(IW,IOBASE,EPS,IREPSPI,NSP,total_mem_new)
!
! executed in parallel or in serial from here.
!
      use interface_to_mpi
      use polprp_cfg     !these data are visible on master node only!
                         !They contain control input for the POLPRP run
                         !and user input for the Davidson diagonalizer
      use qstack         !make stack system available to master node

 
      IMPLICIT NONE
!
!---------------Argument Variables---------------------------------
!
      integer                           :: iw       ! I/O file handle for screen output
      integer                           :: iobase   ! I/O base file handle, all purposes
      real*8, dimension(:)              :: eps
      integer, dimension(:)             :: irepspi  ! spinor definitions
      integer                           :: nsp      ! number of active spinors
      real*8                            :: total_mem_new  ! current max. dyn. memory in bytes!
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
!LTR/LFA are defined here (true, false)
#include  "../relccsd/ccpar.inc"
#include "polprp_servercodes.h"
#include "polprp_stacklines.h"
!
!---------------Local variables---------------------------------
!
      integer                 :: n1,n2
      integer                 :: ialloc
      real*8, allocatable     :: daibj(:)
      integer                 :: irploop,desrep
      CHARACTER*6             :: namematfil = 'PP_ELE'  !name of the main ADC matrix file
      CHARACTER*6             :: namecnffil = 'XCONFG'  !name of the configuration file
      integer                 :: irecl, ireccnt, intbuf
      integer                 :: ladc
      integer                 :: lencouplbl = 0
      integer                 :: nbufs = 0
      integer                 :: signal
      integer                 :: serialindicator
      integer                 :: i
      logical                 :: propertyints
      integer                 :: qst_numofstacks
      integer                 :: qst_verbose
      logical                 :: isthere
!  
!  variables for rdprop_adc integral testing
!
      integer      ::  iprnt = 1
      character*8  :: nameofprop
      real*8, allocatable   :: a_CBUF(:)
      real*8, allocatable   :: a_EVO(:), a_EOV(:), a_EOO(:), a_EVV(:)
      complex*16            :: ephase
!
!-----------------------interface section---------------------
!
      interface

        SUBROUTINE print_input_info(i1,i2)
          integer                  :: i1,i2
        END SUBROUTINE

        SUBROUTINE property_handler(i1,i2,l1,i3,ia1)
          integer                  :: i1,i2,i3
          logical                  :: l1
          integer, dimension(:)    :: ia1
        END SUBROUTINE

        SUBROUTINE propinit(i1,i2,l1,i3,ia1)
          integer                  :: i1,i2,i3
          logical                  :: l1
          integer, dimension(:)    :: ia1
        END SUBROUTINE

        SUBROUTINE BUILD_X_MAIN (ra1,ra2,i1)
          real*8, dimension(:)  :: ra1
          real*8, dimension(:)  :: ra2
          integer               :: i1
        END SUBROUTINE

        SUBROUTINE SOLVE_PHPH_R(ra1,i1,i2)
          real*8, dimension(:)  :: ra1
          integer               :: i1,i2
        END SUBROUTINE

        SUBROUTINE SOLVE_PHPH_C(ra1,i1,i2)
          real*8, dimension(:)  :: ra1
          integer               :: i1,i2
        END SUBROUTINE

        SUBROUTINE BUILD_X_COUSAT_R(h,a,i1,i2,i3,i4,l1,r2,i5)
          real*8, dimension(:)  :: h
          real*8, dimension(:)  :: a
          integer               :: i1,i2,i3,i4
          logical               :: l1
          real*8                :: r2
          integer               :: i5
        END SUBROUTINE

        SUBROUTINE BUILD_X_COUSAT_C(h,a,i1,i2,i3,i4,l1,r2,i5)
          real*8, dimension(:)  :: h
          real*8, dimension(:)  :: a
          integer               :: i1,i2,i3,i4
          logical               :: l1
          real*8                :: r2
          integer               :: i5
        END SUBROUTINE

        SUBROUTINE XKENEDIAG(i1,ra1)
          integer               :: i1
          real*8, dimension(:)  :: ra1
        END SUBROUTINE

        SUBROUTINE mpi_bxcs_r(ra1,i1,ra2,i2,i3,i4,i5,l1,r1,i6,i7)
          real*8, dimension(:)  :: ra1
          integer               :: i1,i2,i3,i4,i5,i6,i7
          real*8, dimension(:)  :: ra2
          real*8                :: r1
          logical               :: l1
        END SUBROUTINE

        SUBROUTINE mpi_bxcs_c(ra1,i1,ra2,i2,i3,i4,i5,l1,r1,i6,i7)
          real*8, dimension(:)  :: ra1
          integer               :: i1,i2,i3,i4,i5,i6,i7
          real*8, dimension(:)  :: ra2
          real*8                :: r1
          logical               :: l1
        END SUBROUTINE

        SUBROUTINE mpi_davidson_master_r(i1,i2,i3,i4,c1,c2,i5,i6,i7)
          integer               :: i1,i2,i3,i4,i5,i6,i7
          character*6           :: c1,c2
        END SUBROUTINE

        SUBROUTINE mpi_davidson_master_ext_r(i1,i2,i3,i4,c1,c2,i5,i6,i7)
          integer               :: i1,i2,i3,i4,i5,i6,i7
          character*6           :: c1,c2
        END SUBROUTINE

        SUBROUTINE mpi_davidson_master_c(i1,i2,i3,i4,c1,c2,i5,i6,i7)
          integer               :: i1,i2,i3,i4,i5,i6,i7
          character*6           :: c1,c2
        END SUBROUTINE

        SUBROUTINE mpi_davidson_master_ext_c(i1,i2,i3,i4,c1,c2,i5,i6,i7)
          integer               :: i1,i2,i3,i4,i5,i6,i7
          character*6           :: c1,c2
        END SUBROUTINE

        SUBROUTINE mpi_block_davidson_r(i1,i2,i3,i4,i5,c1)
          integer               :: i1,i2,i3,i4,i5
          character*6           :: c1
        END SUBROUTINE

        SUBROUTINE mpi_task_server (i1,i2)
          integer               :: i1,i2
        END SUBROUTINE

        SUBROUTINE build_tramo(i1,i2,ra1)
          integer               :: i1,i2
          real*8, dimension(:)  :: ra1
        END SUBROUTINE

        SUBROUTINE PRECONTRACT_PAR(i1,i2,ra1)
          integer               :: i1,i2
          real*8, dimension(:)  :: ra1
        END SUBROUTINE

      end interface
!
!---------------Execution ---------------------------------
!

!  _________________________ INTRO SECTION __________________________
! |
! |
! |
#if defined (VAR_MPI)
      if(MYPROC.EQ.MASTER) then   !in parallel: only master prints info
#endif

        WRITE(IW,*)' ________________________________________________ '
        WRITE(IW,*)'|                                                |'
        WRITE(IW,*)'|           Input control section                |'
        WRITE(IW,*)'|________________________________________________|'
        WRITE(IW,*)
        CALL print_input_info(iw,nrep) 
        WRITE(IW,*)

#if defined (VAR_MPI)
      endif
      WRITE(IW,'(A,I5)') ' Total # of nodes: ',NMPROC
      WRITE(IW,'(A,I5)') ' Current node # :  ',MYPROC
#endif

! |
! |
! |_________________________ END INTRO SECTION __________________________



!  _________________________ CALCULATION SECTION __________________________
! |
! |  in parallel: only master calcs ADC mat (is very fast),
! |  the nodes help with various
! |  tasks. Nodes will massively participate in Davidson.
! |
#if defined (VAR_MPI)
      if(MYPROC.EQ.MASTER) then 
#endif

        WRITE(IW,*)' ________________________________________________ '
        WRITE(IW,*)'|                                                |'
        WRITE(IW,*)'|           Calculation section                  |'
        WRITE(IW,*)'|________________________________________________|'
        WRITE(IW,*)
!
! initialize stack system (only master does this)
!
        qst_numofstacks = NUM_OF_STACKS
        qst_verbose = 0
        if(qstack_init(qst_numofstacks,qst_verbose).ne.NUM_OF_STACKS) 
     &     call quit('Qstack problem')
        WRITE(IW,*)'Stack system initialized with',
     &             NUM_OF_STACKS,' stacks'
        WRITE(IW,*)
!
!  debug purpose: print out IREPSPI array for absolute spinor numbers
!  IREPSPI(relative spinor number,irep,1=occ,2=virt) and
!  analyze the spinor sequence in a 2h2p column of the ADC matrix for all
!  possible final state symmetries (can be a long output)
!___________________________________________________________________
!|
!| check availability of transformed dipole moment integrals
!| required for transition moment calculations.
!|
        if(polprp_dotrmo) then
          write(*,*) 'TMs requested, entering property handler.'
          CALL property_handler(iobase,iobase+1,
     &                          propertyints,NSP,IREPSPI)
          inquire(file='ADCTRMOM',exist=isthere)
          if(isthere) then
            open(unit=iobase,file='ADCTRMOM',status='OLD')
            close(iobase,status='DELETE')
            write(*,*) 'Old ADCTRMOM file cleared.'
          endif
        endif
!|
!|
!|__________________________________________________________________
!
!  do PHPH block first (it can easily be kept in memory)
!
!   all symmetries in the HP-HP block are created at once.
!   Allocate PHPH block memory
!

        n1 = jvovo(nrep+1)
        allocate(daibj(n1*rcw),stat=ialloc)
        if(ialloc.ne.0) stop 'allocation problem daibj'
        daibj=0.0d0

        CALL BUILD_X_MAIN (daibj,eps,n1)
!
! diagonalize main block for Davidson start vectors
!
        if(carith) then
          CALL SOLVE_PHPH_C(daibj,n1,iobase)   !iobase is reusable!
        else
          CALL SOLVE_PHPH_R(daibj,n1,iobase)   !iobase is reusable!
        endif

#if defined (VAR_MPI)
      endif
#endif

!__________________________________________________
!|
!|  in this version the size of the buffer is fixed
!|
        INTBUF = 40*1024
!|
!|_________________________________________________


#if defined (VAR_MPI)
      if(MYPROC.EQ.MASTER) then 
!
!  trigger slaves. Check if all are awake.
!
        call interface_mpi_bcast(SERVER_HELLO,1,MASTER,
     &                           global_communicator)
#endif
!
! ______________ PRECONTRACTION of the VVVV and the VVOO integrals
! |
! |  precontract is only executed by master in case of a par run
! |  Slaves do their share in the service loop.
!
        write(*,*) 'Entering precontraction.'
        CALL PRECONTRACT_PAR(iobase,iw,eps)
        write(*,*) 'Precontraction finished.'

! ______________ SYMMETRY LOOP FOR COUPLING + SAT BLOCK___________
! |
! |   In case of parallel execution 
! |   only master loops through the irreps and the slaves provide 
! |   helper routines.
! |

        DO IRPLOOP = 1,NREP

! check for valid inputs

          desrep = polprp_statesym(IRPLOOP)
          if(desrep.EQ.0) cycle
          if(desrep.GT.NREP) then
             WRITE(IW,*) 'Chosen Symmetry',desrep,
     &       ' does not exist and is skipped!'
             CYCLE
          endif
          if(MVO(desrep).EQ.0) then
            WRITE(iw,*) 'No states in symmetry',desrep,repna(desrep)
            CYCLE
          endif

          write(*,*)
          write(*,*)
          write(*,*)
          write(*,*)
          write(*,*) '-----------------------------------------------',
     &       '-----------------------------------------------'
          write(*,*) '-----------------------------------------------',
     &       '-----------------------------------------------'
          write(*,*) '------  Now doing final state symmetry  ',desrep
          write(*,*) '-----------------------------------------------',
     &       '-----------------------------------------------'
          write(*,*) '-----------------------------------------------',
     &       '-----------------------------------------------'
          write(*,*)


! generate configuration file for this symmetry

          CALL WCONDAT_X(iobase,namecnffil,irecl,desrep,ireccnt)
          WRITE(IW,*) 'Configuration data written to ',namecnffil
          WRITE(IW,*) 'Record length in configuration file:',IRECL
          WRITE(IW,*) 'Number of records (=matrix dimension):',ireccnt
!
! calculate diagonal energies in satellite block for ADC-2 strict calculations
! (only master does it) Energies go on the stack and are fetched by
! XSATMAKE routine. Energies are calculated as real values  (diagonal) and
! are adapted accordingly in the complex routine.
!
          CALL XKENEDIAG(desrep,eps)

! generate coupling and sat block in this symmetry
! together with the helpers simultaneously started.

#if defined (VAR_MPI)
!
          inquire(file=namematfil,exist=isthere)
          if(isthere) then
            open(unit=iobase,file=namematfil,status='OLD')
            close(iobase,status='DELETE')
          endif
          open(unit=iobase,file=namematfil, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='NEW')
          n2 = jvovo(nrep+1)
          write(*,*)
          write(*,*) '**) Entering parallel driver for'
          if(carith) then
            write(*,*) '**) complex coupling and satellite block.'
            write(*,*)
            call mpi_bxcs_c(daibj,
     &                      n2,
     &                      eps,
     &                      desrep,
     &                      lencouplbl,
     &                      intbuf,
     &                      nbufs,
     &                      polprp_doextended,
     &                      polprp_writethr,
     &                      iobase,
     &                      iw)
          else
            write(*,*) '**) real coupling and satellite block.'
            write(*,*)
            call mpi_bxcs_r(daibj,
     &                      n2,
     &                      eps,
     &                      desrep,
     &                      lencouplbl,
     &                      intbuf,
     &                      nbufs,
     &                      polprp_doextended,
     &                      polprp_writethr,
     &                      iobase,
     &                      iw)
            write(*,*) '**) returned from real coupling and sat block.'
            write(*,*) '**) Used buffers of bufsize',intbuf,': ',nbufs
          endif

! ADC matrix is on disk and corresponding file is closed.

          close(iobase)

#else
!
! remind: in the serial case the number of buffers (NBUFS) is known
!         after matrix construction.

          inquire(file=namematfil,exist=isthere)
          if(isthere) then
            open(unit=iobase,file=namematfil,status='OLD')
            close(iobase,status='DELETE')
          endif
          OPEN(UNIT=IOBASE,FILE=namematfil, ACCESS='SEQUENTIAL',
     &       FORM='UNFORMATTED',STATUS='NEW')
          serialindicator = -654321
          Write(iobase) serialindicator
          write(*,*)
          write(*,*) '**) Entering serial driver for'
          IF(CARITH) THEN
            write(*,*) '**) complex coupling and satellite block.'
            write(*,*)
            CALL BUILD_X_COUSAT_C(DAIBJ,EPS,DESREP,LENCOUPLBL,
     &                            INTBUF,NBUFS,polprp_doextended,
     &                            polprp_writethr,iobase)
            CALL XDIAG_INI_C(iobase,INTBUF,NBUFS,
     &                       lencouplbl + MVO(desrep))
          ELSE
            write(*,*) '**) real coupling and satellite block.'
            write(*,*)
            CALL BUILD_X_COUSAT_R(DAIBJ,EPS,DESREP,LENCOUPLBL,
     &                            INTBUF,NBUFS,polprp_doextended,
     &                            polprp_writethr,iobase)
            CALL XDIAG_INI_R(iobase,INTBUF,NBUFS,
     &                       lencouplbl + MVO(desrep))
          ENDIF

! ADC matrix is on disk and corresponding file is closed.

          CLOSE(iobase)

          Write(*,*)
          Write(*,*) 'Serial run information:'
          Write(*,*) 'ADC matrix has',NBUFS,' buffers of length',INTBUF
          Write(*,*)
#endif

!
! check for internal consistency (we use various approaches for that)
! 
          ladc = lencouplbl + MVO(desrep)
          if(ladc.ne.ireccnt) then
            call quit(' *** Internal error. ADC matrix dim. mismatch')
          endif

#if defined (VAR_MPI)
! Master distributes ADC matrix when running in parallel

          write(iw,*)
          write(iw,*) '**** MASTER: distributing ADC matrix.'
          write(iw,*)
          call mpi_master_mat_dist(iobase,namematfil,irploop,
     &                              intbuf,nbufs,ladc,iw)
#endif

!  Davidson diagonalizer driver call. Works for
!  Serial and parallel. In case of parallel run
!  slaves perform partial matrix/vector multiplication.

          if(carith) then

            call mpi_block_davidson_c(iobase,
     &                                intbuf,
     &                                ladc,
     &                                MVO(desrep),
     &                                desrep,
     &                                namematfil)

          else

            call mpi_block_davidson_r(iobase,
     &                                intbuf,
     &                                ladc,
     &                                MVO(desrep),
     &                                desrep,
     &                                namematfil)

          endif

! Main ADC matrix file is deleted

          open(unit=iobase,file=NAMEMATFIL, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='unknown')
          close(iobase,status='DELETE')
!_________________________________________________________________
!|
!| serial or parallel/MASTER: calculate transition moments
!|

          if(polprp_dotrmo) then
            write(*,*) 'User asks for TM calculation',
     &                  ' (symmetry',desrep,')'
            write(*,*)
            if(propertyints) then
              call build_tramo(desrep,iobase,eps)
            else
              write(*,*) 'but required integrals do not exist.'
              write(*,*) 'No transition moments for irep',desrep,
     &                   'can be calculated.'
            endif
          endif
!|
!|________________________________________________________________

        ENDDO   !IRPLOOP

#if defined (VAR_MPI)
        write(iw,*) 'Master indicates finalization.'
        call interface_mpi_bcast(SERVER_DONE,1,MASTER,
     &                           global_communicator)
#endif
        deallocate(daibj)
        call qstack_shutdown()
! |
! |
! |________________ END SYMMETRY LOOP (MASTER) _________________________
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!   |
!   |
!   |
!   |     Slave execution part
!   |
!   |
!   |
#if defined (VAR_MPI)
      else
! _______________________________________________________
!|
!|
!|  ********************************************************
!|  *****  TASK SERVER FOR THE SLAVE NODES  ****************
!|  ********************************************************
!|
!|
        call mpi_task_server(iobase,intbuf)
        write(*,*) 'Node',myproc,' returns to main program'
        return

      endif
!|
!|  ********************************************************
!|  ******************* END OF TASK SERVER *****************
!|  ********************************************************
!|
!|_______________________________________________________


#endif

! |
! |
! |_________________________ END CALCULATION SECTION __________________________

      return
      end 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
