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
!
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PAMADC ()
!
!---------------Description-----------------------------------------------
!
!    Double group symmetry adapted adc code
!    This version is called from within DIRAC
!
!---------------Routines called-------------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : MP adapted from LV
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
#include "dcbmp2no.h"
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
!
!---------------Executable code--------------------------------------
!
      CALL CPUUSED(SEC0)
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      IN = LUCMD
      IW = LUPRI
!
! iface set to DIRAC6 as default, in line with having scheme 6 as default in moltra
! LV 2020: this set up is weird, but appears to give no problems here in contrast to
! polprp. If interface problems show up, look at the fixed code in polprp_driver
!     INTFCE = 'DIRAC     '
      INTFCE = 'DIRAC6    '
      IF( IAOMOSCHEME .eq. 4 ) INTFCE = 'DIRAC     '
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
!     ADC is not yet parallel!
!
!     Summon the slaves, who are waiting in the general menu routine.
!
!     CALL CCMSTI
!
!     Tell them the correct arithmetic (NZ is not initialized on the slaves
!     if you do not run HF !)
!
#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         call interface_mpi_BCAST_l0(CA,1,MASTER,global_communicator)
         call interface_mpi_BCAST_l0(MP2_DENRUN,1,MASTER,
     &                  global_communicator)
      END IF
!
!     if MP2-NO is enabled, tell the slaves about the correct interface
!
      IF( MP2_DENRUN )THEN
        call interface_mpi_BCAST(IAOMOSCHEME,1,MASTER,
     &                 global_communicator)
        IF( IAOMOSCHEME .eq. 4 ) INTFCE = 'DIRAC     '
        call interface_mpi_BCAST(INFO_STRING(1),62,MASTER,
     &                 global_communicator)
        call interface_mpi_BCAST(INFO_STRING(2),62,MASTER,
     &                 global_communicator)
        call interface_mpi_BCAST(INFO_STRING(3),62,MASTER,
     &                 global_communicator)
      END IF
#endif
!
!     Enter the generic ADCMAIN routine (master becomes now primus inter pares)
!

      CALL ADCMAIN_NEW(IN,IW,CA,INTFCE)

!
!     Return the slaves to the general menu routine.
!
!     CALL CCMSTE
!
      CALL CPUUSED(SEC)
      WRITE (IW,1000) SEC-SEC0,SEC0,SEC
      WRITE (IW,'(//80A1)') ('#',I=1,80)
      CLOSE(LUCMD)
!
 1000 FORMAT (//' CPU time (seconds) used in RELADC:',T50,F14.4/
     & ' CPU time (seconds) used before RELADC:',T50,F14.4/
     & ' CPU time (seconds) used in total sofar:',T50,F14.4//
     & '  --- Normal end of RELADC Run ---')
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADCMAIN_NEW (IN123,IW123,CA,INTFCE)
!
      use memory_allocator
      use interface_to_mpi
      use adc_cfg    !  in the general input reader the user data are provided.

      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Driver for Relativistic 1H/2H Propagators (no standalone)
!     RELADC now comprises single/double ionization + Fano-Stieltjes
!     in combination with Lanczos diagonalization.
!
!     Excitations based on the polarization propagator are now 
!     completely independent and are accessible via the
!     new (parallel) module POLPRP (see new manual entry)
!
!---------------Last modified------------------------------------------
!
!     Author : Markus Pernpointner
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
#include "adcinpt.inc"
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
      LOGICAL DOCCSETI
      REAL*8 CPU0,CPU1,CPU2,CPU3
!      integer :: verbose        = 0
!      integer :: ncore_elec_adc = 0
!      integer(kind=klongint) :: maxcor
!      integer(kind=klongint) :: mxcorr
!      integer(kind=klongint) :: mxcorr1
!      integer(kind=klongint) :: mword_limit
!      real(kind=kreal)       :: total_mem_new = 0
!      real(kind=kreal)       :: new_mword_limit = 0
!      integer                :: n,i,ntotsymm

      real(8), allocatable :: CC(:)
      real(8), allocatable :: testoooo(:)

      integer :: ncore_elec_adc = 0  ! new argument for RDSYMI, needed in RELCCSD, not in ADC.

      integer           :: maxcor,mxcorr,mxcorr1,mword_limit
      real*8            :: total_mem_new = 0   !number of bytes, if number of megawords is given in --aw
      real*8            :: new_mword_limit  = 0

      integer     :: MAXCORE    ! contains --mw * 1000 * 1000 (number of words, decimal mega)
      integer     :: maxcore_info(4)
      character*4 :: maxcore_infc(4)

      maxcore_info = 0
      maxcore_infc = (/'SETI','SIPS','DIPS','    '/)

      call legacy_lwork_get(MAXCORE)


!
!---------------Executable code--------------------------------------

!
!     Set the input and output units and a control string
!     In case relccsd is not called beforehand we need to agree
!     on which orbital energies we take (the original ones!)
!     Controlled by inpt.inc from relccsd.
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
!    say hello
!
      WRITE (IW,160)
 160  FORMAT (///10X,"************************************************"/
     &      10X,"****",40X,"****"/
     &      10X,"****       FOUR-COMPONENT         ",10X,"****"/
     &      10X,"****      PROPAGATOR MODULE       ",10X,"****"/
     &      10X,"****",40X,"****"/
     &      10X,"************************************************"///
     &      10X,"****  written by"/
     &      10X,"****  M. Pernpointner
     &          (University of Heidelberg 2004, 2010)"/
     &      10X,"****  For the 1P propagator cite: "/
     &      10X,"****  M. Pernpointner,
     &       J. Chem. Phys. 121, 8782 (2004)."/
     &      10X,"****  For the 2P propagator cite: "/
     &      10X,"****  M. Pernpointner, 
     &       J. Phys. B 43, 205102 (2010)."//)
!
!  first transfer interface variables to internal ADC variables
! ___________________________________________________________
!|
!|
!|
      adclevel = reladc_adclevel
      dosips   = reladc_dosips  
      dodips   = reladc_dodips  
      do i=1,32
         sipreps(i) = reladc_sipreps(i)
         dipreps(i) = reladc_dipreps(i)
      enddo
      readqkl  = reladc_readqkl 
      doconst  = reladc_doconst 
      doadcpop = reladc_doadcpop
      vconv    = reladc_vconv   
      adcthr   = reladc_adcthr  


      dofull   = reladc_dofull
      dolanc   = reladc_dolanc
!     dodavi   = reladc_dodavi

      sipiter  = reladc_sipiter 
      dipiter  = reladc_dipiter 
      do i=1,32
        sipeigv(i) = reladc_sipeigv(i) 
        dipeigv(i) = reladc_dipeigv(i) 
      enddo
      sipprnt  = reladc_sipprnt 
      dipprnt  = reladc_dipprnt 
      doincore = reladc_doincore
      lancmem  = reladc_lancmem
      do i=1,8
        adcprint(i)=0    !fill print level array 
      enddo
!|
!|
!|___________________________________________________________
!
! print overall job specification for the propagator run.
!
      write(IW,*)
      write(IW,'(T20,A)') 'Requested propagators:'
      write(IW,'(T20,A)') '----------------------'
      write(IW,*)
      if(dosips) write(IW,'(10X,A)') 'Single ionization.'
      if(dodips) write(IW,'(10X,A)') 'Double ionization.'
      write(IW,*)
      write(IW,'(T20,A)') 'Requested diagonalizers:'
      write(IW,'(T20,A)') '------------------------'
      write(IW,*)
      if(dofull) write(IW,'(10X,A)') 'Full (dim<5000!).'
      if(dolanc) write(IW,'(10X,A)') 'Lanczos.'
!     if(dodavi) write(IW,'(10X,A)') 'Davidson.'
!
      WRITE (IW,*)
      WRITE (IW,*)
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
      WRITE (IW,*) '*********************************'
      WRITE (IW,*) '***   memory control by user  ***'
      WRITE (IW,*) '*********************************'
      WRITE (IW,*)
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
         PRINT*," INCREASE MAXSP IN ADCMAIN TO",NSP
         CALL QUIT("NSP.GT.MAXSP")
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
!  Report memory in megabytes
!
!
!     Set up symmetry tables and parallel distribution arrays
!     Symtab is required for the pure diagonalization also.
!     The variables NELEC, NFROZ ... are defined in inpt.inc.
!
      write(iw,*)
      write(iw,'(T10,A)') '--------------------------------------------'
      write(iw,'(T10,A)') '     Resolving point group information:'
      write(iw,'(T10,A)') '--------------------------------------------'
      write(iw,*)
      CALL SYMTAB (NELEC,NFROZ,NINACT,NACT,NSP,EPS1,EPS,IREPSP,
     &             IREPSPI,.FALSE.,.TRUE.,NELEC_F1,NELEC_F2)
!
!  determine number of electrons and print general information
!
      NELECT = 0
      DO I = 1, 16
         NELECT = NELECT + NELEC(I)
      ENDDO
!
      WRITE (IW,1004) NELECT,NSP-NELECT,CARITH
1004  FORMAT (//" Number of electrons :",T40,I5
     &        /" Number of virtual spinors :",T40,I5
     &        /" Complex arithmetic mode :",T44,L1/)

!
!  initialize USEOE=T default also for RELADC
      USEOE   = .true.
!  because it relies on CCSETI and needs the orbital
!  energies from the original HF calculation
!  stefan: i am not sure whether you always want to set USEOE = .true.
!          but let us keep it the way the reladc wizards set it...
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

      IF(DOSIPS.OR.DODIPS) then
         DOCCSETI=.TRUE.
      ENDIF
!
!  call integral sorter, determine memory requirements and start right away.
!  
      IF(DOCCSETI) THEN
        ICALC = 1
        CALL CCSETI (ICALC,maxcore_info(1),EPS,NSP,IREPSP,IREPSPI)
        WRITE(IW,*) 'Memory requirements of ',
     &               maxcore_infc(1),maxcore_info(1)
        ICALC = 2
        CALL CCSETI (ICALC,maxcore_info(1),EPS,NSP,IREPSP,IREPSPI)
      ENDIF
! Elke
!      WRITE(iw,*) 'after ccseti'

      IF(DOSIPS) THEN
        ICALC = 2
        CALL ADCSNGL_MAIN (maxcore_info(2),NSP,IREPSPI,EPS)
        WRITE(IW,*) 'Memory requirements of ',
     &               maxcore_infc(2),maxcore_info(2)
      ENDIF

      IF(DODIPS) THEN
        ICALC = 2
        CALL ADCDBLE_MAIN (maxcore_info(3),NSP,IREPSPI,EPS)
        WRITE(IW,*) 'Memory requirements of ',
     &               maxcore_infc(3),maxcore_info(3)
      ENDIF

!<<<<<<< HEAD
!!
!! MXCORR is a cumulative maximum of individual memory requirements
!! ===> if any subroutine would need more core than MAXCORE we quit.
!!
!      IF (MXCORR.GT.MAXCOR) THEN
!         WRITE (IW,1002) MXCORR
!1002     FORMAT (" Not enough core, increase to at least",I20,
!     &        " Words")
!         CALL QUIT("Insufficient memory")
!      ENDIF
!      WRITE (IW,1001) MXCORR
!1001  FORMAT (" Memory used for active modules :",
!     &        T50,I10," 8-byte words")
!!
!!  We have sufficient memory. Now start the calculations
!!  Again, all modules are independent, resetting of mem. pointers
!!
!      ICALC = 2
!      ALLOC_ME = .TRUE.
!      KFREE = 1
!      LFREE = MXCORR
!!
!      write(*,*) 'Providing ',MXCORR,' to the act. modules.'
!      IF(DOSORT) THEN
!        CALL CCSETI (ICALC,MXCORR,EPS,NSP,IREPSP,IREPSPI)
!      ENDIF
!
!      write(*,*) 'MXCORR after allocated CCSETI:',MXCORR
!      write(*,*) 'We give the FULL memory to the ADC modules.'
!
!      LFREE = MAXCOR
!
!
!      IF(DOSIPS) THEN
!        CALL ADCSNGL_MAIN (ICALC,CC,KFREE,LFREE,MXCORR,NSP,IREPSPI,EPS)
!      ENDIF
!
!      IF(DODIPS) THEN
!        CALL ADCDBLE_MAIN (ICALC,CC,KFREE,LFREE,MXCORR,NSP,IREPSPI,EPS)
!      ENDIF
!
!=======
!>>>>>>> master
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADCSNGL_MAIN(MXCORR,NSP,IREPSPI,EPS)
!
!---------------Modules------------------------------------------------
!
      use adc_cfg
      use adc_fano
      use memory_allocator
      use adc_mat
C
      IMPLICIT INTEGER(A-Z)
C
C---------------Description--------------------------------------------
C
C     DHF-ADC(3) Main routine
C
C     Calculates ionization energies including correlation energy 
C     up to a certain order. Currently available: 
C             ADC-2, ADC-2X and ADC-3
C
C     The ADCLEVELS 1-3 correspond to the hierarchy above.
C
C     Author : MP
C
C
C---------------Calling variables--------------------------------------
C
      INTEGER               :: MXCORR,NSP,IREPSPI(NSP,16,2)
      REAL*8                :: EPS(*)
C
C---------------Common Blocks--------------------------------------
C
#include "dgroup.h"
#include  "../relccsd/files.inc"
#include  "adcinpt.inc"
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
C
C---------------Local variables--------------------------------------
C
      REAL(8), PARAMETER    :: autoev = 27.2113957D0

      INTEGER IDUMMY(2)
      REAL*8 RDUMMY(2)
      COMPLEX*16 CDUMMY(2)

      CHARACTER*5 NMSPECF
      PARAMETER(NMSPECF='SSPEC')
      CHARACTER*6 FILEN2,FILEN4,FILEN5
      PARAMETER(FILEN2='ADCCNF',FILEN4='ADCDIA',FILEN5='ADCMPK')
      CHARACTER*6 FILEN9
      PARAMETER(FILEN9='QKLVAL')

      CHARACTER*9 FILE2, FILE3, FILE4, FILE5
      CHARACTER*9 FILE9
      CHARACTER*15 CALTYP(0:3)

      DATA CALTYP/'Direct diagonal','ADC(2) strict  ',
     &            'ADC(2) extended','ADC(3)         '/
      INTEGER IRPLOOP,NBUFS,IBUFPOS
      REAL*8 ASYM
      INTEGER SYMMITER,IONIZLEVEL,IRECL

      real*8 ddot
      logical done

      integer(8) :: ivmax

C ***********************************************************************
C ************* Local arrays to be allocated dynamically  ***************
C ***********************************************************************

      integer, allocatable, dimension (:)    :: BUFIVV,ICRA,ICRB,IOI,
     & IOJ,IPIVOT,OOLT,VVLT

      real*8, allocatable, dimension (:)    :: ADCL, ADIAG1, ADIAG2, 
     & BDIAG1, BDIAG2, BUF1A, BUF1B, BUF1C, BUF1D, BUF4, BUF5A, BUF5B,
     & BUF5C, BUF9A, BUF9B, BUF9C, BUFBOO, BUFBVO2, BUFBVO, BUFBVV,
     & BUFIO, BUFO1, BUFO2, BUFO3, BUFQOO, BUFQVOA, BUFQVOI, BUFQVO,
     & BUFQVV, BUFSIG, BUFU1, BUFU2, BUFU3, BUFU4, BUFVOV,
     & BUFVOW, BUFVVVV, BUFW1, BUFW2, BUFW3,
     & BUFW4, BUFX1, BUFXVO, BUFXVV, BUFY1, BUFYVO, CKKSA, CKKS,
     & CSIGM, DIAG, EAJL, EJAB, FKKS
C
C--------------- description of buffer pointers 
C
      
C     CKKS,CKKSA    : buffer for hole/hole part. one auxiliary array is
C                     needed for ADC(3) since hermitian conjugate has 
C                     to be calculated
C
C     BUF1A - BUF1D : all purpose buffers for
c                     the <VV||OO> and <VO||OO> integrals
C                     must hold the following sorted integral lists:
c                     VV,OO: size IVVOOTT
c                     V,VOO: size IVVOOT
C                     VVO,O: size KVVOOT
C                     OVV,O: size LOOVVT
C                     VOO,V: size LVVOOT
C                     VO,OO: size IVOOOT
C                     VO,OO: size JVOOO  (no triangular OO storage!)
C                     complex storage allowed 
C                     all arrays are of roughly the similar size.
C                     therefore this way of allocation avoids a lot of
C                     different buffers.
C         ==> space : max(IVVOOTT,IVVOOT,KVVOOT,LOOVVT,LVVOOT,IVOOOT,JVOOO) 
C         ==> NDIMT3
C
C     BUF4:  holds the array of epsilon sums for the h/h block in
C            a specific symmetry only. ALWAYS REAL !
C        ==> dimension: NDIMT4
C
C     BUF5A/B  :  same as BUF1A - BUF1D, maximum dimension also NDIMT3
C
C     BUF9A - BUF9C  holds VOVO integrals of length IVOVO(NREP+1) =
C                    JVOVO(NREP+1) !
C
C     BUFO1 - BUFO3  buffers for OOOO integrals.
C                    must hold the following sorting types:
C                    I>J,K>L: size IOOOOTT
C                    I>JK,L:  size KOOOOT
C                    IJ,KL:   size JOOOO
C
C     BUFU1 - U4 : holds coupling coefficients of U affinity/hole index
C                  U^+_i,jab . length: NO * NOVVT
C
C     BUFW1 - W4 : holds coupling coefficients for U ionization/particle
C                  U^-_d,akl . length: NV * NVOOT
C
C     BUFX1:   holds coupling coefficients for U affinity/particle
C              U+_d,jab   length: NV*NOVVT
C
C     BUFY1:   holds coupling coefficients for U ionization/hole
C              U+_i,akl   length: NO*NVOOT
C              needed in *** ALL *** ADC levels
C
C       ** important: all four coupling U matrices are stored
C                     as U_xyz,p in memory !
C
C     EAJL holds the K_ajl,ajl array of dimension MXNVOOT. ALWAYS REAL !
C
C     EJAB holds the K_jab,jab array of dimension MXNOVVT. ALWAYS REAL !
C
C     OOLT holds OOT lookup table for each irrep
C
C     VVLT holds VVT lookup table for each irrep
C
C     ADIAG1,2 hold diagonal elements for the V ionization iteration
C              length: MXNVOOT, complex allowed
C
C     BDIAG1,2 hold diagonal elements for the V affinity iteration
C              length: MXNOVVT, complex allowed
C
C     BUFVVVV holds space for all (!) VVVV integrals. In this implementation
C             an inevitable nuisance. After exiting subroutine this
C             space is, of course, released again.
C  
C     BUFVV1/2 apply in ADC(3) and divide the remaining memory space
C              in two halves, each large enough for the <VO||VV>
C              integrals, and in total large enough for the <VV||VV>
C              integrals.
C
C     BUFQOO   buffer for Q_oo, length MOO(1)
C     BUFQVO   buffer for Q_vo, length MVO(1)
C     BUFQVOI  buffer for Q_vo, length MVO(1) ionization part only
C     BUFQVOA  buffer for Q_vo, length MVO(1) affinity part only
C     BUFQVV   buffer for Q_vv, length MVV(1)
C
C     BUFBOO   buffer for B_oo, length MOO(1)
C     BUFBOV   buffer for B_ov, length MOV(1) : probably not needed !
C     BUFBVO   buffer for B_vo, length MVO(1)
C     BUFBVV   buffer for B_vv, length MVV(1)
C
C  --------- variables for sigma-infinity  code:
C
C     CSIGM    buffer for sigma coefficient matrix in the complex case
C              only. dimensions 2*MVO(1) * 2*MVO(1)
c
c   save the original memory space for the diagonalizers


      IONIZLEVEL=1
c
c  determine necessary dimension for the all-purpose arrays
c  BUF1A...BUF1D, and BUF5A, BUF5B
c
      NDIMT3 = MAX0 (IVVOOTT(NREP+1),IVVOOT(NREP+1))
      NDIMT3 = MAX0 (NDIMT3,KVVOOT(NREP+1))
      NDIMT3 = MAX0 (NDIMT3,LOOVVT(NREP+1))
      NDIMT3 = MAX0 (NDIMT3,LVVOOT(NREP+1))
      NDIMT3 = MAX0 (NDIMT3,IVOOOT(NREP+1))
      NDIMT3 = MAX0 (NDIMT3,JVOOO(NREP+1))

      NDIMTO = MAX0 (IOOOOTT(NREP+1),KOOOOT(NREP+1))
      NDIMTO = MAX0 (NDIMTO,JOOOO(NREP+1))
C
C  determine array dimensions for the  > largest < irrep and
C  combinations of irreps. this is necessary because for sigma(inf) 
C  the program loops over ALL symmetries !!!
C  if we do ADC-3 contributions from ALL irreps come in no matter
C  if we calculate only one specific irrep.
C  We have to allocate arrays such that the largest irrep can fit
C
      MXNO    = 0
      MXNV    = 0
      MXNVOOT = 0
      MXNOVVT = 0

      NDIMT4  = 0
      LADCM   = 0
      MXICRB  = 0

      LBUFUAH = 0
      LBUFUAP = 0
      LBUFUIH = 0
      LBUFUIP = 0

      DO 15 KREP = 1,NREP

        IF(NO(KREP).GT.MXNO)       MXNO=NO(KREP)
        IF(NV(KREP).GT.MXNV)       MXNV=NV(KREP)
        IF(NVOOT(KREP).GT.MXNVOOT) MXNVOOT=NVOOT(KREP)
        IF(NOVVT(KREP).GT.MXNOVVT) MXNOVVT=NOVVT(KREP)


        MX1 = NVVOT(KREP)*NO(KREP)*NO(KREP)
        MX2 = NVOOT(KREP)+NO(KREP)
        MX3 = NOVVT(KREP)+NO(KREP)
        IF(MX1.GT.NDIMT4) NDIMT4 = MX1
        IF(MX2.GT.LADCM)  LADCM = MX2
        IF(MX3.GT.MXICRB) MXICRB = MX3


        LKUAH = NO(KREP) * NOVVT(KREP)
        LKUAP = NV(KREP) * NOVVT(KREP)
        LKUIH = NO(KREP) * NVOOT(KREP)
        LKUIP = NV(KREP) * NVOOT(KREP)
        IF(LKUAH.GT.LBUFUAH) LBUFUAH = LKUAH
        IF(LKUAP.GT.LBUFUAP) LBUFUAP = LKUAP
        IF(LKUIH.GT.LBUFUIH) LBUFUIH = LKUIH
        IF(LKUIP.GT.LBUFUIP) LBUFUIP = LKUIP

 15   CONTINUE

      NCKKS   = MXNO*MXNO
      INTBUF = 5*1024*1024  !output buffer is now 5 MWORDS (40 MByte)
!
!  allocate memory needed in all types of ADC calculations
!
      allocate(CKKS(NCKKS*RCW))
      allocate(BUF1A(NDIMT3*RCW))
      allocate(BUF1B(NDIMT3*RCW))
      allocate(BUF4(NDIMT4))
      allocate(BUF5A(NDIMT3*RCW))
      allocate(BUF5B(NDIMT3*RCW))
      allocate(EAJL(MXNVOOT))
      allocate(ICRA(6*LADCM))
      allocate(ADCL(LADCM*RCW))
      allocate(DIAG(LADCM*RCW))
      allocate(BUFIO(INTBUF*RCW))
      allocate(BUFY1(LBUFUIH*RCW))
      allocate(IOI(INTBUF))
      allocate(IOJ(INTBUF))
!
!  For ADC-2X and ADC-3 we need one more OOOO and VOVO buffer 
!  for the satellite block and the occ/occ lookup table
!
      IF(ADCLEVEL.GE.2) THEN
        OOLTL=MXNO*MXNO*NREP
        allocate(OOLT(OOLTL))
        allocate(BUFO1(NDIMTO*RCW))
        allocate(BUF9A(NV4*RCW))
      ENDIF
!
!  For ADC(3) we need additional buffers.
!
      IF(ADCLEVEL.EQ.3) THEN
        NX1=4*MVO(1)*MVO(1)*RCW
        NX2=MXNV*MXNV*NREP
        allocate(CKKSA(NCKKS*RCW))
        allocate(FKKS(NCKKS*RCW))
        allocate(BUF1C(NDIMT3*RCW))
        allocate(BUF1D(NDIMT3*RCW))
        allocate(BUF5C(NDIMT3*RCW))
        allocate(BUF9B(NV4*RCW))
        allocate(BUF9C(NV4*RCW))
        allocate(BUFO2(NDIMTO*RCW))
        allocate(BUFO3(NDIMTO*RCW))
        allocate(BUFU1(LBUFUAH*RCW))
        allocate(BUFU2(LBUFUAH*RCW))
        allocate(BUFU3(LBUFUAH*RCW))
        allocate(BUFU4(LBUFUAH*RCW))
        allocate(BUFW1(LBUFUIP*RCW))
        allocate(BUFW2(LBUFUIP*RCW))
        allocate(BUFW3(LBUFUIP*RCW))
        allocate(BUFW4(LBUFUIP*RCW))
        allocate(BUFX1(LBUFUAP*RCW))
        allocate(EJAB(MXNOVVT))
        allocate(ADIAG1(MXNVOOT*RCW))
        allocate(ADIAG2(MXNVOOT*RCW))
        allocate(BDIAG1(MXNOVVT*RCW))
        allocate(BDIAG2(MXNOVVT*RCW))
        allocate(BUFQOO(MOO(1)*RCW))
        allocate(BUFQVO(MVO(1)*RCW))
        allocate(BUFQVOI(MVO(1)*RCW))
        allocate(BUFQVOA(MVO(1)*RCW))
        allocate(BUFQVV(MVV(1)*RCW))
        allocate(BUFBOO(MOO(1)*RCW))
        allocate(BUFBVO(MVO(1)*RCW))
        allocate(BUFXVO(MVO(1)*RCW))
        allocate(BUFYVO(MVO(1)*RCW))
        allocate(BUFBVV(MVV(1)*RCW))
        allocate(BUFXVV(MVV(1)*RCW))
        allocate(BUFIVV(MVV(1)))
        allocate(BUFBVO2(2*MVO(1)*RCW))
        allocate(BUFSIG(MOO(1)*RCW))
        allocate(CSIGM(NX1))
        allocate(ICRB(6*MXICRB))
        allocate(IPIVOT(2*MVO(1)))
        allocate(VVLT(NX2))
      ENDIF
C
C  The remaining memory (major part) has to be enough for
C  a) the VV?VV integrals (!). indeed this is not good but
C     a complicated access pattern in the sigma calculations
C     necessitates that.
C
C  ******************* up to here the fixed buffers are all allocated
C  the remaining large chunks of memory are allocated and freed on
C  demand, because it is only needed temporarily. The KFRDYN
C  variable holds the beginning of this area.
C  ******************* 

      IVMAX = 0
      IF(ADCLEVEL.EQ.3) THEN
!
! determine max space for the large <VV|VV>,  <VO|VV> classes
! and inform user about that. Attention! If careless usage is
! done her, machine can start swapping in an unpleasant manner!
!
        IVMAX=MAX0(IVVVVTT(NREP+1),2*IVOVVT(NREP+1))
        IVMAX=MAX0(IVMAX,2*JVVVO(NREP+1))
        IVMAX = IVMAX*RCW
        WRITE(IW,*) ' WORDS/MWORDS dynamically allocated for the'
        WRITE(IW,*) ' <VV|VV> resp. <VO|VV> classes:',IVMAX,
     &              IVMAX/(1024*1024)
      ENDIF
C
C----------------------------------------------------------------------
C------------------------------- START ADC CALCULATION ----------------
C----------------------------------------------------------------------
C
      WRITE(IW,*)
      CALL PST('Parameters for the single ionization run:+')
      WRITE(IW,*)
      write(IW,'(10X,A,T50,I8)') 'Number of spinors:',NSP
      write(IW,'(10X,A,T50,A8)') 'ADC treatment:',CALTYP(adclevel)
      write(IW,'(10X,A,T50,L1)') 'Read previous Qkl data:',readqkl
      write(IW,'(10X,A,T50,L1)') 'Incl. constant diagrams:',doconst
      write(IW,'(10X,A,T50,F18.6)') 'Convergence for CD:',vconv
      write(IW,'(10X,A,T50,F18.6)') 'Threshold for screen output:',
     &      sipprnt
      write(IW,'(10X,A)') 'Requested SIP symmetries (sipreps):'
      n=0
      do i=1,32
        n=n+sipreps(i)
      enddo
      if(n.eq.0) then
        write(IW,'(10X,A)') 'None. Program will do it automatically.'
      else
        write(IW,'(10X,16I4)') (sipreps(i),i=1,16)
        write(IW,'(10X,16I4)') (sipreps(i),i=17,32)
      endif
      write(IW,'(10X,A)') 'Requested SIP eigenvectors:'
      n=0
      do i=1,32
        n=n+sipeigv(i)
      enddo
      if(n.eq.0) then
        write(IW,'(10X,A)') 'None.'
      else
        write(IW,'(10X,16I4)') (sipeigv(i),i=1,16)
        write(IW,'(10X,16I4)') (sipeigv(i),i=17,32)
      endif
C
C
C  next we check if user wants to read in previously calculated Q_kl
C  data. Otherwise we have to calculate them. JSTAT will give
C  information about status of availability of the Q_kl data.
C  if readqkl is false, jstat = 1 and the data will be calculated !
C  if reading is successful jstat will switch to 0 !
C
      JSTAT = 1
      IF( (ADCLEVEL.EQ.3).AND.(READQKL)) THEN
        WRITE(IW,1040)
        CALL RWQKL(BUFQOO,BUFQVO,BUFQVV,LFA,JSTAT,
     &             ITAPADC+9,FILEN9)
        IF(JSTAT.EQ.0) THEN
          WRITE(IW,1041)
        ELSE
          WRITE(IW,1042)
        ENDIF
      ENDIF

 1040 FORMAT (/5X,"Retrieving Q_kl data ...")
 1041 FORMAT (/5X,"Q_kl successfully loaded")
 1042 FORMAT (/5X,"Q_kl not available on disk. Will be computed")
C
C  we have to determine symmiter variable that determines kind of
C  Kramers degeneracy depending on variable NZ
c  If we work in real symmetries we have degenerate
c  Kramers pairs. We can therefore skip one Kramers partner
c  and copy the corresponding Qkl array. 
C  In the case of complex symmetries (C_s C_i C_1) we have to
c  distinguish three cases:
c
c  C_s group: no real MO integrals, inherently complex but
c             existence of Kramers degeneracy due to 
c             TWO degenerate fermion irreps:  symmiter = 2
c  C_i group: no real MO integrals, inherently complex 
c             two fermionic nondegenerate irreps E1g/E1u
c             not related by Kramers symmetry: symmiter =1
c  C_1 group: same as C_i but only one nondegenerate irrep E1
c             NREP=1, symmiter=1
c
c  In case of spinfree or Levy-Leblond calculations we have
c  bosonic irreps exclusively, therefore real arithmetic
c  and an irrep step of 1 (GROUPTYPE = 8)
c
      IF(NZ.eq.1) THEN
        IF(SPFR.eqv..FALSE.) THEN
          SYMMITER = 2
          WRITE(IW,*) 'Real fermionic group (D2h, C2v, linear)'
          WRITE(IW,*) 'Only one Kramers partner calculated!'
        ELSE
          SYMMITER = 1
          WRITE(IW,*) 'Real bosonic group (Spinfree/Lev.-Lebl.)'
          WRITE(IW,*) 'No Kramers skipping'
        ENDIF
      ELSE IF(NZ.eq.2) THEN
        IF(SPFR.eqv..FALSE.) THEN
          SYMMITER = 2
          WRITE(IW,*) 'Complex fermionic group with Kramers pairs'
          WRITE(IW,*) 'Only one Kramers partner calculated!'
        ELSE
          SYMMITER = 1
          WRITE(IW,*) 'Real bosonic group (Spinfree/Lev.-Lebl.)'
          WRITE(IW,*) 'No Kramers skipping'
        ENDIF
      ELSE IF(NZ.eq.4) THEN
        SYMMITER = 1
        WRITE(IW,*) 'Quaternionic fermionic group. No Kramers pairs.'
        WRITE(IW,*) 'Each symmetry treated explicitly.'
      ELSE
        CALL QUIT('NZ unassigned. Should not happen!')
      ENDIF

C
C  if we do ADC-3 we need the Q_kl in  **ALL** symmetries ==> we have to
C  loop over all symmetries for the corresponding entities. This is due
C  to the coupling of symmetries in the b_pq calculation !
C
C  start Q_kl calculation if JSTAT=1 and ADC-3 is requested
C  if DOCONST = FALSE then NO constant diagrams will be calculated
C  only the other ADC-3 contributions will be inserted.
C

      IF( (ADCLEVEL.EQ.3).AND.
     &    (JSTAT.EQ.1)   .AND.
     &    (DOCONST.EQV..TRUE.)   ) THEN

        CALL PST('Starting calculation of constant diagrams+')

        CALL XCOPY(MOO(1),A0,0,BUFQOO,1)
        CALL XCOPY(MVO(1),A0,0,BUFQVO,1)
        CALL XCOPY(MVO(1),A0,0,BUFQVOI,1)
        CALL XCOPY(MVO(1),A0,0,BUFQVOA,1)
        CALL XCOPY(MVV(1),A0,0,BUFQVV,1)

c------------------------------------------------------------
c---------- loop over symmetries for the Qkl calculation ----
c---------- ONLY IF DOCONST = TRUE !! -----------------------
c------------------------------------------------------------

        DO 700 KREP=1,NREP,SYMMITER
          IF( (NO(KREP).EQ.0).AND.(NV(KREP).EQ.0)) THEN
            WRITE(IW,*) 'No Qkl contribution in symmetry',KREP
            GOTO 700
          ENDIF
c
c in the constant diagram loop we need the VVVV integral buffer
c
c we allocate it only as long as we need it !
c
          allocate(BUFVVVV(IVMAX))
c
c -------  prepare V+/- vector iteration
c
          CALL UAFFHO(BUFU1,
     &                BUF1A,BUF1B,BUF1C,
     &                BUF9A,BUF9B,BUF9C,
     &                BUFO1,EPS,KREP)
c
          CALL UIONPA(BUFW1,
     &                BUF1A,BUF1B,BUF1C,
     &                BUF9A,BUF9B,BUF9C,
     &                BUFVVVV,EPS,KREP)
c
c  UAFFHO/UIONPA are for VAFFITER/VIONITER
c  and need to be complex conjugated for the correct BKC !
c  for VAFFITER an incore version is not feasible, we tried....
c
          NX = NO(KREP)*NOVVT(KREP)
          IF(CARITH) CALL CONJUGA(NX,BUFU1,1)
          NX = NV(KREP)*NVOOT(KREP)
          IF(CARITH) CALL CONJUGA(NX,BUFW1,1)

          CALL MAKEEKOVV(EJAB,EPS,KREP,ICRB)
          IF(CARITH) THEN
            CALL CVAFFITER (BDIAG1,BDIAG2,BUFVVVV,BUF9A,
     &                      BUFU1,BUFU2,BUFU3,BUFU4,
     &                 EJAB,EPS,ICRB,VVLT,MXNV,KREP,VCONV)
          ELSE
            CALL VAFFITER (BDIAG1,BDIAG2,BUFVVVV,BUF9A,
     &                     BUFU1,BUFU2,BUFU3,BUFU4,
     &                 EJAB,EPS,ICRB,VVLT,MXNV,KREP,VCONV)
          ENDIF
C
C  after VAFFITER the converged V+ vectors are in BUFU4
C
          CALL MAKEEKVOO(EAJL,EPS,KREP,ICRA)
          IF(CARITH) THEN
            CALL CVIONITER (ADIAG1,ADIAG2,BUFO1,BUF9A,
     &                      BUFW1,BUFW2,BUFW3,BUFW4,
     &                  EAJL,EPS,ICRA,OOLT,MXNO,KREP,VCONV)
          ELSE
            CALL VIONITER (ADCL,DIAG,BUFO1,BUF9A,
     &                     BUFW1,BUFW2,BUFW3,BUFW4,
     &                  EAJL,EPS,ICRA,OOLT,MXNO,KREP,VCONV)
          ENDIF
c
C  after VIONITER the converged V- vectors are in BUFW4 and we can release
c  the huge VVVV buffer BUFVVVV

          deallocate(BUFVVVV)
C
C  next we need to provide VOVV buffers
C
          NX1 = JVVVO (NREP+1) 
          allocate(BUFVOV(NX1*RCW))
          allocate(BUFVOW(NX1*RCW))
C
C  next we generate the remaining two U arrays in order to
C  determine the Q_kl entities and release the VOVV space afterwards
C
          CALL UAFFPA(BUFX1,
     &                BUFVOV,BUFVOW,
     &                BUF1A,BUF1B,BUF1C,BUF9A,
     &                EPS,KREP)

C
          CALL UIONHO(BUFY1,
     &                BUF1A,BUF1B,BUF1C,
     &                BUF9A,BUFVOV,
     &                EPS,KREP)
c
c  complex conjugate UAFFPA/UIONHO for the correct BKC in CMAKE_QKL !
c  because the Q_kl are produced in |><| BKC !
c  afterwards release memory for the VOVV storage.
c
          NX = NV(KREP)*NOVVT(KREP)
          IF(CARITH) CALL CONJUGA(NX,BUFX1,1)
          NX = NO(KREP)*NVOOT(KREP)
          IF(CARITH) CALL CONJUGA(NX,BUFY1,1)

          deallocate(BUFVOV)
          deallocate(BUFVOW)

C
C from this point we have all necessary U coupling blocks
C and the iterated affinity/ionization V vectors for symmetry KREP.
C
C the coupling matrices and V vectors are transferred to the QKL routine
C                   in the following order:
C  ****   Uion,hole // Uaff,part // Vion,part // Vaff,hole ****
C  if the group structure allows for that we duplicate the Q_kl values
C  for KREP+1 symmetry afterwards.
C
          IF(CARITH) THEN
            CALL CMAKE_QKL (BUFY1,BUFX1,BUFW4,BUFU4,
     &                BUFQOO,BUFQVO,BUFQVOI,BUFQVOA,
     &                BUFQVV,EPS,KREP,ADCPRINT)
          ELSE
            CALL RMAKE_QKL (BUFY1,BUFX1,BUFW4,BUFU4,
     &                BUFQOO,BUFQVO,BUFQVOI,BUFQVOA,
     &                BUFQVV,EPS,KREP,ADCPRINT)
          ENDIF
c
c duplicate if possible
c
          IF(SYMMITER.EQ.2) THEN
            IF(CARITH) THEN
              CALL CQKLDUP(BUFQOO, BUFQVO, BUFQVV, KREP)
            ELSE
              CALL QKLDUP(BUFQOO, BUFQVO, BUFQVV, KREP)
            ENDIF
            WRITE(IW,*) ' Q_kl duplication could be performed.'
          ELSE
            WRITE(IW,*) ' No Q_kl duplication possible !'
          ENDIF
c
c  integral buffer is released. nothing additionally allocated.
c
 700    CONTINUE
C
C from here we have the complete Q_kl arrays and we write them out to
C disk
C 
        CALL RWQKL(BUFQOO,BUFQVO,BUFQVV,LTR,JSTAT,
     &             ITAPADC+9,FILEN9)

      ENDIF
C
C  end of Q_kl calculation
C

C
C  next we calculate the B_pq and Sigma(infinity) if we are in ADC-3
C  Of course, these entities are only available if DOCONST was set true.
C
      IF(  (ADCLEVEL.EQ.3).AND.DOCONST)  THEN
C
C in MAKE_BPQ1 we do all calculations requiring arrays up to the VOVV type
C MAKE_BPQ1 is written for real and complex case. The Q_kl were
C generated in |><| BKC ==> NO additional complex conjugation.
C
        NX1 = JVVVO (NREP+1)
        allocate(BUFVOV(NX1*RCW))
        allocate(BUFVOW(NX1*RCW))

        CALL MAKE_BPQ1(BUFQOO,BUFQVO,BUFQVV,
     &                 BUFBOO,BUFBVO,
     &                 BUFO1,BUFO2,BUF1A,BUF1B,
     &                 BUF1C,BUF9A,BUF9B,
     &                 BUFVOV,BUFVOW,ADCPRINT)

        deallocate(BUFVOV)
        deallocate(BUFVOW)

C
C the B_pq are available right now and we can solve the system
C of equations for Sig(infty)
C BKC:
C constant diagrams are created in the <||> order, the ADC matrix,
C however has |><| convention in the hh block. in the complex case we
C therefore complex conjugate the resulting constant diagrams.
C
        IF(CARITH) THEN
          CALL CSIGSOLV(BUFBOO,BUFBVO,BUFXVO,
     &                  BUFYVO,BUFBVO2,
     &                  BUF1A,BUF1B,
     &                  BUF1C,BUF1D,
     &                  BUF9A,BUF9B,BUF9C,
     &                  EPS,IPIVOT,BUFSIG,
     &                  CSIGM,ADCPRINT)
        ELSE
          CALL SIGSOLV(BUFBOO,BUFBVO,
     &               BUF1A,BUF1B,BUF1C,BUF1D,
     &               BUF9A,BUF9B,BUF9C,
     &               EPS,IPIVOT,BUFSIG,ADCPRINT)
        ENDIF
C
C Sigma(infty) is calculated and available in CC(BUFSIG)
C
      ENDIF

C -----------------------------------------------------------------
C
C  from here all the blocks for the ADC(2) and ADC(2x) calculations
C  are generated and combined (if requested) with the previously
C  Sigma(infty) contributions.
C  the loop is over the user-requested symmetries !
C
C  If the user does not enter any specific symmetries,
C  RELADC determines all symmetries there is work to be done for. Only if user
C  specifies symmetries explicitly these ones are then calculated in the ordinary
C  way.
C
C  Diagonalization is done immediately after matrix construction
C  As a consequence, filenames for the diag/matrix file are not
C  labeled by symmetry number anymore.
C
C                  *************************************************
C                  ******* LOOP OVER USER-DEFINED SYMMETRIES *******
C                  ******* LOOP OVER USER-DEFINED SYMMETRIES *******
C                  ******* LOOP OVER USER-DEFINED SYMMETRIES *******
C                  *************************************************
C
c
c determine if user has given explicit symmetries
c
      n = 0
      do i = 1,nrep
        n = n + sipreps(i)
      enddo
      if(n.eq.0) then
        write(iw,*) 'User has not given explicit symmetries.'
        write(iw,*) 'Kramers degeneracy counter:',symmiter
        write(iw,*) 'Program calculates the following SIP symmetries:'
        write(iw,*)
        icount=1
        do i=1,nrep,symmiter
          sipreps(icount) = i
          write(iw,*) '   No: ',icount,'  symmetry: ',i
          icount = icount + 1
        enddo
      else
        write(iw,*) 'User has made explicit symmetry choices.'
        write(iw,*) 'Program calculates the following SIP symmetries:'
        write(iw,*)
        icount=1
        do i=1,nrep
          if(sipreps(i).ne.0) then
            write(iw,*) '   No: ',icount,'  symmetry: ',sipreps(i)
            icount = icount + 1
          endif
        enddo
      endif
      write(iw,*) 'Program calculates',icount-1,' symmetries.'
      write(iw,*)
      
! *********************************
! *****  main symmetry loop  ******
! *********************************

      DO 900 IRPLOOP=1,NREP

! *********************************
! *****  main symmetry loop  ******
! *********************************

        IF(sipreps(IRPLOOP).EQ.0) GOTO 900
        IF(sipreps(IRPLOOP).GT.NREP) THEN
           WRITE(IW,*) 'Chosen Symmetry',sipreps(IRPLOOP),
     &     ' does not exist and is skipped!'
           GOTO 900
        ENDIF
c
c now we actually found a symmetry and start. Further tests
c avoid crashing. KREP is the considered symmetry.
c
        KREP=sipreps(IRPLOOP)
c
        IF(NO(KREP).EQ.0) THEN
          WRITE(IW,290) KREP
          GOTO 900
        ENDIF
 290  FORMAT (" No occupied orbitals in IRREP",I2," ...skipping")

        IF(NVOOT(KREP).EQ.0) THEN
          WRITE(IW,295) KREP
          GOTO 900
        ENDIF
 295  FORMAT (" No VOOT orbitals in IRREP",I2," ...skipping")

        IF(NVVOT(KREP).EQ.0) THEN
          WRITE(IW,298) KREP
          GOTO 900
        ENDIF
 298  FORMAT (" No VVOT orbitals in IRREP",I2," ...skipping")

        WRITE(IW,300) KREP,REPNA(KREP)
 300    FORMAT (/" ******************************************"
     & /" ***"/" ***"/" *** IRREP",I6,3X,"(",A4,")"/
     & " ***"/" ***"/" ******************************************"//)
C
C  Write the ADC files into the corresponding files. 
C  Each symmetry is diagonalized immediately.
C

C
C  set start number of written buffers and initial buffer pointer
C
        NBUFS=0
        IBUFPOS=0
C
C  clear h/h and h/2h1p (Non-Dyson) block
C
        N=NO(KREP)*NO(KREP)
        CALL XCOPY(N,A0,0,CKKS,1)
        N=NO(KREP)*NVOOT(KREP)
        CALL XCOPY(N,A0,0,BUFY1,1)

        LADC = NVOOT(KREP) + NO(KREP)
        WRITE(IW,'(A,2I8)') 
     &      ' Hole space & dimension of ADC matrix in this irrep: ',
     &      NO(KREP),LADC
C
C  ***************************  start h/h block  ****************************
C  ... BKC:  |><|
C
!! Elke:
!!        WRITE(iw,*) 'ckks'
!!        WRITE(iw,'(1X,ES14.5)') (cc(ixx), ixx=ckks,ckks+no(krep)-1)
!
        CALL MAKE_HH2(CKKS,BUF1A,BUF1B,EPS,BUF4,KREP)
C
C  If ADC(3) is activated: compute C(3) and add Sigma(inf).
C  we also need space for the VVVV integrals.   
C  S4P flag is always false because the Sigma(4+) algorithm will not
C  go into this implementation.
C
        IF(ADCLEVEL.EQ.3) THEN
          allocate(BUFVVVV(IVMAX))
          CALL MAKE_HH3(KREP,CKKS,CKKSA,EPS,BUFVVVV,
     &                  BUF1A,BUF1B,BUF1C,BUF1D,
     &                  BUF9A,BUF9B,BUF9C,BUFO1,
     &                  FKKS,.false.)
          deallocate(BUFVVVV)
  
          CALL CKKASYM(CKKS,CKKS,NO(KREP),ASYM)
          WRITE(IW,*) 'h/h hermiticity deviation:',ASYM

          IF(DOCONST) THEN
            IX1=JJOO(KREP,KREP)*RCW + 1
            CALL ADDARR(NO(KREP),BUFSIG(IX1),BUFSIG(IX1),CKKS,CKKS)
          ENDIF
        ENDIF
c
c  check deviation from hermiticity (should be approx. zero)
c
        CALL CKKASYM(CKKS,CKKS,NO(KREP),ASYM)
        WRITE(IW,*) 'h/h hermiticity deviation:',ASYM
        WRITE(IW,*)
C
C  compute 1h/2h-1p block contributions for ADC2 and ADC2X
C
C  The scaling factor of -1 arises because in UIONHO the OSWALD-type
C  formula is implemented that is the negative of the Schirmer/Trofimov
C  formula. One further CC is necessary to comply with the BKC
C  convention chosen for the overall ADC matrix.
C
        IF(ADCLEVEL.EQ.1 .OR. ADCLEVEL.EQ.2) THEN
          CALL MAKE_HHP2(BUF5A,BUF5B,BUFY1,KREP)
        ELSE IF(ADCLEVEL.EQ.3) THEN

          NX1 = JVVVO (NREP+1)
          allocate(BUFVOV(NX1*RCW))

          CALL UIONHO(BUFY1,BUF1A,BUF1B,BUF1C,
     &                BUF9A,BUFVOV,EPS,KREP)
          NX=NO(KREP)*NVOOT(KREP)
          CALL XSCAL(NX,-A1,BUFY1,1)
          IF(CARITH) CALL CONJUGA (NX,BUFY1,1)
 
          deallocate(BUFVOV)

        ELSE
          CALL QUIT('Illegal ADCLEVEL value for coupling block calc.')
        ENDIF
C
C  write the off-diagonal h/h 
C  and the h/2hp buffer to disk but NOT the diagonal entries since
C  there are more diagonal elements to be calculated in MAKEHPHP !
C  diagonal is written out in MAKEHPHP
C  all other elements are finished and can be written out.
C
        WRITE(iw,*) 'CARITH = ', CARITH

!
! This is the big Fano-IF-clause introduced by Elke
!
        IF (reladc_dofano.AND.krep.EQ.reladc_fano_inrep) THEN
          WRITE(iw,*)
          WRITE(iw,*) '**************************************'
          WRITE(iw,*) '*                                    *'
          WRITE(iw,*) '* You selected a FanoADC run         *'
          WRITE(iw,*) '*                                    *'
          WRITE(iw,*) '**************************************'

          CALL MAKEEKVOO(EAJL,EPS,KREP,ICRA)
          WRITE(iw,*) 'CARITH = ', CARITH
          IF (CARITH) THEN
            CALL QUIT ('FanoADC does not run with complex algebra yet')
            IF (ADCLEVEL.EQ. 1) THEN
              CALL fanoadcc(CKKS,BUFY1,LADC,KREP,ITAPADC,IW,
     &                     EAJL,RDUMMY,RDUMMY,MXNO
     &                     )
            ELSE IF (ADCLEVEL.EQ. 2 .OR. ADCLEVEL.EQ. 3) THEN
              CALL fanoadcc(CKKS,BUFY1,LADC,KREP,ITAPADC,IW,
     &                     EAJL,BUFO1,BUF9A,MXNO
     &                     )
            END IF

          ELSE
            IF (ADCLEVEL.EQ. 1) THEN
              CALL fanoadcr(CKKS,BUFY1,LADC,KREP,ITAPADC,IW,
     &                     EAJL,RDUMMY,RDUMMY,MXNO
     &                     )
            ELSE IF (ADCLEVEL.EQ. 2 .OR. ADCLEVEL.EQ. 3) THEN
              CALL fanoadcr(CKKS,BUFY1,LADC,KREP,ITAPADC,IW,
     &                     EAJL,BUFO1,BUF9A,MXNO
     &                     )
            END IF
          END IF
          !CALL MEMREL('Fano-memrel',CC,KFRSAV,KFRDYN,KFREE,LFREE)
        END IF
!
! Do normal ADC calculation
!
        IF (reladc_fanoonly.eqv..false.) THEN
          WRITE(iw,*)
          WRITE(iw,*) '**************************************'
          WRITE(iw,*) '*                                    *'
          WRITE(iw,*) '* Start of normal ADC run            *'
          WRITE(iw,*) '* in symmetry: ', krep, '*'
          WRITE(iw,*) '*                                    *'
          WRITE(iw,*) '**************************************'
          WRITE(iw,*)
          WRITE(IW,*) 'ADC diagonal elements written to  ',FILEN4
          WRITE(IW,*) 'Packed ADC matrix written to      ',FILEN5
  
          OPEN(ITAPADC+4,FILE=FILEN4,FORM='UNFORMATTED',
     &         STATUS='UNKNOWN')
          OPEN(ITAPADC+5,FILE=FILEN5,FORM='UNFORMATTED',
     &         STATUS='UNKNOWN')

C
C  write the off-diagonal h/h 
C  and the h/2hp buffer to disk but NOT the diagonal entries since
C  there are more diagonal elements to be calculated in MAKEHPHP !
C  diagonal is written out in MAKEHPHP
C  all other elements are finished and can be written out.
C
          IF(CARITH) THEN
            CALL CWRITE_HHP(CKKS,BUFY1,ADCL,BUFIO,
     &               INTBUF,LADC,IOI,IOJ,
     &               KREP,NBUFS,IBUFPOS,ITAPADC+5)
          ELSE
            CALL WRITE_HHP(CKKS,BUFY1,ADCL,BUFIO,
     &            INTBUF,LADC,IOI,IOJ,
     &            KREP,NBUFS,IBUFPOS,ITAPADC+5)
          ENDIF

            WRITE(IW,*) 'Number of buffers after 1h/2h-1p:',NBUFS
            WRITE(IW,*) 'Buffer pointer after 1h/2h-1p:   ',IBUFPOS
C
C  ...   ************  start 2h-1p/2h-1p block  ********************
C
C  construct corresp. matrix entries. Arrays not addressed in ADC(2)
C  strict are referenced as dummies in order to avoid dangling
C  pointers. These arrays will not be referenced in the subroutine !
C  due to the out-of-core-character columns of the satellite block
C  and the diagonal are written out immediately in this routine
C
C
C  create K_akl array and lookup tables.
C
          CALL MAKEEKVOO(EAJL,EPS,KREP,ICRA)
c
c  compute the off-diagonal entries if required
!  BKC: C_aij,a'i'j' = |><|<|,<||>|>
c  ATT: the computed ones have complex conjugate BKC due to
c  the Schirmer,Trofimov formula. Therefore all calculated
c  entries will be CC is the routine, if complex case!
c
          IF(ADCLEVEL.EQ.1) THEN
            IF(CARITH) THEN
              CALL CMAKEHPHP(CKKS,DIAG,BUFIO,EAJL,CDUMMY,
     &                      CDUMMY,IOI,IOJ,ICRA,IDUMMY,LADC,
     &                      KREP,MXNO,NBUFS,IBUFPOS,INTBUF,ADCLEVEL,
     &                      ITAPADC+4,ITAPADC+5)
            ELSE
              CALL MAKEHPHP(CKKS,DIAG,BUFIO,EAJL,RDUMMY,
     &                      RDUMMY,IOI,IOJ,ICRA,IDUMMY,LADC,
     &                      KREP,MXNO,NBUFS,IBUFPOS,INTBUF,ADCLEVEL,
     &                      ITAPADC+4,ITAPADC+5)
            ENDIF
          ELSE IF(ADCLEVEL.EQ.2 .or. ADCLEVEL.EQ.3) THEN
            IF(CARITH) THEN
           CALL CMAKEHPHP(CKKS,DIAG,BUFIO,EAJL,BUFO1,
     &                    BUF9A,IOI,IOJ,ICRA,OOLT,
     &                    LADC,KREP,MXNO,NBUFS,IBUFPOS,INTBUF,ADCLEVEL,
     &                    ITAPADC+4,ITAPADC+5)
            ELSE
             CALL MAKEHPHP(CKKS,DIAG,BUFIO,EAJL,BUFO1,
     &                    BUF9A,IOI,IOJ,ICRA,OOLT,
     &                    LADC,KREP,MXNO,NBUFS,IBUFPOS,INTBUF,ADCLEVEL,
     &                    ITAPADC+4,ITAPADC+5)
            ENDIF
          ENDIF

C
C  ... close files in this symmetry and enter diagonalizer
C
        END IF
        CLOSE(ITAPADC+4)
        CLOSE(ITAPADC+5)


C
C  call configuration creator for 1h/2h1p configurations
C  File handle ITAPADC+2 is available. closed after sub call!
C
          WRITE(IW,*) 'ADC hole space configurations written to ',FILEN2
          WRITE(IW,*) 'Number of chunks in ADC matfile:         ',NBUFS
          CALL WCONDAT_S(ITAPADC+2,FILEN2,IRECL,KREP,LADC)
          WRITE(IW,*) 'Record length determined by WCONDAT_S:    ',IRECL
C
C -------------------------------------------------------------
C ---- from here ADC matrix construction has terminated   -----
C ---- and diagonalizers can do their work ...            -----
C ---- For each symmetry diagonalization is done          -----
C ---- immediately.                                       -----
C -------------------------------------------------------------
!
!  Elke: only do a diagonalization, when requested
!        this if construct will become obsolete as soon as the
!        call to the diagonalizer is moved to a separate subroutine
!
          IF (reladc_dodiag) THEN

c
c  Do a full diagonalization of the ADC matrix (for debugging
c  purposes) with eigenvalues and -vectors written to disk.
c  This is only allowed for matrices with a dimension smaller than
c  5000. Otherwise this routine is automatically exited.
c  This feature is especially useful when it comes to eigenvector
c  analysis. Due to simplicity FULLDIAR is wirtten in F90 and
c  retrievable in the separate module adcanal.F90 with dynamic
c  allocation of workspace.
c
            IF(DOFULL) THEN
              IF(CARITH) THEN
                CALL FULLDIAC(ITAPADC+4,FILEN4,FILEN5,INTBUF,LADC,KREP)
              ELSE
                CALL FULLDIAR(ITAPADC+4,FILEN4,FILEN5,INTBUF,LADC,KREP)
              ENDIF
            ENDIF
C
C  remember IRPLOOP counts through the input columns and is not equal to
C  the requested irrep. Therefore SIPEIGV also refers to the number of
C  eigenvectors requested for the KREP entered at the position IRPLOOP!
C
            reladc_md_iobase      = ITAPADC
            reladc_md_ionizl      = IONIZLEVEL
            reladc_md_ioldnew     = 1
            reladc_md_intbuf      = INTBUF
            reladc_md_desrep      = KREP
            reladc_md_rcw         = RCW
            reladc_md_lnzitr      = SIPITER
            reladc_md_matdim      = LADC
            reladc_md_irecl       = IRECL
            reladc_md_nmain       = NO(KREP)
            reladc_md_nbufs       = NBUFS
           reladc_md_eeigv_lower = reladc_md_sip_eeigv(IRPLOOP,1)/autoev
           reladc_md_eeigv_upper = reladc_md_sip_eeigv(IRPLOOP,2)/autoev
            reladc_md_fileadc     = FILEN5
            reladc_md_filediag    = FILEN4 
            reladc_md_filecnf     = FILEN2
            reladc_md_nmspec      = NMSPECF
c
c  call Lanczos diagonalizer. 
c
            CALL DIAG_LANC(IW,DOINCORE)
c
            WRITE(IW,*)
            WRITE(IW,*) '       ------------------------------'
            WRITE(IW,*) '       Finished SIP calculation'
            WRITE(IW,*) '       In symmetry',KREP
            WRITE(IW,*) '       ------------------------------'
            WRITE(IW,*)
            WRITE(IW,*)

          ELSE
            WRITE(IW,*)
            WRITE(IW,*) 'ADC matrix will not be diagonalized.'
            WRITE(IW,*)

          END IF 
!
! Not only Fano requested
!


 900  CONTINUE

      WRITE(IW,*) 'SIP ADC done for all selected symmetries.'
C
C                   ******************************
C                   ***      END USER LOOP     ***
C                   ***      DEALLOCATIONS     ***
C                   ******************************
!
!
      deallocate(CKKS, BUF1A, BUF1B, BUF4, BUF5A, BUF5B, EAJL, 
     &           ICRA, ADCL, DIAG, BUFIO, BUFY1, IOI, IOJ)
!
      IF(ADCLEVEL.GE.2) THEN
        deallocate(OOLT, BUFO1, BUF9A)
      ENDIF
!
      IF(ADCLEVEL.EQ.3) THEN
        deallocate(CKKSA, FKKS, BUF1C, BUF1D, BUF5C, BUF9B, BUF9C,
     &    BUFO2, BUFO3, BUFU1, BUFU2, BUFU3, BUFU4, BUFW1, BUFW2,
     &    BUFW3, BUFW4, BUFX1, EJAB, ADIAG1, ADIAG2, BDIAG1, BDIAG2,
     &    BUFQOO, BUFQVO, BUFQVOI, BUFQVOA, BUFQVV, BUFBOO, BUFBVO,
     &    BUFXVO, BUFYVO, BUFBVV, BUFXVV, BUFIVV, BUFBVO2, BUFSIG,
     &    CSIGM, ICRB, IPIVOT, VVLT)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE UAFFHO(UAH,
     &                  VVOO1,VVOO2,VVOO3,
     &                  VOVO1,VOVO2,VOVO3,
     &                  OOOO,EPS,KREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates coupling elements of affinity/hole range
C     Yields the following output:
C     U_jab,i  in this order and in |j><a|<b| , |i>  BKC
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 UAH(*),VVOO1(*),VVOO2(*),VVOO3(*)
      REAL*8 VOVO1(*),VOVO2(*),VOVO3(*)
      REAL*8 OOOO(*),EPS(*)
      INTEGER KREP
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 ddot,getimag
!     COMPLEX*16 zdotc
C---------------Executable Code --------------------------------------
C
C     CALL PST('Calculating coupling block AFFIN/HOLE+')
C
C  determine actual length of U-buffer in irrep KREP
C  and corresponding offset
C
      LENU = NO(KREP)*NOVVT(KREP)
      OFFS = LOOVVT(KREP)*RCW
C
C  ......................... part A
C
C  read VVOO integrals V_ab,ij and complex conjugate them.
C  Sort to V_jab,i (VVOO2) 
C  pick KREP and store as U_jab,i
C  After part A we have to reconjugate in order to establish
C  original integrals.
C
      CALL GETVVOO(VVOO1)
      CALL SRT22 (NREP,MULTB,LFA,NVVT,NO,NO,NOVVT,LOOVVT,LLOVVT,
     &            VVOO1,VVOO2)
      CALL XCOPY (LENU,VVOO2(OFFS+1),1,UAH,1)

c     if(carith) then
c       write(*,*) 'uaffho part A:',zdotc(lenu,uah,1,uah,1)
c     else
c       write(*,*) 'uaffho part A:',ddot(lenu,uah,1,uah,1)
c     endif
C
C  ......................... part B
C
C  form the v_ab,lm integrals from the V_ab,lm integrals still in VVOO1
C  and contract with V_lm,ij over lm.  ==> W_ab,ij (VVOO2)
C  Due to BKC the resulting W_ab,ij has to be conjugated.
C  sort W_ab,ij (VVOO2) to W_jab,i (VVOO3). Pick symmetry and
C  add **negative** on top of UAH (due to formula)
C
      CALL DENOMVVOO(EPS,VVOO1,VVOO1)
      CALL GETOOOO(OOOO)
      CALL CNTRCT('N','N',NVVT,NOOT,NOOT,A1,VVOO1,OOOO,A0,VVOO2,NREP)
      CALL SRT22 (NREP,MULTB,LFA,NVVT,NO,NO,NOVVT,LOOVVT,LLOVVT,
     &            VVOO2,VVOO3)
      CALL XAXPY (LENU,-A1,VVOO3(OFFS+1),1,UAH,1)

c     if(carith) then
c       write(*,*) 'uaffho part B:',zdotc(lenu,uah,1,uah,1)
c     else
c       write(*,*) 'uaffho part B:',ddot(lenu,uah,1,uah,1)
c     endif

C
C  ......................... part C
C
C  Sort the v_ac,jl (still in VVOO1) to v_aj,cl (VOVO1)
C  fetch the VOVO, conjugate and sort from
C  V_ci,bl (VOVO2) to V_cl,bi (VOVO3)
C  contract the v_aj,cl (VOVO1) and V_cl,bi (VOVO3) over _cl_
C  to A_aj,bi in VOVO2 --> A_a>b,ji via inverse SRT1TS4.
C
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,J2VOVO,JJVO,JJVO,VVOO1,VOVO1)
      CALL GETVOVO (VOVO2)
      IF (CARITH) CALL CONJUGA (NV4,VOVO2,1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            VOVO2,VOVO3)
C
C --- special contraction loop accounting for the different irrep
C     order in 1TT4 and SRT16 sorted arrays:
C
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,VOVO1(OFF1),M,VOVO3(OFF2),K,
     &               A0,VOVO2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
C  A_aj,bi  --> A_a>b,ji (VOVO2 --> VOVO3) only a/b antisymm !
C  A_a>b,ji --> A_ji,a>b (VOVO3 --> VOVO2)
C  A_ji,a>b --> A_ja>b,i (VOVO2 --> VOVO3)
C
C  The result is to be subtracted from the previous terms
C  according to the formulae. Att in the 1C1 sorter complex
C  conjugation is performed and establishes BKC !
C
      CALL SRT1TS4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,J2VOVO,
     &              JJVO,JJVO,VOVO2,VOVO3)
      CALL SRT1C1N (NREP,NVVT,NOO,VOVO3,VOVO2)
      CALL SRT6   (NREP,MULTB,LFA,NOO,NO,NO,NVVT,NOVVT,LOOVVT,LLOVVT,
     &             VOVO2,VOVO3)
      CALL XAXPY (LENU,-A1,VOVO3(OFFS+1),1,UAH,1)

c     if(carith) then
c       write(*,*) '** Coupling UAFFHO (c):',zdotc(lenu,uah,1,uah,1)
c     else
c       write(*,*) '** Coupling UAFFHO (r):',ddot(lenu,uah,1,uah,1)
c     endif
c     if(carith) write(*,*) 
c    &         '** UAFFHO imaginary part:',getimag(lenu,uah)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE UAFFPA(UAP,
     &                  VOVV1,VOVV2,
     &                  VVOO1,VVOO2,VOOO,VOVO,
     &                  EPS,KREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates coupling elements of affinity/particle range
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 UAP(*),VOVV1(*),VOVV2(*)
      REAL*8 VVOO1(*),VVOO2(*),VOOO(*),VOVO(*)
      REAL*8 EPS(*)
      INTEGER KREP
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      integer JVVVO_I48(MXREP+1)   !mxrep+1 is defined in symm.inc !
      LOGICAL DONE
      real*8 ddot,getimag
!     complex*16 zdotc
C
C---------------Executable Code --------------------------------------
C
c     CALL PST('Calculating coupling block AFFIN/PART+')
C
      LENU = NV(KREP)*NOVVT(KREP)
      OFFS = LVOVVT(KREP)*RCW
      DO I=1,NREP+1
        JVVVO_I48(I) = JVVVO(I)
      ENDDO

      IF(IVOVVT(NREP+1) .NE. LVOVVT(NREP+1)) THEN
        WRITE(IW,*) 'Unexpected size differences of arrays'
        CALL QUIT('Subroutine: UAFFPA')
      ENDIF
C
C   ..................... part A
C
C  read ALL VOVV integrals V_dj,ab (VOVV1) and sort to V_d,jab (VOVV2)
C  same number of integrals ! a>b
C  no cc due to BKC
C  then to V_jab,d (VOVV1) and pick corresponding IRREP for UAP
C
c     if(carith) then
c       write(*,*) 'vovv1 chksum (c)',zdotc(nv5,vovv1,1,vovv1,1)
c       write(*,*) 'RCW:',RCW
c     else
c       write(*,*) 'vovv1 chksum (r)',ddot(nv5,vovv1,1,vovv1,1)
c     endif

      NV5_I48 = NV5  ! transfer true I8 variable to the default length variable
      CALL RDVOVV(VOVV1)
      IF(CARITH) CALL CONJUGA(NV5_I48,VOVV1,1)
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NV,NO,NVVT,LVOVVT,LLOVVT,
     &             VOVV1,VOVV2)
      CALL SRT1C1N (NREP,NV,NOVVT,VOVV2,VOVV1)
      CALL XCOPY (LENU,VOVV1(OFFS+1),1,UAP,1)

c     if(carith) then
c       write(*,*) 'uaffpa1:',zdotc(lenu,uap,1,uap,1)
c     else
c       write(*,*) 'uaffpa1:',ddot(lenu,uap,1,uap,1)
c     endif
C
C   ..................... part B
C
C
C get VVOO integrals V_ab,lm and form v_ab,lm (VVOO1)
C and VOOO integrals V_dj,lm and form V*_dj,lm (VOOO)
C then contract over lm to A_dj,ab (VOVV1)
C then sort A_dj,ab (VOVV1) to A_d,jab and then to A_jab,d.
C add add the negative irrep contribution to UAP
C

      CALL GETVVOO (VVOO1)
      CALL DENOMVVOO(EPS,VVOO1,VVOO1)
      CALL GETVOOO (VOOO)
      IF(CARITH) CALL CONJUGA(NV2,VOOO,1)
      CALL CNTRCT ('N','T',NVO,NVVT,NOOT,A1,VOOO,VVOO1,A0,VOVV1,NREP)
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NV,NO,NVVT,LVOVVT,LLOVVT,
     &             VOVV1,VOVV2)
      CALL SRT1C1N (NREP,NV,NOVVT,VOVV2,VOVV1)
      CALL XAXPY (LENU,-A1,VOVV1(OFFS+1),1,UAP,1)

c     if(carith) then
c       write(*,*) 'uaffpa2:',zdotc(lenu,uap,1,uap,1)
c     else
c       write(*,*) 'uaffpa2:',ddot(lenu,uap,1,uap,1)
c     endif
C
C   ..................... part C
C
C  the v_ac,jl are still in VVOO1
C  and sort to v_aj,cl (VOVO)

      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,J2VOVO,JJVO,JJVO,VVOO1,VOVO)
C
C  get VOVV integrals V_bl,cd (VOVV1) and sort to V_bd,cl (VOVV2)
C  no additional complex conjugation ! the transposition is accounted
C  for in the 'T' parameter of the XGEMM routine.
C
      CALL RDVOVV(VOVV2)
      CALL SRT26 (NREP,MULTB,LFA,LFA,NV,NO,NV,NV,MVV,JVVVO_I48,
     &            JJVV,JJVO,VOVV2,VOVV1)
C
C  here again special contraction, with 'T' in XGEMM
C
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVV(IREP)
         K = MVO(JREP)
         OFF2 = J2VVVO(JREP) * RCW + 1
         CALL XGEMM ('N','T',M,N,K,A1,VOVO(OFF1),M,VOVV1(OFF2),N,
     &               A0,VOVV2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
C  next backsort and antisymmetrize  according to the 1XX4 order
C  A_aj,bd --> A_a>b,jd via new SRT1TS4 inverse sorter
C  Then sort A_a>b,jd  (<|<| |>|> structure) to
C  A_jd,a>b and then to A_jab,d via SRT6
C
      CALL SRT1TS4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NV,MVO,J2VOVV,
     &              JJVO,JJVV,VOVV2,VOVV1)
      CALL SRT1C1N (NREP,NVVT,NOV,VOVV1,VOVV2)
      CALL SRT6   (NREP,MULTB,LFA,NOV,NO,NV,NVVT,NOVVT,LVOVVT,LLOVVT,
     &             VOVV2,VOVV1)
      CALL XAXPY (LENU,-A1,VOVV1(OFFS+1),1,UAP,1)

c     if(carith) then
c       write(*,*) '** Coupling UAFFPA (c):',zdotc(lenu,uap,1,uap,1)
c     else
c       write(*,*) '** Coupling UAFFPA (r):',ddot(lenu,uap,1,uap,1)
c     endif
c     if(carith) write(*,*)
c    &    '** UAFFPA imaginary part:',getimag(lenu,uap)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE UIONHO(UIH,
     &                  VOOO1,VOOO2,VVOO,
     &                  VOVO,VOVV,
     &                  EPS,KREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates coupling elements of ionization/hole range
C     The BKC of this output is 
C     U_akl,j = |><|<|,|> and it is stored as U_akl,j
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 UIH(*)
      REAL*8 VOOO1(*),VOOO2(*),VVOO(*)
      REAL*8 VOVO(*),VOVV(*)
      REAL*8 EPS(*)
      INTEGER KREP

C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
      real*8 ddot,getimag
!     complex*16 zdotc
C
C---------------Executable Code --------------------------------------
C

c     CALL PST('Calculating coupling block IONIZ/HOLE+')

C  determine actual length of U-buffer in irrep KREP
C  and corresponding offset

      LENU = NO(KREP)*NVOOT(KREP)
      OFFS = LOVOOT(KREP)*RCW
C
C  ................. part A
C
C
C  get VOOO integrals V_ai,kl (VOOO1), negate/CC and sort to 
C  V_akl,i (VOOO2)
C
      CALL GETVOOO(VOOO1)
      CALL XSCAL(NV2,-A1,VOOO1,1)
      IF(CARITH) CALL CONJUGA(NV2,VOOO1,1)
      CALL SRT6(NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &          VOOO1,VOOO2)
      CALL XCOPY (LENU,VOOO2(OFFS+1),1,UIH,1)

c     if(carith) then
c       write(*,*) 'uionho1:',zdotc(lenu,uih,1,uih,1)
c     else
c       write(*,*) 'uionho1:',ddot(lenu,uih,1,uih,1)
c     endif
C
C  ................. part B
C
C  get the VVOO integrals V_bc.kl and form v*_bc,kl (VVOO)
C  get the VOVV integrals and form V*_ai,bc (VOVV)
C  contract over bc to A_ai,kl (VOOO1) and sort to A_akl,i (VOOO2)
C  
      CALL GETVVOO(VVOO)
      CALL DENOMVVOO(EPS,VVOO,VVOO)
      IF(CARITH) CALL CONJUGA(NV3,VVOO,1)

      NDIMV7=IVOVVT(NREP+1)
      CALL RDVOVV (VOVV)

      IF(CARITH) CALL CONJUGA(NDIMV7,VOVV,1)
      CALL CNTRCT ('N','N',NVO,NOOT,NVVT,A1,VOVV,VVOO,A0,VOOO1,NREP)
      CALL SRT6(NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &          VOOO1,VOOO2)
      CALL XAXPY (LENU,A1,VOOO2(OFFS+1),1,UIH,1)

c     if(carith) then
c       write(*,*) 'uionho2:',zdotc(lenu,uih,1,uih,1)
c     else
c       write(*,*) 'uionho2:',ddot(lenu,uih,1,uih,1)
c     endif
C
C  ................. part C (new nomenclature, see worked out formula !)
C
C  sort v^*_ab,km (still in VVOO) ---> B_ak,bm (VOVO)
C  get the V_bl,im (VOOO1) ----> C_bm,il (VOOO2)
C
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,J2VOVO,JJVO,JJVO,VVOO,VOVO)
      CALL GETVOOO(VOOO1)
      CALL SRT26 (NREP,MULTB,LFA,LFA,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &            VOOO1,VOOO2)
C
C  Perform special contraction loop v^*_ak,bm x V_bm,il ==> A_ak,il
C
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MOO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOOO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,VOVO(OFF1),M,VOOO2(OFF2),K,
     &               A0,VOOO1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
C  Reantisymmetrize A_ak,il --> A_ai,k>l
C  Complex conjugate in order to fulfill the BKC
C  and regroup to A_ai,k>l  ==> A_akl,i
C
      CALL SRT1ST4I(NREP,MULTB,LFA,LTR,NV,NO,NO,NO,MVO,J2VOOO,JJVO,JJOO,
     &              VOOO1,VOOO2)
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           VOOO2,VOOO1)
      CALL XAXPY (LENU,-A1,VOOO1(OFFS+1),1,UIH,1)

c     if(carith) then
c       write(*,*) '** Coupling UIONHO (c):',zdotc(lenu,uih,1,uih,1)
c     else
c       write(*,*) '** Coupling UIONHO (r):',ddot(lenu,uih,1,uih,1)
c     endif
c     if(carith) write(*,*) 
c    &     '** UIONHO imaginary part:',getimag(lenu,uih)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE UIONPA(UIP,
     &                  VVOO1,VVOO2,VVOO3,
     &                  VOVO1,VOVO2,VOVO3,
     &                  VVVV,EPS,KREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates coupling elements of ionization/particle range
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 UIP(*),VVOO1(*),VVOO2(*),VVOO3(*)
      REAL*8 VOVO1(*),VOVO2(*),VOVO3(*)
      REAL*8 VVVV(*),EPS(*)
      INTEGER KREP
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
      real*8 ddot,getimag
!     complex*16 zdotc
C
C---------------Executable Code --------------------------------------
C
c     CALL PST('Calculating coupling block IONIZ/PART+')
C
C  determine actual length of U-buffer in irrep KREP
C  and corresponding offset
C
      LENU = NV(KREP)*NVOOT(KREP)
      OFFS = LVVOOT(KREP)*RCW
C
C ********* part A
C
C
C  read VVOO integrals V_da,kl (VVOO1)
C  sort V_da,kl (VVOO1) to V*kl,da (VVOO2)
C  sort V*kl,da (VVOO2) to V*akl,d (VVOO3)
C  then pick irep KREP and store in UIP
C
      CALL GETVVOO(VVOO1)
      IF(CARITH) CALL CONJUGA(NV3,VVOO1,1)
      CALL SRT1C1N (NREP,NVVT,NOOT,VVOO1,VVOO2)
      CALL SRT22 (NREP,MULTB,LFA,NOOT,NV,NV,NVOOT,LVVOOT,LLVOOT,
     &            VVOO2,VVOO3)
      CALL XCOPY (LENU,VVOO3(OFFS+1),1,UIP,1)
c     if(carith) then
c       write(*,*) 'uionpa1:',zdotc(lenu,uip,1,uip,1)
c     else
c       write(*,*) 'uionpa1:',ddot(lenu,uip,1,uip,1)
c     endif
C
C ********* part B
C
C  form the -v^*_bc,kl integrals from the V_bc,kl (still in VVOO1)
C  negation due to energy denominator definition !
C  complex conjugation due to index pair swapping
C
      CALL GETVVOO(VVOO1)
      CALL DENOMVVOO(EPS,VVOO1,VVOO1)
      IF(CARITH) CALL CONJUGA(NV3,VVOO1,1)
      CALL XSCAL(NV3,-A1,VVOO1,1)
C
C  contract the -v^*_bc,kl (VVOO1) and V_bc,da (VVVV) to
C  A_kl,da (VVOO3)
C  then sort A_kl,da (VVOO3) to A_akl,d (VVOO2) and
C  ADD on top of UIP
C
      CALL RDVVVV(VVVV)
      CALL CNTRCT('T','N',NOOT,NVVT,NVVT,A1,VVOO1,VVVV,A0,VVOO3,NREP)
      CALL SRT22 (NREP,MULTB,LFA,NOOT,NV,NV,NVOOT,LVVOOT,LLVOOT,
     &            VVOO3,VVOO2)
      CALL XAXPY (LENU,A1,VVOO2(OFFS+1),1,UIP,1)

c     if(carith) then
c       write(*,*) 'uionpa2:',zdotc(lenu,uip,1,uip,1)
c     else
c       write(*,*) 'uionpa2:',ddot(lenu,uip,1,uip,1)
c     endif
C
C ********** part C
C
C  sort the v_ab,km (in VVOO1) to v_ak,bm (VOVO1)
C  get VOVO integrals V_bl,dm (VOVO2) and sort to V_bm,dl (VOVO3)
C  contract v_bm,ak and V_bm,dl to A_ak,dl and reantisym to A_ad,k>l.
C  then from A_ad,k>l to A_ak>l,d.
C
      CALL GETVVOO(VVOO1)
      CALL DENOMVVOO(EPS,VVOO1,VVOO1)
      IF(CARITH) CALL CONJUGA(NV3,VVOO1,1)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              VVOO1,VOVO1)
      CALL GETVOVO(VOVO2)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            VOVO2,VOVO3)
C
C  next contract via special loop taking into account the
C  varying bosonic irrep order in VOVO1 and VOVO3 due to the sorting
C  conditions
C
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,VOVO1(OFF1),M,VOVO3(OFF2),K,
     &               A0,VOVO2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT1ST4I(NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              VOVO2,VOVO3)
      CALL SRT6 (NREP,MULTB,LFA,NVV,NV,NV,NOOT,NVOOT,LVVOOT,LLVOOT,
     &           VOVO3,VOVO2)

      CALL XAXPY (LENU,-A1,VOVO2(OFFS+1),1,UIP,1)

c     if(carith) then
c       write(*,*) '** Coupling UIONPA (c):',zdotc(lenu,uip,1,uip,1)
c     else
c       write(*,*) '** Coupling UIONPA (r):',ddot(lenu,uip,1,uip,1)
c     endif
c     if(carith) write(*,*) 
c    6        '** UIONPA imaginary part:',getimag(lenu,uip)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE VAFFITER (DIAG1,DIAG2,VVVV,VOVO,
     &                     UAH,V0,W1,W2,
     &                     EKJAB,EPS,ICRB,VVLT,MXNV,KREP,VCONV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     This subroutine determines V+_k, k occupied
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER MXNV,ICRB(6,*),VVLT(MXNV,MXNV,*),KREP
      REAL*8 DIAG1(*),DIAG2(*),VVVV(*),VOVO(*)
      REAL*8 UAH(*),V0(*),W1(*),W2(*)
      REAL*8 EKJAB(*),EPS(*),VCONV
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER ITERC
      REAL*8 KCVIR,S,SD,EPSK
      REAL*8 SEC0,SEC1,SEC2,SECA,TITERTIME

      REAL*8 DDOT,ECC
C
C---------------Executable code--------------------------------------
C
c     CALL PST('computing V affinity block (real)+')
      WRITE(IW,710) KREP,VCONV
C
C ... set dimensions
C
      LENK = NOVVT(KREP)   ! length of one vector in V_k
      LENV = NO(KREP)*NOVVT(KREP) ! size of V_k array in total
      KOFF = IO(KREP) ! offset for occupied energy
C
C ... clear buffers
C
      CALL XCOPY(LENV,A0,0,V0,1)
      CALL XCOPY(LENV,A0,0,W1,1)
C
C ... read VOVO and VVVV integrals
C
      CALL GETVOVO(VOVO)
      CALL RDVVVV(VVVV)
C
C  construct VVT lookup table for all IRREPS
C
      DO IREP=1,NREP
      DO I=1,MXNV
      DO J=1,MXNV
        VVLT(J,I,IREP)=0
      ENDDO
      ENDDO
      ENDDO
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NV(IREP)
          DO I=J+1,NV(IREP)
            VVLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
C
C ... construct diag( K + C ) affinity block (length of a vector in the
C     V array)
C
      DO I=1,LENK
        DIAG1(I) = KCVIR(I,I,VVVV,VOVO,ICRB,VVLT,MXNV,EPS,ecc)
        DIAG1(I) = DIAG1(I) + EKJAB(I)
c       write(*,*) 'vaffiter diagelement:',i,diag1(i)
      ENDDO

C
C ... generate the V_k^(0) start vectors via
C     V_k^(0) = D_k^(-1) * U_k
C     for each k of the occupied orbitals
C     the specific diagonal D_k has to be created !
C
      OFF = 0
      DO KOCC = 1,NO(KREP)
        EPSK = EPS(KOFF+KOCC)
        WRITE(IW,720) KOCC,EPSK
        DO I=1,LENK
          DIAG2(I) = 1.0D0/(EPSK - DIAG1(I))
        ENDDO
        DO I=1,LENK
          V0(OFF+I)=DIAG2(I)*UAH(OFF+I)
        ENDDO
        OFF = OFF + LENK
      ENDDO
C
C     fill V0 with all k start vectors.
C     att! V0 is still needed is the iteration steps !
C
      CALL XCOPY(LENV,V0,1,W1,1)
C
C ... enter iterations
C ... vector V_k^(n) have to be cleared in each iteration step
C
      WRITE(IW,*)
      ITERC=0
      CALL CPUUSED(SEC0)
      SEC1=SEC0

 111  CONTINUE

      ITERC = ITERC + 1
      CALL XCOPY(LENV,A0,0,W2,1)

C ...      form A' * V_k^(n-1) over ALL k vectors


      DO J = 1,LENK
        DO I = J+1,LENK
          S = KCVIR(I,J,VVVV,VOVO,ICRB,VVLT,MXNV,EPS,ECC)
          OFF = 0
          DO K = 1,NO(KREP)
            W2(OFF+I) = W2(OFF+I) + S*W1(OFF+J)
            W2(OFF+J) = W2(OFF+J) + S*W1(OFF+I)
            OFF = OFF + LENK
          ENDDO
        ENDDO
      ENDDO
C
C  multiply all k vectors in W2 now with the D_k^(-1)
C  that means form D_k^(-1) * A' * V_k^(n-1)
C

      OFF = 0
      DO KOCC = 1,NO(KREP)
        EPSK = EPS(KOFF+KOCC)
        DO I=1,LENK
          DIAG2(I) = 1.0D0/(EPSK - DIAG1(I))
        ENDDO
        DO I=1,LENK
          W2(OFF+I)=DIAG2(I)*W2(OFF+I)
        ENDDO
        OFF = OFF + LENK
      ENDDO
C
C  finally add the zeroth order vectors V0 to W2
C
      CALL XAXPY(LENV,A1,V0,1,W2,1)

C  and form euklidean distance between V_k^n and V_k^(n-1)
C  over all k-vectors

      S=0.0D0 
      DO I=1,LENV
        SD = W2(I) - W1(I)
        S = S + SD*SD
      ENDDO
#if defined (BIT64)
      S=SQRT(S)
#else
      S=DSQRT(S)
#endif
      CALL CPUUSED(SEC2)
      WRITE(IW,730) ITERC,S,SEC2-SEC1
      SEC1=SEC2
      IF(S.GT.VCONV) THEN
C ... not converged yet, reiterate
        CALL XCOPY(LENV,W2,1,W1,1)
        GOTO 111
      ENDIF
c     write(*,*) 'vaffiter chksum',ddot(lenv,w2,1,w2,1)
C
C ... iteration converged, V+ vectors now available in W2
C
      TITERTIME = SEC2-SEC0
      SECA = TITERTIME/DBLE(ITERC)
c     WRITE(IW,740) TITERTIME,SECA
      WRITE(IW,*) ' Convergence reached.'

 710  FORMAT (/1X,'---',39X,'---'
     &        /1X,'--- Inverse Iteration for the hole states ---',
     &        /1X,'---',39X,'---'/1X,'Symmetry:',T10,I2,
     &        /1X,'Threshold (VCONV): ',T20,G15.4//)
 720  FORMAT (5X,'#',I3,' with spinor energy  ',G15.8)
 730  FORMAT (5X,'Iteration  #',I3,' convergence  ',G15.8,
     &        ' Time',F10.3)
 740  FORMAT (//8X,'*** Total time:',F10.3,
     &        ' Average:',F10.3,' ***'//)


      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE VIONITER (DIAG1,DIAG2,OOOO,VOVO,
     &                     UIP,V0,W1,W2,
     &                     EKAKL,EPS,ICRA,OOLT,MXNO,KREP,VCONV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     This subroutine determines V-_k, k virtual
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER MXNO,ICRA(6,*),OOLT(MXNO,MXNO,*),KREP
      REAL*8 DIAG1(*),DIAG2(*),OOOO(*),VOVO(*)
      REAL*8 UIP(*),V0(*),W1(*),W2(*),EKAKL(*),EPS(*),VCONV
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER ITERC
      REAL*8 KCOCC,S,SD,EPSK

      real*8 ddot
C
C---------------Executable code--------------------------------------
C
c     CALL PST('computing V ionization block (real)+')
      WRITE(IW,710) KREP,VCONV
C
C ... clear buffers
C
      NN=NV(KREP)*NVOOT(KREP)
      CALL XCOPY(NN,A0,0,V0,1)
      CALL XCOPY(NN,A0,0,W1,1)
C
C ... get OOOO/VOVO integrals
C
      CALL GETOOOO(OOOO)
      CALL GETVOVO(VOVO)
C
C  construct OOT lookup table for the IRREPS
C
      DO IREP=1,NREP
      DO I=1,MXNO
      DO J=1,MXNO
        OOLT(J,I,IREP)=0
      ENDDO
      ENDDO
      ENDDO
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NO(IREP)
          DO I=J+1,NO(IREP)
            OOLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
C
C ... construct diag( K + C ) ionization block
C
      DO I=1,NVOOT(KREP)
        DIAG1(I) = KCOCC(I,I,OOOO,VOVO,ICRA,OOLT,MXNO)
        DIAG1(I) = DIAG1(I) + EKAKL(I)
      ENDDO
C
C   determine length and energy offset
C
      IMAX = NVOOT(KREP)
      KOFF = IO(NREP+1) + IV(KREP)
C
C ... generate the V_k^(0) start vectors via
C     V_k^(0) = D_k^(-1) * U_k
C     for each k of the virtual orbitals
C     the specific diagonal has to be created !
C
      OFF = 0
      DO KVIR = 1,NV(KREP)
        EPSK = EPS(KOFF+KVIR)
        WRITE(IW,720) KVIR,EPSK
        DO I=1,IMAX
          DIAG2(I) = 1.0D0/(EPSK - DIAG1(I))
        ENDDO
        DO I=1,IMAX
          V0(OFF+I)=DIAG2(I)*UIP(OFF+I)
        ENDDO
        OFF = OFF + IMAX
      ENDDO
c     write(iw,*) ' **** V_ion V_0 start vectors',
c    & ddot(nn,v0,1,v0,1)
C
C     and fill iteration start vector with V0
C     att! V0 is still needed is the iteration steps !
C
      CALL XCOPY(NN,V0,1,W1,1)
C
C ... enter iterations
C
      ITERC=0

 111  CONTINUE

      ITERC = ITERC + 1
      CALL XCOPY(NN,A0,0,W2,1)

C ...      form A' * V_k^(n-1) over ALL k vectors

      DO J = 1,IMAX
        DO I = J+1,IMAX
          S = KCOCC(I,J,OOOO,VOVO,ICRA,OOLT,MXNO)
          OFF = 0
          DO K = 1,NV(KREP)
            W2(OFF+I) = W2(OFF+I) + S*W1(OFF+J)
            W2(OFF+J) = W2(OFF+J) + S*W1(OFF+I)
            OFF = OFF + IMAX
          ENDDO
        ENDDO
      ENDDO

C ...      multiply all k vectors in W2 now with the D_k^(-1)
C ...      that means form D_k^(-1) * A' * V_k^(n-1)

      OFF = 0
      DO KVIR = 1,NV(KREP)
        EPSK = EPS(KOFF+KVIR)
        DO I=1,IMAX
          DIAG2(I) = 1.0D0/(EPSK - DIAG1(I))
        ENDDO
        DO I=1,IMAX
          W2(OFF+I)=DIAG2(I)*W2(OFF+I)
        ENDDO
        OFF = OFF + IMAX
      ENDDO

C ...      finally add the zeroth order vectors V0 to W2

      CALL XAXPY(NN,A1,V0,1,W2,1)

C ...      and form a euklidean distance between V_k^n and V_k^(n-1)
C ...      over all k-vectors

      S=0.0D0 
      DO I=1,NN
        SD = W2(I) - W1(I)
        S = S + SD*SD
      ENDDO
#if defined (BIT64)
      S=SQRT(S)
#else
      S=DSQRT(S)
#endif
      WRITE(IW,730) ITERC,S
      IF(S.GT.VCONV) THEN
C ... not converged yet, reiterate
        CALL XCOPY(NN,W2,1,W1,1)
        GOTO 111
      ENDIF
      WRITE(IW,*) 'Convergence reached.'
c     write(*,*) 'vioniter chksum',ddot(nn,w2,1,w2,1)

C ... iteration converged, V- vectors now available in W2

 710  FORMAT (/1X,'---',43X,'---'
     &        /1X,'--- Inverse Iteration for the particle states ---',
     &        /1X,'---',43X,'---'/1X,'Symmetry:',T30,I2,
     &        /1X,'Threshold (VCONV): ',T20,G15.4//)
 720  FORMAT (5X,'#',I3,' with spinor energy  ',G15.8)
 730  FORMAT (5X,'Iteration  #',I3,' convergence  ',G15.8)
        
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CVAFFITER (DIAG1,DIAG2,VVVV,VOVO,
     &                      UAH,V0,W1,W2,
     &                      EKJAB,EPS,ICRB,VVLT,MXNV,KREP,VCONV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     This subroutine determines V+_k, k occupied in the complex case
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 DIAG1(*),DIAG2(*),VVVV(*),VOVO(*)
      COMPLEX*16 UAH(*),V0(*),W1(*),W2(*)
      REAL*8 EKJAB(*),EPS(*)
      INTEGER ICRB(6,*),VVLT(MXNV,MXNV,*),MXNV,KREP
      REAL*8 VCONV
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER ITERC
      COMPLEX*16 CEPSK,CKCVIR,DII,S_UNST,S_STAR,S
      REAL*8 SEC0,SEC1,SEC2,SECA,TITERTIME

!     COMPLEX*16 ZDOTC
      REAL*8 GETIMAG
C
C---------------Executable code--------------------------------------
C
      CALL PST('computing V affinity block (complex)+')
      WRITE(IW,710) KREP,VCONV
C
C ... set dimensions
C
      LENK = NOVVT(KREP)   ! length of one vector in V_k
      LENV = NO(KREP)*NOVVT(KREP) ! size of V_k array in total
      KOFF = IO(KREP) ! offset for occupied energy
C
C ... clear buffers and diagonals
C
      CALL XCOPY(LENV,A0,0,V0,1)
      CALL XCOPY(LENV,A0,0,W1,1)
      CALL XCOPY(LENK,A0,0,DIAG1,1)
      CALL XCOPY(LENK,A0,0,DIAG2,1)
C
C ... read VOVO and VVVV integrals
C
      CALL GETVOVO(VOVO)
      CALL RDVVVV(VVVV)
C
C  construct VVT lookup table for all IRREPS
C
      DO IREP=1,NREP
      DO I=1,MXNV
      DO J=1,MXNV
        VVLT(J,I,IREP)=0
      ENDDO
      ENDDO
      ENDDO
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NV(IREP)
          DO I=J+1,NV(IREP)
            VVLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
C
C ... construct diag( K + C ) affinity block (length of a vector in the
C     V array)
C
      DO I=1,LENK
        DII = CKCVIR(I,I,VVVV,VOVO,ICRB,VVLT,MXNV)
        IF(DIMAG(DII).GT.1.0D-10) THEN
          WRITE(IW,*) ' *** WARNING ***   Diagelement ',I,
     &                ' has imaginary part of size',DIMAG(DII)
        ENDIF
        DIAG1(I) = DII + DCMPLX(EKJAB(I),0.0D0)
      ENDDO
C
C ... generate the V_k^(0) start vectors via
C     V_k^(0) = D_k^(-1) * U_k
C     for each k of the occupied orbitals
C     the specific diagonal D_k has to be created !
C
      OFF = 0
      DO KOCC = 1,NO(KREP)
        CEPSK = DCMPLX(EPS(KOFF+KOCC),0.0D0)
        WRITE(IW,720) KOCC,DBLE(CEPSK)
        DO I=1,LENK
          DIAG2(I) = (1.0D0,0.0D0)/(CEPSK - DIAG1(I))
        ENDDO
        DO I=1,LENK
          V0(OFF+I) = DIAG2(I) * UAH(OFF+I)
        ENDDO
        OFF = OFF + LENK
      ENDDO
C
C     fill W1 with all k start vectors.
C     att! V0 is still needed in the iteration steps !
C
      CALL XCOPY(LENV,V0,1,W1,1)
C
C ... enter iterations
C ... vector V_k^(n) have to be cleared in each iteration step
C
      WRITE(IW,*)
      ITERC=0
      CALL CPUUSED(SEC0)
      SEC1=SEC0

 111  CONTINUE

      ITERC = ITERC + 1
      CALL XCOPY(LENV,A0,0,W2,1)

C ... form A' * V_k^(n-1) over ALL k vectors
C           *** ATTENTION ***
C ... please note that the upper triangular part is the
C ... complex conjugate of the lower triangular part
C ... and a picked element has to be complex conjugated 
C ... for the upper triangular contribution
C
      DO J = 1,LENK
        DO I = J+1,LENK

          S_UNST = DCONJG(CKCVIR(I,J,VVVV,VOVO,ICRB,VVLT,MXNV))
          S_STAR = DCONJG(S_UNST)

          OFF = 0
          DO K = 1,NO(KREP)
            W2(OFF+I) = W2(OFF+I) + S_UNST*W1(OFF+J)
            W2(OFF+J) = W2(OFF+J) + S_STAR*W1(OFF+I)
            OFF = OFF + LENK
          ENDDO
        ENDDO
      ENDDO
C
C  multiply all k vectors in W2 now with the D_k^(-1)
C  that means form D_k^(-1) * A' * V_k^(n-1)
C

      OFF = 0
      DO KOCC = 1,NO(KREP)
        CEPSK = DCMPLX(EPS(KOFF+KOCC),0.0D0)
        DO I=1,LENK
          DIAG2(I) = (1.0D0,0.0D0)/(CEPSK - DIAG1(I))
        ENDDO
        DO I=1,LENK
          W2(OFF+I) = DIAG2(I) * W2(OFF+I)
        ENDDO
        OFF = OFF + LENK
      ENDDO
C
C  finally add the zeroth order vectors V0 to W2
C
      CALL XAXPY(LENV,A1,V0,1,W2,1)

C  and form euklidean distance between V_k^n and V_k^(n-1)
C  over all k-vectors

      S=A0
      DO I=1,LENV
        DII = W2(I) - W1(I)
        S = S + DCONJG(DII)*DII
      ENDDO
#if defined (BIT64)
      S=SQRT(S)
#else
      S=SQRT(S)
#endif
      CALL CPUUSED(SEC2)
      WRITE(IW,730) ITERC,ABS(S),SEC2-SEC1
      SEC1=SEC2
      IF(ABS(S).GT.VCONV) THEN
C ... not converged yet, reiterate
        CALL XCOPY(LENV,W2,1,W1,1)
        GOTO 111
      ENDIF
c     write(*,*) 'cvaffiter chksum',zdotc(lenv,w2,1,w2,1)
C
C ... iteration converged, V+ vectors now available in W2
C
      TITERTIME = SEC2-SEC0
      SECA = TITERTIME/DBLE(ITERC)
c     WRITE(IW,740) TITERTIME,SECA
      WRITE(IW,*) ' Convergence reached.'

 710  FORMAT (/10X,'---',44X,'---'
     &        /10X,'--- Complex Inv. Iteration for the hole states ---',
     &        /10X,'---',44X,'---'//10X,'Symmetry:',T30,I2,
     &        /10X,'Threshold: ',T30,G15.4//)
 720  FORMAT (5X,'#',I3,' with spinor energy  ',G15.8)
 730  FORMAT (5X,'Iteration  #',I3,' convergence  ',G15.8,
     &        ' Time',F10.3)
 740  FORMAT (//8X,'*** Total AFF iteration time:',F10.3,
     &        ' Average:',F10.3,' ***'//)


      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CVIONITER (DIAG1,DIAG2,OOOO,VOVO,
     &                      UIP,V0,W1,W2,
     &                      EKAKL,EPS,ICRA,OOLT,MXNO,KREP,VCONV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     This subroutine determines V-_k, k virtual for the complex case
C     Due to hermiticity the diagonal entries have to be real!
C     The calling routine therefore provides only arrays of
C     corresponding real length !
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 DIAG1(*),DIAG2(*),OOOO(*),VOVO(*)
      COMPLEX*16 UIP(*),V0(*),W1(*),W2(*)
      INTEGER MXNO,ICRA(6,*),OOLT(MXNO,MXNO,*),KREP
      REAL*8 EKAKL(*),EPS(*),VCONV
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER ITERC
      COMPLEX*16 CKCOCC,DII,CEPSK
      COMPLEX*16 S_UNST,S_STAR
      REAL*8 S

!     COMPLEX*16 ZDOTC
      REAL*8 GETIMAG

C
C---------------Executable code--------------------------------------
C
      CALL PST('computing V ionization block (complex)+')
      WRITE(IW,710) KREP,VCONV
C
C ... clear buffers
C
      LENK=NVOOT(KREP)
      LENV=NV(KREP)*NVOOT(KREP)

      CALL XCOPY(LENV,A0,0,V0,1)
      CALL XCOPY(LENV,A0,0,W1,1)
      CALL XCOPY(LENK,A0,0,DIAG1,1)
      CALL XCOPY(LENK,A0,0,DIAG2,1)
C
C ... get OOOO/VOVO integrals
C
      CALL GETOOOO(OOOO)
      CALL GETVOVO(VOVO)
C
C  construct OOT lookup table for the IRREPS
C
      DO IREP=1,NREP
      DO I=1,MXNO
      DO J=1,MXNO
        OOLT(J,I,IREP)=0
      ENDDO
      ENDDO
      ENDDO
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NO(IREP)
          DO I=J+1,NO(IREP)
            OOLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
C
C ... construct diag( K + C ) ionization block
C
      DO I=1,LENK
        DII = CKCOCC(I,I,OOOO,VOVO,ICRA,OOLT,MXNO)
        IF(DIMAG(DII).GT.1.0D-10) THEN
          WRITE(IW,*) ' W *** Diagelement ',I,' has imaginary part',DII
        ENDIF
        DIAG1(I) = DII + DCMPLX(EKAKL(I),0.0D0)
      ENDDO
C
C   determine length and energy offset
C
      KOFF = IO(NREP+1) + IV(KREP)
C
C ... generate the V_k^(0) start vectors via
C     V_k^(0) = D_k^(-1) * U_k
C     for each k of the virtual orbitals
C     the specific diagonal has to be created !
C
      OFF = 0
      DO KVIR = 1,NV(KREP)
        CEPSK = DCMPLX(EPS(KOFF+KVIR),0.0D0)
        WRITE(IW,720) KVIR,DBLE(CEPSK)
        DO I=1,LENK
          DIAG2(I) = (1.0D0,0.0D0)/(CEPSK - DIAG1(I))
        ENDDO
        DO I=1,LENK
          V0(OFF+I)=DIAG2(I)*UIP(OFF+I)
        ENDDO
        OFF = OFF + LENK
      ENDDO
C
C     and fill iteration start vector with V0
C     att! V0 is still needed is the iteration steps !
C
      CALL XCOPY(LENV,V0,1,W1,1)
C
C ... enter iterations
C
      WRITE(IW,*)
      ITERC=0

 111  CONTINUE

      ITERC = ITERC + 1
      CALL XCOPY(LENV,A0,0,W2,1)

C ... form A' * V_k^(n-1) over ALL k vectors
C ...    *** ATTENTION ***
C ... see comment in CVAFFITER for this multiplication !!


      DO J = 1,LENK
        DO I = J+1,LENK
          S_UNST = CKCOCC(I,J,OOOO,VOVO,ICRA,OOLT,MXNO)
          S_STAR = DCONJG(S_UNST)
          OFF = 0
          DO K = 1,NV(KREP)
            W2(OFF+I) = W2(OFF+I) + S_UNST*W1(OFF+J)
            W2(OFF+J) = W2(OFF+J) + S_STAR*W1(OFF+I)
            OFF = OFF + LENK
          ENDDO
        ENDDO
      ENDDO
   

C ...      multiply all k vectors in W2 now with the D_k^(-1)
C ...      that means form D_k^(-1) * A' * V_k^(n-1)

      OFF = 0
      DO KVIR = 1,NV(KREP)
        CEPSK = DCMPLX(EPS(KOFF+KVIR),0.0D0)
        DO I=1,LENK
          DIAG2(I) = (1.0D0,0.0D0)/(CEPSK - DIAG1(I))
        ENDDO
        DO I=1,LENK
          W2(OFF+I)=DIAG2(I)*W2(OFF+I)
        ENDDO
        OFF = OFF + LENK
      ENDDO

C ...      finally add the zeroth order vectors V0 to W2

      CALL XAXPY(LENV,A1,V0,1,W2,1)

C ...      and form a euklidean distance between V_k^n and V_k^(n-1)
C ...      over all k-vectors

      S=0.0D0 
      DO I=1,LENV
        DII = W2(I) - W1(I)
        S = S + DCONJG(DII)*DII
      ENDDO
#if defined (BIT64)
      S=SQRT(S)
#else
      S=DSQRT(S)
#endif
      WRITE(IW,730) ITERC,S
      IF(S.GT.VCONV) THEN
C ... not converged yet, reiterate
        CALL XCOPY(LENV,W2,1,W1,1)
        GOTO 111
      ENDIF
      WRITE(IW,*) 'Convergence reached.'
c     write(*,*) 'cvioniter chksum',zdotc(LENV,w2,1,w2,1)

C ... iteration converged, V- vectors now available in W2

 710  FORMAT (/10X,'---',48X,'---'
     &    /10X,'--- Complex Inv. Iteration for the particle states ---',
     &        /10X,'---',48X,'---'//10X,'Symmetry:',T30,I2,
     &        /10X,'Threshold: ',T30,G15.4//)
 720  FORMAT (5X,'#',I3,' with spinor energy  ',G15.8)
 730  FORMAT (5X,'Iteration  #',I3,' convergence  ',G15.8)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RMAKE_QKL(UIH,UAP,VIP,VAH,
     &                    QOO,QVO,QVOI,QVOA,
     &                    QVV,EPS,KREP,ADCPRINT)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates the Q_kl values from the corresponding
C     U coupling blocks and iterated V vectors
C     hereby the Q_oo, Q_vo, ... individual arrays are generated
C     as sums of Q+ and Q- already in one step. it is done for the
C     specified symmetry KREP such that Q_kl is totally symmetric
C     when l=krep !
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 UIH(*),UAP(*),VIP(*),VAH(*)
      REAL*8 QOO(*),QVO(*),QVOI(*),QVOA(*)
      REAL*8 QVV(*),EPS(*)
      INTEGER KREP,ADCPRINT(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 E,E1,E2,F
      REAL*8 DDOT
C
C---------------Executable Code --------------------------------------
C
c     CALL PST('Generating real Q_kl+')
c     WRITE(IW,*) ' **** Symmetry',KREP

C
C  LENA is length of affinity blocks (U or V)
C  LENI is length of ionization blocks (U or V)
C
      LENA = NOVVT(KREP)
      LENI = NVOOT(KREP)
C
C calculate Q+_oo + Q-_oo (zero for Q-) for a given symmetry KREP such
C that Q_kl contains only totally symmetric contributions.
C
      OFF1 = JJOO(KREP,KREP) + 1
      DO JX=1,NO(KREP)
        OFFJJ = (JX-1)*LENA + 1
        DO IX=1,NO(KREP)
          OFFII = (IX-1)*LENA + 1
          E=DDOT(LENA,VAH(OFFII),1,VAH(OFFJJ),1)
          QOO(OFF1) = -E
          OFF1 = OFF1 + 1
        ENDDO
      ENDDO
C
C calculate Q+_vo + Q-_vo for a given symmetry KREP such that
C Q_kl contains only totally symmetric contributions.
C

      OFF1 = JJVO(KREP,KREP) + 1
      DO JX=1,NO(KREP)
        JXE = IO(KREP) + JX
        OFFVAFF_J=(JX-1)*LENA + 1
        OFFUION_J=(JX-1)*LENI + 1
        DO AX=1,NV(KREP)
          AXE = IO(NREP+1) + IV(KREP) + AX
          OFFUAFF_A=(AX-1)*LENA + 1
          OFFVION_A=(AX-1)*LENI + 1
          E1=DDOT(LENA,UAP(OFFUAFF_A),1,VAH(OFFVAFF_J),1)
          E2=DDOT(LENI,VIP(OFFVION_A),1,UIH(OFFUION_J),1)
          F=EPS(AXE) - EPS(JXE)
          F=-1.0D0/F
          QVOA(OFF1) = F*(E1)
          QVOI(OFF1) = F*(E2)
          OFF1 = OFF1 + 1
        ENDDO
      ENDDO
C
C  add +/- contributions of QVO in corresponding symmetry
C
      OFF1 = JJVO(KREP,KREP) + 1
      DO JX=1,NO(KREP)
        DO AX=1,NV(KREP)
          QVO(OFF1) = QVOI(OFF1) + QVOA(OFF1)
          OFF1 = OFF1 + 1
        ENDDO
      ENDDO

c  debug
c     nx=no(krep)*nv(krep)
c     off1=jjvo(krep,krep) + 1
c     WRITE(IW,*) 'QVO chksum:',ddot(nx,qvo(off1),1,qvo(off1),1)
c  debug end

C
C calculate Q+_vv + Q-_vv (zero for Q+) for a given symmetry KREP such
C that Q_kl contains only totally symmetric contributions.
C
      OFF1 = JJVV(KREP,KREP) + 1
      DO BX=1,NV(KREP)
        OFFBB = (BX-1)*LENI + 1
        DO AX=1,NV(KREP)
          OFFAA = (AX-1)*LENI + 1
          IF(CARITH) THEN
            STOP 'Should not occur here (rmake_qkl)'
          ELSE
            E=DDOT(LENI,VIP(OFFAA),1,VIP(OFFBB),1)
            QVV(OFF1) = E
          ENDIF
          OFF1 = OFF1 + 1
        ENDDO
      ENDDO

c  debug
c     nx=nv(krep)*nv(krep)
c     off1=jjvv(krep,krep) + 1
c     WRITE(IW,*) 'QVV chksum:',ddot(nx,qvv(off1),1,qvv(off1),1)
c  debug end

C
C  write out Q_kl matrices dependent on PRINT level
C  (will be implemented later)
C  ionization/affinity part is output separately.
C  For the subsequent treatment both contributions have to be added !
C
c     IF(ADCPRINT(5).NE.0) THEN
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ Q_kl (real) in Symmetry:',KREP
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ h/h part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL MATOUT(QOO(JJOO(KREP,KREP)+1),NO(KREP),NO(KREP))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ p/h part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL MATOUT(QVO(JJVO(KREP,KREP)+1),NV(KREP),NO(KREP))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ p/p part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL MATOUT(QVV(JJVV(KREP,KREP)+1),NV(KREP),NV(KREP))
c     ENDIF
      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMAKE_QKL(UIH,UAP,VIP,VAH,
     &                    QOO,QVO,QVOI,QVOA,
     &                    QVV,EPS,KREP,ADCPRINT)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates the Q_kl values from the corresponding
C     complex U coupling blocks and iterated complex V vectors
C     hereby the Q_oo, Q_vo, ... individual arrays are generated
C     as sums of Q+ and Q- already in one step. it is done for the
C     specified symmetry KREP such that Q_kl is totally symmetric
C     when l=krep !
C     important for the complex case:
C     The BKC for the Q_kl is |><| and not <||> !!
C     One sees this when looking at the formula for the b_pq which have
C     the BKC <||> !!
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 UIH(*),UAP(*),VIP(*),VAH(*)
      COMPLEX*16 QOO(*),QVO(*),QVOI(*),QVOA(*),QVV(*)
      REAL*8 EPS(*)
      INTEGER KREP,ADCPRINT(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 E,E1,E2
      COMPLEX*16 ZDOTC
      REAL*8 ASYM,F,GETIMAG
      real*8 evl(200),rwork(200)
      complex*16 cwork(200)
C
C---------------Executable Code --------------------------------------
C
c     CALL PST('Generating complex Q_kl+')
c     WRITE(IW,*) ' **** Symmetry',KREP
C
C  LENA is length of affinity blocks (U or V)
C  LENI is length of ionization blocks (U or V)
C
      LENA = NOVVT(KREP)
      LENI = NVOOT(KREP)
C
C calculate Q+_oo + Q-_oo (zero for Q-) for a given symmetry KREP such
C that Q_kl contains only totally symmetric contributions.
C ATT: we need zdotc because we always multiply a^* and b !
C
      OFF1 = JJOO(KREP,KREP) + 1
      DO JX=1,NO(KREP)
        OFFJJ = (JX-1)*LENA + 1
        DO IX=1,NO(KREP)
          OFFII = (IX-1)*LENA + 1
          E=ZDOTC(LENA,VAH(OFFII),1,VAH(OFFJJ),1)
          QOO(OFF1) = -E
          OFF1 = OFF1 + 1
        ENDDO
      ENDDO
c  debug
c     nx=no(krep)*no(krep)
c     off1=jjoo(krep,krep) + 1
c     WRITE(IW,*) 'QOO chksum:',zdotc(nx,qoo(off1),1,qoo(off1),1)
c  debug end
C
C calculate Q+_vo + Q-_vo for a given symmetry KREP such that
C Q_kl contains only totally symmetric contributions.
C

      OFF1 = JJVO(KREP,KREP) + 1
      DO JX=1,NO(KREP)
        JXE = IO(KREP) + JX
        OFFVAFF_J=(JX-1)*LENA + 1
        OFFUION_J=(JX-1)*LENI + 1
        DO AX=1,NV(KREP)
          AXE = IO(NREP+1) + IV(KREP) + AX

          OFFUAFF_A=(AX-1)*LENA + 1
          OFFVION_A=(AX-1)*LENI + 1
 
          E1=ZDOTC(LENA,UAP(OFFUAFF_A),1,VAH(OFFVAFF_J),1)
          E2=ZDOTC(LENI,VIP(OFFVION_A),1,UIH(OFFUION_J),1)
          F=EPS(AXE) - EPS(JXE)
          F=-1.0D0/F
          QVOA(OFF1) = DCMPLX(F,0.0D0)*E1
          QVOI(OFF1) = DCMPLX(F,0.0D0)*E2
          OFF1 = OFF1 + 1
        ENDDO
      ENDDO
c  debug
c     nx=no(krep)*nv(krep)
c     off1=jjvo(krep,krep) + 1
c     WRITE(IW,*) 'QVOA chksum:',zdotc(nx,qvoa(off1),1,qvoa(off1),1)
c     WRITE(IW,*) 'QVOI chksum:',zdotc(nx,qvoi(off1),1,qvoi(off1),1)
c  debug end
C
C  add +/- contributions of QVO in corresponding symmetry
C
      OFF1 = JJVO(KREP,KREP) + 1
      DO JX=1,NO(KREP)
        DO AX=1,NV(KREP)
          QVO(OFF1) = QVOI(OFF1) + QVOA(OFF1)
          OFF1 = OFF1 + 1
        ENDDO
      ENDDO
c  debug
c     nx=no(krep)*nv(krep)
c     off1=jjvo(krep,krep) + 1
c     WRITE(IW,*) 'QVO chksum:',zdotc(nx,qvo(off1),1,qvo(off1),1)
c  debug end
C
C calculate Q+_vv + Q-_vv (zero for Q+) for a given symmetry KREP such
C that Q_kl contains only totally symmetric contributions.
C
      OFF1 = JJVV(KREP,KREP) + 1
      DO BX=1,NV(KREP)
        OFFBB = (BX-1)*LENI + 1
        DO AX=1,NV(KREP)
          OFFAA = (AX-1)*LENI + 1
          E=ZDOTC(LENI,VIP(OFFAA),1,VIP(OFFBB),1)
          QVV(OFF1) = E
          OFF1 = OFF1 + 1
        ENDDO
      ENDDO
c  debug
c     nx=nv(krep)*nv(krep)
c     off1=jjvv(krep,krep) + 1
c     WRITE(IW,*) 'QVV chksum:',zdotc(nx,qvv(off1),1,qvv(off1),1)
c  debug end

C
C  write out Q_kl matrices dependent on PRINT level
C
c     IF(ADCPRINT(5).NE.0) THEN
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ Q_kl (complex) in Symmetry:',KREP
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ h/h part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL CMATOUT(QOO(JJOO(KREP,KREP) +1),NO(KREP),NO(KREP))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ p/h part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL CMATOUT(QVO(JJVO(KREP,KREP) +1),NV(KREP),NO(KREP))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ p/p part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL CMATOUT(QVV(JJVV(KREP,KREP) +1),NV(KREP),NV(KREP))
c     ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QKLDUP(QOO,QVO,QVV,KREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Duplicates the missing Q_kl entries for symmetry KREP+1
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 QOO(*),QVO(*),QVV(*)
      INTEGER KREP
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C
C---------------Executable Code --------------------------------------
C

      YOO=NO(KREP)*NO(KREP)
      YVO=NV(KREP)*NO(KREP)
      YVV=NV(KREP)*NV(KREP)

      OFFOOS = JJOO(KREP,KREP)+1
      OFFVOS = JJVO(KREP,KREP)+1
      OFFVVS = JJVV(KREP,KREP)+1

      OFFOOD = JJOO(KREP+1,KREP+1)+1
      OFFVOD = JJVO(KREP+1,KREP+1)+1
      OFFVVD = JJVV(KREP+1,KREP+1)+1

      CALL XCOPY(YOO,QOO(OFFOOS),1,QOO(OFFOOD),1)
      CALL XCOPY(YVO,QVO(OFFVOS),1,QVO(OFFVOD),1)
      CALL XCOPY(YVV,QVV(OFFVVS),1,QVV(OFFVVD),1)

c     WRITE(IW,*) ' **** QKL check in QKLDUP symmetry ',KREP+1
c     IF(.TRUE.) THEN
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ Q_kl in Symmetry:',KREP
c       WRITE(IW,*) '     ~~~      (in eV)            ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ h/h part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL MATOUTEV(QOO(JJOO(KREP,KREP)+1),NO(KREP),NO(KREP))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ p/h part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL MATOUTEV(QVO(JJVO(KREP,KREP)+1),NV(KREP),NO(KREP))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ p/h part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL MATOUTEV(QVV(JJVV(KREP,KREP)+1),NV(KREP),NV(KREP))
c     ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CQKLDUP(QOO,QVO,QVV,KREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Duplicates the missing COMPLEX Q_kl entries for symmetry KREP+1
C     ATT !!
C     For the complex case where Kramers partners exist (in this program
C     suite only for C_s symmetry) the KREP+1 Qkl is the complex
C     conjugate of the KREP Q_kl !
C     This was found by actually calculating each symmetry separately.
C     Still, Q_kl duplication can be utilized.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 QOO(*),QVO(*),QVV(*)
      INTEGER KREP
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C
C---------------Executable Code --------------------------------------
C
      IF(.NOT.CARITH) CALL QUIT('Complex Q_kl in real mode?')

      YOO=NO(KREP)*NO(KREP)
      YVO=NV(KREP)*NO(KREP)
      YVV=NV(KREP)*NV(KREP)

      OFFOOS = JJOO(KREP,KREP)+1
      OFFVOS = JJVO(KREP,KREP)+1
      OFFVVS = JJVV(KREP,KREP)+1

      OFFOOD = JJOO(KREP+1,KREP+1)+1
      OFFVOD = JJVO(KREP+1,KREP+1)+1
      OFFVVD = JJVV(KREP+1,KREP+1)+1

      CALL ZCOPY(YOO,QOO(OFFOOS),1,QOO(OFFOOD),1)
      CALL ZCOPY(YVO,QVO(OFFVOS),1,QVO(OFFVOD),1)
      CALL ZCOPY(YVV,QVV(OFFVVS),1,QVV(OFFVVD),1)
c
c  for the related kramers partner which is the CC of the original.
c
      CALL CONJUGA(YOO,QOO(OFFOOD),1)
      CALL CONJUGA(YVO,QVO(OFFVOD),1)
      CALL CONJUGA(YVV,QVV(OFFVVD),1)

c     IF(.FALSE.) THEN
c       WRITE(IW,*) ' **** CQKLDUP check for symmetry ',KREP,'/',KREP+1
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ h/h part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL CMATOUT(QOO(JJOO(KREP,KREP)+1),NO(KREP),NO(KREP))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL CMATOUT(QOO(JJOO(KREP+1,KREP+1)+1),NO(KREP+1),NO(KREP+1))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ p/h part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL CMATOUT(QVO(JJVO(KREP,KREP)+1),NV(KREP),NO(KREP))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL CMATOUT(QVO(JJVO(KREP+1,KREP+1)+1),NV(KREP+1),NO(KREP+1))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       WRITE(IW,*) '     ~~~ p/p part                ~~~'
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL CMATOUT(QVV(JJVV(KREP,KREP)+1),NV(KREP),NV(KREP))
c       WRITE(IW,*) '     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c       CALL CMATOUT(QVV(JJVV(KREP+1,KREP+1)+1),NV(KREP+1),NV(KREP+1))
c     ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKE_BPQ1(QOO,QVO,QVV,
     &                    BOO,BVO,
     &                    OOOO1,OOOO2,VOOO1,VOOO2,
     &                    VVOO,VOVO1,VOVO2,
     &                    VOVV1,VOVV2,ADCPRINT)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates the B_pq values by contraction over Q_kl with the 
C     corresponding integrals
C  ATTENTION ****
C     According to the BKC the Q_kl have to be used as Q_lk. Otherwise
C     a correct contraction over ket/bra pairs is not possible.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 QOO(*),QVO(*),QVV(*)
      REAL*8 BOO(*),BVO(*)
      REAL*8 OOOO1(*),OOOO2(*),VOOO1(*),VOOO2(*)
      REAL*8 VVOO(*),VOVO1(*),VOVO2(*)
      REAL*8 VOVV1(*),VOVV2(*)
      INTEGER ADCPRINT(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
!     COMPLEX*16 zdotc
      REAL*8 ddot,getimag
      INTEGER JVVVO_I48(MXREP+1)   !MXREP is defined in symm.inc !
      INTEGER IVOVVTAUX
C---------------Executable Code --------------------------------------
C
C  preparation of aux array by implicit type conversion
C
      DO I=1,NREP+1
        JVVVO_I48(I) = JVVVO(I)
      ENDDO
C****************
C***** b_oo
C****************
C
C  We do the B_oo/VO first since missing contributions are added by
C  by adding the transpose (hermitian conjugate) to itself.
C  Therefore no other contributions may be there at the beginning !
C  read the VOOO integrals V_ai,lj and sort to V_al,ij. contract with
C  Q_al to b_ij and store in b_oo !  This
C  contribution is not symmetric in ij as can be seen from the original
C  formula. But the OV contribution yiels b*_ji ! ==> we add the
C  hermitian conjugate to this contribution
C  *** ATT ***
C  Obviously, in the schirmer/angonoa paper the contraction over Q_kl
C  violates BKC. we therefore CC the Q_kl obtaining Q_lk  :-(((((
C
c     IF(CARITH) THEN
c       CALL PST('Calculating b_pq inhomogeneities (complex)+')
c     ELSE
c       CALL PST('Calculating b_pq inhomogeneities (real)+')
c     ENDIF

      CALL GETVOOO(VOOO1)
      IF(CARITH) WRITE(IW,*) 'VOOO imag. contr.:',
     &           GETIMAG(IVOOOT(NREP+1),VOOO1)
      CALL SRT1ST4 (NREP,MULTB,LTR,LFA,NV,NO,NO,NO,MVO,J2VOOO,JJVO,JJOO,
     &              VOOO1,VOOO2)
      CALL XGEMV ('T',MVO(1),MOO(1),A1,VOOO2,MVO(1),QVO,1,A0,BOO,1)


      DO KREP=1,NREP
        OFF=JJOO(KREP,KREP)*RCW + 1
        IF(CARITH) THEN
          CALL CMATSYM(BOO(OFF),NO(KREP),NO(KREP))
        ELSE
          CALL MATSYM(BOO(OFF),NO(KREP),NO(KREP))
        ENDIF
      ENDDO
        
c     if(carith) then
c       write(*,*) 'boo/vo checksum(c):',zdotc(moo(1),boo,1,boo,1)
c     else
c       write(*,*) 'boo/vo checksum(r):',ddot(moo(1),boo,1,boo,1)
c     endif

C
C  sort the V_ik,jl (i>k,j>l) to V_ij,kl and take only first (totally symmetric)
C  irrep in kl then contract with Q_kl to b_ij. The resulting b_ij
C  contribution is symmetric in ij  and is added to the b_oo/VO+OV part
C  from above.
C
      CALL GETOOOO(OOOO1)
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NO,NO,NO,NO,MOO,J2OOOO,JJOO,JJOO,
     &              OOOO1,OOOO2)
      CALL XGEMV ('N',MOO(1),MOO(1),A1,OOOO2,MOO(1),QOO,1,A1,BOO,1)
c     if(carith) then
c       write(*,*) 'boo/oo checksum(c):',zdotc(moo(1),boo,1,boo,1)
c     else
c       write(*,*) 'boo/oo checksum(r):',ddot(moo(1),boo,1,boo,1)
c     endif
C
C  sort the V_ai,bj integrals to V_ab,ij and contract with Q_ab to b_ij
C  this contribution is symmetric and is added as b_oo/VV to the part
C  from above.
C
      CALL GETVOVO(VOVO1)
      CALL SRT1SS4F (NREP,MULTB,LTR,NV,NO,NV,NO,MVV,J2VVOO,JJVV,JJOO,
     &               VOVO1,VOVO2)
      CALL XGEMV ('T',MVV(1),MOO(1),A1,VOVO2,MVV(1),QVV,1,A1,BOO,1)

      IF(ADCPRINT(6).NE.0) THEN
        WRITE(IW,*) '   BOO vector:'
        DO IREP=1,NREP
          IF(CARITH) THEN
           CALL CMATOUT(BOO(JJOO(IREP,IREP)*RCW +1),NO(IREP),NO(IREP))
          ELSE
            CALL MATOUT(BOO(JJOO(IREP,IREP)+1),NO(IREP),NO(IREP))
          ENDIF
          WRITE(IW,*) '-------------------'
        ENDDO
      ENDIF
C
C****************
C***** b_vo
C****************
C
C
C  VO/OO
C  read the V_ak,jl and sort to V_aj,kl ==> b_aj
C
      CALL GETVOOO(VOOO1)
      CALL SRT1ST4 (NREP,MULTB,LTR,LFA,NV,NO,NO,NO,MVO,J2VOOO,JJVO,JJOO,
     &              VOOO1,VOOO2)
      CALL XGEMV ('N',MVO(1),MOO(1),A1,VOOO2,MVO(1),QOO,1,A0,BVO,1)
c     if(carith) then
c       write(*,*) 'bvo/oo checksum(c):',zdotc(mvo(1),bvo,1,bvo,1)
c     else
c       write(*,*) 'bvo/oo checksum(r):',ddot(mvo(1),bvo,1,bvo,1)
c     endif
C
C  VO/VO
C  read the V_ab,jl and sort to V_aj,bl. contract over bl with Q_bl
C
      CALL GETVVOO(VVOO)
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,
     &              MVO,J2VOVO,JJVO,JJVO,VVOO,VOVO2)
      CALL XGEMV ('N',MVO(1),MVO(1),A1,VOVO2,MVO(1),QVO,1,A1,BVO,1)
c     if(carith) then
c       write(*,*) 'bvo/vo checksum(c):',zdotc(mvo(1),bvo,1,bvo,1)
c     else
c       write(*,*) 'bvo/vo checksum(r):',ddot(mvo(1),bvo,1,bvo,1)
c     endif
C
C  VO/OV
C  read the V_al,bj and negate them (V_al,jb in the formula)
C  sort to V_aj,bl:blrep=1 and contract with Q_bl to b_aj
C  In the complex case Q_bl has to be conjugated before contraction
C  and is restored after array has been used.
C
      CALL GETVOVO(VOVO1)
      CALL XSCAL(NV4,-A1,VOVO1,1)
      CALL SRT16 (NREP,MULTB,LTR,LFA,NV,NO,NV,NO,MVO,JVOVO,
     &            JJVO,JJVO,VOVO1,VOVO2)
      IF(CARITH) CALL CONJUGA(MVO(1),QVO,1)
      CALL XGEMV ('N',MVO(1),MVO(1),A1,VOVO2,MVO(1),QVO,1,A1,BVO,1)
      IF(CARITH) CALL CONJUGA(MVO(1),QVO,1)
c     if(carith) then
c       write(*,*) 'bvo/ov checksum(c):',zdotc(mvo(1),bvo,1,bvo,1)
c     else
c       write(*,*) 'bvo/ov checksum(r):',ddot(mvo(1),bvo,1,bvo,1)
c     endif
C
C  VO/VV
C  read the V_cj,ba (VOVV1) , conjugate and resort to V_cb,aj (VOVV2)
C  this is done for the VOVV integrals in SRT20M where IMPLICITLY only
C  the first irrep is taken ! 
C  afterwards contract over Q_bc = Q^*_cb. restored after contraction
C
      CALL RDVOVV (VOVV1)
      IVOVVTAUX = IVOVVT(NREP+1)
      IF(CARITH) CALL CONJUGA(IVOVVTAUX,VOVV1,1)
      CALL SRT20M(NREP,MULTB,NVO,NV,NO,NV,NV,MVV,
     &            JVVVO_I48,JJVV,JJVO,VOVV1,VOVV2)
      IF(CARITH) CALL CONJUGA(MVV(1),QVV,1)
      CALL XGEMV ('T',MVV(1),MVO(1),A1,VOVV2,MVV(1),QVV,1,A1,BVO,1)
      IF(CARITH) CALL CONJUGA(MVV(1),QVV,1)
c     if(carith) then
c       write(*,*) 'bvo final checksum(c):',zdotc(mvo(1),bvo,1,bvo,1)
c     else
c       write(*,*) 'bvo final checksum(r):',ddot(mvo(1),bvo,1,bvo,1)
c     endif

c     write(iw,*) '   BVO finished'

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SIGSOLV(BOO,BVO,
     &                   VVOO,VOOO1,VOOO2,VOOO3,
     &                   VOVO1,VOVO2,VOVO3,
     &                   EPS,IPIV,SIGMA,ADCPRINT)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Solves the linear system of equations for the Sigma_pq
C     Since the system is coupling different symmetries of Sigma_pq 
C     the equations yield all symmetries of Sigma_pq
C
C   !! ATT !! all necessary entities for the full Sigma_pq are provided
C   and calculated already in the previous routines ! here only the
C   Sigma_hh part for the non-Dyson ionization part is calculated !
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BOO(*),BVO(*)
      REAL*8 VVOO(*),VOOO1(*),VOOO2(*),VOOO3(*)
      REAL*8 VOVO1(*),VOVO2(*),VOVO3(*)
      REAL*8 EPS(*),SIGMA(*)
      INTEGER IPIV(*),ADCPRINT(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      real*8 one,zero,ddot
C
C---------------Executable Code --------------------------------------
C
c     CALL PST('Solving real equation system for Sigma(inf)+')
C
C  get the V_ab,ij integrals and sort to V_ai,bj
C
      CALL GETVVOO(VVOO)
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,
     &                 MVO,J2VOVO,JJVO,JJVO,VVOO,VOVO1)
C
C  then get the V_aj,bi integrals and negate them. sort to V_ai,bj
C
      CALL GETVOVO(VOVO2)
      CALL XSCAL (NV4,-A1,VOVO2,1)
      CALL SRT16 (NREP,MULTB,LTR,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            VOVO2,VOVO3)
C
C  add the buffer VOVO1 to VOVO3 in symmetry 1 : 
C  B_ai,bj = V_ab,ij - V_aj,bi
C
      N1 = MVO(1) * MVO(1)
      CALL XAXPY(N1,A1,VOVO1,1,VOVO3,1)
C
C  then divide B_ai,bj by (ej - eb) yielding
C  A_ai,bj = (V_ab,ij - V_aj,bi)/(ej - eb)
C  multiply result by (-1) yielding -A_ai,bj and add
C  unity matrix. Result is 1-A_ai,bj in VOVO3
C
      CALL SIGDIV (EPS,VOVO3,VOVO3,MVO(1),MVO(1))
      CALL XSCAL  (N1,-A1,VOVO3,1)
      CALL ADDONE (VOVO3,VOVO3,MVO(1),MVO(1))
C
C  now form the inverse of the (1-A_phph) matrix.
C  Done in a two-step process by first LU factorizing and then
C  backsubstituting via LAPACK routines. the backsubstituting routine
C  needs work space: we provide the VOVO2 buffer which should be
C  plenty and enough. length is NV4 therefore.
C
      N = MVO(1)
      CALL DGETRF(N,N,VOVO3,N,IPIV,INFO)
      IF (INFO.NE.0)
     &  CALL QUIT('Error during LU factorization in SIGSOLV')

      CALL DGETRI(N,VOVO3,N,IPIV,VOVO2,NV4,INFO)
      IF (INFO.NE.0)
     &  CALL QUIT('Error during backsubstitution in SIGSOLV')
C
C  from here we have (1-A_ph,ph)^-1 in VOVO3
C

C  next we form A_hh,ph
C  get V_bm,ji (VOOO1) ==> V_bj,mi (VOOO2)
C  get -V*_bi,mj (VOOO1) ==> -V*_bj,mi (VOOO3)
C  result in VOOO3
C
      CALL GETVOOO(VOOO1)
      CALL SRT1ST4(NREP,MULTB,LTR,LFA,NV,NO,NO,NO,MVO,J2VOOO,JJVO,JJOO,
     &              VOOO1,VOOO2)
      CALL XSCAL (NV2,-A1,VOOO1,1)
      IF (CARITH) CALL CONJUGA(NV2,VOOO1,1)
      CALL SRT26 (NREP,MULTB,LTR,LFA,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &            VOOO1,VOOO3)
      CALL XAXPY(MVO(1)*MOO(1),A1,VOOO2,1,VOOO3,1)
C
C  create MOO(1) * MVO(1) matrix by special transposition routine
C  Now the array in VOOO2 has an MOO(1) * MVO(1) organization.
C  divide second index (p/h pair) by the p/h denominators again.
C  Result A_hh,ph in VOOO2
C
      CALL SIGTRA(VOOO3,VOOO3,VOOO2,VOOO2,MVO(1),MOO(1))
      CALL SIGDIV(EPS,VOOO2,VOOO2,MOO(1),MVO(1))
C
C  next we multiply the A_hh,ph matrix in VOOO2 with (1-A_phph)^(-1)
C  in VOVO3 yielding X_mi,bj. Then we contract X_mi,bj * b_bj and add
C  b_mi finally yielding Sigma(infinity).
C
      M=MOO(1)
      N=MVO(1)
      K=MVO(1)
      CALL XGEMM('N','N',M,N,K,A1,VOOO2,M,VOVO3,K,A0,VOOO3,M)
      CALL XGEMV ('N',M,N,A1,VOOO3,M,BVO,1,A0,SIGMA,1)
      CALL XAXPY (M,A1,BOO,1,SIGMA,1)
C
C  print out sigma(infinity) matrix
C
      IF(ADCPRINT(7).NE.0) THEN
        CALL PST('Constant diagrams (r) for all symmetries (eV)+')
        WRITE(IW,*)
        WRITE(IW,*) (NO(IX),IX=1,NREP)
        WRITE(IW,*)
        DO KREP=1,NREP
          OFF=JJOO(KREP,KREP)*RCW + 1
          IF (CARITH) THEN
            CALL CMATOUT(SIGMA(OFF),NO(KREP),NO(KREP))
          ELSE
            CALL MATOUT(SIGMA(OFF),NO(KREP),NO(KREP))
          ENDIF
        ENDDO 
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CSIGSOLV(BOO,BVO,XVO,YVO,BVO2,
     &                    VVOO,VOOO1,VOOO2,VOOO3,
     &                    VOVO1,VOVO2,VOVO3,
     &                    EPS,IPIV,SIGMA,MC,ADCPRINT)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Solves the COMPLEX linear system of equations for the Sigma_pq
C     Since the system is coupling different symmetries of Sigma_pq
C     the equations yield all symmetries of Sigma_pq
C
C   !! ATT !! all necessary entities for the full Sigma_pq are provided
c   that means p and q can sweep over hole AND particle indices.
C   here only the
C   Sigma_hh part for the non-Dyson ionization part is calculated !
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BOO(*),BVO(*),XVO(*),YVO(*),BVO2(*)
      REAL*8 VVOO(*),VOOO1(*),VOOO2(*),VOOO3(*)
      REAL*8 VOVO1(*),VOVO2(*),VOVO3(*)
      REAL*8 EPS(*),SIGMA(*)
      COMPLEX*16 MC(*)
      INTEGER IPIV(*),ADCPRINT(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 ASYM,GETIMAG
C
C---------------Executable Code --------------------------------------
C
c     CALL PST('Solving complex equation system for Sigma(inf)+')
c     write(iw,*) 'BOO imag contr:',getimag(moo(1),boo)
c     write(iw,*) 'BVO imag contr:',getimag(mvo(1),bvo)
C
C  get the VOVO integrals V_aj,bi and sort to V_ai,bj
C     VOVO2 --> VOVO3: ATT ! only the totally symmetric
C  irrep in bj is selected !!!
C  then divide the V_ai,bj by (e_j - e_b).
C  conjugate for the corresponding lower matrix entries.
C
C               ============> UPPER LEFT/LOWER RIGHT <===============
C
      CALL GETVOVO(VOVO2)
      write(iw,*) 'VOVO imag contr:',getimag(ivovo(nrep+1),vovo2)
      CALL SRT16 (NREP,MULTB,LTR,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            VOVO2,VOVO3)
      CALL SIGDIV (EPS,VOVO3,VOVO3,MVO(1),MVO(1))
      CALL CSIGBLD (MC,VOVO3,2*MVO(1),MVO(1),1)
      CALL CONJUGA (MVO(1)*MVO(1),VOVO3,1)
      CALL CSIGBLD (MC,VOVO3,2*MVO(1),MVO(1),4)
C
C               ============> END UPPER LEFT/LOWER RIGHT <===============
C
C
C               ============> UPPER RIGHT/LOWER LEFT <===============
C
C  get the V_ab,ij integrals, negate and sort to V_ai,bj
C  VVOO --> VOVO1  and divide by (e_j - e_b) ==> B_ai,bj
C  only the totally symmetric irrep is selected !!!
C  conjugate for the corresponding lower matrix entries.
C
      CALL GETVVOO(VVOO)
      write(iw,*) 'VVOO imag contr:',getimag(ivvoott(nrep+1),vvoo)
      CALL XSCAL (NV3,-A1,VVOO,1)
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,
     &                 MVO,J2VOVO,JJVO,JJVO,VVOO,VOVO3)
      CALL SIGDIV (EPS,VOVO3,VOVO3,MVO(1),MVO(1))
      CALL CSIGBLD (MC,VOVO3,2*MVO(1),MVO(1),2)
      CALL CONJUGA (MVO(1)*MVO(1),VOVO3,1)
      CALL CSIGBLD (MC,VOVO3,2*MVO(1),MVO(1),3)
C
C               ============> END UPPER RIGHT/LOWER LEFT <===============
C
C  now complex matrix is ready. We add 1 to the diagonal , check for
C  coefficient relations and call the
C  complex LU-inversion routines.
C
      N=2*MVO(1)
      CALL ADDONE(MC,MC,N,N)
      WRITE(IW,*) 'Coefficient matrix:'
      CALL MCCHK(MC,2*MVO(1),MVO(1),ASYM)

      CALL ZGETRF(N,N,MC,N,IPIV,INFO)
      IF (INFO.NE.0) THEN
        CALL QUIT('Error during LU factorization in CSIGSOLV')
      ENDIF
      CALL ZGETRI(N,MC,N,IPIV,VOVO2,NV4,INFO)
      IF (INFO.NE.0) THEN
        CALL QUIT('Error during backsubstitution in CSIGSOLV')
      ENDIF
      WRITE(IW,*) 'Inverse matrix:'
      CALL MCCHK(MC,2*MVO(1),MVO(1),ASYM)
      WRITE(IW,*)
      WRITE(IW,*)
C
C  next construct the vector B_vo/B*_vo and multiply inverse matrix on
C  it for Sigma_ai. The inverse matrix is in MC and has dimensions
C  MVO(1) x MVO(1). The final Sigma_ai/Sigma*_ai solution vector
C  is temporarily stored in VOVO2 and the unstarred/starred
C  parts then transferred to XVO/YVO
C
C  important: from the system of equations we see that first b_aj is
C  used and then b_ja = b*_aj. We therefore need to resort the b-array
C  via the SRT1C routine. iThis routine is necessary because first come
C  the b_aj of irrep 1x1 then of irrep 2x2 ... nrep x nrep.
C  the auxiliary buffer for that is XVO which has
C  the same dimensions as the BVO buffer. the complete vector then sits
C  in BVO2(*).
C
      N=MVO(1)
      CALL XCOPY (N,BVO(1),1,BVO2(1),1)
      CALL XCOPY (N,BVO(1),1,BVO2(1+N*RCW),1)
      CALL CONJUGA (N,BVO2(1+N*RCW),1)
      CALL CHKVECID(N,BVO2(1),BVO2(1+N*RCW),ASYM)
      WRITE(IW,*) 'B-Vector deviation:',ASYM
      
      M=2*MVO(1)
      N=1
      K=2*MVO(1)
      CALL XGEMM('N','N',M,N,K,A1,MC,M,BVO2,K,A0,VOVO2,M)
      N=MVO(1)
      CALL XCOPY (N,VOVO2(1),1,XVO(1),1)
      CALL XCOPY (N,VOVO2(1+N*RCW),1,YVO(1),1)
      CALL CHKVECID(N,XVO,YVO,ASYM)
      WRITE(*,*) 'SIGMA_ph/SIGMA_hp  Asymmetry:',ASYM
c
c  from here SIGMA_aj and SIGMA_ja are available
c  in the vectors XVO and YVO
c
c   ----- from here form Sigma_mi, the h/h part with Sigma_aj
c         which is the solution vector in XVO
c
C  get VOOO integrals as V_bi,mj, negate and conjugate them
C  and sort to -V*_bj,mi. reorganize and divide by (e_j - e_b)
C  Afterwards its multiplied with Sigma_bj. Initial result is
C  stored in Sigma_hh.
C
      CALL GETVOOO(VOOO1)
      CALL XSCAL (NV2,-A1,VOOO1,1)
      CALL CONJUGA(NV2,VOOO1,1)
      CALL SRT26 (NREP,MULTB,LTR,LFA,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &            VOOO1,VOOO3)
      CALL SIGTRA(VOOO3,VOOO3,VOOO2,VOOO2,MVO(1),MOO(1))
      CALL SIGDIV(EPS,VOOO2,VOOO2,MOO(1),MVO(1))

      M = MOO(1)
      N = MVO(1)
      CALL XGEMV ( 'N', M, N, A1, VOOO2, M, XVO, 1,
     &                   A0, SIGMA, 1)
C
C  next get the V_bm,ji integrals and sort to V_bj,mi
C  transpose matrix and divide by energies.
c
c   ----- from here form Sigma_mi, the h/h part with Sigma_ja
c         which is the solution vector in YVO is in the aj order
C
      CALL GETVOOO(VOOO1)
      CALL SRT1ST4(NREP,MULTB,LTR,LFA,NV,NO,NO,NO,MVO,J2VOOO,JJVO,JJOO,
     &              VOOO1,VOOO3)
      CALL SIGTRA(VOOO3,VOOO3,VOOO2,VOOO2,MVO(1),MOO(1))
      CALL SIGDIV(EPS,VOOO2,VOOO2,MOO(1),MVO(1))
c
c  no resorting because second solution vector is in the right order
c
c     CALL SRT1C1N(NREP,NO,NV,YVO,XVO)

      M = MOO(1)
      N = MVO(1)
      CALL XGEMV ( 'N', M, N, A1, VOOO2, M, YVO, 1,
     &                   A1, SIGMA, 1)
C
C  finally add the b_hh contribution to Sigma.
C  and conjugate them because we need it in |><| convention
C
      CALL XAXPY  (MOO(1),A1,BOO,1,SIGMA,1)
      CALL CONJUGA(MOO(1),SIGMA,1)

      IF(ADCPRINT(7).NE.0) THEN
        WRITE(IW,*) '***'
        WRITE(IW,*) '*** Constant diagrams (c) for all symmetries (eV):'
        WRITE(IW,*) '***'
        WRITE(IW,*)
        WRITE(IW,*)
        DO KREP=1,NREP
          WRITE(IW,*) 'Symmetry ',KREP
          WRITE(IW,*)
          OFF=JJOO(KREP,KREP)*RCW + 1
          CALL CMATOUT(SIGMA(OFF),NO(KREP),NO(KREP))
        ENDDO
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MCCHK(V,LDA,NN,ASYM)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Verifies relations of coefficient matrix. The coefficients
C     have to be complex conjugate to each other.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 V(LDA,LDA),H1,H2,H3
      INTEGER LDA,NN
      REAL*8 ASYM1,ASYM2,ASYM
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/files.inc"
#include  "../relccsd/param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable Code --------------------------------------
C
      IF(LDA.NE.2*NN) THEN
        WRITE(IW,*) 'Illegal dimensions of coeff. matrix !'
        CALL QUIT('MCCHK')
      ENDIF

c
c  compare upper left/lower right and upper right/lower left
C
      H1=A0
      H2=A0
      H3=A0

      DO J=1,NN
        DO I=1,NN
          H1 = H1 + V(I,J) - DCONJG(V(I+NN,J+NN))
          H2 = H2 + V(I,J+NN) - DCONJG(V(I+NN,J))
          H3 = H3 + V(I,J) + DCONJG(V(I,J))
        ENDDO
      ENDDO
      ASYM1 = DSQRT(DBLE(H1*DCONJG(H1)))
      ASYM2 = DSQRT(DBLE(H2*DCONJG(H2)))
      ASYM = ASYM1 + ASYM2

      WRITE(IW,*) 'Asymmetry 1:',ASYM1
      WRITE(IW,*) 'Asymmetry 2:',ASYM2
      WRITE(IW,*) 'Total asymmetry:',ASYM
      WRITE(IW,*) 'Matrix norm:',H3

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION GETIMAG(M,V)
c
c  this subroutine gets the quadratic sum of the imaginary part
c  of a vector
c
      IMPLICIT INTEGER (A-Z)
      INTEGER M
      REAL*8 GETIMAG
      COMPLEX*16 V(*)
c
c   local variables
c
      REAL*8 SUM,H
c
c   start execution
c
      sum=0.0d0
      do i=1,m
        h=dimag(v(i))
        sum = sum + h*h
      enddo

      getimag = sum

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CHKVECID(N,V1,V2,ASYM)
      IMPLICIT INTEGER (A-Z)
      INTEGER N
      COMPLEX*16 V1(N),V2(N),H
      REAL*8 ASYM

      H=(0.0D0,0.0D0)
      DO I=1,N
        H=H + (V1(I) - DCONJG(V2(I)))
      ENDDO
      ASYM = DSQRT(DBLE(DCONJG(H)*H))

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CSIGBLD(DARRAY,SARRAY,LDAD,LDAS,IPLACE)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Constructs the enlarged coefficient matrix for inversion in
C     the complex ADC-matrix case. 
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 DARRAY(LDAD,LDAD),SARRAY(LDAS,LDAS)
      INTEGER LDAD,LDAS,IPLACE
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable Code --------------------------------------
C
      IF(LDAS*2.NE.LDAD) THEN
        WRITE(IW,*) 'Source matrix will not fit in dest. array'
        CALL QUIT('CSIGBLD')
      ENDIF

      IF(IPLACE.EQ.1) THEN
        NROFF=0
        NCOFF=0
      ELSE IF(IPLACE.EQ.2) THEN
        NROFF=0
        NCOFF=LDAS
      ELSE IF(IPLACE.EQ.3) THEN
        NROFF=LDAS
        NCOFF=0
      ELSE IF(IPLACE.EQ.4) THEN
        NROFF=LDAS
        NCOFF=LDAS
      ELSE
        CALL QUIT('Illegal IPLACE value in CSIGBLD')
      ENDIF

      DO J=1,LDAS
        DO I=1,LDAS
          DARRAY(NROFF+I,NCOFF+J)=SARRAY(I,J)
        ENDDO
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SIGDIV(EPS,A,AC,NROW,NCOL)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C   Divides V_ai,bj matrix by corresponding energies in the second
C   index pair p/h. Number of p/h entries: MVO(1) necessarily !
C   Also suitable for rectangular matrices, but columns have to be of
C   ph-type.
C   by the double allocation complex matrices can be handled as well.
C   in either case the other matrix is not referenced
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER NROW,NCOL
      REAL*8 EPS(*),A(NROW,NCOL)
      COMPLEX*16 AC(NROW,NCOL)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 EJ,EB,FAC
C
C---------------Executable code--------------------------------------
C
      JC = 0
      DO I=1,NROW
        DO JREP=1,NREP
          EJOFF = IO(JREP)
          DO 20 BREP=1,NREP
            IF(MULTB(BREP,JREP,2).NE.1) GOTO 20
c      write(iw,*) 'SIGDIV jrep/brep:',jrep,brep
            EBOFF = IO(NREP+1) + IV(BREP)
            DO JFIE=1,NO(JREP)
              DO BFIE=1,NV(BREP)
                EJ=EPS(EJOFF + JFIE)
                EB=EPS(EBOFF + BFIE)
                FAC=1.0D0/(EJ - EB)
                JC = JC+1
                IF(CARITH) THEN
                  AC(I,JC)=AC(I,JC) * DCMPLX(FAC,0.0D0)
                ELSE
                  A(I,JC)=A(I,JC) * FAC
                        ENDIF
              ENDDO
            ENDDO
 20       CONTINUE
        ENDDO
        IF(JC.NE.NCOL) CALL
     &    QUIT('Inconsistency of columns in SIGDIV')
        JC = 0
      ENDDO
        
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SIGTRA(BUF1,CBUF1,BUF2,CBUF2,NR1,NC1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Transposes rectangular matrix BUF1(NR1,NC1) to 
C     BUF2(NC1,NR1). This is preferable to using a sorting
C     algorithm when a row/column organization is evident.
C     ATT: in the calling program the arrays are purely linear !
C     the interpretation of the logical dimension the is important !
C     it is this reinterpretation which allows for easy coding.
C     the reinterpretation of the leading dimensions leads to contiguous
C     storage in memory !
C
C     To common understanding of the purpose of this routine a complex
C     conjugation does not have to occur !
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER NR1,NC1
      REAL*8 BUF1(NR1,NC1),BUF2(NC1,NR1)
      COMPLEX*16 CBUF1(NR1,NC1),CBUF2(NC1,NR1)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      DO I=1,NR1
      DO J=1,NC1
        IF(CARITH) THEN
          CBUF2(J,I)=CBUF1(I,J)
        ELSE
          BUF2(J,I)=BUF1(I,J)
        ENDIF
      ENDDO
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKE_HH2(CKKS,VVOO1,VVOO2,EPS,EABLKKS,KREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes ADC(2) h/h block.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 CKKS(*),VVOO1(*),VVOO2(*),EPS(*),EABLKKS(*)
      INTEGER KREP
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
C  the v_abij integrals (i.e. the VVOO integrals divided
C  by e_a + e_b - e_i - e_j). v_abij and v*_abij are afterwards
C  available as IJK,L:LREP in VVOO1 and VVOO2
C

      CALL GETVVOO(VVOO1)
      CALL DENOMVVOO(EPS,VVOO1,VVOO1)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             VVOO1,VVOO2)
      NX=KVVOOT(NREP+1)
      CALL XCOPY(NX,VVOO2,1,VVOO1,1)
      IF(CARITH) THEN
        CALL CONJUGA(NX,VVOO2,1)
      ENDIF
C
C  Generate the epsilon array of VVOT,O,O structure in the symmetry KREP
C
      CALL EPSARR(EPS,EABLKKS,KREP,ISIZE)
C
C  contract the arrays v_abl,k and v*_ablk' and E_ablkk' in symmetry KREP
C  over the common index range abl and form the C(2)_kk' block
C  Comply to BKC by complex conjugating the VVOO2.
C  The resulting C_kk' block has the BKC row: |>  col: <| =   |><|
C
      CALL CNTRCT3(VVOO1,VVOO2,EABLKKS,CKKS,
     &             NVVOT(KREP),NO(KREP),KVVOOT(KREP))
C  
C  add the occupied spinor energies in symmetry KREP
C  to the diagonal of the h-h block.
C
      CALL ADDEKK(EPS,CKKS,CKKS,KREP)
C
C  write out C_kk matrix
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKE_HH3(KREP,CKKS,CKKSAUX,EPS,VVVV,
     &                   VVOO1,VVOO2,VVOO3,VVOO4,
     &                   VOVO1,VOVO2,VOVO3,OOOO1,
     &                   FKKS,S4P)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates contributions to the h/h block for ADC(3)
C     The resulting BKC has to be |><| for all C(3) contributions.
C
C     Additionally the routine can be used for the calculation of
C     some SIGMA(4+) contributions. If the corresponding flag is set
C     the FKKS array is filled with the third order f_kk' contributions.
C
C     The resulting BKC has to be |><| for all C(3) contributions.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER KREP
      REAL*8 CKKS(*),CKKSAUX(*),EPS(*),VVVV(*)
      REAL*8 VVOO1(*),VVOO2(*),VVOO3(*),VVOO4(*)
      REAL*8 VOVO1(*),VOVO2(*),VOVO3(*),OOOO1(*)
      REAL*8 FKKS(*)
      LOGICAL S4P
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/complex.inc"
#include  "../relccsd/symm.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
      REAL*8 ASYM
C
C---------------Executable Code --------------------------------------
c
      OFFS = JJOO(KREP,KREP)*RCW + 1
      LENC = NO(KREP) * NO(KREP)

C if sig4+ then clear the corresponding position in the f_kk' array

      IF(S4P) CALL XCOPY (LENC,A0,0,FKKS(OFFS),1)

C     get <VV||OO> integrals and form
C          the v_ab,kl   (VVOO1) and
C          the v^*_ab,kl (VVOO2) and
C     these integrals remain there as long as possible since they 
C     are needed multiple times !
C
      CALL GETVVOO(VVOO1)
      CALL DENOMVVOO(EPS,VVOO1,VVOO1)
      CALL XCOPY(NV3,VVOO1,1,VVOO2,1)
      IF(CARITH) THEN
        CALL CONJUGA(NV3,VVOO2,1)
      ENDIF
C
C  .................................................................calculation of C(A)_kk'
C
C
C     contract V_cd,ab with v_ab,kl to A_cd,kl (to buffer VVOO3)
C     in all symmetries ! the prefactor of 1/4 is cancelled due to
C     the double index restriction a>b and c>d
C     In case of S4P the A_cd,kl is divided by the epsilon terms e_cdkl
C    
C     If these routines are used for the generation of SIGMA(4+)
C     some denomiator divisions have to be performed and some hermitian
C     conjugations have to be omitted.
C
      CALL RDVVVV(VVVV)
      CALL CNTRCT('N','N',NVVT,NOOT,NVVT,A1,VVVV,VVOO1,A0,VVOO3,NREP)

      IF(S4P) CALL DENOMVVOO(EPS,VVOO3,VVOO3)
C
C   Sort v^*_cd,k'l to v^*_cdl,k':KSREP   (VVOO2)
C   Sort A_cd,kl to A_cdl,k               (VVOO4)
C   Since the triangular storage mode is lifted all necessary integrals
C   are generated via SRT1T3
C
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             VVOO2,VVOO4)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             VVOO3,VVOO2)
C
C   contract inegrals now over cdl range (twice the index k and l were
C   permuted) and add the result to previous C_kk' matrix
C   the CNTRCT routine must be called with 'T' in the complex
C   case because no ADDITIONAL complex conjugation is allowed !
C  
      CALL CNTRCT('T','N',NO,NO,NVVOT,A1,VVOO2,VVOO4,A0,VVOO3,NREP)
      CALL XAXPY (LENC,A1,VVOO3(OFFS),1,CKKS,1)
      IF (S4P) CALL XAXPY (LENC,A1,VVOO3(OFFS),1,FKKS,1)
C
C  ................................................................calculation of C(B)_kk'
C

C   get VOVO integrals in buffer VOVO1 and negate them
C   (V_lc,bm is written in the mathematical formula !)
C   and the v*_ab,kl integrals to buffer VVOO2
C
      CALL GETVOVO(VOVO1)
      CALL XSCAL (NV4,-A1,VOVO1,1)
      CALL XCOPY(NV3,VVOO1,1,VVOO2,1)
      if(CARITH) THEN
        CALL CONJUGA(NV3,VVOO2,1)
      ENDIF
C
C   v^*_ac,k'm (VVOO2) ==> v^*_ak',cm  (VOVO3)
C   V_cl,bm    (VOVO1) ==> V_cm,bl     (VOVO2)
C   v^*_ak',cm x V_cm,bl  ==> A_ak',bl (VOVO1)
C   A_ak',bl   (VOVO1) ==> A_a>b,k'>l  (VVOO3)
C
      CALL SRT1TT4(NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &             MVO,J2VOVO,JJVO,JJVO,VVOO2,VOVO3)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,
     &            JVOVO,JJVO,JJVO,VOVO1,VOVO2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(JREP)
         K = MVO(JREP)
         OFF2 = JVOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,VOVO3(OFF1),M,VOVO2(OFF2),K,
     &               A0,VOVO1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
C  now sort A_ak',bl (VOVO1) to A_a>b,k'l (VOVO2) only
C  antisymmetrization in the a,b index pair !
C  Then to A_a>bl,k' (VOVO3)
C  Sort the v_a>b,l>k = -v_a>b,k>l to v_a>bl,k (remember minus sign !)
C
      CALL SRT1TS4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,J2VOVO,
     &              JJVO,JJVO,VOVO1,VOVO2)
      CALL SRT19 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &            VOVO2,VOVO3)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             VVOO1,VVOO3)
C
C   now contract v_a>bl,k (VVOO3) and A_a>bl,k' (VOVO3) over abl.
C   Important !! Due to the antisymmetrization in ab for the A array 
C   no additional factor 2 is necessary. Account for the minus sign
C   stemming from the lk permutation in v_a>b,k>l.
C   In the complex case no additional CC must occur ==> 'T'
C   instead of 'C'
C  
      CALL CNTRCT('T','N',NO,NO,NVVOT,-A1,VVOO3,VOVO3,A0,VVOO2,NREP)
      CALL XAXPY (LENC,A1,VVOO2(OFFS),1,CKKS,1)
C
C  ............................................................calculation of C(C)_kk'
C
C     v_ab,lm (VVOO1) x V_lm,jk (OOOO1)  ==> A_ab,jk (VVOO3)
C     A_ab,jk (VVOO3)  ==> A_abj,k (VVOO4)
C     v^*_ab,jk' (VVOO2) ==> v^*_abj,k' (VVOO3)
C     A_abj,k (VVOO4)  x  v^*_abj,k' (VVOO3) ==>  C(C)_kk'
C     create the C_kk' aux array and form the Hermitian conjugate
C     this goes to the final C_kk' array !
C     the factor 1/4 is cancelled by the index restriction
C     a>b and l>m.
C
      CALL GETOOOO(OOOO1)
      CALL CNTRCT('N','N',NVVT,NOOT,NOOT,A1,VVOO1,OOOO1,A0,VVOO3,NREP)
      IF(S4P) CALL  DENOMVVOO(EPS,VVOO3,VVOO3)
      CALL XCOPY(NV3,VVOO1,1,VVOO2,1)
      IF(CARITH) THEN
        CALL CONJUGA(NV3,VVOO2,1)
      ENDIF
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             VVOO3,VVOO4)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             VVOO2,VVOO3)
C
C  again, no additional complex conjugation !
C
      CALL CNTRCT('T','N',NO,NO,NVVOT,A1,VVOO4,VVOO3,A0,VVOO2,NREP)
      CALL XCOPY (LENC,VVOO2(OFFS),1,CKKSAUX,1)
      IF(S4P) CALL XAXPY (LENC,A1,VVOO2(OFFS),1,FKKS,1)
      CALL ADDARR(NO(KREP),CKKSAUX,CKKSAUX,CKKS,CKKS)
      CALL HERMCON(NO(KREP),CKKSAUX,CKKSAUX)
      CALL ADDARR(NO(KREP),CKKSAUX,CKKSAUX,CKKS,CKKS)
C
C  ............................................................calculation of C(D)_kk'
C
C   v_ba,ml (VVOO1) ==> v_bm,al (VOVO1)
C   V^*_ak,cl (VOVO2) ==> V_al,ck (VOVO3)
C   special contraction v_bm,al (VOVO1) * V_al,ck (VOVO3) =>
C   A_bm,ck
C
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,J2VOVO,JJVO,JJVO,VVOO1,VOVO1)
      CALL GETVOVO(VOVO2)
      if(CARITH) THEN
        CALL CONJUGA(NV4,VOVO2,1)
      ENDIF
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,
     &            JVOVO,JJVO,JJVO,VOVO2,VOVO3)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(JREP)
         K = MVO(JREP)
         OFF2 = JVOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,VOVO1(OFF1),M,VOVO3(OFF2),K,
     &               A0,VOVO2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
C  now antisymmetrize A in bc:
C  A_bm,ck (VOVO2) ==> A_b>c,mk (VOVO3) and to A_b>cm,k (VOVO2)
C  Then fetch the v*_b>c,k'>m again to VVOO2
C  and sort v*b>c,m>k' (VVOO2) to v*_b>cm,k' (VVOO3) the minus sign is
C  accounted for at the end of the contraction !
C
      CALL XCOPY(NV3,VVOO1,1,VVOO2,1)
      IF(CARITH) THEN
        CALL CONJUGA(NV3,VVOO2,1)
      ENDIF
      CALL SRT1TS4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,J2VOVO,
     &              JJVO,JJVO,VOVO2,VOVO3)
      CALL SRT1S3  (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &              VOVO3,VOVO2)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             VVOO2,VVOO3)
C
C  again no additional complex conjugation ==> 'T' parameter and not 'C'
C
      CALL CNTRCT('T','N',NO,NO,NVVOT,-A1,VOVO2,VVOO3,A0,VVOO2,NREP)
      CALL XCOPY (LENC,VVOO2(OFFS),1,CKKSAUX,1)
      CALL ADDARR(NO(KREP),CKKSAUX,CKKSAUX,CKKS,CKKS)
      CALL HERMCON(NO(KREP),CKKSAUX,CKKSAUX)
      CALL ADDARR(NO(KREP),CKKSAUX,CKKSAUX,CKKS,CKKS)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKE_HHP2(VOOO1,VOOO2,BUFHP,CREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates 1h/2h-1p contribution for the ADC-2/2x
C     For ADC-3 the UIONHO routine is called that calculates
C     everything at once. This is more organized.
C
C     BKC to be followed: <||>|>,<|
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VOOO1(*),VOOO2(*),BUFHP(*)
      INTEGER CREP
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
#include  "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------

      real*8 ddot,getimag
      complex*16 zdotc
C
C---------------Executable Code --------------------------------------
C
C   storage mode of h/h-2p block: C_ak>l,j corresponds to hermitian
C   conjugate (lower triangular part of the matrix).
C
C  ... read the VOOO integrals and sort according to
C      IJ,K>L:KLREP --> IK>L,J:JREP, Result in VOOO2
C      since we are in lower triangular storage mode we DO NOT
C      additionally cc the VOOO. Only by this we get the correct BKC
C      for the 1h/2h1p block.
C     
      CALL PST('calculating 1h/2h1p coupling block in first order+')

      CALL GETVOOO(VOOO1)

      write(iw,*) 'Imaginary part of the VOOO integrals:',
     &  GETIMAG(ivooot(nrep+1),vooo1)

      CALL SRT6(NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &          VOOO1,VOOO2)
C
C  copy C_aij,k to h/h2p buffer
C

      OFF1 = LOVOOT(CREP)*RCW + 1
      N=NVOOT(CREP)*NO(CREP)
      CALL XCOPY(N,VOOO2(OFF1),1,BUFHP,1)

      if(carith) then
        write(iw,*) 'ADC2 h/h2p chksum (c) in symmetry',CREP,
     &        zdotc(n,bufhp,1,bufhp,1)
      else
        write(iw,*) 'ADC2 h/h2p chksum (r) in symmetry',CREP,
     &        ddot(n,bufhp,1,bufhp,1)
      endif


      RETURN 
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKEHPHP(CKKS,DIAG,BUF,EAJL,OOOO,VOVO,IOI,IOJ,ICRA,
     &                    OOLT,LADC,CREP,MXNO,NBUFS,IB,INTBUF,
     &                    ADCLEVEL,ITAPE1,ITAPE2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates 2h-1p/2h-1p contribution to the ADC matrix.
C     Here the treatment of ADC(2), and ADC(3) is combined due to
C     programming clarity.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 CKKS(*),DIAG(*),BUF(*),EAJL(*),OOOO(*),VOVO(*)
      INTEGER LADC,CREP,MXNO,NBUFS,IB,INTBUF,ADCLEVEL,ITAPE1,ITAPE2
      INTEGER IOI(*),IOJ(*),ICRA(6,*),OOLT(MXNO,MXNO,*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 S,KCOCC
      INTEGER JDUMMY
      INTEGER zeros
C
C---------------Executable Code --------------------------------------
C
      CALL PST('Constructing 2h1p/2h1p block (real)+')
C
      JDUMMY = 0
      NOK=NO(CREP)
      NAJL=NVOOT(CREP)

C
C  fill in the diagonal of the h/h block ...
C
      OFF1 = 1
      DO K=1,NOK
        DIAG(K)=CKKS(OFF1 + K - 1)
        OFF1 = OFF1 + NOK
      ENDDO
C
C  ... and the elements K_akl
C  (this is the only contribution in case of ADC-2 strict)
C
      CALL XCOPY(NAJL,EAJL,1,DIAG(NOK+1),1)
C
C  enter additional code for off-diagonal matrix elements in ADC(2)X and
C  ADC(3) calculations. Only in this case the arrays OOOO and VOVO are
C  addressed.
C
      IF(ADCLEVEL.GE.2) THEN
        CALL GETOOOO(OOOO)
        CALL GETVOVO(VOVO)
C
C  construct OOT lookup table for the IRREPS
C
        DO IREP=1,NREP
        DO I=1,MXNO
        DO J=1,MXNO
          OOLT(J,I,IREP)=0
        ENDDO
        ENDDO
        ENDDO
        DO IREP=1,NREP
          IOFF=1
          DO J=1,NO(IREP)
            DO I=J+1,NO(IREP)
              OOLT(I,J,IREP)=IOFF
              IOFF=IOFF+1
!              WRITE(iw,'(1X,4(A7,I4))') 'OOLT = ',ioff-1,'IREP = ',IREP,
!     &                                  'j = ', j, 'i = ', i
            ENDDO
          ENDDO
        ENDDO

!        WRITE(iw,*) 'size of oolt = ', SIZE(oolt)
C
C  loop over matrix entries and write to buffer immediately
C

        zeros = 0
        DO J=1,NVOOT(CREP)
          DO I=J,NVOOT(CREP)
  
            IROW=NOK+I
            ICOL=NOK+J

            S=KCOCC(I,J,OOOO,VOVO,ICRA,OOLT,MXNO)
            !write(iw,'(1X,A4,ES15.4)') 's = ', S
            IF (S.EQ. 0.0D0) THEN
              zeros = zeros + 1
            END IF
C
C  store matrix element either in the array for diagonal elements
C  or in the write buffer  but then as the negative value ! all diagonal
C  elements will be negated further below.
C
            IF(IROW.EQ.ICOL) THEN
              DIAG(IROW) = DIAG(IROW) + S
            ELSE
              IF(S.NE.0.0D0) THEN
                IB=IB+1
                BUF(IB)=S
                IOI(IB)=IROW
                IOJ(IB)=ICOL
                IF(IB.EQ.INTBUF) THEN
                  NBUFS = NBUFS + 1
                  WRITE(ITAPE2) (-BUF(IXX),IXX=1,INTBUF),
     &                           (IOI(IXX),IXX=1,INTBUF),
     &                           (IOJ(IXX),IXX=1,INTBUF),
     &                           INTBUF,JDUMMY
                  IB = 0
                ENDIF
              ENDIF
            ENDIF

          ENDDO
        ENDDO

        WRITE(iw,*) 'Number of zeros = ', zeros

      ENDIF
C
C  write out remaining elements in matrix buffer
C
      IF(IB.GT.0) THEN
        NBUFS = NBUFS + 1
        WRITE(ITAPE2) (-BUF(IXX),IXX=1,INTBUF),
     &                 (IOI(IXX),IXX=1,INTBUF),
     &                 (IOJ(IXX),IXX=1,INTBUF),
     &                 IB,JDUMMY
        IB = 0
      ENDIF
C
C  write out !negative! diagonal in specific file
C  and the number of records in the off-diagonal file.
C
!      WRITE(iw,*) 'ADC diag' !Elke
!      WRITE(iw,'(1X,ES14.5)') (diag(ixx),ixx=1,ladc)
      WRITE(ITAPE1) (-DIAG(IXX),IXX=1,LADC),NBUFS
      WRITE(IW,1000) INTBUF,NBUFS
C
C  files will be closed in main program
C
 1000 FORMAT (1X,'Length of write buffers:',T30,I7/
     &        1X,'Number of buffers written:',T30,I6)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMAKEHPHP(CKKS,DIAG,BUF,EAJL,OOOO,VOVO,IOI,IOJ,ICRA,
     &                    OOLT,LADC,CREP,MXNO,NBUFS,IB,INTBUF,
     &                    ADCLEVEL,ITAPE1,ITAPE2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates 2h-1p/2h-1p contribution to the ADC matrix.
C     Here the treatment of ADC(2), and ADC(3) is combined due to
C     programming clarity.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EAJL(*)
      COMPLEX*16 CKKS(*),DIAG(*),BUF(*),OOOO(*),VOVO(*)
      INTEGER LADC,CREP,MXNO,NBUFS,IB,INTBUF,ADCLEVEL,ITAPE1,ITAPE2
      INTEGER IOI(*),IOJ(*),ICRA(6,*),OOLT(MXNO,MXNO,*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER JDUMMY,TOTELE
      COMPLEX*16 S,CKCOCC
C
C---------------Executable Code --------------------------------------
C
      CALL PST('Constructing 2h1p/2h1p block (complex)+')
C
      JDUMMY = 0
      NOK=NO(CREP)
      NAJL=NVOOT(CREP)
C
C  fill in the diagonal of the h/h block ...
C
      OFF1 = 1
      DO K=1,NOK
        DIAG(K)=CKKS(OFF1 + K - 1)
        OFF1 = OFF1 + NOK
      ENDDO
C
C  ... and the elements K_akl
C  (this is the only contribution in case of ADC-2 strict)
C
      DO K=1,NAJL
        DIAG(NOK+K) = DCMPLX(EAJL(K),0.0D0)
      ENDDO
C
C  enter additional code for off-diagonal matrix elements in ADC(2)X and
C  ADC(3) calculations. Only in this case the arrays OOOO and VOVO are
C  addressed.
C
      IF(ADCLEVEL.GE.2) THEN
        CALL GETOOOO(OOOO)
        CALL GETVOVO(VOVO)
C
C  construct OOT lookup table for the IRREPS
C
        DO IREP=1,NREP
        DO I=1,MXNO
        DO J=1,MXNO
          OOLT(J,I,IREP)=0
        ENDDO
        ENDDO
        ENDDO
        DO IREP=1,NREP
          IOFF=1
          DO J=1,NO(IREP)
            DO I=J+1,NO(IREP)
              OOLT(I,J,IREP)=IOFF
              IOFF=IOFF+1
            ENDDO
          ENDDO
        ENDDO
C
C  determine additional diagonal and off-diagonal elements
C  ATT* the CKCOCC routine yields matrix entries in the
C  <||>|>,|><|<| convention
C

        DO J=1,NVOOT(CREP)
          DO I=J,NVOOT(CREP)
            IROW=NOK+I
            ICOL=NOK+J
            S=CKCOCC(I,J,OOOO,VOVO,ICRA,OOLT,MXNO)

C
C  store matrix element either in the array for diagonal elements
C  or in the write buffer  but then as the negative value ! all diagonal
C  elements will be negated further below.
C
            IF(IROW.EQ.ICOL) THEN
              DIAG(IROW) = DIAG(IROW) + S
            ELSE
              IF(S.NE.(0.0D0,0.0D0)) THEN
                IB=IB+1
                BUF(IB)=S
                IOI(IB)=IROW
                IOJ(IB)=ICOL
                IF(IB.EQ.INTBUF) THEN
                  NBUFS = NBUFS + 1
                  WRITE(ITAPE2) (-BUF(IXX),IXX=1,INTBUF),
     &                           (IOI(IXX),IXX=1,INTBUF),
     &                           (IOJ(IXX),IXX=1,INTBUF),
     &                           INTBUF,JDUMMY
                  IB = 0
                ENDIF
              ENDIF
            ENDIF

          ENDDO
        ENDDO

      ENDIF
C
C  write out remaining elements in matrix buffer
C  and calculate total number of elements
C  this is possible here because that's the last access to the
C  off-dia file.
C
      TOTELE = 0
      IF(NBUFS.GE.1) THEN
        TOTELE = NBUFS * INTBUF
      ENDIF
      IF(IB.GT.0) THEN
        NBUFS = NBUFS + 1
        WRITE(ITAPE2) (-BUF(IXX),IXX=1,INTBUF),
     &                 (IOI(IXX),IXX=1,INTBUF),
     &                 (IOJ(IXX),IXX=1,INTBUF),
     &                 IB,JDUMMY
        TOTELE = TOTELE + IB
        IB = 0
      ENDIF
C
C  write out !negative! diagonal in specific file
C  and the number of records in the off-diagonal file.
C
      WRITE(ITAPE1) (-DIAG(IXX),IXX=1,LADC),NBUFS
      WRITE(IW,1000) INTBUF,NBUFS,TOTELE
C
C  files will be closed in main program
C
 1000 FORMAT (1X,'Length of write buffers:',T30,I7/
     &        1X,'Number of buffers written:',T30,I6/
     &        1X,'Total number of elements in Offdiag-file:',T45,I8)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION KCOCC(IROW,JCOL,OOOO,VOVO,ICRA,OOLT,MXNO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     This function returns the C(1)_akl,a'k'l' element of the 2h1p/2h1p
C     block in the (N-1) space. The integrals OOOO and VOVO have to be
C     provided and ICRA/OOLT contain translation/lookup tables referring
C     to the irrep under consideration (generated by the caller).
C     The matrix is of dimension NVOOT x NVOOT and the proper range
C     of row I and column J has to be maintained by the caller !
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER IROW,JCOL,MXNO,ICRA(6,*),OOLT(MXNO,MXNO,*)
      REAL*8 OOOO(*),VOVO(*)
      REAL*8 KCOCC
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 S
C
C---------------Executable code--------------------------------------
C
      AS    = ICRA(1,JCOL)
      ASREP = ICRA(2,JCOL)
      KS    = ICRA(3,JCOL)
      KSREP = ICRA(4,JCOL)
      LS    = ICRA(5,JCOL)
      LSREP = ICRA(6,JCOL)

      A     = ICRA(1,IROW)
      AREP  = ICRA(2,IROW)
      K     = ICRA(3,IROW)
      KREP  = ICRA(4,IROW)
      L     = ICRA(5,IROW)
      LREP  = ICRA(6,IROW)


      S=0.0D0

C
C   ..... Part A
C
      IF(AREP.EQ.ASREP) THEN
        IF(A.EQ.AS) THEN
          KLREP=MULTB(KREP,LREP,1)
          IF(KREP.EQ.LREP) THEN
            IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
          ELSE
            IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
          ENDIF
          IOFF=IOFF * NOOT(KLREP) + IOOOOTT(KLREP)
          IF(KSREP.EQ.LSREP) THEN
            IOFF=IOFF + IIOOT(KSREP,LSREP) + OOLT(KS,LS,KSREP)
          ELSE
            IOFF=IOFF + IIOOT(KSREP,LSREP) +
     &        (LS-1)*NO(KSREP) + KS
          ENDIF
          S = S - OOOO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part B
C
      IF(KREP.EQ.KSREP) THEN
        IF(K.EQ.KS) THEN
          ASLREP=MULTB(ASREP,LREP,1)
          IOFF=IIVO(ASREP,LREP) + (L-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASLREP) + IVOVO(ASLREP)
          IOFF=IOFF + IIVO(AREP,LSREP) +
     &        (LS-1)*NV(AREP) + A
          S = S + VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part C
C
      IF(LREP.EQ.LSREP) THEN
        IF(L.EQ.LS) THEN
          ASKREP=MULTB(ASREP,KREP,1)
          IOFF=IIVO(ASREP,KREP) + (K-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASKREP) + IVOVO(ASKREP)
          IOFF=IOFF + IIVO(AREP,KSREP) +
     &        (KS-1)*NV(AREP) + A
          S = S + VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part D
C
      IF(LREP.EQ.KSREP) THEN
        IF(L.EQ.KS) THEN
          ASKREP=MULTB(ASREP,KREP,1)
          IOFF=IIVO(ASREP,KREP) + (K-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASKREP) + IVOVO(ASKREP)
          IOFF=IOFF + IIVO(AREP,LSREP) +
     &        (LS-1)*NV(AREP) + A
          S = S - VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part E
C
      IF(KREP.EQ.LSREP) THEN
        IF(K.EQ.LS) THEN
          ASLREP=MULTB(ASREP,LREP,1)
          IOFF=IIVO(ASREP,LREP) + (L-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASLREP) + IVOVO(ASLREP)
          IOFF=IOFF + IIVO(AREP,KSREP) +
     &        (KS-1)*NV(AREP) + A
          S = S - VOVO(IOFF)
        ENDIF
      ENDIF

      KCOCC = S

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION KCVIR(IROW,JCOL,VVVV,VOVO,ICRB,VVLT,MXNV,EPS,ECC)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     This function returns a C(1)_jab,j'a'b' element of the 2p1h/2p1h
C     block in the (N+1) space. The integrals VVVV and VOVO have to be
C     provided and ICRB/VVLT contain translation/lookup tables referring
C     to the irrep under consideration (generated by the caller).
C     The matrix is of dimension NOVVT x NOVVT and the proper range
C     of row I and column J has to be maintained by the caller !
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER IROW,JCOL,MXNV,ICRB(6,*),VVLT(MXNV,MXNV,*)
      REAL*8 VVVV(*),VOVO(*)
      REAL*8 KCVIR,EPS(*)
      REAL*8 ECC
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 S
C
C---------------Executable code--------------------------------------
C
      JS    = ICRB(1,JCOL)
      JSREP = ICRB(2,JCOL)
      AS    = ICRB(3,JCOL)
      ASREP = ICRB(4,JCOL)
      BS    = ICRB(5,JCOL)
      BSREP = ICRB(6,JCOL)

      J     = ICRB(1,IROW)
      JREP  = ICRB(2,IROW)
      A     = ICRB(3,IROW)
      AREP  = ICRB(4,IROW)
      B     = ICRB(5,IROW)
      BREP  = ICRB(6,IROW)

      S=0.0D0
C
C   ..... Part A
C
      IF(JREP.EQ.JSREP) THEN
        IF(J.EQ.JS) THEN
          ASBSREP=MULTB(ASREP,BSREP,1)
          IF(ASREP.EQ.BSREP) THEN
            IOFF=IIVVT(ASREP,BSREP) + VVLT(AS,BS,ASREP) - 1
          ELSE
            IOFF=IIVVT(ASREP,BSREP) + (BS-1)*NV(ASREP) + AS - 1
          ENDIF
          IOFF=IOFF * NVVT(ASBSREP) + IVVVVTT(ASBSREP)
          IF(AREP.EQ.BREP) THEN
            IOFF=IOFF + IIVVT(AREP,BREP) + VVLT(A,B,AREP)
          ELSE
            IOFF=IOFF + IIVVT(AREP,BREP) +
     &        (B-1)*NV(AREP) + A
          ENDIF
          S = S + VVVV(IOFF)
          if(ioff.gt.ivvvvtt(nrep+1)) stop 'vvvv offset error !'
        ENDIF
      ENDIF
C
C   ..... Part B
C
      IF(AREP.EQ.ASREP) THEN
        IF(A.EQ.AS) THEN
          BSJREP=MULTB(BSREP,JREP,1)
          IOFF=IIVO(BSREP,JREP) + (J-1)*NV(BSREP) + BS - 1
          IOFF=IOFF*NVO(BSJREP) + IVOVO(BSJREP)
          IOFF=IOFF + IIVO(BREP,JSREP) +
     &        (JS-1)*NV(BREP) + B
          S = S - VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part C
C
      IF(BREP.EQ.BSREP) THEN
        IF(B.EQ.BS) THEN
          ASJREP=MULTB(ASREP,JREP,1)
          IOFF=IIVO(ASREP,JREP) + (J-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASJREP) + IVOVO(ASJREP)
          IOFF=IOFF + IIVO(AREP,JSREP) +
     &        (JS-1)*NV(AREP) + A
          S = S - VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part D
C
      IF(AREP.EQ.BSREP) THEN
        IF(A.EQ.BS) THEN
          ASJREP=MULTB(ASREP,JREP,1)
          IOFF=IIVO(ASREP,JREP) + (J-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASJREP) + IVOVO(ASJREP)
          IOFF=IOFF + IIVO(BREP,JSREP) +
     &        (JS-1)*NV(BREP) + B
          S = S + VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part E
C
      IF(BREP.EQ.ASREP) THEN
        IF(B.EQ.AS) THEN
          BSJREP=MULTB(BSREP,JREP,1)
          IOFF=IIVO(BSREP,JREP) + (J-1)*NV(BSREP) + BS - 1
          IOFF=IOFF*NVO(BSJREP) + IVOVO(BSJREP)
          IOFF=IOFF + IIVO(AREP,JSREP) +
     &        (JS-1)*NV(AREP) + A
          S = S + VOVO(IOFF)
        ENDIF
      ENDIF

      KCVIR = S

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION CKCOCC(IROW,JCOL,OOOO,VOVO,ICRA,OOLT,MXNO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     This function returns the C(1)_akl,a'k'l' element of the 2h1p/2h1p
C
C  ATT ! This implementation complies to the BKC notation: 
C               C(1)_akl,a'k'l'=<||>|>,|><|<|
C
C     The matrix is of dimension NVOOT x NVOOT and the proper range
C     of row I and column J has to be maintained by the caller !
C
C  ATT ! One MUST NOT change this by CC the final result because then
C     it will not fit anymore in the chosen overall ADC matrix BKC
C     convention !!
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER IROW,JCOL,MXNO,ICRA(6,*),OOLT(MXNO,MXNO,*)
      COMPLEX*16 OOOO(*),VOVO(*)
      COMPLEX*16 CKCOCC
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 S
C
C---------------Executable code--------------------------------------
C
      AS    = ICRA(1,JCOL)
      ASREP = ICRA(2,JCOL)
      KS    = ICRA(3,JCOL)
      KSREP = ICRA(4,JCOL)
      LS    = ICRA(5,JCOL)
      LSREP = ICRA(6,JCOL)

      A     = ICRA(1,IROW)
      AREP  = ICRA(2,IROW)
      K     = ICRA(3,IROW)
      KREP  = ICRA(4,IROW)
      L     = ICRA(5,IROW)
      LREP  = ICRA(6,IROW)

      S=(0.0D0,0.0D0)

C
C   ..... Part A
C
      IF(AREP.EQ.ASREP) THEN
        IF(A.EQ.AS) THEN
          KLREP=MULTB(KREP,LREP,1)
          IF(KREP.EQ.LREP) THEN
            IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
          ELSE
            IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
          ENDIF
          IOFF=IOFF * NOOT(KLREP) + IOOOOTT(KLREP)
          IF(KSREP.EQ.LSREP) THEN
            IOFF=IOFF + IIOOT(KSREP,LSREP) + OOLT(KS,LS,KSREP)
          ELSE
            IOFF=IOFF + IIOOT(KSREP,LSREP) +
     &        (LS-1)*NO(KSREP) + KS
          ENDIF
          S = S - OOOO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part B
C
      IF(KREP.EQ.KSREP) THEN
        IF(K.EQ.KS) THEN
          ASLREP=MULTB(ASREP,LREP,1)
          IOFF=IIVO(ASREP,LREP) + (L-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASLREP) + IVOVO(ASLREP)
          IOFF=IOFF + IIVO(AREP,LSREP) +
     &        (LS-1)*NV(AREP) + A
          S = S + VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part C
C
      IF(LREP.EQ.LSREP) THEN
        IF(L.EQ.LS) THEN
          ASKREP=MULTB(ASREP,KREP,1)
          IOFF=IIVO(ASREP,KREP) + (K-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASKREP) + IVOVO(ASKREP)
          IOFF=IOFF + IIVO(AREP,KSREP) +
     &        (KS-1)*NV(AREP) + A
          S = S + VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part D
C
      IF(LREP.EQ.KSREP) THEN
        IF(L.EQ.KS) THEN
          ASKREP=MULTB(ASREP,KREP,1)
          IOFF=IIVO(ASREP,KREP) + (K-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASKREP) + IVOVO(ASKREP)
          IOFF=IOFF + IIVO(AREP,LSREP) +
     &        (LS-1)*NV(AREP) + A
          S = S - VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part E
C
      IF(KREP.EQ.LSREP) THEN
        IF(K.EQ.LS) THEN
          ASLREP=MULTB(ASREP,LREP,1)
          IOFF=IIVO(ASREP,LREP) + (L-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASLREP) + IVOVO(ASLREP)
          IOFF=IOFF + IIVO(AREP,KSREP) +
     &        (KS-1)*NV(AREP) + A
          S = S - VOVO(IOFF)
        ENDIF
      ENDIF

c  ** no cc here !

      CKCOCC = S

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION CKCVIR(IROW,JCOL,VVVV,VOVO,ICRB,VVLT,MXNV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     This function returns a C(1)_jab,j'a'b' element of the 2p1h/2p1h
C
C     ATT***!!! This implementation complies to the OSWALD convention
C     and yields: C(1)_jab,j'a'b' = |><|<|,<||>|>
C     The integrals VVVV and VOVO have to be
C     provided and ICRB/VVLT contain translation/lookup tables referring
C     to the irrep under consideration (generated by the caller).
C     The matrix is of dimension NOVVT x NOVVT and the proper range
C     of row I and column J has to be maintained by the caller !
C
C     Due to the **Oswald** convention the value of this function is complex
C     conjugated when used by the VAFFITER routine in order to comply
C     with the overall BKC !
C
C     By convention ALWAYS the lower triangle is accessed. This is also
C     because of the integral storage convention in MDCINT.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER IROW,JCOL,MXNV,ICRB(6,*),VVLT(MXNV,MXNV,*)
      COMPLEX*16 VVVV(*),VOVO(*)
      COMPLEX*16 CKCVIR
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 S
C
C---------------Executable code--------------------------------------
C
      JS    = ICRB(1,JCOL)
      JSREP = ICRB(2,JCOL)
      AS    = ICRB(3,JCOL)
      ASREP = ICRB(4,JCOL)
      BS    = ICRB(5,JCOL)
      BSREP = ICRB(6,JCOL)

      J     = ICRB(1,IROW)
      JREP  = ICRB(2,IROW)
      A     = ICRB(3,IROW)
      AREP  = ICRB(4,IROW)
      B     = ICRB(5,IROW)
      BREP  = ICRB(6,IROW)

      S=(0.0D0,0.0D0)
C
C   ..... Part A
C
      IF(JREP.EQ.JSREP) THEN
        IF(J.EQ.JS) THEN
          ASBSREP=MULTB(ASREP,BSREP,1)
          IF(ASREP.EQ.BSREP) THEN
            IOFF=IIVVT(ASREP,BSREP) + VVLT(AS,BS,ASREP) - 1
          ELSE
            IOFF=IIVVT(ASREP,BSREP) + (BS-1)*NV(ASREP) + AS - 1
          ENDIF
          IOFF=IOFF * NVVT(ASBSREP) + IVVVVTT(ASBSREP)
          IF(AREP.EQ.BREP) THEN
            IOFF=IOFF + IIVVT(AREP,BREP) + VVLT(A,B,AREP)
          ELSE
            IOFF=IOFF + IIVVT(AREP,BREP) +
     &        (B-1)*NV(AREP) + A
          ENDIF
          S = S + VVVV(IOFF)
          if(ioff.gt.ivvvvtt(nrep+1))
     &         stop 'vvvv offset error in CKCVIR!'
        ENDIF
      ENDIF
C
C   ..... Part B
C
      IF(AREP.EQ.ASREP) THEN
        IF(A.EQ.AS) THEN
          BSJREP=MULTB(BSREP,JREP,1)
          IOFF=IIVO(BSREP,JREP) + (J-1)*NV(BSREP) + BS - 1
          IOFF=IOFF*NVO(BSJREP) + IVOVO(BSJREP)
          IOFF=IOFF + IIVO(BREP,JSREP) +
     &        (JS-1)*NV(BREP) + B
          S = S - VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part C
C
      IF(BREP.EQ.BSREP) THEN
        IF(B.EQ.BS) THEN
          ASJREP=MULTB(ASREP,JREP,1)
          IOFF=IIVO(ASREP,JREP) + (J-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASJREP) + IVOVO(ASJREP)
          IOFF=IOFF + IIVO(AREP,JSREP) +
     &        (JS-1)*NV(AREP) + A
          S = S - VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part D
C
      IF(AREP.EQ.BSREP) THEN
        IF(A.EQ.BS) THEN
          ASJREP=MULTB(ASREP,JREP,1)
          IOFF=IIVO(ASREP,JREP) + (J-1)*NV(ASREP) + AS - 1
          IOFF=IOFF*NVO(ASJREP) + IVOVO(ASJREP)
          IOFF=IOFF + IIVO(BREP,JSREP) +
     &        (JS-1)*NV(BREP) + B
          S = S + VOVO(IOFF)
        ENDIF
      ENDIF
C
C   ..... Part E
C
      IF(BREP.EQ.ASREP) THEN
        IF(B.EQ.AS) THEN
          BSJREP=MULTB(BSREP,JREP,1)
          IOFF=IIVO(BSREP,JREP) + (J-1)*NV(BSREP) + BS - 1
          IOFF=IOFF*NVO(BSJREP) + IVOVO(BSJREP)
          IOFF=IOFF + IIVO(AREP,JSREP) +
     &        (JS-1)*NV(AREP) + A
          S = S + VOVO(IOFF)
        ENDIF
      ENDIF

      CKCVIR = S

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CKKASYM(HH,CHH,NSIZE,ASYMHH)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes the maximum deviation from symmetry (hermiticity) 
C     in the given matrix given as the matrix HH or CHH.
C     Dimension of the sqare matrix is NSIZE
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 HH(*)
      COMPLEX*16 CHH(*)
      INTEGER NSIZE
      REAL*8 ASYMHH
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/complex.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER N
      REAL*8 RV
      COMPLEX*16 CV
C
C---------------Executable code--------------------------------------
C
      ASYMHH=0.0D0
      OFFN = 1
      OFFT = 1

      N=NSIZE
      DO KS=1,N
      DO K=1,KS
        OFFNA = OFFN + (KS-1)*N + K  - 1
        OFFTA = OFFT + (K-1)*N  + KS - 1
        IF(CARITH) THEN
          CV = DCONJG(CHH(OFFNA)) - CHH(OFFTA)
#if defined (BIT64)
          ASYMHH = ASYMHH + DBLE(CONJG(CV)*CV)
#else
          ASYMHH = ASYMHH + DBLE(DCONJG(CV)*CV)
#endif
        ELSE
          RV = HH(OFFNA) - HH(OFFTA)
          ASYMHH = ASYMHH + RV*RV
        ENDIF
      ENDDO
      ENDDO

#if defined (BIT64)
      ASYMHH =  SQRT(ASYMHH)
#else
      ASYMHH =  DSQRT(ASYMHH)
#endif

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMVVOO (EPS,T2,CT2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Divide VVOO integrals by denominators ea + eb - ej - ei
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher, adapted by MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 T2(*)
      COMPLEX*16 CT2(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 FAC,FAC1,FAC2,FAC3
      COMPLEX*16 CFAC
C
C---------------Executable code--------------------------------------
C
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = - EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = - EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 + EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 + EPS(AA)
                  ABIJ = ABIJ + 1
                  IF (CARITH) THEN
                     CFAC = DCMPLX(FAC,0.0D0)
                     CT2(ABIJ) = CT2(ABIJ)/CFAC
                  ELSE
                     T2(ABIJ) = T2(ABIJ)/FAC
                  ENDIF
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      IF(ABIJ.ne.IVVOOTT(nrep+1)) THEN
        WRITE(*,*) 'ABIJ/IVVOOTT:',ABIJ,IVVOOTT(NREP+1)
        STOP 'COUNTING ERROR IN DENOMVVOO!'
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EPSARR(EPS,BUF,KREP,ISIZE)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Constructs the epsilon array E(VVOT,O,O)
C     Energies always real ==> this version works for real/complex ADC
C     Matrices
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 BUF(*)
      INTEGER KREP
      INTEGER ISIZE
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 EPSK,EPSKS,EPSL,EPSB,EPSA
C
C---------------Executable code--------------------------------------
C
      DES = 1
      KOFF=IO(KREP)
      DO 55 K=1,NO(KREP)
       EPSK=EPS(KOFF+K)*(-0.5D0)
       DO 50 KS=1,NO(KREP)
        EPSKS=EPS(KOFF+KS)*(-0.5D0)
        DO 45 LREP=1,NREP
         LOFF=IO(LREP)
         ABREP=MULTB(LREP,KREP,1)

C from here ABREP(b) and LREP(f) are known.

         DO LFIE = 1,NO(LREP) 
          EPSL=-EPS(LOFF+LFIE)
          DO 40 BREP=1,NREP
           AREP=MULTB(BREP,ABREP+NREP,2)
           IF(AREP.LT.BREP) GOTO 40
           BOFF=IV(BREP) + IO(NREP+1)
           AOFF=IV(AREP) + IO(NREP+1)
           DO 30 BFIE=1,NV(BREP)
             EPSB=EPS(BOFF+BFIE)
             AMIN=1
             IF(AREP.EQ.BREP) AMIN=BFIE+1
             DO 25 AFIE=AMIN,NV(AREP)
              EPSA=EPS(AOFF+AFIE)
              BUF(DES)=EPSA+EPSB+EPSL+EPSK+EPSKS
              DES=DES+1
 25          CONTINUE
 30         CONTINUE

  40       CONTINUE
         ENDDO                
 45     CONTINUE
 50    CONTINUE
 55   CONTINUE

      ISIZE=DES-1

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CNTRCT3(V1,V2,VE,CKKS,M,N,ELOFFSET)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Contract the V1,V2 and VE array over a common index range
C     dependent on the leading dimension in M. Hereby
C     V1,V2 and CKKS can be complex, VE is always real.
C     OFFKREP :  startoffset into the v_abl,k array since we start with 
C             :  irrep KREP in the main program
C            
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 V1(*),V2(*),VE(*),CKKS(*)
      INTEGER M,N,ELOFFSET
c
c  ATT: arrays declared as real but there can be complex numbers in them
c  therefore we have to stretch the offset by a factor of RCW !
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/complex.inc"
C
C---------------Executable code--------------------------------------
C
      OFF4=1
      OFFE=1
      OFFKREP = ELOFFSET*RCW

      OFF2=OFFKREP + 1
      DO 20 KS=1,N
       OFF1=OFFKREP + 1
       DO 10 K=1,N
         IF(CARITH) THEN
          CALL CCONTR3(V1(OFF1),V2(OFF2),VE(OFFE),CKKS(OFF4),M)
         ELSE
          CALL RCONTR3(V1(OFF1),V2(OFF2),VE(OFFE),CKKS(OFF4),M)
         ENDIF
         OFF1 = OFF1 + M * RCW
         OFFE = OFFE + M
         OFF4 = OFF4 + RCW
 10    CONTINUE
       OFF2 = OFF2 + M * RCW
 20   CONTINUE

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCONTR3 (A,B,C,D,M)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Performs a COMPLEX A,B,C three-array contraction of length M
C     and adds the NEGATIVE result as the FIRST element to COMPLEX array D
C     adding is essential since D could already contain the
C     Sigma(infinity) contributions !
C     The calling routine ensures the correct offsets for all arrays !
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 A(*),B(*),D(*)
      REAL*8 C(*)
      INTEGER M
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 SUM
C
C---------------Executable code--------------------------------------
C
      SUM = (0.0D0,0.0D0)

      DO I=1,M
        SUM = SUM + A(I) * B(I) * DCMPLX(C(I),0.0D0)
      ENDDO

      D(1) = D(1) - SUM

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RCONTR3 (A,B,C,D,M)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Performs a REAL A,B,C three-array contraction of length M
C     and adds the NEGATIVE result as the FIRST element to REAL array D.
C     adding is essential since D could already contain the
C     Sigma(infinity) contributions !
C     The calling routine ensures the correct offsets for all arrays !
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 A(*),B(*),C(*),D(*)
      INTEGER M
C
C---------------Local variables--------------------------------------
C
#include  "../relccsd/files.inc"
      REAL*8 SUM
C
C---------------Executable code--------------------------------------
C
      SUM = 0.0D0

      DO I=1,M
        SUM = SUM + A(I) * B(I) * C(I)
      ENDDO

      D(1) = D(1) - SUM

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADDONE(A,AC,NR,NC)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Adds unit matrix to the quadratic 
C     real/complex matrix A
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER NR,NC
      REAL*8 A(NR,NC)
      COMPLEX*16 AC(NR,NC)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER N
C
C---------------Executable code--------------------------------------
C
      IF(NR.NE.NC) CALL QUIT('ADDONE: Matrix not quadratic !')

      N=NR

      DO I=1,N
        IF(CARITH) THEN
          AC(I,I) = AC(I,I) + (1.0D0,0.0D0)
        ELSE
          A(I,I) = A(I,I) + 1.0D0
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADDEKK(EPS,CKK,CCKK,KREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Adds the occupied spinor energies on the diagonal to C_kk(2).
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 CKK(*)
      COMPLEX*16 CCKK(*),CEPS
      INTEGER KREP
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/complex.inc"
#include  "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
C
C---------------Executable code--------------------------------------
C
      OFFA = 1
      N=NO(KREP)

      DO K=1,N
        OFFS=IO(KREP)+K
        IF(CARITH) THEN
          CEPS=DCMPLX(EPS(OFFS),0.0D0)
          CCKK(OFFA) = CCKK(OFFA) + CEPS
        ELSE
          CKK(OFFA) = CKK(OFFA) + EPS(OFFS)
        ENDIF
        OFFA = OFFA + N + 1
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADDARR(N,A1,CA1,A2,CA2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Adds NxN array A1 (or complex array CA1) to NxN array A2 (CA2 resp.)
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER N
      REAL*8 A1(*),A2(*)
      COMPLEX*16 CA1(*),CA2(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
C
C---------------Executable code--------------------------------------
C
      DO I=1,N*N
        IF(CARITH) THEN
          CA2(I) = CA2(I) + CA1(I)
        ELSE
          A2(I)  = A2(I)  + A1(I)
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HERMCON(N,A,CA)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Forms transpose or hermitian conjugate of the 
C     NxN matrix A1 (or complex array CA1)
C     In the complex case the diagonal will also be conjugated
C     because the routine is called for expressions with a not
C     necessarily vanishing imaginary part on the diagonal due to the
C     asymmetry of the expression with respect to k <==> k' interchange.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER N
      REAL*8 A(N,*)
      COMPLEX*16 CA(N,*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 AX
      COMPLEX*16 CAX
C
C---------------Executable code--------------------------------------
C
      DO J=1,N
        DO I=J+1,N
            IF(CARITH) THEN
            CAX = CA(I,J)
            CA(I,J) = DCONJG(CA(J,I))
            CA(J,I) = DCONJG(CAX)
          ELSE
            AX = A(I,J)
            A(I,J) = A(J,I)
            A(J,I) = AX
          ENDIF
        ENDDO
        IF(CARITH) CA(J,J)=DCONJG(CA(J,J))
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKEEKVOO(E,EPS,AKLREP,ICRA)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Creates the packed energy array -Ea+Ej+El according to the
C     implicit irrep and orbital order as in a NVOOT array for the
C     requested irrep KREP. At the same time it creates the 1p-2h
C     configuration resolution array in this order !
C     structure of ICRA array: #p irep(p) #h (irep h) #h irep(h)
C
C     ATT! energy arrays always real ==> we have to stick to DCOPY !
C     This also holds in the case of CARITH=.TRUE. !!
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 E(*),EPS(*)
      INTEGER AKLREP,ICRA(6,*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/param.inc"
#include "../relccsd/symm.inc"
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 EPSA,EPSK,EPSL
C
C---------------Executable code--------------------------------------
C
C  clear destination array
C
      LEN=NVOOT(AKLREP)
      CALL DCOPY(LEN,AR0,0,E,1)
C
C  loop over irrep structure and form energy array & lookup table
C
      !WRITE(iw,'(1X,3A14)') 'epsa', 'epsk', 'epsl'
      OFF=1
      DO 60 KLREP=1,NREP
        AREP=MULTB(AKLREP,KLREP+NREP,2)
        AOFF=IV(AREP) + IO(NREP+1)
        DO 45 LREP=1,NREP
         LOFF=IO(LREP)
         KREP=MULTB(LREP,KLREP+NREP,2)
!        IF(KREP.LT.LREP) GOTO45
         IF(KREP.LT.LREP) GOTO 45
         KOFF=IO(KREP)
         DO L=1,NO(LREP)
           EPSL=EPS(LOFF+L)
           KMIN=1
           IF(KREP.EQ.LREP) KMIN = L + 1
           DO K=KMIN,NO(KREP)
             EPSK=EPS(KOFF+K)
             DO A=1,NV(AREP)
               EPSA=EPS(AOFF+A)
               E(OFF)=-EPSA+EPSK+EPSL
               ICRA(1,OFF)=A
               ICRA(2,OFF)=AREP
               ICRA(3,OFF)=K
               ICRA(4,OFF)=KREP
               ICRA(5,OFF)=L
               ICRA(6,OFF)=LREP
               !WRITE(iw,'(1X,3ES14.5)') epsa, epsk, epsl
               OFF=OFF+1
             ENDDO
           ENDDO
         ENDDO
! 45     ENDDO
 45      CONTINUE
!60   ENDDO
 60   CONTINUE

      IF( (OFF-1).NE.NVOOT(AKLREP)) THEN
        WRITE(IW,*) ' *** Expected value of counter:',NVOOT(AKLREP)
        WRITE(IW,*) ' *** Actual value of counter:  ',OFF-1
        CALL QUIT('Inconsistency in MAKEEKVOO')
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKEEKOVV(E,EPS,JABREP,ICRB)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Creates the packed energy array -Ej+Ea+Eb according to the
C     implicit irrep and orbital order as in a NOVVT array for the
C     requested irrep KREP. At the same time it creates the 
C     corresponding lookup table for the OVV <> number translation
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 E(*),EPS(*)
      INTEGER JABREP,ICRB(6,*)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/param.inc"
#include "../relccsd/symm.inc"
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 EPSJ,EPSA,EPSB
C
C---------------Executable code--------------------------------------
C
C ************** ATT ! this array is ALWAYS REAL !!!
C that means if we are in complex arithmetic, usage of XCOPY is
C erroneous !!!!!!!! we have to stick to DCOPY
C
C  clear destination array
C
      LEN=NOVVT(JABREP)
      CALL DCOPY(LEN,AR0,0,E,1)
C
C  loop over irrep structure and form energy array & lookup table
C
      OFF=0
      DO 60 ABREP=1,NREP
        JREP=MULTB(JABREP,ABREP+NREP,2)
        JOFF=IO(JREP)
        DO 45 BREP=1,NREP
         BOFF=IV(BREP) + IO(NREP+1)
         AREP=MULTB(BREP,ABREP+NREP,2)
!        IF(AREP.LT.BREP) GOTO45
         IF(AREP.LT.BREP) GOTO 45
         AOFF=IV(AREP) + IO(NREP+1)
c     x=multb(arep,brep,1)
c     write(iw,*) 'MAKEEKOVV irrep',multb(jrep,x+nrep,2)
         DO B=1,NV(BREP)
           EPSB=EPS(BOFF+B)
           AMIN=1
           IF(AREP.EQ.BREP) AMIN = B + 1
           DO A=AMIN,NV(AREP)
             EPSA=EPS(AOFF+A)
             DO J=1,NO(JREP)
               EPSJ=EPS(JOFF+J)
               OFF=OFF+1
               E(OFF)=-EPSJ+EPSA+EPSB
               ICRB(1,OFF)=J
               ICRB(2,OFF)=JREP
               ICRB(3,OFF)=A
               ICRB(4,OFF)=AREP
               ICRB(5,OFF)=B
               ICRB(6,OFF)=BREP
             ENDDO
           ENDDO
         ENDDO
! 45     ENDDO
 45     CONTINUE
!60   ENDDO
 60   CONTINUE

      IF( (OFF).NE.NOVVT(JABREP)) THEN
        WRITE(IW,*) ' *** Expected value of counter:',NOVVT(JABREP)
        WRITE(IW,*) ' *** Actual value of counter:  ',OFF
        CALL QUIT('Inconsistency in MAKEEKOVV')
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ASYMSPA(E,LDE,ASYM,SPARSE)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Determines deviation from symmetric form and sparseness for arbitrary 
C     quadratic matrix E with leading dimension LDE. It distinguishes between
C     E being real or complex.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER LDE
      REAL*8 E(LDE,LDE),ASYM,SPARSE
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER IZE,INZ
      REAL*8 RV
      COMPLEX*16 CV
C
C---------------Executable code--------------------------------------
C
      ASYM=0.0D0
      SPARSE=0.0D0
      IZE=0
      
      DO J=1,LDE
        DO I=J+1,LDE
          IF(CARITH) THEN      
            CV = E(I,J) - E(J,I)
#if defined (BIT64)
            ASYM = ASYM + DBLE(CONJG(CV)*CV)
#else
            ASYM = ASYM + DBLE(DCONJG(CV)*CV)
#endif
            IF(E(I,J).EQ.A0) IZE=IZE+1
          ELSE
            RV = E(I,J) - E(J,I)
            ASYM = ASYM + RV*RV
            IF(E(I,J).EQ.AR0) IZE=IZE+1
          ENDIF
        ENDDO
      ENDDO
C     write(iw,*) ' **** explicit diagonal elements ****'
C     DO J=1,LDE
C       write(iw,*) 'I,I,A:',J,J,E(J,J)
C     ENDDO
C
C  diagonal elements are nonzero and we calculate for the full matrix
C
      SPARSE=DBLE(2*IZE)/DBLE(LDE*LDE)*100.0D0
      INZ=LDE*(LDE-1)/2 - IZE
      WRITE(IW,*) 'Nonzero elements in lower triangle:',INZ
       
C       DO I=1,30
C         write(iw,'(I3,8E13.4)') I,(E(I,J),J=1,8)
C       ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION QMASYMR(E,LDE)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes deviation from symmetry of a quadratic matrix
C
      REAL*8 QMASYMR
      INTEGER LDE
      REAL*8 E(LDE,LDE),ASYM
C
C---------------Executable code--------------------------------------
C
      ASYM=0.0d0
      IF(LDE.eq.1) THEN
       QMASYMR = ASYM
       RETURN
      ENDIF

      DO I=2,LDE
        DO J=1,I-1
          ASYM = ASYM + DABS(E(I,J) - E(J,I))
        ENDDO
      ENDDO

      QMASYMR = ASYM

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION QMASYMQ(E,LDE)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes deviation from symmetry of a quadratic matrix
C
      REAL*8 QMASYMQ
      INTEGER LDE
      REAL*8 ASYM
      COMPLEX*16 E(LDE,LDE)
C
C---------------Executable code--------------------------------------
C
      ASYM=0.0d0
      IF(LDE.eq.1) THEN
       QMASYMQ = ASYM
       RETURN
      ENDIF

      DO I=2,LDE
        DO J=1,I-1
          ASYM = ASYM + ABS(E(I,J) - E(J,I))
        ENDDO
      ENDDO

      QMASYMQ = ASYM

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MATNEG(E,LDE)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Multiplies matrix E of dimension LDE x LDE by (-1)
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER LDE
      REAL*8 E(LDE,LDE)
C
C---------------Common Blocks--------------------------------------
C---------------Local variables--------------------------------------
C---------------Executable code--------------------------------------

      DO I=1,LDE
        DO J=1,LDE
          E(J,I) = - E(J,I)
        ENDDO
      ENDDO
 
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WRITE_HHP(CKKS,BUFHP,COL,BUF,
     &                   INTBUF,LADC,IOI,IOJ,
     &                   CREP,NBUFS,IB,ITAPE2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Writes out the H/H2P buffer to file. The remaining elements will
C     be written out in the next call to this routine.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 CKKS(*),BUFHP(*),COL(*),BUF(*)
      INTEGER INTBUF,LADC,IOI(*),IOJ(*)
      INTEGER CREP,NBUFS,IB,ITAPE2
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/symm.inc"
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER JDUMMY
C---------------Executable Code --------------------------------------
C
C
C  clear output buffer, index arrays and column buffer
C
      DO I=1,INTBUF
        BUF(I)=0.0D0
        IOI(I)=0
        IOJ(I)=0
      ENDDO
      CALL XCOPY(LADC,A0,0,COL,1)
      JDUMMY = 0

C
C  write out h/h and h/2h1p buffers
C
      NOK=NO(CREP)
      NAJL=NVOOT(CREP)

      OFF1 = 1
      OFF3 = 1

      DO 200 K=1,NOK
        ICOL=K
        OFF2=1
        CALL XCOPY(NOK,CKKS(OFF1),1,COL(OFF2),1)
        OFF2 = OFF2 + NOK
        CALL XCOPY(NAJL,BUFHP(OFF3),1,COL(OFF2),1)

C
C  one column is ready. we write out =>> its negative <<= and
C  only THE LOWER TRIANGLE !! 
C  the negative of the ADC matrix has the eigenvalues where we
C  want them
C
        DO IROW = ICOL+1,LADC
          IF(COL(IROW).NE.0.0D0) THEN
            IB=IB+1
            BUF(IB)=COL(IROW)
            IOI(IB)=IROW
            IOJ(IB)=ICOL
            IF(IB.EQ.INTBUF) THEN
              NBUFS = NBUFS + 1
              WRITE(ITAPE2) (-BUF(IXX),IXX=1,INTBUF),
     &                      (IOI(IXX),IXX=1,INTBUF),
     &                      (IOJ(IXX),IXX=1,INTBUF),
     &                      INTBUF,JDUMMY
              IB = 0
            ENDIF
          ENDIF
        ENDDO
C
C  we adjust offsets and go to next column
C
        OFF1 = OFF1 + NOK
        OFF3 = OFF3 + NAJL
 200  CONTINUE

C
C  1h/2h-1p block finished
C  Important: we leave with current IO-buffer position IB
C  and current number of buffers NBUFS for the next part.
C  ==> If there are unwritten elements in the buffer they
C  will be cleared in the final part of the ADC matrix
C
      RETURN 
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CWRITE_HHP(CKKS,BUFHP,COL,BUF,
     &                   INTBUF,LADC,IOI,IOJ,
     &                   CREP,NBUFS,IB,ITAPE2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Writes out the H/H2P buffer to file in the complex case.
C     The remaining elements will
C     be written out in the next call to this routine.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 CKKS(*),BUFHP(*),COL(*),BUF(*)
      INTEGER INTBUF,LADC,IOI(*),IOJ(*)
      INTEGER CREP,NBUFS,IB,ITAPE2
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/symm.inc"
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER JDUMMY
C---------------Executable Code --------------------------------------
C
C
C  clear output buffer, index arrays and column buffer
C
      DO I=1,INTBUF
        BUF(I)=(0.0D0,0.0D0)
        IOI(I)=0
        IOJ(I)=0
      ENDDO
      CALL ZCOPY(LADC,A0,0,COL,1)
      JDUMMY = 0

C
C  write out h/h and h/2h1p buffers
C
      NOK=NO(CREP)
      NAJL=NVOOT(CREP)

      OFF1 = 1
      OFF3 = 1

      DO 200 K=1,NOK
        ICOL=K
        OFF2=1
        CALL ZCOPY(NOK,CKKS(OFF1),1,COL(OFF2),1)
        OFF2 = OFF2 + NOK
        CALL ZCOPY(NAJL,BUFHP(OFF3),1,COL(OFF2),1)
C
C  one column is ready. we write out =>> its negative <<= and
C  only THE LOWER TRIANGLE !!
C  the negative of the ADC matrix has the eigenvalues where we
C  want them
C
        DO IROW = ICOL+1,LADC
          IF(COL(IROW).NE.(0.0D0,0.0D0)) THEN
            IB=IB+1
            BUF(IB)=COL(IROW)
            IOI(IB)=IROW
            IOJ(IB)=ICOL
            IF(IB.EQ.INTBUF) THEN
              NBUFS = NBUFS + 1
              WRITE(ITAPE2) (-BUF(IXX),IXX=1,INTBUF),
     &                      (IOI(IXX),IXX=1,INTBUF),
     &                      (IOJ(IXX),IXX=1,INTBUF),
     &                      INTBUF,JDUMMY
              IB = 0
            ENDIF
          ENDIF
        ENDDO
C
C  we adjust offsets and go to next column
C
        OFF1 = OFF1 + NOK
        OFF3 = OFF3 + NAJL
 200  CONTINUE

C
C  1h/2h-1p block finished
C  Important: we leave with current IO-buffer position IB
C  and current number of buffers NBUFS for the next part.
C  ==> If there are unwritten elements in the buffer they
C  will be cleared in the final part of the ADC matrix
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RDVOVV (VOVV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     This routine reads in all VOVV integrals. accounts for offsets
C     in the complex case.
C
C---------------Routines called----------------------------------------
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VOVV(*)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/symm.inc"
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
      integer*8 OFF1
C
C---------------Executable code--------------------------------------
C
      OFF1 = 1
      DO 30 IRP = 1, NREP

         ISTART = 0
         IF (NVVT(IRP).EQ.0) GOTO 30
         MINT = NVVT(IRP)
         CALL GETVOVV (IRP,ISTART,NINT,DONE,VOVV(OFF1),MINT)

         IF (.NOT.DONE) THEN
            WRITE(IW,*) 'Irrep',IRP,' unfinished'
            CALL QUIT('VOVV reading error in RDVOVV')
         ENDIF
         OFF1 = OFF1 + NVO(IRP)*NVVT(IRP)*RCW
 30   CONTINUE

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RDVVVV (VVVV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     This routine reads in all VVVV integrals (we would not do that if
C     they weren't needed completely in memory anyway due to the random
C     access in the UIONPA routine). Also accounts for complex
C     integrals.
C
C---------------Routines called----------------------------------------
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VVVV(*)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/symm.inc"
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
      integer*8 off1
C
C---------------Executable code--------------------------------------
C
      OFF1 = 1
      DO 30 IRP = 1, NREP
         ISTART = 0
         IF (NVVT(IRP).EQ.0) GOTO 30
         MINT = NVVT(IRP)
         CALL GETVVVV (IRP,ISTART,NINT,DONE,VVVV(OFF1),MINT)
         IF (.NOT.DONE) THEN
            WRITE(IW,*) 'Irrep',IRP,' unfinished'
            CALL QUIT('VVVV reading error in RDVVVV')
         ENDIF
         OFF1 = OFF1 + NVVT(IRP)*NVVT(IRP)*RCW
 30   CONTINUE

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RWQKL(QOO,QVO,QVV,LRW,JSTAT,IOCH,IONAME)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Since the construction of the Q_kl takes most of the time due to
C     the necessary inclusion of all symmetries we provide the
C     possibility to store and reread the qkl entities dependent on the
C     input variable READQKL = T/F.
C
C     The Qkl arrays always come as real*8. they have twice the length
C     if cokmplex arithmetic is turned on. This will be accounted for.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 QOO(*),QVO(*),QVV(*)
      LOGICAL LRW
      INTEGER JSTAT
      CHARACTER*6 IONAME
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/symm.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL LT
      INTEGER SLEN
C
C---------------Executable code--------------------------------------
C
C  first we treat the reading of Q_kl
C
      JSTAT = 0
      IF(.NOT.LRW) THEN
        INQUIRE (FILE=IONAME,EXIST=LT)
C
C  QKL file is not there, we have to generate it
C  code (JSTAT = 1). this is not an error code.
C  otherwise we read in the stored Qkl values
C
        IF(.NOT.LT) THEN
          JSTAT=1
          RETURN
        ENDIF
C
C  files are all there. we read in the Qkl values
C
        OPEN(IOCH,FILE=IONAME,FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IOCH)
        SLEN=MOO(1)*RCW
        READ(IOCH) (QOO(IX),IX=1,SLEN)
        SLEN=MVO(1)*RCW
        READ(IOCH) (QVO(IX),IX=1,SLEN)
        SLEN=MVV(1)*RCW
        READ(IOCH) (QVV(IX),IX=1,SLEN)
        CLOSE(IOCH)
        
      ELSE
C
C  here we treat the writing of Q_kl
C
        OPEN(IOCH,FILE=IONAME,FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IOCH)
        SLEN=MOO(1)*RCW
        WRITE(IOCH) (QOO(IX),IX=1,SLEN)
        SLEN=MVO(1)*RCW
        WRITE(IOCH) (QVO(IX),IX=1,SLEN)
        SLEN=MVV(1)*RCW
        WRITE(IOCH) (QVV(IX),IX=1,SLEN)
        CLOSE(IOCH)
        
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MATSYM(E,NR,NC)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Symmetrizes the REAL quadratic matrix E by adding its transpose
C     to itself. Extension to complex: Hermitian conjugate 
C     Attention: diagonal is treated separately otherwise we end up with
C     three times the original value !
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER NR,NC
      REAL*8 E(NR,NC)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 H
C
C---------------Executable code--------------------------------------
C
      IF(NR.NE.NC) CALL 
     & QUIT('Rectangular matrix given in MATSYM')
      DO I=1,NR
        E(I,I) = E(I,I) + E(I,I)
        DO J=I+1,NR
          H = E(I,J)
          E(I,J) = E(I,J) + E(J,I)
          E(J,I) = E(J,I) + H
        ENDDO
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMATSYM(E,NR,NC)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Adds Hermitian conjugate of E to complex matrix E
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER NR,NC
      COMPLEX*16 E(NR,NC)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 H
C
C---------------Executable code--------------------------------------
C
      IF(NR.NE.NC) CALL
     & QUIT('Rectangular matrix given in CMATSYM')
      DO I=1,NR
        E(I,I) = DCONJG(E(I,I)) + E(I,I)
        DO J=I+1,NR
          H = E(I,J)
          E(I,J) = E(I,J) + DCONJG(E(J,I))
          E(J,I) = E(J,I) + DCONJG(H)
        ENDDO
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MATOUT(E,NR,NC)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Prints out a real quadratic matrix E of dimension N x N 
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER NR,NC
      REAL*8 E(NR,NC)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C
C---------------Executable code--------------------------------------
C
      WRITE(IW,*) 'matrix dimension:',NR,' x ',NC
      WRITE(IW,*)
      DO I=1,NR
        WRITE(IW,'(I6,100E15.7)') I,(E(I,IXX),IXX=1,NC)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MATOUTEV(E,NR,NC)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Prints out a real quadratic matrix E of dimension N x N in eV
C     coming in with entries in atomic units
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER NR,NC
      REAL*8 E(NR,NC)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 AUTOEV
      PARAMETER(AUTOEV = 27.2113957D0)
C
C---------------Executable code--------------------------------------
C
      WRITE(IW,*) 'matrix dimension:',NR,' x ',NC
      WRITE(IW,*)
      DO I=1,NR
        WRITE(IW,'(I6,100E15.7)') I,(E(I,IXX)*AUTOEV,IXX=1,NC)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMATOUT(E,NR,NC)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Prints out a complex quadratic matrix E of dimension NR x NC
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER NR,NC
      COMPLEX*16 E(NR,NC)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      WRITE(IW,*) 'matrix dimension:',NR,' x ',NC
      WRITE(IW,*)
      DO I=1,NR
        WRITE(IW,'(I6,100E15.7)') I,(E(I,IXX),IXX=1,NC)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMATOUTEV(E,NR,NC)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Prints out a hermitian quadratic matrix E of dimension NR x NC in eV
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER NR,NC
      COMPLEX*16 E(NR,NC)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 AUTOEV
      PARAMETER(AUTOEV = (27.2113957D0,0.0D0))
C
C---------------Executable code--------------------------------------
C
      WRITE(IW,*) 'matrix dimension:',NR,' x ',NC
      WRITE(IW,*)
      DO I=1,NR
        WRITE(IW,'(I6,200E15.7)') I,(E(I,IXX)*AUTOEV,IXX=1,NC)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Subroutine Locate(keyword,ilength,ifound)
      Character*10 keyword,line*80
      integer ilength,ifound
C
      ifound = 1
      rewind(5)
 10   read(5,'(A80)',END=15,ERR=15) line
      if (line(3:(ilength+2)).eq.keyword(1:ilength)) then
         backspace(5)
         goto 20
      endif
      goto 10
 15   ifound = -1
 20   Return
      End
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     
      SUBROUTINE PST(A)
c
c  this routine headlines an arbitrary character string with a + delimiter
c  in the input. is is meant to inform the user about program progress.
c
      IMPLICIT INTEGER(A-Z)
      CHARACTER*60 A

#include "../relccsd/files.inc"
      character*60 b,bt

      lim=60
      m=lim
      do i=1,lim
        if(a(i:i).eq.'+') then
           m=i
           goto 44
        endif
      enddo
 44   do i=1,m-1
        bt(i:i) = a(i:i)
        b(i:i) = '#'
      enddo
      do i=m,lim
        bt(i:i)=' '
        b(i:i)=' '
      enddo

      write(iw,*)
      write(iw,*) '####',b
      write(iw,*) '###'
      write(iw,*) '### ',bt
      write(iw,*) '###'
      write(iw,*) '####',b
      write(iw,*)

      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
