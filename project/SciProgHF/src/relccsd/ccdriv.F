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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCHEADER
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Write header for Relativistic Coupled Cluster Program
C     Keep track of all the modifications made in the program.
C     Modification record is now obsolete, consult the CVS history !
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER DATEX*10,TIMEX*8
C
C---------------Executable code--------------------------------------
C
      WRITE (IW,1000)
      WRITE (IW,1002)
      CALL DAYTIME (DATEX,TIMEX)
      WRITE (IW,1001) DATEX,TIMEX
C
 1000 FORMAT (//" Relativistic Coupled Cluster program RELCCSD"
     &//" Written by :",/"    Lucas Visscher"/
     &/"    NASA Ames Research Center    (1994)"
     &/"    Rijks Universiteit Groningen (1995)"
     &/"    Odense Universitet           (1996-1997)"
     &/"    VU University Amsterdam      (1998-present)")
 1001 FORMAT (//" Today is :",T15,A10/" The time is :",T17,A8)
 1002 FORMAT (//" This module is documented in"
     &/"  - Initial implementation :",T30,
     &"L. Visscher, T.J. Lee and K.G. Dyall, ",
     &"J. Chem. Phys. 105 (1996) 8769."
     &/"  - Fock Space (FSCC)      :",T30,
     &"L. Visscher, E. Eliav and U. Kaldor, ",
     &"J. Chem. Phys. 115 (2002) 9720."
     &/"  - Parallelization        :",T30,
     &"M. Pernpointner and L. Visscher, ",
     &"J. Comp. Chem. 24 (2003) 754.",
     &/"  - Intermediate Hamilt. FS:",T30,
     &"E. Eliav, M. J. Vilkas, Y. Ishikawa, and U. Kaldor, ",
     &"J. Chem. Phys. 122 (2005) 224113.",
     &/"  - MP2 expectation values :",T30,
     &"J.N.P. van Stralen, L. Visscher, C.V. Larsen and H.J.Aa. Jensen,",
     &" Chem. Phys. 311 (2005) 81."
     &/"  - CC  expectation values :",T30,
     &"A. Shee, L. Visscher, and T. Saue, ",
     &"J. Chem. Phys. 145 (2016) 184107."
     &/"  - EOM-IP/EA/EE energies  :",T30,
     &"A. Shee, T. Saue, L. Visscher, and A.S.P. Gomes, ",
     &"J. Chem. Phys. 149 (2018) 174113.",
     &/"  - Core spectra (CVS-EOM) :",T30,
     &"L. Halbert, M. L. Vidal, A. Shee, S. Coriani, and A.S.P. Gomes,",
     &" arXiv:2011.08549")
C
      RETURN
C---------------Record of modifications made in the code-------------
C     modifications after 2000 are documented in the cvs system
C--------------------------------------------------------------------
C     July 1999: Implemented Fock space Coupled Cluster
C                Parallelized (part of) the program using MPI
C--------------------------------------------------------------------
C     11-12-97 : MP2 gradients works, CCSD Lambda equations written
C--------------------------------------------------------------------
C     October 1997 : major restructuring
C     Split the code in energy, gradient and second derivative modules
C     Final version should contain energies and 1st and 2nd order
C     properties at all the available levels of theory
C     Improved memory management (thanks to H-J Aa Jensens routines)
C--------------------------------------------------------------------
C     Aug  6 1997 :
C     Implemented excitation energy module.
C--------------------------------------------------------------------
C     Apr 16 1997 :
C     Interfaced with DIRAC, reorganized code.
C--------------------------------------------------------------------
C     Jul 10 1996 :
C     Added RPA module, changed output format,
C     added NOCCSD (MP2) option
C     Note that the RPA is not very optimized : needs complex
C     arithmetics for completely imaginary operators, does not
C     consider symmetry of the operators
C--------------------------------------------------------------------
C     Oct 5 1995 :
C     Corrected errors in complex arithmetic, implemented -T correction
C     Frozen option was also implemented somewhat earlier
C--------------------------------------------------------------------
C     Aug 13 1995 :
C     Implemented triples correction in symmetry-adapted version
C--------------------------------------------------------------------
C     Jun 29 1995 :
C     Dynamic memory allocation & incore files for CRAY
C--------------------------------------------------------------------
C     Jan 10 1995 :
C     Real/complex arithmetic possible
C     Symmetry packed version
C--------------------------------------------------------------------
C     Oct 19 1994 :
C     Writes integrals out on disk
C     Otherwise identical to previous version
C--------------------------------------------------------------------
C     Sep 22 1994 : First numbered version of program
C     Keeps all integrals in core.
C     Does not use any symmetry.
C     But is vectorized by using matrix multiplies where ever possible
C--------------------------------------------------------------------
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCINPT (IN,IW)
      use relcc_cfg
!     All initialization should be done in relcc_cfg, ultimately all common
!     blocks should go, keep them for the moment for compatibility with reladc
C
      IMPLICIT NONE
      INTEGER IN,IW
C
C---------------Description--------------------------------------------
C
C     Reads input for relativistic Coupled Cluster program
C
C   Routine is called twice - first from PSIINP/dirrd.F as part of input reading
C and second time from own RelCC module. To minimize the amount of output we
C need to distinguish between calls.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C              Miro Ilias, July/2015 - add control printouts
C
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
#include "inpt.inc"
#include "ihm.inc"
#include "complex.inc"
#include "symm.inc"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbmp2no.h"
#include "dcbham.h"
C
C---------------Local variables--------------------------------------
C
      LOGICAL REWINP,DOEA,DOIE,DOEA2,DOIE2,DOEXC
      REAL*8 EHMIN,EHMAX,EPMIN,EPMAX
      REAL*8 SHIFTH11,SHIFTH2,SHIFTP11,SHIFTP2,SHIFTH12,SHIFTHP,SHIFTP12
C
!     This namelist is the only one kept for backward compatibility (old HF Hessian).
!     This is to be replaced when a real CC or MP2 Hessian becomes available.
      NAMELIST/CCSOPR/ NAOPER,NBOPER,NFREQ,NAMEA,NAMEB,EFREQ,
     &                 MAXIT,MAXDIM,NTOL
      DOUBLE PRECISION D0
      PARAMETER (D0 = 0.0D00)
      INTEGER NACTH,NACTP,NACTHI,NACTPI,NTOL
      DIMENSION NACTH(16),NACTP(16)
      DIMENSION NACTHI(16),NACTPI(16)
      INTEGER MAXIT00,MAXIT01,MAXIT10,MAXIT11,MAXIT20,MAXIT02,I,IFREQ,
     &        ISECT,IPRP_WF,IOPER,MAXDIM,MAXIT,TSHOLD

      Logical, SAVE :: FirstCall=.TRUE.

C
C---------------Executable code--------------------------------------

C
C     Convert input from module file to old-fashioned keywords
C
      DO I= 1, 16
       if (nelec_input)  then
!         Overwrite information provided by SCF program with user-specified input
          NELEC(I)    = relcc_nelec(I)
       else if (nelec_open_input) then
!         Open shells require user input because the SCF-initialized nelec counts only closed shell electrons
          NELEC(I) = NELEC(I) + relcc_nelec_open(I) 
       end if 
       NFROZ(I)    = relcc_nfroz(I)
       NELEC_F1(I) = relcc_nelec_f1(I)
       NELEC_F2(I) = relcc_nelec_f2(I)
      ENDDO
      IPRNT   = relcc_print
      DOSORT  = relcc_do_sort
      DOENER  = relcc_do_energy
      DOFOPR  = relcc_do_gradient
      DOSOPR  = relcc_do_hessian
      DOFSPC  = relcc_do_fspc
      DOIH    = relcc_fs_do_ih
      DORESTART = relcc_do_restart
      DOMEMCOUNT = relcc_do_count_memory
      DEBUG   = relcc_debug
      TIMING  = relcc_timing
      MWORD   = relcc_memory_in_mw
      INTERFACE = relcc_integral_interface

!     Print section
      CALL PRSYMB(IW,'=',75,0)
      WRITE(IW,'(A)')
     & ' **RELCC: Set-up for Coupled Cluster calculations'
      CALL PRSYMB(IW,'=',75,0)

      If (.NOT.FirstCall)  then
        WRITE(IW,'(1X,A,I3)') '* General print level   : ',IPRNT
        IF (IPRNT.GE.1) THEN
          WRITE(IW,"(1X,A,16I3)") "NEL_F1:",(NELEC_F1(I),I=1,16)
        ENDIF
      EndIf

!MI  deal with memory count demand
      if (DOMEMCOUNT.AND.FirstCall)  then
        print *,
     &'Asked for RelCC memory consumption count. Afterwards leave.'
      endif
C
C     Input for the energy calculation.
C
      MAXDIM  = relcc_ccener_max_dimension_diis
      MAXIT   = relcc_ccener_max_iterations
      NTOL    = relcc_ccener_desired_convergence
      DOMP2   = relcc_do_mp2
      DOCCSD  = relcc_do_ccsd
      DOCCSDT = relcc_do_ccsd_t
      NOCCS   = relcc_no_singles
      NOCCD   = relcc_no_doubles
      MXDIMCC = MAXDIM
      MXITCC = MAXIT
      NTOLCC = NTOL
C
C     Input for the gradient calculation.
C
      DOMP2G    = relcc_do_mp2gradient
      DOMP2GOLD = relcc_do_oldmp2gradient
      DOCCSDG   = relcc_do_ccsdgradient
      DOCCSDTG  = relcc_do_ccsdtgradient
      NOZG      = .not.relcc_do_relaxed
      DONATORB  = relcc_do_naturalorbitals
      NEOPER    = 3
      NAMEE(1)  = 'XDIPLEN '
      NAMEE(2)  = 'YDIPLEN '
      NAMEE(3)  = 'ZDIPLEN '
      MAXDIM    = relcc_fopr_max_dimension_diis
      MAXIT     = relcc_fopr_max_iterations
      NTOL      = relcc_fopr_desired_convergence
C     Switch on evaluation of property output if necessary
      IF (DOMP2G) THEN
          DOEXP = .TRUE.
          IPRP_WF = 0
          DO I = 1, NPRP_WF
             IF (PRP_WF(I).EQ.'MP2 ') IPRP_WF = I
          END DO
          IF (IPRP_WF.EQ.0) THEN
             IF (NPRP_WF .GE. MXPRP_WF) THEN
                IPRP_WF = 1 ! Unlikely that this happens, just skip the HF print out in this case
             ELSE
                NPRP_WF = NPRP_WF + 1
                IPRP_WF = NPRP_WF
             END IF
             PRP_WF(IPRP_WF) = 'MP2 '
          END IF
      END IF

!     following section is for switching on coupled cluster property evaluation.
      IF (DOCCSDG) THEN
          DOEXP = .TRUE.
          IPRP_WF = 0
          DO I = 1, NPRP_WF
             IF (PRP_WF(I).EQ.'CCSD') IPRP_WF = I
          END DO
          IF (IPRP_WF.EQ.0) THEN
             IF (NPRP_WF .GE. MXPRP_WF) THEN
                IPRP_WF = 1 ! Unlikely that this happens, just skip the HF print out in this case
             ELSE
                NPRP_WF = NPRP_WF + 1
                IPRP_WF = NPRP_WF
             END IF
             PRP_WF(IPRP_WF) = 'CCSD'
          END IF
      END IF
! the following section is to activate excited state property evaluation module

      if (relcc_do_eomprop) then
        DOEXP = .TRUE.
          IPRP_WF = 0
          DO I = 1, NPRP_WF
             IF (PRP_WF(I).EQ.'EOM ') IPRP_WF = I
          END DO
          IF (IPRP_WF.EQ.0) THEN
             IF (NPRP_WF .GE. MXPRP_WF) THEN
                IPRP_WF = 1 ! Unlikely that this happens, just skip the HF print out in this case
             ELSE
                NPRP_WF = NPRP_WF + 1
                IPRP_WF = NPRP_WF
             END IF
             PRP_WF(IPRP_WF) = 'EOM '
          END IF
      END IF

      MXDIMGR = MAXDIM
      MXITGR = MAXIT
      NTOLGR = NTOL
C
C     Input for the linear response calculation.
C
      NAOPER = 1
      NBOPER = 1
      NFREQ  = 1
      DO I = 1, 16
         NEXC(I) = 1
      ENDDO
      NAMEA(1) = 'ZDIPLEN '
      NAMEB(1) = 'ZDIPLEN '
      DO IFREQ = 1, MAXFR
         EFREQ(IFREQ) = D0
      ENDDO
      MAXDIM =  relcc_max_dimension_diis
      MAXIT  =  relcc_max_iterations
      NTOL   =  relcc_desired_convergence
      REWIND (IN)
      READ (IN,CCSOPR,END=8,ERR=7)
      GOTO 8
    7 CALL QUIT('Error exit while reading namelist CCSOPR')
    8 CONTINUE
      MXDIMRP = MAXDIM
      MXITRP  = MAXIT
      NTOLRP  = NTOL
C
C     Input for the excitation energy calculation. (no longer supported, corresponded to TD-HF in MO basis)
C
      DO I = 1, 16
         NEXC(I) = 1
      ENDDO
      MAXDIM =  relcc_max_dimension_diis
      MAXIT  =  relcc_max_iterations
      NTOL   =  relcc_desired_convergence
      MXDIMEX = MAXDIM
      MXITEX = MAXIT
      NTOLEX = NTOL
C
C     Input for the Fock-Space CCSD energy calculation.
C
      MAXDIM =  relcc_fs_max_dimension_diis
      MAXIT  =  relcc_fs_max_iterations
      NTOL   =  relcc_fs_desired_convergence
      TSHOLD =  relcc_fs_tshold
      DO I = 1, 16
          NACTH(I)     = relcc_fs_nacth(I)
          NACTP(I)     = relcc_fs_nactp(I)
      ENDDO
      DOEA      = relcc_fs_do_ea
      DOIE      = relcc_fs_do_ie
      DOEXC     = relcc_fs_do_exc
      DOEA2     = relcc_fs_do_ea2
      DOIE2     = relcc_fs_do_ie2
      FSSECT(1) = relcc_fs_fssect(1)
      DO ISECT = 2, 6
          FSSECT(ISECT) = relcc_fs_fssect(ISECT)
      ENDDO

      if ((DOIE .eqv. .true.).and.( ALL(NACTH == 0) ) ) then
        write(IW,*) 'WARNING: you did not specify the set of active 
     &hole orbitals from which electron ionization takes place .NACTH!'
      endif

      if ((DOEA .eqv. .true.).and.( ALL(NACTP == 0) ) ) then
        write(IW,*) 'WARNING: you did not specify the set of active 
     &particle orbitals to which electron attachment takes place .NACTP'
      endif
      
! in case we want to jump to a given sector
      do isect = 1, 6
         FS_SKIP_SECT(isect) = relcc_restart_skipsect(isect)
         FS_REDO_SECT(isect) = relcc_restart_redosect(isect)
         if ((FS_SKIP_SECT(isect) .eqv. .true.).and.
     &       (FS_REDO_SECT(isect) .eqv. .true.)) then
            call quit('Inconsistent input for skip/redo_sect')
         endif
      enddo
      restart_unconverged  = relcc_restart_unconverged
      restart_redo_sort    = relcc_restart_redo_sorting
      restart_redo_ccsd    = relcc_restart_redo_ccsd
      restart_ignore_check = relcc_restart_ignore_check

CMI   ...provide with default values
      MAXIT02=relcc_fs_max02_iterations
      MAXIT01=relcc_fs_max01_iterations
      MAXIT00=relcc_fs_max00_iterations
      MAXIT10=relcc_fs_max10_iterations
      MAXIT11=relcc_fs_max11_iterations
      MAXIT20=relcc_fs_max20_iterations
!     stefan: select state sorted by energy for numerical gradient calculation
      STATE_E_FSCC = relcc_fs_select_state_for_numgrad_energy
      MXDIMFS = MAXDIM
      MXITFS = MAXIT
      DO ISECT = 1, 6
CMI     ...default maximum number iterations per sector
        MAXIT_SECT(ISECT) = MAXIT
      ENDDO
      NTOLFS = NTOL
      LTSHFS = TSHOLD
      DO I = 1, 16
        NINACT(I) = NELEC(I) - NACTH(I)
        NACT(I) = NACTH(I) + NACTP(I)
        IF (.not.FirstCall.and.NINACT(I).LT.0) THEN
           WRITE (IW,*) "more holes than electrons in input"
           CALL QUIT ("Input error")
        END IF
      ENDDO
CMI   ... translate MAXITxy easy input to MAXIT_SECT type input ...
      IF (MAXIT02.NE.-1) MAXIT_SECT(5) = MAXIT02
      IF (MAXIT01.NE.-1) MAXIT_SECT(2) = MAXIT01
      IF (MAXIT10.NE.-1) MAXIT_SECT(3) = MAXIT10
      IF (MAXIT00.NE.-1) MAXIT_SECT(1) = MAXIT00
      IF (MAXIT11.NE.-1) MAXIT_SECT(4) = MAXIT11
      IF (MAXIT20.NE.-1) MAXIT_SECT(6) = MAXIT20
C
C     Translate easy input to FSSECT type input
C
      IF (DOEA) THEN
         FSSECT(1) = 1
         FSSECT(2) = 1
      ENDIF
C
      IF (DOIE) THEN
         FSSECT(1) = 1
         FSSECT(3) = 1
      ENDIF
C
      IF (DOEXC) THEN
         FSSECT(1) = 1
         FSSECT(2) = 1
         FSSECT(3) = 1
         FSSECT(4) = 1
      ENDIF
C
      IF (DOEA2) THEN
         FSSECT(1) = 1
         FSSECT(2) = 1
         FSSECT(5) = 1
      ENDIF
C
      IF (DOIE2) THEN
         FSSECT(1) = 1
         FSSECT(3) = 1
         FSSECT(6) = 1
      ENDIF
C
C     General input for IH procedure
C
      EHMIN    = relcc_ih_eh_min
      EHMAX    = relcc_ih_eh_max
      EPMIN    = relcc_ih_ep_min
      EPMAX    = relcc_ih_ep_max
      DO I = 1, 16
         NACTH(I)     = relcc_fs_nacth(I)! it's 0 for (0,2) sector
         NACTP(I)     = relcc_fs_nactp(I)
      ENDDO
      IHSCHEME = relcc_ih_scheme
C
C     Defaults for scheme 1 that shifts the denominators
C
      SHIFTH11 = relcc_ih_shift_h11
      SHIFTH12 = relcc_ih_shift_h12
      SHIFTP11 = relcc_ih_shift_p11
      SHIFTP12 = relcc_ih_shift_p12
      SHIFTH2  = relcc_ih_shift_h2
      SHIFTP2  = relcc_ih_shift_p2
      SHIFTHP  = relcc_ih_shift_hp
      AIH      = relcc_ih_aih
      NIH      = relcc_ih_nih
C
C     Some processing of the input is necessary
C
      ER_IH(1,1) = EHMIN ! for (1,0) or (2,0) sectors plus (1,1) sector
      ER_IH(2,1) = EHMAX ! the same
      ER_IH(1,2) = EPMIN ! ... limits for Pi space for (0,1)/(0,2) plus (1,1) sectors
      ER_IH(2,2) = EPMAX

C
C     XIH energy shift sign definition.
C     Array with as first dimension 1e/2e shift
C                   second dimension sector number
C
C     Shifts for particles should be negative
      SHIFT_IH(1,2) = -DABS(SHIFTP11)
      SHIFT_IH(2,2) = -DABS(SHIFTP12)
      SHIFT_IH(2,5) = -DABS(SHIFTP2)
CEE     Shifts for holes must be negative (not positive)
!----------------------------------------------------------
      SHIFT_IH(1,3) =  -DABS(SHIFTH11)
      SHIFT_IH(2,3) =  -DABS(SHIFTH12)
      SHIFT_IH(2,6) =  -DABS(SHIFTH2)
C     Shifts for hole-particle excitations may be either negative or positive
      SHIFT_IH(2,4) =  SHIFTHP
C
C     Input for the sorting routine.
C     We can add a finite field here, if desired
C     This assumes that PRTRA is called as well !
C
      USEOE   = relcc_use_orbital_energies
      NORECMP = relcc_no_recompute
      ADD_FINITE_FIELD = relcc_ADD_FINITE_FIELD
      NFFOPER = relcc_NFFOPER
      DO IOPER = 1, MAXOP
         FF_PROP_NAMES(IOPER) = relcc_FF_PROP_NAMES(IOPER)
         FF_PROP_STRENGTHS(1,IOPER) = relcc_FF_PROP_STRENGTHS(1,IOPER)
         FF_PROP_STRENGTHS(2,IOPER) = relcc_FF_PROP_STRENGTHS(2,IOPER)
      ENDDO
!
!     X2Cmmf hamiltonian in use - by default always set USEOE and NORECMP to .true.
!     everything else results in nonsense energies - if overwritten by user input for example...
      if(X2CMMF)then
        if(.not.relcc_use_orbital_energies)
     &          relcc_use_orbital_energies = .true.
        if(.not.relcc_no_recompute)
     &          relcc_no_recompute = .true.
        USEOE   = .true.
        NORECMP = .true.
        write(IW,'(/A,/A,/A,/A,/A/)')
     &  ' ************************************************************',
     &  ' *** molecular-mean-field X2C Hamiltonan (X2Cmmf) active  ***',
     &  ' *** use of SCF orbital energies and no re-computation of ***',
     &  ' *** Fock matrix elements set explicitly                  ***',
     &  ' ************************************************************'
      end if
!
C
C     Some consistency and other checks
C
      IF (DOFSPC.AND..NOT.(DOENER.OR.DOFOPR.OR.DOSOPR)) DOSORT=.FALSE.

      if (DOENER.AND.DOFSPC) then
         relcc_do_energy = .false.
         DOENER = .FALSE.
      endif

      IF (DOFSPC) THEN
!     mi: print out range of intermediate Hamilt.spaces ... later ...
        write(IW,'(/,2X,A,F12.5,1X,A,F12.5)')
     &    'EHMIN:',ER_IH(1,1),'EHMAX:',ER_IH(2,1)
        write(IW,'(2X,A,F12.5,1X,A,F12.5)')
     &     'EPMIN:',ER_IH(1,2),'EPMAX:',ER_IH(2,2)
      ENDIF
C
C     IF( MP2_DENRUN )THEN
C       IN = IN_SAVE
C     END IF
C

!MI   After the first call of this routine reset the control variable
      FirstCall = .FALSE.

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCMAIN (IUNIN,IUNUT,CA,INTFCE)
c
      use memory_allocator
      use allocator_parameters, only : klongint, kreal
      use interface_to_mpi
      use relcc_cfg
      use xmlout
      use eom_driver
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Driver for Relativistic Coupled Cluster program
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "diis.inc"
#include "inpt.inc"
#include "symm.inc"
#include "eqns.inc"
#include "complex.inc"
#include "ccpar.inc"
#include "waio.h"
#include "dcbmp2no.h"
C
C---------------Local variables--------------------------------------
C
C     CC is the static memory  (when called as subroutine)
C     CC is the dynamic memory (when used as stand-alone code)
C
C
      INTEGER MAXSP
      PARAMETER (MAXSP=1000)
      INTEGER IREPSP(MAXSP*4),IREPSPI(MAXSP*MXREP*2)
      REAL*8 EPS(MAXSP),EPS1(MAXSP)
      LOGICAL ALLOC_ME,CA
      REAL*8 CPU0,CPU1,CPU2,CPUTOT,CPUI,CPUS
      REAL*8 CPUE,CPUFO,CPUSO,CPUFS,CPUVAR,CPUEOM
      REAL*8 WALL0,WALL1,WALL2,WALLTOT,WALLI,WALLS,TIME_SYNC_TOT
      REAL*8 WALLE,WALLFO,WALLSO,WALLFS,WALLVAR,WALLEOM,TIME_SYNC
      CHARACTER*10 INTFCE
      integer :: verbose    = 0
      integer :: ncore_elec = 0
      integer(kind=klongint) :: maxcor
      integer(kind=klongint) :: mxcorr
      integer(kind=klongint) :: mxcorr1
      integer(kind=klongint) :: mword_limit
      real(kind=kreal)       :: new_mword_limit = 0
      integer :: word,mycharge,NELECT,NSP
      integer :: jcode,IUNIN,IUNUT,I,ICALC,IERR
C
C---------------Executable code--------------------------------------
C
C     Set the input and output units
C
#if defined (CC_ALONE)
      MP2_DENRUN = .FALSE.
#endif
      INTERFACE = relcc_integral_interface
      CARITH = CA
      IN = IUNIN
      IW = IUNUT
      WRITE (IW,1011) MYPROC+1,NMPROC

C
      CALL CCHEADER
      CALL WAIO_INIT (IW)
C
C     Initialize the timing.
C
      CALL XTIME(0,-1,'                             ')
      CALL GETTIM(CPU0,WALL0)
      CPU1 = CPU0
      CPUI = 0
      CPUS = 0
      CPUE = 0
      CPUFO = 0
      CPUFS = 0
      CPUEOM = 0
      CPUSO = 0
      CPUTOT = 0
      WALL1 = WALL0
      WALLI = 0
      WALLS = 0
      WALLE = 0
      WALLFO = 0
      WALLFS = 0
      WALLEOM = 0
      WALLSO = 0
      WALLTOT = 0
C
C     Memory and I/O setup section
C     ----------------------------
C
      call allocator_get_words_available(maxcor)
      mword_limit = maxcor / (1024 * 1024)
      MXCORR = 0
C
C     Read symmetry information
C
      CALL RDSYMI (NSP,EPS1,IREPSP,ncore_elec)
C
      IF (NSP.GT.MAXSP) THEN
         PRINT*," INCREASE MAXSP IN CCMAIN TO",NSP
         CALL QUIT("NSP.GT.MAXSP")
      ENDIF
C
C     Read input from file (old namelist) / transfer input (new style) from relcc_cfg

!     sync co-workers and master with variables from new input...
#ifdef VAR_MPI
      call relcc_sync_cw()
#endif
C
      CALL CCINPT (IN,IW)
C
C     Number of integer words for REAL*8 and
C     number of REAL*8 words per variable (2 for complex groups,
C     1 for real)
C
#if  defined (INT_STAR8)
      IRW = 1
#else
      IRW = 2
#endif
      IF (CARITH) THEN
         RCW = 2
      ELSE
         RCW = 1
      ENDIF
C
C
C  Check the user supplied memory limit with the static maximum.
C
      IERR = 0
C
C     In case of dynamic allocation :
C     Allow for overhead : subtract 800 kWord
C
#if defined (CC_ALONE)
         MAXCOR = MAXCOR - 800 * 1024
#endif
C
C     Report memory in megaBytes
C
      if ( maxcor .gt. 0 ) then
        WRITE (IW,1000) DFLOAT(MAXCOR)/DFLOAT((128*1024)), 
     &  DFLOAT(MAXCOR)/DFLOAT(128*1024*1024)
      else
        WRITE (IW,1410)
      endif
C
C     General setup
C     -------------
C
      CALL CCDINI(DORESTART)
      EQNS = '      '
C
C     Set up symmetry tables and parallel distribution arrays
C     -------------------------------------------------------
C
C
      CALL SYMTAB (NELEC,NFROZ,NINACT,NACT,NSP,EPS1,EPS,IREPSP,
     &             IREPSPI,.FALSE.,.TRUE.,NELEC_F1,NELEC_F2)
C
      NELECT = 0
      DO I = 1, 16
         NELECT = NELECT + NELEC(I)
      ENDDO
      IF (NELECT.LE.2) DOCCSDT = .FALSE. ! Override default or user set value if this calculation makes no sense and possibly gives errors.

!     total charge of the system
      call get_charge_cc(mycharge,nelect,ncore_elec)
C
C     Write output in xml format
C
      if (doxml) then
         call xml_begin('task','type="RELCCSD"')
         call xml_begin('input')
         call xml_quantity('Number of spinors',NSP)
         call xml_quantity('Number of correlated electrons',NELECT)
         call xml_quantity('Number of virtual spinors',NSP-NELECT)
         call xml_end('input')
         call xml_begin('output')
      endif

C
      if ( maxcor .gt. 0) then
      WRITE (IW,1004) NELECT,mycharge,NSP-NELECT,CARITH,DOSORT,DOENER,
     &                DOFOPR,
     &                DOSOPR,DEBUG,TIMING,IPRNT,mword_limit,INTERFACE
      else
      WRITE (IW,1414) NELECT,mycharge,NSP-NELECT,CARITH,DOSORT,DOENER,
     &                DOFOPR,
     &                DOSOPR,DEBUG,TIMING,IPRNT,INTERFACE
      endif
      
!miro: printout info about the RelCC memory count
      write(IW,"(1X,A,10X,L1)")
     & "Leave after calculating the total memory demand :",
     &  DOMEMCOUNT

!miro: printout of changed DHOLU limit, with warning
      if (relcc_ccener_dholu_limit_set) then
        write(IW,*) 'WARNING: you have altered the DHOLU limit !'
        write(IW,*) 'new value=',relcc_ccener_dholu_limit
        write(IW,*)
     & 'Decreasing DHOLU does not solve CC convergence problems !'
      endif

C
C     Check the memory requirements of the active modules
C     ---------------------------------------------------
C
      ICALC = 1
      ALLOC_ME = .FALSE.
C
      MXCORR1 = 0

      CALL CCSETI (ICALC,MXCORR1,EPS,NSP,IREPSP,IREPSPI)

      IF (DOSORT) MXCORR = MAX(MXCORR,MXCORR1)

      IF (DOENER) THEN
         MXCORR1 = 0
         CALL CCENER (ICALC,MXCORR1,EPS)
         MXCORR = MAX(MXCORR,MXCORR1)
      ENDIF
C
      IF (DOFOPR) THEN
         MXCORR1 = 0
         CALL CCFOPR (ICALC,NSP,IREPSPI,MXCORR1,EPS)
         MXCORR = MAX(MXCORR,MXCORR1)
      ENDIF

C
      IF (DOSOPR) THEN
         MXCORR1 = 0
         CALL CCHESS (ICALC,MXCORR1,EPS,NSP,IREPSPI)
         MXCORR = MAX(MXCORR,MXCORR1)
      ENDIF
C
      IF (DOFSPC) THEN
         MXCORR1 = 0
         CALL CCFSPC(ICALC,MXCORR1,EPS1,EPS,NSP,IREPSP,IREPSPI)
         MXCORR = MAX(MXCORR,MXCORR1)
      ENDIF

!! Activate EOMCC calculation 

      IF (relcc_do_eomcc) THEN
         MXCORR1 = 0
         CALL eom_ee(ICALC,MXCORR1,eps)
         MXCORR = MAX(MXCORR,MXCORR1)
      ENDIF

C
C excitation energy code (e.g. rpa) is disabled
C     IF (DOEXCIT) THEN
C        MXCORR1 = 0
C        CALL CCEXC (ICALC,MXCORR1,EPS,NSP,IREPSPI)
C        MXCORR = MAX(MXCORR,MXCORR1)
C     ENDIF
C
C MXCORR is a cumulative maximum of individual memory requirements
C and MAXCOR is the fixed maximum of core memory available on this
C node. ===> if any subroutine would need more core memory we quit.

      WRITE (IW,1001) MXCORR
C
      CALL GETTIM(CPU2,WALL2)
      CPUI = CPUI + CPU2 - CPU1
      CPU1 = CPU2
      WALLI = WALLI + WALL2 - WALL1
      WALL1 = WALL2
      CALL FLSHFO (IW)

!MI show the total memory demand in MB/GB:
         WRITE(IW,"(/,1X,A,F15.2,A)")
     &    'Predicted RelCC memory demand: ',
     &    DFLOAT(MXCORR)/DFLOAT(2**17)," MB"
         WRITE(IW,"(1X,A,F15.3,A)")
     &    'Predicted RelCC memory demand: ',
     &    DFLOAT(MXCORR)/DFLOAT(2**27)," GB"

!MI leave RelCC calculations after memory demand counted
       if (DOMEMCOUNT) GOTO 77 ! jump to the end section
C
C     Start the calculations
C     ----------------------
C
      ICALC = 2
      ALLOC_ME = .TRUE.
C
      CALL XTIME(5,1,'Sorting of integrals          ')
      IF (DOSORT) CALL CCSETI (ICALC,MXCORR,EPS,
     &                         NSP,IREPSP,IREPSPI)

      CALL XTIME(5,2,'Sorting of integrals          ')
      CALL GETTIM(CPU2,WALL2)
      CPUS = CPUS + CPU2 - CPU1
      CPU1 = CPU2
      WALLS = WALLS + WALL2 - WALL1
      WALL1 = WALL2
      CALL FLSHFO (IW)
C
      IF (DOENER) CALL CCENER (ICALC,MXCORR,EPS)
      CALL GETTIM(CPU2,WALL2)
      CPUE = CPUE + CPU2 - CPU1
      CPU1 = CPU2
      WALLE = WALLE + WALL2 - WALL1
      WALL1 = WALL2
      CALL FLSHFO (IW)
C
      IF (DOFOPR) CALL CCFOPR (ICALC,NSP,IREPSPI,MXCORR,EPS)
      CALL GETTIM(CPU2,WALL2)
      CPUFO = CPUFO + CPU2 - CPU1
      CPU1 = CPU2
      WALLFO = WALLFO + WALL2 - WALL1
      WALL1 = WALL2
      CALL FLSHFO (IW)
C
      IF (DOSOPR) CALL CCHESS (ICALC,MXCORR,EPS,NSP,IREPSPI)
      CALL GETTIM(CPU2,WALL2)
      CPUSO = CPUSO + CPU2 - CPU1
      CPU1 = CPU2
      WALLSO = WALLSO + CPU2 - CPU1
      WALL1 = WALL2
      CALL FLSHFO (IW)
C
      IF (DOFSPC) CALL CCFSPC(ICALC,MXCORR,EPS1,EPS,NSP,IREPSP,IREPSPI)
      CALL GETTIM(CPU2,WALL2)
      CPUFS = CPUFS + CPU2 - CPU1
      CPU1 = CPU2
      WALLFS = WALLFS + CPU2 - CPU1
      WALL1 = WALL2
      CALL FLSHFO (IW)

      IF (relcc_do_eomcc) CALL eom_ee(ICALC,MXCORR,eps)
      CALL GETTIM(CPU2,WALL2)
      CPUEOM = CPUEOM + CPU2 - CPU1
      CPU1 = CPU2
      WALLEOM = WALLEOM + CPU2 - CPU1
      WALL1 = WALL2
      CALL FLSHFO (IW)

C  FINISHED; RETURN AFTER CLOSING FILES
C
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
      IF (DOCCSD.OR.DOCCSDT.OR.DOCCSDG.OR.relcc_do_eomcc) 
     &       CALL WAIO_CLOSE(ITAPT+5,JCODE)
      JCODE = 3
      IF (DOCCSD.OR.DOCCSDG.OR.relcc_do_eomcc)
     &  CALL WAIO_CLOSE(ITAPT+6,JCODE)
      IF (DOFSPC) CLOSE (IUHEFF)
C
      CALL CCSUMM
      CALL FLSHFO (IW)

! MI: jump here if only the RelCC memory count needed

 77   CONTINUE 
C
C     Print timings
C
!!      IF (.NOT.DOCCSDG)
      CALL XTIME(0,4,' RELCCSD                         ')
      CALL FLSHFO (IW)

      CALL GETTIM(CPUTOT,WALLTOT)

      CPUTOT = CPUTOT - CPU0
      CPUVAR = CPUTOT-CPUI-CPUS-CPUE-CPUFO-CPUSO-CPUFS-CPUEOM
      WALLTOT = WALLTOT-WALL0
      WALLVAR = WALLTOT-WALLI-WALLS-WALLE-WALLFO-WALLSO-WALLFS-WALLEOM
      WRITE(IW,1010) WALL0,CPU0,WALLI,CPUI,WALLS,CPUS,WALLE,CPUE,
     & WALLFO,CPUFO,WALLSO,CPUSO,WALLFS,CPUFS,WALLEOM,CPUEOM,
     & WALLVAR,CPUVAR,WALLTOT,CPUTOT


C     Print I/O statistics

      CALL WAIO_PRST (IW)
      CALL FLSHFO (IW)
C
      if (doxml) then
         call xml_end('output')
         call xml_end('task')
      endif

      RETURN
C
1000  FORMAT (/" Total memory available :",F15.2," MB",1X,F8.3," GB")
1410  FORMAT (/" Total memory available has not been determined")
1001  FORMAT (" Memory used for active modules :",
     &        T50,I15," 8-byte words")
1002  FORMAT (" Not enough memory, increase to at least",
     &        T50,I15," 8-byte words"//" Make Mword equal to ",I5)
1003  FORMAT (/" Size of 3- and 4-virtual batches :",T50,I10
     &        /" Number of passes :",T50,I10)
1414  FORMAT (//" Number of electrons :",T40,I5
     &        /" Total charge of the system :",T40,I5
     &        /" Number of virtual spinors :",T40,I5
     &        /" Complex arithmetic mode :",T44,L1
     &        /" Do integral sorting     :",T44,L1
     &        /" Do energy calculation   :",T44,L1
     &        /" Do gradient calculation :",T44,L1
     &        /" Do response calculation :",T44,L1
     &        /" Debug information       :",T44,L1
     &        /" Timing information      :",T44,L1
     &        /" Print level          :",T40,I5
     &        /" Memory limit (MWord) not set in Dirac but by the OS"
     &        /" Interface used       :",T40,A10//)
1004  FORMAT (//" Number of electrons :",T40,I5
     &        /" Total charge of the system :",T40,I5
     &        /" Number of virtual spinors :",T40,I5
     &        /" Complex arithmetic mode :",T44,L1
     &        /" Do integral sorting     :",T44,L1
     &        /" Do energy calculation   :",T44,L1
     &        /" Do gradient calculation :",T44,L1
     &        /" Do response calculation :",T44,L1
     &        /" Debug information       :",T44,L1
     &        /" Timing information      :",T44,L1
     &        /" Print level          :",T40,I5
     &        /" Memory limit (MWord) :",T40,I10
     &        /" Interface used       :",T40,A10//)
1005  FORMAT (" Size of I/O buffers :",T50,I10," kB")
1006  FORMAT (" Too much memory specified, maximum is",
     &        T50,I10," 8-byte words"//" Make Mword equal to ",I5)
1008  FORMAT (" Amount of memory used :",T50,F10.1," megawords")
1010  FORMAT(//' Timing of main modules :',T35,
     &           'Wallclock (s)',T55,'CPU on master (s)',
     &/' Before CC driver :',T30,F12.2,T55,F12.2
     &/' Initialization :',T30,F12.2,T55,F12.2
     &/' Integral sorting :',T30,F12.2,T55,F12.2
     &/' Energy calculation :',T30,F12.2,T55,F12.2
     &/' First order properties :',T30,F12.2,T55,F12.2
     &/' Second order properties :',T30,F12.2,T55,F12.2
     &/' Fock space energies :',T30,F12.2,T55,F12.2
     &/' EOMCC energies :',T30,F12.2,T55,F12.2
     &/' Untimed parts :',T30,F12.2,T55,F12.2 
     &/' Total time in CC driver :',T30,F12.0,T55,F12.2)
 1011 FORMAT (/'---< Process ',I5,' of ',I5,'----<'/)
 1015 FORMAT (/' Checking ',A,'MXCORR1,MXCORR',2I10)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SYMTAB (NELEC,NFROZ,NINACT,NACT,
     &                   NSP,EPS1,EPS,IREPSP,IREPSPI,
     &                   FOCKSP,REPORT,NELEC_F1,NELEC_F2)
C
      use spinor_indexing
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Generate symmetry pointer arrays
C
C---------------Routines called----------------------------------------
C
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "files.inc"
#include "symm.inc"
#include "freeze.inc"
#include "complex.inc"
#include "ihm.inc"
#include "ccpar.inc"
C---------------Calling variables--------------------------------------
C
      INTEGER IREPSP(NSP,4),IREPSPI(NSP,MXREP,2)
      INTEGER NELEC(16),NELEC_F1(16), NELEC_F2(16)
      INTEGER NFROZ(16),NINACT(16),NACT(16)
      REAL*8 EPS1(*),EPS(*)
      LOGICAL FOCKSP,REPORT
      INTEGER NSP
C
C---------------Local variables--------------------------------------
C
      CHARACTER*14 REPNT(16)
      INTEGER NVIRT(16)
      CHARACTER*4 REPNA_ACT(MXREP)
      INTEGER NO_ACT(MXREP),NV_ACT(MXREP)
      INTEGER NIO_ACT(MXREP),NAO_ACT(MXREP)
      INTEGER NAV_ACT(MXREP),NIV_ACT(MXREP)
      !MI: for the IH-FSCC method in linear symmetry
      INTEGER NAOPI_ACT(MXREP),NAVPI_ACT(MXREP)
      LOGICAL ASSREP
      INTEGER IOI(0:3), IVA(0:3)
      integer i,ICURFSS,ii,IJKREP,IJREP,IREP,IRP,ISP,j,JKREP,JREP
      integer KREP,NACT_REP,NREP1,NREP2

C
C---------------Executable code--------------------------------------
C
C     Initialize dimension arrays
C
      DO IREP = 1, NREP
         NE(IREP) = 0
         NO(IREP) = 0
         NV(IREP) = 0
         NC(IREP) = 0
!         NCONT(IREP) = 0
         NIO(IREP) = 0
         NIV(IREP) = 0
         NAO(IREP) = 0
         NAV(IREP) = 0
         NSO(IREP) = 0
         NSV(IREP) = 0
         MOO(IREP) = 0
         MOV(IREP) = 0
         MVO(IREP) = 0
         MVV(IREP) = 0
         NOO(IREP) = 0
         NOV(IREP) = 0
         NVO(IREP) = 0
         NVV(IREP) = 0
         NOOT(IREP) = 0
         NVVT(IREP) = 0
         NOOOT(IREP) = 0
         NVOOT(IREP) = 0
         NVVOT(IREP) = 0
         NOOVT(IREP) = 0
         NOOV(IREP) = 0
         NVVVT(IREP) = 0
         NVOO (IREP) = 0
         NVOV (IREP) = 0
         NVOV2(IREP) = 0
         NOVO (IREP) = 0
         NOVO2 (IREP) = 0
         NVVO (IREP) = 0
         NOVVT(IREP) = 0
         NOVV(IREP) = 0
         NOOOT2(IREP) = 0
         NVVVT2(IREP) = 0
C ** additional arrays for parallel triples
         NT3VVT(IREP) = 0
         NT3VVVT(IREP) = 0
      ENDDO
C
C     Convert input information to abelian subgroup
C
      ASSREP = .FALSE.
      IF (.NOT.FOCKSP) THEN
C
C     Single reference calculation : we need only to distinguish
C     between occupied and virtual, no active orbitals.
C
C     Check if we go through a manual or automatic assignment
C
C
         DO I = 1, 16
            IF (NELEC_F1(I).NE.0.OR.NELEC_F2(I).NE.0) THEN
               ASSREP = .TRUE.
            ENDIF
         ENDDO
         IF (ASSREP) THEN
C
C     Manual assignment for each abelian subgroup (Ivan Infante jan04)
C
C
C...........Loop over parent group
            DO IRP = 1, NSYMRP
               NELEC(IRP) = 0
               NVIRT(IRP) = 0
C..............Loop over Abelian subgroup
               DO IREP = 1, NREP
                  IN = 0
                  DO I = 1, NSP
                     IF (IREPSP(I,2).EQ.IREP
     &                   .AND.IREPSP(I,1).EQ.IRP) THEN
                        IN = IN + 1
                        IF (IREP.LE.(NREP/2)) THEN
C
C     In the case of a molecule with an inversion center (NFSYM=2),
C     we perform two separate assignments according to the gerade
C     (NELEC_F1) and ungerade symmetry (NELEC_F2).
C     For a molecule without inversion centre, we performe one single
C     assignment (NELEC_F1).
C
                           IF (IN.LE.NELEC_F1(IREP)) THEN
                              NELEC(IRP) = NELEC(IRP) + 1
                              NO(IREP) = NO(IREP) + 1
                              IREPSP(I,3) = NO(IREP)
                              IREPSPI(NO(IREP),IREP,1) = I
                           ELSE
                              NVIRT(IRP) = NVIRT(IRP) + 1
                              NV(IREP) = NV(IREP) + 1
                              IREPSP(I,3) = -NV(IREP)
                              IREPSPI(NV(IREP),IREP,2) = I
                           ENDIF
                        ELSE
C
C     Assignment for the ungerade symmetry (only if NFSYM=2)
C
                           IF (IN.LE.NELEC_F2(IREP-(NREP/2))) THEN
                              NELEC(IRP) = NELEC(IRP) + 1
                              NO(IREP) = NO(IREP) + 1
                              IREPSP(I,3) = NO(IREP)
                              IREPSPI(NO(IREP),IREP,1) = I
                           ELSE
                              NVIRT(IRP) = NVIRT(IRP) + 1
                              NV(IREP) = NV(IREP) + 1
                              IREPSP(I,3) = -NV(IREP)
                              IREPSPI(NV(IREP),IREP,2) = I
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
C
         ELSE  ! if (ASSREP) then ... else ... endif
         ! here: Automatic assigment for each abelian subgroup
C
            IF (.not.IFROZ_set_in_input) NFROZEN = 0
            DO IRP = 1, NSYMRP
               NVIRT(IRP) = 0
               IN  = 0
               DO I = 1, NSP
                  IF (IREPSP(I,1).EQ.IRP) THEN
                     IN = IN + 1
                     IREP = IREPSP(I,2)
                     IF (IN.LE.NELEC(IRP)) THEN
                        NO(IREP) = NO(IREP) + 1
                        IREPSP(I,3) = NO(IREP)
                        IREPSPI(NO(IREP),IREP,1) = I
!        write(*,*) 'occ:',NO(IREP),IREP,IREPSPI(NO(IREP),IREP,1)
                        IF (.not.IFROZ_set_in_input) THEN
                           IF (NO(IREP).GT.MAXFROZ)
     &                        CALL QUIT("INCREASE MAXFROZ")
                           IF (IN.LE.NFROZ(IRP)) THEN
                              NC(IREP) = NC(IREP) + 1
                              IFROZ(NO(IREP),IREP) =  NC(IREP)
                              NFROZEN = NFROZEN + 1
                           ELSE
                              IFROZ(NO(IREP),IREP) =  0
                           ENDIF
                        ELSE IF (IFROZ(NO(IREP),IREP) .GT. 0) THEN
                           NC(IREP) = NC(IREP) + 1
                        ENDIF
                     ELSE
                        NVIRT(IRP) = NVIRT(IRP) + 1
                        NV(IREP) = NV(IREP) + 1
                        IREPSP(I,3) = - NV(IREP)
                        IREPSPI(NV(IREP),IREP,2) = I
!        write(*,*) 'vir:',NV(IREP),IREP,IREPSPI(NV(IREP),IREP,2)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDIF  ! if (ASSREP) then ... else ... endif
C
         DO IREP = 1, NREP
C           The number of electrons is now identical to NO, but we need
C           two separate arrays because NO is extended in case of
C           Fock space calculations (see below).
            NE(IREP)  = NO(IREP)
C           Fill in the Fock space arrays in case we call a routine that
C           needs them.
            NIO(IREP) = NO(IREP)
            NAO(IREP) = 0
            NIV(IREP) = NV(IREP)
            NAV(IREP) = 0
            NSO(IREP) = NO(IREP)
            NSV(IREP) = NV(IREP)
         ENDDO
C
      ELSE ! Fockspace
C
C     Active spinors belong to both the occupied and virtual
C     space in Fock space calculations. We need a different sorting
C     that puts all the active occupied at the end of the occupied set
C     and all the active virtuals in the beginning. This may break the
C     original (higher pointgroup) symmetry ordering but keeps the
C     Abelian order.
C
C        Step 1 - Add the inactive occupied to the O space
C               - Add the active virtual to the V space
C
         IF (.not.IFROZ_set_in_input) NFROZEN = 0
         DO IRP = 1, NSYMRP
            NVIRT(IRP) = 0
            IN  = 0
            DO I = 1, NSP
               IF (IREPSP(I,1).EQ.IRP) THEN
                  IN = IN + 1
                  IREP = IREPSP(I,2)
                  IF (IN.LE.NELEC(IRP).AND.
     &                IN.LE.NINACT(IRP)) THEN
                     NIO(IREP) = NIO(IREP) + 1
                     NO(IREP) = NO(IREP) + 1
                     IREPSP(I,3) = NO(IREP)
                     IREPSP(I,4) = IREPSP(I,3)
                     IREPSPI(NO(IREP),IREP,1) = I
C                    Check for frozen core orbitals
                     IF (.not.IFROZ_set_in_input) THEN
                        IF (NO(IREP).GT.MAXFROZ)
     &                  CALL QUIT("INCREASE MAXFROZ")
                        IF (IN.LE.NFROZ(IRP)) THEN
                           NC(IREP) = NC(IREP) + 1
                           IFROZ(NO(IREP),IREP) =  NC(IREP)
                           NFROZEN = NFROZEN + 1
                        ELSE
                           IFROZ(NO(IREP),IREP) =  0
                        ENDIF
                     ELSE IF (IFROZ(NO(IREP),IREP) .GT. 0) THEN
                        NC(IREP) = NC(IREP) + 1
                     ENDIF
                  ENDIF
                  IF (IN.GT.NELEC(IRP).AND.
     &                IN.LE.NACT(IRP)+NINACT(IRP)) THEN
                     NVIRT(IRP) = NVIRT(IRP) + 1
                     NAV(IREP) = NAV(IREP) + 1
                     NV(IREP) = NV(IREP) + 1
                     IREPSP(I,3) = - NV(IREP)
                     IREPSPI(NV(IREP),IREP,2) = I
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
C
C        Step 2 - Add the active occupied to the O space
C               - Add the inactive virtual to the V space
C
         DO IRP = 1, NSYMRP
            IN  = 0
            DO I = 1, NSP
               IF (IREPSP(I,1).EQ.IRP) THEN
                  IN = IN + 1
                  IREP = IREPSP(I,2)
                  IF (IN.LE.NELEC(IRP).AND.
     &                IN.GT.NINACT(IRP)) THEN
                     NAO(IREP) = NAO(IREP) + 1
                     NO(IREP) = NO(IREP) + 1
                     IREPSP(I,3) = NO(IREP)
                     IREPSPI(NO(IREP),IREP,1) = I
                     IF (.not.IFROZ_set_in_input) THEN
                        IFROZ(NO(IREP),IREP) =  0 ! No frozen act. allowed
                     ELSE IF (IFROZ(NO(IREP),IREP) .GT. 0) THEN
                        NC(IREP) = NC(IREP) + 1
                     END IF
                  ENDIF
                  IF (IN.GT.NELEC(IRP).AND.
     &                IN.GT.NACT(IRP)+NINACT(IRP)) THEN
                     NVIRT(IRP) = NVIRT(IRP) + 1
                     NIV(IREP) = NIV(IREP) + 1
                     NV(IREP) = NV(IREP) + 1
                     IREPSP(I,3) = - NV(IREP)
                     IREPSP(I,4) = IREPSP(I,3)
                     IREPSPI(NV(IREP),IREP,2) = I
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
C
C        Step 3 - Add the active occupied also to the V space
C               - Add the active virtual also to the O space
C
         DO IRP = 1, NSYMRP
            IN  = 0
            DO I = 1, NSP
               IF (IREPSP(I,1).EQ.IRP) THEN
                  IN = IN + 1
                  IREP = IREPSP(I,2)
                  IF (IN.LE.NELEC(IRP).AND.
     &                IN.GT.NINACT(IRP)) THEN
                     NV(IREP) = NV(IREP) + 1
                     IREPSP(I,4) = - NV(IREP)
                     IREPSPI(NV(IREP),IREP,2) = I
                  ENDIF
                  IF (IN.GT.NELEC(IRP).AND.
     &                IN.LE.NACT(IRP)+NINACT(IRP)) THEN
                     NO(IREP) = NO(IREP) + 1
                     IREPSP(I,4) = NO(IREP)
                     IREPSPI(NO(IREP),IREP,1) = I
                     IF (.not.IFROZ_set_in_input) THEN
                        IFROZ(NO(IREP),IREP) =  0 ! No frozen act. allowed
                     ELSE IF (IFROZ(NO(IREP),IREP) .GT. 0) THEN
                        NC(IREP) = NC(IREP) + 1
                     END IF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
C
C        The total numbers of occupied/virtual, regardless of act/inact
C
         DO IREP = 1, NREP
            NSO(IREP) = NIO(IREP) + NAO(IREP)
            NSV(IREP) = NIV(IREP) + NAV(IREP)
            NE(IREP)  = NSO(IREP)
         ENDDO
      ENDIF
C
      IF (REPORT.AND.(.NOT.FOCKSP)) THEN
!        NSYMRP - DIRAC initial symmetry (from RECC input)
            WRITE (IW,1000) (REPN(I),I=1,NSYMRP)
            WRITE (IW,1002) 'occupied',(NELEC(I),I=1,NSYMRP)
            WRITE (IW,1002) 'virtual ',(NVIRT(I),I=1,NSYMRP)
         IF (NREP.GT.8) THEN
C           Write only active irreps (NO or NV > 0)
            NACT_REP = 0
!           NREP - number of irreps in abelian subgroup
            DO I = 1, NREP
               IF (NO(I)+NV(I).GT.0) THEN
                  NACT_REP = NACT_REP + 1
                  NO_ACT(NACT_REP) = NO(I)
                  NV_ACT(NACT_REP) = NV(I)
                  REPNA_ACT(NACT_REP) = REPNA(I)
               ENDIF
            ENDDO
            NREP1 = 1
            NREP2 = MIN(NACT_REP,8)
  10        WRITE (IW,1001) (REPNA_ACT(I),I=NREP1,NREP2)
            WRITE (IW,1002) 'occupied',(NO_ACT(I),I=NREP1,NREP2)
            WRITE (IW,1002) 'virtual ',(NV_ACT(I),I=NREP1,NREP2)
            NREP1 = NREP1 + 8
            NREP2 = MIN(NACT_REP,NREP2+8)
            IF (NREP1.LE.NACT_REP) GOTO 10
         ELSE
            WRITE (IW,1001) (REPNA(I),I=1,NREP)
            WRITE (IW,1002) 'occupied',(NO(I),I=1,NREP)
            WRITE (IW,1002) 'virtual ',(NV(I),I=1,NREP)
         ENDIF
      ENDIF
C
      IF (NFROZEN.GT.0.AND.REPORT) THEN
         WRITE (IW,1003)
         WRITE (IW,1000) (REPN(I),I=1,NSYMRP)
         WRITE (IW,1002) 'frozen',(NFROZ(I),I=1,NSYMRP)
         WRITE (IW,1001) (REPNA(I),I=1,NREP)
         WRITE (IW,1002) 'frozen',(NC(I),I=1,NREP)
      ENDIF
C
      IF (FOCKSP.AND.REPORT) THEN
         WRITE (IW,1004)
         WRITE (IW,1000) (REPN(I),I=1,NSYMRP)
         WRITE (IW,1002) 'inactive occupied',(NINACT(I),I=1,NSYMRP)
         WRITE (IW,1002) 'active occupied  ',
     &                   (NELEC(I)-NINACT(I),I=1,NSYMRP)
         WRITE (IW,1002) 'active virtual   ',
     &                   (NINACT(I)+NACT(I)-NELEC(I),I=1,NSYMRP)
         WRITE (IW,1002) 'inactive virtual ',(NELEC(I)+NVIRT(I)-
     &                    NINACT(I)-NACT(I),I=1,NSYMRP)

!        ... write out reprezentations in Abelian subgroup ...
         IF (NREP.GT.8) THEN
C           Write only active irreps (NO or NV > 0)
            NACT_REP = 0
            DO I = 1, NREP
               IF (NIO(I)+NAO(I)+NAV(I)+NIV(I).GT.0) THEN
                !... move nonzero values to new arrays ...
                  NACT_REP = NACT_REP + 1
                  NIO_ACT(NACT_REP) = NIO(I)
                  NAO_ACT(NACT_REP) = NAO(I)
                  NAV_ACT(NACT_REP) = NAV(I)
                  NIV_ACT(NACT_REP) = NIV(I)
                  REPNA_ACT(NACT_REP) = REPNA(I)
               ENDIF
            ENDDO
            NREP1 = 1
            NREP2 = MIN(NACT_REP,8)
  11      WRITE (IW,1001) (REPNA_ACT(I),I=NREP1,NREP2)
          WRITE (IW,1002) 'inactive occupied',(NIO_ACT(I),I=NREP1,NREP2)
          WRITE (IW,1002) 'active occupied',(NAO_ACT(I),I=NREP1,NREP2)
          WRITE (IW,1002) 'active virtual ',(NAV_ACT(I),I=NREP1,NREP2)
          WRITE (IW,1002) 'inactive virtual ',(NIV_ACT(I),I=NREP1,NREP2)
            NREP1 = NREP1 + 8
            NREP2 = MIN(NACT_REP,NREP2+8)
            IF (NREP1.LE.NACT_REP) GOTO 11
          ELSE
            WRITE (IW,1001) (REPNA(I),I=1,NREP)
            WRITE (IW,1002) 'inactive occupied',(NIO(I),I=1,NREP)
            WRITE (IW,1002) 'active occupied  ',(NAO(I),I=1,NREP)
            WRITE (IW,1002) 'active virtual   ',(NAV(I),I=1,NREP)
            WRITE (IW,1002) 'inactive virtual ',(NIV(I),I=1,NREP)
         ENDIF
      ENDIF
C
C     Set up dimension arrays
C
C     Nxx(IREP) gives the total number of |ket>|ket> index pairs that have symmetry IREP
C     Mxx(IREP) gives the total number of <bra |ket> index pairs that have symmetry IREP
C     II arrays provide pointers to the start of irrep combinations for |ket>|ket> or <bra|<bra| index pairs
C     JJ arrays provide pointers to the start of irrep combinations for <bra|ket> index pairs
C
      DO JREP = 1, NREP
         DO IREP = 1, NREP
            IJREP = MULTB(IREP,JREP,1)
            IIOO(IREP,JREP) = NOO(IJREP)
            IIOV(IREP,JREP) = NOV(IJREP)
            IIVO(IREP,JREP) = NVO(IJREP)
            IIVV(IREP,JREP) = NVV(IJREP)
            NOO(IJREP) = NOO(IJREP) + NO(IREP)*NO(JREP)
            NOV(IJREP) = NOV(IJREP) + NO(IREP)*NV(JREP)
            NVO(IJREP) = NVO(IJREP) + NV(IREP)*NO(JREP)
            NVV(IJREP) = NVV(IJREP) + NV(IREP)*NV(JREP)
            IJREP = MULTB(IREP,JREP,2)
            JJOO(IREP,JREP) = MOO(IJREP)
            JJOV(IREP,JREP) = MOV(IJREP)
            JJVO(IREP,JREP) = MVO(IJREP)
            JJVV(IREP,JREP) = MVV(IJREP)
            MOO(IJREP) = MOO(IJREP) + NO(IREP)*NO(JREP)
            MOV(IJREP) = MOV(IJREP) + NO(IREP)*NV(JREP)
            MVO(IJREP) = MVO(IJREP) + NV(IREP)*NO(JREP)
            MVV(IJREP) = MVV(IJREP) + NV(IREP)*NV(JREP)
         ENDDO
      ENDDO
C
C     Same M/N/II/JJ convention as used above, but now for triangular (I>J) arrays
C
      DO JREP = 1, NREP
         IJREP = MULTB(JREP,JREP,1)
         IIOOT(JREP,JREP) = NOOT(IJREP)
         IIVVT(JREP,JREP) = NVVT(IJREP)
         NOOT(IJREP) = NOOT(IJREP) + NO(JREP)*(NO(JREP)-1)/2
         NVVT(IJREP) = NVVT(IJREP) + NV(JREP)*(NV(JREP)-1)/2
         DO IREP = JREP + 1, NREP
            IJREP = MULTB(IREP,JREP,1)
            IIOOT(IREP,JREP) = NOOT(IJREP)
            IIVVT(IREP,JREP) = NVVT(IJREP)
            NOOT(IJREP) = NOOT(IJREP) + NO(IREP)*NO(JREP)
            NVVT(IJREP) = NVVT(IJREP) + NV(IREP)*NV(JREP)
         ENDDO
      ENDDO
C
C     Dimension arrays for 3-index tensors, formed by combining a 2-index array with a 1-index array
C     KK arrays are symmetry pointers for tensors in which the leftmost indices are first grouped together
C     LL arrays are symmetry pointers for tensors in which the rightmost indices are first grouped together
C
      DO KREP = 1, NREP
         DO IJREP = 1, NREP
!
!  discussion with luuk
!  for consistency IJREP is a bra and the order in the MULTB
!  array may be swapped to maintain consistency with the definitions
!
!  original:
!
            IJKREP = MULTB(KREP,IJREP+NREP,2)
!
!  possible alternative in case of totally symmetric products:
!           IJKREP = MULTB(IJREP+NREP,KREP,2)
!
C           Here we use that (i) KLREP = IJREP, (ii) LREP = MULTB(KREP,KLREP,2), (iii) Ordering should be
C           such that we can contract with an array to the right in which the leftmost index is ordered on LREP.

            KKOOOT(IJREP,KREP) = NOOOT(IJKREP)
            KKVVOT(IJREP,KREP) = NVVOT(IJKREP)
            KKOOVT(IJREP,KREP) = NOOVT(IJKREP)
            KKVVVT(IJREP,KREP) = NVVVT(IJKREP)
            KKVOO (IJREP,KREP) = NVOO (IJKREP)
            KKVOV (IJREP,KREP) = NVOV (IJKREP)
            KKOVO (IJREP,KREP) = NOVO2(IJKREP)
            KKOVV (IJREP,KREP) = NOVV (IJKREP)

            NOOOT(IJKREP) = NOOOT(IJKREP) + NOOT(IJREP)*NO(KREP)
            NVVOT(IJKREP) = NVVOT(IJKREP) + NVVT(IJREP)*NO(KREP)
            NOOVT(IJKREP) = NOOVT(IJKREP) + NOOT(IJREP)*NV(KREP)
            NVVVT(IJKREP) = NVVVT(IJKREP) + NVVT(IJREP)*NV(KREP)
            NVOO (IJKREP) = NVOO (IJKREP) + NVO (IJREP)*NO(KREP)
            NOVV (IJKREP) = NOVV (IJKREP) + NOV (IJREP)*NV(KREP)
            NOVO2(IJKREP) = NOVO2(IJKREP) + NOV (IJREP)*NO(KREP)
            NVOV (IJKREP) = NVOV (IJKREP) + NVO (IJREP)*NV(KREP)
         ENDDO
      ENDDO
C
      DO IJREP = 1, NREP
         DO KREP = 1, NREP
            IJKREP = MULTB(KREP,IJREP+NREP,2)
            LLOVO (KREP,IJREP) = NOVO (IJKREP)
            LLOOV (KREP,IJREP) = NOOV (IJKREP)
            LLVVO (KREP,IJREP) = NVVO (IJKREP)
            LLVOOT(KREP,IJREP) = NVOOT(IJKREP)
            LLVOV (KREP,IJREP) = NVOV (IJKREP)
            LLOVVT(KREP,IJREP) = NOVVT(IJKREP)
            LLOOOT(KREP,IJREP) = NOOOT2(IJKREP)
            LLVVVT(KREP,IJREP) = NVVVT2(IJKREP)
            NOVO (IJKREP) = NOVO (IJKREP) + NVO(IJREP)*NO(KREP)
            NOOV (IJKREP) = NOOV (IJKREP) + NOV(IJREP)*NO(KREP)
            NVVO (IJKREP) = NVVO (IJKREP) + NVO(IJREP)*NV(KREP)
            NVOOT(IJKREP) = NVOOT(IJKREP) + NOOT(IJREP)*NV(KREP)
            NVOV2(IJKREP) = NVOV2(IJKREP) + NOV (IJREP)*NV(KREP)
            NOVVT(IJKREP) = NOVVT(IJKREP) + NVVT(IJREP)*NO(KREP)
            NOOOT2(IJKREP) = NOOOT2(IJKREP) + NOOT(IJREP)*NO(KREP)
            NVVVT2(IJKREP) = NVVVT2(IJKREP) + NVVT(IJREP)*NV(KREP)
         ENDDO
      ENDDO
C
C set up number of Fock matrix elements NFOO, NFVO, NFVV
C and maximum number of amplitudes NDIMT1, NDIMT2
C
      NFOO = MOO(1)
      NFVO = MVO(1)
      NFVV = MVV(1)
      NDIMT1 = 0
      NDIMT2 = 0
      DO IREP = 1, NREP
         NDIMT1 = NDIMT1 + NV(IREP)*NO(IREP)
         NDIMT2 = NDIMT2 + NVVT(IREP)*NOOT(IREP)
      ENDDO
C
C     Set up pointer arrays : IXXXX for <IJ,KL:KLREP> ordening
C                             JXXXX for (IK|JL:JLREP) ordening
C                             KXXXX for <IJK,L:LREP> ordening
C                             LXXXX for <I,JKL:JKLREP> ordening
C
      IO(1) = 0
      IV(1) = 0
C
      DO IREP = 1, NREP
         IO(IREP+1) = IO(IREP) + NO(IREP)
         IV(IREP+1) = IV(IREP) + NV(IREP)
      ENDDO
C
      IOO(1) = 0
      IVO(1) = 0
      IVV(1) = 0
      IOOOOTT(1) = 0
      IOOOO(1) = 0
      IVOOOT (1) = 0
      IOVOOT (1) = 0
      IVVOOTT(1) = 0
      IVOVO  (1) = 0
      IOVVO  (1) = 0
      IVOOV  (1) = 0 ! new index
      IVOVVT (1) = 0
      IOVVVT (1) = 0
      IVVVVTT(1) = 0
      IVVOO  (1) = 0
      IVVOOT (1) = 0
C
C setting up the totally symmetric XXXX arrays. IMPORTANT:
C if we look at <VV||VV> the product of all four ireps has to be
C totally symmetric. due to construction the NXX(T) arrays are either
C starred or unstarred and the product NXX(T) * NXX(T) always corresponds
C to GAMMA* * GAMMA = A only in ABELIAN symmetries !!!
C so a further distinction and check for total symmetry is not necessary !!
C ==> therefore the IXXXX(TT) arrays contain the offsets of nonvanishing
C totally symmetric integrals.
C
      DO IREP = 1, NREP
         IOO    (IREP+1) = IOO    (IREP) + NO  (IREP)*NO  (IREP)
         IVO    (IREP+1) = IVO    (IREP) + NV  (IREP)*NO  (IREP)
         IVV    (IREP+1) = IVV    (IREP) + NV  (IREP)*NV  (IREP)
         IOOOOTT(IREP+1) = IOOOOTT(IREP) + NOOT(IREP)*NOOT(IREP)
! new index. will go away.
         IOOOO  (IREP+1) = IOOOO  (IREP) + NOO(IREP)*NOO(IREP)
         IVOOOT (IREP+1) = IVOOOT (IREP) + NVO (IREP)*NOOT(IREP)
         IOVOOT (IREP+1) = IOVOOT (IREP) + NOV (IREP)*NOOT(IREP)
         IVVOOTT(IREP+1) = IVVOOTT(IREP) + NVVT(IREP)*NOOT(IREP)
         IVOVO  (IREP+1) = IVOVO  (IREP) + NVO (IREP)*NVO (IREP)
         IOVVO  (IREP+1) = IOVVO  (IREP) + NOV (IREP)*NVO (IREP)
         IVOOV  (IREP+1) = IVOOV  (IREP) + NVO (IREP)*NOV (IREP) !new index
         IVOVVT (IREP+1) = IVOVVT (IREP) + NVO (IREP)*NVVT(IREP)
         IOVVVT (IREP+1) = IOVVVT (IREP) + NOV (IREP)*NVVT(IREP) !new index
         IVVVVTT(IREP+1) = IVVVVTT(IREP) + NVVT(IREP)*NVVT(IREP)
         IVVOO  (IREP+1) = IVVOO  (IREP) + NVV (IREP)*NOO (IREP)
         IVVOOT (IREP+1) = IVVOOT (IREP) + NVV (IREP)*NOOT(IREP)
      ENDDO
C
      NV1 = IOOOOTT(NREP+1)
      NV2 = IVOOOT (NREP+1)
      NV3 = IVVOOTT(NREP+1)
      NV4 = IVOVO  (NREP+1)
      NV5 = IVOVVT (NREP+1)
      NV6 = IVVVVTT(NREP+1)
C

      JVOXX(1) = 0
      JOOOO(1) = 0
      JOOOV(1) = 0
      JOOVO(1) = 0
      JOOVV(1) = 0
      JOVOO(1) = 0
      JOVOV(1) = 0
      JOVVO(1) = 0
      JOVVV(1) = 0
      JVOOO(1) = 0
      JVVOO(1) = 0
      JVVOV(1) = 0
      JVOVO(1) = 0
      JVOOV(1) = 0
      JVVVO(1) = 0
      JVVVOI4(1)  = 0
C
C ** additional RELADC arrays
C

      J2OOOV(1) = 0
      J2OOVO(1) = 0
      J2OOVV(1) = 0
      J2OVOO(1) = 0
      J2OVOV(1) = 0
      J2OVVO(1) = 0
      J2OVVV(1) = 0

      J2VOVO(1) = 0
      J2VOOV(1) = 0
      J2OOOO(1) = 0
      J2OOVO(1) = 0
      J2VOOO(1) = 0
      J2VVOO(1) = 0
      J2VVOV(1) = 0
      J2VOVV(1) = 0
      J2VVVO(1) = 0
C
C filling
C
      DO IREP = 1, NREP
         JVOXX(IREP+1) = JVOXX(IREP) + MVO(IREP)
         JOOOO(IREP+1) = JOOOO(IREP) + MOO(IREP)*MOO(IREP)
         JOOOV(IREP+1) = JOOOV(IREP) + MOO(IREP)*MOV(IREP)
         JOOVO(IREP+1) = JOOVO(IREP) + MOO(IREP)*MVO(IREP)
         JOOVV(IREP+1) = JOOVO(IREP) + MOO(IREP)*MVV(IREP)
         JOVOO(IREP+1) = JOVOO(IREP) + MOV(IREP)*MOO(IREP)
         JOVOV(IREP+1) = JOVOV(IREP) + MOV(IREP)*MOV(IREP)
         JOVVO(IREP+1) = JOVVO(IREP) + MOV(IREP)*MVO(IREP)
         JOVVV(IREP+1) = JOVVV(IREP) + MOV(IREP)*MVV(IREP)
         JVOOO(IREP+1) = JVOOO(IREP) + MVO(IREP)*MOO(IREP)
         JVVOO(IREP+1) = JVVOO(IREP) + MVV(IREP)*MOO(IREP)
         JVOVO(IREP+1) = JVOVO(IREP) + MVO(IREP)*MVO(IREP)
         JVOOV(IREP+1) = JVOOV(IREP) + MVO(IREP)*MOV(IREP)
         JVVVO(IREP+1) = JVVVO(IREP) + MVV(IREP)*MVO(IREP)
         JVVOV(IREP+1) = JVVOV(IREP) + MVV(IREP)*MOV(IREP)
         JVVVOI4(IREP+1) = JVVVOI4(IREP) + MVV(IREP)*MVO(IREP) !new index to avoid integer8
         JREP = MULTB(IREP+NREP,1+NREP,2)

         J2OOOV(IREP+1) = J2OOOV(IREP) + MOO(IREP)*MOV(JREP)
         J2OOVO(IREP+1) = J2OOVO(IREP) + MOO(IREP)*MVO(JREP)
         J2OOVV(IREP+1) = J2OOVO(IREP) + MOO(IREP)*MVV(JREP)
         J2OVOO(IREP+1) = J2OVOO(IREP) + MOV(IREP)*MOO(JREP)
         J2OVOV(IREP+1) = J2OVOV(IREP) + MOV(IREP)*MOV(JREP)
         J2OVVO(IREP+1) = J2OVVO(IREP) + MOV(IREP)*MVO(JREP)
         J2OVVV(IREP+1) = J2OVVV(IREP) + MOV(IREP)*MVV(JREP)
         J2VOVO(IREP+1) = J2VOVO(IREP) + MVO(IREP)*MVO(JREP)
         J2VOOV(IREP+1) = J2VOOV(IREP) + MVO(IREP)*MOV(JREP)
         J2OOOO(IREP+1) = J2OOOO(IREP) + MOO(IREP)*MOO(JREP)
         J2VOOO(IREP+1) = J2VOOO(IREP) + MVO(IREP)*MOO(JREP)
         J2VVOO(IREP+1) = J2VVOO(IREP) + MVV(IREP)*MOO(JREP)
         J2VVOV(IREP+1) = J2VVOV(IREP) + MVV(IREP)*MOV(JREP)
         J2VOVV(IREP+1) = J2VOVV(IREP) + MVO(IREP)*MVV(JREP)
         J2VVVO(IREP+1) = J2VVVO(IREP) + MVV(IREP)*MVO(JREP)
      ENDDO
      NDIMX = JVOXX(NREP+1)
C
      KVOOV (1) = 0
      KVOVV (1) = 0
      KOVOO(1)  = 0
      KOVOV(1)  = 0
      KOVVOT(1) = 0
      KVOVO (1) = 0
      KVOOO (1) = 0
      KOOOOT(1) = 0
      KVVOOT(1) = 0
      KOOVOT(1) = 0
      KOOOVT(1) = 0
      KOOVVT(1) = 0
      KVVVOT(1) = 0
      KVVOVT(1) = 0
      DO IREP = 1, NREP
         KVOOV  (IREP+1) = KVOOV  (IREP) + NVOO (IREP)*NV(IREP)
         KOVOO  (IREP+1) = KOVOO  (IREP) + NOVO (IREP)*NO(IREP)
         KOVOV  (IREP+1) = KOVOO  (IREP) + NOVO (IREP)*NV(IREP)
         KOVVOT (IREP+1) = KOVVOT (IREP) + NOVVT (IREP)*NO(IREP)
         KVOVV  (IREP+1) = KVOVV  (IREP) + NVOV (IREP)*NV(IREP)
         KVOVO  (IREP+1) = KVOVO  (IREP) + NVOV (IREP)*NO(IREP)
         KVOOO  (IREP+1) = KVOOO  (IREP) + NVOO (IREP)*NO(IREP)
         KOOOOT (IREP+1) = KOOOOT (IREP) + NOOOT(IREP)*NO(IREP)
         KVVOOT (IREP+1) = KVVOOT (IREP) + NVVOT(IREP)*NO(IREP)
         KOOVOT (IREP+1) = KOOVOT (IREP) + NOOVT(IREP)*NO(IREP)
         KOOOVT (IREP+1) = KOOOVT (IREP) + NOOOT(IREP)*NV(IREP)
         KOOVVT (IREP+1) = KOOVVT (IREP) + NVOOT(IREP)*NV(IREP)
         KVVOVT (IREP+1) = KVVOVT (IREP) + NVVOT(IREP)*NV(IREP)
         KVVVOT (IREP+1) = KVVVOT (IREP) + NVVVT(IREP)*NO(IREP)
      ENDDO
C
      LOOVO (1) = 0
      LOVVO (1) = 0
      LOVOV (1) = 0
      LVOVO (1) = 0
      LVVVO (1) = 0
      LOOOOT(1) = 0
      LVVOOT(1) = 0
      LOVOOT(1) = 0
      LVOOOT(1) = 0
      LOVVOT(1) = 0
      LVOOVT(1) = 0 !new index
C ** additional RELADC arrays
      LOOVVT(1) = 0
      LVOVVT(1) = 0

      DO IREP = 1, NREP
         LOOVO  (IREP+1) = LOOVO  (IREP) + NO(IREP)*NOVO (IREP)
         LVOVO  (IREP+1) = LVOVO  (IREP) + NV(IREP)*NOVO (IREP)
         LOVVO  (IREP+1) = LOVVO  (IREP) + NO(IREP)*NVVO(IREP)
         LOVOV  (IREP+1) = LOVOV  (IREP) + NO(IREP)*NVOV(IREP)
         LVVVO  (IREP+1) = LVVVO  (IREP) + NV(IREP)*NVVO(IREP)
         LOOOOT (IREP+1) = LOOOOT (IREP) + NO(IREP)*NOOOT(IREP)
         LVVOOT (IREP+1) = LVVOOT (IREP) + NV(IREP)*NVOOT(IREP)
         LOVOOT (IREP+1) = LOVOOT (IREP) + NO(IREP)*NVOOT(IREP)
         LVOOOT (IREP+1) = LVOOOT (IREP) + NV(IREP)*NOOOT(IREP)
         LOVVOT (IREP+1) = LOVVOT (IREP) + NO(IREP)*NVVOT(IREP)

         LVOOVT (IREP+1) = LVOOVT (IREP) + NV(IREP)*NVOOT(IREP) !new index
         LOOVVT (IREP+1) = LOOVVT (IREP) + NO(IREP)*NOVVT(IREP)
         LVOVVT (IREP+1) = LVOVVT (IREP) + NV(IREP)*NOVVT(IREP)
      ENDDO
C
C     3-index arrays are used in t3corr : make index arrays
C
      DO IJKREP = 1, NREP
       IVVVT(1,IJKREP) = 0
       DO JKREP = 1, NREP
        IREP = MULTB(JKREP+NREP,IJKREP,2)
        IVVVT(JKREP+1,IJKREP)=IVVVT(JKREP,IJKREP)+NV(IREP)*NVVT(JKREP)
       ENDDO
       JVVVT(1,IJKREP) = 0
       DO KREP = 1, NREP
        IJREP = MULTB(KREP,IJKREP,2)
        JVVVT(KREP+1,IJKREP)=JVVVT(KREP,IJKREP)+NVVT(IJREP)*NV(KREP)
       ENDDO
      ENDDO
      IVVVTMX = 1
      DO IJKREP = 1, NREP
      IF (IVVVT(NREP+1,IJKREP).GT.IVVVTMX) IVVVTMX=IVVVT(NREP+1,IJKREP)
      ENDDO

C ** here all symm arrays are set up and we can call the
C ** parallel distribution list generator

      CALL DISLIST

C ** set up the remaining arrays for the parallel triples
C ** Note that these arrays are also used in the serial code !!!

      DO IREP=1,NREP
        NT3VVT(IREP)=IDIST(2,2,IREP)-IDIST(1,2,IREP)
      ENDDO

      DO IJREP = 1,NREP
        DO KREP = 1,NREP
          IJKREP = MULTB(KREP,IJREP+NREP,2)
          LLT3VVVT(KREP,IJREP)=NT3VVVT(IJKREP)
          NT3VVVT(IJKREP)=NT3VVVT(IJKREP) + NT3VVT(IJREP) * NV(KREP)
        ENDDO
      ENDDO
C
      KT3VVVOT(1)=0
      DO IREP=1,NREP
        KT3VVVOT(IREP+1) = KT3VVVOT(IREP) + NT3VVVT(IREP)*NO(IREP)
      ENDDO

      KT3VVVVT(1)=0
      DO IREP=1,NREP
        KT3VVVVT(IREP+1) = KT3VVVVT(IREP) + NT3VVVT(IREP)*NV(IREP)
      ENDDO

C
C----------------------------------------------
C----- F O C K   S P A C E   A R R A Y S ------
C----------------------------------------------
C
C     Make the dimension and pointer arrays need in Fock space calcns.
C     Like in the CCSD(T) calculations we calculate only the arrays
C     needed in the actual sorting steps.
C
      IF (FOCKSP) THEN
C
         DO IREP = 1, NREP
            NIVIO(IREP) = 0
            NAVAO(IREP) = 0
            NAVSO (IREP) = 0
            NSOAV (IREP) = 0
            NAOV (IREP) = 0
            MAOAV(IREP) = 0
            NIOIOT(IREP) = 0
            NAOAOT(IREP) = 0
            NAVAVT(IREP) = 0
            NIVIVT(IREP) = 0
            NSVSVT(IREP) = 0
            NSOSOT(IREP) = 0
            NIVIVIOT(IREP) = 0
            NIOIOIVT(IREP) = 0
            NSVSVOT(IREP) = 0
            NSVSVAVT(IREP) = 0
            NVSOSOT(IREP) = 0
            NAOSOSOT(IREP) = 0
         ENDDO
C
         DO JREP = 1, NREP
            DO IREP = 1, NREP
               IJREP = MULTB(IREP,JREP,1)
               IISOAV(IREP,JREP)=NSOAV(IJREP)
               NIVIO(IJREP) = NIVIO(IJREP) + NIV(IREP)*NIO(JREP)
               NAVAO(IJREP) = NAVAO(IJREP) + NAV(IREP)*NAO(JREP)
               NAVSO(IJREP) = NAVSO(IJREP) + NAV(IREP)*NSO(JREP)
               NSOAV(IJREP) = NSOAV(IJREP) + NSO(IREP)*NAV(JREP)
               NAOV(IJREP)  = NAOV(IJREP)  + NAO(IREP)*NV(JREP)
               IJREP = MULTB(IREP,JREP,2)
               JJAOAV(IREP,JREP) = MAOAV(IJREP)
               MAOAV(IJREP) = MAOAV(IJREP) + NAO(IREP)*NAV(JREP)
            ENDDO
         ENDDO
C
         DO JREP = 1, NREP
            IJREP = MULTB(JREP,JREP,1)
            NIOIOT(IJREP) = NIOIOT(IJREP) + NIO(JREP)*(NIO(JREP)-1)/2
            NAOAOT(IJREP) = NAOAOT(IJREP) + NAO(JREP)*(NAO(JREP)-1)/2
            NAVAVT(IJREP) = NAVAVT(IJREP) + NAV(JREP)*(NAV(JREP)-1)/2
            NIVIVT(IJREP) = NIVIVT(IJREP) + NIV(JREP)*(NIV(JREP)-1)/2
            NSVSVT(IJREP) = NSVSVT(IJREP) + NSV(JREP)*(NSV(JREP)-1)/2
            NSOSOT(IJREP) = NSOSOT(IJREP) + NSO(JREP)*(NSO(JREP)-1)/2
            DO IREP = JREP + 1, NREP
               IJREP = MULTB(IREP,JREP,1)
               NIOIOT(IJREP) = NIOIOT(IJREP) + NIO(IREP)*NIO(JREP)
               NAOAOT(IJREP) = NAOAOT(IJREP) + NAO(IREP)*NAO(JREP)
               NAVAVT(IJREP) = NAVAVT(IJREP) + NAV(IREP)*NAV(JREP)
               NIVIVT(IJREP) = NIVIVT(IJREP) + NIV(IREP)*NIV(JREP)
               NSVSVT(IJREP) = NSVSVT(IJREP) + NSV(IREP)*NSV(JREP)
               NSOSOT(IJREP) = NSOSOT(IJREP) + NSO(IREP)*NSO(JREP)
            ENDDO
         ENDDO
C
         DO KREP = 1, NREP
          DO IJREP = 1, NREP
           IJKREP = MULTB(KREP,IJREP+NREP,2)
           KKIVIVIOT(IJREP,KREP) = NIVIVIOT(IJKREP)
           KKIOIOIVT(IJREP,KREP) = NIOIOIVT(IJKREP)
           KKSVSVAVT(IJREP,KREP) = NSVSVAVT(IJKREP)
           KKSVSVOT(IJREP,KREP) = NSVSVOT(IJKREP)
           NIVIVIOT(IJKREP) = NIVIVIOT(IJKREP) + NIVIVT(IJREP)*NIO(KREP)
           NIOIOIVT(IJKREP) = NIOIOIVT(IJKREP) + NIOIOT(IJREP)*NIV(KREP)
           NSVSVAVT(IJKREP) = NSVSVAVT(IJKREP) + NSVSVT(IJREP)*NAV(KREP)
           NSVSVOT(IJKREP) = NSVSVOT(IJKREP) + NSVSVT(IJREP)*NO(KREP)
          ENDDO
         ENDDO
C
         DO IJREP = 1, NREP
          DO KREP = 1, NREP
           IJKREP = MULTB(KREP,IJREP+NREP,2)
           LLVSOSOT(KREP,IJREP) = NVSOSOT(IJKREP)
           LLAOSOSOT(KREP,IJREP) = NAOSOSOT(IJKREP)
           NVSOSOT(IJKREP) = NVSOSOT(IJKREP)+NSOSOT(IJREP)*NV(KREP)
           NAOSOSOT(IJKREP) = NAOSOSOT(IJKREP)+NSOSOT(IJREP)*NAO(KREP)
          ENDDO
         ENDDO
C
         JAOAVOV(1) = 0
         JAOAVAOAV(1) = 0
         DO IREP = 1, NREP
            JAOAVOV(IREP+1)   = JAOAVOV(IREP)+MAOAV(IREP)*MOV(IREP)
            JAOAVAOAV(IREP+1) = JAOAVAOAV(IREP)+MAOAV(IREP)*MAOAV(IREP)
         ENDDO
C
         IAVSO(1) = 0
         IAVAV(1) = 0
         IAOAO(1) = 0
         KIVIVIOAVT(1) = 0
         KIOIOIVAOT(1) = 0
         KVVOAVT(1) = 0
         KSVSVAVAVT(1) = 0
         DO IREP = 1, NREP
          IAVSO(IREP+1) = IAVSO(IREP) + NAV(IREP) * NSO(IREP)
          IAVAV(IREP+1) = IAVAV(IREP) + NAV(IREP) * NAV(IREP)
          IAOAO(IREP+1) = IAOAO(IREP) + NAO(IREP) * NAO(IREP)
          KIVIVIOAVT(IREP+1) = KIVIVIOAVT(IREP)+NIVIVIOT(IREP)*NAV(IREP)
          KIOIOIVAOT(IREP+1) = KIOIOIVAOT(IREP)+NIOIOIVT(IREP)*NAO(IREP)
          KVVOAVT(IREP+1) = KVVOAVT(IREP)+NVVOT(IREP)*NAV(IREP)
          KSVSVAVAVT(IREP+1) = KSVSVAVAVT(IREP)+NSVSVAVT(IREP)*NAV(IREP)
         ENDDO
C
         LAOVOOT(1) = 0
         LAOAOSOSOT(1) = 0
         DO IREP = 1, NREP
          LAOVOOT(IREP+1) = LAOVOOT(IREP)+NAO(IREP)*NVOOT(IREP)
          LAOAOSOSOT(IREP+1) = LAOAOSOSOT(IREP)+NAO(IREP)*NAOSOSOT(IREP)
         ENDDO

C        Set dimension of effective Hamiltonian for Fock space

         DO IREP = 1, NREP
            NFS_HE(IREP,1) = 0
            NFS_HE(IREP,2) = NAV(IREP)
            NFS_HE(IREP,3) = NAO(IREP)
            NFS_HE(IREP,4) = MAOAV(IREP)
            NFS_HE(IREP,5) = NAVAVT(IREP)
            NFS_HE(IREP,6) = NAOAOT(IREP)
         ENDDO

C        Set offset for addressing symmetry blocks of effective Hamiltonian (IFS_HE)
C        or its eigenvalues (IFS_EV)

         DO ICURFSS = 1, 6
            IFS_HE(IREP,ICURFSS) = 0
            IFS_HE2(IREP,ICURFSS) = 0
            DO IREP = 1, NREP
               IFS_HE(IREP+1,ICURFSS) = IFS_HE(IREP,ICURFSS) +
     &                                  NFS_HE(IREP,ICURFSS)
               IFS_HE2(IREP+1,ICURFSS) = IFS_HE2(IREP,ICURFSS) +
     &                                   NFS_HE(IREP,ICURFSS)**2
            ENDDO
         ENDDO
C
      ENDIF
C
C     Reorder the spinor eigenvalues
C
      DO ISP = 1, NSP
         IREP = IREPSP(ISP,2)
         II = IREPSP(ISP,3)
         IF (II.GT.0) THEN
            I = IO(IREP) + II
         ELSE
            I = IO(NREP+1) + IV(IREP) - II
         ENDIF
         EPS(I) = EPS1(ISP)
      ENDDO
C
      IF (FOCKSP) THEN
          ! ... initialization
          DO IREP = 1, NREP
            NAOPI(IREP) = 0
            NAVPI(IREP) = 0
          ENDDO
!        ENDIF

         DO ISP = 1, NSP
!           ... extract the symmetry of the Abelian subgroup, IREP
            IREP = IREPSP(ISP,2)
            II = IREPSP(ISP,4)
            IF (II.GT.0) THEN
               I = IO(IREP) + II
            ELSE
               I = IO(NREP+1) + IV(IREP) - II
            ENDIF
!           ... place orbital energy into common array; I - particular spinor
            EPS(I) = EPS1(ISP)
C           Always fill this array, used in intermediate Hamiltonian formalism
            IF (I.GT.MAXSPIM)
     &          CALL QUIT ("Increase MAXSPIM in relccsd/ihm.inc")
!        ... 0/1 - intermediate Ham. space
            IPIORB(I)=0
            DO J = 1, 2
               IF(EPS(I).LE.ER_IH(2,J).AND.   ! EPS(I) <= EHMAX,EPMAX
     &            EPS(I).GE.ER_IH(1,J)) THEN  ! EPS(I) >= EHMIN,EPMIN
                 IPIORB(I)=1
               ENDIF
            ENDDO

            IF (IPIORB(I).EQ.1) THEN

             IF ( EPS(I).GE.ER_IH(1,1).AND.EPS(I).LE.ER_IH(2,1) ) THEN
               IF ( NAOPI(IREP).LT.NAO(IREP) ) THEN
                 NAOPI(IREP) = NAOPI(IREP) + 1
               ELSE
               ENDIF
             ENDIF

             IF ( EPS(I).GE.ER_IH(1,2).AND.EPS(I).LE.ER_IH(2,2) ) THEN
               IF (NAVPI(IREP) .LT. NAV(IREP) ) THEN
                 NAVPI(IREP) = NAVPI(IREP) + 1
               ELSE
               ENDIF
             ENDIF
           ENDIF
       ENDDO

       IF (REPORT.AND.DOIH) THEN

         IF (NREP.GT.8) THEN
CMI ..... write only active irreps
            NACT_REP = 0
            DO I = 1, NREP
               IF (NIO(I)+NAO(I)+NAV(I)+NIV(I).GT.0) THEN
                  ! ... move nonzero values to new arrays ...
                  NACT_REP = NACT_REP + 1
                  NAOPI_ACT(NACT_REP) = NAOPI(I)
                  NAVPI_ACT(NACT_REP) = NAVPI(I)
               ENDIF
            ENDDO

            NREP1 = 1
            NREP2 = MIN(NACT_REP,8)
  12      WRITE (IW,1005) (REPNA_ACT(I),I=NREP1,NREP2)
          WRITE (IW,1002) 'inactive occupied',(NIO_ACT(I),I=NREP1,NREP2)
          WRITE (IW,1002) 'active occupied',(NAO_ACT(I),I=NREP1,NREP2)
          WRITE (IW,1002) 'active Pi occup',(NAOPI_ACT(I),I=NREP1,NREP2)
          WRITE (IW,1002) 'active virtual ',(NAV_ACT(I),I=NREP1,NREP2)
          WRITE (IW,1002) 'active Pi virt ',(NAVPI_ACT(I),I=NREP1,NREP2)
          WRITE (IW,1002) 'inactive virtual',(NIV_ACT(I),I=NREP1,NREP2)
            NREP1 = NREP1 + 8
            NREP2 = MIN(NACT_REP,NREP2+8)
            IF (NREP1.LE.NACT_REP) GOTO 12
         ELSE
            WRITE (IW,1005) (REPNA(I),I=1,NREP)
            WRITE (IW,1002) 'inactive occupied',(NIO(I),I=1,NREP)
            WRITE (IW,1002) 'active occupied  ',(NAO(I),I=1,NREP)
            WRITE (IW,1002) 'Pi active occup  ',(NAOPI(I),I=1,NREP)
            WRITE (IW,1002) 'active virtual   ',(NAV(I),I=1,NREP)
            WRITE (IW,1002) 'Pi active virtual',(NAVPI(I),I=1,NREP)
            WRITE (IW,1002) 'inactive virtual ',(NIV(I),I=1,NREP)
         ENDIF

       ENDIF

      ENDIF

!     initialize module needed to reorder matrix elements read from MRCONEE or MDPROP
      call init_index (focksp,nsp,irepsp,irepspi)

C
 1000 FORMAT (/' Configuration in highest pointgroup'//T40,16(1X,A4))
 1001 FORMAT (/' Configuration in abelian subgroup'//T40,8(1X,A4))
 1002 FORMAT (' Spinor class : ',A,T40,16(1X,I4))
 1003 FORMAT (/' List of frozen spinors (amplitudes zeroed out)')
 1004 FORMAT (//' List of inactive and active spinors (Fock space)')
 1005 FORMAT (
     &  //' Configuration in abelian subgroup with the IH-FSCC Pi space'
     &       //T40,8(1X,A4))
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCSETI (ICALC,MXCORR,EPS,NSP,IREPSP,IREPSPI)
C
      use memory_allocator
      use allocator_parameters, only : klongint, kreal

      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Driver to generate the sorted integral files
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Produce desired integral files or rebuild Fock matrix
C          - OOOO, VOOO, VVOO and VOVO files
C          - VOVV file
C          - VVVV file
C          - Fock matrix
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "diis.inc"
#include "inpt.inc"
#include "symm.inc"
#include "eqns.inc"
#include "waio.h"
#include "complex.inc"
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      INTEGER IREPSP(NSP,4),IREPSPI(NSP,MXREP,2),ICALC,NSP
C
C---------------Local variables--------------------------------------
C
      LOGICAL ALLOC_ME,SRTALL
      INTEGER*8 NVT
      integer minbuf,nbuf,nvf,nvtmin
Caspg, replacing memget
      integer,       allocatable, target :: a_kr(:),
     &                              a_indk(:),   a_indl(:),
     &                              a_ibuf(:)
      real (kind=kreal), allocatable, target :: a_cbuf(:),
     &                              a_buf1(:),
     &                              a_buf2(:),   a_foo(:),
     &                              a_fvv(:),    a_v1(:),
     &                              a_v2(:),     a_v2a(:),
     &                              a_v3a(:),    a_vt(:),
     &                              a_fvo(:),    a_v1a(:),
     &                              a_v3(:),     a_vta(:)
      real (kind=kreal), pointer ::     a_v4(:),     a_v4a(:)
c
      integer(kind=8)        :: mxcorr
      integer(kind=8)        :: current_mem_use = 0
      integer(kind=8)        :: start_mem_use   = 0
      integer         :: freemem = 0
C
C
C---------------Executable code--------------------------------------
C
      ALLOC_ME = ICALC.GE.2
C
      call allocator_get_words_inuse(start_mem_use)
C
C  Allocate core for integral sorting step
C
      call alloc(a_kr,2*NSP+2, id="kr")
      call alloc(a_INDK,NSP*NSP    , id="indk")
      call alloc(a_INDL,NSP*NSP    , id="indl")
      call alloc(a_CBUF,NSP*NSP*RCW, id="cbuf")
      IF (EQNS.EQ.'FOCKSP') THEN
         call alloc(a_IBUF,NSP*NSP*7*256, id="ibuf")
      ELSE
         call alloc(a_IBUF,NSP*NSP*7*16, id="ibuf" )
      ENDIF
      call alloc(a_BUF1,  NSP*NSP , id="buf1" )
      call alloc(a_BUF2,  NSP*NSP , id="buf2" )
      call alloc(a_FOO,   NFOO*RCW, id="foo"  )
      call alloc(a_FVO,   NFVO*RCW, id="fvo"  )
      call alloc(a_FVV,   NFVV*RCW, id="fvv"  )
      call alloc(a_V1 ,   NV1*RCW , id="v1"   )
      call alloc(a_V1A,   NV1*RCW , id="v1a"  )
      call alloc(a_V2 ,   NV2*RCW , id="v2"   )
      call alloc(a_V2A,   NV2*RCW , id="v2a"  )
      call alloc(a_V3 ,   NV3*RCW , id="v3"   )
      call alloc(a_V3A,   NV3*RCW , id="v3a"  )

C
C     After sorting and antisymmetrizing the OOOO,VOOO and VVOO integrals
C     parts of the above arrays are not needed anymore and could be used
C     for buffering the VOVV and VVVV integrals as well. but this is not
C     a critical shortage of memory so far.
C
C     The maximum necessary size is
C
      NVT = MAX(NV5,NV6)
      NVT = MAX(NVT,INT(NV4,8))
      NVT = MAX(NVT,INT(JOOOO(2),8))
      NVT = MAX(NVT,INT(JVOOO(2),8))
C
C     The remaining space is
      call allocator_get_maxbuff(freemem,kind(a_VT))
      NVF = freemem / (2 * RCW)
C
C  Determine number of available buffers with respect to free memory
C  Records are defined for integers, calculate how many real or
C  complex words we can store on these in MINBUF
C  -----------------------------------------------------------------
C
      MINBUF = NWORDS_ON_RECORD / (RCW * IRW)
      IF (NVT.LT.NVF) THEN
         NBUF = NVT / MINBUF + 1
      ELSE
         NBUF = NVF / MINBUF
      ENDIF
      IF (NBUF.LT.1) NBUF = 1
      NVT = NBUF * MINBUF
C
C     We must be able to store at least a NV4 size array, check that.
C
      NVTMIN = MAX(NVT,INT(NV4,8))
      NVTMIN = MAX(NVTMIN,JOOOO(2))
      NVTMIN = MAX(NVTMIN,JVOOO(2))
C
      call alloc(a_VT ,NVTMIN*RCW, id="vt" )
      call alloc(a_VTA,NVTMIN*RCW, id="vta" )
      a_V4  => a_VT
      a_V4A => a_VTA

      call allocator_get_words_inuse(current_mem_use)
      mxcorr= current_mem_use - start_mem_use
C
      IF (ICALC.EQ.1) THEN
         WRITE (IW,1001) "reading and sorting integrals",MXCORR
         GOTO 450
      ENDIF
C
C  Read integrals & make fock matrices
C
      SRTALL = .TRUE.
C
C  MP2 calculations need less integrals
C  MP: DHF-ADC also needs all integral classes !
C  MP: The eqns variable hereby turns on full sorting in this case.
C
      IF (EQNS.NE.'FOCKSP') SRTALL = DOCCSD.OR.DOMP2GOLD
      IF (EQNS.EQ.'RELADC') SRTALL = .TRUE.
C
      CALL RDINTS(IPRNT,NSP,a_KR,a_INDK,a_INDL,
     &            a_CBUF,a_IBUF,NVT,a_BUF1,a_BUF2,
     &            a_V1,a_V1A,
     &            a_V2,a_V2A,a_V3,a_V3A,a_V4,a_V4A,
     &            a_VT,a_VTA,a_VT,a_VTA,EPS,a_FOO,a_FVO,
     &            a_FVV,BREIT,SRTALL,INTERFACE)

 450  CONTINUE
      call dealloc(a_kr, id="kr" )
      call dealloc(a_INDK, id="indk" )
      call dealloc(a_INDL, id="indl" )
      call dealloc(a_CBUF, id="cbuf" )
      call dealloc(a_IBUF, id="ibuf" )
      call dealloc(a_BUF1, id="buf1" )
      call dealloc(a_BUF2, id="buf2" )
      call dealloc(a_FOO , id="foo" )
      call dealloc(a_FVO , id="fvo" )
      call dealloc(a_FVV , id="fvv" )
      call dealloc(a_V1  , id="v1" )
      call dealloc(a_V1A , id="v1a" )
      call dealloc(a_V2  , id="v2" )
      call dealloc(a_V2A , id="v2a" )
      call dealloc(a_V3  , id="v3" )
      call dealloc(a_V3A , id="v3a" )

      a_V4  => NULL()
      a_V4A => NULL()
      call dealloc(a_VT  , id="vt" )
      call dealloc(a_VTA , id="vta" )
C
      RETURN
1001  FORMAT (" Memory for ",A," :",T50,I15," 8-byte words")
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCENER (ICALC,MXCORR,EPS)
      use memory_allocator
      use allocator_parameters, only : klongint, kreal
      use xmlout
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Driver for energy calculations
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Calculate energy and wave function at the
C          - SCF level
C          - MP2 level
C          - CCSD level
C          - CCSD(T) level
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
      integer icalc
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "diis.inc"
#include "inpt.inc"
#include "symm.inc"
#include "eqns.inc"
#include "results.inc"
#include "complex.inc"
#include "ccpar.inc" ! not needed, because there is no parallelization in this subroutine
C
C---------------Local variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 TOL
      LOGICAL ALLOC_ME
      INTEGER*8 NV5IRP,NV6IRP,NVIRP
C
      real(kind=kreal), allocatable, target :: a_FVO(:), a_FOO(:),
     &                                     a_FVV(:),
     &                                     a_T1(:),  a_T2(:),
     &                                     a_S1(:),  a_S2(:)
c
      real(kind=kreal), allocatable, target :: a_AZ(:), a_H(:), a_HO(:),
     &           a_HV(:), a_GO(:), a_GV(:), a_HOV(:),a_CCC(:), a_BB1(:),
     &                             a_BB2(:), a_BUF1(:), a_BUF2(:)
c
      real(kind=kreal), allocatable, target :: a_VOOO(:),a_T2S(:),
     &                             a_VVOO(:),a_VVVO(:), a_BUF3(:),
     &                             a_TAU3(:), a_T2T(:)
c
      integer(kind=8)        :: current_mem_use = 0
      integer(kind=8)        :: start_mem_use   = 0
      integer(kind=8)        :: mxcor   = 0
      integer(kind=8)        :: mxcr    = 0
      integer(kind=8)        :: mxcorr
      integer                :: freemem = 0
      integer                :: maxit_save
      logical                :: is_ccsd_converged = .false.
      integer                :: dum,irp,nbuf1,nbuf2,nbuf3,nvovv
C
C---------------Executable code--------------------------------------
C
      ALLOC_ME = ICALC.GE.2
      IF (ICALC.LE.0) RETURN

      call allocator_get_words_inuse(start_mem_use)
C
C Allocate space for the calculation of the amplitudes
C
      call alloc( a_FVO, NFVO*RCW   , id="fvo" )
      call alloc( a_FOO, NFOO*RCW   , id="foo" )
      call alloc( a_FVV, NFVV*RCW   , id="fvv" )
      call alloc( a_T1 , NDIMT1*RCW , id="t1" )
      call alloc( a_T2 , NDIMT2*RCW , id="t2" )
      call alloc( a_S1 , NDIMT1*RCW , id="s1" )
      call alloc( a_S2 , NDIMT2*RCW , id="s2" )

      call allocator_get_words_inuse(current_mem_use)
      MXCOR = current_mem_use - start_mem_use
C
C This is all we need to calculate the MP2 energy
C
      IF (DOMP2) MXCORR = MXCOR
C
      IF (DOCCSD) THEN
         MAXDIM = MXDIMCC
         NBUF1 = MAX(JOOOO(NREP+1),JVOOO(NREP+1),
     &               JVVOO(NREP+1),JVOVO(NREP+1))
         NBUF2 = NBUF1

         call alloc( a_AZ   ,    NV1*RCW              , id="az" )
         call alloc( a_H    ,    NV4*RCW              , id="h" )
         call alloc( a_HO   ,    NFOO*RCW             , id="ho" )
         call alloc( a_HV   ,    NFVV*RCW             , id="hv" )
         call alloc( a_GO   ,    NFOO*RCW             , id="go" )
         call alloc( a_GV   ,    NFVV*RCW             , id="gv" )
         call alloc( a_HOV  ,   NFVO*RCW              , id="hov" )
         call alloc( a_BUF1 ,  NBUF1*RCW              , id="buf1" )
         call alloc( a_BUF2 ,  NBUF2*RCW              , id="buf2" )

C
C        Calculate size of the 3rd work array
C        The last array is used to buffer in the <VV||VV> and <VO||VV>
C        integrals, start by making it as large as possible.
C
         call allocator_get_maxbuff(freemem,kind(a_BUF3))
         NBUF3 = freemem / RCW
C
C        Calculate maximum size necessary for the VOVV and VVVV buffers
C
         NV5IRP = 1
         NV6IRP = 1
         DO IRP = 1, NREP
            NV5IRP = MAX(NV5IRP,INT(NVO (IRP),8)*INT(NVVT(IRP),8))
            NV6IRP = MAX(NV6IRP,INT(NVVT(IRP),8)*INT(NVVT(IRP),8))
         ENDDO
C
         NVIRP = MAX(NV5IRP,NV6IRP)
         NBUF3 = MIN(NVIRP,INT(NBUF3,8))
C
C        Check the minimum size for the general use of this array
C
         NBUF3 = MAX(NBUF1,NBUF3)
C
C        Allocate the work array
C
         call alloc(a_BUF3,NBUF3*RCW, id="buf3" )

         call allocator_get_words_inuse(current_mem_use)
         MXCR = current_mem_use - start_mem_use

         mxcorr = MAX(MXCOR,mxcr)
C
      ENDIF
C
      IF (ICALC.EQ.1) THEN
         WRITE (IW,1001) "calculating amplitudes",mxcorr
         GOTO 420
      ENDIF
C
      WRITE (IW,1007) DOMP2,DOCCSD,DOCCSDT
      EQNS = 'AMPLTD'
C
C     Read the Fock matrix from file
C
      CALL FMTOFILE (.FALSE.,a_FVO,a_FOO,a_FVV)

C     Report Hartree-Fock SCF energy that is used (added to correlation energy)
      if (doxml) then
          call xml_quantity('SCF energy',ESCF,'Hartree')
      end if
C
C  Solve relativistic MP2 equations
C
      IF (DOMP2) THEN
! aspg changed
         IF (IMSTAT(5).EQ.3) THEN
            WRITE (IW,1010)
            GOTO 10
         ENDIF
C
         CALL amplitude_equation_mp2 (EPS,a_FVO,a_T1,a_S1,a_T2,a_S2)
C
C        This is a restart point: save the amplitudes and update status.
C
         CALL PUTAMPT(a_T1,a_T2)
         IMSTAT(5) = 3
         CALL CCDUMP
   10    CONTINUE
         if (doxml) then
            call xml_quantity('MP2 energy',ESCF+EMP2,'Hartree')
         end if
      ENDIF
C
C  Solve relativistic CCSD equations
C
   11 IF (DOCCSD) THEN
! aspg changed
         IF (IMSTAT(6).EQ.3) THEN
            WRITE (IW,1011)
            GOTO 20
         ENDIF
         WRITE (IW,1008) MXITCC,MXDIMCC,10.0D0**(-NTOLCC)
         MAXDIM = MXDIMCC
         MAXIT  = MXITCC


         TOL = 10.0D0**(-NTOLCC)
C
C        Initialize T1 and T2 amplitudes
C
! aspg changed
         IF (IMSTAT(6).EQ.2) THEN
C
C           Use restart amplitudes
C
            CALL GETAMPT (a_T1,a_T2)
C
! aspg changed
         ELSEIF (IMSTAT(5).EQ.3) THEN
C
C           Use MP2 amplitudes
C
            CALL GETAMPT (a_T1,a_T2)
C
         ELSE
C
C           Use standard initialization
C
            CALL XCOPY (NDIMT1,a_FVO,1,a_T1,1)
            CALL GETVVOO (a_T2)
            CALL DENOM (EPS,a_T1,a_T2,a_T1,a_T2)
            CALL ZCORE (a_T1,a_T2)
C
         ENDIF
C
         CALL PUTAMPT(a_T1,a_T2)
! aspg changed
         IMSTAT(6) = 2
         CALL CCDUMP
C
         CALL XTIME(5,1,'CCSD equations                ')

         maxit_save = MAXIT
         CALL cceqn_driver_amplitudes(maxdim,MAXIT_SAVE,EPS,a_FOO,
     &        a_FVO,a_FVV,a_T1,a_S1,a_T2,a_S2,a_AZ,a_H,a_HO,a_HV,a_HOV,
     &        a_GO,a_GV,a_BUF1,a_BUF2,a_BUF3,
     &        TOL,NBUF2,NBUF3)
         CALL XTIME(5,2,'CCSD equations                ')
         if (maxit_save.le.maxit) is_ccsd_converged = .true.

         if (doxml) then
            call xml_quantity('CCSD energy',ESCF+ECCSD,'Hartree')
         end if
C
C        This is a restart point: save the amplitudes and update status.
C
         CALL PUTAMPT(a_T1,a_T2)
! aspg changed

         if (is_ccsd_converged.or.restart_unconverged) then
            IMSTAT(6) = 3
         endif
         CALL CCDUMP
   20    CONTINUE
      ENDIF
C
C  Calculate T3 corrections
C
C
 420  CONTINUE
      call dealloc(a_AZ  )
      call dealloc(a_H   )
      call dealloc(a_HO  )
      call dealloc(a_HV  )
      call dealloc(a_HOV )
      call dealloc(a_GO )
      call dealloc(a_GV )
      call dealloc(a_CCC )
      call dealloc(a_BB1 )
      call dealloc(a_BB2 )
      call dealloc(a_BUF1)
      call dealloc(a_BUF2)
      call dealloc(a_BUF3)

      IF (DOCCSDT) THEN
! aspg changed
         IF (IMSTAT(7).EQ.3) THEN
            WRITE (IW,1012)
            GOTO 30
         ENDIF
C
C        now we start reallocating
C
         call alloc( a_VOOO , NV2*RCW            , id="vooo" )
         call alloc( a_T2S  , IVVOOT(NREP+1)*RCW , id="t2s" )
         call alloc( a_T2T  , LOVVOT(NREP+1)*RCW , id="t2t" )
         call alloc( a_VVOO , NDIMT2*RCW         , id="vvoo" )
C
C        Convert the integer*8 variable into an integer*4
C
         NVOVV = KT3VVVOT(NREP+1)
         call alloc( a_VVVO , NVOVV*RCW          , id="vvvo" )
         call alloc( a_BUF1 , IVVVTMX*RCW        , id="buf1" )
         call alloc( a_BUF2 , IVVVTMX*RCW        , id="buf2" )
         call alloc( a_BUF3 , IVVVTMX*RCW        , id="buf3" )
         call alloc( a_TAU3 , NDIMT2*RCW         , id="tau3" )

C
         IF (ICALC.EQ.1) THEN
            call allocator_get_words_inuse(current_mem_use)
            mxcr = current_mem_use - start_mem_use

            WRITE (IW,1001) "in core evaluation of triples",MXCR
            GOTO 422
         ENDIF
C
C        Check if CCSD amplitudes are available
C
! aspg changed
         IF (IMSTAT(6).EQ.3) THEN
            CALL GETAMPT (a_T1,a_T2)
! aspg changed
         ELSEIF (IMSTAT(5).EQ.3) THEN
            WRITE (IW,1006)
            CALL GETAMPT (a_T1,a_T2)
         ELSE
            DOCCSD = .TRUE.
            GOTO 11
         ENDIF
C
         CALL XTIME(5,1,'CCSD(T) evaluation            ')
         CALL T3CORR (a_FVO,EPS,a_T1,a_T2,a_T2S,a_T2T,
     &                a_BUF1,a_BUF2,a_BUF3,a_TAU3,
     &                a_VVVO,a_VVOO,a_VOOO)
         CALL XTIME(5,2,'CCSD(T) evaluation            ')

         if (doxml) then
       call xml_quantity('CCSD+T energy',ESCF+ECCSD+ET1,'Hartree')
       call xml_quantity('CCSD(T) energy',ESCF+ECCSD+ET1+ET2,'Hartree')
       call xml_quantity('CCSD-T energy',ESCF+ECCSD+ET1+ET3,'Hartree')
         end if
C
C     This is a restart point : update status.
C
! aspg changed
         IF (IMSTAT(6).EQ.3) IMSTAT(7) = 3
         CALL CCDUMP
   30    CONTINUE

 422     CONTINUE
         call dealloc(a_T2S  , id="t2s")
         call dealloc(a_T2T  , id="t2t")
         call dealloc(a_VOOO , id="vooo")
         call dealloc(a_VVOO , id="vvoo")
         call dealloc(a_VVVO , id="vvvo")
         call dealloc(a_BUF1 , id="buf1")
         call dealloc(a_BUF2 , id="buf2")
         call dealloc(a_BUF3 , id="buf3")
         call dealloc(a_TAU3 , id="tau3")
      ENDIF
C
      call dealloc(a_FVO , id="fvo")
      call dealloc(a_FOO , id="foo")
      call dealloc(a_FVV , id="fvv")
      call dealloc(a_T1  , id="t1")
      call dealloc(a_T2  , id="t2")
      call dealloc(a_S1  , id="s1")
      call dealloc(a_S2  , id="s2")
C
      RETURN
1001  FORMAT (" Core used for ",A," :",T50,I15," 8-byte words")
1007  FORMAT (//" Energy calculations"
     &        /" MP2 module active :",T44,L1
     &        /" CCSD module active :",T44,L1
     &        /" CCSD(T) module active :",T44,L1)
1008  FORMAT (//" CCSD options :"
     &        /" Maximum number of iterations :",T40,I5
     &        /" Maximum size of DIIS space :",T40,I5
     &        /" Convergence criterium :",T38,E7.1)
1006  FORMAT (/' WARNING : Calculating CCSD(T) energies with',
     &         ' MP2/unconverged amplitudes')
1010  FORMAT (/' Skipping MP2 calculation')
1011  FORMAT (/' Skipping CCSD calculation')
1012  FORMAT (/' Skipping CCSD(T) calculation')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCFOPR (ICALC,NSP,IREPSPI,MXCORR,EPS)
C
C     CCFOPR - CC First Order PRoperties, aka gradients
C
      use memory_allocator
      use allocator_parameters, only : klongint, kreal
      use interface_to_mpi
      use lambda_equation
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Driver for gradient calculations
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Calculate first order properties at the
C          - SCF level
C          - MP2 level
C          - CCSD level
C          - CCSD(T) level
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      DIMENSION IREPSPI(*)
      integer icalc,nsp,irepspi
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "files.inc"
#include "diis.inc"
#include "inpt.inc"
#include "eqns.inc"
#include "symm.inc"
#include "dgroup.h"
#include "complex.inc"
#include "ccpar.inc"
#include "dcbmp2no.h"
#include "mp1stpr.h"
C
C---------------Local variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 TOL
      COMPLEX*16 EPHASE,CEXPA
      LOGICAL ALLOC_ME
      INTEGER*8 NV5IRP,NV6IRP,NVIRP
C
C
      real(kind=kreal), allocatable, target :: a_FOO(:), a_FVO(:),
     &                             a_T1(:),  a_T2(:),   a_S1(:),
     &                             a_S2(:),  a_DOO(:),  a_DVV(:),
     &                             a_XOO(:), a_XVV(:),  a_DVO(:),
     &                             a_DOV(:), a_D2VO(:), a_D2OV(:),
     &                             a_XVO(:), a_XOV(:),  a_EVO(:),
     &                             a_EOV(:), a_EOO(:),  a_EVV(:),
     &                             a_DG(:),  a_DMO(:),  a_XMO(:),
     &                             a_BUF1(:),a_CBUF(:), a_BUF2(:),
     &                    a_BUF3(:),a_FVV(:),a_AZ(:),a_H(:),a_Hbar(:),
     &                    wovoo(:),a_HO(:), a_HV(:),a_HOV(:),
     &                    wvvvo(:), a_CCC(:),  a_BB1(:),
     &                             a_BB2(:),l1_sort(:),l2_sort(:)


      integer(kind=8)        :: current_mem_use = 0
      integer(kind=8)        :: start_mem_use   = 0
      integer(kind=8)        :: mxcor   = 0
      integer(kind=8)        :: mxcr    = 0
      integer(kind=8)        :: mxcorr
      integer                :: freemem = 0
      integer                :: work_memsize
      integer                :: eoff,eoper,iprint_no,irp,kfree
      integer                :: lfree,nbuf1,nbuf2,nbuf3,NODE,NTEST
C
C---------------Executable code--------------------------------------
C
      call allocator_get_words_inuse(start_mem_use)
      call legacy_lwork_get(work_memsize)
C
      if (doccsdg) then
          WRITE (IW,*) 'Calculate CCSD gradient'
      endif

      if (domp2g) then
         WRITE (IW,*) 'calculate MP2 gradient'
      endif
C
C     here we allocate the variables common to both schemes
C     the specific ones are allocated as needed
C
      call alloc( a_FOO   , NFOO*RCW         , id="foo" )
      call alloc( a_FVO   , NFVO*RCW         , id="fvo" )
      call alloc( a_FVV   , NFVV*RCW         , id="fvv" )
      call alloc( a_T1    , NDIMT1*RCW       , id="t1" )
      call alloc( a_T2    , NDIMT2*RCW       , id="t2" )
      call alloc( a_S1    , NDIMT1*RCW       , id="s1" )
      call alloc( a_S2    , NDIMT2*RCW       , id="s2" )
      call alloc( a_DOO   , NFOO*RCW         , id="doo" )
      call alloc( a_DVV   , NFVV*RCW         , id="dvv" )
      call alloc( a_DVO   , NDIMX*RCW        , id="dvo" )
      call alloc( a_DOV   , NDIMX*RCW        , id="dov" )
      call alloc( a_DMO   , NORBT*NORBT*NZ   , id="dmo" )
      call alloc( a_XMO   , NORBT*NORBT*NZ   , id="xmo" )
      NBUF1 = MAX(JOOOO(NREP+1),JVOOO(NREP+1),
     &            JVVOO(NREP+1),JVOVO(NREP+1))

C     Calculate the size of the array needed for the sorted T2 array that
C     is used in MOLTRA. In spinfree calculations this will be larger than
C     the dimension used in RELCCSD (was a nasty bug in DIRAC08)
C     This initialization does not fix the problem with spinfree calculations,
C     however, various routines (e.g. DOOTODIR and DVVTODIR) do need to be
C     adapted to the fact that in spinfree calculations there is a many-to-1
C     correspondence between the list of irreps used by RELCCSD and the
C     fermion irreps of DIRAC. Probably the best solution is to rewrite these
C     routines using an index vector that gives the correspondence.

      CALL INIT_T1(a_T2,.FALSE.)
      CALL INIT_T2(a_T2,.FALSE.)
      NBUF1 = MAX(NBUF1, NTMATR1, NTMATR2)
      NBUF2 = NBUF1
C
C The DIIS algorithm is used for the Z-vector and Lamda-equations
C Allocate the workarrays only once with the maximum dimensions
C
      MAXIT  = MXITGR
      MAXDIM = MXDIMGR

      call alloc( a_BUF1  , NBUF1*RCW             , id="buf1" )
      call alloc( a_BUF2  , NBUF2*RCW             , id="buf2" )
C
C     Read the Fock matrix from file
C
      if (icalc.ge.2) CALL FMTOFILE (.FALSE.,a_FVO,a_FOO,a_FVV)
C
C     Set equation type for use in MP2EQN
C
      EQNS = 'LAMBDA'
C
C     Get MP2 relaxed density matrix and compute first order properties.
C
C     First comes the old code for the MO implementation that gives
C     an approximate value. This module is mainly intended for tests
C     and backwards compatibility. The deafult is now to use the new
C     AO based module that gives full relaxation and the exact result.
C
      IF (DOMP2GOLD) THEN
C
C allocating stuff specific for this option
C
      call alloc( a_XOO   , NFOO*RCW         , id="xoo" )
      call alloc( a_XVV   , NFVV*RCW         , id="xvv" )
      call alloc( a_D2VO  , NDIMX*RCW        , id="d2vo" )
      call alloc( a_D2OV  , NDIMX*RCW        , id="d2ov" )
      call alloc( a_XVO   , NDIMX*RCW        , id="xvo" )
      call alloc( a_XOV   , NDIMX*RCW        , id="xov" )
      call alloc( a_EVO   , NEOPER*NDIMX*RCW , id="evo" )
      call alloc( a_EOV   , NEOPER*NDIMX*RCW , id="eov" )
      call alloc( a_EOO   , NFOO*RCW         , id="eoo" )
      call alloc( a_EVV   , NFVV*RCW         , id="evv" )
      call alloc( a_DG    , NDIMX*RCW        , id="dg" )

      call alloc( a_CCC   , MAXDIM                , id="ccc" )
      call alloc( a_BB1   , (MAXDIM+1)*(MAXDIM+2) , id="bb1" )
      call alloc( a_BB2   , (MAXDIM+1)*(MAXDIM+2) , id="bb2" )
      call alloc( a_CBUF  , 2*NSP*NSP             , id="cbuf" )
C
C     Calculate size of the 3rd work array that is passed on to the
C     response module as work, we make it as large as possible.
C
      call allocator_get_maxbuff(freemem,kind(a_BUF3))
C aspg: the subtraction of work_memsize here is due to the
C       fact that later on in the call to dirnod2 there will
C       be allocation(s) of the same size for the integral
C       drivers etc; so this is in a way a hack, because of
C       limitations elsewhere. 20120720
!     NBUF3 = (freemem - work_memsize) / RCW
! hjaaj Aug. 2015: modified because this could give NBUF3 < 0
      IF (MYPROC .NE. MASTER) THEN
         NBUF3 = 0
      ELSE
         NBUF3 = freemem / RCW
      END IF

C
C     Allocate now also this remaining work array
C
      call alloc( a_BUF3  , NBUF3*RCW             , id="buf3" )
C
C  now our reference work array is a_BUF3 and KFREE an LFREE relate to that...
C
      KFREE = 1
      LFREE = NBUF3 *RCW

      call allocator_get_words_inuse(current_mem_use)
      mxcorr = current_mem_use - start_mem_use

      IF (ICALC.EQ.1) THEN
         WRITE (IW,1001) "Gradient Calculation",MXCORR
         GOTO 430
      ENDIF

         WRITE (IW,1007) MXITGR,MXDIMGR,10.0D0**(-NTOLGR)
C
C        Start with SCF (reference determinant) expectation values
C
         WRITE (IW,1002) '    SCF'
         CALL XCOPY (NFOO,A0,0,a_DOO,1)
         CALL XCOPY (NFVV,A0,0,a_DVV,1)
         CALL XCOPY (NDIMX,A0,0,a_DVO,1)
         CALL XCOPY (NDIMX,A0,0,a_DOV,1)
C
         DO EOPER = 1, NEOPER
            EOFF = (EOPER-1)*NDIMX*RCW+1
            CALL RDPROP(IPRNT,NAMEE(EOPER),NSP,LTR,a_CBUF,IREPSPI,
     &                  EPHASE,a_EVO,a_EOV,a_EOO,a_EVV)
            CALL EXPVAL(a_DOO,a_DVV,a_DVO,a_DOV,
     &                  a_EOO,a_EVV,a_EVO,a_EOV,
     &                  EPHASE,CEXPA)
            WRITE (IW,1004) NAMEE(EOPER),CEXPA
         ENDDO
C
C        Solve relativistic MP2 lambda equations
C
         CALL lambda_equation_mp2 (EPS,a_FVO,a_T1,a_S1,a_T2,a_S2)
C
C        Get amplitude part of density matrix
C
         CALL CCDENSC (1,a_T1,  a_S1,  a_T2,  a_S2,
     &                   a_BUF1,a_BUF2,a_BUF3,a_DOO,a_DVV,
     &                   a_DVO, a_DOV, a_DMO)
C        Solve Z-vector equations
C
         IF (.NOT. NOZG) THEN
            TOL = 10.0D0**(-NTOLGR)
            MAXIT  = MXITGR
            MAXDIM = MXDIMGR
            CALL CCDENSZ (1, NSP, IREPSPI, EPS,  a_T1, a_S1, a_T2, a_S2,
     &                    a_BUF1, a_BUF2, a_BUF3, NBUF2 , NBUF3, a_CBUF,
     &                    a_XVO , a_XOV , a_DOO , a_DVV , a_DVO, a_DOV,
     &                    a_EVO , a_EOV , a_D2VO, a_D2OV, a_DG,
     &                    a_CCC , a_BB1 , a_BB2 , TOL)
         ENDIF
C
C        Calculate expecation value
C
         IF (NOZG) THEN
            WRITE (IW,1002) 'MP2 noZ'
         ELSE
            WRITE (IW,1002) '    MP2'
         END IF
         DO EOPER = 1, NEOPER
            EOFF = (EOPER-1)*NDIMX*RCW+1
            CALL RDPROP(IPRNT,NAMEE(EOPER),NSP,LTR,a_CBUF,IREPSPI,
     &                  EPHASE,a_EVO,a_EOV,a_EOO,a_EVV)
            CALL EXPVAL(a_DOO,a_DVV,a_DVO,a_DOV,
     &                  a_EOO,a_EVV,a_EVO,a_EOV,
     &                  EPHASE,CEXPA)
            WRITE (IW,1004) NAMEE(EOPER),CEXPA
         ENDDO

 430     CONTINUE
         call dealloc(a_XOO, id="xoo" )
         call dealloc(a_XVV, id="xvv" )
         call dealloc(a_D2VO, id="d2vo" )
         call dealloc(a_D2OV, id="d2ov" )
         call dealloc(a_XVO, id="xvo" )
         call dealloc(a_XOV, id="xov" )
         call dealloc(a_EVO, id="evo" )
         call dealloc(a_EOV, id="eov" )
         call dealloc(a_EOO, id="eoo" )
         call dealloc(a_EVV, id="evv" )
         call dealloc(a_DG, id="dg" )
         call dealloc(a_CCC, id="ccc" )
         call dealloc(a_BB1, id="bb1" )
         call dealloc(a_BB2, id="bb2" )
         call dealloc(a_CBUF, id="cbuf" )
         call dealloc(a_BUF3, id="buf3" )
C
      ENDIF ! End of MO code
C
C     New AO implementation
C
      IF (DOMP2G) THEN
C
C
C     Calculate size of the 3rd work array that is passed on to the
C     response module as work, we make it as large as possible.
C
      call allocator_get_maxbuff(freemem,kind(a_BUF3))
C aspg: the subtraction of work_memsize here is due to the
C       fact that later on in the call to dirnod2 there will
C       be allocation(s) of the same size for the integral
C       drivers etc; so this is in a way a hack, because of
C       limitations elsewhere. 20120720
!     NBUF3 = (freemem - work_memsize) / RCW
! hjaaj Aug. 2015: modified because this could give NBUF3 < 0
      IF (MYPROC .NE. MASTER) THEN
         NBUF3 = 0
      ELSE
         NBUF3 = freemem / RCW
      END IF
C
C     Allocate now also this remaining work array
C
      call alloc( a_BUF3  , NBUF3*RCW             , id="buf3" )
C
      call allocator_get_words_inuse(current_mem_use)
      mxcorr = current_mem_use - start_mem_use

      IF (ICALC.EQ.1) THEN
         WRITE (IW,1001) "Gradient Calculation",MXCORR
         GOTO 431
      ENDIF

          IF (NOZG) THEN
             WRITE (IW,1008)
          ELSE
             WRITE (IW,1018)
          ENDIF
C
C         Slaves will only help in the integral calculation.
C
C         this means that, for the time being, we can deallocate everything?
#if defined (VAR_MPI)
          IF (MYPROC.NE.MASTER) THEN
             IF (ICALC.EQ.2) CALL DIRNOD2()
             goto 431
          ENDIF
#endif
C
C        Get CCSD amplitudes and lambdas
C
         CALL lambda_equation_mp2 (EPS, a_FVO, a_T1, a_S1, a_T2, a_S2)
C
C        Get amplitude part of density matrix
C
         CALL CCDENSC (1, a_T1  , a_S1  , a_T2,   a_S2,
     &                    a_BUF1, a_BUF2, a_BUF3, a_DOO,
     &                    a_DVV , a_DVO , a_DOV , a_DMO)
C
C        Solve Z-vector equations
C
         IF (.NOT.NOZG)
     &      CALL MP2ZVEC( a_T2, a_BUF1, a_BUF2, a_BUF3,
     &                   NBUF3, a_XMO , a_DMO)
C
C        Create upper triangle of density matrix
C
         CALL QHMRST( a_DMO, NORBT,NZ,NORBT,NORBT)
C        Write density matrix to file for later processing
C
         CALL STORE_CC_DENSITY ('MP2 ', a_DMO )
C
C        Generate natural orbitals, if desired
C
         IF( DONATORB )THEN
C
C         MP2 - nat. orb. and nat. orb. occ. num.
C
          IPRINT_NO = 1
          call dealloc(a_BUF3, id="buf3")
C
          CALL MP2_NATORB( a_DMO,IPRINT_NO,NOZG)
         END IF


C
C        Release the slaves (in case of parallel runs)
C
#if defined (VAR_MPI)
         NTEST = -1
         DO NODE = 1,NMPROC-1
            call interface_mpi_SEND(NTEST,1,NODE,27,
     &                     global_communicator)
         ENDDO
#endif
 431  CONTINUE
         if(allocated(a_BUF3)) call dealloc(a_BUF3, id="buf3")
      ENDIF


!-----------------------------------------------------------
**ashee
** new section for lambda equation in coupled cluster level.
!-----------------------------------------------------------
      IF (DOCCSDG) THEN
C
C     Calculate size of the 3rd work array that is passed on to the
C     response module as work, we make it as large as possible.
C
      call allocator_get_maxbuff(freemem,kind(a_BUF3))
C aspg: the subtraction of work_memsize here is due to the
C       fact that later on in the call to dirnod2 there will
C       be allocation(s) of the same size for the integral
C       drivers etc; so this is in a way a hack, because of
C       limitations elsewhere. 20120720
       NBUF3 = freemem / RCW

         NV5IRP = MAX(INT(JVOVV(NREP+1),8),INT(JOVVV(NREP+1),8))
         NV6IRP = 1
         DO IRP = 1, NREP
            NV6IRP = MAX(NV6IRP,INT(NVVT(IRP),8)*INT(NVVT(IRP),8))
         ENDDO

C        Check the minimum size for the general use of this array

         NVIRP = MAX(NV5IRP,NV6IRP)

         NBUF3 = MIN(NVIRP,INT(NBUF3,8))

         IF (NV5IRP.GT.NBUF3) THEN

           write(iw,*) 'The calculation may Fail because of very high //
     &                  memory requirement. Increase mw'
           NBUF3 = NV5IRP
         
         ENDIF

         NBUF3 = MAX(NBUF1,NBUF3)

C     Allocate now also this remaining work array
        call alloc( a_BUF3  , NBUF3*RCW             , id="buf3" )


         MAXDIM = MXDIMCC
         MAXIT  = MXITCC

! arrays specific for the use in the call sequence of cceqns subroutine.

         call alloc( a_AZ   ,    NV1*RCW              , id="az" )
         call alloc( a_H    , ivoov(nrep+1)*RCW              , id="h" )
         call alloc( a_Hbar , iovvo(nrep+1)*RCW           , id="hbar" )
         call alloc( wovoo  , iovoot(nrep+1)*RCW   , id="w1v" )
!        call alloc( wvvvo  ,    NV5*RCW            , id="w3v" )
         allocate(wvvvo(nv5*RCW)) !consumes too much memory
         call alloc( a_HO   ,   NFOO*RCW             , id="ho" )
         call alloc( a_HV   ,   NFVV*RCW             , id="hv" )
         call alloc( a_HOV  ,   NFVO*RCW              , id="hov" )
         call alloc( a_CCC  ,   MAXDIM                , id="ccc" )
         call alloc( a_BB1  ,   (MAXDIM+1)*(MAXDIM+2) , id="bb1" )
         call alloc( a_BB2  ,   (MAXDIM+1)*(MAXDIM+2) , id="bb2" )



      call allocator_get_words_inuse(current_mem_use)
      mxcorr = current_mem_use - start_mem_use

      IF (ICALC.EQ.1) THEN
         WRITE (IW,1001) "Gradient Calculation",MXCORR
         GOTO 451
      ENDIF

          WRITE (IW,1010)
C
C         Slaves will only help in the integral calculation.
C
C         this means that, for the time being, we can deallocate everything?
!#if defined (VAR_MPI)
!         IF (MYPROC.NE.MASTER) THEN
!            IF (ICALC.EQ.2) CALL DIRNOD2()
!         ENDIF
!#endif

!--------------------------------------------------
**ashee
** get coupled cluster amplitudes
!--------------------------------------------------
!      IF (IMSTAT(5).EQ.3) THEN
!        CALL GETAMPT (a_T1,a_T2)
!      ENDIF

!------------------------------------------
** solve coupled cluster lambda equation.
** ashee
!------------------------------------------

!        get coupled cluster lambdas



         CALL CCDUMP

    5   format(10f10.6)

         TOL = 10.0D0**(-NTOLCC)

         CALL cceqn_driver_lambda(MAXIT,EPS,a_FOO,a_FVO,a_FVV,a_T1,a_S1,
     &          a_T2 ,a_S2 ,a_AZ,wvvvo,a_H,a_Hbar,wovoo,a_HO,a_HV,a_HOV,
     &           a_BUF1,a_BUF2,a_BUF3,
     &           a_CCC,a_BB1,a_BB2,
     &           TOL,NBUF2,NBUF3)


       call cc_density (a_T1,a_T2,a_S1,a_S2,a_doo,a_dvv,a_buf1,a_buf2,
     &                   a_dmo)


C        Solve Z-vector equations
         IF (.NOT.NOZG)
     &      CALL MP2ZVEC( a_S2, a_BUF1, a_BUF2, a_BUF3,
     &                   NBUF3, a_XMO , a_DMO)


C        Write density matrix to file for later processing
         if (myproc==0) then 

           CALL STORE_CC_DENSITY ('CCSD', a_DMO )

         endif

 451  CONTINUE
         if(allocated(a_BUF3)) call dealloc(a_BUF3, id="buf3")

         call dealloc( a_AZ )
         call dealloc( a_H  )
         call dealloc( a_Hbar  )
         call dealloc( wovoo  )
!        call dealloc( wvvvo  )
         deallocate(wvvvo)
         call dealloc( a_HO )
         call dealloc( a_HV )
         call dealloc( a_HOV)
         call dealloc( a_CCC)
         call dealloc( a_BB1)
         call dealloc( a_BB2)

      ENDIF
C
      call dealloc( a_FOO, id="foo")
      call dealloc( a_FVO, id="fvo")
      call dealloc( a_FVV, id="fvv")
      call dealloc( a_T1 , id="t1")
      call dealloc( a_T2 , id="t2")
      call dealloc( a_S1 , id="s1")
      call dealloc( a_S2 , id="s2")
      call dealloc( a_DOO, id="doo")
      call dealloc( a_DVV, id="dvv")
      call dealloc( a_DVO, id="dvo")
      call dealloc( a_DOV, id="dov")
      call dealloc( a_DMO, id="dmo")
      call dealloc( a_XMO, id="xmo")
      call dealloc(a_BUF1, id="buf1")
      call dealloc(a_BUF2, id="buf2")

C
      RETURN
1001  FORMAT (" Core used for ",A," :",T50,I15," 8-byte words")
1007  FORMAT (//" Gradient calculations with the old MP2 module"
     &        //" Z-vector equations are done in MO basis:"
     &        /" Maximum number of iterations :",T40,I5
     &        /" Maximum size of DIIS space :",T40,I5
     &        /" Convergence criterium :",T38,E7.1)
1010  FORMAT (//" Relaxed density matrix for CCSD properties")
1008  FORMAT (//" Unrelaxed density matrix for MP2 properties")
1018  FORMAT (//" Relaxed density matrix for MP2 properties")
1002  FORMAT(//,1X,A7,' Gradient results',//3X,'Expectation value',13X,
     & 'Real part',8X,'Imaginary part')
1004  FORMAT('  <',A8,'> ',T30,2F20.13)
1009  FORMAT (/" WARNING : Cannot solve CCSD lambda equations :",
     &         " T1 and T2 amplitudes are not available "/
     &         " Check that you specified DOENER and DOCCSD !")
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCHESS (ICALC,MXCORR,EPS,NSP,IREPSPI)
C
      use memory_allocator
      use allocator_parameters, only : klongint, kreal
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Driver for hessian calculations
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Calculate second order properties at the
C          - SCF (RPA) level
C          - MP2 level
C          - CCSD level
C          - CCSD(T) level
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      INTEGER ICALC,NSP,IREPSPI
      DIMENSION IREPSPI(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "diis.inc"
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 APHASE(MAXOP),BPHASE
      REAL*8 TOL
      LOGICAL BDONE
      LOGICAL ALLOC_ME
C
      real(kind=kreal), allocatable, target :: a_X1(:),   a_X2(:),
     &                                         a_Y1(:),   a_Y2(:),
     &                                         a_BVO(:),  a_BOV(:),
     &                                         a_DG(:),   a_CCC(:),
     &                                         a_BB1(:),  a_BB2(:),
     &                                         a_BUF1(:), a_BUF2(:),
     &                                         a_CBUF(:)

      real(kind=kreal), allocatable, target :: a_AVO(:), a_AOV(:)

      integer(kind=8)        :: current_mem_use = 0
      integer(kind=8)        :: start_mem_use   = 0
      integer(kind=8)        :: mxcor   = 0
      integer(kind=8)        :: mxcr    = 0
      integer(kind=8)        :: mxcorr
      integer                :: freemem = 0
      integer                :: aoff,aoper,boper,dum,ifreq,nbuf1
C
C---------------Executable code--------------------------------------
C
      ALLOC_ME = ICALC.GE.2
      IF (ICALC.LE.0) RETURN
C
      call allocator_get_words_inuse(start_mem_use)
C
C Allocate space for the calculation of the hessian
C
      call alloc( a_X1    ,          NDIMX*RCW    , id="x1"  )
      call alloc( a_X2    ,          NDIMX*RCW    , id="x2"  )
      call alloc( a_Y1    ,          NDIMX*RCW    , id="y1"  )
      call alloc( a_Y2    ,          NDIMX*RCW    , id="y2"  )
      call alloc( a_AVO   , NAOPER * NDIMX*RCW    , id="avo" )
      call alloc( a_AOV   , NAOPER * NDIMX*RCW    , id="aov" )
      call alloc( a_BVO   ,          NDIMX*RCW    , id="bvo" )
      call alloc( a_BOV   ,          NDIMX*RCW    , id="bov" )
      call alloc( a_DG    ,          NDIMX*RCW    , id="dg"  )
      MAXDIM = MXDIMRP
      call alloc( a_CCC   , MAXDIM                , id="ccc" )
      call alloc( a_BB1   , (MAXDIM+1)*(MAXDIM+2) , id="bb1" )
      call alloc( a_BB2   , (MAXDIM+1)*(MAXDIM+2) , id="bb2" )
      NBUF1 = MAX(JVVOO(NREP+1),JVOVO(NREP+1))
      call alloc( a_BUF1  , NBUF1*RCW             , id="buf1"  )
      call alloc( a_BUF2  , NBUF1*RCW             , id="buf2"  )
      call alloc( a_CBUF  , 2*NSP*NSP             , id="cbuf"  )

      call allocator_get_words_inuse(current_mem_use)
      mxcorr = current_mem_use - start_mem_use
C
      IF (ICALC.LE.1) THEN
         WRITE (IW,1001) "RPA calculation",mxcorr
         goto 4441
      ENDIF

      a_X1    = 0.0d0
      a_X2    = 0.0d0
      a_Y1    = 0.0d0
      a_Y2    = 0.0d0
      a_AVO   = 0.0d0
      a_AOV   = 0.0d0
      a_BVO   = 0.0d0
      a_BOV   = 0.0d0
      a_DG    = 0.0d0
      a_CCC   = 0.0d0
      a_BB1   = 0.0d0
      a_BB2   = 0.0d0
      a_BUF1  = 0.0d0
      a_BUF2  = 0.0d0
      a_CBUF  = 0.0d0

C  Solve relativistic Random Phase Approximation equations
C
      WRITE (IW,1009) MXITRP,MXDIMRP,10.0D0**(-NTOLRP)
      MAXIT  = MXITRP
      TOL = 10.0D0**(-NTOLRP)
      MAXDIM = MXDIMRP
      DO AOPER = 1, NAOPER
         AOFF = (AOPER-1)*NDIMX*RCW
         CALL RDPROP(IPRNT,NAMEA(AOPER),NSP,LFA,a_CBUF,IREPSPI,
     &               APHASE(AOPER), a_AVO(AOFF+1),a_AOV(AOFF+1),DUM,DUM)
      ENDDO
      DO BOPER = 1, NBOPER
         BDONE = .FALSE.
         CALL RDPROP(IPRNT,NAMEB(BOPER),NSP,LFA,a_CBUF,IREPSPI,
     &               BPHASE,a_BVO,a_BOV,DUM,DUM)
         DO IFREQ = 1, NFREQ
            IF (IPRNT.GE.1) WRITE (IW,1010) NAMEB(BOPER),EFREQ(IFREQ)
            CALL RPAA
     &       (MAXIT,NAOPER,NAMEA,NAMEB(BOPER),APHASE,BPHASE,BDONE,EPS,
     &        EFREQ(IFREQ), a_AVO, a_AOV, a_BVO, a_BOV, a_X1,
     &        a_X2, a_Y1, a_Y2, a_DG, a_BUF1, a_BUF2,
     &        a_CCC, a_BB1, a_BB2, TOL,IPRNT,DEBUG,
     &        TIMING)
            CALL FLSHFO(IW)
         ENDDO
      ENDDO
C
4441  continue
      call dealloc( a_X1   , id="x1"  )
      call dealloc( a_X2   , id="x2"  )
      call dealloc( a_Y1   , id="y1"  )
      call dealloc( a_Y2   , id="y2"  )
      call dealloc( a_AVO  , id="avo" )
      call dealloc( a_AOV  , id="aov" )
      call dealloc( a_BVO  , id="bvo" )
      call dealloc( a_BOV  , id="bov" )
      call dealloc( a_DG   , id="dg"  )
      call dealloc( a_CCC  , id="ccc" )
      call dealloc( a_BB1  , id="bb1" )
      call dealloc( a_BB2  , id="bb2" )
      call dealloc( a_BUF1 , id="buf1"  )
      call dealloc( a_BUF2 , id="buf2"  )
      call dealloc( a_CBUF , id="cbuf"  )

      RETURN
1001  FORMAT (" Core used for ",A," :",T50,I15," 8-byte words")
1009  FORMAT (//" Propagator (RPA) module"
     &        /" Maximum number of iterations :",T40,I5
     &        /" Maximum size of DIIS space :",T40,I5
     &        /" Convergence criterium :",T38,E7.1)
1010  FORMAT (//" Random Phase Approximation : ",
     &        /" Calculate propagators of ",A8," (B)",
     &        " at frequency :",F15.5)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCEXC (ICALC,MXCORR,EPS,NSP,IREPSPI)
C
      use memory_allocator
      use allocator_parameters, only : klongint, kreal
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Driver for excited state calculations
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Calculate excitation energies at the
C          - SCF (RPA) level
C          - MP2 level
C          - CCSD level
C          - CCSD(T) level
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      INTEGER ICALC,NSP,IREPSPI
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "diis.inc"
#include "inpt.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 APHASE(MAXOP),BPHASE
      REAL*8 TOL
      LOGICAL BDONE
      LOGICAL ALLOC_ME
C
      real(kind=kreal), allocatable, target :: a_X1(:),   a_X2(:),
     &                                         a_Y1(:),   a_Y2(:),
     &                                         a_DG(:),   a_CCC(:),
     &                                         a_BB1(:),  a_BB2(:),
     &                                         a_BB3(:),  a_BB4(:),
     &                                         a_BUF1(:), a_BUF2(:),
     &                                         a_FREQR(:),
     &                                         a_FREQI(:),
     &                                         a_FREQD(:)

      integer(kind=8)        :: current_mem_use = 0
      integer(kind=8)        :: start_mem_use   = 0
      integer(kind=8)        :: mxcor   = 0
      integer(kind=8)        :: mxcr    = 0
      integer(kind=8)        :: mxcorr
      integer                :: freemem = 0
      integer                :: irep,nbuf1
C
C---------------Executable code--------------------------------------
C
      ALLOC_ME = ICALC.GE.2
      IF (ICALC.LE.0) RETURN
C
      call allocator_get_words_inuse(start_mem_use)
C
C Allocate space for the excitation energy calculation
C
      call alloc( a_X1    , NDIMX*RCW     , id="x1"  )
      call alloc( a_X2    , NDIMX*RCW     , id="x2"  )
      call alloc( a_Y1    , NDIMX*RCW     , id="y1"  )
      call alloc( a_Y2    , NDIMX*RCW     , id="y2"  )
      call alloc( a_DG    , NDIMX*RCW     , id="dg"  )
      MAXDIM = MAX(6,2*MXDIMRP+2)
      call alloc( a_FREQR , MAXDIM        , id="freqr" )
      call alloc( a_FREQI , MAXDIM        , id="freqi" )
      call alloc( a_FREQD , MAXDIM        , id="freqd" )
      call alloc( a_CCC   , MAXDIM*MAXDIM , id="ccc" )
      call alloc( a_BB1   , MAXDIM*MAXDIM , id="bb1" )
      call alloc( a_BB1   , MAXDIM*MAXDIM , id="bb2" )
      call alloc( a_BB1   , MAXDIM*MAXDIM , id="bb3" )
      call alloc( a_BB1   , MAXDIM*MAXDIM , id="bb4" )
      NBUF1 = MAX(JVVOO(NREP+1),JVOVO(NREP+1))
      call alloc( a_BUF1  , NBUF1*RCW     , id="buf1"  )
      call alloc( a_BUF2  , NBUF1*RCW     , id="buf2"  )

      call allocator_get_words_inuse(current_mem_use)
      mxcorr = current_mem_use - start_mem_use
C
      IF (ICALC.LE.1) THEN
         WRITE (IW,1001) "excitation energies",mxcorr
         goto 4451
      ENDIF
C
C  Solve relativistic Random Phase Approximation equations
C
      WRITE (IW,1009) MXITRP,MXDIMRP,10.0D0**(-NTOLRP)
      MAXIT   = MXITRP
      MAXDIM  = MXDIMRP
      TOL = 10.0D0**(-NTOLRP)
C
C     This is very old code and..not well tested : better warn
      WRITE (IW,*)
     & " WARNING : This module is obsolete and is not tested anymore"
      WRITE (IW,*) " Use it at your own risk....."
C
C
      DO IREP = 1, NREP
        IF (NEXC(IREP).NE.0) THEN
          WRITE (IW,1011) NEXC(IREP),REPNA(IREP+NREP)
          CALL RPAB (MAXDIM,MAXIT,
     &        EPS, a_X1 , a_X2 , a_Y1 , a_Y2 , a_DG , a_BUF1 ,
     &         a_BUF2 , a_CCC , a_BB1 , a_BB2 , a_BB3 , a_BB4 ,
     &         a_FREQR , a_FREQI , a_FREQD ,
     &        NEXC(IREP),IREP,TOL,IPRNT,DEBUG,TIMING)
          CALL FLSHFO(IW)
        ENDIF
      ENDDO

4451  continue
      call dealloc( a_X1    , id="x1"  )
      call dealloc( a_X2    , id="x2"  )
      call dealloc( a_Y1    , id="y1"  )
      call dealloc( a_Y2    , id="y2"  )
      call dealloc( a_DG    , id="dg"  )
      call dealloc( a_FREQR , id="freqr" )
      call dealloc( a_FREQI , id="freqi" )
      call dealloc( a_FREQD , id="freqd" )
      call dealloc( a_CCC   , id="ccc" )
      call dealloc( a_BB1   , id="bb1" )
      call dealloc( a_BB1   , id="bb2" )
      call dealloc( a_BB1   , id="bb3" )
      call dealloc( a_BB1   , id="bb4" )
      call dealloc( a_BUF1  , id="buf1"  )
      call dealloc( a_BUF2  , id="buf2"  )
C
      RETURN
1001  FORMAT (" Core used for ",A," :",T50,I15," 8-byte words")
1009  FORMAT (//" Excitation energy module (RPA)"
     &        /" Maximum number of iterations :",T40,I5
     &        /" Maximum size of reduced space :",T40,I5
     &        /" Convergence criterium :",T38,E7.1)
1011  FORMAT (//" Random Phase Approximation : ",
     &        /" Calculate",I3," excitation energies for irrep ",A4)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCFSPC (ICALC,MXCORR,EPS1,EPS,NSP,IREPSP,IREPSPI)
C
      use memory_allocator
      use allocator_parameters, only : klongint, kreal

      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Driver for Fock space coupled cluster calculations
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Calculate Fock space coupled cluster energies at the
C          - CCSD level
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS1(*),EPS(*)
      INTEGER ICALC,NSP,IREPSP,IREPSPI
      DIMENSION IREPSP(*),IREPSPI(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "diis.inc"
#include "ihm.inc"
#include "eqns.inc"
#include "results.inc"
#include "inpt.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL ALLOC_ME,FOCKSP,REPORT
      CHARACTER*2 IFSS(6)
      REAL*8 TOL
      INTEGER*8 NV5IRP,NV6IRP,NVIRP
c
      real(kind=kreal), allocatable, target :: a_FVO (:), a_FOO (:),
     &                             a_FVV (:),
     &                             a_T1  (:), a_T2  (:), a_S1  (:),
     &                             a_S2  (:), a_AZ  (:), a_H   (:),
     &                             a_HO  (:), a_HV  (:), a_HOV (:),
     &                             a_GO(:), a_GV(:), a_CCC (:),
     &                             a_BB1 (:), a_BB2 (:), a_BUF1(:),
     &                             a_BUF2(:), a_BUF3(:)
      integer                :: freemem = 0
      integer(kind=8)        :: mxcorr, mxcorr1
      integer                :: maxit_save
      logical                :: is_ccsd_converged = .false.
      integer                :: dum,ifstat,irp,nbuf1,nbuf2,nbuf3
C
C---------------Executable code--------------------------------------
C
      ALLOC_ME = ICALC.GE.2
      IF (ICALC.LE.0) RETURN
      EQNS = 'FOCKSP'
      FOCKSP = .TRUE.
      REPORT = ICALC.GT.1
      IF (REPORT) WRITE (IW,1007)

CMI ... more print out
      IF (FOCKSP.AND.DOIH.AND.REPORT) THEN
        WRITE(IW,'(/A)') ' Intermediate Hamiltonian version !'
        WRITE(IW,'(1X,A,F20.7,A,F20.7,A)')
     &  'Energy range for the occupied Pi space is: <',
     &  ER_IH(1,1),';',ER_IH(2,1),'>'

        WRITE(IW,'(1X,A,F20.7,A,F20.7,A)')
     &  'Energy range for the virtual Pi space is:  <',
     &  ER_IH(1,2),';',ER_IH(2,2),'>'
      ENDIF
C
C Allocate space for the Fock space CCSD energy calculation
C
      MXCORR = 0
C
C  Set up symmetry tables and parallel distribution list
C  for use in Fock space calculations
C
      CALL SYMTAB (NELEC,NFROZ,NINACT,NACT,NSP,EPS1,EPS,IREPSP,IREPSPI,
     &             FOCKSP,REPORT,NELEC_F1,NELEC_F2)
C
C     Ask for the memory requirements of the modules
C
      MXCORR1 = 0
      CALL CCSETI (ICALC,MXCORR1,EPS,NSP,IREPSP,IREPSPI)
      MXCORR = MAX(MXCORR,MXCORR1)
C
C Allocate space for the calculation of the amplitudes
C
      MAXDIM = MXDIMFS
      NBUF1 = MAX(JOOOO(NREP+1),JVOOO(NREP+1),
     &             JVVOO(NREP+1),JVOVO(NREP+1))
      NBUF2 = NBUF1

      call alloc( a_FVO , NFVO*RCW              , id="fvo")
      call alloc( a_FOO , NFOO*RCW              , id="foo")
      call alloc( a_FVV , NFVV*RCW              , id="fvv")
      call alloc( a_T1  , NDIMT1*RCW            , id="t1")
      call alloc( a_T2  , NDIMT2*RCW            , id="t2")
      call alloc( a_S1  , NDIMT1*RCW            , id="s1")
      call alloc( a_S2  , NDIMT2*RCW            , id="s2")
      call alloc( a_AZ  , NV1*RCW               , id="az")
      call alloc( a_H   , NV4*RCW               , id="h")
      call alloc( a_HO  , NFOO*RCW              , id="ho")
      call alloc( a_HV  , NFVV*RCW              , id="hv")
      call alloc( a_HOV , NFVO*RCW              , id="hov")
      call alloc( a_GO  , NFOO*RCW              , id="go")
      call alloc( a_GV  , NFVV*RCW              , id="gv")
      call alloc( a_BUF1,NBUF1*RCW              , id="buf1")
      call alloc( a_BUF2,NBUF2*RCW              , id="buf2")
C
C     Calculate size of the 3rd work array
C     The last array is used to buffer in the <VV||VV> and <VO||VV>
C     integrals, start by making it as large as possible.
C
C     The remaining space is
      call allocator_get_maxbuff(freemem,kind(a_BUF3))
      NBUF3 = freemem / RCW
C
C     Calculate maximum size necessary for the VOVV and VVVV buffers
C
      NV5IRP = 1
      NV6IRP = 1
      DO IRP = 1, NREP
         NV5IRP = MAX(NV5IRP,INT(NVO(IRP),8)*INT(NVVT(IRP),8))
         NV6IRP = MAX(NV6IRP,INT(NVVT(IRP),8)*INT(NVVT(IRP),8))
      ENDDO
C
      NVIRP = MAX(NV5IRP,NV6IRP)
      NBUF3 = MIN(NVIRP,INT(NBUF3,8))
C
C     Check the minimum size for the general use of this array
C
      NBUF3 = MAX(NBUF1,NBUF3)
C
C     Allocate the work array
C
      call alloc( a_BUF3 , NBUF3*RCW , id="buf3")

      call allocator_get_words_available(mxcorr)
C
      IF (ICALC.LE.1) THEN
         WRITE (IW,1001) "Fock space CCSD energies",MXCORR
C
C        Reset symmetry tables to normal (single reference CC use)
C
         EQNS = '      '
         FOCKSP = .FALSE.
         CALL SYMTAB (NELEC,NFROZ,NINACT,NACT,NSP,EPS1,EPS,
     &                IREPSP,IREPSPI,FOCKSP,REPORT,NELEC_F1,NELEC_F2)
         GOTO 440
      ENDIF
C
C     Read the Fock matrix from file
C
      CALL FMTOFILE (.FALSE., a_FVO, a_FOO, a_FVV)
C
C     Solve relativistic Fock space CCSD equations
C
      WRITE (IW,1008) MXITFS,MXDIMFS,10.0D0**(-NTOLFS)
      MAXIT  = MXITFS
      TOL = 10.0D0**(-NTOLFS)
C
C     Initialize T1 and T2 amplitudes
C
! aspg changed
      IF ((IMSTAT(8).ge.2) .or.(IMSTAT(9).ge.2) .or.
     &    (IMSTAT(10).ge.2).or.(IMSTAT(11).ge.2).or.
     &    (IMSTAT(12).ge.2).or.(IMSTAT(13).ge.2)    ) THEN

!        if we already have restartable amplitudes, we read them in
!        and start from there. so far one does not keep track of
!        the full set of amplitudes on a per-sector basis, so the amplitudes
!        recovered here correspond, in fact, to those for the last
!        calculated sector

         write (IW,*)
         write (iw,*) 'INFO: Reading restartable amplitudes!'
         call getampt(a_t1,a_t2)
      ELSE
!        otherwise, we start from mp2 amplitudes
         CALL XCOPY (NDIMT1, a_FVO ,1, a_T1 ,1)
         CALL GETVVOO ( a_T2 )
         CALL DENOMF(EPS, a_T1 , a_T2 , a_T1 , a_T2 ,1)
         CALL DENOMF(EPS, a_T1 , a_T2 , a_T1 , a_T2 ,2)
C     Save the zeroth order wave function
         CALL PUTAMPT ( a_T1 , a_T2 )
         write (IW,*)
         write (iw,*) 'INFO: Initializing amplitudes to MP2 ones!'
      ENDIF
C
 440  CONTINUE
      IFSS(1) = '00'
      IFSS(2) = '01'
      IFSS(3) = '10'
      IFSS(4) = '11'
      IFSS(5) = '02'
      IFSS(6) = '20'
C
C     Solve Fock space equation sector by sector
C
      DO ICURFSS = 1, 6
      IFSTAT=7+ICURFSS
      IF (FSSECT(ICURFSS).NE.0) THEN
        IF (ICALC.LE.1) GOTO 441
        WRITE (IW,1009) IFSS(ICURFSS)
        IF (((IMSTAT(IFSTAT).eq.3).and..not.FS_REDO_SECT(ICURFSS))
     &      .or.FS_SKIP_SECT(ICURFSS)) THEN
           WRITE (IW,*) '... restarting: sector converged or skipped'
        ELSE
           WRITE(IW,'(2X,A,I3)')
     &     '...maximum number of iterations for this sector:',
     &     MAXIT_SECT(ICURFSS)
           if (fs_redo_sect(icurfss).and.(IMSTAT(IFSTAT).eq.3)) then
              write (IW,*)
              write (IW,'(5X,A)') 'INFO: sector previously '//
     &          'marked as converged but user asked to recalculate it'
           endif
        ENDIF


 441    CONTINUE
C
C       skip sectors which have converged, or that we specifically asked to skip
C       conversely calculate sectors we asked to be recalculated even if previously converged
C
        IF ((((IMSTAT(IFSTAT).ne.3).and..not.FS_SKIP_SECT(ICURFSS)).or.
     &      FS_REDO_SECT(ICURFSS)).and.(ICALC.gt.1)) THEN

           maxit_save = MAXIT_SECT(ICURFSS)
          CALL cceqn_driver_amplitudes(maxdim, maxit_save,EPS,
     &             a_FOO , a_FVO , a_FVV  , a_T1   , a_S1   ,
     &             a_T2  , a_S2  , a_AZ   ,a_H,
     &             a_HO  , a_HV  , a_HOV, a_GO, a_GV,
     &             a_BUF1 , a_BUF2 , a_BUF3 ,
     &             TOL, NBUF2, NBUF3)

C         Calculate the energies for this sector
          CALL HEFF1(ICURFSS,IPRNT,a_S1,a_S2,a_FVO,EPS,a_BUF1,
     &         a_BUF2,a_BUF3,NSP,IREPSPI)

          if (maxit_save.gt.maxit_SECT(ICURFSS)) then
              is_ccsd_converged = .false.
          else
             is_ccsd_converged = .true.
          endif

C         Save the updated wave function
          CALL FS_PUTAMPT (ICURFSS, a_T1 , a_T2 , a_BUF1 , a_BUF2 )

          if (is_ccsd_converged.or.restart_unconverged) then
              IMSTAT(IFSTAT) = 3
              write (IW,*)
              write (IW,'(1X,A)') 'INFO:   converged (unconverged'//
     &                            ' if .UNCONVERGED used) amplitudes'//
     &                            ' available for restart'
          else
              IMSTAT(IFSTAT) = 2
              write (IW,*)
              write (IW,'(1X,A)') 'INFO: unconverged amplitudes '//
     &                            'available for restart'
          endif

        ENDIF
C
      ENDIF
      ENDDO
      CALL CCDUMP
C
      call dealloc( a_FVO , id="fvo" )
      call dealloc( a_FOO , id="foo" )
      call dealloc( a_FVV , id="fvv" )
      call dealloc( a_T1  , id="t1" )
      call dealloc( a_T2  , id="t2" )
      call dealloc( a_S1  , id="s1" )
      call dealloc( a_S2  , id="s2" )
      call dealloc( a_AZ  , id="az" )
      call dealloc( a_H   , id="h" )
      call dealloc( a_HO  , id="ho" )
      call dealloc( a_HV  , id="hv" )
      call dealloc( a_HOV , id="hov" )
      call dealloc( a_GO  , id="go" )
      call dealloc( a_GV  , id="gv" )
      call dealloc( a_CCC , id="ccc" )
      call dealloc( a_BB1 , id="bb1" )
      call dealloc( a_BB2 , id="bb2" )
      call dealloc( a_BUF1, id="buf1" )
      call dealloc( a_BUF2, id="buf2" )
      call dealloc( a_BUF3, id="buf3" )
C
      RETURN
1001  FORMAT (" Core used for ",A," :",T50,I20," 8-byte words")
1007  FORMAT (//" Fock space CCSD Energy calculations",
     &        /" Module written by Ephraim Eliav & Luuk Visscher, ",
     &         "July 1999")
1008  FORMAT (//" Fock space CCSD options :"
     &        /" Maximum number of iterations :",T40,I5
     &        /" Maximum size of DIIS space :",T40,I5
     &        /" Convergence criterium :",T38,E7.1)
1009  FORMAT (//" Solving equations for sector ",A2/)
1010  FORMAT (//" Sector ",I2," is not implemented."/)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FS_PUTAMPT (ICURFSS,T1,T2,BUF1,BUF2)
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Augments computed amplitudes with amplitudes from previously
C     calculated sectors and write the updated array to disk
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*),BUF1(*),BUF2(*)
      INTEGER ICURFSS
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*2 DELS1
      CHARACTER*4 DELS2
C
C---------------Executable code--------------------------------------
C
      DELS1 = 'DD'
      IF (ICURFSS.EQ. 2) DELS1 = 'DK'
      IF (ICURFSS.EQ. 3) DELS1 = 'KD'
      IF (ICURFSS.EQ. 4) DELS1 = 'KK'
      IF (ICURFSS.EQ. 5) DELS1 = 'DK'
      IF (ICURFSS.EQ. 6) DELS1 = 'KD'
      DELS2 = 'DDDD'
      IF (ICURFSS.EQ. 2) DELS2 = 'DDDK'
      IF (ICURFSS.EQ. 3) DELS2 = 'KDDD'
      IF (ICURFSS.EQ. 4) DELS2 = 'KDDK'
      IF (ICURFSS.EQ. 5) DELS2 = 'DDKK'
      IF (ICURFSS.EQ. 6) DELS2 = 'KKDD'
C
      IF (ICURFSS.EQ.1) THEN
         CALL PUTAMPT(T1,T2)
      ELSE
         CALL GETAMPT(BUF1,BUF2)
         CALL XAXPY(NDIMT1,A1,BUF1,1,T1,1)
         CALL XAXPY(NDIMT2,A1,BUF2,1,T2,1)
         CALL DELFCK ('VO',DELS1,BUF1)
         CALL DELINT ('VVOO',DELS2,BUF2,0,0,0)
         CALL XAXPY(NDIMT1,-A1,BUF1,1,T1,1)
         CALL XAXPY(NDIMT2,-A1,BUF2,1,T2,1)
         CALL PUTAMPT (T1,T2)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCDINI(do_restart)
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Initialize the dump file and status information.
C
C     Now extended in case of a shared file system: the individual
C     MCCRES files need unique names ! => node number appended to
C     file name. ==> inclusion of ccpar.inc which contains the
C     current MPI data for this specific invocation. We allow for
C     a node number range of 0..999 (three digit node number) which
C     corresponds to a maximum of 1000 nodes.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "results.inc"
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER BUFF,I,IRECL,IRECLEN
      DIMENSION BUFF(RSRCLEN)
      LOGICAL do_restart, file_found
      LOGICAL restart_and_current_differ
      INTEGER RST_NEL(16)
      INTEGER RST_NFRZ(16)
      INTEGER RST_NACT(16)
      INTEGER RST_NINACT(16)
      INTEGER nelec_tot, RST_nelec_tot
      INTEGER act_tot, RST_act_tot
      INTEGER inact_tot, RST_inact_tot
C
C---------------Executable code--------------------------------------
C
C     If MPI is active generate the node-unique file name
C     otherwise use file name without further extension.
C     the global quit routine can be found in gp/gptrygve.F
C
C     set generic RESTART file name
C
      RSBASF='MCCRES'
C
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    !  MPI ID has one digit
         WRITE (RSFILN,'(A6,A1,I1)') RSBASF,'.',MYPROC
         LRSFIL=8
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (RSFILN,'(A6,A1,I2)') RSBASF,'.',MYPROC
         LRSFIL=9
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (RSFILN,'(A6,A1,I3)') RSBASF,'.',MYPROC
         LRSFIL=10
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND CCDINI.F RESTART MODULE")
      ENDIF
#else
      RSFILN=RSBASF
      LRSFIL=6
#endif
C
C     Check the status of the individual restart file
C
      IRECL = IRECLEN (0,RSRCLEN,0)
      INQUIRE (FILE=RSFILN(1:LRSFIL),EXIST=file_found)
C
C     if  the restart file exists ...
C
      IF (file_found) THEN
         if (do_restart) then
C
C        ...and we want to restart we go ahead
C
         write (IW,*)
         write (IW,*) 'INFO: Found restart file(s)!'

         OPEN (MCCRES,FILE=RSFILN(1:LRSFIL),ACCESS='DIRECT',
     &         STATUS='OLD',RECL=IRECL)

C           Read the input information from the previous run
C
             READ(MCCRES,REC=1) RST_NEL,RST_NFRZ,RST_NACT,RST_NINACT

             nelec_tot     = 0
             RST_nelec_tot = 0

             act_tot       = 0
             RST_act_tot   = 0

             inact_tot     = 0
             RST_inact_tot = 0

             restart_and_current_differ = .false.
             DO I = 1, 16
                nelec_tot     = nelec_tot + NELEC(I)
                RST_nelec_tot = RST_nelec_tot + RST_NEL(I)

                act_tot       = act_tot + NACT(I)
                RST_act_tot   = RST_act_tot + RST_NACT(I)

                inact_tot     = inact_tot + NINACT(I)
                RST_inact_tot = RST_inact_tot + RST_NINACT(I)

                IF ((RST_NEL(I)   .NE. NELEC(I)  ).OR.
     &              (RST_NFRZ(I)  .NE. NFROZ(I)  ).OR.
     &              (RST_NACT(I)  .NE. NACT(I)   ).OR.
     &              (RST_NINACT(I).NE. NINACT(I) )) THEN
                   restart_and_current_differ = .true.
                ENDIF
             END DO

             if (restart_and_current_differ) then
                write (IW,*)
                write (IW,*) 'INFO: Current setup and restart differ'
                write (IW,*)
                if (RST_nelec_tot .NE. nelec_tot) then
                write (IW,*) ' # of electrons   (rst):',RST_nelec_tot
                write (IW,*) '                  (now):',nelec_tot
                write (IW,*)
                endif
                if (RST_act_tot .NE. act_tot ) then
                write (IW,*) ' #   act. spinors (rst):',RST_act_tot
                write (IW,*) '                  (now):',act_tot
                write (IW,*)
                endif
                if (RST_inact_tot .NE. inact_tot) then
                write (IW,*) ' # inact. spinors (rst):',RST_inact_tot
                write (IW,*) '                  (now):',inact_tot
                write (IW,*)
                endif

                if (restart_ignore_check) then
                   write (IW,*) 'INFO: Ignoring difference, restarting'
                else
                   call quit('Calculation and restart inconsistent')
                endif
             endif

C
C            Read the status information
C
             READ (MCCRES,REC=2) IMSTAT
C
C            Read the energies (read the whole common block)
C
             READ (MCCRES,REC=3) BUFF
             CALL DCOPY(8,BUFF,1,ETOT,1)

         else
C
C        ... if we don't want to restart we trash any existing files and start over
C
            write (IW,*)
            write (IW,*) 'INFO: Restart file(s) will be overwriten!'
            write (IW,*)
            OPEN (MCCRES,FILE=RSFILN(1:LRSFIL),ACCESS='DIRECT',
     &            STATUS='REPLACE',RECL=IRECL)
            CALL DCOPY(8,A0,0,ETOT,1)
            DO I = 1, NSTATS
               IMSTAT(I) = 0
            ENDDO
         endif
C
C ... if the restart files don't exist, we create them
C
      else
         write (IW,*)
         write (IW,*) 'INFO: No old restart file(s) found!'
         write (IW,*)
         OPEN (MCCRES,FILE=RSFILN(1:LRSFIL),ACCESS='DIRECT',
     &         STATUS='NEW',RECL=IRECL)
         CALL DCOPY(8,A0,0,ETOT,1)
         DO I = 1, NSTATS
            IMSTAT(I) = 0
         ENDDO

      endif
C
C     Initialize the status information
C     When restarting the previous status is retained if the module
C     was on the menu.
C
      IF (DOSORT.AND.IMSTAT(1).EQ.0)             IMSTAT(1) = 1
      IF (DOSORT.AND.IMSTAT(2).EQ.0)             IMSTAT(2) = 1
      IF (DOSORT.AND.IMSTAT(3).EQ.0)             IMSTAT(3) = 1
      IF (DOSORT.AND.IMSTAT(4).EQ.0)             IMSTAT(4) = 1
      IF (DOENER.AND.DOMP2.AND.IMSTAT(1).EQ.0)   IMSTAT(5) = 1
      IF (DOENER.AND.DOCCSD.AND.IMSTAT(1).EQ.0)  IMSTAT(6) = 1
      IF (DOENER.AND.DOCCSDT.AND.IMSTAT(1).EQ.0) IMSTAT(7) = 1
!aspg, restart fs
      IF (DOFSPC.AND.DOMP2.AND.IMSTAT(1).EQ.0)   IMSTAT(5) = 1
      IF (DOFSPC.AND.DOCCSD.AND.IMSTAT(1).EQ.0)  IMSTAT(6) = 1

      IF (DOFSPC.AND.IMSTAT(8).EQ.0)             IMSTAT(8) = 1
      IF (DOFSPC.AND.IMSTAT(9).EQ.0)             IMSTAT(9) = 1
      IF (DOFSPC.AND.IMSTAT(10).EQ.0)            IMSTAT(10) = 1
      IF (DOFSPC.AND.IMSTAT(11).EQ.0)            IMSTAT(11) = 1
      IF (DOFSPC.AND.IMSTAT(12).EQ.0)            IMSTAT(12) = 1
      IF (DOFSPC.AND.IMSTAT(13).EQ.0)            IMSTAT(13) = 1

      IF (restart_redo_sort) then
         DO I = 1, 4
            IMSTAT(I) = 1
         ENDDO
      ENDIF
      if (doccsd.and.restart_redo_ccsd.and.(imstat(6).eq.3)) then
            imstat(6) = 2
            imstat(7) = 1
      endif
!aspg, restart fs
C
C     Initialize the record information
C
      IRECIN = 1 ! Input for the current run
      IRECST = 2 ! Status of the current run
      IRECEN = 3 ! Energies and other results
      IRECT1 = 4 ! T1 amplitudes
C
C     The starting record of the T2 amplitudes is set when the
C     length of the amplitude vectors is known.
C
      CLOSE (MCCRES,STATUS='KEEP')
C
C since a unique restart file name was created in CCDINI we do not
C want to regenerate it in CCDUMP ==> we keep this file name in the
C COMMON BLOCK
C
      CALL CCDUMP
C
      RETURN
 1000 FORMAT (3X,'Step : ',A,' done in previous run; will be skipped')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCDUMP
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Dump all the information in the input and restart block
C     to the begin of the restart file
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "results.inc"
#include "inpt.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER BUFF,IRECL,IRECLEN
      DIMENSION BUFF(RSRCLEN)
C
C---------------Executable code--------------------------------------
C
      IRECL = IRECLEN (0,RSRCLEN,0)
      OPEN (MCCRES,FILE=RSFILN(1:LRSFIL),ACCESS='DIRECT',
     &      STATUS='OLD',RECL=IRECL)
C
C     Dump the input information
C
C     Not implemented yet
      WRITE (MCCRES,REC=1) NELEC, NFROZ, NACT, NINACT
C
C     Dump the status information
C
      WRITE (MCCRES,REC=2) IMSTAT
C
C     Dump the energies (write the whole common block)
C
      CALL DCOPY(8,ETOT,1,BUFF,1)
      WRITE (MCCRES,REC=3) BUFF
C
      CLOSE (MCCRES,STATUS='KEEP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCSUMM
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Write summary of the results, based on the restart file
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "results.inc"
#include "inpt.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER BUFF,I,IRECL,IRECLEN
      DIMENSION BUFF(RSRCLEN)
      CHARACTER DATEX*10,TIMEX*8
      CHARACTER*30 STATCH(0:4)
      DATA STATCH/'Never asked for              ',
     &            'On menu but not yet called   ',
     &            'Started, restartable         ',
     &            'Completed, restartable       ',
     &            'Completed, not restartable   '/
C
C---------------Executable code--------------------------------------
C
      IRECL = IRECLEN (0,RSRCLEN,0)
      OPEN (MCCRES,FILE=RSFILN(1:LRSFIL),ACCESS='DIRECT',STATUS='OLD',
     &      RECL=IRECL)
C
C     Read the status information
C
      READ (MCCRES,REC=2) IMSTAT
C
C     Read the energies (read the whole common block)
C
      READ (MCCRES,REC=3) BUFF
      CALL DCOPY(8,BUFF,1,ETOT,1)
C
      CLOSE (MCCRES,STATUS='KEEP')
C
      CALL DAYTIME (DATEX,TIMEX)
C
      WRITE (IW,1000) ('-',I=1,80)
C
C     Write date & time
C
      WRITE (IW,1001) DATEX,TIMEX
C
C     Write status information
C
      WRITE (IW,1002)
C
      WRITE (IW,1021) STATCH(IMSTAT(1))
      WRITE (IW,1022) STATCH(IMSTAT(2))
      WRITE (IW,1023) STATCH(IMSTAT(3))
      WRITE (IW,1024) STATCH(IMSTAT(4))
      WRITE (IW,1025) STATCH(IMSTAT(5))
      WRITE (IW,1026) STATCH(IMSTAT(6))
      WRITE (IW,1026) STATCH(IMSTAT(7))
C
C     Write energies
C
      if (domp2.or.doccsd) then
         WRITE (IW,1003)
         WRITE (IW,1070) ESCF
         if (domp2)
     &      WRITE (IW,1072) EMP2
         if (doccsd) then
            WRITE (IW,1073) ECCSD
            if (doccsdt) then
               WRITE (IW,1074) ET1
               WRITE (IW,1075) ET2
               WRITE (IW,1076) ET3
            endif
         endif

         if (domp2)
     &      WRITE (IW,1080) ESCF+EMP2

         if (doccsd) then
            WRITE (IW,1081) ESCF+ECCSD
            if (doccsdt) then
               WRITE (IW,1082) ESCF+ECCSD+ET1
               WRITE (IW,1083) ESCF+ECCSD+ET1+ET2
               WRITE (IW,1084) ESCF+ECCSD+ET1+ET3
            endif
         endif
      else
         WRITE (IW,1004)
         if (.not.dofspc) WRITE (IW,1005)
      endif
C
      WRITE (IW,1000) ('-',I=1,80)
C
 1000 FORMAT (//80A)
 1001 FORMAT (//" Today is :",T15,A10/" The time is :",T17,A8)
 1002 FORMAT (/' Status of the calculations')
 1003 FORMAT (/' Overview of calculated energies')
 1004 FORMAT (/' No energies were calculated.')
 1005 FORMAT (/' Check your input!')
 1021 FORMAT (' Integral sort # 1 :',T40,A30)
 1022 FORMAT (' Integral sort # 2 :',T40,A30)
 1023 FORMAT (' Fock matrix build :',T40,A30)
 1024 FORMAT (' MP2 energy calculation :',T40,A30)
 1025 FORMAT (' CCSD energy calculation :',T40,A30)
 1026 FORMAT (' CCSD(T) energy calculation :',T40,A30)
 1070 FORMAT('@ SCF energy :',T40,F25.15)
 1072 FORMAT('@ MP2 correlation energy :',T40,F25.15)
 1073 FORMAT('@ CCSD correlation energy :',T40,F25.15)
 1074 FORMAT('@ 4th order triples correction :',T40,F25.15)
 1075 FORMAT('@ 5th order triples (T) correction :',T40,F25.15)
 1076 FORMAT('@ 5th order triples -T  correction :',T40,F25.15)
 1080 FORMAT('@ Total MP2 energy :',T40,F25.15)
 1081 FORMAT('@ Total CCSD energy :',T40,F25.15)
 1082 FORMAT('@ Total CCSD+T  energy :',T40,F25.15)
 1083 FORMAT('@ Total CCSD(T) energy :',T40,F25.15)
 1084 FORMAT('@ Total CCSD-T  energy :',T40,F25.15)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RDSYMI (NSP,EPS1,IREPSP,ncore_elec)
C
      use interface_to_mpi
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Read symmetry information from MOLFDIR
C
C---------------Routines called----------------------------------------
C
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER IREPSP(*)
      REAL*8 EPS1(*)
      integer ncore_elec,NSP
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "files.inc"
#include "results.inc"
#include "symm.inc"
#include "inpt.inc"
#include "complex.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*14 REPNT(16)
      integer   :: spinor_info(5,2)
      integer   :: i,irep,irp,isp,j,jrep,krep,nfsym,nz
      integer   :: nspc(2,0:2),ncore2
      real(8)   :: toterg
C
C---------------Executable code--------------------------------------
C
C     Read symmetry and spinor information (master node only)
C
      spinor_info(1:5,1:2) = 0
C
      IF (MYPROC.EQ.MASTER) THEN
       OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED')
       READ (MRCONEE,END=10,ERR=10) NSP,BREIT,ECORE,NFSYM,NZ,
     &                              SPFR,NORBT,toterg,nspc,ncore2
       READ (MRCONEE,END=10,ERR=10) NSYMRP,(REPNT(IRP),IRP=1,NSYMRP),
     &                 (NELEC(IRP),IRP=1,NSYMRP),
     &                 (spinor_info(1,irp),irp=1,nfsym),
     &                 (spinor_info(2,irp),irp=1,nfsym),
     &                 (spinor_info(3,irp),irp=1,nfsym),
     &                 (spinor_info(4,irp),irp=1,nfsym),
     &                 (spinor_info(5,irp),irp=1,nfsym)
       ncore_elec = -ncore2
       ncore_elec = ncore_elec+spinor_info(3,1)*2+spinor_info(3,2)*2
       GOTO 11
C---- Keep OLD (MOLFDIR) format here for compatibility
   10     WRITE (IW,*) "Missing information on MRCONEE, I assume"//
     &    " this is a MOLFDIR file and try to continue.."
          INTERFACE='MOLFDIR'
          REWIND(MRCONEE)
          READ (MRCONEE) NSP,BREIT,ECORE
          READ (MRCONEE) NSYMRP,(REPNT(IRP),IRP=1,NSYMRP)
          DO IRP = 1, NSYMRP
             NELEC(IRP) = 0
          ENDDO
   11  CONTINUE
       READ (MRCONEE) NREP,(REPNA(IRP),IRP=1,2*NREP)

c debug, ATT ! these informations are convenient to be printed out
c for debugging purposes. Deactivated in the release.
c      WRITE(IW,*)
c      WRITE(IW,*) '************************************'
c      WRITE(IW,*) '**  MRCONEE CONTROL VARIABLES  *****'
c      WRITE(IW,*) '************************************'
c      WRITE(IW,*)
c      WRITE(IW,*) 'NSP:',NSP
c      WRITE(IW,*) 'BREIT:',BREIT
c      WRITE(IW,*) 'ECORE:',ecore
c      WRITE(IW,*) 'NFSYM:',NFSYM
c      WRITE(IW,*) 'NZ:',NZ
c      WRITE(IW,*) 'SPFR:',SPFR
c      WRITE(IW,*) 'NORBT:',NORBT
c      WRITE(IW,*) 'NSYMRP:',NSYMRP
c      WRITE(IW,*) 'REPNT:',REPNT
c      WRITE(IW,*) 'NELEC:',NELEC
c      WRITE(IW,*) 'NREP:',NREP
c      WRITE(IW,*) 'REPNA:',REPNA
c      WRITE(IW,*)
c      WRITE(IW,*) '************************************'
c      WRITE(IW,*) '**  END CONTROL OUTPUT  ************'
c      WRITE(IW,*) '************************************'
c      WRITE(IW,*)
c debug end

       READ (MRCONEE) ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP)
C      IREPSP is to be used as 2-dimensional array, read it like that
C      IREPSP(ISP)     = IREPSP(ISP,1) contains irrep in parent group
C      IREPSP(NSP+ISP) = IREPSP(ISP,2) contains irrep in Abelian subgroup
       READ (MRCONEE) (IREPSP(ISP),IREPSP(NSP+ISP),
     &                 EPS1(ISP),ISP=1,NSP)
       CLOSE (MRCONEE,STATUS='KEEP')
C
C      Right adjust representation names
C
       DO IRP = 1, NSYMRP
          REPN(IRP) = REPNT(IRP)(1:4)
          CALL BLANKL (4,REPN(IRP))
       ENDDO
       DO IRP = 1, NREP
          CALL BLANKL (4,REPNA(IRP))
       ENDDO
      ENDIF
C
C     A complication using fermion irreps is that 
C     characters are generally complex and that
C     the complex conjugate of an irrep may be another fermion irrep.
C     Only the product of an irrep with its complex conjugate gives 
C     the totally symmetric irrep.
C     
C     Table elements like
C       MULTB(IREP,JREP,1) = KREP
C     provides the symmetry of product functions in either bra of ket positions
C       ket(irep)*ket(jrep) = ket(krep)  OR
C       bra(irep)*bra(jrep) = bra(krep)
C
C     To find the irrep of a product function mixing functions in bra and ket positions
C     we effectively make the inverse of the multiplication table
C     We note that 
C       bra(irep)*ket(irep) = 1 (totally symmetric)
C     This gives
C       bra(irep)*ket(irep)*ket(jrep) = ket(jrep) = bra(irep)*ket(krep)
C       ket(irep)*bra(irep)*bra(jrep) = bra(jrep) = ket(irep)*bra(krep) 
C
      DO JREP = 1, 2*NREP
         DO IREP = 1, 2*NREP
            KREP = MULTB(IREP,JREP,1)
            MULTB(IREP,KREP,2) = JREP
         ENDDO
      ENDDO
C      write(*,*) (REPNA(IREP),IREP=1,2*NREP)

C     write(*,*)
C     write(*,*) '--------------------------------------------'
C     write(*,*) '---  multiplication table (full) -----------'
C     write(*,*) '--------------------------------------------'
C     write(*,*)
C     DO IREP = 1, NREP
C       write(*,'(50I4)') (multb(irep,jrep,1),jrep=1,2*nrep)
C     ENDDO
C     write(*,*) '-------------------------------------'
C     DO IREP = NREP+1, 2*NREP
C       write(*,'(50I4)') (multb(irep,jrep,1),jrep=1,2*nrep)
C     ENDDO
C     write(*,*)
C     write(*,*) '--------------------------------------------'
C     write(*,*) '---  inverse multiplication table (full) ---'
C     write(*,*) '--------------------------------------------'
C     write(*,*)
C     DO IREP = 1, 2*NREP
C       write(*,'(50I4)') (multb(irep,jrep,2),jrep=1,2*nrep)
C     ENDDO

C
C     Get rid of offset for boson reps
C
      DO JREP = 1, NREP
         DO IREP = 1, NREP
            MULTB(IREP,JREP,1) = MULTB(IREP,JREP,1) - NREP
            MULTB(IREP,JREP,2) = MULTB(IREP,JREP,2) - NREP
         ENDDO
      ENDDO
      DO JREP = NREP + 1, 2*NREP
         DO IREP = NREP + 1, 2*NREP
            MULTB(IREP,JREP,1) = MULTB(IREP,JREP,1) - NREP
            MULTB(IREP,JREP,2) = MULTB(IREP,JREP,2) - NREP
         ENDDO
      ENDDO

C     write(*,*)
C     write(*,*) '--------------------------------------------'
C     write(*,*) '---  multiplication table (rebased) --------'
C     write(*,*) '--------------------------------------------'
C     write(*,*)
C     DO IREP = 1, NREP
C       write(*,'(50I4)') (multb(irep,jrep,1),jrep=1,2*nrep)
C     ENDDO
C     write(*,*) '-------------------------------------'
C     DO IREP = NREP+1, 2*NREP
C       write(*,'(50I4)') (multb(irep,jrep,1),jrep=1,2*nrep)
C     ENDDO
C     write(*,*)
C     write(*,*) '-------------------------------------- ------------'
C     write(*,*) '---  inverse multiplication table (rebased) -------'
C     write(*,*) '---------------------------------------------------'
C     write(*,*)
C     DO IREP = 1, 2*NREP
C       write(*,'(50I4)') (multb(irep,jrep,2),jrep=1,2*NREP)
C     ENDDO
C
C     Distribute the information read from MRCONEE
C
#if defined (VAR_MPI)
      IF (NMPROC.GT.1) THEN
         call interface_mpi_BCAST(NSP,1,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST_l0(BREIT,1,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(ECORE,1,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(NSYMRP,1,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(NFSYM,1,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(NZ,1,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST_l0(SPFR,1,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(NORBT,1,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(REPN,NSYMRP*4,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(NELEC,NSYMRP,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(NREP,1,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(REPNA,2*NREP*4,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(MULTB,MXREP*MXREP*8,MASTER,
     &                  global_communicator)
         call interface_mpi_bcast_i1_work_f77(IREPSP,2*NSP,MASTER,
     &                  global_communicator)
         call interface_mpi_bcast_r1_work_f77(EPS1,NSP,MASTER,
     &                  global_communicator)
      END IF
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DISLIST
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     Set up parallel distribution lists for the integral classes VVVV,
C     VOVV and VVOO stored in IDIST(6,6,8) common: ccpar.inc
C     the distribution is performed wrt the second batch
C
C     IOFF counts number of NVVT-sized batches on this node and irrep
C     JOFF points to the starting address in the full integral list
C     KOFF points to the relative starting address in the node integral list
C
C     JOFF and KOFF can be used to go from the absolute indexing to
C     the relative indexing. In both cases the begin is set at zero (last
C     integral of the previous block) and the end to the last integral in
C     the block
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER*8 JOFF,KOFF,I,IOFF,IREP,JREP,K,M,MI
      integer MTOT,mvotot,n,nbtch,no33,nootot,nrest,ntot,nv33
      integer nvoo33,nvotot,nvvo33,nvvtot,off1,off2,off3
C
C---------------Executable code--------------------------------------
C
C ********
C *** VVVV
C ********
C
      KOFF = 0
      DO IREP = 1, NREP
         MTOT = NVVT(IREP)
         NTOT = NVVT(IREP)
         NBTCH = NTOT / NMPROC
         NREST = MOD(NTOT,NMPROC)
         IOFF = 0
         DO I = 1, NMPROC
C
C           Calculate absolute start address of this batch
C
            JOFF = IVVVVTT(IREP) + IOFF * MTOT
C
C           Create even distribution, the higher nodes get the remains
C
            IF (NMPROC-NREST+1.EQ.I) NBTCH = NBTCH + 1
C
C           We need this information for our own process
C
            IF (I.EQ.MYPROC+1) THEN
               IDIST(1,1,IREP) = IOFF
               IDIST(2,1,IREP) = IOFF + NBTCH
               IDIST(3,1,IREP) = JOFF
               IDIST(4,1,IREP) = JOFF + NBTCH * MTOT
               IDIST(5,1,IREP) = KOFF
               IDIST(6,1,IREP) = KOFF + NBTCH * MTOT
               KOFF = KOFF + NBTCH * MTOT
            ENDIF
            IOFF = IOFF + NBTCH
         ENDDO
      ENDDO
C
C     NV6PAR is the number of VVVV integrals which is to handle
C     on THIS node. Can still be too large to be sorted in ONE pass.
C     ==> instead of NV6 the number of passes is calculated from NV6PAR !
C     ==> hierarchy is first determination of individual numbers
C                      second local sorting of VOVV and VVVV classes.
C
      NV6PAR =  KOFF
C
C ********
C *** VOVV
C ********
C
      KOFF = 0
      DO IREP = 1, NREP
         NVOTOT = NVO (IREP)
         NVVTOT = NVVT(IREP)
         NBTCH = NVVTOT / NMPROC
         NREST = MOD(NVVTOT,NMPROC)
         IOFF = 0
         DO I = 1, NMPROC
C
C           Calculate absolute start address of this batch
C
            JOFF = IVOVVT(IREP) + IOFF * NVOTOT
C
C           Create even distribution, the higher nodes get the remains
C
            IF (NMPROC-NREST+1.EQ.I) NBTCH = NBTCH + 1
C
C           We need this information for our own process
C
            IF (I.EQ.MYPROC+1) THEN
               IDIST(1,2,IREP) = IOFF
               IDIST(2,2,IREP) = IOFF + NBTCH
               IDIST(3,2,IREP) = JOFF
               IDIST(4,2,IREP) = JOFF + NBTCH * NVOTOT
               IDIST(5,2,IREP) = KOFF
               IDIST(6,2,IREP) = KOFF + NBTCH * NVOTOT
               KOFF = KOFF + NBTCH * NVOTOT
            ENDIF
            IOFF = IOFF + NBTCH
         ENDDO
      ENDDO
C
C     NV5PAR is the number of VOVV integrals which is to handle
C     on THIS node. Can still be too large to be sorted in ONE pass.
C     ==> instead of NV5 the number of passes is calculated from NV5PAR !
C     ==> hierarchy is first determination of individual numbers
C                      second local sorting of VOVV and VVVV classes.
C
      NV5PAR =  KOFF
C
C
C ****************************************************************
C *** HINTM VVOO contraction, column number in corresponding IRREP
C ****************************************************************
C
      DO IREP = 1, NREP
         JREP=MULTB(IREP+NREP,1+NREP,2)
         MI=MVO(IREP)
         MVOTOT=MVO(JREP)
         NBTCH = MVOTOT / NMPROC
         NREST = MOD(MVOTOT,NMPROC)
         IOFF = J2VOVO(IREP)*RCW + 1
         JOFF = J2VOVO(JREP)*RCW + 1
         DO I = 1, NMPROC

C           Create even distribution, the higher nodes get the remains

            IF (NMPROC-NREST+1.EQ.I) NBTCH = NBTCH + 1

C           Process-local information

            IF (I.EQ.MYPROC+1) THEN
               IDIST(1,4,JREP) = NBTCH
               IDIST(2,4,JREP) = IOFF
               IDIST(3,4,JREP) = JOFF
            ENDIF
            IOFF = IOFF + MI*NBTCH*RCW
            JOFF = JOFF + NBTCH*RCW
         ENDDO
      ENDDO

C
C ****************************************************************
C *** index range distribution for W(A,K,C,I)*T(C,K)
C *** only irrep 1 (totally symmetric ) is needed in the ket space =>
C *** the bosonic bra irrep also has to be 1 and we do not need two
C *** different arrays.
C ****************************************************************
C
      MVOTOT=MVO(1)
      NBTCH = MVOTOT / NMPROC
      NREST = MOD(MVOTOT,NMPROC)
      IOFF = 1
      JOFF = 1
      DO I = 1, NMPROC

C        Create even distribution, the higher nodes get the remains

         IF (NMPROC-NREST+1.EQ.I) NBTCH = NBTCH + 1

C        Process-local information

         IF (I.EQ.MYPROC+1) THEN
            IDIST(1,3,1) = NBTCH
            IDIST(2,3,1) = IOFF
            IDIST(3,3,1) = JOFF
         ENDIF
         IOFF = IOFF + MVOTOT*NBTCH*RCW
         JOFF = JOFF + NBTCH*RCW
      ENDDO

C
C ****************************************************************
C ***  VVOO distr. for AINTM * TAU contraction
C ***  Both arrays have an equal symmetry structure G(i) = G(j)
C ***  suitable for the usage of the CNTRCT call.===> we use this here.
C ***  batches are created according to the individual irreps
C ***  IOFF: source offset    JOFF: start offset of second matrix
C ****************************************************************
C
      DO IREP = 1, NREP
         NOOTOT=NOOT(IREP)
         NVVTOT=NVVT(IREP)
         NBTCH = NOOTOT / NMPROC
         NREST = MOD(NOOTOT,NMPROC)
         IOFF = IVVOOTT(IREP)*RCW + 1
         JOFF = IOOOOTT(IREP)*RCW + 1
         DO I = 1, NMPROC

C           Create even distribution, the higher nodes get the remains

            IF (NMPROC-NREST+1.EQ.I) NBTCH = NBTCH + 1

C           Process-local information

            IF (I.EQ.MYPROC+1) THEN
               IDIST(1,5,IREP) = NBTCH
               IDIST(2,5,IREP) = IOFF
               IDIST(3,5,IREP) = JOFF
            ENDIF
            IOFF = IOFF + NVVTOT*NBTCH*RCW
            JOFF = JOFF + NBTCH*RCW
         ENDDO
      ENDDO

C
C ****************************************************************
C ***  distr. for the contraction TAU * VOOO in T1EQN1.
C ***  the arrays are organised according to
C ***  NV(IREP),NVOOT(IREP) * NVOOT(IREP),NO(IREP)
C ***  Both arrays have an equal symmetry structure G(i) = G(j)
C ***  suitable for the usage of the CNTRCT call.===> we use this here.
C ***  IOFF: source offset for first matrix
C ***  JOFF: entry point for second matrix
C ****************************************************************
C
      DO IREP = 1, NREP
         NV33   = NV(IREP)
         NVOO33 = NVOOT(IREP)
         NBTCH = NVOO33 / NMPROC
         NREST = MOD(NVOO33,NMPROC)
         IOFF = LVVOOT (IREP)*RCW + 1
         JOFF = LOVOOT (IREP)*RCW + 1
         DO I = 1, NMPROC

C           Create even distribution, the higher nodes get the remains

            IF (NMPROC-NREST+1.EQ.I) NBTCH = NBTCH + 1

C           Process-local information

            IF (I.EQ.MYPROC+1) THEN
               IDIST(1,6,IREP) = NBTCH
               IDIST(2,6,IREP) = IOFF
               IDIST(3,6,IREP) = JOFF
            ENDIF
            IOFF = IOFF + NV33*NBTCH*RCW
            JOFF = JOFF + NBTCH*RCW
         ENDDO
      ENDDO

C
C ****************************************************************
C ***  distr. for the contraction T2 * GOINTM in T2EQNS
C ***  the arrays are organised according to
C ***  NVVOT(IREP),NO(IREP) * NO(IREP),NO(IREP)
C ***  Both arrays have an equal symmetry structure G(i) = G(j)
C ***  suitable for the usage of the CNTRCT call.===> we use this here.
C ***  IOFF: source offset for first matrix
C ***  JOFF: entry point for second matrix
C ****************************************************************
C
      DO IREP = 1, NREP
         NVVO33 = NVVOT(IREP)
         NO33   = NO(IREP)
         NBTCH = NO33 / NMPROC
         NREST = MOD(NO33,NMPROC)
         IOFF = KVVOOT(IREP)*RCW + 1
         JOFF = IOO(IREP)*RCW + 1
         DO I = 1, NMPROC

C           Create even distribution, the higher nodes get the remains

            IF (NMPROC-NREST+1.EQ.I) NBTCH = NBTCH + 1

C           Process-local information

            IF (I.EQ.MYPROC+1) THEN
               IDIST(1,7,IREP) = NBTCH
               IDIST(2,7,IREP) = IOFF
               IDIST(3,7,IREP) = JOFF
            ENDIF
            IOFF = IOFF + NVVO33*NBTCH*RCW
            JOFF = JOFF + NBTCH*RCW
         ENDDO
      ENDDO

C ****************************************************************
C ***  distr. for the contraction T1(A,K) * W(KB,IJ)
C ***  the arrays are organised according to
C ***  NV(IREP),NO(IREP) * NO(IREP)*NVOOT(IREP)
C ***  Both arrays have an equal symmetry structure G(i) = G(j)
C ***  suitable for the usage of the CNTRCT call.===> we use this here.
C ***  IOFF: source offset for first matrix
C ***  JOFF: entry point for second matrix
C ****************************************************************
C
      DO IREP = 1, NREP
         NV33 = NV(IREP)
         NO33 = NO(IREP)
         NBTCH = NO33 / NMPROC
         NREST = MOD(NO33,NMPROC)
         IOFF = IVO(IREP)*RCW + 1
         JOFF = LOVOOT(IREP)*RCW + 1
         DO I = 1, NMPROC

C           Create even distribution, the higher nodes get the remains

            IF (NMPROC-NREST+1.EQ.I) NBTCH = NBTCH + 1

C           Process-local information

            IF (I.EQ.MYPROC+1) THEN
               IDIST(1,8,IREP) = NBTCH
               IDIST(2,8,IREP) = IOFF
               IDIST(3,8,IREP) = JOFF
            ENDIF
            IOFF = IOFF + NV33*NBTCH*RCW
            JOFF = JOFF + NBTCH*RCW
         ENDDO
      ENDDO

C ****************************************************************
C ***  distr. for the contraction T1(A,K) * W(BK,IC) * T1(C,J)
C ***  the permutation operators are taken into account by the sorters
C ***  CNTRCT organization: NVOO(IRP),NV(IRP) * NV(IRP),NO(IRP)
C ***  Both arrays have an equal symmetry structure G(i) = G(j)
C ***  suitable for the usage of the CNTRCT call.===> we use this here.
C ***  IOFF: source offset for first matrix
C ***  JOFF: entry point for second matrix
C ****************************************************************
C
      DO IREP = 1, NREP
         NVOO33 = NVOO(IREP)
         NV33   = NV(IREP)
         NBTCH = NV33 / NMPROC
         NREST = MOD(NV33,NMPROC)
         IOFF  = KVOOV(IREP)*RCW + 1
         JOFF  = IVO(IREP)*RCW + 1
         DO I = 1, NMPROC

C           Create even distribution, the higher nodes get the remains

            IF (NMPROC-NREST+1.EQ.I) NBTCH = NBTCH + 1

C           Process-local information

            IF (I.EQ.MYPROC+1) THEN
               IDIST(1,9,IREP) = NBTCH
               IDIST(2,9,IREP) = IOFF
               IDIST(3,9,IREP) = JOFF
            ENDIF
            IOFF = IOFF + NVOO33*NBTCH*RCW
            JOFF = JOFF + NBTCH*RCW
         ENDDO
      ENDDO

C
C
C *******************************
C *** HINTM  contraction in T2EQN
C *******************************
C
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP=MULTB(IREP+NREP,1+NREP,2)
         M=MVO(IREP)
         N=MVO(JREP)
         K=MVO(IREP)
         MVOTOT=MVO(IREP)
         NBTCH = MVOTOT / NMPROC
         NREST = MOD(MVOTOT,NMPROC)
         IOFF = OFF1
         JOFF = OFF2
         DO I = 1, NMPROC

C           Create even distribution, the higher nodes get the remains

            IF (NMPROC-NREST+1.EQ.I) NBTCH = NBTCH + 1

C           Process-local information

            IF (I.EQ.MYPROC+1) THEN
               IDIST(1,10,JREP) = NBTCH
               IDIST(2,10,JREP) = IOFF
               IDIST(3,10,JREP) = JOFF
            ENDIF
            IOFF = IOFF + M*NBTCH*RCW
            JOFF = JOFF + NBTCH*RCW
         ENDDO
C
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
      RETURN
      END

      subroutine get_charge_cc(mycharge,nelect_cc,ncore_elec_cc)
      implicit none
      integer, intent(inout) :: mycharge
      integer, intent(in)    :: nelect_cc
      integer, intent(in)    :: ncore_elec_cc
      integer                :: charge
      charge   = 0
      call rmolchr(charge)
      mycharge = charge - (ncore_elec_cc + nelect_cc)
!     print *, 'charge, ncore_elec_cc, nelect_cc, mycharge',
!    & charge, ncore_elec_cc, nelect_cc, mycharge
      end subroutine


