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
C  /* Deck psiscf */
      SUBROUTINE PSISCF(WFCONV)
C***********************************************************************
C
C     Solve SCF - equations (HF/KS); print results
C
C     This routine is called from:  PAMPSI (main/dirac.F)
C
C     Written by T.Saue November 1994
C
C     Last revisions: Sep 11 1995 - tsaue
C                     2005/6 MI to make BSS2DC/DC2BSS available
C
C***********************************************************************

         use memory_allocator
         use dirac_cfg
         use x2cmod_cfg
         use xmlout

#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbpsi.h"
C
      LOGICAL   WFCONV, DO4C2C_SAVE, EX,
     &          TRIVEC_SAVE, TRIFCK_SAVE, ERGCNV_SAVE,
     &          EVCCNV_SAVE, FCKCNV_SAVE 
      real(8), allocatable :: WORK(:)
C
      CALL QENTER('PSISCF')
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in PSISCF')
C
      WFCONV = .FALSE.
      if (doxml) then
         call xml_begin('task','type="SCF optimization"')
         call xml_begin('input')
         call xml_comment('method, convergence thresholds etc')
         call xml_end('input')
         call xml_begin('output')
      endif
      IF(dirac_cfg_dft_calculation) THEN
        IF (BSS) THEN
          CALL TITLER('Two-component BSS/DKH Kohn-Sham calculation',
     &       '*',125)
        ELSE IF (x2c) THEN
          CALL TITLER('X2C Kohn-Sham calculation',
     &       '*',125)
        ELSE
           CALL TITLER('Kohn-Sham calculation','*',125)
        ENDIF
      ELSEIF(ONESYS) THEN
        CALL TITLER('One-electron system','*',125)
      ELSE
        IF (BSS) THEN
           CALL TITLER(
     &     'Two-component BSS/DKH relativistic HF calculation',
     &     '*',125)
        ELSE IF (x2c) THEN
          CALL TITLER('X2C relativistic HF calculation',
     &       '*',125)
        ELSE
           CALL TITLER('Hartree-Fock calculation','*',125)
        ENDIF
      ENDIF
      CALL FLSHFO(LUPRI)
C
C     Initialize
C     ==========
C
      CALL SETDHF(IPRSCF)
C
C     Memory allocation
C
C     KCMO is N2BBASXQ because it is used for Fock matrix in PREDHF.
C
      CALL MEMGET2('REAL','CMO  ',KCMO ,N2BBASXQ,      WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EIG  ',KEIG ,NTBAS(0),      WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBEIG',KIBE ,NTBAS(0),      WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','FOCK ',KFOCK,N2BBASXQ*NFMAT,WORK,KFREE,LFREE)
C
C     Solve DHF-equations
C     ===================
C
      DO4C2C_SAVE = DO4C2C
      DO4C2C = .FALSE.
      IF (DO2C4C.AND.(BSS.or.x2c)) THEN
CMI    ... prepare for the preliminary BSS-SCF procedure to get
CMI   starting Fock-Dirac MO matrix for the DC-SCF
        ERGCNV_SAVE = ERGCNV
        EVCCNV_SAVE = EVCCNV
        FCKCNV_SAVE = FCKCNV 
        MAXITR_SAVE = MAXITR

        EVCCNV = EVCCNV2
        ERGCNV = ERGCNV2
        FCKCNV = FCKCNV2
        MAXITR = MAXITR2

        SCFCNV1_SAVE = SCFCNV(1)
        SCFCNV2_SAVE = SCFCNV(2)

        SCFCNV(1) = SCFCNV2(1)
        SCFCNV(2) = SCFCNV2(2)

      ENDIF

      CALL DHFSCF(WORK(KCMO),WORK(KEIG),WORK(KIBE),WORK(KFOCK),
     &            WORK(KFREE),LFREE)

      !> in the hypothetical case that we run a property calculation for a one-electron system 
      !> it is useful to define some 2e-common blocks... strange but it is a 'bug' in the post-SCF module setup
      !> which basically always assumes 2e-ints to be present.
  
      if(onesys) call paovec(work(kfree),lfree,0,0)
C
C     Memory deallocation of the Fock matrix only...
      CALL MEMREL('PSISCF.DHFSCF',WORK,KWORK,KFOCK,KFREE,LFREE)
C
C     Output section
C
      CALL DHFOUT(WORK(KCMO),WORK(KEIG),WORK(KIBE),WORK(KFREE),
     &            LFREE,IPRSCF)
CMI  ....  after 2cIOTC-SCF prepare starting data for subsequent 4c SCF
      IF (DO2C4C.AND.(INI2C.EQ.3.OR.INI2C.EQ.4)) THEN
        CALL MEMGET2('REAL','TBUF',K1,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','TMAT',K2,N2BBASXQ,WORK,KFREE,LFREE)
        CALL CMOTR4C(WORK(K1),WORK(K2),WORK(KCMO),WORK(KEIG),
     &       WORK(KIBE),WORK(KFREE),LFREE)
        CALL MEMREL('PSISCF.CMOTR4C',WORK,KWORK,K1,KFREE,LFREE)
      ENDIF
C
C     Memory deallocation of the rest..
C
      CALL MEMREL('PSISCF.DHFOUT',WORK,KWORK,KWORK,KFREE,LFREE)

      call dealloc(work)

C=============================================================
CMI ... after BSS-SCF turn to the four-component picture.....
C=============================================================
      IF (DO2C4C) THEN          ! 2c --> 4c switch
        INTDEF    = INTGEN_SAVE
        INTGEN    = INTDEF 
        SSMTRC    = SSMTRC_SAVE
        ERGCNV = ERGCNV_SAVE
        EVCCNV = EVCCNV_SAVE
        FCKCNV = FCKCNV_SAVE
        SCFCNV(1) = SCFCNV1_SAVE
        SCFCNV(2) = SCFCNV2_SAVE
        MAXITR = MAXITR_SAVE
        CALL GO2C4C()
        DO2C4C = .FALSE.
      ENDIF

      DO4C2C = DO4C2C_SAVE

      IF (DO4C2C) THEN ! 4c --> 2c switch
         CALL GO4C2C()
      END IF 
C
C     Check if DHF/DFT is converged (to hard or soft threshold):
C
      WFCONV = DHFCONV(1) .OR. DHFCONV(2)
C
      CALL FLSHFO(LUPRI)
      if (doxml) then
         call xml_end('output')
         call xml_end('task')
      endif
      if(allocated(work)) call dealloc(WORK)
      CALL QEXIT('PSISCF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck predhf */
      SUBROUTINE PREDHF(CMO,EIG,IBEIG,FOCK,DMAT,BMAT,NPOS,WORK,LWORK)
C*****************************************************************************
C
C     Obtain orthonormal trial vectors for Dirac-Fock
C     either by:
C       i) bare nucleus approximation
C      ii) starting from trial vectors
C     iii) starting from two-electron Fock matrix
C     Also serves as solver module for the Dirac-equation in the
C     finite basis approximation (ONESYS)
C
C     Numerous possible schemes for restart has been implemtented.
C     They are based on the program looking for necessary files.
C     Three modes should be discerned:
C
C     1) One-electron system (ONESYS)
C     2) Many-electron system
C         a) DFCYCL present:
C            The program will restart on as much as possible.
C            (coefficients, convergence acceleration,
C             differential density etc.)
C         b) DFCYCL abscent:
C            Restart is only possible on coefficients/Fock matrix.
C
C     Written by T.Saue November 1994
C     Last revision : Jan 27 1997 - tsaue
C                     July 2005 - MI - adapted for DC2BSS calculations.
C
C*****************************************************************************

      use dirac_cfg
      use interface_mo_specific
      use fde_mod
      use num_grid_gen
      use xmlout
#ifdef HAS_PCMSOLVER
      use pcm_scf
#endif

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0,D0 = 0.0D0, D4 = 4.00D00,
     &          D2 = 2.0D0, DP5=0.5D0)
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbreo.h"
      type(fde_import) :: itmp
      LOGICAL DOSEL,DOLEV,DSCFON,AOCBUF
      CHARACTER CPUTID*12,WALTID*12,DAYTID*24,SECTID*12
      INTEGER NOCCUP(8)
      DIMENSION FOCK(*),EIG(*),IBEIG(*),CMO(*),
     &          BMAT(MXDIIS,MXDIIS),DMAT(*),
     &          NPOS(*),WORK(LWORK),
     &          DVOVLP(MXOPEN)
      integer, allocatable :: IBEIG2(:,:)
#ifdef HAS_PCMSOLVER
      real(8), allocatable :: fock_pcm(:, :, :)
#endif

      real(8) :: TIMSTR(2), TIMEND(2)

C
      CALL QENTER('PREDHF')
C
#include "memint.h"
C     
C     ****************************************************
C     ***  I N I T I A L I Z E  a n d   R E S T A R T ****
C     ****************************************************
C
      TDXC=D0
      CALL DHFTST(DSCFON,ICHKCNV,BMAT,IPRSCF)
CMI    no bare nucleus screening when 2c-4c/4c-2c SCF switch     
      IF (RESTFCK) THEN
         BARNUC  =.FALSE.
         DOHUCKEL=.FALSE.
      END IF

      IF (.NOT.ONESYS) THEN
C
C       Initiate SCF cycle
C
        CALL GTINFO(DAYTID)
        WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +  '########## START ITERATION NO.',NITER,' ##########',
     +  DAYTID
        CALL FLSHFO(LUPRI)
C
C       Open files necessary for Dirac-Fock calculation
C
!        N2TMOTQ = N2TMO(I)*NZ        
!        N2TMOTQ = N2TMO(I)*NZ        
        OPEN (LUEVEC, FILE='DFEVEC', FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*N2TMOTQ,STATUS='UNKNOWN')
        CALL OPNFIL(LUFCK2,'DFFCK2','UNKNOWN','PREDHF')
        CALL OPNFIL(LUFCKT,'DFFCKT','UNKNOWN','PREDHF')
        CALL OPNFIL(LUDENS,'DFDENS','UNKNOWN','PREDHF')
      ENDIF  ! IF (.NOT.ONESYS) THEN
      CALL OPNFIL(LUCMOS,'DFCMOS','UNKNOWN','PREDHF')
      OPEN (LUFOCK, FILE='DFFOCK', FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*N2TMOTQ,STATUS='UNKNOWN')
      CALL OPNFIL(LUCOEF,'DFCOEF','UNKNOWN','PREDHF')
      CALL OPNFIL(LUFCK1,'DFFCK1','OLD','PREDHF')
C
C     *****************************************************************
C     *****  G E N E R A T E   F U L L    F O C K    M A T R I X  *****
C     *****************************************************************
C
      CALL FLSHFO(LUPRI)
      CALL GETTIM(CTIDSTR,WTIDSTR)
      IF (IPRSCF .GE. 3) CALL TIMER2('START',TIMSTR,TIMEND)

C
C     Get 2-electron Fock matrix (save in CMO)
C     ========================================
      DOSEL = .FALSE.
C ... skip the generation of 2-el.Fock matrix when RESTFCK=.true. 
      IF (.NOT.RESTFCK) THEN
      IF(ATHUCK) THEN
        DHF_INTTYP = 'Atom. Huckel'
        CALL ADHOC(CMO,EIG,IBEIG,DHFERG,IPRSCF,WORK,KFREE,LFREE)
        CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,DHFERG,.FALSE.) ! CTROND: check
        GOTO 100
C
C     2) - atomic start
C     -------------------------
      ELSEIF(ATOMST) THEN
C.......In the open-shell case use fractional occupation in first iteration; reset at the end
        DHF_INTTYP = 'Atomic s'
        AOCBUF     = AOC
        AOC        = .FALSE.
        NFMATBUF   = NFMAT
        NFMAT      = 1
        CALL ATOMIC_START(DMAT,IPRSCF,WORK,KFREE,LFREE)
        !fixme: michal for InteRest                            
        CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,DMAT,NFMAT,
     &                NPOS,INTFLG,IPRSCF,WORK(KFREE),LFREE)
C
C     2) - using trial vectors
C     ------------------------
C
      ELSE IF (TRIVEC) THEN
C
C       We calculate a regular Fock matrix so we can check convergence ...
C
        IF (ICHKCNV.EQ.0) ICHKCNV = -1
        IOPT = 15
        IF(ATOMIC) THEN
          CALL ATOMIC_REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,IOPT)           
        ELSE
          CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,IOPT)
        ENDIF
C
C       Generate untransformed coefficients (MO basis)
C       and prepare for overlap selection
C
        CALL PRECMOS(CMO,DMAT,FOCK,DOSEL,IPRSCF,WORK,KFREE,LFREE)
C
C       Determine occupation if this was not given in input
C
        IF (INIOCC) THEN
           CALL DETOCC (EIG,WORK,KFREE,LFREE)
           INIOCC = .FALSE.
        ENDIF
C
C       Generate density matrix/matrices
C       Be aware: All density matrices are normalized to 1,
C                 occupation is in the DF(*) factor
C
        CALL DENMAT(DMAT,CMO,IPRSCF)
C
C       Print Mulliken gross populations, if requested
C
        IF (SCFPOP) CALL SCF_MULPOP(CMO,EIG,IBEIG,WORK(KFREE),LFREE)
C
C       Differential density matrix approach
C
        TDF2 = SECOND()
        IF(DSCFON) THEN
          WRITE(LUPRI,'(I5,A)') NITER,
     &       ' *** Differential density matrix ***'
          CALL READNS(LUDENS,FOCK)
          CALL WRIDNS(LUDENS,DMAT)
          IF (NISHT .GT. 0) THEN
            IF(FIXDIF) THEN
              DCOVLP = D1
            ELSE
              DCOVLP = DDOT(N2BBASXQ,FOCK,1,FOCK,1)
              DCOVLP = DDOT(N2BBASXQ,DMAT,1,FOCK,1) / DCOVLP
            ENDIF
            CALL DAXPY(N2BBASXQ,(-DCOVLP),FOCK,1,DMAT,1)
          END IF
          IF (AOC) THEN
             DO IOPEN = 1,NOPEN
                IF (FIXDIF) THEN
                   DVOVLP(IOPEN) = D1
                ELSE
                   DVOVLP(IOPEN) =
     &                DDOT(N2BBASXQ,FOCK(1+N2BBASXQ*IOPEN),1,
     &                              FOCK(1+N2BBASXQ*IOPEN),1)
                   DVOVLP(IOPEN) =
     &                DDOT(N2BBASXQ,DMAT(1+N2BBASXQ*IOPEN),1,
     &                FOCK(1+N2BBASXQ*IOPEN),1) / DVOVLP(IOPEN)
                END IF
                WRITE(LUPRI,'(I5,A,I2,A,F6.4)') NITER,
     &             ' *** Differential density matrix. DVOVLP(',
     &             IOPEN,') =',DVOVLP(IOPEN)
                CALL DAXPY(N2BBASXQ,(-DVOVLP(IOPEN)),
     &               FOCK(1+N2BBASXQ*IOPEN),1,DMAT(1+N2BBASXQ*IOPEN),1)
             END DO
          END IF
C
C           Note that all Fock matrices are multiplied with 2
C              in TWOFCK.
C
          !fixme: michal for InteRest                            
          CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,DMAT,NFMAT,
     &                NPOS,INTFLG,IPRSCF,WORK(KFREE),LFREE)
C           read old Fock matrix in DMAT to finish new Fock matrix/hj
          CALL REAFCK(LUFCK2,DMAT,.FALSE.,NFMAT)
          IF (NISHT.GT.0)
     &         CALL DAXPY(N2BBASXQ,DCOVLP,DMAT,1,FOCK,1)
          IF (AOC) THEN
            DO IOPEN = 1,NOPEN
              CALL DAXPY(N2BBASXQ,DVOVLP(IOPEN),
     &                   DMAT(1+N2BBASXQ*IOPEN),1,
     &                   FOCK(1+N2BBASXQ*IOPEN),1)
            END DO
          END IF
C         restore density matrix in DMAT/tsaue
          CALL READNS(LUDENS,DMAT)
        ELSE ! DSCFON
          CALL WRIDNS(LUDENS,DMAT)
          !fixme: michal for InteRest                            
          CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,DMAT,NFMAT,
     &                NPOS,INTFLG,IPRSCF,WORK(KFREE),LFREE)
        END IF ! DSCFON
C
C
        IF(DOTSC) CALL TSCORR(FOCK,WORK(KFREE),LFREE,IPRSCF)
C
C
C       Write current two-electron Fock matrix to file
C       Note that possible DFT/PCM contributions are not 
C       included since DFFCK2 is used for completing two-electron
C       Fock matrices generated by differential  densities
C       ==============================================
C
        CALL WRIFCK(LUFCK2,FOCK,NFMAT)
C
        if (dirac_cfg_fde) then
           call fde_get_import_info(itmp)
           if (itmp%im_update_vemb) then
              tdfde = SECOND()
              call interface_mo_write()
              call fde_dirac_set_nz(nz)
#ifdef VAR_MPI
           if (parcal) call dirac_parctl(FDE_PAR)
#endif
              call fde_calculate_emb_pot_mat(ntbas(0),dmat,fock)
              tdfde =  second() - tdfde
           endif
        endif

        tdxc = second()
!       xc contribution
        if (dirac_cfg_dft_calculation) then
           call interface_mo_write()
           call generate_num_grid(dmat)
#ifdef VAR_MPI
           if (parcal) call dirac_parctl( XCINT_PAR )
#endif
           call xcint_potential_rks(ntbas(0),
     &                              dmat,
     &                              fock)
           tdxc = second() - tdxc
        end if

C
C       Add solvent contribution
C
        IF (SOLVEN) THEN
           CALL SOLFCK(FOCK,DMAT,1,ESOLVE,ESOLVN,
     &                 WORK(KFREE),LFREE,IPRSOL)
        ENDIF

C
C     3) - using old 2-electron Fock matrix
C     --------------------------------------
      ELSE IF (TRIFCK) THEN
C
CMI       ... works only for closed-shell systems !
         CALL REAFCK(LUFCK2,FOCK,.TRUE.,NFMAT)
C
C     4) - Huckel start guess
C     --------------------------------------
C
      ELSE IF (DOHUCKEL) THEN
C        ... read 1-electron Fock matrix into FOCK ...
         CALL REAFCK(LUFCK1,FOCK,.FALSE.,1)
         NELECT = NELECT_DHF
         CALL IZERO(NOCCUP,8)  ! non-rel. occupation of boson symmetries
         IF (SPINFR) NOCCUP(1:8) = NISH_BOS(1:8)
         CALL HUCDRV(FOCK,NTBAS(0),NELECT,NOCCUP,
     &               WORK(KFREE),LFREE,IPRSCF)
C
C     5) Get one-electron Fock matrix (Bare Nuclei)
C     --------------------------------------
C
      ELSE IF (BARNUC) THEN
C       ... read 1-electron Fock matrix into FOCK ...
         CALL REAFCK(LUFCK1,FOCK,.FALSE.,1)
        IF(IPRSCF.GE.12) THEN
         CALL HEADER(
     &   'PREDHF: Fresh 1-electron Fock matrix within "BARNUC"',-1)
         CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
!       ... add some kind of screening potential
        IF(BNCRON .AND. .NOT. ONESYS) THEN
          CALL BNCORR(FOCK,WORK(KFREE),LFREE,IPRSCF)
          IF (IPRSCF .GE. 3) call TIMER2('BNCORR old',TIMSTR,TIMEND)
        ELSEIF (BNSPON .AND. .NOT. ONESYS) THEN
          CALL Add_Screening_Potential(FOCK,WORK(KFREE),LFREE,IPRSCF)
          IF (IPRSCF .GE. 3) call TIMER2('BNCORR SCRPOT',TIMSTR,TIMEND)
        ENDIF
        IF(IPRSCF.GE.12) THEN
          CALL HEADER(
     & 'PREDHF: 1-electron Fock matrix + BNCORR within "BARNUC"',-1)
           CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                 NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
      END IF

! =======================================================================

      IF (.NOT. BARNUC) THEN
C     MI     ... get 1-electron Fock matrix into CMO ...
         CALL REAFCK(LUFCK1,CMO,.FALSE.,1)
C
C     Calculate correction according to the One-center app. models 
C     ============================================================
C
        IF (TRIVEC.OR.ATOMST) THEN
           IF (ONECAP.AND. .NOT.SMLV1C) THEN
              CALL CATCORR(CMO,FOCK,DMAT,WORK(KFREE),LFREE,IPRSCF)
           END IF
           CALL ERGCAL(CMO,FOCK,DMAT,WORK(KFREE),LFREE)
        END IF

#ifdef HAS_PCMSOLVER
        if (dirac_cfg_pcm) then
           do ifock = 1, n2bbasxq
              if (fock(ifock) /= fock(ifock)) then
                 call quit('NaN found in Fock matrix')
              end if 
           end do
           WRITE(lupri,'(a)') '* Calling pcm_oper_ao_driver'
           allocate(fock_pcm(ntbas(0), ntbas(0), nz))
           fock_pcm = 0.0d0
! fock is the 2-el Fock matrix
           call pcm_oper_ao_driver(fock_pcm, 'TotASC'//char(0), 
     &          work(kfree), lfree)
           call daxpy(n2bbasx, -1.0d0, fock_pcm, 1, fock, 1)
!     print *, "FOCK_PCM"
!     call prqmat(fock_pcm, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
!                      nz, ipqtoq(1,0), 6)
           deallocate(fock_pcm)
           do ifock = 1, n2bbasxq
              if (fock(ifock) /= fock(ifock)) then
                 call quit('NaN found in Fock matrix after pcm_fock')
              end if 
           end do
        endif 
#endif

CMI      ... add 1-electron Fock matrix into the total Fock matrix ...
        CALL DAXPY(N2BBASXQ,D1,CMO,1,FOCK,1)
C
C       Add active two-electron Fock matrices
C
        IF(AOC) THEN
          DO IOPEN = 1,NOPEN
C
C           The factor in only DF and not D2*DF because of the
C           multiplication with D2 in TWOFCK.
C
            CALL DAXPY(N2BBASXQ,DF(IOPEN),
     &                 FOCK(1+N2BBASXQ*IOPEN),1,FOCK,1)
          END DO
        ENDIF
C
C       Write total Fock matrix
C
        CALL WRIFCK(LUFCKT,FOCK,NFMAT)
        IF(IPRSCF.GE.5) THEN
          IF (NISHT .GT. 0) THEN
            CALL HEADER('PREDHF: Total FD matrix in SA-AO basis',-1)
            CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          END IF
          IF(AOC) THEN
            DO IOPEN = 1,NOPEN
              CALL HEADER('PREDHF: Total FV matrix in SA-AO basis',-1)
              WRITE(LUPRI,'(A,I3)') ' Matrix number',IOPEN
              CALL PRQMAT(FOCK(1+N2BBASXQ*IOPEN),NTBAS(0),NTBAS(0),
     &                    NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
            END DO
          ENDIF
        END IF
      END IF ! (.not. BARNUC)
C
      ENDIF ! (.not. RESTFCK)
C
C     *************************************************
C     *****   SOLVE  DHF  EIGENVALUE  PROBLEM    ******
C     *************************************************
C
      IF(.NOT.TRIVEC) THEN
        DOLEV = .FALSE.
      ELSE
        DOLEV = DOLEVEL
      ENDIF
      CALL DFSOLV(CMO,EIG,IBEIG,FOCK,DMAT,BMAT,ICHKCNV,
     &            DOSEL,DOLEV,WORK(KFREE),LFREE)
      IF(ATOMST) THEN
        AOC      = AOCBUF
        NFMAT    = NFMATBUF
      ENDIF
C
C     Add vacuum polarisation in case of ONESYS
C
#ifdef MOD_UNRELEASED
      IF(ONESYS.AND.VACPOL) THEN
        WRITE(LUPRI,'(/A)') '=> Adding vacuum polarization potential'
        CALL ONEVAC(CMO,EIG,IBEIG,FOCK,WORK(KFREE),LFREE)
      ENDIF
#endif
CMI   ... when (DO4C2C.AND..NOT.CONT2C)=.true.,  !
      IF (DOHUCKEL.OR.BARNUC.OR.(TRIFCK.AND..NOT.DO4C2C)) THEN
        WRITE(LUPRI,'(/A)') '=> Calculating sum of orbital energies'
        CALL ONEERG(EIG,IBEIG)
      END IF

      IF (DO4C2C.AND..NOT.ONESYS) THEN
        IF (.NOT.CONT2C) THEN
         IF (IPRHAM.GE.2.OR.IPRSCF.GE.3) THEN
          WRITE(lupri,'(/2X,A/)')
     &  'PREDHF: Last calculation of the SCF energy'//
     &  ' using pict.change transformed Fock matrixes'//
     &  ' from LUFCK1 and LUFCK2.'
         ENDIF
         CALL MEMGET2('REAL','FCK1',KFCK1,N2BBASXQ,WORK,KFREE,LFREE)
         CALL REAFCK(LUFCK1,WORK(KFCK1),.TRUE.,1)
         CALL REAFCK(LUFCK2,FOCK,.TRUE.,NFMAT)
         CALL DENMAT(DMAT,CMO,IPRSCF)
         CALL WRIDNS(LUDENS,DMAT)
         CALL ERGCAL(WORK(KFCK1),FOCK,DMAT,WORK(KFREE),LFREE)
         CALL MEMREL('PREDHF',WORK,KFCK1,KFCK1,KFREE,LFREE)
        ENDIF
        RESTFCK = .FALSE.
CMI        ... MAXITER is what does the SCF procedure end !!
        IF (CONT2C) THEN
         DHFCONV(1)=.FALSE. 
         DHFCONV(2)=.FALSE. 
            DHFEXIT=.FALSE. 
        ENDIF
      ENDIF
C
 100  CONTINUE
      CALL GETTIM(CTIDEND,WTIDEND)
      CTIDTOT = CTIDEND - CTIDSTR
      WTIDTOT = WTIDEND - WTIDSTR
C
      IF(ONESYS) THEN
        IF(PHCOEF) THEN
          DO I = 1,NFSYM
            CALL PHATRA(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ)
          ENDDO
          IF(ATOMIC) THEN
            CALL ATOMIC_WRICMO(LUCOEF,CMO,EIG,IBEIG,DHFERG)
          ELSE
            CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,DHFERG,.FALSE.)
          ENDIF   
        ENDIF
        CLOSE(LUCOEF,STATUS='KEEP')
        CLOSE(LUFCK1,STATUS='KEEP')
        CLOSE(LUFOCK,STATUS='KEEP')
        CPUTID = SECTID(WTIDTOT)
        WRITE(LUPRI,'(A,A12)')
     &    '* Total time for one-electron system:',CPUTID
      ELSE
        ITRSCF(INTFLG) = ITRSCF(INTFLG) + 1
        SCFTID(INTFLG) = SCFTID(INTFLG) + WTIDTOT
        IF(IPRSCF.GE.1) THEN
          CPUTID = SECTID(TDF2)
          WRITE(LUPRI,'(3X,A,A12)')
     &      '* Two-electron Fock matrix  (CPU):',CPUTID
          IF(dirac_cfg_dft_calculation) THEN
          CPUTID = SECTID(TDXC)
          WRITE(LUPRI,'(3X,A,A12)')
     &      '* DFT xc matrix             (CPU):',CPUTID
          ENDIF
          if (dirac_cfg_fde) then
             cputid = sectid(tdfde)
             WRITE(lupri,'(3X,A,A12)')
     &      '* FDE Fock matrix           (CPU):',CPUTID
          endif
        ENDIF
        CPUTID = SECTID(CTIDTOT)
        CALL GTINFO(DAYTID)
        WRITE(LUCYCL,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,DHF_INTTYP,DAYTID
        CALL FLSHFO(LUCYCL)
        WRITE(LUPRI, 1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,DHF_INTTYP,DAYTID
        if (doxml) then
           call xml_begin('scf-iteration')
           call xml_quantity('SCF energy',DHFERG,'Hartree')
           call xml_end('scf-iteration')
        endif
      ENDIF
C
      IF (TRIFCK) THEN
C
C       Get density matrices (normed to one)
C       and write them to DFDENS.
C
        CALL DENMAT(DMAT,CMO,IPRSCF)
        CALL WRIDNS(LUDENS,DMAT)
      END IF
C
      CALL FLSHFO(LUPRI)
      CALL MEMREL('PREDHF',WORK,KWORK,KWORK,KFREE,LFREE)
C
 666  CONTINUE
      CALL QEXIT('PREDHF')
      RETURN
 1010 FORMAT('It. ',I4,1P,G23.13,3D10.2,3X,A8,3X,A12,3X,A12,3X,A10)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DF_FMO */
      SUBROUTINE DF_FMO(FMO,FAO,DAO,TMAT,SAOMO,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Construct the Fock matrix in the reference MO basis
C
C     On input:
C    ===========
C       FAO   - Fock-matrix in AO basis
C       DAO   - Density matrix in AO-basis
C
C       FAO and DAO are generally destroyed !
C
C     On output:
C    ============
C       FMO    - Fock-matrix in reference MO basis
C       TMAT   - AO -> reference MO transformation matrix
C       SAOMO  - < AO | ref MO > overlap matrix
C
C     Called from: DFSOLV
C                  SCFCYC_1
C
C     Extracted from DFSOLV Feb. 2015 hjaaj
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      REAL*8    FMO(*),FAO(*),DAO(*),TMAT(*),SAOMO(*),WORK(*)
      LOGICAL   FNDLAB
C
      CALL QENTER('DF_FMO')
C
C
C     Transform AO Fock matrix to reference MO basis
C     ============================================== 
C
      CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','DF_FMO')
C         read AO2MO transformation matrix TMAT
      CALL READT(LUTMAT,N2TMT,TMAT)
C         read LUTMAT rec#2 with SAOMO
      CALL READT(LUTMAT,N2TMT,SAOMO)
      CLOSE(LUTMAT,STATUS='KEEP')

CMI ... perform all this only if we do not restart BSS from DC-SCF or DC-SCF from BSS-SCF
      IF (.NOT.RESTFCK) THEN
        IF (NASHT.GT.0.AND.AOC) THEN
          IF (TRIVEC) THEN
           CALL MEMGET2('REAL','CMOSAV',KCMOSAV,NCMOTQ,WORK,KFREE,LFREE)
           CALL MEMGET2('REAL','EIG   ',KEIG   ,NORBT ,WORK,KFREE,LFREE)
           CALL MEMGET2('INTE','IBEIG ',KIBEIG ,NORBT ,WORK,KFREE,LFREE)
           CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMOSAV),
     &                WORK(KEIG),WORK(KIBEIG),TOTERG,6)
          ELSE
           CALL MEMGET2('REAL','CMOSAV',KCMOSAV,0     ,WORK,KFREE,LFREE)
           CALL MEMGET2('REAL','EIG   ',KEIG   ,0     ,WORK,KFREE,LFREE)
           CALL MEMGET2('INTE','IBEIG ',KIBEIG ,0     ,WORK,KFREE,LFREE)
          END IF
           CALL MEMGET2('REAL','DTEMP',KDTEMP,N2BBASXQ*NFMAT,
     &                   WORK,KFREE,LFREE)
           CALL MEMGET2('REAL','DEVEC',KDEVEC,N2BBASXQ,
     &                   WORK,KFREE,LFREE)
           CALL MKMOFK(FMO,FAO,DAO,TMAT,WORK(KCMOSAV),
     &                 SAOMO,WORK(KDTEMP),WORK(KDEVEC),
     &                 WORK(KFREE),LFREE)
        ELSE
          CALL MKMOFK(FMO,FAO,DUMMY,TMAT,DUMMY,
     &                DUMMY,DUMMY,DUMMY,WORK(KFREE),LFREE)
        END IF
      ELSE ! IF  .NOT.RESTFCK  ELSE, i.e. RESTFCK true 

CMI ==================================================================
CMI  ... read the 2comp. Fock-Dirac matrix in Lowdin
CMI        MO basis (from pctra DC-SCF) to START 2c BSS-SCF ...
CMI  ... or the 4comp STARTING Fock-Dirac matrix (from rev.pctra BSS-SCF+H2cPP)
CMI     in RKB_MO basis to restart DC-SCF
CMI ==================================================================
       CALL OPNFIL(LUBSS,'BSSMAT','OLD','DF_FMO')
       IF (DO4C2C) THEN
C ...    read the FD_2c in Lowdin MO basis from the BSSMAT to 
C       provide starting for the postDHF-BSS-SCF procedure
        REWIND LUBSS
        IF (FNDLAB('FD2C_LMO',LUBSS)) THEN 
          ISIZE = (NESH(1)*NESH(1)*NZ) 
          IF (NFSYM.EQ.2) ISIZE=ISIZE+(NESH(2)*NESH(2)*NZ)
          CALL READT(LUBSS,ISIZE,FMO)
          IF (IPRHAM.GE.2) WRITE(LUPRI,'(/2X,A)') 
     &  'DF_FMO: FD_2c elements in Lowdin MO basis were read'// 
     &  ' from the BSSMAT file(FD2C_LMO)-starting guess for BSS-SCF.'
         ELSE
           CALL QUIT('DF_FMO: FD2C_LMO label in BSSMAT not found !') 
         END IF
      ELSE  ! RESTFCK=.false., continuing with four-component DC-SCF
         IF (FNDLAB('FD4C_RKB',LUBSS)) THEN
          ISIZE = (NORB(1)*NORB(1)*NZ) 
          IF (NFSYM.EQ.2) ISIZE=ISIZE+(NORB(2)*NORB(2)*NZ)
          CALL READT(LUBSS,ISIZE,FMO)
          IF (IPRHAM.GE.2) WRITE(LUPRI,'(/2X,A)')
     &    'DF_FMO: FD_4c elements in RKB MO basis were read'// 
     &    ' from the BSSMAT file(FD4C_RKB)-starting guess for DC-SCF.'
         ELSE
           CALL QUIT('DF_FMO: FD4C_RKB  not found in BSSMAT !') 
         END IF
       END IF
        CLOSE(LUBSS,STATUS='KEEP')
        RESTFCK = .FALSE.
      ENDIF 
C
C     Freeze orbitals by zeroing off-diagonal elements of FMO
C     =======================================================
C
      IF (DOMOFREEZE) THEN
         IF (.NOT.DIISAO) THEN
            CALL QUIT('.MOFREEZE Requires .DIISAO')
         ENDIF
         CALL FMO_FREEZE(FMO)
      ENDIF
      CALL QEXIT('DF_FMO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dfsolv */
      SUBROUTINE DFSOLV(FMO,EIG,IBEIG,FAO,DMAT,BMAT,ICHKCNV,
     &                  DOSEL,DOLEV,WORK,LWORK)
C*****************************************************************************
C
C     Solve the (Dirac)-Hartree-Fock eigenvalue equation
C
C     On input:
C    ===========
C       FAO   - Fock-matrix in AO-basis
C       IBEIG - super symmetry information
C       DMAT  - Density matrix in AO-basis
C
C     On output:
C    ============
C       FMO - MO-coefficients from current iteration
C       EIG - orbital eigenvalues
C       DMAT and FAO are generally destroyed !
C
C     Callef from: PREDHF
C                  SCFCYC_1
C
C     Written by T.Saue Mar 09 1995
C
C     Last revisions: Mar 09 1995
C                     Nov. 2005 - MI
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      PARAMETER ( DP5 = 0.50D00 , D0 = 0.00D00, D1 = 1.00D00,
     &            DM1 = -1.00D00, D2 = 2.00D00)
      LOGICAL DOSEL,DOLEV,FNDLAB
      DIMENSION FMO(*),EIG(*),IBEIG(*),FAO(*),DMAT(*),BMAT(*),
     &          WORK(LWORK)

      real(8) :: TIMSTR(2), TIMEND(2), TIMSTR_save(2)
      integer, allocatable :: IBEIG2(:,:)
C
#include "memint.h"

      CALL QENTER('DFSOLV')
      IF (IPRSCF .GE. 3) call TIMER2('START',TIMSTR,TIMEND)
      TIMSTR_save = TIMSTR
C
C
C     Transform AO Fock matrix to reference MO basis
C     ============================================== 
C
C     ..and modify open-shell diagonal blocks for open-shell calculations
C
      CALL MEMGET2('REAL','TMAT' ,KTMAT ,N2TMT,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','SAOMO',KSAOMO,N2TMT,WORK,KFREE,LFREE)
      CALL DF_FMO(FMO,FAO,DMAT,WORK(KTMAT),WORK(KSAOMO),
     &     WORK,KFREE,LFREE)
      IF (IPRSCF .GE. 3) call TIMER2('DF_FMO',TIMSTR,TIMEND)
C
C     Convergence control/acceleration
C     ================================
C
!     Level shift, PSB+MI/Aug2014: move the level shift before calling the
!     convergence acceleration because one needs the full matrix
!     for the convergence accelerations

      IF (DOLEV) CALL LEVSHI(FMO,FAO)

      IF(DOCCNV) THEN
        IF(DIISAO) THEN
          CALL MEMGET2('REAL','EVEC',KEVEC,N2BBASXQ      ,
     &       WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','FBUF',KFBUF,N2BBASXQ*NFMAT,
     &       WORK,KFREE,LFREE)
          CALL DHFCNV(FMO,FAO,DMAT,WORK(KEVEC),WORK(KFBUF),
     &                WORK(KTMAT),BMAT,ICHKCNV,WORK(KSAOMO),
     &                WORK(KFREE),LFREE)
        ELSE
          CALL DHFCNV(FMO,DUMMY,DUMMY,FAO,DMAT,DUMMY,BMAT,ICHKCNV,
     &                DUMMY,WORK(KFREE),LFREE)
        ENDIF
      ENDIF
C
C     If converged, then DHFEXIT is .true. now
C
      IF (DHFEXIT .AND. IPRSCF.GE.6) THEN
        DO IFSYM = 1,NFSYM
C        ... before extracting of EIG print the Fock matrix
           CALL HEADER(
     &    'DFSOLV: The converged Fock matrix, FMO, in MO basis ',-1)
           WRITE(LUPRI,'(2X,A,I1,A,I1)')
     &    'fermion symmetry: ',IFSYM,'/',NFSYM
           CALL PRQMAT(FMO(I2TMOT(IFSYM)+1),NTMO(IFSYM),NTMO(IFSYM),
     &          NTMO(IFSYM),NTMO(IFSYM),NZ,IPQTOQ(1,0),LUPRI)
           CALL FLSHFO(LUPRI)
        ENDDO
      ENDIF
      IF (IPRSCF .GE. 3) call TIMER2('DHFCNV',TIMSTR,TIMEND)
C
C     If converged, exit SCF here
C
      IF (DHFEXIT) GOTO 9999      
C
C     Diagonalize Fock matrix
C     =======================
C
      CALL GETTIM(CTIM1,WTIM1)
      CALL DFDIAG (FMO,EIG,IBEIG,FAO,.TRUE.,WORK,KFREE,LFREE)
      CALL GETTIM(CTIM2,WTIM2)
      TDDG = CTIM2-CTIM1
      WDDG = WTIM2-WTIM1
      IF (IPRSCF .GE. 3) call TIMER2('DFDIAG',TIMSTR,TIMEND)
C
C     sebastian+stefan: subtract level shift from orbital energies
C
      IF(DOLEV) CALL UNDOLEVSHI(EIG)
C
C     Overlap selection
C     =================
C
      IF(DOSEL) CALL SMOSEL(FAO,EIG,IBEIG,WORK(KFREE),LFREE)
C
C     Boson type selection
C     ====================
C
      IF(DOBOSSEL) CALL BOSSEL(EIG,IBEIG,FAO,WORK(KFREE),LFREE)
C
C     Automatic occupation
C     ====================
C
      IF(AUTOCC.OR.INIOCC) CALL DETOCC(EIG,WORK,KFREE,LFREE)      
C
C     
C     Backtransform all spinors (temporary storage in FMO)
C     ====================================================
C     NZT is the NZ value for the transformation matrix
C
      IF(NFROT.GT.0) CALL OPNFIL(LUTMAT,'DFFROZ','OLD','DFSOLV')
      DO 40 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 40
        IF(NFRO(I).EQ.0) THEN
          CALL BCKTR1(FMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &              FAO(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &              NTMO(I),NZ,NTMO(I),1,NFBAS(I,0),
     &              WORK(KTMAT+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,
     &              IPRSCF)
        ELSE
          CALL MEMGET2('INTE','IBUF ',KIND,NFRO(I),WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','BUF  ',KBUF,NFBAS(I,0)*NZ,
     &                 WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','EBUF ',KEBF,NTMO(I),WORK,KFREE,LFREE)
          CALL MEMGET2('INTE','IBBUF',KIBF,NTMO(I),WORK,KFREE,LFREE)
          CALL BCKTRA_FRZ(I,LUTMAT,FMO(ICMOQ(I)+1),FAO(I2TMOT(I)+1),
     &                EIG(IORB(I)+1),IBEIG(IORB(I)+1),
     &                WORK(KTMAT+I2TMT(I)),WORK(KIND),WORK(KBUF),
     &                WORK(KEBF),WORK(KIBF),IPRSCF)
          CALL MEMREL('DFSOLV.froz',WORK,KIND,KIND,KFREE,LFREE)
        ENDIF

        IF (IPRSCF.GE.3) THEN
          WRITE(LUPRI,'(/,2X,A)') 'DFSOLV: Eigenvalues:'
          WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,'/',NFSYM
          WRITE(LUPRI,'(I5,F25.10)') (J,EIG(IORB(I)+J),J=1,NTMO(I))
          IF(IPRSCF.GE.5) THEN
            WRITE(LUPRI,'(/,2X,A)') 
     &       'DFSOLV: Eigenvectors - MO coefficients'
            WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &       '*** Fermion corep ',I,'/',NFSYM
            CALL PRQMAT(FMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &           NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
        ENDIF
 40   CONTINUE
      IF (NFROT.GT.0) CLOSE(LUTMAT,STATUS='KEEP')
      IF (IPRSCF .GE. 3) call TIMER2('BCKTRA',TIMSTR,TIMEND)
C
C     Print eigenvalues 
C     (after backtransformation in case of frozen orbitals)
C
      IF (IPRSCF.GE.2) THEN
        CALL PREIGN(EIG,IBEIG,(IPRSCF .GE. 4),WORK(KFREE),LFREE)
      ELSE IF (IPRSCF.GE.0) THEN
        DO I = 1,NFSYM
           JHOMO = IORB(I)+NPSH(I)+NISH(I)
           JLUMO = JHOMO + NASH(I) + 1
           IF (NISH(I) .EQ. 0) JHOMO = JHOMO + 1 ! e.g. fsym 2, oxygen atom
           IF (NSSH(I) .EQ. 0) JLUMO = JLUMO - 1 ! very small test basis sets
           WRITE(LUPRI,'(A,I2,A,(30X,5(I5,F10.5)))')
     &     'E_HOMO...E_LUMO, symmetry',I,':',
     &     (J, EIG(J), J = JHOMO,JLUMO)
        END DO
      ENDIF
C
C     Write coefficients to file
C     ==========================
C
C     Write untransformed coefficients
      REWIND LUCMOS
      CALL WRITT(LUCMOS,N2TMOTQ,FAO)
C     Write transformed coefficients
      IF(ATOMIC) THEN
        CALL ATOMIC_WRICMO(LUCOEF,FMO,EIG,IBEIG,DHFERG)
      ELSE
        CALL WRICMO(LUCOEF,FMO,EIG,IBEIG,DHFERG,.FALSE.)
      ENDIF   
C
C     Dynamic update of overlap selection matrix, if not converged
C     (note that DMAT is destroyed)
C
      IF (DYNSEL) THEN
        CALL SMOGEN(FAO,DMAT)
      ENDIF
C
 9999 CONTINUE
      IF (IPRSCF .GE. 3) call TIMER2('DFSOLV rest ',TIMSTR,TIMEND)
      IF (IPRSCF .GE. 3) call TIMER2('DFSOLV total',TIMSTR_save,TIMEND)
      CALL MEMREL('DFSOLV',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL FLSHFO(LUPRI)
      CALL QEXIT('DFSOLV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dhfscf */
      SUBROUTINE DHFSCF(CMO,EIG,IBEIG,FOCK,WORK,LWORK)
C*****************************************************************************
C
C
C    Solve SCF
C
C
C   Called from PSISCF
C
C*****************************************************************************

      use interface_ao_specific

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "cbihr1.h"

      DIMENSION FOCK(*),EIG(*),IBEIG(*),CMO(*),
     &          WORK(LWORK)
C
C     Generate one-electron Fock matrix and write to file
C     ===================================================
C
!sk   two-component methods: pct-1-el fock matrix will be read from file (BSSMAT or X2CMAT)
      CALL ONEFCK(FOCK,IPRONE,WORK,LWORK)
!sk   we write 1-el fock matrix here and read it later during the scf cycle from file.

      CALL OPNFIL(LUFCK1,'DFFCK1','UNKNOWN','ONEFC1')
      CALL WRIFCK(LUFCK1,FOCK,1)
      CLOSE(LUFCK1,STATUS='KEEP')

!     write info to files
!     this is then used in dft and visualization
!     if you move this to somewhere else think about
!     geometry optimization!
      call interface_ao_write()

      IF(ONESYS) THEN
C
C     Solve one-electron system ...
C     =============================
C
        OPEN_FAC = 1.0D0
        CALL PREDHF(CMO,EIG,IBEIG,FOCK,DUMMY,DUMMY,IDUMMY,WORK,LWORK)
        DHFCONV(1) = .TRUE.
        DHFCONV(2) = .TRUE.
      ELSE
C
C     ... or perform SCF-cycle
C     ========================
C
        CALL SCFCYC(CMO,EIG,IBEIG,FOCK,WORK,LWORK)
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck scfcyc */
      SUBROUTINE SCFCYC(CMO,EIG,IBEIG,FOCK,WORK,LWORK)
C*****************************************************************************
C
C     Perform Dirac-Fock SCF-cycle
C
C     Written by T.Saue Mar 3 1995
C     Last revision: tsaue Mar 3 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION FOCK(*),EIG(*),IBEIG(*),CMO(*), WORK(LWORK)
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "cbihr2.h"
      CALL QENTER('SCFCYC')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET2('REAL','DMAT',KDMAT,N2BBASXQ*NFMAT,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','BMAT',KBMAT,MXDIIS*MXDIIS,WORK,KFREE,LFREE)

      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. , .TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      if (NPOS.GT.0) THEN
         CALL MEMGET2('INTE','NPOS',KPOS,NPOS,WORK,KFREE,LFREE)
      else
         KPOS = KFREE
      endif
C
C     Initiate/restart SCF-cycle
C
      CALL SCFCYC_1(CMO,EIG,IBEIG,FOCK,WORK(KDMAT),WORK(KBMAT),
     &            WORK(KPOS),WORK(KFREE),LFREE)
      CALL MEMREL('SCFCYC.SCFCYC_1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT ('SCFCYC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck SCFCYC_1 */
      SUBROUTINE SCFCYC_1(CMO,EIG,IBEIG,FOCK,DMAT,BMAT,NPOS,WORK,LWORK)
C*****************************************************************************
C
C     Perform Dirac-Fock SCF-cycle
C
C     Written by T.Saue Mar 3 1995
C     Last revision: tsaue Mar 3 1995
C
C*****************************************************************************

      use dirac_cfg
      use fde_mod
      use num_grid_gen
      use interface_mo_specific
      use xmlout
#ifdef HAS_PCMSOLVER
      use pcm_scf
#endif
#ifdef HAS_PELIB
      use pe_variables, only: peqm
#endif

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D0 = 0.0D0,DM1=-1.0D0,D1=1.0D0, DP5 = 0.50D00,
     &          D2 = 2.00D00)
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbreo.h"
      real(8), allocatable :: aoo2esssoc(:)
      LOGICAL   LBIT, tobe, DODIIS_save
      type(fde_import) :: itmp
      CHARACTER WALTID*12,CPUTID*12,SECTID*12,DAYTID*24
      DIMENSION FOCK(*),EIG(*),IBEIG(*),CMO(*),DMAT(*),
     &          BMAT(*),NPOS(*),WORK(LWORK)
      DIMENSION DVOVLP(MXOPEN)

#ifdef HAS_PCMSOLVER
      real(8)              :: time_pcm_fock
      character(12)        :: pcm_cputid
      real(8), allocatable :: fock_pcm(:, :, :)
#endif
#ifdef HAS_PELIB
      REAL(8) :: PE_TOT_NRG, PE_EL_NRG, PE_TIME
#endif
      CALL QENTER('SCFCYC_1')
#include "memint.h"
C
C     ***********************************
C     ***** Perform SCF-iterations ******
C     ***********************************
C
      DODIIS_save   = DODIIS
C
C     Initiate/restart SCF cycle
C
      DHFEXIT = .FALSE.
      CALL PREDHF(CMO,EIG,IBEIG,FOCK,DMAT,BMAT,NPOS,WORK,LWORK)
C     ... this is a restarted calculation we might already be converged
      IF(DHFEXIT) GOTO 30
C
C     Start cycle
C

      DO 10 ITER = 2,MAXITR
         NITER  = NITER + 1
         KITER  = ITER   ! transfer ITER via COMMON to SUBROUTINE DHFCNV
         CALL GTINFO(DAYTID)
         CALL GETTIM(CTIDSTR,WTIDSTR)
         WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +     '########## START ITERATION NO.',NITER,' ##########',
     &     DAYTID
         CALL FLSHFO(LUPRI)
C
C        Check what integral types contribute
C        ====================================
C
         if (ECPCALC) INTDEF = 1

         CALL INTCON(INTFLG,INTBUF,INTDEF,
     &               ERGVAL,CNVINT,NITER,ITRINT,DHF_INTTYP)
C
C        Construct density matrix
C        ========================
C
C        Get density matrices (normed to one)
C
         CALL DENMAT(DMAT,CMO,IPRSCF)
C
C       Print Mulliken gross populations, if requested
C
        IF (SCFPOP) CALL SCF_MULPOP(CMO,EIG,IBEIG,WORK(KFREE),LFREE)
#ifdef MOD_MATLAB_LOG
         call matexport_text('% In SCFCYC1')
         call matexport_double3('D',DMAT,NTBAS(0),NTBAS(0),
     &        NTBAS(0),NTBAS(0),NZ)
#endif
C
C        CMO not needed any more, may be reused between here and DFSOLV/hj
C
C        Construct two-electron Fock matrix
C        (stored in FOCK)
C        ==================================
C
         TDF2 = SECOND()
         IF(DODSCF.AND.INTFLG.EQ.INTBUF) THEN
C           *** Differential density matrix code ***
C           read old density matrix in FOCK before saving new density matrix/hj
            CALL READNS(LUDENS,FOCK)
            CALL WRIDNS(LUDENS,DMAT)
            IF (NISHT .GT. 0) THEN
               IF(FIXDIF) THEN
                  DCOVLP = D1
               ELSE
                  DCOVLP = DDOT(N2BBASXQ,FOCK,1,FOCK,1)
                  DCOVLP = DDOT(N2BBASXQ,DMAT,1,FOCK,1) / DCOVLP
               END IF
               WRITE(LUPRI,'(I5,A,F7.4)') NITER,
     &            ' *** Differential density matrix. DCOVLP     =',
     &            DCOVLP
               IF (DCOVLP .LT. D0 .OR. DCOVLP .GT. D2) THEN
                  WRITE(LUPRI,'(A)') ' INFO: Useless value for DCOVLP,'
     &               //' value reset to: 1.0D0'
                  DCOVLP = D1
               END IF
               CALL DAXPY(N2BBASXQ,(-DCOVLP),FOCK,1,DMAT,1)
            END IF
            IF (AOC) THEN
               DO IOPEN = 1,NOPEN
                  IF (FIXDIF) THEN
                     DVOVLP(IOPEN) = D1
                  ELSE
                     DVOVLP(IOPEN) =
     &                  DDOT(N2BBASXQ,FOCK(1+N2BBASXQ*IOPEN),1,
     &                                FOCK(1+N2BBASXQ*IOPEN),1)
                     DVOVLP(IOPEN) =
     &                  DDOT(N2BBASXQ,DMAT(1+N2BBASXQ*IOPEN),1,
     &                  FOCK(1+N2BBASXQ*IOPEN),1) / DVOVLP(IOPEN)
                  END IF
                  WRITE(LUPRI,'(I5,A,I2,A,F7.4)') NITER,
     &               ' *** Differential density matrix. DVOVLP(',
     &               IOPEN,') =',DVOVLP(IOPEN)
                  IF (DVOVLP(IOPEN) .LT. D0
     &                  .OR. DVOVLP(IOPEN) .GT. D2) THEN
                     WRITE(LUPRI,'(A,I2,A)')
     &                  ' INFO: Useless value for DVOVLP(',IOPEN,'),'
     &                  //' value reset to: 1.0D0'
                     DVOVLP(IOPEN) = D1
                  END IF
                  CALL DAXPY(N2BBASXQ,(-DVOVLP(IOPEN)),
     &               FOCK(1+N2BBASXQ*IOPEN),1,DMAT(1+N2BBASXQ*IOPEN),1)
               END DO
            END IF
C
C           Note that all Fock matrices are multiplied with 2
C              in TWOFCK.
C
            CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,DMAT,NFMAT,
     &                  NPOS,INTFLG,IPRSCF,WORK(KFREE),LFREE)
C           read old Fock matrix in DMAT to finish new Fock matrix/hj
            CALL REAFCK(LUFCK2,DMAT,.FALSE.,NFMAT)
C
            IF (NISHT .GT. 0)
     &         CALL DAXPY(N2BBASXQ,DCOVLP,DMAT,1,FOCK,1)
            IF (AOC) THEN
               DO IOPEN = 1,NOPEN
                  CALL DAXPY(N2BBASXQ,DVOVLP(IOPEN),
     &               DMAT(1+N2BBASXQ*IOPEN),1,
     &               FOCK(1+N2BBASXQ*IOPEN),1)
               END DO
            END IF
C           restore density matrix in DMAT/hj
            CALL READNS(LUDENS,DMAT)
         ELSE
            CALL WRIDNS(LUDENS,DMAT)
            CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,DMAT,NFMAT,
     &                  NPOS,INTFLG,IPRSCF,WORK(KFREE),LFREE)
         ENDIF
C
         IF(DOTSC) CALL TSCORR(FOCK,WORK(KFREE),LFREE,IPRSCF)
         TDF2 = SECOND() - TDF2
C
C        Write current two-electron Fock matrix to file
C        Note that possible DFT/PCM contributions are not 
C        included since DFFCK2 is used for completing two-electron
C        Fock matrices generated by differential  densities
C        =======================================================
C
         CALL WRIFCK(LUFCK2,FOCK,NFMAT)

        if (dirac_cfg_fde) then
           call fde_get_import_info(itmp)
           if (itmp%im_update_vemb) then
              tdfde= second()
              call interface_mo_write()
              call fde_dirac_set_nz(nz)
#ifdef VAR_MPI
           if (parcal) call dirac_parctl(FDE_PAR)
#endif
              call fde_calculate_emb_pot_mat(ntbas(0),dmat,fock)
              tdfde = second() - tdfde
           endif
        endif
!       xc contribution
        if (dirac_cfg_dft_calculation) then
           tdxc = second()
           call interface_mo_write()
           call generate_num_grid(dmat)
#ifdef VAR_MPI
           if (parcal) call dirac_parctl( XCINT_PAR )
#endif
           call xcint_potential_rks(ntbas(0),
     &                              dmat,
     &                              fock)
           tdxc = second() - tdxc
        end if

C
C        Solvent contribution
C        ====================
C
         IF(SOLVEN) THEN
           WRITE(LUPRI,'(A)') '* Calling SOLFCK'
           NF = 1
           CALL SOLFCK(FOCK,DMAT,NF,ESOLVE,ESOLVN,
     &                 WORK(KFREE),LFREE,IPRSOL)
         ENDIF
C
C        Polarizable embedding contribution
C        ==================================
C
#ifdef HAS_PELIB
         IF (PEQM) THEN
            IF (IPRSCF >= 1) WRITE(LUPRI,'(A)') '* Calling PE library'
            PE_TIME = SECOND()
            CALL PELIB_IFC_FOCK(DMAT, FOCK, PE_TOT_NRG, PE_EL_NRG)
            PE_TIME = SECOND() - PE_TIME
         END IF
#endif
C
C        Read one-electron Fock matrix from file into CMO
C        ================================================
C
         CALL REAFCK(LUFCK1,CMO,.FALSE.,1)
C
C        Calculate correction according to the One-center app. models 
C        ============================================================
C
         IF (ONECAP.AND. .NOT.SMLV1C) THEN
              CALL CATCORR(CMO,FOCK,DMAT,WORK(KFREE),LFREE,IPRSCF)
         END IF
C
C
C        Calculate total energy
C        ======================
         inquire(file='aoo2e-ss-soc-contributions',exist=tobe)
         if(tobe)then
           open(99,file='aoo2e-ss-soc-contributions',status='old',
     &     form='unformatted',
     &     access='sequential',action="readwrite",position='rewind')
           allocate(aoo2esssoc(ntbas(0)**2*nz))
           read(99) aoo2esssoc(1:ntbas(0)**2*nz)
           close(99,status='keep')
           call daxpy(ntbas(0)**2*nz,-1.0d0,aoo2esssoc,1,cmo,1)
           call daxpy(ntbas(0)**2*nz, 1.0d0,aoo2esssoc,1,fock,1)
!          if(aoc) call quit('AOC not working yet for aoo2esssoc...')
           !> stefan: enabling this is wrong but not doing it as well...
           !> if possible the contributions should be added in the MO basis but how to calculate the energy at this point correctly???
!          IF (AOC) THEN
!             DO IOPEN = 1,NOPEN
!                CALL DAXPY(N2BBASXQ,1.0d0,
!    &                      aoo2esssoc,1,
!    &                      FOCK(1+N2BBASXQ*IOPEN),1)
!              END DO
!           END IF
           deallocate(aoo2esssoc)
         end if
C
         ! here F1 is in CMO and F2 is in FOCK
         CALL ERGCAL(CMO,FOCK,DMAT,WORK(KFREE),LFREE)
#ifdef HAS_PCMSOLVER
         if (dirac_cfg_pcm) then
           write(lupri,'(a)') '* Calling pcm_oper_ao_driver'
           time_pcm_fock = second()
           allocate(fock_pcm(ntbas(0), ntbas(0), nz))
           fock_pcm = 0.0d0
! fock is the 2-el Fock matrix
           call pcm_oper_ao_driver(fock_pcm, 'TotASC'//char(0), 
     &          work(kfree), lfree)
           call daxpy(n2bbasx, -1.0d0, fock_pcm, 1, fock, 1)
!     print *, "FOCK_PCM"
!     call prqmat(fock_pcm, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
!                      nz, ipqtoq(1,0), 6)
           deallocate(fock_pcm)
           time_pcm_fock = second() - time_pcm_fock
         endif 
#endif         

C
C        Add one-electron Fock matrix to two-electron matrix
C        ===================================================
C
         CALL DAXPY(N2BBASXQ,D1,CMO,1,FOCK,1)
C
C        Add active two-electron Fock matrices
C
         IF(AOC) THEN
           DO IOPEN = 1,NOPEN
C
C             The factor is only DF and not D2*DF because of the
C                multiplication with D2 in TWOFCK.
C
              CALL DAXPY(N2BBASXQ,DF(IOPEN),
     &                 FOCK(1+N2BBASXQ*IOPEN),1,FOCK,1)
           END DO
         ENDIF
C
C       Write total Fock matrix
C
        CALL WRIFCK(LUFCKT,FOCK,NFMAT)
C
         IF(IPRSCF.GE.5) THEN
           IF (NISHT .GT. 0) THEN
            CALL HEADER(
     &      'SCFCYC_1: Total FD matrix in SA-AO basis',-1)
            CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NZ,IPQTOQ(1,0),LUPRI)
           END IF
           IF(AOC) THEN
            DO IOPEN = 1,NOPEN
             CALL HEADER(
     &       'SCFCYC_1: Total FV matrix in SA-AO basis',-1)
             WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &       'matrix number ',IOPEN,'/',NOPEN
             CALL PRQMAT(FOCK(1+N2BBASXQ*IOPEN),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
            END DO
           ENDIF
         END IF
C
C        Solve DF eigenvalue problem
C        ===========================
C        Input:  FOCK - Total Fock matrix on AO-basis
C                DMAT - Density matrix
C                BMAT - DIIS matrix
C        Output: CMO  - MO-coefficients for current iteration
C
         IF(ATHUCK.AND.(NITER.EQ.2)) THEN
           ICHKCNV = 0
         ELSE
           ICHKCNV = 1
         ENDIF

         CALL DFSOLV(CMO,EIG,IBEIG,FOCK,DMAT,BMAT,ICHKCNV,
     &               OVLSEL,DOLEVEL,WORK(KFREE),LFREE)

        IF(IPRSCF.GE.1) THEN
          CPUTID = SECTID(TDF2)
          WRITE(LUPRI,'(3X,A,A12)')
     &      '* Two-electron Fock matrix  (CPU):',CPUTID
#ifdef HAS_PCMSOLVER
         if (dirac_cfg_pcm) then
           pcm_cputid = sectid(time_pcm_fock) 
           write(lupri,'(3X,A,A12)')
     &      '* PCM Fock matrix contribution  (CPU):', pcm_cputid
         endif 
#endif
#ifdef HAS_PELIB
          IF (PEQM) THEN
            CPUTID = SECTID(PE_TIME)
            WRITE(LUPRI,'(3X,A,A12)')
     &      '* PE Fock matrix            (CPU):',CPUTID
          END IF
#endif
          CPUTID = SECTID(TDDG)
          WALTID = SECTID(WDDG)
          IF(NZ.EQ.1) WRITE(LUPRI,'(3X,4A)')
     &      '* Real diagonalization      (CPU): ',CPUTID,
     &      '  (WALL): ',WALTID
          IF(NZ.EQ.2) WRITE(LUPRI,'(3X,4A)')
     &      '* Complex diagonalization   (CPU): ',CPUTID,
     &      '  (WALL): ',WALTID
          IF(NZ.EQ.4) WRITE(LUPRI,'(3X,4A)')
     &      '* Quaternion diagonalization(CPU): ',CPUTID,
     &      '  (WALL): ',WALTID
        ENDIF
        CALL GETTIM(CTIDEND,WTIDEND)
        CTIDTOT = CTIDEND - CTIDSTR
        WTIDTOT = WTIDEND - WTIDSTR
        ITRSCF(INTFLG) = ITRSCF(INTFLG) + 1
        SCFTID(INTFLG) = SCFTID(INTFLG) + WTIDTOT
        WALTID = SECTID(WTIDTOT)
        CPUTID = SECTID(CTIDTOT)
        WRITE(LUPRI,'(4A)')
     &    '>>> Total wall time: ',WALTID,
     &    ', and total CPU time : ',CPUTID
        CALL GTINFO(DAYTID)
        WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +    '########## END ITERATION NO.',NITER,' ##########',
     &    DAYTID
        WRITE(LUCYCL,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,WALTID,DHF_INTTYP,DAYTID(1:10)
        CALL FLSHFO(LUCYCL)
        WRITE(LUPRI, 1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,WALTID,DHF_INTTYP,DAYTID(1:10)
        if (doxml) then
           call xml_begin('scf-iteration')
           call xml_quantity('SCF energy',DHFERG,'Hartree')
           call xml_end('scf-iteration')
        endif
        CALL FLSHFO(LUPRI)
C
C       switch off onecap approximation ?
        IF (ONECOFF) THEN
           ONECOFF = .FALSE.
C          do continue even if it was converged with ONECAP :
           DHFEXIT = .FALSE.
C          reset intbuf so no old Fock matrices are used :
           INTBUF = 0
           WRITE(LUPRI,'(/A,1P,D15.6/2A/)') 
     &  ' ** INFO : ONECAP. Wave function converged beyond : ONECNV = ',
     &     ONECNV,'          ONECAP will be turned off, and integrals',
     &     ' evaluated on the fly !'
        ELSE
           INTBUF = INTFLG
        END IF
C
        IF(DHFEXIT) GOTO 30
 10   CONTINUE    ! DO 10 ITER = 2,MAXITR

      WRITE(LUPRI,'(/A)')
     &' ** Exit SCF because maximum number of iterations reached.'
C
C     **********************************
C     ***** End of SCF-iterations ******
C     **********************************
C
 30   CONTINUE

      CALL FLSHFO(LUPRI)

      DODIIS   = DODIIS_save
C
C     Close necessary files
C
      CLOSE(LUCOEF,STATUS='KEEP')
      CLOSE(LUCMOS,STATUS='KEEP')
      CLOSE(LUDENS,STATUS='KEEP')
      CLOSE(LUFOCK,STATUS='KEEP')
      CLOSE(LUFCK1,STATUS='KEEP')
      CLOSE(LUFCK2,STATUS='KEEP')
      CLOSE(LUFCKT,STATUS='KEEP')
      CLOSE(LUCYCL,STATUS='KEEP')
      CLOSE(LUEVEC,STATUS='KEEP')

      CALL QEXIT('SCFCYC_1')
      RETURN
 1010 FORMAT('It. ',I4,1P,G23.13,3D10.2,3X,A8,3X,A12,3X,A12,3X,A10)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ergcal */
      SUBROUTINE ERGCAL(F1,F2,DMAT,WORK,LWORK)
C****************************************************************************
C
C     PURPOSE:
C   ==============
C        Calculate the SCF energy using the formula
C
C           E(SCF) = D(IJ)(*)(F1(IJ)+0.5F2(IJ))
C
C     On input: 
C        F1, F2 - one- and two-electron Fock matrixes in SA-AO basis
C        DMAT  - density matrix in SA-AO basis
C
C     On output: internal variables of SCF energy updated
C 
C
C     Last revision: 971022-jth
C                    Sept 2005, Miro Ilias - added control outputs
C
C*****************************************************************************

      use dirac_cfg
      use fde_mod
#ifdef HAS_PCMSOLVER
      use pcm_scf
#endif
#ifdef HAS_PELIB
      use pe_variables, only: peqm
#endif

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D2=2.0D0, DP5 = 0.5D00,
     &          DP25 = 0.25D00)
C
#include "frame.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbpsi.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbgen.h"
      LOGICAL LVCON
      CHARACTER M(4)*2,FMT*6,MXFORM*6
      DIMENSION F1(N2BBASX,NZ),F2(N2BBASX,NZ,NFMAT),
     &          DMAT(N2BBASX,NZ,NFMAT),
     &          E(4,2,0:MXOPEN),WORK(LWORK)
      DATA M/'AR','AI','BR','BI'/
c     besterg is the lowest energy we encountered so far
      DATA INUM /0/
      SAVE INUM
#ifdef HAS_PCMSOLVER
      real(8)      :: pcm_energy
#endif
      REAL(8) :: PE_TOT_NRG, PE_EL_NRG
      REAL(8), DIMENSION(:), ALLOCATABLE :: PE_FMAT
      type(fde_import) :: itmp
      data besterg /1.0D100/
      save besterg
C
      CALL QENTER('ERGCAL')
      INUM = INUM + 1
C
      LVCON = DOLVC.AND.
     &        ( (INTFLG.EQ.3 .or. INTFLG.EQ.11) .OR. 
     &          (ONECAP.AND.INTV1C.EQ.2) )
      CALL DZERO(E,8*(MXOPEN+1))
      ERGBUF  = DHFERG
      E1PART  = D0
      E2PART  = D0
C

      if (dirac_cfg_fde) then
         call fde_get_import_info(itmp)
      end if


      IF(AOC) THEN
        ITOP = NOPEN
      ELSE
        ITOP = 0
      ENDIF
      DO ISHELL = 0,ITOP
C
C        Factor D2 comes from Time-Reversal Symmetry.
C
         DFRAC1 = D2*DF(ISHELL)
C
C        Factor DP25 = DP5*DP5 comes from the formula and the
C           fact that all Fock matrices were multiplied with 2
C           in TWOFCK.
C
         DFRAC2 = DP25*DFRAC1*DFRAC1*DA(ISHELL)
C...     Intra-shell contributions
         DO IZ = 1,NZ
            E(IZ,1,ISHELL) =
     &         DFRAC1*DDOT(N2BBASX,DMAT(1,IZ,ISHELL+1),1,F1(1,IZ),1)
            E(IZ,2,ISHELL) =
     &         DFRAC2*DDOT(N2BBASX,DMAT(1,IZ,ISHELL+1),1,
     &                               F2(1,IZ,ISHELL+1),1)
         END DO
C...     Inter-shell contribution
         DO JSHELL = 0,ISHELL-1
C
C           The factor is 0.5 * (2*f_{o1})*(2*f_{o2})
C           Divide by 2 (from TWOFCK)
C           Total factor is DF(ISHELL)*DF(JSHELL)
C
            DFRAC1 = DF(ISHELL)*DF(JSHELL)
            DO IZ = 1,NZ
               ETEMP =
     &            DFRAC1*DDOT(N2BBASX,DMAT(1,IZ,JSHELL+1),1,
     &                                  F2(1,IZ,ISHELL+1),1)
               E(IZ,2,ISHELL) = E(IZ,2,ISHELL) + ETEMP
               E(IZ,2,JSHELL) = E(IZ,2,JSHELL) + ETEMP
            END DO
         END DO
      END DO
!     correct for contribution from ks potential to two-electron energy
      if (dirac_cfg_dft_calculation) then
!Miro: separate this contribution into own variable, xcint_xc_energy
!     the density for calculation is the same as for xcint_xc_contribution to the 2el Fock matrix
        xcint_xc_energy=xcint_get_xc_mat_energy()
        e(1, 2, 0) = e(1, 2, 0) - xcint_xc_energy
        !print *
        !print *,'ERGCAL: xcint_get_xc_mat_energy()=',xcint_xc_energy
        !print *
      end if
      DO IZ = 1,NZ
        DO IOPEN = 0,ITOP
           E1PART = E1PART + E(IZ,1,IOPEN)
           E2PART = E2PART + E(IZ,2,IOPEN)
        END DO
      END DO
C
      DHFERG = POTNUC
      IF(SOLVEN) THEN
        E2PART = E2PART - ESOLVE
        ERGSOL = ESOLVE + ESOLVN
        DHFERG = DHFERG + ERGSOL
      ENDIF

#ifdef HAS_PCMSOLVER
      if (dirac_cfg_pcm) then
! Add polarization energy to dhferg 
        call pcm_energy_driver(dmat, pcm_energy, work, lwork)
        dhferg = dhferg + pcm_energy
      endif 
#endif
#ifdef HAS_PELIB
      IF (PEQM) THEN
        ! The PE Fock matrix is added to F2 before entering ERGCAL
        ! and the electronic PE energy is therefore added to E2PART.
        ! Therefore we subtract its' contribution from the total PE energy.
        ALLOCATE(PE_FMAT(N2BBASX))
        E_PE = 0.0D0
        CALL PELIB_IFC_FOCK(DMAT, PE_FMAT, PE_TOT_NRG, PE_EL_NRG)
        IF (IPRSCF >= 1) CALL PELIB_IFC_ENERGY(DMAT)
        DHFERG = DHFERG + PE_TOT_NRG - 0.5d0 * PE_EL_NRG
        E_PE = PE_TOT_NRG 
        DEALLOCATE(PE_FMAT)
      END IF
#endif
      IF(LVCON) THEN
         IF (ONECAP .AND. (INTV1C.EQ.2)) THEN
            JPRINT=0
         ELSE
            JPRINT=IPRSCF
         END IF
         CALL LVCORR(DMAT,WORK,LWORK,JPRINT)
         DHFERG = DHFERG + CORRLV
      ENDIF
      ELERGY = E1PART + E2PART

      if (dirac_cfg_dft_calculation) then
        xcint_elergy=xcint_get_xc_energy()
        !elergy = elergy + xcint_get_xc_energy()
        elergy = elergy + xcint_elergy
        IF(IPRSCF.GE.5) THEN
          print *,'ERGCAL: xcint_get_xc_energy()=',xcint_elergy
        ENDIF
      end if

      if (dirac_cfg_fde) then
         call fde_get_import_info(itmp)
         if (itmp%im_update_vemb) then
             elergy = elergy - fde_get_embpot_energy()
         else if (itmp%im_vemb) then
             elergy = elergy - fde_dirac_embpot_static_energy()
         end if
!        deallocate(vemb_mat)
      end if

      DHFERG = DHFERG + ELERGY
      
#ifdef HAS_PELIB
      ! edh: Correct also elergy for printing...
      if (peqm) then 
        elergy = elergy - 0.5d0 * pe_el_nrg
      end if
#endif
C
      IF(IPRSCF.GE.1) THEN
         CALL HEADER('Output from ERGCAL',-1)
CMI      ... print out matrixes ...
         IF (IPRSCF.GE.5) THEN
          CALL HEADER(
     &    'ERGCAL: Entering total Fock F1 matrix in SA-AO',-1)
          CALL PRQMAT(F1,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &               IPQTOQ(1,0),LUPRI)
          DO IFMAT=1,NFMAT
          CALL HEADER(
     &    'ERGCAL: Entering total Fock F2 matrix in SA-AO',-1)
          WRITE(LUPRI,'(2X,A,I2,A,I2)')
     &    'Fock matrix no: ',IFMAT,'/',NFMAT
          CALL PRQMAT(F2(1,1,IFMAT),NTBAS(0),NTBAS(0),
     &           NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
          CALL HEADER(
     &    'ERGCAL: Entering total density matrix in SA-AO',-1)
          WRITE(LUPRI,'(2X,A,I2,A,I2)')
     &    'dens.matrix no: ',IFMAT,'/',NFMAT
          CALL PRQMAT(DMAT(1,1,IFMAT),NTBAS(0),NTBAS(0),
     &           NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
          ENDDO
         ENDIF

         ERGMAX = MAX(ABS(E1PART),ABS(E2PART))*10
         FMT = MXFORM(ERGMAX,20)
            WRITE(LUPRI,'(A,'//FMT//')')
     +       'Total electronic energy            :  ',ELERGY,
     +       'Nuclear potential energy           :  ',POTNUC
         IF(LVCON) WRITE(LUPRI,'(A,'//FMT//')')
     +       'SS Coulombic correction            :  ',CORRLV
         IF(SOLVEN) WRITE(LUPRI,'(A,'//FMT//')')
     +       'Solvent energy                     :  ',ERGSOL
#ifdef HAS_PCMSOLVER
         IF(dirac_cfg_pcm) WRITE(LUPRI,'(A,'//FMT//')')
     +       'Polarization energy                :  ',get_pcm_energy()
#endif
        ! later this will be replaced by a get energy routine
#ifdef HAS_PELIB
        IF (PEQM) WRITE(LUPRI,'(A,'//FMT//')')
     +       'Embedding energy                   :  ',PE_TOT_NRG
#endif
         WRITE(LUPRI,'(A,'//FMT//')')
     +       'Total energy                       :  ',DHFERG
         IF(IPRSCF.GE.2) THEN
            WRITE(LUPRI,'(6X,A,'//FMT//')')
     &         '* 1-electron energy          :  ',E1PART
            IF(NZ.GT.1) THEN
               DO 20 IZ = 1,NZ
               DO IOPEN = 0,ITOP
                  WRITE(LUPRI,'(6X,A,A2,A,I2,A,'//FMT//')')
     +               '            from ',M(IZ),'(Shell ',
     +               IOPEN,'):  ',E(IZ,1,IOPEN)
               END DO
   20          CONTINUE
            ENDIF
            WRITE(LUPRI,'(6X,A,'//FMT//')')
     +         '* 2-electron energy          :  ',E2PART
            IF(NZ.GT.1) THEN
               DO 30 IZ = 1,NZ
               DO IOPEN = 0,ITOP
                  WRITE(LUPRI,'(6X,A,A2,A,I2,A,'//FMT//')')
     +               '            from ',M(IZ),'(Shell ',
     +               IOPEN,'):  ',E(IZ,2,IOPEN)
               END DO
   30          CONTINUE
            ENDIF
         END IF
      ENDIF
c     Check if this is the lowest energy so far, in that case
c     save the DFCOEF file as DFCOEF.BEST
      if (dhferg.lt.besterg) then
#if defined SYS_WINDOWS
         call system('copy DFCOEF DFCOEF.BEST')
#else
         call system('cp DFCOEF DFCOEF.BEST')
#endif
         besterg = dhferg
      endif
      CALL QEXIT('ERGCAL')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck intcon */
      SUBROUTINE INTCON(INTFLG,INTOLD,INTDEF,
     &                  CNV,CNVINT,ITR,ITRINT,INTTYP)
C***********************************************************************
C
C     Check what integrals contribute in a given SCF-iteration
C       CNVINT has precedence over ITRINT !!!
C       The information is bitpacked in INTDEF:
C         0001   LL-contribution
C         0010   SL-contribution
C         0100   SS-contribution
C         1000   GT-contribution
C       and all combinations.....
C
C     Written by T.Saue Mar 09 1995
C     Last revision: Aug 24 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      LOGICAL LBIT
      CHARACTER INTTYP*12
      DIMENSION CNVINT(2),ITRINT(2)
      INTTYP = '            '
C
C     LL-integrals
C     ============
C
      IF(LBIT(INTDEF,1)) THEN
        ILLINT = 1
      ELSE
        ILLINT = 0
      ENDIF
C
C     SL-integrals
C     ============
C
      IF(LBIT(INTDEF,2)) THEN
        ISLINT = 1
        IF(.NOT.LBIT(INTOLD,2)) THEN
          IF(ABS(CNV) == 0.0d0 .and. CNVINT(1) /= dummy) ISLINT = 0
          IF(ABS(CNV).GT.CNVINT(1)) ISLINT = 0
          IF(ITR     .LT.ITRINT(1)) ISLINT = 0
        ENDIF
      ELSE
        ISLINT = 0
      ENDIF
C
C     SS-integrals
C     ============
C
      IF(LBIT(INTDEF,3)) THEN
        ISSINT = 1
        IF(.NOT.LBIT(INTOLD,3)) THEN
          IF(ABS(CNV) == 0.0d0 .and. CNVINT(2) /= dummy) ISSINT = 0
          IF(ABS(CNV).GT.CNVINT(2)) ISSINT = 0
          IF(ITR     .LT.ITRINT(2)) ISSINT = 0
        ENDIF
      ELSE
        ISSINT = 0
      ENDIF
C
C     Gaunt-integrals
C     ===============
C
      IF(LBIT(INTDEF,4)) THEN
        IGTINT = 1
        IF(.NOT.LBIT(INTOLD,4)) THEN
C         Activate at the same point as the SS for the moment
C         One could als put it with SL, comes in that order
          IF(ABS(CNV) == 0.0d0 .and. CNVINT(2) /= dummy) IGTINT = 0
          IF(ABS(CNV).GT.CNVINT(2)) IGTINT = 0
          IF(ITR     .LT.ITRINT(2)) IGTINT = 0
        ENDIF
      ELSE
        IGTINT = 0
      ENDIF
C
      IF(ILLINT.EQ.1) INTTYP(1:3) = 'LL '
      IF(ISLINT.EQ.1) INTTYP(4:6) = 'SL '
      IF(ISSINT.EQ.1) INTTYP(7:9) = 'SS '
      IF(IGTINT.EQ.1) INTTYP(10:12) = 'GT '
      INTFLG= ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bossel */
      SUBROUTINE BOSSEL(EIG,IBEIG,CFU,WORK,LWORK)
C***********************************************************************
C
C     Select occupied orbitals based on (approximate) boson symmetry.
C     Input is coefficients in MO-basis.
C
C     Written by H.J.Aa.Jensen Oct. 2004
C     Modified by A. Sunaga Oct. 2018 for atomic symmetry
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION EIG(*),IBEIG(*),CFU(*),WORK(LWORK)
C     
#include "dcbdhf.h"
#include "dcborb.h"
#include "dgroup.h"
C
      INTEGER IBUF(MAX_BOS_BL)
      INTEGER JOCC(2), KOCC(2)
C
      CALL QENTER('BOSSEL')
C
      CALL ICOPY(MAX_BOS_BL,NISH_BOS,1,IBUF,1)
      DO I = 1,NFSYM
C     JOCC(I) = NPSH(I) + NISH(I) + NACSH(I,1:IOPEN-1)
        KOCC(I) = NPSH(I)
        JOCC(I) = NPSH(I)
        IF(NISH_DHF(I).GT.0) THEN
          KOCC(I) = KOCC(I) + NISH_DHF(I)
          CALL BOSSE2(I,EIG(IORB(I)+1),IBEIG(IORB(I)+1),
     &         CFU(I2ORBT(I)+1),NORB(I),JOCC(I),IBUF)
          IF(JOCC(I).NE.KOCC(I)) THEN
            WRITE(LUPRI,1000)
            CALL QUIT('Check the number of closed-shell electrons !')
          ENDIF
        ENDIF
      ENDDO
C     
      DO IOPEN = 1,NOPEN
        DO I = 1,NFSYM
          IF(NACSH(I,IOPEN).GT.0) THEN
            CALL ICOPY(MAX_BOS_BL,NACSH_BOS(1,IOPEN),1,IBUF,1)
            KOCC(I) = KOCC(I) + NACSH(I,IOPEN)
            CALL BOSSE2(I,EIG(IORB(I)+1),IBEIG(IORB(I)+1),
     &                   CFU(I2ORBT(I)+1),NORB(I),JOCC(I),IBUF)
            IF(JOCC(I).NE.KOCC(I)) THEN
              WRITE(LUPRI,1000)
              write(*,*)"(JOCC(I),KOCC(I)",JOCC(I),KOCC(I)
              CALL QUIT('Check the number of open-shell electrons !')
            ENDIF            
          ENDIF  
        ENDDO
      ENDDO
      CALL QEXIT('BOSSEL')
      RETURN
 1000 FORMAT
     & (/'Please check the below points in input',  
     & /'e.g. The basis set is too small ', 
     &   'compared with the occupation number',
     & /'e.g. The sign of kappa is wrong: 2 is for d_3/2, while ', 
     &   '-2 is for p_3/2')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bosse2 */
      SUBROUTINE BOSSE2(IFSYM,EIG,IBEIG,CFU,NO,JOCC,NB_CNT)
C***********************************************************************
C
C     Select occupied orbitals based on (approximate) boson symmetry.
C     Input is coefficients in MO-basis in a fermion symmetry IFSYM.
C
C     Written by H.J.Aa.Jensen Oct. 2004
C     Modified by A. Sunaga Oct. 2018 for atomic symmetry
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcborb.h"
      DIMENSION EIG(*), CFU(NO,NO,NZ)
      INTEGER   IBEIG(*), NB_CNT(*)

      CALL QENTER('BOSSE2')
C
C hjaaj Oct 2004: must also implement open shells (NASHB*)
      JST = JOCC + 1
      IF(IFSYM.EQ.2)KPMAX = ( -1 + int(sqrt(dble(1 + 8*N_SUB_BL(1)))))/2   
      DO J = JST,NO
         IF (ATOMIC) THEN
            IF(NFSYM.EQ.2)THEN
              JSYM = abs(IBEIG(J))
              IF (IFSYM.EQ.2)THEN
                JSYM  = JSYM + KPMAX*(KPMAX+1)/2
              ENDIF
            ELSE
              INDEX1=INT(SQRT(dble(2*ABS(IBEIG(J)))+0.25D0)+0.4999D0)
              IF(IBEIG(J).LT.0)THEN
                JADD = INDEX1*(INDEX1-1)/2               
              ELSE
                JADD = INDEX1*(INDEX1+1)/2
              ENDIF              
              JSYM = ABS(IBEIG(J)) + JADD
            ENDIF  
C           ... IBEIG = ID. ID = ID_SUB_BL(ISUB,IFRP) in dirscf.F/INIBOS.
         ELSEIF (LINEAR) THEN
            JSYM = (abs(IBEIG(J)) + 1) / 2
C           ... IBEIG = 2 * abs(MJ)
C           IF (IFSYM.EQ.2) JSYM = JSYM + N_SUB_BL(1)
Cayaki      ^ It is no longer neccessary, because Atomic symmetry is used. 
         ELSE
            JSYM = IBEIG(J) + 1
C           ... IBEIG = (approximate) boson irrep
         END IF
         IF (NB_CNT(JSYM) .GT. 0) THEN
            NB_CNT(JSYM) = NB_CNT(JSYM) - 1
            JOCC = JOCC + 1
            IF (J .NE. JOCC) THEN
C              swap:
               IBEIGSWP = IBEIG(J)
               IBEIG(J) = IBEIG(JOCC)
               IBEIG(JOCC) = IBEIGSWP
               EIGSWP   = EIG(J)
               EIG(J)   = EIG(JOCC)
               EIG(JOCC)= EIGSWP
               DO IZ = 1,NZ
                  CALL DSWAP(NO,CFU(1,J,IZ),1,CFU(1,JOCC,IZ),1)
               END DO
            END IF
         END IF
      ENDDO
      CALL QEXIT('BOSSE2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smosel */
      SUBROUTINE SMOSEL(CFU,EIG,IBEIG,WORK,LWORK)
C***********************************************************************
C
C     Select vectors based on overlap.
C     Input is
C       CFU   - coefficients in MO-basis
C       EIG   - eigenvalues
C       IBEIG - supersymmetry info
C
C     Written by T.Saue Nov 11 1995
C     Last revision May 2002 hjaaj: include EIG and IBEIG
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION EIG(*),IBEIG(*),CFU(*),WORK(LWORK)
C Used from COMMON blocks:
C  dcbdhf: NSMOTQ
C
#include "dcbdhf.h"
#include "memint.h"
C     Memory allocation - NSMOTQ initialized in SMOGEN
      CALL MEMGET2('REAL','SMO',KSMO,NSMOTQ,WORK,KFREE,LFREE)

      CALL SMOSE1(EIG,IBEIG,CFU,WORK(KSMO),WORK(KFREE),LFREE)
C     Memory deallocation
      CALL MEMREL('SMOSEL.SMOSE1',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smose1 */
      SUBROUTINE SMOSE1(EIG,IBEIG,CFU,SMO,WORK,LWORK)
C***********************************************************************
C
C     Select vectors based on overlap
C
C     Written by T.Saue Nov 11 1995
C     Last revision May 2002 hjaaj: include EIG and IBEIG
C                      Aug 2007 MI: extended for all electronic shells
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION EIG(*),IBEIG(*),CFU(*),SMO(*),WORK(LWORK)
C
C Used from COMMON blocks:
C  DCBORB: NSMOTQ,NOCC
C  DGROUP: NFSYM
#include "dcbdhf.h"
#include "dcborb.h"
#include "dgroup.h"
      CALL QENTER('SMOSE1')
#include "memint.h"
C
C     Read vectors to select from
C
      CALL OPNFIL(LUSMOS,'DFSMOS','OLD','SMOSE1')
      CALL READT(LUSMOS,NSMOTQ,SMO)
      CLOSE(LUSMOS,STATUS='KEEP')
C
      DO I = 1,NFSYM
      IF(NESHMF(I).NE.0) THEN
C       Memory allocation
        CALL MEMGET2('REAL','EBUF'  ,KEBUF,NESHMF(I),WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBEBUF',KIBEB,NESHMF(I),WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','CBUF'  ,KCBUF,NESHMF(I)*NTMO(I)*NZ,
     &               WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','SBUF'  ,KSBUF,NESHMF(I)*NESHMF(I)*NZ,
     &               WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','RCSEL' ,KRSEL,NESHMF(I),WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','ICSEL' ,KISEL,2*NESHMF(I),WORK,KFREE,LFREE)
        CALL SMOSE2(I,EIG(IORB(I)+1),IBEIG(IORB(I)+1),
     &              CFU(I2ORBT(I)+1),SMO(ISMOQ(I)+1),
     &              WORK(KEBUF),WORK(KIBEB),WORK(KCBUF),WORK(KSBUF),
     &              WORK(KRSEL),WORK(KISEL),
     &              NTMO(I),NOCCMF(I),NESHMF(I),WORK(KFREE),LFREE)
C     Memory deallocation
        CALL MEMREL('SMOSE1.SMOSE2',WORK,KWORK,KWORK,KFREE,LFREE)
      ENDIF
      ENDDO
      CALL QEXIT('SMOSE1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smose2 */
      SUBROUTINE SMOSE2(IFRP,EIG,IBEIG,CFU,SMO,EBUF,IBEBUF,CBUF,SBUF,
     &                  RCSEL,ICSEL,KNORB,KNOCC,KNESH,WORK,LWORK)
C***********************************************************************
C
C     Select vectors  fermion ircop IFRP based on overlap
C
C     Reference vector J has largest overlap with MO-coefficient ICSEL(J,1)
C     MO-coefficient J has largest overlap with refernce vector ICSEL(J,2)
C
C     Written by T.Saue Nov 11 1995
C     Major revisions: Jan 7 1998 - jth NORB(IFRP) etc. is now parameter
C                      May 2002 hjaaj: EIG and IBEIG added
C                      Aug 2007 MI: extended for all electronic shells
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0)
C
C Used from COMMON blocks:
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbdhf.h"
C
C     KNORB = NORB(IFRP)
C     KNOCC = NOCCMF(IFRP)
C     KNESH = NESHMF(IFRP)
C
      CHARACTER FSYM*1
      DIMENSION EIG(KNORB),IBEIG(KNORB),
     &          CFU(KNORB,KNORB,NZ),
     &          SMO(KNORB,KNESH,NZ),
     &          SBUF(KNESH,KNESH,NZ),
     &          EBUF(KNESH),IBEBUF(KNESH),
     &          CBUF(KNORB,KNESH,NZ),
     &          RCSEL(KNESH),ICSEL(KNESH,2),
     &          WORK(LWORK)
C
C            Form product :
C            --------------
C    CFU(KNESH,KNORB)^+ * SMO(KNORB,KNESH) =   SBUF(KNESH,KNESH)
C
      CALL QGEMM(NESHMF(IFRP),NESHMF(IFRP),NTMO(IFRP),D1,
     &       'H','N',IPQTOQ(1,0),CFU(1,NPSH(IFRP)+1,1),
     &                NTMO(IFRP),NTMO(IFRP),NZ,
     &       'N','N',IPQTOQ(1,0),SMO,NTMO(IFRP),NESHMF(IFRP),NZ,
     &       D0,IPQTOQ(1,0),SBUF,NESHMF(IFRP),NESHMF(IFRP),NZ)

      IF(IPRSCF.GE.10) THEN
        CALL HEADER('SMOSE2: SBUF-matrix',-1)
        CALL PRQMAT(SBUF,NESHMF(IFRP),NESHMF(IFRP),NESHMF(IFRP),
     &              NESHMF(IFRP),NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Form squared overlaps
C
      DO J = 1,NESHMF(IFRP)
        DO I = 1,NESHMF(IFRP)
          SBUF(I,J,1) = SBUF(I,J,1)*SBUF(I,J,1)
        ENDDO
      ENDDO
      DO IZ = 2,NZ
        DO J = 1,NESHMF(IFRP)
          DO I = 1,NESHMF(IFRP)
            SBUF(I,J,1) = SBUF(I,J,1) + SBUF(I,J,IZ)*SBUF(I,J,IZ)
          ENDDO
        ENDDO
      ENDDO
      IF(IPRSCF.GE.6) THEN
        CALL HEADER('SMOSE1: Squared-matrix',-1)
        CALL PRQMAT(SBUF,NESHMF(IFRP),NESHMF(IFRP),NESHMF(IFRP),
     &              NESHMF(IFRP),1,IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Perform selection
C
      N2DIM = NTMO(IFRP)*NESHMF(IFRP)
      DO IZ = 1,NZ
        CALL DCOPY(N2DIM,CFU(1,NPSH(IFRP)+1,IZ),1,CBUF(1,1,IZ),1)
      ENDDO
      CALL DCOPY(KNESH,EIG(NPSH(IFRP)+1),1,EBUF,1)
      CALL ICOPY(KNESH,IBEIG(NPSH(IFRP)+1),1,IBEBUF,1)
      CALL ICOPY(2*KNESH,0,0,ICSEL,1)

      ISELM = 0
      DO J = 1, NESHMF(IFRP)
        ISEL = IDAMAX(NESHMF(IFRP),SBUF(1,J,1),1)
        ICSEL(J,1) = ISEL
        ICSEL(ISEL,2) = J
        RCSEL(J) = SBUF(ISEL,J,1)
        ISELM = MAX(ISELM,ISEL)
C       Zero out the selected overlap
        CALL DCOPY(NESHMF(IFRP),D0,0,SBUF(ISEL,1,1),NESH(IFRP))
C       Copy buffer into vector
        JJ = NPSH(IFRP) + J
        DO IZ = 1,NZ
          CALL DCOPY(NTMO(IFRP),CBUF(1,ISEL,IZ),1,CFU(1,JJ,IZ),1)
        ENDDO
C       Copy eigenvalue and (approx.) boson symmetry to new place:
        EIG(JJ) = EBUF(ISEL)
        IBEIG(JJ) = IBEBUF(ISEL)
      ENDDO
C
C     Print section
C
      WRITE(FSYM,'(I1)') IFRP
      CALL HEADER(
     & 'SMOSEL: Reordered (electronic) vectors in fermion ircop '
     &            //FSYM,-1)
      NVEC = 0
      DO J = 1,NESHMF(IFRP)
      IF(ICSEL(J,1).NE.J) THEN
        NVEC = NVEC + 1
        WRITE(LUPRI,'(I5,A,I5,A,F20.10)') ICSEL(J,1),' --> ',J,
     &                                    ' OVERLAP:',RCSEL(J)
      ENDIF
      ENDDO
C
      IF(NVEC.EQ.0) WRITE(LUPRI,'(A)') 'No reordering.'
      IF(IPRSCF.GE.3) THEN
        WRITE(LUPRI,'(/A)') '  I  ICSEL(I,1)  ICSEL(I,2)    RCSEL(I)'
        DO I = 1,ISELM
          WRITE(LUPRI,'(I4,2I12,F12.4)')I,ICSEL(I,1),ICSEL(I,2),RCSEL(I)
        END DO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck denorb */
      SUBROUTINE DENORB(DMAT,JVEC,NVEC,IFRP,CMO,IPRINT)
C***********************************************************************
C
C     Construct density matrix for NVEC orbitals of fermion corep IFRP
C     starting with orbital JVEC (with respect to NPSH(IFRP))
C
C     Written by T.Saue Oct 22 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0,D0=0.0D0)
      DIMENSION DMAT(*),CMO(*)
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
C
CTROND
      CALL DZERO(DMAT,N2BBASXQ)
      JOFF = NPSH(IFRP) + JVEC
      CALL DENST1(DMAT(I2BASX(IFRP,IFRP)+1),NTBAS(0),NTBAS(0),NZ,D1,D0,
     &              CMO(ICMOQ(IFRP)+1),NFBAS(IFRP,0),NORB(IFRP),
     &              JOFF,NVEC,NFBAS(IFRP,0))
C
C     Print section
C
      IF(IPRINT.GE.4) THEN
        CALL TITLER('Output from DENORB','*',103)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',IFRP
        CALL PRQMAT(DMAT(I2BASX(IFRP,IFRP)+1),NFBAS(IFRP,0),
     &              NFBAS(IFRP,0),NTBAS(0),NTBAS(0),
     &              NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck oneerg */
      SUBROUTINE ONEERG(EIG,IBEIG)
C**********************************************************************
C
C     PURPOSE:
C   ===========
C        Calculate one-electron energy as sum of eigenvalues
C
C        Written by T.Saue Jan 28 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D2=2.0D0)
C
#include "frame.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbdhf.h"
      CHARACTER FMT*6,MXFORM*6
      DIMENSION EIG(*)
C
      CALL QENTER('ONEERG')
C
      ERGBUF = DHFERG
      E1PART = D0
      EA1PART = D0
      DO IFRP = 1,NFSYM
        IF(NOCC(IFRP).GT.0) THEN
           IF (NISHT .GT. 0) THEN
              E1PART = E1PART + D2*DSUM(NISH(IFRP),
     &                         EIG(IORB(IFRP)+NPSH(IFRP)+1),1)
              IF (IPRSCF.GE.5) THEN
               write(lupri,*) 'ONEERG: eigenvalues:',
     &          (EIG(IORB(IFRP)+NPSH(IFRP)+i),i=1,NISH(IFRP))
               write(lupri,*) '2*sum: ',D2*DSUM(NISH(IFRP),
     &                         EIG(IORB(IFRP)+NPSH(IFRP)+1),1)
              ENDIF
           END IF
           IF (AOC) THEN
             IDX = IORB(IFRP)+NPSH(IFRP)+NISH(IFRP)+1
             DO IOPEN = 1,NOPEN
                E1PART = E1PART + D2*DF(IOPEN)*DSUM(NACSH(IFRP,IOPEN),
     &                    EIG(IDX),1)
                IDX = IDX + NACSH(IFRP,IOPEN)
             END DO
           END IF
        ENDIF
      ENDDO
      ELERGY = E1PART
      DHFERG = ELERGY + POTNUC
C
      IF(IPRSCF.GE.1) THEN
        Call HEADER('Output from ONEERG',-1)
        ERGMAX = ABS(E1PART+EA1PART)*10
        FMT = MXFORM(ERGMAX,20)
        IF (IPRSCF .GE. 2) THEN
           WRITE(LUPRI,'(A,'//FMT//')')
     &     'Inactive energy           :  ',E1PART
        END IF
        WRITE(LUPRI,'(A,'//FMT//')')
     +     'Electronic energy         :  ',ELERGY,
     +     'Nuclear potential energy  :  ',POTNUC,
     +     'Total energy              :  ',DHFERG
      ENDIF
C
      CALL QEXIT('ONEERG')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reord */
      SUBROUTINE REORD(CMO,EIG,IBEIG,IMOORD,IREORD,MXREORD)
C**********************************************************************
C
C     PURPOSE:
C        Reorder orbitals (either in pre-SCF or post-SCF stage)
C
C        Written by J. Thyssen Nov 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION CMO(*),EIG(*),IBEIG(*)
      DIMENSION IMOORD(MXREORD,2),IREORD(2)

      real(8), allocatable :: cmo_save(:)
      real(8), allocatable :: eig_save(:)
      integer, allocatable :: ibeig_save(:)

#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
C
      CALL QENTER('REORD')
C
      WRITE(LUPRI,'(/A)') ' -> Reordering orbitals as requested'
      DO 10 I = 1,NFSYM
         IF (IREORD(I).EQ.0) GOTO 10
         IF (IPRSCF.GE.2) THEN
            CALL HEADER('REORD: before reordering',-1)
            IF (NFSYM.GT.1) WRITE(LUPRI,'(5X,A,I2/)')
     &         'for fermion ircop',I
            WRITE(LUPRI,'(3X,A)')
     &         '* Electron eigenvalues and boson classification ...'
            WRITE(LUPRI,'(3X,I5,1X,F19.8,I5)')
     &         (J,EIG(J),IBEIG(J),J=1+IORB(I)+NPSH(I),IORB(I)+NORB(I))
         END IF
         IF (IPRSCF.GE.10) THEN
            WRITE(LUPRI,'(3X,A)') '* Eigenvectors...'
            CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &                  NFBAS(I,0),NORB(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
C        Save coefficients
C
         allocate(cmo_save(ncmoq(i)))
         allocate(eig_save(norb(i)))
         allocate(ibeig_save(norb(i)))

         CALL DCOPY(NCMOQ(I),CMO(1+ICMOQ(I)),1,cmo_save,1)
         CALL DCOPY(NORB(I),EIG(1+IORB(I)),1,eig_save,1)
         CALL ICOPY(NORB(I),IBEIG(1+IORB(I)),1,ibeig_save,1)
C
C        Reorder coefficients
C
         CALL REOR1(CMO(1+ICMOQ(I)),cmo_save,
     &              EIG(1+IORB(I)),eig_save,
     &              IBEIG(1+IORB(I)),ibeig_save,
     &              NFBAS(I,0),NORB(I),NZ,
     &              IMOORD(1,I),IREORD(I),NPSH(I),NESH(I))
C
C        Release memory
C
         deallocate(cmo_save)
         deallocate(eig_save)
         deallocate(ibeig_save)
C
         IF (IPRSCF.GE.2) THEN
            CALL HEADER('REORD: after reordering',-1)
            IF (NFSYM.GT.1) WRITE(LUPRI,'(5X,A,I2/)')
     &         'for fermion ircop',I
            WRITE(LUPRI,'(3X,A)')
     &         '* Electron eigenvalues and boson classification ...'
            WRITE(LUPRI,'(3X,I5,1X,F15.8,I5)')
     &         (J,EIG(J),IBEIG(J),J=1+IORB(I)+NPSH(I),IORB(I)+NORB(I))
         END IF
         IF (IPRSCF.GE.10) THEN
            WRITE(LUPRI,'(3X,A)') '* Eigenvectors...'
            CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),LUPRI)
         END IF
 10   CONTINUE
      CALL FLSHFO(LUPRI)
      CALL QEXIT('REORD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reor1 */
      SUBROUTINE REOR1(CMO,CMOBCK,EIG,EIGBCK,IBEIG,IBEIGBCK,
     &                 NBASI,NORBI,NZ,IORDER,NORDER,
     &                 NPOSITRON,NELECTRON)
C**********************************************************************
C
C     PURPOSE:
C        Reorder orbitals
C
C        Written by J. Thyssen Nov 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION CMO(NBASI,NORBI,NZ), CMOBCK(NBASI,NORBI,NZ),
     &          EIG(*),EIGBCK(*), IBEIG(*), IBEIGBCK(*),
     &          IORDER(NORDER)
C
C
      CALL QENTER('REOR1')
C
      ICHK_S = 0
      ICHK_D = 0
      DO I = 1,NORDER
         IORB_S = IORDER(I)+NPOSITRON
         IORB_D = I+NPOSITRON
         IF (IORB_D .GT. NELECTRON+NPOSITRON) THEN
            WRITE(LUPRI,*) 'IORB_D=',IORB_D
            WRITE(LUPRI,*) 'NELECTRON+NPOSITRON=',NELECTRON+NPOSITRON
            CALL QUIT('REOR1 reordering to non-existent orbital')
         ELSE
            ICHK_S = ICHK_S + IORB_S
            ICHK_D = ICHK_D + IORB_D
            if(IORB_S == IORB_D) goto 10
C
C           Copy eigenvalue and boson classification
C
            EIG(IORB_D)   = EIGBCK(IORB_S)
            IBEIG(IORB_D) = IBEIGBCK(IORB_S)
C
C           Copy eigenvector
C
            DO IZ = 1,NZ
               CALL DCOPY(NBASI,CMOBCK(1,IORB_S,IZ),1,
     &                          CMO(1,IORB_D,IZ),1)
            END DO
 10         CONTINUE
         END IF
      END DO
C
C     Make sure we are not copying more than one orbital to the same 
C     destination,
C     or that IORDER(I) does not belong to 1..NORDER
C     (not 100% test: if user specified e.g. IORDER=2,2,2 this would 
C     not catch it,
C     but that is quite unlikely) /hjaaj May 2002
C
      IF (ICHK_S .NE. ICHK_D) THEN
CMI     ... more control print-out necessary
         WRITE(LUPRI,'(/,2x,a,i5)') 'NPOSITRON=',NPOSITRON
         WRITE(LUPRI,'(2x,a,i5,a,i5,/)')
     &     'ICHK_S=',ICHK_S,' ICHK_D=',ICHK_D
         DO I=1,NORDER
           WRITE(LUPRI,'(2x,a,i4,a,i4,a,i4)')
     &     'IORDER(',I,'/',NORDER,')=',IORDER(I)
         ENDDO
         CALL QUIT('REOR1: destination orbs. .ne. source orbitals!')
      ENDIF
C
      CALL QEXIT('REOR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MKMOFK2 */
      SUBROUTINE MKMOFK2(FMO,FAO,TMAT,WORK,LWORK)
C**********************************************************************
C
C     PURPOSE:
C        Construct MO Fock matrix from Fock matrices in AO-basis
C     WARNING:
C        This is a small part of the routine MKMOFK distributed
C           with the open shell code.
C
C        Written by J. Thyssen Dec 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( DP5 = 0.50D00 , D0 = 0.00D00, D1 = 1.00D00,
     &            DM1 = -1.00D00, D2 = 2.00D00)
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION FMO(*),FAO(N2BBASXQ,NFMAT),TMAT(*),WORK(*)
      CHARACTER CTMP*3
C
#include "memint.h"
C
      CALL QENTER('MKMOFK2')
C
C     Transform FD to MO-basis
C
      DO I = 1,NFSYM
      IF(NTMO(I).GT.0) THEN
        IF(SUB_BL) THEN
          CALL TSUBBL(I,FAO,FMO,TMAT,TMAT,WORK(KFREE),LFREE,IPRSCF)
        ELSE
          CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &            NTMO(I),NTMO(I),
     &            FAO(I2BASX(I,I)+1,1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &            TMAT(1+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &            TMAT(1+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRSCF)
        ENDIF
        IF(IPRSCF.GE.10 ) THEN
          WRITE(CTMP,'(I1)') I
          CALL HEADER('MKMOFK2: D Fock matrix in MO-basis'//
     &                ', corep '//CTMP,-1)
          CALL PRQMAT(FMO(1+I2TMOT(I)),NTMO(I),NTMO(I),
     &                NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        END IF
      END IF
      END DO
      CALL QEXIT('MKMOFK2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck EXTRSB */
      SUBROUTINE EXTRSB (IMODE,I,NORB_SUB_L,M_SUB_BL_L,N_SUB_BL_L,
     &                   NORB_SUB_R,M_SUB_BL_R,N_SUB_BL_R,
     &                   ISUB,JSUB,ICOMP,JCOMP,IZS,IZE,
     &                   FSUB,MSUB,NSUB,NZS,
     &                   FFUL,MFUL,NFUL,NZ)
C**********************************************************************
C     PURPOSE:
C        Extract or insert blocks out of (into) a full matrix 
C
C        General parameters :
C        IMODE                 - 1 : extract, 2 : insert
C        I                     Parity (fermion irrep)
C        IZS                   First quaternion unit
C        IZE                   Last quaternion unit
C        NZ                    Number of quaternion units
C        NZS                   Number of active quaternion units
C
C        Left index :
C        NORB_SUB_L            Block dimensions left index
C        M_SUB_BL_L            First dimension of NORB_SUB_L
C        N_SUB_BL_L            Number of blocks left index
C        ISUB                  Active block left (zero : all blocks)
C        ICOMP                 Active component left (zero : all comp.)
C
C        Right index :
C        NORB_SUB_R            Block dimensions right index
C        M_SUB_BL_R            First dimension of NORB_SUB_R
C        N_SUB_BL_R            Number of blocks right index
C        JSUB                  Active block right (zero : all blocks)
C        JCOMP                 Active component right (zero : all comp.)
C
C        Input and output :
C        FSUB                  Sub matrix
C        FFUL                  Full matrix
C
C        Written for for boson symmetry (august 2000
C        Generalized to subblocks (november 2003)
C        L. Visscher
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION FSUB(MSUB,NSUB,NZS),FFUL(MFUL,NFUL,NZ)
      DIMENSION NORB_SUB_L(M_SUB_BL_L,2,0:2)
      DIMENSION NORB_SUB_R(M_SUB_BL_R,2,0:2)
      LOGICAL SUBI, SUBJ
C
C     The full matrix is ordered as (function index, component type, sub block),
C     the submatrix should be ordered as (function index, sub blocK)
C
C     if isub is zero we extract blocks from all sub blocks with parity I
C     if isub has a non-zero value we just extract the block of this type.
C
      IF (ISUB.EQ.0) THEN
         ISS = 1
         ISE = N_SUB_BL_L
      ELSE
         ISS = ISUB
         ISE = ISUB
      ENDIF
C
      IF (JSUB.EQ.0) THEN
         JSS = 1
         JSE =  N_SUB_BL_R
      ELSE
         JSS = JSUB
         JSE = JSUB
      ENDIF
C
C     Extract the subblock by looping over both indices in the right order
C     while using different counting variables for the full and the sub matrices
C     Whenever we encounter a block that fulfills the criteria 
C     we copy it into or out from the subbblock.
C
      DO IZ = IZS, IZE
        JXF = 0
        JXS = 0
C
C       Loop over all blocks for the right index.
C
        DO JS = 1, N_SUB_BL_R
C
C         The small component functions come first, loop from S to L
C
          DO JC = 2, 1, -1
            SUBJ = JC.EQ.JCOMP.AND.JS.GE.JSS.AND.JS.LE.JSE
            DO 100 JX = 1, NORB_SUB_R(JS,I,JC)
              JXF = JXF + 1
              IF (SUBJ) THEN
C                Only when we are in a selected block we update this pointer
                 JXS = JXS + 1
              ELSE
C
C             We are not going to run the loops for the left index 
C             and have already updated JXF so we might as well skip 
C             the inner loops
C
                 GOTO 100
              ENDIF
              IXF = 0
              IXS = 0
C
C             Loop over all blocks for the left index.
C
              DO IS = 1, N_SUB_BL_L
C
C               Loop from S to L for the left index.
C
                DO IC = 2, 1, -1
                  SUBI = IC.EQ.ICOMP.AND.IS.GE.ISS.AND.IS.LE.ISE
                  DO IX = 1, NORB_SUB_L(IS,I,IC)
                    IXF = IXF + 1
C
C                   Check the criteria and fill in
C
                    IF (SUBI.AND.SUBJ) THEN
                       IXS = IXS + 1
                       IF (IMODE.EQ.1) THEN
                          FSUB(IXS,JXS,IZ) = FFUL(IXF,JXF,IZ)
                       ELSEIF (IMODE.EQ.2) THEN
                          FFUL(IXF,JXF,IZ) = FSUB(IXS,JXS,IZ)
                       ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
 100        CONTINUE
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck EXTRSBR */
      SUBROUTINE EXTRSBR (IMODE,I,ISYM,JSYM,ICOMP,JCOMP,IZS,IZE,
     &                    FSUB,MSUB,NSUB,NZS,
     &                    FFUL,MFUL,NFUL)
C**********************************************************************
C     PURPOSE:
C        Extract or insert a sub block out of the full matrix 
C        that is ordered on boson irreps for the left and on small, 
C        large or positron, electron on the right index.
C
C        IMODE                 - 1 : extract, 2 : insert
C
C        Written by L. Visscher
C        September 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION FSUB(MSUB,NSUB,NZS),FFUL(MFUL,NFUL,NZ)
      LOGICAL SUBI, SUBJ
C
C     Calculate number of boson symmetries : the full matrix is ordered
C     as (function index, component type, boson irrep), the submatrix
C     should be ordered as (function index, boson irrep)
C
      NBRP = 4 / NZ
C
C     When isym is zero we extract blocks from all boson irreps that
C     have the right parity, if isym has a non-zero value we just extract
C     the block of this symmetry.
C
      IF (ISYM.EQ.0) THEN
         ISS = 1
         ISE = 4 / NZ
      ELSE
         ISS = ISYM
         ISE = ISYM
      ENDIF
C
      IF (JSYM.EQ.0) THEN
         JSS = 1
         JSE = 4 / NZ
      ELSE
         JSS = JSYM
         JSE = JSYM
      ENDIF
C
C     Extract the subblock by looping over both indices in the 
C     right order while using different counting variables for the 
C     full and the sub matrices. Whenever we encounter a block that 
C     fulfills the criteria we copy it into or out from the subbblock.
C
      DO IZ = IZS, IZE
        JXF = 0
        JXS = 0
C
C       The small component functions come first, loop from S to L
C
        DO JC = 2, 1, -1
C
C       Loop over the boson irreps for the right index.
C
          DO JS = 1, NBRP
            SUBJ = JC.EQ.JCOMP.AND.JS.GE.JSS.AND.JS.LE.JSE
            DO 100 JX = 1, NBORB(JS,I,JC)
              JXF = JXF + 1
              IF (SUBJ) THEN
C                Only when we are in a subblock we update this pointer
                 JXS = JXS + 1
              ELSE
C
C             We are not going to run the loops for the left index 
C             and have already updated JXF so we might as well skip 
C             the inner loops
C
                 GOTO 100
              ENDIF
              IXF = 0
              IXS = 0
C
C             Loop over the boson irreps for the left index.
C
              DO IS = 1, NBRP
C
C               Loop from S to L for the left index.
C
                DO IC = 2, 1, -1
                  SUBI = IC.EQ.ICOMP.AND.IS.GE.ISS.AND.IS.LE.ISE
                  DO IX = 1, NBORB(IS,I,IC)
                    IXF = IXF + 1
C
C                   Check the criteria and fill in
C
                    IF (SUBI.AND.SUBJ) THEN
                       IXS = IXS + 1
                       IF (IMODE.EQ.1) THEN
                          FSUB(IXS,JXS,IZ) = FFUL(IXF,JXF,IZ)
                       ELSEIF (IMODE.EQ.2) THEN
                          FFUL(IXF,JXF,IZ) = FSUB(IXS,JXS,IZ)
                       ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
 100        CONTINUE
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck selbos */
      SUBROUTINE SELBOS (IFSYM,NO,INDX,EIG,IBEIG,EVEC,RDUM,IDUM,IPRINT)
C***********************************************************************
C
C     PURPOSE:
C        Sort eigenvectors from block ordering to energy ordering
C        On entry one should provide an array IBEIG with identification
C        of the spinors.
C
C        Written by L. Visscher
C        April 1999
C
C        Last modifications: MI, august 2005, April 2009 (Tel Aviv)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
      DIMENSION EIG(NO),INDX(NO),IBEIG(NO),IDUM(NO)
      DIMENSION EVEC(NO,NO),RDUM(NO,NO)

      CALL QENTER('SELBOS')

      IF (IPRINT.GE.7) THEN 
       CALL HEADER('*** Output from SELBOS - before resorting ***',-1)
       write(lupri,'(2x,a,i1,a,i1,a,i4)')
     & 'fermion symmetry IFSYM/NFSYM=',IFSYM,'/',NFSYM,
     & ' NORB(IFSYM)=',NO
       write(lupri,*)
     & 'entering eigenvalues and boson irreps of each MO:'
       DO J=1,NO
        write(lupri,'(I3,a,f18.10,a,i2,a)')
     &  J,'. eigval=',EIG(J),' boson irrep=',IBEIG(J),
     &  ' corresponding eigenvector:'
        write(lupri,'(5f16.8)') (EVEC(K,J),K=1,NO)
       ENDDO
      ENDIF
C
C     Establish the order
C
      CALL INDEXX (NO,EIG,INDX)
C
C     Sort the eigenvalues and their IDs
C
      DO J = 1, NO
         RDUM(J,1) = EIG(INDX(J))
         IDUM(J) = IBEIG(INDX(J))
      ENDDO
      CALL DCOPY(NO,RDUM,1,EIG,1)
      CALL ICOPY(NO,IDUM,1,IBEIG,1)

C
C     Sort the eigenvectors
C
      DO J = 1, NO
         CALL DCOPY(NO,EVEC(1,INDX(J)),1,RDUM(1,J),1)
      ENDDO
      CALL DCOPY(NO*NO,RDUM,1,EVEC,1)

      IF (IPRINT.GE.7) THEN
       CALL HEADER('SELBOS: after resorting ',-1)
       write(lupri,'(2x,a,i1,a,i1,a,i4)')
     &  'fermion symmetry, IFSYM/NFSYM=',IFSYM,'/',NFSYM,
     &  ' NORB(IFSYM)=',NO
       write(lupri,*)
     &      'resorted eigenvalues and boson irreps of each MO:'
       DO J=1,NO
        write(lupri,'(I3,a,f18.10,a,i2,a)')
     &  J,'. eigval=',EIG(J),' boson irrep=',IBEIG(J),
     &  ' correspoding eigenvector:'
        write(lupri,'(5f16.8)') (EVEC(K,J),K=1,NO)
       ENDDO
      ENDIF

C
      CALL QEXIT('SELBOS')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck appbos */
      SUBROUTINE APPBOS (IFRP,NO,INDX,IBEIG,EVEC,RDUM)
C***********************************************************************
C
C     PURPOSE:
C        Identify vectors by their largest boson character
C        This approximate labeling can be used in the CI code.
C
C        Works only for the modified Dirac equation
C        (restricted kinetic balance)
C
C     Called from: DFDIAG
C
C        Written by L. Visscher
C        October 2000
C        Last modifications: M.Ilias, Aug 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
      DIMENSION IBEIG(NO),INDX(NO)
      DIMENSION EVEC(NO,NO,NZ),RDUM(NO,4)

CC    CALL QENTER('APPBOS') 
C
      IF (URKBAL) RETURN
C
      NBRP = 4 / NZ
C
C     Calculate the weight of each boson irrep in each vector.
C     Because we work in the orthonormal basis the total weights
C     add up to one.
C
      CALL DZERO(RDUM,4*NO)
      DO I = 1, NO
         IBEIG(I) = -1
      ENDDO
C
      IX = 1
      DO ISYM = 1, NBRP
        NBO = NBORB(ISYM,IFRP,0)
        IF (NBO.GT.0) THEN
          DO IZ = 1, NZ
             DO I = 1, NO
                RDUM(I,ISYM) = RDUM(I,ISYM) +
     &                         DDOT(NBO,EVEC(IX,I,IZ),1,EVEC(IX,I,IZ),1)
             ENDDO
          ENDDO
          IX = IX + NBO
        ENDIF
      ENDDO
C
      DO ISYM = 1, NBRP
         IREP = JFSYM(ISYM,IFRP) - 1
         NBO = NBORB(ISYM,IFRP,0)
         NSKIP = 0
C
C        Establish the order, indexx will give the vectors with the
C        highest occupation for this boson irreps last
C
         CALL INDEXX (NO,RDUM(1,ISYM),INDX)
C
C        Label the vectors, skip one if already claimed by another
C        irrep. This is not fullproof but will at least satisfy the
C        constraint that every boson irrep gets the correct number
C        of vectors.
C
         DO J = NO, NO-NBO+1, -1
            IBOLD = IBEIG(INDX(J))
            IF (IBOLD.EQ.-1) THEN
               IBEIG(INDX(J)) = IREP
            ELSE
               NSKIP = NSKIP + 1
            ENDIF
         ENDDO
C
         ISKIP = 0
         DO J = NO-NBO,1,-1
            IBOLD = IBEIG(INDX(J))
            IF (IBOLD.EQ.-1.AND.ISKIP.LT.NSKIP) THEN
               ISKIP = ISKIP + 1
               IBEIG(INDX(J)) = IREP
            ENDIF
         ENDDO
      ENDDO
C
      IF (IPRHAM.GE.8) THEN
       CALL HEADER('*** Output from APPBOS ***',-1)
       do i = no/2+1, no
       write (LUPRI,'(3(A,I4),A,4F8.2)')
     & 'fermion irrep ',ifrp,' electron function ',i-no/2,' labeled ',
     &  ibeig(i),' percentages ',(rdum(i,ix),ix=1,nbrp)
       enddo
      ENDIF

CC    CALL QEXIT('APPBOS') 
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mvogen */
      SUBROUTINE MVOGEN()
C***********************************************************************
C
C     Generate modified virtual orbitals
C
C     Written by T.Saue Sep 18 1998
C
C***********************************************************************

      use memory_allocator

#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
C
#include "dcbfir.h"
#include "dcbgen.h"
#include "dcbmvo.h"
C
      LOGICAL SAVEFLAGS(4)
      real(8), allocatable :: WORK(:)
C
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in MVOGEN')

      CALL TITLER('Modified virtual orbitals','*',125)
C
      CALL MEMGET2('REAL','CMO  ',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EIG  ',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBEIG',KIBE,NTBAS(0),      WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','FOCK ',KFOCK,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DMAT ',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)
C
      call SaveTaskDistribFlags(saveflags)
      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. , .TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      if (NPOS.GT.0) THEN
         CALL MEMGET2('INTE','NPOS',KPOS,NPOS,WORK,KFREE,LFREE)
      else
         KPOS = KFREE
      endif
C
      CALL MVOGEN_1(WORK(KCMO),WORK(KEIG),WORK(KIBE),
     &            WORK(KFOCK),WORK(KDMAT),
     &            WORK(KPOS),WORK(KFREE),LFREE)
      CALL MEMREL('MVOGEN',WORK,KWORK,KWORK,KFREE,LFREE)
      call dealloc(WORK)
C
      IF(PARCAL) call SetTaskDistribFlags(saveflags)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mvogen_1 */
      SUBROUTINE MVOGEN_1(CMO,EIG,IBEIG,FOCK,DMAT,NPOS,WORK,LWORK)
C***********************************************************************
C     Generate modified virtual orbitals
C
C     Written by T.Saue June 18 1998
C     Last revision : June 18 1998 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,DM1=-1.0D0,DM2=-2.0D0,D2=2.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbmvo.h"
#include "dcbham.h"
      LOGICAL   TOBE
      DIMENSION CMO(*),EIG(*),IBEIG(*),FOCK(*),DMAT(*),NPOS(*),WORK(*)
      DIMENSION IC(2),IE(2),ID(2)
      CHARACTER CTMP*1
      CALL QENTER('MVOGEN')
#include "memint.h"
C
C     Initialization
C
      LUBUF = 22
C     full Fock matrix
      IFCMVO = 1
C     only direct (Coulomb) contributions 
      IF(IF2MVO.EQ.1) IFCMVO = 2
C     only exchange contributions
      IF(IF2MVO.EQ.2) IFCMVO = 3
      ISYMVO = 1
      IHRMVO = 1
C
C     Get all coefficients
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(//A)')
     &   'MVOGEN ERROR: MO coefficient file DFCOEF not present'
         CALL QUIT('MVOGEN: No MO coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','MVOGEN')
         IOPT=14
         IF(ATOMIC) THEN
           CALL ATOMIC_REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,IOPT)           
         ELSE
           CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,IOPT)
         ENDIF
      ENDIF
      IF(IONMVO) THEN
C
C       Generate MVOs by the Bauschlicher method (ionic Fock)
C
        CALL DENMVO(DMAT,CMO,FOCK,WORK(KFREE),LFREE)
        IF (ABS(MVOWEIGHT-1.0D0).GT.1.0D-14) THEN
C     Make a weighted average density D = w*Dmvo + (1 - w)D0, where
C     w = MVOWEIGHT and D0 is the full inactive density matrix.
           CALL GENDEN(FOCK,CMO,1,IPRMVO)
           CALL DSCAL(N2BBASXQ,DMAT,MVOWEIGHT,1)
           CALL DAXPY(N2BBASXQ,1-MVOWEIGHT,FOCK,1,DMAT,1)
        ENDIF
        CALL TWOFCK(ISYMVO,IHRMVO,IFCMVO,FOCK,DMAT,1,NPOS,
     &              INTMVO,IPRMVO,WORK(KFREE),LFREE)
        INQUIRE(FILE='DFFCK1',EXIST=TOBE)
        IF(TOBE) THEN
          CALL OPNFIL(LUBUF,'DFFCK1','OLD','MVOGEN_1')
          CALL REAFCK(LUBUF,DMAT,.TRUE.,1)
          CLOSE(LUBUF,STATUS='KEEP')
        ELSE
          CALL ONEFCK(DMAT,IPRINT,WORK(KFREE),LFREE)
        ENDIF
        CALL DAXPY(N2BBASXQ,D1,DMAT,1,FOCK,1)
      ELSE
C
C       Construct density matrix for chosen vectors
C
        IF(ADDREP) THEN
          CALL DENMVO(FOCK,CMO,DMAT,WORK(KFREE),LFREE)
C         ... generate inactive density matrix
          CALL GENDEN(DMAT,CMO,1,IPRMVO)
          CALL DAXPY(N2BBASXQ,DM2,FOCK,1,DMAT,1)
        ELSE
          CALL DENMVO(DMAT,CMO,FOCK,WORK(KFREE),LFREE)
          CALL DSCAL(N2BBASXQ,DM1,DMAT,1)
        ENDIF
        IF(IFCMVO.EQ.3) CALL DSCAL(N2BBASXQ,DM1,DMAT,1)
C
C       Construct two-electron Fock matrix
C
        CALL TWOFCK(ISYMVO,IHRMVO,IFCMVO,FOCK,DMAT,1,NPOS,
     &              INTMVO,IPRMVO,WORK(KFREE),LFREE)
      ENDIF
C
C       First transform to virtual MO - basis;
C       this liberates FOCK. Transformed matrix in DMAT
C
      ND = 1
      DO I = 1,NFSYM
      IF (NSSH(I).GT.0 ) THEN
        IC(I) = ICMOQ(I)+NFBAS(I,0)*(NOCC(I)+NPSH(I))+1
        IE(I) = IORB(I) + NOCC(I) + NPSH(I) + 1
        ID(I) = ND
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &           NSSH(I),NSSH(I),
     &           FOCK(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           DMAT(ID(I)),NSSH(I),NSSH(I),NZ,IPQTOQ(1,0),
     &           CMO(IC(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           CMO(IC(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRMVO)
        IF (IPRMVO.GE.10 ) THEN
          WRITE(CTMP,'(I1)') I
          CALL HEADER('MVOGEN: Fock matrix in virtual MO-basis'//
     &                 ', corep '//CTMP,-1)
          CALL PRQMAT(DMAT(ID(I)),NSSH(I),NSSH(I),
     &                NSSH(I),NSSH(I),NZ,IPQTOQ(1,0),LUPRI)
        END IF
        ND = ND + NSSH(I)*NFBAS(I,0)*NZ
      ENDIF
      ENDDO
C
C     Diagonalize matrix and backtransform vectors
C
      DO I = 1,NFSYM
      IF (NSSH(I).GT.0 ) THEN
        CALL QDIAG(NZ,NSSH(I),DMAT(ID(I)),NSSH(I),NSSH(I),
     &             EIG(IE(I)),1,FOCK(ID(I)),NSSH(I),
     &             NSSH(I),WORK(KFREE),LFREE,IERR)
        CALL BCKTR1(DMAT(ID(I)),NFBAS(I,0),NSSH(I),
     &              FOCK(ID(I)),NSSH(I),NSSH(I),
     &              NSSH(I),NZ,NSSH(I),1,NFBAS(I,0),
     &              CMO(IC(I)),NFBAS(I,0),NORB(I),NZ,
     &              IPRMVO)
        IF(IPRMVO.GE.0) THEN
          WRITE(LUPRI,'(A,A3)')
     &     '* Virtual eigenvalues of ircop ',FREP(I)
          WRITE(LUPRI,'(5F15.8)') (EIG(IE(I)+J),J=0,(NSSH(I)-1))
        ENDIF
C
C       Insert new coefficients
C
        NDIM = NSSH(I)*NFBAS(I,0)
        ICC = IC(I)
        IDD = ID(I)
        DO IZ = 1,NZ
          CALL DCOPY(NDIM,DMAT(IDD),1,CMO(ICC),1)
          ICC = ICC + NCMO(I)
          IDD = IDD + NFBAS(I,0)*NSSH(I)
        ENDDO
      END IF
      END DO
      IF(ATOMIC) THEN
        CALL ATOMIC_WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG)
      ELSE
        CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG,.FALSE.)
      ENDIF   
      CLOSE(LUCOEF,STATUS='KEEP')
C
      CALL QEXIT('MVOGEN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck denmvo */
      SUBROUTINE DENMVO(DMAT,CMO,BUF,WORK,LWORK)
C***********************************************************************
C
C     Generate density matrix for chosen MVO vectors
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbmvo.h"
      DIMENSION DMAT(*),CMO(*),BUF(*),WORK(*)
      DIMENSION NR(2),IR(2,2),KVEC(2,2),NVEC(2,0:2)
      DIMENSION ICMV(2)
#include "memint.h"
      KFRSAV = KFREE
C
C     Set up index arrays for coefficients contributing
C     =================================================
C
      IC   = 1
      DO I = 1,NFSYM
        NR(I)   = NESH(I)
        IR(1,I) = 0
        IR(2,I) = NESH(I)
      ENDDO
      DO IFRP = 1,NFSYM
         CALL MEMGET2('INTE','ILIST',KVEC(IFRP,1),NR(IFRP),
     &                WORK,KFREE,LFREE)
        NVEC(IFRP,0) = 1
        CALL  NUMLST(VECMVO(IFRP),WORK(KVEC(IFRP,1)),
     &               NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &               IFRP,NVEC(IFRP,0))
        CALL ORBCNT(WORK(KVEC(IFRP,1)),NVEC(IFRP,0),
     &              NPSH(IFRP),NESH(IFRP),
     &              NVEC(IFRP,2),NVEC(IFRP,1))
        NVEC(IFRP,0) = NVEC(IFRP,1) + NVEC(IFRP,2)
        CALL MEMREL('DENMVO',WORK,KFRSAV,KVEC(IFRP,1),KFREE,LFREE)
        CALL MEMGET2('INTE','ILIST',KVEC(IFRP,1),NVEC(IFRP,0),
     &              WORK,KFREE,LFREE)
      ENDDO
C
C     Calculate dimensions of the coefficient array
C
      NCMV = 0
      DO IFRP = 1,NFSYM
        ICMV(IFRP) = NCMV
        NCMV       = NCMV + NFBAS(IFRP,0)*NVEC(IFRP,0)*NZ
      ENDDO
C
C     Select the set that we need
C
      DO IFRP = 1, NFSYM
         CALL SELCFS (CMO(ICMOQ(IFRP)+1),IFRP,BUF(ICMV(IFRP)+1),
     &                NVEC(IFRP,0),
     &                WORK(KVEC(IFRP,1)),NVEC(IFRP,2),NVEC(IFRP,1),
     &                NFBAS(IFRP,0),NORB(IFRP))
      ENDDO
C
C     Get the density matrix
C
CTROND
      CALL DZERO(DMAT,N2BBASXQ)
      DO 10 I = 1,NFSYM
        IF(NVEC(I,0).EQ.0) GOTO 10
        CALL DENST1(DMAT(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,D1,D0,
     &              BUF(ICMV(I)+1),NFBAS(I,0),NVEC(I,0),
     &              1,NVEC(I,0),NFBAS(I,0))
   10 CONTINUE
C
C     Print section
C
      IF(IPRMVO.GE.4) THEN
        CALL TITLER('Output density matrix from DENMVO','*',103)
        DO 20 I = 1,NFSYM
        IF(NVEC(I,0).EQ.0) GOTO 20
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(DMAT(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   20   CONTINUE
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Dhftst */
      SUBROUTINE DHFTST(DSCFON,ICHKCNV,BMAT,IPRINT)
C***********************************************************************
C
C     Restart checks etc. for DHF-cycle
C     Based on old code.
C     Written by T.Saue Oct 1 1998
C    
C     Last modification: MI, July 2005 for BSSpostDC-SCF calculation
C
C***********************************************************************
      use num_grid_cfg

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dcbham.h"
C
      DIMENSION BMAT(MXDIIS*MXDIIS)
      LOGICAL   TOBE,OLDFMO,DSCFON

C     CALL QENTER('DHFTST')
C
      BARNUC  = .TRUE.
      CACC    = '        '
      DSCFON  = .FALSE.
      OLDFMO  = .FALSE.
      ICHKCNV = 0
      KITER   = 1

C
C     1. One-electron system
C     ======================
C
      IF(ONESYS) THEN
        TRIVEC = .FALSE.
        TRIFCK = .FALSE.
        NITER  = 1
        RETURN
      ENDIF
C
C
C     2. Many-electron system
C     =======================
C
C
C     ***************************************
C     **** Looking for old wave function ****
C     ***************************************
C
C
C     * ATOMST: Atomic start
C     ------------------------------------------------
C
      IF(ATOMST) THEN
        BARNUC = .FALSE.
      ENDIF
C
C     * TRIFCK: Looking for old 2-electron Fock matrix
C     ------------------------------------------------
C
      IF(TRIFCK.AND.NASHT.GT.0) THEN
         WRITE(LUPRI,'(/A/A/)') '*** WARNING *** '//
     &      'Trial 2-electron Fock matrix not possible for open shell.',
     &      'Looking for trial vectors instead'
         TRIFCK = .FALSE.
         TRIVEC = .TRUE.
      END IF
      IF(TRIFCK) THEN
        INQUIRE(FILE='DFFCK2',EXIST=TOBE)
        IF(TOBE) THEN
          DHF_INTTYP = 'Old FOCK2   '
          BARNUC = .FALSE.
          DOHUCKEL = .FALSE.
        ELSE
          WRITE(LUPRI,'(/A/)')
     &    '*** WARNING *** No trial 2-electron Fock matrix found.'//
     &    'Looking for trial vectors instead'
          TRIFCK = .FALSE.
          TRIVEC = .TRUE.
        ENDIF
      ENDIF
C
C     * TRIVEC: Looking for old coefficients
C     --------------------------------------
C
CMI/May05 Do not use old DFCOEF when DO4C2C=.TRUE. !
      IF(TRIVEC.AND..NOT.RESTFCK) THEN
        INQUIRE(FILE='DFCOEF',EXIST=TOBE)
        IF(TOBE) THEN
          BARNUC   = .FALSE.
          DOHUCKEL = .FALSE.
        ELSE
          WRITE(LUPRI,'(/A)') '*** INFO *** No trial vectors found.'
          IF(ATHUCK) THEN
            WRITE(LUPRI,'(A)') ' Using atomic Huckel start.'
          ELSEIF(BARNUC) THEN
            IF (BNCRON) WRITE(LUPRI,'(A/)')
     &  " Using bare nucleus plus electronic screening (Slater's"//
     &  " rules) as start potential."
            IF (BNSPON) WRITE(LUPRI,'(A/)')
     &    " Using sum of fitted atomic potentials as start potential.",
     &    " Reference: S. Lehtola, L. Visscher, E. Engel, "//
     &    " J. Chem. Phys. 152 (2020) 144105. (small GRASP fit)"
          ENDIF
          TRIVEC = .FALSE.
!         do not use density-based atomic size adjustments in this case
          num_grid_cfg_estimate_radii = .false.
        ENDIF
      ENDIF
C
C     ***********************************
C     **** Looking for cycle history ****
C     ***********************************
C
      INQUIRE(FILE='DFCYCL',EXIST=TOBE)
CMI   ... do always utilize DFCYCL when continuing with BSS-SCF after DC-SCF
      IF (.NOT.RESTFCK) TOBE = TOBE.AND.(TRIVEC.OR.TRIFCK)

      IF(TOBE) THEN
        OPEN(LUCYCL,FILE ='DFCYCL',STATUS='OLD',
     &              ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
        REWIND LUCYCL
 10     CONTINUE
        READ(LUCYCL,1020,END=20) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &              CACC,WALTID,DHF_INTTYP,DAYTID
        GOTO 10
 20     CONTINUE
#if defined (VAR_MFDS)
C        backspace on multifile systems for correct positioning
C        at the end of the file
        BACKSPACE LUCYCL
#endif
        IF (RESTFCK) THEN
         IF (DO4C2C) THEN
          WRITE(LUPRI,'(/2X,A/)')
     &    '...restarting of BSS-SCF after DC-SCF !'
         ELSE
          WRITE(LUPRI,'(/2X,A/)')
     &    '...restarting of DC-SCF after BSS-SCF !'
         ENDIF
        ENDIF

        WRITE(LUPRI,'(A)') '*** WARNING *** : SCFCYC restart info:'
        WRITE(LUPRI,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &            CACC,WALTID,DHF_INTTYP,DAYTID
        ILLINT = 0
        IF(DHF_INTTYP(1:2).EQ.'LL') ILLINT = 1
        ISLINT = 0
        IF(DHF_INTTYP(4:5).EQ.'SL') ISLINT = 1
        ISSINT = 0
        IF(DHF_INTTYP(7:8).EQ.'SS') ISSINT = 1
        IGTINT = 0
        IF(DHF_INTTYP(10:11).EQ.'GT') IGTINT = 1
        INTBUF = ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT
        NITER = NITER + 1
        CALL INTCON(INTFLG,INTBUF,INTDEF,
     &            ERGVAL,CNVINT,NITER,ITRINT,DHF_INTTYP)
C
C       * Check on convergence acceleration
C       -----------------------------------
C
        INQUIRE(FILE='DFFOCK',EXIST=OLDFMO)
        OLDFMO = OLDFMO.AND.(DODIIS.OR.DODAMP)
        IF(DODIIS) THEN
          DIISON = OLDFMO.AND.INTFLG.EQ.INTBUF.AND.NITER.GT.1
          INQUIRE(FILE='DFDIIS',EXIST=TOBE)
          DIISON = DIISON.AND.TOBE
          INQUIRE(FILE='DFCMOS',EXIST=TOBE)
          DIISON = DIISON.AND.TOBE
          INQUIRE(FILE='DFEVEC',EXIST=TOBE)
          DIISON = DIISON.AND.TOBE
          IF(DIISON) THEN
            CALL OPNFIL(LUDIIS,'DFDIIS','OLD','SCFINT')
            READ(LUDIIS) ITDIIS,MDIIS2,NELMBM,BMAT
            IF(MDIIS2.EQ.MXDIIS) THEN
              OLDFMO = .TRUE.
              WRITE(LUPRI,'(A)') '*** WARNING *** : DIIS RESTART'
              WRITE(LUPRI,'(3X,A,I5)')
     &        '* Previous DIIS iteration        :', ITDIIS,
     &        '* Number of elements in B-matrix :',NELMBM
            ELSE
              ITDIIS = 0
            ENDIF
          ELSE
            ITDIIS = 0
          ENDIF
        ENDIF
C
C       Check for differential density matrix
C
        IF(DODSCF.AND.(INTBUF.EQ.INTFLG)) THEN
          INQUIRE(FILE='DFFCK2',EXIST=DSCFON)
          INQUIRE(FILE='DFDENS',EXIST=TOBE)
          DSCFON = DSCFON.AND.TOBE
          IF(DSCFON) THEN
            WRITE(LUPRI,'(A)')
     &      '*** WARNING *** : Differential density restart'
          ENDIF
        ENDIF
      ELSE
        NITER  = 1
        OPEN(LUCYCL,FILE ='DFCYCL',STATUS='UNKNOWN',
     &     ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
        IF(TRIVEC.OR.ATOMST) CALL INTCON(INTFLG,INTBUF,INTDEF,
     &                  DUMMY,CNVINT,NITER,ITRINT,DHF_INTTYP)
      ENDIF
C
C     If (oldfmo) we can check all convergence criteria
C
      IF (OLDFMO) ICHKCNV = 1
      RETURN
 1010 FORMAT('It. ',I4,1P,G23.13,3D10.2,3X,A8,3X,A12,3X,A12,3X,A10)
 1020 FORMAT(4X,    I4,   G23.13,3D10.2,3X,A8,3X,A12,3X,A12,3X,A10)
C     ... 1020 is used for reading, thus no '1P' !!! /hjaaj
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck fndocc */
      subroutine fndocc(eig, skip_check_for_sym_broken_occ)
C***********************************************************************
C
C     Find DHF occupation based on eigenvalues
C     using the build-up principle.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dgroup.h"
C
      DIMENSION EIG(*)
      DIMENSION IEIGMIN(2)
      logical   symmetry_broken
      integer, save :: icharge       = 0
      integer       :: total_num_e, nish_save(2),nash_save(2)

!     radovan: the check below is skipped if fndocc is called from visual
!              this is to avoid modification of occupation if converged homo-lumo gap
!              is very small
      logical, intent(in) :: skip_check_for_sym_broken_occ
      logical             :: skip_check = .false.

      skip_check = skip_check_for_sym_broken_occ

!     (re-)set NELECT_DHF as it might have changed from last time 
!     because of the detection of symmetry broken solution...
!     -----------------------------------------------------------

!     a. get total number of electrons: sum_of_atomic_numbers - sum_of_charges
      if(icharge.lt.1)then
!       determine total # of electrons
        call rmolchr(icharge)
      end if
      total_num_e = icharge - kcharg

!     b. set nelect_dhf
      nelect_dhf = total_num_e
      naelec_dhf = 0
      if(mod(nelect_dhf,2).eq.1) then
        nelect_dhf = nelect_dhf - 1
        naelec_dhf = 1
      end if
      
!     initialize arrays
      DO I = 1,NFSYM
         nish_save(I)   = NISH(I)
         nash_save(I)   = NASH(I)
         IEIGMIN(I)     = IORB(I) + NPSH(I) + 1
         NISH(I)        = 0
         NASH(I)        = 0
         NACSH(I,1)     = 0
      END DO
C
      NELEC2GO = NELECT_DHF
      IF ((NELECT_DHF+1)/2 .GT. NESHT) THEN
         WRITE(LUPRI,'(I3,A,I3,A)')
     &      NELECT_DHF,' electrons but only ',NESHT,' orbitals.'
         CALL QUIT('ERROR: Not enough eigenvalues for FNDOCC.')
      END IF
!     no closed shell orbitals
      IF(NELEC2GO .eq. 0) GOTO 190
C
C     100 loop is over electron pairs (closed shell orbitals)
C
 100  CONTINUE
C        Find the irrep that has the lowest eigenvalue
         IEIGSYM          = 0
 110     IEIGSYM          = IEIGSYM + 1
         IF (NISH(IEIGSYM) .GE. NESH(IEIGSYM)) GOTO 110
         J                = IEIGSYM
!        stefan: add test if fermion irrep J+1 contains any orbitals
           DO I = J+1,NFSYM
            if(norb(I)+iorb(i).ge.IEIGMIN(I))then
               IF (EIG(IEIGMIN(I)) .LT. EIG(IEIGMIN(IEIGSYM)))
     &         IEIGSYM    = I
            end if
           END DO
C        We found the lowest yet unoccupied MO, occupy it
         NISH(IEIGSYM)    = NISH(IEIGSYM) + 1
         IEIGMIN(IEIGSYM) = IEIGMIN(IEIGSYM) + 1
         NELEC2GO         = NELEC2GO - 2
C        More orbitals to be filled ?
      IF (NELEC2GO .GT. 0) GOTO 100
C
C     Same procedure for the open shell (works only for one e in one open shell)
C     we take care of degeneracies in the following...
C
 190  IF (NAELEC_DHF .gt. 0) THEN
         IEIGSYM          = 0
 200     IEIGSYM          = IEIGSYM + 1
      IF (NISH(IEIGSYM) .GE. NESH(IEIGSYM)) GOTO 200
         J                = IEIGSYM
         DO I = J+1,NFSYM
            IF (EIG(IEIGMIN(I)) .LT. EIG(IEIGMIN(IEIGSYM)))
     &         IEIGSYM    = I
         END DO
         NASH(IEIGSYM)    = 1
         NACSH(IEIGSYM,1) = 1
         IEIGMIN(IEIGSYM) = IEIGMIN(IEIGSYM) + 1
      END IF
C
!     stefan + hjaaj - jan 2011:
!     check for symmetry-broken solutions "caused" by the Aufbau
!     principle and correct it.
!     stop if "fractional occupation" is not in use?
      symmetry_broken = .false.
      if (.not. skip_check) then
         if (.not. aoc) then
           call check_for_symmetry_broken_occupations(eig,nish,nash,
     &                             norb,norbt,ieigsym,nfsym,ieigmin,
     &                             total_num_e,symmetry_broken)
           if (symmetry_broken) then
           !  symmetry broken compared to aufbau, is it also symmetry broken
           !  compared to previous SCF iteration ???
              symmetry_broken = .false.
              do i = 1,nfsym
                 if (nish(i) .ne. nish_save(i)) symmetry_broken = .true.
                 if (nash(i) .ne. nash_save(i)) symmetry_broken = .true.
              end do
           end if
         end if
      end if
!
      NASHMFT = 0
      DO I = 1, NFSYM
         NOCC(I) = NISH(I) + NASH(I)
C
C        LV : CAUTION, POSSIBLE BUG
C
C        After the MCSCF merge we now have two arrays that seem to be
C        used interchangeably at various places, I just copy the info
C        to both so that the auto-occupation feature works for DHF but
C        this needs to be checked by the MCSCF bunch !!!
C
         NISH_DHF(I) = NISH(I)
         NASH_DHF(I) = NASH(I)
      END DO
C     Subtract frozen orbitals (if any) from orbital classes
      CALL SUBFRZ()

      if(symmetry_broken)then
        write(lupri,'(/a)') '  INFO from FNDOCC: following the'//
     &  ' Aufbau principle led to a symmetry_broken solution -'//
     &  ' the program corrected it for you.'
        write(lupri,'(a/)') '  INFO from FNDOCC: the new occupation'//
     &  ' will follow below:'
      end if

      END

C***********************************************************************
      subroutine check_for_symmetry_broken_occupations(eig,nish,nash,
     &                     norb,norbt,ieigsym_last_occ,nfsym,ieigmin,
     &                     total_num_e,symmetry_broken)
!***********************************************************************
!
!     purpose: detect symmetry-broken solutions when using the
!     Aufbau principle to set-up an (initial) electronic configuration.
!
!     In other words: detect if degenerate orbitals not all have same
!     occupation.
!
!     Redistribute electrons to not break symmetry by using
!     average-of-configurations (AOC), if needed.
!
!     written by s. knecht - jan 2011
!
!***********************************************************************
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
#include "dcbdhf.h"
      real(8), parameter   :: dthrs_null     = 1.0d-13
      real(8), parameter   :: degeneracy_thr = 6.0d-02
      real(8), parameter   :: d0             = 0.0d0
      real(8), parameter   :: dp5            = 0.5d0
      real(8), parameter   :: d1             = 1.0d0
      real(8), parameter   :: d2             = 2.0d0
      integer, parameter   :: max_degeneracy = 128

      real(8), intent(in)    :: eig(*)
      integer, intent(inout) :: nish(nfsym), nash(nfsym), norb(nfsym)
      integer, intent(inout) :: ieigmin(nfsym)
      logical, intent(inout) :: symmetry_broken
      integer, intent(in)    :: ieigsym_last_occ,total_num_e

      real(8), allocatable :: eig_deg(:)
      integer, allocatable :: ieig_deg(:)
      real(8)              :: eig_last
      integer              :: orboff, deg_counter, subtract_nish
      logical              :: debug_print

      debug_print      = .false.
!     debug_print      = .true.
      if(debug_print)then 
        orboff           = 0
        do i = 1, nfsym
          print *, 'eigenvalues for isym =',i,' # of orbs =',norb(i)
          call wrtmat(eig(1+orboff),1,norb(i),1,norb(i))
          orboff = orboff + norb(i)
        end do
      end if

!     initialize the last Kramers-pair in which electron(s) was (were) put 
      ilast           = ieigmin(ieigsym_last_occ) - 1

      naelec_dhf_save = naelec_dhf ! =0 if number of electrons is even, =1 if odd

      call alloc(ieig_deg,max_degeneracy,
     &           id='pointer to degenerate orbs')
      call alloc(eig_deg, max_degeneracy,
     &           id='eigenvalues of degenerate orbs')

 400  call izero(ieig_deg,max_degeneracy)
      call dzero(eig_deg, max_degeneracy)

      orboff   = 0
      do i = 1, ieigsym_last_occ-1
        orboff = orboff + norb(i)
      end do


      if(debug_print)then
        print *, '  naelec_dhf is',naelec_dhf
        print *, '  ieigsym_last_occ is',ieigsym_last_occ,
     &           '  the last e- was put into ',ilast
        print *, '  the eigenvalue is: ',eig(ilast)
      end if


      eig_last              = eig(ilast)
      deg_counter           = 1
      ieig_deg(deg_counter) = ilast
      eig_deg( deg_counter) = eig_last
      iup                   = ilast
      idown                 = ilast

!     search up and down for degenerate (or near-degenerate) kramers-pairs
!     a. up
 500  continue
        iup = iup + 1
        if(iup.le.orboff+norb(ieigsym_last_occ))then
          if(abs(eig(iup)-eig_last).le.degeneracy_thr)then
            deg_counter           = deg_counter + 1
            ieig_deg(deg_counter) = iup
            eig_deg( deg_counter) = eig(iup)
            goto 500
          end if
        end if
!       b. down
 600  continue
        idown = idown - 1
        if(idown.gt.orboff)then
          if(abs(eig(idown)-eig_last).le.degeneracy_thr)then
            deg_counter           = deg_counter + 1
            ieig_deg(deg_counter) = idown
            eig_deg( deg_counter) = eig(idown)
            goto 600
          end if
        end if

!     set new # of active e- according to degeneracies
!     ==> (ilast - idown ) * 2 - naelec_dhf_save
!     naelec_dhf_save: 1 for odd numbers of total e-, 0 otherwise
      naelec_dhf =  2*(ilast-idown) - naelec_dhf_save

      if(debug_print)then
        print *, '  active e- initial, revised (for re-distribution):',
     &              naelec_dhf_save, naelec_dhf
      end if
      d_nr_acte_revised = real(naelec_dhf,8)

      dalpha(0) = d0
      df(0)     = d1
      da(0)     = d1
      nashmft   = 0

      nopen     = 1

      if(debug_print)then
        print *,'  found degeneracy in irrep:', ieigsym_last_occ
        do i = 1,deg_counter
          print *,'  orbital #:',ieig_deg(i), '; eigenvalue = ',
     &                                           eig_deg(i)
        end do
      end if

!     start revised automatic occupation taking into account degeneracy

!     case a. number of e- == 2 * number of degenerate orbs; the "closed-shell case"
      if(naelec_dhf .eq. 2*deg_counter) then 
        nopen      = 0
        naelec_dhf = 0
        nasht_dhf  = 0
        nelect_dhf = 0
        call izero(nash_dhf,nfsym)
        call izero(nash    ,nfsym)
        call icopy(nfsym,nish,1,nish_dhf,1)
        do ifrp = 1,nfsym
          nelec_dhf(ifrp) = nish_dhf(ifrp)*2
          nelect_dhf      = nelect_dhf + nelec_dhf(ifrp)
        end do
!       consistency check
        if(nelect_dhf.ne.total_num_e)then
          print *, 'calculated # of e-: ',nelect_dhf
          print *, 'total # of e-     : ',total_num_e
          call quit(' *** error in calculated # of e-')
        end if

      else ! naelec_dhf != 2*deg_counter

!       case b. number of e- < 2*deg_counter and # of degenerate orbs > 1
!       --> p-shell, d-shell, f-shell, or near-degenerate MOs...
        if(deg_counter .gt. 1)then

          if(debug_print)then
            print *,' number of active e- for redistribution ==',
     &                naelec_dhf
          end if
!      
          call izero(nacsh,nfsym*mxopen)
          nacsh(ieigsym_last_occ,1) = deg_counter
          call icopy(nfsym*mxopen,nacsh,1,nacshmf,1)
  
          df(1)             = d_nr_acte_revised/real(deg_counter*2,8)
          dm                = d2*(real(nacsh(1,1)+nacsh(2,1),8))
          dn                = dm*df(1)
          da(1)             = (dm*(dn-d1))/(dn*(dm-d1))
          if(abs(da(1)).le.dthrs_null) da(1) = d0
          dalpha(1)         = (d1-da(1))/(d1-df(1))

          if(debug_print)then
            print *,'dm,dn,df(iopen),da(iopen),dalpha(iopen)',
     &               dm,dn,df(1),da(1),dalpha(1)
            print *,'nish(1), nish(2), nish_dhf(1), nish_dhf(2)',
     &               nish(1), nish(2), nish_dhf(1), nish_dhf(2)
            print *,' ilast, idown =',ilast, idown
            print *,' naelec_dhf_save, naelec_dhf =',
     &                naelec_dhf_save, naelec_dhf
          end if

!         reset nish array - take into account whether an unpaired e-
!         had been placed in the last orb.
          subtract_nish          = ilast - idown
          if(naelec_dhf_save.gt.0) subtract_nish = subtract_nish - 1

          if(debug_print)then
            print *, 'nish(ieigsym_last...), subtract_nish, naelec_dhf',
     &                nish(ieigsym_last_occ),subtract_nish, naelec_dhf
          end if

          nish(ieigsym_last_occ) = nish(ieigsym_last_occ)-subtract_nish
          call icopy(nfsym,nish,1,nish_dhf,1)

          call izero(nash_dhf,nfsym)
          call izero(nash,    nfsym)
          nasht_dhf  = 0
          nelect_dhf = 0
          do ifrp = 1,nfsym
            do iopen = 1,nopen
              nash_dhf(ifrp) = nash_dhf(ifrp) + nacsh(ifrp,iopen)
              nashmft        = nashmft + nacshmf(ifrp,iopen)
            end do
            nelec_dhf(ifrp) = nish_dhf(ifrp)*2
            nasht_dhf       = nasht_dhf  + nash_dhf(ifrp)
            nelect_dhf      = nelect_dhf + nelec_dhf(ifrp)
          end do
          call icopy(nfsym,nash_dhf,1,nash,1)

!         consistency check
          if(nelect_dhf+naelec_dhf.ne.total_num_e)then
            print *, 'calculated # of e-: ',nelect_dhf+naelec_dhf
            print *, 'total # of e-     : ',total_num_e
            call quit(' *** error in calculated # of e-')
          end if

!         initial occupation from the Aufbau principle ended up with a
!         symmetry broken solution - therefore set flag to .true.
          symmetry_broken  = .true.

!       case c. number of e- < 2*deg_counter and # of degenerate orbs == 1
!       --> s-shell or a "lonely" MO...
        else if(deg_counter.eq.1)then

!         consistency check
          if(naelec_dhf.ne.1)then
            print *, 'calculated # of active e- : ',naelec_dhf
            print *, 'total # of degenerate orbs: ',deg_counter
            call quit(' *** error in calculated # of active e-')
          end if

          nasht_dhf  = naelec_dhf
          nelect_dhf = 0
          call izero(nacsh,nfsym*mxopen)
          call izero(nash_dhf,nfsym)
          call izero(nash,    nfsym)

          nacsh(ieigsym_last_occ,1) = deg_counter

          df(1)             = d_nr_acte_revised/real(deg_counter*2,8)
          dm                = d2*(real(nacsh(1,1)+nacsh(2,1),8))
          dn                = dm*df(1)
          da(1)             = (dm*(dn-d1))/(dn*(dm-d1))
          if(abs(da(1)).le.dthrs_null) da(1) = d0
          dalpha(1)         = (d1-da(1))/(d1-df(1))

          call icopy(nfsym,nish,1,nish_dhf,1)
          do ifrp = 1,nfsym
            nash_dhf(ifrp)  = nash_dhf(ifrp) + nacsh(ifrp,1)
            nashmft         = nashmft + nacshmf(ifrp,1)
            nelec_dhf(ifrp) = nish_dhf(ifrp)*2
            nelect_dhf      = nelect_dhf + nelec_dhf(ifrp)
          end do
          call icopy(nfsym,nash_dhf,1,nash,1)

        end if ! deg_counter > 1

        if(debug_print)then
          print *,'nash(1), nash(2), nash_dhf(1), nash_dhf(2)',
     &             nash(1), nash(2), nash_dhf(1), nash_dhf(2)
          print *,'nelect_dhf - end:',nelect_dhf
          print *,'nelec_dhf(1), nelec_dhf(2) - end:',
     &             nelec_dhf(1), nelec_dhf(2)         
          print *,'nish(1), nish(2), nish_dhf(1), nish_dhf(2)',
     &             nish(1), nish(2), nish_dhf(1), nish_dhf(2)  
        end if

      end if ! naelec_dhf .eq. 2*deg_counter ?

      call dealloc( eig_deg)
      call dealloc(ieig_deg)

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dfdiag */
      SUBROUTINE DFDIAG(FMO,EIG,IBEIG,FAO,DODBG,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Diagonalize Dirac-Fock matrix
C
C     On input:
C       FMO - Fock (or other) matrix in MO basis to be diagonalized 
C       DODBG - decide to print out control output
C
C     On output:
C       FAO - eigenvectors of the diagonalized FMO matrix
C       EIG - eigenvalues
C       IBEIG - boson irreps of eigenvectors
C
C     Taken out of DFSOLV and modified by L.Visscher Jan 10 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      PARAMETER (D0 = 0.00D00, D1 = 1.00D00)
      LOGICAL   SPECIAL,DODBG
      DIMENSION FMO(*),EIG(*),IBEIG(*),FAO(*),WORK(*)

      CALL QENTER('DFDIAG')
C
C     Check for special cases :
C     - Spinfree or Levy-Leblond; diagonalize boson irreps (subblocks)
C     - Linear symmetry; diagonalize Omega irreps (subblocks)
C     - ZORA; first diagonalize SS-block, then form and diagonalize
C       modified large component block
C
      SPECIAL = (ZORA.OR.SUB_BL)
      CALL DZERO(FAO,N2BBASXQ)
      CALL ICOPY (NORBT,0,0,IBEIG,1)
      NBRP = 4 / NZ
C
      DO 30 I = 1,NFSYM
C
        IF(NTMO(I).EQ.0) GOTO 30
         iprham_save = iprham
!        iprham = 5
C ... control print out of the matrix to be diagonalized
         IF (DODBG.AND.IPRHAM.GE.5) THEN
          CALL HEADER(
     &   'DFDIAG: Ferm.symm.blocked MO matrix to be diagonalized',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &    '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(FMO(1+I2TMOT(I)),
     &    NTMO(I),NTMO(I),NTMO(I),NTMO(I),
     &    NZ,IPQTOQ(1,0),LUPRI)
         ENDIF 
C
C       Diagonalize Fock-matrix in MO-basis
C         - store coefficients in FAO
C       ===================================
C
C       1 : Normal diagonalization, use matrix as supplied by dfsolv
C
        IF (.NOT.SPECIAL) THEN

          IF(DOJACO) THEN
            CALL RSJACO(NTMO(I),NTMO(I),NTMO(I),FMO(I2TMOT(I)+1),
     &                  EIG(IORB(I)+1),1,1,0,FAO(I2TMOT(I)+1))
          ELSE
            CALL QDIAG(NZ,NTMO(I),FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                 EIG(IORB(I)+1),1,FAO(I2TMOT(I)+1),NTMO(I),
     &                 NTMO(I),WORK(KFREE),LFREE,IERR)
          ENDIF
C
#ifdef MOD_UNRELEASED
          IF(VACPOL) THEN
            CALL HEADER(
     &      'DFDIAG: Eigenvalues (normal diagonalization) :',-1)
            WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,'/',NFSYM
            WRITE(LUPRI,'(I5,F25.10)')
     &      (J,EIG(IORB(I)+J),J=1,NTMO(I)) 
         ENDIF
#endif
         IF (DODBG) THEN
          IF (IPRHAM.GE.3) THEN
           CALL HEADER(
     &     'DFDIAG: Eigenvalues (normal diagonalization) :',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,'/',NFSYM
           WRITE(LUPRI,'(I5,F25.10)')
     &     (J,EIG(IORB(I)+J),J=1,NTMO(I)) 
          ENDIF
C
          IF (IPRHAM.GE.5) THEN
           CALL HEADER(
     &     'DFDIAG: Eigenvectors :',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(FAO(1+I2TMOT(I)),
     &     NTMO(I),NTMO(I),NTMO(I),NTMO(I),
     &     NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
         ENDIF 
C
C
C       2 : Subblock diagonalization (e.g. Levy-Leblond, spinfree, linear symmetry)
C           use matrix as supplied by dfsolv but diagonalize blockwise
C
!       IF (.NOT.SPECIAL) THEN
        ELSEIF (.NOT.ZORA) THEN
C
C         Spinfree calculations can use the boson blocking of the
C         Fock matrix. For quaternion groups we use the fact that
C         the matrix becomes real. For linear symmetry we also have
C         real matrices, but a blocking given by another array.
C         Initialize FAO to zero.
C         --------------------------------------------------------
C
          CALL DZERO(FAO(1+I2TMOT(I)),N2TMO(I)*NZ)

          IBO = 0
          DO ISUB = 1, N_SUB_BL(I)
            NBO = NTMO_SUB(ISUB,I,0)
            IF(NBO.GT.0) THEN
              IMAT = I2TMOT(I) + (NTMO(I)+1)*IBO + 1
              IEIG = IORB(I) + IBO + 1
              IF (DODBG.AND.IPRHAM.GE.5) THEN
                CALL HEADER(
     &          'DFDIAG: MO matrix (subblocks) to be diagonalized',-1)
                WRITE(LUPRI,'(3X,A,I1,A,I1,3X,A,I3,A,I3/)')
     &          '*** Fermion corep ',I,'/',NFSYM,
     &          '*** boson sub-corep ',ISUB,'/',N_SUB_BL(I)
                CALL PRQMAT(FMO(IMAT),NBO,NBO,NTMO(I),NTMO(I),
     &                      1,IPQTOQ(1,0),LUPRI)
              ENDIF
              IF (WRITE_FMO_MATRIX) CALL FMO_WRITE(FMO,IMAT,I,ISUB,NBO)
C
C             Diagonalize the isub'th block of parity i
C
#ifdef MOD_UNRELEASED
              IF (.NOT.DOQJACO) THEN
#endif
                IF (DOJACO) THEN
                  CALL RSJACO(NTMO(I),NBO,NBO,FMO(IMAT),
     &                      EIG(IEIG),1,1,0,FAO(IMAT))
                ELSE
                  CALL QDIAG(1,NBO,FMO(IMAT),NTMO(I),NTMO(I),
     &                   EIG(IEIG),1,FAO(IMAT),NTMO(I),NTMO(I),
     &                   WORK(KFREE),LFREE,IERR)
                ENDIF
#ifdef MOD_UNRELEASED
              ELSE
                CALL MEMGET2('REAL','FMO',KFMO,NBO*NBO,WORK,KFREE,LFREE)
                CALL MEMGET2('REAL','FA ',KFA, NBO*NBO,WORK,KFREE,LFREE)
                IPOS=0
                DO II=1,NBO
                  DO JJ=1,NBO
                    WORK(KFMO+IPOS)=FMO(IMAT+(NTMO(I)*(II-1))+JJ-1)
                    IPOS=IPOS+1
                  ENDDO
                ENDDO
                CALL QJACOBI(WORK(KFMO),WORK(KFA),NBO,1,0,
     &                       IDUMMY,.TRUE.,IPRSCF)
                IPOS=0
                DO II=1,NBO
                  EIG(IEIG+II-1)=WORK(KFMO+(NBO*(II-1)+II-1))
                  DO JJ=1,NBO ! store eigenvectors
                    FAO( IMAT+ (NTMO(I)*(II-1)+JJ-1) ) =
     &              WORK(KFA+(NBO*(II-1))+JJ-1)
                    IPOS=IPOS+1
                  ENDDO
                ENDDO
                CALL MEMREL('DFDIAG',WORK,KFMO,KFMO,KFREE,LFREE)
              ENDIF
#endif

              IF (DODBG) THEN
                IF (IPRHAM.GE.3) THEN
                  CALL HEADER(
     &             'DFDIAG: Eigenvalues (special_subbblock) :',-1)
                   WRITE(LUPRI,'(3X,A,I1,A,I1,3X,A,I3,A,I3/)')
     &             '*** Fermion corep ',I,'/',NFSYM,
     &             '*** boson sub-corep ',ISUB,'/',N_SUB_BL(I)
                   WRITE(LUPRI,'(I5,F25.10)')
     &             (IX,EIG(IX),IX=IEIG,IEIG+NBO-1) 
                ENDIF
                IF (IPRHAM.GE.5) THEN
                  CALL HEADER(
     &            'DFDIAG: Eigenvectors(special_subblock)',-1)
                  WRITE(LUPRI,'(3X,A,I1,A,I1,3X,A,I3,A,I3/)')
     &            '*** Fermion corep ',I,'/',NFSYM,
     &            '*** boson sub-corep ',ISUB,'/',N_SUB_BL(I)
                  CALL PRQMAT(FAO(IMAT), NBO,NBO,NTMO(I),NTMO(I),
     &                        1,IPQTOQ(1,0),LUPRI)
                ENDIF
              ENDIF
              IBO = IBO + NTMO_SUB(ISUB,I,0)
            ENDIF ! IF(NBO.GT.0) THEN
          ENDDO
C
C         Sort the eigenvalues according to energy : keep track of
C         the block identification. (Use FMO as scratch space for the
C         vector sort)
C
          CALL INIBOS(IBEIG(IORB(I)+1),I,.FALSE.,IPRHAM)
          CALL MEMGET2('INTE','INDX',KINDX,NORB(I),WORK,KFREE,LFREE)
          CALL MEMGET2('INTE','IDUM',KIDUM,NORB(I),WORK,KFREE,LFREE)
          CALL SELBOS (I,NTMO(I),WORK(KINDX),EIG(IORB(I)+1),
     &                 IBEIG(IORB(I)+1),FAO(I2TMOT(I)+1),
     &                 FMO(I2TMOT(I)+1),WORK(KIDUM),IPRHAM)
          CALL MEMREL('DFDIAG',WORK,KINDX,KINDX,KFREE,LFREE)

          IF (DODBG.AND.IPRHAM.GE.3) THEN
           write(lupri,'(/,2x,a)') 'DFDIAG(subbl) - IBEIG values'
           write(lupri,'(2x,a,i3,a)') 'IBEIG sym',I,' >: '
           do K=1,NTMO(I)
            write(LUPRI,'(2X,I3,A,I3)') K,'-',IBEIG(IORB(I)+K)
           enddo
          ENDIF
C
C       3 : ZORA; diagonalize modified LL block
C
        ELSEIF (ZORA.AND.(.NOT.SUB_BL)) THEN
          NL = NESH(I)
          NS = NPSH(I)
          NT = NORB(I)
          IEIGP = IORB(I) + 1
          IEIGE = IEIGP + NS
          CALL ZORAEQ(I,0,NZ,NL,NS,NT,EIG,IEIGP,IEIGE,
     &                FMO,FAO,WORK,KFREE,LFREE)
C
C       4 : Spinfree ZORA : Make modified matrix but do it block by block
C
        ELSEIF (ZORA.AND.SPINFR) THEN
          NT = NORB(I)
          IEIGP = IORB(I) + 1
          DO 21 ISUB = 1, N_SUB_BL(I)
            NBO = NTMO_SUB(ISUB,I,0)
            IF (NBO.EQ.0) GOTO 21
            NL = NTMO_SUB(ISUB,I,1)
            NS = NTMO_SUB(ISUB,I,2)
            NZS = 1
            IEIGE = IEIGP + NS
            CALL ZORAEQ(I,ISUB,NZS,NL,NS,NT,EIG,IEIGP,IEIGE,
     &                  FMO,FAO,WORK,KFREE,LFREE)
            IEIGP = IEIGP + NBO
  21      CONTINUE
C
C           Sort the eigenvalues according to energy : keep track of
C           the boson symmetry. (Use FMO as scratch space for the
C           vector sort)
C
          CALL INIBOS(IBEIG(IORB(I)+1),I,.FALSE.,IPRHAM)
          CALL MEMGET2('INTE','INDX',KINDX,NORB(I),WORK,KFREE,LFREE)
          CALL MEMGET2('INTE','IDUM',KIDUM,NORB(I),WORK,KFREE,LFREE)
          CALL SELBOS (I,NORB(I),WORK(KINDX),EIG(IORB(I)+1),
     &                   IBEIG(IORB(I)+1),FAO(I2ORBT(I)+1),
     &                   FMO(I2ORBT(I)+1),WORK(KIDUM),IPRHAM)
          CALL MEMREL('DFDIAG',WORK,KINDX,KINDX,KFREE,LFREE)
C
C       End of different diagonalization cases
C
        ENDIF    !  IF (.NOT.SPECIAL) THEN ...
C
C       Check if any positronic eigenvalues are below -2mc^2
C
        IF(IPRSCF.GE.2) CALL EIGCHK(EIG(IORB(I)+1),NPSH(I))
        iprham = iprham_save
C
 30   CONTINUE
C
      CALL QEXIT('DFDIAG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck zoradg */
      SUBROUTINE ZORADG(NZS,NL,NS,FLL,FLS,FSS,VLL,VSL,VSS,
     &                  EIGE,EIGP,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Diagonalize the Dirac matrix using the ZORA approximation :
C     first fold in the small components and build an effective large
C     component equation.
C
C     Consider the general (modified) Dirac equation:
C
C     [ FLL FLS ] [C(L)]   [ SLL   0   ] [C(L)]
C     [         ] [    ] = [           ] [    ] * E
C     [ FSL FSS ] [C(S)]   [   0   SSS ] [C(S)]
C
C     The ZORA approximation is obtained by setting the SSS block
C     of the metric to zero and then performing an unnormalized
C     elimination of the small components (UESC).
C
C     C(S) = {-FSS}^(-1/2) * FSL * C(L)
C
C     This leads to the following equation for the large components
C
C     [ FLL + FLS*{-FSS}^(-1)*FSL ] C(L) = SLL * C(S) * E
C
C     On input:
C       FLL - LL block of MO Fock matrix
C       FLS - SL block of MO Fock matrix
C       FSS - SS block of MO Fock matrix
C
C     On output:
C       VLL - large component ZORA positive-energy solutions     
C       VSL - small component ZORA positive-energy solutions     
C       VSS - eigenvectors of F(SS)
C
C     Reference:
C       L. Visscher and T. Saue, J. Chem. Phys. 113(2000) 3996
C       "Approximate relativistic electronic structure methods based on 
C       the quaternion modified Dirac equation"
C
C     Taken out of DFDIAG and modified by L.Visscher Aug 28 2000
C     Polish Dec 5 2006 T. Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      PARAMETER (D0 = 0.00D00, D1 = 1.00D00)
      DIMENSION EIGP(NS),EIGE(NL),WORK(*)
      DIMENSION FLL(NL,NL,NZS),FSS(NS,NS,NZS),FLS(NL,NS,NZS)
      DIMENSION VLL(NL,NL,NZS),VSS(NS,NS,NZS),VSL(NS,NL,NZS)
C
C     Create ZORA modified kinetic energy operator
C     via diagonalisation of F(SS)
C     ============================================
C
C     Let U be the matrix that diagonalizes F(SS):
C
C        U^(dagger) F U = f
C
C.....diagonalize F(SS) 
C
      CALL QDIAG(NZS,NS,FSS,NS,NS,EIGP,1,VSS,NS,NS,
     &           WORK(KFREE),LFREE,IERR)
C
C.....form WSS = U * f^(-1/2)
C
      CALL DZERO (FSS,NZS*NS*NS)
      DO IZ = 1, NZS
        DO J = 1, NS
          EIGINV = D1 / SQRT(-EIGP(J))
          CALL DAXPY(NS,EIGINV,VSS(1,J,IZ),1,FSS(1,J,IZ),1)
        ENDDO
      ENDDO
C
C     Form intermediate WLS = FLS*WSS
C     ===============================
C
      CALL QGEMM(NL,NS,NS,D1,
     &     'N','N',IPQTOQ(1,0),FLS,NL,NS,NZS,
     &     'N','N',IPQTOQ(1,0),FSS,NS,NS,NZS,
     &      D0,IPQTOQ(1,0),VSL,NL,NS,NZS)
C
C     Form the ZORA modified K.E. matrix
C     WLL = F(LS)*{-F(SS)}^(-1)*F(SL)
C     and add to FLL
C     ======================================
C
      CALL QGEMM(NL,NL,NS,D1,
     &     'N','N',IPQTOQ(1,0),VSL,NL,NS,NZS,
     &     'H','N',IPQTOQ(1,0),VSL,NL,NS,NZS,
     &      D1,IPQTOQ(1,0),FLL,NL,NL,NZS)
C
C     Diagonalize the ZORA matrix
C     ===========================
C
      CALL QDIAG(NZS,NL,FLL,NL,NL,
     &           EIGE,1,VLL,NL,NL,
     &           WORK(KFREE),LFREE,IERR)
C
C     Form the coefficients for the small components:
C     ===============================================
C
C     C(S) = {-FSS}^(-1/2) * FSL * C(L)
C
C     to be stored in VSL.
C
      CALL QGEMM(NS,NL,NS,D1,
     &     'N','N',IPQTOQ(1,0),FSS,NS,NS,NZS,
     &     'H','N',IPQTOQ(1,0),VSL,NL,NS,NZS,
     &      D0,IPQTOQ(1,0),FLS,NS,NL,NZS)
C
      CALL QGEMM(NS,NL,NL,D1,
     &     'N','N',IPQTOQ(1,0),FLS,NS,NL,NZS,
     &     'N','N',IPQTOQ(1,0),VLL,NL,NL,NZS,
     &      D0,IPQTOQ(1,0),VSL,NS,NL,NZS)
C
C     Determine scaling factors
C
      CALL QGEMM(NL,NL,NS,D1,
     &     'H','N',IPQTOQ(1,0),VSL,NS,NL,NZS,
     &     'N','N',IPQTOQ(1,0),VSL,NS,NL,NZS,
     &      D0,IPQTOQ(1,0),FLL,NL,NL,NZS)
C
C     Normalize the vectors for ZORA-4
C
      IF (ZORA4)  THEN
         DO J = 1, NL
            FACTOR = D1 / SQRT(D1+FLL(J,J,1))
            DO IZ = 1, NZS
               CALL DSCAL (NL,FACTOR,VLL(1,J,IZ),1)
               CALL DSCAL (NS,FACTOR,VSL(1,J,IZ),1)
            ENDDO
         ENDDO
      ENDIF
C
C     Calculate the scaled ZORA eigenvalues
C
      IF (ZORASC)  THEN
         DO J = 1, NL
            FACTOR = D1 / (D1 + FLL(J,J,1))
            EIGE(J) = EIGE(J) * FACTOR
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck inibos */
      SUBROUTINE INIBOS(IBOS,IFRP,SORTLS,IPRINT)
C***********************************************************************
C
C     Fill an array IBOS(*) with boson irrep identification
C     for a given irep IFRP
C
C     SORTLS - flag for the S_L resorting of boson irreps
C
C     A Saue cut and paste job, thanks to Luuk ...April 2003
C     Generalized to general subblocks
C
C     last revison: june, 2005 M.ILIAS
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION IBOS(*)
      LOGICAL SORTLS
C
      NSUB = N_SUB_BL(IFRP)
      IX = 1
      IF (.NOT.SORTLS) THEN
        DO ISUB = 1, NSUB
          ID = ID_SUB_BL(ISUB,IFRP)
          DO IC = 2, 1, -1
            NBO = NTMO_SUB(ISUB,IFRP,IC)
            IF(NBO.GT.0) THEN
              CALL ICOPY(NBO,ID,0,IBOS(IX),1)
              IX = IX + NBO
            ENDIF
          ENDDO
        ENDDO
      ELSE
        DO IC = 2, 1, -1
          DO ISUB = 1, NSUB
            ID = ID_SUB_BL(ISUB,IFRP)
            NBO = NTMO_SUB(ISUB,IFRP,IC)
            IF(NBO.GT.0) THEN
              CALL ICOPY(NBO,ID,0,IBOS(IX),1)
              IX = IX + NBO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
   
      IF (IPRINT.GE.4) THEN
       CALL HEADER('*** Output from INIBOS ***',-1)
       write(lupri,'(2X,A,L1)') 'flag on S_L resorting, SORTLS=',SORTLS
       write(lupri,'(2X,A,I1,A,I3,A,I3,A,I3)')
     & 'IBOS for IFRP=',IFRP,' NTMO(IFRP)=',NTMO(IFRP),
     & ' positr.orbs:',NPSH(IFRP),' electr.orbs:',NESH(IFRP)
       write(lupri,'(2X,A)')
     &      'boson irreps of each MO (from pos. to el.):'
       write(lupri,'(4X,25I3)') (IBOS(J),J=1,NTMO(IFRP))
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck detocc */
      SUBROUTINE DETOCC(EIG,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Handle automatic occupation
C     A cut-and-paste job by T. Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
      LOGICAL TEST
      DIMENSION EIG(*),WORK(*)
      DIMENSION NISHOLD(2),NASHOLD(2)
C
C     Save old occupation
C
      CALL ICOPY(2,NISH,1,NISHOLD,1)
      CALL ICOPY(2,NASH,1,NASHOLD,1)
C
C     Determine new occupation
C
      CALL FNDOCC(EIG, .false.)
C
      TEST = .FALSE.
      DO IFRP = 1,NFSYM
        TEST = TEST.OR.(NISH(IFRP).NE.NISHOLD(IFRP))
        TEST = TEST.OR.(NASH(IFRP).NE.NASHOLD(IFRP))
      END DO
      IF (TEST.OR.INIOCC) THEN
C
C        AUTOCC has changed occupation...
C        Update common block info about N2ISHT,...
C
         CALL SETDC2(0)
         CALL SETDHF(0)
         IF(INIOCC) THEN
           INIOCC = .FALSE.
           WRITE(LUPRI,'(/A/)') '* AUTOCC( 0) : Initial occupation:'
         ELSE
           NEWOCC = NEWOCC + 1
           WRITE(LUPRI,'(/A,I2,A/)')
     &            '* AUTOCC(',NEWOCC,') : New occupation:'
         END IF
C
         CALL PROCC()
C
         IF(NEWOCC .GT. 10) THEN
           WRITE(LUPRI,'(/A,/A,I2,A)') 'Problems with AUTOCC:',
     &            'DHF occupation has changed ',NEWOCC, ' times.'
           CALL QUIT('DHF occupation has changed too many times.')
         END IF
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Subfrz */
      SUBROUTINE SUBFRZ()
C***********************************************************************
C
C     For each orbital class determine the number of frozen orbitals
C
C     Written by T. Saue Apr 16 2004
C
C***********************************************************************

         use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dgroup.h"
      integer, allocatable :: ind(:)
      integer, allocatable :: ibf(:)
C
      IF(NFROT.GT.0) CALL OPNFIL(LUTMAT,'DFFROZ','OLD','SUBFRZ')
      DO IFRP = 1,NFSYM
        NISHMF(IFRP) = NISH(IFRP)
        NOCCMF(IFRP) = NOCC(IFRP)
        DO IOPEN = 1,NOPEN
          NACSHMF(IFRP,IOPEN) = NACSH(IFRP,IOPEN)
        END DO
        IF(NFRO(IFRP).GT.0) THEN
          call alloc(ind, nfro(ifrp))
          call alloc(ibf, norb(ifrp))
          CALL ORBARR(IFRP,IBF)
          CALL READI(LUTMAT,NFRO(IFRP),IND)
          READ(LUTMAT)
          CALL SUBFR1(IFRP,IND,IBF)
          call dealloc(ind)
          call dealloc(ibf)
        ENDIF
      ENDDO
      IF(NFROT.GT.0) CLOSE(LUTMAT,STATUS='KEEP')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Subfr1 */
      SUBROUTINE SUBFR1(IFRP,IND,IBUF)
C***********************************************************************
C
C     For each orbital class determine the number of frozen orbitals
C
C     Written by T. Saue Apr 16 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION IND(*),IBUF(*)
C
      DO I = 1,NFRO(IFRP)
        IP = IBUF(IND(I))
        IF(IP.EQ.0) THEN
          NISHMF(IFRP) = NISHMF(IFRP) - 1
        ELSEIF(IP.GT.0) THEN
          NACSHMF(IFRP,IP) = NACSHMF(IFRP,IP) - 1
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck orbstr */
      SUBROUTINE ORBARR(IFRP,IBUF)
C***********************************************************************
C
C     Make an array to distinguish orbital classes
C     
C     Written by T. Saue Apr 16 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION IBUF(*)
C.....Negative-energy orbitals:
      II = 1
      CALL ICOPY(NPSH(IFRP),-2,0,IBUF(II),1)
C.....Inactive orbitals:
      II = II + NPSH(IFRP)
      CALL ICOPY(NISH(IFRP),0,0,IBUF(II),1)
C.....Active orbitals:
      II = II + NISH(IFRP)
      DO IOPEN = 1,NOPEN
        CALL ICOPY(NACSH(IFRP,IOPEN),IOPEN,0,IBUF(II),1)
        II = II + NACSH(IFRP,IOPEN)
      ENDDO
C.....Secondary positive-energy orbitals:
      CALL ICOPY(NSSH(IFRP),-1,0,IBUF(II),1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Levshi */
      SUBROUTINE LEVSHI(FMO,CMO)
C***********************************************************************
C
C     Level shift
C
C     Written by T. Saue July 4 2006
C
C  MI July 6,2006: Requires DFCOEF with written MO's
!
!     April 2013
!     sknecht: attempt to change level shift for closed-shell cases: i do not understand 
!              the old implementation below but apparently it is ok. 
!              my "intuitive solution" via shift_virtuals that is according
!              to the paper cited below is wrong. 
!              the diagonal elements of the Fock matrix in MO basis do NOT
!              seem to roughly correspond to the orbital energies even when damping/diis/differential density matrix 
!              options are turned off... can anyone explain it to me what is wrong with my idea/approach? 
!              the way i intend to do it is exactly the same as it is done+working 
!              in Turbomole (my reference check) and ReSpect (where i implemented it)
!              for further information check: 
!              V. R. Saunders and I. H. Hillier, 
!              A "Level-Shifting" method for converging closed shell Hartree-Fock wave functions, IJQC, 7, p699-705 (1973).
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D1 = 1.0D0,D0=0.0D0)
#include "dcbdhf.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION FMO(*),CMO(*)
      LOGICAL IS_HERE,IS_OPEN

      CALL QENTER('LEVSHI')
C
      IF (IPRSCF.GE.10) THEN
        INQUIRE(FILE='DFCOEF',EXIST=IS_HERE,OPENED=IS_OPEN)
        write(lupri,*) 'LEVSHI: DFCOEF IS_HERE=',IS_HERE
        write(lupri,*) 'LEVSHI: DFCOEF IS_OPEN=',IS_OPEN
      ENDIF

      REWIND LUCMOS
      CALL READT(LUCMOS,N2TMOTQ,CMO)
      DO I = 1,NFSYM
        ISTART = NPSH(I)+NISHMF(I)+1
        DO IOPEN = 1,NOPEN
          IF(OLEV(IOPEN).NE.D0) THEN
            CALL DENST1(FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,
     &              OLEV(IOPEN),D1,
     &              CMO(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &              ISTART,NACSHMF(I,IOPEN),NTMO(I))
          ENDIF
          ISTART = ISTART + NACSHMF(I,IOPEN)
        ENDDO
        NVEC = NESH(I) - NOCCMF(I)
        IF(DLSHIF.GT.D0) THEN

!         CALL HEADER(
!    &   'shift virt:Ferm.symm.blocked MO matrix to be diagonalized',-1)
!         WRITE(LUPRI,'(3X,A,I1,A,I1/)')
!    &    '*** Fermion corep ',I,'/',NFSYM
!        CALL PRQMAT(FMO(1+I2TMOT(I)),
!    &    NTMO(I),NTMO(I),NTMO(I),NTMO(I),
!    &    NZ,IPQTOQ(1,0),LUPRI)

!         call shift_virtuals(fmo(I2TMOT(I)+1),nz,ntmo(i),
!    &                        istart,dlshif)

          CALL DENST1(FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,
     &              DLSHIF,D1,
     &              CMO(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &              ISTART,NVEC,NTMO(I))
        ENDIF
      ENDDO
C
      CALL QEXIT('LEVSHI')
      RETURN
      END

      subroutine shift_virtuals(fock,nz,lrc,offset,shift)
      implicit none
      !input/output  
      real(8), intent(inout) :: fock(lrc,lrc,nz)
      real(8), intent(in)    :: shift
      integer, intent(in)    :: nz
      integer, intent(in)    :: lrc
      integer, intent(in)    :: offset

      ! local variables
      integer                :: i, j
      
      print *, 'total dimension: lrc + offset',lrc, offset
      do i = 1, nz
        do j = offset, lrc
          print *, 'fock before lshift: j,i ',j,i,fock(j,j,i)
          fock(j,j,i) = fock(j,j,i) + shift
          print *, 'fock after lshift:  j,i ',j,i,fock(j,j,i)
        end do
      end do

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE UNDOLEVSHI(EIG)
C***********************************************************************
C
C     undo Level shift
C
C     Written by S. Hoefener and S. Knecht   Nov 15 2011
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D1 = 1.0D0,D0=0.0D0)
#include "dcbdhf.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION EIG(*)
!#def LEVELSHIFT_DEBUG

      CALL QENTER('UNDOLEVSHI')
C
      ISTART = 0
      DO I = 1,NFSYM
        ISTART = IORB(I) + NPSH(I)+NISHMF(I)+1
#ifdef LEVELSHIFT_DEBUG
        write(lupri,*) ' inactive offset is ',istart
#endif
        DO IOPEN = 1,NOPEN
          IF(OLEV(IOPEN).NE.D0) THEN
            DO J = 1, NACSHMF(I,IOPEN)
              EIG(ISTART+J-1) = EIG(ISTART+J-1) - OLEV(IOPEN)
            END DO
          ENDIF
          ISTART = ISTART + NACSHMF(I,IOPEN)
        END DO
        NVEC = NESH(I) - NOCCMF(I)
#ifdef LEVELSHIFT_DEBUG
        write(lupri,*) ' blabla nvec is ',nvec
        write(lupri,*) ' istart is ',istart
        write(lupri,*) ' undo... i, nesh(i), noccmf(i) ',
     &                           i, nesh(i), noccmf(i)
#endif
        IF(DLSHIF.NE.D0) THEN
          DO J = 1, NVEC
            EIG(ISTART+J-1) = EIG(ISTART+J-1) - DLSHIF
#ifdef LEVELSHIFT_DEBUG
            write(lupri,*) 'orb, old new epsilon ==> ',ISTART+J-1,
     &                      EIG(ISTART+J-1)+DLSHIF,
     &                      EIG(ISTART+J-1)
#undef LEVELSHIFT_DEBUG
#endif
          END DO
        ENDIF
      END DO
C
      CALL QEXIT('UNDOLEVSHI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck precmos */
      SUBROUTINE PRECMOS(COFT,COFU,TMAT,DOSEL,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate untransformed coefficients (MO basis)
C     and prepare for overlap selection
C
C     INPUT:
C       COFT - MO coefficients in SO-basis; they have been read by READCMO
C              and so the number of coefficients is NORB(IFRP) in symmetry IFRP.
C      
C     Written by Trond Saue, July 10 2006
C
C**********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
#include "dcbpsi.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      LOGICAL DOSEL,TOBE
      CHARACTER TEXT*74
      DIMENSION COFT(*),COFU(*),TMAT(*),WORK(*)
      DIMENSION IDIM(3,2)
C.....read SAOMO at rec#2 of LUTMAT
      CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','PRECMOS')
      READ(LUTMAT)
      CALL READT(LUTMAT,N2TMT,TMAT)
      CLOSE(LUTMAT,STATUS='KEEP')
C
C     Transform MO-coefficients to orthonormal basis
C     ==============================================
C
C.....transform coefficients to orthonormal basis --> COFU
C     Frozen orbitals must be projected out, see ADJCMO
      DO IFRP = 1,NFSYM
         IF(NORB(IFRP).EQ.0) CYCLE
         IF(NFRO(IFRP).EQ.0) THEN
          CALL QGEMM(NTMO(IFRP),NTMO(IFRP),NFBAS(IFRP,0),D1,
     &         'H','N',IPQTOQ(1,0),TMAT(I2TMT(IFRP)+1),
     &                 NFBAS(IFRP,0),NTMO(IFRP),NZT,
     &         'N','N',IPQTOQ(1,0),COFT(ICMOQ(IFRP)+1),
     &                 NFBAS(IFRP,0),NORB(IFRP),NZ,
     &              D0,IPQTOQ(1,0),COFU(I2TMOT(IFRP)+1),
     &                 NTMO(IFRP),NTMO(IFRP),NZ)
        ELSE
          CALL ADJCMO(COFT(ICMOQ(IFRP)+1),COFU(I2TMOT(IFRP)+1),
     &                TMAT(I2TMT(IFRP)+1),
     &                NFBAS(IFRP,0),NTMO(IFRP),NORB(IFRP),IFRP,
     &                 WORK,KFREE,LFREE)
        ENDIF
        IF(IPRINT.GE.10) THEN
          WRITE(LUPRI,'(A,I2)') 
     &     'PRECMOS: Coefficients in SO-basis, fermion ircop ',IFRP
          CALL PRQMAT(COFT(ICMOQ(IFRP)+1),NFBAS(IFRP,0),NORB(IFRP),
     &                NFBAS(IFRP,0),NORB(IFRP),NZ,IPQTOQ(1,0),LUPRI)
          WRITE(LUPRI,'(A,I2)') 
     &     'PRECMOS: Coefficients in orthonormal basis, fermion ircop '
     &     ,IFRP
          CALL PRQMAT(COFU(I2TMOT(IFRP)+1),NTMO(IFRP),NTMO(IFRP),
     &                NTMO(IFRP),NTMO(IFRP),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
      ENDDO
C.....write coefficients in orthonormal basis to file
      REWIND LUCMOS
      CALL WRITT(LUCMOS,N2TMOTQ,COFU)
C
C     Prepare for overlap selection
C     =============================
C
      IF(OVLSEL) THEN
        DOSEL = .TRUE.
        INQUIRE(FILE='DFSMOS',EXIST=TOBE)
        IF(.NOT.TOBE .OR. L1ORBM) CALL SMOGEN(COFU,TMAT)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Smogen */
      SUBROUTINE SMOGEN(CMO,SCMO)
C***********************************************************************
C
C     Dump (occupied) coefficients in orthonormal basis to the file DFSMOS
C     to be used for overlap selection.
C
C     Written by T. Saue July 10 2006
C     Last changes: M.Ilias, Aug.2008 - update for all electronic shells
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"

      DIMENSION CMO(*),SCMO(*)
      NSMOTQ = 0
      DO IFRP = 1,NFSYM
        ISMOQ(IFRP) = NSMOTQ
        ISOFF = ISMOQ(IFRP) + 1
        ICOFF = I2TMOT(IFRP) + NPSH(IFRP)*NTMO(IFRP) + 1
        NSDIM = NESHMF(IFRP)*NTMO(IFRP)
        NCDIM = NTMO(IFRP)*NTMO(IFRP)
        DO IZ = 1,NZ
          CALL DCOPY(NSDIM,CMO(ICOFF),1,SCMO(ISOFF),1)
          ISOFF = ISOFF + NSDIM
          ICOFF = ICOFF + NCDIM
        ENDDO
!       ... NSMOTQ to be stored into "dcbdhf.h"
        NSMOTQ = NSMOTQ + NESHMF(IFRP)*NTMO(IFRP)*NZ
      ENDDO
C
C     Write to file
C
      CALL OPNFIL(LUSMOS,'DFSMOS','UNKNOWN','SMOGEN')
      REWIND (LUSMOS)
      CALL WRITT(LUSMOS,NSMOTQ,SCMO)
      CLOSE(LUSMOS,STATUS='KEEP')
C
      IF (IPRSCF.GE.5) THEN
        WRITE(LUPRI,'(2X,A)') 'SMOGEN: cofficients dumped into DFSMOS'
      ENDIF

      RETURN
      END
#ifdef MOD_UNRELEASED
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck onevac */
      SUBROUTINE ONEVAC(CMO,EIG,IBEIG,FOCK,WORK,LWORK)
C***********************************************************************
C
C     For one-electron system add vacuum polarization potential
C     and rediagonalize
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION CMO(*),EIG(*),IBEIG(*),FOCK(*),WORK(*),NB(2),NO(2)
C
#include "memint.h"
C
C     Some memory allocation
C
      CALL MEMGET2('REAL','DMAT',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)

      call SaveTaskDistribFlags(saveflags)
      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. , .TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      if (NPOS.GT.0) THEN
         CALL MEMGET2('INTE','NPOS',KPOS,NPOS,WORK,KFREE,LFREE)
      else
         KPOS = KFREE
      endif
C
C     Modify some variables to make things work...
C
      ONESYS = .FALSE.
      NFMAT  = 1
      NOPEN  = 0
      DO IFRP = 1,NFSYM
        NSSH(IFRP) = NSSH(IFRP)+NOCC(IFRP)
        NB(IFRP)   = NISH(IFRP)
        NISH(IFRP) = 0
        NO(IFRP)   = NOCC(IFRP)
        NOCC(IFRP) = 0
      ENDDO
      DOCCNV = .FALSE.
C
C     Solve one-electron system with vacuum polarisation added.
C
      OPEN_FAC_save = OPEN_FAC
      if (OPEN_FAC .lt. 0.0d0) OPEN_FAC = 1.0D0
      CALL PREDHF(CMO,EIG,IBEIG,FOCK,WORK(KDMAT),DUM,WORK(KPOS),
     &            WORK(KFREE),LFREE)
      OPEN_FAC = OPEN_FAC_save
C
C     Reset variables
C
      ONESYS = .TRUE.
      DO IFRP = 1,NFSYM
        NISH(IFRP) = NB(IFRP)
        NOCC(IFRP) = NO(IFRP)
        NSSH(IFRP) = NSSH(IFRP)-NOCC(IFRP)
      ENDDO
C
C     Release memory
      CALL MEMREL('ONEVAC',WORK,1,KWORK,KFREE,LFREE)
C
      RETURN
      END
#endif
      SUBROUTINE FMO_FREEZE(FMO)
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION FMO(NORBT,NORBT,NZ)
      DO IFR=1,NFSYM
         DO I=1,NMOFREEZE(1)
            IF (IFR.EQ.1) THEN
               IIORB = IMOFREEZE(I) + NPSH(1)
            ELSE
               IIORB = IMOFREEZE(I) + NORB(1) + NPSH(2)
            ENDIF
            DO IZ=1,NZ
               CTMP = FMO(IIORB,IIORB,IZ)
               DO J=1,NORBT
                  FMO(J,IIORB,IZ) = 0.0D0
                  FMO(IIORB,J,IZ) = 0.0D0
               ENDDO
               FMO(IIORB,IIORB,IZ) = CTMP
            ENDDO               
         ENDDO
      ENDDO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Zoraeq */
      SUBROUTINE ZORAEQ(IFRP,ISUB,NZS,NL,NS,NT,EIG,IEIGP,IEIGE,
     &                  FMO,FAO,WORK,KFREE,LFREE)
C***********************************************************************
C         
C         Set up and solve the ZORA equation.
C         A cut and paste job by T. Saue Dec 5 2006
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbham.h"
      DIMENSION EIG(*),FAO(*),FMO(*),WORK(*)
C
      KFRSAV = KFREE
      NBRP = 4 / NZ
C
C     Allocate memory for the ZORA diagonalization. We need to
C     have the LL, LS, SL and SS blocks in separate blocks.
C
      NLL = NL * NL
      NLS = NL * NS
      NSS = NS * NS
      NTT = NT * NT
      CALL MEMGET2('REAL','FLL',KFLL,NZ*NLL,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','FLS',KFLS,NZ*NLS,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','FSS',KFSS,NZ*NSS,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','VLL',KVLL,NZ*NLL,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','VSL',KVSL,NZ*NLS,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','VSS',KVSS,NZ*NSS,WORK,KFREE,LFREE)
C
C     Extract LL block of Fock matrix
      CALL EXTRSB (1,IFRP,NBORB,4,NBRP,NBORB,4,NBRP,
     &             ISUB,ISUB,1,1,1,NZS,
     &             WORK(KFLL),NL,NL,NZS,
     &             FMO(I2ORBT(IFRP)+1),NT,NT,NZ)
C     Extract LS block of Fock matrix
      CALL EXTRSB (1,IFRP,NBORB,4,NBRP,NBORB,4,NBRP,
     &             ISUB,ISUB,1,2,1,NZS,
     &             WORK(KFLS),NL,NS,NZS,
     &             FMO(I2ORBT(IFRP)+1),NT,NT,NZ)
C     Extract SS block of Fock matrix
      CALL EXTRSB (1,IFRP,NBORB,4,NBRP,NBORB,4,NBRP,
     &             ISUB,ISUB,2,2,1,NZS,
     &             WORK(KFSS),NS,NS,NZS,
     &             FMO(I2ORBT(IFRP)+1),NT,NT,NZ)
C
C     Do the diagonalization using the ZORA approximation
C
      CALL ZORADG(NZ,NL,NS,WORK(KFLL),WORK(KFLS),WORK(KFSS),
     &         WORK(KVLL),WORK(KVSL),WORK(KVSS),EIG(IEIGE),
     &         EIG(IEIGP),WORK,KFREE,LFREE)
C
C     Put the MO vectors in the FAO matrix
C
      IF(SUB_BL) THEN
C
C       Insert LL block of Fock matrix
        CALL EXTRSB (2,IFRP,NBORB,4,NBRP,NBORB,4,NBRP,
     &               ISUB,ISUB,1,1,1,NZS,
     &               WORK(KVLL),NL,NL,NZS,
     &               FAO(I2ORBT(IFRP)+1),NT,NT,NZ)
C       Insert SL block of Fock matrix
        CALL EXTRSB (2,IFRP,NBORB,4,NBRP,NBORB,4,NBRP,
     &               ISUB,ISUB,2,1,1,NZS,
     &               WORK(KVSL),NS,NL,NZS,
     &               FAO(I2ORBT(IFRP)+1),NT,NT,NZ)
      ELSE
C
C         Note that we cannot use EXTRSB as this routine will sort the
C         eigenctors assuming a boson block structure which does not
C         apply here. The point is that we now have the solutions
C         on the right index (ordered as positron, electron) while
C         the basis is on the left ordered via boson irreps first.
C         In the spinfree ZORA case we do use EXTRSB because the
C         boson subblocking is now also (enforced) present for the
C         solutions.
C
        CALL EXTRSBR (2,IFRP,0,0,1,1,1,NZ,WORK(KVLL),NL,NL,NZ,
     &                FAO(I2ORBT(IFRP)+1),NT,NT)
        CALL EXTRSBR (2,IFRP,0,0,2,1,1,NZ,WORK(KVSL),NS,NL,NZ,
     &                FAO(I2ORBT(IFRP)+1),NT,NT)
      ENDIF
C
C     Done; release the extra memory used in ZORA
C
      CALL MEMREL('ZORAEQ',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Adjcmo */
      SUBROUTINE ADJCMO(COFT,COFU,TMAT,NB,NT,NO,IFRP,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Project out frozen orbitals from coefficients in orthonormal basis
C     For each orbital class the density matrix in the orthonormal basis
C     is generated, then diagonalized and the null vectors eliminated.
C
C     Written by T. Saue Feb 4 2009
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbdhf.h"
      DIMENSION COFT(NB,NO,NZ),COFU(NT,NT,NZ),TMAT(NB,NT,NZT),WORK(*)
      NF = 0
      IOFT = 1
      IOFU = 1
C.....negative energy shells
      IF(NPSH(IFRP).GT.0) THEN
        CALL ADJCMO1(COFT(1,IOFT,1),COFU(1,IOFU,1),TMAT,
     &               NPSH(IFRP),NPSHMF(IFRP),NB,NO,NT,
     &               WORK,KFREE,LFREE)
        IOFT = IOFT + NPSH(IFRP)
        IOFU = IOFU + NPSHMF(IFRP)
        NF   = NF + NPSH(IFRP) - NPSHMF(IFRP)
      ENDIF
C.....inactive shells
      IF(NISH(IFRP).GT.0) THEN
        CALL ADJCMO1(COFT(1,IOFT,1),COFU(1,IOFU,1),TMAT,
     &               NISH(IFRP),NISHMF(IFRP),NB,NO,NT,
     &               WORK,KFREE,LFREE)
        IOFT = IOFT + NISH(IFRP)
        IOFU = IOFU + NISHMF(IFRP)
        NF   = NF + NISH(IFRP) - NISHMF(IFRP)
      ENDIF
C.....active shells
      DO IOPEN = 1,NOPEN
        IF(NACSH(IFRP,IOPEN).GT.0) THEN
          CALL ADJCMO1(COFT(1,IOFT,1),COFU(1,IOFU,1),TMAT,
     &                 NACSH(IFRP,IOPEN),NACSHMF(IFRP,IOPEN),NB,NO,NT,
     &                 WORK,KFREE,LFREE)
          IOFT = IOFT + NACSH(IFRP,IOPEN)
          IOFU = IOFU + NACSHMF(IFRP,IOPEN)
          NF   = NF + NACSH(IFRP,IOPEN) - NACSHMF(IFRP,IOPEN)
        ENDIF
      ENDDO
C.....virtuals
      IF(NSSH(IFRP).GT.0) THEN
        NSSHMF = NSSH(IFRP) - NFRO(IFRP) + NF
        CALL ADJCMO1(COFT(1,IOFT,1),COFU(1,IOFU,1),TMAT,
     &               NSSH(IFRP),NSSHMF,NB,NO,NT,
     &               WORK,KFREE,LFREE)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck adjcmo1 */
      SUBROUTINE ADJCMO1(COFT,COFU,TMAT,NVEC,NVECMF,NB,NO,NT,
     &                   WORK,KFREE,LFREE)
C***********************************************************************
C
C     Project frozen orbitals out of this shell 
C     (inactive/active or secondary).
C     Coefficients are in orthonormal basis and should therefore
C     be orthonormal, that is
C       S_{ij} = C(dagger)_ik * C_kj} = delta_{ij}
C     The metric is diagonalized and the coefficients rotated to
C     eliminate the null space.
C
C     Written by Trond Saue Feb 6 2009
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
#include "dgroup.h"
#include "dcbbas.h"
      DIMENSION COFT(NB,NO,NZ),COFU(NT,NT,NZ),TMAT(NB,NT,NZT),WORK(*)
C
C     
      KFRSAV = KFREE
      IF(NVECMF.EQ.NVEC) THEN
        CALL QGEMM(NT,NVEC,NB,D1,
     &       'H','N',IPQTOQ(1,0),TMAT,NB,NT,NZT,
     &       'N','N',IPQTOQ(1,0),COFT,NB,NO,NZ,
     &            D0,IPQTOQ(1,0),COFU,NT,NT,NZ)
      ELSEIF(NVECMF.GT.0) THEN
        NDDIM = NVEC*NVEC*NZ
        NCDIM = NT*NVEC*NZ
        CALL MEMGET2('REAL','D',KD,NDDIM,WORK,KFREE,LFREE)      
        CALL MEMGET2('REAL','C',KC,NCDIM,WORK,KFREE,LFREE)      
        CALL QGEMM(NT,NVEC,NB,D1,
     &       'H','N',IPQTOQ(1,0),TMAT,NB,NT,NZT,
     &       'N','N',IPQTOQ(1,0),COFT,NB,NO,NZ,
     &            D0,IPQTOQ(1,0),WORK(KC),NT,NVEC,NZ)
C.......Form metric
        CALL QGEMM(NVEC,NVEC,NT,D1,
     &       'H','N',IPQTOQ(1,0),WORK(KC),NT,NVEC,NZ,
     &       'N','N',IPQTOQ(1,0),WORK(KC),NT,NVEC,NZ,
     &            D0,IPQTOQ(1,0),WORK(KD),NVEC,NVEC,NZ)
        CALL MEMGET2('REAL','V',KV,NDDIM,WORK,KFREE,LFREE)      
        CALL MEMGET2('REAL','E',KE,NVEC,WORK,KFREE,LFREE)      
C.......Diagonalize metric; eigenvalues are in ascending order
        CALL QDIAG(NZ,NVEC,WORK(KD),NVEC,NVEC,WORK(KE),1,
     &             WORK(KV),NVEC,NVEC,WORK(KFREE),LFREE,IERR)
C.....  Rotate orbitals into place
        IOFF = NVEC*(NVEC-NVECMF)
        CALL QGEMM(NT,NVECMF,NVEC,D1,
     &       'N','N',IPQTOQ(1,0),WORK(KC),NT,NVEC,NZ,
     &       'N','N',IPQTOQ(1,0),WORK(KV+IOFF),NVEC,NVEC,NZ,
     &            D0,IPQTOQ(1,0),COFU,NT,NT,NZ)
        CALL MEMREL('ADJCMO1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ATOMIC_START(DMAT,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Routine for atomic start
C     NNEG - number of positronic solutions
C     NPOS - number of electronic solutions
C     NBAS - number of AO-basis functions
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
#include "mxcent.h"
C
#include "dcbdhf.h"
#include "dcbatom.h"
#include "nuclei.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION DMAT(NTBAS(0),NTBAS(0),NZ),WORK(*)
C.....Local variables
      LOGICAL TOBE,FNDLAB
      CHARACTER TEXT*74
      INTEGER NNEG,NPOS,NBAS,IDIM(3)
C
      KFRSAV=KFREE
C.....Say hello
      IF(IPRINT.GE.1) THEN
        CALL HEADER('Atomic start',-1)
      ENDIF
C.....Initialize buffer density matrix
      NDMAT=NTBAS(0)*NTBAS(0)*NZ
      NDBUF=NTBAS(0)*NTBAS(0)*NZC1
      CALL MEMGET2('REAL','DBUF',KDBUF,NDBUF,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KDBUF),NDBUF)
C.....Loop over atomic types
      LUATOM=137
      CALL MEMGET2('INTE','INFO',KINF,NTBAS(0),WORK,KFREE,LFREE)
      ICENT = 1
      DO IATOM = 1,NONTYP
C.......Find the total number of basis functions for this atom type
        CALL LABCOUNT(MBAS,WORK(KINF),NTBAS(0),1,-1,ICENT,1,-1)
        IF (MBAS == 0) CYCLE
C.......Read atomic coefficients in C1 format
        INQUIRE(FILE=ATOMFIL(IATOM),EXIST=TOBE)         
        IF(.NOT.TOBE) THEN
           WRITE(LUPRI,'(A,A6,A)') 
     &     'ATOMIC_START: Atomic coefficient file ',ATOMFIL(IATOM),
     &     ' not found !'
           CALL QUIT
     &     ('ATOMIC_START: Atomic coefficient file not found !')
        ENDIF
        CALL OPNFIL(LUATOM,ATOMFIL(IATOM),'OLD','ATOMIC_START')
        IF(FNDLAB('INFO    ',LUATOM)) THEN
           READ(LUATOM) TEXT,NSYM,NZBUF,NNEG,NPOS,NBAS,TOTERG
        ELSE
           WRITE(LUPRI,*) 'ATOMIC_START: error reading info ',
     &                    ATOMFIL(IATOM) 
          CALL QUIT('ATOMIC START: Error reading info')
        ENDIF
        IF(NBAS.NE.MBAS) THEN
           WRITE(LUPRI,'(A,A,I3/A,I6/A,A6,A,I6)') 
     &     'ATOMIC_START: ',
     &     'Wrong number of basis function for atomic type ',IATOM,
     &     'Correct number          : ', MBAS,
     &     'Number read from ',ATOMFIL(IATOM),' : '  ,NBAS
           CALL QUIT
     &     ('ATOMIC_START: Wrong number of basis functions !')
        ENDIF
        IF(IPRINT.GE.1) THEN
          WRITE(LUPRI,*) NAMN(ICENT),ATOMFIL(IATOM)
          WRITE(LUPRI,*) TEXT,NSYM,NNEG,NPOS,NBAS,TOTERG
        ENDIF
        NAORB   = NNEG+NPOS
        NACMO   = NBAS*NAORB
        NACMOQ  = NACMO*NZC1
        NADEN   = NBAS*NBAS*NZC1
        IDIM(1) = NNEG
        IDIM(2) = NPOS
        IDIM(3) = NBAS
C.......Allocate memory for atomic density matrix and coefficients
        CALL MEMGET2('REAL','ADEN',KADE,NADEN ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','CMO ',KCMO,NACMOQ,WORK,KFREE,LFREE)
        CALL REACMO_coefficients(LUATOM,WORK(KCMO),NACMOQ,1,NZC1,
     &           NAORB,NBAS,NNEG,NPOS,IDIM)         
C.......Loop over closed/open shells and construct contribution to atomic density matrix
        CALL MEMGET2('INTE','VEC ',KVEC,NAORB,WORK,KFREE,LFREE)
        NDIM=MAX(NAORB+1,NBAS)
        CALL MEMGET2('INTE','IBUF',KIBUF,NDIM,WORK,KFREE,LFREE)
        FAC = D0
        DO ISHELL = 1,NVECATOM(IATOM)
C.........starting from orbital string VECATOM return index array WORK(KVEC) of selected orbitals
          NSEL=1
          CALL NUMLS1(VECATOM(ISHELL,IATOM),WORK(KVEC),
     &        NAORB,-NNEG,NPOS,NSEL,WORK(KIBUF))
C.........count number NPVEC and NEVEC of negative- and positive-energy orbitals, respectively
          CALL ORBCNT(WORK(KVEC),NSEL,NNEG,NPOS,NPVEC,NEVEC)
C.........extract selected vectors from CMO and place in CBUF
          CALL MEMGET2('REAL','CBUF',KCBUF,NBAS*NSEL*NZC1,
     &                 WORK,KFREE,LFREE)
          CALL INDEXFILL(WORK(KIBUF),NBAS,1,1)
          CALL SELOWC(NZC1,WORK(KCMO),NBAS,NAORB,NNEG,NPOS,
     &                 WORK(KCBUF),NBAS,NSEL,WORK(KVEC),WORK(KIBUF),
     &                 NPVEC,NEVEC)
          IF(IPRINT.GE.4) THEN
            WRITE(LUPRI,'(A,I3)') '* Atomic coefficients ',IATOM
            CALL PRQMAT(WORK(KCMO),NBAS,NAORB,NBAS,NAORB,
     &                NZ,IPQTOQ(1,0),LUPRI)
            WRITE(LUPRI,'(A,I3)') '* Selected Atomic coefficients',IATOM
            CALL PRQMAT(WORK(KCBUF),NBAS,NSEL,NBAS,NSEL,
     &                NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
C.........Generate contribution to atomic density matrix from this shell
          CALL DENST1(WORK(KADE),NBAS,NBAS,NZC1,OCCATOM(ISHELL,IATOM),
     &                FAC,WORK(KCBUF),NBAS,NSEL,1,NSEL,NBAS)
          IF(IPRINT.GE.4) THEN
            WRITE(LUPRI,'(A,I3)') '* Density matrix from atom ',IATOM
            CALL PRQMAT(WORK(KADE),NBAS,NBAS,NBAS,NBAS,
     &                 NZC1,IPQTOQ(1,0),LUPRI)
          ENDIF
          FAC = D1
          CALL MEMREL('ATOMIC_START.shell',WORK,1,KCBUF,KFREE,LFREE)
        ENDDO
C.......Now we have the atomic density matrix. 
C       Put all contributions from this atomic type into the 
C       molecular C1 density matrix
        CALL MEMGET2('INTE','R1BUF',KR1BUF,NTBAS(0),WORK,KFREE,LFREE)      
        CALL MEMGET2('INTE','C1BUF',KC1BUF,NTBAS(0),WORK,KFREE,LFREE)      
        CALL MEMGET2('INTE','I1BUF',KI1BUF,NBAS,WORK,KFREE,LFREE)      
        CALL MEMGET2('INTE','J1BUF',KJ1BUF,NBAS,WORK,KFREE,LFREE) 
        CALL PUT_AOBLOCK(WORK(KDBUF),NTBAS(0),NTBAS(0),ICENT,
     &                   NONT(IATOM),NUCDEG(ICENT),NZC1,
     &                   WORK(KADE),NBAS,NBAS,-1,NBAS,NBAS,-1,
     &                   WORK(KR1BUF),WORK(KC1BUF),
     &                   WORK(KI1BUF),WORK(KJ1BUF))
C.......Release memory        
        CALL MEMREL('ATOMIC_START',WORK,1,KCMO,KFREE,LFREE)
        CLOSE(LUATOM,STATUS='KEEP')
        ICENT = ICENT+NONT(IATOM)
      ENDDO
C
      IF(NBSYM.EQ.1) THEN
        CALL DCOPY(NDMAT,WORK(KDBUF),1,DMAT,1)
      ELSE
C.......transform atomic start to actual symmetry
        IF(IPRINT.GE.5) THEN 
          WRITE(LUPRI,'(A)') '* C1 Density matrix from atomic guess'
          CALL PRQMAT(WORK(KDBUF),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NZC1,IQDEF,LUPRI)
        ENDIF
        CALL DZERO(DMAT,NDMAT)
        IREP = 0
        IOFF = 0
        NOFF = NTBAS(0)*NTBAS(0)
        DO IZ = 1,NZC1
           IREPD = IRQMAT(IZ,IREP)
           IQ    = IQMULT(1,JQBAS(IREPD,1),IZ)
           IPQ   = IQTOPQ(IQ,IREP)
           CALL DTAOSO(WORK(KDBUF+IOFF),DMAT(1,1,IPQ),NTBAS(0),
     &                 IREPD,IPRINT)
           IOFF = IOFF + NOFF
        ENDDO
        CALL BUTOBS(DMAT,NZ,WORK(KDBUF),NDBUF)
        DO IZ = 1,NZ
           IQ = IPQTOQ(IZ,IREP)
           CALL Q2BPHASE('D',IQ,1,DMAT(1,1,IZ))
        ENDDO
      ENDIF
      IF(IPRINT.GE.4) THEN
        WRITE(LUPRI,'(A)') '* Density matrix from atomic guess'
        DO IFRP = 1,NFSYM
          WRITE(LUPRI,'(/A,A3)')
     &      '* Fermion ircop ',FREP(IFRP)
          CALL PRQMAT(DMAT(I2BASX(IFRP,IFRP)+1,1,1),
     &              NFBAS(IFRP,0),NFBAS(IFRP,0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDDO
      ENDIF
      CALL MEMREL('ATOMIC_START',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADHOC(CMO,EIG,IBEIG,EIGSUM,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Atomic Huckel start
C     The basic theory is given here
C       R. Hoffmann, J. Chem. Phys. 39 (1963) 1397
C     This work was further motivated by 
C       P. Norman and H. J. Aa. Jensen, Chem. Phys. Lett. 531 (2012) 229
C     The particularity of the present approach is to employ 
C     pre-calculated atomic fragments
C
C     Written by Trond Saue June 5 2012
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
#include "mxcent.h"
C
#include "nuclei.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbhoc.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      LOGICAL FNDLAB
      CHARACTER TEXT*74
      DIMENSION IDIM(3)
      DIMENSION NCATOM(0:2,MAXHOC),KCATOM(MAXHOC),NAO(MAXHOC),
     &          NATDEG(MAXHOC),CMO(*),EIG(*),IBEIG(*),WORK(*)
      KFRSAV = KFREE
      CALL MEMGET2('INTE','INFO',KINF,NTBAS(0),WORK,KFREE,LFREE)
      IF(IPRINT.GE.2) THEN
        WRITE(6,'(A)') '** Output from ADHOC (atomic Huckel start)'
      ENDIF
C
C     Find total number of reference orbitals
C     =======================================
C     ININD runs over symmetry independent centers
C     INDEP runs over symmetry dependent centers
      NVECS=0
      INDEP=1
      ININD=1
      LUHOC=137
      DO ITYP = 1,NONTYP
C.......Find the total number of basis functions for this atom type
        CALL LABCOUNT(MBAS,WORK(KINF),NTBAS(0),1,-1,ININD,1,-1)
C.......Read information from coefficient file
        CALL OPNFIL(LUHOC,HOCFIL(ITYP),'OLD','ADHOC')
        IF(FNDLAB('INFO    ',LUHOC)) THEN
          READ (LUHOC,END=10,ERR=20) TEXT,NSYM,NZBUF,
     &          NCATOM(2,ITYP),NCATOM(1,ITYP),NAO(ITYP),DUM
        ELSE
          WRITE(6,*) 'ADHOC: not able to read info from file ',
     &                HOCFIL(ITYP)
          CALL QUIT('ADHOC: Error reading info')
        ENDIF
C.............NPSH............NESH............NBAS...........TOTERG
        CLOSE(LUHOC,STATUS='KEEP')
        IF(NAO(ITYP).NE.MBAS) THEN
           WRITE(LUPRI,'(A/A,I6/A,A6,A,I6)') 
     &     '* Atomic Huckel start: ',
     &     '  Wrong number of basis function for atomic type ', ITYP,
     &     '  Number read from                     ',
     &     HOCFIL(ITYP),' : '  ,NAO(ITYP),
     &     '  Expected                             ',MBAS
           CALL QUIT
     &     ('ADHOC: Wrong number of basis functions !')
        ENDIF
C.......Count the number of orbitals to read
        NCATOM(0,ITYP)=0
        NDIM=NCATOM(2,ITYP)+NCATOM(1,ITYP)+1
        CALL MEMGET2('INTE','ITMP ',KITMP,NDIM,WORK,KFREE,LFREE)
        CALL NUMLS1(VECHOC(ITYP),IDUM,IDUM,
     &       -NCATOM(2,ITYP),NCATOM(1,ITYP),
     &       NCATOM(0,ITYP),WORK(KITMP))
        CALL MEMREL('ADHOC.ref',WORK,1,KITMP,KFREE,LFREE)
        CALL MEMGET2('INTE','CATOM',KCATOM(ITYP),NCATOM(0,ITYP),
     &              WORK,KFREE,LFREE)
        NN = INDEP
        DO J = 0,NONT(ITYP)-1
          INDEP = INDEP + NUCDEG(ININD+J)
        ENDDO
        NATDEG(ITYP) = INDEP - NN
        NVECS = NVECS + NATDEG(ITYP)*NCATOM(0,ITYP)
        ININD = ININD+NONT(ITYP)
      ENDDO
C
C     Calculate dimension of arrays
C     =============================
C
      CALL MEMGET2('REAL','E',KE,NVECS,WORK,KFREE,LFREE)
      NQ = NTBAS(0)*NVECS*NZC1
      CALL MEMGET2('REAL','Q',KQ,NQ      ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KQ),NQ)
C
C     Select atomic coefficients
C     ==========================
C
      JVEC = 1
      ININD= 1
      DO ITYP = 1,NONTYP
C.......Read C1 coefficients and orbitalenergies
        CALL OPNFIL(LUHOC,HOCFIL(ITYP),'OLD','PRJAOB')        
        IF(FNDLAB('INFO    ',LUHOC)) THEN
          READ (LUHOC,END=10,ERR=20) TEXT,NSYM,NZBUF,
     &          NNEG,NPOS,NAO(ITYP),DUM
        ELSE
          WRITE(6,*) 'ADHOC2: not able to read info from file ',
     &                HOCFIL(ITYP)
          CALL QUIT('ADHOC2: Error reading info')
        ENDIF
        NTOT = NNEG+NPOS
        IDIM(1) = NNEG
        IDIM(2) = NPOS
        IDIM(3) = NAO(ITYP)
        NQA  = NAO(ITYP)*NTOT*NZC1
        CALL MEMGET2('REAL','QA',KQA,NQA ,WORK,KFREE,LFREE)
        CALL REACMO_coefficients(LUHOC,WORK(KQA),NQA,1,4,
     &           NTOT,NAO(ITYP),NNEG,NPOS,IDIM)         
        CALL MEMGET2('REAL','EA',KEA,NTOT,WORK,KFREE,LFREE)
        CALL REACMO_eigenvalues(LUHOC,WORK(KEA),NTOT,1,4,
     &           NTOT,NNEG,NPOS,IDIM)         
        CLOSE(LUHOC,STATUS='KEEP')
C.......Make pointer array to selected reference orbitals
        NDIM = NTOT+1
        CALL MEMGET2('INTE','ITMP',KITMP,NDIM,WORK,KFREE,LFREE)
        CALL NUMLS1(VECHOC(ITYP),WORK(KCATOM(ITYP)),NCATOM(0,ITYP),
     &       -NNEG,NPOS,NCATOM(0,ITYP),WORK(KITMP))
        CALL MEMREL('ADHOC.refsel',WORK,1,KITMP,KFREE,LFREE)
        CALL ORBCNT(WORK(KCATOM(ITYP)),NCATOM(0,ITYP),NNEG,NPOS,
     &             NCATOM(2,ITYP),NCATOM(1,ITYP))
C.....Extract coefficients and orbital energies (JVEC updated inside)
        CALL SELATOM(WORK(KQ),WORK(KE),NTBAS(0),NVECS,NZC1,
     &               WORK(KQA),WORK(KEA),NAO(ITYP),NTOT,NNEG,NPOS,
     &               ININD,NONT(ITYP),NUCDEG(ININD),
     &               WORK(KCATOM(ITYP)),NCATOM(1,ITYP),
     &               NCATOM(2,ITYP),JVEC,WORK,KFREE,LFREE)
        CALL MEMREL('ADHOC.atom',WORK,1,KQA,KFREE,LFREE)
        ININD = ININD+NONT(ITYP)
      ENDDO
      IF(IPRINT.GE.6) THEN
        WRITE(6,*) 'ADHOC: Selected coefficients:'
        CALL PRQMAT(WORK(KQ),NTBAS(0),NVECS,NTBAS(0),NVECS,NZC1,
     &                IQDEF,LUPRI)
        WRITE(6,*) 'ADHOC: Selected eigenvalues:'
        CALL OUTPUT(WORK(KE),1,NVECS,1,1,NVECS,1,-1,LUPRI)
      ENDIF
C
C     Get overlap matrix in fragment basis
C     ====================================
C
      CALL MEMGET2('REAL','AMAT',KAMAT,NVECS*NVECS*NZC1,
     &             WORK,KFREE,LFREE)
      CALL OVLFRAG(WORK(KAMAT),NVECS,NONTYP,HOCFIL,NCATOM,NATDEG,
     &             WORK(KQ),NTBAS(0),IPRINT,WORK,KFREE,LFREE)
C
C     Do extended Huckel
C
      CALL EHUCKEL(WORK(KAMAT),WORK(KE),HUCPAR,NVECS,NEFF,NZC1,IQDEF,
     &             NONTYP,HOCFIL,NCATOM,NATDEG,IPRINT,WORK,KFREE,LFREE)
C     Construct extended Huckel molecular orbitals (C1 for now)
      IFRP=1
CTROND
      CALL DZERO(CMO,NCMOQ(IFRP))
      IOFF=NPSH(IFRP)*NFBAS(IFRP,0)+1
      CALL QGEMM(NTBAS(0),NEFF,NVECS,D1,
     &       'N','N',IQDEF,WORK(KQ),NTBAS(0),NVECS,NZC1,
     &       'N','N',IQDEF,WORK(KAMAT),NVECS,NEFF,NZC1,
     &            D0,IQDEF,CMO(IOFF),NFBAS(IFRP,0),NORB(IFRP),NZC1)
      CALL DCOPY(NEFF,WORK(KE),1,EIG(NPSH(IFRP)+1),1)
C     Write coefficients to file
      EIGSUM=D2*DSUM(NVECS,WORK(KE),1)
 137  CONTINUE
      CALL MEMREL('ADHOC',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
 10   CONTINUE
      CALL QUIT('ADHOC: END OF FILE reading TEXT')
 20   CONTINUE
      CALL QUIT('ADHOC: ERROR reading TEXT')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ehuckel */
      SUBROUTINE EHUCKEL(HMAT,EIG,HUCPAR,NVECS,NEFF,NZ,IPQ,
     &           NONTYP,HOCFIL,NCATOM,NATDEG,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Extended Huckel module
C     INPUT:
C       HMAT - overlap matrix in fragment basis
C     OUTPUT:
C       HMAT - eigenvectors of extended Huckel problem
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DP5=0.5D0)
      CHARACTER*6 HOCFIL(NONTYP)
      DIMENSION HMAT(NVECS,NVECS,NZ),EIG(NVECS),IPQ(NZ),
     &          NCATOM(0:2,NONTYP),NATDEG(NONTYP),WORK(*)
      KFRSAV=KFREE
      N2DIM=NVECS*NVECS*NZ
      CALL MEMGET2('REAL','SMAT',KSMAT,N2DIM,WORK,KFREE,LFREE)
      CALL DCOPY(N2DIM,HMAT,1,WORK(KSMAT),1)
C.....Construct extended Huckel Hamiltonian; hermiticity not exploited
      JAT = 0
      JOF = 0
      DO J = 1,NONTYP
      DO JJ = 1,NATDEG(J)
        JAT = JAT + 1
        IAT = 0
        IOF = 0
        DO I = 1,NONTYP
        DO II = 1,NATDEG(I)
          IAT = IAT + 1
          IF(IAT.EQ.JAT) THEN
            DO K = 1,NCATOM(0,I)
              HMAT(IOF+K,IOF+K,1) = EIG(IOF+K)
            ENDDO
          ELSE
            DO L = 1,NCATOM(0,J)
              DO K = 1,NCATOM(0,I)
                FAC = DP5*HUCPAR*(EIG(IOF+K)+EIG(JOF+L))
                DO IZ = 1,NZ
                  HMAT(IOF+K,JOF+L,IZ) = FAC*HMAT(IOF+K,JOF+L,IZ)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IOF = IOF + NCATOM(0,I)
        ENDDO
        ENDDO
        JOF = JOF + NCATOM(0,J)
      ENDDO
      ENDDO
C
      IF(IPRINT.GE.2) THEN
        WRITE(6,*) 'EHUCKEL: Total extended Huckel matrix'
        CALL PRQMAT(HMAT,NVECS,NVECS,NVECS,NVECS,NZ,IPQ,LUPRI)
      ENDIF
C     Solve general eigenvalue problem; coefficients stored in SMAT
      CALL GEIGSOLV(HMAT,WORK(KSMAT),EIG,NVECS,NEFF,NVECS,NVECS,NZ,IPQ,
     &              IPRINT,WORK,KFREE,LFREE)
      CALL DCOPY(N2DIM,WORK(KSMAT),1,HMAT,1)
C
      CALL MEMREL('EHUCKEL',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FMO_WRITE(FMO,IMAT,IFRP,ISUB,NBO)
C***********************************************************************
C      miro: save the matrix for diagonalization, FMO(IMAT), into formatted, unique files
C      -----------------------------------------------------------------------------------
C       only in the linear symmetry !
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"            
#include "dcbdhf.h"      
      character*100 :: hfile_format,hfile_name
      DIMENSION FMO(*)      
! get the file name first, containing symmetry
      luhm=95
      hfile_format="(A8,A7,I1,A1,I1,A7,I1,A1,I1)"
      write(hfile_name,hfile_format)
     & "hmatrix.","fermirp",IFRP,"-",NFSYM,
     & "_bosirp",ISUB,"-",N_SUB_BL(IFRP)
! open the unique file and write the matrix
      open(luhm,file=trim(hfile_name),access="sequential",
     &     form="formatted",status="unknown")
      write(luhm,*) NBO,1 ! this is pure real matrix (NZ=1)
      ic=0
      do ii=1,NBO
        do jj=1,NBO
          write(luhm,*) ii,jj,FMO(IMAT+ic)
          ic = ic + 1
        enddo
        ic = ic + NTMO(IFRP) - NBO
      enddo
      close(luhm,status="keep")
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GO2C4C()
C***********************************************************************
C***********************************************************************
         use memory_allocator
         use dirac_cfg
         use x2cmod_cfg

#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"      
#include "dcbham.h"
#include "dcbbas.h"
        LOGICAL EX
        real(8), allocatable :: WORK(:)

      CALL QENTER('GO2C4C')
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in GO2C4C')
        
      MC        = 2
      !> set integral flags
      ILLINT    = IAND(INTGEN,1)
      ISLINT    = IAND(INTGEN/2,1)
      ISSINT    = IAND(INTGEN/4,1)
      IGTINT    = IAND(INTGEN/8,1)

      print *, 'intdef, ssmtrc ==> ',intdef, ssmtrc

      BSS        = .FALSE.
      twocomp    = .false.
      x2cmod_mmf = .false.
      x2c        = .false.

      CALL SETDC1(0)
      CALL SETDC2(0)

      CALL GMOTRA(.FALSE.)

      IF (INI2C.EQ.1.OR.INI2C.EQ.2) THEN
CMI    ... prepare the starting 4c Fock MO matrix
        CALL MEMGET2('REAL','K1',K1,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','K2',K2,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','K3',K3,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMCHK('FDTR4C 0.0.0',WORK,1)

        CALL FDTR4C(WORK(K1),WORK(K2),WORK(K3),WORK(KFREE),LFREE)

        CALL MEMREL('FDTR4C.GO2C4C',WORK,KWORK,KWORK,KFREE,LFREE)

      ELSE IF (INI2C.EQ.3.OR.INI2C.EQ.4) THEN
C       .... prepare the DFCOEF file with staring MOs
        CALL MEMGET2('REAL','CMO  ',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','TMAT ',K1   ,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','EIG  ',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBEIG',KIBE,NTBAS(0),WORK,KFREE,LFREE)

        CALL CMOINI4C(WORK(KCMO),WORK(K1),WORK(KEIG),
     &                  WORK(KIBE),WORK(KFREE),LFREE)

        CALL MEMREL('FDTR4C',WORK,KWORK,KCMO,KFREE,LFREE)
      ENDIF

      INQUIRE(FILE='DFDIIS',EXIST=EX)
      IF (EX) THEN
#if defined SYS_WINDOWS
        CALL SYSTEM('del DFDIIS')
#else
        CALL SYSTEM('/bin/rm DFDIIS')
#endif
        IF (IPRSCF.GE.0) WRITE(LUPRI,*)
     &  'GO2C4C: DFDIIS file from previous SCF cycle was deleted.'
      ENDIF

      INQUIRE(FILE='DFFOCK',EXIST=EX)
      IF (EX) THEN
        CALL SYSTEM('/bin/rm DFFOCK')
        IF (IPRSCF.GE.0) WRITE(LUPRI,*)
     &  'GO2C4C: DFFOCK file from previous SCF cycle was deleted.'
      ENDIF

      INQUIRE(FILE='DFEVEC',EXIST=EX)
      IF (EX) THEN
        CALL SYSTEM('/bin/rm DFEVEC')
        IF (IPRSCF.GE.0) WRITE(LUPRI,*)
     &  'GO2C4C: DFEVEC file from previous SCF cycle was deleted.'
      ENDIF

      INQUIRE(FILE='DFCOEF',EXIST=EX)
      IF (EX.AND.(INI2C.NE.3.AND.INI2C.NE.4)) THEN
        CALL SYSTEM('/bin/rm DFCOEF')
        IF (IPRSCF.GE.0) WRITE(LUPRI,*)
     &  'GO2C4C: DFCOEF file from previous SCF cycle was deleted.'
      ENDIF

      DHFCONV(1)=.FALSE.
      DHFCONV(2)=.FALSE.
      DHFEXIT   =.FALSE.

      RESTFCK=(INI2C.EQ.1.OR.INI2C.EQ.2) 

      BNCRON = .FALSE.
      BNSPON = .FALSE.

      CALL MEMGET2('REAL','CMO  ',KCMO ,N2BBASXQ,      WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EIG  ',KEIG ,NTBAS(0),      WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBEIG',KIBE,NTBAS(0),      WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','FOCK ',KFOCK,N2BBASXQ*NFMAT,WORK,KFREE,LFREE)

      CALL SETDHF(IPRSCF)
      IF(dirac_cfg_dft_calculation) THEN
         CALL TITLER('Kohn-Shame DFT calculation','*',125)
      ELSE
         CALL TITLER('Hartree-Fock calculation','*',125)
      END IF

CMI       ... when MAXITR.eq.0, we may get pseudo4c MO from previous 2c-SCF
      IF (MAXITR.GT.0) THEN
       CALL DHFSCF(WORK(KCMO),WORK(KEIG),WORK(KIBE),WORK(KFOCK),
     &        WORK(KFREE),LFREE)
       CALL DHFOUT(WORK(KCMO),WORK(KEIG),WORK(KIBE),WORK(KFREE),
     &            LFREE,IPRSCF)
      ELSE
        DHFCONV(1)=.TRUE.
      ENDIF

      CALL MEMREL('GO2C4C.DHFOUT',WORK,KWORK,KWORK,KFREE,LFREE)

      call dealloc(work)
      CALL QEXIT('GO2C4C')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GO4C2C()
C***********************************************************************
C***********************************************************************
      use memory_allocator
      use dirac_cfg
      use x2cmod_cfg

#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"      
#include "dcbham.h"
#include "dcbbas.h"
      LOGICAL EX
      real(8), allocatable :: WORK(:)

      CALL QENTER('GO4C2C')
      if(x2cmod_mmf)then
        x2c      = .true.
      else 
        BSS      = .TRUE.
        CMPEIG   = .TRUE. ! this flag is needed for preparing transformed data
        WRITE(lupri,'(/A)')
     &    '     *** GO4C2C: Descending to two-component'// 
     &    ' mode after previous four-component SCF cycles. ***'
      end if
      TWOCOMPBSS = .FALSE. ! ensure 4comp. mode
      TWOCOMP    = .FALSE.

      START2C = .FALSE.
!
!     call the x2c/bss transformation module (branch 2 inside gmotra)
      CALL GMOTRA(.FALSE.)
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in GO4C2C')
      IF(.NOT.CONT2C) THEN
         DHFCONV(1)=.TRUE. 
         DHFCONV(2)=.TRUE. 
      ELSE                     ! continue with 2c-SCF iterations
        if(BSS)then
          CALL TITLER('Two-component BSS/DKH SCF calculation','*',125)
        end if

        BNCRON = .FALSE.
        BNSPON = .FALSE.
C       ... need to reset them...
        DHFCONV(1)=.FALSE.
        DHFCONV(2)=.FALSE.
        DHFEXIT   =.FALSE.

CMI    ...remove some files from the previous 4c-SCF cycle ...
        INQUIRE(FILE='DFDIIS',EXIST=EX)
        IF (EX) CALL SYSTEM('/bin/rm DFDIIS')
        INQUIRE(FILE='DFEVEC',EXIST=EX)
        IF (EX) CALL SYSTEM('/bin/rm DFEVEC')
        INQUIRE(FILE='DFCOEF',EXIST=EX)
        IF (EX) CALL SYSTEM('/bin/rm DFCOEF')

        CALL MEMGET2('REAL','CMO  ',KCMO ,N2BBASXQ,
     &               WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','EIG  ',KEIG ,NTBAS(0),
     &               WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBEIG',KIBE,NTBAS(0),
     &               WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','FOCK ',KFOCK,N2BBASXQ*NFMAT,
     &               WORK,KFREE,LFREE)

        RESTFCK = .TRUE.
        CALL DHFSCF(WORK(KCMO),WORK(KEIG),WORK(KIBE),WORK(KFOCK),
     &         WORK(KFREE),LFREE)

        CALL DHFOUT(WORK(KCMO),WORK(KEIG),WORK(KIBE),WORK(KFREE),
     &              LFREE,IPRSCF)

        CALL MEMREL('PSISCF.DHFOUT',WORK,KWORK,KWORK,KFREE,LFREE)
      END IF ! continue in 2c-mode switch
      call dealloc(work)
      CALL QEXIT('GO4C2C')
       RETURN
       END
C     --- end of dirac/dirscf.F ---
