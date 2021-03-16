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

C FILE    : krmcopt.F

! define task symbols for CALL DIRAC_PARCTL( task )
#include "dirac_partask.h"
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiopt */
      SUBROUTINE PSIOPT(WFTYP,WFCONV)
C***********************************************************************
C
C     Driver routine to PSIOP1
C     Get second order optimized wave function.
C     Either DHF or KR-MCSCF... (WFTYP = 'SCF   ' or 'MCSCF ')
C
C     Written by J. Thyssen and H. J. Aa. Jensen - Oct 22 1998
C     Last revision :
C
C***********************************************************************

      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
      PARAMETER(D0=0.0D0)
C
      CHARACTER WFTYP*(*)
      LOGICAL   WFCONV
C
#include "maxorb.h"
#include "dcbidx.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbopt.h"
#include "dcbgen.h"
!dcbgen.h needed for DOLVC for LVCORR
#include "dcbham.h"
#include "frame.h"
!frame.h needed for CORRLV for LVCORR
C
      LOGICAL EX
      INTEGER RENAME
      real(8), allocatable :: WORK(:)
C     JVP-2020 modification to save moltra input
C     variables to save 
      CHARACTER*72 SAVTRA(4,2,2)
      INTEGER SAVINT(5)
      LOGICAL SAVLOG(10)
      REAL*8  SAVREAL(2)
C
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in PSIOPT')
C
      CALL QENTER('PSIOPT')
C
C     JVP-2020, save moltra inp info
      CALL SAVE_MOLTRA_INP(SAVTRA,SAVLOG,SAVREAL,SAVINT,IPRINT)
C
C     Determine the wave function type
C     --------------------------------
C
      IF (WFTYP .EQ. 'SCF  ') THEN
         IF (NASHT .EQ. 0) THEN
            MCTYPE = JDHF
         ELSE IF (NASHT .EQ. 1 .AND. NAELEC .EQ. 1) THEN
            MCTYPE = JODHF
         ELSE
            CALL QUIT('2. order optimization of this open shell DHF'//
     &        ' is not implemented')
         END IF
C
         IF (THQKVA .LE. D0) THQKVA = 0.1D0
C        ... default quadratic conv. factor for SCF
C            (not set in KRMCINP because we don't know if SCF there)
         IF (IKRMCCNO .EQ. -3) THEN
C           ... no user input, default for SCF: transform to canonical orbitals
            NATONL = .FALSE.
            FOCKON = .TRUE.
            IKRMCCNO = -1
C           ... code for no RTRACI (no RTRACI for SCF)
         END IF
C
C        Delete DFFCK2 if existing, because we will change CMO here ...
C
         LUBUF = 1
         INQUIRE(FILE='DFFCK2',EXIST=EX)
         IF(EX) THEN
            CALL OPNFIL(LUBUF,'DFFCK2','OLD','PSIOPT')
            WRITE(LUPRI,'(/A)') '* PSIOPT-QCDHF: '//
     &      'deleting DFFCK2 file with old two-electron Fock matrix'
            CLOSE(LUBUF,STATUS='DELETE')
         END IF
         INQUIRE(FILE='DFDENS',EXIST=EX)
         IF(EX) THEN
            CALL OPNFIL(LUBUF,'DFDENS','OLD','PSIOPT')
            WRITE(LUPRI,'(/A)')
     & '* PSIOPT-QCDHF: deleting DFDENS file with old density matrix'
            CLOSE(LUBUF,STATUS='DELETE')
         END IF
C
      ELSE IF (NGAS_DC .EQ. 1) THEN
         MCTYPE = JCAS
      ELSE
         MCTYPE = JGAS
      END IF
      IF (THQKVA .LE. D0) THQKVA = 8.0D00
C     ... default quadratic conv. factor for MCSCF

!     set logical flag on common block in dcbopt.h for spinfree MCSCF
#ifdef MOD_MCSCF_spinfree
!     if(spinfr.or.levyle)then
      if(spinfr)then
        spinfr_krmc = .true.
      else if(levyle)then
!       fixme: already the initial CI using Levy-Leblond yields
!              wrong results compared to Dalton. S. Knecht - Aug 2010
        write(lupri,*) ' *** error in KRMCSCF: MCSCF with
     &            Levy-Leblond Hamiltonian is not working.***'
        call quit('*** error in KRMCSCF: MCSCF with '//
     &            'Levy-Leblond Hamiltonian is not working.***')
      else
        spinfr_krmc = .false.
      end if
#else
        spinfr_krmc = .false.
#endif
cdebug
c     opt_ucibos = .false.
cdebug
C
C     Say Hi!
C     -------
C
      CALL KRMCHI
C
C     ***************************************************
C     *** Open files needed *****************************
C     ***************************************************
C     *** Set logical units ***
C
      CALL SETKRMC(LUPRI)

C
C     ***************************************************
C     *** Initialize configurational and orbital      ***
C     *** parameters.                                 ***
C     ***************************************************
C
C
      CALL SETDCBIDX()
      CALL RSETWOP(WFTYP,WORK,KFREE,LFREE)
      CALL FLSHFO(LUPRI)
C
#if defined (VAR_MPI2)
C
C     ***************************************************
C     *** Initialize common block information         ***
C     ***************************************************
!     1. mpi
      CALL MPIXCALC_TYPESZ()
C
#endif

!     2. lvcorr correction
      if(dolvc.and.((iopt_intdef.eq.3).or.(onecap.and.intv1c.eq.2)))then
        kfree_save = kfree
        call memget2('REAL','lvc',klvc,n2bbasxq*(nopen_mc+1),
     &     work,kfree,lfree)
        call lvcorr(work(klvc),work(kfree),lfree,ipropt)
        call memrel('lvcorr',work,1,kfree_save,kfree,lfree)
!       set lvcorr value (otherwise it was initialized to 0.0d0)
        E_lvcorr_mc = corrlv
      else
        E_lvcorr_mc = 0.0D0
      end if
C
C     ***************************************************
C     *** Memory allocation *****************************
C     ***************************************************
C
C
C hjaaj + sk 04-jul-07: new memory allocation scheme when calling
C                       KR-MCSCF with .CIONLY .or.
C                                     MOPT_MXMACRO .lt. 0
C
      CALL MEMGET2('REAL','CMO',KCMO  ,NCMOTQ        ,WORK,KFREE,LFREE)
      IF (MOPT_MXMACRO .LT. 0) THEN
         CALL MEMGET2('REAL','DV'  ,KDV   ,0          ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','FC'  ,KFC   ,0          ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','FV'  ,KFV   ,0          ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','FQ'  ,KFQ   ,0          ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','PV'  ,KPV   ,0          ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','GCI' ,KGCI  ,0          ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','GOEE',KGOEE ,0          ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','GOEP',KGOEP ,0          ,WORK,KFREE,LFREE)
         IF (OPT_CIPROGRAM .EQ. 'LUCIAREL') THEN
            CALL MEMGET2('REAL','FCAC',KFCAC ,0       ,WORK,KFREE,LFREE)
            CALL MEMGET2('REAL','H2AC',KH2AC ,0       ,WORK,KFREE,LFREE)
            CALL MEMGET2('REAL','CREF',KCREF ,0       ,WORK,KFREE,LFREE)
         ELSE
            CALL MEMGET2('REAL','FCAC',KFCAC ,N2ASHXQ ,WORK,KFREE,LFREE)
            CALL MEMGET2('REAL','H2AC',KH2AC ,NASHT*NASHT*NNASHX*NZ*3,
     &           WORK,KFREE,LFREE)
            CALL MEMGET2('REAL','CREF',KCREF ,NZCONFQ ,WORK,KFREE,LFREE)
         END IF
      ELSE
         CALL MEMGET2('REAL','DV'  ,KDV   ,N2ASHXQ    ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','FC'  ,KFC   ,N2ORBXQ    ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','FV'  ,KFV   ,N2ORBXQ    ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','FQ'  ,KFQ   ,NZ*NASHT*NORBT,
     &        WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','PV'  ,KPV   ,NASHT*NASHT*NNASHX*NZ*3,
     &        WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','FCAC',KFCAC ,N2ASHXQ ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','H2AC',KH2AC ,NASHT*NASHT*NNASHX*NZ*3,
     &        WORK,KFREE,LFREE)
c hjaaj 03-jul-07: TODO : implement buffered i/o of cref and gci
c        where nzconfq_buf is max size in memory
c        for MPI we could perhaps try NZCONFQ_BUF = NZCONFQ / NMPROC
c        CALL MEMGET2('REAL','CREF',KCREF ,NZCONFQ_BUF,WORK,KFREE,LFREE)
c        CALL MEMGET2('REAL','GCI ',KGCI  ,NZCONFQ_BUF,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','CREF',KCREF ,NZCONFQ    ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','GCI ',KGCI  ,NZCONFQ    ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','GOEE',KGOEE ,NZHOPEQ    ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','GOEP',KGOEP ,NZXOPPQ    ,WORK,KFREE,LFREE)
      END IF

!     implement boson symmetry in MCSCF - enabling spinfree MCSCF
!     sk + hjaaj 10-aug-10
      call memget2('INTE','ibeig',kibeig,norbt,work,kfree,lfree)
C
C     ***************************************************
C     *** Call 2nd order optimization *******************
C     ***************************************************
C
      CALL PSIOP1(WORK(KCMO),WORK(KDV),WORK(KPV),WORK(KFC),
     &            WORK(KFV),WORK(KFQ),WORK(KFCAC),WORK(KCREF),
     &            WORK(KH2AC),WORK(KGCI),WORK(KGOEE),WORK(KGOEP),
     &            WORK(KIBEIG),WFCONV,WORK,KFREE,LFREE)
C
C
C     ***************************************************
C     *** Call output routines **************************
C     ***************************************************
C
C
      IF (MOPT_MXMACRO .GT. 0) CALL KRMCOUT(WORK,KFREE,LFREE)
C
C
C
C     ***************************************************
C     *** Memory deallocation ***************************
C     ***************************************************
C
      CALL MEMREL('PSIOPT',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     ***************************************************
C     *** Close files ***********************************
C     ***************************************************
C
C
      CLOSE(LUKRMC,STATUS='KEEP')
      CLOSE(LUKRM1,STATUS='DELETE')
      CLOSE(LUKRM2,STATUS='DELETE')
      CLOSE(LUKRM3,STATUS='KEEP')
c     CLOSE(LUKRM3,STATUS='DELETE')
      CLOSE(LUKRM4,STATUS='DELETE')
      IF (LUKRM5.GT.0) CLOSE(LUKRM5,STATUS='DELETE')
      CLOSE(LUITFO,STATUS='DELETE')
C
C     JVP-2020, restore MOLTRA input info
      CALL RESET_MOLTRA_INP(SAVTRA,SAVLOG,SAVREAL,SAVINT,IPRINT)
C
      call dealloc(WORK)
      CALL QEXIT('PSIOPT')
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiopt */
      SUBROUTINE PSIOP1(CMO,DV,PV,FC,FV,FQ,FCAC,CREF,H2AC,GCI,GOEE,GOEP,
     &                  IBEIG,WFCONV,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get second-order optimized wave function.
C     Either DF or KR-MCSCF...
C
C     Written by J. Thyssen and H. J. Aa. Jensen - Oct 22 1998
C     Last revision :
C
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "maxorb.h"
#include "dcbnrt.h"
#include "dcbpsi.h"
#include "dcbkrmc_itinfo.h"
#include "dcbgen.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "frame.h"
C
      LOGICAL   WFCONV, MOISNO
      DIMENSION WORK(*)
      DIMENSION CMO(*),DV(*),FC(*),FV(*),FQ(*),FCAC(*),CREF(*),
     &          H2AC(*),GCI(*),GOEE(*),GOEP(*), PV(*), IBEIG(*)
C
      LOGICAL CHKSTP,DO4INDEX,TOBE
      LOGICAL NRALW, NEOALW
C
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
      CHARACTER SZTRA(2,4)*200
#ifdef mcscf_fast_debug
      logical blabla
#endif
C
C
      CALL QENTER('PSIOP1')
      CALL GETTIM(CPUTOT1,WALLTOT1)
      KFRSAV = KFREE
C
C
C     ***************************************************
C     *** Initialize some parameters ********************
C     ***************************************************
C
      DO I = 1,2
         IOPT_ISYMOP(I) = 1
         IOPT_IFCKOP(I) = 1
         IOPT_IHRMOP(I) = 1
      END DO
      ITMAC  = 0
      ITBCK  = 0
      DO4INDEX = ( MCTYPE .GE. JMCMIN )
C     ... MCTYPE .lt. JMCMIN is SCF
      EMCOLD = D0
      DEPRED = D0
      STPLEN = D0
      WFCONV = .FALSE.
C
      NRALW = FLAG(1)
      NEOALW = FLAG(2)
      DONR = NRALW
      FROZCI = OPT_NOCI
C
C     KTRLVL = 3 (KTRLVL set in *OPTIMIZE)
      JTRLVL = 0
C
      CALL DZERO(CPUOPT,NOPTTIM)
      CALL DZERO(WLLOPT,NOPTTIM)
C
      CALL DZERO(GNORM,5)
C
C     Integrals in macro iteration zero.
C
      CALL INTCON(IOPT_INTFLG,IOPT_INTBUF,IOPT_INTDEF,GNORM(5),
     &            OPT_CNVINT,ITMAC,IOPT_ITRINT,OPT_INTTYP)
C
C
C     ***************************************************
C     *** Get start vectors *****************************
C     ***************************************************
C
C
#ifdef UNDEF
      write(6,*) '(1) iopt_intflg=',iopt_intflg
      write(6,*) '(1) iopt_intbuf=',iopt_intflg
      write(6,*) '(1) iopt_intdef=',iopt_intflg
      write(6,*) '(1) gnorm(5)=',gnorm(5)
      write(6,*) '(1) opt_cnvint=',opt_cnvint
      write(6,*) '(1) itmac=',itmac
      write(6,*) '(1) iopt_itrint=',iopt_itrint
      write(6,*) '(1) opt_inttyp=',opt_inttyp
#endif
C
C     Get start orbitals and CI vectors:
C
      CALL ROPTST(IOPTST,3,CMO,CREF,IBEIG,IPROPT,WORK,KFREE,LFREE)
C     CALL ROPTST(ISTART,IWRT,CMO,CREF,IBEIG,IPRINT,WORK,KFREE,LFREE)
C
C     Only integral transformation done
C
      IF (KTRLVL.eq.5) GOTO 3000
C
C     Check if only do CI (no MCSCF iterations)
C
      IF (MOPT_MXMACRO.LE.0.OR.NCIROOT.GE.2) THEN
         IF (MOPT_MXMACRO.LE.0) THEN
            WRITE(LUPRI,'(/A)')
     &         ' Stop after initial CI in KR-MCSCF as requested'//
     &         ' with .MAX MACRO'
         ELSE
            WRITE(LUPRI,'(/A,I5,A)')
     &              ' (ROPTST)  Requested multi-root MCSCF aborted.'
         END IF
         GO TO 3000
      END IF
C
C     Initialize some more parameters:
C
C     If DEPRED < 0, we have restarted a calculation.
C
      IF (DEPRED .NE. D0) THEN
         CHKSTP = .TRUE.
         ISTEP  = -1
      ELSE
         CHKSTP = .FALSE.
         ISTEP  = 0
      END IF
C
C     Check for legal values of Beta:
C
      BETA = ABS(BETA)
      IF (BETA .LE. THRZER)  BETA = D1
      IF (BETMIN .LE. D0 .OR.
     &    BETMIN .GT. BETMAX) BETMIN = D1
      IF (BETMAX .LE. D1)     BETMAX = D1
C
C
C     ***************************************************
C     *** Start of macro iteration loop *****************
C     ***************************************************
C
C
  100 CONTINUE
         CALL GETTIM(CPUMAC1,WALLMAC1)
         CALL DZERO(DINFO,LDINFO)
         CALL IZERO(IINFO,LIINFO)
C
C        Write iteration number etc.
C        ---------------------------
C
         IF (ITBCK .GT. 0) THEN
            WRITE(LUPRI,9001) ITMAC, ITBCK
         ELSE IF (ITABS .GT. 0) THEN
            WRITE(LUPRI,9002) ITMAC
         ELSE
            ITMAC = ITMAC+1
            ITMACN = ITMACN + 1
            WRITE(LUPRI,9000) ITMAC
         END IF
C
 9000    FORMAT(//' <<< MACRO ITERATION',I3,' >>>'
     &           /' ---------------------------')
 9001    FORMAT(//' <<< MACRO ITERATION',I3,' BACKSTEP NO.',I3,' >>>'
     &           /' ---------------------------')
 9002    FORMAT(//' <<< MACRO ITERATION',I3,' AFTER ABSORPTION >>>'
     &           /' ---------------------------')
         CALL FLSHFO(LUPRI)
C
C        ***************************************************
C        *** Which integrals to calculate ******************
C        ***************************************************
C
         CALL INTCON(IOPT_INTFLG,IOPT_INTBUF,IOPT_INTDEF,GNORM(5),
     &        OPT_CNVINT,ITMAC,IOPT_ITRINT,OPT_INTTYP)
C
         WRITE(LUPRI,9020) OPT_INTTYP
C
 9020    FORMAT(/' (PSIOPT) Integrals in this macro iteration: ',A)
C
C        ***************************************************
C        *** Get MO-coefficients ***************************
C        ***************************************************
C
C
         JRDCMO = 0
         CALL RREADMO(CMO,JRDCMO,1,LUKRMC)
C        read orbital symmetry info
         call izero(ibeig,norbt)
         if(spinfr_krmc)then
           call ireakrmc(lukrmc,'IBEIG   ',ibeig,norbt)
         else if(opt_chckjz)then
           call ireakrmc(lukrmc,'MJVEC   ',ibeig,norbt)
         end if
C
C
C        ***************************************************
C        *** Transform two-electron integrals **************
C        ***************************************************
C
         CALL GETTIM(CPUITR1,WALLITR1)
         IF ( DO4INDEX ) THEN
C
C           Do 4-index transformation.
C           KTRLVL is the order of the transformation.
C
            IF (.NOT. COMPROT) THEN
C
               JTRLVL = KTRLVL
C
            ELSE
C
C              Orbital compression: we need (ga|aa) integrals to
C              get gradient correct. Later, we'll calculate the
C              rest of the integrals.
C
               JTRLVL = 1
C
            END IF
            CALL INISZT(JTRLVL,SZTRA,LNOROT,NOROT)
C
C           Get all coefficients
C
            CALL MEMGET2('REAL','CMO',KCMO,NCMOTQ,WORK,KFREE,LFREE)
            CALL MEMGET2('INTE','IBEIG',KIBEIG,NORBT,WORK,KFREE,LFREE)
C
            CALL IZERO(WORK(KIBEIG),NORBT)
C
C           Read coefficients from KRMCSCF
C
            INQUIRE(FILE='KRMCSCF',EXIST=TOBE)
            IF (.NOT.TOBE) THEN
               WRITE(LUPRI,'(/A)')
     $        '*** ERROR in PSIOP1 *** File with coefficients not found'
               CALL QUIT('*** ERROR in PSIOP1   ***')
            ELSE
               JRDMO = 1
               CALL RREADMO(WORK(KCMO),JRDMO,1,LUKRMC)

               if(spinfr_krmc)then
                 call ireakrmc(lukrmc,'IBEIG   ',work(kibeig),norbt)
               else if(opt_chckjz)then
                 call ireakrmc(lukrmc,'MJVEC   ',work(kibeig),norbt)
               end if

            END IF
C
C           *** Now do the transformation
C
#ifdef mcscf_fast_debug
            inquire(file="4INDINFO", exist=blabla)
            if(blabla)then
              print*,'skip RTRACTL ... ;)'
            else
#endif
              CALL RTRACTL(WORK(KCMO),work(kibeig),SZTRA,
     &                     WORK(KFREE),LFREE)
#ifdef mcscf_fast_debug
            end if
#endif
C
C           *** Throw away CMO again
C
            CALL MEMREL('PSIOP1.RTRACTL.CMO',WORK,1,KCMO,KFREE,LFREE)
C
         END IF
         CALL GETTIM(CPUITR2,WALLITR2)
C
C
C
C        ***************************************************
C        *** Get Gradient **********************************
C        ***************************************************
C
C
         CALL RGRAD(CMO,CREF,ibeig,
     &              get_orbital_rotation_indices_pp(),
     &              get_orbital_rotation_indices_pn(),
     &              GCI,GOEE,GOEP,EMY,EACTIV,DV,PV,FC,FV,FQ,
     &              FCAC,H2AC,IPROPT,WORK,KFREE,LFREE)
C
         GNRMSV = GNORM(5)
         CALL RGRDINF(GCI,GOEE,GOEP)
C
C
         EMCSCF = POTNUC + EMY + EACTIV + E_LVCORR_MC
C
         if(dabs(E_lvcorr_mc).gt.0.0d0)then
           WRITE(LUPRI,9005) POTNUC,E_LVCORR_MC,EMY,EACTIV,EMCSCF
 9005      FORMAT(/' Nuclear repulsion            : ',F25.15,/,
     &             ' LVCORR correction            : ',F25.15,/,
     &             ' Inactive energy              : ',F25.15,/,
     &             ' Active energy                : ',F25.15,/,
     &             ' Total KR-MCSCF energy        : ',F25.15)
         else
           WRITE(LUPRI,9006) POTNUC,EMY,EACTIV,EMCSCF
 9006      FORMAT(/' Nuclear repulsion            : ',F25.15,/,
     &             ' Inactive energy              : ',F25.15,/,
     &             ' Active energy                : ',F25.15,/,
     &             ' Total KR-MCSCF energy        : ',F25.15)
         end if
C
         WRITE(LUPRI,9010) (ITMAC,GNORM(I),I=1,5)
C
 9010    FORMAT(/' CI gradient        (It. ',I3,') : ',F25.15,/,
     &           ' Orbital e-e grad.  (It. ',I3,') : ',F25.15,/,
     &           ' Orbital e-p grad.  (It. ',I3,') : ',F25.15,/,
     &           ' Total orbital grad.(It. ',I3,') : ',F25.15,/,
     &           ' Total gradient     (It. ',I3,') : ',F25.15,/)
C
         IF (IPROPT .GE. 1) THEN
              IDENSLR_STATE = 1
              NATOLCR       = .TRUE.
              IPRNO = 3*IPROPT
              CALL MEMGET2('REAL','OCCNO',KOCCNO,NORBT  ,
     &             WORK,KFREE,LFREE)
              CALL MEMGET2('REAL','UNO  ',KUNO  ,N2ASHXQ,
     &             WORK,KFREE,LFREE)

C             erase old occupation numbers
              CALL DZERO(WORK(KOCCNO),NORBT)
              CALL DZERO(WORK(KUNO),N2ASHXQ)

              CALL RGETNO(CREF,WORK(KOCCNO),WORK(KUNO),CMO,
     &                    .FALSE.,MOISNO,.FALSE.,IPRNO,WORK,KFREE,LFREE)
              CALL MEMREL('PSIOP1.after RGETNO',WORK,1,KFRSAV,
     &                    KFREE,LFREE)
              IF ( OPT_CIPROGRAM .EQ. 'LUCIAREL' ) THEN
                IF( JKRRUNTYPE .lt. 2 ) NATOLCR = .FALSE.
              END IF
         END IF
C
         CALL FLSHFO(LUPRI)
C
C        ***************************************************
C        *** Converged ??? *********************************
C        ***************************************************
C
         IF (GNORM(5) .LE. OPT_THRGRD) THEN
            ISTEP = -2
            WFCONV = .TRUE.
         END IF
C
C
C        ***************************************************
C        *** Check step ************************************
C        ***************************************************
C
C        Input to DIRSTP:
C        istep = -2: MC converged, get step info for analysis
C        istep = -1: restart prediction (currently unused)
C        istep =  0: normal second-order step check
C
C        DIRSTP returns:
C        istep =  0: Step is OK
C        istep =  1: step too large, do backstep
C        istep =  2: step rejected, no backstep available
C        istep = -1: no check because we are close to convergence.
C
         IF (CHKSTP) THEN
            CALL DIRSTP(ISTEP,CMO,CREF,WORK,KFREE,LFREE)
         ELSE
            EMCOLD = EMCSCF
            ISTEP  = 0
         END IF
C
         DINFO(9) = RTRUST
C
C        Exit if converged
C
         IF (WFCONV) GOTO 1000
C
         IF (ISTEP .EQ. 2) THEN
            WRITE(LUPRI,9050)
            GOTO 1999
         END IF
C
 9050    FORMAT(//' *** WARNING *** ',
     &        'Aborting: step rejected and backstep is NOT available')
C
         IF (ISTEP .GT. 0) THEN
            ISTEP = 0
            DONR  = NRALW
            ITBCK = ITBCK + 1
c           call quit('argh')
            IF (ITBCK .LE. MOPT_MAXBCK) THEN
               WRITE(LUPRI,9060)
            ELSE
               WRITE(LUPRI,9062) MOPT_MAXBCK
            END IF
            GO TO 1000
         ELSE IF (ISTEP .LT. 0) THEN
            IF (GNORM(5) .GT. GNRMSV) THEN
               WRITE(LUPRI,9064)
               GOTO 1999
            END IF
            ITBCK = 0
         ELSE
            ITBCK = 0
         END IF
C
 9060    FORMAT(/' (PSIOP1) Optimization control: ',
     &        'Step was too large, we backstep')
 9062    FORMAT(/' (PSIOP1) Optimization control: ',
     &        '*** ERROR in PSIOP1 *** No convergence',
     &        //' Step was too large and maximum number of ',
     &        'backsteps',I2,' reached'/)
 9064    FORMAT(/' (PSIOP1) Optimization control: ',
     &        '*** ERROR in PSIOP1 *** ',
     &        'Gradient norm increasing in local region.',
     &        /' Probable cause: numerical round-off errors ',
     &        'or too many integrals neglected/screening away')
C
         CHKSTP = .TRUE.
C
C
C        ***************************************************
C        *** Second order step *****************************
C        ***************************************************
C
C
         STPLEN = D0
         STPLCI = D0
         STPLEE = D0
         STPLEP = D0
         DEPRED = D0
c        IF (DONR) THEN
            CALL DIRNR(CREF,GCI,GOEE,GOEP,CMO,IBEIG,DV,FC,FV,FCAC,
     &           H2AC,WORK,KFREE,LFREE)
            LDIRNR = .TRUE.
c        ELSE
c           CALL QUIT('*** ERROR in PSIOP1 ***: NEO not implemented')
C           LDIRNR = .FALSE.
c        END IF
C
C
C
 1000    CONTINUE
         CALL GETTIM(CPUMAC2,WALLMAC2)
         CPUTID = SECTID(CPUMAC2-CPUMAC1)
         WALLTID = SECTID(WALLMAC2-WALLMAC1)
         WRITE(LUPRI,9100) ITMAC,CPUTID,WALLTID
C
C        Save information for final summary output
C
         IINFO(1) = ITMAC
         DINFO(1) = EMY
         DINFO(2) = EACTIV
         DINFO(3) = EMCSCF
         DINFO(16) = CPUMAC2-CPUMAC1
         DINFO(18) = CPUITR2-CPUITR1
C
         WRITE(LUITFO) DINFO,IINFO
C
 9100    FORMAT(/'>>> CPU (WALL) TIME FOR MACRO IT.',
     &        I3,' : ',A,' (',A,')')
         CALL FLSHFO(LUPRI)
         IF (WFCONV) GOTO 2000
         IF ( (ITBCK .EQ. 0 .AND. ITMAC .LT. MOPT_MXMACRO) .OR.
     &        (ITBCK .GT. 0 .AND. ITBCK .LE. MOPT_MAXBCK ) ) THEN
            IOPT_INTBUF = IOPT_INTFLG
            GOTO 100
         END IF
C
C
C        ***************************************************
C        *** End of macro iteration loop *******************
C        ***************************************************
C
C
C     ***************************************************
C     *** Not converged *********************************
C     ***************************************************
C
 1999 CONTINUE
      WRITE(LUPRI,9200)
 9200 FORMAT(/'WARNING: Wave function did not converge...')
      GOTO 3000
 2000 CONTINUE
C
C
C     ***************************************************
C     *** We have converged *****************************
C     ***************************************************
C
C
C
 3000 CONTINUE
      CALL GETTIM(CPUTOT2,WALLTOT2)
      CPUTID = SECTID(CPUTOT2-CPUTOT1)
      WALLTID = SECTID(WALLTOT2-WALLTOT1)
      WRITE(LUPRI,9300) CPUTID,WALLTID

C     calculate 1-particle density to be exported on file - default
      if(.not.no1pdens_save)then
        WRITE(LUPRI,'(//A/)') ' INFO: dumping 1-particle density'//
     &  ' matrix calculated for the (converged) MCSCF wave function'//
     &  ' to file.'
        CALL MEMGET2('REAL','OCCNO',KOCCNO,NORBT  ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','UNO'  ,KUNO  ,N2ASHXQ,WORK,KFREE,LFREE)
        CALL RGETNO(CREF,WORK(KOCCNO),WORK(KUNO),CMO,
     &             .false.,MOISNO,.true.,1,WORK,KFREE,LFREE)
        CALL MEMREL('PSIOP1.after RGETNO',WORK,1,KFRSAV,KFREE,LFREE)
      end if
C
 9300 FORMAT(/'>>> TOTAL CPU (WALL) TIME IN KR-MCSCF: ',A,'(',A,')')
C
      CALL MEMREL('PSIOP1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('PSIOP1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgrad */
      SUBROUTINE RGRAD(CMO,CREF,ibeig,JXOPE,JXOPP,
     &                 GCI,GOEE,GOEP,EMCMY,EMCACT,DV,PV,FC,FV,FQ,FCAC,
     &                 H2AC,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate the gradient.
C
C     Input:
C        CMO     - the orbitals for the CEP
C        CREF    - CI vector
C        ibeig   - orbital symmetry info
C        JXOPE   - orbital e-e rotations
C        JXOPP   - orbital e-p rotations
C
C
C     Output:
C        EMCACT  - active energy
C        GCI     - the CI gradient of the CEP
C        GOEE    - the e-e orbital gradient of the CEP
C        GOEP    - the e-p gradient of the CEP
C        EMCMY   - inactive energy
C        DV      - one-electron active density matrix
C        PV      - two-electron active density matrix
C        H2AC    - two-electron integrals with active indices
C        FC      - the inactive Fock matrix
C        FCAC    - the inactive Fock matrix with active indices
C        FV      - the active Fock matrix
C        FQ      - the Q Fock matrix
C
C
C     Note:
C        default: COMPROT == .FALSE.
C
C     Written by J. Thyssen and H. J. Aa. Jensen - Oct 22 1998
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
      DIMENSION GCI(*),DV(*),FC(*),FV(*),FQ(*),PV(*)
      DIMENSION JXOPE(*),JXOPP(*),GOEE(*),GOEP(*)
      DIMENSION CREF(*), ibeig(*)
      DIMENSION WORK(*)
C
C     Local variable
C
      DIMENSION DFAC(2)
C
      LOGICAL FND,FNDLAB,SAVEDFGEN
      LOGICAL dofq, doh2ac

C
#include "dcbopt.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CALL QENTER('RGRAD')
      CALL RTKTIME(.TRUE.,9)
      KFRSAV = KFREE
C
      WRITE(LUPRI,'(/A)') ' (RGRAD) Calculating gradient.'
C
      DFAC(1) = D1
      DFAC(2) = DM1
C
C     ***************************************************
C     *** Clean file KRMC FOCK **************************
C     ***************************************************
C
      REWIND LUKRM3
      FND =  FNDLAB('SODLABEL',LUKRM3)
      BACKSPACE LUKRM3
      CALL NEWLAB('SODLABEL',LUKRM3,LUPRI)
      REWIND(LUKRM3)
C
C     ***************************************************
C     *** Calculate active density matrices *************
C     ***************************************************
C
      IF (NASHT .NE. 0) THEN
         CALL RMAKDM(.TRUE.,.TRUE.,DV,PV,CREF,WORK,KFREE,LFREE)
      END IF
C
C
C     ***************************************************
C     *** Calculate inactive and active Fock matrices ***
C     ***************************************************
C
      CALL RFCKMAT(DV,CMO,EMCMY,EMCACT,FC,FV,WORK(KFREE),LFREE)
C
C     ***************************************************
C     *** Calculate FQ matrix and two-electron inte-  ***
C     *** grals with active orbital indices.          ***
C     ***************************************************
C
      IF (NASHT .GE. 1) THEN
C
C        Get Fock matrix with both indices active.
C        FCAC is written to KRMC FOCK.
C
         CALL RGETAC(FC,FCAC,IPRINT)
C
      END IF
C
      IF (NASHT .GT. 1 .OR. MCTYPE .GE. JMCMIN) THEN
C
C        Calculate FQ matrix and two-electron integrals
C        with active orbital indices (H2AC).
C        FQ and H2AC are written to KRMC FOCK.
C
         dofq   = .true.
         doh2ac = .true.
         CALL RGETH2(FQ,H2AC,PV,ibeig,dofq,doh2ac,
     &               .true.,WORK(KFREE),LFREE)
C
C        Active energy.
C
         EMCACT = RENRACT(DV,FCAC,FQ,IPROPT)
C
      ELSE IF (NASHT. EQ. 1) THEN
         CALL DZERO(FQ,NZ*NASHT*NORBT)
         EMCACT = RENRACT(DV,FCAC,FQ,IPROPT)
      ELSE
         EMCACT = D0
      END IF
C
C
C     ***************************************************
C     *** Calculate Fock matrix contributions to the  ***
C     *** diagonal of the orbital Hessian             ***
C     ***************************************************
C
C
      IF (NZXOPT .GT. 1) THEN
c        CALL RDIAG
      END IF
C
C
C     ***************************************
C     *** Write matrices to file - part 1 ***
C     ***************************************
C
C     DV & PV
C     -------
C
C     Save active density matrices on file KRMC FOCK
C
      CALL WRTKRMC(LUKRM3,'DVMO    ',DV,N2ASHXQ)
      IF ( NASHT .GT. 1 .OR. MCTYPE .GE. JMCMIN ) THEN
         CALL WRTKRMC(LUKRM3,'PVMO    ',PV,NASHT*NASHT*NNASHX*NZ*3)
      END IF
C
C     ***************************************************
C     *** Calculate the orbital gradient ****************
C     ***************************************************
C
      SAVEDFGEN = COMPROT
      IF (NZHOPE .GT. 0) THEN
         CALL DZERO(GOEE,NZHOPEQ)
         LFPV = NZ*NASHT*NORBT
         CALL MEMGET2('REAL','FPV',KFPV,LFPV,WORK,KFREE,LFREE)
         CALL RORBGRD(.TRUE.,DV,FC,FC,FV,FQ,GOEE,
     &                JXOPE,NZHOPE,WORK(KFPV),DFAC,
     &                IPQTOQ(1,0),IPRINT)
C
C        Save generalized Fock matrix (F_{pv})
C
         IF (.NOT. SAVEDFGEN) THEN
            CALL WRTKRMC(LUKRM3,'FGENPV  ',WORK(KFPV),
     &           NORBT * NASHT * NZ)
            SAVEDFGEN = .TRUE.
         END IF
C
         CALL MEMREL('RORBGRD.KFPV',WORK,KFRSAV,KFPV,KFREE,LFREE)
      END IF
      IF (NZXOPP .GT. 0) THEN
         CALL DZERO(GOEP,NZXOPPQ)
         LFPV = NZ*NASHT*NORBT
         CALL MEMGET2('REAL','FPV',KFPV,LFPV,WORK,KFREE,LFREE)
         CALL RORBGRD(.TRUE.,DV,FC,FC,FV,FQ,GOEP,
     &                JXOPP,NZXOPP,WORK(KFPV),DFAC,
     &                IPQTOQ(1,0),IPRINT)
C
C        Save generalized Fock matrix (F_{pv})
C
         IF (.NOT. SAVEDFGEN)
     &        CALL WRTKRMC(LUKRM3,'FGENPV  ',WORK(KFPV),
     &        NORBT * NASHT * NZ)
         CALL MEMREL('RORBGRD.KFPV',WORK,KFRSAV,KFPV,KFREE,LFREE)
      END IF
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('rGRAD: gradient:',-1)
         CALL RPRKAP(GOEE,GOEP,JXOPE,NZHOPE,
     &               JXOPP,NZXOPP,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C
C     ***************************************************
C     *** Singular-value decomposition of gradient ******
C     ***************************************************
C
      IF (COMPROT) THEN
C
         CALL SVLDEG(GOEE,GOEP,JXOPE,NZHOPE,JXOPP,NZXOPP,
     &        DV,FC,FV,FQ,CMO,IPRINT,WORK(KFREE),LFREE)
C
      END IF
C
C
C     ***************************************************
C     *** Calculate the CI gradient *********************
C     ***************************************************
C
      IF (NZCONF .GT. 0 .AND. .NOT. FROZCI ) THEN
C
C        Calculate CI gradient
C
         CALL DZERO(GCI, NZCONFQ)
         CALL RCIGRAD(CREF,FCAC,H2AC,GCI,EMCACT,WORK,KFREE,LFREE)
      ELSE
         CALL DZERO(GCI, NZCONFQ)
      END IF
C
C
C     ***************************************
C     *** Write matrices to file - part 2 ***
C     ***************************************
C
!     --------------------  start of new part 1 ------------------
!      SK - Jan 2010: matrices are saved on file above
!      because we might need them in the mj-blocking of the
!      orbital rotations.
!
!C     DV & PV
!C     -------
!C
!C     Save active density matrices on file KRMC FOCK
!C
!      CALL WRTKRMC(LUKRM3,'DVMO    ',DV,N2ASHXQ)
!      IF ( NASHT .GT. 1 .OR. MCTYPE .GE. JMCMIN ) THEN
!         CALL WRTKRMC(LUKRM3,'PVMO    ',PV,NASHT*NASHT*NNASHX*NZ*3)
!      END IF
!     --------------------  end of new part 1 --------------------
C
C     FC & FV
C     -------
C
      CALL WRTKRMC(LUKRM3,'FCMO    ',FC,N2ORBXQ)
      IF (NASHT .GT. 0) CALL WRTKRMC(LUKRM3,'FVMO    ',FV,N2ORBXQ)
C
C     FCAC, FQ & H2AC
C     ---------------
C
      IF (MCTYPE .GE. JMCMIN) THEN
         CALL WRTKRMC(LUKRM3,'FCAC    ',FCAC,N2ASHXQ)
         CALL WRTKRMC(LUKRM3,'FQMO    ',FQ  ,NASHT*NORBT*NZ)
         CALL WRTKRMC(LUKRM3,'H2AC    ',H2AC,NASHT*NASHT*NNASHX*NZ*3)
      END IF
C
C     CMO & CREF
C     ----------
C     (we have to read & write CREF again as CREF is placed after NEWORB
C      and is therefore deleted when rewriting CMO
C
      IF (COMPROT) THEN
         IF (NZCONF .GT. 0) CALL REAKRMC(LUKRMC,'CREF    ',CREF,NZCONFQ)
         CALL WRTKRMC(LUKRMC,'NEWORB  ',CMO,NCMOTQ)
         IF (NZCONF .GT. 0) CALL WRTKRMC(LUKRMC,'CREF    ',CREF,NZCONFQ)
      END IF
C
C
      CALL MEMREL('RORBGRD',WORK,1,KFRSAV,KFREE,LFREE)
C
C
      CALL FLSHFO(LUPRI)
      CALL RTKTIME(.FALSE.,9)
      CALL QEXIT('RGRAD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rmakdm */
      SUBROUTINE RMAKDM(CALCDV,CALCPV,DV,PV,CREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Construct one- and two-electron active density matrix
C
C     Input:
C        CMO     - the orbitals for the CEP
C        CREF    - CI vector
C
C     Output:
C        DV      - one electron active density matrix
C        PV      - two electron active density matrix
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision : 24.01.2001, T. Fleig
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
!sk   PARAMETER (D1 = 1.0D00, D0 = 0.0D00, DEPS = 1.0D-8 )
      PARAMETER (D1 = 1.0D00, D0 = 0.0D00, DEPS = 1.0D-7 )
C
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbopt.h"
#include "dgroup.h"
C
      LOGICAL CALCDV, CALCPV
      DIMENSION DV(*), PV(*), CREF(*), WORK(*)
C
C
      DIMENSION DTRCE(4)
C
C
      CALL QENTER('RMAKDM')
      CALL RTKTIME(.TRUE.,8)
      KFRSAV = KFREE
C
      IF ( OPT_CIPROGRAM .EQ. 'GASCIP' .OR.
     &     OPT_CIPROGRAM .EQ. 'LUCIAREL') THEN
C
C        Calculate density matrices in MOLFDIR format
C
         IF (CALCPV) THEN
            LPV = (2*NASHT)**4 * MIN(2,NZ)
         ELSE
            LPV = 0
         END IF
         CALL MEMGET2('REAL','PVMOLF',KPVMOLF,LPV,WORK,KFREE,LFREE)
chjaaj:  DVMOLF is always calculated in this version
c        IF (CALCDV) THEN
            LDV = (2*NASHT)**2 * MIN(2,NZ)
c        ELSE
c           LDV = 0
c        END IF
         CALL MEMGET2('REAL','DVMOLF',KDVMOLF,LDV,WORK,KFREE,LFREE)
C
         IF ( OPT_CIPROGRAM .EQ. 'GASCIP' ) THEN
C        Density matrices generated by KRMC GASCIP
C           (this routine always calculates DVMOLF)
            CALL GASCIP_RMAKDM(.true.,CALCPV,NZCONF,WORK(KZCONF),
     &           CREF,CREF,WORK(KDVMOLF),WORK(KPVMOLF),IPROPT)
C
         ELSE IF (OPT_CIPROGRAM.eq.'LUCIAREL') then
C           CALL LUCI_DENS(CALCR1,CALCR2,RHO1,RHO2,
C    &                     CREFL,CREFR,WORK,KFREE,LFREE)
            CALL LUCI_DENS(CALCDV,CALCPV,WORK(KDVMOLF),WORK(KPVMOLF),
     &                     CREF,CREF,WORK,KFREE,LFREE)
         END IF
C
C        Transform to X+, x++ format for DV and PV
C
!        IPROPT = 40
         IF (CALCDV) CALL MFC2QFC(.TRUE.,WORK(KDVMOLF),DV,IPROPT)
C                    CALL MFC2QFC(KSYMM,DVM,DVD,IPRINT)
!        IPROPT = 00
         IF (CALCPV) THEN
            CALL M2DNZ3(.FALSE.,.TRUE.,.TRUE.,
     &           WORK(KPVMOLF),PV,IPROPT)
C           CALL M2DNZ3(TSYMM,KSYMM,PSYMM,PVM,PVD,IPRINT)
         END IF
CSK      IPROPT = 00
         CALL MEMREL('RMAKDM',WORK,1,KFRSAV,KFREE,LFREE)
C
      ELSE IF ( NASHT .EQ. 1 .AND. NAELEC .EQ. 1) THEN
C
C        Trivial density matrices for single open shell SCF,
C        done without the overhead of general CI routines.
C
         IF (CALCDV) THEN
            DV(1) = D1
            DO IZ = 2, NZ
               DV(IZ) = D0
            END DO
         END IF
         IF (CALCPV) THEN
            CALL QUIT('CALCPV not implemented for NASHT.eq.1')
         END IF
      ELSE
         CALL QUIT('Unknown OPT_CIPROGRAM: '//OPT_CIPROGRAM)
      END IF
C
C
C        Check trace of DV
C
C        be a bit more slack (?) when running NO occ. numbers
C        with LUCIAREL
C
      IF( OPT_CIPROGRAM.eq. 'LUCIAREL' .and. NATOLCR ) THEN
        XDEPS = 1.0D-6
      END IF
C
C     IF (... debug check of DV,PV trace ...) THEN
         IF (CALCDV) THEN
            DO IZ = 1, NZ
               DTRCE(IZ) = D0
               DO I = 1,NASHT
                  DTRCE(IZ) = DTRCE(IZ) +
     &                 DV( I + (I-1)*NASHT + (IZ-1)*N2ASHX )
               END DO
            END DO
#if defined MCSCF_DEBUG
            CALL HEADER('RMAKDM: trace of DV',-1)
            WRITE(LUPRI,'(A,4X,3A)')
     &      'Matrix','    Actual    ','   Expected   ',
     &      '     Diff     '
            DO IZ = 1, NZ
               IF ( IZ .EQ. 1 ) THEN
                  DTREXP = NAELEC
               ELSE
                  DTREXP = D0
               END IF
               WRITE(LUPRI,'(1P,A1,A,6X,D10.4,4X,D10.4,4X,D10.4)')
     &         QUNIT(IPQTOQ(IZ,0)),'-part',
     &         DTRCE(IZ),DTREXP,ABS(DTRCE(IZ)-DTREXP)
            END DO
#endif
 100        CONTINUE
            IF ( IPROPT .GE. 10 ) THEN
               CALL HEADER('RMAKDM: trace of DV',-1)
               WRITE(LUPRI,'(A,4X,3A)')
     &              'Matrix','    Actual    ','   Expected   ',
     &              '     Diff     '
            END IF
            DO IZ = 1, NZ
               IF ( IZ .EQ. 1 ) THEN
                  DTREXP = NAELEC
               ELSE
                  DTREXP = D0
               END IF
               IF (IPROPT .GE. 10) THEN
                  WRITE(LUPRI,'(1P,A1,A,6X,D10.4,4X,D10.4,4X,D10.4)')
     &                 QUNIT(IPQTOQ(IZ,0)),'-part',
     &                 DTRCE(IZ),DTREXP,ABS(DTRCE(IZ)-DTREXP)
               END IF
               IF ( ABS ( DTRCE(IZ) - DTREXP ) .GE. DEPS ) THEN
                  IF(OPT_CIPROGRAM .eq. 'LUCIAREL' .and. NATOLCR ) THEN
                    IF (ABS(DTRCE(IZ) - DTREXP) .GE. XDEPS)THEN
                      GOTO 150
                    END IF
                    GOTO 200
                  END IF
 150              CONTINUE
                  IF (IPROPT .LT. 10) THEN
                     IPROPT = 10
                     GOTO 100
                  END IF
                  WRITE(LUPRI,'(/A)') ' (RMAKDM) DV density matrix'
                  CALL PRQMAT(DV,NASHT,NASHT,NASHT,NASHT,NZ,
     &                        IPQTOQ(1,0),LUPRI)
                  CALL QUIT('*** ERROR in RMAKDM ***')
 200              CONTINUE
               END IF
            END DO
         END IF
         IF (CALCPV) THEN
C
C           TODO FIXME: check trace of PV
C
         END IF
C     END IF
      IF ( IPROPT .GE. 7 ) THEN
         IF (CALCDV) THEN
            WRITE(LUPRI,'(/A)') ' (RMAKDM) DV density matrix'
            CALL PRQMAT(DV,NASHT,NASHT,NASHT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
         IF (CALCPV) THEN
            WRITE(LUPRI,'(/A)') ' (RMAKDM) PV density matrix'
            CALL PRDNZ3(PV,NASHT,NNASHX,NZ,IPQTOQ(1,0),LUPRI)
         END IF
      END IF
C
      CALL RTKTIME(.FALSE.,8)
      CALL QEXIT('RMAKDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rmakdv */
      SUBROUTINE RMAKDV(DV,CREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Construct one-electron active density matrix
C
C     Input:
C        CREF    - CI vector
C
C     Output:
C        DV      - one electron active density matrix
C
C     Written by J. Thyssen - Dec 21 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION DV(*)
      DIMENSION WORK(*)
C
      CALL QENTER('RMAKDV')
C
      CALL RMAKDM(.TRUE.,.FALSE.,DV,DUMMY,CREF,WORK,KFREE,LFREE)
C
      CALL QEXIT('RMAKDV')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rfckmat */
      SUBROUTINE RFCKMAT(DV,CMO,EMCMY,EMCACT,FC,FV,WORK,LWORK)
C***********************************************************************
C
C     Construct active and inactive Fock matrices
C
C     Input:
C        CMO     - the orbitals for the CEP
C        DV      - the active one electron density matrix
C        EMCACT  - the energy of CEP
C
C     Output:
C
C        EMCMY   - inactive energy
C        FC      - the inactive Fock matrix
C        FV      - the active Fock matrix
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision :
C      Oct 24, 2003: Diagonalization of FC + FV matrix,
C                    save orbital energies in ORBENE      T. Fleig
C
C***********************************************************************
      use dirac_cfg ! for MCSCF-srDFT calculations
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "consts.h"
C
#include "dcbfir.h"
#include "blocks.h"
#include "dcborb.h"
#include "dcbkrmc.h"
#include "dcbopt.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "cbihr2.h"
      DIMENSION DV(*),CMO(*),FC(*),FV(*)
      DIMENSION WORK(LWORK)
      LOGICAL NOFV
      LOGICAL FNDLAB
      DIMENSION IQBF(4)
      DIMENSION IOPSST(2,MXGAS+2,2),NOPSST(2,MXGAS+2,2)
C
C
#include "memint.h"
C
      CALL QENTER('RFCKMAT')
      CALL RTKTIME(.TRUE.,12)
C
      CALL MEMGET2('REAL','H1AO',KH1AO,N2BBASXQ,WORK,KFREE,LFREE)
      IF (NASHT .GT. 0) THEN
         IF (NISHT .GT. 0) THEN
C
            CALL MEMGET2('REAL','FCAO',KFCAO,2*N2BBASXQ,
     &         WORK,KFREE,LFREE)
            KFAO  = KFCAO
            KFVAO = KFCAO + N2BBASXQ
C
            CALL MEMGET2('REAL','DCAO',KDCAO,2*N2BBASXQ,
     &         WORK,KFREE,LFREE)
            KDAO  = KDCAO
            KDVAO = KDCAO + N2BBASXQ
C
            NFMAT = 2
C
         ELSE
C
            CALL MEMGET2('REAL','FCAO',KFCAO,2*N2BBASXQ,
     &         WORK,KFREE,LFREE)
            KFVAO = KFCAO + N2BBASXQ
            KFAO  = KFVAO
C
            CALL MEMGET2('REAL','DVAO',KDVAO,N2BBASXQ,WORK,KFREE,LFREE)
            KFCAO = KDVAO
            KDAO  = KDVAO
C
            NFMAT = 1
         END IF
      ELSE
         IF (NISHT .GT. 0) THEN
C
            CALL MEMGET2('REAL','FCAO',KFCAO,N2BBASXQ,WORK,KFREE,LFREE)
            KFAO = KFCAO
C
            CALL MEMGET2('REAL','DCAO',KDCAO,N2BBASXQ,WORK,KFREE,LFREE)
            KDAO = KDCAO
C
            NFMAT = 1
C
         ELSE
            CALL QUIT('no inactive nor active orbitals...')
         END IF
      END IF

      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. , .TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      if (NPOS.GT.0) THEN
         CALL MEMGET2('INTE','POS',KPOS,NPOS,WORK,KFREE,LFREE)
      else
         KPOS = KFREE
      endif
C
C     ***************************************************
C     *** Get one-electron density matrices *************
C     ***************************************************
C
C     Construct DCAO
C
      IF (NISHT .GT. 0) THEN
         CALL DENSTY(WORK(KDCAO),CMO,IPROPT)
      END IF
C
C     Construct DVAO
C
      IF (NASHT .GT. 0) THEN
         CALL DENSDV(WORK(KDVAO),CMO,DV,IPROPT,WORK(KFREE),LFREE)
C
C        Scale DVao with 0.5 due to scaling with 2.0 in TWOFCK.
C        (TWOFCK assumes inactive density matrices normalized to 1)
C
         CALL DSCAL(N2BBASXQ,DP5,WORK(KDVAO),1)
      END IF
C
C     ***************************************************
C     *** Get 2-electron Fock matrices ******************
C     ***************************************************
C
C
C     Manu: if MCSCF-srDFT then lr coulomb+exchange integrals
C           will be calculated (this will be specified in twofck
C           with dirac_cfg_srint=.FALSE.)
C
      CALL TWOFCK(IOPT_ISYMOP,IOPT_IHRMOP,IOPT_IFCKOP,
     &            WORK(KFAO),WORK(KDAO),
     &            NFMAT,WORK(KPOS),IOPT_INTFLG,IPRTWO,WORK(KFREE),LFREE)
C
CCCCCCCCCCCCCCCCC
C     Manu: some old comments based on DALTON implementation ...
C     idea: call twofck twice:
C                             - coulomb{^lr} + exchange^{lr}
C                             - coulomb^{sr}
C     EJKVSR = tr(DAO * J^{sr} * DVAO) see DALTON for more info)
C     now call xcfun... with the (inactive+active)DAO to obtain V_{xc}^{sr}
C     and E_{xc}^{sr}
C     having it at hand we calculate -tr(V_{xc}^{sr} * DVAO) - 1/2 tr(V_{xc}^{sr} * DCAO)
CCCCCCCCCCCCCCCCC
C
C     ***************************************************
C     *** Get 1-electron Fock matrix ********************
C     ***************************************************
C
      CALL DIRH1(WORK(KH1AO),IPROPT,WORK(KFREE),LFREE)
C
C     ***************************************************
C     *** Calculate inactive energy *********************
C     ***************************************************
C
      IF (NISHT .GT. 0) THEN
C
C        1-electron inactive energy
C        (factor D2: DCao normed to 1 (not 2))
C
C                  tr(h1 * DC)
         EMCMY1E = D2*DDOT(N2BBASXQ,WORK(KDCAO),1,WORK(KH1AO),1)
C
C        Add 2-electron inactive energy
C        (factor D2: DCao normed to 1 (not 2))
C
C                1/2 *  tr(DC * G * DC) where G = J -K ! Note that FCAO does not contain h1 yet ...
         EMCMY2E = DP5 * D2 * DDOT(N2BBASXQ,WORK(KDCAO),1,WORK(KFCAO),1)
         EMCMY = EMCMY1E + EMCMY2E
C
C        Note that for MCSCF-srDFT, EMCMY is the purely long-range
C        inactive energy (the srHxc potential has NOT been added to the
C        inactive Fock operator FCAO yet).
C
C
         IF (IPROPT .GE. 4) THEN
            WRITE(LUPRI,9000) 1,EMCMY1E
            WRITE(LUPRI,9000) 2,EMCMY2E
            WRITE(LUPRI,9001) EMCMY
         END IF
      ELSE
         EMCMY = D0
      END IF
 9000 FORMAT (/'  rFCKMAT: ',I1,'-electron inactive energy: ',F25.15)
 9001 FORMAT (/'  rFCKMAT: Total inactive energy     : ',F25.15)
C
C     ***************************************************
C     *** Add one-electron contribution to FCAO *********
C     ***************************************************
C
      IF (NISHT .GT. 0) THEN
         CALL DAXPY(N2BBASXQ,D1,WORK(KH1AO),1,WORK(KFCAO),1)
      ELSE
         CALL DCOPY(N2BBASXQ,WORK(KH1AO),1,WORK(KFCAO),1)
      END IF

#ifdef DEBUG_SRDFT
C
C     If MCSCF-srDFT calculation: add the srH
C     potential Jsr . D  to the inactive Fock operator
C     FCAO, where D = DC + DV is the total density matrix.
C     Note that, in addition, the srH correction
C     1/2 tr (D * Jsr * D) - tr (DV * Jsr * D)
C     is added to the inactive energy. The second term relates
C     to the active energy contribution tr(FC*DV) which will be
C     calculated later on with the MCSCF-srDFT inactive Fock operator
C     that contains Jsr * D.
C
C     call addsrh ...
C
C
C     If MCSCF-srDFT calculation: add the srxc
C     potential V^{sr}_{xc}[D] to the inactive Fock operator
C     FCAO, where D = DC + DV is the total density matrix.
C     Note that, in addition, the srXC correction
C     E^{sr}_{xc}[D] - tr (V^{sr}_xc[D] * dvao)
C     is added to the inactive energy. The second term relates
C     to the active energy contribution tr(FC*DV) which will be
C     calculated later on with the MCSCF-srDFT inactive Fock operator
C     that contains V^{sr}_xc[D].
      if (dirac_cfg_dft_calculation) then
         call addsrxc(ntbas(0),
     &                nz,
     &                nisht,
     &                work(kfcao),
     &                work(kdcao),
     &                work(kdvao),
     &                emcmy) ! inactive energy
      end if
#endif
C


C
      IF ( IPROPT .GE. 10 ) THEN
         WRITE(LUPRI,9005)
         IF ( NISHT .GT. 0 ) THEN
            WRITE(LUPRI,9006) 'FCao'
            CALL PRQMAT(WORK(KFCAO),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     $           NZ,IPQTOQ(1,0),LUPRI)
         END IF
         IF ( NASHT .GT. 0 ) THEN
            WRITE(LUPRI,9006) 'FVao'
            CALL PRQMAT(WORK(KFVAO),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     $           NZ,IPQTOQ(1,0),LUPRI)
         END IF
      END IF
 9005 FORMAT(/'rFCKMAT: AO Fock matrices')
 9006 FORMAT(/'rFCKMAT: ',A,' fock matrix')
C
C     ***************************************************
C     *** Transform to MO basis *************************
C     ***************************************************
C
      CALL DZERO(FC,N2ORBXQ)
      IF (NASHT .GT. 0) CALL DZERO(FV,N2ORBXQ)
      DO 100 I = 1,NFSYM
         IF (NORB(I) .EQ. 0) GOTO 100
         IF (IPROPT .GE. 10) WRITE(LUPRI,9010) I
         CALL QTRANS('AOMO','S',D0,
     &               NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &               WORK(KFCAO+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,
     &               IPQTOQ(1,0),
     &               FC(I2ORBX(I,I)+1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               WORK(KFREE),LFREE,IPROPT)
         IF (IPROPT .GE. 10) THEN
            WRITE(LUPRI,9011) 'FC',I
            CALL PRQMAT(FC(I2ORBX(I,I)+1),NORB(I),NORB(I),
     &                  NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
         IF (NASHT .EQ. 0) GOTO 100
         CALL QTRANS('AOMO','S',D0,
     &               NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &               WORK(KFVAO+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,
     &               IPQTOQ(1,0),
     &               FV(I2ORBX(I,I)+1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               WORK(KFREE),LFREE,IPROPT)
         IF (IPROPT .GE. 10) THEN
            WRITE(LUPRI,9011) 'FC',I
            CALL PRQMAT(FC(I2ORBX(I,I)+1),NORB(I),NORB(I),
     &                  NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
            WRITE(LUPRI,9011) 'FV',I
            CALL PRQMAT(FV(I2ORBX(I,I)+1),NORB(I),NORB(I),
     &                  NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
  100 CONTINUE
      CALL MEMREL('RFCKMAT',WORK,KWORK,KWORK,KFREE,LFREE)
C
c     call getmat('FCMO',FC,10,.TRUE.,WORK(KFREE),LFREE)
c     call prqmat(fc,norbt,norbt,norbt,norbt,nz,ipqtoq(1,0),lupri)
c     call getmat('FVMO',FV,10,.TRUE.,WORK(KFREE),LFREE)
c     call prqmat(fv,norbt,norbt,norbt,norbt,nz,ipqtoq(1,0),lupri)
C
      CALL RTKTIME(.FALSE.,12)
      CALL QEXIT('RFCKMAT')
      RETURN
 9010 FORMAT(/'rFCKMAT: MO Fock matrices, irrep ',I1)
 9011 FORMAT(/'rFCKMAT: ',A,' fock matrix, irrep ',I1)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgetfc */
      SUBROUTINE rGETFC(TOFILE,CMO,FC,ECORE,WORK,LWORK)
C***********************************************************************
C
C     Construct inactive Fock matrices
C
C     Input:
C        CMO     - the orbitals for the CEP
C
C     Output:
C
C        FC      - the inactive Fock matrix
C
C     Written by J. Thyssen - June 30 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "consts.h"
C
#include "dcbfir.h"
#include "blocks.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "frame.h"
#include "cbihr2.h"
      LOGICAL   TOFILE
      DIMENSION CMO(*),FC(*), WORK(LWORK)
      DIMENSION IQBF(4)
C
C
#include "memint.h"
C
      CALL QENTER('rGETFC')

      IPRINT_local = IPROPT
      IF (X2CMMF) IPRINT_local = max(2,IPRINT_local)
C
C
C     ***************************************************
C     *** Get 1-electron Fock matrix ********************
C     ***************************************************
C
      CALL MEMGET2('REAL','H1AO',KH1AO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL DIRH1(WORK(KH1AO),IPRINT_local,WORK(KFREE),LFREE)
C
C     ***************************************************
C     *** Get 2-electron Fock matrix ********************
C     ***************************************************
C
      CALL MEMGET2('REAL','FCAO',KFCAO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DCAO',KDCAO,N2BBASXQ,WORK,KFREE,LFREE)
      IF (NISHT .GT. 0) THEN

         NFMAT = 1

         call SetTaskDistribFlags((/.TRUE., .TRUE., .TRUE., .TRUE. /))
         call SetIntTaskArrayDimension(NPOS,PARCAL)
         CALL MEMGET2('INTE','POS',KPOS,NPOS,WORK,KFREE,LFREE)

C        Construct one-electron density matrix DCAO

         CALL DENSTY(WORK(KDCAO),CMO,IPRINT_local)

         CALL TWOFCK(IOPT_ISYMOP,IOPT_IHRMOP,IOPT_IFCKOP,
     &        WORK(KFCAO),WORK(KDCAO),
     &        NFMAT,WORK(KPOS),IOPT_INTFLG,IPRTWO,WORK(KFREE),LFREE)
      ELSE
         CALL DZERO(WORK(KFCAO),N2BBASXQ)
         CALL DZERO(WORK(KDCAO),N2BBASXQ)
      END IF

C     if x2c_mmf, add the mmf two-electron Fock matrix

      IF (X2CMMF) THEN
         CALL MEMGET2('REAL','F2MMF',KF2MMF,N2BBASXQ,WORK,KFREE,LFREE)
         open(lux2c,file='DFFCK_mmf',status='old',form='unformatted',
     &        access='sequential',action='readwrite',
     &        position='rewind')
         CALL REAFCK(lux2c,WORK(KF2MMF),.TRUE.,1)
         close(lux2c, status='keep')
         CALL DAXPY(N2BBASXQ,-1.0D0,WORK(KH1AO),1,WORK(KF2MMF),1)
         EMMF2E = DDOT(N2BBASXQ,WORK(KDCAO),1,WORK(KF2MMF),1)  ! see "same as" below
         CALL DAXPY(N2BBASXQ,1.0D0,WORK(KF2MMF),1,WORK(KFCAO),1)
      ELSE
         EMMF2E = 0.0D0
      END IF
C
C        1-electron inactive energy
C        (factor D2: DCao normalized to 1 (not 2))
C
         ECORE = D2*DDOT(N2BBASXQ,WORK(KDCAO),1,WORK(KH1AO),1)
         IF (IPRINT_local .GE. 2) WRITE(LUPRI,9000) 1,ECORE
C
C        Calculate 2-electron inactive energy
C        (factor D2: DCao normed to 1 (not 2))
C
         EMCMY2E =            DDOT(N2BBASXQ,WORK(KDCAO),1,WORK(KFCAO),1)
!same as EMCMY2E = DP5 * D2 * DDOT(N2BBASXQ,WORK(KDCAO),1,WORK(KFCAO),1)
         IF (IPRINT_local .GE. 2) THEN
            WRITE(LUPRI,9000) 2,EMCMY2E - EMMF2E
            IF (X2CMMF) WRITE(LUPRI,9007) EMMF2E
         END IF
         ECORE = ECORE + EMCMY2E
C
C        ***************************************************
C        *** Add two-electron contribution to H1AO *********
C        ***************************************************
C
         CALL DAXPY(N2BBASXQ,D1,WORK(KFCAO),1,WORK(KH1AO),1)
C
C
C     add the LVCORR energy correction to the ECORE energy
C     (will be zero if DOLVC .false.) /Sep 2011 hjaaj
      ECORE    = ECORE + POTNUC + E_LVCORR_MC
      ECORE_LR = ECORE
      IF (IPRINT_local .GE. 2) THEN
          WRITE(LUPRI,9001) POTNUC
          IF (E_LVCORR_MC .ne. 0.0D0) WRITE(LUPRI,9002) E_LVCORR_MC
          WRITE(LUPRI,9004) ECORE
      END IF
C
      IF ( IPRINT_local .GE. 10 ) THEN
         WRITE(LUPRI,9006) 'FCao'
         CALL PRQMAT(WORK(KH1AO),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
 9000 FORMAT(/' (rGETFC)',I2,'-electron core energy: ',F25.15)
 9001 FORMAT(/' (rGETFC) Nuclear repulsion     : ',F25.15)
 9002 FORMAT(/' (rGETFC) LVCORR correction     : ',F25.15)
 9004 FORMAT(/' (rGETFC) Total core energy     : ',F25.15)
 9006 FORMAT(/'  rGETFC: ',A,' fock matrix')
 9007 FORMAT(/' (rGETFC) 2-electron mmf energy : ',F25.15)
C
C     ***************************************************
C     *** Transform to MO basis *************************
C     ***************************************************
C
      CALL DZERO(FC,N2ORBXQ)
      DO 100 I = 1,NFSYM
         IF (NORB(I) .EQ. 0) GOTO 100
         CALL QTRANS('AOMO','S',D0,
     &               NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &               WORK(KH1AO+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,
     &               IPQTOQ(1,0),
     &               FC(I2ORBX(I,I)+1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               WORK(KFREE),LFREE,IPRINT_local)
         IF (IPRINT_local .GE. 10) THEN
            WRITE(LUPRI,9011) 'FC',I
            CALL PRQMAT(FC(I2ORBX(I,I)+1),NORB(I),NORB(I),
     &                  NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
  100 CONTINUE
C
C     Write to file KRMC FOCK
C
      IF (TOFILE) THEN
         CALL WRTKRMC(LUKRM3,'FCMO    ',FC,N2ORBXQ)
         IF( CWRTFO_MAT ) CALL WRTKRMC(LUKRM3,'ECORE_LR',ECORE_LR,1)
      END IF
C
      CALL MEMREL('rGETFC',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL FLSHFO(LUPRI)
      CALL QEXIT('rGETFC')
      RETURN
 9011 FORMAT(/' rGETFC: ',A,' fock matrix, irrep ',I1)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirh1 */
      SUBROUTINE DIRH1(H1AO,IPRINT,WORK,LWORK)
C***********************************************************************
C
C     Get one-electron Fock matrix in AO-basis
C
C     Input:
C
C     Output:
C        H1AO    - one-electron Fock matrix
C
C     Written by J. Thyssen - Oct 26 1998
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "cbihr1.h"
C
      DIMENSION H1AO(*), WORK(*)
      LOGICAL   EX
C
      CALL QENTER('DIRH1')
C
C     (1) If file DFFCK1 is present, then read H1AO from DFFCK1.
C     ----------------------------------------------------------
C
      INQUIRE (FILE='DFFCK1', EXIST=EX)
      IF (EX) THEN
         CALL OPNFIL(LUFCK1,'DFFCK1','OLD','DIRH1')
         CALL REAFCK(LUFCK1,H1AO,.TRUE.,1)
         CLOSE(LUFCK1)
      ELSE
C
C     (2) Generate H1AO and save on DFFCK1.
C     -------------------------------------
C
         CALL ONEFCK(H1AO,IPRONE,WORK,LWORK)
         CALL OPNFIL(LUFCK1,'DFFCK1','UNKNOWN','DIRH1')
         CALL WRIFCK(LUFCK1,H1AO,1)
         CLOSE (LUFCK1,STATUS='KEEP')
C
      END IF
      IF ( IPRINT .GE. 10 ) THEN
         WRITE(LUPRI,'(/A,L5)') ' (DIRH1) H1 in AO basis', EX
         CALL PRQMAT(H1AO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     $        NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      CALL QEXIT('DIRH1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsetwop */
      SUBROUTINE RSETWOP(WFTYP,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Initialize configurational and orbital parameters for
C     (KR-MC)SCF calculation.
C
C     Input:
C        OPT_SKIPEE  - skip e-e rotations
C        OPT_SKIPEP  - skip e-e rotations
C        OPT_CHCKJZ  - check rotation vector wrt m_j value of the orbs
C
C     Output:
C        KZCONF   - work address of info about CI determinants
C        LZCONF   - integer length of WORK(KZCONF)
C
C        NZXOPE   - number of ++ orb rotations (non-redundant)
C        NZHOPE   - number of ++ orb rotations  incl. redundant active-active rotations
C        LZXOPE   - integer length of ++ orb rotations
C
C        NZXOPP   - number of +- orb rotations
C        LZXOPP   - integer length of +- orb rotations
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision : July 2002 hjaaj
C                     Jan  2010 sk
C
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
      CHARACTER WFTYP*(*)
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbos.h"
#include "dcbnrt.h"
#include "dcbidx.h"
C
      DIMENSION WORK(*), ILIST(MXCORB)
      LOGICAL   GASWF

      integer, allocatable :: reduced_pp(:)
      integer, allocatable :: reduced_pn(:)

C
      CALL QENTER('RSETWOP')
C
C     ***************************************************
C     *** Approximate boson symmetry ? ******************
C     ***************************************************
C
C     If the use of boson or ``approximate'' boson symmetry
C     is requested calculate these.
C
C     This is needed for the current version of LUCIAREL, but
C     can also be used for spin-free and Levy-Leblond Hamiltonians.
C     (Only if LUCIAREL is run in spinfree formalism. Otherwise
C     boson symmetry array is reset.)
C
C     this step is also done for GASCIP.
C
      IF (OPT_UCIBOS) THEN
C
C        For the current version of LUCIAREL we need the
C        approximate boson symmetries for labelling of
C        orbitals. Call ROPTST to retrieve the orbitals.
C
C        TODO FIXME: the control parameter should be consistent
C        with IOPTST.
C
         CALL MEMGET2('REAL','CMO',KCMO,NCMOTQ,WORK,KFREE,LFREE)
         CALL MEMGET2('INTE','IBEIG',KIBEIG,NORBT,WORK,KFREE,LFREE)
         ISTART = MOD(IOPTST,100)
         CALL ROPTST(ISTART,0,WORK(KCMO),DUMMY,WORK(KIBEIG),-1,
     &               WORK,KFREE,LFREE)
C        CALL ROPTST(ISTART,IWRT,CMO,CREF,IBEIG,IPRINT,WORK,KFREE,LFREE)
         CALL RAPPBOS(WORK(KCMO),IBOSYM,WORK(KIBEIG),IPROPT,
     &                WORK(KFREE),LFREE)
         CALL MEMREL('RSETWOP.LUCIAREL',WORK,1,KCMO,KFREE,LFREE)
C
      END IF
C
C     ***************************************************
C     *** Configurational parameters ********************
C     ***************************************************
C
      CALL RSETCI(OPT_CIPROGRAM,WORK,KFREE,LFREE)
      IF (WFTYP(1:3) .NE. 'SCF' .AND. NZCONF .LE. 0) THEN
         GO TO 9000
      END IF
C
C     ***************************************************
C     *** Orbital parameters ****************************
C     ***************************************************
C
      LNOROT = .FALSE.
      NOROT(1:NORBT) = 0
C
      NDSHT = 0
      NVSHT = 0
C
      WRITE(LUPRI,'(/A)') ' (RSETWOP) Orbital classes:'
      DO I = 1,NFSYM
         IF (NFSYM.GT.1) WRITE(LUPRI,'(3X,A,I1)') 'Fermion irrep ',I
         IF (NISH(I).GT.0) THEN
            WRITE(OPT_INDSTR(1,I),'(I4,A2,I4)')
     &         1,'..',NISH(I)
         ELSE
            OPT_INDSTR(1,I) = ' '
         END IF
         IF (NASH(I).GT.0) THEN
            WRITE(OPT_INDSTR(2,I),'(I4,A2,I4)')
     &         NISH(I)+1,'..',NISH(I)+NASH(I)
         ELSE
            OPT_INDSTR(2,I) = ' '
         END IF
         IF (NPSH(I) .GT. 0 .AND. NSSH(I) .GT. 0) THEN
            WRITE(OPT_INDSTR(3,I),'(I4,A2,I4,A1,I4,A2,I4)')
     &         -NPSH(I),'..',-1,',',NISH(I)+NASH(I)+1,'..',NESH(I)
         ELSE IF (NPSH(I) .GT. 0) THEN
            WRITE(OPT_INDSTR(3,I),'(I4,A2,I4)') -NPSH(I),'..',-1
         ELSE IF (NSSH(I) .GT. 0) THEN
            WRITE(OPT_INDSTR(3,I),'(I4,A2,I4)')
     &         NISH(I)+NASH(I)+1,'..',NESH(I)
         ELSE
            OPT_INDSTR(3,I) = ' '
         END IF
C
         WRITE(LUPRI,'(5X,A,A50)')
     &      'Inactive orbitals   : ',OPT_INDSTR(1,I)
         WRITE(LUPRI,'(5X,A,A50)')
     &      'Active orbitals     : ',OPT_INDSTR(2,I)
         WRITE(LUPRI,'(5X,A,A50)')
     &      'Secondary orbitals  : ',OPT_INDSTR(3,I)
C
C        ***************************************************
C        Special handling of deleted secondary orbitals here.
C        (Information about such deleted orbitals were not
C         available before now, therefore NDSH*, NVSH* and
C         the associated IDX* arrays could not be made in
C         the SETDCBIDX routines.) SK+HJAaJ Aug 08
C        ***************************************************
C
         IF (NFRO(I) .GT. 0) THEN
            WRITE(LUPRI,'(5X,A,I4)')
     &      'Frozen orbitals     : 1..',NFRO(I)
            LNOROT = .TRUE.
            DO J = 1, NFRO(I)
               K = IORB(I) + NPSH(I) + J
               NOROT(K) = 1
            END DO
         END IF
         NVEC = 1
         CALL NUMLST(OPT_DELSTR(I),ILIST,NORB(I),
     &        -NPSH(I),NESH(I),I,NVEC)
         NDSH(I) = 0
C
C        changed to OPT_DELSTR which is more natural and also
C        in accordance with SIRIUS-DALTON code - SK Aug 2008
C
C        we still want to freeze but doubly occupied core
C        orbitals...
C        TODO FIXME: IMPLEMENT ME - SK Aug 2008
C
C        hjaaj comment Feb 2017: we could use OPT_DELSTR
C        also for .FREEZE and remove .DELETE
C        (i.e. let .DELETE be part of .FREEZE);
C        I have today implemented that NDSH(I) only
C        counts deleted secondary orbitals.
C        
C
csk      CALL NUMLST(OPT_FRZSTR(I),ILIST,NORB(I),
csk  &        -NPSH(I),NESH(I),I,NVEC)
         IF (NVEC .GT. 0) THEN
            LNOROT = .TRUE.
csk         WRITE(LUPRI,'(A,A50)')
csk  &           '     Frozen orbitals      : ',OPT_FRZSTR(I)
            WRITE(LUPRI,'(A,A50)')
     &           '     Deleted orbitals     : ',OPT_DELSTR(I)
            NVEC_I = 0
            DO J = 1, NVEC
               K = ILIST(J) + IORB(I) + NPSH(I)
               NOROT(K) = 1
               IF (ILIST(J) .GT. NOCC(I)) NVEC_I = NVEC_I + 1
            END DO
            NDSH(I) = NVEC_I
         END IF
C
         NVSH(I) = NSSH(I) - NDSH(I)
C
         NDSHT = NDSHT + NDSH(I)
         NVSHT = NVSHT + NVSH(I)
      END DO
      NOROTC(1:NORBT) = NOROT(1:NORBT)
C     ... 'active' secondary orbitals
      N2VSHX   = NVSHT   * NVSHT
      N2VSHXQ  = N2VSHX  * NZ
C     ... deleted secondary orbitals
      N2DSHX  = NDSHT  * NDSHT
      N2DSHXQ = N2DSHX * NZ
C     ... mapping for 'active' and deleted secondary orbitals
C     IDXG2D - general to deleted
C        D2G - deleted to general
C        V2G - secondary to general
C        G2V - general to secondary
      IORBIX = 0
      IDSHI  = 0
      IORBDX = 0
      IVSHII = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NPSH(IFSYM)
            IORBIX = IORBIX + 1
            IDXG2D(IORBIX) = -1
            IDXD2G(IORBIX) = -1
            IDXV2G(IORBIX) = -1
            IDXG2V(IORBIX) = -1
         END DO
         DO I = 1, NISH(IFSYM)
            IORBIX = IORBIX + 1
            IDXG2D(IORBIX) = -1
            IDXD2G(IORBIX) = -1
            IDXV2G(IORBIX) = -1
            IDXG2V(IORBIX) = -1
         END DO
         DO I = 1, NASH(IFSYM)
            IORBIX = IORBIX + 1
            IDXG2D(IORBIX) = -1
            IDXD2G(IORBIX) = -1
            IDXV2G(IORBIX) = -1
            IDXG2V(IORBIX) = -1
         END DO
         DO I = 1, NVSH(IFSYM)
            IORBIX = IORBIX + 1
            IVSHII = IVSHII + 1
C
            IDXG2D(IORBIX) = -1
            IDXG2V(IORBIX) = IVSHII
C
            IDXD2G(IORBIX) = -1
            IDXV2G(IVSHII) = IORBIX
         END DO
         DO I = 1, NDSH(IFSYM)
            IORBIX = IORBIX + 1
            IDSHI  = IDSHI  + 1
            IDXG2D(IORBIX) = IDSHI
            IDXG2V(IORBIX) = -1
C
            IDXD2G(IDSHI)  = IORBIX
         END DO
      END DO
      IF (IPROPT .ge. 4 .AND. NDSHT .gt. 0) THEN
         WRITE(LUPRI,'(/A)')
     &        ' (RSETWOP) V and IDXV2G'
         DO I = 1, NVSHT
           WRITE(LUPRI,'(8X,I6,I11)') I, IDXV2G(I)
         END DO
         WRITE(LUPRI,'(/A)')
     &        ' (RSETWOP) D and IDXD2G'
         DO I = 1, NDSHT
           WRITE(LUPRI,'(8X,I6,I11)') I, IDXD2G(I)
         END DO
      END IF
C
      CALL FLSHFO(LUPRI)
C

      IF (WFTYP(1:5) .NE. 'KR-CI') THEN     ! no orbital rotation vectors for CI
         
C
C     Always totally symmetric rotations (=> JOPSY = 1)
C
      JOPSY = 1
      GASWF = NGAS_DC.GT.1
      CALL SETXOP(OPT_INDSTR,OPT_SKIPEE,OPT_SKIPEP,GASWF,JOPSY,NOROT,
     &     NZXOPE,NZHOPE,NZXOPP,
     &     IPROPT)
C
C
C     MCSCF calculation in linear symmetry: remove potentially
C     symmetry-breaking orbital rotations - this solves
C     the problem in MCSCF calculations (observed for U2)
C     with active spaces that contain orbitals of different mj-values.
C
C     similar arguments apply for spinfree MCSCF: remove rotations
C     between orbitals of different boson symmetry.
C
      IF((OPT_CHCKJZ.AND.LINEAR).or. spinfr_krmc)THEN
C
C       set internal control variable
        IREDJOB = 10
C
        if (nzhope > 0) then
          allocate(reduced_pp(nzhope*2))
          call icopy(nzhope*2,
     &               get_orbital_rotation_indices_pp(), 1,
     &               reduced_pp, 1)
        else
          allocate(reduced_pp(1))
        end if

        if (nzxopp > 0) then
          allocate(reduced_pn(nzxopp*2))
          call icopy(nzxopp*2,
     &               get_orbital_rotation_indices_pn(), 1,
     &               reduced_pn, 1)
        else
          allocate(reduced_pn(1))
        end if

        CALL RESETXOP_sym(DUMMY,
     &                    reduced_pp,
     &                    reduced_pn,
     &                    IREDJOB,1)

        call update_orbital_rotation_indices(nzhope,
     &                                       nzxope,
     &                                       nzxopp,
     &                                       reduced_pp,
     &                                       reduced_pn)

         deallocate(reduced_pp, reduced_pn)
C
      END IF

      END IF ! (WFTYP(1:5) .NE. 'KR-CI')
C
C     Define LZXOPE to integer length of ++ orb rotations
C     Define LZXOPP to integer length of +- orb rotations
C
      LZXOPE = 2*NZHOPE
      LZXOPP = 2*NZXOPP
C
C     From now on NOROT is used for the orbitals that should not
C     be included in the 4-index transformation.
C     Reset active indices to 0:
C
      DO I = 1, NFSYM
         DO J = 1, NASH(I)
            NOROT(IORB(I) + NPSH(I) + NISH(I) + J) = 0
         END DO
      END DO
      NOROTC(1:NORBT) = NOROT(1:NORBT)
      IF (IPROPT .GT. 4) THEN
         WRITE(LUPRI,*) 'Orbitals excluded from 4-index transformation:'
         CALL IWRTMA(NOROT,NORBT,1,NORBT,1)
      END IF
C
C     Evaluate all derived NZX* and NZH* info in dcbopt.h;
C     write information to LUKRMC
C
      CALL RSET_KRMC(WORK)
C
 9000 CALL FLSHFO(LUPRI)
      CALL QEXIT('RSETWOP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rset_krmc */
      SUBROUTINE RSET_KRMC(WORK)
C***********************************************************************
C
C     Evaluate all derived NZX* and NZH* info in dcbopt.h;
C     write information to LUKRMC
C
C     Extracted from RSETWOP 9-Jul-2002 hjaaj
C
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbopt.h"
#include "dcborb.h"
C
      DIMENSION WORK(*)
      DIMENSION ITMP(9)
C
      LOGICAL FNDLAB
C
      CALL QENTER('RSET_KRMC')

      NZXOPEQ  = NZXOPE *NZ
      NZXOPPQ  = NZXOPP *NZ
      NZXOPT   = NZXOPE + NZXOPP
      NZXOPTQ  = NZXOPT *NZ
C
      NZHOPEQ  = NZHOPE *NZ
      NZHOPT   = NZHOPE + NZXOPP
      NZHOPTQ  = NZHOPT *NZ
C
C     ***************************************************
C     *** Total number of parameters ********************
C     ***************************************************
C
      NZVAR  = NZCONF  + NZXOPT
      NZVARQ = NZCONFQ + NZXOPTQ
C
C     Write rotations to file
C
      ITMP(1) = NZCONF
      ITMP(2) = NZCONFQ
      ITMP(3) = NZXOPE
      ITMP(4) = NZHOPE
      ITMP(5) = NZXOPP
      ITMP(6) = NZ
      ITMP(7) = LZCONF
      ITMP(8) = LZXOPE
      ITMP(9) = LZXOPP
      CALL IWRTKRMC(LUKRMC,'VAR INFO',ITMP,9)
      IF (LZCONF.GT.0)
     &   CALL IWRTKRMC(LUKRMC,'DET INFO',WORK(KZCONF),LZCONF)
      IF (LZXOPE.GT.0)
     &   CALL IWRTKRMC(LUKRMC,'E-E ROT ',
     &                 get_orbital_rotation_indices_pp(),LZXOPE)
      IF (LZXOPP.GT.0)
     &   CALL IWRTKRMC(LUKRMC,'E-P ROT ',
     &                 get_orbital_rotation_indices_pn(),LZXOPP)
C
C
C     Write DCBORB parameters to KRMCSCF
C     ----------------------------------
C
      CALL LABKRMC(LUKRMC,'DCBORB  ')
      WRITE(LUKRMC)
     &     NFSYM,NGAS_DC,
     &     (NISH(I),I=1,NFSYM),
     &     (NASH(I),I=1,NFSYM),
     &     ((NGSH(I,J),I=1,NFSYM),J=1,NGAS_DC),
     &     (NOCC(I),I=1,NFSYM),
     &     (NSSH(I),I=1,NFSYM),
     &     (NFRO(I),I=1,NFSYM),
     &     NAELEC
C
C     End of reference information, write 'Start Of Data' label
C
      CALL NEWLAB('SODLABEL',LUKRMC,LUPRI)
C
C     Consistency check:
C     ------------------
C
C     Check that content on KRMCOLD = KRMCSCF.
C
C     TODO FIXME: write me!
C     TODO FIXME: also write INTFLG
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('RSET_KRMC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsetci */
      SUBROUTINE RSETCI(USE_CIPROGRAM,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Initialize configurational parameters for
C     (KR-MC)SCF calculation and ESR-CI and ???.
C
C     Input:
C
C     Output:
C
C     Extracted from RSETWOP 9-Jul-2002 by HJAaJ.
C
C***********************************************************************
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
C Used from common blocks:
C  dcbgen: LUKRMC
C  dgroup: NZ
C  dcbopt: NZCONF,NZCONFQ,KZCONF,IOPT_SYMMETRY,LZCONF
C  dcborb: NAELEC,NASHT
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbopt.h"
#include "dcborb.h"
C
      CHARACTER*(*) USE_CIPROGRAM
      DIMENSION WORK(*)
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
      LOGICAL FNDLAB, FND
      character orbital_info_label*8
      integer, allocatable :: orb_sym_vec(:)

      CALL QENTER('RSETCI')
      CALL GETTIM(CPU1,WLL1)


      if(linear.or.spinfr_krmc)then

!       may be oversized if not all orbitals are active
        call alloc(orb_sym_vec,norbt)
        orb_sym_vec = 0

!       get orbital symmetry values - either from a previous MCSCF or HF run
        IF ( LUKRM5 .GT. 0 ) THEN
          REWIND(LUKRM5)
          IF ( FNDLAB('MJVEC   ',LUKRM5) )THEN
            CALL IREAKRMC(LUKRM5,'MJVEC   ',orb_sym_vec,NORBT)
            if(IPROPT .ge. 2)then
              WRITE(LUPRI,'(/A)')
     &   ' (RSETCI) *** orbital mj-vector read from KRMCOLD ***'
            end if
         ELSE IF(FNDLAB('IBEIG   ',LUKRM5))THEN
            CALL IREAKRMC(LUKRM5,'IBEIG   ',orb_sym_vec,NORBT)
            if(IPROPT .ge. 2)then
              WRITE(LUPRI,'(/A)')
     &        ' (RSETCI) *** orbital boson irreps read'//
     &        ' from KRMCOLD ***'
            end if
          END IF
        ELSE
C         read orbital symmetry values from DFCOEF
          CALL OPNFIL(LUCOEF,'DFCOEF','OLD','RSETCI')
          CALL REACMO(LUCOEF,'DFCOEF',DUMMY,DUMMY,orb_sym_vec,TOTERG,8)
          CLOSE(LUCOEF,STATUS='KEEP')
          if(IPROPT .ge. 2)then
            WRITE(LUPRI,'(/A)')
     &      ' (RSETCI) *** orbital symmetry information read'//
     &      ' from DFCOEF ***'
          end if
        END IF
C
C       clean LUKRMC and save orbital symmetry information on file KRMCSCF
        REWIND(LUKRMC)
        FND = FNDLAB('SODLABEL',LUKRMC)
        BACKSPACE LUKRMC
        CALL NEWLAB('SODLABEL',LUKRMC,LUPRI)
        REWIND(LUKRMC)

        if(.not.spinfr_krmc)then
          orbital_info_label = 'MJVEC   '
        else
          orbital_info_label = 'IBEIG   '
        end if

!       save symmetry information
        call iwrtkrmc(lukrmc,orbital_info_label,orb_sym_vec,norbt)

        call dealloc(orb_sym_vec)
      end if ! linear or spinfr_krmc switch

C
C     Check that NAELEC .le. 2 *NASHT
C     (NAELEC: # of active electrons)
C
      IF (NAELEC .GT. 2 * NASHT) THEN
         WRITE(LUPRI,'(//2A,I3,2A,I3,A)')
     &        ' *** ERROR in RSETCI *** ',
     &        'Number of electrons (',NAELEC,') is larger than ',
     &        'the number of spinors (',2 * NASHT,')!'
         CALL QUIT('*** ERROR in RSETCI *** NAELEC > 2 * NASHT')
      END IF
C
C     ***************************************************
C     *** Configurational parameters ********************
C     ***************************************************
C
C
      IF ( NASHT .GT. 0 ) THEN
Chj   IF ( NASHT .GT. 1 ) THEN
C
C        A non-trivial number of CI rotations.
C        hjaaj July 2002: also include NASHT .eq. 1
C                         (giving 1-2 determinants)
C         for debug
         if(fcidump .and. USE_CIPROGRAM .NE. 'GASCIP' ) THEN
           call quit('fcidump option only for GASCIP available')
         end if

         if(fcidump)then
           NZCONF = 1
           LZCONF = 0
           CALL MEMGET2('INT8','ZCONF',KZCONF,LZCONF,WORK,KFREE,LFREE)
           GOTO 999
         end if
C
         IF ( USE_CIPROGRAM .EQ. 'GASCIP' ) THEN
C
C           Generate determinants
C
            KFRSAV = KFREE
C           HJAaJ : re size of MAXSTR ....
C               .... size KISTR ca. sqrt size KIDET ??? No, not when GAS
            MAXSTR = 10 000 000
            MAXDET = MIN(100 000 000, (LFREE-MAXSTR-10 000)/2 )
                         !max 100,000,000 determinants in this version
            CALL MEMGET2('INT8','IDET',KIDET ,2*MAXDET,WORK,KFREE,LFREE)
            CALL MEMGET2('INT8','ISTR',KISTR,   MAXSTR,WORK,KFREE,LFREE)

            CALL GASCIP_GENDET(NZCONF,WORK(KIDET),WORK(KISTR),
     &                        MAXDET,MAXSTR,OPT_UCIBOS,
     &                        IOPT_SYMMETRY,IPROPT)
            CALL MEMREL('RSETCI.GASCIP',WORK,1,KFRSAV,KFREE,LFREE)
C
C           Note: the above memory release should make sure that KZCONF .eq. KIDET,
C           where determinant information is already saved inside GASCIP_GENDET,
C           KZCONF is saved in dcbopt.h and is used in other GASCIP routines
C           to obtain determinant information.
C           (LZCONF = 2 * NCONF because "INTEGER*8 IDET(NZCONF,2)")
C
            LZCONF = 2 * NZCONF
            CALL MEMGET2('INT8','ZCONF',KZCONF,LZCONF,WORK,KFREE,LFREE)
            IF (KZCONF .NE. KIDET) CALL QUIT(
     &         'KZCONF .ne. KIDET -- indicates programming error!')
C
#ifdef MOD_DEBUG
C
C          Generate list of Kramers partners...
C          This is only used for debugging purposes.
C
           IF ( MOD(NAELEC,2) .EQ. 0) THEN
             CALL MEMGET2('INTE','KRLIST',KRLIST,NZCONF,
     &          WORK,KFREE,LFREE)
             CALL GASCIP_KRLIST(NZCONF,WORK(KZCONF),WORK(KRLIST),IPROPT)
             CALL IWRTKRMC(LUKRMC,'KRLIST  ',WORK(KRLIST),NZCONF)
             CALL MEMREL('RSETCI.KRLIST',WORK,1,KRLIST,KFREE,LFREE)
           END IF
#endif
         ELSE IF ( USE_CIPROGRAM .EQ. 'LUCIAREL' ) THEN
C
C           Call driver routine for getting NZCONF
            CALL LUCI_NDET(NZCONF,IOPT_SYMMETRY,WORK(KFREE),LFREE)
C
            LZCONF = 0
            CALL MEMGET2('INTE','ZCONF',KZCONF,LZCONF,WORK,KFREE,LFREE)
C
         ELSE IF ( USE_CIPROGRAM .EQ. 'KRCC' ) THEN
            WRITE(LUPRI,*) 'KRCC: IT is ok'
            NZCONF = 1
            LZCONF = 0
            CALL MEMGET2('INTE','ZCONF',KZCONF,LZCONF,WORK,KFREE,LFREE)
         ELSE
C
            WRITE(LUPRI,'(A/3A)')
     $           '*** ERROR in RSETCI ***',
     $           'Unknown CI program: "',USE_CIPROGRAM,'"'
            CALL QUIT('*** ERROR in RSETCI: Unknown CI program ***')
C
         END IF
         CALL FLSHFO(LUPRI)
      ELSE
C
C        Dirac-Hartree-Fock: only 1 CI coefficient
C
         NZCONF = 0
C
      END IF
C
 999  MZ = MIN(NZ,2)
      NZCONFQ = MZ*NZCONF
C
      CALL GETTIM(CPU2,WLL2)
      CPUTID  = SECTID(CPU2-CPU1)
      WALLTID = SECTID(WLL2-WLL1)
C
      WRITE(LUPRI,'(/A,I12,A,I12,A//5A)')
     &      ' (RSETCI) Number of determinants: ',NZCONF,
     &      ' (dimension ',NZCONFQ,')',
     &      ' CPU (Wall) time for generation of determinants: ',
     &      CPUTID,'(',WALLTID,')'
      IF (NASHT .GT. 0 .AND. NZCONF .LE. 0) THEN
         WRITE(LUPRI,'(/A)') ' INFO: no determinants in this symmetry!'
      END IF
C
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('RSETCI')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RESETXOP_sym(UKAPPA,JXOPE,JXOPP,IREDJOB,IPRINT)
C***********************************************************************
      use memory_allocator
C***********************************************************************
C
C     driver routine for a reset orbital rotation parameters wrt
C
C     - mj-values     --> linear spin-orbit calculation
C     - ml-values     --> linear spinfree   calculation
C     - boson irreps  -->        spinfree   calculation
C
C     Input: (stored on common block in dcbopt.h)
C        NZXOPE  - total number of e-e orbital rotations
C        NZHOPE  - total number of e-e orbital rotations (including
C        redundant rotations)
C        NZXOPP  - total number of e-p orbital rotations
C
C        JXOPE   - e-e orbital rotation vector
C        JXOPP   - e-p orbital rotation vector
C
C
C        strategies:
C        IREDJOB == 9:  delete orbital rotation parameters in the setup.
C        IREDJOB == 10: "purify" kappa matrix in each MACRo iteration.
C
C     Output:
C        updated values of NZXOPE, NZHOPE and NZXOPP
C        reduced lists in JXOPE and JXOPP
C        OR
C        purified kappa matrix
C
C     Written by S. Knecht  - Feb 2010
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dummy.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbopt.h"

C
      DIMENSION UKAPPA(NORBT,NORBT,NZ), JXOPE(2,*), JXOPP(2,*)
      LOGICAL FNDLAB, FND
      character orbital_info_label*8
      integer, allocatable :: orb_sym_vec(:)
      real*8,  allocatable :: ukappa2(:)
C
      CALL QENTER('RESETXOP_sym')
C
      IF(IREDJOB .ge. 9)THEN
C
C...    may be oversized if not all orbitals are active
        call alloc(orb_sym_vec,norbt)
        orb_sym_vec = 0
C
C       get orbital symmetry values - either from a previous MCSCF or HF run
        REWIND(LUKRMC)
        if(.not.spinfr_krmc)then
          orbital_info_label = 'MJVEC   '
        else
          orbital_info_label = 'IBEIG   '
        end if
        IF(FNDLAB(orbital_info_label,LUKRMC))THEN
          CALL IREAKRMC(LUKRMC,orbital_info_label,orb_sym_vec,NORBT)
        ELSE
          WRITE(LUPRI,'(/A)')
     &    ' (RESETXOP_sym) - error: *** orbital symmetry information'//
     &    ' not found on file KRMCSCF. ***'
          call quit('*** error in RESETXOP_sym: orbital symmetry'//
     &              ' information not found on file KRMCSCF.***')
        END IF
C
C       debug print
        IF( IPRINT .ge. 7 )THEN
          WRITE(LUPRI,*) ' DEBUG info: symmetry-info vector of orbitals'
          CALL IWRTMAMN(orb_sym_vec,1,NORBT,1,NORBT,LUPRI)
        END IF
C
        IF( IREDJOB .eq. 10 )THEN
C
C         reduce orbital parameters wrt mj-blocking
          CALL DELPAR(orb_sym_vec,JXOPE,NZXOPE,NZHOPE,JXOPP,NZXOPP,
     &                IPRINT)
        END IF
       call dealloc(orb_sym_vec)
      END IF
C
      CALL QEXIT('RESETXOP_sym')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgetac */
      SUBROUTINE RGETAC(FC,FCAC,IPRINT)
C***********************************************************************
C
C     Extract block with active-active orbital indices out of FC
C
C     Input:
C        FC      - Inactive Fock matrix
C
C     Output:
C        FCAC    - Inactive Fock matrix with active-active indices
C
C     Written by J. Thyssen - Oct 29 1998
C     Last revision :
C                S. Knecht  - Aug 2008
C                simplified to use of MATGAT.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dgroup.h"
#include "maxorb.h"
#include "maxash.h"
#include "dcbidx.h"
C
      DIMENSION FC(*),FCAC(*)
C
      CALL QENTER('RGETAC')
      IF (IPRINT .GE. 7) CALL HEADER('Output from rGETAC',-1)
C
      IF (NASHT .EQ. 0) GOTO 9999
C
C     Print FC
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('rGETAC: FC matrix',-1)
         CALL PRQMAT(FC,NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      CALL DZERO(FCAC,N2ASHXQ)
C
C     extract a-a part of FC
C
      CALL MATGAT(FC,NORBT,NORBT,FCAC,NASHT,NASHT,IDXU2G,NZ)
C
C     Print FCAC
C
      IF (IPRINT .GE. 7) THEN
         CALL HEADER('rGETAC: FCAC matrix',-1)
         CALL PRQMAT(FCAC,NASHT,NASHT,NASHT,NASHT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
 9999 CALL QEXIT('RGETAC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck renract */
      FUNCTION RENRACT(DV,FCAC,FQ,IPRINT)
C***********************************************************************
C
C     Calculate active electronic energy of the KR-MCSCF state
C
C     Input:
C        DV      - one-electron active density matrix
C        FQ      - the Q Fock matrix
C        FCAC    - the inactive Fock matrix with active indices
C
C     Returns:
C        RENRACT - active electronic energy
C
C     Written by J. Thyssen - Oct 29 1998
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
#include "maxorb.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbidx.h"
C
      DIMENSION DV(*),FCAC(*),FQ(NORBT,NASHT,NZ)
C
      CALL QENTER('RENRACT')
C
C     Eactive = \sum_{uv} DV_{uv} FC_{uv} +
C               \sum_{uvxy} P_{uv,xy} (uv|xy)
C             = \sum_{uv} DV_{uv} FC_{uv} +
C               1/2 \sum_{u} FQ_{uu}
C
C
C
C     Print FCAC
C
      IF (IPRINT .GE. 8) THEN
         write(lupri,*) 'n2asht',N2ASHT
         CALL HEADER('rEnract: FCAC matrix',-1)
         CALL PRQMAT(FCAC,NASHT,NASHT,NASHT,NASHT,
     &               NZ,IPQTOQ(1,0),LUPRI)
         CALL HEADER('rEnract: DV matrix',-1)
         CALL PRQMAT(DV,NASHT,NASHT,NASHT,NASHT,
     &               NZ,IPQTOQ(1,0),LUPRI)
         CALL HEADER('rEnract: FQ matrix',-1)
         CALL PRQMAT(FQ,NORBT,NASHT,NORBT,NASHT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C     \sum_{uv} DV_{uv} FC_{uv}
C     -------------------------
C
      EACT = DDOT(N2ASHXQ,DV,1,FCAC,1)
      IF (IPRINT .GE. 4)
     &   WRITE(LUPRI,9000) EACT
C
C     1/2 \sum_{u} FQ_{uu}
C     --------------------
C
      EACTFQ = D0
      DO IV = 1, NASHT
         IGV = IDXU2G(IV)
         EACTFQ = EACTFQ + FQ(IGV,IV,1)
      END DO
      EACTFQ = EACTFQ * DP5
      IF (IPRINT .GE. 4)
     &   WRITE(LUPRI,9001) EACTFQ
C
      EACT = EACT + EACTFQ
C
      IF (IPRINT .GE. 4)
     &   WRITE(LUPRI,9002) EACT
      RENRACT = EACT
      CALL QEXIT('RENRACT')
      RETURN
 9000 FORMAT(/2X,'rENRACT: sum_{uv} DV_{uv} FC_{uv} : ',F25.15)
 9001 FORMAT(/2X,'rENRACT: sum_{u} FQ_{uu}          : ',F25.15)
 9002 FORMAT(/2X,'rENRACT: Total active energy      : ',F25.15)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck renrac1 */
      FUNCTION RENRAC1(FQ,NORBT,INASH,IOFF)
C***********************************************************************
C
C     Calculate \sum_{u} FQ_{u,u}
C
C     Input:
C        FQ      - the Q Fock matrix for a given symmetry
C
C     Returns:
C        RENRAC1 - \sum_{u} FQ_{u,u}
C
C     Written by J. Thyssen - Oct 29 1998
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
      DIMENSION FQ(NORBT,INASH)
C
      CALL QENTER('RENRAC1')
C
      E = D0
      DO I = 1,INASH
         E = E + FQ(I+IOFF,I)
      END DO
      RENRAC1 = E
C
      CALL QEXIT('RENRAC1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rreadmo */
      SUBROUTINE RREADMO(CMO,JRDMO,IOPT,LUNIT)
C***********************************************************************
C
C     Retrieve MO coefficients from file
C
C     Input:
C
C        IOPT = 1: read from label NEWORB
C        IOPT = 2: read from label OLDORB
C        IOPT = 3: read from label NEWNATOB
C        IOPT = 4: read from label MCCINATO
C
C     Output:
C        CMO     - MO coefficients
C
C     Written by J. Thyssen - Nov 3 1998
C     Last revision :
C
C                S. Knecht  - Jul 7 2008
C
C                read from label NEWNATOB (MP2, later also MCSCF...)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
#include "dcbbas.h"
C
      DIMENSION CMO(*)
      LOGICAL NOSTOP
      CHARACTER LABEL*8
C
      CALL QENTER('RREADMO')
      IF (JRDMO .LT. 0) THEN
         JRDMO = 0
         NOSTOP = .TRUE.
      ELSE
         NOSTOP = .FALSE.
      END IF
      IF (IOPT .LT. 1 .OR. IOPT .GT. 4) THEN
          WRITE(LUPRI,'(/A,I10)')
     &       ' rREADMO: ERROR: Unknown option ',IOPT
          CALL QUIT('rREADMO: ERROR: Unknown option.')
      END IF
      GOTO (100,100), IOPT
  100 CONTINUE
C
C     Read orbitals from label 'xxxxXXXX'
C
      REWIND(LUNIT)
      IERR = -1
      IF ( IOPT .EQ. 1 ) THEN
         LABEL = 'NEWORB  '
      ELSE IF( IOPT .EQ. 2 )THEN
         LABEL = 'OLDORB  '
      ELSE IF( IOPT .EQ. 3 )THEN
         LABEL = 'NEWNATOB'
      ELSE IF( IOPT .EQ. 4 )THEN
         LABEL = 'MCCINATO'
      END IF
C
      CALL MOLLAB(LABEL,LUNIT,IERR)
      IF (IERR .GE. 0) THEN
         CALL READT(LUNIT,NCMOTQ,CMO)
      ELSE IF (IERR .EQ. -1) THEN
         IF (NOSTOP) THEN
            JRDMO = JRDMO + 1
         ELSE
            CALL QUIT('rREADMO: ERROR: Label ' // LABEL //
     &                'not found on file KRMCSCF')
         END IF
      ELSE
         IF (NOSTOP) THEN
            JRDMO = JRDMO + 1
         ELSE
            CALL QUIT('rREADMO: ERROR: error reading file "KRMCSCF"')
         END IF
      END IF
      IF ( JRDMO .GT. 0 ) THEN
         GO TO 1010
      ELSE
         GO TO 1000
      END IF
C
C     End:
C
 1000 CONTINUE
      IPROPT_SAVE = IPROPT
csk   IPROPT      = 10
      IF (IPROPT .GE. 10) THEN
         CALL HEADER('Output from rREADMO:',-1)
         DO I = 1,NFSYM
            CALL HEADER('Coefficients for irrep '//CHAR(I+48),-1)
            WRITE(LUPRI,*) ' printing offset ICMOQ(I)+1',ICMOQ(I)+1
            WRITE(LUPRI,*) ' NFBAS(I,0), NORB(I)',NFBAS(I,0),NORB(I)
            CALL PRQMAT(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
      IPROPT = IPROPT_SAVE
 1010 CONTINUE
      CALL QEXIT('RREADMO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck roptst */
      SUBROUTINE ROPTST(ISTART,IWRT,CMO,CREF,IBEIG,IPRINT,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate start vectors
C
C     Input:
C        ISTART
C        IWRT
C
C     Output: None
C
C
C     Written by J. Thyssen - Nov 17 1998
C
C     Last revision :
C                S. Knecht  - Jul 07 2008  start plain CI calculation
C                             from MCSCF/MP2 natural orbitals
C                             (RCI_MCORB = .TRUE.)
C                S. Knecht  - Aug 10 2010  read boson symmetry array
C                             IBEIG (if spinfr_krmc == .true.)
C                             --> reuse mj-related logicals.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "maxorb.h"
#include "dcbbos.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbnrt.h"
C
      LOGICAL EX, FNDLAB, FILEOPEN, MJVEC_AVAIL, readmj_df
      DIMENSION WORK(*)
      DIMENSION CREF(*), CMO(*), IBEIG(*)
      LOGICAL VERBOSE
C
      LOGICAL LBIT, RCI_MCORB, RTRACI, NATOLCR_SAVE, READ_DFCOEF
      CHARACTER SZTRA(2,4)*80
      CHARACTER CNOLBL*6, KRM5_NAME*800, FILELABEL*8
C
      CALL QENTER('ROPTST')
      KFRSAV = KFREE

      RCI_MCORB   = .FALSE.
      MJVEC_AVAIL = .FALSE.
      FILEOPEN    = .FALSE.
      readmj_df   = .false.
C
C     Start timer
C
      CALL RTKTIME(.TRUE.,1)
C
C     IOPTST controls the start guess:
C     ------------------------------
C
      IROST  = MOD(ISTART,100)
      IRCIST = ISTART/100
C
C     IRCIST controls the CI start guess:
C     -----------------------------------
C
C     bit 1 of IRCIST:  read CREF from KRMCSCF
C     bit 2 of IRCIST:  read CREFOLD from KRMCSCF
C     bit 3 of IRCIST:  do MAXCIT CI iterations (MAXCIT<0 ==> full CI)
C     bit 4 of IRCIST:  1/Hdiag
C
C     IROST controls the orbital start guess:
C     ---------------------------------------
C
C     bit 1 of IROST: read NEWNATOB on LUKRM5 (KRMCSCF or KRMCOLD)
C     bit 2 of IROST: read NEWORB   on LUKRM5 (KRMCSCF or KRMCOLD)
C     bit 3 of IROST: read OLDORB   on LUKRM5 (KRMCSCF or KRMCOLD)
C     bit 4 of IROST: read DFCOEF from DHF program
C     bit 5 of IROST: start         on H1DIAG
C
      IF (IPROPT .GE. 2) THEN
         WRITE(LUPRI,'(4(/A,I4))')
     &        ' (ROPTST)  Control parameter   (ISTART) =',ISTART,
     &        '           Orbital start guess (IROST)  =',IROST,
     &        '           Config. start guess (IRCIST) =',IRCIST,
     &        '           Write CMO & CREF    (IWRT)   =',IWRT
      END IF
C
C     For KRMCSCF, LUKRM5 will be KRMCOLD from an old MCSCF calculation,
C         if available.
C     For KR-CI, LUKRM5 will be KRMCSCF from preceeding MCSCF or KRMCOLD
C         from an old MCSCF calculation, if either is available.
C     (Aug 08, SK+HJAaJ)
C

      IF ( LUKRM5 .GT. 0 ) THEN
         INQUIRE( UNIT=LUKRM5, NAME=KRM5_NAME )
         LKRM5_NAME = 0
         LKRM5_NAME = LEN_TRIM(KRM5_NAME)
      ELSE
         KRM5_NAME  = '(file does not exist)'
         LKRM5_NAME = 21
      END IF
C
C     *******************************
C     *** (1) Orbital start guess ***
C     *******************************
C
      IF (IROST .EQ. 0) GOTO 1100
C
C     (1.1) Read MCSCF natural orbitals from label MCCINATO on KRM5_NAME
C
      IF (LBIT(IROST,1)) THEN
         IF ( LUKRM5 .GT. 0 ) THEN
            JRDMO = -1
            CALL RREADMO(CMO,JRDMO,4,LUKRM5)
            IF (JRDMO .EQ. 0) THEN
               IF (IPRINT .GE. 0) THEN
                  WRITE(LUPRI,'(/2A)')
     &                 ' (ROPTST)  Start orbitals read from' //
     &                 ' label MCCINATO on file ',
     &                   KRM5_NAME(1:LKRM5_NAME)
               END IF
C
               IF ( MOPT_MXMACRO .EQ. -1 ) THEN
                 WRITE(LUPRI,'(/A)')
     &                ' (ROPTST)  Starting CI calculation' //
     &                ' from MC/CI natural orbitals! '
                 RCI_MCORB = .TRUE.
               ENDIF
               GO TO 1000
            END IF
         END IF
      END IF
C
C     (1.2) Read MP2 natural orbitals from label NEWNATOB on KRM5_NAME
C
      IF (LBIT(IROST,2)) THEN
         IF ( LUKRM5 .GT. 0 ) THEN
            JRDMO = -1
            CALL RREADMO(CMO,JRDMO,3,LUKRM5)
            IF (JRDMO .EQ. 0) THEN
               IF (IPRINT .GE. 0) THEN
                  WRITE(LUPRI,'(/2A)')
     &                 ' (ROPTST)  Start orbitals read from' //
     &                 ' label NEWNATOB on file ',
     &                   KRM5_NAME(1:LKRM5_NAME)
               END IF
C
               IF( MOPT_MXMACRO .EQ. - 1 ) THEN
                 WRITE(LUPRI,'(/A)')
     &                ' (ROPTST)  Starting CI calculation' //
     &                ' from MC orbitals! '
                 RCI_MCORB = .TRUE.
               ENDIF
               GO TO 1000
            END IF
         END IF
      END IF
C
C     (1.3) Read from label NEWORB on KRM5_NAME
C
      IF (LBIT(IROST,3)) THEN
         IF ( LUKRM5 .GT. 0 ) THEN
            JRDMO = -1
            CALL RREADMO(CMO,JRDMO,1,LUKRM5)
            IF (JRDMO .EQ. 0) THEN
               IF (IPRINT .GE. 0) THEN
                  WRITE(LUPRI,'(/2A)')
     &                 ' (ROPTST)  Start orbitals read from' //
     &                 ' label NEWORB on file ',KRM5_NAME(1:LKRM5_NAME)
               END IF
               IF ( GENFOCK )THEN
                  REWIND(LUKRM5)
                  IF ( FNDLAB('CREF    ',LUKRM5) ) THEN
                     IF (MOPT_MXMACRO .ge. 0 ) THEN
                        CALL REAKRMC(LUKRM5,'CREF    ',CREF,NZCONFQ)
                     END IF
                  END IF
C                 ... generate LUCIAREL block/structure information
                  CALL LUCI_RSTRMC(WORK(KFREE),LFREE)
                  CNOLBL  = 'ONLYFD'
                  RTRACI  = .FALSE.
C                 ... if CI Nat. Orb. module shall be used later
                  NATOLCR_SAVE = NATOLCR
                  NATOLCR = .FALSE.
                  CALL MEMGET2('REAL','EIG',KEIG,NORBT,WORK,KFREE,LFREE)
                  CALL DZERO(WORK(KEIG),NORBT)
C                 ... transformation to Fock-type orbitals
                  CALL KRMCCNO(CNOLBL,RTRACI,CREF,CMO,IPRINT,WORK(KEIG),
     &                         WORK,KFREE,LFREE)
                  NATOLCR = NATOLCR_SAVE
C
                  CALL MEMREL('after FOCK-type',WORK,1,KFRSAV,
     &                         KFREE,LFREE)
               END IF
C
               IF( MOPT_MXMACRO .EQ. - 1 ) THEN
                 WRITE(LUPRI,'(/A)')
     &                ' (ROPTST)  Starting CI calculation' //
     &                ' from MC orbitals! '
                 RCI_MCORB = .TRUE.
               ENDIF
               GO TO 1000
            END IF
         END IF
      END IF
C
C     (1.4) Read from label OLDORB on KRM5_NAME
C
      IF (LBIT(IROST,4)) THEN
         IF ( LUKRM5 .GT. 0 ) THEN
            JRDMO = -1
            CALL RREADMO(CMO,JRDMO,2,LUKRM5)
            IF (JRDMO .EQ. 0) THEN
               IF (IPRINT .GE. 0) THEN
                  WRITE(LUPRI,'(/2A)')
     &                 ' (ROPTST)  Start orbitals read from' //
     &                 ' label OLDORB on file ',KRM5_NAME(1:LKRM5_NAME)
               END IF
               GO TO 1000
            END IF
         END IF
      END IF
C
C     (1.5) From previous DHF calculation.
C           Read orbitals from DFCOEF
C
      IF (LBIT(IROST,5)) THEN
         INQUIRE ( FILE = 'DFCOEF', EXIST = EX )
         IF ( EX ) THEN
            CALL OPNFIL(LUCOEF,'DFCOEF','OLD','RREADMO')
            CALL REACMO(LUCOEF,'DFCOEF',CMO,DUMMY,DUMMY,DUMMY,2)
            CLOSE(LUCOEF,STATUS='KEEP')
            readmj_df = .true.
            IF (IPRINT .GE. 0) THEN
               WRITE(LUPRI,'(/A)')
     &              ' (ROPTST)  Start orbitals read from' //
     &              ' file DFCOEF'
            END IF
            GO TO 1000
         END IF
      END IF
C
C     (1.6) From a usual DHF start guess, i.e. bare nucleus guess.
C           TODO FIXME: to be implemented
C
      IF (LBIT(IROST,6)) THEN
         CALL RH1DIAG(CMO,WORK(KFREE),LFREE)
         IF (IPRINT .GE. 0) THEN
            WRITE(LUPRI,'(/A)')
     &      ' (ROPTST)  Start orbitals from' //
     &      ' diagonalization of H_1'
            IF (BNCRON) WRITE(LUPRI,'(A/)')
     &      "           Improved by an estimate "//
     &      "of the electronic screening (Slater's rules)."
         END IF
         CALL QUIT('*** ERROR in ROPTST: H_1 start guess is
     &             not operational. ***')
         GOTO 1000
      END IF
C
C     Oops, nothing worked:
C
      WRITE(LUPRI,'(//3A,I3)')
     &     '*** ERROR in ROPTST *** ',
     &     'Failed to generate orbital start guess; ',
     &     'value of IROST is ',IROST
      CALL QUIT('*** ERROR in ROPTST *** No orbital start guess')
C
 1000 CONTINUE
C
C     Make sure start guess orbitals are orthonormal !!!
C
      IF (IPROPT .GE. 10) THEN
         DO I = 1, NFSYM
            WRITE(LUPRI,'(A,I2)') ' (ROPTST)  Coefficients, irrep ',I
            CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &           NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
C
      CALL RGMCMO(CMO,WORK(KFREE),LFREE)
C
      IF (IPROPT .GE. 10) THEN
         DO I = 1, NFSYM
            WRITE(LUPRI,'(A,I3)')
     &         ' Coefficients (after Gram-Schmidt), irrep',I
            CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &           NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
!
!     intermediate step:
!     save additional required information (integer arrays) to file KRMCSCF

!     1. store mj-values (MCSCF run in linear symmetry)
!     -------------------------------------------------
      if(.not.spinfr_krmc)then

        IF ( LUKRM5 .GT. 0 ) THEN
          REWIND(LUKRM5)
          FILELABEL = 'MJVEC   '
          MJVEC_AVAIL = FNDLAB(FILELABEL,LUKRM5)
          REWIND(LUKRM5)
        END IF
        if( linear .and. readmj_df ) MJVEC_AVAIL = .true.
        IF( MJVEC_AVAIL )THEN
          CALL MEMGET2('INTE','MJVEC',KMJVEC,NORBT,WORK,KFREE,LFREE)
          CALL IZERO(WORK(KMJVEC),NORBT)
!
!         read mj-vector from DFCOEF or KRMC"SCF"(OLD)
          if( readmj_df )then
            call opnfil(lucoef,'DFCOEF','OLD','RMJVEC ')
            call reacmo(lucoef,'DFCOEF',dummy,dummy,work(kmjvec),
     &                  dummy,8)
            close(lucoef,STATUS='KEEP')
          else
            call ireakrmc(lukrm5,'MJVEC   ',work(kmjvec),norbt)
          end if
C
C         write back on file - if not present
          LUNI = -1
          INQUIRE(FILE='KRMCSCF',EXIST=EX,OPENED=FILEOPEN,NUMBER=LUNI)
          IF(EX.and.FILEOPEN)THEN
            rewind(LUNI)
            if(.not.(FNDLAB(FILELABEL,LUNI)))
     &        CALL IWRTKRMC(LUNI,'MJVEC   ',WORK(KMJVEC),NORBT)
          ELSE IF( EX.and.(.not. FILEOPEN))THEN
            LUNI = 99
            CALL OPNFIL(LUNI,'KRMCSCF','OLD    ','ROPTST')
            rewind(LUNI)
            if(.not.(FNDLAB(FILELABEL,LUNI)))
     &        CALL IWRTKRMC(LUNI,'MJVEC   ',WORK(KMJVEC),NORBT)
            CLOSE(LUNI,STATUS='KEEP')
          ELSE IF(.NOT. EX)THEN ! can only happen if /= MCSCF run
            LUNI = 99
            CALL OPNFIL(LUNI,'KRMCSCF','UNKNOWN','KR CI ')
            CALL NEWLAB('*KRCI   ',LUNI,LUPRI)
            CALL NEWLAB('SODLABEL',LUNI,LUPRI)
            CALL NEWLAB('EOFLABEL',LUNI,LUPRI)
            CALL IWRTKRMC(LUNI,'MJVEC   ',WORK(KMJVEC),NORBT)
            CLOSE(LUNI,STATUS='KEEP')
          END IF

          CALL MEMREL('KRMCSAV.1',WORK,1,KFRSAV,KFREE,LFREE)
          WRITE(LUPRI,'(/A)')
     &    '  Vector containing mj-values saved as MJVEC on KRMCSCF'
        END IF

!     2. store IBEIG array containing boson irreps of each orbital (spinfree MCSCF run)
!     ---------------------------------------------------------------------------------
      else

        IF ( LUKRM5 .GT. 0 ) THEN
          REWIND(LUKRM5)
          FILELABEL = 'IBEIG   '
          MJVEC_AVAIL = FNDLAB(FILELABEL,LUKRM5)
          REWIND(LUKRM5)
        END IF
        if( readmj_df ) MJVEC_AVAIL = .true.
        IF( MJVEC_AVAIL )THEN
!
!         read ibeig (info array on boson irreps) from DFCOEF or KRMC"SCF"(OLD)
          if( readmj_df )then
            call opnfil(lucoef,'DFCOEF','OLD','RIBEIG ')
            call reacmo(lucoef,'DFCOEF',dummy,dummy,ibeig,dummy,8)
!           WRITE(LUPRI,'(/A)') ' IBEIG array read from file DFCOEF'
!           call iwrtmamn(ibeig,1,norbt,1,norbt,lupri)
            close(lucoef,STATUS='KEEP')
          else
            call ireakrmc(lukrm5,'IBEIG   ',ibeig,norbt)
          end if
C
C         write back on file - if not present
          LUNI = -1
          INQUIRE(FILE='KRMCSCF',EXIST=EX,OPENED=FILEOPEN,NUMBER=LUNI)
          IF(EX.and.FILEOPEN)THEN
            rewind(LUNI)
            if(.not.(FNDLAB(FILELABEL,LUNI)))
     &        CALL IWRTKRMC(LUNI,'IBEIG   ',ibeig,NORBT)
          ELSE IF( EX.and.(.not. FILEOPEN))THEN
            LUNI = 99
            CALL OPNFIL(LUNI,'KRMCSCF','OLD    ','ROPTST')
            rewind(LUNI)
            if(.not.(FNDLAB(FILELABEL,LUNI)))
     &        CALL IWRTKRMC(LUNI,'IBEIG   ',ibeig,NORBT)
            CLOSE(LUNI,STATUS='KEEP')
          END IF

          WRITE(LUPRI,'(/A)')
     &    '  Vector containing the boson irreps saved as IBEIG'//
     &    ' on KRMCSCF'
        END IF
      end if
!
!     final step in orbital part: write coefficients back to file as label NEWORB
!     ---------------------------------------------------------------------------
      IF(LBIT(IWRT,1)) THEN
        CALL WRTKRMC(LUKRMC,'NEWORB  ',CMO,NCMOTQ)
         WRITE(LUPRI,'(/A)')
     &   '  Orbitals saved as NEWORB on KRMCSCF'
      END IF
C
 1100 CONTINUE
      CALL FLSHFO(LUPRI)
C
C     *****************************
C     *** (2) Conf. start guess ***
C     *****************************
C
      IF( GENFOCK ) THEN
        IF (LBIT(IWRT,2)) CALL WRTKRMC(LUKRMC,'CREF    ',CREF,NZCONFQ)
CSK     GOTO 2100
      END IF
C
C
C
C     bit 1 of IRCIST:  read CREF from KRMCSCF
C     bit 2 of IRCIST:  read CREFOLD from KRMCSCF
C     bit 3 of IRCIST:  do MAXCIT CI iterations (MAXCIT<0 ==> full CI)
C     bit 4 of IRCIST:  1/Hdiag
C
      IF (IRCIST .EQ. 0) GOTO 2100
C
      IF ( MCTYPE .GE. JMCMIN ) THEN
C
         IF( RCI_MCORB ) GOTO 1250
C
C        (2.1) Read CREF from KRM5_NAME
C              or restart plain LUCIAREL from CI vector
C
         IF (LBIT(IRCIST,1)) THEN
            IF ( LUKRM5 .GT. 0 ) THEN
               REWIND(LUKRM5)
               IF ( FNDLAB('CREF    ',LUKRM5) ) THEN
                  IF (MOPT_MXMACRO.GT.0) THEN
                     CALL REAKRMC(LUKRM5,'CREF    ',CREF,NZCONFQ)
                     IF (IPRINT .GE. 0) THEN
                        WRITE(LUPRI,'(/2A)')
     &                       ' (ROPTST)  CI start vector read from' //
     &                       ' label CREF on file',
     &                         KRM5_NAME(1:LKRM5_NAME)
                     END IF
                  END IF
                  IF (OPT_CIPROGRAM.EQ.'LUCIAREL') THEN
C
                    IF (MOPT_MXMACRO.GT.0) THEN
C
C                     RESTART from KRM5_NAME in MCSCF run
C                     Generate restart and vector blocking information
C                     for LUCIAREL CI vectors:
C
                      WRITE(LUPRI,'(/A)')
     &                      ' (ROPTST)  Restarting KRMC with' //
     &                      ' LUCIAREL CI module'
                      CALL LUCI_RSTRMC(WORK(KFREE),LFREE)
                    ELSE IF (MOPT_MXMACRO.EQ.0) THEN
                      WRITE(LUPRI,*)
     &                '  Starting CI calculation from MC orbitals or '
                      WRITE(LUPRI,*)
     &                '  Restarting CI calculation from vector. '
                      IOPT_INTFLG = IOPT_INTDEF
                      JTRLVL = 0
                      IF( (NATOLCR .or. IANACI .eq. 1) .and. MAXCIT .le.
     &                     0 )THEN
C
C                       RESTART from KRM5_NAME in MCSCF run
C                       Generate restart and vector blocking information
C                       for LUCIAREL CI vectors:
C
                        WRITE(LUPRI,'(/A)')
     &                        ' (ROPTST)  Restarting KRMC with' //
     &                        ' LUCIAREL CI module'
                        CALL LUCI_RSTRMC(WORK(KFREE),LFREE)
                        CALL REAKRMC(LUKRM5,'CREF    ',CREF,NZCONFQ)
                        IF (IPRINT .GE. 0) THEN
                           WRITE(LUPRI,'(/2A)')
     &                         ' (ROPTST)  CI start vector read from' //
     &                         ' label CREF on file ',
     &                           KRM5_NAME(1:LKRM5_NAME)
                        END IF
                      END IF
                      CALL INISZT(JTRLVL,SZTRA,LNOROT,NOROT)
                      CALL RTRACTL(CMO,IBEIG,SZTRA,WORK(KFREE),LFREE)
                      CALL RCIST(CMO,CREF,WORK,KFREE,LFREE)
                    ELSE
                      CALL QUIT(' MOPT_MXMACRO .LT. 0, Stopping. ')
                    END IF
                  END IF
                  GO TO 2000
               END IF
            END IF
         END IF
C
C        (2.2) Read CREFOLD from KRM5_NAME
C
         IF (LBIT(IRCIST,2)) THEN
            IF ( LUKRM5 .GT. 0 ) THEN
               IF ( FNDLAB('CREF    ',LUKRM5) ) THEN
                  REWIND(LUKRM5)
                  CALL REAKRMC(LUKRM5,'CREFOLD ',CREF,NZCONFQ)
                  IF (IPRINT .GE. 0) THEN
                     WRITE(LUPRI,'(/2A)')
     $                    ' (ROPTST)  CI start vector read from' //
     $                    ' label CREFOLD on file ',
     &                      KRM5_NAME(1:LKRM5_NAME)
                  END IF
                  GO TO 2000
               END IF
            END IF
         END IF
C
C        (2.3) Do MAXCIT CI iterations
C
1250     CONTINUE
C
         IF (LBIT(IRCIST,3)) THEN
C
C           Do a 0th order 4-index transformation.
C
            IOPT_INTFLG = IOPT_INTDEF
            JTRLVL = 0
            CALL INISZT(JTRLVL,SZTRA,LNOROT,NOROT)
C
            IF (OPT_CIPROGRAM.EQ.'LUCIAREL') THEN
C
               EX = .FALSE.
C
C              check for existing integral files
C
               IF( CINT_REORD ) THEN
C
                 INQUIRE ( FILE = 'IJKL_REOD', EXIST = EX )
C
                 IF ( EX ) GOTO 1500
C
               END IF
C
               INQUIRE ( FILE = '4INDINFO', EXIST = EX )
               IF ( .NOT. EX ) THEN
                 IF (KTRLVL.EQ.5)
     &           WRITE(LUPRI,*) 'CALL RTRACTL for transformation only !'
                 CALL RTRACTL(CMO,IBEIG,SZTRA,WORK(KFREE),LFREE)
               END IF
            ELSE
               CALL RTRACTL(CMO,IBEIG,SZTRA,WORK(KFREE),LFREE)
            END IF
C
            IF (KTRLVL.EQ.5) GOTO 2100
C
 1500       CONTINUE
C
            CALL RCIST(CMO,CREF,WORK,KFREE,LFREE)
C
            IF (IPRINT .GE. 0) THEN
               IF (MAXCIT .GE. 0) THEN
                  IF( JKRRUNTYPE .gt. 1)THEN
                     WRITE(LUPRI,'(A)')
     &               ' (ROPTST)  KR-CI calculation ended properly.'
                  ELSE
                     WRITE(LUPRI,'(/A,I5,A)')
     &                    ' (ROPTST)  CI start vector generated from',
     &                    MAXCIT,' CI iterations.'
                  END IF
               ELSE
                  IF( JKRRUNTYPE .gt. 1)THEN
                     WRITE(LUPRI,'(/A)')
     &                    ' (ROPTST)  Full-CI calculation succeeded.'
                  ELSE
                     WRITE(LUPRI,'(/A)')
     &                    ' (ROPTST)  CI start vector generated from'//
     &                    ' diag. of full CI matrix.'
                  END IF
               END IF
            END IF
            GO TO 2000
C
         END IF
C
C        (2.4) 1/Hdiag
C
         IF (LBIT(IRCIST,4)) THEN
            CALL QUIT('*** ERROR in ROPTST *** 1/Hdiag CI start guess
     & is not implemented')
         END IF
C
C        Oops, nothing worked:
C
         WRITE(LUPRI,'(//3A,I3)')
     &        '*** ERROR in ROPTST *** ',
     &        'Failed to generate CI start guess; ',
     &        'value of IRCIST is ',IRCIST
         CALL QUIT('*** ERROR in ROPTST *** No CI start guess')
C
 2000    CONTINUE
C
C        Write CREF to file.
C
         IF (IPRINT .GE. 30) THEN
            MZ = MIN(NZ,2)
            WRITE(LUPRI,'(/A)') ' (ROPTST)  CI start vector:'
            CALL OUTPUT(CREF,1,NZCONF,1,MZ,NZCONF,MZ,1,LUPRI)
         END IF
         IF (LBIT(IWRT,2)) THEN
           CALL WRTKRMC(LUKRMC,'CREF    ',CREF,NZCONFQ)
           WRITE(LUPRI,'(/A)')
     &     '  CI vector saved as CREF on KRMCSCF'
         END IF
C
      END IF
C     ... check for MCTYPE .ne. SCF
C
 2100 CONTINUE
C
      CALL MEMREL('ROPTST',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL RTKTIME(.FALSE.,1)
      CALL QEXIT('ROPTST')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rcist */
      SUBROUTINE RCIST(CMO,CREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate CI start vectors.
C
C     Input : CMO - MO coefficents
C
C     Output: CREF - CI vectors
C
C
C     Written by J. Thyssen - Jun 27 2000
C
C***********************************************************************
      use mc_energies
      use qcorr_cfg
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "krmcluci_inf.h"
C
      DIMENSION CMO(*), CREF(*), WORK(*)
      LOGICAL   MOISNO, EX
C
      CALL QENTER('RCIST')
      CALL RTKTIME(.TRUE.,2)
      KFRSAV = KFREE
C
C     Do some (e.g., 3) Davidson CI iterations to get a start vector:
C
C     MAXCIT >= 0:
C
C     (2.1) Contruct diagonal of CI matrix
C     (2.2) Pick start vector with zeros everywhere except for the
C           lowest value of the diagonal of the CI matrix. If 2 or
C           more diagonal values are close, pick more start vectors.
C     (2.3) Do a couple of Davidson CI iterations.
C
C     MAXCIT < 0:
C
C     Diagonalize the full CI matrix now.
C     This is of course only possible for sufficiently small CI's.
C
!     fcidump: generate FCIDUMP file and return
      if(fcidump)then
        call fcidump_driver(CMO,WORK,KFREE,LFREE,IPROPT)
        goto 1000
      end if

      IF (MAXCIT .LT. 0 .AND. NZCONF .GT. 100000) THEN
         WRITE(LUPRI,'(//A,I0,A/A)')
     &      'INFO: Are you mad? Full CI with ',NZCONF,' determinants?',
     &      'INFO: Using MAXCIT = 3 instead!'
         MAXCIT = 3
      END IF
C
      IF ( OPT_CIPROGRAM .EQ. 'GASCIP' ) THEN
C
         IF (MAXCIT .GE. 0) THEN
C
C           MAXCIT Davidson CI iterations.
C
            CALL MEMGET2('REAL','ECI',KECI,NCIROOT,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KECI),NCIROOT)
            THRGCI = 0.1D0 * SQRT(OPT_THRECI)   ! GASCIP uses gradient threshold, not energy threshold
            CALL GASCIP_RCISTD(NCIROOT,THRGCI,MAXCIT,
     &                         WORK(KECI),CREF,NZCONF,WORK(KZCONF),CMO,
     &                         OPT_THR_CVEC,IOPT_STATE,ICONV_TOT,
     &                         THRPCI,IPROPT,WORK(KFREE),LFREE)

C
            ECI = WORK(KECI)
         ELSE
C
C           Full CI.
C
            CALL MEMGET2('REAL','ECI',KECI,NZCONF,WORK,KFREE,LFREE)
            CALL GASCIP_STARTVEC(WORK(KZCONF),CMO,CREF,ECI,
     &                           NZCONF,IOPT_STATE,
     &                           THRPCI,IPROPT, WORK,KFREE,LFREE)
C
         END IF
         CALL MEMREL('RCIST.after GASCIP start vec',WORK,
     &               1,KECI,KFREE,LFREE)
         IF (IANACI.EQ.1) THEN
           CALL GASCIP_ANACI(CREF,NZCONF,WORK(KZCONF),THRPCI,
     &                       WORK(KFREE),LFREE)
         END IF

C
      ELSE IF ( OPT_CIPROGRAM .EQ. 'LUCIAREL' ) THEN
C
         IF( CINT_REORD ) THEN
C
           INQUIRE (FILE='IJKL_REOD',EXIST=EX)
           IF(.NOT.EX)
     &     CALL LUCI_IJKLRO(CMO,WORK(KFREE),LFREE)
C
           IF( MAXCIT .le. 0 ) THEN
C
             IF( NATOLCR .or. ( IANACI .eq. 1 ) .or. luci_cfg_qcorr)then
               GOTO 500
             ELSE
               GOTO 1000
             ENDIF
           ENDIF
C
         END IF
C
         CALL LUCI_RCIST(CREF,CMO,WORK(KFREE),LFREE)
         eci = eci_init
C
         call memchk('after CALL LUCI_RCIST',WORK,KFREE)

 500     if(luci_cfg_qcorr)then
           CALL LUCI_QCORR(CREF,WORK,KFREE,LFREE)
         end if

         IF (IANACI.EQ.1) THEN
           CALL LUCI_ANAL(CREF,WORK,KFREE,LFREE)
         END IF

C
      ELSE IF ( OPT_CIPROGRAM .EQ. 'KRCC' ) THEN
         IF (IPROPT.GT.5)
     &      WRITE(LUPRI,*) 'KRCC: Skipping CI, will do CC ;-)'
      ELSE
         WRITE(LUPRI,'(//3A)')
     &        '*** ERROR in RCIST ***' //
     &        'nothing implemented for CI program "',
     &         OPT_CIPROGRAM,'"'
         CALL QUIT('*** ERROR in RCIST *** Not implemented')
      END IF
C
!     WRITE(LUPRI,'(/A,F20.8)') ' (RCIST) CI energy:',ECI
C
C     General allocation for NO generation
C
      IDENSLR_STATE = 0
      IPRNO = 3*IPROPT
      IF( NATOLCR ) THEN
        CALL MEMGET2('REAL','OCCNO',KOCCNO,NORBT  ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','UNO',  KUNO  ,N2ASHXQ,WORK,KFREE,LFREE)
C
        WRITE(LUPRI,'(//A/A/A)')
     &'  *************************************************************',
     &'  ************* Natural orbital occupation numbers ************',
     &'  *************************************************************'
C
C
        DO ICNON = 1, NCIROOT
C
C         erase old occupation numbers
          CALL DZERO(WORK(KOCCNO),NORBT)
          CALL DZERO(WORK(KUNO),N2ASHXQ)
C
          IDENSLR_STATE = ICNON
C
C         NOTE: IF MOPT_MXMACRO < 0: CREF is not allocated!
C               requires appropriate treatment in LUCI_DENS ...
C
C         SK - Feb 2009: introduced TRA_NATO for transformation to
C                        natural CI/MCSCF orbitals. They can be useful
C                        for further analysis, e.g., visualization.
C
          CALL RGETNO(CREF,WORK(KOCCNO),WORK(KUNO),CMO,
     &                TRA_NATO,MOISNO,.false.,IPRNO,WORK,KFREE,LFREE)
C
        END DO
C
C        set back if KR-MCSCF - otherwise LUCI_DENS call may fail
C        do NOT set back for KR-CI - we might do more than one symmetry
C
        IF( JKRRUNTYPE .lt. 2 ) NATOLCR = .FALSE.
C
      END IF
C
1000  CONTINUE
C
      CALL MEMREL('RCIST.after RGETNO',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL RTKTIME(.FALSE.,2)
      CALL QEXIT('RCIST')
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirnr */
      SUBROUTINE DIRNR(CREF,GCI,GOEE,GOEP,CMO,IBEIG,DV,FC,FV,FCAC,
     &                 H2AC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Perform NR step
C
C     Input:
C        CMO     - the orbitals for the CEP
C        IBEIG   - the boson irreps arry of the orbitals for the CEP
C        CREF    - CI vector
C        G       - the gradient
C         CI       CI
C         OEE      orbital electronic
C         OEP      orbital positronic
C
C     Output:
C
C     Written by J. Thyssen and H. J. Aa. Jensen - Nov 17 1998
C
C***********************************************************************
      use orbital_rotation_indices
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      PARAMETER (D0 = 0.00D00, D1 = 1.00D00)
C
#include "maxorb.h"
#include "dcbnrt.h"
#include "dcbopt.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbkrmc_itinfo.h"
C
      LOGICAL OPTRST,OPTNRM,OPTDIH,OPTSTERN, GOEER
      DIMENSION CREF(*),GCI(*),DV(*),FC(*),FV(*),
     &          H2AC(*),WORK(*),GOEE(*),GOEP(*), IBEIG(*)
      DIMENSION ITRINT_MIC(2),CNVINT_MIC(2)
      DIMENSION IOPTTB(-1:2)
      DATA IOPTTB / 0,0,0,0 /
      integer, allocatable :: save_pp(:)
      integer, allocatable :: save_pn(:)
      integer, allocatable :: compress_pp(:)
      integer, allocatable :: compress_pn(:)
C
      CHARACTER SZTRA(4,2)*200
C
      CALL QENTER('DIRNR')
      KFRSAV = KFREE
C
C     Calculate initial micro it. convergence threshold
C
      RBETA = D0
      THR = RCNVTHR()
C
      WRITE(LUPRI,'(/A,1P,D10.2/)')
     &'(DIRNR) Starting micro iterations. Convergence threshold:',THR
C
C
C     ***************************************************
C     *** Compress orbital rotations ********************
C     ***************************************************
C
C
      IF (COMPROT) THEN
C
C        Save lengths and pointers.
C
         NZXOPE_SAV = NZXOPE
         NZHOPE_SAV = NZHOPE
         NZXOPPQ_SAV = NZXOPPQ
         NZXOPEQ_SAV = NZXOPEQ
         NZHOPEQ_SAV = NZHOPEQ
         NZXOPP_SAV = NZXOPP

C        Compress orbital rotations.
C        (this will destroy GPOE and GPOP but that's OK).
C
         if (nzhope > 0) then
            allocate(save_pp(nzhope*2))
            allocate(compress_pp(nzhope*2))
            call icopy(nzhope*2,
     &                 get_orbital_rotation_indices_pp(), 1,
     &                 save_pp, 1)
            call icopy(nzhope*2,
     &                 get_orbital_rotation_indices_pp(), 1,
     &                 compress_pp, 1)
         else
            allocate(save_pp(1))
            allocate(compress_pp(1))
         end if

         if (nzxopp > 0) then
            allocate(save_pn(nzxopp*2))
            allocate(compress_pn(nzxopp*2))
            call icopy(nzxopp*2,
     &                 get_orbital_rotation_indices_pn(), 1,
     &                 save_pn, 1)
            call icopy(nzxopp*2,
     &                 get_orbital_rotation_indices_pn(), 1,
     &                 compress_pn, 1)
         else
            allocate(save_pn(1))
            allocate(compress_pn(1))
         end if
C
C        Remove rotations.
C
         CALL RMROT(NOROTC,GOEE,compress_pp,
     &              NZXOPE,NZHOPE,
     &              GOEP,compress_pn,
     &              NZXOPP,IPROPT,WORK(KFREE),LFREE)

         call update_orbital_rotation_indices(nzhope,
     &                                        nzxope,
     &                                        nzxopp,
     &                                        compress_pp,
     &                                        compress_pn)

         deallocate(compress_pp)
         deallocate(compress_pn)
C
         NZXOPPQ = NZXOPP*NZ
         NZXOPEQ = NZXOPE*NZ
         NZHOPEQ = NZHOPE*NZ
C
C        Calculate integrals with only needed orbitals.
C
         IF (MCTYPE .GE. JMCMIN) THEN
            JTRLVL = KTRLVL
            CALL INISZT(JTRLVL,SZTRA,.TRUE.,NOROTC)
            CALL RTRACTL(CMO,IBEIG,SZTRA,WORK(KFREE),LFREE)
         END IF
C
      END IF
C
C
C     ***************************************************
C     *** Transfer info to Response common blocks *******
C     ***************************************************
C
C
      OPTRST = .FALSE.
      OPTNRM = .FALSE.
CDEBUG
      OPTNRM = .TRUE.
CEND DEBUG
      OPTDIH = .TRUE.
      OPTSTERN = .FALSE.
C
C     Calculate with the same integrals used for gradient and energy.
C
      INTDEF_MIC = IOPT_INTDEF
      INTDEF_MIC = IOPT_INTFLG
      ITRINT_MIC(1) = 1
      ITRINT_MIC(2) = 1
      CNVINT_MIC(1) = DUMMY
      CNVINT_MIC(2) = DUMMY
C
      IF ( FROZCI ) THEN
C
C        We freeze the CI coefficients this macro iteration.
C        We have to include the otherwise redudant active-active rotations
C        in orbital space.
C
         NZCONF1 = 0
         NZCONFQ1 = 0
C
c        TODO FIXME: this requires many changes:
C        - in GTNORB
C        - in KRMCSAV
C        - possibly other places...
c        NZXOPE1 = NZHOPE
c        NZXOPEQ1 = NZHOPEQ
c        NZXOPT1 = NZHOPT
c        NZXOPTQ1 = NZHOPTQ
         NZXOPE1 = NZXOPE
         NZXOPEQ1 = NZXOPEQ
         NZXOPT1 = NZXOPT
         NZXOPTQ1 = NZXOPTQ
      ELSE
C
C        ``Normal'' microiterations.
C
         NZCONF1 = NZCONF
         NZCONFQ1 = NZCONFQ
C
         NZXOPE1 = NZXOPE
         NZXOPEQ1 = NZXOPEQ
         NZXOPT1 = NZXOPT
         NZXOPTQ1 = NZXOPTQ
      ENDIF
      NZVAR1  = NZCONF1  + NZXOPE1  + NZXOPP
      NZVARQ1 = NZCONFQ1 + NZXOPEQ1 + NZXOPPQ
C
C
C     *****************************
C     *** Reorder GPOE gradient ***
C     *****************************
C
C
C     The GPOE gradient is stored like this:
C
C     non-redundant rotations for IZ = 1
C     redundant for IZ = 1
C     non-redundant for IZ = 2
C     etc...
C
C     This makes it very difficult to do things like:
C
C     X = DNRM2(NZXOPEQ,GPOE,1)
C
C     which should be
C
C     DIMENSION GPOE(NZHOPE,*)
C     X = 0
C     DO IZ = 1, NZ
C        X = X + DNRM2(NZXOPE,GPOE(1,IZ),1)
C     END DO
C
C     ...in a million places in the XRS code.
C
C     Instead, we reorder GPOE to:
C
C     non-redundant for IZ = 1,..,NZ
C     redundant     for IZ = 1,..,NZ
C
C     The only routine where the full GPOE is used is RSIGOOGB where
C     we have to take care.
C     (Note: no reordering if no redundant rotations, e.g. SCF)
C
C
      IF (NZ .GT. 1 .AND. NZHOPE.GT.NZXOPE) THEN
         GOEER = .TRUE.
         CALL MEMGET2('REAL','GOEER',KGOEER,NZHOPEQ,WORK,KFREE,LFREE)
         CALL DCOPY(NZHOPEQ,GOEE,1,WORK(KGOEER),1)
         CALL REOGP(NZXOPE,NZHOPE,NZ,WORK(KGOEER),GOEE)
      ELSE
         GOEER = .FALSE.
      END IF
C
C
      CALL TDCBXRS(OPTRST,OPTNRM,OPTDIH,OPTSTERN,IPROPT,MOPT_MXMICRO,
     &     INTDEF_MIC,ITRINT_MIC,CNVINT_MIC,
     &     NZVAR1, NZVARQ1,NZCONF1,NZCONFQ1,
     &     NZXOPT1,NZXOPTQ1,NZXOPE1,NZXOPEQ1,NZXOPP,NZXOPPQ,
     &     NZHOPE,
     &     KZCONF,THR,OPT_CIPROGRAM,EACTIV,IOPTTB,
     &     MCTYPE,OPT_NOPFQ,OPT_NOFQX,JTRLVL)
C
C
C
      NREDM  = 3 * (MOPT_MXMICRO+2)
      NEXSIM = 1
      CALL MEMGET2('REAL','IBTYP',KIBTYP,2*NREDM     ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','IBCVC',KIBCVC,NREDM       ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','IBEVC',KIBEVC,NREDM       ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','IBPVC',KIBPVC,NREDM       ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EVALR',KEVALR,NEXSIM      ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EVECR',KEVECR,NEXSIM*NREDM,WORK,KFREE,LFREE)
C
C
C     ***************************************************
C     *** Perform micro iterations **********************
C     ***************************************************
C
C     June 2001 hjaaj: CREF is now read from LUKRMC inside XRSCTL
C     for more convenient use in lin.resp. etc.
C     (no more memory required with this change)
C
      WORK(KEVALR) = D0
      CALL XRSCTL(GCI,GOEE,GOEP,
     &        WORK(KIBTYP),WORK(KIBCVC),WORK(KIBEVC),WORK(KIBPVC),
     &        RCNV,WORK(KEVALR),WORK(KEVECR),WORK,KFREE,LFREE)
C
C     Save new orbitals etc.
C
      CALL KRMCSAV(0,CMO,CREF,WORK,KFREE,LFREE)
C
C     Restore GOEE in normal order, if REOGP has been called
C
      IF (GOEER) CALL DCOPY(NZHOPEQ,WORK(KGOEER),1,GOEE,1)
C
C     Save information for final summary output
C
      DINFO(7) = D0
      DINFO(8) = STPLEN
      DINFO(10) = D1
      DINFO(19) = STPLCI
      DINFO(20) = STPLEE
      DINFO(21) = STPLEP
C
      IF (COMPROT) THEN
C
C        Restore lengths and pointers.
C
         NZXOPE = NZXOPE_SAV
         NZHOPE = NZHOPE_SAV
         NZXOPP = NZXOPP_SAV
         NZXOPPQ = NZXOPPQ_SAV
         NZXOPEQ = NZXOPEQ_SAV
         NZHOPEQ = NZHOPEQ_SAV

         call update_orbital_rotation_indices(nzhope,
     &                                        nzxope,
     &                                        nzxopp,
     &                                        save_pp,
     &                                        save_pn)

         deallocate(save_pp)
         deallocate(save_pn)
      END IF
C
      CALL MEMREL('DIRNR.after XRSCTL',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DIRNR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tdcbxlr */
      SUBROUTINE TDCBXRS(OPTRST,OPTNRM,OPTDIH,OPTSTERN,IPROPT,MXMICRO,
     &     INTDEF_MIC,ITRINT_MIC,CNVINT_MIC,
     &     NZVAR_OPT,NZVARQ_OPT,NZCONF_OPT,NZCONFQ_OPT,
     &     NZXOPT_OPT,NZXOPTQ_OPT,NZXOPE_OPT,NZXOPEQ_OPT,
     &     NZXOPP_OPT,NZXOPPQ_OPT,
     &     NZHOPE_OPT,
     &     KZCONF_OPT,THR,CIPROGRAM,
     &     EMC,IOPTTB,MCTYPE_OPT,OPT_NOPFQ,OPT_NOFQX,JOPT_JTRLVL)
C***********************************************************************
C
C     Transfer information to response common blocks (dcbxrs.h)
C
C     Input:
C        Information from optimization module
C
C     Output:
C
C     Written by J. Thyssen - Nov 19 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
C
      LOGICAL OPTRST,OPTNRM,OPTDIH,OPTSTERN,OPT_NOPFQ,OPT_NOFQX
      DIMENSION ITRINT_MIC(2),CNVINT_MIC(2)
      DIMENSION IOPTTB(-1:2)
      CHARACTER CIPROGRAM*8
C
      CALL SETRSP
C
C     Characters...
C
      XRS_CIPROGRAM = CIPROGRAM
C
C     Logicals...
C
      RSREST = OPTRST
      STATIC = .TRUE.
      LINEQ  = .TRUE.
      TKNORM = OPTNRM
      DIAGHE = OPTDIH
      LSVCFG(1) = .TRUE.
      LSVCFG(2) = .TRUE.
      STERNH = OPTSTERN
      FKRMC  = .TRUE.
      XRS_NOPFQ = OPT_NOPFQ
      XRS_NOFQX = OPT_NOFQX
C
C     integers...
C
      JTRLVL = JOPT_JTRLVL
      IPRXRS = IPROPT
      JSYMOP = 1
      JTIMOP = 1
      JOPSY  = 1
      MAXITR = MXMICRO
      NFREQ  = 1
      NEXSIM = NFREQ
      NEXSTV = NFREQ
      NEXCNV = NFREQ
      LOFFTY = 0
      NREDM  = 3 * (MAXITR+2)
      NCRED  = 0
      NERED  = 0
      NPRED  = 0
      NZRED  = NCRED + NERED + NPRED
      N2REDM  = NREDM**2
      INTDEF = INTDEF_MIC
      INTFLG = INTDEF_MIC
      ITRINT(1) = ITRINT_MIC(1)
      ITRINT(2) = ITRINT_MIC(2)
      DO I = -1, 2
         IXRSTB(I) = IOPTTB(I)
      END DO
      NZVAR     = NZVAR_OPT
      NZVARQ    = NZVARQ_OPT
      NZCONF    = NZCONF_OPT
      NZCONFQ   = NZCONFQ_OPT
      NZXOPT    = NZXOPT_OPT
      NZXOPTQ   = NZXOPTQ_OPT
      NZXOPE    = NZXOPE_OPT
      NZXOPEQ   = NZXOPEQ_OPT
C
      NZXOPP    = NZXOPP_OPT
      NZXOPPQ   = NZXOPPQ_OPT
C
      NZHOPE    = NZHOPE_OPT
CDEBUG DEBUG
c      NZXOPP = 0
c      NZXOPPQ = 0
c      NZXOPT = 0
c      NZXOPPT = 0
CDEBUG DEBUG
c     NZCONF = 0
c     NZCONFQ = 0
C
      MCTYPE = MCTYPE_OPT
C
C     reals...
C
      THCXRS = THR
      CNVINT(1) = CNVINT_MIC(1)
      CNVINT(2) = CNVINT_MIC(2)
      RESFAC    = 1.0D3
C     with resfac = 1.01 we only take one new trial vector each micro it.
      RESFAC    = 1.01D00
      ENRGY     = EMC
C
C     pointers...
C
      KZCONF = KZCONF_OPT
C
C     units...
C
      LUXVC = 29
      LURST = 29
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtnorb */
      SUBROUTINE GTNORB(CMO,CREF,IBTYP,IBCVC,IBEVC,IBPVC,
     &                  EVALR,EVECR,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate new orbitals ( CMO_{new} = exp ( -\kappa) CMO_{old} )
C
C     Input:
C        CMO - old coefficients
C        CREF - ref. CI vector
C
C     Output:
C        CMO - new coefficients
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      PARAMETER ( D1 = 1.0D00 )
      PARAMETER ( THRTT = 1.0D-4 )
#include "thrzer.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
C
      DIMENSION CMO(*), WORK(*), IBTYP(2,*), IBCVC(*),
     &          IBEVC(*),IBPVC(*), EVECR(*)
      DIMENSION CREF(*)
C
      CALL QENTER('GTNORB')
      KFRSAV = KFREE
C
C     ***************************************************
C     *** Get solution vectors = kappa ******************
C     ***************************************************
C
      CALL MEMGET2('REAL','KAPPA',KKAPPA,NZXOPTQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DELTACI',KDELTACI,NZCONFQ,WORK,KFREE,LFREE)
C
CDEBUG - jth
CDEBUG DZERO - can safely be removed as XRSXV1 does a DZERO
      CALL DZERO(WORK(KKAPPA),NZXOPTQ)
C
      CALL KRMCVEC(IBTYP,IBCVC,IBEVC,IBPVC,
     &     EVALR,EVECR,WORK(KKAPPA),WORK(KDELTACI),
     &     WORK(KFREE),LFREE)
C
      IF (IPROPT .GE. 5 .AND. NZXOPT .GT. 0) THEN
         WRITE(LUPRI,'(A)') ' (GTNORB) Kappa matrix'
         CALL RPRKAP(WORK(KKAPPA),WORK(KKAPPA+NZXOPEQ),
     &               get_orbital_rotation_indices_pp(),NZXOPE,
     &               get_orbital_rotation_indices_pn(),NZXOPP,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
      IF (IPROPT .GE. 10 .AND. NZCONF .GT. 0) THEN
         WRITE(LUPRI,'(A)') ' (GTNORB) CI delta vector'
         CALL RPRCI(WORK(KDELTACI),NZCONF,NZ,
     &        OPT_CIPROGRAM,THRPCI,LUPRI)
      END IF
C
C     *************************************
C     *** Calculate new CI coefficients ***
C     *************************************
C
      IF ( NZCONF .GT. 0 ) THEN
C
         CALL WRTKRMC(LUKRM2,'CREFOLD ',CREF,NZCONFQ)
C
C        Add delta:
C
         CALL DAXPY(NZCONFQ,D1,WORK(KDELTACI),1,CREF,1)
C
         IF (IPROPT .GE. 20 ) THEN
            WRITE(LUPRI,'(A)') ' (GTNORB) New unnorm. CI vector'
            CALL RPRCI(CREF,NZCONF,NZ,
     &           OPT_CIPROGRAM,THRZER,LUPRI)
         END IF
C
         CALL RCINORM(CREF,10)
C
         IF (IPROPT .GE. 20 ) THEN
            WRITE(LUPRI,'(A)') ' (GTNORB) New norm. CI vector'
            CALL RPRCI(CREF,NZCONF,NZ,
     &           OPT_CIPROGRAM,THRZER,LUPRI)
         END IF
C
      END IF
C
C
C     ***************************************************
C     *** Rotate old coefficients ***********************
C     ***************************************************
C
C     Unpack KAPPA
C
      CALL MEMGET2('REAL','UKAPPA',KUKAPPA,N2ORBXQ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KUKAPPA),N2ORBXQ)
      CALL RUPKWOP(WORK(KKAPPA),WORK(KKAPPA+NZXOPEQ),
     &             WORK(KUKAPPA),NORBT,
     &             get_orbital_rotation_indices_pp(),NZXOPE,
     &             get_orbital_rotation_indices_pn(),NZXOPP,NZ)
      IF (IPROPT .GE. 10) THEN
         WRITE(LUPRI,'(A)') ' (GTNORB) Unpacked kappa'
         CALL PRQMAT(WORK(KUKAPPA),NORBT,NORBT,NORBT,NORBT,NZ,
     &               IPQTOQ(1,0),LUPRI)
      END IF
!
C     Save old coefficients
      CALL WRTKRMC(LUKRM2,'OLDORB  ',CMO,NCMOTQ)
C
C     Get new coefficients
C
      IOFF = 0
      DO I = 1,NFSYM
         IF (NOCC(I) .GT. 0) THEN
CJTH     IF (NORB(I) .GT. 0) THEN
            CALL ROTCMO(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,
     &                  WORK(KUKAPPA+IOFF),
     &                  NORBT,IPQTOQ(1,0),IPROPT,WORK(KFREE),LFREE)
         END IF
         IOFF = IOFF + (NORBT+1)*NORB(I)
      END DO
C
C     Reorthonomalize coefficients:
C
C     inside RGMCMO: order from positron,electron to
C                    electron,positron in CMO changed
C                    for the Gram-Schmidt orthonormalization,
C                    because we want to start with the
C                    occupied orbitals
C
      CALL RGMCMO(CMO,WORK(KFREE),LFREE)
C
C     Write new coefficients
C
      CALL WRTKRMC(LUKRM2,'NEWORB  ',CMO,NCMOTQ)
      CALL WRTKRMC(LUKRM2,'CREF    ',CREF,NZCONFQ)
C
      CALL MEMREL('GTNORB',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('GTNORB')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirstp */
      SUBROUTINE DIRSTP(ISTEP,CMO,CREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Step control
C
C     Input:
C        ISTEP: control parameter
C               istep =  0: normal MC step check
C               istep = -1: restart step (with valid DEPRED)
C               istep = -2: MC converged, generate info for analysis
C        CMO  : current MO orbitals
C
C     Output:
C        ISTEP: return parameter
C               istep = -1: no step check (close to convergence)
C               istep =  0: step is OK
C               istep =  1: step is too large
C               istep =  2: step is too large, backstep NOT possible
C
C     Written by J. Thyssen - Nov 27 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (THDE = 1.0D-10, DP1 = 0.1 D00)
#include "consts.h"
C
#include "dcbopt.h"
#include "dcbkrmc_itinfo.h"
C
      DIMENSION CMO(*), CREF(*), WORK(*)
      LOGICAL   NOREJT
C
      CALL QENTER('DIRSTP')
C
C     We are not allowed to reject restart steps (ISTEP = -1).
C
      NOREJT = .FALSE.
      IF (ISTEP .LT. 0) THEN
         NOREJT = .TRUE.
      END IF
C
C     Calculate actual energy difference, and compare to predicted.
C
      ISTEP = 0
      DEACT = EMCSCF - EMCOLD
      IF ( ABS(DEPRED) .GT. THDE ) THEN
         RATIO = DEACT / DEPRED
         WRITE(LUPRI,9000) DEACT,DEPRED,RATIO
      ELSE
         RATIO = D1
         ISTEP = -1
         WRITE(LUPRI,9002) DEACT,DEPRED
      END IF
 9000 FORMAT(/' (DIRSTP): Energy difference;',
     &       /' actual, predicted, and ratio: ',1P,2D15.5,0P,F12.6)
 9002 FORMAT(/' (DIRSTP): Close to convergence, ratio set to one.'
     &       /' Energy difference; actual and predicted: ',1P,2D15.5)
C
C     Save info for final summary output.
C
      DINFO(4) = DEPRED
      DINFO(5) = DEACT
      DINFO(6) = RATIO
C
C     Normal MCSCF optimization
C
      RATVGD = D1 - DP1 * (D1 * RATGOD)
      RTSAVE = RTRUST
      IF ( IOPT_STATE .EQ. 1 ) THEN
         IF ( RATIO .LT. RATMIN ) THEN
            RTRUST = STPRED * RTRUST
            IF ( RATIO .LT. RATREJ .AND. .NOT. NOREJT ) THEN
               ISTEP = 1
               RTRUST = MIN(RTRUST, STPRED * STPLEN)
               IF (RATIO .LT.  0.0D0) RTRUST = STPRED * RTRUST
               IF (RATIO .LT. -1.0D0) RTRUST = STPRED * RTRUST
            END IF
         ELSE IF (RATIO .GT. RATGOD) THEN
            RTRUST = MIN (STPMAX, STPINC * RTRUST)
            IF (RATIO .GT. RATVGD) RTRUST = MIN(STPMAX, STPINC * RTRUST)
         END IF
      ELSE
         CALL QUIT('*** ERROR in DIRSTP ***: '
     &        //'excited state step control not implemented')
      END IF
C
C     Check if step is accepted
C
C     debug
c     write(6,*) ' (dirstp) itmac, itbck ',itmac,itbck
c     IF (ITMAC .EQ. 2 .AND. ITBCK .EQ. 0) THEN
c        ISTEP = 1
c        RTRUST = 0.5 * RTRUST
c     END IF
c     write(6,*) ' (dirstp) istep ',istep
c     write(6,*) ' (dirstp) donr,ldirnr ',donr,ldirnr
      IF (ISTEP .LE. 0) THEN
         EMCOLD = EMCSCF
         GOTO 1000
      END IF
C
C     ************************
C     *** Step is rejected ***
C     ************************
C
      IF (DONR) THEN
C
C        Backstep for pure Newton-Raphson cannot be implemented
C
         IF ( DEACT .LE. D0 ) THEN
C
C           However, if the energy is decreasing we accept the step
C           anyway.
C
            WRITE(LUPRI,'(/A,F12.6,A)')
     &           ' *** WARNING *** bad ratio (',RATIO,
     &           '), but step accepted because DEACT < 0'
            ISTEP = 0
            EMCOLD = EMCSCF
         ELSE
            ISTEP = 2
         END IF
         GOTO 1000
      END IF
C
      IF (LDIRNR) THEN
C
C        Backstep for NR restricted step algorithm
C
         CALL RNRBSTP(CMO,CREF,WORK,KFREE,LFREE)
C
      ELSE
C
C        Backstep for NEO restricted step algorithm
C
         CALL QUIT('*** NOT IMPLEMENTED ***')
      END IF
C
C
 1000 CONTINUE
C
C     Save information for final summary outout
C
      IF (IPROPT .GE. 2)
     &     WRITE(LUPRI,'(A,2F15.10)')
     &     ' (DIRSTP) Old and new KR-MCSCF trust radius: ',RTSAVE,RTRUST
C
      CALL QEXIT('DIRSTP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmc_getfck */
      SUBROUTINE KRMC_GETFCK(FMO,IPRINT)
C***********************************************************************
C
C     Read Fock matrix from file KRMC FOCK
C
C     Input:
C        IPRINT - print level
C        CRDFO_MAT (in dcbopt.h)
C        
C
C     Output:
C        FMO    - Fock matrix in MO basis
C        ECORE_LR (in dcbopt.h, only if CRDFO_MAT is true)
C
C     Written by J. Thyssen - Nov 27 1998; revised Aug 2015 hjaaj
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C dgroup.h : IPTQOQ
C dcborb.h : NORBT,...
C dcbopt.h : CRDFO_MAT, ECORE_LR
C dcbgen.h : LUKRM3
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbgen.h"
C
      REAL*8   FMO(*)
C
      IF( CRDFO_MAT) THEN
        CALL REAKRMC(LUKRM3,'FCMO    ',FMO,N2ORBXQ)
        CALL REAKRMC(LUKRM3,'ECORE_LR',ECORE_LR,1)
      ELSE
        CALL REAKRMC(LUKRM3,'FC      ',FMO,N2ORBXQ)
      END IF
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('FCMO from KRMC_GETFCK',-1)
         CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsnrprd */
      SUBROUTINE RSNRPRD(REDGP,REDE,EVECR,NREDM,IBTYP,NZRED)
C***********************************************************************
C
C     Calculate second order change in energy
C
C     Written by J. Thyssen - Nov 27 1998
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (DP5 = 0.50D00, D0 = 0.00D00, DP25 = 0.25D00)
C
#include "dcbopt.h"
#include "dcbibn.h"
C
      DIMENSION REDGP(NREDM),REDE(NREDM,NREDM),EVECR(NREDM)
      DIMENSION IBTYP(2,*)
C
      IF ( IPROPT .GE. 5 ) THEN
         WRITE(LUPRI,'(A)')
     &        ' (RSNRPRD) Reduced solution vector '
         CALL OUTPUT(EVECR,1,NZRED,1,1,NREDM,1,1,LUPRI)
C
         WRITE(LUPRI,'(A)')
     &        ' (RSNRPRD) Reduced gradient'
         CALL OUTPUT(REDGP,1,NZRED,1,1,NREDM,1,1,LUPRI)
         WRITE(LUPRI,'(A)')
     &        ' (RSNRPRD) Reduced Hessian'
         CALL OUTPUT(REDE,1,NZRED,1,NZRED,NREDM,NREDM,1,LUPRI)
      END IF
C
C     Calculate step length
C     ---------------------
C
C     Calculate total step length.
C     We multiply with as the trial vectors are normalized to 0.5
C
      STPLEN = DP5*DNRM2(NREDM,EVECR,1)
      STPLCI = D0
      STPLEE = D0
      STPLEP = D0
      DO I = 1, NZRED
         IF ( IBTYP(1,I) .EQ. JBCNDX ) THEN
            STPLCI = STPLCI + EVECR(I) ** 2
         ELSE IF ( IBTYP(1,I) .EQ. JBENDX ) THEN
            STPLEE = STPLEE + EVECR(I) ** 2
         ELSE
            STPLEP = STPLEP + EVECR(I) ** 2
         END IF
      END DO
      STPLCI = DP5 * SQRT(STPLCI)
      STPLEE = DP5 * SQRT(STPLEE)
      STPLEP = DP5 * SQRT(STPLEP)
C
C     Calculate predicted energy change
C     ---------------------------------
C
C     DEPRED = E - E0 = x^t g + 1/2 x^t H x
C
C     where x is the step vector, g is the gradient (note that REDGP is
C     equal to minus the gradient), and H is the Hessian.
C
      SUM = D0
      DO I = 1,NZRED
         SUM1 = D0
         DO J = 1,NZRED
            SUM1 = SUM1 + REDE(I,J)*EVECR(J)
         END DO
         SUM = SUM + EVECR(I)*(-REDGP(I) + DP5 * SUM1)
      END DO
C
C     Multiply with 0.5 as trial vectors are normalized to 0.5
C
      DEPRED = DP5*SUM
C
      IF (IPROPT .GE. 2) WRITE(LUPRI,'(A,F15.8)')
     &     ' (RSNRPRD) Predicted energy change = ',DEPRED
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmcsol */
!     SUBROUTINE KRMCSOL(NZRED,NREDM,REDE,REDS,EVALR,REDGP,
      SUBROUTINE KRMCSOL_old(NZRED,NREDM,REDE,REDS,EVALR,REDGP,
     &        EVECR,IBTYP,NPRED,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C
C     Written by J. Thyssen - Nov 22 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "consts.h"
C
      PARAMETER (DLVLSTEP = 0.2D0)
C
#include "dgroup.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbkrmc_itinfo.h"
C
      DIMENSION REDE(NREDM,NREDM), REDS(NREDM,NREDM)
      DIMENSION REDGP(NREDM), EVECR(NREDM)
      DIMENSION IBTYP(*), INERTIA(3)
      DIMENSION WORK(*)
C
      IF (IOPT_STATE .NE. 1)
     &     CALL QUIT('KMRCSOL: excited state opt. not implemented yet')
C
      RTRMAX = RTRUST*RTTOL
      RTRMIN = RTRUST/RTTOL
C
      CALL DUNIT(REDS,NREDM)
      DLVLMAX = D0
      DLVLMIN = D0
      DLVLSH  = D0
      ITER = 0
C
C
C     Start bisection iterations to find level-shift
C     ----------------------------------------------
C
 100  CONTINUE
      CALL XRSLINRED(NZRED,NREDM,REDE,REDS,DLVLSH,REDGP,
     &     EVECR,IBTYP,INERTIA,IINFO(6),IPROPT,WORK,KFREE,LFREE)
      ITER = ITER + 1
      STPL = DP5 * DNRM2(NREDM,EVECR,1)
      INONPOS = INERTIA(2) + INERTIA(3)
C
C     Save electronic Hessian's index for zero level shift
C     (used in RCNVTHR below)
C
      IF (ITER .EQ. 1) JHESSIAN = INONPOS - NPRED
C
C     Don't iterate if Newton-Raphson is requested.
C
      IF (DONR) GOTO 200
C
      IF (INONPOS .GT. NPRED .OR. STPL .GE. RTRMAX) THEN
         DLVLMAX = DLVLSH
         IF (DLVLMIN .EQ. D0) THEN
            DLVLSH = DLVLSH - DLVLSTEP
         ELSE
            DLVLSH  = DP5*(DLVLMIN+DLVLMAX)
         END IF
#ifdef UNDEF
         write(6,*) 'krmcsol: inonpos, steplen (>?), new level shift',
     &        inonpos,STPL,DLVLSH
#endif
         GOTO 100
      ELSE IF (STPL .LE. RTRMIN .AND. DLVLSH .NE. D0) THEN
         DLVLMIN = DLVLSH
         DLVLSH  = DP5*(DLVLMIN+DLVLMAX)
#ifdef UNDEF
         write(6,*) 'krmcsol: steplen<, new level shift',STPL,DLVLSH
#endif
         GOTO 100
#ifdef UNDEF
      ELSE
         write(6,*) 'krmcsol: steplen accepted, level shift',
     &        STPL,DLVLSH
#endif
      END IF
C
C     Iterations finished.
C
 200  CONTINUE
      RBETA = STPL
      EVALR = DLVLSH
      IPAR = INERTIA(3)
      IF (IPRINT .GE. 2) THEN
         WRITE(LUPRI,'(/A,I3,A,F15.8)')
     &        ' (KRMCSOL) Level shift (after',ITER,' bi-sec. iter.) = ',
     &        DLVLSH
      END IF
      IF(IPRINT.GE.4.OR.IPAR.GT.0) THEN
         IF(IPRINT.GE.3) THEN
            WRITE(LUPRI,'(A,F12.5,A)')
     &           ' (KRMCSOL) Solution vector'
            CALL PRIRVC(EVECR,NZRED,IBTYP)
         ENDIF
         WRITE(LUPRI,'(/3X,A,F12.5)')
     &      '* Reduced matrix. Level shift:',DLVLSH
         WRITE(LUPRI,'(5X,A,I5)')
     &        '* Number of positive eigenvalues:',INERTIA(1),
     &        '* Number of zero     eigenvalues:',INERTIA(3),
     &        '* Number of negative eigenvalues:',INERTIA(2)
      ENDIF
C
      DINFO(22) = DLVLSH
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmcsol */
!     SUBROUTINE KRMCSOL_new_to_come(NZRED,NREDM,REDE,REDS,EVALR,REDGP,
      SUBROUTINE KRMCSOL(NZRED,NREDM,REDE,REDS,EVALR,REDGP,
     &        EVECR,IBTYP,NPRED,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C
C     Written by J. Thyssen - Nov 22 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "consts.h"
C
      PARAMETER (DLVLSTEP = 0.2D0)
C
#include "dgroup.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbkrmc_itinfo.h"
#include "dcbibn.h"
C
      DIMENSION REDE(NREDM,NREDM), REDS(NREDM,NREDM)
      DIMENSION REDGP(NREDM), EVECR(NREDM)
      DIMENSION IBTYP(2,*), INERTIA(3)
      DIMENSION WORK(*)
      PARAMETER (STPL_EP_MAX = 0.0100D0)
C
      IF (IOPT_STATE .NE. 1)
     &     CALL QUIT('KMRCSOL: excited state opt. not implemented yet')
C
      RTRMAX = RTRUST*RTTOL
      RTRMIN = RTRUST/RTTOL
      STPL_EP_MIN = STPL_EP_MAX/RTTOL
C
      IF (IPRINT .GE. 5) THEN
         WRITE(LUPRI,'(/A,2I5,L10/A,3F15.10)')
     &   ' KRMCSOL entry. NZRED, NPRED, DONR:',NZRED,NPRED,DONR,
     &   '             RTTOL, RTRMAX, RTRMIN:',RTTOL,RTRMAX,RTRMIN
      END IF
C
      CALL DUNIT(REDS,NREDM)
      DLVLMAX = D0
      DLVLMIN = D0
      DLVLSH  = D0
      INCREASE_SHIFT = 0 ! 0=unknown, >0=increase, <0=decrease
      ITER = 0
      STPL = D0
C
C
C     Start bisection iterations to find level-shift
C     ----------------------------------------------
C
 100  CONTINUE
      STPL_OLD = STPL
      CALL XRSLINRED(NZRED,NREDM,REDE,REDS,DLVLSH,REDGP,
     &     EVECR,IBTYP,INERTIA,IINFO(6),IPROPT,WORK,KFREE,LFREE)
      INONPOS = INERTIA(2) + INERTIA(3)
      ITER = ITER + 1
      STPL = DP5 * DNRM2(NZRED,EVECR,1)
      STPL_CI = D0
      STPL_EE = D0
      STPL_EP = D0
      IF (NPRED .GT. 0 .OR. IPRINT .GE. 5) THEN
         DO I = 1, NZRED
            IF (IBTYP(1,I) .EQ. JBCNDX) THEN
               STPL_CI = STPL_CI + EVECR(I)**2
            ELSE IF (IBTYP(1,I) .EQ. JBENDX) THEN
               STPL_EE = STPL_EE + EVECR(I)**2
            ELSE IF (IBTYP(1,I) .EQ. JBPNDX) THEN
               STPL_EP = STPL_EP + EVECR(I)**2
            ELSE
               CALL QUIT('unknown IBTYP')
            END IF
         END DO
         STPL_CI = DP5 * SQRT(STPL_CI)
         STPL_EE = DP5 * SQRT(STPL_EE)
         STPL_EP = DP5 * SQRT(STPL_EP)
         IF (IPRINT .GE. 5)
     &      WRITE (LUPRI,'(A,I5,1P,G10.2,I5/A,3G10.2)')
     &      'KRMSOL ITER, STPL, INONPOS',ITER,STPL,INONPOS,
     &      'step length CI, e-e, e-p',STPL_CI,STPL_EE,STPL_EP
      END IF
C
C     Save electronic Hessian's index for zero level shift
C     (used in RCNVTHR below)
C
      IF (ITER .EQ. 1) THEN
         STPL_NR = STPL
         JHESSIAN = INONPOS - NPRED
         IF (IPRINT .GE. 5) WRITE(LUPRI,'(A,I5)')
     &      ' KRMCSOL Newton-Raphson. JHESSIAN',JHESSIAN
C
         IF (DONR) GOTO 200 ! Do not iterate if Newton-Raphson is requested.
C
         IF (STPL_NR .GT. RTRMAX .OR. STPL_EP .GT. STPL_EP_MAX
     &       .OR. JHESSIAN .NE. 0) THEN
            IF (JHESSIAN .GT. 0 .OR. NPRED .EQ. 0) THEN
               INCREASE_SHIFT = -1
            ELSE IF (JHESSIAN .LT. 0) THEN
               INCREASE_SHIFT = 1
            ELSE ! JHESSIAN .eq. 0 .and. NPRED .gt. 0
               INCREASE_SHIFT = 1 ! TODO not always correct!!!! (I think)
            END IF
         ELSE
            GO TO 200 ! Newton-Raphson is optimal
         END IF

      END IF
C
      IF (STPL .GT. RTRMAX .OR. STPL_EP .GT. STPL_EP_MAX .OR.    ! too long step
     &    INONPOS .NE. NPRED ) THEN                              ! wrong number of negative eigenvalues
         IF (ITER .LE. 2 .AND. JHESSIAN .EQ. 0 .AND.
     &       STPL .GT. STPL_OLD .AND. STPL_OLD .GT. RTRMAX) THEN ! reset
            INCREASE_SHIFT = -INCREASE_SHIFT
            DLVLMAX = D0
            DLVLMIN = D0
         END IF
         IF (INCREASE_SHIFT .GT. 0) THEN
            DLVLMIN = DLVLSH
            IF (DLVLMAX .NE. D0) THEN
               DLVLSH  = DP5*(DLVLMIN+DLVLMAX)
            ELSE
               DLVLSH  = DLVLSH + DLVLSTEP
            END IF
         ELSE IF (INCREASE_SHIFT .LT. 0) THEN
            DLVLMAX = DLVLSH
            IF (DLVLMIN .NE. D0) THEN
               DLVLSH  = DP5*(DLVLMIN+DLVLMAX)
            ELSE
               DLVLSH  = DLVLSH - DLVLSTEP
            END IF
         ELSE
            call quit('KRMCSOL: INCREASE_SHIFT=0 not programmed yet')
         END IF
         IF (IPRINT .GE. 5) write(LUPRI,'(A,1P,2G20.10)')
     &      'krmcsol: steplen (>?), new level shift', STPL,DLVLSH
      ELSE IF ( DLVLSH .EQ. 0.0D0 .OR.                                 ! Newton step
     &          (STPL .GE. RTRMIN .AND. INONPOS .EQ. NPRED) .OR.       ! restricted step too short
     &           STPL_EP .GE. STPL_EP_MIN )                            ! e-p in allowed range
     &   THEN ! step accepted
         IF (IPRINT .GE. 5) write(LUPRI,'(A,1P,2G20.10)')
     &      'krmcsol: steplen accepted, level shift', STPL,DLVLSH
         GOTO 200
      ELSE
         IF (INCREASE_SHIFT .GT. 0) THEN
            DLVLMAX = DLVLSH
            DLVLSH  = DP5*(DLVLMIN+DLVLMAX)
         ELSE IF (INCREASE_SHIFT .LT. 0) THEN
            DLVLMIN = DLVLSH
            DLVLSH  = DP5*(DLVLMIN+DLVLMAX)
         ELSE
            call quit('KRMCSOL b: INCREASE_SHIFT=0 not programmed yet')
         END IF
         IF (IPRINT .GE. 5) write(LUPRI,'(A,I5,1P,2G20.10)')
     &      'krmcsol: inonpos, steplen (<), new level shift',
     &      inonpos,STPL,DLVLSH
      END IF
      IF (ITER .GT. 100) CALL QUIT('KRMSOL: iterations not converged')
      GOTO 100
C
C     Iterations finished.
C
 200  CONTINUE
      RBETA = STPL
      EVALR = DLVLSH
      IPAR = INERTIA(3)
      IF (IPRINT .GE. 2) THEN
         WRITE(LUPRI,'(/A,I3,A,F15.8)')
     &        ' (KRMCSOL) Level shift (after',ITER,' bi-sec. iter.) = ',
     &        DLVLSH
      END IF
      IF(IPRINT.GE.4.OR.IPAR.GT.0) THEN
         IF(IPRINT.GE.3) THEN
            WRITE(LUPRI,'(A,F12.5,A)')
     &           ' (KRMCSOL) Solution vector'
            CALL PRIRVC(EVECR,NZRED,IBTYP)
         ENDIF
         WRITE(LUPRI,'(/3X,A,F12.5)')
     &      '* Reduced matrix. Level shift:',DLVLSH
         WRITE(LUPRI,'(5X,A,I5)')
     &        '* Number of positive eigenvalues:',INERTIA(1),
     &        '* Number of zero     eigenvalues:',INERTIA(3),
     &        '* Number of negative eigenvalues:',INERTIA(2)
      ENDIF
C
      DINFO(22) = DLVLSH
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rcnvthr */
      FUNCTION RCNVTHR()
C***********************************************************************
C
C     Calculate threshold for micro iterations.
C
C     Written by J. Thyssen - Dec 12 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "consts.h"
      PARAMETER (THQCXX = 0.10D0, RTTHRG = 4.00D0)
      PARAMETER (THGFAC = 0.98D0, THGFC2 = 0.70D0)
C
#include "dcbopt.h"
C
C     Heuristics ripped from Dalton (CNVTHR, sirneo.F)
C
C     Note that the Dirac norm is 0.5 Dalton norm, so all norms
C     have been replaced with 2 * norm.
C
C
C     If new step (RBETA) is greater than previous step (STPLEN)
C     or the level shift is .NE. 0 we are not in the local region,
C     i.e. use global threshold:
C
      IF ( (STPLEN .GT. D0 .AND. RBETA .GT. STPLEN) .OR.
     &     ( DLVLSH .NE. D0 ) )THEN
         IF (JHESSIAN .NE. 0) THEN
C
C           Hessian has incorrect structure, and it's not worth
C           the effort to converge better than 50% (we think!)
C
C           (TODO FIXME: put DP5 = 50% in common block/input)
C
            THRQN = DP5 * D2 * GNORM(5)
         ELSE
C
C           Hessian has correct structure, we converge a little
C           better (by default: THQLIN = 20%)
C
            THRQN = THQLIN * D2 * GNORM(5)
         END IF
      ELSE
         GCNRM = D2 * GNORM(1)
         GONRM = D2 * GNORM(4)
         THRQN = THQKVA * (GCNRM * GCNRM)
         IF (GCNRM .GE. 5.0D-3) THRQN = THQCXX * THRQN
         THRQN = MAX( THRQN, THQKVA*(GONRM*GONRM) )
         THRQN = MIN( THRQN, THQLIN* D2 * GNORM(5) )
         IF (THRQN .LE. OPT_THRGRD) THEN
            THRQN = THGFAC*D2*OPT_THRGRD
         ELSE IF (THRQN .LE. RTTHRG*OPT_THRGRD) THEN
            THRQN = THGFC2*D2*OPT_THRGRD
         END IF
      END IF
      RCNVTHR = DP5 * THRQN
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rcinorm */
      SUBROUTINE RCINORM(CREF,IPRINT)
C***********************************************************************
C
C     Normalize CI vector.
C
C     Written by J. Thyssen - Dec 14 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
      PARAMETER ( THRTT = 1.0D-4 )
C
#include "dcbopt.h"
#include "dgroup.h"
C
      DIMENSION CREF(NZCONF,*)
C
      LOGICAL CNRM
C
      CALL QENTER('RCINORM')
C
      MZ = MIN(2,NZ)
      CNRM = .FALSE.
      DN = DNORM2(NZCONFQ,CREF,1)
      IF ( DN .NE. D1 ) THEN
         IF ( DN .LE. THRTT ) THEN
            DN = D1 / DN
            IF (IPRINT .GE. 5) WRITE(LUPRI,'(/A,F20.12)')
     &        ' (RCINORM) CI vector renormalized with factor ',DN
            CALL DSCAL(NZCONFQ, DN, CREF, 1)
            DN = DNORM2(NZCONFQ,CREF,1)
         END IF
         DN = D1 / DN
         IF (IPRINT .GE. 5) WRITE(LUPRI,'(/A,F20.12)')
     &        ' (RCINORM) CI vector renormalized with factor ',DN
         CALL DSCAL(NZCONFQ, DN, CREF, 1)
         CNRM = .TRUE.
      END IF
      IF ( IPRINT .GE. 20 .AND. CNRM) THEN
         WRITE(LUPRI,'(/A)') ' (RCINORM) Renormalized CI vector'
         CALL RPRCI(CREF,NZCONF,MZ,OPT_CIPROGRAM, THRZER,LUPRI)
      END IF
C
      CALL QEXIT('RCINORM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rtktime */
      SUBROUTINE RTKTIME(TFLAG,N)
C***********************************************************************
C
C     Start/stop timer no. N
C
C     Used for detailed timing statistics.
C
C     Written by J. Thyssen - Dec 22 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
C
      LOGICAL TFLAG
C
      LOGICAL STATET(NOPTTIM)
      DIMENSION CPU(NOPTTIM), WLL(NOPTTIM)
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      SAVE FIRST, STATET, CPU, WLL
C
C     STATET (state of time; stopped (.false.) or running (.true.)
C
      IF (FIRST) THEN
         DO I = 1, NOPTTIM
            STATET(I) = .FALSE.
         END DO
         FIRST = .FALSE.
      END IF
C
      IF ( N .LT. 1 .OR. N .GT. NOPTTIM ) THEN
         WRITE(LUPRI,'(A,I3)')
     &        ' *** ERROR in RTKTIME ***: parameter N invalid; N = ',N
         CALL QUIT('*** ERROR in RTKTIME ***')
      END IF
      IF ( STATET(N) .EQV. TFLAG ) THEN
         WRITE(LUPRI,'(A)')
     &        ' *** ERROR in RTKTIME ***: state of timer = tflag'
         CALL QUIT('*** ERROR in RTKTIME ***')
      END IF
C
      STATET(N) = TFLAG
      CALL GETTIM(CPU1,WLL1)
C
      IF (TFLAG) THEN
C
C        Start timer no. N
C
         CPU(N) = CPU1
         WLL(N) = WLL1
      ELSE
C
C        Stop timer no. N and add results up.
C
         CPUOPT(N) = CPUOPT(N) + (CPU1 - CPU(N))
         WLLOPT(N) = WLLOPT(N) + (WLL1 - WLL(N))
C
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rnrbstp */
      SUBROUTINE RNRBSTP(CMO,CREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Perform backstep (NR restricted step algorihtm)
C
C     Written by J. Thyssen - Dec 28 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgen.h"
C
      DIMENSION WORK(*)
      DIMENSION CMO(*), CREF(*)
C
      LOGICAL FNDLAB
C
      CALL QENTER('RNRBSTP1')
      KFRSAV = KFREE
C
C     Strategy:
C     (1) read old reduced L (NEO) or H (NR) matrix from prev. macro
C     (2) solve for reduced trust radius
C     (3) calculate new CI ref. vector, orbitals, predictions etc
C
C
C     (1) read reduced L or H from from previous macro
C     ------------------------------------------------
C
C     ...and tons of other information.
C
C     Get old orbitals
C
      CALL REAKRMC(LUKRMC,'OLDORB  ',CMO,NCMOTQ)
      CALL REAKRMC(LUKRMC,'OLDCREF ',CREF,NZCONFQ)
C
C     Get L and other stuff
C
      REWIND(LUKRMC)
      IF ( FNDLAB('LREDUCED',LUKRMC) ) THEN
         READ(LUKRMC) NREDM
      ELSE
         WRITE(LUPRI,'(/A)')
     &        ' *** ERROR in RNRBSTP *** LREDUCED not found on LUKRMC'
         CALL QUIT('*** ERROR in RNRBSTP ***'//
     &        ' LREDUCED not found on LUKRMC')
      END IF
      REWIND(LUKRMC)
      CALL MOLLAB('LREDUCED',LUKRMC,LUPRI)
C
      CALL MEMGET2('REAL','REDE' ,KREDE,NREDM*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','REDS' ,KREDS,NREDM*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','REDGP',KREDGP,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EVECR',KEVECR,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBTYP',KIBTYP,2*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBEVC',KIBEVC,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBPVC',KIBPVC,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBCVC',KIBCVC,NREDM,WORK,KFREE,LFREE)
C
      CALL REALRED(LUKRMC,NZRED,NREDM,WORK(KREDE),WORK(KREDGP),
     &     WORK(KEVECR),WORK(KIBTYP),
     &     WORK(KIBEVC),WORK(KIBPVC),WORK(KIBCVC),NPRED)
C
      IF (IPROPT .GE. 5) THEN
         WRITE(LUPRI,'(/A,3I4)')
     &        ' (RNRBSTP1) NZRED, NREDM, NPRED ',NZRED, NREDM, NPRED
         WRITE(LUPRI,'(A)') ' (RNRBSTP1) Reduced Hessian'
         CALL OUTPUT(WORK(KREDE),1,NZRED,1,NZRED,NREDM,NREDM,1,LUPRI)
         WRITE(LUPRI,'(A)') ' (RNRBSTP1) Minus reduced gradient'
         CALL OUTPUT(WORK(KREDGP),1,NZRED,1,1,NREDM,1,1,LUPRI)
         WRITE(LUPRI,'(A)') ' (RNRBSTP1) IBTYP'
         DO I = 1, NZRED
            WRITE(LUPRI,'(5X,2I4)')
     &           IGETINT(WORK(KIBTYP),1 + (I-1)*2),
     &           IGETINT(WORK(KIBTYP),2 + (I-1)*2)
         END DO
      END IF
C
C     (2) solve for reduced trust radius
C     ----------------------------------
C
C     Trust radius has been reduced in DIRSTP.
C
      CALL KRMCSOL(NZRED,NREDM,WORK(KREDE),WORK(KREDS),
     &     EVALR,WORK(KREDGP),WORK(KEVECR),
     &     WORK(KIBTYP),NPRED,IPROPT,WORK,KFREE,LFREE)
C
      WRITE(LUPRI,'(/A/4X,2F15.8)')
     &     ' (RNRBSTP1) Backstep; new trust radius & new levelshift : ',
     &     RTRUST,DLVLSH
C
C     Write stuff back to KRMC_REDL
C
      REWIND(LUKRM4)
      CALL WRTLRED(LUKRM4,NZRED,NREDM,WORK(KREDE),WORK(KREDGP),
     &     WORK(KEVECR),WORK(KIBTYP),
     &     WORK(KIBEVC),WORK(KIBPVC),WORK(KIBCVC),NPRED)
C
C     Calculate new predicted energy change
C
      CALL RSNRPRD(WORK(KREDGP),WORK(KREDE),WORK(KEVECR),
     &     NREDM,WORK(KIBTYP),NZRED)
C
C     Save on KRMCSCF and find new orbitals and CI ref. vector
C
      CALL MEMREL('RNRBSTP1',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL KRMCSAV(0,CMO,CREF,WORK,KFREE,LFREE)
C
      CALL QEXIT('RNRBSTP1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrtlred */
      SUBROUTINE WRTLRED(LUNIT,NZRED,NREDM,REDE,REDGP,EVECR,
     &     IBTYP,IBEVC,IBPVC,IBCVC,NPRED)
C***********************************************************************
C
C     Dump info needed for backsteps on KRMC_LRED/KRMCSCF
C
C     Written by J. Thyssen - Dec 28 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
C
      DIMENSION REDE(*), REDGP(*)
      DIMENSION IBTYP(2,*), IBEVC(*), IBPVC(*), IBCVC(*)
      DIMENSION EVECR(*)
C
      CALL QENTER('WRTLRED')
C
C     Dump info
C
      WRITE(LUNIT) NREDM,
     &     NZRED,
     &     (REDE(I),I=1,NREDM*NREDM),
     &     (REDGP(I),I=1,NREDM),
     &     (EVECR(I),I=1,NREDM),
     &     (IBTYP(1,I),IBTYP(2,I),I=1,NREDM),
     &     (IBEVC(I),I=1,NREDM),
     &     (IBPVC(I),I=1,NREDM),
     &     (IBCVC(I),I=1,NREDM),
     &     NPRED
C
      CALL QEXIT('WRTLRED')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck realred */
      SUBROUTINE REALRED(LUNIT,NZRED,NREDM,REDE,REDGP,EVECR,
     &     IBTYP,IBEVC,IBPVC,IBCVC,NPRED)
C***********************************************************************
C
C     Read info needed for backsteps
C
C     Written by J. Thyssen - Dec 28 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
C
      DIMENSION REDE(*), REDGP(*)
      DIMENSION IBTYP(2,*), IBEVC(*), IBPVC(*), IBCVC(*)
      DIMENSION EVECR(*)
C
      CALL QENTER('REALRED')
C
C     Dump info
C
      READ(LUNIT) NREDM,
     &     NZRED,
     &     (REDE(I),I=1,NREDM*NREDM),
     &     (REDGP(I),I=1,NREDM),
     &     (EVECR(I),I=1,NREDM),
     &     (IBTYP(1,I),IBTYP(2,I),I=1,NREDM),
     &     (IBEVC(I),I=1,NREDM),
     &     (IBPVC(I),I=1,NREDM),
     &     (IBCVC(I),I=1,NREDM),
     &     NPRED
C
      CALL QEXIT('REALRED')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmchi */
      SUBROUTINE KRMCHI
C***********************************************************************
C
C     Write welcome message.
C
C     Written by J. Thyssen - Jan 18 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "maxorb.h"
#include "dcbidx.h"
#include "dcborb.h"
#include "dcbopt.h"
C
      IF (MCTYPE .EQ. JDHF .OR. MCTYPE .EQ. JODHF) THEN
         CALL TITLER('QC-DHF wave function optimization','*',199)
      ELSE IF (MCTYPE .EQ. JCAS) THEN
         CALL TITLER('KR-CAS MCSCF wave function','*',199)
      ELSE IF (MCTYPE .EQ. JGAS) THEN
         CALL TITLER('KR-GAS MCSCF wave function','*',199)
      ELSE
         WRITE(LUPRI,'(A,I3)')
     &        ' *** ERROR in KRMCHI *** Unknown MCTYPE : ',MCTYPE
         CALL QUIT('*** ERROR in KRMCHI *** Unknown MCTYPE!')
      END IF
C
      WRITE(LUPRI,9000)
 9000 FORMAT(
     &     3X,'This is output from DIRAC KR-MCSCF',/,
     &     3X,' - a relativistic four-component second order ',
     &     'restricted-step ',/,
     &     6X,'MCSCF wave function optimization program.',
     &     ///,
     &     3X,'General structure: ',/,
     &     3X,'  Joern Thyssen and Hans Joergen Aa. Jensen ',//,
     &     3X,'Linear equation and eigenvalue solver: ',/,
     &     3X,'  Hans Joergen Aa. Jensen, Trond Saue, ',
     &     'and Joern Thyssen',//,
     &     3X,'Integral transformation: ',/,
     &     3X,'  Luuk Visscher, Jon K. Laerdahl, and Trond Saue',//,
     &     3X,'Linear symmetry implementation:',/,
     &     3X,'  Stefan Knecht and Hans Joergen Aa. Jensen',//,
     &     3X,'KRMC-GASCIP CI code: ',/,
     &     3X,'  Joern Thyssen and Hans Joergen Aa. Jensen',//,
     &     3X,'LUCIAREL CI code: ',/,
     &     3X,'  Timo Fleig and Jeppe Olsen',//,
     &     3X,'Parallel LUCIAREL CI code: ',/,
     &     3X,'  Stefan Knecht, Hans Joergen Aa. Jensen and Timo Fleig',
     &     ///,
     &     3X,'``Aha!''''                             ',
     &     '- Fysiska poem, Linkoepings Universitet',/,
     &     1X,79('*'))
*
      WRITE(LUPRI,9110)
 9110 FORMAT(/' This module is published in:',//,
     &     4X,'KRMC THEORY:    H J Aa Jensen, K G Dyall, T Saue,',
     &           ' and K Faegri Jr.',/,
     &    20X,'   J Chem Phys 104 (1995) 4083',/,
     &     4X,'KRMC IMPLEMENTATION:',/,
     &     4X,'                J Thyssen, T Fleig, and H J Aa Jensen',/,
     &    20X,'   J Chem Phys 129 (2008) 034109',/,
     &     4X,'DIRAC-LUCIAREL: T Fleig, J Olsen, and L Visscher',/,
     &    20X,'   J Chem Phys, 119 (2003) 2963',/,
     &     4X,'KRMC-LUCIAREL:  T Fleig, H J Aa Jensen, J Olsen,',
     &        ' and L Visscher',/,
     &    20X,'   J Chem Phys, 124 (2006) 104106',/,
     &     4X,'PARALLEL KRMC-LUCIAREL:',/,
     &     4X,'                S. Knecht, H J Aa Jensen, and T Fleig',/,
     &    20X,'   J Chem Phys, 132,1 (2010) 014108.',//,
     &     1X,79('*'),//)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmchi */
      SUBROUTINE SETKRMC(LUPRI)
C***********************************************************************
C
C     Assign logical units and open files.
C
C     Written by T. Fleig - Aug 25 2004
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "dcbdhf.h"
#include "dcbgen.h"
C
      LOGICAL EX
C
C     Needs to be done here as SETDHF is not called in KRMC run:
C
      LUFCK1 = 10
C
C     Units 50-59 are reserved for KRMC.
C     Units 50-56 are opened here, 57-58 in krmctra.F.
C     (hjaaj March 2002).
C     59 is currently used for KRMCOLD, restarting option.
C     (tf August 2003)
C     LUKRMC = 50 is also used in ESR_RESOLVE
C
C     If old KRMCSCF exists rename it to KRMCOLD.
C
C
C     LUKRMC is set globally in dirac/dirset.F
C
      CALL OPNFIL(LUKRMC,'KRMCSCF','UNKNOWN','PSIOPT')
      CALL NEWLAB('*KRMCSCF',LUKRMC,LUPRI)
C
      LUKRM1 = 51
      CALL OPNFIL(LUKRM1,'KRMC_KAPPA','UNKNOWN','PSIOPT')
C
      LUKRM2 = 52
      CALL OPNFIL(LUKRM2,'KRMC_CMO','UNKNOWN','PSIOPT')
      CALL NEWLAB('SODLABEL',LUKRM2,LUPRI)
      REWIND LUKRM2
C
      LUKRM3 = 53
      INQUIRE ( FILE = 'KRMC_FOCK', EXIST = EX )
      IF( EX ) THEN
        CALL OPNFIL(LUKRM3,'KRMC_FOCK','OLD','PSIOPT')
      ELSE
        CALL OPNFIL(LUKRM3,'KRMC_FOCK','UNKNOWN','PSIOPT')
        CALL NEWLAB('SODLABEL',LUKRM3,LUPRI)
      END IF
      REWIND LUKRM3
C
      LUKRM4 = 54
      CALL OPNFIL(LUKRM4,'KRMC_LRED','UNKNOWN','PSIOPT')
      REWIND LUKRM4
C
      LUITFO = 56
      CALL OPNFIL(LUITFO,'KRMC_ITINFO','UNKNOWN','PSIOPT')
C
      INQUIRE ( FILE = 'KRMCOLD', EXIST = EX )
      IF (EX) THEN
         LUKRM5 = 59
         CALL OPNFIL(LUKRM5,'KRMCOLD','OLD','PSIOPT')
         REWIND LUKRM5
      ELSE
         LUKRM5 = -59
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rh1diag */
      SUBROUTINE RH1DIAG(CMO,WORK,LWORK)
C***********************************************************************
C
C     Obtain orbitals from diagonalization of H1.
C
C     Written by J. Thyssen - Mar 3 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
#include "dcbbas.h"
C
      DIMENSION CMO(*)
      DIMENSION WORK(*)
C
#include "memint.h"
C
      CALL QENTER('RH1DIAG')
C
C     Get one-electron Fock matrix
C
      CALL MEMGET2('REAL','H1AO',KH1AO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','H1MO',KH1MO,N2ORBTQ,WORK,KFREE,LFREE)
      CALL DIRH1(WORK(KH1AO),IPROPT,WORK(KFREE),LFREE)
      IF (BNCRON .AND. (2*NISHT+NAELEC).GT.1) THEN
          CALL BNCORR(WORK(KH1AO),WORK(KFREE),LFREE,IPROPT)
      ENDIF
C
C     Transform H1DIAG to MO-basis
C     ----------------------------
C
C     1) Get TMAT (AO-MO transformation)
C
      LUTMAT = 3
      CALL MEMGET2('REAL','TMAT',KTMAT,N2TMT,WORK,KFREE,LFREE)
      CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','RH1DIAG')
      CALL READT(LUTMAT,N2TMT,WORK(KTMAT))
      CLOSE(LUTMAT,STATUS='KEEP')
C
C     2) Transform to MO basis
C
      DO I = 1, NFSYM
         CALL QTRANS('AOMO','S',D0,
     &        NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &        WORK(KH1AO+I2BASX(I,I)),
     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &        WORK(KH1MO+I2ORBT(I)),
     &        NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &        WORK(KTMAT+I2TMT(I)),
     &        NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &        WORK(KTMAT+I2TMT(I)),
     &        NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,IPROPT)
      END DO
C
C
C
C     3) Diagonalize H1DIAG / BNCORR
C
      CALL MEMGET2('REAL','EIG',KEIG,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','IBEIG',KIBEIG,NORBT,WORK,KFREE,LFREE)
      CALL DFDIAG(WORK(KH1MO),WORK(KEIG),WORK(KIBEIG),WORK(KH1AO),
     &            .FALSE.,WORK,KFREE,LFREE)
C
C     4) Backtransform to current MO basis
C
      DO I = 1, NFSYM
         IF(NORB(I).GT.0) THEN
C
            CALL BCKTRA(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &           WORK(KH1AO+I2ORBT(I)),NORB(I),NORB(I),
     &           NORB(I),NZ,NORB(I),1,NFBAS(I,0),
     &           WORK(KTMAT+I2TMT(I)),NFBAS(I,0),NORB(I),NZT,
     &           IPROPT)
C
         END IF
      END DO
C
      CALL MEMREL('RH1DIAG',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RH1DIAG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck RSMPAR */
      SUBROUTINE RSMPAR(GPOE,GPOP,JXOPE,NXOPE,NHOPE,JXOPP,NXOPP,
     &                  THRCOM,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     This routine is a direct rip from REDPAR (../prp/pamset.F),
C     except that the new vectors are returned along with
C     their new sizes. Also, we don't try to compress conf.
C     variables.
C     Also remember to handle redundant e-e rotations (active-active)
C     correctly.
C
C     Reduce number of parameters by deleting parameters
C     for which the corresponding gradient element is
C     below threshold THRCOM
C
C     Written by J. Thyssen - Tue Mar 13 14:52:12 "MET 2001
C
C     NOTE: This routine is not active - S. Knecht/Jan 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,DC = 100.0D0)
C
#include "dgroup.h"
C
      DIMENSION GPOE(*),GPOP(*)
      DIMENSION WORK(*)
      DIMENSION JXOPE(2,*),JXOPP(2,*)
C
      KFRSAV = KFREE
      IF(IPRINT.GE.3) WRITE(LUPRI,'(/A,1P,D15.6)')
     &     ' RSMPAR: compression threshold = ',THRCOM
C
      write(LUPRI,*) 'RSMPAR before'
      CALL RPRKAP(GPOE,GPOP,JXOPE,NHOPE,
     &     JXOPP,NXOPP,NZ,IPQTOQ(1,0),LUPRI)
      IF(NXOPE.GT.0) THEN
         CALL MEMGET2('REAL','BUFe',KBUF,NXOPE*NZ,WORK,KFREE,LFREE)
         II = 0
         III = 0
         DO I = 1,NHOPE
            DTST = D0
            IZOFF = 0
            DO IZ = 1,NZ
               DTST = DTST + GPOE(I+IZOFF)*GPOE(I+IZOFF)
               IZOFF = IZOFF + NXOPE
            ENDDO
            IF(SQRT(DTST).GE.THRCOM) THEN
               II = II + 1
               IF (I .LE. NXOPE) III = III + 1
               JXOPE(1,II) = JXOPE(1,I)
               JXOPE(2,II) = JXOPE(2,I)
               IOFF = KBUF-1
               JOFF = 0
               DO IZ = 1,NZ
                  WORK(IOFF+II) = GPOE(JOFF+I)
                  IOFF = IOFF + NXOPE
                  JOFF = JOFF + NXOPE
               ENDDO
            ENDIF
         ENDDO
         WRITE(LUPRI,9000)
         IF(II.LT.NXOPE) THEN
            IOFF = KBUF
            JOFF = 1
            DO IZ = 1,NZ
               CALL DCOPY(II,WORK(IOFF),1,GPOE(JOFF),1)
               IOFF = IOFF + NXOPE
               JOFF = JOFF + II
            ENDDO
            IF(IPRINT.GE.3) THEN
               DTST = DC*(NXOPE-II)/NXOPE
               WRITE(LUPRI,9001) 'e-e',NHOPE,NXOPE,II,III,DTST
            ENDIF
            NHOPE  = II
            NXOPE  = III
         ELSE
            WRITE(LUPRI,9001) 'e-e',NHOPE,NXOPE,NHOPE,NXOPE,D0
         ENDIF
         CALL MEMREL('REDPAR.ee',WORK,1,KFRSAV,KFREE,LFREE)
      ENDIF
      write(LUPRI,*) 'RSMPAR after'
      CALL RPRKAP(GPOE,GPOP,JXOPE,NHOPE,
     &     JXOPP,NXOPP,NZ,IPQTOQ(1,0),LUPRI)
 9000 FORMAT(3X,'Type',T10,'Before',T30,'After',T50,'Reduction')
 9001 FORMAT(3X,A,T10,I5,T18,'(',I5,')',T30,I5,T38,'(',I5,')',
     &     T50,F5.1,'%')
C
C     e-p orbital rotations
C     =====================
C
      IF(NXOPP.GT.0) THEN
        CALL MEMGET2('REAL','BUFp',KBUF,NXOPP*NZ,WORK,KFREE,LFREE)
        II = 0
        DO I = 1,NXOPP
          DTST = D0
          IZOFF = 0
          DO IZ = 1,NZ
            DTST = DTST + GPOP(I+IZOFF)*GPOP(I+IZOFF)
            IZOFF = IZOFF + NXOPP
          ENDDO
          IF(SQRT(DTST).GE.THRCOM) THEN
            II = II + 1
            JXOPP(1,II) = JXOPP(1,I)
            JXOPP(2,II) = JXOPP(2,I)
            IOFF = KBUF-1
            JOFF = 0
            DO IZ = 1,NZ
              WORK(IOFF+II) = GPOP(JOFF+I)
              IOFF = IOFF + NXOPP
              JOFF = JOFF + NXOPP
            ENDDO
          ENDIF
        ENDDO
        WRITE(LUPRI,9000)
        IF(II.LT.NXOPP) THEN
          IOFF = KBUF
          JOFF = 1
          DO IZ = 1,NZ
            CALL DCOPY(II,WORK(IOFF),1,GPOP(JOFF),1)
            IOFF = IOFF + NXOPP
            JOFF = JOFF + II
          ENDDO
          IF(IPRINT.GE.3) THEN
            DTST = DC*(NXOPP-II)/NXOPP
            III  = II
            WRITE(LUPRI,9001) 'e-p',NXOPP,NXOPP,II,III,DTST
          ENDIF
          NXOPP  = II
        ENDIF
        CALL MEMREL('REDPAR.ep',WORK,1,KFRSAV,KFREE,LFREE)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DELPAR(orb_type_vec,JXOPE,NXOPE,NHOPE,JXOPP,NXOPP,
     &                  IPRINT)
C***********************************************************************
C
C     This routine is inspired by REDPAR (../prp/pamset.F).
C     Note that the new vectors are returned along with their new sizes.
C
C     purpose: delete parameters for which the corresponding orbital
C     --------------------------------------------------------------
C
C     - mj-values     --> linear spin-orbit calculation
C     - ml-values     --> linear spinfree   calculation
C     - boson irreps  -->        spinfree   calculation
C
C     do not match.
C
C     Written by S. Knecht - Feb 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,DC = 100.0D0)
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION JXOPE(2,*), JXOPP(2,*)
      integer orb_type_vec(NORBT)
C
      IF( IPRINT .ge. 5 )THEN
        WRITE(LUPRI,'(/A)')' (DELPAR): vector of orbital types:'//
     &              ' mj-/ml-/boson-irrep values'
        CALL IWRTMAMN(orb_type_vec,1,NORBT,1,NORBT,LUPRI)
        WRITE(LUPRI,'(/A)')' (DELPAR): orbital rotations'
        do i = 1, nhope
          write(lupri,'(a,i0,a,i0,a,i0)') ' rotation ',i,': orbital ',
     &    JXOPE(1,I),' --> ',JXOPE(2,I)
        end do
      END IF

      IF(NXOPE.GT.0) THEN
         II = 0
         III = 0
         DO I = 1,NHOPE
           IF(orb_type_vec(JXOPE(1,I)).eq.orb_type_vec(JXOPE(2,I)))THEN
             II = II + 1
             IF (I .LE. NXOPE) III = III + 1
             JXOPE(1,II) = JXOPE(1,I)
             JXOPE(2,II) = JXOPE(2,I)
           ENDIF
         ENDDO
         if(iprint .ge. 1)then
           WRITE(LUPRI,'(/A)')
     &   ' (DELPAR) ** statistics of deleted orbital rotations: **'
           WRITE(LUPRI,9000)
         end if
         IF(II.LT.NXOPE) THEN
            IF(IPRINT.GE.1) THEN
               DTST = DC*(NXOPE-II)/NXOPE
               WRITE(LUPRI,9001) 'e-e',NHOPE,NXOPE,II,III,DTST
            ENDIF
            NHOPE  = II
            NXOPE  = III
         ELSE
            if(iprint .ge. 1)then
              WRITE(LUPRI,9001) 'e-e',NHOPE,NXOPE,NHOPE,NXOPE,D0
            end if
         ENDIF
      ENDIF
 9000 FORMAT(3X,'Type',T10,'Before',T30,'After',T50,'Reduction')
 9001 FORMAT(3X,A,T10,I5,T18,'(',I5,')',T30,I5,T38,'(',I5,')',
     &     T50,F5.1,'%')

      IF( IPRINT .ge. 5 )THEN
        WRITE(LUPRI,'(/A)')' (DELPAR): orbital rotations (reduced)'
        do i = 1, nhope
          write(lupri,'(a,i0,a,i0,a,i0)') ' rotation ',i,': orbital ',
     &    JXOPE(1,I),' --> ',JXOPE(2,I)
        end do
      END IF
C
C     e-p orbital rotations
C     =====================
C
      IF(NXOPP.GT.0) THEN
        II  = 0
        DO I = 1,NXOPP
          DTST = D0
          IF(orb_type_vec(JXOPP(1,I)).eq.orb_type_vec(JXOPP(2,I)))THEN
            II = II + 1
            JXOPP(1,II) = JXOPP(1,I)
            JXOPP(2,II) = JXOPP(2,I)
          ENDIF
        ENDDO
        if(iprint .ge. 1)then
          WRITE(LUPRI,'(/A)')
     &   ' (DELPAR) ** statistics of deleted orbital rotations: **'
          WRITE(LUPRI,9000)
        end if
        IF(II.LT.NXOPP) THEN
          IF(IPRINT.GE.1) THEN
            DTST = DC*(NXOPP-II)/NXOPP
            III  = II
            WRITE(LUPRI,9001) 'e-p',NXOPP,NXOPP,II,III,DTST
          ENDIF
          NXOPP  = II
        ENDIF
      ENDIF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck fsmrot */
      SUBROUTINE FSMROT(NOROT,THR,GOEE,JXOPE,NXOPE,GOEP,JXOPP,NXOPP,
     &        WORK,LWORK)
C***********************************************************************
C
C     Find orbitals we can freeze.
C
C     Written by J. Thyssen - Tue Mar 13 14:52:12 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION GOEE(*),GOEP(*)
      DIMENSION WORK(*)
      DIMENSION JXOPE(2,*),JXOPP(2,*)
      DIMENSION NOROT(*)
C
#include "memint.h"
C
      CALL QENTER('FSMROT')
C
      CALL MEMGET2('REAL','GRAD',KGRAD,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Scatter gradients into matrix
C
      CALL DZERO(WORK(KGRAD),N2ORBXQ)
      CALL XOPSCT(GOEE,JXOPE,NXOPE,WORK(KGRAD),NORBT,NZ)
      CALL XOPSCT(GOEP,JXOPP,NXOPP,WORK(KGRAD),NORBT,NZ)
C
      CALL FULMAT('A',NORBT,NORBT,WORK(KGRAD))
      DO IZ = 2, NZ
         CALL FULMAT('S',NORBT,NORBT,WORK(KGRAD))
      END DO
C     write(6,*) 'fsmrot full matrix'
C     call prqmat(work(kgrad),norbt,norbt,norbt,norbt,nz,ipqtoq(1,0),6)
C
      NBEFORE = NORBT - ISUM(NORBT,NOROT,1)
C
      CALL FSMRO1(NOROT,THR,WORK(KGRAD))
C
      NAFTER = NORBT - ISUM(NORBT,NOROT,1)
C
      WRITE(LUPRI,'(A,I4,A,I4)')
     &     ' FSMROT: Orbitals before: ',NBEFORE,
     &     ', orbitals after: ',NAFTER
C
      CALL MEMREL('FSMROT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('FSMROT')
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck fsmro1 */
      SUBROUTINE FSMRO1(NOROT,THR,GRAD)
C***********************************************************************
C
C     Find orbitals we can freeze.
C
C     Written by J. Thyssen - Tue Mar 13 14:52:12 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "consts.h"
#include "maxash.h"
C
#include "maxorb.h"
#include "dcbidx.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION NOROT(*),GRAD(NORBT,NORBT,*)
C
C     write(6,*) 'FSMRO1 threshold',thr
      DO 10 I = 1, NORBT
C        ...don't remove occupied orbitals.
         IF (IOBTYP(I) .EQ. JTINAC .OR. IOBTYP(I) .EQ. JTACT) GOTO 10
         DN = D0
         DO IZ = 1, NZ
            DN = DN + DDOT(NORBT,GRAD(1,I,IZ),1,GRAD(1,I,IZ),1)
         END DO
         DN = SQRT(DN)
C        write(6,*) 'orbital ',I,' norm = ',dn
         IF (DN .LE. THR) NOROT(I) = 1
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck rmrot */
      SUBROUTINE RMROT(NOROT,GOEE,JXOPE,NXOPE,NHOPE,
     &        GOEP,JXOPP,NXOPP,IPRINT,WORK,LWORK)
C***********************************************************************
C
C     Remove orbital rotations.
C
C     Written by J. Thyssen - Tue Mar 13 14:52:12 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "consts.h"
      PARAMETER( DC = 1.0D2 )
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION GOEE(*),GOEP(*)
      DIMENSION WORK(*)
      DIMENSION JXOPE(2,*),JXOPP(2,*)
      DIMENSION NOROT(*)
C
#include "memint.h"
C
      CALL QENTER('RMROT')
C
      IF (IPRINT .GE. 3) WRITE(LUPRI,9000)
C
      IF (NHOPE .GT. 0) THEN
C
         CALL MEMGET2('REAL','BUFe',KBUF,NHOPE*NZ,WORK,KFREE,LFREE)
C
C
         II = 0
         III = 0
         DO I = 1, NHOPE
            J1 = JXOPE(1,I)
            J2 = JXOPE(2,I)
            IF (NOROT(J1) .EQ. 0 .AND. NOROT(J2) .EQ. 0 ) THEN
               II = II + 1
               IF (I .LE. NXOPE) III = III + 1
C
               JXOPE(1,II) = JXOPE(1,I)
               JXOPE(2,II) = JXOPE(2,I)
               IOFF = KBUF-1
               JOFF = 0
               DO IZ = 1,NZ
                  WORK(IOFF+II) = GOEE(JOFF+I)
                  IOFF = IOFF + NXOPE
                  JOFF = JOFF + NXOPE
               ENDDO
            END IF
         END DO
         IF(II.LT.NHOPE) THEN
            IOFF = KBUF
            JOFF = 1
            DO IZ = 1,NZ
               CALL DCOPY(II,WORK(IOFF),1,GOEE(JOFF),1)
               IOFF = IOFF + NHOPE
               JOFF = JOFF + II
            ENDDO
            IF(IPRINT.GE.3) THEN
               DTST = DC*(NHOPE-II)/NHOPE
               WRITE(LUPRI,9001) 'e-e',NHOPE,NXOPE,II,III,DTST
            ENDIF
            NHOPE  = II
            NXOPE  = III
         ELSE
            WRITE(LUPRI,9001) 'e-e',NHOPE,NXOPE,NHOPE,NXOPE,D0
         ENDIF
         CALL MEMREL('RMROT.ee',WORK,KWORK,KWORK,KFREE,LFREE)
      ENDIF
C
      IF (NXOPP .GT. 0) THEN
        CALL MEMGET2('REAL','BUFp',KBUF,NXOPP*NZ,WORK,KFREE,LFREE)
C
         II = 0
         DO I = 1, NXOPP
            J1 = JXOPP(1,I)
            J2 = JXOPP(2,I)
            IF (NOROT(J1) .EQ. 0 .AND. NOROT(J2) .EQ. 0 ) THEN
               II = II + 1
C
               JXOPP(1,II) = JXOPP(1,I)
               JXOPP(2,II) = JXOPP(2,I)
               IOFF = KBUF-1
               JOFF = 0
               DO IZ = 1,NZ
                  WORK(IOFF+II) = GOEP(JOFF+I)
                  IOFF = IOFF + NXOPP
                  JOFF = JOFF + NXOPP
               ENDDO
            END IF
         END DO
         WRITE(LUPRI,9000)
         IF(II.LT.NXOPP) THEN
            IOFF = KBUF
            JOFF = 1
            DO IZ = 1,NZ
               CALL DCOPY(II,WORK(IOFF),1,GOEP(JOFF),1)
               IOFF = IOFF + NXOPP
               JOFF = JOFF + II
            ENDDO
            IF(IPRINT.GE.3) THEN
               DTST = DC*(NXOPP-II)/NXOPP
               WRITE(LUPRI,9002) 'e-p',NXOPP,II,DTST
            ENDIF
            NXOPP  = III
         ELSE
            WRITE(LUPRI,9002) 'e-p',NXOPP,NXOPP,D0
         ENDIF
         CALL MEMREL('RMROT.ep',WORK,KWORK,KWORK,KFREE,LFREE)
      END IF
C
 9000 FORMAT(3X,'Type',T10,'Before',T30,'After',T50,'Reduction')
 9001 FORMAT(3X,A,T10,I5,T18,'(',I5,')',T30,I5,T38,'(',I5,')',
     &     T50,F5.1,'%')
 9002 FORMAT(3X,A,T10,I5,T30,I5,T50,F5.1,'%')
C
      CALL QEXIT('RMROT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck svldeg */
      SUBROUTINE SVLDEG(GOEE,GOEP,JXOPE,NXOPE,JXOPP,NXOPP,
     &        DV,FC,FV,FQ,CMO,IPRINT,WORK,LWORK)
C***********************************************************************
C
C     Do singular value decomposition of gradient.
C
C     Written by J. Thyssen - Thu Mar 15 15:15:58 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "maxorb.h"
#include "dcbnrt.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION GOEE(NXOPE,*), GOEP(NXOPP,*), JXOPE(*), JXOPP(*)
      DIMENSION FC(*), FV(*), CMO(*), FQ(NORBT,NASHT,*), DV(*)
      DIMENSION WORK(*)
C
      DIMENSION JMO(2), NMO(2)
      DIMENSION DFAC(2)
C
#include "memint.h"
C
      CALL QENTER('SVLDEG')
C
C     Make canonical orbitals for inactive and virtual orbitals
C     ---------------------------------------------------------
C     write(lupri,*) 'SVLDEG: fq before trans'
C     call prqmat(fq,norbt,nasht,norbt,nasht,nz,
C    &     ipqtoq(1,0),lupri)
C
      CALL MEMGET2('REAL','FD',KFD,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','U',KU,N2ORBXQ,WORK,KFREE,LFREE)
C
      CALL DCOPY(N2ORBXQ,FC,1,WORK(KFD),1)
      CALL DAXPY(N2ORBXQ,D1,FV,1,WORK(KFD),1)
C
      CALL SVL_CANORB('I',WORK(KFD),WORK(KU),WORK(KFREE),LFREE)
C
      CALL ROTFC('I',WORK(KU),FC,WORK(KFREE),LFREE)
      CALL ROTFC('I',WORK(KU),FV,WORK(KFREE),LFREE)
      CALL ROTFQ('I',WORK(KU),FQ,WORK(KFREE),LFREE)
      CALL ROTORB('I',WORK(KU),CMO,WORK(KFREE),LFREE)
C
      CALL SVL_CANORB('S',WORK(KFD),WORK(KU),WORK(KFREE),LFREE)
C
      CALL ROTFC('S',WORK(KU),FC,WORK(KFREE),LFREE)
      CALL ROTFC('S',WORK(KU),FV,WORK(KFREE),LFREE)
      CALL ROTFQ('S',WORK(KU),FQ,WORK(KFREE),LFREE)
      CALL ROTORB('S',WORK(KU),CMO,WORK(KFREE),LFREE)
C
C     write(lupri,*) 'SVLDEG fq after trans'
C     call prqmat(fq,norbt,nasht,norbt,nasht,nz,
C    &     ipqtoq(1,0),lupri)
C
C
C     Calculate trial vector
C     ----------------------
C
      CALL MEMGET2('REAL','BEE',KBEE,NXOPE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','BEP',KBEP,NXOPP*NZ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','HEE',KHEE,NXOPE,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','HEP',KHEP,NXOPP,WORK,KFREE,LFREE)
      LFPV = NASHT*NORBT*NZ
      CALL MEMGET2('REAL','FPV',KFPV,LFPV,WORK,KFREE,LFREE)
C
C     Calculate diagonal part of Hessian
C
      CALL FGEN_PV(FC,FQ,DV,WORK(KFPV),IPQTOQ(1,0),IPRINT)
      CALL E2ODI1(WORK(KHEE),FC,FV,DV,WORK(KFPV),JXOPE,NXOPE)
      CALL E2ODI1(WORK(KHEP),FC,FV,DV,WORK(KFPV),JXOPP,NXOPP)
C
C     Save generalized Fock matrix (F_{pv})
C
      CALL WRTKRMC(LUKRM3,'FGENPV  ',WORK(KFPV),LFPV)
C
C     Trial vector
C     (TODO FIXME: avoid divi by zero)
C
      J = 0
      DO IZ = 1, NZ
         DO I = 1, NXOPE
            J = J + 1
            WORK(KBEE+J-1) = GOEE(I,IZ) / WORK(KHEE+I-1)
         END DO
      END DO
      J = 0
      DO IZ = 1, NZ
         DO I = 1, NXOPP
            J = J + 1
            WORK(KBEP+J-1) = GOEP(I,IZ) / WORK(KHEP+I-1)
         END DO
      END DO
C
C
C     Calculate transformation to SV orbitals
C     ---------------------------------------
C
      NOROTC(1:NORBT) = NOROT(1:NORBT)
C
      IF (NXOPE .GT. 0) THEN
C
         CALL SVLDE1('E',.TRUE.,.TRUE.,
     &        GOEE,WORK(KBEE),JXOPE,NXOPE,
     &        WORK(KU),NOROTC,WORK(KFREE),LFREE)
C
         CALL ROTFC('S',WORK(KU),FC,WORK(KFREE),LFREE)
         CALL ROTFC('S',WORK(KU),FV,WORK(KFREE),LFREE)
         CALL ROTFQ('S',WORK(KU),FQ,WORK(KFREE),LFREE)
         CALL ROTORB('S',WORK(KU),CMO,WORK(KFREE),LFREE)
C
      END IF
      IF (NXOPP .GT. 0) THEN
C
         CALL SVLDE1('P',.TRUE.,.TRUE.,
     &        GOEP,WORK(KBEP),JXOPP,NXOPP,
     &        WORK(KU),NOROT,WORK(KFREE),LFREE)
C
         CALL ROTFC('P',WORK(KU),FC,WORK(KFREE),LFREE)
         CALL ROTFC('P',WORK(KU),FV,WORK(KFREE),LFREE)
         CALL ROTFQ('P',WORK(KU),FQ,WORK(KFREE),LFREE)
         CALL ROTORB('P',WORK(KU),CMO,WORK(KFREE),LFREE)
C
      END IF
C
!     CALL IWRTMAMN(NOROTC,NORBT,1,NORBT,1,LUPRI)
C
C
C     Calculate new gradient
C
      DFAC(1) = D1
      DFAC(2) = DM1
C
      IF (NXOPE .GT. 0) THEN
         CALL DZERO(GOEE,NZ*NXOPE)
         CALL RORBGRD(.TRUE.,DV,FC,FC,FV,FQ,GOEE,
     &        JXOPE,NXOPE,WORK(KFPV),DFAC,
     &        IPQTOQ(1,0),IPRINT)
      END IF
      IF (NXOPP .GT. 0) THEN
         CALL DZERO(GOEP,NZ*NXOPP)
         CALL RORBGRD(.TRUE.,DV,FC,FC,FV,FQ,GOEE,
     &        JXOPP,NXOPP,WORK(KFPV),DFAC,
     &        IPQTOQ(1,0),IPRINT)
      END IF
C
C     write(lupri,*) 'SVLDEG: new gradient'
C     CALL RPRKAP(GOEE,GOEP,JXOPE,NXOPE,JXOPP,NXOPP,
C    &     NZ,IPQTOQ(1,0),LUPRI)
C
      CALL MEMREL('SVLDEG',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SVLDEG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck svlde1 */
      SUBROUTINE SVLDE1(C,DOGRAD,DOTRI,
     &     GO,BO,JXOP,NXOP,U,NOROT,WORK,LWORK)
C***********************************************************************
C
C     Do singular value decomposition of gradient.
C
C     Written by J. Thyssen - Thu Mar 15 15:15:58 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "dummy.h"
C
#include "dcborb.h"
#include "dgroup.h"
C
      CHARACTER*1 C
      LOGICAL   DOGRAD, DOTRI
      DIMENSION GO(*), BO(*), JXOP(2,*), U(*), NOROT(*), WORK(*)
C
#include "memint.h"
C
      CALL QENTER('SVLDE1')
C
      IF (C .EQ. 'E') THEN
         CALL MEMGET2('REAL','VECe',KVEC,2 * (NASHT+NISHT) * NSSHT,
     &        WORK,KFREE,LFREE)
         NDIM = NSSHT
      ELSE IF (C .EQ. 'P') THEN
         CALL MEMGET2('REAL','VECp',KVEC,2 * (NASHT+NISHT) * NPSHT,
     &        WORK,KFREE,LFREE)
         NDIM = NPSHT
      ELSE
         WRITE(LUPRI,'(2A)')
     &        '*** ERROR in SVLDE1 *** Unknown control character: ',C
         CALL QUIT('*** ERROR in SVLDE1 *** Unknown control character')
      END IF
C
      IF (DOGRAD) THEN
         IF (DOTRI) THEN
            CALL SVLDE2(C,DOGRAD,DOTRI,
     &           GO,BO,JXOP,NXOP,U,NOROT,
     &           WORK(KVEC),2*NISHT+2*NASHT,NDIM,WORK(KFREE),LFREE)
         ELSE
            CALL SVLDE2(C,DOGRAD,DOTRI,
     &           GO,DUMMY,JXOP,NXOP,U,NOROT,
     &           WORK(KVEC),NISHT+NASHT,NDIM,WORK(KFREE),LFREE)
         END IF
      ELSE
         IF (DOTRI) THEN
            CALL SVLDE2(C,DOTRI,DOGRAD,
     &           BO,DUMMY,JXOP,NXOP,U,NOROT,
     &           WORK(KVEC),NISHT+NASHT,NDIM,WORK(KFREE),LFREE)
         ELSE
            write(LUPRI,*) 'SVLDE1 ERROR: both DOGRAD and DOTRI false!'
            call quit('*** ERROR in SVLDE1 ***')
         END IF
      END IF
C
      CALL MEMREL('SVLDE2',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SVLDE1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck svlde2 */
      SUBROUTINE SVLDE2(C,DOGRAD,DOTRI,
     &                  GO,BO,JXOP,NXOP,U,NOROT,VEC,NR,NC,WORK,LWORK)
C***********************************************************************
C
C     Do singular value decomposition of gradient.
C
C     Written by J. Thyssen - Thu Mar 15 15:15:58 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "dcborb.h"
#include "maxorb.h"
#include "dcbidx.h"
#include "dgroup.h"
C
      CHARACTER*1 C
      LOGICAL DOGRAD, DOTRI
      DIMENSION GO(NXOP,*), BO(NXOP,*)
      DIMENSION JXOP(2,*)
      DIMENSION U(*)
      DIMENSION NOROT(*)
      DIMENSION VEC(NR,NC,NZ)
      DIMENSION WORK(*)
C
#include "memint.h"
C
      CALL QENTER('SVLDE2')
C
      CALL DZERO(VEC,NR*NC*NZ)
C
      IF (C .EQ. 'E') THEN
C
         DO I = 1, NXOP
            K = JXOP(1,I)
            L = JXOP(2,I)
            IF (IOBTYP(L) .EQ. JTSEC) THEN
               IF (IOBTYP(K) .EQ. JTINAC) THEN
                  KK = IDXG2I(K)
               ELSE
                  KK = IDXG2U(K) + NISHT
               END IF
               LL = IDXG2S(L)
C
               IF (DOGRAD) THEN
                  DO IZ = 1, NZ
                     VEC(KK,LL,IZ) = GO(I,IZ)
                  END DO
               END IF
               IF (DOTRI) THEN
                  DO IZ = 1, NZ
                     VEC(KK+NISHT+NASHT,LL,IZ) = BO(I,IZ)
                  END DO
               END IF
            END IF
         END DO
C
      ELSE
         call quit('SVLDE2: only C="E" implemented yet, write me!')
      END IF
C
      CALL MEMGET2('REAL','VEC',KVEC,NC*NC*NZ,WORK,KFREE,LFREE)
C
C
C     Calculate V^{\dagger} V
C     -----------------------
C
      CALL QGEMM(NC,NC,NR,D1,
     &     'H','N',IPQTOQ(1,0),VEC,NR,NC,NZ,
     &     'N','N',IPQTOQ(1,0),VEC,NR,NC,NZ,
     &     D0,IPQTOQ(1,0),WORK(KVEC),NC,NC,NZ)
C
C
C     Diagonalize matrix
C     ------------------
C
      CALL MEMGET2('REAL','EIG',KEIG,NC,WORK,KFREE,LFREE)
C
      IF (NZ .EQ. 1) THEN
         IJOB = 1
         IORDER = 0
         IPACK = 0
         CALL RSJACO(NC,NC,NC,WORK(KVEC),WORK(KEIG),
     &        IJOB,IORDER,IPACK,U)
      ELSE
         MATZ = 1
         CALL QDIAG(NZ,NC,WORK(KVEC),NC,NC,
     &        WORK(KEIG),MATZ,U,NC,NC,
     &        WORK(KFREE),LFREE,IERR)
         IF (IERR .NE. 0) THEN
            WRITE(LUPRI,'(/2A,I4)')
     &           ' *** ERROR in SVLDE2 ***: ',
     &           'QDIAG failed with error code ',IERR
            CALL QUIT('*** ERROR in SVLDE2 ***')
         END IF
      END IF
C
      write(LUPRI,*) 'SVLDE2: eigenvalues'
      CALL OUTPUT(WORK(KEIG),1,NC,1,1,NC,1,1,LUPRI)
C
      IF (C .EQ. 'E') THEN
         DO I = 1, NSSHT
            IF ( ABS(WORK(KEIG+I-1)) .LE. THRZER) THEN
               II = IDXS2G(I)
               write(LUPRI,*) 'removing virtual ',I,II
               NOROT(II) = 1
            END IF
         END DO
      ELSE
         DO I = 1, NPSHT
            IF ( ABS(WORK(KEIG+I-1)) .LE. THRZER) THEN
               II = IDXP2G(I)
               write(LUPRI,*) 'removing positronic ',I,II
               NOROT(II) = 1
            END IF
         END DO
      END IF
C
      CALL MEMREL('SVLDE2',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SVLDE2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck svl_canorb */
      SUBROUTINE SVL_CANORB(C,FD,U,WORK,LWORK)
C***********************************************************************
C
C     Diagonalize the full FD Fock matrix (*not* symmetry packed!)
C     and return the transformation matrix U.
C     The character "C" controls which part of FD we diagonalize.
C
C     Written by J. Thyssen - Mar 2001
C
C     This routine is used by SVLDEG.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "consts.h"
C
#include "maxorb.h"
#include "dcbidx.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      CHARACTER C*1
      DIMENSION FD(NORBT,NORBT,NZ), U(*), WORK(*)
C
#include "memint.h"
C
      CALL QENTER('SVL_CANORB')
C
C
      IF (C .EQ. 'E') THEN
         CALL MEMGET2('REAL','FDe',KFD,NESHT*NESHT*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFD),NESHT*NESHT*NZ)
         CALL MATGAT(FD,NORBT,NORBT,WORK(KFD),NESHT,NESHT,IDXE2G,NZ)
         NDIM = NESHT
      ELSE IF (C .EQ. 'P') THEN
         CALL MEMGET2('REAL','FDp',KFD,NPSHT*NPSHT*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFD),NPSHT*NPSHT*NZ)
         CALL MATGAT(FD,NORBT,NORBT,WORK(KFD),NPSHT,NPSHT,IDXP2G,NZ)
         NDIM = NPSHT
      ELSE IF (C .EQ. 'I') THEN
         CALL MEMGET2('REAL','FDi',KFD,NISHT*NISHT*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFD),NISHT*NISHT*NZ)
         CALL MATGAT(FD,NORBT,NORBT,WORK(KFD),NISHT,NISHT,IDXI2G,NZ)
         NDIM = NISHT
      ELSE IF (C .EQ. 'A') THEN
         CALL MEMGET2('REAL','FDa',KFD,NASHT*NASHT*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFD),NASHT*NASHT*NZ)
         CALL MATGAT(FD,NORBT,NORBT,WORK(KFD),NASHT,NASHT,IDXU2G,NZ)
         NDIM = NASHT
      ELSE IF (C .EQ. 'S') THEN
         CALL MEMGET2('REAL','FDs',KFD,NSSHT*NSSHT*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFD),NSSHT*NSSHT*NZ)
         CALL MATGAT(FD,NORBT,NORBT,WORK(KFD),NSSHT,NSSHT,IDXS2G,NZ)
         NDIM = NSSHT
      ELSE
         WRITE(LUPRI,'(2A)')
     &        '*** ERROR in MATGAT *** Unknown control character: ',C
         CALL QUIT('*** ERROR in MATGAT *** Unknown control character')
      END IF
C
C
C     Diagonalize matrix
C     ------------------
C
      CALL MEMGET2('REAL','EIG',KEIG,NDIM,WORK,KFREE,LFREE)
C
      IF (NZ .EQ. 1) THEN
         IJOB = 1
         IORDER = 0
         IPACK = 0
         CALL RSJACO(NDIM,NDIM,NDIM,WORK(KFD),WORK(KEIG),
     &        IJOB,IORDER,IPACK,U)
      ELSE
         MATZ = 1
         CALL QDIAG(NZ,NDIM,WORK(KFD),NDIM,NDIM,
     &        WORK(KEIG),MATZ,U,NDIM,NDIM,
     &        WORK(KFREE),LFREE,IERR)
         IF (IERR .NE. 0) THEN
            WRITE(LUPRI,'(/,1X,2A,I4)')
     &           '*** ERROR in SVL_CANORB ***: ',
     &           'QDIAG failed with error code ',IERR
            CALL QUIT('*** ERROR in SVL_CANORB ***')
         END IF
      END IF
C
C     write(6,*) 'eigenvalues'
C     CALL OUTPUT(WORK(KEIG),1,NDIM,1,1,NDIM,1,1,6)
C
c     write(6,*) 'eigenvectors'
c     CALL PRQMAT(U,NDIM,NDIM,NDIM,NDIM,NZ,IPQTOQ(1,0),LUPRI)
C
      CALL MEMREL('SVL_CANORB',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('SVL_CANORB')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRMC_WAKE_NODES
C
C     Call the slaves if we run in parallel for KR-MCSCF
C
C     adaption of the corresp. RELCCSD routine written by Luuk Visscher
C     KR-MCSCF adaption: S. Knecht
C
C     Last revision: S. Knecht           - June 2008, Odense
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "dcbgen.h"
#include "parkrmc.h"
C
C     ITASK = 7 for parallel KR-MCSCF
C
      IF (PARCAL) CALL DIRAC_PARCTL( KRMCSCF_PAR )
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRMC_RELEASE_NODES
C
C     Release the slaves if we run in parallel for KR-MCSCF
C
C     adaption of the corresp. RELCCSD routine written by Luuk Visscher
C     KR-MCSCF adaption: S. Knecht
C
C     Last revision: S. Knecht           - June 2008, Odense
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "dcbgen.h"
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
      SUBROUTINE KRMC_NODE_DRIVER()
C
C     KR-MCSCF driver routine for the slaves
C
C     Last revision: S. Knecht           - June  2008, Odense
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use memory_allocator
C
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
#include "parkrmc.h"
      INTEGER ISTAT(MPI_STATUS_SIZE)
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
      CHARACTER KRMC_NODEFN*17, KRMC_NODEF*13
      LOGICAL KRMCEND

      real(8), allocatable :: WORK(:)
      call legacy_lwork_get(LWORK)
      call alloc(WORK,LWORK,id='WORK in KRMC_NODE_DRIVER')
#include "memint.h"
C
      KFRSAV = KFREE
C     Add the master node, NUMNOD = number of slaves
      NUMCPU = NUMNOD + 1
C     run parallel MCSCF
      PARMCSCF = .TRUE.
C
*
*     create a node-unique filename as output file. Important on
*     shared file systems. Otherwise all the output gets mingled in one
*     file. You don't really want to do this.
*
      KRMC_NODEF="KRMC_NODE.OUT"
*
      IF (MYTID .LT. 10) THEN    ! MPI ID has one digit
         WRITE (KRMC_NODEFN,'(A13,A1,I1)') KRMC_NODEF,'.',MYTID
         IFILELNG=15
      ELSE IF (MYTID .LT. 100) THEN  ! MPI ID has two digits
         WRITE (KRMC_NODEFN,'(A13,A1,I2)') KRMC_NODEF,'.',MYTID
         IFILELNG=16
      ELSE IF (MYTID .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (KRMC_NODEFN,'(A13,A1,I3)') KRMC_NODEF,'.',MYTID
         IFILELNG=17
      ELSE
         CALL Abend2("NUMCPU.GT.1000! EXTEND KRMC_NODE_DRIVER MODULE")
      ENDIF
*
*     open the local input file and the node specific output file.
*     Every access to the local stdout handle then automatically writes
*     to the corresponding output file.
*
*
      OPEN(IKRMCPARF,FILE = KRMC_NODEFN(1:IFILELNG))
C
C     transfer file handle to common block
C
      LMCPRI = IKRMCPARF
C
C
C     LUCIAREL run is finished for the slaves ...
C
      CLOSE (IKRMCPARF,STATUS='KEEP')
C
      call dealloc(WORK)
      END
#endif
