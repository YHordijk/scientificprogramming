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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRPXLR */
      SUBROUTINE PRPXLR(WORK,KLRF,KCNV,KLRO,KFREE,LFREE)
C***********************************************************************
C
C
C     Calculate linear response values defined in /XCBXLR/
C
C     Linear response functions are defined as:
C
C        <<A,;B>>_w = -(E[1]_A)^{+}(E[2] - wS[2])^{-1}(E[1]_B)
C
C        E[1]_A   is gradient of perturbation operator A  (h|A)
C        E[1}_B   is gradient of perturbation operator B  (h|B)
C        E[2]     is the molecular Hessian                (h|H|h)
C        S[2]     is a metric                             (h|h)
C
C     where h is the operator manifold. In RPA h is limited to single
C     excitations.
C
C     The linear response functions are found by iterative techniques.
C     We write
C
C        (E[2] - wS[2])^{-1}(E[1]_B) = X
C
C     which can be rearranged to the linear response equation
C
C        (E[2] - wS[2])X = E[1]_B
C
C     The linear reponse equation is solved by projecting X onto
C     a reduced space of trial vectors:
C
C        X = (b_i)a_i           B = [(b_1)  (b_2)  ...  (b_n)]
C
C     and solving the reduced equation
C
C        (RE[2] - wRS[2])a = RE[1]_B
C
C     in which
C
C       RE[2] = B^{+}E[2]B
C       RS[2] = B^{+}S[2]B
C       RE[1] = B^{+}E[1}
C
C
C     The evaluation of linear reponse functions thus proceeds in
C     three steps:
C
C       1. Generate property gradient      : E[1]_B
C       2. Solve reduced response equation : (RE[2] - wRS[2])a = RE[1]_B
C       3. Generate full solution vector   : X
C
C
C     Each operator A and B is assumed to span a given boson irrep
C     and be either symmetric(+) or antisymmetric(-) with respect
C     to time reversal.
C
C     Written by Trond Saue and Hans Joergen Aa. Jensen July 1996
C     Last revision July 9 1996 - tsaue
C
C***********************************************************************

      use dirac_cfg, only: dirac_cfg_scf_calculation

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)

! we use nopen
#include "dcbdhf.h"

#include "dcbgen.h"
#include "dcborb.h"
#include "dcbxlr.h"
#include "dgroup.h"
      DIMENSION WORK(*)
C
      CALL QENTER('PRPXLR')

      if (dirac_cfg_scf_calculation .and. nopen > 0) then
         call quit('ERROR: open-shell SCF response not implemented')
      end if

      NRSPF  = NLRBPT*NBFREQ
      IF(NRSPF.EQ.0) GOTO 10
      IF(XLR_SKIPEE) THEN
        NEE = 0
      ELSE
        NEE = 1
      ENDIF
      IF(XLR_SKIPEP) THEN
        NEP = 0
      ELSE
        NEP = 1
      ENDIF
      NCI = 0
      NDAMP = 1
      IF(DMPFRLR.NE.D0) NDAMP = 2
      N2LRF  = NLRAPT*NRSPF*(1+NEE+NEP+NCI)*NDAMP
      CALL MEMGET2('REAL','ABLRF',KLRF,N2LRF,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KLRF),N2LRF)
      CALL MEMGET2('REAL','RCNV',KCNV,NRSPF,WORK,KFREE,LFREE)
C
      CALL PRPXL1(WORK(KLRF),WORK(KCNV),WORK,KFREE,LFREE)
C
 10   CONTINUE
      IF (XLRANA) CALL ANAXRS(WORK,KFREE,LFREE,IPRXLR)
C
      CALL QEXIT('PRPXLR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRPXL1 */
      SUBROUTINE PRPXL1(ABLRF,RCNV,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Called by PRPXLR (see PRPXLR for description of task)
C
C     2002/2003, MI & HJAaj(&TEC) added changes for LAO based
C     calculation of magnetic properties
C
C***********************************************************************
      use london_direct_hb
      use london_reorth_hb
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
      use memory_allocator
      use orbital_rotation_indices
      use dirac_cfg

#ifdef HAS_PCMSOLVER
      use pcm_linear_response, only: pcm_linear_response_initialize
#endif


#include "implicit.h"
#include "priunit.h"
      PARAMETER(GPTHRS = 1.0D-9,D0 = 0.0D0,D2=2.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbxlr.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbnmr.h"

C
      LOGICAL   LBUF
      CHARACTER TTYP(-1:1)*2

      integer, allocatable :: ibeig(:)
      real(8), allocatable :: cmo(:)
      logical :: parcal_save
C
C svillaume: NFREQ=0 right now, not defined until call to SETXRS
C it should be ABLRF(NBFREQ,NLRAPT,NLRBPT,*) no?
      DIMENSION ABLRF(NFREQ,NLRAPT,NLRBPT,*),RCNV(NBFREQ,NLRBPT),
     &          WORK(*)
      SAVE TTYP
      DATA TTYP /'T-','T0','T+'/
C
      CALL QENTER('PRPXL1')
      KFRSAV = KFREE
C
C     set control variable for REDPAR (see REDPAR for an explanation of
C     IREDJOB)
      IREDJOB = 6
C
C     ****************************
C     *** I N I T I A L I Z E  ***
C     ****************************
C
      CALL SETXLR(WORK,KFREE,LFREE)

      IF(XLRANA) THEN
        CALL MEMGET2('INTE','IORBCLS',KLRO,4*NORBT,WORK,KFREE,LFREE)
        CALL DEFORB(WORK(KLRO),NCLS_XLR,XLR_INDSTR,
     &              WORK,KFREE,LFREE,IPRXLR)
      ELSE
        KLRO = KFREE
      ENDIF

#ifdef HAS_PCMSOLVER
      if (dirac_cfg_pcm) then
              call pcm_linear_response_initialize(lupri)
        write(lupri, '(//A/)') 'PCMSolver initialized.'
      end if
#endif
C
C     **********************************************************
C     ***                                                    ***
C     *** L I N E A R   R E S P O N S E   F U N C T I O N S  ***
C     ***                                                    ***
C     **********************************************************
      if (london) then
        if (.not.regrad) then

          call alloc(cmo,   ncmotq)
          call alloc(ibeig, ntbas(0))

c.........direct contribution
          call onelon(cmo, ibeig, work, kfree, lfree)
          
          if (.not.notwol) then
            call london_twoel_hb(cmo, ibeig, work(kfree), lfree)
          end if

c.........reorthonormalization contribution
!         gosia: in connection-independent formulation, there is no
!         reorthonormalization contribution in the property gradient.
!         Instead, there is a contribution from magnetic derivatives 
!         of overlap matrix (modifications in london_reorth.F90)
          if (.not. noonei) then
            call london_reort_hb(cmo, ibeig, work(kfree), lfree)
          end if
         
          call dealloc(cmo)
          call dealloc(ibeig)

        endif
      endif

C
C     Memory allocation
C     =================
C
C     Dimensions of reduced system
C
      KXLR = KFREE
      NEVECR = NREDM*NFREQ*NDAMP
      CALL MEMGET2('INTE','IBTYP',KIBTYP,2*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBCVC',KIBCVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBEVC',KIBEVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBPVC',KIBPVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EVALR',KEVALR,NFREQ  ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EVECR',KEVECR,NEVECR ,WORK,KFREE,LFREE)
C
C     Transfer frequencies to EVALR
C     =============================
C
      CALL DCOPY(NFREQ,BFREQ,1,WORK(KEVALR),1)
C
C       Loop over operators
C       ===================
C
      DO 10 JOPSY = 1,NFSYM
        NFC = NBSYM/NFSYM
        DO 20 IS = 1,NFC
          JSYMOP  = JFSYM(IS,JOPSY)
          NOPB    = NLRBP(JSYMOP)
          JOPB    = JLRBP(JSYMOP)
          DO 30 IOPB = 1,NOPB
C
C
C           Initialize configurational and orbital parameters
C           =================================================
C           (Note that this routine can also be shifted one 
C            loop up, but by placing it here local restrictions
C            are possible, see e.g. REDPAR below....)
C
            CALL XRSPAR(XLR_INDSTR,XLR_SKIPEE,XLR_SKIPEP,
     &                  IPRXLR)
            INDBP  = LLBPSU(JOPB+IOPB)
            INDPRB = LLRBPU(INDBP)
            JTIMOP = IPRPTIM(INDPRB)
            WRITE (LUPRI,'(///A/A,A,A3,2X,A2/)')
     &        ' <<<  SOLVING SETS OF LINEAR EQUATIONS '//
     &        'FOR LINEAR RESPONSE PROPERTY >>>',
     &        PRPNAM(INDPRB),
     &        ' - Symmetry: ',REP(JSYMOP-1),TTYP(JTIMOP)
C
C           ******************************************
C           *** P R O P E R T Y    G R A D I E N T ***
C           ******************************************
C
            KGRAD = KFREE
            CALL MEMGET2('REAL','GPCI',KGPCI,NZCONFQ,WORK,KFREE,LFREE)
            CALL MEMGET2('REAL','GPOE',KGPOE,NZXOPEQ,WORK,KFREE,LFREE)
            CALL MEMGET2('REAL','GPOP',KGPOP,NZXOPPQ,WORK,KFREE,LFREE)
            IF(NOSPIB) THEN
              LBUF = NOSPIN
              NOSPIN = .TRUE.
            ENDIF
            CALL PAMGRD(INDPRB,WORK(KGPCI),WORK(KGPOE),WORK(KGPOP),
     &                  JOPSY,
     &                  get_orbital_rotation_indices_pp(),
     &                  get_orbital_rotation_indices_pn(),
     &                  NZCONF,NZXOPE,NZXOPP,WORK,KFREE,LFREE,IPRXLR)
            IF(NOSPIB) NOSPIN = LBUF
C
C           Check norm of property gradient; skip if below threshold
C           ========================================================
C
            GPCNRM = D2*DNRM2(NZCONFQ,WORK(KGPCI),1)
            GPENRM = D2*DNRM2(NZXOPEQ,WORK(KGPOE),1)
            GPPNRM = D2*DNRM2(NZXOPPQ,WORK(KGPOP),1)
            GPTNRM = SQRT(GPCNRM**2+GPENRM**2+GPPNRM**2)
            IF(GPTNRM.LT.GPTHRS) THEN
              WRITE(LUPRI,'(A,A16)')
     &          'WARNING: Skipping B property ',PRPNAM(INDPRB)
              WRITE(LUPRI,'(9X,A,1P,D8.1)')
     &          'Norm of gradient  :',GPTNRM,
     &          'is below threshold:',GPTHRS
              GOTO 30
            ENDIF
C
C           Compress number of orbital rotations
C           ====================================
C
            IF(THRCOM.NE.D0) THEN
              CALL REDPAR(IREDJOB,WORK(KGPCI),WORK(KGPOE),WORK(KGPOP),
     &                    get_orbital_rotation_indices_pp(),
     &                    get_orbital_rotation_indices_pn(),
     &                    NZXOPE,NZXOPEQ,NZXOPP,NZXOPPQ,NZXOPT,NZXOPTQ,
     &                    NZVAR,NZVARQ,NZCONF,NZCONFQ,THRCOM,
     &                    WORK,KFREE,LFREE,IPRXLR)
            ENDIF
C
C           ********************************************************
C           *** S O L V E   R E S P O N S E    E Q U A T I O N S ***
C           ********************************************************
C
            NCRED  = 0
            NERED  = 0
            NPRED  = 0
            NZRED  = NCRED + NERED + NPRED
            RSPLAB = PRPNAM(INDPRB)
            CALL XRSCTL(WORK(KGPCI) ,WORK(KGPOE),WORK(KGPOP),
     &                 WORK(KIBTYP),WORK(KIBCVC),WORK(KIBEVC),
     &                 WORK(KIBPVC),RCNV(1,INDBP),
     &                 WORK(KEVALR),WORK(KEVECR),
     &                 WORK,KFREE,LFREE)
C          ... reuse KGP space for solution vectors.
            CALL MEMREL('PRPXL1.XRSCTL',WORK,1,KGRAD,KFREE,LFREE)
C
C           ************************************************************
C           *** E V A L U A T E   R E S P O N S E   F U N C T I O N  ***
C           ************************************************************
C
            IF(NOSPIA) THEN
              LBUF = NOSPIN
              NOSPIN = .TRUE.
            ENDIF
            CALL GETLRF(ABLRF,WORK(KIBTYP),WORK(KIBCVC),
     &                  WORK(KIBEVC),WORK(KIBPVC),RCNV(1,INDBP),
     &                  WORK(KEVALR),WORK(KEVECR),WORK(KLRO),
     &                  WORK,KFREE,LFREE)
            IF(NOSPIA) NOSPIN = LBUF
 30       CONTINUE
 20     CONTINUE
 10   CONTINUE
      CALL MEMREL('PRPXL1',WORK,1,KXLR,KFREE,LFREE)
C
C     Final print section
C
      CALL XLROUT(ABLRF,RCNV,WORK,KFREE,LFREE)
      CLOSE(LUXVC,STATUS = 'KEEP')
      CALL MEMREL('PRPXL1',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('PRPXL1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xlrinp */
      SUBROUTINE XLRINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for linear response  <<A;B>>
C
C     Called from:  PRPINP
C
C     Written by T.Saue and H.J.Aa.Jensen 1996
C     Last revision: June 14 1996 - tsaue
C
C***********************************************************************
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
#ifdef HAS_PELIB
      use pe_variables, only: peqm
#endif
      use dirac_cfg
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 44,DRESFAC = 1.0D2)
C
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxlr.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
C
      LOGICAL SET, NEWDEF, RESET, FINTER, FLIST, INTFLG_CHANGE
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,
     &          MXFORM*6,FMT*6
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.PRINT ','.A OPER','.B OPER','.OPERAT',
     &            '.B FREQ','.ALLCMB','.THRESH','.MAXITR',
     &            '.SKIPEP','.MAXRED','.TRIPLE','.XLRNRM',
     &            '.ONLYSF','.ONLYSG','.INTFLG','.NOPREC',
     &            '.OCCUP ','.VIRTUA','.ACTIVE','.COMPRE',
     &            '.ITRINT','.CNVINT','.REAXVC','.TRIAB ',
     &            '.RESFAC','.ANALYZ','.STERNH','.SKIPEE',
     &            '.MXLOAD','.A_NOSP','.B_NOSP','.UNCOUP',
     &            '.IMAGIN','.GAUSS-','.E2CHEK','.STERNC',
     &            '.DAMPIN','.FREQ I','.LAOMOD','.EPOLE ',
     &            '.MPOLE ','.NMFREQ','.ANATHR','.XXXXXX'/
      DATA SET/.FALSE./
C
#include "ibtfun.h"
#include "gauss_legendre.h"
C
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
! Miro: satisfy runtimecheck for parallel run
      FINTER     = .FALSE. 
C
C     Local initialization
C
      ILLINT = IBTAND(INTGEN,1)
      ISLINT = IBTAND(INTGEN/2,1)
      ISSINT = IBTAND(INTGEN/4,1)
      IGTINT = IBTAND(INTGEN/8,1)
      INTFLG_CHANGE = .FALSE.
C
C
C     Initialize /XCBLR/
C     ===================
C
      IPRXLR     = 0
      THCLR      = 1.0D-5
      THRCOM     = D0
      RESXLR     = DRESFAC
      CNVXLR(1)  = DUMMY
      CNVXLR(2)  = DUMMY
      MAXRM      = -1
      MXLOAD     = -1
      ITRXLR     = 50
      INTXLR     = ILLINT+2*ISLINT+4*ISSINT+8*IGTINT
CMI ... set up only LL class of integrals 
C        as the default for LevyLe/BSS
      IF (LEVYLE.OR.BSS.or.x2c) THEN
         INTXLR = ILLINT
      ENDIF
      ITRILR(1)  = 1
      ITRILR(2)  = 1
      NLRAPT     = 0
      NLRBPT     = 0
      CALL IZERO(NLRAP,8)
      CALL IZERO(NLRBP,8)
      IF (.NOT. DOVDW) THEN
C     Ordinary linear response calculation
         XLR_IMFREQ = .FALSE.
         NBFREQ     = 1
         BFREQ(1)   = D0
         DAMPFRLR   = D0
         STARTFREQ  = D0
         ENDFREQ    = D0
         STEPFREQ   = D1
         FINTER     = .FALSE.
         FLIST      = .FALSE.
      ELSE
C     Calculation of C6 dispersion coefficient
         XLR_IMFREQ = .TRUE.
         NBFREQ = 11
         BFREQ(1) = 0.0D0
         DO IFREQ=1,10
            BFREQ(IFREQ+1)=FREQ_10PT(IFREQ)
         END DO
      END IF
      ALLCMB     = .FALSE.
      XLR_SKIPEP = .FALSE.
!     if number of positronic shells is zero skip e-p rotations. hjj+sk - aug 2010
      if(x2c.or.bss.or.levyle.or.freepj.or.vextpj) 
     &  XLR_SKIPEP = .TRUE.
#ifdef HAS_PCMSOLVER
      if(dirac_cfg_pcm) 
     &  XLR_SKIPEP = .TRUE.
#endif
      XLR_SKIPEE = .FALSE.
      XLR_LSFG(1)= .TRUE.
      XLR_LSFG(2)= .TRUE.
      XLRNRM     = .FALSE.
      XLRDIH     = .TRUE.
      TRIAB      = .FALSE.
      XLR_XVCFIL = 'NORXVC'
      XLR_UNCOUP = .FALSE.
      XLRANA     = .FALSE.
      XLR_ANATHR = 2.0D0
      XSTERN     = .FALSE.
      XSTERNC    = .FALSE.
      NOSPIA     = .FALSE.
      NOSPIB     = .FALSE.
      XLR_E2CHEK = .FALSE.
      DO I = 1,NFSYM
         XLR_INDSTR(2,I) = ' '
         XLR_INDSTR(1,I) = ' '
         XLR_INDSTR(3,I) = ' '
      ENDDO
      XLR_TRIPLET = .FALSE.
C
C     Initialize user-defined linear response properties
C     ==================================================
C
C     Process input from CBIXLR
C     =========================
C
      NEWDEF = (WORD .EQ. '*LINEAR')
      ICHANG = 0
      INPERR = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
  110       CONTINUE
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24,25,26,27,28,29,30,
     &                     31,32,33,34,35,36,37,38,39,40,
     &                     41,42,43,44), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/3A/)') ' Keyword "',WORD,
     *            '" not recognized in XLRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in XLRINP.')
    1          CONTINUE
C&&&& PRINT:  Print level
                  READ(LUCMD,*) IPRXLR
               GO TO 100
    2          CONTINUE
C&&&& A OPERATOR: Left side operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        IPRXLR)
                  CALL OP1IND('NLRAPT',IND1OP,LLRAPU,NLRAPT,INDXPR,
     &                        MAXLLR)
               GO TO 100
    3          CONTINUE
C&&&& B OPERATOR: Right side operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        IPRXLR)
                  CALL OP1IND('NLRBPT',IND1OP,LLRBPU,NLRBPT,INDXPR,
     &                        MAXLLR)
               GO TO 100
    4          CONTINUE
C&&&& OPERATOR: Both left/right side operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        IPRXLR)
                  CALL OP1IND('NLRAPT',IND1OP,LLRAPU,NLRAPT,INDXPR,
     &                        MAXLLR)
                  CALL OP1IND('NLRBPT',IND1OP,LLRBPU,NLRBPT,INDXPR,
     &                        MAXLLR)
               GO TO 100
    5          CONTINUE
C&&&& B FREQuency
                  IF (DOVDW) THEN
C                 if van der Waals calculation, then frequencies are
C                 set by the choice of quadrature scheme
                     READ (LUCMD,*)
                     READ (LUCMD,*)
                     WRITE (LUPRI,'(/A)') 'Warning: frequencies ignored'
                     GOTO 100
                  END IF
                  FLIST = .TRUE.
                  READ (LUCMD, *) NBFREQ
                  IF (NBFREQ.GT.0) THEN
                    IF (NBFREQ.LE.MAXFLR) THEN
                       READ (LUCMD,*) (BFREQ(J),J=1,NBFREQ)
                    ELSE
                       INPERR = INPERR + 1
                       WRITE (LUPRI,'(3(/A,I5),/)')
     &                 '@ Number of B frequencies specified  :',
     &                 NBFREQ,
     &                 '@ is greater than the allowed number :',
     &                 MAXFLR
                       READ (LUCMD,*) (BFREQ(J),J=1,MAXFLR),
     &                              (FFFF,J=MAXFLR+1,NBFREQ)
                       NBFREQ = MAXFLR
                    END IF
                  END IF
               GO TO 100
    6          CONTINUE
C&&&& ALLCMB: form all possible <<A;B>>_w even if imaginary
               ALLCMB = .TRUE.
               GO TO 100
    7          CONTINUE
C&&&& THRESH - Threshold for convergence:
                  READ (LUCMD,*) THCLR
               GO TO 100
    8          CONTINUE
C&&&& MAXITR - Maximum number of iterations
                  READ (LUCMD, '(I5)') ITRXLR
               GO TO 100
    9          CONTINUE
C&&&& SKIPEP -exclude all electron-positron rotations
                  XLR_SKIPEP = .TRUE.
               GOTO 100
 10            CONTINUE
C&&&& MAXRED: max dimension of matrix in reduced system
                  READ(LUCMD,*) MAXRM
               GOTO 100
 11            CONTINUE
C&&&& TRIPLET: allow triplet excitations in SPINFREE runs
               XLR_TRIPLET = .TRUE.
               GOTO 100
 12            CONTINUE
C&&&& XLRNRM: normalize trial vectors
                 XLRNRM = .TRUE.
               GOTO 100
 13            CONTINUE
C&&&& ONLYSF: Only generate one-index transformed Fock matrix
                  XLR_LSFG(1) = .FALSE.
               GOTO 100
 14            CONTINUE
C&&&& ONLYSG: Only generate G matrix (Fock-matrix of modified density)
                  XLR_LSFG(2) = .FALSE.
               GOTO 100
C&&&& INTFLG: specify what two-electron integrals to include
 15            CONTINUE
                  IF (GAUNT) THEN
                     READ(LUCMD,*) ILLINT,ISLINT,ISSINT,IGTINT
                     INTXLR = ILLINT+2*ISLINT+4*ISSINT+8*IGTINT
                  ELSE
                     READ(LUCMD,*) ILLINT,ISLINT,ISSINT
                     INTXLR = ILLINT+2*ISLINT+4*ISSINT
                  ENDIF
                  INTFLG_CHANGE=.TRUE.
               GOTO 100
 16            CONTINUE
C&&&& NOPREC: no preconditioning of initial trial vectors
                  XLRDIH = .FALSE.
               GOTO 100
 17            CONTINUE
C&&&& OCCUP: String of inactive orbitals in XLR module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') XLR_INDSTR(1,I)
                  ENDDO
               GO TO 100
 18            CONTINUE
C&&&& VIRTUA: String of secondary orbitals in XLR module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') XLR_INDSTR(3,I)
                  ENDDO
               GO TO 100
 19            CONTINUE
C&&&& ACTIVE: String of active orbitals in XLR module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') XLR_INDSTR(2,I)
                  ENDDO
               GO TO 100
 20            CONTINUE
C&&&& COMPRESS: Compress number of orbital variations by checking
C               gradient elements againt a threshold
                  READ(LUCMD,*) THRCOM
               GO TO 100
 21            CONTINUE
C&&& ITRINT: Number of iterations before adding SL- and SS-integrals
                  READ(LUCMD,*) ITRILR(1),ITRILR(2)
               GO TO 100
 22            CONTINUE
C&&& CNVINT: Convergence thresholds for adding SL- and SS-integrals
                  READ(LUCMD,*) CNVXLR(1),CNVXLR(2)
               GO TO 100
 23            CONTINUE
C&&&& REAXVC: Read solution vectors from file XLR_XVCFIL
                  READ(LUCMD,'(A6)') XLR_XVCFIL
               GO TO 100
 24            CONTINUE
C&&&& TRIAB:   Enforce triangularity of <<A;B>>
                  TRIAB = .TRUE.
               GO TO 100
 25            CONTINUE
C&&&& RESFAC:  New trial vectors only for parameters whose residual
C              has a norm that is a fraction RESFAC of the max norm
                  READ(LUCMD,*) RESXLR
               GO TO 100
 26            CONTINUE
C&&&& ANALYZe: Analyze individual contributions to response functions
                  XLRANA = .TRUE.
               GO TO 100
 27            CONTINUE
C&&&& STERNH: The ep,ep-Hessian approximated by -2mc^2 diagonal matrix
                  XSTERN = .TRUE.
               GO TO 100
 28            CONTINUE
C&&&& SKIPEE: Skip all e-e rotations
                  XLR_SKIPEE = .TRUE.
               GO TO 100
 29            CONTINUE
C&&&& MXLOAD: Maximum number of trial vectors to be loaded in each microiteration
                  READ(LUCMD,*) MXLOAD
               GO TO 100
 30            CONTINUE
C&&&& NOSPIA: Eliminate all spin-dependene of property operator A
                 NOSPIA = .TRUE.
               GO TO 100
 31            CONTINUE
C&&&& NOSPIB: Eliminate all spin-dependene of property operator B
                 NOSPIB = .TRUE.
               GO TO 100
 32            CONTINUE
C&&&& UNCOUP: Uncoupled Hartree-Fock
                 XLR_UNCOUP = .TRUE.
               GO TO 100
 33            CONTINUE
C&&& IMAGINARY: Purely imaginary frequencies
                 XLR_IMFREQ = .TRUE.
               GO TO 100
 34            CONTINUE
C%%% GAUSS-LEGENDRE: specify number of quadrature points
               READ(LUCMD,*) NBFREQ
               NBFREQ = NBFREQ + 1
               BFREQ(1) = 0.0D0
               IF (NBFREQ-1.EQ.6) THEN
                  DO IFREQ=1,NBFREQ-1
                     BFREQ(IFREQ+1)=FREQ_6PT(IFREQ)
                  END DO
               ELSE IF (NBFREQ-1.EQ.8) THEN
                  DO IFREQ=1,NBFREQ-1
                     BFREQ(IFREQ+1)=FREQ_8PT(IFREQ)
                  END DO
               ELSE IF (NBFREQ-1.EQ.10) THEN
                  DO IFREQ=1,NBFREQ-1
                     BFREQ(IFREQ+1)=FREQ_10PT(IFREQ)
                  END DO
               ELSE IF (NBFREQ-1.EQ.12) THEN
                  DO IFREQ=1,NBFREQ-1
                     BFREQ(IFREQ+1)=FREQ_12PT(IFREQ)
                  END DO
               ELSE IF (NBFREQ-1.EQ.14) THEN
                  DO IFREQ=1,NBFREQ-1
                     BFREQ(IFREQ+1)=FREQ_14PT(IFREQ)
                  END DO
               ELSE IF (NBFREQ-1.EQ.16) THEN
                  DO IFREQ=1,NBFREQ-1
                     BFREQ(IFREQ+1)=FREQ_16PT(IFREQ)
                  END DO                
               ELSE
                  NBFREQ=11
                  DO IFREQ=1,NBFREQ-1
                     BFREQ(IFREQ+1)=FREQ_10PT(IFREQ)
                  END DO
                  WRITE(LUPRI,'( 3(/A) )') 
     &' Warning! Incorrect number of quadrature points.',
     &' Gauss-Legendre integration only implemented for'//
     &' 6,8,10,12,14,16 points.', 'Default of 10 points will be used.'
               END IF
               GO TO 100
 35            CONTINUE
C&&& E2CHEK: explicit contruction of the electronic Hessian E2
               XLR_E2CHEK = .TRUE.
               XLRDIH = .FALSE.
               GO TO 100
 36            CONTINUE
C&&& STERNC: Sternheimer complement
               XSTERNC = .TRUE.
               GO TO 100
 37            CONTINUE
C&&& DAMPING: specify a damping parameter
               READ(LUCMD,*) DMPFRLR
#ifdef HAS_PELIB
               !edh testing...
               IF (PEQM) CALL QUIT('PE damped response not enabled yet')
#endif
               GO TO 100
 38            CONTINUE
C&&& FREQ INTERVAL: specify frequencies with start, end, step values
               FINTER = .TRUE.               
               READ(LUCMD,*) STARTFREQ,ENDFREQ,STEPFREQ
               GO TO 100
 39            CONTINUE
C&&& LAOMOD: rearrange LR equations with LAOs, so connection matrices are not needed
#ifdef MOD_LAO_REARRANGED
               call set_london_keywords('laomod')
#endif
               GO TO 100
 40            CONTINUE
C&&&& EPOLE : Electric multipoles of order L
               READ(LUCMD,*) IORDER
               NDEG=(IORDER+2)*(IORDER+1)/2
               IF((NLRAPT+NDEG).GT.MAXLLR) THEN
                 WRITE(LUPRI,'(A,A,I3)') 'XLRINP: Pointer array LLR',
     &           'APU out of bounds for electric multipole of order ',
     &           IORDER
                 WRITE(LUPRI,'(A,I5)') 
     &           'Increase MAXLLR to ',(NLRAPT+NDEG)
                 CALL QUIT('XLRINP: MAXLLR out of bounds')
               ENDIF
               IF((NLRBPT+NDEG).GT.MAXLLR) THEN
                 WRITE(LUPRI,'(A,A,I3)') 'XLRINP: Pointer array LLR',
     &           'BPU out of bounds for electric multipole of order ',
     &           IORDER
                 WRITE(LUPRI,'(A,I5)') 
     &           'Increase MAXLLR to ',(NLRBPT+NDEG)
                 CALL QUIT('XLRINP: MAXLLR out of bounds')
               ENDIF
               CALL DEF_EPOLE(IORDER,LLRAPU(NLRAPT+1),IPRXLR)
               CALL ICOPY(NDEG,LLRAPU(NLRAPT+1),1,LLRBPU(NLRBPT+1),1)
               NLRAPT = NLRAPT + NDEG
               NLRBPT = NLRBPT + NDEG
               GOTO 100
 41            CONTINUE
C&&&& MPOLE : Magnetic multipoles of order L 
               READ(LUCMD,*) IORDER
               NDEG=3*(IORDER+1)*IORDER/2
               IF((NLRAPT+NDEG).GT.MAXLLR) THEN
                 WRITE(LUPRI,'(A,A,I3)') 'XLRINP: Pointer array LLR',
     &           'APU out of bounds for magnetic multipole of order ',
     &           IORDER
                 WRITE(LUPRI,'(A,I5)') 
     &           'Increase MAXLLR to ',(NLRAPT+NDEG)
               ENDIF
               IF((NLRBPT+NDEG).GT.MAXLLR) THEN
                 WRITE(LUPRI,'(A,A,I3)') 'XLRINP: Pointer array LLR',
     &           'BPU out of bounds for magnetic multipole of order ',
     &           IORDER
                 WRITE(LUPRI,'(A,I5)') 
     &           'Increase MAXLLR to ',(NLRBPT+NDEG)
               ENDIF
               CALL DEF_MPOLE(IORDER,LLRAPU(NLRAPT+1),IPRXLR)
               CALL ICOPY(NDEG,LLRAPU(NLRAPT+1),1,LLRBPU(NLRBPT+1),1)
               NLRAPT = NLRAPT + NDEG
               NLRBPT = NLRBPT + NDEG
               GOTO 100
 42            CONTINUE
C&&&& NMFREQuency: input frequencies in nanometers
                  IF (DOVDW) THEN
C                 if van der Waals calculation, then frequencies are
C                 set by the choice of quadrature scheme
                     READ (LUCMD,*)
                     READ (LUCMD,*)
                     WRITE (LUPRI,'(/A)') 'Warning: frequencies ignored'
                     GOTO 100
                  END IF
                  FLIST = .TRUE.
                  READ (LUCMD, *) NBFREQ
                  IF (NBFREQ.GT.0) THEN
                    IF (NBFREQ.LE.MAXFLR) THEN
                       READ (LUCMD,*) (BFREQ(J),J=1,NBFREQ)
                       DO J = 1,NBFREQ
                         BFREQ(J)=XTNM/BFREQ(J)
                       ENDDO
                    ELSE
                       INPERR = INPERR + 1
                       WRITE (LUPRI,'(3(/A,I5),/)')
     &                 '@ Number of B frequencies specified  :',
     &                 NBFREQ,
     &                 '@ is greater than the allowed number :',
     &                 MAXFLR
                       READ (LUCMD,*) (BFREQ(J),J=1,MAXFLR),
     &                              (FFFF,J=MAXFLR+1,NBFREQ)
                       NBFREQ = MAXFLR
                    END IF
                  END IF
               GO TO 100
 43            CONTINUE
C&&&& ANATHR: Threshold for analysis of linear response function
               READ(LUCMD,*) XLR_ANATHR
               GO TO 100
 44            CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/3A/)') ' Prompt "',WORD,
     *            '" not recognized in XLRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in XLRINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C     ===============
C
C specific complex response process section
      IF (XLR_IMFREQ .AND. DMPFRLR.NE.D0) THEN
      	XLR_IMFREQ = .FALSE.
      	WRITE(LUPRI,'( 2(/A) )') ' Warning!! .IMAGINARY and .DAMPING'//
     & 	' both specified',
     &	'keyword .IMAGINARY ignored assuming complex frequencies!'
      END IF
      IF (FINTER .AND. FLIST) THEN
        FLIST = .FALSE.
        WRITE(LUPRI,'( 2(/A) )') ' Warning!! .B FREQ and .FREQ I'//
     &  ' both specified',
     &  'keyword .B FREQ ignored!'
      END IF
      IF (FINTER) THEN
        NBFREQ = IFIX(0.5+SNGL((ENDFREQ-STARTFREQ)/STEPFREQ))+1
        IF (NBFREQ .GE. MAXFLR) THEN
     	  WRITE(LUPRI,'(2(A,I3/))') 
     &    'WARNING: Too many frequencies specified : ',NBFREQ,
     &    '         Maximum allowed (MAXFLR) is :    ',MAXFLR
          CALL QUIT('Too many frequencies specified!!!')
        END IF
        BFREQ(1)=STARTFREQ
        DO IFREQ=1,NBFREQ-1
          BFREQ(IFREQ+1)=STARTFREQ+IFREQ*STEPFREQ
        END DO
      END IF
      IF(DMPFRLR.NE.D0) THEN
C.....check frequencies; zeros not allowed
        DO IFREQ = 1,NBFREQ
          IF(BFREQ(IFREQ).EQ.D0) THEN
             WRITE(LUPRI,'(A)')
     &     'FATAL ERROR: do not combine zero frequencies with damping.',
     &     'Purely imaginary frequencies are handled by '//
     &     'keyword IMAGINARY'
           CALL QUIT(
     &     'XLRINP: Do not combine zero frequency with damping')
          ENDIF
        ENDDO
        WRITE(LUPRI,'(A,1P,D12.6)')
     &    ' * Damping frequency   :',DMPFRLR
      ENDIF
C end of complex freq processing
      IF(CNVXLR(1).LT.DUMMY) ITRILR(1) = 1
      IF(CNVXLR(2).LT.DUMMY) ITRILR(2) = 1
C
      IF (MAXRM.LT.0) MAXRM = MAX(200,100+NBFREQ*ITRXLR)
C        ... then calculate default MAXRM
      MXREDZ = MAXRM/2
      MAXRM  = 2*MXREDZ
C     ... also to make sure diamagnetic terms are calculated for magn.prop.
C           (see PRPDEF) /hj jan 2001
      IF(QED .AND. XLR_SKIPEP) THEN
         WRITE(LUPRI,'(/A)') ' FATAL ERROR for *LINEAR: '//
     &   'QED invalid for Levy-Leblond/no-pair/.SKIPEP'
         INPERR = INPERR + 1
      END IF
      IF(XSTERN) THEN
        XLR_UNCOUP = .TRUE.
        XLR_SKIPEE = .TRUE.
        XLR_SKIPEP = .FALSE.
        XLR_XVCFIL = 'NORXVC'
      ENDIF
      IF(XSTERNC) THEN
        XLR_UNCOUP = .TRUE.
        XLR_SKIPEE = .FALSE.
        XLR_SKIPEP = .TRUE.
        XLR_XVCFIL = 'NORXVC'
      ENDIF
      IF (XLR_UNCOUP) THEN
        ITRXLR = 0
        XLRDIH = .TRUE.
      ENDIF
C
C     Print section
C     =============
C
      DOXLR = DOXLR.OR.(NLRAPT.GT.0.AND.NLRBPT.GT.0.AND.NBFREQ.GT.0)
      IF(.NOT.(DOXLR.OR.XLRANA)) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &' XLRINP: Electronic linear reponse properties'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I5)')
     &   ' * Print level                 :',IPRXLR
      IF(XLR_UNCOUP) WRITE(LUPRI,'(A)') 
     &  ' * NB !!! Uncoupled calculation.'
      IF(XLRANA) WRITE(LUPRI,'(A)')
     &   ' * Analysis of solution vectors.'
      IF(XLR_E2CHEK) WRITE(LUPRI,'(A)')
     &   ' * Explicit contruction of the electronic Hessian'
      IF(XSTERN .AND. DOXLR) THEN
        WRITE(LUPRI,'(1X,A)')
     &'* Sternheimer: Approximating Hessian for ep-rotations by -2mc^2',
     &'  and ignoring coupling between e-e and e-p rotations.',
     &'  Note that .STERNH forces .MAXITR 0 and .SKIPEE and',
     &'  neglects any of .SKIPEP, .REAXVC, and .NOPREC'
        IF (LEVYLE.OR.BSS.or.x2c.OR.FREEPJ.OR.VEXTPJ.OR.XLR_SKIPEP)THEN
           WRITE (LUPRI,'(/A)') ' FATAL ERROR: '//
     &     '.STERNH is not valid for Levy-Leblond/no-pair/.SKIPEP'
           CALL QUIT(
     &     '.STERNH is not valid for Levy-Leblond/no-pair/.SKIPEP')
        END IF
      END IF
      IF(XSTERNC .AND. DOXLR) THEN
        WRITE(LUPRI,'(1X,A)')
     &'* Sternheimer complement: ',
     &'  Note that .STERNH forces .MAXITR 0 and .SKIPEP and',
     &'  neglects any of .SKIPEP, .REAXVC, and .NOPREC'
        IF (LEVYLE.OR.BSS.or.x2c.OR.FREEPJ.OR.VEXTPJ.OR.XLR_SKIPEE)THEN
           WRITE (LUPRI,'(/A)') ' FATAL ERROR: '//
     &     '.STERNC is not valid for Levy-Leblond/no-pair/.SKIPEE'
           CALL QUIT(
     &     '.STERNC is not valid for Levy-Leblond/no-pair/.SKIPEE')
        END IF
      END IF
      IF (XLR_SKIPEE) THEN
        WRITE(LUPRI,'(A)') ' * NB!!! No e-e rotations'
      ELSE
        WRITE(LUPRI,'(A)')
     &  ' * e-e rotations are allowed in the linear response function.'
      ENDIF
      IF(XLR_SKIPEP.AND..NOT.(LEVYLE.OR.BSS.or.x2c.OR.FREEPJ.OR.VEXTPJ))
     &     WRITE(LUPRI,'(A)')
     &     ' * NB!!! No rotations between electrons and positrons.'
#ifdef HAS_PCMSOLVER
      if(dirac_cfg_pcm) then
        WRITE(LUPRI,'(A)')
     &  ' * e-p rotations not allowed by default'//
     &  ' when using PCM.'
      endif
#endif
      IF (.NOT.XLR_SKIPEP) WRITE(LUPRI,'(A)')
     &  ' * e-p rotations are allowed in the linear response function.'
      IF(SPINFR.AND.XLR_TRIPLET) WRITE(LUPRI,'(A)')
     &  ' * Triplet excitations allowed.'
      IF (NOSPIA) WRITE(LUPRI,'(A)')
     &  ' * All spin-dependence removed from property operator A'
      IF (NOSPIB) WRITE(LUPRI,'(A)')
     &  ' * All spin-dependence removed from property operator B'
      IF (LEVYLE) WRITE(LUPRI,'(A)')
     &  ' * You are running in the nonrelativistic (Levy-Leblond) mode.'
      IF (BSS) WRITE(LUPRI,'(A)')
     &  ' * You are running in the 2-component BSS relativistic mode.'
      IF (DOPCT) WRITE(LUPRI,'(A)')
     &  '   with activated picture change transformation'//
     &  ' of four-component property operators.'
      IF(QED) WRITE(LUPRI,'(A)')
     &     ' * QED: rotations between positrons and secondary orbitals.'
      IF(.NOT.DOXLR) THEN
        DOXLR = DOXLR.OR.XLRANA
        GOTO 999
      ENDIF
      IF(XLR_XVCFIL.NE.'NORXVC') WRITE(LUPRI,'(A,A6)')
     &    ' * The program will attempt restart from solution vectors '//
     &    'on file ',XLR_XVCFIL
      WRITE(LUPRI,'(A,I5)')
     &   ' * Maximum number of iterations:',ITRXLR
      WRITE(LUPRI,'(A,1P,D12.6)')
     &   ' * Threshold for convergence   :',THCLR
      WRITE(LUPRI,'(A,I5)')
     &   ' * Max. size of reduced matrix :',MAXRM
      IF(TRIAB) WRITE(LUPRI,'(A)')
     &   ' * Assuming <<A;B>> to be symmetric.'
      IF(THRCOM.NE.D0) THEN
        WRITE(LUPRI,'(1X,A/1X,A/1X,A,1P,D12.6)')
     &  '* Reducing orbital parameters by deleting rotations for which',
     &  '  for which the corresponding property gradient is below',
     &  '  threshold ',THRCOM
      ENDIF
C
      IF(ALLCMB) WRITE(LUPRI,'(A)')
     &  ' * All combinations of A and B operators will be formed,'//
     &  ' even if not real.'
      WRITE(LUPRI,'(A)') ' * Include contributions from '//
     +    'the following two-electron integrals:'
      IF(ILLINT.EQ.1) WRITE(LUPRI,'(3X,A)') '- LL-integrals'
      IF ((LEVYLE.OR.BSS.or.x2c).AND.(ISLINT.EQ.1.OR.ISSINT.EQ.1)) THEN
        WRITE(LUPRI,'(3X,A)')
     & 'WARNING: Only LL integrals are allowed for LevyLe,BSS,x2c - '//
     & 'deactivating SL,SS,GT class.'
        ISLINT=0
        ISSINT=0
        IGTINT=0
        INTXLR=1
        INTFLG_CHANGE = .FALSE. ! deactivate indicator of integrals setting
      ENDIF
      IF(ISLINT.EQ.1) THEN
        IF(CNVXLR(1).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,1P,D8.1)')
     +    '- SL-integrals below convergence ',CNVXLR(1)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    '- SL-integrals from iteration ',ITRILR(1)
        ENDIF
      ENDIF
      IF(ISSINT.EQ.1) THEN
        IF(CNVXLR(2).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,1P,D8.1)')
     +    '- SS-integrals below convergence ',CNVXLR(2)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    '- SS-integrals from iteration ',ITRILR(2)
        ENDIF
      ENDIF
      IF (INTFLG_CHANGE) THEN
        WRITE(LUPRI,"(4X,A)")
     & "---> accepted user's setting through .INTFLG"
      ELSE
        WRITE(LUPRI,"(4X,A)")
     & "---> these are default values from Hamiltonian input"
      ENDIF
      IF(.NOT.XLR_LSFG(1)) WRITE(LUPRI,'(A)')
     &     ' * Debug: Only calling FMOLI in sigma-vector routine...'
      IF(.NOT.XLR_LSFG(2)) WRITE(LUPRI,'(A)')
     &     ' * Debug: Only calling GMOLI in sigma-vector routine...'
      IF(XLRNRM) WRITE(LUPRI,'(A)')
     &     ' * Trial vectors will be normalized.'
      IF(XLRDIH) THEN
        WRITE(LUPRI,'(A)')
     &  ' * Initial trial vectors will be preconditioned.'
      ELSE
        WRITE(LUPRI,'(A)')
     &  ' * No preconditioning of initial trial vectors.'
      ENDIF
      IF(MXLOAD.GT.0) WRITE(LUPRI,'(A,I5)')
     &  ' * Maximum number of trial vectors to be loaded'//
     &  ' in each microiteration:', MXLOAD
C
C     User specified A operators
C     ===========================
C
      IF (NLRAPT.GT.0) THEN
       CALL PRSYMB(LUPRI,'-',75,0)
       WRITE(LUPRI,'(A)') ' ** User specified A operators **'
       CALL PRSYMB(LUPRI,'-',75,0)
       DO I = 1,NLRAPT
        INDXPR = LLRAPU(I)
        CALL WRIXPR(I,INDXPR)
       ENDDO
      ENDIF
C
C     User specified B operators
C     ===========================
C
      IF (NLRBPT.GT.0) THEN
       CALL PRSYMB(LUPRI,'-',75,0)
       WRITE(LUPRI,'(A)') ' ** User specified B operators **'
       CALL PRSYMB(LUPRI,'-',75,0)
       DO I = 1,NLRBPT
        INDXPR = LLRBPU(I)
        CALL WRIXPR(I,INDXPR)
       ENDDO
      ENDIF
C
C     Frequencies
C     ===========
C
C     CALL PRSYMB(LUPRI,'-',75,0)
      IF(XLR_IMFREQ) THEN
        WRITE(LUPRI,'(A)')
     &' * Linear response functions evaluated at the following '//
     &   'IMAGINARY frequencies:'
      ELSE IF (DMPFRLR.NE.D0) THEN
      	WRITE(LUPRI,'(A)')
     &' * Linear response functions evaluated at the following '//
     &   'COMPLEX frequencies:'
      ELSE
        WRITE(LUPRI,'(A)')
     &' * Linear response functions evaluated at the following '//
     &   'REAL frequencies:'
      ENDIF
      DO J = 1,NBFREQ
        IF (DMPFRLR.NE.D0) THEN
          FMT = MXFORM(BFREQ(J),12)
          WRITE(LUPRI,'(6X,'//FMT//',A,'//FMT//',A)')
     &      BFREQ(J),' +',DMPFRLR,'i a.u.'
        ELSE
          FMT = MXFORM(BFREQ(J),12)
          WRITE(LUPRI,'(6X,'//FMT//',A)') BFREQ(J),' a.u.'
        ENDIF
      ENDDO
      IF(DOVDW) THEN
        WRITE(LUPRI,'(A)')
     &' * Calculation of C6 dispersion coefficient'
        WRITE(LUPRI,'(A,I2,A)')
     &' * Gauss-Legendre integration with ',NBFREQ-1,' points' 
      ENDIF
  999 CONTINUE
      IF (XLR_SKIPEE .AND. XLR_SKIPEP) THEN
         INPERR = INPERR + 1
         WRITE(LUPRI,'(/A)') ' FATAL ERROR: both SKIPEE and SKIPEP'
      END IF
      IF (INPERR.GT.0) CALL QUIT('Input error in *LINEAR RESPONSE')

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setxlr */
      SUBROUTINE SETXLR(WORK,KFREE,LFREE)
C***********************************************************************
C
C     Linear response module:
C  ================================
C     Initialize variables/Open necessary files etc.
C
C     Written by T.Saue and H.J.Aa.Jensen July 1996
C     Last revision: July 9 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxlr.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbxrs.h"
#include "dcbxpr.h"
C
      LOGICAL TOBE,ORBDEF(3)
      DIMENSION WORK(*),IR(2,2),NR(2)
      CHARACTER  TTYP(-1:1)*2
      DATA TTYP /'T-','T0','T+'/

      CALL QENTER('SETXLR')
      KFRSAV = KFREE

C
C     ************************************
C     *** Set orbital strings          ***
C     ************************************
C
      CALL OSTRING(XLR_INDSTR,XLR_SKIPEE,XLR_SKIPEP,IPRXLR,
     &             WORK,KFREE,LFREE)
C
C     ************************************
C     *** Sort operators on symmetries ***
C     ************************************
C

CMI/TODO heartfull wish July 2006: It would be good to have sorting of
C       A,B operators at the beginning of the run (after the XLRINP)

      CALL OP1SRT('LLAPSU',LLRAPU,LLAPSU,NLRAPT,NLRAP,JLRAP,
     &            WORK(KFREE),LFREE)

      CALL OP1SRT('LLBPSU',LLRBPU,LLBPSU,NLRBPT,NLRBP,JLRBP,
     &            WORK(KFREE),LFREE)

CMI   ... print out on what couples are going into the linear response...
      IF (IPRXLR.GE.7) THEN
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(A)') ' ** SETXLR: List of all A operators **'
        CALL PRSYMB(LUPRI,'-',75,0)
        DO I = 1,NLRAPT
          INDXPR = LLRAPU(I)
          CALL WRIXPR(I,INDXPR)
        ENDDO

        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(A)') ' ** SETXLR: List of all B operators **'
        CALL PRSYMB(LUPRI,'-',75,0)
        DO I = 1,NLRBPT
          INDXPR = LLRBPU(I)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      
        WRITE(LUPRI,'(A)')
     &  ' ** SETXLR: Sorted couples of A,B operators'//
     &  ' for linear response calculations **'

        DO IOPSY = 1,NFSYM
          NBF = NBSYM/NFSYM
          DO IBF = 1,NBF
            ISYM = JFSYM(IBF,IOPSY)
            NOPA = NLRAP(ISYM)
            JOPA = JLRAP(ISYM)
            NOPB = NLRBP(ISYM)
            JOPB = JLRBP(ISYM)
            DO IOPA = 1,NOPA
              INDAP  = LLAPSU(JOPA+IOPA)
              INDPRA = LLRAPU(INDAP)
              IREPA  = IPRPSYM(INDPRA)-1
              ITIMA  = IPRPTIM(INDPRA)
              DO IOPB = 1,NOPB
                INDBP  = LLBPSU(JOPB+IOPB)
                INDPRB = LLRBPU(INDBP)
                IREPB  = IPRPSYM(INDPRB)-1
                ITIMB  = IPRPTIM(INDPRB)
                ITIMAB = ITIMA*ITIMB
                IF    (ITIMAB.EQ. 1) THEN
                  WRITE(LUPRI,'(/A,I2,A,I2,A)')
     &             '<<A(',INDAP,'),B(',INDBP,
     &             ')>> - linear response function (real):'
                ELSEIF(ITIMAB.EQ.-1) THEN
                  WRITE(LUPRI,'(/A,I2,A,I2,A)')
     &             '<<A(',INDAP,'),B(',INDBP,
     &             ')>> - linear response function (imaginary):'
                ENDIF
                CALL PRSYMB(LUPRI,'-',95,0)
                WRITE(LUPRI,'(3X,A,A16,5X,A3,2X,A2)')
     &             'A - ',PRPNAM(INDPRA),REP(IREPA),TTYP(ITIMA),
     &             'B - ',PRPNAM(INDPRB),REP(IREPB),TTYP(ITIMB)
                CALL PRSYMB(LUPRI,'-',95,0)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C
C     **************************************************
C     *** Transfer information from XCBXLR to WRKXRS ***
C     **************************************************
C
      CALL SETRSP
      LINEQ     = .TRUE.
      LSVCFG(1) = XLR_LSFG(1)
      LSVCFG(2) = XLR_LSFG(2)
      TKNORM    = XLRNRM
      DIAGHE    = XLRDIH
      IPRXRS    = IPRXLR
      THCXRS    = THCLR
      RESFAC    = RESXLR
      MAXITR    = ITRXLR
      NFREQ     = NBFREQ
      MAXSIM    = MXLOAD
      NREDM     = MAXRM
      N2REDM    = MAXRM*MAXRM
      LOFFTY    = 0
      STATIC    = (NFREQ.EQ.1.AND.BFREQ(1).EQ.D0)
      UNCOUP    = XLR_UNCOUP
      IMFREQ    = XLR_IMFREQ
      TRIPLET   = XLR_TRIPLET
      CNVINT(1) = CNVXLR(1)
      CNVINT(2) = CNVXLR(2)
      ITRINT(1) = ITRILR(1)
      ITRINT(2) = ITRILR(2)
      INTDEF    = INTXLR
      STERNH    = XSTERN
      STERNC    = XSTERNC
      E2CHEK    = XLR_E2CHEK
      DAMPFREQ  = DMPFRLR
      NDAMP     = 1
      IF(DAMPFREQ.NE.D0) NDAMP = 2
C
C     Restart info used by MCSCF program.
C
      IXRSTB(-1) = 0
      IXRSTB(+1) = 0
      IXRSTB(+2) = 0
C
C     Indicate from vectors are read...
C
      LUBUF = 22
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','EIGEXT')
      CALL REACMO(LUBUF,'DFCOEF',DUM,DUM,DUM,DUM,1)
      CLOSE(LUBUF,STATUS='KEEP')
C
C     Initialize file for solution vectors
C
      LUXVC = 29
      CALL OPNFIL(LUXVC,'PAMXVC','UNKNOWN','PRPXLR')
      WRITE(LUXVC) 'END_OF_THIS_FILE'
C
C     Look for restart file for solution vectors
C
      RSREST = .FALSE.
      IF(XLR_XVCFIL.NE.'NORXVC') THEN
        INQUIRE(FILE=XLR_XVCFIL,EXIST=TOBE)
        RSREST = TOBE
        IF(RSREST) THEN
          LURST = 28
          CALL OPNFIL(LURST,XLR_XVCFIL,'OLD','PRPXLR')
        ELSE
          WRITE(LUPRI,*)
     &    '*** WARNING: Restart file '//XLR_XVCFIL//' not found !'
        ENDIF
      ENDIF

      CALL QEXIT('SETXLR')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETLRF */
      SUBROUTINE GETLRF(ABLRF,IBTYP,IBCVC,IBEVC,IBPVC,RCNV,
     &                  EVALR,EVECR,IORBCL,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Evaluate linear response functions <<A;B>>.
C
C     Written by T.Saue Oct 22 1996
C
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0,D0=0.0D0)
C
#include "dcbxrs.h"
#include "dcbxlr.h"
#include "dcbibn.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
      LOGICAL FIRST
      DIMENSION ABLRF(NFREQ,NLRAPT,NLRBPT,*),RCNV(*),
     &          IBTYP(2,*),IBCVC(*),IBEVC(*),IORBCL(4,*),
     &          IBPVC(*),EVALR(*),EVECR(NEVEC,*),WORK(*)
C
      CALL QENTER('GETLRF')
      KFRSAV = KFREE
      FIRST = .TRUE.
      NRSP = NFREQ*NLRAPT*NLRBPT*NDAMP
C
C     First evaluate solution vectors for all frequencies
C     and write them to direct access file
C     ===================================================
C
      NESIM = 0
      NPSIM = 0
      NCSIM = 0
C
C     Get coefficients
C     ----------------
C
      CALL MEMGET2('REAL','CMO',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      LUBUF = 22
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','GETLRF')
      IF(SPINFR) THEN
        CALL MEMGET2('INTE','IBEIG',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
        CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),DUM,WORK(KBEIG),
     &              TOTERG,10)
      ELSE
        KBEIG = KFREE
        CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),DUM,IDUM,TOTERG,2)
      ENDIF
      CLOSE(LUBUF,STATUS='KEEP')
C
C
C
      IF(XLRANA) THEN
        CALL MEMGET2('INTE','IORBCL',KLRO,4*NORBT,WORK,KFREE,LFREE)
        CALL DEFORB(WORK(KLRO),NCLS_XLR,XLR_INDSTR,
     &              WORK,KFREE,LFREE,IPRXLR)
      ELSE
        KLRO=KFREE
      ENDIF
C
C
C     Solution vector - orbital (e-e) part
C     ====================================
C
      IOFF  = 1*NDAMP
      IF(NZXOPE.GT.0) THEN
        NESIM = NFREQ * NSTAT
        KEE   = KFREE
        CALL MEMGET2('REAL','GPO'  ,KGVEC,NZXOPEQ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','XPO'  ,KXVEC,NZXOPEQ*NFREQ*NDAMP,
     &                            WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','BBUF' ,KBBUF,NZXOPEQ,
     &                              WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IVECS',KIVEC,NESIM,
     &                              WORK,KFREE,LFREE)
        IF(XLRANA) THEN
          NI = 0
          NV = 0
          DO IFRP = 1,NFSYM
C...inactive electron orbitals
            NI = NI + NCLS_XLR(2,1,IFRP)
C...secondary electron orbitals
            NV = NV + NCLS_XLR(2,3,IFRP)
          ENDDO
          NINDX=MAX(NI,NV,NZXOPE)
          CALL MEMGET2('INTE','INDX',KINDX,NINDX,WORK,KFREE,LFREE)
          CALL MEMGET2('INTE','IIBF',KII,NI,WORK,KFREE,LFREE)
          CALL MEMGET2('INTE','IVBF',KIV,NV,WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','RIBF',KRI,NI,WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','RVBF',KRV,NV,WORK,KFREE,LFREE)
        ELSE
          KINDX= KFREE
          KII  = KFREE
          KIV  = KFREE
          KRI  = KFREE
          KRV  = KFREE
        ENDIF
        CALL GETLR1(ABLRF(1,1,1,IOFF+1),JBENDX,EVECR,EVALR,RCNV,
     &              IBTYP,IBEVC,
     &              WORK(KGVEC),WORK(KXVEC),WORK(KBBUF),NZXOPE,
     &              WORK(KCMO),WORK(KBEIG),IORBCL,
     &              get_orbital_rotation_indices_pp(),
     &              WORK(KIVEC),NI,NV,
     &              WORK(KINDX),WORK(KII),WORK(KIV),WORK(KRI),WORK(KRV),
     &              WORK,KFREE,LFREE)
        IF(FIRST) THEN
          CALL DCOPY(NRSP,ABLRF(1,1,1,IOFF+1),1,ABLRF,1)
          FIRST = .FALSE.
        ELSE
          CALL DAXPY(NRSP,D1,ABLRF(1,1,1,IOFF+1),1,ABLRF,1)
        ENDIF
        CALL MEMREL('GETLRF.e-e',WORK,1,KEE,KFREE,LFREE)
        IOFF = IOFF + 1*NDAMP
      ENDIF
C
C     Solution vector - orbital (e-p) part
C     ====================================
C
      IF(NZXOPP.GT.0) THEN
        NPSIM = NFREQ * NSTAT
        KEP   = KFREE
        CALL MEMGET2('REAL','XPO'  ,KXVEC,NZXOPPQ*NFREQ*NDAMP,
     &                              WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','CPO',  KGVEC,NZXOPPQ,
     &                              WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','BBUF' ,KBBUF,NZXOPPQ
     &                              ,WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IVECS',KIVEC,NPSIM
     &                              ,WORK,KFREE,LFREE)
        IF(XLRANA) THEN
          NI = 0
          NV = 0
          DO IFRP = 1,NFSYM
C...inactive electron orbitals
            NI = NI + NCLS_XLR(2,1,IFRP)
C...secondary positron orbitals
            NV = NV + NCLS_XLR(1,3,IFRP)
          ENDDO
          NINDX=MAX(NI,NV,NZXOPP)
          CALL MEMGET2('INTE','INDX',KINDX,NINDX,WORK,KFREE,LFREE)
          CALL MEMGET2('INTE','IIBF',KII,NI,WORK,KFREE,LFREE)
          CALL MEMGET2('INTE','IVBF',KIV,NV,WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','RIBF',KRI,NI,WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','RVBF',KRV,NV,WORK,KFREE,LFREE)
        ELSE
          KINDX= KFREE
          KII  = KFREE
          KIV  = KFREE
          KRI  = KFREE
          KRV  = KFREE
        ENDIF
        CALL GETLR1(ABLRF(1,1,1,IOFF+1),JBPNDX,EVECR,EVALR,RCNV,
     &              IBTYP,IBPVC,
     &              WORK(KGVEC),WORK(KXVEC),WORK(KBBUF),NZXOPP,
     &              WORK(KCMO),WORK(KBEIG),IORBCL,
     &              get_orbital_rotation_indices_pn(),
     &              WORK(KIVEC),NI,NV,
     &              WORK(KINDX),WORK(KII),WORK(KIV),WORK(KRI),WORK(KRV),
     &              WORK,KFREE,LFREE)
        IF(FIRST) THEN
          CALL DCOPY(NRSP,ABLRF(1,1,1,IOFF+1),1,ABLRF,1)
          FIRST = .FALSE.
        ELSE
          CALL DAXPY(NRSP,D1,ABLRF(1,1,1,IOFF+1),1,ABLRF,1)
        ENDIF
        CALL MEMREL('GETLRF.e-p',WORK,1,KEP,KFREE,LFREE)
        IOFF = IOFF + 1*NDAMP
      ENDIF
C
C     Solution vector - configurational part
C     ======================================
C
      IF(NZCONF.GT.0) THEN
C     ..not yet written.....tsaue
      ENDIF
C
      CALL MEMREL('GETLRF.end',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('GETLRF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETLR1 */
      SUBROUTINE GETLR1(ABLRF,NBTYP,EVECR,EVALR,RCNV,IBTYP,IBVEC,
     &                  GPO,XPO,BBUF,NXPAR,CMO,IBEIG,IORBCL,JXOP,
     &                  IVECS,NI,NV,INDX,IIBF,IVBF,RIBF,RVBF,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get linear response function for orbital rotation parameter
C       contribution from given variational parameter type
C       specified by NBTYP
C
C     Written by T.Saue Oct 22 1996
C
C***********************************************************************
#ifdef MOD_LAO_REARRANGED
      use london_helper
      use london_utils 
#endif
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D2 = 2.0D00)
C
#include "dcbibn.h"
#include "dcbxrs.h"
#include "dcbxlr.h"
#include "dgroup.h"
#include "dcbxpr.h"
      CHARACTER FILNAM*6,MEMTYP*6,BLANK16*16,BLANK2*2
      DIMENSION ABLRF(NFREQ,NLRAPT,NLRBPT,NDAMP,*),RCNV(*),
     &          GPO(NXPAR,NZ),XPO(NXPAR,NZ,NDAMP,NFREQ),BBUF(*),
     &          CMO(*),IBEIG(*),IORBCL(4,*),
     &          IBVEC(*),IBTYP(2,*),EVECR(NEVEC,*),EVALR(*),
     &          JXOP(2,*),IVECS(*),IIBF(*),IVBF(*),INDX(*),
     &          RIBF(*),RVBF(*),WORK(*)
      CHARACTER TYP(3)*1
      DATA TYP/'i','a','v'/
      DATA BLANK16/'                '/
      DATA BLANK2/'  '/
C
#include "chrnos.h"
      CALL QENTER('GETLR1')
      KFRSAV = KFREE
      IF(IPRXLR.GE.2) CALL HEADER('Output from GETLR1',-1)
C
      NBDIM = NXPAR*NZ
      IF(NBTYP.EQ.JBENDX) THEN
        NBRED  = NERED
        LUBVEC = LUBOE
        FILNAM = 'PAMBOE'
        MEMTYP = 'GETLRE'
        ITYP   = IPEP
      ELSEIF(NBTYP.EQ.JBPNDX) THEN
        NBRED  = NPRED
        LUBVEC = LUBOP
        FILNAM = 'PAMBOP'
        MEMTYP = 'GETLRP'
        ITYP   = IPPP
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'GETLR1 ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('GETLR1: Unknown NBTYP !')
      ENDIF
C
      INDPRB = LLRBPU(INDBP)
      IF(INDPRB.GT.0) THEN
        IF(INDPRB.LT.10) THEN
          NDIG = 1
        ELSEIF(INDPRB.LT.100) THEN
          NDIG = 2
        ELSEIF(INDPRB.LT.1000) THEN
          NDIG = 3
        ELSE
          WRITE(LUPRI,'(A,I5)') 'GETLR1 ERROR: INDPRB = ',INDPRB
          CALL QUIT('GETLR1: Inappropriate INDPRB !')
        ENDIF
      ELSE
        WRITE(LUPRI,'(A,I5)') 'GETLR1 ERROR: INDPRB = ',INDPRB
        CALL QUIT('GETLR1: Inappropriate INDPRB !')
      ENDIF
C
C     Prepare for orbital analysis
C
      IF(XLRANA) THEN
        CALL IZERO(IIBF,NI)
        CALL IZERO(IVBF,NV)
C       Assign orbitals
        DO J = 1,NXPAR
          I1 = JXOP(1,J)
          J1 = IORBCL(4,I1)
          I2 = JXOP(2,J)
          J2 = IORBCL(4,I2)
          IIBF(J1) = I1
          IVBF(J2) = I2
        ENDDO
      ENDIF
C
C     First get solution vector
C     =========================
C
      OPEN(LUBVEC,FILE=FILNAM,FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NBDIM,STATUS='UNKNOWN')
      ITIM = IPX
      NXDIM = NBDIM*NFREQ*NDAMP
C
C     If static we have one round with ITIM = IPX, if
C     non-zero freq. then we have a second round with ITIM = IMX
C
      IROUND = 0
  100 IROUND = IROUND + 1
      CALL DZERO(XPO,NXDIM)
      CALL XRSXV1(ITIM,NBTYP,XPO,EVECR,NFREQ,IBTYP,IBVEC,IVECS,BBUF)
!gosia:
! if shielding_rearrange then here form: \kappa = \kappa - T:
#ifdef MOD_LAO_REARRANGED
      if (shielding_rearrange) then
        nlength=nxpar*ndamp*nfreq
        call combine_kappa_T(xpo, nlength, nbtyp, rsplab, itim)
      end if
#endif


      ! print section
      if (iprxlr > 2) then
         do ifreq = 1, nfreq
            WRITE(LUPRI,'(/1X,A2,A,I3)') BVTYP(ITYP),
     &     ' part of solution vectors for frequency', ifreq
            do idamp = 1, ndamp
               if (ndamp == 2) then
                  ! we only print this for complex response vectors
                  if (idamp == 1) then
                     write(lupri, '(/1x,a)')                            &
     &   'Contribution to the real part of the complex response vector'
                  else
                     write(lupri, '(/1x,a)')                            &
     &   'Contribution to the imag part of the complex response vector'
                  end if
               end if
               call prbvec(lupri, xpo(1, 1, idamp, ifreq), 1, NXPAR)
            end do
         end do
      end if


C     Next form linear response functions
C     ===================================
C     (for now we form only combinations between
C      operators with same point group symmetry)
C
      IF(ALLCMB) THEN
        WRITE(LUPRI,'(A)')
     &  'GETLR1 WARNING: Sorry, ALLCMB not yet written !'
        ALLCMB = .FALSE.
      ENDIF
      IF(.NOT.ALLCMB) THEN
        NOPA   = NLRAP(JSYMOP)
        JOPA   = JLRAP(JSYMOP)
        DO IOPA = 1,NOPA
          INDAP  = LLAPSU(JOPA+IOPA)
          INDPRA = LLRAPU(INDAP)
          KTIMOP = IPRPTIM(INDPRA)
C
C         Get gradient
C         ------------
C
          IF (KTIMOP .EQ. ITIM) THEN
             CALL GPGET(INDPRA,GPO,JXOP,NXPAR,CMO,IBEIG,
     &                  NBTYP,.TRUE.,WORK,WORK,KFREE,LFREE,IPRXLR)
          END IF
C
C         Form linear reponse function
C         ----------------------------
C
          DO IFREQ = 1,NFREQ
            IF (KTIMOP .EQ. ITIM) THEN
              DO IDAMP = 1,NDAMP
                ABLRF(IFREQ,INDAP,INDBP,IDAMP,1)  =
     &          D2*DDOT(NBDIM,GPO,1,XPO(1,1,IDAMP,IFREQ),1)
                IF(XLRANA.OR.IPRXLR.GE.2) THEN
                  FREQ = EVALR(IFREQ)
                  WRITE(LUPRI,1000)
     &            '<<',INDAP,',',INDBP,'>>:',
     &            ABLRF(IFREQ,INDAP,INDBP,IDAMP,1),
     &            ' a.u.,  ',BVTYP(ITYP),
     &            '-contribution at frequency',FREQ,' a.u.'
                ENDIF
                IF(XLRANA) THEN
                  CALL PRSYMB(LUPRI,'-',95,0)
                  WRITE(LUPRI,'(/A/)') '* ORBITAL ANALYSIS'
                  CALL DZERO(BBUF,NXPAR)
                  DO IZ = 1,NZ
                    DO I = 1,NXPAR
                      BBUF(I) = BBUF(I)+GPO(I,IZ)*XPO(I,IZ,IDAMP,IFREQ)
                    ENDDO
                  ENDDO
                  CALL DSCAL(NXPAR,D2,BBUF,1)
                  TMP = 100.0D0/ABLRF(IFREQ,INDAP,INDBP,IDAMP,1)
                  CALL DZERO(RIBF,NI)
                  CALL DZERO(RVBF,NV)
                  DO J = 1,NXPAR
                    I1 = JXOP(1,J)
                    J1 = IORBCL(4,I1)
                    I2 = JXOP(2,J)
                    J2 = IORBCL(4,I2)
                    RIBF(J1) = RIBF(J1) + BBUF(J)
                    RVBF(J2) = RVBF(J2) + BBUF(J)
                  ENDDO
C.......Give major inactive orbitals
        WRITE(LUPRI,'(3X,A)') 'Dominant inactive orbitals:'
        CALL ABSIND (NI,RIBF,INDX)
        DO J = NI,1,-1
          JJ = INDX(J)
          FAC = RIBF(JJ)*TMP
          IF(ABS(FAC).LE.XLR_ANATHR) GOTO 3
          I1  = IIBF(JJ)
          J1 = IORBCL(2,I1)
          WRITE(LUPRI,'(I6,A1,A3,A1,E16.8,3X,F8.2,A1)')
     &        J1,'(',FREP(IORBCL(3,I1)),')',
     &        RIBF(JJ),FAC,'%'
        ENDDO
 3      CONTINUE
C       Give major virtual orbitals
        WRITE(LUPRI,'(3X,A)') 'Dominant virtual orbitals:'
        CALL ABSIND (NV,RVBF,INDX)
        DO J = NV,1,-1
          JJ = INDX(J)
          FAC = RVBF(JJ)*TMP
          IF(ABS(FAC).LE.XLR_ANATHR) GOTO 4
          I2 = IVBF(JJ)
          J2 = IORBCL(2,I2)
          WRITE(LUPRI,'(I6,A1,A3,A1,E16.8,3X,F8.2,A1)')
     &        J2,'(',FREP(IORBCL(3,I2)),')',
     &        RVBF(JJ),FAC,'%'
        ENDDO
 4      CONTINUE
        CALL ABSIND (NXPAR,BBUF,INDX)
        WRITE(LUPRI,'(3X,A,E10.2/)')
     &   'Amplitudes larger than threshold ',XLR_ANATHR
        DO J = NXPAR,1,-1
          JJ = INDX(J)
          FAC = BBUF(JJ)*TMP
          IF(ABS(FAC).LE.XLR_ANATHR) GOTO 5
          I1 = JXOP(1,JJ)
          I2 = JXOP(2,JJ)
          WRITE(LUPRI,
     &    '(3X,I6,3A1,A3,A1,A,I6,3A1,A3,A1,E16.8,3X,F8.2,A1)')
     &      IORBCL(2,I1),'(',TYP(IORBCL(1,I1)),':',
     &      FREP(IORBCL(3,I1)),')',
     &      ' ---> ',
     &      IORBCL(2,I2),'(',TYP(IORBCL(1,I2)),':',
     &      FREP(IORBCL(3,I2)),')',BBUF(JJ),FAC,'%'
        ENDDO
 5      CONTINUE
            CALL PRSYMB(LUPRI,'-',95,0)
              ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C
C     Write solution vectors to unformatted file for restart
C
CPP...we only write real part to file for now
      DO IFREQ = 1,NFREQ
        FREQ = EVALR(IFREQ)
        IF (IROUND.EQ.1 .OR. FREQ.NE.D0) THEN
          CALL WRTRSP(LUXVC,XPO(1,1,1,IFREQ),JXOP,
     &         RSPLAB,BVTYP(ITYP),FREQ,JSYMOP,JTIMOP,
     &         BLANK16,BLANK2,D0,0,0,RCNV(IFREQ),INTXRS,NBDIM)
        END IF
      ENDDO
C
C     IF not static, then now calculate the part of the
C     B solution vector which is of -JTIMOP symmetry,
C     calculate contribution to properties and write them to file also.
C
      IF(.NOT.STATIC .AND. IROUND.EQ.1) THEN
        ITIM = IMX
        ITYP = ITYP + 1
        GO TO 100
      END IF
      CLOSE(LUBVEC,STATUS = 'KEEP')
C
      CALL QEXIT('GETLR1')
      RETURN
 1000 FORMAT(/A,I2,A,I2,A,E16.8,A,3X,A2,A,F16.10,A)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XLROUT */
      SUBROUTINE XLROUT(ABLRF,RCNV,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Final output from linear response module
C
C     Written by T.Saue Oct 31 1996
C
C***********************************************************************
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
C
      PARAMETER(D0 = 0.0D0)
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxlr.h"
#include "dcbxrs.h"
#include "dcbxpr.h"
C
      CHARACTER MXFORM*6,FMT*6,FMT2*6,TTYP(-1:1)*2,FREQLAB*4,
     &          RESPLAB(3)*20
      DIMENSION ABLRF(NBFREQ,NLRAPT,NLRBPT,NDAMP,*),
     &          RCNV(NBFREQ,NLRBPT),IBUF(2),WORK(*)
      DATA TTYP /'T-','T0','T+'/
      DATA RESPLAB/'Real part           ',
     &             'Imaginary part      ',
     &             '                    '/
C
C     Loop over perturbation operators
C     ================================
C
      FREQLAB = 'real'
      IF(IMFREQ) FREQLAB = 'imag'
      CALL HEADER('Linear response; final output ',-1)
      IF(.NOT.ALLCMB) THEN
        DO IOPSY = 1,NFSYM
          NBF = NBSYM/NFSYM
          DO IBF = 1,NBF
            ISYM = JFSYM(IBF,IOPSY)
            NOPA = NLRAP(ISYM)
            JOPA = JLRAP(ISYM)
            NOPB = NLRBP(ISYM)
            JOPB = JLRBP(ISYM)
            DO IOPA = 1,NOPA
              INDAP  = LLAPSU(JOPA+IOPA)
              INDPRA = LLRAPU(INDAP)
              IREPA  = IPRPSYM(INDPRA)-1
              ITIMA  = IPRPTIM(INDPRA)
              DO IOPB = 1,NOPB
                INDBP  = LLBPSU(JOPB+IOPB)
                INDPRB = LLRBPU(INDBP)
                IREPB  = IPRPSYM(INDPRB)-1
                ITIMB  = IPRPTIM(INDPRB)
                ITIMAB = ITIMA*ITIMB
                IF(DAMPFREQ.NE.D0) THEN
                  WRITE(LUPRI,'(/A,I2,A,I2,A)')
     &               '<<A(',INDAP,'),B(',INDBP,
     &               ')>> - linear response function (complex):'
                  J1 = 1
                  J2 = 2
                  J3 = 1
                ELSE
                  IF    (ITIMAB.EQ. 1) THEN
                    WRITE(LUPRI,'(/A,I2,A,I2,A)')
     &               '<<A(',INDAP,'),B(',INDBP,
     &               ')>> - linear response function (real):'
                    J1 = 1
                    J2 = 3
                    J3 = 2
                  ELSEIF(ITIMAB.EQ.-1) THEN
                    WRITE(LUPRI,'(/A,I2,A,I2,A)')
     &               '<<A(',INDAP,'),B(',INDBP,
     &               ')>> - linear response function (imaginary):'
                    J1 = 2
                    J2 = 3
                    J3 = 1
                  ENDIF
                ENDIF
                CALL PRSYMB(LUPRI,'-',95,0)
                WRITE(LUPRI,'(3X,A,A16,5X,A3,2X,A2)')
     &             'A - ',PRPNAM(INDPRA),REP(IREPA),TTYP(ITIMA),
     &             'B - ',PRPNAM(INDPRB),REP(IREPB),TTYP(ITIMB)
                CALL PRSYMB(LUPRI,'-',95,0)
                WRITE(LUPRI,'(1X,A,A4,A1,5X,2(A,3X),A)')
     &           'Frequency (',FREQLAB,')',(RESPLAB(J),J=J1,J2,J3),
     &           'Convergence'
                CALL PRSYMB(LUPRI,'-',95,0)
                DO I = 1,NBFREQ
                  DO J = J1,J2
                  ENDDO
                  IF(DAMPFREQ.NE.D0) THEN
                    FMT = MXFORM(ABLRF(I,INDAP,INDBP,1,1),15)
                    FMT2 = MXFORM(ABLRF(I,INDAP,INDBP,2,1),15)
                    IF(UNCOUP) THEN
                      WRITE(LUPRI,
     &    '(F12.8,A,3X,'//FMT//',A,3X,'//FMT2//',A,1X,E9.2,1X,A)')
     &                    BFREQ(I),' a.u.',ABLRF(I,INDAP,INDBP,1,1),
     &                    ' a.u.',ABLRF(I,INDAP,INDBP,2,1),
     &                    ' a.u.',RCNV(I,INDBP),
     &                    '--- (uncoupled calculation)'
                    ELSEIF(RCNV(I,INDBP).LT.THCLR) THEN
                      WRITE(LUPRI,
     &    '(F12.8,A,3X,'//FMT//',A,3X,'//FMT2//',A,1X,E9.2,1X,A)')
     &                    BFREQ(I),' a.u.',ABLRF(I,INDAP,INDBP,1,1),
     &                    ' a.u.',ABLRF(I,INDAP,INDBP,2,1),
     &                    ' a.u.',RCNV(I,INDBP),'(converged)'
                    ELSE
                      WRITE(LUPRI,
     &    '(F12.8,A,3X,'//FMT//',A,3X,'//FMT2//',A,1X,E9.2,1X,A)')
     &                    BFREQ(I),' a.u.',ABLRF(I,INDAP,INDBP,1,1),
     &                    ' a.u.',ABLRF(I,INDAP,INDBP,2,1),
     &                    ' a.u.',RCNV(I,INDBP),'(**NOT** converged!)'
                    ENDIF
                  ELSE
                      FMT = MXFORM(ABLRF(I,INDAP,INDBP,1,1),15)
                      IF(UNCOUP) THEN
                        WRITE(LUPRI,
     &                   '(F12.8,A,3X,'//FMT//',A,22X,1P,E9.2,3X,A)')
     &                    BFREQ(I),' a.u.',ABLRF(I,INDAP,INDBP,1,1),
     &                    ' a.u.',RCNV(I,INDBP),
     &                    '--- (uncoupled calculation)'
                      ELSEIF(RCNV(I,INDBP).LT.THCLR) THEN
                        WRITE(LUPRI,
     &                   '(F12.8,A,3X,'//FMT//',A,22X,1P,E9.2,3X,A)')
     &                    BFREQ(I),' a.u.',ABLRF(I,INDAP,INDBP,1,1),
     &                    ' a.u.',RCNV(I,INDBP),'(converged)'
                      ELSE
                        WRITE(LUPRI,
     &                   '(F12.8,A,3X,'//FMT//',A,22X,1P,E9.2,3X,A)')
     &                    BFREQ(I),' a.u.',ABLRF(I,INDAP,INDBP,1,1),
     &                    ' a.u.',RCNV(I,INDBP),'(**NOT** converged!)'
                      ENDIF
                  ENDIF
                ENDDO
                CALL PRSYMB(LUPRI,'-',95,0)
              ENDDO
            ENDDO
            IF(DAMPFREQ.NE.D0) THEN
              DFREQEV = DAMPFREQ * XTEV
              DFREQCM = DAMPFREQ * XTKAYS
              WRITE(LUPRI,'(3X,A,F10.6,A,F10.6,A,F10.1,A)')
     &             'Damping parameter equals',DAMPFREQ,' a.u. =',
     &              DFREQEV,' eV =',DFREQCM,' cm-1'
            ENDIF
          ENDDO
        ENDDO
      ELSE
        CALL QUIT('XLROUT: No code for ALLCMB !')
      ENDIF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Anaxrs */
      SUBROUTINE ANAXRS(WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Analyze solution vectors on file PAMXVC
C
C     Written by Trond Saue July 4 2001
C
C***********************************************************************

      use overlap_diagnostic
      use dft_cfg
      use orbital_rotation_indices


#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcborb.h"
#include "dcbbas.h"
C
      LOGICAL   TOBE,GASWF
      DIMENSION WORK(*)
      INTEGER   NCLS(2,3,2)
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2,LBUF*16,PBUF*1
      integer, allocatable :: norot(:)
C
      CALL QENTER('ANAXRS')
      KFRSAV = KFREE
      CALL HEADER('Analysis of response solution vectors',-1)

#ifdef MOD_UNRELEASED
      if (dft_cfg_overlap_diagnostic) then
         call o_matrix_from_disc()
      end if
#endif

C
C     Check if file of solution vectors is present
C
      INQUIRE(FILE='PAMXVC',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
      LUINP = 1
      CALL OPNFIL(LUINP,'PAMXVC','OLD','ANAXRS')
      IREC = 0
      IERR = 0
      IVEC = 0
      IPAR = 0
      NDIM = 0
      LBUF = '                '
      PBUF = ' '
      TOBE = .FALSE.
      JFRP = -1
C     Set orbital strings and classify orbitals
      CALL OSTRING(INDSTR,SKIPEE,SKIPEP,IPRINT,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IORBCL',KLRO,4*NORBT,WORK,KFREE,LFREE)
      CALL DEFORB(WORK(KLRO),NCLS,INDSTR,
     &              WORK,KFREE,LFREE,IPRINT)
      KLOOP=KFREE
C
C     WHILE DO - loop
C     1) find how many vectors of same LAB1 and PTYPE
C     2) then analyze these vectors
C     3) if not 'END_OF_THIS_FILE' go back to 1)
C
    1 CONTINUE
        READ (LUINP,END=3,ERR=2) LAB1
        IF(LAB1.EQ.'END_OF_THIS_FILE') THEN
          IF(TOBE) THEN
C           Analyze NVEC solution vectors from last parameter class
            NVEC = IVEC+1-IPAR
            NBDIM = NDIM/NZ
            IF(PBUF.EQ.'E') THEN
              IF(SKIPEE) GOTO 20
              IC = 1
              CALL ANAXR1(NBDIM,NVEC,IC,
     &                    get_orbital_rotation_indices_pp(),WORK(KLRO),
     &                    LUINP,LAB1,LAB2,TYP1,TYP2,
     &                    NCLS,WORK,KFREE,LFREE,IPRINT)
            ELSEIF(PBUF.EQ.'P') THEN
              IF(SKIPEP) GOTO 20
              IC = 2
              CALL ANAXR1(NBDIM,NVEC,IC,
     &                    get_orbital_rotation_indices_pn(),WORK(KLRO),
     &                    LUINP,LAB1,LAB2,TYP1,TYP2,
     &                    NCLS,WORK,KFREE,LFREE,IPRINT)
            ENDIF
          ENDIF
          GOTO 20
        ENDIF
        BACKSPACE LUINP
        READ (LUINP,END=3,ERR=2) LAB1,LAB2,TYP1,TYP2,
     &            FREQ1,FREQ2,JSYMOP1,JSYMOP2,
     &            JTIMOP1,JTIMOP2,RNORM,LEN,INTFLG,
     &            ERGRSP,NBAS,NORBS
        IVEC = IVEC + 1
        IREC = IREC + 1
C       New parameter class
        IF(TYP1(1:1).NE.PBUF.OR.LAB1.NE.LBUF) THEN
          IF(TOBE) THEN
C           Analyze NVEC solution vectors from previous parameter class
            NVEC = IVEC-IPAR
            NBDIM = NDIM/NZ
            IF(PBUF.EQ.'E') THEN
              IF(SKIPEE) GOTO 30
              IC = 1
              CALL ANAXR1(NBDIM,NVEC,IC,
     &                    get_orbital_rotation_indices_pp(),WORK(KLRO),
     &                    LUINP,LAB1,LAB2,TYP1,TYP2,
     &                    NCLS,WORK,KFREE,LFREE,IPRINT)
            ELSEIF(PBUF.EQ.'P') THEN
              IF(SKIPEP) GOTO 30
              IC = 2
              CALL ANAXR1(NBDIM,NVEC,IC,
     &                    get_orbital_rotation_indices_pn(),WORK(KLRO),
     &                    LUINP,LAB1,LAB2,TYP1,TYP2,
     &                    NCLS,WORK,KFREE,LFREE,IPRINT)
            ENDIF
            READ (LUINP,END=3,ERR=2) LAB1,LAB2,TYP1,TYP2,
     &              FREQ1,FREQ2,JSYMOP1,JSYMOP2,
     &              JTIMOP1,JTIMOP2,RNORM,LEN,INTFLG,
     &              ERGRSP,NBAS,NORBS
          ENDIF
 30       CONTINUE
C         collect type, first vector, and length for this parameter 
C         class
          PBUF = TYP1(1:1)
          IPAR = IVEC
          NDIM = LEN
          TOBE = .TRUE.
        ENDIF
C       New property
        IF(LAB1.NE.LBUF) THEN
          IREP = JSYMOP1-1
          IFRP = JBTOF(IREP,1)
C         Generate orbital rotation operator
          IF(IFRP.NE.JFRP) THEN
            IF(TOBE) THEN
              CALL MEMREL('ANAXRS',WORK,KLOOP,KLOOP,KFREE,LFREE)
            ENDIF
            allocate(norot(norbt))
            norot = 0
            GASWF = .FALSE.
            CALL SETXOP(INDSTR,SKIPEE,SKIPEP,GASWF,IFRP,NOROT,
     &            NZXOPE,NZHOPE,NZXOPP,
     &            IPRINT)
            deallocate(norot)
            JFRP = IFRP
          ENDIF
C
          LBUF = LAB1
        ENDIF
C.......skip solution vector
        READ(LUINP,END=3,ERR=2)
C.......skip orbital rotations
        READ(LUINP,END=3,ERR=2)
        IREC = IREC + 1
      GO TO 1
C
    2 CONTINUE
      IREC = IREC + 1
      IERR = IERR + 1
      WRITE (LUPRI, '(/A,I5/)') ' ERROR READING RECORD NO.',IREC
      GOTO 20
C
    3 CONTINUE
      WRITE (LUPRI,'(/I10,A)') IREC,
     *   ' records read before EOF on file.'
      GOTO 20
 20   CONTINUE
      CLOSE(LUINP)
      CALL MEMREL('ANAXRS',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('ANAXRS')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(/A)') 'ANAXRS ERROR: No solution vectors found !'
      CALL QUIT('ANAXRS: Solution vectors not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck anaxr1 */
      SUBROUTINE ANAXR1(NBDIM,NVEC,IC,JXOP,IORBCL,
     &                  LUINP,LAB1,LAB2,TYP1,TYP2,
     &                  NCLS,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxlr.h"
      INTEGER   NCLS(2,3,2)
      DIMENSION JXOP(*),IORBCL(4,*),WORK(*)
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2
C
      KFRSAV = KFREE
      CALL MEMGET('REAL',KVEC ,NBDIM*NZ*NVEC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFREQ,NVEC         ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KNORM,NVEC         ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF1,NBDIM*NZ     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF2,NBDIM        ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KBTYP,NVEC         ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFLAG,NVEC         ,WORK,KFREE,LFREE)
      NI = 0
      NV = 0
      IF(IC.EQ.1) THEN
        DO IFRP = 1,NFSYM
C...inactive electron orbitals
          NI = NI + NCLS(2,1,IFRP)
C...secondary electron orbitals
          NV = NV + NCLS(2,3,IFRP)
        ENDDO
      ELSEIF(IC.EQ.2) THEN
        DO IFRP = 1,NFSYM
C...inactive positron orbitals
          NI = NI + NCLS(2,1,IFRP)
C...secondary positron orbitals
          NV = NV + NCLS(1,3,IFRP)
        ENDDO
      ELSE
        CALL QUIT('ANAXR1: Unknown IC !')
      ENDIF
      NINDX=MAX(NI,NV,NBDIM)
      CALL MEMGET('INTE',KINDX,NINDX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KII,NI,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIV,NV,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRI,NI,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRV,NV,WORK,KFREE,LFREE)
      CALL ANAXR2(WORK(KVEC),WORK(KFREQ),WORK(KNORM),
     &            WORK(KBTYP),WORK(KFLAG),WORK(KBUF1),
     &            WORK(KBUF2),WORK(KINDX),JXOP,IORBCL,
     &            NI,NV,WORK(KII),WORK(KIV),WORK(KRI),WORK(KRV),
     &            NBDIM,NVEC,LUINP,LAB1,LAB2,TYP1,TYP2,IPRINT)
      CALL MEMREL('ANAXR1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Anaxr2 */
      SUBROUTINE ANAXR2(BVEC,FREQ,RNORM,IBTYP,LFLAG,BUF1,BUF2,INDX,
     &        JXOP,IORBCL,NI,NV,IIBF,IVBF,RIBF,RVBF,
     &        NBDIM,NVEC,LUINP,LAB1,LAB2,TYP1,TYP2,IPRINT)
C***********************************************************************
C
C     Written by T.Saue July 6 2001
C
C***********************************************************************

      use overlap_diagnostic
      use dft_cfg

#include "implicit.h"
#include "priunit.h"
      PARAMETER(DTOL=1.0D-12,DTOL2=1.0D-4,D0=0.0D0,D1=1.0D0,DM=1.0D2)
C
#include "dgroup.h"
#include "pgroup.h"
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2
      LOGICAL LFLAG(NVEC)
      DIMENSION BVEC(NBDIM,NZ,NVEC),FREQ(NVEC),RNORM(NVEC),IBTYP(NVEC),
     &          BUF1(NBDIM,NZ),BUF2(NBDIM),INDX(*),JXOP(2,NBDIM),
     &          IORBCL(4,*),IIBF(NI),IVBF(NV),RIBF(NI),RVBF(NV)
C     Local variables
      CHARACTER TYP(3)*1
      DATA TYP/'i','a','v'/
      real(8) :: lambda_numerator
      real(8) :: lambda_denominator
C
      CALL IZERO(IIBF,NI)
      CALL IZERO(IVBF,NV)
C     Assign orbitals
      DO J = 1,NBDIM
        I1 = JXOP(1,J)
        J1 = IORBCL(4,I1)
        I2 = JXOP(2,J)
        J2 = IORBCL(4,I2)
        IIBF(J1) = I1
        IVBF(J2) = I2
      ENDDO
C
      BACKSPACE LUINP
      DO I = 1,NVEC
        BACKSPACE LUINP
        BACKSPACE LUINP
        BACKSPACE LUINP
      ENDDO
C     Read vectors
      READ (LUINP) LAB1,LAB2,TYP1,TYP2,
     &        FREQ1,FREQ2,JSYMOP1,JSYMOP2,
     &        JTIMOP1,JTIMOP2,DUM,LEN,INTFLG,
     &        ERGRSP,NBAS,NORBS
      IREP = JSYMOP1 - 1
      WRITE(LUPRI,'(A,A1,A,A16,A,A3,A,I3,A,I10)')
     &     '** ',TYP1(1:1),' solution vectors : ',LAB1,
     &     ' Irrep: ',REP(IREP),' Trev: ',JTIMOP1,
     &     ' Length: ',NBDIM
      BACKSPACE LUINP
      DO I = 1,NVEC
        READ (LUINP) LAB1,LAB2,TYP1,TYP2,
     &        FREQ(I),FREQ2,JSYMOP1,JSYMOP2,
     &        JTIMOP1,JTIMOP2,RNORM(I),LEN,INTFLG,
     &        ERGRSP,NBAS,NORBS
        IF(TYP1(2:2).EQ.'+') THEN
          IBTYP(I) = 1
        ELSEIF(TYP1(2:2).EQ.'-') THEN
          IBTYP(I) = -1
        ELSE
          CALL QUIT('XLRANA: Unknown IBTYP !')
        ENDIF
        CALL READT(LUINP,LEN,BVEC(1,1,I))
        IF(IPRINT.GE.4) THEN
          WRITE(LUPRI,'(3X,A2,A,F15.8)')
     &      TYP1,' vector at frequency: ', FREQ(I)
          CALL PRQMAT(BVEC(1,1,I),NBDIM,1,NBDIM,1,
     &                NZ,IPQTOQ(1,IREP),LUPRI)
        ENDIF
C...    skip orbital rotations
        READ(LUINP)
      ENDDO
C     Analyze major excitations
      CALL LSET(NVEC,.TRUE.,LFLAG)
      DO I = 1,NVEC
      IF(LFLAG(I)) THEN
        LFLAG(I) = .FALSE.
        II = IBTYP(I)
        DO IZ = 1,NZ
          CALL DCOPY(NBDIM,BVEC(1,IZ,I),1,BUF1(1,IZ),1)
        ENDDO
        IF(DABS(FREQ(I)-D0).GT.DTOL) THEN
C       Frequency dependent case; search for partner
          DO J = 1,NVEC
          IF(LFLAG(J)) THEN
            JJ = IBTYP(J)
            IF(DABS(FREQ(I)-FREQ(J)).LT.DTOL.AND.(JJ.EQ.-II)) THEN
              LFLAG(J)=.FALSE.
              DO IZ = 1,NZ
                CALL DAXPY(NBDIM,D1,BVEC(1,IZ,J),1,BUF1(1,IZ),1)
              ENDDO
              GOTO 10
            ENDIF
          ENDIF
          ENDDO
          WRITE(LUPRI,'(A,I5)')
     &      'ANAXR2 WARNING ! No partner found for vector. ',I
        ENDIF
 10     CONTINUE
C       Calculate norm
        D2TOT = D0
        DO J = 1,NBDIM
          BUF2(J) = BUF1(J,1)*BUF1(J,1)
        ENDDO
        DO IZ = 2,NZ
          DO J = 1,NBDIM
            BUF2(J) = BUF2(J) + BUF1(J,IZ)*BUF1(J,IZ)
          ENDDO
        ENDDO
        DO J = 1,NBDIM
          D2TOT = D2TOT + BUF2(J)
        ENDDO
        DNTOT = SQRT(D2TOT)
        WRITE(LUPRI,'(A8,F15.8,A,E16.8,A,E10.2)')
     &        ' Freq.: ',FREQ(I),' Norm: ',DNTOT,
     &        ' Residual norm: ',RNORM(I)
        IF(DNTOT.LT.DTOL2) THEN
          WRITE(LUPRI,'(A)') '--> Norm too small. Skipping.'
          GOTO 20
        ENDIF
C       Distribute squared norm on orbitals
        THR = 1.0D0
        CALL DZERO(RIBF,NI)
        CALL DZERO(RVBF,NV)
        DO J = 1,NBDIM
          I1 = JXOP(1,J)
          J1 = IORBCL(4,I1)
          I2 = JXOP(2,J)
          J2 = IORBCL(4,I2)
          RIBF(J1) = RIBF(J1) + BUF2(J)
          RVBF(J2) = RVBF(J2) + BUF2(J)
        ENDDO
        FAC = DM/D2TOT
        CALL DSCAL(NI,FAC,RIBF,1)
        CALL DSCAL(NV,FAC,RVBF,1)
C       Give major inactive orbitals
        WRITE(LUPRI,'(3X,A)') 'Dominant inactive orbitals:'
        CALL INDEXX (NI,RIBF,INDX)
        DO J = NI,1,-1
          JJ = INDX(J)
          IF(RIBF(JJ).LT.THR) GOTO 3
          I1  = IIBF(JJ)
          J1 = IORBCL(2,I1)
          WRITE(LUPRI,'(I6,A1,A3,A1,F6.2,A1)')
     &        J1,'(',FREP(IORBCL(3,I1)),')',RIBF(JJ),'%'
        ENDDO
 3      CONTINUE
C       Give major virtual orbitals
        WRITE(LUPRI,'(3X,A)') 'Dominant virtual orbitals:'
        CALL INDEXX (NV,RVBF,INDX)
        DO J = NV,1,-1
          JJ = INDX(J)
          IF(RVBF(JJ).LT.THR) GOTO 4
          I2  = IVBF(JJ)
          J2 = IORBCL(2,I2)
          WRITE(LUPRI,'(I6,A1,A3,A1,F6.2,A1)')
     &        J2,'(',FREP(IORBCL(3,I2)),')',RVBF(JJ),'%'
        ENDDO
 4      CONTINUE
C       Give major orbital rotations in excitations
        CALL INDEXX (NBDIM,BUF2,INDX)
        WRITE(LUPRI,'(3X,A,F6.2,A1/)')
     &   'Excitation amplitudes larger than weight ',THR,'%'
        DO J = NBDIM,1,-1
          JJ = INDX(J)
          WT = FAC*BUF2(JJ)
          IF(WT.LT.THR) GOTO 5
          I1 = JXOP(1,JJ)
          I2 = JXOP(2,JJ)
          WRITE(LUPRI,'(3X,I6,3A1,A3,A1,A,I6,3A1,A3,A1,F6.2,A1)')
     &      IORBCL(2,I1),'(',TYP(IORBCL(1,I1)),':',
     &      FREP(IORBCL(3,I1)),')',
     &      ' ---> ',
     &      IORBCL(2,I2),'(',TYP(IORBCL(1,I2)),':',
     &      FREP(IORBCL(3,I2)),')',WT,'%'
        ENDDO
 5      CONTINUE
        IF(IPRINT.GE.4) THEN
          WRITE(LUPRI,'(3X,A/)') 'All excitation amplitudes :'
          DO J = 1,NBDIM
            WRITE(LUPRI,'(3X,I6,E16.8)') J,BUF2(J)
          ENDDO
        ENDIF
#ifdef MOD_UNRELEASED
        if (dft_cfg_overlap_diagnostic) then
           lambda_numerator   = 0.0d0
           lambda_denominator = 0.0d0
           do j = nbdim, 1, -1
              jj = indx(j)
              i1 = jxop(1, jj)
              i2 = jxop(2, jj)
              call sum_lambda(lambda_numerator,  
     &                        lambda_denominator,
     &                        iorbcl(2, i1),            
     &                        iorbcl(2, i2),        
     &                        iorbcl(3, i1),
     &                        iorbcl(3, i2),
     &                        buf2(jj))
           enddo
           write(lupri, '(a)') '   PBHT MO Overlap Diagnostic'
           write(lupri, '(a)') '     ref.: M. J. G. Peach, '//
     &        'P. Benfield, T. Helgaker, and D. J. Tozer'
           write(lupri, '(a44)') 'J. Chem. Phys. 128, 044118 (2008)'
           write(lupri, '(a, f8.5)') '     lambda = ', 
     &        lambda_numerator/lambda_denominator
        end if
#endif /* MOD_UNRELEASED */
      ENDIF
 20   CONTINUE
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRINT_GP_XP(GP,XP,JXOP,NXOP,NZ)
C***********************************************************************
C
C     Printout of proprty gradient and solution vector
C     Written by Trond Saue Nov 26 2019
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION GP(NXOP,NZ),XP(NXOP,NZ),JXOP(2,NXOP)
      CHARACTER M(4)*6
      M(1) = 'A real'
      M(2) = 'A imag'
      M(3) = 'B real'
      M(4) = 'B imag'
      DO IZ = 1,NZ
        WRITE(LUPRI,'(/A,A6,A/)') '*** ',M(IZ),' part ***';
        WRITE(LUPRI,'(A)')
     & '    A    I   Property gradient     Solution vector'
        DO IXOP = 1,NXOP
          WRITE(LUPRI,'(2I5,2ES20.8)')
     &         JXOP(1,IXOP),JXOP(2,IXOP),GP(IXOP,IZ),XP(IXOP,IZ)
        ENDDO
      ENDDO
      RETURN
      END
