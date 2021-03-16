!
! FILE    : dcbopt.h
!
!*** 2nd order opt. ***
!     IPROPT    - print level
!     IOPT_ITRINT   - LS/SS integrals from iteration ITRINT(1)/ITRINT(2)
!     IOPT_INTDEF    - default INTFLG
!     IOPT_INTFLG    - LL/SL/SS integrals (packed in bits 0-2)
!     IOPT_INTBUF    - value of INTFLG from last macro iteration
!     OPT_INTTYP     - string with integral types
!     IOPT_IHMROP    - Hermitian symmetry of Fock matrices
!     IOPT_ISYMOP    - Symmetry of Fock matrices
!     IOPT_IFCKOP    - Type of Fock matrices
!     
!     NZVAR     - total number of variables
!     NZCONF    - number of configurational parameters
!     NZXOPE    - number of orbital electron-electron parameters
!     NZXOPP    - number of orbital electron-positron parameters
!     NZXOPT    - total number of orbital parameters
!     JOPSY     - symmetry of wave-function
!
!     OPT_INDSTR - inactive/active/secondary orbitals
!     
!     MOPT_MXMICRO   - maximum micro iterations
!     MOPT_MXMACRO   - maximum macro iterations
!
!     IOPT_STATE     - converge to state no. ISTATE (0=all NCIROOT roots, 1=ground state, etc)
!
!     FLAG:        (1) Always Newton-Raphson iterations
!                  (2) Always NEO iterations
!                 (19) Tight step contol for ground state also
!                 (20) Skip extra termination tests
!
!     IDENSLR_STATE  - natural orbital occupation numbers for state xxx
!     ISHMEM_TYPE    - type of MPI-2 based shared memory level
!
!     JKRRUNTYPE     - determines KR run type:
!                                             0: second-order SCF
!                                             1: second-order MCSCF
!                                             2: KR-CI
!
!     spinfr_krmc    - flag for running a spinfree MCSCF
!
!
!     LVCORR correction (re-)calculated within the KR-MC/CI framework
!     --> E_LVCORR_MC 
!
!     Parameters for MCTYPE:
!
      INTEGER JDHF, JODHF, JCAS, JMCMIN, JGAS
      PARAMETER (JDHF = 0, JODHF = 1, JCAS = 2, JMCMIN = JCAS, JGAS =   &
     &           3)

      INTEGER IPROPT, IOPT_ITRINT(2), IOPT_INTDEF, IOPT_INTFLG,         &
     &        IOPT_INTBUF, IOPT_IHRMOP(2), IOPT_IFCKOP(2),              &
     &        IOPT_ISYMOP(2), NZVAR, NZVARQ, NZCONF, NZCONFQ, NZXOPT,   &
     &        NZXOPTQ, NZXOPE, NZXOPEQ, NZXOPP, NZXOPPQ, NZHOPT,        &
     &        NZHOPTQ, NZHOPE, NZHOPEQ, MOPT_MXMICRO, MOPT_MXMACRO,     &
     &        MOPT_MAXBCK, NCIROOT, MXCIV, IRESTRK, IOPT_STATE, IANACI, &
     &        IOPT_SYMMETRY, KTRLVL, JTRLVL, LENPV, ITMAC, ITMICT,      &
     &        ITBCK, ITABS, ITMACN, MAXCIT, IOPTST, MCTYPE,             &
     &        IOPT_MK2DEL, IOPT_MK2REF, IOPT_MINMK2, IOPT_MAXMK2,       &
     &        JHESSIAN, IKRMCCNO, IDOPARIO, IMEMFAC,                    &
     &        IDENSLR_STATE, ISHMEM_TYPE, NDSH(2), NDSHT, N2DSHX,       &
     &        N2DSHXQ, NVSH(2), NVSHT, N2VSHX, N2VSHXQ, JKRRUNTYPE,     &
     &        NMAX_SYM, nopen_mc, ISVRONO(2)
      COMMON /DCIOPT/ IPROPT, IOPT_ITRINT, IOPT_INTDEF, IOPT_INTFLG,    &
     &                IOPT_INTBUF, IOPT_IHRMOP, IOPT_IFCKOP,            &
     &                IOPT_ISYMOP, NZVAR, NZVARQ, NZCONF, NZCONFQ,      &
     &                NZXOPT, NZXOPTQ, NZXOPE, NZXOPEQ, NZXOPP, NZXOPPQ,&
     &                NZHOPT, NZHOPTQ, NZHOPE, NZHOPEQ, MOPT_MXMICRO,   &
     &                MOPT_MXMACRO, MOPT_MAXBCK, NCIROOT, MXCIV, IANACI,&
     &                IRESTRK, IOPT_STATE, IOPT_SYMMETRY, KTRLVL,       &
     &                JTRLVL, LENPV, ITMAC, ITMICT, ITBCK, ITABS,       &
     &                ITMACN, MAXCIT, IOPTST, MCTYPE,                   &
     &                IOPT_MK2DEL,IOPT_MK2REF,IOPT_MINMK2,IOPT_MAXMK2,  &
     &                JHESSIAN, IKRMCCNO, IDOPARIO, IMEMFAC,            &
     &                IDENSLR_STATE, ISHMEM_TYPE, NDSH, NDSHT, N2DSHX,  &
     &                N2DSHXQ, NVSH, NVSHT, N2VSHX, N2VSHXQ, JKRRUNTYPE,&
     &                NMAX_SYM, nopen_mc, ISVRONO

      REAL*8 OPT_CNVINT(2), GNORM(5), EMCSCF, EMCOLD, DEPRED, BETA,     &
     &       BETMIN, BETMAX, RTRUST, RTTOL, RATMIN, RATGOD, EMY, EACTIV,&
     &       RATREJ, STPLEN, STPMAX, STPINC, STPRED, THQMIN, THQLIN,    &
     &       THQKVA, GAMMA, RBETA, STPC, STPO, ECORR, OPT_THRGRD,       &
     &       STPLCI, STPLEE, STPLEP, DLVLSH, THRPCI, COMPFAC, ECORE_LR, &
     &       E_LVCORR_MC,OPT_THR_CVEC, opt_threci
      COMMON /DCROPT/ OPT_CNVINT, GNORM, EMCSCF, EMCOLD, DEPRED, BETA,  &
     &                BETMIN, BETMAX, RTRUST, RTTOL, RATMIN, RATGOD,    &
     &                EMY, EACTIV, RATREJ, STPLEN, STPMAX, STPINC,      &
     &                STPRED, THQMIN, THQLIN, THQKVA, GAMMA, RBETA,     &
     &                STPC, STPO, ECORR, OPT_THRGRD, STPLCI, STPLEE,    &
     &                STPLEP, DLVLSH, THRPCI, COMPFAC, ECORE_LR,        &
     &                E_LVCORR_MC,OPT_THR_CVEC, opt_threci

      INTEGER NOPTFLAGS
      PARAMETER (NOPTFLAGS = 20)

      LOGICAL FLAG(NOPTFLAGS)
      LOGICAL OPT_SKIPEE, OPT_SKIPEP, OPT_NOPFQ, OPT_NOFQX
      LOGICAL OPT_NOCI, OPT_UCIBOS, OPT_CHCKJZ
      LOGICAL DONR, FROZCI, LDIRNR, NATONL, FOCKON, COMPROT, NATOLCR
      LOGICAL CRDFO_MAT, CWRTFO_MAT, CSHMEMO, CINT_REORD, CINT_SPLIT
      LOGICAL GENFOCK, CHCKPT_WRT, CPBLCK_FILE, TRA_NATO, CINT_LOWSRT
      LOGICAL spinfr_krmc, cana_mcscf, no1pdens_save, save_reordered_nos
      LOGICAL fcidump
!     cana_mcscf means: analysis of mcscf reference vector in a separate 
!                       ci run (after mcscf convergence)
      COMMON /DCLOPT/ OPT_SKIPEE, OPT_SKIPEP, OPT_NOPFQ, OPT_NOFQX,     &
     &                OPT_NOCI, OPT_UCIBOS,             DONR, FROZCI,   &
     &                LDIRNR, FLAG, NATONL, FOCKON, COMPROT, NATOLCR,   &
     &                CRDFO_MAT, CWRTFO_MAT, CSHMEMO, CINT_REORD,       &
     &                CINT_SPLIT, GENFOCK, CHCKPT_WRT, CPBLCK_FILE,     &
     &                TRA_NATO, OPT_CHCKJZ, CINT_LOWSRT, spinfr_krmc,   &
     &                cana_mcscf, no1pdens_save, save_reordered_nos,    &
     &                fcidump

      CHARACTER*72 OPT_INDSTR(3, 2)
      CHARACTER*12 OPT_INTTYP
      CHARACTER*8  OPT_CIPROGRAM
      CHARACTER*72 OPT_FRZSTR(2), OPT_DELSTR(2)
      COMMON /DCCOPT/ OPT_INTTYP, OPT_INDSTR, OPT_CIPROGRAM, OPT_FRZSTR,&
     &                OPT_DELSTR

      INTEGER         KZCONF, LZCONF, LZXOPE, LZXOPP
      COMMON /DCMOPT/ KZCONF, LZCONF, LZXOPE, LZXOPP

!     lzconf is integer length of WORK(KZCONF) etc.
!
!
      INTEGER NOPTTIM
      PARAMETER (NOPTTIM = 20)

!     Detailed timings are saved in CPUOPT and WLLOPT.
!
!     01: time in ROPTST     (wf. start guess)
!     02: time in RCIST      (CI start guess)
!     03: time in RGETH2     (FQ and H2AC)
!     04: time in RGETH2TX   (FQX, FQT, and H2ACX)
!     05: time in GMOLITX    (FCX, FVX, and FVT)
!     06: time in FMOLI      (one index transforms)
!     07: time in RTRACTL    (MO 4-index transformation)
!     08: time in RMAKDM     (density matrices)
!     09: time in RGRAD      (gradient)
!     10: time in RCIGRAD    (CI gradient)
!     11: time for CI sigma vectors
!     12: time for RFCKMAT   (FC and FV)
!     13: time for XRSSEP    (orbital parts of orbital and conf. sigma vec)
!     14: time for XRSTDM    (transition density matrices)
!     15: time in RTR1H1     (one index transforms)
!     16: time for RSIGOC, RSIGOO (orbital sigma vectors)
!     17: time for RSIGCO    (conf part of orbital sigma vectors)
!     18: time for XRSSVC    (sigma vectors)
!
      REAL*8 CPUOPT(NOPTTIM), WLLOPT(NOPTTIM)
      COMMON /DCTOPT/ CPUOPT, WLLOPT
