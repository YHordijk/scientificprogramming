!
! FILE    : dcbgen.h (Dirac Common Block GENeral information)
!
!*** General information in DIRAC ***
!
!  DOPUT - write MO-coefficients to formatted file
!
!  INTGEN - Give top level integral flag
!  DOLVC  - Treat SS 2-e interaction by point charge model
!  DOTSC  - Treat SS 2-e interaction by TS point charge model
!  DIRSET - Conventional OR direct calculation of 2-e integrals
!           has been specified manually (used for SMLV1C)
!  BNCRON - Improved screening of nucleus for bare nucleus starting guess
!  force_bnc - use input value of BNCRON w/o checking molecular charge (PSB+MI)
!  BNCRON - Use sum of atomic potentials for the starting guess
!
!  TWOCOMP, TWOCOMPBSS - related to the pure two-component mode as follows:
!    TWOCOMP = T means that the program is currently running in 2c mode (either nonrelativistic or 2c-relativistic)
!    TWOCOMPBSS = T means that this is a 2-component relativistic calculation in 2c mode
!    X2Cmod_x2c = T means that this is a 2-component relativistic calculation in 2c mode using the X2C module
!      (now stored in the module x2cmod_cfg)
!
!  SKIP2E  - Skip all two electron steps.
!
!  DOJACO  - do the Jacobi diagonalization of real matrixes
!  DOQJACO - Jacobi diagonalization of quaternion Jacobi matrixes
!
!  These flags should be passed on to every module
!
! NOSFCRASH - No Seg.Fault crash (by call quit). This is important for preserving gmon.out
!
! DIRAC_WF(1:N_WF) the wave function jobs to be done in this dirac run, set in dirac/dirrdn.F
!
      INTEGER LBASDIR
      PARAMETER (LBASDIR = 600)

      CHARACTER*50 TITLE
      CHARACTER*10 DIRAC_WF(100)
      COMMON /CBCPAM/ TITLE, DIRAC_WF

      LOGICAL FAMILY, INPTES, VACUUM, NUCVAC, DOPSI, DOTRA, NOTRA,      &
     &        DOANA, DOPRP,DOHRM, DOJACO,DOPUT,DORKBIMP,                &
     &        OPTIMI, HRINPC, RDINPC, PARCAL, DIRCAL, DOFCK3, SEGBAS,   &
     &        NEWBAS, NEWPRP, RELCAL, NEWSYM, RDMLIN, USRIPR, DOLVC,    &
     &        LVNEW, DOTSC, DIRSET, BNCRON, DOACUT, MP2ORG,             &
     &        TWOCOMP, TWOCOMPBSS, NOSET,                               &
     &        DOACIN, CBKRINI, SKIP2E, LOWJACO, NOSFCRASH, DOQJACO,     &
     &        QM3, QMMM, BARNUC, DOHUCKEL, force_bnc, BNSPON
      LOGICAL OPTWLK, OPTNEW,  NMWALK
!     hjaaj: preparing for implementation of QMMM /Jan 2010

      COMMON /DCBGEN_L_PAM/                                             &
     &                FAMILY, INPTES, VACUUM, DOPSI, NUCVAC, DOTRA,     &
     &                NOTRA, DOPRP, DOANA,    DOHRM, DOJACO,            &
     &                DOFCK3,DORKBIMP,DOPUT, OPTIMI, HRINPC, RDINPC,    &
     &                PARCAL, DIRCAL, SEGBAS, NEWBAS, NEWPRP, RELCAL,   &
     &                NEWSYM, RDMLIN, USRIPR, DOLVC, LVNEW, DOTSC,      &
     &                BNCRON, DOACUT, MP2ORG, TWOCOMP,                  &
     &                TWOCOMPBSS, NOSET,                                &
     &                DOACIN, CBKRINI, SKIP2E, LOWJACO, NOSFCRASH,      &
     &                DOQJACO, QM3, QMMM, BARNUC, DOHUCKEL,             &
     &                OPTWLK, OPTNEW, NMWALK, force_bnc, BNSPON

      INTEGER IPRGEN, ISPHTR, IDFLAG, LENBAS, IPRUSR, KCHARG, INTGEN,   &
     &        ILLDIR, ISLDIR, ISSDIR, IGTDIR, LENDND, N_WF, INTGEN_SAVE
      COMMON /DCBGEN_I_PAM/                                             &
     &                IPRGEN, ISPHTR, IDFLAG, LENBAS, IPRUSR, KCHARG,   &
     &                INTGEN, ILLDIR, ISLDIR, ISSDIR, IGTDIR, DIRSET,   &
     &                LENDND, N_WF, INTGEN_SAVE

!     CVAL   : speed of light used in Dirac
!     STOL(2): thresholds for linear dependence, large and small comp. bf.
!     PANAS  : use Panas correlation correction to 2-electron integrals
      REAL*8 CVAL, STOL(2), PANAS
      COMMON /DCBGEN_R_PAM/ CVAL, STOL, PANAS

      INTEGER LUCOEF, LUOVLP, LUTMAT, LUPMAT, LU1INT, LUKRMC, LUKRM1,   &
     &        LUKRM2, LUKRM3, LUKRM4, LUKRM5, LUITFO, LUBSS,  LUX2C
      COMMON /PAMIOU/ LUCOEF, LUOVLP, LUTMAT, LUPMAT, LU1INT, LUKRMC,   &
     &                LUKRM1, LUKRM2, LUKRM3, LUKRM4, LUKRM5, LUITFO,   &
     &                LUBSS,  LUX2C

      CHARACTER*(LBASDIR) BASDIR, CDIRNOD
      COMMON /DCBGEN_C_PAM/ BASDIR, CDIRNOD
! --- end of dcbgen.h ---
