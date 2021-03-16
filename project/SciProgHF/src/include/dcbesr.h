!*****************************************************************************
!  dcbesr.h:  ESR properties
!
!  IPRESR         : print level
!  THRPCI_ESR     : limit for printing CI coefficients in output
!
!  NESRP          : number of properties to evaluate in CI ESR
!  LESPR(1:NESRP) : pointers to property operators for CI ESR
!  IPSIGMA(i)     : pointer to the spin operator Sigma_i
!
!  THR_CIESR      : convergence threshold for CI
!  N_CIESR        : number of CI roots to converge ( .ge. MULTIPESR)
!  MULTIPESR      : "multiplicity", number of CI roots to use
!                   in quasi-degenerate first order perturbation theory
!
!  ESR CI space definition:
!
!
!    NG1 = NGAS_CIESR(1)         - number of "RAS1-GAS" spaces
!    LVL_CIESR(1:NG1,1)          - max number of holes in each "RAS1-GAS" space
!    NGSH_CIESR(1:NFSYM,1:NG1,1) - number of orbitals in each "RAS1-GAS" space
!
!    NG2 = NGAS_CIESR(2)         - number of "RAS2-GAS" spaces, always = 1 for now
!    "RAS2-GAS" space is the active space from preceeding open-shell SCF or MCSCF
!    LVL_CIESR(1:NG2,2)          - negative = no occupation restrictions
!    NGSH_CIESR(1:NFSYM,1:NG2,2) - number of orbitals in each "RAS2-GAS" space
!
!    NG3 = NGAS_CIESR(3)         - number of "RAS3-GAS" spaces
!    LVL_CIESR(1:NG3,3)          - max number of electrons in each "RAS3-GAS" space
!    NGSH_CIESR(1:NFSYM,1:NG3,3) - number of orbitals in each "RAS3-GAS" space
!
!
!  USE_KRAMERS_CONJ : use Kramers conjugation to generate the partner CI vector
!                     (solve 1 root for doublet instead of solving for 2 roots)
!
!  MAX_CIESR_IT   : max CI iterations
!  THR_CVEC_ESR   : threshold for removing elements in CI trial vectors
!
      INTEGER maxesr, npsigma
      PARAMETER (maxgas_esr = 3, maxesr = 90, npsigma = 3)

      REAL*8  THR_CIESR, THR_CVEC_ESR, THRPCI_ESR
      INTEGER IPRESR,NESRP,LESRP(maxesr),IPSIGMA(npsigma)
      INTEGER N_CIESR, MULTIPESR, MAX_CIESR_IT
      INTEGER NGAS_CIESR(3), LVL_CIESR(maxgas_esr,3)
      INTEGER NGSH_CIESR(2,maxgas_esr,3)
      LOGICAL USE_KRAMERS_CONJ

      COMMON /DCBESR/ THR_CIESR,THR_CVEC_ESR,THRPCI_ESR,                &
     &        IPRESR,   NESRP, LESRP, IPSIGMA,                          &
     &        N_CIESR, MULTIPESR, MAX_CIESR_IT,                         &
     &        NGAS_CIESR, LVL_CIESR, NGSH_CIESR,                        &
     &        USE_KRAMERS_CONJ
! --- end of dcbesr.h ---
