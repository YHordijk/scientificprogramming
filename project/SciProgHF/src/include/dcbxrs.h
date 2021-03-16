!
! FILE    : dcbxrs.h
!
!*** Reduced response equations
!
!    JSYMOP - point group symmetry of operator
!    JTIMOP - time reversal symmetry of operator: 1(+),2(-)
!    JOPSY  - pointer to fermion ircop which spans JSYMOP
!
!    IPRXRS - print level
!    MAXITR  - max number of iterations 
!    NFREQ  - number of frequencies
!    NEXSIM - number of simultaneous solutions
!    NEXSTV - number of start vectors
!
!    FKRMC  - is this a KR-MCSCF optimization
!    ENRGY - energy
!    NSTAT - integer factor =1 in the static case and =2 else
!    NDAMP - integer factor =1 without damping, =2 with damping
!
!    THCXRS - threshold for convergence of reduced response equations
!     
!    Configuration and orbital rotation parameters
!    are held in separate vectors
!
!    Trial vectors have the general structure:
!
!     (b_1)^T = [ Z  Y^*]     (b_2)^T  =  [Y   Z^*]
!
!    
      CHARACTER*12 XRS_INTTYP
      CHARACTER*2 BVTYP(6)
      CHARACTER*16 RSPLAB
      CHARACTER*8 XRS_CIPROGRAM
      CHARACTER*72 INDSTR
      COMMON /XRSCCB/ XRS_INTTYP, BVTYP, RSPLAB, XRS_CIPROGRAM,         &
     &                INDSTR(3,2)
!     Parameters for MCTYPE:
!
      INTEGER JDHF, JODHF, JCAS, JMCMIN, JGAS
      PARAMETER (JDHF = 0, JODHF = 1, JCAS = 2, JMCMIN = JCAS, JGAS =   &
     &           3)

      LOGICAL RSREST, STATIC, LINEQ, LSVCFG(2), TKNORM, DIAGHE, STERNH, &
     &        FKRMC, XRS_NOPFQ, XRS_NOFQX, UNCOUP, TRIPLET, IMFREQ,     &
     &        E2CHEK,STERNC,SKIPEE,SKIPEP
      INTEGER IPRXRS, JSYMOP, JTIMOP, JOPSY, MAXITR, ITMIC, NFREQ,      &
     &        NEXSIM, NEXSTV, NEXCNV, LOFFTY, IPX, IMX, MAXSIM, NZVAR,  &
     &        NZVARQ, NZCONF, NZCONFQ, NTYP, NZXOPT, NZXOPTQ, NZXOPE,   &
     &        NZXOPEQ, NZXOPP, NZXOPPQ, NZHOPE, NVPAR, IPEP, IPEM, IPPP,&
     &        IPPM, IPCP, IPCM, NREDM, N2REDM, NZRED, NCRED, NERED,     &
     &        NPRED, KCONV, NCSIM, NESIM, NPSIM, NOSIM, NTSIM, IZRED,   &
     &        ICRED, IERED, IPRED, INTXRS, INTFLG, INTDEF, INTBUF,      &
     &        ITRINT(2), NFOCKMAT, IXRSTB(-1:2), MCTYPE, JTRLVL,NEVEC,  &
     &        NDAMP, NSTAT
      COMMON /XRSICB/ IPRXRS, JSYMOP, JTIMOP, JOPSY, MAXITR, ITMIC,     &
     &                NFREQ, NEXSIM, NEXSTV, NEXCNV, LOFFTY, IPX, IMX,  &
     &                MAXSIM, NZVAR, NZVARQ, NZCONF, NZCONFQ, NTYP,     &
     &                NZXOPT, NZXOPTQ, NZXOPE, NZXOPEQ, NZXOPP, NZXOPPQ,&
     &                NZHOPE, NVPAR, IPEP, IPEM, IPPP, IPPM, IPCP, IPCM,&
     &                NREDM, N2REDM, NZRED, NCRED, NERED, NPRED, KCONV, &
     &                NCSIM, NESIM, NPSIM, NOSIM, NTSIM, IZRED, ICRED,  &
     &                IERED, IPRED, INTXRS, INTFLG, INTDEF, INTBUF,     &
     &                ITRINT, NFOCKMAT, IXRSTB, MCTYPE, JTRLVL,NEVEC,   &
     &                NDAMP, NSTAT

      REAL*8 THCXRS, CNVINT(2), RESFAC, GPTNRM, ENRGY, DAMPFREQ
      COMMON /XRSRCB/ THCXRS, CNVINT, RESFAC, GPTNRM, ENRGY, DAMPFREQ

      COMMON /WRKFLG/ STATIC, RSREST, LINEQ, LSVCFG, TKNORM, DIAGHE,    &
     &                STERNH, FKRMC, XRS_NOPFQ, XRS_NOFQX, UNCOUP,      &
     &                TRIPLET, IMFREQ, E2CHEK,STERNC,SKIPEE,SKIPEP

      INTEGER         KBTOT, KBCI, KBOE, KBOP, KZCONF
      COMMON /MEMXRS/ KBTOT, KBCI, KBOE, KBOP, KZCONF

      INTEGER LUBCI, LUBOE, LUBOP, LUSCI, LUSOE, LUSOP, LURSP, LUXVC,   &
     &        LURST, LUCYCL
      COMMON /XRSIOU/ LUBCI, LUBOE, LUBOP, LUSCI, LUSOE, LUSOP, LURSP,  &
     &                LUXVC, LURST, LUCYCL
