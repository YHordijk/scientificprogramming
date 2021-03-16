!
! *** PP (Polarization propagator) excitation properties  ***
!
! MAXLPP - max number of A operators specified
!
! NPPAP(NBSYM) - number of A-operators classified 
!      according to boson symmetries
! LPPAPU - pointer to list of A - operators (unsorted)
! LPPAPS - pointer to list of A - operators sorted on symmetries
! JPPAP  - symmetry offsets in LPPAPS
!
! IPRXPP - general print level
! ITRXPP - max number of iterations 
! THCPP  - convergence threshold
! KEXCNV(8) - number of excitations to converge in each boson sym.
! KEXSIM(8) - number of simult. trial vectors in each boson sym.
! KEXSTV(8) - number of start   trial vectors in each boson sym.
! MAXEXC = max(KEXCNV(i), i = 1,nbsym)
! NROTAV = maxL for Lebedev quadrature for rotational averaging
! IEPOFF(IORDER) - offset in list of properties to electric multipoles of IORDER
! IMPOFF(IORDER) - offset in list of properties to magnetic multipoles of IORDER
! IBDOFF(IORDER) - offset in list of properties to truncated light-matter interaction of IORDER (in wave vector)
! GNOISE - turn on noise in coefficients upon gradient evaluation
! GNOISE_FAC - noise level  
! GCLEAN - clean coefficients upon gradient evaluation
! GCLEAN_FAC - relative threshold for setting coefficients to zero
!
!    dcbxpp.h -- Feb. 18, 1999, hjaaj
!
      PARAMETER (MAXLPP = 5316, MAX_KVAL = 20)
      LOGICAL XPP_SKIPEE,XPP_SKIPEP,XPP_LSFG,                           &
     &        XPPNRM,ORBXPP,XPPANA,XPP_TRIPLET,                         &
     &        XPP_E2CHEK,DOBED,ORIENTED,BEDCHK,GNOISE,GCLEAN,           &
     &        DOVELR,DOLENR
      COMMON /XCBLPP/XPP_SKIPEE,XPP_SKIPEP,XPP_LSFG(2),                 &
     &        XPPNRM,ORBXPP,XPPANA,XPP_TRIPLET,XPP_E2CHEK,DOBED,        &
     &        ORIENTED,BEDCHK,GNOISE,GCLEAN,DOVELR,DOLENR
      COMMON /XCBRPP/THCPP,CNVXPP(2),RESXPP,XPPERG,                     &
     &        UWAVE(3),UPOL(3),BEDSTP,GNOISE_FAC,GCLEAN_FAC
      COMMON /XCBIPP/IPRXPP,ITRXPP,MAXEXC,                              &
     &     KEXCNV(8), KEXSIM(8), KEXSTV(8),                             &
     &     LPPAPU(MAXLPP),LPPAPS(MAXLPP),NPPAPT,NPPAP(8),JPPAP(8),      &
     &     INTXPP,ITRIPP(2),NCLS_XPP(2,3,2),KVAL_OSC,                   &
     &     IEPOFF(MAX_KVAL+1),IMPOFF(MAX_KVAL), IBDOFF(0:MAX_KVAL),     &
     &     NROTAV,NBED
      CHARACTER XPP_INDSTR*72, XPPFIL*6
      COMMON/DCCXPP/XPP_INDSTR(3,2),XPPFIL
