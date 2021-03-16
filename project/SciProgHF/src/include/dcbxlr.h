!
! File: dcbxlr.h
! *** Info for linear response  <<A,B>>, label: "XLR"  ***
!
! ----- PARAMETERS:
!
! MAXLLR - max number of operators allowed in XLR
! MAXFLR - max number of frequencies allowed in XLR
! Note: if you change any of these two parameters, then you
!       must rebuild with 'make'.
!
      INTEGER MAXLLR, MAXFLR
      PARAMETER (MAXLLR = 60, MAXFLR = 60)
!
! ----- VARIABLES:
!
! NLRAP/NLRBP(NBSYM,ITIM) - number of A/B-operators classified 
!      according to boson symmetries and time reversal symmetry
! LLRAPU/LLRBPU - pointer to list of A/B - operators (unsorted)
! LLRAPS/LLRBPS - pointer to list of A/B - operators sorted on symmetries
! JLRAP/JLRBP   - symmetry offsets in LLRAPS/LLRBPS
!
! IPRXLR - general print level
! ITRXLR - max number of iterations 
! MAXRM  - max dimension of reduced matrix
! THCLR  - convergence threshold
! NBFREQ - number of frequencies 
! BFREQ  - frequencies
! INDXLR is pointer to current B-operator in LLRAPU
! XLR_LSFG(1)=T   - call GMOLI in sigma-vector generation
! XLR_LSFG(2)=T   - call FMOLI in sigma-vector generation
!
!
      LOGICAL XLR_SKIPEE, XLR_SKIPEP, XLR_LSFG(2), XLRNRM, XLRDIH,      &
     &        ALLCMB, TRIAB, COMPRS, XSTERN, XLRANA, NOSPIA, NOSPIB,    &
     &        XLR_UNCOUP, XLR_TRIPLET, XLR_IMFREQ, XLR_E2CHEK,XSTERNC
      COMMON /XCBLLR/ XLR_SKIPEE, XLR_SKIPEP, XLR_LSFG, XLRNRM, XLRDIH, &
     &                ALLCMB, TRIAB, COMPRS, XSTERN, XLRANA,            &
     &                NOSPIA, NOSPIB, XLR_UNCOUP, XLR_TRIPLET,          &
     &                XLR_IMFREQ,XLR_E2CHEK,XSTERNC

      REAL*8 THCLR, BFREQ(MAXFLR), THRCOM, CNVXLR(2), RESXLR, DMPFRLR,  &
     &       XLR_ANATHR
      COMMON /XCBRLR/ THCLR, BFREQ, THRCOM, CNVXLR, RESXLR, DMPFRLR,    &
     &       XLR_ANATHR

      INTEGER IPRXLR, ITRXLR, NBFREQ, INDBP, MAXRM, MXLOAD,             &
     &        LLRAPU(MAXLLR), LLAPSU(MAXLLR), NLRAPT, NLRAP(8),         &
     &        JLRAP(8), LLRBPU(MAXLLR), LLBPSU(MAXLLR), NLRBPT,         &
     &        NLRBP(8), JLRBP(8), INTXLR, ITRILR(2), NCLS_XLR(2, 3, 2)
      COMMON /XCBILR/ IPRXLR, ITRXLR, NBFREQ, INDBP, MAXRM, MXLOAD,     &
     &                LLRAPU, LLAPSU, NLRAPT, NLRAP, JLRAP, LLRBPU,     &
     &                LLBPSU, NLRBPT, NLRBP, JLRBP, INTXLR, ITRILR,     &
     &                NCLS_XLR

      CHARACTER*72 XLR_INDSTR(3, 2)
      CHARACTER*6 XLR_XVCFIL
      COMMON /DCCXLR/ XLR_INDSTR, XLR_XVCFIL
