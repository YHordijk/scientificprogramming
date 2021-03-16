! dcbxqr.h - control common block for pamxqr.h and pamtpa.h
! **************************************************************
! ***********    Quadratic response  <<A;B,C>>    **************
! **************************************************************
!
! MAXQR  - max number of allowed QR functions
! MAXLQR - max number of operators specified
! MAXFQR - max number of frequencies specified
! MZYEE  - max length of electronic response vectors
! MZYEP  - max length of positronic response vectors
!
! Used to read in USER specified requests
!----------------------------------------
! LAQROP - list of A-operator pointers
! LBQROP - list of B-operator pointers
! LCQROP - list of C-operator pointers
! IPRXQR - general print level
! ITRXQR - max number of iterations
! NBQRFR - number of B frequencies
! BQRFR  - B frequencies
! NCQRFR - number of C frequencies
! CQRFR  - C frequencies
! XQR_SKIPEE - do not include e-e rotations
! XQR_SKIPEP - do not include e-p rotations
! XQR_ALLCMB - calculate all response functions disregarding
!              permutational symmetries (debug option)
! XQR_LEVICI - Include only functions with 3 different cartesian
!              indices (for .VERDET)
! XQR_INDSTR - list of orbitals to include in rotation vectors
!
! Used to store the necessary linear response equations
!------------------------------------------------------
! NQROP  - number of unique operator pointers which are to
!          to be solved linear response for
! LQROP  - list of operator pointers, and number of
!          frequencies for that operator
! QRFREQ - list of frequencies for operators
!
! Used to define nontrivial QR response functions
!------------------------------------------------
! NQRHYP  - number of nontrivial QR functions to be computed
! LQRHYP  - list of A/B/C operator pointers
! QRFRHYP - list of A/B/C frequencies
!
! Used in TPA calculations (EXCIT is also used by EXCPRP)
!--------------------------------------------------------
! IEXCOFF - offset when creating direct access labels for excitations
! TPACNV  - number of states per irrep
! TPASIM  -
! TPASTV  -
! EXCIT   - list of excitation frequencies in each irrep
!
! Used in EXCPRP calculations
! ---------------------------
! LSTATE - list of left excited states
! RSTATE - list of right excited states
!
      INTEGER MAXLQR,MAXFQR,MAXQR,IEXCOFF,LUQRINFO,LUQREE,LUQREP
!
      PARAMETER ( MAXLQR=10,MAXFQR=15,MAXQR=70,IEXCOFF=90 )
      PARAMETER ( LUQRINFO=40,LUQREE=41,LUQREP=42 )
!
      LOGICAL XQRNRM,XQRDIH,XQR_SKIPEE,XQR_SKIPEP,XQR_ALLCMB,           &
     &     XQR_LEVICI
      INTEGER IPRXQR,ITRXQR,NBQRFR,NCQRFR,NQROP,MAXQRM,                 &
     &     LAQROP,NAQROP,LBQROP,NBQROP,                                 &
     &     LCQROP,NCQROP,LQROP,LQRHYP,NQRHYP,INTXQR,ITRIQR,             &
     &     MZYEE,MZYEP,TPACNV,TPASIM,TPASTV,LSTATE,RSTATE,              &
     &     ndamp_qr
      REAL*8 THCQR,BQRFR,CQRFR,QRFREQ,QRFRHYP,CNVXQR,RESXQR,EXCIT,      &
     &     dampfreq_qr
      CHARACTER XQR_INDSTR*72
!
      COMMON /XCBLQR/ XQRNRM,XQRDIH,XQR_SKIPEE,XQR_SKIPEP,XQR_ALLCMB,   &
     &     XQR_LEVICI
      COMMON /XCBRQR/ THCQR,BQRFR(MAXFQR),CQRFR(MAXFQR),                &
     &     QRFREQ(3*MAXFQR,3*MAXLQR),QRFRHYP(MAXQR,3),                  &
     &     CNVXQR(2),RESXQR,EXCIT(MAXFQR,8),                            &
     &     dampfreq_qr
      COMMON /XCBIQR/ IPRXQR,ITRXQR,NBQRFR,NCQRFR,NQROP,MAXQRM,         &
     &     LAQROP(MAXLQR),NAQROP,LBQROP(MAXLQR),NBQROP,                 &
     &     LCQROP(MAXLQR),NCQROP,LQROP(3*MAXLQR,2),                     &
     &     LQRHYP(MAXQR,3),NQRHYP,INTXQR,ITRIQR(2),MZYEE,MZYEP,         &
     &     TPACNV(8),TPASIM(8),TPASTV(8),LSTATE(8),RSTATE(8),           &
     &     ndamp_qr
      COMMON /XCBCQR/XQR_INDSTR(3,2)
