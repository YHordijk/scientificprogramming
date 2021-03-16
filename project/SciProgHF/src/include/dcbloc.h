! ** COMMON block for Molecular orbital localization
!     PRJLOC - localization based on:
!              True  - projection analysis
!              False - Mulliken analysis
!     IPRLOC - print level
!     ITRLOC - maximum number of iterations
!     SELMOS - string which contains information about MOs that will be localized
!
!     ---- Convergence keywords ----
!
!     HESLOC = DIAG - diagonal approximation of the hessian is used
!            = FULL - calculation of the full hessian is used
!            = COMB - there are two stages:
!                     1: in the first stage only diagonal approximation of the hessian is used
!                     2: in the second stage the calculation of the full hessian is used
!     LGFULL - it is used when HESLOC = FULL or in the second stage of the HESLOC = COMB option
!              True  - converge when the difference between the functional value
!                      of two subsequent iterations is lower than threshold (THFULL)
!              False - the gradient criterion or criterion of the number
!                      of negative eigenvalues equals zero will be used
!     LGGRAD - it is used when HESLOC = FULL or HESLOC = DIAG or
!              in the second stage of the HESLOC = COMB option
!              True  - converge when the difference between the functional value
!                      of two subsequent iterations is lower than threshold (THGRAD)
!              False - the functional value criterion or criterion of the number
!                      of negative eigenvalues equals zero will be used
!     LGDIAG - if HESLOC = DIAG it have the the same purpose as LGFULL;
!              for threshold is used THDIAG variable
!            - if HESLOC = COMB it is used to switch from the first to the second
!              stage of the convergence; it will switch
!              True  - when the difference between the functional value 
!                      of two subsequent iterations is lower than threshold (THDIAG)
!              False - if the number of negative eigenvalues equals zero
!     LGCHCK - works with HESLOC = COMB option
!              True  - only one calculation of the full hessian will be performed
!
!#include <dcblab.h>
      LOGICAL PRJLOC,LGFULL,LGDIAG,LGGRAD,LGCHCK
      INTEGER IPRLOC,ITRLOC
      REAL*8  THFULL,THDIAG,THGRAD
      CHARACTER*4 HESLOC
      CHARACTER*72 SELMOS
      COMMON/CBILOC/IPRLOC,ITRLOC
      COMMON/CBRLOC/THFULL,THDIAG,THGRAD
      COMMON/CBCLOC/HESLOC,SELMOS
      COMMON/CBLLOC/PRJLOC,LGFULL,LGDIAG,LGGRAD,LGCHCK
