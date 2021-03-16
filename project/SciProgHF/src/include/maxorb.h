!
! File: maxorb.h
!
! hjaaj July 2000:
! MXSHEL = max number of symmetry independent contracted shells
! MXPRIM = max number of symmetry independent primitive shells
!  (e.g. a p shell contains 3 components)
! MXCORB = max total number of contracted functions = max NBASIS
!
      INTEGER MXSHEL, MXPRIM, MXCORB
      PARAMETER (MXSHEL = 15000, MXPRIM = 15000, MXCORB = 40000)
! hjaaj Jan 2019: fitting is not implemented in Dirac, so
!  to save static memory I introduced the following three
!  variables for the fitting common blocks
      INTEGER MXSHEL_FIT, MXPRIM_FIT, MXCORB_FIT
      PARAMETER (MXSHEL_FIT = 1, MXPRIM_FIT = 1, MXCORB_FIT = 1)
