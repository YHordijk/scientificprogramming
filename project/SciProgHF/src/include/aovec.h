!
!     File: aovec.h
!
!     MXCONT = maximum number of contracted functions in an AO input block
!              (in "MOLECULE.INP" input file)
!     MXAOVC = used for MXCONT in a few places
!
!     IF you change MXCONT you should make sure
!     MAXPRD in cbirea.h is at least the same value and
!     then rebuild the program using the command "make".
!
      INTEGER MXCONT, MXAOVC
      PARAMETER (MXCONT = 40, MXAOVC = MXCONT)
