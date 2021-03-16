#define REALK 8
#define _REALK _8
#define MPI_REALK MPI_DOUBLE_PRECISION

!   hjaaj Aug 2017: NOTE: REALK must be 8 in this version.
!   I had to replace REALK with 8 in predefined constants the fortran files,
!   because pgfortran (and standard cpp) do not replace
!   REALK or _REALK in e.g. 0.0_REALK
!   (gfortran, ifort, ... do change it to 0.0_8)

#define STDOUT 6

#define MAX_LEN_STR 80

#if defined(VAR_MPI)
#define MANAGER 0
#define REQUEST_WORK 0
#define NO_MORE_WORK 0
#endif
