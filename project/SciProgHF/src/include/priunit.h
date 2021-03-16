! priunit.h : LUCMD is for reading DIRAC.INP, LUPRI is for output
! - (LUCMD should not be = 5, this has been seen to cause problems with mpi for gfortran 8.2
!    when all slaves try to close unit 5 stdin to open unit 5 DIRAC.INP in relccsd Feb. 2019)
! - "LUPRI" should not be parameter, it may be temporarily changed to redirect output to a file
!   instead of stdout, especially from MPI slaves.
      integer, parameter :: LUCMD   = 4  ! for reading from DIRAC.INP
      integer, parameter :: LUSTDIN = 5  ! for reading standard input
      integer            :: LUPRI   = 6  ! not defined as parameter, due to redefinition in LUCIAREL_NODE
