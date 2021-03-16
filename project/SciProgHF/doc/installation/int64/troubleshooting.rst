:orphan:
 

DFT runs stop with "programming error in distribution of points in parallel DFT" with 64bit integer Intel MPI
-------------------------------------------------------------------------------------------------------------

This is not a programming error in DIRAC but a bug in 64bit integer Intel MPI
mpi_reduce.  The Intel developers are informed and working on a patch.  In the
meantime, for DFT calculations using Intel MPI you should use 32bit integers.
