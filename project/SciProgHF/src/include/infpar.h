#ifdef COMMENT
! -- infpar.h --
!     my_MPI_INTEGER is used both in .c and .F routines in MPI calls
!        so we can handle "-i8" compilations on 32-bit machines,
!        using INT_STAR8 /Jan-2007 hjaaj
#endif
#if defined (INT_STAR8)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

#if defined(__CVERSION__)
extern struct common_diracinfpar {
    double totval;
#if defined (INT_STAR8)
    long iprpar, ntask, mytid, numnod, mparid;
    long icount, node, nodes, iparex;
    long max_omp_threads;
    long parher, slave, timing;
#else
    int  iprpar, ntask, mytid, numnod, mparid;
    int  icount, node, nodes, iparex;
    int  max_omp_threads;
    int  parher, slave, timing;
#endif 
} diracinfpar_;
#else /* __CVERSION__ */
!  -- File: infpar.h for DIRAC; specific information for parallel calculations
!     .. mpi variables:
!        NUMNOD = number of slaves = number of nodes - 1 (from MPI_Comm_rank)
!        MYTID  = mpi ID of this process (from MPI_Comm_rank)
!
!        PARHER true when MPI parallel Hermit calculation
!                    (only set true in subroutines HERNOD and NODTRA)
!
!     .. openMP variables:
!         MAX_OMP_THREADS
      REAL*8 TOTWAL
      INTEGER IPRPAR, NTASK, MYTID, NUMNOD, MPARID
      INTEGER ICOUNT, NODE,  NODES, IPAREX
      INTEGER MAX_OMP_THREADS ! for openMP
      LOGICAL PARHER, SLAVE, TIMING
      COMMON /DIRACINFPAR/ TOTWAL                                        &
     &       ,IPRPAR, NTASK, MYTID, NUMNOD, MPARID                       &
     &       ,ICOUNT, NODE,  NODES, IPAREX                               &
     &       ,MAX_OMP_THREADS                                            &
     &       ,PARHER, SLAVE, TIMING

! -- end of infpar.h --
#endif
