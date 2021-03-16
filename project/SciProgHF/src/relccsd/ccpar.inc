!
!
! MPI information for RELCCSD
!
!     General MPI information
!
      INTEGER MASTER,MYPROC,NMPROC
      COMMON/CCPAR/ MASTER,MYPROC,NMPROC
!
!     MASTER: contains node number of master (0 in general)
!     MYPROC:  -"-     own process number. According to MPI
!             counting starts at zero.....# proc -1 !
!     NMPROC:  -"-     total number of processes (including master !)
!

! Information for integral distribution
!
!     Distribution of two-electron integrals
!     First index : start, stop of block on node
!     Second index : integral type (VOVV,VVVV)
!     Third index : irrep

      INTEGER*8 NV6PAR,NV5PAR,IDIST(6,10,32)
      COMMON/CCPARV/NV6PAR,NV5PAR,IDIST

! DIRAC-specific preprocessor statement to solve 32/64 bit compatibility problems 

#if defined (INT_STAR8)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

