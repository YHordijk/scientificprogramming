!
! FILE    : typesz_mpi2.h
!
!     MPI-2 type size for memory allocation
!     ======================================
!
#if defined (VAR_MPI2)
! aspg, not sure how to get to use mpi or mpif.h from here, so i took
! the following two lines from mpi.h
!     INTEGER MPI_ADDRESS_KIND
!     PARAMETER (MPI_ADDRESS_KIND=8)
!
!     INTEGER(KIND=MPI_ADDRESS_KIND) ISIZE_dp, ISIZE_my_int
!     INTEGER ISIZE_dp, ISIZE_my_int, ISIZE_intdef, ISIZE_int8,
      INTEGER(KIND=8) ISIZE_dp, ISIZE_my_int, 
     &                ISIZE_intdef, ISIZE_int8,
     &                ISIZE_log
!
      COMMON/TYPE_SZMPI2/ ISIZE_dp, ISIZE_my_int, 
     &                    ISIZE_intdef, ISIZE_int8, ISIZE_log
!
#endif
