!
! This program should not compile if DIRAC integer types don't match linked MPI-library integer types.
! For example, if DIRAC is compiled with 8-bytes integer and corresponding MPI library is integer*4.
!
! Program contains few representative MPI calls which are supposed to be 'sensitive enough' 
! to give compilation error when hitting incompatible integer types.
!
! NOTE: The program may not flush compilation error for ALL MPI libraries, therefore internal MPI check
! of integer*4/8 compatibility must be done inside dirac.x.
!
! Written by Rado Bast
!      extended by Miro Ilias,June 2020,GSI.de
!
Program MPI_int_compat
 implicit none
contains
   subroutine check_mpi_comm_rank()
#ifdef USE_MPI_MOD_F90
      use mpi
#else
#include "mpif.h"
#endif
      integer :: get_my_rank
      integer :: irank
      integer :: ierr
      call mpi_comm_rank(mpi_comm_world, irank, ierr)
      get_my_rank = irank
   end subroutine
!
   subroutine check_mpi_probe()
#ifdef USE_MPI_MOD_F90
      use mpi
#else
#include "mpif.h"
#endif
      integer :: rank
      integer :: message_id
      integer :: communicator
      integer :: status_container(MPI_STATUS_SIZE)
      integer :: ierr
      call mpi_probe(rank,message_id,communicator,status_container,ierr)
   end subroutine
!
   subroutine check_mpi_iprobe()
#ifdef USE_MPI_MOD_F90
      use mpi
#else
#include "mpif.h"
#endif
      logical :: x
      integer :: rank
      integer :: message_id
      integer :: communicator
      integer :: status_container(MPI_STATUS_SIZE)
      integer :: ierr
      call mpi_iprobe(rank,message_id,communicator,x,status_container,ierr)
   end subroutine
!
   subroutine check_mpi_test()
#ifdef USE_MPI_MOD_F90
      use mpi
#else
#include "mpif.h"
#endif
      logical :: x
      integer :: request
      integer :: status_container(MPI_STATUS_SIZE)
      integer :: ierr
      call mpi_test(request,x,status_container,ierr)
   end subroutine
!
   subroutine check_mpi_wait()
#ifdef USE_MPI_MOD_F90
      use mpi
#else
#include "mpif.h"
#endif
      integer :: request
      integer :: status_container(MPI_STATUS_SIZE)
      integer :: ierr
      call mpi_wait(request,status_container,ierr)
   end subroutine
!
   subroutine check_mpi_get_count()
#ifdef USE_MPI_MOD_F90
      use mpi
#else
#include "mpif.h"
#endif
      integer :: status_array(MPI_STATUS_SIZE)
      integer :: datatype
      integer :: elements
      integer :: ierr
      call mpi_get_count(status_array,datatype,elements,ierr)
   end subroutine
end program
