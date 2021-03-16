!      Copyright (c) 2018 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!
!     MPI compatibility test
!
#ifdef VAR_MPI
      program selftest_mpi
!     Call MPI routines with some small tests to see if things look ok
!     implicit real*8 (A-H,O-Z)
      implicit none
#include "mpif.h"
      integer            :: mpierr, ireq
      integer            :: istat(mpi_status_size)
      integer            :: myrank, whatever
      integer            :: my_MPI_INTEGER
      integer, parameter :: master_id = 0

#ifdef INT_STAR8
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

      call mpi_init(mpierr)
      call mpi_comm_rank(mpi_comm_world, myrank, mpierr)
!
!     send dummy message (from MASTER to MASTER...)
!
      call mpi_isend(myrank,1,my_MPI_INTEGER,0,25,               &
                     mpi_comm_world,ireq,mpierr)
      call mpi_probe(mpi_any_source,mpi_any_tag,mpi_comm_world,  &
                     istat,mpierr)
      if( istat(mpi_source) .ne. master_id )then
         print *,' *** ERROR: MPI self test failed! *** '
         print *,' sender of message should have ID 0'
         print *,' but has ID ',istat(mpi_source)
         print *,' '
         print *,' Make sure that you are using the right mpirun for'
         print *,' your particular MPI library version.'
         print *,' '
#ifdef INT_STAR8
         print *,' DIRAC is configured to use integer*8 '        &
              //'(64 bit integers).'
         print *,' You need to link a MPI library that is compiled'
         print *,' using integer*8, make sure that this is the case.'
#else
         print *,' DIRAC is configured to use integer*4 '        &
              //' (32 bit integers).'
         print *,' You need to link a MPI library that is compiled'
         print *,' using integer*4, make sure that this is the case.'
#endif
         print *,' '
         print *,' Cannot continue without a working MPI library.'
         print *,' Quitting.'
         stop ' *** ERROR: MPI self test failed ***'
      end if
!
!     ... it looks ok, receive dummy message
!
      call mpi_irecv(whatever,1,my_MPI_INTEGER,0,mpi_any_tag,    &
                     mpi_comm_world,ireq,mpierr)
      call mpi_wait(ireq,istat,mpierr)

      call mpi_finalize(mpierr)

      end 
#endif
