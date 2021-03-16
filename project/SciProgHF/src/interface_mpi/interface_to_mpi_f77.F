!      Copyright (c) 2019 by the authors of DIRAC.
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
! purpose: provide some interface routines to the MPI library via a MPI-f77 interface
!
!          note: this is only for some notorious difficult cases, e.g.,
!                - logical bcast
!                - work-array data distribution which somehow does not work
!                  properly with the "automatic interface" f90-detection
!
! written by s. knecht - Odense/Santiago de Chile 2011-2012.
!
#ifdef VAR_MPI

      subroutine interface_mpi_bcast_l0(x,ndim,root_proc,communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      logical(kind=integer_kind) :: x
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: root_proc
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: root_proc_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      root_proc_i4    = root_proc
      communicator_i4 = communicator
      call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, 
     &               communicator_i4, ierr_i4)
#else
      if (.not. MPI_INIT_called) return
      call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   ,
     &               communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_bcast_l1(x,ndim,root_proc,communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      logical(kind=integer_kind) :: x(*)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: root_proc
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: root_proc_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      root_proc_i4    = root_proc
      communicator_i4 = communicator
      call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, 
     &               communicator_i4, ierr_i4)
#else
      if (.not. MPI_INIT_called) return
      call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , 
     &               communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! special bcast routines for communicating data on the work array
!-------------------------------------------------------------------------------
    
      subroutine interface_mpi_bcast_i1_work_f77(x,ndim,root_proc,
     &                                           communicator)
    
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      integer(kind=integer_kind) :: x(*)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: root_proc
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: root_proc_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      root_proc_i4    = root_proc
      communicator_i4 = communicator
      call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, 
     &               communicator_i4, ierr_i4)
#else
      if (.not. MPI_INIT_called) return
      call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   ,
     &               communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_bcast_i2_work_f77(x,ndim,root_proc,
     &                                           communicator)
    
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      integer(kind=integer_kind) :: x(:,:)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: root_proc
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: root_proc_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      root_proc_i4    = root_proc
      communicator_i4 = communicator
      call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, 
     &               communicator_i4, ierr_i4)
#else
      if (.not. MPI_INIT_called) return
      call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   ,
     &               communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_bcast_r1_work_f77(x,ndim,root_proc,
     &                                           communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      real(8)                    :: x(*)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: root_proc
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: root_proc_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      root_proc_i4    = root_proc
      communicator_i4 = communicator
      call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4,
     &               communicator_i4, ierr_i4)
#else
      if (.not. MPI_INIT_called) return
      call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   ,
     &               communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_bcast_r2_work_f77(x,ndim,ndim2,root_proc,
     &                                           communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      real(8)                    :: x(ndim,*)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: ndim2
      integer(kind=integer_kind) :: root_proc
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: ndim2_i4
      integer(kind=4)            :: root_proc_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      ndim2_i4        = ndim2
      root_proc_i4    = root_proc
      communicator_i4 = communicator
      call mpi_bcast(x, ndim_i4*ndim2_i4, MPI_REAL8, root_proc_i4,
     &               communicator_i4, ierr_i4)
#else
      if (.not. MPI_INIT_called) return
      call mpi_bcast(x, ndim   *ndim2   , MPI_REAL8, root_proc   ,
     &               communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_reduce_i1_work_f77(x,y,ndim,operation,
     &                                            receiver,communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      integer(kind=integer_kind) :: x(*)
      integer(kind=integer_kind) :: y(*)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: operation
      integer(kind=integer_kind) :: receiver
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: operation_i4
      integer(kind=4)            :: receiver_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      operation_i4    = operation
      receiver_i4     = receiver
      communicator_i4 = communicator
      call mpi_reduce(x, y, ndim_i4, my_MPI_INTEGER, operation_i4, 
     &                receiver_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      call mpi_reduce(x, y, ndim   , my_MPI_INTEGER, operation   ,
     &                receiver   , communicator   , ierr   )
#endif
      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_reduce_r1_work_f77(x,y,ndim,operation,
     &                                            receiver,communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      real(8)                    :: x(*)
      real(8)                    :: y(*)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: operation
      integer(kind=integer_kind) :: receiver
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: operation_i4
      integer(kind=4)            :: receiver_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      operation_i4    = operation
      receiver_i4     = receiver
      communicator_i4 = communicator
      call mpi_reduce(x, y, ndim_i4, MPI_REAL8, operation_i4, 
     &                receiver_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      call mpi_reduce(x, y, ndim   , MPI_REAL8, operation   ,
     &                receiver   , communicator   , ierr   )
#endif
      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_allreduce_i1_work_f77(x,y,ndim,operation,
     &                                               communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      integer(kind=integer_kind) :: x(*)
      integer(kind=integer_kind) :: y(*)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: operation
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: operation_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      operation_i4    = operation
      communicator_i4 = communicator
      call mpi_allreduce(x, y, ndim_i4, my_MPI_INTEGER, operation_i4, 
     &                communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      call mpi_allreduce(x, y, ndim   , my_MPI_INTEGER, operation   ,
     &                communicator   , ierr   )
#endif
      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_allreduce_r1_work_f77(x,y,ndim,operation,
     &                                               communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      real(8)                    :: x(*)
      real(8)                    :: y(*)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: operation
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: operation_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      operation_i4    = operation
      communicator_i4 = communicator
      call mpi_allreduce(x, y, ndim_i4, MPI_REAL8, operation_i4, 
     &                   communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
      call mpi_allreduce(x, y, ndim   , MPI_REAL8, operation   ,
     &                   communicator   , ierr   )
#endif
      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_scatter_r1_work_f77(x,ndim1,y,ndim2,
     &                                             root_proc,
     &                                             communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      real(8)                    :: x(*)
      real(8)                    :: y(*)
      integer(kind=integer_kind) :: ndim1
      integer(kind=integer_kind) :: ndim2
      integer(kind=integer_kind) :: root_proc
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim1_i4
      integer(kind=4)            :: ndim2_i4
      integer(kind=4)            :: root_proc_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim1_i4        = ndim1
      ndim2_i4        = ndim2
      root_proc_i4    = root_proc
      communicator_i4 = communicator
      call mpi_scatter(x, ndim1_i4, MPI_REAL8, y, ndim2_i4, MPI_REAL8,
     &                 root_proc_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      call mpi_scatter(x, ndim1   , MPI_REAL8, y, ndim2   , MPI_REAL8,
     &                 root_proc   , communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_send_r1_work_f77(x,ndim,receiver,
     &                                          message_id,communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      real(8)                    :: x(*)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: receiver
      integer(kind=integer_kind) :: message_id
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: receiver_i4
      integer(kind=4)            :: message_id_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      receiver_i4     = receiver 
      message_id_i4   = message_id
      communicator_i4 = communicator
      call mpi_send(x, ndim_i4, MPI_REAL8, receiver_i4, 
     &              message_id_i4, communicator_i4, ierr_i4)
#else
      if (.not. MPI_INIT_called) return
      call mpi_send(x, ndim   , MPI_REAL8, receiver   , message_id,
     &              communicator   , ierr   )
#endif
      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_recv_r1_work_f77(x,ndim,receiver,
     &                                          message_id,communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      real(8)                    :: x(*)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: receiver
      integer(kind=integer_kind) :: message_id
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: receiver_i4
      integer(kind=4)            :: message_id_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: ierr_i4
      integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      receiver_i4     = receiver 
      message_id_i4   = message_id
      communicator_i4 = communicator
      call mpi_recv(x, ndim_i4, MPI_REAL8, receiver_i4, 
     &              message_id_i4, communicator_i4, 
     &              status_array_i4, ierr_i4)
#else
      integer(kind=integer_kind) :: status_array(MPI_STATUS_SIZE)
      if (.not. MPI_INIT_called) return
      call mpi_recv(x, ndim   , MPI_REAL8, receiver   , message_id,
     &              communicator   , status_array, ierr   )
#endif


      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_issend_i2_work_f77(x,ndim,receiver,
     &                                            message_id,
     &                                            communicator,request)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      real(8)                    :: x(ndim)
      integer(kind=integer_kind) :: ndim
      integer(kind=integer_kind) :: receiver
      integer(kind=integer_kind) :: message_id
      integer(kind=integer_kind) :: communicator
      integer(kind=integer_kind) :: request
      integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
      integer(kind=4)            :: ndim_i4
      integer(kind=4)            :: receiver_i4
      integer(kind=4)            :: message_id_i4
      integer(kind=4)            :: communicator_i4
      integer(kind=4)            :: request_i4
      integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
      if (.not. MPI_INIT_called) return
      ndim_i4         = ndim
      receiver_i4     = receiver 
      message_id_i4   = message_id
      communicator_i4 = communicator
      call mpi_issend(x, ndim_i4, my_MPI_INTEGER, receiver_i4, 
     &                message_id_i4, communicator_i4, request_i4, 
     &                ierr_i4)
      request = request_i4
#else
      if (.not. MPI_INIT_called) return
      call mpi_issend(x, ndim   , my_MPI_INTEGER, receiver, message_id,
     &                communicator   , request, ierr   )
#endif
      end subroutine
!-------------------------------------------------------------------------------
#else
      subroutine dummy_mpi_f77_interface
      end subroutine
#endif
