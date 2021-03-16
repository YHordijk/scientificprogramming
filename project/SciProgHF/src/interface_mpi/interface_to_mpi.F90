#define interface_to_mpi_DEBUG -1
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
! purpose: provide common interface routines from Dirac to MPI based on the
!          default integer type in the respective MPI library
!
!          the library module is by far not complete wrt all functionalities
!          available in MPI but merely resembles those that are at present 
!          used in Dirac.
!
!          note that some of the interface routines (difficult cases)
!          are located in interface_mpi/interface_to_mpi_f77.F.
!
! written by s. knecht - Odense/Santiago de Chile 2011-2012.
!
module interface_to_mpi

#ifdef VAR_MPI

  use integer_kind_mpilib
#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif

#include "interface_to_mpi_f77_typedef.h"

  public interface_mpi_bcast
  public interface_mpi_iprobe
  public interface_mpi_probe
  public interface_mpi_send
  public interface_mpi_isend
  public interface_mpi_recv
  public interface_mpi_irecv
  public interface_mpi_reduce
  public interface_mpi_allreduce
  public interface_mpi_gather
  public interface_mpi_gatherv
  public interface_mpi_allgatherv
  public interface_mpi_allgather
  public interface_mpi_test
  public interface_mpi_testall
  public interface_mpi_wait
  public interface_mpi_comm_rank
  public interface_mpi_comm_size
  public interface_mpi_comm_split
  public interface_mpi_comm_free
  public interface_mpi_barrier
  public interface_mpi_get_processor_name
  public interface_mpi_init
  public interface_mpi_finalize
  public interface_mpi_abort
  public interface_mpi_file_open
  public interface_mpi_file_close
  public interface_mpi_file_set_view
  public interface_mpi_file_seek
  public interface_mpi_file_iwrite
  public interface_mpi_file_iwrite_at_pck
  public interface_mpi_file_write_at
  public interface_mpi_file_write_at_i
  public interface_mpi_file_write_at_r
  public interface_mpi_file_write_at_bi
  public interface_mpi_file_write_at_br
  public interface_mpi_file_read_pck
  public interface_mpi_file_read_at
  public interface_mpi_file_read_at_i
  public interface_mpi_file_read_at_r
  public interface_mpi_file_read_at_bi
  public interface_mpi_file_read_at_br
  public interface_mpi_info_create
  public interface_mpi_info_free
  public interface_mpi_info_set
  public interface_mpi_pack_i0
  public interface_mpi_pack_i2
  public interface_mpi_pack_r
  public interface_mpi_unpack_i0
  public interface_mpi_unpack_i2
  public interface_mpi_unpack_r
  public interface_mpi_get_count
  public interface_mpi_type_extent
  public interface_MPI_WTIME

  interface interface_mpi_bcast
    module procedure interface_mpi_bcast_i0
    module procedure interface_mpi_bcast_i1
    module procedure interface_mpi_bcast_i2
    module procedure interface_mpi_bcast_i3
    module procedure interface_mpi_bcast_i4
    module procedure interface_mpi_bcast_r0
    module procedure interface_mpi_bcast_r1
    module procedure interface_mpi_bcast_r2
    module procedure interface_mpi_bcast_r3
    module procedure interface_mpi_bcast_r4
    module procedure interface_mpi_bcast_cx
    module procedure interface_mpi_bcast_c1
    module procedure interface_mpi_bcast_cmplx2
    module procedure interface_mpi_bcast_cmplx3
  end interface

  interface interface_mpi_gather
    module procedure interface_mpi_gather_i
    module procedure interface_mpi_gather_r
  end interface

  interface interface_mpi_allgather
    module procedure interface_mpi_allgather_i0
    module procedure interface_mpi_allgather_i
    module procedure interface_mpi_allgather_r
  end interface

  interface interface_mpi_send
    module procedure interface_mpi_send_i0
    module procedure interface_mpi_send_i1
    module procedure interface_mpi_send_r0
    module procedure interface_mpi_send_r1
    module procedure interface_mpi_send_r2
    module procedure interface_mpi_send_cx
  end interface

  interface interface_mpi_isend
    module procedure interface_mpi_isend_i0
    module procedure interface_mpi_isend_i1
  end interface

  interface interface_mpi_recv
    module procedure interface_mpi_recv_i0
    module procedure interface_mpi_recv_i0_stat
    module procedure interface_mpi_recv_i1
    module procedure interface_mpi_recv_i1_stat
    module procedure interface_mpi_recv_r0
    module procedure interface_mpi_recv_r1
    module procedure interface_mpi_recv_r1_stat
    module procedure interface_mpi_recv_r2
    module procedure interface_mpi_recv_cx
  end interface

  interface interface_mpi_irecv
    module procedure interface_mpi_irecv_i0
    module procedure interface_mpi_irecv_i1
  end interface

  interface interface_mpi_reduce
    module procedure interface_mpi_reduce_i0
    module procedure interface_mpi_reduce_i1
    module procedure interface_mpi_reduce_i2
    module procedure interface_mpi_reduce_r0
    module procedure interface_mpi_reduce_r1
    module procedure interface_mpi_reduce_r2
  end interface

  interface interface_mpi_allreduce
    module procedure interface_mpi_allreduce_i0
    module procedure interface_mpi_allreduce_i1
#ifndef INT_STAR8
    module procedure interface_mpi_allreduce_i0_i8
#endif
    module procedure interface_mpi_allreduce_r0
    module procedure interface_mpi_allreduce_r1
    module procedure interface_mpi_allreduce_r2
    module procedure interface_mpi_allreduce_r3
  end interface
  
  interface interface_mpi_iprobe
    module procedure interface_mpi_iprobe_l0
  end interface

  interface interface_mpi_file_open
    module procedure interface_mpi_file_open
  end interface

  interface interface_mpi_file_close
    module procedure interface_mpi_file_close
  end interface

  interface interface_mpi_file_set_view
    module procedure interface_mpi_file_set_view
  end interface

  interface interface_mpi_file_seek
    module procedure interface_mpi_file_seek
  end interface

  interface interface_mpi_file_write_at
    module procedure interface_mpi_file_write_at_i0
    module procedure interface_mpi_file_write_at_i1
    module procedure interface_mpi_file_write_at_r0
    module procedure interface_mpi_file_write_at_r1
  end interface

  interface interface_mpi_file_read_at
    module procedure interface_mpi_file_read_at_pck
    module procedure interface_mpi_file_read_at_i0
    module procedure interface_mpi_file_read_at_i1
    module procedure interface_mpi_file_read_at_r0
    module procedure interface_mpi_file_read_at_r1
  end interface

  interface interface_mpi_info_create
    module procedure interface_mpi_info_create
  end interface

  interface interface_mpi_info_free
    module procedure interface_mpi_info_free
  end interface

  interface interface_mpi_info_set
    module procedure interface_mpi_info_set
  end interface

  interface interface_mpi_get_count
    module procedure interface_mpi_get_count
  end interface

  interface interface_mpi_type_extent
    module procedure interface_mpi_type_extent
  end interface

  private

  integer,                        parameter, public :: global_communicator = mpi_comm_world
  integer,                        parameter, public :: self_communicator   = mpi_comm_self
  integer,                        parameter, public :: op_mpi_sum          = mpi_sum
  integer,                        parameter, public :: op_mpi_max          = mpi_max
  integer,                        parameter, public :: op_mpi_min          = mpi_min
  integer,                        parameter, public :: df_mpi_source       = mpi_source
  integer,                        parameter, public :: df_mpi_tag          = mpi_tag
  integer,                        parameter, public :: df_mpi_any_source   = mpi_any_source
  integer,                        parameter, public :: df_mpi_any_tag      = mpi_any_tag
  integer,                        parameter, public :: df_mpi_status_size  = mpi_status_size
  integer,                        parameter, public :: df_mpi_request_null = mpi_request_null
  integer,                        parameter, public :: df_mpi_mode_wronly  = mpi_mode_wronly 
  integer,                        parameter, public :: df_mpi_mode_rdonly  = mpi_mode_rdonly 
  integer,                        parameter, public :: df_mpi_mode_rdwr    = mpi_mode_rdwr   
  integer,                        parameter, public :: df_mpi_mode_create  = mpi_mode_create 
  integer,                        parameter, public :: df_mpi_mode_sequential = mpi_mode_sequential
  integer,                        parameter, public :: df_mpi_mode_excl    = mpi_mode_excl
  integer,                        parameter, public :: df_mpi_mode_delete_on_close  = mpi_mode_delete_on_close 
  integer,                        parameter, public :: df_mpi_file_null    = mpi_file_null
  integer,                        parameter, public :: df_mpi_seek_set     = mpi_seek_set
  integer,                        parameter, public :: df_mpi_seek_cur     = mpi_seek_cur
  integer,                        parameter, public :: df_mpi_packed       = mpi_packed
  integer,                        parameter, public :: df_mpi_info_null    = mpi_info_null
  integer,                        parameter, public :: df_mpi_byte         = mpi_byte
  integer,                        parameter, public :: df_mpi_real8        = mpi_real8
  integer,                        parameter, public :: df_mpi_complex8     = mpi_complex8
  integer,                        parameter, public :: df_mpi_integer      = mpi_integer
  integer,                        parameter, public :: df_mpi_integer8     = mpi_integer8
  integer,                        parameter, public :: df_mpi_logical      = mpi_logical
  integer,                        parameter, public :: df_my_mpi_integer   = my_MPI_INTEGER
  integer(kind=MPI_OFFSET_KIND),  parameter, public :: df_MPI_OFFSET_KIND  = MPI_OFFSET_KIND
  integer(kind=MPI_ADDRESS_KIND), parameter, public :: df_MPI_ADDRESS_KIND = MPI_ADDRESS_KIND

  save

contains

!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_i0(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
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

    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif
  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_i1(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:)
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
    call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_i2(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_i3(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:,:,:)
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
    call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_i4(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:,:,:,:)
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
    call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_r0(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x
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
    call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_r1(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:)
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
    call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_r2(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:)
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
    call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_r3(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:,:)
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
    call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_r4(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:,:,:)
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
    call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_cx(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    character*(*)              :: x
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
    call mpi_bcast(x, ndim_i4, MPI_CHARACTER, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , MPI_CHARACTER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_c1(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    character*(*)              :: x(:)
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
    call mpi_bcast(x, ndim_i4, MPI_CHARACTER, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , MPI_CHARACTER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_cmplx2(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    complex(8)                 :: x(:,:)
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
    call mpi_bcast(x, ndim_i4, MPI_COMPLEX8, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , MPI_COMPLEX8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_cmplx3(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    complex(8)                 :: x(:,:,:)
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
    call mpi_bcast(x, ndim_i4, MPI_COMPLEX8, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, ndim   , MPI_COMPLEX8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_probe(rank,message_id,communicator,status_container)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: status_container(MPI_STATUS_SIZE)
    integer(kind=integer_kind) :: rank
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: status_container_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: rank_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    communicator_i4 = communicator
    message_id_i4   = message_id
    rank_i4         = rank

    call mpi_probe(rank_i4, message_id_i4, communicator_i4, status_container_i4, ierr_i4)

    do i = 1, MPI_STATUS_SIZE
      status_container(i) = status_container_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_probe(rank   , message_id   , communicator   , status_container   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_iprobe_l0(rank,message_id,x,communicator,status_container)

!-------------------------------------------------------------------------------
    logical(kind=integer_kind) :: x
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: status_container(MPI_STATUS_SIZE)
    integer(kind=integer_kind) :: rank
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    logical(kind=4)            :: x_l4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: status_container_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: rank_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    communicator_i4 = communicator
    message_id_i4   = message_id
    rank_i4         = rank

    call mpi_iprobe(rank_i4, message_id_i4, communicator_i4, x_l4, status_container_i4, ierr_i4)

    do i = 1, MPI_STATUS_SIZE
      status_container(i) = status_container_i4(i)
    end do
    x               = x_l4
#else
    if (.not. MPI_INIT_called) return
    call mpi_iprobe(rank   , message_id   , communicator   , x   , status_container   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_test(request,x,status_container)

!-------------------------------------------------------------------------------
    logical(kind=integer_kind) :: x
    integer(kind=integer_kind) :: request
    integer(kind=integer_kind) :: status_container(MPI_STATUS_SIZE)
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    logical(kind=4)            :: x_l4
    integer(kind=4)            :: request_i4
    integer(kind=4)            :: status_container_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    request_i4      = request

    call mpi_test(request_i4, x_l4, status_container_i4, ierr_i4)

    do i = 1, MPI_STATUS_SIZE
      status_container(i) = status_container_i4(i)
    end do
    x               = x_l4
#else
    if (.not. MPI_INIT_called) return
    call mpi_test(request   , x   , status_container   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_testall(ndim,request,x,status_container)

!-------------------------------------------------------------------------------
    logical(kind=integer_kind) :: x
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: request(ndim)
    integer(kind=integer_kind) :: status_container(MPI_STATUS_SIZE,ndim)
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    logical(kind=4)            :: x_l4
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: request_i4(ndim)
    integer(kind=4)            :: status_container_i4(MPI_STATUS_SIZE,ndim)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i,j
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    do i = 1, ndim_i4
      request_i4(i) = request(i)
    end do

    call mpi_testall(ndim_i4, request_i4, x_l4, status_container_i4, ierr_i4)

    do i = 1, MPI_STATUS_SIZE
      do j = 1, ndim_i4
        status_container(i,j) = status_container_i4(i,j)
      end do
    end do
    x               = x_l4
#else
    if (.not. MPI_INIT_called) return
    call mpi_testall(ndim   , request   , x   , status_container   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_wait(request,status_container)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: request
    integer(kind=integer_kind) :: status_container(MPI_STATUS_SIZE)
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: request_i4
    integer(kind=4)            :: status_container_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    request_i4      = request

    call mpi_wait(request_i4, status_container_i4, ierr_i4)

    do i = 1, MPI_STATUS_SIZE
      status_container(i) = status_container_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_wait(request   , status_container   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! mpi_send interfaces
!-------------------------------------------------------------------------------

  subroutine interface_mpi_send_i0(x,ndim,receiver,message_id,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
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
    message_id_i4   = message_id
    receiver_i4     = receiver
    communicator_i4 = communicator
    call mpi_send(x, ndim_i4, my_MPI_INTEGER, receiver_i4, message_id_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_send(x, ndim   , my_MPI_INTEGER, receiver   , message_id   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_send_i1(x,ndim,receiver,message_id,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:)
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
    message_id_i4   = message_id
    receiver_i4     = receiver
    communicator_i4 = communicator
    call mpi_send(x, ndim_i4, my_MPI_INTEGER, receiver_i4, message_id_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_send(x, ndim   , my_MPI_INTEGER, receiver   , message_id   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_send_r0(x,ndim,receiver,message_id,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x
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
    message_id_i4   = message_id
    receiver_i4     = receiver
    communicator_i4 = communicator
    call mpi_send(x, ndim_i4, MPI_REAL8, receiver_i4, message_id_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_send(x, ndim   , MPI_REAL8, receiver   , message_id   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_send_r1(x,ndim,receiver,message_id,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:)
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
    message_id_i4   = message_id
    receiver_i4     = receiver
    communicator_i4 = communicator
    call mpi_send(x, ndim_i4, MPI_REAL8, receiver_i4, message_id_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_send(x, ndim   , MPI_REAL8, receiver   , message_id   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_send_r2(x,ndim,receiver,message_id,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:)
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
    message_id_i4   = message_id
    receiver_i4     = receiver
    communicator_i4 = communicator
    call mpi_send(x, ndim_i4, MPI_REAL8, receiver_i4, message_id_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_send(x, ndim   , MPI_REAL8, receiver   , message_id   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_send_cx(x,ndim,receiver,message_id,communicator)

!-------------------------------------------------------------------------------
    character*(*)              :: x
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
    message_id_i4   = message_id
    receiver_i4     = receiver
    communicator_i4 = communicator
    call mpi_send(x, ndim_i4, MPI_CHARACTER, receiver_i4, message_id_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_send(x, ndim   , MPI_CHARACTER, receiver   , message_id   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_isend_i0(x,ndim,receiver,message_id,communicator,request)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
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
    message_id_i4   = message_id
    receiver_i4     = receiver
    communicator_i4 = communicator
    call mpi_isend(x, ndim_i4, my_MPI_INTEGER, receiver_i4, message_id_i4, communicator_i4, request_i4, ierr_i4)
    request         = request_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_isend(x, ndim   , my_MPI_INTEGER, receiver   , message_id   , communicator   , request   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_isend_i1(x,ndim,receiver,message_id,communicator,request)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:)
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
    message_id_i4   = message_id
    receiver_i4     = receiver
    communicator_i4 = communicator
    call mpi_isend(x, ndim_i4, my_MPI_INTEGER, receiver_i4, message_id_i4, communicator_i4, request_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_isend(x, ndim   , my_MPI_INTEGER, receiver   , message_id   , communicator   , request   , ierr   )
#endif

  end subroutine

!-------------------------------------------------------------------------------
! mpi_recv interfaces
!-------------------------------------------------------------------------------

  subroutine interface_mpi_recv_i0(x,ndim,sender,message_id,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_recv(x, ndim_i4, my_MPI_INTEGER, sender_i4, message_id_i4, communicator_i4, status_array_i4, ierr_i4)
#else
    integer(kind=integer_kind) :: status_array(MPI_STATUS_SIZE)
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , my_MPI_INTEGER, sender   , message_id   , communicator   , status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_recv_i0_stat(x,ndim,sender,message_id,communicator,status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4, i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_recv(x, ndim_i4, my_MPI_INTEGER, sender_i4, message_id_i4, communicator_i4, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , my_MPI_INTEGER, sender   , message_id   , communicator   , status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_recv_i1(x,ndim,sender,message_id,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_recv(x, ndim_i4, my_MPI_INTEGER, sender_i4, message_id_i4, communicator_i4, status_array_i4, ierr_i4)
#else
    integer(kind=integer_kind) :: status_array(MPI_STATUS_SIZE)
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , my_MPI_INTEGER, sender   , message_id   , communicator   , status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------
  subroutine interface_mpi_recv_i1_stat(x,ndim,sender,message_id,communicator,status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4, i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_recv(x, ndim_i4, my_MPI_INTEGER, sender_i4, message_id_i4, communicator_i4, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , my_MPI_INTEGER, sender   , message_id   , communicator   , status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_recv_r0(x,ndim,sender,message_id,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_recv(x, ndim_i4, MPI_REAL8, sender_i4, message_id_i4, communicator_i4, status_array_i4, ierr_i4)
#else
    integer(kind=integer_kind) :: status_array(MPI_STATUS_SIZE)
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , MPI_REAL8, sender   , message_id   , communicator   , status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_recv_r1(x,ndim,sender,message_id,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_recv(x, ndim_i4, MPI_REAL8, sender_i4, message_id_i4, communicator_i4, status_array_i4, ierr_i4)
#else
    integer(kind=integer_kind) :: status_array(MPI_STATUS_SIZE)
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , MPI_REAL8, sender   , message_id   , communicator   , status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_recv_r1_stat(x,ndim,sender,message_id,communicator,status_array)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4, i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_recv(x, ndim_i4, MPI_REAL8, sender_i4, message_id_i4, communicator_i4, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , MPI_REAL8, sender   , message_id   , communicator   , status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_recv_r2(x,ndim,sender,message_id,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_recv(x, ndim_i4, MPI_REAL8, sender_i4, message_id_i4, communicator_i4, status_array_i4, ierr_i4)
#else
    integer(kind=integer_kind) :: status_array(MPI_STATUS_SIZE)
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , MPI_REAL8, sender   , message_id   , communicator   , status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_recv_cx(x,ndim,sender,message_id,communicator)

!-------------------------------------------------------------------------------
    character*(*)              :: x
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_recv(x, ndim_i4, MPI_CHARACTER, sender_i4, message_id_i4, communicator_i4, status_array_i4, ierr_i4)
#else
    integer(kind=integer_kind) :: status_array(MPI_STATUS_SIZE)
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , MPI_CHARACTER, sender   , message_id   , communicator   , status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_irecv_i0(x,ndim,sender,message_id,communicator,request)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: request
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: request_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_irecv(x, ndim_i4, my_MPI_INTEGER, sender_i4, message_id_i4, communicator_i4, request_i4, ierr_i4)
    request         = request_i4
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_irecv(x, ndim   , my_MPI_INTEGER, sender   , message_id   , communicator   , request   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_irecv_i1(x,ndim,sender,message_id,communicator,request)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: sender
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: request
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: sender_i4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: request_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    message_id_i4   = message_id
    sender_i4       = sender
    communicator_i4 = communicator
    call mpi_irecv(x, ndim_i4, my_MPI_INTEGER, sender_i4, message_id_i4, communicator_i4, request_i4, ierr_i4)
    request         = request_i4
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_irecv(x, ndim   , my_MPI_INTEGER, sender   , message_id   , communicator   , request   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! mpi_reduce interfaces
!-------------------------------------------------------------------------------

  subroutine interface_mpi_reduce_i0(x,y,ndim,operation,receiver,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
    integer(kind=integer_kind) :: y
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
    call mpi_reduce(x, y, ndim_i4, my_MPI_INTEGER, operation_i4, receiver_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_reduce(x, y, ndim   , my_MPI_INTEGER, operation   , receiver   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_reduce_i1(x,y,ndim,operation,receiver,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:)
    integer(kind=integer_kind) :: y(:)
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
    call mpi_reduce(x, y, ndim_i4, my_MPI_INTEGER, operation_i4, receiver_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_reduce(x, y, ndim   , my_MPI_INTEGER, operation   , receiver   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_reduce_i2(x,y,ndim,operation,receiver,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:,:)
    integer(kind=integer_kind) :: y(:,:)
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
    call mpi_reduce(x, y, ndim_i4, my_MPI_INTEGER, operation_i4, receiver_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_reduce(x, y, ndim   , my_MPI_INTEGER, operation   , receiver   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_reduce_r0(x,y,ndim,operation,receiver,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x
    real(8)                    :: y
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
    call mpi_reduce(x, y, ndim_i4, MPI_REAL8, operation_i4, receiver_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_reduce(x, y, ndim   , MPI_REAL8, operation   , receiver   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_reduce_r1(x,y,ndim,operation,receiver,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:)
    real(8)                    :: y(:)
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
    call mpi_reduce(x, y, ndim_i4, MPI_REAL8, operation_i4, receiver_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_reduce(x, y, ndim   , MPI_REAL8, operation   , receiver   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_reduce_r2(x,y,ndim,operation,receiver,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:)
    real(8)                    :: y(:,:)
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
    call mpi_reduce(x, y, ndim_i4, MPI_REAL8, operation_i4, receiver_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_reduce(x, y, ndim   , MPI_REAL8, operation   , receiver   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_allreduce_i0(x,y,ndim,operation,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
    integer(kind=integer_kind) :: y
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
    call mpi_allreduce(x, y, ndim_i4, my_MPI_INTEGER, operation_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_allreduce(x, y, ndim   , my_MPI_INTEGER, operation   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_allreduce_i1(x,y,ndim,operation,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:)
    integer(kind=integer_kind) :: y(:)
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
    call mpi_allreduce(x, y, ndim_i4, my_MPI_INTEGER, operation_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_allreduce(x, y, ndim   , my_MPI_INTEGER, operation   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

#ifndef INT_STAR8
  subroutine interface_mpi_allreduce_i0_i8(x,y,ndim,operation,communicator)

!-------------------------------------------------------------------------------
    integer(kind=8)            :: x
    integer(kind=8)            :: y
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
    call mpi_allreduce(x, y, ndim_i4, MPI_INTEGER8, operation_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_allreduce(x, y, ndim   , MPI_INTEGER8, operation   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------
#endif

  subroutine interface_mpi_allreduce_r0(x,y,ndim,operation,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x
    real(8)                    :: y
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
    call mpi_allreduce(x, y, ndim_i4, MPI_REAL8, operation_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_allreduce(x, y, ndim   , MPI_REAL8, operation   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_allreduce_r1(x,y,ndim,operation,communicator)

!-------------------------------------------------------------------------------
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
    call mpi_allreduce(x, y, ndim_i4, MPI_REAL8, operation_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_allreduce(x, y, ndim   , MPI_REAL8, operation   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_allreduce_r2(x,y,ndim,operation,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:)
    real(8)                    :: y(:,:)
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
    call mpi_allreduce(x, y, ndim_i4, MPI_REAL8, operation_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_allreduce(x, y, ndim   , MPI_REAL8, operation   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_allreduce_r3(x,y,ndim,operation,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:)
    real(8)                    :: y(:)
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
    call mpi_allreduce(x, y, ndim_i4, MPI_REAL8, operation_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_allreduce(x, y, ndim   , MPI_REAL8, operation   , communicator   , ierr   )
#endif

  end subroutine

!-------------------------------------------------------------------------------
  subroutine interface_mpi_scatter_r1(x,ndim1,y,ndim2,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:)
    real(8)                    :: y(:)
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
    call mpi_scatter(x, ndim1_i4, MPI_REAL8, y, ndim2_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_scatter(x, ndim1   , MPI_REAL8, y, ndim2   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! mpi_gather interfaces
!-------------------------------------------------------------------------------

  subroutine interface_mpi_gather_i(x,ndim,y,disp,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
    integer(kind=integer_kind) :: y(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: disp
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: disp_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    disp_i4         = disp
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_gather(x, ndim_i4, my_MPI_INTEGER, y, disp_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_gather(x, ndim   , my_MPI_INTEGER, y, disp   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_gather_r(x,ndim,y,disp,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x
    real(8)                    :: y(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: disp
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: disp_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    disp_i4         = disp
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_gather(x, ndim_i4, MPI_REAL8, y, disp_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_gather(x, ndim   , MPI_REAL8, y, disp   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

    subroutine interface_mpi_gatherv(x,ndim,y,counter,disp,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(*)
    real(8)                    :: y(*)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: counter(:)
    integer(kind=integer_kind) :: disp(:)
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4) :: ndim_i4
    integer(kind=4) :: root_proc_i4
    integer(kind=4) :: communicator_i4
    integer(kind=4) :: ierr_i4
    integer(kind=4), allocatable :: disp_i4(:)
    integer(kind=4), allocatable  :: counter_i4(:)
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    allocate (disp_i4(size(disp)))
    allocate (counter_i4(size(counter)))
    ndim_i4         = ndim
    disp_i4         = disp
    counter_i4      = counter
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call MPI_Gatherv(x,ndim_i4,MPI_REAL8,y,counter_i4,disp_i4,MPI_REAL8,root_proc_i4,communicator_i4,ierr_i4)
    deallocate (disp_i4)
    deallocate (counter_i4)
#else
    if (.not. MPI_INIT_called) return
    call MPI_Gatherv(x, ndim   , MPI_REAL8, y,counter,disp,MPI_REAL8,root_proc,communicator,ierr )
#endif

    end subroutine

!-------------------------------------------------------------------------------


    subroutine interface_mpi_allgatherv(x,ndim,y,counter,disp,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(*)
    real(8)                    :: y(*)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: counter(:)
    integer(kind=integer_kind) :: disp(:)
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4) :: ndim_i4
    integer(kind=4) :: communicator_i4
    integer(kind=4) :: ierr_i4
    integer(kind=4), allocatable :: disp_i4(:)
    integer(kind=4), allocatable  :: counter_i4(:)
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    allocate (disp_i4(size(disp)))
    allocate (counter_i4(size(counter)))
    ndim_i4         = ndim
    disp_i4         = disp
    counter_i4      = counter
    communicator_i4 = communicator
    call MPI_AllGatherv(x,ndim_i4,MPI_REAL8,y,counter_i4,disp_i4,MPI_REAL8,communicator_i4,ierr_i4)
    deallocate (disp_i4)
    deallocate (counter_i4)
#else
    if (.not. MPI_INIT_called) return
    call MPI_AllGatherv(x, ndim   , MPI_REAL8, y,counter,disp,MPI_REAL8,communicator,ierr )
#endif

    end subroutine

!--------------------------------------------------------------------------------------------------------------


  subroutine interface_mpi_allgather_i0(x,ndim,y,ndim2,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
    integer(kind=integer_kind) :: y(*)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: ndim2
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: ndim2_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    ndim2_i4        = ndim2
    communicator_i4 = communicator
    call mpi_allgather(x, ndim_i4, my_MPI_INTEGER, y, ndim2_i4, my_MPI_INTEGER, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_allgather(x, ndim   , my_MPI_INTEGER, y, ndim2   , my_MPI_INTEGER, communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------
  subroutine interface_mpi_allgather_i(x,ndim,y,ndim2,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(*)
    integer(kind=integer_kind) :: y(*)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: ndim2
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: ndim2_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    ndim2_i4        = ndim2
    communicator_i4 = communicator
    call mpi_allgather(x, ndim_i4, my_MPI_INTEGER, y, ndim2_i4, my_MPI_INTEGER, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_allgather(x, ndim   , my_MPI_INTEGER, y, ndim2   , my_MPI_INTEGER, communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_allgather_r(x,ndim,y,ndim2,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(*)
    real(8)                    :: y(*)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: ndim2
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: ndim2_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim_i4         = ndim
    ndim2_i4        = ndim2
    communicator_i4 = communicator
    call mpi_allgather(x, ndim_i4, MPI_REAL8, y, ndim2_i4, MPI_REAL8, communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_allgather(x, ndim   , MPI_REAL8, y, ndim2   , MPI_REAL8, communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! interfaces to non-buffer MPI routines
!-------------------------------------------------------------------------------

  subroutine interface_mpi_comm_rank(communicator,rank)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: rank
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: rank_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    communicator_i4 = communicator
    call mpi_comm_rank(communicator_i4, rank_i4, ierr_i4)
    rank            = rank_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_comm_rank(communicator   , rank   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_comm_size(communicator,is_size)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: is_size
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: is_size_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    communicator_i4 = communicator
    call mpi_comm_size(communicator_i4, is_size_i4, ierr_i4)
    is_size         = is_size_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_comm_size(communicator   , is_size      , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_comm_split(old_communicator, color, key, new_communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: old_communicator
    integer(kind=integer_kind) :: new_communicator
    integer(kind=integer_kind) :: color
    integer(kind=integer_kind) :: key
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: old_communicator_i4
    integer(kind=4)            :: new_communicator_i4
    integer(kind=4)            :: color_i4
    integer(kind=4)            :: key_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    old_communicator_i4 = old_communicator
    color_i4            = color
    key_i4              = key
    call mpi_comm_split(old_communicator_i4, color_i4, key_i4, new_communicator_i4, ierr_i4)
    new_communicator    = new_communicator_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_comm_split(old_communicator   , color   , key   , new_communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_comm_free(communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    communicator_i4 = communicator
    call mpi_comm_free(communicator_i4, ierr_i4)
    communicator    = communicator_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_comm_free(communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_barrier(communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    communicator_i4 = communicator
    call mpi_barrier(communicator_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_barrier(communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_get_processor_name(pname,ndim)

!-------------------------------------------------------------------------------
    character(len=255)         :: pname ! len = 255 corresponds to MPI_MAX_PROCESSOR_NAME
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_get_processor_name(pname, ndim_i4, ierr_i4)
    ndim = ndim_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_get_processor_name(pname, ndim   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_init()

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: ierr
    integer(kind=integer_kind) :: mpi_th_provided
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: mpi_th_provided_i4
!-------------------------------------------------------------------------------
!caar    call mpi_init(ierr_i4)
    call mpi_init_thread(MPI_THREAD_MULTIPLE,mpi_th_provided_i4,ierr_i4)
    if (mpi_th_provided_i4 /= MPI_THREAD_MULTIPLE) print*, 'warning: MPI library has no thread support'
#else
!caar    call mpi_init(ierr   )
    call mpi_init_thread(MPI_THREAD_MULTIPLE,mpi_th_provided,ierr)
    if (mpi_th_provided /= MPI_THREAD_MULTIPLE) print*, 'warning: MPI library has no thread support'
#endif
    MPI_INIT_called = .true.

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_finalize()

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    call mpi_finalize(ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_finalize(ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_abort(communicator, error_code)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: error_code
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: error_code_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    communicator_i4 = communicator
    error_code_i4   = error_code
    call mpi_abort(communicator_i4, error_code_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_abort(communicator   , error_code   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_open(communicator, fname, amode, info, file_handle)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=integer_kind)       :: info
    integer(kind=integer_kind)       :: amode
    integer(kind=integer_kind)       :: communicator
    character*(*)                    :: fname
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: info_i4
    integer(kind=4)            :: amode_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    info_i4         = info
    amode_i4        = amode
    communicator_i4 = communicator
    call mpi_file_open(communicator_i4, fname, amode_i4, info_i4, file_handle_i4, ierr_i4)
    file_handle     = file_handle_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_open(communicator   , fname, amode   , info   , file_handle   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_close(file_handle)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4 = file_handle
    call mpi_file_close(file_handle_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_close(file_handle   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_set_view(file_handle, disp, etype, ftype, datarep, info)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=integer_kind)       :: info
    integer(kind=integer_kind)       :: etype
    integer(kind=integer_kind)       :: ftype
    character*(*)                    :: datarep
    integer(kind=df_mpi_offset_kind) :: disp
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: info_i4
    integer(kind=4)            :: etype_i4
    integer(kind=4)            :: ftype_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    info_i4         = info
    etype_i4        = etype
    ftype_i4        = ftype
    file_handle_i4  = file_handle
    call mpi_file_set_view(file_handle_i4, disp, etype_i4, ftype_i4, datarep, info_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_set_view(file_handle   , disp, etype   , ftype   , datarep, info   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_seek(file_handle, offset, whence)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: whence
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: whence_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4 = file_handle
    whence_i4      = whence
    call mpi_file_seek(file_handle_i4, offset, whence_i4, ierr_i4)
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_seek(file_handle   , offset, whence   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_iwrite(file_handle, buf, length, request)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    character*(*)                    :: buf(:)
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: request
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: request_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_iwrite(file_handle_i4, buf, length_i4,    mpi_packed, request_i4, ierr_i4)
    request         = request_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_iwrite(file_handle   , buf, length   , df_mpi_packed, request   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_iwrite_at_pck(file_handle, offset, buf, length, request)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    character*(*)                    :: buf(:)
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: request
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: request_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_iwrite_at(file_handle_i4, offset, buf, length_i4,    mpi_packed, request_i4, ierr_i4)
    request         = request_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_iwrite_at(file_handle   , offset, buf, length   , df_mpi_packed, request   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_write_at_i0(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=integer_kind)       :: buf
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_write_at(file_handle_i4, offset, buf, length_i4, my_MPI_INTEGER, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_write_at(file_handle   , offset, buf, length   , my_MPI_INTEGER, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_write_at_i1(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=integer_kind)       :: buf(:)
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_write_at(file_handle_i4, offset, buf, length_i4, my_MPI_INTEGER, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_write_at(file_handle   , offset, buf, length   , my_MPI_INTEGER, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_write_at_r0(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    real(8)                          :: buf
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_write_at(file_handle_i4, offset, buf, length_i4, MPI_REAL8, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_write_at(file_handle   , offset, buf, length   , MPI_REAL8, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_write_at_r1(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    real(8)                          :: buf(:)
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_write_at(file_handle_i4, offset, buf, length_i4, MPI_REAL8, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_write_at(file_handle   , offset, buf, length   , MPI_REAL8, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_write_at_i(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=integer_kind)       :: buf
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_write_at(file_handle_i4, offset, buf, length_i4, my_MPI_INTEGER, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_write_at(file_handle   , offset, buf, length   , my_MPI_INTEGER, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_write_at_r(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    real(8)                          :: buf(*)
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_write_at(file_handle_i4, offset, buf, length_i4, MPI_REAL8, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_write_at(file_handle   , offset, buf, length   , MPI_REAL8, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_write_at_bi(file_handle, offset, buf, length, byte_tp, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=df_mpi_offset_kind) :: buf
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: byte_tp
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_write_at(file_handle_i4, offset, buf, length_i4, MPI_BYTE, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_write_at(file_handle   , offset, buf, length   , MPI_BYTE, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_write_at_br(file_handle, offset, buf, length, byte_tp, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    real(8)                          :: buf(*)
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: byte_tp
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_write_at(file_handle_i4, offset, buf, length_i4, MPI_BYTE, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_write_at(file_handle   , offset, buf, length   , MPI_BYTE, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_read_pck(file_handle, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    character*(*)                    :: buf(:)
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_read(file_handle_i4, buf, length_i4,    mpi_packed, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_read(file_handle   , buf, length   , df_mpi_packed, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_read_at_pck(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    character*(*)                    :: buf(:)
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_read_at(file_handle_i4, offset, buf, length_i4,    mpi_packed, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_read_at(file_handle   , offset, buf, length   , df_mpi_packed, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_read_at_i0(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=integer_kind)       :: buf
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_read_at(file_handle_i4, offset, buf, length_i4, my_MPI_INTEGER, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_read_at(file_handle   , offset, buf, length   , my_MPI_INTEGER, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_read_at_i1(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=integer_kind)       :: buf(:)
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_read_at(file_handle_i4, offset, buf, length_i4, my_MPI_INTEGER, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_read_at(file_handle   , offset, buf, length   , my_MPI_INTEGER, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_read_at_r0(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    real(8)                          :: buf
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_read_at(file_handle_i4, offset, buf, length_i4, MPI_REAL8, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_read_at(file_handle   , offset, buf, length   , MPI_REAL8, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_read_at_r1(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    real(8)                          :: buf(:)
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_read_at(file_handle_i4, offset, buf, length_i4, MPI_REAL8, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_read_at(file_handle   , offset, buf, length   , MPI_REAL8, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_read_at_i(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=integer_kind)       :: buf
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_read_at(file_handle_i4, offset, buf, length_i4, my_MPI_INTEGER, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_read_at(file_handle   , offset, buf, length   , my_MPI_INTEGER, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_read_at_r(file_handle, offset, buf, length, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    real(8)                          :: buf(*)
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_read_at(file_handle_i4, offset, buf, length_i4, MPI_REAL8, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_read_at(file_handle   , offset, buf, length   , MPI_REAL8, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_read_at_bi(file_handle, offset, buf, length, byte_tp, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    integer(kind=df_mpi_offset_kind) :: buf
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: byte_tp
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_read_at(file_handle_i4, offset, buf, length_i4, MPI_BYTE, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_read_at(file_handle   , offset, buf, length   , MPI_BYTE, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_file_read_at_br(file_handle, offset, buf, length, byte_tp, status_array)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: file_handle
    real(8)                          :: buf(*)
    integer(kind=df_mpi_offset_kind) :: offset
    integer(kind=integer_kind)       :: length
    integer(kind=integer_kind)       :: byte_tp
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: file_handle_i4
    integer(kind=4)            :: length_i4
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    file_handle_i4  = file_handle
    length_i4       = length
    call mpi_file_read_at(file_handle_i4, offset, buf, length_i4, MPI_BYTE, status_array_i4, ierr_i4)
    do i = 1, MPI_STATUS_SIZE
      status_array(i) = status_array_i4(i)
    end do
#else
    if (.not. MPI_INIT_called) return
    call mpi_file_read_at(file_handle   , offset, buf, length   , MPI_BYTE, status_array   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_info_create(info_object)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: info_object
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)                  :: info_object_i4
    integer(kind=4)                  :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    info_object_i4 = info_object
    call mpi_info_create(info_object_i4, ierr_i4)
    info_object    = info_object_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_info_create(info_object   , ierr  )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_info_free(info_object)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: info_object
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)                  :: info_object_i4
    integer(kind=4)                  :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    info_object_i4 = info_object
    call mpi_info_free(info_object_i4, ierr_i4)
    info_object    = info_object_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_info_free(info_object   , ierr  )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_info_set(info_object, key, value)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: info_object
    character*(*)                    :: key
    character*(*)                    :: value
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)                  :: info_object_i4
    integer(kind=4)                  :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    info_object_i4 = info_object
    call mpi_info_set(info_object_i4, key, value, ierr_i4)
    info_object    = info_object_i4
#else
    integer(kind=integer_kind)       :: ierr
    if (.not. MPI_INIT_called) return
    call mpi_info_set(info_object   , key, value, ierr  )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_pack_i0(x, ndim1, y, ndim2, pos, communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: x
    integer(kind=integer_kind)       :: ndim1
    character*(*)                    :: y(:)
    integer(kind=integer_kind)       :: ndim2
    integer(kind=integer_kind)       :: pos
    integer(kind=integer_kind)       :: communicator
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim1_i4
    integer(kind=4)            :: ndim2_i4
    integer(kind=4)            :: pos_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim1_i4        = ndim1
    ndim2_i4        = ndim2
    pos_i4          = pos
    communicator_i4 = communicator
    call mpi_pack(x, ndim1_i4, my_MPI_INTEGER, y, ndim2_i4, pos_i4, communicator_i4, ierr_i4)
    pos             = pos_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_pack(x, ndim1   , my_MPI_INTEGER, y, ndim2   , pos   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_pack_i2(x, ndim1, y, ndim2, pos, communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: x(:,:)
    integer(kind=integer_kind)       :: ndim1
    character*(*)                    :: y(:)
    integer(kind=integer_kind)       :: ndim2
    integer(kind=integer_kind)       :: pos
    integer(kind=integer_kind)       :: communicator
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim1_i4
    integer(kind=4)            :: ndim2_i4
    integer(kind=4)            :: pos_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim1_i4        = ndim1
    ndim2_i4        = ndim2
    pos_i4          = pos
    communicator_i4 = communicator
    call mpi_pack(x, ndim1_i4, my_MPI_INTEGER, y, ndim2_i4, pos_i4, communicator_i4, ierr_i4)
    pos             = pos_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_pack(x, ndim1   , my_MPI_INTEGER, y, ndim2   , pos   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_pack_r(x, ndim1, y, ndim2, pos, communicator)

!-------------------------------------------------------------------------------
    real(8)                          :: x(:)
    integer(kind=integer_kind)       :: ndim1
    character*(*)                    :: y(:)
    integer(kind=integer_kind)       :: ndim2
    integer(kind=integer_kind)       :: pos
    integer(kind=integer_kind)       :: communicator
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim1_i4
    integer(kind=4)            :: ndim2_i4
    integer(kind=4)            :: pos_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim1_i4        = ndim1
    ndim2_i4        = ndim2
    pos_i4          = pos
    communicator_i4 = communicator
    call mpi_pack(x, ndim1_i4, MPI_REAL8, y, ndim2_i4, pos_i4, communicator_i4, ierr_i4)
    pos             = pos_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_pack(x, ndim1   , MPI_REAL8, y, ndim2   , pos   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_unpack_i0(x, ndim1, pos, y, ndim2, communicator)

!-------------------------------------------------------------------------------
    character*(*)                    :: x(:)
    integer(kind=integer_kind)       :: ndim1
    integer(kind=integer_kind)       :: y
    integer(kind=integer_kind)       :: ndim2
    integer(kind=integer_kind)       :: pos
    integer(kind=integer_kind)       :: communicator
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim1_i4
    integer(kind=4)            :: ndim2_i4
    integer(kind=4)            :: pos_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim1_i4        = ndim1
    ndim2_i4        = ndim2
    pos_i4          = pos
    communicator_i4 = communicator
    call mpi_unpack(x, ndim1_i4, pos_i4, y, ndim2_i4, my_MPI_INTEGER, communicator_i4, ierr_i4)
    pos             = pos_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_unpack(x, ndim1   , pos   , y, ndim2   , my_MPI_INTEGER, communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_unpack_i2(x, ndim1, pos, y, ndim2, communicator)

!-------------------------------------------------------------------------------
    character*(*)                    :: x(:)
    integer(kind=integer_kind)       :: ndim1
    integer(kind=integer_kind)       :: y(:,:)
    integer(kind=integer_kind)       :: ndim2
    integer(kind=integer_kind)       :: pos
    integer(kind=integer_kind)       :: communicator
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim1_i4
    integer(kind=4)            :: ndim2_i4
    integer(kind=4)            :: pos_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim1_i4        = ndim1
    ndim2_i4        = ndim2
    pos_i4          = pos
    communicator_i4 = communicator
    call mpi_unpack(x, ndim1_i4, pos_i4, y, ndim2_i4, my_MPI_INTEGER, communicator_i4, ierr_i4)
    pos             = pos_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_unpack(x, ndim1   , pos   , y, ndim2   , my_MPI_INTEGER, communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_unpack_r(x, ndim1, pos, y, ndim2, communicator)

!-------------------------------------------------------------------------------
    character*(*)                    :: x(:)
    integer(kind=integer_kind)       :: ndim1
    real(8)                          :: y(:)
    integer(kind=integer_kind)       :: ndim2
    integer(kind=integer_kind)       :: pos
    integer(kind=integer_kind)       :: communicator
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: ndim1_i4
    integer(kind=4)            :: ndim2_i4
    integer(kind=4)            :: pos_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    ndim1_i4        = ndim1
    ndim2_i4        = ndim2
    pos_i4          = pos
    communicator_i4 = communicator
    call mpi_unpack(x, ndim1_i4, pos_i4, y, ndim2_i4, MPI_REAL8, communicator_i4, ierr_i4)
    pos             = pos_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_unpack(x, ndim1   , pos   , y, ndim2   , MPI_REAL8, communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_get_count(status_array, datatype, elements)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)
    integer(kind=integer_kind)       :: datatype
    integer(kind=integer_kind)       :: elements
    integer(kind=integer_kind)       :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: status_array_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: datatype_i4
    integer(kind=4)            :: elements_i4
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    if (.not. MPI_INIT_called) return
    datatype_i4     = datatype
    do i = 1, MPI_STATUS_SIZE
      status_array_i4(i) = status_array(i)
    end do
    call mpi_get_count(status_array_i4, datatype_i4, elements_i4, ierr_i4)
    elements        = elements_i4
#else
    if (.not. MPI_INIT_called) return
    call mpi_get_count(status_array   , datatype   , elements   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_type_extent(datatype, extent)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind)       :: datatype
    integer(kind=mpi_address_kind)   :: lb, extent 
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4)            :: datatype_i4
    integer(kind=4)            :: ierr_i4
    integer(kind=mpi_address_kind)   :: lb_i4, extent_i4 
!-------------------------------------------------------------------------------
    datatype_i4     = datatype
    lb_i4 = extent 
#if interface_to_mpi_DEBUG > 0
    write (0,*) 'datatype ',datatype
    write (0,*) 'mpi_type_extent called, MPI_INIT ?',MPI_INIT_called
#endif
    if (.not. MPI_INIT_called) return
    call mpi_type_get_extent(datatype_i4, lb_i4, extent_i4, ierr_i4)
    extent          = extent_i4
#else
    integer(kind=integer_kind)       :: ierr
#if interface_to_mpi_DEBUG > 0
    write (0,*) 'datatype ',datatype
    write (0,*) 'mpi_type_extent called, MPI_INIT ?',MPI_INIT_called
#endif
    if (.not. MPI_INIT_called) return
    call mpi_type_get_extent(datatype   , lb, extent   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  real(8) function interface_MPI_WTIME()

!-------------------------------------------------------------------------------
    interface_MPI_WTIME = MPI_WTIME()
  end function
!-------------------------------------------------------------------------------

#else
  implicit none
  integer :: my_dummy_interface_to_mpi
#endif
end module interface_to_mpi
