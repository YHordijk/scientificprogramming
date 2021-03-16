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

module xc_mpi

#ifdef VAR_MPI

   use interface_to_mpi
   use integer_kind_mpilib
#ifdef USE_MPI_MOD_F90
   use mpi
   implicit none
#else
   implicit none
#include "mpif.h"
#endif

#else

   implicit none

#endif

   public xc_mpi_is_master
   public xc_mpi_get_nr_proc
#ifdef VAR_MPI
   public xc_mpi_bcast
   public xc_mpi_reduce
   public xc_mpi_distribute_points
   public xc_mpi_get_rank
   public xc_interface_mpi_barrier
   public xc_interface_mpi_file_open_to_write
   public xc_interface_mpi_file_close
   public xc_interface_mpi_file_set_view_r8
   public xc_interface_mpi_file_write_at
   public xc_interface_mpi_info_create
   public xc_interface_mpi_info_free
#endif

   private

#ifdef VAR_MPI
   interface xc_mpi_get_rank
      module procedure get_my_rank
   end interface
   interface xc_mpi_bcast
      module procedure xc_mpi_bcast_i0
      module procedure xc_mpi_bcast_i1
      module procedure xc_mpi_bcast_i2
      module procedure xc_mpi_bcast_i3
      module procedure xc_mpi_bcast_i4
      module procedure xc_mpi_bcast_r0
      module procedure xc_mpi_bcast_r1
      module procedure xc_mpi_bcast_r2
      module procedure xc_mpi_bcast_r3
      module procedure xc_mpi_bcast_r4
      module procedure xc_mpi_bcast_l0
      module procedure xc_mpi_bcast_l1
      module procedure xc_mpi_bcast_cx
      module procedure xc_mpi_bcast_c1
   end interface

   interface xc_mpi_reduce
      module procedure xc_mpi_reduce_i0
      module procedure xc_mpi_reduce_i1
      module procedure xc_mpi_reduce_r0
      module procedure xc_mpi_reduce_r1
      module procedure xc_mpi_reduce_r2
      module procedure xc_mpi_reduce_r4
   end interface
#endif

#ifdef VAR_MPI_32BIT_INT
   integer(4)         :: ierr
#else
   integer            :: ierr
#endif
   integer, parameter :: irank_master = 0

contains

#ifdef VAR_MPI
   function get_my_rank()
      integer :: get_my_rank
      integer :: irank
      call interface_mpi_comm_rank(global_communicator, irank)
      get_my_rank = irank
   end function
#endif

   function xc_mpi_is_master()
      logical :: xc_mpi_is_master
#include "dcbgen.h"
      xc_mpi_is_master = .true.
#ifdef VAR_MPI
      if (parcal) then
         xc_mpi_is_master = (get_my_rank() == irank_master)
      end if
#endif
   end function

   function xc_mpi_get_nr_proc()
      integer :: xc_mpi_get_nr_proc
      integer :: nr_proc
#include "dcbgen.h"
      xc_mpi_get_nr_proc = 1
#ifdef VAR_MPI
      if (parcal) then
         call interface_mpi_comm_size(global_communicator, nr_proc)
         xc_mpi_get_nr_proc = nr_proc
      end if
#endif
   end function

#ifdef VAR_MPI
   subroutine xc_mpi_bcast_i0(x)
      integer :: x
      call interface_mpi_bcast(x, 1, irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_i1(x)
      integer :: x(:)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_i2(x)
      integer :: x(:, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_i3(x)
      integer :: x(:, :, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_i4(x)
      integer :: x(:, :, :, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_r0(x)
      real(8) :: x
      call interface_mpi_bcast(x, 1,  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_r1(x)
      real(8) :: x(:)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_r2(x)
      real(8) :: x(:, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_r3(x)
      real(8) :: x(:, :, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_r4(x)
      real(8) :: x(:, :, :, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_l0(x)
      logical :: x
      call interface_mpi_bcast_l0(x, 1, irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_l1(x)
      logical :: x(:)
      call interface_mpi_bcast_l1(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_cx(x)
      character*(*) :: x
      call interface_mpi_bcast(x, 1,  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_bcast_c1(x)
      character*(*) :: x(:)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine xc_mpi_reduce_i0(x)
      integer :: x
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, 1, mpi_integer, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, 1, mpi_integer, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine xc_mpi_reduce_i1(x)
      integer :: x(:)
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_integer, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_integer, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine xc_mpi_reduce_r0(x)
      real(8) :: x
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, 1, mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, 1, mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine xc_mpi_reduce_r1(x)
      real(8) :: x(:)
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine xc_mpi_reduce_r2(x)
      real(8) :: x(:, :)
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine xc_mpi_reduce_r4(x)
      real(8) :: x(:, :, :, :)
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine xc_mpi_scatter_r1(x, n)
      real(8)              :: x(:)
      real(8), allocatable :: buffer(:)
      integer              :: n
   end subroutine

   subroutine xc_mpi_distribute_points(rx, ry, rz, rw, nr_points_batch, nr_points_on_this_proc)
      real(8)              :: rx(*)
      real(8)              :: ry(*)
      real(8)              :: rz(*)
      real(8)              :: rw(*)
      integer              :: nr_points_batch
      integer              :: nr_points_on_this_proc
      real(8), allocatable :: buffer(:)
      integer              :: nr_points_rest
      integer              :: i
      integer, allocatable :: dist_overview(:)
      integer              :: isum

      nr_points_on_this_proc = nr_points_batch/xc_mpi_get_nr_proc()
      if (xc_mpi_is_master()) then
         allocate(buffer(nr_points_on_this_proc))
         call interface_mpi_scatter_r1_work_f77(rx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(ry, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(rz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(rw, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         deallocate(buffer)
      else
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, rx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, ry, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, rz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, rw, nr_points_on_this_proc, irank_master, &
                          global_communicator)
      end if

!     if there is a rest of points, send them to the last proc
      nr_points_rest = mod(nr_points_batch, xc_mpi_get_nr_proc())
      if (nr_points_rest > 0) then
         if (xc_mpi_is_master()) then
            i = nr_points_batch - nr_points_rest + 1
            call interface_mpi_send(rx(i), nr_points_rest, xc_mpi_get_nr_proc()-1, 90, global_communicator)
            call interface_mpi_send(ry(i), nr_points_rest, xc_mpi_get_nr_proc()-1, 91, global_communicator)
            call interface_mpi_send(rz(i), nr_points_rest, xc_mpi_get_nr_proc()-1, 92, global_communicator)
            call interface_mpi_send(rw(i), nr_points_rest, xc_mpi_get_nr_proc()-1, 93, global_communicator)
         end if
         if (get_my_rank() == (xc_mpi_get_nr_proc()-1)) then
            i = nr_points_on_this_proc + 1
            call interface_mpi_recv(rx(i), nr_points_rest, irank_master, 90, global_communicator)
            call interface_mpi_recv(ry(i), nr_points_rest, irank_master, 91, global_communicator)
            call interface_mpi_recv(rz(i), nr_points_rest, irank_master, 92, global_communicator)
            call interface_mpi_recv(rw(i), nr_points_rest, irank_master, 93, global_communicator)
            nr_points_on_this_proc = nr_points_on_this_proc + nr_points_rest
         end if
      end if

      allocate(dist_overview(0:xc_mpi_get_nr_proc()-1))
      dist_overview = 0
      dist_overview(get_my_rank()) = nr_points_on_this_proc
      call xc_mpi_reduce(dist_overview)
      if (xc_mpi_is_master()) then
         do i = 0, xc_mpi_get_nr_proc()-1
!           write(*, *) 'proc', i+1, 'will calculate', dist_overview(i), 'points'
         end do
      end if
      deallocate(dist_overview)

!     verify that sum matches nr_points_batch
      isum = nr_points_on_this_proc
      call xc_mpi_reduce(isum)
      if (xc_mpi_is_master()) then
         if (isum /= nr_points_batch) then
            call quit('programming error in distribution of points in parallel DFT')
         end if
      end if
   end subroutine

   subroutine xc_interface_mpi_barrier()
      call interface_mpi_barrier(global_communicator)
   end subroutine

   subroutine xc_interface_mpi_file_open_to_write(file_name, file_handle, file_info)
      integer(kind=integer_kind)       :: file_handle
      integer(kind=integer_kind)       :: file_info
      integer(kind=integer_kind)       :: amode
      character*(*)                    :: file_name

      amode = df_mpi_mode_create + df_mpi_mode_wronly + df_mpi_mode_sequential

      call interface_mpi_file_open(global_communicator, file_name, amode, file_info, file_handle)

   end subroutine

   subroutine xc_interface_mpi_file_set_view_r8(file_handle, disp)
      integer(kind=integer_kind)       :: file_handle
      integer(kind=integer_kind)       :: info
      integer(kind=integer_kind)       :: etype
      integer(kind=integer_kind)       :: ftype
!      character*(*)                    :: datarep
      integer(kind=df_mpi_offset_kind) :: disp

      etype   = df_mpi_real8
      ftype   = df_mpi_real8
!      datarep = "native"
      info    = df_mpi_info_null

      call interface_mpi_file_set_view(file_handle, disp, etype, ftype, "native", info)

   end subroutine

   subroutine xc_interface_mpi_file_close(file_handle)
      integer(kind=integer_kind)       :: file_handle

      call interface_mpi_file_close(file_handle)

   end subroutine

   subroutine xc_interface_mpi_info_create(file_info)
      integer(kind=integer_kind) :: file_info

      call interface_mpi_info_create(file_info)

   end subroutine

   subroutine xc_interface_mpi_info_free(file_info)
      integer(kind=integer_kind) :: file_info

      call interface_mpi_info_free(file_info)

   end subroutine

   subroutine xc_interface_mpi_file_write_at(file_handle, offset, buf, length)
      integer(kind=integer_kind)       :: file_handle
      real(8)                          :: buf(:)
      integer(kind=df_mpi_offset_kind) :: offset
      integer(kind=integer_kind)       :: length
      integer(kind=integer_kind)       :: status_array(MPI_STATUS_SIZE)

      call interface_mpi_file_write_at(file_handle, offset, buf, length, status_array)

   end subroutine

#endif /* ifdef VAR_MPI */

end module
