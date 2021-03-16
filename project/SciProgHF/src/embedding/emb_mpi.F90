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

module fde_mpi

#ifdef VAR_MPI

   use interface_to_mpi
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

   public fde_mpi_is_master
   public fde_mpi_get_nr_proc
#ifdef VAR_MPI
   public fde_mpi_bcast
   public fde_mpi_reduce
   public fde_mpi_distribute_points
   public fde_mpi_distribute_frozen
   public fde_mpi_distribute_frozen_pertden_n
   public fde_mpi_distribute_frozen_pertden_gn
   public fde_mpi_distribute_frozen_pertden_s
   public fde_mpi_distribute_frozen_pertden_gs
   public fde_mpi_distribute_vemb
   public fde_get_my_rank
   public fde_mpi_gather
   public fde_mpi_gatherv
#endif

   private

#ifdef VAR_MPI
   interface fde_get_my_rank
      module procedure get_my_rank
   end interface
   interface fde_mpi_bcast
      module procedure fde_mpi_bcast_i0
      module procedure fde_mpi_bcast_i1
      module procedure fde_mpi_bcast_i2
      module procedure fde_mpi_bcast_i3
      module procedure fde_mpi_bcast_i4
      module procedure fde_mpi_bcast_r0
      module procedure fde_mpi_bcast_r1
      module procedure fde_mpi_bcast_r2
      module procedure fde_mpi_bcast_r3
      module procedure fde_mpi_bcast_r4
      module procedure fde_mpi_bcast_l0
      module procedure fde_mpi_bcast_l1
   end interface

   interface fde_mpi_reduce
      module procedure fde_mpi_reduce_i0
      module procedure fde_mpi_reduce_i1
      module procedure fde_mpi_reduce_r0
      module procedure fde_mpi_reduce_r1
      module procedure fde_mpi_reduce_r2
      module procedure fde_mpi_reduce_r4
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

   function fde_mpi_is_master()
      logical :: fde_mpi_is_master
#include "dcbgen.h"
      fde_mpi_is_master = .true.
#ifdef VAR_MPI
      if (parcal) then
         fde_mpi_is_master = (get_my_rank() == irank_master)
      end if
#endif
   end function

   function fde_mpi_get_nr_proc()
      integer :: fde_mpi_get_nr_proc
      integer :: nr_proc
#include "dcbgen.h"
      fde_mpi_get_nr_proc = 1
#ifdef VAR_MPI
      if (parcal) then
         call interface_mpi_comm_size(global_communicator, nr_proc)
         fde_mpi_get_nr_proc = nr_proc
      end if
#endif
   end function

#ifdef VAR_MPI

    subroutine fde_mpi_gatherv(x, y, n)
       real(8) :: x(:), y(:)
       integer, intent(in) :: n
       integer :: i
       integer, allocatable :: rcounts(:)
       integer, allocatable :: displs(:)
       integer, allocatable :: rcounts_node(:)
       integer :: nproc

! set ncouts and displs arrays
!
       nproc = fde_mpi_get_nr_proc()

       allocate(rcounts_node(nproc))
       allocate(rcounts(nproc))
       allocate(displs(nproc))

       rcounts_node = 0
       rcounts = 0
       rcounts_node(1) = n

       call mpi_gather(rcounts_node, 1, mpi_integer, rcounts, 1, mpi_integer, irank_master, global_communicator, ierr)

       displs(1) = 0
       do i = 2, nproc
             displs(i) = displs(i-1) + rcounts(i-1)
       end do

       call mpi_gatherv(x, n, mpi_double_precision, y, rcounts, displs, mpi_double_precision, &
    &                  irank_master, global_communicator, ierr)

       deallocate(rcounts_node)
       deallocate(rcounts)
       deallocate(displs)

    end subroutine


    subroutine fde_mpi_gather(x,y)
       real(8) :: x(:), y(:)
       integer :: n
       integer :: i

       n = size(x)/fde_mpi_get_nr_proc() 

       call mpi_gather(x, n, mpi_double_precision, y, n, mpi_double_precision, irank_master, global_communicator, ierr)

    end subroutine

   subroutine fde_mpi_bcast_i0(x)
      integer :: x
      call interface_mpi_bcast(x, 1, irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_i1(x)
      integer :: x(:)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_i2(x)
      integer :: x(:, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_i3(x)
      integer :: x(:, :, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_i4(x)
      integer :: x(:, :, :, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_r0(x)
      real(8) :: x
      call interface_mpi_bcast(x, 1,  irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_r1(x)
      real(8) :: x(:)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_r2(x)
      real(8) :: x(:, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_r3(x)
      real(8) :: x(:, :, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_r4(x)
      real(8) :: x(:, :, :, :)
      call interface_mpi_bcast(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_l0(x)
      logical :: x
      call interface_mpi_bcast_l0(x, 1, irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_bcast_l1(x)
      logical :: x(:)
      call interface_mpi_bcast_l1(x, size(x),  irank_master, global_communicator)
   end subroutine

   subroutine fde_mpi_reduce_i0(x)
      integer :: x
      if (fde_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, 1, mpi_integer, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, 1, mpi_integer, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine fde_mpi_reduce_i1(x)
      integer :: x(:)
      if (fde_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_integer, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_integer, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine fde_mpi_reduce_r0(x)
      real(8) :: x
      if (fde_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, 1, mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, 1, mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine fde_mpi_reduce_r1(x)
      real(8) :: x(:)
      if (fde_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine fde_mpi_reduce_r2(x)
      real(8) :: x(:, :)
      if (fde_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine fde_mpi_reduce_r4(x)
      real(8) :: x(:, :, :, :)
      if (fde_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_double_precision, mpi_sum, 0, global_communicator, ierr)
      end if
   end subroutine

   subroutine fde_mpi_scatter_r1(x, n)
      real(8)              :: x(:)
      real(8), allocatable :: buffer(:)
      integer              :: n
   end subroutine

   subroutine fde_mpi_distribute_points(rx, ry, rz, rw, nr_points_batch, nr_points_on_this_proc)
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

      nr_points_on_this_proc = nr_points_batch/fde_mpi_get_nr_proc()
      if (fde_mpi_is_master()) then
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
      nr_points_rest = mod(nr_points_batch, fde_mpi_get_nr_proc())
      if (nr_points_rest > 0) then
         if (fde_mpi_is_master()) then
            i = nr_points_batch - nr_points_rest + 1
            call interface_mpi_send(rx(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 90, global_communicator)
            call interface_mpi_send(ry(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 91, global_communicator)
            call interface_mpi_send(rz(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 92, global_communicator)
            call interface_mpi_send(rw(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 93, global_communicator)
         end if
         if (get_my_rank() == (fde_mpi_get_nr_proc()-1)) then
            i = nr_points_on_this_proc + 1
            call interface_mpi_recv(rx(i), nr_points_rest, irank_master, 90, global_communicator)
            call interface_mpi_recv(ry(i), nr_points_rest, irank_master, 91, global_communicator)
            call interface_mpi_recv(rz(i), nr_points_rest, irank_master, 92, global_communicator)
            call interface_mpi_recv(rw(i), nr_points_rest, irank_master, 93, global_communicator) 
            nr_points_on_this_proc = nr_points_on_this_proc + nr_points_rest
         end if
      end if

      allocate(dist_overview(0:fde_mpi_get_nr_proc()-1))
      dist_overview = 0
      dist_overview(get_my_rank()) = nr_points_on_this_proc
      call fde_mpi_reduce(dist_overview)
      if (fde_mpi_is_master()) then
         do i = 0, fde_mpi_get_nr_proc()-1
!           write(*, *) 'proc', i+1, 'will calculate', dist_overview(i),
!           'points'
         end do
      end if
      deallocate(dist_overview)

!     verify that sum matches nr_points_batch
      isum = nr_points_on_this_proc
      call fde_mpi_reduce(isum)
      if (fde_mpi_is_master()) then
         if (isum /= nr_points_batch) then
            call quit('programming error in distribution of points in parallel FDE')
         end if
      end if
   end subroutine


   subroutine fde_mpi_distribute_frozen(n, gnx, gny, gnz, elpot, nr_points_batch, nr_points_on_this_proc)
      real(8)              :: n(*)
      real(8)              :: gnx(*)
      real(8)              :: gny(*)
      real(8)              :: gnz(*)
      real(8)              :: elpot(*)
      integer              :: nr_points_batch
      integer              :: nr_points_on_this_proc
      real(8), allocatable :: buffer(:)
      integer              :: nr_points_rest
      integer              :: i
      integer, allocatable :: dist_overview(:)
      integer              :: isum

      nr_points_on_this_proc = nr_points_batch/fde_mpi_get_nr_proc()
      if (fde_mpi_is_master()) then
         allocate(buffer(nr_points_on_this_proc))
         call interface_mpi_scatter_r1_work_f77(n,   nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gnx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gny, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gnz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(elpot, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         deallocate(buffer)
      else
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, n,   nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gnx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gny, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gnz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, elpot, nr_points_on_this_proc, irank_master, &
                          global_communicator)
      end if

!     if there is a rest of points, send them to the last proc
      nr_points_rest = mod(nr_points_batch, fde_mpi_get_nr_proc())
      if (nr_points_rest > 0) then
         if (fde_mpi_is_master()) then
            i = nr_points_batch - nr_points_rest + 1
            call interface_mpi_send(n(i),   nr_points_rest, fde_mpi_get_nr_proc()-1, 93, global_communicator)
            call interface_mpi_send(elpot(i),nr_points_rest, fde_mpi_get_nr_proc()-1,94, global_communicator)
            call interface_mpi_send(gnx(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 90, global_communicator)
            call interface_mpi_send(gny(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 91, global_communicator)
            call interface_mpi_send(gnz(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 92, global_communicator)
         end if
         if (get_my_rank() == (fde_mpi_get_nr_proc()-1)) then
            i = nr_points_on_this_proc + 1
            call interface_mpi_recv(n(i),   nr_points_rest, irank_master,  93, global_communicator)
            Call interface_mpi_recv(elpot(i), nr_points_rest, irank_master,94, global_communicator)
            call interface_mpi_recv(gnx(i), nr_points_rest, irank_master,  90, global_communicator)
            call interface_mpi_recv(gny(i), nr_points_rest, irank_master,  91, global_communicator)
            call interface_mpi_recv(gnz(i), nr_points_rest, irank_master,  92, global_communicator)
            nr_points_on_this_proc = nr_points_on_this_proc + nr_points_rest
         end if
      end if

      allocate(dist_overview(0:fde_mpi_get_nr_proc()-1))
      dist_overview = 0
      dist_overview(get_my_rank()) = nr_points_on_this_proc
      call fde_mpi_reduce(dist_overview)
      if (fde_mpi_is_master()) then
         do i = 0, fde_mpi_get_nr_proc()-1
!           write(*, *) 'proc', i+1, 'will calculate', dist_overview(i), 'points'
         end do
      end if
      deallocate(dist_overview)

!     verify that sum matches nr_points_batch
      isum = nr_points_on_this_proc
      call fde_mpi_reduce(isum)
      if (fde_mpi_is_master()) then
         if (isum /= nr_points_batch) then
            call quit('programming error in distribution of frozen in parallel FDE')
         end if
      end if
   end subroutine


   subroutine fde_mpi_distribute_frozen_pertden_n(n_bx, n_by, n_bz, nr_points_batch, nr_points_on_this_proc)
      real(8)              :: n_bx(*)
      real(8)              :: n_by(*)
      real(8)              :: n_bz(*)
      integer              :: nr_points_batch
      integer              :: nr_points_on_this_proc
      real(8), allocatable :: buffer(:)
      integer              :: nr_points_rest
      integer              :: i
      integer, allocatable :: dist_overview(:)
      integer              :: isum

      nr_points_on_this_proc = nr_points_batch/fde_mpi_get_nr_proc()
      if (fde_mpi_is_master()) then
         allocate(buffer(nr_points_on_this_proc))
         call interface_mpi_scatter_r1_work_f77(n_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(n_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(n_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         deallocate(buffer)
      else
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, n_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, n_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, n_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
      end if

!     if there is a rest of points, send them to the last proc
      nr_points_rest = mod(nr_points_batch, fde_mpi_get_nr_proc())
      if (nr_points_rest > 0) then
         if (fde_mpi_is_master()) then
            i = nr_points_batch - nr_points_rest + 1
            call interface_mpi_send(n_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 95, global_communicator)
            call interface_mpi_send(n_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 96, global_communicator)
            call interface_mpi_send(n_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 97, global_communicator)
         end if
         if (get_my_rank() == (fde_mpi_get_nr_proc()-1)) then
            i = nr_points_on_this_proc + 1
            call interface_mpi_recv(n_bx(i), nr_points_rest, irank_master,  95, global_communicator)
            call interface_mpi_recv(n_by(i), nr_points_rest, irank_master,  96, global_communicator)
            call interface_mpi_recv(n_bz(i), nr_points_rest, irank_master,  97, global_communicator)
            nr_points_on_this_proc = nr_points_on_this_proc + nr_points_rest
         end if
      end if

      allocate(dist_overview(0:fde_mpi_get_nr_proc()-1))
      dist_overview = 0
      dist_overview(get_my_rank()) = nr_points_on_this_proc
      call fde_mpi_reduce(dist_overview)
      if (fde_mpi_is_master()) then
         do i = 0, fde_mpi_get_nr_proc()-1
!           write(*, *) 'proc', i+1, 'will calculate', dist_overview(i), 'points'
         end do
      end if
      deallocate(dist_overview)

!     verify that sum matches nr_points_batch
      isum = nr_points_on_this_proc
      call fde_mpi_reduce(isum)
      if (fde_mpi_is_master()) then
         if (isum /= nr_points_batch) then
            call quit('programming error in distribution of pertden n in parallel FDE')
         end if
      end if
   end subroutine

   subroutine fde_mpi_distribute_frozen_pertden_gn(gxn_bx, gxn_by, gxn_bz, &
                                                   gyn_bx, gyn_by, gyn_bz, &
                                                   gzn_bx, gzn_by, gzn_bz, &
                                                   nr_points_batch, nr_points_on_this_proc)
      real(8)              :: gxn_bx(*)
      real(8)              :: gxn_by(*)
      real(8)              :: gxn_bz(*)
      real(8)              :: gyn_bx(*)
      real(8)              :: gyn_by(*)
      real(8)              :: gyn_bz(*)
      real(8)              :: gzn_bx(*)
      real(8)              :: gzn_by(*)
      real(8)              :: gzn_bz(*)
      integer              :: nr_points_batch
      integer              :: nr_points_on_this_proc
      real(8), allocatable :: buffer(:)
      integer              :: nr_points_rest
      integer              :: i
      integer, allocatable :: dist_overview(:)
      integer              :: isum

      nr_points_on_this_proc = nr_points_batch/fde_mpi_get_nr_proc()
      if (fde_mpi_is_master()) then
         allocate(buffer(nr_points_on_this_proc))
         call interface_mpi_scatter_r1_work_f77(gxn_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gxn_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gxn_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gyn_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gyn_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gyn_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzn_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzn_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzn_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         deallocate(buffer)
      else
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxn_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxn_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxn_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gyn_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gyn_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gyn_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzn_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzn_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzn_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
      end if

!     if there is a rest of points, send them to the last proc
      nr_points_rest = mod(nr_points_batch, fde_mpi_get_nr_proc())
      if (nr_points_rest > 0) then
         if (fde_mpi_is_master()) then
            i = nr_points_batch - nr_points_rest + 1
            call interface_mpi_send(gxn_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,101, global_communicator)
            call interface_mpi_send(gxn_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,102, global_communicator)
            call interface_mpi_send(gxn_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,103, global_communicator)
            call interface_mpi_send(gyn_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,104, global_communicator)
            call interface_mpi_send(gyn_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,105, global_communicator)
            call interface_mpi_send(gyn_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,106, global_communicator)
            call interface_mpi_send(gzn_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,107, global_communicator)
            call interface_mpi_send(gzn_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,108, global_communicator)
            call interface_mpi_send(gzn_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,109, global_communicator)
         end if
         if (get_my_rank() == (fde_mpi_get_nr_proc()-1)) then
            i = nr_points_on_this_proc + 1
            call interface_mpi_recv(gxn_bx(i), nr_points_rest, irank_master, 101, global_communicator)
            call interface_mpi_recv(gxn_by(i), nr_points_rest, irank_master, 102, global_communicator)
            call interface_mpi_recv(gxn_bz(i), nr_points_rest, irank_master, 103, global_communicator)
            call interface_mpi_recv(gyn_bx(i), nr_points_rest, irank_master, 104, global_communicator)
            call interface_mpi_recv(gyn_by(i), nr_points_rest, irank_master, 105, global_communicator)
            call interface_mpi_recv(gyn_bz(i), nr_points_rest, irank_master, 106, global_communicator)
            call interface_mpi_recv(gzn_bx(i), nr_points_rest, irank_master, 107, global_communicator)
            call interface_mpi_recv(gzn_by(i), nr_points_rest, irank_master, 108, global_communicator)
            call interface_mpi_recv(gzn_bz(i), nr_points_rest, irank_master, 109, global_communicator)
            nr_points_on_this_proc = nr_points_on_this_proc + nr_points_rest
         end if
      end if

      allocate(dist_overview(0:fde_mpi_get_nr_proc()-1))
      dist_overview = 0
      dist_overview(get_my_rank()) = nr_points_on_this_proc
      call fde_mpi_reduce(dist_overview)
      if (fde_mpi_is_master()) then
         do i = 0, fde_mpi_get_nr_proc()-1
!           write(*, *) 'proc', i+1, 'will calculate', dist_overview(i), 'points'
         end do
      end if
      deallocate(dist_overview)

!     verify that sum matches nr_points_batch
      isum = nr_points_on_this_proc
      call fde_mpi_reduce(isum)
      if (fde_mpi_is_master()) then
         if (isum /= nr_points_batch) then
            call quit('programming error in distribution of pertden gn in parallel FDE')
         end if
      end if
   end subroutine

   subroutine fde_mpi_distribute_frozen_pertden_s(sx_bx, sx_by, sx_bz, &
                                                  sy_bx, sy_by, sy_bz, &
                                                  sz_bx, sz_by, sz_bz, &
                                                  nr_points_batch, nr_points_on_this_proc)
      real(8)              :: sx_bx(*)
      real(8)              :: sx_by(*)
      real(8)              :: sx_bz(*)
      real(8)              :: sy_bx(*)
      real(8)              :: sy_by(*)
      real(8)              :: sy_bz(*)
      real(8)              :: sz_bx(*)
      real(8)              :: sz_by(*)
      real(8)              :: sz_bz(*)
      integer              :: nr_points_batch
      integer              :: nr_points_on_this_proc
      real(8), allocatable :: buffer(:)
      integer              :: nr_points_rest
      integer              :: i
      integer, allocatable :: dist_overview(:)
      integer              :: isum

      nr_points_on_this_proc = nr_points_batch/fde_mpi_get_nr_proc()
      if (fde_mpi_is_master()) then
         allocate(buffer(nr_points_on_this_proc))
         call interface_mpi_scatter_r1_work_f77(sx_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(sx_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(sx_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(sy_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(sy_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(sy_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(sz_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(sz_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(sz_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         deallocate(buffer)
      else
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, sx_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, sx_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, sx_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, sy_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, sy_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, sy_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, sz_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, sz_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, sz_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
      end if

!     if there is a rest of points, send them to the last proc
      nr_points_rest = mod(nr_points_batch, fde_mpi_get_nr_proc())
      if (nr_points_rest > 0) then
         if (fde_mpi_is_master()) then
            i = nr_points_batch - nr_points_rest + 1
            call interface_mpi_send(sx_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,110, global_communicator)
            call interface_mpi_send(sx_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,111, global_communicator)
            call interface_mpi_send(sx_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,112, global_communicator)
            call interface_mpi_send(sy_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,113, global_communicator)
            call interface_mpi_send(sy_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,114, global_communicator)
            call interface_mpi_send(sy_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,115, global_communicator)
            call interface_mpi_send(sz_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,116, global_communicator)
            call interface_mpi_send(sz_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,117, global_communicator)
            call interface_mpi_send(sz_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,118, global_communicator)
         end if
         if (get_my_rank() == (fde_mpi_get_nr_proc()-1)) then
            i = nr_points_on_this_proc + 1
            call interface_mpi_recv(sx_bx(i), nr_points_rest, irank_master, 110, global_communicator)
            call interface_mpi_recv(sx_by(i), nr_points_rest, irank_master, 111, global_communicator)
            call interface_mpi_recv(sx_bz(i), nr_points_rest, irank_master, 112, global_communicator)
            call interface_mpi_recv(sy_bx(i), nr_points_rest, irank_master, 113, global_communicator)
            call interface_mpi_recv(sy_by(i), nr_points_rest, irank_master, 114, global_communicator)
            call interface_mpi_recv(sy_bz(i), nr_points_rest, irank_master, 115, global_communicator)
            call interface_mpi_recv(sz_bx(i), nr_points_rest, irank_master, 116, global_communicator)
            call interface_mpi_recv(sz_by(i), nr_points_rest, irank_master, 117, global_communicator)
            call interface_mpi_recv(sz_bz(i), nr_points_rest, irank_master, 118, global_communicator)
            nr_points_on_this_proc = nr_points_on_this_proc + nr_points_rest
         end if
      end if

      allocate(dist_overview(0:fde_mpi_get_nr_proc()-1))
      dist_overview = 0
      dist_overview(get_my_rank()) = nr_points_on_this_proc
      call fde_mpi_reduce(dist_overview)
      if (fde_mpi_is_master()) then
         do i = 0, fde_mpi_get_nr_proc()-1
!           write(*, *) 'proc', i+1, 'will calculate', dist_overview(i), 'points'
         end do
      end if
      deallocate(dist_overview)

!     verify that sum matches nr_points_batch
      isum = nr_points_on_this_proc
      call fde_mpi_reduce(isum)
      if (fde_mpi_is_master()) then
         if (isum /= nr_points_batch) then
            call quit('programming error in distribution of pertden s in parallel FDE')
         end if
      end if
   end subroutine

   subroutine fde_mpi_distribute_frozen_pertden_gs(gxsx_bx, gxsx_by, gxsx_bz, &
                                                   gxsy_bx, gxsy_by, gxsy_bz, &
                                                   gxsz_bx, gxsz_by, gxsz_bz, &
                                                   gysx_bx, gysx_by, gysx_bz, &
                                                   gysy_bx, gysy_by, gysy_bz, &
                                                   gysz_bx, gysz_by, gysz_bz, &
                                                   gzsx_bx, gzsx_by, gzsx_bz, &
                                                   gzsy_bx, gzsy_by, gzsy_bz, &
                                                   gzsz_bx, gzsz_by, gzsz_bz, &
                                                   nr_points_batch, nr_points_on_this_proc)
      real(8)              :: gxsx_bx(*)
      real(8)              :: gxsx_by(*)
      real(8)              :: gxsx_bz(*)
      real(8)              :: gxsy_bx(*)
      real(8)              :: gxsy_by(*)
      real(8)              :: gxsy_bz(*)
      real(8)              :: gxsz_bx(*)
      real(8)              :: gxsz_by(*)
      real(8)              :: gxsz_bz(*)
      real(8)              :: gysx_bx(*)
      real(8)              :: gysx_by(*)
      real(8)              :: gysx_bz(*)
      real(8)              :: gysy_bx(*)
      real(8)              :: gysy_by(*)
      real(8)              :: gysy_bz(*)
      real(8)              :: gysz_bx(*)
      real(8)              :: gysz_by(*)
      real(8)              :: gysz_bz(*)
      real(8)              :: gzsx_bx(*)
      real(8)              :: gzsx_by(*)
      real(8)              :: gzsx_bz(*)
      real(8)              :: gzsy_bx(*)
      real(8)              :: gzsy_by(*)
      real(8)              :: gzsy_bz(*)
      real(8)              :: gzsz_bx(*)
      real(8)              :: gzsz_by(*)
      real(8)              :: gzsz_bz(*)
      integer              :: nr_points_batch
      integer              :: nr_points_on_this_proc
      real(8), allocatable :: buffer(:)
      integer              :: nr_points_rest
      integer              :: i
      integer, allocatable :: dist_overview(:)
      integer              :: isum

      nr_points_on_this_proc = nr_points_batch/fde_mpi_get_nr_proc()
      if (fde_mpi_is_master()) then
         allocate(buffer(nr_points_on_this_proc))
         call interface_mpi_scatter_r1_work_f77(gxsx_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gxsx_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gxsx_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gxsy_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gxsy_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gxsy_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gxsz_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gxsz_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gxsz_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gysx_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gysx_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gysx_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gysy_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gysy_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gysy_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gysz_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gysz_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gysz_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzsx_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzsx_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzsx_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzsy_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzsy_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzsy_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzsz_bx, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzsz_by, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         call interface_mpi_scatter_r1_work_f77(gzsz_bz, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         deallocate(buffer)
      else
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxsx_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxsx_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxsx_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxsy_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxsy_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxsy_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxsz_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxsz_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gxsz_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gysx_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gysx_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gysx_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gysy_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gysy_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gysy_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gysz_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gysz_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gysz_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzsx_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzsx_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzsx_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzsy_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzsy_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzsy_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzsz_bx, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzsz_by, nr_points_on_this_proc, irank_master, &
                          global_communicator)
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, gzsz_bz, nr_points_on_this_proc, irank_master, &
                          global_communicator)
      end if

!     if there is a rest of points, send them to the last proc
      nr_points_rest = mod(nr_points_batch, fde_mpi_get_nr_proc())
      if (nr_points_rest > 0) then
         if (fde_mpi_is_master()) then
            i = nr_points_batch - nr_points_rest + 1
            call interface_mpi_send(gxsx_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,120, global_communicator)
            call interface_mpi_send(gxsx_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,121, global_communicator)
            call interface_mpi_send(gxsx_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,122, global_communicator)
            call interface_mpi_send(gxsy_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,123, global_communicator)
            call interface_mpi_send(gxsy_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,124, global_communicator)
            call interface_mpi_send(gxsy_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,125, global_communicator)
            call interface_mpi_send(gxsz_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,126, global_communicator)
            call interface_mpi_send(gxsz_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,127, global_communicator)
            call interface_mpi_send(gxsz_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,128, global_communicator)
            call interface_mpi_send(gysx_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,130, global_communicator)
            call interface_mpi_send(gysx_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,131, global_communicator)
            call interface_mpi_send(gysx_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,132, global_communicator)
            call interface_mpi_send(gysy_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,133, global_communicator)
            call interface_mpi_send(gysy_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,134, global_communicator)
            call interface_mpi_send(gysy_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,135, global_communicator)
            call interface_mpi_send(gysz_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,136, global_communicator)
            call interface_mpi_send(gysz_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,137, global_communicator)
            call interface_mpi_send(gysz_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,138, global_communicator)
            call interface_mpi_send(gzsx_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,140, global_communicator)
            call interface_mpi_send(gzsx_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,141, global_communicator)
            call interface_mpi_send(gzsx_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,142, global_communicator)
            call interface_mpi_send(gzsy_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,143, global_communicator)
            call interface_mpi_send(gzsy_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,144, global_communicator)
            call interface_mpi_send(gzsy_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,145, global_communicator)
            call interface_mpi_send(gzsz_bx(i), nr_points_rest, fde_mpi_get_nr_proc()-1,146, global_communicator)
            call interface_mpi_send(gzsz_by(i), nr_points_rest, fde_mpi_get_nr_proc()-1,147, global_communicator)
            call interface_mpi_send(gzsz_bz(i), nr_points_rest, fde_mpi_get_nr_proc()-1,148, global_communicator)
         end if
         if (get_my_rank() == (fde_mpi_get_nr_proc()-1)) then
            i = nr_points_on_this_proc + 1
            call interface_mpi_recv(gxsx_bx(i), nr_points_rest, irank_master, 120, global_communicator)
            call interface_mpi_recv(gxsx_by(i), nr_points_rest, irank_master, 121, global_communicator)
            call interface_mpi_recv(gxsx_bz(i), nr_points_rest, irank_master, 122, global_communicator)
            call interface_mpi_recv(gxsy_bx(i), nr_points_rest, irank_master, 123, global_communicator)
            call interface_mpi_recv(gxsy_by(i), nr_points_rest, irank_master, 124, global_communicator)
            call interface_mpi_recv(gxsy_bz(i), nr_points_rest, irank_master, 125, global_communicator)
            call interface_mpi_recv(gxsz_bx(i), nr_points_rest, irank_master, 126, global_communicator)
            call interface_mpi_recv(gxsz_by(i), nr_points_rest, irank_master, 127, global_communicator)
            call interface_mpi_recv(gxsz_bz(i), nr_points_rest, irank_master, 128, global_communicator)
            call interface_mpi_recv(gysx_bx(i), nr_points_rest, irank_master, 130, global_communicator)
            call interface_mpi_recv(gysx_by(i), nr_points_rest, irank_master, 131, global_communicator)
            call interface_mpi_recv(gysx_bz(i), nr_points_rest, irank_master, 132, global_communicator)
            call interface_mpi_recv(gysy_bx(i), nr_points_rest, irank_master, 133, global_communicator)
            call interface_mpi_recv(gysy_by(i), nr_points_rest, irank_master, 134, global_communicator)
            call interface_mpi_recv(gysy_bz(i), nr_points_rest, irank_master, 135, global_communicator)
            call interface_mpi_recv(gysz_bx(i), nr_points_rest, irank_master, 136, global_communicator)
            call interface_mpi_recv(gysz_by(i), nr_points_rest, irank_master, 137, global_communicator)
            call interface_mpi_recv(gysz_bz(i), nr_points_rest, irank_master, 138, global_communicator)
            call interface_mpi_recv(gzsx_bx(i), nr_points_rest, irank_master, 140, global_communicator)
            call interface_mpi_recv(gzsx_by(i), nr_points_rest, irank_master, 141, global_communicator)
            call interface_mpi_recv(gzsx_bz(i), nr_points_rest, irank_master, 142, global_communicator)
            call interface_mpi_recv(gzsy_bx(i), nr_points_rest, irank_master, 143, global_communicator)
            call interface_mpi_recv(gzsy_by(i), nr_points_rest, irank_master, 144, global_communicator)
            call interface_mpi_recv(gzsy_bz(i), nr_points_rest, irank_master, 145, global_communicator)
            call interface_mpi_recv(gzsz_bx(i), nr_points_rest, irank_master, 146, global_communicator)
            call interface_mpi_recv(gzsz_by(i), nr_points_rest, irank_master, 147, global_communicator)
            call interface_mpi_recv(gzsz_bz(i), nr_points_rest, irank_master, 148, global_communicator)
            nr_points_on_this_proc = nr_points_on_this_proc + nr_points_rest
         end if
      end if

      allocate(dist_overview(0:fde_mpi_get_nr_proc()-1))
      dist_overview = 0
      dist_overview(get_my_rank()) = nr_points_on_this_proc
      call fde_mpi_reduce(dist_overview)
      if (fde_mpi_is_master()) then
         do i = 0, fde_mpi_get_nr_proc()-1
!           write(*, *) 'proc', i+1, 'will calculate', dist_overview(i), 'points'
         end do
      end if
      deallocate(dist_overview)

!     verify that sum matches nr_points_batch
      isum = nr_points_on_this_proc
      call fde_mpi_reduce(isum)
      if (fde_mpi_is_master()) then
         if (isum /= nr_points_batch) then
            call quit('programming error in distribution of pertden gs in parallel FDE')
         end if
      end if
   end subroutine


   subroutine fde_mpi_distribute_vemb(vemb, nr_points_batch, nr_points_on_this_proc)
      real(8)              :: vemb(*)
      integer              :: nr_points_batch
      integer              :: nr_points_on_this_proc
      real(8), allocatable :: buffer(:)
      integer              :: nr_points_rest
      integer              :: i
      integer, allocatable :: dist_overview(:)
      integer              :: isum

      nr_points_on_this_proc = nr_points_batch/fde_mpi_get_nr_proc()
      if (fde_mpi_is_master()) then
         allocate(buffer(nr_points_on_this_proc))
         call interface_mpi_scatter_r1_work_f77(vemb, nr_points_on_this_proc, buffer, &
              nr_points_on_this_proc, irank_master, global_communicator)
         deallocate(buffer)
      else
         call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, vemb, nr_points_on_this_proc, irank_master, &
                          global_communicator)
      end if

!     if there is a rest of points, send them to the last proc
      nr_points_rest = mod(nr_points_batch, fde_mpi_get_nr_proc())
      if (nr_points_rest > 0) then
         if (fde_mpi_is_master()) then
            i = nr_points_batch - nr_points_rest + 1
            call interface_mpi_send(vemb(i), nr_points_rest, fde_mpi_get_nr_proc()-1, 90, global_communicator)
         end if
         if (get_my_rank() == (fde_mpi_get_nr_proc()-1)) then
            i = nr_points_on_this_proc + 1
            call interface_mpi_recv(vemb(i), nr_points_rest, irank_master,  90, global_communicator)
            nr_points_on_this_proc = nr_points_on_this_proc + nr_points_rest
         end if
      end if

      allocate(dist_overview(0:fde_mpi_get_nr_proc()-1))
      dist_overview = 0
      dist_overview(get_my_rank()) = nr_points_on_this_proc
      call fde_mpi_reduce(dist_overview)
      if (fde_mpi_is_master()) then
         do i = 0, fde_mpi_get_nr_proc()-1
!           write(*, *) 'proc', i+1, 'will calculate', dist_overview(i), 'points'
         end do
      end if
      deallocate(dist_overview)

!     verify that sum matches nr_points_batch
      isum = nr_points_on_this_proc
      call fde_mpi_reduce(isum)
      if (fde_mpi_is_master()) then
         if (isum /= nr_points_batch) then
            call quit('programming error in distribution of vemb in parallel FDE')
         end if
      end if
   end subroutine

#endif /* ifdef VAR_MPI */

end module
