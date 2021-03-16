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

module dirac_ao_eval

!   this is based on DIRAC and Dalton
!   routines written by T. Helgaker, P. Salek, T. Saue, and others

! - in order to calculate aos using this module
!   you need to initialize some data first
!   by calling dirac_ao_eval_init from each node
! - dirac_ao_eval_init will set symmetry and basis set information
!   that is needed in each point and then mpi_bcast this info to the nodes
! - we try to not recalculate data which does not change from point to point
!   (like symmetry)
! - *after* the module data is initialized you can allocate
!   the ao vector and a buffer, both with the length nr_ao_slices*nr_ao_cartesian
! - (de)allocate them only once (each iteration), not in every point
!   this is also the reason why the ao vector and buffer are not allocated here - this would
!   be very inefficient
! - aos are not dzeroed for more performance
!   if you change the code make sure that elements of these vectors
!   are not read before they are written for the first time
!   note that in levy-leblond or x2c4 mode some elements will
!   contain garbage but they are not (= should not be) used
! - this module is called dft_ao but should be usable also outside dft
!   this is why some information is passed through dirac_ao_eval_init and not queried
!   from dft_cfg, to keep this module generic, please dont "use dft_cfg" here

!   FAQ
!   ===
!
!   1. how is the ao array sorted on output of get_ao?
!      - it is sorted according to ibbas and nbbas
!        a convenient way to see the sorting of (symmetry adapted) aos
!        in the output is to search for 'Generating Lowdin matrix'
!        as an example, with a d2h calculation i got:
!
!                                 Generating Lowdin matrix:
!                                 -------------------------
!
!        L   Ag    * Deleted:          2(Proj:          2, Lindep:          0)
!        L   B1g   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B2g   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B3g   * Deleted:          0(Proj:          0, Lindep:          0)
!        S   B3u   * Deleted:          2(Proj:          2, Lindep:          0)
!        S   B2u   * Deleted:          2(Proj:          2, Lindep:          0)
!        S   B1u   * Deleted:          6(Proj:          6, Lindep:          0)
!        S   Au    * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B3u   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B2u   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B1u   * Deleted:          1(Proj:          1, Lindep:          0)
!        L   Au    * Deleted:          0(Proj:          0, Lindep:          0)
!        S   Ag    * Deleted:          9(Proj:          9, Lindep:          0)
!        S   B1g   * Deleted:          0(Proj:          0, Lindep:          0)
!        S   B2g   * Deleted:          1(Proj:          1, Lindep:          0)
!        S   B3g   * Deleted:          1(Proj:          1, Lindep:          0)
!
!        so the order is:
!        | L Ag | L B1g | L B2g | L B3g | S B3u | ... | S B3g | and then the derivatives


   use interface_ao
   use xc_max_block_length
   use extra
   use dirac_ao_eval_sub

   implicit none

   public get_ao
   public dirac_ao_eval_init

   integer, public :: ao_off_g1_m0(1:3, 0:0)
   integer, public :: ao_off_g0_m1(0:0, 1:3)
   integer, public :: ao_off_g0_m2(0:0, 1:6)
   integer, public :: ao_off_g2_m0(1:6, 0:0)
   integer, public :: ao_off_g1_m1(1:3, 1:3)
   integer, public :: ao_off_g1_m2(1:3, 1:6)

   private

   real(8), parameter :: bitstring_parity(0:7) = (/1.0d0, -1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0/)

   logical, public    :: need_ao_order(0:2, 0:2) = .false.
!  standard LDA calculation: nr_ao_slices = 1
!                            ao vector composition = | ao_g0 |
!                            meaning there are no geometric derivatives
!  standard GGA calculation: nr_ao_slices = 4
!                            ao vector composition = | ao_g0 | ao_gx | ao_gy | ao_gz |
!  with LAO there are more slices
!  also with higher geometric derivatives (soon) there are more and more slices
!  (higher derivatives and mixed derivatives)
!  and the total length of the ao vector is nr_ao_slices*nr_ao_cartesian
!  buffer has the same length for convenience (one could get away with less for the buffer)
   integer, public    :: nr_ao_slices = 0

contains

   subroutine get_ao(block_length, &
                     cx,           &
                     cy,           &
                     cz,           &
                     ao,           &
                     buffer)

!     --------------------------------------------------------------------------
      integer, intent(in)    :: block_length
      real(8), intent(in)    :: cx(*)
      real(8), intent(in)    :: cy(*)
      real(8), intent(in)    :: cz(*)
      real(8)                :: ao(block_length, *)
      real(8), intent(inout) :: buffer(block_length, *)
!     --------------------------------------------------------------------------
      integer                :: i, k
      integer                :: ishell, isym, ixyz
      real(8)                :: px(max_block_length)
      real(8)                :: py(max_block_length)
      real(8)                :: pz(max_block_length)
      real(8)                :: p2(max_block_length)
      real(8)                :: sx
      real(8)                :: sy
      real(8)                :: sz
!     --------------------------------------------------------------------------

      i = 1
      do ishell = 1, nr_shells
         do isym = 0, nr_boson_ireps - 1
            if (iand(isym, shell_stabilizing_irep(ishell)) == 0) then

               sx = bitstring_parity(iand(irep_of_axes(1), isym))*center_xyz(shell_center(ishell), 1)
               sy = bitstring_parity(iand(irep_of_axes(2), isym))*center_xyz(shell_center(ishell), 2)
               sz = bitstring_parity(iand(irep_of_axes(3), isym))*center_xyz(shell_center(ishell), 3)

!#define DEBUG_GAUGEO
#ifdef DEBUG_GAUGEO
               sx = bitstring_parity(iand(irep_of_axes(1), isym))*(center_xyz(shell_center(ishell), 1) - global_gauge_origin(1))
               sy = bitstring_parity(iand(irep_of_axes(2), isym))*(center_xyz(shell_center(ishell), 2) - global_gauge_origin(2))
               sz = bitstring_parity(iand(irep_of_axes(3), isym))*(center_xyz(shell_center(ishell), 3) - global_gauge_origin(3))
#endif

               do k = 1, block_length
                  px(k) = cx(k) - sx
                  py(k) = cy(k) - sy
                  pz(k) = cz(k) - sz
               end do

               do k = 1, block_length
                  p2(k) = px(k)*px(k) &
                        + py(k)*py(k) &
                        + pz(k)*pz(k)
               end do

               if (need_ao_order(2, 0)) then
                  call get_cartesian_ao_g2(block_length,                  &
                                           ishell,                        &
                                           px,                            &
                                           py,                            &
                                           pz,                            &
                                           p2,                            &
                                           ao(1,                      i), &
                                           ao(1, ao_off_g1_m0(1, 0) + i), &
                                           ao(1, ao_off_g1_m0(2, 0) + i), &
                                           ao(1, ao_off_g1_m0(3, 0) + i), &
                                           ao(1, ao_off_g2_m0(1, 0) + i), &
                                           ao(1, ao_off_g2_m0(2, 0) + i), &
                                           ao(1, ao_off_g2_m0(3, 0) + i), &
                                           ao(1, ao_off_g2_m0(4, 0) + i), &
                                           ao(1, ao_off_g2_m0(5, 0) + i), &
                                           ao(1, ao_off_g2_m0(6, 0) + i))
               else
                  if (need_ao_order(1, 0)) then
                     call get_cartesian_ao_g1(block_length,                  &
                                              ishell,                        &
                                              px,                            &
                                              py,                            &
                                              pz,                            &
                                              p2,                            &
                                              ao(1,                      i), &
                                              ao(1, ao_off_g1_m0(1, 0) + i), &
                                              ao(1, ao_off_g1_m0(2, 0) + i), &
                                              ao(1, ao_off_g1_m0(3, 0) + i))
                  else
                     call get_cartesian_ao_g0(block_length, &
                                              ishell,       &
                                              px,           &
                                              py,           &
                                              pz,           &
                                              p2,           &
                                              ao(1, i))
                  end if
               end if

               if (need_ao_order(0, 2)) then
                  call get_cartesian_ao_m2(block_length,                  &
                                           ishell,                        &
                                           px,                            &
                                           py,                            &
                                           pz,                            &
                                           ao(1,                      i), &
                                           ao(1, ao_off_g0_m1(0, 1) + i), &
                                           ao(1, ao_off_g0_m1(0, 2) + i), &
                                           ao(1, ao_off_g0_m1(0, 3) + i), &
                                           ao(1, ao_off_g0_m2(0, 1) + i), &
                                           ao(1, ao_off_g0_m2(0, 2) + i), &
                                           ao(1, ao_off_g0_m2(0, 3) + i), &
                                           ao(1, ao_off_g0_m2(0, 4) + i), &
                                           ao(1, ao_off_g0_m2(0, 5) + i), &
                                           ao(1, ao_off_g0_m2(0, 6) + i))
                  if (need_ao_order(1, 2)) then
                     do ixyz = 1, 3
                           call get_cartesian_ao_m2(block_length,                     &
                                                    ishell,                           &
                                                    px,                               &
                                                    py,                               &
                                                    pz,                               &
                                                    ao(1, ao_off_g1_m0(ixyz, 0) + i), &
                                                    ao(1, ao_off_g1_m1(ixyz, 1) + i), &
                                                    ao(1, ao_off_g1_m1(ixyz, 2) + i), &
                                                    ao(1, ao_off_g1_m1(ixyz, 3) + i), &
                                                    ao(1, ao_off_g1_m2(ixyz, 1) + i), &
                                                    ao(1, ao_off_g1_m2(ixyz, 2) + i), &
                                                    ao(1, ao_off_g1_m2(ixyz, 3) + i), &
                                                    ao(1, ao_off_g1_m2(ixyz, 4) + i), &
                                                    ao(1, ao_off_g1_m2(ixyz, 5) + i), &
                                                    ao(1, ao_off_g1_m2(ixyz, 6) + i))
                     end do
                  end if
               else
                  if (need_ao_order(0, 1)) then
                     call get_cartesian_ao_m1(block_length,                  &
                                              ishell,                        &
                                              px,                            &
                                              py,                            &
                                              pz,                            &
                                              ao(1,                      i), &
                                              ao(1, ao_off_g0_m1(0, 1) + i), &
                                              ao(1, ao_off_g0_m1(0, 2) + i), &
                                              ao(1, ao_off_g0_m1(0, 3) + i))
                     if (need_ao_order(1, 1)) then
                        do ixyz = 1, 3
                           call get_cartesian_ao_m1(block_length,                     &
                                                    ishell,                           &
                                                    px,                               &
                                                    py,                               &
                                                    pz,                               &
                                                    ao(1, ao_off_g1_m0(ixyz, 0) + i), &
                                                    ao(1, ao_off_g1_m1(ixyz, 1) + i), &
                                                    ao(1, ao_off_g1_m1(ixyz, 2) + i), &
                                                    ao(1, ao_off_g1_m1(ixyz, 3) + i))
                        end do
                     end if
                  end if
               end if

               i = i + cartesian_deg(ishell)
            end if
         end do
      end do

      if ((nr_boson_ireps > 1) .or. (nr_ao_cartesian /= nr_ao)) then
!        transform from AO (Hermit sorting)
!        to SA-AO (DIRAC sorting)
!        and/or transform from cartesian to spherical gaussian basis
         call ao_to_so(block_length, ao, buffer)
      end if

   end subroutine

   subroutine get_cartesian_ao_m1(block_length, &
                                  ishell,       &
                                  px,           &
                                  py,           &
                                  pz,           &
                                  ao_m000,      &
                                  ao_m100,      &
                                  ao_m010,      &
                                  ao_m001)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: block_length
      integer, intent(in)  :: ishell
      real(8), intent(in)  :: px(*)
      real(8), intent(in)  :: py(*)
      real(8), intent(in)  :: pz(*)
      real(8), intent(in)  :: ao_m000(block_length, *)
      real(8), intent(out) :: ao_m100(block_length, *)
      real(8), intent(out) :: ao_m010(block_length, *)
      real(8), intent(out) :: ao_m001(block_length, *)
!     --------------------------------------------------------------------------
      integer              :: k
      integer              :: ic
      real(8)              :: f(max_block_length)
!     --------------------------------------------------------------------------

      do ic = 1, cartesian_deg(ishell)
         do k = 1, block_length
            f(k) = 0.5d0*ao_m000(k, ic)
         end do
         do k = 1, block_length
            ao_m100(k, ic) = f(k)*px(k)
            ao_m010(k, ic) = f(k)*py(k)
            ao_m001(k, ic) = f(k)*pz(k)
         end do
      end do

   end subroutine

   subroutine get_cartesian_ao_m2(block_length, &
                                  ishell,       &
                                  px,           &
                                  py,           &
                                  pz,           &
                                  ao_m000,      &
                                  ao_m100,      &
                                  ao_m010,      &
                                  ao_m001,      &
                                  ao_m200,      &
                                  ao_m110,      &
                                  ao_m101,      &
                                  ao_m020,      &
                                  ao_m011,      &
                                  ao_m002)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: block_length
      integer, intent(in)  :: ishell
      real(8), intent(in)  :: px(*)
      real(8), intent(in)  :: py(*)
      real(8), intent(in)  :: pz(*)
      real(8), intent(in)  :: ao_m000(block_length, *)
      real(8), intent(out) :: ao_m100(block_length, *)
      real(8), intent(out) :: ao_m010(block_length, *)
      real(8), intent(out) :: ao_m001(block_length, *)
      real(8), intent(out) :: ao_m200(block_length, *)
      real(8), intent(out) :: ao_m110(block_length, *)
      real(8), intent(out) :: ao_m101(block_length, *)
      real(8), intent(out) :: ao_m020(block_length, *)
      real(8), intent(out) :: ao_m011(block_length, *)
      real(8), intent(out) :: ao_m002(block_length, *)
!     --------------------------------------------------------------------------
      integer              :: k
      integer              :: ic
      real(8)              :: f(max_block_length)
!     --------------------------------------------------------------------------

      do ic = 1, cartesian_deg(ishell)
         do k = 1, block_length
            f(k) = 0.5d0*ao_m000(k, ic)
         end do
         do k = 1, block_length
            ao_m100(k, ic) = f(k)*px(k)
            ao_m010(k, ic) = f(k)*py(k)
            ao_m001(k, ic) = f(k)*pz(k)

            ao_m200(k, ic) = 0.5d0*f(k)*px(k)*px(k)
            ao_m110(k, ic) = 0.5d0*f(k)*px(k)*py(k)
            ao_m101(k, ic) = 0.5d0*f(k)*px(k)*pz(k)
            ao_m020(k, ic) = 0.5d0*f(k)*py(k)*py(k)
            ao_m011(k, ic) = 0.5d0*f(k)*py(k)*pz(k)
            ao_m002(k, ic) = 0.5d0*f(k)*pz(k)*pz(k)
         end do
      end do

   end subroutine


   subroutine ao_to_so(block_length, &
                       ao,           &
                       buffer)

!     --------------------------------------------------------------------------
      integer, intent(in)    :: block_length
      real(8), intent(inout) :: ao(block_length, nr_ao_slices*nr_ao_cartesian)
      real(8), intent(inout) :: buffer(block_length, nr_ao_slices*nr_ao_cartesian)
!     --------------------------------------------------------------------------
      integer                :: i, k, n, off
!     --------------------------------------------------------------------------

      do i = 1, nr_ao_slices*nr_ao_cartesian
         do k = 1, block_length
            buffer(k, i) = 0.0d0
         end do
      end do
      do n = 1, nr_ao_slices
         off = nr_ao_cartesian*(n - 1)
         do i = 1, nr_so_ao
            do k = 1, block_length
               buffer(k, off + i_so_ao(1, i)) = buffer(k, off + i_so_ao(1, i)) &
                                                  + ao(k, off + i_so_ao(2, i))*f_so_ao(i)
            end do
         end do
      end do
      do i = 1, nr_ao_slices*nr_ao_cartesian
         do k = 1, block_length
            ao(k, i) = buffer(k, i)
         end do
      end do

   end subroutine

   subroutine dirac_ao_eval_init(g_order, &
                           m_order, &
                           use_mpi_if_possible)

!     --------------------------------------------------------------------------
      integer, intent(in) :: g_order
      integer, intent(in) :: m_order
      logical, intent(in) :: use_mpi_if_possible
!     --------------------------------------------------------------------------

      call interface_ao_read(use_mpi_if_possible)

      need_ao_order = .false.
      nr_ao_slices  = 1

      if (g_order > 0) then
         ao_off_g1_m0(1, 0) = (0 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m0(2, 0) = (1 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m0(3, 0) = (2 + nr_ao_slices)*nr_ao_cartesian
         nr_ao_slices = nr_ao_slices + 3
         need_ao_order(1, 0) = .true.
      end if

      if (m_order > 0) then
         ao_off_g0_m1(0, 1) = (0 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g0_m1(0, 2) = (1 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g0_m1(0, 3) = (2 + nr_ao_slices)*nr_ao_cartesian
         nr_ao_slices = nr_ao_slices + 3
         need_ao_order(0, 1) = .true.
      end if

      if (g_order > 1) then
         ao_off_g2_m0(1, 0) = (0 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g2_m0(2, 0) = (1 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g2_m0(3, 0) = (2 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g2_m0(4, 0) = (3 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g2_m0(5, 0) = (4 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g2_m0(6, 0) = (5 + nr_ao_slices)*nr_ao_cartesian
         nr_ao_slices = nr_ao_slices + 6
         need_ao_order(2, 0) = .true.
      end if

      if (m_order > 1) then
         ao_off_g0_m2(0, 1) = (0 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g0_m2(0, 2) = (1 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g0_m2(0, 3) = (2 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g0_m2(0, 4) = (3 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g0_m2(0, 5) = (4 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g0_m2(0, 6) = (5 + nr_ao_slices)*nr_ao_cartesian
         nr_ao_slices = nr_ao_slices + 6
         need_ao_order(0, 2) = .true.
      end if

      if (g_order > 0 .and. m_order > 0) then
         ao_off_g1_m1(1, 1) = (0 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m1(1, 2) = (1 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m1(1, 3) = (2 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m1(2, 1) = (3 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m1(2, 2) = (4 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m1(2, 3) = (5 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m1(3, 1) = (6 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m1(3, 2) = (7 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m1(3, 3) = (8 + nr_ao_slices)*nr_ao_cartesian
         nr_ao_slices = nr_ao_slices + 9
         need_ao_order(1, 1) = .true.
      end if

      if (g_order > 0 .and. m_order > 1) then
         ao_off_g1_m2(1, 1) = (0 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(1, 2) = (1 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(1, 3) = (2 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(1, 4) = (3 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(1, 5) = (4 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(1, 6) = (5 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(2, 1) = (6 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(2, 2) = (7 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(2, 3) = (8 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(2, 4) = (9 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(2, 5) = (10 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(2, 6) = (11 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(3, 1) = (12 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(3, 2) = (13 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(3, 3) = (14 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(3, 4) = (15 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(3, 5) = (16 + nr_ao_slices)*nr_ao_cartesian
         ao_off_g1_m2(3, 6) = (17 + nr_ao_slices)*nr_ao_cartesian
         nr_ao_slices = nr_ao_slices + 18
         need_ao_order(1, 2) = .true.
      end if


   end subroutine

end module
