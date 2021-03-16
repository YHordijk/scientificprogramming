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

module visual_london

  use memory_allocator
  use visual_cfg
  use matrix_defop_old

  implicit none

  public get_d_con_mo
  public get_mdm_con

  private
#include "priunit.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"

  integer :: nr_xyz_comp
  integer :: which_icomp
  logical :: sum_ep_rn

contains

  subroutine get_d_con_mo(D_con_mo, occupation, D_con_mo_ai)
!   ----------------------------------------------------------------------------
    type(matrix)            :: D_con_mo
    real(8)                 :: occupation(*)
    type(matrix), optional  :: D_con_mo_ai

    real(8), allocatable    :: D_ai(:, :, :), D_ji(:, :, :)
    integer                 :: file_unit, icomp
    logical                 :: debug_me, dummy

!   gosia, note for me:
!   if natural connection - T^B matrix was built from 'd|S>/dB' integrals
!   if symmetric connection - T^B matrix was built from 'dS/dB' integrals
!   ----------------------------------------------------------------------------

!   read T^B matrix in MO basis to D_con_mo
!   ---------------------------------------

    file_unit = 71

    open(file_unit,                 &
         file   = 'TBMO',           &
         status = 'old',            &
         form   = 'unformatted',    &
         access = 'sequential',     &
         position = 'rewind',       &
         action = 'read')

1   continue

    read(file_unit) nr_xyz_comp
    read(file_unit) which_icomp
    read(file_unit) D_con_mo%irep
    read(file_unit) sum_ep_rn
    read(file_unit) D_con_mo%elms

    if (nr_xyz_comp > 1) then
!      for example after calculating full shielding tensor in preceding response
!      there will be three matrices on TBMO file: T^Bx, T^By and T^Bz
!      read here the requested one:
       if (visual_cfg_london_component == 'Y') then
          if (which_icomp == 2) then
             go to 2
          else
             go to 1
          end if
       else if (visual_cfg_london_component == 'Z') then
          if (which_icomp == 3) then
             go to 2
          else
             go to 1
          end if
       end if
    end if

2   continue

    close(file_unit, status = 'keep')

    write(lupri, *) 'visual: connection matrix (T^B in MO basis) read from TBMO file'
    write(lupri, *) 'nr_xyz_comp = ', nr_xyz_comp
    write(lupri, *) 'component = ', which_icomp
    write(lupri, *) 'irep = ', D_con_mo%irep
    write(lupri, *) 'sum over e+p shells = ', sum_ep_rn

    debug_me = .false.
    if (debug_me) then
       call print_mat(D_con_mo, label = 'D_con_mo read from TBMO')
    end if

    if (visual_cfg_use_orbital_string) then
       call select_orbitals_reortho(D_con_mo%elms, occupation)
    end if

    if (visual_cfg_london_none) then

!      here discern occ-occ and occ-virt blocks of T^B
!      write occ-occ  to D_con_mo
!      and   occ-virt to D_con_mo_ai
       call alloc(D_ji, norbt, norbt, nz)
       call partition_connection_matrix(D_con_mo%elms, D_con_mo_ai%elms, D_ji, sum_ep_rn)
       D_con_mo%elms = 0.0d0
       call dcopy(norbt*norbt*nz, D_ji, 1, D_con_mo%elms, 1)
       call dealloc(D_ji)

    end if

  end subroutine

   subroutine select_orbitals_reortho(A, occupation)

      real(8), intent(inout)  :: A(norbt, norbt, nz)
      real(8), intent(in)     :: occupation(*)

      real(8), allocatable    :: B(:, :, :)
      integer                 :: i, k, iz

      allocate(B(norbt, norbt, nz))
      B = 0.0d0

      ! radovan: i (re)wrote this but don't quite get it
      !          but it reproduces previous more complicated code
      do iz = 1, nz
         do k = 1, norbt
            do i = 1, norbt
               if (dabs(occupation(k)) > 1.0d-20) then
                  B(i, k, iz) = A(i, k, iz)*occupation(k)
               end if
               if (dabs(occupation(i)) < 1.0d-20) then
                  B(k, i, iz) = A(k, i, iz)*occupation(k)
               end if
            end do
         end do
      end do

      A = B
      deallocate(B)

   end subroutine



  subroutine partition_connection_matrix(T, T_ai, T_ji, sum_ep_rn)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: T(norbt, norbt, nz)
     real(8), intent(inout) :: T_ai(norbt, norbt, nz)
     real(8), intent(inout) :: T_ji(norbt, norbt, nz)
     logical, intent(in)    :: sum_ep_rn

     integer :: ip(2), ap(2), ie(2), ae(2)
     integer :: ifrp, iz, i, j
!    ---------------------------------------------------------------------------

     call dcopy(norbt*norbt*nz, T, 1, T_ji, 1)
     call dcopy(norbt*norbt*nz, T, 1, T_ai, 1)

     do ifrp = 1, nfsym
        ip(ifrp) = iorb(ifrp) + 1
        ap(ifrp) = ip(ifrp) + nish(ifrp)
        ie(ifrp) = ip(ifrp) + npsh(ifrp)
        ae(ifrp) = ie(ifrp) + nish(ifrp)

        write(*, *) 'ifrp, ip = ', ifrp, ip(ifrp)
        write(*, *) 'ifrp, ap = ', ifrp, ap(ifrp)
        write(*, *) 'ifrp, ie = ', ifrp, ie(ifrp)
        write(*, *) 'ifrp, ae = ', ifrp, ae(ifrp)
        write(*, *) 'ifrp, nish = ', ifrp, nish(ifrp)
        write(*, *) 'ifrp, npsh = ', ifrp, npsh(ifrp)
        write(*, *) 'ifrp, norb = ', ifrp, norb(ifrp)
     end do

! todo:
!    if .not. sum_ep_rn - reorthonormalization only over electronic shells
!    (so i can zero-out blocks corresponding to positronic shells)

!    'ai' blocks (zero-out 'ij' blocks)
     do iz = 1, nz
        do ifrp = 1, nfsym

!          zero-out 'ij' blocks
           do j = ie(ifrp), ae(ifrp) - 1
              do i = ie(ifrp), ae(ifrp) - 1
                 T_ai(i, j, iz) = 0.0d0
              end do
           end do

!          zero-out 'aa' blocks
           do j = ip(ifrp), ie(ifrp) - 1
              do i = ip(ifrp), ie(ifrp) - 1
                 T_ai(i, j, iz) = 0.0d0
              end do
              do i = ae(ifrp), norb(ifrp)
                 T_ai(i, j, iz) = 0.0d0
              end do
           end do
           do j = ae(ifrp), norb(ifrp)
              do i = ip(ifrp), ie(ifrp) - 1
                 T_ai(i, j, iz) = 0.0d0
              end do
              do i = ae(ifrp), norb(ifrp)
                 T_ai(i, j, iz) = 0.0d0
              end do
           end do

        end do
     end do

!    'ij' blocks
     do iz = 1, nz
        do ifrp = 1, nfsym

           do j = ip(ifrp), norb(ifrp)
              do i = ip(ifrp), ie(ifrp) - 1
                 T_ji(i, j, iz) = 0.0d0
              end do
              do i = ae(ifrp), norb(ifrp)
                 T_ji(i, j, iz) = 0.0d0
              end do
           end do
           do j = ie(ifrp), ae(ifrp) - 1
              do i = ip(ifrp), ie(ifrp) - 1
                 T_ji(i, j, iz) = 0.0d0
              end do
              do i = ae(ifrp), norb(ifrp)
                 T_ji(i, j, iz) = 0.0d0
              end do
           end do

        end do
     end do

  end subroutine


  subroutine get_mdm_con(D_con,       &
                         D_con_mo,    &
                         C,           &
                         Cig, Csg,    &
                         Ciu, Csu)

!   ----------------------------------------------------------------------------
    type(matrix)            :: D_con, D_con_mo
    type(matrix)            :: C, Cig, Csg
    type(matrix), optional  :: Ciu, Csu

!   ----------------------------------------------------------------------------

    if (visual_cfg_london_none) then

!      here build density matrix modified by first-order overlap
!      (c^+_\mu j)*(S^B)_ij*(c_\nu i) where (S^B)_ij = (T^B)_ij + (T^B+)_ij

       D_con = D_con_mo - dag(D_con_mo)
!      check that it gives the same regardless of connection used
!      due to: T^B + (T^B)^+ = -S^B
!      call print_mat(D_con, label = '-S^B')

       if (nfsym == 2) then
          if (jbtof(D_con_mo%irep, 1) == 2) then
!            ungerade perturbation
             D_con = Cig*(D_con*dag(Ciu)) &
                   + Ciu*(D_con*dag(Cig))
          else
!            gerade perturbation
             D_con = Cig*(D_con*dag(Cig)) &
                   + Ciu*(D_con*dag(Ciu))
          end if
       else
          D_con = Cig*(D_con*dag(Cig))
       end if

    else
!      symmetric or natural connection used in preceding response run
!      "-" is because dag() does only transposition and T^B is antihermitian

       if (nfsym == 2) then
          D_con = C*(D_con_mo*dag(Cig)) &
                + C*(D_con_mo*dag(Ciu))
          D_con = D_con - dag(D_con)
       else

          D_con = C*(D_con_mo*dag(Cig))
          D_con = D_con - dag(D_con)

       end if


    end if

  end subroutine

end module
