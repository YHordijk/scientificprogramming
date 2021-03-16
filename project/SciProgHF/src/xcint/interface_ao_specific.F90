module interface_ao_specific

   use interface_ao
   use file_units

   implicit none

   public interface_ao_write

   private

   interface interface_ao_write
#ifdef PRG_DIRAC
      module procedure interface_ao_write_dirac
#endif
#ifdef VAR_DALTON
      module procedure interface_ao_write_dalton
#endif
   end interface

contains

#ifdef PRG_DIRAC
   subroutine interface_ao_write_dirac()

!     --------------------------------------------------------------------------
      integer                 :: i, j, ic
      integer                 :: ishell, icount, iprim, irep, icent
      integer, allocatable    :: nr_primitives(:)
      character(1), parameter :: ls(2) = (/'L', 'S'/)
!     --------------------------------------------------------------------------

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "shells.h"
#include "symmet.h"
#include "aovec.h"
#include "primit.h"
#include "sphtrm.h"
#include "dcbbas.h"
#include "pincom.h"
#include "dgroup.h"
#include "dcbham.h"
#include "nuclei.h"
#include "orgcom.h"

      call interface_ao_clear()

!     always regenerate this file
      inquire(file = interface_file_name_ao, exist = interface_file_exists)
      if (interface_file_exists) then
         open(interface_file_unit,          &
              file   = interface_file_name_ao, &
              status = 'unknown',              &
              form   = 'formatted',            &
              access = 'sequential')
         close(interface_file_unit, status = 'delete')
      end if
      open(interface_file_unit,          &
           file   = interface_file_name_ao, &
           status = 'new',                  &
           form   = 'formatted',            &
           access = 'sequential')
      rewind(interface_file_unit)

      nr_boson_ireps    = nbsym
      nr_fermion_coreps = nfsym
      nr_shells         = kmax

      allocate(l_quantum_nr(nr_shells))
      allocate(shell_center(nr_shells))
      do ishell = 1, nr_shells
         l_quantum_nr(ishell) = nhkt(ishell) - 1
         shell_center(ishell) = ncent(ishell)
      end do

      allocate(nr_primitives(nr_shells))
      do ishell = 1, nr_shells
         i = jstrt(ishell) + 1
         j = jstrt(ishell) + nuco(ishell)
         icount = 0
         do iprim = i, j
            if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
               icount = icount + 1
               nr_primitives(ishell) = icount
            end if
         end do
      end do

      allocate(primitive_exp(nr_shells,    sum(nr_primitives)))
      allocate(contraction_coef(nr_shells, sum(nr_primitives)))
      do ishell = 1, nr_shells
         i = jstrt(ishell) + 1
         j = jstrt(ishell) + nuco(ishell)
         icount = 0
         do iprim = i, j
            if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
               icount = icount + 1
               contraction_coef(ishell, icount) = priccf(iprim, numcf(ishell))
               primitive_exp(ishell, icount)    = priexp(iprim)
            end if
         end do
      end do

      write(interface_file_unit, '(a)') '*** geometry'

      write(interface_file_unit, '(4x, a)') 'nr_centers'
      write(interface_file_unit, '(4x, i6)') nucind

      write(interface_file_unit, '(4x, a)') 'charges_and_coordinates'
      do icent = 1, nucind
         write(interface_file_unit, '(4f27.16)') charge(icent),  &
                                                    cord(1, icent), &
                                                    cord(2, icent), &
                                                    cord(3, icent)
      end do

      write(interface_file_unit, '(a)') '*** basis'

      basis_is_spherical = .false.
      do ishell = 1, kmax
         if (sphr(ishell)) then
            basis_is_spherical = .true.
         end if
      end do
      write(interface_file_unit, '(4x, a)') 'is_spherical'
      write(interface_file_unit, '(4x, l6)') basis_is_spherical

      write(interface_file_unit, '(4x, a)') 'algebra'
      write(interface_file_unit, '(4x, i6)') nz

      use_only_large = .false.
      if (NOSMLV) then
         use_only_large = .true.
      end if
      write(interface_file_unit, '(4x, a)') 'use_only_large'
      write(interface_file_unit, '(4x, l6)') use_only_large

      write(interface_file_unit, '(4x, a)') 'nr_primitive_exp'
      write(interface_file_unit, '(4x, i6)') sum(nr_primitives)

      write(interface_file_unit, '(4x, 4a6, 2a26)') '# cent', &
                                                       'L/S',    &
                                                       'shell',  &
                                                       'l',      &
                                                       'prim',   &
                                                       'coef'
      do ishell = 1, kmax
         if (ishell > nlrgsh) then
            ic = 2
         else
            ic = 1
         end if
         do i = 1, nr_primitives(ishell)
            write(interface_file_unit, '(4x, i6, a6, 2i6, 2e26.16)') shell_center(ishell),     &
                                                                        ls(ic),                   &
                                                                        ishell,                   &
                                                                        l_quantum_nr(ishell),     &
                                                                        primitive_exp(ishell, i), &
                                                                        contraction_coef(ishell, i)
         end do
      end do
      deallocate(nr_primitives)

      !gosia: also write to interface_file_unit?
      !do i = 1, 3
      !     global_gauge_origin(i) = gagorg(i)
      !end do

      if (nr_boson_ireps > 1) then
         write(interface_file_unit, '(a)') '*** symmetry'

         write(interface_file_unit, '(4x, a)') 'center_degeneracy_and_stabilizing_irep'
         do icent = 1, nucind
            write(interface_file_unit, '(4x, 3i6)') icent, nucdeg(icent), istbnu(icent)
         end do

         write(interface_file_unit, '(4x, a)') 'nr_boson_ireps'
         write(interface_file_unit, '(4x, i6)') nr_boson_ireps

         write(interface_file_unit, '(4x, a)') 'nr_fermion_coreps'
         write(interface_file_unit, '(4x, i6)') nr_fermion_coreps

         write(interface_file_unit, '(4x, a)') 'irep_of_axes'
         write(interface_file_unit, '(4x, i6)') isymax(1, 1)
         write(interface_file_unit, '(4x, i6)') isymax(2, 1)
         write(interface_file_unit, '(4x, i6)') isymax(3, 1)

         write(interface_file_unit, '(4x, a)') 'irep_of_rotations'
         write(interface_file_unit, '(4x, i6)') isymax(1, 2)
         write(interface_file_unit, '(4x, i6)') isymax(2, 2)
         write(interface_file_unit, '(4x, i6)') isymax(3, 2)

         write(interface_file_unit, '(4x, a6, 10a8)') '# irep',  &
                                                         'nr-L',    &
                                                         'nr-S',    &
                                                         'D-off-L', &
                                                         'D-off-S', &
                                                         'jqbas-L', &
                                                         'jqbas-S', &
                                                         'jbtof-L', &
                                                         'jbtof-S', &
                                                         'H-off-L', &
                                                         'H-off-S'
         do irep = 0, nr_boson_ireps - 1
            write(interface_file_unit, '(4x, i6, 10i8)') irep,            &
                                                            nbbas(irep, 1),  &
                                                            nbbas(irep, 2),  &
                                                            ibbas(irep, 1),  &
                                                            ibbas(irep, 2),  &
                                                            jqbas(irep, 1),  &
                                                            jqbas(irep, 2),  &
                                                            jbtof(irep, 1),  &
                                                            jbtof(irep, 2),  &
                                                            icos(irep+1, 1), &
                                                            icos(irep+1, 2)
         end do

         pq_to_uq = ipqtoq
         uq_to_pq = iqtopq
         write(interface_file_unit, '(4x, a)') 'quaternion_packing'
         do i = 1, 4
            write(interface_file_unit, '(8x, 8i2)') pq_to_uq(i, 0:7)
         end do
         do i = 1, 4
            write(interface_file_unit, '(8x, 8i2)') uq_to_pq(i, 0:7)
         end do

         iq(0) = 1
         do i = 1, 3
            iq(i) = jm4pos(4 - i)
            iq(i) = iqtopq(iq(i), isymax(i, 2))
         end do
         write(interface_file_unit, '(4x, a)') 'quaternion_units'
         write(interface_file_unit, '(8x, 4i2)') iq(0:3)
      end if

      close(interface_file_unit, status = 'keep')

      call interface_ao_clear()

   end subroutine
#endif

#ifdef VAR_DALTON
   subroutine interface_ao_write_dalton()

!     --------------------------------------------------------------------------
      integer                 :: i, j, ic
      integer                 :: ishell, icount, iprim, irep, icent
      integer, allocatable    :: nr_primitives(:)
      character(1), parameter :: ls(2) = (/'L', 'S'/)
!     --------------------------------------------------------------------------

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "shells.h"
#include "symmet.h"
#include "aovec.h"
#include "primit.h"
#include "nuclei.h"

      call interface_ao_clear()

!     always regenerate this file
      inquire(file = interface_file_name_ao, exist = interface_file_exists)
      if (interface_file_exists) then
         open(interface_file_unit,          &
              file   = interface_file_name_ao, &
              status = 'unknown',              &
              form   = 'formatted',            &
              access = 'sequential')
         close(interface_file_unit, status = 'delete')
      end if
      open(interface_file_unit,          &
           file   = interface_file_name_ao, &
           status = 'new',                  &
           form   = 'formatted',            &
           access = 'sequential')
      rewind(interface_file_unit)

      nr_boson_ireps    = maxrep + 1
!     nr_fermion_coreps = nfsym
      nr_shells         = kmax

      allocate(l_quantum_nr(nr_shells))
      allocate(shell_center(nr_shells))
      do ishell = 1, nr_shells
         l_quantum_nr(ishell) = nhkt(ishell) - 1
         shell_center(ishell) = ncent(ishell)
      end do

      allocate(nr_primitives(nr_shells))
      do ishell = 1, nr_shells
         i = jstrt(ishell) + 1
         j = jstrt(ishell) + nuco(ishell)
         icount = 0
         do iprim = i, j
            if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
               icount = icount + 1
               nr_primitives(ishell) = icount
            end if
         end do
      end do

      allocate(primitive_exp(nr_shells,    sum(nr_primitives)))
      allocate(contraction_coef(nr_shells, sum(nr_primitives)))
      do ishell = 1, nr_shells
         i = jstrt(ishell) + 1
         j = jstrt(ishell) + nuco(ishell)
         icount = 0
         do iprim = i, j
            if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
               icount = icount + 1
               contraction_coef(ishell, icount) = priccf(iprim, numcf(ishell))
               primitive_exp(ishell, icount)    = priexp(iprim)
            end if
         end do
      end do

      write(interface_file_unit, '(a)') '*** geometry'

      write(interface_file_unit, '(4x, a)') 'nr_centers'
      write(interface_file_unit, '(4x, i6)') nucind

      write(interface_file_unit, '(4x, a)') 'charges_and_coordinates'
      do icent = 1, nucind
         write(interface_file_unit, '(4f27.16)') charge(icent),  &
                                                    cord(1, icent), &
                                                    cord(2, icent), &
                                                    cord(3, icent)
      end do

      write(interface_file_unit, '(a)') '*** basis'

      basis_is_spherical = .false.
      do ishell = 1, kmax
         if (sphr(ishell)) then
            basis_is_spherical = .true.
         end if
      end do
      write(interface_file_unit, '(4x, a)') 'is_spherical'
      write(interface_file_unit, '(4x, l6)') basis_is_spherical

      write(interface_file_unit, '(4x, a)') 'algebra'
      write(interface_file_unit, '(4x, i6)') 1

      write(interface_file_unit, '(4x, a)') 'use_only_large'
      write(interface_file_unit, '(4x, l6)') .true.

      write(interface_file_unit, '(4x, a)') 'nr_primitive_exp'
      write(interface_file_unit, '(4x, i6)') sum(nr_primitives)

      write(interface_file_unit, '(4x, 4a6, 2a26)') '# cent', &
                                                       'L/S',    &
                                                       'shell',  &
                                                       'l',      &
                                                       'prim',   &
                                                       'coef'
      do ishell = 1, kmax
         if (ishell > nlrgsh) then
            ic = 2
         else
            ic = 1
         end if
         do i = 1, nr_primitives(ishell)
            write(interface_file_unit, '(4x, i6, a6, 2i6, 2e26.16)') shell_center(ishell),     &
                                                                        ls(ic),                   &
                                                                        ishell,                   &
                                                                        l_quantum_nr(ishell),     &
                                                                        primitive_exp(ishell, i), &
                                                                        contraction_coef(ishell, i)
         end do
      end do
      deallocate(nr_primitives)

      if (nr_boson_ireps > 1) then
         stop 'program me'
      end if

      close(interface_file_unit, status = 'keep')

      call interface_ao_clear()

   end subroutine
#endif

end module
