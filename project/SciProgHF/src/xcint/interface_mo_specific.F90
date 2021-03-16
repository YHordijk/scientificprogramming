module interface_mo_specific

   use file_units
   use interface_mo

   implicit none

   public interface_mo_write

   private

   interface interface_mo_write
#ifdef PRG_DIRAC
      module procedure interface_mo_write_dirac
#endif
#ifdef VAR_DALTON
      module procedure interface_mo_write_dalton
#endif
   end interface

contains

#ifdef PRG_DIRAC
   subroutine interface_mo_write_dirac()

!     --------------------------------------------------------------------------
      real(8) :: dummy
      integer :: i
      logical :: lucoef_opened_by_dirac
!     --------------------------------------------------------------------------

#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbdhf.h"

      call interface_mo_clear()

      inquire(file = 'DFCOEF', exist = interface_file_exists)
      if (interface_file_exists) then
         inquire(unit=lucoef, opened=lucoef_opened_by_dirac)
         if (lucoef_opened_by_dirac) then
            close(lucoef, status = 'keep')
         end if

         open(interface_file_unit,              &
              file   = 'DFCOEF',      &
              status = 'unknown',     &
              form   = 'unformatted', &
              access = 'sequential')
         rewind(interface_file_unit)
         !miro: increased array size as it's sometimes getting out of bounds: NORBT --> NORBT*2
         allocate(mo_eigenvalues(norbt*2))
         mo_eigenvalues = 0.0d0
         allocate(mo_coef(n2bbasxq))
         mo_coef = 0.0d0

         call reacmo(interface_file_unit,    &
                             'DFCOEF',       &
                             mo_coef,        &
                             mo_eigenvalues, &
                             (/0/),          &
                             dummy,          &
                             6)
         mo_eigenvalues_available = .true.
         mo_coef_available        = .true.

         close(interface_file_unit, status = 'keep')
         if (lucoef_opened_by_dirac) then
            open(lucoef,                 &
                 file   = 'DFCOEF',      &
                 status = 'unknown',     &
                 form   = 'unformatted', &
                 access = 'sequential')
         end if
      end if

!     always regenerate this file
      inquire(file = interface_file_name_mo, exist = interface_file_exists)
      if (interface_file_exists) then
         open(interface_file_unit,             &
              file   = interface_file_name_mo, &
              status = 'unknown',    &
              form   = 'formatted',  &
              access = 'sequential')
         close(interface_file_unit, status = 'delete')
      end if
      open(interface_file_unit,             &
           file   = interface_file_name_mo, &
           status = 'new',        &
           form   = 'formatted',  &
           access = 'sequential')
      rewind(interface_file_unit)

      nr_electrons_inactive = nelect
      nr_electrons_active   = naelec

      nr_mo_gerade_negative_secondary = npsh(1)
      nr_mo_gerade_positive_inactive  = nish(1)
      nr_mo_gerade_positive_active    = nash(1)
      nr_mo_gerade_positive_secondary = nesh(1) - nr_mo_gerade_positive_inactive

      nr_mo_ungerade_negative_secondary = npsh(2)
      nr_mo_ungerade_positive_inactive  = nish(2)
      nr_mo_ungerade_positive_active    = nash(2)
      nr_mo_ungerade_positive_secondary = nesh(2) - nr_mo_ungerade_positive_inactive

      if (aoc) then
         nr_open_shells = nopen
      else
         nr_open_shells = 0
      end if

      write(interface_file_unit, '(4x, a)') 'nr_electrons_inactive_total'
      write(interface_file_unit, '(4x, i6)') nr_electrons_inactive
      write(interface_file_unit, '(4x, a)') 'nr_electrons_active_total'
      write(interface_file_unit, '(4x, i6)') nr_electrons_active
      write(interface_file_unit, '(4x, a)') 'nr_mo_gerade_negative_secondary'
      write(interface_file_unit, '(4x, i6)') nr_mo_gerade_negative_secondary
      write(interface_file_unit, '(4x, a)') 'nr_mo_gerade_positive_inactive'
      write(interface_file_unit, '(4x, i6)') nr_mo_gerade_positive_inactive
      write(interface_file_unit, '(4x, a)') 'nr_mo_gerade_positive_active'
      write(interface_file_unit, '(4x, i6)') nr_mo_gerade_positive_active
      write(interface_file_unit, '(4x, a)') 'nr_mo_gerade_positive_secondary'
      write(interface_file_unit, '(4x, i6)') nr_mo_gerade_positive_secondary
      write(interface_file_unit, '(4x, a)') 'nr_mo_ungerade_negative_secondary'
      write(interface_file_unit, '(4x, i6)') nr_mo_ungerade_negative_secondary
      write(interface_file_unit, '(4x, a)') 'nr_mo_ungerade_positive_inactive'
      write(interface_file_unit, '(4x, i6)') nr_mo_ungerade_positive_inactive
      write(interface_file_unit, '(4x, a)') 'nr_mo_ungerade_positive_active'
      write(interface_file_unit, '(4x, i6)') nr_mo_ungerade_positive_active
      write(interface_file_unit, '(4x, a)') 'nr_mo_ungerade_positive_secondary'
      write(interface_file_unit, '(4x, i6)') nr_mo_ungerade_positive_secondary
      write(interface_file_unit, '(4x, a)') 'nr_open_shells'
      write(interface_file_unit, '(4x, i6)') nr_open_shells
      if (nr_open_shells > 0) then
         do i = 1, nr_open_shells
            write(interface_file_unit, '(4x, i6, e26.16)') i, df(i)
         end do
      end if

      nr_mo_gerade = nr_mo_gerade_negative_secondary &
                   + nr_mo_gerade_positive_inactive  &
                   + nr_mo_gerade_positive_active    &
                   + nr_mo_gerade_positive_secondary
      nr_mo_ungerade = nr_mo_ungerade_negative_secondary &
                     + nr_mo_ungerade_positive_inactive  &
                     + nr_mo_ungerade_positive_active    &
                     + nr_mo_ungerade_positive_secondary
      nr_mo = nr_mo_gerade + nr_mo_ungerade

      if (mo_eigenvalues_available) then
         write(interface_file_unit, '(4x, a)') 'mo_eigenvalues'
         do i = 1, nr_mo
            write(interface_file_unit, '(4x, i6, e26.16)') i, mo_eigenvalues(i)
         end do
      end if

      close(interface_file_unit, status = 'keep')

      if (mo_coef_available) then
!        always regenerate this file
         inquire(file = interface_file_name_mo_coef, exist = interface_file_exists)
         if (interface_file_exists) then
            open(interface_file_unit,                  &
                 file   = interface_file_name_mo_coef, &
                 status = 'unknown',         &
                 form   = 'formatted',       &
                 access = 'sequential')
            close(interface_file_unit, status = 'delete')
         end if
         open(interface_file_unit,                  &
              file   = interface_file_name_mo_coef, &
              status = 'new',             &
              form   = 'formatted',       &
              access = 'sequential')
         rewind(interface_file_unit)
         write(interface_file_unit, *) size(mo_coef)
         do i = 1, size(mo_coef)
            write(interface_file_unit, *) mo_coef(i)
         end do
         close(interface_file_unit, status = 'keep')
      end if

      call interface_mo_clear()

   end subroutine
#endif

#ifdef VAR_DALTON
   subroutine interface_mo_write_dalton()

#include "inforb.h"

      call interface_mo_clear()

!     always regenerate this file
      inquire(file = interface_file_name_mo, exist = interface_file_exists)
      if (interface_file_exists) then
         open(interface_file_unit,             &
              file   = interface_file_name_mo, &
              status = 'unknown',    &
              form   = 'formatted',  &
              access = 'sequential')
         close(interface_file_unit, status = 'delete')
      end if
      open(interface_file_unit,             &
           file   = interface_file_name_mo, &
           status = 'new',        &
           form   = 'formatted',  &
           access = 'sequential')
      rewind(interface_file_unit)

      nr_mo_gerade_negative_secondary = 0
      nr_mo_gerade_positive_inactive  = nisht
      nr_mo_gerade_positive_active    = nasht
      nr_mo_gerade_positive_secondary = nssht

      nr_electrons_inactive = nr_mo_gerade_positive_inactive
      nr_electrons_active   = nr_mo_gerade_positive_active

      nr_mo_ungerade_negative_secondary = 0
      nr_mo_ungerade_positive_inactive  = 0
      nr_mo_ungerade_positive_active    = 0
      nr_mo_ungerade_positive_secondary = 0

      nr_open_shells = 0

      write(interface_file_unit, '(4x, a)') 'nr_electrons_inactive_total'
      write(interface_file_unit, '(4x, i6)') nr_electrons_inactive
      write(interface_file_unit, '(4x, a)') 'nr_electrons_active_total'
      write(interface_file_unit, '(4x, i6)') nr_electrons_active
      write(interface_file_unit, '(4x, a)') 'nr_mo_gerade_negative_secondary'
      write(interface_file_unit, '(4x, i6)') nr_mo_gerade_negative_secondary
      write(interface_file_unit, '(4x, a)') 'nr_mo_gerade_positive_inactive'
      write(interface_file_unit, '(4x, i6)') nr_mo_gerade_positive_inactive
      write(interface_file_unit, '(4x, a)') 'nr_mo_gerade_positive_active'
      write(interface_file_unit, '(4x, i6)') nr_mo_gerade_positive_active
      write(interface_file_unit, '(4x, a)') 'nr_mo_gerade_positive_secondary'
      write(interface_file_unit, '(4x, i6)') nr_mo_gerade_positive_secondary
      write(interface_file_unit, '(4x, a)') 'nr_mo_ungerade_negative_secondary'
      write(interface_file_unit, '(4x, i6)') nr_mo_ungerade_negative_secondary
      write(interface_file_unit, '(4x, a)') 'nr_mo_ungerade_positive_inactive'
      write(interface_file_unit, '(4x, i6)') nr_mo_ungerade_positive_inactive
      write(interface_file_unit, '(4x, a)') 'nr_mo_ungerade_positive_active'
      write(interface_file_unit, '(4x, i6)') nr_mo_ungerade_positive_active
      write(interface_file_unit, '(4x, a)') 'nr_mo_ungerade_positive_secondary'
      write(interface_file_unit, '(4x, i6)') nr_mo_ungerade_positive_secondary
      write(interface_file_unit, '(4x, a)') 'nr_open_shells'
      write(interface_file_unit, '(4x, i6)') nr_open_shells

      call interface_mo_clear()

   end subroutine
#endif

end module
