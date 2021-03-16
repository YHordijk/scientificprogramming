module pcm_gp

   implicit none

   public collect_nctot
   public collect_atoms
   public collect_symmetry_info
   public timing

   private

   contains

      function collect_nctot() result(nr_nuclei)

      use, intrinsic :: iso_c_binding, only: c_int

#include "mxcent.h"
#include "nuclei.h"

      integer(c_int) :: nr_nuclei

      nr_nuclei = nucdep

      end function collect_nctot

      subroutine collect_atoms(atomic_charges, atomic_centers)

      use, intrinsic :: iso_c_binding, only: c_double

#include "mxcent.h"
#include "nuclei.h"

      real(c_double), intent(out) :: atomic_charges(*)
      real(c_double), intent(out) :: atomic_centers(3,*)

      integer :: i, j, k

! Get coordinates
      call getacord(atomic_centers)
! Get charges
      i = 0
      do j = 1, nucind
        do k = 1, nucdeg(j)
          i = i + 1
          atomic_charges(i) = charge(j)
        enddo
      enddo

      end subroutine collect_atoms

      subroutine collect_symmetry_info(symmetry_info)

      use, intrinsic :: iso_c_binding, only: c_int

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"

      integer(c_int), intent(inout) :: symmetry_info(4)

      symmetry_info(1) = 0_c_int !pcm_igen(1)
      symmetry_info(2) = 0_c_int !pcm_igen(2)
      symmetry_info(3) = 0_c_int !pcm_igen(3)
      symmetry_info(4) = 0_c_int !pcm_igen(4)

      end subroutine collect_symmetry_info

      real(8) function timing()

      use interface_to_mpi

#if defined (VAR_MPI)
      timing = interface_MPI_WTIME()
#else
      call cpu_time(timing)
#endif

      end function timing

end module pcm_gp
