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

module num_grid_gen

   implicit none

   public reset_num_grid
   public generate_num_grid
   public stash_zipgrid
   public unstash_zipgrid

   private

   logical :: grid_is_done     = .false.
   logical :: remember_zipgrid = .false.

contains

   subroutine reset_num_grid()

      grid_is_done = .false.

   end subroutine

   subroutine stash_zipgrid(zipgrid)

      logical, intent(inout) :: zipgrid

      remember_zipgrid = zipgrid
      zipgrid          = .false.

   end subroutine

   subroutine unstash_zipgrid(zipgrid)

      logical, intent(out) :: zipgrid

      zipgrid = remember_zipgrid

   end subroutine

   subroutine generate_num_grid(dmat)

      use dirac_cfg
!     --------------------------------------------------------------------------
      real(8), intent(in) :: dmat(*)
!     --------------------------------------------------------------------------

      if (.not. grid_is_done) then
!       calculate grid
        call dftgrd(dmat, 1)
        grid_is_done = .true.
      end if

   end subroutine

end module
