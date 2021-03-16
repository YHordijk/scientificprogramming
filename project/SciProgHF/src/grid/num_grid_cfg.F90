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

module num_grid_cfg

   implicit none

   public report_num_grid

   character(30), public :: num_grid_cfg_gridfile
   integer,       public :: num_grid_cfg_angint                  = 41
   integer,       public :: num_grid_cfg_angmin                  = 15
   integer,       public :: num_grid_cfg_integration_check_level = 0
   logical,       public :: num_grid_cfg_estimate_radii          = .false.
   logical,       public :: num_grid_cfg_force_4c_grid           = .false.
   logical,       public :: num_grid_cfg_import_grid             = .false.
   logical,       public :: num_grid_cfg_multigrid               = .false.
   logical,       public :: num_grid_cfg_no_pruning              = .false.
   logical,       public :: num_grid_cfg_zipgrid                 = .true.
   real(8),       public :: num_grid_cfg_radint                  = 1.0d-13

   private

contains

   subroutine report_num_grid()

      write(*, *) ' ===== Numerical integration grid ====='

      if (num_grid_cfg_import_grid) then

         write(*, '(a)')        '   - integration grid imported from file '//num_grid_cfg_gridfile

      else

         write(*, '(a)')        '   - radial quadrature according to '// &
                                'R. Lindh, P.-Aa. Malmqvist, and L. Gagliardi'
         write(*, '(a, e12.5)') '     precision of radial quadrature set to: ', num_grid_cfg_radint
         write(*, '(a, i2)')    '   - angular quadrature using the Lebedev scheme, '// &
                                'exact up to order L = ', num_grid_cfg_angint

         if (num_grid_cfg_no_pruning) then
            write(*, '(a)')      '   - pruning of angular grid turned off'
         end if

         if (num_grid_cfg_integration_check_level > 0) then
            write(*, '(a)')      '   - numerical integration will be tested'
         end if

         if (num_grid_cfg_estimate_radii) then
            write(*, '(a)')      '   - estimate relative atomic sizes for use in the Becke'
            write(*, '(a)')      '     partitioning scheme from atomic contributions '// &
                                 'to the density'
         end if

      end if

   end subroutine

end module
