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

module dft_cfg

   use xc_mpi

   implicit none

   public sync_dft_cfg

   logical, public :: dft_cfg_alda_hs                 = .false.
   logical, public :: dft_cfg_alda_ha                 = .false.
   logical, public :: dft_cfg_xalda                   = .false.
   logical, public :: dft_cfg_no_sdft                 = .false.
   logical, public :: dft_cfg_sdft_collinear          = .false.
   logical, public :: dft_cfg_grac                    = .false.
   logical, public :: dft_cfg_grac_is_active          = .false.
   logical, public :: dft_cfg_saop                    = .false.
   logical, public :: dft_cfg_saop_is_active          = .false.
   logical, public :: dft_cfg_saop_with_response_part = .false.
   logical, public :: dft_cfg_asymptote_is_lb94       = .false.
   logical, public :: dft_cfg_asymptote_is_lbalpha    = .false.
   logical, public :: dft_cfg_blocked                 = .false.
   logical, public :: dft_cfg_pointwise               = .false.
   logical, public :: dft_cfg_overlap_diagnostic      = .false.

   real(8), public :: dft_cfg_tinydens                = 1.0d-14
   real(8), public :: dft_cfg_screening               = 1.0d-18
   real(8), public :: dft_cfg_ac_ip                   = 0.0d0
   real(8), public :: dft_cfg_ac_threshold            = 0.0d0
   real(8), public :: dft_cfg_grac_alpha              = 0.0d0
   real(8), public :: dft_cfg_grac_beta               = 0.0d0

   private

contains

   subroutine sync_dft_cfg()

#ifdef VAR_MPI
      call xc_mpi_bcast(dft_cfg_alda_hs                 )
      call xc_mpi_bcast(dft_cfg_alda_ha                 )
      call xc_mpi_bcast(dft_cfg_xalda                   )
      call xc_mpi_bcast(dft_cfg_no_sdft                 )
      call xc_mpi_bcast(dft_cfg_sdft_collinear          )
      call xc_mpi_bcast(dft_cfg_grac                    )
      call xc_mpi_bcast(dft_cfg_grac_is_active          )
      call xc_mpi_bcast(dft_cfg_saop                    )
      call xc_mpi_bcast(dft_cfg_saop_is_active          )
      call xc_mpi_bcast(dft_cfg_saop_with_response_part )
      call xc_mpi_bcast(dft_cfg_asymptote_is_lb94       )
      call xc_mpi_bcast(dft_cfg_asymptote_is_lbalpha    )
      call xc_mpi_bcast(dft_cfg_blocked                 )
      call xc_mpi_bcast(dft_cfg_pointwise               )
      call xc_mpi_bcast(dft_cfg_overlap_diagnostic      )
      call xc_mpi_bcast(dft_cfg_tinydens                )
      call xc_mpi_bcast(dft_cfg_screening               )
      call xc_mpi_bcast(dft_cfg_ac_ip                   )
      call xc_mpi_bcast(dft_cfg_ac_threshold            )
      call xc_mpi_bcast(dft_cfg_grac_alpha              )
      call xc_mpi_bcast(dft_cfg_grac_beta               )
#endif /* ifdef VAR_MPI */

   end subroutine

end module
