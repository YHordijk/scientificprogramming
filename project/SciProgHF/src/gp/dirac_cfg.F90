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

module dirac_cfg

!radovan: this module can contain info about the chosen method and hamiltonian
!         the individual methods and hamiltonians can store their configuration
!         in individual modules (like for instance dft_cfg)
!         if slaves need to access the content make sure it is broadcast to them
!         (now it's not) this should be done by a method
!         which should be part of this module

  implicit none

  save

!                    kohn-sham dft
  logical, public :: dirac_cfg_dft_calculation    = .false.
!                    hartree-fock
  logical, public :: dirac_cfg_hf_calculation     = .false.
!                    hf or ks
  logical, public :: dirac_cfg_scf_calculation    = .false.
!                    MCSCF-srDFT
  logical, public :: dirac_cfg_mcsrdft_calculation = .false.
!                    openrsp module
  logical, public :: dirac_cfg_openrsp            = .false.
!                    visualization module
  logical, public :: dirac_cfg_visual             = .false.
!                    exacc module
  logical, public :: dirac_cfg_exacc              = .false.
!                    relcc module
  logical, public :: dirac_cfg_relcc              = .false.
!                    embedding module
  logical, public :: dirac_cfg_fde                = .false.
  logical, public :: dirac_cfg_fde_response       = .false.
  logical, public :: dirac_cfg_fde_export         = .false.
!                    Polarizable Continuum Model module
  logical, public :: dirac_cfg_pcm                = .false.
  private

end module
