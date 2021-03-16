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

module adc_cfg

! This module contains all input set by the user. Used for RELADC


!
!  adc control variables
!
  implicit none

  save    !  if values of variables are changed they are still there after the next call 

  integer, public                  :: reladc_adclevel     = 3
  logical, public                  :: reladc_dosips       = .false.
  logical, public                  :: reladc_dodips       = .false.
  integer, dimension(1:32), public :: reladc_sipreps      = 0  
  integer, public                  :: reladc_no_sipreps   = 0
  integer, dimension(1:32), public :: reladc_dipreps      = 0  
  logical, public                  :: reladc_readqkl      = .false.
  logical, public                  :: reladc_doconst      = .true.
  logical, public                  :: reladc_doadcpop     = .false.
  real*8,  public                  :: reladc_vconv        = 1.0E-06
  real*8,  public                  :: reladc_adcthr       = 0.0
  logical, public                  :: reladc_dodiag       = .true.
!
!  selection of diagonalizer
!
  logical, public                  :: reladc_dolanc     = .false.
  logical, public                  :: reladc_dofull     = .false.
!
!  lanczos control variables
!
  integer, public                  :: reladc_sipiter      = 500
  integer, public                  :: reladc_dipiter      = 500
  integer, dimension(32), public   :: reladc_sipeigv      = 0
  integer, dimension(32), public   :: reladc_dipeigv      = 0
  real*8,  public                  :: reladc_sipprnt      = 20.0
  real*8,  public                  :: reladc_dipprnt      = 20.0
  logical, public                  :: reladc_doincore     = .false.
  integer, public                  :: reladc_lancmem      = 50000000
!
!  davidson control variables
!
  integer, public                  :: reladc_davroots   = 5
  integer, public                  :: reladc_davmaxsp   = 50
  integer, public                  :: reladc_davmaxit   = 5
  real*8,  public                  :: reladc_davconv    = 1.0E-05
  logical, public                  :: reladc_davooc     = .true.
!
!  fanoadc variables
!
  logical, public                  :: reladc_dofano       = .false.
  logical, public                  :: reladc_fanoonly     = .false.
  logical, public                  :: reladc_fano_ovl_in  = .false.
  integer, public                  :: reladc_fano_inrep   = 0
  integer, public                  :: reladc_fano_inrelsp = 0
  integer, public                  :: reladc_fano_nrgroups = 0
  integer, public                  :: reladc_fano_nrchannels = 0
  integer, public, allocatable, dimension(:)    :: reladc_fano_groups
  integer, public, allocatable, dimension(:,:)    :: reladc_fano_channels
  character(4), public, allocatable, dimension(:) :: reladc_fano_labels
  logical, public                  :: reladc_fano_checkfin= .false.
  real(8), public                  :: reladc_fano_fovl = 0.2
  
end module
