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

module polprp_cfg
     use iso_c_binding

! This module contains all input for the POLPRP/DAVIDSON module

  implicit none

  save    !  if values of variables are changed their values are retained

!
!  polarization propagator control variables
!
  integer, dimension(1:32), public :: polprp_statesym     = 0  
  logical, public                  :: polprp_doextended   = .false.
  real*8,  public                  :: polprp_writethr     = 0.0
  logical, public                  :: polprp_dodiag       = .true.
  logical, public                  :: polprp_dotrmo       = .false.
  integer, public                  :: polprp_printlev     = 0
  logical, public                  :: polprp_skipccseti   = .false.
!
!  davidson control variables
!
  integer, public                  :: polprp_davroots   = 5
  integer, public                  :: polprp_davmaxsp   = 200
  integer, public                  :: polprp_davmaxit   = 5
  real*8,  public                  :: polprp_davconv    = 1.0E-05
  logical, public                  :: polprp_davreort   = .false.
  
end module
