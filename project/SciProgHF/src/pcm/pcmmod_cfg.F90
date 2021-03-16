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

module pcmmod_cfg

#ifdef HAS_PCMSOLVER
use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_double, c_int, c_null_char

implicit none

public pcmsolver_input

! logical block
logical, public, save :: pcmmod_host_provides_input = .false.
logical, public, save :: pcmmod_skipss   =  .false.
logical, public, save :: pcmmod_separate =  .false.
logical, public, save :: pcmmod_dospf    =  .false.
logical, public, save :: pcmmod_skipoit  =  .false.
! printlevel
integer, public, save :: pcmmod_print    = 0

! *PCMSOL section
! cavity specification *PCMSOL section
character(kind=c_char, len=8), public, save :: pcmmod_cavity_type = 'gepol  '//c_null_char
integer, public, save  :: pcmmod_patch_level = 2
real(c_double), public, save  :: pcmmod_coarsity = 0.5
real(c_double), public, save  :: pcmmod_cavity_area = 0.3
real(c_double), public, save  :: pcmmod_min_distance = 0.1
integer, public, save  :: pcmmod_der_order = 4
logical(c_bool), public, save :: pcmmod_scaling = .true.
character(kind=c_char, len=8), public, save :: pcmmod_radii_set = 'bondi  '//c_null_char
character(kind=c_char, len=20), public, save :: pcmmod_restart_name = 'cavity.npz         '//c_null_char
real(c_double), public, save  :: pcmmod_min_radius = 100.0
! solver specification *PCMSOL section
character(kind=c_char, len=7), public, save :: pcmmod_solver_type = 'iefpcm'//c_null_char
character(kind=c_char, len=16), public, save :: pcmmod_solvent = '               '//c_null_char
character(kind=c_char, len=11), public, save :: pcmmod_equation_type = 'secondkind'//c_null_char
real(c_double), public, save  :: pcmmod_correction = 0.0
real(c_double), public, save  :: pcmmod_probe_radius = 1.0
! green specification *PCMSOL section
character(kind=c_char, len=7), public, save :: pcmmod_inside_type = 'vacuum'//c_null_char
character(kind=c_char, len=22), public, save :: pcmmod_outside_type = 'uniformdielectric    '//c_null_char
real(c_double), public, save :: pcmmod_outside_epsilon = 1.0

contains

! Performs syntactic checks on PCMSolver input and fills the data structures
! holding input data
function pcmsolver_input() result(host_input)

  use pcmsolver, only: PCMInput

  type(PCMInput) :: host_input

  character(kind=c_char, len=1) :: cavity_type(8)
  character(kind=c_char, len=1) :: radii_set(8)
  character(kind=c_char, len=1) :: restart_name(20)
  character(kind=c_char, len=1) :: solver_type(7)
  character(kind=c_char, len=1) :: solvent(16)
  character(kind=c_char, len=1) :: equation_type(11)
  character(kind=c_char, len=1) :: inside_type(7)
  character(kind=c_char, len=1) :: outside_type(22)

  call pcmsolver_f2c_string(pcmmod_cavity_type, cavity_type, 8_c_int)
  host_input%cavity_type  = cavity_type
  host_input%patch_level  = int(pcmmod_patch_level, kind=c_int)
  host_input%coarsity     = pcmmod_coarsity
  host_input%area         = pcmmod_cavity_area
  host_input%min_distance = pcmmod_min_distance
  host_input%der_order    = int(pcmmod_der_order, kind=c_int)
  host_input%scaling      = pcmmod_scaling
  call pcmsolver_f2c_string(pcmmod_radii_set, radii_set, 8_c_int)
  host_input%radii_set    = radii_set
  call pcmsolver_f2c_string(pcmmod_restart_name, restart_name, 20_c_int)
  host_input%restart_name = restart_name
  host_input%min_radius   = pcmmod_min_radius

  call pcmsolver_f2c_string(pcmmod_solver_type, solver_type, 7_c_int)
  host_input%solver_type   = solver_type
  call pcmsolver_f2c_string(pcmmod_solvent, solvent, 16_c_int)
  host_input%solvent       = solvent
  call pcmsolver_f2c_string(pcmmod_equation_type, equation_type, 11_c_int)
  host_input%equation_type = equation_type
  host_input%correction    = pcmmod_correction
  host_input%probe_radius  = pcmmod_probe_radius

  call pcmsolver_f2c_string(pcmmod_inside_type, inside_type, 7_c_int)
  host_input%inside_type     = inside_type
  host_input%outside_epsilon = pcmmod_outside_epsilon
  call pcmsolver_f2c_string(pcmmod_outside_type, outside_type, 7_c_int)
  host_input%outside_type    = outside_type

end function pcmsolver_input
#endif /* HAS_PCMSOLVER */

end module
