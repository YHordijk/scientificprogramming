!
!   Polarizable Embedding (PE) library
!   Copyright (C) 2013, 2014 The PE library developers. See the CONTRIBUTORS file
!                            in the top-level directory of this distribution.
!
!   This file is part of the PE library.
!
!   The PE library is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as
!   published by the Free Software Foundation, either version 3 of the
!   License, or (at your option) any later version.
!
!   The PE library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the PE library. If not, see <http://www.gnu.org/licenses/>.
!
!   Contact information:
!
!   Jogvan Magnus Haugaard Olsen
!   E-mail: foeroyingur@gmail.com
!
module pe_variables

    use pe_precision

    implicit none

    ! options
    logical, save :: peqm = .false.
    logical, save :: pe_pot = .false.
    logical, save :: pe_iter = .true.
    logical, save :: pe_redthr = .false.
    logical, save :: pe_border = .false.
    logical, save :: pe_ind_damp = .false.
    logical, save :: pe_mul_damp = .false.
    logical, save :: pe_core_damp = .false.
    logical, save :: pe_gspol = .false.
    logical, save :: pe_nomb = .false.
    logical, save :: pe_polar = .false.
    logical, save :: pe_cube = .false.
    logical, save :: pe_restart = .false.
    logical, save :: pe_verbose = .false.
    logical, save :: pe_debug = .false.
    logical, save :: pe_zeromul = .false.
    logical, save :: pe_zeropol = .false.
    logical, save :: pe_isopol = .false.
    logical, save :: pe_lf = .false.
    logical, save :: pe_field = .false.

    ! runtypes
    logical :: fock = .false.
    logical :: energy = .false.
    logical :: response = .false.
    logical :: molgrad = .false.
    logical :: mep = .false.
    logical :: london = .false.
    logical :: effdipole = .false.
    logical :: cube = .false.

    ! in- and output logical, i.e. denmats, fckmats and expvals
    logical :: lden = .false.
    logical :: lfck = .false.
    logical :: lexp = .false.

    ! matrix type
    logical :: trimat = .true.

    ! filenames
    character(len=80) :: potfile = 'POTENTIAL.INP'

    logical, save :: synced = .false.
    logical, save :: initialized = .false.
    integer, save :: master, myid, nprocs, ierr
    integer, save :: site_start, site_finish
    integer, save :: cube_start, cube_finish
    integer, dimension(:), save, allocatable :: siteloops
    integer, dimension(:), save, allocatable :: cubeloops
    integer, dimension(:), save, allocatable :: poldists
    integer, dimension(:), save, allocatable :: sitedists
    integer, dimension(:), save, allocatable :: cubedists
    integer, dimension(:), save, allocatable :: displs

    ! logical unit for output file (default is stdout)
    integer, save :: luout = 6

    ! integer used to check memory allocation status,
    ! i.e. allocate(array(len), stat=mem_stat)
    integer :: mem_stat

    ! dummy variables
    integer :: dummy_int
    real(dp) :: dummy_real
    integer, dimension(1) :: dummy_int_array
    real(dp), dimension(1) :: dummy_real_array

    ! constants, thresholds and stuff
    real(dp), parameter :: zero = 1.0e-8_dp
    integer, save :: scf_cycle = 0
    integer, save :: redlvl = 0
    real(dp), save :: thriter = 1.0e-8_dp
    real(dp), save :: ind_damp = 2.1304_dp
    real(dp), save :: mul_damp = 2.1304_dp
    real(dp), save :: core_damp = 2.1304_dp
    real(dp), save :: Rmin = 2.2_dp
    integer, save :: nredist = 1
    integer, save :: redist_order = 1
    character(len=6), save :: border_type = 'REDIST'
    ! use Cholesky factorization of classical response matrix
    logical, save :: chol = .true.
    ! the order from which multipoles are zeroed
    integer, save :: zeromul_order = 1

    ! polarizable embedding potential info
    ! ------------------------------------
    ! total number of classical sites
    integer, save :: nsites = 0
    ! number of polarizable sites
    integer, save :: npols = 0
    ! exclusion list length
    integer, save :: lexlst = 0
    ! number of density matrices
    integer :: ndens = 0
    ! number of basis functions in core fragment
    integer, save :: nbas
    ! size of packed matrices
    integer, save :: nnbas
    ! size of full matrices
    integer, save :: n2bas
    ! number of nuclei in core region
    integer, save :: nnucs = 0
    ! number of LJ sites in core region
    integer, save :: qmLJsites = 0

    ! Thole qm damping stuff
    ! ----------------------
    integer, save :: nqmpols = 0
    real(dp), dimension(:,:), allocatable, save :: qmpols
    real(dp), save :: qmdamp = 2.1304_dp

    ! specifies what type of parameters are present
    ! lmul(0): monopoles, lmul(1): dipoles etc.
    logical, dimension(0:5), save :: lmul = .false.
    ! lpol(1): dipole-dipole polarizabilities
    logical, dimension(0:2,0:2), save :: lpol = .false.
    ! lhypol(1): dipole-dipole-dipole polarizabilities/1st hyperpolarizability
!    logical, dimension(1), save :: lhypol
!    ! lvdw: LJ parameters
    logical, save :: lvdw = .false.

    ! charges, areas, coordinates, elements and exclusion lists
    ! site nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zs
    ! core fragment nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zm
    ! site coordinates
    real(dp), dimension(:,:), allocatable, save :: Rs
    ! core fragment nuclear coordinates
    real(dp), dimension(:,:), allocatable, save :: Rm
    ! core center-of-mass
    real(dp), dimension(:), allocatable, save :: core_com
    ! core polarizabilities
    real(dp), dimension(:,:), allocatable, save :: core_alphas
    ! site elements
    character(len=2), dimension(:), allocatable, save :: elems
    ! exclusion list
    integer, dimension(:,:), allocatable, save :: exclists

    ! energy contributions
    ! total
    real(dp), dimension(:), allocatable, save :: Epe
    ! electrostatic
    real(dp), dimension(:,:), allocatable, save :: Ees
    ! polarization
    real(dp), dimension(:,:), allocatable, save :: Epol
    ! LJ
    real(dp), dimension(:), allocatable, save :: Elj

    ! multipole moments
    ! monopoles, dipoles, quadrupoles, octopoles, etc.
    real(dp), dimension(:,:), allocatable, save :: M0s, M1s, M2s, M3s, M4s, M5s
    ! polarizabilities
    ! symmetric dipole-dipole polarizabilities
    real(dp), dimension(:,:), allocatable, save :: P11s
    ! .true. if P11 > 0 else .false.
    logical, dimension(:), allocatable, save :: zeroalphas

    ! LJ parameters for MM region - from pot file
    real(dp), dimension(:,:), allocatable, save :: LJs
    ! LJ paramters for QM region - from dal file
    real(dp), dimension(:,:), allocatable, save :: qmLJs

    ! CUBE stuff
    ! ---------

    ! options for CUBE
    ! calculate electric field
    logical, save :: cube_field = .false.

    ! general cube information
    ! number of grid points
    integer, save :: npoints
    ! grid points
    real(dp), dimension(:,:), allocatable, save :: Rg
    ! CUBE file origin and step sizes
    real(dp), dimension(3), save :: origin, step
    ! grid density in x, y and z direction
    integer, save :: xgrid = 6
    integer, save :: ygrid = 6
    integer, save :: zgrid = 6
    ! number of steps in x, y and z direction
    integer, save :: xsteps
    integer, save :: ysteps
    integer, save :: zsteps
    ! box size relative to molecule size
    real(dp), save :: xsize = 8.0_dp
    real(dp), save :: ysize = 8.0_dp
    real(dp), save :: zsize = 8.0_dp

    ! Internal field stuff and locfld stuff
    ! Coordinates on which potential and field are calculated
    real(dp), dimension(:,:), allocatable, save :: crds
    ! Number of coordinates (length of crds) / 3
    integer, save :: ncrds = 0

    ! External field stuff
    ! Field strengths
    real(dp), dimension(3), save :: efields = 0.0_dp

end module pe_variables
