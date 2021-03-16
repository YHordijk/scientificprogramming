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

module x2cmod_cfg

  implicit none

  save

! parameters

  integer, parameter, public :: x2c_funit       =  25    ! file unit for X2CMAT
  integer, parameter, public :: x2c_fscru       =  26    ! file unit for x2cscr
  integer, parameter, public :: x2c_fh1dirty    =  98    ! file unit for x2c_h1_dirty
  integer, parameter, public :: x2c_Sao_scrF    =  88
  integer, parameter, public :: x2c_lowdin_scrF =  99    ! file unit for scratch file containing the 4c-Lowdin matrix
  integer, parameter, public :: x2c_famfi       =  77
  integer, parameter, public :: x2c_mxatom      =  20    ! keep in syc with mxatom in mxcent.h
  integer, parameter, public :: x2c_mxopen      =  10    ! keep in syc with mxopne in dcbdhf.h
  integer, parameter, public :: x2c_mxcent      =  300   ! keep in syc with mxcent in mxcent.h

! character block
  character (len=6), public :: file_name_1e_fock_matrix   = 'undef '
  character (len=6), public :: file_name_2e_fock_matrices = 'undef '

! logical block
  logical, public :: x2cmod_x2c                       =  .false.
  logical, public :: x2cmod_fragment_x2c              =  .false.
  logical, public :: x2cmod_debug                     =  .false.
  logical, public :: x2cmod_h1dirty                   =  .false.
  logical, public :: x2cmod_skip_op_pct               =  .false.
  logical, public :: x2cmod_mmf_restart               =  .false.
  logical, public :: x2c_linsym                       =  .false.
  logical, public :: x2c_prep_bnccorr                 =  .false.
  logical, public :: x2c_do_spherical                 =  .false.
  logical, public :: x2c_free_part_mtx_defh1          =  .false.
  logical, public :: x2c_1e_system                    =  .false.
  logical, public :: x2c_4c_fock_mtx_defh1            =  .false.
  logical, public :: x2c_mdirac                       =  .false. ! .true. if MDIRAC in DIRAC is active
  logical, public :: x2c_spinfree                     =  .false. ! .true. if spinfree in DIRAC is active
  logical, public :: x2cmod_mmf                       =  .false. ! this implies x2c_4c_fock_mtx_defh1 = .true. (mmf-ansatz)
  logical, public :: x2c_2c_mmf_mos                   =  .false. ! this implies x2c_4c_fock_mtx_defh1 = .true. (mmf-ansatz)
                                                                 ! but i think it's much nicer (more intuitive name for programmers) 
                                                                 ! to have this logical on the outside of the X2Cmod - stefan, may 2011
  logical, public :: x2c_fock_saao_basis              =  .false.

! real(8) block
  real(8), public :: x2c_cspeed                       = 0.0d0
  real(8), public :: x2c_dfopen(0:x2c_mxopen)         = 0.0d0

! integer block
  integer, public :: x2c_prt                          =  0 ! print level
  integer, public :: x2c_is_defining_h1mat            = -1 ! no defining h1 hamiltonian: decision is made within the program
  integer, public :: x2c_add_amfi                     =  1 ! current default: use the 'old' AMFI by B. Schimmelpfennig
  integer, public :: x2c_amfi_order                   =  2 ! default (if AMFI is turned on): spin-same orbit contributions
  integer, public :: dim_pshell(2)                    =  0
  integer, public :: dim_eshell(2)                    =  0
  integer, public :: dim_eshell2(2)                   =  0
  integer, public :: dim_e2shell(2)                   =  0
  integer, public :: dim_e2shellq(2)                  =  0
  integer, public :: dim_oshell(2)                    =  0
  integer, public :: dim_o2shellq(2)                  =  0
  integer, public :: fullemat_dim                     =  0
  integer, public :: fullomat_dim                     =  0
  integer, public :: fulleomat_dim                    =  0
  integer, public :: fulllowdmat_dim                  =  0
  integer, public :: fullao2momat_dim                 =  0
  integer, public :: fullmo2momat_dim                 =  0
  integer, public :: nzt_x2c                          =  0
  integer, public :: nr_quat                          =  0
  integer, public :: nr_fsym                          =  0
  integer, public :: mfsym                            =  0
  integer, public :: nr_ao_total_x2c                  =  0
  integer, public :: nr_ao_large_x2c                  =  0
  integer, public :: nr_ao_small_x2c                  =  0
  integer, public :: n2bastq_dim_x2c                  =  0
  integer, public :: nr_cmo_q                         =  0
  integer, public :: nr_ao_all(2)                     =  0
  integer, public :: nr_ao_l(2)                       =  0
  integer, public :: nr_ao_s(2)                       =  0
  integer, public :: nr_mo_lw_all(2)                  =  0
  integer, public :: nr_mo_lw_l(2)                    =  0
  integer, public :: nr_mo_lw_s(2)                    =  0
  integer, public :: ioff_aomat_x(2,2)                =  0
  integer, public :: x2c_cb_pq_to_uq(4, 0:7)          =  0
  integer, public :: x2c_cb_uq_to_pq(4, 0:7)          =  0
  integer, public :: x2c_bs_to_fs(0:7, 2)             =  0
  integer, public :: x2c_qdef(1:4)                    =  0
  integer, public :: x2c_tot_charge_for_mfsum         =  0
  integer, public :: fh_1int_4c                       =  0
  integer, public :: all_prp_op                       =  0
  integer, public :: num_nuclei                       =  0
  integer, public :: type_nuclei                      =  0
  integer, public :: nr_ao_bas_type(1:x2c_mxatom,0:2) =  0
  integer, public :: type_charge(1:x2c_mxatom)        =  0
  integer, public :: nr_symm_indep_cent(1:x2c_mxatom) =  0
  integer, public :: nr_degen_nuc_cent(1:x2c_mxcent)  =  0
  integer, public :: x2c_pointer_quat(0:7, 1:2)       =  0
  integer, public :: x2c_bs_irrep_mat(1:4, 0:7)       =  0
  integer, public :: x2c_iqmult_trip_q(1:4, 1:4, 1:4) =  0
  integer, public :: x2c_pointer_quat_op(0:7)         =  0
  integer, public :: x2c_max_quant_num                =  0
  integer, public :: scf_iter_counter                 =  0
  integer, public :: diis_counter                     =  0
  integer, public :: nr_2e_fock_matrices              =  0

end module
