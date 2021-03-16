module fde_dirac_matrices_integration

!embed modules
   use fde_types
   use fde_cfg
   use fde_io
   use fde_data
   use fde_nadd_derv
   use fde_xcfun_interface
   use fde_max_block_length
   use electrostatic_potential

!gosia: rm later:
   use fde_mag_cfg

! dirac-specific modules
   use interface_ao
   use interface_mo
   use interface_file_io
   use dft_cfg
   use interface_grid
   use dirac_ao_eval
   use density_eval
   use fde_mpi

   use xc_london_c1

   implicit none

   public fde_dirac_emb_matrices_via_integration
   public fde_dirac_emb_pot_energy

   public fde_dirac_interaction_nadd
   public fde_dirac_interaction_el_b_rho_a
   public fde_dirac_interaction_nuc_a_rho_b

   public fde_susc2el_integrated

   private

   type fde_omega_prefactor
      real(8), allocatable :: n(:)
      real(8), allocatable :: gn(:, :)
      real(8), allocatable :: s(:, :)
      real(8), allocatable :: gs(:, :, :)
      real(8), allocatable :: n_b_lao(:)
      real(8), allocatable :: gn_b_lao(:, :)
      real(8), allocatable :: s_b_lao(:, :)
      real(8), allocatable :: gs_b_lao(:, :, :)
      real(8), allocatable :: tn(:)
      real(8), allocatable :: ts(:, :)
   end type

   integer, parameter :: bllen = 1 !later more, now slowdft
   integer            :: block_length

   save

   integer, parameter :: file_unit      = 6

   real(8), allocatable :: rx(:)
   real(8), allocatable :: ry(:)
   real(8), allocatable :: rz(:)
   real(8), allocatable :: rw(:)

   real(8), allocatable :: frozen_n(:)
   real(8), allocatable :: frozen_gnx(:)
   real(8), allocatable :: frozen_gny(:)
   real(8), allocatable :: frozen_gnz(:)
   real(8), allocatable :: frozen_elpot(:)
! -------------------------------------------
! for perturbed density import
   real(8), allocatable :: frozen_n_bx(:)
   real(8), allocatable :: frozen_n_by(:)
   real(8), allocatable :: frozen_n_bz(:)
   real(8), allocatable :: frozen_gxn_bx(:)
   real(8), allocatable :: frozen_gxn_by(:)
   real(8), allocatable :: frozen_gxn_bz(:)
   real(8), allocatable :: frozen_gyn_bx(:)
   real(8), allocatable :: frozen_gyn_by(:)
   real(8), allocatable :: frozen_gyn_bz(:)
   real(8), allocatable :: frozen_gzn_bx(:)
   real(8), allocatable :: frozen_gzn_by(:)
   real(8), allocatable :: frozen_gzn_bz(:)
   real(8), allocatable :: frozen_sx_bx(:)
   real(8), allocatable :: frozen_sx_by(:)
   real(8), allocatable :: frozen_sx_bz(:)
   real(8), allocatable :: frozen_sy_bx(:)
   real(8), allocatable :: frozen_sy_by(:)
   real(8), allocatable :: frozen_sy_bz(:)
   real(8), allocatable :: frozen_sz_bx(:)
   real(8), allocatable :: frozen_sz_by(:)
   real(8), allocatable :: frozen_sz_bz(:)
   real(8), allocatable :: frozen_gxsx_bx(:)
   real(8), allocatable :: frozen_gxsx_by(:)
   real(8), allocatable :: frozen_gxsx_bz(:)
   real(8), allocatable :: frozen_gxsy_bx(:)
   real(8), allocatable :: frozen_gxsy_by(:)
   real(8), allocatable :: frozen_gxsy_bz(:)
   real(8), allocatable :: frozen_gxsz_bx(:)
   real(8), allocatable :: frozen_gxsz_by(:)
   real(8), allocatable :: frozen_gxsz_bz(:)
   real(8), allocatable :: frozen_gysx_bx(:)
   real(8), allocatable :: frozen_gysx_by(:)
   real(8), allocatable :: frozen_gysx_bz(:)
   real(8), allocatable :: frozen_gysy_bx(:)
   real(8), allocatable :: frozen_gysy_by(:)
   real(8), allocatable :: frozen_gysy_bz(:)
   real(8), allocatable :: frozen_gysz_bx(:)
   real(8), allocatable :: frozen_gysz_by(:)
   real(8), allocatable :: frozen_gysz_bz(:)
   real(8), allocatable :: frozen_gzsx_bx(:)
   real(8), allocatable :: frozen_gzsx_by(:)
   real(8), allocatable :: frozen_gzsx_bz(:)
   real(8), allocatable :: frozen_gzsy_bx(:)
   real(8), allocatable :: frozen_gzsy_by(:)
   real(8), allocatable :: frozen_gzsy_bz(:)
   real(8), allocatable :: frozen_gzsz_bx(:)
   real(8), allocatable :: frozen_gzsz_by(:)
   real(8), allocatable :: frozen_gzsz_bz(:)
! -------------------------------------------
 
   real(8), allocatable :: use_vemb(:)
   
! -------------------------------------------
! for electrostatic, density calculation 
   real(8), allocatable :: active_elpot(:)
   real(8), allocatable :: active_nucpot(:)
   real(8), allocatable :: active_n(:)   
   real(8), allocatable :: active_gnx(:) 
   real(8), allocatable :: active_gny(:) 
   real(8), allocatable :: active_gnz(:) 
   real(8), allocatable :: active_hnxx(:)
   real(8), allocatable :: active_hnxy(:)
   real(8), allocatable :: active_hnxz(:)
   real(8), allocatable :: active_hnyy(:)
   real(8), allocatable :: active_hnyz(:)
   real(8), allocatable :: active_hnzz(:)
! -------------------------------------------
! for perturbed density export
   real(8), allocatable :: active_n_bx(:)
   real(8), allocatable :: active_n_by(:)
   real(8), allocatable :: active_n_bz(:)
   real(8), allocatable :: active_gxn_bx(:)
   real(8), allocatable :: active_gxn_by(:)
   real(8), allocatable :: active_gxn_bz(:)
   real(8), allocatable :: active_gyn_bx(:)
   real(8), allocatable :: active_gyn_by(:)
   real(8), allocatable :: active_gyn_bz(:)
   real(8), allocatable :: active_gzn_bx(:)
   real(8), allocatable :: active_gzn_by(:)
   real(8), allocatable :: active_gzn_bz(:)
   real(8), allocatable :: active_sx_bx(:)
   real(8), allocatable :: active_sx_by(:)
   real(8), allocatable :: active_sx_bz(:)
   real(8), allocatable :: active_sy_bx(:)
   real(8), allocatable :: active_sy_by(:)
   real(8), allocatable :: active_sy_bz(:)
   real(8), allocatable :: active_sz_bx(:)
   real(8), allocatable :: active_sz_by(:)
   real(8), allocatable :: active_sz_bz(:)
   real(8), allocatable :: active_gxsx_bx(:)
   real(8), allocatable :: active_gxsx_by(:)
   real(8), allocatable :: active_gxsx_bz(:)
   real(8), allocatable :: active_gxsy_bx(:)
   real(8), allocatable :: active_gxsy_by(:)
   real(8), allocatable :: active_gxsy_bz(:)
   real(8), allocatable :: active_gxsz_bx(:)
   real(8), allocatable :: active_gxsz_by(:)
   real(8), allocatable :: active_gxsz_bz(:)
   real(8), allocatable :: active_gysx_bx(:)
   real(8), allocatable :: active_gysx_by(:)
   real(8), allocatable :: active_gysx_bz(:)
   real(8), allocatable :: active_gysy_bx(:)
   real(8), allocatable :: active_gysy_by(:)
   real(8), allocatable :: active_gysy_bz(:)
   real(8), allocatable :: active_gysz_bx(:)
   real(8), allocatable :: active_gysz_by(:)
   real(8), allocatable :: active_gysz_bz(:)
   real(8), allocatable :: active_gzsx_bx(:)
   real(8), allocatable :: active_gzsx_by(:)
   real(8), allocatable :: active_gzsx_bz(:)
   real(8), allocatable :: active_gzsy_bx(:)
   real(8), allocatable :: active_gzsy_by(:)
   real(8), allocatable :: active_gzsy_bz(:)
   real(8), allocatable :: active_gzsz_bx(:)
   real(8), allocatable :: active_gzsz_by(:)
   real(8), allocatable :: active_gzsz_bz(:)
! -------------------------------------------

   real(8), allocatable :: dmat_sorted(:, :)

   real(8) :: nr_electrons_integrated
   real(8), save :: fde_energy
   real(8), save :: fde_mat_energy
   real(8), save :: non_additive_interaction
   real(8), save :: electrostatic_b_rho_a_interaction 
   real(8), save :: electrostatic_a_rho_b_interaction 

   real(8) :: susc2el_fde_integrated(3, 3)
   real(8) :: susc2el_fde_integrated_sdft(3, 3)

   logical :: use_gga_qr = .false. !quaternion real
   logical :: use_gga_qi = .false. !quaternion imaginary
   logical :: use_gga    = .false.

!             number of points in current batch
   integer :: nr_points_batch

!             number of points currently looped over on this processor
   integer :: nr_points_on_this_proc

!             total nr of points (sum of nr_points_batch)
   integer :: nr_points_total

!             total nr of points that were actually used
!             we skip points with very small density
!             typically smaller than nr_points_total
   integer :: nr_points_used

   integer :: max_ao_g_order
   integer :: max_ao_m_order
   integer :: max_fun_derv

   real(8), external :: second
   real(8)           :: time_integration, time_integration_start
   real(8)           :: time_ao,          time_ao_start
   real(8)           :: time_matrix_dist, time_matrix_dist_start
   real(8)           :: time_density,     time_density_start
   real(8)           :: time_derv,        time_derv_start

   integer                      :: mat_dim
   integer                      :: nz

   real(8), pointer             :: dmat_0(:, :, :)
   real(8), allocatable, target :: dmat_0_container(:, :, :)

   real(8), pointer             :: pertden(:, :)

   real(8), allocatable         :: tmat_1(:, :)
   real(8), allocatable         :: tmat_2(:, :)

   integer                      :: nr_dmat
   real(8), pointer             :: dmat(:, :, :, :)
   real(8), allocatable, target :: dmat_container(:, :, :, :)

   real(8), pointer             :: dmat_pertden_direct(:, :, :)
   real(8), allocatable, target :: dmat_pertden_direct_container(:, :, :)

   real(8), pointer             :: dmat_pertden_reorth(:, :, :, :)
   real(8), allocatable, target :: dmat_pertden_reorth_container(:, :, :, :)

   integer                      :: nr_fmat
   real(8), pointer             :: fmat(:, :, :, :)
   real(8), allocatable, target :: fmat_container(:, :, :, :)

   integer, pointer             :: dmat_pg_sym(:)
   integer, allocatable, target :: dmat_pg_sym_container(:)
   integer, pointer             :: dmat_ih_sym(:)
   integer, allocatable, target :: dmat_ih_sym_container(:)
   integer, pointer             :: fmat_pg_sym(:)
   integer, allocatable, target :: fmat_pg_sym_container(:)

   integer                      :: nr_atoms
   real(8), pointer             :: property_gradient(:)
   real(8), allocatable, target :: property_gradient_container(:)

   type(fde_import) :: itmp

   logical :: do_potential
   logical :: use_potential
   logical :: do_interaction_energy
   logical :: do_electrostatic_pot
   integer :: response_order_mo
   integer :: response_order_ao
   logical :: do_london_rhs_direct_der1
   logical :: do_london_rhs_direct_der2
   logical :: do_london_rhs_reorth_coupling
   logical :: do_london_rhs_direct_coupling
   logical :: do_london_rhs_ro
   logical :: do_london_susc2el_der1
   logical :: do_london_susc2el_der2
   logical :: do_london_susc2el_coupling
   logical :: do_london_suscreo_coupling
   logical :: do_geo
   logical :: do_geo_0
   logical :: do_geo_1
   logical :: do_geo_2
   logical :: do_london_ks
   logical :: do_london_lr
   logical :: do_overlap_diagnostic
   logical :: do_fde_update_vemb
   logical :: do_fde_response
   logical :: do_export_pertden_reorth
   logical :: do_export_pertden_direct

   type(fde_grid) :: active_grid
   type(grid_function):: active_gridfunction

   real(8), allocatable :: elpot(:)

   integer, allocatable, target :: default_fmat_pg_sym(:)
   integer, allocatable, target :: default_dmat_pg_sym(:)
   integer, allocatable, target :: default_dmat_ih_sym(:)

contains

   function fde_dirac_emb_pot_energy()
      real(8) :: fde_dirac_emb_pot_energy 
      fde_dirac_emb_pot_energy = fde_mat_energy
   end function

   function fde_dirac_interaction_nadd()
      real(8) :: fde_dirac_interaction_nadd
      fde_dirac_interaction_nadd = non_additive_interaction
   end function 

   function fde_dirac_interaction_el_b_rho_a()
      real(8) :: fde_dirac_interaction_el_b_rho_a
      fde_dirac_interaction_el_b_rho_a = electrostatic_b_rho_a_interaction
   end function 

   function fde_dirac_interaction_nuc_a_rho_b()
      real(8) :: fde_dirac_interaction_nuc_a_rho_b
      fde_dirac_interaction_nuc_a_rho_b = electrostatic_a_rho_b_interaction
   end function 

   function fde_susc2el_integrated() result(r)
      real(8) :: r(3, 3)
      integer :: i, j
      write(*, *) 'FDE contributions - details:'
      do i = 1, 3
        do j = 1, 3
          write(*, *) 'j, i, susc2el_fde_integrated     ', j, i, susc2el_fde_integrated(j,i)
          if (.not. fde_cfg_no_sdft) then
            write(*, *) 'j, i, susc2el_fde_integrated_sdft', j, i, susc2el_fde_integrated_sdft(j,i)
          end if
        end do
      end do
      r = 0.0d0
      do i = 1,3
        do j = 1,3
          r(j, i) = r(j, i) + susc2el_fde_integrated(j, i)
          if (.not. fde_cfg_no_sdft) then
            r(j, i) = r(j, i) + susc2el_fde_integrated_sdft(j, i)
          end if
        end do
      end do
   end function

!   ----------------------------------------------------------------------------
    subroutine fde_dirac_emb_matrices_via_integration(  &
                             fde_mat_dim,               &
                             fde_nz,                    &
                             fde_dmat_0,                &
                             fde_nr_dmat,               &
                             fde_nr_fmat,               &
                             fde_dmat,                  &
                             fde_dmat_pertden_direct,   &
                             fde_dmat_pertden_reorth,   &
                             fde_fmat,                  &
                             fde_dmat_ih_sym,           &
                             fde_dmat_pg_sym,           &
                             fde_fmat_pg_sym,           &
                             fde_nr_atoms,              &
                             fde_property_gradient,     &
                             fde_do_potential,          &
                             fde_use_potential,          &
                             fde_do_interaction_energy, &
                             fde_do_electrostatic_pot,  &
                             fde_do_geo_0,              &
                             fde_do_geo_1,              &
                             fde_do_geo_2,              &
                             fde_do_london_ks,          &
                             fde_do_london_lr,          &
                             fde_response_order_mo,     &
                             fde_response_order_ao,     &
                             fde_do_london_rhs_direct_der1,   &
                             fde_do_london_rhs_direct_der2,   &
                             fde_do_london_rhs_reorth_coupling, &
                             fde_do_london_rhs_direct_coupling, &
                             fde_do_london_rhs_ro,      &
                             fde_do_london_susc2el_der1,     &
                             fde_do_london_susc2el_der2,     &
                             fde_do_london_susc2el_coupling,   &
                             fde_do_london_suscreo_coupling,   &
                             fde_do_export_pertden_reorth,     &
                             fde_do_export_pertden_direct,     &
                             fde_active_grid,           &
                             fde_active_gridfunction)

! response_order_mo/_ao:  0 - do nothing (no response)
!                         1 - ks lr
!                         2 - ks qr
!                         *_mo - old code, *_ao - openrsp code
!
! do_potential            = .true.  - do embedding potential
! do_interaction_energy   = .true.  - do contributions to interaction energy involving the active density and frozen potential
!
! contributions to differentiated kohn-sham matrix:
! do_london_rhs_direct       = .true.  - xc contributions that require mag dervs of ao-s, "direct" (if london)
! do_london_rhs_ro           = .true.  - xc contributions to reotrhonormalization terms (if london)

!     --------------------------------------------------------------------------
      integer,                   intent(in)    :: fde_mat_dim
      integer,                   intent(in)    :: fde_nz
      real(8),           target, intent(in)    :: fde_dmat_0(fde_mat_dim, fde_mat_dim, fde_nz)
      integer,                   intent(in)    :: fde_nr_dmat
      integer,                   intent(in)    :: fde_nr_fmat
      real(8), optional, target, intent(in)    :: fde_dmat(fde_mat_dim, fde_mat_dim, fde_nz, fde_nr_dmat)
      real(8), optional, target, intent(in)    :: fde_dmat_pertden_direct(fde_mat_dim, fde_mat_dim, fde_nz)
      real(8), optional, target, intent(in)    :: fde_dmat_pertden_reorth(fde_mat_dim, fde_mat_dim, fde_nz, 3)
!     real(8), optional, target, intent(inout) :: fde_fmat(fde_mat_dim, fde_mat_dim, fde_nz, fde_nr_fmat)
      real(8), optional, target                :: fde_fmat(fde_mat_dim, fde_mat_dim, fde_nz, fde_nr_fmat)
      integer, optional, target, intent(in)    :: fde_dmat_ih_sym(fde_nr_dmat)
      integer, optional, target, intent(in)    :: fde_dmat_pg_sym(fde_nr_dmat)
      integer, optional, target, intent(in)    :: fde_fmat_pg_sym(fde_nr_fmat)

      integer, optional,         intent(in)    :: fde_nr_atoms
      real(8), optional, target, intent(inout) :: fde_property_gradient(*)

      logical, optional,         intent(in)    :: fde_do_potential
      logical, optional,         intent(in)    :: fde_use_potential
      logical, optional,         intent(in)    :: fde_do_interaction_energy
      logical, optional,         intent(in)    :: fde_do_electrostatic_pot
      logical, optional,         intent(in)    :: fde_do_geo_0
      logical, optional,         intent(in)    :: fde_do_geo_1
      logical, optional,         intent(in)    :: fde_do_geo_2
      logical, optional,         intent(in)    :: fde_do_london_ks
      logical, optional,         intent(in)    :: fde_do_london_lr
      integer, optional,         intent(in)    :: fde_response_order_mo
      integer, optional,         intent(in)    :: fde_response_order_ao
      logical, optional,         intent(in)    :: fde_do_london_rhs_direct_der1
      logical, optional,         intent(in)    :: fde_do_london_rhs_direct_der2
      logical, optional,         intent(in)    :: fde_do_london_rhs_reorth_coupling
      logical, optional,         intent(in)    :: fde_do_london_rhs_direct_coupling
      logical, optional,         intent(in)    :: fde_do_london_rhs_ro
      logical, optional,         intent(in)    :: fde_do_london_susc2el_der1
      logical, optional,         intent(in)    :: fde_do_london_susc2el_der2
      logical, optional,         intent(in)    :: fde_do_london_susc2el_coupling
      logical, optional,         intent(in)    :: fde_do_london_suscreo_coupling
      logical, optional,         intent(in)    :: fde_do_export_pertden_reorth
      logical, optional,         intent(in)    :: fde_do_export_pertden_direct
!     --------------------------------------------------------------------------
      type(fde_grid),      optional, intent(in):: fde_active_grid
      type(grid_function), optional, intent(in):: fde_active_gridfunction
!     --------------------------------------------------------------------------
      real(8), external :: ddot
      integer, parameter :: max_mat_number = 10
      integer, target    :: default_pg_sym(max_mat_number)
      integer, target    :: default_ih_sym(max_mat_number)
      real(8)             :: save_ddot1, save_ddot2, h
!     --------------------------------------------------------------------------
      integer, parameter :: max_response_order_mo = 5
      integer :: i, j, iz, imat
       
      call fde_get_import_info(itmp)


!#define DEBUG_XC
#ifdef DEBUG_XC
      call matrix_to_file('fde_dmat0', fde_mat_dim*fde_mat_dim, fde_dmat_0)
#endif

!     start timer
      time_integration_start = second()
!     reset timer
      time_ao          = 0.0d0
      time_matrix_dist = 0.0d0
      time_density     = 0.0d0
      time_derv        = 0.0d0

      nr_electrons_integrated = 0.0d0
      fde_energy    = 0.0d0

      non_additive_interaction = 0.0d0
      electrostatic_b_rho_a_interaction = 0.0d0
      electrostatic_a_rho_b_interaction = 0.0d0

      susc2el_fde_integrated = 0.0d0
      susc2el_fde_integrated_sdft = 0.0d0

      parallel_fde = (fde_mpi_get_nr_proc() > 1)
#ifdef VAR_MPI
      if (parallel_fde) then
         call sync_fde_cfg()
      end if
#endif

      if (fde_mpi_is_master()) then
         mat_dim = fde_mat_dim
         nz      = fde_nz
         nr_dmat = fde_nr_dmat
         nr_fmat = fde_nr_fmat
         print *, 'debug emb: mat_dim,nz,nr_dmat,nr_fmat:',mat_dim,nz,nr_dmat,nr_fmat
      end if

!     set defaults
      do_potential          = .false.
      use_potential         = .false.
      do_interaction_energy = .false.
      do_electrostatic_pot = .false.
      response_order_mo     = 0
      response_order_ao     = 0
      do_london_rhs_direct_der1    = .false.
      do_london_rhs_direct_der2    = .false.
      do_london_rhs_reorth_coupling  = .false.
      do_london_rhs_direct_coupling  = .false.
      do_london_rhs_ro             = .false.
      do_london_susc2el_der1       = .false.
      do_london_susc2el_der2       = .false.
      do_london_susc2el_coupling   = .false.
      do_london_suscreo_coupling   = .false.
      do_geo_0              = .false.
      do_geo_1              = .false.
      do_geo_2              = .false.
      do_london_ks          = .false.
      do_london_lr          = .false.
      do_export_pertden_reorth     = .false.
      do_export_pertden_direct     = .false.

      nr_atoms              = 0

     if (fde_mpi_is_master()) then      
!        change defaults
         if (present(fde_do_potential))          do_potential          = fde_do_potential
         if (present(fde_use_potential))         use_potential         = fde_use_potential
         if (present(fde_do_interaction_energy)) do_interaction_energy = fde_do_interaction_energy
         if (present(fde_do_electrostatic_pot))  do_electrostatic_pot  = fde_do_electrostatic_pot
         if (present(fde_response_order_mo))     response_order_mo     = fde_response_order_mo
         if (present(fde_response_order_ao))     response_order_ao     = fde_response_order_ao
         if (present(fde_do_london_rhs_direct_der1))    do_london_rhs_direct_der1    = fde_do_london_rhs_direct_der1
         if (present(fde_do_london_rhs_direct_der2))    do_london_rhs_direct_der2    = fde_do_london_rhs_direct_der2
         if (present(fde_do_london_rhs_reorth_coupling))  do_london_rhs_reorth_coupling  = fde_do_london_rhs_reorth_coupling
         if (present(fde_do_london_rhs_direct_coupling))  do_london_rhs_direct_coupling  = fde_do_london_rhs_direct_coupling
         if (present(fde_do_london_rhs_ro))             do_london_rhs_ro             = fde_do_london_rhs_ro
         if (present(fde_do_london_susc2el_der1))       do_london_susc2el_der1       = fde_do_london_susc2el_der1
         if (present(fde_do_london_susc2el_der2))       do_london_susc2el_der2       = fde_do_london_susc2el_der2
         if (present(fde_do_london_susc2el_coupling))   do_london_susc2el_coupling   = fde_do_london_susc2el_coupling
         if (present(fde_do_london_suscreo_coupling))   do_london_suscreo_coupling   = fde_do_london_suscreo_coupling
         if (present(fde_do_export_pertden_reorth))     do_export_pertden_reorth     = fde_do_export_pertden_reorth
         if (present(fde_do_export_pertden_direct))     do_export_pertden_direct     = fde_do_export_pertden_direct
         if (present(fde_do_geo_0))              do_geo_0              = fde_do_geo_0
         if (present(fde_do_geo_1))              do_geo_1              = fde_do_geo_1
         if (present(fde_do_geo_2))              do_geo_2              = fde_do_geo_2
         if (present(fde_do_london_ks))          do_london_ks          = fde_do_london_ks
         if (present(fde_do_london_lr))          do_london_lr          = fde_do_london_lr
         if (present(fde_nr_atoms))              nr_atoms              = fde_nr_atoms

! aspg, 16/10/2015
!
! the defaults for the fde grid and gridfunction are, respectively, the grid and
! gridfunction variables containing the data for the frozen subsystem.
!
! these defaults should be applicable to all actions having to do with adding
! contributions to (perturbed) fock matrices, such as calculating embedding
! potential, kernel and higher-order contributions, shieldings etc.
!
! one can also give the grid and gridfunctions for the export of density and
! electrostatic potential
!
!  
         if (present(fde_active_grid)) then
            active_grid = fde_active_grid
         else
            active_grid = fde_grid_im
         end if

         if (present(fde_active_gridfunction)) then
            active_gridfunction = fde_active_gridfunction
         else
            active_gridfunction = gf_frozen
         end if


         if (do_potential) then
            do_interaction_energy = .true.
         end if

      end if

!gosia debug
!      if (do_export_pertden_reorth) then
!        write(*, *) 'BLE parallel_fde = ', parallel_fde
!      end if


#ifdef VAR_MPI
      if (parallel_fde) then
!        broadcast non-allocatables
         call fde_mpi_bcast(mat_dim)
         call fde_mpi_bcast(nz)
         call fde_mpi_bcast(nr_dmat)
         call fde_mpi_bcast(nr_fmat)
         call fde_mpi_bcast(nr_atoms)
         call fde_mpi_bcast(do_potential)
         call fde_mpi_bcast(use_potential)
         call fde_mpi_bcast(do_interaction_energy)
         call fde_mpi_bcast(do_electrostatic_pot)
         call fde_mpi_bcast(do_london_rhs_direct_der1)
         call fde_mpi_bcast(do_london_rhs_direct_der2)
         call fde_mpi_bcast(do_london_rhs_reorth_coupling)
         call fde_mpi_bcast(do_london_rhs_direct_coupling)
         call fde_mpi_bcast(do_london_rhs_ro)
         call fde_mpi_bcast(do_london_susc2el_der1)
         call fde_mpi_bcast(do_london_susc2el_der2)
         call fde_mpi_bcast(do_london_susc2el_coupling)
         call fde_mpi_bcast(do_london_suscreo_coupling)
         call fde_mpi_bcast(do_export_pertden_reorth)
         call fde_mpi_bcast(do_export_pertden_direct)
         call fde_mpi_bcast(response_order_mo)
         call fde_mpi_bcast(response_order_ao)
         call fde_mpi_bcast(do_geo_0)
         call fde_mpi_bcast(do_geo_1)
         call fde_mpi_bcast(do_geo_2)
         call fde_mpi_bcast(do_london_ks)
         call fde_mpi_bcast(do_london_lr)
! we don't broadcast the grid, we only keep it in master
      end if
#endif

!     nullify(tmat_1)
!     nullify(tmat_2)
      if (do_geo_1) then
         allocate(tmat_1(mat_dim, mat_dim))
         tmat_1(:, :) = fde_dmat(:, :, 1, 1)
      end if

      if (do_geo_2) then
         allocate(tmat_1(mat_dim, mat_dim))
         allocate(tmat_2(mat_dim, mat_dim))
         tmat_1(:, :) = fde_dmat(:, :, 1, 1)
         tmat_2(:, :) = fde_dmat(:, :, 1, 2)
      end if

     nullify(dmat_0)
      if (.not. fde_mpi_is_master()) then
         allocate(dmat_0_container(mat_dim, mat_dim, nz))
         dmat_0 => dmat_0_container
      else
         dmat_0 => fde_dmat_0
      end if

      nullify(dmat)
      if (nr_dmat > 0) then
         if (.not. fde_mpi_is_master()) then
            allocate(dmat_container(mat_dim, mat_dim, nz, nr_dmat))
            dmat => dmat_container
         else
            if (present(fde_dmat)) then
               dmat => fde_dmat
            end if
         end if
      end if

     ! if (do_export_pertden_direct .and. .not. present(fde_dmat_pertden_direct)) then
     !    call quit('export_pertden_direct error')
     ! end if
      nullify(dmat_pertden_direct)
      if (do_export_pertden_direct) then
         if (.not. fde_mpi_is_master()) then
            allocate(dmat_pertden_direct_container(mat_dim, mat_dim, nz))
            dmat_pertden_direct => dmat_pertden_direct_container
         else
            !if (present(fde_dmat_pertden_direct)) then
               dmat_pertden_direct => fde_dmat_pertden_direct
            !end if
         end if
      end if

     ! if (do_export_pertden_reorth .and. .not. present(fde_dmat_pertden_reorth)) then
     !    call quit('export_pertden_reorth error')
     ! end if
      nullify(dmat_pertden_reorth)
      if (do_export_pertden_reorth) then
         if (.not. fde_mpi_is_master()) then
            allocate(dmat_pertden_reorth_container(mat_dim, mat_dim, nz, 3))
            dmat_pertden_reorth => dmat_pertden_reorth_container
         else
            !if (present(fde_dmat_pertden_reorth)) then
               dmat_pertden_reorth => fde_dmat_pertden_reorth
            !end if
         end if
      end if

      if (use_potential) then
! gosia: here fde_static_vemb seems not to be allocated, even if run with *EMBPOT in the input,
! so i have to read it in again:
         if (fde_mpi_is_master()) then
           call fde_import_static
           write(*, *) 'static potential imported'
         end if
      end if

      if (do_london_rhs_direct_coupling .or. do_london_rhs_reorth_coupling &
       .or. do_london_susc2el_coupling &
       .or. do_london_suscreo_coupling) then
         if (fde_mpi_is_master()) then
            if (.not. import_data_is_initialized) call fde_import_frozen
         end if
      end if

      nullify(fmat)
      if (nr_fmat > 0) then
         if (.not. fde_mpi_is_master()) then
            allocate(fmat_container(mat_dim, mat_dim, nz, nr_fmat))
            fmat_container = 0.0d0
            fmat => fmat_container
         else
            if (present(fde_fmat)) then
               fmat => fde_fmat
            end if
         end if
      end if

      if (nr_fmat > 0) then
         allocate(default_fmat_pg_sym(nr_fmat))
         default_fmat_pg_sym = 1
      end if
      if (nr_dmat > 0) then
         allocate(default_dmat_pg_sym(nr_dmat))
         allocate(default_dmat_ih_sym(nr_dmat))
         default_dmat_pg_sym = 1
         default_dmat_ih_sym = 0
      end if

      nullify(fmat_pg_sym)
      nullify(dmat_pg_sym)
      nullify(dmat_ih_sym)

      if (.not. fde_mpi_is_master()) then
         if (nr_fmat > 0) then
            allocate(fmat_pg_sym_container(nr_fmat))
            fmat_pg_sym => fmat_pg_sym_container
         end if
         if (nr_dmat > 0) then
            allocate(dmat_pg_sym_container(nr_dmat))
            allocate(dmat_ih_sym_container(nr_dmat))
            dmat_pg_sym => dmat_pg_sym_container
            dmat_ih_sym => dmat_ih_sym_container
         end if
      else
         if (present(fde_fmat_pg_sym)) then
            fmat_pg_sym => fde_fmat_pg_sym
         else
            fmat_pg_sym => default_fmat_pg_sym(1:nr_fmat)
         end if
         if (present(fde_dmat_pg_sym)) then
            dmat_pg_sym => fde_dmat_pg_sym
         else
            dmat_pg_sym => default_dmat_pg_sym(1:nr_dmat)
         end if
         if (present(fde_dmat_ih_sym)) then
            dmat_ih_sym => fde_dmat_ih_sym
         else
            dmat_ih_sym => default_dmat_ih_sym(1:nr_dmat)
         end if
      end if

      do_geo = .false.
      if (do_geo_0 .or. do_geo_1 .or. do_geo_2) then
         do_geo = .true.
      end if

      nullify(property_gradient)
      if (.not. fde_mpi_is_master()) then
         if (do_geo) then
            allocate(property_gradient_container(nr_atoms*3))
            property_gradient => property_gradient_container
            property_gradient = 0.0d0
         end if
      else
         if (present(fde_property_gradient)) then
            property_gradient => fde_property_gradient(1:nr_atoms*3)
            property_gradient = 0.0d0
         end if
      end if

      if (response_order_mo > max_response_order_mo) then
          stop 'response_order_mo > max_response_order_mo'
      end if

! functionals are set
      call fde_set_nadd_functionals(max_response_order_mo,parallel_fde) 

!     print *,'info, at process',fde_get_my_rank(),'fde_xcf is gga?',fun_is_gga(fde_xcf)
!     print *,'info, at process',fde_get_my_rank(),'fde_kef is gga?',fun_is_gga(fde_kef)


      use_gga_qr = (fun_is_gga(fde_xcf)) .or. (fun_is_gga(fde_kef))
      use_gga_qi = (fun_is_gga(fde_xcf)) .or. (fun_is_gga(fde_kef))
      if (response_order_mo > 0) then
         if (fde_cfg_alda_hs) use_gga_qr = .false.
         if (fde_cfg_alda_ha) use_gga_qi = .false.
      end if     
      if (use_gga_qr .or. use_gga_qi) then
         use_gga = .true.
      end if

      if (do_electrostatic_pot) then
         use_gga = .true.
      end if

!     in qr it is not so easy to separate h+ and h-
!     contributions since they may mix
!     therefore i quit here if the user wants qr
!     and not the full alda or xalda
!     it is really unfair to quit as late as here
!     this can be checked earlier
!     but this is anyway an exotic option
      if (response_order_mo > 1) then
         if (fde_cfg_alda_hs .and. (.not. fde_cfg_alda_ha)) then
            call fde_quit('qr and partial alda is not implemented')
         end if
         if (fde_cfg_alda_ha .and. (.not. fde_cfg_alda_hs)) then
            call fde_quit('qr and partial alda is not implemented')
         end if
      end if

      max_ao_g_order = 0
      if (use_gga) then
        max_ao_g_order = 1
      end if
      if (do_geo) then
        max_ao_g_order = max_ao_g_order + 1
      end if
     
      max_ao_m_order = 0
      if (do_london_rhs_direct_der1 .or. do_london_rhs_direct_der2 &
         .or. do_london_rhs_reorth_coupling .or. do_london_rhs_direct_coupling) then
        max_ao_m_order = 1
      end if
      if (do_london_susc2el_der2 .or. do_london_susc2el_coupling &
       .or. do_london_suscreo_coupling) then
        max_ao_m_order = 1
      end if
      if (do_london_susc2el_der1) then
        max_ao_m_order = 2
      end if

      call dirac_ao_eval_init(max_ao_g_order, &
                        max_ao_m_order, &
                        .true.)

      if (fde_mpi_is_master()) then
!        initial expectation value of the two-electron matrix
         if (do_potential) then
            save_ddot1 = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
         end if

         print *,'nr_mo_gerade_positive_active #1',nr_mo_gerade_positive_active
         print *,'nr_mo_ungerade_positive_active #1',nr_mo_ungerade_positive_active
!        get average density from open shells
         if ((nr_mo_gerade_positive_active + nr_mo_ungerade_positive_active) > 0) print *,'call gtdoav #1' 
         if ((nr_mo_gerade_positive_active + nr_mo_ungerade_positive_active) > 0) &
            call gtdoav(mat_dim, nz, 1, dmat_0, 1)

         if (.not. do_electrostatic_pot) then
#ifdef PRG_DIRAC
            call scale_density_matrices(2.0d0)
#endif

         !if (do_london_rhs_direct .or. do_london_rhs_ro) then
         !   call mat_lao_sym(mat_dim, nz, nr_fmat, fmat, fmat_pg_sym)
         !end if


         call insert_half_phases()
         end if ! do_electrostatic_pot

      end if ! fde_mpi_is_master

#ifdef VAR_MPI
      if (parallel_fde) then
!        broadcast allocatables
         call fde_mpi_bcast(dmat_0)
         if (nr_fmat > 0) then
!           never bcast fmat
!           only master's fmat contains "HF" contribution
!           otherwise it is multiplied
            call fde_mpi_bcast(fmat_pg_sym)
         end if
         if (nr_dmat > 0) then
            call fde_mpi_bcast(dmat)
            call fde_mpi_bcast(dmat_pg_sym)
            call fde_mpi_bcast(dmat_ih_sym)
         end if
         if (do_export_pertden_reorth) then
            call fde_mpi_bcast(dmat_pertden_reorth)
         end if
         if (do_export_pertden_direct) then
            call fde_mpi_bcast(dmat_pertden_direct)
         end if
      end if
#endif

      call interface_mo_read()

#ifdef DEBUG_XC
      write(*, *) 'debug: XC potential matrix zeroed out prior to integration'
      fmat = 0.0d0
#endif

      if (do_electrostatic_pot) then
         call loop_over_batches_elpot()
      else if (do_export_pertden_reorth .or. do_export_pertden_direct) then
         call loop_over_batches_pertden()
      else
         call loop_over_batches()
      end if
#ifdef VAR_MPI
      if (parallel_fde) then
!        collect results from processors
         call fde_mpi_reduce(nr_electrons_integrated)
         call fde_mpi_reduce(nr_points_used)
         if (do_potential) then
            call fde_mpi_reduce(fde_energy)
         end if
         if (do_interaction_energy) then
            call fde_mpi_reduce(non_additive_interaction)
            call fde_mpi_reduce(electrostatic_b_rho_a_interaction)
            call fde_mpi_reduce(electrostatic_a_rho_b_interaction)
         end if
         if (nr_fmat > 0) then
            call fde_mpi_reduce(fmat)
         end if
         if (do_geo) then
            call fde_mpi_reduce(property_gradient)
         end if
         if (do_london_susc2el_der1 .or. do_london_susc2el_der2 &
        .or. do_london_susc2el_coupling .or. do_london_suscreo_coupling) then
            call fde_mpi_reduce(susc2el_fde_integrated)
            call fde_mpi_reduce(susc2el_fde_integrated_sdft)
         end if
      end if
#endif


      if (fde_mpi_is_master()) then
       if (.not. do_electrostatic_pot) then
#ifdef PRG_DIRAC
         call scale_density_matrices(0.5d0)
#endif

         print *,'nr_mo_gerade_positive_active #2',nr_mo_gerade_positive_active
         print *,'nr_mo_ungerade_positive_active #2',nr_mo_ungerade_positive_active
!        get average density from open shells
         if ((nr_mo_gerade_positive_active + nr_mo_ungerade_positive_active) > 0) print *,'call gtdoav #2'
         if ((nr_mo_gerade_positive_active + nr_mo_ungerade_positive_active) > 0) call gtdoav(mat_dim, nz, 1, dmat_0, 1)

!        symmetrize matrices
         if (fde_use_blas3_scheme .or. (use_gga                         &
                                   .and. .not. do_london_rhs_direct_der1 &
                                   .and. .not. do_london_rhs_direct_der2 &
                                   .and. .not. do_london_rhs_ro &
                                   .and. .not. do_london_rhs_reorth_coupling &
                                   .and. .not. do_london_rhs_direct_coupling &
                                   .and. .not. do_london_susc2el_der1 &
                                   .and. .not. do_london_susc2el_der2 &
                                   .and. .not. do_london_susc2el_coupling &
                                   .and. .not. do_london_ks         &
                                   .and. .not. do_export_pertden_reorth &
                                   .and. .not. do_export_pertden_direct &
                                   .and. .not. do_geo)) then
!gosia: this check should be for nr_fmat, shouldn't it?
            if (nr_dmat > 0) then
               call gga_sym(mat_dim, nz, nr_fmat, fmat, fmat_pg_sym, dmat_ih_sym)
            else
               call gga_sym(mat_dim, nz, 1, fmat, (/1/), (/1/))
            end if
         end if

         if (do_london_rhs_direct_der1 .or. do_london_rhs_ro .or. &
             do_london_rhs_direct_der2 .or.  &
             do_london_rhs_reorth_coupling .or. do_london_rhs_direct_coupling) then
            call mat_lao_sym(mat_dim, nz, nr_fmat, fmat, fmat_pg_sym)
         end if

         call insert_half_phases()

!        final expectation value of the two-electron matrix
         if (do_potential) then
            save_ddot2    = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
            fde_mat_energy = save_ddot2 - save_ddot1
         end if

         call report()
       end if ! not do_electrostatic_pot
      end if

      call nullify_and_release()

#ifdef DEBUG_XC
      write(*, *) 'debug: real part of the XC potential matrix'
      call prqmat(fde_fmat,        &
                  mat_dim,        &
                  mat_dim,        &
                  mat_dim,        &
                  mat_dim,        &
                  1,              &
                  (/1, 2, 3, 4/), &
                  file_unit)
      call matrix_to_file('fde_fmat', mat_dim*mat_dim, fde_fmat)
#endif

!     check for nan in xc matrix
      do imat = 1, fde_nr_fmat
         do iz = 1, nz
            do i = 1, mat_dim
               do j = 1, mat_dim
                  if (fde_fmat(j, i, iz, imat) /= fde_fmat(j, i, iz, imat)) then
                     print *, 'error: found NaN in XC matrix, quitting'
                     stop
                  end if
               end do
            end do
         end do
      end do

   end subroutine


   subroutine find_vemb_ipoint(x, y, z, w, vemb)
      real(8), intent(in)  :: x, y, z, w
      !real(8), intent(out) :: vemb
      real(8), intent(in) :: vemb
      integer :: i
      integer :: nrow
      logical :: found

      do i = 1, fde_grid_sv%npoints
         found = (fde_grid_sv%r(1,i)==x).and.(fde_grid_sv%r(2,i)==y).and.(fde_grid_sv%r(3,i)==z).and.(fde_grid_sv%w(i)==w)
         if (found) then
            nrow = i
            if (vemb .ne. fde_static_vemb(nrow)) then
               write(*, *) 'wrong vemb'
            end if
         end if
      end do

!      vemb = fde_static_vemb(nrow)

   end subroutine

   subroutine find_pertden(x, y, z, w, n_b, gn_b, s_b, gs_b)
      real(8), intent(in) :: x, y, z, w
      real(8), intent(out), optional :: n_b(3)
      real(8), intent(out), optional :: gn_b(3,3)
      real(8), intent(out), optional :: s_b(3,3)
      real(8), intent(out), optional :: gs_b(3,3,3)
      integer :: i, j, k
      integer :: nrow, ncol
      logical :: found


!     first find the corresponding line on pertden:
      do i = 1, fde_pertden_nrow
         found = (pertden(1,i)==x).and.(pertden(2,i)==y).and.(pertden(3,i)==z).and.(pertden(4,i)==w)
         if (found) then
            nrow = i
         end if
      end do

!     gosia todo: add pointers indicating on which column of the file the given property starts
!     for now: we always count from the last column
!.....n_b:
      if (present(n_b)) then
         do i = 1, 3
            n_b(i) = pertden(fde_pertden_ncol-48+i, nrow)
         end do
      end if
!.....gn_b
      if (present(gn_b)) then
         do i = 1, 3
            do j = 1, 3
               ncol = (i-1)*3+j
               gn_b(j, i) = pertden(fde_pertden_ncol-45+ncol, nrow)
            end do
         end do
      end if
!.....s_b:
      if (present(s_b)) then
         do i = 1, 3
            do j = 1, 3
               ncol = (i-1)*3+j
               s_b(j, i) = pertden(fde_pertden_ncol-36+ncol, nrow)
            end do
         end do
      end if
!.....gs_b
      if (present(gs_b)) then
         do i = 1, 3
            do j = 1, 3
               do k = 1, 3
                  ncol = (i-1)*9+(j-1)*3+k
                  gs_b(k, j, i) = pertden(fde_pertden_ncol-27+ncol, nrow)
               end do
            end do
         end do
      end if

   end subroutine


!   ----------------------------------------------------------------------------
   subroutine nullify_and_release()
!   ----------------------------------------------------------------------------

      if (associated(dmat)) then
         nullify(dmat)
      end if
      if (associated(dmat_pertden_direct)) then
         nullify(dmat_pertden_direct)
      end if
      if (associated(dmat_pertden_reorth)) then
         nullify(dmat_pertden_reorth)
      end if
      if (associated(fmat)) then
         nullify(fmat)
      end if
      if (associated(pertden)) then
         nullify(pertden)
      end if
      if (associated(dmat_0)) then
         nullify(dmat_0)
      end if
      if (associated(dmat_pg_sym)) then
         nullify(dmat_pg_sym)
      end if
      if (associated(dmat_ih_sym)) then
         nullify(dmat_ih_sym)
      end if
      if (associated(fmat_pg_sym)) then
         nullify(fmat_pg_sym)
      end if
      if (associated(property_gradient)) then
         nullify(property_gradient)
      end if
      if (do_geo_1) then
         deallocate(tmat_1)
      end if
      if (do_geo_2) then
         deallocate(tmat_1)
         deallocate(tmat_2)
      end if
!     if (associated(tmat_1)) then
!        nullify(tmat_1)
!     end if
!     if (associated(tmat_2)) then
!        nullify(tmat_2)
!     end if

      if (.not. fde_mpi_is_master()) then
         if (allocated(dmat_container)) then
            deallocate(dmat_container)
         end if
         if (allocated(dmat_pertden_direct_container)) then
            deallocate(dmat_pertden_direct_container)
         end if
         if (allocated(dmat_pertden_reorth_container)) then
            deallocate(dmat_pertden_reorth_container)
         end if
         if (allocated(fmat_container)) then
            deallocate(fmat_container)
         end if
         if (allocated(dmat_0_container)) then
            deallocate(dmat_0_container)
         end if
         if (allocated(dmat_pg_sym_container)) then
            deallocate(dmat_pg_sym_container)
         end if
         if (allocated(dmat_ih_sym_container)) then
            deallocate(dmat_ih_sym_container)
         end if
         if (allocated(fmat_pg_sym_container)) then
            deallocate(fmat_pg_sym_container)
         end if
         if (allocated(property_gradient_container)) then
            deallocate(property_gradient_container)
         end if
      end if

      if (allocated(default_fmat_pg_sym)) then
         deallocate(default_fmat_pg_sym)
      end if
      if (allocated(default_dmat_pg_sym)) then
         deallocate(default_dmat_pg_sym)
      end if
      if (allocated(default_dmat_ih_sym)) then
         deallocate(default_dmat_ih_sym)
      end if

   end subroutine

   subroutine scale_density_matrices(f)

!     --------------------------------------------------------------------------
      real(8) :: f
!     --------------------------------------------------------------------------

      call dscal(mat_dim*mat_dim*nz, f, dmat_0, 1)
      if (nr_dmat > 0) then
        call dscal(mat_dim*mat_dim*nz*nr_dmat, f, dmat, 1)
      end if
!gosia todo: check if these matrices should be scaled too
! (now they are made to agree with dft)
! insert_half_phases should be then adapted as well
      if (do_export_pertden_reorth) then
        call dscal(mat_dim*mat_dim*nz*3, f, dmat_pertden_reorth, 1)
      end if
      if (do_export_pertden_direct) then
        call dscal(mat_dim*mat_dim*nz, f, dmat_pertden_direct, 1)
      end if


   end subroutine

   subroutine insert_half_phases()

!     --------------------------------------------------------------------------
      integer :: iz, im, iq
!     --------------------------------------------------------------------------

      if (nz < 4) then
         if (nr_dmat > 0) then
            do iz = 1, nz
               do im = 1, nr_dmat
                  iq = pq_to_uq(iz, dmat_pg_sym(im) - 1)
                  call my_q2bphase('D', iq, 1, dmat(1, 1, iz, im))
               end do
            end do
         end if
         if (nr_fmat > 0) then
            do iz = 1, nz
               do im = 1, nr_fmat
                  iq = pq_to_uq(iz, fmat_pg_sym(im) - 1)
                  call my_q2bphase('F', iq, 1, fmat(1, 1, iz, im))
               end do
            end do
         end if
!gosia try:
         if (do_export_pertden_reorth) then
            do iz = 1, nz
               do im = 1, 3
                  iq = pq_to_uq(iz, dmat_pg_sym(im) - 1)
                  call my_q2bphase('D', iq, 1, dmat_pertden_reorth(1, 1, iz, im))
               end do
            end do
         end if
         if (do_export_pertden_direct) then
            do iz = 1, nz
               iq = pq_to_uq(iz, dmat_pg_sym(1) - 1)
               call my_q2bphase('D', iq, 1, dmat_pertden_direct(1, 1, iz))
            end do
         end if

      end if

   end subroutine


! ------------------------------------------------------------------------------
   subroutine integrate()
! ------------------------------------------------------------------------------
!        unperturbed density
   real(8) :: n_0(max_block_length)
!        unperturbed density gradient
   real(8) :: gn_0_pointwise(3)

   real(8) :: gnn_0(max_block_length)
   real(8) :: n_0_t(max_block_length)
   real(8) :: gnn_0_t(max_block_length)
!        --------------------------------------------------------------------------
   real(8) :: n_0_f(max_block_length)
   real(8) :: gnn_0_f(max_block_length)
   real(8) :: w(max_block_length)
!        --------------------------------------------------------------------------
!     geometric business
   real(8) :: r_t1(max_block_length)
   real(8) :: r_t2(max_block_length)
!     --------------------------------------------------------------------------
!     fde and fde response
!      logical :: fde_resp_verbose =.false.
   real(8) :: z0tt(3)
!     --------------------------------------------------------------------------
   integer :: id, m
   integer :: ii, ji
   real(8) :: temp
   real(8) :: block_threshold
   real(8) :: rho13, rho_outer, rho_resp, f
   integer :: ipoint
   integer :: k
   real(8) :: VT(3)
   integer              :: derv_length, im, i
   real(8), allocatable :: derv_active(:, :)
   real(8), allocatable :: derv_frozen(:, :)
   real(8), allocatable :: derv_total(:, :)

   logical :: print_data_at_point

   real(8) :: n_b(3)
   real(8) :: gnn_b(3)
   real(8) :: gn_b(3, 3)
   real(8) :: n_bb(3, 3)
   real(8) :: n_bb_total(3, 3), n_bb_active(3, 3)
   real(8) :: t, gt(3)

   real(8) :: vemb

   logical :: un_above_threshold
   logical :: us_above_threshold
   logical :: ugn_above_threshold
   logical :: ugs_above_threshold
   logical :: utn_above_threshold
   logical :: uts_above_threshold
   logical :: un_b_lao_above_threshold
   logical :: ugn_b_lao_above_threshold
   logical :: us_b_lao_above_threshold
   logical :: ugs_b_lao_above_threshold

   logical :: nuc_part
   logical :: ele_part
   real(8), allocatable :: active_nucpot(:)
   integer :: np_esp
   integer :: irrep_esp
   integer :: ixyz, jxyz

   real(8), allocatable :: ao(:, :)
   real(8), allocatable :: buffer(:, :)

   type(fde_omega_prefactor) :: u


   if (nr_fmat > 0) then
      allocate(u%n (      nr_fmat))
      allocate(u%gn(   3, nr_fmat))
      allocate(u%s (   3, nr_fmat))
      allocate(u%gs(3, 3, nr_fmat))
      allocate(u%tn(      nr_fmat))
      allocate(u%ts(   3, nr_fmat))
      allocate(u%n_b_lao(nr_fmat)) 
      allocate(u%gn_b_lao(3, nr_fmat)) 
      allocate(u%s_b_lao(3, nr_fmat)) 
      allocate(u%gs_b_lao(3, 3, nr_fmat))
   end if

!  if (fde_use_blas3_scheme) then
!      block_threshold = get_block_threshold(fde_cfg_tinydens, mat_dim, dmat_0)
!   end if

!   check what functional derivatives do we need
!   ============================================

   max_fun_derv = 0

   if (do_potential) then
      max_fun_derv = 1
   end if
   if (do_interaction_energy.and..not.do_potential) then
      max_fun_derv = 0
   end if
   if (response_order_mo > 0) then
      max_fun_derv = response_order_mo + 1
   end if
   if (do_london_rhs_direct_der1) then
      max_fun_derv = 1
   end if
   if (do_london_susc2el_der1) then
      max_fun_derv = 1
   end if
   if (do_london_susc2el_der2 .or. do_london_susc2el_coupling .or. do_london_suscreo_coupling) then
      max_fun_derv = 2
   end if
   if (do_london_rhs_direct_der2 .and. .not. fde_cfg_no_sdft) then
      max_fun_derv = 2
   end if
   if (do_london_rhs_reorth_coupling .or. do_london_rhs_direct_coupling) then
      max_fun_derv = 2
   end if
   if (do_london_rhs_ro) then
      if (max_fun_derv < 2) then
         max_fun_derv = 2
      end if
   end if
   if (do_geo_0) then
    max_fun_derv = 1
   end if
   if (do_geo_1) then
    max_fun_derv = 2
   end if
   if (do_geo_2) then
    max_fun_derv = 3
   end if
   if (do_london_ks) then
    max_fun_derv = 1
   end if
   if (do_london_lr) then
    max_fun_derv = 2
   end if

   derv_length = nr_nonzero_derv
   allocate(derv_active(max_block_length, 0:derv_length))
   allocate(derv_total(max_block_length, 0:derv_length))

!  better safe than sorry
   derv_active = 0.0d0
   derv_total  = 0.0d0

   if (do_interaction_energy) then
      allocate(derv_frozen(max_block_length, 0:derv_length))
      derv_frozen = 0.0d0
   end if

!  set it to zero, otherwise it is undefined with lda
   gnn_0 = 0.0d0

  if (fde_use_blas3_scheme) then
      loop_over_points_block: do ipoint = 1, nr_points_on_this_proc, max_block_length
         block_length = min(max_block_length, nr_points_on_this_proc - ipoint + 1)

         allocate(ao(block_length,     nr_ao_slices*nr_ao_cartesian))
         ao = 0.0d0
         allocate(buffer(block_length, nr_ao_slices*nr_ao_cartesian))
         buffer = 0.0d0

         call get_ao(                                  &
                     block_length,                     &
                     rx(ipoint:ipoint+block_length-1), &
                     ry(ipoint:ipoint+block_length-1), &
                     rz(ipoint:ipoint+block_length-1), &
                     ao,                               &
                     buffer                            &
                    )

!        call fde_contribution_blas3(                          &
!                                   use_gga,                  &
!                                   use_gga_qr,               &
!                                   use_gga_qi,               &
!                                   (response_order_mo == 1), &
!                                   block_length,             &
!                                   block_threshold,          &
!                                   mat_dim,                  &
!                                   nz,                       &
!                                   buffer,                   &
!                                   ao,                       &
!                                   nr_fmat,                  &
!                                   fmat,                     &
!                                   dmat_0,                   &
!                                   dmat_saop_outer,          &
!                                   dmat_saop_resp,           &
!                                   nr_dmat,                  &
!                                   dmat,                     &
!                                   dmat_pg_sym,              &
!                                   dmat_ih_sym,              &
!                                   derv_length,              &
!                                   derv,                     &
!                                   max_fun_derv,             &
!                                   rw(ipoint),               &
!                                   nr_electrons_integrated,  &
!                                   fde_energy                &
!                                  )

         deallocate(ao)
         deallocate(buffer)

      end do loop_over_points_block
   end if 

   if (.not. fde_use_blas3_scheme) then

         allocate(ao(1,     nr_ao_slices*nr_ao_cartesian))
         ao = 0.0d0
         allocate(buffer(1, nr_ao_slices*nr_ao_cartesian))
         buffer = 0.0d0

         print_data_at_point = .true.

         np_esp = 1
         irrep_esp = 1
         if (do_interaction_energy) then
            allocate(active_nucpot(np_esp))
            active_nucpot = 0.0d0
         end if

!  here starts the expensive loop over points
!  do not evaluate anything inside this loop
!  that does not change from point to point
   do ipoint = 1, nr_points_on_this_proc

      time_ao_start = second()
      call get_ao(1,    &
            rx(ipoint), &
            ry(ipoint), &
            rz(ipoint), &
            ao,         &
            buffer)
      time_ao = time_ao + second() - time_ao_start

      time_density_start = second()
      if (need_ao_order(1, 0)) then
!        density and density gradient
         call get_gn(n_0(1), gn_0_pointwise, 0, mat_dim, dmat_0, buffer, ao)
         gnn_0(1) = gn_0_pointwise(1)*gn_0_pointwise(1) &
                  + gn_0_pointwise(2)*gn_0_pointwise(2) &
                  + gn_0_pointwise(3)*gn_0_pointwise(3)
      else
!        density
         call get_n(n_0(1), 0, mat_dim, dmat_0, buffer, ao)
      endif

      time_density = time_density + second() - time_density_start

      if (n_0(1) > fde_cfg_tinydens) then
         nr_points_used = nr_points_used + 1

         if (allocated(u%n)) u%n  = 0.0d0
         if (allocated(u%gn)) u%gn = 0.0d0
         if (allocated(u%s)) u%s  = 0.0d0
         if (allocated(u%gs)) u%gs = 0.0d0
         if (allocated(u%tn)) u%tn = 0.0d0
         if (allocated(u%ts)) u%ts = 0.0d0
         if (allocated(u%n_b_lao)) u%n_b_lao  = 0.0d0
         if (allocated(u%gn_b_lao)) u%gn_b_lao  = 0.0d0
         if (allocated(u%s_b_lao)) u%s_b_lao  = 0.0d0
         if (allocated(u%gs_b_lao)) u%gs_b_lao  = 0.0d0

         w(1) = rw(ipoint)


         if ( frozen_n(ipoint) > 0.0d0 ) then
            n_0_t(1) = n_0(1)            + frozen_n(ipoint)
            z0tt(1)  = gn_0_pointwise(1) + frozen_gnx(ipoint)
            z0tt(2)  = gn_0_pointwise(2) + frozen_gny(ipoint)
            z0tt(3)  = gn_0_pointwise(3) + frozen_gnz(ipoint)

            gnn_0_t(1) = z0tt(1)*z0tt(1) + z0tt(2)*z0tt(2) + z0tt(3)*z0tt(3)

            n_0_f(1)   = frozen_n(ipoint)
            gnn_0_f(1) = frozen_gnx(ipoint)*frozen_gnx(ipoint) &
                       + frozen_gny(ipoint)*frozen_gny(ipoint) &
                       + frozen_gnz(ipoint)*frozen_gnz(ipoint) 
         else
            n_0_t(1)   = n_0(1)
            gnn_0_t(1) = gnn_0(1)

            n_0_f(1)   = 0.0d0
            gnn_0_f(1) = 0.0d0
         endif
!       number of electrons
         nr_electrons_integrated = nr_electrons_integrated + w(1)*n_0(1)

         time_derv_start = second()

! aspg: this way we are forcing the same xc functional for the active and non-additive
!       contributions. this is not desired, but we have to overcome xcfun's limitation
!       on using only up to 5 functionals...
!
         call get_xcke_fun_derv(max_fun_derv, & 
                                bllen,        &
                                w,            &
                                n_0_t,        &
                                gnn_0_t,      &
                                derv_total)

         call get_xcke_fun_derv(max_fun_derv, & 
                                bllen,        &
                                w,            &
                                n_0,          &
                                gnn_0,        &
                                derv_active)

!gosia: we don't need derv_active for the coupling-kernel term, so save time.
         if (.not. do_london_rhs_reorth_coupling .and. .not. do_london_rhs_direct_coupling) then
           call get_xcke_fun_derv(max_fun_derv, & 
                                  bllen,        &
                                  w,            &
                                  n_0,          &
                                  gnn_0,        &
                                  derv_active)
         end if

         time_derv = time_derv + second() - time_derv_start


         time_matrix_dist_start = second()

!       interaction energy components
!       ==============================
        if (do_interaction_energy) then
            
            call get_xcke_fun_derv(max_fun_derv, & 
                                   bllen,        &
                                   w,            &
                                   n_0_f,        &
                                   gnn_0_f,      &
                                   derv_frozen)

            nuc_part = .true.
            ele_part = .false.
            call get_esp(np_esp,        &
                         active_nucpot, &
                         irrep_esp,     &
                         size(dmat, 1), &
                         dmat,          &
                         (/ rx(ipoint), ry(ipoint), rz(ipoint) /), &
                         nuc_part,      &
                         ele_part)

! aspg, 14/10/2015
!       note that the electrostatic, hartree and nuclear potentials in dirac have
!       been multiplied by -1 when exported so that we don't have to multipy the 
!       product potential * density by the electron charge when calculating the
!       energy. with this, the potentials will be in line with those obtained
!       e.g. with adf 

            active_nucpot(1) = -active_nucpot(1)
!
            non_additive_interaction = non_additive_interaction &
                              + (derv_total(1, d0000000) - derv_active(1, d0000000) - derv_frozen(1, d0000000))
            electrostatic_b_rho_a_interaction = electrostatic_b_rho_a_interaction + w(1)*frozen_elpot(ipoint)*n_0(1)
            electrostatic_a_rho_b_interaction = electrostatic_a_rho_b_interaction + w(1)*active_nucpot(1)*n_0_f(1)

        end if ! interaction energy components

!       embedding potential
!       ===================
         if (do_potential) then

            u%n(1) = (derv_total(1, d1000000) - derv_active(1, d1000000)) + w(1)*frozen_elpot(ipoint) 

            if (use_gga) then
!              this factor has nothing to do
!              with symmetrization later
!              it comes from the fact that derivative of z
!              yields two identical terms (one left, one right)
              u%gn(:, 1) = 2.0d0*( derv_total(1, d0010000)*z0tt - derv_active(1, d0010000)*gn_0_pointwise)
            end if

         end if !embedding potential

         if (use_potential) then
            vemb = use_vemb(ipoint)*w(1)
!left for debugging, normally the order of embedding potential should follow the order of grid points, so this is not necessary...
!            call find_vemb_ipoint(rx(ipoint), ry(ipoint), rz(ipoint), w(1), use_vemb(ipoint))
         end if

!       london contributions to differentiated Kohn-Sham matrix, 
!       added to property gradient, closed-shell only

!       (1) direct-lao contributions, 1st functional derivatives
!       --------------------------------------------------------
         if (do_london_rhs_direct_der1) then
             if (use_potential) then
!               use external embedding potential
                call lda_london_distribute_r(vemb,                              &
                                             fmat_pg_sym,                       &
                                             mat_dim,                           &
                                             nz,                                &
                                             fmat,                              &
                                             ao,                                &
                                             rx(ipoint), ry(ipoint), rz(ipoint))
             else
                u%n(1) = (derv_total(1, d1000000) - derv_active(1, d1000000)) + w(1)*frozen_elpot(ipoint)
                if (use_gga) then
                  u%gn(:, 1) = 2.0d0*( derv_total(1, d0010000)*z0tt - derv_active(1, d0010000)*gn_0_pointwise)
                
                  call gga_london_distribute_r(u%n(1),                            &
                                               1.0d0,                             &
                                               u%gn(1:3, 1),                      &
                                               fmat_pg_sym,                       &
                                               mat_dim,                           &
                                               nz,                                &
                                               fmat,                              &
                                               ao, &
                                               rx(ipoint), ry(ipoint), rz(ipoint))
                else
                  call lda_london_distribute_r(u%n(1),                            &
                                               fmat_pg_sym,                       &
                                               mat_dim,                           &
                                               nz,                                &
                                               fmat,                              &
                                               ao,                                &
                                               rx(ipoint), ry(ipoint), rz(ipoint))
                end if

             end if
         
         end if

!       (2) london reorthogonalization contribution, 2nd functional derivatives
!       -----------------------------------------------------------------------
           if (do_london_rhs_ro) then

              call lao_reorth(bllen,          &
                              mat_dim,        &
                              nz,             &
                              ao,             &
                              n_0,            &
                              n_0_t,          &
                              gn_0_pointwise, &
                              z0tt,           &
                              gnn_0,          &
                              gnn_0_t,        &
                              w,              &
                              use_gga,        &
                              dmat,           &
                              dmat_pg_sym,    &
                              buffer,         &
                              derv_length,    &
                              derv_active,    &
                              derv_total,     &
                              u)
           end if

!       (3) direct-lao spin density contribution, 2nd functional derivatives
!       --------------------------------------------------------------------

           if (do_london_rhs_direct_der2) then

              call london_direct_kern(bllen,          &
                                      mat_dim,        &
                                      nz,             &
                                      ao,             &
                                      dmat_0,         &
                                      buffer,         &
                                      n_0,            &
                                      n_0_t,          &
                                      gn_0_pointwise, &
                                      z0tt,           &
                                      use_gga,        &
                                      derv_length,    &
                                      derv_active,    &
                                      derv_total,     &
                                      rx(ipoint), ry(ipoint), rz(ipoint), &
                                      u)
           end if

           if (do_london_rhs_reorth_coupling .or. do_london_rhs_direct_coupling) then
              call london_kern_coupling(ipoint,      &
                                        bllen,       &
                                        mat_dim,     &
                                        nz,          &
                                        ao,          &
                                        n_0_t,       &
                                        z0tt,        &
                                        gnn_0_t,     &
                                        w,           &
                                        use_gga,     &
                                        derv_length, &
                                        derv_total,  &
                                        rx(ipoint), ry(ipoint), rz(ipoint), &
                                        u)
           end if

!       linear response
!       ===============

        if (response_order_mo == 1) then
           call lr_fde_uncoupled(  &
                   bllen,          &
                   mat_dim,        &
                   nz,             &
                   nr_dmat,        &
                   dmat,           &
                   dmat_pg_sym,    &
                   dmat_ih_sym,    &
                   .false.,        &
                   ao,             &
                   n_0,            &
                   n_0_t,          &
                   gn_0_pointwise, &
                   z0tt,           &
                   gnn_0,          &
                   gnn_0_t,        &
                   w,              &
                   use_gga_qr,     &
                   use_gga_qi,     &
                   buffer,         &
                   derv_length,    &
                   derv_active,    &
                   derv_total,     &
                   u)
      endif

#ifdef FDE_PLACEHOLDER
!       quadratic response
!       ==================
! N/A
! what goes below is strictly what the standard dft integrator has, so here it's only a placeholder

      if (response_order_mo == 2) then

         call sdft_qr(bllen,           &
                    use_gga,                &
                    fde_cfg_no_sdft,        &
                    mat_dim, dmat,        &
                    dmat_pg_sym,dmat_ih_sym,ao,&
                    n_0(1),gnn_0(1),gn_0_pointwise,w(1), &
                    buffer,             &
               derv_length,    &
               derv, u)

!         the idea here is that   r^[B]               \omega_{ia}^[C]
!                                 r^[C]               \omega_{ia}^[B]
!                                 r^[BC]              \omega_{ia}
!         and the corresponding   q^[B]  \cdot \nabla \omega_{ia}^[C]
!                                 q^[C]  \cdot \nabla \omega_{ia}^[B]
!                                 q^[BC] \cdot \nabla \omega_{ia}

!         can be evaluated using lr code (below)
!         only the remaining part is done in sdft_qr (above)
!
         call lr_fde_uncoupled( &
                bllen,          &
                mat_dim,        &
                nz,             &
                nr_dmat,        &
                dmat,           &
                dmat_pg_sym,    &
                dmat_ih_sym,    &
                .false.,        &
                ao,             &
                n_0,            &
                n_0_t,          &
                gn_0_pointwise, &
                z0tt,           &
                gnn_0,          &
                gnn_0_t,        &
                w,              &
                use_gga_qr,     &
                use_gga_qi,     &
                buffer,         &
                derv_length,    &
                derv_active,    &
                derv_total,     &
                u)
 
      endif

!
! all the geo, london stuff below are placeholders, roughly taken from the dft integrator
!
      if (do_geo_0) then
         if (use_gga) then
            call xc_geo_0_gga_c1(property_gradient, &
                                 derv_active(1, d1000000), &
                                 derv_active(1, d0010000), &
                                 derv_total(1, d1000000), &
                                 derv_total(1, d0010000), &
                                 gn_0_pointwise,    &
                                 mat_dim,           &
                                 dmat_0,            &
                                 ao)
         else
            call xc_geo_0_lda_c1(property_gradient, &
                                 derv_active(1, d1000000), &
                                 derv_total(1, d1000000), &
                                 mat_dim,           &
                                 dmat_0,            &
                                 ao)
         end if
      end if

      if (do_geo_1) then
         call get_n(r_t1(1), 0, mat_dim, tmat_1, buffer, ao)
         call xc_geo_1_lda_c1(property_gradient, &
                              derv_active(1, d1000000), &
                              derv_active(1, d2000000), &
                              derv_total(1, d1000000), &
                              derv_total(1, d2000000), &
                              r_t1(1),           &
                              mat_dim,           &
                              dmat_0,            &
                              tmat_1,            &
                              ao)
      end if

      if (do_geo_2) then
       call get_n(r_t1(1), 0, mat_dim, tmat_1, buffer, ao)
       call get_n(r_t2(1), 0, mat_dim, tmat_2, buffer, ao)
       call xc_geo_2_lda_c1(property_gradient, &
                            derv_active(1, d2000000), &
                            derv_active(1, d3000000), &
                            derv_total(1, d2000000),  &
                            derv_total(1, d3000000),  &
                            r_t1(1),           &
                            r_t2(1),           &
                            mat_dim,           &
                            dmat_0,            &
                            tmat_1,            &
                            tmat_2,            &
                            ao)
      end if

      if (do_london_ks) then
       if (use_gga) then
          call xc_london_c1_ks_gga( &
                                   derv_active(1, d1000000),        &
                                   derv_active(1, d0010000),        &
                                   derv_total(1, d1000000),         &
                                   derv_total(1, d0010000),         &
                                   gn_0_pointwise,                  &
                                   mat_dim,                         &
                                   nz,                              &
                                   fmat,                            &
                                   ao,                              &
                                   (/rx(ipoint), ry(ipoint), rz(ipoint)/))
       else
          call xc_london_c1_ks_lda( &
                                    derv_active(1, d1000000), &
                                    derv_total(1, d1000000), &
                                   mat_dim,           &
                                   nz,                &
                                   fmat,              &
                                   ao,                &
                                   (/rx(ipoint), ry(ipoint), rz(ipoint)/))
       end if
      end if

      if (do_london_lr) then
       if (use_gga) then
          call xc_london_c1_lr_gga_densities(mat_dim,                         &
                                             nz,                              &
                                             dmat,                            &
                                             ao,                              &
                                           (/rx(ipoint),                      &
                                             ry(ipoint),                      &
                                             rz(ipoint)/),                    &
                                             gn_0_pointwise,                            &
                                             n_b,                             &
                                             gn_b,                            &
                                             gnn_b)
       else
          call xc_london_c1_lr_lda_densities(mat_dim,      &
                                             nz,           &
                                             dmat,         &
                                             ao,           &
                                           (/rx(ipoint),   &
                                             ry(ipoint),   &
                                             rz(ipoint)/), &
                                             n_b)
!radovan:    a short explanation of what is happening here
!            (because this is done differently than in Dalton):
!            the above routines evaluate magnetic field derivatives
!            (there are three: im = 1, 3)
!            of perturbed density variables using the transition density matrix dmat
!            next we form u%n and u%gn and distribute over matrices
!            as if it was a normal linear response run
!            this makes the otherwise very lengthy operation rather compact
             do im = 1, 3
                u%n(im) = derv(1, d2000000)*n_b(im)
                if (use_gga) then
                   u%n(im) = u%n(im) + derv(1, d1010000)*gnn_b(im)
                   u%gn(1:3, im) = 2.0d0*gn_0_pointwise(1:3)*n_b(im)*derv(1, d1010000)   &
                                 + 2.0d0*gn_0_pointwise(1:3)*gnn_b(im)*derv(1, d0020000) &
                                 + 2.0d0*gn_b(1:3, im)*derv(1, d0010000)
                end if
             end do
          end if

#endif


!         distribute over fock matrices
!         =============================

          do im = 1, nr_fmat

             un_above_threshold = .false.
             if (dabs(u%n(im)) > tiny(0.0d0)) then
                un_above_threshold = .true.
             end if

             ugn_above_threshold = .false.
             if (maxval((/dabs(u%gn(1, im)), &
                          dabs(u%gn(2, im)), &
                          dabs(u%gn(3, im))/)) > tiny(0.0d0)) then
                ugn_above_threshold = .true.
             end if

             if (ugn_above_threshold) then
!               F_pq +=   u%n    \chi_p* \chi_q
!                     + 2 u%gn_i \chi_p* \nabla_i \chi_q
                call nabla_omega_real(fmat,                &
                                      im,                  &
                                      fmat_pg_sym(im) - 1, &
                                      mat_dim,             &
                                      nz,                  &
                                      u%n(im),             &
                                      u%gn(1, im),         &
                                      ao)
             else
                if (un_above_threshold) then
!                  F_pq += u%n \chi_p* \chi_q
                   call omega_real(fmat,                &
                                   im,                  &
                                   fmat_pg_sym(im) - 1, &
                                   mat_dim,             &
                                   nz,                  &
                                   u%n(im),             &
                                   ao)
                end if
             end if

! spin-density stuff 

             us_above_threshold = .false.
             if (maxval((/dabs(u%s(1, im)), &
                          dabs(u%s(2, im)), &
                          dabs(u%s(3, im))/)) > tiny(0.0d0)) then
                us_above_threshold = .true.
             end if

             ugs_above_threshold = .false.
             if (maxval((/dabs(u%gs(1, 1, im)), &
                          dabs(u%gs(1, 2, im)), &
                          dabs(u%gs(1, 3, im)), &
                          dabs(u%gs(2, 1, im)), &
                          dabs(u%gs(2, 2, im)), &
                          dabs(u%gs(2, 3, im)), &
                          dabs(u%gs(3, 1, im)), &
                          dabs(u%gs(3, 2, im)), &
                          dabs(u%gs(3, 3, im))/)) > tiny(0.0d0)) then
                ugs_above_threshold = .true.
             end if

             if (ugs_above_threshold) then
!               F_pq +=   u%s_j   \chi_p* \Sigma_j \chi_q
!                     + 2 u%gs_ij \chi_p* \Sigma_j \nabla_i \chi_q
                call nabla_omega_imag(fmat,                &
                                      im,                  &
                                      fmat_pg_sym(im) - 1, &
                                      mat_dim,             &
                                      nz,                  &
                                      u%s(1, im),          &
                                      u%gs(1, 1, im),      &
                                      ao)
             else
                if (us_above_threshold) then
!                  F_pq += u%s_j \chi_p* \Sigma_j \chi_q
                   call omega_imag(fmat,                &
                                   im,                  &
                                   fmat_pg_sym(im) - 1, &
                                   mat_dim,             &
                                   nz,                  &
                                   u%s(1, im),          &
                                   ao)
                end if
             end if

!            london contributions to property gradient
!            gosia: this "if" should not be needed here
!            nevertheless without that the test "fde_response_autodiff" fails on gfortran
             if (do_london_rhs_direct_der2.or.do_london_rhs_reorth_coupling &
                 .or. do_london_rhs_direct_coupling) then

               us_b_lao_above_threshold = .false.
               if (maxval((/dabs(u%s_b_lao(1, im)), &
                            dabs(u%s_b_lao(2, im)), &
                            dabs(u%s_b_lao(3, im))/)) > tiny(0.0d0)) then
                  us_b_lao_above_threshold = .true.
               end if
               
               ugs_b_lao_above_threshold = .false.
               if (maxval((/dabs(u%gs_b_lao(1, 1, im)), &
                            dabs(u%gs_b_lao(1, 2, im)), &
                            dabs(u%gs_b_lao(1, 3, im)), &
                            dabs(u%gs_b_lao(2, 1, im)), &
                            dabs(u%gs_b_lao(2, 2, im)), &
                            dabs(u%gs_b_lao(2, 3, im)), &
                            dabs(u%gs_b_lao(3, 1, im)), &
                            dabs(u%gs_b_lao(3, 2, im)), &
                            dabs(u%gs_b_lao(3, 3, im))/)) > tiny(0.0d0)) then
                  ugs_b_lao_above_threshold = .true.
               end if
               
               if (ugs_b_lao_above_threshold) then
!                 F_pq +=   u%s_j   \chi_p* \Sigma_j \chi_q
!                       + 2 u%gs_ij \chi_p* \Sigma_j \nabla_i \chi_q
                  call nabla_omega_imag(fmat,                &
                                        im,                  &
                                        fmat_pg_sym(im) - 1, &
                                        mat_dim,             &
                                        nz,                  &
                                        u%s_b_lao(1, im),          &
                                        u%gs_b_lao(1, 1, im),      &
                                        ao)
               else
                  if (us_b_lao_above_threshold) then
!                    F_pq += u%s_j \chi_p* \Sigma_j \chi_q
                     call omega_imag(fmat,                &
                                     im,                  &
                                     fmat_pg_sym(im) - 1, &
                                     mat_dim,             &
                                     nz,                  &
                                     u%s_b_lao(1, im),          &
                                     ao)
                  end if
               end if

             end if !if (do_london_rhs_direct_der2.or.do_london_rhs_kern_coupling) then

             if (do_london_rhs_reorth_coupling) then

               un_b_lao_above_threshold = .false.
               if (maxval((/dabs(u%n_b_lao(im)), &
                            dabs(u%n_b_lao(im)), &
                            dabs(u%n_b_lao(im))/)) > tiny(0.0d0)) then
                  un_b_lao_above_threshold = .true.
               end if
               
               ugn_b_lao_above_threshold = .false.
               if (maxval((/dabs(u%gn_b_lao(1, im)), &
                            dabs(u%gn_b_lao(2, im)), &
                            dabs(u%gn_b_lao(3, im))/)) > tiny(0.0d0)) then
                  ugn_b_lao_above_threshold = .true.
               end if
               
               if (ugn_b_lao_above_threshold) then
!                 F_pq +=   u%s_j   \chi_p* \Sigma_j \chi_q
!                       + 2 u%gs_ij \chi_p* \Sigma_j \nabla_i \chi_q
                  call nabla_omega_real(fmat,                &
                                        im,                  &
                                        fmat_pg_sym(im) - 1, &
                                        mat_dim,             &
                                        nz,                  &
                                        u%n_b_lao(im),          &
                                        u%gn_b_lao(1, im),      &
                                        ao)
               else
                  if (un_b_lao_above_threshold) then
!                    F_pq += u%s_j \chi_p* \Sigma_j \chi_q
                     call omega_real(fmat,                &
                                     im,                  &
                                     fmat_pg_sym(im) - 1, &
                                     mat_dim,             &
                                     nz,                  &
                                     u%n_b_lao(im),          &
                                     ao)
                  end if
               end if

             end if !if (do_london_rhs_kern_coupling) then


          end do ! distribute over Fock matrices

! end spin-density stuff
          time_matrix_dist = time_matrix_dist + second() - time_matrix_dist_start

!         fde contributions to expectation value part of magnetizability
!         ==============================================================
!         gosia: closed-shell only, called from pammag.F (susc2el subroutine)
!         these are 'pure-direct' terms
!         now works only in C1 symmetry, because isymop in pammag is hardcoded to 1
!         *c18 subroutines use overlap derivatives and not AO derivs 
!         (please leave it to compare directly with dalton)

          if (do_london_susc2el_der1) then
            n_bb = 0.0d0
!           embpot-dependent term, no SDFT contributions
            t= (derv_total(1, d1000000) - derv_active(1, d1000000)) + w(1)*frozen_elpot(ipoint)
            if (use_gga) then
              gt(:) = 2.0d0*( derv_total(1, d0010000)*z0tt - derv_active(1, d0010000)*gn_0_pointwise)
              call gga_london_sus2el_der1(t,        &
                                          1.0d0,    &
                                          gt,       &
                                          (/1,1,1/),&
                                          mat_dim,  &
                                          nz,       &
                                          dmat_0,   &
                                          ao,       &
                                          rx(ipoint), ry(ipoint), rz(ipoint), &
                                          buffer, &
                                          n_bb)
              !call xc_london_c1_susc_gga_der1(derv(1, d1000000),   &
              !                            2.0d0*derv(1, d0010000),  &
              !                            gn_0_pointwise,           &
              !                            mat_dim,             &
              !                            nz,                  &
              !                            dmat_0,                &
              !                            ao,                  &
              !                            (/rx(ipoint), ry(ipoint), rz(ipoint)/), &
              !                            n_bb)
            else
              call lda_london_sus2el_der1(t,         &
                                          (/1,1,1/), &
                                          mat_dim,   &
                                          nz,        &
                                          dmat_0,    &
                                          ao,        &
                                          rx(ipoint), ry(ipoint), rz(ipoint), &
                                          buffer, &
                                          n_bb)
              !call xc_london_c1_susc_der1(t,                   &
              !                            mat_dim,             &
              !                            nz,                  &
              !                            dmat_0,                &
              !                            ao,                  &
              !                            (/rx(ipoint), ry(ipoint), rz(ipoint)/), &
              !                            n_bb)
            end if
            do ixyz = 1, 3
              do jxyz = 1, 3
                susc2el_fde_integrated(jxyz, ixyz) = susc2el_fde_integrated(jxyz, ixyz) &
                                                  + n_bb(jxyz, ixyz)
              end do
            end do
          end if

!         embker-dependent LAO-direct term, only when SDFT contributions
          if (do_london_susc2el_der2) then
            n_bb = 0.0d0
            n_bb_total = 0.0d0
            n_bb_active = 0.0d0
            t = derv_total(1, d0200000) - derv_active(1, d0200000)
            if (use_gga) then
              call gga_london_sus2el_der2(derv_total(1, d0200000),   &
                                          derv_total(1, d0101000),   &
                                          derv_total(1, d0002000),   &
                                          2.0d0*derv_total(1, d0000100),   &
                                          z0tt, &
                                          (/0,0,0/),           &
                                          mat_dim,             &
                                          nz,                  &
                                          dmat_0,              &
                                          ao,                  &
                                          rx(ipoint), ry(ipoint), rz(ipoint), &
                                          buffer, &
                                          n_bb_total)

              call gga_london_sus2el_der2(derv_active(1, d0200000),   &
                                          derv_active(1, d0101000),   &
                                          derv_active(1, d0002000),   &
                                          2.0d0*derv_active(1, d0000100),   &
                                          gn_0_pointwise, &
                                          (/0,0,0/),           &
                                          mat_dim,             &
                                          nz,                  &
                                          dmat_0,              &
                                          ao,                  &
                                          rx(ipoint), ry(ipoint), rz(ipoint), &
                                          buffer, &
                                          n_bb_active)
              do ixyz = 1, 3
                do jxyz = 1, 3
                  susc2el_fde_integrated_sdft(jxyz, ixyz) = susc2el_fde_integrated_sdft(jxyz, ixyz) &
                                                 + n_bb_total(jxyz, ixyz) - n_bb_active(jxyz, ixyz)
                end do
              end do
            else
              call lda_london_sus2el_der2(t,         &
                                          (/0,0,0/), &
                                          mat_dim,   &
                                          nz,        &
                                          dmat_0,    &
                                          ao,        &
                                          rx(ipoint), ry(ipoint), rz(ipoint), &
                                          buffer, &
                                          n_bb)
              do ixyz = 1, 3
                do jxyz = 1, 3
                  susc2el_fde_integrated_sdft(jxyz, ixyz) = susc2el_fde_integrated_sdft(jxyz, ixyz) &
                                                    + n_bb(jxyz, ixyz)
                end do
              end do
            end if
          end if

          if (do_london_susc2el_coupling) then
            if (.not. fde_cfg_no_sdft) then
              call london_susc2el_coupling(ipoint,   &
                                        bllen,       &
                                        mat_dim,     &
                                        nz,          &
                                        ao,          &
                                        n_0_t,       &
                                        z0tt,        &
                                        gnn_0_t,     &
                                        w,           &
                                        use_gga,     &    
                                        derv_length, &
                                        derv_total,  &
                                        rx(ipoint), ry(ipoint), rz(ipoint),&
                                        buffer)
            end if
          end if !do_london_susc2el_coupling

          if (do_london_suscreo_coupling) then
            n_bb = 0.0d0

            call london_suscreo_coupling(ipoint,   &
                                      bllen,       &
                                      mat_dim,     &
                                      nz,          &
                                      dmat,        &
                                      ao,          &
                                      n_0_t,       &
                                      z0tt,        &
                                      w,           &
                                      use_gga,     &    
                                      derv_length, &
                                      derv_total,  &
                                      rx(ipoint), ry(ipoint), rz(ipoint),&
                                      n_bb,        &
                                      buffer)

          end if !do_london_suscreo_coupling

      end if !if (n_0(1) > fde_cfg_small_density_threshold) then
    end do !do ipoint = 1, nr_points_on_this_proc
    end if

    deallocate(ao)
    deallocate(buffer)
    deallocate(derv_total)
    deallocate(derv_active)

    if (allocated(derv_frozen)) deallocate(derv_frozen)
    if (allocated(active_nucpot)) deallocate(active_nucpot)

    if(allocated(u%n)) deallocate(u%n)
    if(allocated(u%gn)) deallocate(u%gn)
    if(allocated(u%s)) deallocate(u%s)
    if(allocated(u%gs)) deallocate(u%gs)
    if(allocated(u%s_b_lao)) deallocate(u%s_b_lao)
    if(allocated(u%gs_b_lao)) deallocate(u%gs_b_lao)
    if(allocated(u%tn)) deallocate(u%tn)
    if(allocated(u%ts)) deallocate(u%ts)

   end subroutine integrate

   
!   ----------------------------------------------------------------------------
      subroutine report()
!   ----------------------------------------------------------------------------
         real(8) :: error_nr_electrons, fde_nr_electrons = 0

         error_nr_electrons = nr_electrons_integrated - fde_nr_electrons

         write(file_unit, *)
         write(file_unit, *) 'FDE integration'
         write(file_unit, *)
         write(file_unit, '(3x, a, i9)')                            &
            'number of grid points                          = ', &
            nr_points_total
         write(file_unit, '(3x, a, i9)')                            &
            'number of grid points below density threshold  = ', &
            nr_points_total - nr_points_used

      write(file_unit, '(3x, a, f26.16)')                        &
            'number of electrons from numerical integration = ', &
            nr_electrons_integrated

      if (do_interaction_energy) then
         write(file_unit, *)
         write(file_unit, '(3x, a)') 'FDE interaction energy contributions'
         write(file_unit, *)
         write(file_unit,'(6x, a, f26.16)') 'non-additive:',non_additive_interaction
         write(file_unit,'(6x, a, f26.16)') 'E_el^B[nA]  :',electrostatic_b_rho_a_interaction
         write(file_unit,'(6x, a, f26.16)') 'E_nuc^A[nB] :',electrostatic_a_rho_b_interaction
      endif 
      write(file_unit, *)

!     timing
      time_integration = second() - time_integration_start
      call timtxt('  time spent in FDE integration                  =', &
                  time_integration, file_unit)
      write(file_unit, *)

   end subroutine report


!   ----------------------------------------------------------------------------
   real(8) function fde_get_nr_electrons_integrated()
!   ----------------------------------------------------------------------------
      fde_get_nr_electrons_integrated = nr_electrons_integrated
   end function fde_get_nr_electrons_integrated


!   ----------------------------------------------------------------------------
      subroutine loop_over_batches()
!   ----------------------------------------------------------------------------
         real(8) :: dummy
         integer :: idummy
         integer :: ibatch
         integer :: i
         real(8), allocatable :: mpi_buffer(:)

!         if (fde_mpi_is_master()) then
!            nr_points_batch = fde_grid_im%npoints
!         end if   
!#ifdef VAR_MPI
!         if (parallel_fde) then
!            call fde_mpi_bcast(nr_points_batch)
!         end if
!#endif
         ibatch          = 0
         nr_points_total = 0
         nr_points_used  = 0

!        loop over batches
         do
         if (fde_mpi_is_master()) then
            nr_points_batch = active_grid%npoints
         end if   
#ifdef VAR_MPI
         if (parallel_fde) then
            call fde_mpi_bcast(nr_points_batch)
         end if
#endif
            if (nr_points_batch < 0) exit

            allocate(rx(nr_points_batch))
            allocate(ry(nr_points_batch))
            allocate(rz(nr_points_batch))
            allocate(rw(nr_points_batch))

            allocate(frozen_n(nr_points_batch))
            allocate(frozen_gnx(nr_points_batch))
            allocate(frozen_gny(nr_points_batch))
            allocate(frozen_gnz(nr_points_batch))
            allocate(frozen_elpot(nr_points_batch))

            if (use_potential) then
               allocate(use_vemb(nr_points_batch))
            end if

!radovan: grid stuff will be moved outside this module
            if (fde_mpi_is_master()) then
                rx = active_grid%r(1,:)
                ry = active_grid%r(2,:)
                rz = active_grid%r(3,:)
                rw = active_grid%w
 
                if (use_potential) then
                   use_vemb(:) = fde_static_vemb(:)
                end if

                frozen_n     = active_gridfunction%n
                frozen_gnx   = active_gridfunction%gn(1,:)
                frozen_gny   = active_gridfunction%gn(2,:)
                frozen_gnz   = active_gridfunction%gn(3,:)
                frozen_elpot = active_gridfunction%elpot

            end if

            if (do_london_rhs_reorth_coupling &
              .or. do_london_rhs_direct_coupling &
              .or. do_london_susc2el_coupling &
              .or. do_london_suscreo_coupling) then
               allocate(frozen_n_bx(nr_points_batch))
               allocate(frozen_n_by(nr_points_batch))
               allocate(frozen_n_bz(nr_points_batch))
               frozen_n_bx = 0.0d0
               frozen_n_by = 0.0d0
               frozen_n_bz = 0.0d0
               if (use_gga) then
                  allocate(frozen_gxn_bx(nr_points_batch))
                  allocate(frozen_gxn_by(nr_points_batch))
                  allocate(frozen_gxn_bz(nr_points_batch))
                  allocate(frozen_gyn_bx(nr_points_batch))
                  allocate(frozen_gyn_by(nr_points_batch))
                  allocate(frozen_gyn_bz(nr_points_batch))
                  allocate(frozen_gzn_bx(nr_points_batch))
                  allocate(frozen_gzn_by(nr_points_batch))
                  allocate(frozen_gzn_bz(nr_points_batch))
                  frozen_gxn_bx = 0.0d0  
                  frozen_gxn_by = 0.0d0
                  frozen_gxn_bz = 0.0d0
                  frozen_gyn_bx = 0.0d0
                  frozen_gyn_by = 0.0d0
                  frozen_gyn_bz = 0.0d0
                  frozen_gzn_bx = 0.0d0
                  frozen_gzn_by = 0.0d0
                  frozen_gzn_bz = 0.0d0
               end if
               if (.not. fde_cfg_no_sdft) then
                  allocate(frozen_sx_bx(nr_points_batch))
                  allocate(frozen_sx_by(nr_points_batch))
                  allocate(frozen_sx_bz(nr_points_batch))
                  allocate(frozen_sy_bx(nr_points_batch))
                  allocate(frozen_sy_by(nr_points_batch))
                  allocate(frozen_sy_bz(nr_points_batch))
                  allocate(frozen_sz_bx(nr_points_batch))
                  allocate(frozen_sz_by(nr_points_batch))
                  allocate(frozen_sz_bz(nr_points_batch))
                  frozen_sx_bx = 0.0d0  
                  frozen_sx_by = 0.0d0
                  frozen_sx_bz = 0.0d0
                  frozen_sy_bx = 0.0d0
                  frozen_sy_by = 0.0d0
                  frozen_sy_bz = 0.0d0
                  frozen_sz_bx = 0.0d0
                  frozen_sz_by = 0.0d0
                  frozen_sz_bz = 0.0d0
                  if (use_gga) then
                     allocate(frozen_gxsx_bx(nr_points_batch))
                     allocate(frozen_gxsx_by(nr_points_batch))
                     allocate(frozen_gxsx_bz(nr_points_batch))
                     allocate(frozen_gxsy_bx(nr_points_batch))
                     allocate(frozen_gxsy_by(nr_points_batch))
                     allocate(frozen_gxsy_bz(nr_points_batch))
                     allocate(frozen_gxsz_bx(nr_points_batch))
                     allocate(frozen_gxsz_by(nr_points_batch))
                     allocate(frozen_gxsz_bz(nr_points_batch))
                     allocate(frozen_gysx_bx(nr_points_batch))
                     allocate(frozen_gysx_by(nr_points_batch))
                     allocate(frozen_gysx_bz(nr_points_batch))
                     allocate(frozen_gysy_bx(nr_points_batch))
                     allocate(frozen_gysy_by(nr_points_batch))
                     allocate(frozen_gysy_bz(nr_points_batch))
                     allocate(frozen_gysz_bx(nr_points_batch))
                     allocate(frozen_gysz_by(nr_points_batch))
                     allocate(frozen_gysz_bz(nr_points_batch))
                     allocate(frozen_gzsx_bx(nr_points_batch))
                     allocate(frozen_gzsx_by(nr_points_batch))
                     allocate(frozen_gzsx_bz(nr_points_batch))
                     allocate(frozen_gzsy_bx(nr_points_batch))
                     allocate(frozen_gzsy_by(nr_points_batch))
                     allocate(frozen_gzsy_bz(nr_points_batch))
                     allocate(frozen_gzsz_bx(nr_points_batch))
                     allocate(frozen_gzsz_by(nr_points_batch))
                     allocate(frozen_gzsz_bz(nr_points_batch))
                     frozen_gxsx_bx = 0.0d0  
                     frozen_gxsx_by = 0.0d0
                     frozen_gxsx_bz = 0.0d0
                     frozen_gxsy_bx = 0.0d0
                     frozen_gxsy_by = 0.0d0
                     frozen_gxsy_bz = 0.0d0
                     frozen_gxsz_bx = 0.0d0
                     frozen_gxsz_by = 0.0d0
                     frozen_gxsz_bz = 0.0d0
                     frozen_gysx_bx = 0.0d0
                     frozen_gysx_by = 0.0d0
                     frozen_gysx_bz = 0.0d0
                     frozen_gysy_bx = 0.0d0
                     frozen_gysy_by = 0.0d0
                     frozen_gysy_bz = 0.0d0
                     frozen_gysz_bx = 0.0d0
                     frozen_gysz_by = 0.0d0
                     frozen_gysz_bz = 0.0d0
                     frozen_gzsx_bx = 0.0d0
                     frozen_gzsx_by = 0.0d0
                     frozen_gzsx_bz = 0.0d0
                     frozen_gzsy_bx = 0.0d0
                     frozen_gzsy_by = 0.0d0
                     frozen_gzsy_bz = 0.0d0
                     frozen_gzsz_bx = 0.0d0
                     frozen_gzsz_by = 0.0d0
                     frozen_gzsz_bz = 0.0d0
                  end if
               end if
               if (fde_mpi_is_master()) then
                  if (do_london_rhs_direct_coupling.or.do_london_susc2el_coupling) then
                    if (.not. fde_cfg_no_sdft) then
                       frozen_sx_bx(:)  = active_gridfunction%s_b_direct(1,1,:)
                       frozen_sy_bx(:)  = active_gridfunction%s_b_direct(2,1,:)
                       frozen_sz_bx(:)  = active_gridfunction%s_b_direct(3,1,:)
                       frozen_sx_by(:)  = active_gridfunction%s_b_direct(1,2,:)
                       frozen_sy_by(:)  = active_gridfunction%s_b_direct(2,2,:)
                       frozen_sz_by(:)  = active_gridfunction%s_b_direct(3,2,:)
                       frozen_sx_bz(:)  = active_gridfunction%s_b_direct(1,3,:)
                       frozen_sy_bz(:)  = active_gridfunction%s_b_direct(2,3,:)
                       frozen_sz_bz(:)  = active_gridfunction%s_b_direct(3,3,:)
!debug print        
                       !do i = 1, nr_points_batch
                       !   write(*, '(A,9ES30.20E3)') 'test15 import frozen s = ', &
                       !               frozen_sx_bx(i), frozen_sy_bx(i), frozen_sz_bx(i), &
                       !               frozen_sx_by(i), frozen_sy_by(i), frozen_sz_by(i), &
                       !               frozen_sx_bz(i), frozen_sy_bz(i), frozen_sz_bz(i)
                       !end do
                       if (use_gga) then
                          frozen_gxsx_bx(:)  = active_gridfunction%gs_b_direct(1,1,1,:)
                          frozen_gysx_bx(:)  = active_gridfunction%gs_b_direct(2,1,1,:)
                          frozen_gzsx_bx(:)  = active_gridfunction%gs_b_direct(3,1,1,:)
                          frozen_gxsy_bx(:)  = active_gridfunction%gs_b_direct(1,2,1,:)
                          frozen_gysy_bx(:)  = active_gridfunction%gs_b_direct(2,2,1,:)
                          frozen_gzsy_bx(:)  = active_gridfunction%gs_b_direct(3,2,1,:)
                          frozen_gxsz_bx(:)  = active_gridfunction%gs_b_direct(1,3,1,:)
                          frozen_gysz_bx(:)  = active_gridfunction%gs_b_direct(2,3,1,:)
                          frozen_gzsz_bx(:)  = active_gridfunction%gs_b_direct(3,3,1,:)
                          frozen_gxsx_by(:)  = active_gridfunction%gs_b_direct(1,1,2,:)
                          frozen_gysx_by(:)  = active_gridfunction%gs_b_direct(2,1,2,:)
                          frozen_gzsx_by(:)  = active_gridfunction%gs_b_direct(3,1,2,:)
                          frozen_gxsy_by(:)  = active_gridfunction%gs_b_direct(1,2,2,:)
                          frozen_gysy_by(:)  = active_gridfunction%gs_b_direct(2,2,2,:)
                          frozen_gzsy_by(:)  = active_gridfunction%gs_b_direct(3,2,2,:)
                          frozen_gxsz_by(:)  = active_gridfunction%gs_b_direct(1,3,2,:)
                          frozen_gysz_by(:)  = active_gridfunction%gs_b_direct(2,3,2,:)
                          frozen_gzsz_by(:)  = active_gridfunction%gs_b_direct(3,3,2,:)
                          frozen_gxsx_bz(:)  = active_gridfunction%gs_b_direct(1,1,3,:)
                          frozen_gysx_bz(:)  = active_gridfunction%gs_b_direct(2,1,3,:)
                          frozen_gzsx_bz(:)  = active_gridfunction%gs_b_direct(3,1,3,:)
                          frozen_gxsy_bz(:)  = active_gridfunction%gs_b_direct(1,2,3,:)
                          frozen_gysy_bz(:)  = active_gridfunction%gs_b_direct(2,2,3,:)
                          frozen_gzsy_bz(:)  = active_gridfunction%gs_b_direct(3,2,3,:)
                          frozen_gxsz_bz(:)  = active_gridfunction%gs_b_direct(1,3,3,:)
                          frozen_gysz_bz(:)  = active_gridfunction%gs_b_direct(2,3,3,:)
                          frozen_gzsz_bz(:)  = active_gridfunction%gs_b_direct(3,3,3,:)
!debug print        
                       !do i = 1, nr_points_batch
                       !   write(*, '(A,27ES30.20E3)') 'test16 import frozen gs = ', &
                       !               frozen_gxsx_bx(i), frozen_gysx_bx(i), frozen_gzsx_bx(i), &
                       !               frozen_gxsy_bx(i), frozen_gysy_bx(i), frozen_gzsy_bx(i), &
                       !               frozen_gxsz_bx(i), frozen_gysz_bx(i), frozen_gzsz_bx(i), &
                       !               frozen_gxsx_by(i), frozen_gysx_by(i), frozen_gzsx_by(i), &
                       !               frozen_gxsy_by(i), frozen_gysy_by(i), frozen_gzsy_by(i), &
                       !               frozen_gxsz_by(i), frozen_gysz_by(i), frozen_gzsz_by(i), &
                       !               frozen_gxsx_bz(i), frozen_gysx_bz(i), frozen_gzsx_bz(i), &
                       !               frozen_gxsy_bz(i), frozen_gysy_bz(i), frozen_gzsy_bz(i), &
                       !               frozen_gxsz_bz(i), frozen_gysz_bz(i), frozen_gzsz_bz(i)
                       !end do
                       end if
                    end if
                  end if
!
                  if (do_london_rhs_reorth_coupling) then
                    frozen_n_bx  = active_gridfunction%n_b_reorth(1,:)
                    frozen_n_by  = active_gridfunction%n_b_reorth(2,:)
                    frozen_n_bz  = active_gridfunction%n_b_reorth(3,:)
!debug print        
                    !do i = 1, nr_points_batch
                    !  write(*, '(A,3ES30.20E3)') 'test12 import frozen n = ', &
                    !              frozen_n_bx(i), frozen_n_by(i), frozen_n_bz(i)
                    !end do
                    if (use_gga) then
                       frozen_gxn_bx(:)  = active_gridfunction%gn_b_reorth(1,1,:)
                       frozen_gyn_bx(:)  = active_gridfunction%gn_b_reorth(2,1,:)
                       frozen_gzn_bx(:)  = active_gridfunction%gn_b_reorth(3,1,:)
                       frozen_gxn_by(:)  = active_gridfunction%gn_b_reorth(1,2,:)
                       frozen_gyn_by(:)  = active_gridfunction%gn_b_reorth(2,2,:)
                       frozen_gzn_by(:)  = active_gridfunction%gn_b_reorth(3,2,:)
                       frozen_gxn_bz(:)  = active_gridfunction%gn_b_reorth(1,3,:)
                       frozen_gyn_bz(:)  = active_gridfunction%gn_b_reorth(2,3,:)
                       frozen_gzn_bz(:)  = active_gridfunction%gn_b_reorth(3,3,:)
!debug print        
                       !do i = 1, nr_points_batch
                       !   write(*, '(A,9ES30.20E3)') 'test13 import frozen gn = ', &
                       !               frozen_gxn_bx(i), frozen_gyn_bx(i), frozen_gzn_bx(i), &
                       !               frozen_gxn_by(i), frozen_gyn_by(i), frozen_gzn_by(i), &
                       !               frozen_gxn_bz(i), frozen_gyn_bz(i), frozen_gzn_bz(i)
                       !end do
                    end if
                    if (.not. fde_cfg_no_sdft) then
                       frozen_sx_bx(:)  = active_gridfunction%s_b_reorth(1,1,:)
                       frozen_sy_bx(:)  = active_gridfunction%s_b_reorth(2,1,:)
                       frozen_sz_bx(:)  = active_gridfunction%s_b_reorth(3,1,:)
                       frozen_sx_by(:)  = active_gridfunction%s_b_reorth(1,2,:)
                       frozen_sy_by(:)  = active_gridfunction%s_b_reorth(2,2,:)
                       frozen_sz_by(:)  = active_gridfunction%s_b_reorth(3,2,:)
                       frozen_sx_bz(:)  = active_gridfunction%s_b_reorth(1,3,:)
                       frozen_sy_bz(:)  = active_gridfunction%s_b_reorth(2,3,:)
                       frozen_sz_bz(:)  = active_gridfunction%s_b_reorth(3,3,:)
!debug print        
                       !do i = 1, nr_points_batch
                       !   write(*, '(A,9ES30.20E3)') 'test15 import frozen s = ', &
                       !               frozen_sx_bx(i), frozen_sy_bx(i), frozen_sz_bx(i), &
                       !               frozen_sx_by(i), frozen_sy_by(i), frozen_sz_by(i), &
                       !               frozen_sx_bz(i), frozen_sy_bz(i), frozen_sz_bz(i)
                       !end do
                       if (use_gga) then
                          frozen_gxsx_bx(:)  = active_gridfunction%gs_b_reorth(1,1,1,:)
                          frozen_gysx_bx(:)  = active_gridfunction%gs_b_reorth(2,1,1,:)
                          frozen_gzsx_bx(:)  = active_gridfunction%gs_b_reorth(3,1,1,:)
                          frozen_gxsy_bx(:)  = active_gridfunction%gs_b_reorth(1,2,1,:)
                          frozen_gysy_bx(:)  = active_gridfunction%gs_b_reorth(2,2,1,:)
                          frozen_gzsy_bx(:)  = active_gridfunction%gs_b_reorth(3,2,1,:)
                          frozen_gxsz_bx(:)  = active_gridfunction%gs_b_reorth(1,3,1,:)
                          frozen_gysz_bx(:)  = active_gridfunction%gs_b_reorth(2,3,1,:)
                          frozen_gzsz_bx(:)  = active_gridfunction%gs_b_reorth(3,3,1,:)
                          frozen_gxsx_by(:)  = active_gridfunction%gs_b_reorth(1,1,2,:)
                          frozen_gysx_by(:)  = active_gridfunction%gs_b_reorth(2,1,2,:)
                          frozen_gzsx_by(:)  = active_gridfunction%gs_b_reorth(3,1,2,:)
                          frozen_gxsy_by(:)  = active_gridfunction%gs_b_reorth(1,2,2,:)
                          frozen_gysy_by(:)  = active_gridfunction%gs_b_reorth(2,2,2,:)
                          frozen_gzsy_by(:)  = active_gridfunction%gs_b_reorth(3,2,2,:)
                          frozen_gxsz_by(:)  = active_gridfunction%gs_b_reorth(1,3,2,:)
                          frozen_gysz_by(:)  = active_gridfunction%gs_b_reorth(2,3,2,:)
                          frozen_gzsz_by(:)  = active_gridfunction%gs_b_reorth(3,3,2,:)
                          frozen_gxsx_bz(:)  = active_gridfunction%gs_b_reorth(1,1,3,:)
                          frozen_gysx_bz(:)  = active_gridfunction%gs_b_reorth(2,1,3,:)
                          frozen_gzsx_bz(:)  = active_gridfunction%gs_b_reorth(3,1,3,:)
                          frozen_gxsy_bz(:)  = active_gridfunction%gs_b_reorth(1,2,3,:)
                          frozen_gysy_bz(:)  = active_gridfunction%gs_b_reorth(2,2,3,:)
                          frozen_gzsy_bz(:)  = active_gridfunction%gs_b_reorth(3,2,3,:)
                          frozen_gxsz_bz(:)  = active_gridfunction%gs_b_reorth(1,3,3,:)
                          frozen_gysz_bz(:)  = active_gridfunction%gs_b_reorth(2,3,3,:)
                          frozen_gzsz_bz(:)  = active_gridfunction%gs_b_reorth(3,3,3,:)
!debug print        
                       !do i = 1, nr_points_batch
                       !   write(*, '(A,27ES30.20E3)') 'test16 import frozen gs = ', &
                       !               frozen_gxsx_bx(i), frozen_gysx_bx(i), frozen_gzsx_bx(i), &
                       !               frozen_gxsy_bx(i), frozen_gysy_bx(i), frozen_gzsy_bx(i), &
                       !               frozen_gxsz_bx(i), frozen_gysz_bx(i), frozen_gzsz_bx(i), &
                       !               frozen_gxsx_by(i), frozen_gysx_by(i), frozen_gzsx_by(i), &
                       !               frozen_gxsy_by(i), frozen_gysy_by(i), frozen_gzsy_by(i), &
                       !               frozen_gxsz_by(i), frozen_gysz_by(i), frozen_gzsz_by(i), &
                       !               frozen_gxsx_bz(i), frozen_gysx_bz(i), frozen_gzsx_bz(i), &
                       !               frozen_gxsy_bz(i), frozen_gysy_bz(i), frozen_gzsy_bz(i), &
                       !               frozen_gxsz_bz(i), frozen_gysz_bz(i), frozen_gzsz_bz(i)
                       !end do
                       end if
                    end if
                  end if
               end if
            end if

            ibatch = ibatch + 1
            nr_points_total = nr_points_total + nr_points_batch
            nr_points_on_this_proc = nr_points_batch

#ifdef VAR_MPI
         if (parallel_fde) then
            call fde_mpi_distribute_points(rx, ry, rz, rw, nr_points_batch, nr_points_on_this_proc) 
            call fde_mpi_distribute_frozen(frozen_n, frozen_gnx, frozen_gny, frozen_gnz, &
       &                                   frozen_elpot,nr_points_batch, nr_points_on_this_proc) 
            if (use_potential) then
              call fde_mpi_distribute_vemb(use_vemb, nr_points_batch, nr_points_on_this_proc)
            end if
            if (do_london_rhs_reorth_coupling &
              .or. do_london_rhs_reorth_coupling &
              .or. do_london_rhs_direct_coupling &
              .or. do_london_susc2el_coupling &
              .or. do_london_suscreo_coupling) then
              call fde_mpi_distribute_frozen_pertden_n(frozen_n_bx, frozen_n_by, frozen_n_bz, &
        &                                              nr_points_batch, nr_points_on_this_proc)
              if (use_gga) then
              call fde_mpi_distribute_frozen_pertden_gn(frozen_gxn_bx, frozen_gxn_by, frozen_gxn_bz, &
                                                        frozen_gyn_bx, frozen_gyn_by, frozen_gyn_bz, &
                                                        frozen_gzn_bx, frozen_gzn_by, frozen_gzn_bz, &
                                                        nr_points_batch, nr_points_on_this_proc)
              end if
              if (.not. fde_cfg_no_sdft) then
                 call fde_mpi_distribute_frozen_pertden_s(frozen_sx_bx, frozen_sx_by, frozen_sx_bz, &
                                                          frozen_sy_bx, frozen_sy_by, frozen_sy_bz, &
                                                          frozen_sz_bx, frozen_sz_by, frozen_sz_bz, &
                                                          nr_points_batch, nr_points_on_this_proc)
                 if (use_gga) then
                    call fde_mpi_distribute_frozen_pertden_gs(frozen_gxsx_bx, frozen_gxsx_by, frozen_gxsx_bz, &
                                                              frozen_gxsy_bx, frozen_gxsy_by, frozen_gxsy_bz, &
                                                              frozen_gxsz_bx, frozen_gxsz_by, frozen_gxsz_bz, &
                                                              frozen_gysx_bx, frozen_gysx_by, frozen_gysx_bz, &
                                                              frozen_gysy_bx, frozen_gysy_by, frozen_gysy_bz, &
                                                              frozen_gysz_bx, frozen_gysz_by, frozen_gysz_bz, &
                                                              frozen_gzsx_bx, frozen_gzsx_by, frozen_gzsx_bz, &
                                                              frozen_gzsy_bx, frozen_gzsy_by, frozen_gzsy_bz, &
                                                              frozen_gzsz_bx, frozen_gzsz_by, frozen_gzsz_bz, &
                                                              nr_points_batch, nr_points_on_this_proc)
                 end if
              end if
            end if
         end if
#endif
!        if (do_potential .or. (response_order_mo == 1)) then
!            fde_use_blas3_scheme = .true.
!         else
!            fde_use_blas3_scheme = .false.
!         end if
!         if (fun_is_tau_mgga(xc_fun)) then
!            fde_use_blas3_scheme = .false.
!         end if
!         if (fde_cfg_blocked) then
!            fde_use_blas3_scheme = .true.
!         end if
!         if (fde_cfg_pointwise) then
!            fde_use_blas3_scheme = .false.
!         end if

            call integrate()

            deallocate(rx)
            deallocate(ry)
            deallocate(rz)
            deallocate(rw)
            deallocate(frozen_n)
            deallocate(frozen_gnx)
            deallocate(frozen_gny)
            deallocate(frozen_gnz)
            deallocate(frozen_elpot)
            if (allocated(use_vemb)) deallocate(use_vemb)
            if (allocated(frozen_n_bx)) deallocate(frozen_n_bx)
            if (allocated(frozen_n_by)) deallocate(frozen_n_by)
            if (allocated(frozen_n_bz)) deallocate(frozen_n_bz)
            if (allocated(frozen_gxn_bx)) deallocate(frozen_gxn_bx)
            if (allocated(frozen_gxn_by)) deallocate(frozen_gxn_by)
            if (allocated(frozen_gxn_bz)) deallocate(frozen_gxn_bz)
            if (allocated(frozen_gyn_bx)) deallocate(frozen_gyn_bx)
            if (allocated(frozen_gyn_by)) deallocate(frozen_gyn_by)
            if (allocated(frozen_gyn_bz)) deallocate(frozen_gyn_bz)
            if (allocated(frozen_gzn_bx)) deallocate(frozen_gzn_bx)
            if (allocated(frozen_gzn_by)) deallocate(frozen_gzn_by)
            if (allocated(frozen_gzn_bz)) deallocate(frozen_gzn_bz)
            if (allocated(frozen_sx_bx)) deallocate(frozen_sx_bx)
            if (allocated(frozen_sx_by)) deallocate(frozen_sx_by)
            if (allocated(frozen_sx_bz)) deallocate(frozen_sx_bz)
            if (allocated(frozen_sy_bx)) deallocate(frozen_sy_bx)
            if (allocated(frozen_sy_by)) deallocate(frozen_sy_by)
            if (allocated(frozen_sy_bz)) deallocate(frozen_sy_bz)
            if (allocated(frozen_sz_bx)) deallocate(frozen_sz_bx)
            if (allocated(frozen_sz_by)) deallocate(frozen_sz_by)
            if (allocated(frozen_sz_bz)) deallocate(frozen_sz_bz)
            if (allocated(frozen_gxsx_bx)) deallocate(frozen_gxsx_bx)
            if (allocated(frozen_gxsx_by)) deallocate(frozen_gxsx_by)
            if (allocated(frozen_gxsx_bz)) deallocate(frozen_gxsx_bz)
            if (allocated(frozen_gxsy_bx)) deallocate(frozen_gxsy_bx)
            if (allocated(frozen_gxsy_by)) deallocate(frozen_gxsy_by)
            if (allocated(frozen_gxsy_bz)) deallocate(frozen_gxsy_bz)
            if (allocated(frozen_gxsz_bx)) deallocate(frozen_gxsz_bx)
            if (allocated(frozen_gxsz_by)) deallocate(frozen_gxsz_by)
            if (allocated(frozen_gxsz_bz)) deallocate(frozen_gxsz_bz)
            if (allocated(frozen_gysx_bx)) deallocate(frozen_gysx_bx)
            if (allocated(frozen_gysx_by)) deallocate(frozen_gysx_by)
            if (allocated(frozen_gysx_bz)) deallocate(frozen_gysx_bz)
            if (allocated(frozen_gysy_bx)) deallocate(frozen_gysy_bx)
            if (allocated(frozen_gysy_by)) deallocate(frozen_gysy_by)
            if (allocated(frozen_gysy_bz)) deallocate(frozen_gysy_bz)
            if (allocated(frozen_gysz_bx)) deallocate(frozen_gysz_bx)
            if (allocated(frozen_gysz_by)) deallocate(frozen_gysz_by)
            if (allocated(frozen_gysz_bz)) deallocate(frozen_gysz_bz)
            if (allocated(frozen_gzsx_bx)) deallocate(frozen_gzsx_bx)
            if (allocated(frozen_gzsx_by)) deallocate(frozen_gzsx_by)
            if (allocated(frozen_gzsx_bz)) deallocate(frozen_gzsx_bz)
            if (allocated(frozen_gzsy_bx)) deallocate(frozen_gzsy_bx)
            if (allocated(frozen_gzsy_by)) deallocate(frozen_gzsy_by)
            if (allocated(frozen_gzsy_bz)) deallocate(frozen_gzsy_bz)
            if (allocated(frozen_gzsz_bx)) deallocate(frozen_gzsz_bx)
            if (allocated(frozen_gzsz_by)) deallocate(frozen_gzsz_by)
            if (allocated(frozen_gzsz_bz)) deallocate(frozen_gzsz_bz)
            if (allocated(use_vemb))       deallocate(use_vemb)



!          don't try to read another batch
           goto 999

!     end loop over batches
         end do
 999     continue
         
      end subroutine loop_over_batches


!   ----------------------------------------------------------------------------
      subroutine lr_fde_uncoupled(block_length,       &
                    mat_dim,            &
                    nz,                 &
                    nr_mat,             &
                    dmaterturbed,       &
                    isym,               &
                    ih,                 &
                    allow_nonhermitian, &
                    ao,                 &
                    n_0,                &
                    n_0_t,              &
                    gn_0_pointwise,     &
                    gn_0_t_pointwise,   &
                    gnn_0,              &
                    gnn_0_t,            &
                    w,                  &
                    is_gga_qr,          &
                    is_gga_qi,          &
                    buffer,             &
                    n,                  &
                    derv_active,        &
                    derv_total,         &
                    u)
!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    integer, intent(in)    :: nr_mat
    real(8), intent(in)    :: dmaterturbed(*)
    integer, intent(in)    :: isym(nr_mat)
    integer, intent(in)    :: ih(nr_mat)
    logical, intent(in)    :: allow_nonhermitian
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: n_0(*)
    real(8), intent(in)    :: n_0_t(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: gn_0_t_pointwise(*)
    real(8), intent(in)    :: gnn_0(*)
    real(8), intent(in)    :: gnn_0_t(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga_qr
    logical, intent(in)    :: is_gga_qi
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv_active(max_block_length, 0:n)
    real(8), intent(in)    :: derv_total(max_block_length, 0:n)
!   ----------------------------------------------------------------------------
    logical                :: calculate_qr, calculate_qi
    logical, save          :: have_sdft_contrib(100) = .false.
    logical, save          :: have_sdft = .false.
    logical, save          :: print_sdft_contrib = .true.


    real(8)                ::  r_b,     s_b(3)
    real(8)                :: gr_b(3), gs_b(3, 3)
    real(8)                ::  t_b, ts_b(3)

    real(8)                :: z_b
    real(8)                :: z_b_t
    real(8)                :: y_b(3)
    real(8)                :: y_b_t(3)

    integer                :: k_start, k, im, im_off, j
    type(fde_omega_prefactor) :: u
!   ----------------------------------------------------------------------------

    k_start = 1

    do im = 1, nr_mat

      calculate_qr = .false.
      calculate_qi = .false.

      if (ih(im) == +1) calculate_qr = .true.
      if (ih(im) == -1) calculate_qi = .true.

      if (ih(im) == 0 .and. allow_nonhermitian) then
        calculate_qr = .true.
        calculate_qi = .true.
      end if

      if (fde_cfg_no_sdft) then
        calculate_qi = .false.
      end if

      if (calculate_qi) have_sdft_contrib(im) = .true.

      im_off = mat_dim*mat_dim*nz*(im - 1)

!     hermitian contribution (qr = quaternion real)
      if (calculate_qr) then

         if (is_gga_qr) then
            if (allow_nonhermitian) then
               call get_gn_nonhermitian(r_b,                      &
                                        gr_b,                     &
                                        isym(im) - 1,             &
                                        mat_dim,                  &
                                        dmaterturbed(1 + im_off), &
                                        buffer,                   &
                                        ao)
            else
               call get_gn(r_b,                      &
                           gr_b,                     &
                           isym(im) - 1,             &
                           mat_dim,                  &
                           dmaterturbed(1 + im_off), &
                           buffer,                   &
                           ao)
            end if
         else
            call get_n(r_b,                        &
                       isym(im) - 1,               &
                       mat_dim,                    &
                       dmaterturbed(1 + im_off), &
                       buffer,                     &
                       ao)
         end if

        u%n(im) = u%n(im) + (derv_total(1, d2000000) - derv_active(1, d2000000))*r_b

        if (is_gga_qr) then
           z_b = 2.0d0*(gn_0_pointwise(1)*gr_b(1) &
                      + gn_0_pointwise(2)*gr_b(2) &
                      + gn_0_pointwise(3)*gr_b(3))

           z_b_t = 2.0d0*(gn_0_t_pointwise(1)*gr_b(1) &
                        + gn_0_t_pointwise(2)*gr_b(2) &
                        + gn_0_t_pointwise(3)*gr_b(3))

           u%n(im) = u%n(im) + (derv_total(1, d1010000)*z_b_t - derv_active(1, d1010000)*z_b)

           u%gn(1:3, im) = u%gn(1:3, im) + 2.0d0*gn_0_t_pointwise(1:3)*r_b*derv_total(1, d1010000)   &
                                         + 2.0d0*gn_0_t_pointwise(1:3)*z_b_t*derv_total(1, d0020000) &
                                         + 2.0d0*gr_b(1:3)*derv_total(1, d0010000)                   &
                                         - 2.0d0*gn_0_pointwise(1:3)*r_b*derv_active(1, d1010000)    &
                                         - 2.0d0*gn_0_pointwise(1:3)*z_b*derv_active(1, d0020000)    &
                                         - 2.0d0*gr_b(1:3)*derv_active(1, d0010000)
        end if

      end if !end of qr contribution

!     antihermitian contribution (qi = quaternion imaginary)
      if (calculate_qi) then


        if (is_gga_qi) then
          if (allow_nonhermitian) then
            call get_gs_nonhermitian(s_b,                       &
                                     gs_b,                       &
                                     isym(im) - 1,               &
                                     mat_dim,                    &
                                     dmaterturbed(1 + im_off), &
                                     buffer,                     &
                                     ao)
          else
            call get_gs(s_b,                       &
                        gs_b,                       &
                        isym(im) - 1,               &
                        mat_dim,                    &
                        dmaterturbed(1 + im_off), &
                        buffer,                     &
                        ao)
          end if
        else
        call get_s(s_b,                        &
                   isym(im) - 1,               &
                   mat_dim,                    &
                   dmaterturbed(1 + im_off), &
                   buffer,                     &
                   ao)
        end if

        do k = k_start, 3
           u%s(k, im) = u%s(k, im) + (derv_total(1, d0200000) - derv_active(1, d0200000))*s_b(k)

           if (is_gga_qi) then
              y_b(k)  = gn_0_pointwise(1)*gs_b(1, k) &
                      + gn_0_pointwise(2)*gs_b(2, k) &
                      + gn_0_pointwise(3)*gs_b(3, k)

              y_b_t(k) = gn_0_t_pointwise(1)*gs_b(1, k) &
                       + gn_0_t_pointwise(2)*gs_b(2, k) &
                       + gn_0_t_pointwise(3)*gs_b(3, k)

              u%s(k, im) = u%s(k, im) + (derv_total(1, d0101000)*y_b_t(k) - derv_active(1, d0101000)*y_b(k))

              u%gs(1:3, k, im) = u%gs(1:3, k, im) +       gn_0_t_pointwise(1:3)*s_b(k)*derv_total(1, d0101000)   &
                                                  +       gn_0_t_pointwise(1:3)*y_b_t(k)*derv_total(1, d0002000) &
                                                  + 2.0d0*gs_b(1:3, k)*derv_total(1, d0000100)                   &
                                                  -       gn_0_pointwise(1:3)*s_b(k)*derv_active(1, d0101000)    &
                                                  -       gn_0_pointwise(1:3)*y_b(k)*derv_active(1, d0002000)    &
                                                  - 2.0d0*gs_b(1:3, k)*derv_active(1, d0000100)
           end if
        end do

      end if !end of qi contribution

    end do

! now for some output
    if (print_sdft_contrib) then
       do im = 1, nr_mat
          if (have_sdft_contrib(im)) have_sdft = .true.
       end do
 
       if (fde_mpi_is_master()) then
          if (have_sdft) then
             write(*,*) 'SDFT contributions to FDE kernel included'
          else
             write(*,*) 'SDFT contributions to FDE kernel neglected'
          endif
       endif
       print_sdft_contrib = .false.
    endif
!

  end subroutine

  subroutine lao_reorth(block_length,       &
                        mat_dim,            &
                        nz,                 &
                        ao,                 &
                        n_0,                &
                        n_0_t,              &
                        gn_0_pointwise,     &
                        gn_0_t_pointwise,   &
                        gnn_0,              &
                        gnn_0_t,            &
                        w,                  &
                        is_gga,             &
                        tb_dmat,            &
                        isym,               &
                        buffer,             &
                        n,                  &
                        derv_active,        &
                        derv_total,         &
                        u)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: n_0(*)
    real(8), intent(in)    :: n_0_t(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: gn_0_t_pointwise(*)
    real(8), intent(in)    :: gnn_0(*)
    real(8), intent(in)    :: gnn_0_t(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga
    real(8), intent(in)    :: tb_dmat(*)
    integer, intent(in)    :: isym(*)
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv_active(max_block_length, 0:n)
    real(8), intent(in)    :: derv_total(max_block_length, 0:n)
!   ----------------------------------------------------------------------------
    real(8)                ::  r_b,     s_b(3)
    real(8)                :: gr_b(3), gs_b(3, 3)
    real(8)                ::  t_b,    ts_b(3)

    real(8)                :: z_b
    real(8)                :: z_b_t
    real(8)                :: y_b(3)
    real(8)                :: y_b_t(3)

    integer                :: k_start, k, im, im_off, j, irep
    type(fde_omega_prefactor) :: u
    integer :: i
!   ----------------------------------------------------------------------------

    k_start = 1
    if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

    do im = 1, nr_dmat

       im_off = mat_dim*mat_dim*nz*(im - 1)
       irep = isym(im) - 1

!      charge density part
       if (is_gga) then
          call get_gn_nonhermitian(r_b,                 &
                                   gr_b,                &
                                   irep,                &
                                   mat_dim,             &
                                   tb_dmat(1 + im_off), &
                                   buffer,              &
                                   ao)
       else
          call get_n(r_b,                 &
                     irep,                &
                     mat_dim,             &
                     tb_dmat(1 + im_off), &
                     buffer,              &
                     ao)
       end if

       u%n(im) = u%n(im) + (derv_total(1, d2000000) - derv_active(1, d2000000))*r_b

       if (is_gga) then
          z_b = 2.0d0*(gn_0_pointwise(1)*gr_b(1) &
                     + gn_0_pointwise(2)*gr_b(2) &
                     + gn_0_pointwise(3)*gr_b(3))

          z_b_t = 2.0d0*(gn_0_t_pointwise(1)*gr_b(1) &
                       + gn_0_t_pointwise(2)*gr_b(2) &
                       + gn_0_t_pointwise(3)*gr_b(3))


          u%n(im) = u%n(im) + (derv_total(1, d1010000)*z_b_t - derv_active(1, d1010000)*z_b)

          u%gn(1:3, im) = u%gn(1:3, im) + 2.0d0*gn_0_t_pointwise(1:3)*r_b*derv_total(1, d1010000)   &
                                        + 2.0d0*gn_0_t_pointwise(1:3)*z_b_t*derv_total(1, d0020000) &
                                        + 2.0d0*gr_b(1:3)*derv_total(1, d0010000)                   &
                                        - 2.0d0*gn_0_pointwise(1:3)*r_b*derv_active(1, d1010000)    &
                                        - 2.0d0*gn_0_pointwise(1:3)*z_b*derv_active(1, d0020000)    &
                                        - 2.0d0*gr_b(1:3)*derv_active(1, d0010000)

       end if

!      SDFT
       if (.not. fde_cfg_no_sdft) then
          if (is_gga) then
             call get_gs_nonhermitian(s_b,                  &
                                      gs_b,                 &
                                      irep,                 &
                                      mat_dim,              &
                                      tb_dmat(1 + im_off),  &
                                      buffer,               &
                                      ao)
          else
             call get_s(s_b,                 &
                        irep,                &
                        mat_dim,             &
                        tb_dmat(1 + im_off), &
                        buffer,              &
                        ao)
          end if

          do k = k_start, 3

             u%s(k, im) = u%s(k, im) + (derv_total(1, d0200000) - derv_active(1, d0200000))*s_b(k)

             if (is_gga) then
                y_b(k)  = gn_0_pointwise(1)*gs_b(1, k) &
                        + gn_0_pointwise(2)*gs_b(2, k) &
                        + gn_0_pointwise(3)*gs_b(3, k)

                y_b_t(k) = gn_0_t_pointwise(1)*gs_b(1, k) &
                         + gn_0_t_pointwise(2)*gs_b(2, k) &
                         + gn_0_t_pointwise(3)*gs_b(3, k)

                u%s(k, im) = u%s(k, im) + (derv_total(1, d0101000)*y_b_t(k) - derv_active(1, d0101000)*y_b(k))

                u%gs(1:3, k, im) = u%gs(1:3, k, im) +       gn_0_t_pointwise(1:3)*s_b(k)*derv_total(1, d0101000)   &
                                                    +       gn_0_t_pointwise(1:3)*y_b_t(k)*derv_total(1, d0002000) &
                                                    + 2.0d0*gs_b(1:3, k)*derv_total(1, d0000100)                   &
                                                    -       gn_0_pointwise(1:3)*s_b(k)*derv_active(1, d0101000)    &
                                                    -       gn_0_pointwise(1:3)*y_b(k)*derv_active(1, d0002000)    &
                                                    - 2.0d0*gs_b(1:3, k)*derv_active(1, d0000100)
             end if
          end do

      end if
    end do

  end subroutine

  subroutine london_direct_kern(block_length,       &
                                mat_dim,            &
                                nz,                 &
                                ao,                 &
                                dmat_full,          &
                                buffer,             &
                                n_0,                &
                                n_0_t,              &
                                gn_0_pointwise,     &
                                gn_0_t_pointwise,   &
                                is_gga,             &
                                n,                  &
                                derv_active,        &
                                derv_total,         &
                                rx, ry, rz,         &
                                u)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: dmat_full(*)
    real(8), intent(in)    :: buffer(*)
    real(8), intent(in)    :: n_0(*)
    real(8), intent(in)    :: n_0_t(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: gn_0_t_pointwise(*)
    logical, intent(in)    :: is_gga
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv_active(max_block_length, 0:n)
    real(8), intent(in)    :: derv_total(max_block_length, 0:n)
    real(8), intent(in)    :: rx, ry, rz
!   ----------------------------------------------------------------------------
    real(8)                :: s_b_lao(3, 3)
    real(8)                :: gs_b_lao(3, 3, 3)
    real(8)                :: y_b, y_b_t

    integer                :: k_start, k, im
    type(fde_omega_prefactor) :: u
!   ----------------------------------------------------------------------------

    s_b_lao    = 0.0d0
    gs_b_lao   = 0.0d0

    k_start = 1
    if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

    if (is_gga) then

      call get_gs_lao(s_b_lao,     &
                      gs_b_lao,    &
                      (/0,0,0/),   &
                      mat_dim,     &
                      dmat_full,   &
                      ao,          &
                      rx, ry, rz, buffer)
    else
       call get_s_lao(s_b_lao,          &
                      (/0, 0, 0/),      &
                      mat_dim,          &
                      dmat_full,        &
                      ao,               &
                      rx, ry, rz, buffer)
    end if

    do im = 1, 3     !B_{x,y,z}
       do k = k_start, 3  !Sigma_{x,y,z}

          u%s_b_lao(k, im) = u%s_b_lao(k, im)  &
                           + (derv_total(1, d0200000) - derv_active(1, d0200000))*s_b_lao(k, im)
         
          if (is_gga) then

             y_b  = gn_0_pointwise(1)*gs_b_lao(1, k, im) &
                  + gn_0_pointwise(2)*gs_b_lao(2, k, im) &
                  + gn_0_pointwise(3)*gs_b_lao(3, k, im)

             y_b_t  = gn_0_t_pointwise(1)*gs_b_lao(1, k, im) &
                    + gn_0_t_pointwise(2)*gs_b_lao(2, k, im) &
                    + gn_0_t_pointwise(3)*gs_b_lao(3, k, im)

         
             u%s_b_lao(k,im) = u%s_b_lao(k,im) + (derv_total(1, d0101000)*y_b_t &
                                               - derv_active(1, d0101000)*y_b)
         
             u%gs_b_lao(1:3, k, im) = u%gs_b_lao(1:3, k, im) &
                     + s_b_lao(k, im)*( gn_0_t_pointwise(1:3)*derv_total(1, d0101000)    &
                                      - gn_0_pointwise(1:3)*derv_active(1, d0101000))    &
                     + gn_0_t_pointwise(1:3)*y_b_t*derv_total(1, d0002000)           &
                     - gn_0_pointwise(1:3)*y_b*derv_active(1, d0002000)              &
                     + 2.0d0*gs_b_lao(1:3, k, im)*(derv_total(1, d0000100) - derv_active(1, d0000100))
          end if
       end do
    end do


  end subroutine

  subroutine london_kern_coupling(ipoint,            &
                                  block_length,      &
                                  mat_dim,           &
                                  nz,                &
                                  ao,                &
                                  n_0_t,             &
                                  gn_0_t_pointwise,  &
                                  gnn_0_t,           &
                                  w,                 &
                                  is_gga,            &
                                  n,                 &
                                  derv_total,        &
                                  rx, ry, rz,        &
                                  u)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ipoint
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: n_0_t(*)
    real(8), intent(in)    :: gn_0_t_pointwise(*)
    real(8), intent(in)    :: gnn_0_t(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv_total(max_block_length, 0:n)
    real(8), intent(in)    :: rx, ry, rz
!   ----------------------------------------------------------------------------
    real(8)                ::  z_b
    real(8)                ::  n_b_lao(3)
    real(8)                ::  gn_b_lao(3, 3)
    real(8)                ::  s_b_lao(3, 3)
    real(8)                ::  gs_b_lao(3, 3, 3)
    real(8)                ::  gsgn_t_b_lao(3, 3)

    integer                :: k_start, k, a
    type(fde_omega_prefactor) :: u
!   ----------------------------------------------------------------------------
    n_b_lao    = 0.0d0
    gn_b_lao    = 0.0d0
    s_b_lao    = 0.0d0
    gs_b_lao   = 0.0d0
    gsgn_t_b_lao   = 0.0d0

    k_start = 1
    if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

!   charge density contributions - equal 0 when do_london_rhs_direct_coupling, so save time:
    if (.not. do_london_rhs_direct_coupling) then
      n_b_lao(1) = frozen_n_bx(ipoint)
      n_b_lao(2) = frozen_n_by(ipoint)
      n_b_lao(3) = frozen_n_bz(ipoint)
      if (is_gga) then
         gn_b_lao(1,1) = frozen_gxn_bx(ipoint)
         gn_b_lao(1,2) = frozen_gxn_bx(ipoint)
         gn_b_lao(1,3) = frozen_gxn_bx(ipoint)
         gn_b_lao(2,1) = frozen_gyn_bx(ipoint)
         gn_b_lao(2,2) = frozen_gyn_by(ipoint)
         gn_b_lao(2,3) = frozen_gyn_bz(ipoint)
         gn_b_lao(3,1) = frozen_gzn_bx(ipoint)
         gn_b_lao(3,2) = frozen_gzn_by(ipoint)
         gn_b_lao(3,3) = frozen_gzn_bz(ipoint)
      end if
      
!      write(*, *) 'test in coupling propgrad n, reorth, direct ', &
!                   do_london_rhs_reorth_coupling, do_london_rhs_direct_coupling, &
!                   rx, ry, rz, &
!                   n_b_lao(1), n_b_lao(2), n_b_lao(3)

      do a = 1, 3     !B_{x,y,z}
         u%n_b_lao(a) = u%n_b_lao(a) + derv_total(1, d2000000)*n_b_lao(a)
      
         if (is_gga) then
            z_b = 2.0d0*(gn_0_t_pointwise(1)*gn_b_lao(1, a) &
                       + gn_0_t_pointwise(2)*gn_b_lao(2, a) &
                       + gn_0_t_pointwise(3)*gn_b_lao(3, a))
      
            u%n_b_lao(a) = u%n_b_lao(a) + derv_total(1, d1010000)*z_b
      
            u%gn_b_lao(1:3, a) = u%gn_b_lao(1:3, a) &
                       + 2.0d0*gn_0_t_pointwise(1:3)*n_b_lao(a)*derv_total(1, d1010000) &
                       + 2.0d0*gn_0_t_pointwise(1:3)*z_b*derv_total(1, d0020000) &
                       + 2.0d0*gn_b_lao(1:3, a)*derv_total(1, d0010000)
         end if
      end do
    end if

!   spin density contributions

    if (.not. fde_cfg_no_sdft) then
       s_b_lao(1,1) = frozen_sx_bx(ipoint)
       s_b_lao(1,2) = frozen_sx_by(ipoint)
       s_b_lao(1,3) = frozen_sx_bz(ipoint)
       s_b_lao(2,1) = frozen_sy_bx(ipoint)
       s_b_lao(2,2) = frozen_sy_by(ipoint)
       s_b_lao(2,3) = frozen_sy_bz(ipoint)
       s_b_lao(3,1) = frozen_sz_bx(ipoint)
       s_b_lao(3,2) = frozen_sz_by(ipoint)
       s_b_lao(3,3) = frozen_sz_bz(ipoint)
       if (is_gga) then
          gs_b_lao(1,1,1) = frozen_gxsx_bx(ipoint)
          gs_b_lao(1,1,2) = frozen_gxsx_by(ipoint)
          gs_b_lao(1,1,3) = frozen_gxsx_bz(ipoint)
          gs_b_lao(1,2,1) = frozen_gxsy_bx(ipoint)
          gs_b_lao(1,2,2) = frozen_gxsy_by(ipoint)
          gs_b_lao(1,2,3) = frozen_gxsy_bz(ipoint)
          gs_b_lao(1,3,1) = frozen_gxsz_bx(ipoint)
          gs_b_lao(1,3,2) = frozen_gxsz_by(ipoint)
          gs_b_lao(1,3,3) = frozen_gxsz_bz(ipoint)
          gs_b_lao(2,1,1) = frozen_gysx_bx(ipoint)
          gs_b_lao(2,1,2) = frozen_gysx_by(ipoint)
          gs_b_lao(2,1,3) = frozen_gysx_bz(ipoint)
          gs_b_lao(2,2,1) = frozen_gysy_bx(ipoint)
          gs_b_lao(2,2,2) = frozen_gysy_by(ipoint)
          gs_b_lao(2,2,3) = frozen_gysy_bz(ipoint)
          gs_b_lao(2,3,1) = frozen_gysz_bx(ipoint)
          gs_b_lao(2,3,2) = frozen_gysz_by(ipoint)
          gs_b_lao(2,3,3) = frozen_gysz_bz(ipoint)
          gs_b_lao(3,1,1) = frozen_gzsx_bx(ipoint)
          gs_b_lao(3,1,2) = frozen_gzsx_by(ipoint)
          gs_b_lao(3,1,3) = frozen_gzsx_bz(ipoint)
          gs_b_lao(3,2,1) = frozen_gzsy_bx(ipoint)
          gs_b_lao(3,2,2) = frozen_gzsy_by(ipoint)
          gs_b_lao(3,2,3) = frozen_gzsy_bz(ipoint)
          gs_b_lao(3,3,1) = frozen_gzsz_bx(ipoint)
          gs_b_lao(3,3,2) = frozen_gzsz_by(ipoint)
          gs_b_lao(3,3,3) = frozen_gzsz_bz(ipoint)
       end if

!       write(*, *) 'test in coupling propgrad s, reorth, direct ', &
!                    do_london_rhs_reorth_coupling, do_london_rhs_direct_coupling, &
!                    rx, ry, rz, &
!                    s_b_lao(1,1), s_b_lao(1,2), s_b_lao(1,3)

       do a = 1, 3     !B_{x,y,z}
          do k = k_start, 3  !Sigma_{x,y,z}
       
             u%s_b_lao(k, a) = u%s_b_lao(k, a)  &
                             + derv_total(1, d0200000)*s_b_lao(k, a)
            
             if (is_gga) then
            
                gsgn_t_b_lao(k, a) = gs_b_lao(1,k,a)*gn_0_t_pointwise(1) &
                                   + gs_b_lao(2,k,a)*gn_0_t_pointwise(2) &
                                   + gs_b_lao(3,k,a)*gn_0_t_pointwise(3)
       
                u%s_b_lao(k, a) = u%s_b_lao(k, a) + derv_total(1, d0101000)*gsgn_t_b_lao(k, a)
           
                u%gs_b_lao(1:3, k, a) = u%gs_b_lao(1:3, k, a) &
                            + gn_0_t_pointwise(1:3)*s_b_lao(k, a)*derv_total(1, d0101000)      &
                            + gn_0_t_pointwise(1:3)*gsgn_t_b_lao(k, a)*derv_total(1, d0002000) &
                            + 2.0d0*gs_b_lao(1:3, k, a)*derv_total(1, d0000100)
       
             end if
          end do
       end do
    end if


  end subroutine

  subroutine london_susc2el_coupling(ipoint,            &
                                     block_length,      &
                                     mat_dim,           &
                                     nz,                &
                                     ao,                &
                                     n_0_t,             &
                                     gn_0_t_pointwise,  &
                                     gnn_0_t,           &
                                     w,                 &
                                     is_gga,            &
                                     n,                 &
                                     derv_total,        &
                                     rx, ry, rz,        &
                                     buffer)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ipoint
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: n_0_t(*)
    real(8), intent(in)    :: gn_0_t_pointwise(*)
    real(8), intent(in)    :: gnn_0_t(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv_total(max_block_length, 0:n)
    real(8), intent(in)    :: rx, ry, rz
    real(8), intent(in)    :: buffer(*)
!   ----------------------------------------------------------------------------
    real(8)                ::  s_b_lao(3, 3), s_b_lao_active(3, 3)
    real(8)                ::  gs_b_lao(3, 3, 3), gs_b_lao_active(3, 3, 3)
    real(8)                ::  s, d_s, d_sns, d_nsns, d_ss, n_sum
    real(8)                ::  s_1(3), gs_2(3)
    real(8)                ::  gngs(3, 3), gsgs(3, 3, 3), gngs_active(3, 3)
    real(8)                ::  n_bb(3,3)

    integer                :: k_start, k, a
    integer                :: ixyz, jxyz, a1, a2
!   ----------------------------------------------------------------------------
    s_b_lao    = 0.0d0
    gs_b_lao   = 0.0d0
    n_bb = 0.0d0

    k_start = 1
    if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

!   perturbed density of frozen subsystem collected in s_b_lao, gs_b_lao
    s_b_lao(1,1) = frozen_sx_bx(ipoint)
    s_b_lao(1,2) = frozen_sx_by(ipoint)
    s_b_lao(1,3) = frozen_sx_bz(ipoint)
    s_b_lao(2,1) = frozen_sy_bx(ipoint)
    s_b_lao(2,2) = frozen_sy_by(ipoint)
    s_b_lao(2,3) = frozen_sy_bz(ipoint)
    s_b_lao(3,1) = frozen_sz_bx(ipoint)
    s_b_lao(3,2) = frozen_sz_by(ipoint)
    s_b_lao(3,3) = frozen_sz_bz(ipoint)
    if (is_gga) then
       gs_b_lao(1,1,1) = frozen_gxsx_bx(ipoint)
       gs_b_lao(1,1,2) = frozen_gxsx_by(ipoint)
       gs_b_lao(1,1,3) = frozen_gxsx_bz(ipoint)
       gs_b_lao(1,2,1) = frozen_gxsy_bx(ipoint)
       gs_b_lao(1,2,2) = frozen_gxsy_by(ipoint)
       gs_b_lao(1,2,3) = frozen_gxsy_bz(ipoint)
       gs_b_lao(1,3,1) = frozen_gxsz_bx(ipoint)
       gs_b_lao(1,3,2) = frozen_gxsz_by(ipoint)
       gs_b_lao(1,3,3) = frozen_gxsz_bz(ipoint)
       gs_b_lao(2,1,1) = frozen_gysx_bx(ipoint)
       gs_b_lao(2,1,2) = frozen_gysx_by(ipoint)
       gs_b_lao(2,1,3) = frozen_gysx_bz(ipoint)
       gs_b_lao(2,2,1) = frozen_gysy_bx(ipoint)
       gs_b_lao(2,2,2) = frozen_gysy_by(ipoint)
       gs_b_lao(2,2,3) = frozen_gysy_bz(ipoint)
       gs_b_lao(2,3,1) = frozen_gysz_bx(ipoint)
       gs_b_lao(2,3,2) = frozen_gysz_by(ipoint)
       gs_b_lao(2,3,3) = frozen_gysz_bz(ipoint)
       gs_b_lao(3,1,1) = frozen_gzsx_bx(ipoint)
       gs_b_lao(3,1,2) = frozen_gzsx_by(ipoint)
       gs_b_lao(3,1,3) = frozen_gzsx_bz(ipoint)
       gs_b_lao(3,2,1) = frozen_gzsy_bx(ipoint)
       gs_b_lao(3,2,2) = frozen_gzsy_by(ipoint)
       gs_b_lao(3,2,3) = frozen_gzsy_bz(ipoint)
       gs_b_lao(3,3,1) = frozen_gzsz_bx(ipoint)
       gs_b_lao(3,3,2) = frozen_gzsz_by(ipoint)
       gs_b_lao(3,3,3) = frozen_gzsz_bz(ipoint)
    end if

!   perturbed density of active subsystem collected in s_b_lao_active, gs_b_lao_active
    s_b_lao_active = 0.0d0
    gs_b_lao_active = 0.0d0
    if (is_gga) then
       call get_gs_lao(s_b_lao_active,   &
                       gs_b_lao_active,  &
                       (/0, 0, 0/),      &
                       mat_dim,          &
                       dmat_0,           &
                       ao,               &
                       rx, ry, rz,       &
                       buffer)
    else
       call get_s_lao(s_b_lao_active,    &
                      (/0, 0, 0/),       &
                      mat_dim,           &
                      dmat_0,            &
                      ao,                &
                      rx, ry, rz,        &
                      buffer)
    end if

    if (is_gga) then
      d_s    = derv_total(1, d0200000)
      d_sns  = derv_total(1, d0101000)
      d_nsns = derv_total(1, d0002000)
      d_ss   = derv_total(1, d0000100)*2.0d0

      gngs = 0.0d0
      gngs_active = 0.0d0
      gsgs = 0.0d0
      do jxyz = 1, 3
        do a1 = 1, 3
          gngs(jxyz, a1) = gngs(jxyz, a1) &
                         + gn_0_t_pointwise(1)*gs_b_lao(1, jxyz, a1) &
                         + gn_0_t_pointwise(2)*gs_b_lao(2, jxyz, a1) &
                         + gn_0_t_pointwise(3)*gs_b_lao(3, jxyz, a1) 
      
          gngs_active(jxyz, a1) = gngs_active(jxyz, a1) &
                         + gn_0_t_pointwise(1)*gs_b_lao_active(1, jxyz, a1) &
                         + gn_0_t_pointwise(2)*gs_b_lao_active(2, jxyz, a1) &
                         + gn_0_t_pointwise(3)*gs_b_lao_active(3, jxyz, a1) 
           do a2 = 1, 3
            gsgs(jxyz, a1, a2) = gsgs(jxyz, a1, a2) &
                               + gs_b_lao_active(1, jxyz, a1)*gs_b_lao(1, jxyz, a2) &
                               + gs_b_lao_active(2, jxyz, a1)*gs_b_lao(2, jxyz, a2) &
                               + gs_b_lao_active(3, jxyz, a1)*gs_b_lao(3, jxyz, a2)
          end do
        end do
      end do

      do jxyz = 1, 3
        do a1 = 1, 3
          do a2 = 1, 3
            n_bb(a2, a1) = n_bb(a2, a1) &
!                          lda
                         + d_s*s_b_lao_active(jxyz,a1)*s_b_lao(jxyz,a2) &
!                          gradient terms
                         + d_sns*(gngs(jxyz,a1)*s_b_lao_active(jxyz,a2) + gngs_active(jxyz,a2)*s_b_lao(jxyz,a1)) &
                         + d_nsns*(gngs_active(jxyz,a1)*gngs(jxyz,a2)) &
                         + d_ss*gsgs(jxyz, a1, a2)
!g           write(*, *) 'test8', s_b_lao_active(jxyz,a1), s_b_lao(jxyz,a2)
          end do
        end do
      end do

    else
      d_s = derv_total(1, d0200000)
      do jxyz = 1, 3
         do a1 = 1, 3
           do a2 = 1, 3
             n_bb(a2, a1) = n_bb(a2, a1) + d_s*s_b_lao_active(jxyz,a1)*s_b_lao(jxyz,a2)
!g           write(*, *) 'test9', s_b_lao_active(jxyz,a1), s_b_lao(jxyz,a2)
           end do
         end do
      end do
    end if
!g    write(*, *) 'test10', n_bb(1,1), n_bb(2,2), n_bb(3,3)

!   sum to susc2el_xc_integrated_sdft
    do ixyz = 1, 3
      do jxyz = 1, 3
        susc2el_fde_integrated_sdft(jxyz, ixyz) = susc2el_fde_integrated_sdft(jxyz, ixyz) &
                                          + n_bb(jxyz, ixyz)
      end do
    end do

  end subroutine


  subroutine london_suscreo_coupling(ipoint,            &
                                     block_length,      &
                                     mat_dim,           &
                                     nz,                &
                                     dmat_T,            &
                                     ao,                &
                                     n_0,             &
                                     gn_0_pointwise,  &
                                     w,                 &
                                     is_gga,            &
                                     n,                 &
                                     derv_total,        &
                                     rx, ry, rz,        &
                                     n_bb,buffer)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ipoint
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(in)    :: dmat_T(*)
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: n_0(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv_total(max_block_length, 0:n)
    real(8), intent(in)    :: rx, ry, rz
    real(8), intent(inout) :: n_bb(3,3)
    real(8), intent(in)    :: buffer(*)
!   ----------------------------------------------------------------------------
    real(8)                ::  d_1, d_2, d_3, d_4

    real(8)                ::  n_b_lao(3)
    real(8)                ::  gn_b_lao(3, 3)
    real(8)                ::  n_active(3), gn_active(3,3)
    real(8)                ::  gngn_active(3), gngn_frozen(3)
    real(8)                ::  gngn_active_frozen(3,3)

    real(8)                ::  s_b_lao(3, 3)
    real(8)                ::  gs_b_lao(3, 3, 3)
    real(8)                ::  s1(3,3),gs1(3,3,3),gsgn_b_lao(3,3),s2(3,3),gs2(3,3,3)
    real(8)                ::  s_active(3,3), gs_active(3,3,3)
    real(8)                ::  gngs_active(3,3), gngs_frozen(3,3)
    real(8)                ::  gsgs_active_frozen(3,3)

    integer                :: k_start, k, a, ioff
    integer                :: ixyz, jxyz, a1, a2
!   ----------------------------------------------------------------------------
    n_b_lao  = 0.0d0
    gn_b_lao = 0.0d0
    s_b_lao    = 0.0d0
    gs_b_lao   = 0.0d0

!   frozen_* variables read from file correspond to the sum
!   of direct (London phase) and reorthonormalization contributions
!   (LAO + T)^alpha
    n_b_lao(1) = frozen_n_bx(ipoint)
    n_b_lao(2) = frozen_n_by(ipoint)
    n_b_lao(3) = frozen_n_bz(ipoint)
    if (is_gga) then
      gn_b_lao(1, 1) = frozen_gxn_bx(ipoint)
      gn_b_lao(1, 2) = frozen_gxn_by(ipoint)
      gn_b_lao(1, 3) = frozen_gxn_bz(ipoint)
      gn_b_lao(2, 1) = frozen_gyn_bx(ipoint)
      gn_b_lao(2, 2) = frozen_gyn_by(ipoint)
      gn_b_lao(2, 3) = frozen_gyn_bz(ipoint)
      gn_b_lao(3, 1) = frozen_gzn_bx(ipoint)
      gn_b_lao(3, 2) = frozen_gzn_by(ipoint)
      gn_b_lao(3, 3) = frozen_gzn_bz(ipoint)
    end if

!   then we need active densities also to correspond to (LAO+T)^beta
!   the real part of LAO-density is 0, so take only the T-density
    if (is_gga) then
       do a1 = 1,3
         ioff = (mat_dim*mat_dim*nz)*(a1-1) + 1
         call get_gn_nonhermitian(n_active(a1),      &
                                  gn_active(1:3,a1), &
                                  0,                 &
                                  mat_dim,           &
                                  dmat_T(ioff),      &
                                  buffer,            &
                                  ao)
      end do
    else
      do a1 = 1,3
        ioff = (mat_dim*mat_dim*nz)*(a1-1) + 1
        call get_n(n_active(a1),    &
                   0,               &
                   mat_dim,         &
                   dmat_T(ioff),    &
                   buffer,          &
                   ao)
      end do
    end if

!   combine all together:
!   =====================
!   lda
    d_1    = derv_total(1, d2000000)
    do a1 = 1, 3
      do a2 = 1, 3
        n_bb(a2, a1) = n_bb(a2, a1) &
                     + d_1*n_active(a1)*n_b_lao(a2)
      end do
    end do

!   gga
    if (is_gga) then
      d_2    = derv_total(1, d1010000)
      d_3    = derv_total(1, d0020000)
      d_4    = derv_total(1, d0010000)*2.0d0
    
      do a1 = 1, 3! B_{x,y,z}
        gngn_active(a1) = gn_0_pointwise(1)*gn_active(1, a1) &
                        + gn_0_pointwise(2)*gn_active(2, a1) &
                        + gn_0_pointwise(3)*gn_active(3, a1)

        gngn_frozen(a1) = gn_0_pointwise(1)*gn_b_lao(1, a1)  &
                        + gn_0_pointwise(2)*gn_b_lao(2, a1)  &
                        + gn_0_pointwise(3)*gn_b_lao(3, a1)

        do a2 = 1, 3! B_{x,y,z}
          gngn_active_frozen(a2, a1) = gn_active(1, a2)*gn_b_lao(1, a1)  &
                                     + gn_active(2, a2)*gn_b_lao(2, a1)  &
                                     + gn_active(3, a2)*gn_b_lao(3, a1)
        end do
      end do

      do a1 = 1, 3
        do a2 = 1, 3
          n_bb(a2, a1) = n_bb(a2, a1) &
                       + 2.0d0*d_2*(gngn_active(a1)*n_b_lao(a2)+ gngn_frozen(a1)*n_active(a2)) &
                       + 4.0d0*d_3*gngn_active(a1)*gngn_frozen(a1) &
                       + 2.0d0*d_4*gngn_active_frozen(a2, a1)
        end do
      end do
    end if
      
!   sum to susc2el_xc_integrated
    do ixyz = 1, 3
      do jxyz = 1, 3
        susc2el_fde_integrated(jxyz, ixyz) = susc2el_fde_integrated(jxyz, ixyz) &
                                          + n_bb(jxyz, ixyz)
      end do
    end do


!   sdft
    if (.not.fde_cfg_no_sdft) then
     
      k_start = 1
      if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation
      
      s_b_lao(1,1) = frozen_sx_bx(ipoint)
      s_b_lao(1,2) = frozen_sx_by(ipoint)
      s_b_lao(1,3) = frozen_sx_bz(ipoint)
      s_b_lao(2,1) = frozen_sy_bx(ipoint)
      s_b_lao(2,2) = frozen_sy_by(ipoint)
      s_b_lao(2,3) = frozen_sy_bz(ipoint)
      s_b_lao(3,1) = frozen_sz_bx(ipoint)
      s_b_lao(3,2) = frozen_sz_by(ipoint)
      s_b_lao(3,3) = frozen_sz_bz(ipoint)
      if (is_gga) then
         gs_b_lao(1,1,1) = frozen_gxsx_bx(ipoint)
         gs_b_lao(1,1,2) = frozen_gxsx_by(ipoint)
         gs_b_lao(1,1,3) = frozen_gxsx_bz(ipoint)
         gs_b_lao(1,2,1) = frozen_gxsy_bx(ipoint)
         gs_b_lao(1,2,2) = frozen_gxsy_by(ipoint)
         gs_b_lao(1,2,3) = frozen_gxsy_bz(ipoint)
         gs_b_lao(1,3,1) = frozen_gxsz_bx(ipoint)
         gs_b_lao(1,3,2) = frozen_gxsz_by(ipoint)
         gs_b_lao(1,3,3) = frozen_gxsz_bz(ipoint)
         gs_b_lao(2,1,1) = frozen_gysx_bx(ipoint)
         gs_b_lao(2,1,2) = frozen_gysx_by(ipoint)
         gs_b_lao(2,1,3) = frozen_gysx_bz(ipoint)
         gs_b_lao(2,2,1) = frozen_gysy_bx(ipoint)
         gs_b_lao(2,2,2) = frozen_gysy_by(ipoint)
         gs_b_lao(2,2,3) = frozen_gysy_bz(ipoint)
         gs_b_lao(2,3,1) = frozen_gysz_bx(ipoint)
         gs_b_lao(2,3,2) = frozen_gysz_by(ipoint)
         gs_b_lao(2,3,3) = frozen_gysz_bz(ipoint)
         gs_b_lao(3,1,1) = frozen_gzsx_bx(ipoint)
         gs_b_lao(3,1,2) = frozen_gzsx_by(ipoint)
         gs_b_lao(3,1,3) = frozen_gzsx_bz(ipoint)
         gs_b_lao(3,2,1) = frozen_gzsy_bx(ipoint)
         gs_b_lao(3,2,2) = frozen_gzsy_by(ipoint)
         gs_b_lao(3,2,3) = frozen_gzsy_bz(ipoint)
         gs_b_lao(3,3,1) = frozen_gzsz_bx(ipoint)
         gs_b_lao(3,3,2) = frozen_gzsz_by(ipoint)
         gs_b_lao(3,3,3) = frozen_gzsz_bz(ipoint)
      end if
      

!     then we need active densities also to correspond to (LAO+T)^beta
!     first get s_lao and gs_lao:
      if (is_gga) then
        call get_gs_lao(s1,     &
                        gs1,    &
                        (/1,1,1/),   &
                        mat_dim,     &
                        dmat_0,      &
                        ao,          &
                        rx, ry, rz, buffer)
      else
         call get_s_lao(s1,          &
                        (/0, 0, 0/),      &
                        mat_dim,          &
                        dmat_0,           &
                        ao,               &
                        rx, ry, rz, buffer)
      end if
      
!     then get s_T, gs_T:
      if (is_gga) then
        do a1 = 1, 3
          ioff = mat_dim*mat_dim*nz*(a1-1)+1
          call get_gs_nonhermitian(s2(1:3,a1),          &
                                   gs2(1:3,1:3,a1),     &
                                   0,                   &
                                   mat_dim,             &
                                   dmat_T(ioff),        &
                                   buffer,              &
                                   ao)
        end do
      else
        do a1 = 1, 3
          ioff = mat_dim*mat_dim*nz*(a1-1)+1
          call get_s(s2(1:3,a1),     &
                     0,              &
                     mat_dim,        &
                     dmat_T(ioff),   &
                     buffer,         &
                     ao)
        end do
      end if

!     rename & sum
      do a1 = 1, 3
        do a2 = 1, 3
          s_active(a2, a1) = s1(a2, a1) + s2(a2, a1)  !s_active(Sigma, B)
          do k = 1, 3
            gs_active(k, a2, a1) = gs1(k, a2, a1) + gs2(k, a2, a1) !gs_active(dR, Sigma, B)
          end do
        end do
      end do
      
!     combine all together:
!     =====================
!     lda
      d_1    = derv_total(1, d0200000)
      do a1 = 1, 3
        do a2 = 1, 3
          do k = 1, 3
            n_bb(a2, a1) = n_bb(a2, a1) &
                         + d_1*s_active(k, a1)*s_b_lao(k, a2)
          end do
        end do
      end do
      
!     gga
      if (is_gga) then
        d_2    = derv_total(1, d0101000)
        d_3    = derv_total(1, d0002000)
        d_4    = derv_total(1, d0000100)*2.0d0
      
        do a1 = 1, 3! B_{x,y,z}
          do k = 1, 3! Sigma
            gngs_active(k, a1) = gn_0_pointwise(1)*gs_active(1, k, a1) &
                               + gn_0_pointwise(2)*gs_active(2, k, a1) &
                               + gn_0_pointwise(3)*gs_active(3, k, a1)
            
            gngs_frozen(k, a1) = gn_0_pointwise(1)*gs_b_lao(1, k, a1)  &
                               + gn_0_pointwise(2)*gs_b_lao(2, k, a1)  &
                               + gn_0_pointwise(3)*gs_b_lao(3, k, a1)
            
            do a2 = 1, 3! B_{x,y,z}
              gsgs_active_frozen(a2, a1) = gs_active(1, k, a2)*gs_b_lao(1, k, a1)  &
                                         + gs_active(2, k, a2)*gs_b_lao(2, k, a1)  &
                                         + gs_active(3, k, a2)*gs_b_lao(3, k, a1)
            end do
          end do
        end do

        do a1 = 1, 3
          do a2 = 1, 3
            n_bb(a2, a1) = n_bb(a2, a1) &
                         + 2.0d0*d_4*gsgs_active_frozen(a2, a1)
            do k = 1, 3
              n_bb(a2, a1) = n_bb(a2, a1) &
                           + d_2*(gngs_active(k,a1)*s_b_lao(k,a2)+ gngs_frozen(k,a1)*s_active(k,a2)) &
                           + d_3*gngs_active(k,a1)*gngs_frozen(k,a1)
            end do
          end do
        end do
     
      end if

!     sum to susc2el_xc_integrated_sdft
      do ixyz = 1, 3
        do jxyz = 1, 3
          susc2el_fde_integrated_sdft(jxyz, ixyz) = susc2el_fde_integrated_sdft(jxyz, ixyz) &
                                            + n_bb(jxyz, ixyz)
        end do
      end do

    end if ! sdft

  end subroutine



! ----------------------------------------------------------------------------
   subroutine loop_over_batches_elpot()
! ----------------------------------------------------------------------------
     real(8) :: dummy
     integer :: idummy
     integer :: ibatch
     integer :: i
     real(8), allocatable :: mpi_buffer(:)

     if (fde_mpi_is_master()) then
        nr_points_batch = active_grid%npoints
     end if
#ifdef VAR_MPI
     if (parallel_fde) then
        call fde_mpi_bcast(nr_points_batch)
     end if
#endif
     ibatch          = 0
     nr_points_total = 0
     nr_points_used  = 0

!    loop over batches
     do
        if (nr_points_batch < 0) exit

        allocate(rx(nr_points_batch))
        allocate(ry(nr_points_batch))
        allocate(rz(nr_points_batch))
        allocate(rw(nr_points_batch))

        allocate(active_elpot(nr_points_batch))
        allocate(active_nucpot(nr_points_batch))
        allocate(active_n(nr_points_batch))
        allocate(active_gnx(nr_points_batch))
        allocate(active_gny(nr_points_batch))
        allocate(active_gnz(nr_points_batch))
        allocate(active_hnxx(nr_points_batch))
        allocate(active_hnxy(nr_points_batch))
        allocate(active_hnxz(nr_points_batch))
        allocate(active_hnyy(nr_points_batch))
        allocate(active_hnyz(nr_points_batch))
        allocate(active_hnzz(nr_points_batch))

        if (fde_mpi_is_master()) then
           rx = active_grid%r(1,:)
           ry = active_grid%r(2,:)
           rz = active_grid%r(3,:)
           rw = active_grid%w
        end if

        ibatch = ibatch + 1
        nr_points_total = nr_points_total + nr_points_batch
        nr_points_on_this_proc = nr_points_batch

#ifdef VAR_MPI
        if (parallel_fde) then
           call fde_mpi_distribute_points(rx, ry, rz, rw, nr_points_batch, nr_points_on_this_proc)
        end if
#endif
        call get_electrostatic_potential_and_density()

        deallocate(rx)
        deallocate(ry)
        deallocate(rz)
        deallocate(rw)

#ifdef VAR_MPI
! we send the points in the nodes to a final vector in master
        if (parallel_fde) then
           call fde_mpi_gatherv(active_elpot,  active_gridfunction%elpot,   nr_points_on_this_proc)
           call fde_mpi_gatherv(active_nucpot, active_gridfunction%nucpot,  nr_points_on_this_proc)
           call fde_mpi_gatherv(active_n,      active_gridfunction%n,       nr_points_on_this_proc)
           call fde_mpi_gatherv(active_gnx,    active_gridfunction%gn(1,:), nr_points_on_this_proc)
           call fde_mpi_gatherv(active_gny,    active_gridfunction%gn(2,:), nr_points_on_this_proc)
           call fde_mpi_gatherv(active_gnz,    active_gridfunction%gn(3,:), nr_points_on_this_proc)

           if (fde_mpi_is_master()) then

! aspg, 19/10/2015
! here we scale the density, density gradient etc since for the
! calculation of the electrostatic potential we needed the unscaled 
! density matrix

               active_gridfunction%n(:)    = 2.0d0 * active_gridfunction%n(:)

               active_gridfunction%gn(1,:) = 2.0d0 * active_gridfunction%gn(1,:)
               active_gridfunction%gn(2,:) = 2.0d0 * active_gridfunction%gn(2,:)
               active_gridfunction%gn(3,:) = 2.0d0 * active_gridfunction%gn(3,:)

               active_gridfunction%hn(1,:) = 2.0d0 * active_gridfunction%hn(1,:)
               active_gridfunction%hn(2,:) = 2.0d0 * active_gridfunction%hn(2,:)
               active_gridfunction%hn(3,:) = 2.0d0 * active_gridfunction%hn(3,:)
               active_gridfunction%hn(4,:) = 2.0d0 * active_gridfunction%hn(4,:)
               active_gridfunction%hn(5,:) = 2.0d0 * active_gridfunction%hn(5,:)
               active_gridfunction%hn(6,:) = 2.0d0 * active_gridfunction%hn(6,:)
            end if
        else ! parallel_fde
#endif
!
! in the serial version, just update the gridfunction with the contents of
! the arrays 
         active_gridfunction%elpot(:)  = active_elpot(:)
         active_gridfunction%nucpot(:) = active_nucpot(:)

         active_gridfunction%n(:)    = 2.0d0 * active_n(:)

         active_gridfunction%gn(1,:) = 2.0d0 * active_gnx(:)
         active_gridfunction%gn(2,:) = 2.0d0 * active_gny(:)
         active_gridfunction%gn(3,:) = 2.0d0 * active_gnz(:)

         active_gridfunction%hn(1,:) = 2.0d0 * active_hnxx(:)
         active_gridfunction%hn(2,:) = 2.0d0 * active_hnxy(:)
         active_gridfunction%hn(3,:) = 2.0d0 * active_hnxz(:)
         active_gridfunction%hn(4,:) = 2.0d0 * active_hnyy(:)
         active_gridfunction%hn(5,:) = 2.0d0 * active_hnyz(:)
         active_gridfunction%hn(6,:) = 2.0d0 * active_hnzz(:)
 
#ifdef VAR_MPI
        end if ! parallel_fde
#endif
        deallocate(active_elpot)
        deallocate(active_nucpot)
        deallocate(active_n)
        deallocate(active_gnx)
        deallocate(active_gny)
        deallocate(active_gnz)
        deallocate(active_hnxx)
        deallocate(active_hnxy)
        deallocate(active_hnxz)
        deallocate(active_hnyy)
        deallocate(active_hnyz)
        deallocate(active_hnzz)


! don't try to read another batch
        goto 999

!     end loop over batches
     end do
 999 continue

     end subroutine loop_over_batches_elpot


! ----------------------------------------------------------------------------
   subroutine loop_over_batches_pertden()
! ----------------------------------------------------------------------------
     real(8) :: dummy
     integer :: idummy
     integer :: ibatch
     integer :: i
     real(8), allocatable :: mpi_buffer(:)

     ibatch          = 0
     nr_points_total = 0
     nr_points_used  = 0

        if (fde_mpi_is_master()) then
           nr_points_batch = active_grid%npoints
        end if
#ifdef VAR_MPI
        if (parallel_fde) then
          call fde_mpi_bcast(nr_points_batch)
        end if
#endif

!    loop over batches
     do
        if (nr_points_batch < 0) exit

        allocate(rx(nr_points_batch))
        allocate(ry(nr_points_batch))
        allocate(rz(nr_points_batch))
        allocate(rw(nr_points_batch))

        allocate(active_n_bx(nr_points_batch))
        allocate(active_n_by(nr_points_batch))
        allocate(active_n_bz(nr_points_batch))

        if (use_gga) then
           allocate(active_gxn_bx(nr_points_batch))
           allocate(active_gxn_by(nr_points_batch))
           allocate(active_gxn_bz(nr_points_batch))
           allocate(active_gyn_bx(nr_points_batch))
           allocate(active_gyn_by(nr_points_batch))
           allocate(active_gyn_bz(nr_points_batch))
           allocate(active_gzn_bx(nr_points_batch))
           allocate(active_gzn_by(nr_points_batch))
           allocate(active_gzn_bz(nr_points_batch))
        end if
        if (.not. fde_cfg_no_sdft) then
           allocate(active_sx_bx(nr_points_batch))
           allocate(active_sx_by(nr_points_batch))
           allocate(active_sx_bz(nr_points_batch))
           allocate(active_sy_bx(nr_points_batch))
           allocate(active_sy_by(nr_points_batch))
           allocate(active_sy_bz(nr_points_batch))
           allocate(active_sz_bx(nr_points_batch))
           allocate(active_sz_by(nr_points_batch))
           allocate(active_sz_bz(nr_points_batch))
           if (use_gga) then
              allocate(active_gxsx_bx(nr_points_batch))
              allocate(active_gxsx_by(nr_points_batch))
              allocate(active_gxsx_bz(nr_points_batch))
              allocate(active_gxsy_bx(nr_points_batch))
              allocate(active_gxsy_by(nr_points_batch))
              allocate(active_gxsy_bz(nr_points_batch))
              allocate(active_gxsz_bx(nr_points_batch))
              allocate(active_gxsz_by(nr_points_batch))
              allocate(active_gxsz_bz(nr_points_batch))
              allocate(active_gysx_bx(nr_points_batch))
              allocate(active_gysx_by(nr_points_batch))
              allocate(active_gysx_bz(nr_points_batch))
              allocate(active_gysy_bx(nr_points_batch))
              allocate(active_gysy_by(nr_points_batch))
              allocate(active_gysy_bz(nr_points_batch))
              allocate(active_gysz_bx(nr_points_batch))
              allocate(active_gysz_by(nr_points_batch))
              allocate(active_gysz_bz(nr_points_batch))
              allocate(active_gzsx_bx(nr_points_batch))
              allocate(active_gzsx_by(nr_points_batch))
              allocate(active_gzsx_bz(nr_points_batch))
              allocate(active_gzsy_bx(nr_points_batch))
              allocate(active_gzsy_by(nr_points_batch))
              allocate(active_gzsy_bz(nr_points_batch))
              allocate(active_gzsz_bx(nr_points_batch))
              allocate(active_gzsz_by(nr_points_batch))
              allocate(active_gzsz_bz(nr_points_batch))
           end if
        end if

        if (fde_mpi_is_master()) then
           rx = active_grid%r(1,:)
           ry = active_grid%r(2,:)
           rz = active_grid%r(3,:)
           rw = active_grid%w(:)
        end if

        ibatch = ibatch + 1
        nr_points_total = nr_points_total + nr_points_batch
        nr_points_on_this_proc = nr_points_batch

#ifdef VAR_MPI
        if (parallel_fde) then
           call fde_mpi_distribute_points(rx, ry, rz, rw, nr_points_batch, nr_points_on_this_proc)
        end if
#endif

        call get_perturbed_density()

        deallocate(rx)
        deallocate(ry)
        deallocate(rz)
        deallocate(rw)

#ifdef VAR_MPI
        if (parallel_fde) then
           if (do_london_rhs_reorth_coupling) then
             call fde_mpi_gatherv(active_n_bx, active_gridfunction%n_b_reorth(1,:), nr_points_on_this_proc)
             call fde_mpi_gatherv(active_n_by, active_gridfunction%n_b_reorth(2,:), nr_points_on_this_proc)
             call fde_mpi_gatherv(active_n_bz, active_gridfunction%n_b_reorth(3,:), nr_points_on_this_proc)
             
             if (use_gga) then
                call fde_mpi_gatherv(active_gxn_bx, active_gridfunction%gn_b_reorth(1,1,:), nr_points_on_this_proc)
                call fde_mpi_gatherv(active_gxn_by, active_gridfunction%gn_b_reorth(1,2,:), nr_points_on_this_proc)
                call fde_mpi_gatherv(active_gxn_bz, active_gridfunction%gn_b_reorth(1,3,:), nr_points_on_this_proc)
                call fde_mpi_gatherv(active_gyn_bx, active_gridfunction%gn_b_reorth(2,1,:), nr_points_on_this_proc)
                call fde_mpi_gatherv(active_gyn_by, active_gridfunction%gn_b_reorth(2,2,:), nr_points_on_this_proc)
                call fde_mpi_gatherv(active_gyn_bz, active_gridfunction%gn_b_reorth(2,3,:), nr_points_on_this_proc)
                call fde_mpi_gatherv(active_gzn_bx, active_gridfunction%gn_b_reorth(3,1,:), nr_points_on_this_proc)
                call fde_mpi_gatherv(active_gzn_by, active_gridfunction%gn_b_reorth(3,2,:), nr_points_on_this_proc)
                call fde_mpi_gatherv(active_gzn_bz, active_gridfunction%gn_b_reorth(3,3,:), nr_points_on_this_proc)
             end if
           end if
           if (do_london_rhs_reorth_coupling.and..not. fde_cfg_no_sdft) then
              call fde_mpi_gatherv(active_sx_bx, active_gridfunction%s_b_reorth(1,1,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sx_by, active_gridfunction%s_b_reorth(1,2,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sx_bz, active_gridfunction%s_b_reorth(1,3,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sy_bx, active_gridfunction%s_b_reorth(2,1,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sy_by, active_gridfunction%s_b_reorth(2,2,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sy_bz, active_gridfunction%s_b_reorth(2,3,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sz_bx, active_gridfunction%s_b_reorth(3,1,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sz_by, active_gridfunction%s_b_reorth(3,2,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sz_bz, active_gridfunction%s_b_reorth(3,3,:), nr_points_on_this_proc)
              if (use_gga) then
                 call fde_mpi_gatherv(active_gxsx_bx, active_gridfunction%gs_b_reorth(1,1,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsx_by, active_gridfunction%gs_b_reorth(1,1,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsx_bz, active_gridfunction%gs_b_reorth(1,1,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsy_bx, active_gridfunction%gs_b_reorth(1,2,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsy_by, active_gridfunction%gs_b_reorth(1,2,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsy_bz, active_gridfunction%gs_b_reorth(1,2,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsz_bx, active_gridfunction%gs_b_reorth(1,3,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsz_by, active_gridfunction%gs_b_reorth(1,3,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsz_bz, active_gridfunction%gs_b_reorth(1,3,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysx_bx, active_gridfunction%gs_b_reorth(2,1,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysx_by, active_gridfunction%gs_b_reorth(2,1,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysx_bz, active_gridfunction%gs_b_reorth(2,1,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysy_bx, active_gridfunction%gs_b_reorth(2,2,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysy_by, active_gridfunction%gs_b_reorth(2,2,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysy_bz, active_gridfunction%gs_b_reorth(2,2,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysz_bx, active_gridfunction%gs_b_reorth(2,3,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysz_by, active_gridfunction%gs_b_reorth(2,3,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysz_bz, active_gridfunction%gs_b_reorth(2,3,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsx_bx, active_gridfunction%gs_b_reorth(3,1,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsx_by, active_gridfunction%gs_b_reorth(3,1,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsx_bz, active_gridfunction%gs_b_reorth(3,1,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsy_bx, active_gridfunction%gs_b_reorth(3,2,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsy_by, active_gridfunction%gs_b_reorth(3,2,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsy_bz, active_gridfunction%gs_b_reorth(3,2,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsz_bx, active_gridfunction%gs_b_reorth(3,3,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsz_by, active_gridfunction%gs_b_reorth(3,3,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsz_bz, active_gridfunction%gs_b_reorth(3,3,3,:), nr_points_on_this_proc)
              end if
           end if
           if (do_london_rhs_direct_coupling.or.do_london_susc2el_coupling.and..not. fde_cfg_no_sdft) then
              call fde_mpi_gatherv(active_sx_bx, active_gridfunction%s_b_direct(1,1,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sx_by, active_gridfunction%s_b_direct(1,2,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sx_bz, active_gridfunction%s_b_direct(1,3,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sy_bx, active_gridfunction%s_b_direct(2,1,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sy_by, active_gridfunction%s_b_direct(2,2,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sy_bz, active_gridfunction%s_b_direct(2,3,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sz_bx, active_gridfunction%s_b_direct(3,1,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sz_by, active_gridfunction%s_b_direct(3,2,:), nr_points_on_this_proc)
              call fde_mpi_gatherv(active_sz_bz, active_gridfunction%s_b_direct(3,3,:), nr_points_on_this_proc)
              if (use_gga) then
                 call fde_mpi_gatherv(active_gxsx_bx, active_gridfunction%gs_b_direct(1,1,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsx_by, active_gridfunction%gs_b_direct(1,1,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsx_bz, active_gridfunction%gs_b_direct(1,1,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsy_bx, active_gridfunction%gs_b_direct(1,2,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsy_by, active_gridfunction%gs_b_direct(1,2,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsy_bz, active_gridfunction%gs_b_direct(1,2,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsz_bx, active_gridfunction%gs_b_direct(1,3,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsz_by, active_gridfunction%gs_b_direct(1,3,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gxsz_bz, active_gridfunction%gs_b_direct(1,3,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysx_bx, active_gridfunction%gs_b_direct(2,1,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysx_by, active_gridfunction%gs_b_direct(2,1,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysx_bz, active_gridfunction%gs_b_direct(2,1,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysy_bx, active_gridfunction%gs_b_direct(2,2,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysy_by, active_gridfunction%gs_b_direct(2,2,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysy_bz, active_gridfunction%gs_b_direct(2,2,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysz_bx, active_gridfunction%gs_b_direct(2,3,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysz_by, active_gridfunction%gs_b_direct(2,3,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gysz_bz, active_gridfunction%gs_b_direct(2,3,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsx_bx, active_gridfunction%gs_b_direct(3,1,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsx_by, active_gridfunction%gs_b_direct(3,1,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsx_bz, active_gridfunction%gs_b_direct(3,1,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsy_bx, active_gridfunction%gs_b_direct(3,2,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsy_by, active_gridfunction%gs_b_direct(3,2,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsy_bz, active_gridfunction%gs_b_direct(3,2,3,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsz_bx, active_gridfunction%gs_b_direct(3,3,1,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsz_by, active_gridfunction%gs_b_direct(3,3,2,:), nr_points_on_this_proc)
                 call fde_mpi_gatherv(active_gzsz_bz, active_gridfunction%gs_b_direct(3,3,3,:), nr_points_on_this_proc)
              end if
           end if
        end if
#endif
        if (fde_mpi_is_master()) then
           if (do_london_rhs_reorth_coupling) then
             active_gridfunction%n_b_reorth(1,:) = active_n_bx(:)
             active_gridfunction%n_b_reorth(2,:) = active_n_by(:)
             active_gridfunction%n_b_reorth(3,:) = active_n_bz(:)
             if (use_gga) then
                active_gridfunction%gn_b_reorth(1,1,:) = active_gxn_bx(:)
                active_gridfunction%gn_b_reorth(1,2,:) = active_gxn_by(:)
                active_gridfunction%gn_b_reorth(1,3,:) = active_gxn_bz(:)
                active_gridfunction%gn_b_reorth(2,1,:) = active_gyn_bx(:)
                active_gridfunction%gn_b_reorth(2,2,:) = active_gyn_by(:)
                active_gridfunction%gn_b_reorth(2,3,:) = active_gyn_bz(:)
                active_gridfunction%gn_b_reorth(3,1,:) = active_gzn_bx(:)
                active_gridfunction%gn_b_reorth(3,2,:) = active_gzn_by(:)
                active_gridfunction%gn_b_reorth(3,3,:) = active_gzn_bz(:)
             endif
           endif
           if (do_london_rhs_reorth_coupling.and..not. fde_cfg_no_sdft) then
              active_gridfunction%s_b_reorth(1,1,:) = active_sx_bx(:)
              active_gridfunction%s_b_reorth(1,2,:) = active_sx_by(:)
              active_gridfunction%s_b_reorth(1,3,:) = active_sx_bz(:)
              active_gridfunction%s_b_reorth(2,1,:) = active_sy_bx(:)
              active_gridfunction%s_b_reorth(2,2,:) = active_sy_by(:)
              active_gridfunction%s_b_reorth(2,3,:) = active_sy_bz(:)
              active_gridfunction%s_b_reorth(3,1,:) = active_sz_bx(:)
              active_gridfunction%s_b_reorth(3,2,:) = active_sz_by(:)
              active_gridfunction%s_b_reorth(3,3,:) = active_sz_bz(:)
              if (use_gga) then
                 active_gridfunction%gs_b_reorth(1,1,1,:) = active_gxsx_bx(:)
                 active_gridfunction%gs_b_reorth(1,1,2,:) = active_gxsx_by(:)
                 active_gridfunction%gs_b_reorth(1,1,3,:) = active_gxsx_bz(:)
                 active_gridfunction%gs_b_reorth(1,2,1,:) = active_gxsy_bx(:)
                 active_gridfunction%gs_b_reorth(1,2,2,:) = active_gxsy_by(:)
                 active_gridfunction%gs_b_reorth(1,2,3,:) = active_gxsy_bz(:)
                 active_gridfunction%gs_b_reorth(1,3,1,:) = active_gxsz_bx(:)
                 active_gridfunction%gs_b_reorth(1,3,2,:) = active_gxsz_by(:)
                 active_gridfunction%gs_b_reorth(1,3,3,:) = active_gxsz_bz(:)
                 active_gridfunction%gs_b_reorth(2,1,1,:) = active_gysx_bx(:)
                 active_gridfunction%gs_b_reorth(2,1,2,:) = active_gysx_by(:)
                 active_gridfunction%gs_b_reorth(2,1,3,:) = active_gysx_bz(:)
                 active_gridfunction%gs_b_reorth(2,2,1,:) = active_gysy_bx(:)
                 active_gridfunction%gs_b_reorth(2,2,2,:) = active_gysy_by(:)
                 active_gridfunction%gs_b_reorth(2,2,3,:) = active_gysy_bz(:)
                 active_gridfunction%gs_b_reorth(2,3,1,:) = active_gysz_bx(:)
                 active_gridfunction%gs_b_reorth(2,3,2,:) = active_gysz_by(:)
                 active_gridfunction%gs_b_reorth(2,3,3,:) = active_gysz_bz(:)
                 active_gridfunction%gs_b_reorth(3,1,1,:) = active_gzsx_bx(:)
                 active_gridfunction%gs_b_reorth(3,1,2,:) = active_gzsx_by(:)
                 active_gridfunction%gs_b_reorth(3,1,3,:) = active_gzsx_bz(:)
                 active_gridfunction%gs_b_reorth(3,2,1,:) = active_gzsy_bx(:)
                 active_gridfunction%gs_b_reorth(3,2,2,:) = active_gzsy_by(:)
                 active_gridfunction%gs_b_reorth(3,2,3,:) = active_gzsy_bz(:)
                 active_gridfunction%gs_b_reorth(3,3,1,:) = active_gzsz_bx(:)
                 active_gridfunction%gs_b_reorth(3,3,2,:) = active_gzsz_by(:)
                 active_gridfunction%gs_b_reorth(3,3,3,:) = active_gzsz_bz(:)
              end if
           end if
           if (do_london_rhs_direct_coupling.or.do_london_susc2el_coupling.and..not. fde_cfg_no_sdft) then
              active_gridfunction%s_b_direct(1,1,:) = active_sx_bx(:)
              active_gridfunction%s_b_direct(1,2,:) = active_sx_by(:)
              active_gridfunction%s_b_direct(1,3,:) = active_sx_bz(:)
              active_gridfunction%s_b_direct(2,1,:) = active_sy_bx(:)
              active_gridfunction%s_b_direct(2,2,:) = active_sy_by(:)
              active_gridfunction%s_b_direct(2,3,:) = active_sy_bz(:)
              active_gridfunction%s_b_direct(3,1,:) = active_sz_bx(:)
              active_gridfunction%s_b_direct(3,2,:) = active_sz_by(:)
              active_gridfunction%s_b_direct(3,3,:) = active_sz_bz(:)
              if (use_gga) then
                 active_gridfunction%gs_b_direct(1,1,1,:) = active_gxsx_bx(:)
                 active_gridfunction%gs_b_direct(1,1,2,:) = active_gxsx_by(:)
                 active_gridfunction%gs_b_direct(1,1,3,:) = active_gxsx_bz(:)
                 active_gridfunction%gs_b_direct(1,2,1,:) = active_gxsy_bx(:)
                 active_gridfunction%gs_b_direct(1,2,2,:) = active_gxsy_by(:)
                 active_gridfunction%gs_b_direct(1,2,3,:) = active_gxsy_bz(:)
                 active_gridfunction%gs_b_direct(1,3,1,:) = active_gxsz_bx(:)
                 active_gridfunction%gs_b_direct(1,3,2,:) = active_gxsz_by(:)
                 active_gridfunction%gs_b_direct(1,3,3,:) = active_gxsz_bz(:)
                 active_gridfunction%gs_b_direct(2,1,1,:) = active_gysx_bx(:)
                 active_gridfunction%gs_b_direct(2,1,2,:) = active_gysx_by(:)
                 active_gridfunction%gs_b_direct(2,1,3,:) = active_gysx_bz(:)
                 active_gridfunction%gs_b_direct(2,2,1,:) = active_gysy_bx(:)
                 active_gridfunction%gs_b_direct(2,2,2,:) = active_gysy_by(:)
                 active_gridfunction%gs_b_direct(2,2,3,:) = active_gysy_bz(:)
                 active_gridfunction%gs_b_direct(2,3,1,:) = active_gysz_bx(:)
                 active_gridfunction%gs_b_direct(2,3,2,:) = active_gysz_by(:)
                 active_gridfunction%gs_b_direct(2,3,3,:) = active_gysz_bz(:)
                 active_gridfunction%gs_b_direct(3,1,1,:) = active_gzsx_bx(:)
                 active_gridfunction%gs_b_direct(3,1,2,:) = active_gzsx_by(:)
                 active_gridfunction%gs_b_direct(3,1,3,:) = active_gzsx_bz(:)
                 active_gridfunction%gs_b_direct(3,2,1,:) = active_gzsy_bx(:)
                 active_gridfunction%gs_b_direct(3,2,2,:) = active_gzsy_by(:)
                 active_gridfunction%gs_b_direct(3,2,3,:) = active_gzsy_bz(:)
                 active_gridfunction%gs_b_direct(3,3,1,:) = active_gzsz_bx(:)
                 active_gridfunction%gs_b_direct(3,3,2,:) = active_gzsz_by(:)
                 active_gridfunction%gs_b_direct(3,3,3,:) = active_gzsz_bz(:)
              end if
           end if
        endif

        deallocate(active_n_bx)
        deallocate(active_n_by)
        deallocate(active_n_bz)
        if (use_gga) then
           deallocate(active_gxn_bx)
           deallocate(active_gxn_by)
           deallocate(active_gxn_bz)
           deallocate(active_gyn_bx)
           deallocate(active_gyn_by)
           deallocate(active_gyn_bz)
           deallocate(active_gzn_bx)
           deallocate(active_gzn_by)
           deallocate(active_gzn_bz)
        end if
        if (.not. fde_cfg_no_sdft) then
           deallocate(active_sx_bx)
           deallocate(active_sx_by)
           deallocate(active_sx_bz)
           deallocate(active_sy_bx)
           deallocate(active_sy_by)
           deallocate(active_sy_bz)
           deallocate(active_sz_bx)
           deallocate(active_sz_by)
           deallocate(active_sz_bz)
           if (use_gga) then
              deallocate(active_gxsx_bx)
              deallocate(active_gxsx_by)
              deallocate(active_gxsx_bz)
              deallocate(active_gxsy_bx)
              deallocate(active_gxsy_by)
              deallocate(active_gxsy_bz)
              deallocate(active_gxsz_bx)
              deallocate(active_gxsz_by)
              deallocate(active_gxsz_bz)
              deallocate(active_gysx_bx)
              deallocate(active_gysx_by)
              deallocate(active_gysx_bz)
              deallocate(active_gysy_bx)
              deallocate(active_gysy_by)
              deallocate(active_gysy_bz)
              deallocate(active_gysz_bx)
              deallocate(active_gysz_by)
              deallocate(active_gysz_bz)
              deallocate(active_gzsx_bx)
              deallocate(active_gzsx_by)
              deallocate(active_gzsx_bz)
              deallocate(active_gzsy_bx)
              deallocate(active_gzsy_by)
              deallocate(active_gzsy_bz)
              deallocate(active_gzsz_bx)
              deallocate(active_gzsz_by)
              deallocate(active_gzsz_bz)
           end if
        end if

! don't try to read another batch
        goto 999

!     end loop over batches
     end do
 999 continue

     end subroutine loop_over_batches_pertden

! ------------------------------------------------------------------------------
   subroutine get_electrostatic_potential_and_density()
! ------------------------------------------------------------------------------
!        unperturbed density
   real(8) :: n_0(max_block_length)
!        unperturbed density gradient
   real(8) :: gn_0_pointwise(3)
   real(8) :: hn_0_pointwise(6)
!     --------------------------------------------------------------------------
   real(8) :: temp
   real(8) :: block_threshold
   integer :: ipoint
   integer :: k
   integer :: i
   logical :: print_data_at_point
   integer :: np_esp
   integer :: irrep_esp

   real(8), allocatable :: ao(:, :)
   real(8), allocatable :: buffer(:, :)
   real(8), allocatable :: dmat_local(:, :, :)
   real(8), allocatable :: r(:, :)

   if (.not. fde_use_blas3_scheme) then

      allocate(ao(1,     nr_ao_slices*nr_ao_cartesian))
      ao = 0.0d0
      allocate(buffer(1, nr_ao_slices*nr_ao_cartesian))
      buffer = 0.0d0

      allocate(r(3,nr_points_on_this_proc))
      r = 0.0d0

!     allocate(dmat_local(mat_dim,mat_dim,nz))
!     dmat_local = 0.0d0
!     dmat_local(:,:,1) = dmat_0(:,:)

      print_data_at_point = .true.

      np_esp    = 1
      irrep_esp = 1

      active_elpot = 0.0d0
      active_nucpot= 0.0d0
      active_n     = 0.0d0
      active_gnx   = 0.0d0
      active_gny   = 0.0d0
      active_gnz   = 0.0d0
      active_hnxx  = 0.0d0
      active_hnxy  = 0.0d0
      active_hnxz  = 0.0d0
      active_hnyy  = 0.0d0
      active_hnyz  = 0.0d0
      active_hnzz  = 0.0d0

!  here starts the expensive loop over points
!  do not evaluate anything inside this loop
!  that does not change from point to point
      do ipoint = 1, nr_points_on_this_proc
      
         time_ao_start = second()
         call get_ao(1,    &
               rx(ipoint), &
               ry(ipoint), &
               rz(ipoint), &
               ao,         &
               buffer)
!        time_ao = time_ao + second() - time_ao_start

!        time_density_start = second()

         gn_0_pointwise(1) = 0.0d0 
         gn_0_pointwise(2) = 0.0d0 
         gn_0_pointwise(3) = 0.0d0

         hn_0_pointwise(1) = 0.0d0 
         hn_0_pointwise(2) = 0.0d0 
         hn_0_pointwise(3) = 0.0d0
         hn_0_pointwise(4) = 0.0d0 
         hn_0_pointwise(5) = 0.0d0 
         hn_0_pointwise(6) = 0.0d0

         if (need_ao_order(1, 0)) then
!           density and density gradient
            call get_gn(n_0(1), gn_0_pointwise, 0, mat_dim, dmat_0, buffer, ao)
         else
!           density
            call get_n(n_0(1), 0, mat_dim, dmat_0, buffer, ao)
         endif

         active_n(ipoint)   = n_0(1) 

         active_gnx(ipoint) = gn_0_pointwise(1)
         active_gny(ipoint) = gn_0_pointwise(2)
         active_gnz(ipoint) = gn_0_pointwise(3)

         active_hnxx(ipoint)= hn_0_pointwise(1)
         active_hnxy(ipoint)= hn_0_pointwise(2)
         active_hnxz(ipoint)= hn_0_pointwise(3)
         active_hnyy(ipoint)= hn_0_pointwise(4)
         active_hnyz(ipoint)= hn_0_pointwise(5)
         active_hnzz(ipoint)= hn_0_pointwise(6)

         r(1,ipoint) = rx(ipoint)
         r(2,ipoint) = ry(ipoint)
         r(3,ipoint) = rz(ipoint)

    end do !do ipoint = 1, nr_points_on_this_proc
    end if

!   time_density = time_density + second() - time_density_start
    call get_esp(nr_points_on_this_proc,&
              active_elpot,             &
              irrep_esp,                &
              mat_dim,                  &
!             dmat_local,               &
              dmat_0,                   &
              r,                        &
              include_nuc_part=.true.,  &
              include_el_part=.true.)

! aspg, 14/10/2015
!       note that the electrostatic, hartree and nuclear potentials in dirac have
!       been multiplied by -1 when exported so that we don't have to multipy the 
!       product potential * density by the electron charge when calculating the
!       energy. with this, the potentials will be in line with those obtained
!       e.g. with adf 
!
    active_elpot = -active_elpot

    call get_esp(nr_points_on_this_proc,&
              active_nucpot,            &
              irrep_esp,                &
              mat_dim,                  &
!             dmat_local,               &
              dmat_0,                   &
              r,                        &
              include_nuc_part=.true.,  &
              include_el_part=.false.)

    active_nucpot = -active_nucpot
!

    deallocate(ao)
    deallocate(buffer)
!   deallocate(dmat_local)
    deallocate(r)

   end subroutine 

! ------------------------------------------------------------------------------
   subroutine get_perturbed_density()
! ------------------------------------------------------------------------------
!        unperturbed density
   real(8) :: n_0(max_block_length)
!        unperturbed density gradient
   real(8) :: gn_0_pointwise(3)
   real(8) :: gnn_0(max_block_length)
   real(8) :: n_0_t(max_block_length)
   real(8) :: gnn_0_t(max_block_length)
!     --------------------------------------------------------------------------
   real(8) :: w(max_block_length)
   real(8) :: temp
   real(8) :: block_threshold
   integer :: ipoint
   integer :: k
   integer :: i
   logical :: print_data_at_point
   logical :: nuc_part
   logical :: ele_part
   integer :: np_esp
   integer :: irrep_esp

   real(8), allocatable :: ao(:, :)
   real(8), allocatable :: buffer(:, :)

   real(8) :: z0tt(3)

   if (.not. fde_use_blas3_scheme) then

      allocate(ao(1,     nr_ao_slices*nr_ao_cartesian))
      ao = 0.0d0
      allocate(buffer(1, nr_ao_slices*nr_ao_cartesian))
      buffer = 0.0d0


      active_n_bx = 0.0d0
      active_n_by = 0.0d0
      active_n_bz = 0.0d0
      if (use_gga) then
         active_gxn_bx = 0.0d0
         active_gxn_by = 0.0d0
         active_gxn_bz = 0.0d0
         active_gyn_bx = 0.0d0
         active_gyn_by = 0.0d0
         active_gyn_bz = 0.0d0
         active_gzn_bx = 0.0d0
         active_gzn_by = 0.0d0
         active_gzn_bz = 0.0d0
      end if
      if (.not. fde_cfg_no_sdft) then
         active_sx_bx = 0.0d0
         active_sx_by = 0.0d0
         active_sx_bz = 0.0d0
         active_sy_bx = 0.0d0
         active_sy_by = 0.0d0
         active_sy_bz = 0.0d0
         active_sz_bx = 0.0d0
         active_sz_by = 0.0d0
         active_sz_bz = 0.0d0
         if (use_gga) then
            active_gxsx_bx = 0.0d0
            active_gxsx_by = 0.0d0
            active_gxsx_bz = 0.0d0
            active_gxsy_bx = 0.0d0
            active_gxsy_by = 0.0d0
            active_gxsy_bz = 0.0d0
            active_gxsz_bx = 0.0d0
            active_gxsz_by = 0.0d0
            active_gxsz_bz = 0.0d0
            active_gysx_bx = 0.0d0
            active_gysx_by = 0.0d0
            active_gysx_bz = 0.0d0
            active_gysy_bx = 0.0d0
            active_gysy_by = 0.0d0
            active_gysy_bz = 0.0d0
            active_gysz_bx = 0.0d0
            active_gysz_by = 0.0d0
            active_gysz_bz = 0.0d0
            active_gzsx_bx = 0.0d0
            active_gzsx_by = 0.0d0
            active_gzsx_bz = 0.0d0
            active_gzsy_bx = 0.0d0
            active_gzsy_by = 0.0d0
            active_gzsy_bz = 0.0d0
            active_gzsz_bx = 0.0d0
            active_gzsz_by = 0.0d0
            active_gzsz_bz = 0.0d0
         end if
      end if

!  here starts the expensive loop over points
!  do not evaluate anything inside this loop
!  that does not change from point to point
         !do i = 1, 3
         !  write(*, *) 'dmat_pertden_reorth read test, i =', i
         !  call prqmat(dmat_pertden_reorth(1,1,1,i), mat_dim, mat_dim, mat_dim, mat_dim, nz,  &
         !              (/1,2,3,4/), 6)
         !end do

      do ipoint = 1, nr_points_on_this_proc
   
         time_ao_start = second()
         
         call get_ao(1,    &
               rx(ipoint), &
               ry(ipoint), &
               rz(ipoint), &
               ao,         &
               buffer)

         !write(*, *) 'ao dft nr_ao_slices, nr_ao_cartesian = ', nr_ao_slices, nr_ao_cartesian
         !do i = 1, nr_ao_slices*nr_ao_cartesian
         !  write(*, '(A,3E16.8,ES30.20E3)') 'ao in point ', rx(ipoint), ry(ipoint), rz(ipoint), ao(1,i)
         !end do

         if (need_ao_order(1, 0)) then
!           density and density gradient
            call get_gn(n_0(1), gn_0_pointwise, 0, mat_dim, dmat_0, buffer, ao)
            gnn_0(1) = gn_0_pointwise(1)*gn_0_pointwise(1) &
                     + gn_0_pointwise(2)*gn_0_pointwise(2) &
                     + gn_0_pointwise(3)*gn_0_pointwise(3)
         else
!           density
            call get_n(n_0(1), 0, mat_dim, dmat_0, buffer, ao)
         endif 


         if (n_0(1) > fde_cfg_tinydens) then
            nr_points_used = nr_points_used + 1
            w(1) = rw(ipoint)

            n_0_t(1)   = n_0(1)
            gnn_0_t(1) = gnn_0(1)
         
!           lao_reorth
            if (do_export_pertden_reorth) then

               call get_pertden_reorth(ipoint,  &
                                       mat_dim, &
                                       nz,      &
                                       ao,      &
                                       n_0,     &
                                       n_0_t,   &
                                       gn_0_pointwise, &
                                       gnn_0,   &
                                       gnn_0_t, &
                                       w,       &
                                       use_gga, &
                                       dmat_pertden_reorth,    &
                                       dmat_pg_sym,    &
                                       buffer,  &
                       rx(ipoint), ry(ipoint), rz(ipoint))
            end if
!           lao_direct
            if (do_export_pertden_direct .and. .not. fde_cfg_no_sdft) then
               call get_pertden_direct(ipoint,         &
                                      mat_dim,        &
                                      nz,             &
                                      ao,             &
                                      buffer,         &
                                      n_0,            &
                                      n_0_t,          &
                                      gn_0_pointwise, &
                                      gnn_0,          &
                                      gnn_0_t,        &
                                      w,              &
                                      use_gga,        &
                                      dmat_pertden_direct,           &
                                      dmat_pg_sym,    &
                                      rx(ipoint), ry(ipoint), rz(ipoint))
            end if
         
         end if !if (n_0(1) > fde_cfg_tinydens) then
   end do !do ipoint = 1, nr_points_on_this_proc
   end if

    deallocate(ao)
    deallocate(buffer)


   end subroutine

   subroutine get_pertden_direct(ipoint,         &
                                 mat_dim,        &
                                 nz,             &
                                 ao,             &
                                 buffer,         &
                                 n_0,            &
                                 n_0_t,          &
                                 gn_0_pointwise, &
                                 gnn_0,          &
                                 gnn_0_t,        &
                                 w,              &
                                 is_gga,         &
                                 dmat_full,      &
                                 isym,           &
                                 px, py, pz)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ipoint
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: buffer(*)
    real(8), intent(in)    :: n_0(*)
    real(8), intent(in)    :: n_0_t(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: gnn_0(*)
    real(8), intent(in)    :: gnn_0_t(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga
    real(8), intent(in)    :: dmat_full(*)
    integer, intent(in)    :: isym(*)
    real(8), intent(in)    :: px, py, pz
!   ----------------------------------------------------------------------------
    real(8)                ::  s_b(3, 3)
    real(8)                ::  gs_b(3, 3, 3)

    real(8)                :: y_b(3, 3)
    real(8)                :: y_b_t(3, 3)

    integer                :: k_start, k, a, im_off, j
    type(fde_omega_prefactor) :: u
    integer :: i
!   ----------------------------------------------------------------------------

    s_b    = 0.0d0
    gs_b   = 0.0d0

    k_start = 1
    if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

    if (is_gga) then
      call get_gs_lao(s_b,         &
                      gs_b,        &
                      (/1,1,1/),   &
                      mat_dim,     &
                      dmat_full,   &
                      ao,          &
                      px, py, pz, buffer)

    else
       call get_s_lao(s_b,              &
                      (/0, 0, 0/),      &
                      mat_dim,          &
                      dmat_full,        &
                      ao,               &
                      px, py, pz, buffer)
    end if

    !write(*, '(A,3F16.8,9ES30.20E3)') 'test23 get_pertden_direct= ', &
    !            px, py, pz, s_b(1,1), s_b(2,1), s_b(3,1), s_b(1,2), s_b(2,2), s_b(3,2), s_b(1,3), s_b(2,3), s_b(3,3)

! gosia: s_b and gs_b contributions can also come from reorthonormalization part
! which is calculated earlier
! that is why here we sum all up to *active* variables
    active_sx_bx(ipoint) = active_sx_bx(ipoint) + s_b(1,1)
    active_sx_by(ipoint) = active_sx_by(ipoint) + s_b(1,2)
    active_sx_bz(ipoint) = active_sx_bz(ipoint) + s_b(1,3)
    active_sy_bx(ipoint) = active_sy_bx(ipoint) + s_b(2,1)
    active_sy_by(ipoint) = active_sy_by(ipoint) + s_b(2,2)
    active_sy_bz(ipoint) = active_sy_bz(ipoint) + s_b(2,3)
    active_sz_bx(ipoint) = active_sz_bx(ipoint) + s_b(3,1)
    active_sz_by(ipoint) = active_sz_by(ipoint) + s_b(3,2)
    active_sz_bz(ipoint) = active_sz_bz(ipoint) + s_b(3,3)
    if (is_gga) then
       active_gxsx_bx(ipoint) = active_gxsx_bx(ipoint) + gs_b(1,1,1)
       active_gxsx_by(ipoint) = active_gxsx_by(ipoint) + gs_b(1,1,2)
       active_gxsx_bz(ipoint) = active_gxsx_bz(ipoint) + gs_b(1,1,3)
       active_gxsy_bx(ipoint) = active_gxsy_bx(ipoint) + gs_b(1,2,1)
       active_gxsy_by(ipoint) = active_gxsy_by(ipoint) + gs_b(1,2,2)
       active_gxsy_bz(ipoint) = active_gxsy_bz(ipoint) + gs_b(1,2,3)
       active_gxsz_bx(ipoint) = active_gxsz_bx(ipoint) + gs_b(1,3,1)
       active_gxsz_by(ipoint) = active_gxsz_by(ipoint) + gs_b(1,3,2)
       active_gxsz_bz(ipoint) = active_gxsz_bz(ipoint) + gs_b(1,3,3)
       active_gysx_bx(ipoint) = active_gysx_bx(ipoint) + gs_b(2,1,1)
       active_gysx_by(ipoint) = active_gysx_by(ipoint) + gs_b(2,1,2)
       active_gysx_bz(ipoint) = active_gysx_bz(ipoint) + gs_b(2,1,3)
       active_gysy_bx(ipoint) = active_gysy_bx(ipoint) + gs_b(2,2,1)
       active_gysy_by(ipoint) = active_gysy_by(ipoint) + gs_b(2,2,2)
       active_gysy_bz(ipoint) = active_gysy_bz(ipoint) + gs_b(2,2,3)
       active_gysz_bx(ipoint) = active_gysz_bx(ipoint) + gs_b(2,3,1)
       active_gysz_by(ipoint) = active_gysz_by(ipoint) + gs_b(2,3,2)
       active_gysz_bz(ipoint) = active_gysz_bz(ipoint) + gs_b(2,3,3)
       active_gzsx_bx(ipoint) = active_gzsx_bx(ipoint) + gs_b(3,1,1)
       active_gzsx_by(ipoint) = active_gzsx_by(ipoint) + gs_b(3,1,2)
       active_gzsx_bz(ipoint) = active_gzsx_bz(ipoint) + gs_b(3,1,3)
       active_gzsy_bx(ipoint) = active_gzsy_bx(ipoint) + gs_b(3,2,1)
       active_gzsy_by(ipoint) = active_gzsy_by(ipoint) + gs_b(3,2,2)
       active_gzsy_bz(ipoint) = active_gzsy_bz(ipoint) + gs_b(3,2,3)
       active_gzsz_bx(ipoint) = active_gzsz_bx(ipoint) + gs_b(3,3,1)
       active_gzsz_by(ipoint) = active_gzsz_by(ipoint) + gs_b(3,3,2)
       active_gzsz_bz(ipoint) = active_gzsz_bz(ipoint) + gs_b(3,3,3)
    end if
 

   end subroutine

  subroutine get_pertden_reorth(ipoint,     &
                        mat_dim,            &
                        nz,                 &
                        ao,                 &
                        n_0,                &
                        n_0_t,              &
                        gn_0_pointwise,     &
                        gnn_0,              &
                        gnn_0_t,            &
                        w,                  &
                        is_gga,             &
                        tb_dmat,            &
                        isym,               &
                        buffer,             &
                        px, py, pz)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ipoint
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: n_0(*)
    real(8), intent(in)    :: n_0_t(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: gnn_0(*)
    real(8), intent(in)    :: gnn_0_t(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga
    real(8), intent(in)    :: tb_dmat(*)
    integer, intent(in)    :: isym(*)
    real(8), intent(in)    :: buffer(*)
    real(8), intent(in)    :: px, py, pz
!   ----------------------------------------------------------------------------
    real(8)                ::  r_b(3),     s_b(3,3)
    real(8)                :: gr_b(3,3), gs_b(3,3,3)
    real(8)                ::  t_b,    ts_b(3)

    real(8)                :: z_b
    real(8)                :: z_b_t
    real(8)                :: y_b(3)
    real(8)                :: y_b_t(3)

    integer                :: k_start, k, im, im_off, j
    integer :: i
!   ----------------------------------------------------------------------------


    r_b = 0.0d0
    gr_b = 0.0d0
    s_b = 0.0d0
    gs_b = 0.0d0


    k_start = 1
    if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

    do im = 1, 3

       im_off = mat_dim*mat_dim*nz*(im - 1)

!      charge density
       if (is_gga) then
          call get_gn_nonhermitian(r_b(im),                 &
                                   gr_b(1:3,im),            &
                                   !isym(im) - 1,            &
                                   0,                       &
                                   mat_dim,                 &
                                   tb_dmat(1 + im_off),     &
                                   buffer,                  &
                                   ao)
       else
!     write(*, *) 'call get_n in get_pertden_reorth, im, im_off = ', im, im_off
          call get_n(r_b(im),                       &
                     0,              &
                     mat_dim,                   &
                     tb_dmat(1 + im_off),       &
                     buffer,                    &
                     ao)
!    write(*, '(A,3F16.8,I,I,ES30.20E3)') 'test22 get_pertden_reorth= ', &
!                px, py, pz, im, im_off, r_b(im)

       end if

!      spin density
       if (.not. fde_cfg_no_sdft) then
          if (is_gga) then
             call get_gs_nonhermitian(s_b(1:3,im),                 &
                                      gs_b(1:3,1:3,im),              &
                                      !isym(im) - 1,              &
                                      0,                         &
                                      mat_dim,                   &
                                      tb_dmat(1 + im_off),       &
                                      buffer,                    &
                                      ao)
          else
             call get_s(s_b(1:3,im),                  &
                        0,               &
                        mat_dim,                    &
                        tb_dmat(1 + im_off),        &
                        buffer,                     &
                        ao)
          end if
       end if
    end do

!    write(*, '(A,3F16.8,3ES30.20E3)') 'test20  get_pertden_reorth= ', &
!                px, py, pz, r_b(1), r_b(2), r_b(3)

!    write(*, '(A,3F16.8,9ES30.20E3)') 'test21 get_pertden_reorth= ', &
!                px, py, pz, s_b(1,1), s_b(2,1), s_b(3,1), &
!                            s_b(1,2), s_b(2,2), s_b(3,2), &
!                            s_b(1,3), s_b(2,3), s_b(3,3)

    active_n_bx(ipoint) = r_b(1)
    active_n_by(ipoint) = r_b(2)
    active_n_bz(ipoint) = r_b(3)

    if (use_gga) then
       active_gxn_bx(ipoint) = gr_b(1,1)
       active_gxn_by(ipoint) = gr_b(1,2)
       active_gxn_bz(ipoint) = gr_b(1,3)
       active_gyn_bx(ipoint) = gr_b(2,1)
       active_gyn_by(ipoint) = gr_b(2,2)
       active_gyn_bz(ipoint) = gr_b(2,3)
       active_gzn_bx(ipoint) = gr_b(3,1)
       active_gzn_by(ipoint) = gr_b(3,2)
       active_gzn_bz(ipoint) = gr_b(3,3)
    end if
    if (.not. fde_cfg_no_sdft) then
       active_sx_bx(ipoint) = s_b(1,1)
       active_sx_by(ipoint) = s_b(1,2)
       active_sx_bz(ipoint) = s_b(1,3)
       active_sy_bx(ipoint) = s_b(2,1)
       active_sy_by(ipoint) = s_b(2,2)
       active_sy_bz(ipoint) = s_b(2,3)
       active_sz_bx(ipoint) = s_b(3,1)
       active_sz_by(ipoint) = s_b(3,2)
       active_sz_bz(ipoint) = s_b(3,3)
       if (use_gga) then
          active_gxsx_bx(ipoint) = gs_b(1,1,1)
          active_gxsx_by(ipoint) = gs_b(1,1,2)
          active_gxsx_bz(ipoint) = gs_b(1,1,3)
          active_gxsy_bx(ipoint) = gs_b(1,2,1)
          active_gxsy_by(ipoint) = gs_b(1,2,2)
          active_gxsy_bz(ipoint) = gs_b(1,2,3)
          active_gxsz_bx(ipoint) = gs_b(1,3,1)
          active_gxsz_by(ipoint) = gs_b(1,3,2)
          active_gxsz_bz(ipoint) = gs_b(1,3,3)
          active_gysx_bx(ipoint) = gs_b(2,1,1)
          active_gysx_by(ipoint) = gs_b(2,1,2)
          active_gysx_bz(ipoint) = gs_b(2,1,3)
          active_gysy_bx(ipoint) = gs_b(2,2,1)
          active_gysy_by(ipoint) = gs_b(2,2,2)
          active_gysy_bz(ipoint) = gs_b(2,2,3)
          active_gysz_bx(ipoint) = gs_b(2,3,1)
          active_gysz_by(ipoint) = gs_b(2,3,2)
          active_gysz_bz(ipoint) = gs_b(2,3,3)
          active_gzsx_bx(ipoint) = gs_b(3,1,1)
          active_gzsx_by(ipoint) = gs_b(3,1,2)
          active_gzsx_bz(ipoint) = gs_b(3,1,3)
          active_gzsy_bx(ipoint) = gs_b(3,2,1)
          active_gzsy_by(ipoint) = gs_b(3,2,2)
          active_gzsy_bz(ipoint) = gs_b(3,2,3)
          active_gzsz_bx(ipoint) = gs_b(3,3,1)
          active_gzsz_by(ipoint) = gs_b(3,3,2)
          active_gzsz_bz(ipoint) = gs_b(3,3,3)
       end if
    end if


  end subroutine



end module fde_dirac_matrices_integration
