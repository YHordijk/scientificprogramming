
! define task symbols for CALL DIRAC_PARCTL( task )
#include "dirac_partask.h"

module dirac_interface

   use memory_allocator
   use matrix_defop_old
   use dirac_cfg
   use num_grid_gen

   implicit none

   public di_get_gbds
   public di_get_geomderiv_molgrad_dft
   public di_get_geomderiv_fxd_dft
   public di_get_geomderiv_gxd_dft
   public di_get_magderiv_fxd_dft
   public di_get_magderiv_gxd_dft
   public do_dft
   public get_s
   public get_f
   public get_g
   public get_C
   public get_natoms
   public dipnuc_ifc
   public qdrnuc_ifc
   public gradnn_ifc
   public di_set_f77_work
   public di_select_wrk
   public di_deselect_wrk
   public daxpy_iz_jz
   public init_dirac_interface
   public lsquit
   public read_mo_coef
   public read_ibeig
   public dpgnuc_ifc

   save
   private

  integer                           :: unit_input       = 0
  integer, public                   :: unit_output      = 0

  integer, public                   :: nr_atoms         = 0
  integer                           :: nr_sym_dep_nuc   = 0
  integer                           :: nr_sym_indep_nuc = 0
  character(6), public, allocatable :: sym_dep_nuc_name(:)

  logical, public                   :: nonzero_small_metric = .false.

  integer, public                   :: integral_flag        = 0
  logical, public                   :: include_ll_integrals = .false.
  logical, public                   :: include_ls_integrals = .false.
  logical, public                   :: include_ss_integrals = .false.
  logical, public                   :: include_gt_integrals = .false.

  logical, public                   :: include_pp_rotations = .false.
  logical, public                   :: include_pn_rotations = .false.
  logical, public                   :: diamagnetic_via_pn   = .false.

  logical, public :: openrsp_levy_leblond = .false.
  logical, public :: openrsp_spinfree     = .false.

  logical :: dirac_interface_data_allocated = .false.
  logical :: dirac_interface_data_set       = .false.
#ifdef VAR_MPI
  logical :: dirac_interface_data_synced    = .false.
#endif

  interface di_get_gbds
    module procedure get_G
  end interface

#include "priunit.h"
#include "mxcent.h"
#include "../abacus/cbinuc.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbgrd.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dipole.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"

  real(8), pointer   :: f77_work(:)       ! Fortran 77 work array pointer
  integer, public    :: len_f77_work  = 0 ! Length of the work array
  integer            :: memory_usage  = 0 ! Amount of work array used
  integer            :: mem_wrk_count = 0 ! Counts

contains

  !> \brief gets the nuclear contribution to electric dipole gradient
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \param na is the number of atoms
  !> \return DGN contains the nuclear contributions to dipole gradient
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine DPGNUC_ifc( na, DGN )
    ! uses MXCOOR
#include "mxcent.h"
    ! uses DDIPN
#include "dipole.h"
    integer, intent(in) :: na
    real(8), intent(out) :: DGN( 3*na, 3 )
    real(8), parameter :: zero = 0.0D+00
    call DIPNUC(0, .true.)
    DGN = transpose( DDIPN( :, 1:3*na ) )
  end subroutine DPGNUC_ifc

   subroutine lsquit(text, u)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: text
      integer,      intent(in) :: u
!     --------------------------------------------------------------------------

!     dalton has lsquit, we use quit and ignore the unit
!     and use default lupri
      call quit(text)

   end subroutine

   subroutine di_get_geomDeriv_molgrad_DFT(g, nr_atoms, D)

      use xcint_main

!     --------------------------------------------------------------------------
      real(8),      intent(inout) :: g(*)
      integer,      intent(in)    :: nr_atoms
      type(matrix), intent(in)    :: D
!     --------------------------------------------------------------------------

      call generate_num_grid(D%elms)
#ifdef VAR_MPI
      if (parcal) call dirac_parctl( XCINT_PAR )
#endif
      call integrate_xc(xc_mat_dim           = ntbas(0), &
                         xc_nz                = nz,       &
                         xc_dmat_0            = D%elms,   &
                         xc_nr_dmat           = 0,        &
                         xc_nr_fmat           = 0,        &
                         xc_nr_atoms          = nr_atoms, &
                         xc_property_gradient = g,        &
                         xc_do_geo_0          = .true.)

   end subroutine

   subroutine di_get_geomderiv_fxd_dft(g, nr_atoms, D, T1)

      use xcint_main

!     --------------------------------------------------------------------------
      real(8),      intent(inout) :: g(*)
      integer,      intent(in)    :: nr_atoms
      type(matrix), intent(in)    :: D, T1
!     --------------------------------------------------------------------------

      if (iszero(T1)) then
         return
      end if

      call generate_num_grid(D%elms)
#ifdef VAR_MPI
      if (parcal) call dirac_parctl( XCINT_PAR )
#endif
      call integrate_xc(xc_mat_dim           = ntbas(0),    &
                         xc_nz                = nz,          &
                         xc_dmat_0            = D%elms,      &
                         xc_nr_dmat           = 1,           &
                         xc_dmat              = (/T1%elms/), &
                         xc_nr_fmat           = 0,           &
                         xc_nr_atoms          = nr_atoms,    &
                         xc_property_gradient = g,           &
                         xc_do_geo_1          = .true.)

   end subroutine

   subroutine di_get_geomderiv_gxd_dft(g, nr_atoms, D, T1, T2)

      use xcint_main

!     --------------------------------------------------------------------------
      real(8),      intent(inout) :: g(*)
      integer,      intent(in)    :: nr_atoms
      type(matrix), intent(in)    :: D, T1, T2
!     --------------------------------------------------------------------------

      if (iszero(T1) .or. iszero(T2)) then
         return
      end if

      call generate_num_grid(D%elms)
#ifdef VAR_MPI
      if (parcal) call dirac_parctl( XCINT_PAR )
#endif
      call integrate_xc(xc_mat_dim           = ntbas(0),             &
                         xc_nz                = nz,                   &
                         xc_dmat_0            = D%elms,               &
                         xc_nr_dmat           = 2,                    &
                         xc_dmat              = (/T1%elms, T2%elms/), &
                         xc_nr_fmat           = 0,                    &
                         xc_nr_atoms          = nr_atoms,             &
                         xc_property_gradient = g,                    &
                         xc_do_geo_2          = .true.)

   end subroutine

   subroutine di_get_magderiv_fxd_dft(F, D)

      use xcint_main

!     --------------------------------------------------------------------------
      type(matrix)              :: F(3)
      type(matrix), intent(in)  :: D
!     --------------------------------------------------------------------------
      real(8),      allocatable :: mat(:, :, :, :)
!     --------------------------------------------------------------------------

      call alloc(mat, ntbas(0), ntbas(0), nz, 3)
      mat = 0.0d0

      call generate_num_grid(D%elms)
#ifdef VAR_MPI
      if (parcal) call dirac_parctl( XCINT_PAR )
#endif
      call integrate_xc(xc_mat_dim      = ntbas(0),    &
                         xc_nz           = nz,          &
                         xc_dmat_0       = D%elms,      &
                         xc_nr_dmat      = 0,           &
                         xc_nr_fmat      = 3,           &
                         xc_fmat         = mat,         &
                         xc_fmat_pg_sym  = (/1, 1, 1/), &
                         xc_do_london_ks = .true.)

      call dcopy(ntbas(0)*ntbas(0)*nz, mat(1, 1, 1, 1), 1, F(1)%elms, 1)
      call dcopy(ntbas(0)*ntbas(0)*nz, mat(1, 1, 1, 2), 1, F(2)%elms, 1)
      call dcopy(ntbas(0)*ntbas(0)*nz, mat(1, 1, 1, 3), 1, F(3)%elms, 1)
      call dealloc(mat)

   end subroutine

   subroutine di_get_magderiv_gxd_dft(D, T, F)

      use xcint_main

!     --------------------------------------------------------------------------
      type(matrix), intent(in)  :: D
      type(matrix), intent(in)  :: T
      type(matrix)              :: F(3)
!     --------------------------------------------------------------------------
      real(8),      allocatable :: mat(:, :, :, :)
!     --------------------------------------------------------------------------

      call alloc(mat, ntbas(0), ntbas(0), nz, 3)
      mat = 0.0d0

      call generate_num_grid(D%elms)
#ifdef VAR_MPI
      if (parcal) call dirac_parctl( XCINT_PAR )
#endif
      call integrate_xc(xc_mat_dim      = ntbas(0),    &
                         xc_nz           = nz,          &
                         xc_dmat_0       = D%elms,      &
                         xc_nr_dmat      = 1,           &
                         xc_dmat         = T%elms,      &
                         xc_dmat_pg_sym  = (/1/),       &
                         xc_nr_fmat      = 3,           &
                         xc_fmat         = mat,         &
                         xc_fmat_pg_sym  = (/1, 1, 1/), &
                         xc_do_london_lr = .true.)

      call dcopy(ntbas(0)*ntbas(0)*nz, mat(1, 1, 1, 1), 1, F(1)%elms, 1)
      call dcopy(ntbas(0)*ntbas(0)*nz, mat(1, 1, 1, 2), 1, F(2)%elms, 1)
      call dcopy(ntbas(0)*ntbas(0)*nz, mat(1, 1, 1, 3), 1, F(3)%elms, 1)
      call dealloc(mat)

   end subroutine

  subroutine init_dirac_interface()


!   ----------------------------------------------------------------------------
    integer :: ill_int
    integer :: ils_int
    integer :: iss_int
    integer :: igt_int
!   ----------------------------------------------------------------------------

    if (.not. dirac_interface_data_allocated) then

      allocate(sym_dep_nuc_name(mxcent))

      dirac_interface_data_allocated = .true.

    end if

!   if (.not. dirac_interface_data_set) then

      unit_input       = lucmd
      unit_output      = lupri

      nr_atoms         = natoms
      nr_sym_dep_nuc   = nucdep
      nr_sym_indep_nuc = nucind
      sym_dep_nuc_name = namdep

      nonzero_small_metric = .false.
      if (mc == 2) then
        nonzero_small_metric = .true.
      end if
      if (bss .or. twocomp) then
        nonzero_small_metric = .false.
      end if

      ill_int = iand(intgen,   1)
      ils_int = iand(intgen/2, 1)
      iss_int = iand(intgen/4, 1)
      igt_int = iand(intgen/8, 1)

      integral_flag = 0
      if (ill_int == 1) then
        integral_flag = integral_flag + 1
        include_ll_integrals = .true.
      end if
      if (ils_int == 1) then
        integral_flag = integral_flag + 2
        include_ls_integrals = .true.
      end if
      if (iss_int == 1) then
        integral_flag = integral_flag + 4
        include_ss_integrals = .true.
      end if
      if (igt_int == 1) then
        integral_flag = integral_flag + 8
        include_gt_integrals = .true.
      end if


      openrsp_levy_leblond = levyle
      openrsp_spinfree     = spinfr

      dirac_interface_data_set = .true.

!   end if

  end subroutine

  subroutine daxpy_iz_jz(f, A, iz, B, jz)

!   ----------------------------------------------------------------------------
    real(8),      intent(in)    :: f
    integer,      intent(in)    :: iz, jz
    type(matrix), intent(inout) :: A, B
!   ----------------------------------------------------------------------------
    integer                     :: n
!   ----------------------------------------------------------------------------

    n = A%nrow*A%ncol
    call daxpy(n, f, A%elms((iz-1)*n + 1), 1, B%elms((jz-1)*n + 1), 1)

  end subroutine

  SUBROUTINE di_set_f77_work(wrk,lwork)
     real(8),target :: wrk(:)
     integer            :: lwork, ltempwork
     logical            :: ldummy1,ldummy2
     len_f77_work = lwork
     f77_work     => wrk
!    call get_linsca_work(ldummy1,ldummy2,ltempwork,lwork)
!    call di_set_coulomb_mem_req(ltempwork)
!    ! Simen: Refine!
!    call di_set_dft_mem_req(ltempwork)
  END SUBROUTINE di_set_f77_work

  subroutine di_select_wrk(work, need_mem)

!   ----------------------------------------------------------------------------
    real(8), pointer    :: work(:)
    integer, intent(in) :: need_mem
!   ----------------------------------------------------------------------------

    if (need_mem > len_f77_work) then

       call di_create_wrk(work, need_mem)

    else

      work => f77_work

      memory_usage  = need_mem
      mem_wrk_count = mem_wrk_count + 1

      if (mem_wrk_count > 1) call quit('counter > 1 not implemented')

    endif

  end subroutine

  subroutine di_deselect_wrk(work, need_mem)

!   ----------------------------------------------------------------------------
    real(8), pointer    :: work(:)
    integer, intent(in) :: need_mem
!   ----------------------------------------------------------------------------

    if (need_mem > len_f77_work) then

       call di_delete_wrk(work)

    else

      nullify(work)

      memory_usage  = 0
      mem_wrk_count = mem_wrk_count - 1

    endif

  end subroutine

  subroutine di_create_wrk(work, need_mem)

!   ----------------------------------------------------------------------------
    real(8), pointer    :: work(:)
    integer, intent(in) :: need_mem
!   ----------------------------------------------------------------------------
    integer             :: error
!   ----------------------------------------------------------------------------

    allocate(work(need_mem), stat=error)

    if (error /= 0) then
      call quit('unable to allocate work space in di_create_wrk')
    end if

  end subroutine

  subroutine di_delete_wrk(work)

!   ----------------------------------------------------------------------------
    real(8), pointer    :: work(:)
!   ----------------------------------------------------------------------------

    deallocate(work)

  end subroutine

  subroutine dipnuc_ifc(n, f, fg)

    integer       :: n
    real(8) :: f(3), fg(3, 3*n)

    call dipnuc((/0.0d0/), (/0.0d0/), 0, .true.)

    f  = dipmn
    fg = ddipn(:, 1:3*n)

  end subroutine

  subroutine get_natoms(return_nr_atoms)

!   ----------------------------------------------------------------------------
    integer, intent(out) :: return_nr_atoms
!   ----------------------------------------------------------------------------

    return_nr_atoms = nr_atoms

  end subroutine

  subroutine qdrnuc_ifc(q)

!   nuclear contribution to quadrupole moments
!   based on Kenneth's nucqdr in Dalton

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: q(6)
!   ----------------------------------------------------------------------------
    real(8)              :: qn(3, 3), x, y, z, r2
    integer              :: ix, iy, iz, i, isymop, k
!   ----------------------------------------------------------------------------

    qn = 0.0d0

    ix = iptax(1, 1)
    iy = iptax(2, 1)
    iz = iptax(3, 1)

    k = 0
    do i = 1, nucind
      do isymop = 0, maxopr
        if (iand(isymop, istbnu(i)) == 0) then
          k = k + 1

          x  = cord(1, k)
          y  = cord(2, k)
          z  = cord(3, k)

          r2 = x*x + y*y + z*z
          r2 = r2/3.0d0

          qn(ix, ix) = qn(ix, ix) + charge(i)*(x*x - r2)
          qn(iy, iy) = qn(iy, iy) + charge(i)*(y*y - r2)
          qn(iz, iz) = qn(iz, iz) + charge(i)*(z*z - r2)
          qn(ix, iy) = qn(ix, iy) + charge(i)* x*y
          qn(ix, iz) = qn(ix, iz) + charge(i)* x*z
          qn(iy, iz) = qn(iy, iz) + charge(i)* y*z

        end if
      end do
    end do

    qn(iy, ix) = qn(ix, iy)
    qn(iz, ix) = qn(ix, iz)
    qn(iz, iy) = qn(iy, iz)

    qn = 1.5d0*qn

!   ordering is: xx, yx, yy, zx, zy, zx
    q = (/qn(1, 1), qn(1:2, 2), qn(1:3, 3)/)

  end subroutine

  subroutine gradnn_ifc(n, g)

    integer       :: n
    real(8) :: g(3*n)

    iprint = 0
    maxdif = 1

    call nucrep((0.0d0), (0.0d0), (0.0d0))

    g(1:3*n) = gradnn(1:3*n)

  end subroutine

  function do_dft() result(l)

    use dirac_cfg

!   ----------------------------------------------------------------------------
    logical :: l
!   ----------------------------------------------------------------------------

    l = dirac_cfg_dft_calculation

  end function

  subroutine get_F(D, F)

      use xcint_main

    type(matrix), intent(in)    :: D
    type(matrix), intent(inout) :: F
!   ----------------------------------------------------------------------------
    type(matrix)                :: F1, F2
    integer,      parameter     :: file_unit = 70

    call init_mat(F1, ntbas(0), ntbas(0))
    call init_mat(F2, ntbas(0), ntbas(0))

    F1%elms = 0.0d0
    F2%elms = 0.0d0

    call opnfil(file_unit, 'DFFCK1', 'OLD', 'get_F')
    call reafck(file_unit, F1%elms, .true., 1)
    close(file_unit, status = 'keep')

    call opnfil(file_unit, 'DFFCK2', 'OLD', 'get_F')
    call reafck(file_unit, F2%elms, .true., 1)
    close(file_unit, status = 'keep')

    if (.not. isdef(F)) call init_mat(F, ntbas(0), ntbas(0))

    F = F1 + F2

    F1 = 0
    F2 = 0

!   xc contribution to ks matrix
    if (dirac_cfg_dft_calculation) then
        call generate_num_grid(D%elms)
#ifdef VAR_MPI
        if (parcal) call dirac_parctl( XCINT_PAR )
#endif
        call integrate_xc(xc_mat_dim      = ntbas(0), &
                           xc_nz           = nz,       &
                           xc_dmat_0       = D%elms,   &
                           xc_nr_dmat      = 0,        &
                 	    xc_nr_fmat      = 1,        &
                           xc_fmat         = F%elms,   &
                           xc_do_potential = .true.)
    end if

  end subroutine

  subroutine get_G(D, G)

    type(matrix),  intent(in)    :: D
    type(matrix),  intent(inout) :: G
!   ----------------------------------------------------------------------------
    integer                      :: npos_nr
    integer,       allocatable   :: npos(:)

    real(8), pointer   :: work(:)
    integer                      :: lwork
    integer                      :: iprint

    if (.not. isdef(G)) call init_mat(G, ntbas(0), ntbas(0))
    G = 1.0d0*D !inherit symmetry information
    G%elms = 0.0d0

    call get_npos(npos_nr)
    call alloc(npos, npos_nr)

!   fixme:
!   call with proper limit
!   at the moment i don't know the limit
!   simply give entire work and hope for the best
    lwork = len_f77_work
    call di_select_wrk(work, lwork)

    iprint = 0
    call twofck((/G%irep + 1/),  &
!               (/G%ih_sym/),    &
                (/0/),           &
                (/1/),           &
                  G%elms,        &
                  D%elms,        &
                  1,             &
                  npos,          &
                  integral_flag, &
                  iprint,        &
                  work,          &
                  lwork)

    call di_deselect_wrk(work, lwork)

    call dealloc(npos)

  end subroutine

   subroutine read_ibeig(ibeig)

!     --------------------------------------------------------------------------
      integer            :: ibeig(*)
!     --------------------------------------------------------------------------
      integer, parameter :: file_unit = 66
      real(8)            :: energy
!     --------------------------------------------------------------------------

      call opnfil(file_unit, 'DFCOEF', 'OLD', 'read_mo_coef')

      call reacmo(file_unit, &
                          'DFCOEF',  &
                          (/0.0d0/), &
                          (/0.0d0/), &
                          ibeig,     &
                          energy,    &
                          8)

      close(file_unit, status = 'keep')

   end subroutine

   subroutine read_mo_coef(mo_coef)

!     --------------------------------------------------------------------------
      real(8)            :: mo_coef(*)
!     --------------------------------------------------------------------------
      integer, parameter :: file_unit = 66
      real(8)            :: energy
!     --------------------------------------------------------------------------

      call opnfil(file_unit, 'DFCOEF', 'OLD', 'read_mo_coef')

      call reacmo(file_unit, &
                          'DFCOEF',  &
                          mo_coef,   &
                          (/0.0d0/), &
                          (/0/),     &
                          energy,    &
                          2)

      close(file_unit, status = 'keep')

   end subroutine

  subroutine get_C(C, cmo_from_file, i, s, g, u)

!   ----------------------------------------------------------------------------
    type(matrix)        :: C
    real(8), intent(in) :: cmo_from_file(*)
    real(8), intent(in) :: i, s, g, u
!   ----------------------------------------------------------------------------
    integer             :: nr_g_mo_ns
    integer             :: nr_g_mo_pi
    integer             :: nr_g_mo_ps
    integer             :: nr_u_mo_ns
    integer             :: nr_u_mo_pi
    integer             :: nr_u_mo_ps
    integer             :: nr_g_mo
    integer             :: nr_u_mo
    integer             :: nr_g_ao
    integer             :: nr_u_ao
    integer             :: k, l, m, n, ir, iw1, iw2, iz
!   ----------------------------------------------------------------------------

    nr_g_mo_ns = npsh(1)
    nr_g_mo_pi = nish(1)
    nr_g_mo_ps = nesh(1) - nr_g_mo_pi
    nr_g_mo    = nr_g_mo_ns + nr_g_mo_pi + nr_g_mo_ps
    nr_u_mo_ns = npsh(2)
    nr_u_mo_pi = nish(2)
    nr_u_mo_ps = nesh(2) - nr_u_mo_pi
    nr_u_mo    = nr_u_mo_ns + nr_u_mo_pi + nr_u_mo_ps
    nr_g_ao    = nfbas(1, 0)
    nr_u_ao    = nfbas(2, 0)

    if (.not. isdef(C)) call init_mat(C, ntbas(0), norbt)
    C%elms = 0.0d0

    ir = 0

    do iz = 1, nz

      iw1 = (iz - 1)*(nr_g_mo + nr_u_mo)*(nr_g_ao + nr_u_ao)
      k   = 0

      do m = 1, nr_g_mo_ns
        k = k + 1
        do n = 1, nr_g_ao
          ir  = ir + 1
          iw2 = iw1 + (k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = cmo_from_file(ir)*g*s
        end do
      end do

      do m = 1, nr_g_mo_pi
        k = k + 1
        do n = 1, nr_g_ao
          ir  = ir + 1
          iw2 = iw1 + (k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = cmo_from_file(ir)*g*i
        end do
      end do

      do m = 1, nr_g_mo_ps
        k = k + 1
        do n = 1, nr_g_ao
          ir  = ir + 1
          iw2 = iw1 + (k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = cmo_from_file(ir)*g*s
        end do
      end do

    end do

    if (nfsym == 2) then

    do iz = 1, nz

      iw1 = (iz - 1)*(nr_g_mo + nr_u_mo)*(nr_g_ao + nr_u_ao)
      k   = 0

      do m = 1, nr_u_mo_ns
        k = k + 1
        do n = nr_g_ao + 1, nr_g_ao + nr_u_ao
          ir  = ir + 1
          iw2 = iw1 + (nr_g_mo + k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = cmo_from_file(ir)*u*s
        end do
      end do

      do m = 1, nr_u_mo_pi
        k = k + 1
        do n = nr_g_ao + 1, nr_g_ao + nr_u_ao
          ir  = ir + 1
          iw2 = iw1 + (nr_g_mo + k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = cmo_from_file(ir)*u*i
        end do
      end do

      do m = 1, nr_u_mo_ps
        k = k + 1
        do n = nr_g_ao + 1, nr_g_ao + nr_u_ao
          ir  = ir + 1
          iw2 = iw1 + (nr_g_mo + k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = cmo_from_file(ir)*u*s
        end do
      end do

    end do
    end if !if (nfsym == 2) then

  end subroutine

  subroutine get_S(S)

    type(matrix),  intent(inout) :: S
!   ----------------------------------------------------------------------------

    if (.not. isdef(S)) call init_mat(S, ntbas(0), ntbas(0))
    S%elms = 0.0d0

    call gtovlx(S%elms, ssmtrc)

  end subroutine

  subroutine set_orbital_rotation_indices(irep,       &
                                          length,     &
                                          from_block, &
                                          to_block,   &
                                          orbital_rotation_indices)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: irep
    integer, intent(in)  :: length
    integer, intent(in)  :: from_block
    integer, intent(in)  :: to_block
    integer, intent(out) :: orbital_rotation_indices(2, length)
!   ----------------------------------------------------------------------------
    integer, parameter   :: ns = 1
    integer, parameter   :: pi = 2
    integer, parameter   :: ps = 3

    integer, parameter   :: g  = 1
    integer, parameter   :: u  = 2

    logical              :: gu_combine(2, 2)

    integer              :: nr_mo(2, 3)       = 0
    integer              :: range_mo(2, 3, 2) = 0

    integer              :: i, s, k
    integer              :: m, m1, m2
    integer              :: gu, gu1, gu2
    integer              :: block, block1, block2

    logical              :: debug_me = .false.
!   ----------------------------------------------------------------------------

    nr_mo(g, ns) = npsh(g)
    nr_mo(g, pi) = nish(g)
    nr_mo(g, ps) = nesh(g) - nish(g)

    nr_mo(u, ns) = npsh(u)
    nr_mo(u, pi) = nish(u)
    nr_mo(u, ps) = nesh(u) - nish(u)

    k = 0
    do gu = 1, 2
      do block = 1, 3
        do m = 1, nr_mo(gu, block)
          k = k + 1
          if (k > range_mo(gu, block, 1)) range_mo(gu, block, 1) = k
          range_mo(gu, block, 2) = k
        end do
      end do
    end do

    gu_combine = .false.
    if (jbtof(irep, 1) == 2) then
!     ungerade perturbation
      gu_combine(1, 2) = .true.
      gu_combine(2, 1) = .true.
    else
!     gerade perturbation
      gu_combine(1, 1) = .true.
      gu_combine(2, 2) = .true.
    end if

    k = 0
    i = 0
    do gu1 = 1, 2
      do block1 = 1, 3
        do m1 = 1, nr_mo(gu1, block1)
          i = i + 1
          s = 0
          do gu2 = 1, 2
            do block2 = 1, 3
              if (gu_combine(gu1, gu2)         &
                  .and. (block1 == from_block) &
                  .and. (block2 == to_block)) then
                do m2 = 1, nr_mo(gu2, block2)
                  k = k + 1
                  s = s + 1
                  orbital_rotation_indices(1, k) = i
                  orbital_rotation_indices(2, k) = s
                end do
              else
                s = s + nr_mo(gu2, block2)
              end if
            end do
          end do
        end do
      end do
    end do

    if (debug_me) then
      call header('debug orbital_rotation_indices', -1)
      do k = 1, length
        write(*, *) k, orbital_rotation_indices(1, k), &
                 '->', orbital_rotation_indices(2, k)
      end do
    end if

  end subroutine
      SUBROUTINE GET_NPOS(NPOS)
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "cbihr2.h"
#include "blocks.h"
#include "dcbfir.h"
      integer npos
      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. ,.TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      end subroutine

end module
