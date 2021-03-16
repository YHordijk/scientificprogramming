#include "dirac_partask.h"

module fde_evaluators_dirac

! embed modules
   use fde_types
   use fde_cfg
   use fde_data
   use fde_io
   use fde_nadd_derv
   use fde_xcfun_interface
   use fde_max_block_length
   use fde_mpi

! dirac-specific modules
   use numoper_integrals
   use fde_dirac_matrices_integration
   use electrostatic_potential
   use density_eval
   use dirac_ao_eval
   use interface_ao
   use interface_mo_specific

   use memory_allocator
!gosia: temporarily:
   use electrostatic_potential

   implicit none

! these are the interfaces to the dirac-specific code
   public fde_launch_slave_process
   public fde_calculate_elpot
   public fde_calculate_nucpot
   public fde_dirac_set_nz
   public fde_dirac_get_nz
   public fde_dirac_get_isymop
   public fde_dirac_set_isymop
   public fde_dirac_get_ihrmop
   public fde_dirac_set_ihrmop

   public fde_get_elpot_pertden
   public fde_get_pertden

! below we have the routines that actually calculate
   public fde_get_density
   interface fde_get_density
      module procedure fde_dirac_get_density_from_dfcoef
   end interface

!gosia unused   public fde_get_density_matrix

   public fde_qccode_data_interface
   interface fde_qccode_data_interface
      module procedure fde_dirac_data_interface
   end interface

! the matrix for the embedding potential, or its contributions to
! different fock matrices can be calculated in different ways; 
! we take this into account by defining an interface name to the
! code-specific routines
 
   public fde_calculate_emb_pot_mat
   interface fde_calculate_emb_pot_mat
      module procedure fde_dirac_embpot_via_integrator
      module procedure fde_dirac_embpot_via_aoproper
   end interface

   public fde_calculate_emb_linrsp_mat
   interface fde_calculate_emb_linrsp_mat
      module procedure fde_dirac_linrsp_via_integrator
   end interface

   public fde_calculate_interaction_energy
   interface fde_calculate_interaction_energy
      module procedure fde_dirac_emb_interaction_energy
   end interface

   public fde_dirac_embpot_static_energy

   public fde_get_embpot_energy
   interface fde_get_embpot_energy
      module procedure fde_dirac_emb_pot_energy
   end interface

   private 

   integer, save :: nr_boson_irreps
   integer, save :: nba, nba_orb, nba_aux
   integer, save :: dfcoef_unit
   integer, save :: nz = 4
   integer, save :: mat_dim_quat_as_1d
   integer, parameter :: fde_max_nr_mat = 10
   integer, save :: isymop(fde_max_nr_mat)
   integer, save :: ihrmop(fde_max_nr_mat)

   logical, save :: is_dirac_initialized = .false.

   contains

!   ----------------------------------------------------------------------------
      subroutine fde_calculate_elpot(level,grid,gf,myprint,old_scheme)
!   ----------------------------------------------------------------------------
         type(grid_function), intent(inout) :: gf
         type(fde_grid), intent(in)         :: grid
         character(len=4)                   :: level
         real(kind=8), allocatable          :: dmat(:)
         logical, intent(in), optional      :: myprint
         logical, intent(in), optional      :: old_scheme
         logical                            :: printlocal
         logical                            :: old_scheme_local
         type(fde_files)                    :: ftmp
         type(fde_import)                   :: itmp

         if (present(myprint)) then
            printlocal = myprint
         else
            printlocal = .True.
         end if

         if (present(old_scheme)) then
            old_scheme_local = old_scheme
         else
            old_scheme_local = .False.
         end if

         if (old_scheme_local) then
            print *,"NOTE: Using serial code to calculate electrostatic potential"
            call fde_calculate_elpot_via_esp(level,grid,gf,myprint)
            call fde_get_density(level,grid,gf)
         else
            call fde_calculate_elpot_via_integrator(level,grid,gf,myprint)
         endif

      end subroutine fde_calculate_elpot

!   ----------------------------------------------------------------------------
      subroutine fde_calculate_elpot_via_esp(level,grid,gf,myprint)
!   ----------------------------------------------------------------------------
         type(grid_function), intent(inout) :: gf
         type(fde_grid), intent(in)         :: grid
         character(len=4)                   :: level
         real(kind=8), allocatable          :: dmat(:)
         logical, intent(in), optional      :: myprint
         logical                            :: printlocal
         type(fde_files)                    :: ftmp
         type(fde_import)                   :: itmp
         integer :: logunit

         if (present(myprint)) then
            printlocal = myprint
         else
            printlocal = .True.
         endif

         call fde_get_files_info(ftmp)
         call fde_get_import_info(itmp)

         logunit      = ftmp%logfile%unit 

         if (nr_boson_irreps > 1) then
            write (logunit,*) ''
            write (logunit,'(2X,A)') 'FDE electrostatic potential not available'
            write (logunit,'(3X,A)') '    (only C(1) symmetry currently supported)'
            write (logunit,*) ''
            return
         endif

! dmat is allocated and set in the routine below
         call fde_dirac_get_dmat_from_dfcoef(dmat,level,myprint=printlocal)

         call fde_dirac_calculate_elpot(dmat,grid%r,gf%elpot,myprint=printlocal)
         call fde_dirac_calculate_nucpot(dmat,grid%r,gf%nucpot)

! alternatively, we can do things separately 
!        call fde_dirac_calculate_hartree(dmat,grid%r,gf%elpot,myprint=printlocal)
!        call fde_dirac_calculate_nucpot(dmat,grid%r,gf%nucpot,myprint=printlocal)
!        gf%elpot = gf%elpot + gf%nucpot        

         deallocate(dmat)

      end subroutine fde_calculate_elpot_via_esp


!   ----------------------------------------------------------------------------
      subroutine fde_calculate_elpot_via_integrator(level,grid,gf,myprint)
!   ----------------------------------------------------------------------------
         type(grid_function), intent(inout) :: gf
         type(fde_grid), intent(in)         :: grid
         character(len=4)                   :: level
         real(kind=8), allocatable          :: dmat(:)
         logical, intent(in), optional      :: myprint
         logical                            :: printlocal
         type(fde_files)                    :: ftmp
         type(fde_import)                   :: itmp
         integer :: logunit

         if (present(myprint)) then
            printlocal = myprint
         else
            printlocal = .True.
         endif

         call fde_get_files_info(ftmp)
         call fde_get_import_info(itmp)

         if (nr_boson_irreps > 1) then
            write (logunit,*) ''
            write (logunit,'(2X,A)') 'FDE electrostatic potential not available'
            write (logunit,'(3X,A)') '    (only C(1) symmetry currently supported)'
            write (logunit,*) ''
            return
         endif

! dmat is allocated and set in the routine below
         call fde_dirac_get_dmat_from_dfcoef(dmat,level,myprint=printlocal)

         call fde_dirac_set_isymop(isymop)
         call fde_dirac_set_ihrmop(ihrmop)
#ifdef VAR_MPI
         if (fde_mpi_get_nr_proc() > 1) then
            call dirac_parctl(FDE_PAR)
         end if 
#endif
         call fde_dirac_emb_matrices_via_integration(         &
                          fde_mat_dim              = nba_orb, &
                          fde_nz                   = nz,      &
                          fde_dmat_0               = dmat,    &
                          fde_nr_dmat              = 0,       &
                          fde_nr_fmat              = 0,       &
                          fde_do_electrostatic_pot = .true.,  &
                          fde_active_grid          = grid,    &
                          fde_active_gridfunction  = gf)    


      end subroutine fde_calculate_elpot_via_integrator


!   ----------------------------------------------------------------------------
      subroutine fde_calculate_nucpot(level,grid,gf,myprint)
!   ----------------------------------------------------------------------------
         type(grid_function), intent(inout) :: gf
         type(fde_grid), intent(in)         :: grid
         character(len=4)                   :: level
         real(kind=8), allocatable          :: dmat(:)
         logical, intent(in), optional      :: myprint
         logical                            :: printlocal
         type(fde_files)                    :: ftmp
         type(fde_import)                   :: itmp
         integer :: logunit

         if (present(myprint)) then
            printlocal = myprint
         else
            printlocal = .True.
         endif

         call fde_get_files_info(ftmp)
         call fde_get_import_info(itmp)
         call fde_dirac_calculate_nucpot(dmat,grid%r,gf%nucpot,myprint=printlocal)

      end subroutine fde_calculate_nucpot


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_calculate_nucpot(dmat,r,vnuc,myprint)
!   ----------------------------------------------------------------------------
         real(kind=8), intent(in)  :: dmat(:)
         real(kind=8), intent(in)  :: r(:,:)
         real(kind=8), intent(out) :: vnuc(:)
         logical, intent(in), optional :: myprint
         logical :: printlocal
         integer :: i
         integer :: np
         integer :: irrep = 1

         if (present(myprint)) then
            printlocal = myprint
         else
            printlocal = .True.
         endif

         np = size(r,2)
!
! now the nuclear potential
!
         if (printlocal) write (*,*) 'Generating the nuclear potential over the grid'

         call get_esp(np,                      &
                      vnuc,                    &
                      irrep,                   &
                      nba_orb,                 &
                      dmat,                    &
                      r,                       &
                      include_nuc_part=.true., &
                      include_el_part=.false.)

! aspg, 14/10/2015
!       note that the electrostatic, hartree and nuclear potentials in dirac have
!       been multiplied by -1 when exported so that we don't have to multipy the 
!       product potential * density by the electron charge when calculating the
!       energy. with this, the potentials will be in line with those obtained
!       e.g. with adf 

        vnuc = -vnuc

      end subroutine fde_dirac_calculate_nucpot


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_calculate_hartree(dmat,r,vc,myprint)
!   ----------------------------------------------------------------------------
         real(kind=8), intent(in)  :: dmat(:)
         real(kind=8), intent(in)  :: r(:,:)
         real(kind=8), intent(out) :: vc(:)
         logical, intent(in), optional :: myprint
         logical :: printlocal
         integer :: i
         integer :: np
         integer :: irrep = 1

         if (present(myprint)) then
            printlocal = myprint
         else
            printlocal = .True.
         endif

         np = size(r,2)
!
! calculating the hartree potential first
! 
         if (printlocal) write (*,*) 'Generating the Hartree potential over the grid'
!
! aspg, 11/10/2012
! important here: get_esp expects to get a density matrix
!                 that is not scaled wih the occupation numbers
!                 (so like what one would get from calling DENMAT).
!                 
         call get_esp(np,                      &
                      vc,                      &
                      irrep,                   &
                      nba_orb,                 &
                      dmat,                    &
                      r,                       &
                      include_nuc_part=.false.,&
                      include_el_part=.true.)

! aspg, 14/10/2015
!       note that the electrostatic, hartree and nuclear potentials in dirac have
!       been multiplied by -1 when exported so that we don't have to multipy the 
!       product potential * density by the electron charge when calculating the
!       energy. with this, the potentials will be in line with those obtained
!       e.g. with adf 

        vc = -vc

      end subroutine fde_dirac_calculate_hartree


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_calculate_elpot(dmat,r,vc,myprint)
!   ----------------------------------------------------------------------------
         real(kind=8), intent(in)  :: dmat(:)
         real(kind=8), intent(in)  :: r(:,:)
         real(kind=8), intent(out) :: vc(:)
         logical, intent(in), optional :: myprint
         logical :: printlocal
         integer :: i
         integer :: np
         integer :: irrep = 1
         
         if (present(myprint)) then
            printlocal = myprint
         else
            printlocal = .True.
         endif

         np = size(r,2)
!
! calculating the hartree potential first
! 
         if (printlocal) write (*,*) 'Generating the Electrostatic potential over the grid'
!
! aspg, 11/10/2012
! important here: get_esp expects to get a density matrix
!                 that is not scaled wih the occupation numbers
!                 (so like what one would get from calling DENMAT).
!                 
         call get_esp(np,                     &
                      vc,                     &
                      irrep,                  &  
                      nba_orb,                 &
                      dmat,                   &
                      r,                      &
                      include_nuc_part=.true.,&
                      include_el_part=.true.)

! aspg, 14/10/2015
!       note that the electrostatic, hartree and nuclear potentials in dirac have
!       been multiplied by -1 when exported so that we don't have to multipy the 
!       product potential * density by the electron charge when calculating the
!       energy. with this, the potentials will be in line with those obtained
!       e.g. with adf 

        vc = -vc

      end subroutine fde_dirac_calculate_elpot
      
!   ----------------------------------------------------------------------------
      subroutine fde_get_elpot_pertden(ndim, nz, ncomp, fmat)
!     gosia: this subroutine calculates the Coulomb kernel-like coupling contribution
!     to the property gradient
!   ----------------------------------------------------------------------------
         integer                :: ndim
         integer                :: nz
         integer                :: ncomp
         real(8), intent(inout) :: fmat(ndim,ndim,nz,ncomp)
         real(8), pointer       :: pertden(:, :), gridp(:,:), dummy(:, :)
         type(grid_function)    :: gf
         type(fde_grid)         :: grid
         type(fde_files)        :: ftmp
         type(fde_import)       :: itmp
         integer                :: icomp, icol
         integer                :: np
         integer                :: irrep = 1
         integer                :: lwork
         real(8), allocatable   :: work(:)
         integer :: file_unit
         character(len=60) :: file_name
         integer :: iz, i, j, ixyz, ip
         real(8) :: average
         real(8) :: time_start, time_total
         real(8), allocatable :: buffer(:,:), vector(:)
         real(8), external :: second
         logical :: debug_me


         time_start = second()

         call fde_get_files_info(ftmp)
         call fde_get_import_info(itmp)

         np = size(fde_grid_im%r,2)

         nullify(pertden)
         file_name = ftmp%frozen_pert_reorth%name
         file_unit = ftmp%frozen_pert_reorth%unit
         call fde_open_file(file_name,file_unit)
         call read_grid(file_unit,np,48,pertden)
         call fde_close_file(file_unit)

         call legacy_lwork_get(lwork)
         allocate(work(lwork))

         do icomp = 1, ncomp
            allocate(buffer(ndim,ndim))
            buffer = 0.0d0
!           columns correspond only to rho_0 (no spin density, no gradients)
            !icol = 10 + icomp
            icol = icomp
            allocate(vector(np))
            do i = 1, np
              vector(i) = pertden(icol,i)*fde_grid_im%w(i)
            end do
            call get_electronic_mep(np,                &
                                    buffer,            &
                                    ndim,              &
                                    fde_grid_im%r(1:3,1:np),  & 
                                    vector,            &
                                    work, lwork,       &
                                    0,                 &
                                     .true.,  &
                                    .false.)
            deallocate(vector)
            do j = 1, ndim
               do i = 1, j
!                    we have to antisymmetrize it, as fmat ^B is antisymmetric
                     average = 0.5d0*(buffer(i,j) - buffer(j,i))
                     buffer(j,i) = average
                     buffer(i,j) =-average
               end do
            end do

            debug_me = .false.
            if (debug_me) then
               write(*, *) 'Coulomb-coupling term: buffer, icomp = ', icomp
               CALL PRQMAT(buffer,ndim, ndim, ndim, &
                     ndim, 1,(/1, 2, 3, 4/), 6)
               write(*, *) 'Coulomb-coupling term: fmat before daxpy, icomp = ', icomp
               CALL PRQMAT(fmat(1,1,1,icomp),ndim, ndim, ndim, &
                     ndim, nz,(/1, 2, 3, 4/), 6)
            end if

            call daxpy(ndim*ndim, 1.0d0, buffer, 1, fmat(1,1,1,icomp), 1)

            deallocate(buffer)

            if (debug_me) then
               write(*, *) 'Coulomb-coupling term: fmat after daxpy, icomp = ', icomp
               CALL PRQMAT(fmat(1,1,1,icomp),ndim, ndim, ndim, &
                     ndim, nz,(/1, 2, 3, 4/), 6)
            end if

         end do

         deallocate(work)

         if (associated(pertden)) then
            nullify(pertden)
         end if

         time_total = second() - time_start
         call timtxt('  time spent in Coulomb coupling term  =', &
                      time_total, 6)


      end subroutine fde_get_elpot_pertden
 
!   ----------------------------------------------------------------------------
      subroutine fde_dirac_get_dmat_from_dfcoef(dmat,level,myprint)
! the outgoing dmat is not scaled with occupation numbers
!   ----------------------------------------------------------------------------
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "mxcent.h"
#include "nuclei.h"
! Passed/recovered variables
      character(len=4) :: level
      real(kind=8), allocatable, intent(out) :: dmat(:)
      logical, intent(in), optional :: myprint
      logical :: printlocal
! internal variables
      INTEGER :: ierr, i, j
      real(kind=8), allocatable :: cmo(:)
      real(8) :: dummy
      real(kind=8), allocatable :: dmat2(:)
      integer :: ia, ib

      if (present(myprint)) then
         printlocal = myprint
      else
         printlocal = .False.
      endif

      allocate(dmat(mat_dim_quat_as_1d))
      dmat = 0.0d0

      select case(level)
      case('MP2','CCSD')
         if (printlocal) write(lupri,*) 'Reading CC_DENSITY'
         call get_cc_density(level,'AO',DMAT,ierr)
         if (ierr.ne.0) then
            write(lupri,*) 'Error in reading density matrix'
            call fde_quit('error readin density matrix')
         end if
      case default
!     HF density used
         allocate(cmo(ncmotq))
         cmo = 0.0d0

         if (printlocal) write(lupri,*) 'Reading DFCOEF'
         call opnfil(LUCOEF,'DFCOEF','UNKNOWN','EMBDRV')
         call reacmo(lucoef,    &
                            'DFCOEF',  &
                             cmo,       &
                             (/0.0d0/), &
                             (/0/),     &
                             dummy,     &
                             2)
!     SCF construct density matrix
         if (nfmat.gt.1) then 
            allocate(dmat2(mat_dim_quat_as_1d*nfmat))
            dmat2 = 0.0d0
            call denmat(dmat2,cmo,1)
            call gtdoav(nba, nz, nfmat, dmat2, 1)
            dmat = dmat2(1:mat_dim_quat_as_1d)
            deallocate(dmat2)
         else
            call denmat(dmat,cmo,0)
         end if
         deallocate(cmo)
      end select

      end subroutine fde_dirac_get_dmat_from_dfcoef 


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_get_density_from_dfcoef(level,grid,gf)
!   ----------------------------------------------------------------------------
      type(grid_function), intent(inout) :: gf
      type(fde_grid), intent(in)         :: grid
      character(len=4), intent(in)       :: level

      INTEGER :: i
      INTEGER :: NDIMENSION
      INTEGER :: NDER=0, DOLND=0

      REAL(KIND=8),ALLOCATABLE :: GAO(:), GAO1(:),GAO2(:),BUF(:),NCNT(:)
      REAL(KIND=8),ALLOCATABLE :: GAB1(:)

      real(kind=8), allocatable :: dmat(:)

      REAL(KIND=8)::dummy,TOLS,TOLOG
      REAL(KIND=8)::gexp=0.0,ncentc=1.0,factor=1.0

      NDIMENSION=nba_orb

!     DMAT does into account the occupation numbers for scf and mp2 ...
      call fde_dirac_get_dmat_from_dfcoef(dmat,level)
!     so we correct it for scf and mp2...
      dmat = 2.0d0*dmat

      allocate(gao (   ndimension))
      allocate(gao1( 3*ndimension))
      allocate(gao2( 6*ndimension))
      allocate(gab1( 3*ndimension))
      allocate(ncnt(   ndimension))
      allocate(buf ( 4*ndimension))

!     Initialize the matrices to zero.

      GAO  = 0.0
      GAO1 = 0.0
      GAO2 = 0.0
      BUF  = 0.0
      NCNT = 0.0
      GAB1 = 0.0

!     WRITE(*,*) 'Calculating the density and its derivatives.'
      
      do i = 1, size(grid%r,2)

        CALL GETSOS(GAO,GAO1,GAO2,GAB1,NCNT,  &
                grid%r(1,i), grid%r(2,i), grid%r(3,i), &
                    BUF,NDIMENSION,2,DOLND,0)
        CALL GETRHO(gf%n(i),0,GAO,DMAT,BUF)
        CALL GETGRHO(gf%gn(1:3,i),0,GAO,GAO1,DMAT,BUF)
        CALL GETLRHO(gf%hn(1,i),gf%hn(2,i),gf%hn(3,i), &
                     gf%hn(4,i),gf%hn(5,i),gf%hn(6,i), &
                     GAO,GAO1,GAO2,DMAT,BUF)

      END DO

      deallocate (buf)
      deallocate (gao)
      deallocate (gao1)
      deallocate (gao2)
      deallocate (gab1)
      deallocate (ncnt) 
      deallocate (dmat)

      end subroutine fde_dirac_get_density_from_dfcoef

!   ----------------------------------------------------------------------------
      subroutine fde_get_pertden(level,grid,gf)
  use memory_allocator 
#include "priunit.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbham.h"
!   ----------------------------------------------------------------------------
      type(grid_function), intent(inout) :: gf
      type(fde_grid), intent(in)         :: grid
      character(len=4), intent(in)       :: level

      real(8), allocatable :: tbdmat(:,:)
      real(8), allocatable :: dmat(:,:)
      real(8), allocatable :: buffer(:,:)
      real(8), allocatable :: cmo(:)
      integer, allocatable :: ibeig(:)

      integer :: tbdmat_file_unit
      integer :: i, icomp, ndim
      logical :: debug_me
      real(8) :: dum, toterg
      integer :: idum

!     for reorthonormalization part we need the connection matrices
!     there's no easy way to get them, i can either call old subroutines or to use TBMO file
!     it requires that before the shielding calculations are run and this file is saved
!     but it is simpler than anything else.
      
!     read T^B density matrix from TBDMAT file
!     ----------------------------------------
      ndim = nba_orb
      call alloc(tbdmat,  ndim*ndim*nz, 3)
      tbdmat = 0.0d0

      tbdmat_file_unit = 81

      open(tbdmat_file_unit,          &
           file   = 'TBDMAT',         &
           status = 'old',            &
           form   = 'unformatted',    &
           access = 'direct',         &
           recl = 8*ndim*ndim*nz,     &
           action = 'read')

      do i = 1, 3 ! 3 components of B-field
         call readac(tbdmat_file_unit,ndim*ndim*nz,tbdmat(:,i),i)
      end do

      close(tbdmat_file_unit, status='keep')

      debug_me = .false.
      if (debug_me) then
         do i = 1, 3
           write(*, *) 'tbdmat read test, i =', i
           call prqmat(tbdmat(:,i), ndim, ndim, ndim, ndim, nz,  &
                       (/1,2,3,4/), lupri)
         end do
      end if

!     get coefficients from DFCOEF...
      call alloc(cmo,   ncmotq)
      cmo = 0.0d0
      call alloc(ibeig, ndim)
      ibeig = 0

      call opnfil(lucoef,'DFCOEF','OLD','fde_get_pertden')
      if (spinfr) then
         call reacmo(lucoef,'DFCOEF',cmo,dum,ibeig,toterg,10)
      else
         call reacmo(lucoef,'DFCOEF',cmo,dum,idum,toterg,2)
      endif
      close(lucoef,status='keep')

!     ...to construct the density matrix
      call alloc(dmat, ndim*ndim, nz)
      dmat = 0.0d0
      call genden(dmat, cmo, 1, 0)

        
!     then calculate the perturbed density
!     for the moment both contributions to perturbed density (direct and reorthonormalization)
!     have to be done in one go
         call alloc(buffer, ndim*ndim, nz)
         call dcopy(ndim*ndim*nz, dmat, 1, buffer, 1)

#ifdef VAR_MPI
        if (parcal) call dirac_parctl(FDE_PAR)
#endif

          call fde_dirac_emb_matrices_via_integration(     &
                                  fde_mat_dim   = ndim,    &
                                  fde_nz        = nz,      &
                                  fde_dmat_0    = dmat,    & !only real part of dmat
                                  fde_nr_dmat   = 0,       &
                                  fde_nr_fmat   = 0,       &
                                  fde_dmat_pertden_direct      = buffer,  &  !full dmat
                                  fde_dmat_pertden_reorth      = tbdmat,  &
                                  fde_do_export_pertden_reorth = .true.,  &
                                  fde_do_export_pertden_direct = .true.,  &
                                  fde_active_grid = grid,                 &
                                  fde_active_gridfunction = gf)
         call dealloc(buffer)

      call dealloc(tbdmat)
      call dealloc(dmat)
      call dealloc(cmo)
      call dealloc(ibeig)


      end subroutine fde_get_pertden



!   ----------------------------------------------------------------------------
      subroutine fde_dirac_embpot_via_aoproper(comb,iprint)
!   ----------------------------------------------------------------------------
!     the following sets up the matrix elements for ground-state embedding 
!     potential, as a so-called "numerical operator" (that is, we have
!     values of the given operator on a set of gridpoints P, and contract
!     that with the value of the orbitals (or their derivatives, whatever)
!     and the integration weight to get the final integral: 
!
!     M_ab = < phi^n_a | operator | phi^m_b > = \sum^P_i w_i phi^n_a(i) operator(i) phi^m_b(i) 
!
!     where i is the i-th grid point, w_i the integration weight, and
!     phi^{n,m}_{a,b}(i) is the value of the {n,m}-th order derivative of
!     the (contracted) basis function {a,b} at point i.
!
!     information about the ao-basis blocks active for a given operator is
!     carried in the variable COMB, as indicated up. In the case of the
!     embedding potential for FDE, we have a diagonal operator just like
!     for the molecular field.
#include "implicit.h"
#include "priunit.h"
            LOGICAL DOINT(4)
            CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
            character*80 string
#include "dcbbas.h"
#include "mxcent.h"
#include "nuclei.h"

            integer :: file_unit, iprint, i, irrep_vemb
            real(kind=8), allocatable :: a_oneint(:)
            character(len=60) :: file_name
            type(fde_files)   :: ftmp

            call fde_get_files_info(ftmp)
      
            file_unit = ftmp%embpot%unit
            file_name = ftmp%embpot%name
            call fde_open_file(file_name,file_unit)

 
            label = 'FDEVEMB '
            irrep_vemb = 0  ! embedding potential is totally symmetric

            read(comb,'(4L1)',err=1000) (doint(i),i=1,4)

            allocate(a_oneint(nnbbasx)) 
            a_oneint = 0.0d0
            call NumOper_OneElOpMatrix(a_oneint,file_unit,doint,irrep_vemb)

            IF(IPRINT.GE.4) THEN
              write (string,*) 'Integrals over (numerical) operator: '//label
              call HEADER(string,-1)
              call OUTPAK(a_oneint,NTBAS(0),1,6)
            ENDIF
!
!        Generate integral labels
!
            CALL GETDAT(RTNLBL(1),RTNLBL(2))
            RTNLBL(2)(1:2) = 'SY'
            WRITE(RTNLBL(2)(3:4),'(I2)') 1
            RTNLBL(2)(5:8) = COMB
!
!        Write integrals to file
!
            call WRTPRO(a_oneint,NNBBASX,LABEL,RTNLBL,IPRINT)
!
            call fde_close_file(file_unit)
            deallocate(a_oneint)

            return
 1000       continue
!
!     Not able to read DOINT information from COMB
!
            write(LUPRI,'(A,A)') 'FDE_SaveAOPROPER_StaticEmbPot: '// &
                                 'Not able to read COMB =',COMB
            call fde_quit('Vemb2AOPROPER: Not able to read COMB')

         end subroutine fde_dirac_embpot_via_aoproper


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_embpot_via_integrator(mat_dim,dmat_0,fmat)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: mat_dim
         real(kind=8), target, intent(in) :: dmat_0(*)
         real(kind=8), target  :: fmat(*)

         call fde_dirac_emb_matrices_via_integration(      &
                                  fde_mat_dim   = mat_dim, &
                                  fde_nz        = nz,      &
                                  fde_dmat_0    = dmat_0,  &
                                  fde_nr_dmat   = 0,       &
                                  fde_nr_fmat   = 1,       &
                                  fde_fmat      = fmat,    &
                                  fde_do_potential = .true.)
!       call fde_dirac_emb_interaction_energy('DHF ')

      end subroutine fde_dirac_embpot_via_integrator


      subroutine fde_launch_slave_process()
         use fde_dirac_matrices_integration
         implicit none
         call fde_dirac_emb_matrices_via_integration(        &
                                  fde_mat_dim   = 1,         &
                                  fde_nz        = 1,         &
                                  fde_dmat_0    = (/0.0d0/), &
                                  fde_nr_dmat   = 0,         &
                                  fde_nr_fmat   = 0)
      end subroutine


      function fde_dirac_embpot_static_energy()
         real(8)         :: fde_dirac_embpot_static_energy
         real(8)         :: fde_dirac_active_nr_electrons
         integer         :: i

! todo: generalize to other densities
         call fde_dirac_get_density_from_dfcoef('DHF ',fde_grid_sv,gf_active_sv)

         fde_dirac_embpot_static_energy = 0.0d0
         fde_dirac_active_nr_electrons = 0.0d0
         do i = 1, fde_grid_sv%npoints
            fde_dirac_embpot_static_energy = fde_dirac_embpot_static_energy &
                     + fde_grid_sv%w(i) * gf_active_sv%n(i) * fde_static_vemb(i) 
            fde_dirac_active_nr_electrons = fde_dirac_active_nr_electrons + fde_grid_sv%w(i) * gf_active_sv%n(i)
         end do
         print *, "FDE excess energy (<0|v_emb(r)|0>) removed from SCF energy: ",fde_dirac_embpot_static_energy
         print *, "Integrated number of electrons in active subsystem        : ",fde_dirac_active_nr_electrons

      end function


      subroutine fde_dirac_interaction_components(level,verbose,non_additive,el_coulomb_f,el_nuclearAtt_a)
         logical, intent(in)  :: verbose
         character(len=4), intent(in)  :: level
         integer              :: derv_length, max_fun_derv, bllen, i, order
         real(8)              :: threshold
         real(8)              :: df_energy
         real(8), allocatable :: derv_active(:, :)
         real(8), allocatable :: derv_frozen(:, :)
         real(8), allocatable :: derv_total(:, :)
         real(8), intent(out) :: non_additive
         real(8), intent(out) :: el_coulomb_f
         real(8), intent(out) :: el_nuclearAtt_a
         real(8)              :: w(1) 
         real(8)              :: n_0_a(1), gnn_0_a(1)
         real(8)              :: n_0_f(1), gnn_0_f(1)
         real(8)              :: n_0_t(1), gnn_0_t(1)

         threshold    = 1.0e-28
         max_fun_derv = 1
         bllen        = 1
         order        = 1
!        parallel_fde = (fde_mpi_get_nr_proc() > 1) 
         parallel_fde = .false.
         derv_length  = nr_nonzero_derv

         allocate(derv_active(max_block_length, 0:derv_length))
         allocate(derv_frozen(max_block_length, 0:derv_length))
         allocate(derv_total(max_block_length,  0:derv_length))

! aspg: todo for static potential
! problem here is that if using a static potential we don't go over the
! integrator and when broadcasting the functionals there won't be a process
! listening... it's
!         call fde_set_nadd_functionals(order,parallel_fde)

         call fde_dirac_get_density_from_dfcoef(level,fde_grid_im,gf_active)
         call fde_calculate_nucpot(level,fde_grid_im,gf_active,myprint=verbose)

         el_coulomb_f = 0.0d0
         el_nuclearAtt_a = 0.0d0
         non_additive = 0.0d0

         derv_active = 0.0d0
         derv_frozen = 0.0d0
         derv_total  = 0.0d0

         do i = 1, fde_grid_im%npoints

            w(1) = 1.0

            n_0_a(1)   = gf_active%n(i)
            if ( n_0_a(1) .gt. threshold) then
               gnn_0_a(1) = gf_active%gn(1,i) * gf_active%gn(1,i) &
                          + gf_active%gn(2,i) * gf_active%gn(2,i) &
                          + gf_active%gn(3,i) * gf_active%gn(3,i)

               call get_xcke_fun_derv(max_fun_derv, &
                                      bllen,        &
                                      w,            &
                                      n_0_a,        &
                                      gnn_0_a,      &
                                      derv_active)
            else
               n_0_a(1)    = 0.0d0
               gnn_0_a(1)  = 0.0d0
               derv_active = 0.0d0
            endif

            n_0_f(1)   = gf_frozen%n(i)
            if ( n_0_f(1) .gt. threshold) then
               gnn_0_f(1) = gf_frozen%gn(1,i) * gf_frozen%gn(1,i) &
                          + gf_frozen%gn(2,i) * gf_frozen%gn(2,i) &
                          + gf_frozen%gn(3,i) * gf_frozen%gn(3,i)

               call get_xcke_fun_derv(max_fun_derv, &
                                      bllen,        &
                                      w,            &
                                      n_0_f,        &
                                      gnn_0_f,      &
                                      derv_frozen)
            else
               n_0_f(1)    = 0.0d0
               gnn_0_f(1)  = 0.0d0
               derv_frozen = 0.0d0
            endif

            n_0_t(1)   = gf_active%n(i)    + gf_frozen%n(i)
            if ( n_0_t(1) .gt. threshold) then
               gnn_0_t(1) = ((gf_active%gn(1,i) + gf_frozen%gn(1,i))*(gf_active%gn(1,i) + gf_frozen%gn(1,i))) &
                          + ((gf_active%gn(2,i) + gf_frozen%gn(2,i))*(gf_active%gn(2,i) + gf_frozen%gn(2,i))) &
                          + ((gf_active%gn(3,i) + gf_frozen%gn(3,i))*(gf_active%gn(3,i) + gf_frozen%gn(3,i)))

               call get_xcke_fun_derv(max_fun_derv, &
                                      bllen,        &
                                      w,            &
                                      n_0_t,        &
                                      gnn_0_t,      &
                                      derv_total)
            else
               n_0_t(1)    = 0.0d0
               gnn_0_t(1)  = 0.0d0
               derv_total  = 0.0d0
            endif
           df_energy = derv_total(1, d0000000) - derv_active(1, d0000000) - derv_frozen(1, d0000000)
           non_additive     = non_additive    + fde_grid_im%w(i)*df_energy
           el_coulomb_f     = el_coulomb_f    + fde_grid_im%w(i)*gf_frozen%elpot(i)*n_0_a(1)
           el_nuclearAtt_a  = el_nuclearAtt_a - fde_grid_im%w(i)*gf_active%nucpot(i)*n_0_f(1)
         end do

         if (allocated(derv_active)) deallocate(derv_active)
         if (allocated(derv_frozen)) deallocate(derv_frozen)
         if (allocated(derv_total))  deallocate(derv_total)

      end subroutine

      subroutine fde_dirac_emb_interaction_energy(fde_active_energy,level,myprint,recalculate)
         character(len=4), intent(in)  :: level
         logical, intent(in), optional :: myprint
         logical, intent(in), optional :: recalculate
         real(8), intent(in)  :: fde_active_energy 
         logical              :: printlocal
         logical              :: recalculate_local
         integer              :: logunit
         real(8)              :: energy, interaction
         real(8)              :: non_additive
         real(8)              :: electrostatic 
         real(8)              :: el_coulomb_f
         real(8)              :: el_nuclearAtt_a 
         real(8)              :: internuclear 
         type(fde_files)      :: ftmp
         type(fde_import)     :: itmp
         real(8)              :: total_system_energy

         if (present(myprint)) then
            printlocal = myprint
         else
            printlocal = .True.
         endif

         if (present(recalculate)) then
            recalculate_local = recalculate 
         else
            recalculate_local = .False.
         endif

         call fde_get_files_info(ftmp)
         call fde_get_import_info(itmp)

         logunit      = ftmp%logfile%unit

         internuclear = fde_intersub_nucrep

         if (itmp%im_vemb) then 
            WRITE (logunit,'(A)') ''
            write (logunit,'(2X,A)') 'FDE interaction energy contributions not available for static potential'
            WRITE (logunit,'(A)') ''
            return
         endif

         if (.not.itmp%im_frozen) then 
            WRITE (logunit,'(A)') ''
            write (logunit,'(2X,A)') 'FDE interaction energy contributions not available'
            write (logunit,'(3X,A)') '    (frozen density not available)'
            WRITE (logunit,'(A)') ''
            return
         endif

         if (recalculate_local) then
            call fde_dirac_interaction_components(level,printlocal,non_additive,el_coulomb_f,el_nuclearAtt_a)
         else
            non_additive     = fde_dirac_interaction_nadd()
            el_coulomb_f     = fde_dirac_interaction_el_b_rho_a()
            el_nuclearAtt_a  = fde_dirac_interaction_nuc_a_rho_b()
         endif

         electrostatic = el_coulomb_f + el_nuclearAtt_a + internuclear 
         interaction   = electrostatic + non_additive

         WRITE(logunit,'(A)') ''
         write(logunit,'(2X,A)') 'FDE frozen subsystem contributions (imported)'
         write(logunit,'(/3X,A,F22.15)') 'Total energy                             : ',fde_frozen_energy

         write (logunit,'(/,2X,A,/)')    'FDE interaction energy contributions ('//level(1:4)//' active density)'
         write (logunit,'(3X,A,F22.15)') '1. Non-additive  contributions           : ',non_additive
         write (logunit,'(/3X,A)')       '2. Electrostatic contributions'
         write (logunit,'(4X,A,F22.15)') ' a.\int dr [v^B_nuc(r)+v^B_H(r)]rho_A(r): ',el_coulomb_f
         write (logunit,'(4X,A,F22.15)') ' b.           \int dr v^A_nuc(r)rho_B(r): ',el_nuclearAtt_a
         if (internuclear.gt.0.0d0) then
            write (logunit,'(4X,A,F22.15)') ' c. Intersubsystem nuclear repulsion    : ',internuclear
         else
            write (logunit,'(4X,A)')     ' c. Intersubsystem nuclear repulsion    :   not given'
         endif
         write (logunit,'(4X,A,F22.15)') '                                  Total : ',electrostatic
         write (logunit,'(/3X,A,F22.15)') 'Interaction energy (1+2)                 : ',interaction

         if ((internuclear.ne.0.0d0).and.(fde_frozen_energy.ne.0.0d0)) then
            total_system_energy = interaction + fde_frozen_energy + fde_active_energy
            write(logunit,'(A)') ''
            write (logunit,'(3X,A,F22.15)') 'Total system energy                      : ',total_system_energy
         end if

     end subroutine fde_dirac_emb_interaction_energy


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_linrsp_via_integrator(mat_dim,dmat0,ndmat,dmat,fmat)
!   ----------------------------------------------------------------------------
         integer, intent(in)  :: mat_dim
         integer, intent(in)  :: ndmat
         real(kind=8), target, intent(in) :: dmat0(*)
         real(kind=8), target :: dmat(*)
         real(kind=8), target :: fmat(*)

         call fde_dirac_emb_matrices_via_integration(       &
                           fde_mat_dim           = mat_dim, &
                           fde_nz                = nz,      &
                           fde_dmat_0            = dmat0,   &
                           fde_nr_dmat           = ndmat,   &
                           fde_nr_fmat           = ndmat,   &
                           fde_dmat              = dmat,    &
                           fde_fmat              = fmat,    &
                           fde_fmat_pg_sym       = isymop,  &
                           fde_dmat_pg_sym       = isymop,  &
                           fde_dmat_ih_sym       = ihrmop,  &
                           fde_response_order_mo = 1)

      end subroutine fde_dirac_linrsp_via_integrator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! routines to set/recover values for certain dirac-specific variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_set_isymop(value)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: value(fde_max_nr_mat)
         isymop = value
      end subroutine fde_dirac_set_isymop

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_set_ihrmop(value)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: value(fde_max_nr_mat)
         ihrmop = value
      end subroutine fde_dirac_set_ihrmop

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_set_nz(value)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: value
         nz = value
      end subroutine fde_dirac_set_nz


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_get_isymop(value)
!   ----------------------------------------------------------------------------
         integer, intent(out) :: value(fde_max_nr_mat)
         value = isymop
      end subroutine fde_dirac_get_isymop

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_get_ihrmop(value)
!   ----------------------------------------------------------------------------
         integer, intent(out) :: value(fde_max_nr_mat)
         value = ihrmop
      end subroutine fde_dirac_get_ihrmop

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_get_nz(value)
!   ----------------------------------------------------------------------------
        integer, intent(out) :: value
         value = nz
      end subroutine fde_dirac_get_nz


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_data_interface
!   ----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"

         nba_orb = ntbas(0)
         nba = nba_orb
         dfcoef_unit   = LUCOEF
         mat_dim_quat_as_1d = n2bbasxq
         nr_boson_irreps = nbsym

!        call fde_dirac_set_nz(nz)


      end subroutine fde_dirac_data_interface


end module fde_evaluators_dirac
