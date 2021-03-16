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

module fde_data
      
   use fde_types
   use fde_cfg
   use fde_io

   implicit none

   private
!
! routines
!
   public fde_cleanup_data
   public fde_test_frozen_density
   public fde_test_active_density
   public fde_test_export
   public fde_import_frozen
   public fde_import_gridout
   public fde_import_static
   public fde_initialize_import_data

!
! variables 
!
! first, those related to the static embedding potential
!  
   public fde_grid_sv
   public gf_active_sv
   public fde_static_vemb

   integer, save :: fde_static_vemb_grid_np = 0
   logical, save :: static_data_is_initialized = .false.

   type(fde_grid), save :: fde_grid_sv
   type(grid_function), save :: gf_active_sv
   real(kind=8), save, allocatable, target :: fde_static_vemb(:)
!
! importing frozen density, so the active one is needed too
! the grid imported isn't necessarily that of the static 
! potential
!
   public gf_frozen
   public gf_active
   public fde_grid_im
!gosia:
   public import_data_is_initialized

   integer, save :: fde_import_grid_np = 0
   logical, save :: import_data_is_initialized = .false.

   type(grid_function), save :: gf_frozen
   type(grid_function), save :: gf_active
   type(fde_grid), save      :: fde_grid_im
!
! export, using a different grid
!
   public gf_export
   public fde_grid_ex

   logical, save :: export_data_is_initialized = .false.
   integer, save :: fde_export_grid_np = 0

   type(fde_grid), save      :: fde_grid_ex
   type(grid_function), save :: gf_export

   contains


!   ----------------------------------------------------------------------------
   subroutine fde_import_static
!   ----------------------------------------------------------------------------
      integer :: file_unit
      character(len=60) :: file_name
      integer :: i
      real(kind=8), pointer :: fde_vtmp(:), gridp(:,:)
      type(fde_files)  :: ftmp

      call fde_get_files_info(ftmp)
      
      file_unit = ftmp%embpot%unit
      file_name = ftmp%embpot%name

      call fde_open_file(file_name,file_unit)
      call read_grid(file_unit,gridp,fde_vtmp)
      call fde_close_file(file_unit)

      fde_static_vemb_grid_np = size(gridp,2)

      call fde_initialize_static_data

      do i=1, fde_grid_sv%npoints
         fde_static_vemb(i) = fde_vtmp(i)
         fde_grid_sv%r(1,i) = gridp(1,i)
         fde_grid_sv%r(2,i) = gridp(2,i)
         fde_grid_sv%r(3,i) = gridp(3,i)
         fde_grid_sv%w(i)   = gridp(4,i)
      enddo

      write(6,*) "debug, size:",fde_grid_sv%npoints
      write(6,*) "debug, static vemb:",fde_static_vemb(fde_grid_sv%npoints)

      deallocate(fde_vtmp)
      deallocate(gridp)
   end subroutine fde_import_static


   subroutine fde_import_gridout
      integer :: file_unit
      character(len=60) :: file_name
      integer :: i
      real(kind=8), pointer :: gridp(:,:)
      type(fde_files)  :: ftmp

      call fde_get_files_info(ftmp)

      file_unit = ftmp%export%unit
      file_name = ftmp%export%name

      call fde_open_file(file_name,file_unit)
      call read_grid(file_unit,gridp)
      call fde_close_file(file_unit)

      fde_export_grid_np = size(gridp,2)
!aspg print *," np export ",fde_export_grid_np," file "//file_name

      call fde_initialize_export_data

      do i=1, fde_grid_ex%npoints
         fde_grid_ex%r(1,i) = gridp(1,i)
         fde_grid_ex%r(2,i) = gridp(2,i)
         fde_grid_ex%r(3,i) = gridp(3,i)
         fde_grid_ex%w(i)   = gridp(4,i)
      enddo
      deallocate(gridp)
   end subroutine


!   ----------------------------------------------------------------------------
   subroutine fde_import_frozen
!   ----------------------------------------------------------------------------
      integer :: file_unit, file_unit_pert_direct, file_unit_pert_reorth
      character(len=60) :: file_name, file_name_pert_direct, file_name_pert_reorth
      integer :: i
      real(kind=8), pointer :: fde_frozen_prp(:,:), gridp(:,:)
      real(kind=8), pointer :: fde_frozen_prp_pert_direct(:,:), fde_frozen_prp_pert_reorth(:, :)
      type(fde_files)  :: ftmp
      type(fde_import)  :: itmp

      call fde_get_files_info(ftmp)
      call fde_get_import_info(itmp)
      
      file_unit = ftmp%frozen%unit
      file_name = ftmp%frozen%name

      call fde_open_file(file_name,file_unit)
      call read_grid(file_unit,gridp,fde_frozen_prp)
      call fde_close_file(file_unit)

      fde_import_grid_np = size(gridp,2)
      write(*, *) 'fde_import_grid_np 1 = ', fde_import_grid_np

!     "direct" and "reorthonormalization" components of the perturbed density in lao basis
!     are read from respective files
!     number of columns is hardcoded = 4 (for grid) + 48 (for properties)
      if (itmp%im_frozen_pert_direct) then
        file_unit_pert_direct = ftmp%frozen_pert_direct%unit
        file_name_pert_direct = ftmp%frozen_pert_direct%name
        call fde_open_file(file_name_pert_direct,file_unit_pert_direct)
        call read_grid(file_unit_pert_direct,size(gridp,2),48,fde_frozen_prp_pert_direct)
        call fde_close_file(file_unit_pert_direct)
        write(*, *) 'fde: im_frozen_pert_direct file read'
      end if
      if (itmp%im_frozen_pert_reorth) then
        file_unit_pert_reorth = ftmp%frozen_pert_reorth%unit
        file_name_pert_reorth = ftmp%frozen_pert_reorth%name
        call fde_open_file(file_name_pert_reorth,file_unit_pert_reorth)
        call read_grid(file_unit_pert_reorth,size(gridp,2),48,fde_frozen_prp_pert_reorth)
        call fde_close_file(file_unit_pert_reorth)
        write(*, *) 'fde: im_frozen_pert_reorth file read'
      end if


      call fde_initialize_import_data
!g      write(*, *) 'test5 itmp%im_n, itmp%im_n_b = ', itmp%im_n, itmp%im_n_b
      
!aspg print *, "np #1: ",fde_import_grid_np,"  np #2: ",fde_grid_im%npoints
      do i = 1, fde_grid_im%npoints
         gf_frozen%elpot(i)  = fde_frozen_prp(1,i) + fde_frozen_prp(2,i)
!        gf_frozen%nucpot(i) = 0.0
         gf_frozen%n(i)      = fde_frozen_prp(3,i)
         gf_frozen%gn(1,i)   = fde_frozen_prp(4,i)
         gf_frozen%gn(2,i)   = fde_frozen_prp(5,i)
         gf_frozen%gn(3,i)   = fde_frozen_prp(6,i)
         
         if (itmp%im_frozen_pert_direct) then
!        ...s_b(Sigma_xyz, Bxyz, i)
            gf_frozen%s_b_direct(1,1,i)    = fde_frozen_prp_pert_direct(13,i) ! sxBx
            gf_frozen%s_b_direct(2,1,i)    = fde_frozen_prp_pert_direct(14,i) ! syBx
            gf_frozen%s_b_direct(3,1,i)    = fde_frozen_prp_pert_direct(15,i) ! szBx
            gf_frozen%s_b_direct(1,2,i)    = fde_frozen_prp_pert_direct(16,i) ! sxBy
            gf_frozen%s_b_direct(2,2,i)    = fde_frozen_prp_pert_direct(17,i) ! syBy
            gf_frozen%s_b_direct(3,2,i)    = fde_frozen_prp_pert_direct(18,i) ! szBy
            gf_frozen%s_b_direct(1,3,i)    = fde_frozen_prp_pert_direct(19,i) ! sxBz
            gf_frozen%s_b_direct(2,3,i)    = fde_frozen_prp_pert_direct(20,i) ! syBz
            gf_frozen%s_b_direct(3,3,i)    = fde_frozen_prp_pert_direct(21,i) ! szBz
!        ...gs_b(Rxyz, Sigma_xyz, Bxyz, i)
            gf_frozen%gs_b_direct(1,1,1,i) = fde_frozen_prp_pert_direct(22,i) ! dxsxBx
            gf_frozen%gs_b_direct(2,1,1,i) = fde_frozen_prp_pert_direct(23,i) ! dysxBx
            gf_frozen%gs_b_direct(3,1,1,i) = fde_frozen_prp_pert_direct(24,i) ! dzsxBx
            gf_frozen%gs_b_direct(1,2,1,i) = fde_frozen_prp_pert_direct(25,i) ! dxsyBx
            gf_frozen%gs_b_direct(2,2,1,i) = fde_frozen_prp_pert_direct(26,i) ! dysyBx
            gf_frozen%gs_b_direct(3,2,1,i) = fde_frozen_prp_pert_direct(27,i) ! dzsyBx
            gf_frozen%gs_b_direct(1,3,1,i) = fde_frozen_prp_pert_direct(28,i) ! dxszBx
            gf_frozen%gs_b_direct(2,3,1,i) = fde_frozen_prp_pert_direct(29,i) ! dyszBx
            gf_frozen%gs_b_direct(3,3,1,i) = fde_frozen_prp_pert_direct(30,i) ! dzszBx
            gf_frozen%gs_b_direct(1,1,2,i) = fde_frozen_prp_pert_direct(31,i) ! dxsxBy
            gf_frozen%gs_b_direct(2,1,2,i) = fde_frozen_prp_pert_direct(32,i) ! dysxBy
            gf_frozen%gs_b_direct(3,1,2,i) = fde_frozen_prp_pert_direct(33,i) ! dzsxBy
            gf_frozen%gs_b_direct(1,2,2,i) = fde_frozen_prp_pert_direct(34,i) ! dxsyBy
            gf_frozen%gs_b_direct(2,2,2,i) = fde_frozen_prp_pert_direct(35,i) ! dysyBy
            gf_frozen%gs_b_direct(3,2,2,i) = fde_frozen_prp_pert_direct(36,i) ! dzsyBy
            gf_frozen%gs_b_direct(1,3,2,i) = fde_frozen_prp_pert_direct(37,i) ! dxszBy
            gf_frozen%gs_b_direct(2,3,2,i) = fde_frozen_prp_pert_direct(38,i) ! dyszBy
            gf_frozen%gs_b_direct(3,3,2,i) = fde_frozen_prp_pert_direct(39,i) ! dzszBy
            gf_frozen%gs_b_direct(1,1,3,i) = fde_frozen_prp_pert_direct(40,i) ! dxsxBz
            gf_frozen%gs_b_direct(2,1,3,i) = fde_frozen_prp_pert_direct(41,i) ! dysxBz
            gf_frozen%gs_b_direct(3,1,3,i) = fde_frozen_prp_pert_direct(42,i) ! dzsxBz
            gf_frozen%gs_b_direct(1,2,3,i) = fde_frozen_prp_pert_direct(43,i) ! dxsyBz
            gf_frozen%gs_b_direct(2,2,3,i) = fde_frozen_prp_pert_direct(44,i) ! dysyBz
            gf_frozen%gs_b_direct(3,2,3,i) = fde_frozen_prp_pert_direct(45,i) ! dzsyBz
            gf_frozen%gs_b_direct(1,3,3,i) = fde_frozen_prp_pert_direct(46,i) ! dxszBz
            gf_frozen%gs_b_direct(2,3,3,i) = fde_frozen_prp_pert_direct(47,i) ! dyszBz
            gf_frozen%gs_b_direct(3,3,3,i) = fde_frozen_prp_pert_direct(48,i) ! dzszBz
         end if !itmp%im_frozen_pert_direct

         if (itmp%im_frozen_pert_reorth) then
!           n_b(Bxyz, i)
            gf_frozen%n_b_reorth(1,i)      = fde_frozen_prp_pert_reorth(1,i)
            gf_frozen%n_b_reorth(2,i)      = fde_frozen_prp_pert_reorth(2,i)
            gf_frozen%n_b_reorth(3,i)      = fde_frozen_prp_pert_reorth(3,i)
!g          write(*, *) 'test6 inside import', gf_frozen%n(i), gf_frozen%n_b(1, i)
!        ...gn_b(Rxyz, Bxyz, i)
            gf_frozen%gn_b_reorth(1,1,i)   = fde_frozen_prp_pert_reorth(4,i) ! dxBx
            gf_frozen%gn_b_reorth(2,1,i)   = fde_frozen_prp_pert_reorth(5,i) ! dyBx
            gf_frozen%gn_b_reorth(3,1,i)   = fde_frozen_prp_pert_reorth(6,i) ! dzBx
            gf_frozen%gn_b_reorth(1,2,i)   = fde_frozen_prp_pert_reorth(7,i) ! dxBy
            gf_frozen%gn_b_reorth(2,2,i)   = fde_frozen_prp_pert_reorth(8,i) ! dyBy
            gf_frozen%gn_b_reorth(3,2,i)   = fde_frozen_prp_pert_reorth(9,i) ! dzBy
            gf_frozen%gn_b_reorth(1,3,i)   = fde_frozen_prp_pert_reorth(10,i) ! dxBz
            gf_frozen%gn_b_reorth(2,3,i)   = fde_frozen_prp_pert_reorth(11,i) ! dyBz
            gf_frozen%gn_b_reorth(3,3,i)   = fde_frozen_prp_pert_reorth(12,i) ! dzBz
!        ...s_b(Sigma_xyz, Bxyz, i)
            gf_frozen%s_b_reorth(1,1,i)    = fde_frozen_prp_pert_reorth(13,i) ! sxBx
            gf_frozen%s_b_reorth(2,1,i)    = fde_frozen_prp_pert_reorth(14,i) ! syBx
            gf_frozen%s_b_reorth(3,1,i)    = fde_frozen_prp_pert_reorth(15,i) ! szBx
            gf_frozen%s_b_reorth(1,2,i)    = fde_frozen_prp_pert_reorth(16,i) ! sxBy
            gf_frozen%s_b_reorth(2,2,i)    = fde_frozen_prp_pert_reorth(17,i) ! syBy
            gf_frozen%s_b_reorth(3,2,i)    = fde_frozen_prp_pert_reorth(18,i) ! szBy
            gf_frozen%s_b_reorth(1,3,i)    = fde_frozen_prp_pert_reorth(19,i) ! sxBz
            gf_frozen%s_b_reorth(2,3,i)    = fde_frozen_prp_pert_reorth(20,i) ! syBz
            gf_frozen%s_b_reorth(3,3,i)    = fde_frozen_prp_pert_reorth(21,i) ! szBz
!        ...gs_b(Rxyz, Sigma_xyz, Bxyz, i)
            gf_frozen%gs_b_reorth(1,1,1,i) = fde_frozen_prp_pert_reorth(22,i) ! dxsxBx
            gf_frozen%gs_b_reorth(2,1,1,i) = fde_frozen_prp_pert_reorth(23,i) ! dysxBx
            gf_frozen%gs_b_reorth(3,1,1,i) = fde_frozen_prp_pert_reorth(24,i) ! dzsxBx
            gf_frozen%gs_b_reorth(1,2,1,i) = fde_frozen_prp_pert_reorth(25,i) ! dxsyBx
            gf_frozen%gs_b_reorth(2,2,1,i) = fde_frozen_prp_pert_reorth(26,i) ! dysyBx
            gf_frozen%gs_b_reorth(3,2,1,i) = fde_frozen_prp_pert_reorth(27,i) ! dzsyBx
            gf_frozen%gs_b_reorth(1,3,1,i) = fde_frozen_prp_pert_reorth(28,i) ! dxszBx
            gf_frozen%gs_b_reorth(2,3,1,i) = fde_frozen_prp_pert_reorth(29,i) ! dyszBx
            gf_frozen%gs_b_reorth(3,3,1,i) = fde_frozen_prp_pert_reorth(30,i) ! dzszBx
            gf_frozen%gs_b_reorth(1,1,2,i) = fde_frozen_prp_pert_reorth(31,i) ! dxsxBy
            gf_frozen%gs_b_reorth(2,1,2,i) = fde_frozen_prp_pert_reorth(32,i) ! dysxBy
            gf_frozen%gs_b_reorth(3,1,2,i) = fde_frozen_prp_pert_reorth(33,i) ! dzsxBy
            gf_frozen%gs_b_reorth(1,2,2,i) = fde_frozen_prp_pert_reorth(34,i) ! dxsyBy
            gf_frozen%gs_b_reorth(2,2,2,i) = fde_frozen_prp_pert_reorth(35,i) ! dysyBy
            gf_frozen%gs_b_reorth(3,2,2,i) = fde_frozen_prp_pert_reorth(36,i) ! dzsyBy
            gf_frozen%gs_b_reorth(1,3,2,i) = fde_frozen_prp_pert_reorth(37,i) ! dxszBy
            gf_frozen%gs_b_reorth(2,3,2,i) = fde_frozen_prp_pert_reorth(38,i) ! dyszBy
            gf_frozen%gs_b_reorth(3,3,2,i) = fde_frozen_prp_pert_reorth(39,i) ! dzszBy
            gf_frozen%gs_b_reorth(1,1,3,i) = fde_frozen_prp_pert_reorth(40,i) ! dxsxBz
            gf_frozen%gs_b_reorth(2,1,3,i) = fde_frozen_prp_pert_reorth(41,i) ! dysxBz
            gf_frozen%gs_b_reorth(3,1,3,i) = fde_frozen_prp_pert_reorth(42,i) ! dzsxBz
            gf_frozen%gs_b_reorth(1,2,3,i) = fde_frozen_prp_pert_reorth(43,i) ! dxsyBz
            gf_frozen%gs_b_reorth(2,2,3,i) = fde_frozen_prp_pert_reorth(44,i) ! dysyBz
            gf_frozen%gs_b_reorth(3,2,3,i) = fde_frozen_prp_pert_reorth(45,i) ! dzsyBz
            gf_frozen%gs_b_reorth(1,3,3,i) = fde_frozen_prp_pert_reorth(46,i) ! dxszBz
            gf_frozen%gs_b_reorth(2,3,3,i) = fde_frozen_prp_pert_reorth(47,i) ! dyszBz
            gf_frozen%gs_b_reorth(3,3,3,i) = fde_frozen_prp_pert_reorth(48,i) ! dzszBz
         end if !itmp%im_frozen_pert_reorth

         fde_grid_im%r(1,i)  = gridp(1,i)
         fde_grid_im%r(2,i)  = gridp(2,i)
         fde_grid_im%r(3,i)  = gridp(3,i)
         fde_grid_im%w(i)    = gridp(4,i)
      enddo

!gosia comment: these two are allocated inside read_grid...
      if (itmp%im_frozen_pert_direct) deallocate(fde_frozen_prp_pert_direct) 
      if (itmp%im_frozen_pert_reorth) deallocate(fde_frozen_prp_pert_reorth) 
      deallocate(fde_frozen_prp) 
      deallocate(gridp)
   end subroutine fde_import_frozen


!   ----------------------------------------------------------------------------
   subroutine fde_initialize_export_data
!   ----------------------------------------------------------------------------
      if (.not.export_data_is_initialized) then
         export_data_is_initialized = .true.
         if (fde_export_grid_np > 0) then

            call new_fde_grid(fde_grid_ex,fde_export_grid_np)
            call new_grid_function(gf_export,fde_export_grid_np)
         else
            call fde_quit('there are no grid points for fde export!')
         endif
      endif

   end subroutine fde_initialize_export_data


!   ----------------------------------------------------------------------------
   subroutine fde_initialize_static_data
!   ----------------------------------------------------------------------------
      if (.not.static_data_is_initialized) then
         static_data_is_initialized = .true.

         if (fde_static_vemb_grid_np > 0) then

            call new_fde_grid(fde_grid_sv,fde_static_vemb_grid_np)
            call new_grid_function(gf_active_sv,fde_static_vemb_grid_np)

            if (.not.allocated(fde_static_vemb)) then
               allocate(fde_static_vemb(fde_static_vemb_grid_np))
               fde_static_vemb         = 0.0d0
            endif
         else
            call fde_quit('there are no grid points for vemb!')
         endif
      endif
  end subroutine fde_initialize_static_data

!   ----------------------------------------------------------------------------
   subroutine fde_initialize_import_data
!   ----------------------------------------------------------------------------
      if (.not.import_data_is_initialized) then
         import_data_is_initialized = .true.
         if (fde_import_grid_np > 0) then
            call new_fde_grid(fde_grid_im,fde_import_grid_np)
            call new_grid_function(gf_frozen,fde_import_grid_np)
            call new_grid_function(gf_active,fde_import_grid_np)
         else
            call fde_quit('there are no grid points for fde import !')
         endif
      endif
      !write(*, *) 'test3 ', gf_frozen%n_b 
   end subroutine fde_initialize_import_data


!   ----------------------------------------------------------------------------
   subroutine fde_cleanup_data
!   ----------------------------------------------------------------------------
      if (static_data_is_initialized) then
         deallocate(fde_static_vemb)
         call del_fde_grid(fde_grid_sv)
         call del_grid_function(gf_active_sv)
      endif

      if (import_data_is_initialized) then
         call del_grid_function(gf_frozen)
         call del_grid_function(gf_active)
         call del_fde_grid(fde_grid_im)
      endif

      if (export_data_is_initialized) then
         call del_grid_function(gf_export)
         call del_fde_grid(fde_grid_ex)
      endif

   end subroutine fde_cleanup_data


!   ----------------------------------------------------------------------------
   subroutine fde_test_frozen_density
!   ----------------------------------------------------------------------------
      integer      :: i
      real(kind=8) :: frz_particle_nr = 0.0

      do i = 1, fde_grid_im%npoints
         frz_particle_nr = frz_particle_nr + fde_grid_im%w(i)*gf_frozen%n(i) 
      end do
      
      write (*,'(A,F16.6,A)') '  Frozen density integrates to ',frz_particle_nr,' electrons'
      write (*,'(A,I12)') '  Grid points processed ',fde_grid_im%npoints

      frz_particle_nr = 0.0

   end subroutine fde_test_frozen_density


!   ----------------------------------------------------------------------------
   subroutine fde_test_active_density
!   ----------------------------------------------------------------------------
      integer      :: i 
      real(kind=8) :: act_particle_nr 

      act_particle_nr = 0.0d0

      do i = 1, fde_grid_im%npoints
         act_particle_nr = act_particle_nr + fde_grid_im%w(i)*gf_active%n(i) 
      end do
      
      write (*,*) ' Active density integrates to ',act_particle_nr,' electrons'
      write (*,*) ' Grid points processed ',fde_grid_im%npoints

      act_particle_nr = 0.0

   end subroutine fde_test_active_density


!   ----------------------------------------------------------------------------
   subroutine fde_test_export
!   ----------------------------------------------------------------------------
      integer      :: i
      real(kind=8) :: exp_int_n 
      real(kind=8) :: exp_int_ve 
      real(kind=8) :: exp_int_vn 
      real(kind=8) :: exp_int_vh 

      exp_int_n  = 0.0d0
      exp_int_ve = 0.0d0  
      exp_int_vn = 0.0d0  
      exp_int_vh = 0.0d0  

! aspg, 14/10/2015
!       note that the electrostatic, hartree and nuclear potentials in dirac have
!       been multiplied by -1 when exported so that we don't have to multipy the 
!       product potential * density by the electron charge when calculating the
!       energy. with this, the potentials will be in line with those obtained
!       e.g. with adf 
      do i = 1, fde_grid_ex%npoints
         exp_int_n  = exp_int_n  +       fde_grid_ex%w(i)*gf_export%n(i) 
         exp_int_ve = exp_int_ve +       fde_grid_ex%w(i)*gf_export%n(i)*gf_export%elpot(i)
         exp_int_vn = exp_int_vn +       fde_grid_ex%w(i)*gf_export%n(i)*gf_export%nucpot(i) 
      end do
      exp_int_vh = 0.5d0*(exp_int_ve - exp_int_vn)
      
      write (*,*) ' '
      write (*,'(A,I12,A)') ' Exported quantities (over ',fde_grid_ex%npoints,' grid points) integrate to:'
      write (*,'(A,F16.6,A)') ' - electron density    : ',exp_int_n,' electrons'
      write (*,'(A,F16.6,A)') ' - electrostatic energy: ',exp_int_ve,' a.u.'
      write (*,'(A,F16.6,A)') ' - n.-el. attr. energy : ',exp_int_vn,' a.u.'
      write (*,'(A,F16.6,A)') ' - Hartree      energy : ',exp_int_vh,' a.u.'

      exp_int_n  = 0.0d0
      exp_int_ve = 0.0d0  
      exp_int_vn = 0.0d0  
      exp_int_vh = 0.0d0  

   end subroutine fde_test_export

end module fde_data
