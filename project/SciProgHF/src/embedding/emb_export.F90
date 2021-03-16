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

module fde_export_data

   use fde_cfg
   use fde_types
   use fde_data
   use fde_io
   use xml_parser
   use xml_file
   use fde_mag_cfg
   use fde_evaluators_dirac
   
   implicit none

   public fde_export_to_file

   private
   
   contains

!   ----------------------------------------------------------------------------
      subroutine fde_export_to_file(level)
!   ----------------------------------------------------------------------------
         character(len=4) :: level
         character(len=4) :: export_format
         type(fde_export) :: etmp

         call fde_get_export_info(etmp)
         
         if (etmp%do_grid_out .and. (level.eq.etmp%ex_level)) then
            write (*,*) 'Preparing FDE data for export'

! aspg, 14/10/2015
!       note that the electrostatic, hartree and nuclear potentials in dirac have
!       been multiplied by -1 when exported so that we don't have to multipy the 
!       product potential * density by the electron charge when calculating the
!       energy. with this, the potentials will be in line with those obtained
!       e.g. with adf 
            call fde_prepare_export(level,fde_grid_ex,gf_export)

            call fde_get_export_info(etmp) 
            if (etmp%ex_format.eq.'XML') then
               call fde_export_data_as_xml(fde_grid_ex,gf_export)
            else if (etmp%ex_format.eq.'TXT') then
               call fde_export_data_as_txt(fde_grid_ex,gf_export)
            else
               call fde_quit('unrecognized export format')
            endif
         else
               write (*,*) 'inconsistent level of theory in fde_export'
               write (*,*) 'skipping'
         endif
      end subroutine fde_export_to_file


!   ----------------------------------------------------------------------------
      subroutine fde_prepare_export(level,grid,gf)
!   ----------------------------------------------------------------------------
         type(grid_function), intent(inout) :: gf
         type(fde_grid), intent(in)         :: grid
         character(len=4)                   :: level

        call fde_calculate_elpot(level,grid,gf,old_scheme=fde_export_old_esp)

        if (fde_rsp_mag_lao_export) then
           call fde_get_pertden(level,grid,gf)
        end if

      end subroutine fde_prepare_export


!   ----------------------------------------------------------------------------
      subroutine fde_export_data_as_xml(grid,gf)
!   ----------------------------------------------------------------------------
         type(grid_function), intent(in) :: gf
         type(fde_grid), intent(in)      :: grid 

         type(fde_files)  :: ftmp
         type(fde_export) :: etmp

         type(xml_tag),pointer :: tag

         integer :: i
         integer :: export_unit, log_unit
         character(len=60)  :: export_name, log_name

         character(len=132) :: csizept1
         character(len=132) :: csizept2
         character(len=132) :: csizegrdrho1

         real(kind=8), pointer :: gp(:,:)                   ! the grid
         real(kind=8), pointer :: np(:), gnp(:,:), hes(:,:) ! density, its gradient and hessian
         real(kind=8), pointer :: vcp(:), vnp(:)            ! electrostatic (hartree+nuclear), nuclear potentials


         if (associated(grid%r)) then                

            allocate(gp(4,grid%npoints))

            do i = 1, grid%npoints
               gp(1,i) = grid%r(1,i)
               gp(2,i) = grid%r(2,i)
               gp(3,i) = grid%r(3,i)
               gp(4,i) = grid%w(i)
            enddo
         else
            call fde_quit('grid not available in fde export')
         endif

         if (associated(gf%n)) then
            np  => gf%n
         else
            call fde_quit('density not available in fde export')
         endif
         
         if (associated(gf%gn)) then
            gnp => gf%gn
         else
            call fde_quit('density gradient not available in fde export')
         endif   

         if (associated(gf%hn)) then
            hes => gf%hn
         else
            call fde_quit('density gradient not available in fde export')
         endif   
         
         if (associated(gf%elpot)) then
            vcp => gf%elpot
            vnp => gf%nucpot
         else
            call fde_quit('potential not available in fde export')
         endif
         
         call fde_get_files_info(ftmp)
         call fde_get_export_info(etmp)
      
         export_unit = ftmp%export%unit
         export_name = ftmp%export%name

         log_unit    = ftmp%logfile%unit
         log_name    = ftmp%logfile%name
         

         WRITE(log_unit,*) 'Output FDE data to XML file:',trim(export_name)

         call fde_test_export

         write (csizept1,*)     size(gp,1) 
         write (csizept2,*)     size(gp,2)
         write (csizegrdrho1,*) size(gnp,1)

         csizept1     = adjustl(csizept1)
         csizept2     = adjustl(csizept2)
         csizegrdrho1 = adjustl(csizegrdrho1)

         tag=>xml_open(export_name)
           tag=>xml_tag_open(tag,'grid')
             call set_attribute(tag,'size',csizept2)
             tag=>xml_tag_open(tag,'dataset','gridpoints')
               call set_attribute(tag,'size',  csizept2)
               call set_attribute(tag,'width', csizept1)
               call xml_add_data(tag,gp)
             tag=>xml_tag_close(tag)

             if (etmp%ex_n) then
                tag=>xml_tag_open(tag,'dataset','density')
                  call set_attribute(tag,'size',  csizept2)
                  call set_attribute(tag,'width', '1')
                  call xml_add_data(tag,np)
                tag=>xml_tag_close(tag)
             endif

             if (etmp%ex_gn) then
                tag=>xml_tag_open(tag,'dataset','gradient')
                  call set_attribute(tag,'size', csizept2) 
                  call set_attribute(tag,'width', '3')
                  call xml_add_data(tag,gnp)
                tag=>xml_tag_close(tag)
             endif

             tag=>xml_tag_open(tag,'dataset','hessian')
               call set_attribute(tag,'size',  csizept2)
               call set_attribute(tag,'width', '6')    
               call xml_add_data(tag,hes)
             tag=>xml_tag_close(tag)

             if (etmp%ex_coulomb) then
                tag=>xml_tag_open(tag,'dataset','vc')
                  call set_attribute(tag,'size',  csizept2)
                  call set_attribute(tag,'width', '1')
                  call xml_add_data(tag,vcp)
                tag=>xml_tag_close(tag)

                tag=>xml_tag_open(tag,'dataset','nuc')
                  call set_attribute(tag,'size',  csizept2)
                  call set_attribute(tag,'width', '1')
                  call xml_add_data(tag,vnp)
                tag=>xml_tag_close(tag)
             endif

           tag=>xml_tag_close(tag)

           write(log_unit,*) 'done writing xml file'       
           call xml_write(export_name)

           write(log_unit,*) 'closing xml file'
         tag=>xml_close(export_name)

         deallocate(gp)
      end subroutine fde_export_data_as_xml


!   ----------------------------------------------------------------------------
      subroutine fde_export_data_as_txt(grid,gf)
!   ----------------------------------------------------------------------------
         type(grid_function), intent(in) :: gf
         type(fde_grid), intent(in)      :: grid 

         type(fde_files)  :: ftmp
         type(fde_export) :: etmp

         integer :: i, j, k, l
         integer :: nr_variables
         real(kind=8) :: v_pc_array = 0.0d0
         integer :: export_unit, log_unit
         character(len=60)  :: export_name, log_name
       
         nr_variables = 6

         if (.not.associated(grid%r)) &
            call fde_quit('grid not available in fde export')

         if (.not.associated(gf%n)) &
            call fde_quit('density not available in fde export')
         
         if (.not.associated(gf%gn)) &
            call fde_quit('density gradient not available in fde export')

         if (.not.associated(gf%hn)) &
            call fde_quit('density gradient not available in fde export')
         
         if (.not.associated(gf%elpot)) &
            call fde_quit('potential not available in fde export')

!        lao contribs
         !if (fde_rsp_mag_lao_export) then
         !   if (.not.associated(gf%n_b)) then
         !      call fde_quit('n_b not available in fde export')
         !   else
         !      nr_variables = nr_variables + 3
         !   end if
         !   if (.not.associated(gf%gn_b)) then
         !      call fde_quit('gn_b not available in fde export')
         !   else
         !      nr_variables = nr_variables + 9
         !   end if
         !   if (.not. fde_cfg_no_sdft) then
         !      if (.not.associated(gf%s_b)) then
         !         call fde_quit('s_b not available in fde export')
         !      else
         !         nr_variables = nr_variables + 9
         !      end if
         !      if (.not.associated(gf%gs_b)) then
         !         call fde_quit('gs_b not available in fde export')
         !      else
         !         nr_variables = nr_variables + 27
         !      end if
         !   end if
         !end if
         
         call fde_get_files_info(ftmp)
         call fde_get_export_info(etmp)
      
         export_unit = ftmp%export%unit
         export_name = ftmp%export%name

         log_unit    = ftmp%logfile%unit
         log_name    = ftmp%logfile%name

         WRITE(log_unit,*) 'Output FDE data to text file:',trim(export_name)

         call fde_test_export

         call fde_open_file(export_name,export_unit)
         write (export_unit,'(2I8)') grid%npoints, nr_variables

         do i = 1, grid%npoints
!             if (fde_rsp_mag_lao_export.and..not.fde_cfg_no_sdft) then                   ! order of columns:
!               write (export_unit,'(58ES30.20E3)') (grid%r(j,i), j = 1,3), &                 !1-3,   x-y-z
!                                               grid%w(i),             &                      !4
!                                               gf%elpot(i),          &                       !5
!                                               v_pc_array,            &                      !6
!                                               gf%n(i),               &                      !7
!                                               (gf%gn(j,i), j = 1,3), &                      !8-10,  gx-gy-gz
!                                               (gf%n_b(j,i), j = 1,3),             &         !11-13, bx-by-bz
!                                               ((gf%gn_b(j,k,i), j = 1,3), k=1,3), &         !14-22, gxbx-gybx-gzbx - gxby-gyby-gzby - gzbx-gzby-gzbz
!                                               ((gf%s_b(j,k,i), j = 1,3), k=1,3),  &         !23-31, sxbx-sybx-szbx - sxby-syby-szby - szbx-szby-szbz
!                                               (((gf%gs_b(j,k,l,i), j = 1,3), k=1,3), l=1,3) !32-58, gxsxbx-gysxbx-gzsxbx - gxsybx-gysybx-gzsybx - etc.
!             else if (fde_rsp_mag_lao_export.and.fde_cfg_no_sdft) then
!               write (export_unit,'(22ES30.20E3)') (grid%r(j,i), j = 1,3), &
!                                               grid%w(i),             &
!                                               gf%elpot(i),          &
!                                               v_pc_array,            &
!                                               gf%n(i),               &
!                                               (gf%gn(j,i), j = 1,3), &
!                                               (gf%n_b(j,i), j = 1,3), &
!                                               ((gf%gn_b(j,k,i), j = 1,3), k=1,3)
!             else
               write (export_unit,'(10ES30.20E3)') (grid%r(j,i), j = 1,3), &
                                               grid%w(i),             &
                                               gf%elpot(i),          &
                                               v_pc_array,            &
                                               gf%n(i),               &
                                               (gf%gn(j,i), j = 1,3)
!             end if

         end do

         write(log_unit,*) 'done writing plaintext file'       
         write(log_unit,*) 'closing plaintext file'
         call fde_close_file(export_unit)

      end subroutine fde_export_data_as_txt
      

end module fde_export_data 
