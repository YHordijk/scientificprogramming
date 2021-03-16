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
!
!
! module containing the functionality to perform the picture-change
! transformation of a four-component (effective) one-electron operator 
! in SA-AO basis to an infinite-order two-component operator 
! following the X2C formalism. 
!
! written by sknecht april 2010
!
module x2c_h1inf2c_mo2ao

  use x2c_fio
  use x2c_utils, only:                 &
      construct_saomo_overlap_half_tr, &
      print_x2cmat
  use picture_change_operations

  implicit none

  public get_h1_inf_2c_saao_bas

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0

contains

!**********************************************************************
  subroutine get_h1_inf_2c_saao_bas(h1_2cinf_saao,          &
                                    scr_mat1,               &
                                    scr_mat2,               &
                                    scr_mat3,               &
                                    norb_dim,               &
                                    nesh_dim,               &
                                    naosh_L,                &
                                    naosh_all,              &
                                    nr_ao_L,                &
                                    naosh_ls,               &
                                    ioff_aomx,              &
                                    nfsym,                  &
                                    nz,                     &
                                    nzt,                    &
                                    ipq_off_in,             &
                                    is_final_ham_lvl,       &
                                    add_amfi_contrib,       &
                                    is_len_wrk,             &
                                    x2c_file_scr,           &
                                    x2c_file_glb,           &
                                    x2c_file_amf,           &
                                    op_bs_to_fs,            &
                                    spherical_on,           &
                                    fragment_x2c,           &
                                    basic_h1mat,            &
                                    print_lvl)
!**********************************************************************
!
!    purpose: construct the (basic) 2c-h1 matrix in DIRAC-sorted SA-AO basis (h1_infinite_2c_saao) 
!             starting from the operator in matrix form in 2c orthonormal MO basis. 
!
!             save the final h1 matrix in DIRAC-sorted SA-AO basis either 
!
!             - on file X2Camfi_scr if atomic-mean-field contributions are to be added or
!             - on file X2CMAT if no   atomic-mean-field contributions are to be added
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: h1_2cinf_saao(*)
     real(8), intent(inout) :: scr_mat1(*)
     real(8), intent(inout) :: scr_mat2(*)
     real(8), intent(inout) :: scr_mat3(*)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: nesh_dim(nfsym)
     integer, intent(in)    :: naosh_L(nfsym)
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: nr_ao_L
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: nzt
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(inout) :: is_final_ham_lvl
     integer, intent(in)    :: basic_h1mat
     integer, intent(in)    :: add_amfi_contrib
     integer, intent(in)    :: is_len_wrk
     integer, intent(in)    :: x2c_file_scr
     integer, intent(in)    :: x2c_file_glb
     integer, intent(in)    :: x2c_file_amf
     integer, intent(in)    :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: spherical_on
     logical, intent(in)    :: fragment_x2c
!----------------------------------------------------------------------
     real(8), allocatable   :: wrk(:)
     integer                :: lwrk
     integer                :: i
     integer                :: n1esh_f
     integer                :: nbas1_f
     integer                :: nr_sph
     integer                :: ioff_h1_inf_2cmoON
     integer                :: ioff_tramat
     character (len=12)     :: flabel
!**********************************************************************

!      step 1: redirecting data of h1 in ON MO basis (residing on h1_2cinf_saao on entry) to array scr_mat1
!      ------------------------------------------------------------------------------
       ioff_h1_inf_2cmoON = 1
       do i= 1, nfsym
         call dcopy(nesh_dim(i)**2 * nz,h1_2cinf_saao(ioff_h1_inf_2cmoON),1,scr_mat1(ioff_h1_inf_2cmoON),1)
         ioff_h1_inf_2cmoON = ioff_h1_inf_2cmoON + nesh_dim(i)**2 * nz
       end do

!      step 2: get transformation matrix V (constructed "on the fly")
!      ------------------------------------------------------------------------------
       call construct_saomo_overlap_half_tr(scr_mat2,             &
                                            h1_2cinf_saao,        &
                                            scr_mat3,             &
                                            nesh_dim,             &
                                            norb_dim,             &
                                            naosh_all,            &
                                            nz,                   &
                                            nzt,                  &
                                            nfsym,                &
                                            ipq_off_in,           &
                                            x2c_file_scr)

!      initialize pointers for matrices and vectors 
       ioff_tramat        = 1
       ioff_h1_inf_2cmoON = 1

!      step 3: obtain h1_infinite_2c_saao for each fermion ircop
!      ------------------------------------------------------------------------------
       do i = 1, nfsym
        
!        set dimensions for matrices used in the subroutine
         n1esh_f  = nesh_dim(i)
         nbas1_f  = naosh_all(i)

         if(nesh_dim(i) > 0)then

!          h1_2cinf_saao  =    V+     * h1_2cinf_moON *    V
!          h1_2cinf_saao  = scr_mat2+ * scr_mat1      * scr_mat2     
           call qtrans90('MOAO','S',val_d0,                &
                         nbas1_f,                          &
                         nbas1_f,                          &
                         n1esh_f,                          &
                         n1esh_f,                          &
                         h1_2cinf_saao(1+ioff_aomx(i,i)),  &
                         naosh_ls,                         &
                         naosh_ls,                         &
                         nz,                               &
                         ipq_off_in,                       &
                         scr_mat1(ioff_h1_inf_2cmoON),     &
                         n1esh_f,                          &
                         n1esh_f,                          &
                         nz,                               &
                         ipq_off_in,                       &
                         scr_mat2(ioff_tramat),            &
                         nbas1_f,                          &
                         n1esh_f,                          &
                         nz,                               &
                         ipq_off_in,                       &
                         scr_mat2(ioff_tramat),            &
                         nbas1_f,                          &
                         n1esh_f,                          &
                         nz,                               &
                         ipq_off_in,                       &
                         print_lvl)

!          debug print
           if(print_lvl > 1)then
             call print_x2cmat(h1_2cinf_saao(1+ioff_aomx(i,i)),nbas1_f,naosh_ls,nz,ipq_off_in(1,0),          &
                               'x2c - h1_inf_2c-AO',6,.true.)
           end if

         end if

!        update offsets for matrices
         ioff_h1_inf_2cmoON = ioff_h1_inf_2cmoON + n1esh_f * n1esh_f * nz
         ioff_tramat        = ioff_tramat        + n1esh_f * nbas1_f * nz

       end do

!      step 4: store h1_inf_2c in SA-AO basis on file (tag i = 0: marks data storage of "gerade+ungerade" (if nfsym > 1))
!      ------------------------------------------------------------------------------

!      are we going to add atomic-mean-field corrections in a subsequent step?
       if(add_amfi_contrib > 0)then

!        set (intermediate) pointer to Hamiltonian integral type - final adaption after AMFI step
         is_final_ham_lvl = 5

         write(flabel,'(a7,i4,i1)') 'h12cAOa',1,0
!                                          ^
!                                          |
!                                          a ("amfi preparation"); label will be
!                                          changed to "A" after the addition of
!                                          AMFI contributions
         call x2c_write(flabel,h1_2cinf_saao,naosh_ls*naosh_ls*nz,x2c_file_amf)
       else
 
!        set pointer to Hamiltonian integral type 
         if(.not.fragment_x2c)then
           is_final_ham_lvl = 1
         else
           is_final_ham_lvl = 2
         end if

         if(basic_h1mat == 2)then ! basic h1-mat == 4c-fock operator, aka "mmf"

            ! We first save the full 4c Fock matrix (the "F" in "h12cAOF") transformed to
            ! X2C basis in 4c storage mode because we use it later to get the canonical 2c
            ! MOs before all AO and MO index vectors are changed to 2c (i.e. e.g. ntbas(2) = 0 etc.)
            ! Below (label "h12cAOn") we save the 2c storage mode part, i.e. the "++" part, for
            ! use as the X2c one-electron operator and mmf two-electron operator.

            write(flabel,'(a7,i4,i1)') 'h12cAOF',1,0
            call x2c_write(flabel,h1_2cinf_saao,naosh_ls**2*nz,x2c_file_glb)

!           debug print
            if(print_lvl > 1)then
              do i = 1, nfsym
                nbas1_f  = naosh_all(i)
                call print_x2cmat(h1_2cinf_saao(1+ioff_aomx(i,i)),nbas1_f,naosh_ls,nz,ipq_off_in(1,0),      &
                   'x2c - h1_2nd_2c-AO',6,.true.)
              end do
            end if

         end if ! basic_h1mat

         write(flabel,'(a7,i4,i1)') 'h12cAOn',1,0
!                                          ^
!                                          |
!                                          n ("no amfi") == final Hamiltonian integrals

!        select LL-block of h1_inf_2c
           call dzero(scr_mat2,naosh_ls**2*nz)
           call pick_LL_block_saao_bas(h1_2cinf_saao,       &
                                       scr_mat2,            &
                                       naosh_L,             &
                                       naosh_all,           &
                                       naosh_ls,            &
                                       nr_ao_L,             &
                                       nfsym,               &
                                       nz,                  &
                                       0,                   &
                                       op_bs_to_fs)

!          debug print
           if(print_lvl > 1)then
             call print_x2cmat(scr_mat2,nr_ao_L,nr_ao_L,nz,ipq_off_in(1,0),'x2c - h1_inf_2c-AO-LL',6)
           end if

           if(spherical_on)then
  
!            set pointer to WORK array and allow full length for allocations
             lwrk = is_len_wrk
             allocate(wrk(lwrk))
             call sph_iotc(scr_mat2,h1_2cinf_saao,nr_sph,print_lvl,wrk,lwrk)
!            reset pointer to WORK array
             deallocate(wrk)

             call x2c_write(flabel,h1_2cinf_saao,nr_sph**2  * nz,x2c_file_glb)
           else
             call x2c_write(flabel,scr_mat2     ,nr_ao_L**2 * nz,x2c_file_glb)
           end if

       end if ! amfi contributions 

  end subroutine get_h1_inf_2c_saao_bas

!**********************************************************************

end module x2c_h1inf2c_mo2ao
