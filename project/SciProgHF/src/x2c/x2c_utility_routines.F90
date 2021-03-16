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
module x2c_utils


  use x2c_fio
  use x2c_cb_interface, only:  &
      set_x2c_subblock_data

  implicit none

  public construct_pctmat_saao_bas
  public construct_saomo_overlap_half_tr
  public generic_interface_ao2mo_mo2ao
  public store_aomo_trafo_matrices
  public x2c_lowd2c
  public print_x2cmat
  public print_x2cmtx
  public read_1fock_x2c
  public read_4cfock_operator_x2c
  public read_fock_matrices_saao_basis_x2c
  public write_fock_matrices_saao_basis_x2c
  public get_Sao_LL_mat
  public get_saomo_x2c
  public get_boson_irrep_info
  public make_spinfree_h1_onmo
  public make_two_e_fock_so_exclusive
  public get_2e_fock_matrix_x2c
  public get_2e_2c_fock_matrix_x2c
  public hello_dirx2c
  public goodbye_dirx2c

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0

contains

!**********************************************************************
  subroutine construct_pctmat_saao_bas(xmat,               &
                                       scr_mat1,           &
                                       scr_mat2,           &
                                       scr_mat3,           &
                                       norb_dim,           &
                                       nesh_dim,           &
                                       naosh_all,          &
                                       naosh_L,            &
                                       nfsym,              &
                                       nz,                 &
                                       nzt,                &
                                       ipq_off_in,         &
                                       x2c_file_scr,       &
                                       x2c_file_glb,       &
                                       print_lvl)
!**********************************************************************
!
!    purpose: construct the picture-change transformation matrix "U"
!             in DIRAC-sorted SA-AO basis
!             using the AO-to-MO transformation matrix V and the
!             pct-matrix Umat_ON in orthonormal MO basis.
!
!             equation: U = RKB AOtoMO trafo mat ["SL resorted"]  * Umat_ON
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: xmat(*)
     real(8), intent(inout) :: scr_mat1(*)
     real(8), intent(inout) :: scr_mat2(*)
     real(8), intent(inout) :: scr_mat3(*)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: nesh_dim(nfsym)
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: naosh_L(nfsym)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: nzt
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: x2c_file_scr
     integer, intent(in)    :: x2c_file_glb
     integer, intent(in)    :: print_lvl
!----------------------------------------------------------------------
     integer                :: i
     integer                :: n1esh_f
     integer                :: norb1_f
     integer                :: nbas1_f
     integer                :: nbasL_f
     integer                :: ioff_vmat
     integer                :: j,k,l
     character (len=12)     :: flabel
     character (len=28)     :: debug_string
!**********************************************************************

!      initialize offset for matrix
       ioff_vmat = 1

!      step 1: get transformation matrix V (constructed "on the fly")
!      ------------------------------------------------------------------------------
!         out: V -> scr_mat1
       call construct_saomo_overlap_half_tr(scr_mat1,             &
                                            xmat,                 &
                                            scr_mat2,             &
                                            nesh_dim,             &
                                            norb_dim,             &
                                            naosh_all,            &
                                            nz,                   &
                                            nzt,                  &
                                            nfsym,                &
                                            ipq_off_in,           &
                                            x2c_file_scr)

!      step 2: obtain pct matrix in DIRAC-sorted SA-AO basis
!      ------------------------------------------------------------------------------
       do i = 1, nfsym

!        set dimensions for matrices used in the subroutine
         n1esh_f  = nesh_dim(i)
         norb1_f  = norb_dim(i)
         nbas1_f  = naosh_all(i)
         nbasL_f  = naosh_L(i)

         if(norb_dim(i) > 0)then

!          read picture-change transformation matrix in orthonormal MO basis from scratch file
           write(flabel,'(a7,i4,i1)') 'Umat_ON',1,i
           call x2c_read(flabel,scr_mat3,norb1_f*n1esh_f*nz,x2c_file_scr)

!          read RKB AOtoMO transformation matrix ["SL resorted"] from file
           write(flabel,'(a7,i4,i1)') '4cAOMOr',1,i
           call x2c_read(flabel,xmat,nbas1_f*norb1_f*nzt,x2c_file_scr)

!          pctmt (4c-SA-AO -> 2c-ON-RKB) = RKB AOtoMO trafo mat ["SL resorted"]  * Umat_ON
!          scr_mat2                      = xmat                                  * scr_mat3
           call qgemm(nbas1_f,                                  &
                      n1esh_f,                                  &
                      norb1_f,                                  &
                      val_d1,                                   &
                      'N','N',                                  &
                      ipq_off_in,                               &
                      xmat,                                     &
                      nbas1_f,                                  &
                      norb1_f,                                  &
                      nzt,                                      &
                      'N','N',                                  &
                      ipq_off_in,                               &
                      scr_mat3,                                 &
                      norb1_f,                                  &
                      n1esh_f,                                  &
                      nz,                                       &
                      val_d0,                                   &
                      ipq_off_in,                               &
                      scr_mat2,                                 &
                      nbas1_f,                                  &
                      n1esh_f,                                  &
                      nz)

!          pctmt 4c-SA-AO -> 2c-SA-AO    = pctmt (4c-SA-AO -> 2c-ON-RKB)         * V+
!          xmat                          = scr_mat2                              * scr_mat1
           call qgemm(nbas1_f,                                  &
                      nbasL_f,                                  &
                      n1esh_f,                                  &
                      val_d1,                                   &
                      'N','N',                                  &
                      ipq_off_in,                               &
                      scr_mat2,                                 &
                      nbas1_f,                                  &
                      n1esh_f,                                  &
                      nz,                                       &
                      'H','N',                                  &
                      ipq_off_in,                               &
                      scr_mat1(ioff_vmat),                      &
                      nbas1_f,                                  &
                      n1esh_f,                                  &
                      nz,                                       &
                      val_d0,                                   &
                      ipq_off_in,                               &
                      xmat,                                     &
                      nbas1_f,                                  &
                      nbasL_f,                                  &
                      nz)

!          store picture-change transformation matrix in DIRAC-sorted SA-AO basis on global file
           write(flabel,'(a7,i4,i1)') 'pctmtAO',1,i
           call x2c_write(flabel,xmat,nbas1_f*nbasL_f*nz,x2c_file_glb)

!          debug print
           if(print_lvl > 2)then
             write(debug_string,'(a27,i1)') 'x2c - pctmat_saao - fsym = ',i
             call print_x2cmat(xmat,nbas1_f,nbasL_f,nz,ipq_off_in,debug_string,6)
!mathematica write(debug_string,'(a27,i1)') 'x2c - mtxpctmat_saao - fsym = ',i
!mathematica call print_x2cmtx(xmat,nbas1_f,nbasL_f,nz,ipq_off_in,debug_string,6)
           end if

         end if

!        update offset for matrix
         ioff_vmat          = ioff_vmat + n1esh_f * nbas1_f * nz

       end do

  end subroutine construct_pctmat_saao_bas


!**********************************************************************
  subroutine construct_saomo_overlap_half_tr(s_aomo_half_t, &
                                             scr_mat1,      &
                                             scr_mat2,      &
                                             nesh_dim,      &
                                             norb_dim,      &
                                             naosh_all,     &
                                             nz,            &
                                             nzt,           &
                                             nfsym,         &
                                             ipq_off_in,    &
                                             x2c_file_scr)
!**********************************************************************
!
!    purpose:
!       construct half transformed (== right index transformed) overlap matrix:
!       Saomo_Le = Sao(LL) * C(L, electronic shells)
!         (formerly known as "TM2CAOFP")
!       returned in s_aomo_half_t
!
!    requires:  4cAOMOr matrix (present on file "x2cscr")
!----------------------------------------------------------------------
     real(8), intent(inout) :: s_aomo_half_t(*)
     real(8), intent(inout) :: scr_mat1(*)
     real(8), intent(inout) :: scr_mat2(*)
     integer, intent(in)    :: nesh_dim(nfsym)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: nz
     integer, intent(in)    :: nzt
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: x2c_file_scr
!----------------------------------------------------------------------
     integer                :: i
     integer                :: ioff_aomo_half_t
     integer                :: ioff_aomo_sl_r
     integer                :: ioff_temp_sl_r
     integer                :: ioff_ovl
     integer                :: n1esh_f
     integer                :: norb1_f
     integer                :: nbas1_f
     character (len=12)     :: flabel
!**********************************************************************

!      initialize offsets for matrices
       ioff_aomo_half_t = 1
       ioff_temp_sl_r   = 1
       ioff_ovl         = 1

       do  i = 1, nfsym

!        set dimensions for matrices used in the subroutine
         nbas1_f  = naosh_all(i)
         norb1_f  = norb_dim(i)
         n1esh_f  = nesh_dim(i)

         if(nesh_dim(i) > 0)then

!          read (task_id < 0) Sao_LL overlap matrix from file (always read, thus mdirac == .false.)
           call get_Sao_LL_mat(scr_mat1(ioff_ovl),.false.,nbas1_f,i,-1)

!          read RKB AOtoMO transformation matrix ["SL resorted"] from file
           write(flabel,'(a7,i4,i1)') '4cAOMOr',1,i
           call x2c_read(flabel,scr_mat2(ioff_temp_sl_r),nbas1_f*norb1_f*nzt,x2c_file_scr)

!          proper offset in qgemm for the "SL resorted" matrix
           ioff_aomo_sl_r = ioff_temp_sl_r + nbas1_f * n1esh_f

!          Saomo_Le          = Sao_LL    * RKB AOtoMO trafo mat ["SL resorted"]
!          s_aomo_half_t     = scr_mat1  * scr_mat2
           call qgemm(nbas1_f,                                  &
                      n1esh_f,                                  &
                      nbas1_f,                                  &
                      val_d1,                                   &
                      'N','N',                                  &
                      ipq_off_in,                               &
                      scr_mat1(ioff_ovl),                       &
                      nbas1_f,                                  &
                      nbas1_f,                                  &
                      1,                                        &
                      'N','N',                                  &
                      ipq_off_in,                               &
                      scr_mat2(ioff_aomo_sl_r),                 &
                      nbas1_f,                                  &
                      n1esh_f,                                  &
                      nzt,                                      &
                      val_d0,                                   &
                      ipq_off_in,                               &
                      s_aomo_half_t(ioff_aomo_half_t),          &
                      nbas1_f,                                  &
                      n1esh_f,                                  &
                      nz)

         end if

!        update offsets for matrices
         ioff_ovl         = ioff_ovl + nbas1_f**2
         ioff_temp_sl_r   = ioff_temp_sl_r   + nbas1_f*norb1_f*nzt
         ioff_aomo_half_t = ioff_aomo_half_t + nbas1_f*n1esh_f*nz

       end do

  end subroutine construct_saomo_overlap_half_tr
!**********************************************************************

  subroutine generic_interface_ao2mo_mo2ao(scrmat1,                &
                                           scrmat2,                &
                                           scrmat3,                &
                                           ndim_i,                 &
                                           ndim_ii,                &
                                           ndim_iii,               &
                                           nz_1,                   &
                                           nz_2,                   &
                                           nz_3,                   &
                                           ipq_off_in,             &
                                           transform_direction,    &
                                           symmetry_property,      &
                                           print_lvl)
!**********************************************************************
!
!    purpose: generic interface to qtrans90
!             scrmat2 = scrmat3+ * scrmat1 * scrmat3
!
!----------------------------------------------------------------------
     real(8), intent(inout)        :: scrmat1(*)
     real(8), intent(inout)        :: scrmat2(*)
     real(8), intent(inout)        :: scrmat3(*)
     integer, intent(in)           :: ndim_i
     integer, intent(in)           :: ndim_ii
     integer, intent(in)           :: ndim_iii
     integer, intent(in)           :: nz_1
     integer, intent(in)           :: nz_2
     integer, intent(in)           :: nz_3
     integer, intent(in)           :: ipq_off_in(4,0:7)
     integer, intent(in)           :: print_lvl
     character (len=4), intent(in) :: transform_direction
     character (len=1), intent(in) :: symmetry_property
!----------------------------------------------------------------------

!      h1_4cON = V+ * h1_4cao * V
       call qtrans90(transform_direction,        &
                     symmetry_property,          &
                     val_d0,                     &
                     ndim_i,                     &
                     ndim_i,                     &
                     ndim_ii,                    &
                     ndim_ii,                    &
                     scrmat1,                    &
                     ndim_iii,                   &
                     ndim_iii,                   &
                     nz_1,                       &
                     ipq_off_in(1,0),            &
                     scrmat2,                    &
                     ndim_ii,                    &
                     ndim_ii,                    &
                     nz_2,                       &
                     ipq_off_in(1,0),            &
                     scrmat3,                    &
                     ndim_i,                     &
                     ndim_ii,                    &
                     nz_3,                       &
                     ipq_off_in(1,0),            &
                     scrmat3,                    &
                     ndim_i,                     &
                     ndim_ii,                    &
                     nz_3,                       &
                     ipq_off_in(1,0),            &
                     print_lvl)

  end subroutine generic_interface_ao2mo_mo2ao
!**********************************************************************

  subroutine store_aomo_trafo_matrices(mat1,                &
                                       mat2,                &
                                       mat3,                &
                                       norb_dim,            &
                                       naosh_all,           &
                                       nfsym,               &
                                       nz,                  &
                                       nzt,                 &
                                       defining_h1mat,      &
                                       x2c_file_scr,        &
                                       ipq_off_in,          &
                                       print_lvl,           &
                                       linear_sym)
!**********************************************************************
!
!    purpose:
!       save RKB AO to MO transformation matrices to file
!         1. in initial implementation SL_TM_4C:
!            required in general     ["SL resorted" transformation matrix] {ID == "r"}
!         2. in initial implementation TMAT4C  :
!            if(4c Dirac-Fock--> 2c) ["original" transformation matrix]    {ID == "o"}
!         3. in initial implementation TMAT_LSY:
!            if(4c Dirac-Fock--> 2c .and. linear == .true.) [non-linear --> linear symmetry] {ID == "l"}
!
!       re-initialize orbital common blocks in dcbdhf.h if (4c Dirac-Fock--> 2c)
!
!----------------------------------------------------------------------
     real(8), intent(in)    :: mat1(*)
     real(8), intent(in)    :: mat2(*)
     real(8), intent(inout) :: mat3(*)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: nzt
     integer, intent(in)    :: defining_h1mat
     integer, intent(in)    :: x2c_file_scr
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: linear_sym
!----------------------------------------------------------------------
     integer                :: i
     integer                :: ioff_aomo
     integer                :: ioff_momo
     character (len=12)     :: flabel
!**********************************************************************

!      initialize offsets
       ioff_aomo = 1
       ioff_momo = 1

       do i = 1, nfsym

         write(flabel,'(a7,i4,i1)') '4cAOMOr',1,i
         call x2c_write(flabel,mat1(ioff_aomo),naosh_all(i)*norb_dim(i)*nzt,x2c_file_scr)

         if(defining_h1mat == 2)then
           write(flabel,'(a7,i4,i1)') '4cAOMOo',1,i
           call x2c_write(flabel,mat2(ioff_aomo),naosh_all(i)*norb_dim(i)*nzt,x2c_file_scr)

           if(linear_sym)then
             write(flabel,'(a7,i4,i1)') '4cMOMOl',1,i
             call x2c_write(flabel,mat3(ioff_momo),norb_dim(i)**2 * nz,x2c_file_scr)
           end if

         end if

         ioff_aomo = ioff_aomo + naosh_all(i) * norb_dim(i)*nzt
         ioff_momo = ioff_momo + norb_dim(i) *  norb_dim(i)*nz
       end do

       if(defining_h1mat == 2)then
        call set_tmt(0)
       end if

!      redirect ["SL resorted"] matrix to mat3 for later reuse in
!      x2c_get_h1_mat_base (module x2c_h1mat_base)
       ioff_aomo = 1
       do i = 1, nfsym
         mat3(ioff_aomo:naosh_all(i)*norb_dim(i)*nz ) = 0.0d0
         call dcopy(naosh_all(i)*norb_dim(i)*nzt,mat1(ioff_aomo),1,mat3(ioff_aomo),1)
!        debug print
         if(print_lvl > 2)then
           call print_x2cmat(mat3(ioff_aomo),naosh_all(i),norb_dim(i),nzt,ipq_off_in,'x2c - AO2MO resorted',6)
         end if
         ioff_aomo = ioff_aomo + naosh_all(i) * norb_dim(i)*nzt
       end do

  end subroutine store_aomo_trafo_matrices

!**********************************************************************
  subroutine x2c_lowd2c(lowd_mat,            &
                        scr_mat1,            &
                        scr_mat2,            &
                        npsh_dim,            &
                        nesh_dim,            &
                        naosh_L,             &
                        naosh_S,             &
                        naosh_all,           &
                        norb_dim_lwd_L,      &
                        norb_dim_lwd_S,      &
                        norb_dim_lwd_all,    &
                        naosh_ls,            &
                        nfsym,               &
                        nz,                  &
                        ipq_off_in,          &
                        x2c_file_glb,        &
                        print_lvl,           &
                        linear_sym)
!**********************************************************************
!
!    purpose:
!         construct the 2c-Lowdin matrix (only needed in the x2c-mmf approach).
!         in case of linear symmetry we transform the 2c-Lowdin matrix
!         by means of the <j_z>_LL eigenvectors.
!
!         out:                          regular 2c-LOWDIN matrix on file X2CMAT
!              if linear_sym == .true.: linear  2c-LOWDIN matrix on file X2CMAT
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: lowd_mat(*)
     real(8), intent(inout) :: scr_mat1(*)
     real(8), intent(inout) :: scr_mat2(*)
     integer, intent(in)    :: npsh_dim(nfsym)
     integer, intent(in)    :: nesh_dim(nfsym)
     integer, intent(in)    :: naosh_L(nfsym)
     integer, intent(in)    :: naosh_S(nfsym)
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: norb_dim_lwd_L(nfsym)
     integer, intent(in)    :: norb_dim_lwd_S(nfsym)
     integer, intent(in)    :: norb_dim_lwd_all(nfsym)
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: x2c_file_glb
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: linear_sym
!----------------------------------------------------------------------
     integer                :: i
     integer                :: j
     integer                :: nbasL_f
     integer                :: nbasS_f
     integer                :: nbas1_f
     integer                :: ilowdmat
     integer                :: ilowdmat_tmp
     integer                :: nlowdmat
     integer                :: nborb_dim(4,2,0:2)
     character (len=12)     :: flabel
!**********************************************************************

!      initialize offset
       ilowdmat = 1
       nlowdmat = 0

!      step 1: zero out L-S terms
!      --------------------------
       do i = 1, nfsym

!        initialize
         ilowdmat_tmp = ilowdmat
         nbasL_f      = naosh_L(i)
         nbasS_f      = naosh_S(i)
         nbas1_f      = naosh_all(i)

!        upper triangle
         do j = 1, norb_dim_lwd_l(i)
           call dzero(lowd_mat(ilowdmat_tmp+nbasL_f),nbasS_f)
           ilowdmat_tmp = ilowdmat_tmp + nbas1_f
         end do

!        lower triangle
         ilowdmat_tmp   = ilowdmat     + nbas1_f * norb_dim_lwd_L(i)
         call dzero(lowd_mat(ilowdmat_tmp),norb_dim_lwd_S(i)*nbas1_f)

!        update offset in LOWDIN matrix
         ilowdmat = ilowdmat + nbas1_f * norb_dim_lwd_all(i)
       end do

!      step 2: reorder the L blocks for nfsym > 1:
!              L for nfsym==1 followed by L for nfsym==2
!      -------------------------------------------------
       nlowdmat = naosh_all(1) * nesh_dim(1)
       if(nfsym == 2)then
         ilowdmat_tmp =        1 + naosh_all(1)     * norb_dim_lwd_all(1)
         ilowdmat     =        1 + naosh_all(1)     * nesh_dim(1)
         nlowdmat     = nlowdmat + naosh_all(nfsym) * nesh_dim(nfsym)
         call dcopy(naosh_all(2)*nesh_dim(2),lowd_mat(ilowdmat_tmp),1,lowd_mat(ilowdmat),1)
       end if

!      step 3: store 2c-adapted Lowdin matrix in AO2MO basis on global file (tag i = 0: marks data storage of "g+u" (if nfsym > 1))
!                    ^^         ^ ^              ^^ ("g" == general case)
!      label code:   ||         | |              ||   |
       write(flabel,'(a7,i4,i1)') '2cLWAOg',1,0


       call x2c_write(flabel,lowd_mat,nlowdmat,x2c_file_glb)


!      step 4: transform 2c-adapted Lowdin matrix in AO2MO basis to account for linear symmetry
!      ----------------------------------------------------------------------------------------
       if(linear_sym)then

!        a. read the linear-symmetry adapted MO2MO transformation matrix
!        -----------------------------------------------------------------
         open(8,file='JzMOLLb',status='old',form='unformatted', &
              access='sequential',action="readwrite",position='rewind')

!        MO2MO linear symmetry <--> nonlinear symmetry transformation matrix
         ilowdmat = 1
         do i = 1, nfsym

           if(nesh_dim(i) > 0) call readt(8,nesh_dim(i)**2,scr_mat1(ilowdmat))

!          debug print
           if(print_lvl > 2)then
             call print_x2cmat(scr_mat1(ilowdmat),nesh_dim(i),nesh_dim(i),1,ipq_off_in(1,0),'trafo-lin',6)
           end if

           ilowdmat = ilowdmat + nesh_dim(i)**2
         end do

         close(8,status='delete')

!        transformation step: non-linear-symmetry LOWDIN --> linear-symmetry LOWDIN
         ilowdmat     = 1
         ilowdmat_tmp = 1
         do i = 1, nfsym
           call bcktra(scr_mat2(ilowdmat),naosh_all(i),nesh_dim(i),     &
                       scr_mat1(ilowdmat_tmp),nesh_dim(i),nesh_dim(i),  &
                       nesh_dim(i),1,                                   &
                       nesh_dim(i),1,naosh_all(i),                      &
                       lowd_mat(ilowdmat),naosh_all(i),nesh_dim(i),1,   &
                       print_lvl)

!          debug print
           if(print_lvl > 2)then
             call print_x2cmat(scr_mat2(ilowdmat),naosh_all(i),nesh_dim(i),1,ipq_off_in(1,0),'lowd-lin',6)
           end if

           ilowdmat     = ilowdmat     + naosh_all(i)*nesh_dim(i)
           ilowdmat_tmp = ilowdmat_tmp + nesh_dim(i)**2
         end do

!        store adapted 2c-LOWDIN matrix on file
         write(flabel,'(a7,i4,i1)') '2cLWAOl',1,0
         call x2c_write(flabel,scr_mat2,nlowdmat,x2c_file_glb)

       end if


!      we always want to call the common block manipulation routines as the sub-block
!      data (taken from dcborb.h) is written in there
       call set_x2c_subblock_data(nfsym,                &
                                  nz)

  end subroutine x2c_lowd2c

!**********************************************************************
  subroutine print_x2cmat(xmat,        &
                          ndim_x,      &
                          ndim_y,      &
                          nz,          &
                          ipq_off,     &
                          title,       &
                          prntunit,    &
                          lrlc_equal)
!**********************************************************************
!
!    purpose: print a matrix
!
!**********************************************************************
     real(8), intent(in)           :: xmat(*)
     integer, intent(in)           :: ndim_x
     integer, intent(in)           :: ndim_y
     integer, intent(in)           :: nz
     integer, intent(in)           :: ipq_off(4,*)
     integer, intent(in)           :: prntunit
     logical, optional, intent(in) :: lrlc_equal
     character(*), intent(in)      :: title
!**********************************************************************
!
         call header2(title,-1,prntunit)
         if(.not.present(lrlc_equal))then
           call prqmat(xmat,ndim_x,ndim_y,ndim_x,ndim_y,nz,  &
                       ipq_off,prntunit)
         else
           call prqmat(xmat,ndim_x,ndim_x,ndim_y,ndim_y,nz,  &
                       ipq_off,prntunit)
         end if

  end subroutine print_x2cmat

!***********************************************************************
  subroutine print_x2cmtx(xmat,        &
                          ndim_x,      &
                          ndim_y,      &
                          nz,          &
                          ipq_off,     &
                          title,       &
                          prntunit)
!**********************************************************************
!
!    purpose: print a matrix
!
!**********************************************************************
     real(8), intent(in)      :: xmat(ndim_x,ndim_y)
     integer, intent(in)      :: ndim_x
     integer, intent(in)      :: ndim_y
     integer, intent(in)      :: nz
     integer, intent(in)      :: ipq_off(4,0:7)
     integer, intent(in)      :: prntunit
     character(*), intent(in) :: title
!-----------------------------------------------------------------------
     integer                  :: irow, icol
!**********************************************************************
!
         call header(title,-1)
         do icol= 1,ndim_y
           do irow=1,ndim_x
             write(prntunit,*) irow, icol, xmat(irow,icol)
           enddo
         enddo

  end subroutine print_x2cmtx

!**********************************************************************
  subroutine read_1fock_x2c(xmat,      &
                            ndim_tot,  &
                            nz,        &
                            ham_lvl)
!**********************************************************************
!
!    purpose: read 2c-Hamiltonian integrals from file X2CMAT calculated
!             within the X2C module.
!             The parameter "ham_lvl" points to the integrals type,
!             i.e. what kind of exact-two-component Hamiltonian will be read
!             from file.
!
!    ham_lvl  == 1: oo-order 2c-Hamiltonian without AMFI contributions
!    ham_lvl  ==-2: mmf-X2C
!   |ham_lvl| == 3: oo-order 2c-Hamiltonian with AMFI (old == -3 or new == +3) contributions
!
!**********************************************************************
     real(8), intent(inout)   :: xmat(*)
     integer, intent(in)      :: ndim_tot
     integer, intent(in)      :: nz
     integer, intent(in)      :: ham_lvl
!----------------------------------------------------------------------
     integer                  :: abs_ham_lvl
     integer                  :: i
     integer                  :: j
     integer                  :: file_unit
     logical                  :: isopen
     character (len=12)       :: flabel
!**********************************************************************

       if(ham_lvl == -2) return
       abs_ham_lvl = abs(ham_lvl)
       if(abs_ham_lvl <= 0 .or. abs_ham_lvl > 4)           &
         call quit('*** error in read_1fock_x2c: invalid pointer to Hamiltonian integral type. ***')

!      step 1: determine file label according to integral type flag
!      ------------------------------------------------------------

       i = 0
       j = 0
       if(abs_ham_lvl <= 2)then
         write(flabel,'(a7,i4,i1)') 'h12cAOn',1,i
       else
         write(flabel,'(a7,i4,i1)') 'h12cAOA',1,i
       end if

!      step 2: read selected integrals from file X2CMAT
!      ------------------------------------------------

       file_unit = -1
       inquire(file='X2CMAT',opened=isopen,number=file_unit)

       if(.not.isopen)then
!        open file
         file_unit = 99
         open(file_unit,file='X2CMAT',status='old',form='unformatted',         &
              access='sequential',action='read',position='rewind')
       end if

       call x2c_read(flabel,xmat,ndim_tot * nz,file_unit)

       close(file_unit,status='keep')

  end subroutine read_1fock_x2c

!**********************************************************************
  subroutine read_4cfock_operator_x2c(xmat,           &
                                      op_mat_on_file, &
                                      print_lvl)
!**********************************************************************
!
!    purpose: read the 4c-Fock operator matrix elements (in ON-MO basis)
!             from file op_mat_on_file within the X2C module
!             if the basic and/or the defining h1 matrix
!             is the 4c-Fock operator (controlled by defining_h1mat == 2).
!
!**********************************************************************
     use x2cmod_cfg
     real(8),           intent(inout) :: xmat(*)
     character (len=6), intent(in)    :: op_mat_on_file
     integer,           intent(in)    :: print_lvl
#include "dcbdhf.h"
!----------------------------------------------------------------------
     integer                          :: real_bytes    = 8
     integer                          :: i
     integer                          :: record_number
     integer                          :: file_unit
     logical                          :: does_exist
     real(8), allocatable             :: temp_mat(:)
     character (len=28)               :: debug_string
!**********************************************************************

!      THIS IS ONLY A DEBUG OPTION!!!!!!!!!!!!!!!!!
       if(x2c_fock_saao_basis)then

         allocate(temp_mat(nr_ao_total_x2c**2*nr_quat))
         temp_mat = 0
         file_unit = 99

         open(file_unit,file='DFFCK1',status='old',form='unformatted',access='sequential',action='read')
         call read_fock_matrices_saao_basis_x2c(xmat,                       &
                                                nr_ao_total_x2c**2*nr_quat, &
                                                1,                          &
                                                file_unit,                  &
                                                print_lvl)
         close(file_unit, status='keep')
!        debug print
         if(print_lvl > 2)then
           write(debug_string,'(a27,i1)') 'x2c - 4c-fock_AO 1e- fsym = ',0
           call print_x2cmat(xmat,nr_ao_total_x2c,nr_ao_total_x2c,nr_quat,x2c_cb_pq_to_uq,debug_string,6)
         end if

         open(file_unit,file='DFFCK2',status='old',form='unformatted',access='sequential',action='read')
         call read_fock_matrices_saao_basis_x2c(temp_mat,                   &
                                                nr_ao_total_x2c**2*nr_quat, &
                                                1,                          &
                                                file_unit,                  &
                                                print_lvl)
         close(file_unit, status='keep')
!        debug print
         if(print_lvl > 2)then
           write(debug_string,'(a27,i1)') 'x2c - 4c-fock_AO 2e- fsym = ',0
           call print_x2cmat(temp_mat,nr_ao_total_x2c,nr_ao_total_x2c,nr_quat,x2c_cb_pq_to_uq,debug_string,6)
         end if


         call daxpy(nr_ao_total_x2c**2*nr_quat,val_d1,temp_mat,1,xmat,1)
!        debug print
         if(print_lvl > 2)then
           write(debug_string,'(a27,i1)') 'x2c - 4c-fock_AO 12- fsym = ',0
           call print_x2cmat(xmat,nr_ao_total_x2c,nr_ao_total_x2c,nr_quat,x2c_cb_pq_to_uq,debug_string,6)
         end if
         deallocate(temp_mat)

       else ! end of DEBUG ONLY OPTION!

         !> read and dump the F1 Fock matrix in AO basis to file (if needed later)
         file_unit = 99
         open(file_unit,file='DFFCK1',status='old',form='unformatted',access='sequential',action='read')
         call read_fock_matrices_saao_basis_x2c(xmat,                       &
                                                nr_ao_total_x2c**2*nr_quat, &
                                                1,                          &
                                                file_unit,                  &
                                                print_lvl)

         open(81,file='4cf1AO',status='replace',form='unformatted',         &
               access='sequential',action='write',position='rewind')
         write(81) (xmat(i),i=1,nr_ao_total_x2c**2*nr_quat)
         close(81,status="keep")

         close(file_unit, status='keep')
         call dzero(xmat,nr_ao_total_x2c**2*nr_quat)

!      step 1: inquire file and set record number
!      ------------------------------------------
       inquire(file=op_mat_on_file,exist=does_exist)

       if(.not.does_exist) &
         call quit('*** error in read_4cfock_operator_x2c: file with 4c-fock matrix does not exist.***')

!      open file
       file_unit = 99
       open(file_unit,file=op_mat_on_file,status='old',form='unformatted',         &
            access='direct',action='read',recl=real_bytes*n2tmotq)

       record_number = 1
       if(scf_iter_counter > 0)then
         record_number = mod(scf_iter_counter-2,diis_counter-1) + 1
         if(record_number == 0)  &
         record_number = 1
       end if

!      step 2: read matrix elements from file
!      --------------------------------------
       call readac(file_unit,n2tmotq,xmat,record_number)

       close(file_unit,status='keep')
       end if

  end subroutine read_4cfock_operator_x2c

!**********************************************************************
  subroutine read_fock_matrices_saao_basis_x2c(xmat,           &
                                               ndim,           &
                                               nr_entries,     &
                                               fh,             &
                                               print_lvl,      &
                                               loop_factor     &
                                              )
!**********************************************************************
!
!    purpose: read the 4c-Fock operator matrix elements (in SA-AO basis)
!             from file fh within the X2C module
!             if the basic and/or the defining h1 matrix
!             is the 4c-Fock operator (controlled by defining_h1mat == 2).
!
!**********************************************************************
     integer,           intent(in)    :: ndim
     integer,           intent(in)    :: nr_entries
     integer,           intent(in)    :: fh
     integer,           intent(in)    :: print_lvl
     real(8), optional                :: loop_factor(0:nr_entries-1)
     real(8),           intent(inout) :: xmat(ndim * nr_entries)
!----------------------------------------------------------------------
     integer                          :: ndim_check, offset, i,j
     integer                          :: nr_entries_check
     character (len=74)               :: info_text
     character (len=72)               :: file_name
     logical                          :: squash
!**********************************************************************
       squash = .false.
       if(present(loop_factor))then
        squash = .true.
       end if

!      step 1: inquire file name
!      -------------------------
       inquire(unit=fh,name=file_name)

!      step 2: read matrix elements from file
!      --------------------------------------
       rewind fh
       read(fh) info_text
       read(fh) ndim_check, nr_entries_check
       if(print_lvl > 2)then
         print '(/2a)', ' info: fock matrix read from file ',file_name
         print '(/2a)', ' file header :',info_text
       end if
       if(ndim_check /= ndim .or. nr_entries_check /= nr_entries)then
         print *, ' mismatch of fock matrix dimensions: dimension on file = ',ndim_check,' input dimension = ',ndim
         print *, ' mismatch of fock matrix dimensions: # entries on file = ',nr_entries_check,' input # entries = ',&
                                                                              nr_entries
         call quit('*** error in read_fock_matrices_saao_basis_x2c: mismatch of fock matrix dimensions.')
       end if

       read(fh) xmat

       !> possibly squash AO matrices into one with scaling factors applied to each matrix, e.g. AOC factors
       if(nr_entries > 1 .and. squash)then
         offset = 1
         do i = 1, nr_entries
           if(i > 1)then
             call daxpy(ndim,loop_factor(i-1),xmat(offset+(i-1)*(ndim)),1,xmat,1)
           end if
           offset = offset + ndim
         end do
       end if

  end subroutine read_fock_matrices_saao_basis_x2c

!**********************************************************************
  subroutine write_fock_matrices_saao_basis_x2c(xmat,           &
                                                ndim,           &
                                                nr_entries,     &
                                                name,           &
                                                fh,             &
                                                print_lvl)
!**********************************************************************
!
!    purpose: write the 2c-Fock operator matrix elements (in SA-AO basis)
!             to file fh within the X2C module
!             if the basic and/or the defining h1 matrix
!             is the 4c-Fock operator (controlled by defining_h1mat == 2).
!
!**********************************************************************
     integer,           intent(in)    :: ndim
     integer,           intent(in)    :: nr_entries
     integer,           intent(in)    :: fh
     integer,           intent(in)    :: print_lvl
     real(8),           intent(in)    :: xmat(ndim * nr_entries)
!----------------------------------------------------------------------
     character (len=74)               :: info_text
     character (len=72)               :: file_name
     character (len=4)                :: name
!**********************************************************************

!      step 1: inquire file name
!      -------------------------
       inquire(unit=fh,name=file_name)

!      step 2: write matrix elements to file
!      -------------------------------------
       rewind fh
       write(info_text,'(a74)') 'X2C molecular-mean-field: picture-change transformed effective '//name//' matrix.'
       write(fh) info_text
       write(fh) ndim, nr_entries
       if(print_lvl > 2)then
         print '(/2a)', ' info: fock matrix written to file ',file_name
         print '(/2a)', ' file header :',info_text
       end if

       write(fh) xmat

  end subroutine write_fock_matrices_saao_basis_x2c

!**********************************************************************
  subroutine get_Sao_LL_mat(scr_mat1,    &
                            mdirac,      &
                            ndim_mat,    &
                            i,           &
                            task_id)
!**********************************************************************
!
!    purpose: compute, write and read the Sao_LL overlap matrix.
!
!             task_id >= 0: compute and write the Sao_LL overlap matrix
!             task_id <  0: read the Sao_LL overlap matrix (i-th fermion corep part)
!
!**********************************************************************
     real(8), intent(inout)   :: scr_mat1(*)
     integer, intent(in)      :: ndim_mat
     integer, intent(in)      :: i
     integer, intent(in)      :: task_id
     logical, intent(in)      :: mdirac
!----------------------------------------------------------------------
     integer                  :: j
     integer                  :: ndim_mat_tmp
     integer                  :: ioff_ovl
     integer                  :: file_unit
     character (len=12)       :: flabel
!**********************************************************************

       inquire(file='X2C_Sao_LL_scr', number=file_unit)

       ndim_mat_tmp = ndim_mat

       if(task_id >= 0)then

#ifdef PRG_DIRAC
!        compute Sao_LL overlap matrix where the SS part is set to zero
         call gtovlt(scr_mat1,0.0d0,0)
#endif

!        store matrix on file in "gerade/ungerade" format (if applicable)
         ioff_ovl = 1
         do j = 1, i

!          "hidden" code (dirty but helpful)
           if(j == 2) ndim_mat_tmp = task_id

           write(flabel,'(a11,i1)') 'Sao_LL_ovlp',j
           call x2c_write(flabel,scr_mat1(ioff_ovl),ndim_mat_tmp**2,file_unit)

           ioff_ovl = ioff_ovl + ndim_mat_tmp**2
         end do

       else

         write(flabel,'(a11,i1)') 'Sao_LL_ovlp',i
         call x2c_read(flabel,scr_mat1,ndim_mat_tmp**2,file_unit)
       end if

  end subroutine get_Sao_LL_mat
!**********************************************************************

  subroutine get_saomo_x2c(                &
                            saomo,         &
                            vaomo,         &
                            naosh_ls,      &
                            naosh_all,     &
                            norb_dim,      &
                            ioff_aomat_x,  &
                            nfsym,         &
                            nz,            &
                            nzt,           &
                            ipq_off_in,    &
                            mdirac,        &
                            print_lvl)
!**********************************************************************
!
!    purpose: compute the right-index transformed overlap matrix
!
!**********************************************************************
     real(8), intent(inout) :: saomo(*)
     real(8), intent(inout) :: vaomo(*)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: ioff_aomat_x(nfsym,nfsym)
     integer, intent(in)    :: nz
     integer, intent(in)    :: nzt
     integer, intent(in)    :: ipq_off_in(4,0:7)
     logical, intent(in)    :: mdirac
     integer, intent(in)    :: print_lvl
!----------------------------------------------------------------------
     real(8), allocatable   :: saoao(:)
     integer                :: iaomo
     integer                :: i
     integer                :: norb1_f
     integer                :: nbas1_f
!**********************************************************************

       allocate(saoao(nzt*naosh_ls**2))
       saoao(1:nzt*naosh_ls**2) = 0.0d0

!      compute Saoao overlap matrix
#ifdef PRG_DIRAC
       call gtovlx(saoao,1.0d0)
#endif

       iaomo = 1

       do i = 1, nfsym

!        set dimensions for matrices used in the subroutine
         norb1_f  = norb_dim(i)
         nbas1_f  = naosh_all(i)

!        debug print
         if(print_lvl > 2)then
           call print_x2cmat(vaomo(iaomo),nbas1_f,norb1_f,nzt,ipq_off_in,'x2c - v    ',6)
           call print_x2cmat(saoao(ioff_aomat_x(i,i)+1),nbas1_f,naosh_ls,nzt,ipq_off_in,'x2c - saoao',6)
         end if

         if(norb_dim(i) > 0)then

!          Saomo                         = Saoao * RKB AOtoMO trafo mat
!          saomo                         = saoao * vaomo
           call qgemm(nbas1_f,                                  &
                      norb1_f,                                  &
                      nbas1_f,                                  &
                      val_d1,                                   &
                      'N','N',                                  &
                      ipq_off_in,                               &
                      saoao(ioff_aomat_x(i,i)+1),               &
                      naosh_ls,                                 &
                      naosh_ls,                                 &
                      1,                                        &
                      'N','N',                                  &
                      ipq_off_in,                               &
                      vaomo(iaomo),                             &
                      nbas1_f,                                  &
                      norb1_f,                                  &
                      nzt,                                      &
                      val_d0,                                   &
                      ipq_off_in,                               &
                      saomo(iaomo),                             &
                      nbas1_f,                                  &
                      norb1_f,                                  &
                      nzt)

         end if

!        debug print
         if(print_lvl > 2)then
           call print_x2cmat(saomo(iaomo),nbas1_f,norb1_f,nzt,ipq_off_in,'x2c - saomo',6)
         end if

!        keep track of offset
         iaomo = iaomo + nzt * nbas1_f * norb1_f

       end do

       deallocate(saoao)

  end subroutine get_saomo_x2c
!**********************************************************************

  subroutine get_boson_irrep_info(boson_info,  &
                                  nfsym,       &
                                  norb_dim,    &
                                  print_lvl)
!**********************************************************************
!
!    purpose: obtain boson irrep symmetry info.
!
!**********************************************************************
     integer, intent(inout)   :: boson_info(*)
     integer, intent(in)      :: nfsym
     integer, intent(in)      :: norb_dim(nfsym)
     integer, intent(in)      :: print_lvl
!----------------------------------------------------------------------
     integer                  :: i
     integer                  :: ioff_orb
!**********************************************************************

!     initialize orbital subblock arrays norb_sub(ib,ifrp,0) and id_sub_bl(ib,ifrp)
      call inisub

      ioff_orb = 1
      do i = 1, nfsym
        call inibos(boson_info(ioff_orb),i,.true.,print_lvl)
        ioff_orb = ioff_orb + norb_dim(i)
      end do

  end subroutine get_boson_irrep_info
!**********************************************************************

  subroutine make_spinfree_h1_onmo(                            &
                                   h1_onmo_bas,                &
                                   lrow,                       &
                                   lcol,                       &
                                   nrow,                       &
                                   ncol,                       &
                                   boson_info_r,               &
                                   boson_info_c,               &
                                   nz,                         &
                                   print_lvl                   &
                                  )
!**********************************************************************
!
!    purpose: retain only spinfree elements in 1-el Hamiltonian
!
!**********************************************************************
     real(8), intent(inout)   :: h1_onmo_bas(lrow,lcol,nz)
     integer, intent(in)      :: boson_info_r(*)
     integer, intent(in)      :: boson_info_c(*)
     integer, intent(in)      :: lrow
     integer, intent(in)      :: lcol
     integer, intent(in)      :: nrow
     integer, intent(in)      :: ncol
     integer, intent(in)      :: nz
     integer, intent(in)      :: print_lvl
!----------------------------------------------------------------------
     integer                  :: i, j, k
!**********************************************************************

!      part 1: nz = 1
       do i = 1, nrow
         do j = 1, ncol
          if(boson_info_r(i) /= boson_info_c(j)) h1_onmo_bas(i,j,1) = 0.0d0
         end do
       end do

!      part 2: nz >= 1
       if(nz > 1)then
         do k = 2, nz
           do i = 1, nrow
             do j = 1, ncol
               h1_onmo_bas(i,j,k) = 0.0d0
             end do
           end do
         end do
       end if

  end subroutine make_spinfree_h1_onmo
!**********************************************************************

  subroutine make_two_e_fock_so_exclusive(                            &
                                          Gmat,                       &
                                          nrow,                       &
                                          ncol,                       &
                                          nz                          &
                                         )
!**********************************************************************
!
!    purpose: retain only two-electron spin-orbit elements in 2-el Fock matrix
!
!**********************************************************************
     real(8), intent(inout)   :: Gmat(nrow,ncol,4)
     integer, intent(in)      :: nrow
     integer, intent(in)      :: ncol
     integer, intent(in)      :: nz
!----------------------------------------------------------------------
     integer                  :: i, j
!**********************************************************************

!      check dimensionalities
       if(nz /= 4) stop '2-el Fock matrix is not in quaternion format'

!      take out scalar parts
       do j=1,ncol
         do i=1,nrow
           !> quaternion(0)
           Gmat(i     ,j     ,1) = 0.0d0                 ! ll
           Gmat(i     ,j     ,2) = 0.0d0                 ! ll
           Gmat(i     ,j     ,3) = 0.0d0                 ! ll
           Gmat(i     ,j     ,4) = 0.0d0                 ! ll
         enddo
       enddo

  end subroutine make_two_e_fock_so_exclusive
!**********************************************************************

  subroutine get_2e_fock_matrix_x2c(                            &
                                    Gmat,                       &
                                    nrows,                      &
                                    ncols,                      &
                                    nz,                         &
                                    nr_2e_fock_matrices,        &
                                    intflg,                     &
                                    is_len_wrk,                 &
                                    print_lvl,                  &
                                    aoc_factors                 &
                                   )
!**********************************************************************
!
!    purpose: get 2e-fock matrix with integral contributions according to
!             intflg (bit packed integer):
!             (LL|LL)                     == 1
!                       (SS|LL)           == 2
!             (LL|LL) + (SS|LL)           == 3
!                                 (SS|SS) == 4
!             (LL|LL)           + (SS|SS) == 5
!                       (SS|LL) + (SS|SS) == 6
!             (LL|LL) + (SS|LL) + (SS|SS) == 7
!
!             with GAUNT (LS|LS) entering the game (intflg=8) you know what to do by now...
!
!**********************************************************************
#ifdef VAR_MPI
#include "mpif.h"
#endif
     real(8), intent(inout)   :: Gmat(nrows,ncols,nz,nr_2e_fock_matrices)
     integer, intent(in)      :: nrows
     integer, intent(in)      :: ncols
     integer, intent(in)      :: nz
     integer, intent(in)      :: nr_2e_fock_matrices
     integer, intent(in)      :: intflg
     integer, intent(in)      :: is_len_wrk
     integer, intent(in)      :: print_lvl
     real(8), intent(in)      :: aoc_factors(*)
!----------------------------------------------------------------------
     integer                  :: i, npos, local_lwork, lucoef
     real(8)                  :: xx, j
     real(8), allocatable     :: work(:)
     real(8), allocatable     :: dmat(:)
     real(8), allocatable     :: cmo(:)
     integer, allocatable     :: isymop(:), ifckop(:), ihrmop(:), pos(:)
     logical, allocatable     :: saveflags(:)
     logical                  :: parcal
!**********************************************************************

 !   generate density matri(x/ces)
     allocate(dmat(nrows*ncols*nz*nr_2e_fock_matrices))
     allocate(cmo(nrows*ncols*nz))
     dmat = 0.0d0
     cmo  = 0.0d0
!    read coefficients
     lucoef = 99
     call opnfil(lucoef,'DFCOEF','OLD','x2cfck')
     rewind(lucoef)
     call reacmo(lucoef,'DFCOEF',cmo,xx,i,j,2)
     close(lucoef,status='keep')

!    generate density matri(x/ces) (denmat is a general driver routine from Dirac
!    computing all inactive and active density matrices)
     call denmat(dmat,cmo,print_lvl)
     deallocate(cmo)

     allocate(isymop(nr_2e_fock_matrices))
     allocate(ifckop(nr_2e_fock_matrices))
     allocate(ihrmop(nr_2e_fock_matrices))

     ihrmop = 0
     ifckop = 0
     isymop = 0

     do i = 1,nr_2e_fock_matrices
!      totally symmetric operator
       isymop(i) = 1
!      fock matrix type
       ifckop(i) = 1
!      hermitian operator
       ihrmop(i) = 1
     end do

     allocate(saveflags(4))

     call SaveTaskDistribFlags(saveflags)
     call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. ,.TRUE. /))
     parcal = .false.
#ifdef VAR_MPI
     call mpi_initialized(parcal,i)
#endif
     call SetIntTaskArrayDimension(npos,parcal)
     allocate(pos(npos))

     local_lwork = is_len_wrk - 3*nrows*ncols*nz*nr_2e_fock_matrices
     allocate(work(local_lwork))

!    2e-fock matrices (Dirac driver routine)
     call twofck(isymop,ihrmop,ifckop,Gmat,dmat,nr_2e_fock_matrices,   &
                 pos,intflg,print_lvl,work,local_lwork)

!     accumulate active fock matrices (add to inactive fock matrix)
      if(nr_2e_fock_matrices > 1)then
        do i = 1, nr_2e_fock_matrices-1
          call daxpy(nrows*ncols*nz,aoc_factors(i),Gmat(1,1,1,i+1),1,Gmat(1,1,1,1),1)
        end do
     end if

     deallocate(work)
     if(parcal) call SetTaskDistribFlags(saveflags)
     deallocate(pos)
     deallocate(saveflags)
     deallocate(ihrmop)
     deallocate(ifckop)
     deallocate(isymop)
     deallocate(dmat)

  end subroutine get_2e_fock_matrix_x2c
!**********************************************************************

  subroutine get_2e_2c_fock_matrix_x2c(                            &
                                       Gmat,                       &
                                       Dmat,                       &
                                       nrows,                      &
                                       ncols,                      &
                                       nz,                         &
                                       nr_2e_fock_matrices,        &
                                       is_len_wrk,                 &
                                       print_lvl,                  &
                                       aoc_factors                 &
                                       )
!**********************************************************************
!
!    purpose: get 2e-fock matrix
!
!**********************************************************************
#ifdef VAR_MPI
#include "mpif.h"
#endif
     real(8), intent(inout)   :: Gmat(nrows,ncols,nz,nr_2e_fock_matrices)
     real(8), intent(inout)   :: Dmat(nrows,ncols,nz,nr_2e_fock_matrices)
     integer, intent(in)      :: nrows
     integer, intent(in)      :: ncols
     integer, intent(in)      :: nz
     integer, intent(in)      :: nr_2e_fock_matrices
     integer, intent(in)      :: is_len_wrk
     integer, intent(in)      :: print_lvl
     real(8), intent(in)      :: aoc_factors(*)
!----------------------------------------------------------------------
     integer                  :: i, npos, local_lwork, lucoef
     integer, parameter       :: intflg = 1
     real(8)                  :: xx, j
     real(8), allocatable     :: work(:)
     integer, allocatable     :: isymop(:), ifckop(:), ihrmop(:), pos(:)
     logical, allocatable     :: saveflags(:)
     logical                  :: parcal
!**********************************************************************

     allocate(isymop(nr_2e_fock_matrices))
     allocate(ifckop(nr_2e_fock_matrices))
     allocate(ihrmop(nr_2e_fock_matrices))

     ihrmop = 0
     ifckop = 0
     isymop = 0

     do i = 1,nr_2e_fock_matrices
!      totally symmetric operator
       isymop(i) = 1
!      fock matrix type
       ifckop(i) = 1
!      hermitian operator
       ihrmop(i) = 1
     end do

     allocate(saveflags(4))

     call SaveTaskDistribFlags(saveflags)
     call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. ,.TRUE. /))
     parcal = .false.
#ifdef VAR_MPI
     call mpi_initialized(parcal,i)
#endif
     call SetIntTaskArrayDimension(npos,parcal)
     allocate(pos(npos))

     local_lwork = is_len_wrk - 3*nrows*ncols*nz*nr_2e_fock_matrices
     allocate(work(local_lwork))

!    2e-fock matrices (Dirac driver routine)
     call twofck(isymop,ihrmop,ifckop,Gmat,Dmat,nr_2e_fock_matrices,   &
                 pos,intflg,print_lvl,work,local_lwork)

!     accumulate active fock matrices (add to inactive fock matrix)
      if(nr_2e_fock_matrices > 1)then
        do i = 1, nr_2e_fock_matrices-1
          call daxpy(nrows*ncols*nz,aoc_factors(i),Gmat(1,1,1,i+1),1,Gmat(1,1,1,1),1)
        end do
     end if

     deallocate(work)
     if(parcal) call SetTaskDistribFlags(saveflags)
     deallocate(pos)
     deallocate(saveflags)
     deallocate(ihrmop)
     deallocate(ifckop)
     deallocate(isymop)

  end subroutine get_2e_2c_fock_matrix_x2c
!**********************************************************************

  subroutine hello_dirx2c()
!**********************************************************************

    print '(/18x,a)', ' *********************************************************************'
    print '(18x,a )', ' ***   Entering the Exact-Two-Component (X2C) interface in DIRAC   ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** library version:  1.2 (August  2013)                          ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** authors:          - Stefan Knecht                             ***'
    print '(18x,a )', ' ***                   - Trond Saue                                ***'
    print '(18x,a )', ' *** contributors:     - Hans Joergen Aagaard Jensen               ***'
    print '(18x,a )', ' ***                   - Michal Repisky                            ***'
    print '(18x,a )', ' ***                   - Miroslav Ilias                            ***'
    print '(18x,a )', ' *** features:         - X2C                                       ***'
    print '(18x,a )', ' ***                   - X2C-atomic/fragment (X2C-LU)              ***'
    print '(18x,a )', ' ***                   - X2C-spinfree                              ***'
    print '(18x,a )', ' ***                   - X2C-molecular-mean-field (X2Cmmf)         ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' ***                      Universities of                          ***'
    print '(18x,a )', ' ***     Zuerich, Toulouse, Odense, Banska Bystrica and Tromsoe    ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** contact: stefan.knecht@phys.chem.ethz.ch                      ***'
    print '(18x,a/)', ' *********************************************************************'

  end subroutine hello_dirx2c
!**********************************************************************

  subroutine goodbye_dirx2c()
!**********************************************************************

    print '(/18x,a)', ' *********************************************************************'
    print '( 18x,a)', ' ***               X2C transformation ended properly.              ***'
    print '( 18x,a)', ' ***          Calculation continues in two-component mode.         ***'
    print '(18x,a/)', ' *********************************************************************'

  end subroutine goodbye_dirx2c
!**********************************************************************

end module x2c_utils
