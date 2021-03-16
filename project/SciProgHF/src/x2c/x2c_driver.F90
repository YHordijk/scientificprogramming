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

module dirx2c

! stefan: - this module contains all necessary functionality 
!           to perform the transformation from a four-component 
!           to a two-component framework 
!
!           written by sknecht, jan 2010 - july 2010
!           modularization and interface to MDIRAC (InteRest): sknecht and mrepisky aug 2012
!
!           TODO: go through code and carefully check nz vs nzt [a job for TDK aka Bruce Wayne] 
!                 nzt == number of AO2MO transformation matrices 
!                 "MDIRAC" NZT == 1, 
!                 " DIRAC" NZT == NZ
!                 for simplicity we internally use nzt = nz also for MDIRAC (some 0 multiplies will not hurt that much
!                 
! ---------------------------------------------------------------------------------------
!                            "roadmap" of the X2C code
! ---------------------------------------------------------------------------------------
!
!   x2c_driver: top-level routine looping over molecular/atomic fragments
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       1. do fragment = 1 , all_fragments
!
!           --> x2c_tra: for each molecular/atomic fragments:
!    
!               --> a. set up a defining matrix 4c-h1 in orthonormal MO basis
! 
!               --> b. get decoupling matrix "R" in orthonormal MO basis
!                   +  get picture-change transformation matrix U in orthonormal MO basis   
!
!               --> c. construct picture-change transformation matrix U in DIRAC-sorted SA-AO basis
! 
!               --> d. picture-change transform the "basic" 4c-h1 matrix in DIRAC-sorted SA-AO basis
!                      note: the defining matrix 4c-h1 and the "basic" 4c-h1 matrix 
!                            may or may not be identical.
!
!               --> e. (add AMFI contributions to the bare-nucleus 1e Dirac-Hamiltonian in case 
!                      h1 == bare-nucleus 1e Dirac-Hamiltonian (or extended Huckel guess) .and. AMFI is desired)
!
!       end do
!
!       2. (x2c_combine_aba_x2c: approximate a "full" picture-change transformation matrix U from all fragments)
!
!       3. pctrafo_all_op_driver: picture-change transform 
!
!           a. list of 1e property operator 
!
!           b. bare-nucleus corrected full 1e Dirac-Hamiltonian
!
!
!
!   successful X2C procedure: renew common blocks
!   +++++++++++++++++++++++++++++++++++++++++++++
!
! ---------------------------------------------------------------------------------------
!                            end of "roadmap" of the X2C code
! ---------------------------------------------------------------------------------------
!
! x2c modules
  use x2cmod_cfg
  use fragment_x2c_cfg
  use x2c_cb_interface
  use x2c_def_h1_4c
  use x2c_get_h1inf2c_mo
  use x2c_decoupling_mat
  use x2c_h1inf2c_mo2ao
  use x2c_2e_soc_interface
  use x2c_pct_ao
  use x2c_fio
  use x2c_utils
  use x2c_h1mat_base
  use x2c_pct_mo_coefficients
  use x2c_import_exportC1
  use x2c_fragment_x2c_interface

  implicit none

  public x2c_main

  private

  save

  real(8), parameter :: x2c_dummy(2) =  0.0d0


! x2c-mat type definition
! ----------------------------------------------------------------------------
  type x2c_data_type

    real(8), allocatable ::      &
      xmat(:),                   & ! 
      scr1mat(:),                & ! 
      scr2mat(:)                   !

  end type x2c_data_type

! x2c-mat type object
  type(x2c_data_type) :: x2c_data
! ----------------------------------------------------------------------------

contains

!**********************************************************************
  subroutine x2c_main(                        &
                     )
!**********************************************************************
!
!    purpose: main driver routine of the X2C algorithm
!
!             for a detailed overview of the code structure the programmer 
!             may be refered to the "figure" at the beginning of this file.
!             A in-depth theoretical description of the molecular X2C (one-step) 
!             procedure can be found in reference 
!             "M. Ilias and T. Saue, J. Chem. Phys., 126, 064102 (2007)". 
!
!             the driver routine splits into three major parts:
!                 1. loop over (molecular fragments)/ "all centers" for decoupling
!                 2. approximate a combined picture-change transformation matrix U_comb
!                 3. picture-change transform all relevant operator
!
!             possible "defining h1_4c" Hamiltonian:
!
!                 a. x2c_is_defining_h1mat == 0: 1-el Dirac Hamiltonian                           
!                    ==> default for X2C
!
!                 b. x2c_is_defining_h1mat == 2: Fock operator after a (converged) [spinfree==sf] SCF 
!                    ==> default for [sf-]X2Cmmf
!
!                 b. x2c_is_defining_h1mat == 3: free-particle 
!                    ==> test-option
!
!
!**********************************************************************
!----------------------------------------------------------------------
     integer                :: lwrk
     integer                :: final_ham_lvl_x2c        ! level of final hamiltonian stored on global file X2CMAT
     integer                :: dimension_for_2mat
     integer                :: i,j
     character (len=12)     :: flabel                   ! label for matrices on file
     character (len=51)     :: path
!**********************************************************************

!    print header
     call hello_dirx2c()

!    ------------------------------------------------------------------------------
!    step 1.0: initialize X2C dimensions in Dirac-sorted SA-AO basis and symmetry settings 
!              transfer information from (Dirac) common blocks to private X2C variables
!    ------------------------------------------------------------------------------
     call init_x2c_cb(.true.)

!    initialize variables
     final_ham_lvl_x2c = 0
     lwrk              = len_wrk_f77_x2c

!    inform user about the chosen path in X2C module
     path = ''
     select case(fragment_x2c_info%fragment_approach_enabled)
       case(.true.)
                          path = 'fragment/atomic X2C (with spin-orbit contributions)'
         if(x2c_spinfree) path = 'fragment/atomic X2C (spinfree)                     '
       case(.false.)
                          path = 'molecular X2C (with spin-orbit contributions)      '
         if(x2c_spinfree) path = 'molecular X2C (spinfree)                           '
     end select

     print '(/18x,a,a)', ' *** chosen path in X2C module: ',path
 
!    ------------------------------------------------------------------------------
!    step 1.1: initialize and open files
!    ------------------------------------------------------------------------------

!    open file X2CMAT
     open(x2c_funit,file='X2CMAT',status='replace',form='unformatted',    &
          access='sequential',action='readwrite',position='rewind')
     call x2c_write(flabel,x2c_dummy,-1,x2c_funit)

!    open scratch file x2cscr
     open(x2c_fscru,file='x2cscr',status='replace',form='unformatted',  &
          access='sequential',action="readwrite",position='rewind')
     call x2c_write(flabel,x2c_dummy,-1,x2c_fscru)

!    open scratch file X2Camfi_scr
     open(x2c_famfi,file='X2Camfi_scr',status='replace',form='unformatted',  &
          access='sequential',action="readwrite",position='rewind')
     call x2c_write(flabel,x2c_dummy,-1,x2c_famfi)


!    ------------------------------------------------------------------------------
!    step 1.2: obtain for 
!              i. fragment/atomic X2C
!                  a. picture-change transformation matrix "U" in Dirac-sorted SA-AO basis
!                  b. h1_2c_inf (bare-nucleus 1e Dirac-Hamiltonian) in Dirac-sorted SA-AO basis
!             ii. molecular X2C
!                  a. decoupling matrix "R" in orthonormal MO basis
!                  b. picture-change transformation matrix "U" in Dirac-sorted SA-AO basis
!                  c. h1_2c_inf (bare-nucleus 1e Dirac-Hamiltonian) in Dirac-sorted SA-AO basis
!    ------------------------------------------------------------------------------

     if(fragment_x2c_info%fragment_approach_enabled)then

!      construct molecular U matrix from fragments/atoms
!      -------------------------------------------------
       call x2c_fragment_x2c_driver(nr_quat,               &
                                    nr_fsym,               &
                                    x2c_mxatom,            &
                                    type_nuclei,           &
                                    type_charge,           &
                                    nr_symm_indep_cent,    &
                                    nr_degen_nuc_cent,     &
                                    nr_ao_bas_type(1,0),   &
                                    nr_ao_l,               &
                                    nr_ao_all,             &
                                    x2c_funit              &
                                   )

!      get h1_2c_inf in  Dirac-sorted SA-AO basis
!      ------------------------------------------
       call x2c_fragment_get_h1_2c(nr_ao_total_x2c,        &
                                   nr_ao_large_x2c,        &
                                   nr_ao_all,              &
                                   nr_ao_l,                &
                                   nr_fsym,                &
                                   nr_quat,                &
                                   x2c_cb_pq_to_uq,        &
                                   ioff_aomat_x,           &
                                   x2c_add_amfi,           &
                                   final_ham_lvl_x2c,      &
                                   x2c_funit,              &
                                   x2c_famfi,              &
                                   x2c_bs_to_fs,           &
                                   lwrk,                   &
                                   x2c_do_spherical,       &
                                   x2c_mdirac,             &
                                   x2c_prt)

     else

!      molecular (full) X2C 
!      --------------------

!      initialize memory and prepare AO2MO transformation matrices
       call x2c_full_x2c_init(x2c_data,                    &
                              nr_ao_total_x2c,             &
                              nr_ao_all,                   &
                              nr_ao_l,                     &
                              nr_ao_s,                     &
                              dim_eshell,                  &
                              dim_pshell,                  &
                              ioff_aomat_x,                &
                              dimension_for_2mat,          &
                              fullao2momat_dim,            &
                              fullmo2momat_dim,            &
                              fulllowdmat_dim,             &
                              nr_mo_lw_l,                  &
                              nr_mo_lw_s,                  &
                              nr_mo_lw_all,                &
                              nr_fsym,                     &
                              nr_quat,                     &
                              x2c_cb_pq_to_uq,             &
                              x2c_funit,                   &
                              x2c_lowdin_scrF,             &
                              x2c_Sao_scrF,                &
                              x2c_is_defining_h1mat,       &
                              lwrk,                        &
                              x2c_prt,                     &
                              x2c_linsym,                  &
                              x2c_mdirac                   &
                             )

!      get U and h1_2c matrices
       call x2c_full_x2c_driver(x2c_data%xmat,             &
                                x2c_data%scr1mat,          &
                                x2c_data%scr2mat,          &
                                nr_ao_large_x2c,           &
                                nr_ao_total_x2c,           &
                                nr_ao_all,                 &
                                nr_ao_l,                   &
                                dim_oshell,                &
                                dim_eshell,                &
                                dim_pshell,                &
                                fullomat_dim,              &
                                fulleomat_dim,             &
                                ioff_aomat_x,              &
                                nr_fsym,                   &
                                nr_quat,                   &
                                nzt_x2c,                   &
                                x2c_cb_pq_to_uq,           &
                                x2c_funit,                 &
                                x2c_fscru,                 &
                                x2c_famfi,                 &
                                x2c_is_defining_h1mat,     &
                                final_ham_lvl_x2c,         &
                                x2c_add_amfi,              &
                                x2c_bs_to_fs,              &
                                lwrk,                      &
                                x2c_prt,                   &
                                fragment_x2c_info%         &
                                fragment_approach_enabled, &
                                x2c_linsym,                &
                                x2c_mdirac,                &
                                x2c_spinfree,              &
                                x2c_do_spherical)
!      free memory
       deallocate(x2c_data%scr1mat)
       deallocate(x2c_data%scr2mat)
       deallocate(x2c_data%xmat)

!      close temporary file
       close(x2c_Sao_scrF,status="delete")

     end if

      
!    --------------------------------------------------------------------------------------------------
!    step 1.3: add oo-order 2e-SOC or AMFI contributions in case
!              h1 == bare-nucleus 1e Dirac-Hamiltonian .and. 2e-SOC are desired
!    --------------------------------------------------------------------------------------------------
     if(x2c_add_amfi > 0)then

       call x2c_add_2e_so_corrections(final_ham_lvl_x2c,                 &
                                      x2c_add_amfi,                      &
                                      nr_ao_total_x2c,                   &
                                      nr_ao_large_x2c,                   &
                                      nr_ao_all,                         &
                                      nr_ao_l,                           &
                                      nr_fsym,                           &
                                      nr_quat,                           &
                                      x2c_bs_irrep_mat,                  &
                                      x2c_cb_pq_to_uq,                   &
                                      x2c_cb_uq_to_pq,                   &
                                      x2c_iqmult_trip_q,                 &
                                      x2c_bs_to_fs,                      &
                                      x2c_pointer_quat,                  &
                                      x2c_pointer_quat_op,               &
                                      x2c_max_quant_num,                 &
                                      num_nuclei,                        &
                                      type_nuclei,                       &
                                      x2c_tot_charge_for_mfsum,          &
                                      x2c_amfi_order,                    &
                                      lwrk,                              &
                                      x2c_do_spherical,                  &
                                      x2c_famfi,                         &
                                      x2c_funit,                         &
                                      x2c_cspeed,                        &
                                      x2c_prt)
      end if

      close(x2c_famfi,status="delete")

!     ------------------------------------------------------------------------------------------
!     step 2: picture-change transform all (remaining) operator (in DIRAC-sorted SA-AO basis):
!             a. if(appropriate) bare-nucleus corrected 1e Dirac-Hamiltonian
!             b. if(appropriate) 1-electron/2-electron Fock matrices
!             c. list of 1e-property operators

!             in : i.  on file X2CMAT   (fh: x2c_funit):
!                     (combined) picture-change transformation matrix "U"
!                  ii. on file AOPROPER (fh: 1int_4c_file): 
!                      all property integrals in Hermit-sorted AO basis
!
!             out: all operator (listed above) in picture-change transformed form on file
!     -----------------------------------------------------------------------------------
      call pctrafo_all_op_driver(nr_ao_total_x2c,               &
                                 nr_ao_large_x2c,               &
                                 nr_ao_all,                     &
                                 nr_ao_l,                       &
                                 nr_fsym,                       &
                                 nr_quat,                       &
                                 x2c_bs_irrep_mat,              &
                                 x2c_cb_pq_to_uq,               &
                                 x2c_cb_uq_to_pq,               &
                                 x2c_iqmult_trip_q,             &
                                 x2c_pointer_quat,              &
                                 x2c_qdef,                      &
                                 ioff_aomat_x,                  &
                                 x2c_funit,                     &
                                 fh_1int_4c,                    &
                                 x2c_bs_to_fs,                  &
                                 all_prp_op,                    &
                                 lwrk,                          &
                                 x2c_do_spherical,              &
                                 x2c_mdirac,                    &
                                 x2c_prep_bnccorr,              &
                                 x2c_1e_system,                 &
                                 x2c_is_defining_h1mat,         &
                                 x2c_dfopen,                    &
                                 nr_2e_fock_matrices,           &
                                 file_name_1e_fock_matrix,      &
                                 file_name_2e_fock_matrices,    &
                                 x2c_prt)
!      ---------------------------------------------------------------------------------------------------------------
!      step 3: picture-change transform the 4c-MO coefficients (defining
!      operator matrix is the 4c-Dirac Fock operator)
!      ---------------------------------------------------------------------------------------------------------------
       if(x2c_is_defining_h1mat == 2)then ! X2Cmmf
         call x2c_pctrafo_mo_coefficients_driver('DFCOEF',                  &
                                                 2,                         &
                                                 nr_cmo_q,                  &
                                                 nr_ao_total_x2c,           &
                                                 nr_ao_large_x2c,           &
                                                 nr_ao_all,                 &
                                                 nr_ao_l,                   &
                                                 ioff_aomat_x,              &
                                                 dimension_for_2mat,        &
                                                 dim_oshell,                &
                                                 dim_eshell,                &
                                                 dim_pshell,                &
                                                 nr_fsym,                   &
                                                 nr_quat,                   &
                                                 x2c_cb_pq_to_uq,           &
                                                 x2c_funit,                 &
                                                 x2c_linsym,                &
                                                 x2c_prt)
       end if

!      --------------------------------------------------------------------------------------------
!      step 4: export picture-change transformation matrix to AO basis with nosym (C1 double group)
!      --------------------------------------------------------------------------------------------
       if(x2cmod_h1dirty)then
         call x2c_export_pctmat2C1()
       end if

!      ----------------------------------------------------------------------------------------
!      step 5: renew dimensions on common blocks, i.e. "remove" references to positronic shells
!      ----------------------------------------------------------------------------------------
       call renew_x2c_cb_orb_shell_dim(final_ham_lvl_x2c)

       close(x2c_funit,status="keep")

!      print goodbye 
       call goodbye_dirx2c()

  end subroutine x2c_main

!**********************************************************************
  subroutine x2c_full_x2c_init(A,                    &
                               naosh_ls,             &
                               naosh_all,            &
                               naosh_L,              &
                               naosh_S,              &
                               nesh_dim,             &
                               npsh_dim,             &
                               ioff_aomx,            &
                               dimension_for_2mat,   &
                               fullao2momat_dim,     &
                               fullmo2momat_dim,     &
                               fulllowdmat_dim,      &
                               nmo_lw_l,             &
                               nmo_lw_s,             &
                               nmo_lw_all,           &
                               nfsym,                &
                               nz,                   &
                               ipq_off,              &
                               x2c_file_glb,         &
                               lowdin_scrF,          &
                               Sao_scrF,             &
                               defining_h1mat,       &
                               lwrk,                 &
                               print_lvl,            &
                               linear_sym,           &
                               mdirac                &
                              )              
!**********************************************************************
!
!    purpose: full (== molecular approach) X2C transformation initialization routine
!
!             requires: 4c-Lowdin matrix; 
!                       value of "defining_h1mat" determines the defining 4c-h1 matrix
!                       wrt which the decoupling step is performed.
!
!                       defining_h1mat == 0: bare-nucleus one-electron Dirac-Hamiltonian
!                                      == 2: converged 4c-Fock-Dirac operator
!
!             delivers: a. AO2MO transformation matrices
!                       b. allocates matrices with proper dimensions
!
!**********************************************************************
     type(x2c_data_type)    :: A
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: naosh_L(nfsym)
     integer, intent(in)    :: naosh_S(nfsym)
     integer, intent(in)    :: nesh_dim(nfsym)
     integer, intent(in)    :: npsh_dim(nfsym)
     integer, intent(in)    :: nmo_lw_l(nfsym)
     integer, intent(in)    :: nmo_lw_s(nfsym)
     integer, intent(in)    :: nmo_lw_all(nfsym)
     integer, intent(out)   :: dimension_for_2mat
     integer, intent(in)    :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)    :: fullao2momat_dim
     integer, intent(in)    :: fullmo2momat_dim
     integer, intent(in)    :: fulllowdmat_dim
     integer, intent(in)    :: ipq_off(4,0:7)
     integer, intent(in)    :: x2c_file_glb
     integer, intent(in)    :: lowdin_scrF
     integer, intent(in)    :: Sao_scrF
     integer, intent(in)    :: defining_h1mat
     integer, intent(inout) :: lwrk
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: linear_sym
     logical, intent(in)    :: mdirac
!----------------------------------------------------------------------
     real(8), allocatable   :: eigvl(:)
     real(8), allocatable   :: wrk(:)
     integer                :: is_allocated_dyn
     integer                :: local_lwrk
     integer                :: fulllowdmat_dim_tmp
     integer                :: dim1_fsym_Sao
     integer                :: dim2_fsym_Sao
     integer                :: i,j
     character (len=12)     :: flabel             ! label for matrices on file
!**********************************************************************

!    memory allocation - part 1 - allocate necessary matrices and arrays - xmat for step 1.1 + 1.2
!    --------------------------------------------------------------------------------------------- 
     dimension_for_2mat   = naosh_ls**2 * nz
     do i = 1, nfsym
       dimension_for_2mat = max(dimension_for_2mat,ioff_aomx(i,i)+1+naosh_ls*naosh_all(i)*nz)
     end do

     allocate(A%xmat(dimension_for_2mat))

     A%xmat            = 0.0d0
     is_allocated_dyn  = dimension_for_2mat

!    ------------------------------------------------------------------------------
!    step 1: compute 4c-Sao_LL overlap matrix
!    ------------------------------------------------------------------------------

!    open temporary scratch file for 4c-Sao_LL overlap matrix
     open(Sao_scrF,file='X2C_Sao_LL_scr',status='replace',form='unformatted',    &
          access='sequential',action="readwrite",position='rewind')
     call x2c_write(flabel,x2c_dummy,-1,x2c_Sao_scrF)

     dim1_fsym_Sao = naosh_all(1)
     dim2_fsym_Sao = 0
     if(nr_fsym == 2) dim2_fsym_Sao = naosh_all(nr_fsym)

     call get_Sao_LL_mat(A%xmat,mdirac,dim1_fsym_Sao,nfsym,dim2_fsym_Sao)
!                                                                  ^ (compute+write if >= 0)

!    allocate remaining required matrices and arrays: TODO: check if dimensions can actually be reduced...
     j                               = dimension_for_2mat
     if(defining_h1mat < 1) j = naosh_ls**2 * nz
     allocate(A%scr2mat(naosh_ls**2 * nz))
     allocate(A%scr1mat(j))

     allocate(eigvl(naosh_ls))

     A%scr1mat = 0.0d0
     A%scr2mat = 0.0d0
     eigvl     = 0.0d0

!    set maximum allowed length of wrk array allocations
     is_allocated_dyn = is_allocated_dyn + naosh_ls         + naosh_ls**2 * nz + j
     local_lwrk       = lwrk
     local_lwrk       = local_lwrk       - is_allocated_dyn

!    ------------------------------------------------------------------------------
!    step 2: get LOWDIN matrix
!              1.2.a. get the 4c-Lowdin matrix
!              1.2.b. get the 2c-Lowdin matrix (only needed if the defining Hamiltonian is the 4c-Fock operator)
!    ------------------------------------------------------------------------------

!    a. read/get the 4c-Lowdin matrix
!    --------------------------------
     open(lowdin_scrF,file='LOWDMAT',status='old',form='unformatted',    &
          access='sequential',action="readwrite",position='rewind')

     call readi(lowdin_scrF,1,fulllowdmat_dim_tmp)
     call readt(lowdin_scrF,fulllowdmat_dim,A%xmat)

!    consistency check for LOWDIN matrix
     if(fulllowdmat_dim_tmp /= fulllowdmat_dim)then
       call quit('error in x2c-driver: dimension of existing LOWDIN matrix does not match'// & 
                 ' the present LOWDIN dimensions.')
     end if
  
!    b. construct the two-component Lowdin matrix
!    --------------------------------------------
     if(defining_h1mat == 2)then

       call x2c_lowd2c(A%xmat,               &
                       A%scr1mat,            &
                       A%scr2mat,            &
                       npsh_dim,             &
                       nesh_dim,             &
                       naosh_L,              &
                       naosh_S,              &
                       naosh_all,            &
                       nmo_lw_l,             &
                       nmo_lw_s,             &
                       nmo_lw_all,           &
                       naosh_ls,             &
                       nfsym,                &
                       nz,                   &
                       ipq_off,              &
                       x2c_file_glb,         &
                       print_lvl,            &
                       linear_sym)

!      recover the 4c-Lowdin matrix
       A%xmat = 0.0d0
       rewind lowdin_scrF
       call readi(lowdin_scrF,1,fulllowdmat_dim)
       call readt(lowdin_scrF,fulllowdmat_dim,A%xmat)

     end if

     close(lowdin_scrF,status="delete")
!    ------------------------------------------------------------------------------
!    step 3: construct modified dirac equation in ON basis --> retrieve AOtoMO transformation matrices.
!    ------------------------------------------------------------------------------
!
!          in : 4c-LOWDIN matrix -> xmat
!          out: AOtoMO transformation matrices -> xmat, scr1mat
!                                               + (scr2mat, if x2c_is_defining_h1mat == 2 .and. x2c_linsym == .true.) 
!    ------------------------------------------------------------------------------
     select case(defining_h1mat)

       case(2) ! defining h1 hamiltonian is the 4c-fock operator obtained in a preceeding 4c-SCF run

         A%scr1mat = 0.0d0
         A%scr2mat = 0.0d0
         A%xmat    = 0.0d0

!        a. read the (linear-symmetry adapted) AO2MO transformation matrix
!        -----------------------------------------------------------------
         open(11,file='AOMOlin',status='old',form='unformatted', &
              access='sequential',action="readwrite",position='rewind')

         if(linear_sym)then
!          MO2MO linear symmetry --> nonlinear symmetry transformation matrix
           call readt(11,fullmo2momat_dim,A%scr2mat)
!          AO2MO linear symmetry --> nonlinear symmetry transformation matrix
!          skip the record as it is not needed here
           read(11)
         end if

!        AO2MO transformation matrix
         call readt(11,fullao2momat_dim,A%scr1mat)

         close(11,status='delete')

!        b. read the SL-resorted AO2MO transformation matrix
!        ---------------------------------------------------
         open(9,file='AOMOSLR',status='old',form='unformatted', &
              access='sequential',action="readwrite",position='rewind')

         call readt(9,fullao2momat_dim,A%xmat)

         close(9,status='delete')

       case (0, 3) ! 1-el Dirac Hamiltonian / free particle

         allocate(wrk(local_lwrk))
         call modham(A%scr1mat,                                         &
                     A%xmat,                                            &
                     A%scr2mat,                                         &
                     eigvl,                                             &
                     .true.,                                            &
                     .false.,                                           &
                     defining_h1mat,                                    &
                     wrk,                                               &
                     local_lwrk)
         deallocate(wrk)
!        ------------------------------------------------------------------------------------
!        step 4: transfer dimensions on common blocks concerning the orthonormal MO basis -
!                these dimensions have been evaluated in MODHAM
!        ------------------------------------------------------------------------------------
         call init_x2c_cb(.true.)

       case default
     
         stop 'x2c main: defining 4c-Hamiltonian type unknown...'

     end select

!    memory allocation part 2 - free top-level scratch memory and possibly re-allocate scr1mat with appropriate size
!    ---------------------------------------------------------------------------------------------------------------
     deallocate(eigvl)
     if(defining_h1mat < 1)then
       deallocate(A%scr1mat)
!      reallocate with full size
       allocate(A%scr1mat(dimension_for_2mat))
       A%scr1mat = 0.0d0
     end if
!
!    recalculate the new maximum allowed length of wrk array allocations to be used in the subroutines.
     is_allocated_dyn = is_allocated_dyn - naosh_ls         - j + naosh_ls**2 * nz
     lwrk             = lwrk             - is_allocated_dyn
!
  end subroutine x2c_full_x2c_init
!**********************************************************************

  subroutine x2c_full_x2c_driver(xmat,                 &
                                 scr1mat,              &
                                 scr2mat,              &
                                 nr_ao_L,              &
                                 naosh_ls,             &
                                 naosh_all,            &
                                 naosh_L,              &
                                 norb_dim,             &
                                 nesh_dim,             &
                                 npsh_dim,             &
                                 h1mat_dim_4c,         &
                                 momat_orbesh_dim_4c,  &
                                 ioff_aomx,            &
                                 nfsym,                &
                                 nz,                   &
                                 nzt,                  &
                                 ipq_off,              &
                                 x2c_file_glb,         &
                                 x2c_file_scr,         &
                                 x2c_file_amf,         &
                                 defining_h1mat,       &
                                 is_final_ham_lvl,     &
                                 add_amfi_contrib,     &
                                 op_bs_to_fs,          &
                                 is_len_wrk,           &
                                 print_lvl,            &
                                 fragment_x2c,         &
                                 linear_sym,           &
                                 mdirac,               &
                                 spinfree,             &
                                 spherical_on)
!**********************************************************************
!
!    purpose: X2C transformation driver routine
!
!             requires: 4c-Lowdin matrix; 
!                       value of "defining_h1mat" determines the defining 4c-h1 matrix
!                       wrt which the decoupling step is performed.
!
!                       defining_h1mat == 0: bare-nucleus one-electron Dirac-Hamiltonian
!                                      == 1: (extended) Huckel start guess (won't be implemented)
!                                      == 2: converged 4c-Fock-Dirac operator
!
!             delivers: a. decoupling matrix "R" in orthonormal MO basis
!                       b. picture-change transformation matrix "U" in orthonormal MO basis and Dirac-sorted SA-AO basis
!                       c. h1_2c_inf in Dirac-sorted SA-AO basis
!
!    for more information see also the equation (9), (10), (23) and (24) in the
!    reference: M. Ilias and T. Saue, J. Chem. Phys., 126, 064102 (2007)
!
!**********************************************************************
     real(8), intent(inout) :: xmat(*)
     real(8), intent(inout) :: scr1mat(*)
     real(8), intent(inout) :: scr2mat(*)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: nzt
     integer, intent(in)    :: nr_ao_L
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: naosh_L(nfsym)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: nesh_dim(nfsym)
     integer, intent(in)    :: npsh_dim(nfsym)
     integer, intent(in)    :: h1mat_dim_4c
     integer, intent(in)    :: momat_orbesh_dim_4c
     integer, intent(in)    :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)    :: ipq_off(4,0:7)
     integer, intent(in)    :: x2c_file_glb
     integer, intent(in)    :: x2c_file_scr
     integer, intent(in)    :: x2c_file_amf
     integer, intent(in)    :: defining_h1mat
     integer, intent(inout) :: is_final_ham_lvl
     integer, intent(in)    :: add_amfi_contrib
     integer, intent(in)    :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)    :: is_len_wrk
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: fragment_x2c
     logical, intent(in)    :: linear_sym
     logical, intent(in)    :: mdirac
     logical, intent(in)    :: spinfree
     logical, intent(in)    :: spherical_on
!----------------------------------------------------------------------
     real(8), allocatable   :: eigvl(:)
     integer                :: i,j
     integer                :: lwrk
     integer                :: basic_h1mat        ! defining parameter for the basic 4c-h1 which will be transformed to 2c-h1
     character (len=12)     :: flabel             ! label for matrices on file
!**********************************************************************
!

!      initialize maximum allowed length of wrk array allocations
       lwrk = is_len_wrk

!      step 1: save RKB AOtoMO transformation matrices to file x2cscr
!
!              in : AOtoMO transformation matrices
!              out: ["SL resorted"] matrix -> scr2mat
!              on file x2cscr: 
!                 a. AOtoMO transformation matrix ["SL resorted"] under the label 4cAOMOr               --> xmat
!                 b. AOtoMO transformation matrix under the label 4cAOMOo                               --> scr1mat
!                 c. MOtoMO non-linear -> linear symmetry transformation matrix under the label 4cMOMOl --> scr2mat
!      ------------------------------------------------------------------------------
       call store_aomo_trafo_matrices(xmat,                &
                                      scr1mat,             &
                                      scr2mat,             &
                                      norb_dim,            &
                                      naosh_all,           &
                                      nfsym,               &
                                      nz,                  &
                                      nzt,                 &
                                      defining_h1mat,      &
                                      x2c_file_scr,        &
                                      ipq_off,             &
                                      print_lvl,           &
                                      linear_sym)

!      step 2: store the hamiltonian for which the 4c-->2c transformation 
!              should be done on file (== basic "h1_4c" in orthonormal MO basis) h1_4cON
!               in: RKB AOtoMO transformation matrix ["SL resorted"] -> scr2mat
!              out: - (nothing in matrices)
!                   on file x2cscr: basic 4c-h1 in orthonormal basis under label "h1_4cON"
!      ------------------------------------------------------------------------------
!      the basic "h1_4c" will always be the bare-nucleus h1 hamiltonian except 
!      for the "4c-Dirac-Fock operator as basis" option. sknecht+tsaue: march 2011
                               basic_h1mat = 0
       if(defining_h1mat == 2) basic_h1mat = 2
       call x2c_get_h1_mat_base(xmat,                  &
                                scr1mat,               &
                                scr2mat,               &
                                naosh_ls,              &
                                naosh_all,             &
                                norb_dim,              &
                                ioff_aomx,             &
                                nfsym,                 &
                                nz,                    &
                                nzt,                   &
                                ipq_off,               &
                                basic_h1mat,           &
                                lwrk,                  &
                                x2c_file_scr,          &
                                linear_sym,            &
                                mdirac,                &
                                spinfree,              &
                                print_lvl)
       
!      step 3: get defining matrix (== "h1-defining" in orthonormal MO basis)
!              in : ["SL resorted"] matrix -> scr2mat
!              out: defining 4c-h1 in orthonormal basis -> xmat
!      ------------------------------------------------------------------------------
       call x2c_get_h1_defining_mat(xmat,                  &
                                    scr1mat,               &
                                    scr2mat,               &
                                    naosh_ls,              &
                                    naosh_all,             &
                                    norb_dim,              &
                                    ioff_aomx,             &
                                    nfsym,                 &
                                    nz,                    &
                                    nzt,                   &
                                    ipq_off,               &
                                    defining_h1mat,        &
                                    lwrk,                  &
                                    x2c_file_scr,          &
                                    linear_sym,            &
                                    mdirac,                &
                                    spinfree,              &
                                    print_lvl)
 
!      allocate scratch matrix with appropriate dimension
       allocate(eigvl(max(naosh_ls,momat_orbesh_dim_4c,h1mat_dim_4c)))
       eigvl   = 0.0d0

!      set new maximum allowed length of wrk array allocations
       lwrk = lwrk - max(naosh_ls,momat_orbesh_dim_4c,h1mat_dim_4c)

!      step 4: get the transformed basic operator "h1_4c" --> "h1_infinite_2c" in orthonormal basis
!
!              in : defining  4c-h1 in orthonormal basis -> xmat
!                   on file x2cscr in orthonormal basis: the basic operator "h1_4c" (label == h1_4cON)
!              out: the basic 2c-h1 in orthonormal basis -> xmat
!                   on file x2cscr in orthonormal basis: pct-matrix U, "R"-matrix
!      ------------------------------------------------------------------------------
       call get_h1_2c_onmo_bas(xmat,                      &
                               scr1mat,                   &
                               scr2mat,                   &
                               eigvl,                     &
                               norb_dim,                  &
                               nesh_dim,                  &
                               npsh_dim,                  &
                               nfsym,                     &
                               nz,                        &
                               ipq_off,                   &
                               x2c_file_scr,              &
                               defining_h1mat,            &
                               spinfree,                  &
                               print_lvl)

!      step 5: get the transformed operator "h1_infinite_2c" in SA-AO basis
!
!              in : 2c-h1 in orthonormal MO basis -> xmat
!              out: 2c-h1 in SA-AO basis -> xmat + *
!                   
!              * == on file X2CMAT     : under the label h12cAOn (written in h1_2c^++ format)
!                         or
!              * == on file X2Camfi_scr: under the label h12cAOa (written in 4c h1_2c^++ format) 
!                                        in preparation for AMFI contributions.
!                                              
!      ------------------------------------------------------------------------------
       call get_h1_inf_2c_saao_bas(xmat,                    &
                                   eigvl,                   &
                                   scr1mat,                 &
                                   scr2mat,                 &
                                   norb_dim,                &
                                   nesh_dim,                &
                                   naosh_L,                 &
                                   naosh_all,               &
                                   nr_ao_L,                 &
                                   naosh_ls,                &
                                   ioff_aomx,               &
                                   nfsym,                   &
                                   nz,                      &
                                   nzt,                     &
                                   ipq_off,                 &
                                   is_final_ham_lvl,        &
                                   add_amfi_contrib,        &
                                   lwrk,                    &
                                   x2c_file_scr,            &
                                   x2c_file_glb,            &
                                   x2c_file_amf,            &
                                   op_bs_to_fs,             &
                                   spherical_on,            &
                                   fragment_x2c,            &
                                   basic_h1mat,             &
                                   print_lvl)

!      step 6: construct the picture-change transformation matrix in SA-AO basis
!
!              in : -
!              out: -                         
!                   on file X2CMAT: picture-change transformation matrix U in SA-AO basis 
!                                   under the label pctmtAO
!      ------------------------------------------------------------------------------
       call construct_pctmat_saao_bas(xmat,                &
                                      scr1mat,             &
                                      scr2mat,             &
                                      eigvl,               &
                                      norb_dim,            &
                                      nesh_dim,            &
                                      naosh_all,           &
                                      naosh_L,             &
                                      nfsym,               &
                                      nz,                  &
                                      nzt,                 &
                                      ipq_off,             &
                                      x2c_file_scr,        &
                                      x2c_file_glb,        &
                                      print_lvl)

!      free scratch memory
       deallocate(eigvl)

       close(x2c_file_scr,status="delete")

  end subroutine x2c_full_x2c_driver

!**********************************************************************

end module dirx2c
