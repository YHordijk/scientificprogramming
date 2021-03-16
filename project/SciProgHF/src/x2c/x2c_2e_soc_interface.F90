#define DEBUG_SOC
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
! module containing the driver routines for the addition of atomic-mean-field
! two-electron spin-orbit corrections to the hamiltonian in DIRAC-sorted SA-AO basis. 
!
! written by sknecht july 2010
!
module x2c_2e_soc_interface

  use x2c_fio
  use x2c_utils, only:                 &
      print_x2cmat
  use picture_change_operations
#ifdef MOD_AOOSOC
  use atomic_oo_order_so_correction
#endif

  implicit none

  public x2c_add_2e_so_corrections

  private

contains

!**********************************************************************
  subroutine x2c_add_2e_so_corrections(is_final_ham_lvl,        &
                                       x2c_add_amfi,            &
                                       naosh_ls,                &
                                       nr_ao_L,                 &
                                       naosh_all,               &
                                       naosh_L,                 &
                                       nfsym,                   &
                                       nz,                      &
                                       bs_irrep_mat,            &
                                       ipq_off_in,              &
                                       iqp_off_in,              &
                                       iqmult_trip_q,           &
                                       op_bs_to_fs,             &
                                       quat_pointer,            &
                                       quat_pointer_op,         &
                                       max_quant_num,           &
                                       nuclei_totnum,           &
                                       nuclei_type,             &
                                       tot_charge_for_mfsum,    &
                                       amfi_order_contrib,      &
                                       is_len_wrk,              &
                                       spherical_on,            &
                                       x2c_file_amf,            &
                                       x2c_file_glb,            &
                                       cspeed,                  &
                                       print_lvl)
!**********************************************************************
!
!    purpose: driver routine pointing to the appropriate Dirac-2e-SOC/external 
!             AMFI interface routine according to the value x2c_add_amfi.
!             The resulting hamiltonian h1+AMFI is stored on the global file 
!             X2CMAT and the "hamiltonian level" is set to its appropriate 
!             value (see read_1fock_x2c/dirx2c_utils.F90 for a detailed list 
!             of "hamiltonian levels").
!
!----------------------------------------------------------------------
     real(8), intent(in)     :: cspeed
     integer, intent(inout)  :: is_final_ham_lvl
     integer, intent(in)     :: x2c_add_amfi
     integer, intent(in)     :: naosh_ls
     integer, intent(in)     :: nr_ao_L
     integer, intent(in)     :: naosh_all(nfsym)
     integer, intent(in)     :: naosh_L(nfsym)
     integer, intent(in)     :: nfsym
     integer, intent(in)     :: nz
     integer, intent(in)     :: ipq_off_in(4,0:7)
     integer, intent(in)     :: iqp_off_in(4,0:7)
     integer, intent(in)     :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)     :: bs_irrep_mat(4,0:7)
     integer, intent(in)     :: iqmult_trip_q(4,4,4)
     integer, intent(in)     :: quat_pointer(0:7,2)
     integer, intent(in)     :: quat_pointer_op(0:7)
     integer, intent(in)     :: x2c_file_amf
     integer, intent(in)     :: x2c_file_glb
     integer, intent(in)     :: is_len_wrk
     integer, intent(in)     :: max_quant_num
     integer, intent(in)     :: nuclei_totnum
     integer, intent(in)     :: nuclei_type
     integer, intent(in)     :: tot_charge_for_mfsum
     integer, intent(in)     :: amfi_order_contrib
     integer, intent(in)     :: print_lvl
     logical, intent(in)     :: spherical_on
!----------------------------------------------------------------------
     integer                 :: i
!**********************************************************************
 
       if(x2c_add_amfi == 1)then
!        interface to "old" AMFI (by B. Schimmelpfennig)
         call x2c_interface2old_amfi(naosh_ls,             &
                                     nr_ao_L,              &
                                     naosh_all,            &
                                     naosh_L,              &
                                     nfsym,                &
                                     nz,                   &
                                     bs_irrep_mat,         &
                                     ipq_off_in,           &
                                     iqp_off_in,           &
                                     iqmult_trip_q,        &
                                     op_bs_to_fs,          &
                                     quat_pointer,         &
                                     quat_pointer_op,      &
                                     cspeed,               &
                                     max_quant_num,        &
                                     nuclei_totnum,        &
                                     tot_charge_for_mfsum, &
                                     amfi_order_contrib,   &
                                     is_len_wrk,           &
                                     spherical_on,         &
                                     x2c_file_amf,         &
                                     x2c_file_glb,         &
                                     print_lvl)
       else
#ifdef MOD_AOOSOC
         call x2c_interface2aoosoc(naosh_ls,             &
                                   nr_ao_L,              &
                                   naosh_all,            &
                                   naosh_L,              &
                                   nfsym,                &
                                   nz,                   &
                                   bs_irrep_mat,         &
                                   ipq_off_in,           &
                                   iqp_off_in,           &
                                   iqmult_trip_q,        &
                                   op_bs_to_fs,          &
                                   quat_pointer,         &
                                   cspeed,               &
                                   max_quant_num,        &
                                   nuclei_totnum,        &
                                   nuclei_type,          &
                                   tot_charge_for_mfsum, &
                                   amfi_order_contrib,   &
                                   is_len_wrk,           &
                                   spherical_on,         &
                                   x2c_file_amf,         &
                                   x2c_file_glb,         &
                                   print_lvl)

#else
         call quit('unknown 2e-SOC interface in x2c module')
#endif
       end if

       is_final_ham_lvl = -3

  end subroutine x2c_add_2e_so_corrections

!**********************************************************************
  subroutine x2c_interface2old_amfi(naosh_ls,              &
                                    nr_ao_L,               &
                                    naosh_all,             &
                                    naosh_L,               &
                                    nfsym,                 &
                                    nz,                    &
                                    bs_irrep_mat,          &
                                    ipq_off_in,            &
                                    iqp_off_in,            &
                                    iqmult_trip_q,         &
                                    op_bs_to_fs,           &
                                    quat_pointer,          &
                                    quat_pointer_op,       &
                                    cspeed,                &
                                    max_quant_num,         &
                                    nuclei_totnum,         &
                                    tot_charge_for_mfsum,  &
                                    amfi_order_contrib,    &
                                    is_len_wrk,            &
                                    spherical_on,          &
                                    x2c_file_amf,          &
                                    x2c_file_glb,          &
                                    print_lvl)   
!**********************************************************************
!
!    purpose: interface routine to the old AMFI program by B. Schimmelpfennig. 
!             According to the value of "amfi_order_contrib" two-electron
!             spin-orbit contributions are added to the X2C Hamiltonian.
!
!             the outline of this routine may be sketched as follows:
!
!             1. prepare for spherical harmonics
!             2. calculate AMFI contributions
!            /  a. transform from spherical-gaussian AO basis to cartesian-gaussian AO basis
!   "saw"   /   b. transform from simple AO basis to Hermit-sorted SO-AO basis
!  -------  \   c. transform from Hermit-sorted SO-AO basis to Dirac-sorted AO (SA-AO) basis
!            \  d. insert quaternion phases
!             3. read h1_++_saao_cartesian_gaussian from file X2Camfi_scr (in ++ form but 4c-format retained)
!             4. pick out LL block
!             5. save h1_++_LL_saao_cartesian_gaussian to file X2CMAT (only LL-block retained)

!    note: sphgen modifies the common blocks for 
!          cartesian --> spherical gaussian transformations;
!          a call should therefore be followed by a reset call 
!          or rewrite the code such that it works without these common blocks
!          (taken from the notes by Andre Gomes and Luuk Visscher).
!
!----------------------------------------------------------------------
     real(8), intent(in)             :: cspeed
     integer, intent(in)             :: naosh_ls
     integer, intent(in)             :: nr_ao_L
     integer, intent(in)             :: naosh_all(nfsym)
     integer, intent(in)             :: naosh_L(nfsym)
     integer, intent(in)             :: nfsym
     integer, intent(in)             :: nz
     integer, intent(in)             :: ipq_off_in(4,0:7)
     integer, intent(in)             :: iqp_off_in(4,0:7)
     integer, intent(in)             :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)             :: bs_irrep_mat(4,0:7)
     integer, intent(in)             :: iqmult_trip_q(4,4,4)
     integer, intent(in)             :: quat_pointer(0:7,2)
     integer, intent(in)             :: quat_pointer_op(0:7)
     integer, intent(in)             :: x2c_file_amf
     integer, intent(in)             :: x2c_file_glb
     integer, intent(in)             :: is_len_wrk
     integer, intent(in)             :: max_quant_num
     integer, intent(in)             :: nuclei_totnum
     integer, intent(in)             :: tot_charge_for_mfsum
     integer, intent(in)             :: amfi_order_contrib
     integer, intent(in)             :: print_lvl
     logical, intent(in)             :: spherical_on
!--------------------------------------------------------------
     real(8), allocatable            :: wrk(:)
     real(8), allocatable            :: scr1_mat(:)
     real(8), allocatable            :: scr2_mat(:)
     real(8), allocatable            :: rcha_pnuc(:)
     integer, allocatable            :: itmp_mat1(:)
     integer, allocatable            :: itmp_mat2(:)
     integer, allocatable            :: itmp_mat3(:)
     integer, allocatable            :: itmp_mat4(:)
     integer, allocatable            :: itmp_mat5(:)
     integer, allocatable            :: itmp_mat6(:)
     integer, allocatable            :: itmp_mat7(:)
     integer, allocatable            :: itmp_mat8(:)
     integer, allocatable            :: itmp_mat9(:)
     character (len=4), allocatable  :: clab_int(:)
     integer                         :: sph_fun_L
     integer                         :: nr_sph
     integer                         :: i
     integer                         :: j
     integer                         :: iq
     integer                         :: ipq
     integer                         :: iz
     integer                         :: irepd
     integer                         :: nsphcm
     integer                         :: wf_components
     integer                         :: lwrk
     character (len=12)              :: flabel
!**********************************************************************
 
!      ---------------------------------------
!      step 1. prepare for spherical harmonics
!      ---------------------------------------

       lwrk = is_len_wrk
       allocate(wrk(lwrk))

       call sphgen(1,2,.false.,wrk,lwrk,print_lvl)

!      determine the total number of spherical functions for the large (L => 1) component
       sph_fun_L = nsphcm(1)

       allocate(clab_int(max_quant_num**2))

!      get proper labels (including signs) for s,p,d,f,... functions
       call sphlab(max_quant_num-1,clab_int)


!      allocate temporary matrices for AMFI module
       allocate(itmp_mat1(sph_fun_L)    )
       allocate(itmp_mat2(sph_fun_L)    )
       allocate(itmp_mat3(sph_fun_L)    )
       allocate(itmp_mat4(sph_fun_L)    )
       allocate(itmp_mat5(sph_fun_L)    )
       allocate(itmp_mat6(sph_fun_L)    )
       allocate(itmp_mat7(sph_fun_L*2)  )
       allocate(itmp_mat8(nuclei_totnum))
       allocate(itmp_mat9(nuclei_totnum))
       allocate(rcha_pnuc(nuclei_totnum))

!      ------------------------------
!      step 2. add AMFI contributions
!      ------------------------------
!        out: AMFI (spherical) -> scr1_mat

!      allocate first scratch matrix
       allocate(scr1_mat(naosh_ls*naosh_ls*4))
       scr1_mat = 0.0d0

       call amfiin(scr1_mat,               &
                   naosh_ls,               &
                   itmp_mat2,              &
                   itmp_mat1,              &
                   itmp_mat3,              &
                   itmp_mat8,              &
                   itmp_mat4,              &
                   itmp_mat5,              &
                   itmp_mat6,              &
                   itmp_mat9,              &
                   itmp_mat7,              &
                   rcha_pnuc,              &
                   sph_fun_L,              &
                   clab_int,               &
                   .false.,                &
                   .false.,                &
                   tot_charge_for_mfsum,   &
                   amfi_order_contrib,     &
                   quat_pointer_op,        &
                   cspeed,                 &
                   print_lvl,              &
                   wrk,                    &
                   lwrk)


!      release scratch memory
       deallocate(clab_int)
       deallocate(itmp_mat1)
       deallocate(itmp_mat2)
       deallocate(itmp_mat3)
       deallocate(itmp_mat4)
       deallocate(itmp_mat5)
       deallocate(itmp_mat6)
       deallocate(itmp_mat7)
       deallocate(itmp_mat8)
       deallocate(itmp_mat9)
       deallocate(rcha_pnuc)


!      set up the proper factors for the backtransformation
       call sphgen(1,2,.true.,wrk,lwrk,print_lvl)

!      allocate scratch space
       allocate(scr2_mat(naosh_ls*naosh_ls*4))
       scr2_mat = 0

!      a. transform from spherical Gaussian AO basis to cartesian Gaussian AO basis
!      ----------------------------------------------------------------------------
!        in : AMFI contributions -> scr1_mat
!        out: AMFI contributions -> scr2_mat
!             wf_components = 0: large+small
       wf_components = 0 
       do iz = 1, 4
         call mtaosc(scr2_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
                     naosh_ls,                                   &
     &               scr1_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
                     naosh_ls,                                   &
                     wf_components,                              &
                     wrk,                                        &
                     lwrk,                                       &
                     print_lvl)
       end do

!      reset the common blocks (see note in the header)
       call sphgen(1,2,.false.,wrk,lwrk,print_lvl)

       deallocate(wrk)

!      b. transform from AO-basis to SO-AO basis
!      -----------------------------------------
!        in : h1_AO    -> scr2_mat
!        out: h1_SO-AO -> scr1_mat
       scr1_mat = 0.0d0
       i = op_bs_to_fs(0,1)
       do iz = 1, 4
         irepd = bs_irrep_mat(iz,0)
         iq    = iqmult_trip_q(1,quat_pointer(irepd,i),iz)
         ipq   = iqp_off_in(iq,0)
         call mtaoso(scr2_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
                     scr1_mat(1+(naosh_ls*naosh_ls*(ipq-1))),    &
                     naosh_ls,                                   &
                     irepd,                                      &
                     print_lvl)
       end do

!      c. transform from Hermit-sorted AO to DIRAC-sorted AO basis
!      -----------------------------------------------------------
       call butobs_no_work(scr1_mat,nz)

!      d. insert quaternion phase factors
!      ----------------------------------
       if(nz < 4)then
         do iz = 1, nz                                                
           iq = ipq_off_in(iz,0)
           call q2bphase('F',iq,1,scr1_mat(1+(naosh_ls*naosh_ls*(iz-1))))
         end do           
       end if          

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr1_mat,naosh_ls,naosh_ls,nz,ipq_off_in,'x2c - AMFI contribution',6)
       end if

#ifdef DEBUG_SOC
       open(99,file='soc-contributions',status='replace',form='formatted',  &
       access='sequential',action="readwrite",position='rewind')

       call print_x2cmat(scr1_mat,naosh_ls,naosh_ls,nz,ipq_off_in,'x2c - AMFI contribution',99)
       close(99,status='keep')
#endif

!      -------------------------------------------------------------------------------
!      step 3. read h1_2c_saao_cartesian_gaussian from file and add AMFI contributions 
!      -------------------------------------------------------------------------------

       i = 0
       j = 0
       write(flabel,'(a7,i4,i1)') 'h12cAOa',1,i
       scr2_mat = 0.0d0
       call x2c_read(flabel,scr2_mat,naosh_ls*naosh_ls*nz,x2c_file_amf)

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr2_mat,naosh_ls,naosh_ls,nz,ipq_off_in,'x2c - h1 wo AMFI contribution',6)
       end if


!      add amfi contributions to h1_2c
       call daxpy(naosh_ls*naosh_ls*nz,1.0d0,scr2_mat,1,scr1_mat,1)

       scr2_mat = 0.0d0

!      --------------------------------------------------------
!      step 4. select LL-block of h1_++_saao_cartesian_gaussian
!      --------------------------------------------------------

       call pick_LL_block_saao_bas(scr1_mat,            &
                                   scr2_mat,            &
                                   naosh_L,             &
                                   naosh_all,           &
                                   naosh_ls,            &
                                   nr_ao_L,             &
                                   nfsym,               &
                                   nz,                  &
                                   0,                   &
                                   op_bs_to_fs)

!      -------------------------------------------------------------
!      step 5. write h1_++_LL_saao_cartesian_gaussian to file X2CMAT 
!      -------------------------------------------------------------

       i = 0
       j = 0
       write(flabel,'(a7,i4,i1)') 'h12cAOA',1,i

       if(spherical_on)then

         lwrk = is_len_wrk
         allocate(wrk(lwrk))
         call sph_iotc(scr2_mat,scr1_mat,nr_sph,print_lvl,wrk,lwrk)
         deallocate(wrk)

         call x2c_write(flabel,scr1_mat,nr_sph**2  * nz,x2c_file_glb)
       else

!        debug print
         if(print_lvl > 2)then
           call print_x2cmat(scr2_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - h12c+AMFI',6)
         end if

         call x2c_write(flabel,scr2_mat,nr_ao_L**2 * nz,x2c_file_glb)
       end if

!      release scratch space
       deallocate(scr1_mat)
       deallocate(scr2_mat)

  end subroutine x2c_interface2old_amfi

!**********************************************************************

  subroutine x2c_interface2aoosoc(naosh_ls,              &
                                  nr_ao_L,               &
                                  naosh_all,             &
                                  naosh_L,               &
                                  nfsym,                 &
                                  nz,                    &
                                  bs_irrep_mat,          &
                                  ipq_off_in,            &
                                  iqp_off_in,            &
                                  iqmult_trip_q,         &
                                  op_bs_to_fs,           &
                                  quat_pointer,          &
                                  cspeed,                &
                                  max_quant_num,         &
                                  nuclei_totnum,         &
                                  nuclei_type,           &
                                  tot_charge_for_mfsum,  &
                                  amfi_order_contrib,    &
                                  is_len_wrk,            &
                                  spherical_on,          &
                                  x2c_file_amf,          &
                                  x2c_file_glb,          &
                                  print_lvl)   
!**********************************************************************
!
!    purpose: interface routine to th oo-order 2e-SOC module written by S. Knecht 
!             and friends. 
!
!             the outline of this routine may be sketched as follows:
!
!             1. retrieve 2e-SOC from file
!   "saw"   /   a. transform from simple AO basis to Hermit-sorted SO-AO basis
!  -------  \   b. transform from Hermit-sorted SO-AO basis to Dirac-sorted AO (SA-AO) basis
!            \  c. insert quaternion phases
!             2. read h1_++_saao_cartesian_gaussian from file X2Camfi_scr (in ++ form but 4c-format retained)
!             3. pick out LL block
!             4. save h1_++_LL_saao_cartesian_gaussian to file X2CMAT (only LL-block retained)
!
!             partly written on the transatlantic flight from Buenos Aires to Paris ...
!
!----------------------------------------------------------------------
     real(8), intent(in)             :: cspeed
     integer, intent(in)             :: naosh_ls
     integer, intent(in)             :: nr_ao_L
     integer, intent(in)             :: naosh_all(nfsym)
     integer, intent(in)             :: naosh_L(nfsym)
     integer, intent(in)             :: nfsym
     integer, intent(in)             :: nz
     integer, intent(in)             :: ipq_off_in(4,0:7)
     integer, intent(in)             :: iqp_off_in(4,0:7)
     integer, intent(in)             :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)             :: bs_irrep_mat(4,0:7)
     integer, intent(in)             :: iqmult_trip_q(4,4,4)
     integer, intent(in)             :: quat_pointer(0:7,2)
     integer, intent(in)             :: x2c_file_amf
     integer, intent(in)             :: x2c_file_glb
     integer, intent(in)             :: is_len_wrk
     integer, intent(in)             :: max_quant_num
     integer, intent(in)             :: nuclei_totnum
     integer, intent(in)             :: nuclei_type
     integer, intent(in)             :: tot_charge_for_mfsum
     integer, intent(in)             :: amfi_order_contrib
     integer, intent(in)             :: print_lvl
     logical, intent(in)             :: spherical_on
!--------------------------------------------------------------
     real(8), allocatable            :: scr1_mat(:)
     real(8), allocatable            :: scr2_mat(:)
     real(8), allocatable            :: wrk(:)
     integer, allocatable            :: aooso_info(:)
     integer                         :: i
     integer                         :: j
     integer                         :: iq
     integer                         :: ipq
     integer                         :: iz
     integer                         :: irepd
     integer                         :: nsphcm
     integer                         :: nr_sph
     integer                         :: nrows_frag
     integer                         :: ncols_frag
     integer                         :: center_total
     integer                         :: lwrk
     integer                         :: wf_components
     character (len=12)              :: flabel
!**********************************************************************
 
!      allocate scratch space
       allocate(scr1_mat(naosh_ls*naosh_ls*4))
       allocate(scr2_mat(nr_ao_L*nr_ao_L*4))
       scr1_mat = 0.0d0
       scr2_mat = 0.0d0

!      -------------------------------------------------------------------------------
!      step 1. read h1_2c_saao_cartesian_gaussian from file and add aoosoc contributions 
!      -------------------------------------------------------------------------------

       i = 0
       j = 0
       write(flabel,'(a7,i4,i1)') 'h12cAOa',1,i
       call x2c_read(flabel,scr1_mat,naosh_ls*naosh_ls*nz,x2c_file_amf)


!      --------------------------------------------------------
!      step 2. select LL-block of h1_++_saao_cartesian_gaussian
!      --------------------------------------------------------
!          in: scr1_mat
!         out: scr2_mat
       call pick_LL_block_saao_bas(scr1_mat,            &
                                   scr2_mat,            &
                                   naosh_L,             &
                                   naosh_all,           &
                                   naosh_ls,            &
                                   nr_ao_L,             &
                                   nfsym,               &
                                   nz,                  &
                                   0,                   &
                                   op_bs_to_fs)
!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr2_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - h1 wo aoosoc contribution',6)
       end if

!      ---------------------------------------------------------------------------------
!      step 3. read and distribute aoosoc contributions for symmetry independent centers
!      ---------------------------------------------------------------------------------
!        out: aoosoc  -> scr1_mat

       allocate(aooso_info(naosh_ls))

       center_total = 1
       scr1_mat     = 0.0d0

       do i = 1, nuclei_type

!        count the number of basis functions for center_total
         call labcount(ncols_frag,aooso_info,nr_ao_L,1,-1,center_total,1,-1)
         nrows_frag = ncols_frag

#ifdef MOD_AOOSOC
         call put_atomic_oo_order_so_correction(                &
                                                scr1_mat,       &
                                                nr_ao_L,        &
                                                nrows_frag,     &
                                                ncols_frag,     &
                                                4,              &
                                                i,              &
                                                center_total,   &
                                                nuclei_totnum,  &
                                                print_lvl       &
                                               )
#else
         call quit('call put_atomic_oo_order_so_correction inaccessible')
#endif
       end do

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr1_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - aooSOC contribution',6)
       end if

#ifdef DEBUG_SOC
       open(99,file='soc-contributions',status='replace',form='formatted',  &
       access='sequential',action="readwrite",position='rewind')

       call print_x2cmat(scr1_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - aooSOC contribution',99)
       close(99,status='keep')

       open(99,file='aoo2e-ss-soc-contributions',status='replace',form='unformatted',  &
       access='sequential',action="readwrite",position='rewind')

       write(99) scr1_mat(1:nr_ao_L*nr_ao_L*nz)
       close(99,status='keep')
#endif

       deallocate(aooso_info)

!      -----------------------------------------
!      step 4. add aoo-SO contributions to h1_2c
!      -----------------------------------------

       call daxpy(nr_ao_L*nr_ao_L*nz,1.0d0,scr1_mat,1,scr2_mat,1)

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr2_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - h1 + aooSOC contribution',6)
       end if

!      -------------------------------------------------------------
!      step 5. write h1_++_LL_saao_cartesian_gaussian to file X2CMAT 
!      -------------------------------------------------------------

       i = 0
       j = 0
       write(flabel,'(a7,i4,i1)') 'h12cAOA',1,i


       if(spherical_on)then

         lwrk = is_len_wrk
         allocate(wrk(lwrk))
         call sph_iotc(scr2_mat,scr1_mat,nr_sph,print_lvl,wrk,lwrk)
         deallocate(wrk)

         call x2c_write(flabel,scr1_mat,nr_sph**2  * nz,x2c_file_glb)
       else

!        debug print
         if(print_lvl > 2)then
           call print_x2cmat(scr2_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - h12c+ooSO',6)
         end if

         call x2c_write(flabel,scr2_mat,nr_ao_L**2 * nz,x2c_file_glb)
       end if

!      release scratch space
       deallocate(scr1_mat)
       deallocate(scr2_mat)

  end subroutine x2c_interface2aoosoc

!**********************************************************************

end module x2c_2e_soc_interface
