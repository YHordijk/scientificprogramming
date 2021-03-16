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
module x2c_fragment_x2c_interface

! stefan: this module contains all the functionality to carry out an ensuing 
!         transformation
!         ========================
!         h_2c^++ = U^+ * h_4c * U
!         ========================
!         where U was obtained by means of an atom-by-atom X2C approach:
!         U =~ SUM_{all unique atomic center/fragments i} U_i = U^{@}
!
!         for more details see Peng and Reiher, JCP 2012
!
! written by sknecht jan 2013
!

  use fragment_x2c_cfg
  use x2c_fio
  use x2c_pct_ao, only: &
      pctrafo_h1_4cto2c_saao_bas

  implicit none

  public x2c_fragment_x2c_driver
  public x2c_fragment_get_h1_2c

  private


contains

!***********************************************************************
  subroutine x2c_fragment_x2c_driver(nz,                  &
                                     nfsym,               &
                                     mxatom,              &
                                     type_nuclei,         &
                                     type_charge,         &
                                     nr_symm_indep_cent,  &
                                     nr_degen_nuc_cent,   &
                                     nr_ao_bas_type,      &
                                     naosh_L_mol,         &
                                     naosh_all_mol,       &
                                     funit                &
                                    )

!**********************************************************************
!
!    purpose: read pct-matrices U_fragment for atomic types and
!             generate the approximate molecular U
!
!----------------------------------------------------------------------
  integer, intent(in)  :: nz
  integer, intent(in)  :: nfsym
  integer, intent(in)  :: mxatom
  integer, intent(in)  :: naosh_L_mol(nfsym)
  integer, intent(in)  :: naosh_all_mol(nfsym)
  integer, intent(in)  :: nr_ao_bas_type(mxatom,0:2)
  integer, intent(in)  :: type_nuclei
  integer, intent(in)  :: type_charge(*)
  integer, intent(in)  :: nr_symm_indep_cent(*)
  integer, intent(in)  :: nr_degen_nuc_cent(*)
  integer, intent(in)  :: funit
!----------------------------------------------------------------------
  real(8), allocatable :: umol(:)
  integer              :: i
  integer              :: j
  integer, parameter   :: nzc1    = 4
  integer, parameter   :: nfsymc1 = 1
  character(len=12)    :: flabel
!----------------------------------------------------------------------

!     safety check
      if(.not.fragment_x2c_info%fragment_approach_enabled) return

      if(nz /= nzc1)then
        stop ' atomic X2C approach works so far only in C1 symmetry'
      end if

      if(nfsym /= nfsymc1)then
        stop ' atomic X2C approach works so far only in C1 symmetry'
      end if

!     allocate molecular pct-matrix
      j    = 0
      do i = 1, nfsymc1
        j = j + naosh_all_mol(i)*naosh_L_mol(i)*nzc1
      end do
      allocate(umol(j))
      umol = 0.0d0

!     default for the following procedure is to use the atomic fragments
!     ----------------------------------------------------------------------------
!     possible extensions are: unit matrices for light atoms and fragment matrices 
!     for heavy centers being close or even sharing a chemical bond.
!     ----------------------------------------------------------------------------

      call construct_umol_from_atoms(umol,                &
                                     nzc1,                &
                                     nfsymc1,             &
                                     nz,                  &
                                     nfsym,               &
                                     mxatom,              &
                                     type_nuclei,         &
                                     type_charge,         &
                                     nr_symm_indep_cent,  &
                                     nr_degen_nuc_cent,   &
                                     nr_ao_bas_type(1,0), &
                                     naosh_L_mol,         &
                                     naosh_all_mol        &
                                    )

!     put molecular pctmat to file
      write(flabel,'(a7,i4,i1)') 'pctmtAO',1,1

      call x2c_write(flabel,                 &
                     umol,                   &
                     naosh_all_mol(1) *      &
                     naosh_L_mol(1)   *      &
                     nzc1,                   &
                     funit)

      deallocate(umol)

  end subroutine x2c_fragment_x2c_driver
!***********************************************************************

  subroutine construct_umol_from_atoms(umol,                &
                                       nzc1,                &
                                       nfsymc1,             &
                                       nz,                  &
                                       nfsym,               &
                                       mxatom,              &
                                       type_nuclei,         &
                                       type_charge,         &
                                       nr_symm_indep_cent,  &
                                       nr_degen_nuc_cent,   &
                                       nr_ao_bas_type,      &
                                       naosh_L_mol,         &
                                       naosh_all_mol        &
                                      )

!**********************************************************************
!
!    purpose: read pct-matrices U_fragment for atomic types and
!             generate the approximate molecular U
!
!----------------------------------------------------------------------
  real(8), intent(inout) :: umol(*)
  integer, intent(in)    :: nzc1
  integer, intent(in)    :: nfsymc1
  integer, intent(in)    :: nz
  integer, intent(in)    :: nfsym
  integer, intent(in)    :: mxatom
  integer, intent(in)    :: type_nuclei
  integer, intent(in)    :: type_charge(*)
  integer, intent(in)    :: nr_symm_indep_cent(*)
  integer, intent(in)    :: nr_degen_nuc_cent(*)
  integer, intent(in)    :: nr_ao_bas_type(mxatom,0:2)
  integer, intent(in)    :: naosh_all_mol(nfsymc1)
  integer, intent(in)    :: naosh_L_mol(nfsymc1)
!----------------------------------------------------------------------
  integer, allocatable   :: rbuf(:)
  integer, allocatable   :: cbuf(:)
  integer, allocatable   :: ibuf(:)
  integer, allocatable   :: jbuf(:)
  integer                :: icent
  integer                :: iatom
!----------------------------------------------------------------------

!     debug option - read U matrix from molecular calculation
      if(fragment_x2c_info%fragment_approach_ismolecule)then
        call fragment_x2c_init(fragment_x2c_info,          &
                               1,                          &
                               nzc1,                       &
                               0,                          &
                               naosh_all_mol(1),           &
                               naosh_L_mol(1)              &
                              )
        call dcopy(naosh_all_mol(1)*naosh_L_mol(1)*nzc1,   &
                   fragment_x2c_info%pctmat,               &
                   1,                                      & 
                   umol,                                   & 
                   1                                       & 
                   )
        return
      end if

      icent    = 1
      iatom    = 1

      allocate(rbuf(naosh_all_mol(1)))
      allocate(cbuf(naosh_L_mol(1)))

      do ! loop over atom types

        if(iatom > type_nuclei) exit

!       set # of basis functions (L and L+S), allocate pctmat-frag 
!       for this atom type and read atomic pctmat from file
        call fragment_x2c_init(fragment_x2c_info,          &
                               1,                          &
                               nzc1,                       &
                               type_charge(iatom),         &
                               nr_ao_bas_type(iatom,0),    &
                               nr_ao_bas_type(iatom,1)     &
                              )

         allocate(ibuf(fragment_x2c_info%naosh_all(1)))
         allocate(jbuf(fragment_x2c_info%naosh_L(1)))

         rbuf = 0
         cbuf = 0
         ibuf = 0
         jbuf = 0
        
         call put_aoblock(umol,                            &
                          naosh_all_mol(1),                &
                          naosh_L_mol(1),                  &
                          icent,                           &
                          nr_symm_indep_cent(iatom),       &
                          nr_degen_nuc_cent(icent),        &
                          nzc1,                            &
                          fragment_x2c_info%pctmat,        &
                          fragment_x2c_info%naosh_all(1),  &
                          fragment_x2c_info%naosh_all(1),  &
                          -1,                              &
                          fragment_x2c_info%naosh_L(1),    &
                          fragment_x2c_info%naosh_L(1),    &
                          -1,                              &
                          rbuf,                            &
                          cbuf,                            &
                          ibuf,                            &
                          jbuf)

        deallocate(ibuf)
        deallocate(jbuf)

!       update indices
        icent = icent + nr_symm_indep_cent(iatom)
        iatom = iatom + 1
      end do

      deallocate(rbuf)
      deallocate(cbuf)

  end subroutine construct_umol_from_atoms
!***********************************************************************

  subroutine x2c_fragment_get_h1_2c(naosh_ls,             &
                                    nr_ao_L,              &
                                    naosh_all,            &
                                    naosh_L,              &
                                    nfsym,                &
                                    nz,                   &
                                    ipq_off_in,           &
                                    ioff_aomx,            &
                                    add_amfi_contrib,     &
                                    is_final_ham_lvl,     &
                                    x2c_file_glb,         &
                                    x2c_file_amf,         &
                                    op_bs_to_fs,          &
                                    is_len_wrk,           &
                                    spherical_on,         &
                                    mdirac,               &
                                    print_lvl)
!**********************************************************************
!
!    purpose: picture-change transform the bare-nucleus 1e Dirac 
!             Hamiltonian and store it on file.
!
!             The bare-nucleus corrected Hamiltonian is used in the initial 
!             SCF step whereas the default Hamiltonian (including AMFI contributions) 
!             will be considered in all subsequent SCF steps. 
!
!             equation: h_{2c_D}^++ = U+ * h_{4c_D} * U
!
!----------------------------------------------------------------------
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: nr_ao_L
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: naosh_L(nfsym)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(inout) :: is_final_ham_lvl
     integer, intent(in)    :: add_amfi_contrib
     integer, intent(in)    :: is_len_wrk
     integer, intent(in)    :: x2c_file_glb
     integer, intent(in)    :: x2c_file_amf
     integer, intent(in)    :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: spherical_on
     logical, intent(in)    :: mdirac
!----------------------------------------------------------------------
     real(8), allocatable   :: h1_2cinf_saao(:)
     integer                :: funit
     integer                :: nr_aosh_h1
     logical                :: get_ll_block_2c
     logical                :: dobncorr
     character (len=12)     :: flabel
!----------------------------------------------------------------------

      dobncorr = .false.

!     construct h12cAOn or h12cAOa depending on amfi flag
      if(add_amfi_contrib > 0)then

!       set pointer to Hamiltonian integral type - final adaption after AMFI step
        is_final_ham_lvl = 5
        write(flabel,'(a7,i4,i1)') 'h12cAOa',1,0
        funit           = x2c_file_amf
        get_ll_block_2c = .false.

      else

!       set pointer to Hamiltonian integral type 
        is_final_ham_lvl = 2
        write(flabel,'(a7,i4,i1)') 'h12cAOn',1,0
        funit           = x2c_file_glb
        get_ll_block_2c = .true.

      end if ! amfi contributions 

      allocate(h1_2cinf_saao(naosh_ls**2*nz))
      h1_2cinf_saao = 0

      call pctrafo_h1_4cto2c_saao_bas(h1_2cinf_saao,        &
                                      naosh_ls,             &
                                      nr_ao_L,              &
                                      naosh_all,            &
                                      naosh_L,              &
                                      nr_aosh_h1,           &
                                      nfsym,                &
                                      nz,                   &
                                      ipq_off_in,           &
                                      ioff_aomx,            &
                                      x2c_file_glb,         &
                                      op_bs_to_fs,          &
                                      is_len_wrk,           &
                                      dobncorr,             &
                                      get_ll_block_2c,      &
                                      spherical_on,         &
                                      mdirac,               &
                                      flabel,               &
                                      print_lvl)

      call x2c_write(flabel,                                &
                     h1_2cinf_saao,                         &
                     nr_aosh_h1*nr_aosh_h1*nz,              &
                     funit                                  &
                    ) 

      deallocate(h1_2cinf_saao)

  end subroutine x2c_fragment_get_h1_2c
!***********************************************************************

end module x2c_fragment_x2c_interface
