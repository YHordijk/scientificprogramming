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
! module containing the functionality to obtain the exact two-component
! one-electron Hamiltonian h1_infinite_2c in orthonormal restricted-kinetic 
! balance basis.
!
! written by sknecht july 2010
!            sknecht august 2012 - added spinfree Hamiltonian
!
module x2c_get_h1inf2c_mo

  use x2c_fio
  use x2c_utils,        only:         &
      print_x2cmat,                   &
      make_spinfree_h1_onmo,          &
      get_boson_irrep_info
  use x2c_decoupling_mat

  implicit none

  public get_h1_2c_onmo_bas

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0

contains

!**********************************************************************
  subroutine get_h1_2c_onmo_bas(xmat,                     &
                                scr1mat,                  &
                                scr2mat,                  &
                                eigvl,                    &
                                norb_dim,                 &  
                                nesh_dim,                 &
                                npsh_dim,                 &
                                nfsym,                    &
                                nz,                       &
                                ipq_off,                  &
                                x2c_file_scr,             &
                                defining_h1mat,           &
                                spinfree,                 &
                                print_lvl)
!**********************************************************************
!
!    purpose: driver routine for obtaining the two-component infinite
!    order one-electron Hamiltonian h1_infinite_2c in 
!    orthonormal restricted-kinetic balance basis
!
!    the bare-nucleus Hamiltonian h1 is returned (in xmat) in its infinite-order two-component form
!----------------------------------------------------------------------
     real(8), intent(inout) :: xmat(*)           ! IN: defining 4c-h1 , OUT: bare-nucleus 2c-h1, (orthonormal basis)
     real(8), intent(inout) :: scr1mat(*)
     real(8), intent(inout) :: scr2mat(*)
     real(8), intent(inout) :: eigvl(*)
     integer, intent(in)    :: nz
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: nesh_dim(nfsym)
     integer, intent(in)    :: npsh_dim(nfsym)
     integer, intent(in)    :: ipq_off(4,0:7)
     integer, intent(in)    :: x2c_file_scr
     integer, intent(in)    :: defining_h1mat
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: spinfree
!----------------------------------------------------------------------
     integer                :: ioff_osh1
     integer                :: ioff_oeig
     integer                :: ioff_eeig
     integer                :: ioff_esh1
     integer                :: ioff_esh2
     integer                :: norb1_f
     integer                :: n1esh_f
     integer                :: n1psh_f
     integer                :: n2esh_f
     integer                :: nesh2_f
     integer                :: ierr_qdiag
     integer                :: i, j, ix
     character (len=12)     :: flabel
     character (len= 5)     :: flabel_tag
     integer, allocatable   :: boson_info(:)
!**********************************************************************
!

#ifdef spinfree_after_approach
       if(spinfree)then
         allocate(boson_info(norb_dim(1)+norb_dim(2)))
         call get_boson_irrep_info(boson_info,  &
                                   nfsym,       &
                                   norb_dim,    &
                                   print_lvl)
       end if
#endif

!      initialize pointers for matrices and vectors
       ioff_osh1 = 1
       ioff_oeig = 1
       ioff_eeig = 1
       ioff_esh1 = 1
       ioff_esh2 = 1
!
!      obtain h1_infinite_2c for each fermion ircop
       do i = 1, nfsym

!        set dimensions for matrices and vectors used in the subroutines
         norb1_f  = norb_dim(i)
         n1psh_f  = npsh_dim(i)
         n1esh_f  = nesh_dim(i)
         n2esh_f  = nesh_dim(i)**2
         nesh2_f  = 2 * nesh_dim(i)
!        print *,'norb1_f, n1esh_f',norb1_f, n1esh_f

         if( norb_dim(i) > 0 )then

!          set unique last part of the file ID
           write(flabel_tag,'(i4,i1)') 1,i

!          step 1: construct the decoupling matrix R in orthonormal MO basis
!                  this step includes:
!                  ==> diagonalization of the defining h1_4c matrix
!                  ==> setup of the "A" and "B" equations
!                  ==> solving the "A" and "B" equations to obtain the decoupling "R" matrix
!          ---------------------------------------------------------------------------------
!          in : defining 4c-h1 -> scr1mat
!          out: R -> scr1mat
!               on file under the label Rmat_ON

           call construct_r_onmo_bas(xmat,           &
                                     scr1mat,        &
                                     scr2mat,        &
                                     eigvl,          &
                                     ioff_osh1,      &
                                     ioff_esh1,      &
                                     ioff_oeig,      &
                                     ioff_eeig,      &
                                     norb1_f,        &
                                     n1psh_f,        &
                                     n1esh_f,        &
                                     nz,             &
                                     ipq_off(1,0),   &
                                     x2c_file_scr,   &
                                     flabel_tag,     &
                                     print_lvl)

!          step 2: a. calculate the unitary transformation matrix U = w1.w2:
!                     - also known as picture change transformation matrix (in orthonormal MO basis)
!                     - stored on file under the label Umat_ON
!
!                  in : R matrix in orthonormal basis --> scr1mat
!                  out: U -> scr2mat
!          -----------------------------------------------------------------
           call construct_pctmat_onmo_bas(xmat,           &
                                          scr1mat,        &
                                          scr2mat,        &
                                          eigvl,          &
                                          ioff_osh1,      &
                                          ioff_esh1,      &
                                          ioff_esh2,      &
                                          ioff_eeig,      &
                                          norb1_f,        &
                                          n1psh_f,        &
                                          n1esh_f,        &
                                          nesh2_f,        &
                                          nz,             &
                                          ipq_off(1,0),   &
                                          x2c_file_scr,   &
                                          flabel_tag,     &
                                          print_lvl)

!          debug print
           if(print_lvl > 2)then
             call print_x2cmat(scr2mat(ioff_esh2),nesh2_f,n1esh_f,nz,ipq_off(1,0),'x2c - umat_onmo',6)
           end if

!          step 3: perform U+ h1_4c_ON U = h1_infinite_2c_ON
!          -------------------------------------------------

!          read op_4c
           write(flabel,'(a7,a5)') 'h1_4cON',flabel_tag
           call x2c_read(flabel,scr1mat(ioff_osh1),norb1_f**2 * nz,x2c_file_scr)

!          in : h1_4c_ON -> scr1mat; U  -> scr2mat 
!          out: h1_infinite_2c_ON -> xmat
           call pctrafo_op_4cto2c_onmo_bas(xmat(ioff_esh1),    &
                                           scr1mat(ioff_osh1), &
                                           scr2mat(ioff_esh2), &
                                           n1esh_f,            &
                                           nesh2_f,            &
                                           nz,                 &
                                           ipq_off(1,0),       &
                                           print_lvl)
            
#ifdef spinfree_after_approach
           if(spinfree)then
             call make_spinfree_h1_onmo(                               &
                                        xmat(ioff_esh1),               & 
                                        n1esh_f,                       & 
                                        n1esh_f,                       & 
                                        n1esh_f,                       & 
                                        n1esh_f,                       & 
                                        boson_info(ioff_oeig+n1psh_f), &
                                        boson_info(ioff_oeig+n1psh_f), &
                                        nz,                            &
                                        print_lvl                      &
                                       )
           end if
#endif

!          debug print
           if(print_lvl > 1)then
             call print_x2cmat(xmat(ioff_esh1),n1esh_f,n1esh_f,nz,ipq_off(1,0),'x2c - h1_inf_2c',6)
           end if
          
         end if

!        update offsets for matrices and vectors
         ioff_osh1 = ioff_osh1 + norb1_f**2 * nz
         ioff_esh1 = ioff_esh1 + n1esh_f**2 * nz
         ioff_esh2 = ioff_esh2 + n1esh_f**2 * nz * 2

         ioff_oeig = ioff_oeig + norb1_f
         ioff_eeig = ioff_eeig + n1esh_f

       end do

#ifdef spinfree_after_approach
       if(spinfree) deallocate(boson_info)
#endif

  end subroutine get_h1_2c_onmo_bas

!**********************************************************************
  subroutine pctrafo_op_4cto2c_onmo_bas(op_inf_2c,               &
                                        op_4c,                   &
                                        umat,                    &
                                        ndim_o2c,                &
                                        ndim_o4c,                &
                                        nz,                      &
                                        ipq_off_in,              &
                                        print_lvl)
!**********************************************************************
!
!    purpose: perform the decoupling U+ * 4c-op * U = op_infinite_2c^++ (orthonormal MO basis) 
!             by means of the unitary transformation matrix U
!
!    equation (9) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
     real(8), intent(inout) :: op_inf_2c(*)
     real(8), intent(inout) :: op_4c(*)
     real(8), intent(inout) :: umat(*)
     integer, intent(in)    :: ndim_o2c
     integer, intent(in)    :: ndim_o4c
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(*)
     integer, intent(in)    :: print_lvl
!----------------------------------------------------------------------
!**********************************************************************

!      op_infinite_2c^++ = U+ op U
       call qtrans90('AOMO','S',val_d0,        &
                     ndim_o4c,                 & 
                     ndim_o4c,                 &
                     ndim_o2c,                 &
                     ndim_o2c,                 &
                     op_4c,                    &
                     ndim_o4c,                 &
                     ndim_o4c,                 &
                     nz,                       &
                     ipq_off_in,               &
                     op_inf_2c,                &
                     ndim_o2c,                 &
                     ndim_o2c,                 &
                     nz,                       &
                     ipq_off_in,               &
                     umat,                     &
                     ndim_o4c,                 &
                     ndim_o2c,                 &
                     nz,                       &
                     ipq_off_in,               &
                     umat,                     &
                     ndim_o4c,                 &
                     ndim_o2c,                 &
                     nz,                       &
                     ipq_off_in,               &
                     print_lvl)

  end subroutine pctrafo_op_4cto2c_onmo_bas

!**********************************************************************

end module x2c_get_h1inf2c_mo
