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
! module containing the selection of the basic "4c-h1" matrix which 
! will shall later be transformed to "2c-h1" using the pct-matrix "U" obtained 
! on the basis of the "defining 4c-h1".
!
! written by sknecht march 2011
!
module x2c_h1mat_base

  use x2c_fio
  use x2c_cb_interface, only:         &
      reset_x2c_cb_onefck
  use x2c_utils, only:                &
      print_x2cmat,                   &
      make_spinfree_h1_onmo,          &
      get_boson_irrep_info
  use x2c_def_h1_4c, only:            &
      get_h1_4c_saao_basis,           &
      get_h1_4c_on_basis

  implicit none

  public x2c_get_h1_mat_base

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0

contains

!**********************************************************************
  subroutine x2c_get_h1_mat_base(h1_onmo_bas,          &
                                 h1_saao_bas,          &
                                 tra_mat_4cAOMO,       &
                                 naosh_ls,             &
                                 naosh_all,            &
                                 norb_dim,             &
                                 ioff_aomx,            &
                                 nfsym,                &
                                 nz,                   &
                                 nzt,                  &
                                 ipq_off_in,           &
                                 basic_h1mat,          &
                                 is_len_wrk,           &
                                 x2c_file_scr,         &
                                 linear_sym,           &
                                 mdirac,               &
                                 spinfree,             &
                                 print_lvl)
!**********************************************************************
!
!    purpose: driver routine to obtain the basic h1 matrix in orthonormal MO basis
!             which will later be transformed from 4c --> 2c 
!             using the pct-matrix "U" obtained for the defining h1 matrix.
!
!             note: the basic h1 may or may not be identical with the defining h1 matrix.

!             currently implemented: 
!               - basic_h1mat == 0: case a. bare-nucleus 1-el Dirac-Hamiltonian
!               - basic_h1mat == 2: case c. 4c-Fock-Dirac operator matrix
!
!----------------------------------------------------------------------
!
     real(8), intent(inout) :: h1_onmo_bas(*)
     real(8), intent(inout) :: h1_saao_bas(*)
     real(8), intent(inout) :: tra_mat_4cAOMO(*)
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: nzt
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: basic_h1mat
     integer, intent(in)    :: is_len_wrk
     integer, intent(in)    :: x2c_file_scr
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: linear_sym
     logical, intent(in)    :: mdirac
     logical, intent(in)    :: spinfree
!----------------------------------------------------------------------
     integer                :: ioff_osh1
     integer                :: ioff_osh2
     integer                :: ioff_ao_o_sh1   
     integer                :: ioff_ao_ao_sh1
     integer                :: i       
     integer                :: norb1_f       
     integer                :: nbas1_f
     integer                :: ndim_ao 
     integer                :: mz
     character (len=12)     :: flabel
     integer, allocatable   :: boson_info(:)
!**********************************************************************

!      step 1: get the basic 4c-h1 matrix in SA-AO basis according to the value for basic_h1mat
!
!              in : basic_h1mat as integer value: 
!                   0 == case a.; bare-nucleus 4c-h1 
!                   2 == case c.; 4c-Dirac-Fock matrix
!              out: 4c-h1 matrix in SA-AO basis --> h1_saao_bas
!      ------------------------------------------------------------------------------
       call get_h1_4c_saao_basis(h1_saao_bas,            &
                                 h1_onmo_bas,            &
                                 tra_mat_4cAOMO,         &
                                 naosh_ls,               &
                                 naosh_all,              &
                                 norb_dim,               &
                                 nfsym,                  &
                                 nz,                     &
                                 nzt,                    &
                                 ipq_off_in,             &
                                 basic_h1mat,            &
                                 linear_sym,             &
                                 mdirac,                 &
                                 is_len_wrk,             &
                                 x2c_file_scr,           &
                                 print_lvl)

!      step 2: transform the basic 4c-h1 to orthonormal MO basis
!
!              in : 4c h1 matrix in SA-AO basis --> h1_saao_bas; transformation matrix --> tra_mat_4cAOMO
!              out: h1 matrix in 4c-ON MO basis --> h1_onmo_bas
!      ------------------------------------------------------------------------------

       mz = nz
       if (nfsym == 2) then
         allocate(boson_info(norb_dim(1)+norb_dim(2)))
       else
         allocate(boson_info(norb_dim(1)))
       endif
       if(spinfree)then
         call get_boson_irrep_info(boson_info,  &
                                   nfsym,       &
                                   norb_dim,    &
                                   print_lvl)
         mz = 1
       end if

!      initialize matrix offsets
       ioff_osh1      = 1
       ioff_osh2      = 1
       ioff_ao_o_sh1  = 1
       ioff_ao_ao_sh1 = 1

       do i = 1, nfsym
 
!        set dimensions and special offsets for matrices accessed in the
!        subroutine
         nbas1_f = naosh_all(i)
         norb1_f = norb_dim(i)
         ndim_ao = naosh_all(i)

!        the h1 matrices (1el-dirac or free-particle matrix) have a different storage mode
!        in memory in saao basis than the 4c-fock operator
         if(basic_h1mat /= 2)then
           ioff_ao_ao_sh1 = ioff_aomx(i,i) + 1
           ndim_ao        = naosh_ls
         end if

         if( norb_dim(i) > 0 )then

!          debug print
           if (print_lvl >= 3) then
             call print_x2cmat(h1_saao_bas(ioff_ao_ao_sh1),nbas1_f,nbas1_f,nz,ipq_off_in(1,0),      &
                              'x2c - h1_4c-basic-saao',6)
           end if

           call get_h1_4c_on_basis(h1_onmo_bas(ioff_osh1),            &
                                   h1_saao_bas(ioff_ao_ao_sh1),       &
                                   tra_mat_4cAOMO(ioff_ao_o_sh1),     &
                                   ndim_ao,                           &
                                   nbas1_f,                           &
                                   norb1_f,                           &
                                   nz,                                &
                                   mz,                                &
                                   nzt,                               &
                                   ipq_off_in,                        &
                                   spinfree,                          &
                                   boson_info(ioff_osh2),             &
                                   print_lvl)
         end if

!        store the basic 4c-h1 in orthonormal basis on file (label == h1_4cON)
         write(flabel,'(a7,i4,i1)') 'h1_4cON',1,i
         call x2c_write(flabel,h1_onmo_bas(ioff_osh1),norb1_f**2 * nz,x2c_file_scr)

!        update offsets for matrices
         ioff_ao_o_sh1  = ioff_ao_o_sh1  + nbas1_f * norb1_f * nzt
         ioff_ao_ao_sh1 = ioff_ao_ao_sh1 + nbas1_f * nbas1_f * nz
         ioff_osh1      = ioff_osh1      + norb1_f**2 * nz
         ioff_osh2      = ioff_osh2      + norb1_f

       end do

       deallocate(boson_info)

  end subroutine x2c_get_h1_mat_base

!**********************************************************************

end module x2c_h1mat_base
