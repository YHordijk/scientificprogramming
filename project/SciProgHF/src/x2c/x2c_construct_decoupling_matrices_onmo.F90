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
! module containing the functionality for obtaining the decoupling 
! matrix (positive/negative energy solutions) in the X2C formalism
! for any defining 4c-h1 matrix in orthonormal basis
!
! written by sknecht april 2010
!
module x2c_decoupling_mat

  use x2c_fio
  use x2c_utils

  implicit none

  public construct_pctmat_onmo_bas
  public construct_r_onmo_bas

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0
  real(8), parameter :: val_thrnull =  1.0d-10

contains

!**********************************************************************
  subroutine construct_pctmat_onmo_bas(xmat,           &
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
                                       ipq_off_in,     &
                                       funit,          &
                                       mat_id,         &
                                       print_lvl)
!**********************************************************************
!
!    purpose: driver routine for obtaining the unitary transformation matrix 
!    U = w1[R].w2[R] in orthonormal MO basis starting from the decoupling matrix R.

!    (equation (10) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)):
!
!        [ h11 h12 ]     [ h+   0  ]
!    U^+ [         ] U = [         ]
!        [ h21 h22 ]     [  0  h-  ]
!
!    equation (9) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
!    where
!
!            [  1  -R^+ ]                [ N_+^-1    0  ]
!    w1[R] = [          ]   and  w2[R] = [              ]
!            [  R    1  ]                [    0   N_-^-1]
!
!    with N_+^-1 = sqrt(1+R^+.R) and N_-^-1 = sqrt(1+R.R^+)
!
!    see equation (10) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: xmat(*)
     real(8), intent(inout) :: scr1mat(*)
     real(8), intent(inout) :: scr2mat(*)
     real(8), intent(inout) :: eigvl(*)
!    offsets for matrices and vectors
     integer, intent(in)    :: ioff_osh1
     integer, intent(in)    :: ioff_esh1
     integer, intent(in)    :: ioff_esh2
     integer, intent(in)    :: ioff_eeig
!    dimensions for matrices and vectors
     integer, intent(in)    :: norb1_f
     integer, intent(in)    :: n1psh_f
     integer, intent(in)    :: n1esh_f
     integer, intent(in)    :: nesh2_f
!    info about quaternion symmetry
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(*)
!    unique matrix/array ID on file funit
     integer, intent(in)    :: funit
     character (len=5)      :: mat_id
!    debug print level
     integer, intent(in)    :: print_lvl
!----------------------------------------------------------------------
     character (len=12)     :: flabel
!**********************************************************************
!

!      step 1: construct w2 and w1 using the R matrix
!        step 1a: construct w2 
!           in : R -> scr1mat
!           out: 1+R^+.R -> xmat
         call construct_1_plus_r_r(xmat(ioff_esh1),       &
                                   scr1mat(ioff_esh1),    &
                                   n1esh_f,               &
                                   nz,                    &
                                   ipq_off_in)

!           in : 1+R^+.R -> xmat
!           out: X^(-1/2) -> xmat; X -> scr2mat
         call diag_sqrt_1_plus_r_r(xmat(ioff_esh1),       &
                                   scr2mat(ioff_esh1),    &
                                   eigvl(ioff_eeig),      &
                                   n1esh_f,               &
                                   nz)
       
!           in : X^(-1/2) -> xmat; X -> scr2mat
!           out: w2 -> scr1mat
         call generate_w2(scr1mat(ioff_esh1),             &
                          xmat(ioff_esh1),                &
                          scr2mat(ioff_esh1),             &
                          n1esh_f,                        &
                          nz,                             &
                          ipq_off_in)
! For testing purposes you can replace the above w2 matrix with the unit matrix here:
!          call dunit2(scr1mat(ioff_esh1),                  &
!                      n1esh_f,                             &
!                      n1esh_f,n1esh_f,nz)

 
!        INCLUDE HERE A READPOINT FOR the 'R' matrix?
         write(flabel,'(a7,a5)') 'Rmat_ON',mat_id
         call x2c_read(flabel,scr2mat(ioff_esh1),n1esh_f**2 * nz,funit)

!        step 1b: w1
!           in : R  -> scr2mat
!           out: w1 -> xmat
         call generate_w1(xmat(ioff_esh2),                &
                          scr2mat(ioff_esh1),             &
                          n1esh_f,                        &
                          nesh2_f,                        &
                          nz)

!      step 2: obtain the matrix U
!           in : w1 -> xmat; w2 -> scr1mat
!           out: U -> scr2mat
       call generate_u_mat(scr2mat(ioff_esh2),            &
                           xmat(ioff_esh2),               &
                           scr1mat(ioff_esh1),            &
                           nesh2_f,                       &
                           n1esh_f,                       &
                           nz,                            &
                           ipq_off_in)

!      INCLUDE HERE A SAVEPOINT FOR the 'U' matrix?
       write(flabel,'(a7,a5)') 'Umat_ON',mat_id
       call x2c_write(flabel,scr2mat(ioff_esh2),n1esh_f**2 * nz * 2,funit)

  end subroutine construct_pctmat_onmo_bas

!**********************************************************************
  subroutine construct_r_onmo_bas(xmat,           &
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
                                  ipq_off_in,     &
                                  funit,          &
                                  mat_id,         &
                                  print_lvl)
!**********************************************************************
!
!    purpose: obtain the R matrix (decoupling matrix) in orthonormal MO basis 
!    which will then be used to construct the picture-change transformation matrix 
!    U = w1[R].w2[R] (equation (10) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)):
!
!        [ h11 h12 ]     [ h+   0  ]
!    U^+ [         ] U = [         ]
!        [ h21 h22 ]     [  0  h-  ]
!
!    equation (9) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
!    where 
!
!            [  1  -R^+ ]                [ N_+^-1    0  ]
!    w1[R] = [          ]   and  w2[R] = [              ]
!            [  R    1  ]                [    0   N_-^-1]
!
!    with N_+^-1 = sqrt(1+R^+.R) and N_-^-1 = sqrt(1+R.R^+)
!
!    see equation (10) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
!    The R matrix is calculated by solving the "R-equation": A.R = B
!
!    see equations (23+24) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: xmat(*)
     real(8), intent(inout) :: scr1mat(*)
     real(8), intent(inout) :: scr2mat(*)
     real(8), intent(inout) :: eigvl(*)
!    offsets for matrices and vectors
     integer, intent(in)    :: ioff_osh1
     integer, intent(in)    :: ioff_esh1
     integer, intent(in)    :: ioff_eeig
     integer, intent(in)    :: ioff_oeig
!    dimensions for matrices and vectors
     integer, intent(in)    :: norb1_f
     integer, intent(in)    :: n1psh_f
     integer, intent(in)    :: n1esh_f
!    info about quaternion symmetry
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(*)
!    unique matrix/array ID on file funit
     integer, intent(in)    :: funit
     character (len=5)      :: mat_id
!    debug print level
     integer, intent(in)    :: print_lvl
!----------------------------------------------------------------------
     character (len=12)     :: flabel
     integer                :: ierr_qdiag
!**********************************************************************

!      step 1: diagonalize the defining h1_4c in orthonormal basis
!      -----------------------------------------------------------
       ierr_qdiag = 0

!      debug print
       if(print_lvl > 1)then
         call print_x2cmat(xmat(ioff_osh1),norb1_f,norb1_f,nz,ipq_off_in,'x2c - defining-h1_4c',6)
       end if

       call qdiag90(nz,norb1_f,xmat(ioff_osh1),norb1_f,norb1_f,            &
                    eigvl(ioff_oeig),1,scr1mat(ioff_osh1),norb1_f,norb1_f, &
                    ierr_qdiag)

       if( ierr_qdiag /= 0 )                                               &
         call quit(' *** error in get_decoupling_mat: qdiag returned with error code. ***')

!      step 2: set up the 'A' and 'B' equations
!         in : eigenvectors of the defining 4c-h1 -> scr1mat
!         out: A -> xmat; B -> scr2mat
       call setup_a_b(xmat(ioff_esh1),                    &
                      scr2mat(ioff_esh1),                 &
                      scr1mat(ioff_osh1),                 & 
                      n1esh_f,                            &
                      norb1_f,                            &
                      n1psh_f,                            &
                      nz,                                 &
                      ipq_off_in)


!      step 3: solve the R equation "AR = B"
!         in : A -> xmat; B -> scr2mat 
!         out: R -> scr1mat
       call solve_r_eq(scr1mat(ioff_esh1),                &
                       xmat(ioff_esh1),                   &
                       scr2mat(ioff_esh1),                &
                       eigvl(ioff_eeig),                  &
                       n1esh_f,                           &
                       nz,                                &
                       ipq_off_in,                        &
                       print_lvl)

!
!      save the R matrix in orthonormal MO basis on file
       write(flabel,'(a7,a5)') 'Rmat_ON',mat_id
       call x2c_write(flabel,scr1mat(ioff_esh1),n1esh_f**2 * nz,funit)

  end subroutine construct_r_onmo_bas

!**********************************************************************
  subroutine setup_a_b(amat,                      &
                       bmat,                      &
                       scrmat,                    &    
                       ndim_ab,                   &
                       ndim_sc,                   &
                       add_off_sc,                &
                       nz,                        &
                       ipq_off_in)
!**********************************************************************
!
! 
!    purpose: set up the 'A' and 'B' equations
!
!    A :=   Y--.Y--^+
!    B := - Y--.Y+-^+ 
!
!    equation (24) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
     real(8), intent(inout) :: amat(*)
     real(8), intent(inout) :: bmat(*)
     real(8), intent(inout) :: scrmat(*)
     integer, intent(in)    :: ndim_ab
     integer, intent(in)    :: ndim_sc
     integer, intent(in)    :: add_off_sc
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(*)
!**********************************************************************
!
!      A := Y--.Y--^+
       call qgemm(ndim_ab,ndim_ab,ndim_ab,val_d1,              &
                  'N','N',ipq_off_in,                          &
                  scrmat,ndim_sc,ndim_sc,nz,                   &
                  'H','N',ipq_off_in,                          &
                  scrmat,ndim_sc,ndim_sc,nz,                   &
                  val_d0,ipq_off_in,amat,ndim_ab,ndim_ab,nz)

!      B := - Y--.Y+-^+
       call qgemm(ndim_ab,ndim_ab,ndim_ab,val_dm1,             &
                  'N','N',ipq_off_in,                          &
                  scrmat,ndim_sc,ndim_sc,nz,                   &
                  'H','N',ipq_off_in,                          &
                  scrmat(1+add_off_sc),ndim_sc,ndim_sc,nz,     &
                  val_d0,ipq_off_in,bmat,ndim_ab,ndim_ab,nz)

  end subroutine setup_a_b

!**********************************************************************
  subroutine solve_r_eq(rmat,         &
                        amat,         &
                        bmat,         &
                        h1eig,        &
                        ndim_ab,      &
                        nz,           &
                        ipq_off,      &
                        print_lvl)
!**********************************************************************
! 
!    purpose: solve the 'R' equation using Cholesky decomposition
! 
!    IDEA: replace the Cholesky part 2 (qchols) by two calls to 
!          dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
!          performing op(A)*X = alpha*B (for upper/lower triangular matrices A)
!          --> advantage: this would be a true BLAS routine. 
!              what remains unsolved is how to deal with the pivoting...
!
!    'R' equation: (Y--.Y--^+) R = -Y--.Y+-^+
!
!    equation (24) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
     real(8), intent(inout) :: rmat(*)
     real(8), intent(inout) :: amat(*)
     real(8), intent(inout) :: bmat(*)
     real(8), intent(inout) :: h1eig(*)
     integer, intent(in)    :: ndim_ab
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off(*)
     integer, intent(in)    :: print_lvl
!----------------------------------------------------------------------
     integer, allocatable   :: temp_mat1(:)
     real(8), allocatable   :: temp_mat2(:)
     integer                :: job_task
     integer                :: effective_dim
     integer                :: j
!**********************************************************************
       
       allocate(temp_mat1(ndim_ab))
       temp_mat1     = val_d0

       job_task      = 1
       effective_dim = 0
       call qchold(amat,ndim_ab,nz,ndim_ab,ndim_ab,h1eig,            &
                   1.0d-30,effective_dim,job_task,temp_mat1)

!      check for linear dependencies in Cholesky decomposition
       if( effective_dim < ndim_ab )then
 
          print '(/a,/a)', '*** Error in solve_r_eq ***',         &
          ' linear dependencies detected in Cholesky decomposition!'
          print '(a,i0,a,i0)', ' reduced dimensionality: ',ndim_ab,' --> ',effective_dim
          call quit(' *** Error in solve_r_eq: Cholesky decomposition problem!  ***')

       endif

       allocate(temp_mat2(ndim_ab))
       temp_mat2 = val_d0

       call qchols(amat,ndim_ab,effective_dim,ndim_ab,               &
                   nz,ndim_ab,ndim_ab,h1eig,                         &
                   bmat,ndim_ab,ndim_ab,                             &
                   rmat,ndim_ab,ndim_ab,                             &
                   job_task,temp_mat1,temp_mat2)

           if(print_lvl > 2)then
             call print_x2cmtx(rmat,ndim_ab,ndim_ab,nz,ipq_off,'x2c - rmat',6)
           end if

       deallocate(temp_mat2)
       deallocate(temp_mat1)

  end subroutine solve_r_eq

!**********************************************************************
  subroutine construct_1_plus_r_r(xmat,      &
                                  rmat,      &
                                  ndim_r,    & 
                                  nz,        &
                                  ipq_off_in)
!**********************************************************************
!
! 
!    purpose: construct the matrix XMAT containing (1 + R^+.R) 
!
!    equation (10) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
     real(8), intent(inout) :: xmat(*)
     real(8), intent(in)    :: rmat(*)
     integer, intent(in)    :: ndim_r
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(*)
!----------------------------------------------------------------------
     integer                :: j
!**********************************************************************
!

!      step 1: X := 1
       call dunit2(xmat,ndim_r,ndim_r,ndim_r,nz)

!      step 2: X := 1 + R^+.R
       call qgemm(ndim_r,ndim_r,ndim_r,val_d1,               &
                  'H','N',ipq_off_in,                        &
                  rmat,ndim_r,ndim_r,nz,                     &
                  'N','N',ipq_off_in,                        &
                  rmat,ndim_r,ndim_r,nz,                     &
                  val_d1,ipq_off_in,xmat,ndim_r,ndim_r,nz)

  end subroutine construct_1_plus_r_r

!**********************************************************************
  subroutine diag_sqrt_1_plus_r_r(xmat,       &
                                  eigvc,      &
                                  eigvl,      &
                                  ndim_x,     &
                                  nz)
!**********************************************************************
!
! 
!    purpose: 1. diagonalize the matrix X = (1 + R^+.R) 
!             2. take inverse sqrt of eigenvectors of diagonal matrix ((1 + R^+.R))
!
!       note: inverse is taken as we want to derive N_+^{-1} = 1/(1 + R^+.R)^{1/2}
!
!    see also equation (10) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!  
!    for further information on taking the square root of a matrix see, e.g.,
!    http://en.wikipedia.org/wiki/Square_root_of_a_matrix
!
!    X^{1/2} = eigvc . d^{1/2} . eigvc^+
!
     real(8), intent(inout) :: xmat(*)
     real(8), intent(inout) :: eigvc(*)
     real(8), intent(inout) :: eigvl(*)
     integer, intent(in)    :: ndim_x
     integer, intent(in)    :: nz
!----------------------------------------------------------------------
     real(8)                :: scaling_factor
     integer                :: ierr_qdiag
     integer                :: iz
     integer                :: ioff
     integer                :: j
!**********************************************************************
!
       ierr_qdiag =  0

!      step 1: diagonalize matrix X = (1 + R^+.R)
       call qdiag90(nz,ndim_x,xmat,ndim_x,ndim_x,           &
                    eigvl,1,eigvc,ndim_x,ndim_x,ierr_qdiag)
       if( ierr_qdiag /= 0 )                                & 
         call quit(' *** error in diag_sqrt_1_plus_r_r: qdiag returned with error code. ***')

!      copy eigenvectors into xmat
       call dcopy(ndim_x*ndim_x*nz,eigvc,1,xmat,1)

!      step 2: rescale eigenvectors: first step for taking the inverse sqrt
!              of (1 + R^+.R)
       ioff =  1
       do iz = 1, nz
         do j = 1, ndim_x
           scaling_factor = val_d1 / dsqrt(eigvl(j))
!              ....  d^-1/2 . eigvc^+ 
           call dscal(ndim_x,scaling_factor,xmat(ioff),1)
           ioff = ioff + ndim_x
         end do
       end do

  end subroutine diag_sqrt_1_plus_r_r

!**********************************************************************
  subroutine generate_w2(w2mat,xmatd12,xmat,ndim_w,nz,ipq_off_in)
!**********************************************************************
!
! 
!    purpose: generate w2 = (X+ * X)^(-1/2) 
!
!    equation (10) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
     real(8), intent(inout) :: w2mat(*)
     real(8), intent(inout) :: xmatd12(*)
     real(8), intent(inout) :: xmat(*)
     integer, intent(in)    :: ndim_w
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(*)
!----------------------------------------------------------------------
!**********************************************************************
!
!      w2 := X * d^-1/2 * X+
       call qgemm(ndim_w,ndim_w,ndim_w,val_d1,                 &
                  'N','N',ipq_off_in,                          &
                  xmatd12,ndim_w,ndim_w,nz,                    &
                  'H','N',ipq_off_in,                          &
                  xmat,ndim_w,ndim_w,nz,                       &
                  val_d0,ipq_off_in,w2mat,ndim_w,ndim_w,nz)

  end subroutine generate_w2

!**********************************************************************
  subroutine generate_w1(w1mat,rmat,ndim_r,ndim_w,nz)
!**********************************************************************
!
! 
!    purpose: generate w1 = (1 R)
!
!    equation (10) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
     real(8), intent(inout) :: w1mat(*)
     real(8), intent(in)    :: rmat(*)
     integer, intent(in)    :: ndim_r
     integer, intent(in)    :: ndim_w
     integer, intent(in)    :: nz
!----------------------------------------------------------------------
     integer                :: iz
     integer                :: j
     integer                :: ioff_r
     integer                :: ioff_w1
     integer                :: ioff_w1_2
!**********************************************************************
!
!      w1 := (1 R) 
       do iz = 1, nz
         do j = 1, ndim_r
!          copy R into w1mat
           ioff_r  = 1 + (j-1)*ndim_r + (iz-1)*ndim_r**2
           ioff_w1 = 1 + (j-1)*ndim_w + (iz-1)*ndim_r**2 * 2
           call dcopy(ndim_r,rmat(ioff_r),1,w1mat(ioff_w1),1)

!          insert 1 and zero to make the final w1 matrix
           ioff_w1_2 = ioff_w1 + ndim_r
           call dzero(w1mat(ioff_w1_2),ndim_r)
           if(iz == 1) w1mat(ioff_w1_2+j-1) = val_d1

         end do
       end do

  end subroutine generate_w1

!**********************************************************************
  subroutine generate_u_mat(umat,w1mat,w2mat,ndim_w1,ndim_w2,nz,ipq_off_in)
!**********************************************************************
!
! 
!    purpose: obtain the unitary transformation matrix U = w1.w2 
!
!    equation (10) in reference: Ilias and Saue, J. Chem. Phys., 126, 064102 (2007)
!
     real(8), intent(inout) :: umat(*)
     real(8), intent(inout) :: w1mat(*)
     real(8), intent(inout) :: w2mat(*)
     integer, intent(in)    :: ndim_w1
     integer, intent(in)    :: ndim_w2
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(*)
!----------------------------------------------------------------------
!**********************************************************************
!
!      U := w1 . w2
       call qgemm(ndim_w1,         &
                  ndim_w2,         &
                  ndim_w2,         &
                  val_d1,          &
                  'N','N',         &
                  ipq_off_in,      &
                  w1mat,           &
                  ndim_w1,         &
                  ndim_w2,         &
                  nz,              &
                  'N','N',         &
                  ipq_off_in,      &
                  w2mat,           &
                  ndim_w2,         &
                  ndim_w2,         &
                  nz,              &
                  val_d0,          &
                  ipq_off_in,      &
                  umat,            &
                  ndim_w1,         &
                  ndim_w2,         &
                  nz)

  end subroutine generate_u_mat

end module x2c_decoupling_mat
