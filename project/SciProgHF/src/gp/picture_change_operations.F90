module picture_change_operations

  implicit none

  public perform_pct_saao_bas
  public pick_LL_block_saao_bas

  private

contains

!**********************************************************************
  subroutine perform_pct_saao_bas(op_matrix_2c_ll,     &
                                  pctmat,              &
                                  scr_mat1,            &
                                  naosh_ls,            &
                                  nr_ao_L,             &
                                  naosh_all,           &
                                  naosh_L,             &
                                  nfsym,               &
                                  nz,                  &
                                  ipq_off_in,          &
                                  ioff_aomx,           &
                                  op_fer_rep,          &
                                  op_bs_to_fs,         &
                                  op_2c_ll,            &
                                  print_lvl)   
!**********************************************************************
!
!    purpose: perform the picture-change transformation of a four-component
!             (property) operator X_4c using the total pct-matrix U:
!
!             ===================================
!             equation: X_2c^++ = U^+ * X_4c * U
!             ===================================
!
!             if op_2c_ll == .true. : return only the upper op_2c^++ block (remaining parts zeroed out)
!             if op_2c_ll == .false.: return full four-component matrix block containing the
!                                     pct-transformed op_2c^++ block
!
!             input : 4c-operator matrix    --> op_matrix_2c_ll; pct-matrix --> pctmat 
!             output: 2c^++-operator matrix --> op_matrix_2c_ll
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: op_matrix_2c_ll(naosh_ls**2,nz)
     real(8), intent(inout) :: scr_mat1(naosh_ls**2,nz)
     real(8), intent(in)    :: pctmat(*)
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: nr_ao_L
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: naosh_L(nfsym)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)    :: print_lvl
     integer, intent(in)    :: op_fer_rep
     integer, intent(in)    :: op_bs_to_fs(0:7,1:2)
     logical, intent(in)    :: op_2c_ll
!----------------------------------------------------------------------
     integer                :: i
     integer                :: k
     integer                :: ik
     integer                :: nbas1_f_i
     integer                :: nbasL_f_i
     integer                :: nbas1_f_ik
     integer                :: nbasL_f_ik
     integer                :: op_sym
     integer                :: ioff_pct_mat_dag
     integer                :: ioff_pct_mat
!**********************************************************************

!      initialize offset for pctmat+
       ioff_pct_mat_dag = 1

!      step 1: perform pct in DIRAC-sorted SA-AO basis
!      ------------------------------------------------------------------------------

!      operator: map fermion to boson symmetry
       op_sym = op_bs_to_fs(op_fer_rep,1)

       do i = 1, nfsym

         ik = mod(i+op_sym,2) + 1

!        set dimensions for matrices used in the subroutine
         nbas1_f_i  = naosh_all(i)
         nbasL_f_i  = naosh_L(i)
         nbas1_f_ik = naosh_all(ik)
         nbasL_f_ik = naosh_L(ik)

         if(naosh_L(i) > 0 .and. naosh_L(ik) > 0)then

!          determine offset for pctmat
           ioff_pct_mat = 1
           if(ik == 2) ioff_pct_mat = 1 + naosh_L(1) * naosh_all(1) * nz

!          perform picture-change transformation: 
!          X_2c = pctmat+ * X_4c * pctmat
!          scr_mat1 = pctmat+ * op_matrix_2c_ll * pctmat
           call qtrans90('AOMO','S',0.0d0,                       &
                         nbas1_f_i,nbas1_f_ik,nbasL_f_i,nbasL_f_ik,   &
                         op_matrix_2c_ll(1+ioff_aomx(i,ik),1), naosh_ls , naosh_ls , nz, ipq_off_in(1,op_fer_rep),   & ! parent matrix
                         scr_mat1(1+ioff_aomx(i,ik),1)       , naosh_ls , naosh_ls , nz, ipq_off_in(1,op_fer_rep),   & ! transformed matrix
                         pctmat(ioff_pct_mat_dag)            , nbas1_f_i, nbasL_f_i, nz, ipq_off_in(1,0),            & ! left transformation
                         pctmat(ioff_pct_mat),nbas1_f_ik,nbasL_f_ik,nz,ipq_off_in(1,0),                              & ! right transformation
                         print_lvl)
          
         end if

!        update offset for pctmat+
         ioff_pct_mat_dag  = ioff_pct_mat_dag + naosh_L(i) * naosh_all(i) * nz

       end do

!      step 2: return pct-operator in DIRAC-sorted SA-AO basis on array op_matrix_2c_ll
!      ------------------------------------------------------------------------------
!              case a. op_2c_ll == .true. : return upper op_2c++ block (remaining parts zeroed out)
!              case b. op_2c_ll == .false.: return full 4c block with upper op_2c++ block containing the pct-operator

       if(op_2c_ll)then
         op_matrix_2c_ll = 0.0d0
         call pick_LL_block_saao_bas(scr_mat1,            &
                                     op_matrix_2c_ll,     &
                                     naosh_L,             &
                                     naosh_all,           &
                                     naosh_ls,            &
                                     nr_ao_L,             &
                                     nfsym,               &
                                     nz,                  &
                                     op_fer_rep,          &
                                     op_bs_to_fs)
       else 
         call dcopy(naosh_ls**2 * nz,scr_mat1,1,op_matrix_2c_ll,1)
       end if

  end subroutine perform_pct_saao_bas

!**********************************************************************

  subroutine pick_LL_block_saao_bas(opmat_4c,            &
                                    opmat_2c,            &  
                                    naosh_L,             &
                                    naosh_all,           &
                                    naosh_ls,            &
                                    nr_ao_L,             &
                                    nfsym,               &
                                    nz,                  &
                                    op_fer_rep,          &
                                    op_bs_to_fs)          
!**********************************************************************
!
!    purpose: 
!             pick out LL-block of a given operator in full 4c-SA-AO basis
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: opmat_4c(*)
     real(8), intent(inout) :: opmat_2c(*)
     integer, intent(in)    :: naosh_L(nfsym)
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: nr_ao_L
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: op_fer_rep
     integer, intent(in)    :: op_bs_to_fs(0:7,1:2)
!----------------------------------------------------------------------
     integer                :: i
     integer                :: j
     integer                :: k
     integer                :: jk
     integer                :: op_sym
     integer                :: ioff_off_LS
     integer                :: ioff_off_LS_f
     integer                :: ioff_off_L
     integer                :: ioff_off_L_f
!**********************************************************************

!      operator: map fermion to boson symmetry
       op_sym = op_bs_to_fs(op_fer_rep,1)
      
       do i = 1, nz

!        initialize offset
         ioff_off_LS = (i-1) * naosh_ls * naosh_ls + 1
         ioff_off_L  = (i-1) * nr_ao_L  * nr_ao_L  + 1
 
         do j = 1, nfsym

           jk   = mod(j+op_sym,2) + 1

           if(j == 1 .and. jk == 1)then
             ioff_off_LS_f = ioff_off_LS
             ioff_off_L_f  = ioff_off_L
           else if (j == 1 .and. jk == 2)then 
             ioff_off_LS_f = ioff_off_LS + naosh_ls * naosh_all(1)
             ioff_off_L_f  = ioff_off_L  + nr_ao_L  * naosh_L(1)
           else if (j == 2 .and. jk == 1)then 
             ioff_off_LS_f = ioff_off_LS + naosh_all(1)
             ioff_off_L_f  = ioff_off_L  + naosh_L(1)
           else if (j == 2 .and. jk == 2)then 
             ioff_off_LS_f = ioff_off_LS + naosh_ls *naosh_all(1)  + naosh_all(1)
             ioff_off_L_f  = ioff_off_L  + nr_ao_L  *naosh_L(1)    + naosh_L(1)
           end if

           do k = 1, naosh_L(jk)
             call dcopy(naosh_L(j),                                   &
                        opmat_4c(ioff_off_LS_f + (k-1) * naosh_ls),   &
                        1,                                            &
                        opmat_2c(ioff_off_L_f  + (k-1) * nr_ao_L),    &
                        1)
           end do

         end do
       end do

  end subroutine pick_LL_block_saao_bas
!**********************************************************************

end module
