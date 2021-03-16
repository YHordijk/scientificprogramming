module exacorr_eom

#if defined (VAR_MPI)

        use exatensor
        use exacorr_datatypes
        use intermediates
        use, intrinsic:: ISO_C_BINDING

        implicit none
        complex(8), parameter :: MINUS_ONE=(-1.D0,0.D0),ONE_HALF=(0.5D0,0.D0)

        private

        public exacorr_eom_driver

       contains

        subroutine exacorr_eom_driver (exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)

         type(exacc_input), intent(in) :: exa_input
!        Common tensors
         type(exatns_comm_tens),intent(inout)   :: comm_t
         type(exatns_lambda_tens,intent(inout)  :: lambda_t
!        fixed 2-body tensors (two-electron integrals)
         type(exatns_intg_tens),intent(inout) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t),intent(inout) :: l1_tensor, l2_tensor
!        space ids and roots
         type(space_dims),intent(in) :: dims
!        scalar tensor
         type(tens_rcrsv_t), intent(inout) :: one_tensor
!        Type of EOM
         character(len=2) :: flavour

!        tensor ids + roots
         integer                    :: occ_space_id, vir_space_id
         integer(INTL)              :: occ_space_root, vir_space_root
         integer,dimension(2)       :: oo_id, ov_id, vv_id
         integer(INTL),dimension(2) :: oo_root, ov_root, vv_root
         integer,dimension(4)       :: oooo_id, ovoo_id, ooov_id, vovo_id, vvvo_id, vovv_id, vvvv_id
         integer(INTL),dimension(4) :: oooo_root, ovoo_root, ooov_root, vovo_root, vvvo_root, vovv_root, vvvv_root

         logical :: CCD, left=.true.
         integer :: ierr

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         print*, " Starting EOM-CC calculation with exatensor"

         CCD = exa_input%ccd

!        Create all ids and roots
         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root
         oo_id = (/occ_space_id, occ_space_id/)
         ov_id = (/occ_space_id, vir_space_id/)
         vv_id = (/vir_space_id, vir_space_id/)
         oo_root = (/occ_space_root, occ_space_root/)
         ov_root = (/occ_space_root, vir_space_root/)
         vv_root = (/vir_space_root, vir_space_root/)
         oooo_id = (/occ_space_id, occ_space_id, occ_space_id, occ_space_id/)
         ovoo_id = (/occ_space_id, vir_space_id, occ_space_id, occ_space_id/)
         ooov_id = (/occ_space_id, occ_space_id, occ_space_id, vir_space_id/)
         vovo_id = (/vir_space_id, occ_space_id, vir_space_id, occ_space_id/)
         vvvo_id = (/vir_space_id, vir_space_id, vir_space_id, occ_space_id/)
         vovv_id = (/vir_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vvvv_id = (/vir_space_id, vir_space_id, vir_space_id, vir_space_id/)
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
         ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
         ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
         vovo_root = (/vir_space_root, occ_space_root, vir_space_root, occ_space_root/)
         vvvo_root = (/vir_space_root, vir_space_root, vir_space_root, occ_space_root/)
         vovv_root = (/vir_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vvvv_root = (/vir_space_root, vir_space_root, vir_space_root, vir_space_root/)

!        Construct fixed intermediates
         ierr=exatns_tensor_create(lambda_t%fbar_oo,"fbar_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%fbar_ov,"fbar_ov",ov_id,ov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%fbar_vv,"fbar_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_oooo,"w_oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_ovoo,"w_ovoo",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_ooov,"w_ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vovo,"w_vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%wbar_vovo,"wbar_vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vvvo,"w_vvvo",vvvo_id,vvvo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vovv,"w_vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vvvv,"w_vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%goo,"goo_lambda",oo_id,oo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%gvv,"gvv_lambda",vv_id,vv_root,EXA_DATA_KIND_C8)

!        Get fixed intermediates
         call print_date('Calculating fixed intermediates')
         if (flavour.eq."EA".or.flavour.eq."EE") then
            call get_w_vvvv(lambda_t%w_vvvv,comm_t,int_t,CCD,one_tensor)
            call get_wbar_vovo(lambda_t%wbar_vovo,comm_t%t2,int_t,one_tensor)
         end if
         call get_fbar_ov(lambda_t%fbar_ov,comm_t%hov,one_tensor)
         if (flavour.eq."EA".or.flavour.eq."EE") then
            call get_w_vvvo(lambda_t,comm_t,int_t,CCD,one_tensor)
            call get_w_vovv(lambda_t%w_vovv,comm_t%t1,int_t,CCD,one_tensor)
         end if
         call get_w_vovo(lambda_t%w_vovo,comm_t,int_t%oovv,occ_space_id,occ_space_root, &
                         vir_space_id,vir_space_root,CCD,one_tensor)
         if (flavour.eq."IP".or.flavour.eq."EE") then
            call get_w_oooo(lambda_t%w_oooo,comm_t%a_int,one_tensor)
            call get_w_ovoo(lambda_t,comm_t,int_t,CCD,one_tensor)
            call get_w_ooov(lambda_t%w_ooov,comm_t%t1,int_t,CCD,one_tensor)
         end if
         call get_fbar_vv(lambda_t%fbar_vv,comm_t%gvv,one_tensor)
         call get_fbar_oo(lambda_t%fbar_oo,comm_t,CCD,one_tensor)
         if (flavour.eq."EE") then
            call get_gvv_lambda(lambda_t%gvv,l2_tensor,comm_t%t2)
            call get_goo_lambda(lambda_t%goo,l2_tensor,comm_t%t2)
         end if
         call print_date('Finished calculating fixed intermediates')

         if (flavour.eq."IP") then
            call eom_ip_right(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)
            if (left) call eom_ip_left(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)
         else if (flavour.eq."EA") then
            call eom_ea_right(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)
            if (left) call eom_ea_left(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)
         else if (flavour.eq."EE") then
            call eom_ee_right(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)
            if (left) call eom_ee_left(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)
         else
            print*, "You should pick a valid EOM type."
         end if

         ierr=exatns_tensor_destroy(lambda_t%fbar_oo)
         ierr=exatns_tensor_destroy(lambda_t%fbar_ov)
         ierr=exatns_tensor_destroy(lambda_t%fbar_vv)
         ierr=exatns_tensor_destroy(lambda_t%w_oooo)
         ierr=exatns_tensor_destroy(lambda_t%w_ooov)
         ierr=exatns_tensor_destroy(lambda_t%w_ovoo)
         ierr=exatns_tensor_destroy(lambda_t%w_vovo)
         ierr=exatns_tensor_destroy(lambda_t%w_vvvo)
         ierr=exatns_tensor_destroy(lambda_t%w_vovv)
         ierr=exatns_tensor_destroy(lambda_t%w_vvvv)
         ierr=exatns_tensor_destroy(lambda_t%wbar_vovo)
         ierr=exatns_tensor_destroy(lambda_t%goo)
         ierr=exatns_tensor_destroy(lambda_t%gvv)

!        Make some noise so that we know we are leaving
         print*, " Leaving eom_driver routine"
         return

        end subroutine exacorr_eom_driver

        subroutine eom_ip_right(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)

         type(exacc_input), intent(in) :: exa_input
!        Common tensors
         type(exatns_comm_tens),intent(inout)   :: comm_t
         type(exatns_lambda_tens,intent(inout)  :: lambda_t
!        fixed 2-body tensors (two-electron integrals)
         type(exatns_intg_tens),intent(inout) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t),intent(inout) :: l1_tensor, l2_tensor
!        space ids and roots
         type(space_dims),intent(in) :: dims
!        scalar tensor
         type(tens_rcrsv_t), intent(inout) :: one_tensor

         integer                    :: nocc, nvir
         integer, allocatable       :: mo_occ, mo_vir
         !logical                   :: CCD
         real(8)                    :: t_target_precision
         integer                    :: occ_space_id, vir_space_id
         integer(INTL)              :: occ_space_root, vir_space_root
         integer,dimension(3)       :: voo_id, oov_id
         integer(INTL),dimension(3) :: voo_root, oov_root

!        Solution tensors
         type(tens_rcrsv_t) :: r1_coef, r2_coef,l1_coef,l2_coef
         type(tens_rvrsv_t) :: r1_sigma, r2_sigma, l1_sigma, l2_sigma
!        Auxiliary tensor
         type(tens_rcrsv_t) :: v_aux

         integer :: neom=30
         integer :: ierr, i

         !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nkr_occ
         nvir = exa_input%nkr_vir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ(1:nocc) = exa_input%mokr_occ(1:exa_input%nkr_occ)
         mo_vir(1:nvir) = exa_input%mokr_vir(1:exa_input%nkr_vir)
         t_target_precision = exa_input%t_econv
         !CCD=exa_input%ccd

         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root

!        Create all id and root arrays
         voo_id = (/vir_space_id, occ_space_id, occ_space_id/)
         oov_id = (/occ_space_id, occ_space_id, vir_space_id/)
         voo_root = (/vir_space_root, occ_space_root, occ_space_root/)
         oov_root = (/occ_space_root, occ_space_root, vir_space_root/)

!        Create all coefficient and sigma tensors
         ierr=exatns_tensor_create(r1_coef,"r1_coef",(/occ_space_id/),(/occ_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l1_coef,"l1_coef",(/occ_space_id/),(/occ_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r2_coef,"r2_coef",voo_id,voo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l2_coef,"l2_coef",oov_id,oov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r1_sigma,"r1_sigma",(/occ_space_id/),(/occ_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l1_sigma,"l1_sigma",(/occ_space_id/),(/occ_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r2_sigma,"r2_sigma",voo_id,voo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l2_sigma,"l2_sigma",oov_id,oov_root,EXA_DATA_KIND_C8)

!        Create aux tensor
         ierr=exatns_create(v_aux,"v_aux",(/vir_space_id/),(/vir_space_root/),EXA_DATA_KIND_C8)

         do i = 1, neom

!         r1_sigma
          ierr=exatns_tensor_init(r1_sigma)
          ierr=exatns_tensor_contract("S(i)+=F(m,i)*R(m)",r1_sigma,lambda_t%fbar_oo,r1_coef,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i)+=F(m,e)*R(e,i,m)",r1_sigma,lambda_t%fbar_ov,r2_coef)
          ierr=exatns_tensor_contract("S(i)+=W(m,n,i,e)*R(e,m,n)",r1_sigma,lambda_t%w_ooov,r2_coef,ONE_HALF)

!         r2_sigma
          ierr=exatns_tensor_init(r2_sigma)
          ierr=exatns_tensor_contract("S(a,i,j)+=W(m,a,j,i)*R(m)",r2_sigma,lambda_t%w_ovoo,r1_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_create("S(a,i,j)+=F(m,i)*r(a,j,m)",r2_sigma,lambda_t%fbar_oo,r2_coef)
          ierr=exatns_tensor_create("S(a,i,j)+=F(m,j)*r(a,m,i)",r2_sigma,lambda_t%fbar_oo,r2_coef)
          ierr=exatns_tensor_create("S(a,i,j)+=W(m,n,i,j)*r(a,m,n)",r2_sigma,lambda_t%w_oooo,r2_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_create("S(a,i,j)+=W(a,m,e,i)*R(e,j,m)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_create("S(a,i,j)+=W(a,m,e,j)*R(e,m,i)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_create("S(a,i,j)+=F(a,e)*R(e,i,j)",r2_sigma,lambda_t%fbar_vv,r2_coef)
          ierr=exatns_tensor_init(v_aux)
          !ierr=exatns_tensor_create("U(e)+=W(m,n,e,f)*R(f,m,n)",v_aux,w_oovv,r2_coef,ONE_HALF) ==> NO W_OOVV DEFINED ==> V_OOVV??
          ierr=exatns_tensor_create("S(a,i,j)+=U(e)*T(e,a,i,j)",r2_sigma,v_aux,comm_t%t2)

!         l1_sigma
          ierr=exatns_tensor_init(l1_sigma)
          ierr=exatns_tensor_contract("S(i)+=F(i,m)*L(m)",l1_sigma,lambda_t%fbar_oo,l1_coef,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i)+=W(i,e,m,n)*L(m,n,e)",l1_sigma,lambda_t%w_ovoo,l2_coef,ONE_HALF)

!         l2_sigma
          ierr=exatns_tensor_init(l2_sigma)
          ierr=exatns_tensor_contract("S(i,j,a)+=W(j,i,m,a)*L(m)",l2_sigma,lambda_t%w_ooov,l1_coef)
          ierr=exatns_tensor_contract("S(i,j,a)+=W(i,e,m,n)*L(m,n,e)",l2_sigma,lambda_t%w_ovoo,l2_coef,ONE_HALF)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a)+=F(j,m)*L(m,i,a)",l2_sigma,lambda_t%fbar_oo,l2_coef)
          ierr=exatns_tensor_contract("S(i,j,a)+=F(i,m)*L(j,m,a)",l2_sigma,lambda_t%fbar_oo,l2_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a)+=W(e,j,a,m)*L(m,i,e)",l2_sigma,lambda_t%w_vovo,l2_coef)
          ierr=exatns_tensor_contract("S(i,j,a)+=W(e,i,a,m)*L(j,m,e)",l2_sigma,lambda_t%w_vovo,l2_coef)
          ierr=exatns_tensor_contract("S(i,j,a)+=W(i,j,m,n)*L(m,n,a)",l2_sigma,lambda_t%w_oooo,l2_coef,ONE_HALF)
          ierr=exatns_tensor_init(v_aux)
          ierr=exatns_tensor_contract("U(e)+=L(m,n,f)*T(e,f,m,n)",v_aux,l2_coef,comm_t%t2,ONE_HALF)
          ierr=exatns_tensor_contract("S(i,j,a)+=V(i,j,e,a)*U(e)",l2_sigma,int_t%oovv,v_aux)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a)+=L(i)*F(j,a)",l2_sigma,l1_coef,lambda_t%fbar_ov)
          ierr=exatns_tensor_contract("S(i,j,a)+=L(j)*F(i,a)",l2_sigma,l1_coef,lambda_t%fbar_ov,MINUS_ONE)

         end do

         deallocate(mo_occ)
         deallocate(mo_vir)

         ierr=exatns_tensor_destroy(r1_coef)
         ierr=exatns_tensor_destroy(l1_coef)
         ierr=exatns_tensor_destroy(r2_coef)
         ierr=exatns_tensor_destroy(l2_coef)
         ierr=exatns_tensor_destroy(r1_sigma)
         ierr=exatns_tensor_destroy(l1_sigma)
         ierr=exatns_tensor_destroy(r2_sigma)
         ierr=exatns_tensor_destroy(l2_sigma)
         ierr=exatns_tensor_destroy(v_aux)

        end subroutine eom_ip_right

        subroutine eom_ea_right(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)

         type(exacc_input), intent(in) :: exa_input
!        Common tensors
         type(exatns_comm_tens),intent(inout)   :: comm_t
         type(exatns_lambda_tens,intent(inout)  :: lambda_t
!        fixed 2-body tensors (two-electron integrals)
         type(exatns_intg_tens),intent(inout) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t),intent(inout) :: l1_tensor, l2_tensor
!        space ids and roots
         type(space_dims),intent(in) :: dims
!        scalar tensor
         type(tens_rcrsv_t), intent(inout) :: one_tensor

         integer                    :: nocc, nvir
         integer, allocatable       :: mo_occ, mo_vir
         !logical                   :: CCD
         real(8)                    :: t_target_precision
         integer                    :: occ_space_id, vir_space_id
         integer(INTL)              :: occ_space_root, vir_space_root
         integer,dimension(3)       :: vvo_id, ovv_id
         integer(INTL),dimension(3) :: vvo_root, ovv_root

!        Solution tensors
         type(tens_rcrsv_t) :: r1_coef, r2_coef,l1_coef,l2_coef
         type(tens_rvrsv_t) :: r1_sigma, r2_sigma, l1_sigma, l2_sigma
!        Auxiliary tensor
         type(tens_rcrsv_t) :: o_aux

         integer :: neom=30
         integer :: ierr, i

         !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nkr_occ
         nvir = exa_input%nkr_vir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ(1:nocc) = exa_input%mokr_occ(1:exa_input%nkr_occ)
         mo_vir(1:nvir) = exa_input%mokr_vir(1:exa_input%nkr_vir)
         t_target_precision = exa_input%t_econv
         !CCD=exa_input%ccd

         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root

         !Create all id and root arrays, neccessary for exatns_tensor_create
         vvo_id = (/vir_space_id, vir_space_id, occ_space_id/)
         ovv_id = (/occ_space_id, vir_space_id, vir_space_id/)
         vvo_root = (/vir_space_root, vir_space_root, occ_space_root/)
         ovv_root = (/occ_space_root, vir_space_root, vir_space_root/)

!        Create all coefficient and sigma tensors
         ierr=exatns_tensor_create(r1_coef,"r1_coef",(/vir_space_id/),(/vir_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l1_coef,"l1_coef",(/vir_space_id/),(/vir_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r2_coef,"r2_coef",vvo_id,vvo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l2_coef,"l2_coef",ovv_id,ovv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r1_sigma,"r1_sigma",(/vir_space_id/),(/vir_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l1_sigma,"l1_sigma",(/vir_space_id/),(/vir_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r2_sigma,"r2_sigma",vvo_id,vvo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l2_sigma,"l2_sigma",ovv_id,ovv_root,EXA_DATA_KIND_C8)

!        Create aux tensor
         ierr=exatns_create(o_aux,"o_aux",(/occ_space_id/),(/occ_space_root/),EXA_DATA_KIND_C8)

         do i = 1, neom

!         r1_sigma
          ierr=exatns_tensor_init(r1_sigma)
          ierr=exatns_tensor_contract("S(a)+=F(a,e)*R(e)",r1_sigma,lambda_t%fbar_vv,r1_coef)
          ierr=exatns_tensor_contract("S(a)+=F(m,e)*R(e,a,m)",r1_sigma,lambda_t%fbar_ov,r2_coef)
          ierr=exatns_tensor_contract("S(a)+=W(a,m,f,e)*R(e,f,n)",r1_sigma,lambda_t%w_vovv,r2_coef,ONE_HALF)

!         r2_sigma
          ierr=exatns_tensor_init(r2_sigma)
          ierr=exatns_tensor_contract("S(a,b,i)+=W(a,b,e,i)*R(e)",r2_sigma,lambda_t%w_vvvo,r1_coef)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_create("S(a,b,i)+=F(a,e)*r(e,b,i)",r2_sigma,lambda_t%fbar_vv,r2_coef)
          ierr=exatns_tensor_create("S(a,b,i)+=F(b,e)*r(a,e,i)",r2_sigma,lambda_t%fbar_vv,r2_coef)
          ierr=exatns_tensor_create("S(a,b,i)+=W(a,b,e,f)*r(e,f,i)",r2_sigma,lambda_t%w_vvvv,r2_coef,ONE_HALF)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_create("S(a,b,i)+=W(a,m,e,i)*R(b,e,m)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_create("S(a,b,i)+=W(b,m,e,i)*R(e,a,m)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_create("S(a,b,i)+=F(m,i)*R(b,a,m)",r2_sigma,lambda_t%fbar_oo,r2_coef)
          ierr=exatns_tensor_init(o_aux)
          ierr=exatns_tensor_create("U(m)+=V(m,n,e,f)*R(e,f,n)",o_aux,int_t%oovv,r2_coef,ONE_HALF)
          ierr=exatns_tensor_create("S(a,b,i)+=U(m)*T(a,b,m,i)",r2_sigma,o_aux,comm_t%t2)

!         l1_sigma
          ierr=exatns_tensor_init(l1_sigma)
          ierr=exatns_tensor_contract("S(a)+=F(e,a)*L(e)",l1_sigma,lambda_t%fbar_vv,l1_coef)
          ierr=exatns_tensor_contract("S(a)+=W(e,f,a,m)*L(m,f,e)",l1_sigma,lambda_t%w_vvvo,l2_coef,ONE_HALF)

!         l2_sigma
          ierr=exatns_tensor_init(l2_sigma)
          ierr=exatns_tensor_contract("S(i,a,b)+=W(e,i,b,a)*L(e)",l2_sigma,lambda_t%w_vovv,l1_coef)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(i,a,b)+=F(e,a)*L(i,e,b)",l2_sigma,lambda_t%fbar_vv,l2_coef)
          ierr=exatns_tensor_contract("S(i,a,b)+=F(e,b)*L(i,a,e)",l2_sigma,lambda_t%fbar_vv,l2_coef)
          ierr=exatns_tensor_contract("S(i,a,b)+=F(i,m)*L(m,b,a)",l2_sigma,lambda_t$fbar_oo,l2_coef)
          ierr=exatns_tensor_contract("S(i,a,b)+=W(e,f,a,b)*L(i,e,f)",l2_sigma,lambda_t%w_vvvv,l2_coef,ONE_HALF)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(i,a,b)+=W(e,i,b,m)*L(m,e,a)",l2_sigma,lambda_t%w_vovo,l2_coef)
          ierr=exatns_tensor_contract("S(i,a,b)+=W(e,i,a,m)*L(m,b,e)",l2_sigma,lambda_t%w_vovo,l2_coef)
          ierr=exatns_tensor_init(o_aux)
          ierr=exatns_tensor_contract("U(m)+=L(o,e,f)*T(e,f,o,m)",o_aux,l2_coef,comm_t%t2,ONE_HALF)
          ierr=exatns_tensor_contract("S(i,a,b)+=U(e)*V(m,i,b,a)",l2_sigma,o_aux,int_t%oovv)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(i,a,b)+=F(i,a)*L(b)",l2_sigma,lambda_t%fbar_ov,l1_coef)
          ierr=exatns_tensor_contract("S(i,a,b)+=F(i,b)*L(a)",l2_sigma,lambda_t%fbar_ov,l1_coef,MINUS_ONE)

         end do

         deallocate(mo_occ)
         deallocate(mo_vir)

         ierr=exatns_tensor_destroy(r1_coef)
         ierr=exatns_tensor_destroy(l1_coef)
         ierr=exatns_tensor_destroy(r2_coef)
         ierr=exatns_tensor_destroy(l2_coef)
         ierr=exatns_tensor_destroy(r1_sigma)
         ierr=exatns_tensor_destroy(l1_sigma)
         ierr=exatns_tensor_destroy(r2_sigma)
         ierr=exatns_tensor_destroy(l2_sigma)
         ierr=exatns_tensor_destroy(o_aux)

        end subroutine eom_ea_right

        subroutine eom_ee_right(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)

         type(exacc_input), intent(in) :: exa_input
!        Common tensors
         type(exatns_comm_tens),intent(inout)   :: comm_t
         type(exatns_lambda_tens,intent(inout)  :: lambda_t
!        fixed 2-body tensors (two-electron integrals)
         type(exatns_intg_tens),intent(inout) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t),intent(inout) :: l1_tensor, l2_tensor
!        space ids and roots
         type(space_dims),intent(in) :: dims
!        scalar tensor
         type(tens_rcrsv_t), intent(inout) :: one_tensor

         integer                    :: nocc, nvir
         integer, allocatable       :: mo_occ, mo_vir
         !logical                   :: CCD
         real(8)                    :: t_target_precision
         integer                    :: occ_space_id, vir_space_id
         integer(INTL)              :: occ_space_root, vir_space_root
         integer,dimension(2)       :: oo_id, ov_id, vo_id, vv_id
         integer,dimension(4)       :: oovv_id, vvoo_id
         integer(INTL),dimension(2) :: oo_root, ov_root, vo_root, vv_root
         integer(INTL),dimension(4) :: oovv_root, vvoo_root

!        Solution tensors
         type(tens_rcrsv_t) :: r1_coef, r2_coef,l1_coef,l2_coef
         type(tens_rvrsv_t) :: r1_sigma, r2_sigma, l1_sigma, l2_sigma
!        Auxiliary tensor
         type(tens_rcrsv_t) :: oo_aux, vv_aux

         integer :: neom=30
         integer :: ierr, i

         !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nkr_occ
         nvir = exa_input%nkr_vir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ(1:nocc) = exa_input%mokr_occ(1:exa_input%nkr_occ)
         mo_vir(1:nvir) = exa_input%mokr_vir(1:exa_input%nkr_vir)
         t_target_precision = exa_input%t_econv
         !CCD=exa_input%ccd

         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root

         !Create all id and root arrays, neccessary for exatns_tensor_create
         oo_id = (/occ_space_id, occ_space_id/)
         ov_id = (/occ_space_id, vir_space_id/)
         vo_id = (/vir_space_id, occ_space_id/)
         vv_id = (/vir_space_id, vir_space_id/)
         oo_root = (/occ_space_root, occ_space_root/)
         ov_root = (/occ_space_root, vir_space_root/)
         vo_root = (/vir_space_root, occ_space_root/)
         vv_root = (/vir_space_root, vir_space_root/)
         oovv_id = (/occ_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vvoo_id = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
         oovv_root = (/occ_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)

!        Create all coefficient and sigma tensors
         ierr=exatns_tensor_create(r1_coef,"r1_coef",vo_id,vo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l1_coef,"l1_coef",ov_id,ov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r2_coef,"r2_coef",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l2_coef,"l2_coef",oovv_id,oovv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r1_sigma,"r1_sigma",vo_id,vo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l1_sigma,"l1_sigma",ov_id,ov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r2_sigma,"r2_sigma",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l2_sigma,"l2_sigma",oovv_id,oovv_root,EXA_DATA_KIND_C8)

!        Create aux tensors
         ierr=exatns_create(oo_aux,"oo_aux",oo_id,oo_root,EXA_DATA_KIND_C8)
         ierr=exatns_create(vv_aux,"vv_aux",vv_id,vv_root,EXA_DATA_KIND_C8)

         do i = 1, neom

!         r1_sigma
          ierr=exatns_tensor_init(r1_sigma)
          ierr=exatns_tensor_contract("S(a,i)+=F(a,e)*R(e,i)",r1_sigma,lambda_t%fbar_vv,r1_coef)
          ierr=exatns_tensor_contract("S(a,i)+=F(m,i)*R(a,m)",r1_sigma,lambda_t%fbar_oo,r1_coef,MINUS_ONE)
          ierr=exatns_tensor_contract("S(a,i)+=F(m,e)*R(e,a,m,i)",r1_sigma,lambda_t%fbar_ov,r2_coef)
          ierr=exatns_tensor_contract("S(a,i)+=W(a,m,e,i)*R(e,m)",r1_sigma,lambda_t%w_vovo,r1_coef,MINUS_ONE)
          ierr=exatns_tensor_contract("S(a,i)+=W(a,m,e,f)*R(e,f,i,m)",r1_sigma,lambda_t%w_vovv,r2_coef,ONE_HALF)
          ierr=exatns_tensor_contract("S(a,i)+=W(m,n,i,e)*R(e,a,n,m)",r1_sigma,lambda_t%w_ooov,r2_coef,ONE_HALF)

!         r2_sigma
          ierr=exatns_tensor_init(r2_sigma)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(a,b,i,j)+=W(m,b,j,i)*R(a,m)",r2_sigma,lambda_t%w_ovoo,r1_coef)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=W(m,a,i,j)*R(b,m)",r2_sigma,lambda_t%w_ovoo,r1_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(a,b,i,j)+=W(a,b,e,j)*R(e,i)",r2_sigma,lambda_t%w_vvvo,r1_coef)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=W(b,a,e,i)*R(e,j)",r2_sigma,lambda_t%w_vvvo,r1_coef)
          ierr=exatns_tensor_init(vv_aux)
          ierr=exatns_tensor_contract("U(b,f)+=W(b,m,f,e)*R(e,m)",vv_aux,lambda_t%w_vovv,r1_coef)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(a,b,i,j)+=U(b,f)*T(a,f,i,j)",r2_sigma,vv_aux,comm_t%t2)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=U(a,f)*T(f,b,i,j)",r2_sigma,vv_aux,comm_t%t2)
          ierr=exatns_tensor_init(oo_aux)
          ierr=exatns_tensor_contract("U(n,j)+=W(n,m,j,e)*R(e,m)",oo_aux,lambda_t%w_ooov,r1_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(a,b,i,j)+=U(n,j)*T(a,b,i,n)",r2_sigma,oo_aux,comm_t%t2)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=U(n,i)*T(a,b,n,j)",r2_sigma,oo_aux,comm_t%t2)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(a,b,i,j)+=F(b,e)*R(a,e,i,j)",r2_sigma,lambda_t%fbar_vv,r2_coef)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=F(a,e)*R(e,b,i,j)",r2_sigma,lambda_t%fbar_vv,r2_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(a,b,i,j)+=F(m,j)*R(a,b,m,i)",r2_sigma,lambda_t%fbar_vv,r2_coef)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=F(m,i)*R(a,b,j,m)",r2_sigma,lambda_t%fbar_vv,r2_coef)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=W(m,n,i,j)*R(a,b,m,n)",r2_sigma,lambda_t%w_oooo,r2_coef,ONE_HALF)
!         Permute i and j, as well as  a and b, so term appears four times
          ierr=exatns_tensor_contract("S(a,b,i,j)+=W(b,m,e,j)*R(a,e,m,i)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=W(a,m,e,j)*R(b,e,i,m)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=W(b,m,e,i)*R(a,e,j,m)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=W(a,m,e,i)*R(b,e,m,j)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_init(vv_aux)
          ierr=exatns_tensor_contract("U(a,f)+=V(m,n,f,e)*R(e,a,m,n)",vv_aux,int_t%oovv,r2_coef,ONE_HALF)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(a,b,i,j)+=U(a,f)*T(f,b,i,j)",r2_sigma,vv_aux,comm_t%t2)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=U(b,f)*T(a,f,i,j)",r2_sigma,vv_aux,comm_t%t2)
          ierr=exatns_tensor_init(oo_aux)
          ierr=exatns_tensor_contract("U(n,i)+=V(n,m,f,e)*R(f,e,i,m)",oo_aux,int_t%oovv,r2_coef,ONE_HALF)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(a,b,i,j)+=U(n,i)*T(a,b,j,n)",r2_sigma,oo_aux,comm_t%t2)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=U(n,j)*T(b,a,i,n)",r2_sigma,oo_aux,comm_t%t2)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=W(a,b,e,f)*R(e,f,i,j)",r2_sigma,lambda_t%w_vvvv,r2_coef,ONE_HALF)

!         l1_sigma
          ierr=exatns_tensor_init(l1_sigma)
          ierr=exatns_tensor_contract("S(i,a)+=F(i,a)*K()",l1_sigma,lambda_t%fbar_ov,one_tensor)
          ierr=exatns_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",l1_sigma,l1_coef,lambda_t%fbar_vv)
          ierr=exatns_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",l1_sigma,l1_coef,lambda_t%fbar_oo,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",l1_sigma,l2_coef,lambda_t%w_vvvo,ONE_HALF)
          ierr=exatns_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",l1_sigma,lambda_t%gvv,lambda_t%w_vovv)
          ierr=exatns_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",l1_sigma,lambda_t%goo,lambda_t%w_ooov)
          ierr=exatns_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",l1_sigma,l1_coef,lambda_t%w_vovo,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",l1_sigma,l2_coef,lambda_t%w_ovoo,ONE_HALF)

!         l2_sigma
          ierr=exatns_tensor_init(l2_sigma)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,b)*K()",l2_sigma,int_t%oovv,one_tensor)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",l2_sigma,l2_coef,lambda_t%fbar_vv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",l2_sigma,l2_coef,lambda_t%fbar_vv)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,i,a,b)*F(j,m)",l2_sigma,l2_coef,lambda_t%fbar_oo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",l2_sigma,l2_coef,lambda_t%fbar_oo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",l2_sigma,l2_coef,lambda_t%w_oooo,ONE_HALF)
!         Permute i and j, as well as  a and b, so term appears four times
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",l2_sigma,l2_coef,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",l2_sigma,l2_coef,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",l2_sigma,l2_coef,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",l2_sigma,l2_coef,lambda_t%w_vovo)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",l2_sigma,int_t%oovv,lambda_t%gvv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",l2_sigma,int_t%oovv,lambda_t%gvv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,a)*W(i,j,m,b)",l2_sigma,l1_coef,lambda_t%w_ooov)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,m,b,a)*G(j,m)",l2_sigma,int_t%oovv,lambda_t%goo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",l2_sigma,int_t%oovv,lambda_t%goo)
          ierr=exatns_tensor_init(oo_aux)
          ierr=exatns_tensor_contract("U(i,m)+=L(i,e)*T(e,m)",oo_aux,l1_coef,comm_t%t1)
!         Permute i and j, so term appears twice
          ierr=exatans_tensor_contract("S(i,j,a,b)+=V(m,j,a,b)*U(i,m)",l2_sigma,int_t%oovv,oo_aux)
          ierr=exatans_tensor_contract("S(i,j,a,b)+=V(i,m,a,b)*U(j,m)",l2_sigma,int_t%oovv,oo_aux)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,e)*V(e,j,a,b)",l2_sigma,l1_coef,int_t%vovv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,e)*V(e,i,b,a)",l2_sigma,l1_coef,int_t%vovv)
!         Permute i and j, as well as  a and b, so term appears four times
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",l2_sigma,l1_coef,lambda_t%fbar_ov)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",l2_sigma,l1_coef,lambda_t%fbar_ov,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",l2_sigma,l1_coef,lambda_t%fbar_ov,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",l2_sigma,l1_coef,lambda_t%fbar_ov)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",l2_sigma,l2_coef,lambda_t%w_vvvv,ONE_HALF)

         end do

         deallocate(mo_occ)
         deallocate(mo_vir)

        end subroutine eom_ee_right

end module exacorr_eom
