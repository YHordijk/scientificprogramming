module intermediates

!Module to make intermediates needed in CC and Lambda equations

     implicit none
     private

     complex(8), parameter :: ZERO=(0.D0,0.D0), MINUS_ONE=(-1.D0,0.D0), ONE=(1.D0,0.D0), &
                              ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0), &
                              MINUS_ONE_HALF=(-0.5D0,0.D0), &
                              MINUS_TWO=(-2.D0,0.D0), TWO=(2.D0,0.D0)
    real(8), parameter     :: ONE_QUARTER=0.25D0

     interface get_tau
#if defined (VAR_MPI)
        module procedure get_tau_exatns
#endif
        module procedure get_tau_talsh
     end interface get_tau

     interface get_hoo
#if defined (VAR_MPI)
        module procedure get_hoo_exatns
#endif
        module procedure get_hoo_talsh
     end interface get_hoo

     interface get_hov
#if defined (VAR_MPI)
        module procedure get_hov_exatns
#endif
        module procedure get_hov_talsh
     end interface get_hov

     interface get_hvv
#if defined (VAR_MPI)
        module procedure get_hvv_exatns
#endif
        module procedure get_hvv_talsh
     end interface get_hvv

     interface get_goo
#if defined (VAR_MPI)
        module procedure get_goo_exatns
#endif
        module procedure get_goo_talsh
     end interface get_goo

     interface get_gvv
#if defined (VAR_MPI)
        module procedure get_gvv_exatns
#endif
        module procedure get_gvv_talsh
     end interface get_gvv

     interface get_a
#if defined (VAR_MPI)
        module procedure get_a_exatns
#endif
        module procedure get_a_talsh
     end interface get_a

     interface get_h
#if defined (VAR_MPI)
        module procedure get_h_exatns
#endif
        module procedure get_h_talsh
     end interface get_h

     interface get_fbar_oo
#if defined (VAR_MPI)
        module procedure get_fbar_oo_exatns
#endif
        module procedure get_fbar_oo_talsh
     end interface get_fbar_oo

     interface get_fbar_ov
#if defined (VAR_MPI)
        module procedure get_fbar_ov_exatns
#endif
        module procedure get_fbar_ov_talsh
     end interface get_fbar_ov

     interface get_fbar_vv
#if defined (VAR_MPI)
        module procedure get_fbar_vv_exatns
#endif
        module procedure get_fbar_vv_talsh
     end interface get_fbar_vv

     interface get_w_oooo
#if defined (VAR_MPI)
        module procedure get_w_oooo_exatns
#endif
        module procedure get_w_oooo_talsh
     end interface get_w_oooo

     interface get_w_ooov
#if defined (VAR_MPI)
        module procedure get_w_ooov_exatns
#endif
        module procedure get_w_ooov_talsh
     end interface get_w_ooov

     interface get_w_ovoo
#if defined (VAR_MPI)
        module procedure get_w_ovoo_exatns
#endif
        module procedure get_w_ovoo_talsh
     end interface get_w_ovoo

     interface get_w_vovo
#if defined (VAR_MPI)
        module procedure get_w_vovo_exatns
#endif
        module procedure get_w_vovo_talsh
     end interface get_w_vovo

     interface get_wbar_vovo
#if defined (VAR_MPI)
        module procedure get_wbar_vovo_exatns
#endif
        module procedure get_wbar_vovo_talsh
     end interface get_wbar_vovo

     interface get_w_vvvo
#if defined (VAR_MPI)
        module procedure get_w_vvvo_exatns
#endif
        module procedure get_w_vvvo_talsh
     end interface get_w_vvvo

     interface get_w_vovv
#if defined (VAR_MPI)
        module procedure get_w_vovv_exatns
#endif
        module procedure get_w_vovv_talsh
     end interface get_w_vovv

     interface get_w_vvvv
#if defined (VAR_MPI)
        module procedure get_w_vvvv_exatns
#endif
        module procedure get_w_vvvv_talsh
     end interface get_w_vvvv

     interface get_goo_lambda
#if defined (VAR_MPI)
        module procedure get_goo_lambda_exatns
#endif
        module procedure get_goo_lambda_talsh
     end interface get_goo_lambda

     interface get_gvv_lambda
#if defined (VAR_MPI)
        module procedure get_gvv_lambda_exatns
#endif
        module procedure get_gvv_lambda_talsh
     end interface get_gvv_lambda

!! Start CC2
     interface get_tauCC2 !L
#if defined (VAR_MPI)
        module procedure get_tauCC2_exatns
#endif
        module procedure get_tauCC2_talsh
     end interface get_tauCC2

     interface get_fbar_oo_CC2
#if defined (VAR_MPI)
        module procedure get_fbar_oo_CC2_exatns
#endif
        module procedure get_fbar_oo_CC2_talsh
     end interface get_fbar_oo_CC2

     interface get_fbar_vv_CC2
#if defined (VAR_MPI)
        module procedure get_fbar_vv_CC2_exatns
#endif
        module procedure get_fbar_vv_CC2_talsh
     end interface get_fbar_vv_CC2

!fbar_ov_CC2 Not usefulll same as CCSD but here to check equations
     interface get_fbar_ov_CC2
#if defined (VAR_MPI)
        module procedure get_fbar_ov_CC2_exatns
#endif
        module procedure get_fbar_ov_CC2_talsh
     end interface get_fbar_ov_CC2


     interface get_w_vvvv_CC2
#if defined (VAR_MPI)
        module procedure get_w_vvvv_CC2_exatns
#endif
        module procedure get_w_vvvv_CC2_talsh
     end interface get_w_vvvv_CC2

     interface get_w_oooo_CC2
#if defined (VAR_MPI)
        module procedure get_w_oooo_CC2_exatns
#endif
        module procedure get_w_oooo_CC2_talsh
     end interface get_w_oooo_CC2

     interface get_w_vovo_CC2
#if defined (VAR_MPI)
        module procedure get_w_vovo_CC2_exatns
#endif
        module procedure get_w_vovo_CC2_talsh
     end interface get_w_vovo_CC2

!! End CC2


!    CCSD intermediates
     public get_tau
     public get_hoo
     public get_hov
     public get_hvv
     public get_goo
     public get_gvv
     public get_a
     public get_h
!    Lambda intermediates
     public get_fbar_oo
     public get_fbar_ov
     public get_fbar_vv
     public get_w_oooo
     public get_w_ooov
     public get_w_ovoo
     public get_w_vovo
     public get_wbar_vovo
     public get_w_vvvo
     public get_w_vovv
     public get_w_vvvv
     public get_goo_lambda
     public get_gvv_lambda
!    Lambda CC2 Intemediates
     public get_tauCC2  !L
     public get_fbar_oo_CC2
     public get_fbar_ov_CC2 ! Not usefull here just to check equations
     public get_fbar_vv_CC2
     public get_w_vovo_CC2
     public get_w_oooo_CC2
     public get_w_vvvv_CC2


     contains

#if defined (VAR_MPI)
       subroutine get_tau_exatns(comm_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        logical, intent(in) :: CCD
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       make tau ! validated
        ierr=exatns_tensor_init(comm_t%tau,'ZERO')
        ierr=exatns_tensor_contract("T(a,b,i,j)+=R(a,b,i,j)*K()",comm_t%tau,comm_t%t2,one_tensor)
        if (.not.CCD) then
            ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",comm_t%tau,comm_t%t1,comm_t%t1) ! o(2)v(2)
            ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",comm_t%tau,comm_t%t1,comm_t%t1,MINUS_ONE) ! o(2)v(2)
        end if

       end subroutine get_tau_exatns


       subroutine get_hoo_exatns(comm_t,oovv_tensor,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: oovv_tensor
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       oo ! validated
        ierr=exatns_tensor_init(comm_t%hoo,'ZERO')
        ierr=exatns_tensor_contract("H(k,i)+=F(k,i)*K()",comm_t%hoo,comm_t%foo,one_tensor, MINUS_ONE)
        ierr=exatns_tensor_contract("H(k,i)+=V(k,l,c,d)*T(c,d,l,i)",comm_t%hoo,oovv_tensor,comm_t%tau,ONE_HALF) ! o(3)v(2)

       end subroutine get_hoo_exatns

       subroutine get_hov_exatns(comm_t,oovv_tensor,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: oovv_tensor
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       ov ! validated
        ierr=exatns_tensor_init(comm_t%hov,'ZERO')
        ierr=exatns_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",comm_t%hov,oovv_tensor,comm_t%t1) ! o(2)v(2)
        ierr=exatns_tensor_contract("H(k,c)+=F(k,c)*K()",comm_t%hov,comm_t%fov,one_tensor)           ! add fock matrix

       end subroutine get_hov_exatns

       subroutine get_hvv_exatns(comm_t,oovv_tensor,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: oovv_tensor
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       vv ! validated
        ierr=exatns_tensor_init(comm_t%hvv,'ZERO')
        ierr=exatns_tensor_contract("H(a,c)+=F(a,c)*K()",comm_t%hvv,comm_t%fvv,one_tensor)
        ierr=exatns_tensor_contract("H(a,c)+=V(k,l,c,d)*T(d,a,k,l)",comm_t%hvv,oovv_tensor,comm_t%tau,ONE_HALF) ! o(2)v(3)

       end subroutine get_hvv_exatns

       subroutine get_goo_exatns(comm_t,ooov_tensor,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: ooov_tensor
        logical, intent(in) :: CCD
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(comm_t%goo,'ZERO')
        ierr=exatns_tensor_contract("G(k,i)+=H(k,i)*K()",comm_t%goo,comm_t%hoo,one_tensor)
        if (.not.CCD) then 
            ierr=exatns_tensor_contract("G(k,i)+=V(k,l,i,c)*T(c,l)",comm_t%goo,ooov_tensor,comm_t%t1,MINUS_ONE) ! o(3)v(1)
            ierr=exatns_tensor_contract("G(k,i)+=F(k,c)*T(c,i)",comm_t%goo,comm_t%fov,comm_t%t1,MINUS_ONE) ! o(2)v(1)
        end if

       end subroutine get_goo_exatns

       subroutine get_gvv_exatns(comm_t,vovv_tensor,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: vovv_tensor
        logical, intent(in) :: CCD
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(comm_t%gvv,'ZERO')
        ierr=exatns_tensor_contract("G(a,c)+=H(a,c)*K()",comm_t%gvv,comm_t%hvv,one_tensor)
        if (.not.CCD) then
            ierr=exatns_tensor_contract("G(a,c)+=V(a,k,c,d)*T(d,k)",comm_t%gvv,vovv_tensor,comm_t%t1) ! o(1)v(3)
            ierr=exatns_tensor_contract("G(a,c)+=F(k,c)*T(a,k)",comm_t%gvv,comm_t%fov,comm_t%t1,MINUS_ONE) ! o(1)v(2)
        end if

       end subroutine get_gvv_exatns

       subroutine get_a_exatns(comm_t,int_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens,exatns_intg_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(exatns_intg_tens),intent(inout) :: int_t
        logical, intent(in) :: CCD
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       oooo ! validated
        ierr=exatns_tensor_init(comm_t%a_int,'ZERO')
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,j)*K()",comm_t%a_int,int_t%oooo,one_tensor)
        if (.not.CCD) then 
            ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",comm_t%a_int,int_t%ooov,comm_t%t1) ! o(4)v(1)
            ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",comm_t%a_int,int_t%ooov,comm_t%t1,MINUS_ONE) ! o(4)v(1)
        end if
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",comm_t%a_int,int_t%oovv,comm_t%tau,ONE_HALF) ! o(4)v(2)

       end subroutine get_a_exatns

       subroutine get_h_exatns(comm_t,int_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens,exatns_intg_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(exatns_intg_tens),intent(inout) :: int_t
        logical, intent(in) :: CCD
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       vovo ! more tests
        ierr=exatns_tensor_init(comm_t%h_int,'ZERO')
        if (.not.CCD) then
            ierr=exatns_tensor_contract("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",comm_t%h_int,int_t%vovv,comm_t%t1) ! o(2)v(3)
            ierr=exatns_tensor_contract("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",comm_t%h_int,int_t%ooov,comm_t%t1,MINUS_ONE) ! o(3)v(2)
        end if
        ierr=exatns_tensor_contract("H(a,k,c,i)+=V(a,k,c,i)*K()",comm_t%h_int,int_t%vovo,one_tensor)

       end subroutine get_h_exatns

       subroutine get_fbar_oo_exatns(fbar_oo,comm_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: fbar_oo
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        logical,intent(in) :: CCD
        integer :: ierr

        ierr=exatns_tensor_init(fbar_oo,'ZERO')
        ierr=exatns_tensor_contract("F(i,m)+=G(i,m)*K()",fbar_oo,comm_t%goo,one_tensor,MINUS_ONE)
        if (.not.CCD) ierr=exatns_tensor_contract("F(i,m)+=R(i,e)*T(e,m)",fbar_oo,comm_t%fov,comm_t%t1,TWO)

       end subroutine get_fbar_oo_exatns

       subroutine get_fbar_ov_exatns(fbar_ov,hov_tensor,one_tensor)

        use exatensor

        type(tens_rcrsv_t),intent(inout) :: hov_tensor
        type(tens_rcrsv_t),intent(inout) :: fbar_ov
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(fbar_ov,'ZERO')
        ierr=exatns_tensor_contract("F(m,e)+=H(m,e)*K()",fbar_ov,hov_tensor,one_tensor)

       end subroutine get_fbar_ov_exatns

       subroutine get_fbar_vv_exatns(fbar_vv,gvv_tensor,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens

        type(tens_rcrsv_t),intent(inout) :: gvv_tensor
        type(tens_rcrsv_t),intent(inout) :: fbar_vv
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(fbar_vv,'ZERO')
        ierr=exatns_tensor_contract("F(e,a)+=G(e,a)*K()",fbar_vv,gvv_tensor,one_tensor)

       end subroutine get_fbar_vv_exatns

       subroutine get_w_oooo_exatns(w_oooo,a_int,one_tensor)

        use exatensor

        type(tens_rcrsv_t),intent(inout) :: w_oooo
        type(tens_rcrsv_t),intent(inout) :: a_int
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(w_oooo,'ZERO')
        ierr=exatns_tensor_contract("W(i,j,m,n)+=A(i,j,m,n)*K()",w_oooo,a_int,one_tensor)

       end subroutine get_w_oooo_exatns

       subroutine get_w_ooov_exatns(w_ooov,t1_tensor,int_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens

        type(tens_rcrsv_t),intent(inout) :: t1_tensor
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: w_ooov
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        logical, intent(in) :: CCD
        integer :: ierr

!       ooov
        ierr=exatns_tensor_init(w_ooov,'ZERO')
        if (.not.CCD) ierr=exatns_tensor_contract("W(m,n,i,e)+=T(f,i)*V(m,n,f,e)",w_ooov,t1_tensor,int_t%oovv)
        ierr=exatns_tensor_contract("W(m,n,i,e)+=V(m,n,i,e)*K()",w_ooov,int_t%ooov,one_tensor)

       end subroutine get_w_ooov_exatns

       subroutine get_w_ovoo_exatns(lambda_t,comm_t,int_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens,exatns_lambda_tens

        type(exatns_lambda_tens),intent(inout) :: lambda_t
        type(exatns_comm_tens),intent(inout) :: comm_t
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        logical, intent(in) :: CCD
        integer :: ierr

        ierr=exatns_tensor_init(lambda_t%w_ovoo,'ZERO')
        ierr=exatns_tensor_contract("W(i,e,m,n)+=V+(m,n,i,e)*K()",lambda_t%w_ovoo,int_t%ooov,one_tensor)
        ierr=exatns_tensor_contract("W(i,e,m,n)+=F(i,f)*T(f,e,m,n)",lambda_t%w_ovoo,lambda_t%fbar_ov,comm_t%t2)
!       Need to permute m and n, so term appears twice
        ierr=exatns_tensor_contract("W(i,e,m,n)+=V(i,o,m,f)*T(e,f,n,o)",lambda_t%w_ovoo,int_t%ooov,comm_t%t2)
        ierr=exatns_tensor_contract("W(i,e,m,n)+=V(o,i,n,f)*T(e,f,m,o)",lambda_t%w_ovoo,int_t%ooov,comm_t%t2)
!       Need to permute m and n, so term appears twice
        if (.not. CCD) then
            ierr=exatns_tensor_contract("W(i,e,m,n)+=R(e,i,f,n)*T(f,m)",lambda_t%w_ovoo,lambda_t%wbar_vovo,comm_t%t1,MINUS_ONE)
            ierr=exatns_tensor_contract("W(i,e,m,n)+=R(e,i,f,m)*T(f,n)",lambda_t%w_ovoo,lambda_t%wbar_vovo,comm_t%t1)
            ierr=exatns_tensor_contract("W(i,e,m,n)+=R(o,i,m,n)*T(e,o)",lambda_t%w_ovoo,lambda_t%w_oooo,comm_t%t1)
        end if
        ierr=exatns_tensor_contract("W(i,e,m,n)+=V(e,i,f,g)*T(g,f,m,n)",lambda_t%w_ovoo,int_t%vovv,comm_t%tau,ONE_HALF)

       end subroutine get_w_ovoo_exatns

       subroutine get_w_vovo_exatns(w_vovo,comm_t,oovv_tensor,occ_space_id,occ_space_root, &
                                        vir_space_id,vir_space_root,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens

        type(tens_rcrsv_t),intent(inout) :: w_vovo
        type(exatns_comm_tens),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: oovv_tensor
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        type(tens_rcrsv_t) :: tau2_tensor
        integer(INTD), intent(in) :: occ_space_id, vir_space_id
        integer(INTL), intent(in) :: occ_space_root, vir_space_root
        integer(INTD) :: vvoo_id(4) 
        integer(INTL) :: vvoo_root(4)
        logical, intent(in) :: CCD
        integer(INTD) :: ierr

        !Auxiliary tensor
        vvoo_id   = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
        vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)
        ierr=exatns_tensor_create(tau2_tensor,"tau2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)

        ierr=exatns_tensor_init(tau2_tensor,'ZERO')
        ierr=exatns_tensor_init(w_vovo,'ZERO')

        ierr=exatns_tensor_contract("W(b,m,e,j)+=H(b,m,e,j)*K()",w_vovo,comm_t%h_int,one_tensor)
        ierr=exatns_tensor_contract("T(f,b,j,n)+=R(f,b,j,n)*K()",tau2_tensor,comm_t%t2,one_tensor)
        if (.not. CCD) ierr=exatns_tensor_contract("T(f,b,j,n)+=T(f,j)*T(b,n)",tau2_tensor,comm_t%t1,comm_t%t1)
        ierr=exatns_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(f,b,j,n)",w_vovo,oovv_tensor,tau2_tensor)

        ierr=exatns_tensor_destroy(tau2_tensor)

       end subroutine get_w_vovo_exatns

       subroutine get_wbar_vovo_exatns(wbar_vovo,t2_tensor,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens

        type(tens_rcrsv_t),intent(inout) :: wbar_vovo
        type(tens_rcrsv_t),intent(inout) :: t2_tensor
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(wbar_vovo,'ZERO')
        ierr=exatns_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(b,f,n,j)",wbar_vovo,int_t%oovv,t2_tensor)
        ierr=exatns_tensor_contract("W(b,m,e,j)+=V(b,m,e,j)*K()",wbar_vovo,int_t%vovo,one_tensor)

       end subroutine get_wbar_vovo_exatns

       subroutine get_w_vvvo_exatns(lambda_t,comm_t,int_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens,exatns_lambda_tens

        type(exatns_lambda_tens),intent(inout) :: lambda_t
        type(exatns_comm_tens),intent(inout) :: comm_t
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        logical, intent(in) :: CCD
        integer :: ierr

        ierr=exatns_tensor_init(lambda_t%w_vvvo,'ZERO')
        ierr=exatns_tensor_contract("W(e,f,a,m)+=V+(a,m,e,f)*K()",lambda_t%w_vvvo,int_t%vovv,one_tensor)
!       Need to permute e and f, so term appears twice
        ierr=exatns_tensor_contract("W(e,f,a,m)+=V(e,n,a,g)*T(f,g,m,n)",lambda_t%w_vvvo,int_t%vovv,comm_t%t2)
        ierr=exatns_tensor_contract("W(e,f,a,m)+=V(f,n,a,g)*T(g,e,m,n)",lambda_t%w_vvvo,int_t%vovv,comm_t%t2)
        ierr=exatns_tensor_contract("W(e,f,a,m)+=F(n,a)*T(e,f,m,n)",lambda_t%w_vvvo,lambda_t%fbar_ov,comm_t%t2)
        ierr=exatns_tensor_contract("W(e,f,a,m)+=V(n,o,m,a)*T(e,f,o,n)",lambda_t%w_vvvo,int_t%ooov,comm_t%tau,ONE_HALF)
        if (.not.CCD) then
!           Need to permute e and f, so term appears twice
            ierr=exatns_tensor_contract("W(e,f,a,m)+=R(f,n,a,m)*T(e,n)",lambda_t%w_vvvo,lambda_t%wbar_vovo,comm_t%t1)
            ierr=exatns_tensor_contract("W(e,f,a,m)+=R(e,n,a,m)*T(f,n)",lambda_t%w_vvvo,lambda_t%wbar_vovo,comm_t%t1,MINUS_ONE)
            ierr=exatns_tensor_contract("W(e,f,a,m)+=R(e,f,a,g)*T(g,m)",lambda_t%w_vvvo,lambda_t%w_vvvv,comm_t%t1)
        end if

       end subroutine get_w_vvvo_exatns

       subroutine get_w_vovv_exatns(w_vovv,t1_tensor,int_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens

        type(tens_rcrsv_t),intent(inout) :: w_vovv
        type(tens_rcrsv_t),intent(inout) :: t1_tensor
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        logical, intent(in) :: CCD
        integer :: ierr

        ierr=exatns_tensor_init(w_vovv,'ZERO')
        if (.not.CCD) ierr=exatns_tensor_contract("W(a,m,e,f)+=V(m,n,e,f)*T(a,n)",w_vovv,int_t%oovv,t1_tensor)
        ierr=exatns_tensor_contract("W(a,m,e,f)+=V(a,m,e,f)*K()",w_vovv,int_t%vovv,one_tensor)

       end subroutine get_w_vovv_exatns

       subroutine get_w_vvvv_exatns(w_vvvv,comm_t,int_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens

        type(tens_rcrsv_t),intent(inout) :: w_vvvv
        type(exatns_comm_tens),intent(inout) :: comm_t
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        logical, intent(in) :: CCD
        integer :: ierr

        ierr=exatns_tensor_init(w_vvvv,'ZERO')
!       Permute e and f, so term appears twice
        if (.not. CCD) then
            ierr=exatns_tensor_contract("W(e,f,a,b)+=V(f,m,a,b)*T(e,m)",w_vvvv,int_t%vovv,comm_t%t1)
            ierr=exatns_tensor_contract("W(e,f,a,b)+=V(e,m,a,b)*T(f,m)",w_vvvv,int_t%vovv,comm_t%t1,MINUS_ONE)
        end if
        ierr=exatns_tensor_contract("W(e,f,a,b)+=V(m,n,a,b)*T(e,f,m,n)",w_vvvv,int_t%oovv,comm_t%tau,ONE_HALF)
        ierr=exatns_tensor_contract("W(e,f,a,b)+=V(e,f,a,b)*K()",w_vvvv,int_t%vvvv,one_tensor)

       end subroutine get_w_vvvv_exatns

       subroutine get_goo_lambda_exatns(goo_tensor,l2_tensor,t2_tensor)

        use exatensor

        type(tens_rcrsv_t), intent(inout) :: goo_tensor
        type(tens_rcrsv_t), intent(inout) :: l2_tensor, t2_tensor
        integer :: ierr

        ierr=exatns_tensor_init(goo_tensor,'ZERO')
        ierr=exatns_tensor_contract("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,t2_tensor,ONE_HALF)

       end subroutine get_goo_lambda_exatns

       subroutine get_gvv_lambda_exatns(gvv_tensor,l2_tensor,t2_tensor)

        use exatensor

        type(tens_rcrsv_t), intent(inout) :: gvv_tensor
        type(tens_rcrsv_t), intent(inout) :: l2_tensor, t2_tensor
        integer :: ierr

        ierr=exatns_tensor_init(gvv_tensor,'ZERO')
        ierr=exatns_tensor_contract("G(e,a)+=L(m,n,a,f)*T(e,f,n,m)",gvv_tensor,l2_tensor,t2_tensor,ONE_HALF)

       end subroutine get_gvv_lambda_exatns

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Start CC2 Exatens
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine get_tauCC2_exatns(comm_t,one_tensor) !L

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       make tau ! validated
!L        ierr=exatns_tensor_init(comm_t%tau)
!!        ierr=exatns_tensor_contract("T(a,b,i,j)+=R(a,b,i,j)*K()",comm_t%tau,comm_t%t2,one_tensor)
        ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",comm_t%tau,comm_t%t1,comm_t%t1) ! o(2)v(2)
        ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",comm_t%tau,comm_t%t1,comm_t%t1,MINUS_ONE) ! o(2)v(2)
 
       end subroutine get_tauCC2_exatns

       subroutine get_fbar_vv_CC2_exatns(fbar_vv,comm_t_CC2,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens,exatns_lambda_tens

!        type(exatns_lambda_tens),intent(inout) :: lambda_t_CC2
        type(tens_rcrsv_t),intent(inout) :: fbar_vv

        type(exatns_comm_tens),intent(inout) :: comm_t_CC2
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       \Ft ea =\f ea -\SUMs \f ma \ts em - \SUMd mf \V mefa \ts fm - \undemi
!       \Vmnaf \taup efmn (\tau in CC2)  (Eq. A5)

        ierr=exatns_tensor_init(fbar_vv,'ZERO')
        ierr=exatns_tensor_contract("H(e,a)+=F(e,a)*K()",fbar_vv,comm_t_CC2%fvv,one_tensor)
        ierr=exatns_tensor_contract("H(e,a)+=F(m,a)*T(e,m)",fbar_vv,comm_t_CC2%fov,comm_t_CC2%t1,MINUS_ONE)
        ierr=exatns_tensor_contract("H(e,a)+=V(m,e,f,a)*T(f,m)",fbar_vv,int_t%vovv,comm_t_CC2%t1,ONE)
        ierr=exatns_tensor_contract("H(e,a)+=V(m,n,a,f)*T(e,f,m,n)",fbar_vv,int_t%oovv,comm_t_CC2%tau,MINUS_ONE_HALF)

       end subroutine get_fbar_vv_CC2_exatns

       subroutine get_fbar_oo_CC2_exatns(fbar_oo,comm_t_CC2,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens,exatns_lambda_tens

!        type(exatns_lambda_tens),intent(inout) :: lambda_t_CC2
        type(tens_rcrsv_t),intent(inout) :: fbar_oo
        type(exatns_comm_tens),intent(inout) :: comm_t_CC2
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!      \Ft im =\f im +\SUMs e \f ie \ts em +\SUMd en \V inme \ts en +\ undemi
!      \SUMt mef \V inef \taup efmn  (\tau in CC2)  (Eq. A4)

        ierr=exatns_tensor_init(fbar_oo,'ZERO')
        ierr=exatns_tensor_contract("H(i,m)+=F(i,m)*K()",fbar_oo,comm_t_CC2%foo,one_tensor)
        ierr=exatns_tensor_contract("H(i,m)+=F(i,e)*T(e,m)",fbar_oo,comm_t_CC2%fov,comm_t_CC2%t1,ONE)
        ierr=exatns_tensor_contract("H(i,m)+=V(i,n,m,e)*T(e,n)",fbar_oo,int_t%ooov,comm_t_CC2%t1,ONE)
        ierr=exatns_tensor_contract("H(i,m)+=V(i,n,e,f)*T(e,f,m,n)",fbar_oo,int_t%oovv,comm_t_CC2%tau,ONE_HALF)

       end subroutine get_fbar_oo_CC2_exatns

       subroutine get_fbar_ov_CC2_exatns(lambda_t_CC2,comm_t_CC2,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens,exatns_lambda_tens

        type(exatns_lambda_tens),intent(inout) :: lambda_t_CC2
        type(exatns_comm_tens),intent(inout) :: comm_t_CC2
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!      \Ft me =\f me +\SUMd nf \V mnef \ts fn (Eq. A6)
        ierr=exatns_tensor_init(lambda_t_CC2%fbar_ov,'ZERO')
        ierr=exatns_tensor_contract("H(m,e)+=F(m,e)*K()",lambda_t_CC2%fbar_ov,comm_t_CC2%fov,one_tensor)
        ierr=exatns_tensor_contract("H(m,e)+=V(m,n,e,f)*T(f,n)",lambda_t_CC2%fbar_ov,int_t%oovv,comm_t_CC2%t1,ONE)

       end subroutine get_fbar_ov_CC2_exatns

       subroutine get_w_vovo_CC2_exatns(w_vovo,comm_t_CC2,oovv_tensor,vvoo_id,vvoo_root,one_tensor)

        use exatensor
        use tensor_algebra
        use exacorr_datatypes, only : exatns_comm_tens

        type(tens_rcrsv_t),intent(inout) :: w_vovo
        type(exatns_comm_tens),intent(inout) :: comm_t_CC2
        type(tens_rcrsv_t),intent(inout) :: oovv_tensor
        type(tens_rcrsv_t) :: t1t1_tensor
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer(INTD) :: ierr
        integer(INTD), dimension(4) :: vvoo_id
        integer(INTL), dimension(4) :: vvoo_root
!!!!!!!!!!!!!       Auxiliary tensor
         ierr=exatns_tensor_create(t1t1_tensor,"t1t1",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)


!        \W bmje = \V bmje +\SUMs f \V mbef \ts fj - \SUMs n \V mnej \ts bn +
!                  +\SUMd nf \V mnef (\td fbjn - \ts fj \ts bn)

        ierr=exatns_tensor_contract("W(b,m,e,j)+=H(b,m,e,j)*K()",w_vovo,comm_t_CC2%h_int,one_tensor)
        ierr=exatns_tensor_init(t1t1_tensor,'ZERO')
        ierr=exatns_tensor_contract("T(f,b,j,n)+=T(f,j)*T(b,n)",t1t1_tensor,comm_t_CC2%t1,comm_t_CC2%t1,MINUS_ONE)
!ReDo these Lines for Check CC2=CCSD if t2 \neg 0 - That's why t1t1 is negative        !
        ierr=exatns_tensor_contract("T(f,b,j,n)+=R(f,b,j,n)*K()",t1t1_tensor,comm_t_CC2%t2,one_tensor)
!End of
        ierr=exatns_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(f,b,j,n)",w_vovo,oovv_tensor,t1t1_tensor)

       end subroutine get_w_vovo_CC2_exatns

       subroutine get_w_oooo_CC2_exatns(w_oooo_CC2,comm_t_CC2,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens

        type(tens_rcrsv_t),intent(inout) :: w_oooo_CC2
        type(exatns_comm_tens),intent(inout) :: comm_t_CC2
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       oooo ! validated
!        ierr=exatns_tensor_init(w_oooo_CC2)
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,j)*K()",w_oooo_CC2,int_t%oooo,one_tensor)
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",w_oooo_CC2,int_t%ooov,comm_t_CC2%t1) ! o(4)v(1)
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",w_oooo_CC2,int_t%ooov,comm_t_CC2%t1, &
                                    MINUS_ONE) ! o(4)v(1)
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",w_oooo_CC2,int_t%oovv,comm_t_CC2%tau, &
                                    ONE_HALF)

       end subroutine get_w_oooo_CC2_exatns

       subroutine get_w_vvvv_CC2_exatns(w_vvvv_CC2,comm_t_CC2,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens
        
        type(tens_rcrsv_t),intent(inout) :: w_vvvv_CC2
        type(exatns_comm_tens),intent(inout) :: comm_t_CC2
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

! \W efab=\V efab -\Pm ef \Vmfab \ts em + \undemi \V mnab \taup efmn
!        ierr=exatns_tensor_init(w_vvvv_CC2)
        ierr=exatns_tensor_contract("W(e,f,a,b)+=V(e,f,a,b)*K()",w_vvvv_CC2,int_t%vvvv,one_tensor)

        ierr=exatns_tensor_contract("W(e,f,a,b)+=V(f,m,a,b)*T(e,m)",w_vvvv_CC2,int_t%vovv,comm_t_CC2%t1)
        ierr=exatns_tensor_contract("W(e,f,a,b)+=V(e,m,a,b)*T(f,m)",w_vvvv_CC2,int_t%vovv,comm_t_CC2%t1,MINUS_ONE)

        ierr=exatns_tensor_contract("W(e,f,a,b)+=V(m,n,a,b)*T(e,f,m,n)",w_vvvv_CC2,int_t%oovv,comm_t_CC2%tau,ONE_HALF)

       end subroutine get_w_vvvv_CC2_exatns
 
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  End CC2 Exatens
!!!!!!!!!!!!!!!!!!!!!!!

#endif

       subroutine get_tau_talsh(comm_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_comm_tens),intent(inout) :: comm_t
        logical, intent(in) :: CCD
        integer :: ierr

!       tau ! validated
        ierr=talsh_tensor_init(comm_t%tau)
        ierr=talsh_tensor_add("T(a,b,i,j)+=R(a,b,i,j)",comm_t%tau,comm_t%t2)
        if (.not.CCD) then
            ierr=talsh_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",comm_t%tau,comm_t%t1,comm_t%t1)                  ! o(2)v(2)
            ierr=talsh_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",comm_t%tau,comm_t%t1,comm_t%t1,scale=MINUS_ONE)  ! o(2)v(2)
        end if

       end subroutine get_tau_talsh


       subroutine get_hoo_talsh(comm_t,oovv_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: oovv_tensor
        integer :: ierr

!       oo ! validated
        ierr=talsh_tensor_init(comm_t%hoo)
        ierr=talsh_tensor_add("H(k,i)+=F(k,i)",comm_t%hoo,comm_t%foo,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("H(k,i)+=V(k,l,c,d)*T(c,d,l,i)",comm_t%hoo,oovv_tensor,comm_t%tau,scale=ONE_HALF) ! o(3)v(2)

       end subroutine get_hoo_talsh

       subroutine get_hov_talsh(comm_t,oovv_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: oovv_tensor
        integer :: ierr

!       ov ! validated
        ierr=talsh_tensor_init(comm_t%hov)
        ierr=talsh_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",comm_t%hov,oovv_tensor,comm_t%t1) ! o(2)v(2)
        ierr=talsh_tensor_add("H(k,c)+=F(k,c)",comm_t%hov,comm_t%fov)

       end subroutine get_hov_talsh

       subroutine get_hvv_talsh(comm_t,oovv_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: oovv_tensor
        integer :: ierr

!       vv ! validated
        ierr=talsh_tensor_init(comm_t%hvv)
        ierr=talsh_tensor_add("H(a,c)+=F(a,c)",comm_t%hvv,comm_t%fvv)
        ierr=talsh_tensor_contract("H(a,c)+=V(k,l,c,d)*T(d,a,k,l)",comm_t%hvv,oovv_tensor,comm_t%tau,scale=ONE_HALF) ! o(2)v(3)

       end subroutine get_hvv_talsh

       subroutine get_goo_talsh(comm_t,ooov_tensor,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: ooov_tensor
        logical, intent(in) :: CCD
        integer :: ierr

!       oo ! validated
        ierr=talsh_tensor_init(comm_t%goo)
        ierr=talsh_tensor_add("G(k,i)+=H(k,i)",comm_t%goo,comm_t%hoo)
        if (.not.CCD) then
            ierr=talsh_tensor_contract("G(k,i)+=V(k,l,i,c)*T(c,l)",comm_t%goo,ooov_tensor,comm_t%t1,scale=MINUS_ONE)     ! o(3)v(1)
            ierr=talsh_tensor_contract("G(k,i)+=F(k,c)*T(c,i)",comm_t%goo,comm_t%fov,comm_t%t1,scale=MINUS_ONE)                         ! o(2)v(1)
        end if

       end subroutine get_goo_talsh

       subroutine get_gvv_talsh(comm_t,vovv_tensor,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: vovv_tensor
        logical, intent(in) :: CCD
        integer :: ierr

!       vv ! validated
        ierr=talsh_tensor_init(comm_t%gvv)
        ierr=talsh_tensor_add("G(a,c)+=H(a,c)",comm_t%gvv,comm_t%hvv)
        if (.not.CCD) then
            ierr=talsh_tensor_contract("G(a,c)+=V(a,k,c,d)*T(d,k)",comm_t%gvv,vovv_tensor,comm_t%t1)             ! o(1)v(3)
            ierr=talsh_tensor_contract("G(a,c)+=F(k,c)*T(a,k)",comm_t%gvv,comm_t%fov,comm_t%t1,scale=MINUS_ONE) ! o(1)v(2)
        end if

       end subroutine get_gvv_talsh

       subroutine get_a_talsh(comm_t,int_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens

        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_intg_tens),intent(inout) :: int_t
        logical, intent(in) :: CCD
        integer :: ierr

!       oooo ! validated
        ierr=talsh_tensor_init(comm_t%a_int)
        ierr=talsh_tensor_add("A(k,l,i,j)+=V(k,l,i,j)",comm_t%a_int,int_t%oooo)
        if (.not.CCD) then
            ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",comm_t%a_int,int_t%ooov,comm_t%t1) ! o(4)v(1)
            ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",comm_t%a_int,int_t%ooov,comm_t%t1, &
                                    scale=MINUS_ONE) ! o(4)v(1)
        end if
        ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",comm_t%a_int,int_t%oovv,comm_t%tau, &
                                    scale=ONE_HALF)

       end subroutine get_a_talsh

       subroutine get_h_talsh(comm_t,int_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens

        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_intg_tens),intent(inout) :: int_t
        logical, intent(in) :: CCD
        integer :: ierr

!       vovo ! more tests
        ierr=talsh_tensor_init(comm_t%h_int)
        if (.not.CCD) then
            ierr=talsh_tensor_contract("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",comm_t%h_int,int_t%vovv,comm_t%t1)                 ! o(2)v(3)
            ierr=talsh_tensor_contract("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",comm_t%h_int,int_t%ooov,comm_t%t1,scale=MINUS_ONE) ! o(3)v(2)
        end if
        ierr=talsh_tensor_add("H(a,k,c,i)+=V(a,k,c,i)",comm_t%h_int,int_t%vovo)

       end subroutine get_h_talsh

       subroutine get_fbar_oo_talsh(fbar_oo,comm_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_tens_t),intent(inout)   :: fbar_oo
        type(talsh_comm_tens),intent(inout)   :: comm_t
        logical, intent(in) :: CCD
        integer :: ierr

        ierr=talsh_tensor_add("F(i,m)+=G(i,m)",fbar_oo,comm_t%goo,scale=MINUS_ONE)
        if (.not.CCD) ierr=talsh_tensor_contract("F(i,m)+=R(i,e)*T(e,m)",fbar_oo,comm_t%fov,comm_t%t1,scale=TWO)

       end subroutine get_fbar_oo_talsh

       subroutine get_fbar_ov_talsh(fbar_ov,hov_tensor)

        use talsh

        type(talsh_tens_t),intent(inout)   :: fbar_ov
        type(talsh_tens_t),intent(inout)   :: hov_tensor
        integer :: ierr

        ierr=talsh_tensor_add("F(m,e)+=H(m,e)",fbar_ov,hov_tensor)

       end subroutine get_fbar_ov_talsh

       subroutine get_fbar_vv_talsh(fbar_vv,gvv_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_tens_t),intent(inout)   :: fbar_vv
        type(talsh_tens_t),intent(inout)   :: gvv_tensor
        integer :: ierr

        ierr=talsh_tensor_add("F(e,a)+=G(e,a)",fbar_vv,gvv_tensor)

       end subroutine get_fbar_vv_talsh

       subroutine get_w_oooo_talsh(w_oooo,a_int)

        use talsh

        type(talsh_tens_t),intent(inout) :: w_oooo
        type(talsh_tens_t),intent(inout) :: a_int
        integer :: ierr

        ierr=talsh_tensor_add("W(i,j,m,n)+=A(i,j,m,n)",w_oooo,a_int)

       end subroutine get_w_oooo_talsh

       subroutine get_w_ooov_talsh(w_ooov,t1_tensor,int_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens

        type(talsh_tens_t),intent(inout)   :: w_ooov
        type(talsh_intg_tens),intent(inout)   :: int_t
        type(talsh_tens_t),intent(inout) :: t1_tensor
        logical, intent(in) :: CCD
        integer :: ierr

        if (.not.CCD) ierr=talsh_tensor_contract("W(m,n,i,e)+=T(f,i)*V(m,n,f,e)",w_ooov,t1_tensor,int_t%oovv)
        ierr=talsh_tensor_add("W(m,n,i,e)+=V(m,n,i,e)",w_ooov,int_t%ooov)
       
       end subroutine get_w_ooov_talsh

       subroutine get_w_ovoo_talsh(lambda_t,comm_t,int_t,CCD,one_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens,talsh_lambda_tens

        type(talsh_lambda_tens),intent(inout) :: lambda_t
        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_intg_tens),intent(inout) :: int_t
        logical, intent(in) :: CCD
        type(talsh_tens_t) :: one_tensor
        integer :: ierr

        ierr=talsh_tensor_contract("W(i,e,m,n)+=V+(m,n,i,e)*K()",lambda_t%w_ovoo,int_t%ooov,one_tensor)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=F(i,f)*T(f,e,m,n)",lambda_t%w_ovoo,lambda_t%fbar_ov,comm_t%t2)
!       Need to permute m and n, so term appears twice
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(i,o,m,f)*T(e,f,n,o)",lambda_t%w_ovoo,int_t%ooov,comm_t%t2)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(o,i,n,f)*T(e,f,m,o)",lambda_t%w_ovoo,int_t%ooov,comm_t%t2)
        if (.not.CCD) then 
!        Need to permute m and n, so term appears twice
         ierr=talsh_tensor_contract("W(i,e,m,n)+=R(e,i,f,n)*T(f,m)",lambda_t%w_ovoo,lambda_t%wbar_vovo,comm_t%t1,scale=MINUS_ONE)
         ierr=talsh_tensor_contract("W(i,e,m,n)+=R(e,i,f,m)*T(f,n)",lambda_t%w_ovoo,lambda_t%wbar_vovo,comm_t%t1)
         ierr=talsh_tensor_contract("W(i,e,m,n)+=R(o,i,m,n)*T(e,o)",lambda_t%w_ovoo,lambda_t%w_oooo,comm_t%t1)
        end if
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(e,i,f,g)*T(g,f,m,n)",lambda_t%w_ovoo,int_t%vovv,comm_t%tau,scale=ONE_HALF)

       end subroutine get_w_ovoo_talsh

       subroutine get_w_vovo_talsh(w_vovo,comm_t,oovv_tensor,nocc,nvir,CCD)

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_tens_t),intent(inout) :: w_vovo
        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: oovv_tensor
        integer,intent(in) :: nocc, nvir
        logical, intent(in) :: CCD
        type(talsh_tens_t) :: tau2_tensor
        integer(INTD) :: vvoo_dims(4)
        integer(INTD) :: ierr

!       Auxiliary tensor
        vvoo_dims(1:2) = nvir
        vvoo_dims(3:4) = nocc
        ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)

        ierr=talsh_tensor_add("W(b,m,e,j)+=H(b,m,e,j)",w_vovo,comm_t%h_int)
        ierr=talsh_tensor_init(tau2_tensor)
        ierr=talsh_tensor_add("T(f,b,j,n)+=R(f,b,j,n)",tau2_tensor,comm_t%t2)
        if (.not.CCD) ierr=talsh_tensor_contract("T(f,b,j,n)+=T(f,j)*T(b,n)",tau2_tensor,comm_t%t1,comm_t%t1)
        ierr=talsh_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(f,b,j,n)",w_vovo,oovv_tensor,tau2_tensor)

       end subroutine get_w_vovo_talsh

       subroutine get_wbar_vovo_talsh(wbar_vovo,t2_tensor,int_t)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens

        type(talsh_tens_t),intent(inout) :: wbar_vovo
        type(talsh_tens_t),intent(inout) :: t2_tensor
        type(talsh_intg_tens),intent(inout) :: int_t
        integer :: ierr

        ierr=talsh_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(b,f,n,j)",wbar_vovo,int_t%oovv,t2_tensor)
        ierr=talsh_tensor_add("W(b,m,e,j)+=V(b,m,e,j)",wbar_vovo,int_t%vovo)

       end subroutine get_wbar_vovo_talsh

       subroutine get_w_vvvo_talsh(lambda_t,comm_t,int_t,CCD,one_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens,talsh_lambda_tens

        type(talsh_lambda_tens),intent(inout) :: lambda_t
        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_intg_tens),intent(inout) :: int_t
        logical, intent(in) :: CCD
        type(talsh_tens_t) :: one_tensor
        integer :: ierr

        ierr=talsh_tensor_contract("W(e,f,a,m)+=V+(a,m,e,f)*K()",lambda_t%w_vvvo,int_t%vovv,one_tensor)
!       Need to permute e and f, so term appears twice
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(e,n,a,g)*T(f,g,m,n)",lambda_t%w_vvvo,int_t%vovv,comm_t%t2)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(f,n,a,g)*T(g,e,m,n)",lambda_t%w_vvvo,int_t%vovv,comm_t%t2)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=F(n,a)*T(e,f,m,n)",lambda_t%w_vvvo,lambda_t%fbar_ov,comm_t%t2)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(n,o,m,a)*T(e,f,o,n)",lambda_t%w_vvvo,int_t%ooov,comm_t%tau,scale=ONE_HALF)
        if (.not.CCD) then
!        Need to permute e and f, so term appears twice
         ierr=talsh_tensor_contract("W(e,f,a,m)+=R(f,n,a,m)*T(e,n)",lambda_t%w_vvvo,lambda_t%wbar_vovo,comm_t%t1)
         ierr=talsh_tensor_contract("W(e,f,a,m)+=R(e,n,a,m)*T(f,n)",lambda_t%w_vvvo,lambda_t%wbar_vovo,comm_t%t1,scale=MINUS_ONE)
         ierr=talsh_tensor_contract("W(e,f,a,m)+=R(e,f,a,g)*T(g,m)",lambda_t%w_vvvo,lambda_t%w_vvvv,comm_t%t1)
        end if

        end subroutine get_w_vvvo_talsh

       subroutine get_w_vovv_talsh(w_vovv,t1_tensor,int_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens
        
        type(talsh_tens_t),intent(inout) :: w_vovv
        type(talsh_tens_t),intent(inout) :: t1_tensor
        type(talsh_intg_tens),intent(inout) :: int_t
        logical, intent(in) :: CCD
        integer :: ierr

        if (.not.CCD) ierr=talsh_tensor_contract("W(a,m,e,f)+=V(m,n,e,f)*T(a,n)",w_vovv,int_t%oovv,t1_tensor)
        ierr=talsh_tensor_add("W(a,m,e,f)+=V(a,m,e,f)",w_vovv,int_t%vovv)

       end subroutine get_w_vovv_talsh

       subroutine get_w_vvvv_talsh(w_vvvv,comm_t,int_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens
        
        type(talsh_tens_t),intent(inout) :: w_vvvv
        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_intg_tens),intent(inout) :: int_t
        logical, intent(in) :: CCD
        integer :: ierr

        if (.not.CCD) then
!        Permute e and f, so term appears twice
         ierr=talsh_tensor_contract("W(e,f,a,b)+=V(f,m,a,b)*T(e,m)",w_vvvv,int_t%vovv,comm_t%t1)
         ierr=talsh_tensor_contract("W(e,f,a,b)+=V(e,m,a,b)*T(f,m)",w_vvvv,int_t%vovv,comm_t%t1,scale=MINUS_ONE)
        end if
        ierr=talsh_tensor_contract("W(e,f,a,b)+=V(m,n,a,b)*T(e,f,m,n)",w_vvvv,int_t%oovv,comm_t%tau,scale=ONE_HALF)
        ierr=talsh_tensor_add("W(e,f,a,b)+=V(e,f,a,b)",w_vvvv,int_t%vvvv)

       end subroutine get_w_vvvv_talsh

       subroutine get_goo_lambda_talsh(goo_tensor,l2_tensor,t2_tensor)

        use talsh

        type(talsh_tens_t), intent(inout) :: goo_tensor
        type(talsh_tens_t), intent(inout) :: l2_tensor, t2_tensor
        integer :: ierr

        ierr=talsh_tensor_init(goo_tensor)
        ierr=talsh_tensor_contract("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,t2_tensor,scale=ONE_HALF)

       end subroutine get_goo_lambda_talsh

       subroutine get_gvv_lambda_talsh(gvv_tensor,l2_tensor,t2_tensor)

        use talsh

        type(talsh_tens_t), intent(inout) :: gvv_tensor
        type(talsh_tens_t), intent(inout) :: l2_tensor, t2_tensor
        integer :: ierr

        ierr=talsh_tensor_init(gvv_tensor)
        ierr=talsh_tensor_contract("G(e,a)+=L(m,n,a,f)*T(e,f,n,m)",gvv_tensor,l2_tensor,t2_tensor,scale=ONE_HALF)

       end subroutine get_gvv_lambda_talsh

!!!!!!!
!!!!!!!    Start CC2 Talsh
!!!!!!!
       subroutine get_tauCC2_talsh(comm_t_CC2)  !L

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_comm_tens),intent(inout) :: comm_t_CC2
        integer :: ierr

        logical :: shunt_CC2=.false.

!       tau ! validated
        ierr=talsh_tensor_init(comm_t_CC2%tau)


!       Here I do a 'usual tau' to check
        if(shunt_CC2) then
                ierr=talsh_tensor_add("T(a,b,i,j)=H(a,b,i,j)",comm_t_CC2%tau,comm_t_CC2%t2)
                print*, " "
                print*, " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING SHUNT IN TAU_CC2  "
                print*, " "
                print*, "t2 comm CC2 - Subroutine   get_tauCC2 access to t2   "
                print*, " "
        end if

        ierr=talsh_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",comm_t_CC2%tau,comm_t_CC2%t1,comm_t_CC2%t1)                  ! o(2)v(2)
        ierr=talsh_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",comm_t_CC2%tau,comm_t_CC2%t1,comm_t_CC2%t1,scale=MINUS_ONE)  ! o(2)v(2)

       end subroutine get_tauCC2_talsh

       subroutine get_fbar_oo_CC2_talsh(fbar_oo,comm_t_CC2,int_t)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens,talsh_lambda_tens
      
        type(talsh_tens_t),intent(inout)   :: fbar_oo
!        type(talsh_lambda_tens),intent(inout) :: lambda_t_CC2
        type(talsh_comm_tens),intent(inout) :: comm_t_CC2
        type(talsh_intg_tens),intent(inout) :: int_t
        integer :: ierr

!      \Ft im =\f im +\SUMs e \f ie \ts em +\SUMd en \V inme \ts en +\ undemi
!      \SUMt mef \V inef \taup efmn  (\tau in CC2)  (Eq. A4)

        ierr=talsh_tensor_init(fbar_oo)
        ierr=talsh_tensor_add("H(i,m)+=F(i,m)",fbar_oo,comm_t_CC2%foo)
        ierr=talsh_tensor_contract("H(i,m)+=F(i,e)*T(e,m)",fbar_oo,comm_t_CC2%fov,comm_t_CC2%t1,scale=ONE)
        ierr=talsh_tensor_contract("H(i,m)+=V(i,n,m,e)*T(e,n)",fbar_oo,int_t%ooov,comm_t_CC2%t1,scale=ONE)
        ierr=talsh_tensor_contract("H(i,m)+=V(i,n,e,f)*T(e,f,m,n)",fbar_oo,int_t%oovv,comm_t_CC2%tau,scale=ONE_HALF)

       end subroutine get_fbar_oo_CC2_talsh

       subroutine get_fbar_vv_CC2_talsh(fbar_vv,comm_t_CC2,int_t)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens,talsh_lambda_tens

!        type(talsh_lambda_tens),intent(inout) :: lambda_t_CC2
        type(talsh_tens_t),intent(inout)   :: fbar_vv
        type(talsh_comm_tens),intent(inout) :: comm_t_CC2
        type(talsh_intg_tens),intent(inout) :: int_t
        integer :: ierr

!       \Ft ea =\f ea -\SUMs \f ma \ts em - \SUMd mf \V mefa \ts fm - \undemi
!       \Vmnaf \taup efmn (\tau in CC2)  (Eq. A5)

        ierr=talsh_tensor_init(fbar_vv)
        ierr=talsh_tensor_add("H(e,a)+=F(e,a)",fbar_vv,comm_t_CC2%fvv)
        ierr=talsh_tensor_contract("H(e,a)+=F(m,a)*T(e,m)",fbar_vv,comm_t_CC2%fov,comm_t_CC2%t1,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("H(e,a)+=V(e,m,f,a)*T(f,m)",fbar_vv,int_t%vovv,comm_t_CC2%t1,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("H(e,a)+=V(m,n,a,f)*T(e,f,m,n)",fbar_vv,int_t%oovv,comm_t_CC2%tau,scale=MINUS_ONE_HALF)

       end subroutine get_fbar_vv_CC2_talsh

       subroutine get_fbar_ov_CC2_talsh(fbar_ov,comm_t_CC2,int_t)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens,talsh_lambda_tens

        type(talsh_tens_t),intent(inout)   :: fbar_ov
!        type(talsh_lambda_tens),intent(inout) :: lambda_t_CC2
        type(talsh_comm_tens),intent(inout) :: comm_t_CC2
        type(talsh_intg_tens),intent(inout) :: int_t
        integer :: ierr

!      \Ft me =\f me +\SUMd nf \V mnef \ts fn (Eq. A6)
        ierr=talsh_tensor_init(fbar_ov)
        ierr=talsh_tensor_add("H(m,e)+=F(m,e)",fbar_ov,comm_t_CC2%fov)
        ierr=talsh_tensor_contract("H(m,e)+=V(m,n,e,f)*T(f,n)",fbar_ov,int_t%oovv,comm_t_CC2%t1,scale=ONE)

       end subroutine get_fbar_ov_CC2_talsh

       subroutine get_w_vovo_CC2_talsh(w_vovo,comm_t_CC2,oovv_tensor,nocc,nvir)

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_comm_tens

        type(talsh_tens_t),intent(inout) :: w_vovo
        type(talsh_comm_tens),intent(inout) :: comm_t_CC2
        type(talsh_tens_t),intent(inout) :: oovv_tensor
        integer,intent(in) :: nocc, nvir
        type(talsh_tens_t) :: t1t1_aux_tensor
        integer(INTD) :: vvoo_dims(4)
        integer(INTD) :: ierr

        logical :: shunt_CC2 =.false.

!       Auxiliary tensor
        vvoo_dims(1:2) = nvir
        vvoo_dims(3:4) = nocc
        ierr=talsh_tensor_construct(t1t1_aux_tensor,C8,vvoo_dims,init_val=ZERO)

!        \W bmje = \V bmje +\SUMs f \V mbef \ts fj - \SUMs n \V mnej \ts bn +
!                  +\SUMd nf \V mnef (\td fbjn - \ts fj \ts bn)

        ierr=talsh_tensor_add("W(b,m,e,j)+=H(b,m,e,j)",w_vovo,comm_t_CC2%h_int)
        ierr=talsh_tensor_init(t1t1_aux_tensor)
        ierr=talsh_tensor_contract("T(f,b,j,n)+=T(f,j)*T(b,n)",t1t1_aux_tensor,comm_t_CC2%t1,comm_t_CC2%t1,scale=ONE)


!ReDo these Lines for Check CC2=CCSD if t2 \neg 0 - That's why t1t1 is negative        !
        if(shunt_CC2) then
                print*, " "
                print*, " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING SHUNT IN W VOVO_CC2  "
                print*, " "
                ierr=talsh_tensor_add("T(f,b,j,n)+=R(f,b,j,n)",t1t1_aux_tensor,comm_t_CC2%t2,scale=ONE)
        end if
!End of

        ierr=talsh_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(f,b,j,n)",w_vovo,oovv_tensor,t1t1_aux_tensor,scale=ONE)



       end subroutine get_w_vovo_CC2_talsh

       subroutine get_w_oooo_CC2_talsh(w_oooo_CC2,comm_t_CC2,int_t)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens

        type(talsh_comm_tens),intent(inout) :: comm_t_CC2
        type(talsh_tens_t),intent(inout) :: w_oooo_CC2
        type(talsh_intg_tens),intent(inout) :: int_t
        integer :: ierr

!       \W mnij = \V mnij + \Pm mn \SUMs e \V ijen \ts em + \undemi \SUMd ef \V ijef \taup efmn

        ierr=talsh_tensor_add("A(k,l,i,j)+=V(k,l,i,j)",w_oooo_CC2,int_t%oooo)
        ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",w_oooo_CC2,int_t%ooov,comm_t_CC2%t1) ! o(4)v(1)
        ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",w_oooo_CC2,int_t%ooov,comm_t_CC2%t1, &
                                    scale=MINUS_ONE) ! o(4)v(1)
        ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",w_oooo_CC2,int_t%oovv,comm_t_CC2%tau, &
                                    scale=ONE_HALF)

       end subroutine get_w_oooo_CC2_talsh

       subroutine get_w_vvvv_CC2_talsh(w_vvvv_CC2,comm_t_CC2,int_t)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens
        
        type(talsh_tens_t),intent(inout) :: w_vvvv_CC2
        type(talsh_comm_tens),intent(inout) :: comm_t_CC2
        type(talsh_intg_tens),intent(inout) :: int_t
        integer :: ierr

! \W efab=\V efab -\Pm ef \Vmfab \ts em + \undemi \V mnab \taup efmn
        ierr=talsh_tensor_add("W(e,f,a,b)+=V(e,f,a,b)",w_vvvv_CC2,int_t%vvvv)

        ierr=talsh_tensor_contract("W(e,f,a,b)+=V(f,m,a,b)*T(e,m)",w_vvvv_CC2,int_t%vovv,comm_t_CC2%t1)
        ierr=talsh_tensor_contract("W(e,f,a,b)+=V(e,m,a,b)*T(f,m)",w_vvvv_CC2,int_t%vovv,comm_t_CC2%t1,scale=MINUS_ONE)

        ierr=talsh_tensor_contract("W(e,f,a,b)+=V(m,n,a,b)*T(e,f,m,n)",w_vvvv_CC2,int_t%oovv,comm_t_CC2%tau,scale=ONE_HALF)

       end subroutine get_w_vvvv_CC2_talsh
       


end module
