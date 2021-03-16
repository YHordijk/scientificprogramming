module exacorr_gradient

! This module contains routines for solving the Lambda equations

#if defined (VAR_MPI)

  use exatensor
  use exacorr_ao_to_mo,  only : print_tensornorm2
  use exacorr_datatypes
  use exacorr_utils,     only : print_date
  use, intrinsic:: ISO_C_BINDING

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                           MINUS_ONE=(-1.D0,0.D0),ONE=(1.0D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0)

  private

  public exacorr_dm_driver

  contains

    subroutine exacorr_dm_driver(exa_input,dims,t2_tensor,l2_tensor,t1_tensor,l1_tensor)

!     Routine for obtaining the density matrices in the MO basis

!     input variables
      type(exacc_input), intent(in  )  :: exa_input
      type(space_dims),  intent(in   ) :: dims
      type(tens_rcrsv_t),intent(inout) :: t2_tensor, l2_tensor
      type(tens_rcrsv_t),intent(inout), optional :: t1_tensor, l1_tensor
!     density matrices type, output
      type(exatns_dm_tens) :: gamma
!     tensor dimensions
      integer(INTD), dimension(2) :: oo_id,  ov_id,  vo_id,  vv_id
      integer(INTL), dimension(2) :: oo_root,ov_root,vo_root,vv_root
!     CCD switch
      logical :: CCD
!     error code
      integer(INTD) :: ierr

      CCD = exa_input%CCD

!**********************
!     set ids + roots *
!**********************

      oo_id   = (/dims%occ_space_id,   dims%occ_space_id  /)
      oo_root = (/dims%occ_space_root, dims%occ_space_root/)
      vv_id   = (/dims%vir_space_id,   dims%vir_space_id  /)
      vv_root = (/dims%vir_space_root, dims%vir_space_root/)
      if (.not.CCD) then
        ov_id   = (/dims%occ_space_id,   dims%vir_space_id  /)
        vo_id   = (/dims%vir_space_id,   dims%occ_space_id  /)
        ov_root = (/dims%occ_space_root, dims%vir_space_root/)
        vo_root = (/dims%vir_space_root, dims%occ_space_root/)
      end if

!************************************
!     construct gamma prime tensors *
!************************************

      ierr=exatns_tensor_create(gamma%prime_oo,"gprime_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_init(gamma%prime_oo,ZERO)
      if (.not.CCD) then
        ierr=exatns_tensor_create(gamma%prime_ov,"gprime_ov",ov_id,ov_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_create(gamma%symm_ov,"gsymm_ov",ov_id,ov_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_create(gamma%prime_vo,"gprime_vo",vo_id,vo_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_create(gamma%symm_vo,"gsymm_vo",vo_id,vo_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(gamma%prime_ov,ZERO)
        ierr=exatns_tensor_init(gamma%symm_ov, ZERO)
        ierr=exatns_tensor_init(gamma%prime_vo,ZERO)
        ierr=exatns_tensor_init(gamma%symm_vo, ZERO)
      end if
      ierr=exatns_tensor_create(gamma%prime_vv,"gprime_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_oo,"gsymm_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_vv,"gsymm_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_init(gamma%prime_vv,ZERO)
      ierr=exatns_tensor_init(gamma%symm_oo,ZERO)
      ierr=exatns_tensor_init(gamma%symm_vv,ZERO)

!**************************************
!     calculate DM blocks in MO basis *
!**************************************

      if (.not.CCD) then
        call calc_dm_mo(gamma,t2_tensor,l2_tensor,dims,CCD,exa_input%print_level,t1_tensor,l1_tensor)
      else
        call calc_dm_mo(gamma,t2_tensor,l2_tensor,dims,CCD,exa_input%print_level)
      end if
      call print_date("Calculated density matrices in MO basis")

!***********************
!     write DM to file *
!***********************

      call store_dm_mo(exa_input,gamma)

!     cleanup
      ierr=exatns_tensor_destroy(gamma%prime_oo)
      if (.not.CCD) then
        ierr=exatns_tensor_destroy(gamma%prime_ov)
        ierr=exatns_tensor_destroy(gamma%symm_ov)
        ierr=exatns_tensor_destroy(gamma%prime_vo)
        ierr=exatns_tensor_destroy(gamma%symm_vo)
      end if
      ierr=exatns_tensor_destroy(gamma%prime_vv)
      ierr=exatns_tensor_destroy(gamma%symm_oo)
      ierr=exatns_tensor_destroy(gamma%symm_vv)
       
      if (ierr.ne.EXA_SUCCESS) then
        print*," ccdriver: error in evaluation of density matrices, code=",ierr
        return
      end if

    end subroutine exacorr_dm_driver

    subroutine calc_dm_mo(gamma,t2_tensor,l2_tensor,dims,CCD,print_level,t1_tensor,l1_tensor)

!     input variables
      type(exatns_dm_tens), intent(inout) :: gamma
      type(tens_rcrsv_t),   intent(inout) :: t2_tensor, l2_tensor
      type(tens_rcrsv_t),   intent(inout), optional :: t1_tensor, l1_tensor
      type(space_dims),     intent(in   ) :: dims
      logical,              intent(in   ) :: CCD
      integer,              intent(in   ) :: print_level
!     auxiliary tensors
      type(tens_rcrsv_t) :: one_dm_tensor
      type(tens_rcrsv_t) :: tau2_tensor
      type(tens_rcrsv_t) :: vv_aux_tensor,oo_aux_tensor
!     tensor dimensions
      integer(INTD), dimension(2) :: oo_id,   vv_id
      integer(INTL), dimension(2) :: oo_root, vv_root
      integer(INTD), dimension(4) :: vvoo_id
      integer(INTL), dimension(4) :: vvoo_root
!     error code
      integer(INTD) :: ierr

!*********************************
!     construct auxiliary tensor *
!*********************************

      ierr=exatns_tensor_create(one_dm_tensor,"one_dm",EXA_DATA_KIND_C8)
      ierr=exatns_tensor_init(one_dm_tensor,ONE)
      if (.not.CCD) then
        oo_id   = (/dims%occ_space_id,   dims%occ_space_id  /)
        oo_root = (/dims%occ_space_root, dims%occ_space_root/)
        vv_id   = (/dims%vir_space_id,   dims%vir_space_id  /)
        vv_root = (/dims%vir_space_root, dims%vir_space_root/)
        vvoo_id   = (/dims%vir_space_id,   dims%vir_space_id,   dims%occ_space_id,   dims%occ_space_id  /)
        vvoo_root = (/dims%vir_space_root, dims%vir_space_root, dims%occ_space_root, dims%occ_space_root/)

        ierr=exatns_tensor_create(tau2_tensor,"tau2_dm",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_create(vv_aux_tensor,"vv_aux",vv_id,vv_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_create(oo_aux_tensor,"oo_aux",oo_id,oo_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(tau2_tensor,ZERO)
        ierr=exatns_tensor_init(vv_aux_tensor,ZERO)
        ierr=exatns_tensor_init(oo_aux_tensor,ZERO)
      end if

!********************
!     calculate DMs *
!********************

!----------------------------------------------------------------
!     oo: g'(j,i) = - 1/2 t(ef,im) * l(jm,ef) - t(e,i) * l(j,e) |
!----------------------------------------------------------------
      ierr=exatns_tensor_init(gamma%prime_oo,ZERO)
      ierr=exatns_tensor_contract("G(j,i)+=T(e,f,i,m)*L(m,j,e,f)",gamma%prime_oo,t2_tensor,l2_tensor,ONE_HALF)
      if (.not.CCD) ierr=exatns_tensor_contract("G(j,i)+=T(e,i)*L(j,e)",gamma%prime_oo,t1_tensor,l1_tensor,MINUS_ONE)

!--------------------------------------------------------------
!     vv: g'(a,b) = 1/2 t(ae,mn) * l(mn,be) + t(a,m) * l(m,b) |
!--------------------------------------------------------------
      ierr=exatns_tensor_init(gamma%prime_vv,ZERO)
      ierr=exatns_tensor_contract("G(a,b)+=T(a,e,m,n)*L(m,n,b,e)",gamma%prime_vv,t2_tensor,l2_tensor,ONE_HALF)
      if (.not.CCD) ierr=exatns_tensor_contract("G(a,b)+=T(a,m)*L(m,b)",gamma%prime_vv,t1_tensor,l1_tensor)

!     no off-diagonal terms in CCD
      if (.not.CCD) then
!----------------------------
!     vo: g'(a,i) += t(a,i) |
!----------------------------
        ierr=exatns_tensor_init(gamma%prime_vo,ZERO)
        ierr=exatns_tensor_contract("G(a,i)+=T(a,i)*K()",gamma%prime_vo,t1_tensor,one_dm_tensor)

!-----------------------------------------------------------
!     vo: g'(a,i) += l(m,e) * (t(ae,im) - t(e,i) * t(a,m)) |
!-----------------------------------------------------------
        ierr=exatns_tensor_init(tau2_tensor,ZERO)
        ierr=exatns_tensor_contract("T(a,e,i,m)+=R(a,e,i,m)*K()",tau2_tensor,t2_tensor,one_dm_tensor)
        ierr=exatns_tensor_contract("T(a,e,i,m)+=T(e,i)*T(a,m)",tau2_tensor,t1_tensor,t1_tensor,MINUS_ONE)
        ierr=exatns_tensor_contract("G(a,i)+=L(m,e)*T(a,e,i,m)",gamma%prime_vo,l1_tensor,tau2_tensor)

!----------------------------------------------------------------------------
!     vo: g'(a,i) -= 1/2 l(mn,ef) * (t(ef,in) * t(a,m) + t(e,i) * t(af,mn)) |
!----------------------------------------------------------------------------
        ierr=exatns_tensor_init(oo_aux_tensor,ZERO) 
        ierr=exatns_tensor_init(vv_aux_tensor,ZERO) 
        ierr=exatns_tensor_contract("X(m,i)+=L(m,n,e,f)*T(e,f,i,n)",oo_aux_tensor,l2_tensor,t2_tensor)
        ierr=exatns_tensor_contract("G(a,i)+=X(m,i)*T(a,m)",gamma%prime_vo,oo_aux_tensor,t1_tensor,MINUS_ONE_HALF)
        ierr=exatns_tensor_contract("X(a,e)+=L(m,n,e,f)*T(a,f,m,n)",vv_aux_tensor,l2_tensor,t2_tensor)
        ierr=exatns_tensor_contract("G(a,i)+=X(a,e)*T(e,i)",gamma%prime_vo,vv_aux_tensor,t1_tensor,MINUS_ONE_HALF)

!---------------------------
!     ov: g'(i,a) = l(i,a) |
!---------------------------
        ierr=exatns_tensor_init(gamma%prime_ov,ZERO)
        ierr=exatns_tensor_contract("G(i,a)+=L(i,a)*K()",gamma%prime_ov,l1_tensor,one_dm_tensor)

!--------------------------------------------------
!       ov symm: g(i,a) = 1/2 (g'(i,a) + g'(a,i)) |
!--------------------------------------------------
        ierr=exatns_tensor_init(gamma%symm_ov,ZERO)
        ierr=exatns_tensor_contract("G(i,a)+=X(i,a)*K()",gamma%symm_ov,gamma%prime_ov,one_dm_tensor,ONE_HALF)
        ierr=exatns_tensor_contract("G(i,a)+=X+(a,i)*K()",gamma%symm_ov,gamma%prime_vo,one_dm_tensor,ONE_HALF)

!--------------------------------------------------
!       vo symm: g(a,i) = symmetrized c.c. g(i,a)
!--------------------------------------------------
         ierr=exatns_tensor_init(gamma%symm_vo,ZERO)
         ierr=exatns_tensor_contract("G(a,i)+=X+(i,a)*K()",gamma%symm_vo,gamma%symm_ov,one_dm_tensor)

      end if

!------------------------------------------------
!     vv symm: g(b,a) = 1/2 (g'(b,a) + g'(a,b)) |
!------------------------------------------------
      ierr=exatns_tensor_init(gamma%symm_vv,ZERO)
      ierr=exatns_tensor_contract("G(b,a)+=X(b,a)*K()",gamma%symm_vv,gamma%prime_vv,one_dm_tensor,ONE_HALF)
      ierr=exatns_tensor_contract("G(b,a)+=X+(a,b)*K()",gamma%symm_vv,gamma%prime_vv,one_dm_tensor,ONE_HALF)

!------------------------------------------------
!     oo symm: g(j,i) = 1/2 (g'(j,i) + g'(i,j)) |
!------------------------------------------------
      ierr=exatns_tensor_init(gamma%symm_oo,ZERO)
      ierr=exatns_tensor_contract("G(j,i)+=X(j,i)*K()",gamma%symm_oo,gamma%prime_oo,one_dm_tensor,ONE_HALF)
      ierr=exatns_tensor_contract("G(j,i)+=X+(i,j)*K()",gamma%symm_oo,gamma%prime_oo,one_dm_tensor,ONE_HALF)

!     debug_exa
      if (print_level.gt.8) then
        print*, "gprime_oo  = ", print_tensornorm2(gamma%prime_oo)
        if (.not.CCD) then
          print*, "gprime_ov  = ", print_tensornorm2(gamma%prime_ov)
          print*, "gsymm_ov   = ", print_tensornorm2(gamma%symm_ov)
          print*, "gprime_vo  = ", print_tensornorm2(gamma%prime_vo)
          print*, "gsymm_vo   = ", print_tensornorm2(gamma%symm_vo)
        end if
        print*, "gprime_vv  = ", print_tensornorm2(gamma%prime_vv)
        print*, "gsymm_oo   = ", print_tensornorm2(gamma%symm_oo)
        print*, "gsymm_vv   = ", print_tensornorm2(gamma%symm_vv)
      end if

!**************
!     cleanup *
!**************

      ierr=exatns_tensor_destroy(one_dm_tensor)
      if (.not.CCD) then
        ierr=exatns_tensor_destroy(tau2_tensor)
        ierr=exatns_tensor_destroy(vv_aux_tensor)
        ierr=exatns_tensor_destroy(oo_aux_tensor)
      end if

    end subroutine calc_dm_mo

    subroutine store_dm_mo(exa_input,gamma)

      use exacorr_global,         only : get_nsolutions
      use talsh
      use tensor_algebra
      use talsh_common_routines,  only : print_tensor
#if defined (VAR_MPI)
      use interface_to_mpi

!     variables for mpi
      integer(C_INT)    :: read_node
      integer           :: my_MPI_rank
#endif


!     input variables
      type(exacc_input),    intent(in   ) :: exa_input
      type(exatns_dm_tens), intent(inout) :: gamma
!     talsh tensor to get local copy
      type(talsh_tens_t)           :: tensor_slice
      integer(INTD), dimension(2)  :: h_dims
!     matrix dimensions
      integer :: nocc, nvir, nkr_occ, nkr_vir, nesh
!     body pointers
      type(C_PTR)                     :: body_p2
      complex(8), pointer, contiguous :: dm_block(:,:)
      complex(8), allocatable         :: dm_mo_fort(:,:)
      real(8), allocatable            :: dm_diracformat(:,:,:)
!     complete mo list
      integer, allocatable :: mo_list(:)
!     CCD switch
      logical :: CCD
!     error code
      integer :: ierr
!     loop variables
      integer :: i, j, norbt
      integer :: i_unbar, j_unbar, j_bar
!     variables for talsh, density matrix needs to be smaller than buffer size 
      integer(C_SIZE_T)               :: buf_size=1024_8*1024_8*1024_8 
      integer(C_INT)                  :: host_arg_max

      CCD = exa_input%CCD

      call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
      !call interface_mpi_comm_size(global_communicator,read_node)
      read_node=0

!***************
!     set dims *
!***************

      nocc    = exa_input%nocc
      nvir    = exa_input%nvir
      nkr_occ = exa_input%nkr_occ
      nkr_vir = exa_input%nkr_vir

!***************************************
!     copy DM blocks into one large DM *
!***************************************

      allocate (dm_mo_fort(nocc+nvir,nocc+nvir))
      dm_mo_fort = ZERO

      buf_size=exa_input%talsh_buff*buf_size
      ierr=talsh_init(buf_size,host_arg_max)
      call print_date('Initialized talsh library')
      write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

!     -----------
!     | OO | OV |
!     -----------
!     | VO | VV |
!     -----------
!     oo
      h_dims=nocc
      ierr=talsh_tensor_construct(tensor_slice,C8,h_dims,init_val=ZERO)
      ierr=exatns_tensor_get_slice(gamma%symm_oo,tensor_slice)
      ierr=talsh_tensor_get_body_access(tensor_slice,body_p2,C8,read_node,DEV_HOST)
      call c_f_pointer(body_p2,dm_block,h_dims) ! to use <dm_mo_block> as a regular Fortran 2d array
      dm_mo_fort(1:nocc,1:nocc) = dm_block
      ierr=talsh_tensor_destruct(tensor_slice)

!     vo + ov
      if (.not.CCD) then
        h_dims(1)=nvir
        h_dims(2)=nocc
        ierr=talsh_tensor_construct(tensor_slice,C8,h_dims,init_val=ZERO)
        ierr=exatns_tensor_get_slice(gamma%symm_vo,tensor_slice)
        ierr=talsh_tensor_get_body_access(tensor_slice,body_p2,C8,read_node,DEV_HOST)
        call c_f_pointer(body_p2,dm_block,h_dims) ! to use <dm_mo_block> as a regular Fortran 2d array
        dm_mo_fort(1+nocc:nocc+nvir,1:nocc) = dm_block
        ierr=talsh_tensor_destruct(tensor_slice)

        h_dims(1)=nocc
        h_dims(2)=nvir
        ierr=talsh_tensor_construct(tensor_slice,C8,h_dims,init_val=ZERO)
        ierr=exatns_tensor_get_slice(gamma%symm_ov,tensor_slice)
        ierr=talsh_tensor_get_body_access(tensor_slice,body_p2,C8,read_node,DEV_HOST)
        call c_f_pointer(body_p2,dm_block,h_dims) ! to use <dm_mo_block> as a regular Fortran 2d array
        dm_mo_fort(1:nocc,1+nocc:nocc+nvir) = dm_block
        ierr=talsh_tensor_destruct(tensor_slice)
      end if

!     vv
      h_dims=nvir
      ierr=talsh_tensor_construct(tensor_slice,C8,h_dims,init_val=ZERO)
      ierr=exatns_tensor_get_slice(gamma%symm_vv,tensor_slice)
      ierr=talsh_tensor_get_body_access(tensor_slice,body_p2,C8,read_node,DEV_HOST)
      call c_f_pointer(body_p2,dm_block,h_dims) ! to use <dm_mo_block> as a regular Fortran 2d array
      dm_mo_fort(1+nocc:nocc+nvir,1+nocc:nocc+nvir) = dm_block
      ierr=talsh_tensor_destruct(tensor_slice)
      
      ierr = talsh_shutdown()
!*****************************************************************************************************************************
!     Determine where this should go in the full list of MOs (NB: DIRAC assumes restricted densities and uses Kramers pairs) *
!*****************************************************************************************************************************
      
      norbt = nkr_occ + nkr_vir
      allocate(mo_list(norbt))
      mo_list(1:nkr_occ)                 = exa_input%mokr_occ
      mo_list(1+nkr_occ:nkr_occ+nkr_vir) = exa_input%mokr_vir

      nesh = get_nsolutions() / 2
      allocate (dm_diracformat(nesh,nesh,4))
      dm_diracformat = ZERO
      do j = 1, norbt
        j_unbar = 2*j-1
        j_bar   = 2*j
        do i = 1, norbt
          ! We copy the unbar-unbar part of the density matrix, as well as the unbar-bar part
          ! This follows from eq. 27 of Shee, Visscher, Saue JCP 145 184107, 2016
          ! it would be good to test whether the Kramers-restricted assumption holds
          ! Or even to average over the two density matrices and enforce a Kramers-restricted density
          i_unbar = 2*i-1
          dm_diracformat(mo_list(i),mo_list(j),1) = real(dm_mo_fort(i_unbar,j_unbar))
          dm_diracformat(mo_list(i),mo_list(j),2) = imag(dm_mo_fort(i_unbar,j_unbar))
          dm_diracformat(mo_list(i),mo_list(j),3) = real(dm_mo_fort(i_unbar,j_bar))
          dm_diracformat(mo_list(i),mo_list(j),4) = imag(dm_mo_fort(i_unbar,j_bar))
        end do
      end do

      call store_cc_denstiy_exacorr('CCSD',dm_diracformat,nesh,4)

      deallocate (dm_diracformat)
      deallocate (mo_list)
      deallocate (dm_mo_fort)

    end subroutine store_dm_mo


    subroutine store_cc_denstiy_exacorr(wf,dmo,norbt,nz)

    use exacorr_utils,  only : get_free_fileunit

    implicit none 

    character(4)  :: wf
    real(8)       :: dmo(norbt*norbt*nz)
    integer       :: norbt
    integer       :: nz
    integer       :: ierr
    integer       :: mccdens

    character(8)  :: label, eoflabel

    WRITE (*,*) "Storing the ",wf," density matrix on file" 
    label = wf//'DENS'
    eoflabel = 'EOFLABEL'

    call get_free_fileunit(mccdens)

    OPEN(mccdens,FILE='CCDENS',STATUS='NEW',FORM='UNFORMATTED', &
      ACCESS='SEQUENTIAL')

    ierr=0
    REWIND (mccdens)
    CALL NEWLAB(label,mccdens,ierr)
    if (ierr.ne.0) write(*,*) "error in store_cc_denstiy_exacorr: ", ierr
    write (mccdens) dmo
    CALL NEWLAB(eoflabel,mccdens,ierr)
    if (ierr.ne.0) write(*,*) "error in store_cc_denstiy_exacorr: ", ierr
    CLOSE(Unit=mccdens,Status='KEEP')

    end subroutine store_cc_denstiy_exacorr

#endif

end module exacorr_gradient
