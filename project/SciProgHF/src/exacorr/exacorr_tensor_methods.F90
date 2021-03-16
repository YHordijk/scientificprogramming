module exacorr_tensor_methods

!This module contains routines to modify (initialize, scale with denominators, etc.) that can be registered with exatensor
#ifdef VAR_MPI

        use exatensor
        use, intrinsic:: ISO_C_BINDING

        implicit none
        private

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

        !Initializer for the 2e AO integral tensor:
        type, extends(tens_method_uni_t), public:: compute_2e_ao_tensor_t
         contains
          procedure, public:: apply=>compute_2e_ao_tensor
        end type compute_2e_ao_tensor_t
        
        !Initializer for the 2e AO integral tensor, allows smaller segments:
        type, extends(tens_method_uni_t), public:: compute_2e_ao_seg_t
         contains
          procedure, public:: apply=>compute_2e_ao_seg
        end type compute_2e_ao_seg_t

        !Initializer for the mo coefficients tensor:
        type, extends(tens_method_uni_t), public:: init_mocoef_tensor_t
         complex(8), allocatable, private:: mocoef_all(:,:) ! stores mo-coefficients for both indices, both spins
         contains
          procedure, public:: init_mocoef_ctor=>init_mocoef
          procedure, public:: apply=>init_mocoef_tensor
        end type init_mocoef_tensor_t

        !Initializer for the density matrix tensor:
        type, extends(tens_method_uni_t), public:: init_dm_tensor_t
         complex(8), allocatable, private:: mocoef_all_1(:,:,:),mocoef_all_2(:,:,:)
         contains
          procedure, public:: init_dm_ctor=>init_dm
          procedure, public:: apply=>init_dm_tensor
        end type init_dm_tensor_t

        !Initializer for the half-transformed integrals:
        type, extends(tens_method_uni_t), public:: init_ht_tensor_t
         complex(8), allocatable, private:: mocoef_all_1(:,:,:),mocoef_all_2(:,:,:)
         contains
          procedure, public:: init_ht_ctor=>init_ht
          procedure, public:: apply=>init_ht_tensor
        end type init_ht_tensor_t

        !setting elements of tensor to zero (1-real, 2-imaginary, 3-complex)
        type, extends(tens_method_uni_t), public:: set_zero_t
         integer,private  :: mode 
         contains
          procedure, public:: set_zero_init
          procedure, public:: reset=>set_zero_reset
          procedure, public:: pack=>set_zero_pack
          procedure, public:: unpack=>set_zero_unpack
          procedure, public:: apply=>set_zero_apply
        end type set_zero_t

        ! fill delta tensor
        type, extends(tens_method_uni_t), public:: delta_t
         integer(INTL), private :: ind(4)
         contains
          procedure, public:: delta_t_init
          procedure, public:: reset=>delta_t_reset
          procedure, public:: pack=>delta_t_pack
          procedure, public:: unpack=>delta_t_unpack
          procedure, public:: apply=>delta_t_apply
        end type delta_t

        !square tensor
        type, extends(tens_method_uni_t), public:: square_t
         contains
          procedure, public:: apply=>square_t_apply
        end type square_t

        !set diagonal to constant
        type, extends(tens_method_uni_t), public:: set_diagonal_t
         complex(8), private :: constant
         logical, private :: offdiag
         contains
          procedure, public:: set_diagonal_init
          procedure, public:: apply=>set_diagonal_apply
        end type set_diagonal_t

        !set orbital energies on diagonal
        type, extends(tens_method_uni_t), public:: set_energy_t
         real(8), allocatable, private :: eps_occ(:),eps_vir(:)
         integer, private              :: occ_id, vir_id   
         real(8), private              :: level_shift
         integer, private              :: norb(2) 
         contains
          procedure, public:: set_energy_init
          procedure, public:: reset=>set_energy_reset
          procedure, public:: pack=>set_energy_pack
          procedure, public:: unpack=>set_energy_unpack
          procedure, public:: apply=>set_energy_apply
        end type set_energy_t

        !Scale by denominators:
        type, extends(tens_method_uni_t), public:: denom_t
          real(8), allocatable, private :: eps_occ(:),eps_vir(:)
          real(8), private              :: level_shift
          integer(INTD), private        :: occ_id
          integer, private              :: norb(2)
          contains
            procedure, public::  denom_t_init
            procedure, public::  reset=>denom_t_reset
            procedure, public::  pack=>denom_t_pack
            procedure, public::  unpack=>denom_t_unpack
            procedure, public::  apply=>denom_t_apply
        end type denom_t

        !Scale by denominators for delta_triples
        type, extends(tens_method_uni_t), public:: denom3_t
          real(8), allocatable, private :: eps_occ(:),eps_vir(:)
          real(8), private              :: level_shift
          integer, private              :: norb(2)
          real(8), private              :: eps_ijk
          contains
            procedure, public::  denom3_t_init
            procedure, public::  reset=>denom3_t_reset
            procedure, public::  pack=>denom3_t_pack
            procedure, public::  unpack=>denom3_t_unpack
            procedure, public::  apply=>denom3_t_apply
        end type denom3_t

        !Initializer for the one electron integrals
        type, extends(tens_method_uni_t), public:: set_1el_t
         complex(8), allocatable, private :: integrals(:,:)
         contains
          procedure, public:: set_1el_init
          procedure, public:: apply=>set_1el_apply
        end type set_1el_t

        !get property integrals
        type, extends(tens_method_uni_t), public:: set_ff_t
          integer, private              :: occ_id, vir_id 
          integer, private              :: nocc, nvir, min_occ
          integer, allocatable, private :: mo_occ(:), mo_vir(:) 
          integer, private              :: i_prop
         contains
          procedure, public:: set_ff_init
          procedure, public:: reset=>set_ff_reset
          procedure, public:: pack=>set_ff_pack
          procedure, public:: unpack=>set_ff_unpack
          procedure, public:: apply=>set_ff_apply
        end type set_ff_t

        !set projecton matrix
        type, extends(tens_method_uni_t), public:: set_projection_t
         complex(8), private :: constant
         contains
          procedure, public:: set_projection_init
          procedure, public:: apply=>set_projection_apply
        end type set_projection_t

        !Cholesky: Diagonal tensor
        type, extends(tens_method_uni_t), public:: chol_diag
         contains
          procedure, public:: apply=>chol_diag_apply
        end type chol_diag

        !Cholesky: Offdiagonal tensor
        type, extends(tens_method_uni_t), public:: chol_J
         integer(INTL), private :: ind(1:2)
         contains
          procedure, public:: chol_J_init
          procedure, public:: reset=>chol_J_reset
          procedure, public:: pack=>chol_J_pack
          procedure, public:: unpack=>chol_J_unpack
          procedure, public:: apply=>chol_J_apply
        end type chol_J

       contains

        ! --------------------------------------------------------------------------------

        function compute_2e_ao_tensor(this,tensor,scalar) result(ierr)
         implicit none
         integer(INTD):: ierr,n
         class(compute_2e_ao_tensor_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         n=tens%num_dims
         if (n.ne.4) call quit('Error: dimension not 4 in compute_2e_ao_tensor')
         ierr=tens_init_ao_integrals(tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n))
         return
        end function compute_2e_ao_tensor

        integer function tens_init_ao_integrals(tens_body, data_kind, tens_dims, tens_bases) result (ierr)
        !Initialization routine that is used to fill a tensor with AO integrals

         use exacorr_datatypes
         use exacorr_global
         use module_interest_eri

         implicit none

         type(C_PTR), value:: tens_body           !in: C pointer to the tensor body (tensor elements)
         integer(INTD), value:: data_kind            !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(INTL), intent(in):: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(INTL), intent(in):: tens_bases(1:*) !in: base offsets for each tensor dimension)

         integer:: nao(4)
         integer:: basis_angular ! 1=cartesian, 2=spherical
         complex(8),pointer :: ao_tens(:,:,:,:)
         integer :: first_ao(1:4)

         integer :: nshells(4)
         type(basis_func_info_t), allocatable :: gto1(:), gto2(:), gto3(:), gto4(:)  ! arrays with basis function information

         integer :: i, j, k, l, ish, jsh, ksh, lsh, ioff, joff, koff, loff, ijkl
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         !$OMP threadprivate(gout)
         real(8), parameter :: FIJKL = 1.0

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

!        Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
         call interest_initialize(.false.)

         ! map the arguments to names used in the context of integral evaluation
         first_ao(1:4) = tens_bases(1:4)
         nao(1:4) = tens_dims(1:4)
         call c_f_pointer(tens_body,ao_tens,nao(1:4)) !map C pointer to 4-dimensional array

!        Allocate and get shell information for all shells that we will treat
         call get_gtos(first_ao(1),nao(1),gto1,nshells(1))
         call get_gtos(first_ao(2),nao(2),gto2,nshells(2))
         call get_gtos(first_ao(3),nao(3),gto3,nshells(3))
         call get_gtos(first_ao(4),nao(4),gto4,nshells(4))
         basis_angular = get_basis_angular()

!        Intranode parallelization should be possible in the loops below because Interest is thread safe

         loff = 0
         do lsh = 1, nshells(4)
            ll   =  gto4(lsh)%orb_momentum
            el   =  gto4(lsh)%exponent(1)
            xl   =  gto4(lsh)%coord(1)
            yl   =  gto4(lsh)%coord(2)
            zl   =  gto4(lsh)%coord(3)
            cl   =  gto4(lsh)%coefficient(1)
            nl   =  nfunctions(ll,basis_angular)
            koff = 0
            do ksh = 1, nshells(3)
               lk   =  gto3(ksh)%orb_momentum
               ek   =  gto3(ksh)%exponent(1)
               xk   =  gto3(ksh)%coord(1)
               yk   =  gto3(ksh)%coord(2)
               zk   =  gto3(ksh)%coord(3)
               ck   =  gto3(ksh)%coefficient(1)
               nk   =  nfunctions(lk,basis_angular)
               joff = 0
               do jsh = 1, nshells(2)
                  lj   =  gto2(jsh)%orb_momentum
                  ej   =  gto2(jsh)%exponent(1)
                  xj   =  gto2(jsh)%coord(1)
                  yj   =  gto2(jsh)%coord(2)
                  zj   =  gto2(jsh)%coord(3)
                  cj   =  gto2(jsh)%coefficient(1)
                  nj   =  nfunctions(lj,basis_angular)
                  ioff = 0
                  do ish = 1, nshells(1)
                     li   =  gto1(ish)%orb_momentum
                     ei   =  gto1(ish)%exponent(1)
                     xi   =  gto1(ish)%coord(1)
                     yi   =  gto1(ish)%coord(2)
                     zi   =  gto1(ish)%coord(3)
                     ci   =  gto1(ish)%coefficient(1)
                     ni   =  nfunctions(li,basis_angular)
!                    output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
                     call interest_eri('llll',FIJKL,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )
                     do l = 1, nl
                        do k = 1, nk
                           do j = 1, nj
                              do i = 1, ni
                                 ijkl = (l-1)*nk*nj*ni+(k-1)*nj*ni+(j-1)*ni+i
                                 ao_tens(i+ioff,j+joff,k+koff,l+loff) = dcmplx(gout(ijkl),0.D0)
                              end do
                           end do
                        end do
                     end do
                     ioff = ioff + ni
                  end do
                  joff = joff + nj
               end do
               koff = koff + nk
            end do
            loff = loff + nl
         end do

         deallocate(gto1)
         deallocate(gto2)
         deallocate(gto3)
         deallocate(gto4)

        end function tens_init_ao_integrals

        ! --------------------------------------------------------------------------------

        function compute_2e_ao_seg(this,tensor,scalar) result(ierr)
          implicit none
          integer(INTD):: ierr,n
          class(compute_2e_ao_seg_t), intent(in):: this
          class(tens_rcrsv_t), intent(inout):: tensor
          complex(8), intent(inout), optional:: scalar
          type(tens_dense_t):: tens

          ierr=0
          tens=tensor%get_dense_adapter(ierr)
          if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
          n=tens%num_dims
          if (n.ne.4) call quit('Error: dimension not 4 in compute_2e_ao_tensor')
          ierr=tens_init_ao_seg(tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n))
          return

        end function compute_2e_ao_seg

        integer function tens_init_ao_seg(tens_body, data_kind, tens_dims, tens_bases) result (ierr)
          !Initialization routine that is used to fill a tensor with AO integrals

          use exacorr_datatypes
          use exacorr_global
          use module_interest_eri

          implicit none

          type(C_PTR), value:: tens_body           !in: C pointer to the tensorbody (tensor elements)
          integer(INTD), value:: data_kind            !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
          integer(INTL), intent(in):: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
          integer(INTL), intent(in):: tens_bases(1:*) !in: base offsets for each tensor dimension)

          integer:: nao(4)
          integer:: basis_angular ! 1=cartesian, 2=spherical
          complex(8),pointer :: ao_tens(:,:,:,:)
          integer :: first_ao(1:4), off_ao(4)
          integer :: nshells(4)
          type(basis_func_info_t), allocatable :: gto1(:), gto2(:), gto3(:), gto4(:)  ! arrays with basis function information
          integer :: i, j, k, l, ish, jsh, ksh, lsh, ioff, joff, koff, loff, ijkl
          integer :: li,lj,lk,ll,ni,nj,nk,nl,i0,j0,k0,l0
          real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
          real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
          !$OMP threadprivate(gout)
          real(8), parameter :: FIJKL = 1.0


          ierr=0 !set error code to success
          if(data_kind /= EXA_DATA_KIND_C8) then
            ierr=-1 !invalid data kind requested
            return
          endif

!         Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
          call interest_initialize(.false.)

          ! map the arguments to names used in the context of integral evaluation
          first_ao(1:4) = tens_bases(1:4)
          nao(1:4) = tens_dims(1:4)
          call c_f_pointer(tens_body,ao_tens,nao(1:4)) !map C pointer to 4-dimensional array

!         Allocate and get shell information for all shells that we will treat
          call get_gtos(first_ao(1),nao(1),gto1,nshells(1))
          call get_gtos(first_ao(2),nao(2),gto2,nshells(2))
          call get_gtos(first_ao(3),nao(3),gto3,nshells(3))
          call get_gtos(first_ao(4),nao(4),gto4,nshells(4))
          basis_angular = get_basis_angular()

!         get information about shell offset
          do i=1,4
            call get_shell_offset(first_ao(i),j)
            off_ao(i)=tens_bases(i)-j-1
          end do

!         Intranode parallelization should be possible in the loops below because Interest is thread safe

          loff = 0
          do lsh = 1, nshells(4)
            ll   =  gto4(lsh)%orb_momentum
            el   =  gto4(lsh)%exponent(1)
            xl   =  gto4(lsh)%coord(1)
            yl   =  gto4(lsh)%coord(2)
            zl   =  gto4(lsh)%coord(3)
            cl   =  gto4(lsh)%coefficient(1)
            nl   =  nfunctions(ll,basis_angular)
            if(lsh.eq.1) then
              l0=off_ao(4)
            else
              l0=0
            end if
            koff = 0
            do ksh = 1, nshells(3)
              lk   =  gto3(ksh)%orb_momentum
              ek   =  gto3(ksh)%exponent(1)
              xk   =  gto3(ksh)%coord(1)
              yk   =  gto3(ksh)%coord(2)
              zk   =  gto3(ksh)%coord(3)
              ck   =  gto3(ksh)%coefficient(1)
              nk   =  nfunctions(lk,basis_angular)
              if(ksh.eq.1) then
                k0=off_ao(3)
              else
                k0=0
              end if
              joff = 0
              do jsh = 1, nshells(2)
                lj   =  gto2(jsh)%orb_momentum
                ej   =  gto2(jsh)%exponent(1)
                xj   =  gto2(jsh)%coord(1)
                yj   =  gto2(jsh)%coord(2)
                zj   =  gto2(jsh)%coord(3)
                cj   =  gto2(jsh)%coefficient(1)
                nj   =  nfunctions(lj,basis_angular)
                if(jsh.eq.1) then
                  j0=off_ao(2)
                else
                  j0=0
                end if
                ioff = 0
                do ish = 1, nshells(1)
                  li   =  gto1(ish)%orb_momentum
                  ei   =  gto1(ish)%exponent(1)
                  xi   =  gto1(ish)%coord(1)
                  yi   =  gto1(ish)%coord(2)
                  zi   =  gto1(ish)%coord(3)
                  ci   =  gto1(ish)%coefficient(1)
                  ni   =  nfunctions(li,basis_angular)
                  if(ish.eq.1) then
                    i0=off_ao(1)
                  else
                    i0=0
                  end if

!                 output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
                  call interest_eri('llll',FIJKL,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )
                  do l = 1, nl
                    if (l+loff.gt.tens_dims(4)) exit
                    do k = 1, nk
                       if (k+koff.gt.tens_dims(3)) exit
                       do j = 1, nj
                         if (j+joff.gt.tens_dims(2)) exit
                         do i = 1, ni
                           if (i+ioff.gt.tens_dims(1)) exit
                           ijkl = (l+l0-1)*nk*nj*ni+(k+k0-1)*nj*ni+(j+j0-1)*ni+(i+i0)
                           ao_tens(i+ioff,j+joff,k+koff,l+loff) = dcmplx(gout(ijkl),0.D0)
                         end do
                       end do
                     end do
                   end do
                   ioff = ioff + ni
                end do
                joff = joff + nj
              end do
              koff = koff + nk
            end do
            loff = loff + nl
          end do

          deallocate(gto1)
          deallocate(gto2)
          deallocate(gto3)
          deallocate(gto4)

        end function tens_init_ao_seg

        ! --------------------------------------------------------------------------------

        subroutine init_mocoef(this,mo_list,nmo,spin)

         use exacorr_mo
         use exacorr_global

         implicit none
         integer, intent(in)       :: nmo          ! the length of the mo basis
         integer, intent(in)       :: mo_list(:)   ! and their indices
         integer(INTD), intent(in) :: spin
         class(init_mocoef_tensor_t), intent(out):: this

         integer  :: nao
         type(cmo) :: cspinor

!        Copy the coefficients into the alpha or beta array, using the quaternion relations as defined in
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
!        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
!        using Kramers' symmetry.

         nao = get_nao()
         allocate (this%mocoef_all(nao,nmo))
         call get_mo_coefficients (cspinor,mo_list,nmo)
         this%mocoef_all(:,:) = cspinor%coeff(:,:,spin)
         call dealloc_mo(cspinor)

        end subroutine init_mocoef

        function init_mocoef_tensor(this,tensor,scalar) result(ierr)

         implicit none
         integer(INTD):: ierr,n
         class(init_mocoef_tensor_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         n=tens%num_dims
         if (n.ne.2) call quit('Error: dimension not 2 in init_mocoef_tensor')
         ierr=tens_init_mocoef(tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n), &
                               this%mocoef_all)

        end function init_mocoef_tensor

        integer function tens_init_mocoef(tens_body,data_kind,tens_dims,tens_bases, &
                         mocoef_all) result (ierr)
        !Initialization routine that is used to fill a tensor with mo coefficients

         implicit none

         type(C_PTR), value        :: tens_body           !in: C pointer to the tensor body (tensor elements)
         integer(INTD), value      :: data_kind            !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)
         complex(8), intent(in)    :: mocoef_all(:,:)

         integer :: nao, nmo
         integer :: first_ao, last_ao, first_mo, last_mo
         complex(8),pointer :: mocoef_tens(:,:)

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

         nao = tens_dims(1)
         nmo = tens_dims(2)
         call c_f_pointer(tens_body,mocoef_tens,(/nao,nmo/)) !map C pointer to 2-dimensional array

         first_ao = tens_bases(1)
         first_mo = tens_bases(2)
         last_ao  = first_ao + nao - 1
         last_mo  = first_mo + nmo - 1

         mocoef_tens(1:nao,1:nmo) = mocoef_all(first_ao:last_ao,first_mo:last_mo)

        end function tens_init_mocoef

        ! --------------------------------------------------------------------------------

        subroutine init_dm(this,mo_list,nmo)

         use exacorr_mo
         use exacorr_global

         implicit none
         integer, intent(in) :: nmo(2)       ! the length of the mo basis
         integer, intent(in) :: mo_list(:)   ! and their indices
         class(init_dm_tensor_t), intent(out):: this

         integer  :: nao
         type(cmo) :: cspinor

!        Copy the coefficients into the alpha or beta array, using the quaternion relations as defined in
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
!        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
!        using Kramers' symmetry.

         nao = get_nao()
         allocate(this%mocoef_all_1(nao,nmo(1),2))
         allocate(this%mocoef_all_2(nao,nmo(2),2))
         call get_mo_coefficients(cspinor,mo_list(1:nmo(1)+nmo(2)),nmo(1)+nmo(2))
         this%mocoef_all_1 = cspinor%coeff(:,1:nmo(1),:)
         this%mocoef_all_2 = cspinor%coeff(:,nmo(1)+1:nmo(1)+nmo(2),:)
         call dealloc_mo(cspinor)

        end subroutine init_dm

        function init_dm_tensor(this,tensor,scalar) result(ierr)

         implicit none
         integer(INTD):: ierr,n
         class(init_dm_tensor_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         n=tens%num_dims
         if (n.ne.4) call quit('Error: dimension not 4 in init_dm_tensor')
         ierr=tens_init_dm(tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n), &
                               this%mocoef_all_1,this%mocoef_all_2)

        end function init_dm_tensor

        integer function tens_init_dm(tens_body,data_kind,tens_dims,tens_bases, &
                         mocoef_all_1,mocoef_all_2) result (ierr)
        !Initialization routine that is used to fill a tensor with mo coefficients

         implicit none

         type(C_PTR), value        :: tens_body           !in: C pointer to the tensor body (tensor elements)
         integer(INTD), value      :: data_kind            !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)
         complex(8), intent(in)    :: mocoef_all_1(:,:,:),mocoef_all_2(:,:,:)

         integer :: p,q,i,j,ia,ja,pa,qa,nao(2),nmo(2)
         complex(8),pointer :: dm_tens(:,:,:,:)

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

         nao = tens_dims(1:2)
         nmo = tens_dims(3:4)
         call c_f_pointer(tens_body,dm_tens,(/nao(1),nao(2),nmo(1),nmo(2)/)) !map C pointer to 4-dimensional array

         ! carry out the spin integration and get the tensor
         do j = 1, nmo(2)
            ja = tens_bases(4) + j - 1 
            do i = 1, nmo(1)
               ia = tens_bases(3) + i - 1 
               do q = 1, nao(2)
                  qa = tens_bases(2) + q - 1 
                  do p = 1, nao(1)
                     pa = tens_bases(1) + p - 1 
                     dm_tens(p,q,i,j) = &
                       dconjg(mocoef_all_1(pa,ia,1))*mocoef_all_2(qa,ja,1) & ! alpha-alpha
                     + dconjg(mocoef_all_1(pa,ia,2))*mocoef_all_2(qa,ja,2)   ! beta-beta
                  end do
               end do
            end do
         end do

        end function tens_init_dm

        ! --------------------------------------------------------------------------------
        
        subroutine init_ht(this,mo_list,nmo)

         use exacorr_mo
         use exacorr_global

         implicit none
         integer, intent(in) :: nmo(2)       ! the length of the mo basis
         integer, intent(in) :: mo_list(:)   ! and their indices
         class(init_ht_tensor_t), intent(out):: this

         integer  :: nao
         type(cmo) :: cspinor

!        Copy the coefficients into the alpha or beta array, using the quaternion relations as defined in
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
!        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
!        using Kramers' symmetry.

         nao = get_nao()
         allocate(this%mocoef_all_1(nao,nmo(1),2))
         allocate(this%mocoef_all_2(nao,nmo(2),2))
         call get_mo_coefficients(cspinor,mo_list(1:nmo(1)+nmo(2)),nmo(1)+nmo(2))
         this%mocoef_all_1 = cspinor%coeff(:,1:nmo(1),:)
         this%mocoef_all_2 = cspinor%coeff(:,nmo(1)+1:nmo(1)+nmo(2),:)
         call dealloc_mo(cspinor)

        end subroutine init_ht

        function init_ht_tensor(this,tensor,scalar) result(ierr)

         implicit none
         integer(INTD):: ierr,n
         class(init_ht_tensor_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         n=tens%num_dims
         if (n.ne.4) call quit('Error: dimension not 4 in init_dm_tensor')
         ierr=tens_init_ht(tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n), &
                               this%mocoef_all_1,this%mocoef_all_2)

        end function init_ht_tensor

        integer function tens_init_ht(tens_body,data_kind,tens_dims,tens_bases, &
                         mocoef_all_1,mocoef_all_2) result (ierr)
        !Initialization routine that is used to fill a tensor with mo coefficients

         use exacorr_datatypes
         use exacorr_global
         use module_interest_eri

         implicit none

         type(C_PTR), value        :: tens_body           !in: C pointer to the tensor body (tensor elements)
         integer(INTD), value      :: data_kind            !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)
         complex(8), intent(in)    :: mocoef_all_1(:,:,:),mocoef_all_2(:,:,:)

         complex(8),pointer :: ht_tens(:,:,:,:)

         integer :: first_ao(1:4)
         integer :: nshells(4)
         integer :: basis_angular 
         type(basis_func_info_t), allocatable :: gto1(:), gto2(:), gto3(:), gto4(:)  ! arrays with basis function information

         integer :: ioff,joff,koff,loff,kmo,lmo,ka,la,pa,qa,ra,sa,nao(4),nmo(2),ijkl
         integer :: i, j, k, l, ish, jsh, ksh, lsh
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         !$OMP threadprivate(gout)
         complex(8) :: dP
         real(8), parameter :: FIJKL = 1.0

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

!        Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
         call interest_initialize(.false.)

         ! map the arguments to names used in the context of integral evaluation
         first_ao(1:2) = tens_bases(1:2)
         first_ao(3:4) = 1
         nao(1:2)      = tens_dims(1:2)
         nmo(1:2)      = tens_dims(3:4)
         nao(3:4)      = get_nao()
         call c_f_pointer(tens_body,ht_tens,(/nao(1),nao(2),nmo(1),nmo(2)/)) !map C pointer to 4-dimensional array
         ht_tens = ZERO

!        Allocate and get shell information for all shells that we will treat
         call get_gtos(first_ao(1),nao(1),gto1,nshells(1))
         call get_gtos(first_ao(2),nao(2),gto2,nshells(2))
         call get_gtos(first_ao(3),nao(3),gto3,nshells(3))
         call get_gtos(first_ao(4),nao(4),gto4,nshells(4))
         basis_angular = get_basis_angular()

         loff = 0
         do lsh = 1, nshells(4) ! loop over all ao's, this index is transformed to mo
            ll   =  gto4(lsh)%orb_momentum
            el   =  gto4(lsh)%exponent(1)
            xl   =  gto4(lsh)%coord(1)
            yl   =  gto4(lsh)%coord(2)
            zl   =  gto4(lsh)%coord(3)
            cl   =  gto4(lsh)%coefficient(1)
            nl   =  nfunctions(ll,basis_angular)
            koff = 0
            do ksh = 1, nshells(3) ! loop over all ao's, this index is transformed to mo
               lk   =  gto3(ksh)%orb_momentum
               ek   =  gto3(ksh)%exponent(1)
               xk   =  gto3(ksh)%coord(1)
               yk   =  gto3(ksh)%coord(2)
               zk   =  gto3(ksh)%coord(3)
               ck   =  gto3(ksh)%coefficient(1)
               nk   =  nfunctions(lk,basis_angular)
               joff = 0
               do jsh = 1, nshells(2) ! loop over the ao's active in the HT integral tensor block
                  lj   =  gto2(jsh)%orb_momentum
                  ej   =  gto2(jsh)%exponent(1)
                  xj   =  gto2(jsh)%coord(1)
                  yj   =  gto2(jsh)%coord(2)
                  zj   =  gto2(jsh)%coord(3)
                  cj   =  gto2(jsh)%coefficient(1)
                  nj   =  nfunctions(lj,basis_angular)
                  ioff = 0
                  do ish = 1, nshells(1) ! loop over the ao's active in the HT integral tensor block
                     li   =  gto1(ish)%orb_momentum
                     ei   =  gto1(ish)%exponent(1)
                     xi   =  gto1(ish)%coord(1)
                     yi   =  gto1(ish)%coord(2)
                     zi   =  gto1(ish)%coord(3)
                     ci   =  gto1(ish)%coefficient(1)
                     ni   =  nfunctions(li,basis_angular)
                     ! output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
                     call interest_eri('llll',FIJKL,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )

                     do lmo = 1, nmo(2) ! loop over the mo's active in the HT integral tensor block
                        la = tens_bases(4) + lmo - 1 
                        do kmo = 1, nmo(1) ! loop over the mo's active in the HT integral tensor block
                           ka = tens_bases(3) + kmo - 1 
                           ijkl = 0
                           do l = 1, nl
                              sa = loff + l
                              do k = 1, nk
                                 ra = koff + k
                                 dP = dconjg(mocoef_all_1(ra,ka,1))*mocoef_all_2(sa,la,1) & ! alpha-alpha
                                    + dconjg(mocoef_all_1(ra,ka,2))*mocoef_all_2(sa,la,2)   ! beta-beta
                                 do j = 1, nj
                                    qa = joff + j
                                    do i = 1, ni
                                       pa = ioff + i
                                       ijkl = ijkl + 1
                                       ht_tens(pa,qa,kmo,lmo) = ht_tens(pa,qa,kmo,lmo) &
                                                              + dP * gout(ijkl)
                                    end do ! i
                                 end do ! j
                              end do ! k
                           end do ! l
                        end do ! kmo
                     end do ! lmo
                     ioff = ioff + ni
                  end do
                  joff = joff + nj
               end do
               koff = koff + nk
            end do
            loff = loff + nl
         end do

         deallocate(gto1)
         deallocate(gto2)
         deallocate(gto3)
         deallocate(gto4)

        end function tens_init_ht

        ! --------------------------------------------------------------------------------

      subroutine set_zero_init(this,mode)
        implicit none
        class(set_zero_t), intent(out):: this
        integer, intent(in)           :: mode
        
        if (mode.lt.1 .or. mode.gt.3) call quit ('Error in set_zero_init: wrong mode')
        this%mode=mode

      end subroutine set_zero_init

      subroutine set_zero_reset(this,mode)
        implicit none
        class(set_zero_t), intent(inout) :: this
        integer, intent(in)           :: mode
        
        if (mode.lt.1 .or. mode.gt.3) call quit ('Error in set_zero_reset: wrong mode')
        this%mode=mode

      end subroutine set_zero_reset

      subroutine set_zero_pack(this,packet,ierr)
        implicit none
        class(set_zero_t), intent(in)        :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr

        call pack_builtin(packet,this%mode)

        if(present(ierr)) ierr=0  
        return
      end subroutine set_zero_pack

      subroutine set_zero_unpack(this,packet,ierr)
        implicit none
        class(set_zero_t), intent(inout)     :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr

        call unpack_builtin(packet,this%mode)

        if(present(ierr)) ierr=0
        return
      end subroutine set_zero_unpack

      function set_zero_apply(this,tensor,scalar) result(ierr)
        implicit none
        integer(INTD) :: ierr, rank
        class(set_zero_t), intent(in) :: this
        class(tens_rcrsv_t), intent(inout) :: tensor
        complex(8), intent(inout), optional :: scalar

        type(tens_dense_t)         :: tens
        integer(INTL), allocatable :: tens_root(:)
        integer(INTD), allocatable :: tens_id(:)

        ierr=0
        tens=tensor%get_dense_adapter(ierr)
        if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
        rank=tens%num_dims
        
        ierr=setting_zero(this%mode,tens%body_ptr,tens%dims(1:rank),tens%bases(1:rank),rank)

      end function set_zero_apply

      integer function setting_zero(mode,tens_body,tens_dims,tens_bases,rank) result (ierr)

        implicit none

        integer, intent(in)       :: mode
        type(C_PTR), value        :: tens_body       !in: C pointer to the tensor body (tensor elements)
        integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
        integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)
        integer(INTD)             :: rank

        integer                :: p,q,r,s
        complex(8),pointer     :: tens1(:)
        complex(8),pointer     :: tens2(:,:)
        complex(8),pointer     :: tens3(:,:,:)
        complex(8),pointer     :: tens4(:,:,:,:)
        complex(8),pointer     :: tens5(:,:,:,:,:)
        complex(8),pointer     :: tens6(:,:,:,:,:,:)

        ierr=0

        select case (rank)
        case (1)
          call c_f_pointer(tens_body,tens1,tens_dims(1:rank))
          if(mode.eq.3) then
            tens1=ZERO
          else if (mode.eq.1) then
            do p = 1, tens_dims(1)
                tens1(p) = dcmplx(dreal(tens1(p)),0.D0)
            end do
          else if (mode.eq.2) then
            do p = 1, tens_dims(1)
                tens1(p) = dcmplx(0.D0,dimag(tens1(p)))
            end do
          else 
            call quit('Error: wrong mode in set zero')
          end if
        case (2)
          call c_f_pointer(tens_body,tens2,tens_dims(1:rank))
          if(mode.eq.3) then
            tens2=ZERO
          else if (mode.eq.1) then
            do q = 1, tens_dims(2)
              do p = 1, tens_dims(1)
                tens2(p,q) = dcmplx(dreal(tens2(p,q)),0.D0)
              end do
            end do
          else if (mode.eq.2) then
            do q = 1, tens_dims(2)
              do p = 1, tens_dims(1)
                tens2(p,q) = dcmplx(0.D0,dimag(tens2(p,q)))
              end do
            end do
          else 
            call quit('Error: wrong mode in set zero')
          end if
        case (3)
          call c_f_pointer(tens_body,tens3,tens_dims(1:rank))
          if(mode.eq.3) then
            tens3=ZERO
          else if (mode.eq.1) then
            do r = 1, tens_dims(3)
              do q = 1, tens_dims(2)
                do p = 1, tens_dims(1)
                  tens3(p,q,r) = dcmplx(dreal(tens3(p,q,r)),0.D0)
                end do
              end do
            end do
          else if (mode.eq.2) then
            do r = 1, tens_dims(3)
              do q = 1, tens_dims(2)
                do p = 1, tens_dims(1)
                  tens3(p,q,r) = dcmplx(0.D0,dimag(tens3(p,q,r)))
                end do
              end do
            end do
          else 
            call quit('Error: wrong mode in set zero')
          end if
        case (4)
          call c_f_pointer(tens_body,tens4,tens_dims(1:rank))
          if(mode.eq.3) then
            tens4=ZERO
          else if (mode.eq.1) then
            do s = 1, tens_dims(4)
              do r = 1, tens_dims(3)
                do q = 1, tens_dims(2)
                  do p = 1, tens_dims(1)
                    tens4(p,q,r,s) = dcmplx(dreal(tens4(p,q,r,s)),0.D0)
                  end do
                end do
              end do
            end do
          else if (mode.eq.2) then
            do s = 1, tens_dims(4)
              do r = 1, tens_dims(3)
                do q = 1, tens_dims(2)
                  do p = 1, tens_dims(1)
                    tens4(p,q,r,s) = dcmplx(0.D0,dimag(tens4(p,q,r,s)))
                  end do
                end do
              end do
            end do
          else 
            call quit('Error: wrong mode in set zero')
          end if
        case (5)
          call c_f_pointer(tens_body,tens5,tens_dims(1:rank))
          if(mode.eq.3) then
            tens5=ZERO
          else 
            call quit('Error: wrong mode in set zero')
          end if
        case (6)
          call c_f_pointer(tens_body,tens6,tens_dims(1:rank))
          if(mode.eq.3) then
            tens6=ZERO
          else 
            call quit('Error: wrong mode in set zero')
          end if
        case default
          call quit('setting zero: rank not implemented')
        end select
      end function setting_zero

      ! --------------------------------------------------------------------------------

      subroutine delta_t_init(this)
        implicit none
        class(delta_t), intent(out):: this

        this%ind=-9
      end subroutine delta_t_init

      subroutine delta_t_reset(this,val,n)
        implicit none
        class(delta_t), intent(inout) :: this
        integer, intent(in)           :: n
        integer(INTL), intent(in)     :: val(n)

        this%ind=-9
        this%ind(1:n)=val(1:n)

      end subroutine delta_t_reset

      subroutine delta_t_pack(this,packet,ierr)
        implicit none
        class(delta_t), intent(in)           :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr
        integer                              :: i

        do i=1,4
          call pack_builtin(packet,this%ind(i))
        end do

        if(present(ierr)) ierr=0  
        return
      end subroutine delta_t_pack

      subroutine delta_t_unpack(this,packet,ierr)
        implicit none
        class(delta_t), intent(inout)        :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr
        integer                              :: i

        do i=1,4
          call unpack_builtin(packet,this%ind(i))
        end do

        if(present(ierr)) ierr=0
        return
      end subroutine delta_t_unpack

      function delta_t_apply(this,tensor,scalar) result(ierr)

        implicit none
        integer(INTD) :: ierr, rank
        class(delta_t), intent(in) :: this
        class(tens_rcrsv_t), intent(inout) :: tensor
        complex(8), intent(inout), optional :: scalar

        type(tens_dense_t)         :: tens
        integer(INTL), allocatable :: tens_root(:)
        integer(INTD), allocatable :: tens_id(:)

        ierr=0
        tens=tensor%get_dense_adapter(ierr)
        if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
        rank=tens%num_dims
        if (rank.gt.4) call quit('Error: only up to 4 dim. implemented in fill_delta_tensor')
        
        ierr=fill_delta(this%ind,tens%body_ptr,tens%dims(1:rank),tens%bases(1:rank),rank)

      end function delta_t_apply

      integer function fill_delta(ind,tens_body,tens_dims,tens_bases,rank) result (ierr)

        use exacorr_global

        implicit none

        integer(INTL), intent(in) :: ind(4)
        type(C_PTR), value        :: tens_body       !in: C pointer to the tensor body (tensor elements)
        integer(8), intent(in)    :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
        integer(8), intent(in)    :: tens_bases(1:*) !in: base offsets for each tensor dimension)
        integer(INTD)             :: rank

        integer(INTL)          :: hind(1:4)
        integer                :: i
        complex(8),pointer     :: tens1(:)
        complex(8),pointer     :: tens2(:,:)
        complex(8),pointer     :: tens3(:,:,:)
        complex(8),pointer     :: tens4(:,:,:,:)

        ierr=0

        do i=1,rank
          if (ind(i).lt.1) call quit('fill_delta: missing index in fill_delta')
          hind(i)=ind(i)-tens_bases(i)+1
        end do 

        select case (rank)
        case (1)
          call c_f_pointer(tens_body,tens1,tens_dims(1:rank))
          tens1=ZERO
          if (hind(1).gt.0 .and. hind(1).le.tens_dims(1)) tens1(hind(1))=ONE
        case (2)
          call c_f_pointer(tens_body,tens2,tens_dims(1:rank))
          tens2=ZERO
          if (hind(2).gt.0 .and. hind(2).le.tens_dims(2)) then
            if (hind(1).gt.0 .and. hind(1).le.tens_dims(1)) tens2(hind(1),hind(2))=ONE 
          end if
        case (3)
          call c_f_pointer(tens_body,tens3,tens_dims(1:rank))
          tens3=ZERO
          if (hind(3).gt.0 .and. hind(3).le.tens_dims(3)) then
            if (hind(2).gt.0 .and. hind(2).le.tens_dims(2)) then
              if (hind(1).gt.0 .and. hind(1).le.tens_dims(1)) tens3(hind(1),hind(2),hind(3))=ONE
            end if
          end if
        case (4)
          call c_f_pointer(tens_body,tens4,tens_dims(1:rank))
          tens4=ZERO
          if (hind(4).gt.0 .and. hind(4).le.tens_dims(4)) then
            if (hind(3).gt.0 .and. hind(3).le.tens_dims(3)) then
              if (hind(2).gt.0 .and. hind(2).le.tens_dims(2)) then
                if (hind(1).gt.0 .and. hind(1).le.tens_dims(1)) tens4(hind(1),hind(2),hind(3),hind(4))=ONE
              end if
            end if
          end if
        case default
          call quit('fill_delta: wrong rank')
        end select
      end function fill_delta

        ! --------------------------------------------------------------------------------

      function square_t_apply(this,tensor,scalar) result(ierr)

        implicit none
        integer(INTD)                       :: ierr,n
        class(square_t), intent(in)         :: this
        class(tens_rcrsv_t), intent(inout)  :: tensor
        complex(8), intent(inout), optional :: scalar
        type(tens_dense_t)                  :: tens

        ierr=0
        tens=tensor%get_dense_adapter(ierr)
        if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')

        n=tens%num_dims
        if (n.gt.6) call quit('ERROR in square_t_apply: rank too large')
        ierr=square_tensor(tens%num_dims,tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n))

      end function square_t_apply

      integer function square_tensor(rank,tens_body,data_kind,tens_dims,tens_bases) result (ierr)

        implicit none

        integer(INTD),intent(in)    :: rank            !in: number of dimensions of the tensor
        type(C_PTR), value          :: tens_body       !in: C pointer to the tensor body (tensor elements)
        integer(INTD), value        :: data_kind       !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
        integer(INTL), intent(in)   :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
        integer(INTL), intent(in)   :: tens_bases(1:*) !in: base offsets for each tensor dimension)

        integer :: p,q,r,s,t,u
        complex(8),pointer :: t1_tens(:)
        complex(8),pointer :: t2_tens(:,:)
        complex(8),pointer :: t3_tens(:,:,:)
        complex(8),pointer :: t4_tens(:,:,:,:)
        complex(8),pointer :: t6_tens(:,:,:,:,:,:)

        ierr=0 !set error code to success
        if(data_kind /= EXA_DATA_KIND_C8) then
          ierr=-1 !invalid data kind requested
          return
        endif

        select case (rank)
        case (1)
          call c_f_pointer(tens_body,t1_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array
          do p = 1, tens_dims(1)
            t1_tens(p) = conjg(t1_tens(p))*t1_tens(p)
          end do
        case (2)
        call c_f_pointer(tens_body,t2_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array
          do q = 1, tens_dims(2)
            do p = 1, tens_dims(1)
              t2_tens(p,q) = conjg(t2_tens(p,q))*t2_tens(p,q)
            end do
          end do
        case (3)
          call c_f_pointer(tens_body,t3_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array
          do r = 1, tens_dims(3)
            do q = 1, tens_dims(2)
              do p = 1, tens_dims(1)
                t3_tens(p,q,r) = conjg(t3_tens(p,q,r))*t3_tens(p,q,r)
              end do
            end do
          end do
        case (4)
          call c_f_pointer(tens_body,t4_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array
          do s = 1, tens_dims(4)
            do r = 1, tens_dims(3)
              do q = 1, tens_dims(2)
                do p = 1, tens_dims(1)
                  t4_tens(p,q,r,s) = conjg(t4_tens(p,q,r,s))*t4_tens(p,q,r,s)
                end do
              end do
            end do
          end do
        case (6)
          call c_f_pointer(tens_body,t6_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array
          do u = 1, tens_dims(6)
            do t = 1, tens_dims(5)
              do s = 1, tens_dims(4)
                do r = 1, tens_dims(3)
                  do q = 1, tens_dims(2)
                    do p = 1, tens_dims(1)
                      t6_tens(p,q,r,s,t,u) = conjg(t6_tens(p,q,r,s,t,u))*t6_tens(p,q,r,s,t,u)
                    end do
                  end do
                end do
              end do
            end do
          end do
        case default
          call quit ('Error in square: dimension not implemented')
        end select

      end function square_tensor
 
        ! --------------------------------------------------------------------------------

        subroutine set_diagonal_init(this,const,offdiag)

          class(set_diagonal_t), intent(out)   :: this
          complex(8), intent(in) :: const
          logical, intent(in)    :: offdiag

          this%constant=const
          this%offdiag=offdiag

        end subroutine set_diagonal_init

        function set_diagonal_apply(this,tensor,scalar) result(ierr)

         implicit none
         integer(INTD):: ierr,rank
         class(set_diagonal_t), intent(in)    :: this
         class(tens_rcrsv_t), intent(inout)   :: tensor
         complex(8), intent(inout), optional  :: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         rank=tens%num_dims
          
         ierr=diagonal_set(tens%body_ptr,tens%data_kind,tens%dims(1:rank),tens%bases(1:rank), rank, &
                           this%constant,this%offdiag)

        end function set_diagonal_apply

        function diagonal_set(tens_body,data_kind,tens_dims,tens_bases, rank, constant, offdiag) result (ierr)
          ! routine to fill the diagonal of a matrix

          implicit none

          integer(INTD)             :: ierr, rank
          type(C_PTR), value        :: tens_body       !in: C pointer to the tensor body (tensor elements)
          integer(INTD), value      :: data_kind       !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
          integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
          integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)
          complex(8), intent(in)    :: constant
          logical, intent(in)       :: offdiag

          integer            :: po, qo, ro, so, pa, qa, ra, sa
          complex(8),pointer :: h2_tens(:,:), h3_tens(:,:,:), h4_tens(:,:,:,:)
          
          ierr=0

          if (rank==2) then
            call c_f_pointer(tens_body,h2_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array
            if (offdiag) h2_tens=ZERO
            do qo = 1, tens_dims(2) ! loop over columns
              qa = tens_bases(2) + qo - 1
                do ro = 1, tens_dims(1) ! loop over rows
                  ra = tens_bases(1) + ro - 1
                  if (qa.eq.ra) then
                    h2_tens(ro, qo) = constant
                  endif
                end do
            end do

          else if (rank==3) then
            call c_f_pointer(tens_body,h3_tens,tens_dims(1:rank)) !map C pointer to 3-dimensional array
            if (offdiag) h3_tens=ZERO
            do po = 1, tens_dims(3) ! loop over depth
              pa = tens_bases(3) + po - 1
              do qo = 1, tens_dims(2) ! loop over colums
                qa = tens_bases(2) + qo - 1
                if (pa.eq.qa) then
                  do ro = 1, tens_dims(1) ! loop over rows
                    ra = tens_bases(1) + ro - 1
                    if (qa.eq.ra) then
                      h3_tens(ro, qo, po) = constant
                    endif
                  end do
                end if
              end do
            end do

          else if (rank==4) then
            call c_f_pointer(tens_body,h4_tens,tens_dims(1:rank)) !map C pointer to 3-dimensional array
            if (offdiag) h4_tens=ZERO
            do po = 1, tens_dims(4)
              pa = tens_bases(4) + po - 1
              do qo = 1, tens_dims(3)
                qa = tens_bases(3) + qo - 1
                if (pa.eq.qa) then
                  do ro = 1, tens_dims(2)
                    ra = tens_bases(2) + ro - 1
                    if (qa.eq.ra) then
                      do so = 1, tens_dims(1)
                        sa = tens_bases(1) + so - 1
                        if (ra.eq.sa) then
                          h4_tens(ro, qo, po, so) = constant
                        endif
                      end do
                    end if
                  end do
                end if
              end do
            end do

          else
            call quit('Error: only dimensions 2-4 in fill_diagonal')
          end if

        end function diagonal_set

        ! --------------------------------------------------------------------------------

        subroutine set_energy_init(this,mo_occ,nocc,occ_space_id,mo_vir,nvir,vir_space_id,level_shift)
         
         use exacorr_mo
         use exacorr_global
         use exacorr_utils, only : shift_orbital_energy

         implicit none

         class(set_energy_t)      :: this
         integer, intent(in)      :: nocc, nvir
         integer(INTD),intent(in) :: occ_space_id, vir_space_id
         integer, intent(in)      :: mo_occ(:),mo_vir(:)
         real(8), intent(in)      :: level_shift

         type(cmo)   :: cspinor

         this%occ_id=occ_space_id
         this%vir_id=vir_space_id
         this%level_shift=level_shift
         this%norb(1)=nocc
         this%norb(2)=nvir

         call get_mo_coefficients (cspinor,mo_occ,nocc)
         allocate(this%eps_occ(nocc))
         this%eps_occ = cspinor%energy
         call dealloc_mo (cspinor)

         call get_mo_coefficients (cspinor,mo_vir,nvir)
         allocate(this%eps_vir(nvir))
         this%eps_vir = cspinor%energy
         call dealloc_mo (cspinor)

         call shift_orbital_energy(this%eps_vir,this%eps_occ,level_shift)

        end subroutine set_energy_init

      subroutine set_energy_reset(this,eps_occ,eps_vir,nocc,nvir)

        use exacorr_utils, only : shift_orbital_energy

        implicit none

        class(set_energy_t)           :: this
        real(8), intent(in)           :: eps_occ(:),eps_vir(:)
        integer, intent(in)           :: nocc, nvir

        this%eps_occ = eps_occ
        this%eps_vir = eps_vir

        call shift_orbital_energy(this%eps_vir,this%eps_occ,this%level_shift)

      end subroutine set_energy_reset

      subroutine set_energy_pack(this,packet,ierr)
        implicit none
        class(set_energy_t), intent(in)      :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr
        integer                              :: i
        
        do i = 1, this%norb(1)
          call pack_builtin(packet,this%eps_occ(i))
        end do
        do i = 1, this%norb(2)
          call pack_builtin(packet,this%eps_vir(i))
        end do

        if(present(ierr)) ierr=0  
        return
      end subroutine set_energy_pack

      subroutine set_energy_unpack(this,packet,ierr)
        implicit none
        class(set_energy_t), intent(inout)   :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr
        integer                              :: i

        do i = 1, this%norb(1)
          call unpack_builtin(packet,this%eps_occ(i))
        end do
        do i = 1, this%norb(2)
          call unpack_builtin(packet,this%eps_vir(i))
        end do

        if(present(ierr)) ierr=0
        return
      end subroutine set_energy_unpack

        function set_energy_apply(this,tensor,scalar) result(ierr)

          use exacorr_datatypes

          implicit none

          integer(INTD) :: ierr, rank
          class(set_energy_t), intent(in) :: this
          class(tens_rcrsv_t), intent(inout) :: tensor
          complex(8), intent(inout), optional :: scalar

          type(tens_dense_t   )      :: tens
          integer(INTL), allocatable :: dims(:), tens_root(:)
          integer(INTD), allocatable :: tens_id(:)
          integer                    :: space_id

          ierr=0
          tens=tensor%get_dense_adapter(ierr)
          if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
          rank=tens%num_dims
          if (rank.ne.2) call quit('Error: dimension not 2 in set_energy_apply')

          allocate(dims(rank))
          call tensor%get_dims(dims,rank)

          allocate (tens_id (rank))
          allocate (tens_root (rank))
          call tensor%get_space_ids(tens_id,tens_root,ierr)

          if (tens_id(1) == this%occ_id) then
            if (tens_id(2) == this%occ_id) then
              space_id=this%occ_id ! tens_id 1 corresponds to occupied space
            else 
              call quit('incompatible space for second index in set_energy_apply')
            end if
          else if (tens_id(1) == this%vir_id) then
            if (tens_id(2) == this%vir_id) then
             space_id=this%vir_id! tens_id 1 corresponds to occupied space
            else 
              call quit('incompatible space for second index in set_energy_apply')
            end if
          else
            call quit('incompatible space id in set_energy_apply')
          end if

          ierr=energy_set_diag(rank,tens%body_ptr,tens%dims(1:rank),tens%bases(1:rank), space_id, &
                               this%eps_occ, this%eps_vir, this%occ_id, this%vir_id)

          deallocate(dims)
          deallocate(tens_id)
          deallocate(tens_root)
        end function set_energy_apply

        integer function energy_set_diag(rank,tens_body,tens_dims,tens_bases, space_id, &
                              eps_occ, eps_vir, occ_id, vir_id) result (ierr)
          !Initialization routine to fill tensor diagonals with eigenvalues form dfcoef

          use exacorr_datatypes

          implicit none
          integer(INTD)             :: rank
          type(C_PTR), value        :: tens_body       !in: C pointer to the tensor body (tensor elements)
          integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
          integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)

          real(8), intent(in)   :: eps_occ(:), eps_vir(:)
          integer, intent(in)   :: space_id, occ_id, vir_id
          complex(8),pointer    :: i_tens(:,:)
          integer(INTL)         :: ko, lo, kg, lg

          ierr=0
          call c_f_pointer(tens_body,i_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array

          i_tens = ZERO
          do lo = 1, tens_dims(2) 
             lg = tens_bases(2) + lo - 1
             if (lg <  tens_bases(1)) cycle ! we are in an off-diagonal block, kg will always be larger than lg
             if (lg >= tens_bases(1)+tens_dims(1)) exit  ! we are in an off-diagonal block, lg always larger than kg
             do ko = 1, tens_dims(1) 
                kg = tens_bases(1) + ko - 1

                if (kg == lg) then
                  if (space_id.eq.occ_id) i_tens(ko,lo) = eps_occ(lg)
                  if (space_id.eq.vir_id) i_tens(ko,lo) = eps_vir(lg)
                end if
             end do
          end do

        end function energy_set_diag

        ! --------------------------------------------------------------------------------
      
      subroutine denom_t_init(this,mo_occ,nocc,mo_vir,nvir,level_shift,occ_id)

        use exacorr_mo
        use exacorr_global
        use exacorr_utils, only : shift_orbital_energy

        implicit none

        integer, intent(in)          :: nocc, nvir
        integer, intent(in)          :: mo_occ(:),mo_vir(:)
        real(8), intent(in)          :: level_shift
        integer(INTD), intent(in)    :: occ_id
        class(denom_t), intent(out)  :: this
         
        type(cmo)          :: cspinor

        this%occ_id=occ_id
        this%level_shift=level_shift
        this%norb(1)=nocc
        this%norb(2)=nvir

        call get_mo_coefficients (cspinor,mo_occ,nocc)
        allocate(this%eps_occ(nocc))
        this%eps_occ = cspinor%energy
        call dealloc_mo (cspinor)

        call get_mo_coefficients (cspinor,mo_vir,nvir)
        allocate(this%eps_vir(nvir))
        this%eps_vir = cspinor%energy
        call dealloc_mo (cspinor)

        call shift_orbital_energy(this%eps_vir,this%eps_occ,level_shift)

      end subroutine denom_t_init

      subroutine denom_t_reset(this,eps_occ,eps_vir,nocc,nvir)

        use exacorr_utils, only : shift_orbital_energy

        implicit none

        class(denom_t), intent(inout) :: this
        real(8), intent(in)           :: eps_occ(:),eps_vir(:)
        integer, intent(in)           :: nocc, nvir

        if (nocc.ne.this%norb(1)) call quit('denom_t_reset: wrong number of orbitals')
        if (nvir.ne.this%norb(2)) call quit('denom_t_reset: wrong number of orbitals')

        this%eps_occ = eps_occ
        this%eps_vir = eps_vir

        call shift_orbital_energy(this%eps_vir,this%eps_occ,this%level_shift)

      end subroutine denom_t_reset

      subroutine denom_t_pack(this,packet,ierr)
        implicit none
        class(denom_t), intent(in)           :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr
        integer                              :: i
        
        do i = 1, this%norb(1)
          call pack_builtin(packet,this%eps_occ(i))
        end do
        do i = 1, this%norb(2)
          call pack_builtin(packet,this%eps_vir(i))
        end do

        if(present(ierr)) ierr=0  
        return
      end subroutine denom_t_pack

      subroutine denom_t_unpack(this,packet,ierr)
        implicit none
        class(denom_t), intent(inout)        :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr
        integer                              :: i

        do i = 1, this%norb(1)
          call unpack_builtin(packet,this%eps_occ(i))
        end do
        do i = 1, this%norb(2)
          call unpack_builtin(packet,this%eps_vir(i))
        end do

        if(present(ierr)) ierr=0
        return
      end subroutine denom_t_unpack

      function denom_t_apply(this,tensor,scalar) result(ierr)
        implicit none
        integer(INTD)                       :: ierr,n
        class(denom_t), intent(in)          :: this
        class(tens_rcrsv_t), intent(inout)  :: tensor
        complex(8), intent(inout), optional :: scalar
        type(tens_dense_t)                  :: tens
        integer(INTD), allocatable          :: tens_id(:)
        integer(INTL), allocatable          :: tens_root(:)
        logical                             :: lambda=.false.

        ierr=0
        tens=tensor%get_dense_adapter(ierr)
        if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')

        select case (tens%num_dims)
        case (2)
          allocate (tens_id  (2))
          allocate (tens_root(2))
        case (4)
          allocate (tens_id  (4))
          allocate (tens_root(4))
        case (6)
          allocate (tens_id  (6))
          allocate (tens_root(6))
        case default
          call quit ('Error: wrong dimension of tensor in denom_t_apply')
        end select

        call tensor%get_space_ids(tens_id,tens_root,ierr)
        if (tens_id(1) == this%occ_id) lambda = .true. ! tens_id corresponds to occupied space

        ierr=denominate(tens%num_dims,tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n), &
                         this%eps_occ,this%eps_vir,lambda)

        deallocate(tens_id)
        deallocate(tens_root)
      end function denom_t_apply

      integer function denominate(rank,tens_body,data_kind,tens_dims,tens_bases, &
                         eps_occ,eps_vir,lambda) result (ierr)

        !Scale with denominators for 2-d or 4-d tensor for regular CC and for lambda
        !                        and for 3-d and 6-d tensor for triples

        implicit none

        integer(INTD),intent(in)  :: rank            !in: number of dimensions of the tensor
        type(C_PTR), value        :: tens_body       !in: C pointer to the tensor body (tensor elements)
        integer(INTD), value      :: data_kind       !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
        integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
        integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)
        real(8), intent(in)       :: eps_occ(:),eps_vir(:)
        logical, intent(in)       :: lambda

        integer(INTL)      :: i,j,k,a,b,c
        real(8)            :: denominator
        complex(8),pointer :: t1_tens(:,:)
        complex(8),pointer :: t2_tens(:,:,:,:)
        complex(8),pointer :: t3_tens(:,:,:,:,:,:)

        ierr=0 !set error code to success
        if(data_kind /= EXA_DATA_KIND_C8) then
          ierr=-1 !invalid data kind requested
          return
        endif

        select case (rank)
        case (2)
          call c_f_pointer(tens_body,t1_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array
          if (.not.lambda) then 
            do i = 1, tens_dims(2)
              do a = 1, tens_dims(1)
                denominator = eps_occ(i+tens_bases(2)-1) - eps_vir(a+tens_bases(1)-1)
                t1_tens(a,i) = t1_tens(a,i) / denominator
              end do
            end do
          else
            do i = 1, tens_dims(1)
              do a = 1, tens_dims(2)
                denominator = eps_occ(i+tens_bases(1)-1) - eps_vir(a+tens_bases(2)-1)
                t1_tens(i,a) = t1_tens(i,a) / denominator
              end do
            end do
          end if
        case (4)
          call c_f_pointer(tens_body,t2_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array
          if (.not.lambda) then
            do j = 1, tens_dims(4)
              do i = 1, tens_dims(3)
                do b = 1, tens_dims(2)
                  do a = 1, tens_dims(1)
                    denominator = eps_occ(i+tens_bases(3)-1) + eps_occ(j+tens_bases(4)-1) &
                                - eps_vir(a+tens_bases(1)-1) - eps_vir(b+tens_bases(2)-1)
                    t2_tens(a,b,i,j) = t2_tens(a,b,i,j) / denominator
                  end do
                end do
              end do
            end do
          else
            do j = 1, tens_dims(2)
              do i = 1, tens_dims(1)
                do b = 1, tens_dims(4)
                  do a = 1, tens_dims(3)
                    denominator = eps_occ(i+tens_bases(1)-1) + eps_occ(j+tens_bases(2)-1) &
                                - eps_vir(a+tens_bases(3)-1) - eps_vir(b+tens_bases(4)-1)
                    t2_tens(i,j,a,b) = t2_tens(i,j,a,b) / denominator
                  end do
                end do
              end do
            end do
          end if
        case (6)
          call c_f_pointer(tens_body,t3_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array
          if (.not.lambda) then
            do k = 1, tens_dims(6)
              do j = 1, tens_dims(5)
                do i = 1, tens_dims(4)
                  do c = 1, tens_dims(3)
                    do b = 1, tens_dims(2)
                      do a = 1, tens_dims(1)
                        denominator = eps_occ(k+tens_bases(6)-1) + eps_occ(j+tens_bases(5)-1) &
                                    + eps_occ(i+tens_bases(4)-1) - eps_vir(c+tens_bases(3)-1) &
                                    - eps_vir(b+tens_bases(2)-1) - eps_vir(a+tens_bases(1)-1)
                        t3_tens(a,b,c,i,j,k) = t3_tens(a,b,c,i,j,k) / denominator
                      end do 
                    end do
                  end do
                end do
              end do
            end do
          else
            do k = 1, tens_dims(3)
              do j = 1, tens_dims(2)
                do i = 1, tens_dims(1)
                  do c = 1, tens_dims(6)
                    do b = 1, tens_dims(5)
                      do a = 1, tens_dims(4)
                        denominator = eps_occ(k+tens_bases(3)-1) + eps_occ(j+tens_bases(2)-1) &
                                    + eps_occ(i+tens_bases(1)-1) - eps_vir(c+tens_bases(6)-1) &
                                    - eps_vir(b+tens_bases(5)-1) - eps_vir(a+tens_bases(4)-1)
                        t3_tens(i,j,k,a,b,c) = t3_tens(i,j,k,a,b,c) / denominator
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end if
         case default
            call quit ('Error: wrong dimension of tensor in denominate')
         end select

      end function denominate

        ! --------------------------------------------------------------------------------

      subroutine denom3_t_init(this,mo_occ,nocc,mo_vir,nvir,level_shift)

        use exacorr_mo
        use exacorr_global
        use exacorr_utils, only : shift_orbital_energy

        implicit none

        class(denom3_t), intent(out) :: this
        integer, intent(in)          :: nocc, nvir
        integer, intent(in)          :: mo_occ(:),mo_vir(:)
        real(8), intent(in)          :: level_shift
         
        type(cmo)          :: cspinor

        this%level_shift=level_shift
        this%norb(1)=nocc
        this%norb(2)=nvir

        call get_mo_coefficients (cspinor,mo_occ,nocc)
        allocate(this%eps_occ(nocc))
        this%eps_occ = cspinor%energy
        call dealloc_mo (cspinor)

        call get_mo_coefficients (cspinor,mo_vir,nvir)
        allocate(this%eps_vir(nvir))
        this%eps_vir = cspinor%energy
        call dealloc_mo (cspinor)

        call shift_orbital_energy(this%eps_vir,this%eps_occ,level_shift)

        this%eps_ijk=1.0D0

      end subroutine denom3_t_init

      subroutine denom3_t_reset(this,val)
        implicit none
        class(denom3_t), intent(inout) :: this
        real(8), intent(in)           :: val

        this%eps_ijk=val

      end subroutine denom3_t_reset

      subroutine denom3_t_pack(this,packet,ierr)
        implicit none
        class(denom3_t), intent(in)           :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr

        call pack_builtin(packet,this%eps_ijk)

        if(present(ierr)) ierr=0  
        return
      end subroutine denom3_t_pack

      subroutine denom3_t_unpack(this,packet,ierr)
        implicit none
        class(denom3_t), intent(inout)        :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr

        call unpack_builtin(packet,this%eps_ijk)

        if(present(ierr)) ierr=0
        return
      end subroutine denom3_t_unpack

      function denom3_t_apply(this,tensor,scalar) result(ierr)
        implicit none
        integer(INTD)                       :: ierr,n
        class(denom3_t), intent(in)          :: this
        class(tens_rcrsv_t), intent(inout)  :: tensor
        complex(8), intent(inout), optional :: scalar
        type(tens_dense_t)                  :: tens
        integer(INTD), allocatable          :: tens_id(:)
        integer(INTL), allocatable          :: tens_root(:)

        ierr=0
        tens=tensor%get_dense_adapter(ierr)
        if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')

        select case (tens%num_dims)
        case (3)
          allocate (tens_id  (3))
          allocate (tens_root(3))
          if (this%eps_ijk.gt.0) call quit ('Error: eps_ijk not set in denom3_t_apply')
        case default
          call quit ('Error: wrong dimension of tensor in denom3_t_apply')
        end select

        ierr=denominate3(tens%num_dims,tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n), &
                         this%eps_vir,this%eps_ijk)

        deallocate(tens_id)
        deallocate(tens_root)
      end function denom3_t_apply

      integer function denominate3(rank,tens_body,data_kind,tens_dims,tens_bases, &
                         eps_vir,eps_ijk) result (ierr)

        !Scale with denominators for 2-d or 4-d tensor for regular CC and for lambda
        !                        and for 3-d and 6-d tensor for triples

        implicit none

        integer(INTD),intent(in)  :: rank            !in: number of dimensions of the tensor
        type(C_PTR), value        :: tens_body       !in: C pointer to the tensor body (tensor elements)
        integer(INTD), value      :: data_kind       !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
        integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
        integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)
        real(8), intent(in)       :: eps_vir(:)
        real(8), intent(in)       :: eps_ijk

        integer(INTL)      :: a,b,c
        real(8)            :: denominator
        complex(8),pointer :: th_tens(:,:,:)

        ierr=0 !set error code to success
        if(data_kind /= EXA_DATA_KIND_C8) then
          ierr=-1 !invalid data kind requested
          return
        endif

        select case (rank)
        case (3)
          call c_f_pointer(tens_body,th_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array
          do c = 1, tens_dims(3)
            do b = 1, tens_dims(2)
              do a = 1, tens_dims(1)
                denominator = eps_ijk - eps_vir(a+tens_bases(1)-1) &
                            - eps_vir(b+tens_bases(2)-1) - eps_vir(c+tens_bases(3)-1) 
                th_tens(a,b,c) = th_tens(a,b,c) / denominator
              end do
            end do
          end do
         case default
            call quit ('Error: wrong dimension of tensor in denominate')
         end select

      end function denominate3

        ! --------------------------------------------------------------------------------

        subroutine set_1el_init(this,mo1_list,nmo1,mo2_list,nmo2, min_occ)

          use exacorr_global

          implicit none

          class(set_1el_t), intent(out):: this
          integer, intent(in) :: nmo1          ! the length of the mo basis
          integer, intent(in) :: mo1_list(:)   ! and their indices
          integer, intent(in) :: nmo2          ! the length of the mo basis
          integer, intent(in) :: mo2_list(:)   ! and their indices
          integer, intent(in) :: min_occ       ! smallest index

          
          logical        :: one_el_exist  ! we also have one-electron Hamiltonian elements available
          real(8)        :: e_core

          one_el_exist = exist_one_el(e_core)

          if (one_el_exist) then
            ! store the integrals inside this type, so we can extract them later
            allocate(this%integrals(nmo1,nmo2))
            call get_one_el(this%integrals, mo1_list,nmo1,mo2_list,nmo2, min_occ)

          end if

        end subroutine set_1el_init

        function set_1el_apply(this,tensor,scalar) result(ierr)

          use exacorr_datatypes

          implicit none

          integer(INTD)                       :: ierr, rank
          class(set_1el_t), intent(in)        :: this
          class(tens_rcrsv_t), intent(inout)  :: tensor
          complex(8), intent(inout), optional :: scalar

          type(tens_dense_t)         :: tens
          integer(INTL), allocatable :: dims(:)

          ierr=0
          tens=tensor%get_dense_adapter(ierr)
          if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
          rank=tens%num_dims
          if (rank.ne.2) call quit('Error: dimension not 2 in one_el_apply')

          allocate(dims(rank))
          call tensor%get_dims(dims,rank)

          ierr=one_el_set(tens%body_ptr,tens%data_kind,tens%dims(1:rank),tens%bases(1:rank), rank, &
                               this%integrals)
          deallocate(dims)
        end function set_1el_apply

        integer function one_el_set(tens_body,data_kind,tens_dims,tens_bases, rank, &
                         integrals) result (ierr)
          !Initialization routine to fill tensor with one electron integrals

          use exacorr_datatypes

          implicit none

          type(C_PTR), value        :: tens_body           !in: C pointer to the tensor body (tensor elements)
          integer(INTD), value      :: data_kind, rank      !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
          integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
          integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)
          complex(8), intent(in)    :: integrals(:,:)

          complex(8),pointer :: i_tens(:,:)
          integer(INTL)      :: lo, lg, ko, kg


          ierr=0
          call c_f_pointer(tens_body,i_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array


          do lo = 1, tens_dims(2) 
             lg = tens_bases(2) + lo - 1
             do ko = 1, tens_dims(1) 
                kg = tens_bases(1) + ko - 1

                i_tens(ko,lo) = integrals(kg, lg)

             end do
          end do

        end function one_el_set

        ! --------------------------------------------------------------------------------

      subroutine set_ff_init(this,mo_occ,nocc,occ_space_id,mo_vir,nvir,vir_space_id,min_occ)
          implicit none
          class(set_ff_t)          :: this
          integer, intent(in)      :: nocc, nvir, min_occ
          integer(INTD),intent(in) :: occ_space_id, vir_space_id
          integer, intent(in)      :: mo_occ(:),mo_vir(:)

          this%occ_id=occ_space_id
          this%vir_id=vir_space_id
          this%nocc=nocc
          this%nvir=nvir
          this%min_occ=min_occ
          this%i_prop=-1

          allocate(this%mo_occ(nocc))
          this%mo_occ = mo_occ

          allocate(this%mo_vir(nvir))
          this%mo_vir = mo_vir

      end subroutine set_ff_init

      subroutine set_ff_reset(this,i_prop)
        implicit none
        class(set_ff_t)         :: this
        integer, intent(in)     :: i_prop

        this%i_prop=i_prop

      end subroutine set_ff_reset

      subroutine set_ff_pack(this,packet,ierr)
        implicit none
        class(set_ff_t), intent(in)          :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr
        
        call pack_builtin(packet,this%i_prop)

        if(present(ierr)) ierr=0  
        return
      end subroutine set_ff_pack

      subroutine set_ff_unpack(this,packet,ierr)
        implicit none
        class(set_ff_t), intent(inout)       :: this
        class(obj_pack_t), intent(inout)     :: packet
        integer(INTD), intent(out), optional :: ierr

        call unpack_builtin(packet,this%i_prop)

        if(present(ierr)) ierr=0
        return
      end subroutine set_ff_unpack

      function set_ff_apply(this,tensor,scalar) result(ierr)

          use exacorr_global

          implicit none

          integer(INTD)                       :: ierr, rank
          class(set_ff_t), intent(in)         :: this
          class(tens_rcrsv_t), intent(inout)  :: tensor
          complex(8), intent(inout), optional :: scalar

          type(tens_dense_t   )      :: tens
          integer(INTL), allocatable :: dims(:), tens_root(:)
          integer(INTD), allocatable :: tens_id(:)
          integer                    :: space_id
          complex(8), allocatable    :: prop(:,:) 

          ierr=0
          tens=tensor%get_dense_adapter(ierr)
          if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
          rank=tens%num_dims
          if (rank.ne.2) call quit('Error: only 2-dimensional property')

          allocate(dims(rank))
          call tensor%get_dims(dims,rank)

          allocate (tens_id (rank))
          allocate (tens_root (rank))
          call tensor%get_space_ids(tens_id,tens_root,ierr)

          if (tens_id(1) == this%occ_id) then
            if (tens_id(2) == this%occ_id) then
              allocate(prop(this%nocc,this%nocc))
              call get_ff_mat (prop,this%i_prop,this%mo_occ,this%nocc, &
                                                this%mo_occ,this%nocc,this%min_occ)
            else 
              allocate(prop(this%nocc,this%nvir))
              call get_ff_mat (prop,this%i_prop,this%mo_occ,this%nocc, &
                                                this%mo_vir,this%nvir,this%min_occ)
            end if
          else
            if (tens_id(2) == this%occ_id) then
              allocate(prop(this%nvir,this%nocc))
              call get_ff_mat (prop,this%i_prop,this%mo_vir,this%nvir, &
                                                this%mo_occ,this%nocc,this%min_occ)
            else 
              allocate(prop(this%nvir,this%nvir))
              call get_ff_mat (prop,this%i_prop,this%mo_vir,this%nvir, &
                                                this%mo_vir,this%nvir,this%min_occ)
            end if
          end if

          ierr=property_set(rank,tens%body_ptr,tens%dims(1:rank),tens%bases(1:rank), &
                               prop)

          deallocate(dims)
          deallocate(tens_id)
          deallocate(tens_root)
          deallocate(prop)
      end function set_ff_apply

      function property_set(rank,tens_body,tens_dims,tens_bases, &
                              prop) result (ierr)

          implicit none

          integer(INTD)             :: ierr, rank
          type(C_PTR), value        :: tens_body       !in: C pointer to the tensor body (tensor elements)
          integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
          integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)

          complex(8), allocatable   :: prop(:,:) 

          complex(8),pointer    :: p_tens(:,:)
          integer(INTL)         :: po, qo, pg, qg


          ierr=0
          call c_f_pointer(tens_body,p_tens,tens_dims(1:rank)) !map C pointer to 2-dimensional array

          p_tens = ZERO
          do po = 1, tens_dims(2) 
            pg = tens_bases(2) + po - 1
            do qo = 1, tens_dims(1) 
              qg = tens_bases(1) + qo - 1
              p_tens(qo,po) = prop(qg,pg)
            end do
          end do

      end function property_set

        ! --------------------------------------------------------------------------------
       
        subroutine set_projection_init(this, const)

          class(set_projection_t), intent(out)   :: this
          complex(8), intent(in) :: const

          this%constant=const

        end subroutine set_projection_init

        function set_projection_apply(this,tensor,scalar) result(ierr)

         implicit none
         integer(INTD):: ierr,rank
         class(set_projection_t), intent(in)  :: this
         class(tens_rcrsv_t), intent(inout)   :: tensor
         complex(8), intent(inout), optional  :: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         rank=tens%num_dims

         if (rank.ne.4) call quit('Error: omly dimension 4 in set_projection_apply')
          
         ierr=projection_set(tens%body_ptr,tens%data_kind,tens%dims(1:rank),tens%bases(1:rank), rank, this%constant)

        end function set_projection_apply

        function projection_set(tens_body,data_kind,tens_dims,tens_bases, rank, constant) result (ierr)
          ! routine to fill the diagonal of a matrix

          implicit none

          integer(INTD)             :: ierr, rank
          type(C_PTR), value        :: tens_body       !in: C pointer to the tensor body (tensor elements)
          integer(INTD), value      :: data_kind       !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
          integer(INTL), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
          integer(INTL), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)
          complex(8), intent(in)    :: constant

          integer            :: po, qo, ro, so, pa, qa, ra, sa
          complex(8),pointer :: h4_tens(:,:,:,:)
          
          ierr=0

          if (rank==4) then
            call c_f_pointer(tens_body,h4_tens,tens_dims(1:rank)) !map C pointer to 3-dimensional array
          
            h4_tens=ZERO

            do po = 1, tens_dims(4)
              pa = tens_bases(4) + po - 1
              do qo = 1, tens_dims(3)
                qa = tens_bases(3) + qo - 1
                if (pa.eq.qa) then
                  do ro = 1, tens_dims(2)
                    ra = tens_bases(2) + ro - 1
                    do so = 1, tens_dims(1)
                      sa = tens_bases(1) + so - 1
                      if (sa.eq.ra) then
                        h4_tens(so, ro, qo, po) = constant
                      endif
                    end do
                  end do
                end if
              end do
            end do

          else
            call quit('Error: only dimension 4 in projection_set')
          end if

        end function projection_set

        ! --------------------------------------------------------------------------------

        function chol_diag_apply(this,tensor,scalar) result(ierr)
         implicit none
         integer(INTD):: ierr,rank
         class(chol_diag), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error in comp_chol_diag: getting dense_adapter')
         rank=tens%num_dims
         if (rank.ne.2) call quit('Error in comp_chol_diag: dimension not 2')
         ierr=get_chol_diag(tens%body_ptr,tens%data_kind,tens%dims(1:rank),tens%bases(1:rank),rank)
         return
        end function chol_diag_apply

        integer function get_chol_diag(tens_body, data_kind, tens_dims, tens_bases,rank) result (ierr)
          !Initialization routine that is used to fill a tensor with AO integrals

          use exacorr_datatypes
          use exacorr_global
          use module_interest_eri

          implicit none

          type(C_PTR), value        :: tens_body       !in: C pointer to the tensor body (tensor elements)
          integer(4), value         :: data_kind       !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
          integer(8), intent(in)    :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
          integer(8), intent(in)    :: tens_bases(1:*) !in: base offsets for each tensor dimension)
          integer(INTD), intent(in) :: rank            ! number of independent indices

          integer             :: nao(rank)
          integer             :: basis_angular ! 1=cartesian, 2=spherical
          complex(8),pointer  :: ao_tens(:,:)
          integer             :: first_ao(1:rank)

          integer                              :: nshells(rank)
          type(basis_func_info_t), allocatable :: gto1(:), gto2(:)  ! arrays with basis function information

          integer            :: i, j, ish, jsh, ioff, joff, ijij
          integer            :: li,lj,ni,nj
          real(8)            :: ei,ej,ci,cj,xi,yi,zi,xj,yj,zj
          real(8),save       :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
          !$OMP threadprivate(gout)
          real(8), parameter :: FIJKL = 1.0

          ierr=0 !set error code to success
          if(data_kind /= EXA_DATA_KIND_C8) then
            ierr=-1 !invalid data kind requested
            return
          endif

!         Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
          call interest_initialize(.false.)

          ! map the arguments to names used in the context of integral evaluation
          first_ao(1:rank) = tens_bases(1:rank)
          nao(1:rank) = tens_dims(1:rank)
          call c_f_pointer(tens_body,ao_tens,nao(1:rank)) !map C pointer to array

!         Allocate and get shell information for all shells that we will treat
          call get_gtos(first_ao(1),nao(1),gto1,nshells(1))
          call get_gtos(first_ao(2),nao(2),gto2,nshells(2))
          basis_angular = get_basis_angular()

!         Intranode parallelization should be possible in the loops below because Interest is thread safe

          joff = 0
          do jsh = 1, nshells(2)
            lj   =  gto2(jsh)%orb_momentum
            ej   =  gto2(jsh)%exponent(1)
            xj   =  gto2(jsh)%coord(1)
            yj   =  gto2(jsh)%coord(2)
            zj   =  gto2(jsh)%coord(3)
            cj   =  gto2(jsh)%coefficient(1)
            nj   =  nfunctions(lj,basis_angular)
            ioff = 0
            do ish = 1, nshells(1)
              li   =  gto1(ish)%orb_momentum
              ei   =  gto1(ish)%exponent(1)
              xi   =  gto1(ish)%coord(1)
              yi   =  gto1(ish)%coord(2)
              zi   =  gto1(ish)%coord(3)
              ci   =  gto1(ish)%coefficient(1)
              ni   =  nfunctions(li,basis_angular)

!             output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
              call interest_eri('llll',FIJKL,gout,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )

              do j = 1, nj
                do i = 1, ni
                  ijij = (j-1)*ni*nj*ni+(i-1)*nj*ni+(j-1)*ni+i
                  ao_tens(i+ioff,j+joff) = dcmplx(gout(ijij),0.D0)
                end do
              end do

              ioff = ioff + ni
            end do
            joff = joff + nj
          end do

          deallocate(gto1)
          deallocate(gto2)

        end function get_chol_diag

        ! --------------------------------------------------------------------------------

        subroutine chol_J_init(this,val)
         implicit none
         class(chol_J), intent(out):: this
         integer(INTL), intent(in):: val(1:2)

         this%ind(1:2)=val(1:2)
        end subroutine chol_J_init

        subroutine chol_J_reset(this,val)
         implicit none
         class(chol_J), intent(inout):: this
         integer(INTL), intent(in):: val(1:2)

         this%ind(1:2)=val(1:2)
        end subroutine chol_J_reset

        subroutine chol_J_pack(this,packet,ierr)
         implicit none
         class(chol_J), intent(in):: this
         class(obj_pack_t), intent(inout):: packet
         integer(INTD), intent(out), optional:: ierr

         call pack_builtin(packet,this%ind(1))
         call pack_builtin(packet,this%ind(2))
         if(present(ierr)) ierr=0
         return
        end subroutine chol_J_pack

        subroutine chol_J_unpack(this,packet,ierr)
         implicit none
         class(chol_J), intent(inout):: this
         class(obj_pack_t), intent(inout):: packet
         integer(INTD), intent(out), optional:: ierr

         call unpack_builtin(packet,this%ind(1))
         call unpack_builtin(packet,this%ind(2))
         if(present(ierr)) ierr=0
         return
        end subroutine chol_J_unpack

        function chol_J_apply(this,tensor,scalar) result(ierr)
         implicit none
         integer(INTD):: ierr,rank
         class(chol_J), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error in comp_chol_J: getting dense_adapter')
         rank=tens%num_dims
         if (rank.ne.2) call quit('Error in comp_chol_J: dimension not 2')
         
         ierr=get_chol_J(this%ind,tens%body_ptr,tens%data_kind,tens%dims(1:rank),tens%bases(1:rank),rank)
         return
        end function chol_J_apply

        integer function get_chol_J(ind,tens_body, data_kind, tens_dims, tens_bases, rank) result (ierr)
        !Initialization routine that is used to fill a tensor with AO integrals

         use exacorr_datatypes
         use exacorr_global
         use module_interest_eri

         implicit none

         integer(INTL), intent(in) :: ind(1:2)
         type(C_PTR), value        :: tens_body       !in: C pointer to the tensor body (tensor elements)
         integer(4), value         :: data_kind       !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(8), intent(in)    :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(8), intent(in)    :: tens_bases(1:*) !in: base offsets for each tensor dimension)
         integer(INTD), intent(in) :: rank            ! number of independent indices

         integer:: nao(4)
         integer:: basis_angular ! 1=cartesian, 2=spherical
         complex(8),pointer :: ao_tens(:,:)
         integer :: first_ao(1:4)


         integer :: nshells(4)
         type(basis_func_info_t), allocatable :: gto1(:), gto2(:), gto3(:), gto4(:)  ! arrays with basis function information

         integer :: i, j, k, l, ish, jsh, ksh, lsh, ioff, joff, koff, loff, ijkl
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         !$OMP threadprivate(gout)
         real(8), parameter :: FIJKL = 1.0

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

!        Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
         call interest_initialize(.false.)

         ! map the arguments to names used in the context of integral evaluation
         first_ao(1:rank) = tens_bases(1:rank)
         first_ao(3:4)=ind
         nao(1:rank) = tens_dims(1:rank)
         nao(3:4)=1 ! fixed single element (maybe use shell later)
         call c_f_pointer(tens_body,ao_tens,tens_dims(1:rank)) !map C pointer to 4-dimensional array

!        Allocate and get shell information for all shells that we will treat
         call get_gtos(first_ao(1),nao(1),gto1,nshells(1))
         call get_gtos(first_ao(2),nao(2),gto2,nshells(2))
         call get_gtos(first_ao(3),nao(3),gto3,nshells(3))
         call get_gtos(first_ao(4),nao(4),gto4,nshells(4))
         basis_angular = get_basis_angular()

!        Intranode parallelization should be possible in the loops below because Interest is thread safe

         lsh = 1
            ll   =  gto4(lsh)%orb_momentum
            el   =  gto4(lsh)%exponent(1)
            xl   =  gto4(lsh)%coord(1)
            yl   =  gto4(lsh)%coord(2)
            zl   =  gto4(lsh)%coord(3)
            cl   =  gto4(lsh)%coefficient(1)
            nl   =  nfunctions(ll,basis_angular)

            ksh = 1
               lk   =  gto3(ksh)%orb_momentum
               ek   =  gto3(ksh)%exponent(1)
               xk   =  gto3(ksh)%coord(1)
               yk   =  gto3(ksh)%coord(2)
               zk   =  gto3(ksh)%coord(3)
               ck   =  gto3(ksh)%coefficient(1)
               nk   =  nfunctions(lk,basis_angular)
               joff = 0
               do jsh = 1, nshells(2)
                  lj   =  gto2(jsh)%orb_momentum
                  ej   =  gto2(jsh)%exponent(1)
                  xj   =  gto2(jsh)%coord(1)
                  yj   =  gto2(jsh)%coord(2)
                  zj   =  gto2(jsh)%coord(3)
                  cj   =  gto2(jsh)%coefficient(1)
                  nj   =  nfunctions(lj,basis_angular)
                  ioff = 0
                  do ish = 1, nshells(1)
                     li   =  gto1(ish)%orb_momentum
                     ei   =  gto1(ish)%exponent(1)
                     xi   =  gto1(ish)%coord(1)
                     yi   =  gto1(ish)%coord(2)
                     zi   =  gto1(ish)%coord(3)
                     ci   =  gto1(ish)%coefficient(1)
                     ni   =  nfunctions(li,basis_angular)
!                    output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
                     call interest_eri('llll',FIJKL,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )
                     
                     call get_shell_offset(first_ao(4),loff)
                     l = first_ao(4)-loff
                        call get_shell_offset(first_ao(3),koff)
                        k = first_ao(3)-koff
                           do j = 1, nj
                              do i = 1, ni
                                 ijkl = (l-1)*nk*nj*ni+(k-1)*nj*ni+(j-1)*ni+i
                                 ao_tens(i+ioff,j+joff) = dcmplx(gout(ijkl),0.D0)
                              end do
                           end do
                     ioff = ioff + ni
                  end do
                  joff = joff + nj
               end do

         deallocate(gto1)
         deallocate(gto2)
         deallocate(gto3)
         deallocate(gto4)

        end function get_chol_J

        ! --------------------------------------------------------------------------------
        ! --------------------------------------------------------------------------------
#endif

end module exacorr_tensor_methods

