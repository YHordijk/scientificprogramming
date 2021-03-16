module talsh_ao_to_mo

!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.


        use tensor_algebra
        use talsh
        use exacorr_utils
        use mobasis_hartree_fock
        use, intrinsic:: ISO_C_BINDING

        implicit none

        private
        public talsh_compute_mo_integrals
        public compute_ao_integrals
        public get_integral_tensor
        public tensor_norm2
        public ao2mo_ind
        public ao2mo_vec
        public mulliken_to_dirac_sort
        public decompose_cholesky_talsh_vec
        public decompose_cholesky_talsh_all

        interface compute_ao_integrals
           module procedure compute_ao_integrals
        end interface
        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

       contains
!-----------------------------------------------------------
        subroutine talsh_compute_mo_integrals (nmo, mo_list, talsh_buff, print_level, rcw, th_cholesky)

!This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor
!Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code
!The integrals are written to MDCINT, to allow for testing with unmodified DIRAC correlation modules

!        Written by Lucas Visscher, winter 2016/2017

         use exacorr_mo
         use exacorr_datatypes
         use exacorr_global
         use module_interest_interface
         use module_interest_eri
#if defined (VAR_MPI)
         use interface_to_mpi
#endif

         implicit none
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         integer, intent(in) :: talsh_buff   ! buffer size for talsh in MB
         integer, intent(in) :: print_level  
         integer, intent(in) :: rcw          ! real or complex 
         real(8), intent(in) :: th_cholesky

         integer(INTD)       :: ierr
         integer             :: nao, my_MPI_rank
         type(talsh_tens_t)  :: aoint_tensor
         type(talsh_tens_t)  :: moint_tensor
         integer(C_INT)      :: aoint_dims(1:4),moint_dims(1:4)
         integer(C_SIZE_T)   :: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_INT)      :: host_arg_max
         integer             :: m, i
         complex(8), pointer, contiguous :: ao_tens(:,:,:,:),mo_tens(:,:,:,:)
         type(C_PTR)                     :: body_p

         type(talsh_tens_t), allocatable :: chol_tensor(:), temp(:)
         type(talsh_tens_t), allocatable :: upq(:), urs(:)
         integer(C_INT)                  :: dim2(1:2)

         type(talsh_tens_t)              :: chol_tens, temp_tens
         type(talsh_tens_t)              :: upq_tens, urs_tens
         integer(C_INT)                  :: dim3(1:3), dim_chol(1:3)

        ! the additional variables needed for the student projects
         type(one_el_t) :: one_int
         logical        :: one_el_exist
         real(8)        :: e_core

!        Make some noise so that we know we are here
         call print_date('Entered compute_mo_integrals routine')

!        Initialize global data
         call initialize_global_data
#if defined (VAR_MPI)
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
#else
         my_MPI_rank = 0
#endif
         if (my_MPI_rank /= 0) return ! TALSH is not parallelized, leave the work to master
         call print_date('Initialized global data')

!        Initialize libraries
         buf_size=talsh_buff*buf_size
         ierr=talsh_init(buf_size,host_arg_max)
         call print_date('Initialized talsh library')
         write(*,'(" Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max
         call interest_initialize(.true.)
         call print_date('Initialized InteRest library')

!        Retrieve basis set information
         nao = get_nao()

!        Set tensor dimensions for the index transformation (simple brute-force algorithm)
         aoint_dims = nao
         moint_dims = nmo

         if (th_cholesky.le.0) then
!          Get tensor for AO integrals
           ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=ZERO)
!          Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
           ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
           call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
           call compute_ao_integrals (nao,ao_tens_complex=ao_tens)
           call print_date('Calculated AO integrals')

!          Transform to MO basis
           ierr=talsh_tensor_construct(moint_tensor,C8,moint_dims,init_val=ZERO)
           call aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)
           if (print_level.gt.2) then
             print*, "AOint = ",tensor_norm2(aoint_tensor)
             print*, "MOint = ",tensor_norm2(moint_tensor)
             print*, "additional transform to get anti-symmetrized integrals"
           end if

         else if (.FALSE.) then
            call print_date(' Cholesky method 41')
            !determine cholesky vectors for ao integrals
            call decompose_cholesky_talsh_vec (chol_tensor, m, th_cholesky, nao, print_level)
            
            allocate(temp(m))
            allocate(upq(m))
            allocate(urs(m))

            ! transform first part
            dim2=nao
            dim2(1)=nmo(1)
            do i=1,m
              ierr=talsh_tensor_construct(temp(i),C8,dim2,init_val=ZERO)
              ierr=talsh_tensor_construct(upq(i),C8,nmo(1:2),init_val=ZERO)
            end do

            call ao2mo_vec(chol_tensor,temp,m,nmo(1),mo_list(1:nmo(1)),11)
            call ao2mo_vec(temp,upq,m,nmo(2),mo_list(nmo(1)+1:sum(nmo(1:2))),21)
            do i=1,m
              ierr=talsh_tensor_init(temp(i))
            end do
            call ao2mo_vec(chol_tensor,temp,m,nmo(1),mo_list(1:nmo(1)),12)
            call ao2mo_vec(temp,upq,m,nmo(2),mo_list(nmo(1)+1:sum(nmo(1:2))),22)

            call print_date(' upq done ')

            ! transform second part
            dim2=nao
            dim2(1)=nmo(3)
            do i=1,m
              ierr=talsh_tensor_destruct(temp(i))
              ierr=talsh_tensor_construct(urs(i),C8,nmo(3:4),init_val=ZERO)
              ierr=talsh_tensor_construct(temp(i),C8,dim2,init_val=ZERO)
            end do

            call ao2mo_vec(chol_tensor,temp,m,nmo(3),mo_list(sum(nmo(1:2))+1:sum(nmo(1:3))),11)
            call ao2mo_vec(temp,urs,m,nmo(4),mo_list(sum(nmo(1:3))+1:sum(nmo(1:4))),21)
            do i=1,m
              ierr=talsh_tensor_init(temp(i))
            end do
            call ao2mo_vec(chol_tensor,temp,m,nmo(3),mo_list(sum(nmo(1:2))+1:sum(nmo(1:3))),12)
            call ao2mo_vec(temp,urs,m,nmo(4),mo_list(sum(nmo(1:3))+1:sum(nmo(1:4))),22)

            ! delete cholesky vector
            do i=1,m
              ierr=talsh_tensor_destruct(chol_tensor(i))
              ierr=talsh_tensor_destruct(temp(i))
            end do
            deallocate(chol_tensor)
            deallocate(temp)

            call print_date(' urs done ')

            ! get the tensor
            ierr=talsh_tensor_construct(moint_tensor,C8,moint_dims,init_val=ZERO)
            do i=1,m
              ierr=talsh_tensor_contract("T(p,q,r,s)+=K(p,q)*L(r,s)",moint_tensor,upq(i),urs(i))
            end do
            if (ierr.ne.0) stop 'program error: cholesky not working'
            if (print_level.gt.2) write (*,*) "MOint = ",tensor_norm2(moint_tensor)

            ! Clean-up
            do i=1,m
              ierr=talsh_tensor_destruct(upq(i))
              ierr=talsh_tensor_destruct(urs(i))
            end do 
            deallocate(upq)
            deallocate(urs)

         else
            call print_date('Cholesky method 42')
            !determine cholesky vectors for ao integrals
            call decompose_cholesky_talsh_all (chol_tens, th_cholesky, nao, print_level)
            ierr = talsh_tensor_dimensions(chol_tens,i,dim_chol)

            !  transform upq
            dim3=dim_chol
            dim3(1)=nmo(1)
            ierr=talsh_tensor_construct(temp_tens,C8,dim3,init_val=ZERO)
            dim3(1:2)=nmo(1:2)
            ierr=talsh_tensor_construct(upq_tens,C8,dim3,init_val=ZERO)

            call ao2mo_ind(chol_tens,temp_tens,nmo(1),mo_list(1:nmo(1)),11)
            call ao2mo_ind(temp_tens,upq_tens,nmo(2),mo_list(nmo(1)+1:sum(nmo(1:2))),21)
            ierr=talsh_tensor_init(temp_tens)
            call ao2mo_ind(chol_tens,temp_tens,nmo(1),mo_list(1:nmo(1)),12)
            call ao2mo_ind(temp_tens,upq_tens,nmo(2),mo_list(nmo(1)+1:sum(nmo(1:2))),22)

            ierr=talsh_tensor_destruct(temp_tens)

            if (print_level.gt.2) call print_date(' upq done ')

            !  transform urs
            dim3=dim_chol
            dim3(1)=nmo(3)
            ierr=talsh_tensor_construct(temp_tens,C8,dim3,init_val=ZERO)
            dim3(1:2)=nmo(3:4)
            ierr=talsh_tensor_construct(urs_tens,C8,dim3,init_val=ZERO)

            call ao2mo_ind(chol_tens,temp_tens,nmo(3),mo_list(sum(nmo(1:2))+1:sum(nmo(1:3))),11)
            call ao2mo_ind(temp_tens,urs_tens,nmo(4),mo_list(sum(nmo(1:3))+1:sum(nmo)),21)
            ierr=talsh_tensor_init(temp_tens)
            call ao2mo_ind(chol_tens,temp_tens,nmo(3),mo_list(sum(nmo(1:2))+1:sum(nmo(1:3))),12)
            call ao2mo_ind(temp_tens,urs_tens,nmo(4),mo_list(sum(nmo(1:3))+1:sum(nmo)),22)

            ierr=talsh_tensor_destruct(temp_tens)
            ierr=talsh_tensor_destruct(chol_tens)

            if (print_level.gt.2) call print_date(' urs done ')

            ierr=talsh_tensor_construct(moint_tensor,C8,moint_dims,init_val=ZERO)
            ierr=talsh_tensor_contract("T(p,q,r,s)+=K(p,q,u)*L(r,s,u)",moint_tensor,upq_tens,urs_tens)
            if (ierr.ne.0) stop 'program error: cholesky not working'

            ierr=talsh_tensor_destruct(upq_tens)
            ierr=talsh_tensor_destruct(urs_tens)

         end if
         call print_date(' 2-el Integrals done ')

!        Write to old format MO-integrals file to allow for easy debugging of the integral transform
         ierr=talsh_tensor_get_body_access(moint_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,mo_tens,moint_dims(1:4))
         call write_to_mdcint (nmo,mo_tens,rcw)

!        get one electron integrals
         one_el_exist = exist_one_el(e_core)

!        call the two possible drivers, only one needs to be made working
         if (one_el_exist) then
           print*,"all integrals available: we proceed to MO-Hartree-Fock"
           call read_from_mrconee (one_int, ierr)
           call hartree_fock_driver(one_int,mo_tens)
         else
           print*,"no one-electron integrals available"
         end if

         if (print_level.gt.8) then
            print*, "additional transform to get anti-symmetrized integralss"
            call get_integral_tensor (aoint_tensor,moint_tensor,12,nmo,mo_list)
            print*, "MOint-antisym = ",tensor_norm2(moint_tensor)
         end if

!        Destruct remaining tensors and shut down library
         ierr=talsh_tensor_destruct(aoint_tensor)
         ierr=talsh_tensor_destruct(moint_tensor)
!        Get statistics
         ierr=talsh_stats()
         ierr = talsh_shutdown()

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         call print_date('Leaving compute_mo_integrals routine')

        end subroutine talsh_compute_mo_integrals

        subroutine get_mocoef_tensor (nao, nmo, mo_list, mocoefa_tensor, mocoefb_tensor)

!This subroutine reads the desired subset of mo coefficients from DIRAC

!        Written by Lucas Visscher, winter 2016/2017 (but in Goa, India, temperature about 30 C)

         use exacorr_mo
         use exacorr_global

         implicit none
         integer, intent(in) :: nao          ! the length of the ao basis
         integer, intent(in) :: nmo          ! the length of the mo basis
         integer, intent(in) :: mo_list(:)   ! and their indices
         type(talsh_tens_t), intent(inout) :: mocoefa_tensor,mocoefb_tensor

         integer(INTD)       :: ierr, tens_rank
         type(cmo)           :: cspinor
         integer(INTD)       :: mocoef_dims(1:2)

         complex(8), pointer, contiguous:: camo_tens(:,:),cbmo_tens(:,:)
         type(C_PTR):: body_p

         ierr=talsh_tensor_dimensions(mocoefa_tensor,tens_rank,mocoef_dims)
         ierr=talsh_tensor_get_body_access(mocoefa_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,camo_tens,mocoef_dims(1:2)) ! to use camo_tens as a regular Fortran 2d array
         ierr=talsh_tensor_get_body_access(mocoefb_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,cbmo_tens,mocoef_dims(1:2)) ! to use cbmo_tens as a regular Fortran 2d array

         call get_mo_coefficients (cspinor,mo_list,nmo)

!        Copy the coefficients into the alpha and beta tensors, using the quaternion relations as defined in
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
!        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
!        using Kramers' symmetry.

         camo_tens = cspinor%coeff(:,:,1)
         cbmo_tens = cspinor%coeff(:,:,2)

        end subroutine get_mocoef_tensor

        subroutine compute_ao_integrals (nao,ao_tens_complex,ao_tens_real)

!This subroutine drives the generation of AO integrals by InteRest
!
         use exacorr_datatypes
         use exacorr_global
         use module_interest_eri

         implicit none
         integer, intent(in) :: nao           ! number of basis functions
         integer             :: basis_angular ! 1=cartesian, 2=spherical
         complex(8),intent(out),optional :: ao_tens_complex(nao,nao,nao,nao)
         real(8),intent(out),optional    :: ao_tens_real(nao,nao,nao,nao)

         integer :: i, j, k, l, ish, jsh, ksh, lsh, ioff, joff, koff, loff, ijkl, ijkl2
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8) :: gout(21*21*21*21)
         real(8), parameter :: fijkl = 1.0
         integer            :: nshells                    ! number of AO shells
         type(basis_func_info_t), allocatable :: gto(:)   ! array with basis function information
         logical complex_output

         if (present(ao_tens_complex)) then
            complex_output = .true.
         else if (present(ao_tens_real)) then
            complex_output = .false.
         else
            call quit ('no output array given in compute_ao_integrals')
         end if

         basis_angular = get_basis_angular()
         call get_gtos(1,nao,gto,nshells)

         !OpenMP task-based parallelism defined by Eric does not work, deactivated [LV]
         !deactivatedopenmp parallel
         !deactivatedopenmp single
         loff = 0
         do lsh = 1, nshells
            ll   =  gto(lsh)%orb_momentum
            el   =  gto(lsh)%exponent(1)
            xl   =  gto(lsh)%coord(1)
            yl   =  gto(lsh)%coord(2)
            zl   =  gto(lsh)%coord(3)
            cl   =  gto(lsh)%coefficient(1)
            nl   =  nfunctions(ll,basis_angular)
            koff = 0
            do ksh = 1, nshells
               lk   =  gto(ksh)%orb_momentum
               ek   =  gto(ksh)%exponent(1)
               xk   =  gto(ksh)%coord(1)
               yk   =  gto(ksh)%coord(2)
               zk   =  gto(ksh)%coord(3)
               ck   =  gto(ksh)%coefficient(1)
               nk   =  nfunctions(lk,basis_angular)
               joff = 0
               do jsh = 1, nshells
                  lj   =  gto(jsh)%orb_momentum
                  ej   =  gto(jsh)%exponent(1)
                  xj   =  gto(jsh)%coord(1)
                  yj   =  gto(jsh)%coord(2)
                  zj   =  gto(jsh)%coord(3)
                  cj   =  gto(jsh)%coefficient(1)
                  nj   =  nfunctions(lj,basis_angular)
                  ioff = 0
                  do ish = 1, nshells
                     li   =  gto(ish)%orb_momentum
                     ei   =  gto(ish)%exponent(1)
                     xi   =  gto(ish)%coord(1)
                     yi   =  gto(ish)%coord(2)
                     zi   =  gto(ish)%coord(3)
                     ci   =  gto(ish)%coefficient(1)
                     ni   =  nfunctions(li,basis_angular)
!                    call interest and get the (LL|LL) type integrals
!                    output order of eri is (c,d,a,b), so input k,l first and then i,j to get the order that we want

                     ![eric]
                     !deactivatedopenmp task
                     call interest_eri('llll',fijkl,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )

                     ijkl = 0
                     ijkl2 = 1
                     do l = 1, nl
                        do k = 1, nk
                           do j = 1, nj
                              do i = 1, ni
                                 ijkl = ijkl + 1
                                 if (complex_output) then
                                    ao_tens_complex(i+ioff,j+joff,k+koff,l+loff) = dcmplx(gout(ijkl),0.D0)
                                 else
                                    ao_tens_real(i+ioff,j+joff,k+koff,l+loff) = gout(ijkl)
                                 end if
!lv                              ao_tens(i+ioff,j+joff,k+koff,l+loff) = dcmplx(gout2(ijkl2),0.D0)
!lv                              ijkl2 = ijkl2 + 5
                              end do
                           end do
!lv                        ijkl2 = ijkl2 + ni*nj*5*4
                        end do
                     end do
                     !deactivatedopenmp end task

                     ioff = ioff + ni
                  end do
                  joff = joff + nj
               end do
               koff = koff + nk
            end do
            loff = loff + nl
         end do
         !deactivatedopenmp end single nowait
         !deactivatedopenmp end parallel

         deallocate(gto)

        end subroutine compute_ao_integrals

        subroutine aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)

!        Transform integral tensor from ao-basis to mo-basis

         implicit none

         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         type(talsh_tens_t), intent(inout)  :: aoint_tensor, moint_tensor

         type(talsh_tens_t):: mocoefa_tensor(4),mocoefb_tensor(4)
         type(talsh_tens_t):: a1int_tensor,a2int_tensor, a3int_tensor
         integer(INTD) :: aoint_dims(1:4),mocoef_dims(1:2),a1int_dims(1:4)
         integer(INTD) :: a2int_dims(1:4),a3int_dims(1:4),moint_dims(1:4)

         integer(INTD) :: ierr, tens_rank
         integer       :: i, imo, nao

         ! call print_date('entered aomo_transform')
!        Get tensor dimension of input and output tensors
         ierr = talsh_tensor_dimensions(aoint_tensor,tens_rank,aoint_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: aoint tensor corrupted'
         ierr = talsh_tensor_dimensions(moint_tensor,tens_rank,moint_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: moint tensor corrupted'

!        Initialize tensors with MO coefficients
         nao = aoint_dims(1)
         imo = 0
         do i = 1, 4
            mocoef_dims(1) = aoint_dims(i)
            if (mocoef_dims(1) .ne. nao)      stop 'program error: inconsistent dimensions'
            mocoef_dims(2) = moint_dims(i)
            if (mocoef_dims(2) .ne. nmo(i)) stop 'program error: inconsistent dimensions'
            ierr=talsh_tensor_construct(mocoefa_tensor(i),C8,mocoef_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(mocoefb_tensor(i),C8,mocoef_dims,init_val=ZERO)
            call get_mocoef_tensor (nao, nmo(i), mo_list(imo+1:imo+nmo(i)), mocoefa_tensor(i), mocoefb_tensor(i))
            imo = imo + nmo(i)
         end do
         ! call print_date('initialized mo-coefficients')


!        Set tensor dimensions for the index transformation (simple brute-force algorithm)
         a1int_dims(1)   = moint_dims(1)
         a1int_dims(2:4) = aoint_dims(2:4)
         a2int_dims(1:2) = moint_dims(1:2)
         a2int_dims(3:4) = aoint_dims(3:4)
         a3int_dims(1:3) = moint_dims(1:3)
         a3int_dims(4)   = aoint_dims(4)

!        Do 1st and 2nd quarter-transformation
         ierr=talsh_tensor_construct(a1int_tensor,C8,a1int_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(a2int_tensor,C8,a2int_dims,init_val=ZERO)
         ! call print_date('constructed a1 and a2')
         ierr=talsh_tensor_contract("D(i,q,r,s)+=R+(p,i)*L(p,q,r,s)",a1int_tensor,mocoefa_tensor(1),aoint_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=11; return; endif
         ierr=talsh_tensor_contract("D(i,j,r,s)+=L(i,q,r,s)*R(q,j)",a2int_tensor,a1int_tensor,mocoefa_tensor(2))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=21; return; endif
         ierr=talsh_tensor_init(a1int_tensor,ZERO)
         ierr=talsh_tensor_contract("D(i,q,r,s)=R+(p,i)*L(p,q,r,s)",a1int_tensor,mocoefb_tensor(1),aoint_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=12; return; endif
         ierr=talsh_tensor_contract("D(i,j,r,s)+=L(i,q,r,s)*R(q,j)",a2int_tensor,a1int_tensor,mocoefb_tensor(2))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=22; return; endif
         ! call print_date('first ht done')

!        Destruct a1 and a1 integrals:
         ierr=talsh_tensor_destruct(a1int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=25; return; endif

!        Do 3rd and 4th quarter-transformation
         ierr=talsh_tensor_construct(a3int_tensor,C8,a3int_dims,init_val=ZERO)
         ! call print_date('constructed a3')
         ierr=talsh_tensor_contract("D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)",a3int_tensor,mocoefa_tensor(3),a2int_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=31; return; endif
         ierr=talsh_tensor_contract("D(i,j,k,l)+=L(i,j,k,s)*R(s,l)",moint_tensor,a3int_tensor,mocoefa_tensor(4))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=41; return; endif
         ierr=talsh_tensor_init(a3int_tensor,ZERO)
         ierr=talsh_tensor_contract("D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)",a3int_tensor,mocoefb_tensor(3),a2int_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=32; return; endif
         ierr=talsh_tensor_contract("D(i,j,k,l)+=L(i,j,k,s)*R(s,l)",moint_tensor,a3int_tensor,mocoefb_tensor(4))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=42; return; endif
         ! call print_date('second ht done')

!        Destruct a2 and a3 integrals:
         ierr=talsh_tensor_destruct(a2int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=35; return; endif
         ierr=talsh_tensor_destruct(a3int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=45; return; endif
         ! call print_date('leaving aomo_transform')

        end subroutine aomo_transform

        subroutine get_integral_tensor (aoint_tensor,integral_tensor,antisymmetrize,nmo,mo_list)

         implicit none

         type(talsh_tens_t), intent(inout)  :: aoint_tensor,integral_tensor
         integer, intent(in)                :: nmo(4), mo_list(:)
         integer, intent(in)                :: antisymmetrize

         integer(INTD)                      :: ierr, moint_dims(4)
         type(talsh_tens_t)                 :: moint_tensor

         moint_dims = nmo
         ierr=talsh_tensor_construct(moint_tensor,C8,moint_dims,init_val=ZERO)
         if (ierr.ne.0) stop 'get_integral_tensor: tensor construction failed'
!        Transform integrals
         call aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)
!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)
         ierr=talsh_tensor_destruct(moint_tensor)

        end subroutine get_integral_tensor

        subroutine mulliken_to_dirac_sort  (mulliken_tensor,dirac_tensor,antisymmetrize)

         implicit none

         type(talsh_tens_t), intent(inout)  :: mulliken_tensor,dirac_tensor
         integer, intent(in)                :: antisymmetrize

         integer(INTD) :: mulliken_dims(1:4),dirac_dims(1:4)

         integer(INTD) :: ierr, tens_rank, mismatch

!        work-around: addition with permutation is not implemented
         type(talsh_tens_t) :: one_tensor, minusone_tensor
         integer(C_INT)     :: scalar_dims(1)
         scalar_dims(1) = 1
         ierr=talsh_tensor_construct(one_tensor,C8,scalar_dims(1:0),init_val=ONE)
         ierr=talsh_tensor_construct(minusone_tensor,C8,scalar_dims(1:0),init_val=MINUS_ONE)
!        end work-round

!        get dimensions of the tensors and check consistency
         ierr = talsh_tensor_dimensions(mulliken_tensor,tens_rank,mulliken_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input tensor corrupted in mulliken_to_dirac_sort'
         ierr = talsh_tensor_dimensions(dirac_tensor,tens_rank,dirac_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'output tensor corrupted in mulliken_to_dirac_sort'
         mismatch = 0
         if (mulliken_dims(1) .ne. dirac_dims(1)) mismatch = mismatch + 1
         if (mulliken_dims(2) .ne. dirac_dims(3)) mismatch = mismatch + 1
         if (mulliken_dims(3) .ne. dirac_dims(2)) mismatch = mismatch + 1
         if (mulliken_dims(4) .ne. dirac_dims(4)) mismatch = mismatch + 1
         if (antisymmetrize .eq. 12 .and. dirac_dims(1) .ne. dirac_dims(2)) mismatch = mismatch + 1
         if (antisymmetrize .eq. 34 .and. dirac_dims(3) .ne. dirac_dims(4)) mismatch = mismatch + 1
         if (mismatch.gt.0) stop 'dimension mismatch in mulliken_to_dirac_sort'

         ierr=talsh_tensor_init(dirac_tensor)
!shorter ierr=talsh_tensor_add("D(p,q,r,s)+=M(p,r,q,s)",dirac_tensor,mulliken_tensor)
         ierr=talsh_tensor_contract("D(p,q,r,s)+=M(p,r,q,s)*C()",dirac_tensor,mulliken_tensor,one_tensor)

         if (antisymmetrize .eq. 12) then
!shorter     ierr=talsh_tensor_add("D(p,q,r,s)+=M(q,r,p,s)",dirac_tensor,mulliken_tensor,scale=MINUS_ONE)
             ierr=talsh_tensor_contract("D(p,q,r,s)+=M(q,r,p,s)*C()",dirac_tensor,mulliken_tensor,minusone_tensor)
         elseif (antisymmetrize .eq. 34) then
!shorter     ierr=talsh_tensor_add("D(p,q,r,s)+=M(p,s,q,r)",dirac_tensor,mulliken_tensor,scale=MINUS_ONE)
              ierr=talsh_tensor_contract("D(p,q,r,s)+=M(p,s,q,r)*C()",dirac_tensor,mulliken_tensor,minusone_tensor)
         end if

!        work-around: addition with permutation is not implemented
         ierr=talsh_tensor_destruct(one_tensor)
         ierr=talsh_tensor_destruct(minusone_tensor)
!        end work-round

        end subroutine mulliken_to_dirac_sort

        function tensor_norm2(a_tensor) result(a_norm2)

         !compute norm2(complete contraction of a tensor with itself) for debugging purposed

         type(talsh_tens_t), intent(inout):: a_tensor

         integer(INTD) :: a_dims(1:4)
         integer(INTD) :: ierr, a_rank

!        scalars (need to be defined as tensor types)
         type(talsh_tens_t)  :: result_tensor
         integer(C_INT)      :: result_dims(1)
         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         real(8) :: a_norm2

!        Initialize scalars that are to be used as tensors in contractions
         result_dims(1) = 1
         ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,result_dims)

         ierr = talsh_tensor_dimensions(a_tensor,a_rank,a_dims)
         select case (a_rank)
         case (0)
            ierr=talsh_tensor_contract("R()+=A+()*A()",result_tensor,a_tensor,a_tensor)
         case (1)
            ierr=talsh_tensor_contract("R()+=A+(p)*A(p)",result_tensor,a_tensor,a_tensor)
         case (2)
            ierr=talsh_tensor_contract("R()+=A+(p,q)*A(p,q)",result_tensor,a_tensor,a_tensor)
         case (3)
            ierr=talsh_tensor_contract("R()+=A+(p,q,r)*A(p,q,r)",result_tensor,a_tensor,a_tensor)
         case (4)
            ierr=talsh_tensor_contract("R()+=A+(p,q,r,s)*A(p,q,r,s)",result_tensor,a_tensor,a_tensor)
         case default
            call quit ('wrong dimension in tensor_norm2')
         end select

         a_norm2 = real(result_tens(1),8)
         ierr=talsh_tensor_destruct(result_tensor)

        end function tensor_norm2

        subroutine ao2mo_ind(ao_tensor,mo_tensor,nmo,mo_list,ind)
!         Transform integral tensor from ao-basis to mo-basis

          use talsh_common_routines

          implicit none

          type(talsh_tens_t), intent(inout) :: ao_tensor
          type(talsh_tens_t), intent(inout) :: mo_tensor
          integer(INTD), intent(in) :: nmo          ! the size of the mo basis for each of the 2 indices (mulliken ordering)
          integer, intent(in)       :: mo_list(:)   ! and their indices (2 lists consecutively)
          integer, intent(in)       :: ind

          type(talsh_tens_t):: mocoef_alpha, mocoef_beta
          integer(INTD)     :: mocoef_dims(1:2)
          integer(INTD)     :: ao_dim2(1:2), mo_dim2(1:2)
          integer(INTD)     :: ao_dim3(1:3), mo_dim3(1:3)
          integer(INTD)     :: ao_dim4(1:4), mo_dim4(1:4)
          integer(INTD)     :: ierr, rank, i

!         Get rank of ao tensor
          rank=talsh_tensor_rank(ao_tensor)

!         for one electron integrals / matrices
          if (rank.eq.2) then
            ierr = talsh_tensor_dimensions(ao_tensor,rank,ao_dim2)
            if (ierr.ne.0 .or. rank.ne.2) stop 'error in ao2mo_ind: ao tensor corrupted'
!           check output tensor
            ierr = talsh_tensor_dimensions(mo_tensor,rank,mo_dim2)
            if (ierr.ne.0 .or. rank.ne.2) stop 'error in ao2mo_ind: mo tensor corrupted'
            if (mo_dim2(ind).ne.nmo) stop 'error in ao2mo_ind: mo tensor wrong nmo1'

!           Initialize tensors with MO coefficients
            mocoef_dims(1) = ao_dim2(ind)
            mocoef_dims(2) = mo_dim2(ind)
            ierr=talsh_tensor_construct(mocoef_alpha,C8,mocoef_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(mocoef_beta,C8,mocoef_dims,init_val=ZERO)
            call get_mocoef_tensor (mocoef_dims(1),mocoef_dims(2), & 
                        mo_list(1:mocoef_dims(2)), mocoef_alpha, mocoef_beta)

!           transform the index
            if (ind.eq.1) then
              ierr=talsh_tensor_contract("B(p,q)+=A(i,q)*M+(i,p)",mo_tensor,ao_tensor,mocoef_alpha)
              if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_ind: 2a1'
              ierr=talsh_tensor_contract("B(p,q)+=A(i,q)*M+(i,p)",mo_tensor,ao_tensor,mocoef_beta)
              if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_ind: 2b1'
            else if (ind.eq.2) then
                ierr=talsh_tensor_contract("C(p,q)+=B(p,i)*M(i,q)",mo_tensor,ao_tensor,mocoef_alpha)
                if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_ind: 2a2'
                ierr=talsh_tensor_contract("C(p,q)+=B(p,i)*M(i,q)",mo_tensor,ao_tensor,mocoef_beta)
                if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_ind: 2b2'
            else
              stop 'error in ao2mo_ind: index not available'
            end if 

!           clean up
            ierr=talsh_tensor_destruct(mocoef_alpha)
            ierr=talsh_tensor_destruct(mocoef_beta)
            
!         cholesky vectors
          else if (rank.eq.3) then
            ierr = talsh_tensor_dimensions(ao_tensor,rank,ao_dim3)
            if (ierr.ne.0 .or. rank.ne.3) stop 'error in ao2mo_ind: ao tensor corrupted'
!           check output tensor
            ierr = talsh_tensor_dimensions(mo_tensor,rank,mo_dim3)
            if (ierr.ne.0 .or. rank.ne.3) stop 'error in ao2mo_ind: mo tensor corrupted'
            if (mo_dim3(3).ne.ao_dim3(3)) stop 'error in ao2mo_ind: mo tensor wrong chol. dim.'
            if (ind.eq.11 .or. ind.eq.12) then
                i=1
            else if (ind.eq.21 .or. ind.eq.22) then
                i=2
            else
              stop 'error in ao2mo_ind: index not available'
            end if
            if (mo_dim3(i).ne.nmo) stop 'error in ao2mo_ind: mo tensor wrong nmo'

!           Initialize tensors with MO coefficients
            mocoef_dims(1) = ao_dim3(i)
            mocoef_dims(2) = mo_dim3(i)
            ierr=talsh_tensor_construct(mocoef_alpha,C8,mocoef_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(mocoef_beta,C8,mocoef_dims,init_val=ZERO)
            call get_mocoef_tensor (mocoef_dims(1),mocoef_dims(2), & 
                       mo_list(1:mocoef_dims(2)), mocoef_alpha, mocoef_beta)

!           transform index
            if (ind.eq.11) then
                ierr=talsh_tensor_contract("A(p,s,u)+=C(r,s,u)*M+(r,p)",mo_tensor,ao_tensor,mocoef_alpha)
                if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_ind: 3a1'
            else if (ind.eq.12) then
                ierr=talsh_tensor_contract("A(p,s,u)+=C(r,s,u)*M+(r,p)",mo_tensor,ao_tensor,mocoef_beta)
                if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_ind: 3b1'
            else if (ind.eq.21) then
                ierr=talsh_tensor_contract("D(p,q,u)+=A(p,s,u)*M(s,q)",mo_tensor,ao_tensor,mocoef_alpha)
                if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_ind: 3a2'
            else if (ind.eq.22) then
                ierr=talsh_tensor_contract("D(p,q,u)+=A(p,s,u)*M(s,q)",mo_tensor,ao_tensor,mocoef_beta)
                if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_ind: 3b2'
            end if

!           clean up
            ierr=talsh_tensor_destruct(mocoef_alpha)
            ierr=talsh_tensor_destruct(mocoef_beta)

          else 
            stop 'error in ao2mo_ind: wrong rank'
          end if

        end subroutine ao2mo_ind

        subroutine ao2mo_vec(ao_vec,mo_vec,m,nmo,mo_list,ind)
!         Transform integral tensor vector from ao-basis to mo-basis

          use talsh_common_routines

          implicit none
         
          type(talsh_tens_t), allocatable :: ao_vec(:)
          type(talsh_tens_t), allocatable :: mo_vec(:)
          integer(INTD), intent(in)       :: nmo          ! the size of the mo basis for each of the 2 indices (mulliken ordering)
          integer, intent(in)             :: mo_list(:)   ! and their indices (2 lists consecutively)
          integer, intent(in)             :: ind, m

          type(talsh_tens_t) :: mocoef_alpha, mocoef_beta
          integer(INTD)      :: mocoef_dims(1:2)
          integer(INTD)      :: ao_dim2(1:2), mo_dim2(1:2)
          integer(INTD)      :: ierr, rank
          integer            :: i,ii

!         Get info about ao tensor
          ierr = talsh_tensor_dimensions(ao_vec(1),rank,ao_dim2)
          if (ierr.ne.0 .or. rank.ne.2) stop 'error in ao2mo_vec: ao tensor corrupted'

          if (ind.eq.11 .or. ind.eq.12) then
            ii=1
          else if (ind.eq.21 .or. ind.eq.22) then
            ii=2
          else
            stop 'error in ao2mo_ind: index not available'
          end if
          
          mo_dim2=ao_dim2
          mo_dim2(ii)=nmo

!         Initialize tensors with MO coefficients
          mocoef_dims(1) = ao_dim2(ii)
          mocoef_dims(2) = mo_dim2(ii)
          ierr=talsh_tensor_construct(mocoef_alpha,C8,mocoef_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(mocoef_beta,C8,mocoef_dims,init_val=ZERO)
          call get_mocoef_tensor (mocoef_dims(1),mocoef_dims(2), & 
                        mo_list(1:mocoef_dims(2)), mocoef_alpha, mocoef_beta)

          do i=1,m
!           transform the index
            if (ind.eq.11) then
              ierr=talsh_tensor_contract("B(p,q)+=A(i,q)*M+(i,p)",mo_vec(i),ao_vec(i),mocoef_alpha)
              if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_vec: 2a1'
            else if (ind.eq.12) then
              ierr=talsh_tensor_contract("B(p,q)+=A(i,q)*M+(i,p)",mo_vec(i),ao_vec(i),mocoef_beta)
              if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_vec: 2b1'
            else if (ind.eq.21) then
                ierr=talsh_tensor_contract("C(p,q)+=B(p,i)*M(i,q)",mo_vec(i),ao_vec(i),mocoef_alpha)
                if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_vec: 2a2'
            else if (ind.eq.22) then
                ierr=talsh_tensor_contract("C(p,q)+=B(p,i)*M(i,q)",mo_vec(i),ao_vec(i),mocoef_beta)
                if(ierr.ne.TALSH_SUCCESS) stop 'error in ao2mo_vec: 2b2'
            end if 
          end do

!         clean up
          ierr=talsh_tensor_destruct(mocoef_alpha)
          ierr=talsh_tensor_destruct(mocoef_beta)

        end subroutine ao2mo_vec

    subroutine decompose_cholesky_talsh_vec (CD, m, t_cholesky, nao, print_level)

!       Cholesky decomposition based on treating the shells seperately using talsh tensor operations
!       possible additions: 
!                          SVD (How to combine with tensor library)
!                          Additional screening in the computation of integrals

!       written by Johann Pototschnig, Summer 2019

        use exacorr_global
        use exacorr_cholesky
        use talsh_common_routines

        implicit none

        real(8), intent(in)                            :: t_cholesky ! cholesky threshold
        type(talsh_tens_t), allocatable, intent(inout) :: CD(:)
        integer, intent(in)                            :: nao, print_level
        integer(C_INT), intent(out)                    :: m

        integer(C_INT)     :: ierr
        integer            :: lshell, kshell, lsize, ksize
        integer            :: l, k, l0, k0, i
        complex(8)         :: D_max, D_el, factor
        type(talsh_tens_t) :: Mask, Diag, Res, Temp
        type(talsh_tens_t) :: Dlk !, CD_v
        integer(C_INT)     :: dims2(1:2),dims4(1:4), hdims(1:2)
        integer(C_INT)     :: off2(1:2),off4(1:4)
        integer(C_INT)     :: rank
        type(talsh_tens_t) :: one_tensor
        integer(C_INT)     :: one_dims(1)

        ! accuracy for tensor printing
        real(8), parameter :: ACC=1.0E-9

        one_dims(1) = 1
        ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)

        if (print_level.gt.2) write(*,*) "========= entered decompose_cholesky_talsh ========= "

        hdims=1
        dims2=nao
        dims4=nao
        off2=0
        off4=0

        !Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
        call interest_initialize(.false.)

        !alllocate cholesky tensor
        allocate(CD(nao*nao))

        ierr=talsh_tensor_construct(Mask,R4,dims2,init_val=ONE)
        ierr=talsh_tensor_construct(Diag,C8,dims2,init_val=ZERO)
        ierr=talsh_tensor_construct(Res,C8,dims2,init_val=ZERO)
        ierr=talsh_tensor_construct(Temp,C8,dims2,init_val=ZERO)
        !ierr=talsh_tensor_construct(CD_v,C8,hdims,init_val=ZERO)

        call talsh_get_diag(Diag)
        if (print_level.gt.2) write(*,*) "--- diagonal computed "

        call talsh_find_shells(Diag, Mask, lshell, kshell, D_max)
        if (print_level.gt.6) write(*,*) "--- shell found: ", lshell, kshell, lsize, ksize, D_max

        m=0

        do while (abs(D_max)>t_cholesky)
          lsize=get_shell_size(lshell)
          ksize=get_shell_size(kshell)
          if (print_level.gt.2) then
            write (*,'(a,i6,a,a,i6,i6,a,e12.2)') ' cholesky length =',m,' largest element:', &
            ' shell = ', lshell, kshell, ' value = ', abs(D_max)
          end if

          dims4(3)=lsize
          dims4(4)=ksize

          ierr=talsh_tensor_construct(Dlk,C8,dims4,init_val=ZERO)

          !get the values for the current shell
          call talsh_get_shell(Dlk, lshell, kshell)
          if (print_level.gt.15) call print_tensor(Dlk, ACC, "Dlk")

          call talsh_find_el(Diag, Mask, lshell, kshell, l, k, D_el)
          if (print_level.gt.6) write(*,*) '--- largest value found: ', l, k, D_el

          do while(abs(D_el)>abs(D_max/1000))
            call get_shell_offset(l,l0)
            call get_shell_offset(k,k0)

            !create current elment of the cholesky vector
            ierr=talsh_tensor_construct(CD(m+1),C8,dims2,init_val=ZERO)
            ierr=talsh_tensor_init(Res, ZERO)

            !sum up previous CD vectors for the indices
            off2(1)=l
            off2(2)=k
            do i = 1, m
                !ierr=talsh_tensor_init(CD_v, ZERO)
                !ierr=talsh_tensor_slice(CD_v,CD(i),off2,0,DEV_HOST)
                !ierr=talsh_tensor_reshape(CD_v,hdims(1:0))
                !if (print_level.gt.6) call print_tensor(CD_v, ACC, "CD_v")
                !ierr=talsh_tensor_contract("L(p,q)+=X(p,q)*Y()",Res,CD(i),CD_v)
                !ierr=talsh_tensor_reshape(CD_v,hdims)
                factor=get_element(CD(i),off2)
                if (print_level.gt.10) write(*,*) 'factor=',factor
                ierr=talsh_tensor_contract("L(p,q)+=X(p,q)*Y()",Res,CD(i),one_tensor,scale=factor)
                if (print_level.gt.15) call print_tensor(Res, ACC, "Res")
            end do

            !get Dlk for the current l and k
            ierr=talsh_tensor_init(Temp, ZERO)
            dims4(3)=1
            dims4(4)=1
            ierr=talsh_tensor_reshape(Temp,dims4)
            off4=0
            off4(3)=l-l0-1
            off4(4)=k-k0-1
            ierr=talsh_tensor_slice(Temp,Dlk,off4,0,DEV_HOST)
            ierr=talsh_tensor_reshape(Temp,dims2)

            !sum and denominate
            ierr=talsh_tensor_contract("T(p,q)+=V(p,q)",Temp,Res,one_tensor,scale=MINUS_ONE)
            factor=ONE / zsqrt(D_el)
            ierr=talsh_tensor_contract("W(p,q)+=V(p,q)",CD(m+1),Temp,one_tensor,scale=factor)
            
            !update diagonal
            ierr=talsh_tensor_init(Res, ZERO)
            ierr=talsh_tensor_contract("W(p,q)+=V(p,q)",Res,CD(m+1),one_tensor)
            if (print_level.gt.15) call print_tensor(Res, ACC, "CD entry")
            call talsh_square_el(Res)
            ierr=talsh_tensor_contract("W(p,q)+=V(p,q)",Diag,Res,one_tensor,scale=MINUS_ONE)
            if (print_level.gt.15) call print_tensor(Diag, ACC, "Diag")

            !update mask
            call talsh_set_mask_done(Mask,l,k)

            m=m+1

            if (m+1>nao*nao) then
                write(*,*) 'WARNING: Cholesky decomposition did not work'
                write(*,*) 'WARNING: Final error:',D_el
                exit
            end if

            !find next element in this shell
            call talsh_find_el(Diag, Mask, lshell, kshell, l, k, D_el)
            if (print_level.gt.6) write(*,*)  '--- largest value found: ', l, k, D_el

          end do

          !find shell with next largest element
          call talsh_find_shells(Diag, Mask, lshell, kshell, D_max)

          ierr=talsh_tensor_destruct(Dlk)

          if (m+1>nao*nao) then
                exit
          end if

        end do
        if (print_level.gt.-1) write(*,*) ' Cholesky decomposition done, used ',m,' elements'

        !clean up
        ierr=talsh_tensor_destruct(Mask)
        ierr=talsh_tensor_destruct(Diag)
        ierr=talsh_tensor_destruct(Res)
        ierr=talsh_tensor_destruct(Temp)
        ierr=talsh_tensor_destruct(one_tensor)

        if (print_level.gt.2) write(*,*)  "========= leaving decompose_cholesky_talsh ========= "

    end subroutine decompose_cholesky_talsh_vec

    subroutine decompose_cholesky_talsh_all (CD, t_cholesky, nao, print_level)

!       Cholesky decomposition based on treating the shells seperately using talsh tensor operations
!       possible additions: 
!                          SVD (How to combine with tensor library)
!                          Additional screening in the computation of integrals

!       written by Johann Pototschnig, Summer 2019

        use exacorr_global
        use exacorr_cholesky
        use talsh_common_routines

        implicit none

        real(8), intent(in)                  :: t_cholesky ! cholesky threshold
        type(talsh_tens_t), intent(inout)    :: CD
        integer, intent(in)                  :: nao,print_level

        integer(C_INT)     :: ierr
        integer            :: lshell, kshell, lsize, ksize
        integer            :: l, k, l0, k0
        complex(8)         :: D_max, D_el, denom_fac
        type(talsh_tens_t) :: Mask, Diag, Res, Temp
        type(talsh_tens_t) :: Dlk, CD_v, CD_m, CD_int
        integer(C_INT)     :: dims2(1:2)
        integer(C_INT)     :: dims3(1:3), off3(1:3)
        integer(C_INT)     :: dims4(1:4), off4(1:4)
        integer(C_INT)     :: rank
        integer(C_INT)     :: m
        type(talsh_tens_t) :: one_tensor
        integer(C_INT)     :: one_dims(1)

        ! accuracy for tensor printing
        real(8), parameter :: ACC=1.0E-15

        one_dims(1) = 1
        ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)

        if (print_level.gt.2) write(*,*) "========= entered decompose_cholesky_talsh ========= "

        dims2=nao
        dims3=nao
        dims4=nao
        off4=0
        off3=0

        !Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
        call interest_initialize(.false.)

        ierr=talsh_tensor_construct(Mask,R4,dims2,init_val=ONE)
        ierr=talsh_tensor_construct(Diag,C8,dims2,init_val=ZERO)
        ierr=talsh_tensor_construct(Res,C8,dims2,init_val=ZERO)
        ierr=talsh_tensor_construct(Temp,C8,dims2,init_val=ZERO)


        call talsh_get_diag(Diag)
        if (print_level.gt.2) write(*,*) "--- diagonal computed "

        call talsh_find_shells(Diag, Mask, lshell, kshell, D_max)
        if (print_level.gt.6) write(*,*) "--- shell found: ", lshell, kshell, lsize, ksize, D_max

        m=0
        dims3(3)=1
        ierr=talsh_tensor_construct(CD,C8,dims3,init_val=ZERO)
        ierr=talsh_tensor_construct(CD_int,C8,dims3,init_val=ZERO)

        do while (abs(D_max)>t_cholesky)
          lsize=get_shell_size(lshell)
          ksize=get_shell_size(kshell)
          if (print_level.gt.2) then 
            write (*,'(a,i6,a,a,i6,i6,a,e12.2)') 'cholesky length =',m,' largest element:', &
            ' shell = ', lshell, kshell, ' value = ', abs(D_max)
          end if

          dims4(3)=lsize
          dims4(4)=ksize

          ierr=talsh_tensor_construct(Dlk,C8,dims4,init_val=ZERO)

          ierr=talsh_tensor_destruct(CD_int)
          dims3(3)=m+lsize*ksize
          ierr=talsh_tensor_construct(CD_int,C8,dims3,init_val=ZERO)
          if (m>0) then
            ierr = talsh_tensor_dimensions(CD,rank,dims3)
            off3=0
            ierr=talsh_tensor_insert(CD_int,CD,off3)
          end if

          !get the values for the current shell
          call talsh_get_shell(Dlk, lshell, kshell)
          if (print_level.gt.15) call print_tensor(Dlk, ACC, "Dlk")

          call talsh_find_el(Diag, Mask, lshell, kshell, l, k, D_el)
          if (print_level.gt.6) write(*,*) '--- largest value found: ', l, k, D_el

          do while(abs(D_el)>abs(D_max/1000))
            call get_shell_offset(l,l0)
            call get_shell_offset(k,k0)

            !sum up previous CD vectors for the indices
            ierr=talsh_tensor_init(Res, ZERO)
            if (m>0) then
              dims3(1:2)=1
              dims3(3)=m
              ierr=talsh_tensor_construct(CD_v,C8,dims3,init_val=ZERO)
              dims3(1:2)=nao
              off3(1)=l-1
              off3(2)=k-1
              off3(3)=0
              ierr=talsh_tensor_slice(CD_v,CD_int,off3,0,DEV_HOST)
              ierr=talsh_tensor_reshape(CD_v,(/ m /))

              dims3(3)=m
              ierr=talsh_tensor_construct(CD_m,C8,dims3,init_val=ZERO)
              off3=0
              ierr=talsh_tensor_slice(CD_m,CD_int,off3,0,DEV_HOST)
              if (print_level.gt.15) call print_tensor(CD_m, ACC, "CD_m")

              ierr=talsh_tensor_contract("T(p,q)+=C(p,q,u)*D(u)",Res,CD_m,CD_v)

              ierr=talsh_tensor_destruct(CD_m)
              ierr=talsh_tensor_destruct(CD_v)
            end if

            !get Dlk for the current l and k
            ierr=talsh_tensor_init(Temp, ZERO)
            dims4(3)=1
            dims4(4)=1
            ierr=talsh_tensor_reshape(Temp,dims4)
            off4=0
            off4(3)=l-l0-1
            off4(4)=k-k0-1
            ierr=talsh_tensor_slice(Temp,Dlk,off4,0,DEV_HOST)
            ierr=talsh_tensor_reshape(Temp,dims2)

            !sum and denominate
            ierr=talsh_tensor_contract("T(p,q)+=V(p,q)",Temp,Res,one_tensor,scale=MINUS_ONE)
            ierr=talsh_tensor_init(Res, ZERO)
            denom_fac=ONE / zsqrt(D_el)
            ierr=talsh_tensor_contract("W(p,q)=V(p,q)",Res,Temp,one_tensor,scale=denom_fac)
            if (print_level.gt.15) call print_tensor(Res, ACC, "CD entry")

            !add Cholesky vector
            dims3(3)=1
            ierr=talsh_tensor_reshape(Res,dims3)
            off3(3)=m
            ierr=talsh_tensor_insert(CD_int,Res,off3)
            ierr=talsh_tensor_reshape(Res,dims2)

            !update diagonal
            call talsh_square_el(Res)
            ierr=talsh_tensor_contract("W(p,q)+=V(p,q)",Diag,Res,one_tensor,scale=MINUS_ONE)
            if (print_level.gt.15) call print_tensor(Diag, ACC, "Diag")

            !update mask
            call talsh_set_mask_done(Mask,l,k)

            m=m+1

            !update CD
            ierr=talsh_tensor_destruct(CD)
            dims3(3)=m
            ierr=talsh_tensor_construct(CD,C8,dims3,init_val=ZERO)
            off3=0
            ierr=talsh_tensor_slice(CD,CD_int,off3,0,DEV_HOST)
            if (print_level.gt.15) call print_tensor(CD, ACC, "CD")

            if (m>nao*nao) write(*,*) 'Cholesky decomposition did not work'

            !find next element in this shell
            call talsh_find_el(Diag, Mask, lshell, kshell, l, k, D_el)
            if (print_level.gt.6) write(*,*)  '--- largest value found: ', l, k, D_el

          end do

          !find shell with next largest element
          call talsh_find_shells(Diag, Mask, lshell, kshell, D_max)

          ierr=talsh_tensor_destruct(Dlk)

        end do

        if (print_level.gt.-1) write(*,*) 'Cholesky decomposition done, used ',m,' elements'

        !clean up
        ierr=talsh_tensor_destruct(CD_int)

        ierr=talsh_tensor_destruct(Mask)
        ierr=talsh_tensor_destruct(Diag)
        ierr=talsh_tensor_destruct(Res)
        ierr=talsh_tensor_destruct(Temp)
        ierr=talsh_tensor_destruct(one_tensor)

        if (print_level.gt.2) write(*,*)  "========= leaving decompose_cholesky_talsh ========= "

    end subroutine decompose_cholesky_talsh_all

end module talsh_ao_to_mo
