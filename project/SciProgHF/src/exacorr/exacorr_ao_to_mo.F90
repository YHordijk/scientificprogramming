module exacorr_ao_to_mo

!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.

        use exacorr_utils
        use, intrinsic:: ISO_C_BINDING
#if defined (VAR_MPI)
        use exatensor
#endif
        implicit none
        private
        public exacorr_compute_mo_integrals

#if defined (VAR_MPI)

        interface get_integral_tensor
            module procedure get_integral_tensor
        end interface get_integral_tensor

        public mocoef_initlabel
        public dm_initlabel
        public get_integral_tensor
        public get_ht2_integral_tensor
        public print_tensornorm2
        public ao2mo_exat
        public mulliken_to_dirac_sort
        public get_CC_integrals_chol

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

         ! It is not possible to make an array of pointers, so we need a container
         type,public ::  h_spaces_container_t
              class(h_space_t), pointer :: ptr
         end type h_spaces_container_t

       contains

        subroutine exacorr_compute_mo_integrals (nmo, mo_list,print_level)

!        This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor
!        It currently tests two different algorithms as well as the get_integral_tensor routine used by the CC code.
!        Intended final purpose is two run the most efficient algorithm only and write integrals out to a file that can be used by correlation methods.
!        See talsh_compute_mo_integrals to see how this works for the single node implementation.

!        Written by Lucas Visscher, march/september 2018

         use interface_to_mpi
         use exacorr_datatypes
         use exacorr_global
         use exacorr_tensor_methods
         use exatensor
         

         implicit none

         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         integer, intent(in) :: print_level

         integer(INTD)       :: my_role
         integer             :: my_MPI_rank    ! needed to start exatensor, MPI should already be started in threaded mode
         integer             :: nshells                 ! the number of shells (orbitals with the same center and l-value)
         integer             :: nao                     ! the number of atomic orbitals)
         integer             :: basis_angular           ! 1=cartesian, 2=spherical
         type(basis_func_info_t), allocatable :: gto(:) ! array with information about the shells

         ! classes for the definition of spaces (AO and 4 different MO spaces)
         class(h_space_t), pointer                :: ao_space
         type(h_spaces_container_t), dimension(4) :: mo_space

         ! variables for the definition of spaces (AO and 4 different MO spaces)
         type(subspace_basis_t):: basis_ao, basis_mo(4)
         integer(INTL)::          nao_l, ao_space_root
         integer(INTL)::          nmo_l, mo_space_root(4)
         integer(INTD)::          ao_space_id, mo_space_id(4)

         ! variables for the definition of colors inside the spaces
         type(color_symmetry_t), allocatable :: color(:)
         integer(INTL)::                        labs

         ! types of the initializer methods
         type(compute_2e_ao_tensor_t) :: aoint_calculator ! method to initialize the ao-integrals tensor
         type(init_mocoef_tensor_t)   :: mocoef_init(2,4) ! array of initializers for algorithm 1, first index is alpha/beta, second is 1,2,3,4
         type(init_dm_tensor_t)       :: dm_init(2)       ! array of initializers for algorithm 2, index indicates first/second halftransformation
         type(init_ht_tensor_t)       :: ht_init          ! initializers for algorithm 3 in which the HT matrix is made without ao-integral tensor
         type(tens_printer_t)         :: tens_printer     ! generic method, but needs to be defined here in specific form

         ! the two tensors that we need
         type(tens_rcrsv_t):: aoint_tensor
         type(tens_rcrsv_t):: moint_tensor

         ! loop and other auxilliary variables
         integer(INTD) ::  spin,ierr, lsh
         integer:: i,imo,jmo,imo12,jmo12,imo34,jmo34
         integer:: l, nl
         character(len=2)  :: vs_str
         integer(INTD) branch_factor

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         if (my_MPI_rank == 0) print*, " Starting AO-to-MO transform with exatensor"

!        Initialize global data (read from DFCOEF or similar by master)
         call set_talsh_only(.false.)
         call initialize_global_data

         nao           = get_nao()   ! This is the total number of ao's
         nao_l         = int(nao,8)  ! we need this as a long int for exatensor
         basis_angular = get_basis_angular() ! needed for setting basis functions

!        Get the basis set information
         call get_gtos(1,nao,gto,nshells)

         !Create a basis for the ao vector space:
         ierr = 0
         call basis_ao%subspace_basis_ctor(nao_l,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

         !Assign the colors for the ao space (this indicates to exatensor where it may break up into subspaces)
         allocate(color(nshells))
         labs = 0
         do lsh = 1, nshells !set basis functions
            nl   =  nfunctions(gto(lsh)%orb_momentum,basis_angular)
            call color(lsh)%color_symmetry_ctor(lsh) ! take one shell as a color
            do  l = 1, nl
                labs = labs + 1
                call basis_ao%set_basis_func(labs,BASIS_ABSTRACT,ierr,symm=color(lsh))
                if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.set_basis_func() failed!')
            end do
         enddo

         !Finalize the ao space basis
         call basis_ao%finalize(ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.finalize() failed!')

         !Register the ao vector space (block size hardwired to 75 for the moment):
         branch_factor = nao / 75 + 1
         if (branch_factor < 2) branch_factor = 2
         ierr=exatns_space_register('AOspace',basis_ao,ao_space_id,ao_space,branch_factor)
         if (my_MPI_rank == 0) print*," registered AO vector space with nao, nshells:",nao,nshells
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

         !Get the MO-coefficient basis for the four indices that are to be transformed
         ierr = 0
         imo = 1
         do i = 1, 4
            jmo     = imo + nmo(i) - 1 ! last MO for index i, whereas imo is the first
            !nmo_l   = 2*int(nmo(i),8)  ! NB: imo, jmo and nmo count orbitals, nmo_l counts spinors
            nmo_l   = int(nmo(i),8)  ! NB: imo, jmo and nmo count orbitals, nmo_l counts spinors
            call basis_mo(i)%subspace_basis_ctor(nmo_l,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

            !Set the functions inside this basis (without colors, not needed for MOs)
            do labs = 1, nmo_l
               call basis_mo(i)%set_basis_func(labs,BASIS_ABSTRACT,ierr)
               if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.set_basis_func() failed!')
            enddo

            !Finalize the mo space basis
            call basis_mo(i)%finalize(ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.finalize() failed!')

            !Register this mo vector space :
            branch_factor = nmo(i) / 75 + 1
            if (branch_factor < 2) branch_factor = 2
            write (vs_str,'(A,I1)') '_',i
            ierr=exatns_space_register('MOspace'//vs_str,basis_mo(i),mo_space_id(i),mo_space(i)%ptr,branch_factor)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

            imo = imo + nmo(i) ! make imo point to the start of the list for the next index
         end do

         !Define separate initializers for the alpha(A) and beta(B) MO-coefficients for each index
         imo = 1
         do i = 1, 4
            jmo     = imo + nmo(i) - 1 ! last MO for index i, whereas imo is the first
            do spin = 1, 2
               !Register this initializer as a method
               call mocoef_init(spin,i)%init_mocoef_ctor(mo_list(imo:jmo),nmo(i),spin)
               ierr=exatns_method_register(mocoef_initlabel(ao_space_id,mo_space_id(i),spin),mocoef_init(spin,i))
               if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
            end do
         end do

         !Define separate initializers for the first and second density matrix
         imo12 = 1
         jmo12 = nmo(1)+nmo(2)
         imo34 = jmo12 + 1
         jmo34 = jmo12 + nmo(3)+nmo(4)
         !Register the dm initializers as methods
         call dm_init(1)%init_dm_ctor(mo_list(imo12:jmo12),nmo(1:2))
         ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),mo_space_id(1:2)),dm_init(1))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
         call dm_init(2)%init_dm_ctor(mo_list(imo34:jmo34),nmo(3:4))
         ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),mo_space_id(3:4)),dm_init(2))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         if (my_MPI_rank == 0) print*," registered MO vectors spaces and initializers with nmo:",nmo

         !Define separate initializer for the construction of the half-transformed integrals (algorithm 3)
         imo12 = 1
         jmo12 = nmo(1)+nmo(2)
         !Register the ht initializer as method
         call ht_init%init_ht_ctor(mo_list(imo12:jmo12),nmo(1:2))
         ierr=exatns_method_register('HalfTrans',ht_init)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         !Register the ao integral generation as a method
         ierr=exatns_method_register('AointCalculator',aoint_calculator)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         !Register the standard print method with a high treshold to print only important tensor elements
         call tens_printer%reset_thresh(0.01D0)
         ierr=exatns_method_register('TensorPrinter',tens_printer)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register() failed!')

!        All preliminary work is now done on all MPI-nodes, we may proceed to start ExaTensor.
!        After the start of ExaTensor, the driver process (the last MPI process) will be executing
!        the commands and tell the others (that stay inside exatns_process_role) what to do. 
!        These other nodes return only after exatns_stop() is called by the driver.

!        Start ExaTENSOR within MPI_COMM_WORLD (called global_communicator in the MPI interface that DIRAC uses):
         ierr=exatns_start(int(global_communicator,INTD))

         if(ierr.eq.EXA_SUCCESS) then
          ierr=exatns_process_role(my_role) ! only EXA_DRIVER will return immediately
          if(my_role.eq.EXA_DRIVER) then

           call print_date('starting ao integral calculation')

           !Create ao-integral-tensor:
           ao_space_root=ao_space%get_root_id(ierr)
           ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')

           !Initialize ao-tensor:
           ierr=exatns_tensor_init(aoint_tensor,'AointCalculator')
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           call print_date('finished ao integral calculation')

           ! the  next three lines are for debugging,
           print*, "aoint = ", print_tensornorm2 (aoint_tensor) ! this is usually sufficent for checking
           ! when uncommenting the next line: keep in mind that printing an entire tensor yields an enormous output
           ! ierr=exatns_tensor_traverse(aoint_tensor,'TensorPrinter')

           !Create mo-integral-tensor:
           do i = 1, 4
              mo_space_root(i)=mo_space(i)%ptr%get_root_id(ierr)
           end do
           ierr=exatns_tensor_create(moint_tensor,'motensor',mo_space_id,mo_space_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           call print_date('created mo integral tensor')

           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call aomo_transform_1 (aoint_tensor,moint_tensor,print_level)
           call print_date('completed mo integral tensor: algorithm 1')
           print*, "moint = ", print_tensornorm2 (moint_tensor)

           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call aomo_transform_2 (aoint_tensor,moint_tensor,print_level)
           call print_date('completed mo integral tensor: algorithm 2')
           print*, "moint = ", print_tensornorm2 (moint_tensor)

           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call aomo_transform_3 (ao_space_id,ao_space_root,moint_tensor,print_level)
           call print_date('completed mo integral tensor: algorithm 3')
           print*,"moint = ", print_tensornorm2 (moint_tensor)

           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call get_integral_tensor (aoint_tensor,moint_tensor,0,print_level)
           call print_date('completed dirac-sorted mo integral tensor - not anti-symmetrized')
           print*,"moint = ", print_tensornorm2 (moint_tensor)

           ierr=exatns_tensor_destroy(moint_tensor)

           !Recreate mo-integral-tensor with equal index spaces so that we can anti-symmetrize:
           do i = 2, 4
              mo_space_id(i)=mo_space_id(1)
              mo_space_root(i)=mo_space_root(1)
           end do
           ierr=exatns_tensor_create(moint_tensor,'motensor',mo_space_id,mo_space_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           
           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call get_integral_tensor (aoint_tensor,moint_tensor,12,print_level)
           call print_date('completed dirac-sorted mo integral tensor - 12 anti-symmetrized')
           print*,"moint = ", print_tensornorm2 (moint_tensor)
           
           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call get_integral_tensor(aoint_tensor,moint_tensor,34,print_level)
           call print_date('completed dirac-sorted mo integral tensor - 34 anti-symmetrized')
           print*, "moint = ", print_tensornorm2 (moint_tensor)

           !Destroy aoint-tensor:
           ierr=exatns_tensor_destroy(aoint_tensor)

           !Destroy moint-tensor:
           ierr=exatns_tensor_destroy(moint_tensor)

           !Stop ExaTENSOR runtime:
           ierr=exatns_stop()

          end if

         else
          write(*,*) ' Process ',my_MPI_rank,' terminated with error ',ierr
         endif

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         print*, " Process ",my_MPI_rank, " leaving compute_mo_integrals_exa routine"
         return

        end subroutine exacorr_compute_mo_integrals

        subroutine aomo_transform_1 (aoint_tensor,moint_tensor,print_level)

!        Transform integral tensor from ao-basis to mo-basis
!        Algorithm 1: n^5 scaling with index-by-index transform

         implicit none

         ! AO integrals are to be provided on input and are not changed, MO integrals are returned on output
         type(tens_rcrsv_t), intent(inout)  :: aoint_tensor, moint_tensor   ! for technical reasons intent(inout) needs to be used for all tensors
         integer, intent(in)                :: print_level
         ! Auxilliary tensors that are created and destroyed inside this subroutine
         type(tens_rcrsv_t):: a2int_tensor

         integer(INTD)     :: ierr
         integer(INTD)     :: aoint_id(4),moint_id(4),a2int_id(4)
         integer(INTL)     :: aoint_root(4),moint_root(4),a2int_root(4)

         if (print_level >5) call print_date('entered aomo_transform_1')

!        Determine space ids and root
         call aoint_tensor%get_space_ids(aoint_id,aoint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('aoint tensor corrupted')
         call moint_tensor%get_space_ids(moint_id,moint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('moint tensor corrupted')

!        Copy the information about spaces from the starting and target tensors to the auxilliary one and allocate it
         a2int_id(1:2) = moint_id(1:2)
         a2int_id(3:4) = aoint_id(3:4)
         a2int_root(1:2) = moint_root(1:2)
         a2int_root(3:4) = aoint_root(3:4)
         ierr=exatns_tensor_create(a2int_tensor,'a2int',a2int_id,a2int_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(a2int_tensor,'ZERO')

!        Do 1st and 2nd quarter-transformation
         call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
         if (print_level >5) call print_date('first ht done')
!        Do 3rd and 4th quarter-transformation
         call ao2mo_exat (a2int_tensor,moint_tensor,print_level)
         if (print_level >5) call print_date('second ht done')

!        Destruct a2 integrals:
         ierr=exatns_tensor_destroy(a2int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=35; return; endif

         if (print_level >5) call print_date('leaving aomo_transform_1')

        end subroutine aomo_transform_1

        subroutine aomo_transform_2 (aoint_tensor,moint_tensor,print_level)

!        Transform integral tensor from ao-basis to mo-basis
!        Algorithm 2: n^6 scaling with density matrices

         implicit none

         ! AO integrals are to be provided on input and are not changed, MO integrals are returned on output
         type(tens_rcrsv_t), intent(inout)  :: aoint_tensor, moint_tensor
         integer,intent(in)                 :: print_level

         ! Auxilliary tensor needed for this algorithm
         type(tens_rcrsv_t):: dm_tensor
         type(tens_rcrsv_t):: a2int_tensor

         integer(INTD)     :: ierr
         integer(INTD)     :: aoint_id(4),moint_id(4),a2int_id(4),dm_id(4)
         integer(INTL)     :: aoint_root(4),moint_root(4),a2int_root(4),dm_root(4)

         if (print_level >5) call print_date('entered aomo_transform_2')

!        Determine ids and root
         call aoint_tensor%get_space_ids(aoint_id,aoint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('aoint tensor corrupted')
         call moint_tensor%get_space_ids(moint_id,moint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('moint tensor corrupted')

!        Create and initialize the density matrix tensor for the first half-transformation
         dm_id(1:2)    = aoint_id(1:2)
         dm_id(3:4)    = moint_id(1:2)
         dm_root(1:2)  = aoint_root(1:2)
         dm_root(3:4)  = moint_root(1:2)
         ierr=exatns_tensor_create(dm_tensor,'dm12tensor',dm_id,dm_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create(dm1) failed!')
!        Initialize dm-tensor 
         ierr=exatns_tensor_init(dm_tensor,dm_initlabel(dm_id(1:2),dm_id(3:4)))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')

         if (print_level >5) call print_date('initialized first density matrix')

!        Copy the information about spaces from the starting and target tensors to the ht integral tensor
         a2int_id(1:2) = moint_id(1:2)
         a2int_id(3:4) = aoint_id(3:4)
         a2int_root(1:2) = moint_root(1:2)
         a2int_root(3:4) = aoint_root(3:4)

!        Do 1st half-transformation
         ierr=exatns_tensor_create(a2int_tensor,'a2int',a2int_id,a2int_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('ERROR in aomo_transform_2: int. tens. too large ')
         ierr=exatns_tensor_init(a2int_tensor,'ZERO')
         if(ierr.ne.EXA_SUCCESS) call quit('ERROR in aomo_transform_2: initialization failed ')
         if (print_level >5) call print_date('constructed a2')
         ierr=exatns_tensor_contract(a2int_tensor,aoint_tensor,dm_tensor,"D(i,j,r,s)+=L(p,q,r,s)*R(p,q,i,j)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=1; return; endif

!        Destruct dm tensor:
         ierr=exatns_tensor_destroy(dm_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif

         if (print_level >5) call print_date('first ht done')

!        Create and initialize the density matrix tensor for the second half-transformation
         dm_id(1:2)    = aoint_id(3:4)
         dm_id(3:4)    = moint_id(3:4)
         dm_root(1:2)  = aoint_root(3:4)
         dm_root(3:4)   = moint_root(3:4)
         ierr=exatns_tensor_create(dm_tensor,'dm34tensor',dm_id,dm_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create(dm) failed!')
!        Initialize dm-tensor 
         ierr=exatns_tensor_init(dm_tensor,dm_initlabel(dm_id(1:2),dm_id(3:4)))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')

         if (print_level >5) call print_date('initialized second density matrix')

!        Do 2nd half-transformation
         ierr=exatns_tensor_contract(moint_tensor,a2int_tensor,dm_tensor,"D(i,j,k,l)+=L(i,j,r,s)*R(r,s,k,l)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=2; return; endif

         if (print_level >5) call print_date('second ht done')

!        Destruct a2 integrals:
         ierr=exatns_tensor_destroy(a2int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif
!        Destruct dm tensor:
         ierr=exatns_tensor_destroy(dm_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif

         if (print_level >5) call print_date('leaving aomo_transform_2')

        end subroutine aomo_transform_2

        subroutine aomo_transform_3 (ao_space_id,ao_space_root,moint_tensor,print_level)

!        Transform integral tensor from ao-basis to mo-basis
!        Algorithm 3: n^6 scaling with fock matrix builds

         implicit none

         ! MO integrals are returned on output
         type(tens_rcrsv_t), intent(inout)  :: moint_tensor

         ! Auxilliary tensor needed for this algorithm
         type(tens_rcrsv_t):: dm_tensor
         type(tens_rcrsv_t):: a2int_tensor

         integer,intent(in) :: print_level
         integer(INTD)     :: ierr
         integer(INTD)     :: ao_space_id, moint_id(4),a2int_id(4),dm_id(4)
         integer(INTL)     :: ao_space_root,moint_root(4),a2int_root(4),dm_root(4)

         if (print_level >5) call print_date('entered aomo_transform_3')

!        Determine ids and root
         call moint_tensor%get_space_ids(moint_id,moint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('moint tensor corrupted')

!        Copy the information about spaces from the starting and target tensors to the ht integral tensor
         a2int_id(1:2) = ao_space_id
         a2int_id(3:4) = moint_id(3:4)
         a2int_root(1:2) = ao_space_root
         a2int_root(3:4) = moint_root(3:4)

!        Initialize half-transformed integrals
         ierr=exatns_tensor_create(a2int_tensor,'a2int',a2int_id,a2int_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(a2int_tensor,'HalfTrans')
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=1; return; endif

         if (print_level >5) call print_date('first ht done')

!        Create the density matrix tensor for the second half-transformation
         dm_id(1:2)    = ao_space_id
         dm_id(3:4)    = moint_id(1:2)
         dm_root(1:2)  = ao_space_root
         dm_root(3:4)   = moint_root(1:2)
         ierr=exatns_tensor_create(dm_tensor,'dm12tensor',dm_id,dm_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create(dm) failed!')
!        Initialize dm-tensor 
         ierr=exatns_tensor_init(dm_tensor,dm_initlabel(dm_id(1:2),dm_id(3:4)))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')

         if (print_level >5) call print_date('initialized density matrix')

!        Do 2nd half-transformation
         ierr=exatns_tensor_contract(moint_tensor,a2int_tensor,dm_tensor,"D(i,j,k,l)+=L(p,q,k,l)*R(p,q,i,j)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=2; return; endif

         if (print_level >5) call print_date('second ht done')

!        Destruct a2 integrals:
         ierr=exatns_tensor_destroy(a2int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif
!        Destruct dm tensor:
         ierr=exatns_tensor_destroy(dm_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif

         if (print_level >5) call print_date('leaving aomo_transform_3')

        end subroutine aomo_transform_3

        subroutine get_integral_tensor (aoint_tensor,integral_tensor,antisymmetrize,print_level,algorithm)

         implicit none

         type(tens_rcrsv_t), intent(inout)  :: aoint_tensor, integral_tensor
         type(tens_rcrsv_t)                 :: moint_tensor

         integer, intent(in)                :: antisymmetrize
         integer, optional                  :: algorithm
         integer, intent(in)                :: print_level

         integer(INTD)     :: ierr, tens_rank
         integer(INTD)     :: aoint_id(4),integral_id(4),moint_id(4)
         integer(INTL)     :: aoint_root(4),integral_root(4),moint_root(4)

         if (.not.present(algorithm)) algorithm = 1

!        Check integrity of tensors
         tens_rank = aoint_tensor%get_rank(ierr)
         if (ierr.ne.0 .or. tens_rank.ne.4) call quit('input tensor corrupted in aomo_transform')
         tens_rank = integral_tensor%get_rank(ierr)
         if (ierr.ne.0 .or. tens_rank.ne.4) call quit('output tensor corrupted in aomo_transform')
         
!        Determine ids and root
         call aoint_tensor%get_space_ids(aoint_id,aoint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('aoint tensor corrupted')
         call integral_tensor%get_space_ids(integral_id,integral_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('integral tensor corrupted')

!        Copy the information about spaces from the starting and target tensors to the auxilliary ones
!        Use relation (12|34) = <13|24>
         moint_id(1) = integral_id(1)
         moint_id(2) = integral_id(3)
         moint_id(3) = integral_id(2)
         moint_id(4) = integral_id(4)

         moint_root(1) = integral_root(1)
         moint_root(2) = integral_root(3)
         moint_root(3) = integral_root(2)
         moint_root(4) = integral_root(4)

!        Allocate auxilliary tensor to get integrals (without anti-symmetrization) in Mulliken order
         ierr=exatns_tensor_create(moint_tensor,'mulliken',moint_id,moint_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(moint_tensor,'ZERO')
        
!        Transform integrals
         select case (algorithm)
         case (1)
            call aomo_transform_1 (aoint_tensor,moint_tensor,print_level)
         case (2)
            call aomo_transform_2 (aoint_tensor,moint_tensor,print_level)
         case default
            call quit ("asked for non-existing algorithm in get_integral_tensor")
         end select

!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)
         ierr=exatns_tensor_destroy(moint_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=1; return; endif

        end subroutine get_integral_tensor

        subroutine get_ht2_integral_tensor (a2int_tensor,integral_tensor,antisymmetrize,print_level)

         implicit none

         type(tens_rcrsv_t), intent(inout)  :: a2int_tensor, integral_tensor
         type(tens_rcrsv_t)                 :: moint_tensor

         integer, intent(in)                :: antisymmetrize, print_level

         integer(INTD)     :: ierr, tens_rank
         integer(INTD)     :: a2int_id(4),integral_id(4),moint_id(4)
         integer(INTL)     :: a2int_root(4),integral_root(4),moint_root(4)

!        Check integrity of tensors
         tens_rank = a2int_tensor%get_rank(ierr)
         if (ierr.ne.0 .or. tens_rank.ne.4) call quit('input tensor corrupted in aomo_transform')
         tens_rank = integral_tensor%get_rank(ierr)
         if (ierr.ne.0 .or. tens_rank.ne.4) call quit('output tensor corrupted in aomo_transform')

!        Determine ids and root
         call a2int_tensor%get_space_ids(a2int_id,a2int_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('a2int tensor corrupted')
         call integral_tensor%get_space_ids(integral_id,integral_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('integral tensor corrupted')

!        Copy the information about spaces from the starting and target tensors to the auxilliary ones
!        Use relation (12|34) = <13|24>
         moint_id(1) = integral_id(1)
         moint_id(2) = integral_id(3)
         moint_id(3) = integral_id(2)
         moint_id(4) = integral_id(4)

         moint_root(1) = integral_root(1)
         moint_root(2) = integral_root(3)
         moint_root(3) = integral_root(2)
         moint_root(4) = integral_root(4)

!        Allocate auxilliary tensor to get integrals (without anti-symmetrization) in Mulliken order
         ierr=exatns_tensor_create(moint_tensor,'mulliken',moint_id,moint_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(moint_tensor,'ZERO')

!        Transform integrals
         call ao2mo_exat (a2int_tensor,moint_tensor,print_level)
!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)
         ierr=exatns_tensor_destroy(moint_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=1; return; endif

        end subroutine get_ht2_integral_tensor

        subroutine mulliken_to_dirac_sort  (mulliken_tensor,dirac_tensor,antisymmetrize)

         implicit none

         type(tens_rcrsv_t), intent(inout)  :: mulliken_tensor,dirac_tensor
         integer, intent(in)                :: antisymmetrize


         integer(INTD)   :: mulliken_id(4),dirac_id(4)
         integer(INTL)   :: mulliken_root(4),dirac_root(4)

         integer(INTD)   :: ierr, mismatch

         type(tens_rcrsv_t) :: one_tensor

!        get dimensions of the input tensors and check consistency
         call mulliken_tensor%get_space_ids(mulliken_id,mulliken_root,ierr)
         call dirac_tensor%get_space_ids(dirac_id,dirac_root,ierr)
         mismatch = 0
         if (mulliken_id(1) .ne. dirac_id(1)) mismatch = mismatch + 1
         if (mulliken_id(2) .ne. dirac_id(3)) mismatch = mismatch + 1
         if (mulliken_id(3) .ne. dirac_id(2)) mismatch = mismatch + 1
         if (mulliken_id(4) .ne. dirac_id(4)) mismatch = mismatch + 1
         if (antisymmetrize .eq. 12 .and. dirac_id(1) .ne. dirac_id(2)) mismatch = mismatch + 1
         if (antisymmetrize .eq. 34 .and. dirac_id(3) .ne. dirac_id(4)) mismatch = mismatch + 1
         if (mismatch.gt.0) stop 'spaces mismatch in mulliken_to_dirac_sort'

         ! We first need to define a scalar as a tensor to work around the missing tensor addition
         ierr=exatns_tensor_create(one_tensor,'one_tensor',EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create(one_tens) failed!')
         ierr=exatns_tensor_init(one_tensor,ONE)

         ! Now do the actual sort, note that we do not initialize as we may be accumulating in the target tensor
         ierr=exatns_tensor_contract(dirac_tensor,mulliken_tensor,one_tensor,"D(p,q,r,s)+=M(p,r,q,s)*C()")

         if (antisymmetrize .eq. 12) then
             ierr=exatns_tensor_contract(dirac_tensor,mulliken_tensor,one_tensor,"D(p,q,r,s)+=M(q,r,p,s)*C()",prefactor=MINUS_ONE)
         elseif (antisymmetrize .eq. 34) then
             ierr=exatns_tensor_contract(dirac_tensor,mulliken_tensor,one_tensor,"D(p,q,r,s)+=M(p,s,q,r)*C()",prefactor=MINUS_ONE)
         end if

         ierr=exatns_tensor_destroy(one_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=1; return; endif

        end subroutine mulliken_to_dirac_sort

        function print_tensornorm2 (tens) result(scalar_real)

          type(tens_rcrsv_t) :: tens
          type(tens_rcrsv_t) :: tensor_scalar
          integer(INTD)      :: ierr,rank
          complex(8)         :: scalar_complex
          real(8)            :: scalar_real

          ! create a rank-0 tensor
          ierr=exatns_tensor_create(tensor_scalar,'tensor_scalar',EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
          ierr=exatns_tensor_init(tensor_scalar,ZERO)

          ! get the rank and do the contraction
          rank = tens%get_rank(ierr)
          select case (rank)
          case(0)
              ierr=exatns_tensor_contract(tensor_scalar,tens,tens,'E()+=D+()*D()')
          case(1)
              ierr=exatns_tensor_contract(tensor_scalar,tens,tens,'E()+=D+(a)*D(a)')
          case(2)
              ierr=exatns_tensor_contract(tensor_scalar,tens,tens,'E()+=D+(a,b)*D(a,b)')
          case(3)
              ierr=exatns_tensor_contract(tensor_scalar,tens,tens,'E()+=D+(a,b,c)*D(a,b,c)')
          case(4)
              ierr=exatns_tensor_contract(tensor_scalar,tens,tens,'E()+=D+(a,b,c,d)*D(a,b,c,d)')
          case(5)
              ierr=exatns_tensor_contract(tensor_scalar,tens,tens,'E()+=D+(a,b,c,d,e)*D(a,b,c,d,e)')
          case(6)
              ierr=exatns_tensor_contract(tensor_scalar,tens,tens,'E()+=D+(a,b,c,d,e,f)*D(a,b,c,d,e,f)')
          case default
              ierr=1
          end select
          if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_contract() failed!')
          ierr=exatns_tensor_get_scalar(tensor_scalar,scalar_complex)
          if(ierr.ne.EXA_SUCCESS) call quit('print_tensornorm2: get_scalar failed!')
          scalar_real=real(scalar_complex,8)
          ierr=exatns_tensor_destroy(tensor_scalar)
          if(ierr.ne.EXA_SUCCESS) then; ierr=45; return; endif

        end function print_tensornorm2
 
        function mocoef_initlabel (ao_space_id,mo_space_id,spin) result(initlabel)

         !Define a unique name for the initializer of mo coeffiecients for a particular spin
         !we assume that we can use a single digit to indicate the ao and mo spaces
         !if we start using more spaces, we can simply define a longer string here

         integer(INTD),intent(in):: ao_space_id, mo_space_id, spin
         character(len=9)  :: initlabel
         character(len=2)  :: vs_str

         if (ao_space_id>9.or.mo_space_id>9) then
             call quit ('Revise naming scheme for MO coef initializers')
         else
              write (vs_str,'(2I1)') ao_space_id,mo_space_id
         end if
         if (spin == 1) initlabel = 'MOCoef' // vs_str // 'A'
         if (spin == 2) initlabel = 'MOCoef' // vs_str // 'B'

        end function mocoef_initlabel

        function dm_initlabel (ao_space_id,mo_space_id) result(initlabel)

         !Define a unique name for the initializer of density matrices
         !we assume that we can use a single digit to indicate the ao and mo spaces
         !if we start using more spaces, we can simply define a longer string here

         integer(INTD),intent(in):: ao_space_id(2), mo_space_id(2)
         character(len=9)        :: initlabel
         character(len=2)        :: vs_str(2)
         integer                 :: i

         do i = 1, 2
            if (ao_space_id(i)>9.or.mo_space_id(i)>9) then
               call quit ('Revise naming scheme for DM initializers')
            else
              write (vs_str(i),'(2I1)') ao_space_id(i),mo_space_id(i)
           end if
         end do
         initlabel = 'DensM' // vs_str(1) // vs_str(2)

        end function dm_initlabel

      subroutine ao2mo_exat(ao_tensor,mo_tensor,print_level)

!       Transform integral tensor from ao-basis to mo-basis

        implicit none

        type(tens_rcrsv_t), intent(inout) :: ao_tensor, mo_tensor 
        integer,intent(in)                :: print_level

        ! Auxilliary tensors that are created and destroyed inside this subroutine
        type(tens_rcrsv_t):: coef_tensor
        type(tens_rcrsv_t):: t_tensor

        integer           :: i, ind_a, ind_b
        integer(INTD)     :: ierr,spin, aorank, morank
        integer(INTD)     :: ao4_id(4),ao3_id(3),ao2_id(2)
        integer(INTD)     :: mo4_id(4),mo3_id(3),mo2_id(2)
        integer(INTD)     :: t4_id(4), t3_id(3), t2_id(2)
        integer(INTD)     :: coef_id(2)
        integer(INTL)     :: ao4_root(4),ao3_root(3),ao2_root(2)
        integer(INTL)     :: mo4_root(4),mo3_root(3),mo2_root(2)
        integer(INTL)     :: t4_root(4), t3_root(3), t2_root(2)
        integer(INTL)     :: coef_root(2)
        character(len=9)  :: coef_name
        character(len=38) :: info_str

        if (print_level >5) call print_date('=== entered ao2mo_exat')

!       Determine ranks
        aorank = ao_tensor%get_rank(ierr)
        if (ierr.ne.0) call quit('ao2mo_exat: ao tensor corrupted (rank)')
        morank = mo_tensor%get_rank(ierr)
        if (ierr.ne.0 .or. aorank.ne.morank) call quit('ao2mo_exat: mo tensor corrupted (rank) ')

        if (aorank.eq.2) then
          call ao_tensor%get_space_ids(ao2_id,ao2_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: ao tensor corrupted (2) ')
          call mo_tensor%get_space_ids(mo2_id,mo2_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: mo tensor corrupted (2)')

          !create temporary tensor
          t2_id(1)   = mo2_id(1)
          t2_id(2)   = ao2_id(2)
          t2_root(1) = mo2_root(1)
          t2_root(2) = ao2_root(2)

          ierr=exatns_tensor_create(t_tensor,'temp',t2_id,t2_root,EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: temp2 tensor not created')

          do spin = 1, 2
            ierr=exatns_tensor_init(t_tensor,'ZERO')

            i=1 ! first index
            coef_id(1)   = ao2_id(i)
            coef_root(1) = ao2_root(i)
            coef_id(2)   = mo2_id(i)
            coef_root(2) = mo2_root(i)

            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor,coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor,mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>10) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor)

            ierr=exatns_tensor_contract(t_tensor,ao_tensor,coef_tensor,"B(p,s)+=A(r,s)*M+(r,p)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 1 failed')

            ierr=exatns_tensor_destroy(coef_tensor)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

            if (print_level >5) then
              write (info_str,'(A31,I2,A3,I2)') 'half(2) transformation done: s=',spin,',i=',i
              call print_date(info_str)
            end if

            i=2  ! second index
            coef_id(1)   = ao2_id(i)
            coef_root(1) = ao2_root(i)
            coef_id(2)   = mo2_id(i)
            coef_root(2) = mo2_root(i)

            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor,coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor,mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>10) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor)

            ierr=exatns_tensor_contract(mo_tensor,t_tensor,coef_tensor,"C(p,q)+=B(p,s)*M(s,q)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 2 failed')

            ierr=exatns_tensor_destroy(coef_tensor)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

            if (print_level >5) then
              write (info_str,'(A31,I2,A3,I2)') 'quarter transformation done: s=',spin,',i=',i
              call print_date(info_str)
            end if

          end do

          ierr=exatns_tensor_destroy(t_tensor)

        else if (aorank.eq.3) then
          call quit('ao2mo_exat: for now rank 3 not implemented')

        else if (aorank.eq.4) then
          ! Algorithm 1: n^5 scaling with index-by-index transform
          call ao_tensor%get_space_ids(ao4_id,ao4_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: ao tensor corrupted (4)')
          call mo_tensor%get_space_ids(mo4_id,mo4_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: mo tensor corrupted (4)')

          ! create temporary tensor
          if (ao4_id(1).eq.ao4_id(3) .and. ao4_id(2).eq.ao4_id(4)) then
            !first halftransformation
            t4_id(1)     = mo4_id(1)
            t4_id(2:4)   = ao4_id(2:4)
            t4_root(1)   = mo4_root(1)
            t4_root(2:4) = ao4_root(2:4)
            ind_a=1
            ind_b=2
          else
            t4_id(1:3)   = mo4_id(1:3)
            t4_id(4)     = ao4_id(4)
            t4_root(1:3) = mo4_root(1:3)
            t4_root(4)   = ao4_root(4)
            ind_a=3
            ind_b=4
          end if

          ierr=exatns_tensor_create(t_tensor,'temp',t4_id,t4_root,EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: temp4 tensor not created')

          do spin = 1, 2
            ierr=exatns_tensor_init(t_tensor,'ZERO')

            i=ind_a ! first index
            coef_id(1)   = ao4_id(i)
            coef_root(1) = ao4_root(i)
            coef_id(2)   = mo4_id(i)
            coef_root(2) = mo4_root(i)
          
            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor,coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor,mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>10) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor)

            if (ind_a.eq.1) then
              ierr=exatns_tensor_contract(t_tensor,coef_tensor,ao_tensor,"D(i,q,r,s)+=R+(p,i)*L(p,q,r,s)")
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 1 failed')
            else
              ierr=exatns_tensor_contract(t_tensor,coef_tensor,ao_tensor,"D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)")
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 3 failed')
            end if
            
            ierr=exatns_tensor_destroy(coef_tensor)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

            if (print_level >5) then
              write (info_str,'(A31,I2,A3,I2)') 'quarter transformation done: s=',spin,',i=',i
              call print_date(info_str)
            end if

            i=ind_b  ! second index
            coef_id(1)   = ao4_id(i)
            coef_root(1) = ao4_root(i)
            coef_id(2)   = mo4_id(i)
            coef_root(2) = mo4_root(i)

            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor,coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor,mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>10) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor)

            if (ind_b.eq.2) then
              ierr=exatns_tensor_contract(mo_tensor,t_tensor,coef_tensor,"D(i,j,r,s)+=L(i,q,r,s)*R(q,j)")
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 2 failed')
            else
              ierr=exatns_tensor_contract(mo_tensor,t_tensor,coef_tensor,"D(i,j,k,l)+=L(i,j,k,s)*R(s,l)")
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 4 failed')
            end if

            ierr=exatns_tensor_destroy(coef_tensor)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

            if (print_level >5) then
              write (info_str,'(A31,I2,A3,I2)') 'quarter transformation done: s=',spin,',i=',i
              call print_date(info_str)
            end if
          end do
          
          ierr=exatns_tensor_destroy(t_tensor)

        else
          call quit('ao2mo_exat: wrong rank or mode')
        end if

        if (print_level >5) call print_date('=== leaving ao2mo_exat')

      end subroutine ao2mo_exat

      subroutine ao2mo_exat_parallel(ao_tensor,mo_tensor,print_level)

!       Transform integral tensor from ao-basis to mo-basis

        implicit none

        type(tens_rcrsv_t), intent(inout) :: ao_tensor, mo_tensor 
        integer,intent(in)                :: print_level

        ! Auxilliary tensors that are created and destroyed inside this subroutine
        type(tens_rcrsv_t):: coef_tensor(2)
        type(tens_rcrsv_t):: t_tensor(2)

        integer           :: i, ind_a, ind_b
        integer(INTD)     :: ierr,spin, aorank, morank
        integer(INTD)     :: ao4_id(4),ao3_id(3)
        integer(INTD),allocatable :: ao_id(:)
        integer(INTD)     :: mo4_id(4),mo3_id(3),mo2_id(2)
        integer(INTD)     :: t4_id(4), t3_id(3), t2_id(2)
        integer(INTD)     :: coef_id(2)
        integer(INTL)     :: ao4_root(4),ao3_root(3),ao2_root(2)
        integer(INTL)     :: mo4_root(4),mo3_root(3),mo2_root(2)
        integer(INTL)     :: t4_root(4), t3_root(3), t2_root(2)
        integer(INTL)     :: coef_root(2)
        character(len=9)  :: coef_name
        character(len=38) :: info_str

        if (print_level >10) call print_date('=== entered ao2mo_exat')

!       Determine ranks
        aorank = ao_tensor%get_rank(ierr)
        if (ierr.ne.0) call quit('ao2mo_exat: ao tensor corrupted (rank)')
        morank = mo_tensor%get_rank(ierr)
        if (ierr.ne.0 .or. aorank.ne.morank) call quit('ao2mo_exat: mo tensor corrupted (rank) ')

        if (aorank.eq.2) then
          allocate(ao_id(2))

          call ao_tensor%get_space_ids(ao_id,ao2_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: ao tensor corrupted (2) ')
          call mo_tensor%get_space_ids(mo2_id,mo2_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: mo tensor corrupted (2)')

          !create temporary tensor
          t2_id(1)   = mo2_id(1)
          t2_id(2)   = ao_id(2)
          t2_root(1) = mo2_root(1)
          t2_root(2) = ao2_root(2)

          ierr=exatns_tensor_create(t_tensor(1),'temp',t2_id,t2_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_create(t_tensor(2),'temp',t2_id,t2_root,EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: temp2 tensor not created')
          ierr=exatns_tensor_init(t_tensor(1),'ZERO')
          ierr=exatns_tensor_init(t_tensor(2),'ZERO')
          
          i=1 ! first index
          coef_id(1)   = ao_id(i)
          coef_root(1) = ao2_root(i)
          coef_id(2)   = mo2_id(i)
          coef_root(2) = mo2_root(i)

          do spin=1,2
            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor(spin),coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor(spin),mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>14) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor(spin))
          end do

          ierr=exatns_tensor_contract(t_tensor(1),ao_tensor,coef_tensor(1),"B(p,s)+=A(r,s)*M+(r,p)")
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 11 failed')
          ierr=exatns_tensor_contract(t_tensor(2),ao_tensor,coef_tensor(2),"B(p,s)+=A(r,s)*M+(r,p)")
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 12 failed')

          ierr=exatns_tensor_destroy(coef_tensor(1))
          ierr=exatns_tensor_destroy(coef_tensor(2))
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

          if (print_level >12) then
            write (info_str,'(A31,I2)') 'half(2) transformation done: i=',i
            call print_date(info_str)
          end if

          i=2  ! second index
          coef_id(1)   = ao_id(i)
          coef_root(1) = ao2_root(i)
          coef_id(2)   = mo2_id(i)
          coef_root(2) = mo2_root(i)

          do spin=1,2
            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor(spin),coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor(spin),mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>14) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor(spin))
          end do

          ierr=exatns_tensor_contract(mo_tensor,t_tensor(1),coef_tensor(1),"C(p,q)+=B(p,s)*M(s,q)")
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 21 failed')
          ierr=exatns_tensor_contract(mo_tensor,t_tensor(2),coef_tensor(2),"C(p,q)+=B(p,s)*M(s,q)")
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 21 failed')

          ierr=exatns_tensor_destroy(coef_tensor(1))
          ierr=exatns_tensor_destroy(coef_tensor(2))
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

          if (print_level >12) then
            write (info_str,'(A31,I2)') 'quarter transformation done: i=',i
            call print_date(info_str)
          end if

          ierr=exatns_tensor_destroy(t_tensor(1))
          ierr=exatns_tensor_destroy(t_tensor(2))
          deallocate(ao_id)

        else if (aorank.eq.3) then
          call quit('ao2mo_exat: for now rank 3 not implemented')

        else if (aorank.eq.4) then
          ! Algorithm 1: n^5 scaling with index-by-index transform
          call ao_tensor%get_space_ids(ao4_id,ao4_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: ao tensor corrupted (4)')
          call mo_tensor%get_space_ids(mo4_id,mo4_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: mo tensor corrupted (4)')

          ! create temporary tensor
          if (ao4_id(1).eq.ao4_id(3) .and. ao4_id(2).eq.ao4_id(4)) then
            !first halftransformation
            t4_id(1)     = mo4_id(1)
            t4_id(2:4)   = ao4_id(2:4)
            t4_root(1)   = mo4_root(1)
            t4_root(2:4) = ao4_root(2:4)
            ind_a=1
            ind_b=2
          else
            t4_id(1:3)   = mo4_id(1:3)
            t4_id(4)     = ao4_id(4)
            t4_root(1:3) = mo4_root(1:3)
            t4_root(4)   = ao4_root(4)
            ind_a=3
            ind_b=4
          end if

          ierr=exatns_tensor_create(t_tensor(1),'temp',t4_id,t4_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_create(t_tensor(2),'temp',t4_id,t4_root,EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: temp2 tensor not created')
          ierr=exatns_tensor_init(t_tensor(1),'ZERO')
          ierr=exatns_tensor_init(t_tensor(2),'ZERO')

          i=ind_a ! first index
          coef_id(1)   = ao4_id(i)
          coef_root(1) = ao4_root(i)
          coef_id(2)   = mo4_id(i)
          coef_root(2) = mo4_root(i)
          
          do spin=1,2
            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor(spin),coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor(spin),mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>14) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor(spin))
          end do

          if (ind_a.eq.1) then
            ierr=exatns_tensor_contract(t_tensor(1),coef_tensor(1),ao_tensor,"D(i,q,r,s)+=R+(p,i)*L(p,q,r,s)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 11 failed')
            ierr=exatns_tensor_contract(t_tensor(2),coef_tensor(2),ao_tensor,"D(i,q,r,s)+=R+(p,i)*L(p,q,r,s)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 12 failed')
          else
            ierr=exatns_tensor_contract(t_tensor(1),coef_tensor(1),ao_tensor,"D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 31 failed')
            ierr=exatns_tensor_contract(t_tensor(2),coef_tensor(2),ao_tensor,"D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 32 failed')
          end if
            
          ierr=exatns_tensor_destroy(coef_tensor(1))
          ierr=exatns_tensor_destroy(coef_tensor(2))
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

          if (print_level >12) then
            write (info_str,'(A31,I2)') 'quarter transformation done: i=',i
            call print_date(info_str)
          end if

          i=ind_b  ! second index
          coef_id(1)   = ao4_id(i)
          coef_root(1) = ao4_root(i)
          coef_id(2)   = mo4_id(i)
          coef_root(2) = mo4_root(i)

          do spin=1,2
            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor(1),coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor(2),mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>14) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor(2))
          end do

          if (ind_b.eq.2) then
            ierr=exatns_tensor_contract(mo_tensor,t_tensor(1),coef_tensor(1),"D(i,j,r,s)+=L(i,q,r,s)*R(q,j)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 21 failed')
            ierr=exatns_tensor_contract(mo_tensor,t_tensor(2),coef_tensor(2),"D(i,j,r,s)+=L(i,q,r,s)*R(q,j)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 22 failed')
          else
            ierr=exatns_tensor_contract(mo_tensor,t_tensor(1),coef_tensor(1),"D(i,j,k,l)+=L(i,j,k,s)*R(s,l)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 41 failed')
            ierr=exatns_tensor_contract(mo_tensor,t_tensor(2),coef_tensor(2),"D(i,j,k,l)+=L(i,j,k,s)*R(s,l)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 42 failed')
          end if

          ierr=exatns_tensor_destroy(coef_tensor(1))
          ierr=exatns_tensor_destroy(coef_tensor(2))
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

          if (print_level >12) then
            write (info_str,'(A31,I2)') 'quarter transformation done: i=',i
            call print_date(info_str)
          end if
          
          ierr=exatns_tensor_destroy(t_tensor(1))
          ierr=exatns_tensor_destroy(t_tensor(2))

        else
          call quit('ao2mo_exat: wrong rank or mode')
        end if

        if (print_level >10) call print_date('=== leaving ao2mo_exat')

      end subroutine ao2mo_exat_parallel

      subroutine get_CC_integrals_chol (ao_id,ao_space,int_t,moint_scheme, &
                                    exa_blocksize,print_level,t_cholesky,j_chol,f_delta)

!       routine to get 2-el integrals using cholesky

!       Written by Johann Pototschnig, Spring 2020

        use exacorr_tensor_methods
        use exacorr_global
        use exacorr_datatypes
        use exacorr_utils

!       information needed to compute the ao-tensor
        integer(INTD),intent(in)              :: ao_id
        class(h_space_t), pointer, intent(in) :: ao_space
!       target integral tensors to be filled and given back to the caller
        type(exatns_intg_tens), intent(inout) :: int_t
        integer                               :: moint_scheme ! (41: fully pvotised, 42: using shells)
        integer, intent(in)                   :: print_level, exa_blocksize
        real(8),intent(in)                    :: t_cholesky
        type(chol_J), intent(inout)           :: j_chol  ! Cholesky: compute off-diagonal elements
        type(delta_t), intent(inout)          :: f_delta  ! fill delta matrix


        integer(INTD)                   :: ierr
        integer(INTL)                   :: nao, nocc, nvir
        integer(INTD)                   :: o_id, v_id
        integer                         :: i, m
        type(tens_rcrsv_t), allocatable :: chol_tens(:), temp_tens(:)
        integer(INTL),dimension(2)      :: dimA, dimB, DimC
        integer(INTD),dimension(2)      :: idA, idB, idC
        type(tens_rcrsv_t), allocatable :: uoo_tens(:), uov_tens(:), uvo_tens(:), uvv_tens(:)
        integer(INTL),dimension(4)      :: vovo_root, ovvo_root, hroot
        integer(INTD),dimension(4)      :: vovo_id, ovvo_id, hid
        type(tens_rcrsv_t)              :: ovvo_tensor, temp_tensor, one_element
        character(len=12)               :: vs_str


!       get the root space, covers all ao's
        nao=ao_space%get_root_id(ierr)
        call int_t%vovo%get_space_ids(vovo_id,vovo_root,ierr)
        if(ierr.ne.EXA_SUCCESS) call quit('vovo tensor corrupted')
        nocc=vovo_root(2)
        nvir=vovo_root(1)
        o_id=vovo_id(2)
        v_id=vovo_id(1)

        call print_date('-Start- Cholesky Decomposition')
        if (moint_scheme.eq.41) then
            call decompose_cholesky_exatensor_full (chol_tens, m, ao_id, ao_space, t_cholesky, &
                                                    j_chol, f_delta, print_level)
        else
            call decompose_cholesky_exatensor (chol_tens, m, ao_id, ao_space, t_cholesky, &
                                               f_delta,exa_blocksize,print_level)  
        end if
        call print_date('-End- Cholesky Decomposition')

        ierr=exatns_tensor_create(one_element,"one_element",EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(one_element,ONE)

!       allocate cholesky vectors
        allocate(uoo_tens(m))
        allocate(uov_tens(m))
        allocate(uvo_tens(m))
        allocate(uvv_tens(m))
        allocate(temp_tens(m))

        !        transform v part
        dimA=nao
        dimA(1)=nvir
        idA=ao_id
        idA(1)=v_id
        dimB=nocc
        dimB(1)=nvir
        idB=o_id
        idB(1)=v_id
        dimC=nvir
        idC=v_id

        do i=1,m
          write (vs_str, '(A3,I9.9)') 'tem',i
          ierr=exatns_tensor_create(temp_tens(i),vs_str,idA,dimA,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_init(temp_tens(i),'ZERO')
          write (vs_str, '(A3,I9.9)') 'uvo',i
          ierr=exatns_tensor_create(uvo_tens(i),vs_str,idB,dimB,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_init(uvo_tens(i),'ZERO')
          write (vs_str, '(A3,I9.9)') 'uvv',i
          ierr=exatns_tensor_create(uvv_tens(i),vs_str,idC,dimC,EXA_DATA_KIND_C8)  
          ierr=exatns_tensor_init(uvv_tens(i),'ZERO')        
        end do

        call ao2mo_vect(chol_tens,temp_tens,m,11,print_level)
        call ao2mo_vect(temp_tens,uvo_tens,m,21,print_level)
        call ao2mo_vect(temp_tens,uvv_tens,m,21,print_level)
        do i=1,m
          ierr=exatns_tensor_init(temp_tens(i),'ZERO')
        end do
        call ao2mo_vect(chol_tens,temp_tens,m,12,print_level)
        call ao2mo_vect(temp_tens,uvo_tens,m,22,print_level)
        call ao2mo_vect(temp_tens,uvv_tens,m,22,print_level)

        if (print_level.gt.2) call print_date(' uvo/uvv done ')

!       transform o part
        dimA=nao
        dimA(1)=nocc
        idA=ao_id
        idA(1)=o_id        
        dimB=nvir
        dimB(1)=nocc
        idB=v_id
        idB(1)=o_id
        dimC=nocc
        idC=o_id

        do i=1,m
          ierr=exatns_tensor_destroy(temp_tens(i))
          write (vs_str, '(A3,I9.9)') 'uov',i
          ierr=exatns_tensor_create(uov_tens(i),vs_str,idB,dimB,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_init(uov_tens(i),'ZERO')
          write (vs_str, '(A3,I9.9)') 'uoo',i
          ierr=exatns_tensor_create(uoo_tens(i),vs_str,idC,dimC,EXA_DATA_KIND_C8) 
          ierr=exatns_tensor_init(uoo_tens(i),'ZERO')
          write (vs_str, '(A3,I9.9)') 'tem',i
          ierr=exatns_tensor_create(temp_tens(i),vs_str,idA,dimA,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_init(temp_tens(i),'ZERO') 
        end do

        call ao2mo_vect(chol_tens,temp_tens,m,11,print_level)
        call ao2mo_vect(temp_tens,uoo_tens,m,21,print_level)
        call ao2mo_vect(temp_tens,uov_tens,m,21,print_level)
        do i=1,m
          ierr=exatns_tensor_init(temp_tens(i),'ZERO')
        end do
        call ao2mo_vect(chol_tens,temp_tens,m,12,print_level)
        call ao2mo_vect(temp_tens,uoo_tens,m,22,print_level)
        call ao2mo_vect(temp_tens,uov_tens,m,22,print_level)

        if (print_level.gt.2) call print_date(' uoo/uvo done ')

        ! cleanup
        do i=1,m
          ierr=exatns_tensor_destroy(temp_tens(i))
          ierr=exatns_tensor_destroy(chol_tens(i))
        end do
        deallocate(temp_tens)
        deallocate(chol_tens)

!       Determine space ids and root for auxilliary voov tensor from vovo tensor
        ovvo_id = o_id
        ovvo_id(2:3) = v_id
        ovvo_root = nocc
        ovvo_root(2:3) = nvir

!       initialise tensors to zero
        ierr=exatns_tensor_init(int_t%oooo,'ZERO')
        ierr=exatns_tensor_init(int_t%ooov,'ZERO')
        ierr=exatns_tensor_init(int_t%oovv,'ZERO')
        ierr=exatns_tensor_init(int_t%vovo,'ZERO')
        ierr=exatns_tensor_init(int_t%vovv,'ZERO')
        ierr=exatns_tensor_init(int_t%vvvv,'ZERO')

!       Get antisymmetric oooo tensor
        hid=o_id
        hroot=nocc
        ierr=exatns_tensor_create(temp_tensor,'temp',hid,hroot,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(temp_tensor,'ZERO')
        do i=1,m
          ierr=exatns_tensor_contract("T(p,q,r,s)+=K+(q,p)*L(r,s)",temp_tensor,uoo_tens(i),uoo_tens(i))
        end do
        if(ierr.ne.EXA_SUCCESS) call quit('cholesky integrals: oooo')
        call mulliken_to_dirac_sort(temp_tensor,int_t%oooo,12)
        ierr=exatns_tensor_destroy(temp_tensor)
        if (print_level.gt.2) call print_date(' oooo done ')

!       Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
        hid=o_id
        hid(1:2)=v_id
        hroot=nocc
        hroot(1:2)=nvir
        ierr=exatns_tensor_create(temp_tensor,'temp',hid,hroot,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(temp_tensor,'ZERO')
        do i=1,m
          ierr=exatns_tensor_contract("T(p,q,r,s)+=K+(q,p)*L(r,s)",temp_tensor,uvv_tens(i),uoo_tens(i))
        end do
        if(ierr.ne.EXA_SUCCESS) call quit('cholesky integrals: vovo')
        call mulliken_to_dirac_sort(temp_tensor,int_t%vovo,0)
        ierr=exatns_tensor_destroy(temp_tensor)
        ierr=exatns_tensor_create(ovvo_tensor,'ovvo',ovvo_id,ovvo_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(ovvo_tensor,'ZERO')
        ierr=exatns_tensor_create(temp_tensor,'temp',ovvo_id,ovvo_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(temp_tensor,'ZERO')
        do i=1,m
          ierr=exatns_tensor_contract("T(p,q,r,s)+=K+(q,p)*L(r,s)",temp_tensor,uvo_tens(i),uvo_tens(i))
        end do
        if(ierr.ne.EXA_SUCCESS) call quit('cholesky integrals: ovvo')
        call mulliken_to_dirac_sort(temp_tensor,ovvo_tensor,0)
        ierr=exatns_tensor_destroy(temp_tensor)
        ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*C()",int_t%vovo,ovvo_tensor,one_element,MINUS_ONE)
        ierr=exatns_tensor_destroy(ovvo_tensor)
        if (print_level.gt.2) call print_date(' vovo done ')

!       Get ooov tensor
        hid=o_id
        hid(4)=v_id
        hroot=nocc
        hroot(4)=nvir
        ierr=exatns_tensor_create(temp_tensor,'temp',hid,hroot,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(temp_tensor,'ZERO')
        do i=1,m
          ierr=exatns_tensor_contract("T(p,q,r,s)+=K+(q,p)*L(r,s)",temp_tensor,uoo_tens(i),uov_tens(i))
        end do
        if(ierr.ne.EXA_SUCCESS) call quit('cholesky integrals: ooov')
        call mulliken_to_dirac_sort(temp_tensor,int_t%ooov,12)
        ierr=exatns_tensor_destroy(temp_tensor)
        if (print_level.gt.2) call print_date(' ooov done ')

!       uoo no longer needed
        do i=1,m
           ierr=exatns_tensor_destroy(uoo_tens(i))
        end do
        deallocate(uoo_tens)

!       Get antisymmetric oovv tensor
        hid=o_id
        hid(2)=v_id
        hid(4)=v_id
        hroot=nocc
        hroot(2)=nvir
        hroot(4)=nvir
        ierr=exatns_tensor_create(temp_tensor,'temp',hid,hroot,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(temp_tensor,'ZERO')
        do i=1,m
          ierr=exatns_tensor_contract("T(p,q,r,s)+=K+(q,p)*L(r,s)",temp_tensor,uvo_tens(i),uov_tens(i))
        end do
        if(ierr.ne.EXA_SUCCESS) call quit('cholesky integrals: oovv')
        call mulliken_to_dirac_sort(temp_tensor,int_t%oovv,12)
        ierr=exatns_tensor_destroy(temp_tensor)
        if (print_level.gt.2) call print_date(' oovv done ')

!       uvo no longer needed
        do i=1,m
          ierr=exatns_tensor_destroy(uvo_tens(i))
        end do
        deallocate(uvo_tens)

!       Get vovv tensor
        hid=v_id
        hid(3)=o_id
        hroot=nvir
        hroot(3)=nocc
        ierr=exatns_tensor_create(temp_tensor,'temp',hid,hroot,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(temp_tensor,ZERO)
        do i=1,m
          ierr=exatns_tensor_contract("T(p,q,r,s)+=K+(q,p)*L(r,s)",temp_tensor,uvv_tens(i),uov_tens(i))
        end do
        if(ierr.ne.EXA_SUCCESS) call quit('cholesky integrals: vovv')
        if (print_level.gt.6) call print_date(' now sorting ')
        call mulliken_to_dirac_sort(temp_tensor,int_t%vovv,34)
        ierr=exatns_tensor_destroy(temp_tensor)
        if (print_level.gt.2) call print_date(' vovv done ')

!       uov no longer needed
        do i=1,m
          ierr=exatns_tensor_destroy(uov_tens(i))
        end do
        deallocate(uov_tens)

!       Get vvvv tensor
        hid=v_id
        hroot=nvir
        ierr=exatns_tensor_create(temp_tensor,'temp',hid,hroot,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(temp_tensor,'ZERO')
        do i=1,m
          ierr=exatns_tensor_contract("T(p,q,r,s)+=K+(q,p)*L(r,s)",temp_tensor,uvv_tens(i),uvv_tens(i))
        end do
        if(ierr.ne.EXA_SUCCESS) call quit('cholesky integrals: vvvv')
        call mulliken_to_dirac_sort(temp_tensor,int_t%vvvv,12)
        ierr=exatns_tensor_destroy(temp_tensor)
        if (print_level.gt.2) call print_date(' vvvv done ')

!       uvv no longer needed
        do i=1,m
          ierr=exatns_tensor_destroy(uvv_tens(i))
        end do
        deallocate(uvv_tens)

        ierr=exatns_tensor_destroy(one_element)

      end subroutine get_CC_integrals_chol

      subroutine ao2mo_cholesky(ao_tensor,mo_tensor,mo2_dims,mo2_id,h_str,m,print_level)

!       Transform integral tensor from ao-basis to mo-basis

        implicit none

        type(tens_rcrsv_t), intent(inout)     :: ao_tensor(:), mo_tensor(:)
        integer(INTL), intent(in)             :: mo2_dims(2)
        integer(INTD), intent(in)             :: mo2_id(2)
        character(len=3), intent(in)          :: h_str
        integer, intent(in)                   :: m, print_level

        ! Auxilliary tensors that are created and destroyed inside this subroutine
        type(tens_rcrsv_t):: coef_tensor(2,2)
        type(tens_rcrsv_t):: t_tensor
        integer           :: i, ui
        integer(INTD)     :: ierr,spin,aorank
        integer(INTD)     :: ao2_id(2), t2_id(2), coef_id(2)
        integer(INTL)     :: ao2_root(2), t2_root(2), coef_root(2)
        character(len=9)  :: coef_name
        character(len=38) :: info_str
        character(len=12) :: vs_str

        if (print_level >5) call print_date('=== entered ao2mo_cholesky')

!       Determine ranks
        aorank = ao_tensor(1)%get_rank(ierr)
        if (ierr.ne.0) call quit('ao2mo_cholesky: ao tensor corrupted (rank)')

        if (aorank.eq.2) then
          call ao_tensor(1)%get_space_ids(ao2_id,ao2_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: ao tensor corrupted (2) ')

          !create temporary tensor
          t2_id(1)   = mo2_id(1)
          t2_id(2)   = ao2_id(2)
          t2_root(1) = mo2_dims(1)
          t2_root(2) = ao2_root(2)
          ierr=exatns_tensor_create(t_tensor,'temp',t2_id,t2_root,EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: temp2 tensor not created')

          !get coefficent tensors
          do i = 1, 2
            coef_id(1)   = ao2_id(i)
            coef_root(1) = ao2_root(i)
            coef_id(2)   = mo2_id(i)
            coef_root(2) = mo2_dims(i)
            do spin = 1, 2
              if (print_level>-1) write(*,*) 'coef_id:',coef_id,', coef_root:',coef_root
              write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
              ierr=exatns_tensor_create(coef_tensor(spin,i),coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_cholesky: tensor not created')
              ierr=exatns_tensor_init(coef_tensor(spin,i),mocoef_initlabel(coef_id(1),coef_id(2),spin))
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_cholesky: coef_tensor not set')
              if (print_level>-1) write(*,*) "created and initialized coef_tensor ",i,spin, &
                mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor(spin,i))
              if (print_level>15)  ierr=exatns_tensor_traverse(coef_tensor(spin,i),'PrintTensor')
            end do
          end do

          do ui = 1, m
            write (vs_str, '(A3,I9.9)') h_str,ui
            ierr=exatns_tensor_create(mo_tensor(ui),vs_str,mo2_id,mo2_dims,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_cholesky: u tensor creation failed')
            ierr=exatns_tensor_init(mo_tensor(ui),'ZERO')
            do spin = 1, 2
              ierr=exatns_tensor_init(t_tensor,'ZERO')
              i=1 ! first index
              ierr=exatns_tensor_contract(t_tensor,ao_tensor(ui),coef_tensor(spin,i),"B(p,s)+=A(r,s)*M+(r,p)")
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_cholesky: contraction 1 failed')
              i=2  ! second index
              ierr=exatns_tensor_contract(mo_tensor(ui),t_tensor,coef_tensor(spin,i),"C(p,q)+=B(p,s)*M(s,q)")
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_cholesky: contraction 2 failed')
              if (print_level>10) then
                write (info_str,'(A28,I5,A3,I2)') 'full transformation done: u=',ui,',s=',spin
                call print_date(info_str)
              end if
            end do
          end do

          do i = 1, 2
            do spin = 1, 2
              ierr=exatns_tensor_destroy(coef_tensor(spin,i))
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_cholesky: destroying coef tensor')
            end do
          end do
          ierr=exatns_tensor_destroy(t_tensor)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_cholesky: t tensor')

        else
          call quit('ao2mo_cholesky: wrong rank or mode')
        end if

        if (print_level >5) call print_date('=== leaving ao2mo_cholesky')

      end subroutine ao2mo_cholesky

      subroutine ao2mo_vect(ao_vec,mo_vec,m,ind,print_level)
!     Transform integral tensor vector from ao-basis to mo-basis

        implicit none
         
        type(tens_rcrsv_t), allocatable :: ao_vec(:), mo_vec(:)
        integer, intent(in)             :: ind, m, print_level

        type(tens_rcrsv_t)  :: coef_alpha, coef_beta
        integer(INTD)       :: ao_id(2), mo_id(2), coef_id(2)
        integer(INTL)       :: ao_dim(2), mo_dim(2), coef_dim(2)
        integer(INTD)       :: ierr, rank
        integer             :: i, ii

!       Get info about tensors
        rank = ao_vec(1)%get_rank(ierr)
        if (ierr.ne.0 .or. rank.ne.2) call quit('ao2mo_vect: ao tensor corrupted (rank)')
        call ao_vec(1)%get_space_ids(ao_id,ao_dim,ierr)
        if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_vect: ao tensor corrupted (2) ')
        call mo_vec(1)%get_space_ids(mo_id,mo_dim,ierr)
        if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_vect: mo tensor corrupted (2) ')

        if (ind.eq.11 .or. ind.eq.12) then
          ii=1
        else if (ind.eq.21 .or. ind.eq.22) then
          ii=2
        else
          call quit('ao2mo_vect: index not available ')
        end if
       
!       Initialize tensors with MO coefficients
        coef_id(1)=ao_id(ii)
        coef_id(2)=mo_id(ii)
        coef_dim(1)=ao_dim(ii)
        coef_dim(2)=mo_dim(ii)
        ierr=exatns_tensor_create(coef_alpha,'coef_alpha',coef_id,coef_dim,EXA_DATA_KIND_C8)
        if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_vect: tensor not created')
        ierr=exatns_tensor_init(coef_alpha,mocoef_initlabel(coef_id(1),coef_id(2),1))
        if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_vect: coef_tensor not set')
        if (print_level>8) write(*,*) "created and initialized alpha ", coef_id, coef_dim, &
                mocoef_initlabel(coef_id(1),coef_id(2),1),print_tensornorm2(coef_alpha)
        ierr=exatns_tensor_create(coef_beta,'coef_beta',coef_id,coef_dim,EXA_DATA_KIND_C8)
        if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_vect: tensor not created')
        ierr=exatns_tensor_init(coef_beta,mocoef_initlabel(coef_id(1),coef_id(2),2))
        if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_vect: coef_tensor not set')
        if (print_level>8) write(*,*) "created and initialized alpha ", coef_id, coef_dim, &
                mocoef_initlabel(coef_id(1),coef_id(2),2),print_tensornorm2(coef_beta)

!       transform the index
        do i=1,m
          if (ind.eq.11) then
            ierr=exatns_tensor_contract(mo_vec(i),ao_vec(i),coef_alpha,"B(p,q)+=A(i,q)*M+(i,p)")
            if(ierr.ne.EXA_SUCCESS) call quit('error in ao2mo_vect: 2a1')
          else if (ind.eq.12) then
            ierr=exatns_tensor_contract(mo_vec(i),ao_vec(i),coef_beta,"B(p,q)+=A(i,q)*M+(i,p)")
            if(ierr.ne.EXA_SUCCESS) call quit('error in ao2mo_vect: 2b1')
          else if (ind.eq.21) then
            ierr=exatns_tensor_contract(mo_vec(i),ao_vec(i),coef_alpha,"C(p,q)+=B(p,i)*M(i,q)")
            if(ierr.ne.EXA_SUCCESS) call quit('error in ao2mo_vect: 2a2')
          else if (ind.eq.22) then
            ierr=exatns_tensor_contract(mo_vec(i),ao_vec(i),coef_beta,"C(p,q)+=B(p,i)*M(i,q)")
            if(ierr.ne.EXA_SUCCESS) call quit('error in ao2mo_vect: 2b2')
          end if 
        end do

!       clean up
        ierr=exatns_tensor_destroy(coef_alpha)
        ierr=exatns_tensor_destroy(coef_beta)

      end subroutine ao2mo_vect

      subroutine decompose_cholesky_exatensor (CD, m, ao_id, ao_space, t_cholesky, &
                                                 f_delta, exa_blocksize, print_level)

          use exacorr_global
          use exacorr_tensor_methods

          type(tens_rcrsv_t), allocatable, intent(inout) :: CD(:)
          integer                                        :: m
          class(h_space_t), pointer, intent(in)          :: ao_space
          integer(INTD),intent(in)                       :: ao_id
          real(8),intent(in)                             :: t_cholesky
          type(delta_t), intent(inout)                   :: f_delta  ! fill delta matrix
          integer,intent(in)                             :: exa_blocksize, print_level

   
          integer(INTL)              :: ao_root
          integer(INTD)              :: ierr, ao_level
          integer(INTL),dimension(2) :: ao2_root
          integer(INTD),dimension(2) :: ao2_id
          integer(INTL),dimension(4) :: max_root
          integer(INTD),dimension(4) :: max_id

          integer                    :: nao, i, n_el
          type(tens_rcrsv_t)         :: Diag, Res, Factor, Dlk, Delta_lk, Delta, Temp
          type(tens_rcrsv_t)         :: Temp_l, Temp_lk, One_vec, P_l, P_k, one_element
          integer(INTL), allocatable :: ind_max(:)
          complex(8)                 :: scal_fac, D_max, D_el
          character(len=12)          :: vs_str

          if (print_level.gt.0) call print_date(" === entering decompose_cholesky_exatensor === ")

          ierr=exatns_tensor_create(one_element,"one_element",EXA_DATA_KIND_C8)
          ierr=exatns_tensor_init(one_element,ONE)

          !setting up ao information
          nao=get_nao() 
          ao_root=ao_space%get_root_id(ierr)
          ao2_root=ao_root
          ao2_id=ao_id
          max_root=ao_root
          max_id=ao_id

          !figure out the level for segments
          if (.TRUE.) then 
              if (nao/2<exa_blocksize) then
                n_el=nao/2
              else
                n_el=exa_blocksize
              end if
              do i=1,exa_blocksize
                !to do: optimze size of chunks, for now 10
                if (n_el/(2.0**(i-1))<10) then
                  ao_level=i
                  exit
                end if
              end do
          else
            ao_level=1
          end if
          if (print_level.gt.5) then
            write (*,*) ' --- ao_level =     ',ao_level  
            do i=1,nao
              write(*,*) ' i =',i,', id = ',ao_space%get_ancestor_id_of_level(int(i-1,INTL),ao_level,ierr),&
                         ', shell=',get_shell_index(i)
            end do
          end if

          allocate(CD(nao*nao))

          !Initialize
          ierr=exatns_tensor_create(Diag,"Diag",ao2_id,ao2_root,EXA_DATA_KIND_C8)        
          ierr=exatns_tensor_create(Res,"Residual",ao2_id,ao2_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_create(Delta,"Delta",ao2_id,ao2_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_create(Factor,"Factor",EXA_DATA_KIND_C8)
          ierr=exatns_tensor_create(Temp,"Temp",(/ ao_id, ao_id, ao_id /),(/ ao_root, ao_root, ao_root /),EXA_DATA_KIND_C8)
          ierr=exatns_tensor_create(One_vec,"One_vec",ao2_id(1:1),ao2_root(1:1),EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Tensor creation failed!')
          ierr=exatns_tensor_init(One_vec,ONE)
           
          ierr=exatns_tensor_init(Diag,'ZERO') 
          ierr=exatns_tensor_transform(Diag,'CholDiag')
          if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Computation od Diagonal failed!')
          if (print_level.gt.10) write(*,'(A7,E20.10)') "Diag  =", print_tensornorm2 (Diag)

          ! find maximum
          ierr=exatns_tensor_max(Diag,D_el,ind_max)
          if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Getting Maximum failed!')
          if (print_level.gt.4) write(*,'(A27,I8,I8,A9,E18.8,E18.8)') " --- maximum found (all):   ", &
                                    ind_max(1), ind_max(2), ', value: ',D_el
          D_max=D_el          
 
          m=0

          do while (abs(D_max)>t_cholesky)
            max_root(1)=ao_space%get_ancestor_id_of_level(ind_max(1)-1,ao_level,ierr)
            max_root(2)=ao_space%get_ancestor_id_of_level(ind_max(2)-1,ao_level,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Assigning shell failed!')
            if (print_level.gt.6) write(*,'(A20,I8,I8,A9,E18.8,E18.8)') " --- shell found :  ", &
                                  max_root(1), max_root(2), ', value: ',D_max

            ierr=exatns_tensor_create(Dlk,"Dlk",max_id,max_root,EXA_DATA_KIND_C8)  
            ierr=exatns_tensor_init(Dlk,'AointSegment')
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Computing Off-Diagonals failed!')
            if (print_level.gt.10) write(*,'(A7,E20.10)') "Dlk   =",print_tensornorm2 (Dlk)
            ierr=exatns_tensor_create(Delta_lk,"Delta_lk",max_id(1:3),max_root(1:3),EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Creating Delta_lk failed!')
            ierr=exatns_tensor_init(Delta_lk,'ZERO')

            ierr=exatns_tensor_create(P_l,"P_l",(/ ao_id, max_id(1) /),(/ ao_root, max_root(1) /),EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('P_l creation failed!')
            ierr=exatns_tensor_transform(P_l,'Unity')
            ierr=exatns_tensor_create(P_k,"P_k",(/ ao_id, max_id(2), ao_id /),(/ ao_root, max_root(2), ao_root /),EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('P_k creation failed!')
            ierr=exatns_tensor_transform(P_k,'Unity')
            ierr=exatns_tensor_create(Temp_l,"Temp_l",(/ max_id(1), ao_id /),(/ max_root(1), ao_root /),EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(Temp_lk,"Temp_lk",max_id(1:3),max_root(1:3),EXA_DATA_KIND_C8)

            do while (abs(D_el)>abs(D_max)/1000)

              !create current elment of the cholesky vector
              write (vs_str, '(A3,I9.9)') 'CD_',m+1
              ierr=exatns_tensor_create(CD(m+1),vs_str,ao2_id,ao2_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Creating Cholesky vector failed!')
              ierr=exatns_tensor_init(CD(m+1),'ZERO')
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Setting up Cholesky vector failed!')

              !extracting the relevant element
              ierr=exatns_tensor_init(Temp,'ZERO')
              ierr=exatns_tensor_init(Res,'ZERO')
              call f_delta%reset((/ ind_max(1), ind_max(2), int(1,INTL) /),3)
              ierr=exatns_tensor_transform(Delta_lk,f_delta)
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Constructing Delta_lk failed!')
              if (print_level.gt.10) write(*,'(A10,E20.10)') "Delta_lk =",print_tensornorm2 (Delta_lk)
              ierr=exatns_tensor_contract("T(p,q,u)+=V(r,s,p,q)*D(r,s,u)",Temp,Dlk,Delta_lk)
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Writing to temporary!')
              ierr=exatns_tensor_contract("R(p,q)+=T(p,q,u)*O(u)",Res,Temp,One_vec)
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Extracting off-diagonal Elements!')
              if (print_level.gt.10) write(*,'(A7,E20.10)') "Dpq   =",print_tensornorm2 (Res)

              !add previous cholesky elements
              do i=1,m
                call f_delta%reset((/ ind_max(1), ind_max(2) /), 2)
                ierr=exatns_tensor_transform(Delta,f_delta)
                if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Constructing Delta failed!')
                ierr=exatns_tensor_init(Factor,ZERO)
                ierr=exatns_tensor_contract("F()=L(p,q)*D(p,q)",Factor,CD(i),Delta)
                if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Extracting old value failed!')
                ierr=exatns_tensor_get_scalar(Factor,scal_fac)
                ierr=exatns_tensor_contract("R(p,q)+=L+(p,q)*F()",Res,CD(i),Factor,MINUS_ONE)
                if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Getting old contribution failed!')
              end do

              !sum and denominate
              scal_fac=ONE / zsqrt(D_el)
              ierr=exatns_tensor_contract("W(p,q)+=R+(p,q)*O()",CD(m+1),Res,one_element,scal_fac)
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Scaling failed!')
              if (print_level.gt.10) write(*,'(A7,E20.10)') "CD_m  =",print_tensornorm2 (CD(m+1))
            
              !update Diagonal
              ierr=exatns_tensor_init(Res,'ZERO')
              ierr=exatns_tensor_contract("L(p,q)+=R(p,q)*O()",Res,CD(m+1),one_element)
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Copying failed!')
              ierr=exatns_tensor_transform(Res,'Square')
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Squaring failed!')            
              ierr=exatns_tensor_contract("D(p,q)+=L(p,q)*O()",Diag,Res,one_element,MINUS_ONE)
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Updating Diagonal failed!')
              if (print_level.gt.8) write(*,'(A7,E20.10)') "Diag  =",print_tensornorm2(Diag)

              m=m+1

              if (m.ge.nao*nao) then
                write (*,*) 'WARNING: Cholesky decomposition failed, error=',D_max
                exit
              end if

              !find new maximum
              deallocate(ind_max)
              ierr=exatns_tensor_init(Temp_l,'ZERO')
              ierr=exatns_tensor_contract("L(l,q)+=D(p,q)*P(p,l)",Temp_l,Diag,P_l)
              ierr=exatns_tensor_init(Temp_lk,'ZERO')
              ierr=exatns_tensor_contract("D(l,k,u)+=L(l,q)*P(q,k,u)",Temp_lk,Temp_l,P_k)
              ierr=exatns_tensor_max(Temp_lk,D_el,ind_max)
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Getting Maximum failed!')
              if (print_level.gt.6) write(*,'(A27,I8,I8,A9,E18.8,E18.8)') " --- maximum found (shell): ", &
                                    ind_max(1), ind_max(2), ', value: ',D_el
         
              if (abs(D_el)<t_cholesky) exit

            end do

            !find new maximum
            deallocate(ind_max)
            ierr=exatns_tensor_max(Diag,D_el,ind_max)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Getting Maximum failed!')
            if (print_level.gt.4) write(*,'(A27,I8,I8,A9,E18.8,E18.8)') " --- maximum found (all):   ", &
                                    ind_max(1), ind_max(2), ', value: ',D_el
            D_max=D_el

            ierr=exatns_tensor_destroy(P_k)
            ierr=exatns_tensor_destroy(P_l)
            ierr=exatns_tensor_destroy(Temp_l)
            ierr=exatns_tensor_destroy(Temp_lk)
            ierr=exatns_tensor_destroy(Delta_lk)
            ierr=exatns_tensor_destroy(Dlk)

            if (m.ge.nao*nao) exit

          end do

          if (print_level.gt.-1) write(*,*) 'Cholesky decomposition done, used ',m,' elements'
          
          !cleanup
          deallocate(ind_max)
          ierr=exatns_tensor_destroy(Diag)
          ierr=exatns_tensor_destroy(Res)
          ierr=exatns_tensor_destroy(Delta)
          ierr=exatns_tensor_destroy(Factor)
          ierr=exatns_tensor_destroy(Temp)
          ierr=exatns_tensor_destroy(One_vec)
          ierr=exatns_tensor_destroy(one_element)

          !cholesky tensor needs to be destroyed elsewhere

          if (print_level.gt.0) call print_date(" === leaving decompose_cholesky_exatensor === ")

        end subroutine decompose_cholesky_exatensor

        subroutine decompose_cholesky_exatensor_full (CD, m, ao_id, ao_space, t_cholesky, &
                                                 j_chol, f_delta, print_level)

          use exacorr_global
          use exacorr_tensor_methods

          type(tens_rcrsv_t), allocatable, intent(inout) :: CD(:)
          integer                                        :: m
          class(h_space_t), pointer, intent(in)          :: ao_space
          integer(INTD),intent(in)                       :: ao_id
          real(8),intent(in)                             :: t_cholesky
          type(chol_J), intent(inout)                    :: j_chol  ! Cholesky: compute off-diagonal elements
          type(delta_t), intent(inout)                   :: f_delta  ! fill delta matrix
          integer,intent(in)                             :: print_level

   
          integer(INTL)              :: ao_root
          integer(INTD)              :: ierr
          integer(INTL),dimension(2) :: ao2_root
          integer(INTD),dimension(2) :: ao2_id

          integer                    :: nao, i
          type(tens_rcrsv_t)         :: Diag, Res, Delta, Factor, Dlk, one_element
          integer(INTL), allocatable :: ind_max(:)
          complex(8)                 :: scal_fac, D_max
          character(len=12)          :: vs_str

          if (print_level.gt.0) call print_date(" === entering decompose_cholesky_exatensor_full === ")

          ierr=exatns_tensor_create(one_element,"one_element",EXA_DATA_KIND_C8)
          ierr=exatns_tensor_init(one_element,ONE)

          nao  = get_nao()   ! This is the total number of ao's
          ao_root=ao_space%get_root_id(ierr) !get the root space, covers all ao's
          ao2_root=ao_root
          ao2_id=ao_id

          allocate(CD(nao*nao))

          !Initialize
          ierr=exatns_tensor_create(Diag,"Diagonal",ao2_id,ao2_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_create(Dlk,"Dlk",ao2_id,ao2_root,EXA_DATA_KIND_C8)          
          ierr=exatns_tensor_create(Res,"Residual",ao2_id,ao2_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_create(Delta,"Delta",ao2_id,ao2_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_create(Factor,"Factor",EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Tensor creation failed!')
           
          ierr=exatns_tensor_init(Diag,'ZERO') 
          ierr=exatns_tensor_transform(Diag,'CholDiag')
          if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Computation od Diagonal failed!')
          if (print_level.gt.10) write(*,'(A7,E20.10)') "Diag  =", print_tensornorm2 (Diag)

          ! find maximum
          ierr=exatns_tensor_max(Diag,D_max,ind_max)
          if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Finding Maximum failed!')
          call f_delta%reset(ind_max,2)
          ierr=exatns_tensor_init(Delta,'ZERO')
          ierr=exatns_tensor_transform(Delta,f_delta)
          if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_transform(f_delta) failed!')
          if (print_level.gt.10) write(*,'(A7,E20.10)') "Delta =",print_tensornorm2 (Delta)
          
          m=0

          do while (abs(D_max)>t_cholesky)
            if (print_level.gt.4) write(*,*) "--- maximum found: ", ind_max, D_max

            !create current elment of the cholesky vector
            write (vs_str, '(A3,I9.9)') 'CD_',m+1
            ierr=exatns_tensor_create(CD(m+1),vs_str,ao2_id,ao2_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Creating Cholesky vector failed!')
            ierr=exatns_tensor_init(CD(m+1),'ZERO')
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Setting up Cholesky vector failed!')

            call j_chol%reset(ind_max)
            ierr=exatns_tensor_init(Dlk,'ZERO')
            ierr=exatns_tensor_transform(Dlk,j_chol)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Computing Off-Diagonals failed!')
            if (print_level.gt.15) ierr=exatns_tensor_traverse(Dlk,'PrintTensor')
            if (print_level.gt.10) write(*,'(A7,E20.10)') "Dpq   =",print_tensornorm2 (Dlk)

            ierr=exatns_tensor_init(Res,'ZERO')
            do i=1,m
              ierr=exatns_tensor_init(Factor,ZERO)
              ierr=exatns_tensor_contract("F()=L(p,q)*D(p,q)",Factor,CD(i),Delta)
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Extracting old value failed!')
              ierr=exatns_tensor_get_scalar(Factor,scal_fac)
              ierr=exatns_tensor_contract("R(p,q)+=L+(p,q)*F()",Res,CD(i),Factor)
              if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Getting old contribution creation failed!')
            end do
            if (print_level.gt.10) write(*,'(A7,E20.10)') "Res   =",print_tensornorm2 (Res)

            !sum and denominate
            ierr=exatns_tensor_contract("T(p,q)+=R(p,q)*O()",Dlk,Res,one_element,MINUS_ONE)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Computing Residual failed!')
            scal_fac=ONE / zsqrt(D_max)
            ierr=exatns_tensor_contract("W(p,q)+=T+(p,q)*O()",CD(m+1),Dlk,one_element,scal_fac)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Scaling failed!')
            if (print_level.gt.10) write(*,'(A7,E20.10)') "CD_m  =",print_tensornorm2 (CD(m+1))
            
            ierr=exatns_tensor_init(Res,'ZERO')
            ierr=exatns_tensor_contract("L(p,q)+=R(p,q)*O()",Res,CD(m+1),one_element)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Copying failed!')
            ierr=exatns_tensor_transform(Res,'Square')
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Squaring failed!')            
            ierr=exatns_tensor_contract("D(p,q)+=L(p,q)*O()",Diag,Res,one_element,MINUS_ONE)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Updating Diagonal failed!')
            if (print_level.gt.8) write(*,'(A7,E20.10)') "Diag  =",print_tensornorm2(Diag)

            m=m+1

            if (m.ge.nao*nao) then
              !call quit(" Cholesky decomposition failed")
              write (*,*) 'WARNING: Cholesky decomposition failed, error=',D_max
              exit
            end if

            !find new maximum
            deallocate(ind_max)
            ierr=exatns_tensor_max(Diag,D_max,ind_max)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Getting Maximum failed!')
            call f_delta%reset(ind_max,2)
            ierr=exatns_tensor_transform(Delta,f_delta)
            if(ierr.ne.EXA_SUCCESS) call quit('decomp_cholesky: Constructing Delta failed!')
            if (print_level.gt.10) write(*,'(A7,E20.10)') "Delta =",print_tensornorm2 (Delta)
          end do

          if (print_level.gt.-1) write(*,*) 'Cholesky decomposition done, used ',m,' elements'
          
          !cleanup
          deallocate(ind_max)
          ierr=exatns_tensor_destroy(Diag)
          ierr=exatns_tensor_destroy(Dlk)
          ierr=exatns_tensor_destroy(Res)
          ierr=exatns_tensor_destroy(Delta)
          ierr=exatns_tensor_destroy(Factor)
          ierr=exatns_tensor_destroy(one_element)

          !cholesky tensor needs to be destroyed elsewhere

          if (print_level.gt.0) call print_date("=== leaving decompose_cholesky_exatensor_full === ")

        end subroutine decompose_cholesky_exatensor_full

! stub non-MPI routine to satisfy references

#else

       contains

        subroutine exacorr_compute_mo_integrals (nmo, mo_list, print_level)
         implicit none
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         integer, intent(in) :: print_level
         print*,  'called exatensor version of compute_mo_integrals in a serial version'
         print*,  'this functionality only works with MPI, we need to stop here'
         stop ' activate MPI'
        end subroutine exacorr_compute_mo_integrals

#endif

end module exacorr_ao_to_mo
