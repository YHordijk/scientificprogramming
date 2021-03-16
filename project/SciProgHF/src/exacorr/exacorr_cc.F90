module exacorr_cc

!This module contains the routines needed to calculate a coupled cluster wave function
!almost everything requires MPI, but we provide a stub function to make compilation of 
!non-MPI versions of DIRAC possible.

#if defined (VAR_MPI)

        use exatensor
        use exacorr_datatypes
        use exacorr_utils
        use intermediates

        implicit none

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0), &
                                 ONE_QUARTER=(0.25D0,0.D0), MINUS_ONE_QUARTER=(-0.25D0,0.D0),  &
                                 MINUS_TWO=(-2.D0,0.D0), TWO=(2.D0,0.D0), &
                                 THREE=(3.D0,0.D0),SIX=(6.D0,0.D0)
        type(tens_rcrsv_t),public :: one_tensor
        real(8) :: e_core
        logical :: one_el_exist

        private

        public exacorr_cc_driver

        interface diis
           module procedure diis
        end interface

       contains

!-------------------------------------------------------------------------------------------
        subroutine exacorr_cc_driver (exa_input)

!        This subroutine drives the ccsd energy calculation
!        The implementation is based on the parallel ExaTensor library

!        Written by Lucas Visscher, winter 2016/2017

         use interface_to_mpi
         use exacorr_global
         use exacorr_tensor_methods
         use exacorr_ao_to_mo
         use exacorr_gradient
         

         implicit none

         type(exacc_input), intent(in) :: exa_input

         integer                :: nocc,nvir           ! the size of the mo basis for occupied and virtual spinors

         integer(INTD) :: my_role
         integer       :: my_MPI_rank                    ! needed to start exatensor, MPI should already be started in threaded mode

         integer(INTD)                  :: ierr
         complex(8)                     :: result_val

!        Common tensors
         type(exatns_comm_tens) :: comm_t
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(exatns_intg_tens) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t) :: l1_tensor, l2_tensor
!        scalars (need to be defined as tensor types)
         type(tens_rcrsv_t)         :: result_tensor
!        space ids and roots
         type(space_dims) :: dims
!        CCSD control and result variables
         real(8)            :: scf_energy, mp2_energy, cc_energy, t_energy(3) ! triples not implemented
         real(8)            :: t1diag
         real(8)            :: t_target_precision

         ! classes for the definition of spaces (AO and 2 different MO spaces)
         class(h_space_t), pointer                :: ao_space
         class(h_space_t), pointer                :: occ_space, vir_space

         ! variables for the definition of spaces (AO and 2 different MO spaces)
         integer(INTL)              :: ao_space_root, occ_space_root, vir_space_root
         integer(INTL),dimension(2) :: oo_root, ov_root, vo_root, vv_root
         integer(INTL),dimension(4) :: oooo_root, ooov_root, oovv_root, vovo_root, &
                                       vvoo_root, vovv_root, vvvv_root
         integer(INTD),dimension(2) :: oo_id, ov_id, vo_id, vv_id 
         integer(INTD),dimension(4) :: oooo_id, ooov_id, oovv_id, vovo_id, &
                                       vvoo_id, vovv_id, vvvv_id
         integer(INTD)              :: ao_space_id, occ_space_id, vir_space_id

         integer :: ncycles

         !CCD
         logical :: CCD, lambda

         !CC2
         logical :: CC2         !L

         !methods
         type(delta_t)      :: f_delta
         type(denom_t)      :: denom
         type(set_energy_t) :: orb_set
         type(denom3_t)     :: denom3
         type(chol_J)       :: j_chol
         type(set_ff_t)     :: ff_set

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         if (my_MPI_rank == 0) call print_date(" Starting CC calculation with exatensor")

!        Initialize and register AO and MO spaces and get their IDs for further reference
         call exacorr_cc_spaces (exa_input,ao_space,occ_space,vir_space,ao_space_id,occ_space_id,vir_space_id)
         if (my_MPI_rank == 0) call print_date(" Registered AO and MO spaces")

!        Register and initialize user defined methods needed in the CC calculation

         call exacorr_cc_methods (exa_input,ao_space_id,occ_space_id,vir_space_id, &
                                  f_delta,orb_set,denom,denom3,j_chol,ff_set)
         if (my_MPI_rank == 0) call print_date(" Initialized and registered tensor methods")

!        copy input variables to local ones for ease of reference
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         t_target_precision = exa_input%t_econv
         CCD = exa_input%ccd
         CC2= exa_input%CC2     !L
         lambda = exa_input%lambda
         ncycles = exa_input%ncycles

         if (my_MPI_rank == 0) then
            call print_exainput(exa_input,int(ao_space%get_root_id(ierr)))
         end if

         !Create all id and root arrays, neccessary for exatns_tensor_create
         ao_space_root=ao_space%get_root_id(ierr)
         occ_space_root=occ_space%get_root_id(ierr)
         vir_space_root=vir_space%get_root_id(ierr)
         dims%occ_space_id=occ_space_id
         dims%vir_space_id=vir_space_id
         dims%occ_space_root=occ_space_root
         dims%vir_space_root=vir_space_root
         oo_id = (/occ_space_id, occ_space_id/)
         ov_id = (/occ_space_id, vir_space_id/)
         vo_id = (/vir_space_id, occ_space_id/)
         vv_id = (/vir_space_id, vir_space_id/)
         oo_root = (/occ_space_root, occ_space_root/)
         ov_root = (/occ_space_root, vir_space_root/)
         vo_root = (/vir_space_root, occ_space_root/)
         vv_root = (/vir_space_root, vir_space_root/)
         oooo_id = (/occ_space_id, occ_space_id, occ_space_id, occ_space_id/)
         ooov_id = (/occ_space_id, occ_space_id, occ_space_id, vir_space_id/)
         oovv_id = (/occ_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vovo_id = (/vir_space_id, occ_space_id, vir_space_id, occ_space_id/)
         vvoo_id = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
         vovv_id = (/vir_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vvvv_id = (/vir_space_id, vir_space_id, vir_space_id, vir_space_id/)
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
         ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
         oovv_root = (/occ_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vovo_root = (/vir_space_root, occ_space_root, vir_space_root, occ_space_root/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)
         vovv_root = (/vir_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vvvv_root = (/vir_space_root, vir_space_root, vir_space_root, vir_space_root/)

!        All preliminary work is now done on all MPI-nodes, we may proceed to start ExaTensor.
!        After the start of ExaTensor, the driver process (the last MPI process) will be executing
!        the commands and tell the others (that stay inside exatns_process_role) what to do.
!        These other nodes return only after exatns_stop() is called by the driver.

!        Start ExaTENSOR within MPI_COMM_WORLD (called global_communicator in the MPI interface that DIRAC uses):
         !call exatns_ctrl_reset_logging(1,2) !logging for debugging purposes
         ierr=exatns_start(int(global_communicator,INT_MPI))

         if(ierr.eq.EXA_SUCCESS) then
            ierr=exatns_process_role(my_role) ! only EXA_DRIVER will return immediately
            if(my_role.eq.EXA_DRIVER) then

!           Create and initialize scalars that are to be used as tensors in contractions
            ierr=exatns_tensor_create(result_tensor,"result_tensor",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(one_tensor,"another_one_tensor",EXA_DATA_KIND_C8)
            ! Note:  also get_integral_tensors will define a one_tensor and call it "one_tensor", so we needed another name
            ierr=exatns_tensor_init(one_tensor,ONE)
            
!           Get fixed two-body tensors
            call print_date('-Start- Integral Transformation')
            ierr=exatns_tensor_create(int_t%oooo,"oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%ooov,"ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%oovv,"oovv",oovv_id,oovv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%vovo,"vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%vovv,"vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%vvvv,"vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            if (exa_input%moint_scheme.gt.40) then
                call get_CC_integrals_chol (ao_space_id,ao_space,int_t,exa_input%moint_scheme, &
                    exa_input%exa_blocksize,exa_input%print_level,exa_input%t_cholesky,j_chol,f_delta)
            else
                call get_CC_integrals (ao_space_id,ao_space,int_t,exa_input%moint_scheme,exa_input%print_level)
            end if
            call print_date('-End- Integral Transformation')

            ! If we have the 1-body integrals we rebuild the Fock matrix and check consistency, otherwise we use the orbital energies
            ierr=exatns_tensor_create(comm_t%foo,"foo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%fov,"fov",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%fvv,"fvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            call init_fock_exat (comm_t,int_t,oo_id,oo_root,ov_id,ov_root,vv_id,vv_root, &
                            f_delta,denom,orb_set,ff_set,exa_input,scf_energy)
            call print_date('initialized Fock matrix')

!           Initialize t1 amplitudes
            if (.not.CCD) then
               ierr=exatns_tensor_create(comm_t%t1,"t1",vo_id,vo_root,EXA_DATA_KIND_C8)
               ierr=exatns_tensor_init(comm_t%t1,ZERO)
               ierr=exatns_tensor_contract("T(a,i)+=V+(i,a)*K()",comm_t%t1,comm_t%fov,one_tensor)
               ierr=exatns_tensor_transform(comm_t%t1,denom)
            end if

!           Initialize t2 amplitudes with the MP2 values
            ierr=exatns_tensor_create(comm_t%t2,"t2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            ierr=exatns_tensor_init(comm_t%t2,ZERO)
            ierr=exatns_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)*K()",comm_t%t2,int_t%oovv,one_tensor)
            ierr=exatns_tensor_transform(comm_t%t2,denom)

!           Calculate and print the MP2 energy
            ierr=exatns_tensor_init(result_tensor,ZERO)
            ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%t2,int_t%oovv,ONE_QUARTER)
            ierr=exatns_tensor_get_scalar(result_tensor,result_val)
            mp2_energy=real(result_val,8)
            ierr=exatns_tensor_init(result_tensor,ZERO)
            ierr=exatns_tensor_contract("E()+=T(a,i)*V(i,a)",result_tensor,comm_t%t1,comm_t%fov)
            ierr=exatns_tensor_get_scalar(result_tensor,result_val)
            mp2_energy=mp2_energy+real(result_val,8)
         
            write(*,*) ""
            write(*,*) "    MP2 energy = ", mp2_energy
            write(*,*) ""

!           Start CC(S)D iterations
            call print_date('-Start- CC Iterations')
            ierr=exatns_tensor_create(comm_t%goo,"goo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%gvv,"gvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%hoo,"hoo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%hvv,"hvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            if (.not.CCD) ierr=exatns_tensor_create(comm_t%hov,"hov",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%a_int,"aint",oooo_id,oooo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%h_int,"hint",vovo_id,vovo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%tau,"tau",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            if (.not.CCD) then
                call diis (comm_t%t2,vvoo_id,vvoo_root,comm_t%t1,vo_id,vo_root,mode=1)
            else
                call diis (comm_t%t2,vvoo_id,vvoo_root,mode=1)
            end if
            if (ierr.eq.EXA_SUCCESS) then
               write(*,*) " ccdriver: allocated all tensors"
            else
               write(*,*) " ccdriver: error in tensor allocation, code=",ierr
               return
            end if

            if (.not.CC2) then
               call solve_t_equations(comm_t,int_t,dims,CCD,t_target_precision,cc_energy,ncycles,denom,exa_input%print_level)
            else if (CC2) then 
               call solve_CC2_t_equations(comm_t,int_t,dims,CCD,t_target_precision,cc_energy,denom,ncycles,exa_input%print_level)
            end if
            call print_date('-End- CC Iterations')

            if (.not.CCD) then
              ! calculate T1 diagnostic
              t1diag = print_tensornorm2(comm_t%t1)
              t1diag = t1diag/sqrt(dble(nocc))
              if (exa_input%do_triples .and. .not.CC2) then
                call print_date('-Start- Triples')
                !call t3_energy_delta(t_energy,comm_t,int_t,f_delta,denom,exa_input%mo_occ,exa_input%print_level)
                call t3_energy_split(t_energy,comm_t,int_t,occ_space,denom,exa_input)
                call print_date('-End- Triples')
                call print_exaoutput(exa_input%talsh, CCD, CC2, scf_energy, mp2_energy, cc_energy, &
                                    t1diag=t1diag, t_energy=t_energy)
              else
                call print_exaoutput(exa_input%talsh, CCD, CC2, scf_energy, mp2_energy, cc_energy, &
                                    t1diag=t1diag)
              end if
            else
               call print_exaoutput(exa_input%talsh, CCD, CC2, scf_energy, mp2_energy, cc_energy)
            end if

!           Destroy tensors and shut down library
            if (.not.CCD) then
                call diis (comm_t%t2,t1_tensor=comm_t%t1,mode=2)
            else
                call diis(comm_t%t2,mode=2)
            end if

            if (lambda.and..not.CC2) then
             call print_date('-Start- Lambda Equations')
!            construct lambda tensors + initialize with t-amplitudes
             if (.not.CCD) then
                ierr=exatns_tensor_create(l1_tensor,"l1",ov_id,ov_root,EXA_DATA_KIND_C8)
                ierr=exatns_tensor_init(l1_tensor,ZERO)
                ierr=exatns_tensor_contract("L(i,a)+=T+(a,i)*K()",l1_tensor,comm_t%t1,one_tensor)
             end if

             ierr=exatns_tensor_create(l2_tensor,"l2",oovv_id,oovv_root,EXA_DATA_KIND_C8)
             if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
             ierr=exatns_tensor_init(l2_tensor,ZERO)
             ierr=exatns_tensor_contract("L(i,j,a,b)+=T+(a,b,i,j)*K()",l2_tensor,comm_t%t2,one_tensor)

!            solve lambda-equations and caluclate density matrix in MO basis
             if (.not.CCD) then
                if (.not.CC2) then
                  call solve_lambda_equations(comm_t,int_t,l2_tensor,denom,dims,t_target_precision,&
                      & ncycles,exa_input%print_level,l1_tensor)
                else
                  call solve_CC2_lambda_equations(comm_t,int_t,l2_tensor,denom,dims,t_target_precision,&
                      & ncycles,exa_input%print_level,l1_tensor)
                end if
                call exacorr_dm_driver(exa_input,dims,comm_t%t2,l2_tensor,comm_t%t1,l1_tensor)
             else
               call solve_lambda_equations(comm_t,int_t,l2_tensor,denom,dims,t_target_precision, &
                                          ncycles,exa_input%print_level)
               call exacorr_dm_driver(exa_input,dims,comm_t%t2,l2_tensor)
             end if
             
             if (.not.CCD) ierr=exatns_tensor_destroy(l1_tensor)
             ierr=exatns_tensor_destroy(l2_tensor)
             call print_date('-End- Lambda Equations')
            else
               write(*,*) 'Lambda equations with CC2 wavefunctions not supported in this release'
            end if

            ierr=exatns_tensor_destroy(int_t%oooo)
            ierr=exatns_tensor_destroy(int_t%ooov)
            ierr=exatns_tensor_destroy(int_t%oovv)
            ierr=exatns_tensor_destroy(int_t%vovo)
            ierr=exatns_tensor_destroy(int_t%vovv)
            ierr=exatns_tensor_destroy(int_t%vvvv)
            ierr=exatns_tensor_destroy(comm_t%foo)
            ierr=exatns_tensor_destroy(comm_t%fov)
            ierr=exatns_tensor_destroy(comm_t%fvv)
            ierr=exatns_tensor_destroy(comm_t%tau)
            ierr=exatns_tensor_destroy(comm_t%hoo)
            if (.not.CCD) ierr=exatns_tensor_destroy(comm_t%hov)
            ierr=exatns_tensor_destroy(comm_t%hvv)
            ierr=exatns_tensor_destroy(comm_t%goo)
            ierr=exatns_tensor_destroy(comm_t%gvv)
            ierr=exatns_tensor_destroy(comm_t%a_int)
            ierr=exatns_tensor_destroy(comm_t%h_int)
            ierr=exatns_tensor_destroy(result_tensor)
            ierr=exatns_tensor_destroy(one_tensor)
            if (.not.CCD) ierr=exatns_tensor_destroy(comm_t%t1)
            ierr=exatns_tensor_destroy(comm_t%t2)

            !Stop ExaTENSOR runtime
            ierr=exatns_stop()

            end if

         else
            write(*,*) ' Process ',my_MPI_rank,' terminated with error ',ierr
         endif

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         if (my_MPI_rank == 0) then
            write(*,*) ""
            write(*,*) " Leaving cc_driver routine"
         end if
         return

        end subroutine exacorr_cc_driver

        subroutine exacorr_cc_spaces (exa_input,ao_space,occ_space,vir_space,ao_space_id,occ_space_id,vir_space_id)

!        This subroutine defines and registers all spaces for ExaTensor

!        Written by Lucas Visscher, January 2019

         use interface_to_mpi
         use exacorr_global
         use exacorr_ao_to_mo

         implicit none

         type(exacc_input), intent(in) :: exa_input
         ! classes for the definition of spaces (AO and 2 different MO spaces)
         class(h_space_t), pointer     :: ao_space, occ_space, vir_space
         ! ids of these spaces (are not stored in the space class itself ?)
         integer(INTD), intent(out)    :: ao_space_id, occ_space_id, vir_space_id

         integer                :: nocc,nvir       ! the size of the mo basis for occupied and virtual ORBITALS
         integer(INTL)          :: nocc_l, nvir_l  ! the size of the mo basis for occupied and virtual SPINORS

         integer                                :: nao      ! the number of atomic orbitals)
         integer                                :: nshells  ! the number of shells (orbitals with the same center and l-value)
         integer                                :: basis_angular    ! 1=cartesian, 2=spherical
         type(basis_func_info_t), allocatable   :: gto(:)   ! array with information about the shells

         ! variables for the definition of spaces (AO and 2 different MO spaces)
         type(subspace_basis_t)     :: basis_ao, basis_occ, basis_vir
         integer(INTL)              :: nao_l ! (same as nao, as long integer)
         integer(INTD)              :: branch_factor, tavp_mng_depth
         integer                    :: my_MPI_rank, MPI_size

         ! variables for the definition of colors inside the spaces
         type(color_symmetry_t), allocatable    :: color(:)
         integer(INTL)                          :: labs

         ! Loop counters and error codes
         integer(INTD) :: ierr, lsh
         integer :: l, nl
         integer :: num_blocks

         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         call interface_mpi_comm_size(global_communicator,MPI_size)

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

         !Get ExaTENSOR TAVP-MNG depth:
         ierr = exatns_virtual_depth(tavp_mng_depth)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_virtual_depth() failed!')
         if(exa_input%print_level.gt.8 .and. my_MPI_rank == 0) write(*,*) 'ExaCorr: Inferred TAVP-MNG depth = ',tavp_mng_depth

         !Register the ao vector space:
         branch_factor = max(2, (int(get_num_segments(int(nao,8),int(exa_input%exa_blocksize,8)),4) - 1)/tavp_mng_depth + 1)
         ierr=exatns_space_register('AO_space',basis_ao,ao_space_id,ao_space,branch_factor)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')
         if(exa_input%print_level.gt.8 .and. my_MPI_rank == 0) write (*,'(A,I3,A,I3)') & 
           ' ExaCorr: Registered AO vector space of dimension   ',nao, &
           ': Number of segments = ', branch_factor*(2**(tavp_mng_depth-1)) !debug
         if(exa_input%print_level.gt.8 .and. my_MPI_rank == 0) call ao_space%print_it()

         !Copy user input to local arrays for the definition of MO spaces (for standard CC we only distinguish occupied and virtual spinors)
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         !Get the MO basis for the occupied and virtual spaces
         ierr = 0
         nocc_l    = int(nocc,8)  ! NB: imo, jmo and nmo count orbitals, nocc_l counts spinors
         nvir_l    = int(nvir,8)  ! NB: imo, jmo and nmo count orbitals, nvir_l counts spinors
         call basis_occ%subspace_basis_ctor(nocc_l,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')
         call basis_vir%subspace_basis_ctor(nvir_l,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

         !Set the functions inside the occupied basis (without colors, not needed for MOs)
         do labs = 1, nocc_l
            call basis_occ%set_basis_func(labs,BASIS_ABSTRACT,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.set_basis_func() failed!')
         enddo

         !Set the functions inside the virtual basis (without colors, not needed for MOs)
         do labs = 1, nvir_l
            call basis_vir%set_basis_func(labs,BASIS_ABSTRACT,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.set_basis_func() failed!')
         enddo

         !Finalize the mo space basis
         call basis_occ%finalize(ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.finalize() failed!')
         call basis_vir%finalize(ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.finalize() failed!')

         !Register occ mo vector space:
         branch_factor = max(2, (int(get_num_segments(int(nocc,8),int(exa_input%exa_blocksize,8)),4) - 1)/tavp_mng_depth + 1)
         ierr=exatns_space_register('Occupied',basis_occ,occ_space_id,occ_space,branch_factor)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register(occ) failed!')
         if(exa_input%print_level.gt.8 .and. my_MPI_rank == 0) write (*,'(A,I3,A,I3)') &
           ' ExaCorr: Registered OCC vector space of dimension  ',nocc, &
           ': Number of segments = ', branch_factor*(2**(tavp_mng_depth-1))
         num_blocks=(branch_factor*(2**(tavp_mng_depth-1)))**2

         !Register vir mo vector space:
         branch_factor = max(2, (int(get_num_segments(int(nvir,8),int(exa_input%exa_blocksize,8)),4) - 1)/tavp_mng_depth + 1)
         ierr=exatns_space_register('Virtual',basis_vir,vir_space_id,vir_space,branch_factor)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register(vir) failed!')
         if(exa_input%print_level.gt.8 .and. my_MPI_rank == 0) write (*,'(A,I3,A,I3)') &
           ' ExaCorr: Registered VIRT vector space of dimension ',nvir, &
           ': Number of segments = ', branch_factor*(2**(tavp_mng_depth-1)) !debug
         num_blocks=num_blocks*(branch_factor*(2**(tavp_mng_depth-1)))**2
         
         if(exa_input%print_level.gt.-1 .and. my_MPI_rank == 0) write (*,'(A,I9,A,I9)') &
           'number of blocks = ',num_blocks,'; number of MPI processes =',MPI_size
         
         if(num_blocks<MPI_size) then
            write (*,'(A)') 'In the current implementation there is a scaling limit.'
            write (*,'(A)') 'Using your configuration there would be nodes doing no work.'
            write (*,'(A)') 'Therefore decrease the number of MPI processes or .EXA_BLOCKSIZE'
            call quit('Not enough work for MPI processes')
         end if

        end subroutine exacorr_cc_spaces


        subroutine exacorr_cc_methods (exa_input,ao_space_id,occ_space_id,vir_space_id, &
                                       f_delta,orb_set,denom,denom3,j_chol, ff_set)

!        This subroutine initalizes and registers all user-defined methods for ExaTensor

!        Written by Lucas Visscher, January 2019

         use exacorr_global
         use exacorr_tensor_methods
         use exacorr_ao_to_mo

         implicit none

         ! container with user input (we only need the indices of the occupied and virtual MOs here)
         type(exacc_input), intent(in) :: exa_input
         ! identifiers of the spaces that we will use in the methods
         integer(INTD), intent(in)     :: ao_space_id, occ_space_id, vir_space_id
         type(delta_t)                 :: f_delta  ! fill delta matrix
         type(set_energy_t)            :: orb_set  ! place orbital energies on diagonal
         type(denom_t)                 :: denom    ! method to scale CC amplitudes by the inverse of orbital energy differences
         type(denom3_t)                :: denom3   ! scale with fixed occ value
         type(chol_J)                  :: j_chol   ! Cholesky: compute diagonal elements
         type(set_ff_t)                :: ff_set   ! method to add property integrals

         integer                :: nocc,nvir           ! the size of the mo basis for occupied and virtual orbitals (Kramers pairs)
         integer, allocatable   :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
         integer, allocatable   :: mo_oo(:),mo_vo(:),mo_ov(:),mo_vv(:) ! joined lists needed for dm init (ugly code, should be refactored)

         integer(INTD)                :: ierr, spin

         ! types of the initializer methods
         type(tens_printer_t)         :: tens_printer         ! generic method, but needs to be defined here in specific form
         type(compute_2e_ao_tensor_t) :: aoint_calculator     ! method to initialize the ao-integrals tensor
         type(compute_2e_ao_seg_t)    :: aoint_seg            ! only for segment
         type(init_mocoef_tensor_t)   :: occupied_init(2), virtual_init(2) ! methods to initialize MO coeffcient tensors
         type(init_dm_tensor_t)       :: dm_init(4)           ! array of initializers for MO transformation algorithm 2
         type(set_diagonal_t)         :: one_diag             ! identity matrix
         type(set_1el_t)              :: one_el(3)            ! one electron integrals
         type(square_t)               :: t_square             ! square matrix elements
         type(set_projection_t)       :: set_proj             ! projection matrix for triples
         type(set_zero_t)             :: cmplx_zero,imag_zero ! zero out tensors
         type(chol_diag)              :: diag_chol            ! Cholesky: compute diagonal elements

         integer:: min_occ

         !Copy user input to local arrays for the definition of MO spaces (for standard CC we only distinguish occupied and virtual spinors)
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir

         min_occ=minval(mo_occ)

         !Define separate initializers for the alpha(A) and beta(B) occupied and virtual MO-coefficients
         do spin = 1, 2
            !Register this initializer for occupied as a method
            call occupied_init(spin)%init_mocoef_ctor(mo_occ,nocc,spin)
            ierr=exatns_method_register(mocoef_initlabel(ao_space_id,occ_space_id,spin),occupied_init(spin))
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
            !Same procedure for the virtual space
            call virtual_init(spin)%init_mocoef_ctor(mo_vir,nvir,spin)
            ierr=exatns_method_register(mocoef_initlabel(ao_space_id,vir_space_id,spin),virtual_init(spin))
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
         end do

         if (exa_input%moint_scheme == 2) then
           !Define names for the initializers for density matrices
           allocate(mo_oo(nocc+nocc))
           allocate(mo_vo(nvir+nocc))
           allocate(mo_ov(nocc+nvir))
           allocate(mo_vv(nvir+nvir))
           mo_oo(1:nocc) = mo_occ
           mo_ov(1:nocc) = mo_occ
           mo_vo(1:nvir) = mo_vir
           mo_vv(1:nvir) = mo_vir
           mo_oo(nocc+1:nocc+nocc) = mo_occ
           mo_ov(nocc+1:nocc+nvir) = mo_vir
           mo_vo(nvir+1:nvir+nocc) = mo_occ
           mo_vv(nvir+1:nvir+nvir) = mo_vir
           !Register the dm initializers as methods
           call dm_init(1)%init_dm_ctor(mo_oo,(/nocc,nocc/))
           ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),(/occ_space_id,occ_space_id/)),dm_init(1))
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
           call dm_init(2)%init_dm_ctor(mo_vo,(/nvir,nocc/))
           ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),(/vir_space_id,occ_space_id/)),dm_init(2))
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
           call dm_init(3)%init_dm_ctor(mo_ov,(/nocc,nvir/))
           ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),(/occ_space_id,vir_space_id/)),dm_init(3))
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
           call dm_init(4)%init_dm_ctor(mo_vv,(/nvir,nvir/))
           ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),(/vir_space_id,vir_space_id/)),dm_init(4))
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
         end if

         !Register the ao integral generation as a method
         ierr=exatns_method_register('Aoint_calculator',aoint_calculator)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(aoint) failed!')

         !Register the ao integral generation for a segment
         ierr=exatns_method_register('AointSegment',aoint_seg)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(AointSegment) failed!')

         !Register the denominator initializer as a method
         call denom%denom_t_init(mo_occ,nocc,mo_vir,nvir,exa_input%level_shift,occ_space_id)
         ierr=exatns_method_register('Denominate',denom)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(Denominate) failed!')

         !Register denominator of triples using delta
         call denom3%denom3_t_init(mo_occ,nocc,mo_vir,nvir,exa_input%level_shift)
         ierr=exatns_method_register('Denominate3',denom3)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(Denominate3) failed!')

         !Register the standard print method with treshold such that it only prints interesting elements
         call tens_printer%reset_thresh(1.0D-6)
         ierr=exatns_method_register('PrintTensor',tens_printer)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(printer) failed!')

         !Register method to set diagonal to one
         call one_diag%set_diagonal_init(ONE,.TRUE.)
         ierr=exatns_method_register('Unity',one_diag)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(Unity) failed!')

         !Register method to set projection for triples
         call set_proj%set_projection_init(ONE)
         ierr=exatns_method_register('SetProj',set_proj)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(SetProj) failed!')

         one_el_exist = exist_one_el(e_core)
         if (one_el_exist) then
            !Register methods to read one electron integrals
            call one_el(1)%set_1el_init(mo_occ,nocc,mo_occ,nocc,min_occ)
            ierr=exatns_method_register('OE_oo',one_el(1))
            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(OE_oo) failed!')
            call one_el(2)%set_1el_init(mo_occ,nocc,mo_vir,nvir,min_occ)
            ierr=exatns_method_register('OE_ov',one_el(2))
            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(OE_ov) failed!')
            call one_el(3)%set_1el_init(mo_vir,nvir,mo_vir,nvir,min_occ)
            ierr=exatns_method_register('OE_vv',one_el(3))
            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(OE_vv) failed!')

            !Register methods to make matrix with orbital energies on the diagonal
            call orb_set%set_energy_init(mo_occ,nocc,occ_space_id,mo_vir,nvir,vir_space_id,exa_input%level_shift)
            ierr=exatns_method_register('OrbSet',orb_set)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(OrbSet) failed!')
         end if

         !Register method to add property integrals
         call ff_set%set_ff_init(mo_occ,nocc,occ_space_id,mo_vir,nvir,vir_space_id,min_occ)
         ierr=exatns_method_register('FF_set',ff_set)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(FF_set) failed!')

         !Register method to set delta tensor (up to 4 dims)
         ierr=exatns_method_register('Delta_F',f_delta)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(Delta_F) failed!')
         
         !Register method to square tensor
         ierr=exatns_method_register('Square',t_square)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(Square) failed!')
         
         !Register method to set to ZERO
         call cmplx_zero%set_zero_init(3)
         ierr=exatns_method_register('ZERO',cmplx_zero)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(cmplx_zero) failed!')
         
         !Register method to set imaginary part to zero
         call imag_zero%set_zero_init(1)
         ierr=exatns_method_register('ZERO_I',imag_zero)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(imag_zero) failed!')

         !Register computation of diagonals for cholesky
         ierr=exatns_method_register('CholDiag',diag_chol)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(CholDiag) failed!')

         !Register computation of off-diagonals for cholesky
         ierr=exatns_method_register('Chol_J',j_chol)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(Chol_J) failed!')

        end subroutine exacorr_cc_methods

        subroutine solve_t_equations(comm_t,int_t,dims,CCD,t_target_precision,ccsd_energy,ncycles,denom,print_level)

!        Routine for solving the t amplitude equations

         use exacorr_ao_to_mo, only : print_tensornorm2
         use exacorr_tensor_methods

!        Input
         type(exatns_comm_tens), intent(inout) :: comm_t
         type(exatns_intg_tens), intent(inout) :: int_t
         type(space_dims),intent(in) :: dims
         integer(INTD) :: occ_space_id, vir_space_id
         integer(INTL) :: occ_space_root, vir_space_root
         logical, intent(in)  :: CCD
         real(8), intent(in)  :: t_target_precision
         real(8), intent(out) :: ccsd_energy
         type(denom_t)        :: denom
         integer, intent(in)  :: print_level
!        solution tensors
         type(tens_rcrsv_t) :: s1_tensor, s2_tensor
!        intermediate
         type(tens_rcrsv_t) :: tau2_tensor, cint_tensor
!        auxiliary tensors
         type(tens_rcrsv_t) :: hov_aux_tensor, hoo_aux_tensor
!        result tensor
         type(tens_rcrsv_t) :: result_t_tensor
         complex(8)         :: result_val
!        Ids and roots for tensors
         integer(INTD),dimension(2) :: oo_id,ov_id,vo_id
         integer(INTD),dimension(4) :: ovoo_id,vvoo_id
         integer(INTL),dimension(2) :: oo_root,ov_root,vo_root
         integer(INTL),dimension(4) :: ovoo_root,vvoo_root

         real(8) :: t_convergence

         integer :: iccsd
         integer,intent(in) :: ncycles
         integer(INTD) :: ierr

         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root
         if (.not.CCD) then
            oo_id     = (/occ_space_id, occ_space_id/)
            ov_id     = (/occ_space_id, vir_space_id/)
            vo_id     = (/vir_space_id, occ_space_id/)
            oo_root   = (/occ_space_root, occ_space_root/)
            ov_root   = (/occ_space_root, vir_space_root/)
            vo_root   = (/vir_space_root, occ_space_root/)
         end if
         ovoo_id   = (/occ_space_id, vir_space_id, occ_space_id, occ_space_id/)
         vvoo_id   = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
         ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)

!        Initialize scalar that are to be used as tensors in contractions
         ierr=exatns_tensor_create(result_t_tensor,"result_t",EXA_DATA_KIND_C8)

!        Construct intermediate tensors
         if (.not.CCD) then
            ierr=exatns_tensor_create(hoo_aux_tensor,"hoo_aux",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(hov_aux_tensor,"hov_aux",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(s1_tensor,"s1",vo_id,vo_root,EXA_DATA_KIND_C8)
         end if
         ierr=exatns_tensor_create(s2_tensor,"s2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(tau2_tensor,"tau2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(cint_tensor,"cint",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')

         do iccsd = 1, ncycles

          call get_tau(comm_t,CCD,one_tensor)

          ierr=exatns_tensor_init(s2_tensor,'ZERO')

!         The vvvv contribution
          ierr=exatns_tensor_contract("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,int_t%vvvv,comm_t%tau,ONE_HALF) !o(2)v(4)

!         vv
          call get_hvv(comm_t,int_t%oovv,one_tensor)
          call get_gvv(comm_t,int_t%vovv,CCD,one_tensor)

!         vovo
          ierr=exatns_tensor_init(tau2_tensor,'ZERO')
          ierr=exatns_tensor_contract("T(a,d,l,i)+=R(a,d,l,i)*K()",tau2_tensor,comm_t%t2,one_tensor,ONE_HALF)
          if (.not.CCD) ierr=exatns_tensor_contract("T(a,d,l,i)+=T(a,l)*T(d,i)",tau2_tensor,comm_t%t1,comm_t%t1) ! o(2)v(2)
          call get_h(comm_t,int_t,CCD,one_tensor)
          ierr=exatns_tensor_contract("H(a,k,c,i)+=V(k,l,c,d)*T(a,d,l,i)",comm_t%h_int,int_t%oovv,tau2_tensor)! o(3)v(3)

!         This term contributes with four permutations (we need to permute both i and j and a and b), re-use tau2 to store unpermuted result
          ierr=exatns_tensor_init(tau2_tensor,'ZERO')
          ierr=exatns_tensor_contract("R(a,b,i,j)+=H(a,k,c,i)*T(c,b,j,k)",tau2_tensor,comm_t%h_int,comm_t%t2)  ! unpermuted result
          ierr=exatns_tensor_contract("S(a,b,i,j)+=R(a,b,i,j)*K()",s2_tensor,tau2_tensor,one_tensor)           ! add to S
          ierr=exatns_tensor_contract("S(a,b,i,j)+=R(a,b,j,i)*K()",s2_tensor,tau2_tensor,one_tensor,MINUS_ONE) ! add with Pij anti-symmetrization
          ierr=exatns_tensor_contract("S(a,b,i,j)+=R(b,a,i,j)*K()",s2_tensor,tau2_tensor,one_tensor,MINUS_ONE) ! add with Pab anti-symmetrization
          ierr=exatns_tensor_contract("S(a,b,i,j)+=R(b,a,j,i)*K()",s2_tensor,tau2_tensor,one_tensor)           ! add with PijPab anti-symmetrization

!         ovoo

          ierr=exatns_tensor_init(cint_tensor,'ZERO')
          ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,d)*T(c,d,i,j)",cint_tensor,int_t%vovv,comm_t%tau,MINUS_ONE_HALF) ! o(3)v(3)
          if (.not.CCD) then
              ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,i)*T(c,j)",cint_tensor,int_t%vovo,comm_t%t1) ! o(3)v(2)
              ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,j)*T(c,i)",cint_tensor,int_t%vovo,comm_t%t1,MINUS_ONE) ! o(3)v(2)
          end if
          ierr=exatns_tensor_contract("C(k,b,i,j)+=V+(i,j,k,b)*K()",cint_tensor,int_t%ooov,one_tensor) ! o(3)v(1)

!         oo
          call get_hoo(comm_t,int_t%oovv,one_tensor)
          call get_goo(comm_t,int_t%ooov,CCD,one_tensor)

!         ov
          if (.not.CCD) then
            call get_hov(comm_t,int_t%oovv,one_tensor)
            ierr=exatns_tensor_init(hov_aux_tensor,'ZERO')
            ierr=exatns_tensor_contract("L(k,c)+=H(k,c)*K()",hov_aux_tensor,comm_t%hov,one_tensor)        ! copy common term
            ierr=exatns_tensor_contract("H(k,c)+=F(k,c)*K()",hov_aux_tensor,comm_t%fov,one_tensor,MINUS_TWO) ! subtract fock matrix
          end if

!         oooo
          call get_a(comm_t,int_t,CCD,one_tensor)

          if (ierr.ne.EXA_SUCCESS) then
             write(*,*) " ccdriver: error in calculation of common blocks, code=",ierr
             return
          end if

          if (print_level.gt.8) then
             if (.not.CCD) write(*,'(A7,F20.15)') "T1amp = ",print_tensornorm2 (comm_t%t1)
             write(*,'(A7,F20.15)') "T2amp = ",print_tensornorm2 (comm_t%t2)
             write(*,'(A7,F20.15)') "Tau   = ",print_tensornorm2 (comm_t%tau)
             write(*,'(A7,F20.15)') "hoint = ",print_tensornorm2 (comm_t%hoo)
             write(*,'(A7,F20.15)') "goint = ",print_tensornorm2 (comm_t%goo)
             write(*,'(A7,F20.15)') "hvint = ",print_tensornorm2 (comm_t%hvv)
             write(*,'(A7,F20.15)') "gvint = ",print_tensornorm2 (comm_t%gvv)
             if (.not.CCD) write(*,'(A7,F20.15)') "hovint= ",print_tensornorm2 (comm_t%hov)
             write(*,'(A7,E20.15)') "aint  = ",print_tensornorm2 (comm_t%a_int)
             write(*,'(A7,F20.15)') "cint  = ",print_tensornorm2 (cint_tensor)
             write(*,'(A7,F20.15)') "hint  = ",print_tensornorm2 (comm_t%h_int)
          end if

!         This term appears twice as we need to permute i and j
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(k,j)*T(a,b,i,k)",s2_tensor,comm_t%goo,comm_t%t2)                   ! o(3)v(2)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(k,i)*T(a,b,k,j)",s2_tensor,comm_t%goo,comm_t%t2)                   ! o(3)v(2)
!         Again two times, also note the complex conjugation as we use vovv instead of vvvo
          if (.not.CCD) then
              ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(c,j,a,b)*T(c,i)",s2_tensor,int_t%vovv,comm_t%t1)                   ! o(2)v(3)
              ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(c,i,a,b)*T(c,j)",s2_tensor,int_t%vovv,comm_t%t1,MINUS_ONE)         ! o(2)v(3)
          end if

!         This term appears twice as we need to permute a and b
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,comm_t%gvv,comm_t%t2)  ! o(2)v(3)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(b,c)*T(a,c,i,j)",s2_tensor,comm_t%gvv,comm_t%t2)  ! o(2)v(3)
!         Again two times, with the C intermediate (note: in RELCCSD this is done differently)
          if (.not.CCD) then
              ierr=exatns_tensor_contract("S(a,b,i,j)+=C(k,b,i,j)*T(a,k)",s2_tensor,cint_tensor,comm_t%t1,MINUS_ONE) ! o(3)v(2)
              ierr=exatns_tensor_contract("S(a,b,i,j)+=C(k,a,i,j)*T(b,k)",s2_tensor,cint_tensor,comm_t%t1)           ! o(3)v(2)
          end if

!         The oooo contribution
          ierr=exatns_tensor_contract("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,comm_t%a_int,comm_t%tau,ONE_HALF) !o(4)v(2)

          ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)    ! o(2)v(2)

!         calculate s1 tensor
          if (.not.CCD) then
            ierr=exatns_tensor_init(s1_tensor,'ZERO')
            ierr=exatns_tensor_contract("S(a,i)+=F+(i,a)*K()",s1_tensor,comm_t%fov,one_tensor)
            ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,int_t%vovv,comm_t%tau,ONE_HALF) !o(2)v(3)
            ierr=exatns_tensor_contract("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,comm_t%hvv,comm_t%t1)                   !o(1)v(2)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,comm_t%hoo,comm_t%t1)                   !o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,comm_t%hov,comm_t%t2)               !o(2)v(2)
            ierr=exatns_tensor_init(hoo_aux_tensor,'ZERO')
            ierr=exatns_tensor_contract("L(k,i)+=H(k,c)*T(c,i)",hoo_aux_tensor,hov_aux_tensor,comm_t%t1)          !o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoo_aux_tensor,comm_t%t1)               !o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,int_t%ooov,comm_t%tau,ONE_HALF) !o(3)v(2)
            ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,int_t%vovo,comm_t%t1,MINUS_ONE)     !o(2)v(2)
          end if

          if (ierr.ne.EXA_SUCCESS) then
             write(*,*) " ccdriver: error in evaluation of amplitude equations, code=",ierr
             return
          end if

!         prepare for next iteration
          if (.not.CCD) then
              ierr=exatns_tensor_init(comm_t%t1,'ZERO')
              ierr=exatns_tensor_contract("T1(a,i)+=S1(a,i)*K()",comm_t%t1,s1_tensor,one_tensor)
              ierr=exatns_tensor_transform(comm_t%t1,denom)
          end if
          ierr=exatns_tensor_init(comm_t%t2,'ZERO')
          ierr=exatns_tensor_contract("T2(a,b,i,j)+=S2(a,b,i,j)*K()",comm_t%t2,s2_tensor,one_tensor)
          ierr=exatns_tensor_transform(comm_t%t2,denom)

          if (.not.CCD) then
              call diis (comm_t%t2,vvoo_id,vvoo_root,comm_t%t1,vo_id,vo_root,t_convergence)
          else
              call diis (comm_t%t2,vvoo_id,vvoo_root,t_convergence=t_convergence)
          end if

          ierr=exatns_tensor_init(result_t_tensor)
          if (.not.CCD) then
!             Will always be 0 contribution for now, since fov=0 (canonical orbitals)
              ierr=exatns_tensor_contract("E()+=T(a,i)*F(i,a)",result_t_tensor,comm_t%t1,comm_t%fov) ! o(1)v(1)
              ierr=exatns_tensor_get_scalar(result_t_tensor,result_val)
              ccsd_energy=real(result_val,8)
              ierr=exatns_tensor_init(result_t_tensor)
          else
              ccsd_energy=real(ZERO,8)
          end if

!         Update tau with converged values
          call get_tau(comm_t,CCD,one_tensor)

          ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_t_tensor,comm_t%tau,int_t%oovv) ! o(2)v(2)
          ierr=exatns_tensor_get_scalar(result_t_tensor,result_val)
          ccsd_energy = ccsd_energy + real(result_val*ONE_QUARTER,8)

          call print_iteration(iccsd,t_convergence,ccsd_energy,print_level)
          if (t_convergence .lt. t_target_precision) exit
         end do

         write(*,*) "-----------------------------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         write(*,*) ""

         if (.not.CCD) ierr=exatns_tensor_destroy(s1_tensor)
         ierr=exatns_tensor_destroy(s2_tensor)
         ierr=exatns_tensor_destroy(tau2_tensor)
         ierr=exatns_tensor_destroy(cint_tensor)
         if (.not.CCD) then
            ierr=exatns_tensor_destroy(hoo_aux_tensor)
            ierr=exatns_tensor_destroy(hov_aux_tensor)
         end if
         ierr=exatns_tensor_destroy(result_t_tensor)

        end subroutine solve_t_equations

        subroutine solve_CC2_t_equations(comm_t,int_t,dims,CCD,t_target_precision,ccsd_energy,denom,ncycles,print_level)

!        Routine for solving the t amplitude equations

         use exacorr_ao_to_mo, only : print_tensornorm2
         use exacorr_tensor_methods

         logical :: debug_local = .false.

!        Input
         type(exatns_comm_tens), intent(inout) :: comm_t
         type(exatns_intg_tens), intent(inout) :: int_t
         type(space_dims),intent(in) :: dims
         integer(INTD) :: occ_space_id, vir_space_id
         integer(INTL) :: occ_space_root, vir_space_root
         logical, intent(in) :: CCD
         real(8), intent(in) :: t_target_precision
         real(8), intent(out) :: ccsd_energy
         type(denom_t)        :: denom
!        solution tensors
         type(tens_rcrsv_t) :: s1_tensor, s2_tensor

!        intermediary tensor for handling permutation symmetry
         type(tens_rcrsv_t) :: intm_p2_tensor, intm_vvvv_tensor

! Add tensor usefull to avoid tauCC2_tensor 
         type(tens_rcrsv_t) :: Hovvv_CC2_tensor,Hvovv_CC2_tensor
         type(tens_rcrsv_t) :: Hvvov_CC2_tensor,Hvvvo_CC2_tensor
         type(tens_rcrsv_t) :: Hooov_CC2_tensor,Hoovo_CC2_tensor
         type(tens_rcrsv_t) :: Hvooo_CC2_tensor,Hovoo_CC2_tensor

        
! Intermediate for old tauCC2 method
         type(tens_rcrsv_t) :: t1t1_tensor, tauCC2_tensor

         type(tens_rcrsv_t) :: cint_vvvo_tensor, cint_ovoo_tensor
         type(tens_rcrsv_t) :: Wabef_CC2_tensor, Wmnij_CC2_tensor
         type(tens_rcrsv_t) :: Wabej_tensor, Wmbij_tensor
!        auxiliary tensors
         type(tens_rcrsv_t) :: hov_aux_tensor, hoo_aux_tensor
!        result tensor
         type(tens_rcrsv_t) :: result_t_tensor
         complex(8)         :: result_val
!        Ids and roots for tensors
         integer(INTD),dimension(2) :: oo_id,ov_id,vo_id
         integer(INTD),dimension(4) :: ovoo_id,vvoo_id,vvvo_id
         integer(INTD),dimension(4) :: vvvv_id, oooo_id
         integer(INTL),dimension(2) :: oo_root,ov_root,vo_root
         integer(INTL),dimension(4) :: ovoo_root,vvoo_root,vvvo_root
         integer(INTL),dimension(4) :: vvvv_root, oooo_root

! Add new dims for Hxxxx 
         integer(INTD),dimension(4) :: ovvv_id , vovv_id
         integer(INTD),dimension(4) :: vvov_id 
         integer(INTD),dimension(4) :: oovo_id , ooov_id 
         integer(INTD),dimension(4) :: vooo_id

         integer(INTL),dimension(4) :: ovvv_root , vovv_root
         integer(INTL),dimension(4) :: vvov_root 
         integer(INTL),dimension(4) :: oovo_root , ooov_root
         integer(INTL),dimension(4) :: vooo_root

         real(8) :: t_convergence

         integer :: iccsd
         integer,intent(in) :: ncycles,print_level
         integer(INTD) :: ierr

         if (print_level.gt.8) debug_local = .true.

         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root
         if (.not.CCD) then
            oo_id     = (/occ_space_id, occ_space_id/)
            ov_id     = (/occ_space_id, vir_space_id/)
            vo_id     = (/vir_space_id, occ_space_id/)
            oo_root   = (/occ_space_root, occ_space_root/)
            ov_root   = (/occ_space_root, vir_space_root/)
            vo_root   = (/vir_space_root, occ_space_root/)
         end if
         ovoo_id   = (/occ_space_id, vir_space_id, occ_space_id, occ_space_id/)
         vvoo_id   = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
         vvvo_id   = (/vir_space_id, vir_space_id, vir_space_id, occ_space_id/)
         oooo_id   = (/occ_space_id, occ_space_id, occ_space_id, occ_space_id/)
         vvvv_id   = (/vir_space_id, vir_space_id, vir_space_id, vir_space_id/)


         ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)
         vvvo_root = (/vir_space_root, vir_space_root, vir_space_root, occ_space_root/)
         vvvv_root = (/vir_space_root, vir_space_root, vir_space_root, vir_space_root/)
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)

! New dims for Hxxxx
         ovvv_id   = (/occ_space_id,vir_space_id,vir_space_id,vir_space_id/)
         vovv_id   = (/vir_space_id,occ_space_id,vir_space_id,vir_space_id/)
         vvov_id   = (/vir_space_id,vir_space_id,occ_space_id,vir_space_id/)
         oovo_id   = (/occ_space_id,occ_space_id,vir_space_id,occ_space_id/)
         ooov_id   = (/occ_space_id,occ_space_id,occ_space_id,vir_space_id/)
         vooo_id   = (/vir_space_id,occ_space_id,occ_space_id,occ_space_id/)
         ovvv_root   = (/occ_space_root,vir_space_root,vir_space_root,vir_space_root/)
         vovv_root   = (/vir_space_root,occ_space_root,vir_space_root,vir_space_root/)
         vvov_root   = (/vir_space_root,vir_space_root,occ_space_root,vir_space_root/)
         oovo_root   = (/occ_space_root,occ_space_root,vir_space_root,occ_space_root/)
         ooov_root   = (/occ_space_root,occ_space_root,occ_space_root,vir_space_root/)
         vooo_root   = (/vir_space_root,occ_space_root,occ_space_root,occ_space_root/)

         ierr=exatns_tensor_create(Wabej_tensor,"Wabej",vvvo_id,vvvo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(Wmbij_tensor,"Wmbij",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)

!        Initialize scalar that are to be used as tensors in contractions
         ierr=exatns_tensor_create(result_t_tensor,"result_t",EXA_DATA_KIND_C8)

!        Construct intermediate tensors
         if (.not.CCD) then
            ierr=exatns_tensor_create(hoo_aux_tensor,"hoo_aux",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(hov_aux_tensor,"hov_aux",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(s1_tensor,"s1",vo_id,vo_root,EXA_DATA_KIND_C8)
         end if
         ierr=exatns_tensor_create(s2_tensor,"s2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(intm_p2_tensor,"cc2_intm_p2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)

         ierr=exatns_tensor_create(cint_vvvo_tensor,"cint_vvvo",vvvo_id,vvvo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(cint_ovoo_tensor,"cint_ovoo",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(intm_vvvv_tensor,"cc2_intm_vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(Wabef_CC2_tensor,"Wabef_CC2",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(Wmnij_CC2_tensor,"Wmnij_CC2",oooo_id,oooo_root,EXA_DATA_KIND_C8)

! Create intermediaites tensor 
         ierr=exatns_tensor_create(Hovvv_CC2_tensor,"Hovvv_CC2",ovvv_id,ovvv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(Hvovv_CC2_tensor,"Hvovv_CC2",vovv_id,vovv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(Hvvov_CC2_tensor,"Hvvov_CC2",vvov_id,vvov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(Hvvvo_CC2_tensor,"Hvvvo_CC2",vvvo_id,vvvo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(Hooov_CC2_tensor,"Hooov_CC2",ooov_id,ooov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(Hoovo_CC2_tensor,"Hoovo_CC2",oovo_id,oovo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(Hvooo_CC2_tensor,"Hvooo_CC2",vooo_id,vooo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(Hovoo_CC2_tensor,"Hovoo_CC2",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)

! resorted integrals nedded below

         ierr=exatns_tensor_init(cint_vvvo_tensor,'ZERO')
         ierr=exatns_tensor_contract("C(b,c,a,i)+=V+(a,i,b,c)*K()",cint_vvvo_tensor,int_t%vovv,one_tensor)

         ierr=exatns_tensor_init(cint_ovoo_tensor,'ZERO')
         ierr=exatns_tensor_contract("C(k,a,i,j)+=V+(i,j,k,a)*K()",cint_ovoo_tensor,int_t%ooov,one_tensor)

         do iccsd = 1, ncycles

          ierr=exatns_tensor_init(s2_tensor,'ZERO')
          if (debug_local) print*, "norm s2 tensor at the start of CC2 t2",print_tensornorm2(s2_tensor)

          if (iccsd==1) then
                 call get_tau(comm_t,CCD,one_tensor)
                 if (debug_local) print *, " call get_tau "
          end if

          ierr=exatns_tensor_init(s2_tensor,'ZERO')
          if (debug_local) print*, "norm s2 tensor at the start of CC2 t2",print_tensornorm2(s2_tensor)

!         creation of tau like intermediates for CC2

!         create \ts ai \ts bj in intermediate storage
          ierr=exatns_tensor_init(intm_p2_tensor,'ZERO')
          ierr=exatns_tensor_contract("I(a,b,i,j)+=T(a,i)*T(b,j)",intm_p2_tensor,comm_t%t1,comm_t%t1) 

      !         Diagramme 15 \Sd abij = \V abij   ! Inversion 
         ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)
         if (debug_local) print*, "norm s2 tensor after diagram 15 ",print_tensornorm2(s2_tensor)

         ierr=exatns_tensor_init(Wabej_tensor,'ZERO')
         ierr=exatns_tensor_init(Wmbij_tensor,'ZERO')
         ierr=exatns_tensor_contract("W(a,b,e,j)+=C(a,b,e,j)*K()",Wabej_tensor,cint_vvvo_tensor,one_tensor,ONE)   !v3o

         ierr=exatns_tensor_contract("W(a,b,e,j)+=T(a,m)*V(b,m,e,j)",Wabej_tensor,comm_t%t1,int_t%vovo,ONE_HALF)  !v3o2
         ierr=exatns_tensor_contract("W(a,b,e,j)+=T(b,m)*V(a,m,e,j)",Wabej_tensor,comm_t%t1,int_t%vovo,MINUS_ONE_HALF) !v3o2
         if (debug_local) print *, "Wabej : norm " ,print_tensornorm2(Wabej_tensor)

         ierr=exatns_tensor_contract("S(a,b,i,j)=T(e,i)*W(a,b,e,j)",s2_tensor,comm_t%t1,Wabej_tensor,ONE)  !v3o2
         ierr=exatns_tensor_contract("S(a,b,i,j)=T(e,j)*W(a,b,e,i)",s2_tensor,comm_t%t1,Wabej_tensor,MINUS_ONE) !v3o2

         ierr=exatns_tensor_contract("W(m,b,i,j)+=C(m,b,i,j)*K()",Wmbij_tensor,cint_ovoo_tensor,one_tensor,ONE)  !vo3

         ierr=exatns_tensor_contract("W(m,b,i,j)+=T(e,i)*V(b,m,e,j)",Wmbij_tensor,comm_t%t1,int_t%vovo,MINUS_ONE_HALF) !v2o3
         ierr=exatns_tensor_contract("W(m,b,i,j)+=T(e,j)*V(b,m,e,i)",Wmbij_tensor,comm_t%t1,int_t%vovo,ONE_HALF) !v2o3
         if (debug_local) print *, "Wmbij : norm ",print_tensornorm2(Wmbij_tensor)

         ierr=exatns_tensor_contract("S(a,b,i,j)=T(a,m)*W(m,b,i,j)",s2_tensor,comm_t%t1,Wmbij_tensor,MINUS_ONE) !v2o3
         ierr=exatns_tensor_contract("S(a,b,i,j)=T(b,m)*W(m,a,i,j)",s2_tensor,comm_t%t1,Wmbij_tensor,ONE) !v2o3

!         additional CC2 W intermediates
!
!         a.\W mnij = \V mnij +\Pm ij \ts ej \V mnie +\unquart \Pm ij \ts ei \ts fj \V mnef
          ierr=exatns_tensor_init(Wmnij_CC2_tensor,'ZERO')
          ierr=exatns_tensor_contract("Wmnij_CC2(m,n,i,j)+=V(m,n,i,j)*K()",Wmnij_CC2_tensor,int_t%oooo,one_tensor)
          ierr=exatns_tensor_contract("Wmnij_CC2(m,n,i,j)+=V(m,n,i,e)*T(e,j)",Wmnij_CC2_tensor,int_t%ooov,comm_t%t1)
          ierr=exatns_tensor_contract("Wmnij_CC2(m,n,i,j)+=V(m,n,j,e)*T(e,i)",Wmnij_CC2_tensor,int_t%ooov,comm_t%t1,MINUS_ONE)

          ierr=exatns_tensor_init(Hoovo_CC2_tensor,'ZERO')
          ierr=exatns_tensor_contract("Hoovo_CC2(m,n,e,j)+=V(m,n,e,f)*T(f,j)",&
              &       Hoovo_CC2_tensor,int_t%oovv,comm_t%t1)
          ierr=exatns_tensor_contract("Wmnij_CC2(m,n,i,j)+=Hoovo_CC2(m,n,e,j)*T(e,i)",&
              &       Wmnij_CC2_tensor,Hoovo_CC2_tensor,comm_t%t1,ONE_QUARTER)
          ierr=exatns_tensor_init(Hooov_CC2_tensor,'ZERO')
          ierr=exatns_tensor_contract("Hooov_CC2(m,n,j,f)+=V(m,n,e,f)*T(e,j)",&
              &       Hooov_CC2_tensor,int_t%oovv,comm_t%t1)
          ierr=exatns_tensor_contract("Wmnij_CC2(m,n,i,j)+=Hooov_CC2(m,n,j,f)*T(f,i)",&
              &       Wmnij_CC2_tensor,Hooov_CC2_tensor,comm_t%t1,MINUS_ONE_QUARTER)
          if (debug_local) print*, "norm WmnijCC2 new method",print_tensornorm2(Wmnij_CC2_tensor)

!         b.\W abef= \V abef -\Pm ab \ts bm \V amef +\unquart \Pm ab \ts am \ts bn \Vmnef

          ierr=exatns_tensor_init(intm_vvvv_tensor,'ZERO')
          ierr=exatns_tensor_contract("I(a,b,e,f)+=V(a,m,e,f)*T(b,m)",intm_vvvv_tensor,int_t%vovv,comm_t%t1)

          ierr=exatns_tensor_init(Wabef_CC2_tensor,'ZERO')
          ierr=exatns_tensor_contract("Wabef_CC2(a,b,e,f)+=V(a,b,e,f)*K()",Wabef_CC2_tensor,int_t%vvvv,one_tensor)
          ierr=exatns_tensor_contract("Wabef_CC2(a,b,e,f)+=I(a,b,e,f)*K()",Wabef_CC2_tensor,intm_vvvv_tensor,one_tensor,MINUS_ONE)
          ierr=exatns_tensor_contract("Wabef_CC2(a,b,e,f)+=I(b,a,e,f)*K()",Wabef_CC2_tensor,intm_vvvv_tensor,one_tensor)

          ierr=exatns_tensor_init(Hvovv_CC2_tensor,'ZERO')
          ierr=exatns_tensor_init(Hovvv_CC2_tensor,'ZERO')
          ierr=exatns_tensor_contract("Hvovv_CC2(a,n,e,f)+=V(m,n,e,f)*T(a,m)",&
                & Hvovv_CC2_tensor,int_t%oovv,comm_t%t1,ONE)        !Operation cost : v^3*o^2
          ierr=exatns_tensor_contract("Wabef_CC2(a,b,e,f)+=Hvovv_CC2(a,n,e,f)*T(b,n)",&
                & Wabef_CC2_tensor,Hvovv_CC2_tensor,comm_t%t1,ONE_QUARTER)!Operation cost : v^4*o
          ierr=exatns_tensor_contract("Hovvv_CC2(m,a,e,f)+=V(m,n,e,f)*T(a,n)",&
                & Hovvv_CC2_tensor,int_t%oovv,comm_t%t1,ONE)        !Operation cost :v^3*o^2
          ierr=exatns_tensor_contract("Wabef_CC2(a,b,e,f)+=Hovvv_CC2(m,a,e,f)*T(b,m)",&
                & Wabef_CC2_tensor,Hovvv_CC2_tensor,comm_t%t1,MINUS_ONE_QUARTER) ! Operation cost : v^4*o
          if (debug_local) print*, "norm Wabef",print_tensornorm2(Wabef_CC2_tensor)

!         Diagrams using the CC2 W intermediates above
!
!         Diagram 21XX : \Sd abij = \undemi \Pm ab \ts am \ts bn \W mnij

           ierr=exatns_tensor_init(Hvooo_CC2_tensor,'ZERO')
           ierr=exatns_tensor_contract("Hvooo_CC2(a,n,i,j)+=Wmnij_CC2(m,n,i,j)*T(a,m)",&
               & Hvooo_CC2_tensor,Wmnij_CC2_tensor,comm_t%t1,ONE)
           ierr=exatns_tensor_contract("S(a,b,i,j)+=Hvooo_CC2(a,n,i,j)*T(b,n)",&
               & s2_tensor,Hvooo_CC2_tensor,comm_t%t1,ONE_HALF)
           ierr=exatns_tensor_init(Hovoo_CC2_tensor,'ZERO')
           ierr=exatns_tensor_contract("Hovoo_CC2(m,a,i,j)+=Wmnij_CC2(m,n,i,j)*T(a,n)",&
               & Hovoo_CC2_tensor,Wmnij_CC2_tensor,comm_t%t1,ONE)
           ierr=exatns_tensor_contract("S(a,b,i,j)+=Hovoo_CC2(m,a,i,j)*T(b,m)",&
               & s2_tensor,Hovoo_CC2_tensor,comm_t%t1,MINUS_ONE_HALF)
           if (debug_local) print*, "norm s2 tensor after diagram 21",print_tensornorm2(s2_tensor)

!         Diagram 22XX : \Sd abij = \undemi \Pm ij \ts ei \ts fj \W abef 

          ierr=exatns_tensor_init(Hvvov_CC2_tensor,'ZERO')
          ierr=exatns_tensor_init(Hvvvo_CC2_tensor,'ZERO')
          ierr=exatns_tensor_contract("Hvvov_CC2(a,b,i,f)+=Wabef_CC2(a,b,e,f)*T(e,i)",&
                & Hvvov_CC2_tensor,Wabef_CC2_tensor,comm_t%t1,ONE)    ! Operationcost : v^4*o
          ierr=exatns_tensor_contract("S(a,b,i,j)+= Hvvov_CC2(a,b,i,f)*T(f,j)",&
                & s2_tensor,Hvvov_CC2_tensor,comm_t%t1,ONE_HALF)    !Operation cost : v^3*o^2
          ierr=exatns_tensor_contract("Hvvvo_CC2(a,b,e,i)+=Wabef_CC2(a,b,e,f)*T(f,i)",&
                & Hvvvo_CC2_tensor,Wabef_CC2_tensor,comm_t%t1,ONE)  !Operation cost : v^4*o
          ierr=exatns_tensor_contract("S(a,b,i,j)+= Hvvvo_CC2(a,b,e,i)*T(e,j)",&
                & s2_tensor,Hvvvo_CC2_tensor,comm_t%t1,MINUS_ONE_HALF)     !  Operation cost : v^3*o^2
          if (debug_local)  print*, "norm s2 tensor after diagram 22",print_tensornorm2(s2_tensor)

! Always check that the Fock matrix is called before, in comm_t 

!         Diagram 23 : s2_{ij}^{ab} +=  P\_(ab)\sum_{e} f_e^b t_{ij}^{ae} ! Use of the Fock matrix 
          ierr=exatns_tensor_init(intm_p2_tensor,'ZERO')
          ierr=exatns_tensor_contract("I(a,b,i,j)+=F(b,e)*T(a,e,i,j)",intm_p2_tensor,comm_t%fvv,comm_t%t2)

          ierr=exatns_tensor_contract("S(a,b,i,j)+=I(a,b,i,j)*K()",s2_tensor,intm_p2_tensor,one_tensor)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=I(b,a,i,j)*K()",s2_tensor,intm_p2_tensor,one_tensor,MINUS_ONE)
          if (debug_local) print*, "norm s2 tensor after diagram 23",print_tensornorm2(s2_tensor)

!         Diagram 24 : s2_{ij}^{ab} += -P\_(ij)\sum_{k} f_j^k t_{ik}^{ab} ! Use of the Fock matrix
          ierr=exatns_tensor_init(intm_p2_tensor,'ZERO')
          ierr=exatns_tensor_contract("I(a,b,i,j)+=F(m,i)*T(a,b,j,m)",intm_p2_tensor,comm_t%foo,comm_t%t2)

          ierr=exatns_tensor_contract("S(a,b,i,j)+=I(a,b,j,i)*K()",s2_tensor,intm_p2_tensor,one_tensor,MINUS_ONE)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=I(a,b,i,j)*K()",s2_tensor,intm_p2_tensor,one_tensor)
          if (debug_local) print*, "norm s2 tensor after diagram 24",print_tensornorm2(s2_tensor)

!         vv
          call get_hvv(comm_t,int_t%oovv,one_tensor)
          call get_gvv(comm_t,int_t%vovv,CCD,one_tensor)

!         vovo
          call get_h(comm_t,int_t,CCD,one_tensor)

!         oooo
          call get_a(comm_t,int_t,CCD,one_tensor)

!         oo
          call get_hoo(comm_t,int_t%oovv,one_tensor)
          call get_goo(comm_t,int_t%ooov,CCD,one_tensor)

!         ov
          if (.not.CCD) then
            call get_hov(comm_t,int_t%oovv,one_tensor)
            ierr=exatns_tensor_init(hov_aux_tensor,'ZERO')
            ierr=exatns_tensor_contract("L(k,c)+=H(k,c)*K()",hov_aux_tensor,comm_t%hov,one_tensor)        ! copy common term
            ierr=exatns_tensor_contract("H(k,c)+=F(k,c)*K()",hov_aux_tensor,comm_t%fov,one_tensor,MINUS_TWO) ! subtract fock matrix
          end if

          if (ierr.ne.EXA_SUCCESS) then
             print*," ccdriver: error in calculation of common blocks, code=",ierr
             return
          end if

          if (print_level.gt.8) then
             if (.not.CCD) write(*,'(A7,E20.15)') "T1amp = ",print_tensornorm2 (comm_t%t1)
             write(*,'(A7,E20.15)') "T2amp = ",print_tensornorm2 (comm_t%t2)
             write(*,'(A7,E20.15)') "Tau   = ",print_tensornorm2 (comm_t%tau)
             write(*,'(A7,E20.15)') "hoint = ",print_tensornorm2 (comm_t%hoo)
             write(*,'(A7,E20.15)') "goint = ",print_tensornorm2 (comm_t%goo)
             write(*,'(A7,E20.15)') "hvint = ",print_tensornorm2 (comm_t%hvv)
             write(*,'(A7,E20.15)') "gvint = ",print_tensornorm2 (comm_t%gvv)
             if (.not.CCD) write(*,'(A7,E20.15)') "hovint= ",print_tensornorm2 (comm_t%hov)
             write(*,'(A7,E20.15)') "aint  = ",print_tensornorm2 (comm_t%a_int)
             write(*,'(A7,E20.15)') "cint_ovoo  = ",print_tensornorm2 (cint_ovoo_tensor)
             write(*,'(A7,E20.15)') "cint_vvvo  = ",print_tensornorm2 (cint_vvvo_tensor)
             write(*,'(A7,E20.15)') "hint  = ",print_tensornorm2 (comm_t%h_int)
          end if


!         calculate s1 tensor
          if (.not.CCD) then
            ierr=exatns_tensor_init(s1_tensor,'ZERO')
            ierr=exatns_tensor_contract("S(a,i)+=F+(i,a)*K()",s1_tensor,comm_t%fov,one_tensor)
            ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,int_t%vovv,comm_t%tau,ONE_HALF) !o(2)v(3)
            ierr=exatns_tensor_contract("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,comm_t%hvv,comm_t%t1)                   !o(1)v(2)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,comm_t%hoo,comm_t%t1)                   !o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,comm_t%hov,comm_t%t2)               !o(2)v(2)
            ierr=exatns_tensor_init(hoo_aux_tensor,'ZERO')
            ierr=exatns_tensor_contract("L(k,i)+=H(k,c)*T(c,i)",hoo_aux_tensor,hov_aux_tensor,comm_t%t1)          !o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoo_aux_tensor,comm_t%t1)               !o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,int_t%ooov,comm_t%tau,ONE_HALF) !o(3)v(2)
            ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,int_t%vovo,comm_t%t1,MINUS_ONE)     !o(2)v(2)
          end if

          if (ierr.ne.EXA_SUCCESS) then
             print*," ccdriver: error in evaluation of amplitude equations, code=",ierr
             return
          end if

!         prepare for next iteration
          if (.not.CCD) then
              ierr=exatns_tensor_init(comm_t%t1,'ZERO')
              ierr=exatns_tensor_contract("T1(a,i)+=S1(a,i)*K()",comm_t%t1,s1_tensor,one_tensor)
              ierr=exatns_tensor_transform(comm_t%t1,denom)
          end if
          ierr=exatns_tensor_init(comm_t%t2,'ZERO')
          ierr=exatns_tensor_contract("T2(a,b,i,j)+=S2(a,b,i,j)*K()",comm_t%t2,s2_tensor,one_tensor)
          ierr=exatns_tensor_transform(comm_t%t2,denom)

          if (debug_local) then
             if (.not.CCD) write(*,'(A7,F20.15)') "T1amp = ",print_tensornorm2 (comm_t%t1)
             write(*,'(A7,F20.15)') "T1amp = ",print_tensornorm2 (comm_t%t2)
          end if


          if (.not.CCD) then
              call diis (comm_t%t2,vvoo_id,vvoo_root,comm_t%t1,vo_id,vo_root,t_convergence)
          else
              call diis (comm_t%t2,vvoo_id,vvoo_root,t_convergence=t_convergence)
          end if
          if (debug_local) then
             if (.not.CCD) write(*,'(A7,F20.15)') "T1amp = ",print_tensornorm2 (comm_t%t1)
             write(*,'(A7,F20.15)') "T1amp = ",print_tensornorm2 (comm_t%t2)
             write(*,'(A7,F20.15)') "Fov = ",print_tensornorm2 (comm_t%fov)
          end if


          ierr=exatns_tensor_init(result_t_tensor)
          if (.not.CCD) then
!             Will always be 0 contribution for now, since fov=0 (canonical orbitals)
              ierr=exatns_tensor_contract("E()+=T(a,i)*F(i,a)",result_t_tensor,comm_t%t1,comm_t%fov) ! o(1)v(1)
              ierr=exatns_tensor_get_scalar(result_t_tensor,result_val)
              ccsd_energy=real(result_val,8)
              ierr=exatns_tensor_init(result_t_tensor)
          else
              ccsd_energy=real(ZERO,8)
          end if

!         Update tau with converged values
          call get_tau(comm_t,CCD,one_tensor)

          ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_t_tensor,comm_t%tau,int_t%oovv) ! o(2)v(2)
          ierr=exatns_tensor_get_scalar(result_t_tensor,result_val)
          ccsd_energy = ccsd_energy + real(result_val*ONE_QUARTER,8)

          call print_iteration(iccsd,t_convergence,ccsd_energy)
          if (t_convergence .lt. t_target_precision) exit
         end do

         write(*,*) "-----------------------------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) print*, "WARNING: Non-converged amplitudes!"
         write(*,*) ""


         if (.not.CCD) ierr=exatns_tensor_destroy(s1_tensor)
         ierr=exatns_tensor_destroy(s2_tensor)

         ierr=exatns_tensor_destroy(cint_ovoo_tensor)
         ierr=exatns_tensor_destroy(cint_vvvo_tensor)
         ierr=exatns_tensor_destroy(Wmnij_CC2_tensor)
         ierr=exatns_tensor_destroy(Wabef_CC2_tensor)
         ierr=exatns_tensor_destroy(Wmbij_tensor)
         ierr=exatns_tensor_destroy(Wabej_tensor)

         ierr=exatns_tensor_destroy(intm_p2_tensor)
                if (debug_local) print*, "Destroy h",ierr
         ierr=exatns_tensor_destroy(intm_vvvv_tensor)

         ierr=exatns_tensor_destroy(Hooov_CC2_tensor)
         ierr=exatns_tensor_destroy(Hoovo_CC2_tensor)
         ierr=exatns_tensor_destroy(Hovoo_CC2_tensor)
         ierr=exatns_tensor_destroy(Hvooo_CC2_tensor)
         ierr=exatns_tensor_destroy(Hvvov_CC2_tensor)
         ierr=exatns_tensor_destroy(Hvvvo_CC2_tensor)
         ierr=exatns_tensor_destroy(Hovvv_CC2_tensor)
         ierr=exatns_tensor_destroy(Hvovv_CC2_tensor)

         if (.not.CCD) then
            ierr=exatns_tensor_destroy(hoo_aux_tensor)
            ierr=exatns_tensor_destroy(hov_aux_tensor)
         end if
         ierr=exatns_tensor_destroy(result_t_tensor)

        end subroutine solve_CC2_t_equations

        subroutine solve_lambda_equations(comm_t,int_t,l2_tensor,denom,dims,t_target_precision,ncycles,print_level,l1_tensor)

!        Routine for solving the lambda equations

         use exacorr_ao_to_mo, only : print_tensornorm2
         use exacorr_tensor_methods

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(exatns_comm_tens), intent(inout) :: comm_t
         type(exatns_intg_tens), intent(inout) :: int_t
         type(exatns_lambda_tens)              :: lambda_t
!        denominator
         type(denom_t)      :: denom
!        solution tensors
         type(tens_rcrsv_t), intent(inout)           :: l2_tensor
         type(tens_rcrsv_t), intent(inout), optional :: l1_tensor
         type(tens_rcrsv_t) :: s1_tensor, s2_tensor
!        Tensor dimensions
         type(space_dims),intent(in) :: dims
         integer(INTD) :: occ_space_id, vir_space_id
         integer(INTL) :: occ_space_root, vir_space_root
         integer(INTD), dimension(2) :: oo_id, ov_id, vv_id
         integer(INTD), dimension(4) :: oooo_id,ovoo_id,ooov_id,vovo_id,oovv_id, &
                                  vvvo_id,vovv_id,vvvv_id
         integer(INTL), dimension(2) :: oo_root, ov_root, vv_root
         integer(INTL), dimension(4) :: oooo_root,ovoo_root,ooov_root,vovo_root,oovv_root, &
                                        vvvo_root,vovv_root,vvvv_root
         integer :: ilambda
         integer,intent(in) :: ncycles
         integer,intent(in) :: print_level
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence
         logical :: CCD=.false.
         integer(INTD) :: ierr

         write(*,*) ""
         write(*,'(A)') "*******************************************************************************"
         write(*,'(A)') "********************** Solving Left Eigenvector Equation **********************"
         write(*,'(A)') "*******************************************************************************"
         write(*,*) ""

         if (.not.present(l1_tensor)) CCD = .true.

!**********************************
!        create all ids and roots *
!**********************************

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
         oovv_id = (/occ_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vovo_id = (/vir_space_id, occ_space_id, vir_space_id, occ_space_id/)
         vvvo_id = (/vir_space_id, vir_space_id, vir_space_id, occ_space_id/)
         vovv_id = (/vir_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vvvv_id = (/vir_space_id, vir_space_id, vir_space_id, vir_space_id/)
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
         ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
         ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
         oovv_root = (/occ_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vovo_root = (/vir_space_root, occ_space_root, vir_space_root, occ_space_root/)
         vvvo_root = (/vir_space_root, vir_space_root, vir_space_root, occ_space_root/)
         vovv_root = (/vir_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vvvv_root = (/vir_space_root, vir_space_root, vir_space_root, vir_space_root/)

!***************************************
!        construct fixed intermediates *
!***************************************

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
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')

!******************************************
!        obtain intermediates from CC(S)D *
!******************************************

         call get_hvv(comm_t,int_t%oovv,one_tensor)
         call get_gvv(comm_t,int_t%vovv,CCD,one_tensor)
         call get_h(comm_t,int_t,CCD,one_tensor)
         call get_hoo(comm_t,int_t%oovv,one_tensor)
         call get_goo(comm_t,int_t%ooov,CCD,one_tensor)
         call get_hov(comm_t,int_t%oovv,one_tensor)
         call get_a(comm_t,int_t,CCD,one_tensor)
         call print_date("Finished calculating common blocks")

!*********************************
!        get fixed intermediates *
!*********************************

         call get_w_vvvv(lambda_t%w_vvvv,comm_t,int_t,CCD,one_tensor)
         call get_wbar_vovo(lambda_t%wbar_vovo,comm_t%t2,int_t,one_tensor)
         call get_fbar_ov(lambda_t%fbar_ov,comm_t%hov,one_tensor)
         call get_w_vvvo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_vovv(lambda_t%w_vovv,comm_t%t1,int_t,CCD,one_tensor)
         call get_w_vovo(lambda_t%w_vovo,comm_t,int_t%oovv,occ_space_id,occ_space_root, &
                          vir_space_id,vir_space_root,CCD,one_tensor)
         call get_w_oooo(lambda_t%w_oooo,comm_t%a_int,one_tensor)
         call get_w_ovoo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_ooov(lambda_t%w_ooov,comm_t%t1,int_t,CCD,one_tensor)
         call get_fbar_vv(lambda_t%fbar_vv,comm_t%gvv,one_tensor)
         call get_fbar_oo(lambda_t%fbar_oo,comm_t,CCD,one_tensor)
         call print_date('Finished calculating fixed intermediates')

         if (print_level > 0) then
           if (.not.CCD) print*, "T1amp      = ", print_tensornorm2(comm_t%t1)
           print*, "T2amp      = ", print_tensornorm2(comm_t%t2)*0.25D0
           ! print*, "Tau        = ", print_tensornorm2(comm_t%tau)*0.25D0
           ! print*, "Fbar_oo    = ", print_tensornorm2(lambda_t%fbar_oo)*0.5D0
           print*, "Fbar_ov    = ", print_tensornorm2(lambda_t%fbar_ov)
           ! print*, "Fbar_vv    = ", print_tensornorm2(lambda_t%fbar_vv)*0.5D0
           ! print*, "w_vvvv     = ", print_tensornorm2(lambda_t%w_vvvv)*0.25D0
           print*, "wbar_vovo  = ", print_tensornorm2(lambda_t%wbar_vovo)
           print*, "w_vvvo     = ", print_tensornorm2(lambda_t%w_vvvo)*0.5D0
           print*, "w_vovv     = ", print_tensornorm2(lambda_t%w_vovv)*0.5D0
           print*, "w_vovo     = ", print_tensornorm2(lambda_t%w_vovo)
           print*, "w_ovoo     = ", print_tensornorm2(lambda_t%w_ovoo)*0.5D0
           print*, "w_ooov     = ", print_tensornorm2(lambda_t%w_ooov)*0.5D0
           print*, "w_oooo     = ", print_tensornorm2(lambda_t%w_oooo)*0.25D0
         end if

!************************************
!        construct solution tensors *
!************************************

         if (.not.CCD) then
            ierr=exatns_tensor_create(s1_tensor,"s1",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(s1_tensor,'ZERO')
         end if 
         ierr=exatns_tensor_create(s2_tensor,"s2",oovv_id,oovv_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(s2_tensor,'ZERO')

!********************************************
!        construct mixed type intermediates *
!********************************************

         ierr=exatns_tensor_create(lambda_t%goo,"goo_mixed",oo_id,oo_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_create(lambda_t%gvv,"gvv_mixed",vv_id,vv_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(lambda_t%goo,'ZERO')
         ierr=exatns_tensor_init(lambda_t%gvv,'ZERO')

!*****************************************************************
!        re-initialize diis, so it can be used in lambda routine *
!*****************************************************************

         if (.not.CCD) then
           call diis (l2_tensor,oovv_id,oovv_root,l1_tensor,ov_id,ov_root,mode=1)
         else
           call diis (l2_tensor,oovv_id,oovv_root,mode=1)
         end if

!*****************
!        iterate *
!*****************

         do ilambda = 1, ncycles

!***********************************
!         form mixed intermediates *
!***********************************

          call get_gvv_lambda(lambda_t%gvv,l2_tensor,comm_t%t2)
          call get_goo_lambda(lambda_t%goo,l2_tensor,comm_t%t2)

!         debug
          if (print_level.gt.8) then
            if (.not.CCD) print*, "L1  = ", print_tensornorm2(l1_tensor)
            print*, "L2  = ", print_tensornorm2(l2_tensor)*0.25D0
            print*, "goo = ", print_tensornorm2(lambda_t%goo)
            print*, "gvv = ", print_tensornorm2(lambda_t%gvv)
          end if

!******************************
!         calculate s1 tensor *
!******************************

          if (.not.CCD) then
            ierr=exatns_tensor_init(s1_tensor,'ZERO')

!----------------------------------------
!           term 1: S(i,a) += Fbar(i,a) |
!----------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=F(i,a)*K()",s1_tensor,lambda_t%fbar_ov,one_tensor)
!-------------------------------------------------
!           term 2: S(i,a) += L(i,e) * Fbar(e,a) |
!-------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,l1_tensor,lambda_t%fbar_vv)

!-------------------------------------------------
!           term 3: S(i,a) -= L(m,a) * Fbar(i,m) |
!-------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,l1_tensor,lambda_t%fbar_oo,MINUS_ONE)

!------------------------------------------------
!           term 6: S(i,a) -= G(n,m) * W(mi,na) |
!------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,lambda_t%goo,lambda_t%w_ooov)

!------------------------------------------------
!           term 7: S(i,a) += L(m,e) * W(ie,am) |
!------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,l1_tensor,lambda_t%w_vovo,MINUS_ONE)

!---------------------------------------------------------
!           term 8: S(i,a) -= 1/2 * L(mn,ae)  * W(ie,mn) |
!---------------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,l2_tensor,lambda_t%w_ovoo,ONE_HALF)

!------------------------------------------------------
!           term 4: S(i,a) += 1/2 L(im,ef) * W(ef,am) |
!------------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,l2_tensor,lambda_t%w_vvvo,ONE_HALF)

!------------------------------------------------
!           term 5: S(i,a) -= G(f,e) * W(ei,fa) |
!------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,lambda_t%gvv,lambda_t%w_vovv)
          end if

!******************************
!         calculate s2 tensor *
!******************************

          ierr=exatns_tensor_init(s2_tensor,'ZERO')

!---------------------------------------
!         term 1: S(ij,ab) += V(ij,ab) |
!---------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)

!----------------------------------------------------------
!         term 2: S(ij,ab) += P(a,b) L(ij,ae) * Fbar(e,b) |
!----------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,l2_tensor,lambda_t%fbar_vv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,l2_tensor,lambda_t%fbar_vv)

!----------------------------------------------------------
!         term 3: S(ij,ab) -= P(i,j) L(im,ab) * Fbar(j,m) |
!----------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)

!------------------------------------------------------
!         term 4: S(ij,ab) += 1/2 L(mn,ab) * W(ij,mn) |
!------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,l2_tensor,lambda_t%w_oooo,ONE_HALF)

!----------------------------------------------------------------
!         term 5: S(ij,ab) += P(i,j) P(a,b) L(im,ae) * W(je,bm) |
!----------------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)

!-------------------------------------------------------
!         term 6: S(ij,ab) += P(a,b) V(ij,ae) * G(e,b) |
!-------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",s2_tensor,int_t%oovv,lambda_t%gvv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",s2_tensor,int_t%oovv,lambda_t%gvv)

!-------------------------------------------------------
!         term 7: S(ij,ab) -= P(a,b) L(m,a) * W(,j,mb) |
!-------------------------------------------------------

          if (.not.CCD) then
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,l1_tensor,lambda_t%w_ooov)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",s2_tensor,l1_tensor,lambda_t%w_ooov)
          end if

!-------------------------------------------------------
!         term 8: S(ij,ab) -= P(i,j) V(im,ab) * G(j,m) |
!-------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,m,b,a)*G(j,m)",s2_tensor,int_t%oovv,lambda_t%goo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",s2_tensor,int_t%oovv,lambda_t%goo)

!----------------------------------------------------------
!         term 9+10: S(ij,ab) += P(i,j) L(i,e) * W(ej,ab) |
!----------------------------------------------------------

          if (.not.CCD) then
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",s2_tensor,l1_tensor,lambda_t%w_vovv)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",s2_tensor,l1_tensor,lambda_t%w_vovv)
          end if

!----------------------------------------------------------------
!         term 11: S(ij,ab) += P(i,j) P(a,b) L(i,a) * Fbar(j,b) |
!----------------------------------------------------------------

          if (.not.CCD) then
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov,MINUS_ONE)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov,MINUS_ONE)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
          end if

!-----------------------------------------------------
!         term 12: S(ij,ab) += 1/2 L(ij,ef) W(ef,ab) |
!-----------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,l2_tensor,lambda_t%w_vvvv,ONE_HALF)

          if (ierr.ne.EXA_SUCCESS) then
            write(*,*) " ccdriver: error in evaluation of lambda equations, code=",ierr
            return
          end if

!*************************************
!         prepare for next iteration *
!*************************************

          if (.not.CCD) then
            ierr=exatns_tensor_init(l1_tensor,'ZERO')
            ierr=exatns_tensor_contract("L1(i,a)+=S1(i,a)*K()",l1_tensor,s1_tensor,one_tensor)
            ierr=exatns_tensor_transform(l1_tensor,denom)
          end if
          ierr=exatns_tensor_init(l2_tensor,'ZERO')
          ierr=exatns_tensor_contract("L2(i,j,a,b)+=S2(i,j,a,b)*K()",l2_tensor,s2_tensor,one_tensor)
          ierr=exatns_tensor_transform(l2_tensor,denom)

          if (.not.CCD) then
            call diis (l2_tensor,t1_tensor=l1_tensor,t_convergence=t_convergence)
          else
            call diis (l2_tensor,t_convergence=t_convergence)
          end if

          call print_iteration(ilambda,t_convergence)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         if (.not.CCD) write(*,'(A,ES22.16)') "  Final L1amp = ", print_tensornorm2(l1_tensor)
         write(*,'(A,ES22.16)') "  Final L2amp = ", print_tensornorm2(l2_tensor)*0.25D0
         write(*,*) ""

!*****************
!        cleanup *
!*****************

         if (.not.CCD) then
           ierr=exatns_tensor_destroy(s1_tensor)
           call diis(l2_tensor,t1_tensor=l1_tensor,mode=2)
         else
           call diis(l2_tensor,mode=2)
         end if
         ierr=exatns_tensor_destroy(s2_tensor)
         ierr=exatns_tensor_destroy(one_tensor)
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

        end subroutine solve_lambda_equations

        subroutine solve_CC2_lambda_equations(comm_t,int_t,l2_tensor,denom,dims,t_target_precision,ncycles,print_level,l1_tensor)

!        Routine for solving the lambda equations

         use exacorr_ao_to_mo, only : print_tensornorm2
         use exacorr_tensor_methods

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(exatns_comm_tens), intent(inout) :: comm_t
         type(exatns_comm_tens) :: comm_t_CC2
         type(exatns_intg_tens), intent(inout) :: int_t
         type(exatns_lambda_tens)              :: lambda_t, lambda_t_CC2
!         denominator
         type(denom_t)      :: denom
!        solution tensors
         type(tens_rcrsv_t), intent(inout)           :: l2_tensor
         type(tens_rcrsv_t), intent(inout), optional :: l1_tensor
         type(tens_rcrsv_t) :: s1_tensor, s2_tensor
!        Tensor dimensions
         type(space_dims),intent(in) :: dims
         integer(INTD) :: occ_space_id, vir_space_id
         integer(INTL) :: occ_space_root, vir_space_root
         integer(INTD), dimension(2) :: oo_id, ov_id, vv_id, vo_id
         integer(INTD), dimension(4) :: oooo_id,ovoo_id,ooov_id,vovo_id,oovv_id, &
                                  vvvo_id,vovv_id,vvvv_id,vvoo_id
         integer(INTL), dimension(2) :: oo_root, ov_root, vv_root, vo_root
         integer(INTL), dimension(4) :: oooo_root,ovoo_root,ooov_root,vovo_root,oovv_root, &
                                        vvvo_root,vovv_root,vvvv_root,vvoo_root
         integer :: ilambda
         integer,intent(in) :: ncycles
         integer,intent(in) :: print_level
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence
         logical :: CCD=.false.

         logical :: print_CC2_inter=.true.  !Check value

         integer(INTD) :: ierr

         write(*,*) ""
         write(*,'(A)')"*******************************************************************************"
         write(*,'(A)') "********************** Solving Left Eigenvector Equation **********************"
         write(*,'(A)')"*******************************************************************************"
         write(*,*) ""

         if (.not.present(l1_tensor)) CCD = .true.

!**********************************
!        create all ids and roots *
!**********************************

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
         oovv_id = (/occ_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vovo_id = (/vir_space_id, occ_space_id, vir_space_id, occ_space_id/)
         vvvo_id = (/vir_space_id, vir_space_id, vir_space_id, occ_space_id/)
         vovv_id = (/vir_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vvvv_id = (/vir_space_id, vir_space_id, vir_space_id, vir_space_id/)
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
         ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
         ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
         oovv_root = (/occ_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vovo_root = (/vir_space_root, occ_space_root, vir_space_root, occ_space_root/)
         vvvo_root = (/vir_space_root, vir_space_root, vir_space_root, occ_space_root/)
         vovv_root = (/vir_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vvvv_root = (/vir_space_root, vir_space_root, vir_space_root, vir_space_root/)

!***************************************
!        construct fixed intermediates *
!***************************************
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
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')

!***
! List of NEW INTERMEDIATES for the CC2 Method
!***
         ierr=exatns_tensor_create(lambda_t_CC2%fbar_oo,"cc2_fbar_oo",oo_id,oo_root,EXA_DATA_KIND_C8)  
         ierr=exatns_tensor_create(lambda_t_CC2%fbar_vv,"cc2_fbar_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t_CC2%w_oooo,"cc2_w_oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8) 
!!!!!!  Inportant here :  The heart of w_vovo intermediate : h_int is the same
!as in    CCSD and for \Lambda_1 either : 
         ierr=exatns_tensor_create(lambda_t_CC2%w_vovo,"cc2_w_vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)! Usefull in CC2 - Same part as h_int
!!!!!!
         ierr=exatns_tensor_create(lambda_t_CC2%w_vvvv,"cc2_w_vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
!******************************************
!        obtain intermediates from CC(S)D *
!******************************************

         call get_hvv(comm_t,int_t%oovv,one_tensor)
         call get_gvv(comm_t,int_t%vovv,CCD,one_tensor)
         call get_h(comm_t,int_t,CCD,one_tensor)
         call get_hoo(comm_t,int_t%oovv,one_tensor)
         call get_goo(comm_t,int_t%ooov,CCD,one_tensor)
         call get_hov(comm_t,int_t%oovv,one_tensor)
         call get_a(comm_t,int_t,CCD,one_tensor)
         call print_date("Finished calculating common blocks")

!*********************************
!        get fixed intermediates *
!*********************************

         call get_w_vvvv(lambda_t%w_vvvv,comm_t,int_t,CCD,one_tensor)
         call get_wbar_vovo(lambda_t%wbar_vovo,comm_t%t2,int_t,one_tensor)
         call get_fbar_ov(lambda_t%fbar_ov,comm_t%hov,one_tensor)
         call get_w_vvvo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_vovv(lambda_t%w_vovv,comm_t%t1,int_t,CCD,one_tensor)
         call get_w_vovo(lambda_t%w_vovo,comm_t,int_t%oovv,occ_space_id,occ_space_root, &
                          vir_space_id,vir_space_root,CCD,one_tensor)
         call get_w_oooo(lambda_t%w_oooo,comm_t%a_int,one_tensor)
         call get_w_ovoo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_ooov(lambda_t%w_ooov,comm_t%t1,int_t,CCD,one_tensor)
         call get_fbar_vv(lambda_t%fbar_vv,comm_t%gvv,one_tensor)
         call get_fbar_oo(lambda_t%fbar_oo,comm_t,CCD,one_tensor)
         call print_date('Finished calculating fixed intermediates')

         if (print_level > 0) then
           if (.not.CCD) print*, "T1amp      = ", print_tensornorm2(comm_t%t1)
           print*, "T2amp      = ", print_tensornorm2(comm_t%t2)*0.25D0
           print*, "Fbar_ov    = ", print_tensornorm2(lambda_t%fbar_ov)
           print*, "wbar_vovo  = ", print_tensornorm2(lambda_t%wbar_vovo)
           print*, "w_vvvo     = ", print_tensornorm2(lambda_t%w_vvvo)*0.5D0
           print*, "w_vovv     = ", print_tensornorm2(lambda_t%w_vovv)*0.5D0
           print*, "w_vovo     = ", print_tensornorm2(lambda_t%w_vovo)
           print*, "w_ovoo     = ", print_tensornorm2(lambda_t%w_ovoo)*0.5D0
           print*, "w_ooov     = ", print_tensornorm2(lambda_t%w_ooov)*0.5D0
           print*, "w_oooo     = ", print_tensornorm2(lambda_t%w_oooo)*0.25D0
         end if

!****************************************
!        get CC2 Specific Intermediates *
!****************************************
         vo_id = (/vir_space_id, occ_space_id/)
         vo_root = (/vir_space_root, occ_space_root/)
         vvoo_id = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)

         ! comm_t_CC2 t1 and t2 = comm_t  
         ierr=exatns_tensor_create(comm_t_CC2%t1,"cc2_t1",vo_id,vo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(comm_t_CC2%t2,"cc2_t2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(comm_t_CC2%tau,"cc2_tau",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(comm_t_CC2%h_int,"cc2_h_int",vovo_id,vovo_root,EXA_DATA_KIND_C8)

         ierr=exatns_tensor_init(comm_t_CC2%t1,'ZERO')
         ierr=exatns_tensor_init(comm_t_CC2%t2,'ZERO')
         ierr=exatns_tensor_init(comm_t_CC2%tau,'ZERO')
         ierr=exatns_tensor_init(comm_t_CC2%h_int,'ZERO')

         ! Copy of t1, t2, foo, fov and fvv from comm_t to comm_t_CC2   L! :
         ierr=exatns_tensor_contract("H(a,i)+=T(a,i)*K()",comm_t_CC2%t1,comm_t%t1,one_tensor,ONE)
         ierr=exatns_tensor_contract("H(a,b,i,j)+=T(a,b,i,j)*K()",comm_t_CC2%t2,comm_t%t2,one_tensor,ONE)
         ierr=exatns_tensor_contract("H(a,i)+=F(a,i)*K()",comm_t_CC2%fov,comm_t%fov,one_tensor,ONE)
         ierr=exatns_tensor_contract("H(a,b)+=F(a,b)*K()",comm_t_CC2%fvv,comm_t%fvv,one_tensor,ONE)
         ierr=exatns_tensor_contract("H(i,j)+=F(i,j)*K()",comm_t_CC2%foo,comm_t%foo,one_tensor,ONE)
         ierr=exatns_tensor_contract("W(b,m,e,j)+=H(b,m,e,j)*K()",comm_t_CC2%h_int,comm_t%h_int,one_tensor,ONE)

         call get_tauCC2(comm_t_CC2,one_tensor)

         if(print_CC2_inter) then
          print*, "t1 commt         = ", print_tensornorm2(comm_t%t1)
          print*, "t1 commt  CC2    = ", print_tensornorm2(comm_t_CC2%t1)
          print*, "t2 commt         = ", print_tensornorm2(comm_t%t2)
          print*, "t2 commt  CC2    = ", print_tensornorm2(comm_t_CC2%t2)
          print*, "Tau commt         = ", print_tensornorm2(comm_t%tau)
          print*, "Tau commt  CC2    = ", print_tensornorm2(comm_t_CC2%tau)
         end if

         call get_fbar_vv_CC2(lambda_t_CC2%fbar_vv,comm_t_CC2,int_t,one_tensor) ! for CCSD: fbar_vv = gvv  (L! : change ?)
         call get_fbar_oo_CC2(lambda_t_CC2%fbar_oo,comm_t_CC2,int_t,one_tensor)

         call get_w_vovo_CC2(lambda_t_CC2%w_vovo,comm_t_CC2,int_t%oovv,vvoo_id,vvoo_root,one_tensor)
         call get_w_oooo_CC2(lambda_t_CC2%w_oooo,comm_t_CC2,int_t,one_tensor)
         call get_w_vvvv_CC2(lambda_t_CC2%w_vvvv,comm_t_CC2,int_t,one_tensor)

         if(print_CC2_inter) then
          print*, "fvv bar commt         = ", print_tensornorm2(lambda_t%fbar_vv)
          print*, "fvv bar commt  CC2    = ", print_tensornorm2(lambda_t_CC2%fbar_vv)
          print*, "foo bar commt         = ", print_tensornorm2(lambda_t%fbar_oo)
          print*, "foo bar commt  CC2    = ", print_tensornorm2(lambda_t_CC2%fbar_oo)
          print*, "wvovo commt         = ", print_tensornorm2(lambda_t%w_vovo)
          print*, "wvovo commt  CC2    = ", print_tensornorm2(lambda_t_CC2%w_vovo)
          print*, "woooo commt         = ", print_tensornorm2(lambda_t%w_oooo)
          print*, "woooo commt  CC2    = ", print_tensornorm2(lambda_t_CC2%w_oooo)
          print*, "wvvvv commt         = ", print_tensornorm2(lambda_t%w_vvvv)
          print*, "wvvvv commt  CC2    = ", print_tensornorm2(lambda_t_CC2%w_vvvv)
         end if

!************************************
!        construct solution tensors *
!************************************

         if (.not.CCD) then
            ierr=exatns_tensor_create(s1_tensor,"s1",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(s1_tensor,'ZERO')
         end if 

         ierr=exatns_tensor_create(s2_tensor,"s2",oovv_id,oovv_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(s2_tensor,'ZERO')

!********************************************
!        construct mixed type intermediates *
!********************************************

         ierr=exatns_tensor_create(lambda_t%goo,"goo_mixed",oo_id,oo_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_create(lambda_t%gvv,"gvv_mixed",vv_id,vv_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(lambda_t%goo,'ZERO')
         ierr=exatns_tensor_init(lambda_t%gvv,'ZERO')

!*****************************************************************
!        re-initialize diis, so it can be used in lambda routine *
!*****************************************************************

         if (.not.CCD) then
           call diis (l2_tensor,oovv_id,oovv_root,l1_tensor,ov_id,ov_root,mode=1)
         else
           call diis (l2_tensor,oovv_id,oovv_root,mode=1)
         end if

!*****************
!        iterate *
!*****************

         do ilambda = 1, ncycles

!***********************************
!         form mixed intermediates *
!***********************************

          call get_gvv_lambda(lambda_t%gvv,l2_tensor,comm_t%t2)
          call get_goo_lambda(lambda_t%goo,l2_tensor,comm_t%t2)

!         debug
!         if (print_level.gt.8) then
            if (.not.CCD) print*, "L1  = ", print_tensornorm2(l1_tensor)
            print*, "L2  = ", print_tensornorm2(l2_tensor)*0.25D0
            print*, "goo = ", print_tensornorm2(lambda_t%goo)
            print*, "gvv = ", print_tensornorm2(lambda_t%gvv)
!         end if

!******************************
!         calculate s1 tensor *
!******************************

          if (.not.CCD) then
            ierr=exatns_tensor_init(s1_tensor,'ZERO')

!----------------------------------------
!           term 1: S(i,a) += Fbar(i,a) |
!----------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=F(i,a)*K()",s1_tensor,lambda_t%fbar_ov,one_tensor)
!-------------------------------------------------
!           term 2: S(i,a) += L(i,e) * Fbar(e,a) |
!-------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,l1_tensor,lambda_t%fbar_vv)

!-------------------------------------------------
!           term 3: S(i,a) -= L(m,a) * Fbar(i,m) |
!-------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,l1_tensor,lambda_t%fbar_oo,MINUS_ONE)

!------------------------------------------------
!           term 6: S(i,a) -= G(n,m) * W(mi,na) |
!------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,lambda_t%goo,lambda_t%w_ooov)

!------------------------------------------------
!           term 7: S(i,a) += L(m,e) * W(ie,am) |
!------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,l1_tensor,lambda_t%w_vovo,MINUS_ONE)

!---------------------------------------------------------
!           term 8: S(i,a) -= 1/2 * L(mn,ae)  * W(ie,mn) |
!---------------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,l2_tensor,lambda_t%w_ovoo,ONE_HALF)

!------------------------------------------------------
!           term 4: S(i,a) += 1/2 L(im,ef) * W(ef,am) |
!------------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,l2_tensor,lambda_t%w_vvvo,ONE_HALF)

!------------------------------------------------
!           term 5: S(i,a) -= G(f,e) * W(ei,fa) |
!------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,lambda_t%gvv,lambda_t%w_vovv)
          end if

!******************************
!         calculate s2 tensor *
!******************************

          ierr=exatns_tensor_init(s2_tensor,'ZERO')

!---------------------------------------
!         term 1: S(ij,ab) += V(ij,ab) |
!---------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)

!----------------------------------------------------------
!         term 2: S(ij,ab) += P(a,b) L(ij,ae) * Fbar(e,b) |  Fbar_vv - CC2
!----------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,l2_tensor,lambda_t_CC2%fbar_vv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,l2_tensor,lambda_t_CC2%fbar_vv)

!----------------------------------------------------------
!         term 3: S(ij,ab) -= P(i,j) L(im,ab) * Fbar(j,m) |  Fbar_oo - CC2
!----------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,l2_tensor,lambda_t_CC2%fbar_oo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,l2_tensor,lambda_t_CC2%fbar_oo)

!------------------------------------------------------
!         term 4: S(ij,ab) += 1/2 L(mn,ab) * W(ij,mn) |   W_oooo - CC2
!------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,l2_tensor,lambda_t%w_oooo,ONE_HALF)

!----------------------------------------------------------------
!         term 5: S(ij,ab) += P(i,j) P(a,b) L(im,ae) * W(je,bm) |  W vovo - CC2
!----------------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)

!-------------------------------------------------------
!         term 7: S(ij,ab) -= P(a,b) L(m,a) * W(,j,mb) |
!-------------------------------------------------------

          if (.not.CCD) then
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,l1_tensor,lambda_t%w_ooov)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",s2_tensor,l1_tensor,lambda_t%w_ooov)
          end if

!----------------------------------------------------------
!         term 9+10: S(ij,ab) += P(i,j) L(i,e) * W(ej,ab) |
!----------------------------------------------------------

          if (.not.CCD) then
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",s2_tensor,l1_tensor,lambda_t%w_vovv)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",s2_tensor,l1_tensor,lambda_t%w_vovv)
          end if

!----------------------------------------------------------------
!         term 11: S(ij,ab) += P(i,j) P(a,b) L(i,a) * Fbar(j,b) |
!----------------------------------------------------------------

          if (.not.CCD) then
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov,MINUS_ONE)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov,MINUS_ONE)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
          end if

!-----------------------------------------------------
!         term 12: S(ij,ab) += 1/2 L(ij,ef) W(ef,ab) |
!-----------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,l2_tensor,lambda_t%w_vvvv,ONE_HALF)

!!! L! End of S1 construction


          if (ierr.ne.EXA_SUCCESS) then
            write(*,*) " ccdriver: error in evaluation of lambda equations, code=",ierr
            return
          end if

!*************************************
!         prepare for next iteration *
!*************************************

          if (.not.CCD) then
            ierr=exatns_tensor_init(l1_tensor,'ZERO')
            ierr=exatns_tensor_contract("L1(i,a)+=S1(i,a)*K()",l1_tensor,s1_tensor,one_tensor)
            ierr=exatns_tensor_transform(l1_tensor,denom)
          end if
          ierr=exatns_tensor_init(l2_tensor,'ZERO')
          ierr=exatns_tensor_contract("L2(i,j,a,b)+=S2(i,j,a,b)*K()",l2_tensor,s2_tensor,one_tensor)
          ierr=exatns_tensor_transform(l2_tensor,denom)

          if (.not.CCD) then
            call diis (l2_tensor,t1_tensor=l1_tensor,t_convergence=t_convergence)
          else
            call diis (l2_tensor,t_convergence=t_convergence)
          end if

          call print_iteration(ilambda,t_convergence)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         if (.not.CCD) write(*,'(A,ES22.16)') "  Final L1amp = ", print_tensornorm2(l1_tensor)
         write(*,'(A,ES22.16)') "  Final L2amp = ", print_tensornorm2(l2_tensor)*0.25D0
         write(*,*) ""

!*****************
!        cleanup *
!*****************

         if (.not.CCD) then
           ierr=exatns_tensor_destroy(s1_tensor)
           call diis(l2_tensor,t1_tensor=l1_tensor,mode=2)
         else
           call diis(l2_tensor,mode=2)
         end if
         ierr=exatns_tensor_destroy(s2_tensor)
         ierr=exatns_tensor_destroy(one_tensor)
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

        end subroutine solve_CC2_lambda_equations


        subroutine get_CC_integrals (ao_space_id,ao_space,int_t,moint_scheme,print_level)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list. We do assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.

!        Written by Lucas Visscher, summer 2017

         use exacorr_ao_to_mo
         use exacorr_global

!        information needed to compute the ao-tensor
         integer(INTD),intent(in)              :: ao_space_id
         class(h_space_t), pointer, intent(in) :: ao_space
!        target integral tensors to be filled and given back to the caller
         type(exatns_intg_tens), intent(inout) :: int_t
         integer, intent(in)                   :: print_level 

         integer(INTL)              :: ao_space_root,subspace_id(4),num_subspaces
         integer(INTL),allocatable  :: ao_subspaces(:)
         type(tens_rcrsv_t)         :: a2int_tensor, ovvo_tensor ! Auxilliary tensors
         integer(INTL),dimension(4) :: a2int_root, vovo_root, ovvo_root
         integer(INTD),dimension(4) :: a2int_id, vovo_id, ovvo_id
         
         integer(INTD)                     :: ierr
         integer                           :: i, sub3, sub4
         integer                           :: moint_scheme ! (0: put zeroes for testing. 1: n^5 algorithm, 2: n^6 algorithm
                                                           !  3: n^5 with re-use of ints 4: as 3, but with sliced AO and HT)
         character(len=15)                 :: text

!        auxilliary integral tensor (only needed inside this routine)
         type(tens_rcrsv_t) :: aoint_tensor

!        get the root space, covers all ao's
         ao_space_root=ao_space%get_root_id(ierr)

!        Determine space ids and root for auxilliary voov tensor from vovo tensor
         call int_t%vovo%get_space_ids(vovo_id,vovo_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('vovo tensor corrupted')
         ovvo_id(1)   = vovo_id(2)
         ovvo_id(2)   = vovo_id(1)
         ovvo_id(3:4) = vovo_id(3:4)
         ovvo_root(1)   = vovo_root(2)
         ovvo_root(2)   = vovo_root(1)
         ovvo_root(3:4) = vovo_root(3:4)

!        initialise tensors to zero
         ierr=exatns_tensor_init(int_t%oooo,'ZERO')
         ierr=exatns_tensor_init(int_t%ooov,'ZERO')
         ierr=exatns_tensor_init(int_t%oovv,'ZERO')
         ierr=exatns_tensor_init(int_t%vovo,'ZERO')
         ierr=exatns_tensor_init(int_t%vovv,'ZERO')
         ierr=exatns_tensor_init(int_t%vvvv,'ZERO')

!        Get tensor with AO integrals
         if (moint_scheme < 4 ) then
            ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            call print_date('allocated AO integral tensor')
         end if

         !Initialize ao-tensor:
         if (moint_scheme == 0) then
!        test code : initialize integrals to zero to be able to time the integral generation
           ierr=exatns_tensor_init(aoint_tensor)
           call print_date('initalized all ao integrals to zero')
         elseif (moint_scheme == 1 .or. moint_scheme == 2) then
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           if(print_level > 6) write(*,*) " N (aoint) = ", print_tensornorm2(aoint_tensor)
           call print_date('finished ao integral calculation')

!          Get vvvv tensor
           call get_integral_tensor (aoint_tensor,int_t%vvvv,12,print_level,algorithm=moint_scheme)
           if(print_level > 6) write(*,*) " N (vvvv) = ",print_tensornorm2(int_t%vvvv)
           call print_date('finished vvvv integrals')

!          Get vovv tensor
           call get_integral_tensor (aoint_tensor,int_t%vovv,34,print_level,algorithm=moint_scheme)
           if(print_level > 6) write(*,*) " N (vovv) =",print_tensornorm2(int_t%vovv)
           call print_date('finished vovv integrals')

!          Get oovv tensor
           call get_integral_tensor (aoint_tensor,int_t%oovv,12,print_level,algorithm=moint_scheme)
           if(print_level > 6) write(*,*) " N (oovv) =",print_tensornorm2(int_t%oovv)
           call print_date('finished oovv integrals')

!          Get ooov tensor
           call get_integral_tensor (aoint_tensor,int_t%ooov,12,print_level,algorithm=moint_scheme)
           if(print_level > 6) write(*,*) " N (ooov) =",print_tensornorm2(int_t%ooov)
           call print_date('finished ooov integrals')

!          Get oooo tensor
           call get_integral_tensor (aoint_tensor,int_t%oooo,12,print_level,algorithm=moint_scheme)
           if(print_level > 6) write(*,*) " N (oooo) =",print_tensornorm2(int_t%oooo)
           call print_date('finished oooo integrals')

!          Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!          Get vovo tensor
           call get_integral_tensor (aoint_tensor,int_t%vovo,0,print_level,algorithm=moint_scheme)
           ierr=exatns_tensor_create(ovvo_tensor,"ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tens. too large ')
           ierr=exatns_tensor_init(ovvo_tensor,'ZERO')
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: initialization failed ')
           call get_integral_tensor (aoint_tensor,ovvo_tensor,0,print_level,algorithm=moint_scheme)
!          Make anti-symmetrized vovo tensor
           ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*K()",int_t%vovo,ovvo_tensor,one_tensor,MINUS_ONE)
           ierr=exatns_tensor_destroy(ovvo_tensor)
           if(print_level > 6) write(*,*) " N (vovo) =",print_tensornorm2(int_t%vovo)
           call print_date('finished vovo integrals')

         elseif (moint_scheme == 3) then
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           call print_date('finished ao integral calculation')
           call print_date('choosing moint_scheme 3')

           ! Setting dimensions for the half transformed integral tensor, the last two are always AO spaces
           a2int_id(3:4)   = ao_space_id
           a2int_root(3:4) = ao_space_root

           ! get the integrals with vv as first transform
           a2int_id(1:2)   = vovo_id(1)
           a2int_root(1:2) = vovo_root(1)
           ierr=exatns_tensor_create(a2int_tensor,"a2vv",a2int_id,a2int_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ') 
           ierr=exatns_tensor_init(a2int_tensor,'ZERO')
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: init. failed ')
           call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
           call print_date('vv halftransformation completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%vvvv,12,print_level)
           if(print_level > 6) write(*,*) " N (vvvv) =",print_tensornorm2(int_t%vvvv)
           call print_date('vvvv integrals completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%vovv,34,print_level)
           if(print_level > 6) write(*,*) " N (vovv) =",print_tensornorm2(int_t%vovv)
           call print_date('vovv integrals completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%vovo,0,print_level)
           call print_date('vovo integrals completed')
           ierr=exatns_tensor_destroy(a2int_tensor)

           ! get the integrals with ov as first transform
           a2int_id(1:2)   = vovo_id(2:3)
           a2int_root(1:2) = vovo_root(2:3)
           ierr=exatns_tensor_create(a2int_tensor,"a2ov",a2int_id,a2int_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ')
           ierr=exatns_tensor_init(a2int_tensor,'ZERO')
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: init. failed ')
           call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
           call print_date('ov halftransformation completed')
           ierr=exatns_tensor_create(ovvo_tensor,"ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           ierr=exatns_tensor_init(ovvo_tensor,'ZERO')
           call get_ht2_integral_tensor (a2int_tensor,ovvo_tensor,0,print_level)
           call print_date('ovvo integrals completed')
!          Make anti-symmetrized vovo tensor
           ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*K()",int_t%vovo,ovvo_tensor,one_tensor,MINUS_ONE)
           ierr=exatns_tensor_destroy(ovvo_tensor)
           if(print_level > 6) write(*,*) " N (vovo) =",print_tensornorm2(int_t%vovo)
           call print_date('vovo integrals anti-symmetrized')
           call get_ht2_integral_tensor (a2int_tensor,int_t%oovv,12,print_level)
           if(print_level > 6) write(*,*) " N (oovv) =",print_tensornorm2(int_t%oovv)
           call print_date('oovv integrals completed')
           ierr=exatns_tensor_destroy(a2int_tensor)

           ! get the integrals with oo as first transform
           a2int_id(1:2)   = vovo_id(2)
           a2int_root(1:2) = vovo_root(2)
           ierr=exatns_tensor_create(a2int_tensor,"a2oo",a2int_id,a2int_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ')
           ierr=exatns_tensor_init(a2int_tensor,'ZERO')
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: init. failed ')
           call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
           call print_date('oo halftransformation completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%ooov,12,print_level)
           if(print_level > 6) write(*,*) " N (ooov) =",print_tensornorm2(int_t%ooov)
           call print_date('ooov integrals completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%oooo,12,print_level)
           if(print_level > 6) write(*,*) " N (oooo) =",print_tensornorm2(int_t%oooo)
           call print_date('oooo integrals completed')
           ierr=exatns_tensor_destroy(a2int_tensor)

         elseif (moint_scheme == 4) then
           
           call print_date('choosing moint_scheme 4')
           ! in this algorithm we only keep part of the AO and HT integrals in memory
           subspace_id(1:3) = ao_space_root
           call ao_space%get_level_composition(int(1,INTD),ao_subspaces,num_subspaces,ierr)
           do sub4 = 1, num_subspaces
              subspace_id(4)=ao_subspaces(sub4)
              ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),subspace_id,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
              ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')
              if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')

              ! Setting dimensions for the half transformed integral tensor, the last two are always AO spaces
              a2int_id(3:4)   = ao_space_id
              a2int_root(3:4) = subspace_id(3:4)

              ! get the integrals with vv as first transform
              a2int_id(1:2)   = vovo_id(1)
              a2int_root(1:2) = vovo_root(1)
              ierr=exatns_tensor_create(a2int_tensor,"a2vv",a2int_id,a2int_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ')
              ierr=exatns_tensor_init(a2int_tensor,'ZERO')
              call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
              call get_ht2_integral_tensor (a2int_tensor,int_t%vvvv,12,print_level)
              call get_ht2_integral_tensor (a2int_tensor,int_t%vovv,34,print_level)
              call get_ht2_integral_tensor (a2int_tensor,int_t%vovo,0,print_level)
              ierr=exatns_tensor_destroy(a2int_tensor)

              ! get the integrals with ov as first transform
              a2int_id(1:2)   = vovo_id(2:3)
              a2int_root(1:2) = vovo_root(2:3)
              ierr=exatns_tensor_create(a2int_tensor,"a2ov",a2int_id,a2int_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ')
              ierr=exatns_tensor_init(a2int_tensor,'ZERO')
              call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
              ierr=exatns_tensor_create(ovvo_tensor,"ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals
              ierr=exatns_tensor_init(ovvo_tensor,'ZERO')
              call get_ht2_integral_tensor (a2int_tensor,ovvo_tensor,0,print_level)
!             Make anti-symmetrized vovo tensor
              ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*K()",int_t%vovo,ovvo_tensor,one_tensor,MINUS_ONE)
              ierr=exatns_tensor_destroy(ovvo_tensor)
              call get_ht2_integral_tensor (a2int_tensor,int_t%oovv,12,print_level)
              ierr=exatns_tensor_destroy(a2int_tensor)

              ! get the integrals with oo as first transform
              a2int_id(1:2)   = vovo_id(2)
              a2int_root(1:2) = vovo_root(2)
              ierr=exatns_tensor_create(a2int_tensor,"a2oo",a2int_id,a2int_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ')
              ierr=exatns_tensor_init(a2int_tensor,'ZERO')
              call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
              call get_ht2_integral_tensor (a2int_tensor,int_t%ooov,12,print_level)
              call get_ht2_integral_tensor (a2int_tensor,int_t%oooo,12,print_level)
              ierr=exatns_tensor_destroy(a2int_tensor)
              ierr=exatns_tensor_destroy(aoint_tensor)
              write (text,'(I6,A,I6)') sub4,' of',num_subspaces
              call print_date('finished transforming subspace'//text)
           end do

           if(print_level > 6) write(*,*) " N (vvvv) =",print_tensornorm2(int_t%vvvv)
           if(print_level > 6) write(*,*) " N (vovv) =",print_tensornorm2(int_t%vovv)
           if(print_level > 6) write(*,*) " N (vovo) =",print_tensornorm2(int_t%vovo)
           if(print_level > 6) write(*,*) " N (oovv) =",print_tensornorm2(int_t%oovv)
           if(print_level > 6) write(*,*) " N (ooov) =",print_tensornorm2(int_t%ooov)
           if(print_level > 6) write(*,*) " N (oooo) =",print_tensornorm2(int_t%oooo)

         end if

!        Clean-up
         if (moint_scheme < 4) then
            ierr=exatns_tensor_destroy(aoint_tensor)
         end if

        end subroutine get_CC_integrals

        subroutine diis (t2_tensor,vvoo_id,vvoo_root,t1_tensor,vo_id,vo_root,t_convergence,mode)

         integer, parameter             :: maxdim = 5 ! maximum number of tensors that will be stored
         integer, parameter             :: mindim = 3 ! minimum number of tensors to allow extrapolation
         integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
         real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes

         type(tens_rcrsv_t), intent(inout),optional :: t1_tensor
         type(tens_rcrsv_t), intent(inout)          :: t2_tensor
         integer(INTD)                              :: ierr, tens_rank
         integer                                    :: i
         integer(INTD),dimension(2),optional        :: vo_id
         integer(INTD),dimension(4),optional        :: vvoo_id
         integer(INTL),dimension(2),optional        :: vo_root
         integer(INTL),dimension(4),optional        :: vvoo_root
         logical                                    :: init, clean_up, CCD
         integer, save                              :: itens, ntens
         complex(8), save                           :: b_diis(maxdim,maxdim)
         complex(8)                                 :: bb(0:maxdim,0:maxdim),coeff(0:maxdim)
         complex(8)                                 :: scalar
         real(8)                                    :: b_max, b_scale
         type(tens_rcrsv_t), save                   :: t1_inp_tensor, t2_inp_tensor
         type(tens_rcrsv_t), save                   :: t1_res_tensor(maxdim), t2_res_tensor(maxdim)
         type(tens_rcrsv_t), save                   :: t1_err_tensor(maxdim), t2_err_tensor(maxdim)
         integer                                    :: ipiv(0:maxdim), info  ! needed for ZHESV solver
         real(8)                                    :: work(maxdim,maxdim)   ! needed for ZHESV solver
         type(tens_rcrsv_t), save                   :: scalar_tensor

         character(len=1) :: vs_str

         init       = .false.
         clean_up   = .false.
         CCD        = .false.

         if (present(mode)) then
            if (mode.eq.1) init       = .true.
            if (mode.eq.2) clean_up   = .true.
         end if

         if (.not.present(t1_tensor)) CCD = .true.

         if (init) then
!           get dimensions of the tensors
            if (.not.CCD) then
                tens_rank=t1_tensor%get_rank(ierr)
                if (ierr.ne.0 .or. tens_rank.ne.2) stop 'input t1-tensor corrupted in diis'
            end if 
            tens_rank=t2_tensor%get_rank(ierr)
            if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input t2-tensor corrupted in diis'
!           initialize DIIS matrix
            b_diis      = ZERO
!           create storage for input, result and error tensors (in memory for the time being)
            if (.not.CCD) ierr=exatns_tensor_create(t1_inp_tensor,"t1_inp_tensor",vo_id,vo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(t2_inp_tensor,"t2_inp_tensor",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            do i = 1, maxdim
               write (vs_str, '(I1)') i
               if (.not.CCD) then
                   ierr=exatns_tensor_create(t1_res_tensor(i),"t1_res_tensor_" // vs_str,vo_id,vo_root,EXA_DATA_KIND_C8)
                   ierr=exatns_tensor_create(t1_err_tensor(i),"t1_err_tensor_" // vs_str,vo_id,vo_root,EXA_DATA_KIND_C8)
               end if
               ierr=exatns_tensor_create(t2_res_tensor(i),"t2_res_tensor_" // vs_str,vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
               ierr=exatns_tensor_create(t2_err_tensor(i),"t2_err_tensor_" // vs_str,vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            end do
            if(ierr.ne.EXA_SUCCESS) call quit('failed to create storage for input')
            itens = 0
            ntens = 0
!           store the tensors in memory as current input tensors, needed to calculate the error vector later on
            if (.not.CCD) then
                ierr=exatns_tensor_init(t1_inp_tensor,'ZERO')
                ierr=exatns_tensor_contract("T1O(a,i)+=T1N(a,i)*K()",t1_inp_tensor,t1_tensor,one_tensor)
            end if
            ierr=exatns_tensor_init(t2_inp_tensor,'ZERO')
            ierr=exatns_tensor_contract("T2O(a,b,i,j)+=T2N(a,b,i,j)*K()",t2_inp_tensor,t2_tensor,one_tensor)
            ierr=exatns_tensor_create(scalar_tensor,"scalar_tensor",EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            return ! no extrapolation possible
         end if

         if (clean_up) then
!           de-allocate tensors
            if (.not.CCD) ierr=exatns_tensor_destroy(t1_inp_tensor)
            ierr=exatns_tensor_destroy(t2_inp_tensor)
            do i = 1, maxdim
                if (.not.CCD) then
                    ierr=exatns_tensor_destroy(t1_res_tensor(i))
                    ierr=exatns_tensor_destroy(t1_err_tensor(i))
                end if
                ierr=exatns_tensor_destroy(t2_res_tensor(i))
                ierr=exatns_tensor_destroy(t2_err_tensor(i))
            end do
            ierr=exatns_tensor_destroy(scalar_tensor)
            return
         end if

!        update counters
         itens = itens + 1
         if (itens.gt.maxdim) itens = 1       ! cycle because we already have the maximum number of stored amplitudes, will overwrite the oldest one
         ntens = ntens + 1                    ! number of tensors stored after returning
         if (ntens.gt.maxdim) ntens = maxdim  ! we have reached the maximum number of tensors to be used in DIIS

!        store the result tensors
         if (.not.CCD) then
            ierr=exatns_tensor_init(t1_res_tensor(itens),'ZERO')
            ierr=exatns_tensor_contract("T1R(a,i)+=T1(a,i)*K()",t1_res_tensor(itens),t1_tensor,one_tensor)
         end if
         ierr=exatns_tensor_init(t2_res_tensor(itens),'ZERO')
         ierr=exatns_tensor_contract("T2R(a,b,i,j)+=T2(a,b,i,j)*K()",t2_res_tensor(itens),t2_tensor,one_tensor)

!        calculate the error tensors by comparing the result of the t-equations with the tensors used as input
         if (.not.CCD) then
            ierr=exatns_tensor_init(t1_err_tensor(itens),'ZERO')
            ierr=exatns_tensor_contract("T1E(a,i)+=T1R(a,i)*K()",t1_err_tensor(itens),t1_tensor,one_tensor)
            ierr=exatns_tensor_contract("T1E(a,i)+=T1I(a,i)*K()",t1_err_tensor(itens),t1_inp_tensor,one_tensor,prefactor=MINUS_ONE)
         end if
         ierr=exatns_tensor_init(t2_err_tensor(itens),'ZERO')
         ierr=exatns_tensor_contract("T2E(a,b,i,j)+=T2R(a,b,i,j)*K()",t2_err_tensor(itens),t2_tensor,one_tensor)
         ierr=exatns_tensor_contract("T2E(a,b,i,j)+=T2I(a,b,i,j)*K()",t2_err_tensor(itens),t2_inp_tensor,one_tensor, &
                                                                                                prefactor=MINUS_ONE)

!        update the b_diis array that holds the scalar products of the error tensors
         do i = 1, ntens
            ierr=exatns_tensor_init(scalar_tensor)
            if (.not.CCD) ierr=exatns_tensor_contract("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
            ierr=exatns_tensor_contract("E()+=T+(a,b,i,j)*T(a,b,i,j)",scalar_tensor,t2_err_tensor(itens),t2_err_tensor(i))
            ierr=exatns_tensor_get_scalar(scalar_tensor,scalar)
            b_diis(itens,i)=scalar
            b_diis(i,itens)=dconjg(scalar)
         end do

         if (present(t_convergence)) then
            t_convergence = dsqrt(dble(b_diis(itens,itens)))
         end if 

!        transfer to a work array in the way requested by the solver and normalize the matrix
         b_max = maxval (abs(b_diis))
         if (b_max .gt. 1.D-30) then
             b_scale = 1.D0 / b_max
         else
             write(*,*) "WARNING: something is wrong in DIIS: found zero error for non-converged amplitudes"
             return
         end if
         bb = ZERO
         bb(1:ntens,0) = MINUS_ONE
         bb(0,1:ntens) = MINUS_ONE
         bb(1:ntens,1:ntens) = b_diis(1:ntens,1:ntens) * b_scale

!        solve the linear system of equations using the LAPACK solver ZHESV
         coeff(0)       = MINUS_ONE  ! on input this should be the right-hand-side of the matrix equation b coef = rhs
         coeff(1:ntens) = ZERO 
         call ZHESV ('L',ntens+1,1,bb,maxdim+1,ipiv,coeff,maxdim+1,work,maxdim*maxdim,info)

!        make the extrapolated amplitudes if the space is large enough
         if (ntens.ge.mindim) then
            if (.not.CCD) ierr=exatns_tensor_init(t1_tensor,'ZERO')
            ierr=exatns_tensor_init(t2_tensor,'ZERO')
            do i = 1, ntens
               if (.not.CCD) ierr=exatns_tensor_contract("T1(a,i)+=T1O(a,i)*K()",t1_tensor,t1_res_tensor(i), &
                                                                                one_tensor,prefactor=coeff(i))
               ierr=exatns_tensor_contract("T2(a,b,i,j)+=T2O(a,b,i,j)*K()",t2_tensor,t2_res_tensor(i),one_tensor,prefactor=coeff(i))
            end do
         end if
        
!        the extrapolated amplitudes are the new input tensors
         if (.not.CCD) then
            ierr=exatns_tensor_init(t1_inp_tensor,'ZERO')
            ierr=exatns_tensor_contract("T1O(a,i)+=T1(a,i)*K()",t1_inp_tensor,t1_tensor,one_tensor)
         end if
         ierr=exatns_tensor_init(t2_inp_tensor,'ZERO')
         ierr=exatns_tensor_contract("T2O(a,b,i,j)+=T2(a,b,i,j)*K()",t2_inp_tensor,t2_tensor,one_tensor)

        end subroutine diis

        subroutine init_fock_exat (comm_t,int_t,oo_id,oo_root,ov_id,ov_root,vv_id,vv_root,&
                              f_delta,denom,orb_set,ff_set,exa_input,scf_energy)

          use exacorr_global
          use exacorr_tensor_methods
          use exacorr_ao_to_mo, only : print_tensornorm2

          implicit none
  
          type(exatns_comm_tens), intent(inout)  :: comm_t
          type(exatns_intg_tens), intent(inout)  :: int_t
          integer(INTL),dimension(2), intent(in) :: oo_root, ov_root, vv_root
          integer(INTD),dimension(2), intent(in) :: oo_id, ov_id, vv_id 
          type(delta_t)                          :: f_delta
          type(denom_t)                          :: denom
          type(set_energy_t)                     :: orb_set
          type(set_ff_t)                         :: ff_set
          type(exacc_input), intent(in)          :: exa_input
          real(8), intent(out)                   :: scf_energy
          integer                                :: iff=1   !TODO: make a loop for field strngths

          real(8)                    :: ev_energy, off_diag
          !        scalars (need to be defined as tensor types)
          type(tens_rcrsv_t)         :: val_tensor, tensor_one, tensor_h
          !        identity matrix for contraction 
          type(tens_rcrsv_t)         :: id_mat_occ, orb_tensor
          integer(INTD)              :: ierr
          complex(8)                 :: result_val
          !        finite field stuff
          integer                    :: i,i_prop

          if (one_el_exist) then

            ! Allocate two auxilliary tensors needed to extract results or multiply by one
            ierr=exatns_tensor_create(val_tensor,"val_tensor",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(tensor_one,"tensor_one",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(tensor_one, ONE)

            ! Make the density matrix (which is an identity matrix in the MO basis)
            ierr=exatns_tensor_create(id_mat_occ,"id_mat_occ",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_transform(id_mat_occ,'Unity')
            if (exa_input%print_level >0) call print_date('Density matrix initialized')

            ! Initialize the sub blocks of the Fock matrix with the (effective) one-electron integrals
            ierr=exatns_tensor_init(comm_t%foo,'OE_oo')
            ierr=exatns_tensor_init(comm_t%fov,'OE_ov')
            ierr=exatns_tensor_init(comm_t%fvv,'OE_vv')
            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'Creating one electron integrals failed!')
            if (exa_input%print_level >0) call print_date('one electron integrals added to Fock Matrix')

            !adding property integrals
            do i=1,exa_input%nff(1)
              call get_ff_ind (i_prop,exa_input%ff_names(i))
              call ff_set%reset(i_prop)

              ierr=exatns_tensor_create(tensor_h,"tensor_h",oo_id,oo_root,EXA_DATA_KIND_C8)
              write(*,*) 'ff_set create error', ierr
              ierr=exatns_tensor_init(tensor_h, ff_set)
              write(*,*) 'ff_set init error', ierr
              ierr=exatns_tensor_contract("F(p,q)+=P(p,q)*O()",comm_t%foo,tensor_h,tensor_one,exa_input%ff(i,iff))
              write(*,*) 'ff_set contract error', ierr
              write(*,*) " prop =",print_tensornorm2(tensor_h)
              write(*,*) " field =",exa_input%ff(i,iff)
              ierr=exatns_tensor_traverse(tensor_h,'PrintTensor')
              ierr=exatns_tensor_destroy(tensor_h)
              write(*,*) " F_oo+prop =",print_tensornorm2(comm_t%foo)
              write(*,*) 'ff_set destroy error', ierr
              ierr=exatns_tensor_traverse(comm_t%foo,'PrintTensor')
            
            end do


            ! Calculate the Hartree-Fock energy
            ierr=exatns_tensor_init(val_tensor, ZERO)
            ierr=exatns_tensor_contract("R()=F(i,j)*I(i,j)",val_tensor,comm_t%foo,id_mat_occ)
            ierr=exatns_tensor_create(tensor_h,"tensor_h",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_get_scalar(val_tensor,result_val)
            ierr=exatns_tensor_init(tensor_h, 'ZERO')
            ierr=exatns_tensor_contract("H(i,j)+=V(i,k,j,l)*I(k,l)",tensor_h,int_t%oooo,id_mat_occ, ONE_HALF)
            ierr=exatns_tensor_contract("R()+=H(i,j)*I(i,j)",val_tensor,tensor_h,id_mat_occ)
            ierr=exatns_tensor_get_scalar(val_tensor,result_val)
            scf_energy=real(result_val,8) ! Casting into real(8)
            ierr=exatns_tensor_destroy(tensor_h)
            scf_energy=scf_energy+e_core
            write(*,*) ""
            write (*,'(A,F23.15)') "  SCF energy (single determinant, calculated) = ", scf_energy
            write(*,*) ""
            ! Make the complete Fock matrix by addding the valence electron contributions
            ierr=exatns_tensor_contract("F(i,j)+=V(i,k,j,l)*I(k,l)",comm_t%foo,int_t%oooo,id_mat_occ, ONE)
            ierr=exatns_tensor_contract("F(i,a)+=V(k,i,l,a)*I(k,l)",comm_t%fov,int_t%ooov,id_mat_occ, ONE)
            ierr=exatns_tensor_contract("F(a,b)+=V(a,k,b,l)*I(k,l)",comm_t%fvv,int_t%vovo,id_mat_occ, ONE)
            if (exa_input%print_level >0) call print_date('two electron integrals added to Fock Matrix')

            if (exa_input%print_level >10) then
              ierr=exatns_tensor_traverse(comm_t%foo,'PrintTensor')
              ierr=exatns_tensor_traverse(comm_t%fov,'PrintTensor')
              ierr=exatns_tensor_traverse(comm_t%fvv,'PrintTensor')
            end if

            call update_eps(comm_t, f_delta, denom, orb_set, oo_root, vv_root, oo_id, vv_id, exa_input%print_level)

            ! Subtract orbital energies from the diagonal (they are moved to the other side of the CC equations)
            ierr=exatns_tensor_create(orb_tensor,"orb_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            ierr=exatns_tensor_init(orb_tensor,orb_set)
            if(ierr.ne.EXA_SUCCESS) call quit('init_fock: set orbitals failed!')
            ierr=exatns_tensor_contract("F(i,j)+=E(i,j)*O()",comm_t%foo,orb_tensor,tensor_one, MINUS_ONE)
            ierr=exatns_tensor_destroy(orb_tensor)
            ierr=exatns_tensor_create(orb_tensor,"orb_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            ierr=exatns_tensor_init(orb_tensor,orb_set)
            if(ierr.ne.EXA_SUCCESS) call quit('init_fock: set orbitals failed!')
            ierr=exatns_tensor_contract("F(a,b)+=E(a,b)*O()",comm_t%fvv,orb_tensor,tensor_one, MINUS_ONE)
            ierr=exatns_tensor_destroy(orb_tensor)
            if (exa_input%print_level >0) call print_date('orbital energies subtracted from the diagonal')
            
            if (exa_input%print_level >0) then
              ! Compute norms (check whether this single determinant is contructed from optimized orbitals)
              write(*,*) ""
              ierr=exatns_tensor_init(val_tensor, ZERO)
              ierr=exatns_tensor_contract("R()=F(i,j)*F(i,j)",val_tensor,comm_t%foo,comm_t%foo)
              ierr=exatns_tensor_get_scalar(val_tensor,result_val)
              off_diag=real(result_val,8) ! Casting into real(8)
              write(*,*) "    Norm of Fock matrix (OO, no diagonal) = ", off_diag
            
              ierr=exatns_tensor_init(val_tensor, ZERO)
              ierr=exatns_tensor_contract("R()=F(i,a)*F(i,a)",val_tensor,comm_t%fov,comm_t%fov)
              ierr=exatns_tensor_get_scalar(val_tensor,result_val)
              off_diag=real(result_val,8) ! Casting into real(8)
              write(*,*) "    Norm of Fock matrix (OV, all        ) = ", off_diag

              ierr=exatns_tensor_init(val_tensor, ZERO)
              ierr=exatns_tensor_contract("R()=F(a,b)*F(a,b)",val_tensor,comm_t%fvv,comm_t%fvv)
              ierr=exatns_tensor_get_scalar(val_tensor,result_val)
              off_diag=real(result_val,8) ! Casting into real(8)
              write(*,*) "    Norm of Fock matrix (VV, no diagonal) = ", off_diag
              write(*,*) ""
            end if 
            ierr=exatns_tensor_destroy(tensor_h)
            ierr=exatns_tensor_destroy(id_mat_occ)
            ierr=exatns_tensor_destroy(val_tensor)
            ierr=exatns_tensor_destroy(tensor_one)
          else
            ierr=exatns_tensor_init(comm_t%foo,'ZERO')
            ierr=exatns_tensor_init(comm_t%fov,'ZERO')
            ierr=exatns_tensor_init(comm_t%fvv,'ZERO')

            call get_scf ( scf_energy )
            write(*,*) ""
            write (*,'(A,F23.15)') "    SCF energy ( read ) = ", scf_energy
            write(*,*) ""

          end if

        end subroutine init_fock_exat

        subroutine update_eps(comm_t, f_delta, denom, orb_set, oo_root, vv_root, oo_id, vv_id, print_level)
          
          use exacorr_tensor_methods
          use exacorr_utils, only : print_orbital_energy

          implicit none

          type(exatns_comm_tens)                 :: comm_t
          type(delta_t)                          :: f_delta 
          type(denom_t)                          :: denom
          type(set_energy_t)                     :: orb_set 
          integer(INTL),dimension(2), intent(in) :: oo_root, vv_root
          integer(INTD),dimension(2), intent(in) :: oo_id, vv_id 
          integer, intent(in)                    :: print_level

          real(8),allocatable        :: eps_occ(:),eps_vir(:)
          integer                    :: nocc, nvir, i
          integer(INTL),dimension(2) :: delta_set
          type(tens_rcrsv_t)         :: d_mat, res_val
          complex(8)                 :: res
          integer(INTD)              :: ierr

          nocc=oo_root(1)
          nvir=vv_root(1)

          allocate(eps_occ(nocc))
          allocate(eps_vir(nvir))

          ierr=exatns_tensor_create(res_val,"res_val",EXA_DATA_KIND_C8)

          ierr=exatns_tensor_create(d_mat,"d_mat",oo_id,oo_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_init(d_mat,'ZERO')
          do i = 1, nocc
              delta_set=i
              call f_delta%reset(delta_set,2)
              ierr=exatns_tensor_transform(d_mat,f_delta)

              ierr=exatns_tensor_init(res_val,ZERO)
              ierr=exatns_tensor_contract("R()+=W(p,q)*C(p,q)",res_val,comm_t%foo,d_mat)
              ierr=exatns_tensor_get_scalar(res_val,res)
              eps_occ(i)=real(res,8)

          end do
          ierr=exatns_tensor_destroy(d_mat)

          ierr=exatns_tensor_create(d_mat,"d_mat",vv_id,vv_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_init(d_mat,'ZERO')
          do i = 1, nvir
              delta_set=i
              call f_delta%reset(delta_set,2)
              ierr=exatns_tensor_transform(d_mat,f_delta)

              ierr=exatns_tensor_init(res_val,ZERO)
              ierr=exatns_tensor_contract("R()+=W(p,q)*C(p,q)",res_val,comm_t%fvv,d_mat)
              ierr=exatns_tensor_get_scalar(res_val,res)
              eps_vir(i)=real(res,8)

          end do
          ierr=exatns_tensor_destroy(d_mat)
          
          if (print_level.gt.3) call print_orbital_energy(eps_occ,nocc,eps_vir,nvir,2)

          call denom%reset(eps_occ,eps_vir,nocc,nvir)
          call orb_set%reset(eps_occ,eps_vir,nocc,nvir)

          deallocate(eps_occ)
          deallocate(eps_vir)
          ierr=exatns_tensor_destroy(res_val)

        end subroutine update_eps

        subroutine t3_energy_split(t_energy,comm_t,int_t,occ_space,denom,exa_input)
            
            use exacorr_ao_to_mo, only : print_tensornorm2
            use exacorr_tensor_methods
            use interface_to_mpi

            real(8), intent(out)                    :: t_energy(3)
            type(exatns_comm_tens)                  :: comm_t
            type(exatns_intg_tens)                  :: int_t
            class(h_space_t), pointer               :: occ_space
            type(exacc_input), intent(in)           :: exa_input
            type(denom_t)                           :: denom

            integer(INTL),dimension(2) :: ov_root, vo_root, dim2_root
            integer(INTD),dimension(2) :: ov_id, vo_id
            integer(INTL),dimension(4) :: ooov_root, vovv_root, oovv_root, vvoo_root, dim4_root
            integer(INTD),dimension(4) :: ooov_id, vovv_id, vvoo_id, oovv_id, oooo_id
            integer(INTL),dimension(6) :: w_sub
            integer(INTD),dimension(6) :: w_id

            integer(INTD)                :: ierr,rank
            integer(INTD)                :: idocc, idvir
            integer(INTL)                :: nocc, nvir
            integer(INTL)                :: i, j, k, k0, n_subsp
            integer(INTL),allocatable    :: subsp(:)
            type(tens_rcrsv_t)           :: res_tensor
            complex(8)                   :: res_val, ONE_THIRD, prefact
            type(tens_rcrsv_t)           :: w, y, temp, temp2
            type(tens_rcrsv_t)           :: t1_i, t1_j, t1_k
            type(tens_rcrsv_t)           :: t2_i, t2_j, t2_k
            type(tens_rcrsv_t)           :: t2_jk, t2_ki, t2_ij
            type(tens_rcrsv_t)           :: vovv_i, vovv_j, vovv_k
            type(tens_rcrsv_t)           :: ooov_jk, ooov_ki, ooov_ij
            type(tens_rcrsv_t)           :: oovv_jk, oovv_ki, oovv_ij
            type(tens_rcrsv_t)           :: f_i, f_j, f_k
            type(tens_rcrsv_t)           :: P_i, P_j, P_k     ! Projection tensors
            type(tens_rcrsv_t)           :: P_jk, P_ij, P_ki    
            integer                      :: nr_proc
            real(8)                      :: est, nseg
            character(len=40)            :: info_str

            if (exa_input%print_level >4) call print_date('entering t3_energy_split')
            if (exa_input%tripl_block(1).gt.0) then
              k0=exa_input%tripl_block(1)
            else
              k0=1
            end if

            t_energy = ZERO

            ONE_THIRD=ONE/3.00D0

            call comm_t%fov%get_space_ids(ov_id,ov_root,ierr)
            call comm_t%fov%get_dims(ov_root,rank,ierr)
            call comm_t%t1%get_space_ids(vo_id,vo_root,ierr)
            call comm_t%t1%get_dims(vo_root,rank,ierr)
            
            call int_t%ooov%get_space_ids(ooov_id,ooov_root,ierr)
            call int_t%ooov%get_dims(ooov_root,rank,ierr)
            call int_t%vovv%get_space_ids(vovv_id,vovv_root,ierr)
            call int_t%vovv%get_dims(vovv_root,rank,ierr)
            call int_t%oovv%get_space_ids(oovv_id,oovv_root,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: getting ids')
            call int_t%oovv%get_dims(oovv_root,rank,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: getting dims')
            
            nocc=oovv_root(1)
            nvir=oovv_root(3)
            idocc=oovv_id(1)
            idvir=oovv_id(3)
            
            oooo_id=idocc
            vvoo_id=idocc
            vvoo_id(1:2)=idvir
            vvoo_root=nocc
            vvoo_root(1:2)=nvir

            if (exa_input%print_level>8) then
               write(*,*) ' --- ov,   id:',ov_id,', root: ',ov_root
               write(*,*) ' --- vo,   id:',vo_id,', root: ',vo_root
               write(*,*) ' --- ooov, id:',ooov_id,', root: ',ooov_root
               write(*,*) ' --- vovv, id:',vovv_id,', root: ',vovv_root
               write(*,*) ' --- oovv, id:',oovv_id,', root: ',oovv_root
               write(*,*) ' --- vvoo, id:',vvoo_id,', root: ',vvoo_root
               write(*,*) '======================================'
            end if

!           Initialize scalars that are to be used as tensors in contractions
            ierr=exatns_tensor_create(res_tensor,"res_tensor",EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: creation of res_tensor')

            ! set fixed indices
            w_id=idocc
            w_id(1:3)=idvir
            w_sub=nvir

!           In order to get an efficient algorithm the number of indices that are split 
!           has to be adjusted
            call interface_mpi_comm_size(global_communicator, nr_proc)

            !get memory information for splitting
            i=nr_proc*exa_input%talsh_buff
            if (exa_input%tripl_block(2).gt.0) then
              nseg=exa_input%tripl_block(2)
            else
              nseg=max(real(nvir,8)**(1.0/3.0),3.0D0)
            end if
            est=(real(nvir)**3*16.0/1000.0**3+1.0)*nseg**3
            if (exa_input%print_level > 0) then
              write(*,*) ' --- Number of processes (np):      ',nr_proc
              write(*,*) ' --- size of segments (ns):         ',nseg
              write(*,*) ' --- avail. memory (np*TALSH_buff): ',i
              write(*,*) ' --- required memory (ns^3*VVV):    ',est
              write(*,*) '======================================'
            end if

            ! the occupied index are always split
            if (.TRUE.) then
              if (nocc/2<exa_input%exa_blocksize) then
                est=real(nocc)/2.0
              else
                est=real(exa_input%exa_blocksize)
              end if
              do k=1,50
                if (est/(2.0**(k-1))<nseg) then
                  i=k
                  exit
                end if
              end do
              call occ_space%get_level_composition(int(i,INTD),subsp,n_subsp,ierr)
              if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: getting level composition')
              if (exa_input%print_level > 0) then
                write(*,*) ' --- splitting occupied spaces: 1-3 index '
                write(*,*) ' --- using level: ',i
                write(*,*) ' --- number el. : ',n_subsp,', size el.: ',real(nocc)/real(n_subsp)
                write(*,*) ' --- subspace   : ',subsp(1:n_subsp)
                write(*,*) '======================================'
              end if
              
            else
                ! use the whole space - only for testing
                n_subsp=1
                allocate(subsp(n_subsp))
                subsp=nocc

            end if

            do k = k0, n_subsp
              w_sub(6)=subsp(k)
              
              dim2_root=nocc
              dim2_root(2)=subsp(k)
              ierr=exatns_tensor_create(P_k,"P_k",oovv_id(1:2),dim2_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('P_k creation failed!')
              ierr=exatns_tensor_transform(P_k,'Unity')
              if (exa_input%print_level.gt.10) ierr=exatns_tensor_traverse(P_k,'PrintTensor')

              dim4_root=vovv_root
              dim4_root(2)=subsp(k)
              ierr=exatns_tensor_create(vovv_k,"vovv_k",vovv_id,dim4_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('vovv_k creation failed!')
              ierr=exatns_tensor_init(vovv_k,'ZERO')
              ierr=exatns_tensor_contract("I(c,k,a,b)+=V(c,p,a,b)*D(p,k)",vovv_k,int_t%vovv,P_k)
              if(ierr.ne.EXA_SUCCESS) call quit('vovv_k extraction failed!')
              if(exa_input%print_level > 10) write(*,*) " vovv_k =",print_tensornorm2(vovv_k)

              dim4_root=vvoo_root
              dim4_root(3)=subsp(k)
              ierr=exatns_tensor_create(t2_k,"t2_k",vvoo_id,dim4_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('t2_k creation failed!')
              ierr=exatns_tensor_init(t2_k,'ZERO')
              ierr=exatns_tensor_contract("I(a,b,k,i)+=T(a,b,p,i)*D(p,k)",t2_k,comm_t%t2,P_k)
              if(ierr.ne.EXA_SUCCESS) call quit('t2_k extraction failed!')
              if(exa_input%print_level > 10) write(*,*) " t2_k =",print_tensornorm2(t2_k)

              dim2_root=vo_root
              dim2_root(2)=subsp(k)
              ierr=exatns_tensor_create(t1_k,"t1_k",vo_id,dim2_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('t1_k creation failed!')
              ierr=exatns_tensor_init(t1_k,'ZERO')
              ierr=exatns_tensor_contract("I(a,k)+=T(a,p)*D(p,k)",t1_k,comm_t%t1,P_k)
              if(ierr.ne.EXA_SUCCESS) call quit('t1_k extraction failed!')
              if(exa_input%print_level > 10) write(*,*) " t1_k =",print_tensornorm2(t1_k)

              dim2_root=ov_root
              dim2_root(1)=subsp(k)
              ierr=exatns_tensor_create(f_k,"f_k",ov_id,dim2_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('f_k creation failed!')
              ierr=exatns_tensor_init(f_k,'ZERO')
              ierr=exatns_tensor_contract("I(k,a)+=T(p,a)*D(p,k)",f_k,comm_t%fov,P_k)
              if(ierr.ne.EXA_SUCCESS) call quit('f_k extraction failed!')
              if(exa_input%print_level > 10) write(*,*) " f_k =",print_tensornorm2(f_k)

              do j = 1, k
                w_sub(5)=subsp(j)

                dim2_root=nocc
                dim2_root(2)=subsp(j)
                ierr=exatns_tensor_create(P_j,"P_j",oovv_id(1:2),dim2_root,EXA_DATA_KIND_C8)
                if(ierr.ne.EXA_SUCCESS) call quit('P_j creation failed!')
                ierr=exatns_tensor_transform(P_j,'Unity')

                dim4_root=vovv_root
                dim4_root(2)=subsp(j)
                ierr=exatns_tensor_create(vovv_j,"vovv_j",vovv_id,dim4_root,EXA_DATA_KIND_C8)
                if(ierr.ne.EXA_SUCCESS) call quit('vovv_j creation failed!')
                ierr=exatns_tensor_init(vovv_j,'ZERO')
                ierr=exatns_tensor_contract("I(c,j,a,b)+=V(c,p,a,b)*D(p,j)",vovv_j,int_t%vovv,P_j)
                if(ierr.ne.EXA_SUCCESS) call quit('vovv_j extraction failed!')

                dim4_root=vvoo_root
                dim4_root(3)=subsp(j)
                ierr=exatns_tensor_create(t2_j,"t2_j",vvoo_id,dim4_root,EXA_DATA_KIND_C8)
                if(ierr.ne.EXA_SUCCESS) call quit('t2_j creation failed!')
                ierr=exatns_tensor_init(t2_j,'ZERO')
                ierr=exatns_tensor_contract("I(a,b,j,i)+=T(a,b,p,i)*D(p,j)",t2_j,comm_t%t2,P_j)
                if(ierr.ne.EXA_SUCCESS) call quit('t2_j extraction failed!')

                dim2_root=vo_root
                dim2_root(2)=subsp(j)
                ierr=exatns_tensor_create(t1_j,"t1_j",vo_id,dim2_root,EXA_DATA_KIND_C8)
                if(ierr.ne.EXA_SUCCESS) call quit('t1_j creation failed!')
                ierr=exatns_tensor_init(t1_j,'ZERO')
                ierr=exatns_tensor_contract("I(a,j)+=T(a,p)*D(p,j)",t1_j,comm_t%t1,P_j)
                if(ierr.ne.EXA_SUCCESS) call quit('t1_j extraction failed!')

                dim2_root=ov_root
                dim2_root(1)=subsp(j)
                ierr=exatns_tensor_create(f_j,"f_j",ov_id,dim2_root,EXA_DATA_KIND_C8)
                if(ierr.ne.EXA_SUCCESS) call quit('f_j creation failed!')
                ierr=exatns_tensor_init(f_j,'ZERO')
                ierr=exatns_tensor_contract("I(j,a)+=T(p,a)*D(p,j)",f_j,comm_t%fov,P_j)
                if(ierr.ne.EXA_SUCCESS) call quit('f_j extraction failed!')

                ! --------------------------------------- jk

                dim4_root=oovv_root(1)
                dim4_root(1)=subsp(j)
                dim4_root(3)=subsp(k)
                ierr=exatns_tensor_create(P_jk,"P_jk",oooo_id,dim4_root,EXA_DATA_KIND_C8)
                if(ierr.ne.EXA_SUCCESS) call quit('P_jk construction failed!')
                ierr=exatns_tensor_transform(P_jk,'SetProj')
                if(ierr.ne.EXA_SUCCESS) call quit('P_jk settting failed!')

                dim4_root=vvoo_root
                dim4_root(3)=subsp(j)
                dim4_root(4)=subsp(k)
                ierr=exatns_tensor_create(t2_jk,"t2_jk",vvoo_id,dim4_root,EXA_DATA_KIND_C8)
                ierr=exatns_tensor_init(t2_jk,'ZERO')
                ierr=exatns_tensor_contract("I(a,b,j,k)+=T(a,b,p,q)*P(j,p,k,q)",t2_jk,comm_t%t2,P_jk)
                if(ierr.ne.EXA_SUCCESS) call quit('t2_jk extraction failed!')
                if(exa_input%print_level > 10) write(*,*) " t2_jk =",print_tensornorm2(t2_jk)


                dim4_root=ooov_root
                dim4_root(1)=subsp(j)
                dim4_root(2)=subsp(k)
                ierr=exatns_tensor_create(ooov_jk,"ooov_jk",ooov_id,dim4_root,EXA_DATA_KIND_C8)
                if(ierr.ne.EXA_SUCCESS) call quit('ooov_jk creation failed!')
                ierr=exatns_tensor_init(ooov_jk,'ZERO')
                ierr=exatns_tensor_contract("I(j,k,i,a)+=T(p,q,i,a)*P(j,p,k,q)",ooov_jk,int_t%ooov,P_jk)
                if(ierr.ne.EXA_SUCCESS) call quit('ooov_jk extraction failed!')
                if(exa_input%print_level > 10) write(*,*) " ooov_jk =",print_tensornorm2(ooov_jk)

                dim4_root=oovv_root
                dim4_root(1)=subsp(j)
                dim4_root(2)=subsp(k)
                ierr=exatns_tensor_create(oovv_jk,"oovv_jk",oovv_id,dim4_root,EXA_DATA_KIND_C8)
                if(ierr.ne.EXA_SUCCESS) call quit('oovv_jk creation failed!')
                ierr=exatns_tensor_init(oovv_jk,'ZERO')
                ierr=exatns_tensor_contract("I(j,k,a,b)+=T(p,q,a,b)*P(j,p,k,q)",oovv_jk,int_t%oovv,P_jk)
                if(ierr.ne.EXA_SUCCESS) call quit('oovv_jk extraction failed!')
                if(exa_input%print_level > 10) write(*,*) " oovv_jk =",print_tensornorm2(oovv_jk)
                
                ierr=exatns_tensor_destroy(P_jk)

                do i = 1, j
                  w_sub(4)=subsp(i)

                  !determine value of prefactor
                  if (k.eq.j) then
                    if (i.eq.j) then
                      prefact=ONE/SIX
                    else
                      prefact=ONE_HALF
                    end if
                  else
                    if (i.eq.j) then
                      prefact=ONE_HALF
                    else
                      prefact=ONE
                    end if
                 end if


                  dim2_root=nocc
                  dim2_root(2)=subsp(i)
                  ierr=exatns_tensor_create(P_i,"P_i",oovv_id(1:2),dim2_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('P_i creation failed!')
                  ierr=exatns_tensor_transform(P_i,'Unity')

                  dim4_root=vovv_root
                  dim4_root(2)=subsp(i)
                  ierr=exatns_tensor_create(vovv_i,"vovv_i",vovv_id,dim4_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('vovv_i creation failed!')
                  ierr=exatns_tensor_init(vovv_i,'ZERO')
                  ierr=exatns_tensor_contract("I(c,i,a,b)+=V(c,p,a,b)*D(p,i)",vovv_i,int_t%vovv,P_i)
                  if(ierr.ne.EXA_SUCCESS) call quit('vovv_i extraction failed!')

                  dim4_root=vvoo_root
                  dim4_root(3)=subsp(i)
                  ierr=exatns_tensor_create(t2_i,"t2_i",vvoo_id,dim4_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('t2_i creation failed!')
                  ierr=exatns_tensor_init(t2_i,'ZERO')
                  ierr=exatns_tensor_contract("I(a,b,i,j)+=T(a,b,p,j)*D(p,i)",t2_i,comm_t%t2,P_i)
                  if(ierr.ne.EXA_SUCCESS) call quit('t2_i extraction failed!')

                  dim2_root=vo_root
                  dim2_root(2)=subsp(i)
                  ierr=exatns_tensor_create(t1_i,"t1_i",vo_id,dim2_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('t1_i creation failed!')
                  ierr=exatns_tensor_init(t1_i,'ZERO')
                  ierr=exatns_tensor_contract("I(a,i)+=T(a,p)*D(p,i)",t1_i,comm_t%t1,P_i)
                  if(ierr.ne.EXA_SUCCESS) call quit('t1_i extraction failed!')

                  dim2_root=ov_root
                  dim2_root(1)=subsp(i)
                  ierr=exatns_tensor_create(f_i,"f_i",ov_id,dim2_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('f_i creation failed!')
                  ierr=exatns_tensor_init(f_i,'ZERO')
                  ierr=exatns_tensor_contract("I(i,a)+=T(p,a)*D(p,i)",f_i,comm_t%fov,P_i)
                  if(ierr.ne.EXA_SUCCESS) call quit('f_i extraction failed!')

                  ! --------------------------------------- ij

                  dim4_root=oovv_root(1)
                  dim4_root(1)=subsp(i)
                  dim4_root(3)=subsp(j)
                  ierr=exatns_tensor_create(P_ij,"P_ij",oooo_id,dim4_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('P_ij construction failed!')
                  ierr=exatns_tensor_transform(P_ij,'SetProj')
                  if(ierr.ne.EXA_SUCCESS) call quit('P_ij settting failed!')

                  dim4_root=vvoo_root
                  dim4_root(3)=subsp(i)
                  dim4_root(4)=subsp(j)
                  ierr=exatns_tensor_create(t2_ij,"t2_ij",vvoo_id,dim4_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('t2_ij creation failed!')
                  ierr=exatns_tensor_init(t2_ij,'ZERO')
                  ierr=exatns_tensor_contract("I(a,b,i,j)+=T(a,b,p,q)*P(i,p,j,q)",t2_ij,comm_t%t2,P_ij)
                  if(ierr.ne.EXA_SUCCESS) call quit('t2_ij extraction failed!')

                  dim4_root=ooov_root
                  dim4_root(1)=subsp(i)
                  dim4_root(2)=subsp(j)
                  ierr=exatns_tensor_create(ooov_ij,"ooov_ij",ooov_id,dim4_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('ooov_ij creation failed!')
                  ierr=exatns_tensor_init(ooov_ij,'ZERO')
                  ierr=exatns_tensor_contract("I(i,j,k,a)+=T(p,q,k,a)*P(i,p,j,q)",ooov_ij,int_t%ooov,P_ij)
                  if(ierr.ne.EXA_SUCCESS) call quit('ooov_ij extraction failed!')

                  dim4_root=oovv_root
                  dim4_root(1)=subsp(i)
                  dim4_root(2)=subsp(j)
                  ierr=exatns_tensor_create(oovv_ij,"oovv_ij",oovv_id,dim4_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('oovv_ij creation failed!')
                  ierr=exatns_tensor_init(oovv_ij,'ZERO')
                  ierr=exatns_tensor_contract("I(i,j,a,b)+=T(p,q,a,b)*P(i,p,j,q)",oovv_ij,int_t%oovv,P_ij)
                  if(ierr.ne.EXA_SUCCESS) call quit('oovv_ij extraction failed!')
                  
                  ierr=exatns_tensor_destroy(P_ij)

                  ! --------------------------------------- ki
                  
                  dim4_root=oovv_root(1)
                  dim4_root(1)=subsp(k)
                  dim4_root(3)=subsp(i)
                  ierr=exatns_tensor_create(P_ki,"P_ki",oooo_id,dim4_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('P_ki construction failed!')
                  ierr=exatns_tensor_transform(P_ki,'SetProj')
                  if(ierr.ne.EXA_SUCCESS) call quit('P_ki settting failed!')
                  
                  dim4_root=vvoo_root
                  dim4_root(3)=subsp(k)
                  dim4_root(4)=subsp(i)
                  ierr=exatns_tensor_create(t2_ki,"t2_ki",vvoo_id,dim4_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('t2_ki creation failed!')
                  ierr=exatns_tensor_init(t2_ki,'ZERO')
                  ierr=exatns_tensor_contract("I(a,b,k,i)+=T(a,b,p,q)*P(k,p,i,q)",t2_ki,comm_t%t2,P_ki)
                  if(ierr.ne.EXA_SUCCESS) call quit('t2_ki extraction failed!')

                  dim4_root=ooov_root
                  dim4_root(1)=subsp(k)
                  dim4_root(2)=subsp(i)
                  ierr=exatns_tensor_create(ooov_ki,"ooov_ki",ooov_id,dim4_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('ooov_ki creation failed!')
                  ierr=exatns_tensor_init(ooov_ki,'ZERO')
                  ierr=exatns_tensor_contract("I(k,i,j,a)+=T(p,q,j,a)*P(k,p,i,q)",ooov_ki,int_t%ooov,P_ki)
                  if(ierr.ne.EXA_SUCCESS) call quit('ooov_ki extraction failed!')

                  dim4_root=oovv_root
                  dim4_root(1)=subsp(k)
                  dim4_root(2)=subsp(i)
                  ierr=exatns_tensor_create(oovv_ki,"oovv_ki",oovv_id,dim4_root,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('oovv_ki creation failed!')
                  ierr=exatns_tensor_init(oovv_ki,'ZERO')
                  ierr=exatns_tensor_contract("I(k,i,a,b)+=T(p,q,a,b)*P(k,p,i,q)",oovv_ki,int_t%oovv,P_ki)
                  if(ierr.ne.EXA_SUCCESS) call quit('oovv_ki extraction failed!')

                  ierr=exatns_tensor_destroy(P_ki)

                  !--------------------------- start computation of triples

                  ierr=exatns_tensor_create(w,'WW',w_id,w_sub,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: creation of wa')
                  ierr=exatns_tensor_create(y,'Y',w_id,w_sub,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: creation of y')
                  ierr=exatns_tensor_create(temp,'Temp',w_id,w_sub,EXA_DATA_KIND_C8)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: creation of temp')

                  ! contracting terms
                  ierr=exatns_tensor_init(temp,'ZERO')
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: init of wa')
                  ! term g^{abe}_i t^{ec}_jk
                  ierr=exatns_tensor_contract("W(a,b,c,i,j,k)+=V+(e,i,a,b)*T(e,c,j,k)",temp,vovv_i,t2_jk)
                  ierr=exatns_tensor_contract("W(a,b,c,i,j,k)+=V+(e,j,a,b)*T(e,c,k,i)",temp,vovv_j,t2_ki)
                  ierr=exatns_tensor_contract("W(a,b,c,i,j,k)+=V+(e,k,a,b)*T(e,c,i,j)",temp,vovv_k,t2_ij)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 11')
                  ! term g^{abe}_i t^{ec}_jk
                  ierr=exatns_tensor_contract("W(a,b,c,i,j,k)+=V+(i,j,l,a)*T(b,c,k,l)",temp,ooov_ij,t2_k)
                  ierr=exatns_tensor_contract("W(a,b,c,i,j,k)+=V+(j,k,l,a)*T(b,c,i,l)",temp,ooov_jk,t2_i) 
                  ierr=exatns_tensor_contract("W(a,b,c,i,j,k)+=V+(k,i,l,a)*T(b,c,j,l)",temp,ooov_ki,t2_j)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 12')

                  ! antisymmetrization
                  ierr=exatns_tensor_init(w,'ZERO')
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: init of w')
                  ierr=exatns_tensor_contract("C(a,b,c,i,j,k)+=W(a,b,c,i,j,k)*O()",w,temp,one_tensor) 
                  ierr=exatns_tensor_contract("C(a,b,c,i,j,k)+=W(b,c,a,i,j,k)*O()",w,temp,one_tensor)
                  ierr=exatns_tensor_contract("C(a,b,c,i,j,k)+=W(c,a,b,i,j,k)*O()",w,temp,one_tensor)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: copy 1')
                  if(exa_input%print_level > 8) write(*,*) "W[+T ] =",print_tensornorm2(w)

                  ! copy
                  ierr=exatns_tensor_init(y,'ZERO')
                  ierr=exatns_tensor_contract("D(a,b,c,i,j,k)+=W(a,b,c,i,j,k)*O()",y,w,one_tensor)
                  !Denominate
                  ierr=exatns_tensor_transform(y,denom)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: denominate')
                  if(exa_input%print_level > 8) write(*,*) "W[/D ] =",print_tensornorm2(y)

                  ! add to energy
                  ierr=exatns_tensor_init(res_tensor,ZERO)
                  ierr=exatns_tensor_contract("R()+=D(a,b,c,i,j,k)*W+(a,b,c,i,j,k)",res_tensor,y,w) 
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 13')
                  ierr=exatns_tensor_get_scalar(res_tensor,res_val)
                  t_energy(1)=t_energy(1) + res_val*prefact

                  ! Form Y-intermediate and calculate 5th order (T) energy contribution
                  ierr=exatns_tensor_init(temp,'ZERO')
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V+(i,j,a,b)*T(c,k)",temp,oovv_ij,t1_k)
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V+(j,k,a,b)*T(c,i)",temp,oovv_jk,t1_i)
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V+(k,i,a,b)*T(c,j)",temp,oovv_ki,t1_j)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 21')

                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V(a,b,i,j)*T(k,c)",temp,t2_ij,f_k)
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V(a,b,j,k)*T(i,c)",temp,t2_jk,f_i)
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V(a,b,k,i)*T(j,c)",temp,t2_ki,f_j)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 22')

                  ! symmetrization and copy
                  ierr=exatns_tensor_init(y,'ZERO')
                  ierr=exatns_tensor_contract("Y(a,b,c,i,j,k)+=T(a,b,c,i,j,k)*O()",y,temp,one_tensor)
                  ierr=exatns_tensor_contract("Y(a,b,c,i,j,k)+=T(b,c,a,i,j,k)*O()",y,temp,one_tensor)
                  ierr=exatns_tensor_contract("Y(a,b,c,i,j,k)+=T(c,a,b,i,j,k)*O()",y,temp,one_tensor)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: copy 2')
                  if(exa_input%print_level > 8) write(*,*) "Y[(T)] =",print_tensornorm2(y)
                  ierr=exatns_tensor_transform(y,denom)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: denominate')

                  !write value
                  ierr=exatns_tensor_init(res_tensor,ZERO)
                  ierr=exatns_tensor_contract("R()+=Y(a,b,c,i,j,k)*W+(a,b,c,i,j,k)",res_tensor,y,w) 
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 23')
                  ierr=exatns_tensor_get_scalar(res_tensor,res_val)
                  t_energy(2)=t_energy(2) - res_val*prefact

                  !Form Y-intermediate and calculate 5th order -T energy contribution
                  ierr=exatns_tensor_init(temp,'ZERO')
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V(a,b,i,j)*T(c,k)",temp,t2_ij,t1_k)
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V(a,b,j,k)*T(c,i)",temp,t2_jk,t1_i)
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V(a,b,k,i)*T(c,j)",temp,t2_ki,t1_j)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 31')
                  dim4_root=vvoo_root
                  dim4_root(3)=subsp(i)
                  dim4_root(4)=subsp(j)
                  ierr=exatns_tensor_create(temp2,'Temp2',vvoo_id,dim4_root,EXA_DATA_KIND_C8)
                  ierr=exatns_tensor_init(temp2,'ZERO')
                  ierr=exatns_tensor_contract("V(a,b,i,j)+=K(a,i)*L(b,j)",temp2,t1_i,t1_j)
                  ierr=exatns_tensor_contract("V(a,b,i,j)+=K(a,j)*L(b,i)",temp2,t1_j,t1_i,MINUS_ONE)
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V(a,b,i,j)*W(c,k)",temp,temp2,t1_k,ONE_THIRD)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 32')
                  ierr=exatns_tensor_destroy(temp2)
                  dim4_root(3)=subsp(j)
                  dim4_root(4)=subsp(k)
                  ierr=exatns_tensor_create(temp2,'Temp2',vvoo_id,dim4_root,EXA_DATA_KIND_C8)
                  ierr=exatns_tensor_init(temp2,'ZERO')
                  ierr=exatns_tensor_contract("V(a,b,j,k)+=K(a,j)*L(b,k)",temp2,t1_j,t1_k)
                  ierr=exatns_tensor_contract("V(a,b,j,k)+=K(a,k)*L(b,j)",temp2,t1_k,t1_j,MINUS_ONE)
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V(a,b,j,k)*W(c,i)",temp,temp2,t1_i,ONE_THIRD)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 32')
                  ierr=exatns_tensor_destroy(temp2)
                  dim4_root(3)=subsp(k)
                  dim4_root(4)=subsp(i)
                  ierr=exatns_tensor_create(temp2,'Temp2',vvoo_id,dim4_root,EXA_DATA_KIND_C8)
                  ierr=exatns_tensor_init(temp2,'ZERO')
                  ierr=exatns_tensor_contract("V(a,b,k,i)+=K(a,k)*L(b,i)",temp2,t1_k,t1_i)
                  ierr=exatns_tensor_contract("V(a,b,k,i)+=K(a,i)*L(b,k)",temp2,t1_i,t1_k,MINUS_ONE)
                  ierr=exatns_tensor_contract("T(a,b,c,i,j,k)+=V(a,b,k,i)*W(c,j)",temp,temp2,t1_j,ONE_THIRD)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 32')
                  ierr=exatns_tensor_destroy(temp2)

                  ! symmetrization and copy
                  ierr=exatns_tensor_init(y,'ZERO')
                  ierr=exatns_tensor_contract("Y(a,b,c,i,j,k)+=T(a,b,c,i,j,k)*O()",y,temp,one_tensor)
                  ierr=exatns_tensor_contract("Y(a,b,c,i,j,k)+=T(b,c,a,i,j,k)*O()",y,temp,one_tensor)
                  ierr=exatns_tensor_contract("Y(a,b,c,i,j,k)+=T(c,a,b,i,j,k)*O()",y,temp,one_tensor)
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: copy 3')
                  if(exa_input%print_level > 8) write(*,*) "Y[-T ] =",print_tensornorm2(y)

                  !write value
                  ierr=exatns_tensor_init(res_tensor,ZERO)
                  ierr=exatns_tensor_contract("R()+=Y(a,b,c,i,j,k)*W+(a,b,c,i,j,k)",res_tensor,y,w) 
                  if(ierr.ne.EXA_SUCCESS) call quit('ERROR in t3_energy: contraction 35')
                  ierr=exatns_tensor_get_scalar(res_tensor,res_val)
                  t_energy(3)=t_energy(3) - res_val*prefact

                  if (exa_input%print_level > 6) write(*,*) 't3 = ', t_energy

                  ierr=exatns_tensor_destroy(w)
                  ierr=exatns_tensor_destroy(y)
                  ierr=exatns_tensor_destroy(temp)

                  ierr=exatns_tensor_destroy(t2_ij)
                  ierr=exatns_tensor_destroy(ooov_ij)
                  ierr=exatns_tensor_destroy(oovv_ij)

                  ierr=exatns_tensor_destroy(t2_ki)
                  ierr=exatns_tensor_destroy(ooov_ki)
                  ierr=exatns_tensor_destroy(oovv_ki)

                  ierr=exatns_tensor_destroy(vovv_i)
                  ierr=exatns_tensor_destroy(t2_i)
                  ierr=exatns_tensor_destroy(t1_i)
                  ierr=exatns_tensor_destroy(f_i)
                  ierr=exatns_tensor_destroy(P_i)

                  if (exa_input%moint_scheme==5) exit

                end do

                ierr=exatns_tensor_destroy(t2_jk)
                ierr=exatns_tensor_destroy(ooov_jk)
                ierr=exatns_tensor_destroy(oovv_jk)

                ierr=exatns_tensor_destroy(vovv_j)
                ierr=exatns_tensor_destroy(t2_j)
                ierr=exatns_tensor_destroy(t1_j)
                ierr=exatns_tensor_destroy(f_j)
                ierr=exatns_tensor_destroy(P_j)

                if (exa_input%moint_scheme==5) exit

              end do

              ierr=exatns_tensor_destroy(vovv_k)
              ierr=exatns_tensor_destroy(t2_k)
              ierr=exatns_tensor_destroy(t1_k)
              ierr=exatns_tensor_destroy(f_k)
              ierr=exatns_tensor_destroy(P_k)
               
              if (exa_input%print_level>-1) then
                write (info_str,'(A9,I4,A4,I4,A9,I10)') '  block  ',k,' of ',n_subsp,' done, n=',k**2
                call print_date(info_str)
              end if
              if (exa_input%print_level > 0) write (*,*)  'after this block, t3 = ',t_energy

              if (exa_input%moint_scheme==5) exit

            end do

            ! We looped over all (abc) combinations, scale with the factor of 1/6 needed to go to the unique ones.
            t_energy = t_energy / SIX
            
            !clean up
            ierr=exatns_tensor_destroy(res_tensor)
            deallocate(subsp)
            
            if (exa_input%print_level >4) call print_date('leaving t3_energy_split')

        end subroutine t3_energy_split


        subroutine t3_energy_delta(t_energy,comm_t,int_t,f_delta,denom3,mo_occ,print_level,beta_occ)
            
            use exacorr_global
            use exacorr_tensor_methods
            use exacorr_ao_to_mo, only : print_tensornorm2

            real(8), intent(out)                    :: t_energy(3)
            type(exatns_comm_tens)                  :: comm_t
            type(exatns_intg_tens)                  :: int_t
            type(delta_t)                           :: f_delta  ! fill delta matrix
            type(denom3_t)                          :: denom3
            integer                                 :: mo_occ(:)
            integer                                 :: print_level
            logical                                 :: beta_occ

            integer(INTL),dimension(2) :: ov_root, vo_root, subspace_id2
            integer(INTD),dimension(2) :: ov_id, vo_id
            integer(INTL),dimension(4) :: ooov_root, vovv_root, vvoo_root, oovv_root, subspace_id4
            integer(INTD),dimension(4) :: ooov_id,   vovv_id,   vvoo_id,   oovv_id

            integer(INTD)                :: ierr,rank
            integer(INTD)                :: idocc, idvir
            integer(INTL)                :: nocc, nvir
            integer(INTL),dimension(3)   :: abc_root
            integer(INTD),dimension(3)   :: abc_id
            type(tens_rcrsv_t)           :: res_tensor
            complex(8)                   :: result_val, ONE_THIRD
            type(tens_rcrsv_t)           :: w_tensor, y_tensor, we_tensor, temp
            type(tens_rcrsv_t)           :: t1_i, t1_j, t1_k
            type(tens_rcrsv_t)           :: t2_i, t2_j, t2_k
            type(tens_rcrsv_t)           :: t2_jk, t2_ki, t2_ij
            type(tens_rcrsv_t)           :: vovv_i, vovv_j, vovv_k
            type(tens_rcrsv_t)           :: ooov_jk, ooov_ki, ooov_ij
            type(tens_rcrsv_t)           :: oovv_jk, oovv_ki, oovv_ij
            type(tens_rcrsv_t)           :: f_i, f_j, f_k
            type(tens_rcrsv_t)           :: delta_i, delta_j, delta_k, delta_2
            integer(INTL)                :: delta1(1)
            real(8), allocatable         :: eps_occ(:)
            real(8)                      :: eps_ijk
            !complex(8), allocatable      :: y(:,:,:) ! used in form_y_intermediate, but allocated here to have fewer alloc/dealloc calls
            integer                      :: i,j,k
            character(len=40)            :: info_str

            if (print_level > 4) call print_date('entering t3_energy_delta')

            if (beta_occ) call quit('t3_energy_delta: only closed shell')

            ONE_THIRD=ONE/3.00D0

            call comm_t%fov%get_dims(ov_root,rank,ierr)
            call comm_t%fov%get_space_ids(ov_id,subspace_id2,ierr)
            vo_root=ov_root(1)
            vo_id=ov_id(1)
            vo_root(1)=ov_root(2)
            vo_id(1)=ov_id(2)

            call int_t%ooov%get_dims(ooov_root,rank,ierr)
            call int_t%ooov%get_space_ids(ooov_id,subspace_id4,ierr)
            call int_t%vovv%get_dims(vovv_root,rank,ierr)
            call int_t%vovv%get_space_ids(vovv_id,subspace_id4,ierr)
            call int_t%oovv%get_dims(oovv_root,rank,ierr)
            call int_t%oovv%get_space_ids(oovv_id,subspace_id4,ierr)
            vvoo_root(1:2)=oovv_root(3:4)
            vvoo_root(3:4)=oovv_root(1:2)
            vvoo_id(1:2)=oovv_id(3:4)
            vvoo_id(3:4)=oovv_id(1:2)

            t_energy = ZERO
            nocc=vvoo_root(3)
            nvir=vvoo_root(1)
            idocc=vvoo_id(3)
            idvir=vvoo_id(1)

!           get orbital energies
            allocate(eps_occ(nocc))
            call get_eps(eps_occ,mo_occ)

!           Initialize scalars that are to be used as tensors in contractions
            ierr=exatns_tensor_create(res_tensor,"res_tensor",EXA_DATA_KIND_C8)

!           Allocate tensors of dimension nvir**3 to hold intermediates
            abc_root=nvir
            abc_id=idvir
            ierr=exatns_tensor_create(w_tensor,"w_tensor",abc_id,abc_root,EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(we_tensor,"we_tensor",abc_id,abc_root,EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(y_tensor,"y_tensor",abc_id,abc_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(temp,"temp",abc_id,abc_root,EXA_DATA_KIND_C8) 
            

            ierr=exatns_tensor_create(vovv_i,"vovv_i",abc_id,abc_root,EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(vovv_j,"vovv_j",abc_id,abc_root,EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(vovv_k,"vovv_k",abc_id,abc_root,EXA_DATA_KIND_C8) 

            abc_root(3)=nocc
            abc_id(3)=idocc
            ierr=exatns_tensor_create(t2_i,"t2_i",abc_id,abc_root,EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(t2_j,"t2_j",abc_id,abc_root,EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(t2_k,"t2_k",abc_id,abc_root,EXA_DATA_KIND_C8) 

            ierr=exatns_tensor_create(t2_jk,"t2_jk",vvoo_id(1:2),vvoo_root(1:2),EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(t2_ki,"t2_ki",vvoo_id(1:2),vvoo_root(1:2),EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(t2_ij,"t2_ij",vvoo_id(1:2),vvoo_root(1:2),EXA_DATA_KIND_C8) 

            ierr=exatns_tensor_create(ooov_ij,"ooov_ij",ooov_id(3:4),ooov_root(3:4),EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(ooov_jk,"ooov_jk",ooov_id(3:4),ooov_root(3:4),EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(ooov_ki,"ooov_ki",ooov_id(3:4),ooov_root(3:4),EXA_DATA_KIND_C8)

            ierr=exatns_tensor_create(oovv_jk,"oovv_jk",vvoo_id(1:2),vvoo_root(1:2),EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(oovv_ki,"oovv_ki",vvoo_id(1:2),vvoo_root(1:2),EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(oovv_ij,"oovv_ij",vvoo_id(1:2),vvoo_root(1:2),EXA_DATA_KIND_C8) 

            ierr=exatns_tensor_create(t1_i,"t1_i",vvoo_id(1:1),vvoo_root(1:1),EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(t1_j,"t1_j",vvoo_id(1:1),vvoo_root(1:1),EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(t1_k,"t1_k",vvoo_id(1:1),vvoo_root(1:1),EXA_DATA_KIND_C8)

            ierr=exatns_tensor_create(f_i,"f_i",vvoo_id(1:1),vvoo_root(1:1),EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(f_j,"f_j",vvoo_id(1:1),vvoo_root(1:1),EXA_DATA_KIND_C8) 
            ierr=exatns_tensor_create(f_k,"f_k",vvoo_id(1:1),vvoo_root(1:1),EXA_DATA_KIND_C8)    
      
            ierr=exatns_tensor_create(delta_i,"delta_i",ooov_id(1:1),ooov_root(1:1),EXA_DATA_KIND_C8)  
            ierr=exatns_tensor_create(delta_j,"delta_j",ooov_id(1:1),ooov_root(1:1),EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(delta_k,"delta_k",ooov_id(1:1),ooov_root(1:1),EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(delta_2,"delta_2",ooov_id(1:2),ooov_root(1:2),EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(delta_i,'ZERO')
            ierr=exatns_tensor_init(delta_j,'ZERO')
            ierr=exatns_tensor_init(delta_k,'ZERO')
            
   
            do k = 3, nocc
              delta1(1)=k
              call f_delta%reset(delta1,1)
              ierr=exatns_tensor_transform(delta_k,f_delta)
              if (print_level.gt.12) ierr=exatns_tensor_traverse(delta_k,'PrintTensor')
              ierr=exatns_tensor_init(vovv_k,'ZERO')
              ierr=exatns_tensor_contract("I(c,a,b)+=V(c,k,a,b)*D(k)",vovv_k,int_t%vovv,delta_k)
              if(ierr.ne.EXA_SUCCESS) call quit('vovv_k extraction failed!')
              ierr=exatns_tensor_init(t2_k,'ZERO')
              ierr=exatns_tensor_contract("I(a,b,i)+=T(a,b,k,i)*D(k)",t2_k,comm_t%t2,delta_k)
              if(ierr.ne.EXA_SUCCESS) call quit('t2_k extraction failed!')

              ierr=exatns_tensor_init(t1_k,'ZERO')
              ierr=exatns_tensor_contract("I(a)+=T(a,k)*D(k)",t1_k,comm_t%t1,delta_k)
              if(ierr.ne.EXA_SUCCESS) call quit('t1_k extraction failed!')
              ierr=exatns_tensor_init(f_k,'ZERO')
              ierr=exatns_tensor_contract("I(a)+=T(k,a)*D(k)",f_k,comm_t%fov,delta_k)
              if(ierr.ne.EXA_SUCCESS) call quit('f_k extraction failed!')

              do j = 2, k-1
                delta1(1)=j
                call f_delta%reset(delta1,1)
                ierr=exatns_tensor_transform(delta_j,f_delta)
                if (print_level.gt.12) ierr=exatns_tensor_traverse(delta_j,'PrintTensor')
                ierr=exatns_tensor_init(vovv_j,'ZERO')
                ierr=exatns_tensor_contract("I(c,a,b)+=V(c,j,a,b)*D(j)",vovv_j,int_t%vovv,delta_j)
                if(ierr.ne.EXA_SUCCESS) call quit('vovv_j extraction failed!')
                ierr=exatns_tensor_init(t2_j,'ZERO')
                ierr=exatns_tensor_contract("I(a,b,i)+=T(a,b,j,i)*D(j)",t2_j,comm_t%t2,delta_j)
                if(ierr.ne.EXA_SUCCESS) call quit('t2_j extraction failed!')

                ierr=exatns_tensor_init(t1_j,'ZERO')
                ierr=exatns_tensor_contract("I(a)+=T(a,j)*D(j)",t1_j,comm_t%t1,delta_j)
                if(ierr.ne.EXA_SUCCESS) call quit('t1_j extraction failed!')
                ierr=exatns_tensor_init(f_j,'ZERO')
                ierr=exatns_tensor_contract("I(a)+=T(j,a)*D(j)",f_j,comm_t%fov,delta_j)
                if(ierr.ne.EXA_SUCCESS) call quit('f_j extraction failed!')

                ierr=exatns_tensor_init(delta_2,'ZERO')
                ierr=exatns_tensor_contract("D(j,k)+=A(j)*B(k)",delta_2,delta_j,delta_k)
                ierr=exatns_tensor_init(t2_jk,'ZERO')
                ierr=exatns_tensor_contract("I(a,b)+=T(a,b,j,k)*D(j,k)",t2_jk,comm_t%t2,delta_2)
                if(ierr.ne.EXA_SUCCESS) call quit('t2_jk extraction failed!')
                ierr=exatns_tensor_init(ooov_jk,'ZERO')
                ierr=exatns_tensor_contract("I(i,a)+=V(j,k,i,a)*D(j,k)",ooov_jk,int_t%ooov,delta_2)
                if(ierr.ne.EXA_SUCCESS) call quit('ooov_jk extraction failed!')
                ierr=exatns_tensor_init(oovv_jk,'ZERO')
                ierr=exatns_tensor_contract("I(a,b)+=V(j,k,a,b)*D(j,k)",oovv_jk,int_t%oovv,delta_2)
                if(ierr.ne.EXA_SUCCESS) call quit('oovv_jk extraction failed!')

                do i = 1, j-1
                  delta1(1)=i
                  call f_delta%reset(delta1,1)
                  ierr=exatns_tensor_transform(delta_i,f_delta)
                  if (print_level.gt.12) ierr=exatns_tensor_traverse(delta_i,'PrintTensor')
                  ierr=exatns_tensor_init(vovv_i,'ZERO')
                  ierr=exatns_tensor_contract("I(c,a,b)+=V(c,i,a,b)*D(i)",vovv_i,int_t%vovv,delta_i)
                  if(ierr.ne.EXA_SUCCESS) call quit('vovv_i extraction failed!')
                  ierr=exatns_tensor_init(t2_i,'ZERO')
                  ierr=exatns_tensor_contract("I(a,b,j)+=T(a,b,i,j)*D(i)",t2_i,comm_t%t2,delta_i)
                  if(ierr.ne.EXA_SUCCESS) call quit('t2_i extraction failed!')

                  ierr=exatns_tensor_init(t1_i,'ZERO')
                  ierr=exatns_tensor_contract("I(a)+=T(a,i)*D(i)",t1_i,comm_t%t1,delta_i)
                  if(ierr.ne.EXA_SUCCESS) call quit('t1_i extraction failed!')
                  ierr=exatns_tensor_init(f_i,'ZERO')
                  ierr=exatns_tensor_contract("I(a)+=T(i,a)*D(i)",f_i,comm_t%fov,delta_i)
                  if(ierr.ne.EXA_SUCCESS) call quit('f_i extraction failed!')

                  ierr=exatns_tensor_init(delta_2,'ZERO')
                  ierr=exatns_tensor_contract("D(i,j)+=A(i)*B(j)",delta_2,delta_i,delta_j)
                  ierr=exatns_tensor_init(t2_ij,'ZERO')
                  ierr=exatns_tensor_contract("I(a,b)+=T(a,b,i,j)*D(i,j)",t2_ij,comm_t%t2,delta_2)
                  if(ierr.ne.EXA_SUCCESS) call quit('t2_ij extraction failed!')
                  ierr=exatns_tensor_init(ooov_ij,'ZERO')
                  ierr=exatns_tensor_contract("I(k,a)+=V(i,j,k,a)*D(i,j)",ooov_ij,int_t%ooov,delta_2)
                  if(ierr.ne.EXA_SUCCESS) call quit('ooov_ij extraction failed!')
                  ierr=exatns_tensor_init(oovv_ij,'ZERO')
                  ierr=exatns_tensor_contract("I(a,b)+=V(i,j,a,b)*D(i,j)",oovv_ij,int_t%oovv,delta_2)
                  if(ierr.ne.EXA_SUCCESS) call quit('oovv_ij extraction failed!')

                  ierr=exatns_tensor_init(delta_2,'ZERO')
                  ierr=exatns_tensor_contract("D(k,i)+=A(k)*B(i)",delta_2,delta_k,delta_i)
                  ierr=exatns_tensor_init(t2_ki,'ZERO')
                  ierr=exatns_tensor_contract("I(a,b)+=T(a,b,k,i)*D(k,i)",t2_ki,comm_t%t2,delta_2)
                  if(ierr.ne.EXA_SUCCESS) call quit('t2_ki extraction failed!')
                  ierr=exatns_tensor_init(ooov_ki,'ZERO')
                  ierr=exatns_tensor_contract("I(j,a)+=V(k,i,j,a)*D(k,i)",ooov_ki,int_t%ooov,delta_2)
                  if(ierr.ne.EXA_SUCCESS) call quit('ooov_ki extraction failed!')
                  ierr=exatns_tensor_init(oovv_ki,'ZERO')
                  ierr=exatns_tensor_contract("I(a,b)+=V(k,i,a,b)*D(k,i)",oovv_ki,int_t%oovv,delta_2)
                  if(ierr.ne.EXA_SUCCESS) call quit('oovv_ki extraction failed!')

                  ! Form w-intermediate, first without symmetrization and minus sign in the auxilliary tensor array we
                  ierr=exatns_tensor_init(temp,'ZERO')
                  ! For this term we use the permutation symmetry <vv||ov> = - <vo||vv>*, with vvov(abe:i) = vovv*(eab:i)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V+(e,a,b)*T(e,c)",temp,vovv_i,t2_jk) ! o(3)v(4)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V+(e,a,b)*T(e,c)",temp,vovv_j,t2_ki) ! o(3)v(4)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V+(e,a,b)*T(e,c)",temp,vovv_k,t2_ij) ! o(3)v(4)
                  ! For this term we use the permutation symmetry <vo||oo> = - <oo||ov>*, with vooo(al:ij) = ooov*(la:ij)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V+(l,a)*T(b,c,l)",temp,ooov_ij,t2_k) ! o(4)v(3)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V+(l,a)*T(b,c,l)",temp,ooov_jk,t2_i) ! o(4)v(3)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V+(l,a)*T(b,c,l)",temp,ooov_ki,t2_j) ! o(4)v(3)
                  if(print_level > 10) write(*,*) "W[sym] =",print_tensornorm2(temp)

                  ! Copy into w array to symmetrize and add the minus sign omitted in the above contractions
                  ierr=exatns_tensor_init(w_tensor,'ZERO')
                  ! symmetrization via contraction as addition with permutation is not implemented
                  ierr=exatns_tensor_contract("W(a,b,c)+=E(a,b,c)*O()",w_tensor,temp,one_tensor)
                  ierr=exatns_tensor_contract("W(a,b,c)+=E(b,c,a)*O()",w_tensor,temp,one_tensor)
                  ierr=exatns_tensor_contract("W(a,b,c)+=E(c,a,b)*O()",w_tensor,temp,one_tensor)
                  if(print_level > 10) write(*,*) "W[+T ] =",print_tensornorm2(w_tensor)

                  eps_ijk = eps_occ(i) + eps_occ(j) + eps_occ(k)
                  call denom3%reset(eps_ijk)
                  ierr=exatns_tensor_init(we_tensor,'ZERO')
                  ierr=exatns_tensor_contract("E(a,b,c)+=W(a,b,c)*O()",we_tensor,w_tensor,one_tensor)
                  ierr=exatns_tensor_transform(we_tensor,denom3)

                  ierr=exatns_tensor_init(res_tensor,ZERO)
                  ierr=exatns_tensor_contract("R()+=W(a,b,c)*E+(a,b,c)",res_tensor,w_tensor,we_tensor)
                  ierr=exatns_tensor_get_scalar(res_tensor,result_val)
                  t_energy(1)=t_energy(1) + result_val


                  !Form Y-intermediate and calculate 5th order -T energy contribution
                  ierr=exatns_tensor_init(temp,'ZERO')
                  ierr=exatns_tensor_contract("E(a,b,c)+=V(a,b)*T(c)",temp,t2_ij,t1_k)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V(a,b)*T(c)",temp,t2_jk,t1_i)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V(a,b)*T(c)",temp,t2_ki,t1_j)
                  ierr=exatns_tensor_init(t2_ij,'ZERO')
                  ierr=exatns_tensor_contract("E(a,b)+=V(a)*T(b)",t2_ij,t1_i,t1_j)
                  ierr=exatns_tensor_contract("E(a,b)+=V(a)*T(b)",t2_ij,t1_j,t1_i,MINUS_ONE)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V(a,b)*T(c)",temp,t2_ij,t1_k,ONE_THIRD)
                  ierr=exatns_tensor_init(t2_ij,'ZERO')
                  ierr=exatns_tensor_contract("E(a,b)+=V(a)*T(b)",t2_ij,t1_j,t1_k)
                  ierr=exatns_tensor_contract("E(a,b)+=V(a)*T(b)",t2_ij,t1_k,t1_j,MINUS_ONE)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V(a,b)*T(c)",temp,t2_ij,t1_i,ONE_THIRD)
                  ierr=exatns_tensor_init(t2_ki,'ZERO')
                  ierr=exatns_tensor_contract("E(a,b)+=V(a)*T(b)",t2_ki,t1_k,t1_i)
                  ierr=exatns_tensor_contract("E(a,b)+=V(a)*T(b)",t2_ki,t1_i,t1_k,MINUS_ONE)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V(a,b)*T(c)",temp,t2_ki,t1_j,ONE_THIRD)

                  ! symmetrization and copy
                  ierr=exatns_tensor_init(y_tensor,'ZERO')
                  ierr=exatns_tensor_contract("Y(a,b,c)+=E(a,b,c)*O()",y_tensor,temp,one_tensor)
                  ierr=exatns_tensor_contract("Y(a,b,c)+=E(b,c,a)*O()",y_tensor,temp,one_tensor)
                  ierr=exatns_tensor_contract("Y(a,b,c)+=E(c,a,b)*O()",y_tensor,temp,one_tensor)
                  if(print_level > 10) write(*,*) "Y[-T ] =",print_tensornorm2(y_tensor)

                  !write value
                  ierr=exatns_tensor_init(res_tensor,ZERO)
                  ierr=exatns_tensor_contract("R()+=Y(a,b,c)*W+(a,b,c)",res_tensor,y_tensor,w_tensor)
                  ierr=exatns_tensor_get_scalar(res_tensor,result_val)
                  t_energy(3)=t_energy(3) - result_val

                  ! Form Y-intermediate and calculate 5th order (T) energy contribution
                  ierr=exatns_tensor_init(temp,'ZERO')
                  ierr=exatns_tensor_contract("E(a,b,c)+=V+(a,b)*T(c)",temp,oovv_ij,t1_k)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V+(a,b)*T(c)",temp,oovv_jk,t1_i)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V+(a,b)*T(c)",temp,oovv_ki,t1_j)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V(a,b)*T(c)",temp,t2_ij,f_k)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V(a,b)*T(c)",temp,t2_jk,f_i)
                  ierr=exatns_tensor_contract("E(a,b,c)+=V(a,b)*T(c)",temp,t2_ki,f_j)

                  ! symmetrization and copy
                  ierr=exatns_tensor_init(y_tensor,'ZERO')
                  ierr=exatns_tensor_contract("Y(a,b,c)+=E(a,b,c)*O()",y_tensor,temp,one_tensor)
                  ierr=exatns_tensor_contract("Y(a,b,c)+=E(b,c,a)*O()",y_tensor,temp,one_tensor)
                  ierr=exatns_tensor_contract("Y(a,b,c)+=E(c,a,b)*O()",y_tensor,temp,one_tensor)
                  if(print_level > 10) write(*,*) "Y[(T)] =",print_tensornorm2(y_tensor)

                  !write value
                  ierr=exatns_tensor_init(res_tensor,ZERO)
                  ierr=exatns_tensor_contract("R()+=Y(a,b,c)*W+(a,b,c)",res_tensor,y_tensor,we_tensor) 
                  ierr=exatns_tensor_get_scalar(res_tensor,result_val)
                  t_energy(2)=t_energy(2) - result_val
                  
                  if(print_level > 10) write(*,*) "W[/D ]  =",print_tensornorm2(we_tensor)

                end do
              end do
              if (print_level>0) then
                write (info_str,'(A9,I4,A4,I4,A9,I10)') '  block  ',k,' of ',nocc,' done, n=',k**2
                call print_date(info_str)
              end if
              if (print_level >2) write (*,*)  't3_energy:',t_energy
            end do

!           We looped over all (abc) combinations, scale with the factor of 1/6 needed to go to the unique ones.
            t_energy = t_energy / SIX

!           Destruct intermediate tensors
            deallocate(eps_occ)
            !deallocate(y)
            ierr=exatns_tensor_destroy(res_tensor)
            ierr=exatns_tensor_destroy(t2_i)
            ierr=exatns_tensor_destroy(t2_j)
            ierr=exatns_tensor_destroy(t2_k)
            ierr=exatns_tensor_destroy(t2_ij)
            ierr=exatns_tensor_destroy(t2_jk)
            ierr=exatns_tensor_destroy(t2_ki)
            ierr=exatns_tensor_destroy(vovv_i)
            ierr=exatns_tensor_destroy(vovv_j)
            ierr=exatns_tensor_destroy(vovv_k)
            ierr=exatns_tensor_destroy(ooov_ij)
            ierr=exatns_tensor_destroy(ooov_jk)
            ierr=exatns_tensor_destroy(ooov_ki)
            ierr=exatns_tensor_destroy(oovv_ij)
            ierr=exatns_tensor_destroy(oovv_jk)
            ierr=exatns_tensor_destroy(oovv_ki)
            ierr=exatns_tensor_destroy(t1_i)
            ierr=exatns_tensor_destroy(t1_j)
            ierr=exatns_tensor_destroy(t1_k)
            ierr=exatns_tensor_destroy(f_i)
            ierr=exatns_tensor_destroy(f_j)
            ierr=exatns_tensor_destroy(f_k)
            ierr=exatns_tensor_destroy(delta_i)
            ierr=exatns_tensor_destroy(delta_j)
            ierr=exatns_tensor_destroy(delta_k)
            ierr=exatns_tensor_destroy(delta_2)
            ierr=exatns_tensor_destroy(we_tensor)
            ierr=exatns_tensor_destroy(w_tensor)
            ierr=exatns_tensor_destroy(y_tensor)
            ierr=exatns_tensor_destroy(temp)

            if (print_level > 4) call print_date('leaving t3_energy_delta')

        end subroutine t3_energy_delta

        subroutine get_subtensor (tensor,stensor,f_delta,mode,i,j)

            use exacorr_tensor_methods
            use exacorr_ao_to_mo, only : print_tensornorm2

!        i: passive index #1
!        j: passive index #2

            type(tens_rcrsv_t), intent(inout) :: tensor, stensor
            type(delta_t)                     :: f_delta  ! fill delta matrix
            integer, intent(in)               :: i, mode
            integer, intent(in), optional     :: j


            integer(INTD)             :: ierr, t_rank, s_rank
            integer(INTL),allocatable :: t_dim(:),s_dim(:),subspace_id(:)
            integer(INTD),allocatable :: space_id(:)
            type(tens_rcrsv_t)        :: Delta
            integer(INTL)             :: delta2(2),delta1(1)

!           get rank of tensors
            t_rank=tensor%get_rank(ierr)
            s_rank=stensor%get_rank(ierr)

!           check ranks
            if (t_rank .ne. 4) call quit('Error in get_subtensor: Only extraction from 4-dimensional tensors possible')
            if (present(j)) then
                if (s_rank .ne. 2) call quit('Error in get_subtensor: Subtensor to many dimensions')
            else
                if (s_rank .ne. 3) call quit('Error in get_subtensor: Subtensor to many dimensions')
            end if

!           get dimensions of tensors
            allocate(t_dim(t_rank))
            call tensor%get_dims(t_dim,t_rank,ierr)
            allocate(s_dim(s_rank))
            call stensor%get_dims(s_dim,s_rank,ierr)
!           get space id
            allocate(space_id(t_rank),subspace_id(t_rank))
            call tensor%get_space_ids(space_id,subspace_id,ierr)


!           check dims, create htensor
            if (mode==3) then
              if (t_dim(1).ne.s_dim(1) .or. t_dim(2).ne.s_dim(2) .or. t_dim(4).ne.s_dim(3)) then
                call quit('Error in get_subtensor: tensor dimension different')
              end if
              ierr=exatns_tensor_create(Delta,"Delta",space_id(3:3),t_dim(3:3),EXA_DATA_KIND_C8)
            else if (mode==2) then
              if (t_dim(1).ne.s_dim(1) .or. t_dim(3).ne.s_dim(2) .or. t_dim(4).ne.s_dim(3)) then
                call quit('Error in get_subtensor: tensor dimension different')
              end if
              ierr=exatns_tensor_create(Delta,"Delta",space_id(2:2),t_dim(2:2),EXA_DATA_KIND_C8)
            else if (mode==12) then
              if (t_dim(3).ne.s_dim(1) .or. t_dim(4).ne.s_dim(2)) then
                call quit('Error in get_subtensor: tensor dimension different')
              end if
              ierr=exatns_tensor_create(Delta,"Delta",space_id(1:2),t_dim(1:2),EXA_DATA_KIND_C8)
            else if (mode==34) then
              if (t_dim(1).ne.s_dim(1) .or. t_dim(2).ne.s_dim(2)) then
                call quit('Error in get_subtensor: tensor dimension different')
              end if
              ierr=exatns_tensor_create(Delta,"Delta",space_id(3:4),t_dim(3:4),EXA_DATA_KIND_C8)
            else
              call quit('Error in get_subtensor: mode not implemented')
            end if

            ! generate delta tensor
            if (present(j)) then
              delta2(1)=i
              delta2(2)=j
              call f_delta%reset(delta2,2)
            else
              delta1(1)=i
              call f_delta%reset(delta1,1)
            end if

            ierr=exatns_tensor_init(Delta,'ZERO')
            ierr=exatns_tensor_transform(Delta,f_delta)
            !ierr=exatns_tensor_traverse(Delta,'PrintTensor')

            !  extract element
            ierr=exatns_tensor_init(stensor,'ZERO')
            if (mode==3) then
              ierr=exatns_tensor_contract("S(p,q,s)=T(p,q,r,s)*D(r)",stensor,tensor,Delta)
            else if (mode==2) then
              ierr=exatns_tensor_contract("S(p,r,s)=T(p,q,r,s)*D(q)",stensor,tensor,Delta)
            else if (mode==12) then
              ierr=exatns_tensor_contract("S(r,s)=T(p,q,r,s)*D(p,q)",stensor,tensor,Delta)
            else if (mode==34) then
              ierr=exatns_tensor_contract("S(p,q)=T(p,q,r,s)*D(r,s)",stensor,tensor,Delta)
            end if

            ! cleanup
            ierr=exatns_tensor_destroy(Delta)
            deallocate(space_id)
            deallocate(subspace_id)
            deallocate(s_dim)
            deallocate(t_dim)

        end subroutine get_subtensor

#else
        implicit none

        public exacorr_cc_driver

       contains

        subroutine exacorr_cc_driver (exa_input)
         use exacorr_datatypes
         implicit none
         type(exacc_input), intent(in) :: exa_input
         write(*,*) 'called exatensor version of cc_driver in a serial version'
         write(*,*) 'this functionality only works with MPI, we need to stop here'
         stop ' activate MPI'
        end subroutine exacorr_cc_driver

#endif

end module exacorr_cc
