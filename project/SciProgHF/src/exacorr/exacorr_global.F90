module exacorr_global

!This module controls storage and interaction with the global data used in ExaCorr

        use exacorr_datatypes
        use exacorr_mo
        use, intrinsic:: ISO_C_BINDING

        implicit none
        private

!       Global data (should limit this to the minimum required to function)
        integer :: my_MPI_master=0 ! Define the MPI master proces (is always zero in the DIRAC MPI environment)
        type(basis_set_info_t) :: ao_basis ! information about the ao basis
        type(cmo):: cspinor_all ! Stores all MO-coefficients and energies obtained from the parent code
        type(one_el_t)  :: int_1el ! stores 1 electron integrals
        logical :: TALSH_ONLY = .TRUE. ! use single-node TALSH library instead of the full EXATENSOR library
        logical :: one_el_exist=.FALSE. ! logical if MRCONEE is available or not

!       Functions and routines to interact with these data
        public nfunctions
        public get_nao
        public get_basis_angular
        public get_nsolutions
        public get_talsh_only
        public set_talsh_only
        public initialize_global_data
        public delete_global_data
        public get_gtos
        public get_mo_coefficients
        public get_eps
        public exist_one_el
        public get_one_el
        public add_finitefield
        public get_ff_ind
        public get_ff_mat
        public get_shell_index
        public get_shell_size
        public get_shell_offset
        public make_spinor_list
        public get_scf

       contains

         subroutine set_talsh_only(new_value)
             logical, intent(in) :: new_value
             TALSH_ONLY = new_value
         end subroutine set_talsh_only

         logical function get_talsh_only()
             get_talsh_only = TALSH_ONLY
         end function get_talsh_only

         integer function get_nao()
             get_nao = ao_basis%nao
         end function get_nao

         integer function get_basis_angular()
             get_basis_angular = ao_basis%basis_angular
         end function get_basis_angular

         integer function get_nsolutions()
             get_nsolutions = cspinor_all%nmo
         end function get_nsolutions

         subroutine initialize_global_data

             use exacorr_respect
             use exacc_cfg
             use exacorr_utils, only : print_date

#if defined (VAR_MPI)
             use interface_to_mpi
             integer :: my_MPI_rank
#endif
             integer :: ierr, scf_file
             logical :: tobe

#if defined (VAR_MPI)
             call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
             if (my_MPI_rank /= my_MPI_master) then
!               Slaves go to the sync immediately and get all that is necessary from master
                call synchronize_global_data
                return
             end if
#endif
             ! read MO coefficients and ao basis set information
             inquire (file='DFCOEF',exist=tobe)
             if (tobe) then
                scf_file = 0
                call read_from_dirac (ao_basis,cspinor_all,'DFCOEF',ierr)
             else
                inquire (file='RSD_MOS',exist=tobe)
                if (tobe) then
                   scf_file = 1
                   call read_from_respect (ao_basis,cspinor_all,ierr)
                end if
             end if

             if (ierr == 0) then
                call print_date ("retrieved basis set information")
             else
                write(*,*) "Error code",ierr," while reading restart file"
                call quit ('Could not read restart file')
             end if

             ! read one electron integrals
             inquire (file='MRCONEE',exist=tobe)
             if (tobe) then
                call read_from_mrconee (int_1el, ierr)
                if (ierr == 0) then
                  call print_date("read one el. integrals from MRCONEE")
                  one_el_exist = .TRUE.
                else
                  write(*,*) "Error code",ierr," while reading MRCONEE"
                  call quit ('Could not read MRCONEE')
                end if
             else
               one_el_exist = .FALSE.
               if (exa_input%beta_occ .and. exa_input%level_shift.lt.1.0D-14) then
                 if (scf_file == 0) call quit ('One electron integrals required in the open shell case (MOLTRA) or LSHIFT')
               end if
             end if

             ! read property integrals (requires also readding of nsp from mrconee, so testing for this as well)
             if (exa_input%nff(1)>0) then
               inquire (file='MDPROP',exist=tobe)
             else
              exa_input%nff(2)=1
              tobe=.FALSE.
             end if
             if (tobe.and.one_el_exist) then
                call read_mdprop (int_1el, ierr)
                if (ierr == 0) then
                  call print_date("read property integrals from MDPROP")
                else
                  call print_date("could not read property integrals from MDPROP")
                  int_1el%n_prop = 0
                endif
             else
               int_1el%n_prop = 0
             end if

             call synchronize_global_data

         end subroutine initialize_global_data

         subroutine synchronize_global_data

#if defined (VAR_MPI)
             use interface_to_mpi
             integer :: nao, nshells, n_primitives, nsp
             integer :: ish
             integer :: my_MPI_rank

!            The next line is to be commented out when running the stand-alone code
             if (TALSH_ONLY) return
             call interface_mpi_comm_rank (global_communicator,my_MPI_rank)

             call sync_mo(cspinor_all,my_MPI_master)

             call interface_mpi_bcast (ao_basis%nshells,1,my_MPI_master,global_communicator)
             call interface_mpi_bcast (ao_basis%nao,1,my_MPI_master,global_communicator)
             call interface_mpi_bcast (ao_basis%basis_angular,1,my_MPI_master,global_communicator)
             nshells       = ao_basis%nshells
             nao           = ao_basis%nao

             if (my_MPI_rank /= my_MPI_master) nullify (ao_basis%gtos)
             if (my_MPI_rank /= my_MPI_master) allocate(ao_basis%gtos(nshells))
             if (my_MPI_rank /= my_MPI_master) nullify (ao_basis%shell_indices)
             if (my_MPI_rank /= my_MPI_master) allocate(ao_basis%shell_indices(nao))

!            Note that this communication step could be speeded up by packing and unpacking the data for all shells
!            Assuming that the time for this step is insignificant, a simpler option is chosen here.
             do ish = 1, nshells
                call interface_mpi_bcast (ao_basis%gtos(ish)%orb_momentum,1,my_MPI_master,global_communicator)
                call interface_mpi_bcast (ao_basis%gtos(ish)%atom_number,1,my_MPI_master,global_communicator)
                call interface_mpi_bcast (ao_basis%gtos(ish)%n_primitives,1,my_MPI_master,global_communicator)
             end do

!            We now know the size of the shells, the exponents and coefficients are already allocated on master,
!            but not on the slaves, so do this first before starting the next communication loop
             if (my_MPI_rank /= my_MPI_master) then
                do ish = 1, nshells
                   n_primitives = ao_basis%gtos(ish)%n_primitives
                   nullify  (ao_basis%gtos(ish)%exponent)
                   allocate (ao_basis%gtos(ish)%exponent(n_primitives))
                   nullify  (ao_basis%gtos(ish)%coefficient)
                   allocate (ao_basis%gtos(ish)%coefficient(n_primitives))
                end do
             end if

             do ish = 1, nshells
                n_primitives = ao_basis%gtos(ish)%n_primitives
                call interface_mpi_bcast (ao_basis%gtos(ish)%coord,3,my_MPI_master,global_communicator)
                call interface_mpi_bcast (ao_basis%gtos(ish)%exponent,n_primitives,my_MPI_master,global_communicator)
                call interface_mpi_bcast (ao_basis%gtos(ish)%coefficient,n_primitives,my_MPI_master,global_communicator)
             end do

!            We finally broadcast the pointers from a basis function index to its shell that were determined on master
             call interface_mpi_bcast (ao_basis%shell_indices,nao,my_MPI_master,global_communicator)

!            One electron integrals
             call interface_mpi_bcast_l0 (one_el_exist,1,my_MPI_master,global_communicator)
             if(one_el_exist) then
               call interface_mpi_bcast (int_1el%n_spinor,1,my_MPI_master,global_communicator)
               nsp=int_1el%n_spinor
               call interface_mpi_bcast (int_1el%e_core,1,my_MPI_master,global_communicator)

               if (my_MPI_rank /= my_MPI_master) nullify (int_1el%h_core)
               if (my_MPI_rank /= my_MPI_master) allocate(int_1el%h_core(nsp,nsp))

               call interface_mpi_bcast (int_1el%h_core,2*nsp*nsp,my_MPI_master,global_communicator)
             end if

!            Property integrals
             call interface_mpi_bcast (int_1el%n_prop,1,my_MPI_master,global_communicator)
             if (int_1el%n_prop>0) then
                if (my_MPI_rank /= my_MPI_master) then
                   nullify  (int_1el%h_prop)
                   allocate (int_1el%h_prop(nsp,nsp,int_1el%n_prop))
                   nullify  (int_1el%property_labels)
                   allocate (int_1el%property_labels(int_1el%n_prop))
                end if
                call interface_mpi_bcast (int_1el%h_prop,2*nsp*nsp*int_1el%n_prop,my_MPI_master,global_communicator)
                call interface_mpi_bcast(int_1el%property_labels,len(int_1el%property_labels(1))*int_1el%n_prop, &
                                            my_MPI_master,global_communicator)
             end if
#else
             return ! For running a serial version no action is required.
#endif

         end subroutine synchronize_global_data

         subroutine delete_global_data
         !TODO: check for memory leaks here. No so important as we usually stop right after and there is not much data anyhow
             deallocate(ao_basis%shell_indices)
             deallocate(ao_basis%gtos)
             call dealloc_mo (cspinor_all)
             if(one_el_exist) then
               deallocate(int_1el%h_core)
             end if
             if(int_1el%n_prop>0) then
               deallocate(int_1el%h_prop)
               deallocate(int_1el%property_labels)
               int_1el%n_prop=0
             end if
         end subroutine delete_global_data

         integer function nfunctions(l_value,basis_angular)
         ! Computes number of Cartesian/spherical functions for a given l-value
             integer, intent(in) :: l_value
             integer, intent(in) :: basis_angular
             if (basis_angular == 1) then
!               Cartesian
                nfunctions = l_value * ( l_value + 1 ) / 2
             else if (basis_angular == 2) then
!               spherical
                nfunctions = 2 * (l_value - 1) + 1
             end if
         end function nfunctions

         integer function get_shell_index(ifun)
             integer, intent(in) :: ifun
             get_shell_index = ao_basis%shell_indices(ifun)
         end function get_shell_index

         integer function get_shell_size(ish)
             integer, intent(in) :: ish
             integer :: il
             il=ao_basis%gtos(ish)%orb_momentum
             get_shell_size = nfunctions(il,ao_basis%basis_angular)
         end function get_shell_size

         subroutine get_shell_offset(l,l0)

             integer, intent(in)  :: l
             integer, intent(out) :: l0
             integer :: i, lshells

             lshells=get_shell_index(l)
             l0=0
             do i = 1,lshells-1
               l0=l0+get_shell_size(i)
             end do
         end subroutine get_shell_offset

         type(basis_func_info_t) function get_gto(ish)
             integer, intent(in) :: ish
             get_gto = ao_basis%gtos(ish)
         end  function get_gto

         subroutine get_gtos(first_ao,nao,gtos,nshells)

          integer, intent(in):: first_ao, nao
          type(basis_func_info_t), intent(inout), allocatable :: gtos(:)
          integer, intent(out) :: nshells
          integer :: ish, first_shell, last_shell

          first_shell = get_shell_index (first_ao)
          last_shell  = get_shell_index (first_ao+nao-1)
          nshells = last_shell - first_shell + 1

          allocate (gtos(nshells))
          do ish = 1, nshells
             gtos(ish) = get_gto(first_shell+ish-1)
          end do

         end subroutine get_gtos

     subroutine get_MO_coefficients (cspinor,mo_list,nmo)
!***********************************************************************
!
!     Fill array with subset of MO coefficients and related information
!     Retrives this from the qmo_all array that should be synced
!
!     Input
!        nmo     : number of selected MOs
!        mo_list : indices of selected MOs
!
!     Output
!        QMO  : MOs (coefficients, energies, etc.)
!
!
!     Written by L.Visscher Sep 2018
!
!
!***********************************************************************

      use exacorr_mo

      integer, intent(in)     :: nmo,mo_list(nmo)
      type(cmo), intent(inout) :: cspinor

      integer  :: i, nao, nmo_all

      nao     = cspinor_all%nao
      nmo_all = cspinor_all%nmo !WHY IS THIS HERE?

      ! allocate arrays in cmo to hold the active MO's specified in mo_list
      call alloc_mo (cspinor,nao,nmo)

      ! copy the requested mo's
      do i = 1, nmo
         cspinor%index(i)       = mo_list(i)
         cspinor%coeff(:,i,:)   = cspinor_all%coeff(:,mo_list(i),:)
         cspinor%energy(i)      = cspinor_all%energy(mo_list(i))
         cspinor%boson_irrep(i) = cspinor_all%boson_irrep(mo_list(i))
      enddo

      end subroutine get_MO_coefficients

      subroutine get_eps(eps,lmo)

           real(8)  :: eps(:)
           integer  :: lmo(:)

           eps=cspinor_all%energy(lmo)

      end subroutine get_eps


      logical function exist_one_el(e_core)

        real(8), intent(out) :: e_core

        exist_one_el=one_el_exist
        if (one_el_exist) then
           e_core=int_1el%e_core
        else
           e_core = 0.D0
        end if

      end function

      subroutine get_one_el(integrals,mo1_list,nmo1,mo2_list,nmo2,min_occ)

        complex(8), intent(out) :: integrals(:,:) !one lectron integrals
        integer, intent(in)     :: nmo1          ! the length of the mo basis
        integer, intent(in)     :: mo1_list(:)   ! and their indices
        integer, intent(in)     :: nmo2          ! the length of the mo basis
        integer, intent(in)     :: mo2_list(:)   ! and their indices
        integer, intent(in)     :: min_occ       ! smallest index, assumes same closed shell as MOLTRA

        integer :: pmo, qmo

        if((mo1_list(nmo1)-min_occ+1).gt.int_1el%n_spinor) call quit ('Not enough orbitals in MRCONEE (MOLTRA)')
        if((mo2_list(nmo2)-min_occ+1).gt.int_1el%n_spinor) call quit ('Not enough orbitals in MRCONEE (MOLTRA)')

        do qmo = 1, nmo2
          do pmo = 1, nmo1
            integrals(pmo, qmo) = int_1el%h_core(mo1_list(pmo)-min_occ+1,mo2_list(qmo)-min_occ+1)
          end do
        end do

      end subroutine get_one_el

      subroutine add_finitefield (h,exa_input,iff,mo1_list,nmo1,mo2_list,nmo2,min_occ)

        complex(8), intent(inout)     :: h(:,:)      ! one_body matrix to which the property field should be added
        type(exacc_input), intent(in) :: exa_input   ! field strength and label
        integer, intent(in)           :: iff
        integer, intent(in)           :: nmo1          ! the length of the mo basis
        integer, intent(in)           :: mo1_list(:)   ! and their indices
        integer, intent(in)           :: nmo2          ! the length of the mo basis
        integer, intent(in)           :: mo2_list(:)   ! and their indices
        integer, intent(in)           :: min_occ       ! smallest index, assumes same closed shell as MOLTRA

        integer :: pmo, qmo, i_prop, i

        if((mo1_list(nmo1)-min_occ+1).gt.int_1el%n_spinor) call quit ('Not enough orbitals in MRCONEE (MOLTRA)')
        if((mo2_list(nmo2)-min_occ+1).gt.int_1el%n_spinor) call quit ('Not enough orbitals in MRCONEE (MOLTRA)')

        do i=1,exa_input%nff(1)

          call get_ff_ind (i_prop,exa_input%ff_names(i))

          do qmo = 1, nmo2
            do pmo = 1, nmo1
              h(pmo, qmo) = h(pmo, qmo) + &
              int_1el%h_prop(mo1_list(pmo)-min_occ+1,mo2_list(qmo)-min_occ+1,i_prop) * exa_input%ff(i,iff)
            end do
          end do
        end do

      end subroutine add_finitefield

      subroutine get_ff_ind (i_prop,p_label)

        character(8), intent(in)      :: p_label   ! label
        integer,intent(out)           :: i_prop

        logical :: found

        found = .false.
        do i_prop = 1, int_1el%n_prop
          !print*,int_1el%property_labels(i_prop),'=?=',p_label
          if (int_1el%property_labels(i_prop) == p_label) then
            found = .true.
            exit
          end if
        end do

        if (.not.found) then
          print*, " ERROR: property ",p_label," not found on property integrals file"
          stop "property to be added not found"
          i_prop=-1
        end if

      end subroutine get_ff_ind

      subroutine get_ff_mat (h,i_prop,mo1_list,nmo1,mo2_list,nmo2,min_occ)

        complex(8), intent(inout)     :: h(:,:)      ! one_body matrix to which the property field should be added
        integer, intent(in)           :: i_prop
        integer, intent(in)           :: nmo1          ! the length of the mo basis
        integer, intent(in)           :: mo1_list(:)   ! and their indices
        integer, intent(in)           :: nmo2          ! the length of the mo basis
        integer, intent(in)           :: mo2_list(:)   ! and their indices
        integer, intent(in)           :: min_occ       ! smallest index, assumes same closed shell as MOLTRA

        integer :: pmo, qmo

        if((mo1_list(nmo1)-min_occ+1).gt.int_1el%n_spinor) call quit ('Not enough orbitals in MDPROP (MOLTRA)')
        if((mo2_list(nmo2)-min_occ+1).gt.int_1el%n_spinor) call quit ('Not enough orbitals in MDPROP (MOLTRA)')

        do qmo = 1, nmo2
          do pmo = 1, nmo1
            h(pmo, qmo) = h(pmo, qmo) + &
            int_1el%h_prop(mo1_list(pmo)-min_occ+1,mo2_list(qmo)-min_occ+1,i_prop)
          end do
        end do

      end subroutine get_ff_mat

      subroutine string2mo ( inp_str, a_or_b, nmo, mo_list)

        use exacorr_utils, only : quicksort

        implicit none

        character(len = *), intent(in)        :: inp_str
        !indicate if alpha or beta orbitals
        logical, intent(in)                   :: a_or_b
        integer,intent(out),optional          :: nmo
        integer,optional                      :: mo_list(:)

        integer :: str_len, sl, p_len, check, nmo_all
        integer :: mo1, mo2, i, j, k, i_mo, size_mo
        character(len=20) :: sav_str
        character(len=10) :: num_str
        character(len=1)  :: p_str, c_str, e_str
        real(8)                  :: e_min, e_max, e_dist
        real(8), allocatable     :: dist_a(:), eps(:)
        integer, allocatable     :: lmo(:)
        logical                  :: do_nmo, do_mo

        do_nmo=.false.
        do_mo=.false.
        if (present(nmo)) then
          do_nmo = .true.
          if (present(mo_list)) stop 'use only one output in string2mo'
          nmo=0
        else if (present(mo_list)) then
          do_mo = .true.
          mo_list=-9898
          size_mo=size(mo_list)
        else
          stop 'no output in string2mo'
        end if

        str_len=len(inp_str)
        num_str='0123456789'
        p_str='.'
        c_str=','
        e_str=' '
        sl=0
        p_len=0
        check=0
        nmo_all = cspinor_all%nmo/2

        !!!!!!!!!!!!!!!!!!!!!! use energy definition of spinors 
        if (inp_str(1:6) == 'energy') then
          !get orbital energies
          allocate(lmo(nmo_all))
          do i = 1, nmo_all
            if (a_or_b) then
              lmo(i)=i*2-1
            else
              lmo(i)=i*2
            end if
          end do
          allocate(eps(nmo_all))
          eps=cspinor_all%energy(lmo)
          deallocate(lmo)

          !read thresholds
          i=7
          do k=1,100
            if (inp_str(i:i).eq.' ') then
              i=i+1
            else
              EXIT
            end if
          end do
          j=index(inp_str(i:),' ')
          j=i+j
          read(inp_str(i:j-1) , *) e_min
          i=j
          do k=1,100
            if (inp_str(i:i).eq.' ') then
              i=i+1
            else
              EXIT
            end if
          end do
          j=index(inp_str(i:),' ')
          j=i+j
          read(inp_str(i:j-1) , *) e_max
          i=index(inp_str(j+1:),' ')+j
          read(inp_str(j-1:i) , *) e_dist
          if (e_min.ge.e_max) stop 'string2mo use: energy e_min, e_max, e_dist'

          sl=size(eps)
          allocate(dist_a(sl))

          dist_a=eps-e_min
          where(dist_a<0) dist_a=1000000
          j=minloc(abs(dist_a),1)
          mo1=1
          do i=j-1,1,-1
            if (abs(eps(i+1)-eps(i))>e_dist) then
              mo1=i+1
              EXIT
            end if
          end do

          dist_a=eps-e_max
          where(dist_a>0) dist_a=1000000
          j=minloc(abs(dist_a),1)
          mo2=sl
          do i=j,sl-1
            if (abs(eps(i)-eps(i+1))>e_dist) then
              mo2=i
              EXIT
            end if
          end do
          if (do_nmo) then
            nmo=mo2-mo1+1
          end if

          if (do_mo) then
            i=mo2-mo1+1
            if(i.ne.size_mo) stop 'mo_list wrong size in string2mo'
            i_mo=0
            do k=mo1,mo2
              i_mo=i_mo+1
              mo_list(i_mo)=k
            end do
          end if
          deallocate(dist_a)
          deallocate(eps)
        !!!!!!!!!!!!!!!!!!!!!! use orbitals strings for spinors
        else
          i_mo=0
          do i=1,str_len
            if (inp_str(i:i)==e_str) then
              cycle
            else if (inp_str(i:i)==p_str) then
              p_len=p_len+1
              if (p_len==1) then
                read(sav_str(1:sl) , *) mo1
                sl=0
              else if (p_len>2) then
                stop 'wrong number of points in string2mo'
              end if
            else if (inp_str(i:i)==c_str) then
              if(p_len==0) then
                read(sav_str(1:sl) , *) mo1
                if (mo1.le.nmo_all) then
                  if (do_nmo) then
                    nmo=nmo+1
                  end if
                  if (do_mo) then
                    i_mo=i_mo+1
                    mo_list(i_mo)=mo1
                  end if
                end if
                sl=0
              else if (p_len==2) then
                read(sav_str(1:sl) , *) mo2
                if (mo1>mo2) then
                  k=mo1
                  mo1=mo2
                  mo2=k
                end if
                if (mo1<nmo_all) then
                  if (mo2>nmo_all) then
                    mo2=nmo_all
                  end if
                  if (do_nmo) then
                    nmo=nmo+mo2-mo1+1
                  end if
                  if (do_mo) then
                    do k=mo1,mo2
                      i_mo=i_mo+1
                      mo_list(i_mo)=k
                    end do
                  end if
                end if
                p_len=0
                sl=0
              else
                stop 'not enough points in string2mo'
              end if
            else
              check=0
              do j=1,10
                if (inp_str(i:i).eq.num_str(j:j)) then
                  sl=sl+1
                  sav_str(sl:sl)=inp_str(i:i)
                  check=check+1
                end if
              end do
              if (check.ne.1) stop 'wrong character in string2mo'
            end if
          end do
          if (sl>0) then
            if(p_len==0) then
              read(sav_str(1:sl) , *) mo1
              if (mo1.le.nmo_all) then
                if (do_nmo) then
                  nmo=nmo+1
                end if
                if (do_mo) then
                  i_mo=i_mo+1
                  mo_list(i_mo)=mo1
                end if
              end if
              sl=0
            else if (p_len==2) then
              read(sav_str(1:sl) , *) mo2
              if (mo1>mo2) then
                k=mo1
                mo1=mo2
                mo2=k
              end if
              if (mo1<nmo_all) then
                if (mo2>nmo_all) then
                  mo2=nmo_all
                end if
                if (do_nmo) then
                  nmo=nmo+mo2-mo1+1
                end if
                if (do_mo) then
                  do k=mo1,mo2
                    i_mo=i_mo+1
                    mo_list(i_mo)=k
                  end do
                end if
              end if
              p_len=0
              sl=0
            else
              stop 'not enough points in string2mo'
            end if
          end if
          if (do_mo) then
            if (i_mo.ne.size_mo) stop 'mo_list wrong size in string2mo'
            !sort list
            j=1
            call quicksort(mo_list, j, size_mo) ! to do: replace by lapack
            !check for elements occuring twice
            do j=1,size_mo-1
              if (mo_list(j)==mo_list(j+1)) stop 'mo appearing twice in string2mo'
            end do
          end if
        end if
      end subroutine string2mo

      subroutine make_spinor_list ( exa_input, string_occupied, &
        string_occ_beta, string_virtual, string_vir_beta)

        use exacorr_datatypes

        implicit none

        type(exacc_input), intent(inout) :: exa_input
        character(len = *), intent(in)   :: string_occupied
        character(len = *), intent(in)   :: string_virtual
        character(len = *), intent(in)   :: string_occ_beta
        character(len = *), intent(in)   :: string_vir_beta

        integer, allocatable :: mobe_occ(:), mobe_vir(:)
        integer              :: i, j, i0, cont
        integer              :: nmo(4)
        integer              :: imo, ialpha, ibeta

        if (exa_input%beta_occ) then
          if (exa_input%beta_vir) then
            nmo(1:4)=0

            call string2mo(string_occupied,.true.,nmo=nmo(1))
            call string2mo(string_occ_beta,.false.,nmo=nmo(2))
            call string2mo(string_virtual,.true.,nmo=nmo(3))
            call string2mo(string_vir_beta,.false.,nmo=nmo(4))

            exa_input%nkr_occ=nmo(1)
            allocate (exa_input%mokr_occ(exa_input%nkr_occ))
            allocate (mobe_occ(nmo(2)))
            exa_input%nkr_vir=nmo(3)
            allocate (exa_input%mokr_vir(exa_input%nkr_vir))
            allocate (mobe_vir(nmo(4)))

            call string2mo(string_occupied,.true.,mo_list=exa_input%mokr_occ)
            call string2mo(string_occ_beta,.false.,mo_list=mobe_occ)
            call string2mo(string_virtual,.true.,mo_list=exa_input%mokr_vir)
            call string2mo(string_vir_beta,.false.,mo_list=mobe_vir)

            if(exa_input%mokr_occ(1).ne.mobe_occ(1)) then
              write(*,*) "core orbital need to be closed shell"
              call quit ("change .OCC_BETA in EXACORR")
            end if
            
            exa_input%nocc=nmo(1)+nmo(2)
            exa_input%nvir=nmo(3)+nmo(4)
            allocate (exa_input%mo_occ(exa_input%nocc))
            allocate (exa_input%mo_vir(exa_input%nvir))

            imo=0
            do i=1, max(nmo(1), nmo(2))
              if(i.le.nmo(1)) then
                imo=imo+1
                exa_input%mo_occ(imo)=exa_input%mokr_occ(i)*2-1
              end if
              if(i.le.nmo(2)) then
                imo=imo+1
                exa_input%mo_occ(imo)=mobe_occ(i)*2
              end if
            end do

            ialpha=1
            ibeta=1
            imo=0
            do i=1, nmo(3)+nmo(4)
              if(ibeta.gt.nmo(4)) then
                if (ialpha.le.nmo(3)) then
                  imo=imo+1
                  exa_input%mo_vir(imo)=exa_input%mokr_vir(ialpha)*2-1
                  ialpha=ialpha+1
                end if
              else if (ialpha.gt.nmo(3)) then
                if (ibeta.le.nmo(4)) then
                  imo=imo+1
                  exa_input%mo_vir(imo)=mobe_vir(ibeta)*2
                  ibeta=ibeta+1
                end if
           else if(exa_input%mokr_vir(ialpha).eq.mobe_vir(ibeta)) then
                imo=imo+1
                exa_input%mo_vir(imo)=exa_input%mokr_vir(ialpha)*2-1
                ialpha=ialpha+1
                imo=imo+1
                exa_input%mo_vir(imo)=mobe_vir(ibeta)*2
                ibeta=ibeta+1
           else if (exa_input%mokr_vir(ialpha).lt.mobe_vir(ibeta)) then
                imo=imo+1
                exa_input%mo_vir(imo)=exa_input%mokr_vir(ialpha)*2-1
                ialpha=ialpha+1
              else
                imo=imo+1
                exa_input%mo_vir(imo)=mobe_vir(ibeta)*2
                ibeta=ibeta+1
              end if
            end do

            deallocate (mobe_occ)
            deallocate (mobe_vir)
          else
            call quit ("also define .OCC_BETA in EXACORR")
          end if
        else if (exa_input%beta_vir) then
          call quit ("also define .VIR_BETA in EXACORR")
        else
          call string2mo(string_occupied,.true.,nmo=exa_input%nkr_occ)
          allocate (exa_input%mokr_occ(exa_input%nkr_occ))
          call string2mo(string_occupied,.true.,mo_list=exa_input%mokr_occ)
          exa_input%nocc = 2 * exa_input%nkr_occ
          allocate (exa_input%mo_occ(exa_input%nocc))
          do i = 1, exa_input%nkr_occ
            exa_input%mo_occ(2*i-1) = exa_input%mokr_occ(i)*2-1
            exa_input%mo_occ(2*i) = exa_input%mokr_occ(i)*2
          end do

          call string2mo(string_virtual,.true.,nmo=exa_input%nkr_vir)
          allocate (exa_input%mokr_vir(exa_input%nkr_vir))
          call string2mo(string_virtual,.true.,mo_list=exa_input%mokr_vir)
          !remove overlapping orbitals
          allocate (exa_input%mo_vir(exa_input%nkr_vir))
          i0=0
          do i = 1, exa_input%nkr_vir
            cont=0
            do j = 1, exa_input%nkr_occ
              if (exa_input%mokr_vir(i)==exa_input%mokr_occ(j)) cont=cont+1
            end do
            if (cont==0) then
              i0=i0+1
              exa_input%mo_vir(i0)=exa_input%mokr_vir(i)
            end if
          end do
          exa_input%nkr_vir=i0
          deallocate(exa_input%mokr_vir)
          !add spinors
          allocate(exa_input%mokr_vir(exa_input%nkr_vir))
          exa_input%mokr_vir(1:exa_input%nkr_vir)=exa_input%mo_vir(1:exa_input%nkr_vir)
          deallocate(exa_input%mo_vir)
          exa_input%nvir = 2 * exa_input%nkr_vir
          allocate (exa_input%mo_vir(exa_input%nvir))
          do i = 1, exa_input%nkr_vir
            exa_input%mo_vir(2*i-1) = 2*exa_input%mokr_vir(i)-1
            exa_input%mo_vir(2*i) = 2*exa_input%mokr_vir(i)
          end do
        end if
        
        if (exa_input%print_level>1) call print_orbitals(exa_input)

      end subroutine make_spinor_list

      subroutine print_orbitals ( exa_input )

        use exacorr_datatypes

        implicit none

        type(exacc_input), intent(inout) :: exa_input
        integer                          :: i

        write(*,*) ' --- read active occupied spinors --- '
        write(*,*) '   #             E   '
        do i=1,exa_input%nocc
          write(*,'(I8,E18.9)') exa_input%mo_occ(i), cspinor_all%energy(exa_input%mo_occ(i))
        end do 
        
        write(*,*) ' --- read active virtual spinors --- '
        write(*,*) '   #             E   '
        do i=1,exa_input%nvir
          write(*,'(I8,E18.9)') exa_input%mo_vir(i), cspinor_all%energy(exa_input%mo_vir(i))
        end do 

      end subroutine print_orbitals

      subroutine get_scf ( scf_energy )
        real(8), intent(out) :: scf_energy
        scf_energy=cspinor_all%total_energy
      end subroutine get_scf 

end module exacorr_global
