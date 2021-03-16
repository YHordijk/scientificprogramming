!
!    Collection of routines needed for exacorr that use the DIRAC common
!    block infrastructure.
!
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXA_MOLTRA()
C***********************************************************************
C
C     Driver routine for exacorr MOLTRA implementation
C
C     Written by L.Visscher Jan 2017
C
C***********************************************************************

      use exacc_cfg
      use exacorr_ao_to_mo,only     : exacorr_compute_mo_integrals
      use talsh_ao_to_mo,only       : talsh_compute_mo_integrals
      use exacorr_utils, only       : print_date
      use, intrinsic:: ISO_C_BINDING
#if defined (VAR_MPI)
      use interface_to_mpi
#endif

      implicit none

#include "../moltra/dcbtra.h"
! to access tra4_indstr
#include "infpar.h"
! we need to broadcast tra4_indstr
#include "dcborb.h"
! to access nesh (number of positive energy "electronlike" orbitals)
#include "dgroup.h"
#include "dcbham.h"
!information about the group

      character (len=72)   :: orbital_string(4) ! String with range of MO's to be transformed for each index
      integer              :: kvec              ! Dummy variable needed in first call to numlst
      integer, allocatable :: mo_list(:), mo_list_tmp(:)
      integer              :: nmo(4)
      integer              :: i, j, imo, jmo, ntmo, rcw

      exa_input%print_level=IPRTRA

#if defined (VAR_MPI)
      ! parallel version of DIRAC
C     Transfer necessary information from common blocks
      call interface_mpi_comm_rank(global_communicator,mytid)
      if (mytid == 0) then
         orbital_string = tra4_indstr(:,1) ! This is to be changed once we have an exacorr-specific input driver
         ! instruct slaves to go to the current routine (exa_moltra) as well
         call dirac_parctl( 19 )
         call interface_mpi_bcast_l0(exa_input%talsh,1,0,
     &                               global_communicator)
      else
         call interface_mpi_bcast_l0(exa_input%talsh,1,0,
     &                               global_communicator)
         if (exa_input%talsh) return ! TALSH has no use for MPI slaves
      end if
#else
      ! serial version of DIRAC
      if (exa_input%talsh) then
         ! TALSH version works fine with serial version of DIRAC
         orbital_string = tra4_indstr(:,1) ! This is to be changed once we have an exacorr-specific input driver
      else
         ! EXATENSOR version does not work with serial version of DIRAC
         print*, "Attempting to use exatensor with a non MPI version"
         call quit ("error in setup: exatensor requires MPI")
      end if
#endif

! next block is code is only for master
#if defined (VAR_MPI)
      if (mytid == 0) then
#endif
C
C     Loop over indices to convert string information to list of orbitals
C
C     Count number of orbitals to be transformed for each index 
      ntmo = 0
      do i = 1, 4
        nmo(i) = 0
        call numlst (orbital_string(i),kvec,1,1,nesh(1),1,nmo(i))
        ntmo = ntmo + nmo(i)
      end do

C     Fill the array with orbital indices
      allocate (mo_list_tmp(ntmo))
      imo = 1
      do i = 1, 4
       call numlst
     & (orbital_string(i),mo_list_tmp(imo),nmo(i),1,nesh(1),1,nmo(i))
       imo = imo + nmo(i)
      end do

C     Convert to Kramers unrestricted picture
      ntmo = ntmo * 2 ! spinors
      allocate (mo_list(ntmo))
      imo = 1
      jmo = 1
      do j = 1, 4
        do i = 1, 2*nmo(j)
            mo_list(imo) = (mo_list_tmp(jmo)-1)*2 + i
            imo = imo + 1
        end do
        jmo = jmo + nmo(j)
      end do
      deallocate(mo_list_tmp)
      nmo = nmo * 2 ! Kramers unrestricted

#if defined (VAR_MPI)
      end if
#endif
C
      IF (SPINFR.OR.NZ.EQ.1) THEN
         rcw = 1
      ELSE
         rcw = 2
      ENDIF

      if (.true.) then !Always use talsh for now
        call talsh_compute_mo_integrals(nmo,mo_list,
     &      exa_input%talsh_buff,exa_input%print_level,rcw,TH_CHOLESKY)
      else
#if defined (VAR_MPI)
         !slaves need the same info to get started
         call interface_mpi_bcast(ntmo,1,mparid,global_communicator)
         if (mytid /= 0) allocate(mo_list(ntmo))
         call interface_mpi_bcast(nmo,4,mparid,global_communicator)
         call interface_mpi_bcast(mo_list,ntmo,mparid,
     &                            global_communicator)
#endif
         call exacorr_compute_mo_integrals (nmo,mo_list,
     &                            exa_input%print_level)
      end if
C
      !Only master should deallocate, since only master completes talsh
      if (mytid == 0) deallocate (mo_list)

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXA_CC()
C***********************************************************************
C
C     Driver routine for exacorr CC implementation
C
C     Written by L.Visscher Jul 2017
C
C***********************************************************************

      use exacc_cfg
      use talsh_cc,only             : talsh_cc_driver
      use exacorr_cc,only           : exacorr_cc_driver
      use exacorr_global
      use exacorr_utils, only       : print_exacorr_logo
      use exacorr_utils, only       : print_date
#if defined (VAR_MPI)
      use interface_to_mpi
#endif

      implicit none

      integer              :: mytid                    ! MPI rank

#if defined (VAR_MPI)
      ! parallel version of DIRAC
      ! the slaves do not have their input set yet, so call them
      ! and let them know whether we use talsh or exatensor
      call interface_mpi_comm_rank(global_communicator,mytid)
      if (mytid == 0 ) then
         call print_exacorr_logo
         call dirac_parctl( 20 ) 
      end if
      call interface_mpi_bcast_l0(exa_input%talsh,1,0,
     &                               global_communicator)
      if (mytid /= 0 .and. exa_input%talsh) return
#else
      ! serial version of DIRAC
      if (.not. exa_input%talsh) then
         ! EXATENSOR version does not work with serial version of DIRAC
         print*, "Attempting to use exatensor with a non MPI version"
         call quit ("error in setup: exatensor requires MPI")
      end if
      call print_exacorr_logo
#endif
C
      call set_talsh_only(exa_input%talsh)
      call initialize_global_data
      if (mytid == 0) call print_date("Initialized global data")

! next block is code is only for master
#if defined (VAR_MPI)
      if (mytid == 0) then
#endif

       call make_spinor_list ( exa_input, string_occupied,
     & string_occ_beta, string_virtual, string_vir_beta)

#if defined (VAR_MPI)
      end if
#endif

      if (exa_input%talsh) then
         call talsh_cc_driver (exa_input)
      else
#if defined (VAR_MPI)
         call exacc_sync_cw(0,mytid)
#endif
         !if lambda is activated set prop flags
         if (exa_input%lambda) then 
           if (mytid == 0) call set_prop_flags_ccsd()
         end if

         call exacorr_cc_driver (exa_input)
      end if

      call exacc_deallocate_cw()

      end

       subroutine Initialize_Shells (shell,ierr)

!      get information from common blocks about the basis set
!      raise an error for contracted sets that can not yet be handled by Interest

!      NB: this routine is outdated, replaced by read_from_dirac. Just keeping it in case we still need it.

       use exacorr_datatypes

#include "implicit.h"
#include "maxorb.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"

       integer ishell, iPrimStart, iPrimEnd, k, ierr
       type(basis_func_info_t) :: shell(*)

       do ishell = 1, nlrgsh
          shell(ishell)%orb_momentum = nhkt(ishell)
          shell(ishell)%atom_number  = ncent(ishell)
          shell(ishell)%n_primitives = nuco(ishell)
          nullify  (shell(ishell)%exponent)
          allocate (shell(ishell)%exponent(nuco(ishell)))
          nullify  (shell(ishell)%coefficient)
          allocate (shell(ishell)%coefficient(nuco(ishell)))
          iPrimStart = jstrt(ishell) + 1
          iPrimEnd   = jstrt(ishell) + nuco(ishell)
          k          = numcf(ishell)
         shell(ishell)%exponent = priexp(iPrimStart:iPrimEnd)
         shell(ishell)%coefficient = priccf(iPrimStart:iPrimEnd,k)
         shell(ishell)%coord = cent(ishell,1:3,1)
!        In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
         call compress_unc_shell (shell(ishell),ierr)
       end do

       return 
       end

       subroutine compress_unc_shell (shell,ierr)
!If basis functions were grouped together in the input, they appear as quasi-contracted functions
!with a lot of zero coefficients. Remove these zeroes and reduce the length to just one function.
!Routine will produce an error, if the function is a true contraced function (more than one non-zero coefficient)

          use exacorr_datatypes

          integer, intent(out)     :: ierr
          type(basis_func_info_t), intent(inout) :: shell

          integer i, j, non_zero
          real(8), allocatable :: exponents(:), coefficients(:)
          real(8), parameter   :: treshold=1.D-12

          allocate(exponents(shell%n_primitives))
          allocate(coefficients(shell%n_primitives))
          exponents    = shell%exponent
          coefficients = shell%coefficient

          ierr = 0
          non_zero = 0
          do i = 1, shell%n_primitives
             if (abs(coefficients(i)) > treshold) then
                non_zero = non_zero + 1
                j = i
             end if
          end do

          if (non_zero > 1) then
             ierr = 1
             return   ! leave shell untouched and return with error condition
          else if (non_zero == 0) then
             print*, " WARNING: compress_unc_shell found empty shell"
             ierr = 2
             j = 1    ! set j to the first function, to avoid complications with zero length arrays
          end if

          deallocate (shell%exponent)
          deallocate (shell%coefficient)
          allocate   (shell%exponent(1))
          allocate   (shell%coefficient(1))

          shell%n_primitives = 1
          shell%exponent     = exponents(j)
          shell%coefficient  = coefficients(j)

          deallocate (exponents)
          deallocate (coefficients)

       end subroutine compress_unc_shell

       subroutine read_from_dirac (ao_basis,cspinor,file_name,ierr)

       use exacorr_mo
       use exacorr_datatypes
       use exacorr_global,only: nfunctions
       use exacorr_utils

       implicit none

!      read all positive energy mo coefficients present on DFCOEF file
!      check for imcompatibiklities (use of symmetry) as well and continue
!      by reading also the basis set information into gtos

       character(len=*),intent(in) :: file_name
       type(basis_set_info_t), intent(inout) :: ao_basis
       type(cmo), intent(inout):: cspinor
       integer, intent(out)   :: ierr 
       integer, parameter     :: NOERROR  = 0 ! no error
       integer, parameter     :: NOTFOUND = 1 ! file not found
       integer, parameter     :: CORRUPT  = 2 ! file unreadable
       integer, parameter     :: SYMMETRY = 3 ! calculation used symmetry
       integer, parameter     :: NOBASIS  = 4 ! no basis set information could be read (old DFCOEF?)
       integer, parameter     :: CONTRAC  = 5 ! basis set is contracted (Interest can not handle this)
       integer, parameter     :: INCONSIS = 6 ! DFCOEF does not contain consistent information wrt number of basis functions
       character              :: text*74      ! holds text read from DFCOEF, only used for printing this
       integer                :: nsym
       integer                :: idim(3,2)
       real(8)                :: total_energy
       logical                :: tobe
       integer                :: lucoef
       integer                :: i, j, k, ishell, iprimstart, iprimend
       integer                :: nao, nao_bas, nshells, iao, ibas
       integer                :: nmo, npo, nto, npriexp, ncnt, nzbuf
       integer, parameter     :: nz = 4 ! we use no symmetry, so always have this
       real(8), allocatable   :: coeff(:,:,:),eig(:), seig(:)
       integer, allocatable   :: ibeig(:),index(:)
       integer, allocatable   :: ncent(:),nhkt(:),nuco(:),numcf(:)
       integer, allocatable   :: jstrt(:)
       real(8), allocatable   :: cent(:,:),priexp(:),priccf(:,:)
       complex(8), allocatable:: scoeff(:,:,:)
       integer, allocatable   :: ao_map(:,:)
       character(8)           :: label(4)
       logical                :: labeled

       ierr = NOERROR

       ! open file and read dimension info
       inquire (file=file_name,exist=tobe)
       if (.not.tobe) then
          print*, "error: MO coefficient file"//file_name//" not found"
          ierr = NOTFOUND
          return
       else
          call get_free_fileunit(lucoef)
          open (lucoef,file=file_name,status='OLD',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
          read (lucoef) label
          if (label(1) == '********') then
             labeled = .true. ! new format
             read (lucoef,end=10,err=10) text,nsym,nzbuf,
     &         ((idim(i,j),i=1,3),j=1,nsym),total_energy
          else
             labeled = .false.
             backspace(lucoef) ! old format (DIRAC19 and earlier)
             read (lucoef,end=10,err=10) text,nsym,
     &         ((idim(i,j),i=1,3),j=1,nsym),total_energy
          endif
       end if

       ! check symmetry and get dimensions of the coefficient array
       if (nsym > 1) then
          ierr = SYMMETRY
          close (lucoef,status='KEEP')
          return
       else
          nao = idim(3,1) ! number of basis functions
          npo = idim(1,1) ! "positronic" (= negative energy) solutions (can be thrown away for our purposes)
          nmo = idim(2,1) ! "electronic" (= positive energy) solutions (we only need these)
          nto = npo+nmo   ! total number of solutions on this file
          if (npo > 0) print*, "warning: reading from 4C DFCOEF file" ! this warning will go once we have the full index transform implemented
       end if

       allocate(coeff(nao,nto,nz))
       allocate(eig(nto))
       allocate(ibeig(nto))
       allocate(index(2*nmo))
        
       ! read the mo coefficients and orbital eigenvalues
       if (labeled) read (lucoef) label
       read (lucoef,end=10,err=10) coeff
       if (labeled) read (lucoef) label
       read (lucoef,end=10,err=10) eig
       if (labeled) read (lucoef) label
       read (lucoef,end=10,err=10) ibeig
       ! fill in index values
       do i = 1, 2*nmo
          index(i) = i ! we set the index compatible with the fact that we store all mo's
       end do

!      switch to unrestricted spinor picture
       allocate(seig(2*nmo))
       allocate(scoeff(nao,2*nmo,2)) ! nao, 2*nmo, spin
       call convert_quaternion_to_complex(coeff,scoeff,eig,seig,
     &                                              nao,nmo,npo)

       ! put in cmo
       call alloc_mo (cspinor,nao,2*nmo)
       cspinor%total_energy = total_energy
       cspinor%index        = index
       cspinor%coeff        = scoeff
       cspinor%energy       = seig
       cspinor%boson_irrep  = 0 ! because nosym. if there is, include ibeig in converter

       deallocate(coeff)
       deallocate(eig)
       deallocate(ibeig)
       deallocate(index)
       deallocate(seig)
       deallocate(scoeff)

       ! read the basis set information and check for consistency
       if (labeled) read (lucoef) label
       read (LUCOEF,err=20,end=20) nao_bas,nshells,npriexp,ncnt 
       if (nao < nao_bas) then
          ierr = INCONSIS
          close (lucoef,status='KEEP')
          return
       elseif (nao > nao_bas) then
          print*,"WARNING: assuming Levy-Leblond Hamiltonian"
          print*,"Reading fewer basis functions than MO-coefficients"
          allocate (ao_map(2,nao_bas))
          forall (i=1:nao_bas) ao_map(:,i)=i
          call compress_mo (cspinor,ao_map)
       end if

       ! allocate the arrays to read the basis set information in DIRAC format
       allocate(ncent(nshells))
       allocate(nhkt(nshells))
       allocate(nuco(nshells))
       allocate(numcf(nshells))
       allocate(jstrt(nshells))
       allocate(cent(nshells,3))
       allocate(priexp(npriexp))
       allocate(priccf(npriexp,ncnt))

!      orb_momentum, atom_number, n_primives for all shells
       read (LUCOEF,err=20,end=20) nhkt,ncent,nuco
!      start of exponent values, number of contracted for all shells
       read (LUCOEF,err=20,end=20) jstrt,numcf
!      coordinates (x, y, z) of center for all shells
       read (LUCOEF,err=20,end=20) cent
!      values of the exponents for all shells (packed into a 1-d array)
       read (LUCOEF,err=20,end=20) priexp
!      values of the contraction coefficients for all shells (packed into a 1-d array)
       read (LUCOEF,err=20,end=20) priccf

       ! put in gtos
       ao_basis%nshells       = nshells
       ao_basis%nao           = nao_bas
       ao_basis%basis_angular = 1
       nullify(ao_basis%gtos)
       allocate(ao_basis%gtos(nshells))
       do ishell = 1, nshells
          ao_basis%gtos(ishell)%orb_momentum = nhkt(ishell)
          ao_basis%gtos(ishell)%atom_number  = ncent(ishell)
          ao_basis%gtos(ishell)%n_primitives = nuco(ishell)
          nullify  (ao_basis%gtos(ishell)%exponent) 
          allocate (ao_basis%gtos(ishell)%exponent(nuco(ishell)))
          nullify  (ao_basis%gtos(ishell)%coefficient)
          allocate (ao_basis%gtos(ishell)%coefficient(nuco(ishell)))
          iPrimStart = jstrt(ishell) + 1
          iPrimEnd   = jstrt(ishell) + nuco(ishell)
          k          = numcf(ishell)
          ao_basis%gtos(ishell)%exponent = priexp(iPrimStart:iPrimEnd)
          ao_basis%gtos(ishell)%coefficient = 
     &                                   priccf(iPrimStart:iPrimEnd,k)
          ao_basis%gtos(ishell)%coord = cent(ishell,1:3)
!         In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
          call compress_unc_shell (ao_basis%gtos(ishell),ierr)
          if (ierr == 1) ierr = CONTRAC
       end do

       deallocate(ncent)
       deallocate(nhkt)
       deallocate(nuco)
       deallocate(numcf)
       deallocate(jstrt)
       deallocate(cent)
       deallocate(priexp)
       deallocate(priccf)

       ! make array with pointers from a basis function index to its shell
       nullify(ao_basis%shell_indices)
       allocate(ao_basis%shell_indices(nao))
       iao = 0
       do ishell = 1, nshells
          do ibas = 1,
     &       nfunctions(
     &       ao_basis%gtos(ishell)%orb_momentum,ao_basis%basis_angular)
             iao = iao + 1
             if (iao.gt.nao) then
                ierr = INCONSIS
                call quit ('Error in read_from_dirac')
             endif
             ao_basis%shell_indices(iao) = ishell
          end do
       end do
       if (iao.gt.nao) then
          ierr = INCONSIS
          close (lucoef,status='KEEP')
          call quit ('Error in read_from_dirac')
       endif

       close (lucoef,status='KEEP')

       call print_date("read from DFCOEF:"//text)
       return

  10   ierr = CORRUPT
       close (lucoef,status='KEEP')
       return

  20   ierr = NOBASIS
       close (lucoef,status='KEEP')
       return

       end subroutine read_from_dirac

       subroutine write_mos_to_dirac (ao_basis,cspinor,file_name)

       use exacorr_mo
       use exacorr_datatypes
       use exacorr_utils

       implicit none

!      write mo coefficients in dirac format

       character(len=*),intent(in) :: file_name
       type(basis_set_info_t), intent(inout) :: ao_basis
       type(cmo), intent(inout):: cspinor
       type(qumo)             :: qorbital
       character              :: text*74      ! holds text read from DFCOEF, only used for printing this
       integer                :: nsym
       real(8)                :: total_energy
       logical                :: tobe
       integer                :: lucoef
       integer                :: ishell
       integer                :: nao, nshells, iao, ibas
       integer                :: nmo, npo, nto, npriexp, ncnt
       integer, allocatable   :: numcf(:),jstrt(:)
       real(8), allocatable   :: priexp(:),priccf(:)

!      switch to Kramers-restricted format
       text = "MO coefficients produced by geniao "
       npo = 0  ! number of "positronic orbitals is hardwired to 0 for now
       nsym = 1 ! we do not use symmetry, so nsym is always 1
       call convert_complex_to_quaternion(cspinor,qorbital)
       ! print warning in case we overwrite a file and open the file
       call get_free_fileunit(lucoef)
       inquire (file=file_name,exist=tobe)
       if (tobe) then
          print*, "warning: overwriting file "//file_name
          open (lucoef,file=file_name,status='OLD',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
       else
          open (lucoef,file=file_name,status='NEW',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
       end if

       ! convert basis set info to 1-d lists used by DIRAC (assuming uncontracted basis sets)
       ! TODO: change storage on DIRAC side to use the basis set type
       nshells = ao_basis%nshells
       ! we always have uncontracted basis sets, so can hardwire the next variable and equate the number of exponents to the number of shells
       ncnt = 1
       npriexp =  nshells
       allocate(jstrt(nshells))
       allocate(numcf(nshells))
       allocate(priexp(nshells))
       allocate(priccf(nshells))
       do ishell = 1, nshells
          jstrt(ishell) = ishell-1
          numcf(ishell) = 1
          priexp(ishell) = ao_basis%gtos(ishell)%exponent(1)
          priccf(ishell) = ao_basis%gtos(ishell)%coefficient(1)
       end do
 
       ! write the header
       write (lucoef) text,nsym,npo,qorbital%nmo,
     &                qorbital%nao,qorbital%total_energy
       ! write the mo-coefficients
       write (lucoef) qorbital%coeff
       ! write the orbital energies
       write (lucoef) qorbital%energy
       ! write the boson symnmetry information
       write (lucoef) qorbital%boson_irrep
       ! write the ao basis information header
       write (lucoef) ao_basis%nao,ao_basis%nshells,npriexp,ncnt 
!      orb_momentum, atom_number, n_primives for all shells
       write (lucoef) ao_basis%gtos(:)%orb_momentum,
     &                ao_basis%gtos(:)%atom_number,
     &                ao_basis%gtos(:)%n_primitives
!      start of exponent values, number of contracted for all shells
       write (lucoef) jstrt,numcf
!      coordinates (x, y, z) of center for all shells
       write (lucoef) (ao_basis%gtos(ishell)%coord(:),
     &                 ishell=1,nshells)
!      values of the exponents for all shells (packed into a 1-d array)
       write (lucoef) priexp
!      values of the contraction coefficients for all shells (packed into a 1-d array)
       write (lucoef) priccf

       deallocate(jstrt)
       deallocate(numcf)
       deallocate(priexp)
       deallocate(priccf)

       close (lucoef,status='KEEP')

       return

       end subroutine write_mos_to_dirac

       subroutine write_mos_to_dirac_formatted(cspinor,file_name)

       use exacorr_mo
       use exacorr_utils

       implicit none

!      write mo coefficients in dirac format

       character(len=*),intent(in) :: file_name
       type(cmo), intent(inout):: cspinor
       type(qumo)             :: qorbital
       character              :: text*74      ! holds text read from DFCOEF, only used for printing this
       integer                :: nsym
       real(8)                :: total_energy
       logical                :: tobe
       integer                :: lucoef
       integer                :: nao
       integer                :: nmo, npo

!      switch to Kramers-restricted format
       text = "MO coefficients produced by geniao "
       npo = 0  ! number of "positronic orbitals is hardwired to 0 for now
       nsym = 1 ! we do not use symmetry, so nsym is always 1
       call convert_complex_to_quaternion(cspinor,qorbital)
       ! print warning in case we overwrite a file and open the file
       call get_free_fileunit(lucoef)
       inquire (file=file_name,exist=tobe)
       if (tobe) then
          print*, "warning: overwriting file "//file_name
          open (lucoef,file=file_name,status='OLD',FORM='FORMATTED',
     &     access='SEQUENTIAL')
       else
          open (lucoef,file=file_name,status='NEW',FORM='FORMATTED',
     &     access='SEQUENTIAL')
       end if

       ! write the header
       write (lucoef,'(A74)') text
       write (lucoef,'(7(X,I0))') nsym,npo,qorbital%nmo,qorbital%nao
       write (lucoef,'(E24.16)') qorbital%total_energy
       ! write the mo-coefficients
       write (lucoef,'(6F22.16)') qorbital%coeff
       ! write the orbital energies
       write (lucoef,'(6E22.12)') qorbital%energy
       ! write the boson symnmetry information
       write (lucoef,'(66(X,I0))') qorbital%boson_irrep

       close (lucoef,status='KEEP')
       return
       
       end subroutine write_mos_to_dirac_formatted

       subroutine Write_To_Mdcint (nmo,mo_tens,rcw)

!      Utility routine to write to the DIRAC moltra format, makes it easier to test the index transform
 
!      Note that the integral matrix is transposed in moltra, this is why the 3 and 4th indices are written to the left below
!      This routine should therefore be rewritten if the matrix transpose I(ij,kl) .ne. I(kl,ij)

#include "implicit.h"

       integer rcw    ! this is used to toggle between complex (so, rcw=2) and real (scalar, rcw=1) tensors
       integer iclass, nclass
       integer nmo(4)
       real(8) mo_tens(2,nmo(3),nmo(4),nmo(1),nmo(2))
       real(8), allocatable :: fmat(:,:,:),rkl(:,:)
       integer, allocatable :: indk(:),indl(:)
       integer nkr(4),bar(8,4)

       iprint = 0
       irepij = 1
       nkr = nmo / 2 ! go to Kramers restricted
       nclass = 8    ! write the standard classes used in MOLTRA as well

!      Open the MDCINT file
       call mdinti (0) 

!      We only need to write a subset of integrals, the rest will be generated by permutation and/or Kramers' symmetry
 
!      Classes are defined to have a standard subset of integrals, make these consistent with the definition in mdintb (in moltra)
!      Index   :  3    4      1    2
!      Class 1 : (k    l    | i    j   )
!      Class 2 : (kbar lbar | i    j   )
!      Class 3 : (k    lbar | i    jbar)
!      Class 4 : (kbar l    | i    jbar)
!      Class 9 : (kbar l    | i    j   )
!      Class 10: (k    lbar | i    j   )
!      Class 11: (k    l    | ibar j   )
!      Class 12: (k    l    | i    jbar)

!      NB: classes 5 through 8 are not used, we therefore shift the final 4 classes in the array below to get a contiguous index

       bar = 1 ! set everything unbarred to initialize
!      define bars for class 2
       bar(2,3) = 2
       bar(2,4) = 2
!      define bars for class 3
       bar(3,2) = 2
       bar(3,4) = 2
!      define bars for class 4
       bar(4,2) = 2
       bar(4,3) = 2
!      define bars for class 5=9
       bar(5,3) = 2
!      define bars for class 6=10
       bar(6,4) = 2
!      define bars for class 7=11
       bar(7,1) = 2
!      define bars for class 8=12
       bar(8,2) = 2

!      Allocate memory needed by mdintb
       nkl = nkr(3) * nkr(4)
       allocate(rkl(rcw,nkl)) 
       allocate(indk(nkl))
       allocate(indl(nkl)) 
       allocate(fmat(nkr(3),nkr(4),rcw))

!      loop over classes to write out subsets of integrals
       do iclass = 1, nclass
          if (iclass.le.4) then
             jclass = iclass
          else
             jclass = iclass + 4  ! shift iclass to obey the convention used in mdintb (given above)
          end if

!         loop over blocks of integrals, taking i and j as outermost indices like is done in moltra
          do jj = 1, nkr(2)
             j = (jj-1)*2 + bar(iclass,2) ! Spinor index 2
             do ii = 1, nkr(1)
                i = (ii-1)*2 + bar(iclass,1) ! Spinor index 1

!              In moltra the real/imaginary index within a block is last rather than first, so we need to reorganize

!              If the particle/particle (index 12 <-> index 34) symmetry is not fulfilled we should reorganize further
!              this can be done here, but will be quite inefficient and not worth the trouble for this utility routine that
!              is to discarded later anyhow.
     
               do ll = 1, nkr(4)
                  l = (ll-1)*2 + bar(iclass,4) ! Spinor index 4
                  do kk = 1, nkr(3)
                     k = (kk-1)*2 + bar(iclass,3) ! Spinor index 3
                     fmat(kk,ll,1) = mo_tens(1,k,l,i,j)
                     if (rcw.gt.1) fmat(kk,ll,2) = mo_tens(2,k,l,i,j)
                  end do
               end do

!              write a block of integrals to file. mdintb takes care of thresholding and labelling
               call mdintb(iprint,rcw,jclass,irepij,nkr(3),nkr(4),
     &           ii,jj,fmat,nkl,indk,indl,rkl)

             end do
          end do
       end do

!      Close the MDCINT file
       call mdintf 
       deallocate(rkl)
       deallocate(indk)
       deallocate(indl)
       deallocate(fmat)

       return
       end

       subroutine read_from_mrconee (int_1el, ierr)

       use exacorr_datatypes
       use exacorr_utils, only: get_free_fileunit

       implicit none

!      routine to read the one electron integrals from MRCONEE

       type(one_el_t), intent(out) :: int_1el
       integer, intent(out)   :: ierr
       integer, parameter     :: NOERROR  = 0 ! no error
       integer, parameter     :: NOTFOUND = 1 ! file not found
       integer, parameter     :: CORRUPT  = 2 ! file unreadable
       integer, parameter     :: INCONSIS = 6 ! MRCONEE does not match information read from DFCOEF

       integer :: nsp, nkr
       logical :: breit
       real(8) :: e_core
       complex(8), allocatable :: h_core(:,:)

       logical                :: tobe
       integer                :: luone

       ierr=NOERROR

!      open file and read setup infromation
       inquire (file='MRCONEE',exist=tobe)
       if (.not.tobe) then
          print*, "error: file MRCONEE not found"
          ierr = NOTFOUND
          return
       else
          call get_free_fileunit(luone)
          open (luone,file='MRCONEE',status='OLD',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
          read (luone,ERR=10,END=10) nsp,breit,e_core
       end if
!      save number of spinors
       int_1el%n_spinor=nsp
!      save core energy
       int_1el%e_core=e_core
!TODO  check consistency of e_core and nsp with the values read from DFCOEF !

!      skip symmetry information
       read (luone,ERR=10,END=10)
       read (luone,ERR=10,END=10)
       read (luone,ERR=10,END=10)
       read (luone,ERR=10,END=10)

!      allocate array to store the effective 1-body Hamiltonian
       nullify(int_1el%h_core)
       allocate(int_1el%h_core(NSP, NSP))

!      as the ordering on mrconee is different we need an auxilliary array to read in
       allocate(h_core(NSP, NSP))

!      read this "core" Hamiltonian
       read (luone,ERR=20,END=20) h_core

!      reorder and store in the data type
       nkr = nsp / 2 ! number of Kramers pairs is always half the number of spinors
       int_1el%h_core(1:nsp:2,1:nsp:2) = h_core(1:nkr,1:nkr)     ! unbar-unbar
       int_1el%h_core(2:nsp:2,1:nsp:2) = h_core(nkr+1:nsp,1:nkr)   ! bar-unbar
       int_1el%h_core(1:nsp:2,2:nsp:2) = h_core(1:nkr,nkr+1:nsp)   ! unbar-bar
       int_1el%h_core(2:nsp:2,2:nsp:2) = h_core(nkr+1:nsp,nkr+1:nsp) ! bar-bar

       deallocate (h_core)
       close (luone,status='KEEP')
       
       return

  10   ierr = CORRUPT
       return

  20   ierr = INCONSIS
       return

       end subroutine read_from_mrconee

       subroutine read_mdprop (int_1el, ierr)

       use exacorr_datatypes
       use exacorr_utils, only: get_free_fileunit

       implicit none

!      routine to read property integrals and add these to the int_1el data

       type(one_el_t), intent(inout) :: int_1el
       integer, intent(out)   :: ierr
       integer, parameter     :: NOERROR  = 0 ! no error
       integer, parameter     :: NOTFOUND = 1 ! file not found
       integer, parameter     :: CORRUPT  = 2 ! file unreadable

       integer, parameter     :: max_prop = 100 ! limit maximum of arrays to be read, should be high enough for the foreseeable future
       integer :: nsp, nkr, n_prop, i_prop
       complex(8), allocatable :: h_prop(:,:)

       character*32           :: label
       logical                :: tobe
       integer                :: luone

       ierr=NOERROR

!      open file and read setup infromation
       inquire (file='MDPROP',exist=tobe)
       if (.not.tobe) then
          print*, "error: file MDPROP not found"
          ierr = NOTFOUND
          return
       else
          call get_free_fileunit(luone)
          open (luone,file='MDPROP',status='OLD',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
       end if
!      set number of spinors
       nsp=int_1el%n_spinor

!      allocate auxilliary array to read in
       allocate(h_prop(nsp, nsp))

!      count the number of matrices that we need to store
!      MDROP has a very simple format: record with label followed by record containing the data
       n_prop = 0
       do while (n_prop < max_prop)
          read (luone,err=10,end=10) label
          if (label(25:32) == 'EOFLABEL') exit
          read (luone,err=10,end=10) h_prop
          n_prop = n_prop + 1
       end do
       rewind (luone)

!      store property integrals
       nullify(int_1el%property_labels)
       nullify(int_1el%h_prop)
       allocate(int_1el%property_labels(n_prop))
       allocate(int_1el%h_prop(nsp,nsp,n_prop))
       int_1el%n_prop = n_prop
       nkr = nsp / 2 ! number of Kramers pairs is always half the number of spinors
       i_prop = 0
       do while (i_prop < n_prop)
          read (luone,err=10,end=10) label
          i_prop = i_prop + 1
          read (luone,err=10,end=10) h_prop
        ! reorder and store in the data type
          int_1el%property_labels(i_prop) = label(25:32)
          int_1el%h_prop(1:nsp:2,1:nsp:2,i_prop)
     &          = h_prop(1:nkr,1:nkr)     ! unbar-unbar
          int_1el%h_prop(2:nsp:2,1:nsp:2,i_prop)
     &         = h_prop(nkr+1:nsp,1:nkr)   ! bar-unbar
          int_1el%h_prop(1:nsp:2,2:nsp:2,i_prop)
     &         = h_prop(1:nkr,nkr+1:nsp)   ! unbar-bar
          int_1el%h_prop(2:nsp:2,2:nsp:2,i_prop)
     &         = h_prop(nkr+1:nsp,nkr+1:nsp) ! bar-bar
       end do

       close (luone,status='KEEP')
       deallocate (h_prop)
       return

  10   ierr = CORRUPT
       return

       end subroutine read_mdprop

