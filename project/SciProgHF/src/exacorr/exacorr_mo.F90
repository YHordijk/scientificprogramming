module exacorr_mo 

!   Contains all the necessary information for a range of molecular orbitals
!   In particular:  - Orbital energy
!                   - Boson characterization (for spinfree runs)
!                   - Coefficients of AO-functions
!                   - Position in the full list of MO vectors

!   Extended from dirac_mo by Stan Papadopoulos, January 2019

    implicit none

    private

    public :: alloc_mo, dealloc_mo, sync_mo, compress_mo
    public :: convert_quaternion_to_complex,convert_complex_to_quaternion,apply_Kramers_operator

    type  :: mo_root
        integer             :: &
         nao, &               ! number of aos
         nmo                  ! number of mos

        logical             :: &
         allocated=.false.,    &        ! coefficient array is allocated
         restricted=.true.,    &        ! spinors (or spin-orbitals) are paired
         complex_algebra=.true.         ! coefficients are complex

        real(8)             :: &
         total_energy         ! total (SCF) energy as stored on the coefficient file

        real(8), pointer    :: &
         energy (:)           ! orbital energies

        real(8), pointer :: &
         coeff_r (:,:,:)      ! coefficients (AO,MO,spin)

        complex(8), pointer :: &
         coeff_c (:,:,:)      ! coefficients (AO,MO,spin)

        real(8), pointer :: &
         coeff_q (:,:,:)      ! coefficients (AO,MO,IZ)

    end type mo_root

    type, public, extends(mo_root) :: cmo

        integer, pointer    :: &
         index(:), &          ! index of the spinor in the full list of spinors for this fermion symmetry
         boson_irrep(:)       ! for spinfree calculations: boson irrep to which the spinor belongs

        integer             :: &
         fermion_irrep        ! fermion irrep (1=gerade, 2=ungerade) to which the spinors belongs

        complex(8), pointer :: &
         coeff (:,:,:)        ! coefficients (AO,MO,spin)

    end type cmo

    type, public, extends(mo_root) :: qumo

        integer, pointer    :: &
         index(:), &          ! index of the spinor in the full list of spinors for this fermion symmetry
         boson_irrep(:)       ! for spinfree calculations: boson irrep to which the spinor belongs

        integer             :: &
         fermion_irrep        ! fermion irrep (1=gerade, 2=ungerade) to which the spinors belongs

        integer :: &
         nz                   ! number of active quaternion units

        real(8), pointer :: &
         coeff   (:,:,:)      ! coefficients (AO,MO,iz)

    end type qumo

    interface alloc_mo
      module procedure alloc_mo_complex
      module procedure alloc_mo_quaternion
    end interface

    interface dealloc_mo
      module procedure dealloc_mo_complex
      module procedure dealloc_mo_quaternion
    end interface

    interface sync_mo
      module procedure sync_mo_complex
      module procedure sync_mo_quaternion
    end interface

    interface compress_mo
      module procedure compress_mo_complex
    end interface

    contains

     subroutine alloc_mo_complex(M,nao,nmo)

!     Prepare to work with a block of Mos

      type(cmo), intent(inout) :: M
      integer,   intent(in )   :: nao
      integer,   intent(in )   :: nmo

      M%nao = nao  ! Size of AO basis. Is kept here as we may want to use non-global defs
      M%nmo = nmo  ! Size of MO basis. May vary (e.g. in the 4-index when we have different ranges)

!     Reserve the memory needed to work with a block of Mos
      nullify  (M%energy)
      allocate (M%energy(M%nmo))
      M%energy = 0.0

      nullify  (M%coeff)
      allocate (M%coeff(M%nao,M%nmo,2))
      M%coeff = 0.0

      nullify  (M%boson_irrep)
      allocate (M%boson_irrep(M%nmo))
      M%boson_irrep = 0

      nullify  (M%index)
      allocate (M%index(M%nmo))
      M%index = 0
      M%allocated = .true.

     end subroutine alloc_mo_complex

     subroutine alloc_mo_quaternion(M,nao,nmo,nz)

!      Prepare to work with a block of Mos

       type(qumo), intent(inout)           :: M
       integer,    intent(in )             :: nao
       integer,    intent(in )             :: nmo
       integer,    intent(in ), optional   :: nz

       M%nao = nao  ! Size of AO basis. Is kept here as we may want to use non-global defs
       M%nmo = nmo  ! Size of MO basis. May vary (e.g. in the 4-index when we have different ranges)
       if (present(nz)) then
          M%nz = nz
       else
          M%nz = 1  ! We default to one, easy for spinfree applications
       endif

!      Reserve the memory needed to work with a block of Mos
       nullify  (M%energy)
       allocate (M%energy(M%nmo))
       M%energy = 0.0

       nullify  (M%coeff)
       allocate (M%coeff(M%nao,M%nmo,M%nz))
       M%coeff = 0.0

       nullify  (M%boson_irrep)
       allocate (M%boson_irrep(M%nmo))
       M%boson_irrep = 0

       nullify  (M%index)
       allocate (M%index(M%nmo))
       M%index = 0
       M%allocated = .true.

     end subroutine alloc_mo_quaternion

     subroutine dealloc_mo_complex(M)

!     Clean up workspace

      type(cmo), intent(inout) :: M

!     Set the dimensions to zero (not really necessary, but may prevent errors due to abuse of the MO-type)
      M%nao = 0
      M%nmo = 0

!     Free the memory associated with a block of Mos
      deallocate (M%energy)
      deallocate (M%coeff)
      deallocate (M%boson_irrep)
      deallocate (M%index)
      M%allocated = .false.

     end subroutine dealloc_mo_complex

     subroutine dealloc_mo_quaternion(M)

!      Clean up workspace

       type (qumo), intent(inout) :: M

!      Set the dimensions to zero (not really necessary, but may prevent errors due to abuse of the MO-type)
       M%nao = 0
       M%nmo = 0
       M%nz = 0

!      Free the memory associated with a block of Mos
       deallocate (M%energy)
       deallocate (M%coeff)
       deallocate (M%boson_irrep)
       deallocate (M%index)
       M%allocated = .false.

     end subroutine dealloc_mo_quaternion

#if defined (VAR_MPI)

     subroutine sync_mo_complex(M,my_MPI_master)

      use interface_to_mpi

!     Synchronize and allocate (if not on master) a block of MOcoefficients

      type(cmo), intent(inout) :: M
      integer,   intent(in   ) :: my_MPI_master
      integer                  :: my_MPI_rank
      integer                  :: n_coef

      call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
      call interface_mpi_bcast (M%nao,1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%nmo,1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%total_energy,1,my_MPI_master,global_communicator)
      if (my_MPI_rank /= my_MPI_master) then
         call alloc_mo (M,M%nao,M%nmo)
         M%allocated = .true.
      end if
      call interface_mpi_bcast (M%fermion_irrep,    1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%boson_irrep,  M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%index,        M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%energy,       M%nmo,my_MPI_master,global_communicator)
      n_coef = M%nao * M%nmo * 4 ! (AO,MO,spin) * 2 for complex
      call interface_mpi_bcast (M%coeff,n_coef,my_MPI_master,global_communicator)

     end subroutine sync_mo_complex

     subroutine sync_mo_quaternion(M,my_MPI_master)

      use interface_to_mpi

!     Synchronize and allocate (if not on master) a block of MOcoefficients
 
      type(qumo), intent(inout) :: M
      integer,    intent(in   ) :: my_MPI_master
      integer                   :: my_MPI_rank
      integer                   :: n_coef

      call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
      call interface_mpi_bcast (M%nao,1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%nmo,1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%total_energy,1,my_MPI_master,global_communicator) 
      call interface_mpi_bcast (M%nz ,1,my_MPI_master,global_communicator)
      if (my_MPI_rank /= my_MPI_master) then
         call alloc_mo(M,M%nao,M%nmo,M%nz)
         M%allocated = .true.
      end if
      call interface_mpi_bcast (M%fermion_irrep,    1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%boson_irrep,  M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%index,        M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%energy,       M%nmo,my_MPI_master,global_communicator)
      n_coef = M%nao * M%nmo * M%nz
      call interface_mpi_bcast (M%coeff,n_coef,my_MPI_master,global_communicator)

     end subroutine sync_mo_quaternion

#else

     subroutine sync_mo_complex(M,my_MPI_master)

!     do nothing, makes it possible to call sync also in serial runs

      type(cmo), intent(inout) :: M
      integer,   intent(in   ) :: my_MPI_master

     end subroutine sync_mo_complex

     subroutine sync_mo_quaternion(M,my_MPI_master)

!     do nothing, makes it possible to call sync also in serial runs

      type(qumo), intent(inout) :: M
      integer,    intent(in   ) :: my_MPI_master

     end subroutine sync_mo_quaternion

#endif

     subroutine compress_mo_complex(M,ao_map,mo_map)

!     Reduce or reorder MO coefficients based on AO and MO compression maps

      type(cmo), intent(inout) :: M
      integer,   intent(in), optional   :: ao_map(:,:)
      integer,   intent(in), optional   :: mo_map(:,:)

      integer, allocatable :: map_ao(:,:), map_mo(:,:)
      integer :: i, j, nao_new, nmo_new
      real(8), allocatable    :: energy(:)
      complex(8), allocatable :: coeff(:,:,:)
      integer, allocatable    :: index(:), boson_irrep(:)

      ! start by allocating local copies of the ao and mo maps (with a sanity check)
      if (present(ao_map)) then
         nao_new = size(ao_map,2)
         if (size(ao_map,1) /= 2 .or. nao_new > M%nao) stop 'compress_mo: inconsistent ao_map'
         allocate (map_ao(2,nao_new))
         map_ao = ao_map
      else
         nao_new = M%nao
         allocate (map_ao(2,nao_new))
         forall (i=1:nao_new) map_ao(1:2,i) = i
      end if

      if (present(mo_map)) then
         nmo_new = size(mo_map,2)
         if (size(mo_map,1) /= 2 .or. nmo_new > M%nmo) stop 'compress_mo: inconsistent mo_map'
         allocate (map_mo(2,nmo_new))
         map_mo = mo_map
      else
         nmo_new = M%nmo
         allocate (map_mo(2,nmo_new))
         forall (i=1:nmo_new) map_mo(1:2,i) = i
      end if

      ! additional sanity checks
      if (any(map_ao<1))            stop 'compress_mo: corrupt ao_map'
      if (any(map_ao(1,:)>M%nao))   stop 'compress_mo: corrupt ao_map'
      if (any(map_ao(2,:)>nao_new)) stop 'compress_mo: corrupt ao_map'
      if (any(map_mo<1))            stop 'compress_mo: corrupt mo_map'
      if (any(map_mo(1,:)>M%nmo))   stop 'compress_mo: corrupt mo_map'
      if (any(map_mo(2,:)>nmo_new)) stop 'compress_mo: corrupt mo_map'

      ! resort and store the MO information in the object
      ! we only need to change this if we reorder MOs

      if (present(mo_map)) then

         allocate (energy(nmo_new))
         allocate (boson_irrep(nmo_new))
         allocate (index(nmo_new))

         do i = 1, nmo_new
            energy(map_mo(2,i))      = M%energy(map_mo(1,i))
            boson_irrep(map_mo(2,i)) = M%boson_irrep(map_mo(1,i))
            index(map_mo(2,i))       = M%index(map_mo(1,i))
         enddo

         deallocate (M%energy)
         allocate (M%energy(nmo_new))
         M%energy = energy

         deallocate  (M%boson_irrep)
         allocate (M%boson_irrep(nmo_new))
         M%boson_irrep = boson_irrep

         deallocate (M%index)
         allocate (M%index(nmo_new))
         M%index = index

      end if

      ! resort the coefficients
      allocate (coeff(nao_new,nmo_new,2))
      do j = 1, nmo_new
         do i = 1, nao_new
            coeff(map_ao(2,i),map_mo(2,j),1) = M%coeff(map_ao(1,i),map_mo(1,j),1)
            coeff(map_ao(2,i),map_mo(2,j),2) = M%coeff(map_ao(1,i),map_mo(1,j),2)
         end do
      end do

      ! store the resorted coefficients in the spinor object
      deallocate  (M%coeff)
      allocate (M%coeff(nao_new,nmo_new,2))
      M%coeff = coeff

      ! we can now update the dimension information
      M%nao = nao_new
      M%nmo = nmo_new

     end subroutine compress_mo_complex

     subroutine convert_quaternion_to_complex(coeff,scoeff,eig,seig,nao,nmo,npo)

!     switch to Kramer's unrestricted picture

!     input variables
      real(8),    intent(inout) :: coeff(:,:,:), eig(:)
      real(8),    intent(  out) :: seig(:)
      complex(8), intent(  out) :: scoeff(:,:,:)
      integer,    intent(in   ) :: nao, nmo, npo
!     summation and dummy variables
      integer                 :: i, j, spin

!     create new arrays grouped by Kramer's partners + drop "positronic" solutions

      do i = 1, npo+nmo
        seig(2*i-1) = eig(npo+i)
        seig(2*i  ) = eig(npo+i)
      end do

      do spin = 1, 2
        if (spin == 1) then
            do j = 1, nmo
                do i = 1, nao
                scoeff(i,2*j-1,spin) = dcmplx(coeff(i,j+npo,1),coeff(i,j+npo,2))
                scoeff(i,2*j  ,spin) = dcmplx(coeff(i,j+npo,3),coeff(i,j+npo,4))
                end do
            end do
        else if (spin == 2) then
            do j = 1, nmo
                do i = 1, nao
                scoeff(i,2*j-1,spin) = dcmplx(-coeff(i,j+npo,3), coeff(i,j+npo,4))
                scoeff(i,2*j  ,spin) = dcmplx( coeff(i,j+npo,1),-coeff(i,j+npo,2))
                end do
            end do
        end if
      end do

     end subroutine convert_quaternion_to_complex

     subroutine convert_complex_to_quaternion(cm,qm)

!     switch to Kramer's restricted format, taking unbarred spinors (odd indices) as reference
!     NB: will not work for unrestricted calculations in which unbar/bar are not intertwined exactly

      implicit none
      type(cmo),  intent(inout) :: cm
      type(qumo), intent(inout) :: qm
      integer                   :: i, j

      call alloc_mo (qm,cm%nao,cm%nmo/2,4)
    
      qm%total_energy  = cm%total_energy
      qm%fermion_irrep = cm%fermion_irrep

      do j = 1, qm%nmo
        qm%energy(j) = cm%energy(2*j-1)
        qm%boson_irrep(j) = cm%boson_irrep(2*j-1)
      end do

      do j = 1, qm%nmo
         do i = 1, qm%nao
          qm%coeff(i,j,1) =  real(cm%coeff(i,2*j-1,1))
          qm%coeff(i,j,2) = aimag(cm%coeff(i,2*j-1,1))
          qm%coeff(i,j,3) = -real(cm%coeff(i,2*j-1,2))
          qm%coeff(i,j,4) = aimag(cm%coeff(i,2*j-1,2))
         end do
      end do

     end subroutine convert_complex_to_quaternion

     subroutine apply_Kramers_operator(cm,cmbar)

!     apply Kramers operator: returns time-reversed set of spinors

      implicit none
      type(cmo),  intent(inout) :: cm    ! original set of spinors
      type(cmo), intent(inout)  :: cmbar ! time-reversed set of spinors
      integer                   :: i, j

      if (cmbar%allocated) call dealloc_mo(cmbar)
      call alloc_mo (cmbar,cm%nao,cm%nmo)

      cmbar%total_energy  = cm%total_energy
      cmbar%fermion_irrep = cm%fermion_irrep

      do j = 1, cm%nmo
        cmbar%energy(j)      = cm%energy(j)
        cmbar%index(j)       = cm%index(j)
        cmbar%boson_irrep(j) = cm%boson_irrep(j)
      end do

      do j = 1, cm%nmo
         do i = 1, cm%nao
          cmbar%coeff(i,j,1) = - dconjg(cm%coeff(i,j,2)) ! c_alpha <- - c_beta*
          cmbar%coeff(i,j,2) =   dconjg(cm%coeff(i,j,1)) ! c_beta  <- + c_alpha*
         end do
      end do

     end subroutine apply_Kramers_operator


end module exacorr_mo
