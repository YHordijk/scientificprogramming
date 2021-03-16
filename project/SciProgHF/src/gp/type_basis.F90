!==============================================================================!
module type_basis
!------------------------------------------------------------------------------!
!
! This module stores information on basis sets in a derived type
!
!   BASIS_INFO
!
! .  This derived type can created and destroyed by
!
!   ALLOCATE_BASIS     and    DEALLOCATE_BASIS
!
!
!
! kmax                     - number of shells contracted shells
! nuco(1:kmax)             - number of uncontracted functions
! nrco(1:????)             - number of contracted functions
! nhkt(1:kmax)             - angular quantum number (s=1,p=2,d=3 etc.)
! khkt(1:kmax)             - number of spherical (Cartesian) components
! kckt(1:kmax)             - number of Cartesian components
! ncent(1:kmax)            - index of symmetry independent center
! numcf(1:kmax)            - index of shell in AO-block
! jstrt(1:kmax)            - offset for in primitive shell
! ishst(1:kmax)            - offset for in contracted shell
!
! npshel                   - number of contracted shells
! ncont                    - number of contractions
! priexp(1:npshel)         - normalized contraction coefficients
! priccf(1:npshel,1:ncont) - exponents of primitive shells
!
!  Note:
!
!   (1) Please create an instance of the derived type for each kind of basis set
!       (orbital basis, fitting basis (Coulomb fit, MP2 fit ...), 
!        CABS basis  for F12 theories ...)
!
!   (2) The module basically contains information provided on the 
!
!              SHELL.H     and    PRIMIT.H
!
!       common blocks. Those common blocks should be removed at some point
!       and replaced by such containers. The advantage is point (1)
!
!   (3) If you feel that something should be added, please go ahead, but mind
!       that simplicity wins in the long run!
! 
!   (4) So far, this module is used only in the LBM-AOCC and AO-NEVPT2 branch.
!       I would be glad if others would use as well.
!
! 
!  BHP, summer 2016
!
!------------------------------------------------------------------------------!

 implicit none

 private 

 type basis_info

  integer :: kmax   = 0
  integer :: npshel = 0
  integer :: ncont  = 0

  integer, dimension(:), pointer :: nuco(:) => null()
  !integer, dimension(:), pointer :: nrco(:) => null()
  integer, dimension(:), pointer :: nhkt(:) => null()
  integer, dimension(:), pointer :: khkt(:) => null()
  integer, dimension(:), pointer :: kckt(:) => null()
  integer, dimension(:), pointer :: ncent(:) => null()
  integer, dimension(:), pointer :: numcf(:) => null()
  integer, dimension(:), pointer :: jstrt(:) => null()
  integer, dimension(:), pointer :: ishst(:) => null()

  real(8), dimension(:), pointer :: priexp(:) => null()
  real(8), dimension(:), pointer :: priccf(:,:) => null()

 end type basis_info

 public :: basis_info, allocate_basis, deallocate_basis

 contains

 subroutine allocate_basis(basis,kmax,npshel,ncont)

  implicit none

  integer, intent(in) :: kmax,npshel,ncont
  type(basis_info), intent(inout) :: basis

  basis%kmax   = kmax
  basis%npshel = npshel
  basis%ncont  = ncont

  allocate(basis%nuco(kmax),&
           basis%nhkt(kmax),&
           basis%khkt(kmax),&
           basis%kckt(kmax),&
           basis%ncent(kmax),&
           basis%numcf(kmax),&
           basis%jstrt(kmax),&
           basis%ishst(kmax+1),&
           basis%priexp(npshel),&
           basis%priccf(npshel,ncont))

 end subroutine allocate_basis

 subroutine deallocate_basis(basis)

  implicit none

  type(basis_info), intent(inout) :: basis

  if (associated(basis%nuco))    deallocate(basis%nuco)
  if (associated(basis%nhkt))    deallocate(basis%nhkt)
  if (associated(basis%khkt))    deallocate(basis%khkt)
  if (associated(basis%kckt))    deallocate(basis%kckt)
  if (associated(basis%ncent))   deallocate(basis%ncent)
  if (associated(basis%numcf))   deallocate(basis%numcf)
  if (associated(basis%jstrt))   deallocate(basis%jstrt)
  if (associated(basis%ishst))   deallocate(basis%ishst)
  if (associated(basis%priexp))  deallocate(basis%priexp)
  if (associated(basis%priccf))  deallocate(basis%priccf)

  nullify(basis%nuco,basis%nhkt,basis%khkt,basis%kckt,&
          basis%ncent,basis%numcf,basis%jstrt,basis%ishst,&
          basis%priexp,basis%priccf)

  basis%kmax   = 0
  basis%npshel = 0
  basis%ncont  = 0

 end subroutine deallocate_basis

end module type_basis
!==============================================================================!
