!==============================================================================!
module type_laplace
!-------------------------------------------------------------
! Module holding information about Laplace transformation 
! of energy denominators
!-------------------------------------------------------------

 implicit none

 private 

!-------------------------------------------------!
! type laplace
!-------------------------------------------------!
 type laplace_info

  ! Logical variable, that is true, if laplace is used:
  logical :: used = .false.

  ! number of Laplace exponents:
  integer :: num = 0

  ! upper and lower bound of orbital energy denominator
  real(8) :: bounds(2) = (/ 0.0d0, 0.0d0 /)

  ! maximum Laplace error:
  real(8) :: errmax = 0.0d0

  ! RMSD Laplace error:
  real(8) :: errmsd = 0.0d0

  ! Laplace exponents:
  real(8), dimension(:),   pointer :: val => null()

  ! Laplace weights:
  real(8), dimension(:),   pointer :: wei => null()

 end type laplace_info


 public :: laplace_info, allocate_laplace, deallocate_laplace

 contains

!-------------------------------------------------!
! Constructor for type laplace
!-------------------------------------------------!
  subroutine allocate_laplace(laplace,num)

   implicit none

   integer, intent(in) :: num
   type(laplace_info)       :: laplace

   integer             :: ierr

   laplace%num = num

   nullify(laplace%val)
   nullify(laplace%wei)

   allocate(laplace%val(num),laplace%wei(num),stat=ierr)
   if (ierr.ne.0) call quit('Memory allocation error in allocate_laplace')

   laplace%val(:) = 0.0d0
   laplace%wei(:) = 0.0d0

  end subroutine allocate_laplace

!-------------------------------------------------!
! Destructor for type laplace
!-------------------------------------------------!
  subroutine deallocate_laplace(laplace)

   implicit none

   type(laplace_info)  :: laplace
   integer             :: ierr

   ierr=0
   if (associated(laplace%val)) deallocate(laplace%val,stat=ierr)
   if (ierr.eq.0.and.associated(laplace%wei)) deallocate(laplace%wei,stat=ierr)
   if (ierr.ne.0) call quit('Memory deallocation error in deallocate_laplace')

  end subroutine deallocate_laplace


end module type_laplace
!==============================================================================!
