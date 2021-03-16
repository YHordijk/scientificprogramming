module trial_vectors_eom 
  use contraction
  use davidson

 implicit none

 public trial_vectors_ee
 public trial_vectors_ip
 public trial_vectors_ea

 contains

   subroutine trial_vectors_ee (r1,Fbar_mi,Fbar_ae,W_mbej)
#include "symm.inc"
#include "complex.inc"
#include "inpt.inc"
!---------------Description--------------------------------------------
!     generating R1 trial vectors by diagonalizing HBar_SS matrix.
!-----------------------------------------------------------------------

!---------------Calling variables
  real(8), intent(inout)   :: r1(:,:)       ! R-Amplitudes
  real(8), intent(in)  :: Fbar_ae(:),Fbar_mi(:),W_mbej(:)      

!---------------Local ------------------

  real(8), pointer :: sigma1(:,:),evec(:,:), ev(:,:)
  integer          :: i,iroot, nroot,jroot, dimr1

!--------------------------------------------------

    nroot = size(r1,2)
    dimr1 = size(r1,1)

    allocate(sigma1(dimr1,nroot))

    sigma1 = 0.0d0

    do iroot = 1, nroot
! class : (HBar_SS*R)(a,i)
!=========================

!S_ai = S_ai +  F^a_e*R^e_i   
!---------------------------
   
      call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/),Fbar_ae,r1(:,iroot),sigma1(:,iroot), &
  &                       1.0d0,1.0d0,nrep)

!S_ai = S_ai - f_bar_mi(m,i)*r1(a,m) 
!------------------------------------

      call contraction_222 ((/"p1","o2"/),(/"o2","o1"/),(/"p1","o1"/),r1(:,iroot),Fbar_mi, &
  &                         sigma1(:,iroot),-1.0d0,1.0d0,nrep)


!S_ai = S_ai + w_mbej(ma,ei)*r1(e,m)
!------------------------------------ 

      call contraction_422 ((/"o2","p1","p2","o1"/),(/"p2","o2"/),(/"p1","o1"/),r1(:,iroot), &
  &                         sigma1(:,iroot),1.0d0,1.0d0,nrep,W_mbej)

  enddo

   if (iprnt.ge.3) then
       print *,"debug in trial_vectors_ee, nrep:",nrep,dimr1,nroot
       call davidson_print_vectors(sigma1,  label="debug:  sigma1 vector         ")
   end if

   allocate (ev(2,nroot))
   allocate (evec(dimr1,nroot))

   evec = 0.0d0
   ev   = 0.0d0

   call solve_values_vectors(sigma1, ev, evec, evec, .false.)
! note: we are forcing the CCS vectors to be sorted here according to the lowest
! energy, but perhaps we need to perform such sorting at a higher level, at the
! routines generating the different kinds of trial vectors
   call sort_values_vectors(1, ev, evec, evec, .true., .false.,.false., (/ 0.0d0, 0.0d0 /))
!
   if (iprnt.ge.3) then
       print *,'CCS right eigenvalues in (real, imaginary) format'
       do i = 1, nroot
           print *,"(",ev(1,i), ",",ev(2,i),")"
       end do
       call davidson_print_vectors(evec, label="debug: right CCS eigenvectors ")
   end if

! copy the eigenvectors to r1
   r1 = evec

   deallocate(ev)
   deallocate(evec)
   deallocate(sigma1)
   end subroutine


   subroutine trial_vectors_ip (r1,Fbar_mi)
#include "symm.inc"
#include "complex.inc"
#include "inpt.inc"
!---------------Description--------------------------------------------
!     generating R1 trial vectors by diagonalizing HBar_SS matrix.
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------

  real(8), intent(inout)   :: r1(:,:) ! R-Amplitudes
  real(8), intent(in)  :: Fbar_mi(:)      

!---------------Local Variables--------------------

  real(8), pointer :: sigma1(:,:),evec(:,:), ev(:,:)
  integer              :: i,iroot,nroot,jroot,dimr1, j, k, myrep

!--------------------------------------------------

    nroot = size(r1,2)
    dimr1 = size(r1,1)

    allocate(sigma1(dimr1,nroot))

    sigma1 = 0.0d0

   do iroot = 1, nroot

     call contraction_222 ((/"o2","o1"/),(/"c1","o2"/),(/"c1","o1"/),Fbar_mi,r1(:,iroot),sigma1(:,iroot), &
    &                       -1.0d0,1.0d0,nrep)

   enddo
   if (iprnt.ge.3) then
       print *,"debug in trial_vectors_ip, nrep:",nrep,dimr1,nroot
       call davidson_print_vectors(sigma1,  label="debug:  sigma1 vector         ")
   end if


   allocate (ev(2,nroot))
   allocate (evec(dimr1,nroot))

   evec = 0.0d0
   ev   = 0.0d0

   call solve_values_vectors(sigma1, ev, evec, evec, .false.)
! note: we are forcing the CCS vectors to be sorted here according to the lowest
! energy, but perhaps we need to perform such sorting at a higher level, at the
! routines generating the different kinds of trial vectors
   call sort_values_vectors(1, ev, evec, evec, .true., .false.,.false., (/ 0.0d0, 0.0d0 /))
   if (iprnt.ge.3) then
       print *,'CCS eigenvalues in (real, imaginary) format'
       do i = 1, nroot
           print *,"(",ev(1,i), ",",ev(2,i),")"
       end do
       call davidson_print_vectors(evec, label="debug: right CCS eigenvectors ")
   end if

! copy the eigenvectors to r1
   r1 = evec

   deallocate(ev)
   deallocate(evec)
   deallocate(sigma1)

 end subroutine

   subroutine trial_vectors_EA (r1,Fbar_ae)
#include "symm.inc"
#include "complex.inc"
#include "inpt.inc"

!---------------Description--------------------------------------------
!     generating R1 trial vectors by diagonalizing HBar_SS matrix.
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------

  real(8), intent(inout)   :: r1(:,:) ! R-Amplitudes
  real(8), intent(in)  :: Fbar_ae(:)      

!---------------Local Variables--------------------

  real(8), pointer :: sigma1(:,:),evec(:,:), ev(:,:)
  integer          :: i,iroot,nroot,jroot,dimr1

!--------------------------------------------------

    nroot = size(r1,2)
    dimr1 = size(r1,1)

    allocate(sigma1(dimr1,nroot))

    sigma1 = 0.0d0

   do iroot = 1, nroot

    call contraction_222 ((/"p1","p2"/),(/"p2","c1"/),(/"p1","c1"/),Fbar_ae,r1(:,iroot),sigma1(:,iroot), &
   &                    1.0d0,1.0d0,nrep)

   end do

   if (iprnt.ge.3) then
       print *,"debug in trial_vectors_ea, nrep:",nrep,dimr1,nroot
       call davidson_print_vectors(sigma1,  label="debug:  sigma1 vector         ")
   end if


   allocate (ev(2,nroot))
   allocate (evec(dimr1,nroot))

   evec = 0.0d0
   ev   = 0.0d0

   call solve_values_vectors(sigma1, ev, evec, evec, .false.)
! note: we are forcing the CCS vectors to be sorted here according to the lowest
! energy, but perhaps we need to perform such sorting at a higher level, at the
! routines generating the different kinds of trial vectors
! also, we are hardcoding not using any shift
   call sort_values_vectors(1, ev, evec, evec, .true., .false., .false., (/ 0.0d0, 0.0d0 /))
   if (iprnt.ge.3) then
       print *,'CCS eigenvalues in (real, imaginary) format'
       do i = 1, nroot
           print *,"(",ev(1,i), ",",ev(2,i),")"
       end do
       call davidson_print_vectors(evec, label="debug: right CCS eigenvectors ")
   end if

! copy the eigenvectors to r1
   r1 = evec

   deallocate(ev)
   deallocate(evec)
   deallocate(sigma1)


 end subroutine

end module 
