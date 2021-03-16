module sigma_eom

  use contraction
  use interface_to_mpi
  use intermediates_1b_2b
  implicit none

!this module calculates the sigma vectors for IP(2h-1p), EA(1h-2p) and EE(2h-2p)  

!    public eom_ip  
!    public eom_ea
    public sigma_ee_1h1p
    public sigma_IP
    public lambda_IP
    public lambda_EE
    public sigma_EA
    public IM_container
    public create_sigma_vector_left
    public create_sigma_vector_right
    public create_diagonal
    public create_eom_ndet
    public initialize_buffers
    public free_buffers
    public free_intermediates
    public free_eom_ndet
    public get_matrix_size
    public set_eom_ndet

    private
    public intermediates
    type intermediates
      character(8) :: eomtype
      real(8),allocatable,dimension(:) :: Fbar_mi, Fbar_me, Fbar_ae                                         ! all 1-b intermediates 
      real(8),allocatable,dimension(:) :: W_ijmn, W_ejmb, W_mbej, W_iemn, W_efam, W_amef, W_mnie, Wbar_mbej ! all 2-b intermediates   
      real(8),allocatable,dimension(:) :: t1, t2
      real(8),allocatable,dimension(:) :: H_II
    endtype intermediates

    real(kind=8), pointer :: local_buf1(:) 
    real(kind=8), pointer :: local_buf2(:) 
    real(kind=8), pointer :: local_buf3(:) 
    integer :: nbuf3
    integer, pointer :: eom_ndet(:)

   contains


   subroutine create_eom_ndet(max_excitation_level)
       integer, intent(in) :: max_excitation_level

       if (.not.associated(eom_ndet)) allocate(eom_ndet(max_excitation_level))
       eom_ndet = 0
   end subroutine

   subroutine free_eom_ndet()
       if (associated(eom_ndet)) deallocate(eom_ndet)
!       nullify(eom_ndet)
   end subroutine

   subroutine set_eom_ndet(excitation_level, nr_det)
       integer, intent(in) :: excitation_level
       integer, intent(in) :: nr_det
       eom_ndet(excitation_level) = nr_det
   end subroutine 

   function get_matrix_size(A)
      integer :: get_matrix_size
      type(intermediates), intent(in) :: A

      get_matrix_size = eom_ndet(1) + eom_ndet(2)
   end function

   subroutine free_intermediates(A)
      use memory_allocator
      use allocator_parameters, only : klongint
      type(intermediates), intent(inout) :: A

      call dealloc(A%Fbar_mi,   id="Fbar_mi")
      call dealloc(A%Fbar_me,   id="Fbar_me")
      call dealloc(A%Fbar_ae,   id="Fbar_ae")
    
      call dealloc(A%W_ijmn,    id="W_ijmn")
      call dealloc(A%W_ejmb,    id="W_ejmb")
      call dealloc(A%W_mbej,    id="W_mbej")
      call dealloc(A%W_iemn,    id="W_iemn")

      call dealloc(A%W_efam,    id="W_efam")
      call dealloc(A%W_amef,    id="W_amef")

      call dealloc(A%W_mnie,    id="W_mnie")
      call dealloc(A%Wbar_mbej, id="Wbar_mbej")

      call dealloc(A%t1,  id="t1")
      call dealloc(A%t2,  id="t2")

!      if (allocated(A%H_II)) deallocate(A%H_II)

   end subroutine

   subroutine initialize_buffers(buf1, buf2, buf3)
       real(kind=8), allocatable, target ::  buf1(:), buf2(:), buf3(:)
       local_buf1 => buf1 
       local_buf2 => buf2 
       local_buf3 => buf3 
   end subroutine

   subroutine free_buffers(buf1, buf2, buf3)
      use memory_allocator
      use allocator_parameters, only : klongint

      real(kind=8),allocatable,intent(inout) ::  buf1(:), buf2(:), buf3(:)

      call dealloc(buf1, id="buf1")
      if (associated(local_buf1)) then
          nullify(local_buf1)
      end if

      call dealloc(buf2, id="buf2")
      if (associated(local_buf2)) then
          nullify(local_buf2)
      end if

      call dealloc(buf3, id="buf3")
      if (associated(local_buf3)) then
          nullify(local_buf3)
      end if
   end subroutine

   function create_diagonal(A)
      real(kind=8), pointer :: create_diagonal(:)
      type(intermediates) :: A
      integer :: N, N1, N2

      N1 = eom_ndet(1)
      N2 = eom_ndet(2)

      N = N1 + N2

      allocate(create_diagonal(N))
      create_diagonal = 1.0d0

      select case (A%eomtype)
         case ("IP")
             create_diagonal = A%H_II 
         case ("EA")
             create_diagonal = A%H_II 
         case ("EE")
             create_diagonal = A%H_II 
      end select

   end function


   function create_sigma_vector_left(A, r, old_sigma_vector,state_sym)
#include "files.inc"
#include "inpt.inc"
      real(kind=8), pointer :: create_sigma_vector_left(:,:)
      real(kind=8), pointer :: old_sigma_vector(:,:)
      type(intermediates) :: A
      real(kind=8), intent(in) :: r(:,:)
      integer :: L, N, L_already_done
      integer :: N1, N2
      integer :: start_r1, end_r1, start_r2, end_r2
      integer :: iroot
      integer, intent(in) :: state_sym

      N = size(r,1)
      L = size(r,2)

      allocate(create_sigma_vector_left(N,L))

      if (associated(old_sigma_vector)) then 
          L_already_done = size(old_sigma_vector,2) 
          do iroot = 1, L_already_done
              create_sigma_vector_left(:,iroot) = old_sigma_vector(:,iroot)
          end do
          deallocate(old_sigma_vector)
      else
          L_already_done = 0
      endif

      N1 = eom_ndet(1)
      N2 = eom_ndet(2)

      start_r1 = 1
      end_r1   = N1
      start_r2 = end_r1 + 1
      end_r2   = end_r1 + N2

      select case (A%eomtype)
         case ('IP')

!$OMP PARALLEL IF((L_already_done - L) .GT. 1)
            call show_number_of_threads(iw)
!$OMP DO
             do iroot = L_already_done + 1, L 

                create_sigma_vector_left(:,iroot:iroot) = 0.0d0
                if (iprnt.ge.2) write(iw,"(2x,a,i6,a,i6)") "debug, constructing Left EOM-IP sigma vector ",iroot," out of ",L
                call Lambda_IP(r(start_r1:end_r1, iroot), &
                 &             r(start_r2:end_r2, iroot), &
                 &             A,                               &
                 &             create_sigma_vector_left(start_r1:end_r1, iroot), &
                 &             create_sigma_vector_left(start_r2:end_r2, iroot))
             enddo
!$OMP END DO
!$OMP END PARALLEL

         case ("EA")

!$OMP PARALLEL IF((L_already_done - L) .GT. 1)
            call show_number_of_threads(iw)
!$OMP DO
             do iroot = L_already_done + 1, L 

                create_sigma_vector_left(:,iroot:iroot) = 0.0d0
                if (iprnt.ge.2) write(iw,"(2x,a,i6,a,i6)") "debug, constructing Left EOM-EA sigma vector ",iroot," out of ",L
!               call Lambda_EA(r(start_r1:end_r1, iroot:iroot),  &
!                &             r(start_r2:end_r2, iroot:iroot), &
!                &             A,                               &
!                &             create_sigma_vector_left(start_r1:end_r1, iroot:iroot), &
!                &             create_sigma_vector_left(start_r2:end_r2, iroot:iroot), &
!                &             nbuf3) 

             enddo
!$OMP END DO
!$OMP END PARALLEL 

         case ('EE')

!$OMP PARALLEL IF((L_already_done - L) .GT. 1)
            call show_number_of_threads(iw)
!$OMP DO
             do iroot = L_already_done + 1, L

                create_sigma_vector_left(:,iroot:iroot) = 0.0d0
                if (iprnt.ge.2) write(iw,"(2x,a,i6,a,i6)") "debug, constructing Left EOM-EE sigma vector ",iroot," out of ",L
                call Lambda_ee(r(start_r1:end_r1, iroot), &
                 &                  r(start_r2:end_r2, iroot), & 
                 &                  A,                               &
                 &                  create_sigma_vector_left(start_r1:end_r1, iroot), &
                 &                  create_sigma_vector_left(start_r2:end_r2, iroot), state_sym)

             enddo
!$OMP END DO
!$OMP END PARALLEL 
      end select


   end function 
   

   function create_sigma_vector_right(A, r, old_sigma_vector, state_sym)
#include "files.inc"
#include "inpt.inc"
      real(kind=8), pointer :: create_sigma_vector_right(:,:)
      real(kind=8), pointer :: old_sigma_vector(:,:)
      type(intermediates) :: A
      real(kind=8), intent(in) :: r(:,:)
      integer :: L, N, L_already_done, N1, N2
      integer :: start_r1, end_r1, start_r2, end_r2
      integer :: iroot,ierr
      integer, intent(in) :: state_sym

      N = size(r,1)
      L = size(r,2)

      allocate(create_sigma_vector_right(N,L))

      if (associated(old_sigma_vector)) then 
          L_already_done = size(old_sigma_vector,2) 
          do iroot = 1, L_already_done
              create_sigma_vector_right(:,iroot) = old_sigma_vector(:,iroot)
          end do
          deallocate(old_sigma_vector)
      else
          L_already_done = 0
      endif

      N1 = eom_ndet(1)
      N2 = eom_ndet(2)

      start_r1 = 1
      end_r1   = N1
      start_r2 = end_r1 + 1
      end_r2   = end_r1 + N2

      select case (A%eomtype)
         case ('IP')

!$OMP PARALLEL IF((L_already_done - L) .GT. 1)
            call show_number_of_threads(iw)
!$OMP DO
             do iroot = L_already_done + 1, L 

                create_sigma_vector_right(:,iroot:iroot) = 0.0d0
                if (iprnt.ge.2) write(iw,"(2x,a,i6,a,i6)") "debug, constructing EOM-IP sigma vector ",iroot," out of ",L
                call sigma_IP(r(start_r1:end_r1, iroot), &
                 &             r(start_r2:end_r2, iroot),A, &
                 &             create_sigma_vector_right(start_r1:end_r1, iroot), &
                 &             create_sigma_vector_right(start_r2:end_r2, iroot))
             enddo
!$OMP END DO
!$OMP END PARALLEL 

         case ("EA")

!$OMP PARALLEL IF((L_already_done - L) .GT. 1)
            call show_number_of_threads(iw)
!$OMP DO
             do iroot = L_already_done + 1, L 
!               print *,'thread number',omp_get_thread_num(),'root',iroot
                create_sigma_vector_right(:,iroot:iroot) = 0.0d0
                if (iprnt.ge.2) write(iw,"(2x,a,i6,a,i6)") "debug, constructing EOM-EA sigma vector ",iroot," out of ",L
                call sigma_EA(r(start_r1:end_r1, iroot), &
                 &             r(start_r2:end_r2, iroot),A, &
                 &             create_sigma_vector_right(start_r1:end_r1, iroot), &
                 &             create_sigma_vector_right(start_r2:end_r2, iroot)) 

             enddo
!$OMP END DO
!$OMP END PARALLEL 

         case ('EE')

!$OMP PARALLEL IF((L_already_done - L) .GT. 1)
            call show_number_of_threads(iw)
!$OMP DO
             do iroot = L_already_done + 1, L

                create_sigma_vector_right(:,iroot:iroot) = 0.0d0
                if (iprnt.ge.2) write(iw,"(2x,a,i6,a,i6)") "debug, constructing EOM-EE sigma vector ",iroot," out of ",L
                call sigma_ee_1h1p(r(start_r1:end_r1, iroot), &
                 &                  r(start_r2:end_r2, iroot), & 
                 &                  A,                               &
                 &                  create_sigma_vector_right(start_r1:end_r1, iroot), &
                 &                  create_sigma_vector_right(start_r2:end_r2, iroot),state_sym)

             enddo
!$OMP END DO
!$OMP END PARALLEL 
      end select


   end function

   subroutine IM_container(eomtype,A)

#include "symm.inc"
#include "param.inc"
#include "complex.inc"
#include "ccpar.inc"
#include "files.inc"
#include "inpt.inc"

!-----------------------------------------------------------------------------   
!  In this routine we calculate all the intermediates required for a specific
!  problem (IP/EA/EE). 
!-----------------------------------------------------------------------------

!-------------------------------Variables--------------------------------------

!Calling Variables
!----------------- 
character(8), intent (in) :: eomtype  
type(intermediates), intent(inout) :: A

!Local Variables
!----------------
real(8), allocatable :: foo(:),fvo(:),fvv(:) 
integer :: ierr

!-----------------------------------------------------------------------------

     nbuf3 = size(local_buf3)

     allocate (foo(nfoo*rcw))
     allocate (fvo(nfvo*rcw))
     allocate (fvv(nfvv*rcw))

     A%eomtype = eomtype

     write (iw,*) " "
     write (iw,*) " constructing intermediates for EOM-",eomtype

     call fmtofile (.false.,fvo,foo,fvv) ! read one-body integrals from file.

     call f_bar_mi(foo,fvo,A%t1,A%t2,A%Fbar_mi)
     if (iprnt.ge.2) then 
        write (iw,*) "    done with Fbar_mi",dot_product(A%Fbar_mi,A%Fbar_mi)
     else
        write (iw,*) "    done with Fbar_mi"
     end if

     call f_bar_me(fvo,A%t1,A%Fbar_me)
     if (iprnt.ge.2) then 
         write (iw,*) "    done with Fbar_me",dot_product(A%Fbar_me,A%Fbar_me)       
     else
         write (iw,*) "    done with Fbar_me"
     end if

     call f_bar_ae(fvv,fvo,A%t1,A%t2,A%Fbar_ae)
     if (iprnt.ge.2) then 
         write (iw,*) "    done with Fbar_ae",dot_product(A%Fbar_ae,A%Fbar_ae)
     else
         write (iw,*) "    done with Fbar_ae"
     end if

     call w_ijmn(A%t1,A%t2,A%W_ijmn)
     if (iprnt.ge.2) then 
         write (iw,*) "    done with W_ijmn",dot_product(A%W_ijmn,A%W_ijmn)
     else
         write (iw,*) "    done with W_ijmn"
     end if

     call wbar_mbej(A%t2,A%Wbar_mbej)
     if (iprnt.ge.2) then 
         write (iw,*) "    done with Wbar_mbej",dot_product(A%Wbar_mbej,A%Wbar_mbej)
     else
         write (iw,*) "    done with Wbar_mbej"
     end if

     call w_mbej(A%t1,A%t2,A%W_mbej)
     if (iprnt.ge.2) then 
         write (iw,*) "    done with W_mbej",dot_product(A%W_mbej,A%W_mbej)
     else
         write (iw,*) "    done with W_mbej"
     end if

     call w_ejmb(A%t1,A%t2,A%W_ejmb)
     if (iprnt.ge.2) then 
         write (iw,*) "    done with W_ejmb",dot_product(A%W_ejmb,A%W_ejmb)
     else
         write (iw,*) "    done with W_ejmb"
     end if

     call w_iemn(A%t1,A%t2,A%W_iemn,A%Wbar_mbej,A%W_ijmn,A%Fbar_me)
     if (iprnt.ge.2) then 
         write (iw,*) "    done with W_iemn",dot_product(A%W_iemn,A%W_iemn)
     else
         write (iw,*) "    done with W_iemn"
     end if

     if (eomtype /= "IP") then

!      CALL XTIME(5,1,'-- "efam" Intermediate                     ')
         call w_efam(A%t1,A%t2,A%Wbar_mbej,A%Fbar_me,A%W_efam)
         if (iprnt.ge.2) then 
             write (iw,*) "    done with W_efam",dot_product(A%W_efam,A%W_efam)
         else
             write (iw,*) "    done with W_efam"
         end if

         if (myproc.gt.0) then
             A%W_efam(1:nv5*rcw) = 0.0d0
         endif

         call w_efam_diagram2 (A%t1,A%t2,local_buf2,local_buf3,nbuf3,A%W_efam)
         if (iprnt.ge.2) then 
             write (iw,*) "    done with W_efam diagram #2",dot_product(A%W_efam,A%W_efam)
         else
             write (iw,*) "    done with W_efam diagram #2"
         end if

!   synchronize all the w_efam contributions from different nodes to the master 

#if defined (VAR_MPI)
      if (nmproc .gt. 1) then

        call xcopy (nv5,a0,0,local_buf3,1)
        ierr = 0

        call interface_mpi_allreduce_r1_work_f77(A%W_efam,local_buf3(1), &
            rcw*ivovvt(nrep+1), &
            op_mpi_sum,global_communicator)

        if(ierr.gt.0) then
           call quit('mpi_reduce error in collecting w_efam !')
        endif
        call xcopy(ivovvt(nrep+1),local_buf3,1,A%W_efam,1)

      endif
         if (iprnt.ge.2) then
             write (iw,*) "  Sync done   W_efam diagram #2",dot_product(A%W_efam,A%W_efam)
         else
             write (iw,*) "  Sync done   W_efam diagram #2"
         end if

#endif
!      CALL XTIME(5,2,'-- "efam" Intermediate                     ')
     endif 

     call w_mnie(A%t1,A%W_mnie)
     if (iprnt.ge.2) then
         write (iw,*) "    done with W_mnie",dot_product(A%W_mnie,A%W_mnie)
     else
         write (iw,*) "    done with W_mnie"
     end if


     if (eomtype /= "IP") then

         call w_amef(A%t1,A%W_amef)
         if (iprnt.ge.2) then
             write (iw,*) "    done with W_amef",dot_product(A%W_amef,A%W_amef)
         else
             write (iw,*) "    done with W_amef"
         end if

     endif 

   deallocate(foo)
   deallocate(fvo)
   deallocate(fVV)


   end subroutine

   subroutine sigma_ee_1h1p(r1,r2,B,sigma1,sigma2,state_sym)
   use symmetry_offset
   use modified_sorting
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "ccpar.inc"
#include "files.inc"
#include "inpt.inc"

!---------------Description--------------------------------------------
!     Calculates 1h1p & 2h2p type sigma vectors for excitation energy
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------
  real(8), intent(in),contiguous,target   :: r1(:), r2(:) ! R-Amplitudes
  type(intermediates), intent(in)  :: B      
  real(8), intent(inout)           :: sigma1(:), sigma2(:)  ! sigma-vectors
  integer, intent(in) :: state_sym
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
  real(8), allocatable :: int_1b(:),sigma2_temp(:),lambda2_temp(:),w_vvvv(:),local_buffer(:) 
  integer :: off1,off2,irp,jrp,istart,mint,m,n,k,nint
  real(8), allocatable :: sigma2_local(:) 
  integer :: off3,i,j,kstart,ntot,ierr,mrp,nrp,ndimr1,ndimr2, N1, N2
  logical :: done 
  real(8)              :: ddot
  type(Offset) :: f,g
  integer  ::row(nrep), column(nrep)
!--------------------------------------------------

  N1 = size(sigma1,1)
  N2 = size(sigma2,1)

  ndimr1 = N1/rcw
  ndimr2 = N2/rcw

  sigma2(1:N2) = 0.0d0

! class : (HBar_SS*R)(a,i)
!=========================

!S_ai = S_ai +  F^a_e*R^e_i   
!---------------------------
   
    call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/),B%Fbar_ae,r1,sigma1,1.0d0,1.0d0,nrep, &
       irrep_left=1,irrep_right=state_sym)

    if (iprnt.ge.2) &
    write(iw,*)'diagram 1 ', dot_product(sigma1,sigma1)

!S_ai = S_ai - f_bar_mi(m,i)*r1(a,m) 
!------------------------------------

    call contraction_222 ((/"p1","o2"/),(/"o2","o1"/),(/"p1","o1"/),r1,B%Fbar_mi,sigma1,-1.0d0,1.0d0,nrep, &
        irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.2) &
    write(iw,*)'diagram 2 ', dot_product(sigma1,sigma1)

!S_ai = S_ai + w_mbej(ma,ei)*r1(e,m)
!------------------------------------ 

 !!    call contraction_422 ((/"o2","p1","p2","o1"/),(/"p2","o2"/),(/"p1","o1"/),r1,sigma1,1.0d0,1.0d0,nrep,B%W_mbej, &
 !!                       irrep_left=1,irrep_right=state_sym)

       call contraction_422 ((/"p1","o2","o1","p2"/),(/"p2","o2"/),(/"p1","o1"/),r1,sigma1,1.0d0,1.0d0,nrep,B%W_ejmb, &
                  irrep_left=1,irrep_right=state_sym)

!   call contraction_242 ((/"p2","o2"/),(/"o2","p1","p2","o1"/),(/"p1","o1"/),r1,sigma1,1.0d0,1.0d0,nrep,B%W_mbej, &
!                      irrep_left=state_sym,irrep_right=1)

    if (iprnt.ge.2) &
   write(iw,*)'diagram 3 ', dot_product(sigma1,sigma1)

!!class : [HBar_SD*R](a,i)
!!========================

!!S_ai = S_ai + F_bar_me(m,e)*r2(ae,im)
!!-------------------------------------

     call contraction_422 ((/"p1","p2","o1","o2"/),(/"o2","p2"/), &
    & (/"p1","o1"/),B%Fbar_me,sigma1,1.0d0,1.0d0,nrep,LeftTensor=r2,irrep_left=state_sym,irrep_right=1)

    if (iprnt.ge.2) &
    write(iw,*)'diagram 4 ',dot_product(sigma1,sigma1)

!!S_ai = S_ai + w_amef(am,ef) * r2(ef,im)
!!---------------------------------------
 
    call contraction_442((/"p1","o2","p2","p3"/),(/"p2","p3","o1","o2"/),(/"p1","o1"/),sigma1,1.0d0,1.0d0,nrep, &
   &         LeftTensor=B%W_amef,RightTensor=r2,irrep_left=1,irrep_right=state_sym)

    if (iprnt.ge.2) &
    write(iw,*)'diagram 5 ',dot_product(sigma1,sigma1)
!!S_ai = S_ai - W_mnie(mn,ie)*r2(ae,mn)
!!--------------------------------------

     call contraction_442((/"p1","p2","o2","o3"/),(/"o2","o3","o1","p2"/),(/"p1","o1"/),sigma1,-1.0d0,1.0d0,nrep, &
    & RightTensor=B%W_mnie,LeftTensor=r2,irrep_left=state_sym,irrep_right=1)

    if (iprnt.ge.2) &
    write(iw,*)'diagram 6 ', dot_product(sigma1,sigma1)

! class : [HBar_DS*R](ab,ij)
! ==========================

!!S(ab,ij) = S(ab,ij) - P(ab) w_iemn(mb,ij)*r1(a,m)
!!-------------------------------------------------

     call contraction_244 ((/"p1","o3"/),(/"o3","p2","o1","o2"/), &
    &    (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,righttensor=b%w_iemn, &
    &   LeftTensor=r1,irrep_left=state_sym,irrep_right=1) 

    if (iprnt.ge.2) &
    write(iw,*)'diagram 7 ', dot_product(sigma2,sigma2)

!!S(ab,ij) = S(ab,ij) + P(ij) w_efam(ab,ej)*r1(e,i)
!!-------------------------------------------------

    call contraction_424 ((/"p1","p2","p3","o2"/),(/"p3","o1"/), &
   &    (/"p1","p2","o2","o1"/),sigma2,1.0d0,-1.0d0,nrep,LeftTensor=B%W_efam, &
   &    RightTensor=r1,irrep_left=1,irrep_right=state_sym) 

!   call contraction_424 ((/"p1","p2","o2","p3"/),(/"p3","o1"/), &
!  &    (/"p1","p2","o2","o1"/),sigma2,1.0d0,-1.0d0,nrep, &
!  &    RightTensor=r1,irrep_left=1,irrep_right=state_sym) 

!   call contraction_244 ((/"p3","o1"/),(/"p1","p2","p3","o2"/), &
!  &    (/"p1","p2","o2","o1"/),sigma2,1.0d0,-1.0d0,nrep,LeftTensor=r1, &
!  &    RightTensor=B%W_efam,irrep_right=1,irrep_left=state_sym) 

       if (iprnt.ge.2) &
    write(iw,*)'diagram 8 ', dot_product(sigma2,sigma2)

! int_1b(b,f) = int_1b(b,f) + w_amef (bm,fe)*r1(e,m)
! S(ab,ij) = S(ab,ij) + P(ab) int_1b(b,f)*t2(af,ij)
!---------------------------------------------------

      call alloc_array(f,nrep)

      call auto_symmetry_offset(f,nv,nv,.false.,.false.)

    allocate(int_1b(f%oneDirac(state_sym)*rcw))

    int_1b = 0.0d0

   call contraction_422 ((/"p2","o3","p3","p4"/),(/"p3","o3"/),(/"p2","p4"/),r1,int_1b,-1.0d0,1.0d0,nrep,B%W_amef, &
          irrep_left=1,irrep_right=state_sym)

!   call contraction_424 ((/"p1","p4","o1","o2"/),(/"p2","p4"/),(/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0, &
!  &                       nrep,RightTensor=int_1b,LeftTensor=B%t2,irrep_left=1,irrep_right=state_sym)

    call contraction_244 ((/"p2","p4"/),(/"p4","p1","o1","o2"/),(/"p2","p1","o1","o2"/),sigma2,-1.0d0,-1.0d0, &
   &                       nrep,LeftTensor=int_1b,RightTensor=B%t2,irrep_left=state_sym,irrep_right=1)

    deallocate(int_1b)
      call dealloc_array(f)

       if (iprnt.ge.2) &
    write(iw,*)'diagram 9 ',dot_product(sigma2,sigma2)

! int_1b(n,j) = int_1b(n,j) + w_mnie(nm,je)*r1(e,m)
! S(ab,ij) = S(ab,ij) - int_1b(n,j) * t2(ab,in)
!---------------------------------------------------


      call alloc_array(f,nrep)

      call auto_symmetry_offset(f,no,no,.false.,.false.)

    allocate(int_1b(f%oneDirac(state_sym)*rcw))

    int_1b = 0.0d0

    call contraction_422 ((/"o4","o3","o2","p3"/),(/"p3","o3"/),(/"o4","o2"/),r1,int_1b,1.0d0,1.0d0,nrep,LeftTensor=B%W_mnie, &
        irrep_left=1,irrep_right=state_sym)

    call contraction_424 ((/"p1","p2","o1","o4"/),(/"o4","o2"/), &
   & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,RightTensor=int_1b,LeftTensor=B%t2,irrep_left=1,irrep_right=state_sym)

    deallocate(int_1b)

      call dealloc_array(f)

       if (iprnt.ge.2) &
    write(iw,*)'diagram 10 ', dot_product(sigma2,sigma2)

! class : [HBar_DD*R](ab,ij) 
! ==========================

! S(ab,ij) = S(ab,ij) + P(ab) f_bar_ae(b,e) * r2(ae,ij)
! ------------------------------------------------------

 !  call contraction_424 ((/"p1","p3","o1","o2"/),(/"p2","p3"/), &
 ! & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,RightTensor=B%Fbar_ae,LeftTensor=r2,irrep_left=state_sym,irrep_right=1)


    call contraction_244 ((/"p2","p3"/),(/"p1","p3","o1","o2"/), &
   & (/"p2","p1","o1","o2"/),sigma2,1.0d0,-1.0d0,nrep,LeftTensor=B%Fbar_ae,RightTensor=r2,irrep_right=state_sym,irrep_left=1)


       if (iprnt.ge.2) &
    write(iw,*)'diagram 11 ', dot_product(sigma2,sigma2)
     
!!S(ab,ij) = S(ab,ij) - P(ij) f_bar_mi(m,j)*r2(ab,im) 
!!----------------------------------------------------

    call contraction_424 ((/"p1","p2","o1","o3"/),(/"o3","o2"/), &
   & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,RightTensor=B%Fbar_mi,LeftTensor=r2,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.2) &
    write(iw,*)'diagram 12 ',dot_product(sigma2,sigma2)

! S(ab,ij) = S(ab,ij) + w_ijmn(mn,ij) * r2(ab,mn)
!--------------------------------------------------

      call contraction_444((/"p1","p2","o3","o4"/),(/"o3","o4","o1","o2"/), &
     &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=r2,RightTensor=B%W_ijmn,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.2) &
    write(iw,*)'diagram 14 ',dot_product(sigma2,sigma2)

! S(ab,ij) = S(ab,ij) + P(ab)P(ij) w_mbej(mb,ej)*r2(ae,im)
!----------------------------------------------------------

      call contraction_444((/"p1","p3","o1","o3"/),(/"o3","p2","p3","o2"/), &
     &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,0.0d0,nrep,LeftTensor=r2,RightTensor=B%W_mbej,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.2) &
    write(iw,*)'diagram 15 ',dot_product(sigma2,sigma2)

! int_1b(a,f) = int_1b(a,f) - V(nm,fe)*r2(ea,mn)  
! S(ab,ij) = S(ab,ij) + int_1b(a,f)*t2(fb,ij)  
! ----------------------------------------------

    allocate(int_1b(nfvv*rcw))
    int_1b = 0.0d0

     call contraction_442 ((/"p1","p3","o4","o3"/),(/"o4","o3","p4","p3"/),(/"p1","p4"/),int_1b,-1.0d0,1.0d0,nrep,&
    &                     LeftTensor=r2,irrep_left=state_sym,irrep_right=1)

     call contraction_244 ((/"p1","p4"/),(/"p4","p2","o1","o2"/), &
    & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=B%t2,irrep_left=state_sym,irrep_right=1)

    deallocate(int_1b)

       if (iprnt.ge.2) &
    write(iw,*)'diagram 16 ', dot_product(sigma2,sigma2)

! int_1b(n,i) = int_1b(n,i) + V(nm,fe)*r2(fe,im)
! S(ab,ij) = S(ab,ij) - int_1b(n,i)*t(ba,jn)
!!-----------------------------------------------

    allocate(int_1b(nfoo*rcw))
    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"p4","p3","o1","o3"/),(/"o4","o1"/),int_1b,1.0d0,1.0d0,nrep,&
   &                      RightTensor=r2, irrep_left=1,irrep_right=state_sym)

    call contraction_424 ((/"p1","p2","o2","o4"/),(/"o4","o1"/), &
    & (/"p1","p2","o2","o1"/),sigma2,1.0d0,-1.0d0,nrep,RightTensor=int_1b,LeftTensor=B%t2,irrep_left=1,irrep_right=state_sym)

    deallocate(int_1b)

    if (iprnt.ge.2) &
    write(iw,*)'diagram 17 ', dot_product(sigma2,sigma2)

! S(ab,ij) = S(ab,ij) + w_vvvv(ab,ef)*r2(ef,ij)
!------------------------------------------------

      call alloc_array(f,nrep,g)

      call auto_symmetry_offset(f,nvvt,noot,.true.,.true.)
      call auto_symmetry_offset(g,nvvt,noot,.true.,.true.)

      allocate(sigma2_local(N2))
      sigma2_local = 0.0d0 

      do irp = 1, nrep

         mrp = multb(irp+nrep,1+nrep,1)  
         nrp = multb(irp+nrep,state_sym+nrep,1)

         istart = 0
         done =.false.
         if ((nvvt(irp).eq.0).or.(noot(nrp).eq.0)) cycle
         mint = nbuf3/nvvt(irp)
         allocate(w_vvvv(nvvt(irp)*mint*rcw))
 10      kstart = istart
         if (idist(1,1,irp).gt.kstart) kstart=idist(1,1,irp)
         ntot = idist(2,1,irp)
         nint = min0(mint,ntot-kstart)
         call getdz (irp,istart,nint,done,w_vvvv,mint)
         m = nvvt(mrp)
         n = noot(nrp)
         k = nint

         off1 = 1 + f%twoDirac(irp,nrp)*rcw
         off2 = 1 + g%twoDirac(mrp,nrp)*rcw

         call xgemm ('n','n',m,n,k,a1,w_vvvv,m,r2(off1+istart*rcw),m, &
     &               a1,sigma2_local(off2),m)
!    &               a1,sigma2(off2),m)
         if (.not.done) then
            istart = istart + nint
            goto 10
         endif
         deallocate(w_vvvv)
      enddo

#if defined (VAR_MPI)
      if (nmproc .gt. 1) then
        allocate (local_buffer(N2))
        call xcopy (ndimr2,a0,0,local_buffer,1)
        ierr = 0

      if (iprnt.ge.2) &
   write(*,*)'diagram 13a ', ddot(N2,sigma2_local,1,sigma2_local,1),myproc

        call interface_mpi_allreduce(sigma2_local, &
             local_buffer, N2, &
             op_mpi_sum,global_communicator)

      if (iprnt.ge.2) &
   write(*,*)'diagram 13b ', ddot(N2,sigma2_local,1,sigma2_local,1),myproc

        if(ierr.gt.0) then
           call quit('mpi_reduce error in collecting sigma2 !')
        endif
        call xcopy(ndimr2,local_buffer,1,sigma2_local,1)
        deallocate (local_buffer)

     endif
#endif
      if (iprnt.ge.2) &
   write(*,*)'diagram 13c', ddot(N2,sigma2_local,1,sigma2_local,1),myproc

       call XAXPY (ndimr2,A1,sigma2_local,1,sigma2,1)
     deallocate(sigma2_local)

       if (iprnt.ge.2) &
    write(*,*)'diagram 13d', ddot(N2,sigma2,1,sigma2,1),myproc

      call dealloc_array(f,g)

      if (iprnt.ge.2) &
    write(iw,*)'diagram 13 ', dot_product(sigma2,sigma2)

  end subroutine

   subroutine sigma_IP(r1,r2,B,sigma1,sigma2)
   use symmetry_offset
   use modified_sorting

#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "inpt.inc"

!---------------Description--------------------------------------------
!     Calculates 2h1p and 1h type sigma vectors for ionization energy
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------
  real(8), intent(in),contiguous,target   :: r1(:), r2(:) ! R-Amplitudes
  type(intermediates), intent(in)  :: B      
  real(8), intent(inout),contiguous :: sigma1(:),sigma2(:)  ! sigma-vectors
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
  real(8), allocatable :: int_1b(:) 
  real(8)              :: ddot
  integer              :: N1, N2, ndimr1, ndimr2
  type(Offset)         :: e
!--------------------------------------------------
  integer :: off3,off1,off2,mrp,nrp,irp, m, n, k
  type(Offset) :: f,g
  integer  ::row(nrep), column(nrep)
  real*8, allocatable :: some_buf2(:), some_buf3(:)
!--------------------------------------------------

  N1 = size(sigma1,1)
  N2 = size(sigma2,1)

  ndimr1 = N1/rcw
  ndimr2 = N2/rcw

!------------------------------------------
! S(a*,i) = S(a*,i) - Fbar_mi(m,i)*r(a*,m)
!------------------------------------------

  call contraction_222 ((/"c1","o2"/),(/"o2","o1"/),(/"c1","o1"/),r1,B%Fbar_mi,sigma1,-1.0d0,1.0d0,nrep)

  if (iprnt.ge.2) &
  write(iw,*)'diagram 1 ', dot_product(sigma1,sigma1)
!--------------------------------------------
! S(a*,i) = S(a*,i) + Fbar_me(m,e)*r2(a*e,im)
!--------------------------------------------

    call contraction_422 ((/"c1","p2","o1","o2"/),(/"o2","p2"/),(/"c1","o1"/),B%Fbar_me, &
  &                      sigma1,1.0d0,1.0d0,nrep,r2)

  if (iprnt.ge.2) &
  write(iw,*)'diagram 2 ', dot_product(sigma1,sigma1)
!-------------------------------------
!S(a*,i) = S(a*,i) - W_mnie(mn,ie)*r2(a*e,mn)
!-------------------------------------

    call contraction_442((/"c1","p2","o2","o3"/),(/"o2","o3","o1","p2"/),(/"c1","o1"/),sigma1,-1.0d0,1.0d0,nrep, &
   &         RightTensor=B%W_mnie,LeftTensor=r2)

  if (iprnt.ge.2) &
  write(iw,*)'diagram 3 ',dot_product(sigma1,sigma1)

!--------------------------------------------
!  S(a*b,ij) = S(a*b,ij)-W(mb,ij)*r(a*,m)
!--------------------------------------------

   call contraction_244 ((/"c1","o3"/),(/"o3","p2","o1","o2"/), &
  &    (/"c1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,RightTensor=B%W_iemn, &
  &    LeftTensor=r1) 

       if (iprnt.ge.2) &
  write(iw,*)'diagram 4 ',dot_product(sigma2,sigma2)

!----------------------------------------------------
! S(a*b,ij) = S(a*b,ij) - P(ij) f_bar_mi(m,j)*r2(a*b,im) 
!-----------------------------------------------------

    call contraction_424 ((/"c1","p2","o1","o3"/),(/"o3","o2"/), &
   & (/"c1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,RightTensor=B%Fbar_mi,LeftTensor=r2)

!   call contraction_244 ((/"o3","o2"/),(/"c1","p2","o3","o1"/), &
!  & (/"c1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=B%Fbar_mi,RightTensor=r2)

       if (iprnt.ge.2) &
  write(iw,*)'diagram 5 ',dot_product(sigma2,sigma2)

!---------------------------------------------------------
! S(a*b,ij) = S(a*b,ij) + P(ij) w_mbej(mb,ej)*r2(a*e,im)
!----------------------------------------------------------

      call contraction_444((/"c1","p3","o3","o1"/),(/"o3","p2","p3","o2"/), &
     &    (/"c1","p2","o1","o2"/),sigma2,-1.0d0,0.0d0,nrep,LeftTensor=r2,RightTensor=B%W_mbej)

       if (iprnt.ge.2) &
  write(iw,*)'diagram 6 ',dot_product(sigma2,sigma2)

!--------------------------------------------------
! S(a*b,ij) = S(a*b,ij) + w_ijmn(mn,ij) * r2(a*b,mn)
!--------------------------------------------------

      call contraction_444((/"c1","p2","o3","o4"/), (/"o3","o4","o1","o2"/),&
     &    (/"c1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=r2,RightTensor=B%W_ijmn)

       if (iprnt.ge.2) &
  write(iw,*)'diagram 7 ',dot_product(sigma2,sigma2)

!------------------------------------------------------
! S(a*b,ij) = S(a*b,ij) + f_bar_ae(b,e) * r2(a*e,ij)
!------------------------------------------------------

    call contraction_424 ((/"c1","p3","o1","o2"/),(/"p2","p3"/), &
   & (/"c1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,RightTensor=B%Fbar_ae,LeftTensor=r2)

       if (iprnt.ge.2) &
  write(iw,*)'diagram 8 ',dot_product(sigma2,sigma2)

!int_1b(a*,f) = int_1b(a*,f) - V(nm,fe)*r2(a*e,nm)  
!S(a*b,ij) = S(a*b,ij) + int_1b(a*,f)*t2(fb,ij)  
!-----------------------------------------------

    call alloc_array(e,nrep)
    call auto_symmetry_offset(e,nv,ncont,.false.,.false.)
    allocate(int_1b(e%oneNonDirac(1)*rcw))
    int_1b = 0.0d0

    call contraction_442 ((/"c1","p3","o4","o3"/),(/"o4","o3","p4","p3"/),(/"c1","p4"/),int_1b,-1.0d0,1.0d0,nrep,LeftTensor=r2)

    call contraction_244 ((/"c1","p4"/),(/"p4","p2","o1","o2"/), &
    & (/"c1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=B%t2)

    deallocate(int_1b)
   call dealloc_array(e)

 if (iprnt.ge.2) &
  write(iw,*)'diagram 9 ',dot_product(sigma2,sigma2)

 end subroutine

!==============================================================================================================================================================================================================

   subroutine sigma_EA(r1,r2,B,sigma1,sigma2)
   use symmetry_offset
   use modified_sorting

#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "ccpar.inc"
#include "files.inc"
#include "inpt.inc"

!---------------Description--------------------------------------------
!     Calculates 1h2p and 1p type sigma vectors for ionization energy
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------
  real(8), intent(in),contiguous,target   :: r1(:), r2(:) ! R-Amplitudes
  type(intermediates), intent(in)  :: B      
  real(8), intent(inout),contiguous :: sigma1(:),sigma2(:)  ! sigma-vectors
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:),w_vvvv(:)
  real(8), allocatable :: sigma2_local(:)
  real(8), allocatable :: int_1b(:),local_buffer(:) 
  real(8)              :: ddot
  integer :: N1, N2, ndimr1, ndimr2
  integer :: off1,off2,irp,istart,kstart,mint,m,n,k,nint
  integer :: off3,i,j,ntot,ierr
  logical :: done 
  type(Offset) :: f

  integer :: mrp,nrp 
  type(Offset) :: g
  integer  ::row(nrep), column(nrep)
!--------------------------------------------------
  N1 = size(sigma1,1)
  N2 = size(sigma2,1)

  ndimr1 = N1/rcw
  ndimr2 = N2/rcw

  sigma2(1:N2) = 0.0d0

!S(a,i*) = S(a,i*) +  F(a,e)*R(e,i*)   
!---------------------------
   
    call contraction_222 ((/"p1","p2"/),(/"p2","c1"/),(/"p1","c1"/),B%Fbar_ae,r1,sigma1,1.0d0,1.0d0,nrep)

    if (iprnt.ge.2) &
    write(iw,*)'diagram 1 ', ddot(N1,sigma1,1,sigma1,1)

!S(a,i*) = S(a,i*) + F_bar_me(m,e)*r2(ae,i*m)
!--------------------------------------

    call contraction_422 ((/"p1","p2","c1","o2"/),(/"o2","p2"/), &
   &   (/"p1","c1"/),B%Fbar_me,sigma1,1.0d0,1.0d0,nrep,LeftTensor=r2)

!   call contraction_242 ((/"o2","p2"/),(/"p1","p2","c1","o2"/), &
!  &   (/"p1","c1"/),B%Fbar_me,sigma1,1.0d0,1.0d0,nrep,RightTensor=r2)

   if (iprnt.ge.2) &
   write(iw,*)'diagram 2 ', ddot(N1,sigma1,1,sigma1,1)

!S(a,i*) = S(a,i*) + w_amef(am,ef) * r2(ef,i*m)
!----------------------------------------
 
    call contraction_442((/"p1","o2","p2","p3"/),(/"p2","p3","c1","o2"/),(/"p1","c1"/),sigma1,1.0d0,1.0d0,nrep, &
   &         LeftTensor=B%W_amef,RightTensor=r2)

   if (iprnt.ge.2) &
   write(iw,*)'diagram 3 ', ddot(N1,sigma1,1,sigma1,1)

!S(ab,i*j) = S(ab,i*j) +  w_efam(ab,ej)*r1(e,i*)
!--------------------------------------------------

    call contraction_424 ((/"p1","p2","p3","o2"/),(/"p3","c1"/), &
   &    (/"p1","p2","c1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=B%W_efam, &
   &    RightTensor=r1) 

   if (iprnt.ge.2) &
    write(iw,*)'diagram 4 ', ddot(N2,sigma2,1,sigma2,1)

! S(ab,i*j) = S(ab,i*j) + P(ab) f_bar_ae(b,e) * r2(ae,i*j)
!-------------------------------------------------------

    call contraction_424 ((/"p1","p3","c1","o2"/),(/"p2","p3"/), &
   & (/"p1","p2","c1","o2"/),sigma2,1.0d0,1.0d0,nrep,RightTensor=B%Fbar_ae,LeftTensor=r2)

   if (iprnt.ge.2) &
    write(iw,*)'diagram 5 ', ddot(N2,sigma2,1,sigma2,1)

! S(ab,i*j) = S(ab,i*j) - P(ij) f_bar_mi(m,j)*r2(ab,i*m) 
!-----------------------------------------------------

    call contraction_424 ((/"p1","p2","c1","o3"/),(/"o3","o2"/), &
   & (/"p1","p2","c1","o2"/),sigma2,-1.0d0,1.0d0,nrep,RightTensor=B%Fbar_mi,LeftTensor=r2)

       if (iprnt.ge.2) &
    write(iw,*)'diagram 6 ', ddot(N2,sigma2,1,sigma2,1)

! S(ab,i*j) = S(ab,i*j) + P(ab) w_mbej(mb,ej)*r2(ae,i*m)
!----------------------------------------------------------

      call contraction_444((/"p1","p3","c1","o3"/),(/"o3","p2","p3","o2"/), &
     &    (/"p1","p2","c1","o2"/),sigma2,1.0d0,0.0d0,nrep,LeftTensor=r2,RightTensor=B%W_mbej)

!     call contraction_444((/"o3","p2","p3","o2"/),(/"p1","p3","c1","o3"/), &
!    &    (/"p1","p2","c1","o2"/),sigma2,1.0d0,0.0d0,nrep,RightTensor=r2,LeftTensor=B%W_mbej)

!     call contraction_444((/"p2","o3","o2","p3"/),(/"p1","p3","c1","o3"/), &
!    &    (/"p1","p2","c1","o2"/),sigma2,1.0d0,0.0d0,nrep,RightTensor=r2,LeftTensor=B%W_ejmb)

   if (iprnt.ge.2) &
    write(iw,*)'diagram 8 ', ddot(N2,sigma2,1,sigma2,1)

!int_1b(n,i*) = int_1b(n,i*) + V(nm,fe)*r2(fe,i*m)
!S(ab,i*j) = S(ab,i*j) - int_1b(n,i*)*t(ba,jn)
!------------------------------------------------

    call alloc_array(f,nrep)

    call auto_symmetry_offset(f,no,ncont,.false.,.false.)
    allocate(int_1b(f%oneNonDirac(1)*rcw))
    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"p4","p3","c1","o3"/),(/"o4","c1"/),int_1b,1.0d0,1.0d0,nrep,RightTensor=r2)

    call contraction_424 ((/"p2","p1","o2","o4"/),(/"o4","c1"/), &
    & (/"p2","p1","c1","o2"/),sigma2,-1.0d0,-1.0d0,nrep,RightTensor=int_1b,LeftTensor=B%t2)

    deallocate(int_1b)

       if (iprnt.ge.2) &
    write(iw,*)'diagram 9 ', ddot(N2,sigma2,1,sigma2,1)

! S(ab,i*j) = S(ab,i*j) + w_vvvv(ab,ef)*r2(ef,i*j)
!------------------------------------------------

      allocate(sigma2_local(N2))
      sigma2_local = 0.0d0

      off1 = 1
      off2 = 1
      do irp = 1, nrep
         istart = 0
         done =.false.
         if ((nvvt(irp).eq.0).or.(f%oneDirac(irp).eq.0)) cycle
         mint = nbuf3/nvvt(irp)
 10      kstart = istart
         if (idist(1,1,irp).gt.kstart) kstart=idist(1,1,irp)
         ntot = idist(2,1,irp)
         nint = min0(mint,ntot-kstart)
         allocate(w_vvvv(nvvt(irp)*nint*rcw))
         call getdz (irp,istart,nint,done,w_vvvv,mint)
         m = nvvt(irp)
         n = f%oneDirac(irp)
         k = nint
         call xgemm ('n','n',m,n,k,a1,w_vvvv,m,r2(off1+istart*rcw),m, &
     &               a1,sigma2_local(off2),m)
!    &               a1,sigma2(off2),m)

         deallocate(w_vvvv)

         if (.not.done) then
            istart = istart + nint
            goto 10
         endif
         off1 = off1 + n * m * rcw
         off2 = off2 + m * n * rcw
      enddo

      call dealloc_array(f)
      allocate (local_buffer(N2))
#if defined (VAR_MPI)
      if (nmproc .gt. 1) then
        call xcopy (ndimr2,a0,0,local_buffer,1)
        ierr = 0

      if (iprnt.ge.2) &
   write(*,*)'diagram 7a ', ddot(N2,sigma2_local,1,sigma2_local,1),myproc
!  write(*,*)'diagram 7a ', ddot(N2,sigma2,1,sigma2,1),myproc

        call interface_mpi_allreduce(sigma2_local, &
!       call interface_mpi_allreduce(sigma2, &
             local_buffer, N2, &
             op_mpi_sum,global_communicator)
      if (iprnt.ge.2) &
   write(*,*)'diagram 7b ', ddot(N2,sigma2_local,1,sigma2_local,1),myproc
!  write(*,*)'diagram 7b ', ddot(N2,sigma2,1,sigma2,1),myproc

        if(ierr.gt.0) then
           call quit('mpi_reduce error in collecting sigma2 !')
        endif
        call xcopy(ndimr2,local_buffer,1,sigma2_local,1)
!       call xcopy(ndimr2,local_buffer,1,sigma2,1)
      endif
#endif
      deallocate (local_buffer)

      if (iprnt.ge.2) &
   write(*,*)'diagram 7c ', ddot(N2,sigma2_local,1,sigma2_local,1),myproc

       call XAXPY (ndimr2,A1,sigma2_local,1,sigma2,1)
     deallocate(sigma2_local)

       if (iprnt.ge.2) &
    write(*,*)'diagram 7d ', ddot(N2,sigma2,1,sigma2,1),myproc

 end subroutine

   subroutine Lambda_IP(l1,l2,B,lambda1,lambda2)

  use symmetry_offset
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "inpt.inc"

!---------------Description--------------------------------------------
!     Calculates 2h1p and 1h type sigma vectors for ionization energy
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------
  real(8), intent(in),contiguous,target   :: L1(:), L2(:) ! R-Amplitudes
  type(intermediates), intent(in)  :: B      
  real(8), intent(inout),contiguous :: lambda1(:),lambda2(:)  ! sigma-vectors
!---------------Local Variables--------------------
  real(8), allocatable :: int_1b(:), G(:) 
  real(8)              :: ddot
  integer              :: N1, N2, ndiml1, ndiml2
  type(Offset)         :: e
!--------------------------------------------------

  N1 = size(lambda1,1)
  N2 = size(lambda2,1)

  ndiml1 = N1/rcw
  ndiml2 = N2/rcw

!------------------------------------------
! L(i,a*) = L(i,a*) - L(m,a*)*f_bar_mi(i,m)  
!------------------------------------------

     call contraction_222 ((/"o2","c1"/),(/"o1","o2"/),(/"o1","c1"/), &
  &    l1,B%Fbar_mi,lambda1,-1.0d0,1.0d0,nrep)

       if (iprnt.ge.0) &
    write(iw,*)'diagram 1 ', ddot(N1,lambda1,1,lambda1,1)

!-----------------------------------------------------
! L(i,a*) = L(i,a*) -  L2(mn,a*e)*W(ie,mn)
!-----------------------------------------------------

!   call contraction_442((/"o2","o3","c1","p2"/),(/"o1","p2","o2","o3" &
! &          /),(/"c1","p1"/),lambda1,-1.0d0,1.0d0,nrep,LeftTensor=l2,  &
! &          RightTensor=B%W_iemn)


    call contraction_442((/"o1","p2","o2","o3"/), (/"o2","o3","c1","p2"/), &
  &          (/"o1","c1"/),lambda1,-1.0d0,1.0d0,nrep,RightTensor=l2,  &
  &          LeftTensor=B%W_iemn)

       if (iprnt.ge.0) &
    write(iw,*)'diagram 2 ', ddot(N1,lambda1,1,lambda1,1)

!-------------------------------------------------
!  L(ij,a*b) = L(ij,a*b) + L2 (ij,a*e)*f_bar_ae(e,b)
!-------------------------------------------------

     call contraction_424 ((/"o1","o2","c1","p3"/),(/"p3","p2"/), &
   &    (/"o1","o2","c1","p2"/),lambda2,1.0d0,1.0d0,nrep,LeftTensor=l2, &
   &    RightTensor=B%Fbar_ae)

       if (iprnt.ge.0) &
    write(iw,*)'diagram 3 ', ddot(N2,lambda2,1,lambda2,1)

!----------------------------------------------
!  L(ij,a*b) = L(ij,a*b) - P(ij) L2 (im,a*b)*f_bar_im (j,m)
!----------------------------------------------

       call contraction_424 ((/"o1","o3","c1","p2"/),(/"o2","o3"/), &
     &    (/"o1","o2","c1","p2"/),lambda2,-1.0d0,1.0d0,nrep, &
     &    LeftTensor=l2,RightTensor=B%Fbar_mi)

       if (iprnt.ge.0) &
    write(iw,*)'diagram 4 ', ddot(N2,lambda2,1,lambda2,1)

!--------------------------------------------------------
!  L(ij,a*b) = L(ij,a*b) + L2(mn,a*b) * W_ijmn(ij,mn)
!--------------------------------------------------------

    call contraction_444((/"o1","o2","o3","o4"/),(/"o3","o4","c1","p2"/), &
   &    (/"o1","o2","c1","p2"/),lambda2,1.0d0,1.0d0,nrep,RightTensor=l2,LeftTensor=B%W_ijmn)


       if (iprnt.ge.0) &
    write(iw,*)'diagram 5 ', ddot(N2,lambda2,1,lambda2,1)
!--------------------------------------------------------
!  L(ij,a*b) = L(ij,a*b) + P(ij) L2(im,a*e) * W_bmej(je,bm)
!--------------------------------------------------------

     call contraction_444((/"o1","o3","c1","p3"/),(/"o2","p3","p2","o3"/), (/"o1","o2","c1","p2"/),lambda2, &
    &    1.0d0,0.0d0,nrep,righttensor=B%W_mbej,LeftTensor=l2)


       if (iprnt.ge.0) &
    write(iw,*)'diagram 6 ', ddot(N2,lambda2,1,lambda2,1)
!---------------------------------------------------------------------------
!  L(ij,a*b) = L(ij,a*b) + L1(m,a*) * w_mnie(ij,mb)
!---------------------------------------------------------------------------

     call contraction_424 ((/"o1","o2","o3","p2"/),(/"o3","c1"/), &
    &   (/"o1","o2","c1","p2"/),lambda2,-1.0d0,1.0d0,nrep,RightTensor=l1, &
    &   LeftTensor=B%W_mnie)


       if (iprnt.ge.0) &
    write(iw,*)'diagram 7 ', ddot(N2,lambda2,1,lambda2,1)
!---------------------------------------------------------------
!  G (e,a*)  = G(e,a*) + L2(mn,a*f) * t2(ef,mn) 
!  L(ij,a*b) = L(ij,a*b) + V(ij,eb)*G (e,a*)  
!---------------------------------------------------------------

    call alloc_array(e,nrep)
    call auto_symmetry_offset(e,nv,ncont,.false.,.false.)

    allocate (G(e%oneNonDirac(1)*rcw))

     G = 0.0d0

     call contraction_442((/"p3","p4","o3","o4"/),(/"o3","o4","c1","p4"/),(/"p3","c1"/),G,-1.0d0,1.0d0,nrep,LeftTensor=B%t2,  &
    &          RightTensor=l2)

     call contraction_424 ((/"o1","o2","p3","p2"/),(/"p3","c1"/), &
    &        (/"o1","o2","c1","p2"/),Lambda2,1.0d0,1.0d0,nrep,RightTensor=G)

    deallocate(G)

       if (iprnt.ge.0) &
    write(iw,*)'diagram 8 ', ddot(N2,lambda2,1,lambda2,1)

     call dealloc_array(e)

   call lambda2_disconnected_IP(l1,B%Fbar_me,lambda2)

       if (iprnt.ge.0) &
    write(iw,*)'diagram 9 ', ddot(N2,lambda2,1,lambda2,1)

 end subroutine

   subroutine Lambda_EE(l1,l2,B,lambda1,lambda2,state_sym)
  use intermediates_1b_2b 
  use symmetry_offset
  use lambda_equation  
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "inpt.inc"
#include "ccpar.inc"

!---------------Description--------------------------------------------
!     Calculates 2h2p and 1h1p type sigma vectors for Excitation energy
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------
  real(8), intent(in),contiguous,target   :: L1(:), L2(:) ! R-Amplitudes
  type(intermediates), intent(in)  :: B      
  integer, intent(in) :: state_sym
  real(8), intent(inout),contiguous :: lambda1(:),lambda2(:)  ! sigma-vectors
!---------------Local Variables--------------------
  real(8), allocatable :: int_oo(:), G_vv(:), G_oo(:) 
  real(8), allocatable :: w_vvvv(:),tau(:) 
  real(8)              :: ddot
  logical              :: done
  integer     :: jklrep,m2,irep,jrep,abrep,jkloff,ij,ijkl,ijkl1,mint
  integer     :: off2,i,j,istart,m,k,n,irp,ntot,nint,kstart,off1
  integer              :: N1, N2, ndiml1, ndiml2
  type(Offset)         :: e
!--------------------------------------------------

  N1 = size(lambda1,1)
  N2 = size(lambda2,1)

  ndiml1 = N1/rcw
  ndiml2 = N2/rcw

!--------------------------------------------------
! l1(1:ndiml1) => l1
! l2(1:ndiml2) => l2

     allocate(G_vv(nfvv*rcw))
     allocate(G_oo(nfoo*rcw))

      G_vv = 0.0d0
      G_oo = 0.0d0


  call contraction_442 ((/"p2","p3","o1","o2"/),(/"o1","o2","p1","p3"/),(/"p2","p1"/),G_vv,-1.0d0,1.0d0,nrep,&
                        &  LeftTensor=B%t2,RightTensor=l2,irrep_left=1,irrep_right=state_sym)

  call contraction_442 ((/"o1","o3","p1","p2"/),(/"p1","p2","o2","o3"/),(/"o1","o2"/),G_oo,1.0d0,1.0d0,nrep,&
                        & LeftTensor=l2,RightTensor=B%t2,irrep_left=state_sym,irrep_right=1)

       call getoovv (lambda2)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda2 after initialization',dot_product(lambda2,lambda2)
!-------------------------------------------------
!   term1 : lambda2 = lambda2 + P(ab) L2 (ij,ae)*f_bar_ae(e,b)
!-------------------------------------------------

       call contraction_424 ((/"o1","o2","p1","p3"/),(/"p3","p2"/), &
     &    (/"o1","o2","p1","p2"/),lambda2,1.0d0,1.0d0,nrep,LeftTensor=L2, &
     &    RightTensor=B%Fbar_ae,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.1) &

     & write(iw,*)'lambda2 after term1',dot_product(lambda2,lambda2)
!----------------------------------------------
!   term2 : lambda2 = lambda2 - P(ij) L2 (im,ab)*f_bar_im (j,m)
!----------------------------------------------

       call contraction_424 ((/"o1","o3","p1","p2"/),(/"o2","o3"/), &
     &    (/"o1","o2","p1","p2"/),lambda2,-1.0d0,1.0d0,nrep, &
     &    LeftTensor=L2,RightTensor=B%Fbar_mi,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.1) &

     & write(iw,*)'lambda2 after term2',dot_product(lambda2,lambda2)
!--------------------------------------------------------
!   term3 : lambda2 = lambda2 + L2(mn,ab) * W_ijmn(ij,mn)
!--------------------------------------------------------

      call contraction_444((/"o1","o2","o3","o4"/),(/"o3","o4","p1","p2"/), &
     &    (/"o1","o2","p1","p2"/),lambda2, &
     &    1.0d0,1.0d0,nrep,RightTensor=L2,LeftTensor=B%W_ijmn,irrep_right=state_sym,irrep_left=1)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda2 after term3', dot_product(lambda2,lambda2)

!--------------------------------------------------------
!   term4 : lambda2 = lambda2 + P(ij)P(ab) L2(im,ae) * W_bmej(je,bm)
!--------------------------------------------------------

      call contraction_444((/"o1","o3","p1","p3"/),(/"o2","p3","p2","o3"/), (/"o1","o2","p1","p2"/),lambda2, &
     &    1.0d0,0.0d0,nrep,righttensor=B%W_mbej,LeftTensor=L2,irrep_left=state_sym,irrep_right=1)

!-----------------------------------------------------------
!   term5 : lambda2 = lambda2 + P(ab) V(ij,ae)*G_vv(e,b)
!-----------------------------------------------------------

      call contraction_424 ((/"o1","o2","p1","p3"/),(/"p3","p2"/), &
     &        (/"o1","o2","p1","p2"/),lambda2,1.0d0,1.0d0,nrep,RightTensor=G_vv,irrep_left=1,irrep_right=state_sym)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda2 after term5',dot_product(lambda2,lambda2)

!---------------------------------------------------------------------------
!     term6: L2(ij,ab) = L2(ij,ab) + L1(m,a) * w_mnie(ij,mb)
!---------------------------------------------------------------------------

      call contraction_424 ((/"o1","o2","o3","p2"/),(/"o3","p1"/), &
     &   (/"o1","o2","p2","p1"/),lambda2,-1.0d0,-1.0d0,nrep,RightTensor=L1, &
     &   LeftTensor=B%W_mnie,irrep_left=1,irrep_right=state_sym)


       if (iprnt.ge.1) &
     & write(iw,*)'lambda2 after term6',dot_product(lambda2,lambda2)
!--------------------------------------------------------------------
!    term7 : lambda2 = lambda2 - P(ij) V(im,ab)*(G_oo(j,m)+L1(j,e)*t(e,m))
!--------------------------------------------------------------------

      call contraction_424 ((/"o1","o3","p1","p2"/),(/"o2","o3"/), &
     &        (/"o1","o2","p1","p2"/),lambda2,-1.0d0, &
     &        1.0d0,nrep,RightTensor=G_oo,irrep_left=1,irrep_right=state_sym)

      allocate(int_oo(IOO(NREP+1)*rcw))

      int_oo = 0.0d0

      call contraction_222((/"o1","p3"/),(/"p3","o3"/),(/"o1","o3"/),L1, &
     &        B%t1,int_oo,1.0d0,1.0d0,nrep,irrep_left=state_sym,irrep_right=1)

      call contraction_424 ((/"o3","o2","p1","p2"/),(/"o1","o3"/), &
     &        (/"o1","o2","p1","p2"/),lambda2,-1.0d0, &
     &        1.0d0,nrep,RightTensor=int_oo,irrep_left=1,irrep_right=state_sym)

      deallocate(int_oo)

      if (iprnt.ge.1) &

      & write(iw,*)'lambda2 after term7',dot_product(lambda2,lambda2)

!-------------------------------------------------------------------------
!   term 8 : lambda2 = lambda2 + expression to be typed
!-------------------------------------------------------------------------

     call contraction_244 ((/"o1","p3"/),(/"p3","o2","p1","p2"/), &
    &  (/"o1","o2","p1","p2"/),lambda2,1.0d0,1.0d0,nrep,LeftTensor=L1,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda2 after term8',dot_product(lambda2,lambda2)
!-------------------------------------------------------------------------
!   term 9 : lambda2 = lambda2 + P(ij)P(ab) L1(i,a)*f_bar_me(j,b)
!-------------------------------------------------------------------------

!      call lambda2_disconnected(l1,B%Fbar_me,lambda2,state_sym)
      call lambda2_disconnected(l1,B%Fbar_me,lambda2)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda2 after term9',dot_product(lambda2,lambda2)

!--------------------------------------------------------
!    term 10 : lambda2 = lambda2 + L2(ij,ef) * W_vvvv(ef,ab)
!--------------------------------------------------------

      allocate(tau(N2))
      call gettau(B%t1,B%t2,tau)
      off2 = 1
      off1 = 1
      do irp = 1, nrep
         istart = 0
         if (nvvt(irp).eq.0) cycle
         mint = nbuf3/nvvt(irp)
 10      kstart = istart
         if (idist(1,1,irp).gt.kstart) kstart=idist(1,1,irp)
         ntot = idist(2,1,irp)
         nint = min0(mint,ntot-kstart)

         allocate(w_vvvv(nvvt(irp)*nint*rcw))

         call getdz (irp,istart,nint,done,w_vvvv,nbuf3)

         n = nvvt(irp)
         m = noot(irp)
         k = nint
         call xgemm ('n','n',m,k,n,A1,l2(off1),m, &
     &           w_vvvv,n,a1,lambda2(off2+istart*m*rcw),m)

         deallocate(w_vvvv)
         if (.not.done) then
            istart = istart + nint
            goto 10
         endif
         off2 = off2 + m * n * rcw
         off1 = off1 + m * n * rcw
       enddo

       deallocate(tau)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda2 after term10',dot_product(lambda2,lambda2)

      call xcopy (ndiml1,B%Fbar_me,1,lambda1,1)

     if (iprnt.ge.1) &
     &   write(iw,*)'lambda1: init',dot_product(lambda1,lambda1)
!---------------------------------
!    lambda1 = lambda1 + L(i,e) * f_bar_ae(e,a)
!---------------------------------

       call contraction_222 ((/"o1","p2"/),(/"p2","p1"/),(/"o1","p1"/), &
     &    L1,B%Fbar_ae,lambda1,1.0d0,1.0d0,nrep,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda1: after term1',dot_product(lambda1,lambda1)
!----------------------------------------
!    lambda1 = lambda1 - L(m,a)*f_bar_mi(i,m)
!----------------------------------------

        call contraction_222 ((/"o2","p1"/),(/"o1","o2"/),(/"o1","p1"/), &
     &    L1,B%Fbar_mi,lambda1,-1.0d0,1.0d0,nrep,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda1: after term2',dot_product(lambda1,lambda1)
!----------------------------------------------
!      lambda1 = lambda1 - G(m,n) * W(mi,na)
!      w_mina(mi,na) = w_mina(mi,na) + V(mi,fa)*t(f,n)
!----------------------------------------------

      call contraction_242 ((/"o3","o2"/),(/"o2","o1","o3","p1"/), &
     & (/"o1","p1"/),G_oo,lambda1,-1.0d0,1.0d0,nrep,RightTensor=B%W_mnie,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda1: after term3',dot_product(lambda1,lambda1)
!-------------------------------------------------
!      lambda1 = lambda1 + L1(m,e) * W(ie,am)
!-------------------------------------------------

      call contraction_242 ((/"o2","p2"/),(/"o1","p2","p1","o2"/), &
     &   (/"o1","p1"/),L1,lambda1,1.0d0,1.0d0,nrep,RightTensor=B%W_mbej,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda1 after term4',dot_product(lambda1,lambda1)
!-----------------------------------------------------
!       lambda1 = lambda1 -  L2(mn,ae)*W(ie,mn)
!-----------------------------------------------------

      call contraction_442((/"o2","o3","p1","p2"/),(/"o1","p2","o2","o3" &
     &          /),(/"o1","p1"/),lambda1,-1.0d0,1.0d0,nrep,LeftTensor=L2,  &
     &          RightTensor=B%W_iemn,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda1 after term5',dot_product(lambda1,lambda1)

!-------------------------------------------------------
!      lambda1    = lambda1 + L2(im,ef) * W(ef,am)
!-------------------------------------------------------

      call contraction_442((/"o1","o2","p2","p3"/),(/"p2","p3","p1","o2"/),(/"o1","p1"/),lambda1,1.0d0,1.0d0,nrep, &
       &         LeftTensor=L2,RightTensor=B%W_efam,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda1 after term6',dot_product(lambda1,lambda1)

!-------------------------------------------------------
!    lambda1 = lambda1 - G(e,f) * V(ei,fa)
!-------------------------------------------------------

      call contraction_242 ((/"p3","p2"/),(/"p2","o1","p3","p1"/), &
     &  (/"o1","p1"/),G_vv,lambda1,-1.0d0,1.0d0,nrep,RightTensor=B%W_amef,irrep_left=state_sym,irrep_right=1)

       if (iprnt.ge.1) &
     & write(iw,*)'lambda1 after term7',dot_product(lambda1,lambda1)

      deallocate(G_oo)
      deallocate(G_vv)

      end subroutine


  subroutine lambda2_disconnected_IP(l1,f_bar,lambda2)

!---------------Description--------------------------------------------

!Calculating the contribution from the disconnected piece in lambda2 of IP. 

  use symmetry_offset
  implicit none
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------

      real*8, intent(in)    :: l1(*),f_bar(:)
      real*8, intent(inout) :: lambda2(*)

!---------------Local variables--------------------------------------
      integer :: ijaboff,abrep,arep,brep,a,b,amin,imin,i,j,ij,jb,ia,ja, &
     &           ib,ni,irep,jrep
      real*8  :: L1min(2)
      type(Offset)         :: e
!---------------Executable code--------------------------------------

    call alloc_array(e,nrep)
    call auto_symmetry_offset(e,no,ncont,.false.,.false.)

  L1min=(/0.0d0,0.0d0/)

  ijaboff = 0
  do abrep = 1, nrep
  do brep = 1, nrep
  arep = multb(brep,abrep+nrep,2)
  do b = 1, nv(brep)
     do a = 1, ncont(arep)
        do jrep = 1, nrep
        irep = multb(jrep,abrep+nrep,2)
        if (irep.lt.jrep) cycle
!--------------------------------------------------------------
! L2(ij,a*b) = L2(ij,a*b) + L1(i,a*) * f_bar(j,b) 
!--------------------------------------------------------------
         if (jrep.eq.brep) then
           ij = 1
           jb = (ivo(brep) + (b-1) * no(brep)) * rcw + 1
           do j = 1, no(jrep)
              imin = 1
              if (irep.eq.jrep) imin = j + 1
              ia = (e%twoNonDirac(irep,arep) + (a-1) * no(arep) + imin - 1)*rcw + 1 
              ni = no(irep) - imin + 1
              call xaxpy (ni,f_bar(jb),l1(ia),1,lambda2(ijaboff+ij),1)
              ij = ij + ni * rcw
              jb = jb + rcw
           enddo
         endif
!--------------------------------------------------------------
! L2(ij,a*b) = L2(ij,ab) - L1(j,a*) * f_bar(i,b) 
!--------------------------------------------------------------
!        if (jrep.eq.arep) then
        if ((irep.eq.brep).and.(jrep.eq.arep)) then
           ij = 1
           ja = (e%twoNonDirac(jrep,arep) + (a-1) * no(jrep)) * rcw + 1  
           do j = 1, no(jrep)
              imin = 1
              if (irep.eq.jrep) imin = j + 1
              ib = (iivo(irep,brep) + (b-1) * no(brep) + imin - 1)*rcw + 1
!              ib = (ivo(brep) + (b-1) * no(brep) + imin - 1)*rcw + 1
              ni = no(irep) - imin + 1
               L1min(1) = -l1(ja)
               if (carith) L1min(2) = -l1(ja+1)

!               write(*,*)'print when it goes out of bounds', ib

              call xaxpy (ni,L1min,f_bar(ib),1,lambda2(ijaboff+ij),1)
              ij = ij + ni * rcw
              ja = ja + rcw
           enddo
         endif
!--------------------------------------------------------
! update offset and go to next irrep pair
!--------------------------------------------------------
       if (irep.ne.jrep) then
          ijaboff = ijaboff + no(irep) * no(jrep) * rcw
       else
          ijaboff = ijaboff + no(irep) * (no(irep)-1) * rcw / 2
       endif
       enddo
    enddo
 enddo
 enddo
 enddo

  call dealloc_array(e)

 end subroutine

      subroutine show_number_of_threads(write_unit) 
#ifdef HAVE_OPENMP
         use omp_lib, only : omp_get_num_threads, omp_get_num_procs
#endif
         integer :: number_of_threads
         integer :: number_of_procs
         integer :: write_unit

#ifdef HAVE_OPENMP
         number_of_threads = omp_get_num_threads()
         number_of_procs   = omp_get_num_procs()
#else
         number_of_threads = 1
         number_of_procs   = 1
#endif
         write(write_unit,"(2x,a,2i6)") "Number of OMP threads, procs in use:",number_of_threads,number_of_procs
      end subroutine

   end module 
