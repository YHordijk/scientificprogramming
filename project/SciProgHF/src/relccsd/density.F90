  module eom_density
      use memory_allocator
      use allocator_parameters, only : klongint, kreal

  use contraction
  implicit none

  public density_IP
  public density_EE

    private

    character(8) :: eomtype
    integer :: nroot,state_sym

   contains

    subroutine density_IP (l1,l2,r1,r2,t1,t2,dmo)

  use symmetry_offset
#include "symm.inc"
#include "complex.inc"
#include "inpt.inc"
!----------------calling variables-------------------------
   real*8, intent(in),contiguous,target :: l1(:),l2(:),r1(:),r2(:)
   real*8, intent(in) :: t1(:),t2(:)
   real*8, intent(inout) :: dmo(*)
!----------------local variables---------------------------
   real*8, allocatable, target :: doo(:),dvv(:),dvo(:),dov(:)  
   real*8, allocatable :: work1(:),work2(:) 
   real*8, allocatable :: woo_temp(:),wpp_temp(:),wop_temp(:)
   real(8), pointer,contiguous :: l1_ptr(:),l2_ptr(:),r1_ptr(:),r2_ptr(:)
   type(Offset)         :: e
   integer              :: ndiml1,ndiml2
!-----------------------------------------------------------

  ndiml1 = size(l1,1) 
  ndiml2 = size(l2,1) 

  l1_ptr(1:ndiml1) => l1
  l2_ptr(1:ndiml2) => l2
  r1_ptr(1:ndiml1) => r1
  r2_ptr(1:ndiml2) => r2

!   allocate(doo(nfoo*rcw))
!   allocate(dvv(nfvv*rcw))
!   allocate(dvo(nfvo*rcw))
!   allocate(dov(nfvo*rcw))

    call alloc(doo, (nfoo*rcw),  id="doo")
    call alloc(dvv, (nfvv*rcw),  id="dvv")
    call alloc(dvo, (nfvo*rcw),  id="dvo")
    call alloc(dov, (nfvo*rcw),  id="dov")

     doo = 0.0d0
     dvv = 0.0d0
     dvo = 0.0d0
     dov = 0.0d0

!-------------------------------------------------------------------
!   doo(j,i) = doo(j,i) - l2(jm,e*f) * r2(e*f,im)
!-------------------------------------------------------------------

  call contraction_442 ((/"o2","o3","c1","p2"/),(/"c1","p2","o1","o3"/),(/"o2","o1"/),doo,1.0d0,1.0d0,nrep, &
 &  LeftTensor=l2,RightTensor=r2)

    if (iprnt.ge.2) write(*,*)'diagram1',dot_product(doo,doo)

!---------------------------------------------------------------------
!    doo(j,i) = doo(j,i) - l1(j,e*)*r1(e*,i)
!---------------------------------------------------------------------

  call contraction_222 ((/"o2","c1"/),(/"c1","o1"/),(/"o2","o1"/),l1,r1,doo,1.0d0,1.0d0,nrep)

    if (iprnt.ge.2) write(*,*)'diagram2',dot_product(doo,doo) 
!---------------------------------------------------------------------
!    doo(j,i) = doo(j,i) - l2(im,a*e)*r1(a*,m)*t1(e,j) 

!    wop_temp(i,e) = wop_temp(i,e) + l2(im,a*e)*r1(a*,m) 
!    doo(j,i) = doo(j,i) + wop_temp(i,e) * t1(e,j)
!---------------------------------------------------------------------

  call alloc(wop_temp, (nfvo*rcw), id="wop_temp vo") 

  wop_temp = 0.0d0

  call contraction_422 ((/"o1","o3","c1","p2"/),(/"c1","o3"/),(/"o1","p2"/),r1,wop_temp,-1.0d0,1.0d0,nrep,LeftTensor=l2)

  call contraction_222 ((/"p2","o2"/),(/"o1","p2"/),(/"o2","o1"/),t1,wop_temp,doo,1.0d0,1.0d0,nrep)

! call dealloc(wop_temp, id="wop_temp vo")

   if (iprnt.ge.2) write(*,*)'diagram3',dot_product(doo,doo) 
!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l2(mn,ea)*r2(eb,mn) 
!----------------------------------------------------------------------

  call contraction_442 ((/"c1","p2","o1","o2"/),(/"o1","o2","c1","p1"/),(/"p2","p1"/),dvv,1.0d0,1.0d0,nrep, &
 &     RightTensor=l2,LeftTensor=r2)

   if (iprnt.ge.2) write(*,*)'diagram4',dot_product(dvv,dvv)
!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l2(mn,cb)*r1(c,n)*t1(a,m)
!     wop_temp(m,b) = wop_temp(m,b) - l2(mn,cb)*r1(c,n) 
!     dvv(b,a) = dvv(b,a) + wop_temp(m,b)*t1(a,m)
!----------------------------------------------------------------------

! call alloc(wop_temp, (nfvo*rcw), id="wop_temp vo")
  wop_temp = 0.0d0

  call contraction_422 ((/"o1","o2","c1","p2"/),(/"c1","o2"/),(/"o1","p2"/),r1,wop_temp,-1.0d0,1.0d0,nrep,LeftTensor=l2)

  call contraction_222 ((/"o1","p2"/),(/"p1","o1"/),(/"p2","p1"/),wop_temp,t1,dvv,-1.0d0,1.0d0,nrep)

  call dealloc(wop_temp, id="wop_temp vo")

   if (iprnt.ge.2) write(*,*)'diagram5',dot_product(dvv,dvv)
!----------------------------------------------------------------------------
!     dvo(a,i) = t1(a,i)
!----------------------------------------------------------------------------

      call xcopy (nfvo,t1,1,dvo,1)

   if (iprnt.ge.2) write(*,*)'diagram6',dot_product(dvo,dvo)
!--------------------------------------------------------------------------
!     dvo(a,i) = dvo(a,i) + l1(m,b) * r(ba,im)
!--------------------------------------------------------------------------

      call contraction_242 ((/"o2","c1"/),(/"c1","p1","o2","o1"/), & 
     & (/"p1","o1"/),l1,dvo,-1.0d0,1.0d0,nrep,RightTensor=r2)

!     call contraction_422 ((/"c1","p1","o1","o2"/),(/"o2","c1"/), &
!    & (/"p1","o1"/),l1_ptr,dvo,-1.0d0,1.0d0,nrep,LeftTensor=r2_ptr)       

   if (iprnt.ge.2) write(*,*)'diagram7',dot_product(dvo,dvo)
!----------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l1(m,e) * r(e,i)*t(a,m)
!----------------------------------------------------------

      call alloc(woo_temp, (nfoo*rcw), id="wop_temp oo")

      woo_temp = 0.0d0

      call contraction_222 ((/"o2","c1"/),(/"c1","o1"/),(/"o2","o1"/), &
     &              l1_ptr,r1_ptr,woo_temp,1.0d0,1.0d0,nrep)                   


      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/), &
     &              woo_temp,t1,dvo,-1.0d0,1.0d0,nrep)                 

!     call dealloc(woo_temp, id="wop_temp oo")

   if (iprnt.ge.2) write(*,*)'diagram8',dot_product(dvo,dvo)
!-------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) -  l2(mn,ef)*r2(ef,in)*t1(a,m) 
!-------------------------------------------------------------

!     call alloc(woo_temp, (nfoo*rcw), id="wop_temp oo")

      woo_temp = 0.0d0

      call contraction_442((/"o2","o3","c1","p3"/),(/"c1","p3","o1","o3" &
     &       /),(/"o2","o1"/),woo_temp,1.0d0,1.0d0,nrep,LeftTensor=l2_ptr,   &
     &          RightTensor=r2_ptr)

      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),&
     &              woo_temp,t1,dvo,-1.0d0,1.0d0,nrep)                

      call dealloc(woo_temp, id="wop_temp oo")

   if (iprnt.ge.2) write(*,*)'diagram9',dot_product(dvo,dvo)
!-----------------------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l2(mn,ef)*r1(e,i)*t2(af,mn)
!-----------------------------------------------------------------------


    call alloc_array(e,nrep)
    call auto_symmetry_offset(e,nv,ncont,.false.,.false.)

    call alloc(wpp_temp, (e%oneNonDirac(1)*rcw), id="wpp_temp #1")

      wpp_temp = 0.0d0

      call contraction_442((/"o2","o3","c1","p3"/),(/"p1","p3","o2","o3" &
     &       /),(/"p1","c1"/),wpp_temp,-1.0d0,1.0d0,nrep,LeftTensor=l2_ptr, & 
     &          RightTensor=t2)

      call contraction_222 ((/"p1","c1"/),(/"c1","o1"/),(/"p1","o1"/), &
     &              wpp_temp,r1_ptr,dvo,1.0d0,1.0d0,nrep)                  

    call dealloc(wpp_temp, id="wpp_temp #1")

  call dealloc_array(e)

   if (iprnt.ge.2) write(*,*)'diagram10',dot_product(dvo,dvo)
!------------------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l2(mn,ef)*t1(f,i)*r2(ea,mn)
!------------------------------------------------------------------

      call alloc(wpp_temp, (nfvv*rcw), id="wpp_temp vv")

      wpp_temp = 0.0d0

!     call contraction_442((/"c1","p1","o2","o3"/),(/"o2","o3","c1","p3"/), &
!    &       (/"p1","p3"/),wpp_temp,-1.0d0,1.0d0,nrep,LeftTensor=r2_ptr, & 
!    &          RightTensor=l2_ptr)

!     call contraction_222 ((/"p1","p3"/),(/"p3","o1"/),(/"p1","o1"/), &
!    &              wpp_temp,t1,dvo,1.0d0,1.0d0,nrep)                  

      call dealloc(wpp_temp, id="wpp_temp vv")


   if (iprnt.ge.2) write(*,*)'diagram10'
!--------------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) +  l2(mn,ef)*t2(af,in)*r1(e,m) 
!--------------------------------------------------------------------

      call alloc(woo_temp, (nfvo*rcw), id="woo_temp vo")

      woo_temp = 0.0d0

      call contraction_422((/"o2","o3","c1","p3"/),(/"c1","o2"/), &
     &       (/"o3","p3"/),r1_ptr,woo_temp,1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)   

      call contraction_422 ((/"p1","p3","o1","o3"/),(/"o3","p3"/),(/"p1","o1"/),&
     &           woo_temp,dvo,1.0d0,1.0d0,nrep,t2)                

!     call dealloc(woo_temp, id="woo_temp vo")

   if (iprnt.ge.2) write(*,*)'diagram11',dot_product(dvo,dvo)
!--------------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) -  l2(mn,ef)*t1(f,i)*t1(a,n)*r1(e,m) 
!--------------------------------------------------------------------

!     call alloc(woo_temp, (nfvo*rcw), id="woo_temp vo")

      woo_temp = 0.0d0

      call contraction_422((/"o2","o3","c1","p3"/),(/"c1","o2"/), &
     &       (/"o3","p3"/),r1_ptr,woo_temp,1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)   

     call alloc(work1, (nfoo*rcw), id="work1 oo") 

       work1 = 0.0d0

      call contraction_222 ((/"p3","o1"/),(/"o3","p3"/),(/"o3","o1"/),&
     &              t1,woo_temp,work1,1.0d0,1.0d0,nrep)                

      call dealloc(woo_temp, id="woo_temp vo")

      call contraction_222 ((/"p1","o3"/),(/"o3","o1"/),(/"p1","o1"/),&
     &              t1,work1,dvo,-1.0d0,1.0d0,nrep)                

     call dealloc(work1, id="work1 oo")

   if (iprnt.ge.2) write(*,*)'diagram12',dot_product(dvo,dvo)
!----------------------------------------------------------
!       dov(i,a) = dov(i,a) + l2(im,ba) * r1(b,m)
!----------------------------------------------------------

      call contraction_422 ((/"o1","o2","c1","p1"/),(/"c1","o2"/), & 
     & (/"o1","p1"/),r1_ptr,dov,-1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)

   if (iprnt.ge.2) write(*,*)'diagram13',dot_product(dov,dov)
      if (carith) then
         call conjuga (nfoo,doo,1)
         call conjuga (nfvv,dvv,1)
         call conjuga (nfvo,dvo,1)
         call conjuga (nfvo,dov,1)
      endif


     call alloc(work1, (nfvv*rcw), id="work1 vv")
     call alloc(work2, (nfvv*rcw), id="work2 vv")

      call convdm(work1,work2,doo,dvv,dov,dvo,dmo)

     call dealloc(work1, id="work1 vv")
     call dealloc(work2, id="work2 vv")

    call dealloc(doo, id="doo")
    call dealloc(dvv, id="dvv")
    call dealloc(dvo, id="dvo")
    call dealloc(dov, id="dov")

!   deallocate(dov)
!   deallocate(dvo)
!   deallocate(doo)
!   deallocate(dvv)

      end subroutine

  subroutine density_EE (l1,l2,r1,r2,t1,t2,dmo)

  use symmetry_offset
#include "symm.inc"
#include "complex.inc"
#include "inpt.inc"
!----------------calling variables-------------------------
   real*8, intent(in),contiguous,target :: l1(:),l2(:),r1(:),r2(:)
   real*8, intent(in) :: t1(:),t2(:)
   real*8, intent(inout) :: dmo(*)
!----------------local variables---------------------------
   real*8, allocatable :: doo(:),dvv(:),dvo(:),dov(:)  
   real*8, allocatable :: work1(:),work2(:) 
   real*8, allocatable :: woo_temp(:),wpp_temp(:),wop_temp(:)
   real(8), pointer,contiguous :: l1_ptr(:),l2_ptr(:),r1_ptr(:),r2_ptr(:)
   type(Offset)         :: e
   integer              :: ndiml1,ndiml2
!-----------------------------------------------------------

  ndiml1 = size(l1,1) 
  ndiml2 = size(l2,1) 

  l1_ptr(1:ndiml1) => l1
  l2_ptr(1:ndiml2) => l2
  r1_ptr(1:ndiml1) => r1
  r2_ptr(1:ndiml2) => r2

    allocate(doo(nfoo*rcw))
    allocate(dvv(nfvv*rcw))
    allocate(dvo(nfvo*rcw))
    allocate(dov(nfvo*rcw))

     doo = 0.0d0
     dvv = 0.0d0
     dvo = 0.0d0
     dov = 0.0d0

!-------------------------------------------------------------------
!   doo(j,i) = doo(j,i) - l2(jm,ef) * t2(ef,im)
!-------------------------------------------------------------------
        call contraction_442 ((/"o2","o3","p1","p2"/), &
     &    (/"p1","p2","o1","o3"/),(/"o2","o1"/),doo,-1.0d0,1.0d0,nrep, &
     &    LeftTensor=l2,RightTensor=t2)

   if (iprnt.ge.2) write(*,*)'diagram1',dot_product(doo,doo)

!-------------------------------------------------------------------
!   doo(j,i) = doo(j,i) - l2(jm,ef) * r2(ef,im)
!-------------------------------------------------------------------

        call contraction_442 ((/"o2","o3","p1","p2"/),&
     &    (/"p1","p2","o1","o3"/),(/"o2","o1"/),doo,-1.0d0,1.0d0,nrep,&
     &    LeftTensor=l2,RightTensor=r2)


   if (iprnt.ge.2) write(*,*)'diagram2',dot_product(doo,doo)

!---------------------------------------------------------------------
!    doo(j,i) = doo(j,i) - l1(j,e)*t1(e,i)
!---------------------------------------------------------------------

       call contraction_222 ((/"o2","p1"/),(/"p1","o1"/),(/"o2","o1"/),&
     &                         l1,t1,doo,-1.0d0,1.0d0,nrep)



   if (iprnt.ge.2) write(*,*)'diagram3',dot_product(doo,doo)


!---------------------------------------------------------------------
!    doo(j,i) = doo(j,i) - l1(j,e)*r1(e,i)
!---------------------------------------------------------------------

       call contraction_222 ((/"o2","p1"/),(/"p1","o1"/),(/"o2","o1"/),&
     &                         l1,r1,doo,-1.0d0,1.0d0,nrep)


   if (iprnt.ge.2) write(*,*)'diagram4',dot_product(doo,doo)
!---------------------------------------------------------------------
!    doo(j,i) = doo(j,i) - l2(im,ae)*r1(a,m)*t1(e,j) 

!    wop_temp(i,e) = wop_temp(i,e) + l2(im,ae)*r1(a,m) 
!    doo(j,i) = doo(j,i) + wop_temp(i,e) * t1(e,j)
!---------------------------------------------------------------------

  allocate(wop_temp(nfvo*rcw)) 

  wop_temp = 0.0d0

  call contraction_422 ((/"o1","o3","p1","p2"/),(/"p1","o3"/),(/"o1","p2"/),r1,wop_temp,1.0d0,1.0d0,nrep,LeftTensor=l2)

  call contraction_222 ((/"p2","o2"/),(/"o1","p2"/),(/"o2","o1"/),t1,wop_temp,doo,1.0d0,1.0d0,nrep)

 deallocate(wop_temp)


   if (iprnt.ge.2) write(*,*)'diagram5',dot_product(doo,doo)
!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l2(mn,ae)*t2(be,mn) 
!----------------------------------------------------------------------

      call contraction_442 ((/"o1","o2","p1","p3"/),&
     &    (/"p2","p3","o1","o2"/),(/"p2","p1"/),dvv,1.0d0,1.0d0,nrep,&
     &     LeftTensor=l2,RightTensor=t2 )

   if (iprnt.ge.2) write(*,*)'diagram6',dot_product(dvv,dvv)

!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l2(mn,ae)*r2(be,mn) 
!----------------------------------------------------------------------

      call contraction_442 ((/"o1","o2","p1","p3"/),&
     &    (/"p2","p3","o1","o2"/),(/"p2","p1"/),dvv,1.0d0,1.0d0,nrep,&
     &     LeftTensor=l2,RightTensor=r2 )

   if (iprnt.ge.2) write(*,*)'diagram7',dot_product(dvv,dvv)

!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l1(m,a)*t1(b,m)
!----------------------------------------------------------------------

      call contraction_222 ((/"o1","p1"/),(/"p2","o1"/),(/"p2","p1"/),&
     &              l1,t1,dvv,1.0d0,1.0d0,nrep)

   if (iprnt.ge.2) write(*,*)'diagram8',dot_product(dvv,dvv)
!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l1(m,a)*r1(b,m)
!----------------------------------------------------------------------

      call contraction_222 ((/"o1","p1"/),(/"p2","o1"/),(/"p2","p1"/),&
     &              l1,r1,dvv,1.0d0,1.0d0,nrep)


   if (iprnt.ge.2) write(*,*)'diagram9',dot_product(dvv,dvv)
!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l2(mn,cb)*r1(c,n)*t1(a,m)
!     wop_temp(m,b) = wop_temp(m,b) - l2(mn,cb)*r1(c,n) 
!     dvv(b,a) = dvv(b,a) + wop_temp(m,b)*t1(a,m)
!----------------------------------------------------------------------

  allocate(wop_temp(nfvo*rcw))

  wop_temp = 0.0d0

  call contraction_422 ((/"o1","o2","p3","p2"/),(/"p3","o2"/),(/"o1","p2"/),r1,wop_temp,-1.0d0,1.0d0,nrep,LeftTensor=l2)

  call contraction_222 ((/"o1","p2"/),(/"p1","o1"/),(/"p2","p1"/),wop_temp,t1,dvv,-1.0d0,1.0d0,nrep)

  deallocate(wop_temp)


   if (iprnt.ge.2) write(*,*)'diagram10',dot_product(dvv,dvv)


!     dvo(a,i) = t1(a,i)
!----------------------------------------------------------------------------

      call xcopy (nfvo,t1,1,dvo,1)

   if (iprnt.ge.2) write(*,*)'diagram11',dot_product(dvo,dvo)
!--------------------------------------------------------------------------
!     dvo(a,i) = dvo(a,i) + l1(m,e) * t(ae,im)
!--------------------------------------------------------------------------

!     call contraction_242 ((/"o2","p2"/),(/"p1","p2","o1","o2"/), 
!    & (/"p1","o1"/),l1,dvo,1.0d0,1.0d0,nrep,RightTensor=t2)


      call contraction_422 ((/"p1","p2","o1","o2"/),(/"o2","p2"/), &
     & (/"p1","o1"/),l1,dvo,1.0d0,1.0d0,nrep,LeftTensor=t2)

   if (iprnt.ge.2) write(*,*)'diagram12',dot_product(dvo,dvo)


!--------------------------------------------------------------------------
!     dvo(a,i) = dvo(a,i) + l1(m,e) * r(ae,im)
!--------------------------------------------------------------------------

!     call contraction_242 ((/"o2","p2"/),(/"p1","p2","o1","o2"/), 
!    & (/"p1","o1"/),l1,dvo,1.0d0,1.0d0,nrep,RightTensor=t2)


      call contraction_422 ((/"p1","p2","o1","o2"/),(/"o2","p2"/), &
     & (/"p1","o1"/),l1,dvo,1.0d0,1.0d0,nrep,LeftTensor=r2)


   if (iprnt.ge.2) write(*,*)'diagram13',dot_product(dvo,dvo)

!----------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l1(m,e) * t(e,i)*t(a,m)
!----------------------------------------------------------

      allocate(woo_temp(nfoo*rcw))

      woo_temp = 0.0d0

      call contraction_222 ((/"o2","p2"/),(/"p2","o1"/),(/"o2","o1"/), &
     &              l1,t1,woo_temp,1.0d0,1.0d0,nrep)


      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),&
     &              woo_temp,t1,dvo,-1.0d0,1.0d0,nrep)

      deallocate(woo_temp)


   if (iprnt.ge.2) write(*,*)'diagram14',dot_product(dvo,dvo)

!----------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l1(m,e) * r(e,i)*t(a,m)
!----------------------------------------------------------

      allocate(woo_temp(nfoo*rcw))

      woo_temp = 0.0d0

      call contraction_222 ((/"o2","p2"/),(/"p2","o1"/),(/"o2","o1"/),&
     &              l1,r1,woo_temp,1.0d0,1.0d0,nrep)


      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),&
     &              woo_temp,t1,dvo,-1.0d0,1.0d0,nrep)

      deallocate(woo_temp)


   if (iprnt.ge.2) write(*,*)'diagram15',dot_product(dvo,dvo)
!----------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l1(m,e) * t(e,i)*r(a,m)
!----------------------------------------------------------

      allocate(woo_temp(nfoo*rcw))

      woo_temp = 0.0d0

      call contraction_222 ((/"o2","p2"/),(/"p2","o1"/),(/"o2","o1"/),&
     &              l1,t1,woo_temp,1.0d0,1.0d0,nrep)


      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),&
     &              woo_temp,r1,dvo,-1.0d0,1.0d0,nrep)

      deallocate(woo_temp)


   if (iprnt.ge.2) write(*,*)'diagram16',dot_product(dvo,dvo)

!-------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) -  l2(mn,ef)*t2(ef,in)*t1(a,m) 
!-------------------------------------------------------------

      allocate(woo_temp(nfoo*rcw))

      woo_temp = 0.0d0

      call contraction_442((/"o2","o3","p2","p3"/),(/"p2","p3","o1","o3" &
     &       /),(/"o2","o1"/),woo_temp,1.0d0,1.0d0,nrep,LeftTensor=L2, & 
     &          RightTensor=t2)

      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),&
     &              woo_temp,t1,dvo,-1.0d0,1.0d0,nrep)

      deallocate(woo_temp)

   if (iprnt.ge.2) write(*,*)'diagram17',dot_product(dvo,dvo)

!-------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) -  l2(mn,ef)*r2(ef,in)*t1(a,m) 
!-------------------------------------------------------------

      allocate(woo_temp(nfoo*rcw))

      woo_temp = 0.0d0

      call contraction_442((/"o2","o3","p2","p3"/),(/"p2","p3","o1","o3" &
     &       /),(/"o2","o1"/),woo_temp,1.0d0,1.0d0,nrep,LeftTensor=L2, &
     &          RightTensor=r2)

      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),&
     &              woo_temp,t1,dvo,-1.0d0,1.0d0,nrep)

      deallocate(woo_temp)

   if (iprnt.ge.2) write(*,*)'diagram18',dot_product(dvo,dvo)


!-------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) -  l2(mn,ef)*t2(ef,in)*r1(a,m) 
!-------------------------------------------------------------

      allocate(woo_temp(nfoo*rcw))

      woo_temp = 0.0d0

      call contraction_442((/"o2","o3","p2","p3"/),(/"p2","p3","o1","o3" &
     &       /),(/"o2","o1"/),woo_temp,1.0d0,1.0d0,nrep,LeftTensor=L2, &
     &          RightTensor=t2)

      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),&
     &              woo_temp,r1,dvo,-1.0d0,1.0d0,nrep)

      deallocate(woo_temp)


   if (iprnt.ge.2) write(*,*)'diagram19',dot_product(dvo,dvo)

!-----------------------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l2(mn,ef)*t1(e,i)*t2(af,mn)
!-----------------------------------------------------------------------

      allocate(wpp_temp(nfvv*rcw))

      wpp_temp = 0.0d0

      call contraction_442((/"o2","o3","p2","p3"/),(/"p1","p3","o2","o3" &
     &       /),(/"p1","p2"/),wpp_temp,-1.0d0,1.0d0,nrep,LeftTensor=L2, &
     &          RightTensor=t2)

      call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/),&
     &              wpp_temp,t1,dvo,1.0d0,1.0d0,nrep)

      deallocate(wpp_temp)


   if (iprnt.ge.2) write(*,*)'diagram20',dot_product(dvo,dvo)
!-----------------------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l2(mn,ef)*r1(e,i)*t2(af,mn)
!-----------------------------------------------------------------------

      allocate(wpp_temp(nfvv*rcw))

      wpp_temp = 0.0d0

      call contraction_442((/"o2","o3","p2","p3"/),(/"p1","p3","o2","o3" &
     &       /),(/"p1","p2"/),wpp_temp,-1.0d0,1.0d0,nrep,LeftTensor=L2, &
     &          RightTensor=t2)

      call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/), &
     &              wpp_temp,r1,dvo,1.0d0,1.0d0,nrep)

      deallocate(wpp_temp)

   if (iprnt.ge.2) write(*,*)'diagram21',dot_product(dvo,dvo)

!-----------------------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l2(mn,ef)*t1(e,i)*r2(af,mn)
!-----------------------------------------------------------------------

      allocate(wpp_temp(nfvv*rcw))

      wpp_temp = 0.0d0

      call contraction_442((/"o2","o3","p2","p3"/),(/"p1","p3","o2","o3" &
     &       /),(/"p1","p2"/),wpp_temp,-1.0d0,1.0d0,nrep,LeftTensor=L2, &
     &          RightTensor=r2)

      call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/), &
     &              wpp_temp,t1,dvo,1.0d0,1.0d0,nrep)

      deallocate(wpp_temp)


   if (iprnt.ge.2) write(*,*)'diagram22',dot_product(dvo,dvo)
!--------------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) +  l2(mn,ef)*t2(af,in)*r1(e,m) 
!--------------------------------------------------------------------

      allocate(woo_temp(nfvo*rcw))

      woo_temp = 0.0d0

      call contraction_422((/"o2","o3","p4","p3"/),(/"p4","o2"/), &
     &       (/"o3","p3"/),r1_ptr,woo_temp,1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)   

      call contraction_422 ((/"p1","p3","o1","o3"/),(/"o3","p3"/),(/"p1","o1"/),&
     &           woo_temp,dvo,1.0d0,1.0d0,nrep,t2)                

      deallocate(woo_temp)

   if (iprnt.ge.2) write(*,*)'diagram23',dot_product(dvo,dvo)
!--------------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) -  l2(mn,ef)*t1(f,i)*t1(a,n)*r1(e,m) 
!--------------------------------------------------------------------

      allocate(woo_temp(nfvo*rcw))

      woo_temp = 0.0d0

      call contraction_422((/"o2","o3","p4","p3"/),(/"p4","o2"/), &
     &       (/"o3","p3"/),r1_ptr,woo_temp,1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)   

     allocate(work1(nfoo*rcw)) 

       work1 = 0.0d0

      call contraction_222 ((/"p3","o1"/),(/"o3","p3"/),(/"o3","o1"/),&
     &              t1,woo_temp,work1,1.0d0,1.0d0,nrep)                

      deallocate(woo_temp)

      call contraction_222 ((/"p1","o3"/),(/"o3","o1"/),(/"p1","o1"/),&
     &              t1,work1,dvo,-1.0d0,1.0d0,nrep)                

     deallocate(work1)

   if (iprnt.ge.2) write(*,*)'diagram24',dot_product(dvo,dvo)

!----------------------------------------------------------
!       dov(i,a) = l1(i,a)
!----------------------------------------------------------

      call xcopy (nfvo,l1,1,dov,1)

!----------------------------------------------------------
!       dov(i,a) = dov(i,a) + l2(im,ba) * r1(b,m)
!----------------------------------------------------------

      call contraction_422 ((/"o1","o2","p2","p1"/),(/"p2","o2"/), & 
     & (/"o1","p1"/),r1_ptr,dov,-1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)

   if (iprnt.ge.2) write(*,*)'diagram25',dot_product(dov,dov)
      if (carith) then
         call conjuga (nfoo,doo,1)
         call conjuga (nfvv,dvv,1)
         call conjuga (nfvo,dvo,1)
         call conjuga (nfvo,dov,1)
      endif

     allocate(work1(nfvv*rcw))
     allocate(work2(nfvv*rcw))

      work1 = 0.0d0
      work2 = 0.0d0

      call convdm(work1,work2,doo,dvv,dov,dvo,dmo)

     deallocate(work1,work2)  


    deallocate(dov)
    deallocate(dvo)
    deallocate(doo)
    deallocate(dvv)

   end subroutine

   end module 
