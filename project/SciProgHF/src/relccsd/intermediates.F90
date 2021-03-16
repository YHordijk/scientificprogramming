  module intermediates_1b_2b

! Module to make intermediates needed in the CC and Lambda equations and for EOMCC
! Written by Avijit Shee and Lucas Visscher, 2014

  use interface_to_mpi
  use contraction
  use interface_to_mpi
  use wrapper_sorting

    implicit none

    public f_bar_mi
    public f_bar_me
    public f_bar_ae   
    public w_ijmn
    public w_ejmb
    public w_mbej
    public w_iemn
    public w_efam
    public w_amef
    public w_mnie
    public wbar_mbej
    public w_ejmb_out_of_core
    public f_bar_ae_out_of_core
    public w_iemn_out_of_core
    public w_efam_diagram2
    public w_efab_out_of_core
    public G_intm
     private

     contains

   subroutine f_bar_mi(foo,fvo,t1,t2,ho)

#include "symm.inc"
#include "param.inc"
#include "complex.inc"
!--------------Calling variables---------------------------------------
     real(8), intent(in) :: foo(nfoo*rcw), fvo(rcw*nfvo)   ! Fock matrix elements
     real(8), intent(in) :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes
     real(8), intent(inout):: ho(nfoo*rcw)                 ! Intermediate
!---------------Local variables----------------------------------------
     real(8), allocatable :: tau(:) !, vooo(:),vooo_sorted(:),vvoo(:)

!    Initialize  f_bar_mi with oo block of Fock matrix
       ho(1:nfoo*rcw) = foo(1:nfoo*rcw)

!-------------Term-1----------------------
! F_bar_mi = F_bar_mi - V(nm,ef) * T(ef,in) 
!-----------------------------------------

   allocate(tau(ndimt2*rcw))

   call gettau  (t1,t2,tau)

   call contraction_442 ((/"o3","o2","p1","p2"/),(/"p1","p2","o1","o3"/),(/"o2","o1"/),ho,-1.0d0,1.0d0,nrep,RightTensor=tau)

   deallocate(tau) 

!------------------------TERM-2---------------------
!F_bar_mi(m,i) = F_bar_mi(m,i) - V(mn, ei) * T1(e,n) 
!---------------------------------------------------

    call contraction_422 ((/"o2","o3","p1","o1"/),(/"p1","o3"/),(/"o2","o1"/),t1,ho,-1.0d0,1.0d0,nrep)

!-----------------------TERM-3---------
!F_bar_mi = F_bar_mi + F(m,e) * T1(e,i)
!--------------------------------------
   
   call contraction_222 ((/"o2","p1"/),(/"p1","o1"/),(/"o2","o1"/),fvo,t1,ho,1.0d0,1.0d0,nrep)

   end subroutine

   subroutine f_bar_ae(fvv,fvo,t1,t2,hv)

!---------------Common Blocks--------------------------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
      real(8), intent(in)  :: fvv(nfvv*rcw), fvo(nfvo*rcw)
      real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
      real(8), intent(inout) :: hv(nfvv*rcw)
!---------------Local variables--------------------------------------
      real(8), allocatable :: tau(:)
      real(8), allocatable :: fov(:)
!---------------------------------------
! hv(a,e) = fvv(a,e)  
!---------------------------------------

    hv = fvv

!----------------------------------------------------------------------
! hv(a,e) = hv(a,e) - V(mn,ef) * Tau(af,mn)
!----------------------------------------------------------------------

   allocate(tau(ndimt2*rcw))

   call gettau (t1,t2,tau)

   call contraction_442 ((/"o1","o2","p2","p3"/),(/"p1","p3","o1","o2"/),(/"p1","p2"/),hv,-1.0d0,1.0d0,nrep,RightTensor=tau)
  
   deallocate(tau) 

!---------------------------------------------------------------
!F_bar_ae = F_bar_ae + V(am, ef) * T1(f,m) -> V(ae,fm) * T1(fm) 
!---------------------------------------------------------------

  call contraction_422 ((/"p1","o1","p2","p3"/),(/"p3","o1"/),(/"p1","p2"/),t1,hv,1.0d0,1.0d0,nrep)

!-----------------------------------------------------------------
!F_bar_ae = F_bar_ae + F(m,a) * T1(e,m)
!-----------------------------------------------------------------

   allocate (fov(nfvo*rcw))
   call srt1c1(nrep,nv,no,fvo,fov)

   call contraction_222 ((/"o1","p1"/),(/"p2","o1"/),(/"p1","p2"/),fov,t1,hv,-1.0d0,1.0d0,nrep)

   deallocate (fov)

  end subroutine

  subroutine f_bar_me(fvo,t1,hov)

#include "symm.inc"
#include "param.inc"
#include "complex.inc"
!++++++++++++++Calling variables+++++++++++++++++++++++++++++++++++++++
     real(8), intent(in),target  :: fvo(rcw*nfvo)          ! Fock matrix elements
     real(8), intent(in) :: t1(rcw*ndimt1)                 ! Amplitudes
     real(8), intent(out):: hov(rcw*nfvo)                  ! Intermediate
!+++++initialization of hov+++++++++++++
 
  call srt1c1(nrep,nv,no,fvo,hov)

!  hov = fvo

!+++++++++++++Term-1++++++++++++++++++++++
! F_bar_me = F_bar_me + V(mn,ef) * T1(f,n) 
!+++++++++++++++++++++++++++++++++++++++++

 call contraction_422 ((/"o1","o2","p1","p2"/),(/"p2","o2"/),(/"o1","p1"/),t1,hov,1.0d0,1.0d0,nrep)

  end subroutine

  subroutine w_ijmn(t1,t2,w_oooo)

!---------------Description--------------------------------------------
!     Calculates w_oooo intermediate
!-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_oooo(rcw*nv1)                ! Intermediate.
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
!--------------------------------------------------

!--------------------
! initialize w_oooo
!--------------------

   call getoooo (w_oooo)

!*************************************************
! w_oooo(ij,mn) = w_oooo(ij,mn) + (V(ij,en)*T(e,m))
!*************************************************

   call contraction_424 ((/"o1","o2","p1","o4"/),(/"p1","o3"/),(/"o1","o2","o4","o3"/),w_oooo, &
  &      1.0d0,-1.0d0,nrep,RightTensor=t1)

!******************************************************
! w_oooo(ij,mn) = w_oooo(ij,mn) +  V(ij,ef) * Tau(ef,mn)
!******************************************************

   allocate(tau(ndimt2*rcw))
   call gettau (t1,t2,tau)
   call contraction_444 ((/"o1","o2","p1","p2"/),(/"p1","p2","o3","o4"/),(/"o1","o2","o3","o4"/),w_oooo,1.0d0,1.0d0,nrep, &
  &                      RightTensor=tau)

   deallocate(tau)

  end subroutine
 
!===============================================================================================

 subroutine w_ejmb(t1,t2,w_voov)

!---------------Description--------------------------------------------
!     Calculates w_voov intermediate
!-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_voov(rcw*ivoov(nrep+1))                ! Intermediate.
!---------------Local Variables----------------------------------------
 real*8,allocatable      :: array_temp(:)
!--------------------
! initialize w_voov
!--------------------

 call getvoov (w_voov)

!---------------------------------------------------
!   w_voov(ej,mb) = w_voov(ej,mb) - V(nj,mb)*T1(e,n)
!---------------------------------------------------

   call contraction_244 ((/"p2","o3"/),(/"o3","o1","o2","p1"/),(/"p2","o1","o2","p1"/),w_voov, &
            &      -1.0d0,1.0d0,nrep,LeftTensor=t1)

!----------------------------------------------------------------------
!   w_voov(ej,mb) = w_voov(ej,mb) + V(nj,fb)*(T2(ef,mn)-T1(f,m)*T1(e,n)
!----------------------------------------------------------------------

  call contraction_444 ((/"p2","p3","o2","o3"/),(/"o3","o1","p3","p1"/),(/"p2","o1","o2","p1"/),w_voov, &
   &      1.0d0,0.0d0,nrep,LeftTensor=t2)


  allocate(array_temp(iovoot(nrep+1)*rcw))

   array_temp = 0.0d0

  call contraction_424 ((/"o3","o1","p1","p3"/),(/"p3","o2"/),(/"o3","o1","o2","p1"/),array_temp, &
             &      -1.0d0,1.0d0,nrep,RightTensor=t1)

  call contraction_244 ((/"p2","o3"/),(/"o3","o1","o2","p1"/),(/"p2","o1","o2","p1"/),w_voov, &
             &      -1.0d0,1.0d0,nrep,LeftTensor=t1,RightTensor=array_temp)

  deallocate(array_temp)


!--------------------------------------------------
!  w_voov(ej,mb) = w(ej,mb) + V(ej,fb)*t(f,m) 
!--------------------------------------------------

  call contraction_424 ((/"p2","o1","p1","p3"/),(/"p3","o2"/),(/"p2","o1","o2","p1"/),w_voov, &
             &      -1.0d0,1.0d0,nrep,RightTensor=t1)

 end subroutine

!===========================================================================================================

 subroutine w_mbej(t1,t2,w_ovvo)
  
  external srt1l1
  
!---------------Description--------------------------------------------
!     Calculates w_voov intermediate
!-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_ovvo(rcw*ivoov(nrep+1))                ! Intermediate.
!---------------Local Variables----------------------------------------
 real(8), allocatable :: w_vovo(:),array_temp(:)
!--------------------
! initialize w_voov
!--------------------

  allocate (w_vovo(nv4*rcw))

  call getvovo (w_vovo)

  call srt1l1 (nrep,multb,.false.,nvo,nv,no,nvo,iovvo,iiov,w_vovo,w_ovvo) 

  deallocate (w_vovo)

!---------------------------------------------------
!   w_ovvo(mb,ej) = w_voov(mb,ej) - V(mn,ej)*T1(b,n)
!---------------------------------------------------

   call contraction_424 ((/"o2","o3","p2","o1"/),(/"p1","o3"/),(/"o2","p1","p2","o1"/),w_ovvo, &
            &      -1.0d0,1.0d0,nrep,RightTensor=t1)

!-------------------------------------------------------------------------------------
!  w_ovvo(mb,ej) = w_ovvo(mb,ej) - V(mn,ef)*(T2(bf,nj)+T1(f,j)*T1(b,n))
!-------------------------------------------------------------------------------------


  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o3","o1"/),(/"o2","p1","p2","o1"/),w_ovvo, &
             &      -1.0d0,0.0d0,nrep,RightTensor=t2)


  allocate(array_temp(nv2*rcw))

  array_temp = 0.0d0

  call contraction_424 ((/"o2","o3","p2","p3"/),(/"p3","o1"/),(/"o2","o3","p2","o1"/),array_temp, &
             &      1.0d0,1.0d0,nrep,RightTensor=t1)


  call contraction_424 ((/"o2","o3","p2","o1"/),(/"p1","o3"/),(/"o2","p1","p2","o1"/),w_ovvo, &
             &      -1.0d0,1.0d0,nrep,LeftTensor=array_temp,RightTensor=t1)

  deallocate(array_temp)

!--------------------------------------------------
!  w_ovvo(mb,ej) = w(mb,ej) + V(mb,ef)*t(f,j) 
!--------------------------------------------------

  call contraction_424 ((/"o2","p1","p2","p3"/),(/"p3","o1"/),(/"o2","p1","p2","o1"/),w_ovvo, &
             &      1.0d0,1.0d0,nrep,RightTensor=t1)


 end subroutine

!-----=========================--------------

 subroutine w_iemn(t1,t2,w_ovoo,wbar_ovvo,w_oooo,hov)

!---------------Description--------------------------------------------
!     Calculates w_voov intermediate
!-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)                                    ! Amplitudes.
   real(8), intent(out) :: w_ovoo(iovoot(nrep+1)*rcw)                                       ! Intermediate.
   real(8), intent(in)  :: wbar_ovvo(rcw*iovvo(nrep+1)),w_oooo(rcw*nv1),hov(nfvo*rcw)       ! reuse of intermediates.
!---------------Local Variables----------------------------------------
  real*8, allocatable   :: tau(:) 
!--------------------
! initialize wbar_ovvo
!--------------------

  call getovoo (w_ovoo) 

!-----------------------------------------------------
! w_ovoo(ie,mn) = w_ovoo(ie,mn) + F_me(if) * T2(ef,mn)  
!-----------------------------------------------------

  call contraction_244 ((/"o1","p2"/),(/"p2","p1","o2","o3"/),(/"o1","p1","o2","o3"/),w_ovoo, &
           &       1.0d0,1.0d0,nrep,LeftTensor=hov,RightTensor=t2)

!-------------------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) +  V(io,mf) * T2(ef,no)
!-------------------------------------------------------

  call contraction_444 ((/"o1","o4","o2","p2"/),(/"p1","p2","o3","o4"/),(/"o1","p1","o2","o3"/),w_ovoo, &
             &      1.0d0,0.0d0,nrep,RightTensor=t2)

!-------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) +  w_bar_mbej(ie,fn)*T1(fm)
!--------------------------------------------

   call contraction_424 ((/"o1","p1","p2","o3"/),(/"p2","o2"/),(/"o1","p1","o3","o2"/),w_ovoo, &
             &      1.0d0,-1.0d0,nrep,LeftTensor=wbar_ovvo,RightTensor=t1)

!-----------------------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) + w_ijmn(io,mn) * t(e,o)
!-----------------------------------------------------------

   call contraction_424 ((/"o1","o4","o2","o3"/),(/"p1","o4"/),(/"o1","p1","o2","o3"/),w_ovoo, &
             &      -1.0d0,1.0d0,nrep,LeftTensor=w_oooo,RightTensor=t1)

!-----------------------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) + V(ie,fg) * tau(fg,mn) 
!-----------------------------------------------------------
!**out of core term..

  allocate (tau(ndimt2*rcw))

  call gettau (t1,t2,tau)

  call contraction_444 ((/"o1","p1","p2","p3"/),(/"p2","p3","o2","o3"/),(/"o1","p1","o2","o3"/),w_ovoo, &
             &      1.0d0,1.0d0,nrep,RightTensor=tau)

  deallocate(tau)

  end subroutine

 subroutine wbar_mbej(t2,wbar_ovvo)

  external srt1l1
 !---------------Description--------------------------------------------
 !     Calculates w_voov intermediate
 !-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

 !---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t2(rcw*ndimt2)                              ! Amplitudes.
   real(8), intent(out) :: wbar_ovvo(rcw*iovvo(nrep+1))                ! Intermediate.
 !---------------Local Variables----------------------------------------
  real*8, allocatable   :: w_vovo(:) 
 !--------------------
 ! initialize wbar_ovvo
 !--------------------

  allocate (w_vovo(nv4*rcw))

  call getvovo (w_vovo)

  call srt1l1 (nrep,multb,.false.,nvo,nv,no,nvo,iovvo,iiov,w_vovo,wbar_ovvo) 

  deallocate (w_vovo)

!------------------------------------------------------------
!   wbar_ovvo(mb,ej) = wbar_ovvo(mb,ej) - V(mn,ef) * T2(bf,nj)
!------------------------------------------------------------

  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o1","o3"/),(/"o2","p1","p2","o1"/),wbar_ovvo, &
             &      1.0d0,0.0d0,nrep,RightTensor=t2)

 end subroutine

 subroutine w_mnie(t1,int_wooop)

!---------------Description--------------------------------------------
!     Calculates w_mnie intermediate
!-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables----------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1)                                    ! Amplitudes.
  real(8), intent(out) :: int_wooop(:)                                      ! Intermediate.
!------------------------------------------------------------
! initialize int_wooop
!------------------------------------------------------------

  call getooov(int_wooop)

!----------------------------------------------------------------------------
!     w(mn,ie) = V(mn,ie) + t(f,i) * V(mn,fe)
!----------------------------------------------------------------------------

  call contraction_424 ((/"o1","o2","p3","p2"/),(/"p3","o3"/), &
 &        (/"o1","o2","o3","p2"/),int_wooop,1.0d0,1.0d0,nrep,RightTensor=t1)

 end subroutine

 subroutine w_amef(t1,w_v_temp)

!---------------Description--------------------------------------------
!     Calculates w_amef intermediate
!-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables----------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1)                                    ! Amplitudes.
  real(8), intent(out) :: w_v_temp(:)                                      ! Intermediate.
!--------------------------------------------------
! initialize w_v_temp
!--------------------------------------------------

    call getvovv_incore(w_v_temp)

!---------------------------------------------------
! w_amef(am,ef) = w_amef(am,ef) - V(nm,ef) * t(a,n) 
!---------------------------------------------------

    call contraction_244 ((/"p2","o2"/),(/"o2","o1","p3","p1"/), &
   &    (/"p2","o1","p3","p1"/),w_v_temp,-1.0d0,1.0d0,nrep,LeftTensor=t1)

 end subroutine

 subroutine w_efam(t1,t2,wbar_ovvo,hov,w_vvvo)

!--------------------DESCRIPTION---------------------------
!        calculates vvvo type intermediate in-core fashion.      
!----------------------------------------------------------

!----------calling variables----------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "files.inc"

!---------------Calling variables--------------------------------------
   real*8, intent(in)  :: t2(rcw*ndimt2),t1(rcw*ndimt1)         ! Amplitudes.
   real*8, intent(in)  :: wbar_ovvo(rcw*iovvo(nrep+1)),hov(nfvo*rcw)     ! Intermediate.
   real*8, intent(out) :: w_vvvo(rcw*nv5)                ! Intermediate.
!---------------Local Variables----------------------------------------

   integer :: AllocateStatus 
   real*8, allocatable  :: tau(:)                
   real*8, allocatable  :: vovv1(:),temp(:)                

   type(indices) :: e,f,g

!----------------------------------------------------------
!        Initialize
!----------------------------------------------------------

   call getvvvo_incore(w_vvvo)
!-----------------------------------------------------------
!   diagram1 : w_vvvo(ef,am)  = w_vvvo(ef,am) + p(ef) V(en,ag)*t2(fg,mn)  
!-----------------------------------------------------------

  call contraction_444 ((/"p2","o2","p4","p1"/),(/"p3","p4","o1","o2"/),(/"p2","p3","p1","o1"/),w_vvvo, &
             &      -1.0d0,0.0d0,nrep,RightTensor=t2)

!------------------------------------------------------------
!  diagram2 : w_vvvo(ef,am) =  w_vvvo(ef,am) + w (ef,ag)*t(g,m) 
!------------------------------------------------------------


!-------------------------------------------------------------
!  diagram3 : w_vvvo(ef,am) =  w_vvvo(ef,am) + f_bar_me(na) * t(ef,nm)
!-------------------------------------------------------------

      call contraction_424 ((/"p2","p3","o2","o1"/),(/"o2","p1"/),(/"p2","p3","p1","o1"/),w_vvvo, &
     &      -1.0d0,1.0d0,nrep,RightTensor=hov,LeftTensor=t2)


!-----------------------------------------------------------
!   diagram4 : w_vvvo(ef,am) =  w_vvvo(ef,am) + V(no,am) * tau(ef,no) 
!-----------------------------------------------------------

  allocate (tau(ndimt2*rcw))
  call gettau (t1,t2,tau)

   call contraction_444 ((/"p2","p3","o2","o3"/),(/"o2","o3","p1","o1"/),(/"p2","p3","p1","o1"/),w_vvvo, &
  &      1.0d0,1.0d0,nrep,LeftTensor=tau)

  deallocate (tau)

!---------------------------------------------------------------------------
!  diagram5 : w_vvvo(ef,am) =  w_vvvo(ef,am) - w_bar_ovvo(nf,am) * t(e,n)
!---------------------------------------------------------------------------

  call contraction_244 ((/"p2","o2"/),(/"o2","p3","p1","o1"/),(/"p2","p3","p1","o1"/),w_vvvo, &
            &      -1.0d0,1.0d0,nrep,RightTensor=wbar_ovvo,LeftTensor=t1)

 end subroutine

 subroutine w_efam_diagram2(t1,t2,buf2,buf3,nbuf3,w_vvvo)

#include "param.inc"  
#include "complex.inc"
#include "symm.inc"   
#include "ccpar.inc"
#include "files.inc"
!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
   real(8), intent(inout) :: w_vvvo(rcw*nv5)                ! Intermediate.
   real(8), intent(in)    :: buf2(*)                        ! buffer arrays
   real(8), intent(inout) :: buf3(*)                        ! buffer arrays
   integer, intent(in)    :: nbuf3 

!---------------Local variables--------------------------------------
                                                                    
      real*8  :: t1val(2)                                                
      logical :: done
      integer :: cd,cco,istart,nint,akcd,ci,akci,t1di,di,akdi,t1ci,cmin,c,d,i,j
      integer :: crep,drep,irp,mint
      real*8,allocatable  :: tau (:)

!------------------------------------------------------------------
!  w_vvvo(ef,am) =  w_vvvo(ef,am) + w (ef,ag)*t(g,m)
!------------------------------------------------------------------

      t1val = (/0.0d0,0.0d0/)

     allocate (tau(ndimt2*rcw))
     call gettau (t1,t2,tau)

      do 30 irp = 1, nrep
         if ((nvvt(irp).eq.0).or.(nvo(irp).eq.0)) goto 30
         done = .false.
         istart = 0
         cco=0

! ** w_efab is constructed over the batches.
! ** ISTART is updated with the actual batch #
! ** available on the local node, NINT is updated with the number of
! ** actually read integrals from this batch and is a pure output parameter !
! ** We therefore have to count upwards in case of ISTART.NE.0 !!

         call w_efab_out_of_core (irp,istart,done,tau,buf2,t1,nbuf3,buf3,nint)

         call putdz (irp,istart,nint,buf3)

! ** now start SORTING and contraction loop
          
         CD = 0
         DO 20 DREP = 1, NREP
           CREP = MULTB(DREP,IRP+NREP,2)
           IF (CREP.LT.DREP) GOTO 20
           DO 15 D = 1, NV(DREP)
              CMIN = 1
              IF (CREP.EQ.DREP) CMIN = D + 1
              DO 10 C = CMIN, NV(CREP)

                 IF(CCO.LT.ISTART) THEN
                   CCO=CCO+1
                   GOTO 10
                 ENDIF

                 CD = CD + 1
                 IF (CD.GT.NINT) THEN
!                   ---------------------------------------
!                   We need the next buffer in this IRREP !
!                   ---------------------------------------
                    ISTART = ISTART + NINT
                    IF(.NOT.DONE) THEN

                    call w_efab_out_of_core (irp,istart,done,tau,buf2,t1,nbuf3,buf3,nint)
                     
                    call putdz (irp,istart,nint,buf3)
                      CD = 1
                      CCO=ISTART
                    ELSE
                      GOTO 30
                    ENDIF
                 ENDIF
                 AKCD = (CD-1)*NVVT(IRP)*RCW+1
                 DO I = 1, NO(DREP)
                    CI = IIVO(CREP,DREP)+(I-1)*NV(CREP)+C
                    AKCI = (IVOVVT(IRP)+(CI-1)*NVVT(IRP))*RCW+1
                    T1DI = (IVO(DREP)+(I-1)*NV(DREP)+D-1) * RCW + 1
                    T1VAL(1) =  T1(T1DI)
                    IF (CARITH) T1VAL(2) =  T1(T1DI+1)
                    CALL XAXPY(NVVT(IRP),T1VAL,BUF3(AKCD),1,w_vvvo(AKCI),1)
                 ENDDO
                 DO I = 1, NO(CREP)
                    DI = IIVO(DREP,CREP)+(I-1)*NV(DREP)+D
                    AKDI = (IVOVVT(IRP)+(DI-1)*NVVT(IRP))*RCW+1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = -T1(T1CI)
                    IF (CARITH) T1VAL(2) = -T1(T1CI+1)
                    CALL XAXPY(NVVT(IRP),T1VAL,BUF3(AKCD),1,w_vvvo(AKDI),1)
                 ENDDO
 10           CONTINUE
 15        CONTINUE
 20      CONTINUE
 30    CONTINUE

   deallocate (tau)

 end subroutine

 subroutine w_ejmb_out_of_core(t1,t2,buf1,buf2,buf3,nbuf3,w_vovo)

#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
#include "ccpar.inc"                                                  

!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
   real(8), intent(inout) :: w_vovo(rcw,nv4)                ! Intermediate.
   real(8), intent(in)    :: buf1(*),buf2(*),buf3(*)        ! buffer arrays
   integer, intent(in)    :: nbuf3 

!---------------Local variables--------------------------------------
!                                                                    
      real*8  :: T1VAL(2)                                                
      logical :: done
      integer :: cd,cco,istart,nint,akcd,ic,akic,t1di,id,akid,t1ci,cmin,c,d,i
      integer :: crep,drep,irp,mint
!------------------------------------------------------------------
! H(AK,CI) = H(AK,CI) - W(AK,CD) * T(D,I)
! H is kept in BUF1, still ordered VOVO
!------------------------------------------------------------------

      t1val = (/0.0d0,0.0d0/)

      if(myproc.eq.master) then
        CALL xcopy(nv4,w_vovo,1,BUF1,1)
      else 
        call xcopy(nv4,a0,0,buf1,1)
      endif

      CALL XTIME(0,1,'--- HINTM: VOVV*T              ')
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         CCO=0
         MINT = NBUF3/NVO(IRP)

! ** reading from GETVOVV: ISTART is updated with the actual batch #
! ** available on the local node, NINT is updated with the number of
! ** actually read integrals from this batch and is a pure output parameter !
! ** We therefore have to count upwards in case of ISTART.NE.0 !!

         CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,MINT)
         CALL DELINT ('VOVV','KDDD',BUF3,IRP,ISTART,NINT)

! ** now start SORTING and contraction loop

         CD = 0
         DO 20 DREP = 1, NREP
           CREP = MULTB(DREP,IRP+NREP,2)
           IF (CREP.LT.DREP) GOTO 20
           DO 15 D = 1, NV(DREP)
              CMIN = 1
              IF (CREP.EQ.DREP) CMIN = D + 1
              DO 10 C = CMIN, NV(CREP)

                 IF(CCO.LT.ISTART) THEN
                   CCO=CCO+1
                   GOTO 10
                 ENDIF

                 CD = CD + 1
                 IF (CD.GT.NINT) THEN
!                   ---------------------------------------
!                   We need the next buffer in this IRREP !
!                   ---------------------------------------
                    ISTART = ISTART + NINT
                    IF(.NOT.DONE) THEN
                      CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,MINT)
                      CALL DELINT ('VOVV','KDDD',BUF3,IRP,ISTART,NINT)
                      CD = 1
                      CCO=ISTART
                    ELSE
                      GOTO 30
                    ENDIF
                 ENDIF
                 AKCD = (CD-1)*NVO(IRP)*RCW+1
                 DO I = 1, NO(DREP)
!                   CI = IIVO(CREP,DREP)+(I-1)*NV(CREP)+C
                    IC = IIOV(CREP,DREP)+(C-1)*NO(DREP)+I
!                   AKCI = (IVOVO(IRP)+(CI-1)*NVO(IRP))*RCW+1
                    AKIC = (IVOOV(IRP)+(IC-1)*NVO(IRP))*RCW+1
                    T1DI = (IVO(DREP)+(I-1)*NV(DREP)+D-1) * RCW + 1
                    T1VAL(1) =  T1(T1DI)
                    IF (CARITH) T1VAL(2) =  T1(T1DI+1)
                    CALL XAXPY(NVO(IRP),T1VAL,BUF3(AKCD),1,BUF1(AKIC),1)
                 ENDDO
                 DO I = 1, NO(CREP)
!                   DI = IIVO(DREP,CREP)+(I-1)*NV(DREP)+D
                    ID = IIOV(DREP,CREP)+(D-1)*NO(CREP)+I
!                   AKDI = (IVOVO(IRP)+(DI-1)*NVO(IRP))*RCW+1
                    AKID = (IVOOV(IRP)+(ID-1)*NVO(IRP))*RCW+1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = -T1(T1CI)
                    IF (CARITH) T1VAL(2) = -T1(T1CI+1)
                    CALL XAXPY(NVO(IRP),T1VAL,BUF3(AKCD),1,BUF1(AKID),1)
                 ENDDO
 10           CONTINUE
 15        CONTINUE
 20      CONTINUE
 30    CONTINUE
!------------------------------------------------------------------
! Now sort H(AK,CI) to H(AI,CK)  is to be done in any case !!
! The H buffer is not changed any more up to the final XGEMM !
! Every other calculation is done in BUF1 !!
!------------------------------------------------------------------
!   CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO, &
!   &            BUF1,w_vovo)

!  call srt1r1 (nrep,multb,lfa,nvo,nv,no,nvo,ivovo,iiov,buf1,buf2)
!  
!  call xaxpy (nv4,-A1,buf2,1,w_vovo,1) 

   CALL XTIME(0,2,'--- HINTM: VOVV*T              ')

 end subroutine

 subroutine f_bar_ae_out_of_core (t1,buf2,nbuf2,hv)

#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "ccpar.inc"

!---------------Calling variables--------------------------------------
     real*8, intent(in)    :: t1(rcw,nfvo),buf2(*)
     real*8, intent(inout) :: hv(rcw*nfvv)
     integer,intent(in)    :: nbuf2 
!---------------Common Blocks--------------------------------------

!---------------Local variables----
     logical :: usedz,right
     integer :: m,n
!----------------------------------
      if(myproc==master) then
        call xcopy(nfvv*rcw,A0,0,buf2,1)
      else 
        call xcopy(nfvv*rcw,A0,0,buf2,1)
      endif
!------------------------------------
! G(A,C) = G(A,C) + W(AK,CD) * T(D,K)
!------------------------------------
      M = MVV(1)
      N = MVO(1)
     IF (N.EQ.0) RETURN
      USEDZ = .FALSE.
      RIGHT = .TRUE.
      CALL SRT20D(NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV, &
     &            JJVO,BUF2,NBUF2,T1,hv,USEDZ,RIGHT)

   end subroutine

   subroutine w_iemn_out_of_core (t1,t2,buf1,buf2,nbuf3,w_ovoo)

!---------------Common Blocks--------------------------------------   
                                                                     
#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
#include "ccpar.inc"                                                  

!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)! Amplitudes.
   real(8), intent(inout) :: w_ovoo(rcw,iovoot(nrep+1))   ! Intermediate.
   real(8), intent(in)    :: buf1(*),buf2(*)              ! buffer arrays.
   integer, intent(in)    :: nbuf3 

!---------------Local variables----------------
      logical :: done,teq
      integer :: kd,istart,nint,m,n,k
      integer :: irp,mint,off1,off2
      real*8,allocatable :: tau(:)
!----------------------------------------------

!----------------------------------------------
!w(ie,mn) = w(ie,mn) + 1/2 V(ie,fg)*tau(fg,mn)
!----------------------------------------------
   call xcopy(iovoot(nrep+1),A0,0,buf2,1)

   allocate(tau(ndimt2*rcw))
   tau = 0.0d0
   call gettau (t1,t2,tau)

    OFF1 = 1
    OFF2 = 1
   DO 130 IRP = 1, NREP
    ISTART = 0
     IF (NVVT(IRP).EQ.0.OR.NVO(IRP).EQ.0) GOTO 130
      M = NVO (IRP)
         N = NOOT(IRP)
          KD = NVVT(IRP)
           MINT = NBUF3/NVO(IRP)

  110      CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF1,MINT)

           K = NINT

 ! ** now we construct the full NVO x NOOT matrix
 ! ** but the contributions are only from the restricted CD range
 ! ** the method is additive since ZGEMM works according to
 ! **      C=A*B + C  ==> all the already computed contributions
 ! ** will be added in this IRREP.

           CALL XGEMM ('N','N',M,N,K,A1,BUF1,M,tau(off1+istart*rcw), &
            &               KD,A1,BUF2(OFF2),M)
           IF (.NOT.DONE) THEN
           ISTART = ISTART + NINT
           GOTO 110
           ENDIF
           OFF1 = OFF1 + KD * N * RCW
           OFF2 = OFF2 + M * N * RCW
 130  CONTINUE

    deallocate (tau)

      call srt1l1 (nrep,multb,lfa,nvo,nv,no,noot,iovoot,iiov,buf2,buf1)

      call xaxpy (iovoot(nrep+1),A1,buf1,1,w_ovoo,1) 

  end subroutine

   subroutine w_efab_out_of_core (irp,istart,done,tau,buf1,t1,nbuf3,w_vvvv,nint)

!-----------------------description------------------------------------------------
!   W_vvvv intermediates are formed such a manner that it could be used as fetching 
!   VVVV integrals, representationwise. However, we will call them only once to construct
!   term2 of W_efam integrals and transfer it to a file, subsequently, inside that routine.
!   batch size is fixed to the distributed batches of VVVV integrals. I assumed that the batch 
!   size for VOVV and VVVV integrals are the same.
!-----------------------------------------------------------------------------------

#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
#include "ccpar.inc"                                                  

!--------------------calling variables-------------
   real(8), intent(in)    :: tau(ndimt2*rcw),t1(ndimt1*rcw)  ! Amplitudes.
   real(8), intent(out)   :: w_vvvv(*)              ! Intermediate.
   real(8), intent(in)    :: buf1(*)                ! buffer arrays
   integer, intent(in)    :: nbuf3,irp,istart
   integer, intent(out)   :: nint
   logical, intent(inout) :: done
!--------------------------------------------------

!--------------------local variables----------------
   integer :: mint,off1,off2,off3,off4,off5
   integer :: abmn,m,n,k,i,irrep,irep,mm(nrep)
   integer :: xxt3vvvt(nrep,nrep),yt3vvvt(nrep),zt3vvvot(nrep+1),zt3vvvvt(nrep+1)
   integer, allocatable :: counter(:), displ(:)
   real*8,allocatable :: w_vvvv_temp_temp(:),w_vvvv_temp(:),w_v(:)  
!---------------------------------------------------

         off1 = 1
         off2 = 1
         off3 = 0

       allocate(counter(0:nmproc-1))
       allocate(displ(0:nmproc-1))

       counter = 0
       displ = 0 

         mint = nbuf3/nvvt(irp)
         call getvvvv (irp,istart,nint,done,w_vvvv,mint)

!-----------------------------------------------------------------------
! generate offsets to distribute the intermediate on different nodes
!-----------------------------------------------------------------------
#if defined (VAR_MPI)
       call interface_mpi_allgather(NINT,1, &
     &   counter,1,global_communicator)

      do i = 0, (nmproc-1)

       counter(i) = counter(i)*noot(irp)            

      enddo

         displ(0) = 0
         do i = 1, (nmproc-1) 
             displ(i) = displ(i-1) + counter(i-1)
         enddo
#endif

!----------------------------------------------
! w(ef,ab) =  V(ef,ab) + V(mn,ab) * TAU(ef,mn)
!          =  V(ef,ab) + TAU(ef,mn) * V'*(ab,mn)            
!----------------------------------------------
       call getoovv (buf1)
       m = nvvt(irp)
       n = nint
       k = noot(irp)
       abmn = (ivvoott(irp)+displ(myproc))*rcw+1
       off1 = off1 + ivvoott(irp)*rcw  

     if (noot(irp) == 0) then
        call xscal (m*n,a1,w_vvvv,1)
     else 
        call xgemm ('n','n',m,n,k,a1,tau(off1),nvvt(irp),buf1(abmn),k, &
      &               a1,w_vvvv,m)
     endif

!-----------------------------------------------------
! w(ef,ab) = w(ef,ab) - V(em,ab)*t(f,m)
!-----------------------------------------------------

      call offset_outofcore_sorter(irp,xxt3vvvt,yt3vvvt,zt3vvvot,zt3vvvvt)

      if (.not.done) off3 = off3+nint

       allocate(w_v(zt3vvvot(nrep+1)*rcw))  

       w_v = 0.0d0

      call srt33 (irp,istart,nrep,multb,nvo,nv,no,nt3vvt, & 
     &            yt3vvvt,zt3vvvot,xxt3vvvt,w_v,nbuf3,off3,done)

       mm(1:nrep) = 0
       do irrep = 1,nrep
         irep = multb(irrep,irp+nrep,2)
         mm(irrep) = mm(irrep)+ nint*nv(irep)
       enddo

       allocate(w_vvvv_temp_temp(dot_product(mm(1:nrep),nv(1:nrep))*rcw))  

       w_vvvv_temp_temp = 0.0d0

       call cntrct ('N','T',mm,nv,no,-a1,w_v,t1,a0,w_vvvv_temp_temp,nrep)

       deallocate(w_v)

       allocate(w_vvvv_temp(nvvt(irp)*nint*rcw)) 

      w_vvvv_temp = 0.0d0 
      call srt43 (irp,istart,nrep,multb,.true.,nvvt,nv,nv,  &
     &    nt3vvt,yt3vvvt,zt3vvvvt,xxt3vvvt,w_vvvv_temp_temp,&
     &    w_vvvv_temp,nbuf3,off3)

       deallocate(w_vvvv_temp_temp)

       call xaxpy(nint*nvvt(irp),a1,w_vvvv_temp,1,w_vvvv,1)

       deallocate(w_vvvv_temp)

   end subroutine

  subroutine G_intm (t2,l2,G_vv,G_oo)

!---------------Common Blocks--------------------------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
!     Calling variables
!     -----------------  
      real(8), intent(in)  :: l2(rcw*ndimt2),t2(rcw*ndimt2)  ! Amplitudes.
      real(8), intent(inout) :: G_vv(:),G_oo(:)

!----------initialization--------------
! G_vv = 0.0d0
! G_oo = 0.0d0
!--------------------------------------

!----------------------------------------------
!    G_ea(e,a)  = - l2(mn,af) * t2(ef,mn)
!----------------------------------------------

  call contraction_442 ((/"p2","p3","o1","o2"/),(/"o1","o2","p1","p3"/),(/"p2","p1"/),G_vv,-1.0d0,1.0d0,nrep,&
                        &  LeftTensor=t2,RightTensor=l2 )

!---------------------------------------------
!     G_im(i,m) = l2(in,ef) * t2(ef,mn)
!---------------------------------------------

  call contraction_442 ((/"o1","o3","p1","p2"/),(/"p1","p2","o2","o3"/),(/"o1","o2"/),G_oo,1.0d0,1.0d0,nrep,&
                        & LeftTensor=l2,RightTensor=t2)

  end subroutine


      subroutine offset_outofcore_sorter(irp,xxt3vvvt,yt3vvvt,zt3vvvot,zt3vvvvt)

#include "symm.inc"
!------------------------------------------------
        integer, intent(in)  :: irp
        integer, intent(out) :: xxt3vvvt(nrep,nrep),yt3vvvt(nrep),zt3vvvot(nrep+1) 
        integer, intent(out) :: zt3vvvvt(nrep+1) 
!------------------------------------------------

!------------------------------------------------
        integer :: irep,krep,ijkrep,jrep
!------------------------------------------------

        YT3VVVT(1:nrep) = 0
        XXT3VVVT(1:nrep,1:nrep) = 0

        DO KREP = 1,NREP
          IJKREP = MULTB(KREP,IRP+NREP,2)
          XXT3VVVT(KREP,IRP)=YT3VVVT(IJKREP)
          YT3VVVT(IJKREP)=YT3VVVT(IJKREP) + NT3VVT(IRP) * NV(KREP)
        ENDDO

      ZT3VVVOT(1)=0
      DO IREP=1,NREP
        ZT3VVVOT(IREP+1) = ZT3VVVOT(IREP) + YT3VVVT(IREP)*NO(IREP)
      ENDDO

      ZT3VVVVT(1)=0
      DO IREP=1,NREP
        ZT3VVVVT(IREP+1) = ZT3VVVVT(IREP) + YT3VVVT(IREP)*NV(IREP)
      ENDDO

      end subroutine


  end module
