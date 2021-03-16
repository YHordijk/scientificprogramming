  module lambda_equation

!this module consists of subroutines for lambda-equations at different level of theory.

  use interface_to_mpi
  use contraction
  use intermediates_1b_2b
  use wrapper_sorting

    public lambda_equation_mp2
    public lambda_equation_L2
    public lambda_equation_L1
    public lambda2_disconnected


     private

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface


     contains


      subroutine lambda_equation_mp2 (eps,fvo,t1,s1,t2,s2)

      implicit none

!---------------Description--------------------------------------------
!
!     Solves MP2 amplitude or lambda equations.
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

!     Intermediaries etc.

      REAL*8 EPS(*),FVO(*),T1(*),S1(*),T2(*),S2(*)

!---------------Common Blocks--------------------------------------

#include "files.inc"
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "results.inc"
#include "complex.inc"

!---------------Local variables--------------------------------------

      REAL*8 DIAG

!---------------Executable code--------------------------------------

!  INITIALIZE T2 AND T1 AMPLITUDES

      CALL XCOPY (NDIMT1,FVO,1,T1,1)
      CALL GETVVOO (T2)
      CALL DENOM (EPS,T1,T2,T1,T2)
      CALL ZCORE (T1,T2)


!  INITIALIZE LAMDA AMPLITUDES

         CALL XCOPY (NDIMT1,T1,1,S1,1)
         CALL XCOPY (NDIMT2,T2,1,S2,1)
 6071 FORMAT(//,'  MP2 results',/)
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 7000 FORMAT(//' Timing of routine MP2EQN :'      &
     &/' Before MP2EQN :',T30,F12.3,' seconds'    &
     &/' Energy, T1-Diag :',T30,F12.3,' seconds'  &
     &/' Untimed parts :',T30,F12.3,' seconds'    &
     &/' Total time in MP2EQN :',T30,F12.3,' seconds')
      end subroutine lambda_equation_mp2


      subroutine lambda_equation_L2 (t1,t2,l1,l2,ho,hv,hov, &
     &             w_voov,w_oooo,int_wooop,w_v_temp,G_oo,G_vv,buf1,buf2, &
     &             buf3,nbuf3,resL2,iprnt,niter)

      use contraction

      implicit none

!---------------Description--------------------------------------------
!     Evaluate residue for lambda equation.
!----------------------------------------------------------------------

!---intermediates required-------------
!f_bar_ae,f_bar_mi,f_bar_me,w_ejmb,w_ijmn,G_ea,G_im
!--------------------------------------

!---------------Common Blocks--------------------------------------
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "ccpar.inc"
#include "files.inc"

!       calling variables
!       -----------------

      integer,intent(in) :: nbuf3,iprnt,niter
      real(8),intent(in) :: ho(rcw*nfoo),hv(rcw*nfvv),hov(nfvo*rcw)  ! one body intermediates
      real(8),intent(in) :: w_voov(ivoov(nrep+1)*rcw),w_oooo(nv1*rcw)! two body intermediates
      real(8),intent(in) :: w_v_temp(:),int_wooop(:)                 ! two body intermediates
      real(8),intent(in) :: G_vv(nfvv*rcw),G_oo(nfoo*rcw)            ! G-intermediates
      real(8),intent(in) :: t1(ndimt1*rcw),t2(ndimt2*rcw)! converged T1 amplitude
      real(8),intent(in) :: L1(ndimt1*rcw),L2(ndimt2*rcw)! L1 and L2 amplitudes from last iteration.
      real(8),intent(in) :: buf1(*),buf2(*),buf3(*)      ! buffer arrays
      real*8,intent(out) :: ResL2(ndimt2*rcw)            ! residue for L2 amplitude


!     Local variables
!     ---------------

      logical          :: done,teq
      real*8           :: t1val(2)
      real*8,allocatable :: int_vv(:)
      real*8,allocatable :: int_oo(:),w_vvvv(:),w_intm(:)
      real*8,allocatable :: L2_temp(:),ResL2_temp(:),tau(:)
      integer, allocatable :: counter(:), displ(:)
      integer     :: jklrep,m2,irep,jrep,abrep,jkloff,ij,ijkl,ijkl1,mint
      integer     :: off2,i,j,istart,m,k,n,irp,ntot,nint,kstart,off1
      integer     :: jabrep,x(nrep,nrep),y(nrep)
      integer     :: off3,nd,l, totsize
      

      type(indices) ::e,f,g 

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface

!---------------Executable code--------------------------------------

    5   format(10f10.6)

       call getoovv (ResL2)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL2 after initialization',dot_product(ResL2,ResL2)
!-------------------------------------------------
!   term1 : ResL2 = ResL2 + P(ab) L2 (ij,ae)*f_bar_ae(e,b)
!-------------------------------------------------

       call contraction_424 ((/"o1","o2","p1","p3"/),(/"p3","p2"/), &
     &    (/"o1","o2","p1","p2"/),ResL2,1.0d0,1.0d0,nrep,LeftTensor=L2, &
     &    RightTensor=hv)

       if (iprnt.ge.1) &

     & write(iw,*)'RESL2 after term1',dot_product(ResL2,ResL2)
!----------------------------------------------
!   term2 : ResL2 = ResL2 - P(ij) L2 (im,ab)*f_bar_im (j,m)
!----------------------------------------------

       call contraction_424 ((/"o1","o3","p1","p2"/),(/"o2","o3"/), &
     &    (/"o1","o2","p1","p2"/),ResL2,-1.0d0,1.0d0,nrep, &
     &    LeftTensor=L2,RightTensor=ho)

       if (iprnt.ge.1) &

     & write(iw,*)'RESL2 after term2',dot_product(ResL2,ResL2)
!--------------------------------------------------------
!   term3 : ResL2 = ResL2 + L2(mn,ab) * W_ijmn(ij,mn)
!--------------------------------------------------------

      call contraction_444((/"o1","o2","o3","o4"/),(/"o3","o4","p1","p2"/), &
     &    (/"o1","o2","p1","p2"/),ResL2, &
     &    1.0d0,1.0d0,nrep,RightTensor=L2,LeftTensor=w_oooo)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL2 after term3', dot_product(ResL2,ResL2)

!--------------------------------------------------------
!   term4 : ResL2 = ResL2 + P(ij)P(ab) L2(im,ae) * W_bmej(je,bm)
!--------------------------------------------------------

      call contraction_444((/"o1","o3","p1","p3"/),(/"o2","p3","p2","o3"/), (/"o1","o2","p1","p2"/),ResL2, &
     &    1.0d0,0.0d0,nrep,righttensor=w_voov,LeftTensor=L2)

!-----------------------------------------------------------
!   term5 : ResL2 = ResL2 + P(ab) V(ij,ae)*G_vv(e,b)
!-----------------------------------------------------------

      call contraction_424 ((/"o1","o2","p1","p3"/),(/"p3","p2"/), &
     &        (/"o1","o2","p1","p2"/),ResL2,1.0d0,1.0d0,nrep,RightTensor=G_vv)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL2 after term5',dot_product(ResL2,ResL2)

!---------------------------------------------------------------------------
!     term6: L2(ij,ab) = L2(ij,ab) + L1(m,a) * w_mnie(ij,mb)
!---------------------------------------------------------------------------

      call contraction_424 ((/"o1","o2","o3","p2"/),(/"o3","p1"/), &
     &   (/"o1","o2","p2","p1"/),ResL2,-1.0d0,-1.0d0,nrep,RightTensor=L1, &
     &   LeftTensor=int_wooop)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL2 after term6',dot_product(ResL2,ResL2)
!--------------------------------------------------------------------
!    term7 : ResL2 = ResL2 - P(ij) V(im,ab)*G_oo(j,m)
!--------------------------------------------------------------------

      call contraction_424 ((/"o1","o3","p1","p2"/),(/"o2","o3"/), &
     &        (/"o1","o2","p1","p2"/),ResL2,-1.0d0, &
     &        1.0d0,nrep,RightTensor=G_oo)

       if (iprnt.ge.1) &
      & write(iw,*)'RESL2 after term7',dot_product(ResL2,ResL2)

!-------------------------------------------------------------------------
!   term 8 : ResL2 = ResL2 + P(ij)L1(i,e)*W(ej,ab)
!-------------------------------------------------------------------------

     call contraction_244 ((/"o1","p3"/),(/"p3","o2","p1","p2"/), &
    &  (/"o1","o2","p1","p2"/),ResL2,1.0d0,1.0d0,nrep,LeftTensor=L1, &
    &   RightTensor=w_v_temp)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL2 after term8',dot_product(ResL2,ResL2)
!-------------------------------------------------------------------------
!   term 9 : ResL2 = ResL2 + P(ij)P(ab) L1(i,a)*f_bar_me(j,b)
!-------------------------------------------------------------------------

      call lambda2_disconnected(l1,hov,ResL2)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL2 after term9',dot_product(ResL2,ResL2)

!--------------------------------------------------------
!    term 10 : ResL2 = ResL2 + L2(ij,ef) * W_vvvv(ef,ab)
!--------------------------------------------------------

       allocate(counter(0:nmproc-1))
       allocate(displ(0:nmproc-1))

       counter = 0
       displ = 0 

      done = .false.
      allocate(tau(ndimt2*rcw))
      call gettau(t1,t2,tau)
      off2 = 1
      off1 = 1
      do irp = 1, nrep
         istart = 0
         totsize = 0
         if ((nvvt(irp).eq.0).or.(noot(irp).eq.0)) cycle
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

#if defined (VAR_MPI)
       call interface_mpi_allgather(NINT,1, &
     &   counter,1,global_communicator)

      do i = 0, nmproc-1

       counter(i) = counter(i)*noot(irp)*rcw            

      enddo

         displ(0) = 0
         do i = 1, nmproc-1 
             displ(i) = displ(i-1) + counter(i-1)
         enddo

       allocate(ResL2_temp(counter(myproc)*rcw))

       ResL2_temp = 0.0d0    

      call xgemm ('n','n',m,k,n,A1,l2(off1),m, &
     &           w_vvvv,n,a1,ResL2_temp,m)

      call interface_mpi_allgatherv(ResL2_temp,counter(myproc), &
     &    buf2(off2+totsize),counter,displ,global_communicator)

       deallocate(ResL2_temp)
#else
      call xgemm ('n','n',m,k,n,A1,l2(off1),m, &
     &           w_vvvv,n,a1,ResL2(off2+istart*m*rcw),m)
#endif
         totsize = totsize + sum(counter)
         deallocate(w_vvvv) 
         if (.not.done) then
            istart = istart + nint
            goto 10
         endif
         off2 = off2 + m * n * rcw
         off1 = off1 + m * n * rcw
       enddo

       deallocate(tau)

#if defined (VAR_MPI)
          call xaxpy (ndimt2,A1,buf2,1,ResL2,1)
#endif

       if (iprnt.ge.1) &
     & write(iw,*)'RESL2 after term10',dot_product(ResL2,ResL2)

      end subroutine

      subroutine lambda2_disconnected(l1,f_bar,ResL2)

!---------------Description--------------------------------------------

!     Calculating the contribution from the disconnected piece in lambda2 residue

!---------------Common Blocks----------------------------------------

      implicit none
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------

      real*8, intent(in)    :: l1(ndimt1*rcw),f_bar(nfvo*rcw)
      real*8, intent(inout) :: ResL2(ndimt2*rcw)

!---------------Local variables--------------------------------------

      integer :: ijaboff,abrep,arep,brep,a,b,amin,imin,i,j,ij,jb,ia,ja, &
     &           ib,ni,irep,jrep
      real*8  :: T1MIN(2)
!---------------Executable code--------------------------------------

      t1min=(/0.0d0,0.0d0/)

      ijaboff = 0
      DO ABREP = 1, NREP
      DO BREP = 1, NREP
      AREP = MULTB(BREP,ABREP+NREP,2)
      IF (AREP.LT.BREP) CYCLE
      DO B = 1, NV(BREP)
         AMIN = 1
         IF (AREP.EQ.BREP) AMIN = B + 1
         DO A = AMIN, NV(AREP)
            DO JREP = 1, NREP
            IREP = MULTB(JREP,ABREP+NREP,2)
            IF (IREP.LT.JREP) CYCLE
!--------------------------------------------------------------
! Resl2(ij,ab) = Resl2(ij,ab) + l1(i,a) * f_bar(j,b) + f_bar(i,a) * l1(j,b)
!--------------------------------------------------------------
            IF (JREP.EQ.BREP) THEN
               IJ = 1
               JB = (IVO(BREP) + (B-1) * NO(BREP)) * RCW + 1
               DO J = 1, NO(JREP)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = J + 1
                  IA = (IVO(AREP) + (A-1) * NO(AREP) + IMIN - 1)*RCW + 1
                  NI = NO(IREP) - IMIN + 1
                  if (NI > 0) then
                    CALL XAXPY (NI,f_bar(jb),L1(IA),1,RESL2(IJABOFF+IJ),1)
                    CALL XAXPY (NI,L1(JB),f_bar(IA),1,RESL2(IJABOFF+IJ),1)
                  endif
                  IJ = IJ + NI * RCW
                  JB = JB + RCW
               ENDDO
            ENDIF
!--------------------------------------------------------------
! Resl2(ij,ab) = Resl2(ij,ab) - L1(j,a) * f_bar(i,b) - f_bar(j,a) * L1(i,b)
!--------------------------------------------------------------
            IF (JREP.EQ.AREP) THEN
               IJ = 1
               JA = (IVO(AREP) + (A-1) * NO(AREP)) * RCW + 1
               DO J = 1, NO(JREP)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = J + 1
                  IB = (IVO(BREP) + (B-1) * NO(BREP) + IMIN - 1)*RCW + 1
                  NI = NO(IREP) - IMIN + 1
                  T1MIN(1) = -f_bar(JA)
                  IF (CARITH) T1MIN(2) = -f_bar(JA+1)
                  if (NI > 0) then
                  CALL XAXPY (NI,T1MIN,L1(IB),1,ResL2(IJABOFF+IJ),1)
                   T1MIN(1) = -L1(JA)
                   IF (CARITH) T1MIN(2) = -L1(JA+1)
                  CALL XAXPY (NI,T1MIN,f_bar(IB),1,ResL2(IJABOFF+IJ),1)
                  endif
                  IJ = IJ + NI * RCW
                  JA = JA + RCW
               ENDDO
            ENDIF
!--------------------------------------------------------
! UPDATE OFFSET AND GO TO NEXT IRREP PAIR
!--------------------------------------------------------
            IF (IREP.NE.JREP) THEN
               IJABOFF = IJABOFF + NO(IREP) * NO(JREP) * RCW
            ELSE
               IJABOFF = IJABOFF + NO(IREP) * (NO(IREP)-1) * RCW / 2
            ENDIF
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE

      subroutine symm_index_outofcore_sorter(istart,x)

#include "symm.inc"
!------------------------------------------------
        integer, intent(in) :: istart
        integer, intent(out) :: x(nrep,nrep)
!------------------------------------------------

!------------------------------------------------
        integer :: y(nrep), krep, ijrep, ijkrep
!------------------------------------------------

            y = 0

          do ijrep = 1,nrep
            do krep = 1,nrep

            ijkrep = multb(krep,ijrep+nrep,2)
            x(krep,ijrep) = y(ijkrep)
            y(ijkrep) =  y(ijkrep) + (istart+1)*no(krep)

            enddo
         enddo

      end subroutine


      subroutine lambda_equation_L1(t1,t2,l1,l2,ho,hv,hov,w_voov,w_ovoo, &
     &               w_vvvo,w_v_temp,int_wooop,G_oo,G_vv,buf1,buf2,buf3, &
     &               nbuf2,nbuf3,ResL1,iprnt)

       implicit none

!       Description
!       -----------
!     Evaluate       Lambda1 equations
!     Here we have to take care of the fact that certain integral classes
!     are now split over individual nodes. The corresponding contributions
!     which are linear in the integrals are added at the end of the T1
!     routines. Therefore we have to take care of calculating the
!     non-distributed contributions to the T1 amplitudes only on the
!     master node in order to avoid double-counting in the PARTS routine !
!     all direct VOVV/VVVV contributions and VOVV/VVVV contributions stemming
!     from intermediates have to be calculated on EACH node.

!     intermediates required for L1 equation
!     ---------------------------------------
!     f_bar_mi,f_bar_me,f_bar_ae,w_iemn,w_ejmb,G_ea,G_im


!   common blocks
!   -------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "ccpar.inc"
#include "files.inc"

!       calling variables
!       -----------------

      integer,intent(in) :: nbuf2,nbuf3,iprnt
      real(8),intent(in) :: ho(nfoo*rcw),hv(nfvv*rcw),hov(nfvo*rcw)    ! one body intermediates
      real(8),intent(in) :: w_voov(ivoov(nrep+1)*rcw),w_v_temp(:), &
     &                      w_ovoo(iovoot(nrep+1)*rcw),w_vvvo(nv5*rcw)&! two body intermediates
     &                      ,int_wooop(:)
      real(8),intent(in) :: G_vv(nfvv*rcw), G_oo(nfoo*rcw)             ! G-intermediates
      real(8),intent(in) :: t1(ndimt1*rcw),t2(ndimt2*rcw)              ! converged T1 amplitude
      real(8),intent(in) :: L1(ndimt1*rcw),L2(ndimt2*rcw)              ! L1 and L2 amplitudes from last iteration.
      real(8),intent(in) :: buf1(*),buf2(*),buf3(*)                    ! buffer arrays
      real*8,intent(inout) :: ResL1(ndimt1*rcw)                        ! residue for L1 amplitude


!     Local variables
!     ---------------
      logical            :: done,teq,usedz,right
      real*8,allocatable :: int_temp(:), l2_sorted1(:)
      real*8,allocatable :: w_o_temp(:)
      real*8             :: ddot
      integer            :: i,j,ij
!-------------------
! Initialize residue
!-------------------

    5   format(10f10.6)

      call xcopy (ndimt1,hov,1,ResL1,1)

     if (iprnt.ge.1) &
     &   write(iw,*)'RESL1: init',dot_product(resL1,resL1)
!---------------------------------
!    ResL1 = ResL1 + L(i,e) * f_bar_ae(e,a)
!---------------------------------

       call contraction_222 ((/"o1","p2"/),(/"p2","p1"/),(/"o1","p1"/), &
     &    L1,hv,ResL1,1.0d0,1.0d0,nrep)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL1: after term1',dot_product(resL1,resL1)
!----------------------------------------
!    ResL1 = ResL1 - L(m,a)*f_bar_mi(i,m)
!----------------------------------------

        call contraction_222 ((/"o2","p1"/),(/"o1","o2"/),(/"o1","p1"/), &
     &    L1,ho,ResL1,-1.0d0,1.0d0,nrep)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL1: after term2',dot_product(resL1,resL1)
!----------------------------------------------
!      ResL1 = ResL1 - G(m,n) * W(mi,na)
!      w_mina(mi,na) = w_mina(mi,na) + V(mi,fa)*t(f,n)
!----------------------------------------------

      call contraction_242 ((/"o3","o2"/),(/"o2","o1","o3","p1"/), &
     & (/"o1","p1"/),G_oo,ResL1,-1.0d0,1.0d0,nrep,RightTensor=int_wooop)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL1: after term3',dot_product(resL1,resL1)
!-------------------------------------------------
!      ResL1 = ResL1 + L1(m,e) * W(ie,am)
!-------------------------------------------------

      call contraction_242 ((/"o2","p2"/),(/"o1","p2","p1","o2"/), &
     &   (/"o1","p1"/),L1,ResL1,1.0d0,1.0d0,nrep,RightTensor=w_voov)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL1 after term4',dot_product(resL1,resL1)
!-----------------------------------------------------
!       ResL1 = ResL1 -  L2(mn,ae)*W(ie,mn)
!-----------------------------------------------------

      call contraction_442((/"o2","o3","p1","p2"/),(/"o1","p2","o2","o3" &
     &          /),(/"o1","p1"/),ResL1,-1.0d0,1.0d0,nrep,LeftTensor=L2,  &
     &          RightTensor=w_ovoo)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL1 after term5',dot_product(resL1,resL1)

!-------------------------------------------------------
!      ResL1    = ResL1 + L2(im,ef) * W(ef,am)
!-------------------------------------------------------

      call contraction_442((/"o1","o2","p2","p3"/),(/"p2","p3","p1","o2"/),(/"o1","p1"/),ResL1,1.0d0,1.0d0,nrep, &
       &         LeftTensor=L2,RightTensor=w_vvvo)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL1 after term6',dot_product(resL1,resL1)

!-------------------------------------------------------
!    ResL1 = ResL1 - G(e,f) * V(ei,fa)
!-------------------------------------------------------

      call contraction_242 ((/"p3","p2"/),(/"p2","o1","p3","p1"/), &
     &  (/"o1","p1"/),G_vv,ResL1,-1.0d0,1.0d0,nrep,RightTensor=w_v_temp)

       if (iprnt.ge.1) &
     & write(iw,*)'RESL1 after term7',dot_product(resL1,resL1)

      end subroutine

     end module
