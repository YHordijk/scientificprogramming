      SUBROUTINE cceqn_driver_lambda(MAXIT,EPS,FOO,FVO,FVV,L1,S1,L2,S2,
     &                 W_OOOO,CZ,W_VOOV,W_OVVO,W_OVOO,HO,HV,HOV,
     &                 BUF1,BUF2,BUF3,CC1,BB1,BB2,CONV,
     &                 NBUF2,NBUF3)

      use interface_to_mpi
      use xmlout
      use intermediates_1b_2b
      use lambda_equation

      implicit none

C---------------Description--------------------------------------------
C
C     Solve CCSD lambda equations.
C
C
C---------------Calling variables--------------------------------------
C
C     Intermediaries etc.

      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*),L1(*),S1(*),L2(*),S2(*)
      REAL*8 W_OOOO(*),W_VOOV(*),W_OVVO(*),W_OVOO(*),HO(*),HV(*),
     & HOV(*),BUF1(*),BUF2(*),BUF3(*)
      REAL*8 CZ(*) ! is only allocated for Lambda equations
      INTEGER NBUF2,NBUF3,MAXIT

C     for diis ...

      REAL*8 CC1(*)
      REAL*8 BB1(*),BB2(*)

      REAL*8 RMS,CONV,RMS1,ddot,sum1,sum2,zdotc
      REAL*8 DIAG
      INTEGER k,l


C---------------Common Blocks--------------------------------------

#include "files.inc"
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "ihm.inc"
#include "results.inc"
#include "complex.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
      integer ISTAT(df_MPI_STATUS_SIZE)
#endif

C---------------Local variables--------------------------------------

      LOGICAL    :: FOCKSP,PRENE
      real*8,allocatable :: resL1(:),resL2(:),dinv_one(:),dinv_two(:)
      real*8,allocatable :: G_oo(:),G_vv(:),int_wooop(:),w_v_temp(:)
      real*8,allocatable :: l1_temp(:),l2_temp(:)
      integer            :: i !remove
      real*8             :: val  
      integer            :: it, itc, nit,ierr
C---------------Executable code--------------------------------------
C
       buf2(1:nbuf2*rcw) = 0.0d0
       buf1(1:ndimt2*rcw) = 0.0d0


!     Calculate fixed intermediaries for Lambda equations

    5   format(1f10.6)  
                        
        call getampt (s1,s2)

         call f_bar_mi(foo,fvo,s1,s2,ho)

         call w_ijmn(S1,S2,w_oooo)

         call f_bar_me(fvo,S1,hov)

         call wbar_mbej(S2,w_ovvo)

         call w_ijmn(S1,S2,w_oooo)

         call f_bar_me(fvo,S1,hov)

         call f_bar_ae(fvv,fvo,s1,s2,hv)

         call w_mbej(s1,s2,w_voov)

         call w_iemn(s1,s2,w_ovoo,w_ovvo,w_oooo,hov)

      CALL XTIME(5,1,'-- "efam" Intermediate                     ')

        call w_efam(s1,s2,w_ovvo,hov,cz)

        if (myproc.gt.0) then
          cz(1:nv5*rcw) = 0.0d0
        endif

      CALL XTIME(5,2,'-- "efam" Intermediate                     ')

         call w_efam_diagram2 (s1,s2,buf2,buf3,nbuf3,cz)

!   synchronize all the w_efam contributions from different nodes to the master 

#if defined (VAR_MPI)
      if (nmproc .gt. 1) then

!        call xtime(4,1,'-- combining w_efam via mpi_reduce')
        call xcopy (nv5,a0,0,buf3,1)
        ierr = 0

        call interface_mpi_allreduce_r1_work_f77(cz,buf3(1),
     &       rcw*ivovvt(nrep+1),
     &       op_mpi_sum,global_communicator)

        if(ierr.gt.0) then
           call quit('mpi_reduce error in collecting w_efam !')
        endif
        call xcopy(ivovvt(nrep+1),buf3,1,cz,1)
!        call xtime(4,2,'-- combining w_efam via mpi_reduce')

      endif
#endif
        
      allocate (int_wooop(iovoot(nrep+1)*rcw))

         call w_mnie(s1,int_wooop)

      allocate (w_v_temp(nv5*rcw))

         call w_amef(s1,w_v_temp)

! out-of-core section of the intemediates

!      call f_bar_ae_out_of_core (s1,buf2,nbuf2,hv)

!       call w_ejmb_out_of_core (s1,s2,buf1,buf2,buf3,nbuf3,w_voov)

!       call w_iemn_out_of_core (s1,s2,buf1,buf2,nbuf3,w_ovoo)

!       call w_efam_out_of_core (s1,s2,hov,w_voov,w_ovvo,buf1,buf2,
!    &              buf3,nbuf3,nbuf2,CZ)



       CALL TITLER('Solving Left Eigenvector Equation','*',123)

!-------Initialization Step for the Left Eigenvectors-----------

      call srt1c1 (nrep,nv,no,s1,l1) 
      call srt1c1 (nrep,nvvt,noot,s2,l2) 

!     Initialize the work arrays to zero.

       buf1(1:ndimt2*rcw) = 0.0d0

!---------before starting the iterative procedure, we need to initialize the DIIS-------------
      rms = 0.0d0
      itc = 0
      it = -1
      call diisx(ndimt1,ndimt2,L1,L2,it,cc1,bb1,itc,bb2,rms)

      nit = 0

      if (iprnt.ge.1) then
         write(iw,6003)
         write(iw,6001) nit,rms
      endif

!--------------------------------------
!  Start Iterative Procedure
!--------------------------------------

       allocate(resL1(ndimt1*rcw)) 
       allocate(resL2(ndimt2*rcw)) 
       allocate(dinv_one(ndimt1*rcw))
       allocate(dinv_two(ndimt2*rcw))

       call denom_n (eps,dinv_one,dinv_two)

      do nit = 1, maxit 

!     Making non-sequential G_ea, G_im intermediates for lambda equations

      allocate(G_vv(nfvv*rcw))
      allocate(G_oo(nfoo*rcw))

       G_vv = 0.0d0
       G_oo = 0.0d0

      call G_intm (s2,l2,G_vv,G_oo)

!     evaluate L1 equation  

      CALL XTIME(5,1,'- L1 equations                  ')
           
      if (ndimt1 > 0) then

       call xcopy (ndimt1,hov,1,ResL1,1)
       CALL lambda_equation_L1(s1,s2,l1,l2,ho,hv,hov,w_voov,w_ovoo,
     &               cz,w_v_temp,int_wooop,G_oo,G_vv,buf1,buf2,buf3,
     &               nbuf2,nbuf3,resL1,iprnt)

      endif

      CALL XTIME(5,2,'- L1 equations                  ')
                                 
      CALL XTIME(5,1,'-- L2 equations                     ')

       CALL lambda_equation_L2(s1,s2,l1,l2,ho,hv,hov,w_voov,w_oooo,
     &                      int_wooop,w_v_temp,G_oo,G_vv,buf1,buf2,buf3,
     &                      nbuf3,resL2,iprnt,nit)

      CALL XTIME(5,2,'-- L2 equations                     ')

       deallocate (G_oo)
       deallocate (G_vv)


!------------------------------------
!     calculation of RMS value 
!------------------------------------

       rms = dot_product(resL1,resL1) + dot_product(resL2,resL2)

       rms = sqrt((rms)/(ndimt1+ndimt2)*rcw)

      if (iprnt.ge.1) then
       write(iw,*)'maximum value RESL1',maxval(abs(resl1))
       write(iw,*)'maximum value RESL2',maxval(abs(resl2))
      endif
      

!  Make new L1s and L2s by dividing ResL1s and ResL2 by eigenvalues


!     CALL DENOM_L (EPS,RESL1,RESL2,RESL1,RESL2)


!      call xaxpy (ndimt1,-a1,resL1,1,l1,1)
!      call xaxpy (ndimt2,-a2,resL2,1,l2,1)

       call xaxpy (ndimt1,-a1,resL1*dinv_one,1,l1,1)
       call xaxpy (ndimt2,-a2,resL2*dinv_two,1,l2,1)

!!      CALL ZCORE (L1,L2)

      CALL XTIME(5,1,'- DIIS extrapolation            ')
      call diisx(ndimt1,ndimt2,l1,l2,   
     &           it,cc1,bb1,itc,bb2,rms1)
      CALL XTIME(5,2,'- DIIS extrapolation            ')


      CALL PARTS (2,RMS,L1,L2,BUF1,BUF2)


      IF (IPRNT.GE.1) WRITE(IW,6001) NIT,RMS


      allocate(l1_temp(ndimt1*rcw))
      allocate(l2_temp(ndimt2*rcw))

    
      call srt1c1 (nrep,no,nv,l1,l1_temp) 
      call srt1c1 (nrep,noot,nvvt,l2,l2_temp) 

      call gettau(s1,s2,buf2)

      deallocate(l1_temp)
      deallocate(l2_temp)

      CALL FLSHFO(IW)

      if(rms.lt.conv) exit

      enddo  ! end of iteration loop

      write(iw,6079) nit 


      deallocate(dinv_one)
      deallocate(dinv_two)
      deallocate(resL1) 
      deallocate(resL2) 
      deallocate (int_wooop)
      deallocate(w_v_temp)

!     Close and delete DIIS scratch files

      IT = -2
      call diisx(ndimt1,ndimt2,L1,L2,
     &           it,cc1,bb1,itc,bb2,rms)

 111  FORMAT(A20,I4,2f25.15)
 6001 FORMAT(2X,I3,5X,2F20.15,F10.5)
 6002 FORMAT(3X,'MP2   =',F20.15/3X,'EMP2  =',F20.15/)
 6003 FORMAT(//3X,'NIT',10X,'RMS')
 6004 FORMAT(//' WARNING: Requested convergence not reached in CCSD :',
     &         'RMS = ',F20.15)
 6071 FORMAT(//,'  CCSD results',/)
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6074 FORMAT(' CCSD correlation energy :',T40,F25.15)
 6084 FORMAT(' CCD correlation energy :',T40,F25.15)
 6085 FORMAT(' CCS correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6076 FORMAT(' Total CCSD energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 6078 FORMAT(' Convergence :',T40,F25.15)
 6079 FORMAT(' Number of iterations used :',T60,I5)
      RETURN
      END
