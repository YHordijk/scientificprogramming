!      Copyright (c) 2019 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCDENSC(IOPT,T1,S1,T2,S2,
     &                   BUF1,BUF2,BUF3,DOO,DVV,DVO,DOV,DMO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Menu driven routine. Depending on value of IOPT
C
C     Calculate unperturbed density matrix at the
C       - SCF level (0)
C       - MP2 level (1)
C       - CCSD level (2)
C       - CCSD(T) level (3)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER IOPT,nps
      REAL*8 T1(*),S1(*),T2(*),S2(*)  !! T1,T2 are CCSD amplitudes and S1,S2 are Lambdas
      REAL*8 DOO(*),DVV(*),DVO(*),DOV(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
      REAL*8 DMO(*)
      REAL*8 ddot
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (NFOO,A0,0,DOO,1)
      CALL XCOPY (NFVV,A0,0,DVV,1)
      CALL XCOPY (NFVO,A0,0,DVO,1)
      CALL XCOPY (NFVO,A0,0,DOV,1)
      IF (IOPT.LT.1) RETURN
C----------------------------------------------------------------------
C D(I,J) = D(I,J) - T(CD,IL) * L(CD,JL) = D(I,J) - T(CDL,I) * L(CDL,J)
C----------------------------------------------------------------------
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,T2,BUF1)
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,S2,BUF2)
      CALL CNTRCT ('C','N',NO,NO,NVVOT,-A1,BUF1,BUF2,A1,DOO,NREP)
C----------------------------------------------------------------------
C D(I,J) = D(I,J) - T(C,I) * L(C,J) 
C----------------------------------------------------------------------
      CALL CNTRCT ('C','N',NO,NO,NV,-A1,T1,S1,A1,DOO,NREP)
C----------------------------------------------------------------------
C D(A,B) = D(A,B) + L(AD,KL) * T(BD,KL) = H(A,C) + L(A,DKL) * T(B,DKL)
C----------------------------------------------------------------------
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF1)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,S2,BUF2)
      CALL CNTRCT ('N','C',NV,NV,NVOOT,A1,BUF1,BUF2,A1,DVV,NREP)
C----------------------------------------------------------------------
C D(A,B) = D(A,B) + L(A,K) * T(B,K)
C----------------------------------------------------------------------
      CALL CNTRCT ('N','C',NV,NV,NO,A1,T1,S1,A1,DVV,NREP)
C
      IF (CARITH) THEN
         CALL CONJUGA (NFOO,DOO,1)
         CALL CONJUGA (NFVV,DVV,1)
      ENDIF
C
C###########################################
C convert the density matrix to DIRAC format
C###########################################
      
      IF (IOPT.LT.2) CALL CONVDM(BUF1,BUF2,DOO,DVV,DVO,DOV,DMO)
C
      IF (IOPT.LT.2) RETURN
C----------------------------------------------------------------------
C DVO(A,I) = L(A,I)
C----------------------------------------------------------------------
!!      CALL XCOPY (NDIMT1,S1,1,DVO,1)         

      call srt1c1n(nrep,nv,no,S1,DVO)
C----------------------------------------------------------------------
C DOV(A,I) = T(A,I)
C----------------------------------------------------------------------
      CALL XCOPY (NDIMT1,T1,1,DOV,1)         
C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) + T(AC,IK) * L(C,K)
C----------------------------------------------------------------------
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,MVO,JVOVO,JJVO,JJVO,
     &              T2,BUF1)
      CALL XGEMV ('N',MVO(1),MVO(1),A1,BUF1,MVO(1),S1,1,A1,DOV,1)

C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) - T(A,K) * L(C,K) * T(C,I)
C----------------------------------------------------------------------
      CALL CNTRCT ('C','N',NO,NO,NV,-A1,S1,T1,A0,BUF1,NREP)
      CALL CNTRCT ('N','N',NV,NO,NO,A1,T1,BUF1,A1,DOV,NREP)

C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) - T(A,K) * L(CD,LK) * T(CD,LI)
C----------------------------------------------------------------------
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,S2,BUF1)
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,T2,BUF2)
      CALL CNTRCT ('C','N',NO,NO,NVVOT,-A1,BUF1,BUF2,A0,BUF3,NREP)
      CALL CNTRCT ('N','N',NV,NO,NO,A1,T1,BUF3,A1,DOV,NREP)

C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) - L(CD,KL) * T(AD,KL) * T(D,I)
C----------------------------------------------------------------------
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,S2,BUF1)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF2)
      CALL CNTRCT ('N','C',NV,NV,NVOOT,-A1,BUF1,BUF2,A0,BUF3,NREP)
      CALL CNTRCT ('C','N',NV,NO,NV,A1,BUF3,T1,A1,DOV,NREP)

      IF (CARITH) THEN
         CALL CONJUGA (NFVO,DOV,1)
         CALL CONJUGA (NFVO,DVO,1)
      ENDIF

      CALL CONVDM(BUF1,BUF2,DOO,DVV,DVO,DOV,DMO)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

       subroutine cc_density (l1,l2,t1,t2,doo,dvv,buf1,buf2,dmo)

       use contraction

#include "symm.inc"
#include "complex.inc"
!-------------------------calling variables-------------------------
       real*8, intent(in)  :: l1(ndimt1*rcw),l2(ndimt2*rcw),
     &                        t1(ndimt1*rcw),t2(ndimt2*rcw)
       real*8, intent(out) :: doo(nfoo*rcw),dvv(nfvv*rcw),dmo(*) !norbt*
!    &                        norbt*nz) 
       real*8, intent(out) :: buf1(*),buf2(*)
!-------------------------------------------------------------------
       real*8, allocatable :: dvo(:),dov(:)
       real*8, allocatable :: woo_temp(:),wpp_temp(:)
       integer :: iprnt=2, dm_type
!-------------------------------------------------------------------

        if (iprnt.ge.2) write(*,*) 'at start of cc_density'
        if (iprnt.ge.2) write(*,*) 'T1',dot_product(t1,t1)
        if (iprnt.ge.2) write(*,*) 'T2',dot_product(t2,t2)
        if (iprnt.ge.2) write(*,*) 'L1',dot_product(l1,l1)
        if (iprnt.ge.2) write(*,*) 'L2',dot_product(l2,l2)

        doo(1:nfoo*rcw) = 0.0d0
        dvv(1:nfvv*rcw) = 0.0d0

!-------------------------------------------------------------------
!   doo(j,i) = doo(j,i) - l2(jm,ef) * t2(ef,im)
!-------------------------------------------------------------------

        call contraction_442 ((/"o2","o3","p1","p2"/),
     &    (/"p1","p2","o1","o3"/),(/"o2","o1"/),doo,-1.0d0,1.0d0,nrep,
     &    LeftTensor=l2,RightTensor=t2)
        if (iprnt.ge.2) write(*,*)'diagram1a',dot_product(doo,doo)

!---------------------------------------------------------------------
!    doo(j,i) = doo(j,i) - l1(j,e)*t1(e,i)
!---------------------------------------------------------------------

       call contraction_222 ((/"o2","p1"/),(/"p1","o1"/),(/"o2","o1"/),
     &                         l1,t1,doo,-1.0d0,1.0d0,nrep)
        if (iprnt.ge.2) write(*,*)'diagram1b',dot_product(doo,doo)

!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l2(mn,ae)*t2(be,mn) 
!----------------------------------------------------------------------

      call contraction_442 ((/"o1","o2","p1","p3"/),
     &    (/"p2","p3","o1","o2"/),(/"p2","p1"/),dvv,1.0d0,1.0d0,nrep,
     &     LeftTensor=l2,RightTensor=t2 )
        if (iprnt.ge.2) write(*,*)'diagram2a',dot_product(dvv,dvv)


!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l1(m,a)*t1(b,m)
!----------------------------------------------------------------------

      call contraction_222 ((/"o1","p1"/),(/"p2","o1"/),(/"p2","p1"/),
     &              l1,t1,dvv,1.0d0,1.0d0,nrep)
        if (iprnt.ge.2) write(*,*)'diagram2b',dot_product(dvv,dvv)


       allocate(dvo(nfvo*rcw))

!----------------------------------------------------------------------------
!     dvo(a,i) = t1(a,i)
!----------------------------------------------------------------------------

      call xcopy (nfvo,t1,1,dvo,1)
        if (iprnt.ge.2) write(*,*)'diagram3a',dot_product(dvo,dvo)

!--------------------------------------------------------------------------
!     dvo(a,i) = dvo(a,i) + l1(m,e) * t(ae,im)
!--------------------------------------------------------------------------

!     call contraction_242 ((/"o2","p2"/),(/"p1","p2","o1","o2"/), 
!    & (/"p1","o1"/),l1,dvo,1.0d0,1.0d0,nrep,RightTensor=t2)


      call contraction_422 ((/"p1","p2","o1","o2"/),(/"o2","p2"/), 
     & (/"p1","o1"/),l1,dvo,1.0d0,1.0d0,nrep,LeftTensor=t2)
        if (iprnt.ge.2) write(*,*)'diagram3b',dot_product(dvo,dvo)


!----------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l1(m,e) * t(e,i)*t(a,m)
!----------------------------------------------------------

      allocate(woo_temp(nfoo*rcw))

      woo_temp = 0.0d0

      call contraction_222 ((/"o2","p2"/),(/"p2","o1"/),(/"o2","o1"/),
     &              l1,t1,woo_temp,1.0d0,1.0d0,nrep)
        if (iprnt.ge.2) write(*,*)'diagram3c',dot_product(dvo,dvo)


      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),
     &              woo_temp,t1,dvo,-1.0d0,1.0d0,nrep)
        if (iprnt.ge.2) write(*,*)'diagram3d',dot_product(dvo,dvo)

      deallocate(woo_temp)

!-------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) -  l2(mn,ef)*t2(ef,in)*t1(a,m) 
!-------------------------------------------------------------

      allocate(woo_temp(nfoo*rcw))

      woo_temp = 0.0d0

      call contraction_442((/"o2","o3","p2","p3"/),(/"p2","p3","o1","o3"
     &       /),(/"o2","o1"/),woo_temp,1.0d0,1.0d0,nrep,LeftTensor=L2, 
     &          RightTensor=t2)
        if (iprnt.ge.2) write(*,*)'diagram3e',dot_product(dvo,dvo)

      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),
     &              woo_temp,t1,dvo,-1.0d0,1.0d0,nrep)
        if (iprnt.ge.2) write(*,*)'diagram3f',dot_product(dvo,dvo)

      deallocate(woo_temp)

!-----------------------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l2(mn,ef)*t1(e,i)*t2(af,mn)
!-----------------------------------------------------------------------

      allocate(wpp_temp(nfvv*rcw))

      wpp_temp = 0.0d0

      call contraction_442((/"o2","o3","p2","p3"/),(/"p1","p3","o2","o3"
     &       /),(/"p1","p2"/),wpp_temp,-1.0d0,1.0d0,nrep,LeftTensor=L2, 
     &          RightTensor=t2)
        if (iprnt.ge.2) write(*,*)'diagram3g',dot_product(dvo,dvo)

      call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/),
     &              wpp_temp,t1,dvo,1.0d0,1.0d0,nrep)
        if (iprnt.ge.2) write(*,*)'diagram3h',dot_product(dvo,dvo)

      deallocate(wpp_temp)

!----------------------------------------------------------
!       dov(i,a) = l1(i,a)
!----------------------------------------------------------

       allocate(dov(nfvo*rcw))

      call xcopy (nfvo,l1,1,dov,1)
        if (iprnt.ge.2) write(*,*)'diagram4a',dot_product(dov,dov)

      if (iprnt.ge.9) then

         print *,'>>>>> elements, DOV not symmetrized'
         dm_type = 3
         call print_density_matrix(dov,dm_type)
         write(*,*)'dov before c.c.',dot_product(dov,dov)
         print *,'<<<<< elements, DOV not symmetrized'

         print *,'>>>>> elements, DVO not symmetrized'
         dm_type = 4
         call print_density_matrix(dvo,dm_type)
         write(*,*)'dvo before c.c.',dot_product(dvo,dvo)
         print *,'<<<<< elements, DVO not symmetrized'

         print *,'>>>>> elements, DVV not symmetrized'
         dm_type = 2
         call print_density_matrix(dvv,dm_type)
         write(*,*)'dvv before c.c.',dot_product(dvv,dvv)
         print *,'<<<<< elements, DVV not symmetrized'

         print *,'>>>>> elements, DOO not symmetrized'
         dm_type = 1
         call print_density_matrix(doo,dm_type)
         write(*,*)'doo before c.c.',dot_product(doo,doo)
        print *,'<<<<< elements, DOO not symmetrized'

      end if

      if (carith) then
         call conjuga (nfoo,doo,1)
         call conjuga (nfvv,dvv,1)
         call conjuga (nfvo,dvo,1)
         call conjuga (nfvo,dov,1)
      endif

      CALL CONVDM(BUF1,BUF2,DOO,DVV,DOV,DVO,DMO)

       deallocate(dov)
       deallocate(dvo)

      end subroutine


      subroutine print_density_matrix(dm,block_type)
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "dgroup.h"
         real(kind=8), intent(in) :: dm(*)
         integer, intent(in) :: block_type ! 1 : oo, 2: vv, 3: ov, 4: vo
         integer :: pq, irep, p, q

         pq = 0
         if (block_type.eq.1) then
            DO IREP = 1, NREP
               write (*,*) 'DOO for irrep ',irep
               DO p = 1, NO(IREP)
                  DO q = 1, NO(IREP)
                     pq = pq + 1
                     IF (CARITH) THEN
       write (*,'(3X,I4,2x,2F28.10,2I4)') pq,dm(rcw*pq-1),dm(rcw*pq),p,q
                     ELSE
                        write (*,'(3X,2I4,2x,F28.10)') pq,dm(pq),p,q
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         else if (block_type.eq.2) then
            DO IREP = 1, NREP
               write (*,*) 'DVV for irrep ',irep
               DO p = 1, NV(IREP)
                  DO q = 1, NV(IREP)
                     pq = pq + 1
                     IF (CARITH) THEN
       write (*,'(3X,I4,2x,2F28.10,2I4)') pq,dm(rcw*pq-1),dm(rcw*pq),p,q
                     ELSE
                        write (*,'(3X,2I4,2x,F28.10)') pq,dm(pq),p,q
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         else if (block_type.eq.3) then
            DO IREP = 1, NREP
               write (*,*) 'DOV for irrep ',irep
               DO p = 1, NO(IREP)
                  DO q = 1, NV(IREP)
                     pq = pq + 1
                     IF (CARITH) THEN
       write (*,'(3X,I4,2x,2F28.10,2I4)') pq,dm(rcw*pq-1),dm(rcw*pq),p,q
                     ELSE
                        write (*,'(3X,2I4,2x,F28.10)') pq,dm(pq),p,q
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         else if (block_type.eq.4) then
            DO IREP = 1, NREP
               write (*,*) 'DVO for irrep ',irep
               DO p = 1, NV(IREP)
                  DO q = 1, NO(IREP)
                     pq = pq + 1
                     IF (CARITH) THEN
       write (*,'(3X,I4,2x,2F28.10,2I4)') pq,dm(rcw*pq-1),dm(rcw*pq),p,q
                     ELSE
                        write (*,'(3X,2I4,2x,F28.10)') pq,dm(pq),p,q
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         else
            call quit('unknown density matrix block')
         end if

      end subroutine

      SUBROUTINE CONVDM(BUF1,BUF2,DOO,DVV,DVO,DOV,dmo)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Convert the density matrix to DIRAC format
C
C
C---------------Routines called----------------------------------------
C
C     DCOPY
C     DOOTODIR
C     DVVTODIR
C---------------Last modified------------------------------------------
C
C     Author : Joost van Stralen - March 2002 : Extended by Avijit Shee
C     for CC
C
C---------------Calling variables--------------------------------------
C
      REAL*8 DOO(*),DVV(*),DOV(*),DVO(*)
      REAL*8 BUF1(*),BUF2(*)
      REAL*8 DMO(*)
      REAL*8 sym_factor !we will change this factor depending on the
                        !fact whether our desired property opeartor is
                        !hermition or anti-hermition.   
C
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "dgroup.h"
C
C---------------Local variables--------------------------------------
C
       real*8, allocatable :: DOO_TRANS(:),DVV_TRANS(:),D_TRANS(:)
      integer a, ab, ai,b,i,ia,ij,irep,j 
! debug
      logical :: debug=.false.
      integer :: dm_type
C---------------Executable code--------------------------------------
C
C
C     Reorder density matrix
C
C######################################################################
C  STEP 1:
C
C   OO       VO         OV      VV            Re          Im
C +-----+   +-----+   +------+   +-----+    +----------+----------+
C |     |   |     |   |      |   |     |    |          |          |
C |     | , |     | , |      | , |     | -> |          |          |
C +-----+   +-----+   +------+   +-----+    |          |          |
C                                           |          |          |
C                                           +----------+----------+
C                                           <---O+V---->
C                                            unbar + bar
C######################################################################
C

      sym_factor = 1.0d0 !1.0d0 for the hermition part of density matrix.
                         !it will be -1.0d0 for the anti-hermition part. 

      ALLOCATE(D_TRANS(NFVO*RCW))

      if (debug) then
         write(*,*)'dov in symmet. before', 
     &           dot_product(dov(1:NFVO*rcw),dov(1:NFVO*rcw))
      end if 

      CALL SRT1C1 (NREP,NO,NV,DVO,D_TRANS)

      CALL XAXPY (NFVO,sym_factor*A1,D_TRANS,1,DOV,1)

      if (debug) then
         write(*,*)'dov in symmet. after', 
     &           dot_product(dov(1:NFVO*rcw),dov(1:NFVO*rcw))
         write(*,*)'dov_trans in symmet.', 
     &           dot_product(d_trans(1:NFVO*rcw),d_trans(1:NFVO*rcw))
      end if

      CALL XSCAL (NFVO,AP5,DOV,1)

      if (debug) then
         print *,'>>>>> elements, DOV symmetrized'
         dm_type = 3
         call print_density_matrix(dov,dm_type)
         write(*,*)'dov in symmet. after, scaled by 0.5', 
     &              dot_product(dov(1:NFVO*rcw),dov(1:NFVO*rcw))
         print *,'<<<<< elements, DOV symmetrized'
      end if

      DEALLOCATE(D_TRANS)

      CALL DZERO(DMO,NORBT*NORBT*NZ)

      buf1(1:nfvo*rcw) = 0.0d0

      AI = 0
      DO IREP = 1, NREP
         DO I = 1, NO(IREP)
            DO A = 1, NV(IREP)
               AI = AI + 1
               IF (CARITH) THEN
                  BUF1(AI) = DOV(RCW*AI-1)
                  BUF2(AI) = -DOV(RCW*AI)
               ELSE
                  BUF1(AI)= DOV(AI)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF (CARITH) THEN
         CALL DCOPY(NFVO,BUF2,1,BUF1(NFVO+1),1)
      ENDIF

      CALL RELTODIR(DMO,BUF1,NFVO,RCW,NREP,NV,NO)

C##############
C     DOO
C##############

!!! symmetrize the whole array.


      buf1(1:nfoo*rcw) = 0.0d0

      ALLOCATE(DOO_TRANS(NFOO*RCW))

      if (debug) then
         write(*,*)'doo in symmet. before', 
     &           dot_product(doo(1:NFOO*rcw),doo(1:NFOO*rcw))
      end if

      CALL SRT1C1 (NREP,NO,NO,DOO,DOO_TRANS)

      CALL XAXPY (NFOO,sym_factor*A1,DOO_TRANS,1,DOO,1)

      if (debug) then
         write(*,*)'doo in symmet. after', 
     &           dot_product(doo(1:NFOO*rcw),doo(1:NFOO*rcw))
         write(*,*)'doo_trans in symmet.', 
     &          dot_product(doo_trans(1:NFOO*rcw),doo_trans(1:NFOO*rcw))
      end if

      DEALLOCATE(DOO_TRANS)

      CALL XSCAL (NFOO,AP5,DOO,1)

      if (debug) then
         print *,'>>>>> elements, DOO symmetrized'
         dm_type = 1
         call print_density_matrix(doo,dm_type)
         write(*,*)'doo in symmet. after, scaled by 0.5', 
     &           dot_product(doo(1:NFOO*rcw),doo(1:NFOO*rcw))
         print *,'<<<<< elements, DOO symmetrized'
      end if

      IJ = 0
      DO IREP = 1, NREP
         DO J = 1, NO(IREP)
            DO I = 1, NO(IREP)
               IJ = IJ + 1
               IF (CARITH) THEN
                  BUF1(IJ) = DOO(RCW*IJ-1)
                  BUF2(IJ) = -DOO(RCW*IJ)
               ELSE
                  BUF1(IJ)= DOO(IJ)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

C
C Put imaginary part after real part
C
      IF (CARITH) THEN
         CALL DCOPY(NFOO,BUF2,1,BUF1(NFOO+1),1)
      ENDIF

      CALL DOOTODIR(dmo,BUF1,NO,nfoo,NREP,RCW)

C##############
C     DVV
C##############

      buf1(1:nfvv*rcw) = 0.0d0

      ALLOCATE(DVV_TRANS(NFVV*RCW))
      CALL SRT1C1 (NREP,NV,NV,DVV,DVV_TRANS)
      CALL XAXPY (NFVV,sym_factor*A1,DVV_TRANS,1,DVV,1)

      if (debug) then
         write(*,*)'dvv in symmet. after',
     &           dot_product(dvv(1:NFVV*rcw),dvv(1:NFVV*rcw))
         write(*,*)'dvv_trans in symmet.',
     &          dot_product(dvv_trans(1:NFVV*rcw),dvv_trans(1:NFVV*rcw))
      end if

      DEALLOCATE(DVV_TRANS)

      CALL XSCAL (NFVV,AP5,DVV,1)

      if (debug) then
         print *,'>>>>> elements, DVV not symmetrized'
         dm_type = 2
         call print_density_matrix(dvv,dm_type)
         write(*,*)'dvv in symmet. after, scaled by 0.5',
     &           dot_product(dvv(1:NFVV*rcw),dvv(1:NFVV*rcw))
         print *,'<<<<< elements, DVV not symmetrized'
      end if

      AB = 0
      DO IREP = 1, NREP
         DO B = 1, NV(IREP)
            DO A = 1, NV(IREP)
               AB = AB + 1
               IF (CARITH) THEN
                  BUF1(AB) = DVV(RCW*AB-1)
                  BUF2(AB) = -DVV(RCW*AB)
               ELSE
                  BUF1(AB)= DVV(AB)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C Put imaginary part after real part
C
      IF (CARITH) THEN
         CALL DCOPY(NFVV,BUF2,1,BUF1(NFVV+1),1)
      ELSE
      ENDIF

      CALL DVVTODIR(dmo,BUF1,NV,nfvv,NREP,RCW)

C##############
C     DVO
C##############

      buf1(1:nfvo*rcw) = 0.0d0

      CALL SRT1C1 (NREP,NV,NO,DOV,DVO)

      CALL XSCAL (NFVO,sym_factor*A1,DVO,1)

      if (debug) then
         print *,'>>>>> elements, DVO symmetrized'
         dm_type = 4
         call print_density_matrix(dvo,dm_type)
         write(*,*)'dvo :',
     &      dot_product(dvo(1:nfvo*rcw),dvo(1:nfvo*rcw))
         print *,'<<<<< elements, DVO symmetrized'
      end if

      IA = 0
      DO IREP = 1, NREP
         DO A = 1, NV(IREP)
            DO I = 1, NO(IREP)
               IA = IA + 1
               IF (CARITH) THEN
                  BUF1(IA) = DVO(RCW*IA-1)
                  BUF2(IA) = -DVO(RCW*IA)
               ELSE
                  BUF1(IA)= DVO(IA)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF (CARITH) THEN
         CALL DCOPY(NFVO,BUF2,1,BUF1(NFVO+1),1)
      ENDIF

      CALL DOVTODIR(DMO,BUF1,NFVO,RCW,NREP,NO,NV)

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE CCDENSZ(IOPT,NSP,IREPSPI,EPS,T1,S1,T2,S2,
     &                   BUF1,BUF2,BUF3,NBUF2,NBUF3,CBUF,
     &                   XVO,XOV,DOO,DVV,DVO,DOV,
     &                   EVO,EOV,D2VO,D2OV,DG,CC1,BB1,BB2,CONV)
C
      use interface_to_mpi
      implicit none
C
C---------------Description--------------------------------------------
C
C     Menu driven routine. Depending on value of IOPT
C
C     Solve Z-vector equations at the
C       - SCF level (0)
C       - MP2 level (1)
C       - CCSD level (2)
C       - CCSD(T) level (3)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER IOPT,NSP,IREPSPI,NBUF2,NBUF3,CC1,BB1,BB2,CONV
      REAL*8 T1(*),S1(*),T2(*),S2(*)
      REAL*8 DOO(*),DVO(*),DOV(*),DVV(*)
      REAL*8 D2VO(*),D2OV(*),DG(*)
      REAL*8 XVO(*),XOV(*)
      REAL*8 EVO(*),EOV(*)
      REAL*8 EPS(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*),CBUF(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*8 NAME
      COMPLEX*16 EPHASE(MAXOP)
      LOGICAL BDONE,USEDZ,RIGHT
      integer dum,eoff,eoper,jprnt,k,m,maxit,n
C
C---------------Executable code--------------------------------------
C
      IF (IOPT.LT.1) RETURN
C
C     No relaxation if there are no virtuals within the right irrep
C
      IF (NFVO.EQ.0) RETURN
C
C     Calculate X-intermediate
C
      CALL XCOPY (NDIMX,A0,0,XVO,1)
C--------------------------------------------
C X(A,I) = X(A,I) + V(A,K,C,D) * TAU(C,D,I,K)
C--------------------------------------------
      CALL XCOPY (NDIMT2,S2,1,BUF1,1)
      CALL T1EQN2 (BUF1,BUF3,NBUF3,XVO)
C
C------------------------------------
C X(A,I) = X(A,I) + W(CI,DA) * D(C,D)
C------------------------------------
      M = NFVV
      N = NFVO
      USEDZ = .FALSE.
      RIGHT = .FALSE.
      CALL SRT20D(NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
     &            JJVO,BUF2,NBUF2,XVO,DVV,USEDZ,RIGHT)
C
C In parallel runs we add all distributed contribution together here and
C continue on the master.
C
#if defined (VAR_MPI)
      IF(NMPROC.GT.1) THEN
C
        CALL XCOPY(NDIMX,XVO,1,BUF1,1)
        CALL DZERO(XVO,RCW*NDIMX)
        call interface_MPI_allREDUCE_r1_work_f77(BUF1,XVO,RCW*NDIMX,
     &                                           op_mpi_sum,
     &                                           global_communicator)
      ENDIF
#endif
C--------------------------------------------
C X(A,I) = X(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------
      CALL XCOPY (NDIMT2,S2,1,BUF3,1)
      CALL T1EQN1 (BUF3,XVO)
C---------------------------------------------------------------------
C X(A,I) = X(A,I) + V(IK,AL) * D(K,L)
C---------------------------------------------------------------------
      CALL GETVOOO (BUF1)
      CALL SRT1ST4 (NREP,MULTB,LTR,LFA,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &              BUF1,BUF2)
      M = NFVO
      K = NFOO
      CALL XGEMV ('N',M,K,A1,BUF2,M,DOO,1,A1,XVO,1)
C
      CALL XSCAL(NDIMX,-A1,XVO,1)
      CALL XCOPY (NDIMX,XVO,1,XOV,1)
      IF (CARITH) CALL CONJUGA (NDIMX,XOV,1)
C
C     Solve Z-vector equations by calling the RPA routine
C
      NAME = 'MP2-GRAD'
      DO EOPER = 1, NEOPER
         EOFF = (EOPER-1)*NDIMX*RCW+1
         CALL RDPROP(IPRNT,NAMEE(EOPER),NSP,LFA,CBUF,IREPSPI,
     &               EPHASE(EOPER),EVO(EOFF),EOV(EOFF),DUM,DUM)
      ENDDO
C
C     Silence RPA routine when the print level is smaller than one
C
      IF (IPRNT.GE.1) THEN
         JPRNT = IPRNT
      ELSE
         JPRNT = -1
      ENDIF
C
      BDONE = .FALSE.
      MAXIT = MXITGR
      CALL RPAA (MAXIT,NEOPER,NAMEE,NAME,EPHASE,A1,BDONE,EPS,A0,
     &           EVO,EOV,XVO,XOV,DVO,D2VO,DOV,D2OV,DG,BUF1,BUF2,
     &           CC1,BB1,BB2,CONV,JPRNT,DEBUG,TIMING)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MP2ZVEC(T2,BUF1,BUF2,BUF3,NBUF3,XMO,DMO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Solve Z-vector equations at the MP2 level, using the AO based 
C     implementation.
C
C     Called from: CCFOPR/ccdriv.F
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
CMI 
C---------------Calling variables--------------------------------------
C
      REAL*8 T2(*)
      REAL*8 DMO(*)
      REAL*8 XMO(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
      INTEGER NBUF3
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      integer kfree, lfree
C
C---------------Executable code--------------------------------------
C
C     +------------------------------------------------+
C     | Generate the T's in the right order for the    |
C     | contraction with the 3/4 transformed integrals |
C     +------------------------------------------------+
C

      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              T2,BUF1)

      CALL SRTDIR (NREP,MULTB,MVO,J2VOVO,JJVO,JJVO,BUF1,BUF2)
      CALL SRTT1T2(BUF2,BUF1)


C
C     +---------------------------------------------------+
C     | Generate the tau dependent part of the Lagrangian |
C     +---------------------------------------------------+
C
      LFREE=NBUF3*RCW
      CALL TR5INI(XMO,BUF2,BUF1,BUF3,LFREE)

C
C     +-------------------------------------------------------+
C     | Calculate the ji' and b'a parts of the density matrix |
C     +-------------------------------------------------------+
C
      KFREE=1
      LFREE=NBUF3*RCW
      CALL CALCDJI(DMO,XMO,BUF3,KFREE,LFREE)
C
C     +-------------------------------------------------+
C     | Generate the D dependent part of the Lagrangian |
C     +-------------------------------------------------+
C
      KFREE=1 
      LFREE=NBUF3*RCW
      CALL MPGLVC(DMO,XMO,BUF3,KFREE,LFREE,1)

C
C     +----------------------------------------+
C     | Solve Z-vector equations by using XRPA |
C     +----------------------------------------+
C
      KFREE=1
      LFREE=NBUF3*RCW
      CALL PRPXLR_MPG(DMO,XMO,BUF3,KFREE,LFREE)

C 
      CALL FLSHFO(IW)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CONVXM(BUF1,BUF2,XVO,XMO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Convert x-matrix(Lagrangian) to DIRAC format
C
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Joost van Stralen - April 2002
C
C---------------Calling variables--------------------------------------
C
      REAL*8 XVO(*)
      REAL*8 BUF1(*),BUF2(*)
      REAL*8 XMO(*)
C
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "dgroup.h"
C
C---------------Local variables--------------------------------------
C
       integer a, ai, i, irep
C
C---------------Executable code--------------------------------------
C
C
C###########################################
C STEP 1: Put imaginary part after real part
C###########################################
      AI = 0
      DO IREP = 1,NREP
         DO I = 1, NO(IREP)
            DO A = 1, NV(IREP)
               AI = AI + 1
               IF (CARITH) THEN
                  BUF1(AI) = XVO(AI*RCW-1)
                  BUF2(AI) = XVO(AI*RCW)
c                 print*,'real XVO',ai,BUF1(AI)
c                 print*,'imag XVO',ai,BUF2(AI)
               ELSE
                  BUF1(AI) = XVO(AI)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
      IF (CARITH) THEN
         CALL DCOPY(NFVO,BUF2,1,BUF1(NFVO+1),1)
      ELSE
      ENDIF
C######################################################################
C  STEP 2:
C
C           Re      Im
C     -  +-------+-------+
C   V |  |*******|*******|
C  u+b|  |*******|*******|
C     -  +-------+-------+
C        <---O--->
C         unbar+bar
C                     -------------->
C
C          Re        I         J       K
C       +--------+--------+--------+--------+
C     P |        |        |        |        |
C     O |        |        |        |        |
C     V |  ***   |  ***   |  ***   |  ***   |
C       +--------+--------+--------+--------+
C        P  O  V
C######################################################################
C
      CALL DZERO(XMO,NORBT*NORBT*NZ)
      CALL RELTODIR(XMO,BUF1,NFVO,RCW,NREP,NV,NO)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXPVAL (DOO,DVV,DVO,DOV,EOO,EVV,EVO,EOV,EPHASE,CEXPA)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates expectation value of an operator.
C
C---------------Routines called----------------------------------------
C
C     BLAS routines
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 DOO(*),DVV(*),DVO(*),DOV(*)
      REAL*8 EOO(*),EVV(*),EVO(*),EOV(*)
      COMPLEX*16 CEXPA,EPHASE
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 EXPA,EXPA_SUM
      REAL*8 DDOT,ZDOT
      integer ii,irep
C
C---------------Executable code--------------------------------------
C
      EXPA = (0.0, 0.0)
      CALL XDOTU (EXPA,NFVO,EOV,1,DVO,1)
      EXPA_SUM = EXPA
      CALL XDOTU (EXPA,NFVO,EVO,1,DOV,1)
      EXPA_SUM = EXPA_SUM + EXPA
c     print*,' Relaxation contribution ',expa_sum
      CALL XDOTU (EXPA,NFOO,EOO,1,DOO,1)
c     print*,' OO contribution ',EXPA
      EXPA_SUM = EXPA_SUM + EXPA
      CALL XDOTU (EXPA,NFVV,EVV,1,DVV,1)
c     print*,' VV contribution ',EXPA
      EXPA_SUM = EXPA_SUM + EXPA
C
C     This is the correlation contribution, add the reference value
C
      II = 1
      DO IREP = 1, NREP
         CALL XDOTU (EXPA,NO(IREP),EOO(II),NO(IREP)+1,A1,0)
         EXPA_SUM = EXPA_SUM + EXPA
         II = II + NO(IREP)*NO(IREP)*RCW
      ENDDO
C
      CEXPA = EPHASE * EXPA_SUM
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE  MP2_NATORB(DMO,NOPRNT_input,NOZG)
C***********************************************************************
C
C     Calculate MP2 nat. orb. and nat. orb. occ. num.
C
C     Called from: CCFOPR/ccdriv.F
C
C     INPUT:   DMO - second order density matrix
C
C     Written by Stefan Knecht and Hans Joergen Aa. Jensen, Feb. 2008
C
C     Last modifications: 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbmp2no.h"
#include "infpar.h"
C
      DIMENSION DMO(NORBT,NORBT,NZ)
      LOGICAL NOZG
      real(8), allocatable :: CMO_NAT(:), OCC_NAT(:)
      integer, allocatable :: ibeig(:)
C
C---------------Local variables--------------------------------------
C
      CHARACTER*8 LABEL
C
      CALL QENTER('MP2_NATORB')
!define DEBUG_MP2_NO
#ifdef DEBUG_MP2_NO
      noprnt = max(noprnt_input,10)
#else
      noprnt = noprnt_input
#endif
C
C     memory allocation
C     -----------------
      allocate(CMO_NAT(NASHT*NASHT*NZ))
      allocate(OCC_NAT(NESHT))
      allocate(ibeig(norbt))
C
C     initialize
C
      CALL DZERO(CMO_NAT,NASHT*NASHT*NZ)
      CALL DZERO(OCC_NAT,NESHT)
      ibeig(1:norbt) = 0

!     !> read ibeig (info array on boson irreps) from DFCOEF
      if( linear .or. spinfr )then
        call opnfil(lucoef,'DFCOEF','OLD','RIBEIG ')
        call reacmo(lucoef,'DFCOEF',dummy,dummy,ibeig,dummy,8)
        close(lucoef,STATUS='KEEP')
      end if
C
C     build diagonal according to Fock matrix
C     ---------------------------------------
      CALL ADD_TWO_FOCK(DMO,NOZG,NOPRNT)
C
C     calculate NO occ. num. 
C     ----------------------
      CALL MP2_NATORB_GEN(DMO,CMO_NAT,OCC_NAT,ibeig,NOPRNT)
C
C     save nat. orb. occ. num. and transform. matrix to file 
C     ------------------------------------------------------
      IF( MYTID .eq. MPARID )THEN
C
C       ... KRMCOLD file
C
C       sk + hjaaj: we use KRMCOLD and not KRMCSCF for the MP2 nat.orb.
C       coefficients, because the KRMC module starts by deleting the
C       KRMCSCF file !!! (And it uses CMO from KRMCOLD as first priority
C       for starting orbitals.)
C
        CALL OPNFIL(LUKRMC,'KRMCOLD','UNKNOWN','MP2_NA')
        CALL NEWLAB('*KRMCMP2',LUKRMC,LUPRI)
C
        CALL WRTKRMC(LUKRMC, 'MP2NATOB',CMO_NAT,NASHT*NASHT*NZ)
        CALL WRTKRMC(LUKRMC, 'MP2NOOCC',OCC_NAT,NESHT)
        CALL IWRTKRMC(LUKRMC,'IBEIGORI',ibeig,NORBT)
        CLOSE(LUKRMC,STATUS='KEEP')
C
      END IF
C
C     memory deallocation
C     -------------------
      deallocate(OCC_NAT,CMO_NAT,ibeig)
C
      CALL QEXIT('MP2_NATORB')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MP2_NATORB_GEN(DMO,CMO_NAT,OCC_NAT,ibeig,NOPRNT)
C***********************************************************************
C
C     diagonalize second order density matrix DMO
C     and calculate nat. orb. occ. num.
C
C     INPUT:   DMO - second order density matrix
C
C     Written by Stefan Knecht and Hans Joergen Aa. Jensen, Feb. 2008
C
C     Last modifications: 
C
C***********************************************************************
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "maxorb.h"
#include "maxash.h"
#include "dcbidx.h"
#include "../moltra/dcbtr3.h"
C
      DIMENSION OCC_NAT(*), CMO_NAT(NASHT,NASHT,NZ), ibeig(norbt)
      DIMENSION DMO(NORBT,NORBT,*)
      real(8), allocatable :: scratch(:)
      real(8), allocatable :: dv(:,:,:), occ(:)
      real(8), allocatable :: cmo_nat_all(:)
      integer :: lwork = 0
      real(8), allocatable :: work(:)

      call legacy_lwork_get(LWORK)
#include "memint.h"
!      Allocating just the amount we need to avoid memory problems
! miro: guys, you miscalculted the exact amount, this artificial increase by 100 is needed to pass tests ;)
      lwork = norbt + nz*nasht**2 + 100
      call alloc(WORK,LWORK,id='WORK in MP2_NATORB')
      KFRSAV = KFREE

      IF ( NOPRNT .ge. 5 )THEN
         CALL HEADER(
     &               'MP2_NATORB_GEN: DMO :',-1)
         CALL PRQMAT(DMO,
     &               NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF

      allocate(occ(nasht),dv(nasht,nasht,nz),cmo_nat_all(norbt**2*nz))
      dv          = 0
      occ         = 0
      cmo_nat_all = 0
      iashi       = 0
      !> set index array idxu2g
      do i = 1,nfsym
        do j = 1, nash(i)
          iashi         = iashi  + 1
          IDXU2G(IASHI) = iorb(i)+ npsh(i) + niocc(i) + j
        end do
      end do

      
      !> find max j_z value
      max_loop = 1
      if(linear) max_loop = iimax(norbt,ibeig,1)

!     extract a-a part of dmo

      do indx_loop = 1, max_loop

        is_symm = 0
        if(linear)then
          !> get mj-value
          is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
        end if

        CALL MATGAT_symm(dmo,norbt,norbt,dv,nasht,
     &                   nasht,idxu2g,ibeig,is_symm,nz)
      end do

      !> diagonalize
      CALL DIAG_piv2fock('A',dv,occ,cmo_nat,
     &                   idxu2g,ibeig,work,lwork)

!     insert new matrix block into transformation matrix
      do indx_loop = 1, max_loop

        is_symm = 0
        if(linear)then
          !> get mj-value
          is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
        end if
        CALL MATSCT_symm(cmo_nat_all,NORBT,NORBT,cmo_nat,
     &                   NASHT,NASHT,IDXU2G,ibeig,is_symm,NZ)
      end do

      call dealloc(work)
C
C     debug print section
C
      IF ( NOPRNT .ge. 5 ) THEN
         DO IFSYM = 1, NFSYM
           CALL HEADER(
     &     'MP2_NATORB_GEN: Eigenvectors :',-1)
           WRITE(LUPRI,*) ' fermionic irrep', IFSYM
           CALL PRQMAT(CMO_NAT(1+iash(ifsym),1+iash(ifsym),1),
     &                 nasht,nash(IFSYM),nasht,nash(IFSYM),
     &                 NZ,IPQTOQ(1,0),LUPRI)
         END DO
      ENDIF
C
C     MP2 - nat. orb. occ. num.
C
      ioff  = 1
      ioff2 = 1
      do i = 1, nfsym
        occ_nat(ioff:(ioff+niocc(i)-1)) = 2.0d0
        occ_nat((ioff+niocc(i)):(ioff+niocc(i)+nash(i)-1)) = 
     &       occ(ioff2:(ioff2+nash(i)-1))
        ioff  = ioff  + nesh(i) 
        ioff2 = ioff2 + nash(i) 
      end do
C
      IF ( NOPRNT .ge. 3 ) THEN
        IOFF2 = 1
        DO I = 1, NFSYM
           NOCCI = NESH(I)
           IF (NOCCI .EQ. 0) THEN
              WRITE(LUPRI,9001) I
           ELSE
              OCCSUM = DSUM(NOCCI,OCC_NAT(IOFF2),1)
              WRITE(LUPRI,9002) I, (OCC_NAT(IOFF2-1+J),J=1,NESH(I))
              WRITE(LUPRI,9003) OCCSUM
           END IF
           IOFF2 = IOFF2 + NESH(I)
        END DO
      END IF

!     release memory
      deallocate(dv,occ)

C
 9001 FORMAT(/'   Symmetry',I3,//,'   No occupied orbitals')
 9002 FORMAT(/'   Natural orbital occupation numbers, symmetry',
     &            I3,//,(5F14.9))
 9003 FORMAT(/'   Sum =',T15,F14.9)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADD_TWO_FOCK(DMO,NOZG,NOPRNT)
C***********************************************************************
C
C     build proper diagonal part in second order density matrix DMO
C
C     INPUT:   DMO - second order density matrix
C
C     Written by Stefan Knecht, Feb. 2008
C
C     Last modifications: 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "../moltra/dcbtr3.h"
C
      DIMENSION DMO(NORBT,NORBT,*)
      LOGICAL NOZG
C
C     add "alpha" and "beta" by scaling with 2.0
C
      CALL DSCAL(NZ*NORBT*NORBT,2.0D0,DMO,1)
C
      IF ( NOPRNT .ge. 10 )THEN
         CALL HEADER(
     &               'ADD_TWO_FOCK input: DMO :',-1)
         CALL PRQMAT(DMO,
     &               NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      IOFF2 = 1
C
      DO IFSYM = 1, NFSYM  
        IOFF1 = IORB(IFSYM) + NPSH(IFSYM) + 1
C
!       WRITE(LUPRI,*) ' NDMOQC(2,IFSYM,1)',NDMOQC(2,IFSYM,1), nozg
C       zero 'wrong' offdiagonal terms for frozen orbs
        IF( .NOT. NOZG )THEN
          DO IZ = 1, NZ
            CALL DZERO(DMO(IOFF1,IOFF1,IZ),NDMOQC(2,IFSYM,1)*NORBT)
          END DO
        END IF
C
        IOCC_A = 0
!       print *, 'NOCC(IFSYM)',NOCC(IFSYM)
        DO I = 1, NISH_DHF(IFSYM)
           DMO(IOFF1+IOCC_A,IOFF1+IOCC_A,1) = 
     &     DMO(IOFF1+IOCC_A,IOFF1+IOCC_A,1) + 2.0D0
           IOCC_A = IOCC_A + 1
        END DO
        DO I = NISH_DHF(IFSYM)+1,NISH_DHF(IFSYM)+NASH_DHF(IFSYM)
           DMO(IOFF1+IOCC_A,IOFF1+IOCC_A,1) = 
     &     DMO(IOFF1+IOCC_A,IOFF1+IOCC_A,1) + 1.0D0
       ! Setting occupation of open shell orbitals to 1.0, no matter actual value.
       ! This will put them in the right place when nat.orb. are sorted.
           IOCC_A = IOCC_A + 1
        END DO
      END DO
C
      IF ( NOPRNT .ge. 10 )THEN
         CALL HEADER(
     &               'ADD_TWO_FOCK output: DMO :',-1)
         CALL PRQMAT(DMO,
     &               NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_CC_DENSITY (WF,BASIS,DMAT,IERR)
C
C---------------Description--------------------------------------------
C
C     Retrieve relaxed density matrix for processing by other modules
C     of DIRAC (currenty only first-order properties)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Common Blocks--------------------------------------
C
#include "implicit.h"
#include "files.inc"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
C
C---------------Calling variables--------------------------------------
C
      PARAMETER (HALF=0.5D0)
      CHARACTER*4 WF                     ! Wave function type
      CHARACTER*2 BASIS                  ! Representation
      REAL*8 DMAT(*)                     ! Density matrix
      INTEGER IERR                       ! Error code
C
C---------------Local variables--------------------------------------
C
      REAL*8, ALLOCATABLE:: DMO(:,:,:)   ! Density matrix in MO basis
      LOGICAL FND,FNDLAB,ADD_2_FOCK
      CHARACTER*8 LABEL
      CHARACTER*6 FILE_NAME
C
C---------------Executable code--------------------------------------

      ADD_2_FOCK = .FALSE.
      SELECT CASE (WF)

         CASE ("CCSD")
            WRITE (LUPRI,*) " Retrieving CC density matrix from file"
            LABEL        = 'CCSDDENS'
            FILE_NAME    = 'CCDENS'
            ADD_2_FOCK   = .TRUE.    
            IERR         = 0
         CASE ("MP2 ")
            WRITE (LUPRI,*) " Retrieving MP2 density matrix from file"
            LABEL        = 'MP2 DENS'
            FILE_NAME    = 'CCDENS'
            ADD_2_FOCK   = .TRUE.
            IERR         = 0
         CASE ("EOM ")
            WRITE (LUPRI,*) " Retrieving EOMCC density matrix from file"
            LABEL        = 'EOM DENS'
            FILE_NAME    = 'CCDENS'
            ADD_2_FOCK   = .FALSE.
            IERR         = 0
         CASE ("MC  ")
            WRITE (LUPRI,*) " Retrieving MCSCF density matrix from file"
            LABEL        = 'MC  DENS'
            FILE_NAME    = 'MCDENS'
            ADD_2_FOCK   = .FALSE.
            IERR         = 0
            CALL QUIT('MCSCF density read not fully implemented yet!')
         CASE ("CI  ")
            WRITE (LUPRI,*) " Retrieving CI density matrix from file"
            LABEL        = 'CI  DENS'
            FILE_NAME    = 'MCDENS'
            ADD_2_FOCK   = .FALSE.
            IERR         = 0
            CALL QUIT('CI density read not fully implemented yet!')
         CASE ("HF  ")
            WRITE (LUPRI,*) " Generating Hartree-Fock density matrix"
            FILE_NAME    = 'CCDENS'
            ADD_2_FOCK   = .TRUE.
            IERR         = 0
         CASE DEFAULT
            WRITE (LUPRI,*) " Unknown or unimplemented density: ",WF
            WRITE (LUPRI,*) " Will use Hartree-Fock density matrix"
            IERR = 1
      END SELECT

C     Read the appropriate correlated density from file

      ALLOCATE (DMO(NORBT,NORBT,NZ),STAT=IERR)
      IF (IERR.NE.0) THEN
        CALL QUIT('Error in DMO allocation !')
      ENDIF
      CALL OPNFIL(MCCDENS,FILE_NAME,'OLD','GET_CC')
      FND =  FNDLAB(LABEL,MCCDENS)
      CALL READT(MCCDENS,NORBT*NORBT*NZ,DMO)
      CLOSE(Unit=MCCDENS,Status='KEEP')
!      CLOSE(Unit=MCCDENS,Status='DELETE')

      
C     Add the closed shell Hartree-Fock contribution
      IF(ADD_2_FOCK)THEN
C       The add_two routine scales by 2, scale back to the orginal afterwards
        CALL ADD_TWO_FOCK(DMO,.TRUE.,0)
        CALL DSCAL(NZ*NORBT*NORBT,HALF,DMO,1)
      END IF

      SELECT CASE (BASIS)
         CASE ("MO")
            CALL DCOPY (N2ORBTQ,DMO,1,DMAT,1)
         CASE ("AO")
            CALL DENSMO_TO_DENSAO (DMO,DMAT)
         CASE DEFAULT
            WRITE (LUPRI,*) " Unknown or unimplemented basis: ",BASIS
            WRITE (LUPRI,*) " Use density matrix in AO basis "
            IERR = IERR + 2
            CALL DENSMO_TO_DENSAO (DMO,DMAT)
      END SELECT

      DEALLOCATE(DMO)
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_EOM_DENSITY (BASIS,STATE,DMAT,IERR)
C
C---------------Description--------------------------------------------
C
C     Retrieve relaxed density matrix for processing by other modules
C     of DIRAC (currenty only first-order properties)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Common Blocks--------------------------------------
C
#include "implicit.h"
#include "files.inc"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
C
C---------------Calling variables--------------------------------------
C
      PARAMETER (HALF=0.5D0)
      CHARACTER*2 BASIS                  ! Representation
      REAL*8 DMAT(*)                     ! Density matrix
      INTEGER IERR                       ! Error code
      INTEGER STATE                      ! Number of state for which density has been stored  
C
C---------------Local variables--------------------------------------
C
      REAL*8, ALLOCATABLE:: DMO(:,:,:)   ! Density matrix in MO basis
      LOGICAL FND,FNDLAB,ADD_2_FOCK
      CHARACTER*8 LABEL
      CHARACTER*6 FILE_NAME
      INTEGER dlength,dstart  
      INTEGER, PARAMETER :: EOMDEN=299
C
C---------------Executable code--------------------------------------

            ADD_2_FOCK = .FALSE.

            WRITE (LUPRI,*) " Retrieving EOMCC density matrix from file"
           ADD_2_FOCK   = .TRUE.
!              ADD_2_FOCK = .FALSE.
            IERR         = 0

C     Read the appropriate correlated density from file

      ALLOCATE (DMO(NORBT,NORBT,NZ),STAT=IERR)
      IF (IERR.NE.0) THEN
        CALL QUIT('Error in DMO allocation !')
      ENDIF

!      call waio_open(itapt+7)
      
      inquire(iolength=dlength) DMO

      dstart = (state-1)*norbt*norbt*nz+1    
!      dlength = norbt*norbt*nz*8

      open (eomden,file='EOM_DENS',access='direct',status='old',
     &      recl=dlength)

      read (eomden,rec=state) dmo

!      write(*,*)(dmo(i,i,1),i=1,norbt) 
!      call rread (itapt+7,dmo,dlength,dstart)
      
C     Add the closed shell Hartree-Fock contribution
      IF(ADD_2_FOCK)THEN
C       The add_two routine scales by 2, scale back to the orginal afterwards
        CALL ADD_TWO_FOCK(DMO,.TRUE.,0)
        CALL DSCAL(NZ*NORBT*NORBT,HALF,DMO,1)
      END IF


      SELECT CASE (BASIS)
         CASE ("MO")
            CALL DCOPY (N2ORBTQ,DMO,1,DMAT,1)
         CASE ("AO")
            CALL DENSMO_TO_DENSAO (DMO,DMAT)
         CASE DEFAULT
            WRITE (LUPRI,*) " Unknown or unimplemented basis: ",BASIS
            WRITE (LUPRI,*) " Use density matrix in AO basis "
            IERR = IERR + 2
            CALL DENSMO_TO_DENSAO (DMO,DMAT)
      END SELECT

      deallocate(dmo)
      
      close (eomden,status='keep')
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE STORE_CC_DENSITY (WF,DMO)
C
C---------------Description--------------------------------------------
C
C     Retrieve relaxed density matrix for processing by other modules
C     of DIRAC (currenty only first-order properties)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C


C---------------Common Blocks--------------------------------------
C
#include "implicit.h"
#include "files.inc"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbprp.h"
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*4 WF                     ! Wave function type
      REAL*8 DMO(*)                      ! Density matrix (MO-basis)
C
C---------------Local variables--------------------------------------
C
      CHARACTER*8 :: LABEL,EOFLABEL
      LOGICAL, SAVE :: INIT=.FALSE.
      LOGICAL :: FNDLAB
C---------------Executable code--------------------------------------
 
      WRITE (LUPRI,*) "Storing the ",WF," density matrix on file"
      LABEL = WF//'DENS'
      EOFLABEL = 'EOFLABEL'

C     Write the correlated density to file
C     Code should be replaced by proper storage routine, this just adds to the end of a file 
C     without being able to do a replace or check for old data.

      CALL OPNFIL(MCCDENS,'CCDENS','UNKNOWN','STORE_')
      IF (.NOT.INIT) THEN
          CALL NEWLAB(EOFLABEL,MCCDENS,LUPRI)
          INIT = .TRUE.
      END IF
      REWIND (MCCDENS)
      IF (FNDLAB(EOFLABEL,MCCDENS))
     & BACKSPACE MCCDENS
      CALL NEWLAB(LABEL,MCCDENS,LUPRI)
! aspg, debug of properties for exacorr
! keep this code commented out for the moment
!     write (*,'(3X,A,2I6)') '>>>>> debug : norbt, nz',norbt,nz 
!     print *,'>>>>> DMO in dirac format'
!     do i=1, NORBT*NORBT*NZ
!        write(*,'(3X,I4,2x,F28.10)') i, dmo(i)
!     end do
!     print *,'<<<<< DMO in dirac format'
      CALL WRITT(MCCDENS,NORBT*NORBT*NZ,DMO)
      CALL NEWLAB(EOFLABEL,MCCDENS,LUPRI)
      CLOSE(Unit=MCCDENS,Status='KEEP')

C     Switch on evaluation of property output
      DOEXP = .TRUE.
C     Add this density matrix to the list to be evaluated
      IPRP_WF = 0
      DO I = 1, NPRP_WF
         IF (PRP_WF(I).EQ.WF) IPRP_WF = I
      END DO
      IF (IPRP_WF.EQ.0) THEN
         IF (NPRP_WF .GE. MXPRP_WF) THEN
             IPRP_WF = 1 ! Unlikely that this happens, just skip the HF print out in this case
         ELSE
             NPRP_WF = NPRP_WF + 1
             IPRP_WF = NPRP_WF
         END IF
         PRP_WF(IPRP_WF) = WF
      END IF

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE STORE_EOM_DENSITY (WF,STATE,DMO)
C
C---------------Description--------------------------------------------
C     Storing relaxed EOM density matrix for processing by other modules
C     of DIRAC (currenty only first-order properties)
C     This should work with multi-state methods as well...
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Avijit Shee
C
C---------------Common Blocks--------------------------------------
C
#include "implicit.h"
#include "files.inc"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*4 WF                     ! Wave function type
      REAL*8, INTENT(IN):: DMO(*)                      ! Density matrix (MO-basis)
      INTEGER STATE                      ! I designate each state by a Number
C
C---------------Local variables--------------------------------------
C
      CHARACTER*8 :: LABEL,EOFLABEL
      LOGICAL, SAVE :: INIT=.FALSE.
      LOGICAL :: FNDLAB
      logical            :: exist
      INTEGER, PARAMETER :: EOMDEN=299
      integer :: dlength
C---------------Executable code--------------------------------------
 
      WRITE (LUPRI,*) "Storing ",WF," density matrix for state",STATE,
     &                " on file"

C     Write the correlated density to file
C     Code should be replaced by proper storage routine, this just adds to the end of a file 
C     without being able to do a replace or check for old data.

      inquire(iolength=dlength) DMO(1:norbt*norbt*nz)

      inquire(file='EOM_DENS', exist=exist)
      if (exist) then

      open (eomden,file='EOM_DENS',access='direct',status='old',
     &      recl=dlength)

      else

       open (eomden,file='EOM_DENS',access='direct',status='new',
     &      recl=dlength)
     
      endif

      write (eomden,rec=state) dmo(1:norbt*norbt*nz)

      close (eomden,status='keep')

      print *,"done for state ",state
      END
