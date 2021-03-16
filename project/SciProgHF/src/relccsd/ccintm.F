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
      SUBROUTINE AINTM(T1,T2,BUF1,BUF2,AZ)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates A intermediate
C
C---------------Routines called----------------------------------------
C
C     GETTAU
C     GETOOOO, GETVVOO, GETVOOO
C     XGEMM, XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
      REAL*8 AZ(*),BUF1(*),BUF2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      real*8 ddot

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface


C---------------Executable code--------------------------------------
C
C--------------------
C A(KL,IJ) = W(KL,IJ)
C--------------------
      CALL GETOOOO (AZ)
      CALL DELINT ('OOOO','DDKK',AZ,0,0,0)
C------------------------------------------------
C A(KL,IJ) = A(KL,IJ) + PIJ ( W(KL,IC) * T(C,J) )
C------------------------------------------------
      CALL GETVOOO (BUF1)
      CALL DELINT ('VOOO','DKDD',BUF1,0,0,0)
      CALL SRT9 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NOOOT,KOOOVT,KKOOOT,
     &           BUF1,BUF2)
      CALL XSCAL (KOOOVT(NREP+1),-A1,BUF2,1)
      IF (CARITH) CALL CONJUGA (KOOOVT(NREP+1),BUF2,1)
      CALL CNTRCT ('N','N',NOOOT,NO,NV,A1,BUF2,T1,A0,BUF1,NREP)
      CALL SRT1T3 (NREP,MULTB,LTR,NOOT,NO,NO,NOOOT,KOOOOT,KKOOOT,
     &             BUF1,BUF2)
      CALL XAXPY (IOOOOTT(NREP+1),A1,BUF2,1,AZ,1)
C---------------------------------------------
C A(KL,IJ) = A(KL,IJ) +  W(KL,CD) * TAU(CD,IJ)
C---------------------------------------------
      CALL GETVVOO (BUF1)
      CALL DELINT ('VVOO','DDDD',BUF1,0,0,0)
      CALL GETTAU (T1,T2,BUF2)
      CALL CNTRCT ('C','N',NOOT,NOOT,NVVT,A1,BUF1,BUF2,A1,AZ,NREP)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GOINTM (FVO,T1,BUF1,BUF2,GO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates GO intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVO(*)
      REAL*8 T1(*)
      REAL*8 GO(*)
      REAL*8 BUF1(*),BUF2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      integer m,n
C
C---------------Executable code--------------------------------------
C
C----------------------------------
C G(K,I) = H(K,I) + F(K,C) * T(C,I)
C----------------------------------
      CALL XCOPY (NFVO,FVO,1,BUF1,1)
      CALL DELFCK ('VO','DD',BUF1)
      CALL CNTRCT ('C','N',NO,NO,NV,A1,BUF1,T1,A1,GO,NREP)
C-----------------------------------------------------
C G(KI) = G(KI) + W(KL,IC) * T(CL) = G(KI) - W*(CI,KL)
C Order W from CI,KL to CL,KI
C-----------------------------------------------------
      M = MOO(1)
      N = MVO(1)
      IF (N.EQ.0) RETURN
      CALL GETVOOO (BUF1)
      CALL DELINT ('VOOO','DKDD',BUF1,0,0,0)
      CALL SRT26 (NREP,MULTB,LTR,LFA,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &            BUF1,BUF2)
      CALL XGEMV ('C',N,M,-A1,BUF2,N,T1,1,A1,GO,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GVINTM(FVO,T1,BUF3,NBUF3,BUF2,NBUF2,GV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates GV intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVO(*)
      REAL*8 T1(*)
      REAL*8 GV(*)
      REAL*8 BUF3(*),BUF2(*)
      INTEGER NBUF3,NBUF2
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL USEDZ,RIGHT
      real*8 ddot
      integer m,n
C
C---------------Executable code--------------------------------------
C
C ********************
C ** TASK SPLITTING **
C ********************
C ** only master calculates local contr. to G(A,C) !
C ** all nodes calculate nonlocal contr.

      IF(MYPROC.EQ.MASTER) THEN

C----------------------------------
C G(A,C) = H(A,C) - T(A,K) * F(K,C)
C----------------------------------
        CALL XCOPY (NFVO,FVO,1,BUF2,1)
        CALL DELFCK ('VO','DD',BUF2)
        CALL CNTRCT ('N','C',NV,NV,NO,-A1,T1,BUF2,A1,GV,NREP)
  
      ELSE

        CALL XCOPY(NFVV,A0,0,GV,1)

      ENDIF
C ***************************
C ** END OF TASK SPLITTING **
C ***************************

C------------------------------------
C G(A,C) = G(A,C) + W(AK,CD) * T(D,K)
C------------------------------------
      M = MVV(1)
      N = MVO(1)
      IF (N.EQ.0) RETURN
      USEDZ = .FALSE.
      RIGHT = .TRUE.
      CALL SRT20D(NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
     &            JJVO,BUF2,NBUF2,T1,GV,USEDZ,RIGHT)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HINTM(T1,T2,BUF1,BUF2,BUF3,NBUF3,H)
C
      use interface_to_mpi
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates H intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF3
      REAL*8 T1(*),T2(*),H(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 T1VAL(2)
      DATA T1VAL /2*0.0D0/
      LOGICAL DONE,TEQ
      real*8 ddot,xsum
      integer akcd,akci,akdi,c,cco,cd,ci,cmin,crep,d,di,drep,i,irep,irp
      integer istart,jrep,k,kld,m,mint,n,nint,off1,off2,off3,t1ci,t1di
C
C---------------Executable code--------------------------------------

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface

C
      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.
      ELSEIF (EQNS.EQ.'FOCKSP') THEN
         TEQ = .TRUE.
      ELSE
         CALL QUIT('Wrong keyword in call to HINTM')
      ENDIF

C ********************
C ** TASK SPLITTING **
C ********************
C ** master calculates local contr. 
C ** - W(AK,CI)  + T(A,L) * W(LK,CI) 
C ** to the H intermediate which is kept
C ** in BUF1 ! ===> If we want to start from a clear H array we have to
C ** clear BUF1 at the extent of the H length.

      IF(MYPROC.EQ.MASTER) THEN
C---------------------------------------------------------------------
C H(AK,IC) = - W(AK,CI)
C We store h initially as H(AK,CI), ordening and minus sign after next
C step !
C---------------------------------------------------------------------
      CALL GETVOVO (H)
      CALL DELINT ('VOVO','KDDK',H,0,0,0)

C----------------------------------------
C H(AK,CI) = H(AK,CI) + T(A,L) * W(LK,CI)
C----------------------------------------
      CALL GETVOOO (BUF1)
      CALL DELINT ('VOOO','DKDD',BUF1,0,0,0)
      CALL SRT1C1 (NREP,NVO,NOOT,BUF1,BUF2)
      CALL SRT1T2 (NREP,MULTB,LFA,NOOT,NO,NO,NVO,LOOVO,LLOVO,BUF2,BUF1)
      CALL CNTRCT ('N','N',NV,NOVO,NO,A1,T1,BUF1,A0,BUF3,NREP)
      CALL SRT1S2 (NREP,MULTB,LTR,NVO,NV,NO,NVO,LVOVO,LLOVO,BUF3,BUF1)
      CALL XAXPY (IVOVO(NREP+1),-A1,H,1,BUF1,1)

      ELSE
    
        CALL XCOPY(NV4,A0,0,BUF1,1)

      ENDIF

C ***************************
C ** END OF TASK SPLITTING **
C ***************************
C ** here now come the distributed VOVV contributions to be calculated
C ** by each node !

C------------------------------------------------------------------
C H(AK,CI) = H(AK,CI) - W(AK,CD) * T(D,I)
C H is kept in BUF1, still ordered VOVO
C------------------------------------------------------------------
      CALL XTIME(0,1,'--- HINTM: VOVV*T              ')
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         CCO=0
         MINT = NBUF3/NVO(IRP)

C ** reading from GETVOVV: ISTART is updated with the actual batch #
C ** available on the local node, NINT is updated with the number of
C ** actually read integrals from this batch and is a pure output parameter !
C ** We therefore have to count upwards in case of ISTART.NE.0 !!

         CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,MINT)
         CALL DELINT ('VOVV','KDDD',BUF3,IRP,ISTART,NINT)

C ** now start SORTING and contraction loop

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
C                   ---------------------------------------
C                   We need the next buffer in this IRREP !
C                   ---------------------------------------
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
                    CI = IIVO(CREP,DREP)+(I-1)*NV(CREP)+C
                    AKCI = (IVOVO(IRP)+(CI-1)*NVO(IRP))*RCW+1
                    T1DI = (IVO(DREP)+(I-1)*NV(DREP)+D-1) * RCW + 1
                    T1VAL(1) = - T1(T1DI)
                    IF (CARITH) T1VAL(2) = - T1(T1DI+1)
                    CALL XAXPY(NVO(IRP),T1VAL,BUF3(AKCD),1,BUF1(AKCI),1)
                 ENDDO
                 DO I = 1, NO(CREP)
                    DI = IIVO(DREP,CREP)+(I-1)*NV(DREP)+D
                    AKDI = (IVOVO(IRP)+(DI-1)*NVO(IRP))*RCW+1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = T1(T1CI)
                    IF (CARITH) T1VAL(2) = T1(T1CI+1)
                    CALL XAXPY(NVO(IRP),T1VAL,BUF3(AKCD),1,BUF1(AKDI),1)
                 ENDDO
 10           CONTINUE
 15        CONTINUE
 20      CONTINUE
 30    CONTINUE

C------------------------------------------------------------------
C Now sort H(AK,CI) to H(AI,CK)  is to be done in any case !!
C The H buffer is not changed any more up to the final XGEMM !
C Every other calculation is done in BUF1 !!
C------------------------------------------------------------------
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            BUF1,H)
      CALL XTIME(0,2,'--- HINTM: VOVV*T              ')

C------------------------------------------------------------------
C H(AI,CK) = H(AI,CK) + [ 0.5*T(AD,IL) - T(A,L)*T(D,I) ] * W(KL,CD)
C------------------------------------------------------------------
      CALL XCOPY (NDIMT2,T2,1,BUF1,1)
      IF (TEQ) CALL XSCAL (NDIMT2,AP5,BUF1,1)
C----------------------------------------------------------------------
C For the H-intermediate in the lambda equations we should not scale t2
C----------------------------------------------------------------------
      CALL GETTAU (T1,BUF1,BUF2)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,
     &             JJVO,BUF2,BUF1)
C---------------------------------------------------------------
C BUF1 now contains 0.5*T(AD,IL) - T(A,L)*T(D,I) + T(A,I)*T(D,L)
C Subtract T(A,I)*T(D,L)
C---------------------------------------------------------------
      CALL XGEMM ('N','N',NDIMT1,NDIMT1,1,-A1,T1,NDIMT1,T1,1,
     &               A1,BUF1,NDIMT1)
C
      CALL GETVVOO (BUF3)
      CALL DELINT ('VVOO','DDDD',BUF3,0,0,0)
      IF (CARITH) CALL CONJUGA (NV3,BUF3,1)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              BUF3,BUF2)

      CALL XTIME(0,1,'--- HINTM: VVOO contribution      ')

C
C --- new parallelized VOVO contraction
C
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M    = MVO(IREP)
         N    = MVO(IREP)
         KLD  = MVO(JREP)
         K    = IDIST(1,4,JREP)
         OFF1 = IDIST(2,4,JREP)
         OFF2 = IDIST(3,4,JREP)
         CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,BUF2(OFF2),
     &               KLD,A1,H(OFF3),M)
         OFF3 = OFF3 + M * N * RCW
      ENDDO
 
      CALL XTIME(0,2,'--- HINTM: VVOO contribution      ')

#if defined (VAR_MPI)
C *******************
C ** parallel case: the H-intermediate was partially calculated
C ** on the separate slaves. In order to perform the H*T2
C ** contraction we need the complete H array ==> we complete it
C ** here and synchronize all the nodes. They then have the complete
C ** H-array and this lengthy contraction can also be parallelized.
C
      IF (NMPROC .GT. 1) THEN

        IF(NBUF3.LT.IVOVO(NREP+1))
     &    CALL QUIT('BUF3 can not keep the VOVO integrals !')

C ... we collect the H array
        CALL XTIME(4,1,'-- combining HINTM via MPI_ALLREDUCE')
        CALL XCOPY (NBUF3,A0,0,BUF3,1)
        call interface_mpi_allreduce_r1_work_f77(H,BUF3,
     &       RCW*IVOVO(NREP+1),
     &       op_MPI_SUM,global_communicator)
        CALL XCOPY(IVOVO(NREP+1),BUF3,1,H,1)
        CALL XTIME(4,2,'-- combining HINTM via MPI_ALLREDUCE')
      ENDIF
#endif

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HOINTM(FOO,EPS,T1,T2,BUF1,BUF2,BUF3,HO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates HO intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FOO(*),EPS(*)
      REAL*8 T1(*),T2(*)
      REAL*8 HO(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      real*8 ddot
      integer i,ii,irp,n,nd
C---------------Executable code--------------------------------------
      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface
C
C---------------------------------------
C H(K,I) = F(K,I)  - DELTA(K,I) * EPS(I)
C---------------------------------------
      CALL XCOPY (NFOO,FOO,1,HO,1)
      II = 1
      I = 1
      DO IRP = 1, NREP
         N = NO(IRP)
         ND = (N+1) * RCW
         CALL DAXPY(N,-A1,EPS(I),1,HO(II),ND)
         I = I + N
         II = II + N * N * RCW
      ENDDO
      CALL DELFCK ('OO','DK',HO)
C----------------------------------------------------------------------
C H(K,I) = H(K,I) + W(KL,CD) * T(CD,IL) = H(K,I) + W*(CDL,K) * T(CDL,I)
C----------------------------------------------------------------------
      CALL GETVVOO (BUF3)
      CALL DELINT ('VVOO','DDDD',BUF3,0,0,0)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             BUF3,BUF1)
      CALL GETTAU  (T1,T2,BUF3)

      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             BUF3,BUF2)
      CALL CNTRCT ('C','N',NO,NO,NVVOT,A1,BUF1,BUF2,A1,HO,NREP)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HOVINTM(FVO,T1,BUF1,BUF2,HOV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates HOV intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVO(*)
      REAL*8 T1(*)
      REAL*8 HOV(*)
      REAL*8 BUF1(*),BUF2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      integer k,m
C
C---------------Executable code--------------------------------------
C
C--------------------------------
C HOV(C,K) = FOV(K,C) = FVO*(C,K)
C--------------------------------
      CALL XCOPY (NFVO,FVO,1,HOV,1)
      CALL DELFCK ('VO','DD',HOV)
      IF (CARITH) CALL CONJUGA (NFVO,HOV,1)
C-------------------------------------------------------------------
C HOV(CK) = HOV(CK) + W(KL,CD) * T(DL) = HOV(CK) + W*(CD,KL) * T(DL)
C-------------------------------------------------------------------
      M = MVO(1)
      K = MVO(1)
      CALL GETVVOO (BUF1)
      CALL DELINT ('VVOO','DDDD',BUF1,0,0,0)
      CALL XCOPY (JVOVO(2),A0,0,BUF2,1)
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,MVO,JVOVO,JJVO,JJVO,
     &              BUF1,BUF2)
      IF (CARITH) CALL CONJUGA (JVOVO(2),BUF2,1)
      CALL XGEMV ('N',M,K,A1,BUF2,M,T1,1,A1,HOV,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HVINTM(FVV,EPS,T1,T2,BUF1,BUF2,BUF3,HV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates HV intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVV(*),EPS(*)
      REAL*8 T1(*),T2(*)
      REAL*8 HV(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
      REAL*8 ddot
      integer a,aa,irp,n,nd
C
C---------------Executable code--------------------------------------


      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface

C---------------------------------------
C H(A,C) = F(A,C)  - DELTA(A,C) * EPS(A)
C---------------------------------------
      CALL XCOPY (NFVV,FVV,1,HV,1)
      AA = 1
C     Calculate start of virtual orbital energies in eps
      A = 1
      DO IRP = 1, NREP
         A = A + NO(IRP)
      ENDDO
C
      DO IRP = 1, NREP
         N = NV(IRP)
         ND = (N+1) * RCW
         CALL DAXPY(N,-A1,EPS(A),1,HV(AA),ND)
         A = A + N
         AA = AA + N * N * RCW
      ENDDO
      CALL DELFCK ('VV','KD',HV)
C----------------------------------------------------------------------
C H(A,C) = H(A,C) - W(KL,CD) * T(AD,KL) = H(A,C) - T(A,DKL) * W*(C,DKL)
C----------------------------------------------------------------------
      CALL GETTAU (T1,T2,BUF3)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF3,BUF1)
      CALL GETVVOO (BUF3)
      CALL DELINT ('VVOO','DDDD',BUF3,0,0,0)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF3,BUF2)
      CALL CNTRCT ('N','C',NV,NV,NVOOT,-A1,BUF1,BUF2,A1,HV,NREP)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BINTM(T1,T2,BUF1,BUF2,BUF3,NBUF3,S2)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates B intermediate and adds it to new T2's
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF3
      REAL*8 T1(*),T2(*)
      REAL*8 S2(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "complex.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE,TEQ
      REAL*8 T1MIN(2)
      DATA T1MIN /2*0.0D0/
      REAL*8 ddot
      integer irp,istart,k,kd,m,mint,n,nint,off1,off2
C
C---------------Executable code--------------------------------------

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface

      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.
      ELSEIF (EQNS.EQ.'FOCKSP') THEN
         TEQ = .TRUE.
      ELSE
         CALL QUIT('Wrong keyword in call to BINTM')
      ENDIF
C---------------------------------------------
C S(AB,IJ) =  S(AB,IJ) + V(AB,CD) * TAU(CD,IJ)
C---------------------------------------------
      IF (TEQ) THEN
         CALL GETTAU (T1,T2,BUF1)
      ELSE
         CALL XCOPY  (NDIMT2,T2,1,BUF1,1)
      ENDIF
      OFF1 = 1
      OFF2 = 1
      DO 30 IRP = 1, NREP
         ISTART = 0
         IF (NVVT(IRP).EQ.0) GOTO 30
         MINT = NBUF3/NVVT(IRP)
  10     CALL GETVVVV (IRP,ISTART,NINT,DONE,BUF3,MINT)
         CALL DELINT ('VVVV','KKDD',BUF3,IRP,ISTART,NINT)
         M = NVVT(IRP)
         N = NOOT(IRP)
         K = NINT
         CALL XGEMM ('N','N',M,N,K,A1,BUF3,M,BUF1(OFF1+ISTART*RCW),M,
     &               A1,S2(OFF2),M)
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
         OFF1 = OFF1 + M * N * RCW
         OFF2 = OFF2 + M * N * RCW
 30   CONTINUE
C
C--------------------------------------------------------------
C S(AB,IJ) =  S(AB,IJ) - PAB [ T(B,K) * W(AK,CD) * TAU(CD,IJ) ]
C--------------------------------------------------------------
      CALL XCOPY (NV2,A0,0,BUF2,1)
C
      OFF1 = 1
      OFF2 = 1
      DO 130 IRP = 1, NREP
         ISTART = 0
         IF (NVVT(IRP).EQ.0.OR.NVO(IRP).EQ.0) GOTO 130
         M = NVO (IRP)
         N = NOOT(IRP)
         KD = NVVT(IRP)
         MINT = NBUF3/NVO(IRP)

 110     CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,MINT)
         CALL DELINT ('VOVV','KDDD',BUF3,IRP,ISTART,NINT)

         K = NINT

C ** now we construct the full NVO x NOOT matrix
C ** but the contributions are only from the restricted CD range
C ** the method is additive since ZGEMM works according to
C **      C=A*B + C  ==> all the already computed contributions
C ** will be added in this IRREP.

         CALL XGEMM ('N','N',M,N,K,A1,BUF3,M,BUF1(OFF1+ISTART*RCW),
     &               KD,A1,BUF2(OFF2),M)
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 110
         ENDIF
         OFF1 = OFF1 + KD * N * RCW
         OFF2 = OFF2 + M * N * RCW
 130  CONTINUE

C ** now the <AK||IJ> contr. are available in memory sorted according
C ** to combined IRREPS (AK) and (IJ).

C ** next sort resorts according to (KA) and (IJ)

      CALL SRT1L1 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,IVOOOT,IIOV,BUF2,BUF3)

C ** next sort resorts according to K,AIJ which means sorted according to
C ** IRREPS K and not combined IRREPS anymore for subsequent contraction

      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NO,NV,NOOT,LOVOOT,LLVOOT,
     &             BUF3,BUF2)
C---------------------------------------------------------------
C Contract with T1 amplitudes instead of L1 for Lambda equations
C---------------------------------------------------------------
      IF (TEQ) THEN
         CALL CNTRCT ('N','N',NV,NVOOT,NO,-A1,T1,BUF2,A0,BUF3,NREP)
      ELSE
         CALL GETAMPT (BUF1,BUF2)
         CALL CNTRCT ('N','N',NV,NVOOT,NO,-A1,BUF1,BUF2,A0,BUF3,NREP)
      ENDIF
      CALL SRT1T2 (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF3,BUF2)
      CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETTAU (T1,T2,TAU)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates Tau from T1 and T2 amplitudes
C
C---------------Routines called----------------------------------------
C
C     XCOPY, XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8,intent(in) :: T1(*),T2(*)
      REAL*8,intent(out) :: TAU(*)
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8  T1MIN(2)
      DATA    T1MIN /2*0.0D0/
      integer ab,abijoff,ai,aj,amin,arp,b,bi,bj,brp,i,ijrp,imin,irp,j
      integer jrp,na
    
C
C---------------Executable code--------------------------------------
C

      CALL XCOPY (NDIMT2,T2,1,TAU,1)
C
      ABIJOFF = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2) ! Note that IJRP is a boson irrep and thus shifteb by NREP
      IF (IRP.LT.JRP) GOTO 10      ! defines an upper triangular loop on occupied fermion irreps
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)      ! defines an upper triangular loop on occupied spinors: I > J
            DO 20 BRP = 1, NREP
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20 ! defines a triangular loop on virtual fermion irreps
C--------------------------------------------------------
C TAU(AB,IJ) = T2(AB,IJ) + T1(A,I) * T1(B,J)
C--------------------------------------------------------

            IF (BRP.EQ.JRP) THEN
               AB = 1
               BJ = (IVO(JRP) + (J-1) * NV(JRP)) * RCW + 1 ! off-set to (B=1,J)
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1 ! defines a triangular loop A > B
                  AI = (IVO(IRP) + (I-1) * NV(IRP) + AMIN - 1) * RCW +1 ! off-set to (A=1,I)
                  NA = NV(ARP) - AMIN + 1
                  CALL XAXPY (NA,T1(BJ),T1(AI),1,TAU(ABIJOFF+AB),1) ! T1(BJ) is a scalar here
                  AB = AB + NA * RCW
                  BJ = BJ + RCW
               ENDDO
            ENDIF

C--------------------------------------------------------
C TAU(AB,IJ) = T2(AB,IJ) - T1(A,J) * T1(B,I)
C--------------------------------------------------------
            IF (ARP.EQ.JRP) THEN
               AB = 1
               BI = (IVO(IRP) + (I-1) * NV(IRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AJ = (IVO(JRP) + (J-1) * NV(JRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  T1MIN(1) = -T1(BI)
                  IF (CARITH) T1MIN(2) = -T1(BI+1)
                  CALL XAXPY (NA,T1MIN,T1(AJ),1,TAU(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BI = BI + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------
C UPDATE OFFSET AND GO TO NEXT IRREP PAIR
C--------------------------------------------------------
            IF (ARP.NE.BRP) THEN
               ABIJOFF = ABIJOFF + NV(ARP) * NV(BRP) * RCW
            ELSE
               ABIJOFF = ABIJOFF + NV(ARP) * (NV(ARP)-1) * RCW / 2
            ENDIF
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE GETTAU2 (S1,DL,TAU)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates Hamiltonian matrix from S1 and S2 elements
C     Essentially Slater rules for 2-e determinants !
C
C     General formula (indices are active virtual or active occupied)
C     We apply it to the whole array and extract the relevant part
C     later. Note that TAU should already be initialized with S2.
C
C     H2(AB,IJ) = S2(AB,IJ) +
C               { S1(A,I)*DELTA(B,J) + S1(B,J)*DELTA(A,I)
C                - S1(A,J)*DELTA(B,I) - S1(B,I)*DELTA(A,J) }
C
C---------------Routines called----------------------------------------
C
C     XCOPY, XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 S1(*),DL(*),TAU(*)
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 S1MIN(2)
      DATA S1MIN /2*0.0D0/
      integer na,n,jrp,j,irp,imin,ijrp,i,f,e,brp,bj,bi,b,arp,amin,aj
      integer ai,abijoff,ab,a,m
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (NDIMT1,A0,0,DL,1)
      AI = 1
      DO IRP = 1, NREP
      DO I = 1, NO(IRP)
         F = I - NSO(IRP)
         N = I - NIO(IRP)
         DO A = 1, NV(IRP)
             E = A
             M = A - NSV(IRP)
             IF (E.EQ.F) DL(AI) = D1 
             IF (M.EQ.N.AND.M.GT.0) DL(AI) = D1 
             AI = AI + RCW
         ENDDO
      ENDDO
      ENDDO
C
      ABIJOFF = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            DO 20 BRP = 1, NREP
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
C-----------------------------------------------------------------------
C TAU(AB,IJ) = S2(AB,IJ) + S1(A,I) * DELTA (B,J) + DELTA(A,I) * S1 (B,J)
C-----------------------------------------------------------------------
            IF (BRP.EQ.JRP) THEN
               AB = 1
               BJ = (IVO(JRP) + (J-1) * NV(JRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AI = (IVO(IRP) + (I-1) * NV(IRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  CALL XAXPY (NA,DL(BJ),S1(AI),1,TAU(ABIJOFF+AB),1)
                  CALL XAXPY (NA,S1(BJ),DL(AI),1,TAU(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BJ = BJ + RCW
               ENDDO
            ENDIF
C-----------------------------------------------------------------------
C TAU(AB,IJ) = S2(AB,IJ) - S1(A,J) * DELTA (B,I) - DELTA(B,I) * S1 (A,J)
C-----------------------------------------------------------------------
            IF (ARP.EQ.JRP) THEN
               AB = 1
               BI = (IVO(IRP) + (I-1) * NV(IRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AJ = (IVO(JRP) + (J-1) * NV(JRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  S1MIN(1) = -S1(BI)
                  IF (CARITH) S1MIN(2) = -S1(BI+1)
                  CALL XAXPY (NA,S1MIN,DL(AJ),1,TAU(ABIJOFF+AB),1)
                  S1MIN(1) = -DL(BI)
                  IF (CARITH) S1MIN(2) = -DL(BI+1)
                  CALL XAXPY (NA,S1MIN,S1(AJ),1,TAU(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BI = BI + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------
C UPDATE OFFSET AND GO TO NEXT IRREP PAIR
C--------------------------------------------------------
            IF (ARP.NE.BRP) THEN
               ABIJOFF = ABIJOFF + NV(ARP) * NV(BRP) * RCW
            ELSE
               ABIJOFF = ABIJOFF + NV(ARP) * (NV(ARP)-1) * RCW / 2
            ENDIF
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETTAU3 (T1,T2,TAU)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates Tau3 from T1 and T2 amplitudes
C     Tau3 is used in the -T correction : T3 = T2 + 1/3 PIJ T1.T1
C
C---------------Routines called----------------------------------------
C
C     XCOPY, XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*),TAU(*)
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 T1VAL(2)
      DATA T1VAL /2*0.0D0/
      integer ab,abijoff,ai,aj,amin,arp,b,bi,bj,brp,i,ijrp,imin,irp
      integer j,jrp,na
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (NDIMT2,T2,1,TAU,1)
C
      ABIJOFF = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            DO 20 BRP = 1, NREP
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
C--------------------------------------------------------
C TAU(AB,IJ) = T2(AB,IJ) + T1(A,I) * T1(B,J)
C--------------------------------------------------------
            IF (BRP.EQ.JRP) THEN
               AB = 1
               BJ = (IVO(JRP) + (J-1) * NV(JRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AI = (IVO(IRP) + (I-1) * NV(IRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  T1VAL(1) = T1(BJ) * ARP3
                  IF (CARITH) T1VAL(2) = T1(BJ+1) * ARP3
                  CALL XAXPY (NA,T1VAL,T1(AI),1,TAU(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BJ = BJ + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------
C TAU(AB,IJ) = T2(AB,IJ) - T1(A,J) * T1(B,I)
C--------------------------------------------------------
            IF (ARP.EQ.JRP) THEN
               AB = 1
               BI = (IVO(IRP) + (I-1) * NV(IRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AJ = (IVO(JRP) + (J-1) * NV(JRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  T1VAL(1) = -T1(BI) * ARP3
                  IF (CARITH) T1VAL(2) = -T1(BI+1) * ARP3
                  CALL XAXPY (NA,T1VAL,T1(AJ),1,TAU(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BI = BI + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------
C UPDATE OFFSET AND GO TO NEXT IRREP PAIR
C--------------------------------------------------------
            IF (ARP.NE.BRP) THEN
               ABIJOFF = ABIJOFF + NV(ARP) * NV(BRP) * RCW
            ELSE
               ABIJOFF = ABIJOFF + NV(ARP) * (NV(ARP)-1) * RCW / 2
            ENDIF
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CNTRCT (TRANSA,TRANSB,MDIM,NDIM,KDIM,ALPHA,A,B,
     &                   BETA,C,NREP)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generalization of XGEMM to contract symmetry packed matrices
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 A(*),B(*),C(*)
      COMPLEX*16 ALPHA,BETA
      INTEGER MDIM(NREP),NDIM(NREP),KDIM(NREP)
      CHARACTER*1 TRANSA,TRANSB
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      integer nrep,irp,k,lda,ldb,ldc,m,n,off1,off2,off3
C
C---------------Executable code--------------------------------------
C
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = MDIM(IRP)
         N = NDIM(IRP)
         K = KDIM(IRP)
         IF (TRANSA.EQ.'N'.OR.TRANSA.EQ.'n') THEN
            LDA = M
         ELSE
            LDA = K
         ENDIF
         IF (TRANSB.EQ.'N'.OR.TRANSB.EQ.'n') THEN
            LDB = K
         ELSE
            LDB = N
         ENDIF
         LDC = M
         IF (K.EQ.0) THEN
            CALL XSCAL (M*N,BETA,C(OFF3),1)
         ELSE
            CALL XGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A(OFF1),LDA,
     &                  B(OFF2),LDB,BETA,C(OFF3),LDC)
         ENDIF
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
      RETURN
      END

!###############################################################################

      SUBROUTINE CNTRCT_SP (TRANSA,TRANSB,MDIM,NDIM,KDIM,ALPHA,A,B,
     &                   BETA,C)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generalization of XGEMM to contract symmetry packed matrices
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
#include "symm.inc"
      REAL*8 A(*),B(*),C(*)
      COMPLEX*16 ALPHA,BETA
      INTEGER MDIM(NREP),NDIM(NREP),KDIM(NREP)
      CHARACTER*1 TRANSA,TRANSB
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
      INTEGER :: OFFXX(NREP+1),OFFYY(NREP+1),OFFZZ(NREP+1),IRP,JRP
      integer k,lda,ldb,ldc,m,n,off1,off2,off3
C---------------Executable code--------------------------------------
C

       OFFXX(1) = 0
       OFFYY(1) = 0
       OFFZZ(1) = 0

      DO IRP = 1, NREP
         JRP = MULTB(IRP+NREP,1+NREP,2)
         OFFXX(IRP+1) = OFFXX(IRP) + MDIM(IRP)*KDIM(JRP) 
         if (TRANSB.EQ.'N'.OR.TRANSB.EQ.'n') then
          OFFYY(IRP+1) = OFFYY(IRP) + KDIM(IRP)*NDIM(JRP)
         else
          OFFYY(IRP+1) = OFFYY(IRP) + NDIM(IRP)*KDIM(JRP)
         endif
         OFFZZ(IRP+1) = OFFZZ(IRP) + MDIM(IRP)*NDIM(IRP) 
      ENDDO

      DO IRP = 1, NREP

         JRP = MULTB(IRP+NREP,1+NREP,2)

         OFF1 = 1 + OFFXX(JRP) * RCW
         OFF2 = 1 + OFFYY(IRP) * RCW
         OFF3 = 1 + OFFZZ(IRP) * RCW

         M = MDIM(IRP)
         N = NDIM(IRP)
         K = KDIM(JRP)
         IF (TRANSA.EQ.'N'.OR.TRANSA.EQ.'n') THEN
            LDA = M
         ELSE
            LDA = K
         ENDIF
         IF (TRANSB.EQ.'N'.OR.TRANSB.EQ.'n') THEN
            LDB = K
         ELSE
            LDB = N
         ENDIF
         LDC = M
         IF (K.EQ.0) THEN
            CALL XSCAL (M*N,BETA,C(OFF3),1)
         ELSE
            CALL XGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A(OFF1),LDA,
     &                  B(OFF2),LDB,BETA,C(OFF3),LDC)
         ENDIF

          OFF1 = OFF1 + M*K

      ENDDO
C
      RETURN
      END
