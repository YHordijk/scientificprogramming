      SUBROUTINE amplitude_equation_T2 (T1,T2,
     &          GO,GV,HOV,H,AZ,BUF1,BUF2,BUF3,NBUF3,S2)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Evaluate T2 equations
C
C---------------Routines called----------------------------------------
C
C     BLAS routines :                 XGEMM, XCOPY
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
      REAL*8 GO(*),GV(*),HOV(*)
      REAL*8 AZ(*),H(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "eqns.inc"
#include "ihm.inc"
#include "symm.inc"
#include "ccpar.inc"
C
      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface

C---------------Local variables--------------------------------------
C
      LOGICAL DONE,TEQ
      REAL*8 T1VAL(2)
      DATA T1VAL /2*0.0D0/
      complex*16 zdotc
      real*8 ddot
      integer abji,c,cj,crep,i,imin,irep,irp,istart,j,ij,jmin,jrep,k,kld
      integer m,mint,n,nint,off1,off2,off3,t1ci,abij,ji
C
C---------------Executable code--------------------------------------
C
C---- Initialize to zero...
      CALL XCOPY (NDIMT2,A0,0,S2,1)
C---------------------------------------------
C S(AB,IJ) = S(AB,IJ) + TAU(AB,KL) * AZ(KL,IJ)
C---------------------------------------------
      CALL XTIME(0,1,'--- T2EQNS: TAU*AINTM contraction')
      CALL GETTAU (T1,T2,BUF1)
          
      OFF3 = 1
      DO IREP = 1, NREP
        M    = NVVT(IREP)
        N    = NOOT(IREP)
        KLD  = NOOT(IREP)
        K    = IDIST(1,5,IREP)
        OFF1 = IDIST(2,5,IREP)
        OFF2 = IDIST(3,5,IREP)
        CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,AZ(OFF2),KLD,
     &              A1,S2(OFF3),M)
        OFF3 = OFF3 + M * N * RCW
      ENDDO

C ***
C *** original CNTRCT call in serial version
C ***
C *** CALL CNTRCT ('N','N',NVVT,NOOT,NOOT,A1,BUF1,AZ,A1,S2,NREP)
C ***

      CALL XTIME(0,2,'--- T2EQNS: TAU*AINTM contraction')

C ** we will not parallelize a summation over just one occupied index
C ** if the system is small this can lead to zero batches. Only master
C ** does the T2(AB,IK) * GO(K,J) contraction

C ********************
C ** TASK SPLITTING **
C ********************

      IF(MYPROC.EQ.MASTER) THEN
C----------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) - PIJ [ T2(AB,IK) * GO(K,J) ]
C----------------------------------------------------------------
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             T2,BUF1)
      CALL CNTRCT ('N','N',NVVOT,NO,NO,A1,BUF1,GO,A0,BUF2,NREP)
      CALL SRT1T3 (NREP,MULTB,LTR,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             BUF2,BUF1)
      CALL XAXPY (NDIMT2,-A1,BUF1,1,S2,1)

C ***************************
C ** END OF TASK SPLITTING **
C ***************************

      ENDIF
C-----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PIJ [ W(AB,CJ) * T1(C,I) ]
C-----------------------------------------------------------------------
      CALL XTIME(0,1,'--- T2EQNS: VOVV*T1              ')
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         MINT = NBUF3/NVO(IRP)  ! max. # of NVO batches that fit in mem.
  10     CONTINUE
         CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,MINT)
         CALL DELINT ('VOVV','DKKK',BUF3,IRP,ISTART,NINT)
         IF (CARITH) CALL CONJUGA (NVO(IRP)*NINT,BUF3,1)
         IJ = 0
         DO 20 JREP = 1, NREP
           CREP = MULTB(JREP,IRP+NREP,2)
           IREP = CREP
           IF (IREP.LT.JREP) GOTO 20
           DO J = 1, NO(JREP)
              IMIN = 1
              IF (JREP.EQ.IREP) IMIN = J + 1
              DO I = IMIN, NO(IREP)
                 IJ = IJ + 1
                 ABIJ = (IVVOOTT(IRP)+(IJ-1)*NVVT(IRP)+ISTART)*RCW+1
                 DO C = 1, NV(CREP)
                    CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = T1(T1CI)
                    IF (CARITH) T1VAL(2) = T1(T1CI+1)
                    CALL XAXPY(NINT,T1VAL,BUF3(CJ),NVO(IRP),S2(ABIJ),1)
                 ENDDO
              ENDDO
           ENDDO
 20      CONTINUE
         JI = 0
         DO 21 IREP = 1, NREP
           JREP = MULTB(IREP,IRP+NREP,2)
           CREP = IREP
           IF (JREP.LT.IREP) GOTO 21
           DO I = 1, NO(IREP)
              JMIN = 1
              IF (JREP.EQ.IREP) JMIN = I + 1
              DO J = JMIN, NO(JREP)
                 JI = JI + 1
                 ABJI = (IVVOOTT(IRP)+(JI-1)*NVVT(IRP)+ISTART)*RCW+1
                 DO C = 1, NV(CREP)
                    CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = - T1(T1CI)
                    IF (CARITH) T1VAL(2) = - T1(T1CI+1)
                    CALL XAXPY(NINT,T1VAL,BUF3(CJ),NVO(IRP),S2(ABJI),1)
                 ENDDO
              ENDDO
           ENDDO
 21      CONTINUE
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
 30    CONTINUE
      CALL XTIME(0,2,'--- T2EQNS: VOVV*T1              ')

C ********************
C ** TASK SPLITTING **
C ********************
C ** here a VOOV dependent contribution (Term 7 in the paper) is to be
C ** calculated only by the master

      IF(MYPROC.EQ.MASTER) THEN

C--------------------------------
C Terms specific for T2 equations
C--------------------------------
                   
C-----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ [T1(A,K)*W(BK,IC)*T1(C,J)]
C-----------------------------------------------------------------------
         CALL GETVOVO (BUF1)
         CALL DELINT ('VOVO','KDDK',BUF1,0,0,0)
         CALL SRT19 (NREP,MULTB,LFA,NVO,NV,NO,NVOO,KVOOV,KKVOO,
     &               BUF1,BUF2)
         CALL CNTRCT ('N','N',NVOO,NO,NV,A1,BUF2,T1,A0,BUF1,NREP)
         CALL SRT1T3 (NREP,MULTB,LTR,NVO,NO,NO,NVOO,KVOOO,KKVOO,
     &                BUF1,BUF3)
C-------------------------------------------------------------------
C        Integrals are now in VO,OOT format : combine with next term
C-------------------------------------------------------------------
           

      ENDIF

C ****************************
C ** END OF TASK SPLITTING. **
C ****************************
C ** this contribution coming from the intermediate G(A,C) has to be
C ** calculated on each node. The fact that G(A,C) contains LOCAL +
C ** DISTRIBUTED parts as well is accounted for by the original calculation
C ** of the G(A,C) array before entering T2EQNS

C-------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB [ GV(A,C) * T2(CB,IJ)]
C-------------------------------------------------
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF2)
      CALL CNTRCT ('N','N',NV,NVOOT,NV,A1,GV,BUF2,A0,BUF1,NREP)

C ** the contraction of G(A,C) with the T2 amplitudes now is in BUF1 !
C ** these are the partial contributions available on EACH node


C ***
C ****** task splitting
C ***
      IF(MYPROC.EQ.MASTER) THEN

C------------------------------------------------
C S(AB,IJ) = S(AB,IJ) - PAB [ T1(A,K) *W(KB,IJ) ]
C------------------------------------------------
                   
         CALL GETVOOO (BUF2)
         CALL DELINT ('VOOO','KDKK',BUF2,0,0,0)
         CALL XAXPY (NV2,-A1,BUF3,1,BUF2,1)
          
           
      CALL SRT1L1 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,IVOOOT,IIOV,BUF2,BUF3)
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NO,NV,NOOT,LOVOOT,LLVOOT,
     &             BUF3,BUF2)
      CALL CNTRCT ('N','N',NV,NVOOT,NO,-A1,T1,BUF2,A1,BUF1,NREP)

C ***
C *** end of task splitting
C ***
      ENDIF

      CALL SRT1T2 (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF1,BUF2)
      CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)

C ** now the final H*T2 contraction. H is complete since we communicated
C ** it already.

C-----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ [H(AI,CK)*T2(CK,BJ)]
C-----------------------------------------------------------------------

      CALL XTIME(0,1,'--- T2EQNS: HINTM*T2           ')
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              T2,BUF1)
      CALL XCOPY (J2VOVO(NREP+1),A0,0,BUF2,1)

      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M    = MVO(IREP)
         N    = MVO(JREP)
         KLD  = MVO(IREP)
         K    = IDIST(1,10,JREP)
         OFF1 = IDIST(2,10,JREP)
         OFF2 = IDIST(3,10,JREP)
         CALL XGEMM ('N','N',M,N,K,A1,H(OFF1),M,BUF1(OFF2),KLD,
     &                A0,BUF2(OFF3),M)
         OFF3 = OFF3 + M * N * RCW
      ENDDO

      CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              BUF2,BUF1)
      CALL XAXPY (NDIMT2,A1,BUF1,1,S2,1)

      CALL XTIME(0,2,'--- T2EQNS: HINTM*T2           ')
C
C
C  Finally we ADD the B INTERMEDIATE CONTRIBUTION
C
C
      CALL XTIME(5,1,'-- BINTM                       ')
      CALL BINTM (T1,T2,BUF1,BUF2,BUF3,NBUF3,S2)
      CALL XTIME(5,2,'-- BINTM                       ')

      END
