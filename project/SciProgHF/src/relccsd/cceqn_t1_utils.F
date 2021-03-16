      SUBROUTINE T1EQN1 (TAU,S1)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Evaluate part of T1 or Lambda1 equations
C     The same contraction is used in the construction of the
C     X-intermediate.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "ccpar.inc"

      real*8,intent(inout) :: S1(ndimt1*rcw)
      real*8,intent(in)    :: TAU(ndimt2*rcw)

C---------------Local variables--------------------------------------
C
      real*8, allocatable ::  tau_sorted(:), vooo(:), vooo_sorted(:)
      
C---------------Executable code--------------------------------------
C
      allocate (tau_sorted(LVVOOT(NREP+1)*RCW))
      allocate (vooo(IVOOOT(NREP+1)*RCW))
      allocate (vooo_sorted(LOVOOT(NREP+1)*RCW))

C--------------------------------------------
C S(A,I) = S(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------

C ** TAU now reorganised in BUF1: NV*NVOOT

      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             TAU,tau_sorted)
      CALL GETVOOO (vooo)
      CALL DELINT ('VOOO','DKDD',vooo,0,0,0)

C     --------------------------------
C     Sort from 12,34 to 134,2 (VOO,O)
C     --------------------------------
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           vooo,vooo_sorted)

C ** VOOO integrals now reorganised in BUF2: NVOOT * NO

C     ------------------------
C     V(K,L,C,I) = V*(C,I,K,L)
C     ------------------------
      IF (CARITH) CALL CONJUGA (LOVOOT(NREP+1),vooo_sorted,1)

C
      CALL CNTRCT 
     &   ('N','N',NV,NO,NVOOT,A1,tau_sorted,vooo_sorted,A1,S1,NREP)
C

      deallocate (tau_sorted)
      deallocate (vooo)
      deallocate (vooo_sorted)

      RETURN
      END
      SUBROUTINE T1EQN2 (TAU,BUF3,NBUF3,S1)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Evaluate part of T1 or Lambda1 equations
C     The same contraction is used in the construction of the
C     X-intermediate.
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
      REAL*8 S1(*),TAU(*)
      REAL*8 BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
      REAL*8 SDOT,DDOT
      integer a,ai,ak,cdik,i,ik,imin,irep,irp,istart,k,krep,nint
      integer t1ai,t1ak
C
C---------------Executable code--------------------------------------
C
C--------------------------------------------
C S(A,I) = S(A,I) + V(A,K,C,D) * TAU(C,D,I,K)
C--------------------------------------------
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         NINT = NBUF3/NVO(IRP)  ! maximum # of NVO batches

C ** NINT is updated on exit of GETVOVV !

  10     CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,NINT)
         CALL DELINT ('VOVV','KDDD',BUF3,IRP,ISTART,NINT)
         IK = 0
         DO 20 KREP = 1, NREP
           IREP = MULTB(KREP,IRP+NREP,2)
           IF (KREP.GT.IREP) GOTO 20
           DO K = 1, NO(KREP)
              IMIN = 1
              IF (KREP.EQ.IREP) IMIN = K + 1
              DO I = IMIN, NO(IREP)
                 IK = IK + 1
                 CDIK = (IVVOOTT(IRP)+(IK-1)*NVVT(IRP)+ISTART)*RCW+1
                 DO A = 1, NV(IREP)
                    AK = (IIVO(IREP,KREP)+(K-1)*NV(IREP)+A-1) * RCW + 1
                    T1AI = (IVO(IREP)+(I-1)*NV(IREP)+A-1) * RCW + 1
                    IF (CARITH) THEN
                       S1(T1AI) = S1(T1AI)
     &            + DDOT(NINT,TAU(CDIK),RCW,BUF3(AK),RCW*NVO(IRP))
     &            - DDOT(NINT,TAU(CDIK+1),RCW,BUF3(AK+1),RCW*NVO(IRP))
                       S1(T1AI+1) = S1(T1AI+1)
     &            + DDOT(NINT,TAU(CDIK+1),RCW,BUF3(AK),RCW*NVO(IRP))
     &            + DDOT(NINT,TAU(CDIK),RCW,BUF3(AK+1),RCW*NVO(IRP))
                    ELSE
                       S1(T1AI) = S1(T1AI)
     &            + DDOT(NINT,TAU(CDIK),1,BUF3(AK),NVO(IRP))
                    ENDIF
                 ENDDO
                 DO A = 1, NV(KREP)
                    AI = (IIVO(KREP,IREP)+(I-1)*NV(KREP)+A-1) * RCW + 1
                    T1AK = (IVO(KREP)+(K-1)*NV(KREP)+A-1) * RCW + 1
                    IF (CARITH) THEN
                       S1(T1AK) = S1(T1AK)
     &            - DDOT(NINT,TAU(CDIK),RCW,BUF3(AI),RCW*NVO(IRP))
     &            + DDOT(NINT,TAU(CDIK+1),RCW,BUF3(AI+1),RCW*NVO(IRP))
                       S1(T1AK+1) = S1(T1AK+1)
     &            - DDOT(NINT,TAU(CDIK+1),RCW,BUF3(AI),RCW*NVO(IRP))
     &            - DDOT(NINT,TAU(CDIK),RCW,BUF3(AI+1),RCW*NVO(IRP))
                    ELSE
                       S1(T1AK) = S1(T1AK)
     &            - DDOT(NINT,TAU(CDIK),1,BUF3(AI),NVO(IRP))
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
 20      CONTINUE
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
 30   CONTINUE
C
      RETURN
      END
