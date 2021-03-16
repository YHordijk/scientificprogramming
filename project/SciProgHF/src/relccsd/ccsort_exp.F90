  module modified_sorting
  
  use symmetry_offset

   public SRT1SS4_exp
   public SRT16_exp
   public SRT1TT4_exp
   public SRT1TT5_exp
   public SRT1ST4_exp
   public SRT36_exp 
   public SRT46_exp 
   public SRT56_exp 
   public SRT22_exp
   public SRT32_exp
   public SRT26_exp
   public SRT1T3_exp
   public SRT1T2_exp
   public SRT1S2_exp
   public SRT19_exp
   public SRT1S3_exp
   public SRT7_exp
   public SRT6_exp
   public SRT9_exp
   public SRT10_exp
   public SRT1TS4_exp
   public SRT1C1_exp
   public nvo_to_mvo

     private

     contains

      SUBROUTINE SRT1TT4_exp (IRP_INPUT,NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     BUF(I>J,K>L:KLREP) <--> BUF(IK,JL)
!     If (FIRST) take only the totally symmetric (first) irrep.
!     If (DOINV) inverse sort.
!     IRP_INPUT specifies the IRREP to which the incoming tensors belong. 
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT,NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL FIRST,DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      type(Offset) :: d,e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(d,nrep,e,f,g,h)
      call auto_symmetry_offset_triangular(d,IFIE,JFIE)
      call auto_symmetry_offset_triangular(e,KFIE,LFIE)

      call auto_symmetry_offset(f,IFIE,KFIE,.false.,.false.)
      call auto_symmetry_offset(g,JFIE,LFIE,.false.,.false.)
      call auto_symmetry_offset(h,f%oneDirac,g%oneDirac,.true.,.true.)


      IF (.NOT.DOINV) THEN
         ROW = f%oneDirac 
         COLUMN = g%oneDirac
      IF (FIRST) THEN
        CALL XCOPY(f%oneDirac(1)*g%oneDirac(1),A0,0,BUF2,1)
      ELSE
        CALL XCOPY(h%oneNonDirac(IRP_INPUT),A0,0,BUF2,1)
      ENDIF
      ELSE
         ROW    = d%oneNonDirac
         COLUMN = e%oneNonDirac
      ENDIF

!     This loop shows the use of MULTB(*,*,2):
!     For instance, I have
!       ket(k)*ket(l) = KLREP = ket(l)*ket(k)
!     where the latter term comes from the Abelian nature of the point group
!     In the two outermost loops I fix KLREP and LREP and I want to know the
!     corresponding KREP.
!     I obtain this by the operation
!       bra(l)*KLREP = bra(l)*ket(l)*ket(k) = ket(k)
!     which is exactly the information obtained by
!       KREP = MULTB(LREP,KLREP+NREP,2)
!     (the shift by NREP is due to the fact that KLREP is a boson irrep)

      IJKL = 1
      DO KLREP = 1, NREP
         IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,IJREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)

           IKJL0 = h%twoNonDirac(IKREP,JLREP) & 
     &          + (g%twoDirac(JREP,LREP)+(L-1)*JFIE(JREP))*f%oneDirac(IKREP)

           IKJL0 = IKJL0 + f%twoDirac(IREP,KREP) + (K-1)*IFIE(IREP)

           ILJK0 = h%twoNonDirac(ILREP,JKREP) & 
     &          + (g%twoDirac(JREP,KREP)+(K-1)*JFIE(JREP))*f%oneDirac(ILREP)

           ILJK0 = ILJK0 + f%twoDirac(IREP,LREP) + (L-1)*IFIE(IREP)

           JKIL0 = h%twoNonDirac(JKREP,ILREP) & 
     &          + (g%twoDirac(IREP,LREP)+(L-1)*IFIE(IREP)) * f%oneDirac(JKREP)
           JKIL0 = JKIL0 + f%twoDirac(JREP,KREP) + (K-1)*JFIE(JREP)

           JLIK0 = h%twoNonDirac(JLREP,IKREP) & 
     &          + (g%twoDirac(IREP,KREP)+(K-1)*IFIE(IREP)) * f%oneDirac(JLREP)
           JLIK0 = JLIK0 + f%twoDirac(JREP,LREP) + (L-1)*JFIE(JREP)
           DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             N = IFIE(IREP) - IMIN + 1
             IKJL = IKJL0 + IMIN-1
             ILJK = ILJK0 + IMIN-1
             JKIL = JKIL0 + (IMIN-1)*f%oneDirac(JKREP)
             JLIK = JLIK0 + (IMIN-1)*f%oneDirac(JLREP)
             IF (DOINV) THEN
                CALL XCOPY (N,BUF1(IKJL*RCW+1),1,BUF2(IJKL),1)
                CALL XAXPY (N,-A1,BUF1(ILJK*RCW+1),1,BUF2(IJKL),1)
                NP = f%oneDirac(JKREP)
                CALL XAXPY(N,-A1,BUF1(JKIL*RCW+1),NP,BUF2(IJKL),1)
                NP = f%oneDirac(JLREP)
                CALL XAXPY (N,A1,BUF1(JLIK*RCW+1),NP,BUF2(IJKL),1)
             ELSE
                IF ((.NOT.FIRST).OR.(IKREP.EQ.1)) &
     &            CALL XCOPY (N,BUF1(IJKL),1,BUF2(IKJL*RCW+1),1)

                IF ((.NOT.FIRST).OR.(ILREP.EQ.1)) THEN
                  CALL XCOPY (N,BUF1(IJKL),1,BUF2(ILJK*RCW+1),1)
                  CALL XSCAL (N,-A1,BUF2(ILJK*RCW+1),1)
                ENDIF
                IF ((.NOT.FIRST).OR.(JKREP.EQ.1)) THEN
                  NP = f%oneDirac(JKREP)
                  CALL XCOPY(N,BUF1(IJKL),1,BUF2(JKIL*RCW+1),NP)
                  CALL XSCAL(N,-A1,BUF2(JKIL*RCW+1),NP)
                ENDIF
                IF ((.NOT.FIRST).OR.(JLREP.EQ.1))  THEN
                  NP = f%oneDirac(JLREP)
                  CALL XCOPY (N,BUF1(IJKL),1,BUF2(JLIK*RCW+1),NP)
                ENDIF
             ENDIF
             IJKL = IJKL + RCW * N
             IKJL0 = IKJL0 + f%oneDirac(IKREP)
             ILJK0 = ILJK0 + f%oneDirac(ILREP)
             JKIL0 = JKIL0 + 1
             JLIK0 = JLIK0 + 1
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(d,e,f,g,h)
      RETURN
      END SUBROUTINE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1TT5_exp (IRP_INPUT,NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     BUF(I>J,K>L:KLREP) <--> BUF(IK,LJ)
!     If (FIRST) take only the totally symmetric (first) irrep.
!     If (DOINV) inverse sort.
!
!---------------Last modified------------------------------------------
!
!     Author : Avijit Shee, 05/05/15
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP,IRP_INPUT, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL FIRST,DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      type(Offset) :: d,e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(d,nrep,e,f,g,h)
      call auto_symmetry_offset_triangular(d,IFIE,JFIE)
      call auto_symmetry_offset_triangular(e,KFIE,LFIE)

      call auto_symmetry_offset(f,IFIE,KFIE,.false.,.false.)
      call auto_symmetry_offset(g,LFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset(h,f%oneDirac,g%oneDirac,.true.,.true.)


      IF (.NOT.DOINV) THEN
         ROW = f%oneDirac 
         COLUMN = g%oneDirac
         IF (FIRST) THEN
            CALL XCOPY(f%oneDirac(1)*g%oneDirac(1),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
         ENDIF
      ELSE
         ROW    = d%oneNonDirac
         COLUMN = e%oneNonDirac
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP
         IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,IJREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20
           JLREP = MULTB(JREP,LREP,2)
           LJREP = MULTB(LREP,JREP,2) !a bra class symmetry product
           KJREP = MULTB(KREP,JREP,2) !a bra class symmetry product
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           LIREP = MULTB(LREP,IREP,2) !a bra class symmetry product
           KIREP = MULTB(KREP,IREP,2) !a bra class symmetry product
           JKREP = MULTB(JREP,KREP,2)

           IKLJ0 = h%twoDirac(IKREP,LJREP) & 
     &          + (g%twoDirac(LREP,JREP)+L-1)* f%oneDirac(IKREP)

           IKLJ0 = IKLJ0 + f%twoDirac(IREP,KREP) + (K-1)*IFIE(IREP)

           ILKJ0 = h%twoDirac(ILREP,KJREP) & 
     &          + (g%twoDirac(KREP,JREP)+K-1) * f%oneDirac(ILREP)

           ILKJ0 = ILKJ0 + f%twoDirac(IREP,LREP) + (L-1)*IFIE(IREP)

           JKLI0 = h%twoDirac(JKREP,LIREP) &
     &          + (g%twoDirac(LREP,IREP)+L-1) * f%oneDirac(JKREP)

           JKLI0 = JKLI0 + f%twoDirac(JREP,KREP) + (K-1)*JFIE(JREP)

           JLKI0 = h%twoDirac(JLREP,KIREP) &
     &          + (g%twoDirac(KREP,IREP)+K-1) * f%oneDirac(JLREP)
           JLKI0 = JLKI0 + f%twoDirac(JREP,LREP) + (L-1)*JFIE(JREP)

           DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             N = IFIE(IREP) - IMIN + 1
             IKLJ = IKLJ0 + IMIN-1
             ILKJ = ILKJ0 + IMIN-1
             JKLI = (JKLI0+(IMIN-1)*f%oneDirac(JKREP)*LFIE(LREP)+J-1)*RCW+1 
             JLKI = (JLKI0+(IMIN-1)*f%oneDirac(JLREP)*KFIE(KREP)+J-1)*RCW+1 
             IF (DOINV) THEN
                CALL XCOPY (N,BUF1(IKLJ*RCW+1),1,BUF2(IJKL),1)
                CALL XAXPY (N,-A1,BUF1(ILKJ*RCW+1),1,BUF2(IJKL),1)
                NP = f%oneDirac(JKREP)*LFIE(LREP)
                CALL XAXPY(N,-A1,BUF1(JKLI),NP,BUF2(IJKL),1)
                NP = f%oneDirac(JLREP)*KFIE(KREP)
                CALL XAXPY (N,A1,BUF1(JLKI),NP,BUF2(IJKL),1)
             ELSE
                IF ((.NOT.FIRST).OR.(IKREP.EQ.1)) &
     &            CALL XCOPY (N,BUF1(IJKL),1,BUF2(IKLJ*RCW+1),1)

                IF ((.NOT.FIRST).OR.(ILREP.EQ.1)) THEN
                  CALL XCOPY (N,BUF1(IJKL),1,BUF2(ILKJ*RCW+1),1)
                  CALL XSCAL (N,-A1,BUF2(ILKJ*RCW+1),1)

                ENDIF
                IF ((.NOT.FIRST).OR.(JKREP.EQ.1)) THEN
                  NP = f%oneDirac(JKREP)*LFIE(LREP)
                  CALL XCOPY(N,BUF1(IJKL),1,BUF2(JKLI),NP)
                  CALL XSCAL(N,-A1,BUF2(JKLI),NP)

                ENDIF
                IF ((.NOT.FIRST).OR.(JLREP.EQ.1))  THEN
                  NP = f%oneDirac(JLREP)*KFIE(KREP)
                  CALL XCOPY (N,BUF1(IJKL),1,BUF2(JLKI),NP)

                ENDIF
             ENDIF
             IJKL = IJKL + RCW * N
             IKLJ0 = IKLJ0 + f%oneDirac(IKREP)*LFIE(LREP)
             ILKJ0 = ILKJ0 + f%oneDirac(ILREP)*KFIE(KREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(d,e,f,g,h)
      RETURN
      END SUBROUTINE

!=====================================================================================
      SUBROUTINE SRT1ST4_exp (IRP_INPUT,NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(IK,JL:JLREP=1)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT,NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL FIRST,DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      type(Offset) :: d,e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(d,nrep,e,f,g,h)

      call auto_symmetry_offset(d,IFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset_triangular(e,KFIE,LFIE)

      call auto_symmetry_offset(f,IFIE,KFIE,.false.,.false.)
      call auto_symmetry_offset(g,JFIE,LFIE,.false.,.false.)
      call auto_symmetry_offset(h,f%oneDirac,g%oneDirac,.true.,.true.)


      IF (.NOT.DOINV) THEN
         ROW = f%oneDirac 
         COLUMN = g%oneDirac
      IF (FIRST) THEN
        CALL XCOPY(f%oneDirac(1)*g%oneDirac(1),A0,0,BUF2,1)
      ELSE
        CALL XCOPY(h%oneNonDirac(IRP_INPUT),A0,0,BUF2,1)
      ENDIF
      ELSE
         ROW    = d%oneNonDirac
         COLUMN = e%oneNonDirac
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP
        IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,IJREP+NREP,2)
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)
           IKJL0 = h%twoNonDirac(IKREP,JLREP)  & 
     &          + (g%twoDirac(JREP,LREP)+(L-1)*JFIE(JREP)) *f%oneDirac(IKREP)
           IKJL0 = IKJL0 + f%twoDirac(IREP,KREP) + (K-1)*IFIE(IREP)
           IKJL0 = IKJL0 * RCW + 1
           ILJK0 = h%twoNonDirac(ILREP,JKREP) &
     &          + (g%twoDirac(JREP,KREP)+(K-1)*JFIE(JREP)) * f%oneDirac(ILREP)
           ILJK0 = ILJK0 + f%twoDirac(IREP,LREP) + (L-1)*IFIE(IREP)
           ILJK0 = ILJK0 * RCW + 1
           DO J = 1, JFIE(JREP)
            IF (DOINV) THEN
             CALL XCOPY(IFIE(IREP),BUF1(IKJL0),1,BUF2(IJKL),1)
             CALL XAXPY(IFIE(IREP),-A1,BUF1(ILJK0),1,BUF2(IJKL),1)
            ELSE
            IF ((JLREP.EQ.1).OR.(.NOT.FIRST)) &
     &         CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKJL0),1)
            IF ((JKREP.EQ.1).OR.(.NOT.FIRST)) THEN 
               CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(ILJK0),1)
               CALL XSCAL (IFIE(IREP),-A1,BUF2(ILJK0),1)
            ENDIF
            ENDIF
            IJKL = IJKL + IFIE(IREP) * RCW
            IKJL0 = IKJL0 + f%oneDirac(IKREP) * RCW
            ILJK0 = ILJK0 + f%oneDirac(ILREP) * RCW
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(d,e,f,g,h)
      RETURN
      END SUBROUTINE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1SS4_exp (IRP_INPUT,NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,KL:KLREP) to array BUF2(IK,JL:JLREP=1)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL FIRST,DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      type(Offset) :: d,e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(d,nrep,e,f,g,h)
      call auto_symmetry_offset(d,IFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset(e,KFIE,LFIE,.false.,.false.)

      call auto_symmetry_offset(f,IFIE,KFIE,.false.,.false.)
      call auto_symmetry_offset(g,JFIE,LFIE,.false.,.false.)
      call auto_symmetry_offset(h,f%oneDirac,g%oneDirac,.true.,.true.)


      IF (.NOT.DOINV) THEN
         ROW = f%oneDirac 
         COLUMN = g%oneDirac
      IF (FIRST) THEN
        CALL XCOPY(f%oneDirac(1)*g%oneDirac(1),A0,0,BUF2,1)
      ELSE
        CALL XCOPY(h%oneNonDirac(IRP_INPUT),A0,0,BUF2,1)
      ENDIF
      ELSE
         ROW    = d%oneNonDirac
         COLUMN = e%oneNonDirac
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP
        IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,IJREP+NREP,2)
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           IKJL0 = h%twoNonDirac(IKREP,JLREP)  & 
     &          + (g%twoDirac(JREP,LREP)+(L-1)*JFIE(JREP)) * f%oneDirac(IKREP)
           IKJL0 = IKJL0 + f%twoDirac(IREP,KREP) + (K-1)*IFIE(IREP)

            DO J = 1, JFIE(JREP)
            IF (DOINV) THEN
              CALL XCOPY (IFIE(IREP),BUF1(IKJL0*RCW+1),1,BUF2(IJKL),1)
            ELSE
            IF ((JLREP.EQ.1).OR.(.NOT.FIRST))  & 
     &        CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKJL0*RCW+1),1)
            ENDIF
            IJKL = IJKL + IFIE(IREP) * RCW
            IKJL0 = IKJL0 + f%oneDirac(IKREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(d,e,f,g,h)
      RETURN
      END SUBROUTINE
!=======================================================================
      SUBROUTINE SRT1T3_exp (IRP_INPUT,NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     BUF(IJ,K>L:KLREP) <--> BUF(IJ K,L:LREP)
!     If (DOINV) inverse sort.
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      INTEGER IJKREP,IJLREP
      type(Offset) :: e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset_triangular(e,KFIE,LFIE)
      if (all(IFIE == JFIE)) then 
      call auto_symmetry_offset_triangular(f,IFIE,JFIE)
      else
      call auto_symmetry_offset(f,IFIE,JFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,f%onenonDirac,KFIE,.true.,.false.)
      call auto_symmetry_offset(h,g%oneDirac,LFIE,.false.,.false.)

      if (.not.doinv) then 
         CALL XCOPY(h%oneNonDirac(IRP_INPUT),A0,0,BUF2,1)
         row = g%oneDirac
         column = lfie
      else
         row = f%oneNonDirac 
         column = e%oneNonDirac
      endif

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       N = f%oneNonDirac(IJREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IJKREP = MULTB(IJREP+NREP,KREP,2) 
        IJLREP = MULTB(IJREP+NREP,LREP,2) 
        IF (KREP.LT.LREP) GOTO 10
        LOFF = (h%twoNonDirac(IJKREP,LREP) + g%twoDirac(IJREP,KREP)) * RCW
        KOFF = (h%twoNonDirac(IJLREP,KREP) + g%twoDirac(IJREP,LREP)) * RCW
        DO L = 1, LFIE(LREP)
         IJL  = (L-1) * N * RCW + 1
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
            IJK  = (K-1) * N * RCW + 1
            IJKL1 = LOFF + (L-1) * g%oneDirac(IJKREP) * RCW + IJK 
            IJLK1 = KOFF + (K-1) * g%oneDirac(IJLREP) * RCW + IJL 
            IF (DOINV) THEN
               CALL XCOPY (N,BUF1(IJKL1),1,BUF2(IJKL),1)
               CALL XAXPY (N,-A1,BUF1(IJLK1),1,BUF2(IJKL),1)
            ELSE
               CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJKL1),1)
               CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJLK1),1)
               CALL XSCAL (N,-A1,BUF2(IJLK1),1)
            ENDIF
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
      call dealloc_array(e,f,g,h)

      RETURN
      END SUBROUTINE
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SRT1T2_exp (IRP_INPUT, NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     BUF(I>J,KL:KLREP) <--> BUF(I,J KL:JKLREP)
!     If (DOINV) inverse sort.
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      INTEGER JKLREP, IKLREP 
      type(Offset) :: e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset_triangular(e,IFIE,JFIE)
      if (all(KFIE == LFIE)) then !this may be a rstriction since it may happen that the holes and particles are equal in number.
      call auto_symmetry_offset_triangular(f,KFIE,LFIE)
      else
      call auto_symmetry_offset(f,KFIE,LFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,JFIE,f%onenonDirac,.false.,.true.)
      call auto_symmetry_offset(h,IFIE,g%oneDirac,.false.,.false.)

      if (.not.doinv) then 
         CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
         row = jfie
         column = g%oneDirac
      else
         row = e%oneNonDirac 
         column = f%oneNonDirac
      endif

      IJKL = 0
      DO KLREP = 1, NREP
!       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       M = e%onenonDirac(IJREP)
       N = f%onenonDirac(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IF (IREP.LT.JREP) GOTO 10
        JKLREP = MULTB(JREP,KLREP+NREP,2)
        IKLREP = MULTB(IREP,KLREP+NREP,2)
        JKLOFF=(h%twoDirac(IREP,JKLREP)+ &
     &       g%twoDirac(JREP,KLREP)*IFIE(IREP))*RCW
        IKLOFF=(h%twoDirac(JREP,IKLREP)+ &
     &       g%twoDirac(IREP,KLREP)*JFIE(JREP))*RCW
        M2 = IFIE(IREP) * JFIE(JREP)
        DO J = 1, JFIE(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, IFIE(IREP)
            IJKL1 = JKLOFF + ((J-1) * IFIE(IREP) + I - 1) * RCW + 1
            JIKL1 = IKLOFF + ((I-1) * JFIE(JREP) + J - 1) * RCW + 1
            IF (DOINV) THEN
               CALL XCOPY(N,BUF1(IJKL1),M2,BUF2(IJKL+IJ),M)
               CALL XAXPY(N,-A1,BUF1(JIKL1),M2,BUF2(IJKL+IJ),M)
            ELSE
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(IJKL1),M2)
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(JIKL1),M2)
               CALL XSCAL(N,-A1,BUF2(JIKL1),M2)
            ENDIF
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO
      call dealloc_array(e,f,g,h)

      RETURN
      END SUBROUTINE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT6_exp (IRP_INPUT,NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(I,J;KL:KLREP) to array BUF2(I KL, J:JREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
     type(Offset) :: e,f,g,h 
!---------------Executable code--------------------------------------

      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset(e,IFIE,JFIE,.false.,.false.)

      if (all(KFIE == LFIE)) then !this may be a restriction since it may happen that the holes and particles are equal in number.
      call auto_symmetry_offset_triangular(f,KFIE,LFIE)
      else
      call auto_symmetry_offset(f,KFIE,LFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,IFIE,f%onenonDirac,.false.,.true.)
      call auto_symmetry_offset(h,g%oneDirac,JFIE,.false.,.false.)

      if (.not.doinv) then 
      CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
      row = g%oneDirac
      column = jfie
      else
      row = e%oneNonDirac 
      column = f%oneNonDirac
      endif


      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       M = e%onenonDirac(IJREP)
       N = f%onenonDirac(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IKLREP = MULTB(IREP,KLREP+NREP,2)
        JOFF = (h%twoDirac(IKLREP,JREP) + g%twoDirac(IREP,KLREP)) * RCW
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            IKLJ1 = JOFF + ((J-1)*g%oneDirac(IKLREP)+I-1) * RCW + 1
            IF (DOINV) THEN
               CALL XCOPY(N,BUF1(IKLJ1),IFIE(IREP),BUF2(IJKL+IJ),M)
            ELSE
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(IKLJ1),IFIE(IREP))
            ENDIF 
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO

      call dealloc_array(e,f,g,h)
      RETURN
      END SUBROUTINE
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT9_exp (IRP_INPUT,NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(I,J;KL:KLREP) to array BUF2(KL J,I:IREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher ; modified by Avijit Shee, so that now 
!                              it does only the sorting. c.c. and
!                              anti-symmetrization should be done
!                              later, if needed. 
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
     type(Offset) :: e,f,g,h 
!---------------Executable code--------------------------------------

      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset(e,IFIE,JFIE,.false.,.false.)

     if (all(KFIE == LFIE)) then !this may be a restriction since it may happen that the holes and particles are equal in number.
      call auto_symmetry_offset_triangular(f,KFIE,LFIE)
      else
      call auto_symmetry_offset(f,KFIE,LFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,f%onenonDirac,JFIE,.true.,.false.)
      call auto_symmetry_offset(h,g%oneDirac,IFIE,.false.,.false.)

      if (.not.doinv) then 
      CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
      row = g%oneDirac
      column = Ifie
      else
      row = e%oneNonDirac 
      column = f%oneNonDirac
      endif

      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       M = e%onenonDirac(IJREP)
       N = f%onenonDirac(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        KLJREP = MULTB(KLREP+NREP,JREP,2)
        IOFF = (h%twoNonDirac(KLJREP,IREP)+g%twoDirac(KLREP,JREP)) * RCW + 1
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            KLJI = IOFF+((I-1)*g%oneDirac(KLJREP)+(J-1)*N)*RCW
            IF (DOINV) THEN
               CALL XCOPY(N,BUF1(KLJI),1,BUF2(IJKL+IJ),M)
            ELSE
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(KLJI),1)
            ENDIF
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO

      call dealloc_array(e,f,g,h)
      RETURN
      END SUBROUTINE
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT10_exp (IRP_INPUT,NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(I,J;KL:KLREP) to array BUF2(KL I,J:IREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author :  Avijit Shee 
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
     type(Offset) :: e,f,g,h 
!---------------Executable code--------------------------------------

      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset(e,IFIE,JFIE,.false.,.false.)

     if (all(KFIE == LFIE)) then !this may be a restriction since it may happen that the holes and particles are equal in number.
      call auto_symmetry_offset_triangular(f,KFIE,LFIE)
      else
      call auto_symmetry_offset(f,KFIE,LFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,f%onenonDirac,IFIE,.true.,.false.)
      call auto_symmetry_offset(h,g%oneDirac,JFIE,.false.,.false.)

      if (.not.doinv) then 
      CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
      row = g%oneDirac
      column = Jfie
      else
      row = e%oneNonDirac 
      column = f%oneNonDirac
      endif

      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       M = e%onenonDirac(IJREP)
       N = f%onenonDirac(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        KLIREP = MULTB(KLREP+NREP,IREP,2)
        JOFF = (h%twoNonDirac(KLIREP,JREP)+g%twoDirac(KLREP,IREP)) * RCW + 1
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            KLIJ = JOFF+((J-1)*g%oneDirac(KLIREP)+(I-1)*N)*RCW
            IF (DOINV) THEN
               CALL XCOPY(N,BUF1(KLIJ),1,BUF2(IJKL+IJ),M)
            ELSE
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(KLIJ),1)
            ENDIF
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO

      call dealloc_array(e,f,g,h)
      RETURN
      END SUBROUTINE
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


      SUBROUTINE SRT1S2_exp (IRP_INPUT, NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     BUF(I J,KL:KLREP) <--> BUF(I,J KL:JKLREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
       type(Offset) :: e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset(e,IFIE,JFIE,.false.,.false.)
      if (all(KFIE == LFIE)) then !this may be a restriction since it may happen that the holes and particles are equal in number.
      call auto_symmetry_offset_triangular(f,KFIE,LFIE)
      else
      call auto_symmetry_offset(f,KFIE,LFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,JFIE,f%onenonDirac,.false.,.true.)
      call auto_symmetry_offset(h,IFIE,g%oneDirac,.false.,.false.)

      if (.not.doinv) then 
      CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
      row = ifie
      column = g%oneDirac
      else
      row = e%oneNonDirac 
      column = f%oneNonDirac
      endif

      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       M = e%onenonDirac(IJREP)
       N = f%onenonDirac(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        JKLREP = MULTB(JREP,KLREP+NREP,2)
        JKLOFF = (h%twoDirac(IREP,JKLREP)+g%twoDirac(JREP,KLREP)*IFIE(IREP)) * RCW
        M2 = IFIE(IREP) * JFIE(JREP)
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            IJKL1 = JKLOFF + ((J-1) * IFIE(IREP) + I - 1) * RCW + 1
            IF (DOINV) THEN
               CALL XCOPY(N,BUF1(IJKL1),M2,BUF2(IJKL+IJ),M)
            ELSE
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(IJKL1),M2)
            ENDIF
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO

      call dealloc_array(e,f,g,h)
      RETURN
      END SUBROUTINE
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT16_exp (IRP_INPUT,NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,KL:KLREP) to array BUF2(IL,KJ:KJREP)
!     If (FIRST) take only the totally symmetric (first) irrep.
!
!---------------Routines called----------------------------------------

!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL FIRST, DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
       type(Offset) :: d,e,f,g,h
!---------------Executable code--------------------------------------

! 15-8-2003 LV : Removed initialization of buf2, should not be necessary
! because we go from a square array to a square

      call alloc_array(d,nrep,e,f,g,h)
      call auto_symmetry_offset(d,IFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset(e,KFIE,LFIE,.false.,.false.)

      call auto_symmetry_offset(f,IFIE,LFIE,.false.,.false.)
      call auto_symmetry_offset(g,KFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset(h,f%oneDirac,g%oneDirac,.true.,.true.)

      IF (.NOT.DOINV) THEN
         ROW = f%oneDirac 
         COLUMN = g%oneDirac
         IF (FIRST) THEN
            CALL XCOPY(f%oneDirac(1)*g%oneDirac(1),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
         ENDIF
      ELSE
         ROW    = d%oneNonDirac
         COLUMN = e%oneNonDirac
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP

       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)

       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,IJREP+NREP,2)
           KJREP = MULTB(KREP,JREP,2)
           ILREP = MULTB(IREP,LREP,2)
           ILKJ0 = h%twoDirac(ILREP,KJREP) & 
     &          + (g%twoDirac(KREP,JREP)+K-1)*f%oneDirac(ILREP)
           ILKJ0 = ILKJ0 + f%twoDirac(IREP,LREP) + (L-1)*IFIE(IREP)
            DO J = 1, JFIE(JREP)
            IF (DOINV) THEN
               IF ((.NOT.FIRST).OR.(KJREP.EQ.1)) &
     &         CALL XCOPY(IFIE(IREP),BUF1(ILKJ0*RCW+1),1,BUF2(IJKL),1)
            ELSE
               IF ((.NOT.FIRST).OR.(KJREP.EQ.1)) &
     &         CALL XCOPY(IFIE(IREP),BUF1(IJKL),1,BUF2(ILKJ0*RCW+1),1)
            ENDIF
            IJKL = IJKL + RCW * IFIE(IREP)
            ILKJ0 = ILKJ0 + f%oneDirac(ILREP) * KFIE(KREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(d,e,f,g,h)
      RETURN
      END SUBROUTINE
!=============================================================================================================
      SUBROUTINE SRT46_exp (IRP_INPUT, NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,KL:KLREP) to array BUF2(IK,LJ:KJREP)
!     If (FIRST) take only the totally symmetric (first) irrep.
!
!---------------Routines called----------------------------------------

!---------------Last modified------------------------------------------
!
!     Author : Avijit Shee
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL FIRST, DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
       type(Offset) :: d,e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(d,nrep,e,f,g,h)
      call auto_symmetry_offset(d,IFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset(e,KFIE,LFIE,.false.,.false.)

      call auto_symmetry_offset(f,IFIE,KFIE,.false.,.false.)
      call auto_symmetry_offset(g,LFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset(h,f%oneDirac,g%oneDirac,.true.,.true.)

      IF (.NOT.DOINV) THEN
         ROW = f%oneDirac 
         COLUMN = g%oneDirac
         IF (FIRST) THEN
            CALL XCOPY(f%oneDirac(1)*g%oneDirac(1),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
         ENDIF
      ELSE
         ROW    = d%oneNonDirac
         COLUMN = e%oneNonDirac
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,IJREP+NREP,2)

           LJREP = MULTB(LREP,JREP,2)
           IKREP = MULTB(IREP,KREP,2)

           IKLJ0 = h%twoDirac(IKREP,LJREP) & 
     &          + (g%twoDirac(LREP,JREP)+L-1)* f%oneDirac(IKREP)

           IKLJ0 = IKLJ0 + f%twoDirac(IREP,KREP) + (K-1)*IFIE(IREP)

            DO J = 1, JFIE(JREP)
            IF (DOINV) THEN
               IF ((.NOT.FIRST).OR.(IKREP.EQ.1)) &
     &         CALL XCOPY(IFIE(IREP),BUF1(IKLJ0*RCW+1),1,BUF2(IJKL),1)
            ELSE
               IF ((.NOT.FIRST).OR.(IKREP.EQ.1)) &
     &         CALL XCOPY(IFIE(IREP),BUF1(IJKL),1,BUF2(IKLJ0*RCW+1),1)
            ENDIF
            IJKL = IJKL + RCW * IFIE(IREP)
            IKLJ0 = IKLJ0 + f%oneDirac(IKREP) * LFIE(LREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(d,e,f,g,h)
      RETURN

      END SUBROUTINE
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT19_exp (IRP_INPUT,NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)
 
      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,KL:KLREP) to array BUF2(IJ L,K:KREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
      type(Offset) :: e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset(e,KFIE,LFIE,.false.,.false.)
      if (all(IFIE == JFIE)) then 
      call auto_symmetry_offset_triangular(f,IFIE,JFIE)
      else
      call auto_symmetry_offset(f,IFIE,JFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,f%onenonDirac,LFIE,.true.,.false.)
      call auto_symmetry_offset(h,g%oneDirac,KFIE,.false.,.false.)

      if (.not.doinv) then 
      CALL XCOPY(h%oneNonDirac(IRP_INPUT),A0,0,BUF2,1)
      row = g%oneDirac
      column = kfie
      else
      row = f%oneNonDirac 
      column = e%oneNonDirac
      endif

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       N = f%oneNonDirac(IJREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IJLREP = MULTB(IJREP+NREP,LREP,2)
        KOFF = (h%twoNonDirac(IJLREP,KREP) + g%twoDirac(IJREP,LREP)) * RCW
        DO L = 1, LFIE(LREP)
         IJL  = (L-1) * N * RCW + 1
         DO K = 1, KFIE(KREP)
            IJLK1 = KOFF + (K-1) * g%oneDirac(IJLREP) * RCW + IJL
            IF(DOINV) THEN
            CALL XCOPY (N,BUF1(IJLK1),1,BUF2(IJKL),1)
            ELSE 
            CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJLK1),1)
            ENDIF
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
 
      call dealloc_array(e,f,g,h)
      RETURN
      END SUBROUTINE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT22_exp (IRP_INPUT,NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(L IJ,K:KREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      type(Offset) :: e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset_triangular(e,KFIE,LFIE)
      if (all(IFIE == JFIE)) then 
      call auto_symmetry_offset_triangular(f,IFIE,JFIE)
      else
      call auto_symmetry_offset(f,IFIE,JFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,LFIE,f%onenonDirac,.false.,.true.)
      call auto_symmetry_offset(h,g%oneDirac,KFIE,.false.,.false.)

      if (.not.doinv) then 
      row = g%oneDirac
      column = kfie
      else
      row = f%oneNonDirac 
      column = e%oneNonDirac
      endif

      IF (.NOT.DOINV) CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       N = f%oneNonDirac(IJREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        LIJREP = MULTB(LREP,IJREP+NREP,2)
        KIJREP = MULTB(KREP,IJREP+NREP,2) 
        IF (KREP.LT.LREP) GOTO 10
        KOFF = h%twoDirac(LIJREP,KREP) + g%twoDirac(LREP,IJREP)
        LOFF = h%twoDirac(KIJREP,LREP) + g%twoDirac(KREP,IJREP)
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
            LIJK = (KOFF + (K-1) * g%oneDirac(LIJREP) + L - 1) * RCW + 1 
            KIJL = (LOFF + (L-1) * g%oneDirac(KIJREP) + K - 1) * RCW + 1
            IF (DOINV) THEN
               CALL XCOPY (N,BUF1(LIJK),LFIE(LREP),BUF2(IJKL),1)
               CALL XAXPY (N,-A1,BUF1(KIJL),KFIE(KREP),BUF2(IJKL),1)
            ELSE
               CALL XCOPY (N,BUF1(IJKL),1,BUF2(LIJK),LFIE(LREP))
               CALL XCOPY (N,BUF1(IJKL),1,BUF2(KIJL),KFIE(KREP))
               CALL XSCAL (N,-A1,BUF2(KIJL),KFIE(KREP))
            ENDIF  
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(e,f,g,h)
      RETURN
      END SUBROUTINE
!===========================================================================
      SUBROUTINE SRT32_exp (IRP_INPUT,NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,KL:KLREP) to array BUF2(L IJ,K:KREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author :Avijit Shee, from SRT22 routine by Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
     type(Offset) :: e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset(e,KFIE,LFIE,.false.,.false.)
      if (all(IFIE == JFIE)) then 
      call auto_symmetry_offset_triangular(f,IFIE,JFIE)
      else
      call auto_symmetry_offset(f,IFIE,JFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,LFIE,f%onenonDirac,.false.,.true.)
      call auto_symmetry_offset(h,g%oneDirac,KFIE,.false.,.false.)

      if (.not.doinv) then 
      row = g%oneDirac
      column = kfie
      else
      row = f%oneNonDirac 
      column = e%oneNonDirac
      endif

      CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       N = f%oneNonDirac(IJREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        LIJREP = MULTB(LREP,IJREP+NREP,2)
        KOFF = h%twoDirac(LIJREP,KREP) + g%twoDirac(LREP,IJREP)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
            LIJK = (KOFF + (K-1) * g%oneDirac(LIJREP) + L - 1) * RCW + 1 
            IF (DOINV) THEN
               CALL XCOPY (N,BUF1(LIJK),LFIE(LREP),BUF2(IJKL),1)
            ELSE
               CALL XCOPY (N,BUF1(IJKL),1,BUF2(LIJK),LFIE(LREP))
            ENDIF  
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(e,f,g,h)
      RETURN
      END SUBROUTINE

!====================================================================================

      SUBROUTINE SRT1S3_exp (IRP_INPUT,NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     BUF(IJ,K L:KLREP) <--> BUF(IJ K,L:LREP)
!     If (DOINV) inverse sort.
!
!---------------Routines called----------------------------------------

!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT,NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      type(Offset) :: e,f,g,h
!---------------Executable code--------------------------------------


      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset(e,KFIE,LFIE,.false.,.false.)
      if (all(IFIE == JFIE)) then 
      call auto_symmetry_offset_triangular(f,IFIE,JFIE)
      else
      call auto_symmetry_offset(f,IFIE,JFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,f%onenonDirac,KFIE,.true.,.false.)
      call auto_symmetry_offset(h,g%oneDirac,LFIE,.false.,.false.)

      IF (DOINV) THEN
        ROW    = f%onenonDirac
        COLUMN = e%onenonDirac  
      ELSE
        CALL XCOPY(h%oneNonDirac(IRP_INPUT),A0,0,BUF2,1)
        ROW    = g%oneDirac
        COLUMN = lfie 
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       N = f%oneNonDirac(IJREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IJKREP = MULTB(IJREP+NREP,KREP,2)
        LOFF = (h%twoNonDirac(IJKREP,LREP) + g%twoDirac(IJREP,KREP)) * RCW
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
            IJK  = (K-1) * N * RCW + 1
            IJKL1 = LOFF + (L-1) * g%oneDirac(IJKREP) * RCW + IJK 
            IF (DOINV) THEN
                CALL XCOPY (N,BUF1(IJKL1),1,BUF2(IJKL),1)
            ELSE
                CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJKL1),1)
            ENDIF
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(e,f,g,h)
      RETURN
      END SUBROUTINE
      
!==========================================================================================
      SUBROUTINE SRT26_exp (IRP_INPUT,NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     BUF(IJ,K>L:KLREP) <--> BUF(IL,KJ:KJREP=1)
!     If (FIRST) take only the totally symmetric (first) irrep.
!     If (DOINV) inverse sort.
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL FIRST,DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
       type(Offset) :: d,e,f,g,h
       integer :: offx(nrep,nrep)
!---------------Executable code--------------------------------------


      call alloc_array(d,nrep,e,f,g,h) 

      call auto_symmetry_offset(d,IFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset_triangular(e,KFIE,LFIE)

      call auto_symmetry_offset(f,IFIE,LFIE,.false.,.false.)
      call auto_symmetry_offset(g,KFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset(h,f%oneDirac,g%oneDirac,.true.,.true.)

      IF (.NOT.DOINV) THEN
         ROW = f%oneDirac 
         COLUMN = g%oneDirac
         IF (FIRST) THEN
            CALL XCOPY(f%oneDirac(1)*g%oneDirac(1),A0,0,BUF2,1) !have to modify this for non-totally symmetric case.
         ELSE
            CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
         ENDIF
         offx = h%twoDirac 
      ELSE
         ROW    = d%oneNonDirac
         COLUMN = e%oneNonDirac
         offx = h%twoDirac 
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP
         IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,IJREP+NREP,2)
           KJREP = MULTB(KREP,JREP,2)
           LJREP = MULTB(LREP,JREP,2)
           ILREP = MULTB(IREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILKJ0 = offx(ILREP,KJREP) & 
     &          + (g%twoDirac(KREP,JREP)+K-1)*f%oneDirac(ILREP)
           ILKJ0 = ILKJ0 + f%twoDirac(IREP,LREP) + (L-1)*IFIE(IREP)
           IKLJ0 = offx(IKREP,LJREP) &
     &          + (g%twoDirac(LREP,JREP)+L-1)*f%oneDirac(IKREP)
           IKLJ0 = IKLJ0 + f%twoDirac(IREP,KREP) + (K-1)*IFIE(IREP)

            DO J = 1, JFIE(JREP)
            IF (DOINV) THEN
             CALL XCOPY(IFIE(IREP),BUF1(ILKJ0*RCW+1),1,BUF2(IJKL),1)
             CALL XAXPY(IFIE(IREP),-A1,BUF1(IKLJ0*RCW+1),1,BUF2(IJKL),1)
            ELSE
             IF ((.NOT.FIRST).OR.(KJREP.EQ.1)) &
     &          CALL XCOPY(IFIE(IREP),BUF1(IJKL),1,BUF2(ILKJ0*RCW+1),1)
             IF ((.NOT.FIRST).OR.(LJREP.EQ.1)) THEN
                CALL XCOPY(IFIE(IREP),BUF1(IJKL),1,BUF2(IKLJ0*RCW+1),1)
                CALL XSCAL(IFIE(IREP),-A1,BUF2(IKLJ0*RCW+1),1)
             ENDIF
            ENDIF
            IJKL = IJKL + RCW * IFIE(IREP)
            ILKJ0 = ILKJ0 + f%oneDirac(ILREP) * KFIE(KREP)
            IKLJ0 = IKLJ0 + f%oneDirac(IKREP) * LFIE(LREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO


      call dealloc_array(d,e,f,g,h) 

      RETURN
      END SUBROUTINE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT36_exp (IRP_INPUT,NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     BUF(I>J,KL:KLREP) <--> BUF(IL,KJ)
!     If (FIRST) take only the totally symmetric (first) irrep.
!     If (DOINV) inverse sort. 
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Avijit Shee
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL FIRST,DOINV
      INTEGER ROW(nrep),COLUMN(nrep)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------

      integer :: offx(nrep,nrep)
      type(Offset) :: d,e,f,g,h

!---------------Executable code--------------------------------------

      call alloc_array(d,nrep,e,f,g,h)
      call auto_symmetry_offset_triangular(d,IFIE,JFIE)
      call auto_symmetry_offset(e,KFIE,LFIE,.false.,.false.)

      call auto_symmetry_offset(f,IFIE,LFIE,.false.,.false.)
      call auto_symmetry_offset(g,KFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset(h,f%oneDirac,g%oneDirac,.true.,.true.)

      IF (.NOT.DOINV) THEN
         ROW = f%oneDirac 
         COLUMN = g%oneDirac
         offx = h%twoDirac
         IF (FIRST) THEN
            CALL XCOPY((f%oneDirac(1)*g%oneDirac(1)),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
         ENDIF
      ELSE
         ROW = d%oneNonDirac
         COLUMN = e%oneNonDirac
         offx = h%twoDirac
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       DO 10 LREP = 1, NREP
         KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,IJREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20

           KJREP = MULTB(KREP,JREP,2)
           JLREP = MULTB(JREP,LREP,2)
           ILREP = MULTB(IREP,LREP,2)
           KIREP = MULTB(KREP,IREP,2)

           ILKJ0 = OFFX(ILREP,KJREP) & 
     &          + (g%twoDirac(KREP,JREP)+K-1)*f%oneDirac(ILREP)
           ILKJ0 = ILKJ0 + f%twoDirac(IREP,LREP) + (L-1)*IFIE(IREP)

           JLKI0 = OFFX(JLREP,KIREP) &
     &          + (g%twoDirac(KREP,IREP)+K-1)*f%oneDirac(JLREP)
           JLKI0 = JLKI0 + f%twoDirac(JREP,LREP) + (L-1)*JFIE(JREP)
          
           DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             N = IFIE(IREP) - IMIN + 1

             ILKJ = ILKJ0 + IMIN-1             
             JLKI = JLKI0 + (IMIN-1)*f%oneDirac(JLREP)*KFIE(KREP)

             IF (DOINV) THEN
                CALL XCOPY (N,BUF1(ILKJ*RCW+1),1,BUF2(IJKL),1)
                 NP = f%oneDirac(JLREP)*KFIE(KREP)
                CALL XAXPY(N,-A1,BUF1(JLKI*RCW+1),NP,BUF2(IJKL),1)
             ELSE
             IF ((.NOT.FIRST).OR.(KJREP.EQ.1)) THEN
                CALL XCOPY (N,BUF1(IJKL),1,BUF2(ILKJ*RCW+1),1)
             ENDIF
                IF ((.NOT.FIRST).OR.(JLREP.EQ.1)) THEN
                  NP = f%oneDirac(JLREP)*KFIE(KREP)
                  CALL XCOPY(N,BUF1(IJKL),1,BUF2(JLKI*RCW+1),NP)
                  CALL XSCAL(N,-A1,BUF2(JLKI*RCW+1),NP)
                ENDIF
             ENDIF
             IJKL = IJKL + RCW * N
             ILKJ0 = ILKJ0 + f%oneDirac(ILREP)*KFIE(KREP)
             JLKI0 = JLKI0 + 1
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(d,e,f,g,h)

      RETURN
      END subroutine

!====================================================================================

      SUBROUTINE SRT56_exp (IRP_INPUT,NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     BUF(I>J,KL:KLREP) <--> BUF(IK,LJ)
!     If (FIRST) take only the totally symmetric (first) irrep.
!     If (DOINV) inverse sort. 
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Avijit Shee
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL FIRST,DOINV
      INTEGER ROW(nrep),COLUMN(nrep)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------

      integer :: offx(nrep,nrep)
      type(Offset) :: d,e,f,g,h

!---------------Executable code--------------------------------------

      call alloc_array(d,nrep,e,f,g,h)
      call auto_symmetry_offset_triangular(d,IFIE,JFIE)
      call auto_symmetry_offset(e,KFIE,LFIE,.false.,.false.)

      call auto_symmetry_offset(f,IFIE,KFIE,.false.,.false.)
      call auto_symmetry_offset(g,LFIE,JFIE,.false.,.false.)
      call auto_symmetry_offset(h,f%oneDirac,g%oneDirac,.true.,.true.)

      IF (.NOT.DOINV) THEN
         ROW = f%oneDirac 
         COLUMN = g%oneDirac
         offx = h%twoDirac
         IF (FIRST) THEN
            CALL XCOPY((f%oneDirac(1)*g%oneDirac(1)),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
         ENDIF
      ELSE
         ROW = d%oneNonDirac
         COLUMN = e%oneNonDirac
         offx = h%twoDirac
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       DO 10 LREP = 1, NREP
         KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,IJREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20

           IKREP = MULTB(IREP,KREP,2)
           LJREP = MULTB(LREP,JREP,2) !a bra class symmetry product
           JKREP = MULTB(JREP,KREP,2)
           LIREP = MULTB(LREP,IREP,2) !a bra class symmetry product

           IKLJ0 = OFFX(IKREP,LJREP) & 
     &          + (g%twoDirac(LREP,JREP)+L-1)* f%oneDirac(IKREP)

           IKLJ0 = IKLJ0 + f%twoDirac(IREP,KREP) + (K-1)*IFIE(IREP)

           JKLI0 = OFFX(JKREP,LIREP) &
     &          + (g%twoDirac(LREP,IREP)+L-1) * f%oneDirac(JKREP)

           JKLI0 = JKLI0 + f%twoDirac(JREP,KREP) + (K-1)*JFIE(JREP)
          
           DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             N = IFIE(IREP) - IMIN + 1

             IKLJ = IKLJ0 + IMIN-1

!             JKLI = (JKLI0+(IMIN-1)*f%oneDirac(JKREP)*LFIE(LREP)+J-1)*RCW+1 
             JKLI = JKLI0+(IMIN-1)*f%oneDirac(JKREP)*LFIE(LREP) 

             IF (DOINV) THEN
                CALL XCOPY (N,BUF1(IKLJ*RCW+1),1,BUF2(IJKL),1)
                NP = f%oneDirac(JKREP)*LFIE(LREP)
                CALL XAXPY(N,-A1,BUF1(JKLI*RCW+1),NP,BUF2(IJKL),1)
             ELSE

                IF ((.NOT.FIRST).OR.(LJREP.EQ.1)) &
     &            CALL XCOPY (N,BUF1(IJKL),1,BUF2(IKLJ*RCW+1),1)
                IF ((.NOT.FIRST).OR.(LIREP.EQ.1)) THEN
                  NP = f%oneDirac(JKREP)*LFIE(LREP)
                  CALL XCOPY(N,BUF1(IJKL),1,BUF2(JKLI*RCW+1),NP)
                  CALL XSCAL(N,-A1,BUF2(JKLI*RCW+1),NP)
                ENDIF
             ENDIF
             IJKL = IJKL + RCW * N
             IKLJ0 = IKLJ0 + f%oneDirac(IKREP)*LFIE(LREP)
             JKLI0 = JKLI0 + 1
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO

      call dealloc_array(d,e,f,g,h)

      RETURN

      END SUBROUTINE
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SRT1TS4_exp (IRP_INPUT,NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     BUF(I>J,KL:KLREP) <--> BUF(IK,JL)
!     If (FIRST) take only the totally symmetric (first) irrep.
!     If (DOINV) inverse sort.
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : MP modified 1TT4 by LV
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP)
      LOGICAL FIRST,DOINV
      INTEGER ROW(nrep),COLUMN(nrep)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      type(Offset) :: d,e,f,g,h
      integer :: offx(nrep,nrep)
!---------------Executable code--------------------------------------


      call alloc_array(d,nrep,e,f,g,h)

      call auto_symmetry_offset_triangular(d,IFIE,JFIE)
      call auto_symmetry_offset(e,KFIE,LFIE,.false.,.false.)

      call auto_symmetry_offset(f,IFIE,KFIE,.false.,.false.)
      call auto_symmetry_offset(g,JFIE,LFIE,.false.,.false.)
      call auto_symmetry_offset(h,f%oneDirac,g%oneDirac,.true.,.true.)

      IF (.NOT.DOINV) THEN
         ROW    = f%oneDirac  
         COLUMN = g%oneDirac  
         offx   = h%twoNonDirac 
         IF (FIRST) THEN
            CALL XCOPY(f%oneDirac(1)*g%oneDirac(1),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(h%oneNonDirac(IRP_INPUT),A0,0,BUF2,1)
         ENDIF
      ELSE
         ROW = d%oneNonDirac 
         COLUMN = e%oneNonDirac 
         offx = h%twoNonDirac 
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP
        IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,IJREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20
           IKREP = MULTB(IREP,KREP,2)
           JKREP = MULTB(JREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JLREP = MULTB(JREP,LREP,2)
           IKJL0 = offx(IKREP,JLREP) & 
     &          + (g%twoDirac(JREP,LREP)+(L-1)*JFIE(JREP)) * f%oneDirac(IKREP)
           IKJL0 = IKJL0 + f%twoDirac(IREP,KREP) + (K-1)*IFIE(IREP)
           JKIL0 = offx(JKREP,ILREP) &
     &          + (g%twoDirac(IREP,LREP)+(L-1)*IFIE(IREP)) * f%oneDirac(JKREP)
           JKIL0 = JKIL0 + f%twoDirac(JREP,KREP) + (K-1)*JFIE(JREP)
           DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             N = IFIE(IREP) - IMIN + 1
             IKJL = IKJL0 + IMIN-1
             JKIL = JKIL0 + (IMIN-1)*f%oneDirac(JKREP)
             IF (DOINV) THEN
                CALL XCOPY (N,BUF1(IKJL*RCW+1),1,BUF2(IJKL),1)
                NP = f%oneDirac(JKREP)
                CALL XAXPY(N,-A1,BUF1(JKIL*RCW+1),NP,BUF2(IJKL),1)
             ELSE
                IF ((.NOT.FIRST).OR.(IKREP.EQ.1)) THEN
                  CALL XCOPY (N,BUF1(IJKL),1,BUF2(IKJL*RCW+1),1)
                ENDIF
                IF ((.NOT.FIRST).OR.(JKREP.EQ.1)) THEN
                  NP = f%oneDirac(JKREP)
                  CALL XCOPY(N,BUF1(IJKL),1,BUF2(JKIL*RCW+1),NP)
                  CALL XSCAL(N,-A1,BUF2(JKIL*RCW+1),NP)
                ENDIF
             ENDIF
             IJKL = IJKL + RCW * N
             IKJL0 = IKJL0 + f%oneDirac(IKREP)
             JKIL0 = JKIL0 + 1
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO


      call dealloc_array(d,e,f,g,h)

      RETURN
      END SUBROUTINE

!=======================================================================

      SUBROUTINE SRT7_exp (IRP_INPUT,NREP,MULTB,DOINV,IFIE,JFIE,KFIE,LFIE,BUF1,BUF2,ROW,COLUMN)

      IMPLICIT INTEGER (A-Z)

!---------------Description--------------------------------------------
!
!     Sort array BUF1(I>J;KL:KLREP) to array BUF2(I KL, J:JREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : LV, extended for the case I>J by MP
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT, NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP) 
      LOGICAL DOINV
      INTEGER ROW(NREP),COLUMN(NREP)

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      type(Offset) :: e,f,g,h 
!---------------Executable code--------------------------------------

      call alloc_array(e,nrep,f,g,h)

      call auto_symmetry_offset_triangular(e,IFIE,JFIE)
      if (all(KFIE == LFIE)) then !this may be a restriction since it may happen that the holes and particles are equal in number.
      call auto_symmetry_offset_triangular(f,KFIE,LFIE)
      else
      call auto_symmetry_offset(f,KFIE,LFIE,.false.,.false.)
      endif

      call auto_symmetry_offset(g,IFIE,f%onenonDirac,.false.,.true.)
      call auto_symmetry_offset(h,g%oneDirac,JFIE,.false.,.false.)

      if (.not.doinv) then 
      row = g%oneDirac
      column = jfie
      else
      row = e%oneNonDirac 
      column = f%oneNonDirac
      endif

      IF (.NOT.DOINV) CALL XCOPY(h%oneDirac(IRP_INPUT),A0,0,BUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       M = e%onenonDirac(IJREP)
       N = f%onenonDirac(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IKLREP = MULTB(IREP,KLREP+NREP,2)
        JKLREP = MULTB(JREP,KLREP+NREP,2)
        IF(IREP.LT.JREP) GOTO 10
        JOFF = (h%twoDirac(IKLREP,JREP) + g%twoDirac(IREP,KLREP)) * RCW
        IOFF = (h%twoDirac(JKLREP,IREP) + g%twoDirac(JREP,KLREP)) * RCW
        DO J = 1, JFIE(JREP)
         IMIN=1
         IF(IREP.EQ.JREP) IMIN=J+1
         DO I = IMIN, IFIE(IREP)
            IKLJ1 = JOFF + ((J-1)*g%oneDirac(IKLREP)+I-1) * RCW + 1
            JKLI1 = IOFF + ((I-1)*g%oneDirac(JKLREP)+J-1) * RCW + 1
            IF (DOINV) THEN
               CALL XCOPY(N,BUF1(IKLJ1),IFIE(IREP),BUF2(IJKL+IJ),M)
               CALL XAXPY(N,-A1,BUF1(JKLI1),JFIE(JREP),BUF2(IJKL+IJ),M)
            ELSE
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(IKLJ1),IFIE(IREP))
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(JKLI1),JFIE(JREP))
               CALL XSCAL(N,-A1,BUF2(JKLI1),JFIE(JREP))
            ENDIF 
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO

      call dealloc_array(e,f,g,h)
      RETURN
      END SUBROUTINE

      SUBROUTINE SRT1C1_exp (IRP_INPUT,SHIFT_BOSON,NPAIR1,NPAIR2,BUF1,BUF2, &
    &    IRREP_PROD_TYPE)
      IMPLICIT INTEGER (A-Z)
#include "symm.inc"
!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,KL:KLREP) to array BUF2*(KL,IJ:KLREP)
!     BUF1 and BUF2 belong to bosonic type irreps 
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER NPAIR1(NREP),NPAIR2(NREP),IRP_INPUT
      LOGICAL SHIFT_BOSON
      INTEGER, OPTIONAL :: IRREP_PROD_TYPE
      INTEGER  :: IRREP_PROD_TYPE_LOC

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------

!---------------Executable code--------------------------------------

!      CALL XCOPY(dot_product(npair1,npair2),A0,0,BUF2,1)

      IRREP_PROD_TYPE_LOC = 1

      IF (PRESENT(IRREP_PROD_TYPE)) IRREP_PROD_TYPE_LOC=IRREP_PROD_TYPE

      IJKL = 1
      KLIJ = 1
      DO KLREP = 1, NREP
      SELECT CASE(IRREP_PROD_TYPE_LOC)
      CASE(1)
       IF (SHIFT_BOSON) THEN
        IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
       ELSE
        IJREP = MULTB(IRP_INPUT+NREP,KLREP,2)
!        IJREP = MULTB(KLREP,IRP_INPUT+NREP,1)
       ENDIF
      CASE(2)
       IF (SHIFT_BOSON) THEN
        IJREP = MULTB(IRP_INPUT+NREP,KLREP+NREP,2)
 !       IJREP = MULTB(KLREP+NREP,IRP_INPUT+NREP,2)
       ELSE
 !       IJREP = MULTB(KLREP,IRP_INPUT+NREP,2)
       IJREP = MULTB(IRP_INPUT+NREP,KLREP,2)
       ENDIF
      END SELECT

       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       DO IJ = 1, M
         CALL XCOPY (N,BUF1(IJKL),M,BUF2(KLIJ),1)
         IJKL = IJKL + RCW
         KLIJ = KLIJ + N * RCW
       ENDDO
       IJKL = KLIJ
      ENDDO

      RETURN
      END SUBROUTINE


     SUBROUTINE nvo_to_mvo (IRP_INPUT,DOINV,IFIE,JFIE,BUF1,BUF2)
      IMPLICIT INTEGER (A-Z)
#include "symm.inc"

!---------------Description--------------------------------------------
!
!     BUF(I>J,K>L:KLREP) <--> BUF(IK,JL)
!     If (FIRST) take only the totally symmetric (first) irrep.
!     If (DOINV) inverse sort.
!     IRP_INPUT specifies the IRREP to which the incoming tensors belong. 
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

      REAL*8 BUF1(*),BUF2(*)
      INTEGER IRP_INPUT
      INTEGER IFIE(NREP),JFIE(NREP)
      LOGICAL DOINV

!---------------Common Blocks--------------------------------------

#include "complex.inc"
#include "param.inc"

!---------------Local variables--------------------------------------
      type(Offset) :: d,e,f,g,h
!---------------Executable code--------------------------------------

      call alloc_array(d,nrep)

      call auto_symmetry_offset(d,IFIE,JFIE,.false.,.false.)

      BUF2(1:d%oneDirac(IRP_INPUT)) = 0.0d0   

      IJ = 1
      DO JREP = 1, NREP
         IREP = MULTB(JREP,IRP_INPUT+NREP,1)

           DO J = 1, JFIE(JREP)
               IJ1 = d%twoDirac(IREP,JREP)+(J-1)*IFIE(IREP)

             IF (DOINV) THEN
                CALL XCOPY (IFIE(IREP),BUF1(IJ1*RCW+1),1,BUF2(IJ),1)
             ELSE
                CALL XCOPY (IFIE(IREP),BUF1(IJ),1,BUF2(IJ1*RCW+1),1)

             ENDIF
             IJ = IJ + RCW * IFIE(IREP)
           ENDDO
       ENDDO

      call dealloc_array(d)
      RETURN
      END SUBROUTINE


  end module      
