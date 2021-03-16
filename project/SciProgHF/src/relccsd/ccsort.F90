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

!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
! This file contains the integral sorters.
! 4 Types of symmetry-packed integral orderings can be defined
! All integrals are assumed antisymmetrized and in Dirac notation.
!
! Order 1 (default and on file) ..... < IJ | KL> as (IJ,KL,KLREP)
! Order 2 ........................... < IJ | KL> as (I,JKL,JKLREP)
! Order 3 ........................... < IJ | KL> as (IJK,L,LREP)
! Order 4 ........................... < IJ | KL> as (IK,JL,JLREP)
!
! When I and J belong to the same class (occupied or virtual), 
! triangular packing (I>J) can be used in IJ or IJK.
!
! 9-12-97 : Introduced systematic names for most of the routines,
!           some special sorters remain.
! 1-3-01 : Note that MXREP is still hardwired, better to have a small
!          common block that just contains NREP and MULTB ?
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1TT4 (NREP,MULTB,FIRST,DOINV,                       &
     &                    IFIE,JFIE,KFIE,LFIE,NPAIR,                    &
     &                    OFF,OFF1,OFF2,BUF1,BUF2)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     BUF(I>J,K>L:KLREP) <--> BUF(IK,JL)
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST,DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer ijkl,ikjl,ikjl0,ikrep,iljk,iljk0,ilrep,imin,irep,j,jkil
      integer jkil0,jkrep,jlik,jlik0,jlrep,jrep,k,klrep,kmin,krep,l
      integer lrep,n,np
!
!---------------Executable code--------------------------------------
!
      IF (.NOT.DOINV) THEN
         IF (FIRST) THEN
            CALL XCOPY(OFF(2),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
         ENDIF
      ENDIF
!
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
!
      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)
           IKJL0 = OFF(IKREP)                                           &
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           ILJK0 = OFF(ILREP)                                           &
     &          + (OFF2(JREP,KREP)+(K-1)*JFIE(JREP)) * NPAIR(ILREP)
           ILJK0 = ILJK0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           JKIL0 = OFF(JKREP)                                           &
     &          + (OFF2(IREP,LREP)+(L-1)*IFIE(IREP)) * NPAIR(JKREP)
           JKIL0 = JKIL0 + OFF1(JREP,KREP) + (K-1)*JFIE(JREP)
           JLIK0 = OFF(JLREP)                                           &
     &          + (OFF2(IREP,KREP)+(K-1)*IFIE(IREP)) * NPAIR(JLREP)
           JLIK0 = JLIK0 + OFF1(JREP,LREP) + (L-1)*JFIE(JREP)
           DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             N = IFIE(IREP) - IMIN + 1
             IKJL = IKJL0 + IMIN-1
             ILJK = ILJK0 + IMIN-1
             JKIL = JKIL0 + (IMIN-1)*NPAIR(JKREP)
             JLIK = JLIK0 + (IMIN-1)*NPAIR(JLREP)
             IF (DOINV) THEN
                CALL XCOPY (N,BUF1(IKJL*RCW+1),1,BUF2(IJKL),1)
                CALL XAXPY (N,-A1,BUF1(ILJK*RCW+1),1,BUF2(IJKL),1)
                NP = NPAIR(JKREP)
                CALL XAXPY(N,-A1,BUF1(JKIL*RCW+1),NP,BUF2(IJKL),1)
                NP = NPAIR(JLREP)
                CALL XAXPY (N,A1,BUF1(JLIK*RCW+1),NP,BUF2(IJKL),1)
             ELSE
                IF ((.NOT.FIRST).OR.(IKREP.EQ.1))                       &
     &            CALL XCOPY (N,BUF1(IJKL),1,BUF2(IKJL*RCW+1),1)

                IF ((.NOT.FIRST).OR.(ILREP.EQ.1)) THEN
                  CALL XCOPY (N,BUF1(IJKL),1,BUF2(ILJK*RCW+1),1)
                  CALL XSCAL (N,-A1,BUF2(ILJK*RCW+1),1)
                ENDIF
                IF ((.NOT.FIRST).OR.(JKREP.EQ.1)) THEN
                  NP = NPAIR(JKREP)
                  CALL XCOPY(N,BUF1(IJKL),1,BUF2(JKIL*RCW+1),NP)
                  CALL XSCAL(N,-A1,BUF2(JKIL*RCW+1),NP)
                ENDIF
                IF ((.NOT.FIRST).OR.(JLREP.EQ.1))  THEN
                  NP = NPAIR(JLREP)
                  CALL XCOPY (N,BUF1(IJKL),1,BUF2(JLIK*RCW+1),NP)
                ENDIF
             ENDIF
             IJKL = IJKL + RCW * N
             IKJL0 = IKJL0 + NPAIR(IKREP)
             ILJK0 = ILJK0 + NPAIR(ILREP)
             JKIL0 = JKIL0 + 1
             JLIK0 = JLIK0 + 1
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      RETURN
      END

!======================================================================================================================
      SUBROUTINE SRT1TT5 (NREP,MULTB,FIRST,DOINV,                       &
     &                    IFIE,JFIE,KFIE,LFIE,NPAIR,                    &
     &                    OFF,OFF1,OFF2,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST,DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
       integer ijkl,iklj,iklj0,ikrep,ilkj,ilkj0,ilrep,imin,irep,j,jkli
       integer jkli0,jkrep,jlki,jlki0,jlrep,jrep,k,klrep,kmin,krep,l
       integer ljrep,lrep,n,np
!
!---------------Executable code--------------------------------------
!

      IF (.NOT.DOINV) THEN
         IF (FIRST) THEN
            CALL XCOPY(OFF(2),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
         ENDIF
      ENDIF
!
      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20
           JLREP = MULTB(JREP,LREP,2)
           LJREP = MULTB(LREP,JREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)
           IKLJ0 = OFF(IKREP)                                           &
     &          + (OFF2(LREP,JREP)+L-1)* NPAIR(IKREP)
           IKLJ0 = IKLJ0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)

           ILKJ0 = OFF(ILREP)                                           &
     &          + (OFF2(KREP,JREP)+K-1) * NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)

           JKLI0 = OFF(JKREP)                                           &
     &          + (OFF2(LREP,IREP)+L-1) * NPAIR(JKREP)
           JKLI0 = JKLI0 + OFF1(JREP,KREP) + (K-1)*JFIE(JREP)
           JLKI0 = OFF(JLREP)                                           &
     &          + (OFF2(KREP,IREP)+K-1) * NPAIR(JLREP)
           JLKI0 = JLKI0 + OFF1(JREP,LREP) + (L-1)*JFIE(JREP)

           DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             N = IFIE(IREP) - IMIN + 1
             IKLJ = IKLJ0 + IMIN-1
             ILKJ = ILKJ0 + IMIN-1
             JKLI = (JKLI0+(IMIN-1)*NPAIR(JKREP)*LFIE(LREP)+J-1)*RCW+1 
             JLKI = (JLKI0+(IMIN-1)*NPAIR(JLREP)*KFIE(KREP)+J-1)*RCW+1 
             IF (DOINV) THEN
                CALL XCOPY (N,BUF1(IKLJ*RCW+1),1,BUF2(IJKL),1)
                CALL XAXPY (N,-A1,BUF1(ILKJ*RCW+1),1,BUF2(IJKL),1)
                NP = NPAIR(JKREP)*LFIE(LREP)
                CALL XAXPY(N,-A1,BUF1(JKLI),NP,BUF2(IJKL),1)
                NP = NPAIR(JLREP)*KFIE(KREP)
                CALL XAXPY (N,A1,BUF1(JLKI),NP,BUF2(IJKL),1)
             ELSE
                IF ((.NOT.FIRST).OR.(IKREP.EQ.1))                       &
     &            CALL XCOPY (N,BUF1(IJKL),1,BUF2(IKLJ*RCW+1),1)

                IF ((.NOT.FIRST).OR.(ILREP.EQ.1)) THEN
                  CALL XCOPY (N,BUF1(IJKL),1,BUF2(ILKJ*RCW+1),1)
                  CALL XSCAL (N,-A1,BUF2(ILKJ*RCW+1),1)

                ENDIF
                IF ((.NOT.FIRST).OR.(JKREP.EQ.1)) THEN
                  NP = NPAIR(JKREP)*LFIE(LREP)
                  CALL XCOPY(N,BUF1(IJKL),1,BUF2(JKLI),NP)
                  CALL XSCAL(N,-A1,BUF2(JKLI),NP)

                ENDIF
                IF ((.NOT.FIRST).OR.(JLREP.EQ.1))  THEN
                  NP = NPAIR(JLREP)*KFIE(KREP)
                  CALL XCOPY (N,BUF1(IJKL),1,BUF2(JLKI),NP)

                ENDIF
             ENDIF
             IJKL = IJKL + RCW * N
             IKLJ0 = IKLJ0 + NPAIR(IKREP)*LFIE(LREP)
             ILKJ0 = ILKJ0 + NPAIR(ILREP)*KFIE(KREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SRT1TT6 (NREP,MULTB,FIRST,DOINV,                       &
     &                    IFIE,JFIE,KFIE,LFIE,NPAIR,                    &
     &                    OFF,OFF1,OFF2,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST,DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
       integer ijkl,iklj,iklj0,ikrep,ilkj,ilkj0,ilrep,imin,irep,j,jkli
       integer jkli0,jkrep,jlki,jlki0,jlrep,jrep,k,klrep,kmin,krep,l
       integer ljrep,lrep,n,np,i
!
!---------------Executable code--------------------------------------
!

      IF (.NOT.DOINV) THEN
         IF (FIRST) THEN
            CALL XCOPY(OFF(2),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
         ENDIF
      ENDIF
!
      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20
           JLREP = MULTB(JREP,LREP,2)
           LJREP = MULTB(LREP,JREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)
           IKLJ0 = OFF(IKREP)                                           &
     &          + (OFF2(LREP,JREP)+L-1)* NPAIR(IKREP)
           IKLJ0 = IKLJ0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)

           ILKJ0 = OFF(ILREP)                                           &
     &          + (OFF2(KREP,JREP)+K-1) * NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)

           JKLI0 = OFF(JKREP)                                           &
     &          + (OFF2(LREP,IREP)+L-1) * NPAIR(JKREP)
           JKLI0 = JKLI0 + OFF1(JREP,KREP) + (K-1)*JFIE(JREP)
           JLKI0 = OFF(JLREP)                                           &
     &          + (OFF2(KREP,IREP)+K-1) * NPAIR(JLREP)
           JLKI0 = JLKI0 + OFF1(JREP,LREP) + (L-1)*JFIE(JREP)

           DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
              DO I = IMIN,IFIE(IREP) 

             IKLJ = IKLJ0 + I
             ILKJ = ILKJ0 + I
             JKLI = (JKLI0+(I-1)*NPAIR(JKREP)*LFIE(LREP)+J-1)*RCW+1 
             JLKI = (JLKI0+(I-1)*NPAIR(JLREP)*KFIE(KREP)+J-1)*RCW+1 



             IF (DOINV) THEN

                N = 1
                CALL XCOPY (N,BUF1(IKLJ*RCW+1),1,BUF2(IJKL),1)
                CALL XAXPY (N,-A1,BUF1(ILKJ*RCW),1,BUF2(IJKL),1)
                NP = 1
                CALL XAXPY(N,-A1,BUF1(JKLI*RCW),NP,BUF2(IJKL),1)
                NP = 1
                CALL XAXPY (N,A1,BUF1(JLKI*RCW),NP,BUF2(IJKL),1)
             ELSE
                IF ((.NOT.FIRST).OR.(IKREP.EQ.1))                       &

     &           BUF2(IKLJ*RCW) = BUF1(IJKL)

                IF ((.NOT.FIRST).OR.(ILREP.EQ.1)) THEN
                  BUF2(ILKJ*RCW) = -BUF1(IJKL)

                ENDIF
                IF ((.NOT.FIRST).OR.(JKREP.EQ.1)) THEN

                   BUF2(JKLI) = -BUF1(IJKL)

                ENDIF
                IF ((.NOT.FIRST).OR.(JLREP.EQ.1))  THEN

                  BUF2(JLKI) = BUF1(IJKL)
                 
                ENDIF
             ENDIF
             IJKL = IJKL + 1 
           ENDDO
             IKLJ0 = IKLJ0 + NPAIR(IKREP)*LFIE(LREP)
             ILKJ0 = ILKJ0 + NPAIR(ILREP)*KFIE(KREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


      SUBROUTINE SRT1ST4 (NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,   &
     &                    NPAIR,OFF,OFF1,OFF2,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST,DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
       integer ijkl,iklj,iklj0,ikrep,ilkj,ilkj0,ilrep,imin,irep,j,jkli
       integer jkli0,jkrep,jlki,jlki0,jlrep,jrep,k,klrep,kmin,krep,l
       integer ljrep,lrep,n,np,ikjl0,iljk0
!
!---------------Executable code--------------------------------------
!
      IF (.NOT.DOINV) THEN
      IF (FIRST) THEN
        CALL XCOPY(OFF(2),A0,0,BUF2,1)
      ELSE
        CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      ENDIF
      ENDIF
      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)
           IKJL0 = OFF(IKREP)                                           &
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           IKJL0 = IKJL0 * RCW + 1
           ILJK0 = OFF(ILREP)                                           &
     &          + (OFF2(JREP,KREP)+(K-1)*JFIE(JREP)) * NPAIR(ILREP)
           ILJK0 = ILJK0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           ILJK0 = ILJK0 * RCW + 1
           DO J = 1, JFIE(JREP)
            IF (DOINV) THEN
             CALL XCOPY(IFIE(IREP),BUF1(IKJL0),1,BUF2(IJKL),1)
             CALL XAXPY(IFIE(IREP),-A1,BUF1(ILJK0),1,BUF2(IJKL),1)
            ELSE
            IF ((JLREP.EQ.1).OR.(.NOT.FIRST))                           &
     &         CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKJL0),1)
            IF ((JKREP.EQ.1).OR.(.NOT.FIRST)) THEN 
               CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(ILJK0),1)
               CALL XSCAL (IFIE(IREP),-A1,BUF2(ILJK0),1)
            ENDIF
            ENDIF
            IJKL = IJKL + IFIE(IREP) * RCW
            IKJL0 = IKJL0 + NPAIR(IKREP) * RCW
            ILJK0 = ILJK0 + NPAIR(ILREP) * RCW
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1SS4 (NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,   &
     &              NPAIR,OFF,OFF1,OFF2,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST,DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer ijkl,ikjl0,ikrep,irep,j,jlrep,jrep,k,klrep,krep,l,lrep
!
!---------------Executable code--------------------------------------
!
      IF (FIRST) THEN
        CALL XCOPY(OFF(2),A0,0,BUF2,1)
      ELSE
        CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      ENDIF

      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           IKJL0 = OFF(JLREP)                                           &
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
            DO J = 1, JFIE(JREP)
            IF (DOINV) THEN
              CALL XCOPY (IFIE(IREP),BUF1(IKJL0*RCW+1),1,BUF2(IJKL),1)
            ELSE
            IF ((JLREP.EQ.1).OR.(.NOT.FIRST))                           &
     &        CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKJL0*RCW+1),1)
            ENDIF
            IJKL = IJKL + IFIE(IREP) * RCW
            IKJL0 = IKJL0 + NPAIR(IKREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1T3 (NREP,MULTB,DOINV,NPAIR,KFIE,LFIE,NTRIPL,       &
     &                   OFF,OFF1,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR(NREP),KFIE(NREP),LFIE(NREP),NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer ijk,ijkl,ijkl1,ijl,ijlk1,ijrep,k,klrep,kmin,koff,krep,l
      integer loff,lrep,n
!
!---------------Executable code--------------------------------------
!

      IF (.NOT.DOINV) THEN
         CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      ENDIF
      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = MULTB(KLREP+NREP,1+NREP,1)
       N = NPAIR(IJREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        LOFF = (OFF(LREP) + OFF1(IJREP,KREP)) * RCW
        KOFF = (OFF(KREP) + OFF1(IJREP,LREP)) * RCW
        DO L = 1, LFIE(LREP)
         IJL  = (L-1) * N * RCW + 1
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
            IJK  = (K-1) * N * RCW + 1
            IJKL1 = LOFF + (L-1) * NTRIPL(LREP) * RCW + IJK 
            IJLK1 = KOFF + (K-1) * NTRIPL(KREP) * RCW + IJL 
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
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1T2 (NREP,MULTB,DOINV,NPAIR1,IFIE,JFIE,NPAIR2,OFF,  &
     &                   OFF2,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF2(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer i,ij,ijkl,ijkl1,ijrep,ikloff,iklrep,imin,irep,j,jikl1
      integer jkloff,jklrep,jrep,klrep,m,m2,n
!
!---------------Executable code--------------------------------------
!
      IF (.NOT.DOINV) THEN
         CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      ENDIF
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = MULTB(KLREP+NREP,1+NREP,1)
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IF (IREP.LT.JREP) GOTO 10
        JKLREP = IREP
        IKLREP = JREP
        JKLOFF = (OFF(JKLREP)+OFF2(JREP,KLREP)*IFIE(IREP)) * RCW
        IKLOFF = (OFF(IKLREP)+OFF2(IREP,KLREP)*JFIE(JREP)) * RCW
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
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT6 (NREP,MULTB,DOINV,NPAIR1,IFIE,JFIE,NPAIR2,NTRIPL, &
     &                 OFF,OFF1,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
       integer i,ij,ijkl,ijrep,iklj1,irep,j,joff,jrep,klrep,m,n
!
!---------------Executable code--------------------------------------
!

      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = MULTB(KLREP+NREP,1+NREP,1)
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        JOFF = (OFF(JREP) + OFF1(IREP,KLREP)) * RCW
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            IKLJ1 = JOFF + ((J-1)*NTRIPL(JREP)+I-1) * RCW + 1
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
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT9 (NREP,MULTB,DOINV,NPAIR1,IFIE,JFIE,NPAIR2,NTRIPL, &
     &                 OFF,OFF1,BUF1,BUF2)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     Sort array BUF1(I,J;KL:KLREP) to array BUF2(KL J,I:IREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher ; modified by Avijit Shee, so that now
!                              it can do only the sorting, c.c. and
!                              anti-symmetrization should be done
!                              later if needed, that is a separate task.
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
       integer i,ij,ijkl,ijrep,ioff,irep,j,jrep,klji,klrep,m,n
!
!---------------Executable code--------------------------------------
!

      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = MULTB(KLREP+NREP,1+NREP,1)
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IOFF = (OFF(IREP)+OFF1(KLREP,JREP)) * RCW + 1
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            KLJI = IOFF+((I-1)*NTRIPL(IREP)+(J-1)*N)*RCW
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
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1L1 (NREP,MULTB,DOINV,NPAIR1,IFIE,JFIE,NPAIR2,      &
     &                   OFF,OFF1,BUF1,BUF2)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,KL:KLREP) to array -BUF2(JI,KL:KLREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer i,ijkl,ijrep,irep,j,jikl,jikloff,jrep,klrep,m,n,nis,njs
!
!---------------Executable code--------------------------------------
!

      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)

      JIKLOFF = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJKL = OFF(KLREP) * RCW + 1
       JIKLOFF = OFF(KLREP)
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        NIS = IFIE(IREP)
        NJS = JFIE(JREP)
        DO J = 1, NJS
         DO I = 1, NIS
           JIKL = (JIKLOFF+OFF1(JREP,IREP)+(I-1)*NJS+J-1)*RCW + 1
           IF (DOINV) THEN
              CALL XCOPY (N,BUF1(JIKL),M,BUF2(IJKL),M)
           ELSE
              CALL XCOPY (N,BUF1(IJKL),M,BUF2(JIKL),M)
           ENDIF
           IJKL = IJKL + RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      CALL XSCAL (OFF(NREP+1),-A1,BUF2,1)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1S2 (NREP,MULTB,DOINV,NPAIR1,IFIE,JFIE,NPAIR2,OFF,  &
     &                   OFF2,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF2(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer i,ij,ijkl,ijkl1,ijrep,irep,j,jkloff,jklrep,jrep,klrep,m,m2
      integer n
!
!---------------Executable code--------------------------------------
!

      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        JKLREP = IREP
        JKLOFF = (OFF(JKLREP)+OFF2(JREP,KLREP)*IFIE(IREP)) * RCW
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
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1C1 (NREP,NPAIR1,NPAIR2,BUF1,BUF2)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,KL:KLREP) to array BUF2*(KL,IJ:KLREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP,NPAIR1(NREP),NPAIR2(NREP)
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer ij,ijkl,klij,klrep,m,n
!
!---------------Executable code--------------------------------------
!

      CALL XCOPY(dot_product(npair1,npair2),A0,0,BUF2,1)

      IJKL = 1
      KLIJ = 1
      DO KLREP = 1, NREP
       M = NPAIR1(KLREP)
       N = NPAIR2(KLREP)
       DO IJ = 1, M
         CALL XCOPY (N,BUF1(IJKL),M,BUF2(KLIJ),1)
         IJKL = IJKL + RCW
         KLIJ = KLIJ + N * RCW
       ENDDO
       IJKL = KLIJ
      ENDDO
!
      N = (KLIJ - 1) / RCW
      IF (CARITH) CALL CONJUGA (N,BUF2,1)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT16 (NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,     &
     &                  NPAIR,OFF,OFF1,OFF2,BUF1,BUF2)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,KL:KLREP) to array BUF2(IL,KJ:KJREP)
!     If (FIRST) take only the totally symmetric (first) irrep.
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST, DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer ijkl,ilkj0,ilrep,irep,j,jrep,k,kjrep,klrep,krep,l,lrep
!
!---------------Executable code--------------------------------------
!
! 15-8-2003 LV : Removed initialization of buf2, should not be necessary
! because we go from a square array to a square
!

      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           KJREP = MULTB(KREP,JREP,2)
           ILREP = MULTB(IREP,LREP,2)
           ILKJ0 = OFF(KJREP)                                           &
     &          + (OFF2(KREP,JREP)+K-1)*NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
            DO J = 1, JFIE(JREP)
            IF (DOINV) THEN
               IF ((.NOT.FIRST).OR.(KJREP.EQ.1))                        &
     &         CALL XCOPY(IFIE(IREP),BUF2(ILKJ0*RCW+1),1,BUF1(IJKL),1)
            ELSE
               IF ((.NOT.FIRST).OR.(KJREP.EQ.1))                        &
     &         CALL XCOPY(IFIE(IREP),BUF1(IJKL),1,BUF2(ILKJ0*RCW+1),1)
            ENDIF
            IJKL = IJKL + RCW * IFIE(IREP)
            ILKJ0 = ILKJ0 + NPAIR(ILREP) * KFIE(KREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT19(NREP,MULTB,DOINV,NPAIR,KFIE,LFIE,NTRIPL,OFF,OFF1,&
     &                  BUF1,BUF2)
 
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR(NREP),KFIE(NREP),LFIE(NREP),NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer ijkl,ijl,ijlk1,ijrep,k,klrep,koff,krep,l,lrep,n
!
!---------------Executable code--------------------------------------
!

      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       N = NPAIR(KLREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        KOFF = (OFF(KREP) + OFF1(IJREP,LREP)) * RCW
        DO L = 1, LFIE(LREP)
         IJL  = (L-1) * N * RCW + 1
         DO K = 1, KFIE(KREP)
            IJLK1 = KOFF + (K-1) * NTRIPL(KREP) * RCW + IJL
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
 
      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SRT1LS1 (NREP,MULTB,NPAIR1A,NPAIR1B,IFIE,JFIE,NPAIR2,  &
     &                    OFF,OFF1,BUF1,BUF2)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     Sort array BUF1(I>J,KL:KLREP) to array BUF2(I J,KL:KLREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR1A(NREP),NPAIR1B(NREP),NPAIR2(NREP)
      INTEGER IFIE(NREP),JFIE(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer i,ijkl1,ijkl2,ijrep,imin,irep,j,jikl2,jrep,klrep,m1,m2,n
!
!---------------Executable code--------------------------------------
!
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL1 = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       M1 = NPAIR1A(IJREP)
       M2 = NPAIR1B(IJREP)
       N = NPAIR2(KLREP)
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IF (IREP.LT.JREP) GOTO 10
        DO J = 1, JFIE(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, IFIE(IREP)
            IJKL2 = (OFF(KLREP)+OFF1(IREP,JREP)                         &
     &               +(J-1)*IFIE(IREP)+I-1)*RCW + 1
            JIKL2 = (OFF(KLREP)+OFF1(JREP,IREP)                         &
     &               +(I-1)*JFIE(JREP)+J-1)*RCW + 1
            CALL XCOPY(N,BUF1(IJKL1),M1,BUF2(IJKL2),M2)
            CALL XCOPY(N,BUF1(IJKL1),M1,BUF2(JIKL2),M2)
            CALL XSCAL(N,-A1,BUF2(JIKL2),M2)
            IJKL1 = IJKL1 + RCW
         ENDDO
        ENDDO
 10    CONTINUE
      IJKL1 = IJKL1 + M1 * (N-1) * RCW
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT22 (NREP,MULTB,DOINV,NPAIR,KFIE,LFIE,NTRIPL,OFF,OFF1&
     &                  ,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR(NREP),KFIE(NREP),LFIE(NREP),NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer ijkl,ijrep,k,kijl,klrep,kmin,koff,krep,l,lijk,loff,lrep,n 
!
!---------------Executable code--------------------------------------
!
      IF (.NOT.DOINV)                                                   &
     &  CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       N = NPAIR(KLREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        KOFF = OFF(KREP) + OFF1(LREP,IJREP)
        LOFF = OFF(LREP) + OFF1(KREP,IJREP)
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
            LIJK = (KOFF + (K-1) * NTRIPL(KREP) + L - 1) * RCW + 1 
            KIJL = (LOFF + (L-1) * NTRIPL(LREP) + K - 1) * RCW + 1
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
!
      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT32 (NREP,MULTB,DOINV,NPAIR,KFIE,LFIE,NTRIPL,OFF,OFF1&
     &                  ,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR(NREP),KFIE(NREP),LFIE(NREP),NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer ijkl,ijrep,k,klrep,koff,krep,l,lijk,lrep,n
!
!---------------Executable code--------------------------------------
!
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)

      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       N = NPAIR(KLREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        KOFF = OFF(KREP) + OFF1(LREP,IJREP)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
            LIJK = (KOFF + (K-1) * NTRIPL(KREP) + L - 1) * RCW + 1 
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
!
      RETURN
      END


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT23 (NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,NTRIPL,      &
     &                  OFF,OFF1,BUF1,BUF2,NBUF1)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     Sort array VOVV(I,J;KL:KLREP) to array BUF2(I KL, J:JREP)
!     Read first array from file, representationwise
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2),NBUF1
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER NTRIPL(NREP)
      INTEGER OFF1(32,32)
      INTEGER*8 OFF(NREP+1)
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
#include "files.inc"
!
!---------------Local variables--------------------------------------
!
      LOGICAL DONE
      INTEGER ISTART,NINT,MINT
      INTEGER*8 IKLJ1,JOFF
      integer i,ij,ijrep,irep,j,jrep,kl1,klrep,m,n,nvvvo
!
!---------------Executable code--------------------------------------
!
      NVVVO = OFF(NREP+1)
      CALL XCOPY(NVVVO,A0,0,BUF2,1)
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
!
!      We can only fetch integrals that are present on our node, ISTART will be 
!      modified by GETVOVV to represent the first value of KL that is present on
!      our node. Likewise we address the target array such that we count from this
!      first value. This offset is kept in KL
!
       ISTART = 0
       KL1 = 0
       MINT=NBUF1/(MAX0(M,1))
   1   CONTINUE
       CALL GETVOVV (IJREP,ISTART,NINT,DONE,BUF1,MINT)
       CALL DELINT ('VOVV','DKKK',BUF1,IJREP,ISTART,NINT)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        JOFF = (OFF(JREP) + OFF1(IREP,KLREP) + KL1*IFIE(IREP) ) * RCW
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            IKLJ1 = JOFF + ((J-1)*NTRIPL(JREP)+I-1) * RCW + 1
            CALL XCOPY(NINT,BUF1(IJ),M,BUF2(IKLJ1),IFIE(IREP))
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IF (.NOT.DONE) THEN
          ISTART = ISTART + NINT
          KL1 = KL1 + NINT
          GOTO 1
       ENDIF
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT33 (IRP,ISTART,NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,  &
     &              NTRIPL,OFF,OFF1,BUF2,NBUF1,OFF3,DONE)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     Sort array VOVV(I,J;KL:KLREP) to array BUF2(I KL, J:JREP)
!     Read first array from file, representationwise
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Avijit Shee
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF2(*)
      INTEGER IRP,NREP,MULTB(64,64,2),OFF3,ISTART
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER NTRIPL(NREP),NBUF1
      INTEGER OFF1(NREP,NREP)
      INTEGER OFF(NREP+1)
      LOGICAL DONE
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
#include "files.inc"
!
!---------------Local variables--------------------------------------
!
!      LOGICAL DONE
      INTEGER MINT,NINT
      INTEGER IKLJ1,JOFF
      REAL*8, ALLOCATABLE :: BUF1(:)
      integer i,ij,ijrep,irep,j,jrep,klrep,m,n,nvvvo
!
!---------------Executable code--------------------------------------

       NVVVO = OFF(NREP+1)
       buf2(1:NVVVO*rcw) = 0.0d0
       KLREP = IRP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       if (M*N*rcw.eq.0) return ! no integrals on this node, we're done
       allocate(buf1(M*N*rcw))

!
!      We can only fetch integrals that are present on our node, ISTART will be 
!      modified by GETVOVV to represent the first value of KL that is present on
!      our node. Likewise we address the target array such that we count from this
!      first value. This offset is kept in KL
!
       MINT=NBUF1/(MAX0(M,1))
       CALL GETVOVV (IJREP,ISTART,NINT,DONE,BUF1,MINT)
       CALL DELINT ('VOVV','DKKK',BUF1,IJREP,ISTART,NINT)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        JOFF = (OFF(JREP) + OFF1(IREP,KLREP) + OFF3*IFIE(IREP) ) * RCW
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            IKLJ1 = JOFF + ((J-1)*NTRIPL(JREP)+I-1) * RCW + 1
            CALL XCOPY(NINT,BUF1(IJ),M,BUF2(IKLJ1),IFIE(IREP))
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE

      deallocate(buf1)
!
      RETURN
      END

!====================================================================================
      SUBROUTINE SRT43 (IRP,ISTART,NREP,MULTB,DOINV,NPAIR1,IFIE,JFIE,   &
     &         NPAIR2,NTRIPL,OFF,OFF1,BUF1,BUF2,NBUF1,OFF3)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     Sort array BUF1(I>J;KL:KLREP) to array BUF2(I KL, J:JREP)
!     Read first array from file, representationwise
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Avijit Shee
!
!---------------Calling variables--------------------------------------
!
      INTEGER IRP,NREP,MULTB(64,64,2),ISTART
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER NTRIPL(NREP),NBUF1
      INTEGER OFF1(NREP,NREP),OFF3
      INTEGER OFF(NREP+1)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
#include "files.inc"
!
!---------------Local variables--------------------------------------
!
      LOGICAL DONE
      INTEGER MINT,NINT
      INTEGER IKLJ1,JOFF,JKLI1,IOFF
      integer i,ij,ijrep,imin,irep,j,jrep,klrep,m,n,nvvvv
!
!---------------Executable code--------------------------------------
      NVVVV = OFF(IRP+1)
      IF (.NOT.DOINV) CALL XCOPY(NVVVV,A0,0,BUF2,1)
       KLREP = IRP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
!
!      We can only fetch integrals that are present on our node, ISTART will be 
!      modified by GETVOVV to represent the first value of KL that is present on
!      our node. Likewise we address the target array such that we count from this
!      first value. This offset is kept in KL
!
       MINT=NBUF1/(MAX0(M,1))
       IF (.NOT.DOINV) THEN
       CALL GETVVVV (IJREP,ISTART,NINT,DONE,BUF1,MINT)
       CALL DELINT ('VVVV','DKKK',BUF1,IJREP,ISTART,NINT)
       ENDIF
       IJ = 1
       DO JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IF(IREP.LT.JREP) CYCLE
        JOFF = (OFF(JREP) + OFF1(IREP,KLREP) + OFF3*IFIE(IREP) ) * RCW
        IOFF = (OFF(IREP) + OFF1(JREP,KLREP) + OFF3*JFIE(JREP) ) * RCW
        DO J = 1, JFIE(JREP)
         IMIN=1
         IF(IREP.EQ.JREP) IMIN=J+1
         DO I = IMIN, IFIE(IREP)
            IKLJ1 = JOFF + ((J-1)*NTRIPL(JREP)+I-1) * RCW + 1
            JKLI1 = IOFF + ((I-1)*NTRIPL(IREP)+J-1) * RCW + 1
            IF(DOINV) THEN
            CALL XCOPY(N,BUF1(IKLJ1),IFIE(IREP),BUF2(IJ),M)
            CALL XAXPY(N,-A1,BUF1(JKLI1),JFIE(JREP),BUF2(IJ),M)
            ELSE
            CALL XCOPY(N,BUF1(IJ),M,BUF2(IKLJ1),IFIE(IREP))
            CALL XCOPY(N,BUF1(IJ),M,BUF2(JKLI1),JFIE(JREP))
            CALL XSCAL(N,-A1,BUF2(JKLI1),JFIE(JREP))
            ENDIF
            IJ = IJ + RCW
         ENDDO
        ENDDO
       ENDDO
!
      RETURN
      END

!====================================================================================

      SUBROUTINE SRT1S3 (NREP,MULTB,DOINV,NPAIR,KFIE,LFIE,NTRIPL,       &
     &                   OFF,OFF1,BUF1,BUF2)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     BUF(IJ,K L:KLREP) <--> BUF(IJ K,L:LREP)
!     If (DOINV) inverse sort.
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR(NREP),KFIE(NREP),LFIE(NREP),NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
      integer ijk,ijkl,ijkl1,ijrep,k,klrep,krep,l,loff,lrep,n
!
!---------------Executable code--------------------------------------
!
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       N = NPAIR(KLREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        LOFF = (OFF(LREP) + OFF1(IJREP,KREP)) * RCW
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
            IJK  = (K-1) * N * RCW + 1
            IJKL1 = LOFF + (L-1) * NTRIPL(LREP) * RCW + IJK 
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
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT26 (NREP,MULTB,FIRST,DOINV,                         &
     &                  IFIE,JFIE,KFIE,LFIE,NPAIR,                      &
     &                  OFF,OFF1,OFF2,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST,DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
       integer ijkl,iklj0,ijrep,ilkj0,ilrep,irep,j,jrep,k,kjrep,klrep
       integer kmin,krep,l,ljrep,lrep,ikrep
!
!---------------Executable code--------------------------------------
!
      IF (.NOT.DOINV) THEN
         IF (FIRST) THEN
            CALL XCOPY(OFF(2),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
         ENDIF
      ENDIF
!
      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           KJREP = MULTB(KREP,JREP,2)
           LJREP = MULTB(LREP,JREP,2)
           ILREP = MULTB(IREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILKJ0 = OFF(KJREP)                                           &
     &          + (OFF2(KREP,JREP)+K-1)*NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           IKLJ0 = OFF(LJREP)                                           &
     &          + (OFF2(LREP,JREP)+L-1)*NPAIR(IKREP)
           IKLJ0 = IKLJ0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
            DO J = 1, JFIE(JREP)
            IF (DOINV) THEN
             CALL XCOPY(IFIE(IREP),BUF1(ILKJ0*RCW+1),1,BUF2(IJKL),1)
             CALL XAXPY(IFIE(IREP),-A1,BUF1(IKLJ0*RCW+1),1,BUF2(IJKL),1)
            ELSE
             IF ((.NOT.FIRST).OR.(KJREP.EQ.1))                          &
     &          CALL XCOPY(IFIE(IREP),BUF1(IJKL),1,BUF2(ILKJ0*RCW+1),1)
             IF ((.NOT.FIRST).OR.(LJREP.EQ.1)) THEN
                CALL XCOPY(IFIE(IREP),BUF1(IJKL),1,BUF2(IKLJ0*RCW+1),1)
                CALL XSCAL(IFIE(IREP),-A1,BUF2(IKLJ0*RCW+1),1)
             ENDIF
            ENDIF
            IJKL = IJKL + RCW * IFIE(IREP)
            ILKJ0 = ILKJ0 + NPAIR(ILREP) * KFIE(KREP)
            IKLJ0 = IKLJ0 + NPAIR(IKREP) * LFIE(LREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT36 (NREP,MULTB,FIRST,DOINV,                         &
     &                    IFIE,JFIE,KFIE,LFIE,NPAIR,                    &
     &                    OFF,OFF1,OFF2,BUF1,BUF2)
!
      implicit none
!
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
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST,DOINV
!
!q---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
#include "files.inc"
!
!---------------Local variables--------------------------------------
!
      integer ijkl,ilkj,ilkj0,ilrep,imin,irep,j,jlki,jlki0,jlrep,k
      integer kirep,kjrep,klrep,krep,l,lrep,n,np,jrep
!
!---------------Executable code--------------------------------------
!

      IF (.NOT.DOINV) THEN
         IF (FIRST) THEN
            CALL XCOPY(OFF(2),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
         ENDIF
      ENDIF
!
      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
         KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20

           KJREP = MULTB(KREP,JREP,2)
           JLREP = MULTB(JREP,LREP,2)
           ILREP = MULTB(IREP,LREP,2)
           KIREP = MULTB(KREP,IREP,2)


           ILKJ0 = OFF(KJREP)                                           &
     &          + (OFF2(KREP,JREP)+K-1)*NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)

           JLKI0 = OFF(JLREP)                                           &
     &          + (OFF2(KREP,IREP)+K-1)*NPAIR(JLREP)
           JLKI0 = JLKI0 + OFF1(JREP,LREP) + (L-1)*JFIE(JREP)

          
           DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             N = IFIE(IREP) - IMIN + 1

             ILKJ = ILKJ0 + IMIN-1             
             JLKI = JLKI0 + (IMIN-1)*NPAIR(JLREP)*KFIE(KREP)

             IF (DOINV) THEN
                CALL XCOPY (N,BUF1(ILKJ*RCW+1),1,BUF2(IJKL),1)
                 NP = NPAIR(JLREP)*KFIE(KREP)
                CALL XAXPY(N,-A1,BUF1(JLKI*RCW+1),NP,BUF2(IJKL),1)
             ELSE
             IF ((.NOT.FIRST).OR.(KJREP.EQ.1)) THEN
                CALL XCOPY (N,BUF1(IJKL),1,BUF2(ILKJ*RCW+1),1)
             ENDIF
                IF ((.NOT.FIRST).OR.(JLREP.EQ.1)) THEN
                  NP = NPAIR(JLREP)*KFIE(KREP)
                  CALL XCOPY(N,BUF1(IJKL),1,BUF2(JLKI*RCW+1),NP)
                  CALL XSCAL(N,-A1,BUF2(JLKI*RCW+1),NP)
                ENDIF
             ENDIF
             IJKL = IJKL + RCW * N
             ILKJ0 = ILKJ0 + NPAIR(ILREP)*KFIE(KREP)
             JLKI0 = JLKI0 + 1
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

       SUBROUTINE SRT1R1 (NREP,MULTB,DOINV,NPAIR1,KFIE,LFIE,NPAIR2,     &
     &                   OFF,OFF1,BUF1,BUF2)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,KL:KLREP) to array -BUF2(IJ,LK:KLREP)
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR1(NREP),KFIE(NREP),LFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
!
!---------------Local variables--------------------------------------
!
       integer ijkl,ijlk,ijlkoff,ijrep,k,klrep,krep,l,lrep,m,n,nks,nls
!
!---------------Executable code--------------------------------------
!

      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)

      IJLKOFF = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJKL = OFF(KLREP) * RCW + 1
       IJLKOFF = OFF(KLREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        NKS = KFIE(KREP)
        NLS = LFIE(LREP)
        DO L = 1, NLS
         DO K = 1, NKS
           IJLK = (IJLKOFF+(OFF1(LREP,KREP)+(K-1)*NLS+L-1)*M)*RCW + 1
           IF (DOINV) THEN
              CALL XCOPY (M,BUF1(IJLK),1,BUF2(IJKL),1)
           ELSE
              CALL XCOPY (M,BUF1(IJKL),1,BUF2(IJLK),1)
           ENDIF
           IJKL = IJKL + M*RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      CALL XSCAL (OFF(NREP+1),-A1,BUF2,1)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT20D(NREP,MULTB,IJPAIR,IFIE,JFIE,KFIE,LFIE,          &
     &                  NPAIR,OFF,OFF1,OFF2,BUF1,NBUF1,T1,GV,USEDZ,     &
     &                  RIGHT)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(IK,LJ:LJREP=1)
!     Do contraction GV(LJ) = BUF2(IK,LJ) * T(LJ)
!     OR do contraction T(LJ) = BUF2*(IK,LJ) * GV(LJ)
!     Integrals are fetched from disk and BUF2 is never formed.
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------
!
      LOGICAL USEDZ,RIGHT
      INTEGER NBUF1
      REAL*8 BUF1(*),T1(*),GV(*)
      INTEGER NREP,MULTB(64,64,2),IJPAIR(NREP)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
#include "ccpar.inc"
!
!---------------Local variables--------------------------------------
!
      LOGICAL DONE
      REAL*8 TVAL(2)
      DATA TVAL /2*0.0D0/
! ** control variable for counting upwards to the value of KLSTAR
      INTEGER KCO,i,ijkl,ijkl1,ik,ikrep,il,ilrep,irep,j,jrep,k,kj,kj0
      integer kjrep,kl,klincr,klmax,klrep,klstar,kmin,krep,l,lj,lj0
      integer ljrep,lrep,ik0,il0
!
!---------------Executable code--------------------------------------
!
      DONE=.FALSE.
      DO 10 KLREP = 1, NREP
       IF (IJPAIR(KLREP).EQ.0) GOTO 10
       IJKL = 1
       KL = 0
       KCO = 0
! ** IJPAIR contains the length of the NVO batch
       KLMAX = NBUF1 / IJPAIR(KLREP)
       KLSTAR = 0
!
! ** KLINCR is updated on exit of GETVOVV. pure output parameter !
!
       IF (.NOT.USEDZ) THEN
          CALL GETVOVV (KLREP,KLSTAR,KLINCR,DONE,BUF1,KLMAX)
          CALL DELINT ('VOVV','KDDD',BUF1,KLREP,KLSTAR,KLINCR)
       ELSE
          CALL GETDZ   (KLREP,KLSTAR,KLINCR,DONE,BUF1)
       ENDIF
!
!      If we do a left hand contraction the integrals need to be
!      conjugated 
!
       IF (CARITH.AND..NOT.RIGHT)                                       &
     &    CALL CONJUGA (KLINCR*IJPAIR(KLREP),BUF1,1)
       DO 20 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 20
        DO 30 L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO 40 K = KMIN, KFIE(KREP)
!
          IF(KCO.LT.KLSTAR) THEN
            KCO = KCO + 1
            GOTO 40
          ENDIF
!
          KL = KL + 1
! 
          IF (KL.EQ.KLINCR+1) THEN
! -----------------------
! We need the next buffer
! -----------------------
             KLSTAR = KLSTAR + KLINCR
             KLINCR = KLMAX
             IF (.NOT.USEDZ) THEN
                IF(.NOT.DONE) THEN
                  CALL GETVOVV (KLREP,KLSTAR,KLINCR,DONE,BUF1,KLMAX)
                  CALL DELINT ('VOVV','KDDD',BUF1,KLREP,KLSTAR,KLINCR)
                ELSE
                  GOTO 10
                ENDIF
             ELSE
                CALL GETDZ   (KLREP,KLSTAR,KLINCR,DONE,BUF1)
             ENDIF
             IF (CARITH.AND..NOT.RIGHT)                                 &
     &          CALL CONJUGA (KLINCR*IJPAIR(KLREP),BUF1,1)
             KL = 1
             IJKL = 1
             KCO = KLSTAR
          ENDIF
!
          IF (RIGHT) THEN
!
             DO 50 JREP = 1, NREP
              IREP = MULTB(JREP,KLREP+NREP,2)
              LJREP = MULTB(LREP,JREP,2)
              KJREP = MULTB(KREP,JREP,2)
              LJ    = OFF2(LREP,JREP) + L - 1
              LJ    = LJ * RCW + 1
              KJ    = OFF2(KREP,JREP) + K - 1
              KJ    = KJ * RCW + 1
              IK0   = OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
              IK0   = IK0 * RCW + 1
              IL0   = OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
              IL0   = IL0 * RCW + 1
              DO 60 J = 1, JFIE(JREP)
               IF (LJREP.EQ.1) THEN
                  CALL XAXPY (IFIE(IREP),T1(LJ),BUF1(IJKL),1,GV(IK0),1)
               ENDIF
               IF (KJREP.EQ.1) THEN 
                  CALL XCOPY (1,T1(KJ),1,TVAL,1)
                  CALL XSCAL (1,-A1,TVAL,1)
                  CALL XAXPY (IFIE(IREP),TVAL,BUF1(IJKL),1,GV(IL0),1)
               ENDIF
               IJKL = IJKL + IFIE(IREP) * RCW
               LJ   = LJ   + LFIE(LREP) * RCW
               KJ   = KJ   + KFIE(KREP) * RCW
 60           CONTINUE
 50          CONTINUE
!
          ELSE
!
             DO 51 JREP = 1, NREP
              IREP = MULTB(JREP,KLREP+NREP,2)
              IKREP = MULTB(IREP,KREP,2)
              ILREP = MULTB(IREP,LREP,2)
              DO 61 I = 1, IFIE(IREP)
                 IK    = OFF1(IREP,KREP) + (K-1)*IFIE(IREP) + I - 1
                 IK    = IK * RCW + 1
                 IL    = OFF1(IREP,LREP) + (L-1)*IFIE(IREP) + I - 1
                 IL    = IL * RCW + 1
                 LJ0   = OFF2(LREP,JREP) + L - 1
                 LJ0   = LJ0 * RCW + 1
                 KJ0   = OFF2(KREP,JREP) + K - 1
                 KJ0   = KJ0 * RCW + 1
                 IJKL1 = IJKL + (I-1) * RCW 
                 IF (IKREP.EQ.1) THEN
                  CALL XAXPY (JFIE(JREP),GV(IK),BUF1(IJKL1),IFIE(IREP), &
     &                        T1(LJ0),LFIE(LREP))
                 ENDIF
                 IF (ILREP.EQ.1) THEN 
                  CALL XCOPY (1,GV(IL),1,TVAL,1)
                  CALL XSCAL (1,-A1,TVAL,1)
                  CALL XAXPY (JFIE(JREP),TVAL,BUF1(IJKL1),IFIE(IREP),   &
     &                        T1(KJ0),KFIE(KREP))
                 ENDIF
 61           CONTINUE
              IJKL = IJKL + IFIE(IREP) * JFIE(JREP) * RCW
 51          CONTINUE
!
          ENDIF
!
 40      CONTINUE
 30     CONTINUE
 20    CONTINUE
 10   CONTINUE
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRTDIR(NREP,MULTB,NPAIR,OFF,OFF1,OFF2,BUF1,BUF2)
!
      implicit none
!
!---------------Description--------------------------------------------
!
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Joost van Stralen
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP,MULTB(64,64,2)
      INTEGER NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      INTEGER INFOT(6,4,4,2)
!
!---------------Common Blocks--------------------------------------
!
#include "complex.inc"
#include "param.inc"
#include "ccpar.inc"
#include "dcborb.h"
#include "dgroup.h"
!
!---------------Local variables--------------------------------------
!
       integer c,c1,c2,cirep1,cirep2,crep,crep1,crep2,d,d1t,d2t,drep
       integer drep1,drep2,i,i1,i2,icidk,iq1,iq2,irep,irep1,irep2,irepc
       integer irepci,irepd,irepi,irepk,ireprs,itdir,iz1,iz2,k,k1,k2
       integer krep,krep1,krep2
!
!---------------Executable code--------------------------------------
!
!
!     Initialize the array with reordered tau's
!
      CALL INIT_T1(BUF2,.TRUE.)
!
!     Information needed to convert the relccsd tau's to the moltra 
!     T's                           
!
      CALL TAUINFO(INFOT)

!
      ITDIR = 1 ! the pointer for the T needed in moltra
!
      DO IREPCI = 1, NFSYM
        DO IREPRS = 0, NBSYM -1
          IF(IREPCI.EQ.JBTOF(IREPRS,1)) THEN
            DO IREPI = 1, NFSYM
              IREP = 2*IREPI - 1
              IREPC = MOD(IREPI + IREPCI,2) + 1
              CREP = 2*IREPC - 1
              DO IZ2 = 1, NZ
                IQ2 = IPQTOQ(IZ2,IREPRS)
                DO I = 1, NAOCC(IREPI)
                  DO C = 1, NAVIR(IREPC)
                    DO IZ1 = 1, NZ
                      IQ1 = IPQTOQ(IZ1,IREPRS)
                      IF(NZ.LT.4) THEN
                        CREP1 = CREP + INFOT(1,IQ1,IQ2,1)
                        IREP1 = IREP + INFOT(2,IQ1,IQ2,1)
                        CREP2 = CREP + INFOT(1,IQ1,IQ2,2)
                        IREP2 = IREP + INFOT(2,IQ1,IQ2,2)
                      ELSE
                         C1 = C + INFOT(1,IQ1,IQ2,1)*NAVIR(IREPC)
                         I1 = I + INFOT(2,IQ1,IQ2,1)*NAOCC(IREPI)
                         C2 = C + INFOT(1,IQ1,IQ2,2)*NAVIR(IREPC)
                         I2 = I + INFOT(2,IQ1,IQ2,2)*NAOCC(IREPI)
                      ENDIF
                      DO IREPD = 1, NFSYM
                        IREPK = MOD(IREPD + IREPCI,2) + 1
!
                        IF(NZ.LT.4) THEN
!
                          DREP = 2*IREPD - 1
                          KREP = 2*IREPK - 1
                          DREP1 = DREP + INFOT(3,IQ1,IQ2,1)
                          KREP1 = KREP + INFOT(4,IQ1,IQ2,1)
                          CIREP1 = MULTB(CREP1,IREP1,2)
                          DREP2 = DREP + INFOT(3,IQ1,IQ2,2)
                          KREP2 = KREP + INFOT(4,IQ1,IQ2,2)
                          CIREP2 = MULTB(CREP2,IREP2,2)
!
!
!     Real and complex groups
!     =======================
!
                          DO D = 1, NAVIR(IREPD)
                            DO K = 1, NAOCC(IREPK)
!           The first tau
                            ICIDK = OFF(CIREP1) + (OFF2(DREP1,KREP1) +  &
     &                              (K-1)*NAVIR(IREPD))*NPAIR(CIREP1)   &
     &                            + (D-1)*NPAIR(CIREP1)                 &
     &                            + OFF1(CREP1,IREP1)                   &
     &                            + (I-1)*NAVIR(IREPC) + C - 1
                            IF(INFOT(6,IQ1,IQ2,1).EQ.1) THEN
      BUF2(ITDIR) = BUF2(ITDIR) - BUF1(ICIDK*RCW+1+INFOT(5,IQ1,IQ2,1))
                            ELSE
      BUF2(ITDIR) = BUF2(ITDIR) + BUF1(ICIDK*RCW+1+INFOT(5,IQ1,IQ2,1))
                            ENDIF
!           The second tau
                            ICIDK = OFF(CIREP2) + (OFF2(DREP2,KREP2) +  &
     &                              (K-1)*NAVIR(IREPD))*NPAIR(CIREP2)   &
     &                            + (D-1)*NPAIR(CIREP2)                 &
     &                            + OFF1(CREP2,IREP2)                   &
     &                            + (I-1)*NAVIR(IREPC) + C - 1
                            IF(INFOT(6,IQ1,IQ2,2).EQ.1) THEN
      BUF2(ITDIR) = BUF2(ITDIR) - BUF1(ICIDK*RCW+1+INFOT(5,IQ1,IQ2,1))
                            ELSE
      BUF2(ITDIR) = BUF2(ITDIR) + BUF1(ICIDK*RCW+1+INFOT(5,IQ1,IQ2,1))
                            ENDIF
                            ITDIR = ITDIR + 1
                            ENDDO !k
                          ENDDO !d
!
                        ELSE !nz
!
!     Quaternion groups
!     =================
!
                          DO D = 1, NAVIR(IREPD)
                            D1T = D + INFOT(3,IQ1,IQ2,1)*NAVIR(IREPD)
                            D2T = D + INFOT(3,IQ1,IQ2,2)*NAVIR(IREPD)
                            DO K = 1, NAOCC(IREPK)
!           The first tau
                            K1 = K + INFOT(4,IQ1,IQ2,1)*NAOCC(IREPK)
                            ICIDK = OFF(IREPCI) + (OFF2(IREPD,IREPK) +  &
     &                              (K1-1)*2*NAVIR(IREPD))*NPAIR(IREPCI)&
     &                            + (D1T-1)*NPAIR(IREPCI)               &
     &                            + OFF1(IREPC,IREPI)                   &
     &                            + (I1-1)*2*NAVIR(IREPC) + C1 -1
                            IF(INFOT(6,IQ1,IQ2,1).EQ.1) THEN
      BUF2(ITDIR) = BUF2(ITDIR) - BUF1(ICIDK*RCW+1+INFOT(5,IQ1,IQ2,1))
                            ELSE
      BUF2(ITDIR) = BUF2(ITDIR) + BUF1(ICIDK*RCW+1+INFOT(5,IQ1,IQ2,1))
                            ENDIF
!           The second tau
                            K2 = K + INFOT(4,IQ1,IQ2,2)*NAOCC(IREPK)
                            ICIDK = OFF(IREPCI) + (OFF2(IREPD,IREPK) +  &
     &                              (K2-1)*2*NAVIR(IREPD))*NPAIR(IREPCI)&
     &                            + (D2T-1)*NPAIR(IREPCI)               &
     &                            + OFF1(IREPC,IREPI)                   &
     &                            + (I2-1)*2*NAVIR(IREPC) + C2 - 1
                            IF(INFOT(6,IQ1,IQ2,2).EQ.1) THEN
      BUF2(ITDIR) = BUF2(ITDIR) - BUF1(ICIDK*RCW+1+INFOT(5,IQ1,IQ2,1))
                            ELSE
      BUF2(ITDIR) = BUF2(ITDIR) + BUF1(ICIDK*RCW+1+INFOT(5,IQ1,IQ2,1))
                            ENDIF
                            ITDIR = ITDIR + 1
!
                            ENDDO !k
                          ENDDO !d
!
                        ENDIF !nz
!
                      ENDDO !irepd
                    ENDDO !iz1
                  ENDDO !c
                ENDDO !i
              ENDDO !iz2
            ENDDO !irepi
          ENDIF
        ENDDO !ireprs
      ENDDO !irepci
!
      RETURN                
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

