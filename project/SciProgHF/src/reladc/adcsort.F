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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1C1N (NREP,NPAIR1,NPAIR2,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,KL:KLREP) to array BUF2(KL,IJ:KLREP)
C     Hereby no complex conjugation is performed
C
C---------------Routines called----------------------------------------
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP,NPAIR1(NREP),NPAIR2(NREP)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
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
C
      N = (KLIJ - 1) / RCW
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1L1F (NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,
     &                    OFF,OFF1,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,KL:KLREP) to array -BUF2(JI,KL:KLREP=1)
C     necessary if the array contains only contributions from the first
C     irrep ! The scaling therefore must NOT exceed the number given in
C     OFF(2) which contains the number of elements in irrep 1
C     No inverse sorting.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher, modified by MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      JIKLOFF = 0
      KLREP = 1
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
             CALL XCOPY (N,BUF1(IJKL),M,BUF2(JIKL),M)
           IJKL = IJKL + RCW
         ENDDO
        ENDDO
 10    CONTINUE
C
      CALL XSCAL (OFF(2),-A1,BUF2,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1SS4F (NREP,MULTB,FIRST,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                     OFF,OFF1,OFF2,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,KL:KLREP) to array BUF2(IK,JL:JLREP=1)
C     now it is possible to pick the first IRREP in JLREP
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : LV modified by MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(2),A0,0,BUF2,1)
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
           IKJL0 = OFF(JLREP)
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           DO J = 1, JFIE(JREP)
            IF (JLREP.EQ.1.OR.(.NOT.FIRST))
     &       CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKJL0*RCW+1),1)
            IJKL = IJKL + IFIE(IREP) * RCW
            IKJL0 = IKJL0 + NPAIR(IKREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1ST4I (NREP,MULTB,FIRST,DOINV,IFIE,JFIE,KFIE,LFIE,
     &                     NPAIR,OFF,OFF1,OFF2,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(IK,JL:JLREP=1)
C     *** ATT *** !!  addressing mode of the master offset array OFF
C     has been changed in order to reproduce the correct bosonic irrep
C     order in the resulting sorted/backsorted array !
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST,DOINV
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
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
           IKJL0 = OFF(IKREP)
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           IKJL0 = IKJL0 * RCW + 1
           ILJK0 = OFF(ILREP)
     &          + (OFF2(JREP,KREP)+(K-1)*JFIE(JREP)) * NPAIR(ILREP)
           ILJK0 = ILJK0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           ILJK0 = ILJK0 * RCW + 1
           DO J = 1, JFIE(JREP)
            IF (DOINV) THEN
               CALL XCOPY (IFIE(IREP),BUF1(IKJL0),1,BUF2(IJKL),1)
               CALL XAXPY (IFIE(IREP),-A1,BUF1(ILJK0),1,BUF2(IJKL),1)
            ELSE
              IF (JLREP.EQ.1.OR.(.NOT.FIRST))
     &           CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKJL0),1)
              IF (JKREP.EQ.1.OR.(.NOT.FIRST)) THEN
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
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1TS4 (NREP,MULTB,FIRST,DOINV,
     &                    IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                    OFF,OFF1,OFF2,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     BUF(I>J,KL:KLREP) <--> BUF(IK,JL)
C     If (FIRST) take only the totally symmetric (first) irrep.
C     If (DOINV) inverse sort.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : MP modified 1TT4 by LV
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      LOGICAL FIRST,DOINV
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (.NOT.DOINV) THEN
         IF (FIRST) THEN
            CALL XCOPY(OFF(2),A0,0,BUF2,1)
         ELSE
            CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
         ENDIF
      ENDIF
C
      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20
           IKREP = MULTB(IREP,KREP,2)
           JKREP = MULTB(JREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JLREP = MULTB(JREP,LREP,2)
           IKJL0 = OFF(IKREP) 
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           JKIL0 = OFF(JKREP) 
     &          + (OFF2(IREP,LREP)+(L-1)*IFIE(IREP)) * NPAIR(JKREP)
           JKIL0 = JKIL0 + OFF1(JREP,KREP) + (K-1)*JFIE(JREP)
           DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             N = IFIE(IREP) - IMIN + 1
             IKJL = IKJL0 + IMIN-1
             JKIL = JKIL0 + (IMIN-1)*NPAIR(JKREP)
             IF (DOINV) THEN
                CALL XCOPY (N,BUF1(IKJL*RCW+1),1,BUF2(IJKL),1)
                NP = NPAIR(JKREP)
                CALL XAXPY(N,-A1,BUF1(JKIL*RCW+1),NP,BUF2(IJKL),1)
             ELSE
                IF ((.NOT.FIRST).OR.(IKREP.EQ.1)) THEN
                  CALL XCOPY (N,BUF1(IJKL),1,BUF2(IKJL*RCW+1),1)
c           write(iw,*) IKJL,ikrep,jlrep
                ENDIF
                IF ((.NOT.FIRST).OR.(JKREP.EQ.1)) THEN
                  NP = NPAIR(JKREP)
                  CALL XCOPY(N,BUF1(IJKL),1,BUF2(JKIL*RCW+1),NP)
                  CALL XSCAL(N,-A1,BUF2(JKIL*RCW+1),NP)
c           write(iw,*) JKIL,jkrep,ilrep
                ENDIF
             ENDIF
             IJKL = IJKL + RCW * N
             IKJL0 = IKJL0 + NPAIR(IKREP)
             JKIL0 = JKIL0 + 1
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT20M (NREP,MULTB,IJPAIR,IFIE,JFIE,KFIE,LFIE,
     &                   NPAIR,OFF,OFF1,OFF2,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(IK,LJ:LJREP)
C     or the inverse direction.
C     This is a modification of the original SRT20 where
C     LJREP = 1 is selected and VOVV integrals are fetched from disk.
C     Here all irreps are sorted in the above way and buffers are
C     transferred, inverse sorting with antisymmetrization enabled.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : LV modified by MP
C
C---------------Calling variables--------------------------------------
C
      INTEGER NREP,MULTB(64,64,2),IJPAIR(NREP)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      REAL*8 BUF1(*),BUF2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 P
      PARAMETER(P=555.5555554D0)
C
C---------------Executable code--------------------------------------
C
C
      CALL XCOPY (OFF(2),A0,0,BUF2,1)

      IJKL = 1
      DO 10 KLREP = 1, NREP
       IF (IJPAIR(KLREP).EQ.0) GOTO 10
       DO 20 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 20
        DO 30 L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO 40 K = KMIN, KFIE(KREP)
          DO 50 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           LJREP = MULTB(LREP,JREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           KJREP = MULTB(KREP,JREP,2)
C          IF(OFF2(LREP,JREP).eq.0) write(iw,*) 'off2 alarm LREP/JREP'
           IKLJ0 = OFF(LJREP) + (OFF2(LREP,JREP)+L-1) * NPAIR(IKREP)
           IKLJ0 = IKLJ0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
C          IF(OFF1(IREP,KREP).eq.0) write(iw,*) 'off1 alarm IREP/KREP',K
           IKLJ0 = IKLJ0 * RCW + 1
           ILKJ0 = OFF(KJREP) + (OFF2(KREP,JREP)+K-1) * NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           ILKJ0 = ILKJ0 * RCW + 1
           DO 60 J = 1, JFIE(JREP)
             IF (LJREP.EQ.1) THEN
                CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKLJ0),1)
             ENDIF
             IF (KJREP.EQ.1) THEN
                CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(ILKJ0),1)
                CALL XSCAL (IFIE(IREP),-A1,BUF2(ILKJ0),1)
             ENDIF
C            IF(DOINV) THEN
C              CALL XCOPY (IFIE(IREP),BUF1(IKLJ0),1,BUF2(IJKL),1)
C              CALL XAXPY (IFIE(IREP),-A1,BUF1(ILKJ0),1,BUF2(IJKL),1)
C            ELSE
C              CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKLJ0),1)
C              CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(ILKJ0),1)
C              CALL XSCAL (IFIE(IREP),-A1,BUF2(ILKJ0),1)
C            ENDIF
             IJKL = IJKL + IFIE(IREP) * RCW
             IKLJ0 = IKLJ0 + NPAIR(IKREP) * LFIE(LREP) * RCW
             ILKJ0 = ILKJ0 + NPAIR(ILREP) * KFIE(KREP) * RCW
 60        CONTINUE
 50       CONTINUE
 40      CONTINUE
 30     CONTINUE
 20    CONTINUE
 10   CONTINUE
C
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT20_NTS (NREP,MULTB,IJPAIR,IFIE,JFIE,KFIE,LFIE,
     &                      NPAIR,OFF,OFF1,OFF2,BUF1,BUF2)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(IK,LJ:LJREP).
!     This routine is a modification of the original SRT20/SRT20M
!     tree where also nontotally symmetric destination ireps LJREP and
!     KJREP are included. This is necessary due to the occurrence of
!     nontotally symmetric property integrals. No inverse sort is
!     provided here since it is not needed.
!
!---------------Calling variables--------------------------------------
!
      INTEGER NREP,MULTB(64,64,2),IJPAIR(NREP)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32),OFF2(32,32)
      REAL*8 BUF1(*),BUF2(*)
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
#include "../relccsd/files.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8 P
      PARAMETER(P=555.5555554D0)
!
!---------------Executable code--------------------------------------
!
!
      CALL XCOPY (OFF(NREP+1),A0,0,BUF2,1)

      IJKL = 1
      DO 10 KLREP = 1, NREP
       IF (IJPAIR(KLREP).EQ.0) GOTO 10
       DO 20 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 20
        DO 30 L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO 40 K = KMIN, KFIE(KREP)
          DO 50 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           LJREP = MULTB(LREP,JREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           KJREP = MULTB(KREP,JREP,2)
!          IF(OFF2(LREP,JREP).eq.0) 
!    &        write(iw,*) 'off2 alarm LREP/JREP in SRT20_NTS'
           IKLJ0 = OFF(LJREP) + (OFF2(LREP,JREP)+L-1) * NPAIR(IKREP)
           IKLJ0 = IKLJ0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
!          IF(OFF1(IREP,KREP).eq.0)
!    &        write(iw,*) 'off1 alarm IREP/KREP in SRT20_NTS'
           IKLJ0 = IKLJ0 * RCW + 1
           ILKJ0 = OFF(KJREP) + (OFF2(KREP,JREP)+K-1) * NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           ILKJ0 = ILKJ0 * RCW + 1
           DO 60 J = 1, JFIE(JREP)
!
!  copy to destination for all ireps and resolve antisymmetric storage
!
             CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKLJ0),1)
             CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(ILKJ0),1)
             CALL XSCAL (IFIE(IREP),-A1,BUF2(ILKJ0),1)

             IJKL = IJKL + IFIE(IREP) * RCW
             IKLJ0 = IKLJ0 + NPAIR(IKREP) * LFIE(LREP) * RCW
             ILKJ0 = ILKJ0 + NPAIR(ILREP) * KFIE(KREP) * RCW
 60        CONTINUE
 50       CONTINUE
 40      CONTINUE
 30     CONTINUE
 20    CONTINUE
 10   CONTINUE

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT7 (NREP,MULTB,DOINV,NPAIR1,IFIE,JFIE,NTRIPL,
     &                 OFF,OFF1,BUF1,BUF2,NPAIR2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(I>J;KL:KLREP) to array BUF2(I KL, J:JREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : LV, extended for the case I>J by MP
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(64,64,2)
      INTEGER NPAIR1(NREP)
      INTEGER NPAIR2(NREP)
      INTEGER IFIE(NREP),JFIE(NREP) 
      INTEGER NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(32,32)
      LOGICAL DOINV
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C

      IF (.NOT.DOINV) CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IF(IREP.LT.JREP) GOTO 10
        JOFF = (OFF(JREP) + OFF1(IREP,KLREP)) * RCW
        IOFF = (OFF(IREP) + OFF1(JREP,KLREP)) * RCW
        DO J = 1, JFIE(JREP)
         IMIN=1
         IF(IREP.EQ.JREP) IMIN=J+1
         DO I = IMIN, IFIE(IREP)
            IKLJ1 = JOFF + ((J-1)*NTRIPL(JREP)+I-1) * RCW + 1
            JKLI1 = IOFF + ((I-1)*NTRIPL(IREP)+J-1) * RCW + 1
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

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
