      SUBROUTINE BFLOPS(CPU4V,CPUH1,CPUH2,NIT)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculate performance of largest contraction
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 CPU4V,CPUH1,CPUH2
      INTEGER NIT
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 MFLOPS,DRCW2,DRCW,DNVVT,DNOOT,DNVO,DNVOOT,DNV,DNO
      REAL*8 MFLOPB,MFLOPH,DMVO
      CHARACTER*30 CONTR
      integer i,irp
C
C---------------Executable code--------------------------------------
C
C     Prevent division by zero in case somebody could not figure out
C     how the cpu timer works on the system.
C
      IF (CPU4V.EQ.0.0.OR.CPUH1.EQ.0.0.OR.CPUH2.EQ.0.0) RETURN
      MFLOPB = 0.0
      MFLOPH = 0.0
C
      DRCW2 = DBLE(2*RCW**2)
      DRCW  = DBLE(RCW)
      DO IRP = 1, NREP

         DNVVT = DBLE(NVVT(IRP))
         DNOOT = DBLE(NOOT(IRP))
         DNVO  = DBLE(NVO(IRP))
         DNVOOT= DBLE(NVOOT(IRP))
         DNV   = DBLE(NV(IRP))
         DNO   = DBLE(NO(IRP))
         DMVO  = DBLE(MVO(IRP))

C        Compute number of operations for BINTM
         MFLOPB = MFLOPB + DNVVT  * DNVVT * DNOOT * DRCW2
         MFLOPB = MFLOPB + DNVVT  * DNVO  * DNOOT * DRCW2
         MFLOPB = MFLOPB + DNVOOT * DNV   * DNO   * DRCW2
         MFLOPB = MFLOPB + DNVVT  * DNOOT * DRCW

C        Compute number of operations for HINTM
         MFLOPH = MFLOPH + DMVO  * DMVO * DMVO * DRCW2

      ENDDO
C
C     Convert to MFlops
C
      MFLOPB = MFLOPB * 0.000001
      MFLOPH = MFLOPH * 0.000001
C
      WRITE (IW,1000)
C
      DO I = 1, 3
         IF (I.EQ.1) THEN
            MFLOPS = MFLOPB * DBLE(NIT) / CPU4V
            CONTR = 'VVVV+VOVV (in B: includes I/O)'
         ELSEIF (I.EQ.2) THEN
            MFLOPS = MFLOPH * DBLE(NIT) / CPUH1
            CONTR = 'VOVO (in H: only XGEMM)       '
         ELSE
            MFLOPS = MFLOPH * DBLE(NIT) / CPUH2
            CONTR = 'VOVO (in T2EQN: includes sort)'
         ENDIF
C
         IF (MFLOPS.GE.1.0E6) THEN
            WRITE (IW,1001) CONTR,MFLOPS/1.0E6,'T'
         ELSEIF (MFLOPS.GE.1.0E5) THEN
            WRITE (IW,1003) CONTR,MFLOPS/1.0E3,'G'
         ELSEIF (MFLOPS.GE.1.0E4) THEN
            WRITE (IW,1002) CONTR,MFLOPS/1.0E3,'G'
         ELSEIF (MFLOPS.GE.1.0E3) THEN
            WRITE (IW,1001) CONTR,MFLOPS/1.0E3,'G'
         ELSEIF (MFLOPS.GE.1.0E2) THEN
            WRITE (IW,1003) CONTR,MFLOPS,'M'
         ELSEIF (MFLOPS.GE.1.0E1) THEN
            WRITE (IW,1002) CONTR,MFLOPS,'M'
         ELSE
            WRITE (IW,1001) CONTR,MFLOPS,'M'
         ENDIF
C
      ENDDO
C
 1000 FORMAT (//' Performance of BLAS GEMM in the largest contractions',
     &//' Contraction type',T45,'Performance')
 1001 FORMAT (1X,A,T40,F10.3,3X,A1,'flop/s')
 1002 FORMAT (1X,A,T40,F10.2,3X,A1,'flop/s')
 1003 FORMAT (1X,A,T40,F10.1,3X,A1,'flop/s')
      RETURN
      END
      SUBROUTINE CCENG(T1,T2,FVO,BUF1,BUF2,ECC,ECCIM)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates CC energy. Checks result (should be real).
C
C---------------Routines called----------------------------------------
C
C     BLAS routines
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C              Miro Ilias,July 2007 - added imaginary part of CCSD energy
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
      REAL*8 FVO(*)
      REAL*8 BUF1(*),BUF2(*)
      REAL*8 ECC,ECCIM
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "files.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CECC1(2),CECC2(2)
C
C---------------Executable code--------------------------------------
C

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface


      CALL XCOPY (NDIMT1,FVO,1,BUF1,1)
      CALL DELFCK ('VO','DD',BUF1) ! Only used in Fock space calculations
C
      CALL XDOTC(CECC1,NFVO,BUF1,1,T1,1) ! Non-zero only in non-canonical case
      ECC = CECC1(1)
      IF (CARITH) THEN
         IF (ABS(CECC1(2)).GT.ACCUR)
     &      WRITE (IW,1000) 'T1',CECC1(2)
      ENDIF
C
      CALL GETVVOO (BUF1)
      CALL DELINT ('VVOO','DDDD',BUF1,0,0,0) ! Only used in Fock space calculations
      CALL GETTAU  (T1,T2,BUF2)
C
      CALL XDOTC(CECC2,NDIMT2,BUF1,1,BUF2,1)
      ECC = ECC + CECC2(1)
      IF (CARITH) THEN
         IF (ABS(CECC2(2)).GT.ACCUR)
     &      WRITE (IW,1000) 'T2',CECC2(2)
      ENDIF

      IF (CARITH) THEN
        ECCIM = CECC1(2) + CECC2(2)
      ENDIF
C
      RETURN
1000  FORMAT (/' WARNING : Imaginary contribution to ',A2,
     &' part of CC energy :',F20.15)
      END
      SUBROUTINE DELFCK (CLASS,DODEL,FXX)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Delete integrals that contribute to unphysical excitations
C     in the Fock space calculations.
C
C     The spinors space is arranged as
C
C
C      #       O         V
C     -----------------------
C      1      IO        AV
C      2      AO        IV
C      3      AV        AO
C
C      This routine deletes the AV integrals in the occupied space if
C      the integrals is to be used in a contraction over truly O-type
C      spinor, while in the virtual space it will delete the AO
C      class, if desired. The parameter CLASS specifies the type
C      of the incoming integrals, the parameter DODEL specifies the
C      indices that are to be cleaned.
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FXX(*)
      CHARACTER*2 CLASS,DODEL
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
      integer ic,idel,ij,imin,index,irep,j,jc,jrep,lngth,nindx
      DIMENSION NINDX(0:2,MXREP,2)
C
C---------------Executable code--------------------------------------
C
C     Return when we do not do Fock space calculations
C
      IF (EQNS.NE.'FOCKSP') RETURN
C
C     NINDX(1,IREP,INDEX) : Spinors to be kept
C     NINDX(2,IREP,INDEX) : Additional spinors to be zeroed.
C
      DO INDEX = 1, 2
         IF (CLASS(INDEX:INDEX).EQ.'O') THEN
            DO IREP = 1, NREP
               IF (DODEL(INDEX:INDEX).EQ.'D') THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSO(IREP)
                  NINDX(2,IREP,INDEX) = NO(IREP)
               ELSE
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NO(IREP)
                  NINDX(2,IREP,INDEX) = NO(IREP)
               ENDIF
            ENDDO
         ELSEIF (CLASS(INDEX:INDEX).EQ.'V') THEN
            DO IREP = 1, NREP
               IF (DODEL(INDEX:INDEX).EQ.'D') THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSV(IREP)
                  NINDX(2,IREP,INDEX) = NV(IREP)
               ELSE
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NV(IREP)
                  NINDX(2,IREP,INDEX) = NV(IREP)
               ENDIF
            ENDDO
         ELSE
            CALL QUIT ('DELFCK integral class ')
         ENDIF
      ENDDO
C
      IJ = 1
C
C     Loop over J : irrep, spinor class, spinors
C
      DO 10 JREP = 1, NREP
         IREP = JREP
         DO JC =  1, 2
         DO J =  NINDX(JC-1,JREP,2)+1, NINDX(JC,JREP,2)
C
C    Loop over I : (irrep fixed via JREP), spinor class
C    The loop over spinors is done within XCOPY
C
            DO IC = 1, 2
               IMIN = NINDX(IC-1,IREP,1)+1
               LNGTH=NINDX(IC,IREP,1)-IMIN+1
               IDEL = IC + JC
C
C              Only when IDEL = 2 all indices are allowed
C
               IF (LNGTH.GT.0.AND.IDEL.GT.2)
     &            CALL XCOPY (LNGTH,A0,0,FXX(IJ),1)
               IF (LNGTH.GT.0) IJ = IJ + RCW * LNGTH
            ENDDO    ! IC
C
         ENDDO    ! J
         ENDDO    ! JC
 10   CONTINUE ! JREP
C
      RETURN
      END
      SUBROUTINE DELINT (CLASS,DODEL,VXXXX,IREPVV,ISTART,NINT)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Delete integrals that contribute to unphysical excitations
C     in the Fock space calculations.
C
C     The spinors space is arranged as
C
C
C      #       O         V
C     -----------------------
C      1      IO        AV
C      2      AO        IV
C      3      AV        AO
C
C      This routine deletes the AV integrals in the occupied space if
C      the integrals is to be used in a contraction over truly O-type
C      spinor, while in the virtual space it will delete the AO
C      class, if desired. The parameter CLASS specifies the type
C      of the incoming integrals, the parameter DODEL specifies the
C      indices that are to be cleaned.
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VXXXX(*)
      CHARACTER*4 CLASS,DODEL
      INTEGER IREPVV,ISTART,NINT
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
      DIMENSION NINDX(0:2,MXREP,4)
      LOGICAL TRIANG(2),ALLOWED,ALLOW1,ALLOW2
      integer ic,ijkl,imin,index,indexe,indexp,indexs,irep,j,jc,jrep,k
      integer kc,kl,klrep,kmin,krep,l,lc,lngth,lrep,nindx
C
C---------------Executable code--------------------------------------
C
C     Return when we do not do Fock space calculations
C
      IF (EQNS.NE.'FOCKSP') RETURN
C
      TRIANG(1) = CLASS(1:1).EQ.CLASS(2:2)
      TRIANG(2) = CLASS(3:3).EQ.CLASS(4:4)
C
C     NINDX(1,IREP,INDEX) : Spinors to be kept
C     NINDX(2,IREP,INDEX) : Additional spinors to be zeroed.
C
      INDEXS = 1
      INDEXE = 2
      DO INDEXP = 1, 2
      DO INDEX = INDEXS, INDEXE
         IF (CLASS(INDEX:INDEX).EQ.'O') THEN
            DO IREP = 1, NREP
               IF (DODEL(INDEX:INDEX).EQ.'D') THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSO(IREP)
                  NINDX(2,IREP,INDEX) = NO(IREP)
               ELSE
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NO(IREP)
                  NINDX(2,IREP,INDEX) = NO(IREP)
               ENDIF
C
C              Triangular arrays need special treatment, we can have the
C              case that the to be deleted index is related by anti-
C              symmetry to the index to be kept.
C              We treat this as the condition 'DD'to have the proper
C              loop structure but modify the delete condition somewhat
C
               IF (TRIANG(INDEXP).AND.
     &             (DODEL(INDEXS:INDEXE).EQ.'DK'.OR.
     &              DODEL(INDEXS:INDEXE).EQ.'KD')) THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSO(IREP)
                  NINDX(2,IREP,INDEX) = NO(IREP)
               ENDIF
            ENDDO
         ELSEIF (CLASS(INDEX:INDEX).EQ.'V') THEN
            DO IREP = 1, NREP
               IF (DODEL(INDEX:INDEX).EQ.'D') THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSV(IREP)
                  NINDX(2,IREP,INDEX) = NV(IREP)
               ELSE
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NV(IREP)
                  NINDX(2,IREP,INDEX) = NV(IREP)
               ENDIF
C
C              The same treatment for the triangular condition as above
C
               IF (TRIANG(INDEXP).AND.
     &             (DODEL(INDEXS:INDEXE).EQ.'DK'.OR.
     &              DODEL(INDEXS:INDEXE).EQ.'KD')) THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSV(IREP)
                  NINDX(2,IREP,INDEX) = NV(IREP)
               ENDIF
            ENDDO
         ELSE
            CALL QUIT ('DELINT integral class ')
         ENDIF
      ENDDO
      INDEXS = INDEXS + 2
      INDEXE = INDEXE + 2
      ENDDO
C
      IJKL = 1
C
C     Loop over compound irreps
C
      DO 40 KLREP = 1, NREP
         KL = 0
C
C        Loop over L : irrep, spinor classes, spinors
C
         DO 30 LREP = 1, NREP
         KREP = MULTB(LREP,KLREP+NREP,2)
         IF (TRIANG(2).AND.KREP.LT.LREP) GOTO 30

         DO LC = 1, 2
         DO L =  NINDX(LC-1,LREP,4)+1, NINDX(LC,LREP,4)
C
C           Loop over K: irrep(fixed by KLREP,JREP),spinor class,spinors
C
            DO KC =  1, 2
            KMIN = NINDX(KC-1,KREP,3)+1
            IF (KREP.EQ.LREP.AND.TRIANG(2)) KMIN = MAX0(KMIN,L+1)
C
C           Condition for keep/delete of this index pair
C
            IF (KC+LC.EQ.2) THEN
               ALLOW2 = .TRUE.
            ELSEIF (KC+LC.EQ.3.AND.TRIANG(2).AND.
     &             (DODEL(3:4).EQ.'DK'.OR.
     &              DODEL(3:4).EQ.'KD')) THEN
               ALLOW2 = .TRUE.
            ELSE
               ALLOW2 = .FALSE.
            ENDIF
C
            DO 20 K = KMIN, NINDX(KC,KREP,3)
C
C           For the VOVV and VVVV integrals we have only a batch in
C           memory, jump out if we are outside this area
C
            KL = KL + 1
            IF (CLASS.EQ.'VOVV'.OR.CLASS.EQ.'VVVV') THEN
               IF (KLREP.NE.IREPVV) GOTO 40
               IF (KL.LE.ISTART.OR.KL.GT.ISTART+NINT) GOTO 20
            ENDIF
C
C    Loop over J : irrep, spinor class, spinors
C
               DO 10 JREP = 1, NREP
               IREP = MULTB(JREP,KLREP+NREP,2)
               IF (IREP.LT.JREP.AND.TRIANG(1)) GOTO 10
C
               DO JC =  1, 2
               DO J =  NINDX(JC-1,JREP,2)+1, NINDX(JC,JREP,2)
C
C    Loop over I : (irrep fixed via KLREP and JREP), spinor class
C    The loop over spinors is done within XCOPY
C
                  DO IC = 1, 2
                     IMIN = NINDX(IC-1,IREP,1)+1
                     IF (IREP.EQ.JREP.AND.TRIANG(1)) IMIN=MAX0(IMIN,J+1)
                     LNGTH=NINDX(IC,IREP,1)-IMIN+1
C
C                    Condition for keep/delete of this index pair
C
                     IF (IC+JC.EQ.2) THEN
                        ALLOW1 = .TRUE.
                     ELSEIF (IC+JC.EQ.3.AND.TRIANG(1).AND.
     &                      (DODEL(1:2).EQ.'DK'.OR.
     &                       DODEL(1:2).EQ.'KD')) THEN
                        ALLOW1 = .TRUE.
                     ELSE
                        ALLOW1 = .FALSE.
                     ENDIF
C
                     ALLOWED = ALLOW1.AND.ALLOW2
                     IF (LNGTH.GT.0.AND.(.NOT.ALLOWED))
     &                  CALL XCOPY (LNGTH,A0,0,VXXXX(IJKL),1)
                     IF (LNGTH.GT.0) IJKL = IJKL + RCW * LNGTH
                  ENDDO    ! IC
C
               ENDDO    ! J
               ENDDO    ! JC
 10            CONTINUE ! JREP
C
 20         CONTINUE ! K
            ENDDO ! KC
C
         ENDDO    ! L
         ENDDO    ! LC
 30      CONTINUE ! LREP
C
 40   CONTINUE ! KLREP
C
      RETURN
      END


c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      SUBROUTINE DENOM (EPS,T1,T2,CT1,CT2)
C
      use relcc_cfg, only: relcc_ccener_dholu_limit
      implicit none
C
C---------------Description--------------------------------------------
C
C     Divide T1 and T2 by denominators
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 T1(*),T2(*)
      COMPLEX*16 CT1(*),CT2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 FAC,FAC1,FAC2,FAC3,ATVAL,ABS,DHOLU
      integer a,aa,abij,ai,amin,aoff,arp,b,bb,brp,i,ii,ijrp,imin,ioff,j
      integer irp,jj,jrp
C
C---------------Executable code--------------------------------------
C
      DHOLU = relcc_ccener_dholu_limit

      II = 0
      AI = 0
      DO IRP = 1, NREP
      DO I = 1, NO(IRP)
         II = II + 1
         FAC1 = EPS(II)
         AOFF = IO(NREP+1) + IV(IRP)
         DO A = 1, NV(IRP)
             AA = AOFF + A
             FAC = FAC1 - EPS(AA)
             AI = AI + 1
C
C            Check for small denominators. All denominators should be 
C            larger than the minimum HOMO-LUMO gap DHOLU to avoid 
C            convergence problems. An exception is the case in which
C            the numerator is zero (e.g. by symmetry), then we allow
C            execution to be continued. 
C
             IF (-FAC.GT.DHOLU) THEN
                IF (CARITH) THEN
                   CT1(AI) = CT1(AI)/FAC
                ELSE
                   T1(AI) = T1(AI)/FAC
                ENDIF
             ELSE
                IF (CARITH) THEN
                   ATVAL = ABS(CT1(AI))
                ELSE
                   ATVAL = ABS(T1(AI))
                ENDIF
                IF (ATVAL.GT.DHOLU) THEN
                  WRITE(*,*) '1.IRP,NOcc(IRP),NVirt(IRP)=',
     &                        IRP,NO(IRP),NV(IRP)
                  WRITE(*,*) 'II=',II,' AA=',AA,' AI=',AI
                  WRITE(*,*) 'eps_II=',EPS(II),' eps_AA=',EPS(AA)
                  WRITE(*,*) 'EPS(II)-EPS(AA)=FAC=',FAC,
     &          ' ABS(T1(AI)/FAC)=',ATVAL,' DHOLU=',DHOLU
                  CALL QUIT ('Zero or negative HOMO-LUMO gap')
                ELSE
                  IF (CARITH) THEN
                     CT1(AI) = A0
                  ELSE
                     T1(AI) = AR0
                  ENDIF
                ENDIF
             ENDIF
         ENDDO
      ENDDO
      ENDDO
C
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      IOFF = IO(IRP) 
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 - EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 - EPS(AA)
                  ABIJ = ABIJ + 1
C                 Check for small denominators, see comment above.
                  IF (-FAC.GT.DHOLU) THEN
                     IF (CARITH) THEN
                        CT2(ABIJ) = CT2(ABIJ)/FAC
                     ELSE
                        T2(ABIJ) = T2(ABIJ)/FAC
                     ENDIF
                  ELSE
                     IF (CARITH) THEN
                        ATVAL = ABS(CT2(ABIJ))
                     ELSE
                        ATVAL = ABS(T2(ABIJ))
                     ENDIF
                     IF (ATVAL.GT.DHOLU) THEN
                       WRITE(*,*) 'II,JJ,AA,BB=',II,JJ,AA,BB
                       WRITE(*,*) 'eps_II=',EPS(II),'eps_JJ=',EPS(JJ),
     &                             'eps_AA=',EPS(AA),'eps_BB=',EPS(BB)
                       WRITE(*,*) 'eII+eJJ-eAA-eBB=',FAC,' ATVAL=',
     &                 ATVAL,' DHOLU=',DHOLU
                        CALL QUIT ('Zero or negative HOMO-LUMO gap')
                     ELSE
                       IF (CARITH) THEN
                          CT2(ABIJ) = A0
                       ELSE
                          T2(ABIJ) = AR0
                       ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      

      SUBROUTINE DENOM_L (EPS,T1,T2,CT1,CT2)
C
      use relcc_cfg, only: relcc_ccener_dholu_limit
      implicit none
C
C---------------Description--------------------------------------------
C
C     Divide T1 and T2 by denominators
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 T1(*),T2(*)
      COMPLEX*16 CT1(*),CT2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 FAC,FAC1,FAC2,FAC3,ATVAL,DHOLU
      integer a,aa,abij,abrp,amin,aoff,arp,b,bb,brp,i,ia,ii,imin,ioff
      integer irp,iw,j,jj,jrp
C
C---------------Executable code--------------------------------------
C

      DHOLU = relcc_ccener_dholu_limit
      IA = 0
      DO ARP = 1, NREP
         AA = IO(NREP+1) + IV(ARP)
      DO A = 1, NV(ARP)
         AA = AA + 1
         FAC1 = EPS(AA)
         IOFF = IO(ARP)
         DO I = 1, NO(ARP)
             II = IOFF + I
             FAC = FAC1 - EPS(II)
             IA = IA + 1
C
C            Check for small denominators. All denominators should be
C            larger than the minimum HOMO-LUMO gap DHOLU to avoid
C            convergence problems. An exception is the case in which
C            the numerator is zero (e.g. by symmetry), then we allow
C            execution to be continued.
C
             IF (FAC.GT.DHOLU) THEN
                IF (CARITH) THEN
                   CT1(IA) = CT1(IA)/FAC
                ELSE
                   T1(IA) = T1(IA)/FAC
                ENDIF
             ELSE
                IF (CARITH) THEN
                   ATVAL = ABS(CT1(IA))
                ELSE
                   ATVAL = ABS(T1(IA))
                ENDIF
                IF (ATVAL.GT.DHOLU) THEN
                  WRITE(IW,*) 'IRP,NOcc(IRP),NVirt(IRP)=',IRP,NO(ARP),
     &                        NV(ARP)
                  WRITE(IW,*) 'I,eps_I=',I,EPS(II),' A,eps_A=',A,EPS(AA)
                  CALL QUIT ('Zero or negative HOMO-LUMO gap')
                ELSE
                  IF (CARITH) THEN
                     CT1(IA) = A0
                  ELSE
                     T1(IA) = AR0
                  ENDIF
                ENDIF
             ENDIF
         ENDDO
      ENDDO
      ENDDO
C
      ABIJ = 0
      DO ABRP = 1, NREP
      DO 10 BRP = 1, NREP
      BB = IV(BRP) + IO(NREP+1)
      ARP = MULTB(BRP,ABRP+NREP,2)
      IF (ARP.LT.BRP) GOTO 10
      DO B = 1, NV(BRP)
         BB = BB + 1
         FAC1 = EPS(BB)
            AOFF = IV(ARP) + IO(NREP+1)
           AMIN = 1
           IF (ARP.EQ.BRP) AMIN = B + 1
         DO A = AMIN, NV(ARP)
            AA = AOFF + A
            FAC2 = EPS(AA) + FAC1
            DO 20 JRP = 1, NREP
            JJ = IO(JRP)
            IRP = MULTB(JRP,ABRP+NREP,2)
            IF (IRP.LT.JRP) GOTO 20
            IOFF = IO(IRP)
            DO J = 1, NO(JRP)
               JJ = JJ + 1
               FAC3 = FAC2 - EPS(JJ)
               IMIN = 1
               IF (IRP.EQ.JRP) IMIN = J + 1
               DO I = IMIN, NO(IRP)
                  II = IOFF + I
                  FAC = FAC3 - EPS(II)
                  ABIJ = ABIJ + 1
C                 Check for small denominators, see comment above.
                  IF (FAC.GT.DHOLU) THEN
                     IF (CARITH) THEN
                        CT2(ABIJ) = CT2(ABIJ)/FAC
                     ELSE
                        T2(ABIJ) = T2(ABIJ)/FAC
                     ENDIF
                  ELSE
                     IF (CARITH) THEN
                        ATVAL = ABS(CT2(ABIJ))
                     ELSE
                        ATVAL = ABS(T2(ABIJ))
                     ENDIF
                     IF (ATVAL.GT.DHOLU) THEN
             WRITE(IW,*) 'IRP,NOcc(IRP),NVirt(IRP)=',IRP,NO(IRP),NV(IRP)
             WRITE(IW,*) 'I,eps_I=',I,EPS(II),' A,eps_A=',A,EPS(AA)
                        CALL QUIT ('Zero or negative HOMO-LUMO gap')
                     ELSE
                       IF (CARITH) THEN
                          CT2(ABIJ) = A0
                       ELSE
                          T2(ABIJ) = AR0
                       ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE DENOM_N (EPS,DINV_ONE,DINV_TWO)
C
      use relcc_cfg, only: relcc_ccener_dholu_limit

      implicit none
C
C---------------Description--------------------------------------------
C
C     Collect the inverse of orbital energy differences in an array. 
C
C---------------Last modified------------------------------------------
C
C     Author : Avijit Shee
C
C---------------Calling variables--------------------------------------
C

      real*8,intent(in)    :: EPS(*)
      real*8,intent(out) :: DINV_ONE(*),DINV_TWO(*)

C---------------Common Blocks--------------------------------------

#include "symm.inc"
#include "complex.inc"
#include "param.inc"

C---------------Local variables--------------------------------------

      REAL*8 FAC,FAC1,FAC2,FAC3,ATVAL,DHOLU,CONST
      integer a,aa,abij,abrp,amin,aoff,arp,b,bb,brp,i,ia,ii,imin,j,jj
      integer ioff,irp,jrp

C---------------Executable code--------------------------------------


      DHOLU = relcc_ccener_dholu_limit
      IA = 1
      DO ARP = 1, NREP
         AA = IO(NREP+1) + IV(ARP)
      DO A = 1, NV(ARP)
         AA = AA + 1
         FAC1 = EPS(AA)
         IOFF = IO(ARP)
         DO I = 1, NO(ARP)
             II = IOFF + I
             FAC = FAC1 - EPS(II)

             IF (FAC.GT.DHOLU) THEN
                   CONST = AR1/FAC
                   DINV_ONE(IA:(IA+(RCW-1))) = CONST
             ELSE
                   DINV_ONE(IA:(IA+(RCW-1))) = AR0
             ENDIF

             IA = IA + RCW
         ENDDO
      ENDDO
      ENDDO

      ABIJ = 1
      DO ABRP = 1, NREP
      DO 10 BRP = 1, NREP
      BB = IV(BRP) + IO(NREP+1)
      ARP = MULTB(BRP,ABRP+NREP,2)
      IF (ARP.LT.BRP) GOTO 10
      DO B = 1, NV(BRP)
         BB = BB + 1
         FAC1 = EPS(BB)
            AOFF = IV(ARP) + IO(NREP+1)
           AMIN = 1
           IF (ARP.EQ.BRP) AMIN = B + 1
         DO A = AMIN, NV(ARP)
            AA = AOFF + A
            FAC2 = EPS(AA) + FAC1
            DO 20 JRP = 1, NREP
            JJ = IO(JRP)
            IRP = MULTB(JRP,ABRP+NREP,2)
            IF (IRP.LT.JRP) GOTO 20
            IOFF = IO(IRP)
            DO J = 1, NO(JRP)
               JJ = JJ + 1
               FAC3 = FAC2 - EPS(JJ)
               IMIN = 1
               IF (IRP.EQ.JRP) IMIN = J + 1
               DO I = IMIN, NO(IRP)
                  II = IOFF + I
                  FAC = FAC3 - EPS(II)
                  IF (FAC.GT.DHOLU) THEN
                        CONST = AR1/FAC
                        DINV_TWO(ABIJ:(ABIJ+(RCW-1))) = CONST
                  ELSE
                    DINV_TWO(ABIJ:(ABIJ+(RCW-1))) = AR0
                  ENDIF
                  ABIJ = ABIJ + RCW
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      RETURN
      END

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      

      SUBROUTINE PAIREN (T1,T2,FVO,BUF1,BUF2)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Attribute correlation energy to electron pairs
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      REAL*8 T1(*),T2(*)
      INTEGER FVO
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "param.inc"
#include "files.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 EI,EIJ
      integer a,abij,ai,amin,arp,b,brp,i,ijrp,imin,irp,j,jrp
C
C---------------Executable code--------------------------------------
C
      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface


      WRITE (IW,1000)
C
      CALL XCOPY (NDIMT1,FVO,1,BUF1,1)
      CALL DELFCK ('VO','DD',BUF1)
C
      WRITE (IW,2000)
      AI = 1
      DO IRP = 1, NREP
      DO I = 1, NO(IRP)
         EI = AR0
         DO A = 1, NV(IRP)
             EI = EI + BUF1(AI) * T1(AI)
             IF (CARITH) THEN
                EI = EI + BUF1(AI+1) * T1(AI+1)
             ENDIF
             AI = AI + RCW
         ENDDO
         IF (ABS(EI).GT.1.D-6) WRITE (IW,2001) IRP,I,EI
      ENDDO
      ENDDO
C
      CALL GETVVOO (BUF1)
      CALL DELINT ('VVOO','DDDD',BUF1,0,0,0)
      CALL GETTAU  (T1,T2,BUF2)
C
      WRITE (IW,3000)
      ABIJ = 1
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            EIJ = AR0
            DO 20 BRP = 1, NREP
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            DO B = 1, NV(BRP)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  EIJ = EIJ + BUF1(ABIJ) * BUF2(ABIJ)
                  IF (CARITH) THEN
                     EIJ = EIJ + BUF1(ABIJ+1) * BUF2(ABIJ+1)
                  ENDIF
                  ABIJ = ABIJ + RCW
               ENDDO
            ENDDO
 20         CONTINUE
            IF (ABS(EIJ).GT.1.D-6) WRITE (IW,3001) IRP,I,JRP,J,EIJ
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      WRITE (IW,4000)
C
 1000 FORMAT (//' Individual terms in the correlation energy')
 2000 FORMAT (/'   Irrep Spinor',5X,'Relaxation energy (T1-amplitudes)')
 2001 FORMAT (I7,I8,F30.15)
 3000 FORMAT (/'   Irrep1 Spinor1 Irrep2 Spinor2',5X,
     &'Correlation energy (T2-amplitudes)')
 3001 FORMAT (2(I7,I8),F30.15)
 4000 FORMAT (//)
C
      RETURN
      END
      SUBROUTINE PARTS(IMODE,RMS,T1,T2,BUF1,BUF2)
C
      use interface_to_mpi
      implicit none
C
C---------------Description--------------------------------------------
C
C     IMODE :
C     1) Sums partial contributions to T1 and T2
C     2) Synchronizes T1 and T2
C
C     BUF1 and BUF2 have the length of T1 and T2 respectively and are
C     not used with imode 2. RMS is also updated because it controls
C     the iterative process
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*),BUF1(*),BUF2(*)
      INTEGER IMODE,RMS
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
      integer ISTAT(df_MPI_STATUS_SIZE)
#endif
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
#if defined (VAR_MPI)
      IF (NMPROC.EQ.1) RETURN
      IF (IMODE.EQ.1) THEN
          call interface_mpi_reduce_r1_work_f77(T1,BUF1,NDIMT1*RCW,
     &                    op_MPI_SUM,MASTER,global_communicator)
          call interface_mpi_reduce_r1_work_f77(T2,BUF2,NDIMT2*RCW,
     &                    op_MPI_SUM,MASTER,global_communicator)
          IF (MYPROC.EQ.MASTER) CALL XCOPY (NDIMT1,BUF1,1,T1,1)
          IF (MYPROC.EQ.MASTER) CALL XCOPY (NDIMT2,BUF2,1,T2,1)
      ENDIF
C
      call interface_mpi_bcast_r1_work_f77(T1,NDIMT1*RCW,MASTER,
     &               global_communicator)
      call interface_mpi_bcast_r1_work_f77(T2,NDIMT2*RCW,MASTER,
     &               global_communicator)
      call interface_mpi_bcast_r1_work_f77(RMS,1,MASTER,
     &               global_communicator)
#endif
      RETURN
      END
      SUBROUTINE S1S2_ADDFO (FVO,BUF1,BUF2,S1,S2)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Add the first order wavefunction contribution to S1 and S2.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2,FVO(*)
      REAL*8 S1(*),S2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C---------------------------------------
C Only the master gets this contribution
C---------------------------------------
      IF (MYPROC.NE.MASTER) RETURN
C-------------------------------
C S(A,I) = S(A,I) + F(A,I)
C-------------------------------
      IF (EQNS.EQ.'FOCKSP') THEN
C       Delete first order Hamiltonian parts
C       first for Fockspace
        CALL XCOPY (NDIMT1,FVO,1,BUF1,1)
        CALL DENOMF (BUF1,BUF1,BUF2,BUF1,BUF2,2)
        CALL XAXPY (NDIMT1,A1,BUF1,1,S1,1)
      ELSE
        CALL XAXPY (NDIMT1,A1,FVO,1,S1,1)
      ENDIF
C-------------------------------
C S(AB,IJ) = S(AB,IJ) + W(AB,IJ)
C-------------------------------
      CALL GETVVOO (BUF2)
      CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)
C
      RETURN
      END
      SUBROUTINE ZCORE(T1,T2)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Zeroes out core contributions (frozen spinors)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "freeze.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      integer i,ijrp,imin,irp,j,jrp,n,off
C
C---------------Executable code--------------------------------------
C
      OFF = 1
      DO IRP = 1, NREP
         DO I = 1, NO(IRP)
            IF (IFROZ(I,IRP).NE.0)
     &      CALL XCOPY (NV(IRP),A0,0,T1(OFF),1)
            OFF = OFF + NV(IRP) * RCW
         ENDDO
      ENDDO
C
      OFF = 1
      DO IJRP = 1, NREP
      N = NVVT(IJRP)
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            IF (IFROZ(I,IRP).NE.0.OR.IFROZ(J,JRP).NE.0)
     &         CALL XCOPY (N,A0,0,T2(OFF),1)
            OFF = OFF + N * RCW
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
C
      RETURN
      END

      SUBROUTINE FCORE(T1,T2,EPS)
      use relcc_cfg  
      use projectors 
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Zeroes out core contributions (frozen core approximation)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*),EPS(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      integer :: idet, j, k
      logical, pointer :: projector(:) => NULL()
C
C---------------Executable code--------------------------------------
C
      if (.not.relcc_projectors_frozen_core) RETURN

      projector => create_t_projector(NDIMT1+NDIMT2,EPS)
       
      do idet=1, (NDIMT1+NDIMT2)
        if (projector(idet)) then
          if (idet.le.NDIMT1) then
            if(CARITH) then
              j = 2*(idet-1) + 1
              k = j + 1
              T1(j) = 0d0
              T1(k) = 0d0
            else
              T1(idet) = 0d0
            endif
          else 
            if(CARITH) then
              j = 2*(idet-1) + 1
              k = j + 1
              T2(j-NDIMT1*2) = 0d0
              T2(k-NDIMT1*2) = 0d0
            else
              T2(idet-NDIMT1) = 0d0
            endif
          endif
        endif
      end do
C
      RETURN
      END


      SUBROUTINE MP2ENG(T1,T2,FVO,BUF1,EMP2)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculates MP2 energy. Checks result (should be real).
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
      REAL*8 T1(*),T2(*)
      REAL*8 FVO(*)
      REAL*8 BUF1(*)
      REAL*8 EMP2
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "files.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CEMP21(2),CEMP22(2)
C
C---------------Executable code--------------------------------------
C
      CALL XDOTC(CEMP21,NFVO,FVO,1,T1,1)
      EMP2 = CEMP21(1)
      IF (CARITH) THEN
         IF (ABS(CEMP21(2)).GT.ACCUR)
     &      WRITE (IW,1000) 'T1',CEMP21(2)
      ENDIF
C
c     print*,'T1 contribution to MP2 energy',emp2
C
      CALL GETVVOO(BUF1)
      CALL XDOTC(CEMP22,NDIMT2,BUF1,1,T2,1)
      EMP2 = EMP2 + CEMP22(1)
      IF (CARITH) THEN
         IF (ABS(CEMP22(2)).GT.ACCUR)
     &      WRITE (IW,1000) 'T2',CEMP22(2)
      ENDIF
C
      RETURN
1000  FORMAT (/' WARNING : Imaginary contribution to ',A2,
     &' part of MP2 energy :',F20.15)
      END
      SUBROUTINE T1DIAG(DIAG,T1,BUF1)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Calculate T1-diagnostic
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),BUF1(*),DIAG
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER SUM,I
      REAL*8 RSUM,DDOT,SDOT
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (NDIMT1,T1,1,BUF1,1)
      CALL DELFCK ('VO','DD',BUF1)
C
      SUM = 0
      DO I = 1, NREP
         SUM = SUM + NE(I)
      ENDDO
      RSUM = DBLE(SUM)
C
      DIAG = DDOT (RCW*NDIMT1,BUF1,1,BUF1,1)
      DIAG =  SQRT(DIAG/RSUM)
C
      RETURN
      END
