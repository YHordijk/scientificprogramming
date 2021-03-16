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
C
#ifdef UNDEF
/* Comdeck log */
961004 - tsaue INTFCK has undergone major surgery:
               - Inner loop now over integrals, prevents cache misses
               - Redefinition of IFCTYP
               - Possibility of screening on Coulomb/exchange
               - Indices unpacked in buffers
               - Deletetion of all "zero" integrals
941019-hjaaj
changed ISYMDM, IFCTYP to ISYMDM(NDMAT), IFCTYP(NDMAT)
implemented all six IFCTYPes (two was already implemented)
improved efficiency of INTFCK by moving GINT outside and
 moving common factors to SKLFCK_1
halved no. of  operations for IFCTYP=2 in INTFCK
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck intfck */
      SUBROUTINE INTFCK(FMAT,AOINT,DMAT,NDMAT,NCCINT,NINTYP,
     &                  ICORBA,ICORBB,ICORBC,ICORBD,SYMFAC,
     &                  IPRINT,NODV,NINDAB,NINDCD,SUSCEP,IFCTYP,
     &                  DNSBUF,DINTSKP,HFXFAC,WORK,LWORK)
C*****************************************************************************
C
C  Direct calculation of Fock matrices in AO-basis
C
C  Major surgery Aug 26 1996 by T.Saue
C
C  Written by Trygve Helgaker and Henrik Koch  93-93 ??
C
C  All IFCTYPs added and partly optimized Nov. 1994 by
C  Hans Joergen Aa. Jensen
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER(D0 = 0.0D0)
      LOGICAL NODV, SUSCEP
      DIMENSION ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBD),
     &          AOINT(NCCINT,NINTYP), FMAT(NBASIS,NBASIS,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          DMAT(NBASIS,NBASIS,NDMAT), IFCTYP(NDMAT),
     &          DNSBUF(2,NDMAT),DINTSKP(*),WORK(LWORK)
#include "twocom.h"
#include "nuclei.h"
C
      CALL QENTER('INTFCK')
#include "memint.h"
C
      IF(DOSCRN) THEN
        FCKTOL = MAX((SCRTHR/DNSMAX),1.00D-15)
      ELSE
        FCKTOL = 1.00D-15
      ENDIF
      CALL MEMGET2('INTE','IND',KIND ,4*NCCINT,WORK,KFREE,LFREE)
C
      CALL INTFC1(FMAT,AOINT,DMAT,NDMAT,NCCINT,NINTYP,
     &                  ICORBA,ICORBB,ICORBC,ICORBD,FCKTOL,SYMFAC,
     &                  IPRINT,NINDAB,NINDCD,SUSCEP,IFCTYP,
     &                  DNSBUF,DINTSKP,HFXFAC,WORK(KIND))
C
      CALL MEMREL('INTFCK',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('INTFCK')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck intfc1 */
      SUBROUTINE INTFC1(FMAT,AOINT,DMAT,NDMAT,NCCINT,NINTYP,
     &                  ICORBA,ICORBB,ICORBC,ICORBD,FCKTOL,SYMFAC,
     &                  IPRINT,NINDAB,NINDCD,SUSCEP,IFCTYP,
     &                  DNSBUF,DINTSKP,HFXFAC,IND)
C*****************************************************************************
C
C  Direct calculation of Fock matrices in AO-basis
C
C  Revised Aug 26 1996 by T.Saue
C
C  Written by Trygve Helgaker and Henrik Koch  93-93 ??
C
C  All IFCTYPs added and partly optimized Nov. 1994 by
C  Hans Joergen Aa. Jensen
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER(DMP25 = -0.25D00,DMP5=-0.5D0,D1=1.0D0,D2=2.0D0,D0=0.0D0,
     &          DM1=-1.D0, DM3=-3.D0)
      INTEGER A,B,C,D
      LOGICAL SUSCEP, NOTEST,LBIT,DOC,DOE
      DIMENSION ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBD),
     &          AOINT(NCCINT), FMAT(NBASIS,NBASIS,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          DMAT(NBASIS,NBASIS,NDMAT), IFCTYP(NDMAT),
     &          DNSBUF(2,NDMAT),IND(4,*),DINTSKP(2,4)
#include "twocom.h"
#include "nuclei.h"
C
      IF (IPRINT .GT. 9) THEN
         CALL HEADER('Output from INTFC1',-1)
         WRITE (LUPRI, '(A,2L5)') ' SUSCEP', SUSCEP
         WRITE (LUPRI, '(A,4I5)') ' NORB ', NORBA,NORBB,NORBC,NORBD
         WRITE (LUPRI, '(A,2L5)') ' DIAGAB/CD ', DIAGAB,DIAGCD
         WRITE (LUPRI, '(A,3L5)') ' SHAEQB, SHCEQD, SHABAB ',
     *                              SHAEQB, SHCEQD, SHABAB
         WRITE (LUPRI, '(A,4I5)') ' NHKTA', NHKTA,NHKTB,NHKTC,NHKTD
         WRITE (LUPRI, '(A,4I5)') ' KHKTA', KHKTA,KHKTB,KHKTC,KHKTD
         WRITE (LUPRI, '(A,F12.6)') ' FCKTOL ', FCKTOL
         WRITE (LUPRI, '(A,F12.6)') ' SYMFAC ', SYMFAC
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals A '
         WRITE (LUPRI, '(20I5)') (ICORBA(I),I=1, NORBA)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals B '
         WRITE (LUPRI, '(20I5)') (ICORBB(I),I=1, NORBB)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals C '
         WRITE (LUPRI, '(20I5)') (ICORBC(I),I=1, NORBC)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals D '
         WRITE (LUPRI, '(20I5)') (ICORBD(I),I=1, NORBD)
      END IF
C
C     ****************************************************
C     ***** Generate indices and permutation factors *****
C     ****************************************************
C
      DOC = DOSCRN.AND.LBIT(ICEFLG,1)
      DOE = DOSCRN.AND.LBIT(ICEFLG,2)
      NOTEST = .NOT.(DIAGAB.OR.DIAGCD.OR.TCONAB.OR.TCONCD)
      SFAC = SYMFAC
      IF (.NOT.SHABAB) SFAC = SFAC+SFAC
      IF (.NOT.SHAEQB) SFAC = SFAC+SFAC
      IF (.NOT.SHCEQD) SFAC = SFAC+SFAC
      TIND = SECOND()
      NBUF  = 0
      IOFF  = 0
      AOMAX = D0
        IF(NOTEST) THEN
          DO ICOMPA = 1,KHKTA
            DO ICOMPB = 1,KHKTB
              DO ICOMPC = 1,KHKTC
                DO ICOMPD = 1,KHKTD
                  DO IORBAB = 1,NORBAB
                    IF (SUSCEP) THEN
                      IORBA = 1 + NINDAB(IORBAB,1)/KHKTA
                      IORBB = 1 + NINDAB(IORBAB,2)/KHKTB
                    ELSE
                      IORBA = NINDAB(IORBAB,1)
                      IORBB = NINDAB(IORBAB,2)
                    END IF
                    A = ICORBA(IORBA) + ICOMPA
                    B = ICORBB(IORBB) + ICOMPB
                    DO IORBCD = 1, NORBCD
                      AINT  = AOINT(IOFF+IORBCD)
                      IF(ABS(AINT).GT.FCKTOL) THEN
                        NBUF  = NBUF + 1
                        IF (SUSCEP) THEN
                          IORBC = 1 + NINDCD(IORBCD,1)/KHKTC
                          IORBD = 1 + NINDCD(IORBCD,2)/KHKTD
                        ELSE
                          IORBC = NINDCD(IORBCD,1)
                          IORBD = NINDCD(IORBCD,2)
                        END IF
                        AOMAX       = MAX(AOMAX,ABS(AINT))
                        AOINT(NBUF) = AINT
                        IND(1,NBUF) = A
                        IND(2,NBUF) = B
                        IND(3,NBUF) = ICORBC(IORBC) + ICOMPC
                        IND(4,NBUF) = ICORBD(IORBD) + ICOMPD
                      ENDIF
                    ENDDO
                    IOFF = IOFF + NORBCD
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF(SFAC.NE.D1) CALL DSCAL(NBUF,SFAC,AOINT,1)
        ELSE
          DO ICOMPA = 1,KHKTA
                        KHKTBB = KHKTB
            IF (DIAGAB) KHKTBB = ICOMPA
            DO ICOMPB = 1,KHKTBB
              FACAB = SFAC
              IF (DIAGAB.AND.ICOMPA.NE.ICOMPB)
     &          FACAB = FACAB+FACAB
              DO ICOMPC = 1,KHKTC
                            KHKTDD = KHKTD
                IF (DIAGCD) KHKTDD = ICOMPC
                DO ICOMPD = 1,KHKTDD
                  FCABCD = FACAB
                  IF (DIAGCD.AND.ICOMPC.NE.ICOMPD)
     &              FCABCD = FCABCD+FCABCD
                  DO IORBAB = 1,NORBAB
                    IF (SUSCEP) THEN
                      IORBA = 1 + NINDAB(IORBAB,1)/KHKTA
                      IORBB = 1 + NINDAB(IORBAB,2)/KHKTB
                    ELSE
                      IORBA = NINDAB(IORBAB,1)
                      IORBB = NINDAB(IORBAB,2)
                    END IF
                    FAB = FCABCD
                    IF (TCONAB.AND.IORBA.NE.IORBB)
     &                FAB = FAB+FAB
                    A = ICORBA(IORBA) + ICOMPA
                    B = ICORBB(IORBB) + ICOMPB
                    DO IORBCD = 1, NORBCD
                      AINT  = AOINT(IOFF+IORBCD)
                      IF(ABS(AINT).GT.FCKTOL) THEN
                        NBUF  = NBUF + 1
                        IF (SUSCEP) THEN
                          IORBC = 1 + NINDCD(IORBCD,1)/KHKTC
                          IORBD = 1 + NINDCD(IORBCD,2)/KHKTD
                        ELSE
                          IORBC = NINDCD(IORBCD,1)
                          IORBD = NINDCD(IORBCD,2)
                        END IF
                        FABCD = FAB
                        IF (TCONCD.AND.IORBC.NE.IORBD)
     &                     FABCD = FABCD+FABCD
                        AOMAX       = MAX(AOMAX,ABS(AINT))
                        AOINT(NBUF) = FABCD*AINT
                        IND(1,NBUF) = A
                        IND(2,NBUF) = B
                        IND(3,NBUF) = ICORBC(IORBC) + ICOMPC
                        IND(4,NBUF) = ICORBD(IORBD) + ICOMPD
                      ENDIF
                    ENDDO
                    IOFF = IOFF + NORBCD
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
C     Statistics
C
      NDEL = NCCINT - NBUF
      DINTSKP(1,2) = DINTSKP(1,2) + NCCINT
      DINTSKP(2,2) = DINTSKP(2,2) + NDEL
      DINTSKP(1,3) = DINTSKP(1,3) + NCM*NBUF
      DINTSKP(1,4) = DINTSKP(1,4) + NEM*NBUF
      TIND = SECOND() - TIND
C
C     **********************************
C     ***** Print section  *************
C     **********************************
C
      IF(IPRINT.GE.15) THEN
        DO INT = 1,NBUF
          WRITE (LUPRI,'(1X,A,F16.12,4I5)')
     &      'GINT,A,B,C,D',AOINT(INT),
     &      IND(1,INT),IND(2,INT),IND(3,INT),IND(4,INT)
        ENDDO
      ENDIF
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 General: no symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C           X = 3 Mixed: some blocks symmetric and some antisymmetric
C         G indicates whether we build Coulomb or Gaunt
C           G = 0 Coulomb
C           G = 1 Gaunt
C
C         Y indicates contributions (Direct and Exchange)
C
C         For Coulomb they have the following meaning:
C           Y = 0 Nothing (neither direct nor exchange)
C           Y = 1 Direct
C           Y = 2 Exchange
C           Y = 3 Direct + Exchange
C
C         For Gaunt we need to distinguish between real and imag. matrices
C         See the note on the implementation of Gaunt:
C         There is never a direct contribution for IZ=1.
C         The exchange part for IZ=1 scales by 3, not -1.
C
C           Y = 4 Nothing           for IZ = 1
C           Y = 5 Exchange          for IZ = 1
C           Y = 6 Nothing           for IZ > 1
C           Y = 7 Direct            for IZ > 1
C           Y = 8 Exchange          for IZ > 1
C           Y = 9 Direct + Exchange for IZ > 1
C
      TFCK = SECOND()
      DO 400 I = 1,NDMAT
        FCM = DNSBUF(1,I)*AOMAX
        FEM = DNSBUF(2,I)*AOMAX
C       Unpack the bitcode
        IX = IFCTYP(I)/10
        IY = MOD(IFCTYP(I),10)
        IG = IY/4
C       We can build Coulomb and Gaunt type matrices by specifying the
C       appropriate multiplication factors.
        IF (IG.EQ.0) THEN
C..........Coulomb interaction
           IC = MOD(IY,2)
           IE = IY - IC
C          ... no direct term for antisymmetric density matrix
           IF (IX.EQ.2) IC = 0
C
C          Exchange may be switched off in DFT calculations
C
           IF (HFXFAC.EQ.D0) IE = 0
        ELSE
C..........Gaunt interaction
           IF (IY.LE.5) THEN
              IZ = 1
              IC = 0
              IE = 2*MOD(IY-4,2)
           ELSE
              IZ = 2 ! It could also be 3 or 4, makes no difference
              IY = IY-6
              IC = MOD(IY,2)
              IE = IY - IC
           ENDIF
C          ... no direct term for symmetric density matrix
           IF (IX.EQ.1) IC = 0
C
C          Exchange may be switched off in DFT calculations
C
           IF (LSCALE_GAUNT .AND. HFXFAC.EQ.D0) IE = 0
        ENDIF
C       Screening on direct contributions
        IF(DOC.AND.IC.EQ.1.AND.FCM.LT.SCRTHR) THEN
          IC = 0
          DINTSKP(2,3) = DINTSKP(2,3) + NBUF
        ENDIF
C       Screening on exchange contributions
        IF(DOE.AND.IE.EQ.2.AND.FEM.LT.SCRTHR) THEN
          IE = 0
          DINTSKP(2,4) = DINTSKP(2,4) + NBUF
        ENDIF
        IY = IC + IE
        IF(IY.EQ.0) GOTO 400
        IF(IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
          IF (IG.EQ.0) THEN
C............Coulomb multiplication factors
             DIRFAC = D1
             EXCFAC = DMP25*HFXFAC
          ELSE
C............Gaunt multiplication factors
             IF (IZ.EQ.1) THEN
                DIRFAC = D0
                IF (LSCALE_GAUNT) THEN
                   EXCFAC = DM3*DMP25*HFXFAC
                ELSE
                   EXCFAC = DM3*DMP25
                ENDIF
             ELSE
                DIRFAC = DM1
                IF (LSCALE_GAUNT) THEN
                   EXCFAC = DMP25*HFXFAC
                ELSE
                   EXCFAC = DMP25
                ENDIF
             ENDIF
          ENDIF
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (1/4) * (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
            DO INT = 1,NBUF
              A = IND(1,INT)
              B = IND(2,INT)
              C = IND(3,INT)
              D = IND(4,INT)
              DINT = DIRFAC*AOINT(INT)
              EINT = EXCFAC*AOINT(INT)
C...............Direct contributions
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(C,D,I)
              FMAT(C,D,I) = FMAT(C,D,I) + DINT*DMAT(A,B,I)
C...............Exchange (LL or SS for Gaunt)
              FMAT(C,A,I) = FMAT(C,A,I) + EINT*DMAT(D,B,I)
              FMAT(D,B,I) = FMAT(D,B,I) + EINT*DMAT(C,A,I)
C...............Exchange (LS or SL for Gaunt)
              FMAT(D,A,I) = FMAT(D,A,I) + EINT*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) + EINT*DMAT(D,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (1/4) (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
            DO INT = 1,NBUF
              A = IND(1,INT)
              B = IND(2,INT)
              C = IND(3,INT)
              D = IND(4,INT)
              EINT = EXCFAC*AOINT(INT)
C...............Exchange (LL or SS for Gaunt)
              FMAT(C,A,I) = FMAT(C,A,I) + EINT*DMAT(D,B,I)
              FMAT(D,B,I) = FMAT(D,B,I) + EINT*DMAT(C,A,I)
C...............Exchange (LS or SL for Gaunt)
              FMAT(D,A,I) = FMAT(D,A,I) + EINT*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) + EINT*DMAT(D,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
            DO INT = 1,NBUF
              A = IND(1,INT)
              B = IND(2,INT)
              C = IND(3,INT)
              D = IND(4,INT)
              DINT = DIRFAC*AOINT(INT)
C...............Direct contributions
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(C,D,I)
              FMAT(C,D,I) = FMAT(C,D,I) + DINT*DMAT(A,B,I)
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'INTFC1 ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in INTFC1: specified IFCTYP not implemented.')
          ENDIF
        ELSE IF (IX.EQ.0) THEN
          DIRFAC = D1
          EXCFAC = DMP5*HFXFAC
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          IF(IY.EQ.3) THEN
            DO INT = 1,NBUF
              A = IND(1,INT)
              B = IND(2,INT)
              C = IND(3,INT)
              D = IND(4,INT)
              DINT = DIRFAC*AOINT(INT)
C...............Direct contributions
              GCD  = DINT*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GCD
              FMAT(B,A,I) = FMAT(B,A,I) + GCD
              GAB  = DINT*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I) = FMAT(C,D,I) + GAB
              FMAT(D,C,I) = FMAT(D,C,I) + GAB
              EINT = EXCFAC*AOINT(INT)
C...............Exchange (LL or SS for Gaunt)
              FMAT(C,A,I) = FMAT(C,A,I) + EINT*DMAT(D,B,I)
              FMAT(A,C,I) = FMAT(A,C,I) + EINT*DMAT(B,D,I)
              FMAT(D,B,I) = FMAT(D,B,I) + EINT*DMAT(C,A,I)
              FMAT(B,D,I) = FMAT(B,D,I) + EINT*DMAT(A,C,I)
C...............Exchange (LS or SL for Gaunt)
              FMAT(D,A,I) = FMAT(D,A,I) + EINT*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) + EINT*DMAT(D,A,I)
              FMAT(A,D,I) = FMAT(A,D,I) + EINT*DMAT(B,C,I)
              FMAT(B,C,I) = FMAT(B,C,I) + EINT*DMAT(A,D,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
            DO INT = 1,NBUF
              A = IND(1,INT)
              B = IND(2,INT)
              C = IND(3,INT)
              D = IND(4,INT)
              EINT = EXCFAC*AOINT(INT)
C...............Exchange (LL or SS for Gaunt)
              FMAT(C,A,I) = FMAT(C,A,I) + EINT*DMAT(D,B,I)
              FMAT(A,C,I) = FMAT(A,C,I) + EINT*DMAT(B,D,I)
              FMAT(D,B,I) = FMAT(D,B,I) + EINT*DMAT(C,A,I)
              FMAT(B,D,I) = FMAT(B,D,I) + EINT*DMAT(A,C,I)
C...............Exchange (LS or SL for Gaunt)
              FMAT(D,A,I) = FMAT(D,A,I) + EINT*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) + EINT*DMAT(D,A,I)
              FMAT(A,D,I) = FMAT(A,D,I) + EINT*DMAT(B,C,I)
              FMAT(B,C,I) = FMAT(B,C,I) + EINT*DMAT(A,D,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
            DO INT = 1,NBUF
              A = IND(1,INT)
              B = IND(2,INT)
              C = IND(3,INT)
              D = IND(4,INT)
              DINT = DIRFAC*AOINT(INT)
C...............Direct contributions
              GCD  = DINT*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GCD
              FMAT(B,A,I) = FMAT(B,A,I) + GCD
              GAB  = DINT*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I) = FMAT(C,D,I) + GAB
              FMAT(D,C,I) = FMAT(D,C,I) + GAB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'INTFC1 ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in INTFC1: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'INTFC1 ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in INTFC1: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
      TFCK = SECOND() - TFCK
C
      IF(IPRINT.GE.4) THEN
        WRITE(LUPRI,'(A)')       '* INTFC1 - statistics:'
        WRITE(LUPRI,'(A,I10)')   '*   Integrals in     :',NCCINT
        WRITE(LUPRI,'(A,I10)')   '*   Integrals out    :',NBUF
        WRITE(LUPRI,'(A,F15.6)') '*   Time indices     :',TIND
        WRITE(LUPRI,'(A,F15.6)') '*   Time Fock matrix :',TFCK
        CALL FLSHFO(LUPRI)
      ENDIF
      RETURN
C
      END
C  /* Deck dsotao */
      SUBROUTINE DSOTAO(DSO,DAO,NBAST,IREPDM,IPRINT)
C
C     Take density matrix in symmetry orbital basis and generate
C     density matrix over distinct pairs of AOs
C
C                                          880418  PRT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION DSO(NBAST,NBAST), DAO(NBAST,NBAST)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
      IF (IPRINT .GT. 10) CALL HEADER('Subroutine DSOTAO',-1)
C
C     Loop over all irreps in molecule
C
      ISTRA = 1
      CALL DZERO(DAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         NORBA = NAOS(IREPA+1)
         DO 200 I = ISTRA,ISTRA + NORBA - 1
            IA   = IAND(ISHFT(IPIND(I),-16),65535)
            NA   = IAND(ISHFT(IPIND(I), -8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            DO 300 ISYMA = 0, MAXOPR
            IF (IAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = PT(IAND(ISYMA,IEOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
                  NORBB = NAOS(IREPB+1)
                  IF (IEOR(IREPA,IREPB).EQ.IREPDM) THEN
                  DO 500 J = ISTRB,ISTRB + NORBB - 1
                     IB   = IAND(ISHFT(IPIND(J),-16),65535)
                     NB   = IAND(ISHFT(IPIND(J), -8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IAND(ISYMB,
     &                             IEOR(IREPB,ISYMAO(NHKTB,NB))))
                        DAO(INDA,INDB) = DAO(INDA,INDB)
     &                                 + FACA*FACB*DSO(I,J)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NORBB
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NORBA
  100 CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Total density matrix in SO basis',-1)
         CALL OUTPUT(DSO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         CALL HEADER('Total density matrix in AO basis',-1)
         CALL OUTPUT(DAO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck sklfck */
      SUBROUTINE SKLFCK(FMAT,WORK,LWORK,IPRINT,DIRFCK,DDFOCK,EXPECT,
     &                  PERTUR,NODV,MAXDER,LONDON,NDMAT,ISYMDM,IFCTYP,
     &                  IATOME)
C*****************************************************************************
C
C      Revised Oct 11 1996 by T.Saue
C      - changed IFCTYP !
C      - symmetrization/anti-symmetrization should be done aftewards !
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "iratdef.h"
      LOGICAL   DIRFCK, DDFOCK, NODV, EXPECT, LONDON, PERTUR
      REAL*8    FMAT(NBASIS*NBASIS,*), WORK(LWORK) 
      INTEGER   ISYMDM(NDMAT), IFCTYP(NDMAT)
#include "abainf.h"
#include "ccom.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "nuctap.h"
CTROND: Take out inforb and abainf if possible !
      CALL QENTER('SKLFCK')
#include "memint.h"
      N2BASX = NBASIS*NBASIS
      NNBASX = (NBASIS*(NBASIS+1))/2
C
C     Undifferentiated Fock matrix
C     ============================
C
      IF (DIRFCK) THEN
         NMAT = NDMAT
         CALL MEMGET2('REAL','SKLTON',KSKLTN,NMAT*N2BASX,
     &      WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','FAC',   KFAC  ,NBASIS*(MAXREP+1),
     &      WORK,KFREE,LFREE)
         NINDEX = KMAX*(MAXREP + 1)*KHK(NHTYP)
         CALL MEMGET2('INTE','IINDEX',KINDEX,NINDEX,
     &      WORK,KFREE,LFREE)
         CALL MEMGET2('INTE','IPOINT',KPOINT,NBASIS*(MAXREP+1),
     &      WORK,KFREE,LFREE)
         IXYZ = 0
         CALL SKLFCK_1(FMAT,WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),
     &                 WORK(KPOINT),DUM,DUM,DUM,IDUM,IXYZ,
     &                 NMAT,ISYMDM,IFCTYP,IPRINT)
      END IF
C
C     Derivative Fock matrices
C     ========================
C
      IF (DDFOCK .AND. .NOT.LONDON) THEN
         IF (PERTUR) THEN
            NMATMX = NUCDEG(IATOME)
         ELSE
            NMATMX = 0
            DO 100 IATOM = 1, NUCIND
               NMATMX = MAX(NMATMX,NUCDEG(IATOM))
  100       CONTINUE
         END IF
C
         CALL MEMGET('REAL',KSKLTN,3*NMATMX*N2BASX  ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KFAC  ,NBASIS*(MAXREP+1),WORK,KFREE,LFREE)
         NINDEX = KMAX*(MAXREP + 1)*KHK(NHTYP)
         CALL MEMGET('INTE',KINDEX,NINDEX           ,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KPOINT,NBASIS*(MAXREP+1),WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KNFAC ,NMATMX*(MAXREP+1),WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KNINDX,(MAXREP+1)       ,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KNPINT,NMATMX*(MAXREP+1),WORK,KFREE,LFREE)
         IOFF = 0
         IFMAT = 0
         DO 200 ITYPE = 1, 2
         IF (ITYPE .EQ. 2 .AND. NODV) GO TO 200
            DO 300 IATOM = 1, NUCIND
            IF (.NOT.PERTUR .OR. IATOM.EQ.IATOME) THEN
               NMAT = NUCDEG(IATOM)
               CALL DCOPY(3*NMAT*N2BASX,FMAT(1,IOFF+1),1,WORK(KSKLTN),1)
               IOFFX = IOFF
               IOFFY = IOFF +   NMAT
               IOFFZ = IOFF + 2*NMAT
               DO 310 I = 1, NMAT
                  IADR = KSKLTN + 3*N2BASX*(I - 1)
                  CALL DCOPY(N2BASX,WORK(IADR),1,
     &                              FMAT(1,IOFFX + I),1)
                  CALL DCOPY(N2BASX,WORK(IADR+N2BASX),1,
     &                              FMAT(1,IOFFY + I),1)
                  CALL DCOPY(N2BASX,WORK(IADR+2*N2BASX),1,
     &                              FMAT(1,IOFFZ + I),1)
  310          CONTINUE
C     We ought also to reorder IFCTYP and ISYMDM, but as
C     all values are 1 and 0, resp., never mind. /941221-hjaaj
               DO 400 IXYZ = 1, 3
                  IADR = IOFF + NMAT*(IXYZ - 1) + 1
                  CALL SKLFCK_1(FMAT(1,IADR),WORK(KSKLTN),WORK(KFAC),
     &                        WORK(KINDEX),WORK(KPOINT),WORK(KNFAC),
     &                        WORK(KNINDX),WORK(KNPINT),IATOM,IXYZ,NMAT,
     &                        ISYMDM(IADR),IFCTYP(IADR),IPRINT)
                  CALL DSITSP(NBASIS,FMAT(1,IADR),WORK(KSKLTN))
C
                  IF (EXPFCK) THEN
                     IF (NODV) THEN
                        IREC = 3*(IATOM - 1) + IXYZ
                     ELSE
                        IREC = 6*(IATOM - 1) + 3*(ITYPE - 1) + IXYZ
                     END IF
                     CALL WRITDX(LUDFCK,LRDFCK,IREC,IRAT*NNBASX,
     &                           WORK(KSKLTN))
                  ELSE
                     IFMAT = IFMAT + 1
                     CALL DCOPY(NNBASX,WORK(KSKLTN),1,FMAT(1,IFMAT),1)
                  END IF
  400          CONTINUE
               IOFF = IOFF + 3*NMAT
            END IF
  300       CONTINUE
  200    CONTINUE
      END IF
C
C     Expectation values
C     ==================
C
      IF (EXPECT) THEN
        CALL MEMGET('REAL',KSKLTN,MXCOOR*2           ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFAC  ,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDEX,3*NUCIND*(MAXREP+1),WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KPOINT,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
         CALL SKLEXP(WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),WORK(KPOINT),
     &               MAXDER,IPRINT)
      END IF
C
C     London Fock matrics
C     ===================
C
      IF (DDFOCK .AND. LONDON) THEN
        CALL MEMGET('REAL',KSKLTN,3*N2BASX         ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFAC  ,NBASIS*(MAXREP+1),WORK,KFREE,LFREE)
        NINDEX = KMAX*(MAXREP + 1)*KHK(NHTYP)
        CALL MEMGET('INTE',KINDEX,NINDEX           ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KPOINT,NBASIS*(MAXREP+1),WORK,KFREE,LFREE)
C
        CALL SKLOND(FMAT(1,1),WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),
     &               WORK(KPOINT),IPRINT)
        IF (.NOT.NODV) THEN
           CALL SKLOND(FMAT(1,4),WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),
     &                  WORK(KPOINT),IPRINT)
        END IF
C
        NTYPE = 3
        IF (.NOT.NODV) NTYPE = 6
        DO 510 N = 1, NTYPE
           CALL DGETAP(NBASIS,FMAT(1,N),WORK(KSKLTN))
           CALL DCOPY(NNBASX,WORK(KSKLTN),1,FMAT(1,N),1)
  510   CONTINUE
      END IF
      CALL MEMREL('SKLFCK',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('SKLFCK')
      RETURN
      END
C  /* Deck SKLFCK_1 */
      SUBROUTINE SKLFCK_1(FMAT,SKLTON,FACTOR,INDEX,IPOINT,NFACTR,NINDEX,
     &                    NPOINT,IATOM,IXYZ,NMAT,ISYMDM,IFCTYP,IPRINT)
C*****************************************************************************
C
C Trygve Helgaker and Henrik Koch 93-94 ??
C
C Implemented all IFCTYPs Oct-Nov 94 Hans Joergen Aa. Jensen
C Revised Oct 11 1996 by T.Saue:
C   a) changed definition of IFCTYP
C   b) symmetrization/anti-symmetrization is taken out
C   c) operations containing inversion provided with minus sign
C      when acting on small components (RELCAL)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (D1 = 1.0D0,D2I = 0.5D0)
      PARAMETER (D4I = 0.25D0, D8I = 0.125D0, D16I = 0.0625D0)
      INTEGER A, B, AP, BP
      DIMENSION FMAT(NBASIS,NBASIS,*), SKLTON(NBASIS,NBASIS,*),
     &          INDEX(KMAX,0:MAXREP,*), FACTOR(NBASIS,0:MAXREP),
     &          IPOINT(NBASIS,0:MAXREP),
     &          NINDEX(0:MAXREP), NPOINT(NMAT,0:MAXREP),
     &          NFACTR(NMAT,0:MAXREP), ISYMDM(NMAT), IFCTYP(NMAT)
#include "nuclei.h"
#include "pincom.h"
#include "shells.h"
#include "symmet.h"
C
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Subroutine SKLFCK_1',-1)
         CALL HEADER("Input skeleton Fock matrices in SKLFCK_1",1)
         DO 10 IMAT = 1, NMAT
            WRITE (LUPRI,*) 'Input skeleton Fock matrix no.',IMAT,
     &         ' out of',NMAT
            CALL OUTPUT(FMAT(1,1,IMAT),1,NBASIS,1,NBASIS,
     &                  NBASIS,NBASIS,-1,LUPRI)
   10    CONTINUE
      END IF
C
C     Finish skeleton Fock matrices
C     ===================================
C     CTROND: symmetrization/anti-symmetrization is removed !!!!
C     (b) multiply with factors omitted in INTFCK
C         (see comments in INTFCK)
C     /941019-hjaaj
C
      IF (IPRINT .GT. 5) THEN
        CALL HEADER("Final scaled skeleton Fock matrices in SKLFCK_1",1)
      END IF
      DO 100 IMAT = 1, NMAT
         IY = MOD(IFCTYP(IMAT),10)
         IX = IFCTYP(IMAT)/10
         FAC = D2I
         IF(IX.EQ.0) FAC = D4I
         CALL DSCAL(NBASIS*NBASIS,FAC,FMAT(1,1,IMAT),1)
         IF (IPRINT .GT. 5) THEN
            WRITE (LUPRI,*) 'Skeleton Fock matrix no.',IMAT,
     &         ' out of',NMAT
            CALL OUTPUT(FMAT(1,1,IMAT),1,NBASIS,1,NBASIS,
     &                  NBASIS,NBASIS,-1,LUPRI)
         END IF
  100 CONTINUE
C
C     Form AO Fock matrix from symmetric skeleton matrix
C     ==================================================
C
      IF (MAXREP .GT. 0) THEN
C
C        First construct index and symmetry factor arrays
C
         IORB = 0
         DO 200 ISHELL = 1, KMAX
            DO 210 ISYMOP = 0, MAXREP
            IF (IAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 220 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  INDEX(ISHELL,ISYMOP,ICOMP) = IORB
  220          CONTINUE
            ELSE
               ISYMF = IEOR(ISYMOP,IAND(ISYMOP,ISTBAO(ISHELL)))
               DO 230 ICOMP  = 1, KHKT(ISHELL)
                  INDEX(ISHELL,ISYMOP,ICOMP) = INDEX(ISHELL,ISYMF,ICOMP)
  230          CONTINUE
            END IF
  210       CONTINUE
  200    CONTINUE
C
         IORB = 0
         DO 300 ISHELL = 1, NLRGSH
            NHKTA = NHKT(ISHELL)
            DO 310 ISYMOP = 0, MAXREP
            IF (IAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 320 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  DO 330 ISYM = 0, MAXREP
                     ISYMPR = IEOR(ISYMOP,ISYM)
                     IPOINT(IORB,ISYM) = INDEX(ISHELL,ISYMPR,ICOMP)
                     FACTOR(IORB,ISYM) =
     &                        PT(IAND(ISYM,ISYMAO(NHKTA,ICOMP)))
  330             CONTINUE
  320          CONTINUE
            END IF
  310       CONTINUE
  300    CONTINUE
#ifdef PRG_DIRAC
C
C        The next paragraph concerns relativistic calculations only !
C        IPARity is invoked since the operation of inversion
C        receives a minus sign when acting on small components....
C
         DO 340 ISHELL = NLRGSH+1,KMAX
            NHKTA = NHKT(ISHELL)
            DO 350 ISYMOP = 0, MAXREP
            IF (IAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 360 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  DO 370 ISYM = 0, MAXREP
                     ISYMPR = IEOR(ISYMOP,ISYM)
                     IPOINT(IORB,ISYM) = INDEX(ISHELL,ISYMPR,ICOMP)
                     FACTOR(IORB,ISYM) = IPAR(ISYM)*
     &                        PT(IAND(ISYM,ISYMAO(NHKTA,ICOMP)))
  370             CONTINUE
  360          CONTINUE
            END IF
  350       CONTINUE
  340    CONTINUE
#endif
C
C        NPOINT and NFACTR
C
         IF (IXYZ .GT. 0) THEN
            IMAT = 0
            DO 215 ISYMOP = 0, MAXREP
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               IMAT = IMAT + 1
               NINDEX(ISYMOP) = IMAT
            ELSE
               ISYMF = IEOR(ISYMOP,IAND(ISYMOP,ISTBNU(IATOM)))
               NINDEX(ISYMOP) = NINDEX(ISYMF)
            END IF
  215       CONTINUE
C
            IMAT = 0
            DO 315 ISYMOP = 0, MAXREP
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               IMAT = IMAT + 1
               DO 335 ISYM = 0, MAXREP
                  ISYMPR = IEOR(ISYMOP,ISYM)
                  NPOINT(IMAT,ISYM) = NINDEX(ISYMPR)
                  NFACTR(IMAT,ISYM) = PT(IAND(ISYM,ISYMAX(IXYZ,1)))
  335          CONTINUE
            END IF
  315       CONTINUE
         END IF
C
C        Construct full matrix from skeleton
C
         CALL DCOPY(NMAT*(NBASIS**2),FMAT,1,SKLTON,1)
         IF (IXYZ .EQ. 0) THEN
            DO 400 ISYMOP = 1, MAXREP
               DO 435 N = 1, NMAT
                  FP = PT(IAND(ISYMOP,ISYMDM(N)))
                  DO 420 B = 1, NBASIS
                     BP = IPOINT(B,ISYMOP)
                     FBFP = FACTOR(B,ISYMOP)*FP
                     DO 410 A = 1, NBASIS
                        AP = IPOINT(A,ISYMOP)
                        FA = FACTOR(A,ISYMOP)
                        FMAT(AP,BP,N) = FMAT(AP,BP,N)
     &                                + FA*FBFP*SKLTON(A,B,N)
  410                CONTINUE
  420             CONTINUE
  435          CONTINUE
  400       CONTINUE
         ELSE
            DO 430 ISYMOP = 1, MAXREP
               DO 440 N = 1, NMAT
                  NP = NPOINT(N,ISYMOP)
                  FP = NFACTR(N,ISYMOP)
                  DO 450 A = 1, NBASIS
                     AP = IPOINT(A,ISYMOP)
                     FA = FACTOR(A,ISYMOP)
                     DO 460 B = 1, NBASIS
                        BP = IPOINT(B,ISYMOP)
                        FB = FACTOR(B,ISYMOP)
                        FMAT(AP,BP,NP) = FMAT(AP,BP,NP)
     &                                 + FA*FB*FP*SKLTON(A,B,N)
  460                CONTINUE
  450             CONTINUE
  440          CONTINUE
  430       CONTINUE
         END IF
C
         FAC = D1/dble(MAXREP + 1)
         CALL DSCAL(NMAT*(NBASIS**2),FAC,FMAT,1)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("AO Fock matrices in SKLFCK_1",1)
            DO 40 N = 1, NMAT
               WRITE (LUPRI,*) 'AO Fock matrix no.',N,
     &         ' out of',NMAT
               CALL OUTPUT(FMAT(1,1,N),1,NBASIS,1,NBASIS,
     &                     NBASIS,NBASIS,-1,LUPRI)
   40       CONTINUE
         END IF
C
C
#ifndef PRG_DIRAC
C
CCTOND:     Matrices must be symmetrized/anti-symmetrized first !
C
C
C        Transformation to symmetry basis
C        ================================
C
         CALL DCOPY(NMAT*NBASIS**2,FMAT,1,SKLTON,1)
         CALL DZERO(FMAT,NMAT*NBASIS*NBASIS)
         IF (IXYZ .EQ. 0) THEN
            ISTRA = 1
            DO 500 IREPA = 0, MAXREP
               NORBA = NAOS(IREPA+1)
               DO 510 I = ISTRA,ISTRA + NORBA - 1
                  IA   = IAND(ISHFT(IPIND(I),-16),65535)
                  NA   = IAND(ISHFT(IPIND(I), -8),  255)
                  NHKTA  = NHKT(IA)
                  KHKTA  = KHKT(IA)
                  MULA   = ISTBAO(IA)
                  INDA   = KSTRT(IA) + NA - KHKTA
                  DO 520 ISYMA = 0, MAXOPR
                  IF (IAND(ISYMA,MULA) .EQ. 0) THEN
                     INDA = INDA + KHKTA
                     FACA = PT(IAND(ISYMA,
     &                                IEOR(IREPA,ISYMAO(NHKTA,NA))))
                     DO 550 N = 1, NMAT
                        ISTRB = 1
                        DO 525 IREPB = 0, MAXREP
                           NORBB = NAOS(IREPB+1)
                           IF (IEOR(IREPA,IREPB) .EQ. ISYMDM(N)) THEN
                           DO 530 J = ISTRB,ISTRB + NORBB - 1
                              IB   = IAND(ISHFT(IPIND(J),-16),65535)
                              NB   = IAND(ISHFT(IPIND(J), -8),  255)
                              NHKTB  = NHKT(IB)
                              KHKTB  = KHKT(IB)
                              MULB   = ISTBAO(IB)
                              INDB   = KSTRT(IB) + NB - KHKTB
                              DO 540 ISYMB = 0, MAXOPR
                              IF (IAND(ISYMB,MULB) .EQ. 0) THEN
                                 INDB = INDB + KHKTB
                                 FACB = PT(IAND(ISYMB,
     &                                  IEOR(IREPB,ISYMAO(NHKTB,NB))))
                                 FMAT(I,J,N) = FMAT(I,J,N)
     &                                 + FACA*FACB*SKLTON(INDA,INDB,N)
                              END IF
540                           CONTINUE
530                        CONTINUE
                           END IF
                           ISTRB = ISTRB + NORBB
  525                   CONTINUE
  550                CONTINUE
                  END IF
520               CONTINUE
510            CONTINUE
               ISTRA = ISTRA + NORBA
500         CONTINUE
         ELSE
            MULN = ISTBNU(IATOM)
            ISTRA = 1
            DO 600 IREPA = 0, MAXREP
               NORBA = NAOS(IREPA+1)
               DO 610 I = ISTRA,ISTRA + NORBA - 1
                  IA   = IAND(ISHFT(IPIND(I),-16),65535)
                  NA   = IAND(ISHFT(IPIND(I), -8),  255)
                  NHKTA = NHKT(IA)
                  KHKTA = KHKT(IA)
                  MULA  = ISTBAO(IA)
                  INDA  = KSTRT(IA) + NA - KHKTA
                  IVARA = IEOR(IREPA,ISYMAO(NHKTA,NA))
                  DO 620 ISYMA = 0, MAXOPR
                  IF (IAND(ISYMA,MULA) .EQ. 0) THEN
                     FA    = PT(IAND(ISYMA,IVARA))
                     INDA  = INDA + KHKTA
                     ISTRB = 1
                     DO 630 IREPB = 0, MAXREP
                        IVARN=IEOR(IEOR(IREPA,IREPB),ISYMAX(IXYZ,1))
                        NORBB = NAOS(IREPB+1)
                        IF (IAND(MULN,IVARN).EQ.0) THEN
                           DO 640 J = ISTRB,ISTRB + NORBB - 1
                              IB   = IAND(ISHFT(IPIND(J),-16),65535)
                              NB   = IAND(ISHFT(IPIND(J), -8),  255)
                              NHKTB = NHKT(IB)
                              KHKTB = KHKT(IB)
                              MULB  = ISTBAO(IB)
                              INDB  = KSTRT(IB) + NB - KHKTB
                              IVARB = IEOR(IREPB,ISYMAO(NHKTB,NB))
                              DO 650 ISYMB = 0, MAXOPR
                              IF (IAND(ISYMB,MULB) .EQ. 0) THEN
                                 INDB = INDB + KHKTB
                                 FAB = FA*PT(IAND(ISYMB,IVARB))
                                 DO 660 ISYMN = 0, MAXOPR
                                 IF (IAND(ISYMN,MULN) .EQ. 0) THEN
                                    INDN = NINDEX(ISYMN)
                                    FAC  = FAB*PT(IAND(ISYMN,IVARN))
                                    FMAT(I,J,1) = FMAT(I,J,1)
     &                                      + FAC*SKLTON(INDA,INDB,INDN)
                                 END IF
  660                            CONTINUE
                              END IF
  650                         CONTINUE
  640                      CONTINUE
                           END IF
                        ISTRB = ISTRB + NORBB
  630                CONTINUE
                  END IF
  620             CONTINUE
  610          CONTINUE
               ISTRA = ISTRA + NORBA
  600       CONTINUE
         END IF
#endif
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("SO Fock matrices in SKLFCK_1",1)
            DO 700 N = 1, NMAT
               WRITE (LUPRI,*) 'SO Fock matrix no.',N,
     &         ' out of',NMAT
               CALL OUTPUT(FMAT(1,1,N),1,NBASIS,1,NBASIS,NBASIS,NBASIS,
     &                     -1,LUPRI)
  700       CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck sklexp */
      SUBROUTINE SKLEXP(SKLTON,FACTOR,INDEX,IPOINT,MAXDER,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      INTEGER A, B, AP, BP
      DIMENSION SKLTON(MXCOOR,MXCOOR),
     &          INDEX(NUCIND,3,0:MAXREP), FACTOR(3*NUCDEP,0:MAXREP),
     &          IPOINT(3*NUCDEP,0:MAXREP)
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
#include "nuclei.h"
#include "symmet.h"
C
      IF (IPRINT .GT. 5) CALL HEADER('Subroutine SKLEXP',-1)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER("Input skeleton molecular gradient",1)
         CALL OUTPUT(GRADEE,1,1,1,3*NUCDEP,1,MXCOOR,-1,LUPRI)
      END IF
      IF (MAXDER .EQ. 2 .AND. IPRINT .GT. 5) THEN
         CALL HEADER("Input skeleton molecular Hessian",1)
         CALL OUTPUT(HESSEE,1,3*NUCDEP,1,3*NUCDEP,MXCOOR,MXCOOR,
     &               -1,LUPRI)
      END IF
C
#ifndef PRG_DIRAC
      IF (MAXDER .EQ. 2) THEN
         DO 10 I = 1, 3*NUCDEP
         DO 10 J = 1, I - 1
            HESSEE(I,J) = HESSEE(I,J) + HESSEE(J,I)
            HESSEE(J,I) = HESSEE(I,J)
   10    CONTINUE
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("Skeleton molecular Hessian symmetrized",1)
            CALL OUTPUT(HESSEE,1,3*NUCDEP,1,3*NUCDEP,MXCOOR,MXCOOR,
     &                  -1,LUPRI)
         END IF
      END IF
#endif
C
      IF (MAXREP .GT. 0) THEN
         ICOOR = 0
         DO 100 IATOM = 1, NUCIND
            DO 110 ISYMOP = 0, MAXREP
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               DO 120 ICART = 1, 3
                  ICOOR = ICOOR + 1
                  INDEX(IATOM,ICART,ISYMOP) = ICOOR
  120          CONTINUE
            ELSE
               ISYMF = IEOR(ISYMOP,IAND(ISYMOP,ISTBNU(IATOM)))
               DO 130 ICART = 1, 3
                  INDEX(IATOM,ICART,ISYMOP) = INDEX(IATOM,ICART,ISYMF)
  130          CONTINUE
            END IF
  110       CONTINUE
  100    CONTINUE
C
         ICOOR = 0
         DO 200 IATOM = 1, NUCIND
            DO 210 ISYMOP = 0, MAXREP
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               DO 220 ICART = 1, 3
                  ICOOR = ICOOR + 1
                  DO 230 ISYM = 0, MAXREP
                     ISYMPR = IEOR(ISYMOP,ISYM)
                     IPOINT(ICOOR,ISYM)=INDEX(IATOM,ICART,ISYMPR)
                     FACTOR(ICOOR,ISYM)=PT(IAND(ISYM,ISYMAX(ICART,1)))
  230             CONTINUE
  220          CONTINUE
            END IF
  210       CONTINUE
  200    CONTINUE
C
C        Construct full gradient from skeleton
C
         CALL DCOPY(MXCOOR,GRADEE,1,SKLTON,1)
         DO 300 ISYMOP = 1, MAXREP
            DO 310 A = 1, 3*NUCDEP
               AP = IPOINT(A,ISYMOP)
               GRADEE(AP) = GRADEE(AP) + FACTOR(A,ISYMOP)*SKLTON(A,1)
  310       CONTINUE
  300    CONTINUE
         FAC = D1/dble(MAXREP + 1)
         CALL DSCAL(MXCOOR,   FAC,GRADEE,1)
#ifndef PRG_DIRAC
C
C        Construct full Hessian from skeleton
C
         IF (MAXDER .EQ. 2) THEN
            CALL DCOPY(MXCOOR**2,HESSEE,1,SKLTON,1)
            DO 400 ISYMOP = 1, MAXREP
               DO 410 A = 1, 3*NUCDEP
               DO 410 B = 1, 3*NUCDEP
                  AP  = IPOINT(A,ISYMOP)
                  BP  = IPOINT(B,ISYMOP)
                  FAC = FACTOR(A,ISYMOP)*FACTOR(B,ISYMOP)
                  HESSEE(AP,BP) = HESSEE(AP,BP) + FAC*SKLTON(A,B)
  410          CONTINUE
  400       CONTINUE
            FAC = D1/dble(MAXREP + 1)
            CALL DSCAL(MXCOOR**2,FAC,HESSEE,1)
         END IF
#endif
C
         IF (IPRINT .GT. 5) THEN
           CALL HEADER("Full molecular gradient (non-symmetry basis)",1)
           CALL OUTPUT(GRADEE,1,1,1,3*NUCDEP,1,MXCOOR,-1,LUPRI)
         END IF
#ifndef PRG_DIRAC
         IF (MAXDER .EQ. 2 .AND. IPRINT .GT. 5) THEN
           CALL HEADER("Full molecular Hessian (non-symmetry basis)",1)
           CALL OUTPUT(HESSEE,1,3*NUCDEP,1,3*NUCDEP,MXCOOR,MXCOOR,
     &                 -1,LUPRI)
         END IF
#endif
C
C        Transformation to symmetry basis
C        ================================
C
         CALL DCOPY(MXCOOR,GRADEE,1,SKLTON,1)
         CALL DZERO(GRADEE,MXCOOR)
         ICOOR = 0
         DO 500 IATOM = 1, NUCIND
         DO 500 ICART = 1, 3
            IF (IAND(ISTBNU(IATOM),ISYMAX(ICART,1)).EQ.0) THEN
               ICOOR = ICOOR + 1
               GRD = D0
               DO 510 ISYMOP = 0, MAXOPR
               IF (IAND(ISTBNU(IATOM),ISYMOP).EQ.0) THEN
                  FAC = PT(IAND(ISYMOP,ISYMAX(ICART,1)))
                  GRD = GRD + FAC*SKLTON(INDEX(IATOM,ICART,ISYMOP),1)
               END IF
  510          CONTINUE
               GRADEE(ICOOR) = GRD
            END IF
  500    CONTINUE
C
#ifndef PRG_DIRAC
         IF (MAXDER .EQ. 2) THEN
            CALL DCOPY(MXCOOR**2,HESSEE,1,SKLTON,1)
            CALL DZERO(HESSEE,MXCOOR**2)
            A = 0
            DO 600 IREPA  = 0, MAXREP
            DO 600 IATOMA = 1, NUCIND
            DO 600 ICARTA = 1, 3
               IVARA = IEOR(IREPA,ISYMAX(ICARTA,1))
               IF (IAND(ISTBNU(IATOMA),IVARA).EQ.0) THEN
                  A = A + 1
                  B = 0
                  DO 700 IREPB  = 0, MAXREP
                  DO 700 IATOMB = 1, NUCIND
                  DO 700 ICARTB = 1, 3
                     IVARB = IEOR(IREPB,ISYMAX(ICARTB,1))
                     IF (IAND(ISTBNU(IATOMB),IVARB).EQ.0) THEN
                        B = B + 1
                        HES = D0
                        DO 800 ISYMA = 0, MAXOPR
                        IF (IAND(ISTBNU(IATOMA),ISYMA).EQ.0) THEN
                           FA = PT(IAND(IVARA,ISYMA))
                           AP = INDEX(IATOMA,ICARTA,ISYMA)
                           DO 810 ISYMB = 0, MAXOPR
                           IF (IAND(ISTBNU(IATOMB),ISYMB).EQ.0) THEN
                              FB = PT(IAND(IVARB,ISYMB))
                              BP = INDEX(IATOMB,ICARTB,ISYMB)
                              HES = HES + FA*FB*SKLTON(AP,BP)
                           END IF
  810                      CONTINUE
                        END IF
  800                   CONTINUE
                        HESSEE(A,B) = HES
                     END IF
  700             CONTINUE
               END IF
  600       CONTINUE
         END IF
#endif
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER("Symmetrized molecular gradient",1)
         CALL PRIGRD(GRADEE)
      END IF
#ifndef PRG_DIRAC
      IF (MAXDER .EQ. 2 .AND. IPRINT .GT. 5) THEN
         CALL HEADER("Symmetrized molecular Hessian",1)
         CALL PRIHES(HESSEE,'CENTERS')
      END IF
#endif
      RETURN
      END
C  /* Deck sklond */
      SUBROUTINE SKLOND(FMAT,SKLTON,FACTOR,INDEX,IPOINT,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0)
      INTEGER A, B, AP, BP, X
      DIMENSION FMAT(NBASIS,NBASIS,3), SKLTON(NBASIS,NBASIS,3),
     &          INDEX(KMAX,0:MAXREP,*), FACTOR(NBASIS,0:MAXREP),
     &          IPOINT(NBASIS,0:MAXREP)
#include "nuclei.h"
#include "pincom.h"
#include "shells.h"
#include "symmet.h"
C
      IF (IPRINT .GT. 5) CALL HEADER('Subroutine SKLOND',-1)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER("Input skeleton Fock matrices in SKLOND",1)
         DO 100 X = 1, 3
            WRITE(LUPRI,*) 'Input skeleton Fock matrix no.',X,NFMAT
            CALL OUTPUT(FMAT(1,1,X),1,NBASIS,1,NBASIS,
     &                  NBASIS,NBASIS,-1,LUPRI)
  100    CONTINUE
      END IF
C
C     Form symmetric skeleton Fock matrix
C     ===================================
C
      DO 200 X = 1, 3
         DO 210 I = 1, NBASIS
         DO 210 J = 1, I
            FMATIJ = DP5*(FMAT(I,J,X) - FMAT(J,I,X))
            FMAT(I,J,X) =   FMATIJ
            FMAT(J,I,X) = - FMATIJ
  210    CONTINUE
  200 CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER("Symmetric skeleton Fock matrix in SKLOND",1)
         DO 300 X = 1, 3
            WRITE(LUPRI,*) 'Symmetric skeleton Fock matrix no.',X,NFMAT
            CALL OUTPUT(FMAT(1,1,X),1,NBASIS,1,NBASIS,NBASIS,NBASIS,
     &                  -1,LUPRI)
  300    CONTINUE
      END IF
C
C     Form AO Fock matrix from symmetric skeleton matrix
C     ==================================================
C
      IF (MAXREP .GT. 0) THEN
C
C        First construct index and symmetry factor arrays
C
         IORB = 0
         DO 400 ISHELL = 1, KMAX
            DO 410 ISYMOP = 0, MAXREP
            IF (IAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 420 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  INDEX(ISHELL,ISYMOP,ICOMP) = IORB
  420          CONTINUE
            ELSE
               ISYMF = IEOR(ISYMOP,IAND(ISYMOP,ISTBAO(ISHELL)))
               DO 430 ICOMP  = 1, KHKT(ISHELL)
                  INDEX(ISHELL,ISYMOP,ICOMP) = INDEX(ISHELL,ISYMF,ICOMP)
  430          CONTINUE
            END IF
  410       CONTINUE
  400    CONTINUE
C
         IORB = 0
         DO 500 ISHELL = 1, KMAX
            NHKTA = NHKT(ISHELL)
            DO 510 ISYMOP = 0, MAXREP
            IF (IAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 520 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  DO 530 ISYM = 0, MAXREP
                     ISYMPR = IEOR(ISYMOP,ISYM)
                     IPOINT(IORB,ISYM) = INDEX(ISHELL,ISYMPR,ICOMP)
                     FACTOR(IORB,ISYM) =
     &                        PT(IAND(ISYM,ISYMAO(NHKTA,ICOMP)))
  530             CONTINUE
  520          CONTINUE
            END IF
  510       CONTINUE
  500    CONTINUE
C
C        Construct full matrix from skeleton
C
         CALL DCOPY(3*(NBASIS**2),FMAT,1,SKLTON,1)
         DO 600 ISYMOP = 1, MAXREP
            DO 610 X = 1, 3
               FP = PT(IAND(ISYMOP,ISYMAX(X,2)))
               DO 620 A = 1, NBASIS
                  AP = IPOINT(A,ISYMOP)
                  FA = FACTOR(A,ISYMOP)
                  DO 630 B = 1, NBASIS
                     BP = IPOINT(B,ISYMOP)
                     FB = FACTOR(B,ISYMOP)
                     FMAT(AP,BP,X)=FMAT(AP,BP,X)+FA*FB*FP*SKLTON(A,B,X)
  630             CONTINUE
  620          CONTINUE
  610       CONTINUE
  600    CONTINUE
C
         FAC = D1/dble(MAXREP + 1)
         CALL DSCAL(3*(NBASIS**2),FAC,FMAT,1)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("AO Fock matrices in SKLOND",1)
            DO 700 X = 1, 3
               WRITE(LUPRI,*) 'AO skeleton Fock matrix no.',X,NFMAT
               CALL OUTPUT(FMAT(1,1,X),1,NBASIS,1,NBASIS,
     &                     NBASIS,NBASIS,-1,LUPRI)
  700       CONTINUE
         END IF
C
C        Transformation to symmetry basis
C        ================================
C
         CALL DCOPY(3*(NBASIS**2),FMAT,1,SKLTON,1)
         CALL DZERO(FMAT,3*NBASIS*NBASIS)
         ISTRA = 1
         DO 800 IREPA = 0, MAXREP
            NORBA = NAOS(IREPA+1)
            DO 810 I = ISTRA,ISTRA + NORBA - 1
               IA    = IAND(ISHFT(IPIND(I),-16),65535)
               NA    = IAND(ISHFT(IPIND(I), -8),  255)
               NHKTA = NHKT(IA)
               KHKTA = KHKT(IA)
               MULA  = ISTBAO(IA)
               INDA  = KSTRT(IA) + NA - KHKTA
               IVARA = IEOR(IREPA,ISYMAO(NHKTA,NA))
               DO 820 ISYMA = 0, MAXOPR
               IF (IAND(ISYMA,MULA) .EQ. 0) THEN
                  FA    = PT(IAND(ISYMA,IVARA))
                  INDA  = INDA + KHKTA
                  ISTRB = 1
                  DO 830 IREPB = 0, MAXREP
                     NORBB = NAOS(IREPB+1)
                     DO 840 J = ISTRB,ISTRB + NORBB - 1
                        IB   = IAND(ISHFT(IPIND(J),-16),65535)
                        NB   = IAND(ISHFT(IPIND(J), -8),  255)
                        NHKTB = NHKT(IB)
                        KHKTB = KHKT(IB)
                        MULB  = ISTBAO(IB)
                        INDB  = KSTRT(IB) + NB - KHKTB
                        IVARB = IEOR(IREPB,ISYMAO(NHKTB,NB))
                        DO 850 ISYMB = 0, MAXOPR
                        IF (IAND(ISYMB,MULB) .EQ. 0) THEN
                           INDB = INDB + KHKTB
                           FAB  = FA*PT(IAND(ISYMB,IVARB))
                           DO 860 X = 1, 3
                           IF (IEOR(IREPA,IREPB).EQ.ISYMAX(X,2)) THEN
                              FMAT(I,J,X) = FMAT(I,J,X)
     &                                    + FAB*SKLTON(INDA,INDB,X)
                           END IF
  860                      CONTINUE
                        END IF
  850                   CONTINUE
  840                CONTINUE
                     ISTRB = ISTRB + NORBB
  830             CONTINUE
               END IF
  820          CONTINUE
  810       CONTINUE
            ISTRA = ISTRA + NORBA
  800    CONTINUE
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("SO Fock matrix in SKLOND",1)
            DO 900 X = 1, 3
               CALL OUTPUT(FMAT(1,1,X),1,NBASIS,1,NBASIS,
     &                     NBASIS,NBASIS,-1,LUPRI)
  900       CONTINUE
         END IF
      END IF
      RETURN
      END
