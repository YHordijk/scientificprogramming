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

C FILE    : diropen.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MKMOFK */
      SUBROUTINE MKMOFK(FMO,FAO,DMAT,TMAT,CMO,TDINV,DTEMP1,DTEMP2,
     &                  WORK,LWORK)
C****************************************************************************
C
C     PURPOSE:
C  =============
C        Construct MO Fock matrix from Fock matrices in AO-basis
C
C     On input:
C  FAO(*,1)       - FD-matrix, that is FD = h + {sum(S) Q_s}, see (3.41)
C                   of the thesis of J. Thyssen
C  FAO(*,1+IOPEN) - FV-matrix, two-electron Fock matrix constructed from
C                   the normalized density matrix of open shell IOPEN
C  TMAT
C  DMAT
C  CMO
C  TDINV
C
C     On output:
C  FMO
C
C        Written by J. Thyssen
C        Last revision: Dec 30 1997
C          Nov.2005 - MI - only print-out info
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( DP5 = 0.50D00 , D0 = 0.00D00, D1 = 1.00D00,
     &            DM1 = -1.00D00, D2 = 2.00D00, DFTHR = 1.0D-12)
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION FMO(*),FAO(N2BBASXQ,NFMAT),TMAT(*),CMO(*),
     &          DTEMP1(N2BBASXQ,*),DTEMP2(N2BBASXQ),
     &          DMAT(N2BBASXQ,NFMAT),TDINV(*),WORK(*)
      CHARACTER CTMP*3
C
#include "memint.h"
C
      CALL QENTER('MKMOFK')
C
C     The various blocks of the Fock matrix in MO-basis is defined in (3.40)
C     They are constructed in AO-basis using density matrices as projection
C     operators.
C
C     In AO-basis the Fock matrix can be expressed as in (3.41) where
C     the first two terms correspond to FD and the next two terms
C     are the coupling operator L and its Hermitian conjugate.
C
C     The coupling operator L is defined in (3.42) with the active-active
C     part between two different open shells given in Eqs. (3.43) - (3.45)
C
C     One or more open shells that do NOT couple:
C     -------------------------------------------
C
C     Note: t stands for hermitian conjungation
C     Fmo = FDmo + DVmo * FVmo * (DCmo - DSmo)t + H.C.
C         = FDmo + Ct Sao [DVao FVao (beta DCao + gamma DSao)t] Sao C + h.c.
C     where beta  = alpha*f
C     and   gamma = a-1
C
C     Eg. for a single open shell (NACSH = 1, NAELEC = 1)
C        beta  = 1 and
C        gamma = -1
C
C     Open shells that do couple:
C     ---------------------------
C
C     See thesis of J. Thyssen ( can be found on http://www.diracprogram.org )
C
C     Transform FD to MO-basis
C     ========================
C
      DO I = 1,NFSYM
      IF(NTMO(I).GT.0) THEN
        IF(SUB_BL) THEN
          CALL TSUBBL(I,FAO,FMO,TMAT,TMAT,WORK(KFREE),LFREE,IPRSCF)
        ELSE
          CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &            NTMO(I),NTMO(I),
     &            FAO(I2BASX(I,I)+1,1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &            TMAT(1+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &            TMAT(1+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRSCF)
        ENDIF

        IF(IPRSCF.GE.10) THEN
         CALL HEADER(
     &    'MKMOFK: Entering Fock AO matrix - closed shell part',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')
     &    '* Fermion ircop no.',I,'/',NFSYM
         CALL PRQMAT(FAO(I2BASX(I,I)+1,1),NFBAS(I,0),NFBAS(I,0),
     &          NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'MKMOFK: Entering AO2MO transformation matrix,TMAT(NZT)',-1)
         WRITE(LUPRI,'(A,I1,A,I1)')
     &   '* Fermion ircop no.',I,'/',NFSYM
         CALL PRQMAT(TMAT(1+I2TMT(I)),NFBAS(I,0),NTMO(I),
     &           NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &    'MKMOFK: D Fock matrix transformed into MO-basis',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')
     &    '* Fermion ircop no.',I,'/',NFSYM
         CALL PRQMAT(FMO(1+I2TMOT(I)),NTMO(I),NTMO(I),
     &          NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        END IF
      END IF
      END DO
C
C     If not .trivec. then we don't have coefficients to form a density
C     matrix, so return.
C
      IF (.NOT. TRIVEC) THEN
C
C        In next iteration we'll have some coefficients
C
         TRIVEC = .TRUE.
         GOTO 30
      END IF
C
C     Construct coupling operator L
C     =============================
C
      IF(.NOT.AOC) GOTO 30
C
      DO 10 IOPEN = 1,NOPEN
C
         BETA  = DALPHA(IOPEN)*DF(IOPEN)
         GAMMA = DA(IOPEN)-D1
         IF (IPRSCF.GE.7) THEN
           WRITE(LUPRI,'(A,I2,A,F10.6)')
     &     'MKMOFK: BETA  factor, open shell',IOPEN,' : ',BETA,
     &     'MKMOFK: GAMMA factor, open shell',IOPEN,' : ',GAMMA,
     &     'MKMOFK: .OPENFAC    , open shell',IOPEN,' : ',abs(OPEN_FAC)
         END IF
C
C        Get secondary density matrices as projection operator
C        in AO-basis
C        IFLAG = 4 + 8 = 12
C          4 - 000100: Secondary electron density matrix
C          8 - 001000: Secondary positron density matrix
C
         IF(VACPOL) THEN
C..........negative energy solutions are occupied
           IFLAG = 4
         ELSE
           IFLAG = 12
         ENDIF
         CALL GENDEN(DTEMP1,CMO,IFLAG,IPRSCF)
C        This is the inactive projector P^Z of (3.42)
C
C        Construct the projection operator:
C         BETA DCao + GAMMA DSao
C           = alpha*f DCao + (a_o-1)*DSao
C
         CALL DSCAL(N2BBASXQ,GAMMA,DTEMP1(1,1),1) ! (a_U-1)*P^Z
C        Next we add the contribution from the inactive projector
C        alpha_U f_U P^I to the secondary projector
         IF(VACPOL) THEN
C..........negative energy solutions are occupied
           IFLAG = 9
           CALL GENDEN(DTEMP1(1,2),CMO,IFLAG,IPRSCF)
           CALL DAXPY(N2BBASXQ,BETA,DTEMP1(1,2),1,DTEMP1(1,1),1)
         ELSE
           CALL DAXPY(N2BBASXQ,BETA,DMAT(1,1),1,DTEMP1(1,1),1)
         ENDIF
C        Finally we add the contribution from the active projector
C        of this open shell to complete the first term of (3.42)
C        From (3.42) it should be noted that this projector is used
C        to generate the active-active block with U=V.
C        Since this diagonal block is also generated by the
C        Hermitian conjugate of the coupling operator, we must multiply
C        with a factor half to avoid double counting
         FAC = GAMMA*DP5 * abs(OPEN_FAC)
         IF (FAC.NE.D0)
     &      CALL DAXPY(N2BBASXQ,FAC,DMAT(1,1+IOPEN),1,DTEMP1(1,1),1)
C        Note that there is freedom in the choice of diagonal blocks
C        of the Fock matrix. The choice for active-active diagonal
C        blocks of (3.40) allows identification of the resulting
C        orbital eigenvalues as ionization energies, by extension
C        of Koopmans theorem to the average-of-configuration case.
C
C        We now have the complete right-hand projector of the first
C        term of (3.42) for shell U contained in DTEMP1
C
         IOFFTDINV = 1
         DO 20 I = 1,NFSYM
C
            IF (NACSH(I,IOPEN).GT.0) THEN
C
C              We next multiply the right-hand projector with Q^U
C              and store the result in DTEMP1(1,2)
C              QVao (beta DCao + gamma DSao) (store in DTEMP1(1,2))
C
C              QVao = 2*f/2 * Fock(TWOINT) = f * Fock(TWOINT)
C
               CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),DF(IOPEN),
     &            'N','N',IPQTOQ(1,0),FAO(1+I2BASX(I,I),1+IOPEN),
     &            NTBAS(0),NTBAS(0),NZ,
     &            'N','N',IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),1),
     &            NTBAS(0),NTBAS(0),NZ,
     &            D0,IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),2),
     &            NTBAS(0),NTBAS(0),NZ)
C
C              Now we start constructing the active-active part
C              of the coupling operator, that is the second term
C              of (3.42)
C
               DO 25 JOPEN = IOPEN+1,NOPEN
                  IF (NACSH(I,JOPEN).EQ.0) GOTO 25
C
C                 We use DTEMP1(xxx,3) as temporary storage
C
C                 We walk through the three cases described in (3.43),
C                 (3.44) and (3.45) by comparing fractional occupations
C                 of open shells I and J, that is DF(IOPEN) vs. DF(JOPEN)
C
                  CALL DZERO(DTEMP1(1,3),N2BBASXQ)
C
                  IF (ABS(DF(IOPEN)-DF(JOPEN)).LT.DFTHR) THEN
C
C                    This is a tough case...
C                    The gradient is correct, but there might
C                    be serious trouble converging the wave-function,
C                    because differences of diagonal elements are not
C                    proportional to the orbital hessian.
C
C                    If fi = fj = 1 then it corresponds to inactive shells,
C                    corresponding to case 3, that is, (3.45),
C                    so don't make the correction below.
C
C                    Otherwise, we have case 2 and (3.44)
                     IF (ABS(DF(IOPEN)-D1).GE.DFTHR) THEN
C                       First we subtract FD
                        CALL DAXPY(N2BBASXQ,-D1,FAO,1,DTEMP1(1,3),1)
C
C                       add two-electron Fock matrices...
C
C                       ( + (ai-1) fi QViao
C                         - (aj-1) fj QVjao ) DVjao
C
                        DSCALE = DF(IOPEN)*
     &                           (DA(IOPEN)-D1)
C
C                       QViao = 2*fi * 1/2 * Fock(TWOINT)
C
                        DSCALE = DSCALE * DF(IOPEN)
C
                        CALL DAXPY(N2BBASXQ,DSCALE,
     &                     FAO(1,1+IOPEN),1,DTEMP1(1,3),1)
                        DSCALE = -DF(JOPEN)*
     &                           (DA(JOPEN)-D1)
C
C                       QVjao = 2*fj * 1/2 * Fock(TWOINT)
C
                        DSCALE = DSCALE * DF(JOPEN)
C
                        CALL DAXPY(N2BBASXQ,DSCALE,
     &                     FAO(1,1+JOPEN),1,DTEMP1(1,3),1)
C                       Here we multiply the combination of the
C                       to open-shell Fock matrices with the right-hand
C                       projector P^V
                        CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),D1,
     &                     'N','N',IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),3),
     &                     NTBAS(0),NTBAS(0),NZ,
     &                     'N','N',IPQTOQ(1,0),
     &                        DMAT(1+I2BASX(I,I),1+JOPEN),
     &                     NTBAS(0),NTBAS(0),NZ,
     &                     D1,IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),2),
     &                     NTBAS(0),NTBAS(0),NZ)
C
                     END IF
C
C
                  ELSE
C                    This is case 1, corresponding to (3.43)
C
C                    ( + (ai-1)/(fi-fj) fi QViao
C                      - (aj-1)/(fi-fj) fj QVjao ) DVjao
C
                     DSCALE = DF(IOPEN)*
     &                        (DA(IOPEN)-D1)/(DF(IOPEN)-DF(JOPEN))
C
C                    QViao = 2*fi * 1/2 * Fock(TWOINT)
C
                     DSCALE = DSCALE * DF(IOPEN)
C
                     CALL DAXPY(N2BBASXQ,DSCALE,
     &                  FAO(1,1+IOPEN),1,DTEMP1(1,3),1)
                     DSCALE = -DF(JOPEN)*
     &                        (DA(JOPEN)-D1)/(DF(IOPEN)-DF(JOPEN))
C
C                    QVjao = 2*fj * 1/2 * Fock(TWOINT)
C
                     DSCALE = DSCALE * DF(JOPEN)
C
                     CALL DAXPY(N2BBASXQ,DSCALE,
     &                  FAO(1,1+JOPEN),1,DTEMP1(1,3),1)
C                    Here we multiply with the right-hand projector P^V
                     CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),D1,
     &                  'N','N',IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),3),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  'N','N',IPQTOQ(1,0),DMAT(1+I2BASX(I,I),1+JOPEN),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  D1,IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),2),
     &                  NTBAS(0),NTBAS(0),NZ)
C
                  END IF
C
 25            CONTINUE
C
C              Multiply the lot with DViao (store in DTEMP2),
C              that is, we multiply from the left woth P^U
C
C
               CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),D1,
     &            'N','N',IPQTOQ(1,0),DMAT(1+I2BASX(I,I),1+IOPEN),
     &            NTBAS(0),NTBAS(0),NZ,
     &            'N','N',IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),2),
     &            NTBAS(0),NTBAS(0),NZ,
     &            D0,IPQTOQ(1,0),DTEMP2(1+I2BASX(I,I)),
     &            NTBAS(0),NTBAS(0),NZ)
C
               IF ( IPRSCF .GE. 15 ) THEN
                  WRITE(CTMP,'(I1)') I
                  CALL HEADER('MKMOFK: coupling operator in AO-basis'//
     &               ', corep '//CTMP,-1)
                  CALL PRQMAT(DTEMP2(1+I2BASX(I,I)),
     &               NFBAS(I,0),NFBAS(I,0),
     &               NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
               END IF
C
C              Transform with inverse T from left and right
C              Store the result in DTEMP1(1,1)
C
               CALL QTRANS('AOMO','S',D0,
     &           NFBAS(I,0),NFBAS(I,0),NTMO(I),NTMO(I),
     &           DTEMP2(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           DTEMP1(1+I2TMOT(I),1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &           TDINV(IOFFTDINV),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &           TDINV(IOFFTDINV),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRSCF)
               IF ( IPRSCF .GE. 15 ) THEN
                  WRITE(CTMP,'(I1)') I
                  CALL HEADER('MKMOFK: Correction matrix in MO-basis'//
     &               ' step 1, corep '//CTMP,-1)
                  CALL PRQMAT(DTEMP1(1+I2TMOT(I),1),NTMO(I),NTMO(I),
     &               NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
               END IF
C
C              Calculate A + At
C
               CALL SUMMMH(DTEMP1(1+I2TMOT(I),1),NTMO(I),NZ,
     &                  NTMO(I),NTMO(I))
               IF ( IPRSCF .GE. 10 ) THEN
                  WRITE(CTMP,'(I1)') I
                  CALL HEADER('MKMOFK: Correction matrix in MO-basis'//
     &               ' step 2, corep '//CTMP,-1)
                  CALL PRQMAT(DTEMP1(1+I2TMOT(I),1),NTMO(I),NTMO(I),
     &               NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
               END IF
C
C              Sum up FD + A
C
               CALL QAXPY(N2TMO(I),'N','N',D1,IPQTOQ(1,0),1,
     &            DTEMP1(1+I2TMOT(I),1),N2TMO(I),1,
     &            IPQTOQ(1,0),NZ,
     &            FMO(1+I2TMOT(I)),N2TMO(I),1,
     &            IPQTOQ(1,0),NZ)

               IF ( IPRSCF .GE. 7 ) THEN
                  CALL HEADER(
     &            'MKMOFK: Total Fock matrix (output) in MO-basis',-1)
                   WRITE(LUPRI,'(A,I1,A,I1)')
     &             '* Fermion ircop no.',I,'/',NFSYM
                  CALL PRQMAT(FMO(1+I2TMOT(I)),NTMO(I),NTMO(I),
     &               NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
               END IF
            END IF
            IOFFTDINV = IOFFTDINV + NZT*NTMO(I)*NFBAS(I,0)
 20      CONTINUE
 10   CONTINUE
C
 30   CONTINUE
      CALL QEXIT('MKMOFK')
      CALL FLSHFO(LUPRI)
      RETURN
      END
! -- end of diropen.F --
