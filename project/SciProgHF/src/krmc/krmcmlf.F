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
C     Subroutines related to the Molfdir format:
C
C     QFC2MFC, MFC2QFC: quaternion active-active matrix to
C                       Molfdir active-active matrix and vice versa.
C     M2DNZ3, DNZ32M  : active Molfdir matrix to Dirac (NZ,3) format and
C                       vice versa.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck QFC2MFC */
      SUBROUTINE QFC2MFC(DQFC,DMFC,ISYM,ITIM,IPRINT)
C***********************************************************************
C
C
C     Input : DQFC  - FCAC matrix in quaternion basis
C                     (FCAC may also be a hermitian or an antihermitian
C                      property matrix)
C             ISYM  - Boson symmetry of FCAC
C             ITIM  - Time reversal symmetry of FCAC
C
C     Output: DMFC  - FCAC matrix in molfdir basis
C
C
C     Written by J. Thyssen - Jun 27 2000
C
C     Last revision : S. Knecht - Nov 2008  
C                     removed nasty bug in offset calculation 
C                     for DCOPY if ITIM == -1 and MZ == 2.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION DQFC(NASHT,NASHT,NZ)
      DIMENSION DMFC(2*NASHT,2*NASHT,2)
C
C
      CALL QENTER('QFC2MFC')
      IF (ITIM.NE.1 .AND. ITIM.NE.-1) THEN
         CALL QUIT('Only ITIM equal to 1 or -1 implemented')
      END IF
      IREP  = ISYM-1
      IOPSY = JBTOF(IREP,1)
CSK   IPRINT = 30
      IF ( IPRINT .GE. 30 ) THEN
         CALL HEADER('QFC2MFC: Input matrix:',-1)
         WRITE(LUPRI,'(A,3I5)')
     &   ' - Boson, Fermion, and Time reversal symmetry :',
     &   ISYM,IOPSY,ITIM
         CALL PRQMAT(DQFC,NASHT,NASHT,NASHT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      IF ( NZ .EQ. 1 ) THEN
C        DMFC matrix is pure real or pure imaginary, in this
C        special case we store in DMFC(*,*,1) whether real or imaginary.
         MZ = 1
         CALL DZERO(DMFC,4 * NASHT**2)
      ELSE
C        In this case real part of DMFC matrix is stored in DMFC(*,*,1)
C        and imaginary part of DMFC matrix is stored in DMFC(*,*,2).
         MZ = 2
         CALL DZERO(DMFC,8 * NASHT**2)
      END IF
      DO IFSY = 1, NFSYM
         JFSY = MOD(IFSY+IOPSY,2) + 1
         IOFFQ = IASH(IFSY)
         IOFFA = 2*IASH(IFSY)
         IOFFB = IOFFA + NASH(IFSY)
         JOFFQ = IASH(JFSY)
         JOFFA = 2*IASH(JFSY)
         JOFFB = JOFFA + NASH(JFSY)
      DO IPQZ = 1,NZ
         IQZ = IPQTOQ(IPQZ,IREP)
CSK      WRITE(LUPRI,*) ' IQZ is',IQZ
      IF (ITIM .EQ. 1) THEN
C     ... ITIM .eq.  1, time invariant property
C        DQFC = + AR + AI qi + BR qj + BI qk
C        DMFC(a,a) =  AR + AI i = (+DQFC(1)) + (+DQFC(2)) i
C        DMFC(a,b) =  BR + BI i = (+DQFC(3)) + (+DQFC(4)) i
C        DMFC(b,a) = -BR + BI i = (-DQFC(3)) + (+DQFC(4)) i
C        DMFC(b,b) =  AR - AI i = (+DQFC(1)) + (-DQFC(2)) i
         IF (IQZ .EQ. 1) THEN
C           WRITE(LUPRI,*) ' copy AR'
C           copy  AR
            DO J = 1, NASH(JFSY)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFA,J+JOFFA,1),1)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFB,J+JOFFB,1),1)
            ENDDO
         END IF
         IF (IQZ .EQ. 2) THEN
C           copy  AI
C           WRITE(LUPRI,*) ' copy AI'
            DO J = 1, NASH(JFSY)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFA,J+JOFFA,MZ),1)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFB,J+JOFFB,MZ),1)
               CALL DSCAL(NASH(IFSY),DM1,DMFC(1+IOFFB,J+JOFFB,MZ),1)
            ENDDO
         END IF
         IF (IQZ .EQ. 3) THEN
C           copy  BR
            DO J = 1, NASH(JFSY)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFA,J+JOFFB,1),1)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFB,J+JOFFA,1),1)
               CALL DSCAL(NASH(IFSY),DM1,DMFC(1+IOFFB,J+JOFFA,1),1)
            ENDDO
         ENDIF
         IF (IQZ .EQ. 4) THEN
C           copy  BI
            DO J = 1, NASH(JFSY)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFA,J+JOFFB,MZ),1)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFB,J+JOFFA,MZ),1)
            ENDDO
         ENDIF
      ELSE
C     ... ITIM .eq. -1, time variant property
C        remember that the quaternion matrix DQFC is for '-i P'
C        where P is the current property operator;
C        symbolically: DMFC = QFC2MFC(i * DQFC)
C        DQFC = i * ( + AR + AI qi + BR qj + BI qk )
C             = - AI + AR qi - BI qj + BR qk
C        DMFC(a,a) =  AR + AI i = (+DQFC(2)) + (-DQFC(1)) i
C        DMFC(a,b) =  BR + BI i = (+DQFC(4)) + (-DQFC(3)) i
C        DMFC(b,a) =  BR - BI i = (+DQFC(4)) + (+DQFC(3)) i
C        DMFC(b,b) = -AR + AI i = (-DQFC(2)) + (-DQFC(1)) i
         IF (IQZ .EQ. 1) THEN
C           copy -AI
C           WRITE(LUPRI,*) ' copy -AI'
            DO J = 1, NASH(JFSY)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFA,J+JOFFA,MZ),1)
Cbug &                               DMFC(1+IOFFA,J+JOFFA,1),MZ)
               CALL DSCAL(NASH(IFSY),DM1,DMFC(1+IOFFA,J+JOFFA,MZ),1)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFB,J+JOFFB,MZ),1)
               CALL DSCAL(NASH(IFSY),DM1,DMFC(1+IOFFB,J+JOFFB,MZ),1)
            ENDDO
         END IF
         IF (IQZ .EQ. 2) THEN
C           copy  AR
C           WRITE(LUPRI,*) ' copy AR'
            DO J = 1, NASH(JFSY)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFA,J+JOFFA,1),1)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFB,J+JOFFB,1),1)
               CALL DSCAL(NASH(IFSY),DM1,DMFC(1+IOFFB,J+JOFFB,1),1)
            ENDDO
         END IF
         IF (IQZ .EQ. 3) THEN
C           copy -BI
            DO J = 1, NASH(JFSY)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFA,J+JOFFB,MZ),1)
               CALL DSCAL(NASH(IFSY),DM1,DMFC(1+IOFFA,J+JOFFB,MZ),1)
               CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                               DMFC(1+IOFFB,J+JOFFA,MZ),1)
            ENDDO
         ENDIF
         IF (IQZ .EQ. 4) THEN
C           copy  BR
            DO J = 1, NASH(JFSY)
                  CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                                  DMFC(1+IOFFA,J+JOFFB,1),1)
                  CALL DCOPY(NASH(IFSY),DQFC(1+IOFFQ,J+JOFFQ,IPQZ),1,
     &                                  DMFC(1+IOFFB,J+JOFFA,1),1)
            ENDDO
         ENDIF
      END IF
C     ... end if ITIM check
      END DO
C     ... end do IZ = 1,NZ
      END DO
C     ... end do IFSY = 1,NFSYM
C
!#define DMRG_DEBUG
#ifdef DMRG_DEBUG

      do IZ = 1, min(nz,2)
      DO L = 1,2*NASHT
         DO K = 1,2*NASHT
           if(k/=l) DMFC(k,l,iz) = 0.0d0
         END DO
      END DO
      end do

#undef DMRG_DEBUG
#endif
      IF ( IPRINT .GE. 30 ) THEN
         CALL HEADER('QFC2MFC: Output matrix:',-1)
         CALL PRQMAT(DMFC,2*NASHT,2*NASHT,2*NASHT,2*NASHT,MIN(NZ,2),
     $        IPQTOQ(1,0),LUPRI)
      END IF
CSK   IPRINT = 00

C
      CALL QEXIT('QFC2MFC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mfc2qfc */
      SUBROUTINE MFC2QFC(KSYMM,DMFC,DQFC,IPRINT)
      use memory_allocator
C***********************************************************************
C
C
C     Input : DMFC  - FCAC matrix in molfdir basis
C             KSYMM - Do Kramers symmetrization
C
C     Output: DQFC  - FCAC matrix in quaternion basis
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION DQFC(NASHT,NASHT,NZ)
      DIMENSION DMFC(2*NASHT,2*NASHT,2)
      LOGICAL KSYMM, ex, fileopen
      integer, allocatable :: mjvec(:)
C
C
csk   IPRINT = 30
      IF ( IPRINT .GE. 30 ) THEN
         CALL HEADER('MFC2QFC: Input matrix:',-1)
         CALL PRQMAT(DMFC,2*NASHT,2*NASHT,2*NASHT,2*NASHT,MIN(NZ,2),
     &        IPQTOQ(1,0),LUPRI)
      END IF
C
C             -----------
C             | A' | B' |
C     DMFC =  |----|----|   DQFC = DQFC^+ = [ A+ | B+ ]
C             | B" | A" |
C             -----------
C
C     where A+ = A' + (A")* and B+ = B' - (B")*
C     (if KSYMM is false then we assume that A' = (A")* and
C     B' = -(B")*).
C
C
      CALL DZERO(DQFC, NASHT**2 * NZ )
      IOFF1 = 0
      IOFF2 = 0
      IOFF3 = 0
      DO IFSYM = 1, NFSYM
         IOFF3 = IOFF3 + NASH(IFSYM)
         DO J = 1, NASH(IFSYM)
            CALL DCOPY(NASH(IFSYM), DMFC(1+IOFF2,J+IOFF2,1), 1,
     &           DQFC (1+IOFF1,J+IOFF1,1),1 )
            IF (KSYMM) CALL DAXPY(NASH(IFSYM),D1,
     &           DMFC(1+IOFF3,J+IOFF3,1),1,
     &           DQFC(1+IOFF1,J+IOFF1,1),1)
         ENDDO
         IF (NZ.GE.2) THEN
            DO J = 1, NASH(IFSYM)
               CALL DCOPY(NASH(IFSYM), DMFC(1+IOFF2,J+IOFF2,2), 1,
     &              DQFC (1+IOFF1,J+IOFF1,2),1 )
               IF (KSYMM) CALL DAXPY(NASH(IFSYM),DM1,
     &              DMFC(1+IOFF3,J+IOFF3,2),1,
     &              DQFC(1+IOFF1,J+IOFF1,2),1)
            ENDDO
         END IF
         IF (NZ.GE.4) THEN
            DO J = 1, NASH(IFSYM)
               CALL DCOPY(NASH(IFSYM), DMFC(1+IOFF2,J+IOFF3,1), 1,
     &              DQFC (1+IOFF1,J+IOFF1,3), 1)
               CALL DCOPY(NASH(IFSYM), DMFC(1+IOFF2,J+IOFF3,2), 1,
     &              DQFC (1+IOFF1,J+IOFF1,4), 1)
               IF (KSYMM) THEN
                  CALL DAXPY(NASH(IFSYM),DM1,
     &                 DMFC(1+IOFF3,J+IOFF2,1),1,
     &                 DQFC(1+IOFF1,J+IOFF1,3),1)
                  CALL DAXPY(NASH(IFSYM),D1,
     &                 DMFC(1+IOFF3,J+IOFF2,2),1,
     &                 DQFC(1+IOFF1,J+IOFF1,4),1)
               END IF
            ENDDO
         ENDIF
         IOFF1 = IOFF1 + NASH(IFSYM)
         IOFF2 = IOFF2 + 2*NASH(IFSYM)
         IOFF3 = IOFF3 + NASH(IFSYM)
      ENDDO

#ifdef LINSYM_DEBUG
!     
!     linear symmetry - make sure only proper <mj|mj> off-diagonal terms
!     enter the KR-MCSCF machinery... (get rid off numerical noise)
      if(linear)then

        call alloc(mjvec,norbt)
        mjvec =  0
        luni  = -1
        inquire(file='KRMCSCF',exist=ex,opened=fileopen,number=luni)
        if(ex.and.fileopen)then
           call ireakrmc(luni,'MJVEC   ',mjvec,norbt)
        else if( ex.and.(.not. fileopen))then
          luni = 99
          call opnfil(luni,'KRMCSCF','OLD    ','MFC2QF')
          call ireakrmc(luni,'MJVEC   ',mjvec,norbt)
          close(luni,status='keep')
        else if(.not. ex)then
          call quit('*** error in mfc2qfc: no mj-vector available
     &              or block-diagonalization. ***')
        end if
!       call iwrtmamn(mjvec,1,norbt,1,norbt,lupri)
        call block_one_dens(mjvec,dqfc)
        call dealloc(mjvec)

      end if
#endif
!
!     print section
!     IPRINT = 30
      IF ( IPRINT .GE. 30 ) THEN
         CALL HEADER('MFC2QFC: Output matrix:',-1)
         CALL PRQMAT(DQFC,NASHT,NASHT,NASHT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
!     IPRINT = 00
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck m2dnz3 */
      SUBROUTINE M2DNZ3(TSYMM,KSYMM,PSYMM,VMUUUU,VDUUUU,IPRINT)
C***********************************************************************
C
C     Transform a UUUU matrix in Molfdir format to a UUUU matrix
C     in Dirac (NZ,3) format.
C
C     Input:
C        TSYMM: transition density matrix symmetrization:
C           do <0|...|B> + <B|...0> symmetrization
C        KSYMM: do Kramers symmetrization (i.e., calculate P++)
C        PSYMM: particle symmetrization: pqrs + rspq
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
C
      DIMENSION VDUUUU(NASHT,NASHT,NNASHX,NZ,*)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      DIMENSION NSTR(2,0:2,4)
      DIMENSION KR(-300:300)
      DIMENSION DINT(8,2)
C
      LOGICAL PSYMM, KSYMM, TSYMM
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
      CALL DZERO(VDUUUU,NASHT * NASHT * NNASHX * NZ * 3)
C
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO
C
 9001 FORMAT(4I5,5X,4F20.10)
C
      IF (IPRINT .GE. 30) THEN
         THR_PRINT = 1.0D-12
         WRITE(LUPRI,'(/A,3(1X,L1)//A)')
     &        ' (M2DNZ3) TSYMM, KSYMM, PSYMM = ',TSYMM,KSYMM,PSYMM,
     &        ' (M2DNZ3) Input VUUUU matrix in Molfdir format'
         DO L = 1,2*NASHT
            DO K = 1,2*NASHT
               DO J = 1,2*NASHT
                  DO I = 1,2*NASHT
                     IF ( (ABS(VMUUUU(I,J,K,L,1)) .GT. THR_PRINT 
     &                .OR. ABS(VMUUUU(I,J,K,L,2)) .GT. THR_PRINT))
     &                    WRITE(LUPRI,9001)
     &                    I,J,K,L,(VMUUUU(I,J,K,L,IZ),IZ=1,NZ_in_CI)
                  END DO
               END DO
            END DO
         END DO
      END IF
C
C
      IJ = 0
      DO II = 1, NASHT
         DO JJ = 1, II
            IJ = IJ + 1
            DO LL = 1, NASHT
               DO KK = 1, NASHT
                  I = KR(II)
                  J = KR(JJ)
                  K = KR(KK)
                  L = KR(LL)
                  ITR = KR(-II)
                  JTR = KR(-JJ)
                  KTR = KR(-KK)
                  LTR = KR(-LL)
C
C                 Calculate integral classes.
C
                  IF (TSYMM) THEN
                     IF (PSYMM) THEN
                        IF (KSYMM) THEN
C
C                          (1) Symm. transition density matrix
C                          (2) Do particle symmetrization
C                          (3) Do Kramer's symmetrization
C
                           VDUUUU(KK,LL,IJ,1,1) =
     &                          VMUUUU(K,L,I,J,1) + VMUUUU(L,K,J,I,1)
     &                          + VMUUUU(I,J,K,L,1) + VMUUUU(J,I,L,K,1)
     &                          + VMUUUU(LTR,KTR,I,J,1)
     &                          + VMUUUU(KTR,LTR,J,I,1)
     &                          + VMUUUU(I,J,LTR,KTR,1)
     &                          + VMUUUU(J,I,KTR,LTR,1)
     &                          + VMUUUU(K,L,JTR,ITR,1)
     &                          + VMUUUU(L,K,ITR,JTR,1)
     &                          + VMUUUU(JTR,ITR,K,L,1)
     &                          + VMUUUU(ITR,JTR,L,K,1)
     &                          + VMUUUU(LTR,KTR,JTR,ITR,1)
     &                          + VMUUUU(KTR,LTR,ITR,JTR,1)
     &                          + VMUUUU(JTR,ITR,LTR,KTR,1)
     &                          + VMUUUU(ITR,JTR,KTR,LTR,1)
                           VDUUUU(KK,LL,IJ,1,2) =
     $                          VMUUUU(K,LTR,I,JTR,1)
     $                          + VMUUUU(I,JTR,K,LTR,1)
     $                          + VMUUUU(LTR,K,JTR,I,1)
     $                          + VMUUUU(JTR,I,LTR,K,1)
     $                          - VMUUUU(L,KTR,I,JTR,1)
     $                          - VMUUUU(I,JTR,L,KTR,1)
     $                          - VMUUUU(KTR,L,JTR,I,1)
     $                          - VMUUUU(JTR,I,KTR,L,1)
     $                          - VMUUUU(K,LTR,J,ITR,1)
     $                          - VMUUUU(J,ITR,K,LTR,1)
     $                          - VMUUUU(LTR,K,ITR,J,1)
     $                          - VMUUUU(ITR,J,LTR,K,1)
     $                          + VMUUUU(L,KTR,J,ITR,1)
     $                          + VMUUUU(J,ITR,L,KTR,1)
     $                          + VMUUUU(KTR,L,ITR,J,1)
     $                          + VMUUUU(ITR,J,KTR,L,1)
                           VDUUUU(KK,LL,IJ,1,3) =
     $                          VMUUUU(KTR,L,I,JTR,1)
     $                          + VMUUUU(I,JTR,KTR,L,1)
     $                          + VMUUUU(L,KTR,JTR,I,1)
     $                          + VMUUUU(JTR,I,L,KTR,1)
     $                          - VMUUUU(LTR,K,I,JTR,1)
     $                          - VMUUUU(I,JTR,LTR,K,1)
     $                          - VMUUUU(K,LTR,JTR,I,1)
     $                          - VMUUUU(JTR,I,K,LTR,1)
     $                          - VMUUUU(KTR,L,J,ITR,1)
     $                          - VMUUUU(J,ITR,KTR,L,1)
     $                          - VMUUUU(L,KTR,ITR,J,1)
     $                          - VMUUUU(ITR,J,L,KTR,1)
     $                          + VMUUUU(LTR,K,J,ITR,1)
     $                          + VMUUUU(J,ITR,LTR,K,1)
     $                          + VMUUUU(K,LTR,ITR,J,1)
     $                          + VMUUUU(ITR,J,K,LTR,1)
                        ELSE ! not KSYMM
C
C                          (1) Symm. transition density matrix
C                          (2) Do particle symmetrization
C
                           VDUUUU(KK,LL,IJ,1,1) =
     &                          VMUUUU(K,L,I,J,1) + VMUUUU(L,K,J,I,1)
     &                          + VMUUUU(I,J,K,L,1) + VMUUUU(J,I,L,K,1)
                           VDUUUU(KK,LL,IJ,1,2) =
     &                          VMUUUU(K,LTR,I,JTR,1)
     &                          + VMUUUU(LTR,K,JTR,I,1)
     &                          + VMUUUU(I,JTR,K,LTR,1)
     &                          + VMUUUU(JTR,I,LTR,K,1)
                           VDUUUU(KK,LL,IJ,1,3) =
     &                          VMUUUU(KTR,L,I,JTR,1)
     &                          + VMUUUU(L,KTR,JTR,I,1)
     &                          + VMUUUU(I,JTR,KTR,L,1)
     &                          + VMUUUU(JTR,I,L,KTR,1)
                        END IF
                     ELSE ! not PSYMM
                        IF (KSYMM) THEN
C
C                          (1) Symm. transition density matrix
C                          (3) Kramers symmetrization
C
                           VDUUUU(KK,LL,IJ,1,1) =
     &                          VMUUUU(K,L,I,J,1)
     $                          + VMUUUU(L,K,J,I,1)
     &                          + VMUUUU(LTR,KTR,I,J,1)
     &                          + VMUUUU(KTR,LTR,J,I,1)
     &                          + VMUUUU(K,L,JTR,ITR,1)
     &                          + VMUUUU(L,K,ITR,JTR,1)
     &                          + VMUUUU(LTR,KTR,JTR,ITR,1)
     &                          + VMUUUU(KTR,LTR,ITR,JTR,1)
                           VDUUUU(KK,LL,IJ,1,2) =
     $                          VMUUUU(K,LTR,I,JTR,1)
     $                          + VMUUUU(LTR,K,JTR,I,1)
     $                          - VMUUUU(L,KTR,I,JTR,1)
     $                          - VMUUUU(KTR,L,JTR,I,1)
     $                          - VMUUUU(K,LTR,J,ITR,1)
     $                          - VMUUUU(LTR,K,ITR,J,1)
     $                          + VMUUUU(L,KTR,J,ITR,1)
     $                          + VMUUUU(KTR,L,ITR,J,1)
                           VDUUUU(KK,LL,IJ,1,3) =
     $                          VMUUUU(KTR,L,I,JTR,1)
     $                          + VMUUUU(L,KTR,JTR,I,1)
     $                          - VMUUUU(LTR,K,I,JTR,1)
     $                          - VMUUUU(K,LTR,JTR,I,1)
     $                          - VMUUUU(KTR,L,J,ITR,1)
     $                          - VMUUUU(L,KTR,ITR,J,1)
     $                          + VMUUUU(LTR,K,J,ITR,1)
     $                          + VMUUUU(K,LTR,ITR,J,1)
                        ELSE ! not KSYMM
C
C                          (1) Symm. transition density matrix
C
                           VDUUUU(KK,LL,IJ,1,1) =
     &                          VMUUUU(K,L,I,J,1) + VMUUUU(L,K,J,I,1)
                           VDUUUU(KK,LL,IJ,1,2) =
     &                          VMUUUU(K,LTR,I,JTR,1)
     &                          + VMUUUU(LTR,K,JTR,I,1)
                           VDUUUU(KK,LL,IJ,1,3) =
     &                          VMUUUU(KTR,L,I,JTR,1)
     &                          + VMUUUU(L,KTR,JTR,I,1)
                        END IF
                     END IF
                  ELSE ! not TSYMM
                     IF (PSYMM) THEN
                        IF (KSYMM) THEN
C
C                          (2) Particle symmetrization
C                          (3) Kramers symmetrization
C
                           VDUUUU(KK,LL,IJ,1,1) =
     $                          VMUUUU(K,L,I,J,1)
     $                          + VMUUUU(I,J,K,L,1)
     $                          + VMUUUU(LTR,KTR,I,J,1)
     $                          + VMUUUU(I,J,LTR,KTR,1)
     $                          + VMUUUU(K,L,JTR,ITR,1)
     $                          + VMUUUU(JTR,ITR,K,L,1)
     $                          + VMUUUU(LTR,KTR,JTR,ITR,1)
     $                          + VMUUUU(JTR,ITR,LTR,KTR,1)
                           VDUUUU(KK,LL,IJ,1,2) =
     $                          VMUUUU(K,LTR,I,JTR,1)
     $                          + VMUUUU(I,JTR,K,LTR,1)
     $                          - VMUUUU(L,KTR,I,JTR,1)
     $                          - VMUUUU(I,JTR,L,KTR,1)
     $                          - VMUUUU(K,LTR,J,ITR,1)
     $                          - VMUUUU(J,ITR,K,LTR,1)
     $                          + VMUUUU(L,KTR,J,ITR,1)
     $                          + VMUUUU(J,ITR,L,KTR,1)
                           VDUUUU(KK,LL,IJ,1,3) =
     $                          VMUUUU(KTR,L,I,JTR,1)
     $                          + VMUUUU(I,JTR,KTR,L,1)
     $                          - VMUUUU(LTR,K,I,JTR,1)
     $                          - VMUUUU(I,JTR,LTR,K,1)
     $                          - VMUUUU(KTR,L,J,ITR,1)
     $                          - VMUUUU(J,ITR,KTR,L,1)
     $                          + VMUUUU(LTR,K,J,ITR,1)
     $                          + VMUUUU(J,ITR,LTR,K,1)
                        ELSE ! not KSYMM
C
C                          (2) Particle symmetrization
C
                           VDUUUU(KK,LL,IJ,1,1) =
     $                          VMUUUU(K,L,I,J,1)
     $                          + VMUUUU(I,J,K,L,1)
                           VDUUUU(KK,LL,IJ,1,2) =
     $                          VMUUUU(K,LTR,I,JTR,1)
     $                          + VMUUUU(I,JTR,K,LTR,1)
                           VDUUUU(KK,LL,IJ,1,3) =
     $                          VMUUUU(KTR,L,I,JTR,1)
     $                          + VMUUUU(I,JTR,KTR,L,1)
                        END IF
                     ELSE ! not PSYMM
                        IF (KSYMM) THEN
C
C                          (3) Kramers symmetrization
C
                           VDUUUU(KK,LL,IJ,1,1) =
     $                          VMUUUU(K,L,I,J,1)
     $                        + VMUUUU(LTR,KTR,I,J,1)
     $                        + VMUUUU(K,L,JTR,ITR,1)
     $                        + VMUUUU(LTR,KTR,JTR,ITR,1)
                           VDUUUU(KK,LL,IJ,1,2) =
     $                          VMUUUU(K,LTR,I,JTR,1)
     $                        - VMUUUU(L,KTR,I,JTR,1)
     $                        - VMUUUU(K,LTR,J,ITR,1)
     $                        + VMUUUU(L,KTR,J,ITR,1)
                           VDUUUU(KK,LL,IJ,1,3) =
     $                          VMUUUU(KTR,L,I,JTR,1)
     $                        - VMUUUU(LTR,K,I,JTR,1)
     $                        - VMUUUU(KTR,L,J,ITR,1)
     $                        + VMUUUU(LTR,K,J,ITR,1)
                        ELSE ! not KSYMM
C
C                          No symmetrization
C
                           VDUUUU(KK,LL,IJ,1,1) = VMUUUU(K,L,I,J,1)
                           VDUUUU(KK,LL,IJ,1,2) = VMUUUU(K,LTR,I,JTR,1)
                           VDUUUU(KK,LL,IJ,1,3) = VMUUUU(KTR,L,I,JTR,1)
                        END IF
                     END IF
                  END IF
C
C                 Complex point groups
C                 --------------------
C
                  IF (.NOT.spinfr .AND. .NOT.levyle) THEN
                  ! all IZ =2,3,4 terms are zero
                  IF (NZ .GE. 2) THEN
                     IF (TSYMM) THEN
                        IF (PSYMM) THEN
                           IF (KSYMM) THEN
C
C                             (1) Symm. transition density matrix
C                             (2) Do particle symmetrization
C                             (3) Do Kramer's symmetrization
C
                              VDUUUU(KK,LL,IJ,2,1) =
     &                             VMUUUU(K,L,I,J,2)
     &                             + VMUUUU(I,J,K,L,2)
     $                             - VMUUUU(L,K,J,I,2)
     $                             - VMUUUU(J,I,L,K,2)
     &                             + VMUUUU(LTR,KTR,I,J,2)
     &                             + VMUUUU(I,J,LTR,KTR,2)
     &                             - VMUUUU(KTR,LTR,J,I,2)
     &                             - VMUUUU(J,I,KTR,LTR,2)
     &                             + VMUUUU(K,L,JTR,ITR,2)
     &                             + VMUUUU(JTR,ITR,K,L,2)
     &                             - VMUUUU(L,K,ITR,JTR,2)
     &                             - VMUUUU(ITR,JTR,L,K,2)
     &                             + VMUUUU(LTR,KTR,JTR,ITR,2)
     &                             + VMUUUU(JTR,ITR,LTR,KTR,2)
     &                             - VMUUUU(KTR,LTR,ITR,JTR,2)
     &                             - VMUUUU(ITR,JTR,KTR,LTR,2)
                              VDUUUU(KK,LL,IJ,2,2) =
     $                             VMUUUU(K,LTR,I,JTR,2)
     $                             + VMUUUU(I,JTR,K,LTR,2)
     $                             - VMUUUU(LTR,K,JTR,I,2)
     $                             - VMUUUU(JTR,I,LTR,K,2)
     $                             - VMUUUU(L,KTR,I,JTR,2)
     $                             - VMUUUU(I,JTR,L,KTR,2)
     $                             + VMUUUU(KTR,L,JTR,I,2)
     $                             + VMUUUU(JTR,I,KTR,L,2)
     $                             - VMUUUU(K,LTR,J,ITR,2)
     $                             - VMUUUU(J,ITR,K,LTR,2)
     $                             + VMUUUU(LTR,K,ITR,J,2)
     $                             + VMUUUU(ITR,J,LTR,K,2)
     $                             + VMUUUU(L,KTR,J,ITR,2)
     $                             + VMUUUU(J,ITR,L,KTR,2)
     $                             - VMUUUU(KTR,L,ITR,J,2)
     $                             - VMUUUU(ITR,J,KTR,L,2)
                              VDUUUU(KK,LL,IJ,2,3) =
     $                             VMUUUU(KTR,L,I,JTR,2)
     $                             + VMUUUU(I,JTR,KTR,L,2)
     $                             - VMUUUU(L,KTR,JTR,I,2)
     $                             - VMUUUU(JTR,I,L,KTR,2)
     $                             - VMUUUU(LTR,K,I,JTR,2)
     $                             - VMUUUU(I,JTR,LTR,K,2)
     $                             + VMUUUU(K,LTR,JTR,I,2)
     $                             + VMUUUU(JTR,I,K,LTR,2)
     $                             - VMUUUU(KTR,L,J,ITR,2)
     $                             - VMUUUU(J,ITR,KTR,L,2)
     $                             + VMUUUU(L,KTR,ITR,J,2)
     $                             + VMUUUU(ITR,J,L,KTR,2)
     $                             + VMUUUU(LTR,K,J,ITR,2)
     $                             + VMUUUU(J,ITR,LTR,K,2)
     $                             - VMUUUU(K,LTR,ITR,J,2)
     $                             - VMUUUU(ITR,J,K,LTR,2)
                           ELSE
C
C                            (1) Symm. transition density matrix
C                            (2) Do particle symmetrization
C
                              VDUUUU(KK,LL,IJ,2,1) =
     &                             VMUUUU(K,L,I,J,2)
     $                             - VMUUUU(L,K,J,I,2)
     &                             + VMUUUU(I,J,K,L,2)
     $                             - VMUUUU(J,I,L,K,2)
                              VDUUUU(KK,LL,IJ,2,2) =
     &                             VMUUUU(K,LTR,I,JTR,2)
     &                             - VMUUUU(LTR,K,JTR,I,2)
     &                             + VMUUUU(I,JTR,K,LTR,2)
     &                             - VMUUUU(JTR,I,LTR,K,2)
                              VDUUUU(KK,LL,IJ,2,3) =
     &                             VMUUUU(KTR,L,I,JTR,2)
     &                             - VMUUUU(L,KTR,JTR,I,2)
     &                             + VMUUUU(I,JTR,KTR,L,2)
     &                             - VMUUUU(JTR,I,L,KTR,2)
                           END IF
                        ELSE
                           IF (KSYMM) THEN
C
C                            (1) Symm. transition density matrix
C                            (2) Kramers symmetrization
C
                              VDUUUU(KK,LL,IJ,2,1) =
     &                             VMUUUU(K,L,I,J,2)
     $                             - VMUUUU(L,K,J,I,2)
     &                             + VMUUUU(LTR,KTR,I,J,2)
     &                             - VMUUUU(KTR,LTR,J,I,2)
     &                             + VMUUUU(K,L,JTR,ITR,2)
     &                             - VMUUUU(L,K,ITR,JTR,2)
     &                             + VMUUUU(LTR,KTR,JTR,ITR,2)
     &                             - VMUUUU(KTR,LTR,ITR,JTR,2)
                              VDUUUU(KK,LL,IJ,2,2) =
     $                             VMUUUU(K,LTR,I,JTR,2)
     $                             - VMUUUU(LTR,K,JTR,I,2)
     $                             - VMUUUU(L,KTR,I,JTR,2)
     $                             + VMUUUU(KTR,L,JTR,I,2)
     $                             - VMUUUU(K,LTR,J,ITR,2)
     $                             + VMUUUU(LTR,K,ITR,J,2)
     $                             + VMUUUU(L,KTR,J,ITR,2)
     $                             - VMUUUU(KTR,L,ITR,J,2)
                              VDUUUU(KK,LL,IJ,2,3) =
     $                             VMUUUU(KTR,L,I,JTR,2)
     $                             - VMUUUU(L,KTR,JTR,I,2)
     $                             - VMUUUU(LTR,K,I,JTR,2)
     $                             + VMUUUU(K,LTR,JTR,I,2)
     $                             - VMUUUU(KTR,L,J,ITR,2)
     $                             + VMUUUU(L,KTR,ITR,J,2)
     $                             + VMUUUU(LTR,K,J,ITR,2)
     $                             - VMUUUU(K,LTR,ITR,J,2)
                           ELSE
C
C                            (1) Symm. transition density matrix
C
                              VDUUUU(KK,LL,IJ,2,1) =
     &                             VMUUUU(K,L,I,J,2)
     $                             - VMUUUU(L,K,J,I,2)
                              VDUUUU(KK,LL,IJ,2,2) =
     &                             VMUUUU(K,LTR,I,JTR,2)
     &                             - VMUUUU(LTR,K,JTR,I,2)
                              VDUUUU(KK,LL,IJ,2,3) =
     &                             VMUUUU(KTR,L,I,JTR,2)
     &                             - VMUUUU(L,KTR,JTR,I,2)
                           END IF
                        END IF
                     ELSE
                        IF (PSYMM) THEN
                           IF (KSYMM) THEN
C
C                            (1) Particle symmetrization
C                            (2) Kramers symmetrization
C
                              VDUUUU(KK,LL,IJ,2,1) =
     $                             VMUUUU(K,L,I,J,2)
     $                             + VMUUUU(I,J,K,L,2)
     $                             + VMUUUU(LTR,KTR,I,J,2)
     $                             + VMUUUU(I,J,LTR,KTR,2)
     $                             + VMUUUU(K,L,JTR,ITR,2)
     $                             + VMUUUU(JTR,ITR,K,L,2)
     $                             + VMUUUU(LTR,KTR,JTR,ITR,2)
     $                             + VMUUUU(JTR,ITR,LTR,KTR,2)
                              VDUUUU(KK,LL,IJ,2,2) =
     $                             VMUUUU(K,LTR,I,JTR,2)
     $                             + VMUUUU(I,JTR,K,LTR,2)
     $                             - VMUUUU(L,KTR,I,JTR,2)
     $                             - VMUUUU(I,JTR,L,KTR,2)
     $                             - VMUUUU(K,LTR,J,ITR,2)
     $                             - VMUUUU(J,ITR,K,LTR,2)
     $                             + VMUUUU(L,KTR,J,ITR,2)
     $                             + VMUUUU(J,ITR,L,KTR,2)
                              VDUUUU(KK,LL,IJ,2,3) =
     $                             VMUUUU(KTR,L,I,JTR,2)
     $                             + VMUUUU(I,JTR,KTR,L,2)
     $                             - VMUUUU(LTR,K,I,JTR,2)
     $                             - VMUUUU(I,JTR,LTR,K,2)
     $                             - VMUUUU(KTR,L,J,ITR,2)
     $                             - VMUUUU(J,ITR,KTR,L,2)
     $                             + VMUUUU(LTR,K,J,ITR,2)
     $                             + VMUUUU(J,ITR,LTR,K,2)
                           ELSE
C
C                            (1) Particle symmetrization
C
                              VDUUUU(KK,LL,IJ,2,1) =
     $                             VMUUUU(K,L,I,J,2)
     $                           + VMUUUU(I,J,K,L,2)
                              VDUUUU(KK,LL,IJ,2,2) =
     $                             VMUUUU(K,LTR,I,JTR,2)
     $                           + VMUUUU(I,JTR,K,LTR,2)
                              VDUUUU(KK,LL,IJ,2,3) =
     $                             VMUUUU(KTR,L,I,JTR,2)
     $                           + VMUUUU(I,JTR,KTR,L,2)
                           END IF
                        ELSE
                           IF (KSYMM) THEN
C
C                            (1) Kramers symmetrization
C
                              VDUUUU(KK,LL,IJ,2,1) =
     $                             VMUUUU(K,L,I,J,2)
     $                           + VMUUUU(LTR,KTR,I,J,2)
     $                           + VMUUUU(K,L,JTR,ITR,2)
     $                           + VMUUUU(LTR,KTR,JTR,ITR,2)
                              VDUUUU(KK,LL,IJ,2,2) =
     $                             VMUUUU(K,LTR,I,JTR,2)
     $                           - VMUUUU(L,KTR,I,JTR,2)
     $                           - VMUUUU(K,LTR,J,ITR,2)
     $                           + VMUUUU(L,KTR,J,ITR,2)
                              VDUUUU(KK,LL,IJ,2,3) =
     $                             VMUUUU(KTR,L,I,JTR,2)
     $                           - VMUUUU(LTR,K,I,JTR,2)
     $                           - VMUUUU(KTR,L,J,ITR,2)
     $                           + VMUUUU(LTR,K,J,ITR,2)
                           ELSE
C
C                             No symmetrization
C
                              VDUUUU(KK,LL,IJ,2,1) =
     $                             VMUUUU(K,L,I,J,2)
                              VDUUUU(KK,LL,IJ,2,2) =
     $                             VMUUUU(K,LTR,I,JTR,2)
                              VDUUUU(KK,LL,IJ,2,3) =
     $                             VMUUUU(KTR,L,I,JTR,2)
                           END IF
                        END IF
                     END IF
                  END IF ! NZ .GE. 2
C
                  IF (NZ .EQ. 4) THEN
                     CALL QUIT('*** ERROR in M2DNZ3 ***'
     $                    //' NZ=4 not implemented - '
     $                    //'we are seeking volunteers!')
                     IF (TSYMM) THEN
                        VDUUUU(KK,LL,IJ,3,1) =
     &                       VMUUUU(K,L,I,JTR,1) + VMUUUU(L,K,JTR,I,1)
                        VDUUUU(KK,LL,IJ,4,1) =
     &                       VMUUUU(K,L,I,JTR,2) - VMUUUU(L,K,JTR,I,2)
                        VDUUUU(KK,LL,IJ,3,2) =
     &                       VMUUUU(K,LTR,I,J,1) + VMUUUU(LTR,K,J,I,1)
                        VDUUUU(KK,LL,IJ,4,2) =
     &                       VMUUUU(K,LTR,I,J,2) - VMUUUU(LTR,K,J,I,2)
                        VDUUUU(KK,LL,IJ,3,3) =
     &                       VMUUUU(KTR,L,I,J,1) + VMUUUU(L,KTR,J,I,1)
                        VDUUUU(KK,LL,IJ,4,3) =
     &                       VMUUUU(KTR,L,I,J,2) - VMUUUU(L,KTR,J,I,2)
                     ELSE
                        VDUUUU(KK,LL,IJ,3,1) = VMUUUU(K,L,I,JTR,1)
                        VDUUUU(KK,LL,IJ,4,1) = VMUUUU(K,L,I,JTR,2)
                        VDUUUU(KK,LL,IJ,3,2) = VMUUUU(K,LTR,I,J,1)
                        VDUUUU(KK,LL,IJ,4,2) = VMUUUU(K,LTR,I,J,2)
                        VDUUUU(KK,LL,IJ,3,3) = VMUUUU(KTR,L,I,J,1)
                        VDUUUU(KK,LL,IJ,4,3) = VMUUUU(KTR,L,I,J,2)
                     ENDIF
                  END IF ! NZ .eq. 4
                  END IF ! (.NOT.spinfr .AND. .NOT.levyle)
                         ! for which all IZ = 2,3,4 terms are zero
               END DO
            END DO
         END DO
      END DO
C
!     remove integrals that should be zero in spinfree calculations but
!     are not due to numerical noise
      if(spinfr.or.levyle)then
        do i = 2,3
          call dzero(VDUUUU(1,1,1,1,i),NASHT*NASHT*NNASHX*NZ)
        end do
      end if
C
      IF (IPRINT .GE. 30) THEN
         CALL HEADER('M2DNZ3: VUUUU matrix in Dirac (NZ,3) format',-1)
         CALL PRDNZ3(VDUUUU,NASHT,NNASHX,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dnz32m */
      SUBROUTINE DNZ32M(VDUUUU,VMUUUU,IPRINT)
C***********************************************************************
C
C     Transform a UUUU matrix in Dirac (NZ,3) format to a UUUU matrix
C     in Molfdir format.
C
C     Written by J. Thyssen - Nov 20 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
      PARAMETER (THROUT = 1.0D-15)
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION VDUUUU(NASHT,NASHT,NNASHX,NZ,*)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      DIMENSION NSTR(2,0:2,4)
      DIMENSION KR(-300:300)
      DIMENSION DINT(8,2)
      real*8, allocatable :: vmuuuu_tmp(:,:,:,:,:)
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
      CALL DZERO(VMUUUU, (2*NASHT)**4 * MIN(NZ,2))
C
      iprint_local = iprint
#ifdef DMRG_DEBUG
!     debug
      iprint_local = 30
#endif
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO
C
      IF (iprint_local .GE. 30) THEN
         WRITE(LUPRI,'(/1X,A)')
     &        '(DNZ32M) Input VUUUU matrix in Dirac (NZ,3) format'
         CALL PRDNZ3(VDUUUU,NASHT,NNASHX,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      DO II = 1, NASHT
         DO JJ = 1, II
            IJ = (II * (II - 1))/2 + JJ
            DO LL = 1, NASHT
               DO KK = 1, NASHT
                  I = KR(II)
                  J = KR(JJ)
                  K = KR(KK)
                  L = KR(LL)
                  ITR = KR(-II)
                  JTR = KR(-JJ)
                  KTR = KR(-KK)
                  LTR = KR(-LL)
C
C                 Canonical integrals (I>=J)
C                 (kl|ij), (kL|iJ), (Kl|iJ)
C
                  VMUUUU(K,L,I,J,1)     = VDUUUU(KK,LL,IJ,1,1)
                  VMUUUU(K,LTR,I,JTR,1) = VDUUUU(KK,LL,IJ,1,2)
                  VMUUUU(KTR,L,I,JTR,1) = VDUUUU(KK,LL,IJ,1,3)
C
C                 The rest:
C
C                 (Kl|Ij) = - (Kl|Ji) = - (lK|iJ)^*
                  VMUUUU(KTR,L,ITR,J,1) = - VDUUUU(LL,KK,IJ,1,2)
C                 (KL|ij) = (lk|ij)
                  VMUUUU(KTR,LTR,I,J,1) = VDUUUU(LL,KK,IJ,1,1)
C                 (KL|IJ) = (lk|ji) = (kl|ij)^*
                  VMUUUU(KTR,LTR,ITR,JTR,1) = VDUUUU(KK,LL,IJ,1,1)
C                 (kl|IJ) = (kl|ji) = (lk|ij)
                  VMUUUU(K,L,ITR,JTR,1) = VDUUUU(LL,KK,IJ,1,1)
C                 (kL|Ij) = - (kL|Ji) = - (Lk|iJ)^*
                  VMUUUU(K,LTR,ITR,J,1) = - VDUUUU(LL,KK,IJ,1,3)
C
                  IF (II .NE. JJ) THEN
C
                     VMUUUU(L,K,J,I,1)     = VDUUUU(KK,LL,IJ,1,1)
                     VMUUUU(LTR,K,JTR,I,1) = VDUUUU(KK,LL,IJ,1,2)
                     VMUUUU(L,KTR,JTR,I,1) = VDUUUU(KK,LL,IJ,1,3)
C
                     VMUUUU(L,KTR,J,ITR,1) = - VDUUUU(LL,KK,IJ,1,2)
                     VMUUUU(LTR,KTR,J,I,1) =   VDUUUU(LL,KK,IJ,1,1)
                     VMUUUU(LTR,KTR,JTR,ITR,1) = VDUUUU(KK,LL,IJ,1,1)
                     VMUUUU(L,K,JTR,ITR,1) =   VDUUUU(LL,KK,IJ,1,1)
                     VMUUUU(LTR,K,J,ITR,1) = - VDUUUU(LL,KK,IJ,1,3)
                  END IF
C
C
                  IF (NZ .GE. 2) THEN
C
C                    Canonical integrals (I>=J)
C                    (kl|ij), (kL|iJ), (Kl|iJ)
C
                     VMUUUU(K,L,I,J,2)     = VDUUUU(KK,LL,IJ,2,1)
                     VMUUUU(K,LTR,I,JTR,2) = VDUUUU(KK,LL,IJ,2,2)
                     VMUUUU(KTR,L,I,JTR,2) = VDUUUU(KK,LL,IJ,2,3)
C
C                    The rest:
C
C                    (Kl|Ij) = - (Kl|Ji) = - (lK|iJ)^*
                     VMUUUU(KTR,L,ITR,J,2) = VDUUUU(LL,KK,IJ,2,2)
C                    (KL|ij) = (lk|ij)
                     VMUUUU(KTR,LTR,I,J,2) = VDUUUU(LL,KK,IJ,2,1)
C                    (KL|IJ) = (lk|ji) = (kl|ij)^*
                     VMUUUU(KTR,LTR,ITR,JTR,2) = - VDUUUU(KK,LL,IJ,2,1)
C                    (kl|IJ) = (kl|ji) = (lk|ij)^*
                     VMUUUU(K,L,ITR,JTR,2) = - VDUUUU(LL,KK,IJ,2,1)
C                    (kL|Ij) = - (kL|Ji) = - (Lk|iJ)^*
                     VMUUUU(K,LTR,ITR,J,2) = VDUUUU(LL,KK,IJ,2,3)
C
                     IF (II .NE. JJ) THEN
C
C                       (lk|ji) = (kl|ij)^*
                        VMUUUU(L,K,J,I,2)     = - VDUUUU(KK,LL,IJ,2,1)
C                       (Lk|Ji) = (kL|iJ)^*
                        VMUUUU(LTR,K,JTR,I,2) = - VDUUUU(KK,LL,IJ,2,2)
C                       (lK|Ji) = (Kl|iJ)^*
                        VMUUUU(L,KTR,JTR,I,2) = - VDUUUU(KK,LL,IJ,2,3)
C
C                       (lK|jI) = (kL|iJ)
                        VMUUUU(L,KTR,J,ITR,2) =   VDUUUU(KK,LL,IJ,2,2)
C                       (LK|ji) = (kl|ji) = (lk|ij)^*
                        VMUUUU(LTR,KTR,J,I,2) = - VDUUUU(LL,KK,IJ,2,1)
C                       (LK|JI) = (kl|ij)
                        VMUUUU(LTR,KTR,JTR,ITR,2) = VDUUUU(KK,LL,IJ,2,1)
C                       (lk|JI) = (lk|ij)
                        VMUUUU(L,K,JTR,ITR,2) = VDUUUU(LL,KK,IJ,2,1)
C                       (Lk|jI) = (Kl|iJ)
                        VMUUUU(LTR,K,J,ITR,2) = VDUUUU(KK,LL,IJ,2,3)
                     END IF
C
                  END IF
C
                  IF (NZ .EQ. 4) THEN
C
C                    (Kl|ij) (canonical)
C
                     VMUUUU(KTR,L,I,J,1) = VDUUUU(KK,LL,IJ,3,3)
                     VMUUUU(KTR,L,I,J,2) = VDUUUU(KK,LL,IJ,4,3)
                     IF (II .NE. JJ) THEN
C                       (lK|ji) = (Kl|ij)^*
                        VMUUUU(L,KTR,J,I,1) =   VDUUUU(KK,LL,IJ,3,3)
                        VMUUUU(L,KTR,J,I,2) = - VDUUUU(KK,LL,IJ,4,3)
                     END IF
C
C                    (kL|ij) (canonical)
C
                     VMUUUU(K,LTR,I,J,1) = VDUUUU(KK,LL,IJ,3,2)
                     VMUUUU(K,LTR,I,J,2) = VDUUUU(KK,LL,IJ,4,2)
                     IF (II .NE. JJ) THEN
C                       (Lk|ji) = (kL|ij)^*
                        VMUUUU(LTR,K,J,I,1) =   VDUUUU(KK,LL,IJ,3,2)
                        VMUUUU(LTR,K,J,I,2) = - VDUUUU(KK,LL,IJ,4,2)
                     END IF
C
C                    (kl|Ij) = - (kl|Ji) = - (lk|iJ)^*
C
                     VMUUUU(K,L,ITR,J,1) = - VDUUUU(LL,KK,IJ,3,1)
                     VMUUUU(K,L,ITR,J,2) =   VDUUUU(LL,KK,IJ,4,1)
                     IF (II .NE. JJ) THEN
C                       (lk|jI) = - (lk|iJ)
                        VMUUUU(L,K,J,ITR,1) = - VDUUUU(LL,KK,IJ,3,1)
                        VMUUUU(L,K,J,ITR,2) = - VDUUUU(LL,KK,IJ,4,1)
                     END IF
C
C                    (kl|iJ) (canonical)
C
                     VMUUUU(K,L,I,JTR,1) = VDUUUU(KK,LL,IJ,3,1)
                     VMUUUU(K,L,I,JTR,2) = VDUUUU(KK,LL,IJ,4,1)
                     IF (II .NE. JJ) THEN
C                       (lk|Ji) = (kl|iJ)^*
                        VMUUUU(L,K,JTR,I,1) =   VDUUUU(KK,LL,IJ,3,1)
                        VMUUUU(L,K,JTR,I,2) = - VDUUUU(KK,LL,IJ,4,1)
                     END IF
C
C
C                    (3 bars)
C
C                    (kL|IJ) = -(Kl|ij)^*
C
                     VMUUUU(K,LTR,ITR,JTR,1) = - VDUUUU(KK,LL,IJ,3,3)
                     VMUUUU(K,LTR,ITR,JTR,2) =   VDUUUU(KK,LL,IJ,4,3)
                     IF (II .NE. JJ) THEN
C                       (Lk|JI) = (kL|IJ)^* = -(Kl|ij)
                        VMUUUU(LTR,K,JTR,ITR,1) = - VDUUUU(KK,LL,IJ,3,3)
                        VMUUUU(LTR,K,JTR,ITR,2) = - VDUUUU(KK,LL,IJ,4,3)
                     END IF
C
C                    (Kl|IJ) = -(kL|ij)^*
C
                     VMUUUU(KTR,L,ITR,JTR,1) = - VDUUUU(KK,LL,IJ,3,2)
                     VMUUUU(KTR,L,ITR,JTR,2) =   VDUUUU(KK,LL,IJ,4,2)
                     IF (II .NE. JJ) THEN
C                       (lK|JI) = (Kl|IJ)^* = -(kL|ij)
                        VMUUUU(L,KTR,JTR,ITR,1) = - VDUUUU(KK,LL,IJ,3,2)
                        VMUUUU(L,KTR,JTR,ITR,2) = - VDUUUU(KK,LL,IJ,4,2)
                     END IF
C
C                    (KL|Ij) = -(kl|iJ)^*
C
                        VMUUUU(KTR,LTR,ITR,J,1) = - VDUUUU(KK,LL,IJ,3,1)
                        VMUUUU(KTR,LTR,ITR,J,2) =   VDUUUU(KK,LL,IJ,4,1)
                     IF (II .NE. JJ) THEN
C                       (KL|jI) = -(kl|Ji)^* = -(lk|iJ)
                        VMUUUU(KTR,LTR,J,ITR,1) = - VDUUUU(LL,KK,IJ,3,1)
                        VMUUUU(KTR,LTR,J,ITR,2) = - VDUUUU(LL,KK,IJ,4,1)
                     END IF

                     IF (KK .GE. LL) THEN
                        KL = (KK*(KK-1))/2 + LL
C
C                    (KL|iJ) = -(kl|Ij)^* = -(Ij|kl)^*
C
                        VMUUUU(KTR,LTR,I,JTR,1) = - VDUUUU(II,JJ,KL,3,3)
                        VMUUUU(KTR,LTR,I,JTR,2) =   VDUUUU(II,JJ,KL,4,3)
                       IF (II .NE. JJ) THEN
C                       (KL|Ji) = -(kl|jI)^* = -(jI|kl)^*
                        VMUUUU(KTR,LTR,JTR,I,1) = - VDUUUU(JJ,II,KL,3,2)
                        VMUUUU(KTR,LTR,JTR,I,2) =   VDUUUU(JJ,II,KL,4,2)
                       END IF
                     ELSE ! KK .lt. LL
                        LK = (LL*(LL-1))/2 + KK
C
C                    (KL|iJ) = -(kl|Ij)^* = -(Ij|kl)^* = -(jI|lk)
C
                        VMUUUU(KTR,LTR,I,JTR,1) = - VDUUUU(JJ,II,LK,3,2)
                        VMUUUU(KTR,LTR,I,JTR,2) = - VDUUUU(JJ,II,LK,4,2)
                       IF (II .NE. JJ) THEN
C                       (KL|Ji) = -(kl|jI)^* = -(jI|kl)^* = -(Ij|lk)
                        VMUUUU(KTR,LTR,JTR,I,1) = - VDUUUU(II,JJ,LK,3,3)
                        VMUUUU(KTR,LTR,JTR,I,2) = - VDUUUU(II,JJ,LK,4,3)
                       END IF
                     END IF
C
                  END IF
               END DO
            END DO
         END DO
      END DO
C

!#define DMRG_DEBUG
#ifdef DMRG_DEBUG
      allocate(vmuuuu_tmp(2*nasht,2*nasht,2*nasht,2*nasht,nz_in_ci))
      vmuuuu_tmp = 0
      do IZ = 1, NZ_in_CI
      DO L = 1,2*NASHT
         DO K = 1,2*NASHT
            if(k/=l)cycle
            DO J = 1,2*NASHT
               if(j/=l) cycle
               DO I = 1,2*NASHT
                 if(i/=l) cycle
                  print *, 'debug ints... ',i,j,k,l,VMUUUU(I,J,K,L,IZ)
                  VMUUUU_tmp(I,J,K,L,IZ) = VMUUUU(I,J,K,L,IZ)
               END DO
            END DO
         END DO
      END DO
      end do
      call dcopy(2*nasht*2*nasht*2*nasht*2*nasht*nz_in_ci,
     &           vmuuuu_tmp,1,vmuuuu,1)
      deallocate(vmuuuu_tmp)
#undef DMRG_DEBUG
#endif

C
!9001 FORMAT(1P,1X,4I3,4X,2(D14.7,1X))
 9001 FORMAT(4I5,5X,4F20.10)
      IF (iprint_local .GE. 30) THEN
         WRITE(LUPRI,'(A,1P,D10.2)')
     &   ' (DNZ32M) Output VUUUU matrix in Molfdir format,'//
     &   ' print threshold =',THROUT
#ifdef old_code
! hjaaj mar 2018: new code below groups the 2-electron integrals,
! making it easier to spot some problems (particle symmetry, hermiticity)
! but not Kramers conjugation problems - that requires 2**4=16 times lot more
! lines for each quadruplet, but with "DO I = 1,NASHT" instead.
         DO L = 1,2*NASHT
            DO K = 1,2*NASHT
               DO J = 1,2*NASHT
                  DO I = 1,2*NASHT
                     IF ( ABS(VMUUUU(I,J,K,L,1)) .GT. THROUT .OR.
     &                    ABS(VMUUUU(I,J,K,L,NZ_in_CI)) .GT. THROUT)
     &                    WRITE(LUPRI,9001)
     &                    I,J,K,L,(VMUUUU(I,J,K,L,IZ),IZ=1,NZ_in_CI)
                  END DO
               END DO
            END DO
         END DO
      ELSE IF (iprint_local .GE. 10) THEN
#endif
         DO I = 1, 2*NASHT
            DO J = 1, I
               DO L = 1, I
                  DO K = 1, L
                  WRITE(LUPRI,'(A)')
     & 'New quadruplet re(ij|kl), re(kl|ij), im(ij|kl), im(kl|ij)'
                  WRITE(LUPRI,9001)
     & I,J,K,L,(VMUUUU(I,J,K,L,IZ),VMUUUU(K,L,I,J,IZ),IZ=1,NZ_in_CI)
                  WRITE(LUPRI,9001)
     & J,I,L,K,(VMUUUU(J,I,L,K,IZ),VMUUUU(L,K,J,I,IZ),IZ=1,NZ_in_CI)
                  WRITE(LUPRI,9001)
     & I,J,L,K,(VMUUUU(I,J,L,K,IZ),VMUUUU(L,K,I,J,IZ),IZ=1,NZ_in_CI)
                  WRITE(LUPRI,9001)
     & J,I,K,L,(VMUUUU(J,I,K,L,IZ),VMUUUU(K,L,J,I,IZ),IZ=1,NZ_in_CI)
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck chkmlf */
      SUBROUTINE CHKNZ3R(VDUUUU)
C***********************************************************************
C
C     Check relations for Dirac (NZ,3) matrix.
C
C     Input:
C       None
C
C     Output:
C       None
C
C     Written by J. Thyssen - Jan 21 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "thrzer.h"
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION VDUUUU(NASHT,NASHT,NNASHX,NZ,3)
C
      IJ = 0
      DO I = 1, NASHT
         DO J = 1, I
            IJ = IJ + 1
            DO K = 1, NASHT
               DO L = 1, NASHT
                  KL = MAX(K,L) * ( MAX(K,L) - 1 ) / 2 + MIN(K,L)
C
C                 (pQ|rS) = - (qQ|rS)
C
                  CALL CHKNZ3R2('K(1,2)',I,J,K,L,
     &                 VDUUUU(K,L,IJ,1,2),
     &                 - VDUUUU(L,K,IJ,1,2))
                  IF (NZ .GE. 2) THEN
                     CALL CHKNZ3R2('K(2,2)',I,J,K,L,
     &                    VDUUUU(K,L,IJ,2,2),
     &                    - VDUUUU(L,K,IJ,2,2))
                  END IF
C
C                 (Pq|rS) = - (Qp|rS)
C
                  CALL CHKNZ3R2('K(1,3)',I,J,K,L,
     &                 VDUUUU(K,L,IJ,1,2),
     &                 - VDUUUU(L,K,IJ,1,2))
                  IF (NZ .GE. 2) THEN
                     CALL CHKNZ3R2('K(2,3)',I,J,K,L,
     &                    VDUUUU(K,L,IJ,2,2),
     &                    - VDUUUU(L,K,IJ,2,2))
                  END IF
C
C                 Particle symmetry
C
                  IF (K .GE. L) THEN
                     CALL CHKNZ3R2('P(1,2)',I,J,K,L,
     &                    VDUUUU(K,L,IJ,1,2),
     &                    VDUUUU(I,J,KL,1,2))
                  END IF
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
C
      SUBROUTINE CHKNZ3R2(A,I,J,K,L,DA,DB)
#include "implicit.h"
#include "priunit.h"
#include "thrzer.h"
#include "consts.h"
C
      CHARACTER A*(*)
C
      DM = MAX(ABS(DA),ABS(DB))
      IF (DM .NE. D0) THEN
         IF ( (ABS(DA-DB))/DM .GT. THRZER) THEN
            WRITE(LUPRI,'(1X,A10,4I4,1P,2D25.15)') A,I,J,K,L,DA,DB
         END IF
      END IF
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck chkmlf */
      SUBROUTINE CHKNZ3(VDUUUU,WORK,LWORK)
C***********************************************************************
C
C     Check DNZ32M and M2DNZ3 transformations.
C
C     Input:
C       None
C
C     Output:
C       None
C
C     Written by J. Thyssen - Jan 21 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "thrzer.h"
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION VDUUUU(*)
      DIMENSION WORK(*)
C
#include "memint.h"
C
      CALL MEMGET('REAL',KVDUUUU,NASHT*NASHT*NNASHX*NZ*3,
     &     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVMUUUU,(2*NASHT)**4 * MIN(2,NZ),
     &     WORK,KFREE,LFREE)
      CALL DNZ32M(VDUUUU,WORK(KVMUUUU),0)
      CALL M2DNZ3(.FALSE.,.FALSE.,.FALSE.,WORK(KVMUUUU),WORK(KVDUUUU),0)
C
      write(LUPRI,*) 'chknz3 pass one...'
      DO I = 1, NASHT*NASHT*NNASHX*NZ*3
         IF (ABS(VDUUUU(I) - WORK(KVDUUUU-1+I)) .GT. THRZER) THEN
            WRITE(LUPRI,'(1X,A,I8,1P,2D20.10)')
     &           'chknz3: deviation in element ',I,
     &           vDuuuu(i),WORK(KVDUUUU-1+I)
         END IF
      END DO
C
      write(LUPRI,*) 'chknz3 pass two...'
      CALL DNZ32M(VDUUUU,WORK(KVMUUUU),0)
      CALL M2DNZ3(.TRUE.,.FALSE.,.FALSE.,WORK(KVMUUUU),WORK(KVDUUUU),0)
      DO I = 1, NASHT*NASHT*NNASHX*NZ*3
         IF (ABS(VDUUUU(I) - DP5 * WORK(KVDUUUU-1+I)) .GT. THRZER) THEN
            WRITE(LUPRI,'(1X,A,I8,1P,2D20.10)')
     &           'chknz3: deviation in element ',I,
     &           vDuuuu(i),WORK(KVDUUUU-1+I)
         END IF
      END DO
C
      call memrel('bla',work,kwork,kwork,kfree,lfree)
C
      RETURN
      END
#ifdef UNDEF
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dnz32ml */
      SUBROUTINE DNZ32ML(VDUUUU,VMUUUU,KR,INDL,INDK,IOFF12,NUMINT,
     &                   NUMINTT,IPRINT)
C***********************************************************************
C
C     Transform a UUUU matrix in Dirac (NZ,3) format to a UUUU matrix
C     in Molfdir format.
C
C     Input:
C        VDUUUU: the integrals in Dirac (NZ,3) format
C
C     Output:
C        INDL, INDK: (KL) indeces
C        IOFF12: offset for (IJ) into INDL and INDK
C        NUMINT: number of integrals for (IJ)
C        NUMINTT: total number of integrals
C        KR: KR(I)  gives unbarred orbital number
C            KR(-I) gives barrred orbital number
C
C     Written by J. Thyssen - Jan 18 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
      PARAMETER (THROUT = 1.0D-10)
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION VDUUUU(NASHT,NASHT,NNASHX,NZ,*)
      DIMENSION VMUUUU( (2*NASHT)**4, *)
      DIMENSION KR(-2 * NASHT: 2 * NASHT)
      DIMENSION NUMINT(2*NASHT,2*NASHT), IOFF12(2*NASHT,2*NASHT)
      DIMENSION INDL(*), INDK(*)
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
C
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO
C
      IF (IPRINT .GE. 30) THEN
         WRITE(LUPRI,'(1X,A)')
     &        '(DNZ32M) Input VUUUU matrix in Dirac (NZ,3) format'
         CALL PRDNZ3(VDUUUU,NASHT,NN2ASHX,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      NTOT = 0
      DO II = 1, NASHT
         DO JJ = 1, II
            IJ = MAX(II,JJ) * ( MAX(II,JJ) - 1 ) / 2 + MIN(II,JJ)
C
            I = KR(II)
            J = KR(JJ)
            ITR = KR(-II)
            JTR = KR(-JJ)
C
C           (kl|ij) integrals
C
            NUMINT(I,J) = N2ASHX
            IOFF12(I,J) = NTOT + 1
            DO LL = 1, NASHT
               L = KR(LL)
               DO KK = 1, NASHT
                  K = KR(KK)
                  NTOT = NTOT + 1
C
                  INDK(NTOT) = K
                  INDL(NTOT) = L
                  VMUUUU(NTOT,1) = VDUUUU(KK,LL,IJ,1,1)
                  IF (NZ .GE. 2) VMUUUU(NTOT,2) = VDUUUU(KK,LL,IJ,2,1)
               END DO
            END DO
C
C           (Kl|ij)
C           (kL|ij)
C           (KL|ij)
C
C
C           (kl|ij) integrals
C
            NUMINT(I,J) = N2ASHX
            IOFF12(I,J) = NTOT + 1
            DO LL = 1, NASHT
               L = KR(LL)
               DO KK = 1, NASHT
                  K = KR(KK)
                  NUMINTT = NUMINTT + 1
C
                  INDK(NUMINTT) = K
                  INDL(NUMINTT) = L
                  VMUUUU(NUMINTT,1) = VDUUUU(KK,LL,IJ,1,1)
                  IF (NZ.GE.2) VMUUUU(NUMINTT,2) = VDUUUU(KK,LL,IJ,2,1)
               END DO
            END DO



                  KTR = KR(-KK)
                  LTR = KR(-LL)

C
C                 Canonical integrals (I>=J)
C                 (kl|ij), (kL|iJ), (Kl|iJ)
C
                  VMUUUU(K,L,I,J,1)     = VDUUUU(KK,LL,IJ,1,1)
                  VMUUUU(K,LTR,I,JTR,1) = VDUUUU(KK,LL,IJ,1,2)
                  VMUUUU(KTR,L,I,JTR,1) = VDUUUU(KK,LL,IJ,1,3)
C
C                 The rest:
C
C                 (Kl|Ij) = - (Kl|Ji) = - (lK|iJ)^*
                  VMUUUU(KTR,L,ITR,J,1) = - VDUUUU(LL,KK,IJ,1,2)
C                 (KL|ij) = (lk|ij)
                  VMUUUU(KTR,LTR,I,J,1) =   VDUUUU(LL,KK,IJ,1,1)
C                 (KL|IJ) = (lk|ji) = (kl|ij)^*
                  VMUUUU(KTR,LTR,ITR,JTR,1) = VDUUUU(KK,LL,IJ,1,1)
C                 (kl|IJ) = (kl|ji) = (lk|ij)
                  VMUUUU(K,L,ITR,JTR,1) =   VDUUUU(LL,KK,IJ,1,1)
C                 (kL|Ij) = - (kL|Ji) = - (Lk|iJ)^*
                  VMUUUU(K,LTR,ITR,J,1) = - VDUUUU(LL,KK,IJ,1,3)
C
                  IF (II .NE. JJ) THEN
C
                     VMUUUU(L,K,J,I,1)     = VDUUUU(KK,LL,IJ,1,1)
                     VMUUUU(LTR,K,JTR,I,1) = VDUUUU(KK,LL,IJ,1,2)
                     VMUUUU(L,KTR,JTR,I,1) = VDUUUU(KK,LL,IJ,1,3)
C
                     VMUUUU(L,KTR,J,ITR,1) = - VDUUUU(LL,KK,IJ,1,2)
                     VMUUUU(LTR,KTR,J,I,1) =   VDUUUU(LL,KK,IJ,1,1)
                     VMUUUU(LTR,KTR,JTR,ITR,1) = VDUUUU(KK,LL,IJ,1,1)
                     VMUUUU(L,K,JTR,ITR,1) =   VDUUUU(LL,KK,IJ,1,1)
                     VMUUUU(LTR,K,J,ITR,1) = - VDUUUU(LL,KK,IJ,1,3)
                  END IF
C
C
                  IF (NZ .GE. 2) THEN
C
C                    Canonical integrals (I>=J)
C                    (kl|ij), (kL|iJ), (Kl|iJ)
C
                     VMUUUU(K,L,I,J,2)     = VDUUUU(KK,LL,IJ,2,1)
                     VMUUUU(K,LTR,I,JTR,2) = VDUUUU(KK,LL,IJ,2,2)
                     VMUUUU(KTR,L,I,JTR,2) = VDUUUU(KK,LL,IJ,2,3)
C
C                    The rest:
C
C                    (Kl|Ij) = - (Kl|Ji) = - (lK|iJ)^*
                     VMUUUU(KTR,L,ITR,J,2) = VDUUUU(LL,KK,IJ,2,2)
C                    (KL|ij) = (lk|ij)
                     VMUUUU(KTR,LTR,I,J,2) = VDUUUU(LL,KK,IJ,2,1)
C                    (KL|IJ) = (lk|ji) = (kl|ij)^*
                     VMUUUU(KTR,LTR,ITR,JTR,2) = - VDUUUU(KK,LL,IJ,2,1)
C                    (kl|IJ) = (kl|ji) = (lk|ij)^*
                     VMUUUU(K,L,ITR,JTR,2) = - VDUUUU(LL,KK,IJ,2,1)
C                    (kL|Ij) = - (kL|Ji) = - (Lk|iJ)^*
                     VMUUUU(K,LTR,ITR,J,2) = VDUUUU(LL,KK,IJ,2,3)
C
                     IF (II .NE. JJ) THEN
C
C                       (lk|ji) = (kl|ij)^*
                        VMUUUU(L,K,J,I,2) = - VDUUUU(KK,LL,IJ,2,1)
C                       (Lk|Ji) = (kL|iJ)^*
                        VMUUUU(LTR,K,JTR,I,2) = - VDUUUU(KK,LL,IJ,2,2)
C                       (lK|Ji) = (Kl|iJ)^*
                        VMUUUU(L,KTR,JTR,I,2) = - VDUUUU(KK,LL,IJ,2,3)
C
C                       (lK|jI) = (kL|iJ)
                        VMUUUU(L,KTR,J,ITR,2) = VDUUUU(KK,LL,IJ,2,2)
C                       (LK|ji) = (kl|ji) = (lk|ij)^*
                        VMUUUU(LTR,KTR,J,I,2) = - VDUUUU(LL,KK,IJ,2,1)
C                       (LK|JI) = (kl|ij)
                        VMUUUU(LTR,KTR,JTR,ITR,2) = VDUUUU(KK,LL,IJ,2,1)
C                       (lk|JI) = (lk|ij)
                        VMUUUU(L,K,JTR,ITR,2) = VDUUUU(LL,KK,IJ,2,1)
C                       (Lk|jI) = (Kl|iJ)
                        VMUUUU(LTR,K,J,ITR,2) = VDUUUU(KK,LL,IJ,2,3)
                     END IF
C
                  END IF
C
                  IF (NZ .EQ. 4) THEN
C
C                    (Kl|ij) (canonical)
C
                     VMUUUU(KTR,L,I,J,1) = VDUUUU(KK,LL,IJ,3,3)
                     VMUUUU(KTR,L,I,J,2) = VDUUUU(KK,LL,IJ,4,3)
                     IF (II .NE. JJ) THEN
C                       (lK|ji) = (Kl|ij)^*
                        VMUUUU(L,KTR,J,I,1) = VDUUUU(KK,LL,IJ,3,3)
                        VMUUUU(L,KTR,J,I,2) = - VDUUUU(KK,LL,IJ,4,3)
                     END IF
C
C                    (kL|ij) (canonical)
C
                     VMUUUU(K,LTR,I,J,1) = VDUUUU(KK,LL,IJ,3,2)
                     VMUUUU(K,LTR,I,J,2) = VDUUUU(KK,LL,IJ,4,2)
                     IF (II .NE. JJ) THEN
C                       (Lk|ji) = (kL|ij)^*
                        VMUUUU(L,KTR,J,I,1) = VDUUUU(KK,LL,IJ,3,2)
                        VMUUUU(L,KTR,J,I,2) = - VDUUUU(KK,LL,IJ,4,2)
                     END IF
C
C                    (kl|Ij) = - (kl|Ji) = - (lk|iJ)^*
C
                     VMUUUU(K,L,ITR,J,1) = - VDUUUU(LL,KK,IJ,3,2)
                     VMUUUU(K,L,ITR,J,2) = VDUUUU(LL,KK,IJ,4,2)
                     IF (II .NE. JJ) THEN
C                       (lk|jI) = - (lk|iJ)
                        VMUUUU(L,KTR,J,I,1) = - VDUUUU(KK,LL,IJ,3,2)
                        VMUUUU(L,KTR,J,I,2) = - VDUUUU(KK,LL,IJ,4,2)
                     END IF
C
C                    (kl|iJ) (canonical)
C
                     VMUUUU(K,L,I,JTR,1) = VDUUUU(KK,LL,IJ,3,1)
                     VMUUUU(K,L,I,JTR,2) = VDUUUU(KK,LL,IJ,4,1)
                     IF (II .NE. JJ) THEN
C                       (lk|Ji) = (kl|iJ)^*
                        VMUUUU(L,K,JTR,I,1) = VDUUUU(KK,LL,IJ,3,1)
                        VMUUUU(L,K,JTR,I,2) = - VDUUUU(KK,LL,IJ,4,1)
                     END IF
C
C
C                    (3 bars)
C
C                    (kL|IJ) = (kL|ji) = (Lk|ij)^*
C
                     VMUUUU(K,LTR,ITR,JTR,1) = VDUUUU(LL,KK,IJ,3,3)
                     VMUUUU(K,LTR,ITR,JTR,2) = - VDUUUU(LL,KK,IJ,4,3)
                     IF (II .NE. JJ) THEN
C                       (Lk|JI) = (Lk|ij)
                        VMUUUU(LTR,K,JTR,ITR,1) = VDUUUU(KK,LL,IJ,3,3)
                        VMUUUU(LTR,K,JTR,ITR,2) = VDUUUU(KK,LL,IJ,4,3)
                     END IF
C
C                    (Kl|IJ) = (Kl|ji) = (lK|ij)^*
C
                     VMUUUU(KTR,L,ITR,JTR,1) = VDUUUU(LL,KK,IJ,3,2)
                     VMUUUU(KTR,L,ITR,JTR,2) = - VDUUUU(LL,KK,IJ,4,2)
                     IF (II .NE. JJ) THEN
C                       (lK|JI) = (lK|ij)
                        VMUUUU(L,KTR,J,I,1) = VDUUUU(KK,LL,IJ,3,2)
                        VMUUUU(L,KTR,J,I,2) = VDUUUU(KK,LL,IJ,4,2)
                     END IF
C
C                    (KL|iJ) = (lk|iJ)
C
                     VMUUUU(K,L,ITR,J,1) = VDUUUU(KK,LL,IJ,3,1)
                     VMUUUU(K,L,ITR,J,2) = VDUUUU(KK,LL,IJ,4,1)
                     IF (II .NE. JJ) THEN
C                       (LK|Ji) = (kl|Ji) = (lk|iJ)^*
                        VMUUUU(L,KTR,J,I,1) = VDUUUU(LL,KK,IJ,3,1)
                        VMUUUU(L,KTR,J,I,2) = - VDUUUU(LL,KK,IJ,4,1)
                     END IF
C
C                    (KL|Ij) = - (lk|Ji) = - (kl|iJ)^*
C
                     VMUUUU(K,L,ITR,J,1) = - VDUUUU(KK,LL,IJ,3,1)
                     VMUUUU(K,L,ITR,J,2) = VDUUUU(KK,LL,IJ,4,1)
                     IF (II .NE. JJ) THEN
C                       (LK|jI) = - (kl|iJ)
                        VMUUUU(L,K,JTR,I,1) = - VDUUUU(KK,LL,IJ,3,1)
                        VMUUUU(L,K,JTR,I,2) = - VDUUUU(KK,LL,IJ,4,1)
                     END IF
C
                  END IF
               END DO
            END DO
         END DO
      END DO
C
C
C
!9001 FORMAT(1P,1X,4I3,4X,2(D14.7,1X))
 9001 FORMAT(4I5,5X,2F20.10)
C
      IF (IPRINT .GE. 30) THEN
         WRITE(LUPRI,'(A,1P.D10.2)')
     &   ' (DNZ32ML) Output VUUUU matrix in Molfdir format,'//
     &   ' print threshold =',THROUT
         DO L = 1,2*NASHT
            DO K = 1,2*NASHT
               DO J = 1,2*NASHT
                  DO I = 1,2*NASHT
                     IF ( ABS(VMUUUU(I,J,K,L,1)) .GT. THROUT .OR.
     &                    ABS(VMUUUU(I,J,K,L,2)) .GT. THROUT)
     &                    WRITE(LUPRI,9001)
     &                    I,J,K,L,(VMUUUU(I,J,K,L,IZ),IZ=1,MIN(2,NZ))
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      RETURN
      END
#endif
!***********************************************************************

#ifdef MOD_UNRELEASED
      subroutine write_fcidump_file(vmuuuu,hcore,ecore)
!***********************************************************************
!
!     
!     
!     purpose: write one- and two-electron integrals for a given active
!              space to disk
!
!     written by S. Knecht
!
!     based on an analogous routine written for Dalton by
!     Katharina Boguslawski and Pawel Tecmer
!
!***********************************************************************
      use mospinor_info
      use symmetry_setup_krci
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dcbham.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbkrci.h"
C
      DIMENSION  VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
      DIMENSION  hcore(2*NASHT,2*NASHT,*)
      real(8), intent(in) :: ecore
      real(8), parameter  :: THR_VMUUUU = 1.d-15
      DIMENSION  KR(-300:300)
      integer :: LUFCID = -1
      character(len=30) :: form1
      character(len=30) :: form2
      character(len=30) :: form3
      integer, allocatable  :: orb_sym_vec(:)
      integer, allocatable  :: ibosym(:)
      logical :: ex_mcdf, ex_mcol, isopen

!     set up symmetry information
      IF (GROUP.eq.'C1 ') THEN
         DOUGRP      = 8
         NIRR_DG     = 2
      ELSE IF (GROUP.eq.'Ci ') THEN
         DOUGRP      = 7
         NIRR_DG     = 4
      ELSE IF (GROUP.eq.'Cs ') THEN
         DOUGRP      = 6
         NIRR_DG     = 4
      ELSE IF (GROUP.eq.'C2 ') THEN
         DOUGRP      = 5
         NIRR_DG     = 4
      ELSE IF (GROUP.eq.'D2 ') THEN
         DOUGRP      = 5
         NIRR_DG     = 4
      ELSE IF (GROUP.eq.'C2h') THEN
         DOUGRP      = 4
         NIRR_DG     = 8
      ELSE IF (GROUP.eq.'C2v') THEN
         if(linear)then
!          Cinfv (known as C2v + linear) treated as C32 inside KR-CI
           DOUGRP      = 10
           NIRR_DG     = 128
         else
!          C2v treated as C2 inside KR-CI
           DOUGRP      = 5
           NIRR_DG     = 4
         end if
      ELSE IF (GROUP.eq.'D2h') THEN
         if(linear)then
!          Dinfh (known as D2h + linear) treated as C16h inside KR-CI
           DOUGRP      = 11
           NIRR_DG     = 128
         else if(spinfr.and..not.linear)then
!          D2h spinfree
           DOUGRP      = 9
           NIRR_DG     = 8
           print *,
     & '*** error in fcidump interface: no spinfree formalism yet. ***'
           call quit(
     & '*** error in fcidump interface: no spinfree formalism yet. ***')
         else
!          D2h treated as C2h inside KR-CI
           DOUGRP      = 4
           NIRR_DG     = 8
         end if
      ELSE
         call quit(
     & '*** error in fcidump interface: unknown double group 
     &      specified! ***')
      END IF

!     setup double group multiplication table
      call symmetry_setup_init(dougrp,nirr_dg)

!     initialize mo-spinor information arrays
      call mospinor_info_init()

      allocate(orb_sym_vec(norbt))
!     distribute orbital symmetry information (if applicable)
      if(linear.or.spinfr)then

        orb_sym_vec = 0

!       read data from master file krmcscf or dfcoef
        inquire(file="KRMCSCF",opened=isopen,exist=ex_mcdf)
        if(ex_mcdf)then
          if(.not.isopen) open(lukrmc,status="old",
     &                         form="unformatted")
          rewind(lukrmc)
          if(linear)then
            call ireakrmc(lukrmc,'MJVEC   ',orb_sym_vec,norbt)
          else if(spinfr)then
            call ireakrmc(lukrmc,'IBEIG   ',orb_sym_vec,norbt)
          end if
          if(.not.isopen) close(lukrmc,status="keep")
        else
          inquire(file="KRMCOLD",opened=isopen,exist=ex_mcol,
     &            number=luxxxx)
          if(ex_mcol)then
            if(.not.isopen) then 
              luxxxx = 99
              open(luxxxx,status="old",form="unformatted")
            end if
            rewind(luxxxx)
            if(linear)then
              call ireakrmc(luxxxx,'MJVEC   ',orb_sym_vec,norbt)
            else if(spinfr)then
              call ireakrmc(luxxxx,'IBEIG   ',orb_sym_vec,norbt)
            end if
            if(.not.isopen) close(luxxxx,status="keep")
          else
            call opnfil(lucoef,'DFCOEF','OLD','mjlist ')
            call reacmo(lucoef,'DFCOEF',dummy,dummy,
     &                          orb_sym_vec,dummy,8)
            close(lucoef,status='KEEP')
          end if
        end if
      end if
      
!     place spinors (as MOs) in appropriate boson irreps  for further processing (see next step)
!     generalized routine for handling of all double groups including
!     linear symmetry (dougrp == 10 or 11)
      allocate(ibosym(norbt))
      call match_spinor_2_boson_irrep(ngsob,ngsob2,ibosym,ngsh,
     &                                orb_sym_vec,iorb,npsh,nish,
     &                                nirr_dg,1,nfsym,1,
     &                                dougrp,nirr_dg,mj2rep,
     &                                imosp_dirac_counter1,
     &                                imosp_dirac_counter2,
     &                                imosp_dirac_mjub,
     &                                imosp_dirac_mjb)
      deallocate(ibosym)

      deallocate(orb_sym_vec)

!     define printing format
      form1="(e24.14, 4X, I6, I6, I6, I6)"
      form2="(A11,I3,A7,I2,A5,I2,A1)"
      form3="(e24.14, 4X, I6, I6, I6, I6)"

!     FCIDUMP info
!     write(lupri,*) 
!     write(lupri,*)" FCIDUMP file details "
!     write(lupri,*)"======================" 
!     write(lupri,*) 

!     Print header of FCIDUMP file using MOLPRO format
      if (LUFCID .le. 0) then
         CALL GPOPEN(LUFCID,'FCIDUMP','new',' ',
     &                'FORMATTED',IDUMMY,.false.)
      end if

      write(LUFCID,form2) ' &FCI NORB=', nasht*2 , 
     & ',NELEC=', NKRCIAELEC, ',MS2=',
     &  0, ','
      write(LUFCID,"(A)",advance='no') '  ORBSYM=' 
      do IFSYM=1,NFSYM,1
         if(NASH(IFSYM) /= 0)then
           do i=1,nash(ifsym),1
             write(LUFCID,"(I3,A1)",advance='no') 
     &       imosp_dirac_counter1(i),','
           end do
           do i=nash(ifsym)+1,2*nash(ifsym),1
             write(LUFCID,"(I3,A1)",advance='no') 
     &       imosp_dirac_counter2(i-nash(ifsym)),','
           end do
         end if
      end do

      write(LUFCID,*) 
      write(LUFCID,"(A7,I1)") '  ISYM=',1
      write(LUFCID,"(A5)") ' &END'

      NCLASS = NZ * NZ * NBSYM / NFSYM
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO

      !> 2e-integrals
      DO L = 1,2*NASHT
         DO K = 1,2*NASHT
            DO J = 1,2*NASHT
               DO I = 1,2*NASHT

                  if (abs(VMUUUU(I,J,K,L,1))<=THR_VMUUUU) CYCLE
                  if (nz > 1 .and. abs(VMUUUU(I,J,K,L,2))<= THR_VMUUUU)
     &            CYCLE

!     only real groups needed for now ==> only even number of bars: 0, 2, 4

!     linear symmetry case: To determine the MJ value we need to know the classes :
!     Class 1 : (i    j    | k    l   )
!     Class 2 : (i    j    | kbar lbar)
!     Class 3 : (i    jbar | k    lbar)
!     Class 4 : (ibar j    | k    lbar)

      if(l<=nasht .and. k<=nasht .and. j<=nasht .and. i<=nasht)then !(ij|kl)
        continue
         isym = imosp_dirac_mjub(i)
         jsym = imosp_dirac_mjub(j)
         ksym = imosp_dirac_mjub(k)
         lsym = imosp_dirac_mjub(l)
         iclass = 1
      else if(l> nasht .and. k> nasht .and. j> nasht .and. i> nasht)then!(IJ|KL)
        continue
         isym = imosp_dirac_mjub(i-nasht)
         jsym = imosp_dirac_mjub(j-nasht)
         ksym = imosp_dirac_mjub(k-nasht)
         lsym = imosp_dirac_mjub(l-nasht)
         iclass = 1
      else if(l<=nasht .and. k<=nasht .and. j> nasht .and. i> nasht)then!(IJ|kl)
         isym = imosp_dirac_mjub(i-nasht)
         jsym = imosp_dirac_mjub(j-nasht)
         ksym = imosp_dirac_mjub(k)
         lsym = imosp_dirac_mjub(l)
         iclass = 2
        continue
      else if(l> nasht .and. k> nasht .and. j<=nasht .and. i<=nasht)then!(ij|KL)
        continue
         isym = imosp_dirac_mjub(i)
         jsym = imosp_dirac_mjub(j)
         ksym = imosp_dirac_mjub(k-nasht)
         lsym = imosp_dirac_mjub(l-nasht)
         iclass = 2
      else if(l> nasht .and. k<=nasht .and. j> nasht .and. i<=nasht)then!(iJ|kL)
        continue
         isym = imosp_dirac_mjub(i)
         jsym = imosp_dirac_mjub(j-nasht)
         ksym = imosp_dirac_mjub(k)
         lsym = imosp_dirac_mjub(l-nasht)
         iclass = 3
      else if(l> nasht .and. k<=nasht .and. j<=nasht .and. i> nasht)then!(Ij|kL)
        continue
         isym = imosp_dirac_mjub(i-nasht)
         jsym = imosp_dirac_mjub(j)
         ksym = imosp_dirac_mjub(k)
         lsym = imosp_dirac_mjub(l-nasht)
         iclass = 4
      else if(l<=nasht .and. k> nasht .and. j<=nasht .and. i> nasht)then!(Ij|Kl)
        continue
         isym = imosp_dirac_mjub(i-nasht)
         jsym = imosp_dirac_mjub(j)
         ksym = imosp_dirac_mjub(k-nasht)
         lsym = imosp_dirac_mjub(l)
         iclass = 3
      else if(l<=nasht .and. k> nasht .and. j> nasht .and. i<=nasht)then!(iJ|Kl)
        continue
         isym = imosp_dirac_mjub(i)
         jsym = imosp_dirac_mjub(j-nasht)
         ksym = imosp_dirac_mjub(k-nasht)
         lsym = imosp_dirac_mjub(l)
         iclass = 4
      else
        cycle ! odd number of bars
      end if
                 IF (ICLASS.LE.2) THEN
                    MIJ = - isym + jsym
                 ELSEIF (ICLASS.EQ.3) THEN
                    MIJ = - isym - jsym
                 ELSE
                    MIJ =  isym + jsym 
                 ENDIF
C
                 IF (ICLASS.EQ.1) THEN
                    MKL = - ksym + lsym
                 ELSEIF (ICLASS.EQ.2) THEN
                    MKL =   ksym - lsym
                 ELSE
                    MKL = - ksym - lsym
                 ENDIF

!     write(6,'(2x,a,2i6)') '(linzero): MIJ, MKL are:', MIJ,MKL
                 IF(linear .and. (MIJ+MKL /= 0)) cycle

!     IF (MIJ+MKL /= 0)
!    &write(6,'(2x,a,4i3,a,1p,d10.2,a,i4)')
!    &'(linzero): removed integral(',i,j,k,l,')',VMUUUU(I,J,K,L,1),
!    & ' class: ',iclass
!     IF (MIJ+MKL ==0 .and. dabs(VMUUUU(I,J,K,L,1)) > 1.0d-99)
!    &write(6,'(2x,a,4i3,a,1p,d10.2,a,i4)')
!    &'(linzero): kept integral(',i,j,k,l,')',VMUUUU(I,J,K,L,1), 
!    & ' class: ',iclass

                 write(LUFCID,form1) VMUUUU(I,J,K,L,1),I,J,K,L
               END DO
            END DO
         END DO
      END DO
      !> 1e-integrals
      DO K = 1,NASHT
         DO L = 1,K
           if(linear .and. 
     &        (imosp_dirac_mjub(l) /= imosp_dirac_mjub(k))) cycle
           write(LUFCID,form1) hcore(K,L,1),K,L,0,0
!          write(LUFCID,form1) hcore(K+nasht,L+nasht,1),
!    &                                      K+nasht,L+nasht,0,0
         END DO
      END DO

      !> core energy
      write(LUFCID,form3) ecore, 0,0 ,0,0
 
      call GPCLOSE(LUFCID,'KEEP')

      !> finalize
      call mospinor_info_delete()
      call symmetry_setup_delete()
C
      END
#endif
