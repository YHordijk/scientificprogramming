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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRAPRP */
      SUBROUTINE TRAPRP(WORK,KFREE,LFREE)
C***********************************************************************
C
C     Driver routine for property integral transformations
C
C     Written by L. Visscher Jul 1997
C     Revised Feb 17 2000 hjaaj (new parameter list to PRPMAT)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "dcbham.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "mxcent.h"
#include "dcbtrp.h"
      LOGICAL TOBE,TOBEK,TRASAM
      DIMENSION WORK(*),NSTR(2,0:2,2),NSTRT(2),KVEC(2,2),KQ(2,2),
     &          IR(2,2),NR(2),NQ(2,2),NQT(2),
     &          NSPC(2,0:2),KQC(2),KE(2,2),KIBE(2,2)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2),ICMOQC(2)
      CHARACTER*32 ACHAR
C
      CALL QENTER('TRAPRP')
      KFRSAV = KFREE
      CALL TRAHI(IPRTRP,1)
C
C     =======================================      
C     Set up index arrays for active orbitals
C     =======================================      
C
      IF(NOPAIR) THEN
        IC   = 1
        DO I = 1,NFSYM
          NR(I)   = NESH(I)
          IR(1,I) = 1
          IR(2,I) = NESH(I)
        ENDDO
      ELSE
        IC   = 0
        DO I = 1,NFSYM
          NR(I)   =  NORB(I)
          IR(1,I) = -NPSH(I)
          IR(2,I) =  NESH(I)
        ENDDO
      ENDIF
!
!.s/sya,2007.0907,add/ - SK - 30-11-2009
      IF( PRPSYA ) THEN
        WRITE(LUPRI,"(
     &    '* Setting MO range in property integral transformation')")
        DO IFRP = 1, NFSYM
          TRA2_INDSTR(1,IFRP) = TRA2_INDPRP(IFRP)
          TRA2_INDSTR(2,IFRP) = TRA2_INDPRP(IFRP)
        END DO
      END IF
!.q
C
      DO I = 1,2
        NSTRT(I) = 0
        DO IFRP = 1,NFSYM
          CALL MEMGET('INTE',KVEC(IFRP,I),NR(IFRP),WORK,KFREE,LFREE)
          NSTR(IFRP,0,I) = 1
          CALL  NUMLST(TRA2_INDSTR(I,IFRP),WORK(KVEC(IFRP,I)),
     &                 NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &                 IFRP,NSTR(IFRP,0,I))
          CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I),
     &              NPSH(IFRP),NESH(IFRP),
     &              NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NSTRT(I) = NSTRT(I) +  NSTR(IFRP,0,I)
          CALL MEMREL('TRAPRP',WORK,KFRSAV,KVEC(IFRP,I),KFREE,LFREE)
          CALL MEMGET('INTE',KVEC(IFRP,I),NSTR(IFRP,0,I),
     &                WORK,KFREE,LFREE)
        ENDDO
      ENDDO      
C
C     Print section
C     =============
C
      DO I = 1, NFSYM
         CALL TRAPRI(2,I,WORK(KVEC(I,1)),WORK(KVEC(I,2)),0,0,NSTR)
      ENDDO
C
C     Get all coefficients and eigenvalues
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
      INQUIRE(FILE='KRMCSCF',EXIST=TOBEK)
      IF (TOBEK) THEN 
         CALL OPNFIL(LUKRMC,'KRMCSCF','OLD','TRAPRP')
         JRDMO = -1
         CALL RREADMO(WORK(KCMO),JRDMO,1,LUKRMC)
         IF (JRDMO .EQ. 0) THEN
            IF (IPRTRP .GE. 0) THEN
               WRITE(LUPRI,'(/A)')
     &              ' (TRAPRP)  Orbitals read from' //
     &              ' label NEWORB on file KRMCSCF'
            END IF
         END IF
         CLOSE (LUKRMC)
      ELSE
         INQUIRE(FILE='DFCOEF',EXIST=TOBE)
         IF (.NOT.TOBE) THEN
            WRITE(LUPRI,'(A)') 'Coefficient file not present'
            CALL QUIT('TRAPAM: No coefficients !')
         ELSE
            CALL OPNFIL(LUCOEF,'DFCOEF','OLD','TRAPRP')
            CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),
     &                  WORK(KIBEIG),TOTERG,14)
         ENDIF
      END IF
C
C     Get the boson symmetries of the transformed spinors
C
      DO I = 1,2
       DO IFRP = 1, NFSYM
        CALL MEMGET('INTE',KIBE(IFRP,I),NSTR(IFRP,0,I),WORK,KFREE,LFREE)
        CALL SELIBEIG(WORK(KIBEIG),IORB(IFRP),IFRP,
     &                WORK(KIBE(IFRP,I)),WORK(KVEC(IFRP,I)),
     &                NSTR(IFRP,2,I),NSTR(IFRP,1,I))
       ENDDO
       IF (NFSYM.EQ.1) CALL MEMGET('INTE',KIBE(2,I),0,WORK,KFREE,LFREE)
C      ... always define KIBE(2,*) for PRPMLF call below
C          so we don't get undefined address
      ENDDO
C
C
C     Transform property matrices to MO basis
C     ---------------------------------------
C
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
      OPEN (LUMLF1,FORM='UNFORMATTED',FILE='MDPROP')
C
      CALL MEMGET('REAL',KPROPMO,NSTRT(1)*NSTRT(2)*NZ,WORK,KFREE,LFREE)
C
      DO I = 1,NTRPP
         INDXPR = LTRPP(I)
         ISYM  = IPRPSYM(INDXPR)
         IOPSY = JBTOF(ISYM-1,1)
         DO IX = 1, 32
            ACHAR(IX:IX) = '*'
         ENDDO
         ACHAR(25:32) = PRPNAM(INDXPR)(1:8)
         CALL MEMGET('REAL',KPRPMO,N2ORBXQ,WORK,KFREE,LFREE)
         CALL PRPMAT(INDXPR,IOPSY,WORK(KPRPMO),.TRUE.,WORK,WORK(KCMO),
     &               WORK(KIBEIG),ICMOQ,NORB,WORK,KFREE,LFREE,IPRTRP)
C
C        Select the set that we need
C
         CALL DZERO(WORK(KPROPMO),NSTRT(1)*NSTRT(2)*NZ)
C
         DO IFRP = 1, NFSYM
C
C          Determine fermion symmetry of indices for the bar and the ket.
C          We fix the ket index, so that the bra index follows from the 
C          combination of the symmetry of the ket and that of the operator.
C
           IFRP2 = IFRP
           IFRP1 = MOD(IFRP+IOPSY,2) + 1
C
C          Calculate offset in gathering matrix
C          Necessary when we have two fermion irreps
C
           IOFF = I2ORBX(IFRP1,IFRP2)
           JOFF = (IFRP2-1)*NSTRT(1)*NSTR(1,0,2)
     &          + (IFRP1-1)*NSTR(1,0,1)
           CALL SELMES(WORK(KPRPMO+IOFF),IFRP1,IFRP2,WORK(KPROPMO+JOFF),
     &                 NSTRT(1),NSTRT(2),
     &                 WORK(KVEC(IFRP1,1)),WORK(KVEC(IFRP2,2)),
     &                 NSTR(IFRP1,2,1),NSTR(IFRP2,2,2),
     &                 NSTR(IFRP1,1,1),NSTR(IFRP2,1,2))
         ENDDO
         CALL MEMREL('TRAPRP',WORK,KPRPMO,KPRPMO,KFREE,LFREE)
C
         CALL MEMGET('REAL',KPRCMO,NSTRT(1)*NSTRT(2)*8,WORK,KFREE,LFREE)
         CALL PRPMLF (WORK(KPRCMO),WORK(KPROPMO),IPRTRP,NSTR,
     &                NSTRT(1),NSTRT(2),
     &                WORK(KIBE(1,1)),WORK(KIBE(2,1)),
     &                WORK(KIBE(1,2)),WORK(KIBE(2,2)),
     &                ISYM,ACHAR,LUMLF1)
         CALL MEMREL('TRAPRP',WORK,KPRCMO,KPRCMO,KFREE,LFREE)
      ENDDO
      CLOSE(LU1INT,STATUS = 'KEEP')
C
      DO IX = 1, 32
         ACHAR(IX:IX) = '*'
      ENDDO
      ACHAR(25:32) = 'EOFLABEL'
      WRITE(LUMLF1) ACHAR
      CLOSE(LUMLF1,STATUS = 'KEEP')
C
      CALL MEMREL('TRAPRP',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRAPRP')
C
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SELMES(PRPMO,IFRP1,IFRP2,PBUF,NDIM1,NDIM2,JVEC1,JVEC2,
     &                  NPVEC1,NPVEC2,NEVEC1,NEVEC2)

C***********************************************************************
C     Pick out a set of matrix elements from PRPMO according to arrays
C     JVEC1 and JVEC2
C
C     Written by L.Visscher 1997
C     Last revision Jul 23 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION PRPMO(NORBT,NORBT,NZ),
     &          PBUF(NDIM1,NDIM2,NZ),
     &          JVEC1(*),JVEC2(*)
C
      DO IZ = 1,NZ
C
C       P-P matrix elements
C
        DO J = 1,NPVEC2
          JJ = NPSH(IFRP2)+1+JVEC2(J)
          DO I = 1,NPVEC1
             II = NPSH(IFRP1)+1+JVEC1(I)
             PBUF(I,J,IZ) = PRPMO(II,JJ,IZ)
          ENDDO
        ENDDO
C
C       E-P matrix elements
C
        DO J = 1,NPVEC2
          JJ = NPSH(IFRP2)+1+JVEC2(J)
          DO I = NPVEC1+1,NPVEC1+NEVEC1
             II = NPSH(IFRP1)+JVEC1(I)
             PBUF(I,J,IZ) = PRPMO(II,JJ,IZ)
          ENDDO
        ENDDO
C
C       P-E matrix elements
C
        DO J = NPVEC2+1,NPVEC2+NEVEC2
          JJ = NPSH(IFRP2)+JVEC2(J)
          DO I = 1,NPVEC1
             II = NPSH(IFRP1)+1+JVEC1(I)
             PBUF(I,J,IZ) = PRPMO(II,JJ,IZ)
          ENDDO
        ENDDO
C
C       E-E matrix elements
C
        DO J = NPVEC2+1,NPVEC2+NEVEC2
          JJ = NPSH(IFRP2)+JVEC2(J)
          DO I = NPVEC1+1,NPVEC1+NEVEC1
             II = NPSH(IFRP1)+JVEC1(I)
             PBUF(I,J,IZ) = PRPMO(II,JJ,IZ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C
C     ==============================================================
      SUBROUTINE PRPMLF (PRCMO,PROPMO,IPRINT,NSTR,NSTRT1,NSTRT2,
     &                   IBE11,IBE21,IBE12,IBE22,ISYM,
     &                   ACHAR,LUMLF1)
C     Express property integrals in complex form and write them to
C     file
C     Instead of array NSTRT now NSTRT1 and NSTRT2 are parameters
C     ==============================================================
C
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
      PARAMETER (DM1=-1.D0,ZEROTHRS=1.0D-11)
      CHARACTER*32 ACHAR
      DIMENSION NSTR(2,0:2,2),PRCMO(NSTRT1*2,NSTRT2*2,2)
      DIMENSION PROPMO(NSTRT1,NSTRT2,NZ)
      DIMENSION IOFF1(2,2),IOFF2(2,2),IOFF3(2,2)
      DIMENSION IBE11(*),IBE21(*),IBE12(*),IBE22(*)
!
C
      IF (IPRINT.GE.5) THEN
         WRITE(LUPRI,'(/3X,3A/)') '*** Property matrix :',ACHAR,' ***'
         CALL PRQMAT(PROPMO,NSTRT1,NSTRT2,NSTRT1,NSTRT2,
     &               NZ,IPQTOQ(1,ISYM-1),LUPRI)
      ENDIF
C
C     Fill in complete property matrix in complex form
C
      IOPSY = JBTOF(ISYM-1,1)
C
      IOFF1(1,1) = 0
      IOFF1(1,2) = 0
      IOFF2(1,1) = 0
      IOFF2(1,2) = 0
      IOFF3(1,1) = NSTR(1,0,1)
      IOFF3(1,2) = NSTR(1,0,2)
      DO IFRP = 2, NFSYM
         IOFF1(IFRP,1) = IOFF1(IFRP-1,1) + NSTR(IFRP-1,0,1)
         IOFF1(IFRP,2) = IOFF1(IFRP-1,2) + NSTR(IFRP-1,0,2)
         IOFF2(IFRP,1) = IOFF2(IFRP-1,1) + 2*NSTR(IFRP-1,0,1)
         IOFF2(IFRP,2) = IOFF2(IFRP-1,2) + 2*NSTR(IFRP-1,0,2)
         IOFF3(IFRP,1) = IOFF2(IFRP,1) + NSTR(IFRP,0,1)
         IOFF3(IFRP,2) = IOFF2(IFRP,2) + NSTR(IFRP,0,2)
      ENDDO
C
      CALL DZERO(PRCMO,NSTRT1*NSTRT2*8)
C
      DO IZ = 1, NZ
        IQ = IPQTOQ(IZ,ISYM-1)
        DO JFRP = 1, NFSYM
          IFRP  = MOD(JFRP+IOPSY,2) + 1
          DO J = 1, NSTR(JFRP,0,2)
            IF (IQ.EQ.1) THEN
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF2(IFRP,1),J+IOFF2(JFRP,2),1),1)
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF3(JFRP,2),1),1)
            ELSEIF (IQ.EQ.2) THEN
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF2(IFRP,1),J+IOFF2(JFRP,2),2),1)
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF3(JFRP,2),2),1)
              CALL DSCAL(NSTR(IFRP,0,1),DM1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF3(JFRP,2),2),1)
            ELSEIF (IQ.EQ.3) THEN
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF2(IFRP,1),J+IOFF3(JFRP,2),1),1)
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF2(JFRP,2),1),1)
              CALL DSCAL(NSTR(IFRP,0,1),DM1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF2(JFRP,2),1),1)
            ELSEIF (IQ.EQ.4) THEN
              CALL DCOPY(NSTR(IFRP,0,1),
     &                 PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                 PRCMO (1+IOFF2(IFRP,1),J+IOFF3(JFRP,2),2),1)
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF2(JFRP,2),2),1)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
C     The symmetry adapted functions were scaled by quaternion phases
C     to obtain as few active quaternion units as possible in double
C     group calculations (make NZ less than 4).
C     This scaling is, however, not convenient in spinfree calculations
C     because it makes real spinor products quaternion and thereby
C     hampers interfacing to nonrelativistic codes. We therefore
C     eliminate this scaling by multiplying the transformed property
C     matrix with the conjugate phase factors. The precise phase factor
C     depends on the boson identification of the spinor because it
C     should correspond to the factor on the non-zero coefficents.
C
C     A drawback is now that transformed property integrals are
C     real or imaginary depending on the BOSON symmetry identification
C     instead of the FERMION as was the case in the double group
C     adapted choice. This means that spinfree calculations with
C     RELCCSD that use non totally symmetric property matrices, i.e.
C     second order properties, should now be carried out in
C     complex arithmetic. Since these calculations will usually be
C     done by the dedicated response code instead of the RPA module
C     in RELCCSD this should be a minor problem.
C
C     A similar procedure is found in the 2-electron transformation,
C     but for the transformation of the core Fock matrix (and other
C     totally symmetric matrices) it not necessary to do this because
C     the combination with totally symmetric spinor products 
C     makes that the phase factors cancel. 
C
      IF (SPINFR) THEN
      CALL DZERO(PRCMO,NSTRT1*NSTRT2*8)
      DO IZ = 1, NZ
C       Determine the boson character of the operator and the
C       active quaternion unit.
        IQPRP = IPQTOQ(IZ,ISYM-1)
        DO JFRP = 1, NFSYM
          IFRP  = MOD(JFRP+IOPSY,2) + 1
          JJ1 = IOFF1(JFRP,2)
          JJ2 = IOFF2(JFRP,2)
          JJ3 = IOFF3(JFRP,2)
          DO J = 1, NSTR(JFRP,0,2)
             JJ1 = JJ1 + 1
             JJ2 = JJ2 + 1
             JJ3 = JJ3 + 1
             IF (JFRP.EQ.1) IREPJB = IBE12(J)
             IF (JFRP.EQ.2) IREPJB = IBE22(J)
             IQJ = IPQTOQ(1,IREPJB)
             II1 = IOFF1(IFRP,1)
             II2 = IOFF2(IFRP,1)
             II3 = IOFF3(IFRP,1)
             DO I = 1, NSTR(IFRP,0,1)
                II1 = II1 + 1
                II2 = II2 + 1
                II3 = II3 + 1
                IF (IFRP.EQ.1) IREPIB = IBE11(I)
                IF (IFRP.EQ.2) IREPIB = IBE21(I)
                IQI = IPQTOQ(1,IREPIB)
                IQ  = IQMULT (IQI,IQJ,IQPRP)
! LV
!   first amendment to make MO transformation of property integrals
!   work for the spinfree case:
! LV
! LV            IQS = IQPHASE(IQI,IQJ,IQPRP)*IQSIGN(IQJ,2,1)
                IQS = IQPHASE(IQI,IQJ,1)*IQSIGN(IQJ,2,1)
C
C The code is to be debugged. The MOLTRA part is probably OK but the
C phase transformed integrals do not give the correct answer when fed
C into RELCCSD. If you use this transformation with another code, you 
C can proceed at your own risk by uncommenting the stop.
C Stop now when we find non-symmetric operators.
!
! LV second amendment: we do not stop here, since it should work.
!
C  
! LV                IF (IQ.NE.IQS)
! LV      &          CALL QUIT ("This does not work with spinfree !")
! LV            if (iq.eq.1)
! LV &          write (*,'(7i4,f20.10)') i,j,iqi,iqj,iqprp,iq,iqs,
! LV &                                   propmo(ii1,jj1,iz)
C
C               Scale the matrix element
                PRPME = PROPMO(II1,JJ1,IZ) * IQS
C               Translate to the complex (bar/unbar) notation.
                IF (IQ.EQ.1) THEN
                   PRCMO(II2,JJ2,1) = PRPME
                   PRCMO(II3,JJ3,1) = PRPME
                ELSEIF (IQ.EQ.2) THEN
                   PRCMO(II2,JJ2,2) = PRPME
                   PRCMO(II3,JJ3,2) = - PRPME
                ELSEIF (IQ.EQ.3) THEN
                   PRCMO(II2,JJ3,1) = PRPME
                   PRCMO(II3,JJ2,1) = - PRPME
                ELSEIF (IQ.EQ.4) THEN
                   PRCMO(II2,JJ3,2) = PRPME
                   PRCMO(II3,JJ2,2) = PRPME
                ENDIF
             ENDDO
          ENDDO
        ENDDO
      ENDDO
      ENDIF
C
! Printout of complex matrix elements of the transformed operator
         IF (IPRINT.GE.2) THEN
           ! print out all NONZERO matrix elements
            WRITE (LUPRI,1012) ACHAR(25:32)
            IF ( IPRINT.EQ.2) THEN
              ! print nonzero elemenys
              WRITE(LUPRI,'(4x,a,d9.3)') 'zero threshold=',ZEROTHRS
            ELSE
              WRITE(LUPRI,"(A)")
     &        "Printout of all property elements"
            ENDIF
            DO I = 1, NSTRT1*2
               DO J = 1, NSTRT2*2
                IF (IPRINT.EQ.2) THEN
                  IF ( (DABS(PRCMO(I,J,1)).GT.ZEROTHRS).OR.
     &                 (DABS(PRCMO(I,J,2)).GT.ZEROTHRS) ) THEN
                    WRITE (LUPRI,1013) I,J,PRCMO(I,J,1),PRCMO(I,J,2)
                  ENDIF
                ELSE  ! print all elements
                    WRITE (LUPRI,1013) I,J,PRCMO(I,J,1),PRCMO(I,J,2)
                ENDIF
               ENDDO
            ENDDO
            WRITE(LUPRI,"(A)")
     &      "End of property elements printout."
            WRITE (6,*)
         ENDIF
C
C     WRITE PRCMO TO FILE
C
      WRITE (LUMLF1) ACHAR
      WRITE (LUMLF1) ((PRCMO(I,J,1),PRCMO(I,J,2),
     &                I=1,NSTRT1*2),J=1,NSTRT2*2)
C
 1000 FORMAT (/' Transforming ',A8,' evaluation date ',A8,' label ',A8)
 1001 FORMAT (/' Finished transformation : integrals on ',A8)
 1010 FORMAT (/'Nonzero diagonal elements of ',A8,
     & ' operator in Kramers paired molecular spinor basis:'
     & //12X,'Real part',11X,'Imaginary Part')
 1011 FORMAT (I6,2F20.10)
 1012 FORMAT (/'Property elements of the "',A8,'" operator',
     &' in Kramers paired molecular spinor basis:'/)
 1013 FORMAT (2I5,1P,2E15.6)
      RETURN
      END
