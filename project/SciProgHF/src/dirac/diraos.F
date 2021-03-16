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

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck getrho */
      SUBROUTINE GETRHO(RHO,IREPOP,GAO,DMAT,BUF)
C***********************************************************************
C
C     Get density in point (PX,PY,PZ).
C     Gets values of AOs in this point from GAO
C
C     Written by T.Saue and T.Helgaker
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),DMAT(NBAST,NBAST),BUF(*)
#else
      DIMENSION GAO(NTBAS(0)),DMAT(NTBAS(0),NTBAS(0)),BUF(*)
#endif
      RHO = D0
      DO IC = 1,MC
         CALL GETRH1(RHO,IC,IREPOP,GAO,DMAT,BUF)      
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getrh1 */
      SUBROUTINE GETRH1(RHO,IC,IREPOP,GAO,DMAT,BUF)      
C***********************************************************************
C
C     Get density in point (PX,PY,PZ).
C     Gets values of AOs in this point from GAO
C
C     Written by T.Saue and T.Helgaker
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),DMAT(NBAST,NBAST),BUF(*)
#else
      DIMENSION GAO(NTBAS(0)),DMAT(NTBAS(0),NTBAS(0)),BUF(*)
#endif
C
#include "ibtfun.h"
C      
      DO 10 IRP1 = 0,NBSYM-1
        NB1 = NBBAS(IRP1,IC)
        IF(NB1.EQ.0) GOTO 10
        JB1 = IBBAS(IRP1,IC)+1
        IRP2 = IBTXOR(IRP1,IREPOP)
        NB2 = NBBAS(IRP2,IC)
        IF(NB2.EQ.0) GOTO 10
        JB2 = IBBAS(IRP2,IC)+1
        CALL DGEMV('N',NBBAS(IRP1,IC),NBBAS(IRP2,IC),D1,
     &             DMAT(JB1,JB2),NTBAS(0),GAO(JB2),1,D0,BUF,1)
        RHO = RHO + DDOT(NBBAS(IRP1,IC),GAO(JB1),1,BUF,1)
   10 CONTINUE
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getgrho */
      SUBROUTINE GETGRHO(RHG,IREPOP,GAO,GAO1,DMAT,BUF)      
C***********************************************************************
C
C     Get gradient of density in point (PX,PY,PZ).
C     Gets values of AOs in this point from GAO
C
C     WARNING: This routine assumes Hermitian density matrix !
C
C     Written by T.Saue and T.Helgaker
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),DMAT(NBAST,NBAST),
     &          BUF(NBAST),GAO1(NBAST,3)
#else
      DIMENSION GAO(NTBAS(0)),DMAT(NTBAS(0),NTBAS(0)),
     &          BUF(NTBAS(0)),GAO1(NTBAS(0),3)
#endif
      DIMENSION RHG(3)
C      
#include "ibtfun.h"
C      
      DO IAX = 1,3
        RHG(IAX) = D0
      ENDDO
      DO IC = 1,MC
        DO 10 IRP1 = 0,NBSYM-1
          NB1  = NBBAS(IRP1,IC)
          IF(NB1.EQ.0) GOTO 10
          JB1  = IBBAS(IRP1,IC)+1
          IRP2 = IBTXOR(IRP1,IREPOP)
          NB2  = NBBAS(IRP2,IC)
          IF(NB2.EQ.0) GOTO 10
          JB2 = IBBAS(IRP2,IC)+1
          CALL DGEMV('N',NBBAS(IRP1,IC),NBBAS(IRP2,IC),D1,
     &                 DMAT(JB1,JB2),NTBAS(0),GAO(JB2),1,D0,BUF,1)
          DO IAX = 1,3
            RHG(IAX) = RHG(IAX) 
     &             + DDOT(NBBAS(IRP1,IC),GAO1(JB1,IAX),1,BUF,1)
          ENDDO
   10   CONTINUE
      ENDDO
      DO IAX = 1,3
        RHG(IAX) = RHG(IAX)+RHG(IAX)
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getlrho */
      SUBROUTINE GETLRHO(RXX,RXY,RXZ,RYY,RYZ,RZZ,GAO,GAO1,GAO2,
     &                   DMAT,BUF)      
C***********************************************************************
C
C     Get Hessian of density in point (PX,PY,PZ).
C     Gets values of AOs in this point from GAO
C
C     Written by T.Saue and T.Helgaker
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
C
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION RHO(0:2),GAO(NBAST),DMAT(NBAST,NBAST),
     &          BUF(NBAST,4),GAO1(NBAST,3),GAO2(NBAST,6)
#else
      DIMENSION RHO(0:2),GAO(NTBAS(0)),DMAT(NTBAS(0),NTBAS(0)),
     &          BUF(NTBAS(0),4),GAO1(NTBAS(0),3),GAO2(NTBAS(0),6)
#endif
C      
      RXX = D0
      RXY = D0
      RXZ = D0
      RYY = D0
      RYZ = D0
      RZZ = D0
C
      NBRP = 4/NZ
      DO IFRP = 1,NFSYM
        DO IC = 1,MC
          IP = MOD(IFRP+IC,NFSYM) + 1
          DO JSYM = 1,NBRP
            IREP = JFSYM(JSYM,IP) - 1
            IF(NBBAS(IREP,IC).GT.0) THEN
              JJ = IBBAS(IREP,IC)
              DO J = 1, NBBAS(IREP,IC)
                BUF(J,1) = DDOT(NBBAS(IREP,IC),
     &                          DMAT(JJ+1,JJ+J),1,GAO(JJ+1),1)
                BUF(J,2) = DDOT(NBBAS(IREP,IC),
     &                          DMAT(JJ+1,JJ+J),1,GAO1(JJ+1,1),1)
                BUF(J,3) = DDOT(NBBAS(IREP,IC),
     &                          DMAT(JJ+1,JJ+J),1,GAO1(JJ+1,2),1)
                BUF(J,4) = DDOT(NBBAS(IREP,IC),
     &                          DMAT(JJ+1,JJ+J),1,GAO1(JJ+1,3),1)
              ENDDO
              RXX = RXX 
     &            + DDOT(NBBAS(IREP,IC),BUF(1,1),1,GAO2(JJ+1,1),1)
     &            + DDOT(NBBAS(IREP,IC),BUF(1,2),1,GAO1(JJ+1,1),1)
              RXY = RXY 
     &            + DDOT(NBBAS(IREP,IC),BUF(1,1),1,GAO2(JJ+1,2),1)
     &            + DDOT(NBBAS(IREP,IC),BUF(1,2),1,GAO1(JJ+1,2),1)
     &            + DDOT(NBBAS(IREP,IC),BUF(1,3),1,GAO1(JJ+1,1),1)
              RXZ = RXZ 
     &            + DDOT(NBBAS(IREP,IC),BUF(1,1),1,GAO2(JJ+1,3),1)
     &            + DDOT(NBBAS(IREP,IC),BUF(1,2),1,GAO1(JJ+1,3),1)
     &            + DDOT(NBBAS(IREP,IC),BUF(1,4),1,GAO1(JJ+1,1),1)
              RYY = RYY 
     &            + DDOT(NBBAS(IREP,IC),BUF(1,1),1,GAO2(JJ+1,4),1)
     &            + DDOT(NBBAS(IREP,IC),BUF(1,3),1,GAO1(JJ+1,2),1)
              RYZ = RYZ 
     &            + DDOT(NBBAS(IREP,IC),BUF(1,1),1,GAO2(JJ+1,5),1)
     &            + DDOT(NBBAS(IREP,IC),BUF(1,3),1,GAO1(JJ+1,3),1)
     &            + DDOT(NBBAS(IREP,IC),BUF(1,4),1,GAO1(JJ+1,2),1)
              RZZ = RZZ 
     &            + DDOT(NBBAS(IREP,IC),BUF(1,1),1,GAO2(JJ+1,6),1)
     &            + DDOT(NBBAS(IREP,IC),BUF(1,4),1,GAO1(JJ+1,3),1)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RXX = RXX + RXX
      RXY = RXY + RXY
      RXZ = RXZ + RXZ
      RYY = RYY + RYY
      RYZ = RYZ + RYZ
      RZZ = RZZ + RZZ
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getsos */
      SUBROUTINE GETSOS(GAO,GAO1,GAO2,GAB1,NCNT,CORPX,CORPY,CORPZ,
     &                  BUF,NBAST,NDER,DOLND,IPRINT)
C***********************************************************************
C
C     Calculate value of symmetry-adapted basis functions in a given
C     point (CORPX,CORPY,CORPZ)
C
C     Written by T.Saue May 12 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
      DIMENSION GAO(NBAST), GAO1(NBAST,3), GAO2(NBAST,6),
     &          GAB1(NBAST,3), BUF(NBAST), NCNT(NBAST)
C
C     First get values for the Ao-basis
      CALL GETAOS(GAO,GAO1,GAO2,GAB1,NCNT,CORPX,CORPY,CORPZ,
     &            NBAST,NDER,DOLND,IPRINT)
      IF(NBSYM.EQ.1) RETURN
C
C     Functions
      CALL WTAOSO(GAO,BUF,NBAST,1,IPRINT)
      CALL VTBUBS(BUF,GAO,1)
C
C     First derivatives
      IF(NDER.GE.1) THEN
        DO J = 1,3
          CALL WTAOSO(GAO1(1,J),BUF,NBAST,1,IPRINT)
          CALL VTBUBS(BUF,GAO1(1,J),1)
        ENDDO
      ENDIF
C
C     Second derivatives
      IF(NDER.GE.2) THEN
        DO J = 1,6
          CALL WTAOSO(GAO2(1,J),BUF,NBAST,1,IPRINT)
          CALL VTBUBS(BUF,GAO2(1,J),1)
        ENDDO
      ENDIF
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      END
C  /* Deck getaos */
      SUBROUTINE GETAOS(GAO,GAO1,GAO2,GAB1,NCNT,CORPX,CORPY,CORPZ,
     &                  NBAST,NDER,DOLND,IPRINT)
C
C     T. Helgaker sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C   
      LOGICAL DOLND
      DIMENSION GAO(NBAST), GAO1(NBAST,3), GAO2(NBAST,6),
     &          GAB1(NBAST,3), NCNT(NBAST)
C
#include "dcbham.h"
#include "onecom.h"
#include "lmns.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "ibtfun.h"
      XYZ(I,J) = PT(IBTAND(ISYMAX(I,1),J))
C
      CALL DZERO(GAO,NBAST)
      IF (NDER.GT.0) CALL DZERO(GAO1,3*NBAST)
      IF (NDER.GT.1) CALL DZERO(GAO2,6*NBAST)
      IADR = 1
      ISHBEGIN=1
      ISHEND  =KMAX

!     if Levy-Leblond or IOTC
!     don't work on small component shells

!radovan:
!     comment out this check (= include all shells) if you wish
!     to calculate the relativistic current with Levy-Leblond
      IF (MC .EQ. 1) ISHEND = NLRGSH
C
      DO 100 ISHELA = ISHBEGIN,ISHEND
         NHKTA = NHKT(ISHELA)
         KHKTA = KHKT(ISHELA)
         KCKTA = KCKT(ISHELA)
         SPHRA = SPHR(ISHELA)
         IF (NDER.GT.1 .OR. NHKTA.GT.2 .OR. DOLND) THEN
            CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
         END IF
         NUCA   = NUCO(ISHELA)
         NUMCFA = NUMCF(ISHELA)
         JSTA   = JSTRT(ISHELA)
C
Ctsaue: Loop over symmetry dependent centers   
C
         MULA   = ISTBAO(ISHELA)
         DO ISYMOP = 0,MAXOPR
         IF (IBTAND(ISYMOP,MULA) .EQ. 0) THEN
           PAX = CORPX-PT(IBTAND(ISYMAX(1,1),ISYMOP))*CENT(ISHELA,1,1)
           PAY = CORPY-PT(IBTAND(ISYMAX(2,1),ISYMOP))*CENT(ISHELA,2,1)
           PAZ = CORPZ-PT(IBTAND(ISYMAX(3,1),ISYMOP))*CENT(ISHELA,3,1)
           IF (NDER.EQ.0) CALL GETGAO(GAO(IADR),PAX,PAY,PAZ)
           IF (NDER.GT.0) THEN
             CALL GETGA1(GAO(IADR),
     &                   GAO1(IADR,1),GAO1(IADR,2),GAO1(IADR,3),
     &                   PAX,PAY,PAZ)
           END IF
           IF (NDER.GT.1) THEN
             CALL GETGA2(GAO2(IADR,1),GAO2(IADR,2),GAO2(IADR,3),
     &                  GAO2(IADR,4),GAO2(IADR,5),GAO2(IADR,6),
     &                  PAX,PAY,PAZ)
           END IF
           IF (.false.) THEN
             CALL GETGB1(GAO(IADR),
     &                   GAB1(IADR,1),GAB1(IADR,2),GAB1(IADR,3),
     &                   PAX,PAY,PAZ,CORPX,CORPY,CORPZ)
           END IF
           IF (NDER.GT.0 .OR. DOLND) THEN
              NCNTA = NCENT(ISHELA)
              DO 300 ICMP = 0, KHKTA - 1
                NCNT(IADR + ICMP) = NCNTA
  300         CONTINUE
           END IF
           IADR = IADR + KHKTA
        ENDIF
        ENDDO
  100 CONTINUE
C
C     Print section
C
      IF (IPRINT .GE. 50) THEN
          CALL HEADER('Output from GETAOS',-1)
          WRITE (LUPRI,'(A,3F12.6)') ' CORP ', CORPX,CORPY,CORPZ 
          IF (IPRINT .GT. 200) THEN
             CALL HEADER('Undifferentiated integrals at this point:',-1)
             CALL OUTPUT(GAO,1,1,1,NBAST,1,NBAST,-1,LUPRI)
             IF (NDER.GT.0) THEN
               CALL HEADER('1st derivative integrals at this point:',-1)
               CALL OUTPUT(GAO1,1,NBAST,1,3,NBAST,3,-1,LUPRI)
             END IF
             IF (NDER.GT.1) THEN
               CALL HEADER('2nd derivative integrals at this point:',-1)
               CALL OUTPUT(GAO2,1,NBAST,1,6,NBAST,6,-1,LUPRI)
             END IF
             IF (DOLND) THEN
               CALL HEADER('London der. integrals at this point:',-1)
               CALL OUTPUT(GAB1,1,NBAST,1,3,NBAST,3,-1,LUPRI)
             END IF
          END IF
          WRITE (LUPRI,'(1X,A)') ' '
      END IF
      RETURN
      END
C  /* Deck getgao */
      SUBROUTINE GETGAO(GAO,PAX,PAY,PAZ)
C
C     T. Helgaker Sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "maxaqn.h"
C
      DIMENSION GAO(KCKTA), CAO(MXAQN)
C
#include "lmns.h"
#include "onecom.h"
#include "primit.h"
#include "sphtrm.h"
C
      ISTR = JSTA + 1
      IEND = JSTA + NUCA
      PA2 = PAX**2 + PAY**2 + PAZ**2
      IF (KCKTA .EQ. 1) THEN
         GAO(1) = PRICCF(ISTR,NUMCFA)*DEXP(-PRIEXP(ISTR)*PA2)
         DO 100 IPRIMA = ISTR + 1, IEND
            GAO(1) = GAO(1) + PRICCF(IPRIMA,NUMCFA)
     &                   *DEXP(-PRIEXP(IPRIMA)*PA2)
  100    CONTINUE
      ELSE IF (KCKTA .EQ. 3) THEN
         GA = PRICCF(ISTR,NUMCFA)*DEXP(-PRIEXP(ISTR)*PA2)
         GAO(1) = PAX*GA
         GAO(2) = PAY*GA
         GAO(3) = PAZ*GA
         DO 200 IPRIMA = ISTR + 1, IEND
            GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-PRIEXP(IPRIMA)*PA2)
            GAO(1) = GAO(1) + PAX*GA
            GAO(2) = GAO(2) + PAY*GA
            GAO(3) = GAO(3) + PAZ*GA
  200    CONTINUE
      ELSE
         IF (SPHRA) THEN
            GA = PRICCF(ISTR,NUMCFA)*DEXP(-PRIEXP(ISTR)*PA2)
            DO 300 ICOMPA = 1,KCKTA
               CAO(ICOMPA) = (PAX**LVALUA(ICOMPA))*
     &                       (PAY**MVALUA(ICOMPA))*
     &                       (PAZ**NVALUA(ICOMPA))*GA
  300       CONTINUE
            DO 400 IPRIMA = ISTR + 1, IEND
               GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-PRIEXP(IPRIMA)*PA2)
               DO 500 ICOMPA = 1,KCKTA
                  CAO(ICOMPA) = CAO(ICOMPA) + (PAX**LVALUA(ICOMPA))*
     &                                        (PAY**MVALUA(ICOMPA))*
     &                                        (PAZ**NVALUA(ICOMPA))*GA
  500          CONTINUE
  400       CONTINUE
C           CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAO,KCKTA,GAO,1)
            CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAO,1,KCKTA,GAO,1,
     &                KHKTA,KHKTA,KCKTA,1)
         ELSE
            GA = PRICCF(ISTR,NUMCFA)*DEXP(-PRIEXP(ISTR)*PA2)
            DO 600 ICOMPA = 1,KCKTA
               GAO(ICOMPA) = (PAX**LVALUA(ICOMPA))*
     &                       (PAY**MVALUA(ICOMPA))*
     &                       (PAZ**NVALUA(ICOMPA))*GA
  600       CONTINUE
            DO 700 IPRIMA = ISTR + 1, IEND
               GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-PRIEXP(IPRIMA)*PA2)
               DO 800 ICOMPA = 1,KCKTA
                  GAO(ICOMPA) = GAO(ICOMPA) + (PAX**LVALUA(ICOMPA))*
     &                                        (PAY**MVALUA(ICOMPA))*
     &                                        (PAZ**NVALUA(ICOMPA))*GA
  800          CONTINUE
  700       CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck getga1 */
      SUBROUTINE GETGA1(GAO,GAOX,GAOY,GAOZ,PAX,PAY,PAZ)
C
C     T. Helgaker Sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (D2 = 2.0D0) 
C
      LOGICAL FIRST
      DIMENSION GAO(KCKTA), GAOX(KCKTA), GAOY(KCKTA), GAOZ(KCKTA),
     &          CAO(MXAQN), CAOX(MXAQN), CAOY(MXAQN), CAOZ(MXAQN)
C
#include "lmns.h"
#include "onecom.h"
#include "primit.h"
#include "sphtrm.h"
C
C
      PA2 = PAX**2 + PAY**2 + PAZ**2
      FIRST = .TRUE.
      DO 100 IPRIMA = JSTA + 1, JSTA + NUCA
         ALPHA = PRIEXP(IPRIMA)
         TALPH = -D2*ALPHA
         GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-ALPHA*PA2)
         IF (KCKTA .EQ. 1) THEN
            TGA = TALPH*GA
            IF (FIRST) THEN
               GAO(1)  = GA
               GAOX(1) = PAX*TGA
               GAOY(1) = PAY*TGA
               GAOZ(1) = PAZ*TGA
            ELSE
               GAO(1)  = GAO(1)  + GA
               GAOX(1) = GAOX(1) + PAX*TGA
               GAOY(1) = GAOY(1) + PAY*TGA
               GAOZ(1) = GAOZ(1) + PAZ*TGA
            END IF
         ELSE IF (KCKTA .EQ. 3) THEN
            TGA = TALPH*GA 
            TGX = PAX*TGA 
            TGY = PAY*TGA 
            TGZ = PAZ*TGA
            IF (FIRST) THEN
               GAO(1)  = PAX*GA
               GAO(2)  = PAY*GA
               GAO(3)  = PAZ*GA
               GAOX(1) = PAX*TGX + GA
               GAOX(2) = PAY*TGX
               GAOX(3) = PAZ*TGX
               GAOY(1) = PAX*TGY
               GAOY(2) = PAY*TGY + GA
               GAOY(3) = PAZ*TGY
               GAOZ(1) = PAX*TGZ
               GAOZ(2) = PAY*TGZ
               GAOZ(3) = PAZ*TGZ + GA
            ELSE
               GAO(1)  = GAO(1)  + PAX*GA
               GAO(2)  = GAO(2)  + PAY*GA
               GAO(3)  = GAO(3)  + PAZ*GA
               GAOX(1) = GAOX(1) + PAX*TGX + GA
               GAOX(2) = GAOX(2) + PAY*TGX
               GAOX(3) = GAOX(3) + PAZ*TGX
               GAOY(1) = GAOY(1) + PAX*TGY
               GAOY(2) = GAOY(2) + PAY*TGY + GA
               GAOY(3) = GAOY(3) + PAZ*TGY
               GAOZ(1) = GAOZ(1) + PAX*TGZ
               GAOZ(2) = GAOZ(2) + PAY*TGZ
               GAOZ(3) = GAOZ(3) + PAZ*TGZ + GA
            END IF
         ELSE
            DO 200 ICOMPA = 1,KCKTA
               L = LVALUA(ICOMPA)
               M = MVALUA(ICOMPA)
               N = NVALUA(ICOMPA)
               P0  = (PAX**L)*(PAY**M)*(PAZ**N)
C
C              up contributions 
C
               TP0 = TALPH*P0
               GAX = PAX*TP0
               GAY = PAY*TP0
               GAZ = PAZ*TP0
C
C              down contributions
C
               IF(L.GT.0)GAX=GAX+dble(L)*(PAX**(L-1))*(PAY**M)*(PAZ**N)
               IF(M.GT.0)GAY=GAY+dble(M)*(PAX**L)*(PAY**(M-1))*(PAZ**N)
               IF(N.GT.0)GAZ=GAZ+dble(N)*(PAX**L)*(PAY**M)*(PAZ**(N-1))
C
               IF (SPHRA) THEN
                  IF (FIRST) THEN
                     CAO (ICOMPA) =  P0*GA
                     CAOX(ICOMPA) = GAX*GA
                     CAOY(ICOMPA) = GAY*GA
                     CAOZ(ICOMPA) = GAZ*GA
                  ELSE
                     CAO (ICOMPA) = CAO (ICOMPA) +  P0*GA
                     CAOX(ICOMPA) = CAOX(ICOMPA) + GAX*GA
                     CAOY(ICOMPA) = CAOY(ICOMPA) + GAY*GA
                     CAOZ(ICOMPA) = CAOZ(ICOMPA) + GAZ*GA
                  END IF
               ELSE
                  IF (FIRST) THEN
                     GAO (ICOMPA) =   P0*GA
                     GAOX(ICOMPA) =  GAX*GA
                     GAOY(ICOMPA) =  GAY*GA
                     GAOZ(ICOMPA) =  GAZ*GA
                  ELSE
                     GAO (ICOMPA) = GAO (ICOMPA) +  P0*GA
                     GAOX(ICOMPA) = GAOX(ICOMPA) + GAX*GA
                     GAOY(ICOMPA) = GAOY(ICOMPA) + GAY*GA
                     GAOZ(ICOMPA) = GAOZ(ICOMPA) + GAZ*GA
                  END IF
               END IF
  200       CONTINUE
         END IF
         FIRST = .FALSE.
  100 CONTINUE
      IF (SPHRA) THEN
C         CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAO ,KCKTA,GAO ,1)
C         CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAOX,KCKTA,GAOX,1)
C         CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAOY,KCKTA,GAOY,1)
C         CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAOZ,KCKTA,GAOZ,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAO ,1,KCKTA,GAO ,1,
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOX,1,KCKTA,GAOX,1,
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOY,1,KCKTA,GAOY,1,
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOZ,1,KCKTA,GAOZ,1,
     &              KHKTA,KHKTA,KCKTA,1)
      END IF
      RETURN
      END
C  /* Deck getga2 */
      SUBROUTINE GETGA2(GAOXX,GAOXY,GAOXZ,GAOYY,GAOYZ,GAOZZ,PAX,PAY,PAZ)
C
C     T. Helgaker Sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "maxaqn.h"
      LOGICAL FIRST
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0) 
C
      DIMENSION GAOXX(KCKTA), GAOXY(KCKTA), GAOXZ(KCKTA),
     &          GAOYY(KCKTA), GAOYZ(KCKTA), GAOZZ(KCKTA),
     &          CAOXX(MXAQN), CAOXY(MXAQN), CAOXZ(MXAQN),
     &          CAOYY(MXAQN), CAOYZ(MXAQN), CAOZZ(MXAQN)
C
#include "lmns.h"
#include "onecom.h"
#include "primit.h"
#include "sphtrm.h"
C
C
      PA2 = PAX**2 + PAY**2 + PAZ**2
      FIRST = .TRUE.   
      DO 100 IPRIMA = JSTA + 1, JSTA + NUCA
         ALPHA = PRIEXP(IPRIMA)
         TALPH = -D2*ALPHA
         TAPAX = TALPH*PAX
         TAPAY = TALPH*PAY
         TAPAZ = TALPH*PAZ
         GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-ALPHA*PA2)
         DO 200 ICOMPA = 1,KCKTA
            L = LVALUA(ICOMPA)
            M = MVALUA(ICOMPA)
            N = NVALUA(ICOMPA)
C
            PXD = D0
            PYD = D0
            PZD = D0
            IF (L.GT.1) PXD = dble(L*(L-1))*(PAX**(L-2))
            IF (M.GT.1) PYD = dble(M*(M-1))*(PAY**(M-2))
            IF (N.GT.1) PZD = dble(N*(N-1))*(PAZ**(N-2))
            PXM = D0
            PYM = D0
            PZM = D0
            IF (L.GT.0) PXM = dble(L)*(PAX**(L-1))
            IF (M.GT.0) PYM = dble(M)*(PAY**(M-1))
            IF (N.GT.0) PZM = dble(N)*(PAZ**(N-1))
            PX0 = PAX**L
            PY0 = PAY**M
            PZ0 = PAZ**N
            PXP = TAPAX*PX0
            PYP = TAPAY*PY0
            PZP = TAPAZ*PZ0
            P000 = PX0*PY0*PZ0 
C
            GAXX = (TAPAX**2 + TALPH*dble(2*L + 1))*P000 + PXD*PY0*PZ0
            GAYY = (TAPAY**2 + TALPH*dble(2*M + 1))*P000 + PX0*PYD*PZ0
            GAZZ = (TAPAZ**2 + TALPH*dble(2*N + 1))*P000 + PX0*PY0*PZD
            GAXY = TAPAX*TAPAY*P000 + (PXP*PYM + PXM*PYP + PXM*PYM)*PZ0
            GAXZ = TAPAX*TAPAZ*P000 + (PXP*PZM + PXM*PZP + PXM*PZM)*PY0
            GAYZ = TAPAY*TAPAZ*P000 + (PYP*PZM + PYM*PZP + PYM*PZM)*PX0
C
            IF (SPHRA) THEN
               IF (FIRST) THEN
                  CAOXX(ICOMPA) = GAXX*GA 
                  CAOXY(ICOMPA) = GAXY*GA 
                  CAOXZ(ICOMPA) = GAXZ*GA
                  CAOYY(ICOMPA) = GAYY*GA
                  CAOYZ(ICOMPA) = GAYZ*GA
                  CAOZZ(ICOMPA) = GAZZ*GA
               ELSE
                  CAOXX(ICOMPA) = CAOXX(ICOMPA) + GAXX*GA 
                  CAOXY(ICOMPA) = CAOXY(ICOMPA) + GAXY*GA 
                  CAOXZ(ICOMPA) = CAOXZ(ICOMPA) + GAXZ*GA
                  CAOYY(ICOMPA) = CAOYY(ICOMPA) + GAYY*GA
                  CAOYZ(ICOMPA) = CAOYZ(ICOMPA) + GAYZ*GA
                  CAOZZ(ICOMPA) = CAOZZ(ICOMPA) + GAZZ*GA
               END IF
            ELSE
               IF (FIRST) THEN
                  GAOXX(ICOMPA) = GAXX*GA 
                  GAOXY(ICOMPA) = GAXY*GA 
                  GAOXZ(ICOMPA) = GAXZ*GA
                  GAOYY(ICOMPA) = GAYY*GA
                  GAOYZ(ICOMPA) = GAYZ*GA
                  GAOZZ(ICOMPA) = GAZZ*GA
               ELSE
                  GAOXX(ICOMPA) = GAOXX(ICOMPA) + GAXX*GA 
                  GAOXY(ICOMPA) = GAOXY(ICOMPA) + GAXY*GA 
                  GAOXZ(ICOMPA) = GAOXZ(ICOMPA) + GAXZ*GA
                  GAOYY(ICOMPA) = GAOYY(ICOMPA) + GAYY*GA
                  GAOYZ(ICOMPA) = GAOYZ(ICOMPA) + GAYZ*GA
                  GAOZZ(ICOMPA) = GAOZZ(ICOMPA) + GAZZ*GA
               END IF
            END IF
  200    CONTINUE
         FIRST = .FALSE.
  100 CONTINUE
      IF (SPHRA) THEN
C         CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAOXX,KCKTA,GAOXX,1)
C         CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAOXY,KCKTA,GAOXY,1)
C         CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAOXZ,KCKTA,GAOXZ,1)
C         CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAOYY,KCKTA,GAOYY,1)
C         CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAOYZ,KCKTA,GAOYZ,1)
C         CALL MXM(CSP(ISPADR(NHKTA)),KHKTA,CAOZZ,KCKTA,GAOZZ,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOXX,1,KCKTA,GAOXX,1,
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOXY,1,KCKTA,GAOXY,1,
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOXZ,1,KCKTA,GAOXZ,1,
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOYY,1,KCKTA,GAOYY,1,
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOYZ,1,KCKTA,GAOYZ,1,
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOZZ,1,KCKTA,GAOZZ,1,
     &              KHKTA,KHKTA,KCKTA,1)
      END IF
      RETURN
      END
C  /* Deck getgb1 */
      SUBROUTINE GETGB1(GAO,GAOX,GAOY,GAOZ,PAX,PAY,PAZ,PX,PY,PZ)
C
C     T. Helgaker Oct 99
C
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "maxaqn.h"
C
      PARAMETER (DP5 = 0.5D0)
      DIMENSION GAOX(KCKTA), GAOY(KCKTA), GAOZ(KCKTA), GAO(KCKTA)
C
#include "onecom.h"
#include "primit.h"
#include "orgcom.h"
C
      POX = PX - ORIGIN(1)
      POY = PY - ORIGIN(2)
      POZ = PZ - ORIGIN(3)
      FCX = -DP5*(PAY*POZ - PAZ*POY)
      FCY = -DP5*(PAZ*POX - PAX*POZ)
      FCZ = -DP5*(PAX*POY - PAY*POX)
      DO 100 ICOMPA = 1,KCKTA
         GA = GAO(ICOMPA)
         GAOX(ICOMPA) = FCX*GA
         GAOY(ICOMPA) = FCY*GA
         GAOZ(ICOMPA) = FCZ*GA
  100 CONTINUE
      RETURN
      END
