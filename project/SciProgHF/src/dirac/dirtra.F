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
C  /* Deck getlab */
      SUBROUTINE GETLAB(IPRINT)
C***********************************************************************
C
C      Generate primitive labels for basis functions in dcblab.h
C
C      Called from PAMINP
C
C      Written by T.Saue - November 1994
C      Last revision: Nov 23 1994 - tsaue
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "dcblab.h"
#include "symmet.h"
#include "aosotr.h"
#include "dgroup.h"
#include "pgroup.h"
#include "nuclei.h"
#include "shells.h"
#include "dcbbas.h"
#include "ccom.h"
      CHARACTER CLS(2)*1,DEG*2
      integer, allocatable :: kp(:)
C
#include "dcbibt.h"
      CALL QENTER('GETLAB')
      CLS(1) = 'L'
      CLS(2) = 'S'
C
C     **********************************************
C     ***** L A B E L S  for   A O - B A S I S *****
C     **********************************************
C
C
C     Make labels in the order they appear in HERMIT
C     ==============================================
C
      NLAB  = 0
      NAORB = 0
      DO 10 ISHELL = 1,NLRGSH
        LVAL  = NHKT(ISHELL)
        ICENT = NCENT(ISHELL)
        ICLS  = LCLASS(ISHELL)
        NDEG  = NUCDEG(ICENT)
        DO 20 IDEG = 1,NDEG
          IOFF = NHKOFF(LVAL)
          DO 30 ICOMP = 1,KHK(LVAL)
            ITYP  = IOFF + ICOMP
            NAORB = NAORB + 1
            JLAB  = IPACK(ICLS,ICENT,IDEG,ITYP)
C
C           Check if label is already defined
C           =================================
            DO 40 ILAB = 1,NLAB
            IF(JLAB.EQ.IATTR(ILAB,1)) THEN
              IPLAB(NAORB,1) = ILAB
              GOTO 30
            ENDIF
C
 40         CONTINUE
            NLAB = NLAB + 1
            IF(NLAB.GT.MAXLAB) GOTO 999
            WRITE(DEG,'(I2)') IDEG
            IATTR(NLAB,1) = JLAB
            PLABEL(NLAB,1) = CLS(ICLS)//' '//NAMN(ICENT)(1:3)//
     &          DEG//' '//GTOTYP(ITYP)
            IPLAB(NAORB,1) = NLAB
 30       CONTINUE
 20     CONTINUE
 10   CONTINUE
      NPLAB(1) = NLAB
      DO 50 ISHELL = NLRGSH+1,KMAX
        LVAL  = NHKT(ISHELL)
        ICENT = NCENT(ISHELL)
        ICLS  = LCLASS(ISHELL)
        NDEG  = NUCDEG(ICENT)
        DO 60 IDEG = 1,NDEG
          IOFF = NHKOFF(LVAL)
          DO 70 ICOMP = 1,KHK(LVAL)
            ITYP  = IOFF + ICOMP
            NAORB = NAORB + 1
            JLAB  = IPACK(ICLS,ICENT,IDEG,ITYP)
C
C           Check if label is already defined
C           =================================
            DO 80 ILAB = 1,NLAB
            IF(JLAB.EQ.IATTR(ILAB,1)) THEN
              IPLAB(NAORB,1) = ILAB
              GOTO 70
            ENDIF
 80         CONTINUE
C
C           Define new label
C           ================
            NLAB = NLAB + 1
            IF(NLAB.GT.MAXLAB) GOTO 999
            WRITE(DEG,'(I2)') IDEG
            IATTR(NLAB,1) = JLAB
            PLABEL(NLAB,1) = CLS(ICLS)//' '//NAMN(ICENT)(1:3)//
     &         DEG//' '//GTOTYP(ITYP)
            IPLAB(NAORB,1) = NLAB
 70       CONTINUE
 60     CONTINUE
 50   CONTINUE
      NPLAB(0) = NLAB
      NPLAB(2) = NPLAB(0) - NPLAB(1)
C
C     Print section
C     =============
C
      IF(IPRINT.GE.1) THEN
        CALL HEADER('GETLAB: AO-labels',-1)
        WRITE(LUPRI,'(3X,A,I5)') '* Large components:',NPLAB(1)
        WRITE(LUPRI,'(6(I6,2X,A12))') (I,PLABEL(I,1),I=1,NPLAB(1))
        IOFF = NPLAB(1)
        WRITE(LUPRI,'(3X,A,I5)') '* Small components:',NPLAB(2)
        WRITE(LUPRI,'(6(I6,2X,A12))')
     &      (I,PLABEL(I,1),I=1+IOFF,NPLAB(0))
        IF(IPRINT.GE.7) THEN
          CALL PRIVEC('IPLAB ',IPLAB(1,1),NAORB)
        ENDIF
      ENDIF
C
C     **********************************************
C     ***** L A B E L S  for   S O - B A S I S *****
C     **********************************************
C
      NLAB = 0
      DO IC = 1,2
        DO IRP = 0,MAXREP
          DO 100 I = 1,NBBAS(IRP,IC)
            IU    = ICOS(IRP+1,IC) + I
            IS    = IBBAS(IRP,IC)  + I
            ICENT = IPCEN(IU)
            ITYP  = IPTYP(IU)
            JLAB  = IPACK(IC,ICENT,IRP,ITYP)
            DO ILAB = 1,NLAB
            IF(JLAB.EQ.IATTR(ILAB,2)) THEN
              IPLAB(IS,2) = ILAB
              GOTO 100
            ENDIF
            ENDDO
            NLAB = NLAB + 1
            IATTR(NLAB,2)  = JLAB
            PLABEL(NLAB,2) = CLS(IC)//' '//REP(IRP)
     &              //NAMN(ICENT)(1:3)//GTOTYP(ITYP)
            IPLAB(IS,2)    = NLAB
 100      CONTINUE
        ENDDO
      ENDDO      
C
C     Print section
C     =============
C
      IF((IPRINT.GE.1).AND.(NBSYM.GT.1)) THEN
        CALL HEADER('GETLAB: SO-labels',-1)
        WRITE(LUPRI,'(3X,A,I5)') '* Large components:',NPLAB(1)
        WRITE(LUPRI,'(6(I6,2X,A12))') (I,PLABEL(I,2),I=1,NPLAB(1))
        IOFF = NPLAB(1)
        WRITE(LUPRI,'(3X,A,I5)') '* Small components:',NPLAB(2)
        WRITE(LUPRI,'(6(I6,2X,A12))')
     &      (I,PLABEL(I,2),I=1+IOFF,NPLAB(0))
        IF(IPRINT.GE.7) THEN
           CALL PRIVEC('IPLAB ',IPLAB(1,2),NAORB)
        ENDIF
      ENDIF
C
C     Use label information to list SO-basis:
C
      NDIM = 2*NPLAB(0)
      allocate(kp(ndim))
      CALL LABLST(KP,iprint)
      deallocate(kp)
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GETLAB')
      RETURN
C
  999 CONTINUE
      WRITE (LUPRI,'(/A,I6/A,I6)')
     *   ' >>> PLABEL error, '//
     *   'no. of labels',NLAB,
     *   '                  '//
     *   ' current maximum number                  ',MAXLAB
      CALL QUIT('GETLAB:Too many labels !!!')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qnphase */
      SUBROUTINE QNPHASE(IM,JSYMOP,ONEINT,IDOINT)
C***********************************************************************
C
C     Quaternion phase insertion in lower triangular matrix
C       Phase: (e_K1*)(e_IM)(e_K2)
C     + overall factor from IDOINT (e.g. from beta factor)
C
C     Theory: T.Saue &  H.J.Aa.Jensen
C     Written by T.Saue October 1994
C     Last revision: Nov 7 2001 - hjaaj
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dgroup.h"
#include "dcbbas.h"
      LOGICAL DIAG,SAME
      DIMENSION ONEINT(NNBBASX), IDOINT(2,2)
      PARAMETER (IL = 2, IR = 1)
C
#include "ibtfun.h"
C
      I1 = 0
      DO 10 IRP1 = 0,NBSYM-1
        DO 20 IC1 = 1,2
          NB1 = NBBAS(IRP1,IC1)
          K1  = JQBAS(IRP1,IC1)
          I2 = 0
          DO 30 IRP2 = 0,IRP1
            SAME = IRP1.EQ.IRP2
            ICMX = 2
            IF(SAME) ICMX = IC1
            DO 40 IC2 = 1,ICMX
              NB2 = NBBAS(IRP2,IC2)
              IF((IDOINT(IC1,IC2) .NE. 0) .AND.
     &           (IBTXOR(IRP1,IRP2).EQ.JSYMOP)) THEN
                K2   = JQBAS(IRP2,IC2)
                KFAC = IQPHASE(K1,K2,IM)*IQPH(K1,IL)*IQPH(K2,IR)
     &                 *IDOINT(IC1,IC2)
                IF(KFAC.EQ.-1) THEN
                  DIAG = SAME.AND.(IC1.EQ.IC2)
                  DO 50 J1 = 1,NB1
                    JMX = NB2
                    IF(DIAG) JMX = J1
                    JOFF = (I1+J1)*(I1+J1-1)/2 + I2
                    DO 60 J2 = 1,JMX
                      ONEINT(JOFF+J2) = -ONEINT(JOFF+J2)
   60               CONTINUE
   50             CONTINUE
                ENDIF
              ENDIF
              I2 = I2 + NB2
   40       CONTINUE
   30     CONTINUE
          I1   = I1 + NB1
   20   CONTINUE
   10 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck q2phase */
      SUBROUTINE Q2PHASE(TYP,IM,JFAC,AMAT)
C***********************************************************************
C
C     Quaternion phase insertion in square matrix
C       TYP = F:  (e_K1)(e_IM)(e_K2*)   as Fock matrix
C       TYP = D:  (e_K1*)(e_IM)(e_K2)   as density matrix
C
C     Written by T.Saue, January 1995
C     Last revision: June 25 1996 - tsaue
C
C     radovan: in case you ask yourself where the difference is
C              between q2phase and q2bphase
C              q2phase  (this routine) is for Hermit sorted matrices
C              q2bphase (next routine) is for DIRAC  sorted matrices
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
C
#include "symmet.h"
#include "dcbbas.h"
#include "dgroup.h"
      CHARACTER TYP*1
      DIMENSION AMAT(*)
C
      IF(TYP.NE.'F'.AND.TYP.NE.'D') THEN
        CALL QUIT('Q2PHASE: Unknown TYP !')
      ENDIF
      DO 10 IRP1 = 0,NBSYM-1
        ISY1 = IRP1 + 1
        DO 20 IC1 = 1,2
          NB1 = NBBAS(IRP1,IC1)
          IF(NB1.EQ.0) GOTO 20
          K1  = JQBAS(IRP1,IC1)
          DO 30 IRP2 = 0,NBSYM-1
            ISY2 = IRP2 + 1
            DO 40 IC2 = 1,2
              NB2 = NBBAS(IRP2,IC2)
              IF(NB2.EQ.0) GOTO 40
              K2   = JQBAS(IRP2,IC2)
              IF    (TYP.EQ.'F') THEN
                KFAC = JFAC*IQPHASE(K1,K2,IM)*IQPH(K1,1)*IQPH(K2,2)
              ELSEIF(TYP.EQ.'D') THEN
                IQ   = IQMULT(K1,K2,IM)
                KFAC = JFAC*IQPHASE(K1,K2,IQ)*IQPH(K1,2)*IQPH(K2,1)
              ENDIF
              IF(KFAC.EQ.-1) THEN
                IOFF = I2COSX(ISY1,ISY2,IC1,IC2)
                DO 50 ICOL = 1,NB2
                  CALL DSCAL(NB1,DM1,AMAT(IOFF+1),1)
                  IOFF = IOFF + NTBAS(0)
   50           CONTINUE
              ENDIF
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck q2bphase */
      SUBROUTINE Q2BPHASE(TYP,IM,JFAC,AMAT)
C***********************************************************************
C
C     Quaternion phase insertion in square matrix
C       TYP = F:  (e_K1)(e_IM)(e_K2*)   as Fock matrix
C       TYP = D:  (e_K1*)(e_IM)(e_K2)   as density matrix
C
C     Written by T.Saue, January 1995
C     Last revision: June 25 1996 - tsaue
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
C
#include "symmet.h"
#include "dcbbas.h"
#include "dgroup.h"
      CHARACTER TYP*1
      DIMENSION AMAT(*)
C
      IF(TYP.NE.'F'.AND.TYP.NE.'D') THEN
        CALL QUIT('Q2BPHASE: Unknown TYP !')
      ENDIF
      DO 10 IRP1 = 0,NBSYM-1
        ISY1 = IRP1 + 1
        DO 20 IC1 = 1,2
          NB1 = NBBAS(IRP1,IC1)
          IF(NB1.EQ.0) GOTO 20
          K1  = JQBAS(IRP1,IC1)
          DO 30 IRP2 = 0,NBSYM-1
            ISY2 = IRP2 + 1
            DO 40 IC2 = 1,2
              NB2 = NBBAS(IRP2,IC2)
              IF(NB2.EQ.0) GOTO 40
              K2   = JQBAS(IRP2,IC2)
              IF    (TYP.EQ.'F') THEN
                KFAC = JFAC*IQPHASE(K1,K2,IM)*IQPH(K1,1)*IQPH(K2,2)
              ELSEIF(TYP.EQ.'D') THEN
                IQ   = IQMULT(K1,K2,IM)
                KFAC = JFAC*IQPHASE(K1,K2,IQ)*IQPH(K1,2)*IQPH(K2,1)
              ENDIF
              IF(KFAC.EQ.-1) THEN
CTROND                IOFF = I2COSX(ISY1,ISY2,IC1,IC2)
                IOFF = I2BBASX(IRP1,IRP2,IC1,IC2)
                DO 50 ICOL = 1,NB2
                  CALL DSCAL(NB1,DM1,AMAT(IOFF+1),1)
                  IOFF = IOFF + NTBAS(0)
   50           CONTINUE
              ENDIF
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      END

      subroutine q2bphase_nonquadmat_LS_L(typ,im,jfac,amat)
C***********************************************************************
C
C     Quaternion phase insertion in non-quadratic matrices amat(LS,L)
C       TYP = F:  (e_K1)(e_IM)(e_K2*)   as Fock matrix
C       TYP = D:  (e_K1*)(e_IM)(e_K2)   as density matrix
C
C     based on q2bphase; adaption for non-quadratic matrices by S. Knecht
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
C
#include "symmet.h"
#include "dcbbas.h"
#include "dgroup.h"
      CHARACTER TYP*1
      DIMENSION AMAT(*)
C
      do irp1 = 0, nbsym-1
        do ic1 = 1, 2

          nb1  = nbbas(irp1,ic1)
          k1   = jqbas(irp1,ic1)

          if(nb1.gt.0) then

            do irp2 = 0, nbsym-1

              nb2 = nbbas(irp2,1)
              k2  = jqbas(irp2,1)

              if(nb2.gt.0)then

                if     (typ.eq.'F')then
                  kfac = jfac*iqphase(k1,k2,im)*iqph(k1,1)*iqph(k2,2)
                else if(typ.eq.'D')then
                  iq   = iqmult(k1,k2,im)
                  kfac = jfac*iqphase(k1,k2,iq)*iqph(k1,2)*iqph(k2,1)
                end if

                if(kfac.eq.-1)then
                  ioff = i2bbasx(irp1,irp2,ic1,1)
                  do icol = 1, nb2
                    call dscal(nb1,dm1,amat(ioff+1),1)
                    ioff = ioff + ntbas(0)
                  end do
                endif
              end if ! NB2 > 0
            end do ! loop over boson irreps in columns
          end if ! NB1 > 0
        end do ! loop over L,S
      end do ! loop over boson irreps in rows
C
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rpgen */
      SUBROUTINE RPGEN(IREPA)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
      DIMENSION IREPA(*)
      DO 10 ISYM = 1,NBSYM
        IREP = ISYM - 1
        DO 20 IC = 1,2
          CALL ICOPY(NCOS(ISYM,IC),IREP,0,IREPA(ICOS(ISYM,IC)+1),1)
 20     CONTINUE
 10   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtqoso */
      SUBROUTINE MTQOSO(AQO,ASO,IPRINT)
C***********************************************************************
C
C     Matrix Transform: QO-basis to SO-basis
C     Written by Trond Saue
C     Last revision: Jan 7 1998 - jth LINUX
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "dcblab.h"
#include "dgroup.h"
#include "dcbbas.h"
C  Used from COMMON blocks:
C     MXORB :  MXCORB
C     LABELS:  IPLAB,IATTR
C     DGROUP:  NZ,IQDEF
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AQO(NBAST,NBAST,NZ),ASO(NBAST,NBAST,4)
#else
      DIMENSION AQO(NTBAS(0),NTBAS(0),NZ),ASO(NTBAS(0),NTBAS(0),4)
#endif
#include "dcbibt.h"
      CALL DZERO(ASO,N2BBASX*4)
      DO K = 1,NZ
        IK = IPQTOQ(K,0)
        DO J = 1,NTBAS(0)
          IQJ  = JQBAS(IGET(IATTR(IPLAB(J,1),1)),
     &                 JGET(IATTR(IPLAB(J,1),1)))
          DO I = 1,NTBAS(0)
             IQI  = JQBAS(IGET(IATTR(IPLAB(I,1),1)),
     &                    JGET(IATTR(IPLAB(I,1),1)))
             IM   = IQMULT(IQI,IQJ,IK)
             ASO(I,J,IM) = (IQPHASE(IQI,IQJ,IK)*
     &          IQPH(IQI,1)*IQPH(IQJ,2))*AQO(I,J,K)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtsoqo */
      SUBROUTINE MTSOQO(ASO,AQO,IPRINT)
C***********************************************************************
C
C     Matrix transform: SO-basis to QO-basis
C     Written by Trond Saue
C     Last revision: Jan 7 1998 -jth LINUX
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "dcblab.h"
#include "dgroup.h"
#include "dcbbas.h"
C  Used from COMMON blocks:
C     MXORB :  MXCORB
C     LABELS:  IPLAB,IATTR
C     DGROUP:  NZ,IQDEF
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AQO(NBAST,NBAST,NZ),ASO(NBAST,NBAST,4)
#else
      DIMENSION AQO(NTBAS(0),NTBAS(0),NZ),ASO(NTBAS(0),NTBAS(0),4)
#endif
#include "dcbibt.h"
      CALL DZERO(AQO,N2BBASX*NZ)
      DO M = 1,NZ
        DO J = 1,NTBAS(0)
          IQJ  = JQBAS(IGET(IATTR(IPLAB(J,1),1)),
     &                 JGET(IATTR(IPLAB(J,1),1)))
          DO I = 1,NTBAS(0)
             IQI  = JQBAS(IGET(IATTR(IPLAB(I,1),1)),
     &                    JGET(IATTR(IPLAB(I,1),1)))
             IM    = IPQTOQ(M,0)
             IK    = IQMULT(IQI,IQJ,IM)
             AQO(I,J,M) = (IQPHASE(IQI,IQJ,IK)*IQPH(I,2)*IQPH(J,1))
     &                    *ASO(I,J,IK)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtqoso */
      SUBROUTINE VTQOSO(AQO,ASO,IPRINT)
C***********************************************************************
C
C     Vector Transform: QO-basis to SO-basis
C     Last revision: Jan 7 1998 - jth LINUX
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "dcblab.h"
#include "dgroup.h"
#include "dcbbas.h"
C  Used from COMMON blocks:
C     MXORB :  MXCORB
C     LABELS:  IPLAB,IATTR
C     DGROUP:  NZ,IQDEF
C     DCBBAS:  NTBAS(0),N2BBASX
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AQO(NBAST,NBAST,NZ),ASO(NBAST,NBAST,4)
#else
      DIMENSION AQO(NTBAS(0),NTBAS(0),NZ),ASO(NTBAS(0),NTBAS(0),4)
#endif
#include "dcbibt.h"
      CALL DZERO(ASO,N2BBASX*4)
      DO K = 1,NZ
        IK = IPQTOQ(K,0)
        DO J = 1,NTBAS(0)
          IQJ  = JQBAS(IGET(IATTR(IPLAB(J,1),1)),
     &                 JGET(IATTR(IPLAB(J,1),1)))
          DO I = 1,NTBAS(0)
             IQI  = JQBAS(IGET(IATTR(IPLAB(I,1),1)),
     &                    JGET(IATTR(IPLAB(I,1),1)))
             IM   = IQMULT(IQI,IQJ,IK)
             ASO(I,J,IM) = (IQPHASE(IQI,IQJ,IK)*
     &          IQPH(IQI,1)*IQPH(IQJ,2))*AQO(I,J,K)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtsoqo */
      SUBROUTINE VTSOQO(ASO,AQO,IPRINT)
C***********************************************************************
C
C     Vector transform: SO-basis to QO-basis
C     Written by Trond Saue
C     Last revision: Jan 7 1998 - jth LINUX
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "dcblab.h"
#include "dgroup.h"
#include "dcbbas.h"
C  Used from COMMON blocks:
C     MXORB :  MXCORB
C     LABELS:  IPLAB,IATTR
C     DGROUP:  NZ,IQDEF
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AQO(NBAST,NBAST,NZ),ASO(NBAST,NBAST,4)
#else
      DIMENSION AQO(NTBAS(0),NTBAS(0),NZ),ASO(NTBAS(0),NTBAS(0),4)
#endif
#include "dcbibt.h"
      CALL DZERO(AQO,N2BBASX*NZ)
      DO M = 1,NZ
        DO J = 1,NTBAS(0)
          IQJ  = JQBAS(IGET(IATTR(IPLAB(J,1),1)),
     &                 JGET(IATTR(IPLAB(J,1),1)))
          DO I = 1,NTBAS(0)
             IQI  = JQBAS(IGET(IATTR(IPLAB(I,1),1)),
     &                    JGET(IATTR(IPLAB(I,1),1)))
             IM    = IPQTOQ(M,0)
             IK    = IQMULT(IQI,IQJ,IM)
             AQO(I,J,M) = (IQPHASE(IQI,IQJ,IK)*IQPH(I,2)*IQPH(J,1))
     &                   *ASO(I,J,IK)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtsoao */
      SUBROUTINE MTSOAO(ASO,AAO,NBAST,IREPDM,IPRINT)
C***********************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine MTSOAO',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
      CALL DZERO(AAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            AN     = D1/FMULT(MULA)
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = AN*
     &              PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
CTOND               FACA = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     BN     = D1/FMULT(MULB)
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = BN*PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        AAO(INDA,INDB) = AAO(INDA,INDB)
     &                                 + FACA*FACB*ASO(I,J)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         CALL HEADER('OUTPUT: Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
      RETURN
      END

C***********************************************************************
      SUBROUTINE MTSOAO_nonquadmatL(ASO,AAO,NBAST1,NBAST2,IREPDM,IPRINT)

C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST1,NBAST2), AAO(NBAST1,NBAST2)
      DIMENSION IR(2)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF(IPRINT .GE. 10) 
     &  CALL HEADER('Subroutine MTSOAO_nonquadmatL',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
      CALL DZERO(AAO,NBAST1*NBAST2)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            AN     = D1/FMULT(MULA)
            print *, 'NA, IA, KSTRT(IA), NHKTA, KHKTA',
     &                NA, IA, KSTRT(IA), NHKTA, KHKTA
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = AN*
     &              PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB_tmp   = 1
               icount_S_ao = 0
               print *,'new INDA --> II',INDA,II
               DO 400 IREPB = 0, MAXREP
                DO 450 JJ = 1,1
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  print *, 'IREPA,IREPB',
     &                      IREPA,IREPB
                  print *, 'component is',
     &            LCLASS(IBTAND(IBTSHR(IPIND(ISTRB_tmp),16),65535))
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB_tmp,ISTRB_tmp + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     BN     = D1/FMULT(MULB)
                     print *, 'class of a, b = ',LCLASS(ia),lclass(ib)
                     print *, 'NB, IB, KSTRT(IB), NHKTB, KHKTB',
     &                         NB, IB, KSTRT(IB), NHKTB, KHKTB
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = BN*PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        AAO(INDA,INDB) = AAO(INDA,INDB)
     &                                 + FACA*FACB*ASO(I,J-icount_S_ao)
!    &                                 + FACA*FACB*ASO(I,J)
                        print *, 
     &                  'ISYMB, ASO(',
     &                   I,J-icount_S_ao,') ==> AAO(',INDA,INDB,')',
     &                   ISYMB, ASO(I,J-icount_S_ao)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB_tmp = ISTRB_tmp + NBB
                  print *, 'NBB',NBB,'added: --> ISTRB_tmp =',ISTRB_tmp
  450             CONTINUE
                  ISTRB       = 1
                  icount_S_ao = 0
                  do irepb_tmp = 0,IREPB
                    do ic2_tmp = 1, 2
                      ISTRB = ISTRB + NCOS(irepb_tmp+1,ic2_tmp)
                      print '(a,i3,a,i3,a,i4,a,i5)', 
     &                          ' extra computing for component:',
     &                          ic2_tmp,' and boson irrep ',irepb_tmp,
     &                          ' yields with',
     &                          NCOS(irepb_tmp+1,ic2_tmp),
     &                          ' basis func. the new offset:',ISTRB
                      if(ic2_tmp.eq.2) icount_S_ao =
     &                icount_S_ao + ncos(irepb_tmp+1,ic2_tmp)
                      print '(a,i4)', 'total # of S-AOs so far:',
     &                                         icount_S_ao
                    end do
                  end do
                  ISTRB_tmp = ISTRB
                  print *, 'reset ISTRB_tmp (L+S) --> ',ISTRB_tmp
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST1,1,NBAST2,NBAST1,NBAST2,-1,LUPRI)
         CALL HEADER('OUTPUT: Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST1,1,NBAST2,NBAST1,NBAST2,-1,LUPRI)
      END IF
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtaoso */
      SUBROUTINE MTAOSO(AAO,ASO,NBAST,IREPDM,IPRINT)
C***********************************************************************
C
C     Matrix transform: AO-basis to SO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"

      CALL QENTER('MTAOSO')

      IF (IPRINT .GE. 10) THEN
        CALL HEADER('Subroutine MTAOSO',-1)
        WRITE(LUPRI,'(2X,A,I2)') 'entering symmetry, IREPDM=',IREPDM
      ENDIF
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            AN     = D1/FMULT(MULA)
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     BN     = D1/FMULT(MULB)
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        ASO(I,J) = ASO(I,J)
     &                                 + FACA*FACB*AAO(INDA,INDB)
          IF (IPRINT.GE.28) THEN 
CMI        ... control print-out to track the AO-SO matrix construction
             WRITE(LUPRI,
     &'(2X,A,I3,A,I3,A,D12.6,A,D10.4,A,D10.4,A,I3,A,I3,A,D12.6)')
     &'ASO(',I,',',J,') =',ASO(I,J),
     &'  from :  ',
     & FACA,'*',FACB,'*AAO(',INDA,',',INDB,')=',AAO(INDA,INDB)
          ENDIF
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE

      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT(MTAOSO):  Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         CALL HEADER('OUTPUT(MTAOSO): Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF

      CALL QEXIT('MTAOSO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dtsoao */
      SUBROUTINE DTSOAO(ASO,AAO,NBAST,IREPDM,IPRINT)
C***********************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine DTSOAO',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
      CALL DZERO(AAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        AAO(INDA,INDB) = AAO(INDA,INDB)
     &                                 + FACA*FACB*ASO(I,J)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         CALL HEADER('OUTPUT:Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dtaoso */
      SUBROUTINE DTAOSO(AAO,ASO,NBAST,IREPDM,IPRINT)
C***********************************************************************
C
C     Matrix transform: AO-basis to SO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine DTAOSO',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            AN     = D1/FMULT(MULA)
            INDA   = KSTRT(IA) + NA - KHKTA
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = AN*
     &              PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     BN     = D1/FMULT(MULB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = BN*PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        ASO(I,J) = ASO(I,J)
     &                                 + FACA*FACB*AAO(INDA,INDB)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         CALL HEADER('OUTPUT:Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ntaoso */
      SUBROUTINE NTAOSO(AAO,ASO,NBAST,IREPM,IPRINT)
C***********************************************************************
C
C     Matrix Transform:  AO-basis to SO-basis
C     This routine is based on HERMIT routine DAOTSO (880418  PRT)
C     Note that the SO-matrix is NOT initialized !!!
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sept 10 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "aosotr.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine MTAOSO',-1)
C
C     Loop over all irreps in molecule
C
      IFXYZ  = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(1) = IREPM
      IR(2) = IBTXOR(IREPM,IFXYZ)
      ISTRA = 1
C      CALL DZERO(ASO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            JKBA   = JTRAN(I)
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            MULA   = ISTBAO(IA)
            AN     = D1/FMULT(MULA)
            DO 300 IKBA = 1,JKBA
               INDA = ITRAN(I,IKBA)
               FACA = CTRAN(I,IKBA)
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     JKBB = JTRAN(J)
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     MULB   = ISTBAO(IB)
                     BN     = D1/FMULT(MULB)
                     DO 600 IKBB = 1,JKBB
                        INDB = ITRAN(J,IKBB)
                        FACB = CTRAN(J,IKBB)
                        ASO(I,J) = ASO(I,J)
     &                                 + FACA*FACB*AAO(INDA,INDB)
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         CALL HEADER('Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtsoao */
      SUBROUTINE VTSOAO(ASO,AAO,NBAST,NVEC,IPRINT)
C***********************************************************************
C
C     Vector transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NVEC), AAO(NBAST,NVEC)
      DIMENSION IR(2)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine VTSOAO',-1)
C
C     Loop over all irreps in molecule
C
      ISTRA = 1
      CALL DZERO(AAO,NBAST*NVEC)
      DO 100 IREPA = 0, MAXREP
        NORBA = NAOS(IREPA+1)
        DO 200 I = ISTRA,ISTRA + NORBA - 1
          IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
          NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
          NHKTA  = NHKT(IA)
          KHKTA  = KHKT(IA)
          MULA   = ISTBAO(IA)
          INDA   = KSTRT(IA) + NA - KHKTA
C
C         Loop over symmetry dependent centers
C
          DO 300 ISYMA = 0, MAXOPR
          IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
            INDA = INDA + KHKTA
            FAC  = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
            CALL DAXPY(NVEC,FAC,ASO(I,1),NBAST,AAO(INDA,1),NBAST)
          ENDIF
  300     CONTINUE
  200   CONTINUE
        ISTRA = ISTRA + NORBA
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Vectors in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NVEC,NBAST,NVEC,-1,LUPRI)
         CALL HEADER('OUTPUT: Vectors in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NVEC,NBAST,NVEC,-1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wtsoao */
      SUBROUTINE WTSOAO(ASO,AAO,IFRP,KNFBAS,NBAST,NVEC,IPRINT)
C***********************************************************************
C
C     Vector transform: SO-basis to AO-basis (in sorted basis !)
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Jan 7 1998 (jth: NFBAS(IFRP,0) now parameter)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
C
C     KNFBAS = NFBAS(IFRP,0)
C     NBAST  = NTBAS(0)
C
      DIMENSION ASO(KNFBAS,NVEC), AAO(NBAST,NVEC)
      DIMENSION IR(2)
C
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine WTSOAO',-1)
C.....Loop over all irreps in molecule
      CALL DZERO(AAO,NTBAS(0)*NVEC)
      DO 100 IREPA = 0, MAXREP
        DO IC = 1,MC
        NORBA = NCOS(IREPA+1,IC)
        IF(JBTOF(IREPA,IC).EQ.IFRP) THEN
          DO 200 I = 1,NORBA
            IU     = ICOS(IREPA+1,IC) + I
            IS     = IBBAS(IREPA,IC) - IBAS(IFRP) + I
            IA     = IBTAND(IBTSHR(IPIND(IU),16),65535)
            NA     = IBTAND(IBTSHR(IPIND(IU), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
C...........Loop over symmetry dependent centers
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
              INDA = INDA + KHKTA
              FAC  = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
              CALL DAXPY(NVEC,FAC,ASO(IS,1),NFBAS(IFRP,0),
     &                            AAO(INDA,1),NTBAS(0))
            ENDIF
  300       CONTINUE
  200     CONTINUE
        ENDIF
        ENDDO
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Vectors in SO basis',-1)
         CALL OUTPUT(ASO,1,NFBAS(IFRP,0),1,NVEC,
     &               NFBAS(IFRP,0),NVEC,-1,LUPRI)
         CALL HEADER('INPUT: Vectors in AO basis',-1)
         CALL OUTPUT(AAO,1,NTBAS(0),1,NVEC,NTBAS(0),NVEC,-1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wtsoao */
      SUBROUTINE WTSOAO_mod(ASO,AAO,IFRP,KNFBAS,NBAST,NVEC,IROW,IPRINT)
C***********************************************************************
C
C     Vector transform: SO-basis to AO-basis (in sorted basis !)
C     Based on WTSOAO, but with row ordering
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
C
C     KNFBAS = NFBAS(IFRP,0)
C     NBAST  = NTBAS(0)
C
      DIMENSION ASO(KNFBAS,NVEC), AAO(NBAST,NVEC),IROW(NBAST)
      DIMENSION IR(2)
C
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine WTSOAO_mod',-1)
C.....Loop over all irreps in molecule
      CALL DZERO(AAO,NTBAS(0)*NVEC)
      DO 100 IREPA = 0, MAXREP
        DO IC = 1,MC
        NORBA = NCOS(IREPA+1,IC)
        IF(JBTOF(IREPA,IC).EQ.IFRP) THEN
          DO 200 I = 1,NORBA
            IU     = ICOS(IREPA+1,IC) + I
            IS     = IBBAS(IREPA,IC) - IBAS(IFRP) + I
            IA     = IBTAND(IBTSHR(IPIND(IU),16),65535)
            NA     = IBTAND(IBTSHR(IPIND(IU), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
C...........Loop over symmetry dependent centers
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
              INDA = INDA + KHKTA
              IA   = IROW(INDA)
              FAC  = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
              CALL DAXPY(NVEC,FAC,ASO(IS,1),NFBAS(IFRP,0),
     &                            AAO(IA,1),NTBAS(0))
            ENDIF
  300       CONTINUE
  200     CONTINUE
        ENDIF
        ENDDO
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Vectors in SO basis',-1)
         CALL OUTPUT(ASO,1,NFBAS(IFRP,0),1,NVEC,
     &               NFBAS(IFRP,0),NVEC,-1,LUPRI)
         CALL HEADER('INPUT: Vectors in AO basis',-1)
         CALL OUTPUT(AAO,1,NTBAS(0),1,NVEC,NTBAS(0),NVEC,-1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wtaoso */
      SUBROUTINE WTAOSO(AAO,ASO,NBAST,NVEC,IPRINT)
C***********************************************************************
C
C     Vector Transform:  AO-basis to SO-basis
C     VTAOSO/WTAOSO are covariant/contarvariant versions
C       VTAOSO is good for vectors and WTAOSO for functions
C     This routine is based on HERMIT routine DAOTSO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sept 10 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NVEC), AAO(NBAST,NVEC)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine WTAOSO',-1)
C
C     Loop over all irreps in molecule
C
      ISTRA = 1
      CALL DZERO(ASO,NBAST*NVEC)
      DO 100 IREPA = 0, MAXREP
        NORBA = NAOS(IREPA+1)
        DO 200 I = ISTRA,ISTRA + NORBA - 1
          IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
          NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
          NHKTA  = NHKT(IA)
          KHKTA  = KHKT(IA)
          MULA   = ISTBAO(IA)
          INDA   = KSTRT(IA) + NA - KHKTA
C
C         Loop over symmetry dependent centers
C
          DO 300 ISYMA = 0, MAXOPR
          IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
            INDA = INDA + KHKTA
            FAC  = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
            CALL DAXPY(NVEC,FAC,AAO(INDA,1),NBAST,ASO(I,1),NBAST)
          ENDIF
  300     CONTINUE
  200   CONTINUE
        ISTRA = ISTRA + NORBA
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Vectors in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NVEC,NBAST,NVEC,-1,LUPRI)
         CALL HEADER('OUTPUT: Vectors in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NVEC,NBAST,NVEC,-1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtaoso */
      SUBROUTINE VTAOSO(AAO,ASO,NBAST,NVEC,IPRINT)
C***********************************************************************
C
C     Vector Transform:  AO-basis to SO-basis
C     This routine is based on HERMIT routine DAOTSO (880418  PRT)
C     Note that the SO-matrix is NOT initialized !!!
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sept 10 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NVEC), AAO(NBAST,NVEC)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine VTAOSO',-1)
C
C     Loop over all irreps in molecule
C
      ISTRA = 1
      CALL DZERO(ASO,NBAST*NVEC)
      DO 100 IREPA = 0, MAXREP
        NORBA = NAOS(IREPA+1)
        DO 200 I = ISTRA,ISTRA + NORBA - 1
          IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
          NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
          NHKTA  = NHKT(IA)
          KHKTA  = KHKT(IA)
          MULA   = ISTBAO(IA)
          AN     = D1/FMULT(MULA)
          INDA   = KSTRT(IA) + NA - KHKTA
C
C         Loop over symmetry dependent centers
C
          DO 300 ISYMA = 0, MAXOPR
          IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
            INDA = INDA + KHKTA
            FAC  = AN*PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
            CALL DAXPY(NVEC,FAC,AAO(INDA,1),NBAST,ASO(I,1),NBAST)
          ENDIF
  300     CONTINUE
  200   CONTINUE
        ISTRA = ISTRA + NORBA
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Vectors in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NVEC,NBAST,NVEC,-1,LUPRI)
         CALL HEADER('OUTPUT: Vectors in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NVEC,NBAST,NVEC,-1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bstobu */
      SUBROUTINE BSTOBU(AMAT,MZ,WORK,LWORK)
C***********************************************************************
C
C     Matrix transform:
C     from class/symmetry sorted basis ("sorted basis" = Dirac sorted)
C     to symmetry/class sorted basis ("unsorted basis" = Hermit sorted)
C
C     Written by T.Saue Mar 2 1995
C     Last revision : May 6 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
      DIMENSION AMAT(*),WORK(LWORK)
      CALL QENTER('BSTOBU')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET2('REAL','BUF',KBUF,N2BBASX,WORK,KFREE,LFREE)
      IZOFF = 0
      DO IZ = 1,MZ
        CALL DCOPY(N2BBASX,AMAT(IZOFF+1),1,WORK(KBUF),1)
        DO JC = 1,2
          DO JREP = 0,MAXREP
            DO J = 0,NBBAS(JREP,JC)-1
              JSOFF =         (IBBAS(JREP  ,JC)+J)*NTBAS(0)
              JUOFF = IZOFF + (ICOS (JREP+1,JC)+J)*NTBAS(0)
              DO IC = 1,2
                DO IREP = 0,MAXREP
                  IUOFF = JUOFF + ICOS (IREP+1,IC)
                  ISOFF = JSOFF + IBBAS(IREP  ,IC)
                  IF(NBBAS(IREP,IC).GT.0)
     &              CALL DCOPY(NBBAS(IREP,IC),WORK(KBUF+ISOFF),1,
     &                                        AMAT(IUOFF+1)   ,1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IZOFF = IZOFF + N2BBASX
      ENDDO
C
C     Memory deallocation
C
      CALL MEMREL('BSTOBU',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT ('BSTOBU')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck butobs */
      SUBROUTINE BUTOBS(AMAT,MZ,WORK,LWORK)
C***********************************************************************
C
C     Matrix transform:
C     from symmetry/class sorted basis ("unsorted basis" = Hermit sorted)
C     to class/symmetry sorted basis ("sorted basis" = Dirac sorted)
C
C     Written by T.Saue Mar 2 1995
C     Last revision : May 6 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
      DIMENSION AMAT(N2BBASX,MZ),WORK(LWORK)
      CALL QENTER('BUTOBS')
#include "memint.h"
C

      CALL MEMGET2('REAL','BUF',KBUF,N2BBASX,WORK,KFREE,LFREE)
      IZOFF = 0
      DO IZ = 1,MZ
        CALL DCOPY(N2BBASX,AMAT(1,IZ),1,WORK(KBUF),1)
        DO JC = 1,2
          DO JREP = 0,MAXREP
            DO J = 0,NBBAS(JREP,JC)-1
              JSOFF = (IBBAS(JREP  ,JC)+J)*NTBAS(0)
              JUOFF = (ICOS (JREP+1,JC)+J)*NTBAS(0)
              DO IC = 1,2
                DO IREP = 0,MAXREP
                  IUOFF = JUOFF + ICOS (IREP+1,IC)
                  ISOFF = JSOFF + IBBAS(IREP  ,IC)
!                  write(*,*)""
!                  write(*,*)"IREP", IREP
!                  write(*,*)"IUOFF",IUOFF
!                  write(*,*)"ISOFF",ISOFF
                  IF(NBBAS(IREP,IC).GT.0) 
     &              CALL DCOPY(NBBAS(IREP,IC),WORK(KBUF+IUOFF),1,
     &                                        AMAT(ISOFF+1,IZ)   ,1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      CALL MEMREL('BUTOBS',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT ('BUTOBS')
      RETURN
C
      END

      SUBROUTINE BSTOBU_no_work(AMAT,MZ)
C***********************************************************************
C
C     Matrix transform:
C     from class/symmetry sorted basis ("sorted basis" = Dirac sorted)
C     to symmetry/class sorted basis ("unsorted basis" = Hermit sorted)
C
C     Written by T.Saue Mar 2 1995
C     Last revision : May 6 1996 - tsaue
C
C***********************************************************************
         use memory_allocator
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
      DIMENSION AMAT(N2BBASX,MZ)
      real(8), allocatable :: buf(:)
      CALL QENTER('BSTOBU_no_work')
C
C     Memory allocation
C
      call alloc(buf, n2bbasx)
      DO IZ = 1,MZ
        CALL DCOPY(N2BBASX,AMAT(1,IZ),1,BUF,1)
        DO JC = 1,2
          DO JREP = 0,MAXREP
            DO J = 0,NBBAS(JREP,JC)-1
              JSOFF = (IBBAS(JREP  ,JC)+J)*NTBAS(0)
              JUOFF = (ICOS (JREP+1,JC)+J)*NTBAS(0)
              DO IC = 1,2
                DO IREP = 0,MAXREP
                  IUOFF = JUOFF + ICOS (IREP+1,IC)
                  ISOFF = JSOFF + IBBAS(IREP  ,IC)
                  IF(NBBAS(IREP,IC).GT.0)
     &              CALL DCOPY(NBBAS(IREP,IC),buf(1+ISOFF),1,
     &                                        AMAT(IUOFF+1,IZ)   ,1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     Memory deallocation
C
      call dealloc(buf)
      CALL QEXIT ('BSTOBU_no_work')
      END

      subroutine bstobu_no_work_nonquadmat_LS_L(amat,mz)
!***********************************************************************
!
!     non-quadratic matrix amat transformation 
!     row dimension   : L+S
!     column dimension: L
!     from class/symmetry sorted basis (  "sorted basis" =  Dirac sorted)
!     to symmetry/class   sorted basis ("unsorted basis" = Hermit sorted)
!
!     adaption of bstobu for non-quadratic matrices: S. Knecht, Feb 2011
!
!***********************************************************************
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
      real(8), intent(inout) :: amat(ntbas(0)*ntbas(1),mz)
      real(8), allocatable   :: bufx(:)
      integer                :: ndim_mat

      call qenter('bstobu_no_work_nonquadmat_LS_L')
      ! component L (1) S (2)
      do i = 1, 2
        do j = 0, maxrep
          print '(a,i3,a,i3,a,i5)', 'jbtof(',j,',',i,') =',jbtof(j,i)
          print '(a,i3,a,i3,a,i5)', 'nbbas(',j,',',i,') =',nbbas(j,i)
          print '(a,i3,a,i3,a,i5)', 'nfbas(',jbtof(j,i),',',i,') =',
     &                               nfbas(jbtof(j,i),i)
        end do
      end do
      do i = 1, nfsym
        print '(a,i3,a,i3,a,i5)', 'nfbas(',i,',',0,') =',
     &                             nfbas(i,0)
      end do
      call flshfo(6)
C
C     Memory allocation
      call alloc(bufx, ntbas(1)*ntbas(0), id='bufx LS-L Dao --> Hao')


!     IREP+1 for icos because loop indices in herrdn.F run from 
!     1 --> maxrep+1 for icos!

      do iz = 1,mz
        call dcopy(ntbas(0)*ntbas(1),amat(1,iz),1,bufx,1)
        call dzero(amat(1,iz),ntbas(0)*ntbas(1))
        do jrep = 0,maxrep
          do j = 0,nbbas(jrep,1)-1
            joff = (ibbas(jrep,1)+j)*ntbas(0)
            print *,'joff ==> column in L part', joff
            do ic = 1, 2
              do irep = 0,maxrep
                print *,'IC, IREP', IC, IREP
                isoff = joff + ibbas(irep  ,ic)
                iuoff = joff + icos (irep+1,ic)
                if(nbbas(irep,ic).gt.0)then

                   call dcopy(nbbas(irep,ic),bufx(1+isoff),1,
     &                                       amat(1+iuoff,iz),1)

                end if

              end do
            end do
          end do
        end do
      end do
C
C     Memory deallocation
      call dealloc(bufx)

      call qexit('bstobu_no_work_nonquadmat_LS_L')
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck butobs */
      SUBROUTINE BUTOBS_no_work(AMAT,MZ)
C***********************************************************************
C
C     Matrix transform:
C     from symmetry/class sorted basis ("unsorted basis" = Hermit sorted)
C     to class/symmetry sorted basis ("sorted basis" = Dirac sorted)
C
C     Written by T.Saue Mar 2 1995
C     Last revision : May 6 1996 - tsaue
C
C***********************************************************************

         use memory_allocator
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
      DIMENSION AMAT(N2BBASX,MZ)
      real(8), allocatable             :: buf(:)
      CALL QENTER('BUTOBS')
C
      call alloc(buf, n2bbasx, id='buf in butobs')
      DO IZ = 1,MZ
        CALL DCOPY(N2BBASX,AMAT(1,IZ),1,buf,1)
        DO JC = 1,2
          DO JREP = 0,MAXREP
            DO J = 0,NBBAS(JREP,JC)-1
              JSOFF = (IBBAS(JREP  ,JC)+J)*NTBAS(0)
              JUOFF = (ICOS (JREP+1,JC)+J)*NTBAS(0)
              DO IC = 1,2
                DO IREP = 0,MAXREP
                  IUOFF = JUOFF + ICOS (IREP+1,IC)
                  ISOFF = JSOFF + IBBAS(IREP  ,IC)
                  if(NBBAS(IREP,IC).GT.0)then
                    CALL DCOPY(NBBAS(IREP,IC),buf(1+IUOFF),1,
     &                                        AMAT(ISOFF+1,IZ)   ,1)
                  end if
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      call dealloc(buf)
      CALL QEXIT ('BUTOBS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtbubs */
      SUBROUTINE VTBUBS(VBU,VBS,NVEC)
C***********************************************************************
C
C     Vector transform:  from unsorted basis to sorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : Jan 7 1998 - jth LINUX
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION VBU(NBAST,NVEC),VBS(NBAST,NVEC)
#else
      DIMENSION VBU(NTBAS(0),NVEC),VBS(NTBAS(0),NVEC)
#endif
      CALL QENTER('VTBUBS')
C
      DO J = 1,NVEC
        DO IC = 1,2
          DO IREP = 0,MAXREP
            IUOFF = ICOS (IREP+1,IC)
            ISOFF = IBBAS(IREP  ,IC)
            IF(NBBAS(IREP,IC).GT.0) 
     &        CALL DCOPY(NBBAS(IREP,IC),VBU(IUOFF+1,J),1,
     &                                  VBS(ISOFF+1,J),1)
          ENDDO
        ENDDO
      ENDDO
C
      CALL QEXIT ('VTBUBS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtbsbu */
      SUBROUTINE VTBSBU(VBS,VBU,NVEC)
C***********************************************************************
C
C     Vector transform:  from sorted basis to unsorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : Jan 7 1998 - jth LINUX
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION VBS(NBAST,NVEC),VBU(NBAST,NVEC)
#else
      DIMENSION VBS(NTBAS(0),NVEC),VBU(NTBAS(0),NVEC)
#endif
      CALL QENTER('VTBSBU')
      DO J = 1,NVEC
        DO IC = 1,2
          DO IREP = 0,MAXREP
            IUOFF = ICOS (IREP+1,IC)
            ISOFF = IBBAS(IREP  ,IC)
            IF(NBBAS(IREP,IC).GT.0)
     &        CALL DCOPY(NBBAS(IREP,IC),VBS(ISOFF+1,J),1,
     &                                  VBU(IUOFF+1,J),1)
          ENDDO
        ENDDO
      ENDDO
      CALL QEXIT ('VTBSBU')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtbsbu */
      SUBROUTINE MTBSBU(AMS,AMU)
C***********************************************************************
C
C     Matrix transform:  from sorted basis to unsorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : July 17 1997 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
      DIMENSION AMS(*),AMU(*)
      CALL QENTER('MTBSBU')
C
      DO JC = 1,2
        DO JREP = 0,MAXREP
          DO J = 0,NBBAS(JREP,JC)-1
            JSOFF = (IBBAS(JREP  ,JC)+J)*NTBAS(0)
            JUOFF = (ICOS (JREP+1,JC)+J)*NTBAS(0)
            DO IC = 1,2
              DO IREP = 0,MAXREP
                IUOFF = JUOFF + ICOS (IREP+1,IC)
                ISOFF = JSOFF + IBBAS(IREP  ,IC)
                IF(NBBAS(IREP,IC).GT.0)
     &            CALL DCOPY(NBBAS(IREP,IC),AMS(ISOFF+1),1,
     &                                      AMU(IUOFF+1),1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL QEXIT ('MTBSBU')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtbubs */
      SUBROUTINE MTBUBS(AMU,AMS)
C***********************************************************************
C
C     Matrix transform:  from unsorted basis to sorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : July 17 1997 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
      DIMENSION AMU(*),AMS(*)
      CALL QENTER('MTBUBS')
      DO JC = 1,2
        DO JREP = 0,MAXREP
          DO J = 0,NBBAS(JREP,JC)-1
            JSOFF = (IBBAS(JREP  ,JC)+J)*NTBAS(0)
            JUOFF = (ICOS (JREP+1,JC)+J)*NTBAS(0)
            DO IC = 1,2
              DO IREP = 0,MAXREP
                IUOFF = JUOFF + ICOS (IREP+1,IC)
                ISOFF = JSOFF + IBBAS(IREP  ,IC)
                IF(NBBAS(IREP,IC).GT.0)
     &            CALL DCOPY(NBBAS(IREP,IC),AMU(IUOFF+1),1,
     &                                      AMS(ISOFF+1),1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL QEXIT ('MTBUBS')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sphgen */
      SUBROUTINE SPHGEN(MORDER,MINTEG,CSCADAPT,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate coefficients for transforming from Cartesian to
C     spherical components
C
C     This is just a driver exploiting all the goodies in SPHCOM of
C     the HERMIT package
C
C     Written by T.Saue Sep 19 1995
C     Last revision: Sep 19 1995
C
CMI   MORDER=1, MINTEG=2 ... previous settings
CMI   CSCADAPT ... if true, adapt CSP matrix for C-S-C transformation
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0)
#include "maxaqn.h"
C
C Used from COMMON blocks
C  SPHCOM: CSP,ISPADR
C
#include "dcbsph.h"

      DIMENSION WORK(LWORK)
      LOGICAL CSCADAPT
C
      CALL QENTER('SPHGEN')
#include "memint.h"

      NCSPT = (MXQN*(MXQN + 1)/2)*(2*MXQN - 1)
      CALL DZERO(CSP,NCSP)
C
      CALL MEMGET('REAL',KCSPT,NCSPT,WORK,KFREE,LFREE)
C
C     MORDER eq 0 : M_l order 0, +1, -1, +2, -2, ..., +LVAL, -LVAL
C     MORDER ne 0 : M_l order -LVAL, ..., -1, 0, 1, ..., +LVAL
C
C     MINTEG eq 1 : the spherical components are normalized such
C                   that all coefficients are integers.
CMI   MORDER = 1
CMI    ... only this works ...
CMI   MINTEG = 2
CMI ... try for CSC tra..
CMI   MINTEG = 1
      IOFF = 1
      DO 10 LVAL = 0, MXQN - 1
        NLM            = 2*LVAL + 1
        NXYZ           = (LVAL+1)*(LVAL+2)/2
        ISPADR(LVAL+1) = IOFF
        IOFF           = IOFF + NXYZ*NXYZ
        IF (LVAL .EQ. 0) THEN
          CSP(1)       = D1
        ELSE IF (LVAL .EQ. 1) THEN
          CSP (2)      = D1
          CSP (6)      = D1
          CSP(10)      = D1
        ELSE
          CALL SPHCOM(LVAL,WORK(KCSPT),NLM,NXYZ,MORDER,MINTEG,
     &                WORK(KFREE),LFREE,IPRINT)
          CALL MTRSP(NXYZ,NLM,WORK(KCSPT),NXYZ,CSP(ISPADR(LVAL+1)),NLM)
        END IF
 10   CONTINUE

C
C     Generate backtransformation as well
C
      DO L = 1,MXQN
        IOFF = ISPADR(L)
        NLM  = 2*L- 1
        NXYZ = (L*(L+1))/2
        CALL MTRSP(NLM,NXYZ,CSP(IOFF),NLM,SPC(IOFF),NXYZ)
      ENDDO

CMI ... adapt the SPC matrixes for S-C correct transformations...
      IF (CSCADAPT) THEN
       CALL MEMGET('REAL',KSCCS,NCSPT,WORK,KFREE,LFREE)
       CALL MEMGET('REAL',KSCCS1,NCSPT,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KIPIV,NCSPT,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KINDXR,NCSPT,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KINDXC,NCSPT,WORK,KFREE,LFREE)

CMI     ... start from d-functions...
        DO L = 3,MXQN
          IOFF = ISPADR(L)
          NLM  = 2*L- 1
          NXYZ = (L*(L+1))/2
C        ... do the multiplication CSP*SPC = WORK(KSCCS)
          CALL DZERO(WORK(KSCCS),NCSPT)
          CALL DZERO(WORK(KSCCS1),NCSPT)
          CALL IZERO(WORK(KIPIV),NCSPT)
          CALL IZERO(WORK(KINDXR),NCSPT)
          CALL IZERO(WORK(KINDXC),NCSPT)

          CALL DGEMM('N','N',NLM,NLM,NXYZ,D1,CSP(IOFF),NLM,
     &          SPC(IOFF),NXYZ,D0,WORK(KSCCS),NLM)

          IF (IPRINT.GE.4) THEN
            WRITE (LUPRI,'(/1X,A,I5)')
     &      'CSP*SPC matrix for L=',L-1
            WRITE (LUPRI,'(1X,A,2I5)') ' size NLM,NLM',NLM,NLM
            CALL OUTPUT(WORK(KSCCS),1,NLM,1,NLM,NLM,NLM,-1,LUPRI)
          ENDIF
C        ... do the matrix inverse, WORK(KSCCS)
          CALL MATRIX_INVERSION(WORK(KSCCS),NLM,WORK(KIPIV),
     &                          WORK(KINDXR),WORK(KINDXC))
          IF (IPRINT.GE.4) THEN
            WRITE (LUPRI,'(1X,A,I5)')
     &      'CSP*SPC inversion matrix for L=',L-1
            WRITE (LUPRI,'(1X,A,2I5)') ' size NLM,NLM',NLM,NLM
            CALL OUTPUT(WORK(KSCCS),1,NLM,1,NLM,NLM,NLM,-1,LUPRI)
          ENDIF
C        ... do the matrix multiplication, WORK(KSCCS)*CSP = CSP
           NLX = NLM*NXYZ
          CALL DCOPY(NLX,CSP(IOFF),1,WORK(KSCCS1),1) 
          CALL DGEMM('N','N',NLM,NXYZ,NLM,D1,WORK(KSCCS),NLM,
     &          WORK(KSCCS1),NLM,D0,CSP(IOFF),NLM)
          IF (IPRINT.GE.4) THEN
            WRITE (LUPRI,'(1X,A,I5)')
     &      'inv(CSP*SPC)*CSP for L=',L-1
            WRITE (LUPRI,'(1X,A,2I5)') ' size NLM, NXYZ:',NLM,NXYZ
            CALL OUTPUT(CSP(IOFF),1,NLM,1,NXYZ,NLM,NXYZ,-1,LUPRI)
          ENDIF
        ENDDO
      ENDIF

C
C     Print section
C
      IF(IPRINT.GE.4) THEN
        CALL HEADER('SPHGEN: Cartesian transformation matrix:',-1)
        DO L = 1,MXQN
          LVAL = L-1
          NLM            = 2*LVAL + 1
          NXYZ           = (LVAL+1)*(LVAL+2)/2
          WRITE (LUPRI,'(1X,A,I5)')
     &     ' Coefficients for angular quantum number ',LVAL
          CALL OUTPUT(CSP(ISPADR(L)),1,NLM,1,NXYZ,NLM,NXYZ,-1,LUPRI)
        ENDDO
        CALL HEADER('SPHGEN: Spherical transformation matrix:',-1)
        DO L = 1,MXQN
          LVAL = L-1
          NLM            = 2*LVAL + 1
          NXYZ           = (LVAL+1)*(LVAL+2)/2
          WRITE (LUPRI,'(1X,A,I5)')
     &     ' Coefficients for angular quantum number ',LVAL
          CALL OUTPUT(SPC(ISPADR(L)),1,NXYZ,1,NLM,NXYZ,NLM,-1,LUPRI)
        ENDDO
      ENDIF
      CALL MEMREL('SPHGEN.SPHGE1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SPHGEN')
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtaocs */
      SUBROUTINE VTAOCS(AOV,LDA,CSV,LDC,NVEC,LCOMP,IPRINT)
C***********************************************************************
C
C     Transform NVEC vectors from Cartesian AO basis to Spherical
C     Components
C
C     Options:
C     LCOMP = 0   Transform all components
C     LCOMP = 1   Transform only large components
C     LCOMP = 2   Transform only small components
C
C     NSPH - number of spherical components
C
C       LDC <= LDA
C
C     Written by T.Saue Sep 19 1995
C     Last revision Sep 19 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
C Used from COMMON blocks:
C    SHELLS: KMAX
C    SPHCOM:
C    NUCLEI: NLARGE,NSMALL,NBASIS
#include "nuclei.h"
#include "shells.h"
#include "dcbsph.h"
      DIMENSION AOV(LDA,NVEC),CSV(LDC,NVEC)

      IF(IPRINT.GE.5) THEN
        CALL HEADER('VTAOCS: Input vectors:',-1)
        WRITE(LUPRI,'(2X,A,2I4)') 'NROWS, NCOLS',LDA,NVEC
        CALL OUTPUT(AOV,1,LDA,1,NVEC,LDA,NVEC,-1,LUPRI)
      ENDIF
C
C     Set limits
C
      NSPH = NSPHCM(LCOMP)
      IF    (LCOMP.EQ.0) THEN
C       All components
        ISH1 = 1
        ISH2 = KMAX
      ELSEIF(LCOMP.EQ.1) THEN
C       Large components
        ISH1 = 1
        ISH2 = NLRGSH
      ELSEIF(LCOMP.EQ.2) THEN
C       Small components
        ISH1 = NLRGSH+1
        ISH2 = KMAX
      ELSE
        WRITE(LUPRI,'(A,I5)') '* Unknown LCOMP = ',LCOMP
        CALL QUIT('VTAOCS: Unknown component LCOMP !!!')
      ENDIF
C
C     Transform to cartesian components
C
      IOFF = 1
      JOFF = 1
      DO ISHELL = ISH1,ISH2
        L     = NHKT(ISHELL)
        NLM   = 2*(L-1) + 1
        NXYZ  = L*(L+1)/2
        NDEG   = NUCDEG(NCENT(ISHELL))
        DO IDEG = 1, NDEG
          CALL DGEMM('N','N',NLM,NVEC,NXYZ,D1,CSP(ISPADR(L)),NLM,
     &               AOV(IOFF,1),LDA,D0,CSV(JOFF,1),LDC)
          IOFF = IOFF + NXYZ
          JOFF = JOFF + NLM
        ENDDO
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
CMI     CALL HEADER('VTAOCS: Transformed vectors(full size):',-1)
CMI     CALL OUTPUT(CSV,1,LDC,1,NVEC,LDC,NVEC,-1,LUPRI)

        CALL HEADER('VTAOCS: Transformed vectors:',-1)
        WRITE(LUPRI,'(2X,A,2I4)') 'NROWS, NCOLS',NSPH,NVEC
        CALL OUTPUT(CSV,1,NSPH,1,NVEC,LDC,NVEC,-1,LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtaocs */
      SUBROUTINE MTAOCS(AOM,LDA,CSM,LDC,LCOMP,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Transform NVEC vectors from Cartesian AO basis to Spherical
C     Components
C
C     Options:
C     LCOMP = 0   Transform all components
C     LCOMP = 1   Transform only large components
C     LCOMP = 2   Transform only small components
C
C     NSPH - total number of spherical components
C
C     Written by T.Saue Sep 19 1995
C     Last revision Sep 19 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include "mxcent.h"
C
#include "nuclei.h"
      DIMENSION AOM(LDA,*),CSM(LDC,*),WORK(LWORK)
      CALL QENTER('MTAOCS')

      IF(IPRINT.GE.5) THEN
        CALL HEADER('MTAOCS: Input AO (C) matrix:',-1)
        WRITE(LUPRI,'(2X,A,2I4)') 'NROWS, NCOLS',LDA,LDA
        CALL OUTPUT(AOM,1,LDA,1,LDA,LDA,LDA,-1,LUPRI)
      ENDIF
C
C     Set limits
C
      NSPH = NSPHCM(LCOMP)
      IF    (LCOMP.EQ.0) THEN
C       All components
        NVEC = NBASIS
      ELSEIF(LCOMP.EQ.1) THEN
C       Large components
        NVEC = NLARGE
      ELSEIF(LCOMP.EQ.2) THEN
C       Small components
        NVEC = NSMALL
      ELSE
        WRITE(LUPRI,'(A,I5)') '* Unknown LCOMP = ',LCOMP
        CALL QUIT('MTAOCS: Unknown component LCOMP !!!')
      ENDIF
C
C     Allocate work space
C
#include "memint.h"
      NWRK = NVEC*NSPH
      CALL MEMGET('REAL',KW1,NWRK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KW2,NWRK,WORK,KFREE,LFREE)
C
C     First transform rows
C
      CALL VTAOCS(AOM,LDA,WORK(KW1),NSPH,NVEC,LCOMP,IPRINT)
      IF(IPRINT.GE.6) THEN
        CALL HEADER('MTAOCS: Half-transformed matrix:',-1)
        WRITE(LUPRI,'(2X,A,2I4)') 'NSPH,NVEC:',NSPH,NVEC
        CALL OUTPUT(WORK(KW1),1,NSPH,1,NVEC,NSPH,NVEC,-1,LUPRI)
      ENDIF
C
C     Transpose and transform columns
C
      CALL MTRSP(NSPH,NVEC,WORK(KW1),NSPH,WORK(KW2),NVEC)
      CALL VTAOCS(WORK(KW2),NVEC,WORK(KW1),NSPH,NSPH,LCOMP,IPRINT)
C
C     Backtranspose to obtain final matrix
C
      CALL MTRSP(NSPH,NSPH,WORK(KW1),NSPH,CSM,LDC)
C
C     Memory deallocation
      CALL MEMREL('MTAOCS',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('MTAOCS: Transformed (S) matrix:',-1)
        WRITE(LUPRI,'(2X,A,2I4)') 'NSPH,NSPH:',NSPH,NSPH
        CALL OUTPUT(CSM,1,NSPH,1,NSPH,LDC,NSPH,-1,LUPRI)
      ENDIF
      CALL QEXIT ('MTAOCS')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtaocs */
      SUBROUTINE MTAOSC(AOM,LDA,CSM,LDC,LCOMP,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Transform NVEC vectors from Spherical AO basis to Cartesian
C     AO basis
C
C     Options:
C     LCOMP = 0   Transform all components
C     LCOMP = 1   Transform only large components
C     LCOMP = 2   Transform only small components
C
C     NSPH - total number of spherical components
C
C     On input:  CSM in spherical components (NB,NB)
C     On output: AOM in cartesian components (NB,NB)
C
C
C  Written by Miro Ilias (after T.Saue's MTAOCS), Strasburg, 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include "mxcent.h"
C
#include "nuclei.h"
      DIMENSION AOM(LDA,*),CSM(LDC,*)
C     DIMENSION AOM(*),CSM(*)

      DIMENSION WORK(LWORK)
      CALL QENTER('MTAOSC')

C
C     Set limits
C
      NSPH = NSPHCM(LCOMP)

      IF    (LCOMP.EQ.0) THEN
C       All components
        NVEC = NBASIS
      ELSEIF(LCOMP.EQ.1) THEN
C       Large components
        NVEC = NLARGE
      ELSEIF(LCOMP.EQ.2) THEN
C       Small components
        NVEC = NSMALL
      ELSE
        WRITE(LUPRI,'(A,I5)') '* Unknown LCOMP = ',LCOMP
        CALL QUIT('MTAOCS: Unknown component LCOMP !!!')
      ENDIF

      IF(IPRINT.GE.5) THEN
        CALL HEADER('MTAOSC: Entering matrix in S AO:',-1)
        WRITE(LUPRI,'(2X,A,2I4)') 'size: NROWS, NCOLS',LDC,LDC
        WRITE(LUPRI,'(2X,A,2I4,A,2I4)')
     &   'desired:  (S) NSPH,NSPH:',NSPH,NSPH,
     &   ' -> (C) NVEC,NVEC:',NVEC,NVEC
        CALL OUTPUT(CSM,1,LDC,1,LDC,LDC,LDC,-1,LUPRI)
      ENDIF

CMI   NVEC >= NSPH !!

C
C     Allocate work space
C
#include "memint.h"
      NWRK = NVEC*NSPH
      NWRK2 = NVEC*NVEC

CMI   CALL MEMGET('REAL',KW1,NWRK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KW1,NWRK2,WORK,KFREE,LFREE)
CMI   CALL MEMGET('REAL',KW2,NWRK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KW2,NWRK2,WORK,KFREE,LFREE)
C
C     First transform rows
C
CMI   CALL VTAOCS(AOM,LDA,WORK(KW1),NSPH,NVEC,LCOMP,0)
CMI   CALL VTCSAO(CSM,LDC,WORK(KW1),NVEC,NSPH,LCOMP,IPRINT)
C ... in calling of this routine is the origin of the memmory fault !!!
C     CALL VTCSAO(CSM,LDC,WORK(KW1),NSPH,NVEC,LCOMP,IPRINT)
C     CALL VTCSAO(CSM,LDC,WORK(KW1),NSPH,NVEC,LCOMP,IPRINT)
C1    CALL VTCSAO(CSM,LDC,WORK(KW1),NVEC,NSPH,LCOMP,IPRINT)
      CALL VTCSAO(CSM,LDC,WORK(KW1),NVEC,NVEC,LCOMP,IPRINT)

CMI   ... we get WORK(KW1) ... NVECxNSPH
      IF(IPRINT.GE.6) THEN
        CALL HEADER('MTAOSC: Half-transformed matrix, WORK(KW1):',-1)
        WRITE(LUPRI,'(2X,A,2I4)') 'NVEC, NSPH:',NVEC,NSPH
C       CALL OUTPUT(WORK(KW1),1,NVEC,1,NSPH,NVEC,NSPH,-1,LUPRI)
        CALL OUTPUT(WORK(KW1),1,NVEC,1,NSPH,NVEC,NVEC,-1,LUPRI)
      ENDIF

C      WORK(KW1),NVECxNSPH -> WORK(KW2),NSPHxNVEC
C
C     Transpose and transform columns
C
CMI   CALL MTRSP(NSPH,NVEC,WORK(KW1),NSPH,WORK(KW2),NVEC)
C     CALL MTRSP(NVEC,NSPH,WORK(KW1),NVEC,WORK(KW2),NSPH)
C1    CALL MTRSP(NVEC,NSPH,WORK(KW1),NVEC,WORK(KW2),NVEC)
      CALL MTRSP(NVEC,NVEC,WORK(KW1),NVEC,WORK(KW2),NVEC)

      IF(IPRINT.GE.6) THEN
        CALL HEADER(
     &    'MTAOSC: Half-transformed transposed matrix, WORK(KW2):',-1)
        WRITE(LUPRI,'(2X,A,2I4)') 'NSPH, NVEC:',NSPH,NVEC
        CALL OUTPUT(WORK(KW2),1,NSPH,1,NVEC,NVEC,NVEC,-1,LUPRI)
      ENDIF


CMI ... we have WORK(KW2) NSPHxNVEC -> WORK(KW1) ... NVECxNVEC
CMI   CALL VTAOCS(WORK(KW2),NVEC,WORK(KW1),NSPH,NSPH,LCOMP,0)
C     CALL VTCSAO(WORK(KW2),NVEC,WORK(KW1),NVEC,NSPH,LCOMP,IPRINT)
      CALL VTCSAO(WORK(KW2),NVEC,WORK(KW1),NVEC,NVEC,LCOMP,IPRINT)

      IF(IPRINT.GE.6) THEN
        CALL HEADER(
     &    'MTAOSC: Transformed transposed matrix, WORK(KW1):',-1)
        WRITE(LUPRI,'(2X,A,2I4)') 'NVEC, NVEC:',NVEC,NVEC
        CALL OUTPUT(WORK(KW1),1,NVEC,1,NVEC,NVEC,NVEC,-1,LUPRI)
      ENDIF

C
C     Backtranspose to obtain final matrix WORK(KW1) ->  AOM
C
CMI   CALL MTRSP(NSPH,NSPH,WORK(KW1),NSPH,CSM,LDC)
CMI   CALL MTRSP(NVEC,NVEC,WORK(KW1),NVEC,AOM,LDA)
C1    CALL MTRSP(NVEC,NSPH,WORK(KW1),NVEC,AOM,LDA)
      CALL MTRSP(NVEC,NVEC,WORK(KW1),NVEC,AOM,LDA)
C
C     Memory deallocation
C  ... MI makes error, removed...
C   ... origin unknown...
      CALL MEMREL('MTAOSC',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('MTAOSC: Final transformed matrix in C AO:',-1)
        CALL OUTPUT(AOM,1,LDA,1,LDA,LDA,LDA,-1,LUPRI)
      ENDIF

      CALL QEXIT ('MTAOSC')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtcsao */
      SUBROUTINE VTCSAO(CSV,LDC,AOV,LDA,NVEC,LCOMP,IPRINT)
C***********************************************************************
C
C     Transform NVEC vectors from Spherical Components to Cartesian
C     AO basis.
C     NOTE: The transformation matrix is generally rectanuglar,
C           so the extra components are lost !!!!
C
C     Options:
C     LCOMP = 0   Transform all components
C     LCOMP = 1   Transform only large components
C     LCOMP = 2   Transform only small components
C
C     NSPH - number of spherical components
C
C     Written by T.Saue Sep 19 1995
C     Last revision Sep 19 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
C Used from COMMON blocks:
C    SHELLS: KMAX
C    SPHCOM:
C    NUCLEI: NLARGE,NSMALL,NBASIS
#include "nuclei.h"
#include "shells.h"
#include "dcbsph.h"
      DIMENSION AOV(LDA,NVEC),CSV(LDC,NVEC)
C
C     Set limits
C
      NSPH = NSPHCM(LCOMP)
      IF    (LCOMP.EQ.0) THEN
C       All components
        ISH1 = 1
        ISH2 = KMAX
      ELSEIF(LCOMP.EQ.1) THEN
C       Large components
        ISH1 = 1
        ISH2 = NLRGSH
      ELSEIF(LCOMP.EQ.2) THEN
C       Small components
        ISH1 = NLRGSH+1
        ISH2 = KMAX
      ELSE
        WRITE(LUPRI,'(A,I5)') '* Unknown LCOMP = ',LCOMP
        CALL QUIT('VTCSAO: Unknown component LCOMP !!!')
      ENDIF
C
C     Transform to spherical components
C
      IOFF = 1
      JOFF = 1
      DO ISHELL = ISH1,ISH2
        L     = NHKT(ISHELL)
        NLM   = 2*(L-1) + 1
        NXYZ  = L*(L+1)/2
        NDEG   = NUCDEG(NCENT(ISHELL))
        DO IDEG = 1,NDEG
          CALL DGEMM('T','N',NXYZ,NVEC,NLM,D1,CSP(ISPADR(L)),NLM,
     &               CSV(JOFF,1),LDC,D0,AOV(IOFF,1),LDA)
          IOFF = IOFF + NXYZ
          JOFF = JOFF + NLM
        ENDDO
      ENDDO
      IF(IPRINT.GE.5) THEN
        NDIM = IOFF-1
        CALL HEADER('VTCSAO: Transformed vectors:',-1)
        WRITE(LUPRI,'(2X,A,2I4)') 'NROWS, NCOLS',NDIM,NVEC
        CALL OUTPUT(AOV,1,NDIM,1,NVEC,LDA,NVEC,-1,LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nsphcm */
      FUNCTION NSPHCM(LCOMP)
C***********************************************************************
C
C     Count the number of spherical functions
C
C     Options:
C     LCOMP = 0   All components
C     LCOMP = 1   Only large components
C     LCOMP = 2   Only small components
C
C     Written by T.Saue Sep 19 1995
C     Last revision Sep 19 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include "maxorb.h"
#include "mxcent.h"
C
C Used from COMMON blocks:
C    SHELLS: KMAX,NHKT
C    NUCLEI: NLARGE,NSMALL,NBASIS
#include "nuclei.h"
#include "shells.h"
C
C     Set limits
C
      IF    (LCOMP.EQ.0) THEN
C       All components
        ISH1 = 1
        ISH2 = KMAX
      ELSEIF(LCOMP.EQ.1) THEN
C       Large components
        ISH1 = 1
        ISH2 = NLRGSH
      ELSEIF(LCOMP.EQ.2) THEN
C       Small components
        ISH1 = NLRGSH+1
        ISH2 = KMAX
      ELSE
        WRITE(LUPRI,'(A,I5)') '* Unknown LCOMP = ',LCOMP
        CALL QUIT('NSPHCM: Unknown component LCOMP !!!')
      ENDIF
C
C     Count the number of spherical functions
C
      NSPHCM = 0
      DO ISHELL = ISH1,ISH2
        L      = NHKT(ISHELL)
        NDEG   = NUCDEG(NCENT(ISHELL))
        NSPHCM = NSPHCM + (2*(L-1) + 1)*NDEG
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sphtst */
      SUBROUTINE SPHTST(WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Check normalization of tranformation matrix
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
C
#include "dcbsph.h"
#include "dcbbas.h"
      DIMENSION CSPSUM(((MXQN+1)*MXQN)/2),WORK(LWORK)
#include "memint.h"
C
      CALL HEADER('Output from SPHTST',-1)
CMI   CALL SPHGEN(WORK,LWORK,IPRINT)
      CALL SPHGEN(1,2,.FALSE.,WORK,LWORK,IPRINT)
C
C     Check normalization
C
      DO L = 1,MXQN
        LVAL = L-1
        NLM            = 2*LVAL + 1
        NXYZ           = (LVAL+1)*(LVAL+2)/2
        IOFF = ISPADR(L)
        DO K = 1,NXYZ
          CSPSUM(K) = SQRT(DDOT(NLM,CSP(IOFF),1,CSP(IOFF),1))
          IOFF = IOFF + NLM
        ENDDO
        WRITE(LUPRI,'(I5,5X,40F9.4)') (L-1),(CSPSUM(I),I=1,NXYZ)
      ENDDO
      NSPH = NSPHCM(0)
      WRITE(LUPRI,'(A,I5)') 'NSPH: ',NSPH
      CALL MEMGET('REAL',KA,NTBAS(0)*NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KB,NSPH*NTBAS(0),WORK,KFREE,LFREE)
      CALL SPHTS1(WORK(KA),WORK(KB),NSPH,WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('SPHTST',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sphts1 */
      SUBROUTINE SPHTS1(AMAT,BMAT,NSPH,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AMAT(NBAST,NBAST),BMAT(NSPH,NBAST),WORK(LWORK)
#else
      DIMENSION AMAT(NTBAS(0),NTBAS(0)),BMAT(NSPH,NTBAS(0)),WORK(LWORK)
#endif
C
C     GENERATE MATRIX
C
      OFF = 0.0D0
      DO J = 1,NTBAS(0)
        DO I = 1,NTBAS(0)
          OFF = OFF + 0.1D0
          AMAT(I,J) = OFF
        ENDDO
      ENDDO
      CALL HEADER('SPHTS1: Test matrix',-1)
      CALL OUTPUT(AMAT,1,NTBAS(0),1,NTBAS(0),NTBAS(0),NTBAS(0),
     &            -1,LUPRI)
      CALL VTAOCS(AMAT,NTBAS(0),BMAT,NSPH,NTBAS(0),0,IPRINT)
      CALL VTCSAO(BMAT,NSPH,AMAT,NTBAS(0),NTBAS(0),0,IPRINT)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mmsoao */
      SUBROUTINE MMSOAO(ASO,AAO,NBAST,IPRINT)
C***********************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine MMSOAO',-1)
C
C     Loop over all irreps in molecule
C
      ISTRA = 1
      CALL DZERO(AAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            AN     = D1/FMULT(MULA)
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA =
     &              PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
C                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     BN     = D1/FMULT(MULB)
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        AAO(INDA,INDB) = AAO(INDA,INDB)
     &                                 + FACA*FACB*ASO(I,J)
                     END IF
  600                CONTINUE
  500             CONTINUE
C                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         CALL HEADER('OUTPUT: Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bcktra */
      SUBROUTINE BCKTRA(COFT,LRT,LCT,COFU,LRU,LCU,NORB,NZC,
     &                  NVEC,JVEC,NBAS,TMAT,LR1,LC1,NZT,IPRINT)
C***********************************************************************
C
C     PURPOSE: Transform vectors (basis) from one-representation 
C              to the other
C
C     METHOD: Uses blas routine QGEMM
C
C     INPUT:
C	NORB - number of vectors to be tranformed
C	COFU - array of untransformed vectors of dimension 
C              (NORB,NORB,NZC)
C	COFT - array of backtransformed vectors of dimension 
C              (NBAS,NVEC,NZC)
C       TMAT - transformation matrix of dimension (NBAS,NORB,NZT)
C
C	Written by T.Saue January 1995
C	LAST VERSION: June 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0)
C
#include "dgroup.h"
      DIMENSION COFU(LRU,LCU,NZC),COFT(LRT,LCT,NZC),
     &          TMAT(LR1,LC1,NZT)

C   ...Define internal print level 
      IBPRNT = 12
C
      IF(NVEC.EQ.0) RETURN
      IF(IPRINT.GE.IBPRNT) THEN
        CALL HEADER('BCKTRA: Nontransformed coefficients',-1)
        WRITE(LUPRI,'(2(3X,A,I5))')
     +    'Start index: ',JVEC,'No. of vectors:',NVEC
        IF(ATOMIC .eqv. .FALSE.)THEN   !ayaki  Error appear in PRQMAT for ATOMIC
          CALL PRQMAT(COFU,NORB,NORB,LRU,LCU,NZC,IPQTOQ(1,0),LUPRI)
        ENDIF   
        CALL HEADER('BCKTRA: Transformation matrix:',-1)
        CALL PRQMAT(TMAT,NBAS,NORB,LR1,LC1,NZT,IPQTOQ(1,0),LUPRI)
      ENDIF
      CALL QGEMM(NBAS,NVEC,NORB,D1,
     &           'N','N',IPQTOQ(1,0),TMAT,LR1,LC1,NZT,
     &           'N','N',IPQTOQ(1,0),COFU(1,JVEC,1),LRU,LCU,NZC,
     &                D0,IPQTOQ(1,0),COFT,LRT,LCT,NZC)
      IF(IPRINT.GE.IBPRNT) THEN
        CALL HEADER('BCKTRA: Backtransformed coefficients',-1)
        WRITE(LUPRI,'(2(3X,A,I5))')
     +    'Start index: ',JVEC,'No. of vectors:',NVEC
        CALL PRQMAT(COFT,NBAS,NVEC,LRT,LCT,NZC,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck delsph */
      SUBROUTINE DELSPH(AMAT,LCOMP,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Delete totally symmetric contributions of spherical transformation
C     from LL or SS-block of matrix
C
C     LCOMP  = 1   LL-block
C     LCOMP  = 2   SS-block
C
C     Written June 20 1996 - Trond Saue
C     Last revision June 20 1996 - tsaue
C   
C     Used in BSS mode, M.ILIAS, Strasbourg, 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
C
#include "nuclei.h"
C
      DIMENSION AMAT(*),WORK(LWORK)
C
      CALL QENTER('DELSPH')
#include "memint.h"
      IF    (LCOMP.EQ.1) THEN
C     LL-block
        NDIM = NLARGE
        IOFF = 1
      ELSEIF(LCOMP.EQ.2) THEN
C     SS-block
        NDIM = NSMALL
        IOFF = NBASIS*NLARGE+NLARGE+1
      ELSE
        WRITE(LUPRI,'(A,I5)') 
     &    'DELSPH: Unknown option ',LCOMP
        CALL QUIT('DELSPH: Unknown option LCOMP !')
      ENDIF
C
      NSPH = NSPHCM(LCOMP)
      NWRK = NDIM*NSPH

      IF (IPRINT.GE.5) THEN
       WRITE(LUPRI,'(A,I6)') '* DELSPH: NSPH',NSPH
       WRITE(LUPRI,'(A)') '* DELSPH: Entering AMAT matrix '
       CALL OUTPUT(AMAT,1,NLARGE,1,NLARGE,NLARGE,NLARGE,-1,LUPRI)
      ENDIF
C     Allocate buffers
      CALL MEMGET('REAL',KW1,NWRK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KW2,NWRK,WORK,KFREE,LFREE)
C     First transform rows
      CALL VTAOCS(AMAT(IOFF),NBASIS,WORK(KW1),NSPH,
     &              NDIM,LCOMP,IPRINT)
C     ...then transpose
      CALL MTRSP(NSPH,NDIM,WORK(KW1),NSPH,WORK(KW2),NDIM)
C     ...and transform columns
      CALL VTAOCS(WORK(KW2),NDIM,WORK(KW1),NSPH,
     &              NSPH,LCOMP,IPRINT)
C     Backtransform columns
      CALL VTCSAO(WORK(KW1),NSPH,WORK(KW2),NDIM,
     &              NSPH,LCOMP,IPRINT)
C     ...then backtranspose
      CALL MTRSP(NDIM,NSPH,WORK(KW2),NDIM,WORK(KW1),NSPH)
C     ...and finally backtransform rows
      CALL VTCSAO(WORK(KW1),NSPH,AMAT(IOFF),NBASIS,
     &              NDIM,LCOMP,IPRINT)
C
C
C     Print section
C     =============
C
C      IF(IPRINT.GE.5) THEN
C        CALL HEADER('Output from DELSPH',-1)
C        CALL PRQMAT(AMAT(IOFF),NBASIS,NBASIS,NBASIS,NBASIS,1,
C     &              IPQTOQ(1,0,1),LUPRI)
C      ENDIF
      CALL QEXIT('DELSPH')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mksphc */
      SUBROUTINE MKSPHC(WORK,LWORK,IPRINT)
C***********************************************************************
C
C     This routine generates an orthogonal transformation from Cartesian
C     to solid harmonics
C
C     Written by T.Saue Aug 15 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxaqn.h"
C
      DIMENSION WORK(LWORK)
      CALL QENTER('MKSPHC')
#include "memint.h"
      NIJK = MXQN*MXQN*MXQN
      MXYZ = (MXQN*(MXQN+1))/2      
      NFUN = (MXQN*(MXQN+1)*(MXQN+2))/6
      CALL MEMGET('INTE',KINADR,3*NFUN,WORK,KFREE,LFREE)            
      CALL MEMGET('INTE',KINIJK,  NIJK,WORK,KFREE,LFREE)            
      CALL MEMGET('INTE',KIX,     MXYZ,WORK,KFREE,LFREE)            
      CALL MEMGET('INTE',KJY,     MXYZ,WORK,KFREE,LFREE)            
      CALL MEMGET('INTE',KKZ,     MXYZ,WORK,KFREE,LFREE)            
      CALL MEMGET('INTE',KLOFF,   MXQN,WORK,KFREE,LFREE)            
      CALL MEMGET('REAL',KDFAC,   MXQN,WORK,KFREE,LFREE)            
      CALL MKSPH1(WORK(KINADR),WORK(KINIJK),
     &            WORK(KIX),WORK(KJY),WORK(KKZ),WORK(KLOFF),
     &            WORK(KDFAC),
     &            WORK(KFREE),LFREE,IPRINT)
      CALL QEXIT('MKSPHC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mksph1 */
      SUBROUTINE MKSPH1(INADR,INIJK,IX,JY,KZ,LOFF,DFAC,
     &                  WORK,LWORK,IPRINT)
C***********************************************************************
C
C     This routine generates an orthogonal transformation from Cartesian
C     to solid harmonics
C
C     Written by T.Saue Aug 15 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
#include "maxaqn.h"
C
#include "sphtrm.h"
      DIMENSION INADR(3,*),INIJK(0:MXQN-1,0:MXQN-1,0:MXQN-1),
     &          IX(*),JY(*),KZ(*),LOFF(*),DFAC(*),WORK(*)
C
      CALL QENTER('MKSPH1')
C
      IND  = 0
      IOFF = 1
      FAC  = D1
      CALL DZERO(CSP,NCSP)
      DO LVAL = 0,MXQN-1
        NLM            = 2*LVAL + 1
        NXYZ           = (LVAL+1)*(LVAL+2)/2
        N2XYZ          = NXYZ*NXYZ
        ISPADR(LVAL+1) = IOFF
        IOFF           = IOFF + N2XYZ
        LOFF(LVAL+1)   = IND
        FAC            = FAC*NLM
        DFAC(LVAL+1)   = FAC
        CALL LMNVAL(LVAL+1,NXYZ,IX,JY,KZ)
        DO K = 1,NXYZ
          IND = IND + 1
          INADR(1,IND) = IX(K)
          INADR(2,IND) = JY(K)
          INADR(3,IND) = KZ(K)
        ENDDO
        IF (LVAL .EQ. 0) THEN
          CSP(1)       = D1
        ELSE IF (LVAL .EQ. 1) THEN
          CSP (2)      = D1
          CSP (6)      = D1
          CSP(10)      = D1
        ELSE
          CALL MKSPH2(LVAL,INADR,INIJK,LOFF,DFAC,WORK,LWORK,IPRINT)
        ENDIF
      ENDDO
C
      CALL QEXIT('MKSPH1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mksph2 */
      SUBROUTINE MKSPH2(LVAL,INADR,INIJK,LOFF,DFAC,
     &                  WORK,LWORK,IPRINT)
C***********************************************************************
C
C     This routine generates an orthogonal transformation from Cartesian
C     to solid harmonics for a given l-value.
C
C     Written by T.Saue Aug 15 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
#include "maxaqn.h"
C
      DIMENSION INADR(3,*),INIJK(0:MXQN-1,0:MXQN-1,0:MXQN-1),
     &          LOFF(*),DFAC(*),WORK(*)
C
      CALL QENTER('MKSPH2')
#include "memint.h"
      NXYZ           = (LVAL+1)*(LVAL+2)/2
      NNXYZ          = (NXYZ*(NXYZ+1))/2
      N2XYZ          = NXYZ*NXYZ
      CALL MEMGET('REAL',KSXYZ,NNXYZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOEF,N2XYZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF ,NXYZ ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,NXYZ ,WORK,KFREE,LFREE)
      CALL MKSPH3(LVAL,INADR,INIJK,LOFF,DFAC,WORK(KSXYZ),
     &            WORK(KCOEF),WORK(KBUF),WORK(KIBUF),
     &                  WORK(KFREE),LFREE,IPRINT)
      CALL QEXIT('MKSPH2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mksph3 */
      SUBROUTINE MKSPH3(LVAL,INADR,INIJK,LOFF,DFAC,
     &                  SXYZ,COEF,BUF,IBUF,
     &                  WORK,LWORK,IPRINT)
C***********************************************************************
C
C     This routine generates an orthogonal transformation from Cartesian
C     to solid harmonics for a given l-value.
C
C     Written by T.Saue Aug 15 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
#include "maxaqn.h"
C
#include "ccom.h"
      DIMENSION INADR(3,*),INIJK(0:MXQN-1,0:MXQN-1,0:MXQN-1),
     &          LOFF(*),DFAC(*),SXYZ(*),COEF(*),BUF(*),
     &          IBUF(*),WORK(*)
C
      NXYZ = (LVAL+1)*(LVAL+2)/2
C
C     Make overlap matrix in minimal basis
C
      CALL MKSXYZ(LVAL,SXYZ,INADR,LOFF,DFAC,IPRINT)      
C
C     Diagonalize minimal overlap matrix
C
      CALL DUNIT(COEF,NXYZ)
      CALL JACO(SXYZ,COEF,NXYZ,NXYZ,NXYZ,BUF,IBUF)
      IF(IPRINT.GE.4) THEN
        WRITE(LUPRI,'(A)') 
     &    '* Eigenvalues of minimal overlap matrix'
        II = 0
        DO I = 1,NXYZ
          II = (I*(I+1))/2
          WRITE(LUPRI,'(A,I5,F10.5)') '*EIG :',II,SXYZ(II)
        ENDDO
        WRITE(LUPRI,'(A)') 
     &    '* Eigenvectors of minimal overlap matrix'
        CALL OUTPUT(COEF,1,NXYZ,1,NXYZ,NXYZ,NXYZ,-1,LUPRI)
      
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mksxyz */
      SUBROUTINE MKSXYZ(LVAL,SXYZ,INADR,LOFF,DFAC,IPRINT)
C***********************************************************************
C
C     Make lower triangular row-packed overlap matrix in minimal 
C     Cartesian Gaussian basis for angular momentum value L
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1 = 1.0D0)
C
      DIMENSION INADR(3,*),LOFF(*),DFAC(*),SXYZ(*)
C
      NXYZ = (LVAL+1)*(LVAL+2)/2
      IOFF = LOFF(LVAL+1)
      IND  = 0
      DO I = 1,NXYZ
        DO J = 1,I
          IND = IND + 1
          SBUF = D1
          DO K = 1,3
            II = INADR(K,IOFF+I) + INADR(K,IOFF+J)
            IF(MOD(II,2).EQ.1) THEN
              SBUF = D0
            ELSEIF(II.GT.0) THEN
              II = II/2
              SBUF = SBUF*DFAC(II)
            ENDIF
          ENDDO
          SXYZ(IND) = SBUF
        ENDDO
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.4) THEN
        CALL HEADER('Output from MKSXYZ',-1)
        WRITE(LUPRI,'(A,I5)') '* L -value:',LVAL
        WRITE(LUPRI,'(A)') '* Minimal overlap matrix'
        CALL OUTPAK(SXYZ,NXYZ,NXYZ,LUPRI)
      ENDIF      
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck keepll */
      SUBROUTINE KEEPLL(AMAT,MZ,BASIS)
C***********************************************************************
C
C     Delete all parts of the matrix, except the LL block.
C
C     Written by L.Visscher, 6-1-2000
C     Last revision : 6-1-2000 - lvisscher
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
      DIMENSION AMAT(*)
      CHARACTER*3 BASIS
C
      IZOFF = 0
      DO IZ = 1,MZ
        DO JC = 1,2
          DO JREP = 0,MAXREP
            DO J = 0,NBBAS(JREP,JC)-1
              JSOFF = IZOFF + (IBBAS(JREP  ,JC)+J)*NTBAS(0)
              JUOFF = IZOFF + (ICOS (JREP+1,JC)+J)*NTBAS(0)
              DO IC = 1,2
                DO IREP = 0,MAXREP
                  IUOFF = JUOFF + ICOS (IREP+1,IC)
                  ISOFF = JSOFF + IBBAS(IREP  ,IC)
                  IF (IC+JC.GT.2.AND.NBBAS(IREP,IC).GT.0) THEN
                    IF (BASIS.EQ.'UNS') THEN
                       CALL DZERO(AMAT(IUOFF+1),NBBAS(IREP,IC))
                    ELSEIF (BASIS.EQ.'SRT') THEN
                       CALL DZERO(AMAT(ISOFF+1),NBBAS(IREP,IC))
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IZOFF = IZOFF + N2BBASX
      ENDDO
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck tsubbl */
      SUBROUTINE TSUBBL(IFRP,FAO,FMO,TM1,TM2,BUF,LBUF,IPRINT)
C***********************************************************************
C 
C     Perform AO->MO transformation only on diagonal blocks of a given irrep
C
C     Generalized routine TBOSON that was written by Trond Saue,
C     LV November 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C      
#include "dcbbas.h"
#include "dcborb.h"  
#include "dcbdhf.h"  
#include "dgroup.h" 
      DIMENSION FAO(*),FMO(*),TM1(*),TM2(*),BUF(LBUF)
C
      NDIM = NTMO(IFRP)*NTMO(IFRP)*NZ
      CALL DZERO(FMO(I2TMOT(IFRP)+1),NDIM)
      IBORB = 0    
      DO ISUB = 1,N_SUB_BL(IFRP)
      IF(NTMO_SUB(ISUB,IFRP,0).GT.0) THEN
        IAO = I2BASX(IFRP,IFRP)+1
        IMO = I2TMOT(IFRP) + (NTMO(IFRP)+1)*IBORB + 1
        ITM = I2TMT(IFRP)  + NFBAS(IFRP,0)*IBORB + 1

        CALL QTRANS('AOMO','S',0.0d0,NFBAS(IFRP,0),NFBAS(IFRP,0),
     &         NTMO_SUB(ISUB,IFRP,0),NTMO_SUB(ISUB,IFRP,0),
     &         FAO(IAO),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &         FMO(IMO),NTMO(IFRP),NTMO(IFRP),1,IPQTOQ(1,0),
     &         TM1(ITM),NFBAS(IFRP,0),NTMO(IFRP),NZT,IPQTOQ(1,0),
     &         TM2(ITM),NFBAS(IFRP,0),NTMO(IFRP),NZT,IPQTOQ(1,0),
     &         BUF,LBUF,IPRINT)
          IBORB = IBORB + NTMO_SUB(ISUB,IFRP,0)
      ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck tsbbl2 */
      SUBROUTINE TSBBL2(IFRP,IC,ISUB,FAO,FMO,TM1,TM2,BUF,LBUF,IPRINT)
C***********************************************************************
C 
C     Perform AO->MO transformation separatly Large and Small
C     (TSUBBL can transoform only total  )
C
C     Modified from TSUBBL
C     A Sunaga 2018
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C      
#include "dcbbas.h"
#include "dcborb.h"  
#include "dcbdhf.h"  
#include "dgroup.h" 
      DIMENSION FAO(*),FMO(*),TM1(*),TM2(*),BUF(LBUF)
C
      IF(ISUB.GE.1)THEN
        IBORB = 0 
        DO IS = 1, ISUB-1
           IBORB = IBORB + NTMO_SUB(IS,IFRP,0)
        ENDDO   
      ELSE
        IBORB = 0
      ENDIF
C      
      IF(NTMO_SUB(ISUB,IFRP,IC).GT.0) THEN
        IF(IC.EQ.1)THEN
          NTMO_BUF = NTMO_SUB(ISUB,IFRP,2)
        ELSEIF(IC.EQ.2)THEN
          NTMO_BUF = 0
        ENDIF
        IAO = I2BASX(IFRP,IFRP)+1
        IMO = 1
        ITM = I2TMT(IFRP)  +
     &        NFBAS(IFRP,0)*(IBORB+NTMO_BUF) + 1
C
        CALL QTRANS('AOMO','S',0.0d0,NFBAS(IFRP,0),NFBAS(IFRP,0),
     &         NTMO_SUB(ISUB,IFRP,IC),NTMO_SUB(ISUB,IFRP,IC),
     &         FAO(IAO),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &         FMO(IMO),NTMO_SUB(ISUB,IFRP,IC),NTMO_SUB(ISUB,IFRP,IC),
     &         1,IPQTOQ(1,0),
     &         TM1(ITM),NFBAS(IFRP,0),NTMO(IFRP),NZT,IPQTOQ(1,0),
     &         TM2(ITM),NFBAS(IFRP,0),NTMO(IFRP),NZT,IPQTOQ(1,0),
     &         BUF,LBUF,IPRINT)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MKSMOAO */
      SUBROUTINE MKSMOAO(SMOAO,TMAT,IPRINT,WORK,LWORK)
C**********************************************************************
C     PURPOSE:
C        Old: Construct half transformed overlap matrix, Smoao = Ct Sao
C        New: Construct half transformed overlap matrix, Smoao = Sao C
C
C        Written by J. Thyssen
C        Modified by T.Saue Aug 2000
C        Last revision: Jan 2 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION SMOAO(*),TMAT(*),WORK(*)
C
#include "memint.h"
C
      CALL QENTER('MKSMOAO')
C
      ISIZE = 0
      DO I = 1,NFSYM
         ISIZE = ISIZE + NZ*NFBAS(I,0)*NORB(I)
      END DO
      CALL DZERO(SMOAO,ISIZE)
C
C     Get AO overlap matrix
C
      CALL MEMGET('REAL',KSAO,N2BBASX,WORK,KFREE,LFREE)
      CALL GTOVLX(WORK(KSAO),SSMTRC)
C
      IOFF = 1
      DO I = 1,NFSYM
         IF ( NORB(I) .GT. 0 ) THEN
            CALL QGEMM(NORB(I),NFBAS(I,0),NFBAS(I,0),D1,
     &         'H','N',IPQTOQ(1,0),TMAT(I2TMT(I)+1),
     &         NFBAS(I,0),NORB(I),NZT,
     &         'N','N',IPQTOQ(1,0),WORK(KSAO+I2BASX(I,I)),
     &         NTBAS(0),NTBAS(0),1,
     &         D0,IPQTOQ(1,0),SMOAO(IOFF),
     &         NORB(I),NFBAS(I,0),NZ)
C
            IF ( IPRINT .GE. 10 ) THEN
               CALL HEADER('MKSMOAO: Left index transf. overlap'//
     &            ' matrix Smoao',-1)
               CALL PRQMAT(SMOAO(IOFF),NORB(I),NFBAS(I,0),
     &            NORB(I),NFBAS(I,0),NZ,IPQTOQ(1,0),LUPRI)
               CALL FLSHFO(LUPRI)
            END IF
         END IF
      IOFF = IOFF + NZ*NORB(I)*NFBAS(I,0)
      END DO
C
      CALL MEMREL('MKSMOAO',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MKSMOAO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MKSAOMO */
      SUBROUTINE MKSAOMO(SAOMO,TMAT,TWOC,IPRINT,WORK,LWORK)
C**********************************************************************
C
C     PURPOSE:
C   =============
C        Construct half transformed overlap matrix, Saomo = Sao T
C        From the relation T^(dagger) Sao T = 1 follows that
C        Sao T is the inverse of T^(dagger)
C
C    On input:  TMAT - transformation matrix, Y
C               TWOC - flag for setting the two-component scheme
C               note: after MAKE_BSS in GMOTRA we may be anyhow the 2c picture,
C                      so we use TWOC=.false.
C
C   On output:  SAOMO - half transformed overlap matrix
C
C   Note(MI,Aug05): It uses NZT,SSMTRC, I2TMT and NTMO variables !
C
C   Written by L. Visscher (based on MKSMOAO by J. Thyssen).
C   Last modifications: M.Ilias, aug 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION SAOMO(*),TMAT(*),WORK(*)
      LOGICAL TWOC
C
#include "memint.h"
C
      CALL QENTER('MKSAOMO')
C

      IF (IPRINT.GE.5) THEN
       CALL HEADER('*** Output from MKSAOMO ***',-1)
       write(lupri,'(2X,A,I1,A,I1)')
     & 'NZ=',NZ,' >>>> Dirac uses NZT=',NZT
       IF (TWOC) THEN
         write(lupri,'(2X,A)')
     &   'TWOC=.true. => Hard used two-component mode !'
       ELSE
         write(lupri,'(2X,A)')
     &  'TWOC=.false. => Used four-component mode'
       ENDIF
       IF (TWOCOMPBSS.AND.(BSS.or.x2c)) THEN
        WRITE(LUPRI,'(2X,A)') 
     &  'X2C/IOTC in two-component mode!' 
       ELSE IF (.NOT.TWOCOMPBSS.AND.(BSS.or.x2c)) THEN
        WRITE(LUPRI,'(2X,A)') 
     &  'X2C/IOTC in four-component mode!'
       ENDIF
      ENDIF
      
      ISIZE = 0
      DO I = 1,NFSYM
        IF (TWOC) THEN
          NREF = NESH(I) 
        ELSE
          NREF = NTMO(I)
        ENDIF
        ISIZE = ISIZE + NZT*NFBAS(I,0)*NREF
        IF (IPRINT.GE.5) THEN
          write(lupri,'(2X,A,I1,A,I1,A,I3)')
     &    'ifsym=',I,'/',NFSYM, ' NFBAS(I,0)=',NFBAS(I,0)
          write(lupri,'(2X,A,I3,A,I3,A,I3,A,I3,A,I3)')
     &    '>>>> NPSH(I)=',NPSH(I),
     &    ' NESH(I)=',NESH(I),' NTMO(I)=',NTMO(I),
     &    ' - NORB(I)=',NORB(I),'>>> used NREF=',NREF
        ENDIF
      ENDDO

      CALL DZERO(SAOMO,ISIZE)
C
C     Get AO overlap matrix
C
      CALL MEMGET('REAL',KSAO,N2BBASX,WORK,KFREE,LFREE)
!     SSMTRC is 1.0d0 when doing 4c-->2c for h1 == 4c-fock operator
      CALL GTOVLX(WORK(KSAO),SSMTRC)
C
      IOFF = 1
      ITMAT = 0
      DO I = 1,NFSYM
        IF (TWOC) THEN
          NREF = NESH(I) 
        ELSE
          NREF = NTMO(I)
          ITMAT = I2TMT(I)
        ENDIF
        IF ( NREF .GT. 0 ) THEN

            IF ( IPRINT .GE. 10 ) THEN
               CALL HEADER('MKSAOMO: Overlap matrix'//
     &            ' (from GTOVLX)',-1)
               WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &         '*** Fermion corep ',I,'/',NFSYM 
               write(lupri,*) ' SSMTRC=',SSMTRC
               CALL PRQMAT(WORK(KSAO+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
     &                     NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)

               CALL HEADER('MKSAOMO: Entering transformation matrix,'//
     &            ' TMAT (symm.adapted,NZT)',-1)
               WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &         '*** Fermion corep ',I,'/',NFSYM 
               CALL PRQMAT(TMAT(ITMAT+1),NFBAS(I,0),NREF,
     &            NFBAS(I,0),NREF,NZT,IPQTOQ(1,0),LUPRI)
               CALL FLSHFO(LUPRI)
            END IF

            CALL QGEMM(NFBAS(I,0),NREF,NFBAS(I,0),D1,
     &         'N','N',IPQTOQ(1,0),WORK(KSAO+I2BASX(I,I)),
     &         NTBAS(0),NTBAS(0),1,
     &         'N','N',IPQTOQ(1,0),TMAT(ITMAT+1),
     &         NFBAS(I,0),NREF,NZT,
     &         D0,IPQTOQ(1,0),SAOMO(IOFF),
     &         NFBAS(I,0),NREF,NZT)
C
            IF ( IPRINT .GE. 10 ) THEN
              CALL HEADER('MKSAOMO: Right index transf. overlap'//
     &          ' matrix SaoCmo, with NZT',-1)
              WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM 
              CALL PRQMAT(SAOMO(IOFF),NFBAS(I,0),NREF,
     &           NFBAS(I,0),NREF,NZT,IPQTOQ(1,0),LUPRI)
              CALL FLSHFO(LUPRI)
            END IF
         END IF
         IOFF = IOFF + NZT*NFBAS(I,0)*NREF
         IF (TWOC) ITMAT = ITMAT + (NFBAS(I,0)*NREF*NZT)
      END DO
C
      CALL MEMREL('SAOMO',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MKSAOMO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck LABLST */
      SUBROUTINE LABLST(ILAB,iprint)
C***********************************************************************
C
C     Use information in SO-labels to generate SO-basis
C
C     Written by T.Saue March 12 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "chrsgn.h"
C
#include "dcbbas.h"
#include "dcblab.h"
#include "symmet.h"
#include "nuclei.h"
#include "aosotr.h"
#include "pgroup.h"
      DIMENSION ILAB(2,*)
#include "dcbibt.h"
C
C     Scan SO-labels
C
      DO I = 1,NPLAB(0)
        ILAB(1,I)=0
      ENDDO
      DO I = 1,NBASIS
        JLAB = IPLAB(I,2)
        ILAB(1,JLAB) = ILAB(1,JLAB) + 1
        ILAB(2,JLAB) = I
      ENDDO
      if(iprint > 0)then
        CALL HEADER('Symmetry Orbitals',1)
        WRITE (LUPRI,'(1X,A,6X,8I6)')
     *     ' Number of orbitals in each symmetry: ',
     &      (NAOS(I),I=1,MAXREP+1)
          WRITE (LUPRI,'(1X,A,8I6)')
     *     ' Number of large orbitals in each symmetry: ',
     &      (NCOS(I,1),I=1,MAXREP+1)
          WRITE (LUPRI,'(1X,A,8I6)')
     *     ' Number of small orbitals in each symmetry: ',
     &      (NCOS(I,2),I=1,MAXREP+1)
      end if
C
C     Loop over large component labels
C
      JC   = -1
      DO I = 1,NPLAB(0)
        IC   = IGET(IATTR(I,2))
        IF(IC.NE.JC) THEN
          JC = IC
          JREP = -1
          IF(IC.EQ.1) THEN
            if(iprint > 0)then
              WRITE(LUPRI,'(/A)') '* Large component functions'
            end if
          ELSE
            if(iprint > 0)then
              WRITE(LUPRI,'(/A)') '* Small component functions'
            end if
          ENDIF
        ENDIF
        IREP = KGET(IATTR(I,2))
        IF(IREP.NE.JREP) THEN
          JREP = IREP
          if(iprint > 0)then
            WRITE (LUPRI,'(/2X,A,2X,A3,A1,I2,A1/)')
     &            'Symmetry',REP(IREP),'(',(IREP+1),')'
          end if
        ENDIF
        INUC = JGET(IATTR(I,2))
        NDEG = NUCDEG(INUC)
        II = ILAB(2,I)+ICOS(IREP+1,IC)-IBBAS(IREP,IC)
        IF(NDEG.GT.1) THEN
          if(iprint > 0)then
            WRITE(LUPRI,'(3X,I6,A,3X,A7,I1,7(A1,I1))')
     &        ILAB(1,I),' functions: ',
     &        PLABEL(I,2)(6:12),1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG)
          end if
        ELSE
          if(iprint > 0)then
            WRITE(LUPRI,'(3X,I6,A,3X,A7)') ILAB(1,I),' functions: ',
     &          PLABEL(I,2)(6:12)
          end if
        ENDIF
      ENDDO
C
      END
C***********************************************************************
C  /* Deck indaoc */
      SUBROUTINE INDAOC (INDAO,INDAOS,INDNAO)
C***********************************************************************
C
C      Sort AO-matrix from the default shell ordering to one where all 
C      functions of a given center are consecutive (but keep L, S overall
C      structure).
C
C      Written by L. Visscher august 2001
C      Last revision: Aug 15 2001 - LV 
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "symmet.h"
#include "aosotr.h"
#include "pgroup.h"
#include "nuclei.h"
#include "shells.h"
#include "dcbbas.h"
#include "ccom.h"
      DIMENSION INDAO(*),INDAOS(*)
      DIMENSION INDNAO(2,*)
C
C     Create index arrays from shell order to center order
C     ====================================================
C
C     NAORB  gives the position of the function in the default order
C     NAORBC gives the position of the function in the new order
C     The index arrays can be used to sort (INDAO and INDAOS) or
C     to give the start of the set of functions for a given center
C     (INDNAO). 
C
      NAORBC = 0
C
      JCENTD = 0
      DO JCENT = 1, NUCIND
        NDEGJ  = NUCDEG(JCENT)
        DO JDEG = 1, NDEGJ
          JCENTD = JCENTD + 1
          INDNAO(1,JCENTD) = NAORBC
          NAORB = 0
          DO ISHELL = 1,NLRGSH
            ICENT = NCENT(ISHELL)
            LVAL  = NHKT(ISHELL)
            NDEGI  = NUCDEG(ICENT)
            DO IDEG = 1,NDEGI
              DO ICOMP = 1,KHK(LVAL)
                NAORB = NAORB + 1
                IF (ICENT.EQ.JCENT.AND.IDEG.EQ.JDEG) THEN
                   NAORBC = NAORBC + 1
                   INDAOS(NAORB) = NAORBC
                   INDAO(NAORBC) = NAORB
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INDNAO(1,JCENTD+1) = NAORBC
      NAORBL = NAORB
C
      JCENTD = 0
      DO JCENT = 1, NUCIND
        NDEGJ = NUCDEG(JCENT)
        DO JDEG = 1, NDEGJ
          JCENTD = JCENTD + 1
          INDNAO(2,JCENTD) = NAORBC
          NAORB = NAORBL
          DO ISHELL = NLRGSH+1,KMAX
            ICENT = NCENT(ISHELL)
            LVAL  = NHKT(ISHELL)
            NDEGI  = NUCDEG(ICENT)
            DO IDEG = 1,NDEGI
              DO ICOMP = 1,KHK(LVAL)
                NAORB = NAORB + 1
                IF (ICENT.EQ.JCENT.AND.IDEG.EQ.JDEG) THEN
                   NAORBC = NAORBC + 1
                   INDAOS(NAORB) = NAORBC
                   INDAO(NAORBC) = NAORB
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INDNAO(2,JCENTD+1) = NAORBC
C
C     Sanity check : every function belongs to a center, which
C     means that the numbers NAORB and NAORBC should match at the end.
C
      IF (NAORB.NE.NAORBC) CALL QUIT ("Error in INDAOC")
C
      RETURN
C
      END
C***********************************************************************
      SUBROUTINE DTSOUAO(DMAT,SDMAT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Density transform: Dirac sorted SO-basis to Hermit sorted AO-basis
C     Can be made more efficient by combining routines BSTOBU and DSOTAO 
C
C     Written by L. Visscher november 2003 (based on part of LSQMAT)
C
C***********************************************************************
#include "implicit.h"
#include "maxorb.h"
#include "dcbbas.h"
      DIMENSION DMAT(*), SDMAT(*)
      DIMENSION WORK(*)
#include "dgroup.h"
C
      IPRINT = -10
      IF(NBSYM.GT.1) THEN
        CALL MEMGET('REAL',KBUF ,N2BBASX ,WORK,KFREE,LFREE)
        CALL DCOPY(N2BBASX,DMAT,1,WORK(KBUF),1)
        CALL BSTOBU(WORK(KBUF),NZ,WORK(KFREE),LFREE)
        CALL DTSOAO(WORK(KBUF),SDMAT,NTBAS(0),0,IPRINT)
        CALL MEMREL('DSOUAO',WORK,KBUF,KBUF,KFREE,LFREE)
      ELSE
        CALL DCOPY(N2BBASX,DMAT,1,SDMAT,1)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qbtrans */
      SUBROUTINE QBTRANS(KREP,TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,
     &                  FAO,LRAO,LCAO,NZAO,IQAO,
     &                  FMO,LRMO,LCMO,NZMO,IQMO,
     &                  TM1,LR1,LC1,NZTM1,IQTM1,IB1,
     &                  TM2,LR2,LC2,NZTM2,IQTM2,IB2,
     &                  BUF,LBUF,IPRINT)
C***********************************************************************
C
C     This routine performs unitary transformations indicated by TYP:
C
C       TYP = AOMO : AO-to-MO-transformation  FMO = (C+)FAOC
C       TYP = MOAO : MO-to-AO-transformation  FAO = CFMO(C+)
C
C     FADD gives the possibility of adding the results
C
C     TREV indicates symmetry of F-matrix under timereversal
C
C       TREV = 'S' - symmetric     (t =  1)
C       TREV = 'A' - anti-symmtric (t = -1)
C
C     The unitary C-matrix has a time-symmetric structure.
C
C     The AO-to-MO-transformation can be set up quaternionically as
C
C       FMO = [(Ca+)-t(CbT)j][(FAOa)+(FAOb)j][(Ca)+(Cb)j]
C
C       F(NBAS,NBAS,NZ) --->  F(NORB,NORB,NZ)
C
C     and the MO-to-Ao-transformation can be set up as:
C
C       FAO = [(Ca)+t(Cb)j][(FMOa)+(FMOb)j][(Ca+)-(Cb)j]
C
C       F(NORB,NORB,NZ) --->  F(NBAS,NBAS,NZ)
C
C     Written by T.Saue,January 1995
C     Last revision: 29 January 1997 : Luuk Visscher
C**********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
      CHARACTER TREV*1,TRA*1,TYP*4
      DIMENSION FAO(LRAO,LCAO,*),IQAO(*),
     &          FMO(LRMO,LCMO,*),IQMO(*),
     &          TM1(LR1,LC1,*),IQTM1(*),IB1(*),
     &          TM2(LR2,LC2,*),IQTM2(*),IB2(*),
     &          BUF(LBUF)
      DIMENSION IQW(4)
C
#include "ibtfun.h"
      CALL QENTER('QBTRANS')
C
C     Add the results on to the contents of the MO matrix ?
C
#ifdef DEBUG
      print*, 'QBTRANS: The AO matrix'
      CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IQAO,LUPRI)
      print*, 'QBTRANS: The transformation matrix TM1'
      CALL PRQMAT(TM1,NRAO,NCMO,LR1,LC1,NZTM1,IQTM1,LUPRI)
      print*, 'QBTRANS: The transformation matrix TM2'
      CALL PRQMAT(TM2,NRAO,NCMO,LR2,LC2,NZTM2,IQTM2,LUPRI)
#endif
      IF((NCMO.EQ.0).OR.(NRMO.EQ.0)) GOTO 999
C
C     Time reversal symmetry of F-matrix
C
      IF    (TREV.EQ.'S') THEN
        TRA = 'N'
      ELSEIF(TREV.EQ.'A') THEN
        TRA = 'I'
      ELSE
        CALL QUIT('QBTRANS: Unknown keyword TREV '//TREV)
      ENDIF
C
C     *******************************************************
C     *****  AO - to - MO -transformation : F' = (C+)FC *****
C     *******************************************************
C
      IF(TYP.EQ.'AOMO') THEN
C
C       Look for the most efficient transformation
C       ------------------------------------------
C
        IF (NCAO.GE.NRAO) THEN
C
C         Perform the two steps:
C           1. W  = FC
C           2. F' = (C+)W
C
          CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
          NBUF = NRAO*NCMO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QBTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF
C
C         First part of transformation: W = FC
C         ------------------------------------
C
          CALL QGEMM(NRAO,NCMO,NCAO,D1,
     &              'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &              'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &               D0,IQW,BUF,NRAO,NCMO,NZW)
C
C         Second part of transformation: H = (C+)W
C         ----------------------------------------
C
          IF(KREP.LT.0) THEN
C...     ...general case (spin-orbit)
            CALL QGEMM(NRMO,NCMO,NRAO,D1,
     &                 'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,
     &                 'N','N',IQW,BUF,NRAO,NCMO,NZW,
     &                 FADD,IQMO,FMO,LRMO,LCMO,NZMO)
          ELSE
C...     ...spin-free case: transform only symmetry blocks
C           of the real part corresponding to operator symmetry
            JREP = IBTXOR(KREP,IB2(1))
            NJ = 0
            IJ = 1
            DO J = 1,NCMO
              IF(J.EQ.NCMO) THEN
                IBJ = -1
              ELSE
                IBJ = IBTXOR(KREP,IB2(J+1))
              ENDIF
              NJ = NJ + 1
              IF(JREP.NE.IBJ) THEN
                IREP = IB1(1)
                NI = 0
                II = 1
                DO I = 1,NRMO
                  IF(I.EQ.NRMO) THEN
                   IBI = -1
                  ELSE
                   IBI = IB1(I+1)
                  ENDIF  
                  NI = NI + 1
                  IF(IREP.NE.IBI) THEN
                    IF(IREP.EQ.JREP) THEN
C...             ...transform this block (only real part)
                      IW = NRAO*(IJ-1)+1
                      CALL QGEMM(NI,NJ,NRAO,D1,
     &                   'H',TRA,IQTM1,TM1(1,II,1),LR1,LC1,NZTM1,
     &                   'N','N',IQW,BUF(IW),NRAO,NCMO,NZW,
     &                   FADD,IQMO,FMO(II,IJ,1),LRMO,LCMO,1)
                    ENDIF
                    IREP = IBI
                    NI = 0
                    II = I + 1
                  ENDIF
                ENDDO
                JREP = IBJ
                NJ = 0
                IJ = J + 1
              ENDIF
            ENDDO
          ENDIF
        ELSE
C
C         Perform the two steps:
C           1. W  = (C+)F
C           2. F' = WC
C
          CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
          NBUF = NRMO*NCAO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QBTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF

C
C         First part of transformation: W = (C+)F
C         ----------------------------------------
C
          CALL QGEMM(NRMO,NCAO,NRAO,D1,
     &       'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,
     &       'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &        D0,IQW,BUF,NRMO,NCAO,NZW)
C
C         Second part of transformation: H = WC
C         -------------------------------------
C
          IF(KREP.LT.0) THEN
C...     ...general case (spin-orbit)
            CALL QGEMM(NRMO,NCMO,NCAO,D1,
     &         'N','N',IQW ,BUF,NRMO,NCAO,NZW,
     &         'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &         FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
          ELSE
C...     ...spin-free case: transform only symmetry blocks
C           of the real part corresponding to operator symmetry
            JREP = IBTXOR(KREP,IB2(1))
            NJ = 0
            IJ = 1
            DO J = 1,NCMO
              IF(J.EQ.NCMO) THEN
                IBJ = -1
              ELSE
                IBJ = IBTXOR(KREP,IB2(J+1))
              ENDIF
              NJ = NJ + 1
              IF(JREP.NE.IBJ) THEN
                IREP = IB1(1)
                NI = 0
                II = 1
                DO I = 1,NRMO
                  IF(I.EQ.NRMO) THEN
                   IBI = -1
                  ELSE
                   IBI = IB1(I+1)
                  ENDIF  
                  NI = NI + 1
                  IF(IREP.NE.IBI) THEN
                    IF(IREP.EQ.JREP) THEN
C...             ...transform this block (only real part)
                      CALL QGEMM(NI,NJ,NCAO,D1,
     &                 'N','N',IQW ,BUF(II),NRMO,NCAO,NZW,
     &                 'N','N',IQTM2,TM2(1,IJ,1),LR2,LC2,NZTM2,
     &                 FADD,IQMO,FMO(II,IJ,1),LRMO,LCMO,1)
                    ENDIF
                    IREP = IBI
                    NI = 0
                    II = I + 1
                  ENDIF
                ENDDO
                JREP = IBJ
                NJ = 0
                IJ = J + 1
              ENDIF
            ENDDO
          ENDIF
C
        ENDIF
        IF(IPRINT.GE.5) THEN
          CALL HEADER('QBTRANS: AO-to-MO-transformed Fock matrix',-1)
          CALL PRQMAT(FMO,NRMO,NCMO,LRMO,LCMO,NZMO,IQMO,LUPRI)
        ENDIF
C
C     *******************************************************
C     *****  MO - to - AO -transformation : F' = CF(C+) *****
C     *******************************************************
C
      ELSEIF(TYP.EQ.'MOAO') THEN
C
C       Look for the most efficient transformation
C       ------------------------------------------
C
        IF(NCMO.GE.NRMO) THEN
          CALL IQPACK(IQMO,NZMO,IQTM2,NZTM2,IQW,NZW)
          NBUF = NRMO*NCAO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QBTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF
C
C         First part of transformation: W = F(C+)
C
          CALL QGEMM(NRMO,NCAO,NCMO,D1,
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       D0,IQW,BUF,NRMO,NCAO,NZW)
C
C         Second part of transformation: H = CW
C
          CALL QGEMM(NRAO,NCAO,NRMO,D1,
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQW ,BUF,NRMO,NCAO,NZW,
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
        ELSE
          CALL IQPACK(IQMO,NZMO,IQTM1,NZTM1,IQW,NZW)
          NBUF = NRAO*NCMO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QBTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF
C
C         First part of transformation: W = CF
C
          CALL QGEMM(NRAO,NCMO,NRMO,D1,
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,
     &       D0,IQW,BUF,NRAO,NCMO,NZW)
C
C         Second part of transformation: H = CW
C
          CALL QGEMM(NRAO,NCAO,NCMO,D1,
     &       'N','N',IQW ,BUF,NRAO,NCMO,NZW,
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
        ENDIF
       IF(IPRINT.GE.5) THEN
          CALL HEADER('QBTRANS: MO-to-AO-transformed Fock matrix',-1)
          CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IQAO,LUPRI)
        ENDIF
      ELSE
        CALL QUIT('QBTRANS: Unknown keyword TYP '//TYP)
      ENDIF
C
 999  CONTINUE
      CALL QEXIT('QBTRANS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bcktr1 */
      SUBROUTINE BCKTR1(COFT,LRT,LCT,COFU,LRU,LCU,NORB,NZC,
     &                  NVEC,JVEC,NBAS,TMAT,LR1,LC1,NZT,IPRINT)
C*****************************************************************************
C
C     PURPOSE: Transform vectors (basis) from one-representation to the other
C
C     METHOD: Uses blas routine QGEMM
C
C     INPUT:
C	NORB - number of vectors to be tranformed
C	COFU - array of untransformed vectors of dimension (NORB,NORB,NZC)
C       TMAT - transformation matrix of dimension (NBAS,NORB,NZT)
C     On output:
C	COFT - array of backtransformed vectors of dimension (NBAS,NVEC,NZC)
C
C	Written by T.Saue January 1995
C	LAST VERSION: June 27 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0)
C
#include "dgroup.h"
      DIMENSION COFU(LRU,LCU,NZC),COFT(LRT,LCT,NZC),
     &          TMAT(LR1,LC1,NZT)

C   ...Define internal print level 
      IBPRNT = 12 
C
      IF(NVEC.EQ.0) RETURN
      IF(IPRINT.GE.IBPRNT) THEN
        CALL HEADER('BCKTR1: Nontransformed coefficients',-1)
        WRITE(LUPRI,'(2(3X,A,I5))')
     +    'Start index: ',JVEC,'No. of vectors:',NVEC
        CALL PRQMAT(COFU,NORB,NORB,LRU,LCU,NZC,IPQTOQ(1,0),LUPRI)
        CALL HEADER('BCKTR1: Transformation matrix:',-1)
        CALL PRQMAT(TMAT,NBAS,NORB,LR1,LC1,NZT,IPQTOQ(1,0),LUPRI)
      ENDIF
      CALL QGEMM(NBAS,NVEC,NORB,D1,
     &           'N','N',IPQTOQ(1,0),TMAT,LR1,LC1,NZT,
     &           'N','N',IPQTOQ(1,0),COFU(1,JVEC,1),LRU,LCU,NZC,
     &                D0,IPQTOQ(1,0),COFT,LRT,LCT,NZC)
      IF(IPRINT.GE.IBPRNT) THEN
        CALL HEADER('BCKTR1: Backtransformed coefficients',-1)
        WRITE(LUPRI,'(2(3X,A,I5))')
     +    'Start index: ',JVEC,'No. of vectors:',NVEC
        CALL PRQMAT(COFT,NBAS,NVEC,LRT,LCT,NZC,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Bcktr2 */
      SUBROUTINE BCKTRA_FRZ(IFRP,IUNIT,COFT,COFU,EIG,IBEIG,
     &                  TMAT,IBUF,BUF,EBUF,IBBUF,IPRINT)
C***********************************************************************
C
C     Backtransform coefficients and insert the frozen orbitals
C     Written by T. Saue May 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION COFT(*),COFU(*),EIG(*),IBEIG(*),TMAT(*),
     &          BUF(*),EBUF(*),IBBUF(*),IBUF(*)
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dgroup.h"
C   ...Define internal print level 
      IBPRNT = 12 
      IF(IPRINT.GE.IBPRNT) THEN
        CALL HEADER('BCKTRA_FRZ: Nontransformed coefficients',-1)
        CALL PRQMAT(COFU,NTMO(IFRP),NTMO(IFRP),
     &            NTMO(IFRP),NTMO(IFRP),NZ,IPQTOQ(1,0),LUPRI)
        CALL HEADER('BCKTRA_FRZ: Transformation matrix:',-1)
        CALL PRQMAT(TMAT,NFBAS(IFRP,0),NTMO(IFRP),
     &              NFBAS(IFRP,0),NTMO(IFRP),NZT,IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Read pointer to positions of frozen orbitals
C
      CALL READI(IUNIT,NFRO(IFRP),IBUF)
      NBUF = NFBAS(IFRP,0)*NZ
      NOFF = NFBAS(IFRP,0)*NORB(IFRP)
      IFRO = 1
      JOFT = 1
      JOFU = 1
      CALL DCOPY(NTMO(IFRP),EIG,1,EBUF,1)        
      CALL ICOPY(NTMO(IFRP),IBEIG,1,IBBUF,1)        
      DO IFO=1,NFRO(IFRP)
        JFRO = IBUF(IFO)
        JEND = JFRO-1
        NVEC = JFRO-JOFT
C       Backtransform coefficients inbetween
        IF(NVEC.GT.0) THEN
          ITO = NFBAS(IFRP,0)*(JOFT-1)+1
          CALL BCKTR1(COFT(ITO),NFBAS(IFRP,0),NORB(IFRP),
     &                COFU,NTMO(IFRP),NTMO(IFRP),
     &                NTMO(IFRP),NZ,NVEC,JOFU,NFBAS(IFRP,0),
     &                TMAT,NFBAS(IFRP,0),NTMO(IFRP),NZT,0)
          CALL DCOPY(NVEC,EBUF(JOFU),1,EIG(JOFT),1)
          CALL ICOPY(NVEC,IBBUF(JOFU),1,IBEIG(JOFT),1)
          JOFU = JOFU+NVEC
          JOFT = JOFT+NVEC
        ENDIF
C       Insert frozen orbital
        CALL READT(IUNIT,NBUF,BUF)
        IB = 1
        IC = NFBAS(IFRP,0)*JEND+1
        DO IZ = 1,NZ
          CALL DCOPY(NFBAS(IFRP,0),BUF(IB),1,COFT(IC),1)        
          IB = IB + NFBAS(IFRP,0)
          IC = IC + NOFF
        ENDDO
        READ(IUNIT) EIG(JFRO),IBEIG(JFRO)
        JOFT = JOFT + 1
      ENDDO
C
C     Backtransform remaining coefficients
C
      NVEC = NORB(IFRP)+1-JOFT
      IF(NVEC.GT.0) THEN
        ITO = NFBAS(IFRP,0)*(JOFT-1)+1
        CALL BCKTR1(COFT(ITO),NFBAS(IFRP,0),NORB(IFRP),
     &                COFU,NTMO(IFRP),NTMO(IFRP),
     &                NTMO(IFRP),NZ,NVEC,JOFU,NFBAS(IFRP,0),
     &                TMAT,NFBAS(IFRP,0),NTMO(IFRP),NZT,0)
        CALL DCOPY(NVEC,EBUF(JOFU),1,EIG(JOFT),1)
        CALL ICOPY(NVEC,IBBUF(JOFU),1,IBEIG(JOFT),1)
      ENDIF
C
      IF(IPRINT.GE.IBPRNT) THEN
        CALL HEADER('BCKTRA_FRZ: Backtransformed coefficients',-1)
        CALL PRQMAT(COFT,NFBAS(IFRP,0),NORB(IFRP),
     &            NFBAS(IFRP,0),NORB(IFRP),NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENSMO_TO_DENSAO (DMO,DAO)
C***********************************************************************
C
C     Transform density matrix in MO form to AO format using coefficients
C     taken from DFCOEF
C     Written by L. Visscher, March 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
C
      REAL*8  DMO(NORBT,NORBT,NZ), DAO(N2BBASX,NZ)
      REAL*8, ALLOCATABLE :: CMO(:),WORK(:)
      LOGICAL FILEX
C
C     Allocate the memory that is required by QTRANS (REACMO will not
C     use the work memory). QTRANS needs to store the halftransformed
C     matrix that has the same dimension as the MO coefficients
C
      KFREE = 1
      LFREE = NCMOTQ
      ALLOCATE (WORK(LFREE))
      ALLOCATE (CMO(NCMOTQ))
     
C
C     Read the MO coefficients
      INQUIRE(FILE='DFCOEF',EXIST=FILEX)
      IF(.NOT.FILEX) CALL QUIT ('DFCOEF missing')
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPEXP')
      CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,IDUM,DUM,2)
      CLOSE(LUCOEF,STATUS='KEEP')

C     Perform the transformation
      DO IFSYM = 1,NFSYM
          CALL QTRANS('MOAO','S',D0,
     &       NFBAS(IFSYM,0),NFBAS(IFSYM,0),
     &       NORB(IFSYM),NORB(IFSYM),
     &       DAO(I2BASX(IFSYM,IFSYM)+1,1),NTBAS(0),NTBAS(0),
     &       NZ,IPQTOQ(1,0),
     &       DMO(IORB(IFSYM)+1,IORB(IFSYM)+1,1),NORBT,NORBT,
     &       NZ,IPQTOQ(1,0),
     &       CMO(ICMOQ(IFSYM)+1),NFBAS(IFSYM,0),NORB(IFSYM),
     &       NZ,IPQTOQ(1,0),
     &       CMO(ICMOQ(IFSYM)+1),NFBAS(IFSYM,0),NORB(IFSYM),
     &       NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,1)
      END DO
   
C     Free the memory
      DEALLOCATE (WORK)
      DEALLOCATE (CMO)

      END
C***********************************************************************
C/* Deck Geigsolv */
      SUBROUTINE GEIGSOLV(AMAT,SMAT,EIG,N,NEFF,LRM,LCM,NZ,IPQ,
     &                    IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Solve general Hermitian eigenvalue problem: Ac = eSc
C     INPUT:
C        AMAT - Hermitian A
C        SMAT - Hermitian metric S
C     OUTPUT:  
C        AMAT - scramble
C        EIG  - eigenvalues
C        SMAT - eigenvectors
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(STOL=1.0D-6,D0=0.0D0,D1=1.0D0)
      DIMENSION AMAT(LRM,LCM,NZ),SMAT(LRM,LCM,NZ),EIG(N),IPQ(*),WORK(*)
C
      KFRSAV=KFREE
C.....Get canonical orthogonalization matrix
      CALL MEMGET('REAL',KVMAT,N*N*NZ,WORK,KFREE,LFREE)
      CALL QLOWDIN(WORK(KVMAT),N,N,SMAT,LRM,LCM,EIG,N,NEFF,NZ,STOL,
     &             WORK,KFREE,LFREE)
C.....Transform AMAT to orthonormal basis; store result in SMAT
      CALL QTRANS('AOMO','S',D0,N,N,NEFF,NEFF,
     &            AMAT,LRM,LCM,NZ,IPQ,SMAT,LRM,LCM,NZ,IPQ,
     &            WORK(KVMAT),N,N,NZ,IPQ,WORK(KVMAT),N,N,NZ,IPQ,
     &            WORK(KFREE),LFREE,IPRINT)
C.....Diagonalize AMAT(SMAT) in orthonormal basis; vectors stored in AMAT
      CALL QDIAG(NZ,NEFF,SMAT,LRM,LCM,EIG,1,AMAT,LRM,LCM,
     &           WORK(KFREE),LFREE,IERR)
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,*) 'GEIGSOLV: Eigenvalues...'
        DO I = 1,NEFF
          WRITE(6,*) I,EIG(I)
        ENDDO
      ENDIF
C.....Backtransform eigenvectors -> store in SMAT
      CALL QGEMM(N,NEFF,NEFF,D1,
     &           'N','N',IPQ,WORK(KVMAT),N,N,NZ,
     &           'N','N',IPQ,AMAT,LRM,LCM,NZ,
     &                D0,IPQ,SMAT,LRM,LCM,NZ)
      IF(IPRINT.GE.5) THEN
        WRITE(LUPRI,*) 'GEIGSOLV: Eigenvectors....'
        CALL PRQMAT(SMAT,N,NEFF,LRM,LCM,NZ,IPQ,LUPRI)
      ENDIF
      CALL MEMREL('GEIGSOLV',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Qlowdin */
      SUBROUTINE QLOWDIN(VMAT,LRV,LCV,SMAT,LRS,LCS,EIG,N,NV,NZ,STOL,
     &                   WORK,KFREE,LFREE)
C***********************************************************************
C
C     Routine for the generation of Lowdin canonical orthogonalization
C     Written by Trond Saue June 10 2012
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0,DM1=-1.0D0)
      DIMENSION VMAT(LRV,LCV,NZ),SMAT(LRS,LCS,NZ),EIG(N),WORK(*)
      NSDIM=LRS*LCS*NZ
C.....Scale overlap matrix to get eigenvalues in the correct order
C     QDIAG gives ascending and we want descending....
      CALL DSCAL(NSDIM,DM1,SMAT,1)
      CALL QDIAG(NZ,N,SMAT,LRS,LCS,EIG,1,VMAT,LRV,LCV,
     &           WORK(KFREE),LFREE,IERR)
      NV = 0
      DO I = 1,N
        TMP = DABS(EIG(I))
        IF(TMP.LE.STOL) GOTO 10
        NV  = NV + 1
        FAC = D1/SQRT(TMP)
        DO IZ = 1,NZ
          CALL DSCAL(N,FAC,VMAT(1,I,IZ),1)
        ENDDO
      ENDDO
 10   CONTINUE
C
      RETURN
      END

