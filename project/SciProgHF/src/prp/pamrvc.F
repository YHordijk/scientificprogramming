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

! define task symbols for CALL DIRAC_PARCTL( task )
#include "dirac_partask.h"

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrssvc */
      SUBROUTINE XRSSVC(BCVEC,BEVEC,BPVEC,
     &                  SCVEC,SEVEC,SPVEC,IBTYP,
     &                  IBCVC,IBEVC,IBPVC,NPOS,
     &                  GPCI,GPOE,GPOP,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate sigma vector SVEC:
C     molecular Hessian E[2] times trial vector B(I)
C
C        SVEC(I) = E[2]*B(I)
C
C     The trial vectors come in three types:
C        C-type: contains only configurational parameters
C        E-type: contains only electron-eelectron (e-e) rotations
C        P-type: contains only electron-positron  (e-p) rotations
C
C     The sigma vectors are split up into a C, E and P -part.
C
C     Written by H.J.Aa.Jensen and Trond Saue 1996
C     Last revision July 31 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
C
      DIMENSION BCVEC(*),BEVEC(*),BPVEC(*)
      DIMENSION SCVEC(*),SEVEC(*),SPVEC(*)
      DIMENSION IBCVC(*),IBEVC(*),IBPVC(*)
      DIMENSION IBTYP(2,*),NPOS(*),WORK(*)
      DIMENSION GPCI(*), GPOE(*), GPOP(*)
C
C Used from common blocks:
C   DCBXRS : LZYVAR,NCSIM,NESIM,NPSIM,NOSIM
C
#include "dcbxrs.h"
C
C
      CALL QENTER('XRSSVC')
      CALL RTKTIME(.TRUE.,18)
C
C     Initialize sigma vectors
C     ========================
C
      IF(NZCONF.GT.0) CALL DZERO(SCVEC,NTSIM*NZCONFQ)
      IF(NZXOPE.GT.0) CALL DZERO(SEVEC,NTSIM*NZXOPEQ)
      IF(NZXOPP.GT.0) CALL DZERO(SPVEC,NTSIM*NZXOPPQ)
C
C     Construct orbital part of sigma vectors
C     ====================================================
C
      IF (NZXOPT .GT. 0) THEN
         CALL XRSSEP(BCVEC,BEVEC,BPVEC,
     &        SCVEC,SEVEC,SPVEC,
     &        IBCVC,IBEVC,IBPVC,
     &        IBTYP,GPOE,GPOP,NPOS,WORK,KFREE,LFREE)
      END IF
C
C     Construct configurational part of sigma vectors
C     for the configurational trial vectors
C     ===============================================
C
      IF (NCSIM.GT.0) THEN
         CALL XRSSCI(BCVEC,SCVEC,IBCVC,WORK,KFREE,LFREE)
      END IF
C
      CALL RTKTIME(.FALSE.,18)
      CALL QEXIT('XRSSVC')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrssci */
      SUBROUTINE XRSSCI(BCVEC,SCVEC,IBCVC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate conf. contribution to conf. sigma vector.
C     (Equation 3.29 in the KRMC paper (JCP 104, 4083 (1996)).
C
C     Input:
C       BCVEC: the CI trial vectors.
C       IBCVC: index of the CI trial vectors.
C     from common blocks:
C       FKRMC: Wave function opt. (.TRUE.) or linear response calc?
C       NCSIM: number of CI trial vectors.
C
C     Output:
C       SCVEC: the CI sigma vectors.
C
C     Written by J. Thyssen - Nov 23 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION BCVEC(*),SCVEC(*),WORK(*),IBCVC(*)
C
C
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcborb.h"
C
      CALL QENTER('XRSSCI')
      CALL RTKTIME(.TRUE.,11)
      KFRSAV = KFREE
C
      IF (JSYMOP .NE. 1 .OR. JTIMOP .NE. 1) THEN
C        hjaaj March 2002: we need to implement in
C        the CI routines for BCVEC not of same symmetry
C        as CREF, and for time-variant BCVECs. Maybe
C        IBCVC is useful for that. TODO
         CALL QUIT('Only JSYMOP.eq.1 and JTIMOP.eq.1 implemented')
      END IF
C
C     Get FCAC and H2AC
C
      LH2AC = NASHT * NASHT * NNASHX * NZ * 3
      CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2AC,LH2AC,WORK,KFREE,LFREE)
C
      CALL GETMAT('FCAC',WORK(KFCAC),IPRXRS,FKRMC,WORK(KFREE),LFREE)
      CALL GETMAT('H2AC',WORK(KH2AC),IPRXRS,FKRMC,WORK(KFREE),LFREE)
C
      IF (XRS_CIPROGRAM .EQ. 'GASCIP' .OR.
     &    XRS_CIPROGRAM .EQ. 'LUCIAREL') THEN
C        use routine with MOLFDIR format for integrals
         CALL XRSSCIM(BCVEC,SCVEC,
     &      WORK(KFCAC),WORK(KH2AC),WORK,KFREE,LFREE)
      ELSE
         WRITE(LUPRI,'(/A/2A)')
     &        '*** ERROR in XRSSCI ***',
     &        'No sigma vector implemented for CI program ',
     &        XRS_CIPROGRAM
         CALL QUIT('*** ERROR in XRSSCI ***')
      END IF
      CALL MEMREL('XRSSCI',WORK,1,KFRSAV,KFREE,LFREE)
      CALL RTKTIME(.FALSE.,11)
      CALL QEXIT ('XRSSCI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrsscim */
      SUBROUTINE XRSSCIM(BCVEC,SCVEC,FCAC,H2AC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     XRSSCIm: Calculate conf. contribution to conf. sigma vector.
C     (Equation 3.29 in the KRMC paper (JCP 104, 4083 (1996)).
C
C     xrssciM: This routine converts integrals to MOLFDIR format before
C     calling the CI program.
C
C     Input:
C       BCVEC: the CI trial vectors.
C       FCAC : active-active part of FC fock matrix
C       H2AC : two-electron integrals.
C       NCSIM: number of CI trial vectors.
C
C     Output:
C       SCVEC(1,NOSIM+1:NOSIM+NCSIM): the CI sigma vectors.
C
C     Written by J. Thyssen - Nov 23 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
      DIMENSION BCVEC(NZCONFQ,*),SCVEC(NZCONFQ,*)
      DIMENSION FCAC(*),H2AC(*),WORK(*)
C
#include "dcbxrs.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CALL QENTER('XRSSCIM')
      KFRSAV = KFREE
      MZ = MIN(NZ,2)
C
C     Transform quarternion FCAC (Dirac MO-basis) to Molfdir MO-basis
C     ---------------------------------------------------------------
C
      CALL MEMGET('REAL',KMOLFFC,(2*NASHT)*(2*NASHT)*MZ,
     &             WORK,KFREE,LFREE)
      CALL QFC2MFC(FCAC,WORK(KMOLFFC),1,1,IPRXRS)
C
      IF (XRS_CIPROGRAM .EQ. 'GASCIP') THEN
C
C        Transform integrals from Dirac (NZ,3) basis to Molfdir basis
C        ------------------------------------------------------------
C
         CALL MEMGET('REAL',KMUUUU,(2*NASHT)**4 * 2, WORK,KFREE,LFREE)
         CALL DZERO(WORK(KMUUUU),(2*NASHT)**4 * 2)
C
         CALL DNZ32M(H2AC,WORK(KMUUUU),IPRXRS)
C
      END IF
C
C     Loop over trial vectors
C     -----------------------
C
      DO I  = 1, NCSIM
C
C        Call CI sigma vector module for calculating H * B
C        -------------------------------------------------
C
         IF (XRS_CIPROGRAM .EQ. 'GASCIP') THEN
            CALL GASCIP_SVC(.FALSE.,DUMMY,SCVEC(1,I+NOSIM),
     &                     NZCONF,WORK(KZCONF),BCVEC(1,I),
     &                     WORK(KMOLFFC),WORK(KMUUUU),IPRXRS)
         ELSE IF (XRS_CIPROGRAM .EQ. 'LUCIAREL') THEN
            CALL LUCI_SIGMA(.FALSE.,DUMMY,BCVEC(1,I),SCVEC(1,I+NOSIM),
     &                      WORK(KMOLFFC),H2AC,WORK,KFREE,LFREE)
         ELSE
            CALL QUIT('Unknown CI program requested!')
         END IF
C
C        Print ?
C
         IF ( IPRXRS .GE. 30 ) THEN
            WRITE(LUPRI,'(A,I2,A)') '( XRSSCIM)   B(',I,')>'
            IF (MZ.EQ.2) WRITE(LUPRI,'(T12,A)') '(real and imag parts)'
            CALL OUTPUT(BCVEC(1,I),1,NZCONF,1,MZ,
     &                  NZCONF,MZ,1,LUPRI)
            WRITE(LUPRI,'(A,I2,A)') '( XRSSCIM) H|B(',I,')>'
            IF (MZ.EQ.2) WRITE(LUPRI,'(T12,A)') '(real and imag parts)'
            CALL OUTPUT(SCVEC(1,I+NOSIM),1,NZCONF,1,MZ,
     &                  NZCONF,MZ,1,LUPRI)
         END IF
C
C        Finish sigma vector: - E[0] b_{\mu}
C
         CALL DAXPY(NZCONFQ,-ENRGY,
     &              BCVEC(1,I),1,SCVEC(1,I+NOSIM),1)
C
      END DO
C
      CALL MEMREL('XRSSCIM',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT ('XRSSCIM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrsphp */
      SUBROUTINE XRSPHP(SCVEC,WORK,LWORK)
C***********************************************************************
C
C     Project out CREF from conf. sigma vectors.
C
C     Input:
C       SCVEC: the CI sigma vectors.
C
C     Output:
C       SCVEC: the CI sigma vectors with CREF projected out.
C
C     Written by J. Thyssen - Nov 23 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER ( THREPS = 1.0D-12 )
#include "consts.h"
C
#include "dgroup.h"
#include "dcbxrs.h"
C
      DIMENSION SCVEC(NZCONFQ,*),WORK(*)
C
#include "memint.h"
C
      CALL QENTER('XRSPHP')
C
C     Get CI ref. vector from record 1 on LUBCI.
C
      CALL MEMGET('REAL',KCREF,NZCONFQ,WORK,KFREE,LFREE)
!     CALL READAC(LUBCI,NZCONFQ,WORK(KCREF),1)
      CALL READ_DX(LUBCI,1,NZCONFQ,WORK(KCREF))
C
C     Do projection.
C
      KCREFI = KCREF + NZCONF
C     ... pointer to imaginary part of CREF
      DO ICSIM = 1, NTSIM
         TTR = DDOT(NZCONF,WORK(KCREF),1,SCVEC(1,ICSIM),1)
         IF (NZ .GE. 2) THEN
            TTR = TTR
     &          - DDOT(NZCONF,WORK(KCREFI),1,SCVEC(1+NZCONF,ICSIM),1)
            TTI =
     &           DDOT(NZCONF,WORK(KCREF ),1,SCVEC(1+NZCONF,ICSIM),1)
     &         + DDOT(NZCONF,WORK(KCREFI),1,SCVEC(1,ICSIM),1)
         ELSE
            TTI = D0
         END IF
C
C        CI ref. vector is normalized to 0.5
C
         TTR = D4 * TTR
         TTI = D4 * TTI
C
         IF ( SQRT(TTR*TTR+TTI*TTI) .GT. THREPS) THEN
            CALL DAXPY(NZCONF,-TTR,WORK(KCREF),1,SCVEC(1,ICSIM),1)
            IF (NZ .GE. 2) THEN
               CALL DAXPY(NZCONF, TTI,WORK(KCREFI),1,SCVEC(1,ICSIM),1)
               CALL DAXPY(NZCONF,-TTR,WORK(KCREFI),1,
     &              SCVEC(1+NZCONF,ICSIM),1)
               CALL DAXPY(NZCONF,-TTI,WORK(KCREF ),1,SCVEC(1,ICSIM),1)
            END IF
            IF (IPRXRS .GE. 5) WRITE(LUPRI,9000) ICSIM,TTR,TTI
         END IF
      END DO
C
 9000 FORMAT(/' CI ref. vector projected out of sigma vector no.',I3,
     &       /'    - overlap was ',1P,D20.12,' + ',D20.12,'i')
C
      IF(IPRXRS.GE.30) THEN
         WRITE(LUPRI,'(/A,I5,A)') ' Real configuration part of',
     &        NCSIM,' configurational sigma vector(s)'
         CALL OUTPUT(SCVEC,1,NZCONF,1,NCSIM,NZCONFQ,1,1,LUPRI)
         IF (NZ .GE. 2) THEN
            WRITE(LUPRI,'(/A,I5,A)') ' Imaginary configuration part of',
     &        NCSIM,' configurational sigma vector(s)'
            CALL OUTPUT(SCVEC(1+NZCONF,1),1,NZCONF,1,NCSIM,
     &           NZCONFQ,1,1,LUPRI)
         END IF
      END IF
C
      CALL MEMREL('XRSPHP',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT ('XRSPHP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrssep */
      SUBROUTINE XRSSEP(BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,
     &     IBCVC,IBEVC,IBPVC,
     &     IBTYP,GPOE,GPOP,NPOS,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Construct orbital part of sigma vectors
C
C     SA_(pq) = E[2](o*o)_(pq,rs)Z_(rs) + E[2](o*o*)_(pq,rs)Y^*_(rs)
C             = -<0|[Qp,[Rs,H]]|0>Z_(rs) + <0|[Qp,[Sr,H]]|0>Y^*_(rs)
C             = -<0|[Qp,[Rs,H]]|0>Z_(rs) - <0|[Qp,[Sr,H]]|0>{-Y+}_(sr)
C             = -<0|[Qp,[Rs,H]]|0>WA_(r,s)
C
C     where the matrix WA has the structure:
C
C           P  I   A   S             P  I  A  S
C          --------------         ---------------
C        P | 0  0   0   0       P |  0  Z  Z  0
C        I | 0  0  -Y+ -Y+      I | -Y+ 0  0  0
C        A | 0  Z   X  -Y+      A | -Y+ 0  0  0
C        S | 0  Z   Z   0       S |  0  0  0  0
C
C           e-e rotations           e-p rotations
C
C     with the convention:
C        P - positronic spinors
C        I - inactive spinors
C        A - active spinors
C        S - secondary spinors
C
C     (X has -Y+ in upper diagonal and Z in lower diagonal)
C
C     Similarly
C
C     SB_(pq) = E[2](o*o)_(pq,rs)Y_(rs) + E[2](o*o*)_(pq,rs)Z^*_(rs)
C             = -<0|[Qp,[Rs,H]]|0>WB_(r,s)
C
C     where WB is identical to WA with Y and Z interchanged.
C
C-----------------------------------------------------------------------
C
C     The key to efficient evaluation of SVEC is to reformulate it
C     as a gradient element:
C
C     SA_(pq) = -<0|[Qp,HA]]|0>
C
C     where HA is a one-index transformed Hamiltonian
C
C     HA = [Tu,H] WA_(t,u) = ha_(pq)Pq + ga_(pq,rs)PRsq
C
C     with
C
C     ha_(pq)    = WA_(pt)h_(tq)    - h_(pt)WA_(tq)
C
C     ga_(pq,rs) = WA_(pt)g_(tq,rs) - g_(pt,rs)WA_(tq)
C                + WA_(rt)g_(pq,ts) - g_(pq,rt)WA_(ts)
C
C     In RPA we need to evaluate
C
C     SA_(si) = -FA_(si) - GA_(si)
C     SB_(si) = -FB_(si) - GB_(si)
C
C     with
C
C     FA_(si) = WA_(st)F_(ti) - F_(st)WA_(ti)
C
C     GA_(si) = [g_(si,xy)-g_(sy,xi)] DA_(yx)
C
C     where x and y are AO-indices and
C
C     DA_(yx) = DA1_(yx) - DA2_(yx)
C
C     DA1_(yx) = C_(yj)(-Y+)_(js)(C+)_(sx)
C     DA2_(yx) = C_(ys)( Z )_(sj)(C+)_(jx)
C
C-----------------------------------------------------------------------
C
C     For efficient evaluation we shall exploit the relations
C
C     Time - symmetric operators:      DA1 = -[DA2^+]
C     Time - anti symmetric operators: DA1 = -[-i(DA2^+)i] (i-transform)
C
C
C-----------------------------------------------------------------------
C
C     The algorithm proceeds in the following steps:
C
C     1. Generate standard and modified two-electron Fock matrices.
C     2. Accumulate GA/GB into SVEC
C     3. Accumulate FA/FB into SVEC
C
C     Step 1 is the most memory and CPU-intensive step, but after
C     step 2 most of the memory can be released.
C
C     The b-vectors are scatterd into W matrices twice to
C     minimize memory.
C
C
C     Written by T.Saue and H.J.Aa.Jensen 1996
C     Last revision Aug 2 1996 - tsaue
C
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0)
      DIMENSION BCVEC(*),BEVEC(*),BPVEC(*),
     &          SCVEC(*),SEVEC(*),SPVEC(*),
     &          IBCVC(*),IBEVC(*),IBPVC(*),IBTYP(2,*),NPOS(*),
     &          WORK(*)
      DIMENSION GPOE(*),GPOP(*)
C
#include "dcborb.h"
#include "dcbham.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxrs.h"
C
      CALL QENTER('XRSSEP')
      CALL RTKTIME(.TRUE.,13)
C
      KFRSAV = KFREE
C
C     ************************
C     *** Density matrices ***
C     ************************
C
C
C     Get CI ref. vector from record 1 on LUBCI.
C
C
C     Transition densities
C     ---------------------
C
      CALL MEMGET('REAL',KDVT,NCSIM * N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVT,
     &     NCSIM * NASHT * NASHT * N2ASHX * NZ * 3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCREF,NZCONFQ,WORK,KFREE,LFREE)
!     IF ( NZCONF .GT. 0) CALL READAC(LUBCI,NZCONFQ,WORK(KCREF),1)
      IF ( NZCONF .GT. 0) CALL READ_DX(LUBCI,1,NZCONFQ,WORK(KCREF))
C
      IF ( NCSIM .GT. 0 ) THEN
         CALL XRSTDM(WORK(KCREF),BCVEC,WORK(KDVT),WORK(KPVT),
     &        WORK,KFREE,LFREE)
      END IF
      CALL MEMREL('after XRSTDM',WORK,1,KCREF,KFREE,LFREE)
C
C     Wf. densities
C     -------------
C
      CALL MEMGET('REAL', KDV, N2ASHXQ, WORK, KFREE, LFREE)
      IF ( NESIM + NPSIM .GT. 0 .AND. NASHT .GT. 0) THEN
         CALL GETMAT('DVMO',WORK(KDV),IPRXRS,FKRMC,WORK(KFREE),LFREE)
      END IF
C
C
C
C     ***************************************************
C     *** Generate Fock matrices (FVT, FCT, FVX, FCX) ***
C     ***************************************************
C
C
C     Allocate memory for FCX, FVT, and FVX
C     -------------------------------------
C
C     (allocate matrices with dimension N2BBASXQ so they
C     can be used as AO matrices in Fock matrix construction).
C
C     The matrices are allocated consequtive, so they can
C     be used as density matrices in the call to TWOFCK.
C
C     NCSIM FVT's + NOSIM FCX's + NOSIM FVX's
C
C
C     Calculate two-electron part of FCX, FVT, and FVX
C     ------------------------------------------------
C
      NDIM = NCSIM + NOSIM
      IF (NASHT .GT. 0) NDIM = NDIM + NOSIM
      CALL MEMGET2('REAL','F',KF,NDIM*N2BBASXQ,WORK,KFREE,LFREE)
C
      IF (LSVCFG(1)) THEN
        CALL GMOLITX(WORK(KF),WORK(KDVT),WORK(KDV),
     &     BEVEC,BPVEC,IBCVC,IBEVC,IBPVC,IBTYP,NPOS,WORK,KFREE,LFREE)
      ELSE
        CALL DZERO(WORK(KF),N2BBASXQ*NDIM)
      ENDIF
C
C
C     Add one-electron part to FCX and FVX
C     ------------------------------------
C
      IF ( NOSIM .GT. 0 .AND. LSVCFG(2)) THEN
         CALL FMOLI(WORK(KF + NCSIM * N2ORBXQ),BEVEC,BPVEC,
     &        IBEVC,IBPVC,IBTYP,WORK,KFREE,LFREE)
      END IF
C
C
C     ********************************
C     *** Generate FQT, FQX, H2XAC ***
C     ********************************
C
      IF (NASHT.GT.1) THEN
         CALL MEMGET('REAL',KFQT,NCSIM * NASHT * NORBT * NZ,
     &        WORK,KFREE,LFREE)

         CALL MEMGET('REAL',KFQX,NOSIM * NASHT * NORBT * NZ,
     &        WORK,KFREE,LFREE)

         CALL MEMGET('REAL',KH2ACX,
     &        NOSIM * NASHT * NASHT * NNASHX * NZ * 3,
     &        WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KPV,
     &        NASHT * NASHT * N2ASHX * NZ * 3,WORK,KFREE,LFREE)
         CALL GETMAT('PVMO',WORK(KPV),IPRXRS,FKRMC,WORK(KFREE),LFREE)
C
         CALL MEMGET('INTE',kibeig,norbt,WORK,KFREE,LFREE)
         call izero(work(kibeig),norbt)
C        read orbital symmetry info
         if(spinfr.or.levyle)then
           call ireakrmc(lukrmc,'IBEIG   ',work(kibeig),norbt)
         else if(linear)then
           call ireakrmc(lukrmc,'MJVEC   ',work(kibeig),norbt)
         end if
C
         CALL RGETH2TX(WORK(KFQX),WORK(KFQT),WORK(KH2ACX),
     &        work(kibeig),WORK(KPV),WORK(KPVT),
     &        BEVEC,BPVEC,
     &        get_orbital_rotation_indices_pp(),
     &        get_orbital_rotation_indices_pn(),
     &        IBEVC,IBPVC,IBTYP,IPRXRS,
     &        WORK(KFREE),LFREE)
C
C        Approximate positronic part of FQX or all of FQX:
C
         IF (NOSIM.GT.0 .AND. (XRS_NOPFQ .OR. XRS_NOFQX) ) THEN
            CALL XRSAPPFQ(WORK(KFQX),WORK(KF+NCSIM*N2ORBXQ),WORK(KDV),
     &           WORK(KFREE),LFREE)
         END IF
C
C        Dispose PV
C
         CALL MEMREL('after RGETH2TX',WORK,1,KPV,KFREE,LFREE)
      ELSE
         KFQT   = KFREE
         KFQX   = KFREE
         KH2ACX = KFREE
      END IF
C
C
C     *******************************************************
C     *** Calculate orbital part of conf. sigma vector(s) ***
C     *******************************************************
C
C
      IF ( NCSIM .GT. 0 ) THEN
         CALL RSIGOC(SEVEC,SPVEC,WORK(KDVT),WORK(KF),WORK(KFQT),
     &        get_orbital_rotation_indices_pp(),NZXOPE,
     &        get_orbital_rotation_indices_pn(),NZXOPP,
     &        IBCVC,IBTYP,WORK(KFREE),LFREE)
      END IF
C
C
C     *********************************************************
C     *** Calculate orbital part of orbital sigma vector(s) ***
C     *********************************************************
C
C
      IF ( NOSIM .GT. 0 ) THEN
         CALL RSIGOO(SEVEC,SPVEC,WORK(KF+NCSIM * N2ORBXQ),
     &        WORK(KFQX),WORK(KDV),
     &        get_orbital_rotation_indices_pp(),NZXOPE,
     &        get_orbital_rotation_indices_pn(),NZXOPP,
     &        IBEVC,IBPVC,IBTYP,WORK(KFREE),LFREE)
         IF (FKRMC) THEN
            CALL RSIGOOGB(SEVEC,SPVEC,
     &           GPOE,GPOP,BEVEC,BPVEC,
     &           get_orbital_rotation_indices_pp(),
     &           get_orbital_rotation_indices_pn(),
     &           WORK(KFREE),LFREE)
         END IF
      END IF
C
C
C     *******************************************************
C     *** Calculate conf. part of orbital sigma vector(s) ***
C     *******************************************************
C
C
      IF ( NOSIM .GT. 0 .AND. NZCONF .GT. 0) THEN
         CALL RSIGCO(SCVEC,IBCVC,WORK(KF),WORK(KH2ACX),
     &        WORK,KFREE,LFREE)
      END IF
C
      CALL RTKTIME(.FALSE.,13)
C
      CALL MEMREL('XRSSEP',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT ('XRSSEP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dmoli3 */
      SUBROUTINE DMOLI3(ITYP,IH,DMAT,CMO,BMAT)
C***********************************************************************
C
C     Generate modified density matrix for use in generation of
C     sigma vectors
C
C     Bit 1 of ITYP = 0   e-e  rotations in BMAT
C                     1   e-p rotations in BMAT
C     Bit 2 of ITYP = 0   inactive density matrix
C                     1   active density matrix
C
C     We define: DB_{ab} = c_{as}B_{si}c_{bi}^*
C                a,b - AO-indices
C                s   - all indices (electronic or positronic)
C                i   - inactive index
C
C     IH   = 1        Hermitian density matrix:  DMAT = -[DB + DB^+]
C                       based on Hermitian BMAT
C     IH   =-1   anti-Hermitian density matrix:  DMAT = -[DB - DB^+]
C                       based on anti-Hermitian BMAT
C
C     For time-symmetric operators DB^+ refers to the Hermite conjugate
C     For time-antisymmetric operators DB^+ refers to the Hermite
C     conjugate of the i-transform
C
C     Written by T.Saue Sep 17 1996
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, DP5 = 0.5D00)
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxrs.h"
C
      DIMENSION DMAT(N2BBASX,NZ),CMO(*),BMAT(NORBT,NORBT,NZ)
      LOGICAL   LBIT, LBIT_1
C
C
C     Initialization
C     ==============
C
      IREP = JSYMOP-1
      CALL DZERO(DMAT,N2BBASXQ)
C
C     Print section
C
      IF(IPRXRS.GE.8) THEN
        WRITE(LUPRI,'(A)') 'DMOLI3: Input B-matrix:'
        CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                IPQTOQ(1,IREP),LUPRI)
      ENDIF
C
C     Generate DB
C     ============
C
!     ITYP_TEMP = ITYP
!     ITYP = ITYP_TEMP

      DO 10 I1 = 1, NFSYM
        I2  = MOD(I1+JOPSY,2) + 1
C
C       Parameters for first index
C    
        IF (.NOT. LBIT(ITYP,1)) THEN
C
C         electronic indices...
C
          NS = NESH(I1)
          IS = IORB(I1) + NPSH(I1) + 1
          JS = ICMOQ(I1) + NPSH(I1)*NFBAS(I1,0) + 1
        ELSE
C
C         or positronic indices...
C
          NS = NPSH(I1)
          IS = IORB(I1) + 1
          JS = ICMOQ(I1) + 1
        ENDIF
C
C       Parameters for inactive index
C
         IF (.NOT. LBIT(ITYP,2)) THEN
C           write(lupri,*) 'inactive DMOLI3'
C
C           inactive density matrix...
C
 
          ! hjaaj Feb 2017: do not exclude frozen orbitals here:
          !   simpler code, and does not matter for performance.
          ! NI  = NISH(I2) - NFRO(I2)
          ! II  = IORB(I2) + NPSH(I2) + NFRO(I2) + 1
          ! JI  = ICMOQ(I2) + (NPSH(I2)+ NFRO(I2))*NFBAS(I2,0) + 1
            NI  = NISH(I2)
            II  = IORB(I2) + NPSH(I2) + 1
            JI  = ICMOQ(I2) + NPSH(I2)*NFBAS(I2,0) + 1
C
C          Generate modified density matrix as transformation MO to AO
C
           IF(NI.EQ.0.OR.NS.EQ.0) GOTO 10
CMI         bug trap
            ITYP_TEMP = ITYP
            CALL QTRANS90('MOAO','S',D0,
     &                  NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &                  DMAT(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),
     &                  NZ,IPQTOQ(1,IREP),
     &                  BMAT(IS,II,1),NORBT,NORBT,
     &                  NZ,IPQTOQ(1,IREP),
     &                  CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &                  CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &                  IPRXRS)
CMI        'empirical'  bug trap for gfortran in OpenMPI, ILP64, optimized 
           ITYP_TEMP = ITYP
C
         ELSE
C           write(lupri,*) 'active DMOLI3'
C
C           or active density matrix...
C
            NI  = NASH(I2)
            II  = IORB(I2) + NPSH(I2) + NISH(I2) + 1
            JI  = ICMOQ(I2) +
     &               (NPSH(I2)+ NISH(I2))*NFBAS(I2,0) + 1
C
C
C          Generate modified density matrix as transformation MO to AO
C
           IF(NI.EQ.0.OR.NS.EQ.0) GOTO 10
            CALL QTRANS90('MOAO','S',D0,
     &                  NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &                  DMAT(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),
     &                  NZ,IPQTOQ(1,IREP),
     &                  BMAT(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &                  CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,IREP),
     &                  CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,IREP),
     &                  IPRXRS)
C
C
         END IF
 10   CONTINUE
C
      IF (LBIT(ITYP,2)) THEN
C
C        For active density matrices scale the density matrix by
C        1/2 because all Fock matrices are multiplied with 2 in TWOFCK.
C
         CALL DSCAL(N2BBASXQ,DP5,DMAT,1)
      END IF
C
      IF(IPRXRS.GE.8) THEN
        CALL HEADER(
     &     'DMOLI3: Non-symmetrized AO modified density matrix',-1)
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
C
C
C     Take symmetric/antisymmetric combinations to form DMAT
C     ======================================================
C
      DO IZ = 1,NZ
        IQ = IPQTOQ(IZ,IREP)
        IM = IHQMAT(IQ,IH)
        IF(IM.EQ.1) THEN
          CALL FULMAT('S',NTBAS(0),NTBAS(0),DMAT(1,IZ))
        ELSEIF(IM.EQ.2) THEN
          CALL FULMAT('A',NTBAS(0),NTBAS(0),DMAT(1,IZ))
        ENDIF
      ENDDO
C
C     Print section
C
      IF(IPRXRS.GE.8) THEN
        CALL HEADER('DMOLI3: Symmetrized AO modified density matrix',-1)
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Fmoli */
      SUBROUTINE FMOLI(FMO,BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &                 WORK,KFREE,LFREE)
C***********************************************************************
C
C     FMO contains G-matrix; add one-index transformed Fock-matrix
C
C     Written by T.Saue Sep 17 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbxrs.h"
      DIMENSION FMO(*),BEVEC(*),BPVEC(*),IBTYP(*),
     &          IBEVC(*),IBPVC(*),WORK(*)
C
      CALL QENTER('FMOLI')
      CALL RTKTIME(.TRUE.,6)
      KFRSAV = KFREE
      CALL MEMGET('REAL',KBMAT,N2BBASXQ,WORK,KFREE,LFREE)
      IFAC = 1
      IF ( NASHT .GT. 0) IFAC = 2
      CALL MEMGET('REAL',KFBUF,N2BBASXQ*IFAC,WORK,KFREE,LFREE)
C
C     Get Fock matrices in MO-basis (FCmo and FVmo)
C     =============================================
C
      IF (NASHT .GT. 0) THEN
         CALL GETMAT('FCMO+FVMO',WORK(KFBUF),IPRXRS,FKRMC,
     $        WORK(KFREE),LFREE)
      ELSE
         CALL GETMAT('FCMO',WORK(KFBUF),IPRXRS,FKRMC,WORK(KFREE),LFREE)
      END IF
c     CALL GETFCK(WORK(KFBUF),IPRXRS,WORK,KFREE,LFREE)
C
C     Do one-index transformation
C
      CALL FMOLI2(FMO,WORK(KFBUF),WORK(KBMAT),BEVEC,BPVEC,
     &            IBEVC,IBPVC,IBTYP,WORK,KFREE,LFREE)
C
      CALL MEMREL('FMOLI',WORK,1,KFRSAV,KFREE,LFREE)
      CALL RTKTIME(.FALSE.,6)
      CALL QEXIT('FMOLI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Fmoli2 */
      SUBROUTINE FMOLI2(FMO,FBUF,BMAT,BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     FMO contains G-matrix; add one-index transformed Fock-matrix
C     FBUF contains original Fock-matrix
C
C     Written by T.Saue Sep 17 1996
C
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION FMO(N2ORBXQ,*),FBUF(*),BMAT(*),BEVEC(*),
     &          BPVEC(*),IBEVC(*),IBPVC(*),IBTYP(2,*),WORK(*)
C
      DIMENSION JOCCMO(2), NOCCMO(2), JGENMO(2), NGENMO(2)
C
      DO I = 1, NFSYM
         JOCCMO(I) = IORB(I) + 1 + NPSH(I)
         NOCCMO(I) = NOCC(I)
         JGENMO(I) = IORB(I) + 1
         NGENMO(I) = NORB(I)
      END DO
      IM = 0
      IVEC = 0
      IREP = JSYMOP-1
C
C
C     *************************************************************
C     *** Add one-index tranformed Fock-matrices; e-e rotations ***
C     *************************************************************
C
C
      IBOFF = 1
      DO I = 1,NESIM
C
         IVEC = IVEC + 1
         IH   = IBEVC(NERED+IVEC)
         IH = IBTYP(2,IH)
C
C        Inactive part (FC)
C        ------------------
C
C        There is always an inactive part.
C        If there is no inactive orbitals then FMO has been DZERO'ed earlier
C        so we can safely add the FC B - B FC part.
C
         ITYP  = 0*1 + 0*2
         IM = IM + 1
C
C        Scatter BVEC into BMAT
C        ----------------------
C
         CALL DZERO(BMAT,N2ORBXQ)
         CALL XOPSCT(BEVEC(IBOFF),get_orbital_rotation_indices_pp(),
     &               NZXOPE,BMAT,NORBT,NZ)
C
C        Construct full BMAT
C        -------------------
C
         IOFF = 1
         DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            IS = IHQMAT(IQ,IH)
            IF(IS.EQ.1) THEN
               CALL FULMAT('A',NORBT,NORBT,BMAT(IOFF))
            ELSEIF(IS.EQ.2) THEN
               CALL FULMAT('S',NORBT,NORBT,BMAT(IOFF))
            ENDIF
            IOFF = IOFF + N2ORBX
         ENDDO
C
         IF(IPRXRS.GE.6) THEN
            WRITE(LUPRI,'(//A,I5)')
     &           '* FMOLI: Scattered (e-e) trial vector no.:',I
            CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,'(/A,I5)')
     &           '* FMOLI: QXC no.:',I
            CALL PRQMAT(FMO(1,IM),NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,'(/A,I5)')
     &           '* FMOLI: FC:',I
            CALL PRQMAT(FBUF,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,IREP),LUPRI)
         ENDIF
C
C
C        Add one-index transformed Fock-matrix
C        -------------------------------------
C
C
         CALL FMOLI3(ITYP,FMO(1,IM),BMAT,FBUF,
     &        WORK,KFREE,LFREE)
C
         IF (NASHT .GT. 0) THEN
C
C           active Fock matrices (FV)
C           -------------------------
C
C
            ITYP  = 0*1 + 0*2
            IM = IM + 1
C
C           Scatter BVEC into BMAT
C           ----------------------
C
            CALL DZERO(BMAT,N2ORBXQ)
            CALL XOPSCT(BEVEC(IBOFF),get_orbital_rotation_indices_pp(),
     &                  NZXOPE,BMAT,NORBT,NZ)
C
C           Construct full BMAT
C           -------------------
C
            IOFF = 1
            DO IZ = 1,NZ
               IQ = IPQTOQ(IZ,IREP)
               IS = IHQMAT(IQ,IH)
               IF(IS.EQ.1) THEN
                  CALL FULMAT('A',NORBT,NORBT,BMAT(IOFF))
               ELSEIF(IS.EQ.2) THEN
                  CALL FULMAT('S',NORBT,NORBT,BMAT(IOFF))
               ENDIF
               IOFF = IOFF + N2ORBX
            ENDDO
C
            IF(IPRXRS.GE.6) THEN
               WRITE(LUPRI,'(//A,I5)')
     &              '* FMOLI: Scattered (e-e) trial vector no.:',
     &              I
               CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,'(/A,I5)')
     &              '* FMOLI: FVX no.:',I
               CALL PRQMAT(FMO(1,IM),NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,'(/A,I5)')
     &              '* FMOLI: FV:',I
               CALL PRQMAT(FBUF(1+N2ORBXQ),NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
            ENDIF
C
C           Add one-index transformed Fock-matrix
C           -------------------------------------
C
            CALL FMOLI3(ITYP,FMO(1,IM),BMAT,FBUF(1+N2ORBXQ),
     &           WORK,KFREE,LFREE)
C
         END IF
         IBOFF = IBOFF + NZXOPEQ
      ENDDO
C
C     *************************************************************
C     *** Add one-index tranformed Fock-matrices; e-p rotations ***
C     *************************************************************
C
      IBOFF = 1
      IVEC = 0
      DO I = 1,NPSIM
C
         IVEC = IVEC + 1
         IH   = IBPVC(NPRED+IVEC)
         IH = IBTYP(2,IH)
C
C        Inactive
C        --------
C
         ITYP  = 1*1 + 0*2
         IM = IM + 1
C
C        Scatter BVEC into BMAT
C        ----------------------
C
         CALL DZERO(BMAT,N2ORBXQ)
         CALL XOPSCT(BPVEC(IBOFF),get_orbital_rotation_indices_pn(),
     &               NZXOPP,BMAT,NORBT,NZ)
C
C        Construct full BMAT
C        -------------------
C
         IOFF = 1
         DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            IS = IHQMAT(IQ,IH)
            IF(IS.EQ.1) THEN
               CALL FULMAT('A',NORBT,NORBT,BMAT(IOFF))
            ELSEIF(IS.EQ.2) THEN
               CALL FULMAT('S',NORBT,NORBT,BMAT(IOFF))
            ENDIF
            IOFF = IOFF + N2ORBX
         ENDDO
C
         IF(IPRXRS.GE.6) THEN
            WRITE(LUPRI,'(//A,I5)')
     &           '* FMOLI: Scattered inactive (e-p) trial vector no.:',I
            CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,'(/A,I5)')
     &           '* FMOLI: QCX no.:',I
            CALL PRQMAT(FMO(1,IM),NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,'(/A,I5)')
     &           '* FMOLI: FC:',I
            CALL PRQMAT(FBUF,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,IREP),LUPRI)
         ENDIF
C
C        Add one-index transformed Fock-matrix
C        -------------------------------------
C
         CALL FMOLI3(ITYP,FMO(1,IM),BMAT,FBUF,
     &               WORK,KFREE,LFREE)
C
         IF (NASHT .GT. 0) THEN
C
C           Active
C           ------
C
            ITYP  = 1*1 + 0*2
            IM = IM + 1
C
C           Scatter BVEC into BMAT
C           ----------------------
C
            CALL DZERO(BMAT,N2ORBXQ)
            CALL XOPSCT(BPVEC(IBOFF),get_orbital_rotation_indices_pn(),
     &                  NZXOPP,BMAT,NORBT,NZ)
C
C           Construct full BMAT
C           -------------------
C
            IOFF = 1
            DO IZ = 1,NZ
               IQ = IPQTOQ(IZ,IREP)
               IS = IHQMAT(IQ,IH)
               IF(IS.EQ.1) THEN
                  CALL FULMAT('A',NORBT,NORBT,BMAT(IOFF))
               ELSEIF(IS.EQ.2) THEN
                  CALL FULMAT('S',NORBT,NORBT,BMAT(IOFF))
               ENDIF
               IOFF = IOFF + N2ORBX
            ENDDO
C
            IF(IPRXRS.GE.6) THEN
               WRITE(LUPRI,'(A,I5)')
     &              '* FMOLI: Scattered active (e-p) trial vector no.:',
     &              I
               CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,'(/A,I5)')
     &              '* FMOLI: FVX no.:',I
               CALL PRQMAT(FMO(1,IM),NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,'(/A,I5)')
     &              '* FMOLI: FV:',I
               CALL PRQMAT(FBUF(1+N2ORBXQ),NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
            ENDIF
C
C           Add one-index transformed Fock-matrix
C           -------------------------------------
C
            CALL FMOLI3(ITYP,FMO(1,IM),BMAT,FBUF(1+N2ORBXQ),
     &           WORK,KFREE,LFREE)
C
         END IF
         IBOFF = IBOFF + NZXOPPQ
      ENDDO
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rtr1h1 */
      SUBROUTINE RTR1H1(J1MO,N1MO,J2MO,N2MO,
     &     H1X,IPH1X,BMAT,IPBMAT,H1,IPH1,IH1SYM,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate one-index transformed Fock-matrices
C     for use in generation of sigma vectors
C
C     H1X_{pq} = H1X_{pq} +
C        \sum_r BMAT_{pr} H1_{rq} - H1_{pr} BMAT{rq}
C
C     J{1,2}MO is the start orbital for the transformation (p & q above)
C     N{1,2}MO is the number of orbitals to transform
C     IH1SYM is the symmetry of H1.
C
C     Written by J. Thyssen - Dec 6 2000
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0,DM1 = -1.0D0)
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxrs.h"
C
      DIMENSION H1X(NORBT,NORBT,NZ)
      DIMENSION H1(NORBT,NORBT,NZ)
      DIMENSION BMAT(NORBT,NORBT,NZ)
      DIMENSION J1MO(2), J2MO(2), N1MO(2), N2MO(2)
      DIMENSION WORK(*)
      DIMENSION IPH1X(*), IPH1(*), IPBMAT(*)
C
C
      CALL RTKTIME(.TRUE.,15)
C
C     Generate one-index transformed Fock matrix
C     ------------------------------------------
C
      DO 10 I1 = 1,NFSYM
         IRSYM = MOD(I1+JOPSY,2) + 1
         I2 = MOD(IRSYM+IH1SYM,2) + 1
C
         IF (N1MO(I1) .NE. 0 .AND. N2MO(I2) .NE. 0) THEN
C
C           First part: H1X_{pq} = H1X_{pq} + BMAT_{pr}H1_{rq}
C
            CALL QGEMM(N1MO(I1),N2MO(I2),NORB(IRSYM),D1,
     &           'N','N',IPBMAT,
     &           BMAT(J1MO(I1),1+IORB(IRSYM),1),NORBT,NORBT,NZ,
     &           'N','N',IPH1,
     &           H1(1+IORB(IRSYM),J2MO(I2),1),NORBT,NORBT,NZ,
     &           D1,IPH1X,H1X(J1MO(I1),J2MO(I2),1),
     &           NORBT,NORBT,NZ)
         END IF
C
C
C        Second part: H1X_{pq} = H1X_{pq} - H1^t_{pr}BMAT_{rq}
C
         IRSYM = MOD(I1+IH1SYM,2) + 1
         I2 = MOD(IRSYM+JOPSY,2) + 1
         IF (N1MO(I1) .NE. 0 .AND. N2MO(I2) .NE. 0) THEN
            CALL QGEMM(N1MO(I1),N2MO(I2),NORB(IRSYM),DM1,
     &           'N','N',IPH1,
     &           H1(J1MO(I1),1+IORB(IRSYM),1),NORBT,NORBT,NZ,
     &           'N','N',IPBMAT,
     &           BMAT(1+IORB(IRSYM),J2MO(I2),1),NORBT,NORBT,NZ,
     &           D1,IPH1X,H1X(J1MO(I1),J2MO(I2),1),
     &           NORBT,NORBT,NZ)
         END IF
C
 10   CONTINUE
C
C     Print section
C
      IF(IPRXRS.GE.30) THEN
         CALL HEADER('RTR1H1: one-index transformed matrix',-1)
         CALL PRQMAT(H1X,NORBT,NORBT,NORBT,
     &        NORBT,NZ,IPH1X,LUPRI)
      ENDIF
      CALL RTKTIME(.FALSE.,15)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Gmoli3 */
      SUBROUTINE GMOLI3(ITYP,FMO,CMO,FAO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Transform the G-matrix to MO-basis, but only the blocks needed
C
C     Bit 1 of ITYP = 0   e-e  rotations in BMAT
C                     1   e-p rotations in BMAT
C     Bit 2 of ITYP = 0   inactive density matrix
C                     1   active density matrix
C
C     Written by T.Saue Sep 18 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbxrs.h"
#include "dgroup.h"
      LOGICAL EPMIX,LBIT
      DIMENSION FMO(NORBT,NORBT,NZ),CMO(*),
     &          FAO(N2BBASXQ),WORK(*)
C
C
C     Transform only blocks of right symmetry combination
C     ===================================================
C
      IF(.NOT. LBIT(ITYP,1))     THEN
        EPMIX = NZXOPP.GT.0
      ELSE
        EPMIX = NZXOPE.GT.0
      ENDIF
      IREP  = JSYMOP-1
      DO 10 I1 = 1,NFSYM
        I2  = MOD(I1+JOPSY,2) + 1
C
C       Parameters for first index
C
        IF (.NOT. LBIT(ITYP,1)) THEN
C
C         electronic indices...
C
          NS = NESH(I1)
          IS = IORB(I1) + NPSH(I1) + 1
          JS = ICMOQ(I1) + NPSH(I1)*NFBAS(I1,0) + 1
        ELSE
C
C         or positronic indices...
C
          NS = NPSH(I1)
          IS = IORB(I1) + 1
          JS = ICMOQ(I1) + 1
        ENDIF
C
C       Parameters for second index.
C
         NI = NOCC(I2)
         II = IORB(I2) + NPSH(I2) + 1
         JI = ICMOQ(I2) + NPSH(I2) * NFBAS(I2,0) + 1
C
        IF(NS.NE.0.AND.NI.NE.0) THEN
          CALL QTRANS('AOMO','S',D0,
     &                NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &                FAO(I2BASX(I1,I2)+1),NTBAS(0),NTBAS(0),NZ,
     &                IPQTOQ(1,IREP),
     &                FMO(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &                CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,IPRXRS)
        ENDIF
C
C       Transform other block as well,
C       if both electronic and postitronic rotations exist
C
        IF(EPMIX) THEN
C
C         Parameters for secondary index(note the switch...)
C
C
C       Parameters for first index
C
        IF (.NOT. LBIT(ITYP,1)) THEN
C
C         positronic indices...
C
          NS = NPSH(I1)
          IS = IORB(I1) + 1
          JS = ICMOQ(I1) + 1
        ELSE
C
C         or electronic indices...
C
          NS = NESH(I1)
          IS = IORB(I1) + NPSH(I1) + 1
          JS = ICMOQ(I1) + NPSH(I1)*NFBAS(I1,0) + 1
        ENDIF
          IF(NS.NE.0.AND.NI.NE.0) THEN
            CALL QTRANS('AOMO','S',D0,
     &                NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &                FAO(I2BASX(I1,I2)+1),NTBAS(0),NTBAS(0),NZ,
     &                IPQTOQ(1,IREP),
     &                FMO(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &                CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,IPRXRS)
          ENDIF
        ENDIF
 10   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck setgmt */
      SUBROUTINE SETGMT(ISYMOP,IHRMOP,IFCKOP,IBCVC,IBEVC,IBPVC,IBTYP,
     &                  IPRINT)
C***********************************************************************
C
C     Set arrays IHRMOP, IFCKOP, ISYMOP in preparation
C     for construction of G-matrices
C
C     Written by T.Saue Oct 15 1996
C
C***********************************************************************
       use dirac_cfg
       use dft_cfg

#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dgroup.h"
#include "dcborb.h"
C
      DIMENSION IHRMOP(*),IFCKOP(*),ISYMOP(*)
      DIMENSION IBCVC(*),IBEVC(*),IBPVC(*),IBTYP(2,*)
C
      IM = 0
C
      NFMAT = 0
      IF (NISHT .GT. 0) NFMAT = NFMAT + 1
      IF (NASHT .GT. 0) NFMAT = NFMAT + 1
C
C     conf.
C     -----
C
C     Only FVT Fock matrices, so no loop over NFMAT.
C
      DO I = 1, NCSIM
         IM = IM + 1
         IS = IBCVC(NCRED+I)
         IHRMOP(IM) = IBTYP(2,IS)
         IFCKOP(IM) = 1
!        turn off exchange for .ALDA
!        unless user wants .XALDA
         if (dirac_cfg_dft_calculation .and.
     &       (dft_cfg_alda_hs .or. dft_cfg_alda_ha)
     &       .and. .not. dft_cfg_xalda) ifckop(im) = 2
         ISYMOP(IM) = JSYMOP
      END DO
C
C     e-e rotations
C     -------------
C
      DO I = 1,NESIM
         DO IMAT = 1,NFMAT
           IM = IM + 1
           IS = IBEVC(NERED+I)
           IHRMOP(IM) = IBTYP(2,IS)
           IFCKOP(IM) = 1
!        turn off exchange for .ALDA
!        unless user wants .XALDA
         if (dirac_cfg_dft_calculation .and.
     &       (dft_cfg_alda_hs .or. dft_cfg_alda_ha)
     &       .and. .not. dft_cfg_xalda) ifckop(im) = 2
           ISYMOP(IM) = JSYMOP
         ENDDO
      END DO
C
C
C
C     e-p rotations
C     -------------
C
      DO I = 1,NPSIM
         DO IMAT = 1,NFMAT
           IM = IM + 1
           IS = IBPVC(NPRED+I)
           IHRMOP(IM) = IBTYP(2,IS)
           IFCKOP(IM) = 1
!        turn off exchange for .ALDA
!        unless user wants .XALDA
         if (dirac_cfg_dft_calculation .and.
     &       (dft_cfg_alda_hs .or. dft_cfg_alda_ha)
     &       .and. .not. dft_cfg_xalda) ifckop(im) = 2
           ISYMOP(IM) = JSYMOP
         ENDDO
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrstvc */
      SUBROUTINE XRSTVC(BEVEC,IBEVC,NBEVC,IEOFF,
     &                  BPVEC,IBPVC,NBPVC,IPOFF,
     &                  BCVEC,IBCVC,NBCVC,ICOFF,
     &                  TEVEC,TPVEC,TCVEC,NBTOT,ITOFF,
     &                  IBTYP,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate tau vectors : metric S[2] times trial vector B(I)
C
C        TVEC(I) = S[2]*B(I)
C
C     The trial vectors come in three types:
C        C-type: contains only configurational parameters
C        E-type: contains only electron-eelectron (e-e) rotations
C        P-type: contains only electron-positron  (e-p) rotations
C
C     The tau vectors are correspondingly split up into a C, E and P -part.
C     The tau vectors are accumulated, so that one tau vector may
C     receive contributions from more than one trial vector.
C
C     Written by H.J.Aa.Jensen and Trond Saue 1996
C     Last revision Sep 19 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
C
      DIMENSION BCVEC(*),BEVEC(*),BPVEC(*),
     &          TCVEC(*),TEVEC(*),TPVEC(*),
     &          IBCVC(*),IBEVC(*),IBPVC(*),
     &          IBTYP(*),
     &          WORK(*)
C
#include "dcbxrs.h"
C
C
      CALL QENTER('XRSTVC')
C
C     Initialize tau vectors
C     ========================
C
      IF(NZCONF.GT.0) CALL DZERO(TCVEC,NBTOT*NZCONFQ)
      IF(NZXOPE.GT.0) CALL DZERO(TEVEC,NBTOT*NZXOPEQ)
      IF(NZXOPP.GT.0) CALL DZERO(TPVEC,NBTOT*NZXOPPQ)
C
C     Construct orbital part of tau vectors
C     =====================================
C
      IF (NESIM + NPSIM .GT. 0) THEN
          CALL XRSTEP(BEVEC,IBEVC,NBEVC,IEOFF,
     &                BPVEC,IBPVC,NBPVC,IPOFF,
     &                BCVEC,IBCVC,NBCVC,ICOFF,
     &                TEVEC,TPVEC,NBTOT,ITOFF,
     &                IBTYP,WORK,KFREE,LFREE)
       END IF
C
C     Construct configurational part of linear transformed vectors
C     for the configurational trial vectors
C     ============================================================
C
      IF (NCSIM.GT.0) THEN
         CALL XRSTCI(BCVEC,TCVEC,TEVEC,TPVEC,WORK,KFREE,LFREE)
      END IF
C
      CALL QEXIT('XRSTVC')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrstci */
      SUBROUTINE XRSTCI(BCVEC,TCVEC,TEVEC,TPVEC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate cc block of tau vector
C
C     Written by H.J.Aa.Jensen and T.Saue 1996
C     Last revision Oct 15 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION BCVEC(*),TCVEC(*),TEVEC(*), TPVEC(*),WORK(*)
C
C Used from common blocks:
C     XCBPRI: LUPRI
C     WRKXRS: NCSIM
C
#include "dcbxrs.h"
      CALL QENTER('XRSTCI')
      WRITE(LUPRI,'(A)')
     &     '*** Sorry, XRSTCI called, but is not implemented yet!'
      CALL QUIT('XRSCLI called, but is not implemented yet!')
      CALL QEXIT ('XRSTCI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrstov */
      SUBROUTINE XRSTEP(BEVEC,IBEVC,NBEVC,IEOFF,
     &                  BPVEC,IBPVC,NBPVC,IPOFF,
     &                  BCVEC,IBCVC,NBCVC,ICOFF,
     &                  TEVEC,TPVEC,NBTOT,ITOFF,
     &                  IBTYP,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Construct upper orbital part of tau vectors
C
C     TVEC_(mn) = S[2](o*o)_(mn,pq)U_(pq) + hS[2](o*o*)_(mn,pq)U^*_(pq)
C               = <0|[Nm,Pq]|0>U_(pq) - h<0|[Nm,Qp]|0>U^*_(pq)
C               = <0|[Nm,Pq]|0>W_(p,q)
C
C     where the matrix WA has the structure:
C
C           P  I   A   S             P  I  A  S
C          --------------         ---------------
C        P | 0  0   0   0       P |  0  U  U  0
C        I | 0  0 -hU+-hU+      I |-hU+ 0  0  0
C        A | 0  U   X -hU+      A |-hU+ 0  0  0
C        S | 0  U   U   0       S |  0  0  0  0
C
C           e-e rotations           e-p rotations
C
C     with the convention:
C        P - positronic spinors
C        I - inactive spinors
C        A - active spinors
C        S - secondary spinors
C
C     (X has -Y+ in upper diagonal and Z in lower diagonal)
C
C-----------------------------------------------------------------------
C     In RPA we need to evaluate
C
C     TVEC_(si) = W_(si)
C
C     Written by T.Saue and H.J.Aa.Jensen 1996
C     Last revision Oct 15 1996 - tsaue
C
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, DM1 = -1.0D00)
      DIMENSION BEVEC(NZXOPEQ,*),IBEVC(*),
     &          BPVEC(NZXOPPQ,*),IBPVC(*),
     &          BCVEC(NZCONFQ,*),IBCVC(*),
     &          TEVEC(NZXOPEQ,*),TPVEC(NZXOPPQ,*),
     &          IBTYP(2,*),WORK(*)
      CHARACTER SSTR(2)*1
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxrs.h"
C
      DATA SSTR /'S','A'/
C
C
C     Hermicity of operator
C     =====================
C
C
C     Allocate memory and get active density matrix
C     =============================================
C
      IREP = JSYMOP - 1
      KFRSAV = KFREE
      MNZXOP = MAX(NZXOPEQ,NZXOPPQ)
      CALL MEMGET('REAL',KUBVEC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KOP   ,MNZXOP ,WORK,KFREE,LFREE)
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDV,N2ASHXQ,WORK,KFREE,LFREE)
         CALL GETMAT('DVMO',WORK(KDV),IPRXRS,FKRMC,WORK(KFREE),LFREE)
      ELSE
         CALL MEMGET('REAL',KDV,0,WORK,KFREE,LFREE)
      END IF
C
C     Calculate tau vectors
C     =====================
C
C     Electronic...
C
      DO I = 1,NBEVC
         IH = IBEVC(I+IEOFF)
         IH = IBTYP(2,IH)
         IF (IH .EQ. -1) THEN
            IH  = 1
            IHH = 1
         ELSE IF (IH .EQ. 1) THEN
            IH  = -1
            IHH =  2
         ELSE
            CALL QUIT('argh ralle ralle')
         END IF
C
C        Construct full matrix
C        ---------------------
C
         CALL DZERO(WORK(KUBVEC),N2ORBXQ)
         CALL XOPSCT(BEVEC(1,I),get_orbital_rotation_indices_pp(),
     &               NZXOPE,
     &               WORK(KUBVEC),NORBT,NZ)

         IF (IPRXRS .GE. 100) THEN
            WRITE(LUPRI,'(1X,A,I3)')
     &         'XRSTEP: Scattered trial vector no. ',I
            CALL PRQMAT(WORK(KUBVEC),NORBT,NORBT,NORBT,NORBT,
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
         IZOFF = 0
         DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            IS = IHQMAT(IQ,IH)
            CALL FULMAT(SSTR(IS),NORBT,NORBT,WORK(KUBVEC+IZOFF))
            IZOFF = IZOFF + N2ORBX
         END DO
         IF (IPRXRS .GE. 100) THEN
            WRITE(LUPRI,'(1X,A,I3)')
     &         'XRSTEP: Full trial vector no. ',I
            CALL PRQMAT(WORK(KUBVEC),NORBT,NORBT,NORBT,NORBT,
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
C        Calculate contribution to tau vector
C        ------------------------------------
C
         CALL DZERO(WORK(KOP),NZXOPEQ)
         CALL PRPORB(WORK(KOP),get_orbital_rotation_indices_pp(),
     &               NZXOPE,WORK(KUBVEC),
     &               NORBT,WORK(KDV),NZ,IHH,IPRXRS)
C
C
         J  = IBEVC(I+IEOFF) - ITOFF
         IF (IPRXRS .GE. 20) THEN
             WRITE(LUPRI,'(A,I2)') 'XRSTEP: contribution to J = ',J
             WRITE(LUPRI,'(I5,3X,1P,E18.10)')
     &         (K,WORK(KOP+K-1), K = 1,NZXOPE)
         END IF
C
C
         CALL DAXPY(NZXOPEQ,D1,WORK(KOP),1,TEVEC(1,J),1)
C
      END DO
C
C
C     Positronic...
C
      DO I = 1,NBPVC
         IH = IBPVC(I+IPOFF)
         IH = IBTYP(2,IH)
         IF (IH .EQ. -1) THEN
            IH  = 1
            IHH = 1
         ELSE IF (IH .EQ. 1) THEN
            IH  = -1
            IHH =  2
         ELSE
            CALL QUIT('argh ralle ralle')
         END IF
C
C        Construct full matrix
C        ---------------------
C
         CALL DZERO(WORK(KUBVEC),N2ORBXQ)
         CALL XOPSCT(BPVEC(1,I),get_orbital_rotation_indices_pn(),
     &               NZXOPP,
     &               WORK(KUBVEC),NORBT,NZ)

         IF (IPRXRS .GE. 100) THEN
            WRITE(LUPRI,'(1X,A,I3)')
     &         'XRSTEP: Scattered trial vector no. ',I
            CALL PRQMAT(WORK(KUBVEC),NORBT,NORBT,NORBT,NORBT,
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
         IZOFF = 0
         DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            IS = IHQMAT(IQ,IH)
            CALL FULMAT(SSTR(IS),NORBT,NORBT,WORK(KUBVEC+IZOFF))
            IZOFF = IZOFF + N2ORBX
         END DO
         IF (IPRXRS .GE. 100) THEN
            WRITE(LUPRI,'(1X,A,I3)')
     &         'XRSTEP: Full trial vector no. ',I
            CALL PRQMAT(WORK(KUBVEC),NORBT,NORBT,NORBT,NORBT,
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
C        Calculate contribution to tau vector
C        ------------------------------------
C
         CALL DZERO(WORK(KOP),NZXOPPQ)
         CALL PRPORB(WORK(KOP),get_orbital_rotation_indices_pn(),
     &               NZXOPP,WORK(KUBVEC),
     &               NORBT,WORK(KDV),NZ,IHH,IPRXRS)
C
C
         J  = IBPVC(I+IPOFF) - ITOFF
         IF (IPRXRS .GE. 20) THEN
             WRITE(LUPRI,'(A,I2)') 'XRSTEP: contribution to J = ',J
             WRITE(LUPRI,'(I5,3X,1P,E18.10)')
     &         (K,WORK(KOP+K-1), K = 1,NZXOPP)
         END IF
C
C
         CALL DAXPY(NZXOPPQ,D1,WORK(KOP),1,TPVEC(1,J),1)
C
      END DO
      CALL MEMREL('XRSTEP',WORK,1,KFRSAV,KFREE,LFREE)

Cnew code
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Xrsxvc */
      SUBROUTINE XRSXVC(ITYP,XEVEC,XPVEC,XCVEC,EVECR,MFREQ,
     &           IBTYP,IBCVC,IBEVC,IBPVC,BEBUF,BPBUF,BCBUF,
     &           WORK,KFREE,LFREE)
C***********************************************************************
C
C     Form solution vector in response equations.
C     The solution vector is split into a (+) and a (-) part.
C     In the static case the (-) part is automatically zero.
C
C     Written by T.Saue and H.J.Aa.Jensen Oct 16 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dcbibn.h"
      DIMENSION XEVEC(*),XPVEC(*),XCVEC(*),IBTYP(2,*),
     &          IBCVC(*),IBEVC(*),IBPVC(*),EVECR(*),
     &          BEBUF(*),BPBUF(*),BCBUF(*),WORK(*)
C
      CALL QENTER('XRSXVC')
      IF(NZXOPE.GT.0) CALL DZERO(XEVEC,MFREQ*NZXOPEQ)
      IF(NZXOPP.GT.0) CALL DZERO(XPVEC,MFREQ*NZXOPPQ)
      IF(NZCONF.GT.0) CALL DZERO(XCVEC,MFREQ*NZCONFQ)
C
C     e-e rotations
C
      IF(NERED.GT.0) THEN
        CALL XRSXV1(ITYP,JBENDX,XEVEC,EVECR,MFREQ,
     &              IBTYP,IBEVC,IVECS,BEBUF)
      ENDIF
C
C     e-p rotations
C
      IF(NPRED.GT.0) THEN
        CALL XRSXV1(ITYP,JBPNDX,XPVEC,EVECR,MFREQ,
     &              IBTYP,IBPVC,IVECS,BPBUF)
      ENDIF
C
C     Configurational parameters
C
      IF(NCRED.GT.0) THEN
        CALL XRSXV1(ITYP,JBCNDX,XCVEC,EVECR,MFREQ,
     &              IBTYP,IBCVC,IVECS,BCBUF)
      ENDIF
C
      CALL QEXIT('XRSXVC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Xrsxv1 */
      SUBROUTINE XRSXV1(ITYP,NBTYP,XVEC,EVECR,MFREQ,IBTYP,IBVEC,
     &                  IVECS,BVBUF)
C***********************************************************************
C
C     Form electronic/positronic/configurational
C     part of solution vector in linear response
C
C     Written by T.Saue Oct 16 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dcbibn.h"
      PARAMETER(DTOL=1.0D-10,D0=0.0D0)
C
      DIMENSION XVEC(*),EVECR(NEVEC,MFREQ),
     &          IBTYP(2,*),IBVEC(*),IVECS(*),BVBUF(*)
C
C
      IF    (NBTYP.EQ.JBCNDX) THEN
        NBDIM  = NZCONFQ
        NBRED  = NCRED
        LUBVEC = LUBCI
        IRECOFF = 1
      ELSEIF(NBTYP.EQ.JBENDX) THEN
        NBDIM  = NZXOPEQ
        NBRED  = NERED
        LUBVEC = LUBOE
        IRECOFF = 0
      ELSEIF(NBTYP.EQ.JBPNDX) THEN
        NBDIM  = NZXOPPQ
        NBRED  = NPRED
        LUBVEC = LUBOP
        IRECOFF = 0
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'XRSXVC ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('XRSXVC: Unknown NBTYP !')
      ENDIF
C
      CALL IZERO(IVECS,MFREQ)
      CALL DZERO(XVEC,NBDIM*MFREQ)
      IF(DAMPFREQ.EQ.D0) THEN
        DO IREC = 1,NBRED
          IND = IBVEC(IREC)
          IH  = IBTYP(2,IND)
          IF(IH.EQ.ITYP) THEN
            CALL READ_DX(LUBVEC,IREC+IRECOFF,NBDIM,BVBUF)
            IBOFF = 1
            DO IFREQ = 1,MFREQ
              IF(ABS(EVECR(IND,IFREQ)).GT.DTOL) THEN
                IVECS(IFREQ) = IVECS(IFREQ) + 1
                CALL DAXPY(NBDIM,EVECR(IND,IFREQ),
     &                       BVBUF,1,XVEC(IBOFF),1)
              ENDIF
              IBOFF = IBOFF + NBDIM
            ENDDO
          ENDIF
        ENDDO
      ELSE
        DO IREC = 1,NBRED
          IND = IBVEC(IREC)
          INDR = 2*IND - 1
          INDI = 2*IND
          IH  = IBTYP(2,IND)
          IF(IH.EQ.ITYP) THEN
            CALL READ_DX(LUBVEC,IREC+IRECOFF,NBDIM,BVBUF)
            IBOFF = 1
            DO IFREQ = 1,MFREQ
              IF(ABS(EVECR(INDR,IFREQ)).GT.DTOL) THEN
                IVECS(IFREQ) = IVECS(IFREQ) + 1
                CALL DAXPY(NBDIM,EVECR(INDR,IFREQ),
     &                       BVBUF,1,XVEC(IBOFF),1)
              ENDIF
              IBOFF = IBOFF + NBDIM
              IF(ABS(EVECR(INDI,IFREQ)).GT.DTOL) THEN
                IVECS(IFREQ) = IVECS(IFREQ) + 1
                CALL DAXPY(NBDIM,EVECR(INDI,IFREQ),
     &                       BVBUF,1,XVEC(IBOFF),1)
              ENDIF
              IBOFF = IBOFF + NBDIM
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ACCS2X */
      SUBROUTINE ACCS2X(BEPVC,BEMVC,BPPVC,BPMVC,BCPVC,BCMVC,
     &                  IBTYP,IBCVC,IBEVC,IBPVC,
     &                  EVECR,EVALR,MFREQ,BEBUF,BPBUF,BCBUF,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     In the construction of the residual accumulate -freq*S[2]X terms
C     for all frequencies W+Y*i.
C       W -> real part of frequency
C       Y -> imaginary part of frequency
C
C     Input:
C       B(+)-vec contains current X(-) solution vectors
C       B(-)-vec contains current X(+) solution vectors
C
C     Output:
C       For REAL frequencies:
C         B(+)-vec contains T(+) = -W*S[2]*X(-)
C         B(-)-vec contains T(-) = -W*S[2]*X(+)
C       For IMAGINARY frequencies:     
C         B(+)-vec contains T(+) = +Y*S[2]*X(-)
C         B(-)-vec contains T(-) = -Y*S[2]*X(+)
C       For COMPLEX frequencies:
C see JCP133,064105,2010 eq(21)
C /!\ eq(21) first line: the first "+" should be a "-"
C         B(+)-vec contains T(+) = -W*S[2]*X(-)  +Y*S[2]*X(-)'
C                      and  T(+)'= -W*S[2]*X(-)' -Y*S[2]*X(-)
C         B(-)-vec contains T(-) = -W*S[2]*X(+)  +Y*S[2]*X(+)'
C                      and  T(-)'= -W*S[2]*X(+)' -Y*S[2]*X(+)
C
C     Written by T.Saue Oct 17 1996
c     Updated Dec10 - Jan11 S.Villaume & T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dcbxrs.h"
#include "dcbibn.h"
      DIMENSION BEPVC(NZXOPEQ,*),BEMVC(NZXOPEQ,*),BEBUF(NZXOPEQ,NDAMP),
     &    BPPVC(max(1,NZXOPPQ),*),BPMVC(max(1,NZXOPPQ),*),
     &    BPBUF(max(1,NZXOPPQ),NDAMP),
     &    BCPVC(max(1,NZCONFQ),*),BCMVC(max(1,NZCONFQ),*),
     &    BCBUF(max(1,NZCONFQ)),
     &          IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),
     &          EVECR(NEVEC,MFREQ),EVALR(MFREQ),WORK(*)
C
      CALL QENTER('ACCS2X')
C
      NBEVC = MIN(NZXOPE,1)
      NBPVC = MIN(NZXOPP,1)
      NBCVC = MIN(NZCONF,1)
C
      II = 1
      DO I = 1,MFREQ
        FREQ = EVALR(I)
        IF(FREQ.NE.D0) THEN
C
C         Accumulate T(-) = S[2]X(+) in buffer
C         ------------------------------------
C
!radovan: bepvc, bppvc, bcpvc can go out of bounds
!         if corresponding degrees of freedom are zero
          CALL XRSTVC(BEPVC(1,II),1,NBEVC,0,
     &                BPPVC(1,II),1,NBPVC,0,
     &                BCPVC(1,I) ,1,NBCVC,0,
     &                BEBUF(1,1),BPBUF(1,1),BCBUF,1,0,
     &                IBTYP,WORK,KFREE,LFREE)
          IF(DAMPFREQ.NE.D0) THEN
            CALL XRSTVC(BEPVC(1,II+1),1,NBEVC,0,
     &                  BPPVC(1,II+1),1,NBPVC,0,
     &                  BCPVC(1,I)   ,1,NBCVC,0,
     &                  BEBUF(1,2),BPBUF(1,2),BCBUF,1,0,
     &                  IBTYP,WORK,KFREE,LFREE)
          ENDIF
C
C         Accumulate T(+) = S[2]X(-) in b(+)
C         ------------------------------------
C
!radovan: bepvc, bppvc, bcpvc can go out of bounds
!         if corresponding degrees of freedom are zero
          CALL XRSTVC(BEMVC(1,II),1,NBEVC,0,
     &                BPMVC(1,II),1,NBPVC,0,
     &                BCMVC(1,I) ,1,NBCVC,0,
     &                BEPVC(1,II),BPPVC(1,II),BCPVC(1,I),1,0,
     &                IBTYP,WORK,KFREE,LFREE)
          IF(DAMPFREQ.NE.D0) THEN
            CALL XRSTVC(BEMVC(1,II+1),1,NBEVC,0,
     &                  BPMVC(1,II+1),1,NBPVC,0,
     &                  BCMVC(1,I)   ,1,NBCVC,0,
     &                  BEPVC(1,II+1),BPPVC(1,II+1),BCPVC(1,I),1,0,
     &                  IBTYP,WORK,KFREE,LFREE)
          ENDIF
C
C         Scale T = S[2]X with negative frequency
C         ---------------------------------------
C
          IF(NZXOPE.GT.0) THEN
            CALL DCOPY(NDAMP*NZXOPEQ,BEBUF,1,BEMVC(1,II),1)
            CALL ACCS2XX(BEPVC(1,II),BEMVC(1,II),BEBUF,FREQ,NZXOPEQ)
          ENDIF
          IF(NZXOPP.GT.0) THEN
            CALL DCOPY(NDAMP*NZXOPPQ,BPBUF,1,BPMVC(1,II),1)
            CALL ACCS2XX(BPPVC(1,II),BPMVC(1,II),BPBUF,FREQ,NZXOPPQ)
          ENDIF
          IF(NZCONF.GT.0) THEN
            CALL DCOPY(NZCONFQ,BCBUF,1,BCMVC(1,I),1)
            CALL DSCAL(NZCONFQ,FAC,BCMVC(1,I),1)
            CALL DSCAL(NZCONFQ,FAC,BCPVC(1,I),1)
          ENDIF
        ELSE
C
C       IF FREQ = 0, then R(-) is automatically zero
C       --------------------------------------------
C
          IF(NZXOPE.GT.0) THEN
            CALL DZERO(BEPVC(1,I),NZXOPEQ)
            CALL DZERO(BEMVC(1,I),NZXOPEQ)
          ENDIF
          IF(NZXOPP.GT.0) THEN
            CALL DZERO(BPPVC(1,I),NZXOPPQ)
            CALL DZERO(BPMVC(1,I),NZXOPPQ)
          ENDIF
          IF(NZCONF.GT.0) THEN
            CALL DZERO(BCPVC(1,I),NZCONFQ)
            CALL DZERO(BCMVC(1,I),NZCONFQ)
          ENDIF
        ENDIF
        II = II + NDAMP
      ENDDO
C
      CALL QEXIT('ACCS2X')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ADDE2X */
      SUBROUTINE ADDE2X(ITYP,BEVEC,BPVEC,BCVEC,EVECR,MFREQ,
     &            IBTYP,BEBUF,BPBUF,BCBUF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     In the construction of the residual vector, add E[2}X.
C
C     Written by T.Saue Oct 17 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dcbibn.h"
      DIMENSION BEVEC(*),BPVEC(*),BCVEC(*),IBTYP(2,*),
     &          EVECR(*),BEBUF(*),BPBUF(*),BCBUF(*),WORK(*)
C
      CALL QENTER('ADDE2X')
C
C     e-e rotations
C
      IF(NZXOPE.GT.0) THEN
        CALL ADDSV1(ITYP,JBENDX,BEVEC,EVECR,MFREQ, IBTYP,BEBUF)
      ENDIF
C
C     e-p rotations
C
      IF(NZXOPP.GT.0) THEN
        CALL ADDSV1(ITYP,JBPNDX,BPVEC,EVECR,MFREQ, IBTYP,BPBUF)
      ENDIF
C
C     Configurational parameters
C
      IF(NZCONF.GT.0) THEN
        CALL ADDSV1(ITYP,JBCNDX,BCVEC,EVECR,MFREQ, IBTYP,BCBUF)
      ENDIF
C
      CALL QEXIT('ADDE2X')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ADDSV1 */
      SUBROUTINE ADDSV1(ITYP,NBTYP,BVEC,EVECR,MFREQ,IBTYP,BVBUF)
C***********************************************************************
C
C     Form electronic/positronic/configurational
C     part of solution vector in linear response
C
C     Written by T.Saue Oct 16 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0)
#include "dcbxrs.h"
#include "dcbibn.h"
      DIMENSION BVEC(*),EVECR(NEVEC,MFREQ),
     &               IBTYP(2,*),BVBUF(*)
C
C
      IF (NBTYP.EQ.JBCNDX) THEN
        NBDIM  = NZCONFQ
        LUBVEC = LUSCI
        IRECOFF = 1
      ELSEIF(NBTYP.EQ.JBENDX) THEN
        NBDIM  = NZXOPEQ
        LUBVEC = LUSOE
        IRECOFF = 0
      ELSEIF(NBTYP.EQ.JBPNDX) THEN
        NBDIM  = NZXOPPQ
        LUBVEC = LUSOP
        IRECOFF = 0
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'ADDSV1 ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('ADDSV1: Unknown NBTYP !')
      ENDIF
C
      IF(DAMPFREQ.EQ.D0) THEN
        DO IREC = 1,NZRED
          IH  = IBTYP(2,IREC)
          IF(IH.EQ.ITYP) THEN
 !           CALL READAC(LUBVEC,NBDIM,BVBUF,IREC+IRECOFF)
          CALL READ_DX(LUBVEC,IREC+IRECOFF,NBDIM,BVBUF)
            IBOFF = 1
            DO IFREQ = 1,MFREQ
              CALL DAXPY(NBDIM,EVECR(IREC,IFREQ),BVBUF,1,BVEC(IBOFF),1)
              IBOFF = IBOFF + NBDIM
            ENDDO
          ENDIF
        ENDDO
      ELSE
        DO IREC = 1,NZRED
          IH  = IBTYP(2,IREC)
          INDR = 2*IREC - 1
          INDI = 2*IREC
          IF(IH.EQ.ITYP) THEN
!            CALL READAC(LUBVEC,NBDIM,BVBUF,IREC+IRECOFF)
          CALL READ_DX(LUBVEC,IREC+IRECOFF,NBDIM,BVBUF)
            IBOFF = 1
            DO IFREQ = 1,MFREQ
              CALL DAXPY(NBDIM,EVECR(INDR,IFREQ),BVBUF,1,BVEC(IBOFF),1)
              IBOFF = IBOFF + NBDIM
              CALL DAXPY(NBDIM,EVECR(INDI,IFREQ),BVBUF,1,BVEC(IBOFF),1)
              IBOFF = IBOFF + NBDIM
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XVCNRM */
      SUBROUTINE XVCNRM(XNRM,XVEC,NDIM,MFREQ)
C***********************************************************************
C
C     Take norm of resultant vector X
C
C     Written by T.Saue Oct 17 1996
C
C     Last revision Oct 17 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C

      DIMENSION XVEC(NDIM,*),XNRM(*)
C
      DO I = 1,MFREQ
        XNRM(I) = DNRM2(NDIM,XVEC(1,I),1)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/*Deck xrsne1 */
      SUBROUTINE XRSNE1(BEVEC,BPVEC,BCVEC,MFREQ,EVALR,EVECR,
     &                  GPOE,GPOP,GPCI,IBTYP,IBCVC,IBEVC,IBPVC,
     &                  BEBUF,BPBUF,BCBUF,RCNV,RNORM,XNORM,IVECS,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     This routine will form new trial vectors for the
C     next microiteration.
C
C     Written by T.Saue Oct 18 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "dcbxrs.h"
#include "dcbibn.h"
      DIMENSION BEVEC(*),BPVEC(*),BCVEC(*),EVALR(*),
     &          EVECR(NEVEC,*),GPOE(*),GPOP(*),GPCI(*),
     &          IBTYP(2,NREDM),IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),
     &          BEBUF(*),BPBUF(*),BCBUF(*),IVECS(MFREQ,NVPAR),
     &          RNORM(NVPAR,MFREQ*NDAMP),XNORM(MFREQ*NDAMP,NVPAR),
     &          RCNV(*),WORK(*)
C
      IF (IPRXRS.GE.3) THEN
         CALL HEADER('Output from XRSNE1',-1)
         WRITE(LUPRI,'(A,2L10,I5)') 'STATIC, FKRMC, MFREQ',
     &      STATIC,FKRMC,MFREQ
      END IF
C
C     ***********************************
C     ***** Form solution vectors X *****
C     ***********************************
C
      IND = NZRED
C
C     Static case
C     ===========
C
      IF(STATIC) THEN
C
C       Orbital (e-e) part
C       ------------------
C
        IF(NZXOPE.GT.0) THEN
          NESIM = MFREQ
          IF((IND+NESIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          DO I = 1,NESIM
            IND = IND + 1
            IBTYP(1,IND)   = JBENDX
            IBTYP(2,IND)   = IPX
            IBEVC(NERED+I) = IND
          ENDDO
          CALL DZERO(BEVEC,NESIM*NZXOPEQ)
          IF(NERED.GT.0) THEN
            CALL XRSXV1(IPX,JBENDX,BEVEC,EVECR,MFREQ,
     &                  IBTYP,IBEVC,IVECS(1,IPEP),BEBUF)
            CALL XVCNRM(XNORM(1,IPEP),BEVEC,NZXOPEQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)')
     &          'Norm of orbital(e-e) X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,D8.2))') (XNORM(I,IPEP),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &          ' orbital (e-e) X(+) solution vectors'
              CALL PRBVEC(LUPRI,BEVEC,MFREQ,NZXOPE)
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPEP),MFREQ)
          ENDIF
        ENDIF
C
C       Orbital (e-p) part
C       ------------------
C
        IF(NZXOPP.GT.0) THEN
          NPSIM = MFREQ
          IF((IND+NPSIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          DO I = 1,NPSIM
            IND = IND + 1
            IBTYP(1,IND)   = JBPNDX
            IBTYP(2,IND)   = IPX
            IBPVC(NPRED+I) = IND
          ENDDO
          CALL DZERO(BPVEC,NPSIM*NZXOPPQ)
          IF(NPRED.GT.0) THEN
            CALL XRSXV1(IPX,JBPNDX,BPVEC,EVECR,MFREQ,
     &                  IBTYP,IBPVC,IVECS(1,IPPP),BPBUF)
            CALL XVCNRM(XNORM(1,IPPP),BPVEC,NZXOPPQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)')
     &          'Norm of orbital(e-p) X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPPP),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &          ' orbital (e-p) X(+) solution vectors'
              CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPPP),MFREQ)
          ENDIF
        ENDIF
C
C       Configurational part
C       --------------------
C
        IF(NZCONF.GT.0) THEN
          NCSIM = MFREQ
          IF((IND+NCSIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          DO I = 1,NCSIM
            IND = IND + 1
            IBTYP(1,IND)   = JBCNDX
            IBTYP(2,IND)   = IPX
            IBCVC(NCRED+I) = IND
          ENDDO
          CALL DZERO(BCVEC,NCSIM*NZCONFQ)
          IF(NCRED.GT.0) THEN
            CALL XRSXV1(IPX,JBCNDX,BCVEC,EVECR,MFREQ,
     &                  IBTYP,IBCVC,IVECS(1,IPCP),BCBUF)
            CALL XVCNRM(XNORM(1,IPCP),BCVEC,NZCONFQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)')
     &          'Norm of configurational X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPCP),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &          ' configurational X(+) solution vectors'
              CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPCP),NCSIM)
          ENDIF
        ENDIF
      ELSE ! not (STATIC)
C
C       In frequency-dependent case, the resultant X-vectors
C       split into X(+) and X(-).
C
        KBE2 = 1 + MFREQ*NZXOPEQ*NDAMP
        KBP2 = 1 + MFREQ*NZXOPPQ*NDAMP
        KBC2 = 1 + MFREQ*NZCONFQ
C
C       Orbital (e-e) part
C       ------------------
C
        IF(NZXOPE.GT.0) THEN
          NESIM = MFREQ*NSTAT*NDAMP
          IF((IND+NESIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          NBPREV = NERED
          DO IP = 1,NFREQ*NDAMP,NDAMP
            IM = IP + NFREQ*NDAMP
            DO IDAMP = 0,NDAMP-1
              IBTYP(1,IND+IP+IDAMP)     = JBENDX
              IBTYP(2,IND+IP+IDAMP)     = IPX
              IBEVC(NBPREV+IP+IDAMP)    = IND+IP+IDAMP
              IBTYP(1,IND+IM+IDAMP)     = JBENDX
              IBTYP(2,IND+IM+IDAMP)     = IMX
              IBEVC(NBPREV+IM+IDAMP)    = IND+IM+IDAMP
            ENDDO
          ENDDO
          IND = IND + NESIM         
          CALL DZERO(BEVEC,NESIM*NZXOPEQ)
          IF(NERED.GT.0) THEN
            CALL XRSXV1(IPX,JBENDX,BEVEC      ,EVECR,MFREQ,
     &                  IBTYP,IBEVC,IVECS(1,IPEP),BEBUF)
            CALL XVCNRM(XNORM(1,IPEP),BEVEC      ,NZXOPEQ,MFREQ*NDAMP)
            CALL XRSXV1(IMX,JBENDX,BEVEC(KBE2),EVECR,MFREQ,
     &                  IBTYP,IBEVC,IVECS(1,IPEM),BEBUF)
            CALL XVCNRM(XNORM(1,IPEM),BEVEC(KBE2),NZXOPEQ,MFREQ*NDAMP)
C           Print section
            IF(IPRXRS.GE.3) THEN
              IF(DAMPFREQ.NE.D0) THEN
                WRITE(LUPRI,'(A)')
     &            "Norm of orbital(e-e) X(+)-solution vectors:"
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &            (XNORM(I,IPEP),I=1,MFREQ*NDAMP,NDAMP)
                WRITE(LUPRI,'(A)')
     &            "Norm of orbital(e-e) X'(+)-solution vectors:"
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &            (XNORM(I,IPEP),I=2,MFREQ*NDAMP,NDAMP)
                WRITE(LUPRI,'(A)')
     &            "Norm of orbital(e-e) X(-)-solution vectors:"
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &            (XNORM(I,IPEM),I=1,MFREQ*NDAMP,NDAMP)
                WRITE(LUPRI,'(A)')
     &            "Norm of orbital(e-e) X'(-)-solution vectors:"
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &            (XNORM(I,IPEM),I=2,MFREQ*NDAMP,NDAMP)
              ELSE
                WRITE(LUPRI,'(A)')
     &            'Norm of orbital(e-e) X(+)-solution vectors:'
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &            (XNORM(I,IPEP),I=1,MFREQ*NDAMP)
                WRITE(LUPRI,'(A)')
     &            'Norm of orbital(e-e) X(-)-solution vectors:'
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &            (XNORM(I,IPEM),I=1,MFREQ*NDAMP)
              ENDIF
            ENDIF
            IF(IPRXRS.GE.5) THEN
              IF(DAMPFREQ.NE.D0) THEN
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           " orbital (e-e) X(+) solution vectors"
               CALL PRBVCC(LUPRI,BEVEC,1,MFREQ*NDAMP,NDAMP,NZXOPE)
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           " orbital (e-e) X'(+) solution vectors"
               CALL PRBVCC(LUPRI,BEVEC,2,MFREQ*NDAMP,NDAMP,NZXOPE)
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           " orbital (e-e) X(-) solution vectors"
               CALL PRBVCC(LUPRI,BEVEC(KBE2),1,MFREQ*NDAMP,NDAMP,NZXOPE)
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           " orbital (e-e) X'(-) solution vectors"
               CALL PRBVCC(LUPRI,BEVEC(KBE2),2,MFREQ*NDAMP,NDAMP,NZXOPE)
              ELSE
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           ' orbital (e-e) X(+) solution vectors'
               CALL PRBVEC(LUPRI,BEVEC      ,MFREQ,NZXOPE)
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           ' orbital (e-e) X(-) solution vectors'
               CALL PRBVEC(LUPRI,BEVEC(KBE2),MFREQ,NZXOPE)
              ENDIF
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPEP),MFREQ*NDAMP)
            CALL DZERO(XNORM(1,IPEM),MFREQ*NDAMP)
          ENDIF
        ENDIF
C
C       Orbital (e-p) part
C       ------------------
C
        IF(NZXOPP.GT.0) THEN
          NPSIM = MFREQ*NSTAT*NDAMP
          IF((IND+NPSIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          NBPREV = NPRED
          DO IP = 1,NFREQ*NDAMP,NDAMP
            IM = IP + NFREQ*NDAMP
            DO IDAMP = 0,NDAMP-1
              IBTYP(1,IND+IP+IDAMP)     = JBPNDX
              IBTYP(2,IND+IP+IDAMP)     = IPX
              IBPVC(NBPREV+IP+IDAMP)    = IND+IP+IDAMP
              IBTYP(1,IND+IM+IDAMP)     = JBPNDX
              IBTYP(2,IND+IM+IDAMP)     = IMX
              IBPVC(NBPREV+IM+IDAMP)    = IND+IM+IDAMP
            ENDDO
          ENDDO
          IND = IND + NPSIM         
          CALL DZERO(BPVEC,NPSIM*NZXOPPQ)
          IF(NPRED.GT.0) THEN
            CALL XRSXV1(IPX,JBPNDX,BPVEC      ,EVECR,MFREQ,
     &                  IBTYP,IBPVC,IVECS(1,IPPP),BPBUF)
            CALL XVCNRM(XNORM(1,IPPP),BPVEC      ,NZXOPPQ,MFREQ*NDAMP)
            CALL XRSXV1(IMX,JBPNDX,BPVEC(KBP2),EVECR,MFREQ,
     &                  IBTYP,IBPVC,IVECS(1,IPPM),BPBUF)
            CALL XVCNRM(XNORM(1,IPPM),BPVEC(KBP2),NZXOPPQ,MFREQ*NDAMP)
C           Print section
            IF(IPRXRS.GE.3) THEN
              IF(DAMPFREQ.NE.D0) THEN
                WRITE(LUPRI,'(A)')
     &            "Norm of orbital(e-p) X(+)-solution vectors:"
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &            (XNORM(I,IPPP),I=1,MFREQ*NDAMP,NDAMP)
                WRITE(LUPRI,'(A)')
     &            "Norm of orbital(e-p) X'(+)-solution vectors:"
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &            (XNORM(I,IPPP),I=2,MFREQ*NDAMP,NDAMP)
                WRITE(LUPRI,'(A)')
     &            "Norm of orbital(e-p) X(-)-solution vectors:"
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &            (XNORM(I,IPPM),I=1,MFREQ*NDAMP,NDAMP)
                WRITE(LUPRI,'(A)')
     &            "Norm of orbital(e-p) X'(-)-solution vectors:"
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &            (XNORM(I,IPPM),I=2,MFREQ*NDAMP,NDAMP)
              ELSE
                WRITE(LUPRI,'(A)')
     &            'Norm of orbital(e-p) X(+)-solution vectors:'
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &           (XNORM(I,IPPP),I=1,MFREQ)
                WRITE(LUPRI,'(A)')
     &            'Norm of orbital(e-p) X(-)-solution vectors:'
                WRITE(LUPRI,'(1P,8(3X,E8.2))')
     &           (XNORM(I,IPPM),I=1,MFREQ)
              ENDIF
            ENDIF
            IF(IPRXRS.GE.5) THEN
              IF(DAMPFREQ.NE.D0) THEN
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           " orbital (e-p) X(+) solution vectors"
               CALL PRBVCC(LUPRI,BPVEC,1,MFREQ*NDAMP,NDAMP,NZXOPP)
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           " orbital (e-p) X'(+) solution vectors"
               CALL PRBVCC(LUPRI,BPVEC,2,MFREQ*NDAMP,NDAMP,NZXOPP)
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           " orbital (e-p) X(-) solution vectors"
               CALL PRBVCC(LUPRI,BPVEC(KBP2),1,MFREQ*NDAMP,NDAMP,NZXOPP)
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           " orbital (e-p) X'(-) solution vectors"
               CALL PRBVCC(LUPRI,BPVEC(KBP2),2,MFREQ*NDAMP,NDAMP,NZXOPP)
              ELSE
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           ' orbital (e-p) X(+) solution vectors'
               CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
               WRITE(LUPRI,'(/I5,A)') MFREQ,
     &           ' orbital (e-p) X(-) solution vectors'
               CALL PRBVEC(LUPRI,BPVEC(KBP2),MFREQ,NZXOPP)
              ENDIF
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPPP),MFREQ*NDAMP)
            CALL DZERO(XNORM(1,IPPM),MFREQ*NDAMP)
          ENDIF
        ENDIF
C
C       Configurational part
C       --------------------
C
        IF(NZCONF.GT.0) THEN
          NCSIM = MFREQ + MFREQ
          IF((IND+NCSIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          IOFF = NCRED
          DO I = 1,MFREQ
            IND = IND + 1
            IBTYP(1,IND)  = JBCNDX
            IBTYP(2,IND)  = IPX
            IBCVC(IOFF+I) = IND
          ENDDO
          IOFF = NCRED + MFREQ
          DO I = 1,MFREQ
            IND = IND + 1
            IBTYP(1,IND) = JBCNDX
            IBTYP(2,IND) = IMX
            IBCVC(IOFF)  = IND
          ENDDO
          CALL DZERO(BCVEC,NCSIM*NZCONFQ)
          IF(NCRED.GT.0) THEN
            CALL XRSXV1(IPX,JBCNDX,BCVEC      ,EVECR,MFREQ,
     &                  IBTYP,IBCVC,IVECS(1,IPCP),BCBUF)
            CALL XRSXV1(IMX,JBCNDX,BCVEC      ,EVECR,MFREQ,
     &                  IBTYP,IBCVC,IVECS(1,IPCM),BCBUF)
            CALL XVCNRM(XNORM(1,IPCP),BCVEC      ,NZCONFQ,MFREQ)
            CALL XVCNRM(XNORM(1,IPCM),BCVEC(KBC2),NZCONFQ,MFREQ)
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)')
     &          'Norm of configurational X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPCP),I=1,MFREQ)
              WRITE(LUPRI,'(A)')
     &          'Norm of configurational X(-)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPCP),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &          ' configurational X(+) solution vectors'
              CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &          ' configurational X(-) solution vectors'
              CALL PRBVEC(LUPRI,BCVEC(KBC2),MFREQ,NZCONF)
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPCP),MFREQ)
            CALL DZERO(XNORM(1,IPCM),MFREQ)
          ENDIF
        ENDIF
      ENDIF ! from IF(STATIC) THEN ... ELSE ...
C
      NTSIM = NESIM + NPSIM + NCSIM
C
C     Print section
C
      IF(IPRXRS.GT.0) THEN
        WRITE(LUPRI,'(/A)') '* XRSNEX: Norm of resultant vectors:'
        IF(DAMPFREQ.EQ.D0) THEN
          WRITE(LUPRI,'(A15,2X,A8,9(2X,A2,A6))')
     &      '* Frequency    ',
     &       'Total   ',(BVTYP(I),' part ',I=1,NVPAR)
          DO IFREQ = 1,MFREQ
CPP check XTNORM
            XTNORM = DNRM2(NVPAR,XNORM(IFREQ,1),MFREQ)
            WRITE(LUPRI,'(F15.8,1P,10(2X,E8.2))')
     &          EVALR(IFREQ),XTNORM,(XNORM(IFREQ,I),I=1,NVPAR)
          ENDDO
        ELSE
          DO IFREQ = 1,MFREQ
            INDR = 2 * IFREQ - 1
            INDI = 2 * IFREQ
            XTNORMR = DNRM2(NVPAR,XNORM(INDR,1),MFREQ*NDAMP)
            XTNORMI = DNRM2(NVPAR,XNORM(INDR,1),MFREQ*NDAMP)
            XTNORM = DSQRT(XTNORMR**2 + XTNORMI**2)
CPP Need to be changed with correct complex label ++,+-,--,-+!
            WRITE(LUPRI,'(A15,2X,A8,9(2X,A2,A6))')
     &       '* Frequency    ',
     &       'Total   ',(BVTYP(I),' part ',I=1,NVPAR),
     &                  (BVTYP(I),' part ',I=1,NVPAR)
            WRITE(LUPRI,'(F15.8,1P,16(2X,E8.2))')
     &          EVALR(IFREQ),XTNORM,
     &         (XNORM(INDR,I),I=1,NVPAR),
     &         (XNORM(INDI,I),I=1,NVPAR)
            WRITE(LUPRI,'(A,F15.8)') '  Damping :',DAMPFREQ
          ENDDO
        ENDIF
      ENDIF
      CALL FLSHFO(LUPRI)
C
C     ****************************
C     ***** Form residuals R *****
C     ****************************
C
C     For KRMCSCF form -EVALR(1)*X for restricted step optimization
C     (the level shift EVALR(1) is calculated in KRMCSOL)
C
C     In frequency-dependent case, form -W*S[2]X
C     (in static case W = 0 )
C     ==========================================
C
      IF ( FKRMC ) THEN
         CALL DSCAL(NZXOPEQ,-EVALR(1),BEVEC,1)
         CALL DSCAL(NZXOPPQ,-EVALR(1),BPVEC,1)
         CALL DSCAL(NZCONFQ,-EVALR(1),BCVEC,1)
      ELSE IF(STATIC) THEN
         IF(NESIM.GT.0) CALL DZERO(BEVEC,NZXOPEQ*NESIM)
         IF(NPSIM.GT.0) CALL DZERO(BPVEC,NZXOPPQ*NPSIM)
         IF(NCSIM.GT.0) CALL DZERO(BCVEC,NZCONFQ*NCSIM)
      ELSE
        CALL ACCS2X(BEVEC,BEVEC(KBE2),BPVEC,BPVEC(KBP2),
     &              BCVEC,BCVEC(KBC2),IBTYP,IBCVC,IBEVC,IBPVC,
     &              EVECR,EVALR,MFREQ,BEBUF,BPBUF,BCBUF,
     &              WORK,KFREE,LFREE)
      ENDIF
C
C     Print section
      IF(IPRXRS.GE.5) THEN
        IF(NZXOPE.GT.0) THEN
          IF(DAMPFREQ.NE.D0) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       " orbital (e-e) -WS[2]X(-)+YS[2]X'(-) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC,1,MFREQ*NDAMP,NDAMP,NZXOPE)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       " orbital (e-e) -WS[2]X'(-)-YS[2]X(-) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC,2,MFREQ*NDAMP,NDAMP,NZXOPE)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       " orbital (e-e) -WS[2]X(+)+YS[2]X'(+) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC(KBE2),1,MFREQ*NDAMP,NDAMP,NZXOPE)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       " orbital (e-e) -WS[2]X'(+)-YS[2]X(+) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC(KBE2),2,MFREQ*NDAMP,NDAMP,NZXOPE)
          ELSE
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       ' orbital (e-e) -WS[2]X(-) residual vectors'
            CALL PRBVEC(LUPRI,BEVEC,MFREQ,NZXOPE)
            IF(.NOT.STATIC) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &         ' orbital (e-e) -WS[2]X(+) residual vectors'
              CALL PRBVEC(LUPRI,BEVEC(KBE2),MFREQ,NZXOPE)
            ENDIF
          ENDIF
        ENDIF
        IF(NZXOPP.GT.0) THEN
          IF(DAMPFREQ.NE.D0) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       " orbital (e-p) -WS[2]X(-)+YS[2]X'(-) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC,1,MFREQ*NDAMP,NDAMP,NZXOPP)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       " orbital (e-p) -WS[2]X'(-)-YS[2]X(-) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC,2,MFREQ*NDAMP,NDAMP,NZXOPP)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       " orbital (e-p) -WS[2]X(+)+YS[2]X'(+) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC(KBP2),1,MFREQ*NDAMP,NDAMP,NZXOPP)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       " orbital (e-p) -WS[2]X'(+)-YS[2]X(+) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC(KBP2),2,MFREQ*NDAMP,NDAMP,NZXOPP)
          ELSE
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       ' orbital (e-p) -WS[2]X(-) residual vectors'
            CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
            IF(.NOT.STATIC) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &         ' orbital (e-p) -WS[2]X(+) residual vectors'
              CALL PRBVEC(LUPRI,BPVEC(KBP2),MFREQ,NZXOPP)
            ENDIF
          ENDIF
        ENDIF
        IF(NZCONF.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ,
     &     ' configurational -WS[2]X(-) residual vectors'
          CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       ' configurational -WS[2]X(+) residual vectors'
            CALL PRBVEC(LUPRI,BCVEC(KBC2),MFREQ,NZCONF)
          ENDIF
        ENDIF
      ENDIF
C
C     Add E[2]X
C     =========
C
      CALL ADDE2X(IPX,BEVEC,BPVEC,BCVEC,EVECR,MFREQ,
     &            IBTYP,BEBUF,BPBUF,BCBUF,WORK,KFREE,LFREE)
      IF(.NOT.STATIC) THEN
        CALL ADDE2X(IMX,BEVEC(KBE2),BPVEC(KBP2),BCVEC(KBC2),EVECR,MFREQ,
     &              IBTYP,BEBUF,BPBUF,BCBUF,WORK,KFREE,LFREE)

      ENDIF
C
C     Print section
C
      IF(IPRXRS.GE.5) THEN
        IF(NZXOPE.GT.0) THEN
          IF(DAMPFREQ.NE.D0) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-e) E[2]X(+)-WS[2]X(-)+YS[2]X'(-) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC,1,MFREQ*NDAMP,NDAMP,NZXOPE)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-e) E[2]X'(+)-WS[2]X'(-)-YS[2]X(-) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC,2,MFREQ*NDAMP,NDAMP,NZXOPE)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-e) E[2]X(-)-WS[2]X(+)+YS[2]X'(+) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC(KBE2),1,MFREQ*NDAMP,NDAMP,NZXOPE)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-e) E[2]X'(-)-WS[2]X'(+)-YS[2]X(+) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC(KBE2),2,MFREQ*NDAMP,NDAMP,NZXOPE)
          ELSE
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       ' orbital (e-e) (E[2]X(+)-WS[2]X(-)) residual vectors'
            CALL PRBVEC(LUPRI,BEVEC,MFREQ,NZXOPE)
            IF(.NOT.STATIC) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &         ' orbital (e-e) (E[2]X(-)-WS[2]X(+)) residual vectors'
              CALL PRBVEC(LUPRI,BEVEC(KBE2),MFREQ,NZXOPE)
            ENDIF
          ENDIF
        ENDIF
        IF(NZXOPP.GT.0) THEN
          IF(DAMPFREQ.NE.D0) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-p) E[2]X(+)-WS[2]X(-)+YS[2]X'(-) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC,1,MFREQ*NDAMP,NDAMP,NZXOPP)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-p) E[2]X'(+)-WS[2]X'(-)-YS[2]X(-) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC,2,MFREQ*NDAMP,NDAMP,NZXOPP)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-p) E[2]X(-)-WS[2]X(+)+YS[2]X'(+) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC(KBP2),1,MFREQ*NDAMP,NDAMP,NZXOPP)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-p) E[2]X'(-)-WS[2]X'(+)-YS[2]X(+) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC(KBP2),2,MFREQ*NDAMP,NDAMP,NZXOPP)
          ELSE
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       ' orbital (e-p) (E[2]X(+)-WS[2]X(-)) residual vectors'
            CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
            IF(.NOT.STATIC) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &         ' orbital (e-p) (E[2]X(-)-WS[2]X(+)) residual vectors'
              CALL PRBVEC(LUPRI,BPVEC(KBP2),MFREQ,NZXOPP)
            ENDIF
          ENDIF
        ENDIF
        IF(NZCONF.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ,
     &     ' configurational (E[2]X(+)-WS[2]X(-)) residual vectors'
          CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       'configurational (E[2]X(-)-WS[2]X(+)) residual vectors'
            CALL PRBVEC(LUPRI,BCVEC(KBC2),MFREQ,NZCONF)
          ENDIF
        ENDIF
      ENDIF ! from IF(IPRXRS.GE.5) THEN
C
C     If (LINEQ) add gradient
C     =======================
C
      IF(LINEQ) THEN
        IF(NZXOPE.GT.0) THEN
          CALL GRDADD(BEVEC,GPOE,NZXOPEQ,MFREQ*NDAMP)
        ENDIF
        IF(NZXOPP.GT.0) THEN
          CALL GRDADD(BPVEC,GPOP,NZXOPPQ,MFREQ*NDAMP)
        ENDIF
        IF(NZCONF.GT.0) THEN
          CALL GRDADD(BCVEC,GPCI,NZCONFQ,MFREQ*NDAMP)
        ENDIF
C
C       Print section
C
CPP...need to take care of complex case printing
      IF(IPRXRS.GE.5) THEN
        IF(NZXOPE.GT.0) THEN
          IF(DAMPFREQ.NE.D0) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     & " orbital (e-e) E[2]X(+)-WS[2]X(-)+YS[2]X'(-)+G residual vectors"
            CALL PRBVCC(LUPRI,BEVEC,1,MFREQ*NDAMP,NDAMP,NZXOPE)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-e) E[2]X'(+)-WS[2]X'(-)-YS[2]X(-) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC,2,MFREQ*NDAMP,NDAMP,NZXOPE)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-e) E[2]X(-)-WS[2]X(+)+YS[2]X'(+) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC(KBE2),1,MFREQ*NDAMP,NDAMP,NZXOPE)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-e) E[2]X'(-)-WS[2]X'(+)-YS[2]X(+) residual vectors"
            CALL PRBVCC(LUPRI,BEVEC(KBE2),2,MFREQ*NDAMP,NDAMP,NZXOPE)
          ELSE
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       ' orbital (e-e) (E[2]X(+)-WS[2]X(-)+G) residual vectors'
            CALL PRBVEC(LUPRI,BEVEC,MFREQ,NZXOPE)
            IF(.NOT.STATIC) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &         ' orbital (e-e) (E[2]X(-)-WS[2]X(+)) residual vectors'
              CALL PRBVEC(LUPRI,BEVEC(KBE2),MFREQ,NZXOPE)
            ENDIF
          ENDIF
        ENDIF
        IF(NZXOPP.GT.0) THEN
          IF(DAMPFREQ.NE.D0) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     & " orbital (e-p) E[2]X(+)-WS[2]X(-)+YS[2]X'(-)+G residual vectors"
            CALL PRBVCC(LUPRI,BPVEC,1,MFREQ*NDAMP,NDAMP,NZXOPP)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-p) E[2]X'(+)-WS[2]X'(-)-YS[2]X(-) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC,2,MFREQ*NDAMP,NDAMP,NZXOPP)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-p) E[2]X(-)-WS[2]X(+)+YS[2]X'(+) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC(KBP2),1,MFREQ*NDAMP,NDAMP,NZXOPP)
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &  " orbital (e-p) E[2]X'(-)-WS[2]X'(+)-YS[2]X(+) residual vectors"
            CALL PRBVCC(LUPRI,BPVEC(KBP2),2,MFREQ*NDAMP,NDAMP,NZXOPP)
          ELSE
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       ' orbital (e-p) (E[2]X(+)-WS[2]X(-)+G) residual vectors'
            CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
            IF(.NOT.STATIC) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ,
     &         ' orbital (e-p) (E[2]X(-)-WS[2]X(+)) residual vectors'
              CALL PRBVEC(LUPRI,BPVEC(KBP2),MFREQ,NZXOPP)
            ENDIF
          ENDIF
        ENDIF
        IF(NZCONF.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ,
     &     ' configurational (E[2]X(+)-WS[2]X(-)+G) residual vectors'
          CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ,
     &       ' configurational (E[2]X(-)-WS[2]X(+)) residual vectors'
            CALL PRBVEC(LUPRI,BCVEC(KBC2),MFREQ,NZCONF)
            ENDIF
          ENDIF
        ENDIF
      ENDIF ! from IF(LINEQ) THEN
C
      CALL FLSHFO(LUPRI)
C
C     **********************************
C     ***** Form new trial vectors *****
C     **********************************
C
      CALL XRSBVC(BEVEC,BEBUF,IBEVC,
     &            BPVEC,BPBUF,IBPVC,
     &            BCVEC,BCBUF,IBCVC,
     &            MFREQ,IBTYP,IVECS,EVALR,EVECR,RCNV,RNORM,XNORM,
     &            WORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GRDADD */
      SUBROUTINE GRDADD(BVEC,GP,NBDIM,MFREQ)
C***********************************************************************
C
C     Add property gradient to residual
C
C     Written by T.Saue Oct 17 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
#include "dcbxrs.h"
C
      DIMENSION BVEC(NBDIM,MFREQ),GP(NBDIM)
C

      DO I = 1,MFREQ,NDAMP
        CALL DAXPY(NBDIM,D1,GP,1,BVEC(1,I),1)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XRSBVC */
      SUBROUTINE XRSBVC(BEVEC,BEBUF,IBEVC,
     &                  BPVEC,BPBUF,IBPVC,
     &                  BCVEC,BCBUF,IBCVC,
     &                  MFREQ,IBTYP,IVECS,EVALR,EVECR,RCNV,RNORM,XNORM,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     Construct new trial vectors for next microiteration
C
C     The selection of trial vectors proceed in three steps:
C       1) If |R| .LE. (THRS*|X|)
C             where |R| is norm of residual
C             and   |X| is norm of current solution vector
C          then no new trial vectors are generated from this root.
C       2) If RESFAC*|R*| .LE. |RMAX|
C             where |R*| is the norm of a partial residual
C             (E-,P-,C-type)
C             and |RMAX| is the largest partial residual
C          then no new trial vectors of this type is generated from
C             this root.
C          If RESFAC = 1, then only one trial vector is generated
C             from each root.
C          If RESFAC very big, then all possible trial vectors are
C             generated from the root.
C
C     Input:
C       BEVEC - residual : orbital (e-e) part
C       BPVEC - residual : orbital (e-p) part
C       BCVEC - residual : configurational part
C
C       XECNV - norm of resultant vector: orbital (e-e) part
C       XPCNV - norm of resultant vector: orbital (e-p) part
C       XCCNV - norm of resultant vector: configurational part
C
C     Output:
C
C       BEVEC - NESIM trial vectors: orbital (e-e) type
C       BPVEC - NPSIM trial vectors: orbital (e-p) type
C       BCVEC - NCSIM trial vectors : configurational part
C
C       RNORM - residual norm
C
C     Written by T.Saue, Oct 18 1996
C     Last modification Nov 26 1997
C
C***********************************************************************

      use dirac_cfg

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,D5=5.0D0)
C
#include "dcbxrs.h"
#include "dcbibn.h"
      LOGICAL   USEXVC
      DIMENSION BEVEC(NZXOPEQ,*),BEBUF(*),IBEVC(NREDM),
     &          BPVEC(NZXOPPQ,*),BPBUF(*),IBPVC(NREDM),
     &          BCVEC(NZCONFQ,*),BCBUF(*),IBCVC(NREDM),
     &          IBTYP(2,NREDM),IVECS(MFREQ*NDAMP,NVPAR),EVALR(*),
     &          EVECR(*),RCNV(*),RNORM(NVPAR,MFREQ*NDAMP),
     &          XNORM(MFREQ*NDAMP,NVPAR),WORK(*)
C
      CALL QENTER('XRSBVC')
      CNVMAX = D0
      IF(DAMPFREQ.NE.D0) THEN
        WRITE(LUPRI,'(/A)') '* Norm of residual vectors:'
        WRITE(LUPRI,'(A12,2X,A9,1X,A8,18(2X,A2,A6))')
     &        '  Frequency ','Threshold ','Total   ',
     &       (BVTYP(I),' part ',I=1,NVPAR),
     &       (BVTYP(I),' part ',I=1,NVPAR)
      ELSE
        WRITE(LUPRI,'(/A)') '* Norm of residual vectors:'
        WRITE(LUPRI,'(A15,2X,A9,1X,A8,9(2X,A2,A6))')
     &        '  Frequency    ','Threshold ','Total   ',
     &        (BVTYP(I),' part ',I=1,NVPAR)
      ENDIF
      DO IFREQ = 1,MFREQ
         INDR = 2*IFREQ-1
         INDI = 2*IFREQ
C
C        ************************************
C        ***** Get norms of residual vectors
C        ************************************
C
C        Orbital (e-e) part
C        ------------------
C
         IF(NZXOPE.GT.0) THEN
           IF(DAMPFREQ.NE.D0) THEN
             IVECS(INDR,IPEP) = 1
             RNORM(IPEP,INDR) = DNRM2(NZXOPEQ,BEVEC(1,INDR),1)
             IVECS(INDI,IPEP) = 1
             RNORM(IPEP,INDI) = DNRM2(NZXOPEQ,BEVEC(1,INDI),1)
             IVECS(INDR,IPEM) = 1
             IOFF = INDR + MFREQ*NDAMP
             RNORM(IPEM,INDR) = DNRM2(NZXOPEQ,BEVEC(1,IOFF),1)
             IVECS(INDI,IPEM) = 1
             IOFF = INDI + MFREQ*NDAMP
             RNORM(IPEM,INDI) = DNRM2(NZXOPEQ,BEVEC(1,IOFF),1)
           ELSE
             IVECS(IFREQ,IPEP) = 1
             RNORM(IPEP,IFREQ) = DNRM2(NZXOPEQ,BEVEC(1,IFREQ),1)
             IF(.NOT.STATIC) THEN
                IVECS(IFREQ,IPEM) = 1
                IOFF = IFREQ + MFREQ
                RNORM(IPEM,IFREQ) = DNRM2(NZXOPEQ,BEVEC(1,IOFF),1)
             ENDIF
           ENDIF
         ENDIF
C
C        Orbital (e-p) part
C        ------------------
C
         IF(NZXOPP.GT.0) THEN
           IF(DAMPFREQ.NE.D0) THEN
             IVECS(INDR,IPPP) = 1
             RNORM(IPPP,INDR) = DNRM2(NZXOPPQ,BPVEC(1,INDR),1)
             IVECS(INDI,IPPP) = 1
             RNORM(IPPP,INDI) = DNRM2(NZXOPPQ,BPVEC(1,INDI),1)
             IVECS(INDR,IPPM) = 1
             IOFF = INDR + MFREQ*NDAMP
             RNORM(IPPM,INDR) = DNRM2(NZXOPPQ,BPVEC(1,IOFF),1)
             IVECS(INDI,IPPM) = 1
             IOFF = INDI + MFREQ*NDAMP
             RNORM(IPPM,INDI) = DNRM2(NZXOPPQ,BPVEC(1,IOFF),1)
           ELSE
             IVECS(IFREQ,IPPP) = 1
             RNORM(IPPP,IFREQ) = DNRM2(NZXOPPQ,BPVEC(1,IFREQ),1)
             IF(.NOT.STATIC) THEN
                IOFF = IFREQ + MFREQ
                IVECS(IFREQ,IPPM) = 1
                RNORM(IPPM,IFREQ) = DNRM2(NZXOPPQ,BPVEC(1,IOFF),1)
             ENDIF
           ENDIF
         ENDIF
C
C        Configurational part
C        --------------------
C
         IF(NZCONF.GT.0) THEN
            IVECS(IFREQ,IPCP) = 1
            RNORM(IPCP,IFREQ) = DNRM2(NZCONFQ,BCVEC(1,IFREQ),1)
            IF(.NOT.STATIC) THEN
               IOFF = IFREQ + MFREQ
               IVECS(IFREQ,IPCM) = 1
               RNORM(IPCM,IFREQ) = DNRM2(NZCONFQ,BCVEC(1,IOFF),1)
            ENDIF
         ENDIF
C
C        *****************************
C        ***** Check convergence *****
C        *****************************
C
         IF(DAMPFREQ.NE.D0) THEN
           RNRMAX = D0
           DO I = 1,NVPAR
              RNRMAX = MAX(RNRMAX,RNORM(I,INDR))
              RNRMAX = MAX(RNRMAX,RNORM(I,INDI))
           ENDDO
           RTNORMR = DNRM2(NVPAR,RNORM(1,INDR),1)
           RTNORMI = DNRM2(NVPAR,RNORM(1,INDI),1)
           RTNORM = DSQRT(RTNORMR**2 + RTNORMI**2)
           XTNORMR = DNRM2(NVPAR,XNORM(INDR,1),MFREQ*NDAMP)
           XTNORMI = DNRM2(NVPAR,XNORM(INDI,1),MFREQ*NDAMP)
           XTNORM = DSQRT(XTNORMR**2 + XTNORMI**2)
         ELSE
           RNRMAX = D0
           DO I = 1,NVPAR
              RNRMAX = MAX(RNRMAX,RNORM(I,IFREQ))
           ENDDO
           RTNORM = DNRM2(NVPAR,RNORM(1,IFREQ),1)
           XTNORM = DNRM2(NVPAR,XNORM(IFREQ,1),MFREQ)
         ENDIF
         IF ( .NOT. FKRMC ) THEN
            IF(XTNORM.EQ.D0) THEN
               if (dirac_cfg_openrsp) then
                  do i = 1, nvpar
                     ivecs(ifreq, i) = 0
                  end do
               end if
               WRITE(LUPRI,'(A)')
     &         '*** WARNING! *** Zero norm of solution vector !'
            ELSE
               RCNV(IFREQ) = RTNORM/GPTNRM
               CNVMAX = MAX(CNVMAX,RCNV(IFREQ))
            ENDIF
C
C           We are calculating properties.
C           Use relative norm as conv. measure
C
            IF(RTNORM.LE.(GPTNRM*THCXRS)) THEN
C              ... this root converged
              IF(DAMPFREQ.NE.D0) THEN
                DO I = 1,NVPAR
                   IVECS(INDR,I) = 0
                   IVECS(INDI,I) = 0
                ENDDO
              ELSE
                DO I = 1,NVPAR
                   IVECS(IFREQ,I) = 0
                ENDDO
              ENDIF
            ELSE
C              hjj+patrick:
C              Check each type separately whether residual norm is less than
C              a fraction of the maximun residual norm, or if less than
C              a fraction (0.2) of convergence threshold.
               IF (LINEQ) THEN
                 IF(DAMPFREQ.NE.D0) THEN
                   DO I = 1,NVPAR
                     IF(RESFAC*RNORM(I,INDR).LT.RNRMAX .OR.
     &                  D5*RNORM(I,INDR).LE.GPTNRM*THCXRS)
     &                  IVECS(INDR,I) = 0
                     IF(RESFAC*RNORM(I,INDI).LT.RNRMAX .OR.
     &                  D5*RNORM(I,INDI).LE.GPTNRM*THCXRS)
     &                  IVECS(INDI,I) = 0
                   ENDDO
                 ELSE
                   DO I = 1,NVPAR
                     IF(RESFAC*RNORM(I,IFREQ).LT.RNRMAX .OR.
     &                  D5*RNORM(I,IFREQ).LE.GPTNRM*THCXRS)
     &                  IVECS(IFREQ,I) = 0
                   ENDDO
                 ENDIF
               ELSE
C                hjaaj jan 2001: for PP excitation energies we must include
C                both the X(+) and the X(-) vectors to keep the paired structure ...
                 DO I = 1,NVPAR,2
                   RNTEST = RNORM(I,IFREQ) + RNORM(I+1,IFREQ)
                   IF(RESFAC*RNTEST.LT.RNRMAX.OR.
     &                D5*RNTEST.LE.GPTNRM*THCXRS) THEN
                     IVECS(IFREQ,I) = 0
                     IVECS(IFREQ,I+1) = 0
                   END IF
                 ENDDO
               END IF
            END IF
         ELSE
C
C           **** KRMCSCF code **** (FKRMC is .true.)
C
            RCNV(IFREQ) = RTNORM
            CNVMAX = MAX(CNVMAX,RCNV(IFREQ))
C
C           Update conv. threshold
C
            THCXRS = RCNVTHR()
C
C           We are doing a MCSCF opt.
C           Use abs. norm as conv. measure.
C
            IF ( RTNORM .LE. THCXRS ) THEN
               DO I = 1,NVPAR
                  IVECS(IFREQ,I) = 0
               ENDDO
            ELSE
C              hjj+patrick:
C              Check each type separately whether residual norm is
C              less than a fraction of the maximun residual norm, or
C              if less than a fraction (0.2) of convergence threshold.
               DO I = 1,NVPAR
                  IF(RESFAC*RNORM(I,IFREQ).LT.RNRMAX.OR.
     &                 D5*RNORM(I,IFREQ).LE.THCXRS) IVECS(IFREQ,I) = 0
               ENDDO
            ENDIF
         ENDIF
         IF(DAMPFREQ.NE.D0) THEN
           WRITE(LUPRI,'(F12.8,1P,16(2X,E8.2))')
     &        EVALR(IFREQ),THCXRS,RTNORM,
     &        (RNORM(I,INDR),I=1,NVPAR),
     &        (RNORM(I,INDI),I=1,NVPAR)
         ELSE
           WRITE(LUPRI,'(F15.8,1P,11(2X,E8.2))')
     &        EVALR(IFREQ),THCXRS,RTNORM,(RNORM(I,IFREQ),I=1,NVPAR)
         ENDIF
      ENDDO
C
C     Check what integral classes contribute
C     ======================================
C
      ITNEX  = ITMIC + 1
      INTBUF = INTFLG
      CALL INTCON(INTFLG,INTBUF,INTDEF,
     &            CNVMAX,CNVINT,ITNEX,ITRINT,XRS_INTTYP)
      USEXVC = INTFLG.NE.INTBUF
C
C     Check whether calculation has converged
C     =======================================
C
      IF(.NOT.USEXVC) THEN
        KCONV = ISUM(NTSIM,IVECS,1)
        IF(KCONV.EQ.0) GO TO 9999
      ENDIF
C
C     If maximum number of iterations passed, return
C     ==============================================
C
      IF(ITMIC.GE.MAXITR) GO TO 9999
C
C     ***************************************
C     ***** Generate new trial vectors  *****
C     ***** (update IBTYP and IB*VC)    *****
C     ***************************************
C
C
C     If requested, allocate space for diagonal elements of Hessian (approximate)
C
      KFRSAV = KFREE
      IF(DIAGHE.AND.(.NOT.USEXVC)) THEN
         CALL MEMGET('REAL',KEDIAH,NZXOPE,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEDIAS,NZXOPE,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KPDIAH,NZXOPP,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KPDIAS,NZXOPP,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KCDIAH,NZCONF,WORK,KFREE,LFREE)
         CALL E2DIAG(WORK(KEDIAH),WORK(KPDIAH),WORK(KCDIAH),
     &               WORK(KEDIAS),WORK(KPDIAS),
     &               WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KEDIAH,0,WORK,KFREE,LFREE)
         KPDIAH = KEDIAH
         KCDIAH = KEDIAH
         KEDIAS = KEDIAH
         KPDIAS = KEDIAH
      ENDIF
      IF(USEXVC) THEN
         IBOFF  = 0
         NEPREV = 0
         NPPREV = 0
         NCPREV = 0
      ELSE
         IBOFF  = NZRED
         NEPREV = NERED
         NPPREV = NPRED
         NCPREV = NCRED
      ENDIF
C
C     Orbital (e-e) type trial vectors
C
      IF(NZXOPE.GT.0) THEN
        CALL MEMGET('REAL',KPREC,2*NDAMP*NZXOPE,WORK,KFREE,LFREE)
        CALL XBVGEN(BEVEC,USEXVC,WORK(KEDIAH),WORK(KEDIAS),EVALR,EVECR,
     &              IBTYP,IBEVC,IVECS(1,IPEP),
     &              MFREQ,JBENDX,IBOFF,NEPREV,NZXOPE,
     &              BEBUF,WORK(KPREC),WORK,KFREE,LFREE)
        CALL XRSORT(BEVEC,JBENDX,NESIM,NEPREV,
     &              IBEVC,IBTYP,IVECS(1,IPEP),WORK,KFREE,LFREE)
        CALL XRSAVE(BEVEC,JBENDX,NZXOPEQ,NESIM,NEPREV,
     &            IBEVC,IBTYP,IVECS(1,IPEP),WORK,KFREE,LFREE)
        IBOFF = IBOFF + NESIM
      ENDIF
C
C     Orbital (e-p) type trial vectors
C
      IF(NZXOPP.GT.0) THEN
        CALL MEMGET('REAL',KPREC,2*NDAMP*NZXOPP,WORK,KFREE,LFREE)
        CALL XBVGEN(BPVEC,USEXVC,WORK(KPDIAH),WORK(KPDIAS),EVALR,EVECR,
     &              IBTYP,IBPVC,IVECS(1,IPPP),
     &              MFREQ,JBPNDX,IBOFF,NPPREV,NZXOPP,
     &              BPBUF,WORK(KPREC),WORK,KFREE,LFREE)
        CALL XRSORT(BPVEC,JBPNDX,NPSIM,NPPREV,
     &              IBPVC,IBTYP,IVECS(1,IPPP),WORK,KFREE,LFREE)
        CALL XRSAVE(BPVEC,JBPNDX,NZXOPPQ,NPSIM,NPPREV,
     &            IBPVC,IBTYP,IVECS(1,IPPP),WORK,KFREE,LFREE)
        IBOFF = IBOFF + NPSIM
      ENDIF
C
C     Configuration type trial vectors
C
      IF(NZCONF.GT.0) THEN
        KPREC  = KFREE
        KCDIAS = KFREE
        CALL XBVGEN(BCVEC,USEXVC,WORK(KCDIAH),WORK(KCDIAS),EVALR,EVECR,
     &              IBTYP,IBCVC,IVECS(1,IPCP),
     &              MFREQ,JBCNDX,IBOFF,NCPREV,NZCONF,
     &              BCBUF,WORK(KPREC),WORK,KFREE,LFREE)
        CALL XRSORT(BCVEC,JBCNDX,NCSIM,NCPREV,
     &              IBCVC,IBTYP,IVECS(1,IPCP),WORK,KFREE,LFREE)
        CALL XRSAVE(BCVEC,JBCNDX,NZCONFQ,NCSIM,NCPREV,
     &            IBCVC,IBTYP,IVECS(1,IPCP),WORK,KFREE,LFREE)
        IBOFF = IBOFF + NCSIM
      ENDIF
C
C     Check final number of new trial vectors
C     =======================================
C
      NTSIM = NESIM + NPSIM + NCSIM
      IF(NTSIM.EQ.0) THEN
        KCONV = -1
        GO TO 9999
      ENDIF
C
C     Check dimension of reduced space
C     ================================
C
      IF(USEXVC) GO TO 9999
      MREDM = NZRED + NTSIM
      IF(MREDM.GT.NREDM) THEN
        KCONV = -2
        GO TO 9999
      ENDIF
C
 9999 CALL QEXIT('XRSBVC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck XBVGEN */
      SUBROUTINE XBVGEN(BVEC,USEXVC,HDIAG,SDIAG,EVALR,EVECR,IBTYP,IBVEC,
     &            IVECS,MFREQ,JBTYP,IBOFF,NBPREV,NVAR,BBUF,PREC,
     &            WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate new trial vectors b_k from residual vectors Res_k
C
C     b_k = (P')inv Res_k
C
C     where P' is an approximation to E[2] - w S[2].
C
C     Written by T.saue Dec 1 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "thrzer.h"
#include "dcbibn.h"
C
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,DM1 = -1.0D0)
      PARAMETER(DLIM = 1.0D-6)
!     PARAMETER(DLIM = 1.0D-3) sk+hjaaj test
C
#include "dcbxrs.h"
#include "dgroup.h"
C
      LOGICAL USEXVC,TEST
      DIMENSION BVEC(NVAR,NZ,*),HDIAG(*),SDIAG(*),EVECR(*),
     &          EVALR(MFREQ),IBTYP(2,*),IBVEC(*),IVECS(*),
     &          BBUF(NVAR,NZ,*),PREC(NVAR,2,*),WORK(*)
C
      CALL QENTER('XBVGEN')
C
      TEST = DIAGHE.AND.(.NOT.USEXVC)
      FFAC = D1
      IF(IMFREQ) FFAC = DM1
      IF (JBTYP .EQ. JBCNDX) THEN
CHJ-jun2004      prepare for NZ.eq.4 (with MZ)
        MZ = MIN(2,NZ)
      ELSE
         MZ = NZ
      END IF
      NVARQ = NVAR*MZ
C
      IF(USEXVC) THEN
C       New integral class added;
C       start from present solution vector
        CALL XRSXV1(IPX,JBTYP,BVEC,EVECR,MFREQ,
     &        IBTYP,IBVEC,IVECS,BBUF)
        DO IP = 1,MFREQ
          IBTYP(1,IBOFF+IP) = JBTYP
          IBTYP(2,IBOFF+IP) = IPX
          IBVEC(NBPREV +IP) = IBOFF + IP
          IF(IVECS(IP).GT.0) IVECS(IP) = 1
        ENDDO
        IF(.NOT.STATIC) THEN
          CALL XRSXV1(IMX,JBTYP,BVEC(1,1,(MFREQ+1)),EVECR,MFREQ,
     &           IBTYP,IBVEC,IVECS(MFREQ+1),BBUF)
          DO IP = 1,MFREQ
            IM = IP + MFREQ
            FREQ = EVALR(IP)
            IBTYP(1,IBOFF+IM) = JBTYP
            IBTYP(2,IBOFF+IM) = IMX
            IBVEC(NBPREV +IM) = IBOFF + IM
            IF(FREQ.EQ.D0) THEN
               IVECS(IM) = 0
            ELSE
              IF(IVECS(IM).GT.0) IVECS(IM) = 1
            ENDIF
          ENDDO
        ENDIF
        GOTO 10
      ENDIF
C
C     If .not. USEXVC:
C
      IFREQ = 0
      DO IP = 1,MFREQ*NDAMP,NDAMP
        IM = IP + MFREQ*NDAMP
        IFREQ = IFREQ + 1
        FREQ = EVALR(IFREQ)
C
C       For KR-MCSCF optimizations FREQ can be !=0, since it
C       contains the level shift.
C
        IF(FREQ.EQ.D0 .OR. FKRMC) THEN
C         X(+) - vector
          IBTYP(1,IBOFF+IP) = JBTYP
          IBTYP(2,IBOFF+IP) = IPX
          IBVEC(NBPREV +IP) = IBOFF + IP
          IF(TEST.AND.IVECS(IP).EQ.1) THEN
            BNRMP = DNRM2(NVARQ,BVEC(1,1,IP),1)
            IF (JBTYP .EQ. JBCNDX) THEN
              IVAR = IDAMAX(NVARQ,BVEC(1,1,IP),1)
              THRBVEC = ABS(BVEC(IVAR,1,IP))*1.0D-1
CHJ           THRBVEC = BNRMP/SQRT(1.0D2*NVAR)
CHJ           THRBVEC = THCXRS/SQRT(1.0D6*NVAR)
CHJ 16-jun-2004: try use BNRMP instead of THCXRS to zero BCVEC elements
CHJ could this give problems for MZ.eq.2 (e.g. real part included, imag. not)?
            ELSE
              THRBVEC = THRZER
!sk+hjaaj testTHRBVEC = THCXRS*1.0D-4
            END IF
            DO IZ = 1,MZ
              DO IVAR = 1,NVAR
C
C             For MCSCF optimizations FREQ (EVALR) is the
C             level shift. For LR it's zero.
C
                IF ( ABS(BVEC(IVAR,IZ,IP)) .GT. THRBVEC ) THEN
!miro: throwing "forrtl: error (65): floating invalid" for ifort & VAR_MPI
                  IF ( ABS(HDIAG(IVAR)-FREQ*SDIAG(IVAR)) .GT. DLIM) THEN
                       BVEC(IVAR,IZ,IP) =
     &                -BVEC(IVAR,IZ,IP) / (HDIAG(IVAR)-FREQ*SDIAG(IVAR))
                  ELSE
                    BVEC(IVAR,IZ,IP) = -BVEC(IVAR,IZ,IP)
!sk+hjaaj test      BVEC(IVAR,IZ,IP) = -BVEC(IVAR,IZ,IP) / DLIM
                  END IF
                ELSE
                  BVEC(IVAR,IZ,IP) = D0
                END IF
              ENDDO
            ENDDO
            BNRMP = BNRMP/DNRM2(NVARQ,BVEC(1,1,IP),1)
C           retain norm for numerical stability !
            CALL DSCAL(NVARQ,BNRMP,BVEC(1,1,IP),1)
          ENDIF
C         No X(-) - vector
          IF(.NOT.STATIC) THEN
            IBTYP(1,IBOFF+IM) = JBTYP
            IBTYP(2,IBOFF+IM) = IMX
            IBVEC(NBPREV+IM)  = IBOFF + IM
            IVECS(IM)          = 0
          ENDIF
        ELSE
CPP need modification for cpp
          DO IDAMP = 0,NDAMP-1
            IBTYP(1,IBOFF+IP+IDAMP)     = JBTYP
            IBTYP(2,IBOFF+IP+IDAMP)     = IPX
            IBVEC(NBPREV+IP+IDAMP)      = IBOFF + IP + IDAMP
            IBTYP(1,IBOFF+IM+IDAMP)     = JBTYP
            IBTYP(2,IBOFF+IM+IDAMP)     = IMX
            IBVEC(NBPREV+IM+IDAMP)      = IBOFF + IM + IDAMP
          ENDDO
          IF(TEST) THEN
C ------
C           hj-nov99
C           Preconditioning: 
C           For REAL frequencies multiply with minus the elements of
C
C           /    hdiag  -freq*sdiag \ -1      /    hdiag*fac   freq*sdiag*fac \
C           |                       |     =   |                               |
C           \ -freq*sdiag   hdiag   /         \  freq*sdiag*fac     hdiag*fac /
C
C           where fac = (hdiag**2 - freq**2*sdiag2)**(-1)
C           for the res(ip), res(im) pair
C
C           For IMAGINARY frequencies multiply with minus the elements of
C
C           /    hdiag  +freq*sdiag \ -1      /    hdiag*fac  -freq*sdiag*fac \
C           |                       |     =   |                               |
C           \ -freq*sdiag   hdiag   /         \ +freq*sdiag*fac     hdiag*fac /
C
C           where fac = (hdiag**2 + freq**2*sdiag2)**(-1)
C           for the res(ip), res(im) pair
C
C           For COMPLEX frequencies, it is more complicated.
C
C           see eq.(23) and eq.(24) of JCP 133, 064105 (2010)
C
C ------
C           build preconditioner
            CALL GENPREC(PREC,HDIAG,SDIAG,FREQ,NVAR)
C           Note that all residuals are employed, independent of norm
            IF(DAMPFREQ.NE.D0) THEN
              CALL MEMGET('REAL',KPRBUF,2*NDAMP*NVARQ,WORK,KFREE,LFREE)
              CALL XRSPREC(PREC,BVEC(1,1,IP),IVECS(IP),WORK(KPRBUF),
     &                     NVAR,WORK,KFREE,LFREE)
              CALL MEMREL('XRSPREC.cpp',WORK,1,KPRBUF,KFREE,LFREE)
            ELSE
              IF(IVECS(IP).EQ.1) THEN
                DO IZ = 1,MZ
                  DO IVAR = 1,NVAR
                   BBUF(IVAR,IZ,1)=PREC(IVAR,1,1)*BVEC(IVAR,IZ,IP)+
     &                        FFAC*PREC(IVAR,2,1)*BVEC(IVAR,IZ,IM)
                   ENDDO
                ENDDO
CPP what about the norm retained? the new trial vector ++ keeps the norm
C      of the Residual ++?
                FAC=DNRM2(NVARQ,BVEC(1,1,IP),1)/
     &              DNRM2(NVARQ,BBUF(1,1,1),1)
                CALL DSCAL(NVARQ,FAC,BBUF(1,1,1),1)
              ENDIF
              IF(IVECS(IM).EQ.1) THEN
                DO IZ = 1,MZ
                  DO IVAR = 1,NVAR
                    BBUF(IVAR,IZ,2)=PREC(IVAR,1,1)*BVEC(IVAR,IZ,IM)+
     &                              PREC(IVAR,2,1)*BVEC(IVAR,IZ,IP)
                  ENDDO
                ENDDO
                FAC=DNRM2(NVARQ,BVEC(1,1,IM),1)/
     &              DNRM2(NVARQ,BBUF(1,1,2),1)
                CALL DSCAL(NVARQ,FAC,BBUF(1,1,2),1)
              ENDIF
              IF(IVECS(IP).EQ.1) THEN
                CALL DCOPY(NVARQ,BBUF(1,1,1),1,BVEC(1,1,IP),1)
              ENDIF
              IF(IVECS(IM).EQ.1) THEN
                CALL DCOPY(NVARQ,BBUF(1,1,2),1,BVEC(1,1,IM),1)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C
 10   CONTINUE
      CALL QEXIT('XBVGEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrsprec */
      SUBROUTINE XRSPREC(PREC,BVEC,IVECS,PRBUF,NVAR,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate new trial vectors b_k from residual vectors Res_k
C
C     for the complex frequency case.
C
C     Written by S.Villaume Jan 2011
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dgroup.h"
C
      DIMENSION BVEC(NVAR,NZ,*),IVECS(*),PREC(NVAR,2,*),
     &          PRBUF(NVAR,NZ,*),WORK(*)
C
      CALL QENTER('XRSPREC')
C
      NVARQ = NVAR*NZ
      IOFF = NFREQ*NDAMP
C
      CALL DZERO(PRBUF,4*NVAR*NZ)
C
C     build b(+)
      IF(IVECS(1).EQ.1) THEN
        DO IZ = 1,NZ
          DO IVAR = 1,NVAR
            PRBUF(IVAR,IZ,1)=PREC(IVAR,1,1)*BVEC(IVAR,IZ,1)+
     &                       PREC(IVAR,2,1)*BVEC(IVAR,IZ,1+IOFF)-
     &                       PREC(IVAR,2,2)*BVEC(IVAR,IZ,2+IOFF)-
     &                       PREC(IVAR,1,2)*BVEC(IVAR,IZ,2)
          ENDDO
        ENDDO
        FAC=DNRM2(NVARQ,BVEC(1,1,1),1)/
     &      DNRM2(NVARQ,PRBUF(1,1,1),1)
        CALL DSCAL(NVARQ,FAC,PRBUF(1,1,1),1)
      ENDIF
C
C     build b'(+)
      IF(IVECS(2).EQ.1) THEN
        DO IZ = 1,NZ
          DO IVAR = 1,NVAR
            PRBUF(IVAR,IZ,2)=PREC(IVAR,1,2)*BVEC(IVAR,IZ,1)+
     &                       PREC(IVAR,2,2)*BVEC(IVAR,IZ,1+IOFF)+
     &                       PREC(IVAR,2,1)*BVEC(IVAR,IZ,2+IOFF)+
     &                       PREC(IVAR,1,1)*BVEC(IVAR,IZ,2)
          ENDDO
        ENDDO
        FAC=DNRM2(NVARQ,BVEC(1,1,2),1)/
     &      DNRM2(NVARQ,PRBUF(1,1,2),1)
        CALL DSCAL(NVARQ,FAC,PRBUF(1,1,2),1)
      ENDIF
C
C     build b(-)
      IF(IVECS(1+IOFF).EQ.1) THEN
        DO IZ = 1,NZ
          DO IVAR = 1,NVAR
            PRBUF(IVAR,IZ,3)=PREC(IVAR,2,1)*BVEC(IVAR,IZ,1)+
     &                       PREC(IVAR,1,1)*BVEC(IVAR,IZ,1+IOFF)-
     &                       PREC(IVAR,1,2)*BVEC(IVAR,IZ,2+IOFF)-
     &                       PREC(IVAR,2,2)*BVEC(IVAR,IZ,2)
          ENDDO
        ENDDO
        FAC=DNRM2(NVARQ,BVEC(1,1,1+IOFF),1)/
     &      DNRM2(NVARQ,PRBUF(1,1,3),1)
        CALL DSCAL(NVARQ,FAC,PRBUF(1,1,3),1)
      ENDIF
C
C     build b'(-)
      IF(IVECS(2+IOFF).EQ.1) THEN
        DO IZ = 1,NZ
          DO IVAR = 1,NVAR
            PRBUF(IVAR,IZ,4)=PREC(IVAR,2,2)*BVEC(IVAR,IZ,1)+
     &                       PREC(IVAR,1,2)*BVEC(IVAR,IZ,1+IOFF)+
     &                       PREC(IVAR,1,1)*BVEC(IVAR,IZ,2+IOFF)+
     &                       PREC(IVAR,2,1)*BVEC(IVAR,IZ,2)
          ENDDO
        ENDDO
        FAC=DNRM2(NVARQ,BVEC(1,1,2+IOFF),1)/
     &      DNRM2(NVARQ,PRBUF(1,1,4),1)
        CALL DSCAL(NVARQ,FAC,PRBUF(1,1,4),1)
      ENDIF
C
C     copy new vectors from buffer to BVEC
      IF(IVECS(1).EQ.1) THEN
        CALL DCOPY(NVARQ,PRBUF(1,1,1),1,BVEC(1,1,1),1)
      ENDIF
      IF(IVECS(2).EQ.1) THEN
        CALL DCOPY(NVARQ,PRBUF(1,1,2),1,BVEC(1,1,2),1)
      ENDIF
      IF(IVECS(1+IOFF).EQ.1) THEN
        CALL DCOPY(NVARQ,PRBUF(1,1,3),1,BVEC(1,1,1+IOFF),1)
      ENDIF
      IF(IVECS(2+IOFF).EQ.1) THEN
        CALL DCOPY(NVARQ,PRBUF(1,1,4),1,BVEC(1,1,2+IOFF),1)
      ENDIF
C
      CALL QEXIT('XRSPREC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrstdm */
      SUBROUTINE XRSTDM(CREF,BCVEC,DVT,PVT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate transtion density matrices.
C
C     Input:
C       CREF    - reference CI vector
C       BCVEC   - CI trial vectors.
C
C     Output:
C
C       DVT     - one-electron transition densities
C       PVT     - two-electron transition densities
C
C     Written by J. Thyssen - Nov 24 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcborb.h"
C
      DIMENSION WORK(*)
      DIMENSION CREF(*)
      DIMENSION BCVEC(NZCONFQ,NCSIM)
      DIMENSION DVT(N2ASHX,NZ,NCSIM)
      DIMENSION PVT(NASHT*NASHT*NNASHX*NZ*3,NCSIM)
C
      CALL QENTER('XRSTDM')
      CALL RTKTIME(.TRUE.,14)
C
      DO ICSIM = 1, NCSIM
C
C        Generate density matrix
C
         CALL XRSDM(.TRUE., DVT(1,1,ICSIM), PVT(1,ICSIM),
     &        CREF, BCVEC(1,ICSIM),WORK,KFREE,LFREE)
C
C        Symmetrize <0|X|B> + <B|X|0>
C
         CALL FULMAT('S',NASHT,NASHT,DVT(1,1,ICSIM))
         DO IZ = 2, NZ
            CALL FULMAT('A',NASHT,NASHT,DVT(1,IZ,ICSIM))
         END DO
C
         IF ( IPRXRS .GE. 10 ) THEN
            WRITE(LUPRI,'(1X,A,I3)')
     &           '(XRSTDM) DV density matrix no. ',ICSIM
            CALL PRQMAT(DVT(1,1,ICSIM),NASHT,NASHT,NASHT,NASHT,NZ,
     &           IPQTOQ(1,0),LUPRI)
            WRITE(LUPRI,'(1X,A,I3)')
     &           '(XRSTDM) PV density matrix no. ',ICSIM
            CALL PRDNZ3(PVT(1,ICSIM),NASHT,NNASHX,NZ,IPQTOQ(1,0),LUPRI)
         END IF
      END DO
C
      CALL RTKTIME(.FALSE.,14)
      CALL QEXIT('XRSTDM')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrsdm */
      SUBROUTINE XRSDM(TSYMM,DV,PV,CLREF,CRREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate density matrices.
C
C     Input:
C       TSYMM   - symmetrize <L|..|R> + <R|..|L>
C       CLREF, CRREF - left and right CI vectors.
C
C     Output:
C
C       DV      - one-electron transition densities
C       PV      - two-electron transition densities
C
C     Written by J. Thyssen - Nov 24 2000
C     Last revision :
C
C***********************************************************************
      use mcscf_routines
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (D1 = 1.0D00, D0 = 0.0D00)
C
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcborb.h"
#include "dcbham.h"
C
      LOGICAL TSYMM
      DIMENSION WORK(*)
      DIMENSION CLREF(*), CRREF(*)
      DIMENSION DV(*), PV(*)
C
      CALL QENTER('XRSDM')
      KFRSAV = KFREE
C
C
      IF ( NASHT .EQ. 1 ) THEN
         DV(1) = D1
         PV(1) = D1
         DO IZ = 2, NZ
            DV(IZ) = D0
            PV(IZ) = D0
         END DO
         GO TO 100
      END IF
C
      IF ( XRS_CIPROGRAM .EQ. 'GASCIP' .OR.
     &     XRS_CIPROGRAM .EQ. 'LUCIAREL' ) THEN
C
C        Allocate memory
C        ---------------
C
         CALL MEMGET('REAL',KPVMOLF,(2*NASHT)**4 * MIN(2,NZ)
     &        ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KDVMOLF,(2*NASHT)**2 * MIN(2,NZ),
     &        WORK,KFREE,LFREE)
!
         CALL DZERO(WORK(KPVMOLF),(2*NASHT)**4 * MIN(2,NZ))
         CALL DZERO(WORK(KDVMOLF),(2*NASHT)**2 * MIN(2,NZ))
C
C        Calculate density matrices
C        --------------------------
C
         IF ( XRS_CIPROGRAM .EQ. 'GASCIP' ) THEN
            CALL GASCIP_RMAKDM(.TRUE.,.TRUE.,NZCONF,WORK(KZCONF),
     &           CLREF,CRREF,WORK(KDVMOLF),WORK(KPVMOLF),IPRXRS)
         ELSE
C           LUCIAREL:
            CALL LUCI_DENS(.TRUE.,.TRUE.,WORK(KDVMOLF),WORK(KPVMOLF),
     &                     CLREF,CRREF,WORK,KFREE,LFREE)
         END IF
C
C        Transform matrices from Molfdir format to Dirac format
C        ------------------------------------------------------
C
C        One-electron density matrix
C
         CALL DZERO(DV,NASHT**2 * MIN(2,NZ))
!
         CALL MFC2QFC(.TRUE.,WORK(KDVMOLF),DV,IPRXRS)
C        CALL MFC2QFC(KSYMM,DVM,DV,IPRINT)
C
C        Two-electron density matrix
C
         CALL DZERO(PV,NASHT**4 * MIN(2,NZ))
!
         NM5 = MIN(2,NZ)
         CALL M2DNZ390(TSYMM,.TRUE.,.TRUE.,WORK(KPVMOLF),PV,IPRXRS,
     &                 IPQTOQ,NASH,NASHT,NNASHX,NZ,NZ_in_CI,NM5,
     &                 NBSYM,NFSYM,LUPRI,SPINFR,LEVYLE)
C        CALL M2DNZ3(TSYMM,KSYMM,PSYMM,PVM,PV,IPRINT)
C
         CALL MEMREL('XRSDM',WORK,1,KFRSAV,KFREE,LFREE)
C
      END IF
C
 100  CONTINUE
      IF ( IPRXRS .GE. 10 ) THEN
         CALL HEADER('XRSDM: DV density matrix',-1)
         CALL PRQMAT(DV,NASHT,NASHT,NASHT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
         CALL HEADER('XRSDM: PV density matrix',-1)
         CALL PRDNZ3(PV,NASHT,NNASHX,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      CALL QEXIT('XRSDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmolitx */
      SUBROUTINE GMOLITX(FMO,DVT,DV,
     &                   BEVEC,BPVEC,IBCVC,IBEVC,IBPVC,IBTYP,
     &                   NPOS,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Driver routine.
C
C     Calculate myriads of 2-electron Fock matrices "G-matrices".
C
C     Input:
C        DVT    : transition density matrices.
C        DV     : wf. density matrix
C        BEVEC  : electronic trial vectors
C        BPVEC  : positronic trial vectors
C        IBCVC, IBEVC, IBPVC: index stuff
C        NPOS   : something needed for the call to TWOFCK
C
C     Output:
C        FMO    : Fock matrices in MO basis: FMO(N2ORBXQ,NDSIM)
C                 NOTE: FMO is used as scratch space for AO matrices and
C                 it should therefore be allocated with NDSIM*NBBASXQ.
C
C     Written by J. Thyssen - Nov 24 2000
C     Replaces GMOLI by T.Saue Sep 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
C
#include "dcbxrs.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
C
      DIMENSION FMO(*), DV(*), DVT(*)
      DIMENSION BEVEC(*), BPVEC(*)
      DIMENSION IBCVC(*), IBEVC(*), IBPVC(*), IBTYP(*), NPOS(*)
      DIMENSION WORK(*)
C
      LOGICAL TOBE
C
C
      CALL QENTER('GMOLITX')
      CALL RTKTIME(.TRUE.,5)
      KFRSAV = KFREE
C
C
C     Allocate space
C     --------------
C
      NDIM = NTSIM
      IF (NASHT .GT. 0) NDIM = NDIM + NOSIM
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      IF(SPINFR) THEN
        CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
      ELSE
        KBEIG = KFREE
      ENDIF
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAO,NDIM * N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KSYM,NTSIM * NDIM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHRM,NTSIM * NDIM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KFCK,NTSIM * NDIM,WORK,KFREE,LFREE)
C
C     Get coefficients
C     ----------------
C
      IF ( .NOT. FKRMC) THEN
C
C        We're not doing a KRMCSCF optimization
C
         INQUIRE(FILE='DFCOEF',EXIST=TOBE)
         IF (.NOT.TOBE) THEN
            WRITE(LUPRI,'(/A)')
     &      '*** ERROR in GMOLITX : Coefficient file DFCOEF not found!'
            CALL QUIT('GMOLITX: coeff. file DFCOEF file not found!')
         END IF
         LUBUF = 1
         CALL OPNFIL(LUBUF,'DFCOEF','OLD','GMOLI1')
         IF(SPINFR) THEN
           IOPT = 10
         ELSE
           IOPT = 2
         ENDIF
         IF (IPRXRS .GE. 10) IOPT = IOPT + 1
         CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),WORK(KEIG),WORK(KBEIG),
     &               TOTERG,IOPT)
         CLOSE(LUBUF,STATUS='KEEP')
      ELSE
C
C        We are doing a KRMCSCF opt. Read orbitals from KRMCSCF
         JRDMO = 1
         CALL RREADMO(WORK(KCMO),JRDMO,1,LUKRMC)
         call izero(work(kbeig),norbt)
         if(spinfr) call ireakrmc(lukrmc,'IBEIG   ',work(kbeig),norbt)
      ENDIF
C
C     Print MO coefficients ?
C
      IF (IPRXRS .GE. 10) THEN
         CALL HEADER('GMOLITX: MO coefficients',-1)
         DO I = 1, NFSYM
            IF (NFSYM.GT.1) WRITE(LUPRI,'(/A,I2/)')
     &         ' *** MO coefficients of fermion symmetry',I
            CALL PRQMAT(WORK(KCMO+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &           NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
C
C     Calculate modified GX Fock matrix
C     ---------------------------------
C
      CALL GMOLITX1(FMO,DVT,DV,
     &     WORK(KCMO),WORK(KBEIG),BEVEC,BPVEC,IBCVC,IBEVC,IBPVC,IBTYP,
     &     NPOS,WORK(KFAO),WORK(KSYM),WORK(KHRM),WORK(KFCK),
     &     WORK,KFREE,LFREE)
C
      CALL MEMREL('GMOLITX',WORK,1,KFRSAV,KFREE,LFREE)
      CALL RTKTIME(.FALSE.,5)
      CALL QEXIT('GMOLITX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmolitx1 */
      SUBROUTINE GMOLITX1(FMO,DVT,DV,
     &     CMO,IBEIG,BEVEC,BPVEC,IBCVC,IBEVC,IBPVC,IBTYP,
     &     NPOS,FAO,ISYMOP,IHRMOP,IFCKOP,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate myriads of Fock matrices, called by GMOLITX.
C
C     Input:
C        DVT    : transition density matrices.
C        DV     : wf. density matrix
C        BEVEC  : electronic trial vectors
C        BPVEC  : positronic trial vectors
C        IBCVC, IBEVC, IBPVC: index stuff
C        NPOS   : something needed for the call to TWOFCK
C
C        FAO    : scratch used for AO Fock matrices.
C
C     Output:
C        FMO    : Fock matrices.
C
C     Written by J. Thyssen - Nov 24 2000
C
C***********************************************************************
#ifdef MOD_LAO_REARRANGED
        use london_helper
#endif
        use dirac_cfg
        use dft_cfg
        use xcint_main
        use fde_mod
        use fde_dirac_matrices_integration
        use num_grid_gen
#ifdef HAS_PCMSOLVER
        use pcm_linear_response, only: pcm_oit_asc
        use pcmmod_cfg
#endif

#ifdef HAS_PELIB
        use polarizable_embedding, only: pe_master
        use pe_variables, only: peqm, pe_gspol
#endif
        use orbital_rotation_indices


#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
C
#include "cbihr2.h"
#include "dcbxrs.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbgen.h"
C
      DIMENSION FMO(N2ORBXQ,*), DV(*), DVT(*), CMO(*),IBEIG(*)
      DIMENSION BEVEC(*), BPVEC(*), FAO(N2BBASXQ,*)
      DIMENSION IBCVC(*), IBEVC(*), IBPVC(*), IBTYP(*), NPOS(*)
      DIMENSION ISYMOP(*), IHRMOP(*), IFCKOP(*)
      DIMENSION WORK(*)
      type(fde_import) :: itmp
      real(8), allocatable ::mep_oit_asc(:, :, :)

C
C
      CALL QENTER('GMOLITX1')
      IREP = JSYMOP - 1
C
C
C     ************************
C     *** Get AO densities ***
C     ************************
C
C
C     Set IHRMOP, IFCKOP and ISYMOP
C     -----------------------------
C
      CALL SETGMT(ISYMOP,IHRMOP,IFCKOP,IBCVC,IBEVC,IBPVC,IBTYP,IPRXRS)
C
C     Calculate AO density matrices
C     -----------------------------
C
C     Note: reuse memory and save in FMO:
C        DVT(1..NCSIM), [DCX, DVX](1..NESIM), [DCX, DVX](1..NPSIM)
C     Only the necessary matrices are saved.
C
      CALL DMOLITX(FMO,BEVEC,BPVEC,CMO,DVT,DV,IHRMOP,NDMAT,
     &     WORK,KFREE,LFREE)
C
      IF (IPRXRS.GE. 10) WRITE(LUPRI,'(/A,I3)')
     &     ' (GMOLITX1) Number of density matrices:',NDMAT
C
C
C     ****************************
C     *** Get AO Fock matrices ***
C     ****************************
C
C
C     Calculate AO Fock matrices
C     --------------------------
      CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FAO,FMO,NDMAT,NPOS,INTFLG,
     &     IPRTWO,WORK(KFREE),LFREE)

!       construct kohn-sham contribution
!       ================================

        if (dirac_cfg_dft_calculation) then

          if (ncsim > 0) then
            call quit('gmolitx1: mcscf-dft not implemented')
          end if

!         generate standard density matrix
          call memget('REAL', kbuf, n2bbasxq, work, kfree, lfree)
          call genden(work(kbuf), cmo, 1, iprxrs)

!                     buf holds the unperturbed density matrix
!                     fmo is actually the perturbed ao density matrix
!                     fao is really fao (confusing naming reuses memory)

          call generate_num_grid(work(kbuf))
#ifdef VAR_MPI
          if (parcal) call dirac_parctl( XCINT_PAR )
#endif
          call integrate_xc (xc_mat_dim            = ntbas(0),
     &                       xc_nz                 = nz,
     &                       xc_dmat_0             = work(kbuf),
     &                       xc_nr_dmat            = ndmat,
     &                       xc_nr_fmat            = ndmat,
     &                       xc_dmat               = fmo,
     &                       xc_fmat               = fao,
     &                       xc_fmat_pg_sym        = isymop,
     &                       xc_dmat_pg_sym        = isymop,
     &                       xc_dmat_ih_sym        = ihrmop,
     &                       xc_response_order_mo = 1)

          call memrel('gmoli1.dft', work, 1, kbuf, kfree, lfree)

        end if

        if (dirac_cfg_fde_response) then
           call fde_get_import_info(itmp)
           if (itmp%im_frozen) then

              if (.not.dirac_cfg_dft_calculation) then
!         generate standard density matrix
                 call memget('REAL', kbuf, n2bbasxq, work, kfree, lfree)
                 call genden(work(kbuf), cmo, 1, iprxrs)
              endif

              call fde_dirac_set_isymop(isymop)
              call fde_dirac_set_ihrmop(ihrmop)
#ifdef VAR_MPI
          if (parcal) call dirac_parctl( FDE_PAR )
#endif
              call fde_dirac_emb_matrices_via_integration(       
     &                       fde_mat_dim            = ntbas(0),
     &                       fde_nz                 = nz,
     &                       fde_dmat_0             = work(kbuf),
     &                       fde_nr_dmat            = ndmat,
     &                       fde_nr_fmat            = ndmat,
     &                       fde_dmat               = fmo,
     &                       fde_fmat               = fao,
     &                       fde_fmat_pg_sym        = isymop,
     &                       fde_dmat_pg_sym        = isymop,
     &                       fde_dmat_ih_sym        = ihrmop,
     &                       fde_response_order_mo = 1)

              call memrel('gmoli1.dft', work, 1, kbuf, kfree, lfree)

              write (*,*) 'FDE response contributions included'
           else if (dirac_cfg_fde.and..not.dirac_cfg_fde_response) then
              write (*,*) 'FDE response contributions not included'
           endif
        endif
C
C     Add solvent contributions
C     -------------------------
C
      IF(SOLVEN) THEN
         IF (NCSIM .GT. 0) THEN
            WRITE(LUPRI,'(/A)') '*** ERROR in GMOLITX *** '//
     &           'Solvent not implemented for MCSCF'
            CALL QUIT('GMOLITX1: solvent not implemented for MCSCF.')
         END IF
         CALL SOLFCK(FAO,FMO,NDMAT,ESOLVE,ESOLVN,
     &        WORK(KFREE),LFREE,IPRSOL)
      END IF

! Insert safeguards against symmetry and X2C?
#ifdef HAS_PCMSOLVER
        if (dirac_cfg_pcm) then
           if (pcmmod_skipoit) then
              write(lupri,'(a)') '* Skipping formation of '//
     &          'One-Index Transformed electronic polarization charges.'
           else
              write(lupri,'(a)') '* Calling pcm_oit_asc to get '//
     &          'One-Index Transformed electronic polarization charges.'
! fmo is actually the perturbed ao density matrix
! fao is really fao (confusing naming reuses memory)
              allocate(mep_oit_asc(ntbas(0), ntbas(0), nz))
              mep_oit_asc = 0.0d0
              call pcm_oit_asc(mep_oit_asc, fmo, ndmat, 
     &         work(kfree), lfree)
              call daxpy(n2bbasx, -1.0d0, mep_oit_asc, 1, fao, 1)
!           print *, "Contraction of MEP and pertubed ASC: MEP_OIT_ASC"
!           print *, "mep_oit_asc"
!           call prqmat(mep_oit_asc, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
!                       nz, ipqtoq(1,0), 6 
              deallocate(mep_oit_asc)
           endif
        endif 
#endif      
#ifdef HAS_PELIB
        if (peqm .and. .not. pe_gspol) then
           if (IPRXRS .gt. 3) then
              write(lupri,'(/A)') ' Polarizable embedding'//
     &        ' contributions included in linear response.'
           endif
           call pelib_ifc_response(fao, fmo, ndmat, ihrmop)
        endif
#endif

C
C
C     *******************************************
C     *** Transform Fock matrices to MO basis ***
C     *******************************************
C
 9010 FORMAT(/' (GMOLITX1) FVT (AO-basis) no.: ',I2/)
 9013 FORMAT(/' (GMOLITX1) FVT (MO-basis) no.: ',I2/)
 9011 FORMAT(/' (GMOLITX1) ',A,' (AO-basis) no.: ',I2/)
 9012 FORMAT(/' (GMOLITX1) ',A,' (MO-basis) no.: ',I2/)
C
C     FVT Fock matrices
C     -----------------
C
!stefan & radovan: we decide to set jrep to zero
!                  otherwise undefined
      JREP = 0
      IAOMAT = 0
      IMOMAT = 0
      DO I = 1, NCSIM
         IAOMAT = IAOMAT + 1
         IMOMAT = IMOMAT + 1
         CALL DZERO(FMO(1,IMOMAT),N2ORBXQ)
         CALL GMOLITX2(-1,JREP,FMO(1,IMOMAT),FAO(1,IAOMAT),CMO,IBEIG,
     &        WORK(KFREE),LFREE)
         IF ( IPRXRS .GE. 10) THEN
            WRITE(LUPRI,9010) I
            CALL PRQMAT(FAO(1,IAOMAT),NTBAS(0),NTBAS(0),
     &           NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,9013) I
            CALL PRQMAT(FMO(1,IMOMAT),NORBT,NORBT,
     &           NORBT,NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
         END IF
      END DO
C
C     e-e FCX- and FVX- Fock matrices
C     -------------------------------
C
C     Increment IMAT for FCX's since there is always a FCX
C     independent of NISHT.
C
      JZ    = -1
      IBOFF =  1
      DO I = 1, NESIM
         IF(SPINFR .AND. (JTIMOP .EQ. 1)) THEN
           CALL SPFEXC(BEVEC(IBOFF),
     &                 get_orbital_rotation_indices_pp(),
     &                 NZXOPE,IBEIG,JZ,JREP)
           IBOFF = IBOFF + NZXOPEQ
         ENDIF
          IMOMAT = IMOMAT + 1
         CALL DZERO(FMO(1,IMOMAT),N2ORBXQ)
         IF (NISHT .GT. 0) THEN
            IAOMAT = IAOMAT + 1
!radovan: jrep seems to be irep of kappa
            CALL GMOLITX2(JZ,JREP,FMO(1,IMOMAT),FAO(1,IAOMAT),CMO,IBEIG,
     &           WORK(KFREE),LFREE)
            IF ( IPRXRS .GE. 10) THEN
               WRITE(LUPRI,9011) 'e-e inactive QCX', I
               CALL PRQMAT(FAO(1,IAOMAT),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,9012) 'e-e inactive QCX', I
               CALL PRQMAT(FMO(1,IMOMAT),NORBT,NORBT,
     &              NORBT,NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
            END IF
         END IF
         IF (NASHT .GT. 0) THEN
            IMOMAT = IMOMAT + 1
            IAOMAT = IAOMAT + 1
            CALL DZERO(FMO(1,IMOMAT),N2ORBXQ)
            CALL GMOLITX2(JZ,JREP,FMO(1,IMOMAT),FAO(1,IAOMAT),CMO,IBEIG,
     &           WORK(KFREE),LFREE)
            IF ( IPRXRS .GE. 8) THEN
               WRITE(LUPRI,9011) 'e-e active FVX', I
               CALL PRQMAT(FAO(1,IAOMAT),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,9012) 'e-e active FVX', I
               CALL PRQMAT(FMO(1,IMOMAT),NORBT,NORBT,
     &              NORBT,NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
            END IF
         END IF
      END DO
C
C     e-p FCX- and FVX- Fock matrices
C     -------------------------------
C
      JZ    = -1
      IBOFF =  1
      DO I = 1, NPSIM
         IF(SPINFR .AND. (JTIMOP .EQ. 1)) THEN
           CALL SPFEXC(BPVEC(IBOFF),get_orbital_rotation_indices_pn(),
     &                 NZXOPP,IBEIG,JZ,JREP)
           IBOFF = IBOFF + NZXOPPQ
         ENDIF
         IMOMAT = IMOMAT + 1
         CALL DZERO(FMO(1,IMOMAT),N2ORBXQ)
         IF (NISHT .GT. 0) THEN
            IAOMAT = IAOMAT + 1
            CALL GMOLITX2(JZ,JREP,FMO(1,IMOMAT),FAO(1,IAOMAT),CMO,IBEIG,
     &           WORK(KFREE),LFREE)
            IF ( IPRXRS .GE. 10) THEN
               WRITE(LUPRI,9011) 'e-p inactive QCX', I
               CALL PRQMAT(FAO(1,IAOMAT),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,9012) 'e-p inactive QCX', I
               CALL PRQMAT(FMO(1,IMOMAT),NORBT,NORBT,
     &              NORBT,NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
            END IF
         END IF
         IF (NASHT .GT. 0) THEN
            IMOMAT = IMOMAT + 1
            IAOMAT = IAOMAT + 1
            CALL DZERO(FMO(1,IMOMAT),N2ORBXQ)
            CALL GMOLITX2(-1,JREP,FMO(1,IMOMAT),FAO(1,IAOMAT),CMO,IBEIG,
     &           WORK(KFREE),LFREE)
            IF ( IPRXRS .GE. 8) THEN
               WRITE(LUPRI,9011) 'e-p active FVX', I
               CALL PRQMAT(FAO(1,IAOMAT),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,9012) 'e-p active FVX', I
               CALL PRQMAT(FMO(1,IMOMAT),NORBT,NORBT,
     &              NORBT,NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
            END IF
         END IF
      END DO
C
      CALL QEXIT('GMOLITX1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dmolitx */
      SUBROUTINE DMOLITX(DMAT,BEVEC,BPVEC,CMO,DVT,DV,IHRMOP,
     &     NDMAT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Driver routine.
C
C     Calculate various density matrices.
C
C     Written by J. Thyssen - Nov 24 2000
C     Replaces DKMOLI by T.Saue Aug 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dcborb.h"
#include "dcbbas.h"
C
      DIMENSION DMAT(*), BEVEC(*), BPVEC(*), CMO(*), DVT(*), DV(*)
      DIMENSION IHRMOP(*), WORK(*)
C
      CALL QENTER('DMOLITX')
      KFRSAV = KFREE
C
      CALL MEMGET('REAL',KBMAT,N2ORBXQ,WORK,KFREE,LFREE)
      CALL DMOLITX1(DMAT,BEVEC,BPVEC,CMO,DVT,DV,WORK(KBMAT),
     &            IHRMOP,NDMAT,WORK,KFREE,LFREE)
C
      CALL MEMREL('DMOLITX',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DMOLITX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dmolitx1 */
      SUBROUTINE DMOLITX1(DMAT,BEVEC,BPVEC,CMO,DVT,DV,BMAT,IHRMOP,
     &     NDMAT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate various backtransformed density matrices
C     for use in generation of orbital part of sigma vectors.
C
C     Input:
C
C     Output:
C        NDMAT - number of density matrices
C        DMAT  - the density matrices
C
C     Written by J. Thyssen - Nov 24 2000
C     Replaces DMOLI2 by T. Saue Aug 1996
C
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcbxrs.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION DMAT(N2BBASXQ,*)
      DIMENSION BEVEC(*), BPVEC(*), CMO(*)
      DIMENSION DVT(*), DV(NASHT,NASHT,NZ), BMAT(*)
      DIMENSION IHRMOP(*), WORK(*)
C
      IM = 0
      IREP = JSYMOP - 1
C
      KFRSAV = KFREE
C
C
C     ***************************************************
C     *** Generate backtransformed DVT density matrix ***
C     ***************************************************
C
      IDVTOFF = 1
      DO I = 1, NCSIM
C
         ITYP  = 0*1+1*2
         IM = IM + 1
         CALL DENSDV(DMAT(1,IM),CMO,DVT(IDVTOFF),IPRXRS,
     &        WORK(KFREE),LFREE)
         CALL DSCAL(N2BBASXQ,DP5,DMAT(1,IM),1)
#ifdef UNDEF
C
C        Expand DVT into BMAT
C        --------------------
C
         CALL DZERO(BMAT, N2ORBXQ)
         CALL DVSCT(DVT(IDVTOFF),BMAT)
         IF (IPRXRS .GT. 10) THEN
            WRITE(LUPRI,'(/A,I5)')
     &           '* DMOLITX1: Full DVT matrix no.:',I
            CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,IREP),LUPRI)
         END IF
C
C        Construct backtransformed DVT
C        -----------------------------
C
C        The matrix is scaled with 0.5 in DMOLI3.
C
         CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT)
C
#endif
         IF (IPRXRS .GT. 10) THEN
            WRITE(LUPRI,'(/A,I5)')
     &           '* DMOLITX1: DVT(AO) matrix no.:',I
            CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &           NTBAS(0),NTBAS(0),NZ,
     &           IPQTOQ(1,IREP),LUPRI)
         END IF
C
         IDVTOFF = IDVTOFF + N2ASHXQ
      END DO
C
C
C     **********************************************************
C     *** Generate modified density matrices ; e-e rotations ***
C     **********************************************************
C
C
      IBOFF = 1
      DO I = 1,NESIM
C
C        Inactive density matrices
C        -------------------------
C
C        DCX_{pi} = \sum_r D_{ri} B_{rp}
C                 = \sum_r \delta_{ri} B_{rp}
C                 = 2 B_{ip}
C
C        The factor 2 is multiplied onto the Fock matrix in TWOFCK.
C
C        The projection onto \delta_{ri} is done implicitly by
C        DMOLI3 in the MO-AO transformation.
C
C
         IF (NISHT .GT. 0) THEN
            ITYP  = 0*1+0*2
            IM = IM + 1
C
C           Scatter BVEC into BMAT
C
            CALL DZERO(BMAT,N2ORBXQ)
            CALL XOPSCT(BEVEC(IBOFF),
     &                  get_orbital_rotation_indices_pp(),
     &                  NZXOPE,
     &           BMAT,NORBT,NZ)
C
C           Take symmetric/antisymmetric combinations to form BMAT
C
            IH = IHRMOP(IM)
            IZOFF = 1
            DO IZ = 1,NZ
               IQ = IPQTOQ(IZ,IREP)
               IS = IHQMAT(IQ,IH)
               IF(IS.EQ.1) THEN
                  CALL FULMAT('A',NORBT,NORBT,BMAT(IZOFF))
               ELSEIF(IS.EQ.2) THEN
                  CALL FULMAT('S',NORBT,NORBT,BMAT(IZOFF))
               ENDIF
               IZOFF = IZOFF + N2ORBX
            ENDDO
            CALL DSCAL(N2ORBXQ,DM1,BMAT,1)
C
            IF(IPRXRS.GE.10) THEN
               WRITE(LUPRI,'(A,I3,A,I3)')
     &              '* DMOLITX1: Scattered inactive (e-e)'//
     &              ' trial vector no.:',I,' of',NESIM
               CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
            END IF
C
C           Backtransform to MO basis and (anti)symmetrize.
C
            CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT)
C
            IF(IPRXRS.GE.10) THEN
               WRITE(LUPRI,'(3(A,I3))')
     &              '* DMOLITX1: Modified inactive (e-e) matrix no.:',I,
     &              ' of',NESIM,' - hermiticity:',IH
               CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
            END IF
         END IF
C
C
C        Active density matrices
C        -----------------------
C
C
C        DVX_{pv} = \sum_u B_{up} D_{uv}
C                 = \sum_u (B^T)_{pu} D_{uv}
C
C        DMOLI3 has conveniently divided D with the factor of 2 that TWOFCK
C        multiplies onto Fock matrices.
C
         IF (NASHT .GT. 0) THEN
C
            ITYP  = 0*1+1*2
            IM = IM + 1
C
C           Scatter BVEC into scratch area
C
            CALL MEMGET('REAL',KBTMP,N2ORBXQ,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KDVF,N2ORBXQ,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KBTMP),N2ORBXQ)
            CALL XOPSCT(BEVEC(IBOFF),get_orbital_rotation_indices_pp(),
     &                  NZXOPE,WORK(KBTMP),
     &           NORBT,NZ)
C
C           Take symmetric/antisymmetric combinations to form BMAT
C
            IH = IHRMOP(IM)
            IZOFF = 0
            DO IZ = 1,NZ
               IQ = IPQTOQ(IZ,IREP)
               IS = IHQMAT(IQ,IH)
               IF(IS.EQ.1) THEN
                  CALL FULMAT('A',NORBT,NORBT,WORK(KBTMP+IZOFF))
               ELSEIF(IS.EQ.2) THEN
                  CALL FULMAT('S',NORBT,NORBT,WORK(KBTMP+IZOFF))
               ENDIF
               IZOFF = IZOFF + N2ORBX
            ENDDO
C
            IF(IPRXRS.GE.10) THEN
               WRITE(LUPRI,'(A,I5)')
     &              '* DMOLITX1: Scattered active (e-e)'//
     &              ' trial vector no.:',I
               CALL PRQMAT(WORK(KBTMP),NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
            END IF
C
C           Contract with active density matrix.
C
C           FIXME: The formulae need to be checked... (jth)
C           FIXME: run over symmetries instead of transforming the
C           whole matrix.
C
            CALL DVSCT(DV,WORK(KDVF))
#ifdef UNDEF
            CALL QGEMM(NORBT,NORBT,NORBT,D1,
     &           'T','N',IPQTOQ(1,0),WORK(KDVF),NORBT,NORBT,NZ,
     &           'N','N',IPQTOQ(1,0),WORK(KBTMP),NORBT,NORBT,NZ,
     &           D0,IPQTOQ(1,0),BMAT,NORBT,NORBT,NZ)
            CALL MEMREL('active e-e',WORK,1,KFRSAV,KFREE,LFREE)
#else
            CALL QGEMM(NORBT,NORBT,NORBT,D1,
     &           'T','N',IPQTOQ(1,0),WORK(KBTMP),NORBT,NORBT,NZ,
     &           'N','N',IPQTOQ(1,0),WORK(KDVF),NORBT,NORBT,NZ,
     &           D0,IPQTOQ(1,0),BMAT,NORBT,NORBT,NZ)
            CALL MEMREL('active e-e',WORK,1,KFRSAV,KFREE,LFREE)
#endif
C
C           Take symmetric/antisymmetric combinations to form BMAT
C
#ifdef UNDEF
            IH = IHRMOP(IM)
            IZOFF = 1
            DO IZ = 1,NZ
               IQ = IPQTOQ(IZ,IREP)
               IS = IHQMAT(IQ,IH)
               IF(IS.EQ.1) THEN
                  CALL FULMAT('A',NORBT,NORBT,BMAT(IZOFF))
               ELSEIF(IS.EQ.2) THEN
                  CALL FULMAT('S',NORBT,NORBT,BMAT(IZOFF))
               ENDIF
               IZOFF = IZOFF + N2ORBX
            ENDDO
#endif
C
            IF(IPRXRS.GE.10) THEN
               WRITE(LUPRI,'(A,I5)')
     &              '* DMOLI: DV * B (e-e) vector no.:',I
               CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
            END IF
C
C           Backtransform to AO basis
C
            CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT)
C
C           Print section
C
            IF(IPRXRS.GE.10) THEN
               WRITE(LUPRI,'(A,I5)')
     &              '* DMOLI: Modified active (e-e) matrix no.:',I
               CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
            END IF
         END IF
         IBOFF = IBOFF + NZXOPEQ
C
      ENDDO
C
C
C     **********************************************************
C     *** Generate modified density matrices ; e-p rotations ***
C     **********************************************************
C
C
      IBOFF = 1
      DO I = 1,NPSIM
C
C        Inactive density matrices
C        -------------------------
C
         IF (NISHT .GT. 0) THEN
C
            ITYP = 1*1 + 0*2
            IM = IM + 1
C
C           Scatter BVEC into BMAT
C
            CALL DZERO(BMAT,N2ORBXQ)
            CALL XOPSCT(BPVEC(IBOFF),get_orbital_rotation_indices_pn(),
     &                  NZXOPP,
     &           BMAT,NORBT,NZ)
C
C           Take symmetric/antisymmetric combinations to form BMAT
C
            IH = IHRMOP(IM)
            IZOFF = 1
            DO IZ = 1,NZ
               IQ = IPQTOQ(IZ,IREP)
               IS = IHQMAT(IQ,IH)
               IF(IS.EQ.1) THEN
                  CALL FULMAT('A',NORBT,NORBT,BMAT(IZOFF))
               ELSEIF(IS.EQ.2) THEN
                  CALL FULMAT('S',NORBT,NORBT,BMAT(IZOFF))
               ENDIF
               IZOFF = IZOFF + N2ORBX
            ENDDO
            CALL DSCAL(N2ORBXQ,DM1,BMAT,1)
C
            IF(IPRXRS.GE.10) THEN
               WRITE(LUPRI,'(A,I3,A,I3)')
     &              '* DMOLITX1: Scattered inactive'//
     &              ' (e-p) trial vector no.:',I,' of',NPSIM
               CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
            END IF
C
C           Backtransform to AO basis
C
            CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT)
C
            IF(IPRXRS.GE.10) THEN
               WRITE(LUPRI,'(A,I5)')
     &              '* DMOLITX1: Modified inactive (e-p) matrix no.:',I
               CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
            END IF
         END IF
C
C
C        Active density matrices
C        -----------------------
C
         IF (NASHT .GT. 0) THEN
C
            ITYP  = 1*1+1*2
            IM = IM + 1
C
C           Scatter BVEC into scratch area
C
            CALL MEMGET('REAL',KBTMP,N2ORBXQ,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KDVF,N2ORBXQ,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KBTMP),N2ORBXQ)
            CALL XOPSCT(BPVEC(IBOFF),get_orbital_rotation_indices_pn(),
     &                  NZXOPP,WORK(KBTMP),
     &           NORBT,NZ)
C
C           Take symmetric/antisymmetric combinations to form BMAT
C
            IH = IHRMOP(IM)
            IZOFF = 0
            DO IZ = 1,NZ
               IQ = IPQTOQ(IZ,IREP)
               IS = IHQMAT(IQ,IH)
               IF(IS.EQ.1) THEN
                  CALL FULMAT('A',NORBT,NORBT,WORK(KBTMP+IZOFF))
               ELSEIF(IS.EQ.2) THEN
                  CALL FULMAT('S',NORBT,NORBT,WORK(KBTMP+IZOFF))
               ENDIF
               IZOFF = IZOFF + N2ORBX
            ENDDO
C
            IF(IPRXRS.GE.10) THEN
               WRITE(LUPRI,'(A,I5)')
     &              '* DMOLITX1: Scattered (e-p)'//
     &              ' trial vector no.:',I
               CALL PRQMAT(WORK(KBTMP),NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
            END IF
C
C           Contract with active density matrix.
C
C           FIXME: The formulae need to be checked... (jth)
C           FIXME: run over symmetries instead of transforming the
C           whole matrix.
C
            CALL DVSCT(DV,WORK(KDVF))
            CALL QGEMM(NORBT,NORBT,NORBT,D1,
     &           'T','N',IPQTOQ(1,0),WORK(KBTMP),NORBT,NORBT,NZ,
     &           'N','N',IPQTOQ(1,0),WORK(KDVF),NORBT,NORBT,NZ,
     &           D0,IPQTOQ(1,0),BMAT,NORBT,NORBT,NZ)
            CALL MEMREL('active e-p',WORK,1,KFRSAV,KFREE,LFREE)
C
            IF(IPRXRS.GE.10) THEN
               WRITE(LUPRI,'(A,I5)')
     &              '* DMOLITX1: DV * B (e-p) vector no.:',I
               CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
            END IF
C
C
C           Backtransform to AO basis
C
            CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT)
C
C           Print section
C
            IF(IPRXRS.GE.10) THEN
               WRITE(LUPRI,'(A,I5)')
     &              '* DMOLITX1: Modified active (e-p) matrix no.:',I
               CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
            END IF
         END IF
         IBOFF = IBOFF + NZXOPPQ
C
C
      ENDDO
C
      NDMAT = IM
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dvsct */
      SUBROUTINE DVSCT(DV,DVF)
C***********************************************************************
C
C     Expand N2ASHX DV into N2ORBXQ DVF
C
C     Input:
C        DV - active density matrix with dimension N2ASHX
C
C     Output:
C        DVF - active density matrix with dimenion N2ORBX
C
C     Written by J. Thyssen - Nov 24 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "maxorb.h"
#include "dcbidx.h"
C
      DIMENSION DV(NASHT,NASHT,NZ)
      DIMENSION DVF(NORBT,NORBT,NZ)
C
      CALL DZERO(DVF,N2ORBXQ)
      DO IZ = 1, NZ
         DO J = 1, NASHT
            JG = IDXU2G(J)
            DO I = 1, NASHT
               IG = IDXU2G(I)
               DVF(IG,JG,IZ) = DV(I,J,IZ)
            END DO
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmolitx2 */
      SUBROUTINE GMOLITX2(IOPT,JREP,FMO,FAO,CMO,IBEIG,WORK,LWORK)
C***********************************************************************
C
C     Transform the G-matrix to MO-basis.
C
C     Input:
C     - FAO:  the Fock matrix in AO basis
C     - CMO:  the MO-AO transformation matrix
C     - IOPT: in the spinfree time-symmetric case IOPT.GT.0
C             indicates that QBTRANS should be used, that is
C             SO effects are quenched by only transforming
C             to specific symmetry-blocks of the Fock
C             matrix in MO-basis. IOPT then indicates what
C             component of FMO to transform, whereas
C     - JREP: indicates what symmetry-blocks to transform to
C
C     Output:
C     - FMO: the Fock matrix in MO basis
C
C     Written by J. Thyssen - Nov 27 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcbham.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbxrs.h"
#include "dgroup.h"
C
      DIMENSION FMO(NORBT,NORBT,NZ), CMO(*), IBEIG(*), FAO(*), WORK(*)
C
#include "memint.h"
C
      IREP  = JSYMOP-1
      ITIM  = JTIMOP
      DO 10 I1 = 1,NFSYM
         I2  = MOD(I1+JOPSY,2) + 1
C
C        Parameters for first index
C
C        For historical (read: copy'n'paste) reasons the first
C        index is named NI as in the original version the first
C        index was over inactive orbitals. Today, it stands for number
C        of occupIed orbitals.
C
         NI = NOCC(I2)
         II = IORB(I2) + NPSH(I2) + 1
         JI = ICMOQ(I2) + NPSH(I2) * NFBAS(I2,0) + 1
C
C        Analogous for the second index, which is named NS for number
C        of secondary orbitals. Today, it stands for Sotal number of orbitals.C
C
         NS = NORB(I1)
         IS = IORB(I1) + 1
         JS = ICMOQ(I1) + 1
C
C        Transform the d*** thing:
C
C        In the original routine (GMOLI3) you would transform
C        the electronic part and the positronic part separately.
C        This saves some FLOPs when there is electronic rotations, but
C        no positronic rotations (and vice versa). However, it's just simpler
C        to let the second index run over all orbitals.
C
         IF(NS.NE.0.AND.NI.NE.0) THEN
           iopt_spinfree = IOPT
#ifdef MOD_MCSCF_spinfree
!radovan:  the following line breaks .SPINFREE and .LEVY-LEBLOND
!          magnetic response
!          if (SPINFR) iopt_spinfree = 1
!          rather do:
           if (SPINFR.AND.((ITIM.EQ.1).OR.NOSPIN)) iopt_spinfree = 1
#endif
          IF(IOPT.GT.0) THEN
!           IF(iopt_spinfree.GT.0) THEN
             CALL QBTRANS(JREP,'AOMO','S',D0,
     &           NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &           FAO(I2BASX(I1,I2)+1),NTBAS(0),NTBAS(0),NZ,
     &           IPQTOQ(1,IREP),
C     &           FMO(IS,II,iopt_spinfree),NORBT,NORBT,1,
C     &           IPQTOQ(iopt_spinfree,IREP),
     &           FMO(IS,II,iopt),NORBT,NORBT,1,
     &           IPQTOQ(iopt,IREP),
     &           CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &           IBEIG(IS),
     &           CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &           IBEIG(II),
     &           WORK(KFREE),LFREE,IPRXRS)
             DO IZ = 1,NZ
CTROND             IF(IZ.NE.iopt_spinfree) THEN
             IF(IZ.NE.iopt) THEN
               CALL DZERO(FMO(1,1,IZ),N2ORBX)
             ENDIF  
             ENDDO
           ELSE
             CALL QTRANS('AOMO','S',D0,
     &           NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &           FAO(I2BASX(I1,I2)+1),NTBAS(0),NTBAS(0),NZ,
     &           IPQTOQ(1,IREP),
     &           FMO(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &           CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &           CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRXRS)
           ENDIF
         ENDIF
 10   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck rsigoc */
      SUBROUTINE RSIGOC(SEVEC,SPVEC,DVT,FVT,FQT,
     &     JXOPE,NXOPE,JXOPP,NXOPP,
     &     IBCVC,IBTYP,WORK,LWORK)
C***********************************************************************
C
C     Calculate orbital part of conf. sigma vectors.
C
C     Input:
C        FVT, FQT    : transition Fock matrices
C        DVT         : one-electron active transition density matrix
C        JXOPE, JXOPP: orbital rotations.
C                      Note that JXOPE should include the redundant
C                      active-active rotations.
C        IB*         : index vectors
C
C     Output:
C        SEVEC       : e-e orbital part of conf. sigma vector(s)
C        SPVEC       : e-p orbital part of conf. sigma vector(s)
C
C     Written by J. Thyssen - Nov 28 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbxrs.h"
      DIMENSION SEVEC(NXOPE,NZ,*),SPVEC(NXOPP,NZ,*)
      DIMENSION FVT(N2ORBXQ,*), FQT(NORBT*NASHT*NZ,*)
      DIMENSION DVT(N2ASHXQ,*)
      DIMENSION JXOPE(*), JXOPP(*)
      DIMENSION WORK(*)
      DIMENSION IBCVC(*),IBTYP(2,*)
C
      DIMENSION DFAC(2)
C
#include "memint.h"
C
      CALL QENTER('RSIGOC')
      CALL RTKTIME(.TRUE.,16)
C
      LFPV = NZ*NASHT*NORBT
      CALL MEMGET2('REAL','FPV',KFPV,LFPV,WORK,KFREE,LFREE)
C
      IREP = JSYMOP - 1
      IM = 0
C
      CALL MEMGET2('REAL','FC',KFC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL GETMAT('FCMO',WORK(KFC),IPRXRS,FKRMC,WORK(KFREE),LFREE)
C
C
      DO I = 1, NCSIM
C
         IM = IM + 1
C
C        Hermiticity (used in RORBGRD)
C        -----------------------------
C
         IS = IBCVC(NCRED+I)
         IS = IBTYP(2,IS)
         DFAC(1) = DM1 * D2 * DM1
         IF (IS .EQ. -1) THEN
            DFAC(2) = DM1 * D2 * DM1
         ELSE
            DFAC(2) = D1 * D2 * DM1
         END IF
C
C        ***********************************
C        *** Calculate e-e sigma vectors ***
C        ***********************************
C
         IF (NXOPE .GT. 0) THEN
C
C           Calculate e-e sigma vector as gradient
C           --------------------------------------
C
            CALL DZERO(WORK(KFPV),LFPV)
            CALL RORBGRD(.FALSE.,DVT(1,IM),WORK(KFC),DUMMY,
     &           FVT(1,IM),FQT(1,IM),
     &           SEVEC(1,1,I+NOSIM),JXOPE,NXOPE,WORK(KFPV),DFAC,
     &           IPQTOQ(1,IREP),IPRXRS)
C
         END IF
C
C        ***********************************
C        *** Calculate e-p sigma vectors ***
C        ***********************************
C
         IF (NXOPP .GT. 0) THEN
C
C           Calculate e-p sigma vector as gradient
C           --------------------------------------
C
            CALL DZERO(WORK(KFPV),LFPV)
            CALL RORBGRD(.FALSE.,DVT(1,IM),WORK(KFC),DUMMY,
     &           FVT(1,IM),FQT(1,IM),
     &           SPVEC(1,1,I+NOSIM),JXOPP,NXOPP,WORK(KFPV),DFAC,
     &           IPQTOQ(1,IREP),IPRXRS)
C
         END IF
      END DO
C
      IF(IPRXRS.GE.5) THEN
         CALL HEADER('Output from RSIGOC',-1)
      IF (NXOPE .GT. 0) THEN
         WRITE(LUPRI,'(/I5,A)') NCSIM,
     &        ' e-e orbital part of conf. sigma vectors'
         CALL PRBVEC(LUPRI,SEVEC(1,1,1+NOSIM),NCSIM,NXOPE)
      END IF
      IF (NXOPP .GT. 0) THEN
         WRITE(LUPRI,'(/I5,A)') NCSIM,
     &        ' e-p orbital part of conf. sigma vectors'
         CALL PRBVEC(LUPRI,SPVEC(1,1,1+NOSIM),NCSIM,NXOPP)
      END IF
      ENDIF
C
      CALL MEMREL('RSIGOC',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL RTKTIME(.FALSE.,16)
      CALL QEXIT('RSIGOC')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck rsigoo */
      SUBROUTINE RSIGOO(SEVEC,SPVEC,FX,FQX,DV,
     &                 JXOPE,NXOPE,JXOPP,NXOPP,
     &                 IBEVC,IBPVC,IBTYP,WORK,LWORK)
C***********************************************************************
C
C     Calculate orbital part of orbital sigma vector(s).
C
C     Written by J. Thyssen - Nov 29 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbxrs.h"
#include "dummy.h"
C
      DIMENSION SEVEC(NXOPE,NZ,*), SPVEC(NXOPP,NZ,*)
      DIMENSION FX(N2ORBXQ,*), FQX(NASHT*NORBT*NZ,*), DV(*)
      DIMENSION JXOPE(*), JXOPP(*), IBEVC(*),IBPVC(*),IBTYP(2,*)
      DIMENSION WORK(*)
C
      DIMENSION DFAC(2)
C
#include "memint.h"
C
      CALL QENTER('RSIGOO')
      CALL RTKTIME(.TRUE.,16)
      CALL MEMGET('REAL',KFPV,NZ*NASHT*NORBT,WORK,KFREE,LFREE)
C
C
C     *******************************
C     *** Calculate sigma vectors ***
C     *******************************
C
      IREP = JSYMOP - 1
      IM = 0
C
C
C     ***************************
C     *** e-e trial vector(s) ***
C     ***************************
C
C
      DO I = 1, NESIM
C
         IM = IM + 1
C
C        Hermicity (used in RORBGRD)
C        ---------------------------
C
         IS = IBEVC(NERED+I)
         IS = IBTYP(2,IS)
         DFAC(1) = -1.0D00
         IF (IS .EQ. -1) THEN
            DFAC(2) = -1.0D00
         ELSE
            DFAC(2) = +1.0D00
         END IF
C
C        e-e part of sigma vector
C        ------------------------
C
         IF (NXOPE .GT. 0) THEN
           IF (NASHT.NE.0) THEN
            CALL RORBGRD(.TRUE.,DV,FX(1,IM),FX(1,IM),FX(1,IM+1),
     &           FQX(1,I),SEVEC(1,1,I),JXOPE,NXOPE,WORK(KFPV),DFAC,
     &           IPQTOQ(1,IREP),IPRXRS)
           ELSE
            CALL RORBGRD(.TRUE.,DV,FX(1,IM),FX(1,IM),FX(1,IM+1),
     &           DUMMY,SEVEC(1,1,I),JXOPE,NXOPE,WORK(KFPV),DFAC,
     &           IPQTOQ(1,IREP),IPRXRS)
           ENDIF
         END IF
C
C        e-p part of sigma vector
C        ------------------------
C
         IF (NXOPP .GT. 0) THEN
           IF (NASHT.NE.0) THEN
            CALL RORBGRD(.TRUE.,DV,FX(1,IM),FX(1,IM),FX(1,IM+1),
     &           FQX(1,I),SPVEC(1,1,I),JXOPP,NXOPP,WORK(KFPV),DFAC,
     &           IPQTOQ(1,IREP),IPRXRS)
           ELSE
            CALL RORBGRD(.TRUE.,DV,FX(1,IM),FX(1,IM),FX(1,IM+1),
     &           DUMMY,SPVEC(1,1,I),JXOPP,NXOPP,WORK(KFPV),DFAC,
     &           IPQTOQ(1,IREP),IPRXRS)
           ENDIF
         END IF
C
         IF (NASHT .GT. 0) IM = IM + 1
C
      END DO
C
C
C     ***************************
C     *** e-p trial vector(s) ***
C     ***************************
C
C
      DO I = 1, NPSIM
C
         IM = IM + 1
C
C        Hermicity (used in RORBGRD)
C        ---------------------------
C
         IS = IBPVC(NPRED+I)
         IS = IBTYP(2,IS)
         DFAC(1) = -1.0D00
         IF (IS .EQ. -1) THEN
            DFAC(2) = -1.0D00
         ELSE
            DFAC(2) = +1.0D00
         END IF
C
C        e-e part of sigma vector
C        ------------------------
C
         IF (NXOPE .GT. 0) THEN
           IF (NASHT.NE.0) THEN
            CALL RORBGRD(.TRUE.,DV,FX(1,IM),FX(1,IM),FX(1,IM+1),
     &           FQX(1,I+NESIM),
     &           SEVEC(1,1,I+NESIM),JXOPE,NXOPE,WORK(KFPV),DFAC,
     &           IPQTOQ(1,IREP),IPRXRS)
           ELSE
            CALL RORBGRD(.TRUE.,DV,FX(1,IM),FX(1,IM),FX(1,IM+1),
     &           DUMMY,
     &           SEVEC(1,1,I+NESIM),JXOPE,NXOPE,WORK(KFPV),DFAC,
     &           IPQTOQ(1,IREP),IPRXRS)
           ENDIF
         END IF
C
C        e-p part of sigma vector
C        ------------------------
C
         IF (NXOPP .GT. 0) THEN
           IF (NASHT.NE.0) THEN
            CALL RORBGRD(.TRUE.,DV,FX(1,IM),FX(1,IM),FX(1,IM+1),
     &           FQX(1,I+NESIM),
     &           SPVEC(1,1,I+NESIM),JXOPP,NXOPP,WORK(KFPV),DFAC,
     &           IPQTOQ(1,IREP),IPRXRS)
           ELSE
            CALL RORBGRD(.TRUE.,DV,FX(1,IM),FX(1,IM),FX(1,IM+1),
     &           DUMMY,
     &           SPVEC(1,1,I+NESIM),JXOPP,NXOPP,WORK(KFPV),DFAC,
     &           IPQTOQ(1,IREP),IPRXRS)
           ENDIF
         ENDIF
         IF (NASHT .GT. 0) IM = IM + 1
C
      END DO
C
C
      CALL MEMREL('RORBGRD.KFPV',WORK,KFPV,KFPV,KFREE,LFREE)
C
      IF(IPRXRS.GE.5) THEN
        IF (NXOPE .GT. 0) THEN
         WRITE(LUPRI,'(/I5,A)') NOSIM,
     &        ' (e-e) orbital part of orbital (e-e) sigma vectors'
         CALL PRBVEC(LUPRI,SEVEC(1,1,1),NOSIM,NXOPE)
        ENDIF
        IF (NXOPP .GT. 0) THEN
         WRITE(LUPRI,'(/I5,A)') NOSIM,
     &        ' (e-p) orbital part of orbital (e-p) sigma vectors'
         CALL PRBVEC(LUPRI,SPVEC(1,1,1),NOSIM,NXOPP)
        ENDIF
      ENDIF
C
      CALL RTKTIME(.FALSE.,16)
      CALL QEXIT('RSIGOO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck rsigoogb */
      SUBROUTINE RSIGOOGB(SEVEC,SPVEC,GPOE,GPOP,BEVEC,BPVEC,
     &                    JXOPE,JXOPP,WORK,LWORK)
C***********************************************************************
C
C     For KR-MCSCF wave function optimizations:
C     Add the G * B term.
C
C     Written by J. Thyssen - Nov 30 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION SEVEC(*), SPVEC(*), GPOE(*), GPOP(*), BEVEC(*), BPVEC(*)
      DIMENSION JXOPE(*), JXOPP(*)
      DIMENSION WORK(*)
C
#include "memint.h"
C
      CALL QENTER('RSIGOOGB')
      CALL RTKTIME(.TRUE.,16)
C
C     Allocate memory for temporary matrix.
C     -------------------------------------
C
      CALL MEMGET('REAL',KUBO,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUGO,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUGB,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Reorder gradient back
C     ---------------------
C
C     Transform gradient from (NZXOPE,NZ) + (NZHOPE-NZXOPE,NZ) format to
C     (NZHOPE,NZ) format:
C
      CALL MEMGET('REAL',KGPOE,NZHOPE*NZ,WORK,KFREE,LFREE)
      IF (NZ .GT. 1 .AND. NZHOPE.GT.NZXOPE) THEN
         CALL INVREOGP(NZXOPE,NZHOPE,NZ,WORK(KGPOE),GPOE)
      ELSE
         CALL DCOPY(NZHOPE*NZ,GPOE,1,WORK(KGPOE),1)
      END IF
C
C     Add G * B term
C     --------------
C
      CALL RSIGOOGB1(SEVEC,SPVEC,WORK(KGPOE),GPOP,BEVEC,BPVEC,
     &        JXOPE,JXOPP,
     &        WORK(KUBO),WORK(KUGO),WORK(KUGB),
     &        WORK(KFREE),LFREE)
C
      CALL MEMREL('RSIGOOGB',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RSIGOOGB')
      CALL RTKTIME(.FALSE.,16)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck rsigoogb1 */
      SUBROUTINE RSIGOOGB1(SEVEC,SPVEC,GPOE,GPOP,BEVEC,BPVEC,
     &                     JXOPE,JXOPP, UBO,UGO,UGB,WORK,LWORK)
C***********************************************************************
C
C     For KR-MCSCF wave func optimizations:
C     Add the G * B term.
C
C     Written by J. Thyssen - Nov 30 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbxrs.h"
C
      DIMENSION SEVEC(NZXOPE,NZ,*), SPVEC(NZXOPP,NZ,*)
      DIMENSION GPOE(*),GPOP(*)
      DIMENSION BEVEC(NZXOPE,NZ,*),BPVEC(NZXOPP,NZ,*)
      DIMENSION JXOPE(2,*), JXOPP(2,*)
      DIMENSION UBO(N2ORBX,*),UGO(N2ORBX,*),UGB(NORBT,NORBT,NZ)
      DIMENSION WORK(*)
C
#include "memint.h"
C
      CALL QENTER('RSIGOOGB1')
C
C     *******************************
C     *** Calculate sigma vectors ***
C     *******************************
C
      IF (JSYMOP .NE. 1) THEN
         CALL QUIT('*** ERROR in RSIGOOGB1 : JSYMOP .ne. 1 ***')
      END IF
C
C     Scatter GPOE & GPOP into UGO
C     ----------------------------
C
      CALL DZERO(UGO,N2ORBXQ)
      CALL XOPSCT(GPOE,JXOPE,NZHOPE,UGO,NORBT,NZ)
      CALL XOPSCT(GPOP,JXOPP,NZXOPP,UGO,NORBT,NZ)
C
C     Construct full UGO
C     ------------------
C
      CALL FULMAT('A',NORBT,NORBT,UGO(1,1))
      DO IZ = 2, NZ
         CALL FULMAT('S',NORBT,NORBT,UGO(1,IZ))
      END DO
C
C     ***************************
C     *** e-e sigma vector(s) ***
C     ***************************
C
C
      DO I = 1, NESIM
C
C        Scatter BEVEC into UBO
C        ----------------------
C
         CALL DZERO(UBO,N2ORBXQ)
         CALL XOPSCT(BEVEC(1,1,I),JXOPE,NZXOPE,UBO,NORBT,NZ)
C
C        Construct full UBO
C        ------------------
C
         CALL FULMAT('A',NORBT,NORBT,UBO(1,1))
         DO IZ = 2, NZ
            CALL FULMAT('S',NORBT,NORBT,UBO(1,IZ))
         END DO
C
C        Calculate G * b^t
C        -----------------
C
         CALL QGEMM(NORBT,NORBT,NORBT,D1,
     &        'N','N',IPQTOQ(1,0),UGO,NORBT,NORBT,NZ,
     &        'T','N',IPQTOQ(1,0),UBO,NORBT,NORBT,NZ,
     &           D0,IPQTOQ(1,0),UGB,NORBT,NORBT,NZ)
C
         IF ( IPRXRS .GE. 30 ) THEN
            WRITE(LUPRI,'(1X,A)')
     &           '(RSIGOOGB1) G * b^t matrix (e-e trial vector)'
            CALL PRQMAT(UGO,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,0),LUPRI)
            CALL PRQMAT(UBO,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,0),LUPRI)
            CALL PRQMAT(UGB,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,0),LUPRI)
         END IF
         DO IOP = 1, NZXOPE
            K = JXOPE(1,IOP)
            L = JXOPE(2,IOP)
            SEVEC(IOP,1,I) = SEVEC(IOP,1,I) -
     &          DP5 * (UGB(K,L,1) - UGB(L,K,1) )
            DO IZ = 2, NZ
               SEVEC(IOP,IZ,I) = SEVEC(IOP,IZ,I) -
     &              DP5 * (UGB(K,L,IZ) - UGB(L,K,IZ))
            END DO
         END DO
         DO IOP = 1, NZXOPP
            K = JXOPP(1,IOP)
            L = JXOPP(2,IOP)
            SPVEC(IOP,1,I) = SPVEC(IOP,1,I) -
     &           DP5*(UGB(K,L,1) - UGB(L,K,1))
            DO IZ = 2, NZ
               SPVEC(IOP,IZ,I) = SPVEC(IOP,IZ,I) -
     &              DP5*(UGB(K,L,IZ) - UGB(L,K,IZ))
            END DO
         END DO
      END DO
C
C     ***************************
C     *** e-p sigma vector(s) ***
C     ***************************
C
C
      DO I = 1, NPSIM
C
C        Scatter BPVEC into UBO
C        ----------------------
C
         CALL DZERO(UBO,N2ORBXQ)
         CALL XOPSCT(BPVEC(1,1,I),JXOPP,NZXOPP,UBO,NORBT,NZ)
C
C        Construct full UBO
C        ------------------
C
         CALL FULMAT('A',NORBT,NORBT,UBO(1,1))
         DO IZ = 2, NZ
            CALL FULMAT('S',NORBT,NORBT,UBO(1,IZ))
         END DO
C
C        Calculate G * b^t
C        -----------------
C
         CALL QGEMM(NORBT,NORBT,NORBT,D1,
     &        'N','N',IPQTOQ(1,0),UGO,NORBT,NORBT,NZ,
     &        'T','N',IPQTOQ(1,0),UBO,NORBT,NORBT,NZ,
     &           D0,IPQTOQ(1,0),UGB,NORBT,NORBT,NZ)
C
         IF ( IPRXRS .GE. 30 ) THEN
            WRITE(LUPRI,'(1X,A)')
     &           '(RSIGOOGB1) G * b^t matrix (e-p trial vector)'
            CALL PRQMAT(UGO,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,0),LUPRI)
            CALL PRQMAT(UBO,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,0),LUPRI)
            CALL PRQMAT(UGB,NORBT,NORBT,NORBT,NORBT,NZ,
     &           IPQTOQ(1,0),LUPRI)
         END IF
C
         DO IOP = 1, NZXOPE
            K = JXOPE(1,IOP)
            L = JXOPE(2,IOP)
            SEVEC(IOP,1,I+NESIM) = SEVEC(IOP,1,I+NESIM) -
     &           DP5 * (UGB(K,L,1) - UGB(L,K,1))
            DO IZ = 2, NZ
               SEVEC(IOP,IZ,I+NESIM) = SEVEC(IOP,IZ,I+NESIM) -
     &              DP5*(UGB(K,L,IZ) - UGB(L,K,IZ))
            END DO
         END DO
         DO IOP = 1, NZXOPP
            K = JXOPP(1,IOP)
            L = JXOPP(2,IOP)
            SPVEC(IOP,1,I+NESIM) = SPVEC(IOP,1,I+NESIM) -
     &           DP5*(UGB(K,L,1) - UGB(L,K,1))
            DO IZ = 2, NZ
               SPVEC(IOP,IZ,I+NESIM) = SPVEC(IOP,IZ,I+NESIM) -
     &              DP5*(UGB(K,L,IZ) - UGB(L,K,IZ))
            END DO
         END DO
      END DO
C
C
      IF(IPRXRS.GE.30) THEN
         WRITE(LUPRI,'(/1X,A,I5,A)')
     &        '(RSIGOOGB1) ',NOSIM,' (e-e) orbital part(s) of orbital'//
     &        ' sigma vector(s) after G * B'
         CALL PRBVEC(LUPRI,SEVEC(1,1,1),NOSIM,NZXOPE)
         WRITE(LUPRI,'(/1X,A,I5,A)')
     &        '(RSIGOOGB1) ',NOSIM,' (e-p) orbital part(s) of orbital'//
     &        ' sigma vector(s) after G * B'
         CALL PRBVEC(LUPRI,SPVEC(1,1,1),NOSIM,NZXOPP)
      ENDIF
C
      CALL QEXIT('RSIGOOGB1')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsigco */
      SUBROUTINE RSIGCO(SCVEC,IBCVC,FX,H2ACX,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate orbital contribution to conf. sigma vector.
C     (Equation 3.30 in the KRMC paper (JCP 104, 4083 (1996)).
C
C     Input:
C       IBCVC: index of the CI trial vectors.
C       FX   : one index transformed Fock matrices
C              (we only need the FC matrices)
C       H2ACX: one-index transformed two-electron integrals
C
C     Output:
C       SCVEC: the CI sigma vectors.
C
C     Written by J. Thyssen - Dec 5 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION SCVEC(*)
      DIMENSION FX(*),H2ACX(*)
      DIMENSION IBCVC(*)
      DIMENSION WORK(*)
C
#include "dcbxrs.h"
#include "dcborb.h"
C
      CALL QENTER('RSIGCO')
      CALL RTKTIME(.TRUE.,17)
      KFRSAV = KFREE
C
      IF (XRS_CIPROGRAM .EQ. 'GASCIP' .OR.
     &    XRS_CIPROGRAM .EQ. 'LUCIAREL') THEN
         CALL MEMGET('REAL',KCREF,NZCONFQ,WORK,KFREE,LFREE)
!        CALL READAC(LUBCI,NZCONFQ,WORK(KCREF),1)
         CALL READ_DX(LUBCI,1,NZCONFQ,WORK(KCREF))
C        use routine with MOLFDIR format for integrals
         CALL RSIGCOM(WORK(KCREF),SCVEC,FX,H2ACX,IBCVC,
     &        WORK,KFREE,LFREE)
         CALL MEMREL('RSIGCO',WORK,1,KFRSAV,KFREE,LFREE)
      ELSE
         WRITE(LUPRI,'(/A/2A)')
     &        '*** ERROR in RSIGCO ***',
     &        'No sigma vector implemented for CI program ',
     &        XRS_CIPROGRAM
         CALL QUIT('*** ERROR in RSIGCO ***')
      END IF
      CALL RTKTIME(.FALSE.,17)
      CALL QEXIT ('RSIGCO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsigcom */
      SUBROUTINE RSIGCOM(CREF,SCVEC,FX,H2ACX,IBCVC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate orbital contribution to conf. sigma vector.
C     (Equation 3.30 in the KRMC paper (JCP 104, 4083 (1996)).
C
C     This routine converts integrals to MOLFDIR format before
C     calling the CI program.
C
C     Input:
C       IBCVC: index of the CI trial vectors.
C       FX   : one index transformed Fock matrices
C              (we only need the FC matrices)
C       H2ACX: one-index transformed two-electron integrals
C       CREF : reference CI vector.
C
C     Output:
C       SCVEC: the CI sigma vectors.
C
C     Written by J. Thyssen - Dec 05 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
C
      DIMENSION CREF(NZCONF,*), SCVEC(NZCONFQ,*)
      DIMENSION FX(N2ORBXQ,*), H2ACX(NASHT*NASHT*NNASHX*3*NZ,*)
      DIMENSION IBCVC(*), WORK(*)
C
C
#include "dcbxrs.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CALL QENTER('RSIGCOM')
      KFRSAV = KFREE
C
C     *******************************
C     *** Calculate sigma vectors ***
C     *******************************
C
      CALL MEMGET('REAL',KMOLFFC,(2*NASHT)*(2*NASHT)*MIN(NZ,2),
     $     WORK,KFREE,LFREE)
      IOFF1 = 1 + NCSIM
      DO I = 1, NOSIM
C
C        Get active-active part of FCX Fock matrices.
C
C        Remember that the ordering in the FX array is
C
C        FVT[1..NCSIM], {FCX, FVX}[1..NOSIM]
C
         CALL MEMGET('REAL',KFCACX,N2ASHXQ,WORK,KFREE,LFREE)
         CALL RGETAC(FX(1,IOFF1),WORK(KFCACX),IPRXRS)
         IOFF1 = IOFF1 + 2
C
C        Calculate CI sigma vector
C        -------------------------
C
C        Transform quarternion FCAC (Dirac MO-basis) to Molfdir MO-basis
C
         CALL QFC2MFC(WORK(KFCACX),WORK(KMOLFFC),JSYMOP,JTIMOP,IPRXRS)
C
         IF (XRS_CIPROGRAM .EQ. 'GASCIP') THEN
C
C           Transform integrals from Dirac (NZ,3) basis to Molfdir basis
C
C           debug
c           CALL CHKNZ3R(H2ACX(1,I))
C           end debug
C
            CALL MEMGET('REAL',KMUUUU,(2*NASHT)**4 * 2,
     &                  WORK,KFREE,LFREE)
            CALL DNZ32M(H2ACX(1,I),WORK(KMUUUU),IPRXRS)
C
         END IF
C
C        Calculate sigma vector
C
         IF (XRS_CIPROGRAM .EQ. 'GASCIP') THEN
            CALL GASCIP_SVC(.FALSE.,DUMMY,SCVEC(1,I),NZCONF,
     &         WORK(KZCONF),CREF,WORK(KMOLFFC),WORK(KMUUUU),IPRXRS)
         ELSE IF (XRS_CIPROGRAM .EQ. 'LUCIAREL') THEN
            CALL LUCI_SIGMA(.FALSE.,DUMMY,CREF,SCVEC(1,I),
     &                      WORK(KMOLFFC),H2ACX(1,I),WORK,KFREE,LFREE)
         ELSE
            CALL QUIT('Unknown CI PROGRAM requested!')
         END IF
C
C        FIXME: document this factor -2
C        (the 2 is easy (cref normalized to 0.5))
C
         CALL DSCAL(NZCONFQ,-D2,SCVEC(1,I),1)
         IF ( IPRXRS .GE. 30 ) THEN
            WRITE(LUPRI,'(/1X,A,I2,A)')
     &           '(RSIGCOM) conf. part of orbital sigma vector',I
            CALL RPRCI(SCVEC(1,I),NZCONF,NZ,
     &           XRS_CIPROGRAM,THRZER,LUPRI)
         END IF
      END DO
      CALL MEMREL('RSIGCOM',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT ('RSIGCOM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Fmoli3 */
      SUBROUTINE FMOLI3(ITYP,FMO,BMAT,FBUF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate one-index transformed Fock-matrices
C     for use in generation of sigma vectors
C
C
C     Bit 1 of ITYP = 0   e-e rotations in BMAT
C                     1   e-p rotations in BMAT
C     Bit 2 of ITYP = 0   inactive density matrix
C                     1   active density matrix
C
C     FMO_{si} = BMAT_{sj}FBUF_{ji} - FBUF^t_{sj}Z_{ji}
C
C     For time-antisymmetric operators FBUF^t is i-transform
C
C     Written by T.Saue Sep 17 1996
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0,DM1 = -1.0D0)
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxrs.h"
C
      LOGICAL EPMIX,LBIT
      DIMENSION FMO(NORBT,NORBT,NZ),FBUF(NORBT,NORBT,NZ),
     &          BMAT(NORBT,NORBT,NZ),WORK(*)
C
C
C     Initialization
C     ==============
C
      IF(.NOT. LBIT(ITYP,1))     THEN
         EPMIX = NZXOPP.GT.0
      ELSE
         EPMIX = NZXOPE.GT.0
      ENDIF
      IREP = JSYMOP-1
C
C     Generate one-index transformed Fock matrix
C     ==========================================
C
      DO 10 I1 = 1,NFSYM
         I2  = MOD(I1+JOPSY,2) + 1
C
C        Do the whole transformation.
C
         II = IORB(I2) + 1
         NI = NORB(I2)
C
         IS = IORB(I1) + 1
         NS = NORB(I1)
C
         CALL QGEMM(NS,NI,NI,D1,
     &        'N','N',IPQTOQ(1,IREP),
     &        BMAT(IS,II,1),NORBT,NORBT,NZ,
     &        'N','N',IPQTOQ(1,0),
     &        FBUF(II,II,1),NORBT,NORBT,NZ,
     &        D1,IPQTOQ(1,IREP),FMO(IS,II,1),NORBT,NORBT,NZ)
C
         CALL QGEMM(NS,NI,NS,DM1,
     $        'N','N',IPQTOQ(1,0),
     $        FBUF(IS,IS,1),NORBT,NORBT,NZ,
     $        'N','N',IPQTOQ(1,IREP),
     $        BMAT(IS,II,1),NORBT,NORBT,NZ,
     $        D1,IPQTOQ(1,IREP),FMO(IS,II,1),NORBT,NORBT,NZ)
C
C
C
#ifdef UNDEF
C
C     There is a problem with this code for active-active rotations.
C
C       Parameters for first index
C
        NI = NOCC(I2) - NFRO(I2)
        II = IORB(I2) + NPSH(I2) + NFRO(I2) + 1
C
C
C       Parameters for second index
C
        IF    (.NOT. LBIT(ITYP,1)) THEN
C
C         Electronic...
C
          NS = NESH(I1)
          IS = IORB(I1) + NPSH(I1) + 1
C
        ELSE
C
C         Positronic...
C
          NS = NPSH(I1)
          IS = IORB(I1) + 1
C
        ENDIF
C
c       NS = NORB(I1)
c       IS = IORB(I1) + 1
c       EPMIX = .FALSE.
        IF (NI.NE.0 .AND. NS.NE.0) THEN
C
C         First part: FMO_{si} = FMO_{si} + BMAT_{sj}FBUF_{ji}
C
          CALL QGEMM(NS,NI,NI,D1,
     &          'N','N',IPQTOQ(1,IREP),
     &          BMAT(IS,II,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,0),
     &          FBUF(II,II,1),NORBT,NORBT,NZ,
     &          D1,IPQTOQ(1,IREP),FMO(IS,II,1),NORBT,NORBT,NZ)
C
C
C         Second part: FMO_{si} = FMO_{si} - FBUF^t_{st}BMAT_{ti}
C
          CALL QGEMM(NS,NI,NS,DM1,
     &          'N','N',IPQTOQ(1,0),
     &          FBUF(IS,IS,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,IREP),
     &          BMAT(IS,II,1),NORBT,NORBT,NZ,
     &          D1,IPQTOQ(1,IREP),FMO(IS,II,1),NORBT,NORBT,NZ)
          IF(EPMIX) THEN
C
C           Parameters for secondary index(note the switch...)
C
            IF    (.NOT. LBIT(ITYP,1)) THEN
              NT = NPSH(I1)
              IT = IORB(I1) + 1
            ELSE
              NT = NESH(I1)
              IT = IORB(I1) + NPSH(I1) + 1
            ENDIF
            IF(NT.NE.0) CALL QGEMM(NT,NI,NS,DM1,
     &          'N','N',IPQTOQ(1,0),
     &          FBUF(IT,IS,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,IREP),
     &          BMAT(IS,II,1),NORBT,NORBT,NZ,
     &          D1,IPQTOQ(1,IREP),FMO(IT,II,1),NORBT,NORBT,NZ)
          ENDIF
        ENDIF
#endif
C
 10   CONTINUE
C
C     Print section
C
      IF(IPRXRS.GE.8) THEN
        CALL HEADER('FMOLI3: F-G-matrix',-1)
        CALL PRQMAT(FMO,NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck spfexc   */
      SUBROUTINE SPFEXC(BVEC,JXOP,NXOP,IBEIG,JZ,JREP)
C***********************************************************************
C
C     Analyze spinfree trial vector to identify spin symmetry
C     Only one symmetry and one component IZ should be non-zero.
C
C     Written by T. Saue 2009
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)    
C
#include "dgroup.h"
      DIMENSION BVEC(NXOP,NZ),JXOP(2,NXOP),IBEIG(*)
      LOGICAL JJ_JZ_SET
C
#include "ibtfun.h"
      JZ  = 0
      TMP = D0
      JJ_JZ_SET = .FALSE.
      DO IZ = 1,NZ
        DO J = 1,NXOP
          IF(ABS(BVEC(J,IZ)).GT.TMP) THEN
            TMP = ABS(BVEC(J,IZ))
            JZ = IZ
            JJ = J
            JJ_JZ_SET = .TRUE.
          ENDIF
        ENDDO
      ENDDO
      IF (.NOT.JJ_JZ_SET) THEN
       WRITE(LUPRI,*) 'SPFEXC: indexes JZ,JJ not set !'
       WRITE(LUPRI,*) 'NZ,NXOP: ',NZ,NXOP
       CALL QUIT('SPFEXC: indexes JZ,JJ not set !')
      ENDIF
      I1 = JXOP(1,JJ)
      I2 = JXOP(2,JJ)
      JREP = IBTXOR(IBEIG(I1),IBEIG(I2))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrsappfq */
      SUBROUTINE XRSAPPFQ(FQX,FX,DV,WORK,LWORK)
C***********************************************************************
C
C     Approximate positronic part of FQ.
C
C     Written by J. Thyssen - Feb 18 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbxrs.h"
#include "dummy.h"
C
      DIMENSION FQX(NORBT*NASHT*NZ,*), FX(N2ORBXQ,*), DV(*)
      DIMENSION WORK(*)
C
      CALL QENTER('XRSAPPFQ')
C
      IF (XRS_NOFQX) THEN
C
C        Approximate the whole FQX with DV * FVX
C
         DO I = 1, NOSIM
            CALL APPFQ('G',FQX(1,I),DV,FX(1,2*I),WORK,LWORK)
         END DO
C
C
      ELSE IF (XRS_NOPFQ) THEN
C
C        Approximate positronic part of FQX matrices.
C
         IF (IPRXRS .GT. 1) WRITE(LUPRI,'(/A)')
     &      ' approximating positronic part of FQX ...'
         DO I = 1, NOSIM
            CALL APPFQ('P',FQX(1,I),DV,FX(1,2*I),WORK,LWORK)
         END DO
C
      END IF
C
      CALL QEXIT('XRSAPPFQ')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ACCS2XX */
      SUBROUTINE ACCS2XX(TPVEC,TMVEC,BBUF,FREQ,NVARQ)
C***********************************************************************
C
C     Constructing contribution from S[2] to residuals,
C     see Eq. (21) in J. Chem. Phys. 133 (2010) 064105
C
C     On input
C       TPVEC contains T(+)
C       TMVEC contains T(-)
C
C     T. Saue and S. Villaume Dec 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,DM1=-1.0D0)
C
#include "dcbxrs.h"
      DIMENSION BBUF(NVARQ,NDAMP),TPVEC(NVARQ,NDAMP),TMVEC(NVARQ,NDAMP)
      FFAC = D1
      IF(IMFREQ) FFAC = DM1
      IF(DAMPFREQ.NE.D0) THEN
C  Take care of Complex case
        CALL DCOPY(NVARQ*NDAMP,          TMVEC(1,1),1,BBUF,1)
        CALL DSCAL(NVARQ*NDAMP,    -FREQ,TMVEC(1,1),1)
        CALL DAXPY(NVARQ      , DAMPFREQ, BBUF(1,2),1,TMVEC(1,1),1)
        CALL DAXPY(NVARQ      ,-DAMPFREQ, BBUF(1,1),1,TMVEC(1,2),1)
CPP: error in the first line of eq (21) of JCP133,064105,2010
        CALL DCOPY(NVARQ*NDAMP,          TPVEC(1,1),1,BBUF,1)
        CALL DSCAL(NVARQ*NDAMP,    -FREQ,TPVEC(1,1),1)
        CALL DAXPY(NVARQ      , DAMPFREQ, BBUF(1,2),1,TPVEC(1,1),1)
        CALL DAXPY(NVARQ      ,-DAMPFREQ, BBUF(1,1),1,TPVEC(1,2),1)
      ELSE
C  Take care of Real or Imaginary case
        CALL DSCAL(NVARQ,-FREQ,TMVEC(1,1),1)
        CALL DSCAL(NVARQ,-FFAC*FREQ,TPVEC(1,1),1)
      ENDIF
C
      RETURN
      END

