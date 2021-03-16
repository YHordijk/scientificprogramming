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

C FILE    : pamgrd.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamgrd */
      SUBROUTINE PAMGRD(INDXPR,GPCI,GPOE,GPOP,
     &                  JOPSY,JXOPE,JXOPP,NCONF,NXOPE,NXOPP,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct gradient for property of index INDXPR in /XCBPRP/
C     GPC - configurational part
C     GPO - orbital part
C
C     Written by T.Saue - July 1996
C     Last revision: Jul 10 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxpr.h"
C
      DIMENSION GPCI(*),GPOE(*),GPOP(*),JXOPE(*),JXOPP(*),
     &          WORK(*)
      CALL QENTER('PAMGRD')
C
C     Consistency check
C
      ISYM  = IPRPSYM(INDXPR)
      IOPSY = JBTOF(ISYM-1,1)
      IF(IOPSY.NE.JOPSY) THEN
        WRITE(LUPRI,'(A)') 
     &    '*** ERROR(PAMGRD): Incorrect fermion ircop of property ',
     &    PRPNAM(INDXPR)
        CALL QUIT('*** ERROR(PAMGRD) Incorrect fermion ircop')
      ENDIF
C
C     Configurational part
C     ====================
C
      IF(NCONF.GT.0) THEN
C        CALL GRDCI(INDXPR,GPC,WORK,KFREE,LFREE,IPRINT)
      ENDIF      
C
C     Orbital part
C     ============
C
      NXOPT = NXOPE + NXOPP
      IF(NXOPT.GT.0) THEN
        CALL GRDORB(INDXPR,GPOE,GPOP,JXOPE,JXOPP,NXOPE,NXOPP,IOPSY,
     &              WORK,KFREE,LFREE,IPRINT)
      ENDIF
C
      CALL QEXIT('PAMGRD')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grdorb */
      SUBROUTINE GRDORB(INDXPR,GPOE,GPOP,JXOPE,JXOPP,NXOPE,NXOPP,IOPSY,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct orbital part of gradient for property INDXPR 
C     listed in /XCBXPR/
C
C     The property gradient is defined by
C
C              | g  | 
C     E[1]_X = | g* | ; 
C
C     where
C
C     g_{pq} = <0|[-q(+)p,H_X]|0> = -h_{X;ps}D_{qs} + h_{X;rq}D_{rp}
C
C     Active-inactive block:
C
C     g_{ai} = -h_{X;ai}
C
C     Orbital classes:
C
C       Inactive orbitals  :  i,j,k,l
C       Active orbitals    :  v,w,x,y,z
C       Secondary orbitals :  a,b,c,d,e
C       General orbitals   :  m,n,o,p,q,r,s,t,u,v
C
C     GRDORB constructs g_{pq}.
C     It employs the orbital rotation index JXOP(p,q)
C     with the following ordering:
C
C          IA : inactive-active   
C          IS : inactive-secondary
C          AS : active-secondary
C          AB : active-active
C
C     Indices are therefore permuted...
C
C     Written by T.Saue - July 1996
C     Last revision: Jul 10 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
C
      DIMENSION GPOE(*),GPOP(*),JXOPE(*),JXOPP(*),WORK(*)
C
      CALL QENTER('GRDORB')
      KFRSAV = KFREE
C
      CALL MEMGET2('REAL','PRPMO',KPRPMO,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','CMO'  ,KCMO  ,NCMOTQ ,WORK,KFREE,LFREE)
      IF(SPINFR) THEN
        CALL MEMGET2('INTE','IBTYP',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
      ELSE
        KIBRP = KFREE
      ENDIF
C
      CALL GRDOR1(INDXPR,GPOE,GPOP,JXOPE,JXOPP,NXOPE,NXOPP,IOPSY,
     &            WORK(KPRPMO),WORK(KCMO),WORK(KIBRP),
     &            WORK,KFREE,LFREE,IPRINT) 
C
      CALL MEMREL('GRDORB',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('GRDORB')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grdor1 */
      SUBROUTINE GRDOR1(INDXPR,GPOE,GPOP,JXOPE,JXOPP,NXOPE,NXOPP,IOPSY,
     &                  PRPMO,CMO,IBEIG,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct orbital part of gradient for property INDXPR 
C     listed in /XCBXPR/
C
C     The property gradient is defined by
C
C              | g  | 
C     E[1]_X = | g* | ; 
C
C     where
C
C     g_{pq} = <0|[-q(+)p,H_X]|0> = -h_{X;ps}D_{qs} + h_{X;rq}D_{rp}
C
C     Active-inactive block:
C
C     g_{ai} = -h_{X;ai}
C
C     Orbital classes:
C
C       Inactive orbitals  :  i,j,k,l
C       Active orbitals    :  v,w,x,y,z
C       Secondary orbitals :  a,b,c,d,e
C       General orbitals   :  m,n,o,p,q,r,s,t,u,v
C
C     GRDORB constructs g_{pq}.
C     It employs the orbital rotation index JXOP(p,q)
C     with the following ordering:
C
C          IA : inactive-active   
C          IS : inactive-secondary
C          AS : active-secondary
C          AB : active-active
C
C     Indices are therefore permuted...
C
C     Written by T.Saue July 1996
C     Last revision: Feb 17 2000 (hjaaj, new parameter list for PRPMAT)
C
C     MI,HJAaJ(&TEC), march 2003: added new routine calling
C     that gets the total property gradient dH(B)/dB|B=0 when using
C     London atomic orbitals
C
C     S. Knecht, Jan 2010: added read-in of coefficients from a 
C                          previous MCSCF run.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbham.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbpsi.h"
C
      DIMENSION GPOE(*),GPOP(*),JXOPE(*),JXOPP(*),
     &          PRPMO(NORBT,NORBT,NZ),CMO(*),IBEIG(*),WORK(*)
      LOGICAL TOBEK
C
      CALL QENTER('GRDOR1')
C
      TOBEK =.FALSE.
C
C     Read coefficients
C     =================
C
C
C     check for coefficients from the current KRMCSCF run
      INQUIRE(FILE='KRMCSCF',EXIST=TOBEK)
C
      IF(TOBEK.and.DOKRMC)THEN
C
        IF(SPINFR) THEN
          WRITE(LUPRI,'(A)') ' *** ERROR in GRDOR1: reading'//
     &                       ' of MCSCF coefficients in a spinfree'//
     &                       ' calculation. ***'
          CALL QUIT('Spinfree MCSCF calculations are not possible.')
        END IF
        CALL OPNFIL(LUKRMC,'KRMCSCF','UNKNOWN','GRDOR1')
        JRDMO = -1
        CALL RREADMO(CMO,JRDMO,1,LUKRMC)
C       CLOSE (LUKRM5,STATUS='KEEP')
C       ... set to zero
        TOTERG = 0.0D0
      ELSE
C
C       read HF coefficients 
        CALL OPNFIL(LUCOEF,'DFCOEF','OLD','GRDOR1')
        IF(SPINFR) THEN
          CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,IBEIG,TOTERG,10)
        ELSE
          CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,DUM,TOTERG,2)
        ENDIF
        CLOSE(LUCOEF,STATUS='KEEP')
      END IF
C
C     Get property matrix in MO-basis and scale by factor -1
C     ======================================================
C
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
CMI... distinguish between standard and LAO-based atomic orbitals
CMI NOTE/27.10'04 The new IBEIG parameters for spin-free properties !
CMI This has to be adapted also for LAO properties !
      IF (PRPNAM(INDXPR)(1:3).EQ.'LAO') THEN
        IB = 7 ! 0111 = full H(B) prop. gradient
!gosia: fac2lo is not set? i'm setting it to 1.0
        fac2lo = 1.0d0
        CALL GETDHB(IB,FAC2LO,PRPMO,INDXPR,IOPSY,
     &                    WORK,KFREE,LFREE,IPRINT)
      ELSE
        CALL PRPMAT(INDXPR,IOPSY,PRPMO,.TRUE.,WORK,CMO,IBEIG,ICMOQ,NORB,
     &              WORK,KFREE,LFREE,IPRINT)
      ENDIF

      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Hermicity of operator
C     =====================
C
      ITIM = IPRPTIM(INDXPR)
      IREP = IPRPSYM(INDXPR) - 1
      IQ   = IPQTOQ(1,IREP)
      IH   = IHQMAT(IQ,ITIM)
C
C
C
C     Construct gradient - electronic part
C     ====================================
C
      KFRSAV = KFREE
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDV,N2ASHXQ,WORK,KFREE,LFREE)
         CALL GETMAT('DVMO',WORK(KDV),IPRINT,TOBEK,WORK(KFREE),LFREE)
      ELSE
         CALL MEMGET('REAL',KDV,0,WORK,KFREE,LFREE)
      END IF
      IF(NXOPE.GT.0) THEN
        CALL DZERO(GPOE,NZ*NXOPE)
        CALL PRPORB(GPOE,JXOPE,NXOPE,PRPMO,NORBT,WORK(KDV),NZ,IH,IPRINT)
c       no scaling
        CALL DSCAL(NZ*NXOPE,DM1,GPOE,1)
      ENDIF
C
C     Construct gradient - positronic part
C     ====================================
C
      IF(NXOPP.GT.0) THEN
        CALL DZERO(GPOP,NZ*NXOPP)
        CALL PRPORB(GPOP,JXOPP,NXOPP,PRPMO,NORBT,WORK(KDV),NZ,IH,IPRINT)
c       no scaling
        CALL DSCAL(NZ*NXOPP,DM1,GPOP,1)
      ENDIF
      CALL MEMREL('GRDOR1',WORK,1,KFRSAV,KFREE,LFREE)
C
C     Print section
C     =============
C
#ifndef ANALYZE_PROPERTY_GRADIENT
      IF(IPRINT.GE.4) THEN
#endif
        CALL HEADER('Output from GRDORB',-1)
        IREP = IPRPSYM(INDXPR)-1
        IF(NXOPE.GT.0) THEN
          WRITE(LUPRI,'(A,A16)') 
     &    '* Gradient (e-e) of property ',PRPNAM(INDXPR)
          CALL PRQMAT(GPOE,NXOPE,1,NXOPE,1,NZ,IPQTOQ(1,IREP),
     &                LUPRI)
        ENDIF
        IF(NXOPP.GT.0) THEN
          WRITE(LUPRI,'(A,A16)') 
     &    '* Gradient (e-p) of property ',PRPNAM(INDXPR)
          CALL PRQMAT(GPOP,NXOPP,1,NXOPP,1,NZ,IPQTOQ(1,IREP),
     &                LUPRI)
        ENDIF
#ifndef ANALYZE_PROPERTY_GRADIENT
      ENDIF
#endif
      CALL QEXIT ('GRDOR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xopgat */
      SUBROUTINE XOPGAT(QVEC,JXOP,NXOP,QMAT,NDIM,NZ)
C***********************************************************************
C
C     Gather matrix elements into vector using JXOP
C
C     Written by T.Saue August 2 1996
C     Last revision Aug 2 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION QVEC(NXOP,NZ),JXOP(2,NXOP),QMAT(NDIM,NDIM,NZ)
C
      DO IZ = 1,NZ
#include "ivdep.h"
        DO IXOP = 1,NXOP
          IA = JXOP(1,IXOP)
          IM = JXOP(2,IXOP)
          QVEC(IXOP,IZ) = QMAT(IM,IA,IZ)
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xopsct */
      SUBROUTINE XOPSCT(QVEC,JXOP,NXOP,QMAT,NDIM,NZ)
C***********************************************************************
C
C     Scatter vector elements into matrix using JXOP
C
C     Written by T.Saue August 2 1996
C     Last revision Aug 2 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION QVEC(NXOP,NZ),JXOP(2,NXOP),QMAT(NDIM,NDIM,NZ)
C
      DO IZ = 1,NZ
#include "ivdep.h"
        DO IXOP = 1,NXOP
          INDI = JXOP(1,IXOP)
          INDA = JXOP(2,IXOP)
          QMAT(INDA,INDI,IZ) = QVEC(IXOP,IZ)
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpmat */
      SUBROUTINE PRPMAT(INDXPR,IOPSY,PRPMO,ONFILE,PRPINT,CMO,
     &                  IBEIG,IPCMO,NPMO,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct MO-matrix for property INDXPR listed in /XCBXPR/
C
C     -- NPMO(i) is dimension of output property matrix 
C                in fermion symmetry i
C     -- IPCMO(i) is off-set to first MO-coeff in fermion symmetry i
C        NPMO = NORB and IPCMO = ICMOQ gives full PRPMO matrix
C        NPMO = NASH and IPCMO = xxxxx gives      PRPMO matrix over 
C                                active orbitals
C
C   On input:      INDXPR, IOPSY
C ============     CMO - molecular orbitals ..
C                  IBEIG - boson info on orbitals...
C                  IPCMO, NPMO
C               
C
C  On output:  PRPMO(N2ORBXQ) - resulting operator in MO basis
C ===========  (taken care of the spin-sfree case!)
C
C
C     Written by T.Saue July 1996
C     Revisions:
C       Feb 17 2000 (hjaaj, IPCMO(2), NPMO(2) parameters for 
C                    partial MO transf.)
C       Dec 2005, MI - extended for SQUARE property matrixes
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxpr.h"
C
      LOGICAL ONFILE
      DIMENSION PRPMO(*),PRPINT(*),CMO(*),IBEIG(*),IPCMO(*),NPMO(*),
     &          WORK(*)
C
      CALL QENTER('PRPMAT')
      KFRSAV = KFREE
C     Determine whether integrals are fetched on triangular or square form;
C     check if first scalar property operator is has diagonal symmetry or not
      NPMOT = NPMO(1)
      IF (NFSYM.EQ.2) NPMOT = NPMOT + NPMO(2)
      N2PMOX = NPMOT*NPMOT
      IF (IPRINT .GE. 6) THEN
        CALL HEADER('Output from PRPMAT',-1)
        WRITE(LUPRI,'(/A,4I10)')
     &  ' Input to PRPMAT: IPCMO, NPMO:',
     &  (IPCMO(I),NPMO(I),I = 1,NFSYM)
        WRITE(LUPRI,*)' NPMOT,N2PMOX:',NPMOT,N2PMOX
      ENDIF
C
C     Memory allocation
C
      CALL MEMGET('REAL',KPRPAO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ      ,WORK,KFREE,LFREE)
      CALL PRPMA1(INDXPR,PRPMO,ONFILE,PRPINT,CMO,IBEIG,IPCMO,NPMO,
     &            N2PMOX,WORK(KPRPAO),WORK(KFIRST),
     &            WORK,KFREE,LFREE,IPRINT)
      CALL MEMREL('PRPMAT',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('PRPMAT')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpma1 */
      SUBROUTINE PRPMA1(INDXPR,PRPMO,ONFILE,PRPINT,CMO,IBEIG,IPCMO,NPMO,
     &                  N2PMOX,PRPAO,FIRST,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct MO-matrix for property INDXPR listed in /XCBXPR/
C
C     Written by T.Saue July 1996
C     Last revision: Feb 17 2000 (hjaaj, see PRPMAT)
C     Dec 2005, MI - adapting for both triangular and square operators
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dummy.h"
C
      LOGICAL QBTR,TOBE,FIRST(NZ),ONFILE
      DIMENSION PRPMO(N2PMOX,NZ),CMO(*),IBEIG(*),IPCMO(2),NPMO(2),
     &          PRPAO(N2BBASX,NZ),PRPINT(*),WORK(*)
      DIMENSION J2PMOX(2,2)
C
      CALL QENTER('PRPMA1')
C
C     Set up index information for
C     PRPMO(N2PMOX,NZ) [equivalent to "PRPMO(NPMOT,NPMOT,NZ)"]
C     such that PRPMO(J2PMOX(1,1),NZ) gives "PRPMO(g1,g1,NZ)"
C               PRPMO(J2PMOX(2,1),NZ) gives "PRPMO(u1,g1,NZ)"
C               PRPMO(J2PMOX(1,2),NZ) gives "PRPMO(g1,u1,NZ)"
C               PRPMO(J2PMOX(2,2),NZ) gives "PRPMO(u1,u1,NZ)"
C     where g1 is index of first gerade orbital (always g1 = 1)
C       and u1 is index of first ungerade orbital (u1 = g1 + NPMO(1)).
C
      NPMOT = NPMO(1)
      IF (NFSYM.EQ.2) NPMOT = NPMOT + NPMO(2)
      JOFF = 1
      DO J = 1,NFSYM
        IJOFF = JOFF
        DO I = 1,NFSYM
          J2PMOX(I,J) = IJOFF
          IJOFF = IJOFF + NPMO(I)
        END DO
        JOFF = JOFF + NPMO(J)*NPMOT
      END DO
C
C     Extract symmetry information - ISYM, IREP, ITIM
C     =============================
C
      ISYM = IPRPSYM(INDXPR)
      IREP = ISYM - 1
      ITIM = IPRPTIM(INDXPR)
C
C     Get property matrix in AO-basis
C     ===============================
      CALL PRPMSAO(INDXPR,PRPAO,ONFILE,PRPINT,FIRST,
     &             WORK(KFREE),LFREE,IPRINT)
C .. control print out of 2c/4c-component
C    operator in the SA-AO basis
      IF (IPRINT.GE.15) THEN
        WRITE(LUPRI,'(/A,A16,1X,A,I1,A,I2)')
     &'PRPMA1: Total 4c or 2c matrix (in SA-AO basis) of property ',
     & PRPNAM(INDXPR),'IREP=',IREP,' ITIM=',ITIM
        CALL PRQMAT(PRPAO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &            IPQTOQ(1,IREP),LUPRI)
      ENDIF
C
C     Transform to MO-basis only the blocks of correct fermion symmetry
C     =================================================================
C         
CTROND
!      IF (PRPNAM(INDXPR)(1:5).EQ.'BEDcZ'.AND.
!     &    PRPNAM(INDXPR)(8:10).EQ.'X00' .AND.
!     &    PRPNAM(INDXPR)(14:16).EQ.'Z00') THEN
!          DO I = 1,NXPAR
!             WRITE(6,*) 'JXOP..',JXOP(1,I),JXOP(2,I)
!          ENDDO
!          I1=1
!          I2=2
CTROND    Radium
!         II = 165
!         IA = 213
CTROND     Magnesium
C          II = 30
C          IA = 40
!          CALL GPOTEST(PRPNAM(INDXPR),II,IA,NZ,
!     &           PRPAO(I2BASX(I1,I2)+1,1),NTBAS(0),
!     &           CMO(IPCMO(I1)+1),NFBAS(I1,0),NORB(I1),
!     &           CMO(IPCMO(I2)+1),NFBAS(I2,0),NORB(I2))
!      ENDIF
      CALL DZERO(PRPMO,N2PMOX*NZ)
      IOPSY = JBTOF(IREP,1)
!     write(*,*) 'N2PMOX,NZ:', N2PMOX,NZ
!     write(*,*) 'IOPSY:', IOPSY
!     write(*,*) 'IREP:', IREP
      QBTR = .FALSE.
      DO I1 = 1,NFSYM
        I2  = MOD(I1+IOPSY,2) + 1
        IF (NFBAS(I1,0).GT.0.AND.NFBAS(I2,0).GT.0) THEN
          IF (IPRINT.GE.20) THEN
            CALL HEADER('PRPMA1: entering symm.blocked AO-matrix',-1)
            WRITE(LUPRI,'(A,I1,1X,I1,A,I1)')
     &      '* Fermion irps no.',I1,I2,' of ',NFSYM
            CALL PRQMAT(PRPAO(I2BASX(I1,I2)+1,1),
     &                  NFBAS(I1,0),NFBAS(I2,0),NTBAS(0),NTBAS(0),NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
            IF (IPRINT.GE.25) THEN
              CALL HEADER(
     &        'PRPMA1: entering symm.blocked MO-coefficients',-1)
              WRITE(LUPRI,'(A,I1,1X,I1,A,I1)')
     &        '* Fermion ircop no.',I1,I2,' of ',NFSYM
              CALL PRQMAT(CMO(IPCMO(I1)+1),NFBAS(I1,0),NORB(I1),
     &         NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,IREP),LUPRI)
              CALL FLSHFO(LUPRI)
            ENDIF
          ENDIF

          IF(SPINFR.AND.((ITIM.EQ.1).OR.NOSPIN)) THEN
            CALL QBTRANS(IREP,'AOMO','S',D0,
     &           NFBAS(I1,0),NFBAS(I2,0),NPMO(I1),NPMO(I2),
     &           PRPAO(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),NZ,
     &           IPQTOQ(1,IREP),
     &           PRPMO(J2PMOX(I1,I2),1),NPMOT,NPMOT,NZ,IPQTOQ(1,IREP),
     &           CMO(IPCMO(I1)+1),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &           IBEIG(IORB(I1)+1),
     &           CMO(IPCMO(I2)+1),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &           IBEIG(IORB(I2)+1),
     &           WORK(KFREE),LFREE,IPRINT)
           QBTR = .TRUE.
           ELSE
             CALL QTRANS('AOMO','S',D0,
     &           NFBAS(I1,0),NFBAS(I2,0),NPMO(I1),NPMO(I2),
     &           PRPAO(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),NZ,
     &           IPQTOQ(1,IREP),
     &           PRPMO(J2PMOX(I1,I2),1),NPMOT,NPMOT,NZ,IPQTOQ(1,IREP),
     &           CMO(IPCMO(I1)+1),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &           CMO(IPCMO(I2)+1),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
           ENDIF ! IF(SPINFR.AND.((ITIM.EQ.1).OR.NOSPIN)) THEN
         ENDIF ! IF (NFBAS(I1,0).GT.0.AND.NFBAS(I2,0).GT.0) THEN
  
        IF (IPRINT.GE.15) THEN
          CALL HEADER('PRPMA1: Symm.blocked MO-matrix',-1)
          WRITE(LUPRI,'(A,I1,1X,I1,A,I1)')
     &    '* Fermion irps no.',I1,I2,' of ',NFSYM
          CALL PRQMAT(PRPMO(J2PMOX(I1,I2),1),NPMO(I1),NPMO(I2),
     &           NPMOT,NPMOT,NZ,IPQTOQ(1,IREP),LUPRI)
          CALL FLSHFO(LUPRI)
        ENDIF
!
      ENDDO
C
C     Print section
C     =============
C
      IF(IPRINT.GE.10) THEN
        CALL HEADER('PRPMA1: Total MO matrix of property',-1)
        IF (QBTR) THEN
         WRITE(LUPRI,'(4X,A)')
     &   'Applied BOSON symmetry based AO_to_MO transformation !'
        ELSE
         WRITE(LUPRI,'(4X,A)')
     &   'Applied FERMION symmetry based AO_to_MO transformation !'
        ENDIF
        WRITE(LUPRI,'(A,A16)') 
     &  '*MO-matrix of property ',PRPNAM(INDXPR)
        WRITE(LUPRI,'(2X,A,I3,A,I2,A,I2,A,I2)')
     &  'with parameters INDXPR=',INDXPR,' ISYM=',ISYM,
     &  ' IREP=',IREP,' ITIM=',ITIM
        CALL PRQMAT(PRPMO,NPMOT,NPMOT,NPMOT,NPMOT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
C
      CALL QEXIT('PRPMA1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gpget */
      SUBROUTINE GPGET(INDXPR,GP,JXPAR,NXPAR,CMO,IBEIG,
     &                 ITYP,ONFILE,PRPINT,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct orbital part of gradient for property INDXPR 
C     listed in /XCBXPR/
C
C     The property gradient is defined by
C
C
C              | (q+|Q) |   | -<0|[q ,Q]|Q> |   | e[1]_Q* |
C     E[1]_Q = | (q |Q) | = |  <0|[q+,Q]|Q> | = | e[1]_Q  |
C
C     In RPA q = a(+)m, where a is occupied and m secondary.
C     This routine will construct e[1]_Q*
C
C     e[1]_Q*(ma) = -<0|[a(+)m,Q]|Q> 
C                 = -<0|a(+)[m,Q]|0> - <0|[a(+),Q]m|0>
C                 = -<0|a(+)[m,Q]|0> - <0|,m(+)[a,Q]m|0>(*) = -Q_(ma)
C
C     GRDORB employes the orbital rotation index JXOP(a,m).
C     Indices are therefore permuted...
C
C     Written by T.Saue - July 1996
C     Last revision: Jul 10 1996 - tsaue
C
C  MI,HJ (march 2003) - added routine calling for cumulationg the dH/dB prop.gradient with LAO
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
       PARAMETER(DM1 = -1.0D0)
C
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbibn.h"
       LOGICAL ONFILE
       DIMENSION GP(*),CMO(*),IBEIG(*),JXPAR(*),PRPINT(*),WORK(*)
C
      CALL QENTER('GPGET')
      KFRSAV = KFREE
C
      ISYM  = IPRPSYM(INDXPR)
      IOPSY = JBTOF(ISYM-1,1)
      IF(ITYP.EQ.JBENDX) THEN
        IF(NXPAR.GE.0) THEN
          CALL MEMGET('REAL',KPRPMO,N2BBASXQ,WORK,KFREE,LFREE)
          CALL GPOGET(INDXPR,GP,JXPAR,NXPAR,WORK(KPRPMO),
     &                CMO,IBEIG,ONFILE,PRPINT,WORK,KFREE,LFREE,IPRINT)
          CALL MEMREL('GPGET.e-e',WORK,1,KFRSAV,KFREE,LFREE)
#ifndef ANALYZE_PROPERTY_GRADIENT
          IF(IPRINT.GE.4) THEN
#endif
            CALL HEADER('Output from GPGET',-1)
            IREP = ISYM - 1
            WRITE(LUPRI,'(A,A16)') 
     &      '* Gradient (e-e) of property ',PRPNAM(INDXPR)
            CALL PRQMAT(GP,NXPAR,1,NXPAR,1,NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
#ifndef ANALYZE_PROPERTY_GRADIENT
          ENDIF
#endif
        ENDIF
      ELSEIF(ITYP.EQ.JBPNDX) THEN
        IF(NXPAR.GE.0) THEN
          CALL MEMGET('REAL',KPRPMO,N2BBASXQ,WORK,KFREE,LFREE)
          CALL GPOGET(INDXPR,GP,JXPAR,NXPAR,WORK(KPRPMO),
     &                CMO,IBEIG,ONFILE,PRPINT,WORK,KFREE,LFREE,IPRINT)
          CALL MEMREL('GPGET.e-p',WORK,1,KFRSAV,KFREE,LFREE)
          IF(IPRINT.GE.4) THEN
            CALL HEADER('Output from GPGET',-1)
            IREP = ISYM-1
            WRITE(LUPRI,'(A,A16)') 
     &      '* Gradient (e-p) of property ',PRPNAM(INDXPR)
              CALL PRQMAT(GP,NXPAR,1,NXPAR,1,NZ,
     &                    IPQTOQ(1,IREP),LUPRI)
          ENDIF
        ENDIF
      ELSEIF(ITYP.EQ.JBCNDX) THEN
        WRITE(LUPRI,'(A)') 
     &    'GPGET ERROR: Configurational part not yet written !'
        CALL QUIT('GPGET: Configurational part not yet written !')
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'GPGET ERROR: ',
     &    'unknown ITYP = ',ITYP
        CALL QUIT('GPGET: Unknown ITYP !')
      ENDIF
C
C
      CALL QEXIT('GPGET')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck GPOGET */
      SUBROUTINE GPOGET(INDXPR,GPO,JXOP,NXOP,PRPMO,CMO,
     &                  IBEIG,ONFILE,PRPINT,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct orbital part of gradient for property INDXPR 
C     listed in /XCBXPR/
C
C     The property gradient is defined by
C
C              | (q+|Q) |   | -<0|[q ,Q]|Q> |   | e[1]_Q* |
C     E[1]_Q = | (q |Q) | = |  <0|[q+,Q]|Q> | = | e[1]_Q  |
C
C     In RPA q = a(+)m, where a is occupied and m secondary.
C     This routine will construct e[1]_Q*
C
C     e[1]_Q*(ma) = -<0|[a(+)m,Q]|Q> 
C                 = -<0|a(+)[m,Q]|0> - <0|[a(+),Q]m|0>
C                 = -<0|a(+)[m,Q]|0> - <0|,m(+)[a,Q]m|0>(*) = -Q_(ma)
C
C     GRDORB employes the orbital rotation index JXOP(a,m).
C     Indices are therefore permuted...
C
C     Written by T.Saue July 1996
C     Last revision: Feb 17 2000 hjaaj (new parameter list for PRPMAT)
C
C     MI/2003,2004- add dH/dB prop.gradient accumulation when using LAO
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "mxcent.h"
#include "dcbprp.h"
C
      LOGICAL ONFILE
      DIMENSION GPO(*),PRPMO(NORBT,NORBT,NZ),CMO(*),
     &          IBEIG(*),JXOP(NXOP),PRPINT(*),WORK(*)
C
      CALL QENTER('GPOGET')
C
C     Get property matrix and scale by factor -1
C     ==========================================
C
      ISYM  = IPRPSYM(INDXPR)
      IOPSY = JBTOF(ISYM-1,1)
      IF(ONFILE) THEN
        OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &           FILE='AOPROPER')
      ENDIF
      IF (PRPNAM(INDXPR)(1:3).EQ.'LAO') THEN
        IB = 7 ! 111
!gosia: fac2lo is not set? i'm setting it to 1.0
        fac2lo = 1.0d0
        CALL GETDHB(IB,FAC2LO,PRPMO,INDXPR,IOPSY,
     &                    WORK,KFREE,LFREE,IPRINT)
      ELSE
         CALL PRPMAT(INDXPR,IOPSY,PRPMO,ONFILE,PRPINT,
     &               CMO,IBEIG,ICMOQ,NORB,
     &               WORK,KFREE,LFREE,IPRINT)
      ENDIF

      IF(ONFILE) CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Hermicity of operator
C     =====================
C
      ITIM = IPRPTIM(INDXPR)
      IREP = IPRPSYM(INDXPR) - 1
      IQ   = IPQTOQ(1,IREP)
      IH   = IHQMAT(IQ,ITIM)
C
C     Construct gradient - electronic part
C     ====================================
C
      KFRSAV = KFREE
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDV,N2ASHXQ,WORK,KFREE,LFREE)
         CALL GETMAT('DVMO',WORK(KDV),IPRINT,.FALSE.,WORK(KFREE),LFREE)
      ELSE
         CALL MEMGET('REAL',KDV,0,WORK,KFREE,LFREE)
      END IF
      CALL DZERO(GPO,NZ*NXOP)
      CALL PRPORB(GPO,JXOP,NXOP,PRPMO,NORBT,WORK(KDV),NZ,IH,IPRINT)
C
      CALL MEMREL('GPOGET',WORK,1,KFRSAV,KFREE,LFREE)
c       no scaling
      CALL DSCAL(NZ*NXOP,DM1,GPO,1)
      CALL QEXIT ('GPOGET')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prporb */
      SUBROUTINE PRPORB(QVEC,JXOP,NXOP,QMAT,NDIM,DV,NZ,IH,IPRINT)
C***********************************************************************
C
C     Construct 1-electron property gradient
C
C     Input:
C        JXOP    - list of orbital rotations
C        QMAT    - one-electron property matrix
C        UDV     - active density matrix (symmetry packed)
C        IH      - hermicity
C           IH = 1: g(k,l) = sum_m Dlm Qkm - ( sum_m Dmk Qml )^*
C           IH = 2: g(k,l) = sum_m Dlm Qkm + ( sum_m Dmk Qml )^*
C
C     Written by J. Thyssen May 7 1999
C     Last revision 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
      PARAMETER ( D2 = 2.0D00, D1 = 1.0D00, DM1 = -1.0D00 )
C
#include "dcborb.h"
#include "maxorb.h"
#include "dcbidx.h"
C
      DIMENSION QVEC(NXOP,NZ),JXOP(2,NXOP),QMAT(NDIM,NDIM,NZ),
     &          DV(*)

      CALL QENTER('PRPORB')
C
C
      IF (IH .EQ. 1) THEN
         DS = D1
      ELSE IF (IH .EQ. 2) THEN
         DS = DM1
      ELSE
         CALL QUIT('PRPORB: Unknown IH')
      END IF
C
      DO IG = 1,NXOP
         K = JXOP(1,IG)
         L = JXOP(2,IG)
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         ISYMK = IFSMO(K)
         ISYML = IFSMO(L)
C
C        First index:
C
         IF (ITYPK .EQ. JTINAC) THEN
C
C        ...inactive
C
            DO IZ = 1,NZ
               QVEC(IG,IZ) = QVEC(IG,IZ) + D2 * QMAT(L,K,IZ)
            END DO
C
         ELSE IF (ITYPK .EQ. JTACT) THEN
C
C        ...active
C           TODO : this code is written for MCSCF, it is not tested
C                  for open shell RPA /Trond 2001, hjaaj July 2002
C
            KDV = K - (NPSH(ISYMK) + NISH(ISYMK)) - IORB(ISYMK)
            JDV = I2ASHT(ISYMK) + (KDV-1)*NASH(ISYMK) + 1
            DO IZ = 1,NZ
               QVEC(IG,IZ) = QVEC(IG,IZ) 
     &            + DDOT(NASH(ISYMK),DV(JDV),1,QMAT(L,K,IZ),1)
               JDV = JDV + N2ASH(ISYMK)
            END DO
C
         END IF
C
C        Second index:
C
         IF (ITYPL .EQ. JTACT) THEN
C
C           ...active
C           TODO : this code is written for MCSCF, it is not tested
C                  for open shell RPA /Trond 2001, hjaaj July 2002
C
            KDV = L - (NPSH(ISYML) + NISH(ISYML)) - IORB(ISYML)
            JDV = I2ASHT(ISYML) + (KDV-1)*NASH(ISYML) + 1
            DO IZ = 1,NZ
               QVEC(IG,IZ) = QVEC(IG,IZ) 
     &            + DS * DDOT(NASH(ISYML),DV(JDV),1,QMAT(K,L,IZ),1)
               JDV = JDV + N2ASH(ISYML)
            END DO
C
         END IF
C

      END DO

      CALL QEXIT('PRPORB')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GPOTEST(PRPNAM,I1,I2,NZ,PRPAO,NB,
     &           CMO1,NF1,NO1,CMO2,NF2,NO2)
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION PRPAO(NB,NB),CMO1(NF1,NO1),CMO2(NF2,NO2)
      CHARACTER*16 PRPNAM
      WRITE(LUPRI,*) 'GPOTEST: Checking element (',I1,',',I2,
     &   ') of property ',PRPNAM
      IF(NZ.GT.1) CALL QUIT('GPOTEST: NZ. gt. 1. Write me!')
      GPO = 0.0D0
      DO J = 1,NF2
        DO I = 1,NF1
          TMP = CMO1(I,I1)*PRPAO(I,J)*CMO2(J,I2)
          GPO = GPO + TMP
          IF(ABS(TMP).GT.1e-20) THEN
            WRITE(LUPRI,'(ES16.8,2I4,3ES16.8)') 
     &        ABS(TMP),I,J,
     &        CMO1(I,I1),PRPAO(I,J),CMO2(J,I2)
          ENDIF
        ENDDO
      ENDDO
      WRITE(LUPRI,'(A,A16,A,ES16.8)') 
     &   'GPOTEST: Total ',PRPNAM,' : ',GPO
      RETURN
      END

