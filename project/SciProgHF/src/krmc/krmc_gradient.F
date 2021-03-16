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
C  /* Deck rgrdinf */
      SUBROUTINE RGRDINF(GCI,GOEE,GOEP)
C***********************************************************************
C
C     Calculate norm of gradient
C
C     Input:
C        GCI     - the CI gradient of the CEP
C        GOEE    - the e-e orbital gradient of the CEP
C        GOEP    - the e-p gradient of the CEP
C
C     Output:
C        GNORM   - (1) CI-gradient norm
C                  (2) orbital e-e gradient norm
C                  (3) orbital e-p gradient norm
C                  (4) total orbital gradient norm
C                  (5) total norm
C
C     Written by J. Thyssen - Nov 3 1998
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (D0 = 0.0D00)
C
#include "dcbopt.h"
#include "dcbkrmc_itinfo.h"
#include "dgroup.h"
C
      DIMENSION GCI(NZCONF,*),GOEE(NZHOPE,*),GOEP(NZXOPP,*)
C
      CALL QENTER('RGRDINF')
C
      IF (NZCONF .GT. 0) THEN
         GCINRM = DNORM2(NZCONFQ,GCI,1)
      ELSE
         GCINRM = D0
      END IF
      IF (NZHOPE .EQ. NZXOPE) THEN
         GOBEENRM = DNORM2(NZXOPEQ,GOEE(1,1),1)
      ELSE
         GOBEENRM = D0
         DO IZ = 1, NZ
            GOBEENRM = GOBEENRM + DDOT(NZXOPE,GOEE(1,IZ),1,GOEE(1,IZ),1)
         END DO
         GOBEENRM = SQRT(GOBEENRM)
      END IF
C
      GOBEPNRM = DNORM2(NZXOPPQ,GOEP,1)
C
      GNORM(1) = GCINRM
      GNORM(2) = GOBEENRM
      GNORM(3) = GOBEPNRM
      GNORM(4) = SQRT(GOBEENRM**2+GOBEPNRM**2)
      GNORM(5) = SQRT(GCINRM**2+GNORM(4)**2)
C
C     Save gradient information for final summary output
C
C
      CALL DCOPY(5,GNORM,1,DINFO(11),1)
C
      CALL QEXIT('RGRDINF')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rcigrad */
      SUBROUTINE RCIGRAD(CREF,FCAC,H2AC,GCI,EMCACT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate CI gradient
C
C     Input:
C        CREF    - CI vector
C        FCAC    - FC Fock matrix with active indices
C        H2AC    - Two-electron integrals with active indices
C        EMCACT  - active energy
C
C     Output:
C        GCI     - CI gradient
C
C     Written by J. Thyssen - Nov 20 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER ( THRSML = 1.0D-8 )
      PARAMETER ( THRZER = 1.0D-12 )
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION CREF(*), FCAC(*), H2AC(*), GCI(*), WORK(*)
C
      CALL QENTER('RCIGRAD')
      CALL RTKTIME(.TRUE.,10)
C
C     Calculate CI gradient
C     ---------------------
C
C
      CALL RCIGRD1(CREF,WORK(KZCONF),FCAC,H2AC,GCI,EMCACX,
     &             IPROPT,WORK,KFREE,LFREE)
C
C     Consistency check of active energy
C     ----------------------------------
C
      DNZCONF = NZCONF
      DIFFOK = SQRT(DNZCONF) * THRZER
      IF ( ABS (EMCACT) .GT. THRZER)
     &     DIFFOK = DIFFOK * ABS(EMCACT)
      IF (NZXOPT .GT. 0 .AND. ABS(EMCACT-EMCACX) .GT. DIFFOK) THEN
         IF (ABS(EMCACT-EMCACX) .GT. THRSML) THEN
            WRITE(LUPRI,9000)
     &           '*** ERROR in RCIGRAD ***',
     &           EMCACT,EMCACX
            CALL QUIT('*** ERROR in RCIGRAD ***')
         ELSE
            WRITE(LUPRI,9000)
     &           '*** WARNING in RCIGRAD ***',
     &           EMCACT,EMCACX
         END IF
      END IF
 9000 FORMAT(1P,//,A,
     &     /,1X,'(RCIGRAD) Active energy from DV and PV:',D25.15,
     &     /,1X,'(RCIGRAD) Active energy from RCIGRAD:  ',D25.15,/)
C
      CALL RTKTIME(.FALSE.,10)
      CALL QEXIT('RCIGRAD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rcigrd1 */
      SUBROUTINE RCIGRD1(CREF,IDET,FCAC,H2AC,GCI,EMCACT,
     &                   IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate CI gradient.
C
C     Input:
C        CREF    - CI vector
C        FCAC    - FC Fock matrix with active indices
C        H2AC    - Two-electron integrals with active indices
C        EMCACT  - active energy
C
C     Output:
C        GCI     - CI gradient
C
C     Written by J. Thyssen - Nov 20 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
C
      DIMENSION CREF(NZCONF,*), GCI(NZCONF,*), FCAC(*), H2AC(*)
      DIMENSION WORK(*), EMC(2)
      DIMENSION IDET(*)
C
      CALL QENTER('RCIGRD1')
      KFRSAV = KFREE
C
      MZ = MIN(2,NZ)
C
C     ***************************************************
C     *** Calculate CI matrix times reference vector. ***
C     ***************************************************
C
      CALL MEMGET('REAL',KDIAG,NZCONF,WORK,KFREE,LFREE)
C
C     Transform quarternion FCAC (Dirac MO-basis) to Molfdir MO-basis
C     ---------------------------------------------------------------
C
      CALL MEMGET('REAL',KMOLFFC,(2*NASHT)*(2*NASHT)*MIN(NZ,2),
     &            WORK,KFREE,LFREE)
      CALL QFC2MFC(FCAC,WORK(KMOLFFC),1,1,IPROPT)
C
C     Calculate CI matrix times reference vector
C     ------------------------------------------
      IF ( OPT_CIPROGRAM .EQ. 'GASCIP') THEN
C
C        Transform integrals from Dirac (NZ,3) basis to Molfdir basis
C        ------------------------------------------------------------
C
         CALL MEMGET('REAL',KMUUUU,(2*NASHT)**4 * 2, WORK,KFREE,LFREE)
         CALL DZERO(WORK(KMUUUU),(2*NASHT)**4 * 2)
C
         CALL DNZ32M(H2AC,WORK(KMUUUU),IPROPT)
C
C        Calculate CI gradient with CI program GASCIP.
         CALL GASCIP_SVC(.TRUE.,WORK(KDIAG),GCI,NZCONF,WORK(KZCONF),
     &        CREF,WORK(KMOLFFC),WORK(KMUUUU),IPROPT)
      ELSE IF ( OPT_CIPROGRAM .EQ. 'LUCIAREL') THEN
C
         CALL LUCI_SIGMA(.TRUE.,WORK(KDIAG),CREF,GCI,
     &                   WORK(KMOLFFC),H2AC,WORK,KFREE,LFREE)
C
      ELSE
         WRITE(LUPRI,'(/A/2A)')
     $        '*** ERROR in RCIGRD1 ***',
     $        'gradient not implemented for CI program ',
     $        OPT_CIPROGRAM
         CALL QUIT('*** ERROR in RCIGRD1 ***')
      END IF
C
      IF ( IPRINT .GE. 5 ) THEN
         WRITE(LUPRI,'(/A)')
     &        ' (RCIGRD1) H|0>'
         CALL OUTPUT(GCI(1,1),1,NZCONF,1,MZ,NZCONF,MZ,1,LUPRI)
         WRITE(LUPRI,'(/A)')
     &        ' (RCIGRD1) Diagonal of H'
         CALL OUTPUT(WORK(KDIAG),1,NZCONF,1,1,NZCONF,1,1,LUPRI)
      END IF
C
C
C     Calculate active energy
C     -----------------------
C
C     EACTIVE = <0|H|0> = <0 * H|0>
C                       = CREF^{\dagger} GCI
C     (since GCI is the sigma vector)
C
C     Real part of energy:
C
CSK      WRITE(LUPRI,*) 'C vector'
CSK      CALL WRTMATMN(CREF,1,NZCONF,1,NZCONF,LUPRI)
CSK      WRITE(LUPRI,*) 'sigma vector'
CSK      CALL WRTMATMN(GCI,1,NZCONF,1,NZCONF,LUPRI)
C
      EMC(1) = DDOT(NZCONF,CREF(1,1),1,GCI(1,1),1)
C
      IF (NZ .GE. 2) THEN
         EMC(1) = EMC(1) + DDOT(NZCONF,CREF(1,2),1,GCI(1,2),1)
         EMC(2) =
     &        - DDOT(NZCONF,CREF(1,2),1,GCI(1,1),1)
     &        + DDOT(NZCONF,CREF(1,1),1,GCI(1,2),1)
      ELSE
         EMC(2) = D0
      END IF
      IF (IPRINT .GE. 0 .OR. EMC(2) .NE. D0) THEN
         WRITE(LUPRI,9000)
     &        '(RCIGRD1) Active energy from CI gradient module:',
     &        'Real part: ', EMC(1),
     &        'Imag part: ', EMC(2)
         IF ( ABS(EMC(2)) .GT. THRZER) THEN
            WRITE(LUPRI,'(//,2A,//)')
     &           '*** WARNING ***',
     &           'Active energy is imaginary!'
         END IF
      END IF
C
 9000 FORMAT(1P,//,1X,A,2(/,5X,A,D25.15),//)
#ifdef UNDEF
C
C     Allocate memory for matrices:
C
      NZCONF2 = NZCONF ** 2
      CALL MEMGET('REAL',KHR,NZCONF2 * 2,WORK,KFREE,LFREE)
      KHI = KHR + NZCONF2
C
C     Construct Hamiltonian:
C     ----------------------
C
      CALL DZERO(WORK(KHR),NZCONF2)
      CALL DZERO(WORK(KHI),NZCONF2)
CTROND..check this call !!!
      CALL GASCIP_MAKEH(.TRUE.,1,NZCONF,IDET,WORK(KHR),WORK(KHI),
     &     WORK(KMOLFFC),WORK(KMUUUU))
      write(6,*) 'debug hamil'
      call prqmat(work(khr),nzconf,nzconf,nzconf,nzconf,
     &     min(nz,2),ipqtoq(1,0),lupri)
C
      mz = min(nz,2)
      call qgemm(nzconf,1,nzconf,1.0d00,
     &     'N','N',ipqtoq(1,0),work(khr),nzconf,nzconf,mz,
     &     'N','N',ipqtoq(1,0),cref,nzconf,1,mz,
     &     0.0d00,ipqtoq(1,0),gci,nzconf,1,mz)
      write(6,*) 'debug H|0>'
      call output(gci,1,nzconf,1,min(nz,2),nzconf,2,1,6)
      write(6,*) 'sigma vector'
      do i = 1, nzconf
         write(6,*) gci(I,1),gci(i,2)
      end do
C
C     Calculate active energy
C     -----------------------
C
C     EACTIVE = <0|H|0> = <0 * H|0>
C                       = CREF^{\dagger} GCI
C     (since GCI is the sigma vector)
C
C     Real part of energy:
C
      EMC(1) = DDOT(NZCONF,CREF(1,1),1,GCI(1,1),1)
C
      IF (NZ .GE. 2) THEN
         EMC(1) = EMC(1) + DDOT(NZCONF,CREF(1,2),1,GCI(1,2),1)
         EMC(2) =
     &        - DDOT(NZCONF,CREF(1,2),1,GCI(1,1),1)
     &        + DDOT(NZCONF,CREF(1,1),1,GCI(1,2),1)
      ELSE
         EMC(2) = D0
      END IF
      IF (IPRINT .GE. 0 .OR. EMC(2) .NE. D0) THEN
         WRITE(LUPRI,9000)
     &        '(RCIGRD1) Active energy from CI gradient module:',
     &        'Real part: ', EMC(1),
     &        'Imag part: ', EMC(2)
         IF ( ABS(EMC(2)) .GT. THRZER) THEN
            WRITE(LUPRI,'(//,2A,//)')
     &           '*** WARNING ***',
     &           'Active energy is imaginary!'
         END IF
      END IF
#endif
C
C
      EMCACT = EMC(1)
C
C     Finish the CI gradient
C     ----------------------
C
C     GCI -= EMCACT * CREF
C
      DO IZ = 1, MIN0(NZ,2)
         DO I = 1, NZCONF
            GCI(I,IZ) = GCI(I,IZ) - EMCACT * CREF(I,IZ)
!           SK - Aug 2010: suggestion for improvement
!           call daxpy(NZCONF,-EMCACT,CREF(1,IZ),1,GCI(1,IZ),1)
         END DO
      END DO
C
CSK      IPRINT = 20
C
      IF (IPRINT .GE. 5) THEN
         WRITE(LUPRI,'(A)')
     &        ' (RCIGRD1) CI gradient'
         CALL OUTPUT(GCI(1,1),1,NZCONF,1,MZ,NZCONF,MZ,1,LUPRI)
      END IF
C
C     Save diagonal of Hessian
C     ------------------------
C
C     E^{[2]c^*c}_{\mu,\mu} = <\mu|H|\mu> - E[0]
C        - ( g_{\mu}^* C_{\mu}^* + C_{\mu} g_{\mu} ) =
C
C     E^{[2]c^*c}_{\mu,\mu} = <\mu|H|\mu> - E[0]
C        - 2 Re ( g_{\mu} C_{\mu} )
C
      DO I = 1, NZCONF
         WORK(KDIAG + I - 1) = WORK(KDIAG + I - 1) - EMCACT
     &        - D2 * GCI(I,1) * CREF(I,1)
         IF (NZ .GE. 2) THEN
            WORK(KDIAG + I - 1) = WORK(KDIAG + I - 1)
     &           + D2 * GCI(I,2) * CREF(I,2)
         END IF
      END DO
      IF (IPRINT .GE. 5) THEN
         WRITE(LUPRI,'(A)')
     &        ' (RCIGRD1) configurational E[2] diagonal'
         CALL OUTPUT(WORK(KDIAG),1,NZCONF,1,1,NZCONF,1,1,LUPRI)
      END IF
      CALL WRTKRMC(LUKRM3,'E2CDIAG ',WORK(KDIAG),NZCONF)
C
      CALL MEMREL('E2 diag',WORK,1,KFRSAV,KFREE,LFREE)
C
CSK      IPRINT = 00
      CALL QEXIT('RCIGRD1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reogp */
      SUBROUTINE REOGP(NZXOPE,NZHOPE,NZ,G,GR)
C***********************************************************************
C
C     Reorder gradient (see comments in DIRNR).
C
C     Input:
C       G - gradient
C       NZ{X,H}OP and NZ - dimensions
C
C     Output:
C       GR - reordered gradient
C
C     Written by J. Thyssen - Jan 19 2001
C     Last revision :
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION G(NZHOPE,*), GR(*)
C
C     write(6,*) 'gradient on entry'
C     call output(g,1,nzhope,1,nz,nzhope,nz,1,6)
C
      II = 1
      DO IZ = 1, NZ
         CALL DCOPY(NZXOPE,G(1,IZ),1,GR(II),1)
         II = II + NZXOPE
      END DO
      DO IZ = 1, NZ
         CALL DCOPY(NZHOPE - NZXOPE,G(1+NZXOPE,IZ),1,GR(II),1)
         II = II + NZHOPE - NZXOPE
      END DO
C
C     write(6,*) 'gradient on exit, non-red'
C     call output(gr,1,nzxope,1,nz,nzxope,nz,1,6)
C     write(6,*) 'gradient on exit, non-red'
C     call output(gr(1+nz*nzxope),1,nzhope-nzxope,1,nz,nzhope-nzxope,
C    &     nz,1,6)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck invreogp */
      SUBROUTINE INVREOGP(NZXOPE,NZHOPE,NZ,G,GR)
C***********************************************************************
C
C     Inverse of REOGP.
C
C     Input:
C       GR - reduced gradient
C       NZ{X,H}OP and NZ - dimensions
C
C     Output:
C       G - gradient
C
C     Written by J. Thyssen - Jan 19 2001
C     Last revision :
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION G(NZHOPE,*), GR(*)
C
C     write(6,*) 'gradient on entry, non-red'
C     call output(gr,1,nzxope,1,nz,nzxope,nz,1,6)
C     write(6,*) 'gradient on entry, non-red'
C     call output(gr(1+nz*nzxope),1,nzhope-nzxope,1,nz,nzhope-nzxope,
C    &     nz,1,6)
C
      II = 1
      DO IZ = 1, NZ
         CALL DCOPY(NZXOPE,GR(II),1,G(1,IZ),1)
         II = II + NZXOPE
      END DO
      DO IZ = 1, NZ
         CALL DCOPY(NZHOPE - NZXOPE,GR(II),1,G(1+NZXOPE,IZ),1)
         II = II + NZHOPE - NZXOPE
      END DO
C
C     write(6,*) 'gradient on exit'
C     call output(g,1,nzhope,1,nz,nzhope,nz,1,6)
      RETURN
      END
