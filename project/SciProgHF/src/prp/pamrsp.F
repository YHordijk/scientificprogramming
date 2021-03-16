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
C  /* Deck xrsctl */
      SUBROUTINE XRSCTL(GPCI,GPOE,GPOP,IBTYP,IBCVC,IBEVC,IBPVC,
     &                  RCNV,EVALR,EVECR,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Module for solving the linear response equations.
C     Based on RSPCTL from RESPONS.
C
C
C     LINEQ.EQ.FALSE  direct the solution the generalized response
C     eigenvalue problem
C
C        (E[2] - w(I)S[2]) X(I) = 0
C
C     LINEQ.EQ.TRUE   direct the solution the generalized response
C     linear equations
C
C        (E[2] - w(I)S[2]) X(I) = E[1]_B
C
C        w(I)     is frequency no. I
C        E[2]     is the molecular Hessian                (h|H|h)
C        S[2]     is a metric                             (h|h)
C        E[1}_B   is gradient of perturbation operator B  (h|B)
C
C     where h is the operator manifold. In RPA h is limited to single
C     excitations.
C
C-----------------------------------------------------------------------
C     The linear reponse equation is solved by projecting X onto
C     a reduced space of trial vectors:
C
C        X = (b_i)a_i           B = [(b_1)  (b_2)  ...  (b_n)]
C
C     and solving the corresponding reduced equations
C
C        (RE[2] - wRS[2])a = 0            (LINEQ.EQ.FALSE.)
C
C        (RE[2] - wRS[2])a = RE[1]_B      (LINEQ.EQ.TRUE.)
C
C     where
C
C        RE[2] = B^{+}E[2]B
C        RS[2] = B^{+}S[2]B
C        RE[1] = B^{+}E[1}
C
C-----------------------------------------------------------------------
C
C     Variables:
C        GPCI    - property gradient , configurational part
C        GPOE    - property gradient , orbital part(e-e rotations)
C        GPOP    - property gradient , orbital part(e-p rotations)
C        REDGP   - reduced property gradient
C        REDE    - reduced Hessian matrix E[2]
C        REDS    - reduced metric         S[2]
C        IBTYP   - trial vectors will be of three kinds:
C                  C-type: contains only configurational parameters
C                  E-type: contains only e-e rotations
C                  P-type: contains only e-p rotations
C        IBCVC   - pointer from C-type trial vectors to IBTYP
C        IBEVC   - pointer from E-type trial vectors to IBTYP
C        IBPVC   - pointer from P-type trial vectors to IBTYP
C        EVALR   - eigenvalues of reduced system
C        EVECR   - eigenvectors of reduced system
C
C     Files:
C
C        LUBCI   - direct access file of C-type trial vectors
C        LUBOE   - direct access file of E-type trial vectors
C        LUBOP   - direct access file of B-type trial vectors
C
C        LUSCI   - direct access file of C-part of sigma vectors
C        LUSOE   - direct access file of E-part of sigma vectors
C        LUSOP   - direct access file of P-part of sigma vectors
C
C        LURSP   - unformatted file with restart information
C-----------------------------------------------------------------------
C
C     In the frequency-dependent case, the eigenvectors come in pairs
C
C     b_1^T = [Z  Y^*]        b_2^T = [Y  Z^*]
C
C     The paired structure is exploited in that only the upper half
C     (Z and Y respectively) of the pair of trial vectors is stored.
C
C-----------------------------------------------------------------------
C
C MAXSIM: desired number of simultaneous roots in a microiteration
C In each microiteration the reduced matrices are extended by NTSIM
C trial vectors.
C    NTSIM = NCSIM + NESIM + NPSIM
C            NCSIM - number of configurational trial vectors
C            NESIM - number of orbital trial vectors (e-e rotations)
C            NPSIM - number of orbital trial vectors (e-p rotations)
C
C IBTYP is an array running over ALL trial vectors:
C   IBTYP(1,*) - indicates type of this trial vector
C                JBCNDX =-1: CI-type
C                JBENDX =+1: electron orb-type
C                JBPNDX =+2: positron orb-type
C   IBTYP(2,*) - indicates Hermiticity of trial vector
C                +1 - [Z-Y  Y*-Z*] - Hermitian combination
C                -1 - [Z+Y  Y*+Z*] - anti-Hermitian combination
C                Note that for time-antisymmetric operators, an
C                imaginary phase is extracted, thus reversing
C                hermiticity during the calculation.
C
C  ( TRIAL VECTORS ARE ONLY ADDED FOR THE NON CONVERGED VECTORS AND
C    LINEAR DEPENDENCE IS REMOVED)
C
C MAXRIT: MAXIMUM NUMBER OF MICROITERATIONS
C
C     Written by Hans Joergen Aa. Jensen and T.Saue 1996
C     Last revision: July 30 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
      PARAMETER(D0=0.0D0)
C
      DIMENSION GPCI(*), GPOE(*),GPOP(*),
     &          IBTYP(2,NREDM), IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),
     &          RCNV(*),EVALR(*),EVECR(*),WORK(*)
C
#include "dcbgen.h"
#include "dcbxrs.h"
#include "blocks.h"
#include "dcbfir.h"
C
      CHARACTER*1 TTYPEG, TTYPSG
C
      CALL QENTER('XRSCTL')
      CALL ICOPY(2*NREDM,0,0,IBTYP,1)
      KFRSAV = KFREE
      IF(LINEQ) THEN
        NREDGP = NREDM
      ELSE
        NREDGP = 0
      ENDIF
C.....NSTAT is set to one in the static case, two otherwise
      IF (E2CHEK) STATIC = .FALSE.
      NSTAT = 2
      IF(STATIC) NSTAT = 1
C.....NDAMP is set to two in the complex case since the solution vector now
C     now has complex coefficients
      NDAMP = 1
      IF(DAMPFREQ.NE.D0) NDAMP = 2
      NEVEC=NDAMP*NREDM
      CALL MEMGET2('REAL','REDE' ,KREDE ,N2REDM,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','REDS' ,KREDS ,N2REDM,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','REDGP',KREDGP,NREDGP,WORK,KFREE,LFREE)
C
C     New NPOS(*) array for HER_PARDRV, reset dcbfir.h to true
C

      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. , .TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      if (NPOS.GT.0) THEN
         CALL MEMGET2('INTE','POS',KPOS,NPOS,WORK,KFREE,LFREE)
      else
         KPOS = KFREE
      endif
C
C
C     Fix the number of variational parameters
C
C     TTYPEG: time reversal type of E[2] G
C     TTYPSG: time reversal type of S[2] G
C
      IF (JTIMOP .GT. 0) THEN
         TTYPEG = '+'
         TTYPSG = '-'
      ELSE
         TTYPEG = '-'
         TTYPSG = '+'
      END IF
C
      NVPAR = 0
      IF(NZXOPE.GT.0) THEN
        NVPAR        = NVPAR + 1
        IPEP         = NVPAR
        BVTYP(NVPAR) = 'E'//TTYPEG
        IF(.NOT.STATIC) THEN
          NVPAR        = NVPAR + 1
          IPEM         = NVPAR
          BVTYP(NVPAR) = 'E'//TTYPSG
        ENDIF
      ENDIF
      IF(NZXOPP.GT.0) THEN
        NVPAR        = NVPAR + 1
        IPPP         = NVPAR
        BVTYP(NVPAR) = 'P'//TTYPEG
        IF(.NOT.STATIC) THEN
          NVPAR      = NVPAR + 1
          IPPM       = NVPAR
          BVTYP(NVPAR) = 'P'//TTYPSG
        ENDIF
      ENDIF
      IF(NZCONF.GT.0) THEN
        NVPAR        = NVPAR + 1
        IPCP         = NVPAR
        BVTYP(NVPAR) = 'C'//TTYPEG
        IF(.NOT.STATIC) THEN
          NVPAR        = NVPAR + 1
          IPCM         = NVPAR
          BVTYP(NVPAR) = 'C'//TTYPSG
        ENDIF
      ENDIF
      CALL MEMGET2('REAL','RNORM',KRNORM,NVPAR*NFREQ*NDAMP,
     &   WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','ISIM' ,KISIM ,NVPAR            ,
     &   WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IVECS',KIVECS,NVPAR*NFREQ*NDAMP,
     &   WORK,KFREE,LFREE)
C
      CALL XRSCT1(GPCI,GPOE,GPOP,RCNV,
     &            IBTYP,IBCVC,IBEVC,IBPVC,EVALR,EVECR,
     &            WORK(KREDE),WORK(KREDS),WORK(KREDGP),
     &            WORK(KRNORM),WORK(KIVECS),WORK(KISIM),
     &            WORK(KPOS),WORK,KFREE,LFREE)
C
      CALL MEMREL('XRSCTL.XRSCT1',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('XRSCTL')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XRSCT1 */
      SUBROUTINE XRSCT1(GPCI,GPOE,GPOP,RCNV,
     &                  IBTYP,IBCVC,IBEVC,IBPVC,EVALR,EVECR,
     &                  REDE,REDS,REDGP,RNORM,IVECS,ISIM,
     &                  NPOS,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Module for solving the linear response equations.
C     Based on RSPCTL from RESPONS.
C
C
C     LINEQ.EQ.FALSE  direct the solution the generalized response
C     eigenvalue problem
C
C        (E[2] - w(I)S[2]) X(I) = 0
C
C     LINEQ.EQ.TRUE   direct the solution the generalized response
C     linear equations
C
C        (E[2] - w(I)S[2]) X(I) = E[1]_B
C
C        w(I)     is frequency no. I
C        E[2]     is the molecular Hessian                (h|H|h)
C        S[2]     is a metric                             (h|h)
C        E[1}_B   is gradient of perturbation operator B  (h|B)
C
C     where h is the operator manifold. In RPA h is limited to single
C     excitations.
C
C-----------------------------------------------------------------------
C     The linear reponse equation is solved by projecting X onto
C     a reduced space of trial vectors:
C
C        X = (b_i)a_i           B = [(b_1)  (b_2)  ...  (b_n)]
C
C     and solving the corresponding reduced equations
C
C        (RE[2] - wRS[2])a = 0            (LINEQ.EQ.FALSE.)
C
C        (RE[2] - wRS[2])a = RE[1]_B      (LINEQ.EQ.TRUE.)
C
C     where
C
C        RE[2] = B^{+}E[2]B
C        RS[2] = B^{+}S[2]B
C        RE[1] = B^{+}E[1}
C
C-----------------------------------------------------------------------
C
C     Variables:
C        GPCI    - property gradient , configurational part
C        GPOE    - property gradient , orbital part(e-e rotations)
C        GPOP    - property gradient , orbital part(e-p rotations)
C        REDGP   - reduced property gradient
C        REDE    - reduced Hessian matrix E[2]
C        REDS    - reduced metric         S[2]
C        IBTYP   - trial vectors will be of three kinds:
C                  C-type: contains only configurational parameters
C                  E-type: contains only e-e rotations
C                  P-type: contains only e-p rotations
C        IBCVC   - pointer from C-type trial vectors to IBTYP
C        IBEVC   - pointer from E-type trial vectors to IBTYP
C        IBPVC   - pointer from P-type trial vectors to IBTYP
C        EVALR   - eigenvalues of reduced system
C        EVECR   - eigenvectors of reduced system
C
C     Files:
C
C        LUBCI   - direct access file of C-type trial vectors
C        LUBOE   - direct access file of E-type trial vectors
C        LUBOP   - direct access file of B-tpye trial vectors
C
C        LUSCI   - direct access file of C-part of sigma vectors
C        LUSOE   - direct access file of E-part of sigma vectors
C        LUSOP   - direct access file of P-part of sigma vectors
C
C        LURSP   - unformatted file with restart information
C-----------------------------------------------------------------------
C
C     In the frequency-dependent case, the eigenvectors come in pairs
C
C     b_1^T = [Z  Y^*]        b_2^T = [Y  Z^*]
C
C     The paired structure is exploited in that only the upper half
C     (Z and Y respectively) of the pair of trial vectors is stored.
C
C-----------------------------------------------------------------------
C
C MAXSIM: desired number of simultaneous roots in a microiteration
C In each microiteration the reduced matrices are extended by NTSIM
C trial vectors.
C    NTSIM = NCSIM + NESIM + NPSIM
C            NCSIM - number of configurational trial vectors
C            NESIM - number of orbital trial vectors (e-e rotations)
C            NPSIM - number of orbital trial vectors (e-p rotations)
C
C IBTYP is an array running over ALL trial vectors:
C   IBTYP(1,*) - indicates type of this trial vector
C                JBCNDX =-1: CI-type
C                JBENDX =+1: electron orb-type
C                JBPNDX =+2: positron orb-type
C   IBTYP(2,*) - indicates Hermiticity of trial vector
C                +1 - [Z-Y  Y*-Z*] - Hermitian combination
C                -1 - [Z+Y  Y*+Z*] - anti-Hermitian combination
C                Note that for time-antisymmetric operators, an
C                imaginary phase is extracted, thus reversing
C                hermiticity during the calculation.
C
C  ( TRIAL VECTORS ARE ONLY ADDED FOR THE NON CONVERGED VECTORS AND
C    LINEAR DEPENDENCE IS REMOVED)
C
C MAXITR: maximum number of microiterations
C
C     Written by Hans Joergen Aa. Jensen and T.Saue 1996
C     Last revision: July 30 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D1 = 1.0D0)
C
#include "dcbgen.h"
#include "dcbxrs.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcborb.h"
#include "cbihr2.h"
#include "dcbbas.h"
#include "dcbkrmc_itinfo.h"
C
      CHARACTER CPUTID*12,SECTID*12,DAYTID*24
      DIMENSION GPCI(*),GPOE(*),GPOP(*),
     &          IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),RCNV(*),
     &          EVALR(*),EVECR(NEVEC,NFREQ),
     &          RNORM(NVPAR,NFREQ*NDAMP),IVECS(NFREQ,NVPAR),
     &          ISIM(NVPAR),NPOS(*),WORK(*)
C
#include "dcbibn.h"
      LOGICAL BVLOAD
      CHARACTER LAB123(3)*8,RESTAR*8
C
      DATA LAB123/'********','********','XRSPONPP'/
      DATA RESTAR/'RESTART '/
C
#include "ibtfun.h"
C
      CALL QENTER('XRSCT1')
      KFRSAV = KFREE
C
C     ***************************
C     *** I n i t i a l i z e ***
C     ***************************
C
      IF ( FKRMC ) CALL GETTIM(CPUTOT1,CPUWALL1)
C
      CALL XRSINI(REDGP,REDS,REDE,GPCI,GPOE,GPOP,RCNV,EVECR,RNORM,
     &            WORK,KFREE,LFREE)
      KCONV = 1
C
C     Set up initial trial vectors
C     ============================
C     NCSIM - number of configurational trial vectors
C     NESIM - number of orbital trial vectors (e-e rotations)
C     NPSIM - number of orbital trial vectors (e-p rotations)
C     The vectors are written to direct access files.
C
      KFRITR = KFREE
      CALL XRSST(IBTYP,IBCVC,IBEVC,IBPVC,
     &           GPCI,GPOE,GPOP,EVALR,EVECR,
     &           WORK,KFREE,LFREE)
      IF (MAXSIM.LE.0) THEN
         MAXSIM1 = MAX(NTSIM,NVPAR*NFREQ)
         MXBVEC  = MAX(NZCONFQ,NZXOPEQ,NZXOPPQ)
         MAXSIM2 = (LFREE-100000) / (MXBVEC + NZVARQ + 4*N2BBASXQ)
C        ... based on rough estimate of memory per vector
C            (for b-vector, sigma-vector, some Fock/density matrices),
C            have subtracted 100000 arbitrarily for all other small matrices
C            /Sep 2005 hjaaj
         MAXSIM = MIN(MAXSIM1,MAXSIM2)
         MAXSIM = MAX(1,MAXSIM)
C        ... always attempt at least one vector
      END IF
C
C     ============================
C     hjj+panor:
C     Save screening parameter SCRFCK that is used in the construction
C     of sigma-vectors. The parameter is before change equal to the
C     value used in the wave function optimization and after the change
C     to that which should be used for properties. Both these values
C     can be user specified. Furthermore, if the trial vectors are NOT
C     normalized then the screening threshold is modified by the norm
C     of the gradient further down.
C
      SCRSAV = SCRFCK
C
C     If   (not user specified SCRFCK)
C     then SCRFCK may be dynamically determined
C          (both for LINEQ and for PP) /hjaaj jan 2001
C
      IF (E2CHEK) THEN
         SCRFCK = 1.D-14 ! we need accuracy for symmetry check of E2 matrix
      ELSE IF (.NOT. USRSCR) THEN
         IDTEST = IBTAND(INTDEF, IDFLAG)
         IF (IDTEST .NE. 0) THEN
C
C           Some Fock matrices will be calculated with direct code :
C           Set screening threshold based on SCFCNV
C           Coulomb contribution in AO basis F_ij = sum(kl) (ij|kl) D_kl
C           will contain max N2BBAST terms (N2BBAST is the total number
C           of totally symmetric "k,l" elements in D_kl; other boson
C           symmetries will have similar size)   /hjaaj 17 Jan 2001
C
C           The SQRT is based on statistical arguments for errors in
C           a sum when you omit elements abs value below a certain threshold.
C
            SCRFAC = N2BBAST
            SCRFAC = 10 * SQRT(SCRFAC)
            SCRFCK = THCXRS / SCRFAC
            WRITE(LUPRI,'(/A,1P,D10.2/A,D9.2)')
     &         '* XRSP INFO * Screening threshold in direct Fock'//
     &         ' matrix construction reset to:',SCRFCK,
     &         '         which is response convergence threshold'//
     &         ' divided by 10 * sqrt(N2BBAST):',SCRFAC
         END IF
      END IF
C
      IF (LINEQ .AND. .NOT.TKNORM) THEN
C       ... asking for relative convergence
C           thus screening threshold must also be scaled
C           otherwise we may screen (far) too much away!
        SCRFCK = SCRFCK*GPTNRM
        WRITE (LUPRI,'(A,1P,D12.4)')
     &     '* Screening threshold in LR solver scaled with norm of'//
     &     ' property vector:',GPTNRM,
     &     '  and changed to:',SCRFCK
      ENDIF
C
C     Handle case of no new start vectors
C     ===================================
      IF(NTSIM.LE.0) THEN
C
C       Restart: go to equation solver
C
        IF(NZRED.GT.0) THEN
          CALL GTINFO(DAYTID)
          WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +    '########## RESTART in MICROITERATION NO.',ITMIC,' ##########'
     +    ,DAYTID
          TIMITR = SECOND()
          JZRED = NZRED
          GO TO 30
        ELSE
C
C       No restart: exit routine.
C
          KCONV = -3
          GOTO 40
        ENDIF
      ENDIF
C
C     Handle case of zero microiterations
C     ===================================
C     Set first coefficient for each frequency equal to one
C     for getting non-zero response function with zero
C     microiterations
C     (used for Sternheimer and/or uncoupled approx.)
C
      IF(MAXITR.EQ.0.OR.UNCOUP) THEN
        NZRED  = NTSIM
        NCRED  = NCSIM
        NERED  = NESIM
        NPRED  = NPSIM
        IF (.NOT. UNCOUP) THEN
C         ... EVECR set in subroutine XRSST2 for UNCOUP
          DO I = 1,NFREQ
            EVECR(1,I) = D1
          ENDDO
        END IF
        IF(IPRXRS.GE.4) THEN
          IF(NCSIM.GT.0) THEN
            WRITE(LUPRI,'(/I5,A)')NCSIM,' Config. trial vectors'
            CALL PRBVEC(LUPRI,WORK(KBCI),NCSIM,NZCONF)
          ENDIF
          IF(NESIM.GT.0) THEN
            WRITE(LUPRI,'(/I5,A)')NESIM,
     &      ' Electronic (++) trial vectors'
            CALL PRBVEC(LUPRI,WORK(KBOE),NESIM,NZXOPE)
          ENDIF
          IF(NPSIM.GT.0) THEN
            WRITE(LUPRI,'(/I5,A)')NPSIM,
     &      ' Electronic (-+) trial vectors'
            CALL PRBVEC(LUPRI,WORK(KBOP),NPSIM,NZXOPP)
          ENDIF
         END IF
      ENDIF
C
C     **********************************************
C     ***                                        ***
C     *** M I C R O I T E R A T I O N    L O O P ***
C     ***  - add new set of trial vectors to     ***
C     ***    reduced system                      ***
C     ***                                        ***
C     **********************************************
 10   CONTINUE
        IF(ITMIC.EQ.MAXITR) GOTO 40
        TIMITR = SECOND()
        CALL FLSHFO(LUPRI)
        ITMIC = ITMIC + 1
        CALL GTINFO(DAYTID)
        WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +  '########## START MICROITERATION NO.',ITMIC,' ##########',
     +  DAYTID
        IF (NTSIM.GT.MAXSIM) THEN
          BVLOAD = .TRUE.
          NBCSIM = NCSIM
          NBESIM = NESIM
          NBPSIM = NPSIM
          CALL MEMREL('XRSCT1.BVLOAD',WORK,1,KFRITR,KFREE,LFREE)
        ELSE
          BVLOAD = .FALSE.
        ENDIF
        MTSIM = NTSIM
        ILOAD = 0
        JZRED = NZRED
        DO ITLOAD = 1,MTSIM,MAXSIM
          ILOAD = ILOAD + 1
          IF (BVLOAD) THEN
C            Trial vectors not in memory; read from file
             CALL REABLD(NBCSIM,NBESIM,NBPSIM,WORK,KFREE,LFREE)
          ENDIF
          NOSIM = NESIM + NPSIM
          NTSIM = NCSIM + NOSIM
C         Print section
          IF(BVLOAD.AND.IPRXRS.GE.1) THEN
            CALL BVINFO(ISIM,IBTYP,NZRED,NTSIM)
            WRITE(LUPRI,'(A,I3,A,6(2X,A2,A1,I3))')
     &     '* Trial vectors in load ',ILOAD,':',
     &     (BVTYP(I),':',ISIM(I),I=1,NVPAR)
          ENDIF
C
          IF(IPRXRS.GE.4) THEN
            WRITE(LUPRI,'(/A)')' Before XRSSVC'
            IF(NCSIM.GT.0) THEN
              WRITE(LUPRI,'(/I5,A)')NCSIM,' Config. trial vectors'
              CALL PRBVEC(LUPRI,WORK(KBCI),NCSIM,NZCONF)
            ENDIF
            IF(NESIM.GT.0) THEN
              WRITE(LUPRI,'(/I5,A)')NESIM,
     &        ' Electronic (++) trial vectors'
              CALL PRBVEC(LUPRI,WORK(KBOE),NESIM,NZXOPE)
            ENDIF
            IF(NPSIM.GT.0) THEN
              WRITE(LUPRI,'(/I5,A)')NPSIM,
     &        ' Electronic (-+) trial vectors'
              CALL PRBVEC(LUPRI,WORK(KBOP),NPSIM,NZXOPP)
            ENDIF
          END IF
          CALL FLSHFO(LUPRI)
C
C         Calculate sigma-vectors (E[2]B)
C           Configurational part in WORK(KSCI)
C           e-e part in WORK(KSOE)
C           e-p part in WORK(KSOP)
C         ===============================================
C
          CALL MEMGET2('REAL','SCI',KSCI,(NTSIM*NZCONFQ),
     &       WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','SOE',KSOE,(NTSIM*NZXOPEQ),
     &       WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','SOP',KSOP,(NTSIM*NZXOPPQ),
     &       WORK,KFREE,LFREE)
          CALL XRSSVC(WORK(KBCI),WORK(KBOE),WORK(KBOP),
     &         WORK(KSCI),WORK(KSOE),WORK(KSOP),
     &         IBTYP,IBCVC,IBEVC,IBPVC,NPOS,
     &         GPCI,GPOE,GPOP, WORK,KFREE,LFREE)
C
C         On file we save the Hessian PHP = PLP (i.e. beta=0), so
C         we must project CREF out.
C
          IF ( NZCONF .GT. 0 ) CALL XRSPHP(WORK(KSCI),WORK(KFREE),LFREE)
C
          CALL FLSHFO(LUPRI)
C
C         Write sigma-vectors to file
C         ============================
C
          ICOUNT = NZRED
          CALL WRSVEC(ICOUNT,WORK(KSCI),WORK(KSOE),WORK(KSOP))
C
C         Update restart info file
C         ========================
C
          REWIND (LURSP)
          KTOT = NZRED + LOFFTY
          IF (IPRXRS.GE.7) THEN
            WRITE(LUPRI,'(/A,I5)')' Trial vectors for restart:',KTOT
          ENDIF
          WRITE(LURSP) LAB123,RESTAR
          WRITE(LURSP) JSYMOP,JTIMOP,KTOT,
     &                 (IBTYP(1,I),IBTYP(2,I),I=1,KTOT)
          REWIND (LURSP)
C
C         Update matrices and vectors of reduced RSP equation
C         ===================================================
C
          IF (NZRED.GT.NREDM) THEN
             WRITE(LUPRI,'(//A/A,I5,/A,I5)')
     *         ' >>> ERROR IN XRSCTL >>>',
     *         ' Dimension of reduced space is  ',NZRED,
     *         ' which exceeds allowed dimension',NREDM
            CALL QUIT('XRSCTL: Too large dimension of reduced space')
          END IF
          IZRED  = NZRED
          NZRED  = NZRED + NTSIM
          ICRED  = NCRED
          NCRED  = NCRED + NCSIM
          IERED  = NERED
          NERED  = NERED + NESIM
          IPRED  = NPRED
          NPRED  = NPRED + NPSIM
          CALL FLSHFO(LUPRI)
          CALL XRSRED (REDGP,REDE,REDS,GPCI,GPOE,GPOP,
     &                 IBCVC,IBEVC,IBPVC,IBTYP,
     &                 WORK(KBCI),WORK(KBOE),WORK(KBOP),
     &                 WORK(KSCI),WORK(KSOE),WORK(KSOP),
     &                 WORK,KFREE,LFREE)
          CALL MEMREL('XRSCT1.XRSRED',WORK,1,KFRITR,KFREE,LFREE)
          CALL FLSHFO(LUPRI)
C
C       Perform a new load if some trial vectors still remains
C       ======================================================
          IF (BVLOAD) THEN
            NBCSIM = NBCSIM - NCSIM
            NBESIM = NBESIM - NESIM
            NBPSIM = NBPSIM - NPSIM
          ENDIF
        ENDDO
C
 30     CONTINUE
C
C       Info on current trial vectors
C       -----------------------------
C
        CALL BVINFO(ISIM,IBTYP,JZRED,MTSIM)
C
C       Write info to PAMCYC
C       ====================
C
        WRITE(LUCYCL,'(I3,2X,A9,6(2X,A2,A1,I3))')
     &     ITMIC,XRS_INTTYP,(BVTYP(I),':',ISIM(I),I=1,NVPAR)
        WRITE(LUPRI,'(A,6(2X,A2,A1,I3))')
     &     '* Trial vectors in this micro iteration:',
     &     (BVTYP(I),':',ISIM(I),I=1,NVPAR)
C
C       Solve reduced RSP equation
C       ==========================
C
        CALL XRSSOL (REDGP,REDE,REDS,EVALR,EVECR,RCNV,
     &               IBTYP,IBEVC,WORK,KFREE,LFREE)
        CALL FLSHFO(LUPRI)
C
C        CALCULATE TRANSITION MOMENT
C
C        CALL XRSMOM(IBTYP,EVALR,EVECR,BVECS,SVECS,PRVEC,NSIM)
C
C SET UP REDUCED MATRIX EXPLICLY AND CHECK BLOCK STRUCTURE IN E2 AND S2
C
C         IF (ABSYM) THEN
C            CALL XE2SYM(IBTYP,WORK,KFREE,LFREE)
C         END IF
C
C        Generate MAXSIM linear independent trial vectors for
C        next microiteration using solutions of reduced RSP equation
C        ===========================================================
C
         CALL XRSNEX(IBTYP,IBCVC,IBEVC,IBPVC,RCNV,
     &               EVALR,EVECR,GPCI,GPOE,GPOP,
     &               RNORM,IVECS,WORK,KFREE,LFREE)
C
C
C        Zero reduced space if new integrals contribute
C        ==============================================
C
         IF(INTFLG.NE.INTBUF) THEN
           N2RED = NZRED*NREDM
           IF(LINEQ) CALL DZERO(REDGP,NZRED)
           CALL DZERO(REDE,N2RED)
           CALL DZERO(REDS,N2RED)
           NERED = 0
           NPRED = 0
           NCRED = 0
           NZRED = 0
         ENDIF
 137     CONTINUE
         TIMITR = SECOND() - TIMITR
         CPUTID = SECTID(TIMITR)
C
C        Write info to PAMCYC
C        ====================
C
         WRITE(LUCYCL,'(1P,6D10.2)')
     &                        ((RNORM(I,J),I=1,NVPAR),J=1,NFREQ*NDAMP)
         WRITE(LUCYCL,'(A12)') CPUTID
         IF(IPRXRS.GE.1) THEN
           WRITE(LUPRI,'(/A,I5,A,A12)')
     &       '* End of microiteration no.',ITMIC,' after ',CPUTID
         ENDIF
         CALL FLSHFO(LUPRI)
C         GOTO 138
CTROND         OPTORB = OPTSAV
C
C        Check for state of the calculation
C        ==================================
C
         IF(KCONV.GT.0.AND.ITMIC.LT.MAXITR) THEN
           INTBUF = INTFLG
           GOTO 10
         ENDIF
C
C****************************************************
C****                                            ****
C**** E N D  O F   M I C R O I T E R A T I O N S ****
C****                                            ****
C****************************************************
C
 40   CONTINUE
CPP
 138  CONTINUE
      CALL XRSOUT(RCNV,RNORM,ISIM,EVALR,IBTYP)
      IF(NZCONF.GT.0) THEN
        CLOSE(LUBCI,STATUS = 'KEEP')
        CLOSE(LUSCI,STATUS = 'DELETE')
      ENDIF
      IF(NZXOPE.GT.0) THEN
        CLOSE(LUBOE,STATUS = 'KEEP')
        IF (E2CHEK) THEN
           CLOSE(LUSOE,STATUS = 'KEEP')
        ELSE
           CLOSE(LUSOE,STATUS = 'DELETE')
        END IF
      ENDIF
      IF(NZXOPP.GT.0) THEN
        CLOSE(LUBOP,STATUS = 'KEEP')
        IF (E2CHEK) THEN
           CLOSE(LUSOP,STATUS = 'KEEP')
        ELSE
           CLOSE(LUSOP,STATUS = 'DELETE')
        END IF
      ENDIF
      CLOSE(LUCYCL,STATUS='DELETE')
      CALL MEMREL('XRSCT1',WORK,1,KFRSAV,KFREE,LFREE)
C
C     Restore the screening parameter SCRFCK
      SCRFCK = SCRSAV
C
C
C     For KR-MCSCF iterations, calculate predicted energy
C     change and update some info elements.
C
      IF ( FKRMC ) THEN
C
C        Calculate predicted energy change
C
         CALL RSNRPRD(REDGP,REDE,EVECR,NREDM,IBTYP,NZRED)
C
C        Dump info needed for backsteps
C
         REWIND(LUKRM4)
         CALL WRTLRED(LUKRM4,NZRED,NREDM,REDE,REDGP,EVECR,
     &        IBTYP,IBEVC,IBPVC,IBCVC,NPRED)
C
C        Save info for final summary output
C
         CALL GETTIM(CPUTOT2,CPUWALL2)
         DINFO(17) = CPUTOT2-CPUTOT1
         IINFO(2) = ITMIC
         IINFO(3) = NZRED
         IINFO(4) = NCRED
         IINFO(5) = NERED
         IINFO(6) = NPRED
      END IF
C
C
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('XRSCT1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XRSST */
      SUBROUTINE XRSST(IBTYP,IBCVC,IBEVC,IBPVC,
     &                 GPCI,GPOE,GPOP,EVALR,EVECR,
     &                 WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate start vectors in RSP
C
C     Written by H.J.Aa.Jensen and T.Saue 1996
C     Last revision : July 30 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),
     &          GPCI(*),GPOE(*),GPOP(*),
     &          EVALR(*),EVECR(*),WORK(*)
C
#include "dcbibn.h"
C
#include "dcbxrs.h"
#include "dgroup.h"
C
      CALL QENTER('XRSST')
      KFRSAV = KFREE
C
      NESIM = 0
      NPSIM = 0
      NCSIM = 0
C
C.....Generate complete set of trial functions
      IF(E2CHEK) THEN
        IFAC = NZ*NSTAT
        NESIM = NZXOPE*IFAC
        NPSIM = NZXOPP*IFAC
        NCSIM = NZCONF*IFAC
      ELSEIF(LINEQ) THEN
C       ... linear response
C.......DIAGHE: diagonal Hessian calculated => preconditioning
C       RSREST: restart on solution vectors from previous run
C       UNCOUP: uncoupled response; solution vectors generated here
        IF(DIAGHE.OR.RSREST.OR.UNCOUP) THEN
          NVEC = NFREQ
C.......no preconditioning; use gradient only for trial vectors
        ELSE
          NVEC = 1
        ENDIF
        IF(NZXOPE.GT.0) NESIM = NVEC*NSTAT*NDAMP
        IF(NZXOPP.GT.0) NPSIM = NVEC*NSTAT*NDAMP
        IF(NZCONF.GT.0) NCSIM = NVEC*NSTAT*NDAMP
      ELSE
C       ... linear response single residue; excitation energies
CHJ TODO nov.2000:
C       In this version prefer conf start vectors,
C       if not available then use e-e start vectors (i.e. for RPA)
        NEXST2 = NEXSTV + 2
C       ... we include two extra trial vectors to check if
C           elements NEXSTV and NEXSTV+1/+2 are (nearly) degenerate
        IF (NZCONF.GT.0) THEN
           NEXST2 = MIN(NEXST2,NZCONF)
           NCSIM = 2*NEXST2
        ELSE IF (NZXOPE.GT.0) THEN
           NEXST2 = MIN(NEXST2,NZXOPEQ)
           NESIM = 2*NEXST2
        ELSE
           CALL QUIT('XRSST-EXCITA: No conf or e-e variables')
        END IF
      ENDIF
C
      NTSIM = NESIM + NPSIM + NCSIM
C.....Check for out of bounds
C reminder NSTAT = 1 in static case, = 2 otherwise
C          NDAMP = 1 in real or imag case, = 2 in complex case
      NSIM = NFREQ * NSTAT * NDAMP
C
      NBUF = 0
      IF(NZXOPE.GT.0) NBUF = NBUF + 1
      IF(NZXOPP.GT.0) NBUF = NBUF + 1
      IF(NZCONF.GT.0) NBUF = NBUF + 1
      NTOT = NZRED+NTSIM+NSIM*NBUF
      IF(NTOT.GT.NREDM) THEN
        WRITE(LUPRI,'(A,I12)')
     &    '* Current number of trial vectors       :',NZRED,
     &    '* Requested number of new trial vectors :',NTSIM,
     &    '* Solution vectors generated in XRSNEX  :',NSIM*NBUF,
     &    '* Maximum number of trial vectors NREDM :',NREDM,
     &    '* Increase MAXRED to                    :',NTOT
        CALL QUIT('* XRSST: Increase size of reduced space MAXRED !')
      ENDIF
      NBTOT = NESIM*NZXOPEQ + NPSIM*NZXOPPQ + NCSIM*NZCONFQ
      CALL MEMGET2('REAL','BTOT',KBTOT,NBTOT,WORK,KFREE,LFREE)
      KFIRST = KBTOT
C
      KBOE   = KFIRST
      KFIRST = KFIRST + NZXOPEQ*NESIM
C
      KBOP   = KFIRST
      KFIRST = KFIRST + NZXOPPQ*NPSIM
C
      KBCI   = KFIRST
      KFIRST = KFIRST + NZCONFQ*NCSIM
C
      KFRBUF = KFREE
C     ... set proper dimension for integer arrays. default. SK - Sep 08
!     mi : +1 memleak fix
      LEVEC  = NESIM +1
      LPVEC  = NPSIM +1
      LCVEC  = NCSIM

!     krmc optimization
      if(fkrmc)then
!       ... SK - Sep 08: if XRSST is called as part from the MCSCF
!       optimization procedure we need to allocate at least 2* NxSIM where
!       x = E (e-e), P (e-p) or C (CI gradient).
!       In subroutine XRSST2 the corresponding integer array is accessed
!       by
!       IVEC(IP) = x
!       IVEC(IM) = x
!       where IP = NFREQ + IM and IM = 1, NFREQ
!
!       example for bad access to IVEC which caused the crash for GOSCI
!       (and in principle also for LUCIAREL):
!
!       NESIM = 1 ==> LEVEC = 1 --> DIMENSION IVEC(1)
!       NFREQ = 1 ==> IM = 1, IP = 2
!
!       allowed access:       IVEC(IM) = ...
!       resulting bad access: IVEC(IP) = ...
!
!       FIXME: does this also apply to other cases with STATIC = .TRUE.?
!       Do we need a "if(fkrmc)then" statement?
!       SK - Sep 08
!
        levec  = nesim * 2
        lpvec  = npsim * 2
        lcvec  = ncsim * 2
      end if

      CALL MEMGET2('INTE','EVEC',KEVEC,LEVEC,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','PVEC',KPVEC,LPVEC,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','CVEC',KCVEC,LCVEC,WORK,KFREE,LFREE)
C     ... initialize
      CALL IZERO(WORK(KEVEC),LEVEC)
      CALL IZERO(WORK(KPVEC),LPVEC)
      CALL IZERO(WORK(KCVEC),LCVEC)
C
C     If requested, allocate space for diagonal elements of Hessian
C     (approximate)
C
      IF(DIAGHE) THEN
         CALL MEMGET2('REAL','EDIAH',KEDIAH,NZXOPE,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','EDIAS',KEDIAS,NZXOPE,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','PDIAH',KPDIAH,NZXOPP,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','PDIAS',KPDIAS,NZXOPP,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','CDIAH',KCDIAH,NZCONF,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET2('REAL','EDIAH',KEDIAH,0,WORK,KFREE,LFREE)
         KEDIAS = KEDIAH
         KPDIAH = KEDIAH
         KPDIAS = KEDIAH
         KCDIAH = KEDIAH
      ENDIF
      CALL XRSST1(IBTYP,IBCVC,IBEVC,IBPVC,
     &            WORK(KEVEC),WORK(KPVEC),WORK(KCVEC),
     &            GPCI,GPOE,GPOP,EVALR,EVECR,
     &            WORK(KEDIAH),WORK(KPDIAH),WORK(KCDIAH),
     &            WORK(KEDIAS),WORK(KPDIAS),
     &            WORK,KFREE,LFREE)

      NOSIM = NESIM + NPSIM
      NTSIM = NCSIM + NOSIM
C
C     Compress trial vectors in memory
C     ================================
C
      MBTOT = NESIM*NZXOPEQ + NPSIM*NZXOPPQ + NCSIM*NZCONFQ
      IF(MBTOT.NE.NBTOT) THEN
C     ...compress vectors in memory
        KFIRST = KBTOT
C
        JBOE   = KBOE
        KBOE   = KFIRST
        NEDIM  = NZXOPEQ*NESIM
        KFIRST = KFIRST + NEDIM
        IF(KBOE.NE.JBOE.AND.NEDIM.GT.0) THEN
Chj       CALL DCOPY(NEDIM,WORK(JBOE),1,WORK(KBOE),1)
Chj       Oct 2000: dcopy may cause problems because of overlap
chj       between WORK(KBOE:KBOE+NEDIM-1) and WORK(JBOE:JBOE+NEDIM-1)
          DO I = 0,NEDIM-1
             WORK(KBOE+I) = WORK(JBOE+I)
          END DO
        END IF
C
        JBOP   = KBOP
        KBOP   = KFIRST
        NPDIM  = NZXOPPQ*NPSIM
        KFIRST = KFIRST + NPDIM
        IF(KBOP.NE.JBOP.AND.NPDIM.GT.0) THEN
Chj       CALL DCOPY(NPDIM,WORK(JBOP),1,WORK(KBOP),1)
          DO I = 0,NPDIM-1
             WORK(KBOP+I) = WORK(JBOP+I)
          END DO
        END IF
C
        JBCI   = KBCI
        KBCI   = KFIRST
        NCDIM  = NZCONFQ*NCSIM
        KFIRST = KFIRST + NCDIM
        IF(KBCI.NE.JBCI.AND.NCDIM.GT.0) THEN
Chj       CALL DCOPY(NCDIM,WORK(JBCI),1,WORK(KBCI),1)
          DO I = 0,NCDIM-1
             WORK(KBCI+I) = WORK(JBCI+I)
          END DO
        END IF
C
cmi      deallocate all
        CALL MEMREL('XRSINI.CMP',WORK,1,KFRSAV,KFREE,LFREE)
        CALL MEMGET2('REAL','BTOT',KBTOT,NBTOT,WORK,KFREE,LFREE)
      ELSE
        CALL MEMREL('XRSINI.BUF',WORK,1,KFRBUF,KFREE,LFREE)
      ENDIF
C
      CALL QEXIT ('XRSST')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XRSST1 */
      SUBROUTINE XRSST1(IBTYP,IBCVC,IBEVC,IBPVC,
     &                 IEVEC,IPVEC,ICVEC,GPCI,GPOE,GPOP,EVALR,EVECR,
     &                 HEDIAG,HPDIAG,HCDIAG,SEDIAG,SPDIAG,
     &                 WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate start vectors in RSP
C
C     Written by H.J.Aa.Jensen and T.Saue 1996
C     Last revision : July 30 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "dcbibn.h"
C
#include "dcbxrs.h"
#include "dgroup.h"
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),EVALR(*),EVECR(*),
     &          GPCI(*),GPOE(*),GPOP(*),IEVEC(NDAMP,*),IPVEC(NDAMP,*),
     &          ICVEC(*),HEDIAG(*),HPDIAG(*),HCDIAG(*),
     &          SEDIAG(*),SPDIAG(*),WORK(*)
C
      CALL QENTER('XRSST1')
      KFRSAV = KFREE
C
C     If requested, get diagonal elements of Hessian (approximate)
C        and metric.
C
      IF(DIAGHE) CALL E2DIAG(HEDIAG,HPDIAG,HCDIAG,SEDIAG,SPDIAG,
     &                       WORK,KFREE,LFREE)
C
C     L I N E A R   R E S P O N S E   F U N C T I O N S
C     =================================================
C     (or single residue for excitations)
C
C
      IND = NZRED
C
C     Orbital (e-e) type trial vector
C     -------------------------------
C
      IF(NESIM.GT.0) THEN
CPP.....test for which case you need this allocation
        CALL MEMGET2('REAL','EPREC',KPREC,2*NDAMP*NZXOPE,
     &     WORK,KFREE,LFREE)
        CALL XRSST2('E',WORK(KBOE),GPOE,HEDIAG,SEDIAG,EVALR,EVECR,
     &              IBTYP,IBEVC,IEVEC,JBENDX,NESIM,IND,NERED,NZXOPE,NZ,
     &              WORK(KPREC),WORK,KFREE,LFREE)
        CALL MEMREL('XRSST1.eprec',WORK,1,KPREC,KFREE,LFREE)
        IND = IND + NESIM
      END IF
C
C     Orbital (e-p) type trial vector
C     -------------------------------
C
      IF(NPSIM.GT.0) THEN
CPP.....test for which case you need this allocation
        CALL MEMGET2('REAL','PPREC',KPREC,2*NDAMP*NZXOPP,
     &     WORK,KFREE,LFREE)
        CALL XRSST2('P',WORK(KBOP),GPOP,HPDIAG,SPDIAG,EVALR,EVECR,
     &              IBTYP,IBPVC,IPVEC,JBPNDX,NPSIM,IND,NPRED,NZXOPP,NZ,
     &              WORK(KPREC),WORK,KFREE,LFREE)
        CALL MEMREL('XRSST1.pprec',WORK,1,KPREC,KFREE,LFREE)
        IND = IND + NPSIM
      END IF
C
C     Configurational type trial vector
C     ---------------------------------
C
      IF(NCSIM.GT.0) THEN
        CALL MEMGET2('REAL','CIDIAS',KCDIAS,0,WORK,KFREE,LFREE)
        KPREC=KFREE
        CALL XRSST2('C',WORK(KBCI),GPCI,HCDIAG,WORK(KCDIAS),
     &              EVALR,EVECR,IBTYP,IBCVC,ICVEC,
     &              JBCNDX,NCSIM,IND,NCRED,NZCONF,NZ_in_CI,
     &              WORK(KPREC),WORK,KFREE,LFREE)
        IND = IND + NCSIM
      END IF

      CALL MEMREL('XRSST1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT ('XRSST1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XRSORT */
      SUBROUTINE XRSORT(BVECS,NBTYP,NBNEW,NBPREV,
     &                  IBVEC,IBTYP,IVECS,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Purpose:
C     Orthogonalize new b-vectors against all previous b-vectors
C     and among themselves, and, if requested, renormalize.
C
C     b-vectors generally have the structure
C
C              |  B  |
C       b(h) = | hB* |
C
C     where h = +/- 1 and only the upper part B is explicitly stored.
C
C     b-vectors of opposite h-parameter are automatically
C     orthogonal, provided that the have identical time reversal symmetry
C     structure, so only b-vectors of the same h need to be orthogonbalized
C     against each others.
C
C     Input:
C       BVECS,  non-orthogonal b-vectors
C       NBNEW,  number of b-vectors in BVECS
C       NBPREV, number of previous b-vectors on file
C       IBTYP(1,*) = NBTYP   b-vector type (CI, e-obital, p-orbital)
C       IBTYP(2,*) = H       h-parameter
C       IBVEC  - pointer to IBTYP
C
C     Output:
C       BVECS,  orthogonal b-vectors (also written to file)
C       NBNEW,  number of b-vectors written to file
C
C     Scratch:
C       BVBUF(*) - buffer for old b-vectors
C
C     Written by H.J.Aa.Jensen and T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbibn.h"
#include "thrldp.h"
      DIMENSION BVECS(*),IBVEC(*),IBTYP(2,*),IVECS(*),WORK(*)
C
C
C Used from common blocks:
C  WRKXRS: LZCONF,LZXOPT,JSYMOP
C
#include "dgroup.h"
#include "dcbxrs.h"
C
      CALL QENTER('XRSORT')
C
C     Set dimensions;
C     Set b vector input/output unit;
C     Set linear dependency limit =
C         sqrt( number of elements * machine # of digits squared ) =
C         sqrt( number of elements ) * machine # of digits
C
      KFRSAV = KFREE
CTROND : Remember LOFFTY !!!!!
      IF    (NBTYP.EQ.JBCNDX) THEN
        NBDIM  = NZCONFQ
        LUBVEC = LUBCI
      ELSEIF(NBTYP.EQ.JBENDX) THEN
        NBDIM  = NZXOPEQ
        LUBVEC = LUBOE
      ELSEIF(NBTYP.EQ.JBPNDX) THEN
        NBDIM  = NZXOPPQ
        LUBVEC = LUBOP
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'XRSORT ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('XRSORT: Unknown NBTYP !')
      ENDIF
      IF (LINEQ) THEN
Chj-nov99: If differential screening (TKNORM false) then
C     only keep trial vectors which norm has been reduced by less than
C     a factor of THRLDV under Gram-Schmidt orthogonalization
Chj-jan01: what has this to do with TKNORM ??? now always
C          THRLDV2 = 5.0D-3
C     -- the THRLDV1 is used to remove for true linear dependency
C     -- the THRLDV2 is used to reduce the number of trial vectors
C        for linear response calculations with many close-lying frequencies.
chj   IF (TKNORM) THEN
chj      THRLDV = NBDIM*THRLDP
chj      THRLDV = SQRT(THRLDV)
         THRLDV1 = NBDIM*THRLDP
         THRLDV1 = SQRT(THRLDV1)
chj   ELSE
chj      THRLDV = 5.0D-3
         THRLDV2 = 5.0D-3
chj   END IF
      ELSE
C     ... hjaaj Oct 2004: for excitations we don't want
C         to use the THRLDV2 test (we don't have close-lying solutions
C         vectors as for LINEQ with close-lying frequencies,
C         and it has been seen to destroy the paired structure
C         by removing a P(-) trial vector but not the
C         corresponding P(+) trial vector.
C
         THRLDV1 = NBDIM*THRLDP
         THRLDV1 = SQRT(THRLDV1)
         THRLDV2 = 0.0D0
      END IF
C
C     Memory allocation
C
      CALL MEMGET2('REAL','BVBUF',KBVBUF,NBDIM,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','BVNRM',KBVNRM,NBNEW,WORK,KFREE,LFREE)
      CALL XRSOR1(BVECS,NBDIM,NBNEW,NBPREV,THRLDV1,THRLDV2,IBVEC,IBTYP,
     &            WORK(KBVBUF),WORK(KBVNRM),LUBVEC,IVECS,NBTYP)
      CALL MEMREL('XRSORT',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('XRSORT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck XRSOR1 */
      SUBROUTINE XRSOR1(BVECS,NBDIM,NBNEW,NBPREV,THRLDV1,THRLDV2,IBVEC,
     &                  IBTYP,BVBUF,BVNRM,LUBVEC,IVECS,NBTYP)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (T1MIN = 1.D-8, THRRND = 1.D-4)
C
#include "dcbibn.h"
#include "dcbxrs.h"
      DIMENSION BVECS(NBDIM,NBNEW),IBTYP(2,NREDM),IBVEC(NREDM),
     &     BVBUF(NBDIM),IVECS(NBNEW),BVNRM(NBNEW)
      CHARACTER BTYP(-1:2)*1,XTYP(-1:1)*3
      DATA BTYP/'C',' ','E','P'/,XTYP/'(-)','   ','(+)'/
!     real(8), parameter :: thr_vec = 1.0d-03
      real(8), parameter :: thr_vec = 0.0d+00
      integer            :: iround_screening, izero, i_removed

      IROUND_screening = 0
  5   CONTINUE
      IROUND_screening = IROUND_screening + 1
      I_ZERO           = 0
      I_REMOVED        = 0
C
C
C     Take norm of new trial vectors
C     ==============================
C
      DO IBN = 1,NBNEW
         IF(IVECS(IBN).EQ.0) THEN
            BVNRM(IBN) = D0
         ELSE
            BVNRM(IBN) = D2*DNRM2(NBDIM,BVECS(1,IBN),1)
C
C           Check for zero norm
C
            IF(BVNRM(IBN).EQ.D0) THEN
               IVECS(IBN) = 0
            ENDIF
         ENDIF
      ENDDO
C
C     IBSTART:
C        For orbital trial vectors IBSTART = 1, but for
C        conf. trial vectors we also need to orthogonalize against
C        CREF. CREF is conviently written to LUBCI in record 1.
C
      IF ( NBTYP .EQ. JBCNDX) THEN
         IBSTART = 0
      ELSE
         IBSTART = 1
      END IF
C
      IROUND = 0
      ITURN  = 0
 10   CONTINUE
      IREC   = 0
      ITURN  = ITURN + 1
C
C     Orthogonalize new b-vectors against previous b-vectors
C     ======================================================
C
      DO IB = IBSTART,NBPREV
C
C        Read old vector from LUBVEC
C
         IREC = IREC + 1
         CALL READ_DX(LUBVEC,IREC,NBDIM,BVBUF)
         IF (IB .NE. 0) THEN
            IH   = IBTYP(2,IBVEC(IB))
         ELSE
            IH = 1
         END IF
         IF(.NOT.TKNORM) THEN
            BNBN = D4*DNRM2(NBDIM,BVBUF,1)**2
         ELSE
            BNBN = D1
         ENDIF
C
         DO 100 JBN = 1,NBNEW
            JB = JBN + NBPREV
            JH = IBTYP(2,IBVEC(JB))
            IF(IVECS(JBN).EQ.0.OR.(IH.NE.JH)) GOTO 100
C
C           Form -<BVECS|BVBUF>/<BVBUF|BVBUF> :
C
            BNBP = -D4*DDOT(NBDIM,BVECS(1,JBN),1,BVBUF,1)/BNBN
C
C           Form BVECS = BVECS -<BVECS|BVBUF>/<BVBUF|BVBUF>BVBUF
C
            CALL DAXPY(NBDIM,BNBP,BVBUF,1,BVECS(1,JBN),1)
 100     CONTINUE
      ENDDO
C
C     Orthogonalize new vectors against each others
C     =============================================
C
      DO 200 IBN = 1,NBNEW
        IB = IBN + NBPREV
        IH = IBTYP(2,IBVEC(IB))
        IF(IVECS(IBN).EQ.0) THEN
          IF(IPRXRS.GE.2) THEN
            WRITE(LUPRI,'(A1,A3,A,I3,A)')
     &      BTYP(NBTYP),XTYP(IH),'-vector no.',IBN,' not accepted'
          ENDIF
          GOTO 200
        ENDIF
C
        BNRM1 = D2*DNRM2(NBDIM,BVECS(1,IBN),1)
        DO 300 JBN = 1,(IBN-1)
          JB = JBN + NBPREV
          JH = IBTYP(2,IBVEC(JB))
          IF(IVECS(JBN).EQ.0.OR.(IH.NE.JH)) GOTO 300
          IF(.NOT.TKNORM) THEN
            BNBN = D4*DNRM2(NBDIM,BVECS(1,JBN),1)**2
          ELSE
            BNBN = D1
          ENDIF
C
C         Form -<BVECS1|BVECS2>/<BVECS2|BVECS2> :
C
          BNBP = -D4*DDOT(NBDIM,BVECS(1,JBN),1,BVECS(1,IBN),1)/BNBN
C
C         Form BVECS1 = BVECS1 -<BVECS1|BVECS2>/<BVECS2|BVECS2>BVECS2
C
          CALL DAXPY(NBDIM,BNBP,BVECS(1,JBN),1,BVECS(1,IBN),1)
 300    CONTINUE
C
C       Take norm of vector number IBN;
C
        BNRM = D2*DNRM2(NBDIM,BVECS(1,IBN),1)
C
C       Check size of norm
C       ==================
C
        IF ( BNRM.LT.THRLDV1*BVNRM(IBN) ) THEN
C          ... see comment in XRSORT /hjaaj Jan 2001
C
          WRITE(LUPRI,'(A1,A3,A,I3,A,1P,D13.5)')
     &      BTYP(NBTYP),XTYP(IH), '-vector no.',IBN,
     &      ' removed (linear dependence). Norm :',BNRM
          IVECS(IBN) = 0
        ELSE IF ( BNRM.LT.THRLDV2*BNRM1 ) THEN
          WRITE(LUPRI,'(A1,A3,A,I3,A,1P,D13.5)')
     &      BTYP(NBTYP),XTYP(IH), '-vector no.',IBN,
     &      ' removed (close to previous accepted vector). Norm :',BNRM
          IVECS(IBN) = 0
        ELSE
          IF(BNRM.LT.THRRND*BVNRM(IBN).AND.ITURN.EQ.1) IROUND=IROUND+1
          IF(TKNORM) THEN
            IF(BNRM.LT.T1MIN) THEN
              BNRM = D1 / BNRM
              CALL DSCAL(NBDIM,BNRM,BVECS(1,IBN),1)
              BNRM = D2*DNRM2(NBDIM,BVECS(1,IBN),1)
            END IF
            BNRM = D1 / BNRM
            CALL DSCAL(NBDIM,BNRM,BVECS(1,IBN),1)
            IF(IPRXRS.GE.2) THEN
              WRITE(LUPRI,'(A1,A3,A,I3,A,1P,D12.5)')
     &          BTYP(NBTYP),XTYP(IH),
     &          '-vector no.',IBN,' accepted. Norm : ',D1
            ENDIF
          ELSE
C           hj+patrick:
C           In case of differential screening: keep the original norm
C           (before Gram-Schmidt)
            BNRM = BVNRM(IBN)/BNRM
            CALL DSCAL(NBDIM,BNRM,BVECS(1,IBN),1)
            IF(IPRXRS.GE.2) THEN
               WRITE(LUPRI,'(A1,A3,A,I3,A,1P,D12.5)')
     &              BTYP(NBTYP),XTYP(IH),
     &              '-vector no.',IBN,' accepted. Norm : ',BVNRM(IBN)
            ENDIF
          ENDIF
        ENDIF
        IF (THR_VEC .GT. 0.0D0 .AND. IROUND_screening .EQ. 1)THEN
          DO I = 1, NBDIM
            IF(ABS(BVECS(i,IBN)) .LT. THR_VEC)THEN
              IF(BVECS(i,IBN) .EQ. 0.0d0) THEN
                I_ZERO = I_ZERO + 1
              ELSE
                BVECS(i,IBN) = 0.0d0
                I_REMOVED    = I_REMOVED + 1
              END IF
            END IF
          END DO
        END IF
 200  CONTINUE
C
C     Check whether a new round of orthonormalization is
C     necessary
C     ==================================================
C
      IF(THR_VEC .GT. 0.0D0 .AND. IROUND_screening .EQ. 1) THEN
         IF(I_REMOVED .GT. 0) THEN
!#ifdef LUCI_DEBUG
            WRITE (lupri,'(/A,I12,A,A,1P,D10.2,I14)')
     &       'info: Removed',I_REMOVED,
     &       ' elements in new trial vector',
     &       '; threshold & zeroes',THR_VEC, I_ZERO
!#endif
            ITURN     = 0
            I_REMOVED = 0
            I_ZERO    = 0
            GO TO 5
         END IF
      ELSE
        IF(IROUND.GT.0) THEN
           IROUND = 0
           IF (ITURN.EQ.1) THEN
              WRITE(LUPRI,'(A)')
     &             '* XRSORT: Second round of orthonormalization !'
              GO TO 10
           ENDIF
           CALL QUIT(
     &          'XRSORT error: IROUND .gt. 0 and ITURN .gt. 1')
        END IF
      END IF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck Reabld */
      SUBROUTINE REABLD(NBCSIM,NBESIM,NBPSIM,
     &     WORK,KFREE,LFREE)
C***********************************************************************
C
C     Read a load of B-vectors, limited by maximum number MAXSIM
C
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
      DIMENSION WORK(*)
C
      CALL QENTER('REABLD')
      KFIRST = KFREE
      NREST  = MAXSIM
C
C     Electronic type trial vectors
C     -----------------------------
C
      KBOE = KFIRST
      NESIM = MIN(NREST,NBESIM)
      IF(NESIM.GT.0) THEN
         NBUF   = NESIM*NZXOPEQ
         CALL MEMGET2('REAL','BOE',KBOE,NBUF,WORK,KFREE,LFREE)
         IREC = NERED
         JOFF = KBOE
         DO I = 1,NESIM
            IREC = IREC + 1
            CALL READ_DX(LUBOE,IREC,NZXOPEQ,WORK(JOFF))
            JOFF = JOFF + NZXOPEQ
         ENDDO
         NREST  = NREST - NESIM
         KFIRST = KFIRST + NZXOPEQ*NESIM
      ENDIF
C
C     Positronic type trial vectors
C     -----------------------------
C
      KBOP = KFIRST
      NPSIM = MIN(NREST,NBPSIM)
      IF(NPSIM.GT.0) THEN
         NBUF   = NPSIM*NZXOPPQ
         CALL MEMGET2('REAL','BOP',KBOP,NBUF,WORK,KFREE,LFREE)
         IREC = NPRED
         JOFF = KBOP
         DO I = 1,NPSIM
            IREC = IREC + 1
            CALL READ_DX(LUBOP,IREC,NZXOPPQ,WORK(JOFF))
            JOFF = JOFF + NZXOPPQ
         ENDDO
         NREST  = NREST - NPSIM
         KFIRST = KFIRST + NZXOPPQ*NPSIM
      ENDIF
C
C     Configurational type trial vectors
C     ----------------------------------
C
      KBCI = KFIRST
      NCSIM = MIN(NREST,NBCSIM)
      IF(NCSIM.GT.0) THEN
         NBUF  = NCSIM*NZCONFQ
         CALL MEMGET2('REAL','BCI',KBCI,NBUF,WORK,KFREE,LFREE)
         IREC = NCRED
         JOFF = KBCI
         DO I = 1,NCSIM
            IREC = IREC + 1
            CALL READ_DX(LUBCI,IREC,NZCONFQ,WORK(JOFF))
            JOFF = JOFF + NZCONFQ
         ENDDO
         NREST  = NREST - NCSIM
         KFIRST = KFIRST + NZCONFQ*NCSIM
      ENDIF
C
      CALL QEXIT('REABLD')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck wrsvec */
      SUBROUTINE WRSVEC(ICOUNT,SCVEC,SEVEC,SPVEC)
C***********************************************************************
C
C     Write sigma-vectors to file
C
C     Written by T.Saue Sep 18 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
      DIMENSION SCVEC(*),SEVEC(*),SPVEC(*)
C
C     Configurational part
C     --------------------
C
      IF(NZCONF.GT.0) THEN
         IREC = ICOUNT
         JOFF = 1
         DO ISIM = 1,NTSIM
            IREC = IREC + 1
!           CALL WRTDAC(LUSCI,NZCONFQ,SCVEC(JOFF),IREC + 1)
            CALL WRIT_DX(LUSCI,IREC+1,NZCONFQ,SCVEC(JOFF))
            JOFF = JOFF + NZCONFQ
         ENDDO
      ENDIF
C
C     Electronic part
C     ---------------
C
      IF(NZXOPE.GT.0) THEN
         IREC = ICOUNT
         JOFF = 1
         DO ISIM = 1,NTSIM
            IREC = IREC + 1
!           CALL WRTDAC(LUSOE,NZXOPEQ,SEVEC(JOFF),IREC)
            CALL WRIT_DX(LUSOE,IREC,NZXOPEQ,SEVEC(JOFF))
            JOFF = JOFF + NZXOPEQ
         ENDDO
      ENDIF
C
C     Positronic part
C     ---------------
C
      IF(NZXOPP.GT.0) THEN
         IREC = ICOUNT
         JOFF = 1
         DO ISIM = 1,NTSIM
            IREC = IREC + 1
!           CALL WRTDAC(LUSOP,NZXOPPQ,SPVEC(JOFF),IREC)
            CALL WRIT_DX(LUSOP,IREC,NZXOPPQ,SPVEC(JOFF))
            JOFF = JOFF + NZXOPPQ
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrsred */
      SUBROUTINE XRSRED(REDGP,REDE,REDS,GPCI,GPOE,GPOP,
     &                  IBCVC,IBEVC,IBPVC,IBTYP,BCVEC,BEVEC,BPVEC,
     &                  SCVEC,SEVEC,SPVEC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     PURPOSE:
C     Extend the reduced matrices with the NTSIM new trial vectors.
C
C     The reduced matrices are
C       Reduced Hessian:  RE[2] = B^D*E[2]*B     (D means dagger)
C       Reduced metric :  RS[2] = B^D*s[2]*B
C       E[2]B is already computed and given by the sigma-vector.
C     In the linear equation (LINEQ.EQ.TRUE) there is in addition the
C     reduced property gradient
C       Reduced gradient: RE[1] = B^D*E[1]
C
C     In the frequency-dependent case, trial b-vectors and sigma vectors
C     come in pairs
C
C        | Z  |    | Y  |
C        | Y* |  , | Z* |
C
C     Only the parts Z and Y are stored. Likewise the gradient has the
C     structure
C
C        | G* |
C        | G  |
C
C     Here only G* is stored.
C
C     In the frequency-independent case (STATIC) Z = Y, so the
C     dimension of the problem is reduced by a factor 2.
C
C   Input:
C     REDGP - reduced property gradient
C     REDE  - old reduced E[2]-matrix (dimension: NZRED-NTSIM)
C     REDS  - old reduced S[2]-matrix (dimension: NZRED-NTSIM)
C     BCVEC - new trial vectors corresponding to configurational
C             parameters
C     BEVEC - new trial vectors corresponding to electron-electron
C             rotations
C     BPVEC - new trial vectors corresponding to electron-positron
C             rotations
C     SCVEC - configurational part of new sigma-vectors
C     SEVEC - electron-electron part of new sigma vectors
C     SPVEC - electron-postrion part of new sigma vectors
C     S times the N b-vectors are created in XRSSLI and stored in
C             sigma-vectors
C     after extending the reduced Hessian
C
C   Output:
C     REDE, the new, extended reduced E[2]-matrix (dimension: NZRED)
C     REDS, rhe new, extended reduced S[2]-matrix (dimension: NZRED)
C
C   Scratch:
C
C     WORK, real scratch array dimension KZYVAR
C
C     Written by H.J.Aa.Jensen and T.Saue
C     Last revision Sep 19 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
c#include "ibndxdef.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "dcbxrs.h"
#include "dgroup.h"
      DIMENSION REDGP(NREDM),REDE(NREDM,NREDM),REDS(NREDM,NREDM),
     &          GPCI(*),GPOE(*),GPOP(*),
     &          IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),IBTYP(2,NREDM),
     &          BCVEC(NZCONFQ,*),BEVEC(NZXOPEQ,*),BPVEC(NZXOPPQ,*),
     &          SCVEC(NZCONFQ,*),SEVEC(NZXOPEQ,*),SPVEC(NZXOPPQ,*),
     &          WORK(*)
      LOGICAL   E2FULL
C
      CALL QENTER('XRSRED')
C
C     Exit routine if no new trial vectors
C     ====================================
C
      IF (NTSIM.LT.1) GO TO 9999
C
      E2FULL = .FALSE.
      IF (E2CHEK) E2FULL = .TRUE.
      IF (IPRXRS .GE. 15) THEN
         CALL HEADER('(XRSRED) Sigma vectors',-1)
         IF (NZCONF .GT. 1) THEN
         DO I = 1, NTSIM
            WRITE(LUPRI,'(/I5,A,I5)')
     &        I,'. conf. sigma vector out of',NTSIM
            CALL RPRCI(SCVEC(1,I),NZCONF,NZ,
     &           XRS_CIPROGRAM,THRZER,LUPRI)
         END DO
         END IF
         IF (NZXOPE .GT. 0) THEN
         WRITE(LUPRI,'(/I5,A)') NTSIM,' orbital e-e sigma vectors'
         CALL PRBVEC(LUPRI,SEVEC,NTSIM,NZXOPE)
         END IF
         IF (NZXOPP .GT. 0) THEN
         WRITE(LUPRI,'(/I5,A)') NTSIM,' orbital e-p sigma vectors'
         CALL PRBVEC(LUPRI,SPVEC,NTSIM,NZXOPP)
         END IF
      END IF
C
C     Memory check
C     ============
C
      KFRSAV = KFREE
C
C     ******************************************************************
C     ***** If linear equations (LINEQ), then extend reduced gradients *
C     ***** (we make the negative gradient)                            *
C     ******************************************************************
C
C     The gradient is a +-type vector and only gets contributions for
C     trial vector of (+)-type, that is for which IBTYP(2,*) equals 1
C
      IF(LINEQ)THEN
        CALL EXTRVC(IPX,REDGP,GPCI,GPOE,GPOP,BCVEC,BEVEC,BPVEC,
     &              IBTYP,IBCVC,IBEVC,IBPVC)
        CALL DSCAL(NTSIM,DM1,REDGP(IZRED+1),1)
        IF(IPRXRS.GE.4) THEN
          WRITE(LUPRI,'(/A)') '* XRSRED: Minus reduced gradient:'
          CALL PRIRVC(REDGP,NZRED,IBTYP)
        ENDIF
      ENDIF
C
C     **************************************
C     ****  Extend reduced E[2]-matrix  ****
C     **************************************
C
C     We need only make the upper triangle explicitly
C     since the matrices are real symmetric
C
C     The E[2}-matrix conserves the IH-parameter,
C     thus it only couples trial vectors of the same kind
C       (IH1 = IH2)
C
C     Allocate work space for previous trial vectors
      NBVEC = MAX(NZCONFQ,NZXOPEQ,NZXOPPQ)
      CALL MEMGET2('REAL','B1',KB1,NBVEC,WORK,KFREE,LFREE)
      CALL EXTRDE(REDE,E2FULL,IBTYP,IBCVC,IBEVC,IBPVC,
     &            BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,WORK(KB1))
C
C     Check for anti-Hermicity
C
      IF (E2FULL) THEN
         DAHS = D0
         DMAH = D0
         MAHI = 0
         MAHJ = 0
         DO I = 1, NZRED
            DO J = 1, NZRED
               AH = ABS(REDE(I,J) - REDE(J,I))
               IF ( AH .GE. DMAH ) THEN
                  DMAH = AH
                  MAHI = I
                  MAHJ = J
               END IF
               DAHS = DAHS + AH
            END DO
         END DO
         IF (DAHS .GT. 1.0D-14 .OR. IPRXRS.GE.4) THEN
            WRITE(LUPRI,'(/A,1P,D10.2/A,D10.2,5X,A,0P,I3,A,I3,A/)')
     &        '(XRSRED)  1-norm of anti-symmetry of REDE = ',DAHS,
     &        '(XRSRED) oo-norm of anti-symmetry of REDE = ',
     &        DMAH,' (element (',MAHI,',',MAHJ,'))'
         END IF
      END IF
C
      CALL MEMREL('EXTRE2',WORK,1,KB1,KFREE,LFREE)
C
C     Print section
C     =============
C
      IF(IPRXRS.GE.4) THEN
        WRITE(LUPRI,'(/A,I5)')
     &    'XRSRED: Reduced Hessian. Dimension:',NZRED
        IF (E2FULL) WRITE(LUPRI,'(A)')
     &    'XRSRED: E2FULL: E2(i,j) = s(i)*b(j)'
        CALL PRIRMT(REDE,NREDM,NZRED,NZRED,IBTYP)
      ENDIF
C
C     **************************************
C     ****  Extend reduced S[2]-matrix  ****
C     **************************************
C
C     We need only make the upper triangle explicitly
C     since the matrices are real symmetric
C
C     The S[2}-matrix switches the IH-parameter,
C     thus it only couples trial vectors of the opposite kind
C       (IH1 = - IH2)
C
C     In the static case, we only have the (+)-combination,
C     thus there can be no contribution from the
C     generalized metric.
C
C     Metric times trial vectors is stored in SVECS...
C
      IF(.NOT.STATIC) THEN
C
C       Generate tau-vectors: S[2]b
C       ---------------------------
C
        CALL XRSTVC(BEVEC,IBEVC,NESIM,IERED,
     &              BPVEC,IBPVC,NPSIM,IPRED,
     &              BCVEC,IBCVC,NCSIM,ICRED,
     &              SEVEC,SPVEC,SCVEC,NTSIM,IZRED,
     &              IBTYP,WORK,KFREE,LFREE)
C
C       Extend reduced S-matrix
C       -----------------------
C
        NBVEC = MAX(NZCONFQ,NZXOPEQ,NZXOPPQ)
        CALL MEMGET2('REAL','B1',KB1,NBVEC,WORK,KFREE,LFREE)
        CALL EXTRMT(-1,REDS,IBTYP,IBCVC,IBEVC,IBPVC,
     &              BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,WORK(KB1))
        CALL MEMREL('EXTRS2',WORK,1,KB1,KFREE,LFREE)
C
C       Print section
C       =============
C
        IF(IPRXRS.GE.4) THEN
          WRITE(LUPRI,'(/A,I5)')
     &      '* XRSRED: Reduced metric. Dimension:',NZRED
          CALL PRIRMT(REDS,NREDM,NZRED,NZRED,IBTYP)
        ENDIF
      ENDIF
C
C *** End of subroutine XRSRED
C
 9999 CONTINUE
      CALL QEXIT('XRSRED')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck xrsini */
      SUBROUTINE XRSINI(REDGP,REDS,REDE,GPCI,GPOE,GPOP,RCNV,EVECR,
     &                  RNORM,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Initialize XRS-module
C
C     Written by T.Saue Oct 15 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "consts.h"
C
#include "dcbibn.h"
C
#include "dcbfir.h"
#include "dcbxrs.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcborb.h"
      DIMENSION WORK(*),REDGP(*), REDS(*), REDE(*),
     &          GPCI(*),GPOE(*),GPOP(*),
     &          EVECR(NEVEC,NFREQ),RNORM(*),RCNV(*)
      logical old_dx
C
      CALL QENTER('XRSINI')
      KFRSAV = KFREE
C
C     File control
C     ============
C
      LURSP  = 30
      LUBCI  = 31
      LUBOE  = 32
      LUBOP  = 33
      LUSCI  = 34
      LUSOE  = 35
      LUSOP  = 36
      LUCYCL = 37
C
C     LURSP - file of restart info
C
      CALL OPNFIL(LURSP,'PAMRSP','UNKNOWN','PAMRSP')
C
C     LUB* - files of b-vectors
C     LUS* - files of sigma-vectors
C
      IF(NZCONF.GT.0) THEN
        call open_dx(lubci,'PAMBCI',nzconfq,'UNKNOWN',old_dx)
        call open_dx(lusci,'PAMSCI',nzconfq,'UNKNOWN',old_dx)
      ENDIF
      IF(NZXOPE.GT.0) THEN
        call open_dx(luboe,'PAMBOE',NZXOPEQ,'UNKNOWN',old_dx)
        call open_dx(lusoe,'PAMSOE',NZXOPEQ,'UNKNOWN',old_dx)
      ENDIF
      IF(NZXOPP.GT.0) THEN
        call open_dx(lubop,'PAMBOP',NZXOPPQ,'UNKNOWN',old_dx)
        call open_dx(lusop,'PAMSOP',NZXOPPQ,'UNKNOWN',old_dx)
      ENDIF
      OPEN(LUCYCL,FILE ='PAMCYC',STATUS='UNKNOWN',
     &            ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
C
C     Initialize
C     ==========
C
      IPX  =  JTIMOP
      IMX  = -JTIMOP
C svillaume: NTYP is defined but not used, and it looks like our NSTAT
C      NTYP = 2
C      IF(STATIC) NTYP = 1
      ITMIC = 0
      CALL DCOPY(NFREQ,DUMMY,0,RCNV,1)
      CALL DZERO(EVECR,NEVEC*NFREQ)
      CALL DZERO(REDE,N2REDM)
      CALL DZERO(REDS,N2REDM)
      IF (LINEQ) CALL DZERO(REDGP,NREDM)
      CALL DZERO(RNORM,NVPAR*NFREQ)
C
      NFOCKMAT = 0
      IF (NISHT .GT. 0) NFOCKMAT = NFOCKMAT + 1
      IF (NASHT .GT. 0) NFOCKMAT = NFOCKMAT + 1
      IF (NFOCKMAT .EQ. 0) THEN
         CALL QUIT('XRSINI: Neither inactive nor active orbitals!')
      END IF
C
C     Calculate GP norm for elimination of trial vectors with small
C     norm.
C
      IF (LINEQ) THEN
C svillaume: GPTNRM is already evaluated in PRPXL1 and is stored in a
C            common block (dbcxrs), so we don't need to do it now
C         GPTNRM = (D2*DNRM2(NZCONFQ,GPCI,1))**2
C     &          + (D2*DNRM2(NZXOPEQ,GPOE,1))**2
C     &          + (D2*DNRM2(NZXOPPQ,GPOP,1))**2
C         GPTNRM = SQRT(GPTNRM)
         GPTNRM = DNRM2(NZCONFQ,GPCI,1)**2 + DNRM2(NZXOPEQ,GPOE,1)**2 +
     &     DNRM2(NZXOPPQ,GPOP,1)**2
         GPTNRM = SQRT(GPTNRM)
         IF (IPRXRS.GE.3) WRITE(LUPRI,'(1X,A,1P,D12.5)')
     &     'Norm of gradient in XRSINI:',GPTNRM
      ELSE
C     ... excitation energies
         GPTNRM = D1
      END IF
C
C     Determine integral contribution
C     ===============================
C
      INTBUF = 0
      ITNEX  = ITMIC + 1
      CALL INTCON(INTFLG,INTBUF,INTDEF,
     &            DUMMY,CNVINT,ITNEX,ITRINT,XRS_INTTYP)
C
      IF (NZCONF .GT. 0) THEN
C
C       Write CI ref. vector to first record of PAMBCI
C       ----------------------------------------------
C
C       Note, that we have to renormalize CREF to 0.5, as all
C       other trial vectors.
C
         CALL MEMGET2('REAL','CREF',KCREF,NZCONFQ,WORK,KFREE,LFREE)
C
         CALL REAKRMC(LUKRMC,'CREF    ',WORK(KCREF),NZCONFQ)
         CALL DSCAL(NZCONFQ,DP5,WORK(KCREF),1)
!        stefan+hjaaj - fix for large-scale mcscf
!        CALL WRTDAC(LUBCI,NZCONFQ,WORK(KCREF),1)
         CALL WRIT_DX(LUBCI,1,NZCONFQ,WORK(KCREF))
         CALL MEMREL('XRSINI',WORK,1,KFRSAV,KFREE,LFREE)
C
C        Write CI gradient to first record of PAMSCI
C        -------------------------------------------
C
!        stefan+hjaaj - fix for large-scale mcscf
!        CALL WRTDAC(LUSCI,NZCONFQ,GPCI,1)
         CALL WRIT_DX(LUSCI,1,NZCONFQ,GPCI)
      END IF
C
      CALL QEXIT ('XRSINI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Extrmt */
      SUBROUTINE EXTRMT(IPAR,REDMAT,IBTYP,IBCVC,IBEVC,IBPVC,
     &                  BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,BVBUF)
C***********************************************************************
C
C     Extend reduced matrices in response equations
C     IPAR =  1	Hessian-type :
C               Coupling of trial vectors of same type (IH = JH)
C     IPAR = -1 Metric type:
C               Coupling of trial vectors of opposite type (IH = -JH)
C
C     Written by T.Saue Oct 15 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4=4.0D0,DM1=-1.0D0)
C
#include "dcbxrs.h"
C
      DIMENSION IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),IBTYP(2,NREDM),
     &          BCVEC(NZCONFQ,*),BEVEC(NZXOPEQ,*),BPVEC(NZXOPPQ,*),
     &          SCVEC(NZCONFQ,*),SEVEC(NZXOPEQ,*),SPVEC(NZXOPPQ,*),
     &          BVBUF(*),REDMAT(NREDM,NREDM)
      FFAC = D1
      IF(IMFREQ) FFAC = DM1
C
C     C-type trial vectors
C     ====================
      DO I = 1,NCSIM
CTROND If LOFFTY, remember to skip reference vector
CTROND    Nothing written yet.....
         call quit('argh ralle ralle: not implemented yet!')
      ENDDO
C
C     E-type trial vectors
C     ====================
C
C     Read in old vectors and multiply with all new sigma-vectors
C     -----------------------------------------------------------
C
      DO IREC = 1,IERED
        CALL READ_DX(LUBOE,IREC,NZXOPEQ,BVBUF)
        I  = IBEVC(IREC)
        IH = IBTYP(2,I)
        DO JN = 1,NTSIM
          J  = JN + IZRED
          JH = IBTYP(2,J)
          IF(IH*JH.EQ.IPAR) THEN
            REA = D4*DDOT(NZXOPEQ,BVBUF,1,SEVEC(1,JN),1)
            IF(JH .EQ. -1) THEN
               REDMAT(I,J) = FFAC * REA
               REDMAT(J,I) = REA
            ELSE
               REDMAT(I,J) = REA
               REDMAT(J,I) = FFAC * REA
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
C     For new vectors only upper triangular loop is needed
C     ----------------------------------------------------
C
      DO IN = 1,NESIM
        I  = IBEVC(IN+IERED)
        IH = IBTYP(2,I)
        DO J = I,NZRED
          JN = J - IZRED
          JH = IBTYP(2,J)
          IF(IH*JH.EQ.IPAR) THEN
            REA = D4*DDOT(NZXOPEQ,BEVEC(1,IN),1,SEVEC(1,JN),1)
            IF(JH .EQ. -1) THEN
               REDMAT(I,J) = FFAC * REA
               REDMAT(J,I) = REA
            ELSE
               REDMAT(I,J) = REA
               REDMAT(J,I) = FFAC * REA
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
C     P-type trial vectors
C     ====================
C
C
C     Read in old vectors and multiply with all new sigma-vectors
C     -----------------------------------------------------------
C
      DO IREC = 1,IPRED
        CALL READ_DX(LUBOP,IREC,NZXOPPQ,BVBUF)
        I  = IBPVC(IREC)
        IH = IBTYP(2,I)
        DO JN = 1,NTSIM
          J  = JN + IZRED
          JH = IBTYP(2,J)
          IF(IH*JH.EQ.IPAR) THEN
            REA = D4*DDOT(NZXOPPQ,BVBUF,1,SPVEC(1,JN),1)
            IF(JH .EQ. -1) THEN
               REDMAT(I,J) = FFAC * REA
               REDMAT(J,I) = REA
            ELSE
               REDMAT(I,J) = REA
               REDMAT(J,I) = FFAC * REA
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
C     For new vectors only upper triangular loop is needed
C     ----------------------------------------------------
C
      DO IN = 1,NPSIM
        I  = IBPVC(IN+IPRED)
        IH = IBTYP(2,I)
        DO J = I,NZRED
          JN = J - IZRED
          JH = IBTYP(2,J)
          IF(IH*JH.EQ.IPAR) THEN
            REA = D4*DDOT(NZXOPPQ,BPVEC(1,IN),1,SPVEC(1,JN),1)
            IF(JH .EQ. -1) THEN
               REDMAT(I,J) = FFAC * REA
               REDMAT(J,I) = REA
            ELSE
               REDMAT(I,J) = REA
               REDMAT(J,I) = FFAC * REA
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Extrde */
      SUBROUTINE EXTRDE(REDMAT,E2FULL,IBTYP,IBCVC,IBEVC,IBPVC,
     &                  BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,BVBUF)
C***********************************************************************
C
C     Extend reduced Hessian in response equations
C
C     Written by T.Saue and H.J.Aa.Jensen Sep 29 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4=4.0D0)
C
#include "dcbxrs.h"
C
      DIMENSION IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),IBTYP(2,NREDM),
     &          BCVEC(NZCONFQ,*),BEVEC(NZXOPEQ,*),BPVEC(NZXOPPQ,*),
     &          SCVEC(NZCONFQ,*),SEVEC(NZXOPEQ,*),SPVEC(NZXOPPQ,*),
     &          BVBUF(*),REDMAT(NREDM,NREDM)
      LOGICAL   E2FULL
C
      IPAR = 1
C
C     C-type trial vectors
C     ====================
C
C
      IF (NCSIM .GT. 0) THEN
C
         DO I = 1, IZRED
            IH = IBTYP(2,I)
            CALL READ_DX(LUSCI,I+1,NZCONFQ,BVBUF)
            DO JN = 1, NCSIM
               J = IBCVC(JN+ICRED)
               JH = IBTYP(2,J)
               IF ( IH * JH .EQ. IPAR ) THEN
                  REA = D4 * DDOT(NZCONFQ,BVBUF,1,BCVEC(1,JN),1)
                  REDMAT(I,J) = REA
                  REDMAT(J,I) = REA
               END IF
            END DO
         END DO
         IF (E2FULL) THEN
            DO IN = 1,ICRED
               I = IBCVC(IN)
               IH = IBTYP(2,I)
               CALL READ_DX(LUBCI,IN+1,NZCONFQ,BVBUF)
               DO JNC = 1, NCSIM
                  J = IBCVC(JNC+ICRED)
                  JH = IBTYP(2,J)
                  JN = J - IZRED
                  IF ( IH * JH .EQ. IPAR ) THEN
                     REA = D4 * DDOT(NZCONFQ,BVBUF,1,SCVEC(1,JN),1)
                     REDMAT(J,I) = REA
                  END IF
               END DO
            END DO
            DO IN = 1,IERED
               I = IBEVC(IN)
               IH = IBTYP(2,I)
               CALL READ_DX(LUBOE,IN,NZXOPEQ,BVBUF)
               DO JNC = 1, NCSIM
                  J = IBCVC(JNC+ICRED)
                  JH = IBTYP(2,J)
                  JN = J - IZRED
                  IF ( IH * JH .EQ. IPAR ) THEN
                     REA = D4 * DDOT(NZXOPEQ,BVBUF,1,SEVEC(1,JN),1)
                     REDMAT(J,I) = REA
                  END IF
               END DO
            END DO
            DO IN = 1,IPRED
               I = IBPVC(IN)
               IH = IBTYP(2,I)
               CALL READ_DX(LUBOP,IN,NZXOPPQ,BVBUF)
               DO JNC = 1, NCSIM
                  J = IBCVC(JNC+ICRED)
                  JH = IBTYP(2,J)
                  JN = J - IZRED
                  IF ( IH * JH .EQ. IPAR ) THEN
                     REA = D4 * DDOT(NZXOPPQ,BVBUF,1,SPVEC(1,JN),1)
                     REDMAT(J,I) = REA
                  END IF
               END DO
            END DO
         END IF
C
         DO IN = 1, NCSIM
            I = IBCVC(IN + ICRED)
            IH = IBTYP(2,I)
            IF (E2FULL) THEN
               J1 = IZRED+1
            ELSE
               J1 = I
            END IF
            DO J = J1,NZRED
               JN = J - IZRED
               JH = IBTYP(2,J)
               IF(IH*JH.EQ.IPAR) THEN
                  REA = D4 * DDOT(NZCONFQ,BCVEC(1,IN),1,SCVEC(1,JN),1)
                  REDMAT(J,I) = REA
                  IF (.NOT.E2FULL) REDMAT(I,J) = REA
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
C     E-type trial vectors
C     ====================
C
      IF(NESIM.GT.0) THEN
C
C        Read in all old sigma-vectors and multiply with new E-type
C        trial-vectors
C        ---------------------------------------------------------------
C
         DO I = 1,IZRED
            CALL READ_DX(LUSOE,I,NZXOPEQ,BVBUF)
            IH = IBTYP(2,I)
            DO JN = 1,NESIM
               J  = IBEVC(JN+IERED)
               JH = IBTYP(2,J)
               IF(IH*JH.EQ.IPAR) THEN
                  REA = D4*DDOT(NZXOPEQ,BVBUF,1,BEVEC(1,JN),1)
                  REDMAT(I,J) = REA
                  REDMAT(J,I) = REA
               ENDIF
            ENDDO
         ENDDO
         IF (E2FULL) THEN
            DO IN = 1,ICRED
               I = IBCVC(IN)
               IH = IBTYP(2,I)
               CALL READ_DX(LUBCI,IN+1,NZCONFQ,BVBUF)
               DO JNE = 1, NESIM
                  J = IBEVC(JNE+IERED)
                  JH = IBTYP(2,J)
                  JN = J - IZRED
                  IF ( IH * JH .EQ. IPAR ) THEN
                     REA = D4 * DDOT(NZCONFQ,BVBUF,1,SCVEC(1,JN),1)
                     REDMAT(J,I) = REA
                  END IF
               END DO
            END DO
            DO IN = 1,IERED
               I = IBEVC(IN)
               IH = IBTYP(2,I)
               CALL READ_DX(LUBOE,IN,NZXOPEQ,BVBUF)
               DO JNE = 1, NESIM
                  J = IBEVC(JNE+IERED)
                  JH = IBTYP(2,J)
                  JN = J - IZRED
                  IF ( IH * JH .EQ. IPAR ) THEN
                     REA = D4 * DDOT(NZXOPEQ,BVBUF,1,SEVEC(1,JN),1)
                     REDMAT(J,I) = REA
                  END IF
               END DO
            END DO
            DO IN = 1,IPRED
               I = IBPVC(IN)
               IH = IBTYP(2,I)
               CALL READ_DX(LUBOP,IN,NZXOPPQ,BVBUF)
               DO JNE = 1, NESIM
                  J = IBEVC(JNE+IERED)
                  JH = IBTYP(2,J)
                  JN = J - IZRED
                  IF ( IH * JH .EQ. IPAR ) THEN
                     REA = D4 * DDOT(NZXOPPQ,BVBUF,1,SPVEC(1,JN),1)
                     REDMAT(J,I) = REA
                  END IF
               END DO
            END DO
         END IF
C
C        For new vectors only upper triangular loop is needed
C        ----------------------------------------------------
C
         DO IN = 1,NESIM
            I  = IBEVC(IN+IERED)
            IH = IBTYP(2,I)
            IF (E2FULL) THEN
               J1 = IZRED+1
            ELSE
               J1 = I
            END IF
            DO J = J1,NZRED
               JN = J - IZRED
               JH = IBTYP(2,J)
               IF(IH*JH.EQ.IPAR) THEN
                  REA = D4*DDOT(NZXOPEQ,BEVEC(1,IN),1,SEVEC(1,JN),1)
                  REDMAT(J,I) = REA
                  IF (.NOT. E2FULL) REDMAT(I,J) = REA
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
C     P-type trial vectors
C     ====================
C
      IF(NPSIM.GT.0) THEN
C
C        Read in old sigma-vectors and multiply with all new
C        trial-vectors
C        ----------------------------------------------------
C
         DO I = 1,IZRED
            CALL READ_DX(LUSOP,I,NZXOPPQ,BVBUF)
            IH = IBTYP(2,I)
            DO JN = 1,NPSIM
               J  = IBPVC(JN+IPRED)
               JH = IBTYP(2,J)
               IF(IH*JH.EQ.IPAR) THEN
                  REA = D4*DDOT(NZXOPPQ,BVBUF,1,BPVEC(1,JN),1)
                  REDMAT(I,J) = REA
                  REDMAT(J,I) = REA
               ENDIF
            ENDDO
         ENDDO
         IF (E2FULL) THEN
            DO IN = 1,ICRED
               I = IBCVC(IN)
               IH = IBTYP(2,I)
               CALL READ_DX(LUBCI,IN+1,NZCONFQ,BVBUF)
               DO JNP = 1, NPSIM
                  J = IBPVC(JNP+IPRED)
                  JH = IBTYP(2,J)
                  JN = J - IZRED
                  IF ( IH * JH .EQ. IPAR ) THEN
                     REA = D4 * DDOT(NZCONFQ,BVBUF,1,SCVEC(1,JN),1)
                     REDMAT(J,I) = REA
                  END IF
               END DO
            END DO
            DO IN = 1,IERED
               I = IBEVC(IN)
               IH = IBTYP(2,I)
               CALL READ_DX(LUBOE,IN,NZXOPEQ,BVBUF)
               DO JNP = 1, NPSIM
                  J = IBPVC(JNP+IPRED)
                  JH = IBTYP(2,J)
                  JN = J - IZRED
                  IF ( IH * JH .EQ. IPAR ) THEN
                     REA = D4 * DDOT(NZXOPEQ,BVBUF,1,SEVEC(1,JN),1)
                     REDMAT(J,I) = REA
                  END IF
               END DO
            END DO
            DO IN = 1,IPRED
               I = IBPVC(IN)
               IH = IBTYP(2,I)
               CALL READ_DX(LUBOP,IN,NZXOPPQ,BVBUF)
               DO JNP = 1, NPSIM
                  J = IBPVC(JNP+IPRED)
                  JH = IBTYP(2,J)
                  JN = J - IZRED
                  IF ( IH * JH .EQ. IPAR ) THEN
                     REA = D4 * DDOT(NZXOPPQ,BVBUF,1,SPVEC(1,JN),1)
                     REDMAT(J,I) = REA
                  END IF
               END DO
            END DO
         END IF
C
C        For new vectors only upper triangular loop is needed
C        ----------------------------------------------------
C
         DO IN = 1,NPSIM
            I  = IBPVC(IN+IPRED)
            IH = IBTYP(2,I)
            IF (E2FULL) THEN
               J1 = IZRED+1
            ELSE
               J1 = I
            END IF
            DO J = J1,NZRED
               JN = J - IZRED
               JH = IBTYP(2,J)
               IF(IH*JH.EQ.IPAR) THEN
                  REA = D4*DDOT(NZXOPPQ,BPVEC(1,IN),1,SPVEC(1,JN),1)
                  REDMAT(J,I) = REA
                  IF (.NOT. E2FULL) REDMAT(I,J) = REA
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Extrvc */
      SUBROUTINE EXTRVC(IPAR,REDVC,FCVEC,FEVEC,FPVEC,
     &                  BCVEC,BEVEC,BPVEC,IBTYP,IBCVC,IBEVC,IBPVC)
C***********************************************************************
C
C     Extend reduced vectors in response equations
C     IPAR =  1	(+)-vector (gradient-type)
C     IPAR = -1 (-)-vector
C
C     Written by T.Saue Oct 15 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4=4.0D0)
C
#include "dcbxrs.h"
#include "dgroup.h"
C
      DIMENSION IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),IBTYP(2,NREDM),
     &          BCVEC(NZCONFQ,*),BEVEC(NZXOPEQ,*),BPVEC(NZXOPPQ,*),
     &          FCVEC(NZCONFQ,*),FEVEC(NZXOPEQ,*),FPVEC(NZXOPPQ,*),
     &          REDVC(NREDM)
C
C     C-type trial vectors
C     ====================
      DO I = 1,NCSIM
         IND = IBCVC(I + ICRED)
         IF ( IBTYP(2,IND) .EQ. IPAR ) THEN
            REDVC(IND) = D4 * DDOT(NZCONFQ,BCVEC(1,I),1,FCVEC,1)
         END IF
      ENDDO
C
C     E-type trial vectors
C     ====================
C
      DO I = 1,NESIM
        IND        = IBEVC(I+IERED)
        IF(IBTYP(2,IND).EQ.IPAR) THEN
          REDVC(IND) = D4*DDOT(NZXOPEQ,BEVEC(1,I),1,FEVEC,1)
        ENDIF
      ENDDO
C
C     P-type trial vectors
C     ====================
C
      DO I = 1,NPSIM
        IND        = IBPVC(I+IPRED)
        IF(IBTYP(2,IND).EQ.IPAR) THEN
          REDVC(IND) =  D4*DDOT(NZXOPPQ,BPVEC(1,I),1,FPVEC,1)
        ENDIF
      ENDDO
C
C     Print section
C     =============
C
      IF(IPRXRS.GE.200) THEN
        IREP = JSYMOP-1
        IF(NZXOPE.GT.0) THEN
          WRITE(LUPRI,'(A)')'EXTRVC : Full vector (e-e) part:'
          CALL PRQMAT(FEVEC,NZXOPE,1,NZXOPE,1,NZ,
     &                IPQTOQ(1,IREP),LUPRI)
        ENDIF
        IF(NZXOPP.GT.0) THEN
          WRITE(LUPRI,'(A)')'EXTRVC : Full vector (e-p) part'
          CALL PRQMAT(FPVEC,NZXOPP,1,NZXOPP,1,NZ,
     &                IPQTOQ(1,IREP),LUPRI)
        ENDIF
      ENDIF
      IF(IPRXRS.GE.50) THEN
        WRITE(LUPRI,'(A,I3)')
     &    'EXTRVC: Reduced vector Dim.: NZRED = ',NZRED
        CALL PRIRVC(REDVC,NZRED,IBTYP)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrssol */
      SUBROUTINE XRSSOL (REDGP,REDE,REDS,EVALR,EVECR,
     *                   RCNV,IBTYP,IBEVC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     This routine solves the response equation in reduced subspace
C     of dimension NZRED.
C
C     LINEQ = T:       (REDE - wREDS) X = REDG
C                      (linear system of equations)
C
C     LINEQ = F:       (REDE - wREDS) X = 0     (eigenvalue equation)
C
C     Input variables:
C       REDE(NREDM,NZRED) - reduced Hessian
C       REDS(NREDM,NZRED) - reduced metric
C       REDG (NFREQ)      - reduced property gradient (if LINEQ)
C       EVALR(NFREQ)      - frequencies
C
C     Output variables:
C       EVALR(NZRED)         - real part of eigenvalues
C                              (excitation energies)
C                              (if. NOT.LINEQ)
C       EVECR(NEVEC,NZRED,*) - solutions vectors in reduced space
C
C
C     Written by T.Saue and Hans Joergen Aa. Jensen. 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
C
      DIMENSION REDGP(NREDM),REDE(NREDM,NREDM),REDS(NREDM,NREDM),
     &          EVALR(NFREQ),EVECR(NEVEC,NFREQ),IBTYP(2,*),RCNV(*),
     &          IBEVC(*),WORK(*)
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , D2 = 2.0D0 , DP5 = 0.5D0)
      PARAMETER ( DP75 = 0.75D00 )
      PARAMETER ( DRTEST = 1.0D-5 )
C
      DIMENSION INERT(3)
C
      CALL QENTER('XRSSOL')
C
C     Check dimension of reduced problem
C     ==================================
C
      IF (NZRED.GT.NREDM) THEN
         WRITE(LUPRI,'(//A/A,I5,/A,I5)')
     *   ' >>> ERROR IN XRSSOL >>>',
     *   ' Dimension of reduced space is  ',NZRED,
     *   ' which exceeds allowed dimension',NREDM
         CALL QUIT('XRSSOL: Too large dimension of reduced space')
      END IF
C
C     Print section
C     =============
C
      IF(IPRXRS.GE.4) THEN
        WRITE(LUPRI,'(/A,I5)')
     &    'XRSSOL: Reduced Hessian. Dimension:',NZRED
        CALL PRIRMT(REDE,NREDM,NZRED,NZRED,IBTYP)
        IF(.NOT.STATIC) THEN
          WRITE(LUPRI,'(/A,I5)')
     &      'XRSSOL: Reduced metric. Dimension:',NZRED
          CALL PRIRMT(REDS,NREDM,NZRED,NZRED,IBTYP)
        ENDIF
      ENDIF
C
C     Check for instabilies in reduced Hessian
C
      CALL XRSSTAB(REDE,NREDM,NZRED,NPRED,IPRXRS,INERT,
     &             IBTYP,IBEVC,WORK,KFREE,LFREE)
C
      IF (FKRMC) THEN
         CALL KRMCSOL(NZRED,NREDM,REDE,REDS,EVALR,REDGP,
     &        EVECR(1,1),IBTYP,NPRED,IPRXRS,WORK,KFREE,LFREE)
      ELSE IF (LINEQ) THEN
C
C     Solve reduced linear response problem in subspace
C     =================================================
C
         DO IFREQ = 1,NFREQ
           AFREQ = EVALR(IFREQ)
           CALL XRSLINRED(NZRED,NREDM,REDE,REDS,AFREQ,REDGP,
     &          EVECR(1,IFREQ),IBTYP,INERT,NPRED,IPRXRS,
     &          WORK,KFREE,LFREE)
         ENDDO
C
C        Print section
C        -------------
C
         IF(IPRXRS.GE.1) THEN
           WRITE(LUPRI,'(/A,I4)') '<<B;B>> at iteration: ',ITMIC
           CALL PRSYMB(LUPRI,'-',44,0)
           DO IFREQ = 1,NFREQ
           VLRF = DP5*DDOT(NZRED,REDGP,1,EVECR(1,IFREQ),1)
           IF(RCNV(IFREQ).GT.THCXRS) THEN
             WRITE(LUPRI,'(3X,A,F15.8,5X,F15.8)')
     &         'Freq.:',EVALR(IFREQ),VLRF
           ELSE
              WRITE(LUPRI,'(3X,A,F15.8,5X,F15.8,3X,A)')
     &         'Freq.:',EVALR(IFREQ),VLRF,'(converged)'
           ENDIF
           ENDDO
        ENDIF
      ELSE
C     ... i.e. .NOT. LINEQ
C
C     Solve reduced general problem in subspace
C     =========================================
C
        CALL EIGRED(NFREQ,NZRED,NREDM,REDE,REDS,EVALR,EVECR,
     &          IERR,IPRXRS,WORK,KFREE,LFREE)
        IF (IPRXRS .GE. 1) THEN
          WRITE(LUPRI,'(/A,I4,A,I4/,(5F15.5))')
     &      'XRSSOL:',NFREQ,
     &      ' excitation frequencies (in a.u.) at iteration:',ITMIC,
     &      (EVALR(I),I = 1,NFREQ)
C         END IF
C         IF (IPRXRS .GE. 3) THEN
          WRITE(LUPRI,'(/A)')'XRSSOL: and the reduced solution vectors:'
          CALL PRIRMT (EVECR,NREDM,NZRED,NFREQ,IBTYP)
        ENDIF
      ENDIF
C
      CALL QEXIT('XRSSOL')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XRSLinred */
      SUBROUTINE XRSLINRED(N,LDA,AE,AS,W,G,X,IBTYP,INERT,
     &                  NPVEC,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Solve the reduced linear reponse equation
C
C       [AE - W*AS]X = G
C
C     Written by T.Saue Oct 16 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
#include "dcbxrs.h"
C
      DIMENSION AE(LDA,N),AS(LDA,N),G(N),X(*),
     &          IBTYP(2,*),INERT(3),DET(2),WORK(*)
C
      CALL QENTER('XRSLINRED')
      KFRSAV = KFREE
C
C     Memory allocation
C
      N2A = LDA*N
      IF(DAMPFREQ.NE.D0) THEN
        CALL MEMGET2('REAL','A2-cpp',KA2,2*N2A,WORK,KFREE,LFREE)
        CALL XRSLINRE2(N,LDA,N2A,AE,AS,WORK(KA2),W,G,X,IBTYP,
     &                  NPVEC,IPRINT,WORK,KFREE,LFREE)
      ELSE
        CALL MEMGET2('REAL','A2',KA2,N2A,WORK,KFREE,LFREE)
        CALL XRSLINRE1(N,LDA,N2A,AE,AS,WORK(KA2),W,G,X,IBTYP,
     &                  NPVEC,IPRINT,INERT,WORK,KFREE,LFREE)
      ENDIF
      CALL MEMREL('XRSLINRED',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('XRSLINRED')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XRSLinre1 */
      SUBROUTINE XRSLINRE2(N,LDA,N2A,AE,AS,A2,W,G,X,IBTYP,NPVEC,IPRINT,
     &                     WORK,KFREE,LFREE)
C***********************************************************************
C
C     Solve the reduced linear *complex* reponse equation
C
C       [AE - (W+i*DAMPFREQ)*AS]X = G
C
C     Written by S.Villaume & T.Saue 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D2 = 2.0D0, DP5 = 0.5D0)
#include "dcbxrs.h"
C
      DIMENSION AE(LDA,N),AS(LDA,N),A2(*),G(N),X(*),
     &          IBTYP(2,*),WORK(*)
C
      CALL QENTER('XRSLINRE2')
      KFRSAV = KFREE
C
C put the reduced gradient in X
      CALL DCOPY(N,G,1,X,2)
      CALL DCOPY(N,D0,0,X(2),2)
C build the reduced matrices in A2
      CALL DCOPY(N2A,AE,1,A2,2)
      CALL DCOPY(N2A,D0,0,A2(2),2)
      CALL DAXPY(N2A,-W,AS,1,A2,2)
      CALL DAXPY(N2A,-DAMPFREQ,AS,1,A2(2),2)
C solve A2*X=G
      CALL MEMGET2('INTE','PVT',KPVT,N,WORK,KFREE,LFREE)
      CALL ZSYSV('L',N,1,A2,LDA,WORK(KPVT),X,LDA,
     &           WORK(KFREE),LFREE,INFO)
CPP remove after debug
!     ALPHAR = DP5 * DDOT(N,G,1,X,2)
!     ALPHAI = DP5 * DDOT(N,G,1,X(2),2)
C
      CALL MEMREL('XRSLINRE2',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('XRSLINRE2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XRSLinre1 */
      SUBROUTINE XRSLINRE1(N,LDA,N2A,AE,AS,A2,W,G,X,IBTYP,
     &                  NPVEC,IPRINT,INERT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Solve the reduced linear *real or imaginary* reponse equation
C
C       [AE - W*AS]X = G
C
C     Written by T.Saue Dec 1996
C     Updated S.Villaume & T.Saue 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
#include "dcbxrs.h"
C
      DIMENSION AE(LDA,N),AS(LDA,N),A2(*),G(N),X(NEVEC),
     &          DET(2),INERT(3),IBTYP(2,*),WORK(*)
C
      CALL QENTER('XRSLINRE1')
      KFRSAV = KFREE
C
C     Form matrix [AE - W*AS]
C
      CALL DCOPY(N2A,AE,1,A2,1)
      IF (W.NE.D0) THEN
        CALL DAXPY(N2A,-W,AS,1,A2,1)
      END IF
C
C     Place G in buffer X
C
      CALL DCOPY(N,G,1,X,1)
      IF(IPRINT.GE.8) THEN
        WRITE(LUPRI,'(A)') 'XRSLINRED: Reduced E[2] - WS[2] matrix'
        CALL PRIRMT(A2,LDA,N,N,IBTYP)
        WRITE(LUPRI,'(A)') 'XRSLINRED: Reduced gradient'
        CALL PRIRVC(X,N,IBTYP)
      ENDIF
      IF(IMFREQ) THEN
        CALL MEMGET2('INTE','PVT',KPVT,N,WORK,KFREE,LFREE)
        CALL DGESOL(1,N,LDA,LDA,A2,X,WORK(KPVT),INFO)
      ELSE
        NNA = (N*(N+1))/2
        CALL MEMGET2('REAL','AN',KAN,NNA,WORK,KFREE,LFREE)
        CALL DAMATR(N,A2,LDA,WORK(KAN))
        CALL DSPSLI(N,1,WORK(KAN),X,A2,INFO,DET,INERT)
C
C       Print section
C
      ENDIF
      IF(IPRINT.GE.8) THEN
        WRITE(LUPRI,'(A,F12.5,A)')
     &      'XRSLINRED: Solution vector at frequency:',W,' a.u.'
        CALL PRIRVC(X,N,IBTYP)
      ENDIF

      CALL MEMREL('XRSLINRE1',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('XRSLINRE1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrsnex */
      SUBROUTINE XRSNEX(IBTYP,IBCVC,IBEVC,IBPVC,RCNV,EVALR,EVECR,
     &                  GPCI,GPOE,GPOP,RNORM,IVECS,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Purpose:
C       1) Construct residual:
C
C             LINEQ = T:  R(I) = (E[2] - W(I)*S[2])*X(I) + G
C             LINEQ = F:  R(I) = (E[2] - W(I)*S[2])*X(I)
C
C       2) Test for convergence of NEXCNV eigenvectors
C          Convergence criterium: |R(I)| .LE. THCXRS * |X(I)|
C
C       3) Generate new trial vectors
C
C     Written by T.Saue and H.J.Aa.Jensen 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),RCNV(*),
     &          EVALR(*),EVECR(NEVEC,*),GPCI(*),GPOE(*),GPOP(*),
     &          RNORM(NVPAR,NFREQ*NDAMP),IVECS(NFREQ,NVPAR),WORK(*)
      CALL QENTER('XRSNEX')
      KFRSAV = KFREE
      NESIM = 0
      NPSIM = 0
      NCSIM = 0
C
C       Memory allocation
C
      IFREQ = 1
      MFREQ = NFREQ
      NSIM  = MFREQ*NSTAT*NDAMP
      IF((NSIM+NZRED).GT.NREDM) THEN
        WRITE(LUPRI,'(A,I6)')
     &    '* Current number of trial vectors       :',NZRED,
     &    '* Requested number of new trial vectors :',NSIM,
     &    '* Maximum number of trial vectors NREDM :',NREDM
        CALL QUIT('* XRSNEX:Increase size of reduced space MAXRED !')
      ENDIF
      NBTOT = NSIM*(NZXOPEQ+NZXOPPQ+NZCONFQ)
      CALL MEMGET2('REAL','BTOT',KBTOT,NBTOT,WORK,KFREE,LFREE)
      KFIRST = KBTOT
      KBOE   = KFIRST
      KFIRST = KFIRST + NZXOPEQ*NSIM
      KBOP   = KFIRST
      KFIRST = KFIRST + NZXOPPQ*NSIM
      KBCI   = KFIRST
      KFIRST = KFIRST + NZCONFQ*NSIM
C
C     Allocate buffers
C
      KFRBUF = KFREE
      CALL MEMGET2('REAL','BEBUF',KBEBUF,NZXOPEQ*NDAMP*NSTAT,
     &   WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','BPBUF',KBPBUF,NZXOPPQ*NDAMP*NSTAT,
     &   WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','BCBUF',KBCBUF,NZCONFQ,
     &   WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','XNORM',KXNORM,NVPAR*MFREQ*NDAMP,
     &   WORK,KFREE,LFREE)
C
C     Generate new trial vectors
C     ==========================
C
      CALL XRSNE1(WORK(KBOE),WORK(KBOP),WORK(KBCI),
     &            MFREQ,EVALR(IFREQ),EVECR(1,IFREQ),
     &            GPOE,GPOP,GPCI,IBTYP,IBCVC,IBEVC,IBPVC,
     &            WORK(KBEBUF),WORK(KBPBUF),WORK(KBCBUF),RCNV,
     &            RNORM(1,IFREQ),WORK(KXNORM),IVECS,
     &            WORK,KFREE,LFREE)
C
C     Compress trial vectors in memory
C     ================================
C
      MBTOT = NBTOT
      NBTOT = NESIM*NZXOPEQ + NPSIM*NZXOPPQ + NCSIM*NZCONFQ
      IF(MBTOT.NE.NBTOT) THEN
C     ...compress vectors in memory
        KFIRST = KBTOT
C
        JBOE   = KBOE
        KBOE   = KFIRST
        NEDIM  = NZXOPEQ*NESIM
        KFIRST = KFIRST + NEDIM
        IF(KBOE.NE.JBOE.AND.NEDIM.GT.0)
     &   CALL DCOPY(NEDIM,WORK(JBOE),1,WORK(KBOE),1)
C
        JBOP   = KBOP
        KBOP   = KFIRST
        NPDIM  = NZXOPPQ*NPSIM
        KFIRST = KFIRST + NPDIM
        IF(KBOP.NE.JBOP.AND.NPDIM.GT.0)
     &    CALL DCOPY(NPDIM,WORK(JBOP),1,WORK(KBOP),1)
C
        JBCI   = KBCI
        KBCI   = KFIRST
        NCDIM  = NZCONFQ*NCSIM
        KFIRST = KFIRST + NCDIM
        IF(KBCI.NE.JBCI.AND.NCDIM.GT.0)
     &    CALL DCOPY(NCDIM,WORK(JBCI),1,WORK(KBCI),1)
C
        CALL MEMREL('XRSNEX.CMP',WORK,1,KFRSAV,KFREE,LFREE)
        CALL MEMGET2('REAL','BTOT',KBTOT,NBTOT,WORK,KFREE,LFREE)
      ELSE
        CALL MEMREL('XRSNEX.BUF',WORK,1,KFRBUF,KFREE,LFREE)
      ENDIF
      CALL QEXIT('XRSNEX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRBVEC */
      SUBROUTINE PRBVEC(IUNIT,BVEC,NBSIM,NBDIM)
C***********************************************************************
C     Print trial vector
C
C     Written by T.Saue Oct 24 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxrs.h"
      DIMENSION BVEC(NBDIM,NZ,NBSIM)
C
      IREP = JSYMOP - 1
      DO I = 1,NBSIM
        WRITE(LUPRI,'(A,I4)') '* Vector nr. ',I
        CALL PRQMAT(BVEC(1,1,I),NBDIM,1,NBDIM,1,NZ,
     &              IPQTOQ(1,IREP),IUNIT)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRBVCC */
      SUBROUTINE PRBVCC(IUNIT,BVEC,NBSTA,NBSIM,NBSTP,NBDIM)
C***********************************************************************
C     Print trial vector in cpp
C
C     Adapted from PRBVEC
C     Written by T.Saue & S.Villaume Jan2011
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxrs.h"
      DIMENSION BVEC(NBDIM,NZ,NBSIM)
C
      IREP = JSYMOP - 1
      ICOUNT = 1
      DO I = NBSTA,NBSIM,NBSTP
        WRITE(LUPRI,'(A,I4)') '* Vector nr. ',ICOUNT
        CALL PRQMAT(BVEC(1,1,I),NBDIM,1,NBDIM,1,NZ,
     &              IPQTOQ(1,IREP),IUNIT)
        ICOUNT = ICOUNT + 1
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck BVINFO */
      SUBROUTINE BVINFO(ISIM,IBTYP,IOFF,NBVEC)
C***********************************************************************
C
C     Give information on number of trial vectors
C
C     Written by T.Saue Nov 8 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dcbibn.h"
      DIMENSION IBTYP(2,*),ISIM(*)
C
C     Initialize
C
      CALL IZERO(ISIM,NVPAR)
C
C     Count vectors
C
      DO I = 1,NBVEC
        IBV = IBTYP(1,IOFF+I)
        IBH = IBTYP(2,IOFF+I)
        IF    (IBV.EQ.JBENDX) THEN
          IF(IBH.EQ.IPX) THEN
            ISIM(IPEP) = ISIM(IPEP) + 1
          ELSE
            ISIM(IPEM) = ISIM(IPEM) + 1
          ENDIF
        ELSEIF(IBV.EQ.JBPNDX) THEN
          IF(IBH.EQ.IPX) THEN
            ISIM(IPPP) = ISIM(IPPP) + 1
          ELSE
            ISIM(IPPM) = ISIM(IPPM) + 1
          ENDIF
        ELSEIF(IBV.EQ.JBCNDX) THEN
          IF(IBH.EQ.IPX) THEN
            ISIM(IPCP) = ISIM(IPCP) + 1
          ELSE
            ISIM(IPCM) = ISIM(IPCM) + 1
          ENDIF
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PRIRMT */
      SUBROUTINE PRIRMT (REDMAT,LDA,NRED,NCOL,IBTYP)
C***********************************************************************
C
C     Print reduced matrix with information
C
C     Written by T.Saue Nov 11 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0)
      INTEGER BEGIN
      DIMENSION REDMAT(LDA,NCOL),IBTYP(2,NRED)
      CHARACTER BTP(-1:2)*1,HTP(-1:1)*1
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA BTP/'C',' ','E','P'/,HTP/'-',' ','+'/
C
      NCTL = 1
C
C     First check if zero matrix
C
      AMAX = D0
      DO 10 J = 1,NCOL
        DO 10 I = 1,NRED
          AMAX = MAX(AMAX,ABS(REDMAT(I,J)))
   10 CONTINUE
      IF (AMAX .EQ. D0) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      KCOL = 4
C
      LAST = MIN(NCOL,KCOL)
      DO 2 BEGIN = 1,NCOL,KCOL
         WRITE (LUPRI,'(/18X,4(1X,A3,I4,1X,4A1,2X))')
     &       ('Col',I,'(',BTP(IBTYP(1,I)),HTP(IBTYP(2,I)),')',
     &       I = BEGIN,LAST)
         DO 1 K = 1,NRED
            DO 4 I = BEGIN,LAST
               IF (REDMAT(K,I).NE.D0) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LUPRI,'(4X,A3,I4,1X,4A1,2X,1P,8D15.7)')
     &        'Row',K,'(',BTP(IBTYP(1,K)),HTP(IBTYP(2,K)),')',
     &        (REDMAT(K,I), I = BEGIN,LAST)
    1    CONTINUE
         LAST = MIN(LAST+KCOL,NCOL)
    2 CONTINUE
    3 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PRIRVC */
      SUBROUTINE PRIRVC (REDVEC,NRED,IBTYP)
C***********************************************************************
C
C     Print reduced vector with information
C
C     Written by T.Saue Nov 11 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0)
      DIMENSION REDVEC(NRED),IBTYP(2,NRED)
      CHARACTER BTP(-1:2)*1,HTP(-1:1)*1,PFMT*26
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA BTP/'C',' ','E','P'/,HTP/'-',' ','+'/
C
C
C     First check if zero matrix
C
      AMAX = D0
      DO 10 I = 1,NRED
        AMAX = MAX(AMAX,ABS(REDVEC(I)))
   10 CONTINUE
      IF (AMAX .EQ. D0) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero vector.'
         GO TO 3
      END IF
C
C     Set print format
C
      IF (FFMIN.LE.AMAX .AND. AMAX.LE.FFMAX) THEN
C        use F output format
         PFMT = '(4X,A3,I4,1X,4A1,2X,F15.8)'
      ELSE
C        use 1PD output format
         PFMT = '(4X,A3,I4,1X,4A1,1P,D17.6)'
      END IF
C
      DO 2 I = 1,NRED
        IF(REDVEC(I).EQ.D0) GOTO 2
        WRITE (LUPRI,PFMT)
     &        'Row',I,'(',BTP(IBTYP(1,I)),HTP(IBTYP(2,I)),')',
     &        REDVEC(I)
    2 CONTINUE
    3 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XRSOUT */
      SUBROUTINE XRSOUT(RCNV,RNORM,ISIM,EVALR,IBTYP)
C***********************************************************************
C
C     Final output from XRSCTL
C       RCNV  - total norm of residual relative to norm of solution vector
C       RNORM - residual norm
C
C     Written by T.Saue Nov 19 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "dcbxrs.h"
C
      CHARACTER CPUTID*12,WORD*3
      DIMENSION RNORM(NVPAR,NFREQ*NDAMP),EVALR(NFREQ),ISIM(*),
     &          RCNV(NFREQ),IBTYP(2,NREDM)
C
      IF(UNCOUP) RETURN
      CALL HEADER('Final output from XRSCTL',-1)
C
C     Check for state of the calculation
C     ==================================
C
      IF (KCONV.EQ.-3) THEN
C
C     Start vectors rejected
C     ----------------------
C
        WRITE (LUPRI,'(/A/4X,A)')
     &      '*** WARNING - Start vectors rejected !',
     &      'No microiterations.'
        RETURN
      ELSE IF (KCONV.EQ.-2) THEN
C
C     Maximum dimension of reduced space exceeded
C     -------------------------------------------
C
        WRITE (LUPRI,'(/A,I3,A/4X,A)')
     &      '*** WARNING - Microiterations stopped after',
     &      ITMIC,' iterations',
     &      'because maximum dimension of reduced space exceeded.'
      ELSE IF (KCONV.EQ.-1) THEN
C
C     Linear dependence between new trial vectors
C     -------------------------------------------
C
        WRITE (LUPRI,'(/A,I3,A/4X,A)')
     &      '*** WARNING - Microiterations stopped after',
     &      ITMIC,' iterations',
     &      'because of linear dependence between new trial vectors'
      ELSE IF(KCONV.EQ.0)THEN
C
C     Convergence of iterative solution
C     ---------------------------------
C
        WRITE(LUPRI,'(/A,I3,A)')
     &     '*** Microiterations converged after',ITMIC,
     &     ' microiterations'
      ELSE
C
C           No convergence
C           --------------
C
        IF (ITMIC.GE.MAXITR) THEN
C
C         Max no of microiterations reached
C         ---------------------------------
C
          WRITE(LUPRI,'(/A,I4,A)')
     *      ' *** WARNING-Maximum number of microiterations,',
     *      ITMIC,', reached'
        END IF
      END IF
C
C     Give summary of calculations
C
      IF(MAXITR.GT.0) THEN
        IF(DAMPFREQ.NE.D0) THEN
          WRITE(LUPRI,'(/12X,A9,3X,A10,2X,12(2X,A2,A6))')
     &      'Freq.(au)','Total norm',
     &      (BVTYP(I),' norm ',I=1,NVPAR),
     &      (BVTYP(I),' norm ',I=1,NVPAR)
        ELSE
          WRITE(LUPRI,'(/12X,A9,3X,A10,2X,6(2X,A2,A6))')
     &      'Freq.(au)','Total norm',(BVTYP(I),' norm ',I=1,NVPAR)
        ENDIF
        NCOL = 36 + NDAMP*NVPAR*10
        REWIND LUCYCL
        DO ITR = 1,ITMIC
          CALL PRSYMB(LUPRI,'=',NCOL,0)
          READ(LUCYCL,'(I3,2X,A9,6(2X,A3,I3))')
     &       JTR,XRS_INTTYP,(WORD,ISIM(I),I=1,NVPAR)
          READ(LUCYCL,'(6D10.2)')
     &       ((RNORM(I,J),I=1,NVPAR),J=1,NFREQ*NDAMP)
          READ(LUCYCL,'(A12)') CPUTID
          WRITE(LUPRI,'(A5,I3,2X,A12,1X,A,A9,A,6(5X,A1,I3,A1))')
     &        '* It.',ITR,CPUTID,' (',XRS_INTTYP,')',
     &        ('[',ISIM(I),']',I=1,NVPAR)
          CALL PRSYMB(LUPRI,'-',NCOL,0)
          DO IFREQ = 1,NFREQ
            IF(DAMPFREQ.NE.D0) THEN
              INDR = 2*IFREQ-1
              INDI = 2*IFREQ
              RTNORMR = DNRM2(NVPAR,RNORM(1,INDR),1)
              RTNORMI = DNRM2(NVPAR,RNORM(1,INDI),1)
              RTNORM = DSQRT(RTNORMR**2 + RTNORMI**2)
              WRITE(LUPRI,'(10X,F10.6,A3,1P,D10.2,A3,12(D10.2))')
     &        EVALR(IFREQ),' : ',RTNORM,' | ',
     &        (RNORM(I,INDR),I=1,NVPAR),
     &        (RNORM(I,INDI),I=1,NVPAR)
            ELSE
              RTNORM = DNRM2(NVPAR,RNORM(1,IFREQ),1)
              WRITE(LUPRI,'(10X,F10.6,A3,1P,D10.2,A3,6(D10.2))')
     &        EVALR(IFREQ),' : ',RTNORM,' | ',(RNORM(I,IFREQ),I=1,NVPAR)
            ENDIF
          ENDDO
        ENDDO
        CALL PRSYMB(LUPRI,'=',NCOL,0)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck E2DIAG */
      SUBROUTINE E2DIAG(HEDIAG,HPDIAG,HCDIAG,SEDIAG,SPDIAG,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get diagonal elements of Hessian(approximate)
C
C     Written by T.Saue 1996
C
C***********************************************************************

      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
      DIMENSION HEDIAG(*),HPDIAG(*),HCDIAG(*),
     &          SEDIAG(*),SPDIAG(*),WORK(*)
      CALL QENTER('E2DIAG')
C
C     Orbital part
C
      IF((NZXOPE.GT.0).OR.(NZXOPP.GT.0)) THEN
        CALL E2ODIA(HEDIAG,HPDIAG,SEDIAG,SPDIAG,
     &              get_orbital_rotation_indices_pp(),
     &              get_orbital_rotation_indices_pn(),
     &              WORK,KFREE,LFREE)
      ENDIF
C
C     Configurational part
C
      IF(NZCONF.GT.0) THEN
         CALL E2CDIA(HCDIAG,WORK,KFREE,LFREE)
      ENDIF
C
      CALL QEXIT('E2DIAG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck e2cdia */
      SUBROUTINE E2CDIA(HCDIAG,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get configurational part of diagonal hessian elements.
C
C     Written by J. Thyssen - Nov 21 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dcbgen.h"
C
      DIMENSION HCDIAG(*)
      CALL QENTER('E2CDIA')
C
C
C     Read Hessian diagonal from file
C
      CALL REAKRMC(LUKRM3,'E2CDIAG ',HCDIAG,NZCONF)
C
      CALL QEXIT('E2CDIA')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck E2ODIA */
      SUBROUTINE E2ODIA(HEDIAG,HPDIAG,SEDIAG,SPDIAG,
     &                  JXOPE,JXOPP,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get orbital part of diagonal hessian elements
C
C     In RPA the diagonal elements are given by:
C
C       E[2}_{si,si} = F_{ss} - F_{ii} - [g_{ss,ii} - g_{si,is}]
C
C     which is approximated by
C
C       E[2}_{si,si} ~ F_{ss} - F_{ii}
C
C     Written by T.Saue Dec 1 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM4 = -4.0D0)
C
#include "dcbgen.h"
#include "dcbxrs.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbibn.h"
#include "dgroup.h"
C
      DIMENSION HEDIAG(*),HPDIAG(*),JXOPE(2,*),JXOPP(2,*),
     &          SEDIAG(*),SPDIAG(*),WORK(*)
      integer :: lfmo
C
      CALL QENTER('E2ODIA')

      KFRSAV = KFREE
      lfmo = max(NFOCKMAT*N2BBASXQ,2*N2ORBXQ)
      CALL MEMGET2('REAL','FMO',KFMO,lfmo,WORK,KFREE,LFREE)
      KFC = KFMO
      KFV = KFMO + N2ORBXQ
C
C     Get total Fock matrix in MO-basis
C
      IF (NASHT .GT. 0) THEN
         CALL MEMGET2('REAL','DV',KDV,N2ASHXQ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','FGENPV',KFGENPV,NORBT*NASHT*NZ,
     &      WORK,KFREE,LFREE)
         CALL GETMAT('FCMO+FVMO',WORK(KFMO),IPRXRS,FKRMC,
     $        WORK(KFREE),LFREE)
         CALL GETMAT('DVMO'     ,WORK(KDV) ,IPRXRS,FKRMC,
     $        WORK(KFREE),LFREE)
         CALL GETMAT('FGENPV',WORK(KFGENPV),IPRXRS,FKRMC,
     &        WORK(KFREE),LFREE)
      ELSE
         CALL MEMGET2('REAL','DV 0',KDV,0,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','FGENPV 0',KFGENPV,0,WORK,KFREE,LFREE)
         CALL GETMAT('FCMO',WORK(KFMO),IPRXRS,FKRMC,WORK(KFREE),LFREE)
      END IF
C
C     Electronic part
C
      IF(NZXOPE.GT.0) THEN
         IF(STERNC) THEN
C           Sternheimer approx.: diagonal of e-p Hessian
C           approximated with -2mc^2
C           An extra factor two (see E2ODI1 for details)
            FAC = -DM4*CVAL*CVAL
            CALL DCOPY(NZXOPE,FAC,0,HEDIAG,1)
         ELSE
           CALL E2ODI1(HEDIAG,WORK(KFC),WORK(KFV),
     &             WORK(KDV),WORK(KFGENPV),JXOPE,NZXOPE)
         ENDIF
         IF (FKRMC) THEN
            SEDIAG(1:NZXOPE) = 1.0D0
         ELSE
            CALL S2ODI1(SEDIAG,WORK(KDV),JXOPE,NZXOPE)
         END IF
         IF(IPRXRS.GE.5) THEN
            WRITE(LUPRI,'(/A)')
     &         '** E2ODIA: Diagonal Hessian and metric (e-e):'
            WRITE(LUPRI,'(2I5,1P,2D20.10)')
     &         (JXOPE(1,I),JXOPE(2,I),HEDIAG(I),SEDIAG(I), I = 1,NZXOPE)
         ENDIF
      ENDIF
C
C     Positronic part
C
      IF(NZXOPP.GT.0) THEN
         IF(STERNH) THEN
C           Sternheimer approx.: diagonal of e-p Hessian
C           approximated with -2mc^2
C           An extra factor two (see E2ODI1 for details)
            FAC = DM4*CVAL*CVAL
            CALL DCOPY(NZXOPP,FAC,0,HPDIAG,1)
         ELSE
            CALL E2ODI1(HPDIAG,WORK(KFC),WORK(KFV),
     &           WORK(KDV),WORK(KFGENPV),JXOPP,NZXOPP)
         ENDIF
         IF (FKRMC) THEN
            SPDIAG(1:NZXOPP) = 1.0D0
         ELSE
            CALL S2ODI1(SPDIAG,WORK(KDV),JXOPP,NZXOPP)
         END IF
         IF(IPRXRS.GE.5) THEN
            WRITE(LUPRI,'(/A)')
     &         '** E2ODIA: Diagonal Hessian and metric (e-p):'
            WRITE(LUPRI,'(2I5,1P,2D20.10)')
     &         (JXOPP(1,I),JXOPP(2,I),HPDIAG(I),SPDIAG(I), I = 1,NZXOPP)
         ENDIF
      ENDIF
C
C
      CALL MEMREL('E2ODIA',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('E2ODIA')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck E2ODI1 */
      SUBROUTINE E2ODI1(HODIAG,FC,FV,DV,FGENPV,JXOP,NXOP)
C***********************************************************************
C
C     Get orbital part of diagonal hessian elements
C
C     Hopefully this routine will be documented in my Ph.D. thesis.
C
C     H_{iu,iu} ~= D_{uu} FC_{ii} + 2 FD_{uu} - 2 FD_{ii} - F_{uu}
C
C     H_{ia,ia} ~= 2 FD_{aa} - 2 FD_{ii}
C
C     H_{uv,uv} ~= DV_{vv} FC_{uu} + DV_{uu} FC_{vv} - 2 D_{uv} FC_{uv}
C                  - F_{uu} - F_{vv}
C
C     H_{ua,ua} ~= D_{uu} FC_{aa} - F_{uu}
C
C     The neglected two-electron terms can be approximated by
C     D_{uu} FV_{pp}.
C
C     Written by J. Thyssen - Nov 22 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
      PARAMETER(D2 = 2.0D0, DP5 = 0.5D0)
C
#include "dcborb.h"
#include "dgroup.h"
#include "maxorb.h"
#include "dcbidx.h"
C
      DIMENSION HODIAG(NXOP)
      DIMENSION FC(NORBT,NORBT,NZ)
      DIMENSION FV(NORBT,NORBT,NZ)
      DIMENSION FGENPV(NORBT,NASHT,NZ)
      DIMENSION DV(NASHT,NASHT,NZ)
      DIMENSION JXOP(2,NXOP)
C
C
      CALL DZERO(HODIAG,NXOP)
      DO IXOP = 1,NXOP
         IS = JXOP(2,IXOP)
         II = JXOP(1,IXOP)
         ITYPS = IOBTYP(IS)
         ISYMS = IFSMO(IS)
         ITYPI = IOBTYP(II)
         ISYMI = IFSMO(II)
C
         IF (ITYPI .EQ. JTINAC) THEN
C
C           First index inactive:
C
C           - 2 FD_{ii} = - 2 FC_{ii} - 2 FV_{ii}
C
            HODIAG(IXOP) = - D2 * FC(II,II,1)
            IF (NASHT .GT. 0) THEN
               HODIAG(IXOP) = HODIAG(IXOP) - D2 * FV(II,II,1)
            END IF
C
C           + 2 FD_{pp} = 2 FC_{pp} + 2 FV_{pp}
C
            HODIAG(IXOP) = HODIAG(IXOP) + D2 * FC(IS,IS,1)
            IF (NASHT .GT. 0)
     &           HODIAG(IXOP) = HODIAG(IXOP) + D2 * FV(IS,IS,1)
C
            IF (ITYPS .EQ. JTACT) THEN
C
C              DV_{uu} FC_{ii} - F_{uu}
C
               ISACT = IDXG2U(IS)
               HODIAG(IXOP) = HODIAG(IXOP)
     &              + DV(ISACT,ISACT,1) * FC(II,II,1)
     &              - FGENPV(IS,ISACT,1)
C
C              Approximation to integrals
C
               HODIAG(IXOP) = HODIAG(IXOP) +
     &              DV(ISACT,ISACT,1) * FV(II,II,1)
C
            END IF
         ELSE IF (ITYPI .EQ. JTACT) THEN
C
C           First index active
C
C           DV_{uu} FC_{vv} - F_{uu}
C           + DV_{uu} FV_{vv} (approx. of integrals.)
C
            IIACT = IDXG2U(II)
C
            HODIAG(IXOP) = HODIAG(IXOP)
     &           + DV(IIACT,IIACT,1) * FC(IS,IS,1)
     &           + DV(IIACT,IIACT,1) * FV(IS,IS,1)
     &           - FGENPV(II,IIACT,1)
C
            IF (ITYPS .EQ. JTACT) THEN
C
C              DV_{vv} FC_{uu} - F_{vv}
c              + DV_vv} FV_{uu} (approx. of integrals.)
C              - 2 Re ( DV_{uv} FC_{uv} )
C              - 2 Re ( DV_{uv} FV_{uv} ) (approx. of integrals)
C
               ISACT = IDXG2U(IS)
C
               HODIAG(IXOP) = HODIAG(IXOP)
     &              + DV(ISACT,ISACT,1) * FC(II,II,1)
     &              + DV(ISACT,ISACT,1) * FV(II,II,1)
     &              - FGENPV(IS,ISACT,1)
     &              - D2 * DV(IIACT,ISACT,1) * FC(II,IS,1)
     &              - D2 * DV(IIACT,ISACT,1) * FV(II,IS,1)
C
               DO IZ = 2, NZ
                  HODIAG(IXOP) = HODIAG(IXOP)
     &                 + D2 * DV(IIACT,ISACT,IZ) * FC(II,IS,IZ)
     &                 + D2 * DV(IIACT,ISACT,IZ) * FV(II,IS,IZ)
               END DO
            END IF
         END IF
C
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck S2ODI1 */
      SUBROUTINE S2ODI1(SODIAG,DV,JXOP,NXOP)
C***********************************************************************
C
C     Get orbital part of diagonal metric elements
C
C     The diagonal elements are given by:
C
C       S[2]_{pq,pq} = <0| [Qp,Pq] |0>
C                    = Dpp - Dqq
C
C     Only real part is stored. What about MCSCF?
C
C
C     Written by J. Thyssen - Jun 1 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
      PARAMETER(D2 = 2.0D0)
C
#include "dcborb.h"
#include "dgroup.h"
#include "maxorb.h"
#include "dcbidx.h"
C
      DIMENSION SODIAG(NXOP),JXOP(2,NXOP),DV(NASHT,NASHT,NZ)
C
C
      DO IXOP = 1,NXOP
         IS = JXOP(2,IXOP)
         II = JXOP(1,IXOP)
         ITYPS = IOBTYP(IS)
         ISYMS = IFSMO(IS)
         ITYPI = IOBTYP(II)
         ISYMI = IFSMO(II)
C
C
         IF (ITYPI .EQ. JTINAC) THEN
            SODIAG(IXOP) = D2
         ELSE IF (ITYPI .EQ. JTACT) THEN
            IIACT = IDXG2U(II)
            SODIAG(IXOP) = DV(IIACT,IIACT,1)
         END IF
C
         IF (ITYPS .EQ. JTACT) THEN
            ISACT = IDXG2U(IS)
            SODIAG(IXOP) = SODIAG(IXOP) - DV(ISACT,ISACT,1)
         END IF
C
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck XRSST2 */
      SUBROUTINE XRSST2(TYP,BVEC,GP,HDIAG,SDIAG,
     &     EVALR,EVECR,IBTYP,IBVEC,IVEC,NBTYP,NBVEC,IBOFF,
     &     NBPREV,NDIM,NZL,PREC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate start vectors, called from XRSST1.
C
C     Written by T.Saue Dec 1 1996
C
C     TYP='E' : NDIM = NZXOPE, NZL = NZ
C     TYP='P' : NDIM = NZXOPP, NZL = NZ
C     TYP='C' : NDIM = NCONF,  NZL = NZ_in_CI
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, D2 = 2.0D0)
      PARAMETER (DLIM = 1.0D-6)
#include "dcbibn.h"
C
C Used from include files:
C  dcborb.h : NASHT
C  dcbxrs.h : several
C
#include "dcborb.h"
#include "dcbxrs.h"
C
      CHARACTER*1 TYP
C
      DIMENSION BVEC(NDIM,NZL,*),GP(NDIM,NZL),HDIAG(NDIM),SDIAG(NDIM),
     &          EVALR(NFREQ),IBTYP(2,*),IBVEC(*),IVEC(*),
     &          EVECR(NEVEC,NFREQ),PREC(NDIM,2,*),WORK(*)
      LOGICAL FNDRSP
#include "chrnos.h"
C
      CALL QENTER('XRSST2')
C
      FFAC = D1
      IF(IMFREQ) FFAC = DM1
      NBDIM = NDIM*NZL
      IF (NBTYP.EQ.JBENDX) THEN
         ITYP    = IPEP
         LUBVEC  = LUBOE
         IRECOFF = 0
      ELSE IF (NBTYP.EQ.JBPNDX) THEN
         ITYP    = IPPP
         LUBVEC  = LUBOP
         IRECOFF = 0
      ELSE IF (NBTYP .EQ. JBCNDX) THEN
         ITYP    = IPCP
         LUBVEC  = LUBCI
C        For CI record 1 is used for the CI ref. vector.
         IRECOFF = 1
      ELSE
         WRITE(LUPRI,'(A,A,I5)') 'XRSST2 ERROR: ',
     &        'unknown NBTYP = ',NBTYP
         CALL QUIT('XRSST2: Unknown NBTYP !')
      ENDIF
      IF(E2CHEK) THEN
        IB = 0
        DO IVAR = 1,NDIM
          DO IZ = 1,NZL
            IXX  = IPX
            DO IH = 1,2
              IB = IB + 1
              CALL DZERO(BVEC(1,1,IB),NBDIM)
              BVEC(IVAR,IZ,IB) = D1
              IBTYP(1,IBOFF+IB)     = NBTYP
              IBTYP(2,IBOFF+IB)     = IXX
              IBVEC(NBPREV+IB)      = IBOFF + IB
              IVEC(IB)              = 1
              IXX = IMX
            ENDDO
          ENDDO
        ENDDO
      ELSEIF(LINEQ) THEN
C
C     L I N E A R   R E S P O N S E   F U N C T I O N S
C     =================================================
C
C
        IF(RSREST) THEN
C
C       Restart on solution vectors; scan through file
C
          REWIND LURST
          NVEC = 0
          IXX  = IPX
          IROUND = 1
 10       CONTINUE
          IF(FNDRSP(LURST,BVEC(1,1,NVEC+1),RSPLAB,BVTYP(ITYP),
     &           '                ','  ',NBDIM)) THEN
            NVEC = NVEC + 1
            IBTYP(2,IBOFF+NVEC) = IXX
            IBTYP(1,IBOFF+NVEC) = NBTYP
            IBVEC(NBPREV+NVEC)  = IBOFF + NVEC
            IVEC(NVEC)          = 1
            IF(NVEC.EQ.NBVEC) GOTO 20
            GOTO 10
          ENDIF
          IF(.NOT.STATIC.AND.IROUND.EQ.1) THEN
            REWIND LURST
            ITYP = ITYP + 1
            IXX  = IMX
            IROUND = IROUND + 1
            GOTO 10
          ENDIF
          IF(NVEC.GT.0) THEN
            DO I = (NVEC+1),NBVEC
               IBTYP(1,IBOFF+I) = NBTYP
               IBTYP(2,IBOFF+I) = IXX
               IBVEC(NBPREV+I)  = IBOFF + I
               IVEC(I)          = 0
            ENDDO
            GOTO 20
          ELSE
            WRITE(LUPRI,'(A)')
     &           '*** WARNING: No restart vectors found !'
          ENDIF
        ENDIF
        IF(DIAGHE) THEN
          DO IP = 1,NFREQ*NDAMP,NDAMP
            IM = IP + NFREQ*NDAMP
            DO IDAMP = 0,NDAMP-1
              IBTYP(1,IBOFF+IP+IDAMP)     = NBTYP
              IBTYP(2,IBOFF+IP+IDAMP)     = IPX
              IBVEC(NBPREV+IP+IDAMP)      = IBOFF+IP+IDAMP
              IBTYP(1,IBOFF+IM+IDAMP)     = NBTYP
              IBTYP(2,IBOFF+IM+IDAMP)     = IMX
              IBVEC(NBPREV+IM+IDAMP)      = IBOFF+IM+IDAMP
              IF(UNCOUP) EVECR(IBOFF+IP,IP) = D1
            ENDDO
          ENDDO
C.........If the contribution to the norm of the property gradient
C         from this parameter type is below threshold, then no
C         trial vectors of this kind is generated
          GPTMP = DNRM2(NBDIM,GP,1)
          IF (IPRXRS.GE.3) WRITE(LUPRI,'(1X,A,1P,D12.5)')
     &         'Norm of gradient ratio in XRSST2:',GPTMP/GPTNRM
          IF (GPTMP.LT.0.3D0*GPTNRM.AND.(.NOT.UNCOUP)) THEN
            NMDIM = NFREQ*NSTAT*NDAMP
            CALL IZERO(IVEC,NMDIM)
            GOTO 20
          ENDIF
          IFREQ = 0
          DO IP = 1,NFREQ*NDAMP,NDAMP
            IM = IP + NFREQ*NDAMP
            IFREQ = IFREQ + 1
C
C           hjj+patrick:
C           If partial gradient (e-e, e-p, or conf) is small compared to the
C           total gradient we assume that it is unnecessary to provide a
C           corresponding trial vector now. If not then create a trial vector
C           from the gradient vector and the diagonal of the Hessian.
C
            FREQ = EVALR(IFREQ)
            IF(FREQ.EQ.D0) THEN
              DO IZ = 1,NZL
                DO IVAR = 1,NDIM
                  IF(ABS(HDIAG(IVAR)).GT.DLIM ) THEN
                    BVEC(IVAR,IZ,IP) = -GP(IVAR,IZ)/HDIAG(IVAR)
                  ELSE
                    BVEC(IVAR,IZ,IP) = -GP(IVAR,IZ)
                  END IF
                ENDDO
              ENDDO
C             tsaue:         retain norm for numerical stability !
              IF(.NOT.UNCOUP) THEN
                BNRM = GPTMP/DNRM2(NBDIM,BVEC(1,1,IP),1)
                CALL DSCAL(NBDIM,BNRM,BVEC(1,1,IP),1)
              ENDIF
              IVEC(IP) = 1
              IVEC(IM) = 0
            ELSE
              CALL GENPREC(PREC,HDIAG,SDIAG,FREQ,NDIM)
              DO ID = 0,NDAMP-1
                ID2 = ID+1
                DO IZ = 1,NZL
                  DO IVAR = 1,NDIM
                    BVEC(IVAR,IZ,IP+ID) = PREC(IVAR,1,ID2)*GP(IVAR,IZ)
                    BVEC(IVAR,IZ,IM+ID) = PREC(IVAR,2,ID2)*GP(IVAR,IZ)
                  ENDDO
                ENDDO
C             tsaue:         retain norm for numerical stability !
                IF(.NOT.UNCOUP) THEN
                  BNRM = GPTMP/DNRM2(NBDIM,BVEC(1,1,IP+ID),1)
                  CALL DSCAL(NBDIM,BNRM,BVEC(1,1,IP+ID),1)
                  BNRM = GPTMP/DNRM2(NBDIM,BVEC(1,1,IM+ID),1)
                  CALL DSCAL(NBDIM,BNRM,BVEC(1,1,IM+ID),1)
                ENDIF
                IVEC(IP+ID) = 1
                IVEC(IM+ID) = 1
              ENDDO
            END IF
          ENDDO
        ELSE IF (IXRSTB(NBTYP) .GT. 0) THEN
C
C         Restart on IXRSTB trial vectors read from trial vector file
C         (this is a useful interface; e.g., let the MCSCF program
C         write some trial vectors on the trial vector file, and let
C         the solver program start from that).
C
          DO I = 1, IXRSTB(NBTYP)
            NVEC = NVEC + 1
C           Assume it is a `+' type trial vector
            IBTYP(2,IBOFF+NVEC) = IPX
            IBTYP(1,IBOFF+NVEC) = NBTYP
            IBVEC(NBPREV+NVEC)  = IBOFF + NVEC
            IVEC(NVEC)          = 1
            CALL READ_DX(LUBVEC,I+IRECOFF,NBDIM,BVEC(1,1,NVEC))
          END DO
          WRITE(LUPRI,'(/,1X,A,I2,3A)')
     &         '(XRSST2)',IXRSTB(NBTYP), ' trial vectors of type ',
     &         BVTYP(ITYP),' read.'
C
        ELSE
C         start with GP vectors as start vectors
          IF(STATIC) THEN
            CALL DCOPY(NBDIM,GP,1,BVEC(1,1,1),1)
            IBTYP(1,IBOFF+1) = NBTYP
            IBTYP(2,IBOFF+1) = IPX
            IBVEC(NBPREV+1)  = IBOFF + 1
            IVEC(1)         = 1
          ELSE
            CALL DCOPY(NBDIM,GP,1,BVEC(1,1,1),1)
            IBTYP(1,IBOFF+1) = NBTYP
            IBTYP(2,IBOFF+1) = IPX
            IBVEC(NBPREV+1)  = IBOFF + 1
            IVEC(1)          = 1
            CALL DCOPY(NBDIM,GP,1,BVEC(1,1,2),1)
            IBTYP(1,IBOFF+2) = NBTYP
            IBTYP(2,IBOFF+2) = IMX
            IBVEC(NBPREV+2)  = IBOFF + 2
            IVEC(2)          = 1
          ENDIF
        ENDIF
 20     CONTINUE
      ELSE
C     ... .not. LINEQ (implemented by HJAaJ Nov.2000)
C
C     E X C I T A T I O N  E N E R G I E S etc.
C     ========================================
C
Ctsaue: divide also by NZL to get all possible singlets and triplets
         NEXST2 = NBVEC/2/NZL
         IF ( TYP .NE. 'C' .AND. NASHT .GT. 0) THEN
C        ... SDIAG(i) = 1 for CI-vector and = 2 for RPA
            DO IVAR = 1,NDIM
               HDIAG(IVAR) = HDIAG(IVAR) / SDIAG(IVAR)
            END DO
         END IF
C
C        Find lowest NEXST2 elements in HDIAG,
C        IVEC is temporarily used as pointer
C
         CALL FNDMIN(NEXST2,IVEC,HDIAG,NDIM,WORK(KFREE),LFREE)

         IF (IPRXRS.GE.3) WRITE(LUPRI,'(2A/,(2I10,F20.5))')
     &         ' Lowest HDIAG elements in XRSST2, type ',TYP,
     &         (I,IVEC(I),HDIAG(IVEC(I)),I=1,NEXST2)
C
Chj:     Here NEXST2 may be greater than NEXSTV from XRSST;
C        number of start trial vectors is increased with
C        the vectors which are nearly degenerate with
C        the last requested, i.e. number NEXSTV.
C

!Miro: Subscript #1 of the array HDIAG has value 0 which is less than the lower bound of 1 !
! here is naive fix:
         IF (IVEC(NEXSTV).GE.1) HDGSTV = HDIAG(IVEC(NEXSTV))
!
  210    IF (NEXST2 .GT. NEXSTV) THEN
            HDGDIF = HDIAG(IVEC(NEXST2)) - HDGSTV
            IF ( HDGDIF .GT. 0.01D0) THEN
               NEXST2 = NEXST2 - 1
               GO TO 210
            END IF
         END IF
         NEXST2Q=NEXST2*NZL
         NBVEC = 2*NEXST2Q
         CALL DZERO(BVEC,NBVEC*NBDIM)
C
Ctsaue:  One should run over all values of NZL to get all
C        singlets and triplets allowed in the selected symmetry
C        In C1 symmetry the spin states go as
C           1 - S0
C           2 - Tz
C           3 - Ty
C           4 - Tx
         IP = 0
         DO I = 1,NEXST2
           IVAR = IVEC(I)
           DO IZ = 1,NZL
             IP = IP + 1
             IM = IP + NEXST2Q
             IBTYP(1,IBOFF+IP)     = NBTYP
             IBTYP(2,IBOFF+IP)     = IPX
             IBVEC(NBPREV+IP)      = IBOFF + IP
             IBTYP(1,IBOFF+IM)     = NBTYP
             IBTYP(2,IBOFF+IM)     = IMX
             IBVEC(NBPREV+IM)      = IBOFF + IM
             BVEC(IVAR,IZ,IP) = D1
             BVEC(IVAR,IZ,IM) = D1
           ENDDO
         ENDDO
         CALL ICOPY(NBVEC,1,0,IVEC,1)
      END IF
C
C     Orthogonalize new b-vectors against all previous b-vectors
C     and among themselves, and, if requested, renormalize.
C
      IF(.NOT.UNCOUP) THEN
        CALL XRSORT(BVEC,NBTYP,NBVEC,NBPREV,
     &            IBVEC,IBTYP,IVEC,WORK,KFREE,LFREE)
      ENDIF
C
C
C     Save new vectors on file and update IBCVEC and IBTYP
C     ====================================================
C     (is is assumed that all new BVECS are consecutively
C     ordered on IBTYP)
C
      CALL XRSAVE(BVEC,NBTYP,NBDIM,NBVEC,NBPREV,
     &            IBVEC,IBTYP,IVEC,WORK,KFREE,LFREE)
C
      CALL QEXIT('XRSST2')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck GENPREC */
      SUBROUTINE GENPREC(PREC,HDIAG,SDIAG,FREQ,NVAR)
C***********************************************************************
C
C     Purpose:
C     compute the preconditioning element P, A, B, C, D for the
C     complex response solver
C     see Eqs.(24),(25) and (26) of J. Chem. Phys. 133(2010) 064105
C
C     Input
C     PREC: were to store the preconditioner element
C     HDIAG: diagonal of the hessian
C     SDIAG: diagonal of the metric
C     FREQ: omega, imaginary if IMFREQ
C     NVAR: total number of element in each vector
C
C     written by S.Villaume & T.Saue Nov-Dec 10
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,DM1 = -1.0D0, D2 = 2.0D0, D4 = 4.0D0)
#include "dcbibn.h"
C
#include "dcbxrs.h"
#include "dgroup.h"
C
      DIMENSION PREC(NVAR,2,*),HDIAG(NVAR),SDIAG(NVAR)
C
      CALL QENTER('GENPREC')
CTROND: this one probably be removed
      CALL DZERO(PREC,2*NVAR*NDAMP)

      SIGNFAC = D1
      IF(IMFREQ) SIGNFAC=DM1
C.....Eq.(24)
      IF(DAMPFREQ.NE.D0) THEN
        FREQ2 = FREQ * FREQ
        DAMP  = DAMPFREQ
        DAMP2 = DAMPFREQ * DAMPFREQ
C  compute 4 * omega**2 * gamma**2
        FREQ2TDAMP2T4 = FREQ2 * DAMP2 * D4
C  compute 2 * omega * gamma
        FREQTDAMPT2 = FREQ * DAMP * D2
C  compute omega**2 - gamma**2 and omega**2 + gamma**2
        FREQ2MDAMP2 = FREQ2 - DAMP2
        FREQ2PDAMP2 = FREQ2 + DAMP2
        DO IVAR=1,NVAR
C  compute few common expressions
          H2DIAG = HDIAG(IVAR) * HDIAG(IVAR)
          S2DIAG = SDIAG(IVAR) * SDIAG(IVAR)
          FAC1 = (H2DIAG - (FREQ2MDAMP2 * S2DIAG))
          FAC2 = FREQ2TDAMP2T4 * S2DIAG * S2DIAG
          FAC3 = FREQ2PDAMP2 * S2DIAG
C.........denominator
          FAC = DM1/((FAC1 * FAC1) + FAC2)
C.........factor A/P
          PREC(IVAR,1,1) = HDIAG(IVAR) * FAC1 * FAC
C.........factor B/P
          PREC(IVAR,2,1) = FREQ * SDIAG(IVAR) * (H2DIAG - FAC3) * FAC
C.........factor D/P
          PREC(IVAR,1,2) = FREQTDAMPT2 * HDIAG(IVAR) * S2DIAG * FAC
C.........factor C/P
          PREC(IVAR,2,2) = DAMP * SDIAG(IVAR) * (H2DIAG + FAC3) * FAC
        END DO
C.....Eq.(25,26)
      ELSE
        FREQ2 = FREQ * FREQ
        DO IVAR=1,NVAR
C  compute few common expressions
          H2DIAG = HDIAG(IVAR) * HDIAG(IVAR)
          S2DIAG = SDIAG(IVAR) * SDIAG(IVAR)
C.......denominator
          FAC = DM1/(H2DIAG - (SIGNFAC * FREQ2 * S2DIAG))
C.......factor A/P
          PREC(IVAR,1,1) = HDIAG(IVAR) * FAC
C.......factor B/P
          PREC(IVAR,2,1) = FREQ * SDIAG(IVAR) * FAC
        END DO
      ENDIF
C
      CALL QEXIT('GENPREC')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fndxlb */
      LOGICAL FUNCTION FNDXLB(RTNLBL,LU)
C
C  5-Aug-1986 hjaaj
C  Feb 1997 - tsaue
C
      CHARACTER*8 RTNLBL(3), B(4), STAR8
      PARAMETER (STAR8 = '********')
Ctsaue      REWIND LU
      IRDERR = 0
    1 READ (LU,END=3,ERR=6) B
      IRDERR = 0
      IF (B(1).NE.STAR8)     GO TO 1
      IF (B(4).NE.RTNLBL(3)) GO TO 1
      IF (B(3)(1:1).NE.RTNLBL(2)(1:1)) GO TO 1
      FNDXLB    = .TRUE.
      RTNLBL(1) = B(2)
      RTNLBL(2) = B(3)
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      GO TO 8
    3 CONTINUE
#if defined (SYS_CRAY) || defined (VAR_MFDS) || defined (SYS_T3D)
C 880916-hjaaj -- attempt to fix an IBM problem
C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
C     of FTxxF001), backspace makes LU ready for append.
C 940510-hjaaj: same change for Cray's multifile datasets
      BACKSPACE LU
#endif
    8 FNDXLB = .FALSE.
C
   10 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Wrtrsp*/
      SUBROUTINE WRTRSP(LU,RSPVEC,JXOP,
     &                  LAB1,TYP1,FREQ1,JSYMOP1,JTIMOP1,
     &                  LAB2,TYP2,FREQ2,JSYMOP2,JTIMOP2,
     &                     RNORM,INTFLG,LEN)
C***********************************************************************
C
C   Write response vectors
C   Written by T.Saue and P.Norman Aug 23 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      LOGICAL FNDLAB16
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2
      DIMENSION RSPVEC(LEN),JXOP(2,LEN)
      REWIND LU
      IF(FNDLAB16('END_OF_THIS_FILE',LU)) THEN
        BACKSPACE LU
        ERGRSP = 0
CHJ TODO: note that ERGRSP was not defined, what is the idea ????
        WRITE(LU) LAB1,LAB2,TYP1,TYP2,
     &          FREQ1,FREQ2,JSYMOP1,JSYMOP2,
     &          JTIMOP1,JTIMOP2,RNORM,LEN,INTFLG,
     &          ERGRSP,NTBAS(0),NORBT
        CALL WRITT(LU,LEN,RSPVEC)
!radovan: this was previously LEN but this goes out of bounds!
!         the length of the orbital rotation index vector is 2*LEN/NZ
!       CALL WRITI(LU,2*LEN,JXOP)
        CALL WRITI(LU,2*LEN/NZ,JXOP)
        WRITE(LU) 'END_OF_THIS_FILE'
      ELSE
        WRITE(LUPRI,'(A)') 'WRTRSP: EOF flag not found !'
        CALL QUIT('WRTRSP: EOF flag not found !')
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Rearsp */
      LOGICAL FUNCTION REARSP(LU,RSPVEC,
     &                  LAB1,TYP1,FREQ1,JSYMOP1,JTIMOP1,
     &                  LAB2,TYP2,FREQ2,JSYMOP2,JTIMOP2,
     &                     RNORM,LEN,INTFLG)
C***********************************************************************
C
C   Read response vectors
C   Written by T.Saue and P.Norman Aug 23 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
      CHARACTER*16 LAB1,LAB2,XLAB1,XLAB2
      CHARACTER*2 TYP1,TYP2,XTYP1,XTYP2
      DIMENSION RSPVEC(LEN)
C
      REWIND LU
 5    CONTINUE
      READ(LU,ERR=10,END=20) XLAB1
      IF(XLAB1.EQ.'END_OF_THIS_FILE') THEN
        REARSP = .FALSE.
        RETURN
      ELSE
        BACKSPACE LU
        READ(LU,ERR=10,END=20) XLAB1,XLAB2,XTYP1,XTYP2,
     &          XFREQ1,XFREQ2,ISYMOP1,ISYMOP2,
     &          ITIMOP1,ITIMOP2,RNORM,LEN,JNTFLG,ERGRSP,NB,NO
        IF(XLAB1.EQ.LAB1.AND.XLAB2.EQ.LAB2.AND.
     &   XFREQ1.EQ.FREQ1.AND.XFREQ2.EQ.FREQ2.AND.
     &   XTYP1.EQ.TYP1.AND.XTYP2.EQ.TYP2.AND.
     &   ISYMOP1.EQ.JSYMOP1.AND.ISYMOP2.EQ.JSYMOP2.AND.
     &   ITIMOP1.EQ.JTIMOP2.AND.ITIMOP2.EQ.JTIMOP2.AND.
     &   INTFLG.EQ.JNTFLG) THEN
          CALL READT(LU,LEN,RSPVEC)
          REARSP = .TRUE.
        ELSE
          READ(LU)
          GOTO 5
        ENDIF
      ENDIF
      RETURN
 10   CONTINUE
      IF(XLAB1.EQ.'END_OF_THIS_FILE') THEN
        REARSP = .FALSE.
        RETURN
      ELSE
        WRITE(LUPRI,'(A,I2)')
     &   '* REARSP: Error reading unit:',LU
        CALL QUIT('* REARSP: Error reading rspvec')
      ENDIF
 20   CONTINUE
      WRITE(LUPRI,'(A,I2)')
     &   '* REARSP: EOF reading unit:',LU
      CALL QUIT('* REARSP: EOF reading rspvec')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck fndrsp */
      LOGICAL FUNCTION FNDRSP(LU,RSPVEC,LAB1,TYP1,LAB2,TYP2,LEN)
C***********************************************************************
C
C     Scan file for response vectors
C     Written by T.Saue and P.Norman Aug 23 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "pgroup.h"
      CHARACTER*16 LAB1,LAB2,XLAB1,XLAB2
      CHARACTER*2 TYP1,TYP2,XTYP1,XTYP2
      DIMENSION RSPVEC(LEN)
C
 5    CONTINUE
      READ(LU,ERR=10,END=20) XLAB1
      IF(XLAB1.EQ.'END_OF_THIS_FILE') THEN
         FNDRSP = .FALSE.
         RETURN
      ELSE
         BACKSPACE LU
         READ(LU,ERR=10,END=20) XLAB1,XLAB2,XTYP1,XTYP2,
     &        XFREQ1,XFREQ2,ISYMOP1,ISYMOP2,
     &        ITIMOP1,ITIMOP2,RNORM,LEN,INTFLG,ERGRSP,NB,NO
         IF(XLAB1.EQ.LAB1.AND.XLAB2.EQ.LAB2.AND.
     &        XTYP1.EQ.TYP1.AND.XTYP2.EQ.TYP2) THEN
            CALL READT(LU,LEN,RSPVEC)
C...........skip orbital rotations
            READ(LU)
            WRITE(LUPRI,'(A)')
     &           '*** Restart on response vector:'
            WRITE(LUPRI,'(A,A16,2X,A2,A,A3,A,I2,A,F7.4)')
     &           'A: ',LAB1,TYP1,' boson irrep: ',REP(ISYMOP1-1),
     &           ' timrev:',ITIMOP1,' frequency: ',XFREQ1
            IF(LAB2.NE.'                ') THEN
               WRITE(LUPRI,'(A,A16,2X,A2,A,A3,A,I2,A,F7.4)')
     &              'B: ',LAB2,TYP2,' boson irrep: ',REP(ISYMOP2-1),
     &              ' timrev:',ITIMOP2,' frequency: ',XFREQ2
            ENDIF
            FNDRSP = .TRUE.
            RETURN
         ELSE
            READ(LU)
            READ(LU)
            GOTO 5
         ENDIF
      ENDIF
      RETURN
 10   CONTINUE
      WRITE(LUPRI,'(A,I2)')
     &     '* FNDRSP: Error reading unit:',LU
      CALL QUIT('* FNDRSP: Error reading rspvec')
 20   CONTINUE
      WRITE(LUPRI,'(A,I2)')
     &     '* FNDRSP: EOF reading unit:',LU
      CALL QUIT('* FNDRSP: EOF reading rspvec')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fndlab */
      LOGICAL FUNCTION FNDLAB16(A,LU)
C
C     Based on FNDLAB
C     T.Saue & P.Norman Aug 23 2000
C
      CHARACTER*16 A,B
      IRDERR = 0
    1 READ(LU,END=3,ERR=6)B
      IRDERR = 0
      IF (B.NE.A) GO TO 1
      FNDLAB16 = .TRUE.
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      GO TO 8
    3 CONTINUE
#if defined (SYS_CRAY) || defined (VAR_MFDS) || defined (SYS_T3D)
C 880916-hjaaj -- attempt to fix an IBM problem
C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
C     of FTxxF001), backspace makes LU ready for append.
C 940510-hjaaj: same change for Cray's multifile datasets
      BACKSPACE LU
#endif
C
    8 FNDLAB16 = .FALSE.
C
   10 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Xrsave */
      SUBROUTINE XRSAVE(BVEC,NBTYP,NBDIM,NBVEC,NBPREV,
     &                  IBVEC,IBTYP,IVEC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Save trial vectors to file and update IBVEC and IBTYP
C
C     Written by T. Saue Jan 2 2006
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbibn.h"
#include "dcbxrs.h"
      DIMENSION BVEC(NBDIM,*),IBTYP(2,*),IBVEC(*),IVEC(*),WORK(*)
C
      IF (NBTYP.EQ.JBENDX) THEN
         LUBVEC  = LUBOE
         IRECOFF = 0
      ELSE IF (NBTYP.EQ.JBPNDX) THEN
         LUBVEC  = LUBOP
         IRECOFF = 0
      ELSE IF (NBTYP .EQ. JBCNDX) THEN
         LUBVEC  = LUBCI
C
C     For CI record 1 is used for the CI ref. vector.
C
         IRECOFF = 1
      ELSE
         WRITE(LUPRI,'(A,A,I5)') 'XRSAVE ERROR: ',
     &        'unknown NBTYP = ',NBTYP
         CALL QUIT('XRSST2: Unknown NBTYP !')
      ENDIF
C
      IREC = NBPREV
C
      DO 10 IBNB = 1,NBVEC
         IB   = IBNB + NBPREV
         INDB = IBVEC(IB)
         IF(IVEC(IBNB).EQ.0) GOTO 10
         IREC = IREC + 1
!        CALL WRTDAC(LUBVEC,NBDIM,BVEC(1,IBNB),IREC+IRECOFF)
         CALL WRIT_DX(LUBVEC,IREC+IRECOFF,NBDIM,BVEC(1,IBNB))
         IF(IB.NE.IREC) THEN
            IND          = IBVEC(IREC)
            IBTYP(1,IND) = IBTYP(1,INDB)
            IBTYP(2,IND) = IBTYP(2,INDB)
            IBN          = IREC - NBPREV
            CALL DCOPY(NBDIM,BVEC(1,IBNB),1,BVEC(1,IBN),1)
         ENDIF
 10   CONTINUE
      NBVEC = IREC - NBPREV
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck setrsp */
      SUBROUTINE SETRSP
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Set default values for a number of control variables of the response
C     module such that they need to be modified only in modules where
C     they are used.
C
C***********************************************************************
C
      implicit none
#include "dcbxrs.h"
      FKRMC  = .FALSE.
      UNCOUP = .FALSE.
      IMFREQ = .FALSE.
      E2CHEK = .FALSE.
      STERNC = .FALSE.
      DAMPFREQ = 0.0D0
      NEVEC  = NREDM
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck xrsstab */
      SUBROUTINE XRSSTAB(A,LDA,N,NEPROT,IPRINT,INERT,IBTYP,IBEVC,
     &                   WORK,KFREE,LFREE)
C***********************************************************************
C
C     Check for instability by looking for negative eigenvalues of
C     the reduced Hessian.
C
C        A - reduced Hessian of dimension N
C
C        INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.
C        INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.
C        INERT(3)  =  NUMBER OF ZERO EIGENVALUES.
C
C     In the general relativistic case the ground state corresponds to
C     a saddle point and so one should subtract the number of e-p
C     rotations given by NEPROT
C
C     Written by T. Saue - updated Dec 9 2010
C
C***********************************************************************
      use orbital_rotation_indices
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
      DIMENSION A(LDA,N),IBTYP(*),IBEVC(*),WORK(*)
      DIMENSION INERT(3),DET(2)
C
      KFRSAV=KFREE
      NNA = (N*(N+1))/2
      CALL MEMGET2('INTE','PVT',KPVT,N,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','AN', KAN,NNA,WORK,KFREE,LFREE)
      CALL DAMATR(N,A,LDA,WORK(KAN))                                ! Row-pack lower triangle of matrix AMAT into ATR
      CALL DSPFA (WORK(KAN),N,WORK(KPVT),INFO)                      ! Factorization: A = U*D*TRANS(U)
      CALL DSPDI(WORK(KAN),N,WORK(KPVT),DET,INERT,WORK(KFREE),100)  ! Calculate the inertia
      IPAR = INERT(2) - NEPROT                                      
!     stefan: reduce output for mcscf - FIXME: find a better solution to handle the quantity of output...
      if(nasht .le. 0)then
        IF(IPRINT.GE.2 .OR. IPAR.GT.0 .OR. INERT(3).ne.0) THEN
          WRITE(LUPRI,'(3X,A,4(A,I5))')
     &        '* Reduced Hessian: ',
     &        'Number of positive/zero/negative eigenvalues : (',
     &        INERT(1),'/',INERT(3),'/',INERT(2),
     &        ');   +- (ep) trial vectors:',NEPROT
          IF(IPAR.GT.0) THEN
            WRITE(LUPRI,'(//A,I3/A/A/)')
     &      '*** WARNING - If your reference state is not the'//
     &         ' excited state with index',IPAR+1,
     &      '    WARNING - then you probably have an instability in'//
     &        ' your reference state !!!',
     &      '    WARNING - (e.g. unrestricted lower energy than'//
     &            ' restricted, or ...)'
!           Analyze instability in more detail
            CALL MEMREL('XRSSTAB',WORK,1,KFRSAV,KFREE,LFREE)
            CALL MEMGET2('REAL','EIG ',KEIG ,N    ,WORK,KFREE,LFREE)
            CALL MEMGET2('REAL','EVEC',KEVEC,LDA*N,WORK,KFREE,LFREE)
            CALL MEMGET2('INTE','IORBCL',KLRO,4*NORBT,WORK,KFREE,LFREE)
            CALL INSTABANA(LDA,N,A,WORK(KEVEC),WORK(KEIG),IPAR,IBTYP,
     &           IBEVC,get_orbital_rotation_indices_pp(),
     &           WORK(KLRO),WORK,KFREE,LFREE)
          ENDIF
        ENDIF
      ENDIF
      CALL MEMREL('XRSSTAB',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/*   Deck instabana */
      SUBROUTINE INSTABANA(LDA,N,REDE,EVEC,EIG,IPAR,IBTYP,IBEVC,JXOPE,
     &                     IORBCL,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Upon detection of instability in reduced Hessian, analyze corresponding solution vector
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxrs.h"
#include "dgroup.h"
#include "dcbibn.h"
#include "dcborb.h"
      PARAMETER(D1=1.0D0)
      DIMENSION REDE(LDA,N),EVEC(LDA,N),EIG(N),IBTYP(2,NREDM),
     &     IBEVC(NREDM),JXOPE(2,NZXOPE),IORBCL(4,NORBT),
     &     WORK(*)
      real*8, allocatable  :: XVEC(:,:),BBUF(:,:)
      integer, allocatable :: INDX(:)
      INTEGER   NCLS(2,3,2)
      CHARACTER*1 TYP(3)
      DATA TYP/'i','a','v'/
C.....First diagonalize reduced Hessian      
      CALL RSJACO(LDA,N,N,REDE,EIG,1,1,0,EVEC)
C      WRITE(LUPRI,*) '* Eigenvalues of the reduced Hessian. NDIM=',N,
C     &           ' NPRED= ',NPRED,NZXOPE
C      DO I = 1,N
C        WRITE(LUPRI,*) I, EIG(I)
C      ENDDO
C.....We next form solution e-e vectors
      allocate(XVEC(NZXOPE,NZ))
      allocate(BBUF(NZXOPEQ,NSTAT))
      allocate(INDX(NZXOPE))
C     Set orbital strings and classify orbitals
      CALL OSTRING(INDSTR,SKIPEE,SKIPEP,-1,WORK,KFREE,LFREE)
      CALL DEFORB(IORBCL,NCLS,INDSTR,WORK,KFREE,LFREE,IPRXRS)
C.....The calls to XRSXV1 assumes that the file PAMBOE exists and is open
      DO I = 1,IPAR
        II = NPRED + I
        CALL XRSXV1(IPX,JBENDX,XVEC,EVEC(1,II),1,IBTYP,IBEVC,
     &              J,BBUF)
        IF(.NOT.STATIC) THEN
          CALL XRSXV1(IMX,JBENDX,BBUF(1,2),EVEC(1,II),1,IBTYP,IBEVC,
     &                J,BBUF)
          CALL DAXPY(NZXOPEQ,D1,BBUF(1,2),1,XVEC,1)          
        ENDIF
C       Calculate norm
        DO J = 1,NZXOPE
          BBUF(J,1) = XVEC(J,1)*XVEC(J,1)
        ENDDO
        DO IZ = 2,NZ
          DO J = 1,NZXOPE
            BBUF(J,1) = BBUF(J,1) + XVEC(J,IZ)*XVEC(J,IZ)
          ENDDO
        ENDDO
        D2TOT=DSUM(NZXOPE,BBUF,1)
        THR   = SQRT(D2TOT/NZXOPE)
        DNTOT = SQRT(D2TOT)
        CALL INDEXX (NZXOPE,BBUF,INDX)
C        WRITE(LUPRI,*) '* XVEC no. ',I
C        DO J = 1,NZXOPE
C          WRITE(LUPRI,*) JXOPE(1,J),JXOPE(2,J),J,BBUF(J,1)
C        ENDDO
        WRITE(LUPRI,'(/3X,A,I3)')
     &    '* Hessian electronic eigenstate nr. ',I
        WRITE(LUPRI,'(3X,A,ES16.8)')
     &    '    Eigenvalue : ',EIG(II)
        WRITE(LUPRI,'(3X,A,ES10.2/)')
     &    '    Excitation amplitudes larger than threshold ',THR
        THR = THR*THR
        DO J = NZXOPE,1,-1
          JJ = INDX(J)
          IF(ABS(BBUF(JJ,1)).LT.THR) GOTO 5
          I1 = JXOPE(1,JJ)
          I2 = JXOPE(2,JJ)
          WRITE(LUPRI,'(6X,I6,3A1,A3,A1,A,I6,3A1,A3,A1,E16.8)')
     &      IORBCL(2,I1),'(',TYP(IORBCL(1,I1)),':',
     &      FREP(IORBCL(3,I1)),')',
     &      ' ---> ',
     &      IORBCL(2,I2),'(',TYP(IORBCL(1,I2)),':',
     &      FREP(IORBCL(3,I2)),')',SQRT(BBUF(JJ,1))
        ENDDO
 5      CONTINUE
      ENDDO
      deallocate(BBUF)
      deallocate(XVEC)
      deallocate(INDX)
      RETURN
      END
      
      
