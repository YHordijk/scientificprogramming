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
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rmolgrd */
      SUBROUTINE RMOLGRD(WORK,LWORK)
C*****************************************************************************
C
C     Calculate molecular gradient.
C
C     Written by J. Thyssen - 1997/05
C
C     Last revision: jth 1997/07/09
C
C*****************************************************************************
         use dirac_cfg
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#include "mxcent.h"
C
      DIMENSION WORK(LWORK)
      LOGICAL   LSAVE(3)
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "blocks.h"
#include "dcbgrd.h"
#include "dcbpsi.h"
#include "symmet.h"
C
      CALL QENTER('RMOLGRD')
#include "memint.h"
C
C
C     Start by writing information to file RMOLGRD:
C
C     - Non-folded and Folded contravariant density matrices in AO basis
C     - Folded contravariant real part of Fock matrix
C
C     80 is the the integer part of 50 * log(5)...
C     (don't drink beer while programming)
C
      LUMGRD = 80
      CALL OPNFIL(LUMGRD,'RMOLGRD','NEW','RMOLGRD')
C
      if (dirac_cfg_hf_calculation .or. dirac_cfg_dft_calculation) then
C
C        Initialize some variables
C
C        If calculation of gradient requested without
C        calculating the wavefunction SETDHF must be called -
C        so we call it anyway.
C
C        HJAAJ July 2001: I have discovered that some of these
C        dirgrd.F routines (e.g. CONMAT_DHF) assume
C        NOPEN = 1 and DF(1) = 0.5 if NASHT .gt. 0 !!!!
C
         IF (NOPEN .EQ. 1) THEN
            IF (NASHT .GT. 1 .OR. ABS(DF(1)-0.5D0) .GT. 1.0D-14) THEN
              WRITE(LUPRI,'(/A/A/A)')
     &        '*** ERROR in RMOLGRD ***',
     &        'Molecular gradient not implemented for open-shell DHF',
     &        'with more than one single occupied open-shell orbital.'
              CALL QUIT('*** ERROR in RMOLGRD: Molecular gradient' // 
     &        ' not implemented for this open-shell DHF')
            END IF
         ELSE IF (NOPEN .GT. 1) THEN
            WRITE(LUPRI,'(/A/A/A)')
     &        '*** ERROR in RMOLGRD ***',
     &        'Molecular gradient not implemented for open-shell DHF',
     &        'with more than one open shell.'
            CALL QUIT('*** ERROR in RMOLGRD: Molecular gradient' // 
     &        ' not implemented for this open-shell DHF')
         END IF

         if (dirac_cfg_dft_calculation) then
            if (nopen > 0) then
               call quit('error in molecular gradient:' // 
     &        ' not implemented for open-shell DFT')
            end if
            if (maxrep > 0) then
               call quit('error in molecular gradient:' // 
     &        ' not implemented for DFT with symmetry')
            end if
         end if

         if (bss .or. x2c) then
            call quit('error in molecular gradient:' // 
     &      ' not implemented for X2C')
         end if

         CALL SETDHF(IPRGRD)
C
C        Density and Fock matrices for (Dirac-)Hartree-Fock
C
         CALL CONMAT_DHF(IPRGRD,WORK(KFREE),LFREE)
C
      ELSE
C
         WRITE(LUPRI,'(/A/A)')
     &        '*** ERROR in RMOLGRD ***',
     &        'Molecular gradient only implemented for DHF'
         CALL QUIT('*** ERROR in RMOLGRD: ' // 
     &             'Molecular gradient only implemented for DHF')
C
      END IF
C
C
      CALL TITLER('Analytical calculation of molecular gradient',
     &            '*',110)
C
C     Calculate gradient
C
      IGRD_INTBUF = IGRD_INTFLG
      CALL RMOLGRD1(WORK(KFREE),LFREE,IPRGRD)
C
      CLOSE (LUMGRD, STATUS = 'DELETE' )
C
C     Output of gradient
C
      CALL GRDOUT(WORK(KFREE),LFREE,IPRGRD)
C
C
C
C     Memory deallocation
      CALL MEMREL('RMOLGRD.DHFOUT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('RMOLGRD')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck conmat_dhf */
      SUBROUTINE CONMAT_DHF(IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Calculate various contravariant density and Fock matrices.
C
C     Memory usage is as follows:
C
C     Fock matrix contruction:
C
C     1. Covariant AO Fock matrices:
C        3 * NZ * N2BBASX  (3 Fock matrices)
C     2. Transformation to MO basis:
C     3. Transformation to contravariant AO basis
C     4. Fold matrix
C     5. Insert half phases
C     6. Transform from QO to AO basis
C        2 * NZ   N2BBASX  +  2 * NZ * N2ORBT (4 Fock matrices)
C     
C     Density matrix contructions:
C
C     7. Total folded covariant density matrix:
C        3 * NZ * N2BBASX  (3 Fock matrices)
C     8. Individual non-folded covariant density matrices:
C        3 * NZ * N2BBASX  (3 Fock matrices)
C
C     Density matrix construction
C
C     J. Thyssen, 2000/8/27
C
C*****************************************************************************
      use quaternion_algebra
      use xcint_main
      use dirac_cfg
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
      PARAMETER ( D1 = 1.00D00, DP5 = 0.5D00, D0 = 0.0D00,
     $     D2 = 2.00D00 )
C
#include "dcbham.h"
#include "nuclei.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbgrd.h"
#include "symmet.h"
C
      DIMENSION WORK(*)
      LOGICAL TOBE
      CHARACTER CTMP*1, M(4)*2

      real(8), allocatable :: dmat(:)
      real(8), allocatable :: cmo(:)
C
      DATA M /'AR','AI','BR','BI'/
C
#include "ibtfun.h"
#include "memint.h"
C
      CALL QENTER('CONMAT_DHF')
      IF ( IPRINT .GT. 5 )
     $     CALL HEADER('Output from CONMAT_DHF',-1)
C
      LUBUF = 81
C
C     1. Start by constructing AO Fock matrices in covariant basis.
C     =============================================================
C
      CALL MEMGET('REAL',KFAO,N2BBASXQ * 2,WORK,KFREE,LFREE)
      KFCAO = KFAO
      KFVAO = KFAO + N2BBASXQ
C     ... KFVAO is also used for CMO later
      CALL MEMGET('REAL',KH1AO,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Two-electron Fock matrix/matrices
C
      INQUIRE(FILE='DFFCK2',EXIST=TOBE)
      IF(TOBE) THEN
         CALL OPNFIL(LUBUF,'DFFCK2','OLD','CONMAT_DHF')
         WRITE(LUPRI,'(/A)')
     &   '* Reading two-electron Fock matrix for DHF molecular gradient'
         CALL REAFCK(LUBUF,WORK(KFCAO),.TRUE.,NFMAT)
         CLOSE(LUBUF,STATUS='KEEP')
      ELSE
         WRITE(LUPRI,'(2(/A))')
     &        '*** ERROR in CONMAT_DHF ***',
     $        'Gradient calc. without DFFCK2 is not implemented yet'
         CALL QUIT('*** ERROR in CONMAT_DHF ***')
      END IF
C
C     One electron Fock Matrix
C
      INQUIRE(FILE='DFFCK1',EXIST=TOBE)
      IF (TOBE) THEN
         CALL OPNFIL(LUBUF,'DFFCK1','OLD','CONMAT_DHF')
         WRITE(LUPRI,'(/A)')
     &   '* Reading one-electron Fock matrix for DHF molecular gradient'
         CALL REAFCK(LUBUF,WORK(KH1AO),.TRUE.,1)
         CLOSE(LUBUF,STATUS='KEEP')
      ELSE
         WRITE(LUPRI,'(/A)')
     &    '*** WARNING *** CONMAT_DHF: No 1-Fock found. Regenerating.'
         CALL ONEFCK(WORK(KH1AO),IPRINT,WORK(KFREE),LFREE)
      END IF
C
C     Correct Fock matrices for ONECAP
C
      IF (ONECAP) THEN
         CALL MEMGET('REAL',KSCQ,NUCIND,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KDSO,N2BBASXQ,WORK,KFREE,LFREE)
c        ... Get density
         CALL GETDSO(WORK(KDSO),WORK(KFREE),LFREE,IPRINT)
c        ... Calculate and add correction.
         CALL CATCORR(WORK(KH1AO),WORK(KFCAO),WORK(KDSO),
     &                WORK(KFREE),LFREE,IPRINT)
      END IF

      gradxc = 0.0d0
      if (dirac_cfg_dft_calculation) then

         allocate(cmo(n2bbasxq))
         cmo = 0.0d0
         call opnfil(lucoef, 'DFCOEF', 'OLD', 'xc molecular gradient')
         call reacmo(lucoef, 'DFCOEF',cmo, dummy, dummy, dummy, 3)
         close(lucoef, status='keep')
       
         allocate(dmat(n2bbasxq))
         dmat = 0.0d0
         call genden(dmat, cmo, 1, 0)
         deallocate(cmo)

#ifdef VAR_MPI
         if (parcal) call dirac_parctl( XCINT_PAR )
#endif
!fixme: there is some issue with one-step - investigate later
!       not correct for gga
!        in one go integrate xc potential
!        and xc contribution to molecular gradient
!        call integrate_xc(xc_mat_dim           = ntbas(0),
!    &                     xc_nz                = nz,
!    &                     xc_dmat_0            = dmat,
!    &                     xc_nr_dmat           = 0, 
!    &                     xc_nr_fmat           = 1,
!    &                     xc_fmat              = work(kh1ao),
!    &                     xc_do_potential      = .true.,
!    &                     xc_nr_atoms          = nucdep, 
!    &                     xc_property_gradient = gradxc,
!    &                     xc_do_geo_0          = .true.)
         call integrate_xc(xc_mat_dim           = ntbas(0),
     &                     xc_nz                = nz,
     &                     xc_dmat_0            = dmat,
     &                     xc_nr_dmat           = 0, 
     &                     xc_nr_fmat           = 1,
     &                     xc_fmat              = work(kh1ao),
     &                     xc_do_potential      = .true.)
#ifdef VAR_MPI
         if (parcal) call dirac_parctl( XCINT_PAR )
#endif
         call integrate_xc(xc_mat_dim           = ntbas(0),
     &                     xc_nz                = nz,
     &                     xc_dmat_0            = dmat,
     &                     xc_nr_dmat           = 0, 
     &                     xc_nr_fmat           = 0,
     &                     xc_nr_atoms          = nucdep, 
     &                     xc_property_gradient = gradxc,
     &                     xc_do_geo_0          = .true.)
         deallocate(dmat)

      end if



C     FCao = H1ao + QCao
C
      CALL DAXPY(N2BBASXQ,D1,WORK(KH1AO),1,WORK(KFCAO),1)
      CALL MEMREL('CONMAT_DHF 1',WORK,KWORK,KH1AO,KFREE,LFREE)
C
C     2. Transform to MO basis
C     ========================
C
      CALL MEMGET('REAL',KFMO,MAX(N2ORBTQ,NCMOTQ),WORK,KFREE,LFREE)
      KCMO = KFMO
      CALL MEMGET('REAL',KFMOTMP,N2ORBTQ,WORK,KFREE,LFREE)
C
      IF (NASHT.GT.0) THEN
C
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CONMAT_DHF')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUMMY,DUMMY,DUMMY,3)
         CLOSE(LUCOEF,STATUS='KEEP')
C
C
C        Transform QV to MO basis
C
         DO 10 I = 1,NFSYM
            IF (NASH(I).EQ.0) GOTO 10
            CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &           NORB(I),NORB(I),
     &           WORK(KFVAO + I2BASX(I,I)),NTBAS(0),NTBAS(0),
     &           NZ,IPQTOQ(1,0),
     &           WORK(KFMOTMP +I2ORBT(I)),NORB(I),NORB(I),
     &           NZ,IPQTOQ(1,0),
     &           WORK(KCMO+ICMOQ(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KCMO+ICMOQ(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
            IF (IPRINT.GE.10) THEN
               WRITE(LUPRI,*) 'QV in MO-basis, irrep ',I
               CALL PRQMAT(WORK(KFMOTMP+I2ORBT(I)),NORB(I),NORB(I),
     &              NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
 10      CONTINUE
C
         CALL DZERO(WORK(KFMO),N2ORBTQ)
C
         DO 11 I = 1,NFSYM
            IF (NASH(I).EQ.0) GOTO 11
            IADR = I2ORBT(I)+(NISH(I)+NPSH(I))*NORB(I)
            CALL DCOPY(NASH(I)*NORB(I),WORK(KFMOTMP+IADR),1,
     &                                 WORK(KFMO+IADR),1)
            IF (IPRINT.GE.10) THEN
               WRITE(LUPRI,*) 'Active part of F, irrep ',I
               CALL PRQMAT(WORK(KFMO+I2ORBT(I)),NORB(I),NORB(I),
     &              NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
 11      CONTINUE
C
C        FDao = H1ao + QCao + QVao
C             = FCao + QVao
C
         CALL DAXPY(N2BBASXQ,DP5,WORK(KFVAO),1,WORK(KFCAO),1)
C
      ELSE
         CALL DZERO(WORK(KFMO),N2ORBTQ)
      END IF
C
C
      IF (NISHT.GT.0) THEN
         DO 20 I = 1,NFSYM
C
            KCMO = KFVAO
            CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CONMAT_DHF')
            CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUMMY,DUMMY,DUMMY,3)
            CLOSE(LUCOEF,STATUS='KEEP')
C
            IF (NISH(I).EQ.0) GOTO 20
            CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &           NORB(I),NORB(I),
     &           WORK(KFCAO+I2BASX(I,I)),NTBAS(0),NTBAS(0),
     &           NZ,IPQTOQ(1,0),
     &           WORK(KFMOTMP+I2ORBT(I)),NORB(I),NORB(I),
     &           NZ,IPQTOQ(1,0),
     &           WORK(KCMO+ICMOQ(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KCMO+ICMOQ(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
            IF (IPRINT.GE.10) THEN
               WRITE(LUPRI,*) 'FD in MO-basis, irrep ',I
               CALL PRQMAT(WORK(KFMOTMP+I2ORBT(I)),NORB(I),NORB(I),
     &              NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
            IADR = I2ORBT(I)+NPSH(I)*NORB(I)
            CALL DSCAL(NISH(I)*NORB(I),D2,WORK(KFMOTMP+IADR),1)
            CALL DCOPY(NISH(I)*NORB(I),WORK(KFMOTMP+IADR),1,
     &           WORK(KFMO+IADR),1)
            IF (IPRINT.GE.10) THEN
               WRITE(LUPRI,*) 'Total Fock matrix, irrep ',I
               CALL PRQMAT(WORK(KFMO+I2ORBT(I)),NORB(I),NORB(I),
     &                    NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
 20      CONTINUE
      END IF
C
C
C     3. Construct contravariant Fock matrix in AO-basis
C     ===============================================
C
C     FOCK = CCtFCCt = DFDt = CF(MO)Ct
C
      CALL DZERO(WORK(KFCAO),N2BBASXQ)
      KCMO = KFVAO
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CONMAT_DHF')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUMMY,DUMMY,DUMMY,3)
      CLOSE(LUCOEF,STATUS='KEEP')
      IF ( IPRINT .GE. 10 ) THEN
         DO I=1,NFSYM
           IF (NORB(I) .NE. 0) THEN
           WRITE(CTMP,'(I1)') I
           CALL HEADER('Coefficient matrix, irrep '//CTMP,-1)
           CALL PRQMAT(WORK(KCMO+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &                 NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
           END IF
         END DO
      END IF
C
      DO I= 1,NFSYM
         IF (NORB(I) .NE. 0) THEN
            NZA = NZ
            CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &           NORB(I),NORB(I),
     &           WORK(KFAO+I2BASX(I,I)),NTBAS(0),NTBAS(0),
     &           NZ,IPQTOQ(1,0),
     &           WORK(KFMO+I2ORBT(I)),NORB(I),NORB(I),
     &           NZ,IPQTOQ(1,0),
     &           WORK(KCMO + ICMOQ(I)),
     $           NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KCMO + ICMOQ(I)),
     $           NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
         END IF
      END DO
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('CONMAT_DHF: FOCK After QTRANS to AO-basis',-1)
         CALL FLSHFO(LUPRI)
         CALL PRQMAT(WORK(KFAO),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Transform to unsorted basis
C
      CALL BSTOBU(WORK(KFAO),NZ,WORK(KFREE),LFREE)
      IF(IPRINT.GE.5) THEN
         CALL HEADER('CONMAT_DHF: FOCK After BSTOBU to unsorted basis',
     &        -1)
         CALL PRQMAT(WORK(KFAO),NTBAS(0),NTBAS(0),
     &               NTBAS(0),NTBAS(0),NZ,
     &               IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     4. Fold Fock matrix
C     ================
C
      DO IZ = 1, NZ
         CALL RMOLGRD_FOLDMAT(WORK(KFCAO + ( IZ - 1 ) * N2BBASX ),
     &        NTBAS(0), (IZ .EQ. 1) )
      END DO
C
      IF ( IPRINT .GT. 5 ) THEN
         CALL HEADER('CONMAT_DHF: FOCK After folding',-1)
         CALL PRQMAT(WORK(KFCAO),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &        NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C     5. Insert half-phases
C     ==================
C
      IF(NZ.LT.4) THEN
         DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,0)
            CALL Q2PHASE('F',IQ,1,WORK(KFCAO + N2BBASX * (IZ - 1) ) )
         ENDDO
         IF(IPRINT.GE.5) THEN
            CALL HEADER('CONMAT_DHF: '//
     $           'Unsorted FOCK with half-phases inserted',-1)
            CALL PRQMAT(WORK(KFCAO),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NZ,
     &                  IPQTOQ(1,0),LUPRI)
         ENDIF
      ENDIF
C
C     6. Transform fock matrix from QO to AO
C     ======================================
C
C     We are only interested in the real part, as it is used in the product:
C
C     F * dS/dX, 
C
C     where S is the overlap matrix, which is real.
C
      CALL DZERO(WORK(KFVAO),N2BBASX)
      CALL DSYM1R(WORK(KFCAO),WORK(KFVAO),NTBAS(0),0,IPRINT - 20)
C
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('Final folded contravariant AO Fock matrix',-1)
         CALL PRQMAT(WORK(KFVAO),NTBAS(0),NTBAS(0),
     &        NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Write to file
C     =============
C
      REWIND (LUMGRD)
      CALL NEWLAB('FKFLDCT1',LUMGRD,LUPRI)
      CALL WRITT(LUMGRD,N2BBASX,WORK(KFVAO))
C
      CALL MEMREL('CONMAT_DHF.after fock',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     7. Folded contravariant total AO density matrix
C     ===============================================
C
      CALL MEMGET('REAL',KDAO,3 * N2BBASXQ,WORK,KFREE,LFREE)
      KDCAO = KDAO
      KDVAO = KDAO + N2BBASXQ
      KCMO  = KDVAO + N2BBASXQ
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CONMAT_DHF')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUMMY,DUMMY,DUMMY,3)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Get density matrix
C
      CALL GENDEN(WORK(KDCAO),WORK(KCMO),1,IPRINT)
      IF (NASHT.GT.0) THEN
         CALL GENDEN(WORK(KDVAO),WORK(KCMO),2,IPRINT)
C
C        DVao + DCao
C
         CALL DAXPY(N2BBASXQ,DP5,WORK(KDVAO),1,
     &                          WORK(KDCAO),1)
      END IF
C
C     Since Trond insists in DC_ij = \sum_a C_ia C_ja^*
C     and .NOT. DC_ij = 2 \sum C_ia C_ja^*, we multiply with 2.0
C
      CALL DSCAL(N2BBASXQ,D2,WORK(KDCAO),1)
C
C     Transform to unsorted basis
C
      CALL BSTOBU(WORK(KDCAO),NZ,WORK(KFREE),LFREE)
      IF(IPRINT.GE.10) THEN
         CALL HEADER('CONMAT_DHF: DMAT After BSTOBU to unsorted basis',
     &               -1)
         CALL PRQMAT(WORK(KDCAO),NTBAS(0),NTBAS(0),
     &        NTBAS(0),NTBAS(0),NZ, IPQTOQ(1,0),LUPRI)
      END IF
C
C     Fold the damn thing
C
c      CALL RMOLGRD_FOLDMAT(WORK(KDCAO),NTBAS(0),NZ)
C
C     Insert half-phases
C
      IREP = ISYMOP(1)-1
      IF(NZ.LT.4) THEN
         DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            CALL Q2PHASE('D',IQ,1,WORK(KDCAO + N2BBASX * (IZ - 1) ) )
         ENDDO
         IF(IPRINT.GE.5) THEN
            CALL HEADER('Unsorted DMAT with half-phases inserted',-1)
            CALL PRQMAT(WORK(KDCAO),NTBAS(0),NTBAS(0),NTBAS(0),
     &           NTBAS(0),NZ, IPQTOQ(1,IREP),LUPRI)
         ENDIF
      END IF
C
C
C     Note ad Imaginary part of the density matrix:
C
C     In the AVEMOM routine X, Y-, and Z-derivates are multiplied 
C     with the k-part, j-part, and the i-part, respectively.
C     
C
C     The real part:
C     --------------
C
      CALL DCOPY(N2BBASX, WORK(KDCAO), 1, WORK(KCMO), 1)
      CALL DZERO(WORK(KDVAO),N2BBASX)
      CALL RMOLGRD_FOLDMAT( WORK(KCMO), NTBAS(0), .TRUE. )
      CALL DSYM1R(WORK(KCMO),WORK(KDVAO),NTBAS(0),0,IPRINT - 20)
      CALL NEWLAB('DTFLDCT1',LUMGRD,LUPRI)
      CALL WRITT(LUMGRD,N2BBASX,WORK(KDVAO))
C
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('CONMAT_DHF: Final DMAT:',-1)
         WRITE(LUPRI,'(3A)') '*** ',M(1), ' part ***'
         CALL OUTPUT(WORK(KDVAO),1,NTBAS(0),1,NTBAS(0),
     $        NTBAS(0),NTBAS(0),-1,LUPRI)
      END IF
C
C     The imaginary parts:
C     --------------------
C
      IC = 2
C
      DO I = 2,4
C
C        Direction:
C        IDIR = 1(X), 2(Y), or 3(Z)
C
         IDIR = 5-I
C
C
C        Symmetry of axis:
C
         ISYM = ISYMAX(IDIR,1)
C
C        Output goes to...
C
         IZ = JM4POS(I-1)
C
C        from...
C
         IQ = IQMULT(1,JQBAS(ISYM,IC),IZ)
         IQP = IQTOPQ(IQ,0)
         IREPD = IRQMAT(IZ,0)
C
         CALL DZERO(WORK(KDVAO),N2BBASX)
C
         CALL DSYM1R(WORK(KDCAO + (IQP - 1) * N2BBASX),WORK(KDVAO),
     &              NTBAS(0),IREPD,IPRINT - 20)
C
         CALL RMOLGRD_FOLDMAT( WORK(KDVAO), NTBAS(0), .FALSE. )
C
         CALL NEWLAB('DTFLDCT'//CHAR(I+48),LUMGRD,LUPRI)
         CALL WRITT(LUMGRD,N2BBASX,WORK(KDVAO))
C
         IF ( IPRINT .GE. 10) THEN
            WRITE(LUPRI,'(3A)') '*** ',M(I), ' part ***'
            CALL OUTPUT(WORK(KDVAO),1,NTBAS(0),1,NTBAS(0),
     $           NTBAS(0),NTBAS(0),-1,LUPRI)
         END IF
C     
      END DO
C
C
C
C     8. Non-folded contravariant individual AO density matrices
C     ==========================================================
C
C
C     Contruct density matrix
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CONMAT_DHF')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUMMY,DUMMY,DUMMY,3)
      CLOSE(LUCOEF,STATUS='KEEP')
      DO IMAT = 1,NFMAT
         IBIT = IBTSHL(1,IMAT-1)
         CALL GENDEN(WORK(KDCAO),WORK(KCMO),IBIT,IPRINT)
         IF (IPRINT .GE. 5) THEN
            CALL HEADER('CONMAT_DHF: DMAT from GENDEN',-1)
            CALL PRQMAT(WORK(KDCAO),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  IPQTOQ(1,0),LUPRI)
         END IF
C
C        Scale inactive density matrix with D2 (= occupation)
C
         IF (IMAT.EQ.1)
     &      CALL DSCAL(N2BBASXQ,D2,WORK(KDCAO),1)
C
C        Transform to unsorted basis
C
         IREP = ISYMOP(IMAT)-1
         CALL BSTOBU(WORK(KDCAO),NZ,WORK(KFREE),LFREE)
         IF(IPRINT.GE.5) THEN
            CALL HEADER('CONMAT_DHF:'//
     $           ' DMAT After BSTOBU to unsorted basis', -1)
            WRITE(LUPRI,'(A,I5)') '*** Matrix no. ',IMAT
            CALL PRQMAT(WORK(KDCAO),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
         END IF
C
C        Transform density matrix from QO to AO
C        ======================================
C
C
C        The real part:
C        --------------
C
         CALL DZERO(WORK(KDVAO),N2BBASX)
         CALL DTSOAO(WORK(KDCAO),WORK(KDVAO),NTBAS(0),0,IPRINT)
         CALL NEWLAB('D'//CHAR(IMAT+48)//'   CT1',LUMGRD,LUPRI)
         CALL WRITT(LUMGRD,N2BBASX,WORK(KDVAO))
C
         IF ( IPRINT .GE. 10 ) THEN
            CALL HEADER('CONMAT_DHF: Non-folded DMAT '//
     $           CHAR(IMAT+48)//' :',-1)
            WRITE(LUPRI,'(3A)') '*** ',M(1), ' part ***'
            CALL OUTPUT(WORK(KDVAO),1,NTBAS(0),1,NTBAS(0),
     $           NTBAS(0),NTBAS(0),-1,LUPRI)
         END IF
C
C        The imaginary parts:
C        --------------------
C
         IC = 2
C
         DO I = 2,4
C
C           Direction:
C           IDIR = 1(X), 2(Y), or 3(Z)
C
            IDIR = 5-I
C
C           Symmetry of axis:
C
            ISYM = ISYMAX(IDIR,1)
C
C           Output goes to...
C
            IZ = JM4POS(I-1)
C
C           from...
C
            IQ = IQMULT(1,JQBAS(ISYM,IC),IZ)
            IQP = IQTOPQ(IQ,0)
            IREPD = IRQMAT(IZ,0)
C     
            CALL DZERO(WORK(KDVAO),N2BBASX)
C     
            CALL DTSOAO(WORK(KDCAO + (IQP - 1) * N2BBASX),WORK(KDVAO),
     &           NTBAS(0),IREPD,IPRINT - 20)
C     
            CALL NEWLAB('D'//CHAR(IMAT+48)//'   CT'//CHAR(I+48),
     &           LUMGRD,LUPRI)
            CALL WRITT(LUMGRD,N2BBASX,WORK(KDVAO))
C     
            IF ( IPRINT .GE. 10) THEN
               WRITE(LUPRI,'(3A)') '*** ',M(I), ' part ***'
               CALL OUTPUT(WORK(KDVAO),1,NTBAS(0),1,NTBAS(0),
     $              NTBAS(0),NTBAS(0),-1,LUPRI)
            END IF
C     
         END DO
      ENDDO
C
      CALL MEMREL('CONMAT_DHF.after dmat',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('CONMAT_DHF')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rmolgrd_foldmat */
      SUBROUTINE RMOLGRD_FOLDMAT(A,N,FLAG)
C*****************************************************************************
C
C
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER ( D0 = 0.0D00 )
C
      DIMENSION A(N,N)
      LOGICAL FLAG
C
C     FLAG: symmetric folding ( FLAG = true )
C         : anti-symmetric folding ( FLAG = false )
C
      IF (FLAG) THEN
         DO I = 2, N
            DO J = 1, I - 1
               A(I,J) = A(I,J) + A(J,I)
               A(J,I) = D0
            END DO
         END DO
      ELSE
         DO I = 2, N
            DO J = 1, I - 1
               A(I,J) = A(I,J) - A(J,I)
               A(J,I) = D0
            END DO
         END DO
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rmolgrd1 */
      SUBROUTINE RMOLGRD1(WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Calculate molecular gradient
C
C     J. Thyssen, 1997/18/04
C    
C     Last revision: Aug 28 2000 jth, general cleanup
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      PARAMETER(D0 = 0.0D0)
C
      DIMENSION WORK(LWORK)
      DIMENSION DNTSKP(3,14,3)
      DIMENSION IQBF(4)
      DIMENSION IFCTYP(4*2)
      CHARACTER CTMP*1, STR*3, GABLAB*8, OMITVNUC(2)*4
      LOGICAL   DOINT(2,2)
      LOGICAL   DOLS,DOSS,LBIT,DOGRADLV
      LOGICAL   NODV,NOPV,FIRST,NOCONT, NOINFO
      LOGICAL   FNDLAB
C     Timing variables
      DOUBLE PRECISION DPRTIM(4,2)
      SAVE DPRTIM
      DATA DPRTIM /D0,D0,D0,D0,D0,D0,D0,D0/
C
#include "dcbbas.h"
#include "nuclei.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dorps.h"
C two el. part
#include "ccom.h"
#include "cbihr2.h"
#include "blocks.h"
#include "dcbgrd.h"
C
#include "exeinf.h"
#include "symmet.h"
#include "chrxyz.h"
C
#include "ibtfun.h"
C
      CALL QENTER('RMOLGRD1')
C
#include "memint.h"
C
      CALL GETTIM(CPUTO1,WLLTO1)
      MAXDER = 1
C  
C     Initialisation of timings
C
      CPUONE = D0
      CPUTLL = D0
      CPUTLS = D0
      CPUTSS = D0
      CPUTOT = D0
      WLLONE = D0
      WLLTLL = D0
      WLLTLS = D0
      WLLTSS = D0
      WLLTOT = D0
      CPUGAB = D0
      WLLGAB = D0
C
      NOINFO = .TRUE.
      GRD_NOSMLV = .FALSE.
      GRD_NOSMLS = .FALSE.
      OMITVNUC(1) = 'FFFF'
      OMITVNUC(2) = 'FFFF'
      DOGRADLV  = (ONECAP .AND. INTV1C .EQ.2) .OR. DOLVC .OR.
     &   ( LBIT( IGRD_INTBUF, 2 ) .AND. LBIT( IGRD_INTFLG, 3 ) )
C    
      DO J = 1,2
         DO I = 1,2
            DOINT(I,J) = .TRUE.
         END DO
      END DO
C
      WRITE(LUPRI,'(A)') 
     &   ' Special information about calculation '//
     &   'of the molecular gradient:'
      IF (LEVYLE) THEN
         WRITE(LUPRI,'(/A,2(/,4X,A))') 
     &      ' Levy-Leblond Hamiltonian: ',
     &      '- No small-small nuclear attraction integrals.',
     &      '- No small-small reorthonormalization gradient'
         NOINFO = .FALSE.
         GRD_NOSMLS = .TRUE.
         GRD_NOSMLV = .TRUE.
         IF (.NOT. NOSMLV) THEN
            WRITE(LUPRI,'(/2A)')
     &         ' WARNING: small-small nuclear attraction integrals AND',
     &         ' Levy-Leblond Hamiltonian: Inconsistent...'
         END IF
         IGRD_INTBUF = IBTAND( 1, IGRD_INTBUF )
      ELSE IF (NOSMLV) THEN
         WRITE(LUPRI,'(/A)') 
     &      ' - No small-small nuclear attraction integrals.'
         NOINFO = .FALSE.
         GRD_NOSMLV = .TRUE.
      END IF
C
      IF (N1OPER .GT. 0) THEN
         WRITE(LUPRI,'(/A,I3)') 
     &      ' - Number of finite field operators:',N1OPER
C TODO   insert list of finite fields here ??? /hjaaj
         NOINFO = .FALSE.
      END IF
C
      IF ( DOLVC .AND. .NOT.(ONECAP.AND.INTV1C.EQ.2) .AND.
     &   LBIT( IGRD_INTBUF, 3 ) ) THEN
         IGRD_INTBUF = IBTAND ( 1+2, IGRD_INTBUF )
         WRITE(LUPRI,'(/A/A/)')
     &     ' INFO: LVCORR used in wavefunction, but not in gradient =>',
     &     ' LVCORR gradient used (SS-gradient is not calculated)'
         NOINFO = .FALSE.
      END IF
      IF (ONECAP) THEN
         IF (INTV1C.EQ.1) THEN
            WRITE(LUPRI,'(/A)')
     &      ' - Only one-center contributions to d/dXA Vss.'
            NOINFO = .FALSE.
         ELSE IF (INTV1C.EQ.2) THEN
            WRITE(LUPRI,'(/A/A)') 
     &      ' - ONECAP type 2 : Only atomic contributions to '//
     &      ' d/dXA Vss (thus no LS ','   and SS 2-electron gradient'//
     &      ' and no small-small nuclear attraction gradient).'
            NOINFO = .FALSE.
         ELSE IF (INTV1C.EQ.3) THEN
            WRITE(LUPRI,'(/A)') 
     &      ' - ONECAP type 3 : Only atomic contributions to '//
     &      ' d/dXA Vss (thus no LS and SS 2-electron gradient '
            NOINFO = .FALSE.
         END IF
      END IF
      IF (NOINFO) WRITE(LUPRI,'(/1X,A)')
     &   'No special information in this calculation!'
C
      IF (N1OPER .GT. 0) THEN
         WRITE (LUPRI,'(/A)')
     &   ' ERROR: molecular gradient has not been implemented yet'//
     &   ' for finite field calculations.'
         CALL QUIT('Molecular grad. not implemented for finite fields')
      END IF
C
C     Do totally symmetric contributions (RSETSYM(.FALSE.,...)
C     and calculate translation/rotation invariance
C
      DO I = 1, 3*MXCENT
         DOPERT(I,1) = .TRUE.
         DOPERT(I,2) = .TRUE.
      END DO
      FTRONV = .TRUE.
      RELGRD = .TRUE.
      CALL RSETSYM(.FALSE.,RELGRD,.FALSE.,.FALSE.)
      CALL RTROINI(.FALSE.)
      CALL TROINV(WORK(KFREE),LFREE)
C
C     Nuclear repulsion
C     =================
C
      CALL NUCINI
      CALL MEMGET('REAL',KHESNN,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCSTRA,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCTRA,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
      CALL NUCREP(WORK(KHESNN),WORK(KCSTRA),WORK(KSCTRA))
      CALL MEMREL('RMOLGRD1.NUCREP',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Calculated GRADLV:
C     "SS-integral contribution by classical repulsion of 
C      small component atomic charges"-gradient
C      or ONECAP model 2 to gradient
C
C     GRADLV will only be used if relevant, even if calculated here.
C
      IF (DOGRADLV) THEN
C           ... get small and large component charges
C               for one-center model 2 (and maybe type3!)
         CALL MEMGET('REAL',KSCQ,NUCIND,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KDSO,N2BBASXQ,WORK,KFREE,LFREE)
         CALL GETDSO(WORK(KDSO),WORK(KFREE),LFREE,IPRINT)
         CALL SCDENS(WORK(KSCQ),WORK(KDSO), WORK(KFREE),LFREE,IPRINT)
         CALL MEMREL('RMOLGRD1.dolvc',WORK,KWORK,KDSO,KFREE,LFREE)
         CALL LVCGRD(IPRINT,WORK(KSCQ),WORK(KFREE),LFREE)
      END IF
C
C     One-center model 2 correction to SS-nuclear attration
C     ================================================
C     
C     if one-center model 1 : LVCGRD is called with modified nuclear charges.
C     The correction is therefore stored in  GRADLV as :   
C     GRADLV = GRADLV - GRADNN
C     
      IF (ONECAP .AND. INTV1C .EQ. 2) THEN
         DO I = 1,3*NUCDEP
            GRADLV(I) = GRADLV(I) - GRADNN(I)
         END DO
      END IF
C
C     Set IFCTYP
C     ==========
C
      IOFF = 0
      DO IMAT = 1,NFMAT
C        Density matrix 1: Coulomb and Exchange
         IFCTYP(1+IOFF) = 13
C        Density matrix 2-4: Only Exchange
         IFCTYP(2+IOFF) = 22
         IFCTYP(3+IOFF) = 22
         IFCTYP(4+IOFF) = 22
         IOFF = IOFF + 4
      END DO
C
C     One-center corrections to the nuclear-attraction gradient.
C     ===========================================================
C
      IF (ONECAP) THEN
         IF (INTV1C .EQ. 1) THEN
c           ... one-center model 1 (Luuk's)
            WRITE(LUPRI,'(/1X,A)') 
     &      '- Only one-center contributions to d/dXA Vss.'
            OMITVNUC(1) = 'TFFF'
            OMITVNUC(2) = 'TFTT'
         ELSE IF (INTV1C .EQ. 2) THEN
            IGRD_INTBUF = IBTAND ( 1, IGRD_INTBUF )
C           omit all SS, use corrected charges for LL
            GRD_NOSMLV = .TRUE.
            OMITVNUC(1)   = 'TFFF'
            OMITVNUC(2)   = 'TTTT'
            CALL MEMGET('REAL',KNQ,NUCIND,WORK,KFREE,LFREE)
C
C           Store the original nuclear charges for later use,
C           and make the one-center model 2 corrected nuclear charges.
C
            DO I=1,NUCIND
               WORK(KNQ-1+I) = CHARGE(I)
               CHARGE(I) = WORK(KSCQ-1+I) + CHARGE(I)
            END DO
         ELSE IF (INTV1C .EQ. 3) THEN
            IGRD_INTBUF = IBTAND ( 1, IGRD_INTBUF )
C           use corrected charges for LL and SS 
C           (Corrections now added to two-electron Fock-matrix - the factor
C            one half needed is tranferred to the charges here.)   
C            --- ONCE SS-PART IS IMPLEMENTED !!!
c            OMITVNUC   = 'TFFF'
c            CALL MEMGET('REAL',KNQ,NUCIND,WORK,KFREE,LFREE)
C
C           Store the original nuclear charges for later use,
C           and make the one-center model 3 corrected nuclear charges.
C
            WRITE(LUPRI,'(/1X,A,A)') 
     &      ' Warning !!!  Molecuar gradient not implemented for ',
     &      'one-center model 3 yet.'
             CALL QUIT('*** ERROR in RMOLGRD1 ***')
         ELSE 
            WRITE(LUPRI,'(/1X,A)') 
     &      '-- Molecular gradient not implemented for specified model!'
            CALL QUIT('*** ERROR in RMOLGRD1 ***')
         END IF
         NOINFO = .FALSE.
      ELSE
         OMITVNUC(1) = 'FFFF'
         OMITVNUC(2) = 'FFFF'
      END IF
C    
C     One electron part of gradient
C     =============================
C
      CALL HEADER('RMOLGRD1: One-electron part',-1)
      IF ( DPRTIM(1,1) .GT. D0 ) THEN
         CALL TSTAMP(' ',LUPRI)
         WRITE(LUPRI,'(/1X,A,F10.2,A)') 
     &               'Estimated CPU time:  ',DPRTIM(1,1),' seconds.'
         WRITE(LUPRI,'(1X,A,F10.2,A/)') 
     &               'Estimated WALL time: ',DPRTIM(1,2),' seconds.'
      END IF
C
C
      CALL GETTIM(CPU1,WALL1)
C
C     Get contravariant Fock matrix in symmetry distinct AO-basis
C
      CALL MEMGET('REAL',KFAO,N2BBASX,WORK,KFREE,LFREE)
C
      REWIND (LUMGRD)
      IF ( FNDLAB('FKFLDCT1',LUMGRD) ) THEN
         CALL READT(LUMGRD,N2BBASX,WORK(KFAO))
      ELSE
         WRITE(LUPRI,'(/A/A)')
     $        '*** ERROR in RMOLGRD1 ***',
     $        'label "FKFLDCT1" not found on file RMOLGRD'
         CALL QUIT('*** ERROR in RMOLGRD1 ***')
      END IF
C
C     Get folded total contravariant density matrix
C     in symmetry distinct AO-basis
C
      CALL MEMGET('REAL',KDTAO,4 * N2BBASX,WORK,KFREE,LFREE)
C
      IOFF = 0
      DO I = 1, 4
         REWIND (LUMGRD)
         IF ( FNDLAB('DTFLDCT'//CHAR(I+48),LUMGRD) ) THEN
            CALL READT(LUMGRD,N2BBASX,WORK(KDTAO+IOFF))
         ELSE
            WRITE(LUPRI,'(/A/A)')
     $           '*** ERROR in RMOLGRD1 ***',
     $           'label "DTFLDCT'//CHAR(I+48)//'"'//
     &           ' not found on file RMOLGRD'
            CALL QUIT('*** ERROR in RMOLGRD1 ***')
         END IF
         IOFF = IOFF + N2BBASX
      END DO
C
C
C     Call ONEDRV for construction of gradient
C
      CALL MEMGET('REAL',KSTHMA,3*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFACIN,2*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
C
C     PROPTY = .TRUE.
C     MAXDIF = 1
C     DIFDIP = .FALSE.
C     SECDER = .FALSE.
C     NOPV   = .FALSE.
C     NODV   = .FALSE.
C
C     Make one-electron part of gradient
C
      CALL ONEDR1(WORK(KSTHMA),WORK(KDTAO),WORK(KFAO),WORK(KFACIN),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KNCENT),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),DUMMY,WORK(KFREE),LFREE,
     &            IPRGRD,.TRUE.,1,.TRUE.,.TRUE.,.FALSE.,NTBAS(0),
     &            NNBBASX,NNBBAST,DOINT,OMITVNUC,.FALSE.)
C
C     Restore original nuclear charges changed by ONECAP
C     for calculating corrected nuclear attraction integrals.
C
      IF (ONECAP .AND. INTV1C .EQ. 2) THEN
         DO I=1,NUCIND
            CHARGE(I) = WORK(KNQ-1+I)
         END DO
      END IF
C
C     Dispose temp. memory
C
      CALL MEMREL('RMOLGRD1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
      CALL GETTIM(CPU2,WALL2)
C
      CPUONE = CPU2 - CPU1
      WLLONE = WALL2 - WALL1
      WRITE(LUPRI,'(/A,F10.2,A/)') 
     &      'Time for one-electron part of gradient    :',
     &      CPUONE,' seconds.'
      CALL FLSHFO(LUPRI)
C
C
C     Two electron part of gradient
C     =============================
C
      CALL HEADER('RMOLGRD1: Two-electron part',-1)
      CALL GETTIM(CPU11,WALL11)
C
C
      CALL DZERO(GRADEE,MXCOOR)
      CALL DZERO(GRADER(1,1),MXCOOR)
      CALL DZERO(GRADER(1,2),MXCOOR)
      CALL DZERO(GRADER(1,3),MXCOOR)
C
C     Get contravariant density matrix, no folding
C
C
C     Get non-folded total contravariant density matrix
C     in symmetry distinct AO-basis
C
C     Matrices 1..4: DCao + DVao
C     Matrces  5..8: DVao
C
      CALL MEMGET('REAL',KDAO,4 * NFMAT * N2BBASX,WORK,KFREE,LFREE)
C
      IOFF = 0
      DO IMAT = 1, NFMAT
         DO I = 1, 4
            REWIND (LUMGRD)
            IF ( FNDLAB('D'//CHAR(IMAT+48)//'   CT'//CHAR(I+48),
     &           LUMGRD) ) THEN
               CALL READT(LUMGRD,N2BBASX,WORK(KDAO+IOFF))
            ELSE
               WRITE(LUPRI,'(/A/A)')
     $              '*** ERROR in RMOLGRD1 ***',
     $              'label "D'//CHAR(IMAT+48)//'   CT'//
     &              CHAR(I+48)//'"'//
     &              ' not found on file RMOLGRD'
               CALL QUIT('*** ERROR in RMOLGRD1 ***')
            END IF
            IOFF = IOFF + N2BBASX
         END DO
      END DO
      IF ( NFMAT .EQ. 2 ) THEN
C        ... add DVao to DCao
         CALL DAXPY(4*N2BBASX,D1,WORK(KDAO+4*N2BBASX),1,WORK(KDAO),1)
      END IF
C
      CALL GETTIM(CPU1,WALL1)
      IF (SCRGRD .GT. D0) THEN
C
C     Prepare screening
C     =================
C
C     CAll RTROINI with COMPAR true to get necessary gab-integrals.
C
         DO I = 1, 3*MXCENT
            DOPERT(I,1) = .TRUE.
            DOPERT(I,2) = .TRUE.
         END DO
         FTRONV = .TRUE.
         CALL RSETSYM(.TRUE.,RELGRD,.FALSE.,.FALSE.)
         CALL RTROINI(.TRUE.)
         CALL TROINV(WORK(KFREE),LFREE)
C
C        Read/generate undiff. GAB matrices
C
         N2GAB   = NSYMBL*NSYMBL
         CALL MEMGET ('REAL',KGABRAO, 4 * N2GAB, WORK, KFREE, LFREE )
C
         IJOB   = 0
         ITYPE  = 0
         IGTYP  = 1
         MAXDF = 0
         CALL GETGAB(IJOB,ITYPE,IGTYP,MAXDF,
     &        IPRINT,WORK(KGABRAO),WORK(KFREE),LFREE)
C
C        Read/generate diff. GAB matrices
C
         IJOB   = 0
         ITYPE  = 10
         IGTYP  = 1
         MAXDF  = 2
         CALL GETGAB(IJOB,ITYPE,IGTYP,MAXDF,
     &        IPRINT,WORK(KGABRAO+N2GAB),WORK(KFREE),LFREE)
         IF (IPRINT.GE.10) THEN
            DO I = 1, 3
               CALL HEADER('RMOLGRD1: GABRAO-matrix ('//CHRXYZ(I)//')',
     &              -1)
               CALL OUTPUT(WORK(KGABRAO+N2GAB*I),1,NSYMBL,1,NSYMBL,
     &                     NSYMBL,NSYMBL,-1,LUPRI)
            END DO
         END IF
C
C        Make reduced density matrices
C
         CALL MEMGET ('REAL',KDMRAO,4 * NFMAT * N2GAB, WORK,
     &        KFREE, LFREE )
C
         CALL MKDRAO(WORK(KDAO), WORK(KDMRAO), 4*NFMAT,
     &        WORK(KFREE),LFREE,IPRINT)
C
C        Restore translation/rotation invariance for gradient calc.
C
         DO I = 1, 3*MXCENT
            DOPERT(I,1) = .TRUE.
            DOPERT(I,2) = .TRUE.
         END DO
         FTRONV = .TRUE.
         CALL RSETSYM(.FALSE.,RELGRD,.FALSE.,.FALSE.)
         CALL RTROINI(.FALSE.)
         CALL TROINV(WORK(KFREE),LFREE)
C
      ELSE
C     ... no screening on gradient
         CALL MEMGET ('REAL',KGABRAO, 0, WORK, KFREE, LFREE )
         CALL MEMGET ('REAL',KDMRAO, 0, WORK, KFREE, LFREE )
      END IF
      CALL GETTIM(CPU2,WALL2)
      CPUGAB = CPU2 - CPU1
      WLLGAB = WALL2 - WALL1
      WRITE(LUPRI,'(/1X,A,F10.2,A/)') 
     &      'Time for GAB matrix                       :',
     &      CPUGAB,' seconds.'
C
C
C
      CALL HEADER('RMOLGRD1: LL-Two-electron part',-1)
      CALL CYESNO(STR,LBIT( IGRD_INTBUF, 1 ))
      WRITE(LUPRI,'(1X,2A/)') 'Calculate LL-Two-electron gradient: ',
     &                        STR
C
C
      IF ( LBIT( IGRD_INTBUF, 1 ) ) THEN
         IF ( DPRTIM(2,1) .GT. D0 ) THEN
            CALL TSTAMP(' ',LUPRI)
            WRITE(LUPRI,'(/1X,A,F10.2,A)') 
     &                  'Estimated CPU time:  ',DPRTIM(2,1),' seconds.'
            WRITE(LUPRI,'(1X,A,F10.2,A/)') 
     &                  'Estimated WALL time: ',DPRTIM(2,2),' seconds.'
         END IF
      ELSE
         CALL HEADER('RMOLGRD1: Skipping LL-Two-el.',-1)
         GOTO 100
      END IF
C
C
      CALL GETTIM(CPU1,WALL1)
C
C     LL-gradient
C     IREPDM = IDUMMY
      ITYPE  = 2
C     GMAT   = DUMMY
C     INDX   = IDUMMY
C     INDXAB = IDUMMY
      MAXDIF = 1
      JATOM  = 0
      NODV   = (NFMAT.EQ.1)
      NOPV   = .TRUE.
      NOCONT = .FALSE.
      FIRST  = .TRUE.
C     NPOS   = IDUMMY
      MTOTTK = NLRGBL*(NLRGBL+1)/2
      I2TYP  = 1
      INTCLASS  = 0
C     GABRAO = DUMMY
C     DMRSO  = DUMMY
C
      NEWBAS = .TRUE. 
      NDMAT   = NFMAT*4
      CALL DZERO(DNTSKP,3*14*3)
      CALL FLSHFO(LUPRI)
#if defined (VAR_MPI)
      IF (PARCAL) THEN
         CALL DIRAC_PARCTL( HERFCK_PAR )
         CALL MEMGET('INTE',KNPOS,MTOTTK,WORK,KFREE,LFREE)
         CALL HER_PARDRV(WORK(KFREE),LFREE,DUMMY,WORK(KDAO),NDMAT,
     &        IDUMMY,IFCTYP,ITYPE,MAXDIF,JATOM,NODV,NOPV,NOCONT,TKTIME,
     &        RTNTWO,FIRST,WORK(KNPOS),MTOTTK,I2TYP,ICEDIF,
     &        SCRGRD,WORK(KGABRAO),WORK(KDMRAO),DUMMY,DNTSKP)
      ELSE
#endif
         CALL TWOINT(WORK(KFREE),LFREE,DUMMY,WORK(KDAO),NDMAT,IDUMMY,
     &        IFCTYP,DUMMY,IDUMMY,IDUMMY,ITYPE,MAXDIF,JATOM,NODV,
     &        NOPV,NOCONT,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &        IPRNTD,RTNTWO,IDUMMY,I2TYP,ICEDIF,SCRGRD,
     &        WORK(KGABRAO),WORK(KDMRAO),DUMMY,DNTSKP,
     &        RELCAL,.false.,DUMMY,DUMMY)
#if defined (VAR_MPI)
      END IF
#endif
C
C     Symmetrization of GRADEE
C
      CALL MEMGET('REAL',KSKLTN,MXCOOR*2           ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAC  ,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDEX,3*NUCIND*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KPOINT,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL SKLEXP(WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),WORK(KPOINT),
     &               MAXDER,0)
C
C     Copy LL-GRADEE gradient to GRADER(1,1)
C
      CALL DCOPY(MXCOOR,GRADEE,1,GRADER(1,1),1)
      CALL DZERO(GRADEE,MXCOOR)
C
C     Memory deallocation
C
      CALL MEMREL('RMOLGRD1',WORK,KWORK,KSKLTN,KFREE,LFREE)
C
      CALL GETTIM(CPU2,WALL2)
C
      CPUTLL = CPU2 - CPU1
      WLLTLL = WALL2 - WALL1
      IF ( PARCAL ) THEN
         CPUTLL = MAX(CPUTLL,DMXCPT)
         WLLTLL = MAX(WLLTLL,DMXWLT)
      END IF
      IF (SCRGRD .GT. D0) CALL SCRGDS('SCR:LL',DNTSKP,IPRINT)
      WRITE(LUPRI,'(/1X,A,F10.2,A/)') 
     &      'Time for LL-two-electron part of gradient :',
     &      CPUTLL,' seconds.'
C
      CALL FLSHFO(LUPRI)
 100  CONTINUE
C
C
C     Calculate estimated LS- and SS- gradient and skip if
C     they seem unnecessary.
C
C
      CALL HEADER('RMOLGRD1: LS-Two-electron part',-1)
      IF ( DOTRCK .AND. LBIT ( IGRD_INTBUF, 2 ) ) THEN
         CALL ESTGLS(DOLS)
         IF ( .NOT. DOLS ) THEN
            IGRD_INTBUF = IBTAND ( 1, IGRD_INTBUF )
C           ... skip LS and SS gradient
         END IF
      END IF
      CALL CYESNO(STR, LBIT ( IGRD_INTBUF, 2 ) )
      WRITE(LUPRI,'(/1X,2A/)') 'Calculate LS-Two-electron gradient: ',
     &                        STR
C
      IF ( LBIT ( IGRD_INTBUF, 2 ) ) THEN
         IF ( DPRTIM(3,1) .GT. D0 ) THEN
            CALL TSTAMP(' ',LUPRI)
            WRITE(LUPRI,'(/1X,A,F10.2,A)') 
     &                  'Estimated CPU time:  ',DPRTIM(3,1),' seconds.'
            WRITE(LUPRI,'(1X,A,F10.2,A/)') 
     &                  'Estimated WALL time: ',DPRTIM(3,2),' seconds.'
         END IF
      ELSE
         CALL HEADER('RMOLGRD1: Skipping LS-Two-el.',-1)
         GOTO 200
      END IF
      CALL GETTIM(CPU1,WALL1)
C
C
C     LS-gradient
C     IREPDM = IDUMMY
      ITYPE  = 2
C     GMAT   = DUMMY
C     INDX   = IDUMMY
C     INDXAB = IDUMMY
      MAXDIF = 1
      JATOM  = 0
      NODV   = (NFMAT.EQ.1)
      NOPV   = .TRUE.
      NOCONT = .FALSE.
      FIRST  = .TRUE.
C     NPOS   = IDUMMY
      MTOTTK = NLRGBL*(NLRGBL+1)/2
      I2TYP  = 2
      INTCLASS  = 1
C     GABRAO = DUMMY
C     DMRSO  = DUMMY
C
      NEWBAS = .TRUE.
      NDMAT   = NFMAT*4
      CALL DZERO(DNTSKP,3*14*3)
      CALL FLSHFO(LUPRI)
#if defined (VAR_MPI)
      IF (PARCAL) THEN
         CALL DIRAC_PARCTL( HERFCK_PAR )
         CALL MEMGET('INTE',KNPOS,MTOTTK,WORK,KFREE,LFREE)
         CALL HER_PARDRV(WORK(KFREE),LFREE,DUMMY,WORK(KDAO),NDMAT,
     &        IDUMMY,IFCTYP,ITYPE,MAXDIF,JATOM,NODV,NOPV,NOCONT,TKTIME,
     &        RTNTWO,FIRST,WORK(KNPOS),MTOTTK,I2TYP,ICEDIF,
     &        SCRGRD,WORK(KGABRAO),WORK(KDMRAO),DUMMY,DNTSKP)
      ELSE
#endif
         CALL TWOINT(WORK(KFREE),LFREE,DUMMY,WORK(KDAO),NDMAT,IDUMMY,
     &        IFCTYP,DUMMY,IDUMMY,IDUMMY,ITYPE,MAXDIF,JATOM,NODV,
     &        NOPV,NOCONT,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &        IPRNTD,RTNTWO,IDUMMY,I2TYP,ICEDIF,SCRGRD,
     &        WORK(KGABRAO),WORK(KDMRAO),DUMMY,DNTSKP,
     &        RELCAL,.false.,DUMMY,DUMMY)
#if defined (VAR_MPI)
      END IF
#endif
      CALL FLSHFO(LUPRI)
C
C     Symmetrization of GRADEE
C
      CALL MEMGET('REAL',KSKLTN,MXCOOR*2           ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAC  ,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDEX,3*NUCIND*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KPOINT,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL SKLEXP(WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),WORK(KPOINT),
     &               MAXDER,0)
C
C     Copy LS-GRADEE gradient to GRADER(1,2)
C
      CALL DCOPY(MXCOOR,GRADEE,1,GRADER(1,2),1)
      CALL DZERO(GRADEE,MXCOOR)
C
C     Dispose memory
C
      CALL MEMREL('RMOLGRD1',WORK,KWORK,KSKLTN,KFREE,LFREE)
C
      CALL GETTIM(CPU2,WALL2)
C
      CPUTLS = CPU2 - CPU1
      WLLTLS = WALL2 - WALL1
      IF ( PARCAL ) THEN
         CPUTLS = MAX(CPUTLS,DMXCPT)
         WLLTLS = MAX(WLLTLS,DMXWLT)
      END IF
      IF (SCRGRD .GT. D0) CALL SCRGDS('SCR:LS',DNTSKP,IPRINT)
      WRITE(LUPRI,'(/1X,A,F10.2,A/)') 
     &      'Time for LS-two-electron part of gradient :',
     &      CPUTLS,' seconds.'
C
C
C
C
      CALL FLSHFO(LUPRI)
 200  CONTINUE
C
      CALL HEADER('RMOLGRD1: SS-Two-electron part',-1)
      IF ( LBIT( IGRD_INTBUF, 3 ) .AND.
     &     .NOT. LBIT( IGRD_INTBUF, 2 ) ) THEN
         WRITE(LUPRI,'(/1X,A/)') 'WARNING: SS-gradient .NOT. '//
     &       'calculated since LS-gradient is not calculated.'
         IGRD_INTBUF = IBTAND ( IGRD_INTBUF, 1+2 )
      END IF
C
      IF ( DOTRCK .AND. LBIT ( IGRD_INTBUF, 3 )
     &     .AND. .NOT. DOLVC ) THEN
CHJAaJ July 2001: Idea: this can be done much better based on GRADLV !!! TODO
         CALL ESTGSS(DOSS)
         IF ( .NOT. DOSS ) THEN
            IGRD_INTBUF = IBTAND ( 1+2, IGRD_INTBUF )
         END IF
      END IF
C
      CALL CYESNO(STR, LBIT( IGRD_INTBUF, 3 ) )
      WRITE(LUPRI,'(1X,2A/)') 'Calculate SS-Two-electron gradient: ',
     &                        STR
C
      IF ( LBIT( IGRD_INTBUF, 3 ) ) THEN
         IF ( DPRTIM(4,1) .GT. D0 ) THEN
            CALL TSTAMP(' ',LUPRI)
            WRITE(LUPRI,'(/1X,A,F10.2,A)') 
     &                  'Estimated CPU time:  ',DPRTIM(4,1),' seconds.'
            WRITE(LUPRI,'(1X,A,F10.2,A/)') 
     &                  'Estimated WALL time: ',DPRTIM(4,2),' seconds.'
         END IF
      ELSE
         CALL HEADER('RMOLGRD1: Skipping SS-Two-el.',-1)
         GOTO 300
      END IF
      CALL GETTIM(CPU1,WALL1)
C
C
C     SS-gradient
C
C     IREPDM = IDUMMY
      ITYPE  = 2
C     GMAT   = DUMMY
C     INDX   = IDUMMY
C     INDXAB = IDUMMY
      MAXDIF = 1
      JATOM  = 0
      NODV   = (NFMAT.EQ.1)
      NOPV   = .TRUE.
      NOCONT = .FALSE.
      FIRST  = .TRUE.
C     NPOS   = IDUMMY
      MTOTTK = NSMLBL*(NSMLBL+1)/2
      I2TYP  = 3
      INTCLASS  = 2
C     GABRAO = DUMMY
C     DMRSO  = DUMMY
C
      NEWBAS = .TRUE.
      NDMAT   = NFMAT*4
      CALL DZERO(DNTSKP,3*14*3)
      CALL FLSHFO(LUPRI)
#if defined (VAR_MPI)
      IF (PARCAL) THEN
         CALL DIRAC_PARCTL( HERFCK_PAR )
         CALL MEMGET('INTE',KNPOS,MTOTTK,WORK,KFREE,LFREE)
         CALL HER_PARDRV(WORK(KFREE),LFREE,DUMMY,WORK(KDAO),NDMAT,
     &        IDUMMY,IFCTYP,ITYPE,MAXDIF,JATOM,NODV,NOPV,NOCONT,TKTIME,
     &        RTNTWO,FIRST,WORK(KNPOS),MTOTTK,I2TYP,ICEDIF,
     &        SCRGRD,WORK(KGABRAO),WORK(KDMRAO),DUMMY,DNTSKP)
      ELSE
#endif
         CALL TWOINT(WORK(KFREE),LFREE,DUMMY,WORK(KDAO),NDMAT,IDUMMY,
     &        IFCTYP,DUMMY,IDUMMY,IDUMMY,ITYPE,MAXDIF,JATOM,NODV,
     &        NOPV,NOCONT,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &        IPRNTD,RTNTWO,IDUMMY,I2TYP,ICEDIF,SCRGRD,
     &        WORK(KGABRAO),WORK(KDMRAO),DUMMY,DNTSKP,
     &        RELCAL,.false.,DUMMY,DUMMY)
#if defined (VAR_MPI)
      END IF
#endif
      CALL FLSHFO(LUPRI)
C
C     Dispose memory
C
      CALL MEMREL('RMOLGRD1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Symmetrization of GRADEE
C
      CALL MEMGET('REAL',KSKLTN,MXCOOR*2           ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAC  ,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDEX,3*NUCIND*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KPOINT,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL SKLEXP(WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),WORK(KPOINT),
     &               MAXDER,0)
C
C     Copy SS-GRADEE gradient to GRADER(1,3)
C
      CALL DCOPY(MXCOOR,GRADEE,1,GRADER(1,3),1)
      CALL DZERO(GRADEE,MXCOOR)
C
C
C     End of gradient calculation
C     ===========================
C
C     Memory deallocation
C
      CALL MEMREL('RMOLGRD1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL GETTIM(CPU2,WALL2)
C
      CPUTSS = CPU2 - CPU1
      WLLTSS = WALL2 - WALL1
      IF ( PARCAL ) THEN
         CPUTSS = MAX(CPUTSS,DMXCPT)
         WLLTSS = MAX(WLLTSS,DMXWLT)
      END IF
      IF (SCRGRD .GT. D0) CALL SCRGDS('SCR:SS',DNTSKP,IPRINT)
      WRITE(LUPRI,'(/1X,A,F10.2,A/)') 
     &      'Time for SS-two-electron part of gradient :',
     &      CPUTSS,' seconds.'
C
      CALL FLSHFO(LUPRI)
 300  CONTINUE
C
C     Output of timings
C
      CALL HEADER('Timings for RMOLGRD1',-1)
      WRITE(LUPRI,'(/1X,A,F10.2,A,F10.2,A)') 
     &   'One-electron part of gradient    CPU  :',CPUONE,'s    WALL :',
     &   WLLONE,'s'
      WRITE(LUPRI,'(/1X,A,F10.2,A,F10.2,A,F10.2,A)') 
     &   'GAB matrix for screening         CPU  :',CPUGAB,'s    WALL :',
     &   WLLGAB,'s'
      IF ( LBIT( IGRD_INTBUF, 1 ) ) THEN
         WRITE(LUPRI,'(1X,A,F10.2,A,F10.2,A)') 
     &   'LL-Two-electron part of gradient CPU  :',CPUTLL,'s    WALL :',
     &   WLLTLL,'s'
      ELSE
         WRITE(LUPRI,'(1X,A)') 
     &   'LL-Two-electron part SKIPPED!'
      END IF
      IF ( LBIT( IGRD_INTBUF, 2 ) ) THEN
         WRITE(LUPRI,'(1X,A,F10.2,A,F10.2,A)') 
     &   'LS-Two-electron part of gradient CPU  :',CPUTLS,'s    WALL :',
     &   WLLTLS,'s'
      ELSE
         WRITE(LUPRI,'(1X,A)') 
     &    'LS-Two-electron part SKIPPED!'
      END IF
      IF ( LBIT( IGRD_INTBUF, 3 ) ) THEN
         WRITE(LUPRI,'(1X,A,F10.2,A,F10.2,A)') 
     &   'SS-Two-electron part of gradient CPU  :',CPUTSS,'s    WALL :',
     &   WLLTSS,'s'
      ELSE
         WRITE(LUPRI,'(1X,A)') 
     &   'SS-Two-electron part SKIPPED!'
      END IF
      WRITE(LUPRI,'(1X,A,F10.2,A,F10.2,A)') 
     &    'Total two-electron part of grad. CPU  :',
     &      CPUTSS+CPUTLS+CPUTLL+CPUGAB,'s    WALL :',
     &      WLLTSS+WLLTLS+WLLTLL+WLLGAB,'s'
      CALL GETTIM(CPUTO2,WLLTO2)
      CPUTOT = CPUTO2 - CPUTO1
      WLLTOT = WLLTO2 - WLLTO1
      WRITE(LUPRI,'(1X,A,F10.2,A,F10.2,A/)') 
     &   'Total time in RMOLGRD1           CPU  :',CPUTOT,'s    WALL :',
     &    WLLTOT,'s'
C
C     Save timings for later use
C
      DPRTIM(1,1) = CPUONE
      DPRTIM(2,1) = CPUTLL
      DPRTIM(3,1) = CPUTLS
      DPRTIM(4,1) = CPUTSS
      DPRTIM(1,2) = WLLONE
      DPRTIM(2,2) = WLLTLL
      DPRTIM(3,2) = WLLTLS
      DPRTIM(4,2) = WLLTSS
C
C     Finished gradient, reset RELGRD to .FALSE.
C     also reset MOLGRD in abainf.h to false
C
      RELGRD = .FALSE.
      CALL RSETSYM(.FALSE.,RELGRD,.FALSE.,.FALSE.)
C
      CALL QEXIT('RMOLGRD1')
      CALL FLSHFO(LUPRI)
C
      RETURN
      END
C
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsetsym */
      SUBROUTINE RSETSYM(ALLSYM,RELGRD,RELHES,RELDIPD)
C*****************************************************************************
C
C     Temporary subroutine: sets DOSYM-array from abainf.h
C
C     jth, 1997/24/04; revised 010626-hjaaj: also set MOLGRD,MOLHES,DIPDER
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "abainf.h"
      LOGICAL ALLSYM,RELGRD,RELHES,RELDIPD
      DOSYM(1) = .TRUE.
      DO I=2,8
         DOSYM(I) = ALLSYM
      END DO
      MOLGRD = RELGRD
      MOLHES = RELHES
      DIPDER = RELDIPD
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dsym1r */
      SUBROUTINE DSYM1R(ASO,AAO,NBAST,IREPDM,IPRINT)
C*****************************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Slightly changed to DSYM1R by jth 1997/05/09
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C*****************************************************************************
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
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine DSYM1R',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
C      CALL DZERO(AAO,NBAST*NBAST)
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
Cjth - insert
                     MAB    = IBTOR(MULA,MULB)
                     KAB    = IBTAND(MULA,MULB)
                     HKAB   = FMULT(KAB)
Cjth - endinsert
                     INDB   = KSTRT(IB) + NB - KHKTB
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        AAO(INDA,INDB) = AAO(INDA,INDB)
     &                                 + HKAB*FACA*FACB*ASO(I,J)
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
         CALL HEADER('DSYM1R.INPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         CALL HEADER('DSYM1R.OUTPUT:Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grdout */
      SUBROUTINE GRDOUT(WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     (c) 1997/05/28 by J. Thyssen
C
C
C     Last revision: 1997/06/11 jth
C
C*****************************************************************************

      use dirac_cfg

#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
#include "taymol.h"
#include "dummy.h"
C
#include "dcbpsi.h"
#include "dcbgen.h"
#include "dcbgrd.h"
#include "dcbham.h"
C
      LOGICAL LBIT
      CHARACTER*(20) LABEL(3)
      DIMENSION WORK(*)
C
      CALL QENTER('GRDOUT')
C
#include "memint.h"
C
      CALL MEMGET('REAL',KCSTRA,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCTRA,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
C
      CALL TITLER('Final output from PSIGRD','*',110)
C
C     Sum up gradient
C
      DO I = 1,3*NUCDEP
         GRADER(I,0) = GRADER(I,1) + GRADER(I,2) + GRADER(I,3)
         GRADRO(I,0) = GRADRO(I,1) + GRADRO(I,2)
         GRADNU(I,0) = GRADNU(I,1) + GRADNU(I,2)
         GRADKN(I,0) = GRADKN(I,1) + GRADKN(I,2)
      END DO
C
      IF ( IPRINT .GE. 1 ) THEN

        CALL HEADER('Kinetic energy integral gradient',-1)
        LABEL(1) = 'Momentum'
        LABEL(2) = 'Beta mtx.'
        CALL PRIGDR(GRADKN,2,WORK(KCSTRA),WORK(KSCTRA),LABEL,' ')

        CALL HEADER('Nuclear attraction integral gradient',-1)
        LABEL(1) = 'LL'
        LABEL(2) = 'SS'
        CALL PRIGDR(GRADNU,2,WORK(KCSTRA),WORK(KSCTRA),LABEL,' ')

        CALL HEADER('Reorthonormalization gradient',-1)
        LABEL(1) = 'LL'
        LABEL(2) = 'SS'
        CALL PRIGDR(GRADRO,2,WORK(KCSTRA),WORK(KSCTRA),LABEL,' ')

        CALL HEADER('Two-electron integral gradient',-1)
        LABEL(1) = 'LL'
        LABEL(2) = 'LS'
        LABEL(3) = 'SS'
        CALL PRIGDR(GRADER,3,WORK(KCSTRA),WORK(KSCTRA),LABEL,' ')

        if (dirac_cfg_dft_calculation) then
          call header('XC contribution to molecular gradient',-1)
          call prigdr(gradxc, 0, work(kcstra), work(ksctra), label, ' ')
        end if
C
C       "SS-integral contribution by classical repulsion of
C        small component atomic charges"-gradient
C       ===================================================
C
        IF (ONECAP .AND. INTV1C.EQ.2) THEN
           CALL HEADER(
     &      'One-center model 2 SS correction to integral gradient',-1)
           CALL PRIGDR(GRADLV,0,WORK(KCSTRA),WORK(KSCTRA),LABEL,' ')
        ELSE IF ( LBIT( IGRD_INTBUF, 2 ) .AND. LBIT( IGRD_INTFLG, 3 ) )
     &     THEN
           CALL HEADE3('"SS-integral contribution by classical ',
     &                 'repulsion of small component atomic ',
     &                 'charges"-gradient',-1)
           IF ( LBIT( IGRD_INTBUF, 3 ) ) THEN
              WRITE(LUPRI,'(/A)') 'SCC gradient is NOT used, since '//
     &             'SS-gradient has been calculated.'
           ELSE IF (.NOT.DOLVC) THEN
              WRITE(LUPRI,'(/A)') 'SCC gradient is used, since '//
     &             'SS-gradient was NOT calculated for this geometry.'
           END IF 
           CALL PRIGDR(GRADLV,0,WORK(KCSTRA),WORK(KSCTRA),LABEL,' ')
        END IF
C       Nuclear repulsion
C       =================
C
        CALL HEADER('Nuclear repulsion gradient',-1)
        CALL PRIGDR(GRADNN,0,WORK(KCSTRA),WORK(KSCTRA),LABEL,' ')
C
C
      END IF
C
C
      CALL ZERGRD
      CALL ADDGRD(GRADKN(1,0))
      CALL ADDGRD(GRADNU(1,0))
      CALL ADDGRD(GRADRO(1,0))
      CALL ADDGRD(GRADER(1,0))
      CALL ADDGRD(GRADNN)
      if (dirac_cfg_dft_calculation) then
        call addgrd(gradxc)
      end if
      IF ( (ONECAP .AND. INTV1C.EQ.2) .OR. 
     &   ( LBIT( IGRD_INTBUF, 2 ) .AND. LBIT( IGRD_INTFLG, 3 ) .AND.
     &     .NOT. LBIT ( IGRD_INTBUF, 3 ) ) ) THEN
         CALL ADDGRD(GRADLV)
      END IF
      CALL HEADE2('Total molecular gradient','@',-1)
      CALL PRIGDR(GRDMOL,0,WORK(KCSTRA),WORK(KSCTRA),LABEL,'@')
C
      IF (IPRINT .GE. 3) THEN
C
C        Print norm of molecular gradient
C
         CALL GRDNRM(GRDMOL,WORK(KCSTRA),WORK(KSCTRA),IPRINT)
      ENDIF
C
      CALL MEMREL('GRDOUT',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GRDOUT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck estgls */
      SUBROUTINE ESTGLS(DOLS)
C*****************************************************************************
C
C     (c) 1997/06/27 by J. Thyssen
C
C     On basis of one-electron gradient,g(1), and two-electron LL-gradient, 
C     g(2,LL), estimate if two-electron LS- gradient is necessary. 
C
C     If ||g(2,LS)|| > 0.1 ||g(1)+g(2,LL)|| then we calculate g(2,LS).
C     As an estimate of ||g(2,LS)|| we use 0.005 ||g(2,LL)|| based
C     on empirical facts.
C     Combined: DOLS if ||g(2,LL)|| > 20. ||g(1)+g(2,LL)||
C
C     Since the LS-gradient is based on (LL|SS) or (SS|LL) integrals
C     it scales with 1/(c^2). 
C
C     Last revision: 1997/07/02 jth
C
C*****************************************************************************
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "consts.h"
C
      DIMENSION TGRAD(MXCOOR)
      LOGICAL DOLS
      CHARACTER*3 STR
#include "nuclei.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbgrd.h"
C
      CALL QENTER('ESTGLS')
C
      DOLS = .FALSE.
C
      DO I = 1, 3*NUCDEP
C
C        Calculate "total" gradient so far
C
         TGRAD(I)  = GRADNN(I)   + 
     &               GRADRO(I,1) + GRADRO(I,2) +
     &               GRADKN(I,1) + GRADKN(I,2) +
     &               GRADNU(I,1) + GRADNU(I,2) +
     &               GRADER(I,1)
C
      END DO
C
C     Calculate norms
C
      TOTNRM = DDOT(3*NUCDEP,TGRAD,      1,TGRAD,      1)
      TOTNRM = SQRT(TOTNRM)
      DLLNRM = DDOT(3*NUCDEP,GRADER(1,1),1,GRADER(1,1),1)
      DLLNRM = SQRT(DLLNRM)
C
C     If ||g(2-el,LL)|| / ||g(1-el)+g(2-el,LL)|| > 20 then calculate LS-gradient
C
      CFAC = (CVAL / CVEL) ** 2
C     user may have changed CVAL, the velocity of light (!)
      IF ( TOTNRM .GT. D0 ) THEN
         DL = DLLNRM / TOTNRM
         DR = 20.D0 * CFAC
         IF ( DL .GE. DR ) THEN
            DOLS = .TRUE.
         END IF
      ELSE
         DOLS = .TRUE.
      END IF
      WRITE(LUPRI,'(/A)') ' Output from ESTGLS'
      WRITE(LUPRI,1000) ' Norm of g(1)+g(2,LL): ',TOTNRM
      WRITE(LUPRI,1000) ' Norm of g(2,LL)     : ',DLLNRM
      WRITE(LUPRI,1000) ' Est. norm of g(2,LS): ',0.005D0*DLLNRM
      CALL CYESNO(STR,DOLS)
      WRITE(LUPRI,'(A,A3/)') ' Calculate g(2,LS)   : ',STR
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('ESTGLS')
      RETURN
 1000 FORMAT(A,F16.9)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck estgss */
      SUBROUTINE ESTGSS(DOSS)
C*****************************************************************************
C
C     (c) 1997/07/02 by J. Thyssen
C
C     Do something similar as in ESTGLS
C
C     Last revision: 1997/07/02 jth
C
C*****************************************************************************
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "consts.h"
C
      PARAMETER(DP1 = 0.1D00)
      DIMENSION TGRAD(MXCOOR)
      LOGICAL DOSS
      CHARACTER*3 STR
#include "nuclei.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbgrd.h"
#include "optinf.h"
C
      CALL QENTER('ESTGSS')
C
      DOSS = .FALSE.
C
      DO I = 1, 3*NUCDEP
C
C        Calculate "total" gradient so far
C
         TGRAD(I)  = GRADNN(I)   + 
     &               GRADRO(I,1) + GRADRO(I,2) +
     &               GRADKN(I,1) + GRADKN(I,2) +
     &               GRADNU(I,1) + GRADNU(I,2) +
     &               GRADER(I,1) + GRADER(I,2)
C
      END DO
C
C     Calculate norms
C
      TOTNRM = DDOT(3*NUCDEP,TGRAD,      1,TGRAD,      1)
      TOTNRM = SQRT(TOTNRM)
      DLLNRM = DDOT(3*NUCDEP,GRADER(1,1),1,GRADER(1,1),1)
      DLLNRM = SQRT(DLLNRM)
      DLSNRM = DDOT(3*NUCDEP,GRADER(1,2),1,GRADER(1,2),1)
      DLSNRM = SQRT(DLSNRM)
C
C     As estimate for ||g(2,SS)|| use ||g(2,LS)||^2 / ||g(2,LL)||
C
CHJ   CFAC = (CVAL / CVEL) ** 4
C     user may have changed CVAL, the velocity of light (!)
C     hjaaj aug 2001: no, this factor is already in DLSNRM !
      IF ( TOTNRM .GT. D0 ) THEN
         DR = ( DLSNRM ** 2 ) / DLLNRM
CHJ      DL = CFAC * DP1 * MAX( TOTNRM, GRDTHR )
         DL = DP1 * MAX( TOTNRM, GRDTHR )
         IF ( DR .GE. DL ) THEN
            DOSS = .TRUE.
         END IF
      ELSE
         DOSS = .TRUE.
      END IF
C
      WRITE(LUPRI,'(/1X,A)') 'Output from ESTGSS'
      WRITE(LUPRI,1000) 'Norm of g(1)+g(2,LL)+g(2,LS): ',TOTNRM
      WRITE(LUPRI,1000) 'Norm of g(2,LL)             : ',DLLNRM
      WRITE(LUPRI,1000) 'Norm of g(2,LS)             : ',DLSNRM
      WRITE(LUPRI,1000) 'Estimated norm of g(2,SS)   : ',DR    
      CALL CYESNO(STR,DOSS)
      WRITE(LUPRI,'(1X,A,A3/)') 'Calculate g(2,SS)           : ',
     &                              STR
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('ESTGSS')
      RETURN
 1000 FORMAT(1X,A,F16.9)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck cyesno */
      SUBROUTINE CYESNO(STR,BOOL)
C*****************************************************************************
C
C     (c) 1997/07/08 by J. Thyssen
C
C     Return yes or no depending on BOOL
C
C     Last revision: 1997/07/08 jth
C
C*****************************************************************************
#include "implicit.h"
      CHARACTER*(*) STR
      LOGICAL BOOL
C
      IF ( BOOL ) THEN
         STR = 'yes'
      ELSE
         STR = 'no'
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck delfile */
      SUBROUTINE DELFILE(IUNIT,FILENAME)
C*****************************************************************************
C
C     (c) 1997/07/09 by J. Thyssen
C
C     Delete file FILENAME
C
C     Last revision: 1997/07/09 jth
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) FILENAME
      LOGICAL EX
C
      CALL QENTER('DELFILE')
C
      INQUIRE(FILE=FILENAME,EXIST=EX)
      IF ( EX ) THEN
         CALL OPNFIL(IUNIT,FILENAME,'UNKNOWN','DELFILE')
         CLOSE(IUNIT,STATUS='DELETE')
         INQUIRE(FILE=FILENAME,EXIST=EX)
         IF ( EX ) THEN
            WRITE(LUPRI,'(//1X,3A//)') 'WARNING: ',
     &         FILENAME,' could .NOT. be deleted!' 
            CALL QUIT('DELFILE: file could not be deleted.')
         END IF
      END IF
C
C
      CALL QEXIT('DELFILE')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lvcgrd */
      SUBROUTINE LVCGRD(IPRINT,SCQ,WORK,LWORK)
C*****************************************************************************
C
C     (c) 1997/09/23 by J. Thyssen
C
C     "SS-integral contribution by classical repulsion of
C      small component atomic charges"-gradient
C
C     Last revision: 1997/09/23 jth
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0,
     *           THREE = 3.0D0)
#include "nuclei.h"
#include "frame.h"
#include "dorps.h"
#include "symmet.h"
C
      DIMENSION SCQ(NUCIND),WORK(LWORK)
#include "dcbgrd.h"
#include "dcbham.h"
#include "ibtfun.h"
C
C
      CALL QENTER('LVCGRD')
C
      IF (ONECAP .AND. INTV1C .EQ. 2) THEN
         DO I = 1,NUCIND
            SCQ(I) = SCQ(I) + CHARGE(I)
         END DO
      END IF
C
C     This is copied from ABACUS/NUCREP
C
      NCDEP3 = 3*NUCDEP
      CALL DZERO(GRADLV,NCDEP3)
C
C
C     Run over symmetry-independent nuclei A
C
      DO 200 NCENTA = 1, NUCIND
         CHARGA = SCQ(NCENTA)
         IF (ABS(CHARGA) .GT. ZERO) THEN
            NAX    = 3*NCENTA - 2
            NAY    = 3*NCENTA - 1
            NAZ    = 3*NCENTA
            CORDXA = CORD(1,NCENTA)
            CORDYA = CORD(2,NCENTA)
            CORDZA = CORD(3,NCENTA)
            MULA   = ISTBNU(NCENTA)
C
C           Run over symmetry-independent nuclei B
C
            DO 300 NCENTB =  1, NCENTA
               CHARGB = SCQ(NCENTB)
               IF (ABS(CHARGB) .GT. ZERO) THEN
                  NBX    = 3*NCENTB - 2
                  NBY    = 3*NCENTB - 1
                  NBZ    = 3*NCENTB
                  CORBX0 = CORD(1,NCENTB)
                  CORBY0 = CORD(2,NCENTB)
                  CORBZ0 = CORD(3,NCENTB)
                  MULB   = ISTBNU(NCENTB)
C
                  MAB    = IBTOR (MULA,MULB)
                  KAB    = IBTAND(MULA,MULB)
                  HKAB   = FMULT(KAB)
                  CROSS = ONE
                  IF (NCENTA .EQ. NCENTB) THEN
                     HKAB = HALF*HKAB
                     CROSS = TWO
                  END IF
C
C                 Run over symmetry independent charge distributions
C
                  DO 400 ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MAB) .EQ. 0) THEN
                     ICENTA = NUCNUM(NCENTA,1)
                     ICENTB = NUCNUM(NCENTB,ISYMOP+1)
                  IF (ICENTA .EQ. ICENTB) GO TO 400
C
                     SBX = PT(IBTAND(ISYMAX(1,1),ISYMOP))
                     SBY = PT(IBTAND(ISYMAX(2,1),ISYMOP))
                     SBZ = PT(IBTAND(ISYMAX(3,1),ISYMOP))
                     XAB = CORDXA - SBX*CORBX0
                     YAB = CORDYA - SBY*CORBY0
                     ZAB = CORDZA - SBZ*CORBZ0
C
                     XAB2   = XAB*XAB
                     YAB2   = YAB*YAB
                     ZAB2   = ZAB*ZAB
                     RAB2   = XAB2 + YAB2 + ZAB2
                     RAB1   = SQRT(RAB2)
                     ZZR3IN = - HKAB*CHARGA*CHARGB/(RAB1*RAB2)
                     VNUCX  = XAB*ZZR3IN
                     VNUCY  = YAB*ZZR3IN
                     VNUCZ  = ZAB*ZZR3IN
C
C                    ********************
C                    ***** Gradient *****
C                    ********************
C
C                    Totally symmetric contributions only
C
                     IF (DOREPS(0)) THEN
                        IAX  = IPTCNT(NAX,0,1)
                        IAY  = IPTCNT(NAY,0,1)
                        IAZ  = IPTCNT(NAZ,0,1)
                        IBX  = IPTCNT(NBX,0,1)
                        IBY  = IPTCNT(NBY,0,1)
                        IBZ  = IPTCNT(NBZ,0,1)
                        IF (IAX.NE.0) GRADLV(IAX) =GRADLV(IAX)+VNUCX
                        IF (IAY.NE.0) GRADLV(IAY) =GRADLV(IAY)+VNUCY
                        IF (IAZ.NE.0) GRADLV(IAZ) =GRADLV(IAZ)+VNUCZ
                        IF (IBX.NE.0) GRADLV(IBX) =GRADLV(IBX)-SBX*VNUCX
                        IF (IBY.NE.0) GRADLV(IBY) =GRADLV(IBY)-SBY*VNUCY
                        IF (IBZ.NE.0) GRADLV(IBZ) =GRADLV(IBZ)-SBZ*VNUCZ
                     END IF
                  END IF
  400             CONTINUE
               END IF
  300       CONTINUE
         END IF
  200 CONTINUE
C
C      Restore small component charges
C
      IF (ONECAP .AND. INTV1C .EQ. 2) THEN
         DO I = 1,NUCIND
            SCQ(I) = SCQ(I) - CHARGE(I)
         END DO
      END IF
      CALL QEXIT('LVCGRD')
C
      RETURN
      END
C
C
      SUBROUTINE HEADE3(HEAD1,HEAD2,HEAD3,IN)
#include "implicit.h"
      CHARACTER*(*) HEAD1,HEAD2,HEAD3
#include "priunit.h"
C
      IF (IN .GE. 0) THEN
         I1 = IN + 1
         I2 = I1
         I3 = I1
      ELSE
         L1 = LEN(HEAD1)
         L2 = LEN(HEAD2)
         L3 = LEN(HEAD3)
         LM = MAX(L1,L2,L3)
         I1 = (72 - L1) /2 + 1
         I2 = (72 - L2) /2 + 1
         I3 = (72 - L3) /2 + 1
         IM = (72 - LM) /2 + 1
      END IF
      WRITE(LUPRI, '(//,80A)') (' ',I=1,I1), HEAD1
      WRITE(LUPRI, '(80A)')    (' ',I=1,I2), HEAD2
      WRITE(LUPRI, '(80A)')    (' ',I=1,I3), HEAD3
      WRITE(LUPRI, '(80A)') (' ',I=1,IM), ('-',I=1,LM)
      WRITE(LUPRI, '()')
      RETURN
      END
C  /* Deck prigdr */
      SUBROUTINE PRIGDR(GRAD,ICOMPO,CSTRA,SCTRA,LABEL,CFIRST)
C
C     Written by Joern Thyssen 1997/06/20
C     Last revision: 1997/06/20 jth
C
C     Print relativistic gradient
C     Based on PRIGRD
C     GRAD is asummed to be an array GRAD(MXCOOR,ICOMPO), where
C     (1) GRAD(0) is total gradient
C     (2) GRAD(1) is LABEL(1) component of gradient
C     (3) GRAD(2) is LABEL(2) c... etc.
C     CFIRST is the first character printed on each line
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "abainf.h"
#include "nuclei.h"
#include "symmet.h"
      DIMENSION GRAD(MXCOOR,0:ICOMPO), CGRAD(MXCOOR,0:3), 
     &          CSTRA(*), SCTRA(*)
      CHARACTER*(*) LABEL(ICOMPO)
      CHARACTER*10  TLABEL
      CHARACTER*1   CFIRST
C
      DIMENSION DNORM(0:3)
C
C
      IF (.NOT.DOSYM(1)) THEN
         WRITE (LUPRI,'(2X,A)')
     &  ' Gradient not calculated - '//
     &  ' totally symmetric distortions not requested.'
      ELSE IF (MAXREP .EQ. 0) THEN
         IOFF = 0
         DO 100 IATOM = 1, NUCDEP
            TLABEL = 'Total     '
            IF ( ICOMPO .EQ. 0 ) TLABEL = '          '
            WRITE (LUPRI,1000) CFIRST,NAMDEP(IATOM),TLABEL,
     &                         (GRAD(IOFF+J,0),J=1,3)
            DO K = 1, ICOMPO
               WRITE (LUPRI,1000) CFIRST,'      ',LABEL(K),
     &                            (GRAD(IOFF+J,K),J=1,3)
            END DO
            IF ( ICOMPO .GT. 0 ) WRITE(LUPRI,'()')
            IOFF = IOFF + 3
  100    CONTINUE
      ELSE
         TLABEL = 'Total     '
         IF ( ICOMPO .EQ. 0 ) TLABEL = '          '
         DO 200 I = 1, NCRREP(0,1)
            WRITE (LUPRI,'(16X,A6,3X,A10,F17.10)') 
     &        NAMEX(IPTCOR(I,1)),TLABEL,GRAD(I,0)
            DO J = 1, ICOMPO
               WRITE (LUPRI,'(16X,A6,3X,A10,F17.10)') 
     &           '      ',LABEL(J),GRAD(I,J)
            END DO
  200    CONTINUE
         WRITE (LUPRI,'(//)')
         DO K = 0, ICOMPO
            CALL TRAGRD(GRAD(1,K),CGRAD(1,K),CSTRA,SCTRA,NCRREP(0,1),
     &                  3*NUCDEP)
         END DO
         IOFF = 0
         DO 300 IATOM = 1, NUCDEP
            WRITE (LUPRI,1000) CFIRST,NAMDEP(IATOM),TLABEL,
     &                         (CGRAD(IOFF+J,0),J=1,3)
            DO K = 1, ICOMPO
               WRITE (LUPRI,1000) CFIRST,'      ',LABEL(K),
     &                            (CGRAD(IOFF+J,K),J=1,3)
            END DO
            IF ( ICOMPO .GT. 0 ) WRITE(LUPRI,'()')
            IOFF = IOFF + 3
  300    CONTINUE
C
      END IF
C 
C     Calculate weight statistics
C
      IF ( ICOMPO .GT. 0 ) THEN
         DNORM(0) = 0.0D00
         DO K = 1,ICOMPO
            DNORM(K) = DDOT(3*NUCDEP,GRAD(1,K),1,GRAD(1,K),1)
            DNORM(K) = SQRT(DNORM(K))
            DNORM(0) = DNORM(0) + DNORM(K)
         END DO
         WRITE(LUPRI,'(A)') 'Weights:'
         DO K = 1,ICOMPO
            WRITE(LUPRI,'(1X,A10,5X,F13.9,A)') LABEL(K),
     &            1.0D02 * DNORM(K) / DNORM(0),'%'
         END DO
      END IF
      WRITE (LUPRI,'(//)')
      CALL FLSHFO(LUPRI)
      RETURN
 1000 FORMAT (A1,A6,1X,A10,F17.10,2F21.10)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck heade2 */
      SUBROUTINE HEADE2(HEAD,CFIRST,IN)
C*****************************************************************************
C
C     Print header like SUBR HEADER
C     except for the first character being CFIRST
C
      CHARACTER HEAD*(*), CFIRST*1
#include "priunit.h"
C
      LENGTH = LEN(HEAD)
      IF (IN .GE. 0) THEN
         INDENT = IN + 1
      ELSE
         INDENT = (72 - LENGTH)/2 + 1
      END IF
      WRITE (LUPRI, '(//,80A)') CFIRST,(' ',I=1,INDENT), HEAD
      WRITE (LUPRI, '(80A)') CFIRST,(' ',I=1,INDENT), ('-',I=1,LENGTH)
      WRITE (LUPRI, '()')
      RETURN
      END
      SUBROUTINE RPRIGEO(COORD,CFIRST)
C
C     Jun 29 1988 tuh
C     Slight modified Feb 17 1998 jth
C     Copied from abacus/abatro.F (TROINI)
C     and modified for use in Dirac./jth 1998-02-04
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      CHARACTER CFIRST*1
      DIMENSION COORD(3,*)
C
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
C
C
      IATOM = 0
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 200 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               WRITE (LUPRI,'(A1,A,F17.10,2F24.10)') 
     *           CFIRST,NAMDEP(IATOM),
     *           (PT(IBTAND(ISYMAX(I,1),ISYMOP))*COORD(I,ICENT),I=1,3)
            END IF
  200    CONTINUE
  100 CONTINUE
      WRITE (LUPRI,'(//)')
      RETURN
      END
C  /* Deck troini */
      SUBROUTINE RTROINI(LCOMPA)
C
C     Initialize /CBITRO/
C
C     Copied from abacus/abatro.F (TROINI)
C     and modified for use in Dirac./jth 1998-02-04
C
#include "implicit.h"
#include "mxcent.h"
#include "cbitro.h"
#include "abainf.h"
#include "dcbgen.h"
      LOGICAL LCOMPA
C
      IPRINT = IPRGEN
      THRESH = 0.1D00
      COMPAR = LCOMPA
c     SKIP   = .NOT. ((MOLGRD .OR. MOLHES .OR. DIPDER) .AND. DOSYM(1))
      SKIP   = .FALSE.
      CUT    = .FALSE.
      NOROT1 = .FALSE.
      HESTRO = .TRUE.
      GDALL  = .FALSE.
      GDTRO  = .FALSE.
      RDTRO  = .FALSE.
c     TROGRD = MOLGRD .AND. .NOT. (NFIELD .GT. 0)
      TROGRD = .TRUE.
c     TROHES = MOLHES
      TROHES = .FALSE.
c     TRODIP = DIPDER
      TRODIP = .FALSE.
      RETURN
      END
C
C
      SUBROUTINE GETDSO(DSO,WORK,LWORK,IPRINT)
C
C     Get total density matrix in SO basis (used for SCDENS & LCDENS)
C
C     31-July-2001 HJAaJ
C
#include "implicit.h"
      DIMENSION DSO(*), WORK(LWORK)
C
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcbdhf.h"
#include "dcborb.h"
C
#include "memint.h"
C
      CALL QENTER('GETDSO')
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','GETDSO')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUMMY,DUMMY,DUMMY,3)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Get total density matrix
C
      CALL GENDEN(DSO,WORK(KCMO),1,IPRINT)
      IF (NASHT.GT.0) THEN
         IF (NOPEN > 1) THEN
            CALL QUIT('GETDSO only implemented for one open shell')
         END IF
         CALL MEMGET('REAL',KDVSO,N2BBASXQ,WORK,KFREE,LFREE)
         CALL GENDEN(WORK(KDVSO),WORK(KCMO),2,IPRINT)
         CALL DAXPY(N2BBASXQ,DF(1),WORK(KDVSO),1,DSO,1)
      END IF
      CALL MEMREL('GETDSO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('GETDSO')
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      RETURN
      END
C
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GRDNRM(GRAD,CSTRA,SCTRA,IPRINT)
C*****************************************************************************
C
C     Written by Jesper Pedersen 4/3-2002
C
C     Transform molecular gradient and print calculate the norm -
C     for printing norm of gradient in single-point calculations.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "abainf.h"
#include "nuclei.h"
#include "symmet.h"
      DIMENSION GRAD(MXCOOR),CGRAD(MXCOOR),CSTRA(*),SCTRA(*)
C
C     Transform gradient to Cartesian coordinates
C
      CALL TRAGRD(GRAD,CGRAD,CSTRA,SCTRA,NCRREP(0,1),3*NUCDEP)
      IF(IPRINT .GE. 5) THEN 
         JOFF=0
         DO 100 I = 1,NUCDEP
            WRITE(LUPRI,'(3X,A10,3F21.10)') 
     &            NAMDEP(I),(CGRAD(JOFF+J),J=1,3)
            JOFF=JOFF+3
100      CONTINUE
         WRITE(LUPRI,'(A,//)') ' '
      ENDIF
C
C     Calculate norm
C
      TMP = DDOT(3*NUCDEP,CGRAD,1,CGRAD,1)
      GRDNORM=SQRT(TMP)
      CALL HEADER('Norm of total molecular gradient',-1)
      WRITE(LUPRI,'(17X,A22,F17.10,///)') 
     &      'Norm of gradient  :  ',GRDNORM
C
C
      RETURN
      END
C*****************************************************************************
