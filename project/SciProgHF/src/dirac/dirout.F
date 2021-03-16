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
C  /* Deck dhfout */
      SUBROUTINE DHFOUT(CMO,EIG,IBEIG,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Output module for DIRAC
C
C***********************************************************************
      use dirac_cfg
      use fde_evaluators_dirac, only: fde_calculate_interaction_energy

#ifdef HAS_PCMSOLVER
      use pcm_scf, only: get_pcm_energy, pcm_scf_finalize
#endif      

#ifdef HAS_PELIB
      use pe_variables, only: peqm
#endif

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "nuclei.h"
#include "frame.h"
#include "dcblab.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbpsi.h"
      LOGICAL LBIT,LVCON
      CHARACTER CPUTID*12,FMT*6,MXFORM*6,SECTID*12,DAYTID*10
      DIMENSION CMO(*),EIG(*),IBEIG(*),WORK(LWORK)
      INTEGER ITER
#include "memint.h"
C
C     Summary of calculation
C     ======================
C
      LVCON = DOLVC.AND.(INTFLG.EQ.3 .or. INTFLG.EQ.11)
      IF(.NOT.ONESYS) THEN
        CALL HEADER('SCF - CYCLE',-1)
        IF(ERGCNV) THEN
          WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &      '* Convergence on total energy.',
     &      '  Desired convergence :',SCFCNV(1),
     &      '  Allowed convergence:',SCFCNV(2)
        ELSEIF(FCKCNV) THEN
          WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &      '* Convergence on total Fock matrix.',
     &      '  Desired convergence:',SCFCNV(1),
     &      '  Allowed convergence:',SCFCNV(2)
        ELSEIF(EVCCNV) THEN
          WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &  '* Convergence on norm of error vector (gradient).',
     &  '  Desired convergence:',SCFCNV(1),
     &  '  Allowed convergence:',SCFCNV(2)
        ENDIF
        WRITE(LUPRI,'(A)')
     &   '* ERGVAL - convergence in total energy',
     &   '* FCKVAL - convergence in maximum change'//
     &     ' in total Fock matrix',
     &   '* EVCVAL - convergence in error vector (gradient)'
        CALL PRSYMB(LUPRI,'-',128,0)
        WRITE(LUPRI,1000)
        CALL PRSYMB(LUPRI,'-',128,0)
        OPEN(LUCYCL,FILE ='DFCYCL',STATUS='OLD',ACCESS='SEQUENTIAL',
     &              FORM = 'FORMATTED')
        REWIND LUCYCL
        DO 10 I = 1,NITER
          READ(LUCYCL,1020) ITER,DHFTMP,ERGVAL,FCKVAL,EVCVAL,
     &                 CACC,CPUTID,DHF_INTTYP,DAYTID
          WRITE(LUPRI,1010) ITER,DHFTMP,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,DHF_INTTYP,DAYTID
   10   CONTINUE
        CLOSE(LUCYCL,STATUS='KEEP')
        CALL PRSYMB(LUPRI,'-',128,0)
        IF(DHFCONV(1)) THEN
          WRITE(LUPRI,'(A,I5,A)')
     +       '* Convergence after',NITER,' iterations.'
        ELSEIF (DHFCONV(2)) THEN
          WRITE(LUPRI,'(A,I5,A)')
     +       '* Desired convergence limit not reached after',NITER,
     +       ' iterations but the current convergence is acceptable.'
        ELSE
          WRITE(LUPRI,'(A,I5,A)')
     +       '* No convergence after',NITER,' iterations'
        ENDIF
        WRITE(LUPRI,'(A)') '* Average elapsed time per iteration: '
        DHF_INTTYP = 'No 2-ints   '
        IF(ITRSCF(0).GT.0) THEN
          SCFTID(0) = SCFTID(0)/ITRSCF(0)
          CPUTID    = SECTID(SCFTID(0))
          WRITE(LUPRI,'(6X,A,A,3X,A12)') DHF_INTTYP,' :',CPUTID
        ENDIF
        DO I = 1,7
        IF(ITRSCF(I).GT.0) THEN
          SCFTID(I) = SCFTID(I)/ITRSCF(I)
          CPUTID    = SECTID(SCFTID(I))
          DHF_INTTYP= '            '
          IF(LBIT(I,1)) DHF_INTTYP(1:2) = 'LL'
          IF(LBIT(I,2)) DHF_INTTYP(4:5) = 'SL'
          IF(LBIT(I,3)) DHF_INTTYP(7:8) = 'SS'
          IF(LBIT(I,4)) DHF_INTTYP(10:11) = 'GT'
          WRITE(LUPRI,'(6X,A,A,3X,A12)') DHF_INTTYP,' :',CPUTID
        ENDIF
        ENDDO
        CALL HEADER('TOTAL ENERGY',-1)
        ERGMAX = MAX(ABS(ELERGY),ABS(DHFERG))*10
        FMT = MXFORM(ERGMAX,20)
C       Start by writing the computed electronic energy
        if (dirac_cfg_fde) then
           WRITE(LUPRI,'(A)') ''
           WRITE(LUPRI,'(2X,A)') 'FDE active subsystem contributions'
           WRITE(LUPRI,'(A)') ''
        endif
        WRITE(LUPRI,'(3X,A,'//FMT//')')
     +      'Electronic energy                        :   ',ELERGY
C       The total energy consists of the electronic energy plus 
C       additional terms
        WRITE(LUPRI,'(/3X,A)') 'Other contributions to the total energy'
C       We will always want to add the repulsion between the nuclei
        WRITE(LUPRI,'(3X,A,'//FMT//')')
     +      'Nuclear repulsion energy                 :   ',POTNUC
C       In case we invoked a solvent model we need to add the 
C       solvation energy as well
C       (LV: I have not changed the name into solvation energy but 
C            think this would be more appropriate than solvent energy, 
C            leave this to authors of this part) Done, hjaaj.
        IF(SOLVEN) THEN
          ERGSOL = ESOLVN + ESOLVE
          WRITE(LUPRI,'(3X,A,'//FMT//')')
     +      'Solvation energy                         :   ',ERGSOL
        ENDIF
#ifdef HAS_PCMSOLVER        
        if (dirac_cfg_pcm) then
          write(lupri,'(3X,A,'//FMT//')')
     +  'Solvation energy                         :   ',get_pcm_energy()
        end if 
#endif
#ifdef HAS_PELIB
! edh should be replaced with energy routine
        IF (PEQM) THEN
            write(lupri,'(3X,A,'//FMT//')')
     +  'Embedding energy                         :   ', E_PE
        END IF
#endif
C       A one-center approximation may be used that gives a correction 
C       to the electronic energy
        IF (ONECAP .AND. INTV1C.EQ.2) THEN
          WRITE(LUPRI,'(3X,A,'//FMT//')')
     +      'ONECAP model 2 SS correction             :   ',CORRLV
C       The SCC model may be used in cases where no SS integrals were 
C       computed
        ELSE IF(LVCON) THEN
          WRITE(LUPRI,'(3X,A,'//FMT//')')
     +      'SS Coulombic correction                  :   ',CORRLV
        ENDIF
C       We now write the sum, the actual summation is done in the 
C       ERGCAL routine 
        WRITE(LUPRI,'(/3X,A)') 'Sum of all contributions to the energy'
        WRITE(LUPRI,'(3X,A,'//FMT//')')
     +      'Total energy                             :   ',DHFERG

        if (dirac_cfg_fde) then
          call fde_calculate_interaction_energy(DHFERG,'DHF ')
        endif
      ENDIF
C
C     Check that there are no eigenvalues in the gap
C
      DO IFRP = 1,NFSYM
         CALL EIGCHK(EIG(IORB(IFRP)+1),NPSH(IFRP))
      END DO
C
C     Print eigenvalues
C
      IF(LBIT(IPREIG,1))
     &   CALL PREIGN(EIG,IBEIG,LBIT(IPREIG,2),WORK,LWORK)
C
C     Check HOMO - LUMO gap, and if there are any negative
C     energy virtuals.
C
      CALL CHKHLG(EIG,IBEIG)
C
#ifdef HAS_PCMSOLVER      
      if (dirac_cfg_pcm) then
              call pcm_scf_finalize
       write(lupri, '(//A/)') 'PCMSolver interface correctly finalized'
      end if
#endif      
      RETURN
C
 1000 FORMAT(11X,'Energy',15X,'ERGVAL',4X,'FCKVAL',4X,'EVCVAL',6X,
     &       'Conv.acc',4X,'CPU',10X,'Integrals',3X,'Time stamp')
 1010 FORMAT('It. ',I4,1P,G23.13,3D10.2,3X,A8,3X,A12,3X,A12,3X,A10)
 1020 FORMAT(4X,    I4,   G23.13,3D10.2,3X,A8,3X,A12,3X,A12,3X,A10)
C     ... 1020 is used for reading, thus no '1P' factor !
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prcoef */
      SUBROUTINE PRCOEF(CMO,EIG,VECPRI,IPREIG,IPRCMP,ILABDF,IPRINT,
     &                  WORK,LWORK)
C***********************************************************************
C
C     Print coefficients in
C       AO-basis  (ILABDF = 1)
C       SO-basis  (ILABDF = 2)
C     VECPRI is an character array that contains information about which
C     vectors to print.
C
C     Written by T.Saue 1996
C     Last revision Jan 9 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      CHARACTER VECPRI*72
      DIMENSION CMO(*),EIG(*),VECPRI(NFSYM),WORK(LWORK)
C
      CALL QENTER('PRCOEF')
#include "memint.h"
C
C     Loop over fermion ircops
C     ========================
C
      DO 10 IFRP = 1,NFSYM
C
C       Find number of electronic/positronic vectors
C       ============================================
C
        NVEC = 1
        CALL MEMGET('INTE',KJVEC,NORB(IFRP),WORK,KFREE,LFREE)
        CALL NUMLST(VECPRI(IFRP),WORK(KJVEC),NORB(IFRP),
     &              -NPSH(IFRP),NESH(IFRP),IFRP,NVEC)
        CALL ORBCNT(WORK(KJVEC),NVEC,NPSH(IFRP),NESH(IFRP),NPVEC,NEVEC)
        IF(NVEC.EQ.0) GOTO 20
        CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
        NCBAS=NTBAS(0)
        IF(ILABDF.EQ.2) NCBAS=NFBAS(IFRP,0)
        CALL MEMGET('REAL',KCBF,NCBAS*NVEC*4,WORK,KFREE,LFREE)
        CALL SELCMO(ILABDF,IFRP,CMO,WORK(KJVEC),NVEC,NPVEC,NEVEC,
     &              WORK(KCBF),NCBAS,NVEC,NCBAS,
     &              IPRINT,WORK,KFREE,LFREE)
C
C       Print coefficients
C       ==================
C
        CALL PRCOE1(WORK(KCBF),EIG(IORB(IFRP)+1),IFRP,WORK(KJVEC),
     &              NPVEC,NEVEC,NCBAS,IPREIG,IPRCMP,ILABDF)
 20     CONTINUE
        CALL MEMREL('PRCOEF',WORK,1,KWORK,KFREE,LFREE)
 10   CONTINUE
C
C       Memory deallocation
C
      CALL MEMREL('PRCOEF.PRCOE1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('PRCOEF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prcoe1 */
      SUBROUTINE PRCOE1(CMO,EIG,IFRP,JVEC,NPVEC,NEVEC,NCBAS,
     &                  IPREIG,IPRCMP,ILABDF)
C***********************************************************************
C
C     Print vectors (AO/SO-basis)
C     The coefficients, whether in AO- or SO-basis, are provided on 
C     quaternion (NZ=4) form, which corresponds to
C  
C       C^Q = C^a - C^b* = C^a_R + C^a_I - C^b_R + C^b_I
C
C     It can be seen that the part C^b_R should be multiplied by -1.
C     
C     Written by T.Saue May 6 1996
C
C***********************************************************************
#include "implicit.h"
      PARAMETER(D0 = 0.0D0,DM10 = 1.0D-10,DM1 = -1.0D0)
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcblab.h"
#include "dgroup.h"
#include "symmet.h"
#include "dcbbas.h"
#include "dcborb.h"
      LOGICAL LBIT
      CHARACTER FMT*6,MXFORM*6
      DIMENSION CMO(NCBAS,NPVEC+NEVEC,4),EIG(*),JVEC(*)
C
      NVEC = NPVEC + NEVEC
      NTOT = NVEC*NCBAS
C
C.....Properly scale C^b_R
      CALL DSCAL(NTOT,DM1,CMO(1,1,3),1)
C
C     Print vectors in AO-basis
C     =========================
C
      IF(ILABDF.EQ.1) THEN
        DO I = 1,NVEC
          IF(I.GT.NPVEC) THEN
            IF (IPREIG .EQ. 1) THEN
              II = NPSH(IFRP)+JVEC(I)
              FMT = MXFORM(EIG(II),16)
              WRITE(LUPRI,'(/A,I3,A,'//FMT//')')
     +         '* Electronic eigenvalue no.',JVEC(I),': ', EIG(II)
            ELSE
              WRITE(LUPRI,'(/A,I3)')
     +         '* Electronic orbital no.',JVEC(I)
            END IF
          ELSE
            IF (IPREIG .EQ. 1) THEN
              II = NPSH(IFRP)+1+JVEC(I)
              FMT = MXFORM(EIG(II),16)
              WRITE(LUPRI,'(/A,I3,A,'//FMT//')')
     +         '* Positronic eigenvalue no.',II,': ', EIG(II)
            ELSE
              WRITE(LUPRI,'(/A,I3)')
     +         '* Positronic orbital no.',JVEC(I)
            END IF
          ENDIF
          CALL PRSYMB(LUPRI,'=',52,0)
          IOFF = 0
          DO IC = 1,2
          IF(LBIT(IPRCMP,IC)) THEN
            DO J = 1,NTBAS(IC)
              JJ  = IOFF + J
              DVEC = ABS(CMO(JJ,I,1)) + ABS(CMO(JJ,I,2))
     &             + ABS(CMO(JJ,I,3)) + ABS(CMO(JJ,I,4))
              IF(DVEC.GT.DM10) WRITE(LUPRI,100)
     &          JJ,PLABEL(IPLAB(JJ,1),1),(CMO(JJ,I,IZ),IZ = 1,4)
            ENDDO
          ENDIF
          IOFF = IOFF + NTBAS(IC)
          ENDDO
        ENDDO
C
C     Print vectors in SO-basis
C     =========================
C
      ELSEIF(ILABDF.EQ.2) THEN
        NVEC = NPVEC+NEVEC
        NBRP = 4/NZ
        DO I = 1,NVEC
          IF(I.GT.NPVEC) THEN
            IF (IPREIG .EQ. 1) THEN
              II = NPSH(IFRP)+JVEC(I)
              FMT = MXFORM(EIG(II),16)
              WRITE(LUPRI,'(/A,I3,A,'//FMT//')')
     +         '* Electronic eigenvalue no.',JVEC(I),': ', EIG(II)
            ELSE
              WRITE(LUPRI,'(/A,I3)')
     +         '* Electronic orbital no.',JVEC(I)
            END IF
          ELSE
            IF (IPREIG .EQ. 1) THEN
              II = NPSH(IFRP)+1+JVEC(I)
              FMT = MXFORM(EIG(II),16)
              WRITE(LUPRI,'(/A,I3,A,'//FMT//')')
     +         '* Positronic eigenvalue no.',II,': ',EIG(II)
            ELSE
              WRITE(LUPRI,'(/A,I3)')
     +         '* Positronic orbital no.',JVEC(I)
            ENDIF
          ENDIF
          CALL PRSYMB(LUPRI,'=',52,0)
          DO IC = 1,2
          IF(LBIT(IPRCMP,IC)) THEN
C           IP is gerade(1)/ungerade(2) !
            IP    = MOD(IFRP+IC,2) + 1
            DO ISYM = 1,NBRP
              IBSYM = JFSYM(ISYM,IP)
              IBREP = IBSYM - 1
              DO J = 1,NBBAS(IBREP,IC)
                JU    = ICOS(IBSYM,IC)  + J
                JS    = IBBAS(IBREP,IC) - IBAS(IFRP) + J
                DVEC = ABS(CMO(JS,I,1)) + ABS(CMO(JS,I,2))
     &               + ABS(CMO(JS,I,3)) + ABS(CMO(JS,I,4))
!                IF(DVEC.GT.DM10) WRITE(LUPRI,100)
!     &            JU,PLABEL(IPLAB(IBAS(IFRP)+JS,2),2),
!     &               (CMO(JS,I,IZ),IZ = 1,4)
                WRITE(LUPRI,100)
     &            JS,PLABEL(IPLAB(IBAS(IFRP)+JS,2),2),
     &               (CMO(JS,I,IZ),IZ = 1,4)
              ENDDO
            ENDDO
          ENDIF
          ENDDO
        ENDDO
      ENDIF
C.....Rescale C^b_R
      CALL DSCAL(NTOT,DM1,CMO(1,1,3),1)
C
      RETURN
 100  FORMAT(3X,I5,2X,A12,2X,4F14.10)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck cfsort */
      SUBROUTINE QEXPAND(CBUF,CMO,NVEC,IFRP,KNFBAS)
C***********************************************************************
C
C     In the case of NZ.LT.4 this routine expands real/imaginary matrix
C     to quaternionic
C     KNFBAS = NFBAS(IFRP,0)
C
C     Written by T.Saue November 1994
C     Last revision: Jan 7 1998 - jth NFBAS(IFRP,0) is now parameter
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcborb.h"
      DIMENSION CMO(KNFBAS,NVEC,NZ),CBUF(KNFBAS,NVEC,4)
C
      IF(NVEC.EQ.0) RETURN
      NDIM  = NFBAS(IFRP,0)
      NTOT  = NVEC*NDIM
      NTOTZ = NTOT*NZ
      IF(NZ.NE.4) CALL DZERO(CBUF,NTOT*4)
      NBRP = 4/NZ
      DO IC = 1,MC
        IP     = MOD(IFRP+IC,2) + 1
        DO ISYM = 1,NBRP
          IBRP = JFSYM(ISYM,IP) - 1
          IOFF = IBBAS(IBRP,IC) - IBAS(IFRP)
          DO IZ = 1,NZ
            IQ = IQFROMPQ(ISYM,IFRP,IC,IZ)
            DO I = 1,NBBAS(IBRP,IC)
              CALL DCOPY(NVEC,CMO (I+IOFF,1,IZ),NFBAS(IFRP,0),
     &                        CBUF(I+IOFF,1,IQ),NFBAS(IFRP,0))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck preigv */
      SUBROUTINE PREIGV(EIG,NVAL,IVAL,WORK,LWORK)
C***********************************************************************
C
C     Prints eigenvalues and their degeneracies
C
C     Written by T. Saue
C     LAST VERSION: July 23 1994
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION EIG(*),WORK(LWORK)
C
      CALL QENTER('PREIGV')
#include "memint.h"
C
      IF(IVAL.LT.0) IVAL = 0
C
C     Memory allocation
      CALL MEMGET('INTE',KFIRST,NVAL,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNMULT ,NVAL,WORK,KFREE,LFREE)
      CALL PREIG1(EIG(IVAL+1),NVAL,WORK(KFIRST),WORK(KNMULT))
C
C     Memory deallocation
      CALL MEMREL('PREIGV.PREIG1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('PREIGV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck preig1 */
      SUBROUTINE PREIG1(EIG,NVAL,IFIRST,NMULT)
C***********************************************************************
C
C     Prints eigenvalues and their degeneracies
C
C     Written by T. Saue
C     LAST VERSION: July 23 1994
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0 = 0.0D0)
      CHARACTER FMT*6,MXFORM*6
      DIMENSION EIG(NVAL),IFIRST(NVAL),NMULT(NVAL)
C
      CALL QENTER('PREIG1')
C
C     Find max value
C     ==============
C
      EIGMAX = D0
      DO J = 1,NVAL
        EIGMAX = MAX(EIGMAX,ABS(EIG(J)))
      ENDDO
      FMT = MXFORM(EIGMAX,18)
C
C     Find degeneracies
C     =================
C
      CALL DEGANA(NDEG,NVAL,EIG,IFIRST,NMULT)
C
C     Print section
C     =============
C
      WRITE(LUPRI,'(5('//FMT//',2X,A1,I2,A1,3X))')
     &     (EIG(IFIRST(I)),'(',2*NMULT(I),')',I=1,NDEG)
C
      CALL QEXIT('PREIG1')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck degana */
      SUBROUTINE DEGANA(NDEG,NVAL,EIG,IFIRST,NMULT)
C***********************************************************************
C
C     Given a list EIG of ordered eigenvalues, this routine
C     finds degeneracies. For each disctinct value, it will list
C     its degeneracy(NMULT) and first occurrence(IFIRST)
C
C     Written by T. Saue
C     LAST VERSION: July 23 1994
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DEGTOL = 1.0D-6)
      DIMENSION EIG(NVAL),IFIRST(NVAL),NMULT(NVAL)
C
      CALL QENTER('DEGANA')
C
      NDEG = 1
      DEGVAL    = EIG(1)
      IFIRST(1) = 1
      NMULT(1)  = 1
      DO 10 I = 2,NVAL
        DEGTST = ABS(EIG(I)-DEGVAL)
C       We better take the relative difference for large eigenvalues
        AEIG = ABS(EIG(I))
        IF (AEIG.GT.1.0D3) DEGTST = DEGTST / AEIG
C
C       Degeneracy
C       ==========
C
        IF(DEGTST.LT.DEGTOL) THEN
          NMULT(NDEG) = NMULT(NDEG)+1
        ELSE
C
C       First distinct value
C       ====================
C
          NDEG         = NDEG + 1
          DEGVAL       = EIG(I)
          IFIRST(NDEG) = I
          NMULT(NDEG)  = 1
        ENDIF
   10 CONTINUE
C
      CALL QEXIT('DEGANA')
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck eigchk */
      SUBROUTINE EIGCHK(EIG,NVAL)
C***********************************************************************
C
C     This routine checks for intruder states
C
C     Written by T.Saue - November 1994
C     Last revision : Nov 15 1994 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
      DIMENSION EIG(NVAL)
C
      NGAP = 0
      GAP = -2*CVAL*CVAL
      DO 10 I = NVAL,1,-1
        IF(EIG(I).LE.GAP) GOTO 20
        NGAP = NGAP + 1
   10 CONTINUE
   20 CONTINUE
      IF(NGAP.GT.0) THEN
        CALL HEADER('WARNING from EIGCHK',-1)
        WRITE(LUPRI,'(3X,I5,A)') NGAP,' positron states intruding'
        WRITE(LUPRI,'(6X,I5,2F18.8)')
     *     ((NVAL-I),EIG(NVAL-I),EIG(NVAL-I)-GAP,I=0,(NGAP-1))
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ORBCNT(JVEC,NVEC,NNEG,NPOS,NPVEC,NEVEC)
C***********************************************************************
C     "Orbital Count"
C     Find NEVEC/NPVEC, the number of electronic/positronic orbital vectors,
C     in the list of orbitals in JVEC(1:NVEC).
C     Orbitals outside the range -NNEG:NPOS are ignored and eliminated
C     from JVEC.
C
C     In   : NNEG, NPOS
C     InOut: NVEC; JVEC(1:NVEC)
C     Out  : NEVEC,NPVEC
C
C     Written by T.Saue Jan 7 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION JVEC(NVEC)
C
      NPVEC = 0
      NEVEC = 0
      II = 0
      DO 10 I = 1,NVEC
        J = JVEC(I)
C
C       Positronic orbital in range [-NNEG(IFRP),-1]
C
        IF    (J.LT.0.AND.J.GE.(-NNEG)) THEN
          II       = II + 1
          JVEC(II) = J
          NPVEC    = NPVEC + 1
C
C       Electronic orbital in range [1,NPOS]
C
        ELSEIF (J.GT.0.AND.J.LE.NPOS) THEN
          II       = II + 1
          JVEC(II) = J
          NEVEC    = NEVEC + 1
        ENDIF
 10   CONTINUE
C
      NVEC = NPVEC+NEVEC
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SELCFS(CMO,IFRP,CBUF,LCB,JVEC,NPVEC,NEVEC,
     &                 KNFBAS,KNORB)
C***********************************************************************
C     Pick out a set of vectors from CMO according to array JVEC
C
C     Written by T.Saue 1997
C     Last revision Jan 8 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
C     KNFBAS = NFBAS(IFRP,0)
C     KNORB  = NORB(IFRP)
C
      DIMENSION CMO(KNFBAS,KNORB,NZ),CBUF(KNFBAS,LCB,NZ),
     &          JVEC(NPVEC+NEVEC)
C
      NVEC = NPVEC + NEVEC
      DO IZ = 1,NZ
C
C       Positronic vectors
C
        DO I = 1,NPVEC
          II = NPSH(IFRP)+1+JVEC(I)
          CALL DCOPY(NFBAS(IFRP,0),CMO(1,II,IZ),1,CBUF(1,I,IZ),1)
        ENDDO
C
C       Electronic vectors
C
        DO I = NPVEC+1,NVEC
          II = NPSH(IFRP)+JVEC(I)
          CALL DCOPY(NFBAS(IFRP,0),CMO(1,II,IZ),1,CBUF(1,I,IZ),1)
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PHATRA */
      SUBROUTINE PHATRA(VEC,NDIM,NVEC,NZ)
C***********************************************************************
C
C     Phase transform of vectors
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0, DM1 = -1.0D0)
      DIMENSION VEC(NDIM,NVEC,NZ),A(4),B(4)
C
      IF    (NZ.EQ.1) THEN
C
C     Real vectors
C     ============
C
        DO J = 1,NVEC
C
C         Find largest absolute value
C
          IB = 1
          DB = ABS(VEC(1,J,1))
          DO I = 2,NDIM
            DM = ABS(VEC(I,J,1))
            IF(DM.GT.DB) THEN
              DB = DM
              IB = I
            ENDIF
          ENDDO
C
C         Scale with minus one if largest element is negative
C
          IF(DM.LT.D0) CALL DSCAL(NVEC,DM1,VEC,1)
        ENDDO
      ELSEIF(NZ.EQ.2) THEN
C
C     Complex vectors
C     ===============
C
        DO J = 1,NVEC
C
C         Find largest absolute value
C
          IB = 1
          DB = PYTHAG(VEC(1,J,1),VEC(1,J,2))
          DO I = 2,NDIM
            DM = PYTHAG(VEC(I,J,1),VEC(I,J,2))
            IF(DM.GT.DB) THEN
              DB = DM
              IB = I
            ENDIF
          ENDDO
C
          DO IZ = 1,NZ
            A(IZ) =  VEC(IB,J,IZ)/DB
          ENDDO
          DO I = 1,NDIM
            B(1)       = VEC(I,J,2)*A(1) - VEC(I,J,1)*A(2)
            VEC(I,J,1) = VEC(I,J,1)*A(1) + VEC(I,J,2)*A(2)
            VEC(I,J,2) = B(1)
          ENDDO
        ENDDO
      ELSEIF(NZ.EQ.4) THEN
C
C     Quaternion vectors
C     ==================
C
        DO J = 1,NVEC
C
C         Find largest absolute value
C
          IB = 1
          DB = PYTHAG(PYTHAG(VEC(1,J,1),VEC(1,J,2)),
     &                PYTHAG(VEC(1,J,3),VEC(1,J,4)))
          DO I = 2,NDIM
            DM = PYTHAG(PYTHAG(VEC(I,J,1),VEC(I,J,2)),
     &                  PYTHAG(VEC(I,J,3),VEC(I,J,4)))
            IF(DM.GT.DB) THEN
              DB = DM
              IB = I
            ENDIF
          ENDDO
C
          DO IZ = 1,NZ
            A(IZ) =  VEC(IB,J,IZ)/DB
          ENDDO
          DO I = 1,NDIM
            B(1)       = VEC(I,J,2)*A(1) - VEC(I,J,1)*A(2)
     &                 + VEC(I,J,4)*A(3) - VEC(I,J,3)*A(4)
            B(2)       = VEC(I,J,3)*A(1) - VEC(I,J,4)*A(2)
     &                 - VEC(I,J,1)*A(3) + VEC(I,J,2)*A(4)
            B(3)       = VEC(I,J,4)*A(1) + VEC(I,J,3)*A(2)
     &                 - VEC(I,J,2)*A(3) - VEC(I,J,1)*A(4)
            VEC(I,J,1) = VEC(I,J,1)*A(1) + VEC(I,J,2)*A(2)
     &                 + VEC(I,J,3)*A(3) + VEC(I,J,4)*A(4)
            VEC(I,J,2) = B(1)
            VEC(I,J,3) = B(2)
            VEC(I,J,4) = B(3)
          ENDDO
        ENDDO
      ELSE
        WRITE(LUPRI,'(A,I5)') 'PHATRA: * Unknown NZ = ',NZ
        CALL QUIT('PHATRA: Unknown NZ !')
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SELEIG(EIG,IFRP,EBUF,JVEC,NPVEC,NEVEC)
C***********************************************************************
C     Pick out a set of eigenvalues from EIG according to array JVEC
C
C     Written by L. Visscher 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION EIG(*),EBUF(NPVEC+NEVEC),JVEC(*)
C
      NVEC = NPVEC + NEVEC
C
C     Positronic eigenvalues
C
      DO I = 1,NPVEC
        II = NPSH(IFRP)+1+JVEC(I)
        EBUF(I) = EIG(II)
      ENDDO
C
C     Electronic eigenvalues
C
      DO I = NPVEC+1,NVEC
        II = NPSH(IFRP)+JVEC(I)
        EBUF(I) = EIG(II)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SELIBEIG(IBEIG,IOFF,IFRP,IBBUF,JVEC,NPVEC,NEVEC)
C***********************************************************************
C     Pick out a set of integers from IBEIG according to array JVEC
C
C     Written by L. Visscher 1999
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION IBEIG(*),IBBUF(*),JVEC(*)
C
      NVEC = NPVEC + NEVEC
C
C     Positronic eigenvalues
C
      DO I = 1,NPVEC
        II = NPSH(IFRP)+1+JVEC(I)
        IBBUF(I) = IBEIG(IOFF+II)
      ENDDO
C
C     Electronic eigenvalues
C
      DO I = NPVEC+1,NVEC
        II = NPSH(IFRP)+JVEC(I)
        IBBUF(I) = IBEIG(IOFF+II)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck preigb */
      SUBROUTINE PREIGB(EIG,IBEIG,IVAL,NVAL)
C***********************************************************************
C
C     Prints eigenvalues and their degeneracies
C
C     Written by T. Saue
C     LAST VERSION: July 23 1994
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
#include "pgroup.h"
C
      DIMENSION EIG(*),IBEIG(*)
      CHARACTER FMT*6,MXFORM*6
C
      EIGMAX = D0
      DO J = IVAL+1,IVAL+NVAL
        EIGMAX = MAX(EIGMAX,ABS(EIG(J)))
      ENDDO
      FMT = MXFORM(EIGMAX,18)
C
C     Print section
C     =============
C
      WRITE(LUPRI,'(5('//FMT//',1X,A1,A3,A1,2X))')
     &(EIG(J),'(',REP(IBEIG(J)),')',J=IVAL+1,IVAL+NVAL)
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck preign */
      SUBROUTINE PREIGN(EIG,IBEIG,PRPOS,WORK,LWORK)
C***********************************************************************
C
C     Prints eigenvalues and their degeneracies
C
C     Written by J. Thyssen - 26 Aug 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
C
      DIMENSION EIG(NORBT),IBEIG(NORBT),WORK(LWORK)
      LOGICAL PRPOS
C
C
      CALL QENTER('PREIGN')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET2('INTE','IFIRST',KFIRST,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','NMULT' ,KNMULT,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EIGBUF',KEIGBF,NORBT,WORK,KFREE,LFREE)
      CALL PREIGN1(EIG,WORK(KEIGBF),IBEIG,
     &   WORK(KFIRST),WORK(KNMULT),PRPOS)
C
C     Memory deallocation
      CALL MEMREL('PREIGN',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('PREIGN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck preign */
      SUBROUTINE PREIGN1(EIG,EIGBUF,IBEIG,IFIRST,NMULT,PRPOS)
C***********************************************************************
C
C     Prints eigenvalues - new version.
C
C     Written by J. Thyssen - 26 Aug 1999
C     Modified by A. Sunaga - 20 Sep 2018 for Atomic output
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER(NCOL=16)
#include "dgroup.h"
#include "pgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbibt.h"            
C     
      DIMENSION EIG(NORBT),IBEIG(NORBT),EIGBUF(NORBT),
     &          IFIRST(NORBT),NMULT(NORBT)
      CHARACTER FMTE*6,FMTP*6,MXFORM*6, TLINE*(5*NCOL+5),
     &          TLINE2*(6*NCOL+5),TLINE3*(6*NCOL+5)
      CHARACTER*5 NAME_SUB2(-MAX_SUB_BL:MAX_SUB_BL),
     &            NAME_SUB3(-MAX_SUB_BL:MAX_SUB_BL)! for MJ of ATOMIC
      CHARACTER SPDCAR*1 ! type declaration for external function
      DIMENSION ISHEL_SUB(3+MXOPEN,MAX_SUB_BL,2)
      LOGICAL   PRPOS
C     
C
C     Closed shell eigenvalues
C
      CALL HEADER('Eigenvalues',-1)
C.....Common format for all electron eigenvalues
      EIGMAX = D0
C
      DO IFSYM = 1,NFSYM
        IF(NESH(IFSYM).GT.0)THEN 
          DO I = 1,NESH(IFSYM)
            EIGMAX = MAX(EIGMAX,ABS(EIG(IORB(IFSYM)+NPSH(IFSYM)+I)))
          ENDDO
        ENDIF
      END DO
      FMTE = MXFORM(EIGMAX,18)
C.....Common format for all positron eigenvalues
      EIGMAX = D0
      DO IFSYM = 1,NFSYM
        IF(NPSH(IFSYM).GT.0)THEN 
          DO I = 1,NPSH(IFSYM)
            EIGMAX = MAX(EIGMAX,ABS(EIG(IORB(IFSYM)+I)))
          ENDDO
        ENDIF
      END DO
      FMTP = MXFORM(EIGMAX,18)
C
C     Standard case
C     =============
C
      IF (.NOT. SUB_BL ) THEN
C
C       Print eigenvalues
C
        DO IFSYM = 1,NFSYM
          WRITE(LUPRI,'(/2A)') '* Fermion symmetry ',FREP(IFSYM)
C
C         Closed shell eigenvalues
C
          IOFF = IORB(IFSYM)+NPSH(IFSYM)
          IF (NISH(IFSYM) .GT. 0) THEN
            WRITE(LUPRI,'(2X,A,F6.4)') '* Closed shell, f = ',DF(0)
            CALL PREIGN2(EIG(IOFF+1),IFIRST,NMULT,NISH(IFSYM),FMTE)
          END IF
C
C         Open shell eigenvalues
C
          IOFF = IOFF+NISH(IFSYM)
          DO IOPEN = 1,NOPEN
            IF (NACSH(IFSYM,IOPEN) .GT. 0) THEN
              WRITE(LUPRI,'(2X,A,I1,A,F6.4)')
     &           '* Open shell #',IOPEN,', f = ',DF(IOPEN)
              CALL PREIGN2(EIG(IOFF+1),IFIRST,NMULT,
     &                     NACSH(IFSYM,IOPEN),FMTE)
              IOFF = IOFF + NACSH(IFSYM,IOPEN)
            END IF
          END DO
C
C         Virtual eigenvalues
C
          IF (NSSH(IFSYM) .GT. 0)  THEN
            WRITE(LUPRI,'(2X,A,F6.4)')
     &            '* Virtual eigenvalues, f = ',D0
            CALL PREIGN2(EIG(IOFF+1),IFIRST,NMULT,NSSH(IFSYM),FMTE)
          END IF
C
C         Positronic eigenvalues
C
          IF (PRPOS) THEN
            IF (NPSH(IFSYM) .GT. 0) THEN
              WRITE(LUPRI,'(2X,A,F6.4)')
     &               '* Positronic eigenvalues, f = ',D0
              CALL PREIGN2(EIG(IORB(IFSYM)+1),
     &                     IFIRST,NMULT,NPSH(IFSYM),FMTP)
            END IF
          END IF
        END DO
      ELSE
C
C     Supersymmetry: linear/atomic symmetry and/or spinfree case
C     ===================================================
C
        DO IFSYM = 1, NFSYM
          IF(N_SUB_BL(IFSYM).GT.0)THEN 
          DO ISUB = 1, N_SUB_BL(IFSYM)
C
C         Write identification of the blocks. 
C
            ID = ID_SUB_BL(ISUB,IFSYM)
            IF (SPINFR) THEN     ! Spinfree case
              WRITE(LUPRI,'(/2A)') '* Boson symmetry ',REP(ID)
              NAME_SUB2(ID) = '  '//REP(ID)
            ELSEIF(ATOMIC) THEN  ! Atomic supersymmetry
              CALL ATOMIC_ID(ID,KP,J,MJ,LL)
              WRITE(LUPRI,'(/A,I4,4A,I2,2A,I3,A)')
     &            '* Block',ISUB,' in ',FREP(IFSYM),
     &            ':  ',SPDCAR(LL),J,'/2',';',MJ,'/2'
              WRITE (NAME_SUB2(ID),'(A,I2,"/2")') SPDCAR(LL),J
              WRITE (NAME_SUB3(ID),'(I3,"/2")')MJ
            ELSE                 ! Linear supersymmetry       
              WRITE(LUPRI,'(/A,I4,3A,I2,A)')
     &            '* Block',ISUB,' in ',FREP(IFSYM),
     &            ':  Omega = ',ABS(ID),'/2'
              WRITE (NAME_SUB2(ID),'(I3,"/2")') ABS(ID)
            ENDIF
C
            NTOT = 0
C...........Closed shell eigenvalues: pick eigenvalues of symmetry ID
            IOFF = IORB(IFSYM)+NPSH(IFSYM)
            CALL PCKEIG(EIG,IOFF+1,IOFF+NISH(IFSYM),
     &                IBEIG,ID,NEIG,EIGBUF)
            ISHEL_SUB(1,ISUB,IFSYM) = NEIG
C
            IF (NEIG .GT. 0) THEN
              NTOT = NTOT + NEIG
              WRITE(LUPRI,'(2X,A,F6.4)') '* Closed shell, f = ',DF(0)
              CALL PREIGN2(EIGBUF,IFIRST,NMULT,NEIG,FMTE)
            END IF
C...........Open shell eigenvalues: pick eigenvalues of symmetry ID
            IOFF = IOFF+NISH(IFSYM)            
            DO IOPEN = 1,NOPEN
              CALL PCKEIG(EIG,IOFF+1,IOFF+NACSH(IFSYM,IOPEN),
     &                    IBEIG,ID,NEIG,EIGBUF)
              ISHEL_SUB(1+IOPEN,ISUB,IFSYM) = NEIG
C
              IF (NEIG .GT. 0) THEN
                NTOT = NTOT + NEIG
                WRITE(LUPRI,'(2X,A,I1,A,F6.4)')
     &              '* Open shell #',IOPEN,', f = ',DF(IOPEN)
                CALL PREIGN2(EIGBUF,IFIRST,NMULT,NEIG,FMTE)
              END IF
              IOFF = IOFF + NACSH(IFSYM,IOPEN)
            END DO
C...........Virtual eigenvalues: pick eigenvalues of symmetry ID
            CALL PCKEIG(EIG,IOFF+1,
     &                  IORB(IFSYM)+NORB(IFSYM),IBEIG,ID,NEIG,EIGBUF)
            ISHEL_SUB(1+MXOPEN+1,ISUB,IFSYM) = NEIG
C          
            IF (NEIG .GT. 0) THEN
              NTOT = NTOT + NEIG
              WRITE(LUPRI,'(2X,A,F6.4)')
     &            '* Virtual eigenvalues, f = ',D0
              CALL PREIGN2(EIGBUF,IFIRST,NMULT,NEIG,FMTE)
            END IF  
C...........Positronic eigenvalues: pick eigenvalues of symmetry ID
            ISHEL_SUB(1+MXOPEN+2,ISUB,IFSYM) = 0
            IF (PRPOS .AND. NPSH(IFSYM) .GT.0) THEN
              CALL PCKEIG(EIG,IORB(IFSYM)+1,IORB(IFSYM)+NPSH(IFSYM),
     &                    IBEIG,ID,NEIG,EIGBUF)
              ISHEL_SUB(1+MXOPEN+2,ISUB,IFSYM) = NEIG
C
              IF (NEIG .GT. 0) THEN
                NTOT = NTOT + NEIG
                WRITE(LUPRI,'(2X,A,F6.4)')
     &              '* Negative energy eigenvalues, f = ',DF(0)
                CALL PREIGN2(EIGBUF,IFIRST,NMULT,NEIG,FMTP)
              END IF
            END IF
C...........No eigenvalues
            IF(NTOT.EQ.0) THEN
              WRITE(LUPRI,'(2X,A)') '* No eigenvalues.'
            ENDIF
          END DO
          ENDIF
        END DO
C     
C       Print occupation string for each fermion symmetry
C
        IF (NBSYM.GT.1) THEN
          DO IFSYM = 1, NFSYM
            WRITE(LUPRI,'(/2A)')
     &          '* Occupation in fermion symmetry ',FREP(IFSYM)
C...........Closed shell
            JST = IORB(IFSYM) + NPSH(IFSYM)
            IF (NISH(IFSYM) .GT. 0) THEN
              WRITE(LUPRI,'(2X,A)') '* Inactive orbitals'
              TMP = NISH(IFSYM)
              DO I = JST,JST+NISH(IFSYM)-1,18
                IF(TMP.GT.18)THEN
                  JMAX = 18
                ELSE
                  JMAX = TMP
                ENDIF   
                IF(ATOMIC)THEN
                  WRITE(LUPRI,'(4X,18(A5,1X))')
     &                      (NAME_SUB2(IBEIG(J)),J=I+1,I+JMAX)
                  WRITE(LUPRI,'(2X,A,18(A5,1X))')
     &                 'Mj',(NAME_SUB3(IBEIG(J)),J=I+1,I+JMAX)
                ELSE
                  WRITE(LUPRI,'(4X,18(A5))')
     &                 (NAME_SUB2(IBEIG(J)),J=I+1,I+JMAX)
                ENDIF
                TMP = TMP - 18
              ENDDO
            END IF
C...........Open shells
            JST = JST + NISH(IFSYM)
            IF (NASH(IFSYM) .GT. 0) THEN
              WRITE(LUPRI,'(2X,A)') '* Active orbitals'
              TMP = NASH(IFSYM)
              DO I = JST,JST+NASH(IFSYM)-1,18
                IF(TMP.GT.18)THEN
                  JMAX = 18
                ELSE
                  JMAX = TMP
                ENDIF                   
                IF(ATOMIC)THEN
                  WRITE(LUPRI,'(4X,18(A5,1X))')
     &                      (NAME_SUB2(IBEIG(J)),J=I+1,I+JMAX)
                  WRITE(LUPRI,'(2X,A,18(A5,1X))')
     &                 'Mj',(NAME_SUB3(IBEIG(J)),J=I+1,I+JMAX)
                ELSE  
                  WRITE(LUPRI,'(4X,18(A5))')
     &                 (NAME_SUB2(IBEIG(J)),J=I+1,I+JMAX)
                ENDIF
                TMP = TMP - 18                    
              ENDDO
            END IF
C...........Virtual shells
            JST = JST + NASH(IFSYM)
            IF (NSSH(IFSYM) .GT. 0) THEN
              WRITE(LUPRI,'(2X,A)') '* Virtual orbitals'
              TMP = NSSH(IFSYM)
              DO I = JST,JST+NSSH(IFSYM)-1,18
                IF(TMP.GT.18)THEN
                  JMAX = 18
                ELSE
                  JMAX = TMP
                ENDIF                
                IF(ATOMIC)THEN
                  WRITE(LUPRI,'(4X,18(A5,1X))')
     &                 (NAME_SUB2(IBEIG(J)),J=I+1,I+JMAX)
                  WRITE(LUPRI,'(2X,A,18(A5,1X))')
     &                 'Mj',(NAME_SUB3(IBEIG(J)),J=I+1,I+JMAX)
                ELSE 
                  WRITE(LUPRI,'(4X,18(A5))')
     &                 (NAME_SUB2(IBEIG(J)),J=I+1,I+JMAX)
                ENDIF
                TMP = TMP - 18                   
              ENDDO              
            END IF
C...........Negative energy shells
            IF (PRPOS .AND. NPSH(IFSYM) .GT. 0) THEN
              WRITE(LUPRI,'(2X,A)') '* Negative energy orbitals'
              JST = IORB(IFSYM)
              TMP = NPSH(IFSYM)          
              DO I = JST,JST+NPSH(IFSYM)-1,18
                IF(TMP.GT.18)THEN
                  JMAX = 18
                ELSE
                  JMAX = TMP
                ENDIF 
                IF(ATOMIC)THEN
                  WRITE(LUPRI,'(4X,18(A5,1X))')
     &                      (NAME_SUB2(IBEIG(J)),J=I+1,I+JMAX)
                  WRITE(LUPRI,'(2X,A,18(A5,1X))')
     &                 'Mj',(NAME_SUB3(IBEIG(J)),J=I+1,I+JMAX)
                ELSE
                  WRITE(LUPRI,'(4X,18(A5))')
     &                 (NAME_SUB2(IBEIG(J)),J=I+1,I+JMAX)
                ENDIF
                TMP = TMP - 18                
              ENDDO              
            END IF !(PRPOS .AND. NPSH(IFSYM) .GT. 0)
          END DO ! IFSY
        END IF !  IF (NBSYM.GT.1) THEN
C
C     Print out occupation of subblocks
C
        WRITE(LUPRI,'(/A)') '* Occupation of subblocks' 
        DO IFSYM=1,NFSYM
          TLINE = ' '
          IF(NORB(IFSYM).GT.0)THEN
          IF(ATOMIC)THEN
             WRITE(TLINE2,'(A3,A2)') FREP(IFSYM),': '
             WRITE(TLINE3,'(A)') ' Mj: '
          ELSE   
             WRITE(TLINE ,'(A3,A2)') FREP(IFSYM),': '
          ENDIF   
          ILOW = 0
          ITOP = N_SUB_BL(IFSYM)
 10       CONTINUE
          ITOP = MIN(ITOP,NCOL)+ILOW
          ILOW = ILOW+1
          IT = 6
          DO ISUB = ILOW, ITOP
            ID = ID_SUB_BL(ISUB,IFSYM)
            IF (SPINFR) THEN
              WRITE(TLINE(IT:),'(1X,A)') REP(ID)
            ELSEIF(ATOMIC) THEN
              WRITE(TLINE2(IT:),'(1X,A5)') NAME_SUB2(ID)
              WRITE(TLINE3(IT:),'(1X,A5)') NAME_SUB3(ID)
              IT = IT + 1
            ELSE
              WRITE(TLINE(IT:),'(A5)') NAME_SUB2(ID)
            ENDIF  
            IT = IT + 5
          ENDDO
          IF(ATOMIC)THEN
            WRITE (LUPRI,'(23X,A)') TLINE2
            WRITE (LUPRI,'(23X,A)') TLINE3
          ELSE  
            WRITE (LUPRI,'(23X,A)') TLINE
          ENDIF
          IF(ATOMIC)THEN
            WRITE(LUPRI,'(A,F6.4,A,16I6)') 
     &       '  closed shells (f=',DF(0),'):',
     &       (ISHEL_SUB(1,ISUB,IFSYM),ISUB=ILOW,ITOP)
            DO IOPEN=1,NOPEN
              WRITE(LUPRI,'(A,I1,A,F6.4,A,16I6)')
     &        '  open shell #',IOPEN,' (f=',DF(IOPEN),'):',
     &        (ISHEL_SUB(1+IOPEN,ISUB,IFSYM),ISUB=ILOW,ITOP)
            ENDDO
            WRITE(LUPRI,'(A,F6.4,A,16I6)')
     &       ' virtual shells (f=',D0,'):',
     &       (ISHEL_SUB(1+MXOPEN+1,ISUB,IFSYM),ISUB=ILOW,ITOP)
            IF (PRPOS.AND.(.NOT.LEVYLE)) THEN
              WRITE(LUPRI,'(A,F6.4,A,16I6)')
     &        ' neg.erg shells (f=',D1,'):',
     &       (ISHEL_SUB(1+MXOPEN+2,ISUB,IFSYM),ISUB=ILOW,ITOP)
            ENDIF
C
            WRITE(LUPRI,'(A,16I6)')
     &        'tot.num. of pos.erg shells:',
     &      (NORB_SUB(ISUB,IFSYM,1),ISUB=ILOW,ITOP)
            IF (PRPOS.AND.(.NOT.LEVYLE)) THEN
            WRITE(LUPRI,'(A,16I6)')
     &        'tot.num. of neg.erg shells:',
     &      (NORB_SUB(ISUB,IFSYM,2),ISUB=ILOW,ITOP)
            ENDIF
            ILOW = ITOP
            ITOP = N_SUB_BL(IFSYM)-ITOP
            IF(ITOP.GT.0) GOTO 10
          ELSE
            WRITE(LUPRI,'(A,F6.4,A,16I5)') 
     &       '  closed shells (f=',DF(0),'):',
     &       (ISHEL_SUB(1,ISUB,IFSYM),ISUB=ILOW,ITOP)
            DO IOPEN=1,NOPEN
              WRITE(LUPRI,'(A,I1,A,F6.4,A,16I5)')
     &        '  open shell #',IOPEN,' (f=',DF(IOPEN),'):',
     &        (ISHEL_SUB(1+IOPEN,ISUB,IFSYM),ISUB=ILOW,ITOP)
            ENDDO
            WRITE(LUPRI,'(A,F6.4,A,16I5)')
     &       ' virtual shells (f=',D0,'):',
     &       (ISHEL_SUB(1+MXOPEN+1,ISUB,IFSYM),ISUB=ILOW,ITOP)
            IF (PRPOS.AND.(.NOT.LEVYLE)) THEN
              WRITE(LUPRI,'(A,F6.4,A,16I5)')
     &        ' neg.erg shells (f=',D1,'):',
     &       (ISHEL_SUB(1+MXOPEN+2,ISUB,IFSYM),ISUB=ILOW,ITOP)
            ENDIF
C
            WRITE(LUPRI,'(A,16I5)')
     &        'tot.num. of pos.erg shells:',
     &      (NORB_SUB(ISUB,IFSYM,1),ISUB=ILOW,ITOP)
            IF (PRPOS.AND.(.NOT.LEVYLE)) THEN
            WRITE(LUPRI,'(A,16I5)')
     &        'tot.num. of neg.erg shells:',
     &      (NORB_SUB(ISUB,IFSYM,2),ISUB=ILOW,ITOP)
            ENDIF
            ILOW = ITOP
            ITOP = N_SUB_BL(IFSYM)-ITOP
            IF(ITOP.GT.0) GOTO 10
          ENDIF !ATOMIC
          ENDIF !NORB(IFSYM).GT.0
        ENDDO !IFSYM
        WRITE(LUPRI,'(/)')

        IC =1
        IF (PRPOS) IC=0
        DO IFSYM=1, NFSYM
          IF(NORB(IFSYM).GT.0)THEN 
          DO ISUB=1,N_SUB_BL(IFSYM)
C.........closed shells
            ISUM = ISHEL_SUB(1,ISUB,IFSYM)
C.........add open shells
            DO IOPEN=1,NOPEN
              ISUM = ISUM + ISHEL_SUB(1+IOPEN,ISUB,IFSYM)
            ENDDO
C.........add virtual shells
            ISUM = ISUM + ISHEL_SUB(1+MXOPEN+1,ISUB,IFSYM)
C.........add positronic shells
            ISUM = ISUM + ISHEL_SUB(1+MXOPEN+2,ISUB,IFSYM)
C.........note: IC=0 when print of positr. shells desired
            IF (ISUM.NE.NORB_SUB(ISUB,IFSYM,IC)) THEN
              WRITE(LUPRI,'(2X,A,A,A,I3)') 
     &          'WARNING! PREIGN: Inconsistency for ifsym=',
     &        FREP(IFSYM),' and  mj_id=',ID_SUB_BL(ISUB,IFSYM)
              WRITE(LUPRI,'(A,I4,A,I4)')
     &          'ISUM=',ISUM,'<> NORB_SUB(ISUB,IFSYM,IC)=',
     &          NORB_SUB(ISUB,IFSYM,IC)
           ENDIF
          ENDDO
          ENDIF !NORB(IFSYM).GT.0        
        ENDDO
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck preign2 */
      SUBROUTINE PREIGN2(EIG,IFIRST,NMULT,NSH,FMT)
C***********************************************************************
C
C     Prints eigenvalues - new version.
C
C     Written by J. Thyssen - 26 Aug 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
      DIMENSION EIG(*),IFIRST(*),NMULT(*)
      CHARACTER FMT*6
C.s/sya,2006.12.02 - SK - 30-11-2009
      CHARACTER*26 FMT5
C.q
C
C     Find degeneracies
C
      CALL DEGANA(NDEG,NSH,EIG,IFIRST,NMULT)
C
C     Print section
C
C.s/sya,2006.11.30 - SK - 30-11-2009
C#    WRITE(LUPRI,'(5('//FMT//',2X,A1,I2,A1,3X))')
C#   &   (EIG(IFIRST(I)),'(',
C#   &   2*NMULT(I),')',I=1,NDEG)
      FMT5 = '(5(' // FMT // ',2X,A1,I2,A1,3X))'
      WRITE(LUPRI,FMT5)
     &   (EIG(IFIRST(I)),'(',
     &   2*NMULT(I),')',I=1,NDEG)
C.q
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pckeig */
      SUBROUTINE PCKEIG(EIG,ISTART,IEND,IBEIG,ISYM,NEIG,EIGBUF)
C***********************************************************************
C
C     Pick all eigenvalues in the from the array 
C     EIG(ISTART),...,EIG(IEND) and put them in EIGBUF
C
C     Written by J. Thyssen - 26 Aug 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION EIG(*),IBEIG(*),EIGBUF(*)
C
      NEIG = 0
      DO I = ISTART,IEND
         IF (IBEIG(I) .EQ. ISYM) THEN
            NEIG = NEIG + 1
            EIGBUF(NEIG) = EIG(I)
         END IF
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck eighlc */
      SUBROUTINE CHKHLG(EIG,IBEIG)
C***********************************************************************
C
C     Check HOMO-LUMO gap
C
C     Written by J. Thyssen - 26 Aug 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (GAPMIN = 0.1D00)
C
      DIMENSION EIG(*),IBEIG(*)
      CHARACTER CLUMO*3,CHOMO*3
C
#include "dgroup.h"
#include "pgroup.h"
#include "dcborb.h"
#include "dcbham.h"
C

!     initialize EHOMO and ELUMO
!     Below my pragmatic solution that
!     should work in most conceivable cases. LV, 1-9-2000.
      ELUMO = 1.D15
      EHOMO = -1.D15

C     Initialize (return of no HOMO or LUMO)
C
      IHOMO = -1
      DO I = 1,NFSYM
        IF(NOCC(I).GT.0) THEN
          IF (SPINFR) THEN
            IHOMO = IBEIG(IORB(I)+NPSH(I)+1)
          ELSE
           IHOMO = I
          ENDIF
          GOTO 10
        ENDIF
      ENDDO
 10   CONTINUE
      IF(IHOMO.LT.0) RETURN
      ILUMO = -1
      DO I = 1,NFSYM
        IF(NSSH(I).GT.0) THEN
          IF (SPINFR) THEN
            ILUMO = IBEIG(IORB(I)+NPSH(I)+NOCC(I)+1)
          ELSE
            ILUMO = I
          ENDIF
          GOTO 20
        ENDIF
      ENDDO
 20   CONTINUE
      IF(ILUMO.LT.0) RETURN
      DO I = 1,NFSYM
         DO J = 1,NSSH(I)
            IF (ELUMO .GT. EIG(IORB(I)+NPSH(I)+NOCC(I)+J)) THEN
               ELUMO = EIG(IORB(I)+NPSH(I)+NOCC(I)+J)
               IF (SPINFR) THEN
                  ILUMO = IBEIG(IORB(I)+NPSH(I)+NOCC(I)+J)
               ELSE
                  ILUMO = I
               END IF
            END IF
         END DO
         DO J = 1,NOCC(I)
            IF (EHOMO .LT. EIG(IORB(I)+NPSH(I)+J)) THEN
               EHOMO = EIG(IORB(I)+NPSH(I)+J)
!radovan: but if levyle == .true. also spinfr == .true.
               IF (LEVYLE .OR. SPINFR) THEN
                  IHOMO = IBEIG(IORB(I)+NPSH(I)+J)
               ELSE
                  IHOMO = I
               END IF
            END IF
         END DO
      END DO
C
C     HOMO - LUMO gap
C
      IF (SPINFR) THEN
         CHOMO = REP(IHOMO)
         CLUMO = REP(ILUMO)
      ELSE
         CHOMO = FREP(IHOMO)
         CLUMO = FREP(ILUMO)
      END IF
      WRITE (LUPRI,'(A)') '* HOMO - LUMO gap:'
      WRITE (LUPRI,'(2(/A,F15.8,A,A,A)/A/A,F15.8,A/)')
     &   '    E(LUMO) :',ELUMO,' au (symmetry ',CLUMO,')',
     &   '  - E(HOMO) :',EHOMO,' au (symmetry ',CHOMO,')',
     &   '  ------------------------------------------',
     &   '    gap     :',ELUMO-EHOMO,' au'
      IF (ELUMO - EHOMO .LT. GAPMIN) THEN
         WRITE(LUPRI,'(/A)')
     &   '* INFO: E(LUMO) - E(HOMO) small or negative.'
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck selown */
      SUBROUTINE SELOWN(IOPT,CMO,EIG,IBEIG,INUC,NNUC,IUNIT,CMOFIL,
     &                  CSEL,ESEL,IBSEL,NVECS,IOFF,KVEC,NSTR,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     Select coeffiecients for symmetry independent nucleus INUC
C     calculated in its own basis
C     Read option is provided from bit-packed IOPT:
C       0001 - give restart info
C       0010 - read coefficients
C       0100 - read eigenvalues
C       1000 - read boson irrep identification
C
C     The array IDIM contains the following information:
C       IDIM(1,IFRP) = NPSH(IFRP)    : number of positronic solutions
C       IDIM(2,IFRP) = NESH(IFRP)    : number of electronic solutions
C       IDIM(3,IFRP) = NFBAS(IFRP,0) : number of AO-basis functions
C
C     modified for MCSCF coefficients, occupation numbers and mj-values 
C     by S. Knecht - April 2010
C
C     Written by T.Saue Sep 25 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "nuclei.h"
      DIMENSION EIG(*),CMO(*),IBEIG(*),NSTR(2,0:2),KVEC(2),NVECS(2),
     &          CSEL(*),ESEL(*),IBSEL(*),IOFF(2),WORK(*)
C     Local variables
      LOGICAL LBIT,FNDLAB
      CHARACTER TEXT*74,CMOFIL*6,FMT*6,MXFORM*6
      DIMENSION IDIM(3,2),JDIM(2,0:3)
      REAL*8 CMCSCF_MAGIC
!
!     flag for MCSCF coefficients on file CMOFIL
      CMCSCF_MAGIC = 137.0D0
C
      KFRSAV = KFREE
      REWIND IUNIT
C
C     Read title line
C     ===============
C
      IF(FNDLAB('INFO    ',IUNIT)) THEN      
        READ (IUNIT,END=10,ERR=20) TEXT,NSYM,NZBUF,
     &          ((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG
      ELSE
         CALL QUIT('SELOWN: Error reading info')
      ENDIF
      
!
!     set local IOPT_TMP
      IOPT_TMP = IOPT
!
!     found MCSCF coefficients on file CMOFIL?
      IF( TOTERG .eq. CMCSCF_MAGIC ) IOPT_TMP = 14
C
C     File info
C     =========
C
      IF(LBIT(IOPT_TMP,1)) THEN
        FMT = MXFORM(TOTERG,20)
        WRITE(LUPRI,'(/A,A6,3X,A,'//FMT//')')
     &   'SELOWN: Coefficients read from file ',CMOFIL,
     &   '- Total energy: ',TOTERG
        WRITE(LUPRI,'(2A)') '* Heading :',TEXT
        WRITE(LUPRI, '(2A, I8)') '- number of symmetry independent',  
     &   ' nuclei in this fragment : ', NNUC      
        WRITE(LUPRI,'(A,2I8)')
     &   '- Positrons : ',(IDIM(1,I),I=1,NFSYM)
        WRITE(LUPRI,'(A,2I8)')
     &   '- Electrons : ',(IDIM(2,I),I=1,NFSYM)
        WRITE(LUPRI,'(A,2I8)')
     &   '- SO-basis  : ',(IDIM(3,I),I=1,NFSYM)
      ENDIF
C
C     Calculate dimensions
C     ====================
C
      NCDIM = 0
      NEDIM = 0
      DO IFRP = 1,NFSYM
        JDIM(IFRP,1) = IDIM(1,IFRP)
        JDIM(IFRP,2) = IDIM(2,IFRP)
        JDIM(IFRP,0) = IDIM(1,IFRP) + IDIM(2,IFRP)
        JDIM(IFRP,3) = IDIM(3,IFRP)                 
        NEDIM = NEDIM + (IDIM(1,IFRP)+IDIM(2,IFRP))
        NCDIM = NCDIM + (IDIM(1,IFRP)+IDIM(2,IFRP))*IDIM(3,IFRP)
      ENDDO
      NCDIM = NCDIM*NZ
C
C     Read coefficients
C     =================
C    
      IF(LBIT(IOPT_TMP,2)) THEN
        CALL REACMO_coefficients(IUNIT,CMO,NCDIM,NFSYM,NZ,
     &           JDIM(1,0),JDIM(1,3),JDIM(1,1),JDIM(1,2),IDIM)         
      ENDIF
C
C     Read eigenvalues
C     ================
C
      IF(LBIT(IOPT_TMP,3)) THEN
        CALL REACMO_eigenvalues(IUNIT,EIG,NEDIM,NFSYM,NZ,
     &           JDIM(1,0),JDIM(1,1),JDIM(1,2),IDIM)         
      ENDIF
C
C     Read boson irrep info
C     =====================
C
      IF(LBIT(IOPT_TMP,4)) THEN
        CALL REACMO_supersymmetry(IUNIT,'SUPERSYM',IBEIG,NEDIM,NFSYM,NZ,
     &           JDIM(1,0),JDIM(1,1),JDIM(1,2),IDIM)         
      ENDIF
C
C     Select and adjust format to molecular basis
C     ===========================================
C
      ICMO  = 1
      IEIG  = 1
      ICOFF = 1
      IEOFF = 1
      DO IFRP = 1,NFSYM
      IF(NSTR(IFRP,0).GT.0) THEN
        NCMO = IDIM(1,IFRP)+IDIM(2,IFRP)
C       Generate pointer array for bases
        CALL MEMGET('INTE',KBAS,IDIM(3,IFRP),WORK,KFREE,LFREE)
        CALL SELOWI(INUC,NNUC,IFRP,WORK(KBAS),NSBAS)
        IF(NSBAS.NE.IDIM(3,IFRP)) THEN
          WRITE(LUPRI,'(A,I3/A,A4,A,I3)')
     &      'SELOWN: Error in selection of coefficients in ircop ',
     &      IFRP,
     &      'for symmetry independent center ',NAMN(INUC),' no. ',INUC
          WRITE(LUPRI,'(A,I8)')
     &      'Number of basis functions for this center is:',NSBAS,
     &      'Number of basis functions in coefficient file:',
     &       IDIM(3,IFRP)
           CALL QUIT('SELOWN: Error in cf.selection !')
        ENDIF
C       Select coefficients
        IF(LBIT(IOPT_TMP,2)) THEN
          ICSEL = ICOFF + NFBAS(IFRP,0)*IOFF(IFRP)
          CALL SELOWC(NZ,CMO(ICMO),IDIM(3,IFRP),NCMO,
     &              IDIM(1,IFRP),IDIM(2,IFRP),
     &              CSEL(ICSEL),NFBAS(IFRP,0),NVECS(IFRP),
     &              WORK(KVEC(IFRP)),WORK(KBAS),
     &              NSTR(IFRP,2),NSTR(IFRP,1))
          ICMO = ICMO + NCMO*IDIM(3,IFRP)*NZ
        ENDIF
        IESEL = IEOFF + IOFF(IFRP)
C       Select eigenvalues
        IF(LBIT(IOPT_TMP,3)) THEN
          CALL SELOWE(EIG(IEIG),IDIM(1,IFRP),IDIM(2,IFRP),
     &                ESEL(IESEL),WORK(KVEC(IFRP)),
     &                NSTR(IFRP,2),NSTR(IFRP,1))
        ENDIF
C       Select boson irrep information
        IF(LBIT(IOPT_TMP,4)) THEN
          CALL SELOWB(IBEIG(IEIG),IDIM(1,IFRP),IDIM(2,IFRP),
     &                IBSEL(IESEL),WORK(KVEC(IFRP)),
     &                NSTR(IFRP,2),NSTR(IFRP,1))
        ENDIF
        IEIG = IEIG + NCMO
        CALL MEMREL('SELOWN.loop',WORK,KBAS,KBAS,KFREE,LFREE)
        ICOFF = ICOFF + NFBAS(IFRP,0)*NVECS(IFRP)*NZ
        IEOFF = IEOFF + NVECS(IFRP)
      ENDIF
      ENDDO
      RETURN
 10   CONTINUE
      CALL QUIT('SELOWN: END reading TEXT')
 20   CONTINUE
      CALL QUIT('SELOWN: ERROR reading TEXT')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck selow2 */
      SUBROUTINE SELOWC(NZ,CMO,KBAS,KORB,KPSH,KESH,
     &                     CBUF,MBAS,MORB,
     &                     JVEC,JBAS,NPVEC,NEVEC)
C***********************************************************************
C     Pick out a set of vectors from CMO according to array JVEC
C
C     Written by T.Saue 1997
C     Last revision Jan 8 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION CMO(KBAS,KORB,NZ),CBUF(MBAS,MORB,NZ),JVEC(*),
     &          JBAS(*)
C
      NVEC = NPVEC + NEVEC
      NTOT = NVEC*MBAS
      DO IZ = 1,NZ
        CALL DZERO(CBUF(1,1,IZ),NTOT)
C
C       Positronic vectors
C
        DO I = 1,NPVEC
          II = KPSH+1+JVEC(I)
          DO J = 1,KBAS
            CBUF(JBAS(J),I,IZ) = CMO(J,II,IZ)
          ENDDO
        ENDDO
C
C       Electronic vectors
C
        DO I = NPVEC+1,NVEC
          II = KPSH+JVEC(I)
          DO J = 1,KBAS
            CBUF(JBAS(J),I,IZ) = CMO(J,II,IZ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck selowe */
      SUBROUTINE SELOWE(EIG,KPSH,KESH,EBUF,JVEC,NPVEC,NEVEC)
C***********************************************************************
C     Pick out a set of eigenvalues from EIG according to array JVEC
C
C     Written by L. Visscher 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION EIG(*),EBUF(NPVEC+NEVEC),JVEC(*)
C
      NVEC = NPVEC + NEVEC
C
C
C     Positronic eigenvalues
C
      DO I = 1,NPVEC
        II = KPSH+1+JVEC(I)
        EBUF(I) = EIG(II)
      ENDDO
C
C     Electronic eigenvalues
C
      DO I = NPVEC+1,NVEC
        II = KPSH+JVEC(I)
        EBUF(I) = EIG(II)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck selowb */
      SUBROUTINE SELOWB(IBEIG,KPSH,KESH,IBBUF,JVEC,NPVEC,NEVEC)
C***********************************************************************
C     Pick out a set of from boson irrep info array IBEIG 
C     according to array JVEC
C
C     Written by L. Visscher 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION IBEIG(*),IBBUF(NPVEC+NEVEC),JVEC(*)
C
      NVEC = NPVEC + NEVEC
C
C
C     Positronic eigenvalues
C
      DO I = 1,NPVEC
        II = KPSH+1+JVEC(I)
        IBBUF(I) = IBEIG(II)
      ENDDO
C
C     Electronic eigenvalues
C
      DO I = NPVEC+1,NVEC
        II = KPSH+JVEC(I)
        IBBUF(I) = IBEIG(II)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck selowi */
      SUBROUTINE SELOWI(INUC,NNUC,IFRP,JBAS,NSBAS)
C***********************************************************************
C
C     Make index array pointing from basis functions of
C     symmetry independent centers INUC..INUC+NNUC-1 to the full
C     molecular basis
C
C     Written by T.Saue Sep 25 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "dcbbas.h"
#include "dcblab.h"
      DIMENSION JBAS(*)
C
#include "dcbibt.h"
      NSBAS = 0
      DO I = 1,NFBAS(IFRP,0)
        II = IBAS(IFRP) + I
        ILAB = IPLAB(II,2)
        ICENT = JGET(IATTR(ILAB,2))
        IF ((ICENT.GE.INUC).AND.(ICENT.LT.(INUC+NNUC))) THEN
          NSBAS = NSBAS + 1
          JBAS(NSBAS) = I
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Pripol */
      SUBROUTINE PRIPOL(CPOL,EPOL,NBAS,NPOL,NZ,LRPOL,LCPOL)
C***********************************************************************
C
C     Print polarization contributions from projection analysis
C
C     Written by T. Saue Oct 25 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,DTOL = 1.0D-5)
C
#include "maxorb.h"
#include "dcblab.h"
      CHARACTER FMT*6,MXFORM*6
      DIMENSION CPOL(LRPOL,LCPOL,NZ),EPOL(NPOL)
C
      DO I = 1,NPOL
        FMT = MXFORM(EPOL(I),16)
        WRITE(LUPRI,'(/A,I3,A,'//FMT//')')
     +    '* Electronic eigenvalue nr.',I,': ',EPOL(I)
        CALL PRSYMB(LUPRI,'=',52,0)
        BIG=D0
        IBIG=0
        DO J = 1,NBAS
          DVEC = D0
          DO IZ = 1,NZ
            DVEC = DVEC + ABS(CPOL(J,I,IZ))
          ENDDO
          IF(DVEC.GT.DTOL) THEN
            WRITE(LUPRI,100)
     &          J,PLABEL(IPLAB(J,2),2),(CPOL(J,I,IZ),IZ = 1,NZ)
            IF(DVEC.GT.BIG) THEN
              BIG  = DVEC
              IBIG = J
            ENDIF
          ENDIF
        ENDDO
        WRITE(6,*) ' * Max. ',BIG
        IF(IBIG.GT.0) THEN
          WRITE(LUPRI,100)
     &          IBIG,PLABEL(IPLAB(IBIG,2),2),(CPOL(IBIG,I,IZ),IZ = 1,NZ)
        ENDIF
      ENDDO
      CALL PRSYMB(LUPRI,'=',52,0)
      RETURN
 100  FORMAT(3X,I5,2X,A12,2X,4F14.10)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ORBCN2(JVEC,NVEC,IFRP,NPVEC,NEVEC)
C***********************************************************************
C     Find number of electronic/positronic vectors
C
C     Written by T.Saue Jan 7 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION JVEC(NVEC)
#include "dcborb.h"
C
      NPVEC = 0
      NEVEC = 0
      II = 0
      DO 10 I = 1,NVEC
        J = JVEC(I)
C
C       Positronic orbital in range [-NPSH(IFRP),-1]
C
        IF    (J.LT.0.AND.J.GE.(-NPSH(IFRP))) THEN
          JVEC(I)  = NPSH(IFRP)+1+J
          NPVEC    = NPVEC + 1
C
C       Electronic orbital in range [1,NESH(IFRP)]
C
        ELSEIF (J.GT.0.AND.J.LE.NESH(IFRP)) THEN
          JVEC(I)  = NPSH(IFRP)+J
          NEVEC    = NEVEC + 1
        ENDIF
 10   CONTINUE
C
      NVEC = NPVEC+NEVEC
      RETURN
      END
