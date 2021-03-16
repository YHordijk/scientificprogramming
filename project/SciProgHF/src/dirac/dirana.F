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
C  /* Deck pamana */
      SUBROUTINE PAMANA()
C***********************************************************************
C
C     Analyze Dirac-Fock wavefunction
C
C     Written by T.Saue Sep 11 1995
C     Last revision : KApril 27 1996
C
C***********************************************************************

      use dirac_cfg
      use second_order_minimization
      use memory_allocator
      use xmlout

#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbbas.h"
#include "dcbrho.h"
      real(8), allocatable :: WORK(:)
C Used from COMMON blocks:
C  CBIANA: DOPOP,DOANA,DOEXP
C  DCBBAS: N2BBASX
C
C
      CALL QENTER('PAMANA')

      call legacy_lwork_get(LWORK)
      call alloc(WORK,LWORK,id='WORK in PAMANA')
C
C     Write output in xml format
C
      if (doxml) then
         call xml_begin('task','type="ANALYSIS"')
         call xml_begin('output')
      endif
C
C***********************************************************************
C*****  V E C T O R    P R I N T ***************************************
C***********************************************************************
C
      IF(DOVEC) CALL ANAVEC(WORK,LWORK)
C
C***********************************************************************
C*****  M U L L I K E N    P O P U L A T I O N     A N A L Y S I S  ****
C***********************************************************************
C
      IF(DOPOP) CALL MULPOP(WORK,LWORK)
C
C***********************************************************************
C*****  P L O T   1 - D I M E N S I O N A L     D E N S I T Y   ********
C***********************************************************************
C
      IF(DO1RHO) CALL GT1RHO(WORK,LWORK)
C
C***********************************************************************
C*****  P L O T   3 - D I M E N S I O N A L     D E N S I T Y   ********
C***********************************************************************
C
      IF (DO3RHO) THEN
         IF (NCUBORB(1).EQ.0.AND.NCUBORB(2).EQ.0) THEN
C     Total density
            CALL GT3RHO(WORK,LWORK)
         ELSE
C     Density of individual orbitals
            CALL WRITE_ORB_CUBES(NCUBORB(1),NCUBORB(2),
     &           ICUBORB,WORK,LWORK)
         ENDIF
      ENDIF
C
C***********************************************************************
C*****  P R O J E C T I O N   O N T O    A N O T H E R   S O L U T I O N
C***********************************************************************
C
      IF(DOPRJ) CALL ANAPRJ(WORK,LWORK)
C
C
C***********************************************************************
C*****  P L O T  I N T E G R A T I O N                   ***************
C*****  W E I G H T S   I N   O N E   D I M E N S I O N  ***************
C***********************************************************************

      IF(DO1WT) CALL GT1WT(WORK,LWORK)

C
C***********************************************************************
C*************  O R B I T A L   L O C A L I S A T I O N  ***************
C***********************************************************************
C
      IF(DOLOC)THEN
        CALL LOCALIZATION_INTERFACE
        CALL MINIMIZATION_DRV
      ENDIF
C
C***********************************************************************
C
C
#ifdef ANALYZE_PROPERTY_GRADIENT
      CALL ANALYZE_PROPERTY_GRAD(WORK,LWORK)
#endif
C
      if (doxml) then
         call xml_end('output')
         call xml_end('task')
      endif

      call dealloc(WORK)
C
      CALL QEXIT('PAMANA')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck scf_mulpop */
      SUBROUTINE SCF_MULPOP(CMO,EIG,IBEIG,WORK,LWORK)
C***********************************************************************
C
C     Mulliken population analysis in each SCF iteration
C
C     Written by H.J.Aa.Jensen and P.Norman, Jan. 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
C Used from COMMON blocks:
C  DCBBAS: N2BBASX
C  LABELS: NLAB
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbpop.h"
#include "dgroup.h"
C
      DIMENSION CMO(*), EIG(*), IBEIG(*), WORK(LWORK)
      LOGICAL   DONETP_M,LABDEF_M
C
      CALL QENTER('SCF_MULPOP')
#include "memint.h"
C
C     Save possible user input (for default Mull. pop) and set
C     SCF_MULPOP default values
C
      IPRPOP_M = IPRPOP
      ILABDF_M = ILABDF
      DONETP_M = DONETP
      LABDEF_M = LABDEF
C
      IPRPOP   = -1 ! suppress normal MULPOP output
      ILABDF   =  1 ! select AO basis for easy summing for an atom
      DONETP   = .FALSE. ! no net populations
      LABDEF   = .FALSE. ! labels not already defined
      DO_SCFPOP = .TRUE. ! select SCFPOP special output;
                         ! implies only total charges for each atom, omit point charges
C
      KRMC_FLG = 0
C
      DO I = 1, NFSYM
         VECPOP_SAVE(I) = VECPOP(I)
         NOCCI = NOCC(I)
         IF (NOCCI .eq. 0) THEN
            VECPOP(I) = '   '
         ELSE IF( NOCCI .eq. 1 )THEN
            VECPOP(I) = '1  '
         ELSE
            WRITE(VECPOP(I),'(A3,I4)') '1..',NOCCI
         END IF
      END DO
C
C     Dummy memory allocation
      CALL MEMGET('REAL',KOCC ,0       ,WORK,KFREE,LFREE)
C
C     Get labels for population analysis
C
      CALL GETPLB(IPRPOP)
C
C     Memory allocation
      NNLAB = (NPOPLAB*(NPOPLAB+1))/2
      N2LAB = NPOPLAB*NPOPLAB
      CALL MEMGET('INTE',KIPGR ,NTBAS(0)     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT ,N2LAB*2      ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT ,N2BBASX      ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPOPG ,NPOPLAB*2    ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIEFF ,NPOPLAB      ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPGTOT,NPOPLAB*NFSYM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPNTOT,N2LAB*NFSYM  ,WORK,KFREE,LFREE)
C
      CALL MULPOP_1(WORK(KPMAT),WORK(KSMAT),CMO,EIG,
     &            WORK(KOCC),IBEIG,WORK(KIPGR),WORK(KPOPG),
     &            WORK(KIEFF),WORK(KPGTOT),WORK(KPNTOT),WORK(KFREE),
     &            LFREE,KRMC_FLG)
C
C     Memory deallocation
      CALL MEMREL('SCF_MULPOP',WORK,KWORK,KWORK,KFREE,LFREE)
C     restore user input for Mull. pop. ana.
      IPRPOP   = IPRPOP_M
      ILABDF   = ILABDF_M
      DONETP   = DONETP_M
      LABDEF   = LABDEF_M
      DO I = 1, NFSYM
         VECPOP(I) = VECPOP_SAVE(I)
      END DO
      DO_SCFPOP = .FALSE.
C
      CALL QEXIT('SCF_MULPOP')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mulpop */
      SUBROUTINE MULPOP(WORK,LWORK)
C***********************************************************************
C
C     Mulliken population analysis
C
C     Written by T.Saue - Sep 11 1995
C     Last revision: April 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
C Used from COMMON blocks:
C  DCBBAS: N2BBASX
C  LABELS: NLAB
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbpop.h"
#include "dgroup.h"
#include "dcbmp2no.h"
C
      LOGICAL TOBE,TOBEK,DONETP_M,ADDSML_M,LABDEF_M
      LOGICAL FNDLAB,CINATORB
      LOGICAL MJVEC_AVAIL
      CHARACTER FILELABEL*8
      DIMENSION WORK(LWORK)
C
      CALL QENTER('MULPOP')
#include "memint.h"
C
      CALL TITLER('Mulliken population analysis','*',122)
C
      DO_SCFPOP = .FALSE.
C
      IF( MP2NATPOP )THEN
C       set-up analysis string
        DO I = 1, NFSYM
          VECPOP_SAVE(I) = VECPOP(I)
          NOCCI = NACTMP2NO(I)
          IF( SELPOPNAT .eq. 'all' ) NOCCI = NESH(I)
          IF (NOCCI .eq. 0) THEN
            VECPOP(I) = '   '
          ELSE IF( NOCCI .eq. 1 )THEN
            VECPOP(I) = '1  '
          ELSE
            WRITE(VECPOP(I),'(A3,I4)') '1..',NOCCI
          END IF
        END DO
C       save possible user input (for default Mull. pop) and set
C       MP2-NO default values
        IPRPOP_M = IPRPOP
        ILABDF_M = ILABDF
        DONETP_M = DONETP
        ADDSML_M = ADDSML
        LABDEF_M = LABDEF
C
        IPRPOP   = IPRMP2NO
        ILABDF   = 2
        DONETP   = .FALSE.
        ADDSML   = .TRUE.
        LABDEF   = .FALSE.
C
      ELSE
        DO I = 1, NFSYM
          IF (VECPOP(I)(1:9).eq.'UNDEFINED') THEN
            NOCCI = NOCC(I)
            IF (NOCCI .eq. 0) THEN
              VECPOP(I) = '   '
            ELSE IF( NOCCI .eq. 1 )THEN
              VECPOP(I) = '1  '
            ELSE
              WRITE(VECPOP(I),'(A3,I4)') '1..',NOCCI
            END IF
          END IF
        END DO
      END IF
C
C     Check if coefficients are on file
C
      IF( MP2NATPOP ) THEN
         INQUIRE(FILE='KRMCOLD',EXIST=TOBEK)
      ELSE
         INQUIRE(FILE='KRMCSCF',EXIST=TOBEK)
      END IF
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF ((.NOT.TOBE).AND.(.NOT.TOBEK)) GOTO 1000
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KBEIG,NORBT   ,WORK,KFREE,LFREE)
C...  may be oversized when not all orbitals are investigated
      CALL MEMGET('REAL',KOCC ,NORBT   ,WORK,KFREE,LFREE)
C
      CALL IZERO(WORK(KBEIG),NORBT)
      CALL DZERO(WORK(KCMO),N2BBASXQ)
      CALL DZERO(WORK(KOCC),NORBT)
C
C     Read coefficients, eigenvalues and
C     irrep identification + restart info
C
      IF (TOBEK) THEN
         IF( MP2NATPOP )THEN
           CALL OPNFIL(LUKRMC,'KRMCOLD','OLD','MULPOP')
           KRMC_FLG = 2
           CALL MEMGET('REAL',KEIG_NAT,NESHT,WORK,KFREE,LFREE)
C          ... initialize
           CALL DZERO(WORK(KEIG_NAT),NESHT)
C          ... read from file
           CALL READMP2NAT(WORK(KCMO),WORK(KOCC),WORK(KEIG_NAT),
     &                     WORK(KBEIG))
           WRITE(LUPRI,'(/A)')
     &        ' (MULPOP) MP2 natural orbitals and occupancies read'//
     &        ' from label NEWNATOB/NEWNATOC on file KRMCOLD'
C
           CLOSE(LUKRMC,STATUS='KEEP')
         ELSE
           CINATORB = .FALSE.
           MJVEC_AVAIL = .FALSE.
           CALL OPNFIL(LUKRMC,'KRMCSCF','OLD','MULPOP')
           REWIND(LUKRMC)
           FILELABEL = 'MCCINATO'
           CINATORB = FNDLAB(FILELABEL,LUKRMC)
           JRDMO = -1
           REWIND(LUKRMC)
           FILELABEL = 'MJVEC   '
           MJVEC_AVAIL = FNDLAB(FILELABEL,LUKRMC)
           IF(CINATORB)THEN
             KRMC_FLG = 3
C            coefficients
             CALL RREADMO(WORK(KCMO),JRDMO,4,LUKRMC)
             IF (JRDMO .EQ. 0) THEN
                WRITE(LUPRI,'(/A)')
     &               ' (MULPOP)  Orbitals read from' //
     &               ' label MCCINATO on file KRMCSCF'
             END IF
C            occupancies
             CALL REAKRMC(LUKRMC,'MCNATOCC',WORK(KOCC),NORBT)
             WRITE(LUPRI,'(/A)')
     &             ' (MULPOP)  Natural orbital occupancies read from'//
     &             ' label MCNATOCC on file KRMCSCF'
           ELSE
             KRMC_FLG = 1
             CALL RREADMO(WORK(KCMO),JRDMO,1,LUKRMC)
             IF (JRDMO .EQ. 0) THEN
                WRITE(LUPRI,'(/A)')
     &               ' (MULPOP)  Orbitals read from' //
     &               ' label NEWORB on file KRMCSCF'
             END IF
           END IF
           IF( MJVEC_AVAIL )THEN
             CALL IREAKRMC(LUKRMC,'MJVEC   ',WORK(KBEIG),NORBT)
             WRITE(LUPRI,'(/A)')
     &            ' (MULPOP)  mj-values read from' //
     &            ' label MJVEC on file KRMCSCF'
           END IF
         END IF
         CLOSE (LUKRMC)
      ELSE
         KRMC_FLG = 0
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMANA')
         IOPT=14
         IF(ATOMIC) THEN
           CALL ATOMIC_REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),
     &          WORK(KBEIG),TOTERG,IOPT)           
         ELSE
           CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),
     &                 WORK(KBEIG),TOTERG,IOPT)
         ENDIF
         CLOSE(LUCOEF,STATUS='KEEP')
C
C        check for energy 137.0D0 - if found we are doing a Mulliken
C        population analysis of MCSCF/CI natural orbitals in C1
C        symmetry - SK Feb 2009 - this is only a temporary fix. In the
C        long run, occupancies should also be stored on DFCOEF. FIXME!
         TEST_ENERGY = 137.0D0
         IF(TOTERG .eq. TEST_ENERGY )THEN
           CALL DCOPY(NORBT,WORK(KEIG),1,WORK(KOCC),1)
           KRMC_FLG = 3
         END IF
      END IF
C
C     Get labels for population analysis
C
      CALL GETPLB(IPRPOP)
C
C     Memory allocation
      NNLAB = (NPOPLAB*(NPOPLAB+1))/2
      N2LAB = NPOPLAB*NPOPLAB
      CALL MEMGET('INTE',KIPGR ,NTBAS(0)       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT ,N2LAB*2        ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT ,N2BBASX        ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPOPG ,NPOPLAB*2      ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIEFF ,NPOPLAB        ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPGTOT,NPOPLAB*NFSYM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPNTOT,N2LAB*NFSYM  ,WORK,KFREE,LFREE)
C
      CALL MULPOP_1(WORK(KPMAT),WORK(KSMAT),WORK(KCMO),WORK(KEIG),
     &            WORK(KOCC),WORK(KBEIG),WORK(KIPGR),WORK(KPOPG),
     &            WORK(KIEFF),WORK(KPGTOT),WORK(KPNTOT),WORK(KFREE),
     &            LFREE,KRMC_FLG)
C
C     Memory deallocation
      CALL MEMREL('MULPOP',WORK,KWORK,KWORK,KFREE,LFREE)
C     restore user input for Mull. pop. ana.
      IF (MP2NATPOP )THEN
        IPRPOP   = IPRPOP_M
        ILABDF   = ILABDF_M
        DONETP   = DONETP_M
        ADDSML   = ADDSML_M
        LABDEF   = LABDEF_M
        DO I = 1, NFSYM
          VECPOP(I) = VECPOP_SAVE(I)
        END DO
      END IF
      CALL QEXIT('MULPOP')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'MULPOP: Coefficient file not found !'
      CALL QUIT('MULPOP: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MULPOP_1 */
      SUBROUTINE MULPOP_1(PMAT,SMAT,CMO,EIG,OCC,IBEIG,IPGRS,POPG,IEFF,
     &                  PGTOT,PNTOT,WORK,LWORK,KRMC_FLG)
C***********************************************************************
C
C     This routine generates net and gross Mulliken populations
C     separately for each component of a 2-spinor. The array
C     IPLAB determines the grouping of basis functions into labels.
C
C     PNET - net population(triangularly row-packed)
C     PGRS - gross population
C     SMAT - overlap matrix
C     COEF - coefficients
C
C     Written by T.Saue July 24 1994
C     LAST VERSION: Jan 9 1997
C
C***********************************************************************
      use xmlout
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "maxorb.h"
C Used from COMMON blocks:
C  DCBORB: NSPH
C  DCBBAS:
C  DGROUP: NZ
C
#include "dcbgen.h"
#include "dcbpop.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dgroup.h"
C
      DIMENSION PMAT(NPOPLAB,NPOPLAB,2),SMAT(*),CMO(*),
     &          EIG(NORBT),IBEIG(NORBT),
     &          IPGRS(*),POPG(NPOPLAB,2),IEFF(NPOPLAB),OCC(*),
     &          PNTOT(NPOPLAB,NPOPLAB,NFSYM),
     &          PGTOT(NPOPLAB,NFSYM),WORK(LWORK)
C
      CALL QENTER('MULPOP_1')
#include "memint.h"
C
      N2LAB = NPOPLAB*NPOPLAB
      CALL DZERO(PNTOT,N2LAB*NFSYM)
      CALL DZERO(PGTOT,NPOPLAB*NFSYM)
C
C     Make pointers to group labels
C     =============================
C
      DO I = 1,NTBAS(0)
        IPGRS(I) = IPOPLAB(IPLAB(I,ILABDF))
      ENDDO
C
C     Get overlap matrix
C     ==================
C
      CALL GTOVLX(SMAT,SSMTRC)
C
C     If necessary, backtransform overlap matrix to AO-basis
C     ======================================================
C
      IF(ILABDF.EQ.1.AND.NBSYM.GT.1) THEN
        CALL MEMGET('REAL',KBUF ,N2BBASX ,WORK,KFREE,LFREE)
        CALL MTBSBU(SMAT,WORK(KBUF))
        CALL MTSOAO(WORK(KBUF),SMAT,NTBAS(0),0,IPRPOP)
        CALL MEMREL('MULPOP.MULPOP_1',WORK,KWORK,KWORK,KFREE,LFREE)
      ENDIF
C
C     Loop over fermion ircops
C     ========================
C
      DO IFRP = 1,NFSYM
        IF (IPRPOP .GE. 0)
     &     CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
        IF (IPRPOP .GE. 2) WRITE(LUPRI,*) 'vecpop : ',vecpop(ifrp)
C
C       Find number of electronic/positronic vectors
C       ============================================
C
        NVEC = 1
        CALL MEMGET('INTE',KJVEC,NORB(IFRP),WORK,KFREE,LFREE)
        CALL NUMLST(VECPOP(IFRP),WORK(KJVEC),NORB(IFRP),
     &             -NPSH(IFRP),NESH(IFRP),IFRP,NVEC)
        CALL ORBCNT(WORK(KJVEC),NVEC,NPSH(IFRP),NESH(IFRP),
     &              NPVEC,NEVEC)
        IF(NVEC.EQ.0) GOTO 20
        CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
        NCBAS=NTBAS(0)
        IF(ILABDF.EQ.2) NCBAS=NFBAS(IFRP,0)
        CALL MEMGET('REAL',KCBF,NCBAS*NVEC*4,WORK,KFREE,LFREE)
        CALL SELCMO(ILABDF,IFRP,CMO,WORK(KJVEC),NVEC,NPVEC,NEVEC,
     &              WORK(KCBF),NCBAS,NVEC,NCBAS,
     &              IPRPOP,WORK,KFREE,LFREE)
C
C       Perform Mulliken population analysis
C       ====================================
C
        CALL MEMGET('REAL',KOCC2,NVEC,WORK,KFREE,LFREE)
        CALL MULPOP_2(PMAT,IFRP,SMAT,WORK(KCBF),EIG(IORB(IFRP)+1),
     &              IBEIG(IORB(IFRP)+1),
     &              IPGRS,POPG,IEFF,PGTOT,PNTOT,NPVEC,NEVEC,NCBAS,
     &              WORK(KJVEC),OCC,WORK(KOCC2),WORK(KFREE),LFREE,
     &              KRMC_FLG)
        CALL MEMREL('MULPOP_1.occ',WORK,KOCC2,KOCC2,KFREE,LFREE)
 20     CONTINUE
        CALL MEMREL('MULPOP_1',WORK,1,KWORK,KFREE,LFREE)
      END DO ! IFRP = 1, NFSYM
C
C     Write out total populations
C     ===========================
C
      IF (NFSYM.GT.1) CALL DAXPY(NPOPLAB,D1,PGTOT(1,2),1,PGTOT(1,1),1)
      IF (DO_SCFPOP) THEN
!       calculate net atomic charge
        DO I = 1, NPOPLAB
           PGTOT(I,1) = CHRGNUC(I) - PGTOT(I,1)
        END DO
        CALL WR_MULCHRG(PGTOT(1,1),POPLAB,NPOPLAB)
      ELSE
        WRITE(LUPRI,'(//A)') '*** Total gross population ***'
        CALL WRPOPT(PGTOT(1,1),POPLAB,NPOPLAB,IEFF)
        if (doxml) then
           do i = 1, npoplab
              call xml_quantity('Mulliken gross population '//POPLAB(i),
     &                           PGTOT(i,1),'electron')
           end do
        end if
        IF(DONETP) THEN
! Miro: PNTOT's 3rd index is 2, fix needed wrt NFSYM
          IF (NFSYM.GT.1) THEN
            CALL DAXPY(N2LAB,D1,PNTOT(1,1,2),1,PNTOT(1,1,1),1)
          ENDIF
          WRITE(LUPRI,'(//A)') '*** Total net population ***'
          CALL WRPOPN(PNTOT(1,1,1),PGTOT(1,1),
     &                POPLAB,NPOPLAB,IEFF,'total')
        ENDIF
      END IF
C     Memory deallocation
      CALL MEMREL('MULPOP_1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MULPOP_1')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MULPOP_2 */
      SUBROUTINE MULPOP_2(PMAT,IFRP,SMAT,CMO,EIG,IBEIG,IPGRS,POPG,IEFF,
     &           PGTOT,PNTOT,NPVECS,NEVECS,NCBAS,JVEC,OCC,OCC2,WORK,
     &           LWORK,KRMC_FLG)
C***********************************************************************
C
C     This routine generates net and gross Mulliken populations
C     separately for each component of a 2-spinor. The array
C     IPLAB determines the grouping of basis functions into labels.
C
C     PNET - net population(triangularly row-packed)
C     PGRS - gross population
C     SMAT - overlap matrix
C     COEF - coefficients
C
C     KRMC_FLG = 2: MP2 natural orbitals and occupancies read from
C                   label NEWNATOB/NEWNATOC on file KRMCOLD
C     KRMC_FLG = 1: Orbitals read from
C                   label NEWORB on file KRMCSCF
C     KRMC_FLG = 0  Orbitals read from DFCOEF
C
C     Written by T.Saue July 24 1994
C     Revised       : Jan 7 1998 - jth LINUX
C     Last revision : Aug 20 2002 - jkp (take fractional occ. into
C                                        account)
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D2 = 2.0D0,DHALF = 0.5D0)
#include "dummy.h"
#include "maxorb.h"
C Used from COMMON blocks:
C  DCBORB: NSPH
C  DCBBAS:
C  DGROUP: NZ
C  DCBDHF: NOPEN
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbpop.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbdhf.h"
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION SMAT(NBAST,NBAST)
#else
      DIMENSION SMAT(NTBAS(0),NTBAS(0))
#endif
      CHARACTER FMT*6,MXFORM*6,SYMLAB*9,SYMLAB2*11
      CHARACTER SPDCAR*1 ! type declaration for external function      
      DIMENSION PMAT(NPOPLAB,NPOPLAB,2),JVEC(*),OCC(*),OCC2(*),
     &          CMO(NCBAS,NPVECS+NEVECS,NZ),EIG(*),IBEIG(*),IPGRS(*),
     &          POPG(NPOPLAB,2),IEFF(NPOPLAB),
     &          PNTOT(NPOPLAB,NPOPLAB,NFSYM),
     &          PGTOT(NPOPLAB,NFSYM),WORK(LWORK)
      DATA SYMLAB/'         '/
C
      CALL QENTER('MULPOP_2')
#include "memint.h"
C
      N2LAB = NPOPLAB*NPOPLAB
      NVEC = NPVECS + NEVECS
      CALL DZERO(OCC2,NVEC)
      DO I = 1,NVEC
c
c       Find fractional occupation for each requested orbital string
c
        IF (KRMC_FLG.EQ.0) THEN
C       ... SCF-NOs
          IF (JVEC(I) .GT. 0 .AND. JVEC(I) .LE. NISH(IFRP)) THEN
c         ... closed shell
            OCC2(I) = DF(0)
          ELSE
c         ... test open shell(s)
            NVORB =  NISH(IFRP)
            DO IOPEN=1,NOPEN
              NVORB =  NVORB + NACSH(IFRP,IOPEN)
              IF (JVEC(I).GT.NISH(IFRP).AND.JVEC(I).LE.NVORB) THEN
                OCC2(I) = DF(IOPEN)
                GO TO 100
              ENDIF
            END DO
  100       CONTINUE
          END IF
        ELSE IF(KRMC_FLG.EQ.2) THEN
C       ... MP2-NOs
           OCC2(I) = OCC(NPSH(IFRP)+IORB(IFRP)+JVEC(I))
        ELSE IF(KRMC_FLG.EQ.3) THEN
C       ... CI-NOs
           OCC2(I) = OCC(NPSH(IFRP)+IORB(IFRP)+JVEC(I))
        END IF
CSKdebug CALL WRTMATMN(OCC,1,NORBT,1,NORBT,LUPRI)
c
       IF (IPRPOP .GE. 0) THEN
        IF(I.GT.NPVECS) THEN ! Electronic eigenvalues
          II = NPSH(IFRP)+JVEC(I)
          FMT = MXFORM(EIG(II),16)
          IF (KRMC_FLG.gt.0 )THEN
            ID = IBEIG(II)
            WRITE(SYMLAB,'("m_j=",I3,"/2")') ID
            IF (KRMC_FLG.EQ.1) THEN
              IF( ID .eq. 0) THEN 
                WRITE(LUPRI,'(/A,I4,A)')
     &          '* Electronic eigenvalue no.',JVEC(I),': '
              ELSE
                WRITE(LUPRI,'(/A,I4,A,2X,A9)')
     &          '* Electronic eigenvalue no.',JVEC(I),': ',SYMLAB
              END IF
            ELSE
C             ... CI-NOs, MP2-NOs,...
              IF( ID .eq. 0) THEN 
                WRITE(LUPRI,'(/A,I4,A,A,F12.8)')
     &          '* Electronic eigenvalue no.',JVEC(I),': ',
     &          '      Occupation :',OCC2(I)
              ELSE
                WRITE(LUPRI,'(/A,I4,A,A,F12.8,2X,A9)')
     &          '* Electronic eigenvalue no.',JVEC(I),': ',
     &          '      Occupation :',OCC2(I),SYMLAB
              END IF
            END IF
          ELSE
            IF(SPINFR) THEN
              ID = IBEIG(II)
              WRITE(SYMLAB,'("sym= ",A3,1X)') REP(ID)
            ELSEIF(ATOMIC) THEN
              ID = IBEIG(II)
              CALL ATOMIC_ID(ID,KP,J,MJ,LL)
              WRITE(SYMLAB2,'(A,I2,"/2",";",I3,"/2")') SPDCAR(LL),J,MJ
            ELSEIF(LINEAR) THEN
              ID = IBEIG(II)
              WRITE(SYMLAB,'("m_j=",I3,"/2")') ID
            ENDIF
C           
            IF(ATOMIC)THEN
              WRITE(LUPRI,'(/A,I4,A,'//FMT//',A,F6.4,A,2X,A11)')
     &           '* Electronic eigenvalue no.',JVEC(I),': ',
     &           EIG(II),'   (Occupation : f = ',OCC2(I),')',
     &           SYMLAB2
            ELSE
              WRITE(LUPRI,'(/A,I4,A,'//FMT//',A,F6.4,A,2X,A9)')
     &           '* Electronic eigenvalue no.',JVEC(I),': ',
     &           EIG(II),'       (Occupation : f = ',OCC2(I),')',
     &           SYMLAB
            ENDIF   
          END IF   !  IF (KRMC_FLG.gt.0 )
        ELSE ! Positronic eigenvalues
          II  = NPSH(IFRP)+1+JVEC(I)
          FMT = MXFORM(EIG(II),16)
          IF (KRMC_FLG.GT.0) THEN
             WRITE(LUPRI,'(/A,I3)')
     &       '* Positronic eigenvalue no.',II
          ELSE
            IF(SPINFR) THEN
              ID = IBEIG(II)
              WRITE(SYMLAB,'("sym= ",A3,1X)') REP(ID)
            ELSEIF(ATOMIC) THEN
              ID = IBEIG(II)
              CALL ATOMIC_ID(ID,KP,J,MJ,LL)
              WRITE(SYMLAB2,'(A,I2,"/2",";",I3,"/2")') SPDCAR(LL),J,MJ
            ELSEIF(LINEAR) THEN
              ID = IBEIG(II)
              WRITE(SYMLAB,'("m_j=",I3,"/2")') ID
            ENDIF
            IF(ATOMIC)THEN
              WRITE(LUPRI,'(/A,I4,A,'//FMT//',A,F6.4,A,2X,A11)')
     &          '* Positronic eigenvalue no.',II,': ',
     &          EIG(II),'       (Occupation : f = ',OCC2(I),')',
     &          SYMLAB2
            ELSE
              WRITE(LUPRI,'(/A,I4,A,'//FMT//',A,F6.4,A,2X,A9)')
     &          '* Positronic eigenvalue no.',II,': ',
     &          EIG(II),'       (Occupation : f = ',OCC2(I),')',
     &          SYMLAB
            ENDIF   
          END IF
        END IF
        CALL PRSYMB(LUPRI,'=',92,0)
       END IF ! (IPRPOP .GE. 0)
C
C       Calculate SO population matrix
C       ==============================
C
        CALL POPMAT(IFRP,I,PMAT,SMAT,CMO,NVEC,
     &              IPGRS,NPOPLAB,NCBAS,WORK(KFREE),LFREE,IPRPOP)
C
C       Calculate gross populations
C       ===========================
C
        CALL GRSPOP(POPG,PMAT,NPOPLAB)
        CALL WRPOPG(POPG,POPLAB,NPOPLAB,IEFF,IPRPOP)
C
C       Multiply fractional occupation on population.
C       =============================================
C
        IF(KRMC_FLG.EQ.0) THEN
           FACTOR  = D2*OCC2(I)
        ELSE IF(KRMC_FLG.EQ.2) THEN
C       ... MP2-NOs
           FACTOR = OCC2(I)
        ELSE IF(KRMC_FLG.EQ.3) THEN
C       ... CI-NOs
           FACTOR = OCC2(I)
        ELSE
           GO TO 1001
        END IF
        CALL DAXPY(NPOPLAB,FACTOR,POPG(1,1),1,PGTOT(1,IFRP),1)
        CALL DAXPY(NPOPLAB,FACTOR,POPG(1,2),1,PGTOT(1,IFRP),1)
C
C       Calculate net populations
C       ===========================
C
        IF(DONETP) THEN
          IF (IPRPOP.GT.0) THEN
            CALL WRPOPN(PMAT(1,1,1),POPG(1,1),POPLAB,NPOPLAB,
     &           IEFF,'alpha')
            CALL WRPOPN(PMAT(1,1,2),POPG(1,2),POPLAB,NPOPLAB,
     &           IEFF,'beta ')
          END IF
          CALL DAXPY(N2LAB,D2,PMAT(1,1,1),1,PNTOT(1,1,IFRP),1)
          CALL DAXPY(N2LAB,D2,PMAT(1,1,2),1,PNTOT(1,1,IFRP),1)
        ENDIF
 1001   CONTINUE
      ENDDO ! I = 1,NVEC
C
C     Write out total populations
C     ==========================
C
      IF(NFSYM.GT.1 .AND. IPRPOP .GE. 0) THEN
        WRITE(LUPRI,'(//A,A3,A)')
     &    '** Total gross population of fermion ircop ',
     &    FREP(IFRP),' **'
        CALL WRPOPT(PGTOT(1,IFRP),POPLAB,NPOPLAB,IEFF)
        IF(DONETP) THEN
          WRITE(LUPRI,'(//A,A3,A)')
     &      '** Total net population of fermion ircop ',
     &      FREP(IFRP),' **'
          CALL WRPOPN(PNTOT(1,1,IFRP),PGTOT(1,IFRP),
     &                POPLAB,NPOPLAB,IEFF,'total')
        ENDIF
      ENDIF

 9999 CONTINUE
      CALL MEMREL('MULPOP_2.POPMAT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MULPOP_2')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck popmat */
      SUBROUTINE POPMAT(IFRP,IVEC,PMAT,SMAT,CMO,NVECS,IGLAB,NGLAB,
     &                  NCBAS,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Construct population matrix
C         M(i,j) = Re(CONJ[CMO(i,IVEC)]*SMAT(i,j)*CMO(j,IVEC))
C     The final population matrix is compressed based on a labeling system.
C
C     Written by T.Saue Oct 9 1995
C     Last revision: Oct 9 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C Used from COMMON blocks
C  DCBAS: N2BBASX
#include "dcbbas.h"
      DIMENSION PMAT(*),SMAT(*),CMO(*),IGLAB(*),WORK(LWORK)
C
      CALL QENTER('POPMAT')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('REAL',KPBUF,2*N2BBASX,WORK,KFREE,LFREE)
C
      CALL POPMAT_1(IFRP,IVEC,PMAT,SMAT,CMO,NVECS,WORK(KPBUF),
     &     IGLAB,NGLAB,NCBAS,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('POPMAT.POPMAT_1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('POPMAT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck POPMAT_1 */
      SUBROUTINE POPMAT_1(IFRP,IVEC,PMAT,SMAT,CMO,NVECS,PBUF,
     &                  IGLAB,NGLAB,NCBAS,IPRINT)
C***********************************************************************
C
C     Construct population matrix
C         M(i,j) = Re(CONJ[CMO(i,IVEC)]*SMAT(i,j)*CMO(j,IVEC))
C     The final population matrix is compressed based on a labeling system.
C
C     Written by T.Saue Oct 9 1995
C     Last revision: Jan 7 1998 - jth LINUX
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
#include "maxorb.h"
C
C Used from COMMON blocks
C  DCBBAS: NTBAS,NFBAS
C  DCBORB: NORB
C  DGROUP: NZ
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbpop.h"
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION SMAT(NBAST,NBAST),PBUF(NBAST,NBAST,*)
#else
      DIMENSION SMAT(NTBAS(0),NTBAS(0)),PBUF(NTBAS(0),NTBAS(0),*)
#endif
      DIMENSION PMAT(NGLAB,NGLAB,2),CMO(NCBAS,NVECS,4),IGLAB(*)
C
C     Generate population matrix for basis
C     ====================================
C
      SUMA = D0
      SUMB = D0
      CALL DZERO(PBUF,N2BBASX*2)
      IF(ILABDF.EQ.1) THEN
        DO J = 1,NTBAS(0)
          DO I = 1,NTBAS(0)
            PBUF(I,J,1) = (CMO(I,IVEC,1)*CMO(J,IVEC,1)+
     &                     CMO(I,IVEC,2)*CMO(J,IVEC,2))*SMAT(I,J)
            SUMA = SUMA + PBUF(I,J,1)
            PBUF(I,J,2) = (CMO(I,IVEC,3)*CMO(J,IVEC,3)+
     &                     CMO(I,IVEC,4)*CMO(J,IVEC,4))*SMAT(I,J)
            SUMB = SUMB + PBUF(I,J,2)
          ENDDO
        ENDDO
        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,'(3X,A)') 'POPMAT: AO-basis'
          WRITE(LUPRI,'(6X,A,F9.4)') 'alpha-part:',SUMA
          IF(IPRINT.GE.4) THEN
            CALL OUTPUT(PBUF(1,1,1),1,NTBAS(0),
     &                1,NTBAS(0),NTBAS(0),NTBAS(0),1,LUPRI)
          ENDIF
          WRITE(LUPRI,'(6X,A,F9.4)') 'beta-part :',SUMB
          IF(IPRINT.GE.4) THEN
            CALL OUTPUT(PBUF(1,1,2),1,NTBAS(0),
     &                  1,NTBAS(0),NTBAS(0),NTBAS(0),1,LUPRI)
          ENDIF
        ENDIF
      ELSE
        IOFF = IBAS(IFRP)
        DO J = 1,NFBAS(IFRP,0)
          JJ = IOFF + J
          DO I = 1,NFBAS(IFRP,0)
            II = IOFF + I
            PBUF(II,JJ,1) = (CMO(I,IVEC,1)*CMO(J,IVEC,1)+
     &                       CMO(I,IVEC,2)*CMO(J,IVEC,2))*SMAT(II,JJ)
            SUMA = SUMA + PBUF(II,JJ,1)
            PBUF(II,JJ,2) = (CMO(I,IVEC,3)*CMO(J,IVEC,3)+
     &                       CMO(I,IVEC,4)*CMO(J,IVEC,4))*SMAT(II,JJ)
            SUMB = SUMB + PBUF(II,JJ,2)
          ENDDO
        ENDDO
        IF(IPRINT.GE.2) THEN
          IOFF = IBAS(IFRP) + 1
          WRITE(LUPRI,'(3X,A)') 'POPMAT: SO-basis'
          WRITE(LUPRI,'(6X,A,F9.4)') 'alpha-part:',SUMA
          IF(IPRINT.GE.4) THEN
            CALL OUTPUT(PBUF(IOFF,IOFF,1),1,NFBAS(IFRP,0),
     &                1,NFBAS(IFRP,0),NTBAS(0),NTBAS(0),1,LUPRI)
          ENDIF
          WRITE(LUPRI,'(6X,A,F9.4)') 'beta-part :',SUMB
          IF(IPRINT.GE.4) THEN
            CALL OUTPUT(PBUF(IOFF,IOFF,2),1,NFBAS(IFRP,0),
     &                  1,NFBAS(IFRP,0),NTBAS(0),NTBAS(0),1,LUPRI)
          ENDIF
        ENDIF
      ENDIF
C
C     Compress population matrix for labels
C     =====================================
C
      N2PMAT = NGLAB*NGLAB*2
      CALL DZERO(PMAT,N2PMAT)
      CALL GATMAT(1,NTBAS(0),NTBAS(0),PBUF(1,1,1),PMAT(1,1,1),
     &            IGLAB,NGLAB)
      CALL GATMAT(1,NTBAS(0),NTBAS(0),PBUF(1,1,2),PMAT(1,1,2),
     &            IGLAB,NGLAB)
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,'(3X,A,I5)') 'POPMAT: labels',NGLAB
        SUMA = D0
        DO J = 1,NGLAB
          DO I = 1,NGLAB
            SUMA = SUMA + PMAT(I,J,1)
          ENDDO
        ENDDO
        WRITE(LUPRI,'(6X,A,F9.4)') 'Sum alpha-part:',SUMA
        IF(IPRINT.GE.3) THEN
          CALL OUTPUT(PMAT(1,1,1),1,NGLAB,
     &                1,NGLAB,NGLAB,NGLAB,1,LUPRI)
        ENDIF
        SUMB = D0
        DO J = 1,NGLAB
          DO I = 1,NGLAB
            SUMB = SUMB + PMAT(I,J,2)
          ENDDO
        ENDDO
        WRITE(LUPRI,'(6X,A,F9.4)') 'Sum beta-part :',SUMB
        IF(IPRINT.GE.3) THEN
          CALL OUTPUT(PMAT(1,1,2),1,NGLAB,
     &                1,NGLAB,NGLAB,NGLAB,1,LUPRI)
        ENDIF
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grspop */
      SUBROUTINE GRSPOP(POPG,PMAT,NLAB)
C***********************************************************************
C
C     From Mulliken population matrix generate gross populations for
C     2-spinor.
C
C     Written by T.Saue - July 24 1994
C     LAST VERSION: Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION POPG(NLAB,2),PMAT(NLAB,NLAB,2)
C
      DO 10 J = 1,NLAB
        POPG(J,1) = DSUM(NLAB,PMAT(1,J,1),1)
        POPG(J,2) = DSUM(NLAB,PMAT(1,J,2),1)
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrpopg */
      SUBROUTINE WRPOPG(POPG,LAB,NLAB,IEFF,IPRINT)
C***********************************************************************
C
C     Output routine for Mulliken gross population
C
C     Written by T.Saue - Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DTOL=0.0001D0,D0 = 0.0D0)
      CHARACTER LAB(NLAB)*12,FMT*7
      DIMENSION POPG(NLAB,2),IEFF(*)
C
C     Find nonzero contributions
C     ==========================
C
      SUMA = D0
      SUMB = D0
      NEFF = 0
      DO I = 1,NLAB
        SUMA = SUMA + POPG(I,1)
        SUMB = SUMB + POPG(I,2)
        POPT = ABS(POPG(I,1))+ABS(POPG(I,2))
        IF(POPT.GE.DTOL) THEN
          NEFF = NEFF + 1
          IEFF(NEFF) = I
        ENDIF
      ENDDO
C
C     Write out nonzero contributions
C     ===============================
C
      IF (IPRINT .LT. 0) GO TO 9999
C
      NROW = (NEFF+1)/8
      IF(MOD((NEFF+1),8).GT.0) NROW = NROW + 1
      NPRI = MIN(NEFF+1,8)
      NFMT = 8 + 15*NPRI
      NPRI = NPRI - 1
      WRITE(LUPRI,'(/A,F8.5)') '* Gross populations greater than',DTOL
      WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
      WRITE(LUPRI,'(/A8,2X,A6,2X,A1,4X,7(A12,3X))')
     &     'Gross   ','Total ','|',(LAB(IEFF(I)),I=1,NPRI)
      WRITE(LUPRI,FMT) ('-',I=1,NFMT)
      WRITE(LUPRI,'(A8,F8.4,2X,A1,4X,7(F8.4,7X))')
     &     ' alpha  ',SUMA,'|',(POPG(IEFF(I),1),I=1,NPRI)
      WRITE(LUPRI,'(A8,F8.4,2X,A1,4X,7(F8.4,7X))')
     &     ' beta   ',SUMB,'|',(POPG(IEFF(I),2),I=1,NPRI)
      NOFF = NPRI
      DO J = 2,NROW
        NPRI = MIN(NEFF-NOFF,8)
        NFMT = 8 + 15*NPRI
        WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
        WRITE(LUPRI,'(/A8,8(1X,A12,2X))')
     &       'Gross  |',(LAB(IEFF(NOFF+I)),I=1,NPRI)
        WRITE(LUPRI,FMT) ('-',I=1,NFMT)
        WRITE(LUPRI,'(A8,8(1X,F8.4,6X))')
     &       ' alpha |',(POPG(IEFF(NOFF+I),1),I=1,NPRI)
        WRITE(LUPRI,'(A8,8(1X,F8.4,6X))')
     &       ' beta  |',(POPG(IEFF(NOFF+I),2),I=1,NPRI)
        NOFF = NOFF + NPRI
      ENDDO
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/A)') '* Gross populations (all):'
        DO I = 1,NLAB
          WRITE(LUPRI,'(A12,1P,2E16.8)') LAB(I),POPG(I,1),POPG(I,2)
        ENDDO
      ENDIF
C
 9999 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck netpop */
      SUBROUTINE NETPOP(POPN,PMAT,NLAB)
C***********************************************************************
C
C     Mulliken net populations (upper triangular columnwise packed)
C
C     Written  by T.Saue Oct 11 1995
C     Last revision: Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D2 = 2.0D0)
      DIMENSION POPN((NLAB*(NLAB+1)/2),2),PMAT(NLAB,NLAB,2)
C
C     Generate net populations
C
      IND = 0
      DO J = 1,NLAB
        IND = IND + 1
        POPN(IND,1) = PMAT(J,J,1)
        POPN(IND,2) = PMAT(J,J,2)
        DO I = 1,J
          IND = IND + 1
          POPN(IND,1) = D2*PMAT(I,J,1)
          POPN(IND,2) = D2*PMAT(I,J,2)
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrpopn */
      SUBROUTINE WRPOPN(PMAT,POPG,LAB,NLAB,IEFF,COMP)
C***********************************************************************
C
C     Output routine for Mulliken net populations
C
C     Written by T.Saue Oct 11 1995
C     Last revision Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DTOL=0.0001D0,D2 = 2.0D0)
      CHARACTER LAB(NLAB)*12,FMT*7,PFMT*23,COMP*5
      DIMENSION PMAT(NLAB,NLAB),POPG(NLAB),IEFF(NLAB)
C
C     Find nonzero contributions
C     ==========================
C
      NEFF = 0
      DO I = 1,NLAB
        POPT = ABS(POPG(I))
        IF(POPT.GE.DTOL) THEN
          NEFF = NEFF + 1
          IEFF(NEFF) = I
        ENDIF
      ENDDO
      IF(NEFF.EQ.0) RETURN
C
      CALL HEADER('Net population -- '//COMP,0)
C
C     Write net populations
C     =====================
C
      NBATCH = NEFF/8
      IF(MOD(NEFF,8).GT.0) NBATCH = NBATCH + 1
      NOFF = 0
      DO J = 1,NBATCH
        NEND = MIN(NEFF-NOFF,8)
        NFMT = 8 + 15*NEND
        WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
        WRITE(LUPRI,'(/A8,8(A12,3X))')
     &       'Net.pop.',(LAB(IEFF(NOFF+K)),K=1,NEND)
        WRITE(LUPRI,FMT) ('-',I=1,NFMT)
        NT = NOFF+NEND
        DO I = 1,NT
          NBEG = MAX(NOFF+1,I)
          NPRI = NT - NBEG + 1
          NX   = NEND-NPRI
          IF(NX.GT.0) THEN
            WRITE(PFMT,'(A,I4,A,I2,A)')
     &         '(A12,',15*NX,'X,',NPRI,'(F8.4,7X))'
          ELSE
            WRITE(PFMT,'(A,I2,A)') '(A12,',NPRI,'(F8.4,7X))      '
          ENDIF
          IF(I.EQ.NBEG) THEN
            WRITE(LUPRI,PFMT) LAB(IEFF(I)),
     &         PMAT(IEFF(I),IEFF(I)),
     &         (D2*PMAT(IEFF(NOFF+K),IEFF(I)),K=NX+2,NEND)
          ELSE
            WRITE(LUPRI,PFMT) LAB(IEFF(I)),
     &       (D2*PMAT(IEFF(NOFF+K),IEFF(I)),K=NX+1,NEND)
          ENDIF
        ENDDO
        NOFF = NOFF + NEND
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck anaprj */
      SUBROUTINE ANAPRJ(WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
C
#include "dcbprj.h"
#include "dgroup.h"
      DIMENSION NVECS(2,0:2),KVEC(2,0:MAXREF),NSTR(2,0:2,0:MAXREF)
      DIMENSION WORK(LWORK)
#include "memint.h"
      CALL TITLER('Projection analysis','*',122)
C
C     Prepare selected molecular and fragment orbitals
C     ================================================
C
      NDIM=2*3*(MAXREF+1)
      CALL IZERO(NSTR,NDIM)
      IOPT=14
      IF(PATOMS) THEN
C.....Use atomic fragments in C1 basis
        CALL PRJAOB(KCSEL,KESEL,KBSEL,KBVEC,NVECS,KVEC,NSTR,IOPT,
     &              KRMC_FLG,WORK,KFREE,LFREE)
      ELSE
C.....Use symmetry-adapted fragments
        CALL PRJSOB(KCSEL,KESEL,KBSEL,KBVEC,NVECS,KVEC,NSTR,IOPT,
     &              KRMC_FLG,WORK,KFREE,LFREE)
      ENDIF
C
C     Do projection analysis
C     ======================
C
C     The polarization contribution will be appended to the solution vector
C
      CALL ANAPR1(NVECS,KVEC,NSTR,WORK(KCSEL),WORK(KESEL),WORK(KBSEL),
     &            WORK(KBVEC),IOPT,KRMC_FLG,WORK,KFREE,LFREE)
C
C     If requested, do phase transformation of molecule coefficients
C     using a selected reference orbital for alignment
C     ==============================================================
C
C      IF(PRJADJ) THEN
C        CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRJANA')
C        CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
C        CALL MEMGET('REAL',KEIG,NORBT   ,WORK,KFREE,LFREE)
C        CALL MEMGET('INTE',KIBE,NORBT   ,WORK,KFREE,LFREE)
C        CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),WORK(KIBE),
C     &              TOTERG,14)
C
C        DO IFRP = 1,NFSYM
C          CALL PRJPHA(IREFADJ(IFRP),IFRP,
C     &                WORK(KCMO+ICMOQ(IFRP)),NFBAS(IFRP,0),NORB(IFRP),
C     &                WORK(KBV(IFRP)),NVECS(IFRP,0),NVECS(IFRP,2),
C     &                WORK(KVEC(IFRP,0)),NSTR(IFRP,2,0),NSTR(IFRP,1,0))
C        ENDDO
C        REWIND LUCOEF
C        CALL WRICMO(LUCOEF,WORK(KCMO),WORK(KEIG),WORK(KIBE),TOTERG,.FALSE.)
C        CLOSE(LUCOEF,STATUS='KEEP')
C      ENDIF
C
      CALL MEMREL('ANAPRJ',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck anasob */
      SUBROUTINE PRJSOB(KCSEL,KESEL,KBSEL,KBVEC,NVECS,KVEC,NSTR,IOPT,
     &                  KRMC_FLG,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0,DM1=-1.0D0,D1=1.0D0)
C
#include "dcbgen.h"
#include "dcbprj.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dummy.h"
C
      LOGICAL TOBE,FNDLAB
      DIMENSION WORK(*)
      DIMENSION NSTR(2,0:2,0:MAXREF),KVEC(2,0:MAXREF),
     &          NVECS(2,0:2),IOFF(2),IDIM(3,2)
      CHARACTER TEXT*74
      KFRSAV=KFREE
C
C
C     If requested, do Mayer polarization of reference orbitals
C     =========================================================
C
CTROND      IF(MAYER) THEN
CTROND        CALL MAYER_POLARIZE(WORK,KFREE,LFREE)
CTROND      ENDIF
C
C     Find total numbers of reference orbitals;
C     first scan coefficient files for NPSH/NESH;
C     this is needed in the case when 'all' or 'oo'
C     is used in the orbital strings
C
      DO I = 1,NREFS
        CALL OPNFIL(LUCOEF,REFFIL(I),'OLD','ANAPRJ')
        IF(FNDLAB('INFO    ',LUCOEF)) THEN         
          READ (LUCOEF,END=10,ERR=20) TEXT,NSYM,NZBUF,
     &       (NSTR(J,2,I),NSTR(J,1,I),IDUM,J=1,NSYM),DUM
C.............NPSH........NESH........NBAS...........TOTERG
       ELSE
          WRITE(6,*) 'PRJSOB: Error reading first record on ',REFFIL(I)     
          CALL QUIT('PRJSOB: Error reading info')
        ENDIF
        CLOSE(LUCOEF,STATUS='KEEP')
      ENDDO
      DO IFRP = 1,NFSYM
        NVECS(IFRP,1) = 0
        DO I = 1,NREFS
          CALL ORBNUM(VECREF(IFRP,I),IFRP,KVEC(IFRP,I),NSTR(1,0,I),
     &                WORK,KFREE,LFREE)
          NVECS(IFRP,1) = NVECS(IFRP,1) + NSTR(IFRP,0,I)
        ENDDO
      ENDDO
C
C     Find number of molecular orbitals to analyze
C
      DO IFRP = 1,NFSYM
        NSTR(IFRP,1,0) = NESH(IFRP)
        NSTR(IFRP,2,0) = NPSH(IFRP)
        CALL ORBNUM(VECPRJ(IFRP),IFRP,KVEC(IFRP,0),NSTR(1,0,0),
     &              WORK,KFREE,LFREE)
        NVECS(IFRP,2) = NSTR(IFRP,0,0)
      ENDDO
C
C     Selected eigenvalues/syminfo
C
      NDIM = 0
      DO IFRP = 1, NFSYM
        NVECS(IFRP,0) = NVECS(IFRP,1) + NVECS(IFRP,2)
        NDIM = NDIM + NVECS(IFRP,0)
      ENDDO
      CALL MEMGET('REAL',KESEL,NDIM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KBSEL,NDIM,WORK,KFREE,LFREE)
C     Calculate dimensions of coefficient arrays:
C       reference orbitals are followed by molecular ones
C
      NCDIM = 0
      DO IFRP = 1,NFSYM
        NCDIM = NCDIM + NFBAS(IFRP,0)*NVECS(IFRP,0)*NZ
      ENDDO
      CALL MEMGET('REAL',KCSEL,NCDIM,WORK,KFREE,LFREE)
C
C     COEFFICIENT SELECTION
C     =====================
C
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBE,NORBT   ,WORK,KFREE,LFREE)
      CALL IZERO(WORK(KIBE),NORBT)
C
C
C     Select reference orbitals and eigenvalues
C
      DO IFRP = 1,NFSYM
        IOFF(IFRP) = 0
      ENDDO
      CALL SELFRAG(WORK(KCSEL),NVECS(1,0),WORK(KESEL),WORK(KBSEL),
     &             NREFS,NSTR(1,0,1),REFFIL,OWNBAS,NPROJNUC,IOPT,
     &             KVEC(1,1),WORK(KCMO),WORK(KEIG),WORK(KIBE),
     &             IOFF,KDUM,WORK,KFREE,LFREE)
C
C     Get molecular coefficients and eigenvalues
C     Note that integer arrays IOFF is used to get correct
C     positioning...
C
      CALL SELFRAG(WORK(KCSEL),NVECS(1,0),WORK(KESEL),WORK(KBSEL),
     &             1,NSTR(1,0,0),'DFCOEF',.FALSE.,1,IOPT,
     &             KVEC(1,0),WORK(KCMO),WORK(KEIG),WORK(KIBE),
     &             IOFF,KRMC_FLG,WORK,KFREE,LFREE)
C
C     Throw away the full coefficients and eigenvalues
C     ================================================
C
      CALL MEMREL('PRJSOB',WORK,KFRSAV,KCMO,KFREE,LFREE)
      NBDIM = 0
      DO IFRP = 1,NFSYM
        NBDIM = NBDIM + NVECS(IFRP,0)*NVECS(IFRP,2)*NZ
      ENDDO
      CALL MEMGET('REAL',KBVEC,NBDIM,WORK,KFREE,LFREE)
      IF(IPRPRJ.GE.2) THEN
        NDIM = 0
        DO IFRP = 1,NFSYM
          NDIM = NDIM + NVECS(IFRP,0)
        ENDDO
        WRITE(6,*) 'All selected eigenvalues'
        CALL OUTPUT(WORK(KESEL),1,NDIM,1,1,NDIM,1,-1,LUPRI)
      ENDIF
      RETURN
 10   CONTINUE
      CALL QUIT('PRJSOB: END OF FILE reading TEXT')
 20   CONTINUE
      CALL QUIT('PRJSOB: ERROR reading TEXT')
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjan1 */
      SUBROUTINE PRJAN1(CMO,EIG,IBEIG,OCC,BVEC,NTOT,NREF,NMOL,NBAS,
     &                  NFRAGV,FREP,IPQ,IFRP,NZP,JMOL,JREF,NSTR,SMAT,
     &                  GROSS,POL,IOPT,KRMC_FLG,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Projection analysis:
C
C     Molecular orbitals |Y(mol;i> are expanded as
C
C     |Y(mol;i> = SUM(B,n)|Y(ref;B,n)> c(B,n,i) + |Y(pol,i)>
C
C     where |Y(ref;B,n)> are fragment orbitals of center B with index n.
C     The orthogonal complement is the polarization contribution |Y(pol,i)>
C
C     The expansion coefficients c(B,n,i) are found from the linear equation
C
C     SUM(B,n) <Y(ref;A,m)|Y(ref;B,n> c(B,n,i) = <Y(ref;A,m)|Y(mol;i)>
C
C     The first NREF columns of CMO contain the coefficients of the 
C     fragment orbitals. On input the next NMOL columns contain the 
C     coefficients of the molecular orbitals to be analyzed.
C     On ouput these are replaced by the polarization contributions.
C
C     KRMC_FLG = 1: Orbitals read from
C                   label NEWORB on file KRMCSCF
C     KRMC_FLG = 0  Orbitals read from DFCOEF
C
C    Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,DC=1.0D-2,D2=2.0D0,D1=1.0D0,DM1=-1.0D0)
C
#include "dcbprj.h"
C
      LOGICAL DOPOP
      CHARACTER*3 FREP(2)
      DIMENSION GROSS(2,0:2,MAXREF)
      DIMENSION CMO(NBAS,NTOT,NZP),EIG(NTOT),IBEIG(NTOT),OCC(NMOL),
     &          BVEC(NTOT,NMOL,NZP),SMAT(NBAS,NBAS),NFRAGV(NREFS)
      DIMENSION NSTR(2,0:2,0:MAXREF),JREF(3,NREF),JMOL(3,NMOL),IPQ(NZP),
     &          WORK(LFREE)
C
      KFRSAV = KFREE
      IMOL   = NREF+1
      DOPOP  = NREFS.GT.1
      PTOL   = PROTHR*DC
      IF(NREFS.GT.1.OR.POLREF) THEN
        NDIM=NTOT*NTOT*NZP
        CALL MEMGET('REAL',KAMAT,NTOT*NTOT*NZP,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KAMAT1,NTOT*NTOT,WORK,KFREE,LFREE)
      ELSE
        KAMAT  = KFREE
        KAMAT1 = KFREE
      ENDIF
C
C     If requested, polarize reference using scheme of Knizia
C
      IF(POLREF) THEN
        CALL PRJAN2(CMO(1,IMOL,1),EIG(IMOL),NMOL,NBAS,NTOT,
     &              CMO          ,EIG      ,NREF,NBAS,NTOT,
     &              NTOT,NBAS,NZP,IPQ,WORK(KAMAT),WORK(KAMAT1),
     &              BVEC,SMAT,PTOL,
     &              NREFS,REFFIL,NFRAGV,IPRINT,WORK,KFREE,LFREE)
        CALL KNIZIA(CMO,BVEC,SMAT,WORK(KAMAT),NTOT,NMOL,
     &              NBAS,NZP,IPQ,NREFS,REFFIL,NFRAGV,
     &              IPRPRJ,WORK,KFREE,LFREE)
      ENDIF
C
C     Get fragment expansion coefficients
C
      CALL PRJAN2(CMO(1,IMOL,1),EIG(IMOL),NMOL,NBAS,NTOT,
     &            CMO          ,EIG      ,NREF,NBAS,NTOT,
     &            NTOT,NBAS,NZP,IPQ,WORK(KAMAT),WORK(KAMAT1),
     &            BVEC,SMAT,PTOL,
     &            NREFS,REFFIL,NFRAGV,IPRINT,WORK,KFREE,LFREE)
C
C     Replace molecular coefficients by polarization contribution 
C
      CALL DUNIT2(BVEC(IMOL,1,1),NMOL,NTOT,NMOL,NZP)
      CALL QGEMM(NBAS,NMOL,NREF,DM1,
     &       'N','N',IPQ,CMO ,NBAS,NTOT,NZP,
     &       'N','N',IPQ,BVEC,NTOT,NMOL,NZP,
     &       D1,IPQ,CMO(1,IMOL,1),NBAS,NTOT,NZP)
      IF(IPRINT.GE.2) THEN
        CALL HEADER('Solution vectors',0)
        CALL PRQMAT(BVEC,NTOT,NMOL,NTOT,NMOL,NZP,IPQ,LUPRI)        
      ENDIF
      IF(IPRINT.GE.5) THEN
        CALL HEADER('Polarization vectors',0)
        CALL PRIPOL(CMO(1,IMOL,1),EIG(IMOL),NBAS,NMOL,NZP,NBAS,NTOT)
      ENDIF
C
C     Make overlap matrix AMAT = <\psi^B_j1\psi^A_i>
C     but now with polarization contributions included
C
      IF(NREFS.GT.1) THEN
        CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NTOT,NTOT,
     &           SMAT       ,NBAS,NBAS, 1,IPQ,
     &           WORK(KAMAT),NTOT,NTOT,NZP,IPQ,
     &           CMO,NBAS,NTOT,NZP,IPQ,
     &           CMO,NBAS,NTOT,NZP,IPQ,
     &           WORK(KFREE),LFREE,IPRINT)
      ENDIF
C
C     Mulliken-type population analysis
C
      !Miro-fix out-of-bounds, array size increased 
      CALL MEMGET('REAL',KGVEC,(NREF+2)*2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT,NREF*NREF ,WORK,KFREE,LFREE)
      CALL ANAPR2(DOPOP,NREF,NMOL,NTOT,WORK(KPMAT),WORK(KGVEC),
     &            WORK(KAMAT),BVEC,NZP,
     &            JMOL,EIG(IMOL),IBEIG(IMOL),NSTR,
     &            JREF,EIG,GROSS,POL,OCC,IOPT,KRMC_FLG,IPQ)
      CALL MEMREL('PRJAN1',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRJAN2(CMOL,EMOL,NMOL,LRMOL,LCMOL,
     &                  CREF,EREF,NREF,LRREF,LCREF,
     &                  NTOT,NBAS,NZP,IPQ,AMAT,AMAT1,BVEC,SMAT,PTOL,
     &                  NFRAG,FRAGLAB,NFRAGV,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get fragment expansion coefficients
C     Written by Trond Saue
C
C***********************************************************************
#include "priunit.h"
C
      real(8),parameter       :: D0=0.0D0,DM1=-1.0D0,D1=1.0D0
      character(6),intent(in) :: FRAGLAB(NFRAG)
      integer,intent(in)      :: LRMOL,LCMOL,LRREF,LCREF
      integer,intent(in)      :: NREF,NMOL,NTOT,NBAS,NZP,NFRAG
      integer,intent(in)      :: NFRAGV(NFRAG),IPQ(NZP),IPRINT
      real(8),intent(in)      :: CMOL(LRMOL,LCMOL,NZP),EMOL(NMOL)
      real(8),intent(in)      :: CREF(LRREF,LCREF,NZP),EREF(NREF)
      real(8),intent(in)      :: SMAT(NBAS,NBAS),PTOL
      integer,intent(inout)   :: KFREE,LFREE
      real(8),intent(out)     :: AMAT(NTOT,NTOT,NZP),BVEC(NTOT,NMOL,NZP)
      real(8),intent(out)     :: AMAT1(NTOT,NTOT,1)
      real(8),intent(out)     :: WORK(*)
      real(8), allocatable    :: BUF1(:),BUF2(:)
      integer, allocatable    :: IPIVOT(:)
      real(8)                 :: AMAT_SUM_NZP(1:4),AMAT_SUM,AMAT_SUM_RE

C
C
C     Make projection vector  BVEC = <\psi^B_j1\psi^{MO}_i>
C
      CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NREF,NMOL,
     &            SMAT,NBAS ,NBAS ,1 ,IPQ,
     &            BVEC,NTOT ,NMOL ,NZP,IPQ,
     &            CREF,LRREF,LCREF,NZP,IPQ,
     &            CMOL,LRMOL,LCMOL,NZP,IPQ,
     &            WORK(KFREE),LFREE,IPRINT)
      IF(IPRINT.GE.2) THEN
        WRITE(6,*) '* Projection vector BVEC: '
        CALL PRQMAT(BVEC,NREF,NMOL,NTOT,NMOL,NZP,IPQ,LUPRI)
      ENDIF
      IF(NFRAG.GT.1) THEN
C
C       Make overlap matrix AMAT = <\psi^B_j1\psi^A_i>
C
        CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NREF,NREF,
     &             SMAT,NBAS ,NBAS ,1 ,IPQ,
     &             AMAT,NTOT ,NTOT ,NZP,IPQ,
     &             CREF,LRREF,LCREF,NZP,IPQ,
     &             CREF,LRREF,LCREF,NZP,IPQ,
     &             WORK(KFREE),LFREE,IPRINT)

!Miro: get real numbers matrix AMAT1 of quaternion matrix AMAT
        CALL NZ_MTX_AMPLITUDES(AMAT,AMAT1,NTOT,NZP)
C
C       Print overlaps
C       ==============
C
        IF(IPRINT.GE.1) THEN
          JOFF = 1
          DO J = 1,NFRAG
            IOFF = 1
            DO I = 1,(J-1)
              WRITE(LUPRI,'(/,A,2(3X,A6))') '*Overlaps :',
     &             FRAGLAB(I),FRAGLAB(J)
              CALL PRQMAT(AMAT(IOFF,JOFF,1),NFRAGV(I),NFRAGV(J),
     &                    NTOT,NTOT,NZP,IPQ,LUPRI)
              WRITE(LUPRI,'(A,2(3X,A6))') '*Overlaps - amplitudes :',
     &             FRAGLAB(I),FRAGLAB(J)
              CALL PRQMAT(AMAT1(IOFF,JOFF,1),NFRAGV(I),NFRAGV(J),
     &                    NTOT,NTOT,1,IPQ,LUPRI)

!Miro: sum all overlaps for given fragments
              CALL SUM_NZ_MATRIX(
     &             AMAT(IOFF,JOFF,1),NFRAGV(I),NFRAGV(J),
     &             NTOT,NTOT,NZP,AMAT_SUM_NZP,AMAT_SUM)

              CALL SUM_RE_MATRIX(AMAT1(IOFF,JOFF,1),
     &             NFRAGV(I),NFRAGV(J),NTOT,NTOT,AMAT_SUM_RE)

              WRITE(LUPRI,'(A,2(3X,A6),//)')
     & '*Overlaps - summed per all projected orbitals, per fragments :',
     &                FRAGLAB(I),FRAGLAB(J)
              WRITE(LUPRI,'(2F12.4)') AMAT_SUM,AMAT_SUM_RE

              IOFF = IOFF + NFRAGV(I)
            ENDDO
            JOFF = JOFF + NFRAGV(J)
          ENDDO
        ENDIF
C
C       Solve linear system by Cholesky decomposition
C       with full pivoting
C
        JOB = 1
        allocate(BUF1(NREF))
        allocate(BUF2(NREF))
        allocate(IPIVOT(NREF))
        CALL QCHOLD(AMAT,NREF,NZP,NTOT,NTOT,BUF1,PTOL,NEFF,JOB,IPIVOT)
        CALL QCHOLS(AMAT,NREF,NEFF,NMOL,NZP,NTOT,NTOT,BUF1,BVEC,
     &         NTOT,NMOL,BVEC,NTOT,NMOL,JOB,IPIVOT,BUF2)
        deallocate(BUF1)
        deallocate(BUF2)
        deallocate(IPIVOT)
        IF(NEFF.LT.NREF) THEN
          WRITE(LUPRI,'(A,A)') '* WARNING: ',
     &     ' linear dependencies detected in Cholesky decomposition'
          WRITE(LUPRI,'(11X,A,I5,A,I5)')
     &     ' Reduced dimensionality: ',NREF,' --> ',NEFF
        ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjinf */
      SUBROUTINE PRJINF(JREF,IND,IFRAG,IVEC,NSTR,KREF,NPOFF)
C***********************************************************************
C
C     Make pointer vector of reference orbital information for use
C     in projection analysis.
C
C     Wrintten by T.Saue March 31 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbprj.h"
      DIMENSION JREF(3,*),NSTR(2,0:2,0:MAXREF),KREF(*)
      IFLAG = IFRAG
      IF(IFRAG.EQ.0) IFLAG = IND
C
C     Positrons
C
      IP = 2
      DO J = 1,NSTR(IND,2,IFRAG)
        JREF(1,IVEC+J) = NPOFF+1+KREF(J)
        JREF(2,IVEC+J) = IP
        JREF(3,IVEC+J) = IFLAG
      ENDDO
C
C     Electrons
C
      IP = 1
      DO J = NSTR(IND,2,IFRAG)+1,NSTR(IND,0,IFRAG)
        JREF(1,IVEC+J) = KREF(J)
        JREF(2,IVEC+J) = IP
        JREF(3,IVEC+J) = IFLAG
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck anapr2 */
      SUBROUTINE ANAPR2(DOPOP,NREF,NMOL,NTOT,PMAT,GVEC,AMAT,BVEC,NZP,
     &                  JMOL,EMOL,IBMOL,NSTR,JREF,EREF,GROSS,POL,OCC,
     &                  IOPT,KRMC_FLG,IPQ)
C***********************************************************************
C
C     Using the coefficients obtained from projection analysis,
C     carry out a Mulliken-like analysis.
C     PROTHR is threshold for print of projection coeffiencients.
C
C     Written by Trond Saue March 13 1996
C          revision Jan  7 1998 - jth NFBAS(IFRP,0) is now parameter
C     Last revision Apr 16 2010 - SK introduced proper handling of 
C                                    MCSCF-NO coefficients, occupation
C                                    numbers (and mj-values)
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dcbprj.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbham.h"      
C     
      LOGICAL DOPOP
      CHARACTER FMT*6,MXFORM*6,OTYP*1,SYMLAB*9
      DIMENSION GROSS(2,0:2,MAXREF),IPQ(*)
      DIMENSION PMAT(NREF,NREF),GVEC(NREF+2,2),AMAT(NTOT,NTOT,NZP),
     &          BVEC(NTOT,NMOL,NZP),JMOL(3,*),EMOL(*),IBMOL(*),
     &          NSTR(2,0:2,0:MAXREF),JREF(3,*),EREF(*),OCC(*)
      DIMENSION P(4),OTYP(2)
      REAL*8 Dno_occ
      DATA SYMLAB/'         '/
      DATA OTYP/'E','P'/
      SAVE OTYP
C
C     Initialization
C
      CALL DZERO(GVEC,(NREF+2))
      CALL DCOPY(4,D0,0,P,1)
C
C     Loop over molecular orbitals
C
      Dno_occ = D2
      IF(KRMC_FLG.EQ.1) Dno_occ = D1
C
      IFRP = JMOL(3,1)
      CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
      DO 10 K = 1,NMOL
        FMT = MXFORM(EMOL(K),16)
        IP = JMOL(2,K)
        IFRP_TMP = JMOL(3,K)
        IF(IFRP_TMP.NE.IFRP) THEN
          IFRP = IFRP_TMP
          CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
        ENDIF
        IF(IP.EQ.1) THEN
C...    Positive-energy (electronic) orbital
          IF(KRMC_FLG.EQ.0)THEN
            IF(SPINFR) THEN
              ID = IBMOL(K)
              WRITE(SYMLAB,'("sym= ",A3,1X)') REP(ID)
            ELSEIF(LINEAR) THEN
              ID = IBMOL(K)
              WRITE(SYMLAB,'("m_j=",I3,"/2")') ID
            ENDIF
            WRITE(LUPRI,'(/A,I4,A,'//FMT//',A,F6.4,A,2X,A9)')
     &         '* Electronic eigenvalue nr.',JMOL(1,K),': ',
     &         EMOL(K),'       (Occupation : f = ',OCC(K),')',
     &         SYMLAB
         ELSE ! MCSCF
            ID = IBMOL(K)
            WRITE(SYMLAB,'("m_j=",I3,"/2")') ID
            IF( ID .eq. 0) THEN
              WRITE(LUPRI,'(/A,I4,A,A,F12.8)')
     &         '* Electronic eigenvalue nr.',JMOL(1,K),': ',
     &            '      Occupation :',OCC(K)
            ELSE
                  WRITE(LUPRI,'(/A,I4,A,A,F12.8,2X,A9)')
     &            '* Electronic eigenvalue no.',JMOL(1,K),': ',
     &            '      Occupation :',OCC(K),SYMLAB
            END IF
          END IF
        ELSE
C...    Negative-energy (positronic) orbital
          IF(SPINFR) THEN
            ID = IBMOL(K)
            WRITE(SYMLAB,'("sym= ",A3,1X)') REP(ID)
          ELSEIF(LINEAR) THEN
            ID = IBMOL(K)
            WRITE(SYMLAB,'("m_j=",I3,"/2")') ID
          ENDIF
          WRITE(LUPRI,'(/A,I4,A,'//FMT//',A,F6.4,A,2X,A9)')
     &       '* Positronic eigenvalue nr.',JMOL(1,K),': ',
     &       EMOL(K),'       (Occupation : f = ',OCC(K),')',
     &       SYMLAB
        ENDIF
        CALL PRSYMB(LUPRI,'=',92,0)
        WRITE(LUPRI,100) 'Orbital','Total','Eigenvalue',
     &       'Kramers partner 1','Kramers partner 2'
C
C       Print coefficients
C
        NPRI = 0
C
        DO 20 I = 1,NREF
          CW = D0
          DO IZ = 1,NZP
            CW = CW + BVEC(I,K,IZ)*BVEC(I,K,IZ)
          ENDDO
          CW = SQRT(CW)
          IF(CW.LT.PROTHR) GOTO 20
          NPRI = NPRI + 1
          KORB = JREF(1,I)
          KTYP = JREF(2,I)
          KREF = JREF(3,I)
          DO IZ = 1,NZP
            P(IZ) = BVEC(I,K,IZ)/CW
          ENDDO
C
          WRITE(LUPRI,110) REFFIL(KREF),OTYP(KTYP),KORB,CW,EREF(I),
     &     '(',P(1),',',P(2),')','(',P(3),',',P(4),')'
   20   CONTINUE
        IF(NPRI.EQ.0) THEN
          WRITE(LUPRI,'(A)') ' --> No contributions !'
          GOTO 10
        ENDIF

C
C       Population analysis
C
        IF(DOPOP) THEN
          N2REF = NREF*NREF
          CALL DZERO(PMAT,N2REF)
C
C         Calculate net population matrix
C
          DO JZ = 1,NZP
            DO IZ = 1,NZP
              KZ  = IQMULT(IZ,JZ,1)
              FAC = IQSIGN(IZ,2,1)*IQSIGN(JZ,1,1)*IQPHASE(IZ,JZ,KZ)
              DO J = 1,NREF
                DO I = 1,NREF
                  PMAT(I,J) = PMAT(I,J) + FAC*
     &                        BVEC(I,K,IZ)*AMAT(I,J,KZ)*BVEC(J,K,JZ)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C
C         Calculate gross population vector
C
          IF(LWGPOP) THEN
            DO J = 1,NREF
              GVEC(J,2) = PMAT(J,J)
              DO I = 1,(J-1)
                CW = PMAT(I,I) + PMAT(J,J)
                IF(CW.GT.D0) THEN
                 GVEC(J,2) = GVEC(J,2) + D2*(PMAT(J,J)/CW)*PMAT(I,J)
                 GVEC(I,2) = GVEC(I,2) + D2*(PMAT(I,I)/CW)*PMAT(I,J)
                ENDIF
              ENDDO
            ENDDO
          ELSE
            DO J = 1,NREF
              GVEC(J,2) = DSUM(NREF,PMAT(1,J),1)
            ENDDO
          ENDIF
         IF(IPRPRJ.GE.2) THEN         
            WRITE(LUPRI,'(A)') '* Net population'
            CALL OUTPUT(PMAT,1,NREF,1,NREF,NREF,NREF,1,LUPRI)
            WRITE(LUPRI,'(A)') '* Gross population'
            CALL OUTPUT(GVEC(1,2),1,NREF,1,1,NREF,1,1,LUPRI)
          ENDIF
        ELSE
C
C         Calculate gross population vector
C
          DO J = 1,NREF
            GVEC(J,2) = BVEC(J,K,1)*BVEC(J,K,1)
            DO JZ = 2,NZP
              GVEC(J,2) = GVEC(J,2) + BVEC(J,K,JZ)*BVEC(J,K,JZ)
            ENDDO
          ENDDO
          IF(IPRPRJ.GE.2) THEN          
            WRITE(LUPRI,'(A)') '* Gross population'
            CALL OUTPUT(GVEC(1,2),1,NREF,1,1,NREF,1,1,LUPRI)
          ENDIF
        ENDIF
        WRITE(LUPRI,'(A)') '* Gross contributions:'
        CALL GROSSC(GROSS,NREFS,IFRP,NREF,D1,GVEC(1,2),JREF,REFFIL)
        FACTOR  = OCC(K)*Dno_occ
        CALL DAXPY((NREF+1),FACTOR,GVEC(1,2),1,GVEC(1,1),1)
 10   CONTINUE
      CALL PRSYMB(LUPRI,'=',48,0)
      WRITE(LUPRI,'(/A)')
     &    '* Total reference orbital contributions:'
      DO I = 1,NREF
        KORB = JREF(1,I)
        KTYP = JREF(2,I)
        KREF = JREF(3,I)
        WRITE(LUPRI,110)
     &   REFFIL(KREF),OTYP(KTYP),KORB,GVEC(I,1),EREF(I)
      ENDDO
      WRITE(LUPRI,'(/A)') '* Total gross contributions:'
      TOTM = Dno_occ*DSUM(NMOL,OCC,1)
      CALL GROSSC(GROSS,NREFS,IFRP,NREF,TOTM,GVEC(1,1),JREF,REFFIL)
      POL = GVEC(NMOL+1,1)
C
      RETURN
 100  FORMAT(4X,A7,12X,A5,11X,A10,4X,2(3X,A17))
 110  FORMAT(4X,A6,2X,A1,I6,3X,F11.5,5X,E15.8,
     &       2(3X,A1,F7.4,A1,F7.4,A1))
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getplb */
      SUBROUTINE GETPLB(IPRINT)
C***********************************************************************
C
C     Generate labels for use in population analysis
C
C     Written by Trond Saue May 1996
C     Last revision: tsaue - May 2 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "dcbpop.h"
C
#include "mxcent.h"
#include "nuclei.h"
C
#include "dcbibt.h"
C
      IF (DO_SCFPOP) THEN
C- combine all large and small component types for each center
C- omit point charges
        NPOPLAB = 0
        DO I = 1, NPLAB(0)
          ICENT      = JGET(IATTR(I,ILABDF))
          IDEG       = KGET(IATTR(I,ILABDF))
          ICENT_DEP  = NUCPRE(ICENT) + IDEG
          IPOPLAB(I)         = ICENT_DEP
          CHRGNUC(ICENT_DEP) = CHARGE(ICENT)
          POPLAB(ICENT_DEP)  = '    '//NAMDEP(ICENT_DEP)
          IF (.NOT. NOORBT(ICENT)) THEN   ! omit point charges
            NPOPLAB = MAX(NPOPLAB, ICENT_DEP)
          END IF
        END DO
        GO TO 8000
      ELSE IF (ADD_ALL) THEN
C- combine all large and small component types for each center
        DO I = 1, NPLAB(0)
          ICENT = JGET(IATTR(I,ILABDF))
          IDEG  = KGET(IATTR(I,ILABDF))
          ICENT_DEP = NUCPRE(ICENT) + IDEG
          IPOPLAB(I) = ICENT_DEP
        END DO
        NPOPLAB = NUCDEP
        DO I = 1,NPOPLAB
           POPLAB(I) = '    '//NAMDEP(I)
        END DO
        GO TO 8000
      END IF
C- large component types
      IF(.NOT.LABDEF) THEN
        DO I = 1,NPLAB(1)
          POPLAB(I) = PLABEL(I,ILABDF)
          IPOPLAB(I)  = I
        ENDDO
        NPOPLAB = NPLAB(1)
      ENDIF
      IF(ADDSML) THEN
C- combine all small component types for each center
        JCENT = 0
        JDEG  = 0
        DO I = NPLAB(1)+1,NPLAB(0)
          ICENT = JGET(IATTR(I,ILABDF))
          IDEG  = KGET(IATTR(I,ILABDF))
          IF(ICENT.NE.JCENT.OR.IDEG.NE.JDEG) THEN
            NPOPLAB = NPOPLAB + 1
            POPLAB(NPOPLAB) = PLABEL(I,ILABDF)(3:8)//'_small'
            JCENT = ICENT
            JDEG  = IDEG
          ENDIF
          IPOPLAB(I) = NPOPLAB
        ENDDO
      ELSE
C- individual small component types
        IF(.NOT.LABDEF) THEN
          DO I = NPLAB(1)+1,NPLAB(0)
            POPLAB(I) = PLABEL(I,ILABDF)
            IPOPLAB(I)  = I
          ENDDO
          NPOPLAB = NPOPLAB + NPLAB(2)
        ENDIF
      ENDIF
C
C     Print section
C     =============
C
 8000 IF(IPRINT.GE.2) THEN
        CALL HEADER('Output from GETPLB: Group labels',-1)
        WRITE(LUPRI,'(A,I5)') '* Total number : ',NPOPLAB
        WRITE(LUPRI,'(6(I6,2X,A12))') (I,POPLAB(I),I=1,NPOPLAB)
        IF(IPRINT.GE.5) THEN
          CALL PRIVEC('IPOPLAB ',IPOPLAB,NPLAB(0))
        ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck anavec */
      SUBROUTINE ANAVEC(WORK,LWORK)
C***********************************************************************
C
C     Print vectors
C
C     Written by Trond Saue May 1996
C     Last revision: tsaue - Jan 9 1997
C     Vector analysis for KRMC added: T. Fleig, Oct 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "dcbgen.h"
#include "dcbvec.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "consts.h"
C
      LOGICAL TOBE
      DIMENSION WORK(*)
      DIMENSION IOFPFS(2,2)
#include "memint.h"
C
      NTESTL = 000
C
      IF(IPRCMP.LE.0) RETURN
      CALL TITLER('Vector print','*',130)
C
C     Check if MCSCF coefficients are on file
C
      INQUIRE (FILE ='KRMCSCF',EXIST=TOBE)
C
      IF(.NOT.TOBE) GOTO 990
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
C
      CALL OPNFIL(LUKRMC,'KRMCSCF','OLD','PAMANA')
      IOPT = 1       ! Read from NEWORB label
      JRDMO = -1
      CALL RREADMO(WORK(KCMO),JRDMO,IOPT,LUKRMC)
      CLOSE(LUKRMC,STATUS='KEEP')
C
      IF (NTESTL.GE.20) THEN
C
C     Print all coefficients (including positronic and S component)
C
        DO I = 1, NFSYM
           WRITE(LUPRI,'(A,I2)') ' (ROPTST)  Coefficients, irrep ',I
           CALL PRQMAT(WORK(KCMO+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &          NFBAS(I,0),NORB(I),NZP,IPQTOQ(1,0),LUPRI)
        END DO
      END IF
C
      CALL HEADER('Coefficients from KRMCSCF',2)
      KEIG = KFREE
      IPREIG = 0
      CALL PRCOEF(WORK(KCMO),WORK(KEIG),VECPRI,
     &            IPREIG,IPRCMP,ILABDF,IPRVEC,WORK(KFREE),LFREE)
C     Memory deallocation
      CALL MEMREL('ANAVEC',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
C
C     No MCSCF, then check if RHF on DFCOEF :
C
 990  CONTINUE
      INQUIRE (FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KBEIG,NORBT   ,WORK,KFREE,LFREE)
C
C     Read coefficients
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMANA')
      IOPT=14
      IF(ATOMIC) THEN
        CALL ATOMIC_REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),
     &                     WORK(KBEIG),TOTERG,IOPT)           
      ELSE
        CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),WORK(KBEIG),
     &              TOTERG,IOPT)
      ENDIF
      CLOSE(LUCOEF,STATUS='KEEP')
C
      IF (NTESTL.GE.10) THEN
C
C     Check orbital energies.
C
        DO J=1,NORBT,1
          WRITE(6,*) 'ORBENE(',J,') = ',WORK(KEIG+J-1)
        END DO
      END IF
C
      CALL HEADER('Coefficients from DFCOEF',2)
      IPREIG = 1
      CALL PRCOEF(WORK(KCMO),WORK(KEIG),VECPRI,
     &            IPREIG,IPRCMP,ILABDF,IPRVEC,WORK(KFREE),LFREE)
C     Memory deallocation
      CALL MEMREL('ANAVEC',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(//A)') 'ANAVEC: Coefficient file not found !'
      CALL QUIT('ANAVEC: Coefficient file not found !')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrpopt */
      SUBROUTINE WRPOPT(POPG,LAB,NLAB,IEFF)
C***********************************************************************
C
C     Output routine for Mulliken gross population
C
C     Written by T.Saue - Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DTOL=0.0001D0,D0 = 0.0D0)
      CHARACTER LAB(NLAB)*12,FMT*7
      DIMENSION POPG(NLAB),IEFF(*)
C
C     Find nonzero contributions
C     ==========================
C
      SUMA = D0
      NEFF = 0
      DO I = 1,NLAB
        SUMA = SUMA + POPG(I)
        POPT = ABS(POPG(I))
        IF(POPT.GE.DTOL) THEN
          NEFF = NEFF + 1
          IEFF(NEFF) = I
        ENDIF
      ENDDO
C
C     Write out nonzero contributions
C     ===============================
C
      NROW = (NEFF+1)/8
      IF(MOD((NEFF+1),8).GT.0) NROW = NROW + 1
      NPRI = MIN(NEFF+1,8)
      NFMT = 8 + 15*NPRI
      NPRI = NPRI - 1
      WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
      WRITE(LUPRI,'(/A8,3X,A6,3X,A1,4X,7(A12,3X))')
     &     'Gross   ','Total ','|',(LAB(IEFF(I)),I=1,NPRI)
      WRITE(LUPRI,FMT) ('-',I=1,NFMT)
      WRITE(LUPRI,'(A8,F10.5,2X,A1,4X,7(F9.5,6X))')
     &     ' total  ',SUMA,'|',(POPG(IEFF(I)),I=1,NPRI)
      NOFF = NPRI
      DO J = 2,NROW
        NPRI = MIN(NEFF-NOFF,8)
        NFMT = 8 + 15*NPRI
        WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
        WRITE(LUPRI,'(/A8,8(1X,A12,2X))')
     &       'Gross  |',(LAB(IEFF(NOFF+I)),I=1,NPRI)
        WRITE(LUPRI,FMT) ('-',I=1,NFMT)
        WRITE(LUPRI,'(A8,8(1X,F9.5,5X))')
     &       ' total |',(POPG(IEFF(NOFF+I)),I=1,NPRI)
        NOFF = NOFF + NPRI
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck WR_MULCHRG */
      SUBROUTINE WR_MULCHRG(POPG,LAB,NLAB)
C***********************************************************************
C
C     Output routine for Mulliken atomic charges
C
C     Written by H. J. Aa. Jensen, Jan 2010, based on WRPOPT
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER LAB(NLAB)*12,FMT*10
      DIMENSION POPG(NLAB)
      PARAMETER (NCOL = 12, LCOL=10)
C
C     ===============================
C
      WRITE(LUPRI,'(A)') '* Mulliken charges:'
      NROW = 1 + (NLAB-1)/NCOL
      NOFF = 0
      DO J = 1,NROW
        NPRI = MIN(NLAB-NOFF,NCOL)
        NFMT = LCOL*NPRI - 2
        WRITE(FMT,'(A4,I3,A3)') '(2X,',NFMT,'A1)'
        IF (J .EQ. 1) WRITE(LUPRI,FMT) ('-',I=1,NFMT)
        WRITE(LUPRI,'(20A10)') (LAB(NOFF+I)(1:10),I=1,NPRI)
        WRITE(LUPRI,'(20F10.3)') (POPG(NOFF+I),I=1,NPRI)
        WRITE(LUPRI,FMT) ('-',I=1,NFMT)
        NOFF = NOFF + NPRI
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck grossc */
      SUBROUTINE GROSSC(GROSS,NFRAG,IFRP,NREF,TOTM,GVEC,JREF,REFFIL)
C***********************************************************************
C
C     Determine gross contribution in projection analysis
C
C     Written by T.Saue Oct 5 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
      CHARACTER REFFIL*6
      DIMENSION GROSS(2,0:2,NFRAG)
      DIMENSION GVEC(*),REFFIL(*),JREF(3,NREF)
      CALL DZERO(GROSS,6*NFRAG)
      GTOT = D0
      DO I = 1,NREF
        IP    = JREF(2,I)
        IFRAG = JREF(3,I)
        GROSS(IFRP,IP,IFRAG) = GROSS(IFRP,IP,IFRAG) + GVEC(I)
        GROSS(IFRP,0 ,IFRAG) = GROSS(IFRP,0 ,IFRAG) + GVEC(I)
        GTOT  = GTOT + GVEC(I)
      ENDDO
      DO I = 1,NFRAG
        WRITE(LUPRI,'(3X,A6,8X,3X,F10.4,2(3X,A4,F10.4))')
     &      REFFIL(I),GROSS(IFRP,0,I),
     &      'E - ',GROSS(IFRP,1,I),'P - ',GROSS(IFRP,2,I)
      ENDDO
      GVEC(NREF+1) = TOTM - GTOT
      WRITE(LUPRI,'(3X,A,F10.4)') 'Polarization: ',GVEC(NREF+1)
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gt1rho  */
      SUBROUTINE GT1RHO(WORK,LWORK)
C***********************************************************************
C
C     Get density in one dimension, that is along bonds
C     (ghost centers may be defined to get other lines)
C
C     Written by Trond Saue May 9 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
      PARAMETER (D2=2.0D0)
      LOGICAL TOBE
      DIMENSION WORK(LWORK)
#include "dcbgen.h"
#include "dcbrho1.h"
#include "dcbbas.h"
C
C for NFMAT:
#include "dcbdhf.h"
C
      CALL QENTER('GT1RHO')
#include "memint.h"
C
      CALL TITLER('Density along bonds','*',122)
      IF (NFMAT .GT. 1) THEN
        WRITE(LUPRI,'(/A)')
     &  '--> Not implemented for open shell yet, nothing done.'
        GO TO 9999
      END IF
C
C     Check if coefficients are on file
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) THEN
        WRITE(LUPRI,'(/A)')'GT1RHO: Coefficient file DFCOEF not found !'
        CALL QUIT('GT1RHO: Coefficients not found !')
      END IF
C     Memory allocation
CTROND: Consider open shell ! (NFMAT...)
      CALL MEMGET('REAL',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAOS,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOOR,6       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Read coefficients + give info
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','GT1RHO')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),TEMP,TEMP,
     &            TOTERG,3)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Generate density matrix
C     Be aware: Both density matrices are normed to 1
C     Coeffients may now be deallocated.....
C
      CALL DENMAT(WORK(KDMAT),WORK(KCMO),IPRHO1)
      IF(IOPTDD.NE.0)
     &  CALL ATTACHMENT_DENSITY(WORK(KDMAT),IOPTDD,IPRHO1,
     &                          WORK,KFREE,LFREE)
      CALL DSCAL(N2BBASXQ,D2,WORK(KDMAT),1)
      CALL GT1RH1(WORK(KDMAT),WORK(KGAOS),WORK(KCOOR),WORK(KBUF),
     &            IPRHO1,WORK,KFREE,LFREE)
C
C     Memory deallocation
      CALL MEMREL('GT1RHO',WORK,KWORK,KWORK,KFREE,LFREE)
C
 9999 CALL QEXIT('GT1RHO')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gt1rh1 */
      SUBROUTINE GT1RH1(DMAT,GAO,COORD,BUF,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get density in one dimension, that is along bonds
C     (ghost centers may be defined to get other lines)
C
C     Written by Trond Saue May 9 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "nuclei.h"
      DIMENSION DMAT(*),GAO(*),COORD(3,2),BUF(*),WORK(*)
      CHARACTER FILNAM*10
#include "ibtfun.h"
C
      LURHO = 1
C
C     Loop over pairs of nuclear centers
C
      DO ICENTA = 1, NUCIND
         MULA  = ISTBNU(ICENTA)
         IDEGA = 1
         COORD(1,1) = CORD(1,ICENTA)
         COORD(2,1) = CORD(2,ICENTA)
         COORD(3,1) = CORD(3,ICENTA)
         DO ICENTB = ICENTA+1, NUCIND
           MULB = ISTBNU(ICENTB)
           IDEGB = 0
           DO ISYMOP = 0, MAXOPR
             IF (IBTAND(ISYMOP,MULB) .EQ. 0) THEN
               IDEGB = IDEGB + 1
               COORD(1,2)=PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENTB)
               COORD(2,2)=PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENTB)
               COORD(3,2)=PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENTB)
               WRITE(FILNAM,'(A2,2(A3,I1))')
     &           'RH',NAMN(ICENTA),IDEGA,NAMN(ICENTB),IDEGB
               DO I = 1,10
                 IF(FILNAM(I:I).EQ.' ') FILNAM(I:I)='_'
               ENDDO
               OPEN(LURHO,FILE =FILNAM,STATUS='NEW',
     &              ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
               CALL GT1RH2(DMAT,GAO,COORD,BUF,LURHO,
     &             IPRINT,WORK,KFREE,LFREE)
               CLOSE(LURHO,STATUS='KEEP')
             ENDIF
           ENDDO
         ENDDO
         IF (MULT(MULA).GT.1) THEN
           IDEGB = 1
           DO ISYMOP = 1, MAXOPR
             IF (IBTAND(ISYMOP,MULA) .EQ. 0) THEN
               IDEGB = IDEGB + 1
               COORD(1,2)=PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENTA)
               COORD(2,2)=PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENTA)
               COORD(3,2)=PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENTA)
               WRITE(FILNAM,'(A2,2(A3,I1))')
     &           'RH',NAMN(ICENTA),IDEGA,NAMN(ICENTA),IDEGB
               DO I = 1,10
                 IF(FILNAM(I:I).EQ.' ') FILNAM(I:I)='_'
               ENDDO
               OPEN(LURHO,FILE =FILNAM,STATUS='NEW',
     &              ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
               CALL GT1RH2(DMAT,GAO,COORD,BUF,LURHO,
     &                     IPRINT,WORK,KFREE,LFREE)
               CLOSE(LURHO,STATUS='KEEP')
              ENDIF
           ENDDO
         ENDIF
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Gt1rh2 */
      SUBROUTINE GT1RH2(DMAT,GAO,COORD,BUF,LURHO,
     &                  IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get density in one dimension, that is along bonds
C     (ghost centers may be defined to get other lines)
C
C     Written by Trond Saue May 9 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
C
      DIMENSION DMAT(*),GAO(*),COORD(3,2),RHO(0:2),BUF(*),WORK(*)
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbrho1.h"
#include "codata.h"
C
C     Calculate distance between centers
C
      DX = (COORD(1,2)-COORD(1,1))
      DY = (COORD(2,2)-COORD(2,1))
      DZ = (COORD(3,2)-COORD(3,1))
      DISTAN = SQRT(DX*DX+DY*DY+DZ*DZ)
      DX     = DX/DISTAN
      DY     = DY/DISTAN
      DZ     = DZ/DISTAN
      ASTEP  = DSTEP/XTANG
      NPOINT = INT(DISTAN/ASTEP)
      IF(MOD(NPOINT,2).EQ.1) NPOINT = NPOINT+1
      ASTEP  = DISTAN/dble(NPOINT)
C
C     Calculate density along between centers:
C     If bond distance is 2R, then points are calculated between (-R,3R)
C
      IP1 = -NPOINT/2
      IP2 = NPOINT - IP1
      DO IP = IP1,IP2
        PP  = IP*ASTEP
        PX  = COORD(1,1) + PP*DX
        PY  = COORD(2,1) + PP*DY
        PZ  = COORD(3,1) + PP*DZ
C       Get value of SOs at (PX,PY,PZ)
        CALL GETSOS(GAO,TEMP,TEMP,TEMP,TEMP,PX,PY,PZ,
     &              BUF,NTBAS(0),0,.FALSE.,IPRINT)
C       Get density in this point
        RHO(0) = D0
        DO IC = 1,MC
          RHO(IC) = D0
          CALL GETRH1(RHO(IC),IC,0,GAO,DMAT,BUF)
          RHO(0) = RHO(0) + RHO(IC)
        ENDDO
        PP = XTANG*PP
        WRITE(LURHO,'(F10.6,3E20.10)') PP,(RHO(I),I=0,MC)
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gt3rho  */
      SUBROUTINE GT3RHO(WORK,LWORK)
C***********************************************************************
C
C     Prepare for 3-dimensional density plot
C
C     Written by Trond Saue Nov 5 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0)
      LOGICAL TOBE,LBIT
      DIMENSION WORK(LWORK),CUBE(3,0:4),LUCUBE(2)
#include "dcbgen.h"
#include "dcbrho.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
C
      CALL QENTER('GT3RHO')
c
c     Start timer
c
      TIMSTR = SECOND()
#include "memint.h"
C
      CALL TITLER('3-dimensional density plot','*',122)
C
      IF(LEVYLE) ILSCUB = 1
      IF(LBIT(ILSCUB,1)) THEN
         WRITE(LUPRI,'(A)')
     &   ' GT3RHO : Gaussian cube file generated for large component.'
         LUCUBE(1) = 11
         OPEN(LUCUBE(1),FILE='rhoL.cube',STATUS='NEW',
     &        ACCESS='SEQUENTIAL',FORM='FORMATTED')
      ENDIF
      IF (LBIT(ILSCUB,2)) THEN
         WRITE(LUPRI,'(A)')
     &   ' GT3RHO : Gaussian cube file generated for small component.'
         LUCUBE(2) = 12
         OPEN(LUCUBE(2),FILE='rhoS.cube',STATUS='NEW',
     &        ACCESS='SEQUENTIAL',FORM='FORMATTED')
      ENDIF
      CALL MEMGET('REAL',KATOM,4*NUCDEP,WORK,KFREE,LFREE)
C     Get cube parameters
      CALL CUBPAR(CUBE,WORK(KATOM))
C     Write header for cube file
      DO IC = 1,2
        IF(LBIT(ILSCUB,IC)) CALL CUBHED(LUCUBE(IC),CUBE,WORK(KATOM))
      ENDDO
C
C     Check if coefficients are on file
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
C     Memory allocation
CTROND: Consider open shell ! (NFMAT...)
      CALL MEMGET('REAL',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAOS ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF  ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Read coefficients + give restart info
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),TEMP,TEMP,
     &            TOTERG,3)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Generate density matrix
C     Be aware: Both density matrices are normed to 1
      CALL GENDEN(WORK(KDMAT),WORK(KCMO),1,IPRRHO)
C     Coeffients may now be deallocated.....
      CALL MEMREL('GT3RHO.cmo',WORK,KCMO,KCMO,KFREE,LFREE)
C
      CALL MEMGET('REAL',KRHO,NCUBE(3)*2,WORK,KFREE,LFREE)
      CALL GT3RH1(LUCUBE,CUBE,WORK(KRHO),
     &            WORK(KDMAT),WORK(KGAOS),WORK(KBUF),
     &            WORK,KFREE,LFREE)
C
      DO IC = 1,2
        IF(LBIT(ILSCUB,IC)) CLOSE(LUCUBE(IC),STATUS='KEEP')
      ENDDO
C     Memory deallocation
      CALL MEMREL('GT3RHO',WORK,KWORK,KWORK,KFREE,LFREE)
C
c
c     Stop timer
c
      TIMEND = SECOND()
      TIME   = TIMEND - TIMSTR
      CALL TIMTXT('>>> Time used in GT3RHO is',TIME,LUPRI)
      CALL QEXIT('GT3RHO')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'GT3RHO: Coefficient file not found !'
      CALL QUIT('GT1RHO: Coefficients not found !')
C
      END
      SUBROUTINE WRITE_ORB_CUBES(N1,N2,ORBS,WORK,LWORK)
c     Write cube files for gerade orbitals orbs(1,1)..orbs(N1,1)
c     and ungerade orbitals orbs(1,2)..orbs(N2,2).
c     orbital indices are defined so that 1 is first occupied.
c     Based on GT3RHO, this one by ulfek.
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      integer n1,n2,orbs
      PARAMETER(D0=0.0D0,D1=1.0D0)
      LOGICAL TOBE,LBIT
      DIMENSION WORK(LWORK),CUBE(3,0:4),LUCUBE(2),orbs(1000,2)
      character*16 filename
#include "dcbgen.h"
#include "dcbrho.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
C
      CALL QENTER('GT3RHO')
c
c     Start timer
c
      TIMSTR = SECOND()
#include "memint.h"
C
      CALL TITLER('3-dimensional density plot','*',122)
C
      IF(LEVYLE) ILSCUB = 1
      CALL MEMGET('REAL',KATOM,4*NUCDEP,WORK,KFREE,LFREE)
C
C     Check if coefficients are on file
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
C     Memory allocation
CTROND: Consider open shell ! (NFMAT...)
      CALL MEMGET('REAL',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAOS ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF  ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Read coefficients + give restart info
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),TEMP,TEMP,
     &            TOTERG,3)
      CLOSE(LUCOEF,STATUS='KEEP')

C     Get cube parameters
      CALL CUBPAR(CUBE,WORK(KATOM))

      CALL MEMGET('REAL',KRHO,NCUBE(3)*2,WORK,KFREE,LFREE)

      do i=1,N1+N2
      IF(LBIT(ILSCUB,1)) THEN
         if (i.le.n1) then
            write(filename,'(A,I4.4,A)') 'rhoL',orbs(i,1),'g.cube'
         else
            write(filename,'(A,I4.4,A)') 'rhoL',orbs(i-N1,2),'u.cube'
         endif
         WRITE(LUPRI,'(A)')
     &   ' GT3RHO : Gaussian cube file generated for large component.'
         LUCUBE(1) = 11
         OPEN(LUCUBE(1),FILE=filename,STATUS='NEW',
     &        ACCESS='SEQUENTIAL',FORM='FORMATTED')
      ENDIF
      IF (LBIT(ILSCUB,2)) THEN
         if (i.le.n1) then
            write(filename,'(A,I4.4,A)') 'rhoS',orbs(i,1),'g.cube'
         else
            write(filename,'(A,I4.4,A)') 'rhoS',orbs(i-N1,2),'u.cube'
         endif
         WRITE(LUPRI,'(A)')
     &   ' GT3RHO : Gaussian cube file generated for small component.'
         LUCUBE(2) = 12
         OPEN(LUCUBE(2),FILE=filename,STATUS='NEW',
     &        ACCESS='SEQUENTIAL',FORM='FORMATTED')
      ENDIF
C     Write header for cube file
      DO IC = 1,2
        IF(LBIT(ILSCUB,IC)) CALL CUBHED(LUCUBE(IC),CUBE,WORK(KATOM))
      ENDDO

C     Generate density matrix
      if (i.gt.n1) then
         ifrp = 2
         iorb = i - n1
      else
         ifrp = 1
         iorb = i
      endif
C     Be aware: Both density matrices are normed to 1
      call denorb(work(kdmat),orbs(iorb,ifrp),1,ifrp,work(kcmo),IPRRHO)
      CALL GT3RH1(LUCUBE,CUBE,WORK(KRHO),
     &            WORK(KDMAT),WORK(KGAOS),WORK(KBUF),
     &            WORK,KFREE,LFREE)
C
      DO IC = 1,2
        IF(LBIT(ILSCUB,IC)) CLOSE(LUCUBE(IC),STATUS='KEEP')
      ENDDO
      enddo
C     Coeffients may now be deallocated.....
      CALL MEMREL('GT3RHO.cmo',WORK,KCMO,KCMO,KFREE,LFREE)
C     Memory deallocation
      CALL MEMREL('GT3RHO',WORK,KWORK,KWORK,KFREE,LFREE)
C
c
c     Stop timer
c
      TIMEND = SECOND()
      TIME   = TIMEND - TIMSTR
      CALL TIMTXT('>>> Time used in GT3RHO is',TIME,LUPRI)
      CALL QEXIT('GT3RHO')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'GT3RHO: Coefficient file not found !'
      CALL QUIT('GT1RHO: Coefficients not found !')
C
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Gt3rh1 */
      SUBROUTINE GT3RH1(LUCUBE,CUBE,RHO,DMAT,GAO,BUF,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     3-dimensional plot of density
C
C     Written by Trond Saue Nov 5 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
      LOGICAL LBIT
      DIMENSION LUCUBE(*),CUBE(3,0:3),RHO(2,*),DMAT(*),GAO(*),
     &          BUF(*),WORK(*)
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbrho.h"
C
C
      DO I1 = 1,NCUBE(1)
        DO I2 = 1,NCUBE(2)
          DO I3 = 1,NCUBE(3)
            PX = CUBE(1,0)+(I1-1)*CUBE(1,1)
     &                    +(I2-1)*CUBE(1,2)
     &                    +(I3-1)*CUBE(1,3)
            PY = CUBE(2,0)+(I1-1)*CUBE(2,1)
     &                    +(I2-1)*CUBE(2,2)
     &                    +(I3-1)*CUBE(2,3)
            PZ = CUBE(3,0)+(I1-1)*CUBE(3,1)
     &                    +(I2-1)*CUBE(3,2)
     &                    +(I3-1)*CUBE(3,3)
C           Get value of SOs at (PX,PY,PZ)
            CALL GETSOS(GAO,TEMP,TEMP,TEMP,TEMP,PX,PY,PZ,
     &              BUF,NTBAS(0),0,.FALSE.,IPRRHO)
C           Get density in this point
            DO IC = 1,2
              IF(LBIT(ILSCUB,IC)) THEN
                RHO(IC,I3) = D0
                CALL GETRH1(RHO(IC,I3),IC,0,GAO,DMAT,BUF)
              ENDIF
            ENDDO
          ENDDO
C         Write buffer
          DO IC = 1,2
            IF(LBIT(ILSCUB,IC))
     &      WRITE(LUCUBE(IC),'(6E13.5)') (RHO(IC,I),I=1,NCUBE(3))
          ENDDO
        ENDDO
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Cubpar */
      SUBROUTINE CUBPAR(CUBE,ATOMS)
C***********************************************************************
C
C     Generate parameters for the cube. All coordinates in Angstroms.
C
C     NCUBE(0) - Number of atoms
C     NCUBE(1) - Number of points in the "X"-direction
C     NCUBE(2) - Number of points in the "Y"-direction
C     NCUBE(3) - Number of points in the "Z"-direction
C
C     CUBE(*,0) - Coordinates of initial point
C     CUBE(*,1) - Step vector for the "X"-direction
C     CUBE(*,2) - Step vector for the "Y"-direction
C     CUBE(*,3) - Step vector for the "Z"-direction
C
C     ATOMS(1,I) - x-coordinate of atom I
C     ATOMS(2,I) - y-coordinate of atom I
C     ATOMS(3,I) - z-coordinate of atom I
C     ATOMS(4,I) - Charge of atom I
C
C     Written by Trond Saue Nov 6 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "symmet.h"
#include "nuclei.h"
#include "dcbrho.h"
#include "codata.h"
      DIMENSION CUBE(3,0:3),ATOMS(4,NUCDEP)
C
#include "ibtfun.h"
C
      NCUBE(0) = NUCDEP
C     Initialize
      CALL DZERO(CUBE,12)
      DO I = 1,3
        TEMP = CORD(I,1)
        CUBE(I,0) = TEMP
        CUBE(I,I) = TEMP
      ENDDO
C     Loop over all centers
      JCENT = 0
      DO ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
C        Loop over all centers
         DO ISYMOP = 0, MAXOPR
         IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
           JCENT = JCENT + 1
           DO I = 1,3
             ATOMS(I,JCENT) =
     &             PT(IBTAND(ISYMAX(I,1),ISYMOP))*CORD(I,ICENT)
             CUBE(I,0) = MIN(CUBE(I,0),ATOMS(I,JCENT))
             CUBE(I,I) = MAX(CUBE(I,I),ATOMS(I,JCENT))
           ENDDO
           ATOMS(4,JCENT)=CHARGE(ICENT)
         ENDIF
         ENDDO
      ENDDO
C     Add some space, calculate step sizes and convert step vectors to Angstrom
      DO I = 1,3
        CUBE(I,0) = CUBE(I,0)  - CUBADJ(1)
        CUBE(I,I) = (CUBE(I,I) + CUBADJ(2))*XTANG
        CUBE(I,I) = (CUBE(I,I)-CUBE(I,0))/(NCUBE(I)-1)
      ENDDO
C      CALL DSCAL(12,XTANG,CUBE,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Cubhed */
      SUBROUTINE CUBHED(LUCUBE,CUBE,ATOMS)
C***********************************************************************
C
C     Write header for formatted Gaussian cube file
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbgen.h"
#include "nuclei.h"
#include "dcbrho.h"
C
      DIMENSION CUBE(3,0:3),ATOMS(4,NUCDEP)
      WRITE(LUCUBE,10) TITLE
      WRITE(LUCUBE,10) 'Total density'
      DO I = 0,3
        WRITE(LUCUBE,20) NCUBE(I),(CUBE(J,I),J=1,3)
      ENDDO
      DO I = 1,NUCDEP
        WRITE(LUCUBE,20) NINT(ATOMS(4,I)),ATOMS(4,I),
     &                   (ATOMS(J,I),J=1,3)
      ENDDO
C
      RETURN
   10 FORMAT(1X,A)
   20 FORMAT(I5,4F12.6)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Anapr3 */
      SUBROUTINE ANAPR3(AMAT,BVEC,NMOL,NREF,NTOT,EREF,JREF,
     &                  PSUM,WORK,KFREE,LFREE)
C***********************************************************************
C
C      Do Loewdin symmetric orthonormalization
C      Written by T. Saue Oct 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dgroup.h"
#include "dcbprj.h"
      CHARACTER OTYP(2)*1
      DIMENSION AMAT(NTOT,NTOT,NZ),BVEC(NTOT,NMOL,NZ),
     &          JREF(3,*),EREF(*),PSUM(NREFS),WORK(*)
      DATA OTYP/'E','P'/
      SAVE OTYP
C
      KFRSAV = KFREE
C
C     Get transformation matrix
C
      CALL SQSTRA(1,AMAT,NREF,NTOT,NTOT,NZ,WORK(KFREE),LFREE)
C
C     Transform coefficients
C
      NDIM = NTOT*NTOT*NZ
      CALL MEMGET('REAL',KCVEC,NDIM,WORK,KFREE,LFREE)
      CALL QGEMM(NREF,NMOL,NREF,D1,
     &           'N','N',IPQTOQ(1,0),AMAT,NTOT,NTOT,NZ,
     &           'N','N',IPQTOQ(1,0),BVEC,NTOT,NMOL,NZ,
     &           D0,IPQTOQ(1,0),WORK(KCVEC),NTOT,NMOL,NZ)
C
C     Form density matrix in symmetry orthonormalized basis
C
      CALL QGEMM(NTOT,NTOT,NMOL,D1,
     &           'N','N',IPQTOQ(1,0),WORK(KCVEC),NTOT,NMOL,NZ,
     &           'H','N',IPQTOQ(1,0),WORK(KCVEC),NTOT,NMOL,NZ,
     &           D0,IPQTOQ(1,0),AMAT,NTOT,NTOT,NZ)
      IF(IPRPRJ.GE.2) THEN
        WRITE(LUPRI,'(/A/)')
     &   'Coefficients in symmetry orthonormalized basis:'
        CALL PRQMAT(WORK(KCVEC),NTOT,NMOL,NTOT,NMOL,NZ,
     &              IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/A/)')
     &  'Density matrix in symmetry orthonormalized basis:'
        CALL PRQMAT(AMAT,NTOT,NTOT,NTOT,NTOT,NZ,
     &            IPQTOQ(1,0),LUPRI)
      ENDIF
      WRITE(LUPRI,'(/A)')
     &   '* Occupations in symmetry orthonormalized basis:'
      CALL DZERO(PSUM,NREFS)
      DO I = 1,NREF
        KORB = JREF(1,I)
        KTYP = JREF(2,I)
        KREF = JREF(3,I)
        WRITE(LUPRI,110)
     &   REFFIL(KREF),OTYP(KTYP),KORB,D2*AMAT(I,I,1),EREF(I)
        PSUM(KREF) = PSUM(KREF) + D2*AMAT(I,I,1)
      ENDDO
      WRITE(LUPRI,'(/A)')
     &   '* Fragment occupations in symmetry orthonormalized basis:'
      DO I = 1,NREFS
        WRITE(LUPRI,'(3X,A6,8X,3X,F7.4)')
     &      REFFIL(I),PSUM(I)
      ENDDO
      CALL MEMREL('ANAPR3',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      RETURN
 110  FORMAT(4X,A6,2X,A1,I6,3X,F11.5,5X,E15.8,
     &       2(3X,A1,F5.2,A1,F5.2,A1,3X))
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjpha */
      SUBROUTINE PRJPHA(ISEL,IFRP,CMO,NBAS,NVEC,BVEC,NTOT,NMOL,
     &                  JVEC,NPVEC,NEVEC)
C***********************************************************************
C
C     Phase adjustment using a selected reference orbital
C
C     Written by T. Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION CMO(NBAS,NVEC,NZ),BVEC(NTOT,NMOL,NZ),JVEC(NMOL)
C
C     Positronic vectors
C
      DO I = 1,NPVEC
        II = NPSH(IFRP)+1+JVEC(I)
        CALL PRJPH1(ISEL,CMO,II,NBAS,NVEC,BVEC,I,NTOT,NMOL)
      ENDDO
C
C     Electronic vectors
C
      ITOP = NEVEC+NPVEC
      DO I = NPVEC+1,ITOP
        II = NPSH(IFRP)+JVEC(I)
        CALL PRJPH1(ISEL,CMO,II,NBAS,NVEC,BVEC,(I+NPVEC),NTOT,NMOL)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjph1 */
      SUBROUTINE PRJPH1(ISEL,CMO,IC,NBAS,NVEC,BVEC,IB,NTOT,NMOL)
C***********************************************************************
C
C     Phase adjustment using a selected reference orbital
C
C     Written by T. Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0, DM1 = -1.0D0)
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION CMO(NBAS,NVEC,NZ),BVEC(NTOT,NMOL,NZ),JVEC(NMOL)
      DIMENSION A(4),B(4)
      IF(NZ.EQ.1) THEN
        IF(BVEC(ISEL,IB,1).LT.D0) CALL DSCAL(NBAS,DM1,CMO(1,IC,1),1)
      ELSEIF(NZ.EQ.2) THEN
        DB = PYTHAG(BVEC(ISEL,IB,1),BVEC(ISEL,IB,2))
        IF(DB.GT.D0) THEN
          DO IZ = 1,NZ
            A(IZ) =  BVEC(ISEL,IB,IZ)/DB
          ENDDO
          DO I = 1,NBAS
            B(1)        = CMO(I,IC,2)*A(1) - CMO(I,IC,1)*A(2)
            CMO(I,IC,1) = CMO(I,IC,1)*A(1) + CMO(I,IC,2)*A(2)
            CMO(I,IC,2) = B(1)
          ENDDO
        ENDIF
      ELSEIF(NZ.EQ.4) THEN
        DB = PYTHAG(PYTHAG(BVEC(ISEL,IB,1),BVEC(ISEL,IB,2)),
     &              PYTHAG(BVEC(ISEL,IB,3),BVEC(ISEL,IB,4)))
        IF(DB.GT.D0) THEN
          DO IZ = 1,NZ
            A(IZ) =  BVEC(ISEL,IB,IZ)/DB
          ENDDO
          DO I = 1,NBAS
            B(1)        = CMO(I,IC,2)*A(1) - CMO(I,IC,1)*A(2)
     &                  + CMO(I,IC,4)*A(3) - CMO(I,IC,3)*A(4)
            B(2)        = CMO(I,IC,3)*A(1) - CMO(I,IC,4)*A(2)
     &                  - CMO(I,IC,1)*A(3) + CMO(I,IC,2)*A(4)
            B(3)        = CMO(I,IC,4)*A(1) + CMO(I,IC,3)*A(2)
     &                  - CMO(I,IC,2)*A(3) - CMO(I,IC,1)*A(4)
            CMO(I,IC,1) = CMO(I,IC,1)*A(1) + CMO(I,IC,2)*A(2)
     &                  + CMO(I,IC,3)*A(3) + CMO(I,IC,4)*A(4)
            CMO(I,IC,2) = B(1)
            CMO(I,IC,3) = B(2)
            CMO(I,IC,4) = B(3)
          ENDDO
        ENDIF
      ELSE
        WRITE(LUPRI,'(A,I5)') 'PRJPH1: * Unknown NZ = ',NZ
        CALL QUIT('PRJPH1: Unknown NZ !')
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mayer_polarize */
      SUBROUTINE MAYER_POLARIZE(WORK,KFREE,LFREE)
C***********************************************************************
C
C     Polarize reference orbitals
C
C     Written by T. Saue Sep 2013
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbprj.h"
#include "dgroup.h"
      LOGICAL TOBE
      CHARACTER TEXT*74,FMT*6,MXFORM*6
      DIMENSION WORK(*)
      DIMENSION NSTR(2,0:2),KVEC(2),NVECS(2),IOFF(2),NQ(2),
     &          NPSH_FRAG(2),NESH_FRAG(2),NFBAS_FRAG(2)
C
      KFRSAV = KFREE    
C
C     Find number of MO-orbitals to analyse
C
      NCDIM = 0
      DO IFRP = 1,NFSYM
        NSTR(IFRP,1) = NESH(IFRP)
        NSTR(IFRP,2) = NPSH(IFRP)
        CALL ORBNUM(VECPRJ(IFRP),IFRP,KVEC(IFRP),NSTR(1,0),
     &              WORK,KFREE,LFREE)
        IOFF(IFRP)   = 0
        NVECS(IFRP)  = NSTR(IFRP,0)
        NQ(IFRP)     = NFBAS(IFRP,0)*NVECS(IFRP)*NZ
        NCDIM        = NCDIM + NQ(IFRP)
      ENDDO
C
C     Get selected molecular coefficients
C     
      CALL MEMGET('REAL',KCSEL,NCDIM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      IOPT = 2
      IF(IPRPRJ.GE.2) IOPT = 3
      CALL SELFRAG(WORK(KCSEL),NVECS,DUM,DUM,
     &             1,NSTR(1,0),'DFCOEF',.FALSE.,1,IOPT,
     &             KVEC,WORK(KCMO),DUM,DUM,
     &             IOFF,KRMC_FLG,WORK,KFREE,LFREE)
      CALL MEMREL('MAYERPOL.cfs',WORK,KFRSAV,KCMO,KFREE,LFREE)
C
C     Get overlap matrix packed on fermion ircops
C
      CALL MEMGET('REAL',KSMAT,N2BAST,WORK,KFREE,LFREE)
      CALL GTOVLT(WORK(KSMAT),SSMTRC)
C
C     Loop over fragments
C
      INUC = 1
      DO IFRAG = 1,NREFS
        NNUC = NPROJNUC(IFRAG)
        INQUIRE(FILE=REFFIL(IFRAG),EXIST=TOBE)
        IF (.NOT.TOBE) THEN
           WRITE(LUPRI,'(A,A6,A)') 
     &     'Fragment coefficient file ',REFFIL(IFRAG),' not present'
           CALL QUIT('MAYERPOL: No fragment coefficients !')
        ENDIF
        CALL OPNFIL(LUCOEF,REFFIL(IFRAG),'OLD','MAYERPOL')
        READ (LUCOEF,END=10,ERR=20) TEXT,NSYM,
     &       (NPSH_FRAG(I),NESH_FRAG(I),NFBAS_FRAG(I),I=1,NSYM),TOTERG
        FMT = MXFORM(TOTERG,20)
        WRITE(LUPRI,'(/A,A6,3X,A,'//FMT//')')
     &   'MAYERPOL: Coefficients read from file ',REFFIL(IFRAG),
     &   '- Total energy: ',TOTERG
        WRITE(LUPRI,'(2A)') '* Heading :',TEXT
        WRITE(LUPRI, '(2A, I8)') '- number of symmetry independent',  
     &   ' nuclei in this fragment : ', NPROJNUC(IFRAG)      
        WRITE(LUPRI,'(A,2I8)')
     &   '- Positrons : ',(NPSH_FRAG(I),I=1,NFSYM)
        WRITE(LUPRI,'(A,2I8)')
     &   '- Electrons : ',(NESH_FRAG(I),I=1,NFSYM)
        WRITE(LUPRI,'(A,2I8)')
     &   '- SO-basis  : ',(NFBAS_FRAG(I),I=1,NFSYM)
C
C       Calculate dimensions
C       ====================
C
        NCFRAG = 0
        NEFRAG = 0
        DO IFRP = 1,NFSYM
          NORB2  = NPSH_FRAG(IFRP)+NESH_FRAG(IFRP)
          NEFRAG = NEFRAG + NORB2
          NCFRAG = NCFRAG + NORB2*NFBAS_FRAG(IFRP)
        ENDDO
        NCFRAG = NCFRAG*NZ
        CALL MEMGET('REAL',KCMO2,NCFRAG,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIG2,NEFRAG,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBE2,NEFRAG,WORK,KFREE,LFREE)
        CALL READT(LUCOEF,NCFRAG,WORK(KCMO2))
        CALL READT(LUCOEF,NEFRAG,WORK(KEIG2))
        CALL READI(LUCOEF,NEFRAG,WORK(KIBE2))
        CLOSE(LUCOEF,STATUS='KEEP')
        IC1 = KCSEL
        IC2 = KCMO2
        IS1= KSMAT
        DO IFRP = 1,NFSYM
C.........Generate pointer array for from fragment basis to molecular basis
          CALL MEMGET('INTE',KBAS,NFBAS_FRAG(IFRP),WORK,KFREE,LFREE)
          CALL SELOWI(INUC,NNUC,IFRP,WORK(KBAS),NSBAS)
          IF(NSBAS.NE.NFBAS_FRAG(IFRP)) THEN
            WRITE(LUPRI,'(A,A,I3/A,A4,A,I3)')
     &        'POLARIZE)_REF: ',
     &        'Error in selection of coefficients in ircop ',
     &        IFRP,
     &        'for fragment ',REFFIL(IFRAG),' no. ',IFRAG
            WRITE(LUPRI,'(A,I8)')
     &        'Number of basis functions for this fragment is:',NSBAS,
     &        'Number of basis functions in coefficient file:',
     &         NFBAS_FRAG(IFRP)
           CALL QUIT('MAYERPOL: Error in cf.selection !')
          ENDIF
C.........Generate polaried fragment orbitals: restrict to posivite-energy orbitals
          IOFF2  = NPSH_FRAG(IFRP) 
          NVEC2  = NESH_FRAG(IFRP)
          NORB2  = NPSH_FRAG(IFRP)+NESH_FRAG(IFRP)
          NSCOMP = NVEC2*NVECS(IFRP)*NZ
          CALL MEMGET('REAL',KSCOMP,NSCOMP,WORK,KFREE,LFREE)
          NCHUNK = NFBAS_FRAG(IFRP)*NVECS(IFRP)*NZ
          CALL MEMGET('REAL',KCHUNK,NCHUNK,WORK,KFREE,LFREE)
          CALL MPOLFRAG(WORK(IC1),NFBAS(IFRP,0),NVECS(IFRP),
     &                 WORK(IC2),NFBAS_FRAG(IFRP),NORB2,IOFF2,NVEC2,
     &                 WORK(IS1),WORK(KSCOMP),WORK(KCHUNK),WORK(KBAS),
     &                 WORK,KFREE,LFREE,IPRPRJ)
          IC1 = IC1 + NQ(IFRP)
          IC2 = IC2 + NORB2*NFBAS_FRAG(IFRP)*NZ
          IS1 = IS1 + NFBAS(IFRP,0)*NFBAS(IFRP,0)
          CALL MEMREL('MAYERPOL.if',WORK,KFRSAV,KBAS,KFREE,LFREE)
        ENDDO ! End loop over fermion irreps IFRP
C.......Dump polarized coefficients to file
        REFFIL(IFRAG)(1:1)='P'
        CALL OPNFIL(LUCOEF,REFFIL(IFRAG),'NEW','MAYERPOL')
        CALL GTINFO(TEXT(51:74))
        WRITE (LUCOEF) TEXT,NSYM,
     &       (NPSH_FRAG(I),NESH_FRAG(I),NFBAS_FRAG(I),I=1,NSYM),TOTERG
        CALL WRITT(LUCOEF,NCFRAG,WORK(KCMO2))
        CALL WRITT(LUCOEF,NEFRAG,WORK(KEIG2))
        CALL WRITI(LUCOEF,NEFRAG,WORK(KIBE2))
        CLOSE(LUCOEF,STATUS='KEEP')
        INUC = INUC + NNUC
      ENDDO ! End loop over REFS
      CALL MEMREL('MAYERPOL',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
 10   CONTINUE
      CALL QUIT('MAYERPOL: END reading TEXT')
 20   CONTINUE
      CALL QUIT('MAYERPOL: ERROR reading TEXT')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mpolfrag */
      SUBROUTINE MPOLFRAG(CMO1,NBAS1,NORB1,CMO2,NBAS2,NORB2,IOFF2,NVEC2,
     &                   SMAT1,SCOMP,CHUNK,IBAS,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,DM1=-1.0D0,D0 = 0.0D0)
#include "dgroup.h"
C
      DIMENSION CMO1(NBAS1,NORB1,NZ),CMO2(NBAS2,NORB2,NZ),
     &          SMAT1(NBAS1,NBAS1),SCOMP(NVEC2,NORB1,NZ),
     &          CHUNK(NBAS2,NORB1,NZ),IBAS(NBAS2),WORK(*)
C
      KFRSAV=KFREE      
C.....Extract fragment chunk from selected molecular coefficients
      NCOL=NORB1*NZ
      CALL EXTVEC(CMO1,NBAS1,CHUNK,NBAS2,NBAS2,NCOL,IBAS)
      IF(IPRINT.GE.4) THEN
        WRITE(LUPRI,*) '* Fragment chunk'
        CALL PRQMAT(CHUNK,NBAS2,NORB1,NBAS2,NORB1,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C.....Extract fragment chunk of overlap matrix
      NSMAT2 = NBAS2*NBAS2
      CALL MEMGET('REAL',KSMAT2,NSMAT2,WORK,KFREE,LFREE)
      CALL EXTELM('N',SMAT1,NBAS1,NBAS1,WORK(KSMAT2),NBAS2,NBAS2,NBAS2,
     &            IBAS,IBAS)
C.....Transform to overlap in fragment chunk basis
      NSMATC  = NORB1*NORB1*NZ
      CALL MEMGET('REAL',KSMATC,NSMATC,WORK,KFREE,LFREE)
      CALL QTRANS('AOMO','S',D0,NBAS2,NBAS2,NORB1,NORB1,
     &           WORK(KSMAT2),NBAS2,NBAS2,1,IPQTOQ(1,0),
     &           WORK(KSMATC),NORB1,NORB1,NZ,IPQTOQ(1,0),
     &           CHUNK,NBAS2,NORB1,NZ,IPQTOQ(1,0),
     &           CHUNK,NBAS2,NORB1,NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,*) '* Overlap matrix in fragment chunk basis'
        CALL PRQMAT(WORK(KSMATC),NORB1,NORB1,NORB1,NORB1,
     &               NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
      IF(NORB1.EQ.1) THEN
C.......Renormalization
        NEFF = 1
        FAC = D1/SQRT(WORK(KSMATC))
        CALL DSCAL(NCHUNK,FAC,CHUNK,1)
      ELSE
C.......Canonical orthonormalization + realigment
        CALL MEMGET('REAL',KVMATC,NSMATC,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KSEIG,NORB1,WORK,KFREE,LFREE)
        STOL=1.0D-2
        CALL LOWGENQ(NZ,NORB1,NEFF,
     &               WORK(KSMATC),NORB1,NORB1,
     &               WORK(KVMATC),NORB1,NORB1,
     &               STOL,IPRINT,WORK(KSEIG),IPQTOQ(1,0),
     &               WORK(KFREE),LFREE)
        IF(IPRINT.GE.2) THEN
          WRITE(6,*) 'Eigenvalues of overlap matrix',NEFF
          CALL OUTPUT(WORK(KSEIG),1,NORB1,1,1,NORB1,1,-1,LUPRI)
        ENDIF
        WRITE(6,*) 'Number of polarized atomic orbitals: ',NEFF
        NCHUNK=NBAS2*NORB1*NZ
        CALL MEMGET('REAL',KBUF,NCHUNK,WORK,KFREE,LFREE)
        CALL QGEMM(NBAS2,NEFF,NORB1,D1,
     &             'N','N',IPQTOQ(1,0),CHUNK,NBAS2,NORB1,NZ,
     &             'N','N',IPQTOQ(1,0),WORK(KVMATC),NORB1,NORB1,NZ,
     &             D0,IPQTOQ(1,0),WORK(KBUF),NBAS2,NORB1,NZ)
        CALL DCOPY(NCHUNK,WORK(KBUF),1,CHUNK,1)
        IF(IPRINT.GE.2) THEN
          CALL QTRANS('AOMO','S',D0,NBAS2,NBAS2,NEFF,NEFF,
     &             WORK(KSMAT2),NBAS2,NBAS2,1,IPQTOQ(1,0),
     &             WORK(KSMATC),NEFF,NEFF,NZ,IPQTOQ(1,0),
     &             CHUNK,NBAS2,NORB1,NZ,IPQTOQ(1,0),
     &             CHUNK,NBAS2,NORB1,NZ,IPQTOQ(1,0),
     &             WORK(KFREE),LFREE,IPRINT)
          WRITE(LUPRI,*) 
     &       '* Overlap matrix in orthonormalized chunk basis'
          CALL PRQMAT(WORK(KSMATC),NEFF,NEFF,NEFF,NEFF,
     &            NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C.......Align free and polarized orbitals based on corresponding orbitals
        CALL QTRANS('AOMO','S',D0,NBAS2,NBAS2,NVEC2,NEFF,
     &             WORK(KSMAT2),NBAS2,NBAS2,1,IPQTOQ(1,0),
     &             SCOMP,NVEC2,NORB1,NZ,IPQTOQ(1,0),
     &             CMO2(1,IOFF2+1,1),NBAS2,NORB2,NZ,IPQTOQ(1,0),
     &             CHUNK,NBAS2,NORB1,NZ,IPQTOQ(1,0),
     &             WORK(KFREE),LFREE,IPRINT)
        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,*) 
     &    '* Overlap free/polarized fragment orbitals before alignment:'
          CALL PRQMAT(SCOMP,NVEC2,NEFF,NVEC2,NORB1,NZ,IPQTOQ(1,0),
     &            LUPRI)
        ENDIF
        CALL QGEMM(NEFF,NEFF,NEFF,D1,
     &             'H','N',IPQTOQ(1,0),SCOMP,NVEC2,NORB1,NZ,
     &             'N','N',IPQTOQ(1,0),SCOMP,NVEC2,NORB1,NZ,
     &             D0,IPQTOQ(1,0),WORK(KSMATC),NEFF,NEFF,NZ)
        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,*) '* SCOMP squared ...'
          CALL PRQMAT(WORK(KSMATC),NEFF,NEFF,NEFF,NEFF,
     &                NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C.......Scale matrix with -1 to get most overlapping orbitals first
        NDIM=NEFF*NEFF*NZ
        CALL DSCAL(NDIM,DM1,WORK(KSMATC),1)
        CALL QDIAG(NZ,NEFF,WORK(KSMATC),NEFF,NEFF,
     &             WORK(KSEIG),1,WORK(KVMATC),NEFF,
     &             NEFF,WORK(KFREE),LFREE,IERR)
        IF(IPRINT.GE.2) THEN
          WRITE(6,*) 'Eigenvalues of SCOMP squared'
          CALL OUTPUT(WORK(KSEIG),1,NEFF,1,1,NEFF,1,-1,LUPRI)
        ENDIF
        CALL QGEMM(NBAS2,NEFF,NEFF,D1,
     &             'N','N',IPQTOQ(1,0),CHUNK,NBAS2,NORB1,NZ,
     &             'N','N',IPQTOQ(1,0),WORK(KVMATC),NEFF,NEFF,NZ,
     &             D0,IPQTOQ(1,0),WORK(KBUF),NBAS2,NORB1,NZ)
        CALL DCOPY(NCHUNK,WORK(KBUF),1,CHUNK,1)
C.......Calculate overlap between free- and polarized fragment orbitals
        CALL QTRANS('AOMO','S',D0,NBAS2,NBAS2,NVEC2,NEFF,
     &             WORK(KSMAT2),NBAS2,NBAS2,1,IPQTOQ(1,0),
     &             SCOMP,NVEC2,NORB1,NZ,IPQTOQ(1,0),
     &             CMO2(1,IOFF2+1,1),NBAS2,NORB2,NZ,IPQTOQ(1,0),
     &             CHUNK,NBAS2,NORB1,NZ,IPQTOQ(1,0),
     &             WORK(KFREE),LFREE,IPRINT)
        IF(IPRINT.GE.2) THEN
          CALL QGEMM(NEFF,NEFF,NEFF,D1,
     &               'H','N',IPQTOQ(1,0),SCOMP,NVEC2,NORB1,NZ,
     &               'N','N',IPQTOQ(1,0),SCOMP,NVEC2,NORB1,NZ,
     &               D0,IPQTOQ(1,0),WORK(KSMATC),NEFF,NEFF,NZ)
          WRITE(LUPRI,*) '* SCOMP squared again ...'
          CALL PRQMAT(WORK(KSMATC),NEFF,NEFF,NEFF,NEFF,
     &              NZ,IPQTOQ(1,0),LUPRI)                
        ENDIF
      ENDIF
      CALL MEMREL('MPOLFRAG',WORK,KFRSAV,KSMATC,KFREE,LFREE)
C.....Calculate overlap between free- and polarized fragment orbitals
      CALL QTRANS('AOMO','S',D0,NBAS2,NBAS2,NVEC2,NEFF,
     &           WORK(KSMAT2),NBAS2,NBAS2,1,IPQTOQ(1,0),
     &           SCOMP,NVEC2,NORB1,NZ,IPQTOQ(1,0),
     &           CMO2(1,IOFF2+1,1),NBAS2,NORB2,NZ,IPQTOQ(1,0),
     &           CHUNK,NBAS2,NORB1,NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,*) 
     &     '* Overlap between free/polarized fragment orbitals '
        CALL PRQMAT(SCOMP,NVEC2,NEFF,NVEC2,NORB1,NZ,IPQTOQ(1,0),
     &              LUPRI)
      ELSE
        WRITE(LUPRI,*) 
     &     '* Overlap between free/polarized fragment orbitals '
        CALL PRQMAT(SCOMP,NEFF,NEFF,NVEC2,NORB1,NZ,IPQTOQ(1,0),
     &              LUPRI)
      ENDIF
C.....Use overlap selection to put polarized orbitals in corresponding free position
      DO J = 1, NEFF
        DO I = 1,NVEC2
           SCOMP(I,J,1) = SCOMP(I,J,1)*SCOMP(I,J,1)
        ENDDO
      ENDDO
      DO IZ = 2,NZ
        DO J = 1, NEFF
          DO I = 1,NVEC2
             SCOMP(I,J,1) = SCOMP(I,J,1)+SCOMP(I,J,IZ)*SCOMP(I,J,IZ)
          ENDDO
      ENDDO
      ENDDO
      DO J = 1,NEFF
        ISEL = IDAMAX(NVEC2,SCOMP(1,J,1),1)
C       Zero out the selected overlap
        CALL DCOPY(NORB1,D0,0,SCOMP(ISEL,1,1),NVEC2)
C       Copy polarized orbitals into corresponding free orbital position
        DO IZ = 1,NZ
          CALL DCOPY(NBAS2,CHUNK(1,J,IZ),1,CMO2(1,IOFF2+ISEL,IZ),1)
        ENDDO
      ENDDO
      CALL MEMREL('MPOLFRAG',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck knizia */
      SUBROUTINE KNIZIA(CMO,BVEC,SMAT,AMAT,NTOT,NMOL,NBAS,NZ,IPQ,
     &                  NREFS,REFFIL,NFRAGV,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Polarize fragment orbitals according to the intrinsic AO scheme
C     propozed by Knizia, see
C
C     G. Knizia , J. Chem. Theory Comp. 
C     "Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum Theory and Chemical Concepts
C
C     Written by T. Saue - September 2013
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,DM1=-1.0D0,D1=1.0D0)
C
      CHARACTER*6 REFFIL(NREFS)
      DIMENSION CMO(NBAS,NTOT,NZ),BVEC(NTOT,NMOL,NZ),SMAT(NBAS,NBAS),
     &          NFRAGV(NREFS),AMAT(NTOT,NTOT,NZ),WORK(*)
C
      KFRSAV=KFREE
      NREF=NTOT-NMOL
      IMOL=NREF+1
C
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,*) '* Starting reference coefficients'
        CALL PRQMAT(CMO,NBAS,NREF,NBAS,NTOT,
     &          NZ,IPQ,LUPRI)
        CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NREF,NREF,
     &           SMAT,NBAS,NBAS,1,IPQ,
     &           AMAT,NREF,NREF,NZ,IPQ,
     &           CMO,NBAS,NTOT,NZ,IPQ,
     &           CMO,NBAS,NTOT,NZ,IPQ,
     &           WORK(KFREE),LFREE,IPRINT)
        WRITE(LUPRI,*) 
     &     '* Overlap matrix of starting coefficients'
        CALL PRQMAT(AMAT,NREF,NREF,NREF,NREF,
     &          NZ,IPQ,LUPRI)
      ENDIF
C.....Generate depolarized MOs (stored in WORK(KBUF))
      NDIM=NBAS*NMOL*NZ
      CALL MEMGET('REAL',KBUF,NDIM,WORK,KFREE,LFREE)
      CALL QGEMM(NBAS,NMOL,NREF,D1,
     &       'N','N',IPQ,CMO,NBAS,NTOT,NZ,
     &       'N','N',IPQ,BVEC,NTOT,NMOL,NZ,
     &       D0,IPQ,WORK(KBUF),NBAS,NMOL,NZ)
C.....Orthogonalize depolarized MOs
      CALL ORTHO_SELVEC(WORK(KBUF),NBAS,NMOL,1,NMOL,NEFF,SMAT,AMAT,
     &                  IPQ,NZ,IPRINT,WORK,KFREE,LFREE)
      IF(NEFF.NE.NMOL) THEN
        WRITE(6,*) '* WARNING ! KNIZIA: Linear dependencies removed ',
     &    'in depolarized MOs. NMOL = ',NMOL,' NEFF = ',NEFF
      ENDIF
      IF(IPRINT.GE.5) THEN
        CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NEFF,NEFF,
     &           SMAT,NBAS,NBAS,1,IPQ,
     &           AMAT,NMOL,NMOL,NZ,IPQ,
     &           WORK(KBUF),NBAS,NMOL,NZ,IPQ,
     &           WORK(KBUF),NBAS,NMOL,NZ,IPQ,
     &           WORK(KFREE),LFREE,IPRINT)
        WRITE(LUPRI,*) 
     &     '* Overlap matrix in orthonormalized depolarized MOs'
        CALL PRQMAT(AMAT,NEFF,NEFF,NMOL,NMOL,
     &          NZ,IPQ,LUPRI)
      ENDIF
C
C.....Construct projectors
C     ====================
      NDIM = NBAS*NREF
      NVEC = NDIM*NZ
      NMAT = NBAS*NBAS*NZ
      CALL MEMGET('REAL',KM1,NMAT,WORK,KFREE,LFREE)            
      CALL MEMGET('REAL',KV1,NVEC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KM2,NMAT,WORK,KFREE,LFREE)
C.....Projector onto depolarized space P' = D'S -> store in M1
      CALL DENST1(WORK(KM2),NBAS,NBAS,NZ,D1,D0,
     &            WORK(KBUF),NBAS,NMOL,1,NEFF,NBAS)
      CALL QGEMM(NBAS,NBAS,NBAS,D1,
     &           'N','N',IPQ,WORK(KM2),NBAS,NBAS,NZ,
     &           'N','N',IPQ,SMAT,NBAS,NBAS,1,
     &           D0,IPQ,WORK(KM1),NBAS,NBAS,NZ)
      CALL MEMREL('KNIZIA.proj1',WORK,KFRSAV,KM2,KFREE,LFREE)
      IF(IPRINT.GE.5) THEN
        WRITE(LUPRI,*) '* Projector onto depolarized space'
        CALL PRQMAT(WORK(KM1),NBAS,NBAS,NBAS,NBAS,
     &          NZ,IPQ,LUPRI)
      ENDIF
C.....apply first projector to fragment coefficients: b11 = P'b0; stored in V1
      CALL QGEMM(NBAS,NREF,NBAS,D1,
     &           'N','N',IPQ,WORK(KM1),NBAS,NBAS,NZ,
     &           'N','N',IPQ,CMO,NBAS,NTOT,NZ,
     &           D0,IPQ,WORK(KV1),NBAS,NREF,NZ)
      IF(IPRINT.GE.8) THEN
        WRITE(LUPRI,*) '* Intermediate b11 '
        CALL PRQMAT(WORK(KV1),NBAS,NREF,NBAS,NREF,
     &          NZ,IPQ,LUPRI)
      ENDIF
C.....construct orthogonal complement: b12 = b0 - b11; stored in v2
      CALL MEMGET('REAL',KV2,NVEC,WORK,KFREE,LFREE)
      IOFF = 0
      DO IZ = 1,NZ
        CALL DCOPY(NDIM,CMO(1,1,IZ),1,WORK(KV2+IOFF),1)
        IOFF = IOFF + NDIM
      ENDDO
      CALL DAXPY(NDIM*NZ,DM1,WORK(KV1),1,WORK(KV2),1)
      IF(IPRINT.GE.8) THEN
        WRITE(LUPRI,*) '* Intermediate b12 '
        CALL PRQMAT(WORK(KV2),NBAS,NREF,NBAS,NREF,
     &          NZ,IPQ,LUPRI)
      ENDIF
C.....Projector onto starting MO space P = DS -> store in M1
      CALL MEMGET('REAL',KM2,NMAT,WORK,KFREE,LFREE)
      CALL DENST1(WORK(KM2),NBAS,NBAS,NZ,D1,D0,
     &            CMO,NBAS,NTOT,IMOL,NMOL,NBAS)
      CALL QGEMM(NBAS,NBAS,NBAS,D1,
     &           'N','N',IPQ,WORK(KM2),NBAS,NBAS,NZ,
     &           'N','N',IPQ,SMAT,NBAS,NBAS,1,
     &           D0,IPQ,WORK(KM1),NBAS,NBAS,NZ)
      CALL MEMREL('KNIZIA.proj',WORK,KFRSAV,KM2,KFREE,LFREE)
      IF(IPRINT.GE.5) THEN
        WRITE(LUPRI,*) '* Projector onto starting MO space'
        CALL PRQMAT(WORK(KM1),NBAS,NBAS,NBAS,NBAS,
     &          NZ,IPQ,LUPRI)
      ENDIF
      CALL MEMGET('REAL',KV3,NVEC,WORK,KFREE,LFREE)
C.....apply: b21 = P*b11; store in v3
      CALL QGEMM(NBAS,NREF,NBAS,D1,
     &           'N','N',IPQ,WORK(KM1),NBAS,NBAS,NZ,
     &           'N','N',IPQ,WORK(KV1),NBAS,NREF,NZ,
     &           D0,IPQ,WORK(KV3),NBAS,NREF,NZ)
      IF(IPRINT.GE.8) THEN
        WRITE(LUPRI,*) '* Intermediate b21 '
        CALL PRQMAT(WORK(KV3),NBAS,NREF,NBAS,NREF,
     &          NZ,IPQ,LUPRI)
      ENDIF
C.....copy b12 to V1
      CALL DCOPY(NDIM*NZ,WORK(KV2),1,WORK(KV1),1)
C.....subtract P*b12 from b12 -> v1
      CALL QGEMM(NBAS,NREF,NBAS,DM1,
     &           'N','N',IPQ,WORK(KM1),NBAS,NBAS,NZ,
     &           'N','N',IPQ,WORK(KV2),NBAS,NREF,NZ,
     &           D1,IPQ,WORK(KV1),NBAS,NREF,NZ)
      IF(IPRINT.GE.8) THEN
        WRITE(LUPRI,*) '* Intermediate b22 '
        CALL PRQMAT(WORK(KV1),NBAS,NREF,NBAS,NREF,
     &          NZ,IPQ,LUPRI)
      ENDIF
      IF(IPRINT.GE.5) THEN
        WRITE(LUPRI,*) '* Polarized reference coefficients - P'
        CALL PRQMAT(WORK(KV3),NBAS,NREF,NBAS,NREF,
     &          NZ,IPQ,LUPRI)
        WRITE(LUPRI,*) '* Polarized reference coefficients - Q'
        CALL PRQMAT(WORK(KV1),NBAS,NREF,NBAS,NREF,
     &          NZ,IPQ,LUPRI)
      ENDIF
C.....Assemble and normalize polarized fragment coefficients
      CALL DAXPY(NDIM*NZ,D1,WORK(KV3),1,WORK(KV1),1)
      CALL MEMREL('KNIZIA.ass',WORK,KFRSAV,KV2,KFREE,LFREE)
C
C     Calculate overlap between original and polarized coefficients
C
      CALL ABSOLUTE_OVERLAP(CMO,NBAS,NTOT,1,NREF,
     &                      WORK(KV1),NREF,1,SMAT,AMAT,
     &                      IPQ,NZ,IPRINT,WORK,KFREE,LFREE)
      WRITE(LUPRI,*) 
     &  '* Squared overlap between original and polarized orbitals',
     &  ' before orthonormalization.'
      CALL PRQMAT(AMAT,NREF,NREF,NREF,NREF,1,IPQ,LUPRI)
C
C     Orthonormalization and align within each fragment
C
      IVEC= 1
      DO IREF = 1,NREFS
        NVEC = NFRAGV(IREF)
        CALL ORTHO_SELVEC(WORK(KV1),NBAS,NREF,IVEC,NVEC,NEFF,SMAT,AMAT,
     &                    IPQ,NZ,IPRINT,WORK,KFREE,LFREE)
        IVEC = IVEC + NVEC
      ENDDO
      CALL ABSOLUTE_OVERLAP(CMO,NBAS,NTOT,1,NREF,
     &                      WORK(KV1),NREF,1,SMAT,AMAT,
     &                      IPQ,NZ,IPRINT,WORK,KFREE,LFREE)
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,*) 
     &  '* Squared overlap between original and polarized orbitals',
     &  ' after orthonormalization.'
        CALL PRQMAT(AMAT,NREF,NREF,NREF,NREF,1,IPQ,LUPRI)
      ENDIF
C.....Replace original polarization coefficients
      IOFF = 0
      DO IZ = 1,NZ
        CALL DCOPY(NDIM,WORK(KV1+IOFF),1,CMO(1,1,IZ),1)
        IOFF = IOFF + NDIM
      ENDDO
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,*) '* Polarized reference coefficients'
        CALL PRQMAT(CMO,NBAS,NREF,NBAS,NTOT,
     &          NZ,IPQ,LUPRI)
      ENDIF
      CALL MEMREL('KNIZIA',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Ortho_selvec */
      SUBROUTINE ORTHO_SELVEC(CMO,NBAS,NVEC,ISEL,NSEL,NEFF,
     &                    SMAT,AMAT,IPQ,NZ,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Orthogonalize selected vectors
C     Written by Trond Saue Sep 2013
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      DIMENSION CMO(NBAS,NVEC,NZ),SMAT(NBAS,NBAS),
     &          AMAT(NSEL,NSEL,NZ),IPQ(NZ),WORK(*)
C
      KFRSAV = KFREE
      CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NSEL,NSEL,
     &         SMAT,NBAS,NBAS,1,IPQ,
     &         AMAT,NSEL,NSEL,NZ,IPQ,
     &         CMO(1,ISEL,1),NBAS,NVEC,NZ,IPQ,
     &         CMO(1,ISEL,1),NBAS,NVEC,NZ,IPQ,
     &         WORK(KFREE),LFREE,IPRINT)
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,*) '* ORTHO_SELVEC: Overlap matrix'
        CALL PRQMAT(AMAT,NSEL,NSEL,NSEL,NSEL,NZ,IPQ,LUPRI)
      ENDIF
      IF(NSEL.EQ.1) THEN
        NEFF = 1
        FAC = D1/SQRT(AMAT(1,1,1))
        DO IZ = 1,NZ
          CALL DSCAL(NBAS,FAC,CMO(1,ISEL,IZ),1)
        ENDDO
      ELSE
        NDIM=NSEL*NSEL*NZ
        CALL MEMGET('REAL',KEIG ,NSEL,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVMAT,NDIM,WORK,KFREE,LFREE)
        STOL=1.0D-2
        CALL SYM_ORTHO(NZ,NSEL,NEFF,
     &             AMAT       ,NSEL,NSEL,
     &             WORK(KVMAT),NSEL,NSEL,
     &             STOL,IPRINT,WORK(KEIG),IPQ,
     &             WORK(KFREE),LFREE)
CTROND        CALL LOWGENQ(NZ,NSEL,NEFF,
CTROND     &             AMAT       ,NSEL,NSEL,
CTROND     &             WORK(KVMAT),NSEL,NSEL,
CTROND     &             STOL,IPRINT,WORK(KEIG),IPQ,
CTROND     &             WORK(KFREE),LFREE)
        IF(IPRINT.GE.2) THEN
          WRITE(6,*) 'ORTHO_SELVEC: Eigenvalues of overlap matrix',NEFF
          CALL OUTPUT(WORK(KEIG),1,NSEL,1,1,NSEL,1,-1,LUPRI)
          WRITE(6,*) 'ORTHO_SELVEC: Span of selected coefficients ',NEFF
        ENDIF
        NDIM=NBAS*NSEL
        CALL MEMGET('REAL',KBUF,NDIM*NZ,WORK,KFREE,LFREE)
        IOFF = 0
        DO IZ = 1,NZ
          CALL DCOPY(NDIM,CMO(1,ISEL,IZ),1,WORK(KBUF+IOFF),1)
          IOFF = IOFF + NDIM
        ENDDO
        CALL QGEMM(NBAS,NEFF,NSEL,D1,
     &             'N','N',IPQ,WORK(KBUF),NBAS,NSEL,NZ,
     &             'N','N',IPQ,WORK(KVMAT),NSEL,NSEL,NZ,
     &             D0,IPQ,CMO(1,ISEL,1),NBAS,NVEC,NZ)
        IF(IPRINT.GE.2) THEN
          CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NEFF,NEFF,
     &         SMAT,NBAS,NBAS,1,IPQ,
     &         AMAT,NEFF,NEFF,NZ,IPQ,
     &         CMO(1,ISEL,1),NBAS,NVEC,NZ,IPQ,
     &         CMO(1,ISEL,1),NBAS,NVEC,NZ,IPQ,
     &         WORK(KFREE),LFREE,IPRINT)
          WRITE(LUPRI,*) '* ORTHO_SELVEC: Overlap matrix after'
          CALL PRQMAT(AMAT,NEFF,NEFF,NEFF,NEFF,NZ,IPQ,LUPRI)
        ENDIF
        CALL MEMREL('ORTHO_SELVEC',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Align_selvec */
      SUBROUTINE ALIGN_SELVEC(CMO,NBAS,NVEC,ISEL,NSEL,
     &                        CREF,NCOL,IREF,SMAT,AMAT,
     &                        IPQ,NZ,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Align selected vectors of CMO with reference vectors CREF
C     Written by Trond Saue Sep 2013
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,DM1=-1.0D0,D0 = 0.0D0)
      DIMENSION CMO(NBAS,NVEC,NZ),CREF(NBAS,NCOL,NZ),
     &          SMAT(NBAS,NBAS),AMAT(NSEL,NSEL,NZ),IPQ(NZ),WORK(*)
C
      KFRSAV = KFREE
      IF(IPRINT.GE.5) THEN
        WRITE(LUPRI,*) '* ALIGN_SELVEC: Starting coefficients'
        CALL PRQMAT(CMO(1,ISEL,1),NBAS,NSEL,NBAS,NVEC,NZ,IPQ,LUPRI)
        WRITE(LUPRI,*) '* ALIGN_SELVEC: Reference coefficients'
        CALL PRQMAT(CREF(1,IREF,1),NBAS,NSEL,NBAS,NCOL,NZ,IPQ,LUPRI)
      ENDIF
C.....Calculate overlap between reference and selected orbitals
      CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NSEL,NSEL,
     &         SMAT,NBAS,NBAS,1,IPQ,
     &         AMAT,NSEL,NSEL,NZ,IPQ,
     &         CREF(1,IREF,1),NBAS,NCOL,NZ,IPQ,
     &         CMO(1,ISEL,1),NBAS,NVEC,NZ,IPQ,
     &         WORK(KFREE),LFREE,IPRINT)
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,*) '* ALIGN_SELVEC: Overlap matrix ref/sel'
        CALL PRQMAT(AMAT,NSEL,NSEL,NSEL,NSEL,NZ,IPQ,LUPRI)
      ENDIF
C.....Use overlap selection to put polarized orbitals in corresponding free position
      DO J = 1, NSEL
        DO I = 1,NSEL
           AMAT(I,J,1) = AMAT(I,J,1)*AMAT(I,J,1)
        ENDDO
      ENDDO
      DO IZ = 2,NZ
        DO J = 1, NSEL
          DO I = 1,NSEL
             AMAT(I,J,1) = AMAT(I,J,1)+AMAT(I,J,IZ)*AMAT(I,J,IZ)
          ENDDO
        ENDDO
      ENDDO
      DO J = 1,NSEL
        IMAX = IDAMAX(NSEL,AMAT(1,J,1),1)
        J1 = ISEL + J - 1
        J2 = IREF + IMAX - 1
C       Zero out the selected overlap
        CALL DCOPY(NSEL,D0,0,AMAT(IMAX,1,1),NSEL)
C       Copy polarized orbitals into corresponding free orbital position
        DO IZ = 1,NZ
          CALL DCOPY(NBAS,CMO(1,J1,IZ),1,CREF(1,J2,IZ),1)
        ENDDO
      ENDDO
      CALL MEMREL('ALIGN_SELVEC',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Norm_selvec */
      SUBROUTINE NORM_SELVEC(CMO,NBAS,NVEC,ISEL,NSEL,SMAT,AMAT,
     &                        IPQ,NZ,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Normalize selected vectors
C     Written by Trond Saue Sep 2013
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,DM1=-1.0D0,D0 = 0.0D0)
      DIMENSION CMO(NBAS,NVEC,NZ),SMAT(NBAS,NBAS),AMAT(NSEL,NSEL,NZ),
     &          IPQ(NZ),WORK(*)
C
      KFRSAV = KFREE
C.....Calculate overlap matrix
      CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NSEL,NSEL,
     &         SMAT,NBAS,NBAS,1,IPQ,
     &         AMAT,NSEL,NSEL,NZ,IPQ,
     &         CMO(1,ISEL,1),NBAS,NVEC,NZ,IPQ,
     &         CMO(1,ISEL,1),NBAS,NVEC,NZ,IPQ,
     &         WORK(KFREE),LFREE,IPRINT)
      DO IVEC = ISEL,NSEL-1
        FAC = D1/SQRT(AMAT(1,1,1))
        DO IZ = 1,NZ
          CALL DSCAL(NBAS,FAC,CMO(1,IVEC,IZ),1)
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Anapr1 */
      SUBROUTINE ANAPR1(NVECS,KVEC,NSTR,CSEL,ESEL,IBSEL,BVEC,
     &                  IOPT,KRMC_FLG,WORK,KFREE,LFREE)
C***********************************************************************
C
C    Driver for projection analysis
C
C    NVECS(IFRP,1) - number of fragments orbitals in ircop IFRP
C    NVECS(IFRP,2) - number of MOs in ircop IFRP
C    NVECS(IFRP,0) - total number of orbitals in ircop IFRP
C
C    KVEC(IFRP,0)  - address for index array of indexes for MOs
C    KVEC(IFRP,IFRAG) - address for index array of indexes of orbitals of fragment IFRAG
C
C    NSTR(IFRP,0,0) - number of MOs in ircop IFRP
C    NSTR(IFRP,1,0) - number of positive-energy MOs in ircop IFRP
C    NSTR(IFRP,1,0) - number of negative-energy MOs in ircop IFRP
C    NSTR(IFRP,0,IFRAG) - number of orbitals of fragment IFRAG in ircop IFRP
C    NSTR(IFRP,1,IFRAG) - number of positive-energy orbitals of fragment IFRAG in ircop IFRP
C    NSTR(IFRP,1,IFRAG) - number of negative-energy orbitals of fragment IFRAG in ircop IFRP
C
C    CSEL, ESEL and IBSEL are arrays of coefficients, orbital eigenvalues and symmetry information
C    For each fermion ircop data for reference orbitals are followed by those of MOs
C    
C    When symmetry-adapted orbitals are NOT used, the analysis is effectively carried
C    out in C1 symmetry and so only one ircop is defined.
C      
C     KRMC_FLG = 1: Orbitals read from
C                   label NEWORB on file KRMCSCF
C     KRMC_FLG = 0  Orbitals read from DFCOEF
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D0=0.0D0)
C
#include "nuclei.h"
#include "dcbgen.h"
#include "dcbprj.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION NVECS(2,0:2),KVEC(2,0:MAXREF),NSTR(2,0:2,0:MAXREF),
     &          CSEL(*),ESEL(*),IBSEL(*),BVEC(*),
     &          GROSS(2,0:2,MAXREF),POL(2),WORK(*)
      integer, allocatable :: NFRAGV(:),JMOL(:,:),JREF(:,:)
      real(8), allocatable :: OCC(:)
C
      CALL QENTER('ANAPR1')
      KFRSAV=KFREE
      POLT = D0
      IF(PATOMS) THEN
C.....Analysis in molecular symmetry in terms of individual atoms
        KBUF   = KFREE
        NMOL   = NVECS(1,2)
        NREF   = NVECS(1,1)
        IMOL   = NREF + 1
C.......ntot   = nmol + nref
        NTOT   = NVECS(1,0)
        NBAS   = NTBAS(0)
        allocate(NFRAGV(NREFS))
        allocate(OCC(NMOL))
        allocate(JMOL(3,NMOL))
        allocate(JREF(3,NREF))
        CALL MEMGET('REAL',KSMAT,N2BBASX,WORK,KFREE,LFREE)
        CALL GTOVLX(WORK(KSMAT),SSMTRC)
        CALL MEMGET('REAL',KBBUF ,N2BBASX ,WORK,KFREE,LFREE)
        CALL MTBSBU(WORK(KSMAT),WORK(KBBUF))
        CALL MTSOAO(WORK(KBBUF),WORK(KSMAT),NTBAS(0),0,IPRPRJ)
        CALL MEMREL('ANAPR1.aoovl',WORK,1,KBBUF,KFREE,LFREE)
        IVEC = 0
        NPTOT=0
        DO IFRP = 1,NFSYM
          CALL PRJINF(JMOL,IFRP,0,IVEC,NSTR,WORK(KVEC(IFRP,0)),
     &                NPSH(IFRP))
          NPTOT = NPTOT + NPSH(IFRP)
          IVEC=IVEC + NSTR(IFRP,0,0)
        ENDDO
        IF(KRMC_FLG.GT.0) THEN
          CALL DCOPY(NMOL,ESEL(IMOL),1,OCC,1)
        ELSE
          CALL SETOCC(NMOL,OCC,JMOL)
        ENDIF
        IVEC  = 0
        DO I = 1,NREFS
          CALL PRJINF(JREF,1,I,IVEC,NSTR,WORK(KVEC(1,I)),NPTOT)
          NFRAGV(I) = NSTR(1,0,I)
          IVEC      = IVEC + NSTR(1,0,I)
        ENDDO
        CALL PRJAN1(CSEL,ESEL,IBSEL,OCC,
     &              BVEC,NTOT,NREF,NMOL,NBAS,NFRAGV,FREP,
     &              IQDEF,1,NZC1,JMOL,JREF,
     &              NSTR,WORK(KSMAT),GROSS,POL,IOPT,KRMC_FLG,IPRPRJ,
     &              WORK,KFREE,LFREE)
        CALL MEMREL('ANAPR1.aoloop',WORK,KBUF,KBUF,KFREE,LFREE)
        deallocate(NFRAGV)
        deallocate(OCC)
        deallocate(JMOL)
        deallocate(JREF)
      ELSE
C.....Analysis in terms of symmetry-adapted fragments
        CALL MEMGET('REAL',KSMAT,N2BAST,WORK,KFREE,LFREE)
        CALL GTOVLT(WORK(KSMAT),SSMTRC,0)
        ISOFF = KSMAT
        ICOFF = 1
        IVOFF = 1
        IBOFF = 1
        NBDIM = 0
        NCDIM = 0
        DO IFRP = 1,NFSYM
          NMOL   = NVECS(IFRP,2)
          NREF   = NVECS(IFRP,1)
          NTOT   = NVECS(IFRP,0)
          NBAS   = NFBAS(IFRP,0)
          IF(NMOL*NREF.NE.0) THEN
            KBUF   = KFREE
C...........ntot   = nmol + nref
            NBDIM  = NBDIM + NTOT*NMOL*NZ
            NCDIM  = NCDIM + NBAS*NTOT*NZ
            allocate(NFRAGV(NREFS))
            allocate(OCC(NMOL))
            allocate(JMOL(3,NMOL))
            allocate(JREF(3,NREF))
            CALL PRJINF(JMOL,IFRP,0,0,NSTR,WORK(KVEC(IFRP,0)),
     &                  NPSH(IFRP))
            IF(KRMC_FLG.GT.0) THEN
              CALL DCOPY(NMOL,ESEL(IVOFF+NREF),1,OCC(IOCC),1)
            ELSE
              CALL SETOCC(NMOL,OCC,JMOL)
            ENDIF
            IVEC = 0
            DO I = 1,NREFS
              CALL PRJINF(JREF,IFRP,I,IVEC,NSTR,WORK(KVEC(IFRP,I)),
     &                    NPSH(IFRP))
              NFRAGV(I) = NSTR(IFRP,0,I)
              IVEC = IVEC + NSTR(IFRP,0,I)
            ENDDO
            CALL PRJAN1(CSEL(ICOFF),ESEL(IVOFF),IBSEL(IVOFF),OCC,
     &                  BVEC(IBOFF),NTOT,NREF,NMOL,NBAS,NFRAGV,FREP,
     &                  IPQTOQ(1,0),IFRP,NZ,JMOL,JREF,
     &                  NSTR,WORK(ISOFF),GROSS,POL(IFRP),IOPT,KRMC_FLG,
     &                  IPRPRJ,WORK,KFREE,LFREE)
            CALL MEMREL('ANAPR1.soloop',WORK,KBUF,KBUF,KFREE,LFREE)
            deallocate(NFRAGV)
            deallocate(OCC)
            deallocate(JMOL)
            deallocate(JREF)
          ENDIF
          POLT  = POLT + POL(IFRP)
          ISOFF = ISOFF + N2BAS(IFRP)
          ICOFF = ICOFF + NBAS*NTOT*NZ
          IVOFF = IVOFF + NTOT
          IBOFF = IBOFF + NTOT*NMOL*NZ
        ENDDO
 5      IF(NFSYM.EQ.2) THEN
          CALL HEADER('Total gross contributions:',-1)
          GTOT = D0
          DO I = 1,NREFS
            WRITE(LUPRI,'(3X,A6,8X,3X,F7.4,2(3X,A4,F7.4))')
     &          REFFIL(I),(GROSS(1,0,I)+GROSS(2,0,I)),
     &          'E - '   ,(GROSS(1,2,I)+GROSS(2,2,I)),
     &          'P - '   ,(GROSS(1,1,I)+GROSS(2,1,I))
            GTOT = GTOT + GROSS(1,0,I)+GROSS(2,0,I)
          ENDDO
          WRITE(LUPRI,'(3X,A,F10.4)') 'Polarization: ',POLT
        ENDIF
C
C       Dump reference and projection coefficients to file
C       ==================================================
C
        CALL OPNFIL(LUCOEF,'DFPRJC','UNKNOWN','ANAPR1')
        WRITE(LUCOEF) NREFS
        DO IFRP = 1,NFSYM
C         Dimensions....NREF..........NMOL.........
          WRITE(LUCOEF) NVECS(IFRP,1),NVECS(IFRP,2),
     &                 (NSTR(IFRP,0,I),I=1,NREFS)
        ENDDO
C.......Reference + polarization coefficients
        CALL WRITT(LUCOEF,NCDIM,CSEL)
C.......Projection coefficients
        CALL WRITT(LUCOEF,NBDIM,BVEC)
        CLOSE(LUCOEF,STATUS='KEEP')
      ENDIF
C
      CALL QEXIT('ANAPR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Normpol */
      SUBROUTINE NORMPOL(POLVEC,AMAT,BVEC,NBAS,NTOT,NMOL,IPRINT)
C***********************************************************************
C
C     Normalize polarization contributions
C     POLVEC contains
C       NREF reference orbitals + NMOL polarization contributions
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
#include "dgroup.h"
      DIMENSION POLVEC(NBAS,NTOT,NZ),AMAT(NTOT,NTOT,NZ),
     &          BVEC(NTOT,NMOL,NZ)
C
      NREF = NTOT-NMOL
      IF(IPRINT.GE.5) THEN
        WRITE(6,*) '* NORMPOL: Overlaps..',NTOT,NMOL,NREF
        CALL PRQMAT(AMAT,NTOT,NTOT,NTOT,NTOT,NZ,
     &             IPQTOQ(1,0),LUPRI)
      ENDIF
      DO J = 1,NMOL
        FAC = D1/SQRT(AMAT(NREF+J,NREF+J,1))
        DO IZ = 1,NZ
          CALL DSCAL(NBAS,FAC,POLVEC(1,NREF+J,IZ),1)
        ENDDO
        BVEC(NREF+J,J,1)=D1/FAC
      ENDDO
      RETURN
      END
#ifdef ANALYZE_PROPERTY_GRADIENT
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck ANALYZE_PROPERTY_GRAD */
      SUBROUTINE ANALYZE_PROPERTY_GRAD(WORK,LWORK)
C=======================================================================
C     analyze property gradient
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      PARAMETER (NZ   = 4)
!     PARAMETER (NREF = 38) !se
      PARAMETER (NREF = 88) !po
      PARAMETER (NMOL = 2)
      PARAMETER (NTOT = NREF + NMOL)

#include "inforb.h"

      DIMENSION BVEC(NTOT,NMOL,NZ),
     &          PMAT(NTOT,NTOT,NZ),
     &          TMAT(NTOT,NTOT,NZ),
     &          QMAT(NTOT,NTOT,NZ),
     &          ZMAT(NTOT,NTOT),
     &          ZVEC(NTOT*NTOT),
     &          IND(NTOT*NTOT),
     &          IQ(4),
     &          GRADIENT(4),
     &          WORK(LWORK)

      CHARACTER PROP_NAME*16

#include "memint.h"

!     IHOMO = ...
!     ILUMO = ...

      DO I = 1,4
        IQ(I) = I
      ENDDO

!     read projection coefficients

      LU_PROJ = 13
      CALL OPNFIL(LU_PROJ,'DFPRJC','OLD','MAINNN')
      READ(LU_PROJ)
      READ(LU_PROJ)
      READ(LU_PROJ)
      READ(LU_PROJ) BVEC
      CLOSE(LU_PROJ,STATUS = 'KEEP')

!     read property matrix in fragment basis

      LU_PROP = 14
      CALL OPNFIL(LU_PROP,'PMAT','OLD','MAINNN')
 10   CONTINUE
      READ(LU_PROP,END = 20,ERR = 20) PROP_NAME
      READ(LU_PROP) PMAT

!     get property gradient

      CALL QTRANS('AOMO','S',D0,
     &            NTOT,NTOT,1,1,
     &            PMAT,NTOT,NTOT,NZ,IQ,
     &            GRADIENT,1,1,NZ,IQ,
     &            BVEC(1,1,1),NTOT,NMOL,NZ,IQ,
     &            BVEC(1,2,1),NTOT,NMOL,NZ,IQ,
     &            WORK(KFREE),LFREE,0)
      CALL DSCAL(4,D2,GRADIENT,1)

      WRITE(LUPRI,*) 'GRADANA: gradient for property ', PROP_NAME
      DO IZ = 1,NZ
        WRITE(LUPRI,*) GRADIENT(IZ)
      ENDDO

      CALL DZERO(QMAT,NTOT*NTOT*NZ)
      DO K = 1,NTOT
        DO L = 1,NTOT
          CALL QTRANS('AOMO','S',D0,
     &                1,1,1,1,
     &                PMAT(K,L,1),NTOT,NTOT,NZ,IQ,
     &                QMAT(K,L,1),NTOT,NTOT,NZ,IQ,
     &                BVEC(K,1,1),NTOT,NMOL,NZ,IQ,
     &                BVEC(L,2,1),NTOT,NMOL,NZ,IQ,
     &                WORK(KFREE),LFREE,0)
        ENDDO
      ENDDO
      CALL DSCAL(NTOT*NTOT*NZ,D2,QMAT,1)

      WRITE(LUPRI,*) 'GRADANA: compare'
      DO IZ = 1,NZ
        WRITE(LUPRI,*) DSUM(NTOT*NTOT,QMAT(1,1,IZ),1)
      ENDDO

!     symmetrize QMAT -> triangular TMAT

      CALL DZERO(TMAT,NTOT*NTOT*NZ)
      DO IZ = 1,NZ
        DO L = 1,NTOT
          DO K = 1,(L-1)
            TMAT(K,L,IZ) = QMAT(K,L,IZ) + QMAT(L,K,IZ)
          ENDDO
          TMAT(L,L,IZ) = QMAT(L,L,IZ)
        ENDDO
      ENDDO

!     square TMAT -> ZMAT or ZVEC

      CALL DZERO(ZMAT,NTOT*NTOT)
      CALL DZERO(ZVEC,NTOT*NTOT)
      DO K = 1,NTOT
        DO L = 1,NTOT
          DO IZ = 1,NZ
            ZMAT(K,L) = ZMAT(K,L) + TMAT(K,L,IZ)*TMAT(K,L,IZ)
          ENDDO
        ENDDO
      ENDDO
      DO K = 1,NTOT
        DO L = 1,NTOT
          ZVEC((K-1)*NTOT-1+L) = ZMAT(K,L)
        ENDDO
      ENDDO

      CALL INDEXX(NTOT*NTOT,ZVEC,IND)

      WRITE(LUPRI,*) 'GRADANA: ZMAT'
      CALL PRQMAT(ZMAT,NTOT,NTOT,NTOT,NTOT,1,IQ,LUPRI)
      WRITE(LUPRI,*) 'GRADANA: index'
      WRITE(LUPRI,'(2A6,6A16)') 'INDK', 'INDL', 'ZVEC',
     &                          'ZMAT',
     &                  'TMAT(IZ1)', 'TMAT(IZ2)','TMAT(IZ3)','TMAT(IZ4)'
      DO IDUMMY = 1,50
        I = NTOT*NTOT - IDUMMY + 1
        INDK = INT((IND(I)-0.1D0)/NTOT)
        INDL = IND(I) - INDK*NTOT + 1
        INDK = INDK + 1
        WRITE(LUPRI,'(2I6,6F16.8)') INDK, INDL, ZVEC(IND(I)),
     &                              ZMAT(INDK,INDL),
     &                             (TMAT(INDK,INDL,IZ),IZ=1,NZ)
!       IF((INDK.EQ.IHOMO).AND.(INDL.EQ.ILUMO)) THEN
!         WRITE(LUPRI,*) 'GRADANA: QMAT homo-lumo'
!         DO MYZ = 1,NZ
!           WRITE(LUPRI,*) MYZ, QMAT(INDK,INDL,MYZ)
!         ENDDO
!         WRITE(LUPRI,*) 'GRADANA: QMAT lumo-homo'
!         DO MYZ = 1,NZ
!           WRITE(LUPRI,*) MYZ, QMAT(INDL,INDK,MYZ)
!         ENDDO
!       ENDIF
      ENDDO

      GOTO 10
 20   CONTINUE

      RETURN
      END
#endif

      SUBROUTINE GET_AO_ATOM_MAP(IAOATOM)
C     Return an array with the index of the atom carrying a certain
C     AO. In case of SO combinations from multiple atoms give the index
C     of the first of these atoms. /UE 2008 (adapted from INDAOC)
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
      DIMENSION IAOATOM(NTBAS(0))
      NAORBC = 0
      JCENTD = 0
      DO JCENT = 1, NUCIND
        NDEGJ  = NUCDEG(JCENT)
        DO JDEG = 1, NDEGJ
          JCENTD = JCENTD + 1
          NAORB = 0
          DO ISHELL = 1,NLRGSH
            ICENT  = NCENT(ISHELL)
            LVAL   = NHKT(ISHELL)
            NDEGI  = NUCDEG(ICENT)
            DO IDEG = 1,NDEGI
              DO ICOMP = 1,KHKT(ISHELL)
                NAORB = NAORB + 1
                IF (ICENT.EQ.JCENT.AND.IDEG.EQ.JDEG) THEN
                   NAORBC = NAORBC + 1
                   IAOATOM(NAORBC) = JCENT
c                   print *,naorb,'<->',naorbc
                endif
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      NAORBL = NAORB
C
      JCENTD = 0
      DO JCENT = 1, NUCIND
        NDEGJ = NUCDEG(JCENT)
        DO JDEG = 1, NDEGJ
          JCENTD = JCENTD + 1
          NAORB  = NAORBL
          DO ISHELL = NLRGSH+1,KMAX
            ICENT  = NCENT(ISHELL)
            LVAL   = NHKT(ISHELL)
            NDEGI  = NUCDEG(ICENT)
            DO IDEG = 1,NDEGI
              DO ICOMP = 1,KHKT(ISHELL)
                NAORB = NAORB + 1
                IF (ICENT.EQ.JCENT.AND.IDEG.EQ.JDEG) THEN
                   NAORBC = NAORBC + 1
                   IAOATOM(NAORBC) = JCENT
                endif
              ENDDO
             ENDDO
          ENDDO
        ENDDO
      ENDDO
      END
c
C     Calculate a capped density matrix of the atoms
c     listed in icent.
c
      subroutine make_cap_dmat(icent,ncent,cmo,D,thres)
      implicit none
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      double precision cmo(*), D(norbt,norbt,nz),thres,ovsum
      integer icent, ncent, i, j, k
      dimension icent(ncent)
      print *,'making fragment projector'
      call make_fragment_projector(icent,ncent,D,cmo)
c     keep only occupied block of the projector
      do k=1,nz
         do i=nish(1)+1,norbt
            do j=1,norbt
               D(j,i,k) = 0
               D(i,j,k) = 0
            enddo
         enddo
      enddo
c     purify to get orbitals of the cap
      call purify(D,norbt,nz,thres)
      call matexport_double3('Dcap',D,norbt,norbt,norbt,norbt,nz)
      end

      subroutine make_fragment_projector(icent,ncent,P,cmo)
c     Make a projector on the basis of centers icent(1)
c     .. icent(ncent) (and their symmetry related centers). The
c     projector is created in MO basis (given by cmo) -ulfek
      implicit none
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      integer icent, ncent, aomap, i, j, ifind, n
      double precision P, tmp, cmo(*)
      double precision, allocatable :: Pao(:,:), Pso(:,:), vec(:,:),
     &     val(:)
      dimension P(norbt,norbt,nz), icent(ncent), aomap(ntbas(0))
c     Set up projector in AO basis
      call get_ao_atom_map(aomap)
      allocate( Pao(ntbas(0),ntbas(0)) )
      allocate( Pso(ntbas(0),ntbas(0)) )
      call dzero(Pao(1,1),ntbas(0)**2)
      n = 0
      do i=1,ntbas(0)
         if (ifind(aomap(i),icent,ncent).gt.0) then
            Pao(i,i) = 1
            n = n + 1
         endif
      enddo
      print *,'Marking',n,' orbitals as belonging to system A'
      call dzero(Pso(1,1),ntbas(0)**2)
c     go to so basis TODO: is irepdmp=0 correct?
      call mtaoso(Pao,Pso,ntbas(0),0,0)
      n = 0
      do i=1,ntbas(0)
         if (abs(Pso(i,i)).gt.1e-14) then
            n = n + 1
         endif
      enddo
c     to mo basis
      call qtrans90('AOMO','S',0.0D0,
     &     ntbas(0),ntbas(0),norbt,norbt,
     &     Pso,ntbas(0),ntbas(0),1,ipqtoq(1,0),
     &     P,norbt,norbt,nz,ipqtoq(1,0),
     &     cmo,ntbas(0),norbt,nz,ipqtoq(1,0),
     &     cmo,ntbas(0),norbt,nz,ipqtoq(1,0),0)
      deallocate( Pso )
      deallocate( Pao )
      call purify(P,norbt,nz,1.0D-8)
      end
C
      subroutine purify(P,n,mz,thres)
c     purify quaternion matrix P by setting eigenvalues > thres to 1 and
c     eigenvalues <= thres to 0.
      implicit none
#include "dgroup.h"
#include "dcbbas.h"
      double precision P(n,n,mz), thres
      double precision, allocatable :: val(:), vec(:,:,:)
      integer ierr, kveg, keig, nvec, n, mz, i, j, k
      allocate( vec(n,n,mz) )
      allocate( val(n) )
      call qdiag90(mz,n,P,n,n,val,1,vec,n,n,ierr)
      if (ierr.ne.0) then
         call quit('qdiag failed in purify()')
      endif
c     see how many nonzero eigenvalues we have
      nvec = 0
      print *,'In purify() - eigenvalues:'
      do i=1,n
         print *,val(i)
         if (val(i).gt.thres) then
            nvec = nvec + 1
         endif
      enddo
c     P = vec*vec' for vectors n-nvec+1 .. n
      CALL qgemm(n,n,nvec,1.0D0,
     &     'N','N',IPQTOQ(1,0),vec(1,n-nvec+1,1),n,n,mz,
     &     'H','N',IPQTOQ(1,0),vec(1,n-nvec+1,1),n,n,mz,
     &          0.0D0,IPQTOQ(1,0),P,n,n,mz)
      deallocate( vec )
      deallocate( val )
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prjaob */
      SUBROUTINE PRJAOB(KCSEL,KESEL,KBSEL,KBVEC,NVECS,KVEC,NSTR,IOPT,
     &              KRMC_FLG,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Carry out projection analysis in AO-basis (C1 atomic fragments)
C     Written by T.Saue May 14 2012
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbprj.h"
#include "nuclei.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION NVECS(2,0:2),KVEC(2,0:MAXREF),NSTR(2,0:2,0:MAXREF),
     &          NAO(MAXREF),WORK(*),IDIM(3)
      DIMENSION IOFF(2)
      CHARACTER TEXT*74
      LOGICAL FNDLAB
C
      KFRSAV = KFREE
      IF(IPRPRJ.GE.2) WRITE(6,'(A/)') '** Output from PRJAOB:'
      CALL MEMGET('INTE',KINF,NTBAS(0),WORK,KFREE,LFREE)
C
C     Find total number of reference orbitals
C     =======================================
C     ININD runs over symmetry independent centers
C     IREFS runs over symmetry dependent centers, exluding centers without basis functions
      NVECS(1,1)=0
      IREFS=1
      ININD=1
      DO ITYP = 1,NONTYP
        IF(.NOT.NOORBT(ININD)) THEN ! Skip centers with no orbitals
C.........Find the total number of basis functions for this atom type
          CALL LABCOUNT(MBAS,WORK(KINF),NTBAS(0),1,-1,ININD,1,-1)
C.........Read information from coefficient file
          CALL OPNFIL(LUCOEF,REFFIL(IREFS),'OLD','PRJAOB')        
          IF(FNDLAB('INFO    ',LUCOEF)) THEN         
            READ (LUCOEF,END=10,ERR=20) TEXT,NSYM,NZBUF,
     &          NSTR(1,2,IREFS),NSTR(1,1,IREFS),NAO(IREFS),DUM
C...............NPSH............NESH............NBAS.......TOTERG
          ELSE
            WRITE(6,*) 'PRJAOB: Error reading first record on ',
     &        REFFIL(IREFS)     
            CALL QUIT('PRJAOB: Error reading info')
          ENDIF
          CLOSE(LUCOEF,STATUS='KEEP')
          IF(NAO(IREFS).NE.MBAS) THEN
             WRITE(LUPRI,'(A,A,I3,/,A,A6,A,I6)') 
     &       'Atomic projection analysis: ',
     &       '  Wrong number of basis function for atomic type ', ITYP,
     &       '  Number read from                     ',
     &       REFFIL(IREFS),' : '  ,MBAS
             write(LUPRI,*) 'NAO(IREFS), MBAS:',NAO(IREFS), MBAS
             CALL QUIT
     &       ('PRJAOB: Wrong number of basis functions !')
          ENDIF
C.........Count the number of orbitals to read: 
C         NCATOM(1,0,IREFS) is the total number of orbitals per atom of this type
C         The number is multiplied with total number of atoms of this type and
C         added to NVECS(1,1)
          NSTR(1,0,IREFS)=0
          NDIM=NSTR(1,2,IREFS)+NSTR(1,1,IREFS)+1
          CALL MEMGET('INTE',KITMP,NDIM,WORK,KFREE,LFREE)
          CALL NUMLS1(VECREF(1,IREFS),IDUM,IDUM,
     &         -NSTR(1,2,IREFS),NSTR(1,1,IREFS),
     &         NSTR(1,0,IREFS),WORK(KITMP))
          CALL MEMREL('PRJAOB.ref',WORK,1,KITMP,KFREE,LFREE)
          CALL MEMGET('INTE',KVEC(1,IREFS),NSTR(1,0,IREFS),
     &         WORK,KFREE,LFREE)
          NDEP = 0
          DO J = 0,NONT(ITYP)-1
            NDEP = NDEP + NUCDEG(ININD+J)
          ENDDO
          NVECS(1,1) = NVECS(1,1) + NDEP*NSTR(1,0,IREFS)
          IREFS=IREFS+NDEP
        ENDIF  
        ININD = ININD + NONT(ITYP)
      ENDDO
C
C     Find number of molecular orbitals to analyze
C     ============================================
C
      NVECS(1,2) = 0
      DO IFRP = 1,NFSYM
        NSTR(IFRP,1,0) = NESH(IFRP)
        NSTR(IFRP,2,0) = NPSH(IFRP)
        CALL ORBNUM(VECPRJ(IFRP),IFRP,KVEC(IFRP,0),NSTR(1,0,0),
     &              WORK,KFREE,LFREE)
        NVECS(1,2) = NVECS(1,2) + NSTR(IFRP,0,0)
      ENDDO
C
C     Calculate dimension of arrays
C     =============================
C
      NVECS(1,0) = NVECS(1,1) + NVECS(1,2)
      CALL MEMGET('REAL',KESEL,NVECS(1,0),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KBSEL,NVECS(1,0),WORK,KFREE,LFREE)
      CALL IZERO(WORK(KBSEL),NVECS(1,1))
      NCDIM = NTBAS(0)*NVECS(1,0)*NZC1
      CALL MEMGET('REAL',KCSEL,NCDIM   ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KCSEL),NCDIM)
C
C     Select atomic coefficients
C     ==========================
      JVEC  = 1
      IREFS = 0
      ININD = 1
      IVOFF = 0
      DO ITYP = 1,NONTYP
        IF(.NOT.NOORBT(ININD)) THEN
          IREFS = IREFS + 1
C.........Read C1 coefficients and orbitalenergies
          CALL OPNFIL(LUCOEF,REFFIL(IREFS),'OLD','PRJAOB')        
          IF(FNDLAB('INFO    ',LUCOEF)) THEN
            READ (LUCOEF,END=10,ERR=20) TEXT,NSYM,NZBUF,
     &        NNEG,NPOS,NAO(IREFS),DUM
          ELSE
            WRITE(LUPRI,*) 'PRJAOB: error reading coefficients ',
     &                    REFFIL(IREFS) 
            CALL QUIT('PRJAOB: Error reading coefficients')
          ENDIF
          NTOT = NNEG+NPOS
          NQA  = NAO(IREFS)*NTOT*NZC1
          IDIM(1) = NNEG
          IDIM(2) = NPOS
          IDIM(3) = NAO(IREFS)
          CALL MEMGET('REAL',KQA,NQA ,WORK,KFREE,LFREE)
          CALL REACMO_coefficients(LUCOEF,WORK(KQA),NQA,1,NZC1,
     &           NTOT,NAO(IREFS),NNEG,NPOS,IDIM)         
          CALL MEMGET('REAL',KEA,NTOT,WORK,KFREE,LFREE)
          CALL REACMO_eigenvalues(LUCOEF,WORK(KEA),NTOT,1,NZC1,
     &           NTOT,NNEG,NPOS,IDIM)         
          CLOSE(LUCOEF,STATUS='KEEP')
C.........Make pointer array to selected reference orbitals
          NDIM = NTOT+1
          CALL MEMGET('INTE',KITMP,NDIM,WORK,KFREE,LFREE)
          CALL NUMLS1(VECREF(1,IREFS),WORK(KVEC(1,IREFS)),
     &         NSTR(1,0,IREFS),-NNEG,NPOS,NSTR(1,0,IREFS),WORK(KITMP))
          CALL MEMREL('PRJAOB.refsel',WORK,1,KITMP,KFREE,LFREE)
          CALL ORBCNT(WORK(KVEC(1,IREFS)),NSTR(1,0,IREFS),NNEG,NPOS,
     &               NSTR(1,2,IREFS),NSTR(1,1,IREFS))
C.......Extract coefficients and orbital energies (JVEC updated inside)
          CALL SELATOM(WORK(KCSEL),WORK(KESEL),NTBAS(0),NVECS(1,0),NZC1,
     &                 WORK(KQA),WORK(KEA),NAO(IREFS),NTOT,NNEG,NPOS,
     &                 ININD,NONT(ITYP),NUCDEG(ININD),
     &                 WORK(KVEC(1,IREFS)),NSTR(1,1,IREFS),
     &                 NSTR(1,2,IREFS),JVEC,WORK,KFREE,LFREE)
          CALL MEMREL('PRJAOB.atom',WORK,1,KQA,KFREE,LFREE)
          IF(IPRPRJ.GE.2) THEN
            WRITE(6,'(A,I5,A,3X,A6/A,2I6/A,I6)') 
     &       '  Atomic type ',ITYP,'  ',REFFIL(IREFS),
     &       '  NPSH/NESH       :',NSTR(1,2,IREFS),NSTR(1,1,IREFS),
     &       '  Basis set size  :',NAO(IREFS)
            WRITE(6,*) 'Selected coefficients:'
            KCOFF = KCSEL+NTBAS(0)*IVOFF
            NVDIM = NSTR(1,0,IREFS)
            CALL PRQMAT(WORK(KCOFF),NTBAS(0),NVDIM,NTBAS(0),
     &                  NVECS(1,0),NZC1,IQDEF,LUPRI)
            WRITE(6,*) 'Selected eigenvalues:'
            KEOFF = KESEL+IVOFF
            CALL OUTPUT(WORK(KEOFF),1,NVDIM,1,1,NVDIM,1,-1,LUPRI)
          ENDIF
C.........Duplicate information for other atoms of this type
          NDEP=0
          DO J = 0,NONT(ITYP)-1
            NDEP = NDEP + NUCDEG(ININD+J)
          ENDDO
          NN = IREFS
          WRITE(REFFIL(NN)(5:6),'(I2)') 1
          DO IDEP = 2,NDEP
            IREFS = IREFS + 1
            WRITE(REFFIL(IREFS),'(A4,I2)') REFFIL(NN)(1:4),IDEP
            NSTR(1,0,IREFS) = NSTR(1,0,NN)
            NSTR(1,1,IREFS) = NSTR(1,1,NN)
            NSTR(1,2,IREFS) = NSTR(1,2,NN)
            CALL MEMGET('INTE',KVEC(1,IREFS),NSTR(1,0,IREFS),
     &                WORK,KFREE,LFREE)
            CALL ICOPY(NSTR(1,0,IREFS),WORK(KVEC(1,NN)),1,
     &                                 WORK(KVEC(1,IREFS)),1)
          ENDDO
          IVOFF = IVOFF + NDEP*NSTR(1,0,IREFS)
        ENDIF 
        ININD = ININD+NONT(ITYP)
      ENDDO
C
C     Select molecular coefficients
C     =============================
C     Selected MOs etc. in full symmetry first goes into buffers
      NCDIM = 0
      DO IFRP = 1,NFSYM
        NCDIM = NCDIM + NFBAS(IFRP,0)*NSTR(IFRP,0,0)*NZ
      ENDDO
      CALL MEMGET('REAL',KCBUF,NCDIM     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEBUF,NVECS(1,2),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,NVECS(1,2),WORK,KFREE,LFREE)
C
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBE,NORBT   ,WORK,KFREE,LFREE)
C
      DO IFRP = 1,NFSYM
        IOFF(IFRP) = 0
      ENDDO
      CALL SELFRAG(WORK(KCBUF),NSTR(1,0,0),WORK(KEBUF),WORK(KIBUF),
     &             1,NSTR(1,0,0),'DFCOEF',.FALSE.,1,IOPT,
     &             KVEC(1,0),WORK(KCMO),WORK(KEIG),WORK(KIBE),
     &             IOFF,KRMC_FLG,WORK,KFREE,LFREE)
      IF(IPRPRJ.GE.2) THEN
        KCOFF = KCBUF
        DO IFRP = 1,NFSYM
        WRITE(6,*) 'Selected molecular coefficients in SO-basis',IFRP
        CALL PRQMAT(WORK(KCOFF),NFBAS(IFRP,0),NSTR(IFRP,0,0),
     &              NFBAS(IFRP,0),NSTR(IFRP,0,0),NZ,IPQTOQ(1,0),LUPRI)
        KCOFF = KCOFF + NFBAS(IFRP,0)*NSTR(IFRP,0,0)*NZ
        ENDDO
      ENDIF
      CALL MEMREL('PRJAOB.mol',WORK,1,KCMO,KFREE,LFREE)
C.....Then transform coefficients to AO-basis
      CALL C1COEF(WORK(KCBUF),WORK(KCSEL),NSTR(1,0,0),
     &            NVECS(1,1),NVECS(1,0),0)
      CALL DCOPY(NVECS(1,2),WORK(KEBUF),1,WORK(KESEL+NVECS(1,1)),1)
      CALL IICOPY(NVECS(1,2),WORK(KIBUF),1,1,
     &                       WORK(KBSEL),1,NVECS(1,1)+1)
      IF(IPRPRJ.GE.2) THEN
        WRITE(6,*) 'Selected molecular coefficients'
        KCOFF = KCSEL+NTBAS(0)*NVECS(1,1)
        NVDIM = NVECS(1,2)
        CALL PRQMAT(WORK(KCOFF),NTBAS(0),NVDIM,NTBAS(0),
     &                NVECS(1,0),NZC1,IQDEF,LUPRI)
          WRITE(6,*) 'Selected molecular eigenvalues:'
          KEOFF = KESEL+NVECS(1,1)
          CALL OUTPUT(WORK(KEOFF),1,NVDIM,1,1,NVDIM,1,-1,LUPRI)
      ENDIF
      IF(IPRPRJ.GE.2) THEN
        WRITE(6,*) 'All selected coefficients'
        CALL PRQMAT(WORK(KCSEL),NTBAS(0),NVECS(1,0),NTBAS(0),
     &                NVECS(1,0),NZC1,IQDEF,LUPRI)
        WRITE(6,*) 'All selected eigenvalues'
        CALL OUTPUT(WORK(KESEL),1,NVECS(1,0),1,1,NVECS(1,0),1,-1,LUPRI)
      ENDIF
      CALL MEMREL('PRJAOB.mol2',WORK,1,KCBUF,KFREE,LFREE)
      NBDIM = NVECS(1,0)*NVECS(1,2)*NZC1
      CALL MEMGET('REAL',KBVEC,NBDIM,WORK,KFREE,LFREE)
      RETURN
 10   CONTINUE
      CALL QUIT('PRJAOB: END OF FILE reading TEXT')
 20   CONTINUE
      CALL QUIT('PRJAOB: ERROR reading TEXT')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck selatom */
      SUBROUTINE SELATOM(CMO,EIG,NRM,NCM,NZ,
     &                   ACMO,AEIG,NRA,NCA,NNEG,NPOS,
     &                   INUC,NONT,NDEG,
     &                   IVEC,NEVEC,NPVEC,JVEC,
     &                   WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by T. Saue May 14 2012
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION CMO(NRM,NCM,NZ),EIG(*),ACMO(NRA,NCA,NZ),AEIG(*),
     &          NDEG(*),IVEC(*),WORK(*)
      CALL MEMGET('INTE',KMOL ,NRM,WORK,KFREE,LFREE)      
      CALL MEMGET('INTE',KATOM,NRA,WORK,KFREE,LFREE)
      NVEC=NEVEC+NPVEC
      ININD = INUC      
      DO ITYP = 1,NONT
        DO IDEG = 1,NDEG(ITYP)
C.........Create pointer mol->atom for rows
          CALL LABCOUNT(MBAS,WORK(KMOL),NRM,1,-1,ININD,IDEG,-1)
C.........Create reverse pointer atom->mol
          CALL REVINDEX(WORK(KMOL),NRM,WORK(KATOM),NRA)
C.........Extract coefficients
          CALL SELOWC(NZ,ACMO,NRA,NCA,NNEG,NPOS,
     &                CMO(1,JVEC,1),NRM,NCM,
     &                IVEC,WORK(KATOM),NPVEC,NEVEC)
C.........Extract orbital energies
          CALL SELOWE(AEIG,NNEG,NPOS,
     &                EIG(JVEC),IVEC,NPVEC,NEVEC)
          JVEC = JVEC + NVEC
        ENDDO
        ININD = ININD + 1
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjana_mod */
      SUBROUTINE PRJANA_MOD(CSEL,IPQSEL,NZC,AMAT,BVEC,IPQBVC,NZB,
     &                      SMAT,LRSMT,LCSMT,IPQSMT,NZS,
     &                      NBAS,NTOT,NREF,NMOL,PTOL,
     &                      IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Projection analysis for fermion ircop IFRP
C
C     INPUT:
C       CSEL(NBAS,NTOT,NZC) contains reference and molecular coeffients
C     OUTPUT:
C       AMAT(NTOT,NTOT,NZB) contains overlap matrix in fragment basis
C       BVEC(NTOT,NMOL) contains solution vector
C       CSEL(NBAS,NTOT,NZC) contains reference and polarization coefficients
C
C     Written by T.Saue June 5 2012
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,DM1=-1.0D0)
C
      DIMENSION CSEL(NBAS,NTOT,NZC),SMAT(LRSMT,LCSMT,NZS),
     &          AMAT(NTOT,NTOT,NZB),BVEC(NTOT,NMOL,NZB),
     &          WORK(*)
C
      KFRSAV = KFREE
      KREF   = 1
      KMOL   = KREF + NREF
C
C     Make projection vector  BVEC = <\psi^B_j1\psi^{MO}_i>
C
      CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NREF,NMOL,
     &            SMAT,LRSMT,LCSMT,NZS,IPQSMT,
     &            BVEC,NTOT ,NMOL ,NZB,IPQBVC,
     &            CSEL(KREF,1,1),NBAS,NTOT,NZC,IPQSEL,
     &            CSEL(KMOL,1,1),NBAS,NTOT,NZC,IPQSEL,
     &            WORK(KFREE),LFREE,IPRINT)
C
C     Make overlap matrix AMAT = <\psi^B_j1\psi^A_i>
C
      CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NREF,NMOL,
     &            SMAT,LRSMT,LCSMT,NZS,IPQSMT,
     &            AMAT,NTOT ,NTOT ,NZB,IPQBVC,
     &            CSEL(KREF,1,1),NBAS,NTOT,NZC,IPQSEL,
     &            CSEL(KREF,1,1),NBAS,NTOT,NZC,IPQSEL,
     &            WORK(KFREE),LFREE,IPRINT)
C
C     Print overlaps
C
C     ...write me...
C
C     Solve linear system by Cholesky decomposition
C     with full pivoting
C
      JOB = 1
      CALL MEMGET('REAL',KB1,NREF,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KB2,NREF,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIP,NREF,WORK,KFREE,LFREE)
      CALL QCHOLD(AMAT,NREF,NZB,NTOT,NTOT,WORK(KB1),
     &            PTOL,NEFF,JOB,WORK(KIP))
      CALL QCHOLS(AMAT,NREF,NEFF,NMOL,NZB,NTOT,NTOT,WORK(KB1),
     &            BVEC,NTOT,NMOL,BVEC,NTOT,NMOL,
     &            JOB,WORK(KIP),WORK(KB2))
      IF(NEFF.LT.NREF) THEN
        WRITE(LUPRI,'(A,A)') '* WARNING: ',
     &   ' linear dependencies detected in Cholesky decomposition'
        WRITE(LUPRI,'(11X,A,I5,A,I5)') 
     &   ' Reduced dimensionality: ',NREF,' --> ',NEFF
      ENDIF
      CALL MEMREL('PRJANA_MOD',WORK,1,KFRSAV,KFREE,LFREE)
C
C     Calculate polarization contribution;
C     overwrites molecular coefficients
C
      CALL DUNIT2(BVEC(NREF+1,1,1),NMOL,NTOT,NMOL,NZB)
      CALL QGEMM(NBAS,NMOL,NREF,DM1,
     &           'N','N',IPQSEL,CSEL(KREF,1,1),NBAS,NTOT,NZC,
     &           'N','N',IPQBVC,BVEC,NTOT,NMOL,NZB,
     &                D1,IPQSEL,CSEL(KMOL,1,1),NBAS,NTOT,NZC)
      IF(IPRINT.GE.5) THEN
        CALL HEADER('Polarization vectors',0)
      ENDIF      
C
C     Make overlap matrix AMAT = <\psi^B_j1\psi^A_i>
C     but now with polarization contributions included
C
      CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NTOT,NTOT,
     &           SMAT,LRSMT,LCSMT,NZS,IPQSMT,
     &           AMAT,NTOT,NTOT,NZB,IPQBVC,
     &           CSEL,NBAS,NTOT,NZC,IPQSEL,
     &           CSEL,NBAS,NTOT,NZC,IPQSEL,
     &           WORK(KFREE),LFREE,IPRINT)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck OVLFRAG */
      SUBROUTINE OVLFRAG(AMAT,NVECS,NONTYP,REFFIL,NCATOM,NATOMS,
     &                   CSEL,NBAS,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Starting from selected atomic coefficients in C1 symmetry,
C     generate overlap matrix in fragment basis
C
C     Written by T. Saue June 6 2012
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,DM1=-1.0D0)
#include "dcbham.h"
#include "dgroup.h"
      CHARACTER*6 REFFIL(NONTYP)
      DIMENSION AMAT(NVECS,NVECS,NZC1),NCATOM(0:2,NONTYP),
     &          NATOMS(NONTYP),CSEL(NBAS,NVECS,NZC1),WORK(*)
      KFRSAV = KFREE
      CALL MEMGET('REAL',KSMAT,NBAS*NBAS,WORK,KFREE,LFREE)
      CALL GET_OVERLAP_AO(WORK(KSMAT),SSMTRC,IPRINT,WORK,KFREE,LFREE)
      CALL QTRANS('AOMO','S',D0,NBAS,NBAS,NVECS,NVECS,
     &           WORK(KSMAT),NBAS,NBAS,1,1,
     &           AMAT,NVECS,NVECS,NZC1,IQDEF,
     &           CSEL,NBAS,NVECS,NZC1,IQDEF,
     &           CSEL,NBAS,NVECS,NZC1,IQDEF,
     &           WORK(KFREE),LFREE,IPRINT)
      IF(IPRINT.GE.4) THEN
        WRITE(LUPRI,*) '** Output from OVLFRAG **'
        KJJ = 1
        DO J = 1,NONTYP
          DO JJ = 1,1,NATOMS(J)
            KII = 1
            DO I = 1,J
              IF(I.EQ.J) THEN
                IIEND = JJ-1
              ELSE
                IIEND = NATOMS(I)
              ENDIF
              DO II = 1,IIEND
                WRITE(LUPRI,*) ' * Overlap ',REFFIL(I),II,REFFIL(J),JJ
                CALL PRQMAT(AMAT(KII,KJJ,1),NCATOM(0,I),NCATOM(0,J),
     &                      NVECS,NVECS,NZC1,IQDEF,LUPRI)
                KII = KII + NCATOM(0,I)
              ENDDO
            ENDDO
            KJJ = KJJ + NCATOM(0,J)
          ENDDO
        ENDDO
        WRITE(LUPRI,*) '* Total overlap matrix'
        CALL PRQMAT(AMAT,NVECS,NVECS,NVECS,NVECS,NZC1,IQDEF,LUPRI)
      ENDIF
      CALL MEMREL('OVLFRAG',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
! --- end of dirac/dirana.F ---

