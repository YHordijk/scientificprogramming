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
C MAIN AUTHOR : Ulf Ekstroem, 2005-2006.
C ORIGINAL PUBLICATION OF THIS WORK :
C           U.Ekstroem, P. Norman, V. Carraveta,
C           Phys.Rev.A 73(2006)022501
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE STEX_INPUT(WORD,RESET,WORK,LWORK)
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbstex.h"
      PARAMETER (NTABLE = 6)
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,
     &          PNAME(STEX_NR_OP)*16
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.HOLES ','.PRINT ','.BATCH ','.CUBE  ',
     &            '.CUTOFF','.SCREEN'/
      DATA SET/.FALSE./ 
C
#include "ibtfun.h"
C
      WRITE (LUPRI,*) 'STEX input section'
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF

      DOSTEX = .TRUE.
c     cofactor if separate ref state only, checked in STEX()
      STEX_COFACTOR = .TRUE.
C     Should we neglect all 2-electron integrals?
      STEX_NO2E = .FALSE.
C     Should we keep only the diagonal of the Hessian?
      STEX_DIAG = .FALSE.
C     Couple the excitations from different holes
      STEX_NOCOUPLING = .FALSE.
C     Nr of simultanous Fock matrices to build.
C     What is the optimal fock matrix calculation order?
C     speed is cut in half if the matrices are not
C     (anti)hermitian, so batch size 4 is maybe good.
      STEX_BATCHSIZE = 4
C     Don't write cube files by default.
      STEX_NCUBE = 0
C     Set default cutoff in eV
      STEX_CUTOFFEV = 100.0D0
      STEX_CUTOFF = STEX_CUTOFFEV/XTEV
C     Screening (-1 means not set here)
      STEX_SCRFCK = -1.0D0
C     Print level
      STEX_PRINT=0

C     Add diplen and dipvel to DCBPRP
      PNAME(1) = 'XDIPLEN '
      PNAME(2) = 'YDIPLEN '
      PNAME(3) = 'ZDIPLEN '
      PNAME(4) = 'XDIPVEL '
      PNAME(5) = 'YDIPVEL '
      PNAME(6) = 'ZDIPVEL '
      NSTXOP = 0
      DO I=1,STEX_NR_OP
         CALL XPRIND(PNAME(I),1,1,1.0D0,PNAME(I),INDXPR,STEX_OP_SYM(I),
     &        STEX_OP_TIM(I),0)
         CALL OP1IND('STXDOP',IND1OP,STEX_OP_IND,NSTXOP,INDXPR
     $        ,STEX_NR_OP)
      END DO
C
C     Process input
C     =========================
C
      NEWDEF = (WORD .EQ. '*STEX  ')
      ICHANG = 0
      INPERR = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
  110       CONTINUE
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1, 2, 3, 4, 5, 6), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in STXINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in STXINP.')
C     .HOLES
C     Expects format like (NFSYM=2)
C     2    - nr of holes in fs 1.
C     5 9  - holes in fs 1
C     1    - nr of holes in fs 2.
C     1    - holes in fs 2
C
    1          CONTINUE
               STEX_NHOLES = 0
               DO IFSYM=1,NFSYM
                  READ (LUCMD,*) NHOLES
                  IF (NHOLES+STEX_NHOLES.GT.STEX_MAX_HOLES) THEN
                     WRITE (LUPRI,*)
     $                    'Error in STEX_INPUT(): Nr holes too large.'
                     WRITE (LUPRI,*) 'STEX_MAX_HOLES =',STEX_MAX_HOLES
                     CALL QUIT('STEX input error')
                  ENDIF
                  READ (LUCMD,*)
     &                 (STEX_HOLES(II+STEX_NHOLES),II=1,NHOLES)
                  DO II=1,NHOLES
                     STEX_IFSYM(II+STEX_NHOLES) = IFSYM
                  ENDDO
                  STEX_NHOLES = STEX_NHOLES + NHOLES
               ENDDO

               GOTO 100
C     .PRINT
    2          CONTINUE
               READ (LUCMD,*) STEX_PRINT
               GOTO 100
C     .BATCH SIZE
    3          CONTINUE
               READ (LUCMD,*) STEX_BATCHSIZE
               GOTO 100
C     TODO: move .CUBE to analyze
C     .CUBE
    4          CONTINUE
               READ (LUCMD,*) STEX_NCUBE
               GOTO 100
C     .CUTOFF
    5          CONTINUE
               READ (LUCMD,*) STEX_CUTOFFEV
               STEX_CUTOFF = STEX_CUTOFFEV/XTEV
               GOTO 100
C     .SCREEN
    6          CONTINUE
               READ (LUCMD,*) STEX_SCRFCK
               GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in STXINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in STXINP.')
            END IF
      END IF
  300 CONTINUE
  999 CONTINUE           
      WRITE (LUPRI,*) 'STEX CUTOFF:',STEX_CUTOFFEV
      WRITE (LUPRI,*) 'STEX SCRFCK:',STEX_SCRFCK
      WRITE (LUPRI,*) 'STEX NHOLES:',STEX_NHOLES
      WRITE (LUPRI,*) 'Hole orbital   Fermion symmetry'
      DO I=1,STEX_NHOLES
         WRITE (LUPRI,*) STEX_HOLES(I),STEX_IFSYM(I)
      ENDDO
      IF (INPERR.GT.0) CALL QUIT('Input error in *STEX')      
      WRITE (LUPRI,*) 'End of STEX input section'
      END

      SUBROUTINE STEX(WORK,LWORK)
C
C     At this point we expect a converged groundstate DFCOEF in
C     the working directory, plus an optional DFCOEF.ION for the
C     core-relaxed determinant.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcbdhf.h"
#include "dcbstex.h"
C#include "dcbrho.h"
      DOUBLE COMPLEX REFGSOVL
      DIMENSION WORK(LWORK)
      LOGICAL LION,LTOBE
#include "memint.h"
      CALL QENTER('STEX')
      WRITE(LUPRI,*) 'Start of STEX'
C
C     Sanity check. Unfortunately we cannot do this in
C     the input section since the environment is not complete
C     there.
C
      DO IFSYM=1,NFSYM
         IF (NASH(IFSYM).GT.0) THEN
            CALL QUIT
     $         ('ERROR: Active shells not supported in STEX, quitting.')
         ENDIF
      ENDDO
      IF (NZ.LT.4) THEN
         CALL QUIT('NZ < 4 not yet supported in STEX, run in C1')
      ENDIF
      IF (NFSYM.GT.1) THEN
         CALL QUIT('NFSYM > 1 not yet supported in STEX, run in C1')
      ENDIF

C     Save ground state energy
      GSERG = DHFERG

C     Check for existence of DFCOEF (should always be present)
      INQUIRE(FILE='DFCOEF',EXIST=LTOBE)
      IF(.NOT.LTOBE) THEN
         WRITE(LUPRI,*) 'Missing DFCOEF'
         CALL QUIT('File DFCOEF missing')
      ENDIF
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KBREP,NORBT,WORK,KFREE,LFREE)
      LUBUF=22
      IOPT = 2+4+8
C     See if a separate reference state will be used and read reference orbitals and eigenvalues
      INQUIRE(FILE='DFCOEF.ION',EXIST=LION)
      IF (LION) THEN
         WRITE(LUPRI,*) 'DFCOEF.ION found, using it as reference state'
         WRITE(LUPRI,*) 'for STEX calculation'
         STEX_COFACTOR = .TRUE.
         CALL OPNFIL(LUBUF,'DFCOEF.ION','OLD','STEX  ')
         CALL REACMO(LUBUF,'DFCOEF.ION',WORK(KCMO),WORK(KEIG),
     &        WORK(KBREP),TOTERG,IOPT)
      ELSE
         WRITE(LUPRI,*) 'No separate reference state (DFCOEF.ION) found'
         WRITE(LUPRI,*) 'using groundstate as STEX reference state.'
         STEX_COFACTOR = .FALSE.
         CALL OPNFIL(LUBUF,'DFCOEF','OLD','STEX  ')
         CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),WORK(KEIG),
     &        WORK(KBREP),TOTERG,IOPT)
      ENDIF
      CLOSE(LUBUF,STATUS='KEEP')

C     Set up cutoff and indices
      STEX_NVIR(2) = 0
      print *,'Virtual cutoff',STEX_CUTOFF,' Hartree'
      DO ISYM=1,NFSYM
         JEIG = KEIG + IORB(ISYM) + NPSH(ISYM) + NISH(ISYM) - 1
         I = 1
         DO WHILE (I.LE.NSSH(ISYM).AND.
     &        (WORK(JEIG+I).LT.STEX_CUTOFF
     &        .OR.WORK(JEIG+I)-WORK(JEIG+I-1).LT.0.05D0)) 
            I = I + 1
         ENDDO
         STEX_NVIR(ISYM) = I - 1
      ENDDO
      WRITE (LUPRI,*) 'STEX Nr virtuals:',(STEX_NVIR(I),I=1,NFSYM)
      STEX_NVIRT = STEX_NVIR(1) + STEX_NVIR(2)
      IF (STEX_NVIRT .LE. 0) THEN
         WRITE (LUPRI,*) 'STEX aborting, no virtuals below cutoff'
         RETURN
      END IF
      STEX_N2VIRX = STEX_NVIRT*STEX_NVIRT
      STEX_N2VIRXQ = NZ*STEX_N2VIRX
      STEX_I2VIRX(1,1) = 0
      STEX_I2VIRX(2,1) = STEX_NVIR(1)
      STEX_I2VIRX(1,2) = STEX_NVIR(1)*STEX_NVIRT
      STEX_I2VIRX(2,2) = STEX_I2VIRX(1,2) + STEX_NVIR(1)
c      WRITE(LUPRI,*) 'STEX: Dimension of A matrix is',STEX_NVIRT,'*'
c     $     ,STEX_NHOLES,'*4'
      WRITE(LUPRI,*) 'STEX holes in orbitals (total orbital index)'
      DO II=1,STEX_NHOLES
         CALL STEX_ORB_INDEX(IORBT,STEX_HOLES(II),STEX_IFSYM(II),2)
         STEX_HOLES(II) = IORBT
         WRITE(LUPRI,*) STEX_HOLES(II)
      ENDDO

      STEX_NPARAM=STEX_NVIRT*STEX_NHOLES*4
      WRITE(LUPRI,*) 'STEX: Nr of variational parameters is ',
     &     STEX_NPARAM
      NHOLECMB = (STEX_NHOLES*(STEX_NHOLES+1))/2

      IF (LFREE.LT.2*STEX_NPARAM**2) THEN
         WRITE (LUPRI,*)
     &        'Error: Not enough memory for STEX Hamiltonian and'
     &        ,' eigenvectors. '
         WRITE (LUPRI) 'Needed words (approx): ',
     &        2*STEX_NPARAM**2,', available:',LFREE
         CALL QUIT('Out of memory in STEX')
      ENDIF

C     Allocate memory for STEX_TM. (transition moments)
C     property gradients are complex, so *2 size of double.
      CALL MEMGET('REAL',KPRP,STEX_NPARAM*STEX_NR_OP*2,WORK,KFREE,LFREE)

      IF (STEX_COFACTOR) THEN         
         ERGION = TOTERG
         CALL MEMGET('REAL',KCMOGS,NCMOTQ,WORK,KFREE,LFREE)
         CALL OPNFIL(LUBUF,'DFCOEF','OLD','STEX  ')
         IOPT = 2+8
         CALL REACMO(LUBUF,'DFCOEF',WORK(KCMOGS),DUM,
     &        WORK(KBREP),TOTERG,IOPT)
         CALL COF_PRPGRAD(WORK(KPRP),WORK(KCMOGS),WORK(KBREP)
     $        ,WORK(KCMO),REFGSOVL,WORK(KFREE),LFREE)
      ELSE
         ERGION = 0.0D0
         REFGSOVL = 1.0D0
C     Calculate the property gradient the normal way
         CALL FROZEN_PRPGRAD(WORK(KPRP),WORK(KCMO),WORK(KBREP)
     $        ,WORK(KFREE),LFREE)
      ENDIF
C       IF (STEX_NCUBE.GT.0) THEN
C          WRITE (LUPRI,*) 'Writing cube files'
C          CALL FLUSH(LUPRI)
C C     Write cube files for the hole and the first STEX_NCUBE virtuals
C          ILSCUB = 3
C          IPRRHO = 0
C          NCUBE(1) = 80
C          NCUBE(2) = 80
C          NCUBE(3) = 80
C          CUBADJ(1) = 4.0D0
C          CUBADJ(2) = 8.0D0
C          DO II=1,NORB(1)
C             DO JJ=1,STEX_NHOLES
C                IF (II.EQ.STEX_HOLES(JJ)) THEN
C                   CALL GT3RHO(II,WORK(KFREE),LFREE)
C                ENDIF
C             ENDDO
C             IVIRT = II-NPSH(1)-NISH(1)
C             IF (IVIRT.GE.1.AND.IVIRT.LE.STEX_NCUBE) THEN
C                CALL GT3RHO(II,WORK(KFREE),LFREE)
C             ENDIF
C          ENDDO
C       ENDIF

C     Construct hole densities & calculate integrals
C     Write integrals to disk (to save memory)
      WRITE (LUPRI,*) 'STEX needs',NHOLECMB*4,' Fock matrices',
     &     ' that will be computed in batches of at most',STEX_BATCHSIZE
      STEX_LUINT = 22
      OPEN(STEX_LUINT,FILE='STEXINT',FORM='UNFORMATTED',
     +     ACCESS='DIRECT',RECL=8*STEX_N2VIRXQ,STATUS='UNKNOWN')
      CALL STEX_CALC_INTEGRALS(WORK(KCMO),WORK,KFREE,LFREE)
      CALL FLUSH(LUPRI)
      IF (STEX_COFACTOR) THEN
C     If we use a separate reference state the 2e Fock matrices
C     on disk have to be re-evaluated. This is done automatically
C     if we delete the old ones.
         I = SYSTEM('rm -f DFFCK2')
         I = SYSTEM('rm -f DFFOCK')
      ENDIF
      CALL MEMGET('REAL',KFOCKREF,N2BBASXQ,WORK,KFREE,LFREE)
      CALL FLUSH(LUPRI)
C     Add the reference state Fock matrix to the integrals
      WRITE(LUPRI,*) 'Calculating reference Fock matrix'
C     Calculate Fock matrix and energy for the reference determinant
      CALL STEX_GETFCK(WORK(KCMO),WORK(KFOCKREF),REFERG,
     &     WORK(KFREE),LFREE)
      WRITE (LUPRI,*)
      WRITE (LUPRI,*) '--------------- STEX Final Output --------------'
      WRITE (LUPRI,*)
      WRITE (LUPRI,*) 'STEX reference state energy (au)        :',REFERG
      IF (STEX_COFACTOR) THEN
         WRITE (LUPRI,*) 'Ground state energy (au)                :'
     $        ,GSERG
         WRITE (LUPRI,*) 'Ground/reference state overlap (abs val):',
     &     ABS(REFGSOVL)
         WRITE (LUPRI,*) 'DSCF Ionization potential (eV)          :',
     &        (ERGION-GSERG)*27.211383
      ENDIF
      WRITE (LUPRI,*)
      CALL STEX_MOD_INTEGRALS(WORK(KFOCKREF),WORK(KFREE),LFREE)
c     Build A matrix
      CALL MEMGET('REAL',KAMAT,STEX_NPARAM**2,WORK,KFREE,LFREE)
      CALL STEX_BUILD_A(WORK(KAMAT) ,WORK(KFREE),LFREE)
      CLOSE(STEX_LUINT,STATUS='DELETE')
      CALL MEMGET('REAL',KVEC,STEX_NPARAM**2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVAL,STEX_NPARAM,WORK,KFREE,LFREE)
c      WRITE(LUPRI,*) 'Diagonalizing the STEX Hamiltonian'
      CALL FLUSH(LUPRI)
      CALL QDIAG(1,STEX_NPARAM,WORK(KAMAT),STEX_NPARAM,STEX_NPARAM
     $     ,WORK(KVAL),1,WORK(KVEC),STEX_NPARAM,STEX_NPARAM,WORK(KFREE)
     $     ,LFREE,IERR)
      IF (IERR.NE.0) THEN
         PRINT *,'QDIAG IERR=',IERR
         CALL QUIT('QDIAG failed')
      ENDIF
      CALL STEX_ANALYZE(STEX_NPARAM,WORK(KVEC),WORK(KVAL),WORK(KPRP),
     &     GSERG,REFERG,REFGSOVL,WORK,LWORK)
      PRINT *,'STEX End'
      CALL QEXIT('STEX')
      END

      SUBROUTINE STEX_GETFCK(CMO,FMAT,ENERGY,WORK,LWORK)
C     Calculate the Fock matrix and total energy of the CMO
C     determinant (only inactive part). Assumes 1-e part is
C     on disk in DFFCK1
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcbdhf.h"
#include "cbihr2.h"
#include "dcbstex.h"
      DIMENSION WORK(LWORK),FMAT(N2BBASXQ),CMO(NCMOTQ)
#include "memint.h"

      CALL MEMGET('REAL',KF1,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KF2,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDEN,N2BBASXQ,WORK,KFREE,LFREE)
C     Read one electron fock matrix from disk
      CALL OPNFIL(LUFCK1,'DFFCK1','OLD','STEXGF')
      CALL REAFCK(LUFCK1,WORK(KF1),.FALSE.,1)
      CLOSE(LUFCK1,STATUS='KEEP')
C     density
      CALL GENDEN(WORK(KDEN),CMO,1,0)
C     Calculate 2e part
      ISYMOP(1) = 1
      IFCKOP(1) = 1
      IHRMOP(1) = 1
      IPRINT = 0
      CALL STEX_GET_NPOS(NPOS)
      CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
      SCRFCK_SAVE = SCRFCK
      print *,'SCRFCK:',SCRFCK
      IF (STEX_SCRFCK.GE.0.0D0) THEN
         SCRFCK = STEX_SCRFCK
         print *,'SCRFCK:',SCRFCK
      ENDIF
      CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,WORK(KF2),WORK(KDEN),1,
     &     WORK(KPOS),INTFLG,IPRINT,WORK(KFREE),LFREE)
      SCRFCK = SCRFCK_SAVE 
      print *,'SCRFCK:',SCRFCK
C     Calculate energy
      CALL ERGCAL(WORK(KF1),WORK(KF2),WORK(KDEN),WORK(KFREE),LFREE)
C     Add 2e part to 1e part
      DO II=0,N2BBASXQ-1
         WORK(KF1+II) = WORK(KF1+II) + WORK(KF2+II)
      ENDDO
C     Go to MO basis
      CALL QRFAOMO(FMAT,WORK(KF1),1,1,1,0,CMO,WORK(KFREE),LFREE)
      ENERGY = DHFERG
      END

      SUBROUTINE QIMUL(A,NDIM,NZ,FAC)
c     A = FAC*i*A
      IMPLICIT NONE
      INTEGER I,J,NDIM,NZ
      DOUBLE PRECISION A,SWAP,FAC
      DIMENSION A(NDIM,NDIM,NZ)
      DO I=1,NDIM
         DO J=1,NDIM
            SWAP = A(I,J,1)
            A(I,J,1) = -FAC*A(I,J,2)
            A(I,J,2) = FAC*SWAP
         ENDDO
      ENDDO
      DO I=1,NDIM
         DO J=1,NDIM
            SWAP = A(I,J,3)
            A(I,J,3) = -FAC*A(I,J,4)
            A(I,J,4) = FAC*SWAP
         ENDDO
      ENDDO
      END

      SUBROUTINE STEX_ORB_INDEX(IORBT, IORBI, IFSYM, IORBTYP)
C     Give the "total" orbital index IORBT from fermion symmetry of
C     the orbital (IFSYM = 1,2), and orbital type:
C     IORBTYP
C     1 Positronic
C     2 Inactive
C     3 Active
C     4 Virtual
C     IORBI is the number of the orbital "in its symmetry and type"
C     Example: CALL STEX_ORB_INDEX(IORBT,3,1,4) will put the index
C     of the third gerade virtual orbital in IORBT.
#include "implicit.h"
#include "dgroup.h"
#include "dcborb.h"
      IF (IFSYM.GT.NFSYM.OR.IORBTYP.LT.1.OR.IORBTYP.GT.4.OR.IORBI.LT.1)
     $     GOTO 666
      IORBT = IORBI
      IF (IORBTYP.GT.1) IORBT = IORBT + NPSH(IFSYM)
      IF (IORBTYP.GT.2) IORBT = IORBT + NISH(IFSYM)
      IF (IORBTYP.GT.3) IORBT = IORBT + NASH(IFSYM)
      IF (IFSYM.EQ.2) IORBT = IORBT + NORB(IFSYM)
      IF (IORBTYP.EQ.1.AND.IORBI.GT.NPSH(IFSYM)) GOTO 666
      IF (IORBTYP.EQ.2.AND.IORBI.GT.NISH(IFSYM)) GOTO 666
      IF (IORBTYP.EQ.3.AND.IORBI.GT.NASH(IFSYM)) GOTO 666
      IF (IORBTYP.EQ.4.AND.IORBI.GT.NSSH(IFSYM)) GOTO 666
      RETURN
 666  CONTINUE
         PRINT *,'STEX_ORB_INDEX: Invalid orbital: ',IORBI,IFSYM,IORBTYP
         CALL QUIT('Error in STEX_ORB_INDEX')
      END

      SUBROUTINE STEX_ORB_CLASSIFY(IORBT, IORBI, IFSYM, IORBTYP)
C     Given a "total" orbital index IORBT determine fermion symmetry of
C     the orbital (IFSYM = 1,2), and orbital type:
C     IORBTYP
C     1 Positronic
C     2 Inactive
C     3 Active
C     4 Virtual
C     Also determine the index within the orbitals symmetry and class (IORBI).
#include "implicit.h"
#include "dcborb.h"
      IF (IORBT.LE.NORB(1)) THEN
         IFSYM = 1
         IOFF = 0
         IORBI = IORBT
      ELSE
         IFSYM = 2
         IOFF = NORB(1)
         IORBI = IORBT - NORB(1)
      ENDIF
      IF (IORBT.LE.IOFF+NPSH(IFSYM)) THEN
         IORBTYP = 1
      ELSE
         IORBI = IORBI - NPSH(IFSYM)
         IF (IORBT.LE.IOFF+NPSH(IFSYM)+NISH(IFSYM)) THEN
            IORBTYP = 2
         ELSE
            IORBI = IORBI - NISH(IFSYM)
            IF (IORBT.LE.IOFF+NPSH(IFSYM)+NISH(IFSYM)+NASH(IFSYM))
     &           THEN
               IORBTYP = 3
            ELSE
               IORBI = IORBI - NASH(IFSYM)
               IORBTYP = 4
            ENDIF
         ENDIF
      ENDIF
      END


      SUBROUTINE STEX_PAIRDENS(CMO,IIORB,JJORB,DAO,IJK,WORK,LWORK)
C
C
C     Calculate psi_i Q psi_j^\dag, where Q is a quaternion phase {1,i,j,k}
C     IJK determines the quaternion phase, leading to
C     IJK=1 -> -0.5*(  |i><j| + |i'><j'| )
C         2 -> -0.5*(  |i><j| - |i'><j'| )
C         3 -> -0.5*(  |i><j'| - |i'><j| )
C         4 -> -0.5*(  |i><j'| + |i'><j| )
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
      DIMENSION WORK(LWORK),CMO(NTBAS(0),NORBT,NZ),DAO(N2BBASXQ)
      DIMENSION CORB1(NTBAS(0),NZ),CORB2(NTBAS(0),NZ)
      CHARACTER OP
#include "memint.h"
C     FIXME: NZ < 4 and IFSYM = 2 not handled.
c     maybe unpack, run and pack?

      DO IZ=1,NZ
         DO II=1,NTBAS(0)
            CORB1(II,IZ) = CMO(II,IIORB,IZ)
         ENDDO
      ENDDO
c     insert the quaternion phase according to ijk
      IF (IJK.LE.2) THEN
         DO IZ=1,NZ
            DO II=1,NTBAS(0)
               CORB2(II,IZ) = CMO(II,JJORB,IZ)
            ENDDO
         ENDDO
      ELSE
c     Multiply with -j from the right
         DO IZ=1,2
            DO II=1,NTBAS(0)
               CORB2(II,IZ) = CMO(II,JJORB,IZ+2)
            ENDDO
         ENDDO
         DO IZ=3,4
            DO II=1,NTBAS(0)
               CORB2(II,IZ) = -CMO(II,JJORB,IZ-2)
            ENDDO
         ENDDO
      ENDIF

      IF (IJK.EQ.2.OR.IJK.EQ.4) THEN
         OP = 'I'
      ELSE
         OP = 'N'
      ENDIF

      CALL QGEMM(NTBAS(0),NTBAS(0),1,-0.5D0,
     &     'N',OP,IPQTOQ(1,0),CORB1,NTBAS(0),1,NZ,
     &     'H','N',IPQTOQ(1,0),CORB2,NTBAS(0),1,NZ,
     &     0.0D0,IPQTOQ(1,0),DAO,NTBAS(0),NTBAS(0),NZ)
      END


      SUBROUTINE STEX_GET_NPOS(NPOS)
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "cbihr2.h"
#include "blocks.h"
#include "dcbfir.h"
      integer npos
      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. ,.TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      END

      SUBROUTINE STEX_CALC_INTEGRALS(CMO,WORK,KFREE,LFREE)
C
C     Calculate the needed 2e integrals and transform them to
C     mo basis. put the result on unit STEX_LUINT in mo basis
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcbstex.h"
#include "dcbdhf.h"
c#include "cbihr2.h"
      DIMENSION WORK(*),CMO(NCMOTQ),
     &     ISYMO(STEX_BATCHSIZE),
     &     IHRMO(STEX_BATCHSIZE),
     &     IFCKO(STEX_BATCHSIZE)
      CALL QENTER('STEX_CALC_INTEGRALS')
      CALL MEMCHK('INT0',WORK,1)
C     NPOS is for parallel stuff
      CALL STEX_GET_NPOS(NPOS)
      CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
      IPRINT=1
      ITWOFCKPRINT = 0
      CALL MEMGET('REAL',KDEN,STEX_BATCHSIZE*N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFMAT,STEX_BATCHSIZE*N2BBASXQ,WORK,KFREE,LFREE)
      WRITE (LUPRI,*) 'Calculating STEX 2e integrals'
      CALL FLUSH(LUPRI)
      IDMAT = 0
      KF = KFMAT
      ITOT = 0
      IBATCH = 1
c     put the density in the FMAT matrix
      DO II=1,STEX_NHOLES
         DO JJ=II,STEX_NHOLES
            DO IJK=1,4
               IDMAT = IDMAT + 1
               CALL STEX_PAIRDENS(CMO,STEX_HOLES(II),STEX_HOLES(JJ),
     &              WORK(KF),IJK,
     &              WORK(KFREE),LFREE)
               IF (IJK.EQ.2.OR.IJK.EQ.4) THEN
C     Multiply in an imaginary i to get a trev matrix.
                  CALL QIMUL(WORK(KF),NTBAS(0),NZ,1.0D0)
               ENDIF
C     Set properties of the corresponding Fock matrices
               ISYMO(IDMAT) = 1
               IF (II.EQ.JJ) THEN
                  IF (IJK.EQ.1) THEN
                     IHRMO(IDMAT) = 1
                  ELSE
                     IHRMO(IDMAT) = -1
                  ENDIF
               ELSE
                  IHRMO(IDMAT) = 0
               ENDIF
               IF (II.EQ.JJ) THEN
                  IF (IJK.EQ.1) THEN
                     IFCKO(IDMAT) = 1
                  ELSE
C     Coulomb term vanishes for i = j
                     IFCKO(IDMAT) = 3
                  ENDIF
               ELSE
                  IF (IJK.LE.2) THEN
                     IFCKO(IDMAT) = 1
                  ELSE
c     hmm..
                     IFCKO(IDMAT) = 1
                  ENDIF
               ENDIF
               KF = KF + N2BBASXQ
C     Calculate Fock matrices when the batch is full, or at the end of
C     the run.
               IF (IDMAT.GE.STEX_BATCHSIZE.OR.(II.EQ.STEX_NHOLES.AND.
     &              JJ.EQ.STEX_NHOLES.AND.IJK.EQ.4)) THEN
c     put the AO integrals in WORK(KDEN)
                  WRITE (LUPRI,*) 'Integral batch, i, j, herm. ',
     &                 IBATCH,II,JJ,IHRMO(IDMAT)
                  CALL FLUSH(LUPRI)
c                  SCRFCK_SAVE = SCRFCK
c                  IF (STEX_SCRFCK.GE.0.0D0) THEN
c                     SCRFCK = STEX_SCRFCK
c                  ENDIF
                  CALL TWOFCK(ISYMO,IHRMO,IFCKO,WORK(KDEN),WORK(KFMAT)
     $                 ,IDMAT,WORK(KPOS),INTFLG,ITWOFCKPRINT,
     &                 WORK(KFREE),LFREE)
c                  SCRFCK = SCRFCK_SAVE 
c     go to mo basis (put mo integrals on disk)
                  ISYM = 1
                  IFER = 1
                  ITIM = 1
                  DO KK=0,IDMAT-1
                     CALL STEXFAOMO(WORK(KFMAT),
     &                    WORK(KDEN+KK*N2BBASXQ),
     &                    ISYM,ITIM,IFER,0,CMO,WORK(KFREE),LFREE)
                     IMOD = MOD(KK+ITOT,2)
                     IF (IMOD.EQ.1) THEN
                        CALL QIMUL(WORK(KFMAT),
     &                       STEX_NVIRT,NZ,-1.0D0)
                     ENDIF
                     IREC = ITOT+KK+1
c     TODO: optimize write below (one big block is better)
                     WRITE (STEX_LUINT, REC=IREC)
     &                    (WORK(KFMAT+IJ),IJ=0,STEX_N2VIRXQ-1)
                  ENDDO
                  ITOT = ITOT + IDMAT
                  IDMAT = 0
                  KF = KFMAT
                  IBATCH = IBATCH + 1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      WRITE (LUPRI,*) 'Integrals done'
      CALL FLUSH(LUPRI)
      CALL QEXIT('STEX_CALC_INTEGRALS')
      END

      SUBROUTINE STEX_MOD_INTEGRALS(FMAT,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcbstex.h"
      DIMENSION FMAT(NORBT,NORBT,NZ), WORK(LWORK)
#include "memint.h"
      CALL MEMGET('REAL',KTMP,STEX_N2VIRXQ,WORK,KFREE,LFREE)
      N = 1
      IVIROFF = NPSHT + NISHT
      DO II=1,STEX_NHOLES
         DO JJ=II,STEX_NHOLES
            DO KK=0,3
C     Read integrals to modify
               READ (STEX_LUINT,REC=N)
     &         (WORK(KTMP+IJ),IJ=0,STEX_N2VIRXQ-1)
C     Add d_ij F_ab part
               IF (KK.EQ.0.AND.II.EQ.JJ.AND.MOD(N-1,4).EQ.0) THEN
                  DO IZ=1,4
                     DO I=1,STEX_NVIRT
                        K = KTMP + (I-1)*STEX_NVIRT - 1 + 
     &                       (IZ-1)*STEX_N2VIRX
                        DO J=1,STEX_NVIRT
                           WORK(K+J) = WORK(K+J) + 
     &                          FMAT(IVIROFF+J,IVIROFF+I,IZ)
                        ENDDO
                     ENDDO                       
                  ENDDO
               ENDIF
C     Add d_ab F_ij part
               IF (KK.EQ.0.OR.KK.EQ.2) THEN
                  CSIGN=-1.0D0
               ELSE
                  CSIGN=1.0D0
               ENDIF
               DO IJ=0,STEX_NVIRT-1
                  WORK(KTMP+IJ+STEX_NVIRT*IJ) = 
     &                 WORK(KTMP+IJ+STEX_NVIRT*IJ) +
     &                 CSIGN*FMAT(STEX_HOLES(II),STEX_HOLES(JJ),KK+1)
               ENDDO
C     Write out modified integrals
               WRITE (STEX_LUINT,REC=N)
     &         (WORK(KTMP+IJ),IJ=0,STEX_N2VIRXQ-1)
               N = N + 1
            ENDDO
         ENDDO
      ENDDO
      END

      SUBROUTINE STEX_ADD_VPART(AMAT,FMAT,IZ,FAC)
C     Add the IZ part of FMAT, times FAC, to AMAT.
C     perhaps this could be done with BLAS with the odense
C     simplification (only virtual part stored on disk)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcbstex.h"
      DIMENSION AMAT(STEX_NPARAM,*), FMAT(STEX_NVIRT,STEX_NVIRT,NZ)
      DO II=1,STEX_NVIRT
         CALL DAXPY(STEX_NVIRT,FAC,FMAT(1,II,IZ),1,AMAT(1,II),1)
      ENDDO
      END

      SUBROUTINE STEX_BUILD_A(AMAT,WORK,LWORK)
C
C     Build the A matrix.
C     Vectors are ordered as V(NVIRT,4,NHOLE)
C     AMAT - Output CI Hamiltonian (A matrix)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcbstex.h"
      DIMENSION WORK(LWORK),
     $     AMAT(STEX_NVIRT,4,STEX_NHOLES,STEX_NVIRT,4,STEX_NHOLES)
#include "memint.h"
      CALL QENTER('STEX_BUILD_A')
      CALL DZERO(AMAT,STEX_NPARAM**2)
      CALL MEMGET('REAL',KTMP,STEX_N2VIRXQ,WORK,KFREE,LFREE)
      IFMAT = 0
      DO II=1,STEX_NHOLES
         DO JJ=II,STEX_NHOLES
            IFMAT = IFMAT + 1
            READ (STEX_LUINT, REC=IFMAT)
     &      (WORK(KTMP+KK),KK=0,STEX_N2VIRXQ-1)
C     FIXME: expand to quaternion format if NZ < 4 
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,1,JJ),WORK(KTMP),1,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,2,JJ),WORK(KTMP),2,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,3,JJ),WORK(KTMP),3,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,4,JJ),WORK(KTMP),4,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,2,JJ),WORK(KTMP),1,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,3,JJ),WORK(KTMP),4,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,4,JJ),WORK(KTMP),3,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,3,JJ),WORK(KTMP),1,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,4,JJ),WORK(KTMP),2,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,4,JJ),WORK(KTMP),1,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,1,JJ),WORK(KTMP),2,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,1,JJ),WORK(KTMP),3,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,2,JJ),WORK(KTMP),4,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,1,JJ),WORK(KTMP),4,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,2,JJ),WORK(KTMP),3,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,3,JJ),WORK(KTMP),2,-1.0D0)

            IFMAT = IFMAT + 1
            READ (STEX_LUINT, REC=IFMAT)
     &      (WORK(KTMP+KK),KK=0,STEX_N2VIRXQ-1)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,1,JJ),WORK(KTMP),1,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,2,JJ),WORK(KTMP),2,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,3,JJ),WORK(KTMP),3,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,4,JJ),WORK(KTMP),4,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,2,JJ),WORK(KTMP),1,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,3,JJ),WORK(KTMP),4,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,4,JJ),WORK(KTMP),3,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,3,JJ),WORK(KTMP),1,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,4,JJ),WORK(KTMP),2,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,4,JJ),WORK(KTMP),1,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,1,JJ),WORK(KTMP),2,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,1,JJ),WORK(KTMP),3,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,2,JJ),WORK(KTMP),4,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,1,JJ),WORK(KTMP),4,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,2,JJ),WORK(KTMP),3,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,3,JJ),WORK(KTMP),2,1.0D0)

            IFMAT = IFMAT + 1
            READ (STEX_LUINT, REC=IFMAT)
     &      (WORK(KTMP+KK),KK=0,STEX_N2VIRXQ-1)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,1,JJ),WORK(KTMP),3,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,2,JJ),WORK(KTMP),4,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,3,JJ),WORK(KTMP),1,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,4,JJ),WORK(KTMP),2,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,2,JJ),WORK(KTMP),3,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,3,JJ),WORK(KTMP),2,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,4,JJ),WORK(KTMP),1,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,3,JJ),WORK(KTMP),3,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,4,JJ),WORK(KTMP),4,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,4,JJ),WORK(KTMP),3,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,1,JJ),WORK(KTMP),4,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,1,JJ),WORK(KTMP),1,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,2,JJ),WORK(KTMP),2,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,1,JJ),WORK(KTMP),2,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,2,JJ),WORK(KTMP),1,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,3,JJ),WORK(KTMP),4,-1.0D0)

            IFMAT = IFMAT + 1
            READ (STEX_LUINT, REC=IFMAT)
     &      (WORK(KTMP+KK),KK=0,STEX_N2VIRXQ-1)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,1,JJ),WORK(KTMP),3,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,2,JJ),WORK(KTMP),4,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,3,JJ),WORK(KTMP),1,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,1,II,1,4,JJ),WORK(KTMP),2,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,2,JJ),WORK(KTMP),3,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,3,JJ),WORK(KTMP),2,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,4,JJ),WORK(KTMP),1,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,3,JJ),WORK(KTMP),3,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,4,JJ),WORK(KTMP),4,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,4,JJ),WORK(KTMP),3,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,2,II,1,1,JJ),WORK(KTMP),4,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,1,JJ),WORK(KTMP),1,-1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,3,II,1,2,JJ),WORK(KTMP),2,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,1,JJ),WORK(KTMP),2,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,2,JJ),WORK(KTMP),1,1.0D0)
            CALL STEX_ADD_VPART(AMAT(1,4,II,1,3,JJ),WORK(KTMP),4,1.0D0)

        ENDDO
      ENDDO
      IF (STEX_NOCOUPLING) THEN
C     Zero off-diagonal blocks.
         DO IH=1,STEX_NHOLES
            DO JH=IH+1,STEX_NHOLES
               PRINT *,'Zeroing coupling block',IH,JH
               DO IB=1,4
                  DO JB=1,4
                     DO II=1,STEX_NVIRT
                        DO JJ=1,STEX_NVIRT
                           AMAT(JJ,JB,IH,II,IB,JH) = 0.0D0
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF (STEX_DIAG) THEN
c     Zero all non-diagonal elements. Mainly for testing.
         CALL FDIAG(AMAT,STEX_NPARAM,STEX_NPARAM)
      ELSE
C     Make symmetric
         CALL USYM(AMAT,STEX_NPARAM,STEX_NPARAM)
      ENDIF
      CALL QEXIT('STEX_BUILD_A')
      END

      SUBROUTINE STEX_TM(XPRPTM,AMO,ITIM,IHOLE)
c
c     Calculate the transition moment elements
c
c     sqrt(2)*
c     <0|A[q+_ai + q+_a'i']|0>
c     <0|A[q+_ai - q+_a'i']|0>
c     <0|A[q+_a'i - q+_ai']|0>
c     <0|A[q+_a'i + q+_ai']|0>
c
C     The index a runs over all included virtual orbitals
C
C     In:
C     AMO - Property matrix in mo basis
C     ITIM - Timer reversal symmetry of A
C     IHOLE - the hole index "i".
C
C     Out:
C     XPRPTM(STEX_NVIRT,4) - in the order written above
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcbstex.h"
      DOUBLE COMPLEX XPRPTM(STEX_NVIRT,4)
      DIMENSION AMO(NORBT,NORBT,NZ)
      CALL QENTER('STEX_TM')
      IF (ITIM.NE.1) THEN
         CALL QUIT('ITIM=-1 not implemented, sorry.')
      ENDIF
c     FIXME: symmetry
      IF (NFSYM .EQ. 2) THEN
         CALL QUIT('NFSYM.eq.2 not implemented, sorry.')
      END IF
      IF (NZ .NE. 4) THEN
         CALL QUIT('NZ != 4 not implemented, sorry.')
      END IF
      DSQRT2 = SQRT(2.0D0)
      DO IVIRT=1,STEX_NVIR(1)
         XPRPTM(IVIRT,1) =  DSQRT2*AMO(STEX_HOLES(IHOLE),
     &        NPSH(1) + NISH(1) + IVIRT,1)
      ENDDO
      DO IVIRT=1,STEX_NVIRT
         XPRPTM(IVIRT,2) = -DSQRT2*AMO(STEX_HOLES(IHOLE),
     &        NPSH(1) + NISH(1) + IVIRT,2)
      ENDDO
      DO IVIRT=1,STEX_NVIRT
         XPRPTM(IVIRT,3) =  DSQRT2*AMO(STEX_HOLES(IHOLE),
     &        NPSH(1) + NISH(1) + IVIRT,3)
      ENDDO
      DO IVIRT=1,STEX_NVIRT
         XPRPTM(IVIRT,4) = -DSQRT2*AMO(STEX_HOLES(IHOLE),
     &        NPSH(1) + NISH(1) + IVIRT,4)
      ENDDO
      CALL QEXIT('STEX_TM')
      END

      SUBROUTINE FROZEN_PRPGRAD(XPRPGRAD,CMO,IBEIG,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcbstex.h"
      DOUBLE COMPLEX XPRPGRAD
      DIMENSION CMO(NCMOTQ),WORK(LWORK),IBEIG(*),
     &     XPRPGRAD(STEX_NVIRT,4,STEX_NHOLES,STEX_NR_OP)
      LOGICAL LOPEN
#include "memint.h"
      CALL QENTER('FROZEN_PRPGRAD')
      CALL MEMGET('REAL',KPROP,N2ORBXQ,WORK,KFREE,LFREE)
      INQUIRE(FILE='AOPROPER',OPENED=LOPEN)
      IF (.NOT. LOPEN) THEN
         OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &        FILE='AOPROPER')
      ENDIF
      DO IPRP=1,STEX_NR_OP
         CALL WRIXPR(IPRP,STEX_OP_IND(IPRP))
         call flush(lupri)
         CALL PRPMAT(STEX_OP_IND(IPRP),IDUM,WORK(KPROP),.TRUE.,WORK,CMO,
     &               IBEIG,ICMOQ,NORB,WORK,KFREE,LFREE,1)
c         CALL HEADER('FROZEN_PRPGRAD: Property matrix in MO-basis',-1)
c         CALL PRQMAT(WORK(KPROP),NORBT,NORBT,NORBT,NORBT,NZ,
c     &        IPQTOQ(1,0),LUPRI)
         IF (IPRP.LE.3) THEN
            ITIM=1
         ELSE
C     Check this for dipvel (p?)
            ITIM=1
         ENDIF
         DO IHOLE=1,STEX_NHOLES
            CALL STEX_TM(XPRPGRAD(1,1,IHOLE,IPRP),WORK(KPROP), ITIM,
     $           IHOLE)
         ENDDO
      ENDDO
      CALL QEXIT('FROZEN_PRPGRAD')
      END

      DOUBLE COMPLEX FUNCTION DZDOT(N,DVEC,ZVEC)
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION DVEC
      DOUBLE COMPLEX ZVEC,SUM
      DIMENSION DVEC(N),ZVEC(N)
      SUM = 0
      DO I=1,N
         SUM = SUM + DVEC(I)*ZVEC(I)
      ENDDO
      DZDOT = SUM
      END

      SUBROUTINE STEX_ANALYZE(NVEC,EVEC,EVAL,PRPGRAD,GSERG,REFERG
     $     ,REFGSOVL,WORK,LWORK)
C
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcbstex.h"
      DOUBLE COMPLEX XMOM,PRPGRAD,DZDOT,REFGSOVL
      DIMENSION EVEC(STEX_NVIRT,4,STEX_NHOLES,NVEC),EVAL(NVEC),
     $     PRPGRAD(STEX_NVIRT,4,STEX_NHOLES,STEX_NR_OP),WORK(LWORK)
     $     ,OSCL(3), OSCV(3),DOMCHANW(50),
     &     IDOMCHANH(50),IDOMCHANV(50), HOLEWEIGHT(STEX_NHOLES)
#include "memint.h"
      CALL QENTER('STEX_ANALYZE')
      WRITE(LUPRI,'(1A)') 'Excitation energies and oscillator strengths'
      WRITE(LUPRI,'(1A)') ' State    Energy (Ev) '//
     &     ' <f_l>/3   <f_v>/3     f_lx      f_ly  '//
     &     '    f_lz      f_vx      f_vy      f_vz'
      WRITE(LUPRI,'(160A1)') ('-',I=1,160)
      DERG = REFERG - GSERG
C     Check sum rule for oscillator strengths
      SUML = 0.0D0
      SUMV = 0.0D0
      DO IVEC=1,NVEC
c     transition moments
         DO IPRP=1,3
c     Length gauge
            XMOM = DZDOT(STEX_NPARAM,EVEC(1,1,1,IVEC),
     &           PRPGRAD(1,1,1,IPRP))
            OSCL(IPRP) = 2.0D0*ABS(XMOM)**2*(EVAL(IVEC) + DERG)
C     Velocity gauge
            XMOM = DZDOT(STEX_NPARAM,EVEC(1,1,1,IVEC),
     &           PRPGRAD(1,1,1,IPRP+3))
            OSCV(IPRP) = 2.0D0*ABS(XMOM)**2/(EVAL(IVEC) + DERG)
         ENDDO
         WRITE(LUPRI,'(1I5,1F16.8,8(3X,E14.6))')
     &        IVEC,(EVAL(IVEC) + DERG)*27.211383,
     &        (OSCL(1)+OSCL(2)+OSCL(3))/3.0D0,
     &        (OSCV(1)+OSCV(2)+OSCV(3))/3.0D0,
     &        OSCL(1),OSCL(2),OSCL(3),OSCV(1),OSCV(2),OSCV(3)
         SUML = SUML + (OSCL(1)+OSCL(2)+OSCL(3))/3.0D0
         SUMV = SUMV + (OSCV(1)+OSCV(2)+OSCV(3))/3.0D0
      ENDDO
      WRITE(LUPRI,*) '-------------------------------------'/
     $     /'---------------------------------------'
      WRITE(LUPRI,'(1A,2F10.6)')
     $     '   Sum of oscillator strengths (length, velocity gauge): '
     $     ,SUML ,SUMV

      WRITE (LUPRI,'(//A)')
     &   ' Dominant excitation channels for each state '//
     &   '[hole -> virtual (weight)]'
      DO IVEC=1,NVEC
C     Find the dominant channels
         NDOMCHAN = 0
         DO IH=1,STEX_NHOLES
            DO IV=1,STEX_NVIRT
               SUM = 0.0D0
               DO IT=1,4
                  SUM = SUM + EVEC(IV,IT,IH,IVEC)**2
               ENDDO
               IF (SUM.GE.0.05D0) THEN
                  NDOMCHAN = NDOMCHAN + 1
                  DOMCHANW(NDOMCHAN) = SUM
                  IDOMCHANH(NDOMCHAN) = IH
                  IDOMCHANV(NDOMCHAN) = IV
               ENDIF
            ENDDO
         ENDDO
         WRITE(LUPRI,'(1I5,50(1I3,1A,1I4,1A,1F4.2,1A))') IVEC,
     &        (IDOMCHANH(II),' ->',IDOMCHANV(II),' (',DOMCHANW(II),
     &        ')',II=1,NDOMCHAN)
      ENDDO
      IF (STEX_NHOLES .GT. 1) THEN
        WRITE (LUPRI,'(/A)') ' Hole contribution to each state'
        DO IVEC=1,NVEC
         DO IH=1,STEX_NHOLES
            HOLEWEIGHT(IH) = 0.0D0
            DO IV=1,STEX_NVIRT
               DO IT=1,4
                  HOLEWEIGHT(IH) =
     &                 HOLEWEIGHT(IH) + EVEC(IV,IT,IH,IVEC)**2
               ENDDO
            ENDDO
         ENDDO
         WRITE(LUPRI,'(1I5,50(1I3,1A,1F6.3,1A))') IVEC,
     &        (IH,' ',HOLEWEIGHT(IH),' ',IH=1,STEX_NHOLES)
        ENDDO
      END IF
      CALL QEXIT('STEX_ANALYZE')
      END


      SUBROUTINE COF_ADJUGATE(N,LDA,A,TMP)
C     Compute the adjugate of the NxN double complex matrix A by
C     using a SVD factorization and perturbing the smallest singular
C     value if it is too close to zero.
C     TMP should be of size N+2*N**2+3*N*4+5*N (double complex)
      IMPLICIT NONE
      INTEGER LDA,N,LWORK
      DOUBLE COMPLEX A,TMP
      DIMENSION A(LDA,N),TMP(*)
c     assume that a double complex is twice the size of a double
C     LWORK should be at least 3*N double complex numbers
      LWORK = 3*N *4
      CALL COF_ADJUGAT1(N,LDA,A,TMP(1),TMP(N+1),TMP(N+N**2+1),
     & TMP(N+2*N**2+1),LWORK,TMP(N+2*N**2+LWORK+1))
      END

      SUBROUTINE COF_ADJUGAT1(N,LDA,A,SIGMA,U,VT,WORK,LWORK,RWORK)
C     A = USV'
C     adj(A) = det(U)det(V)det(S)VS^-1U'
      IMPLICIT NONE
      INTEGER I,J,N,LWORK,INFO,LDA
      DOUBLE PRECISION SIGMA,RWORK
      DOUBLE COMPLEX A,U,VT,WORK,DETS,ZDET,ZZERO,ZONE
      DIMENSION A(LDA,N),SIGMA(N),U(N,N),VT(N,N),WORK(LWORK),
     &     RWORK(5*N)
      ZZERO = 0.0
      ZONE = 1.0
      CALL ZGESVD('A','A',N,N,A,LDA,SIGMA,U,N,VT,N,WORK,LWORK,RWORK
     &     ,INFO)
      IF (INFO.NE.0) THEN
         PRINT *,'COF_ADJUGATE: ZGESVD returned ',INFO
      ENDIF
C     Compute det(sigma)
c      print *,'sigma',sigma
      DETS = 1.0D0
      DO I=1,N
         IF (ABS(SIGMA(I)).LT.1e-14) THEN
            SIGMA(I) = 1e-14
         ENDIF
         DETS = DETS*SIGMA(I)
      ENDDO
C     Compute VS^-1, i.e. divide each column of V by sigma(i)
C     this is equivalent to dividing each row of VT.
      DO J=1,N
         DO I=1,N
            VT(I,J) = VT(I,J)/SIGMA(I)
         ENDDO
      ENDDO
C     Compute VS^-1U' = (S^-1V')'U' = VT'*U' (VT = S^-1V' here)
C     store the result back in A
      CALL ZGEMM('C','C',N,N,N,ZONE,VT,N,U,N,ZZERO,A,LDA)
C     det(U) and det(V)
C     VT is now multiplied by S^-1, undo that operation
      DO J=1,N
         DO I=1,N
            VT(I,J) = VT(I,J)*SIGMA(I)
         ENDDO
      ENDDO
      CALL ZDETER(U,N,N,ZDET)
      DETS = DETS*ZDET
      CALL ZDETER(VT,N,N,ZDET)
c      DETS = DETS*CONJG(ZDET) ??
      DETS = DETS*ZDET
      DO J=1,N
         DO I=1,N
            A(I,J) = A(I,J)*DETS
         ENDDO
      ENDDO
      END

      SUBROUTINE ZDETER(A,LDA,N,DET)
C     Compute determinant of general double complex matrix A.
      IMPLICIT NONE
      INTEGER N,IWORK,INFO,LDA
      DOUBLE COMPLEX A,DET,WORK,VDET
      DIMENSION A(LDA,N),IWORK(N),WORK(N),VDET(2)
      CALL ZGEFA(A,LDA,N,IWORK,INFO)
      CALL ZGEDI(A,LDA,N,IWORK,VDET,WORK,10)
      DET = VDET(1)*10.0**vdet(2)
      END

      SUBROUTINE ZINVER(A,N,DET)
C     Compute determinant and inverse of general double complex matrix A.
      IMPLICIT NONE
      INTEGER N,IWORK,INFO
      DOUBLE COMPLEX A,DET,WORK,VDET
      DIMENSION A(N,N),IWORK(N),WORK(N),VDET(2)
      DET = 0.0
      CALL ZGEFA(A,N,N,IWORK,INFO)
      IF (INFO.NE.0) RETURN
      CALL ZGEDI(A,N,N,IWORK,VDET,WORK,11)
      DET = VDET(1)*10.0**vdet(2)
      END

      SUBROUTINE QTOC(QMAT,CMAT,ITIM,NROW,NCOL,NLDR,NLDC,MZ)
C     Convert a quaternion matrix
C     QMAT = A+Bj
C     to complex form
C     CMAT = ( A1    B1   A2    B2  ..)
C            (-tB1* tA1* -tB2* tA2* ..)
C     where t = ITIM is the time reversal character
C     NCOL and NROW is the dimension of A (and B).
C     CMAT is stored in interlaced format, so that columns C_ai and C_ai' are
C     next to each other. The rows are also interlaced, so we have
C     Q_ai -> ( C_ai  C_ai'  )
C             ( C_a'i C_a'i' )
C
      IMPLICIT NONE
      INTEGER I,J,NROW,NCOL,ITIM,NLDR,NLDC,MZ
      DOUBLE PRECISION QMAT
      DOUBLE COMPLEX CMAT
      DIMENSION QMAT(NLDR,NLDC,MZ), CMAT(NROW*2, NCOL*2)
C     FIXME: implement MZ < 4
      IF (MZ.NE.4) THEN
         CALL QUIT('QTOC: NZ != 4')
      ENDIF
      DO J = 1,NCOL
         DO I = 1,NROW
C     A
            CMAT(I*2-1,J*2-1) = DCMPLX(QMAT(I,J,1), QMAT(I,J,2))
C     -tB*
            CMAT(I*2,J*2-1)=-ITIM*DCMPLX(QMAT(I,J,3),-QMAT(I,J,4))
         ENDDO
         DO I = 1,NROW
C     B
            CMAT(I*2-1,J*2) = DCMPLX(QMAT(I,J,3), QMAT(I,J,4))
C     tA*
            CMAT(I*2,J*2) = ITIM*DCMPLX(QMAT(I,J,1),-QMAT(I,J,2))
         ENDDO
      ENDDO
      END


      SUBROUTINE STEX_FAOMO2(FMO,FAO,MZ,ISYM,ITIM,IFER,CMO1,CMO2,
     &     WORK,LWRK)
C     Compute <i|F|j>, where i are orbitals from CMO1 and j are obitals
C     from CMO2
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcbxqr.h"
C
      PARAMETER ( D0 = 0.0D0 )
C
      CHARACTER T*1
C
      DIMENSION FAO(N2BBASX,MZ),FMO(N2ORBX,MZ)
      DIMENSION T(-1:1),CMO1(*),CMO2(*),WORK(LWRK)
C
      CALL QENTER('STEX_FAOMO2')
C
      KFREE = 1
      LFREE = LWRK
C
      T(1)='S'
      T(-1)='A'
C
      CALL DZERO(FMO,N2ORBXQ)
C
      IPRINT=0
      DO I1=1,NFSYM
         I2 = MOD(I1+IFER,2) + 1
         IF(NORB(I1).GT.0 .AND. NORB(I2).GT.0) THEN
            CALL QTRANS('AOMO',T(ITIM),
     &           D0,NFBAS(I1,0),NFBAS(I2,0),NORB(I1),NORB(I2),
     &           FAO(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),MZ,
     &           IPQTOQ(1,ISYM-1),
     &           FMO(I2ORBX(I1,I2)+1,1),NORBT,NORBT,NZ,
     &           IPQTOQ(1,ISYM-1),
     &           CMO1(ICMOQ(I1)+1),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &           CMO2(ICMOQ(I2)+1),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
         END IF
      END DO
      CALL QEXIT('STEX_FAOMO2')
C
      RETURN
      END

      SUBROUTINE ZCOLSWAP(A,ROWS,COLS,I,J)
      IMPLICIT NONE
      DOUBLE COMPLEX A,TMP
      INTEGER ROWS,COLS,I,J,K
      DIMENSION A(ROWS,COLS)
      DO K=1,ROWS
         TMP = A(K,I)
         A(K,I) = A(K,J)
         A(K,J) = TMP
      ENDDO
      END

      SUBROUTINE ZROWSWAP(A,ROWS,COLS,I,J)
      IMPLICIT NONE
      DOUBLE COMPLEX A,TMP
      INTEGER ROWS,COLS,I,J,K
      DIMENSION A(ROWS,COLS)
      DO K=1,COLS
         TMP = A(I,K)
         A(I,K) = A(J,K)
         A(J,K) = TMP
      ENDDO
      END

      SUBROUTINE STEX_PRPMO12(AQ,AC,OVL,WORK,LWORK)
C     Given a property matrix in mo basis (AQ), and an overlap
C     matrix OVL = <i1|j2>,
C     Place a complexified version of <i1|A|j2> in AC. i runs over
C     all inactive orbitals, while j runs over inactive and virtual.
C     AQ12 = sum j <i1|A|j1><j1|k2> = AQ*S
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
      DOUBLE COMPLEX AC
      DIMENSION AC(*),AQ(N2ORBXQ),OVL(N2ORBXQ),WORK(LWORK)
#include "memint.h"
      CALL MEMGET('REAL',KAQ12, N2ORBXQ, WORK,KFREE,LFREE)
      CALL QGEMM(NORBT,NORBT,NORBT,1.0D0,
     &     'N','N',IPQTOQ(1,0),AQ,NORBT,NORBT,NZ,
     &     'N','N',IPQTOQ(1,0),OVL,NORBT,NORBT,NZ,
     $     0.0D0,IPQTOQ(1,0),WORK(KAQ12),NORBT,NORBT,NZ)
      NOCCT = NISH(1)
      NOCCVIRT = NORB(1) - NPSH(1)
      CALL QTOC(WORK(KAQ12+NPSH(1)*(1+NORBT)),AC,1,NOCCT
     $     ,NOCCVIRT,NORBT,NORBT,NZ)
      END


      SUBROUTINE STEX_OVL12(QOVL,CMO1,CMO2,WORK,LWORK)
C     Calculate the overlap between the mo's in cmo1 and cmo2
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
      DIMENSION CMO1(NCMOTQ),CMO2(NCMOTQ),QOVL(N2ORBXQ),WORK(LWORK)
#include "memint.h"
      CALL MEMGET('REAL',KOVL, N2BBASX, WORK,KFREE,LFREE)
      CALL GTOVLX(WORK(KOVL),SSMTRC)
      CALL STEX_FAOMO2(QOVL,WORK(KOVL),1,1,1,1,CMO1,CMO2,WORK(KFREE)
     $     ,LFREE)
      END

      FUNCTION ZCOF(A,N,I,J)
C     Compute the cofactor (i,j) of A.
C     A is destroyed in the process.
      IMPLICIT NONE
      INTEGER N,I,J,II,JJ
      DOUBLE COMPLEX A,ZCOF
      DIMENSION A(N,N)
      DO II=1,I-1
         DO JJ=1,J-1
            A(II,JJ) = A(II,JJ)
         ENDDO
         DO JJ=J,N-1
            A(II,JJ) = A(II,JJ+1)
         ENDDO
      ENDDO
      DO II=I,N-1
         DO JJ=1,J-1
            A(II,JJ) = A(II+1,JJ)
         ENDDO
         DO JJ=J,N-1
            A(II,JJ) = A(II+1,JJ+1)
         ENDDO
      ENDDO
      CALL ZDETER(A,N,N-1,ZCOF)
      ZCOF = (-1)**(I+J)*ZCOF
      END

      SUBROUTINE COF_SUM(NORB,NOVORB,PRP,NPRP,S,SINV,SDET
     $     ,NUCDIP,MOM,IHOLE,IVIRT,WORK,LWORK)
C     Given a property matrix <i1|A|j2>, and the overlap matrix
C     <i1|j2>, between two sets of orbitals, compute the
C     transition moments <1|A|2> through cofactor expansion.
C     <1|2> = det(S)
C     <1|A|2> = sum_ij <i1|A|j2> * cof_ij S
C     cof_ij S = det S * (S^-1)_ji (if det S != 0)
C     SINV should contain the inverse of S, if SDET != 0.
C     Else SINV will be used as scratch space in cofactor
C     computations.
C     NOVORB is number of occ+virt orbitals
C     NORB is number of occupied
      IMPLICIT NONE
      DOUBLE COMPLEX PRP,MOM,S,SINV,SDET,ZCOF,COF
      INTEGER NORB,NOVORB,I,J,K,L,IPRP,NPRP,IHOLE,IVIRT,LWORK
      DOUBLE PRECISION NUCDIP,WORK
      LOGICAL LCOF,LSVD
      DIMENSION PRP(NORB,NOVORB,NPRP),SINV(NORB,NORB),S(NORB,NORB),
     &     MOM(NPRP), NUCDIP(3)
      DO IPRP=1,NPRP
         MOM(IPRP) = 0.0
      ENDDO
C     check if we have to calculate cof. explicitly.
      IF (ABS(SDET).LT.1e-8) THEN
         LCOF = .TRUE.
      ELSE
         LCOF = .FALSE.
      ENDIF
c     new cofactor evaluation by svd
      LSVD = .FALSE.
      IF (LCOF) THEN
c         print *,'calling adjugate, det = ',ABS(SDET)
         DO K=1,NORB
            DO L=1,NORB
               SINV(L,K) = S(L,K)
            ENDDO
         ENDDO
         CALL COF_ADJUGATE(NORB,NORB,SINV,WORK)
         LSVD = .TRUE.
      ENDIF
c     TODO: can we always use the SVD approach now?
      DO I=1,NORB
         DO J=1,NORB
            IF (LCOF) THEN
               IF (LSVD) THEN
                  COF = SINV(J,I)
               ELSE
C     Copy S to SINV, because zcof destroys its argument
                  DO K=1,NORB
                     DO L=1,NORB
                        SINV(L,K) = S(L,K)
                     ENDDO
                  ENDDO
                  COF = ZCOF(SINV,NORB,I,J)
               ENDIF
            ELSE
               COF = SDET*SINV(J,I)
            ENDIF
            DO IPRP=1,NPRP
               IF (J.EQ.IHOLE) THEN
                  MOM(IPRP)=MOM(IPRP)+COF*PRP(I,IVIRT,IPRP)
               ELSE
                  MOM(IPRP)=MOM(IPRP)+COF*PRP(I,J,IPRP)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C     add nuclear dipole moment term for the length gauge
      DO IPRP=1,3
         MOM(IPRP) = MOM(IPRP) - SDET*NUCDIP(IPRP)
      ENDDO
      END


      SUBROUTINE ZMATNORM(NAME,A,ROWS,COLS)
      IMPLICIT NONE
      DOUBLE COMPLEX A
      DOUBLE PRECISION NOR
      INTEGER ROWS,COLS,II,JJ
      CHARACTER NAME*(*)
      DIMENSION A(ROWS,COLS)
      NOR = 0.0D0
      DO JJ=1,COLS
         DO II=1,ROWS
            NOR = NOR + ABS(A(II,JJ))**2
         ENDDO
      ENDDO
      PRINT *,'ZMATNORM ',NAME,' ',SQRT(NOR)
      END

      SUBROUTINE COF_PRPGRAD(XPRPGRAD,CMO1,IBREP,CMO2,REFGSOVL,WORK
     $     ,LWORK)
C     Calculate the overlap between the determinants formed from
C     the occupied orbitals of CMO1 and the occupied orbitals of
C     CMO2, where the IHOLE orbital in CMO2 has been replace with
C     the IVIRT orbital. The indices start at electronic orbitals.
C     IHBAR and IVBAR indicate if it is the "primary" (0) orbital or its
C     Kramers partner (1) that is ment.
C     Input:
C     CMO1 - the groundstate
C     CMO2 - the reference state
C     Output:
C     REFGSOVL - overlap between CMO1 and CMO2
C     XPRPGRAD - the property gradient
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dipole.h"
#include "dcbstex.h"
      LOGICAL LOPEN
      DOUBLE COMPLEX XPRPGRAD,DETS,MOM,REFGSOVL
      DIMENSION XPRPGRAD(STEX_NVIRT,4,STEX_NHOLES,STEX_NR_OP),
     $     WORK(LWORK),IBREP(*),CMO1(NCMOTQ),CMO2(NCMOTQ),MOM(STEX_NR_OP
     $     )
#include "memint.h"
      CALL QENTER('COF_PRPGRAD')
      WRITE(LUPRI,*) 'Calculating property gradient by cofactors'
      CALL FLUSH(LUPRI)
      NOCCT = NISH(1)
      NOCCVIRT = NORB(1) - NPSH(1)
C     NCSIZE - nr of doubles for the nocc*2 x noccvirt*2 dcmplx matrix
      NCSIZE = NOCCT*NOCCVIRT*8
      CALL MEMGET('REAL',KCOVL, NCSIZE, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQOVL, N2ORBXQ, WORK,KFREE,LFREE)
C     Get overlap between CMO1 and CMO2
      CALL STEX_OVL12(WORK(KQOVL),CMO1,CMO2,WORK(KFREE),LFREE)
C     Put the complex overlap matrix in WORK(KCOVL)
      CALL QTOC(WORK(KQOVL+NPSH(1)*(1+NORBT)),WORK(KCOVL),1,NOCCT,
     &     NOCCVIRT,NORBT,NORBT,NZ)
C     Load the property matrices and convert them to complex form.
      INQUIRE(FILE='AOPROPER',OPENED=LOPEN)
      IF (.NOT. LOPEN) THEN
         OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &        FILE='AOPROPER')
      ENDIF
      CALL MEMGET('REAL',KCPRP, NCSIZE*STEX_NR_OP, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQA, N2ORBXQ, WORK,KFREE,LFREE)

      DO IPRP=1,STEX_NR_OP
         CALL PRPMAT(STEX_OP_IND(IPRP),IDUM,WORK(KQA),.TRUE.,WORK,CMO1,
     &               IBREP,ICMOQ,NORB,WORK,KFREE,LFREE,0)
         CALL STEX_PRPMO12(WORK(KQA),WORK(KCPRP+(IPRP-1)*NCSIZE),
     &        WORK(KQOVL),WORK(KFREE),LFREE)
      ENDDO
c     Get <gs|ref>
      CALL DCOPY(NOCCT**2*8,WORK(KCOVL),1,WORK(KQA),1)
      CALL ZDETER(WORK(KQA),NOCCT*2,NOCCT*2,REFGSOVL)
C     Calculate nuclear dipole moment
      CALL DIPNUC(0,.FALSE.)
      WRITE (LUPRI,*) 'Nuclear dipole moment:',DIPMN
c     prepare for cofactors
      CALL DZERO(XPRPGRAD,STEX_NVIRT*4*STEX_NHOLES*STEX_NR_OP*2)
      DSQRT2I = 1.0D0/SQRT(2.0D0)
C     Loop over excitations. Indices starting with I are over holes
C     and those starting with J are over virtuals
      DO II=1,STEX_NHOLES
         DO IBAR=0,1
            IHOLE = (STEX_HOLES(II)-NPSH(1))*2 - 1 + IBAR
            DO JJ=1,STEX_NVIRT
               DO JBAR=0,1
c                  PRINT *,'from ',II,IBAR,' to ',JJ,JBAR
                  IVIRT = NOCCT*2 + JJ*2 - 1 + JBAR
C     Swap in columns to make a compact matrix for the excited state
C     determinant
                  CALL ZCOLSWAP(WORK(KCOVL),NOCCT*2,NOCCVIRT*2,
     &                 IHOLE,IVIRT)
C     Try to invert the overlap matrix to get cofactors.
C     reuse KQOVL for this, since it's no longer needed.
                  CALL DCOPY(NOCCT**2*8,WORK(KCOVL),1,WORK(KQOVL),1)
                  CALL ZINVER(WORK(KQOVL),NOCCT*2,DETS)
                  CALL COF_SUM(NOCCT*2,NOCCVIRT*2,WORK(KCPRP),
     $                 STEX_NR_OP,WORK(KCOVL),WORK(KQOVL),DETS,
     $                 DIPMN,MOM,IHOLE,IVIRT,WORK(KFREE),LFREE)
c                  PRINT *,II,IBAR,JJ,JBAR,ABS(DETS)**2,MOM(1)
C     Now we have a moment <gs|A a+_j a_i |ref>, where
C     j is given by JJ and JBAR and i is given by II and IBAR.
C     Transform to the order given in the header of this subroutine.
                  DO IPRP=1,STEX_NR_OP
                     IF (IBAR.EQ.0.AND.JBAR.EQ.0) THEN
                        XPRPGRAD(JJ,1,II,IPRP) =
     &                       XPRPGRAD(JJ,1,II,IPRP) + MOM(IPRP)*DSQRT2I
                        XPRPGRAD(JJ,2,II,IPRP) =
     &                       XPRPGRAD(JJ,2,II,IPRP) +
     &                       DCMPLX(0,1)*MOM(IPRP)*DSQRT2I
                     ELSE IF (IBAR.EQ.1.AND.JBAR.EQ.1) THEN
                        XPRPGRAD(JJ,1,II,IPRP) =
     &                       XPRPGRAD(JJ,1,II,IPRP) + MOM(IPRP)*DSQRT2I
                        XPRPGRAD(JJ,2,II,IPRP) =
     &                       XPRPGRAD(JJ,2,II,IPRP) -
     &                       DCMPLX(0,1)*MOM(IPRP)*DSQRT2I
                     ELSE IF (IBAR.EQ.0.AND.JBAR.EQ.1) THEN
                        XPRPGRAD(JJ,3,II,IPRP) =
     &                       XPRPGRAD(JJ,3,II,IPRP) + MOM(IPRP)*DSQRT2I
                        XPRPGRAD(JJ,4,II,IPRP) =
     &                       XPRPGRAD(JJ,4,II,IPRP) +
     &                       DCMPLX(0,1)*MOM(IPRP)*DSQRT2I
                     ELSE
                        XPRPGRAD(JJ,3,II,IPRP) =
     &                       XPRPGRAD(JJ,3,II,IPRP) - MOM(IPRP)*DSQRT2I
                        XPRPGRAD(JJ,4,II,IPRP) =
     &                       XPRPGRAD(JJ,4,II,IPRP) +
     &                       DCMPLX(0,1)*MOM(IPRP)*DSQRT2I
                     ENDIF
                  ENDDO
C     Swap back columns
                  CALL ZCOLSWAP(WORK(KCOVL),NOCCT*2,NOCCVIRT*2,
     &                 IHOLE,IVIRT)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL QEXIT('COF_PRPGRAD')
      END


      SUBROUTINE PRINTGRAD(GRAD,NHOLES,NVIRT)
      IMPLICIT NONE
      DOUBLE COMPLEX GRAD
      DOUBLE PRECISION NRM
      INTEGER I,J,K,IPRP,NHOLES,NVIRT
      DIMENSION GRAD(NVIRT,4,NHOLES,6)
      DO IPRP=1,6
         PRINT *,'Property gradient ',IPRP
         DO K=1,4
            NRM = DCMPLX(0,0)
            DO J=1,NHOLES
               DO I=1,NVIRT
                  PRINT *,K,J,I,GRAD(I,K,J,IPRP)
                  NRM = NRM + ABS(GRAD(I,K,J,IPRP))**2
               ENDDO
            ENDDO
            PRINT *,'Norm ',K,' = ',NRM
         ENDDO
      ENDDO
      END

      SUBROUTINE STEXFAOMO(FMO,FAO,ISYM,ITIM,IFER,IPRINT,CMO,WORK,LWRK)
C*****************************************************************************
C
C     Transform the virtual part of a (Fock) matrix from AO to MO basis using QTRANS
C
C     Written by panor 1999/ulfek&hjj 2006
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcbstex.h"

      PARAMETER ( D0 = 0.0D0 )

      CHARACTER T*1

      DIMENSION FAO(N2BBASXQ),FMO(STEX_N2VIRXQ)
      DIMENSION T(-1:1),CMO(*),WORK(LWRK)
      INTEGER   JCMO(2)

      CALL QENTER('STEXFAOMO')
      IF (IPRINT.GE.10) THEN
         CALL HEADER('STEXFAOMO: Fock matrix in AO-basis',-1)
         CALL PRQMAT(FAO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &        IPQTOQ(1,ISYM-1),LUPRI)
      ENDIF

      KFREE = 1 
      LFREE = LWRK

      JCMO(1) = ICMOQ(1) + NFBAS(1,0)*(NPSH(1)+NISH(1)) + 1
      JCMO(2) = ICMOQ(2) + NFBAS(2,0)*(NPSH(2)+NISH(2)) + 1      

      CALL DZERO(FMO,STEX_N2VIRXQ)
      !IPRINT=0 ! mi-commented out (fixed openmpi gfortran 4.4.5 prob)
      DO I1=1,NFSYM
         I2 = MOD(I1+IFER,2) + 1
         IF(STEX_NVIR(I1).GT.0 .AND. STEX_NVIR(I2).GT.0) THEN
            CALL QTRANS('AOMO','S',
     &           D0,NFBAS(I1,0),NFBAS(I2,0),STEX_NVIR(I1),STEX_NVIR(I2),
     &           FAO(I2BASX(I1,I2)+1),NTBAS(0),NTBAS(0),NZ,
     &           IPQTOQ(1,ISYM-1),
     &           FMO(STEX_I2VIRX(I1,I2)+1),STEX_NVIRT,STEX_NVIRT,NZ,
     &           IPQTOQ(1,ISYM-1),
     &           CMO(JCMO(I1)),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &           CMO(JCMO(I2)),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
         END IF
      END DO
C
C     Print
C
      IF (IPRINT.GE.10) THEN
         CALL HEADER('STEXFAOMO: Fock matrix in MO-basis',-1)
         CALL PRQMAT(FMO,STEX_NVIRT,STEX_NVIRT,STEX_NVIRT,STEX_NVIRT,NZ,
     &        IPQTOQ(1,ISYM-1),LUPRI)
      ENDIF
C
      CALL QEXIT('STEXFAOMO')
C     
      RETURN
      END

#ifdef NOTUSED
      SUBROUTINE ORBDENS(CMO,IIORB,DAO,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
      DIMENSION WORK(LWORK),CMO(NTBAS(0),NORBT,NZ),DAO(N2BBASXQ)
      DIMENSION CORB1(NTBAS(0),NZ),CORB2(NTBAS(0),NZ)
      CHARACTER OP
#include "memint.h"

      DO IZ=1,NZ
         DO II=1,NTBAS(0)
            CORB1(II,IZ) = CMO(II,IIORB,IZ)
         ENDDO
      ENDDO

      DO IZ=1,NZ
         DO II=1,NTBAS(0)
            CORB2(II,IZ) = CMO(II,IIORB,IZ)
         ENDDO
      ENDDO

      OP = 'N'

      CALL QGEMM(NTBAS(0),NTBAS(0),1,1.0D0,
     &     'N',OP,IPQTOQ(1,0),CORB1,NTBAS(0),1,NZ,
     &     'H','N',IPQTOQ(1,0),CORB2,NTBAS(0),1,NZ,
     &     0.0D0,IPQTOQ(1,0),DAO,NTBAS(0),NTBAS(0),NZ)
      END
#endif
