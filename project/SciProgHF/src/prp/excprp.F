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

C***********************************************************************
C
C     Driver routine for computing matrix norm.
C
C     Written by Ulf Ekstrom
C***********************************************************************
      SUBROUTINE QMATNORM(NAME,A,ROWS,COLS,NZ)
      IMPLICIT NONE
#include "priunit.h"
      DOUBLE PRECISION A,NOR
      INTEGER ROWS,COLS,II,JJ,IZ,NZ
      CHARACTER NAME*(*)
      DIMENSION A(ROWS,COLS,NZ)
      NOR = 0.0D0
      DO IZ=1,NZ
         DO JJ=1,COLS
            DO II=1,ROWS
               NOR = NOR + ABS(A(II,JJ,IZ))**2
            ENDDO
         ENDDO
      ENDDO
      WRITE(LUPRI,'(3A,F14.8)') ' QMATNORM: Norm of ',NAME,' is ',
     &     SQRT(NOR)
      END
      SUBROUTINE EXCPRP(WORK,LWRK)
C***********************************************************************
C
C     Driver routine for computing excited state properties and
C     transition moments of the form < LEFT | OP | RIGHT >
C
C     Each column of the matrix tmeqlist is a specification of a 
C     transition moment. The data contained in a column is:
C     Oper., Leftsym., Left excit., Rightsym., Right excit.
C
C     Written by Erik Tellgren, aug 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxpr.h"
#include "dcbxqr.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
C
      INTEGER LUCMO
      PARAMETER ( LUCMO = 22 )
      REAL*8 RES_LR(3,MAXLQR,MAXLQR,MAXFQR)
      REAL*8 RES_EXCPRP(3,NBSYM,MAXFQR,NBSYM,MAXFQR,MAXLQR)
      REAL*8 WORK(*)
      INTEGER tmeqlist(5,MAXQR)
      INTEGER numleftst, numrightst, I
      INTEGER nexcsym(8),nexcsim(8),nexcstv(8)
      INTEGER LWRK,KFREE,LFREE,KCMO,KIBRP,IDUM
      INTEGER KVECAEE,KVECAEP,KVECBEE,KVECBEP,KVECCEE,KVECCEP
      REAL*8 DUM,TOTERG
      LOGICAL prponly
C
      CALL QENTER('EXCPRP')
      KFREE = 1
      LFREE = LWRK
!     Determine if the left and right states are allowed to be different
      numleftst = lstate(1)
      numrightst = rstate(1)
      DO I = 2,NBSYM
         numleftst = numleftst + lstate(I)
         numrightst = numrightst + rstate(I)
      END DO
      IF (numrightst .le. 0) THEN
         prponly = .TRUE.
      ELSE
         prponly = .FALSE.
      END IF
!     Print some info
      WRITE(LUPRI,'(/A)') '=================EXCPRP====================='
      WRITE(LUPRI,*) ' Calculation of transition moments for excited '//
     &     'states.',' In an exact state formalism, '//
     &     'these transition moments',
     &     ' are formally given as residues ',
     &     'of the response function.'
      WRITE(LUPRI,'(/A)') ' User input:'
      WRITE(LUPRI,'(I4,A)') NAQROP, ' operator(s): '
      WRITE(LUPRI,*) (PRPNAM(LAQROP(I)),I=1,NAQROP)
      IF (.not. prponly) THEN
         WRITE(LUPRI,'(I4,A)') numleftst, ' left state(s): '
         WRITE(LUPRI,'(8I4)') (lstate(I),I=1,NBSYM)
         WRITE(LUPRI,'(I4,A)') numrightst, ' right state(s): '
         WRITE(LUPRI,'(8I4)') (rstate(I),I=1,NBSYM)
      ELSE
         WRITE(LUPRI,'(I4,A)') numleftst, ' state(s): '
         WRITE(LUPRI,'(8I4)') (lstate(I),I=1,NBSYM)
      END IF
      WRITE(LUPRI,*) '==================EXCPRP====================='
      WRITE(LUPRI,'(///2A/)') 
     &     '*** ENTERING EXCITED STATE PROPERTIES ',
     &     'CALCULATION'
!
!     Get the MO coefficients
!
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL OPNFIL(LUCMO,'DFCOEF','OLD','PRPTPA')
      IF(SPINFR) THEN
        CALL MEMGET('INTE',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
        CALL REACMO(LUCMO,'DFCOEF',WORK(KCMO),DUM,WORK(KIBRP),TOTERG,11)
      ELSE
        KIBRP = KFREE
        CALL REACMO(LUCMO,'DFCOEF',WORK(KCMO),DUM,IDUM,TOTERG,3)
      ENDIF
      CLOSE(LUCMO,STATUS='KEEP')
!
!     Determine which symmetry classes that give contributions
      CALL EXCLISTSETUP(prponly,nexcsym)
      CALL ICOPY(8,nexcsym,1,nexcsim,1)
      CALL ICOPY(8,nexcsym,1,nexcstv,1)
!     Solve for the excitation vectors and energies in these symmetries
      CALL CALCEXCVEC(WORK(KFREE),LFREE,nexcsym,nexcsim,nexcstv)
!     Determine for which symmetries and frequencies to solve
!     the linear response equations
      CALL IZERO(tmeqlist,5*MAXQR)
      CALL FREQLISTSETUP(prponly, tmeqlist)

      CALL PRSYMB(LUPRI,'-',70,1)
      WRITE(LUPRI,*) 'EXCPRP: Number of linear response',
     &     ' eqs. to solve: ', NQRHYP
      CALL PRSYMB(LUPRI,'-',70,1)
      IF (NQRHYP .LE. 0) THEN
         CALL QUIT('No nontrivial response functions.')
      END IF
!     Solve linear response equations
      CALL PRSYMB(LUPRI,'=',70,1)
      WRITE(LUPRI,'(2A)') 
     &     ' >>>>>>>>              L I N E A R  R E S P O N S E',
     &     '            <<<<<<<<'
      CALL PRSYMB(LUPRI,'=',70,1)
!
      CALL QRVEC(WORK(KFREE),LFREE)
!     Calculate linear response function (not necessary, but done
!     because the computational cost is negligible)
      CALL MEMGET('REAL',KVECAEE,MZYEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECAEP,MZYEP*NZ,WORK,KFREE,LFREE)
!      IF (.not. prponly) THEN
         CALL LRCALC(RES_LR,WORK(KVECAEE),WORK(KVECAEP),WORK(KCMO),
     &        WORK(KIBRP),WORK(KFREE),LFREE)
!      END IF
!     
!     Calculate quadratic response function values
      WRITE(LUPRI,'(/A)') ' '
      CALL PRSYMB(LUPRI,'=',70,1)
      WRITE(LUPRI,'(2A)') 
     &     ' >>>>>>>>     E X C I T E D   S T A T E   ',
     &     'P R O P E R T I E S  <<<<<<<<'
      WRITE(LUPRI,'(2A)') 
     &     ' >>>>>>>>     ( A N D  M O M E N T S )    ',
     &     '                     <<<<<<<<'

      CALL PRSYMB(LUPRI,'=',70,1)
C     
      CALL MEMGET('REAL',KVECBEE,MZYEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECBEP,MZYEP*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECCEE,MZYEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECCEP,MZYEP*NZ,WORK,KFREE,LFREE)
      CALL EXCPRPCALC(prponly,tmeqlist,RES_EXCPRP,WORK(KVECAEE),
     &     WORK(KVECBEE),WORK(KVECCEE),WORK(KVECAEP),WORK(KVECBEP),
     &     WORK(KVECCEP),WORK(KCMO),WORK(KIBRP),
     &     WORK(KFREE),LFREE)
!     Print final results of the calculation
      CALL EXCPRP_PRINT(prponly,tmeqlist,RES_LR,RES_EXCPRP,WORK,LFREE)
      CALL QEXIT('EXCPRP')
      RETURN
      END
C
C
      SUBROUTINE EXCLISTSETUP(prponly,exclist)
C***********************************************************************
C
C     Routine for determining which excitation vectors
C     that should be computed.
C
C     Written by Erik Tellgren, aug 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
!
      LOGICAL prponly
      LOGICAL nonzeroflag
      INTEGER I,exclist(8),excsym
      INTEGER ifer,ibos,NFC
!
      CALL QENTER('EXCLISTSETUP')
!
      DO I=1,NBSYM
         exclist(I) = max(lstate(I),rstate(I))
      END DO
!
!     number of boson symmetries per fermion symmetry
      NFC = NBSYM / NFSYM
!     check if there's a matching right state and operator
      DO ifer=1,NFSYM
         DO ibos=1,NFC
!           store symmetry of the excited state
            excsym = JFSYM(ibos,ifer)
!           check if any non-zero transition moment is possible
            CALL EXCPRPSTCHK(prponly,excsym,nonzeroflag)
            IF (.not. nonzeroflag) THEN
               exclist(I) = 0
            END IF
         END DO
      END DO
      CALL QEXIT('EXCLISTSETUP')
      RETURN
      END
C
      SUBROUTINE EXCPRPCALC(prponly,tmomlist,RESULT,VECAEE,VECBEE,
     &     VECCEE,VECAEP,VECBEP,VECCEP,CMO,IBEIG,WORK,LWRK)
C***********************************************************************
C
C     Take the previously calculated excitation vectors and
C     response vectors, and assemble all the excited state properties or
C     matrix elements. The array RESULT will contain the results.
C
C     Written by Erik Tellgren, sep 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dcbgen.h"
#include "dgroup.h"
C
      REAL*8 D0
      PARAMETER ( D0=0.0D0 )
      LOGICAL prponly
      INTEGER LWRK,KFREE,LFREE
      INTEGER tmomlist(5,MAXQR)
      REAL*8 VECAEE(MZYEE,NZ),VECBEE(MZYEE,NZ),VECCEE(MZYEE,NZ),
     &          VECAEP(MZYEP,NZ),VECBEP(MZYEP,NZ),VECCEP(MZYEP,NZ)
      REAL*8 CMO(*),WORK(LWRK)
      INTEGER IBEIG(*)
      REAL*8 RESULT(3,NBSYM,MAXFQR,NBSYM,MAXFQR,MAXLQR)
      INTEGER ifer1,ifer2,ibos1,ibos2,iop
      INTEGER leftsym,rightsym,opbsym
      INTEGER NFC,IQRF,I,J,K,JF,JB
C
#include "ibtfun.h"
      MULD2H(I,J) = IBTXOR(I-1,J-1) + 1
C
      CALL QENTER('EXCPRPCALC')
C
      KFREE = 1
      LFREE = LWRK
C
C     Open the property file
C     ======================
C
      OPEN(LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &            FILE='AOPROPER')
C     
      CALL DZERO(RESULT,3*NBSYM*MAXFQR*NBSYM*MAXFQR*MAXLQR)
C
      NFC  = NBSYM/NFSYM
      DO K = 1,min(NQRHYP,MAXQR)
!        extract data from tmomlist
         iop = tmomlist(1,K)
         leftsym = tmomlist(2,K)
         I = tmomlist(3,K)
         rightsym = tmomlist(4,K)
         J = tmomlist(5,K)
!        determine fermion symmetries
         DO JF = 1,NFSYM
            DO JB = 1,NFC
               IF (leftsym .eq. JFSYM(JB,JF)) THEN
                  ifer1 = JF
                  ibos1 = JB
               END IF
               IF (rightsym .eq. JFSYM(JB,JF)) THEN
                  ifer2 = JF
                  ibos2 = JB
               END IF
            END DO
         END DO
!        determine operator symmetry
         opbsym = IPRPSYM(LAQROP(iop))
!        assemble transition moment
         IQRF = K
         CALL EXCPRP_TMASS(IQRF,
     &        RESULT(1,leftsym,I,rightsym,J,iop),
     &        ibos1,ifer1,
     &        ibos2,ifer2,iop,I,J,VECAEE,VECBEE,
     &        VECCEE,VECAEP,VECBEP,VECCEP,CMO,
     &        IBEIG,WORK(KFREE),LFREE)
      END DO
C
C     Close the property file
C     =======================
C
      CLOSE(LU1INT,STATUS='KEEP')
C
      CALL QEXIT('EXCPRPCALC')
      RETURN
      END
C
      SUBROUTINE EXCPRP_TMASS(IQRF,RESULT,ibos1,ifer1,ibos2,ifer2,iop,
     &     exc1,exc2,VECAEE,VECBEE,VECCEE,VECAEP,
     &     VECBEP,VECCEP,CMO,IBEIG,WORK,LWRK)
C***********************************************************************
C
C     Routine for assembling an excited state property/transition moment
C     from stored excitation and response vectors.
C
C     Written by Erik Tellgren, aug 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dcbxrs.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
C
      REAL*8 D1
      PARAMETER ( D1=1.0D0 )
      INTEGER IQRF,ibos1,ifer1,ibos2,ifer2,iop,exc1,exc2,LWRK
      INTEGER IBEIG(*)
      REAL*8 RESULT(3)
      REAL*8 VECAEE(MZYEE,NZ),VECBEE(MZYEE,NZ),VECCEE(MZYEE,NZ),
     &          VECAEP(MZYEP,NZ),VECBEP(MZYEP,NZ),VECCEP(MZYEP,NZ)
      REAL*8 CMO(*),WORK(LWRK)
      INTEGER KFREE, LFREE
      CHARACTER optag*8, lefttag*8, righttag*8
      INTEGER leftsym,rightsym,indlexc,indrexc,itim1,itim2
      INTEGER opbsym,opfsym,optim,oppointer
      REAL*8 freq1,freq2,lrfreq
      INTEGER NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP
!
      KFREE = 1
      LFREE = LWRK
!
      RESULT(1) = 0.0
      RESULT(2) = 0.0
      RESULT(3) = 0.0
!     determine symmetry properties and indices to excitation eqs.
      leftsym = JFSYM(ibos1,ifer1)
      indlexc = IEXCOFF + leftsym
      rightsym = JFSYM(ibos2,ifer2)
      indrexc = IEXCOFF + rightsym
      opbsym = IPRPSYM(LAQROP(iop))
!     read left excited state from file (store as "B vector")
      CALL EXCINIT(indlexc,exc1,freq1,leftsym,ifer1,
     &     itim1,lefttag,NZYBEE,NZYBEP,VECBEE,VECBEP,
     &     WORK(KFREE),LFREE)
!     read right excited state from file (store as "C vector")
      CALL EXCINIT(indrexc,exc2,freq2,rightsym,ifer2,
     &     itim2,righttag,NZYCEE,NZYCEP,VECCEE,
     &     VECCEP,WORK(KFREE),LFREE)
!     read response vector from file (store as "A vector")
      CALL HYPINIT(IQRF,1,lrfreq,oppointer,opbsym,
     &     opfsym,optim,optag,NZYAEE,NZYAEP,VECAEE,
     &     VECAEP,WORK(KFREE),LFREE)
!
      WRITE(LUPRI,'(/A)') 'Assembling excited state '//
     &     'property/moment from:'
      WRITE(LUPRI,*) ' Response vector ',
     &     IQRF,' out of ',NQRHYP
      WRITE(LUPRI,'(A,A12,I5,F16.8)') 
     &     ' A operator, boson symmetry, frequency: ',
     &     optag,opbsym,lrfreq
      WRITE(LUPRI,'(A,A,I5,F16.8)')
     &     ' Left state, boson symmetry, frequency: ',
     &     lefttag,leftsym,freq1
      WRITE(LUPRI,'(A,A,I5,F16.8)')
     & ' Right state, boson symmetry, frequency: ',
     &     righttag,rightsym,freq2
      CALL FLSHFO(LUPRI)
C     -----------------------
C     Assemble matrix element
C     -----------------------
!     add Nb_j A[2]_jk Nc_k to RESULT
!     the Nb_j vector should be transformed by QSWAP
      CALL QSWAP(NZYBEE,VECBEE,leftsym,D1,.true.)
      CALL QSWAP(NZYBEP,VECBEP,leftsym,D1,.true.)
      CALL A2DRV(oppointer,VECBEE,VECBEP,VECCEE,
     &     VECCEP,opbsym,leftsym,rightsym,optim,
     &     itim1,itim2,opfsym,ifer1,ifer2,
     &     NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,
     &     RESULT,CMO,IBEIG,WORK(KFREE),LFREE)
!     add Nc_j A[2]_jk Nb_k to RESULT
      CALL A2DRV(oppointer,VECCEE,VECCEP,VECBEE,
     &     VECBEP,opbsym,rightsym,leftsym,optim,
     &     itim2,itim1,opfsym,ifer2,ifer1,
     &     NZYAEE,NZYAEP,NZYCEE,NZYCEP,NZYBEE,NZYBEP,
     &     RESULT,CMO,IBEIG,WORK(KFREE),LFREE)
!     add Na_j (E[3]_jkl + E[3]_jlk + ...) Nc_k Nb_l
      CALL T3DRV(VECAEE,VECAEP,VECCEE,VECCEP,
     &     VECBEE,VECBEP,opbsym,rightsym,leftsym,
     &     optim,itim2,itim1,opfsym,ifer2,ifer1,
     &     NZYAEE,NZYAEP,NZYCEE,NZYCEP,
     &     NZYBEE,NZYBEP,freq2,-freq1,
     &     RESULT,CMO,WORK(KFREE),LFREE)
      WRITE(LUPRI,'(A,E16.8)') ' Value of electronic part:', RESULT(2)
      WRITE(LUPRI,'(A,E16.8)') ' Value of positronic part:', RESULT(3)
      WRITE(LUPRI,'(A,E16.8)') ' Value of total moment   :', RESULT(1)
      WRITE(LUPRI,'(A,A)') '(N.B. <g|op|f> - delta_gf <0|op|0>, not',
     &     ' <g|op|f>, is given here.)'
      RETURN
      END
!
      SUBROUTINE EXCPRPSTCHK(prponly,excstsym,stateOK)
C***********************************************************************
C
C     Routine for checking if an excited state will contribute
C     to any non-zero transition moments.
C
C     Written by Erik Tellgren, aug 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dcbxrs.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
C
      LOGICAL stateOK, prponly
      INTEGER excstsym
      INTEGER I,iop,ifer,ibos
      INTEGER lstsym,rstsym
      INTEGER NFC
#include "ibtfun.h"
      MULD2H(I,J) = IBTXOR(I-1,J-1) + 1
!     number of boson symmetries per fermion symmetry
      NFC = NBSYM / NFSYM
!     are different left and right states allowed?
      IF (.NOT. prponly) THEN
!        check if there's a matching right state and operator
         DO ifer=1,NFSYM
            DO ibos=1,NFC
!              store left state symmetry
               rstsym = JFSYM(ibos,ifer)
               IF (rstate(rstsym) .gt. 0) THEN
!                 loop over operators
                  DO iop=1,NAQROP
                     IF (MULD2H(rstsym,MULD2H(IPRPSYM(LAQROP(iop)),
     &                    excstsym)) .eq. 1) THEN
                        stateOK = .true.
                        RETURN
                     END IF
                  END DO
               END IF
            END DO
         END DO
!        check if there's a matching left state and operator
         DO ifer=1,NFSYM
            DO ibos=1,NFC
!              store left state symmetry
               lstsym = JFSYM(ibos,ifer)
               IF (lstate(lstsym) .gt. 0) THEN
!                 loop over operators
                  DO iop=1,NAQROP
                     IF (MULD2H(lstsym,MULD2H(IPRPSYM(LAQROP(iop)),
     &                    excstsym)) .eq. 1) THEN
                        stateOK = .true.
                        RETURN
                     END IF
                  END DO
               END IF
            END DO
         END DO
!     identical left and right states should be used
      ELSE
!        check if there's a matching operator
         DO iop=1,NAQROP
            IF (IPRPSYM(LAQROP(iop)) .eq. 1) THEN
               stateOK = .true.
               RETURN
            END IF
         END DO
      END IF
      stateOK = .false.
      RETURN
      END
!
!
      SUBROUTINE EXCPRPINP(KEYWID, DOFLAG, SETTINGS)
C*****************************************************************************
C
C     Input handling for excited state properties calculation
C
C     Written by Erik Tellgren, aug 2004
C
C     Input: pointers DOFLAG and SETTINGS
C            KEYWID indicates which keyword or request that
C            triggered the call to EXCPRPINP
C            (1-.EXCDIP, 2-.EXCPRP, 4-consistency check)
C     Output: Stored in common var. lstate and rstate,
C             info about the calculation is stored in SETTINGS
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
!
      REAL*8 D1
      PARAMETER ( D1=1.0D0 )
      INTEGER KEYWID, SETTINGS, LWRK
      LOGICAL DOFLAG
      LOGICAL lstateOK
      INTEGER LFREE, KFREE, DUMMY, DUMMY2, OPINDEX
      INTEGER I
      CHARACTER OPNAME*16
!
      CALL QENTER('EXCPRPINP')
!     Do consistency check if KEYWID == 4
      IF (KEYWID .eq. 4) THEN
         DOFLAG = .false.
         IF ((SETTINGS.ne.1).and.(SETTINGS.ne.2)) THEN
            CALL QEXIT('EXCPRPINP')
            RETURN
         END IF
         ! check for B and C operators and frequencies
         IF ((NBQROP .gt. 0) .or. (NCQROP .gt. 0)) THEN
            NBQROP = 0
            NCQROP = 0
            WRITE(LUPRI,*) 'EXCPRPINP: Warning, the specified B and C ',
     &           'operators will be ignored.'
         END IF
         IF ((NBQRFR .gt. 0) .or. (NCQRFR .gt. 0)) THEN
            NBQRFR = 0
            NCQRFR = 0
            WRITE(LUPRI,*) 'EXCPRPINP: Warning, the specified B and C ',
     &           'frequencies will be ignored.'
         END IF
         ! check that an operator has been defined
         IF (NAQROP .lt. 1) THEN
            WRITE(LUPRI,*) 'EXCPRPINP: Error, no operator has been ',
     &           'specified for the excited state property calculation.'
            CALL QUIT('EXCPRPINP: Error, no operator has been '//
     &           'specified for the excited state property '//
     &           'calculation.')
            RETURN
         END IF
         ! if SETTINGS==1, check that the dipole op is selected
         IF ((SETTINGS .eq. 1) .and. (NAQROP .ne. 3)) THEN
            WRITE(LUPRI,*) 'EXCPRPINP: Error, only the dipole ',
     &           'operator should be used.'
            CALL QUIT('EXCPRPINP: Error, only the dipole '//
     &           'operator should be used.')
            RETURN
         END IF
         DOFLAG = .true.
         CALL QEXIT('EXCPRPINP')
         RETURN
      END IF

!     Check if SETTINGS has already been set
      IF (SETTINGS .ne. 0) THEN
         WRITE(LUPRI,*) 'EXCPRPINP: Error, multiple specifications ',
     &        'of excited state property calculations.'
         CALL QUIT('EXCPRPINP: Error, multiple specifications '//
     &        'of excited state property calculations.')
         RETURN
      END IF
!     Check if NBSYM > 8
      IF (NBSYM .gt. 8) THEN
         WRITE(LUPRI,*) 'EXCPRPINP: Error, too large symmetry ',
     &        'group (NBSYM > 8).'
         CALL QUIT('EXCPRPINP: Error, too large symmetry '//
     &        'group (NBSYM > 8).')
         RETURN
      END IF
!
!     initialize left and right states
      CALL IZERO(lstate,8)
      CALL IZERO(rstate,8)
!     keyword: .EXCDIP
      IF (KEYWID .eq. 1) THEN
         SETTINGS = 1
!        read excited states (only left states are read, leaving
!        rstate to contain only zeros; this is as it should be)
         READ(LUCMD,*,END=999,ERR=999) (lstate(I),I=1,NBSYM)
!        set the A operators to the x,y,z-dipole operators
         OPNAME = 'XDIPLEN'
         CALL XPRIND(OPNAME,1,1,D1,OPNAME(1:8),OPINDEX,DUMMY,
     &        DUMMY2,IPRXQR)
         CALL OP1IND('NAQROP',DUMMY,LAQROP,NAQROP,OPINDEX,MAXLQR)
         OPNAME = 'YDIPLEN'
         CALL XPRIND(OPNAME,1,1,D1,OPNAME(1:8),OPINDEX,DUMMY,
     &        DUMMY2,IPRXQR)
         CALL OP1IND('NAQROP',DUMMY,LAQROP,NAQROP,OPINDEX,MAXLQR)
         OPNAME = 'ZDIPLEN'
         CALL XPRIND(OPNAME,1,1,D1,OPNAME(1:8),OPINDEX,DUMMY,
     &        DUMMY2,IPRXQR)
         CALL OP1IND('NAQROP',DUMMY,LAQROP,NAQROP,OPINDEX,MAXLQR)
!     keyword: .EXCPRP
      ELSE IF (KEYWID .eq. 2) THEN
         SETTINGS = 2
!        Read left excited states
         lstateOK = .false.
         READ(LUCMD,*,END=999,ERR=999) (lstate(I),I=1,NBSYM)
         lstateOK = .true.
!        Read right excited states
         READ(LUCMD,*,END=999,ERR=999) (rstate(I),I=1,NBSYM)
      END IF
!
      CALL QEXIT('EXCPRPINP')
      RETURN
C     Handle erroneous input
 999  CONTINUE
      IF (KEYWID .eq. 1) THEN
         IF (lstateOK) THEN
            WRITE(LUPRI,*) 'EXCPRPINP: Error while reading right',
     &           'states (',NBSYM,' integers expected).'
         ELSE
            WRITE(LUPRI,*) 'EXCPRPINP: Error while reading left',
     &        'states (',NBSYM,' integers expected).'
         END IF
         CALL QUIT('EXCPRPINP: Error, '//
     &        'incorrect input data following .EXCPRP.')
      ELSE IF (KEYWID .eq. 2) THEN
         WRITE(LUPRI,*) 'EXCPRPINP: Error while reading',
     &        'states (',NBSYM,' integers expected).'
         CALL QUIT('EXCPRPINP: Error, '//
     &        'incorrect input data following .EXCDIP.')
      END IF
      RETURN
      END
!
!
      SUBROUTINE ADDTOFREQLIST(tmomlist,iop,sym1,exc1,sym2,exc2,numeqs)
C***********************************************************************
C
C     Routine for determining which linear response equations
C     that should be solved. It also prepares for the actual solving
!     by setting certain common blocked variables.
C
C     Written by Erik Tellgren, aug 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxpr.h"
#include "dcbxqr.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
!
      INTEGER tmomlist(5,MAXQR)
      INTEGER iop,sym1,sym2,exc1,exc2,numeqs
      INTEGER leftsym,rightsym,leftexc,rightexc
      INTEGER K
      LOGICAL unique
!
#include "ibtfun.h"
      MULD2H(I,J) = IBTXOR(I-1,J-1) + 1
!
      CALL QENTER('ADDTOFREQLIST')
!     copy some data
      leftsym = sym1
      leftexc = exc1
      rightsym = sym2
      rightexc = exc2
      IF (exc1 .gt. exc2) THEN
         leftsym = sym2
         leftexc = exc2
         rightsym = sym1
         rightexc = exc1
      ELSEIF ((exc1 .eq. exc2) .and. (sym1 .gt. sym2)) THEN
         leftsym = sym2
         leftexc = exc2
         rightsym = sym1
         rightexc = exc1         
      ENDIF
!     check if there already is an equivalent entry in the list
      unique = .true.
      DO K=1,numeqs
         IF ((tmomlist(1,K) .eq. iop) .and.
     &        (tmomlist(2,K) .eq. leftsym) .and.
     &        (tmomlist(3,K) .eq. leftexc) .and.
     &        (tmomlist(4,K) .eq. rightsym) .and.
     &        (tmomlist(5,K) .eq. rightexc)) THEN
            unique = .false.
         END IF
      END DO
!     add entry to list, if unique
      IF (unique .and. (numeqs .lt. MAXQR)) THEN
         numeqs = numeqs + 1
         tmomlist(1,numeqs) = iop
         tmomlist(2,numeqs) = leftsym
         tmomlist(3,numeqs) = leftexc
         tmomlist(4,numeqs) = rightsym
         tmomlist(5,numeqs) = rightexc
      END IF
      CALL QEXIT('ADDTOFREQLIST')
      RETURN
      END
C
      SUBROUTINE FREQLISTSETUP(prponly, tmomlist)
C***********************************************************************
C
C     Routine for determining which linear response equations
C     that should be solved. It also prepares for the actual solving
!     by setting certain common blocked variables.
C
C     Written by Erik Tellgren, aug 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbxpr.h"
#include "dcbxqr.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
!
      LOGICAL prponly, toomanyeqs
      INTEGER tmomlist(5,MAXQR)
      INTEGER ifer1,ifer2,ibos1,ibos2,iop
      INTEGER leftsym,rightsym,opbsym,oppointer
      INTEGER recnum, NFC, I, J, K, numeqs
      REAL*8 freq1,freq2,lrfreq
!
#include "ibtfun.h"
      MULD2H(I,J) = IBTXOR(I-1,J-1) + 1
!
      CALL QENTER('FREQLISTSETUP')
      toomanyeqs = .false.
      numeqs = 0
!     number of boson symmetries per fermion symmetry
      NFC = NBSYM / NFSYM
!     check if the the left and right states must be the same
      IF (prponly) THEN
         DO ifer1=1,NFSYM
         DO ibos1=1,NFC
            leftsym = JFSYM(ibos1,ifer1)
            DO iop = 1,NAQROP
               oppointer = LAQROP(iop)
               opbsym = IPRPSYM(oppointer)
               IF (MULD2H(leftsym,MULD2H(opbsym,leftsym)).eq.1) THEN
!                 loop over states in the symmetry classes
                  DO I = 1,lstate(leftsym)
                     IF (numeqs .lt. MAXQR) THEN
                        rightsym = leftsym
                        J = I
                        CALL ADDTOFREQLIST(tmomlist,iop,leftsym,I,
     &                       rightsym,J,numeqs)
                     ELSE
                        numeqs = MAXQR
                        toomanyeqs = .true.
                     END IF
                  END DO
               END IF
            END DO
         END DO
         END DO
      ELSE
!        Different left and right states are allowed.
!        loop over left state fermion and boson symmetries
         DO ifer1=1,NFSYM
         DO ibos1=1,NFC
!        loop over right state fermion and boson symmetries
         DO ifer2=1,NFSYM
         DO ibos2=1,NFC
            leftsym = JFSYM(ibos1,ifer1)
            rightsym = JFSYM(ibos2,ifer2)
!           loop over operators
            DO iop = 1,NAQROP
               oppointer = LAQROP(iop)
               opbsym = IPRPSYM(oppointer)
               IF (MULD2H(leftsym,MULD2H(opbsym,rightsym)).eq.1) THEN
!                 loop over states in the symmetry classes
                  DO I = 1,lstate(leftsym)
                  DO J = 1,rstate(rightsym)
                     IF (numeqs .lt. MAXQR) THEN
                        CALL ADDTOFREQLIST(tmomlist,iop,leftsym,I,
     &                       rightsym,J,numeqs)
                     ELSE
                        numeqs = MAXQR
                        toomanyeqs = .true.
                     END IF
                  END DO
                  END DO
               END IF
            END DO
         END DO
         END DO
         END DO
         END DO
      END IF
!     Print warning if too many transition moments were requested
      IF (toomanyeqs) THEN
         WRITE(LUPRI,*) '**** Warning! Too many linear ',
     &        'response equations specified.',
     &        'Maximum is ',MAXQR
         WRITE(LUPRI,*) 'The calculation will continue, but some ',
     &        'requested transition moments will be omitted.'
      END IF
!     Use INDQR to insert linear response equation into a list
!     of equations to be solved
      DO K = 1,numeqs
!        extract transition moment specification
         iop = tmomlist(1,K)
         oppointer = LAQROP(iop)
         leftsym = tmomlist(2,K)
         I = tmomlist(3,K)
         rightsym = tmomlist(4,K)
         J = tmomlist(5,K)
!        insert into list
         freq1 = EXCIT(I,leftsym)
         freq2 = EXCIT(J,rightsym)
         lrfreq = freq2 - freq1
         CALL INDQR(oppointer,abs(lrfreq),recnum)
!        set up other lists
         LQRHYP(K,1) = oppointer
         QRFRHYP(K,1) = -lrfreq
         QRFRHYP(K,2) = -lrfreq
      END DO
!
      NQRHYP = numeqs
!
      CALL QEXIT('FREQLISTSETUP')
      RETURN
      END
C
      SUBROUTINE EXCPRP_PRINT(prponly,tmomlist,RES_LR,RES_EXCPRP,
     &     WORK,LWRK)
C***********************************************************************
C
C     Printing of results aquired during EXCPRPCALC.
C
C     Written by Erik Tellgren
C
C***********************************************************************
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "dcbxqr.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "pgroup.h"
C
      PARAMETER ( LUCMO = 22, D0 = 0.0D0, D2R3 = (2.0D0/3.0D0),
     &            D2 = 2.0D0, D4 = 4.0D0 )
      LOGICAL prponly
      INTEGER tmomlist(5,MAXQR)
      CHARACTER opname*8,opname2*8,excopname*8,PFMT*6
      REAL*8 RES_LR(3,MAXLQR,MAXLQR,MAXFQR)
      REAL*8 RES_EXCPRP(3,NBSYM,MAXFQR,NBSYM,MAXFQR,MAXLQR)
      REAL*8 WORK(*)
      INTEGER LWRK, KFREE, LFREE, IQRF, NFC, cnt
      INTEGER ifer1,ifer2,ibos1,ibos2,leftsym,rightsym
      INTEGER EXCOPFROMSYM
      INTEGER iop,iop2,opsym,opsym2,opptr,opptr2,
     &     excop,excsym
      INTEGER I,J,K,NFR
      REAL*8 freq1,freq,rsp,oscstr
C
#include "ibtfun.h"
      MULD2H(I,J) = IBTXOR(I-1,J-1) + 1
C
      CALL QENTER('EXCPRP_PRINT')
C
      KFREE = 1
      LFREE = LWRK
C
      CALL TITLER('RESULTS OF EXCITED STATE PROPERTIES CALCULATION',
     &            '*',116)
      CALL HEADER('Linear response properties',-1)
C
C     ***********************************
C     *** Print excitation properties ***
C     ***********************************
C
      CALL HEADER('Excitation properties',0)
C
      NFC = NBSYM/NFSYM
C
C     Loop over excitation symmetries
C     ===============================
C
      DO ifer1 = 1,NFSYM
      DO ibos1 = 1,NFC
         excsym = JFSYM(ibos1,ifer1)
         CALL PRSYMB(LUPRI,'-',80,0)
         WRITE(LUPRI,'(A,I3,A3,A)') ' *** Excited state boson symmetry',
     &        excsym,' : ',REP(excsym-1)
         IF (prponly) THEN
            WRITE(LUPRI,'(A,I3)')
     &           'No. of excitations in the symmetry: ',lstate(excsym)
         ELSE
            WRITE(LUPRI,'(A,I3,A7,I3,A8)')
     &           'No. of excitations in the symmetry: ',
     &           lstate(excsym), ' (left), ',
     &           rstate(excsym), ' (right)'
         END IF
C        loop over excitations (same left and right state)
         DO I = 1,max(lstate(excsym),rstate(excsym))
            freq1 = EXCIT(I,excsym)
            WRITE(LUPRI,'(A,I3,A,F14.6,A)') '* Excitation energy no. ',
     &           I,' : ',freq1, ' a.u.'
            WRITE(LUPRI,'(F14.4,A)')
     &           freq1*XTEV, ' eV'
C           determine oscillator strength
            oscstr = D0
            DO iop2 = 1,NQROP
               opptr2 = LQROP(iop2,1)
               opname2  = PRPNAM(opptr2)(1:8)
               opsym2 = IPRPSYM(opptr2)
C     check that this is an A operator
               IF (opptr2 .lt. IEXCOFF) THEN
                  IF (MULD2H(opsym2,excsym) .EQ. 1) THEN
                     excop = EXCOPFROMSYM(excsym)
                     rsp = RES_LR(1,iop2,excop,I)
                     WRITE(LUPRI,'(/A,A10,F14.6,A)')
     &                    ' Operator, transition moment: ',
     &                    opname2,rsp,' a.u.'
                     oscstr = oscstr + D2R3*rsp*rsp*freq1
                  END IF
               END IF
C     end loop over A operators
            END DO
C     print oscillator strength
            IF (oscstr .ge. 1.0D-10) THEN
               PFMT = 'F14.6'   !MXFORM(oscstr,12)
               WRITE(LUPRI,'(/5X,A,'//PFMT//')')
     &              ' Oscillator strength: ',oscstr
            ELSE
               WRITE(LUPRI,'(/5X,A)')
     &              'Forbidden transition'
            END IF
C        end loop over excitation frequencies
         END DO
      END DO
      END DO
      CALL PRSYMB(LUPRI,'=',80,0)
      CALL FLSHFO(LUPRI)
C
C     *****************************
C     *** Print linear response ***
C     *****************************
C     
      IF (.not. prponly) THEN
         CALL HEADER('Linear response',0)
         CALL PRSYMB(LUPRI,'-',80,0)
         WRITE(LUPRI,*)
     &        ' Oper. symmetry    A-operator    B-operator    '//
     &        'Frequency [a.u]      Response [a.u.]'
         CALL PRSYMB(LUPRI,'-',80,0)
C
C     Loop over A- and B-operators
C     ============================
C     Check that neither A- nor B-operator is excitation operator and
C     check if symmetry allowed.
C
         DO iop = 1,NQROP
            opptr = LQROP(iop,1)
            opname  = PRPNAM(opptr)(1:8)
            opsym = IPRPSYM(opptr)
            DO iop2 = iop,NQROP
               opptr2  = LQROP(iop2,1)
               opname2  = PRPNAM(opptr2)(1:8)
               opsym2 = IPRPSYM(opptr2)
!              added conjunct to IF statement / ET 2005-03-29
               IF ((opptr .lt. IEXCOFF) .and.
     &              (opptr2 .lt. IEXCOFF) .and.
     &              (MULD2H(opsym,opsym2) .EQ. 1)) THEN
!                 Loop over frequencies and print
                  DO I = 1,min(LQROP(iop,2),MAXFQR)
                     freq = QRFREQ(I,iop)
                     rsp  = RES_LR(1,iop2,iop,I)
                     WRITE(LUPRI,
     &                    '(7X,A3,10X,A8,6X,A8,6X,F14.6,4X,F14.6)')
     &                    REP(opsym2-1),opname2,opname,freq,rsp
                  END DO
               END IF
            END DO
         END DO
         CALL PRSYMB(LUPRI,'-',80,0)
      END IF
C     
C     ****************************
C     *** Print EXCPRP results ***
C     ****************************
C
      CALL HEADER('Quadratic response properties',-1)
C
      CALL PRSYMB(LUPRI,'-',80,0)
      WRITE(LUPRI,*) 'Excited state properties and moments that are',
     &     ' not trivially zero due to symmetry follows:'
      CALL PRSYMB(LUPRI,'-',80,0)
      WRITE(LUPRI,*)
      DO cnt = 1,min(NQRHYP,MAXQR)
!        extract data from tmomlist
         iop = tmomlist(1,cnt)
         leftsym = tmomlist(2,cnt)
         I = tmomlist(3,cnt)
         rightsym = tmomlist(4,cnt)
         J = tmomlist(5,cnt)
!        extract more operator data
         opptr = LAQROP(iop)
         opsym = IPRPSYM(opptr)
!        print transition moment
         IF (leftsym.eq.rightsym.and.I.eq.J) THEN
            WRITE(LUPRI,'(I3,A,F14.6)') cnt,' Energy: ',EXCIT(I,leftsym)
            WRITE(LUPRI,'(A2,I3,A1,A3,A2,A9,A1,I3,A1,'//
     &           'A3,A8,A9,A6,F14.6)')
     &           '< ',I,' ',REP(leftsym-1),
     &           '| ',PRPNAM(opptr)(1:8),'|',I,' ',
     &           REP(rightsym-1),'> - <0|',
     &           PRPNAM(opptr)(1:8),'|0> = ',
     &           RES_EXCPRP(1,leftsym,I,rightsym,J,iop)
         ELSE
            WRITE(LUPRI,'(I3,A,F14.6,F14.6)')
     &           cnt,' Left energy, right energy: ',
     &           EXCIT(I,leftsym),EXCIT(J,rightsym)
            WRITE(LUPRI,'(A2,I3,A1,A3,A2,A9,A1,I3,A1,A3,A4,F14.6)')
     &           '<',I,' ',REP(leftsym-1),
     &           '| ',PRPNAM(opptr)(1:8),'|',J,' ',
     &           REP(rightsym-1),'> = ',
     &           RES_EXCPRP(1,leftsym,I,rightsym,J,iop)
         END IF
         WRITE(LUPRI,*)
      END DO
!     also print the transition moments that are zero due to symmetry
      CALL PRSYMB(LUPRI,'-',80,0)
      WRITE(LUPRI,*) 'Excited state properties and moments that are',
     &     ' trivially zero due to symmetry follows:'
      CALL PRSYMB(LUPRI,'-',80,0)
      WRITE(LUPRI,*)
      NFC  = NBSYM/NFSYM
!     check if the left and right states must be the same
      IF (prponly) THEN
!        loop over left states (used also for the right states)
         DO ifer1=1,NFSYM
         DO ibos1=1,NFC
            ifer2 = ifer1
            ibos2 = ibos1
            leftsym = JFSYM(ibos1,ifer1)
            rightsym = leftsym
!           loop over operators
            DO iop = 1,NAQROP
               opptr = LAQROP(iop)
               opsym = IPRPSYM(opptr)
               IF (MULD2H(leftsym,MULD2H(opsym,rightsym)).ne.1 .and.
     &              (lstate(leftsym) .gt. 0)) THEN
                  WRITE(LUPRI,*) '< * ',REP(leftsym-1),
     &                 '| ',PRPNAM(opptr)(1:8),'| * ',
     &                 REP(rightsym-1),'> = ',
     &                 '0     (by symmetry)'
                  WRITE(LUPRI,*)
               END IF
               CALL FLSHFO(LUPRI)
            END DO
         END DO
         END DO
         CALL QEXIT('EXCPRP_PRINT')
         RETURN         
      END IF
!     --------------------------------------------
!     Different left and right states are allowed.
!     --------------------------------------------
!     loop over left state fermion and boson symmetries
      DO ifer1=1,NFSYM
      DO ibos1=1,NFC
         leftsym = JFSYM(ibos1,ifer1)
!        loop over right state fermion and boson symmetries
         DO ifer2=1,NFSYM
         DO ibos2=1,NFC
            rightsym = JFSYM(ibos2,ifer2)
!           loop over operators
            DO iop = 1,NAQROP
               opptr = LAQROP(iop)
               opsym = IPRPSYM(opptr)
               IF (MULD2H(leftsym,MULD2H(opsym,rightsym)).ne.1 .and.
     &              (lstate(leftsym) .gt. 0) .and.
     &              (rstate(rightsym) .gt. 0)) THEN
                  WRITE(LUPRI,*) '< * ',REP(leftsym-1),
     &                 '| ',PRPNAM(opptr)(1:8),'| * ',
     &                 REP(rightsym-1),'> = ',
     &                 '0   (by symmetry)'
                  WRITE(LUPRI,*)
               END IF
            END DO
         END DO
         END DO
      END DO
      END DO
      CALL PRSYMB(LUPRI,'*',80,0)
C
      CALL QEXIT('EXCPRP_PRINT')
      RETURN
      END
C
      INTEGER FUNCTION EXCOPFROMSYM(sym)
C***********************************************************************
C
C     This function determines the index to the excitation operator
C     (it is assumed that there is only one) that has symmetry SYM.
C     This function is needed here because RES_LR is indexed by operator
C     indices and not symmetries.
C
C     Written by Erik Tellgren
C
C***********************************************************************
      IMPLICIT NONE
#include "dcbxqr.h"
C
      integer sym
      integer excop, excopptr, excopsym
      do excop = 1,NQROP
         excopptr = LQROP(excop,1)
         if (excopptr .ge. IEXCOFF) then
            excopsym = excopptr - IEXCOFF
            if (excopsym .eq. sym) then
               EXCOPFROMSYM = excop
               return
            end if
         end if
      end do

      CALL QUIT('EXCOPFROMSYM: No exc. operator found.')
      end

      
