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

C FILE    : dirset.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamset */
      SUBROUTINE PAMSET()
C***********************************************************************
C
C     Subroutine for setting dimensions and generate overlap and Lowdin
C     matrices
C
C     Called from:  DIRCTL (after PAMINP)
C                   EXEDIT ( -||-) - optimization
C
C***********************************************************************
      use memory_allocator
      use visual_cfg
      use x2cmod_cfg
#ifdef MOD_CAP
      use memory
#endif

#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "dcbham.h"
#include "dcblab.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "cbihr1.h"
#include "cbirea.h"
#include "infpar.h"

#ifdef BUILD_GEN1INT
#include "gen1int.h"
#endif
C
      LOGICAL TWOBUF,GOFAST,TOBE
      real(8) :: TIMSTR(2), TIMEND(2)
      real(8), allocatable :: WORK(:)

      CALL QENTER('PAMSET')
      call legacy_lwork_get(LWORK)
      call alloc(WORK,LWORK,id='WORK in PAMSET - 1')

      CALL TIMER2('START ',TIMSTR,TIMEND)
C
C     Attempt a restart unless geometry optimization
C     ==============================================
C
      IF(OPTIMI) THEN
        GOFAST = .FALSE.
      ELSE
        CALL GOTEST(GOFAST,IPREAD)
      ENDIF
C
C ... if 2c-mode IOTC, read basis in 4c-mode here
C
      IF (TWOCOMPBSS .or. x2c) THEN
        TWOCOMP = .FALSE.
        RDINPC = .FALSE.
        CNTMAT = .FALSE.
!       IPREAD = IPREAD - 5
        CALL READIN(.FALSE.)
        CALL SETDC1(IPREAD)
!       IPREAD = IPREAD + 5
      ENDIF
C
C     Generate all one-electron integrals
C     ===================================
C     exit if DOHRM: .ONLY INTEGRALS -- only calc. 1-el. integrals
C
C     Start RECP integral     
C     ===================
      IF (ECPCALC) CALL RECP

#ifdef BUILD_GEN1INT
!     test suite of gen1int interface
      if (test_gen1int) then
         call gen1int_host_test(lwork, work, lupri, iprdef)
      end if
#endif

      INQUIRE(FILE='AOPROPER',EXIST=TOBE)
!     radovan: don't restart if you do a geometry optimization!
!     gosia: don't restart and don't generate integrals if visual_cfg_london
      if (tobe .and. visual_cfg_london) then
        goto 11
      endif
      IF(TOBE .and. .not. optimi) THEN
        WRITE(LUPRI,'(A)')
     &  '* WARNING: DIRAC found the AOPROPER file '//
     &  '(one-electron integrals) ...will attempt restart.'
        IF (NOSMLV) CALL NOSSNUCATT
      ELSE 
        CALL ONEGEN(WORK,LWORK)
      ENDIF
 11   continue
#ifdef MOD_CAP
        IF (CAP) CALL CAP_STUFF(WORK,LWORK)
#endif

      IF(DOHRM) THEN
        CALL TIMER2('PAMSET',TIMSTR,TIMEND)
        RETURN
      ENDIF
C
C     Generate MO-transformation matrix and write to file
C     ===================================================
C
C
C     DO4C2C=T : after 4c SCF perform IOTC transformation
C     DO2C4C=T : jump to 4c level after some 2c SCF iterations
C
      call dealloc(work) ! release memory - particularly because of X2C...

      CALL GMOTRA(GOFAST)
C
      call alloc(WORK,LWORK,id='WORK in PAMSET - 2')
C
C     Set DCB common blocks - second time
C     ===================================
C
      CALL SETDC2(IPRGEN)
C
C     If semiclassical QED, write density of bare vacuum to file
C     ==========================================================
C
#ifdef MOD_UNRELEASED
      IF(VACPOL.AND.VACREF) CALL BARVAC(WORK,LWORK,IPRGEN)      
#endif
C
C     Prepare for two-electron integrals
C     ==================================
C
      IF(.NOT.ONESYS) CALL TWOPRP(WORK,LWORK)
C
C     Prepare for parallel calculation (set inforb)
C     ==============================================
C
      CALL RELINF
C
      CALL TIMER2('PAMSET',TIMSTR,TIMEND)
      call dealloc(WORK)
      CALL QEXIT('PAMSET')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setdc1 */
      SUBROUTINE SETDC1(IPRINT)
C***********************************************************************
C
C     Set DIRAC common blocks (DCB)
C
C     The general variable syntax is:
C
C     prefix1+prefix2+body+suffix1+suffix2 (irp,ic)
C
C     prefix1:
C       N   - length of vector
C       N2  - length of full matrix
C       NN  - length of row-packed matrix(lower triangle)
C       I   - offset of vector
C       I2  - offset of full matrix
C       IN  - offset of row-packed matrix(lower triangle)
C       J   - pointer of vector
C       J2  - pointer of full matrix
C       JJ  - pointer of row-packed matrix(lower triangle)
C
C     prefix2:
C       B - boson irreps
C       F - fermion irreps
C
C     body:
C       BAS - AO-basis
C       ORB - MO-basis
C       ISH - inactive spinors
C
C     suffix1:
C       T - symmetry-packing
C       X - no symmetry-packing
C       Q - quaternionic symmetry-packed matrix (NZ matrices adjoined)
C
C     indices:
C       irp - irrep (boson or fermion)
C       ic  - component:
C              0 - total
C              1 - large
C              2 - small
C
C     Written by T.Saue - October 1994
C     Last revision: Nov 23 1994 - tsaue
C
CMI/jan.2006 All quantities are constructed from NAOS, NCOS (symmet.h) variables, which
CMI are calculated in READIN
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "symmet.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbfir.h"
#include "inftap.h"
      CHARACTER COMP(2)*1
C
      IF(IPRINT.GE.2) CALL TITLER('Output from SETDC1','*',103)
C
      COMP(1) = 'L'
      COMP(2) = 'S'
C
C     *********************************************
C     ***** COMMON/DCBBAS/ - Basis dimensions *****
C     *********************************************
C
C     NTBAS(ic)      - total number of basis functions for component ic
C     NBBAS(ibrp,ic) - number of SO basis functions
C                      for boson irrep ibrp and component ic
C     NFBAS(ifrp,ic) - number of SO basis functions
C                      for fermion irrep ifrp and component ic
C     MXBBAS         - max.number of SO basis functions in a boson irrep
C
C     Matrix dimensions:
C        N2BBAST,NNBBAST,N2BBASX,N2BBASXQ
C        N2BAS(ifrp),N2BASQ(ifrp),N2BAST,
C     Offsets:
C        IBBAS(ibrp,ic),
C        I2BBASX(ibrp,jbrp,ic,jc)
C        IBAS(ifrp),I2BASX(ifrp,ifrp),I2BAST(ifrp)
C
      NNBBAST = 0
      N2BBAST = 0
      DO IBREP = 0,MAXREP
        NBBAS(IBREP,0) = NAOS(IBREP+1)
        NNBBAST        = NNBBAST + NBBAS(IBREP,0)*(NBBAS(IBREP,0)+1)/2
        N2BBAST        = N2BBAST + NBBAS(IBREP,0)*NBBAS(IBREP,0)
      ENDDO
      CALL IZERO(NTBAS,3)
      CALL IZERO(NFBAS,6)
      MXBBAS   = 0
      IBASI    = 0
      NBRP     = 4/NZ
C
      N2BAST  = 0
      N2BASTQ = 0
C
      IF(NFSYM.EQ.1)     THEN
        IBAS(1)    = IBASI
        DO IC = 1,2
          DO ISYM = 1,NBRP
            IBSYM           = JFSYM(ISYM,1)
            IBREP           = IBSYM - 1
            IBBAS(IBREP,IC) = IBASI
            NBBAS(IBREP,IC) = NCOS(IBSYM,IC)
            IBASI           = NCOS(IBSYM,IC) + IBASI
            NFBAS(1,IC)     = NCOS(IBSYM,IC) + NFBAS(1,IC)
            MXBBAS          = MAX(MXBBAS,NBBAS(IBREP,IC))
          ENDDO
          NTBAS(IC) = NTBAS(IC) + NFBAS(1,IC)
        ENDDO
        NFBAS(1,0) = NFBAS(1,1)  + NFBAS(1,2)
        N2BAS(1)   = NFBAS(1,0)*NFBAS(1,0)
        N2BAST     = N2BAS(1)
        N2BASQ(1)  = N2BAS(1)*NZ
        I2BAST(1)  = 0
        N2BASTQ    = N2BASQ(1)
      ELSEIF(NFSYM.EQ.2) THEN
        DO IFRP = 1,NFSYM
          IBAS(IFRP)    = IBASI
          DO IC = 1,2
            IF(MDIRAC) THEN
              IP = IFRP
            ELSE
C             IP is gerade(1)/ungerade(2) !
              IP    = MOD(IFRP+IC,2) + 1
            ENDIF
            DO ISYM = 1,NBRP
              IBSYM           = JFSYM(ISYM,IP)
              IBREP           = IBSYM - 1
              IBBAS(IBREP,IC) = IBASI
              NBBAS(IBREP,IC) = NCOS(IBSYM,IC)
              IBASI           = NCOS(IBSYM,IC) + IBASI
              NFBAS(IFRP,IC)  = NCOS(IBSYM,IC) + NFBAS(IFRP,IC)
              MXBBAS          = MAX(MXBBAS,NBBAS(IBREP,IC))
            ENDDO
            NTBAS(IC)       = NTBAS(IC) + NFBAS(IFRP,IC)
          ENDDO
          NFBAS(IFRP,0) = NFBAS(IFRP,1)  + NFBAS(IFRP,2)
          N2BAS(IFRP)   = NFBAS(IFRP,0)*NFBAS(IFRP,0)
          N2BAST        = N2BAST  + N2BAS(IFRP)
          N2BASQ(IFRP)  = N2BAS(IFRP)*NZ
          I2BAST(IFRP)  = N2BASTQ
          N2BASTQ       = N2BASTQ + N2BASQ(IFRP)
        ENDDO
      ENDIF
      NTBAS(0)   = NTBAS(1) + NTBAS(2)
      NNBBASX    = NTBAS(0)*(NTBAS(0)+1)/2
      N2BBASX    = NTBAS(0)*NTBAS(0)
      N2BBASXQ   = N2BBASX*NZ
      MXFBAS = MAX(NFBAS(1,0),NFBAS(2,0))
C
      DO JC = 1,2
        DO IC = 1,2
          DO JBRP = 0,NBSYM-1
            DO IBREP = 0,NBSYM-1
              I2BBASX(IBREP,JBRP,IC,JC) =
     &               NTBAS(0)*IBBAS(JBRP,JC)+IBBAS(IBREP,IC)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C
      IOFF = 0
      DO J = 1,NFSYM
        DO I = 1,NFSYM
          I2BASX(I,J) = IOFF + IBAS(I)
        ENDDO
        IOFF = IOFF + NTBAS(0)*NFBAS(J,0)
      ENDDO
C
C     Print section
C     =============
C
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,'(2X,6(1X,A6))')
     &        'Irrep ','NTBAS ','NBASL ','NBASS ','IBASL ','IBASS '
        WRITE(LUPRI,'(3X,41A1)') ('-',I=1,41)
        WRITE(LUPRI,'(3X,I1,2X,A3,5I7)')
     &    ((I+1),REP(I),(NBBAS(I,J),J=0,2),(IBBAS(I,J),J=1,2),
     &     I=0,NBSYM-1)
        WRITE(LUPRI,'(3X,41A1)') ('-',I=1,41)
        WRITE(LUPRI,'(3X,A4,I2,3I7)')
     &       ('Frep',I,(NFBAS(I,J),J=0,2),I = 1,NFSYM)
        WRITE(LUPRI,'(3X,41A1)') ('-',I=1,41)
        WRITE(LUPRI,'(3X,A6,3I7,4X,A3,I7)') 'Total ',(NTBAS(J),J=0,2),
     &        'Max',MXBBAS
        WRITE(LUPRI,'(/3X,A,2I9))') 'NNBBAST, NNBASX:',NNBBAST,NNBBASX
        WRITE(LUPRI,'(/A)')
     &      '* Matrix offsets to symmetry blocks (I2BBASX):'
        WRITE(LUPRI,'(4(3X,2(A3,2X),2(A1,2X),I10))')
     &     ((((REP(I),REP(J),COMP(IC),COMP(JC),
     &        I2BBASX(I,J,IC,JC),IC = 1,2),JC = 1,2)
     &        ,I=0,NBSYM-1),J=0,NBSYM-1)
      ENDIF
C
C     *****************************************
C     ***** COMMON/PAMIOU/ - Unit numbers *****
C     *****************************************
C
      LUCOEF = 1
      LUOVLP = 2
      LUTMAT = 3
      LU1INT = 19
      LUPMAT = 20
      LUINTR = 0
      LUKRMC = 50 
C
C     *****************************************
C     ***** COMMON/FIRSTN/ -              *****
C     *****************************************
C
C     
      FIRST1 = .TRUE.
      FIRST2 = .TRUE.
      FIRST3 = .TRUE.
      FIRST4 = .TRUE.
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setdc2 */
      SUBROUTINE SETDC2(IPRINT)
C***********************************************************************
C
C     Set DIRAC common blocks (DCB) - second run
C
C     The general variable syntax is:
C
C     prefix1+prefix2+body+suffix1+suffix2 (irp,ic)
C
C     prefix1:
C       N   - length of vector
C       N2  - length of full matrix
C       NN  - length of row-packed matrix(lower triangle)
C       I   - offset of vector
C       I2  - offset of full matrix
C       IN  - offset of row-packed matrix(lower triangle)
C       J   - pointer of vector
C       J2  - pointer of full matrix
C       JJ  - pointer of row-packed matrix(lower triangle)
C
C     prefix2:
C       B - boson irreps
C       F - fermion irreps
C
C     body:
C       BAS - AO-basis
C       ORB - MO-basis
C       ISH - inactive spinors
C
C     suffix1:
C       T - symmetry-packing
C       X - no symmetry-packing
C       Q - quaternionic symmetry-packed matrix (NZ matrices adjoined)
C
C     indices:
C       irp - irrep (boson or fermion)
C       ic  - component:
C              0 - total
C              1 - large
C              2 - small
C
C     Written by T.Saue - November 1994
C     Last revision: Nov 23 1994 - tsaue
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "pgroup.h"
#include "dgroup.h"
      CHARACTER COMP(2)*1
C
C
      COMP(1) = 'L'
      COMP(2) = 'S'
C
C     *********************************************
C     ***** COMMON/DCBBAS/ - Basis dimensions *****
C     *********************************************
C
C     The following dimensions are set in LOWDIN:
C       NFORB(ifrp,ic) - number of MO basis functions
C                        for fermion irrep ifrp and component ic
C     Here we set:
C       NTORB(ic)      - number of MO basis functions for a component ic
C
      IF(.NOT.URKBAL) NTORB(2) = NTORB(1)
C
C     **********************************************
C     ***** COMMON/DCBORB - Orbital dimensions *****
C     **********************************************
C
C     A shell is occupied by a Kramers pair
C
C     NORB(ifrp) - total number of shells
C     NISH(ifrp) - inactive (electron) shells
C     NASH(ifrp) - active (electron) shells
C     NOCC(ifrp) - total number of partially/fully occupied shells
C     NSSH(ifrp) - secondary (electron) shells
C     NESH(ifrp) - electron shells
C     NPSH(ifrp) - positron shells
C
C     NCMO(ifrp) - dimension of MO - coefficient matrix
C     NCAO(ifrp) - dimension of occupied part of MO-coef. matrix
C
      NORBT  = 0
      N2ORBT = 0
      NFROT  = 0
      NISHT  = 0
      N2ISHT = 0
      NASHT  = 0
      N2ASHT = 0
      NOCCT  = 0
      N2OCCT = 0
      NSSHT  = 0
      N2SSHT = 0
      NESHT  = 0
      N2ESHT = 0
      NPSHT  = 0
      N2PSHT = 0
      NCMOT  = 0
      DO 100 IFRP = 1,NFSYM
C
        IORB(IFRP)  = NORBT
        IASH(IFRP)  = NASHT
        IOCC(IFRP)  = NOCCT
        NOCC(IFRP)  = NISH(IFRP)+NASH(IFRP)
C
        N2ORB(IFRP) = NORB(IFRP)*NORB(IFRP)
        NORBT       = NORBT  + NORB(IFRP)
        I2ORBT(IFRP)= N2ORBT*NZ
        N2ORBT      = N2ORBT + N2ORB(IFRP)
        N2ISH(IFRP) = NISH(IFRP)*NISH(IFRP)
        NFROT       = NFROT  + NFRO(IFRP)
        NISHT       = NISHT  + NISH(IFRP)
        N2ISHT      = N2ISHT + N2ISH(IFRP)
C
        N2ASH(IFRP) = NASH(IFRP)*NASH(IFRP)
        NASHT       = NASHT  + NASH(IFRP)
        I2ASHT(IFRP)= N2ASHT*NZ
        N2ASHT      = N2ASHT + N2ASH(IFRP)
C
        N2OCC(IFRP) = NOCC(IFRP)*NOCC(IFRP)
        NOCCT       = NOCCT  + NOCC(IFRP)
        I2OCCT(IFRP)= N2OCCT * NZ
        N2OCCT      = N2OCCT + N2OCC(IFRP)
C
        NSSH(IFRP)  = NESH(IFRP) - NOCC(IFRP)
        N2SSH(IFRP) = NSSH(IFRP)*NSSH(IFRP)
        NSSHT       = NSSHT  + NSSH(IFRP)
        N2SSHT      = N2SSHT + N2SSH(IFRP)
        N2ESH(IFRP) = NESH(IFRP)*NESH(IFRP)
        NESHT       = NESHT  + NESH(IFRP)
        N2ESHT      = N2ESHT + N2ESH(IFRP)
        N2PSH(IFRP) = NPSH(IFRP)*NPSH(IFRP)
        NPSHT       = NPSHT  + NPSH(IFRP)
        N2PSHT      = N2PSHT + N2PSH(IFRP)
        ICMO(IFRP)  = NCMOT
        ICMOQ(IFRP) = NCMOT*NZ
        NCMO(IFRP)  = NFBAS(IFRP,0)*NORB(IFRP)
        NCMOQ(IFRP) = NCMO(IFRP)*NZ
        NCMOT       = NCMOT + NCMO(IFRP)
  100 CONTINUE
      N2ORBTQ = N2ORBT*NZ
      NCMOTQ  = NCMOT*NZ
      IOFF = 0
      DO 200 J = 1,NFSYM
        DO 300 I = 1,NFSYM
          I2ORBX(I,J) = IOFF + IORB(I)
 300    CONTINUE
        IOFF = IOFF + NORB(J)*NORBT
 200  CONTINUE
C
C     Calculate I2ASHX
C
      IOFF = 0
      DO J = 1,NFSYM
         DO I = 1,NFSYM
            I2ASHX(I,J) = IOFF + IASH(I)
         END DO
         IOFF = IOFF + NASH(J)*NASHT
      END DO
C
C     Calculate I2OCCX
C
      IOFF = 0
      DO J = 1,NFSYM
         DO I = 1,NFSYM
            I2OCCX(I,J) = IOFF + IOCC(I)
         END DO
         IOFF = IOFF + NOCC(J)*NOCCT
      END DO
      N2ORBX  = NORBT*NORBT
      N2ORBXQ = NZ*N2ORBX
      N2ISHX  = NISHT*NISHT
      N2ASHX  = NASHT*NASHT
      N2ASHXQ = N2ASHX * NZ
      NNASHX  = NASHT * (NASHT+1) / 2
      N2OCCX  = NOCCT*NOCCT
      N2SSHX  = NSSHT*NSSHT
      N2ESHX  = NESHT*NESHT
      N2PSHX  = NPSHT*NPSHT
C
      IF(IPRINT.GE.1) THEN
        CALL HEADER('Output from SETDC2',-1)
        WRITE(LUPRI,'(A7,9X,A5,14X,A8)') 
     &       'Fermion','Basis','Orbitals'
        WRITE(LUPRI,'(A6,2(5X,3A5))')
     &       'ircop  ','Total',' Lbas',' Sbas',
     &                 'Total',' Eorb',' Porb'
        WRITE(LUPRI,'(2X,A3,1X,5X,3I5,5X,3I5)')
     &    (FREP(I),(NFBAS(I,J),J = 0,2),(NFORB(I,J),J = 0,2),
     &     I = 1,NFSYM)
      ENDIF
C
C     Sanity check if NOCC > NORB
C
      DO I = 1,NFSYM
         IF (NOCC(I) .GT. NORB(I)) THEN
            WRITE(LUPRI,*) 
     &         'Number of occ. orbitals (',NOCC(I),
     &         ') greater than orbitals (',NORB(I),') for symmetry ',I
            CALL QUIT('NOCC>NORB')
         END IF
      END DO
C
C     Set a variable so that we may check whether the orbital information
C     is already defined.
C
      I_DCBORB_SET = 1
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setdhf */
      SUBROUTINE SETDHF(IPRINT)
C***********************************************************************
C
C     Set dimensions for COMMON blocks in SCF-module
C
C     Written by T.Saue November 1994
C     Last revision : Nov 23 1994 - tsaue
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
C
C     *****************************************
C     ***** COMMON/DHFIOU/ - Unit numbers *****
C     *****************************************
C
      LUCYCL    =  9
      LUFCK1    = 10
      LUFCK2    = 11
      LUFOCK    = 12
      LUDIIS    = 13
      LUFCKT    = 14
      LU2HER    = 20
      LUDENS    = 21
      LUEVEC    = 22
      LUCMOS    = 23
      LUSMOS    = 24
C
C     Closed shell Hartree - Fock
C     ===========================
C
      IF (.NOT.AOC) THEN
         NFMAT = 1
      ELSEIF (NOPEN.GT.0) THEN
         NFMAT = NOPEN+1
      ELSE
         IF (NASHT.GT.0) THEN
            NFMAT = 2
         ELSE
            NFMAT = 1
         END IF
      END IF
      DO I = 1,NFMAT
C        Totally symmetric operator
         ISYMOP(I) = 1
C        Fock matrix type
         IFCKOP(I) = 1
C        Hermitian operator
         IHRMOP(I) = 1
      END DO
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SCDENS_2 */
      SUBROUTINE SCDENS_2(SCQ,ISTART,ISTOP,IPRINT)
C***********************************************************************
C
C     INPUT:  Charge of atom
C     OUTPUT: Small component charge
C
C     Small component charges have been obtained from calculations
C     using the relativistic atomic code GRASP.
C     see T.SAUE, Ph.D. thesis, University of Oslo 1997
C
C     Quadratic fit:
C       A[0] = 0.0299562717
C       A[1] = -0.00349803581
C       A[2] = 0.000122746761
C       R square = 0.9985796
C     Cubic fit:
C       A[0] = -0.000873131732
C       A[1] = -2.4427764e-05
C       A[2] = 3.96486383e-05
C       A[3] = 5.32680272e-07
C       R square = 0.9999936
C     4th degree fit:
C       A[0] = 0.00137797077
C       A[1] = -0.000439315949
C       A[2] = 5.73988753e-05
C       A[3] = 2.68047049e-07
C       A[4] = 1.27227511e-09
C       R square = 0.999999
C     5th degree fit:
C       A[0] = 0.000348425837
C       A[1] = -0.000161165557
C       A[2] = 3.90375579e-05
C       A[3] = 7.35391124e-07
C       A[4] = -3.77065743e-09
C       A[5] = 1.93958944e-11
C       R square = 0.9999999
C     Written by T.Saue June 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
      PARAMETER(DTOL = 1.0D-3,D0=0.0D0)
      DIMENSION SCCHRG(103),BFIT(6),SCQ(ISTART:ISTOP)
      SAVE SCCHRG,BFIT
C     
      DATA (SCCHRG(I),I = 1,103)
     +/ 0.000013D0, 0.0001D0, 0.0002D0, 0.0004D0, 0.0007D0, 
     1  0.0010D0, 0.0014D0, 0.0020D0, 0.0026D0, 0.0034D0, 
     +  0.0043D0, 0.0053D0, 0.0065D0, 0.0077D0, 0.0091D0, 
     2  0.0106D0, 0.0123D0, 0.0141D0, 0.0160D0, 0.0181D0, 
     +  0.0203D0, 0.0227D0, 0.0252D0, 0.0279D0, 0.0308D0, 
     3  0.0339D0, 0.0371D0, 0.0405D0, 0.0440D0, 0.0478D0, 
     +  0.0517D0, 0.0558D0, 0.0602D0, 0.0647D0, 0.0694D0, 
     4  0.0743D0, 0.0793D0, 0.0846D0, 0.0901D0, 0.0958D0, 
     +  0.1017D0, 0.1078D0, 0.1141D0, 0.1206D0, 0.1273D0, 
     5  0.1343D0, 0.1415D0, 0.1489D0, 0.1566D0, 0.1644D0, 
     +  0.1725D0, 0.1809D0, 0.1895D0, 0.1983D0, 0.2073D0, 
     6  0.2166D0, 0.2261D0, 0.2359D0, 0.2459D0, 0.2563D0, 
     +  0.2668D0, 0.2777D0, 0.2888D0, 0.3003D0, 0.3118D0, 
     7  0.3238D0, 0.3360D0, 0.3485D0, 0.3614D0, 0.3746D0, 
     +  0.3879D0, 0.4017D0, 0.4157D0, 0.4301D0, 0.4447D0, 
     8  0.4597D0, 0.4750D0, 0.4907D0, 0.5067D0, 0.5230D0, 
     +  0.5396D0, 0.5566D0, 0.5739D0, 0.5916D0, 0.6097D0, 
     9  0.6281D0, 0.6468D0, 0.6659D0, 0.6854D0, 0.7053D0, 
     +  0.7256D0, 0.7464D0, 0.7675D0, 0.7888D0, 0.8107D0, 
     O  0.8332D0, 0.8560D0, 0.8790D0, 0.9027D0, 0.9268D0, 
     +  0.9514D0, 0.9765D0, 1.0020D0 /
      DATA BFIT/0.000348425837D0,
     &         -0.000161165557D0,
     &          3.90375579D-05,
     &          7.35391124D-07,
     &         -3.77065743D-09,
     &          1.93958944D-11/
C
C
C     Check that charge is on list
C            
      DO 10 I = ISTART,ISTOP
        IF(NOORBT(I)) GOTO 10
        IQ   = NINT(CHARGE(I))
        IF(IQ.GT.103) THEN
          SCQ(I) = - POLVAL(5,BFIT,CHARGE(I))
        ELSEIF(IQ.GT.0) THEN
          SCQ(I) = - SCCHRG(IQ)
        ELSE
C       ... this takes care of floating orbitals ! /HJAaJ
          SCQ(I) = D0
        ENDIF
        IF(IPRINT.GE.1) THEN
          IF    (IQ.GT.103) THEN
            WRITE(LUPRI,'(A,F9.4,A,F12.8,A)')
     &            'Charge : ',CHARGE(I),
     &            ' --> Small component charge : ',SCQ(I),
     &            '(extrapolation)'
          ELSE
            WRITE(LUPRI,'(A,F9.4,A,F12.8)')
     &            'Charge : ',CHARGE(I),
     &            ' --> Small component charge : ',SCQ(I)
          ENDIF
        ENDIF            
 10   CONTINUE
      RETURN
C
      ENTRY SCQSET(IZ,SCQVAL)
C
C     entry SCQSET added 30-Mar-2001 HJAaJ
C
      IF ((IZ .GT. 0 .AND. IZ .LE. 103) .AND.
     &    (SCQVAL .GE. D0)) THEN
         WRITE(LUPRI,'(/A,F15.10,A,F15.10/A,I5)')
     &      ' * Tabulated small component charge changed from',
     &      SCCHRG(IZ),'   to',SCQVAL, '   for nuclei no.',IZ
         SCCHRG(IZ) = SCQVAL
      ELSE
         WRITE(LUPRI,'(/A/A,I5,F15.10)')
     &      ' ERROR: Illegal input to .SCQSET :',
     &      ' Z, small component charge:',IZ,SCQVAL
         CALL QUIT('Illegal input to .SCQSET')
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Lvcorr */
      SUBROUTINE LVCORR(DMAT,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Calculate classic repulsive potential of 
C     small component charges
C
C     Modified by O. Fossgaard July 2000 for SCDENS
C
C***********************************************************************

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "frame.h"
C
#include "symmet.h"
#include "nuclei.h"
#include "shells.h"
#include "dcbbas.h"
#include "dcbdhf.h"     
#include "dcbgen.h"
#include "dcbham.h"
      DIMENSION DMAT(N2BBASXQ,*),WORK(*)
#include "memint.h"
      IF(LEVYLE) THEN
        CORRLV = D0
      ELSE
        CALL MEMGET('REAL',KSCQ,NUCIND,WORK,KFREE,LFREE)
C
C     Get small component charges
C     Two options (decided in SCDENS) :
C     - Old method uses a table lookup with to get the DHF limit
C       small component charges
C     - New method uses a Mulliken analysis of the density
C       THIS WILL NOT WORK FOR PARTIAL DENSITIES AS IN MOLTRA !!!!
C
        IF(AOC) THEN
          CALL MEMGET('REAL',KBUF,N2BBASX,WORK,KFREE,LFREE)
          CALL DCOPY(N2BBASX,DMAT,1,WORK(KBUF),1)           
          DO ISHELL = 1, NOPEN
            CALL DAXPY(N2BBASX,DF(ISHELL),DMAT(1,(ISHELL+1)),1,
     &                 WORK(KBUF),1)
          ENDDO
          CALL SCDENS(WORK(KSCQ),WORK(KBUF),WORK(KFREE),LFREE,IPRINT)
        ELSE
          CALL SCDENS(WORK(KSCQ),DMAT,WORK(KFREE),LFREE,IPRINT)
        ENDIF
C
C       Calculate electrostatic correction to nuclear repulsion energy
C     
        CALL LVCORR_1(WORK(KSCQ),IPRINT)
        CALL MEMREL('LVCORR',WORK,KWORK,KWORK,KFREE,LFREE)  
      ENDIF
      RETURN
      
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck scdens */
      SUBROUTINE SCDENS(SCQ,DMAT,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Using density matrix, compute small component charges
C     based on a Mulliken population analysis
C
C     Written by O. Fossgaard July 2000     
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
      LOGICAL TABSCQ
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "symmet.h"
#include "nuclei.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbgen.h"
#include "dcbdhf.h"
      DIMENSION SCQ(*),DMAT(*),WORK(LWORK)      
C
#include "memint.h"
      IF(IPRINT.GE.1) CALL HEADER('Output from SCDENS',-1)
      CALL DZERO(SCQ,NUCIND)
c     ... when to use tabulated charges :
          TABSCQ = (.NOT.LVNEW.AND..NOT.ONECAP).OR.
     &             (ONECAP.AND.(INTV1C .EQ. 2 .OR. INTV1C .EQ. 4)
     &             .AND. ICTLV1C(2) .EQ. 1)
      IF(.NOT.LEVYLE) THEN
       IF(TABSCQ) THEN
C       ... use tabulated small component charges
         IF(IPRINT.GE.1) WRITE(LUPRI,'(A/)')
     &      ' * Using tabulated small component atomic charges:'
         CALL SCDENS_2(SCQ,1,NUCIND,IPRINT)
       ELSE
C       ... use Mulliken small component charges
         CALL MEMGET('REAL',KSMAT,NTBAS(0)*NTBAS(0),WORK,KFREE,LFREE)
         CALL GTOVLX(WORK(KSMAT),SSMTRC)
         CALL MEMGET('INTE',KICLAB,NTBAS(2),WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KGROSS,NTBAS(2),WORK,KFREE,LFREE)
         CALL SCDENS_1(SCQ,DMAT,WORK(KSMAT),WORK(KGROSS),
     &               WORK(KICLAB),WORK(KFREE),LFREE)
C        Adjust for nuclear multiplicities
         DO I = 1,NUCIND
           SCQ(I) = D2*SCQ(I)/FMULT(ISTBNU(I))
         ENDDO
         CALL MEMREL('SCDENS',WORK,KWORK,KWORK,KFREE,LFREE)  
         IF(IPRINT.GE.1) THEN
           WRITE(LUPRI,'(A/)')
     &       ' * Using Mulliken small component atomic charges:'
           WRITE(LUPRI,'(4X,A4,F12.8)') (NAMN(I),SCQ(I),I = 1,NUCIND)
         ENDIF
       ENDIF
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck lvcorr_1 */
      SUBROUTINE LVCORR_1(SCQ,IPRINT)
C***********************************************************************
C
C     Calculate classic repulsive potential of 
C     small component charges (given in SCQ on input)
C
C     Written by T.Saue June 1997
C     Corrected for ONECAP March 2001 HJAaJ
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
C
#include "nuclei.h"
#include "frame.h"
#include "dcbham.h"
      DIMENSION SCQ(*)
C
      IPRINT_local = IPRINT - 3
      IF (ONECAP .AND. INTV1C .EQ. 2) THEN
         DO I = 1,NUCIND
            SCQ(I) = SCQ(I) + CHARGE(I)
         END DO
         CORRLV = GETPOT(SCQ,IPRINT_local) - POTNUC
         IF (IPRINT .GT. 0) WRITE(LUPRI,'(/A,F18.12)')
     &     ' * ONECAP model 1 classical repulsion from small component',
     &     CORRLV
      ELSE
C
C        Calculate repulsion energy for small component atomic charges
C
         CORRLV = GETPOT(SCQ,IPRINT_local)
C
C        Print section
C
         IF (IPRINT .GT. 0) WRITE(LUPRI,'(/A,F18.12)')
     &     ' * Classical repulsion of small component atomic charges:',
     &     CORRLV
      END IF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SCDENS_1 */
      SUBROUTINE SCDENS_1(SCQ,DMAT,SMAT,GROSS,ICLAB,WORK,LWORK)
C***********************************************************************
C
C Calculate classic repulsive potential of 
C     small component charges
C
C     Written by O. Fossgaard July 2000     
C
C***********************************************************************

#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "frame.h"
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
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION SMAT(NBAST,*),DMAT(NBAST,*)
#else
      DIMENSION SMAT(NTBAS(0),NTBAS(0)),DMAT(NTBAS(0),NTBAS(0))
#endif
      DIMENSION ICLAB(*),GROSS(*),SCQ(*)
#include "dcbibt.h"
C
C     Get array ICLAB which maps small component basis functions onto
C     atomic centers.
C
      I = 0
      DO IFRP = 1,NFSYM
         IOFF = IBAS(IFRP) + NFBAS(IFRP,1)
         DO ISBAS  = IOFF+1,IOFF+NFBAS(IFRP,2)
            I = I + 1
            IOLAV  = IATTR(IPLAB(ISBAS,2),2)
            ICENT  = JGET(IOLAV)
            ICLAB(I) = ICENT
         ENDDO   
      ENDDO   
C
C     Calculate the gross population of each small component basis
C     function
C
      IC = 2      
      DO IBRP = 0,MAXREP
         IFRP = JBTOF(IBRP,IC)
         DO I = 1,NBBAS(IBRP,IC)
            IOFF = IBBAS(IBRP,IC) + I
            JOFF = IOFF
            DO J = 1,IFRP
               JOFF = JOFF - NFBAS(J,1)
            ENDDO
            GROSS(JOFF) = DDOT(NBBAS(IBRP,IC),
     &         SMAT(IBBAS(IBRP,IC) + 1,IOFF),1 ,
     &         DMAT(IBBAS(IBRP,IC) + 1,IOFF),1)              
         ENDDO
      ENDDO  
C
C     Calculate the small component charge density on each atom
C
      DO I = 1, NTBAS(2)
         SCQ(ICLAB(I)) = SCQ(ICLAB(I)) - GROSS(I)
      ENDDO
C
      RETURN
C      
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck lcdens */
      SUBROUTINE LCDENS(QLC,DMAT,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Using density matrix, compute large component charges
C     based on a Mulliken population analysis
C
C     Written by Jesper K. Pedersen Aug. 2001
C     (based on SCDENS by O. Fossgaard)
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "symmet.h"
#include "nuclei.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbdhf.h"
      DIMENSION QLC(*),DMAT(*),WORK(LWORK)      
C
#include "memint.h"
      CALL DZERO(QLC,NUCIND)
      CALL MEMGET('REAL',KSMAT,NTBAS(0)*NTBAS(0),WORK,KFREE,LFREE)
      CALL GTOVLX(WORK(KSMAT),SSMTRC)
      CALL MEMGET('INTE',KICLAB,NTBAS(1),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGROSS,NTBAS(1),WORK,KFREE,LFREE)
      CALL LCDEN1(QLC,DMAT,WORK(KSMAT),WORK(KGROSS),
     &               WORK(KICLAB),WORK(KFREE),LFREE)
C        Adjust for nuclear multiplicities
      DO I = 1,NUCIND
        QLC(I) = D2*QLC(I)/FMULT(ISTBNU(I))
      ENDDO
      CALL MEMREL('LCDENS',WORK,KWORK,KWORK,KFREE,LFREE)  
      IF(IPRINT.GE.1) THEN
        CALL HEADER('Output from LCDENS',-1)
        WRITE(LUPRI,'(/1X,A)')
     &       '* Large component Mulliken atomic charges:'
      WRITE(LUPRI,'(4X,A4,F12.8)') (NAMN(I),QLC(I),I = 1,NUCIND)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck lcden1 */
      SUBROUTINE LCDEN1(QLC,DMAT,SMAT,GROSS,ICLAB,WORK,LWORK)
C***********************************************************************
C
C     Using density matrix, compute large component charges
C     based on a Mulliken population analysis
C
C     Written by Jesper K. Pedersen Aug. 2001
C     (based on SCDENS by O. Fossgaard)
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "frame.h"
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
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION SMAT(NBAST,*),DMAT(NBAST,*)
#else
      DIMENSION SMAT(NTBAS(0),NTBAS(0)),DMAT(NTBAS(0),NTBAS(0))
#endif
      DIMENSION ICLAB(*),GROSS(*),QLC(*)
#include "dcbibt.h"

C
C     Get array ICLAB which maps large component basis functions onto
C     atomic centers.
C
      I = 0
      DO IFRP = 1,NFSYM
         IOFF = IBAS(IFRP)
         DO ILBAS  = IOFF+1,IOFF+NFBAS(IFRP,1)
            I = I + 1
            IOLAV  = IATTR(IPLAB(ILBAS,2),2)
            ICENT  = JGET(IOLAV)
            ICLAB(I) = ICENT
         ENDDO
      ENDDO
C
C     Calculate the gross population of each large component basis
C     function
C
      IC = 1      
      DO IBRP = 0,MAXREP
         IFRP = JBTOF(IBRP,IC)
         DO I = 1,NBBAS(IBRP,IC)
            IOFF = IBBAS(IBRP,IC) + I
            JOFF = IOFF
            IF (IFRP.EQ.2) THEN
               JOFF = IOFF - NFBAS(1,2)
            ENDIF
            GROSS(JOFF) = DDOT(NBBAS(IBRP,IC),
     &         SMAT(IBBAS(IBRP,IC) + 1,IOFF),1 ,
     &         DMAT(IBBAS(IBRP,IC) + 1,IOFF),1)              
         ENDDO
      ENDDO  
C
C     Calculate the large component charge density on each atom
C
      DO I = 1, NTBAS(1)
         QLC(ICLAB(I)) = QLC(ICLAB(I)) - GROSS(I)
      ENDDO
C
      RETURN
C      
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck lsqmat */
      SUBROUTINE LSQMAT(QMAT,DMAT,NMPOLOP,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Using density and overlap matrix, compute matrix of large 
C     component charges, placing the charges on the nuclei and
C     at points between the nuclei.
C
C     Written by Jesper K. Pedersen sep. 2001
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "symmet.h"
#include "nuclei.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbgen.h"
      DIMENSION QMAT(NUCDEP,NUCDEP,2,NMPOLOP),DMAT(*),WORK(LWORK)
      DIMENSION INDDIP(3)
      LOGICAL   DIPDONE
      CHARACTER*8  PLABEL(4), RTNLBL(2)
      CHARACTER*17 OPNAME(4)
      DATA PLABEL/'OVERLAP ', 'XDIPLEN ', 'YDIPLEN ', 'ZDIPLEN '/
      DATA OPNAME/'Charges          ', 'Dipole moment (x)', 
     &            'Dipole moment (y)', 'Dipole moment (z)'/
      SAVE DIPDONE
      DATA DIPDONE /.FALSE./
C
#include "memint.h"
      CALL QENTER('LSQMAT')
      CALL TIMER('START ',TIMSTR,TIMEND)
      IF (NMPOLOP .NE. 1 .AND. NMPOLOP .NE. 4) THEN
         CALL QUIT('LSQMAT: NMPOLOP .ne. 1 nor 4')
      END IF
      CALL DZERO(QMAT,NMPOLOP*2*NUCDEP*NUCDEP)
C
C     ... Allocate memory for DMAT in AO basis
C         If necessary, backtransform 
C         real part of density matrix to AO-basis
C
      CALL MEMGET('REAL',KDAO ,NTBAS(0)*NTBAS(0),WORK,KFREE,LFREE)
C
C     ... transform density matrix
C
      CALL DTSOUAO(DMAT,WORK(KDAO),WORK,KFREE,LFREE)
C
C     ... Calculate dipole moment integrals (stored on LU1INT)
C
      IF (NMPOLOP .GT. 1 .AND. .NOT. DIPDONE) THEN
         CALL DEF_DIPOLE(INDDIP,IPRINT)
C
C        ... Open AOPROPER integral file for ONEGEN
C            (is closed again inside ONEGEN)
C
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
         CALL ONEGEN(WORK(KFREE),LFREE)
         DIPDONE = .TRUE.
      END IF
C
C     ... Open AOPROPER integral file (again) for PRPREA
C
      OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
C
C     ... Allocate memory for RLMMAT in AO basis
C         (L=0: Overlap, L=1: Dipole)
C         and for buffer (using NNBBASX for PPRREA,
C         N2BBASX for NBSYM.gt.1)
C
      CALL MEMGET('REAL',KRLMMAT,NTBAS(0)*NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF,N2BBASX,WORK,KFREE,LFREE)
C
C     ... Make list of type (L or S) of each basis function
C
      CALL MEMGET('INTE',KICLAB,NTBAS(0),WORK,KFREE,LFREE)
      CALL LSQMA1(WORK(KICLAB))
C
      DO IOP = 1,NMPOLOP
C
C     ... Get property matrix
C
C        Read integral from file
C        =======================
C
         CALL PRPREA(LU1INT,PLABEL(IOP),RTNLBL,'TFFT',WORK(KBUF),
     &               NNBBASX,IPRINT)
C
C        Make full matrix
C        ================
C
         CALL DSPTSI(NTBAS(0),WORK(KBUF),WORK(KRLMMAT))
C
C     ... If necessary, backtransform overlap to AO-basis
C
         IF(NBSYM.GT.1) THEN
           READ(RTNLBL(2)(4:4),'(I1)') IPRPSY
           IF (IPRPSY .GT. 8 .OR. IPRPSY .LT. 0)
     &     CALL QUIT('LSQMAT: Illegal symmetry on label '//RTNLBL(2))
           IPRPSY = IPRPSY - 1
           CALL DCOPY(N2BBASX,WORK(KRLMMAT),1,WORK(KBUF),1)
           CALL MTSOAO(WORK(KBUF),WORK(KRLMMAT),NTBAS(0),IPRPSY,IPRINT)
         ENDIF
C
C     ... Calculate the matrix of charges/dipole components.
C
         CALL LSQMA2(QMAT(1,1,1,IOP),WORK(KDAO),WORK(KRLMMAT),
     &               WORK(KICLAB))
C
C    Output from LSQMAT
C
         IF(IPRINT .GE. 3) CALL HEADER('Output from LSQMAT',-1)
C
C    1) Write total atomic charges
C
         IF(IPRINT .GE. 1 .AND. IOP.EQ.1) THEN
            CALL MEMGET('REAL',KACHRG,2*NUCIND,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KACHRG),2*NUCIND)
            IOFF = 1
            QLT  = 0.0D0
            QST  = 0.0D0
            DO I = 1,2*NUCIND,2
               IF (IOFF .GT. NUCDEP)
     &            CALL QUIT('LSQMAT : IOFF .GT. NUCDEP')
                  DO J = 1,NUCDEP
                     WORK(KACHRG-1+I)=WORK(KACHRG-1+I)+QMAT(J,IOFF,1,1)
                     WORK(KACHRG  +I)=WORK(KACHRG  +I)+QMAT(J,IOFF,2,1)
                  END DO
               IOFF = IOFF + FMULT(ISTBNU((I+1)/2))
               QLT = QLT + WORK(KACHRG-1+I)*FMULT(ISTBNU((I+1)/2))
               QST = QST + WORK(KACHRG  +I)*FMULT(ISTBNU((I+1)/2))
            ENDDO
            WRITE(LUPRI,'(/17X,A/)')
     &      '** Total atomic charges (Mulliken) **'
            WRITE(LUPRI,'(A/)') 
     &      '       Atom(s)               Q_l                     Q_s'
C           ... first line
            WRITE(LUPRI,'(4X,A3,1P,I2,A3,A4,1P,2D25.12)') 
     &      '   ',INT(FMULT(ISTBNU(1))),' x ',NAMN(1),
     &            WORK(KACHRG),WORK(KACHRG+1)
C           ... the rest
            WRITE(LUPRI,'(4X,A3,1P,I2,A3,A4,1P,2D25.12)') 
     &      (' + ',INT(FMULT(ISTBNU((I+1)/2))),' x ',NAMN((I+1)/2),
     &      WORK(KACHRG-1+I),WORK(KACHRG+I),I = 3,2*NUCIND,2)
            WRITE(LUPRI,'(5X,A)')
     &  '--------------------------------------------------------------'
            WRITE(LUPRI,'(5X,A,1P,2D25.12)') 
     &      '=          ',QLT,QST
         END IF
C
C     2) Write the full matrix of charges or dipole moment integrals.
C
         IF(IPRINT.GE.4) THEN
           DO IC = 1,2
              IF ( IC .EQ. 1) THEN
                 WRITE(LUPRI,'(/A,A/)') 
     &           '  ***   Matrix of large component ',
     &                    OPNAME(IOP)
              ELSE
                 WRITE(LUPRI,'(/A,A/)') 
     &           '  ***   Matrix of small component ',
     &                    OPNAME(IOP)
              END IF
              CALL OUTPUT(QMAT(1,1,IC,IOP),1,NUCDEP,1,NUCDEP,
     &                    NUCDEP,NUCDEP,-1,LUPRI)
           END DO
         ENDIF
      END DO
      CALL MEMREL('LCQMAT',WORK,KWORK,KWORK,KFREE,LFREE)  
      CLOSE(LU1INT,STATUS = 'KEEP')
C
      CALL TIMER('LSQMAT',TIMSTR,TIMEND)
      CALL QEXIT('LSQMAT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck lsqma1 */
      SUBROUTINE LSQMA1(ICLAB)
C***********************************************************************
C
C     Make array ICLAB(I) which maps each basis function (L or S) to
C     a specific nucleus.
C
C     Written by Jesper K. Pedersen sep. 2001
C     
C***********************************************************************

#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxorb.h"
C
#include "dcblab.h"
#include "dgroup.h"
#include "nuclei.h"
#include "dcbbas.h"
      DIMENSION ICLAB(*)
      integer, allocatable :: IDEGN(:,:)
#include "dcbibt.h"
      allocate(IDEGN(NUCIND,NBSYM))
      DO ICENT = 1,NUCIND
        IDEG = 0
        DO ISYM = 1,NBSYM
          IF (NUCNUM(ICENT,ISYM) .NE. 0) THEN
            IDEG = IDEG + 1
            IDEGN(ICENT,IDEG) = NUCNUM(ICENT,ISYM)
          ENDIF
        ENDDO
      ENDDO
      DO 100 I = 1, NTBAS(0)
        ITEMP = IATTR(IPLAB(I,1),1)
        ICENT = JGET(ITEMP)
        IDEG  = KGET(ITEMP)
        ICLAB(I) = IDEGN(ICENT,IDEG)
  100 CONTINUE
      deallocate(IDEGN)
C
      RETURN
C      
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck lsqma2 */
      SUBROUTINE LSQMA2(QMAT,DMAT,RLMMAT,ICLAB)
C***********************************************************************
C
C     Calculate
C
C     Written by Jesper K. Pedersen sep. 2001
C     
C***********************************************************************

#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dcbbas.h"
      PARAMETER(D2=2.0D0)
#include "nuclei.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION RLMMAT(NBAST,*),DMAT(NBAST,*)
#else
      DIMENSION RLMMAT(NTBAS(0),NTBAS(0)),DMAT(NTBAS(0),NTBAS(0))
#endif
      DIMENSION ICLAB(*),QMAT(NUCDEP,NUCDEP,2)
#include "dcbibt.h"
C
C     Calculate the population of each large and small component basis
C     function and add to appropiate element of QMAT(NUCDEP,NUCDEP,IC).
C
      JTBEND = 0
      DO IC = 1,2
C
         JTBSTA = JTBEND + 1
         JTBEND = (JTBSTA - 1) + NTBAS(IC)
         DO I = JTBSTA, JTBEND
            ICENTB = ICLAB(I)
            DO J = JTBSTA, JTBEND
               ICENTA = ICLAB(J)
               QMAT(ICENTA,ICENTB,IC) = QMAT(ICENTA,ICENTB,IC) -
     &                                  D2*RLMMAT(J,I)*DMAT(J,I)
            ENDDO
         ENDDO  
C
      END DO
C
      RETURN
C      
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck barvac */
      SUBROUTINE BARVAC(WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate and write bare vacuum density to file
C
C     Written by T. Saue Nov 7 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
      LOGICAL TOBE
      DIMENSION WORK(LWORK)
#include "memint.h"
      INQUIRE(FILE='DFBVAC',EXIST=TOBE)
      IF(TOBE) THEN
        WRITE(LUPRI,'(A)')
     &   '* WARNING: BARVAC found file DFBVAC. '//
     &   'No reference vacuum generated.'
      ELSE
        CALL MEMGET('REAL',KFAO ,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFMO ,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KTMAT,N2TMT   ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
        CALL BARVA1(WORK(KFAO),WORK(KFMO),WORK(KTMAT),
     &            WORK(KEIG),WORK(KIBRP),IPRINT,WORK,KFREE,LFREE)
        CALL MEMREL('BARVAC',WORK,KWORK,KWORK,KFREE,LFREE)
      ENDIF
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck barva1 */
      SUBROUTINE BARVA1(FAO,FMO,TMAT,EIG,IBEIG,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate and write bare vacuum density to file
C
C     Written by T. Saue Nov 7 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION FAO(*),FMO(*),TMAT(*),EIG(*),IBEIG(*),WORK(*)
C
      LUBUF  = 10
C
C     Generate free particle matrix
C
      CALL FREEMT(FAO,IPRINT,WORK(KFREE),LFREE)      
C
C     Transform to MO-basis FAO --> FMO
C
      CALL OPNFIL(LUBUF,'AOMOMAT','OLD','BARVAC')
      CALL READT(LUBUF,N2TMT,TMAT)
      CLOSE(LUBUF,STATUS='KEEP')
      CALL MKMOFK2(FMO,FAO,TMAT,WORK(KFREE),LFREE)
C
C     Diagonalize
C
      CALL DFDIAG (FMO,EIG,IBEIG,FAO,.FALSE.,WORK,KFREE,LFREE)
C
C     Inspect eigenvalues
C
      DO I = 1,NFSYM
        WRITE(LUPRI,*) '* BARVAC fermion ircop :',I
        DO J = 1,NORB(I)
          WRITE(6,'(I5,E20.10)') J,EIG(IORB(I)+J)
        ENDDO
      ENDDO
C     Backtransform solutions to FMO
C
      DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        CALL BCKTRA(FMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &              FAO(I2ORBT(I)+1),NORB(I),NORB(I),
     &              NORB(I),NZ,NORB(I),1,NFBAS(I,0),
     &              TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZT,
     &              IPRINT)
 10   CONTINUE
C
C     Construct bare vacuum density matrix and write to file
C
      CALL GENDEN(FAO,FMO,8,IPRINT)
      CALL OPNFIL(LUBUF,'DFBVAC','UNKNOWN','BARVAC')
      CALL WRITT(LUBUF,N2BBASXQ,FAO)
      CLOSE(LUBUF,STATUS='KEEP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gofast */
      SUBROUTINE GOTEST(GOFAST,IPRINT)
C***********************************************************************
C
C     Look for MO-transformation matrices and other necessary
C     information in order to skip GMOTRA which may be lengthy
C
C     Written by T. Saue Jun 15 2007
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbbas.h"
C
      LOGICAL GOFAST
      GOFAST = .FALSE.
      INQUIRE(FILE='AOMOMAT',EXIST=GOFAST)
      IF(GOFAST) THEN
        WRITE(LUPRI,'(2A)')
     &  '* WARNING: DIRAC found file AOMOMAT (transformation matrices)',
     &  '  ...will attempt restart.'
      ENDIF
C
C     The necessary files have been found;
C     try to get the necessary information from
C     file AOMOMAT
C
      IF(GOFAST) THEN
        NBRP = 4/NZ
        CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','GOFAST')
C.......skip AO2MO transformation matrix
        READ(LUTMAT)
C.......skip inverse transformation matrix
        READ(LUTMAT)
C.......read records containing orbital information
        READ(LUTMAT,END=10,ERR=10) NZT,URKBAL,SUB_BL,I2COFK
        N2TMT   = 0
        N2TMOTQ = 0
        DO IFRP = 1,NFSYM
          READ(LUTMAT,END=10,ERR=10) 
     &       NORB(IFRP),(NFORB(IFRP,J),J=0,2),
     &       ((NBORB(JSYM,IFRP,IC),JSYM=1,NBRP),IC=1,2),
     &        NISH(IFRP),NESH(IFRP),NPSH(IFRP),
     &        NASH(IFRP),(NACSH(IFRP,IOPEN),IOPEN=1,NOPEN),
     &        NFRO(IFRP),NTMO(IFRP),NISHMF(IFRP),NOCCMF(IFRP),
     &        (NACSHMF(IOPEN,IFRP),IOPEN=1,NOPEN),
     &        NESHMF(IFRP),NPSHMF(IFRP)
          I2TMT(IFRP) = N2TMT
          N2TMT = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZT
          I2TMOT(IFRP) = N2TMOTQ
          N2TMO(IFRP) = NTMO(IFRP)*NTMO(IFRP)
          N2TMOTQ = N2TMOTQ + N2TMO(IFRP)*NZ
        ENDDO
        IF(SUB_BL) THEN
          DO IFRP = 1,NFSYM
            READ(LUTMAT,END=10,ERR=10) N_SUB_BL(IFRP),
     &       ((NORB_SUB(IB,IFRP,IC),IC=0,2),
     &       ID_SUB_BL(IB,IFRP),IB=1,N_SUB_BL(IFRP))
          ENDDO
        ENDIF
        CLOSE(LUTMAT,STATUS='KEEP')
        GOTO 20
C.......did not find the necessary records. 
C.......Look for a formatted file with same information
 10     CONTINUE
        CLOSE(LUTMAT,STATUS='KEEP')
        INQUIRE(FILE='GOFAST',EXIST=GOFAST)
        IF(GOFAST) THEN
          OPEN(LUTMAT,FILE ='GOFAST',STATUS='OLD',
     &             ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
          READ(LUTMAT,*,END=30,ERR=30) NZT,URKBAL,SUB_BL,I2COFK
          N2TMT   = 0
          N2TMOTQ = 0
          DO IFRP = 1,NFSYM
            READ(LUTMAT,*,END=30,ERR=30) 
     &      NORB(IFRP),(NFORB(IFRP,J),J=0,2),
     &      ((NBORB(JSYM,IFRP,IC),JSYM=1,NBRP),IC=1,2),
     &        NISH(IFRP),NESH(IFRP),NPSH(IFRP),
     &        NASH(IFRP),(NACSH(IFRP,IOPEN),IOPEN=1,NOPEN),
     &        NFRO(IFRP),NTMO(IFRP),NISHMF(IFRP),NOCCMF(IFRP),
     &        (NACSHMF(IOPEN,IFRP),IOPEN=1,NOPEN),
     &        NESHMF(IFRP),NPSHMF(IFRP)
            I2TMT(IFRP) = N2TMT
            N2TMT = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZT
            I2TMOT(IFRP) = N2TMOTQ
            N2TMO(IFRP) = NTMO(IFRP)*NTMO(IFRP)
            N2TMOTQ = N2TMOTQ + N2TMO(IFRP)*NZ
          ENDDO
          IF(SUB_BL) THEN
            DO IFRP = 1,NFSYM
              READ(LUTMAT,*,END=30,ERR=30) N_SUB_BL(IFRP),
     &        ((NORB_SUB(IB,IFRP,IC),IC=0,2),
     &        ID_SUB_BL(IB,IFRP),IB=1,N_SUB_BL(IFRP))
            ENDDO
          ENDIF
          GOTO 20
 30       CONTINUE
          GOFAST = .FALSE.
        ENDIF
      ENDIF
 20   CONTINUE
      IF(GOFAST) THEN
        WRITE(LUPRI,'(A)') '*** Restart successful !'
        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,*) 'NZT,URKBAL,SUB_BL,I2COFK...',
     &                    NZT,URKBAL,SUB_BL,I2COFK
          WRITE(LUPRI,*) 'NORB...',(NORB(I),I=1,NFSYM)
          WRITE(LUPRI,*) 'NFORB..',
     &       ((NFORB(IFRP,J),J=0,2),IFRP=1,NFSYM)
          WRITE(LUPRI,*) 'NISH...',(NISH(I),I=1,NFSYM)
          WRITE(LUPRI,*) 'NESH...',(NESH(I),I=1,NFSYM)
          WRITE(LUPRI,*) 'NPSH...',(NPSH(I),I=1,NFSYM)
          WRITE(LUPRI,*) 'NASH...',(NASH(I),I=1,NFSYM)
          WRITE(LUPRI,*) 'NACSH..',
     &       ((NACSH(IFRP,IOPEN),IOPEN=1,NOPEN),IFRP=1,NFSYM)
          WRITE(LUPRI,*) 'NFRO...',(NFRO(I),I=1,NFSYM)
          WRITE(LUPRI,*) 'NTMO...',(NTMO(I),I=1,NFSYM)
          WRITE(LUPRI,*) 'NISHMF.',(NISHMF(I),I=1,NFSYM)
          WRITE(LUPRI,*) 'NOCCMF.',(NOCCMF(I),I=1,NFSYM)
          WRITE(LUPRI,*) 'NACSHMF',
     &       ((NACSHMF(IOPEN,IFRP),IOPEN=1,NOPEN),IFRP=1,NFSYM)
          IF(SUB_BL) THEN
            WRITE(LUPRI,*) 'SUB_BL information:'
            DO IFRP = 1,NFSYM
              WRITE(LUPRI,*) N_SUB_BL(IFRP),
     &        ((NORB_SUB(IB,IFRP,IC),IC=0,2),
     &        ID_SUB_BL(IB,IFRP),IB=1,N_SUB_BL(IFRP))
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Orbdat*/
      SUBROUTINE ORBDAT
C***********************************************************************
C
C     Dump information from GMOTRA as rec 3 on AOMOMAT
C
C     Written by T. Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcbham.h"
C
      NBRP = 4/NZ
      WRITE(LUTMAT) NZT,URKBAL,SUB_BL,I2COFK
      DO IFRP = 1,NFSYM
        WRITE(LUTMAT) NORB(IFRP),(NFORB(IFRP,J),J=0,2),
     &       ((NBORB(JSYM,IFRP,IC),JSYM=1,NBRP),IC=1,2),
     &        NISH(IFRP),NESH(IFRP),NPSH(IFRP),
     &        NASH(IFRP),(NACSH(IFRP,IOPEN),IOPEN=1,NOPEN),
     &        NFRO(IFRP),NTMO(IFRP),NISHMF(IFRP),NOCCMF(IFRP),
     &        (NACSHMF(IOPEN,IFRP),IOPEN=1,NOPEN),
     &        NESHMF(IFRP),NPSHMF(IFRP)
      ENDDO
      IF(SUB_BL) THEN
        DO IFRP = 1,NFSYM
          NN = N_SUB_BL(IFRP)
          WRITE(LUTMAT) 
     &     NN,((NORB_SUB(IB,IFRP,IC),IC=0,2),
     &     ID_SUB_BL(IB,IFRP),IB=1,NN)
        ENDDO
      ENDIF
C
      RETURN
      END
! -- end of dirac/dirset.F --
