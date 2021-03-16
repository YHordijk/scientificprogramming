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

***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ALLOC_LUCITA
*
* Dimensions and
* Allocation of static memory
*
* =====
* Input
* =====
*
* KFREE : Pointer to first element of free space
* Information in /LUCINP/,/ORBINP/,/CSYM/
*
* ======
* Output
* ======
* KFREE : First array of free space after allocation of
*         static memory
* /GLBBAS/,/CDIM/
*
*
* =======
* Version
* =======
*
* Modified Jan 1997
*           Fall 97 (KPGINT1 added )

*. Input
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cstate.inc"
#include "csm.inc"
#include "crun.inc"
#include "cprnt.inc"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

*./CINTFO/
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*.Output
#include "glbbas.inc"
*.1 : One electron integrals( Complete matrix allocated )
      CALL MEMMAN(KINT1,NTOOB ** 2,'ADDS  ',2,'INT1  ')
*. A copy of the original UNMODIFIED 1-elecs ints
      CALL MEMMAN(KINT1O,NTOOB ** 2,'ADDS  ',2,'INT1O ')
*. Zero to avoid problems with elements that will not
*. be initialized
      ZERO = 0.0D0
      CALL DZERO(WORK(KINT1),NTOOB**2)
      CALL DZERO(WORK(KINT1O),NTOOB**2)
*.1.1 : Inactive fock matrix
      CALL MEMMAN(KFI  ,NTOOB ** 2,'ADDS  ',2,'FI    ')
      CALL MEMMAN(KFIO ,NTOOB ** 2,'ADDS  ',2,'FIO   ')
*.2 : Two electron integrals
      IF(NOINT.EQ.0.AND.INCORE.EQ.1) THEN
         CALL MEMMAN(KINT2,NINT2,'ADDS  ',2,'INT2  ')
      END IF
*. Pointers to symmetry block of integrals
      CALL MEMMAN(KPINT1,NBINT1,'ADDS  ',2,'PINT1 ')
      CALL MEMMAN(KPINT2,NBINT2,'ADDS  ',2,'PINT2 ')
*. Pointers to nonsymmetric one-electron integrals
      DO ISM = 1, NSMOB
*. triangular packed
        CALL MEMMAN(KPGINT1(ISM),NSMOB,'ADDS  ',2,'PGINT1')
*. no packing
        CALL MEMMAN(KPGINT1A(ISM),NSMOB,'ADDS  ',2,'PGIN1A')
      END DO
*. Symmetry of last index as a function of initial index
      CALL MEMMAN(KLSM1,NBINT1,'ADDS  ',2,'LSM1   ')
      CALL MEMMAN(KLSM2,NBINT2,'ADDS  ',2,'LSM2   ')
*.3 One-body density
      CALL MEMMAN(KRHO1,NTOOB ** 2,'ADDS  ',2,'RHO1  ')
*.4  Two-body density matrix
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2
      IF(IDENSI.EQ.2)
     &CALL MEMMAN(KRHO2,LRHO2     ,'ADDS  ',2,'RHO2  ')
      KPNIJ = 1
      KIJKK = 1
      IF(MXP1+MXP2+MXQ. NE. 0 .OR. IPROCC.NE.0 ) THEN
*. Allocate memory for explicit hamiltonian and roots
        NSBDET = MXP1 + MXP2 + MXQ
        NSBDETP = MAX(NSBDET,IPROCC)
        MXP = MXP1 + MXP2
        CALL MEMMAN(KSBEVC,NSBDET*NROOT,'ADDS  ',2,'KSBEVC')
        CALL MEMMAN(KSBEVL,       NROOT,'ADDS  ',2,'KSBEVL')
*. KSBIDT must be able to hold list of subspace dets / print dets
        LSCR = NSBDETP+1
        CALL MEMMAN(KSBIDT,LSCR,'ADDS  ',1,'KSBIDT')
*. ( NSBDET + 1 elements allocated since NSBDET+1 values are
* obtained in order to check for degenerencies )
CADDB
        CALL MEMMAN(KSBCNF,LSCR        ,'ADDS  ',1,'KSBCNF')
CADDE
        CALL MEMMAN(KSBIA ,LSCR        ,'ADDS  ',1,'KSBIA ')
        CALL MEMMAN(KSBIB ,LSCR        ,'ADDS  ',1,'KSBIB ')
*
        LH0 = MXP*(MXP+1)/2 + MXP1*MXQ + MXQ
        CALL MEMMAN(KH0,LH0,'ADDS  ',2,'KH0   ')
*. Scratch space for manipulating with H0
        LH0SCR =
     &  MAX(6*NSBDETP,4*NSBDETP+4*NOCOB,
     &      MXP1*(MXP1+1)/2+MXP1**2)
        CALL MEMMAN(KH0SCR,LH0SCR      ,'ADDS  ',2,'KH0SCR')
      ELSE
        KSBEVC = 0
        KSBEVL = 0
        KSBIDT = 0
        KSBCNF = 0
        KH0    = 0
        KH0SCR = 0
      END IF
*. indeces for pair of orbitals symmetry ordered
*
      CALL MEMMAN(KINH1,NTOOB*NTOOB,'ADDS  ',2,'KINH1  ')
*
*. Length of MO-MO and MO_AO expansion file
      LMOMO = 0
      LMOAO = 0
      DO ISM = 1, NSMOB
        LMOMO = LMOMO + NMOS_ENV(ISM)*NMOS_ENV(ISM)
        LMOAO = LMOAO + NMOS_ENV(ISM)*NAOS_ENV(ISM)
      END DO
      CALL MEMMAN(KMOMO  ,LMOMO ,'ADDS  ',2,'MOMO  ')
      IF(NOMOFL.EQ.0) THEN
        WRITE(6,*)
     &  ' Alloc : space for MO-AO matrices will be allocated '
        CALL MEMMAN(KMOAOIN,LMOAO ,'ADDS  ',2,'MOAOIN')
        CALL MEMMAN(KMOAOUT,LMOAO ,'ADDS  ',2,'MOAOUT')
      ELSE
        KMOAOIN = 0
        KMOAOUT = 0
      END IF
*. (more to come )
      CALL MEMMAN(IDUMMY,IDUMMU,'CHECK',IDUMMY,'Dummy ')
      RETURN
      END
***********************************************************************
      SUBROUTINE LMEMCHK(TEXT)

      CHARACTER*(*) TEXT
#include "priunit.h"
*
* Check memory allocated  with the memory manager
      write (lupri,'(/2A)') 'Memory check - ',TEXT
      CALL MEMMAN(IDUM,IDUM,'CHECK ',IDUM,'IDUM  ')
      write (lupri,'(/3A/)') 'Memory check - ',TEXT,' passed'
*
      END
***********************************************************************
      SUBROUTINE MEMMAN(KBASE,KADD,TASK,IR,IDENT)
*
* Memory manager routine
*
* KBASE : New base address
*         If TASK = INI, KBASE is offset for memory to be controlled
*         by MEMMAN
* KADD  : Dimension of array to be added
*         If TASK = INI, KADD is total length of array
* TASK  : = INI  : Initialize                 Character*6
*         = ADDS : Add static memory
*         = ADDL : Add Local memory
*         = FLUSH : Flush local memory
*         = CHECK : Check memory paddings
*         = FREE  : Return first Free word in KBASE
*         = MARK  : Set a mark at current free adress
*         = FLUSM : Flush local memory to previous mark
*         = SFREEM: Return the free memory currently available in real*8 words
* IR    : 1 => integer , 2 => real,
*         ratio between integer and real is RtoI
* IDENT : identifier of memory slice,Character*6
*
* Local Memory not flushed before allocation of additional static memory
* is tranferred to static memory
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER   KBASE, KADD
      INTEGER*8 KFREES, KFREEL, NWORD, IMARK, IBASE, JBASE,
     &          JBASEN, IBASE_L
!               for addressing of WORK
      CHARACTER*6 TASK,IIDENT,MARKC
      CHARACTER IDENT*(*)
      integer, parameter :: npad = 1
      PARAMETER(MAXLVL = 1024 )
      PARAMETER(MAXMRK = 1024)
      COMMON/CMEMO8/NWORD,KFREES,KFREEL,IBASE(MAXLVL),
     &             IMARK(MAXMRK)
      COMMON/CMEMO4/NS,NL,NM,
     &             LENGTH(MAXLVL),IIDENT(MAXLVL),
     &             MARKL(MAXMRK),MARKS(MAXMRK),MARKC(MAXMRK)
*. Two real*8 words, one added NPAD times before each array, another
*. added NPAD times after each array
      DATA PAD1/0.123456789D0/
      DATA PAD2/0.987654321D0/
*
#include "priunit.h"
#include "mxpdim.inc"
#include "wrkspc.inc"
#ifdef INT_STAR8
      integer, parameter :: rtoi = 1
#else
      integer, parameter :: rtoi = 2
#endif
*
      ISTOP   = 0
      ITSOK   = 1
      IPRNTMP = 0
*
      IF(TASK(1:3).EQ.'INI') THEN
*
**.Initialize
*
        NS     = 0
        NL     = 0
        NSNLI  = 0
        KFREES = KBASE
        KFREEL = KBASE
*
        NM        = 1
        MARKC(NM) = '-INI--'
        IMARK(NM) = KFREEL
        MARKL(NM) = NL
        MARKS(NM) = NS
*
*       Compute pointer to end element
*
        NWORD = KBASE + MXPWRD-1
        IPRNTMP = 0
        ISTOP = 0

!       write(lupri,*)'RtoI',RtoI,'for IR',IR
*
** First free word
*
      ELSE IF (TASK(1:4).EQ.'FREE') THEN
       KBASE = KFREEL
       NSNLI = NS+NL
*
**. Static memory
*
      ELSE IF(TASK(1:4).EQ.'ADDS') THEN
!Miro: KFREEL is too big, causing out-of-bounds !
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREES = KFREEL +(KADD+1)/RtoI + 2*NPAD
        ELSE
          KFREES = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREES-1 .GT. NWORD ) THEN
          WRITE(lupri,*)
          WRITE(lupri,*) ' You can''t always get what you want'
          WRITE(lupri,*) ' No, you can''t always get what you want'
          WRITE(lupri,*) ' But if you try sometime, you may find '
          WRITE(lupri,*) ' you get what you need '
          WRITE(lupri,*) '                       Jagger/Richard '
*
          WRITE(lupri,*) ' MEMMAN : work array too short '
          WRITE(lupri,*) ' current and required length ',NWORD,KFREES-1
*
          WRITE(lupri,*) ' Trying to allocate : identifer,length'
          WRITE(lupri,'(20X,A,I9)')IDENT,KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NS = NS + NL + 1
        NL = 0
        NSNLI = NS + NL
        IF(NS.GT.MAXLVL) THEN
          WRITE(lupri,*) ' Too many levels in MEMMAN '
          WRITE(lupri,*) ' Increase MAXLVL from ', MAXLVL
          call quit('*** error in MEMMAN: too many levels. ***')
        END IF
        IIDENT(NS)     = IDENT
        KFREEL         = KFREES
        IBASE(NS)      = KBASE

!     Miro: indexes "KBASE-1","KFREEL-1" are out-of-bounds , see
!     wrkspc.inc; runtime check flags report problem
        work(KBASE-1)  = PAD1
        work(KFREEL-1) = PAD2
        call setvec(work(KBASE-1), PAD1,1)
        call setvec(work(KFREEL-1),PAD2,1)
!       WRITE(lupri,*) ' global tag set at: ',KBASE-1,KFREEL-1
*
**. Local memory
*
      ELSE IF(TASK(1:4).EQ.'ADDL') THEN
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREEL = KFREEL +(KADD+1)/RtoI + 2*NPAD
        ELSE
          KFREEL = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREEL-1 .GT. NWORD ) THEN
          WRITE(lupri,*)
          WRITE(lupri,*) ' You can''t always get what you want'
          WRITE(lupri,*) ' No, you can''t always get what you want'
          WRITE(lupri,*) ' But if you try sometime, you may find '
          WRITE(lupri,*) ' you get what you need '
          WRITE(lupri,*) '                       Jagger/Richard '

          WRITE(lupri,*) ' MEMMAN : work array too short '
          WRITE(lupri,*) ' current and required length ',NWORD,KFREEL-1
          WRITE(lupri,*) ' Trying to allocate : identifer,offset,length'
          WRITE(lupri,'(24X,A,2I9)')IDENT,KBASE,KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NL    =  NL + 1
        NSNLI = NS+NL
        IF(NS+NL.GT.MAXLVL) THEN
          WRITE(lupri,*) ' Too many levels in MEMMAN '
          WRITE(lupri,*) ' Increase MAXLVL from ', MAXLVL
          call quit('*** error in MEMMAN: too many levels. ***')
        END IF
        IIDENT(NSNLI)  = IDENT
        IBASE(NSNLI)   = KBASE

!       work(KBASE-1)  = PAD1
!       work(KFREEL-1) = PAD2
        call setvec(work(KBASE-1), PAD1,1)
        call setvec(work(KFREEL-1),PAD2,1)
!       WRITE(lupri,*) ' local tag set at: ',KBASE-1,KFREEL-1
*
** Flush local memory
*
      ELSE IF(TASK(1:5).EQ.'FLUSH') THEN
        NSNLI = NS+NL
        KFREEL = KFREES
        NL = 0
      ELSE IF(TASK(1:4).EQ.'MARK') THEN
*. Set a mark at current free address
        NM = NM + 1
        IF(NM.GT.MAXMRK) THEN
          WRITE(lupri,*) ' Too many marks  in MEMMAN '
          WRITE(lupri,*) ' Increase MAXMRK from ', MAXMRK
          call quit('*** error in MEMMAN: too many marks. ***')
        END IF
        MARKC(NM) = IDENT
        IMARK(NM) = KFREEL
        MARKL(NM) = NL
        MARKS(NM) = NS
        NSNLI = NS + NL
      ELSE IF (TASK(1:5).EQ.'FLUSM') THEN
        NSNLI = NS+NL
        KFREEL = IMARK(NM)
        IF(KFREES.GT.IMARK(NM)) KFREES = IMARK(NM)
        NL = MARKL(NM)
        NS = MARKS(NM)
        NM = NM - 1
      ELSE IF( TASK(1:5).EQ.'CHECK') THEN
        NSNLI = NS+ NL
      ELSE IF( TASK(1:6) .EQ. 'SFREEM') THEN
*. compute current free memory
        KBASE = NWORD - KFREEL - 1
        NSNLI = NS+ NL
      ELSE
          WRITE(lupri,*) ' MEMMAN : Unknown task parameter ',TASK
          WRITE(lupri,*) ' Too confused to continue  '
          call quit('*** error in MEMMAN: unknown task. ***')
      END IF
*
**. Check paddings
*
      ICHECK = 1
      ITSOK  = 1
!     IF(ICHECK == 1)THEN
      IF(TASK(1:5) == 'CHECK' .OR. ICHECK == 1)THEN

        DO 100 IL = 1, NSNLI

          JBASE  = IBASE(IL)

          IF(IL.NE.NSNLI) THEN
            JBASEN = IBASE(IL+1)
          ELSE
           JBASEN = KFREEL + 1
          END IF

          L1OK = 1
          L2OK = 1
          IF(WORK(JBASE-1).NE.PAD1.OR. WORK(JBASEN-2).NE.PAD2)THEN
            ITSOK = 0
            WRITE(lupri,*) ' Memory problem for : '
            WRITE(lupri,*) '   Level (IL) ',IL

            IF(WORK(JBASE-1)  .NE.PAD1) L1OK = 0
            IF(WORK(JBASEN-2) .NE.PAD2) L2OK = 0
!           write(lupri,*) 'bla bla: JBASE-1, JBASEN-2 ==> ',
!    &      JBASE-1, JBASEN-2
!           write(lupri,*) 'WORK(JBASE-1), WORK(JBASEN-2)',           
!    &      WORK(JBASE-1), WORK(JBASEN-2)
!           do i = 1, il
!            write(lupri,*) 'place: start and end tags',IBASE(I)-1,
!    &       IBASE(I+1)-2
!           end do

            IF(L1OK.EQ.1.AND.L2OK.EQ.1) THEN
              WRITE(lupri,'(4X,A,I20,4X,A)')
     &        IIDENT(IL),IBASE(IL),'    OKAY     OKAY '
            ELSE IF(L1OK.EQ.1.AND.L2OK.EQ.0) THEN
              WRITE(lupri,'(4X,A,I20,4X,A)')
     &        IIDENT(IL),IBASE(IL),'    OKAY       -  '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
              WRITE(lupri,'(4X,A,I20,4X,A)')
     &        IIDENT(IL),IBASE(IL),'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
              WRITE(lupri,'(4X,A,I20,4X,A)')
     &        IIDENT(IL),IBASE(IL),'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.0) THEN
              WRITE(lupri,'(4X,A,I20,4X,A)')
     &        IIDENT(IL),IBASE(IL),'     -       -    '
            END IF
          END IF
  100   CONTINUE
      END IF
 1001 CONTINUE
*
        IF(ITSOK == 0 .OR. IPRNTMP > 0) THEN
          WRITE(lupri,'(A,A,1X,A)') ' Current task : ', TASK, IDENT
          WRITE(lupri,*) ' NS, NL, NSNLI',NS,NL,NSNLI
          IF(ITSOK.EQ.0)
     &    WRITE(lupri,*) '  Sorry to say it , but memory is CORRUPTED '
          WRITE(lupri,*) '  Memory map : '
          WRITE(lupri,*) 
     &    '  Identifier    Offset     Pad1 okay Pad2 okay '
          WRITE(lupri,*) 
     &    '  ========== ============  ========= ========= '
          DO 200 IL = 1, NSNLI
            JBASE = IBASE(IL)
            IF(IL.NE.NSNLI) THEN
             JBASEN = IBASE(IL+1)
            ELSE
             JBASEN = KFREEL + 1
            END IF
            L1OK = 1
            L2OK = 1
            IF(WORK(JBASE-1).NE.PAD1) L1OK = 0
            IF(WORK(JBASEN-2).NE.PAD2) L2OK = 0
C           find local off-set to WORK for print
            IBASE_L = IBASE(IL) - IMARK(1)
            IF(L1OK.EQ.1.AND.L2OK.EQ.1) THEN
               WRITE(lupri,'(4X,A,I12,4X,A)')
     &         IIDENT(IL),IBASE_L,'    OKAY     OKAY '
            ELSE IF(L1OK.EQ.1.AND.L2OK.EQ.0) THEN
               WRITE(lupri,'(4X,A,I12,4X,A)')
     &         IIDENT(IL),IBASE_L,'    OKAY       -  '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(lupri,'(4X,A,I12,4X,A)')
     &         IIDENT(IL),IBASE_L,'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(lupri,'(4X,A,I12,4X,A)')
     &         IIDENT(IL),IBASE_L,'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.0) THEN
               WRITE(lupri,'(4X,A,I12,4X,A)')
     &         IIDENT(IL),IBASE_L,'     -       -    '
            END IF
  200     CONTINUE
*
* Marks
*
          WRITE(lupri,*)
          WRITE(lupri,*) '======='
          WRITE(lupri,*) ' Marks '
          WRITE(lupri,*) '======='
          WRITE(lupri,*)
*
          WRITE(lupri,*) ' Identifier  Start of free memory '
          WRITE(lupri,*) ' ================================='
          DO JMARK = 1, NM
            WRITE(lupri,'(3X,A6,10X,I18)') MARKC(JMARK),IMARK(JMARK)
          END DO
*
        IF (ITSOK.EQ.0)
     &    call quit('*** error in MEMMAN: memory corrupted. ***')
        END IF
*
      IF (ISTOP.NE.0)
     &   call quit(
     & '*** error in MEMMAN: not enough memory available. ***')

      END
