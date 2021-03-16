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

      SUBROUTINE RMEMCHK_opt(check_label)
#include "ipoist8.inc"
#include "priunit.h"
      character (len=20) :: check_label
      write(lupri,'(a,a20)') ' memcheck for label ==> ',check_label
*
* Check memory allocated  with the memory manager
*
      call memmar(KDUM,IDUM,'CHECK ',IDUM,'IDUM  ')
*
      END
***********************************************************************

      SUBROUTINE RMEMCHK
#include "ipoist8.inc"
*
* Check memory allocated  with the memory manager
*
      call memmar(KDUM,IDUM,'CHECK ',IDUM,'IDUM  ')
*
      END
***********************************************************************

      SUBROUTINE MEMMAR(KBASE,KADD,TASK,IR,IDENT)
C
C     Memory manager routine for LUCIAREL "R"
C
C     KBASE : New base address
C             If TASK = INI, KBASE is offset for memory to be controlled
C             by MEMMAR
C     KADD  : Dimension of array to be added - KADD is default INTEGER
C             if TASK = INI, KADD is total length of array
C     TASK  : 
C             = INI    : Initialize                 Character*6
C             = ADDS   : Add static memory
C             = ADDL   : Add Local memory
C             = FLUSH  : Flush local memory
C             = CHECK  : Check memory paddings
C             = FREE   : Return first Free word in KBASE
C             = MARK   : Set a mark at current free adress
C             = FLUSM  : Flush local memory to previous mark
C             = PRINT  : Print memory allocations
C             = SFREEM : Print memory allocations
C     IR    : 1 => integer , 2 => real, 3 => integer*8
C     IDENT : identifier of memory slice, CHARACTER*6
C
C     NOTE:    Local Memory not flushed before allocation of additional 
C              static memory
C              is tranferred to static memory
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C     memory common block
C
#include "memory_r.inc"
C
C     Two REAL*8 words, one added NPAD times before each array, another
C     added NPAD times after each array
C
      DATA PAD1/0.123456789D0/
      DATA PAD2/0.987654321D0/
      INTEGER*8 KBASE_LOCAL
      REAL*8    XBASE_LOCAL
      logical test1,test2
      integer*8 ibase1,ibase2
      REAL*8 work1,work2
C
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "parluci.h"
C
      KBASE_LOCAL = 0
      XBASE_LOCAL = 0.0D0
      ILPRNT      = 0
      IF(ILPRNT .NE. 0 ) THEN
        WRITE(LUWRT,*) ' KBASE KADD, IR ', KBASE,KADD,IR
        WRITE(LUWRT,'(3A)') ' TASK , IDENT ', TASK,IDENT
        WRITE(LUWRT,*) '  NS NL NM ',NS,NL,NM
      WRITE(LUWRT,*) ' MEMMAR: allocate array of length', KADD
      END IF
C
      ISTOP = 0
      ITSOK  = 1
      IPRNTMP = 0
C
      IF(TASK(1:5).EQ.'PRINT') THEN
        IPRNTMP = 1
        NSNLI = NS + NL
C
      ELSE IF(TASK(1:3).EQ.'INI') THEN
C
C      Initialize
C
        NS = 0
        NL = 0
        NSNLI = 0
        NM = 0
        KFREES = KBASE
        KFREEL = KBASE
C
C       Compute pointer to end element
C
        NWORD = KBASE + MXPWRD-1
        IPRNTMP = 0
        ISTOP = 0
C
C     First free word
C
      ELSE IF (TASK(1:4).EQ.'FREE') THEN
        KBASE = KFREEL
        NSNLI = NS+NL
C
C     Static memory
C
      ELSE IF(TASK(1:4).EQ.'ADDS') THEN
        KBASE = KFREEL+NPAD
#ifdef LUCI_DEBUG
        write(luwrt,*) 'bla adds:  KBASE ==> ',KBASE
        write(luwrt,*) 'bla adds: KFREEL ==> ',KFREEL
        write(luwrt,*) 'bla adds: NPAD   ==> ',NPAD
#endif
        IF( IR .EQ. 1 ) THEN
          KFREES = KFREEL +(KADD+1) + 2*NPAD
        ELSE IF( IR .EQ. 2) THEN
          KFREES = KFREEL + KADD + 2*NPAD
        ELSE
          KFREES = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREES-1 .GT. NWORD ) THEN
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' You can''t always get what you want'
          WRITE(LUWRT,*) ' No, you can''t always get what you want'
          WRITE(LUWRT,*) ' But if you try sometime, you may find '
          WRITE(LUWRT,*) ' you get what you need '
          WRITE(LUWRT,*) '                       Jagger/Richard '
C
          WRITE(LUWRT,*) ' MEMMAR : work array too short '
          WRITE(LUWRT,*) ' current and required length ',NWORD,KFREES-1
          WRITE(LUWRT,*) ' difference in double words: ',KFREES-1-NWORD
C
          WRITE(LUWRT,*) ' Trying to allocate : identifer,offset,length'
          WRITE(LUWRT,'(20X,A,2X,2I15)')IDENT, KBASE, KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NS = NS + NL + 1
        NL = 0
        NSNLI = NS + NL
        IF(NS.GT.MAXLVL) THEN
          WRITE(LUWRT,*) ' Too many levels in MEMMAR '
          write(LUWRT,*) ' NS has value  ',NS
          WRITE(LUWRT,*) ' Increase MAXLVL from ', MAXLVL, 
     &    ' to at least ',MAXLVL+1
          Call quit(' *** error in memmar: too many memory access'//
     & ' pointers. increase MAXLVL in luciarel/memory_r.inc and'//
     & ' re-compile. ***')
        END IF
        IIDENT(NS) = IDENT
        KFREEL = KFREES
        IBASE(NS) = KBASE
        
#ifdef LUCI_DEBUG
        write(luwrt,*) 'bla adds: IBASE(NS) ==> ',IBASE(NS)
        write(luwrt,*) 'bla adds: IIDENT(NS)==> ',IIDENT(NS)
        write(luwrt,*) 'bla adds: KFREEL    ==> ',KFREEL
        write(luwrt,*) 'bla adds: KBASE     ==> ',KBASE
#endif

#ifdef OLD_CODE_CAUSING_PROBLEMS
        DO 10 IPAD = 1, NPAD
          WORK(KBASE-NPAD-1+IPAD) = PAD1
          WORK(KFREEL-NPAD-1+IPAD) =  PAD2
   10   CONTINUE
#endif
        WORK(KBASE-1 ) = PAD1
        WORK(KFREEL-1) = PAD2

        call setvec(work(KBASE-1), PAD1,1)
        call setvec(work(KFREEL-1),PAD2,1)

#ifdef LUCI_DEBUG
        write(luwrt,*) 'bla adds: work(pad1)    ==> ',
     &                  WORK(KBASE-1),KBASE-1
        write(luwrt,*) 'bla adds: work(pad2)    ==> ',
     &                  WORK(KFREEL-1),KFREEL-1
#endif
C
C     Local memory
C
      ELSE IF(TASK(1:4).EQ.'ADDL') THEN
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREEL = KFREEL +(KADD+1) + 2*NPAD
        ELSE IF( IR .EQ. 2) THEN
          KFREEL = KFREEL + KADD + 2*NPAD
        ELSE
          KFREEL = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREEL-1 .GT. NWORD ) THEN
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' You can''t always get what you want'
          WRITE(LUWRT,*) ' No, you can''t always get what you want'
          WRITE(LUWRT,*) ' But if you try sometime, you may find '
          WRITE(LUWRT,*) ' you get what you need '
          WRITE(LUWRT,*) '                       Jagger/Richard '

          WRITE(LUWRT,*) ' MEMMAR : work array too short '
          WRITE(LUWRT,*) ' current and required length ',NWORD,KFREEL-1
          WRITE(LUWRT,*) ' difference in double words: ',KFREEL-1-NWORD
          WRITE(LUWRT,*) ' Trying to allocate : identifer,offset,length'
          WRITE(LUWRT,'(24X,A,2X,2I15)')IDENT,KBASE,KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NL =  NL + 1
        NSNLI = NS+NL
        IF(NS+NL.GT.MAXLVL) THEN
          WRITE(LUWRT,*) ' Too many levels in MEMMAR '
          WRITE(LUWRT,*) ' Increase MAXLVL from ', MAXLVL, 
     &    ' to at least ',MAXLVL+1
          Call quit(' *** error in memmar: too many memory access'//
     & ' pointers. increase MAXLVL in luciarel/memory_r.inc and'//
     & ' re-compile. ***')
        END IF
        IIDENT(NS+NL) = IDENT
        IBASE(NS+NL) = KBASE
#ifdef OLD_CODE_CAUSING_PROBLEMS
        DO 20 IPAD = 1, NPAD
          WORK(KBASE-NPAD-1+IPAD) = PAD1
          WORK(KFREEL-NPAD-1+IPAD) =  PAD2
   20   CONTINUE
#endif
        call setvec(work(KBASE-1), PAD1,1)
        call setvec(work(KFREEL-1),PAD2,1)
C
C     Flush local memory
C
      ELSE IF(TASK(1:5).EQ.'FLUSH') THEN
        NSNLI = NS+NL
        KFREEL = KFREES
        NL = 0
C       Flush output unit
        LU6 = 6
        CALL GFLUSH(LU6)
      ELSE IF(TASK(1:4).EQ.'MARK') THEN
C       Set a mark at current free address
        NM = NM + 1
        IF(NM.GT.MAXMRK) THEN
          WRITE(LUWRT,*) ' Too many marks  in MEMMAR '
          WRITE(LUWRT,*) ' Increase MAXMRK from ', MAXMRK
          Call QUIT( ' *** ERROR in MEMMAR: too many memory marks. ***')
        END IF
        MARKC(NM) = IDENT
        MARK_ACT = IDENT 
        IMARK(NM) = KFREEL
        MARKL(NM) = NL
        MARKS(NM) = NS
        NSNLI = NS + NL
      ELSE IF (TASK(1:5).EQ.'FLUSM') THEN
        I_AM_INTOLERANT = 1
        IF(I_AM_INTOLERANT.EQ.1.AND.IDENT(1:6).NE.MARK_ACT(1:6)) THEN
          WRITE(LUWRT,*) ' Error in Flushing :  MARKS not consistent '
          WRITE(LUWRT,'(A,A,A,A)')
     &    ' Actual MARK and MARK to be flushed ',
     &     IDENT,' versus ',MARK_ACT
          Call QUIT( ' *** ERROR in MEMMAR: mark mismatch 
     &                 while flushing***')
        END IF
C
        NSNLI = NS+NL
        KFREEL = IMARK(NM)
        IF(KFREES.GT.IMARK(NM)) KFREES = IMARK(NM)
        IF(NM.GT.1) MARK_ACT = MARKC(NM-1)
        NL = MARKL(NM)
        NS = MARKS(NM)
        NM = NM - 1
C       Flush output unit
        LU6 = 6
        CALL GFLUSH(LU6)
      ELSE IF( TASK(1:5).EQ.'CHECK') THEN
        NSNLI = NS+ NL
      ELSE IF( TASK(1:6) .EQ. 'SFREEM') THEN
C       compute current free memory
        KBASE_LOCAL = NWORD - KFREEL - 1
        NSNLI = NS+ NL
        XBASE_LOCAL  = REAL(KBASE_LOCAL) 
        LMEMFREE_PTR = INT(XBASE_LOCAL)
csk     WRITE(LUWRT,*) 
csk  &          ' MEMMAR: KBASE_LOCAL, XBASE_LOCAL, LMEMFREE_PTR',
csk  &                    KBASE_LOCAL, XBASE_LOCAL, LMEMFREE_PTR
      ELSE
        WRITE(LUWRT,*) ' MEMMAR : Unknown task parameter ',TASK
        WRITE(LUWRT,*) ' Too confused to continue  '
        Call QUIT( ' *** ERROR in MEMMAR: unknown task parameter.*** ')
      END IF
C
C     Check paddings
C
cjp circumvention of gfotran -O bug
      ICHECK = 1
      ITSOK = 1
      IF(TASK(1:5).EQ.'CHECK'.OR.ICHECK.EQ.1) THEN
        DO 100 IL = 1, NSNLI
          JBASE = IBASE(IL)
          IF(IL.NE.NSNLI) THEN
           JBASEN = IBASE(IL+1)
          ELSE
           JBASEN = KFREEL + 1
          END IF
          L1OK = 1
          L2OK = 1
          ibase1=JBASE-NPAD-1
          ibase2=JBASEN-2*NPAD-1
          DO IPAD = 1, NPAD
            work1=WORK(ibase1+IPAD)
            test1=work1.ne.pad1
            IF(test1) THEN
               ITSOK = 0
               WRITE(LUWRT,*) 'memory problem1: ',IIDENT(IL),IBASE(IL)
            END IF
           END DO
           DO IPAD = 1, NPAD
            work2=WORK(ibase2+IPAD)
            test2=work2.ne.pad2
            IF(test2) THEN
               ITSOK = 0
               WRITE(LUWRT,*) 'memory problem2: ',IIDENT(IL),IBASE(IL)
            END IF
           END DO
  100   CONTINUE
      END IF
 1001 CONTINUE
C
        IF(ITSOK.EQ.0.OR.IPRNTMP.NE.0) THEN
          WRITE(LUWRT,'(A,A)') ' Current task : ', TASK
          WRITE(LUWRT,*) ' NS, NL, NSNLI',NS,NL,NSNLI
          WRITE(LUWRT,'(A,A)') ' Current Identifier : ', IDENT
          IF (ITSOK.EQ.0) THEN
            WRITE(LUWRT,*) '  Sorry to say it , but memory is CORRUPTED'
          END IF
          WRITE(LUWRT,*) '  Memory map : '
          WRITE(LUWRT,*) '  Identifier   Offset    Pad1 okay Pad2 okay '
          WRITE(LUWRT,*) '  ========== ==========  ========= ========= '
          DO 200 IL = 1, NSNLI
            JBASE = IBASE(IL)
            IF(IL.NE.NSNLI) THEN
             JBASEN = IBASE(IL+1)
            ELSE
             JBASEN = KFREEL + 1
            END IF
            L1OK = 1
            L2OK = 1
            DO 40 IPAD = 1, NPAD
              IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1) L1OK = 0
   40       CONTINUE
            DO 50 IPAD = 1, NPAD
              IF(WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) L2OK = 0
   50       CONTINUE
            IF(L1OK.EQ.1.AND.L2OK.EQ.1) THEN
               WRITE(LUWRT,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'    OKAY     OKAY '
            ELSE IF(L1OK.EQ.1.AND.L2OK.EQ.0) THEN
               WRITE(LUWRT,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'    OKAY       -  '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(LUWRT,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(LUWRT,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.0) THEN
               WRITE(LUWRT,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       -    '
            END IF
  200     CONTINUE
C
C         Marks
C
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '======='
          WRITE(LUWRT,*) ' Marks '
          WRITE(LUWRT,*) '======='
          WRITE(LUWRT,*)
C
          WRITE(LUWRT,*) ' Identifier  Start of free memory '
          WRITE(LUWRT,*) ' ================================='
          DO JMARK = 1, NM
            WRITE(LUWRT,'(3X,A6,10X,I20)') MARKC(JMARK),IMARK(JMARK)
          END DO
C
          IF(ITSOK.EQ.0) THEN
            WRITE(LUWRT,*) ' *** error in memmar: array with'//
     &      ' identification tag ',IDENT,' was accessed out-of-bounds'//
     &      ' prior to the memory check. ***'
            call quit(
     &      ' *** quit in MEMMAR: out-of-bounds access detected. *** ')
          ENDIF
        END IF
C
      IF(ISTOP.NE.0) THEN
            WRITE(LUWRT,*) ' *** error in memmar: array with'//
     &      ' identification tag ',IDENT,' was accessed out-of-bounds'//
     &      ' prior to the memory check. ***'
        Call QUIT(
     &      ' *** quit in MEMMAR: out-of-bounds access detected. *** ')
      END IF
      END
***********************************************************************

      SUBROUTINE MEMMAR_I8(KBASE,KADD,TASK,IR,IDENT)
C
C     Memory manager routine for LUCIAREL "R"
C
C     KBASE : New base address
C             If TASK = INI, KBASE is offset for memory to be controlled
C             by MEMMAR
C     KADD  : Dimension of array to be added - KADD is INTEGER*8
C             if TASK = INI, KADD is total length of array
C     TASK  :
C             = INI    : Initialize                 Character*6
C             = ADDS   : Add static memory
C             = ADDL   : Add Local memory
C             = FLUSH  : Flush local memory
C             = CHECK  : Check memory paddings
C             = FREE   : Return first Free word in KBASE
C             = MARK   : Set a mark at current free adress
C             = FLUSM  : Flush local memory to previous mark
C             = PRINT  : Print memory allocations
C             = SFREEM : Print memory allocations
C     IR    : 1 => integer , 2 => real, 3 => integer*8
C     IDENT : identifier of memory slice, CHARACTER*6
C
C     NOTE:    Local Memory not flushed before allocation of additional
C              static memory
C              is tranferred to static memory
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C     memory common block
C
#include "memory_r.inc"
C
C     MEMMAR_I8 specific
C
      INTEGER*8 KADD
C
C     Two real*8 words, one added NPAD times before each array, another
C     added NPAD times after each array
C
      DATA PAD1/0.123456789D0/
      DATA PAD2/0.987654321D0/
C
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "parluci.h"
C
cjp
      logical test1,test2
      integer*8 ibase1,ibase2
      REAL*8 work1,work2

      ILPRNT = 0
      IF(ILPRNT .NE. 0 ) THEN
        WRITE(LUWRT,*) ' KBASE KADD, IR ', KBASE,KADD,IR
        WRITE(LUWRT,'(3A)') ' TASK , IDENT ', TASK,IDENT
        WRITE(LUWRT,*) '  NS NL NM ',NS,NL,NM
      WRITE(LUWRT,*) ' MEMMAR_I8: allocate array of length', KADD
      END IF
C
      ISTOP = 0
      ITSOK  = 1
      IPRNTMP = 0
C
C     Static memory
C
      IF(TASK(1:4).EQ.'ADDS') THEN
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREES = KFREEL +(KADD+1) + 2*NPAD
        ELSE IF( IR .EQ. 2) THEN
          KFREES = KFREEL + KADD + 2*NPAD
        ELSE
          KFREES = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREES-1 .GT. NWORD ) THEN
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' You can''t always get what you want'
          WRITE(LUWRT,*) ' No, you can''t always get what you want'
          WRITE(LUWRT,*) ' But if you try sometime, you may find '
          WRITE(LUWRT,*) ' you get what you need '
          WRITE(LUWRT,*) '                       Jagger/Richard '
C
          WRITE(LUWRT,*) ' MEMMAR : work array too short '
          WRITE(LUWRT,*) ' current and required length ',NWORD,KFREES-1
          WRITE(LUWRT,*) ' difference in double words: ',KFREES-1-NWORD
C
          WRITE(LUWRT,*) ' Trying to allocate : identifer,offset,length'
          WRITE(LUWRT,'(20X,A,2X,2I15)')IDENT, KBASE, KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NS = NS + NL + 1
        NL = 0
        NSNLI = NS + NL
        IF(NS.GT.MAXLVL) THEN
          WRITE(LUWRT,*) ' Too many levels in MEMMAR '
          write(LUWRT,*) ' NS has value  ',NS
          WRITE(LUWRT,*) ' Increase MAXLVL from ', MAXLVL, 
     &    ' to at least ',MAXLVL+1
          Call quit(' *** ERROR in memmar_i8: too many memory access'//
     & ' pointers. Increase MAXLVL in luciarel/memory_r.inc and'//
     & ' re-compile. ***')
        END IF
        IIDENT(NS) = IDENT
        KFREEL = KFREES
        IBASE(NS) = KBASE
#ifdef OLD_CODE_CAUSING_PROBLEMS
        DO 10 IPAD = 1, NPAD
          WORK(KBASE-NPAD-1+IPAD) = PAD1
          WORK(KFREEL-NPAD-1+IPAD) =  PAD2
   10   CONTINUE
#endif
        call setvec(work(KBASE-1), PAD1,1)
        call setvec(work(KFREEL-1),PAD2,1)

C
C     Local memory
C
      ELSE IF(TASK(1:4).EQ.'ADDL') THEN
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREEL = KFREEL +(KADD+1) + 2*NPAD
        ELSE IF( IR .EQ. 2) THEN
          KFREEL = KFREEL + KADD + 2*NPAD
        ELSE
          KFREEL = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREEL-1 .GT. NWORD ) THEN
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' You can''t always get what you want'
          WRITE(LUWRT,*) ' No, you can''t always get what you want'
          WRITE(LUWRT,*) ' But if you try sometime, you may find '
          WRITE(LUWRT,*) ' you get what you need '
          WRITE(LUWRT,*) '                       Jagger/Richard '

          WRITE(LUWRT,*) ' MEMMAR : work array too short '
          WRITE(LUWRT,*) ' current and required length ',NWORD,KFREEL-1
          WRITE(LUWRT,*) ' difference in double words: ',KFREEL-1-NWORD
          WRITE(LUWRT,*) ' Trying to allocate : identifer,offset,length'
          WRITE(LUWRT,'(24X,A,2X,2I15)')IDENT,KBASE,KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NL =  NL + 1
        NSNLI = NS+NL
        IF(NS+NL.GT.MAXLVL) THEN
          WRITE(LUWRT,*) ' Too many levels in MEMMAR '
          WRITE(LUWRT,*) ' Increase MAXLVL from ', MAXLVL, 
     &    ' to at least ',MAXLVL+1
          Call quit(' *** ERROR in memmar_i8: too many memory access'//
     & ' pointers. Increase MAXLVL in luciarel/memory_r.inc and'//
     & ' re-compile. ***')
        END IF
        IIDENT(NS+NL) = IDENT
        IBASE(NS+NL) = KBASE
#ifdef OLD_CODE_CAUSING_PROBLEMS
        DO 20 IPAD = 1, NPAD
          WORK(KBASE-NPAD-1+IPAD) = PAD1
          WORK(KFREEL-NPAD-1+IPAD) =  PAD2
   20   CONTINUE
#endif
        call setvec(work(KBASE-1), PAD1,1)
        call setvec(work(KFREEL-1),PAD2,1)

      ELSE
        WRITE(LUWRT,*) ' MEMMAR_I8 : Unknown task parameter ',TASK
        WRITE(LUWRT,*) ' Too confused to continue  '
        call quit('*** ERROR in MEMMAR_I8: unknown task parameter.***')
      END IF
C
C     Check paddings
C
cjp circumvention of gfortran -O bug
      ICHECK = 1
      ITSOK = 1
      IF(TASK(1:5).EQ.'CHECK'.OR.ICHECK.EQ.1) THEN
        DO 100 IL = 1, NSNLI
          JBASE = IBASE(IL)
          IF(IL.NE.NSNLI) THEN
           JBASEN = IBASE(IL+1)
          ELSE
           JBASEN = KFREEL + 1
          END IF
          L1OK = 1
          L2OK = 1
          ibase1=JBASE-NPAD-1
          ibase2=JBASEN-2*NPAD-1
          DO IPAD = 1, NPAD
            work1=WORK(ibase1+IPAD)
            test1=work1.ne.pad1
            IF(test1) THEN
               ITSOK = 0
               WRITE(LUWRT,*) 'memory problem3: ',IIDENT(IL),IBASE(IL)
            END IF
           END DO
           DO IPAD = 1, NPAD
            work2=WORK(ibase2+IPAD)
            test2=work2.ne.pad2
            IF(test2) THEN
               ITSOK = 0
               WRITE(LUWRT,*) 'memory problem4: ',IIDENT(IL),IBASE(IL)
            END IF
           END DO
  100   CONTINUE
      END IF
 1001 CONTINUE
C
        IF(ITSOK.EQ.0.OR.IPRNTMP.NE.0) THEN
          WRITE(LUWRT,'(A,A)') ' Current task : ', TASK
          WRITE(LUWRT,*) ' NS, NL, NSNLI',NS,NL,NSNLI
          WRITE(LUWRT,'(A,A)') ' Current Identifier : ', IDENT
          IF (ITSOK.EQ.0) THEN
            WRITE(LUWRT,*) '  Sorry to say it , but memory is CORRUPTED'
          END IF
          WRITE(LUWRT,*) '  Memory map : '
          WRITE(LUWRT,*) '  Identifier   Offset    Pad1 okay Pad2 okay '
          WRITE(LUWRT,*) '  ========== ==========  ========= ========= '
          DO 200 IL = 1, NSNLI
            JBASE = IBASE(IL)
            IF(IL.NE.NSNLI) THEN
             JBASEN = IBASE(IL+1)
            ELSE
             JBASEN = KFREEL + 1
            END IF
            L1OK = 1
            L2OK = 1
            DO 40 IPAD = 1, NPAD
              IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1) L1OK = 0
   40       CONTINUE
            DO 50 IPAD = 1, NPAD
              IF(WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) L2OK = 0
   50       CONTINUE
            IF(L1OK.EQ.1.AND.L2OK.EQ.1) THEN
               WRITE(LUWRT,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'    OKAY     OKAY '
            ELSE IF(L1OK.EQ.1.AND.L2OK.EQ.0) THEN
               WRITE(LUWRT,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'    OKAY       -  '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(LUWRT,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(LUWRT,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.0) THEN
               WRITE(LUWRT,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       -    '
            END IF
  200     CONTINUE
C
C         Marks
C
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '======='
          WRITE(LUWRT,*) ' Marks '
          WRITE(LUWRT,*) '======='
          WRITE(LUWRT,*)
C
          WRITE(LUWRT,*) ' Identifier  Start of free memory '
          WRITE(LUWRT,*) ' ================================='
          DO JMARK = 1, NM
            WRITE(LUWRT,'(3X,A6,10X,I20)') MARKC(JMARK),IMARK(JMARK)
          END DO
C
          IF(ITSOK.EQ.0) THEN
            WRITE(LUWRT,*) ' *** error in memmar_i8: array with'//
     &    ' identification tag ',IDENT,' was accessed out-of-bounds'//
     &    ' prior to the memory check. ***'
            Call QUIT(
     &    ' *** quit in MEMMAR_i8: out-of-bounds access detected. *** ')
          END IF
        END IF
C
      IF(ISTOP.NE.0) THEN       
        WRITE(LUWRT,*) ' *** error in memmar_i8: array with'//
     &  ' identification tag ',IDENT,' was accessed out-of-bounds'//
     &  ' prior to the memory check. ***'
        Call QUIT(
     &  ' *** quit in MEMMAR_i8: out-of-bounds access detected. *** ')
      END IF
      END
