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
C
!/* Comdeck memlog */
!060316-hjaaj: for DBG_LENID check memcb1.h
!910110-hjaaj: MEMREL: removed non-standard char. concat. in CALL MEMCHK
!901004-hjaaj: NWNMEM introduced, no final check on NWNMEM (yet)
!900420-hjaaj: MEMPKG.U: corrected errors for IRAT=1 (e.g. Cray)
!891208-hjaaj: MEMINI: new routine to define output units in MEM* routines
!              (makes MEM* routines indep. of sirius, uses LUWMEM not LUWSIR)
!891114-hjaaj: MEMREL: new routine to release memory

C  /* Deck memini */
      SUBROUTINE MEMINI(LUWME1,LUEME1)
C
C Copyright 8-Dec-1989,4-Oct-1990 Hans Joergen Aa. Jensen
C
C     Define output units for memory allocation routines
C     Initialize mempkg warning count NWNMEM to zero
C
#include "implicit.h"
#include "memcb1.h"
C
      LUWMEM = LUWME1
      LUEMEM = LUEME1
      NWNMEM = 0
      IDENT8 = '** KFREE'
      READ (IDENT8,'(A8)') WIDENT_KFREE
      RETURN
      END
C  /* Deck memget */
      SUBROUTINE MEMGET(DATA_TYPE,KBASE,LENGTH,WORK,KFREE,LFREE)
C
C Copyright 9-Jan-1988/27-Jul-1993 Hans Joergen Aa. Jensen
C
C l.r. 980824-hjaaj: QENTER/QEXIT calls use too much time because MEMGET
C is called often. Now only call QENTER/QEXIT if warning or error.
C Restructured code (e.g. changed errors from IF () THEN to IF () GO TO)
C to get fewer logical tests and fewer jumps after logical tests.
C
C     Memory allocation
C
!-----------------------------------------------------------------------
!  N.B. (sya, 2007.05.08)
!     The length of array is counted by REAL*8.
!
!     INTENT(IN)    :: TYPE, LENGTH
!     INTENT(OUT)   :: KBASE
!     INTENT(INOUT) :: KFREE, LFREE
!
!     Before call
!       1                 KFREE[old]                     SIZE(WORK)
!       +-----------------+------------------------------+
!                          <------------LFREE----------->
!                          <--LENGTH/IRAT-->
!
!     After call
!       1                 KBASE             KFREE[new]   SIZE(WORK)
!       +-----------------+-----------------+------------+
!                                            <---LFREE-->
!                         KBASE=KFREE[old]+1
!
!-----------------------------------------------------------------------
#include "implicit.h"
      CHARACTER*(*) DATA_TYPE
      DIMENSION     WORK(*)
      CALL MEMGET2(DATA_TYPE,'-MEMGET-',KBASE,LENGTH,WORK,KFREE,LFREE)
      RETURN
      END
C  /* Deck memget */
      SUBROUTINE MEMGET2(DATA_TYPE,IDENT,KBASE,LENGTH,WORK,KFREE,LFREE)
C
C Copyright 9-Jan-1988/27-Jul-1993 Hans Joergen Aa. Jensen
C
C l.r. 980824-hjaaj: QENTER/QEXIT calls use too much time because MEMGET
C is called often. Now only call QENTER/QEXIT if warning or error.
C Restructured code (e.g. changed errors from IF () THEN to IF () GO TO)
C to get fewer logical tests and fewer jumps after logical tests.
C
C     Memory allocation
C
#include "implicit.h"
      CHARACTER*(*) DATA_TYPE, IDENT
      DIMENSION     WORK(*)
C
#include "iratdef.h"
#include "memcb1.h"
#include "logging.h"
C
Chj1  CALL QENTER('MEMGET2 ')
      IF (DATA_TYPE(1:4) .EQ. 'REAL' .OR.
     &    DATA_TYPE(1:4) .EQ. 'INT8') THEN
         LREAL = LENGTH
      ELSE IF (DATA_TYPE(1:4) .EQ. 'COMP') THEN
         LREAL = 2*LENGTH
      ELSE IF (DATA_TYPE(1:4) .EQ. 'INTE') THEN
         LREAL = (LENGTH-1)/IRAT + 1
      ELSE IF (DATA_TYPE(1:4) .EQ. 'INT4') THEN
         LREAL = (LENGTH-1)/2 + 1
      ELSE IF (DATA_TYPE(1:4) .EQ. 'INT2') THEN
         LREAL = (LENGTH-1)/4 + 1
      ELSE IF (DATA_TYPE(1:4) .EQ. 'LOGI') THEN
         LREAL = (LENGTH-1)/LRAT + 1
      ELSE IF (DATA_TYPE(1:4) .EQ. 'WORK') THEN
         LENGTH = LFREE - 2*LENID
C                         ^- make sure OK for KFREE .eq. 1
         LREAL  = LENGTH
      ELSE
         CALL QENTER('MEMGET2 ')
         WRITE (LUWMEM,'(/2A/2A)')
     *      ' MEMGET2 ERROR, illegal data type : ',DATA_TYPE,
     *      ' MEMGET2 ERROR, for ident         : ',IDENT
         CALL QTRACE(LUWMEM)
         IF (LUEMEM .NE. LUWMEM) THEN
            WRITE (LUEMEM,'(/2A/2A)')
     *      ' MEMGET2 ERROR, illegal data type : ',DATA_TYPE,
     *      ' MEMGET2 ERROR, for ident         : ',IDENT
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMGET2: programming ERROR, illegal data type')
      END IF
      IF (LOGMEMGET) THEN
         WRITE (LUPRI,'(/A,2I10)') '@memget:',LREAL,LFREE
      ENDIF
C
      IDENT8 = IDENT
C
      IF (KFREE .NE. 1) THEN

C        length error check:
         IF (LREAL+LENID .GT. LFREE) GO TO 8001

C        Error check, memory check of WORK(KFREE) OK ? :

         WMEMCK = WORK(KFREE-1)
         IF (MEMCK .NE. MEMID) GO TO 8002
C        ... comparing integers, to be sure no renormalization of float
#if DBG_LENID
         WMEMCK = WORK(KFREE-LENID+2)
         IF (MEMCK .NE. MEMID) GO TO 8002
#endif

C        Warning check, is length of WORK(KFREE) .ne. 0 ? :

         WMEMCK = WORK(KFREE-LENID+1)
         IF (MEMCK .NE. LZERO) GO TO 8003
C        ... go issue warning and continue at 7003
      ELSE
C        KFREE = 1 case, initialize memory vector info:
         IF (LREAL+2*LENID .GT. LFREE) GO TO 8001
#if DBG_LENID
         do i = 3,lenid
           WORK(i) = WMEMID
         end do
#else
         WORK(3) = WMEMID
#endif
         KFREE = KFREE + LENID
         LFREE = LFREE - LENID
      END IF

C     save ident and length for new allocation

 7003 CONTINUE
      READ (IDENT8,'(A8)') WORK(KFREE-LENID)
      WORK(KFREE-LENID+1) = WLREAL

      KBASE  = KFREE
      KFREE  = KBASE + (LREAL + LENID)
      LFREE  = LFREE - (LREAL + LENID)
      WORK(KFREE-LENID)   = WIDENT_KFREE
      WORK(KFREE-LENID+1) = WLZERO
#ifdef DBG_LENID
      do i = 2,lenid-1
         WORK(KFREE-LENID+i) = WMEMID
      end do
#else
      WORK(KFREE-1) = WMEMID
#endif
Chj1  CALL QEXIT('MEMGET2 ') ! too expensive as MEMGET is called many times
      RETURN
C
C     error branches:
C
C        error branch  IF (LREAL+LENIDS .GT. LFREE) :
C
 8001    CALL QENTER('MEMGET2 ')
         IF (KFREE .GT. 1) THEN
            LENIDS = LENID
         ELSE
            LENIDS = 2*LENID
         END IF
         WRITE (LUWMEM,1010) DATA_TYPE,LREAL+LENIDS,LFREE
         CALL MEM_SHORTAGE_INFO(LREAL+LENIDS,LFREE)
         WRITE (LUWMEM,'(/2A/2A/A,I10)')
     &      ' MEMGET2 INFO, type  : ',DATA_TYPE,
     &      ' MEMGET2 INFO, ident : ',IDENT,
     &      ' MEMGET2 INFO, length: ',LENGTH
         IF (KFREE.GT.1)
     &      CALL MEMCHK('MEMGET ERROR (insuff. memory)',WORK,1)
         CALL QTRACE(LUWMEM)
         IF (LUEMEM .NE. LUWMEM) THEN
            WRITE (LUEMEM,1010) DATA_TYPE,LREAL+LENIDS,LFREE
            CALL MEM_SHORTAGE_INFO(LREAL+LENIDS,LFREE)
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMGET ERROR, insufficient work space in memory')
 1010 FORMAT(
     &/' MEMGET2 ERROR, insufficient free space for next allocation ',A,
     &/T17,'( Need:',I10,', available (LFREE):',I10,' )')
C
C        error branch for invalid MEMID for WORK(KFREE)
C
 8002    CALL QENTER('MEMGET2 ')
         WRITE (LUWMEM,'(/A/A,I12,2(/A,I25))')
     &      ' MEMGET2 ERROR',
     &      ' KFREE =',KFREE,
     &      ' found memory checks:',MEMCK,
     &      ' expected           :',MEMID
         WRITE (LUWMEM,'(/2A/2A/A,I10)')
     &      ' MEMGET2 INFO, type  : ',DATA_TYPE,
     &      ' MEMGET2 INFO, ident : ',IDENT,
     &      ' MEMGET2 INFO, length: ',LENGTH
         CALL MEMCHK('MEMGET2 call of MEMCHK:',WORK,1)
         CALL QUIT('MEMGET2 ERROR, not a valid memget '//
     &             'id in work(kfree)')
C
C        warning branch for non-zero length for WORK(KFREE):
C
 8003    NWNMEM = NWNMEM + 1
         CALL QENTER('MEMGET2 ')
         WRITE (LUWMEM,'(/A/A,I10,/A,I25)')
     *      ' MEMGET2 WARNING, nonzero allocation work(kfree)',
     *      ' KFREE =',KFREE,
     *      ' found memory allocation :',MEMCK
         WRITE (LUWMEM,'(/2A/2A/A,I10)')
     &      ' MEMGET2 INFO, type  : ',DATA_TYPE,
     &      ' MEMGET2 INFO, ident : ',IDENT,
     &      ' MEMGET2 INFO, length: ',LENGTH
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMGET2 WARNING, nonzero allocation work(kfree)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QEXIT('MEMGET2 ')
C        Warning issued, continue and do the assignment now
         GO TO 7003
C
C     end of MEMGET2
C
      END

      SUBROUTINE MEM_SHORTAGE_INFO (MEM_NEED, MEM_DISPOSAL)
!------------------------------------------------------------------------
!
!  Print out worth info on how much the memory has to be increased
! in the its shortage
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"

      integer*8 mem_need

      ! calculate the difference in megawords
      DIFF_MEM=(DFLOAT(MEM_NEED)-DFLOAT(MEM_DISPOSAL))/1000000.0D0
      DIFF_MEM=DIFF_MEM+1.0D0

      WRITE(LUPRI,'(/,/,8X,A,/,A,F6.0,A)')
     & "**** To overcome this memory crash: ****",
     & "...please increase the current memory allocation "//
     & "(pam 'mw' or 'nw' parameter) by at least", 
     & DIFF_MEM," megawords."

#ifdef VAR_MPI
      WRITE(LUPRI,'(2X,A,/,A)') "Well, this is parallel run.",
     & "Check in the error output whether"//
     &" the memory increase is for the master"//
     &" (pam 'mw' parameter) or for the slave (pam 'nw' parameter)."
#else
      WRITE(LUPRI,'(2X,A,A)') "Well, this is serial run. ",
     &"Increase the pam 'mw' parameter."
#endif
      RETURN
      END

C  /* Deck memrel */
      SUBROUTINE MEMREL(TEXT,WORK,KFIRST,KREL,KFREE,LFREE)
C
C 14-Mar-1989 Hans Joergen Aa. Jensen
C
C Check memory allocation from WORK(KFIRST)
C then free memory from WORK(KREL)
C
!-----------------------------------------------------------------------
!  N.B. (sya, 2007.05.08)
!     Usually, KFIRST is set to equal to KFREE.  However, another value
!     is possible.
!
!     INTENT(IN)    :: TEXT, KFIRST, KREL
!     INTENT(OUT)   ::
!     INTENT(INOUT) :: KFREE, LFREE
!
!     Before call
!       1                 KREL              KFREE        SIZE(WORK)
!       +-----------------+-----------------+------------+
!                                            <---LFREE-->
!
!     After call
!       1                 KFREE                          SIZE(WORK)
!       +-----------------+------------------------------+
!                          <------------LFREE----------->
!
!-----------------------------------------------------------------------
#include "implicit.h"
      DIMENSION     WORK(*)
      CHARACTER*(*) TEXT
C
#include "memcb1.h"
C
      CHARACTER     TMPTXT*71
C
      CALL QENTER('MEMREL  ')
C
C     First check if memory allocation is intact
C     (check if anything has been out of bounds)
C
      LENTXT = LEN(TEXT)
      LENTXT = MIN(LENTXT,50)
      TMPTXT = TEXT(1:LENTXT) // ' (called from MEMREL)'
      LENTXT = LENTXT + 21
      CALL MEMCHK(TMPTXT(1:LENTXT),WORK,KFIRST)
C
C     Check if WORK(krel) has legal memid
C
      IF (KREL .GT. KFREE) GO TO 8001
      IF (KREL .NE. 1) THEN
         KFREE1 = KREL
      ELSE
         KFREE1 = 1 + LENID
      END IF

      WMEMCK = WORK(KFREE1-1)
      IF (MEMCK .NE. MEMID) GO TO 8002
C
C     release ...
C
      WORK(KFREE1-LENID  ) = WIDENT_KFREE
      WORK(KFREE1-LENID+1) = WLZERO
      LFREE = LFREE + (KFREE-KREL)
      KFREE = KREL
C
      CALL QEXIT('MEMREL  ')
      RETURN
C
C     error branch IF (KREL .GT. KFREE) THEN
C
 8001    WRITE (LUWMEM,'(/A/2A/A,3I10)')
     *      ' MEMREL ERROR, krel .gt. kfree',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KREL,KFREE =',KFIRST,KREL,KFREE
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMREL ERROR, krel .gt. kfree',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMREL ERROR, krel .gt. kfree')
C
C     error branch IF (MEMCK .NE. MEMID) THEN
C
 8002    WRITE (LUWMEM,'(/A/2A/A,2I10,2(/A,I25))')
     *      ' MEMREL ERROR, not a valid memget id in work(krel)',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KREL =',KFIRST,KREL,
     *      ' found memory check at work(krel):',MEMCK,
     *      ' expected                        :',MEMID
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMREL ERROR, not a valid memget id in work(krel)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMREL ERROR, not a valid memget id in work(krel)')
C
C     end of MEMREL
C
      END
C  /* Deck memchk */
      SUBROUTINE MEMCHK(TEXT,WORK,KFIRST)
C
C 17-Mar-1989 Hans Joergen Aa. Jensen
C
C Check if memory allocation from WORK(KFIRST) is intact
C (check if anything has been out of bounds).
C
#include "implicit.h"
      DIMENSION     WORK(*)
      CHARACTER*(*) TEXT
C
#include "memcb1.h"
C
Chj1  CALL QENTER('MEMCHK  ')
C
      KALLOC = KFIRST
      IF (KALLOC .EQ. 1) KALLOC = 1 + LENID
C
      IALLOC = 0
  100 CONTINUE
         IALLOC = IALLOC + 1
#ifdef DBG_LENID
         do i = 2,lenid-1
            WMEMCK = WORK(KALLOC-LENID+i)
            IF (MEMCK .NE. MEMID) GO TO 8001
         end do
#else
         WMEMCK = WORK(KALLOC-1)
         IF (MEMCK .NE. MEMID) GO TO 8001
#endif
         WLREAL = WORK(KALLOC-LENID+1)
         IF (LREAL .GE. 0) THEN
            KALLOC = KALLOC + LENID + LREAL
            GO TO 100
         END IF
         IF (LREAL .NE. LZERO) GO TO 8002
C
Chj1  CALL QEXIT('MEMCHK  ')
      RETURN
C
C        error branch IF (MEMCK .NE. MEMID) THEN
C
 8001    WRITE (LUWMEM,'(/A/2A/A,3I10,2(/A,I25))')
     *      ' MEMCHK ERROR, not a valid memget id in work(kalloc-1)',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KALLOC,IALLOC =',KFIRST,KALLOC,IALLOC,
     *      ' found memory checks:',MEMCK,
     *      ' expected           :',MEMID
         WRITE (LUWMEM,'(/A/A)')
     &      ' Dump of allocations before invalid memget id:',
     &      ' IALLOC #    IDENT     KBASE         Length'
         KALLOC = KFIRST
         IF (KALLOC .EQ. 1) KALLOC = 1 + LENID
         DO I = 1, IALLOC
            WRITE (IDENT8,'(A8)') WORK(KALLOC-LENID)
            WLREAL = WORK(KALLOC-LENID+1)
            WRITE (LUWMEM,'(I5,2X,A8,2I15)') I,IDENT8,KALLOC,LREAL
            KALLOC = KALLOC + LENID + LREAL
         END DO
#ifdef DBG_LENID
         do i = 2,lenid-1
            WMEMCK = WORK(KALLOC-LENID+i)
            if (memck .ne. memid) write (luwmem,*)
     &         'memchk error in debug element: ',i+1,memck
         end do
#endif
         CALL QENTER('MEMCHK  ')
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMCHK ERROR, not a valid memget id in work(kalloc-1)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMCHK ERROR, not a valid memget '//
     &             'id in work(kalloc-1)')
C
C        error branch IF (LREAL .LT. 0) THEN
C
 8002    WRITE (LUWMEM,'(/A/2A/A,4I10)')
     *      ' MEMCHK ERROR, next allocation LENGTH is negative',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KALLOC,IALLOC,LENGTH =',KFIRST,KALLOC,IALLOC,LREAL
         CALL QENTER('MEMCHK  ')
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMCHK ERROR, next allocation LENGTH is negative',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMCHK ERROR, next '//
     &             'allocation has negative length.')
C
C     end of MEMCHK
C
      END
C  /* Deck memdump */
      SUBROUTINE MEMDUMP(TEXT,WORK,KFIRST)
C
C June 2015, based on MEMCHK (Hans Joergen Aa. Jensen)
C
C Check if memory allocation from WORK(KFIRST) is intact
C (check if anything has been out of bounds).
C
#include "implicit.h"
      DIMENSION     WORK(*)
      CHARACTER*(*) TEXT
C
#include "memcb1.h"
C
Chj1  CALL QENTER('MEMDUMP ')
C
      KALLOC = KFIRST
      IF (KALLOC .EQ. 1) KALLOC = 1 + LENID
C
      IALLOC = 0
  100 CONTINUE
         IALLOC = IALLOC + 1
#ifdef DBG_LENID
         do i = 2,lenid-1
            WMEMCK = WORK(KALLOC-LENID+i)
            IF (MEMCK .NE. MEMID) GO TO 8001
         end do
#else
         WMEMCK = WORK(KALLOC-1)
         IF (MEMCK .NE. MEMID) CALL MEMCHK(TEXT,WORK,KFIRST)
#endif
         WLREAL = WORK(KALLOC-LENID+1)
         IF (LREAL .GE. 0) THEN
            KALLOC = KALLOC + LENID + LREAL
            GO TO 100
         END IF
         IF (LREAL .NE. LZERO) CALL MEMCHK(TEXT,WORK,KFIRST)
C
         WRITE (LUWMEM,'(/A/2A//A)')
     &      ' *** Dump of memory allocations ***',
     &      '     Text from calling routine: ',TEXT,
     &      ' IALLOC #    IDENT     KBASE         Length'
         KALLOC = KFIRST
         IF (KALLOC .EQ. 1) KALLOC = 1 + LENID
         DO I = 1, IALLOC
            WRITE (IDENT8,'(A8)') WORK(KALLOC-LENID)
            WLREAL = WORK(KALLOC-LENID+1)
            WRITE (LUWMEM,'(I5,2X,A8,2I15)') I,IDENT8,KALLOC,LREAL
            KALLOC = KALLOC + LENID + LREAL
         END DO

Chj1  CALL QEXIT('MEMDUMP ')
      RETURN
C
C     end of MEMDUMP
C
      END
C  /* Deck memopt */
      SUBROUTINE MEMOPT(DATA_TYPE,KBASE,LENGTH,WORK,KFREE,LFREE,ALLOC)
C
C     Almost indentical to MEMGET. Luuk Visscher 21-10-1997.
C     Updated to optimized version. J. Thyssen Aug 25 1998
C
C     Memory allocation
C
#include "implicit.h"
      CHARACTER*(*) DATA_TYPE
      DIMENSION     WORK(*)
      LOGICAL       ALLOC
C
#include "iratdef.h"
#include "memcb1.h"
C
      IF (DATA_TYPE(1:4) .EQ. 'REAL') THEN
         LREAL = LENGTH
      ELSE IF (DATA_TYPE(1:4) .EQ. 'INTE') THEN
         LREAL = (LENGTH-1)/IRAT + 1
      ELSE IF (DATA_TYPE(1:4) .EQ. 'LOGI') THEN
         LREAL = (LENGTH-1)/LRAT + 1
      ELSE IF (DATA_TYPE(1:4) .EQ. 'INT2') THEN
         LREAL = (LENGTH-1)/IRAT2 + 1
      ELSE IF (DATA_TYPE(1:4) .EQ. 'WORK') THEN
         LENGTH = LFREE - 2*LENID
C                         ^- make sure OK for KFREE .eq. 1
         LREAL  = LENGTH
      ELSE
         CALL QENTER('MEMOPT  ')
         WRITE (LUWMEM,'(/2A)')
     *      ' MEMOPT ERROR, illegal type :',DATA_TYPE
         CALL QTRACE(LUWMEM)
         IF (LUEMEM .NE. LUWMEM) THEN
            WRITE (LUEMEM,'(/2A)')
     *         ' MEMOPT ERROR, illegal type :',DATA_TYPE
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMOPT ERROR, illegal type')
      END IF
C
      IDENT8 = '-MEMOPT-'
C
      IF (ALLOC) THEN
         IF (KFREE .NE. 1) THEN

C           length error check:
            IF (LREAL+LENID .GT. LFREE) GO TO 8001

C            Error check, memory check of WORK(KFREE) OK ? :

            WMEMCK = WORK(KFREE-1)
            IF (MEMCK .NE. MEMID) GO TO 8002
C           ... comparing integers, to be sure no renormalization of float

C           Warning check, is length of WORK(KFREE) .ne. 0 ? :

            WMEMCK = WORK(KFREE-LENID+1)
            IF (MEMCK .NE. LZERO) GO TO 8003
C           ... go issue warning and continue at 8013
         ELSE
C           KFREE = 1 case, initialize memory vector info:
            IF (LREAL+2*LENID .GT. LFREE) GO TO 8001
            WORK(3) = WMEMID
            KFREE = KFREE + LENID
            LFREE = LFREE - LENID
         END IF
      ELSE
         IF (KFREE .EQ. 1) THEN
            KFREE = KFREE + LENID
            LFREE = LFREE - LENID
         END IF
      END IF

C     save ident and length for new allocation

 8013 CONTINUE
      IF (ALLOC) THEN
         READ (IDENT8,'(A8)') WORK(KFREE-LENID)
         WORK(KFREE-LENID+1) = WLREAL
      END IF

      KBASE  = KFREE
      KFREE  = KBASE + (LREAL + LENID)
      LFREE  = LFREE - (LREAL + LENID)

      IF (ALLOC) THEN
         WORK(KFREE-LENID)   = WIDENT_KFREE
         WORK(KFREE-LENID+1) = WLZERO
         WORK(KFREE-1)       = WMEMID
      END IF
      RETURN
C
C     error branches:
C
C        error branch  IF (LREAL+LENIDS .GT. LFREE) :
C
 8001    CALL QENTER('MEMOPT  ')
         IF (KFREE .GT. 1 .AND. ALLOC) THEN
            LENIDS = LENID
            WRITE (LUWMEM,1010) LREAL+LENIDS,LFREE
            CALL MEMCHK('MEMOPT ERROR (insuff. memory)',WORK,1)
         ELSE
            LENIDS = 2*LENID
            WRITE (LUWMEM,1010) LREAL+LENIDS,LFREE
         END IF
         CALL QTRACE(LUWMEM)
         IF (LUEMEM .NE. LUWMEM) THEN
            WRITE (LUEMEM,1010) LREAL+LENIDS,LFREE
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMOPT ERROR, insufficient work space in memory')
 1010 FORMAT(
     &   /' MEMOPT ERROR, insufficient free space for next allocation',
     &   /T16,'( Need:',I10,', available (LFREE):',I10,' )')
C
C        error branch IF (MEMCK .NE. MEMID) :
C
 8002    CALL QENTER('MEMOPT  ')
         WRITE (LUWMEM,'(/A/A,I10,2(/A,I25))')
     *      ' MEMOPT ERROR, not a valid MEMOPT id in work(kfree)',
     *      ' KFREE =',KFREE,
     *      ' found memory check :',MEMCK,
     *      ' expected           :',MEMID
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMOPT ERROR, not a valid MEMOPT id in work(kfree)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMOPT ERROR, not a valid MEMOPT '//
     &             'id in work(kfree)')
C
C        warning branch IF (MEMCK(2) .NE. 0) :
C
 8003    NWNMEM = NWNMEM + 1
         WRITE (LUWMEM,'(/A/A,I10,/A,I25)')
     *      ' MEMOPT WARNING, nonzero allocation work(kfree)',
     *      ' KFREE =',KFREE,
     *      ' found memory allocation :',MEMCK
         CALL QENTER('MEMOPT  ')
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMOPT WARNING, nonzero allocation work(kfree)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QEXIT('MEMOPT  ')
C        Warning issued, continue and do the assignment now
         GO TO 8013
C
C     end of MEMOPT
C
      END
C --- end of mempgk.F ---
