C
C luciarel/memory_r.inc : memory information for LUCIAREL
C only used in luciarel/memory_r.F, which contains the memory manager routines for LUCIAREL
C
      INTEGER     NPAD, MAXLVL, MAXMRK

C MAXLVL : maximum number of simultaneous memory allocations in LUCIAREL,
C          can become quite big if you ask for many properties for many states.
C MAXMRK : maximum number of memory marks for releasing blocks of memory,
C          generally called at most once in each subroutine in traceback stack, thus
C          needed number is much smaller than MAXLVL, and 100 is probably always enough.
C ---> To user: increase MAXLVL or MAXMRK in next line if Dirac output tells you to do so.
      PARAMETER(MAXLVL = 5000, MAXMRK = 100)
C
      PARAMETER(NPAD = 1 )
      INTEGER*8   NWORD,KFREES,KFREEL,IBASE,IMARK
      INTEGER     LENGTH, NS, NL, NM
      INTEGER     MARKL, MARKS
      CHARACTER*6 IIDENT, MARKC, MARK_ACT

      COMMON/CMEMOR/NWORD,KFREES,KFREEL,IBASE(MAXLVL),IMARK(MAXMRK),   ! integer*8
     &              LENGTH(MAXLVL),NS,NL,NM,                           ! default integer
     &              MARKL(MAXMRK),MARKS(MAXMRK),                       ! default integer
     &              IIDENT(MAXLVL),MARKC(MAXMRK),MARK_ACT              ! character*6

! local variables for use in subroutines using /CMEMOR/
      INTEGER*8   KBASE,JBASE,JBASEN
      CHARACTER*6 TASK, IDENT

C -- end of luciarel/memory_r.inc
