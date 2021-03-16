! File: luciarel/wrkspc.inc
#ifdef LUCIAREL_STANDALONE
      COMMON/BIGGY/WORK(MXPWRD)
#else
!     WORK will be referenced out of bounds, this is not an error.
!     In Dirac the arrays are really allocated in the
!     global WRK array allocated outside LUCIAREL by using
!     the off-set between the two memory locations. See
!     subroutine COMPOFF for more information.

!     Note: all Ksomething used for WORK(Ksomething) in LUCIAREL
!     should be declared integer*8 in luciarel/ipoist8.inc

      COMMON/BIGGY/WORK(1)
#endif
