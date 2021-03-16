! File: lucita/wrkspc.inc
!
!     WORK will be referenced out of bounds, this is not an error.
!     In Dirac the arrays are really allocated in the
!     global WRK array allocated outside LUCITA by using
!     the off-set between the two memory locations. See
!     subroutine COMPOFF for more information.

      COMMON/BIGGY/WORK(1)
