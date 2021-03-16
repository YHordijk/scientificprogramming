#if defined (SYS_ALLIANT)
!VD$ NOSYNC
#endif
#if !defined (SYS_ALLIANT)
!parallelization note: no synchronization problems
#include <ivdep.h>
#endif
