!*** Unpack four-indices ***
      integer ibit2,IBIT2M,IBIT4, IBIT4M
#if defined (SYS_CRAY) || defined (INT_STAR8)
      PARAMETER (IBIT2=65536,IBIT2M=65535,IBIT4=4294967296,             &
     & IBIT4M=4294967295)
#else
      PARAMETER (IBIT2 = 256,IBIT2M = 255,IBIT4 = 65536,IBIT4M = 65535)
#endif
#include <ibtfun.h>
#if defined (SYS_CRAY) || defined (INT_STAR8)
      INDGET(I,J) = IBTAND(IBTSHR(I,J),IBIT2M)
      LGET(I) = INDGET(I,0)
      KGET(I) = INDGET(I,16)
      JGET(I) = INDGET(I,32)
      IGET(I) = INDGET(I,48)
      IPACK(I,J,K,L) = IBTOR(IBTOR(IBTOR(IBTSHL(L,0),                   &
     &IBTSHL(K,16)),IBTSHL(J,32)),IBTSHL(I,48))
#else
      INDGET(I,J) = IBTAND(IBTSHR(I,J),IBIT2M)
      LGET(I) = INDGET(I,0)
      KGET(I) = INDGET(I,8)
      JGET(I) = INDGET(I,16)
      IGET(I) = INDGET(I,24)
      IPACK(I,J,K,L) = IBTOR(IBTOR(IBTOR(IBTSHL(L,0),                   &
     &IBTSHL(K,8)),IBTSHL(J,16)),IBTSHL(I,24))
#endif
