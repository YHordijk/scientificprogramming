! dcbmp2no.h -- info for MP2 natural orbitals
      CHARACTER VECMP2NO*72
      CHARACTER VECNOSEL*72
      CHARACTER INFO_STRING*62
      CHARACTER MP2NOFN*14 
      CHARACTER SELPOPNAT*6
!     max. number of virt. spinors in NO-MP2 run
      INTEGER JMAXVSP,IMO_NAT,IOCC_NAT,NACTMP2NO,
     &        LMP2_DEN_INP,LUNATORB,IAOMOSCHEME,
     &        IF2MP2NO,IMAXVSP,INTMP2NO,IPRMP2NO,
     &        MP2NOFN_L
      PARAMETER (JMAXVSP = 200)
!     dimension of transform. matrix : IMO_NAT
!
      LOGICAL         DEF_EN,RUN_CCMOD,MP2_DENRUN,SKIP_MP2DEN,          &
     &                SKIP_ITRAFO,MP2NATPOP,MP2NATPOP2
      COMMON/CBLMP2NO/DEF_EN,RUN_CCMOD,MP2_DENRUN,SKIP_MP2DEN,          &
     &                SKIP_ITRAFO,MP2NATPOP,MP2NATPOP2
      COMMON/CBIMP2NO/IPRMP2NO,INTMP2NO,IF2MP2NO,IMAXVSP,IAOMOSCHEME,   &
     &                LMP2_DEN_INP,MP2NOFN_L,IMO_NAT,IOCC_NAT,          &
     &                LUNATORB,NACTMP2NO(2)
      COMMON/CBCMP2NO/VECMP2NO(2),VECNOSEL(2),MP2NOFN,INFO_STRING(3),   &
     &                SELPOPNAT
! -- end of dcbmp2no.h --

