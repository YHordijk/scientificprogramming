! FILE : dcbkrmc_itinfo.h
!
!     i     IINFO      DINFO
!    ---   -------    -------
!     01    ITMAC      EMY
!     02    ITMIC      EACTIV
!     03    NZRED      EMCSCF
!     04    NCRED      DEPRED
!     05    NERED      DEACT
!     06    NPRED      RATIO
!     07               BETA
!     08               STPLNG
!     09               RTRUST
!     10               GAMMA
!     11               GCINRM ( GNORM(1) )
!     12               GEENRM ( GNORM(2) )
!     13               GEPNRM ( GNORM(3) )
!     14               GOBNRM ( GNORM(4) )
!     15               GTOTNRM( GNORM(5) )
!     16               TIMMAC
!     17               TIMMIC
!     18               TIMITR
!     19               CI STPLNG
!     20               E STPLNG
!     21               P STPLNG
!     22               LEVEL SHIFT
!
      INTEGER LDINFO, LIINFO
      PARAMETER (LDINFO = 22, LIINFO = 6)

      INTEGER IINFO(LIINFO)
      COMMON /IITINFO/ IINFO

      REAL*8 DINFO(LDINFO)
      COMMON /RITINFO/ DINFO
