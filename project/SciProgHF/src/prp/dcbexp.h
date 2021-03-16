!C
!C
!C*****************************************************************************
!C  DCBEXP:  Expectation values
!C
      INTEGER       MAXEXP
      PARAMETER(MAXEXP=350)

      INTEGER       IPREXP,NEXPP,LEXPP,IPDIPOLE,IPQUADRU,IPEFG,         &
     &              IPPVC,IPEFT,IPEFTNTL,IPDSUSC,IPEFFDEN,IPEFN
      LOGICAL ORBANA,PRJANA,PRPCAN

      COMMON/DCBEXP/IPREXP,NEXPP,LEXPP(MAXEXP),                         &
     &              IPDIPOLE(3),IPQUADRU(6),IPEFG(6*MXCENT),            &
     &              IPPVC(MXCENT),IPEFT(15*MXCENT),IPEFTNTL(15*MXCENT), &
     &              IPDSUSC(39),IPEFFDEN(MXCENT),                       &
     &              IPEFN(3*MXCENT),ORBANA,PRJANA,PRPCAN

