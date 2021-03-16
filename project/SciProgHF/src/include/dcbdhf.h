!
! FILE    : dcbdhf.h
!
!*** Dirac common blocks for Dirac-Hartree-Fock calculations ***
!
!     PROJEC     - any projection, defined by a set of coefficients
!     ... (descriptions are to be added)
!
!     ERGCNV2, EVCCNV2, FCKCNV2, MAXITR2 - related to the preliminary SCF method before turning to 4c mode
!
!     BSS-SCF (...MI/sept05 - can be extended by more SCF control flags )
!


      INTEGER MFMAT, MXOPEN, MXMOFREEZE
      PARAMETER (MFMAT = 10, MXOPEN = 10, MXMOFREEZE = 100)

      INTEGER MFRAG
      PARAMETER (MFRAG = 6)

      INTEGER MAX_BOS_BL
      PARAMETER (MAX_BOS_BL = 72)
!     MAX_BOS_BL must be equal to MAX_SUB_BL in dcborb.h /hjaaj oct 2004

      CHARACTER*12 DHF_INTTYP
      CHARACTER*8  CACC
      CHARACTER*6  PRJFIL(MFRAG)
      CHARACTER*72 VCPROJ(2, MFRAG), VCFROZ(2)

      COMMON /DCCDHF/ DHF_INTTYP, CACC, PRJFIL, VCPROJ, VCFROZ

      LOGICAL DHFCONV(2), DHFEXIT,    TRIVEC, TRIFCK, ERGCNV, EVCCNV,   &
     &        FCKCNV, DODSCF, DODAMP, DODIIS, DIISON, FIXDIF, SIRIFC,   &
     &        DIISAO, DIISMO, OVLSEL, DYNSEL, DOBOSSEL,AUTOCC,INIOCC,   &
     &        DHF_SKIPEE, DHF_SKIPEP, NOQCDHF,AOC,    DOCCNV, SCFPOP,   &
     &        ERGCNV2, EVCCNV2, FCKCNV2, DOMOFREEZE,  DOLEVEL,ATOMST,   &
     &        WRITE_FMO_MATRIX,ATHUCK,NOSWIT

      COMMON /DCLDHF/ ERGCNV, EVCCNV, FCKCNV, DHFCONV, DHFEXIT, TRIVEC, &
     &                TRIFCK, DODSCF, DODAMP, DODIIS,  DIISON,  OVLSEL, &
     &                FIXDIF, DIISAO, DIISMO, DYNSEL,  DOBOSSEL,AUTOCC, &
     &                INIOCC, DHF_SKIPEE, DHF_SKIPEP,  NOQCDHF, AOC,    &
     &                ERGCNV2,EVCCNV2,FCKCNV2,DOMOFREEZE,       DOCCNV, &
     &                DOLEVEL,ATOMST, SIRIFC, SCFPOP, WRITE_FMO_MATRIX, &
     &                ATHUCK,NOSWIT

      INTEGER NITER, MAXITR, IPRSCF, KITER, MXDIIS, ITDIIS, NELMBM,     &
     &        INTDEF, INTDEF_SAVE, INTFLG,                              &
     &        INTBUF, ITRINT(2), IPREIG,NSMOTQ, ISMOQ(2),               &
     &        ITRSCF(0:15), NFMAT, ISYMOP(MFMAT), IFCKOP(MFMAT),        &
     &        IHRMOP(MFMAT), NEWOCC, MXMACRO, MXMICRO, MAXBCK, NFRAG,   &
     &        NTMO(2), I2TMOT(2), N2TMO(2),N2TMOTQ, NISHMF(2),NESHMF(2),&
     &        NACSHMF(2, MXOPEN), NASHMFT, NOCCMF(2), NPSHMF(2),        &
     &        NPRJNUC(MFRAG), MAXITR2,                                  &
     &        NMOFREEZE(2),  IMOFREEZE(MXMOFREEZE)

      COMMON /DCIDHF/ NITER, MAXITR, IPRSCF, KITER, MXDIIS, ITDIIS,     &
     &                NELMBM, INTDEF, INTDEF_SAVE, INTFLG, INTBUF,      &
     &                ITRINT, IPREIG, NSMOTQ, ISMOQ, ITRSCF, NFMAT,     &
     &                ISYMOP, IFCKOP, IHRMOP, NEWOCC, MXMACRO, MXMICRO, &
     &                MAXBCK, NFRAG, NTMO, I2TMOT, N2TMO, N2TMOTQ,      &
     &                NISHMF, NACSHMF, NASHMFT, NOCCMF, NESHMF,NPSHMF,  &
     &                NPRJNUC, MAXITR2, NMOFREEZE, IMOFREEZE

!..number of closed and open shell orbitals with frozen ones subtracted

      REAL*8 SCFCNV(2), ERGVAL, EVCVAL, FCKVAL, DHFERG, ELERGY, E1PART, &
     &       E2PART, ERGBUF, CONVRG, TDF2, TDDG, DNSDHF, DAMPFC, DIISTH,&
     &       BMCOND, CNVINT(2), SCFTID(0:15), ESOLVE, ESOLVN, PRJTHR,   &
     &       SCFCNV2(2), DLSHIF, OPEN_FAC

      COMMON /DCRDHF/ SCFCNV, ERGVAL, EVCVAL, FCKVAL, DHFERG, ELERGY,   &
     &                E1PART, E2PART, ERGBUF, CONVRG, TDF2, TDDG,       &
     &                DNSDHF, DAMPFC, DIISTH, BMCOND, CNVINT, SCFTID,   &
     &                ESOLVE, ESOLVN, PRJTHR, SCFCNV2, DLSHIF, OPEN_FAC

      REAL*8 DF(0:MXOPEN), DA(0:MXOPEN), DALPHA(0:MXOPEN), OLEV(MXOPEN)

      INTEGER NISH_DHF(2), NASH_DHF(2), NSSH_DHF(2), NESH_DHF(2),       &
     &        NPSH_DHF(2), NFRO_DHF(2), NELEC_DHF(2), NELECT_DHF,       &
     &        NAELEC_DHF, NACSH(2, MXOPEN), NOPEN,                      &
     &        NISH_BOS(MAX_BOS_BL), NACSH_BOS(MAX_BOS_BL,MXOPEN)
      COMMON /DCODHF/ DF, DA, DALPHA, OLEV, NISH_DHF, NASH_DHF,         &
     &                NSSH_DHF, NESH_DHF, NPSH_DHF, NFRO_DHF,           &
     &                NELEC_DHF, NELECT_DHF, NAELEC_DHF, NACSH, NOPEN,  &
     &                NISH_BOS, NACSH_BOS

      INTEGER         LUCYCL, LUFCK1, LUFCK2, LUFOCK, LUDIIS, LUDENS,   &
     &                LUEVEC, LUCMOS, LUSMOS, LUFCKT
      COMMON /DHFIOU/ LUCYCL, LUFCK1, LUFCK2, LUFOCK, LUDIIS, LUDENS,   &
     &                LUEVEC, LUCMOS, LUSMOS, LUFCKT
! --- end of dcbdhf.h ---
