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

#ifndef PRG_DIRAC
C  /* Deck tstdir */
      SUBROUTINE TSTDIR(FMAT,DMAT,NDMAT,WORK,LWORK,IPRINT)
C
C     Written by Henrik Koch and Trygve Helgaker 27-November-1991.
C
C     PURPOSE : Driver routine for the test calculation of the
C               two-electron part of the fock matrices.
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      PARAMETER (MAXDIS=500)
      DIMENSION DMAT(*), FMAT(*), WORK(LWORK)
#include "inforb.h"
C
      CALL TIMER('START ',TIMSTR,TIMEND)
      KGMAT = 1
      KINDX = KGMAT + NBAST*NBAST*MAXDIS
      KEND  = KINDX + (MAXDIS + 1)/IRAT
      LFREE = LWORK - KEND
      IF (KEND .GT. LWORK) CALL STOPIT('TSTDIR',' ',KEND,LWORK)
C
      NUMDIS = -1
    1 CONTINUE
C
C     Get distributions
C
      CALL GETDIS(WORK(KGMAT),WORK(KINDX),NUMDIS,MAXDIS,WORK(KEND),
     &            LFREE)
      IF (NUMDIS .EQ. -1) GOTO 2
C
C     Add contributions to Fock matrix
C
      CALL DISFCK(FMAT,DMAT,NDMAT,WORK(KGMAT),WORK(KINDX),NUMDIS)
C
      GOTO 1
C
    2 CONTINUE
C
C     Fock matrix has now been constructed
C
C------------------------------------------------------
C     Write out densities and associated fock matrices.
C------------------------------------------------------
C
      IF (IPRINT.GT.2) THEN
         CALL HEADER('Density and Fock matrices from '//
     &               'distributions in TSTDIR',-1)
         DO 300 I = 1, NDMAT
            ISTR = NBAST*NBAST*(I - 1) + 1
            WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
            CALL OUTPUT(DMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
      CALL TIMER('GETDIS',TIMSTR,TIMEND)
      RETURN
      END
C  /* Deck disfck */
      SUBROUTINE DISFCK(FMAT,DMAT,NDMAT,GMAT,INDEX,NUMDIS)
C
C     Written by Henrik Koch and Trygve Helgaker 27-November-1991.
C
C
C     PURPOSE : Calculate fock matrices using distributions.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (ONE=1.0D00, HALF=0.5D00, FOURTH=0.25D00)
      INTEGER   P,Q,R,S
      DIMENSION FMAT(NBAST,NBAST,NDMAT), DMAT(NBAST,NBAST,NDMAT)
      DIMENSION GMAT(NBAST,NBAST,NUMDIS), INDEX(NUMDIS)
#include "inforb.h"
C
      DO 100 I = 1,NUMDIS
C
         CALL UNPKIJ(INDEX(I),R,S)
C
         IF (R .EQ. S) THEN
            FACCOU = HALF
            FACEXC = FOURTH
         ELSE
            FACCOU = ONE
            FACEXC = HALF
         END IF
C
         DO 110 J = 1,NDMAT
C
            VALUE = FACCOU*DDOT(NBAST*NBAST,DMAT(1,1,J),1,GMAT(1,1,I),1)
            FMAT(R,S,J) = FMAT(R,S,J) + VALUE
            FMAT(S,R,J) = FMAT(S,R,J) + VALUE
C
            DO 120 Q = 1,NBAST
C
               FMAT(R,Q,J) = FMAT(R,Q,J)
     &                  - FACEXC*DDOT(NBAST,DMAT(1,S,J),1,GMAT(1,Q,I),1)
               FMAT(S,Q,J) = FMAT(S,Q,J)
     &                  - FACEXC*DDOT(NBAST,DMAT(1,R,J),1,GMAT(1,Q,I),1)
C
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck getdis */
      SUBROUTINE GETDIS(GMAT,INDXAB,NUMDIS,MAXDIS,WORK,LWORK)
C
C     Written by Henrik Koch and Trygve Helgaker 26-November-1991.
C
C     PURPOSE : Driver routine for the calculation of the two-electron
C               distributions (**|cd).
C               The distributions are stored as full squares and without
C               symmetry reduction.
C
C                      *CALL CBITWO   from TWOINP
C                      *CALL INFORB   from sirius
C                      *CALL GAMCOM   from READIN
C                      *CALL CCOM     from READIN
C                      *CALL PINCOM   from READIN
C                      *CALL PRIMIT   from READIN
C                      *CALL XYZPOW   from READIN
C                      *CALL SHELLS   from READIN
C                      *CALL NUCLEI   from READIN
C                      *CALL SYMMET   from READIN
C                      *CALL DORPS    from abatro.u
C                       ( only DOREPS and DOCOOR used)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "dummy.h"
      LOGICAL ABA
      DIMENSION GMAT(NBAST,NBAST,MAXDIS), INDXAB(MAXDIS)
      DIMENSION WORK(LWORK)
#include "inforb.h"
#include "cbitwo.h"
#include "dorps.h"
C
      CALL ABARUN(ABA)
C
      I2TYP = 0
C
C-------------------------------------------------------------
C     Setup information for the two-electron integralroutines.
C-------------------------------------------------------------
C
      IF (.NOT.ABA) THEN
         NODV = .FALSE.
         NOPV = .FALSE.
         IF (NASHT .EQ. 0) NODV = .TRUE.
         IF (NASHT .LT. 2) NOPV = .TRUE.
         IPRALL =  0
         DO 100 I = 0,7
            DOREPS(I) = .TRUE.
  100    CONTINUE
         DO 110 I = 1,MXCENT
            DOCOOR(1,I) = .TRUE.
            DOCOOR(2,I) = .TRUE.
            DOCOOR(3,I) = .TRUE.
  110    CONTINUE
         IPRNTA = 0
         IPRNTB = 0
         IPRNTC = 0
         IPRNTD = 0
         RETUR  = .FALSE.
         NOCONT = .FALSE.
         TKTIME = .FALSE.
      END IF
C
      MAXDIF = 0
      ITYPE  = 4
      KJSTRS = 1
      KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KIORBS = KNCONT + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KJORBS = KIORBS + (MXSHEL*MXAOVC + 1)/IRAT
      KKORBS = KJORBS + (MXSHEL*MXAOVC + 1)/IRAT
      KLAST  = KKORBS + (MXSHEL*MXAOVC + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('GETDIS','PAOVEC',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,.FALSE.,IPRALL)
      KLAST = KJORBS
      LWRK  = LWORK - KLAST + 1
C
C----------------------------
C     Calculate distributions
C----------------------------
C
      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL TWOINT(WORK(KLAST),LWRK,DUMMY,DUMMY,DUMMY,NDMAT,
     &            IDUMMY,IDUMMY,GMAT,
     &            INDXAB,NUMDIS,MAXDIS,ITYPE,MAXDIF,0,NODV,NOPV,NOCONT,
     &            TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &            IDUMMY,I2TYP,WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &            WORK(KIORBS),
     &            IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,.FALSE.)
      CALL TIMER('TWOINT',TIMSTR,TIMEND)
C
C--------------------------------------
C     Write out integrals distributions
C--------------------------------------
C
      IF (IPRINT.GT.3) THEN
         CALL HEADER('Distribution matrices in GETDIS',-1)
         DO 300 I = 1, NUMDIS
            CALL UNPKIJ(INDXAB(I),IA,IB)
            WRITE(LUPRI,'(//,1X,A,I4)')' Distribution matrix No.',I
            WRITE(LUPRI,'(1X,A,2I4)')  ' Orbital indices:       ',IA,IB
            CALL OUTPUT(GMAT(1,1,I),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck fcktes */
#if !defined (VAR_GETFCK)
      SUBROUTINE FCKTES(WORK,LWORK,MAXDIF,NODV,NOPV,NOCONT,TKTIME,
     &                  IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR)
C     April 1995: FCKTES below must be changed to call new Fock matrix
C     routines instead of GETFCK.
      CALL QUIT('FCKTES not implemented in this version!')
      END
#else
      SUBROUTINE FCKTES(WORK,LWORK,MAXDIF,NODV,NOPV,NOCONT,TKTIME,
     &                  IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR)
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D1 = 1.0D0)
      LOGICAL NODV, NOPV, NOCONT, TKTIME, RETUR, DIRECT
      DIMENSION WORK(LWORK)
#include "inforb.h"
C
      IF (NODV) THEN
         NDMAT = 1
      ELSE
         NDMAT = 2
      END IF
      KFMAT = 1
      KDMAT = KFMAT + N2BASX*NDMAT
      KCMO  = KDMAT + N2BASX*NDMAT
      KDV   = KCMO  + NCMOT
      KDTSO = KDV   + NNASHX
      KDASO = KDTSO + NNBAST
      KLAST = KDASO + NNBAST
      IF (KLAST .GT. LWORK) CALL STOPIT('FCKTES',' ',KLAST,LWORK)
      CALL ONEDSF(WORK(KCMO),WORK(KDV),WORK(KDMAT),
     &            WORK(KDMAT+N2BASX),IPRINT,NODC,NODV)
      CALL DAXPY(N2BASX,-D1,WORK(KDMAT+N2BASX),1,WORK(KDMAT),1)
      CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KDTSO),NNBAST)
      CALL SQDENS(WORK(KDTSO),WORK(KFMAT),IPRINT)
      IF (.NOT.NODV) CALL DZERO(WORK(KFMAT+N2BASX),N2BASX)
C
      DIRECT = .TRUE.
      KLAST = KCMO
      LWRK  = LWORK - KLAST + 1
      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL GETFCK(WORK(KFMAT),WORK(KDMAT),NDMAT,WORK(KLAST),LWRK,DIRECT,
     *          .FALSE.)
      CALL TIMER('GETFCK',TIMSTR,TIMEND)
      if (.true.) return
C
C     Test distributions
C
      KFMAT = 1
      KDMAT = KFMAT + N2BASX*NDMAT
      KCMO  = KDMAT + N2BASX*NDMAT
      KDV   = KCMO  + NCMOT
      KDTSO = KDV   + NNASHX
      KDASO = KDTSO + NNBAST
      KLAST = KDASO + NNBAST
      IF (KLAST .GT. LWORK) CALL STOPIT('FCKTES',' ',KLAST,LWORK)
      CALL ONEDSF(WORK(KCMO),WORK(KDV),WORK(KDMAT),WORK(KDMAT+N2BASX),
     &            IPRINT,NODC,NODV)
      CALL DAXPY(N2BASX,-D1,WORK(KDMAT+N2BASX),1,WORK(KDMAT),1)
      CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KDTSO),NNBAST)
      CALL SQDENS(WORK(KDTSO),WORK(KFMAT),IPRINT)
      IF (.NOT.NODV) CALL DZERO(WORK(KFMAT+N2BASX),N2BASX)
C
      KLAST = KCMO
      LWRK  = LWORK - KLAST + 1
      CALL TSTDIR(WORK(KFMAT),WORK(KDMAT),NDMAT,WORK(KLAST),LWRK,IPRINT)
      RETURN
      END
#endif
C  /* Deck sqdens */
      SUBROUTINE SQDENS(DTSO,DMAT,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION DTSO(*), DMAT(NBAST,NBAST)
#include "inforb.h"
#include "symmet.h"
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from SQDENS','*',103)
      CALL DZERO(DMAT,NBAST*NBAST)
      ISOFF = 0
      DO 100 IREP = 1, NSYM
         DO 200 I = 1, NBAS(IREP)
            INDA = IBAS(IREP) + I
            DO 300 J = 1, I
               INDB = IBAS(IREP) + J
               ISOFF = ISOFF + 1
               DMAT(INDA,INDB) = DMAT(INDA,INDB) + DTSO(ISOFF)
               DMAT(INDB,INDA) = DMAT(INDA,INDB)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Squared matrix in SQDENS',-1)
         CALL OUTPUT(DMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
#endif /* ifndef PRG_DIRAC */
C  /* Deck sether */
      SUBROUTINE SETHER(JPRINT,NEWGEO,RELCAL,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  CBIREA: LUMLCL
C  SYMMET: MAXREP,NAOS(8))
C  INFORB: NSYM,NBAS(8)
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "cbirea.h"
#include "symmet.h"
#include "inforb.h"
C
      LOGICAL SFIRST, RELCAL, TSTINP, HERMIT, NEWGEO
      SAVE    SFIRST
      DATA    SFIRST /.TRUE./
C
C     Initialize /CBIREA/ and call readin
C
      IF (SFIRST .OR. NEWGEO) THEN
         TSTINP = .FALSE.
         CALL REAINI(JPRINT,RELCAL,TSTINP)
         LUMLCL = 56
C
         IPREAD = JPRINT
         HERMIT = .FALSE.
         CALL READIN(HERMIT)
C
C        consistency check
C
         IF(.NOT.RELCAL) THEN
           IF (MAXREP+1 .NE. NSYM) THEN
            WRITE (LUPRI,'(2(/A,I5))')
     &      ' SETHER fatal error: NSYM  from SIRIUS is',NSYM,
     &      '                  MAXREP+1 from READIN is',MAXREP+1
            CALL QUIT('SETHER error: NSYM .ne. MAXREP+1')
           END IF
           NERR = 0
           DO 100 ISYM = 1,NSYM
              IF (NAOS(ISYM) .NE. NBAS(ISYM)) NERR = NERR + 1
  100      CONTINUE
           IF (NERR .GT. 0) THEN
              WRITE (LUPRI,'(/A/A/)')
     &        ' SETHER fatal error: NBAS(:) from LUONEL',
     &        '   is different from NAOS(:) from READIN.'
              WRITE (LUPRI,'(A,8I5)') ' NBAS(:) =',(NBAS(I),I=1,NSYM)
              WRITE (LUPRI,'(A,8I5)') ' NAOS(:) =',(NAOS(I),I=1,NSYM)
              CALL QUIT('SETHER error: NBAS(:) .ne. NAOS(:)')
           END IF
         ELSE
C
C        In relativistic calculation we must generate what is needed in INFORB
C
           CALL RELINF
         ENDIF
         SFIRST = .FALSE.
      END IF
      RETURN
      END
C  /* Deck aindex */
#ifdef PRG_DIRAC
      SUBROUTINE AINDEX(ISHELA,NAINTS,INDEXA,DOINDX,IPRINT)
#else /* PRG_DIRAC */
      SUBROUTINE AINDEX(ISHELA,NAINTS,INDEXA,DOINDX,IORBSH,IPRINT)
#endif /* PRG_DIRAC */
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#include "mxcent.h"
      LOGICAL DOINDX
#ifdef PRG_DIRAC
      DIMENSION INDEXA(*)
#else /* PRG_DIRAC */
      DIMENSION INDEXA(*), IORBSH(MXSHEL,MXAOVC)
#endif /* PRG_DIRAC */
#include "blocks.h"
#include "symmet.h"
#include "ibtfun.h"
C
      IORB = 0
      MULA   = ISTBSH(ISHELA)
      DO 100 IORBA = 1, NORBSH(ISHELA)
         DO 200 ICOMPA = 1, KHKTSH(ISHELA)
            ITYNA = ISYMAO(NHKTSH(ISHELA),ICOMPA)
            IADR  = IORBSB(IORBSH(ISHELA,1))
     &            + KHKTSH(ISHELA)*(IORBA - 1) + ICOMPA
            DO 300 IREPA = 0, MAXREP
               IF (IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
                  IORB = IORB + 1
                  IF (DOINDX) INDEXA(IORB) = IPTSYM(IADR,IREPA)
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      NAINTS = IORB
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from AINDEX',2)
         IF (.NOT.DOINDX) THEN
            WRITE (LUPRI,'(2X,A,I5)') ' Number of a orbitals:',NAINTS
         ELSE
            WRITE (LUPRI,'(2X,I3,A,8I5/,(40X,8I5))')
     &          NAINTS,' distributions in this TWOINT call:',
     &            (INDEXA(I),I=1,NAINTS)
         END IF
      END IF
      RETURN
      END
#ifndef PRG_DIRAC
C  /* Deck abarun */
      SUBROUTINE ABARUN(RUNABA)
#include "implicit.h"
      LOGICAL ABA, RUNABA
      SAVE ABA
      DATA ABA /.FALSE./
      RUNABA = ABA
      RETURN
      ENTRY ABASET
      ABA = .TRUE.
      RETURN
      END
#endif /* ifndef PRG_DIRAC */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck RELINF */
      SUBROUTINE RELINF
C*****************************************************************************
C
C     Generate INFORB for DIRAC
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "nuclei.h"
#include "inforb.h"
#include "symmet.h"
C
      NSYM  = MAXREP + 1
      NBAST = NBASIS
      N2BASX = NBAST*NBAST
      DO I = 1,NSYM
        NBAS(I) = NAOS(I)
        IBAS(I) = ICOS(I,1)
      ENDDO
C
      RETURN
      END
