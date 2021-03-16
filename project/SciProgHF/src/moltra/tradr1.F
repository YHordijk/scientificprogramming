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

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDR1T */
      SUBROUTINE TRDR1T(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,
     &                  NSTR,ANTIS,LMP2,KINDX,KQ,KE,KIBE,DINTSKP)
C
C     Written by Luuk Visscher, Februari 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "dgroup.h"
      PARAMETER (MAXTSK=1)
      LOGICAL ANTIS,LMP2
      DIMENSION NSTR(2,0:2,4),NSTR1(2,0:2,4,MAXTSK)
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4),KQ(2,4),KE(2,4),KIBE(2,4)
      DIMENSION NSTRS(2,4,MAXTSK),ICMOQR1(2,4)
      DIMENSION WORK(*),DINTSKP(*)
C
      CALL QENTER('TRDR1T')
      KFRSAV = KFREE
C
C     Divide index transformation in task with different batches
C     of spinors
C
C     CALL DRV1T (WORK,KFREE,LFREE,MAXTSK,NSTR,NTASK,NSTR1,NSTRS)
C
C     NOTE 1 TODO :
C     DRV1T remains to be written, number of tasks will depend on memory
C     availabe and required for the three-index transformed integrals
C     At the moment : assume that we need only one run through.
C     
      ntask = 1
      call icopy (24,NSTR(1,0,1),1,NSTR1(1,0,1,1),1)
      call icopy (8,1,0,NSTRS,1)
C
C     NOTE 2 TODO :
C     Obviously the best place to parallelize this scheme is over this
C     loop, if the number of task is sufficiently large compared to the 
C     number of processors. More fine-grain parallelization is possible
C     by parallelizing over the calls to CALDIS and MS4IN1 in the actual
C     driver TRDRV1
C
      DO ITASK = 1, NTASK
C
C        Give relative offset for the coefficients
C
         DO INDEX = 1, 4
            DO I = 1, NFSYM
               ICMOQR1(I,INDEX) = ICMOQR(I,INDEX)
     &                          + (NSTRS(I,INDEX,ITASK)-1)*NFBAS(I,0)
            ENDDO
         ENDDO
C
         CALL TRDRV1(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR1,
     &               NSTR1(1,0,1,ITASK),ANTIS,LMP2,WORK(KINDX),KQ,KE,
     &               KIBE,DINTSKP)
C
      ENDDO
C
      CALL MEMREL('TRDR1T',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDR1T')
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDRV1 */
      SUBROUTINE TRDRV1(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,
     &                  NSTR,ANTIS,LMP2,INDX,KQ,KE,KIBE,DINTSKP)
C
C     Written by Luuk Visscher, december 1996
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      LOGICAL ANTIS, LMP2, NOPV, NODV, TRIAN(2)
      DIMENSION NSTR(2,0:2,4),INDX(3,*),KQ(2,4),KE(2,4),KIBE(2,4)
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4)
      DIMENSION IJP12(0:7)
      DIMENSION IJP34(2)
      DIMENSION IHTM(0:7),IHM(0:7),IFM(0:7),NHTM(2)
      DIMENSION WORK(*),DINTSKP(*)
C
      CALL QENTER('TRDRV1')
      KFRSAV = KFREE
C
C     Define the packing of the 2-index transformed integrals
C
CTROND      TRIAN(1) = ISAME(2).EQ.ISAME(1) - can we use this ? TODO
      TRIAN(1) = .FALSE.
      TRIAN(2) = .FALSE.
      CALL PCK2IN(NSTR,TRIAN,IPRINT)
C
C     Make index array for half-transformed integrals
C
      NHTM(1) = 0
      NHTM(2) = 0
      NHTMT = 0
      DO IREPIJ = 1, NFSYM
         DO IREPAB = 0, NBSYM-1
            IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
               IHTM(IREPAB) = NHTM(IREPIJ) + NHTMT
               NHTM(IREPIJ) = NHTM(IREPIJ) + NSPCK(IREPAB,0)*NZ
            ENDIF
         ENDDO
         NHTMT = NHTMT + NHTM(IREPIJ)*NFPCK12(IREPIJ)
      ENDDO
C
C     Allocate memory for half transformed integrals
C
      IF (IPRINT .GE. 0) WRITE (LUPRI,1000) NHTMT*8./(1024.*1024.)
      CALL MEMGET('REAL',KHTMAT,NHTMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KHTMAT),NHTMT)
C
C     Prepare for screening
C
      IF(SCRTRA.GT.D0) THEN
        CALL PR4SCR(KGAB,KDRIJ,WORK,KFREE,LFREE,DINTSKP,
     &              WORK(KQ(1,1)),WORK(KQ(1,2)),
     &              WORK(KQ(1,3)),WORK(KQ(1,4)),
     &              NDMOQR,ICMOQR,NSTR,ISAME,IPRINT)
      ELSE
        KDRIJ = KFREE
        KGAB  = KFREE
      ENDIF
C
C     Loop over (XX|LL) and (XX|SS) integrals
C
C     ----------------
      DO 100 IC = 1, 2
C     ----------------
C
         CALL SHRNGE (IC,INTFLG,I2TYP,IASTRT,IBSTRT,IASMAX,IBSMAX)
         IF (I2TYP.LT.0) THEN
            GOTO 100
         ELSEIF (I2TYP.EQ.1) THEN
            ICS = 1
            ICF = 1
         ELSEIF (I2TYP.EQ.2) THEN
            ICS = 2
            ICF = 2
         ELSE
            ICS = 1
            ICF = 2
         ENDIF
C
C        Loop over the shells
C
         DO ISHLA = IASTRT, IASMAX
         DO ISHLB = IBSTRT, ISHLA
C
            KFRSAV3 = KFREE
            NINSHA  = NINSH(ISHLA,-1,INDX,IPRINT)
            NINSHB  = NINSH(ISHLB,-1,INDX,IPRINT)
C
C           Set up gather array for distributions, 
C           calculate memory requirements for this batch.
C
C           NOTE TODO : The actual dimension of INDXAB is kept in the
C           common block dcbtra.h. When we want to parallelize it may
C           be better to allocate it always as MXINSH*MXINSH*4,
C           or (even better) to send it also down to HERMIT.
C           Luuk.
C
            LDXAB = 5*NINSHA*NINSHB
            CALL MEMGET('INTE',KDXAB,LDXAB,WORK,KFREE,LFREE)
C           Get memory needed for Hermit
            MWHER = MXMEMHER(IC,I2TYP)
            NSIZG = LFREE - MWHER
            NUMDIS = NDISTRN(I2TYP,IPRINT,.TRUE.,.FALSE.,.FALSE.,.TRUE.,
     &                       INDX,WORK(KDXAB),IJP12,IDUM,NSIZG,IDUM)
            NSIZH = 0
            DO IREPAB = 0, NBSYM-1
               IHM(IREPAB) = NSIZH
               IREPIJ = JBTOF(IREPAB,1)
               NSIZH = NSIZH + IJP12(IREPAB)*NFPCK12(IREPIJ)
            ENDDO
C
            IF (IPRINT .GE. 3) THEN
               WRITE (LUPRI,'(//,4(1X,A,2I4,/))')
     &         ' Calculation of integral distribution set:',ISHLA,ISHLB
               WRITE (LUPRI,'(//,4(1X,A,I16,/))')
     &         ' Number of distributions in this set:     ',NUMDIS,
     &         ' Size of scalar integral batch:           ',NSIZG,
     &         ' Size of half-transformed integral batch: ',NSIZH
            END IF
C
            CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KGMAT),NSIZG)
C
C           Get distributions
C
            NODV = NASHT.EQ.0
            NOPV = NASHT.LT.2
            CALL CALDIS(I2TYP,WORK(KGMAT),INDX,WORK(KDXAB),
     &               NODV,NOPV,WORK(KGAB),WORK(KDRIJ),DINTSKP,
     &               SCRTRA,WORK(KFREE),LFREE,IPRINT)
C
C           Do first step of 4-index transformation :
C           Transform first pair of indices
C 
            CALL MEMGET('REAL',KHMAT,NSIZH*NZ,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KHMAT),NSIZH*NZ)
            CALL MS4IN1 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,.FALSE.,
     &                   NSTR(1,0,1),NSTR(1,0,2),
     &                   NDMOQR,ICMOQR,IJP12,
     &                   WORK(KQ(1,1)),WORK(KQ(1,2)),
     &                   WORK(KGMAT),WORK(KHMAT))
C
C           Scatter half-transformed integrals to the right position
C
            CALL MS4INS (WORK,KFREE,LFREE,IPRINT,IC,IHM,NHTM,IHTM,INDX,
     &                   WORK(KDXAB),WORK(KHMAT),WORK(KHTMAT))
C
            CALL MEMREL('TRDRV1',WORK,1,KFRSAV3,KFREE,LFREE)
C
         ENDDO
         ENDDO
C
C----------------
 100     CONTINUE
C----------------
C
C
C     Screening statistics
C
      IF(SCRTRA.GT.D0) THEN
        CALL ST4SCR(DINTSKP)
      ENDIF
      NSIZF = 0
C
      DO IREPAB = 0, NBSYM-1
         IREPIJ = JBTOF(IREPAB,1)
         NSIZF = NSIZF + NFPCK34(IREPIJ)*NZ*NZ*NFPCK12(IREPIJ)
      ENDDO
C
      IF (IPRINT .GE. 0) WRITE (LUPRI,1001) NSIZF*8./(1024.*1024.)
      CALL MEMGET('REAL',KFMAT,NSIZF,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KFMAT),NSIZF)
C
      icount1 = 0
      icount2 = 0
      do ix = 1, nhtmt
c        print*,ix,work(khtmat+ix-1)
         if (abs(work(khtmat+ix-1)).gt.1.0d-16) then
            icount1 = icount1 + 1
         else
            icount2 = icount2 + 1
         endif
      enddo
      print*,' Non-zero half-transformed integrals ',icount1
      print*,' Zero half-transformed integrals     ',icount2
      print*,' Total half-transformed integrals    ',icount1+icount2
      ICS = 1
      ICF = 2
      CALL MS4IN2 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,
     &             NSTR(1,0,3),NSTR(1,0,4),
     &             NDMOQR(1,1,3),ICMOQR(1,3),
     &             NFPCK12,WORK(KQ(1,3)),WORK(KQ(1,4)),
     &             WORK(KHTMAT),WORK(KFMAT))
c     write (LUPRI,'(//A)') ' The fully transformed integrals'
c     do ix = 1, nsizf
c        print*,ix,work(kfmat+ix-1)
c     enddo
C
c     IF (ANTIS) WRITE(LUPRI,'(A)') 'Antisymmetrizing'
c     IF (LMP2) THEN
c        WRITE(LUPRI,'(A)') 'DOING MP2'
C         CALL MP2SEC(WORK(KE(1,1)),WORK(KE(1,2)),NSTR)
c     ENDIF
C
C     Initialize MOLFDIR integral file
C
      CALL MDINTI(IPRINT)
C
C     Number of boson symmetry for each fermion irrep (parity)
C
      NBSYMP = NBSYM/NFSYM
C
C     Complete the integrals to the desired classes
C
      CALL MEMGET('INTE',KDXKR12,NFPCK12T*2,WORK,KFREE,LFREE)
      CALL MKINDXKR (NSTR,.FALSE.,WORK(KDXKR12))
C
      CALL MEMGET('INTE',KDXB12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MKINDXB (NSTR,WORK(KIBE(1,1)),WORK(KIBE(2,1)),
     &              WORK(KIBE(1,2)),WORK(KIBE(2,2)),
     &              .FALSE.,WORK(KDXB12))
      CALL MEMGET('INTE',KDXB34,2*NFPCK34T,WORK,KFREE,LFREE)
      CALL MKINDXB (NSTR,WORK(KIBE(1,3)),WORK(KIBE(2,3)),
     &              WORK(KIBE(1,4)),WORK(KIBE(2,4)),
     &              .FALSE.,WORK(KDXB34))
C
      NIJOFF = 0
      NKLOFF = 0
      IOFF = 0
      DO IREPIJ = 1, NFSYM
         NIJ = NFPCK12(IREPIJ)
         NKL = NFPCK34(IREPIJ)
         NF = NFPCK34(IREPIJ)*NZ*NZ*NBSYMP
         CALL DRV1SYM(WORK,KFREE,LFREE,IPRINT,NSTR(1,0,3),NSTR(1,0,4),
     &                IREPIJ,NIJ,WORK(KDXKR12),NIJOFF,
     &                WORK(KDXB12),WORK(KDXB34),NKLOFF,
     &                WORK(KFMAT+IOFF),NF)
         NIJOFF = NIJOFF + NIJ
         NKLOFF = NKLOFF + NKL
         IOFF = IOFF + NF*NIJ
      ENDDO
C
C     Close integral file
C
      CALL MDINTF
C
C     Release all memory and exit
C
      CALL MEMREL('TRDRV1',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDRV1')
C
      RETURN
C
 1000 FORMAT (/' Storing all half-transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
 1001 FORMAT (/' Storing all fully-transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MS4INS*/
      SUBROUTINE MS4INS (WORK,KFREE,LFREE,IPRINT,IC,IHM,NHTM,IHTM,
     &                   INDX,INDXAB,HMAT,HTMAT)
C
C     Sort half-transformed integrals
C
C     Luuk Visscher
C
C
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbtra.h"
#include "dcbbas.h"
C
      DIMENSION INDX(3,*),INDXAB(NINSHA,NINSHB,4)
      DIMENSION IHTM(0:7),IHM(0:7),NHTM(2)
      DIMENSION HMAT(*),HTMAT(*)
      DIMENSION WORK(*)
      DIMENSION IND(2)
C
#include "ibtfun.h"
C
C     Loop over the functions inside the shell-block.
C     We do it this way rather than just over the boson irreps of the
C     density and than over the pairs within that irrep, because we
C     need to access the information in INDXAB. The alternative gives
C     less jumping through memory, but requires introduction of a new 
C     n^2 type array in which indab is stored.
C
      DO 20 INDAR = 1, NINSHA
         DO 10 INDBR = 1, NINSHB
            INDAB  = INDXAB(INDAR,INDBR,1)
            IF (INDAB.EQ.0) GOTO 10
            IREPAB = INDXAB(INDAR,INDBR,3)
            IPAB   = INDXAB(INDAR,INDBR,4)
            IREPIJ = JBTOF(IREPAB,1)
            IOFF = (IHM(IREPAB)+(IPAB-1)*NFPCK12(IREPIJ))*NZ+1
            CALL IUNPCK(INDAB,2,IND)
            INDA = IND(1)
            INDB = IND(2)
            IREPA = INDX(2,INDA)
            IREPB = INDX(2,INDB)
            INDA1 = INDA - ICOS(IREPA+1,IC)
            INDB1 = INDB - ICOS(IREPB+1,IC)
C
C           The offset within the untransformed indexes. These were
C           the right-hand side (slowest varying) and are sorted to
C           become the left hand side. Symmetrization is also done
C           in this step (the batch indices were triangular).
C           
            JOFFAB = ISPCK(IREPA,IREPB,IC)
     &             + (INDB1-1)*NBBAS(IREPA,IC)+INDA1
            JOFFBA = ISPCK(IREPB,IREPA,IC)
     &             + (INDA1-1)*NBBAS(IREPB,IC)+INDB1
            JOFF1 = IHTM(IREPAB) + JOFFAB
            JOFF2 = IHTM(IREPAB) + JOFFBA
            DO IZ = 1, NZ
C
C              We copy all transformed indices that now become slowest
C              varying (right-hand). This gives a stride.
C
               CALL DCOPY(NFPCK12(IREPIJ),HMAT(IOFF),1,
     &                   HTMAT(JOFF1),NHTM(IREPIJ))
               IF (JOFF1.NE.JOFF2)
     &         CALL DCOPY(NFPCK12(IREPIJ),HMAT(IOFF),1,
     &                   HTMAT(JOFF2),NHTM(IREPIJ))
               IOFF  = IOFF  + NFPCK12(IREPIJ)
               JOFF1 = JOFF1 + NSPCK(IREPAB,0)
               JOFF2 = JOFF2 + NSPCK(IREPAB,0)
            ENDDO
  10     CONTINUE
  20  CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MS4IN2 */
      SUBROUTINE MS4IN2 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,NSTR1,NSTR2,
     & NDMOQR,ICMOQR,IJPAIRS,Q1,Q2,GMAT,HMAT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher January 1997.
C     A twin brother of MS4IN1. Don't get confused by the name of the
C     indices. We DO transform the second index pair.
C
C     PURPOSE : Do 4-index transformation to molecular spinor basis
C               Second index-pair transformation.
C   
C
C     Input :
C
C     - IPRINT       Print flag
C     - ICS          First class of integrals 1 : (XX|LL), 2 : (XX|SS)
C     - ICL          Last class of integrals
C     - NSTR1        Number of active spinors for index 1
C     - NSTR2        Number of active spinors for index 2
C     - NFPCK        Number of spinor pairs for each compound symmetry
C     - IFPCK        Pointers to spinor pairs
C     - NDMOQR       Dimensions of the coefficient array
C     - ICMOQR       Pointers to coefficients
C     - IJPAIRS      Number of blocks for each boson symmetry
C     - Q1           Coefficients for index 1
C     - Q2           Coefficients for index 2
C     - GMAT         Symmetry packed scalar integrals
C
C     Output :
C
C     - HMAT         Symmetry packed half-transformed integrals
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
      INTEGER   NSTR1(2),NSTR2(2)
      DIMENSION WORK(*)
      DIMENSION Q1(*),Q2(*),GMAT(*),HMAT(*)
      DIMENSION IJPAIRS(2)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2)
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbibt.h"
C
      CALL QENTER('MS4IN2')
C
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
C
C     Transform righthand indices to molecular spinor basis
C
      IOFF = 1
      JOFF = 1
      DO IREPIJ = 1, NFSYM
        CALL MEMGET ('REAL',KHTMP,NFPCK34(IREPIJ)*NZ,WORK,KFREE,LFREE)
        CALL DZERO (WORK(KHTMP),NFPCK34(IREPIJ)*NZ)
        DO IJ = 1, IJPAIRS(IREPIJ)
          DO IREPPQ = 0, NBSYM-1
          IF (IREPIJ.EQ.JBTOF(IREPPQ,1)) THEN
            DO IZ = 1, NZ
               DO IC = ICS, ICF
                  IF (ICS.EQ.2) IOFF = IOFF + NSPCK(IREPPQ,1)
                  DO IREPQ = 0, NBSYM-1
                     IREPP = IBTXOR(IREPPQ,IREPQ)
                     IREPI =  JBTOF(IREPP,IC)
                     IREPJ =  JBTOF(IREPQ,IC)
                     KOFF = ICMOQR(IREPI,1) + IBBAS(IREPP,IC) 
     &                    - IBAS(IREPI)
                     LOFF = ICMOQR(IREPJ,2) + IBBAS(IREPQ,IC) 
     &                    - IBAS(IREPJ)
                     NP = NBBAS(IREPP,IC)
                     NQ = NBBAS(IREPQ,IC)
                     NI = NSTR1(IREPI)
                     NJ = NSTR2(IREPJ)
                     NRQ1 = NDMOQR(1,IREPI,1)
                     NCQ1 = NDMOQR(2,IREPI,1)
                     NRQ2 = NDMOQR(1,IREPJ,2)
                     NCQ2 = NDMOQR(2,IREPJ,2)
                     IF (IC.EQ.1) THEN
                        IREPPI = IREPP
                        IREPQI = IREPQ
                     ELSE
                        IREPPI = IBTXOR(IXYZ,IREPP)
                        IREPQI = IBTXOR(IXYZ,IREPQ)
                     ENDIF
                     IF ((NP*NQ.NE.0).AND.(NI*NJ.NE.0)) THEN
                       CALL QTRANS('AOMO','S',D0,NP,NQ,NI,NJ,
     &                      GMAT(IOFF),NP,NQ,1,IPQTOQ(1,0),
     &                      WORK(KHTMP),NI,NJ,NZ,IPQTOQ(1,IREPPQ),
     &                      Q1(KOFF),NRQ1,NCQ1,NZ,IPQTOQ(1,IREPPI),
     &                      Q2(LOFF),NRQ2,NCQ2,NZ,IPQTOQ(1,IREPQI),
     &                      WORK(KFREE),LFREE,IPRINT)
C
C                      Order the integrals with NKL as first and IZ2
C                      as second index
C
                       DO IZ2 = 1, NZ
                          JOFF1 = JOFF + (IZ2-1)*NFPCK34(IREPIJ) +
     &                                   IFPCK34(IREPI,IREPJ)
                          JOFF2 = KHTMP + (IZ2-1)*NI*NJ
                          CALL DAXPY(NI*NJ,D1,WORK(JOFF2),1,
     &                               HMAT(JOFF1),1)
                       ENDDO
                     ENDIF
                     IOFF = IOFF + NP * NQ
                  ENDDO
               ENDDO
               JOFF = JOFF + NFPCK34(IREPIJ)*NZ
            ENDDO
          ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      CALL QEXIT('MS4IN2')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DRV1SYM*/
      SUBROUTINE DRV1SYM(WORK,KFREE,LFREE,IPRINT,NSTR3,NSTR4,
     & IREPIJ,NIJ,INDXKR,NIJOFF,INDXB12,INDXB34,NKLOFF,FMAT,NF)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher May 1999.
C
C     PURPOSE : Do 4-index transformation to molecular spinor basis
C               Symmetrization of integrals (which are in core).
C
C     Input :
C
C     - IPRINT       Print flag
C     - NSTR3        Number of active spinors for index 3
C     - NSTR4        Number of active spinors for index 4
C     - IREPIJ       Fermion symmetry (parity)
C     - NIJ          Number of blocks for this fermion symmetry
C     - Q3           Coefficients for index 3
C     - Q4           Coefficients for index 4
C
C     Output is written directly in the MOLFDIR-type file MDCINT
C
#include "implicit.h"
#include "priunit.h"
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION WORK(*)
      DIMENSION INDXKR(2,*),INDXB12(2,*),INDXB34(2,*)
      DIMENSION FMAT(NF,*)
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbtr3.h"
C
      NKL = NFPCK34(IREPIJ)
      DO IJ = 1, NIJ
         IKR = INDXKR(1,IJ+NIJOFF)
         JKR = INDXKR(2,IJ+NIJOFF)
         CALL SYMFINT(IPRINT,IREPIJ,
     &                INDXB12(1,IJ+NIJOFF),INDXB34(1,1+NKLOFF),
     &                IJ,IKR,JKR,NSTR3,NSTR4,NKL,FMAT(1,IJ))
      ENDDO
C
      RETURN
      END
