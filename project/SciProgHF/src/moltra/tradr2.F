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

! define task symbols for CALL DIRAC_PARCTL( task )
#include "dirac_partask.h"
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDR2T */
      SUBROUTINE TRDR2T(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOBF,ICMOIJ,
     &                  NSTR,ANTIS,LMP2,EMP2,KINDX,KQ,KE,KIBE,DINTSKP)
C
C     Written by Luuk Visscher, Februari 1997
C     Modifications: Jon K. Laerdahl, June 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
#include "maxorb.h"
#include "aovec.h"
C
#if defined (VAR_MPI)
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#endif
#include "infpar.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "blocks.h"
#include "dcbtra.h"
#include "dcbtr3.h"
      LOGICAL ANTIS,LMP2
      DIMENSION NSTRT(4),KQ(2,4),KE(2,4),KIBE(2,4)
      DIMENSION NSTR(2,0:2,4),NDMOBF(2,2,4)
      DIMENSION NSTTIJ(4),NSTRIJ(2,0:2,4),ICMOIJ(2,4)
      DIMENSION WORK(*),DINTSKP(*)
C
      CALL QENTER('TRDR2T')
      KFRSAV = KFREE
C
C     TODO (/hjaaj aug 2001):
      IF (.NOT. LMP2)
     &     CALL QUIT('TRDR2T: Only MP2 is implemented for strategy 2')
C
C     Check that no negative energy solutions are active.
C     Active negative energy solutions is not implemented
C     and is NOT a good idea unless one knows exactly
C     what one is doing...
C                                     Jon
C
      DO IND = 1,4
         NSTRT(IND) = 0
         DO IFRP = 1,NFSYM
            IF (NSTR(IFRP,2,IND).NE.0) 
C     TODO (/hjaaj aug 2001): remember to modify mp2.tex 
C          (.ACTIVE, .POSITR) if this is implemented
     &           CALL QUIT
     &           ('TRDR2T: Negative energy solutions in active space!')
            NSTRT(IND) = NSTRT(IND) + NSTR(IFRP,0,IND)
         ENDDO
      ENDDO
      IF (NSTRT(1).NE.NSTRT(3)) CALL QUIT
     &     ('TRD2T: Different number of active for index 1 and 3!')
      IF (NSTRT(2).NE.NSTRT(4)) CALL QUIT
     &     ('TRD2T: Different number of active for index 2 and 4!')
C
C     Initialize
C
      CALL ICOPY(16,NDMOBF,1,NDMOQR,1)
      CALL ICOPY( 8,ICMOIJ,1,ICMOQR,1)
C
C     Zero MP2 energy and start looping over tasks
C
      IF (LMP2) THEN
         WRITE(LUPRI,'(/,15X,A,/)') '**** Calculating MP2 energy. ****'
         EMP2 = D0
      ENDIF
C
C     Get Cauchy-Schwartz integrals if screening
C
      IF(SCRTRA.GT.D0) THEN
         N2GAB   = NSYMBL*NSYMBL
         CALL MEMGET('REAL',KGAB,N2GAB,WORK,KFREE,LFREE)
         IJOB   = 0
         ITYPE  = 0
         IGTYP  = 1
         MAXDF = 0
         CALL GETGAB(IJOB,ITYPE,IGTYP,MAXDF,
     &        IPRINT,WORK(KGAB),WORK(KFREE),LFREE)
         CALL DZERO(DINTSKP,8)
      ELSE
         KGAB = KFREE
      ENDIF
C
C     *****************************************
C     ****** P A R A L L E L    C O D E  ******
C     *****************************************
C
#if defined (VAR_MPI)
      IF (PARCAL .AND. (.NOT. MP2ORG)) THEN
C
C     *** New parallel MP2 algorithm ***
C
         KDRIJ = KFREE
         CALL MEMGET('INTE',KDSLST,MAXSHL*6,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KTIMDS,MAXSHL,WORK,KFREE,LFREE)
         CALL TR2MST(EMP2,WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,
     &        KINDX,KGAB,KDRIJ,INTFLG,ICMOIJ,NSTR,NSTRT,NSTRIJ,NSTTIJ,
     &        ANTIS,LMP2,MAXSHL,WORK(KDSLST),WORK(KTIMDS),DINTSKP)
      ELSE IF (PARCAL .AND. MP2ORG) THEN
C
C     *** Original parallel MP2 algorithm ***
C
C     Get hold of the slaves
C     ( ITASK = 2 for parallel integral transformation )
C
         CALL DIRAC_PARCTL( MOLTRA_PAR )
C     
C     Initialize
C     
         NTIM = 7*NODES
         CALL MEMGET ('REAL',KTIM,NTIM,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KTIM),NTIM)
         KDRIJ = KFREE
         CALL TRAPARI(WORK,KFREE,LFREE,KQ,KE,KIBE,
     &        KINDX,KGAB,KDRIJ,INTFLG,NSTR,ANTIS,LMP2)
C     
C     Divide index transformation in tasks with different batches
C     of spinors. First find number of tasks.
C     
         CALL TSKCHK(MAXTSK,LENTSK,LFREE,NODES,NSTRT(1),IPRINT)
C     
C     Start loop over tasks
C     
         DO IJTSK = 1, MAXTSK
#if defined (MPE)
            CALL MPE_LOG_EVENT(3, 0, "start wait")
#endif
            call interface_mpi_RECV(ITEST,1,df_MPI_ANY_SOURCE,20,
     &           global_communicator,ISTAT)
#if defined (MPE)
            CALL MPE_LOG_EVENT(4, 0, "end wait")
#endif
C     .. ITEST .eq. 0: code for "send me new task"
            IF (ITEST.NE.0) CALL QUIT(' Illegal ITEST in MPI_RECV')
            NODE = ISTAT(df_MPI_SOURCE)
            IF (IPRINT.GE.1) WRITE(LUPRI,'(3(A,I6))') 
     &           '* Sending to node:',NODE,
     &           '; IJTSK  =',IJTSK,', LENTSK =',LENTSK
#if defined (MPE)
            CALL MPE_LOG_EVENT(5, 0, "start sndtsk")
#endif
            call interface_mpi_SEND(IJTSK,1,NODE,30,
     &           global_communicator)
            call interface_mpi_SEND(LENTSK,1,NODE,30,
     &           global_communicator)
#if defined (MPE)
            CALL MPE_LOG_EVENT(6, 0, "end sndtsk")
#endif
         ENDDO
C
C     Get the slaves out of the transformation loop
C     and receive final results
C
         WALLBEF = 0.0
         CALL TRARES(WORK,KFREE,LFREE,DINTSKP,LMP2,EMP2,
     &        WORK(KTIM),WALLBEF)
C
C     Timing information
C
         CALL TR2TIM(WORK(KTIM),NODES,LMP2)
C     
C     Release the slaves...
C     (ITASK = -1)
C     
         CALL DIRAC_PARCTL( EXIT_NODEMENU )
#if defined (VAR_PFS)
C
C       close global file (opened in TRAPARI)
C
        CALL GLOBAL_FILE_CLOSE_READ
#endif
      ELSE
#endif
C
C     *************************************
C     ****** S E R I A L    C O D E  ******
C     *************************************
C     
C       Initialize timing
C
        CALL MEMGET('REAL',KTIM,7,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KTIM),7)
C
C       Divide index transformation in tasks with different batches
C       of spinors. First find number of tasks.
C
        CALL TSKCHK(MAXTSK,LENTSK,LFREE,1,NSTRT(1),IPRINT)
C
C       Start loop over tasks
C
        DO IJTSK = 1, MAXTSK
           CALL UNPKIJ(IJTSK,ISPNR,JSPNR)
           CALL TSKARR(ISPNR,JSPNR,LENTSK,NSTRT,NSTR,NSTRIJ,
     &          ICMOQR,ICMOIJ,IPRINT) 
C
C     *** New MP2 algorithm ***
C
           IF (.NOT. MP2ORG) THEN
              IF (IPRINT.GE.1) WRITE(LUPRI,'(/A/)') 
     &             '**** New MP2 algorithm selected'
              CALL MEMGET('INTE',KDSLST,MAXSHL*6,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KTIMDS,MAXSHL,WORK,KFREE,LFREE)
              CALL TD23V2(IJTSK,MAXTSK,WORK,KFREE,LFREE,IPRINT,EMP2,
     &             ISPNR,JSPNR,LENTSK,INTFLG,NDMOQR,ICMOIJ,NSTRIJ,
     &             WORK(KINDX),KQ,KE,DINTSKP,WORK(KGAB),WORK(KTIM),
     &             MAXSHL,WORK(KDSLST),WORK(KTIMDS))
C
C     *** Original MP2 algorithm ***
C
           ELSE
              IF (IPRINT.GE.1) WRITE(LUPRI,'(/A/)') 
     &             '**** Original MP2 algorithm selected'
              CALL TRDR23(WORK,KFREE,LFREE,IPRINT,EMP2,ISPNR,JSPNR,
     &             LENTSK,INTFLG,NDMOQR,ICMOIJ,NSTRIJ,
     &             ANTIS,LMP2,WORK(KINDX),KQ,KE,
     &             DINTSKP,WORK(KGAB),WORK(KTIM))
           END IF
        ENDDO
C
C       Timing information
C
        CALL TR2TIM(WORK(KTIM),1,LMP2)
#if defined (VAR_MPI)
      ENDIF
#endif
C
C     Screening statistics
C
      IF(SCRTRA.GT.D0) THEN
        CALL ST4SCR(DINTSKP)
      ENDIF
C
      CALL MEMREL('TRDR2T',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDR2T')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TR2TIM */
      SUBROUTINE TR2TIM(TIM1,NDIM,LMP2)
C***********************************************************************
C
C     Timing information from scheme 2.
C
C     Written by T.Saue Sep 18 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
      LOGICAL LMP2
      DIMENSION TIM1(7,NDIM)
C
      WRITE(LUPRI,*)
      IF(PARCAL) THEN
        DO I = 1,NDIM
          WRITE(LUPRI,'(A,I4)') '** Node :',I
          CALL TIMTXT('Time used in CALDIS is: ',TIM1(1,I),LUPRI)
          CALL TIMTXT('Time used in MS4IN1 is: ',TIM1(2,I),LUPRI)
          CALL TIMTXT('Time used in TQTRNS is: ',TIM1(3,I),LUPRI)
          CALL TIMTXT('Time used in LQTRNS is: ',TIM1(4,I),LUPRI)
          IF(LMP2) THEN
            CALL TIMTXT('Time used in MP2CAL is: ',TIM1(5,I),LUPRI)
          ENDIF
        ENDDO
      ELSE
          CALL TIMTXT('Time used in CALDIS is: ',TIM1(1,1),LUPRI)
          CALL TIMTXT('Time used in MS4IN1 is: ',TIM1(2,1),LUPRI)
          CALL TIMTXT('Time used in TQTRNS is: ',TIM1(3,1),LUPRI)
          CALL TIMTXT('Time used in LQTRNS is: ',TIM1(4,1),LUPRI)
          IF(LMP2) THEN
            CALL TIMTXT('Time used in MP2CAL is: ',TIM1(5,1),LUPRI)
          ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDR23 */
      SUBROUTINE TRDR23(WORK,KFREE,LFREE,IPRINT,
     &     EMP2,ISPNR,JSPNR,LENTSK,INTFLG,NDMOQR,ICMOQR,
     &     NSTR,ANTIS,LMP2,INDX,KQ,KE,DINTSKP,GABRAO,TIM1)
C
C     Written by Luuk Visscher, December 1996
C     Modifications by Jon Laerdahl, June 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbmpt.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
      LOGICAL ANTIS,LMP2,NODV,NOPV,L2TYP(0:3),IBEQJB,TRIAN(2)
      DIMENSION KQ(2,4),KE(2,2),NSTR(2,0:2,4),INDX(3,*),
     &          NSTRT(4),ITQMT(0:7,0:7,2),
     &          NDMOQR(2,2,4),ICMOQR(2,4),
     &          IJP12(0:7),IJP34(2),IHM(0:7),
     &          TIM1(7),DINTSKP(*),GABRAO(*),WORK(*)
      DIMENSION ISAMEIJ(4)
C
#include "ibtfun.h"
C
C     Notice: Range for index 1, 2, 3 and 4 are not
C     the same. E1(3) is NOT necessarily = E3(3)
C     and in MP2 E1(3).NE.E2(3).
C
C     To implement MSOUT: Get absolute position of
C     active occupied Kramers pairs from outside.
C     Apart from this Luuk has already the output
C     routine for MOLFIR integrals ready. Use this..
C
      CALL QENTER('TRDR23')
      KFRSAV = KFREE
      DO I = 1,4
        NSTRT(I) = 0
        DO IFRP = 1, NFSYM
          NSTRT(I) = NSTRT(I) + NSTR(IFRP,0,I)
        ENDDO
      ENDDO
C
C     Define the packing of the 2-index transformed integrals
C
      TRIAN(1) = ISAME(2).EQ.ISAME(1)
      TRIAN(2) = .FALSE.
      CALL PCK2IN(NSTR,TRIAN,IPRINT)
C
C     Define the packing of the 3-index transformed integrals
C
      CALL PCK3IN(ITQMT,NTQMT,NSTR,IPRINT)
C
      IF (IPRINT .GE. 1) WRITE (LUPRI,1000) NTQMT*8./(1024.*1024.)
C
      IF (LFREE.LT.NTQMT) THEN
            WRITE(LUPRI,'(/A,/A,/A)') 
     &       '### Not enough memory for 3-index transf. matrix.   ###',
     &       '### Reduce memory requirements by setting           ###',
     &       '### .IJTSK to a lower value!                        ###'
C
      ENDIF
      CALL MEMGET('REAL',KTQMAT,NTQMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KTQMAT),NTQMT)
C
C     Loop over (XX|LL) and (XX|SS) integrals
C
      IF (TESTLS) WRITE(LUPRI,'(/A)') 'YES: Doing (LL|SS) ONLY.'
      IF (TESTSL) WRITE(LUPRI,'(/A)') 'YES: Doing (SS|LL) ONLY.'
C
C     Get screening matrices
C
      IF(SCRTRA.GT.0) THEN
CTROND  Not using sameness relations...
CTROND  ...this may be improved TODO
        DO I = 1,4
          ISAMEIJ(I) = I
        ENDDO
        CALL PR4SC2(KDRIJ,WORK,KFREE,LFREE,GABRAO,
     &              WORK(KQ(1,1)),WORK(KQ(1,2)),
     &              WORK(KQ(1,3)),WORK(KQ(1,4)),
     &              NDMOQR,ICMOQR,NSTR,ISAMEIJ,IPRINT)
      ELSE 
        KDRIJ = KFREE
      ENDIF
C
C     ----------------
      DO 100 IC = 1, 2 
C     ----------------
C
         IF ((IC.EQ.1).AND.TESTLS) GOTO 100
         IF ((IC.EQ.2).AND.TESTSL) GOTO 100
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
            NINSHA = NINSH(ISHLA,-1,INDX,IPRINT)
            NINSHB = NINSH(ISHLB,-1,INDX,IPRINT)
C
C           Set up gather array for distributions, 
C           calculate memory requirements for this batch.
C
C           NOTE : The actual dimension of INDXAB is kept in the
C           common block dcbtra.h. When we want to parallelize it may
C           be better to allocate it always as MXINSH*MXINSH*4,
C           or (even better) to send it also down to HERMIT.
C           Luuk. TODO
C
C           First a dummy run to calculate the memory requirements
C           for GMAT. Max size of GMAT is set in MAXSCL from
C           .SCLMEM in input or chosen based on MWHER.
C
C     
            IF(MAXSCL.EQ.0) THEN
               MWHER = MXMEMHER(IC,I2TYP)
               MEMSV1 = MIN(LFREE/10, 1 000 000) + MWHER
               NSIZG = LFREE - MEMSV1
C              subtract max mem needed in Hermit and 10%/1Mw for anything else
            ELSE
               NSIZG = MAXSCL
               IF (MAXSCL.GT.LFREE/100*80) WRITE(LUPRI,'(2(1X,A,I10))')
     &              '.GT. 80%.., LFREE:',LFREE,'MAXSCL:',MAXSCL
               IF (MAXSCL.GT.LFREE/100*95) THEN
                  WRITE(LUPRI,'(2(1X,A,I10))')
     &                 '.GT. 95%.., LFREE:',LFREE,'MAXSCL:',MAXSCL
                  CALL QUIT('GMAT too large...')
               ENDIF
            ENDIF
C
            NUMDIS = NDISTRN(I2TYP,IPRINT,.FALSE.,.FALSE.,.FALSE.,
     &           .TRUE.,INDX,IDUM,IJP12,IDUM,NSIZG,IDUM)
C
C     If several passes are needed, memory requirements
C     must be reconsidered
C
            IF (NPASS.GT.1) THEN
               NGBTMP = (IRAT*NSIZG)/((1+IRAT)*NPASS)
               MEMSV2 = 2*NGBTMP
               NSIZG = LFREE - MAX(MEMSV1,MEMSV2)
               NUMDIS = NDISTRN(I2TYP,IPRINT,.FALSE.,.FALSE.,.FALSE.,
     &              .TRUE.,INDX,IDUM,IJP12,IDUM,NSIZG,IDUM)
            END IF
C
            LDXAB = 5*NINSHA*NINSHB
            CALL MEMGET('INTE',KDXAB,LDXAB,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KIJPASS,8*NPASS,WORK,KFREE,LFREE)
C
C           Set INDXAB and INDPASS
C
            NUMDIS = NDISTRN(I2TYP,IPRINT,.TRUE.,.TRUE.,.FALSE.,.TRUE.,
     &              INDX,WORK(KDXAB),IJP12,IDUM,NSIZG,WORK(KIJPASS))
C
            NSIZH = 0
            DO IREPAB = 0, NBSYM-1
               IHM(IREPAB) = NSIZH
               IREPIJ      = JBTOF(IREPAB,1)
               NSIZH       = NSIZH + IJP12(IREPAB)*NFPCK12(IREPIJ)
            ENDDO
C
            CALL MEMGET('REAL',KHMAT,NSIZH*NZ,WORK,KFREE,LFREE)
            CALL DZERO (WORK(KHMAT),NSIZH*NZ)
C
            IF (IPRINT .GE. 3) THEN
               WRITE (LUPRI,
     &              '(//,1X,A,2I16/,2(1X,A,I16/),2(1X,A,I16,F10.3/))')
     &         ' Calculation of integral distribution set:',
     &              ISHLA,ISHLB,
     &         ' Number of distributions in this set:     ',NUMDIS,
     &         ' Number of passes:                        ',NPASS,
     &         ' Size of scalar integral batch:           ',NSIZG,
     &              NSIZG*8./(1024.*1024.),
     &         ' Size of half-transformed integral batch: ',NSIZH,
     &              NSIZH*8./(1024.*1024.)
            END IF
C
            IF (NPASS.GT.1) THEN
C
C              Initialize scalar integral buffer file and arrays
C              Use the space reserved for GMAT now as buffer storage
C
               LGFIL = LUTRA2
               NGBFSZ = (IRAT*NSIZG)/((1+IRAT)*NPASS)
C
C              Allocate the buffer arrays in WORK and set KGMAT to the
C              start of this allocation section. We now only need to
C              send down GMAT into HERMIT to find this section.
C
               CALL MEMGET('INTE',KALGREC,NPASS,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KALGBUF,NPASS,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KAIGBUF,NPASS*NGBFSZ,WORK,KFREE,LFREE)
               CALL MEMGET('REAL',KARGBUF,NPASS*NGBFSZ,WORK,KFREE,LFREE)
C
C              As we refer to the places inside the section we need
C              to subtract KGMAT.
C
               KGMAT  = KALGREC
               KLGREC = KALGREC - KGMAT + 1
               KLGBUF = KALGBUF - KGMAT + 1
               KIGBUF = KAIGBUF - KGMAT + 1
               KRGBUF = KARGBUF - KGMAT + 1
            ELSE
C
C              Write directly into the array GMAT
C
               CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
               CALL DZERO(WORK(KGMAT),NSIZG)
               KLGBUF = 1
            ENDIF
C
C           Get distributions
C
            NODV = NASHT.EQ.0
            NOPV = NASHT.LT.2
C
            CALL TIMEC('START ',TIMSTR,TIMEND)
            CALL CALDIS(I2TYP,WORK(KGMAT),INDX,WORK(KDXAB),
     &                  NODV,NOPV,GABRAO,WORK(KDRIJ),DINTSKP,SCRTRA,
     &                  WORK(KFREE),LFREE,IPRINT)
C
            CALL TIMEC('CALDIS',TIMSTR,TIMEND) 
            TIM1(1) = TIM1(1) + TIMEND - TIMSTR
C
C           Do the transformation of the first index pair and sort the
C           half-transformed integrals to memory
C
C           De-allocate buffers and allocate GMAT if necessary.
C
            IF (NPASS.GT.1) THEN
               CALL MEMREL('TRDRV3.1',WORK,1,KALGBUF,KFREE,LFREE)
               CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
               CALL DZERO(WORK(KGMAT),NSIZG)
            ENDIF
C
            CALL TIMEC('START ',TIMSTR,TIMEND)
            CALL MS4IN1E(WORK,KFREE,LFREE,IPRINT,WORK(KALGREC),ICS,ICF,
     &                   NSTR(1,0,1),NSTR(1,0,2),
     &                   NDMOQR,ICMOQR,WORK(KIJPASS),
     &                   WORK(KQ(1,1)),WORK(KQ(1,2)),
     &                   WORK(KGMAT),NSIZG,WORK(KHMAT)) 
            CALL TIMEC('MS4IN1',TIMSTR,TIMEND) 
            TIM1(2) = TIM1(2) + TIMEND - TIMSTR
C
            CALL MEMREL('TRDRV3.2',WORK,1,KGMAT,KFREE,LFREE)
C
C           Delete scalar integral buffer file
C
            IF (NPASS.GT.1) THEN
               CALL DELGBUF(LGFIL)
            ENDIF
C
C
C     3/4 transformation
C
            CALL TIMEC('START ',TIMSTR,TIMEND)
            CALL TQTRNS (WORK,KFREE,LFREE,IPRINT,IC,IHM,
     &           INDX,WORK(KDXAB),WORK(KHMAT),ITQMT,
     &           ICMOQR(1,3),NDMOQR(1,1,3),WORK(KTQMAT),NSTR(1,0,3),
     &           WORK(KQ(1,3)))
C
            CALL TIMEC('TQTRNS',TIMSTR,TIMEND) 
            TIM1(3) = TIM1(3) + TIMEND - TIMSTR
            CALL MEMREL('TRDR23',WORK,1,KFRSAV3,KFREE,LFREE)
C
         ENDDO
         ENDDO
C
C----------------
 100     CONTINUE
C----------------
C
c      write (LUPRI,'(//A)') ' The 3/4 transformed integrals:'
c      do ix = 1, NTQMT
c         print*,ix,work(KTQMAT+ix-1)
c      enddo
C
      N4QMT = 0
      DO IREPIJ= 1,NFSYM
         DO IREPRS = 0, NBSYM-1
            IF(IREPIJ.EQ.JBTOF(IREPRS,1)) THEN
               NIJ = NFPCK12(IREPIJ) ! Number of ij pairs 
               NKL = NFPCK34(IREPIJ) ! Number of kl pairs 
               N4QMT = N4QMT + NIJ*NKL*NZ**2
            ENDIF
         ENDDO
      ENDDO
C
C
C     Allocate memory for 4/4 Transformed integrals
C
      IF (IPRINT .GE. 1) WRITE (LUPRI,1001) N4QMT*8./(1024.*1024.)
      CALL MEMGET('REAL',K4QMAT,N4QMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(K4QMAT),N4QMT)
C     
C     4/4 transformation
C
      ICS2 = 1
      ICF2 = 2
      L2TYP(1) = MOD(INTFLG,2).EQ.1
      L2TYP(2) = MOD((INTFLG/2),2).EQ.1
      L2TYP(3) = MOD(((INTFLG/2)/2),2).EQ.1
      IF ((.NOT.L2TYP(1)).AND.(.NOT.L2TYP(2))) ICS2 = 2
      IF ((.NOT.L2TYP(2)).AND.(.NOT.L2TYP(3))) ICF2 = 1
      IF (TESTLS) ICS2 = 2
      IF (TESTSL) ICF2 = 1
C
      CALL TIMEC('START ',TIMSTR,TIMEND)
      CALL LQTRNS (WORK,KFREE,LFREE,IPRINT,
     &     ITQMT,ICMOQR(1,4),NDMOQR(1,1,4),
     &     WORK(KTQMAT),NSTR(1,0,3),NSTR(1,0,4),
     &     WORK(KQ(1,4)),WORK(K4QMAT),ICS2,ICF2)
      CALL TIMEC('LQTRNS',TIMSTR,TIMEND) 
      TIM1(4) = TIM1(4) + TIMEND - TIMSTR
C     
c      write (LUPRI,'(//A)') ' The fully transformed integrals'
c      do ix = 1, n4qmt
c         print*,ix,work(k4qmat+ix-1)
c      enddo
C
      CALL TIMEC('START ',TIMSTR,TIMEND)
C
      IF (NZ.EQ.1) THEN
         NCLASS = 4
      ELSEIF (NZ.EQ.2) THEN
         NCLASS = 8
      ELSE
         NCLASS = 16
      ENDIF
C
      CALL DEFCL (NCLASS)
C
C     Allocate memory for naive antisymmetrize
C     -> Fully transformed integrals are
C     few compared to other arrays...
C     Storing final integrals in
C     TQMAT if there is enough space.
C
      LSIZ=16
      IF (NZ.LT.4) LSIZ=8
      DO I=1,4
         LSIZ=LSIZ*NSTRT(I)
      ENDDO
      IF(NTQMT.GE.LSIZ) THEN
         KANT = KTQMAT
         IF(IPRINT.GE.1)
     &        WRITE(LUPRI,'(A)') 
     &        '## Storing final integs. in 3/4-matrix. ##'
      ELSE
         IF(IPRINT.GE.1)
     &        WRITE(LUPRI,'(A)') 
     &        '## Final integs. .GT. 3/4-matrix. ##'
         CALL MEMGET('REAL',KANT,LSIZ,WORK,KFREE,LFREE)
      ENDIF
      IF (IPRINT .GE. 1) WRITE (LUPRI,1002) LSIZ*8./(1024.*1024.)
      CALL DZERO(WORK(KANT),LSIZ)
C
C     FQMAT stored as (ij,iz1,kl,iz2,Pk,Grs,Pkl) and the order
C     of ij (and kl) is gg,uu,ug,gu or as (i,j,iz1,k,l,iz2,Pk,Grs,Pkl)
C
C     F.ex. for Pkl=g, iz2 loops faster than loop over gg/uu for kl,
C     but slower than kl within gg/uu
C
      IOFF3 = 0
      DO IREPKL = 1,NFSYM
         CALL MEMGET('REAL',KFTMAT,2*NFPCK12(IREPKL),WORK,KFREE,LFREE)
         LKR = 0
         DO IREPL = 1, NFSYM
            IREPK = MOD(IREPL+IREPKL,2) + 1
C     Offset for k and l symmetry (i.e. for Pkl=g, offset for uu)
            IOFF2 = IFPCK34(IREPK,IREPL)*NFPCK12(IREPKL)*NZ*NZ
            DO L = 1, NSTR(IREPL,0,4)
               LKR = LKR + 1
               DO K = 1, NSTR(IREPK,0,3)
                  IF (IREPK.EQ.1) THEN
                     KKR = K
                  ELSE
                     KKR = NSTR(1,0,3) + K
                  ENDIF
C     
C     Real and complex classes
C     
C     Offset for K and L
                  IOFF1 = NFPCK12(IREPKL)*NZ*
     &                 (NSTR(IREPK,0,3)*(L-1)+(K-1))
                  DO JCLASS = 1, 4
                     CALL DZERO(WORK(KFTMAT),2*NFPCK12(IREPKL))
                     ICNT1 = 0
                     DO IREPRS = 0, NBSYM-1
                        IF (IREPKL.EQ.JBTOF(IREPRS,1)) THEN
C     Counter for size of IREPKL=1 contribution and offset for
C     IREPKL
                           IOFF6 = ICNT1
                           ICNT1=ICNT1
     &                          +NFPCK12(IREPKL)*NFPCK34(IREPKL)*NZ*NZ
C     Offset for Grs
                           IRPCN1 = 0
                           DO IZ2 = 1, NZ
C     Offset for iz2 
                              IOFF4 = NFPCK12(IREPKL)*NZ*IRPCN1*
     &                             NSTR(IREPK,0,3)*NSTR(IREPL,0,4)
C
                              IRPCN1 = IRPCN1 + 1
                              IQ2 = IPQTOQ(IZ2,IREPRS)
                              IRPCN2 = 0
                              DO IZ1 = 1, NZ
C     Offest for iz1
                                 IOFF5 = NFPCK12(IREPKL)*IRPCN2
                                 IRPCN2 = IRPCN2 + 1
                                 IQ1 = IPQTOQ(IZ1,IREPRS)
                                 IFAC=INTCL(JCLASS,IQ1,IQ2)
                                 IOFF=IOFF1+IOFF2+IOFF3
                                 IOFF=IOFF+IOFF4+IOFF5+IOFF6
C
                                 IF (IFAC.NE.0) THEN
                                    FACTR=IFAC
                                    CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                   WORK(K4QMAT+IOFF),1,
     &                                   WORK(KFTMAT),1)
                                 ENDIF
C 
C     Non-real groups
C
                                 IF (NZ.GT.1) THEN
                                    IFAC=INTCL(JCLASS+4,IQ1,IQ2)
                                    IF (IFAC.NE.0) THEN
                                       FACTR=IFAC
                                       CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                      WORK(K4QMAT+IOFF),1,
     &                                      WORK(KFTMAT+
     &                                      NFPCK12(IREPKL)),1)
                                    ENDIF
                                 ENDIF
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDDO
                     CALL MP2FI1(JCLASS,IREPKL,NSTR,KKR,LKR,
     &                    WORK(KFTMAT),WORK(KANT),
     &                    NFPCK12(IREPKL),NSTRT(1),NSTRT(2),
     &                    NSTRT(3),NSTRT(4),LSIZ)
                  ENDDO
C     
C     Quaternion groups
C     
                  IF(NZ.EQ.4) THEN
                     DO JCLASS = 9, 12
                        CALL DZERO(WORK(KFTMAT),2*NFPCK12(IREPKL))
                        ICNT1 = 0
                        DO IREPRS = 0, NBSYM-1
                           IF (IREPKL.EQ.JBTOF(IREPRS,1)) THEN
C     Counter for size of IREPKL=1 contribution and offset for
C     IREPKL
                              IOFF6 = ICNT1
                              ICNT1=ICNT1
     &                             +NFPCK12(IREPKL)
     &                             *NFPCK34(IREPKL)*NZ*NZ
C     Offset for Grs
                              IRPCN1 = 0
                              DO IZ2 = 1, NZ
C     Offset for iz2 
                                 IOFF4 = NFPCK12(IREPKL)*NZ*IRPCN1*
     &                                NSTR(IREPK,0,3)*NSTR(IREPL,0,4)
C     
                                 IRPCN1 = IRPCN1 + 1
                                 IQ2 = IPQTOQ(IZ2,IREPRS)
                                 IRPCN2 = 0
                                 DO IZ1 = 1, NZ
C     Offest for iz1
                                    IOFF5 = NFPCK12(IREPKL)*IRPCN2
                                    IRPCN2 = IRPCN2 + 1
                                    IQ1 = IPQTOQ(IZ1,IREPRS)
                                    IFAC=INTCL(JCLASS,IQ1,IQ2)
                                    IOFF=IOFF1+IOFF2+IOFF3
                                    IOFF=IOFF+IOFF4+IOFF5+IOFF6
C
                                    IF (IFAC.NE.0) THEN
                                       FACTR=IFAC
                                       CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                      WORK(K4QMAT+IOFF),1,
     &                                      WORK(KFTMAT),1)
                                    ENDIF
                                    IFAC=INTCL(JCLASS+4,IQ1,IQ2)
                                    IF (IFAC.NE.0) THEN
                                       FACTR=IFAC
                                       CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                      WORK(K4QMAT+IOFF),1,
     &                                      WORK(KFTMAT
     &                                      +NFPCK12(IREPKL)),1)
                                    ENDIF
                                 ENDDO
                              ENDDO
                           ENDIF
                        ENDDO
                        CALL MP2FI1(JCLASS,IREPKL,NSTR,KKR,LKR,
     &                       WORK(KFTMAT),WORK(KANT),
     &                       NFPCK12(IREPKL),NSTRT(1),NSTRT(2),
     &                       NSTRT(3),NSTRT(4),LSIZ)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CALL MEMREL('TRDR23',WORK,1,KFTMAT,KFREE,LFREE)
         IOFF3=ICNT1
      ENDDO
C
      I1 = KE(1,1) + (ISPNR-1)*LENTSK
      I2 = KE(1,2)
      I3 = KE(1,1) + (JSPNR-1)*LENTSK
      I4 = KE(1,2)
      CALL CONMP2(EMP2,WORK(KANT),
     &            WORK(I1),WORK(I2),WORK(I3),WORK(I4),
     &            NSTRT(1),NSTRT(2),NSTRT(3),NSTRT(4),
     &            NSTR,ISPNR,JSPNR,IPRINT)
      CALL TIMEC('CONMP2',TIMSTR,TIMEND) 
      TIM1(5) = TIM1(5) + TIMEND - TIMSTR
C     
C     Release all memory and exit
C
      CALL MEMREL('TRDR23',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('TRDR23')
C
      RETURN
C
 1000 FORMAT (/' Storing all 3-index transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
 1001 FORMAT (/' Storing all fully-transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
 1002 FORMAT (/' Storing all final integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TD23V2 */
      SUBROUTINE TD23V2(IJBTCH,MXBTCH,WORK,KFREE,LFREE,IPRINT,ERGMP2,
     &     ISPNR,JSPNR,LENTSK,INTFLG,NDMOQR,ICMOQR,NSTR,INDX,
     &     KQ,KE,DINTSKP,GABRAO,TIM1,NDSLST,IDSLST,TIMDIS)
C
C     Sequential MP2 driver (modified version of TRDR23)
C
C     Vebjorn Bakken, Summer 2003
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dcbtra.h"
#include "dcborb.h"
#include "dcbmpt.h"
#include "dcbmp2.h"
#include "dgroup.h"
      DIMENSION KQ(2,4), KE(2,2), NSTR(2,0:2,4), INDX(3,*),
     &     NSTRT(4), ITQMT(0:7,0:7,2), IJP12(0:7), NDMOQR(2,2,4),
     &     ICMOQR(2,4), TIM1(7), DINTSKP(*), GABRAO(*), WORK(*),
     &     IDSLST(6,NDSLST), TIMDIS(NDSLST), IDSNFO(6)
      DIMENSION ISAMEIJ(4)
      CHARACTER*(24) TMSTMP
      LOGICAL NODV, NOPV, L2TYP(0:3), IBEQJB, TRIAN(2)
C
      CALL QENTER('TD23V2')
      CALL TIMEC('START ',TIMSHS,TIMSHE)
C
C     Initialize
C        Define packing of 2-index and 3-index transformed integrals
C        and determine memory for 3/4 and 4/4 transformed matrices
C
      KFRSAV = KFREE
      CALL MKNSTT(NSTR,NSTRT)
      CALL PCKINI(NSTR,ITQMT,NTQMT,N4QMT,IPRINT)
C
C     Allocate memory for 3/4 transformed matrix
C
      IF (IPRINT .GE. 1) WRITE (LUPRI,1000) NTQMT*8./(1024.*1024.)
      IF (LFREE .LT. NTQMT) THEN
            WRITE(LUPRI,'(/A,/A,/A,/A,/A,I12,A,/A,I12,A,/A)') 
     &       '########################################################',
     &       '###  Not enough memory for 3-index transf. matrix.   ###',
     &       '###   Increase available memory or reduce memory     ###',
     &       '### requirements by setting .IJTSK to a lower value! ###',
     &       '###          Available memory: ', LFREE,  '          ###',
     &       '###          Required memory : ', NTQMT,  '          ###',
     &       '########################################################'
            CALL QUIT('TD23V2: Not enough memory for 3-index matrix')
      ENDIF
      KFRSV1 = KFREE
      CALL MEMGET('REAL',KTQMAT,NTQMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KTQMAT),NTQMT)
C
C     Get screening matrices
C
      IF (SCRTRA .GT. 0) THEN
         DO I = 1, 4
            ISAMEIJ(I) = I
         ENDDO
         CALL PR4SC2(KDRIJ,WORK,KFREE,LFREE,GABRAO,
     &        WORK(KQ(1,1)),WORK(KQ(1,2)),WORK(KQ(1,3)),WORK(KQ(1,4)),
     &        NDMOQR,ICMOQR,NSTR,ISAMEIJ,IPRINT)
      ELSE
         KDRIJ = KFREE
      ENDIF
C
C     Set up list of shells and sort them
C
      CALL DZERO(TIMDIS,NDSLST)
      CALL SHLLST(NDSLST,IDSLST,INTFLG,INDX,LFREE,IPRINT,SRTSHL)
      IF ((IPRINT .EQ. 0) .AND. (IJBTCH .EQ. 1))
     &     CALL MP2HDR(LUPRI,MXBTCH,NDSLST,-1)
      CALL TIMEC('TD23V2',TIMSHS,TIMSHE)
C ----------------------------------------------------------------------
C     Inner loop over R shells (for all requested integral types)
C        Results are accumulated in the 3/4 transformed matrix
C ----------------------------------------------------------------------
      DO 100 ISHLND = 1, NDSLST
         CALL TIMEC('START ',TIMSHS,TIMSHE)
         IF (IPRINT .GE. 1) THEN
            IF (MXBTCH .GT. 1) THEN
               WRITE(LUPRI,'(2(A,I3,A,I3))')
     &              'Batch ', IJBTCH, '/', MXBTCH,
     &              ' -=- Calculating shell ', ISHLND, '/', NDSLST
            ELSE
               WRITE(LUPRI,'(A,I3,A,I3)')
     &              'Calculating shell ', ISHLND, '/', NDSLST
            END IF
         ELSE IF (IPRINT .EQ. 0) THEN
            CALL GTINFO(TMSTMP)
            WRITE(LUPRI,'(I6,I9,8X,A1,8X,A1,5X,A19)') 
     &           IJBTCH, ISHLND, '-', '-', TMSTMP(1:19)
         END IF
         CALL FLSHFO(LUPRI)
         IF (IDSLST(1,ISHLND) .EQ. 0) GOTO 100
         CALL SHLTRA(IDSLST(1,ISHLND),WORK,KFREE,LFREE,IPRINT,KQ,
     &        INDX,INTFLG,NSTR,NSTRT,DINTSKP,GABRAO,WORK(KDRIJ),
     &        NDMOQR,ICMOQR,ITQMT,NTQMT,KTQMAT,TIM1)
C
C     Print integrals if requested (above threshold TPRI34):
C
         IF (TPRI34 .GT. 0.D0)
     &        CALL PRIINT('The accumulated 3/4 transformed integrals',
     &        NTQMT, WORK(KTQMAT), TPRI34, -1)
C
C     End loop over R shell index
C     
         CALL TIMEC('TD23V2',TIMSHS,TIMSHE) 
         TIMDIS(ISHLND) = TIMSHE - TIMSHS
 100  CONTINUE
C ----------------------------------------------------------------------
C     End inner loop over tasks
C ----------------------------------------------------------------------
      IF (IPRINT .EQ. 0)
     &     WRITE(LUPRI,'(58A1)') ('-', I=1,58)
C
C     Print timings for each R shell
C
      IF (IPRINT .GE. 1) THEN
         WRITE(LUPRI,'(/A/A/A/A)')
     &        '=============================================',
     &        '===          CPU usage per shell          ===',
     &        '=============================================',
     &        ' ISHLA  IC/I2TYP   Tot. dist.   CPU time (s)'
         DO I = 1, NDSLST
            WRITE(LUPRI,'(I6,I4,A1,I2,I15,F13.1)') 
     &           IDSLST(1,I), IDSLST(4,I), '/',
     &           IDSLST(5,I), IDSLST(6,I), TIMDIS(I)
         END DO
         WRITE(LUPRI,'(A/)')
     &        '============================================='
      END IF
C     
C     Allocate memory for 4/4 transformation
C
      IF (IPRINT .GE. 1) WRITE (LUPRI,1010) N4QMT*8./(1024.*1024.)
      CALL MEMGET('REAL',K4QMAT,N4QMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(K4QMAT),N4QMT)
C     
C     Carry out 4/4 transformation
C
      CALL TIMEC('START ',TIMSTR,TIMEND)
      CALL T44INI(INTFLG,ICS2,ICF2)
      CALL LQTRNS(WORK,KFREE,LFREE,IPRINT,ITQMT,ICMOQR(1,4),
     &     NDMOQR(1,1,4),WORK(KTQMAT),NSTR(1,0,3),NSTR(1,0,4),
     &     WORK(KQ(1,4)),WORK(K4QMAT),ICS2,ICF2)
      CALL TIMEC('LQTRNS',TIMSTR,TIMEND) 
      TIM1(4) = TIM1(4) + TIMEND - TIMSTR
C
C     Print integrals if requested (above threshold TPRI44):
C       First the fully transformed integrals from this shell,
C       then the accumulated fully transformed integrals.
C
      IF (TPRI44 .GT. 0.D0)
     &     CALL PRIINT('The fully transformed integrals',
     &     N4QMT, WORK(K4QMAT), TPRI44, -1)
C
C     Allocate memory for naive antisymmetrize
C       Store final integrals in TQMAT if there is enough space,
C       fully transformed integrals are usually few compared to
C       other arrays.
C
      CALL TIMEC('START ',TIMSTR,TIMEND)
      CALL LSZINI(NSTRT,LSIZ)
      IF (NTQMT .GE. LSIZ) THEN
         KANT = KTQMAT
         IF (IPRINT .GE. 1) WRITE(LUPRI,'(/A)') 
     &        ' Storing final integrals in 3/4-matrix'
      ELSE
         IF(IPRINT .GE. 1) WRITE(LUPRI,'(/A)') 
     &        ' Final integrals too large to store in 3/4-matrix'
         CALL MEMGET('REAL',KANT,LSIZ,WORK,KFREE,LFREE)
      ENDIF
      IF (IPRINT .GE. 1) WRITE (LUPRI,1020) LSIZ*8./(1024.*1024.)
      CALL DZERO(WORK(KANT),LSIZ)
C
C     Anti-symmetrize final integrals
C
      CALL NANTSY(WORK,KFREE,LFREE,IPRINT,NSTR,NSTRT,NCLASS,LSIZ,
     &     K4QMAT,KANT)
C
      I1 = KE(1,1) + (ISPNR-1)*LENTSK
      I2 = KE(1,2)
      I3 = KE(1,1) + (JSPNR-1)*LENTSK
      I4 = KE(1,2)
      CALL CONMP2(ERGMP2,WORK(KANT),WORK(I1),WORK(I2),WORK(I3),
     &     WORK(I4),NSTRT(1),NSTRT(2),NSTRT(3),NSTRT(4),
     &     NSTR,ISPNR,JSPNR,IPRINT)
      CALL TIMEC('CONMP2',TIMSTR,TIMEND) 
      TIM1(5) = TIM1(5) + TIMEND - TIMSTR
CTROND
      WRITE(6,*) 'After CONMP2...'
C
      CALL MEMREL('TD23V2',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TD23V2')
      RETURN
C
 1000 FORMAT (/' Storing all 3-index transformed integrals to memory',
     &     /' Memory requirements ',F10.3,' Megabytes')
 1010 FORMAT (/' Storing all fully-transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
 1020 FORMAT (/' Storing all final integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SHLLST */
      SUBROUTINE SHLLST(NDSLST,IDSLST,INTFLG,INDX,LFREE,IPRINT,SORTDS)
C
C     Set up list of shells to be calculated. If requested,
C     the list will be sorted by estimated computational cost.
C
C     Vebjorn Bakken, Summer 2003
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dcbmpt.h"
      DIMENSION IDSLST(6,NDSLST), INDX(3,*), IJP12(0:7)
      LOGICAL SORTDS, DONE
C
C     Count total number of distributions for each R shell
C     Information is stored in IDSLST(n,shell), where:
C
C     n = 1 : ISHLA
C     n = 2 : IBSTRT
C     n = 3 : IBSMAX
C     n = 4 : IC
C     n = 5 : I2TYP
C     n = 6 : Total number of distributions to calculate for ISHLA
C
C     IC:     1 = SS, 2 = LL
C     I2TYP:  1 = SS, 2 = LL, 12 = SS|LL
C
      CALL IZERO(IDSLST,6*NDSLST)
C
      IF (TESTLS) WRITE(LUPRI,'(/A)') 'YES: Doing (LL|SS) ONLY.'
      IF (TESTSL) WRITE(LUPRI,'(/A)') 'YES: Doing (SS|LL) ONLY.'
      DO 100 IC = 1, 2
         IF ((IC.EQ.1).AND.TESTLS) GOTO 100
         IF ((IC.EQ.2).AND.TESTSL) GOTO 100
         CALL SHRNGE(IC,INTFLG,I2TYP,IASTRT,IBSTRT,IASMAX,IBSMAX)
         IF (I2TYP.LT.0) GOTO 100
         DO ISHLA = IASTRT, IASMAX
            IDSLST(1,ISHLA) = ISHLA
            IDSLST(2,ISHLA) = IBSTRT
C     IBSMAX = ISHLA, making use of RS-symmetry
            IDSLST(3,ISHLA) = ISHLA
            IDSLST(4,ISHLA) = IC
            IDSLST(5,ISHLA) = I2TYP
            ITOTDS = 0
            NINSHA  = NINSH(ISHLA,-1,INDX,IPRINT)
            DO ISHLB = IBSTRT, ISHLA
               NINSHB  = NINSH(ISHLB,-1,INDX,IPRINT)
               IF(MAXSCL.EQ.0) THEN
                  MWHER = MXMEMHER(IC,I2TYP)
                  NSIZG = LFREE - MIN(LFREE/5, 3 000 000) - MWHER
               ELSE
                  NSIZG = MAXSCL
               ENDIF
               ITOTDS = ITOTDS + NDISTRN(I2TYP,IPRINT,.FALSE.,.FALSE.,
     &              .FALSE.,.TRUE.,INDX,IDUM,IJP12,IDUM,NSIZG,IDUM)
            END DO
            IDSLST(6,ISHLA) = ITOTDS
         END DO
 100  CONTINUE
C
C     Sort R shells by the total number of distributions to be calculated
C     (includes loop over S shells)
C
      IF (SORTDS) THEN
 200     CONTINUE
         DONE = .TRUE.
         DO ISHL = 1, NDSLST-1
C
C     Simple bubble sort, but weigthed so that I2TYP = SS|LL is
C     multiplied by a factor of 4 compared to I2TYP = SS or LL
C
            IWGHT1 = NINT(LOG(1.D0*IDSLST(5,ISHL))/LOG(10.D0)) + 1
            IWGHT2 = NINT(LOG(1.D0*IDSLST(5,ISHL+1))/LOG(10.D0)) + 1
            IF (IDSLST(6,ISHL+1)*IWGHT2*IWGHT2 .GT.
     &           IDSLST(6,ISHL)*IWGHT1*IWGHT1) THEN
               DO I = 1, 6
                  ITMP = IDSLST(I,ISHL)
                  IDSLST(I,ISHL)   = IDSLST(I,ISHL+1)
                  IDSLST(I,ISHL+1) = ITMP
               END DO
               DONE = .FALSE.
            END IF
         END DO
         IF (.NOT. DONE) GOTO 200
      END IF
      IF (IPRINT .GE. 1) THEN
         WRITE(LUPRI,'(/A/A/A/A)')
     &        '===============================================',
     &        '===              List of shells             ===',
     &        '===============================================',
     &        ' ISHLA IBSTRT IBSMAX   IC I2TYP   No. distrib.'
         DO J = 1, NDSLST
            IF (IDSLST(1,J) .GT. 0) 
     &           WRITE(LUPRI,'(I6,2I7,2I5,I13)') (IDSLST(I,J), I = 1,6)
         END DO
         WRITE(LUPRI,'(A/)')
     &        '==============================================='
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SHLTRA */
      SUBROUTINE SHLTRA(IDSNFO,WORK,KFREE,LFREE,IPRINT,KQ,
     &     INDX,INTFLG,NSTR,NSTRT,DINTSKP,GABRAO,DRIJ,
     &     NDMOQR,ICMOQR,ITQMT,NTQMT,KTQMAT,TIM1)
C
C     Calculate contribution to 3/4-transformed matrix from _one_
C     R shell as specified in IDSNFO
C
C     Vebjorn Bakken, Summer 2003
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dcbtra.h"
#include "dcborb.h"
#include "dcbmpt.h"
#include "dcbmp2.h"
#include "dgroup.h"
      DIMENSION IDSNFO(6), KQ(2,4), KE(2,2), INDX(3,*),
     &     NSTR(2,0:2,4), NSTRT(4), DINTSKP(*), GABRAO(*), WORK(*),
     &     NDMOQR(2,2,4), ICMOQR(2,4), TIM1(7), DRIJ(*),
     &     ITQMT(0:7,0:7,2), IJP12(0:7), IHM(0:7), ISAMEIJ(4)
      CHARACTER*(24) TMSTMP
      LOGICAL NODV, NOPV, L2TYP(0:3)
C
      CALL QENTER('SHLTRA')
      KFRSV1 = KFREE
C
C     Do a specific R shell (for all requested integral types)
C
      ISHLA  = IDSNFO(1)
      IBSTRT = IDSNFO(2)
      IBSMAX = IDSNFO(3)
      IC     = IDSNFO(4)
      I2TYP  = IDSNFO(5)
      ITOTDS = IDSNFO(6)
C
      IF (ITOTDS .LE. 0) GOTO 100
      IF (I2TYP .LT. 0) THEN
         GOTO 100
      ELSEIF (I2TYP .EQ. 1) THEN
         ICS = 1
         ICF = 1
      ELSEIF (I2TYP .EQ. 2) THEN
         ICS = 2
         ICF = 2
      ELSE
         ICS = 1
         ICF = 2
      ENDIF
C
C     Loop over S shells
C
      DO ISHLB = IBSTRT, IBSMAX
         KFRSV2 = KFREE
         NINSHA  = NINSH(ISHLA,-1,INDX,IPRINT)
         NINSHB  = NINSH(ISHLB,-1,INDX,IPRINT)
C
C     First a dummy run to calculate the memory requirements
C     for GMAT. Max size of GMAT is set in MAXSCL from
C     .SCLMEM in input or chosen based on MWHER.
C
         IF(MAXSCL.EQ.0) THEN

            MWHER = MXMEMHER(IC,I2TYP)
            MEMSV1 = MIN(LFREE/10, 1 000 000) + MWHER
C     Subtract max mem needed in Hermit and 10%/1Mw for anything else
            NSIZG = LFREE - MEMSV1
         ELSE
            NSIZG = MAXSCL
            IF (MAXSCL.GT.LFREE/100*80) WRITE(LUPRI,'(2(1X,A,I10))')
     &           '.GT. 80%.., LFREE:',LFREE,'MAXSCL:',MAXSCL
            IF (MAXSCL.GT.LFREE/100*95) THEN
               WRITE(LUPRI,'(2(1X,A,I10))')
     &              '.GT. 95%.., LFREE:',LFREE,'MAXSCL:',MAXSCL
               CALL QUIT('GMAT too large...')
            ENDIF
         ENDIF
C
C     Determine number of passes based on available memory
C
         NUMDIS = NDISTRN(I2TYP,IPRINT,.FALSE.,.FALSE.,
     &        .FALSE.,.TRUE.,INDX,IDUM,IJP12,IDUM,NSIZG,IDUM)
C
C     If several passes are needed, memory requirements
C     must be reconsidered
C
         IF (NPASS.GT.1) THEN
            NGBTMP = (IRAT*NSIZG)/((1+IRAT)*NPASS)
            MEMSV2 = 2*NGBTMP
            NSIZG = LFREE - MAX(MEMSV1,MEMSV2)
            NUMDIS = NDISTRN(I2TYP,IPRINT,.FALSE.,.FALSE.,
     &           .FALSE.,.TRUE.,INDX,IDUM,IJP12,IDUM,NSIZG,IDUM)
         END IF
C
         LDXAB = 5*NINSHA*NINSHB
         CALL MEMGET('INTE',KDXAB,LDXAB,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KIJPASS,8*NPASS,WORK,KFREE,LFREE)
C
C     Set INDXAB and INDPASS
C
         NUMDIS = NDISTRN(I2TYP,IPRINT,.TRUE.,.TRUE.,.FALSE.,.TRUE.,
     &        INDX,WORK(KDXAB),IJP12,IDUM,NSIZG,WORK(KIJPASS))
         NSIZH = 0
         DO IREPAB = 0, NBSYM-1
            IHM(IREPAB) = NSIZH
            IREPIJ      = JBTOF(IREPAB,1)
            NSIZH       = NSIZH + IJP12(IREPAB)*NFPCK12(IREPIJ)
         ENDDO
         CALL MEMGET('REAL',KHMAT,NSIZH*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KHMAT),NSIZH*NZ)
         IF (IPRINT .GE. 3) THEN
            CALL GTINFO(TMSTMP)
            WRITE (LUPRI,
     &           '(/A/A,A/2(A,2I10/),2(A,I10/),2(A,I10,F10.3/))')
     &           '   SHLTRA',
     &           '   Time stamp:                      ', TMSTMP(1:19),
     &           '   Calculation of integral type:   ', IC,I2TYP,
     &           '   Integral distribution set:      ', ISHLA,ISHLB,
     &           '   Number of distributions in set: ', NUMDIS,
     &           '   Number of passes:               ', NPASS,
     &           '   Size of scalar integral batch:  ', NSIZG,
     &           NSIZG*8./(1024.*1024.),
     &           '   Size of half-transf. int. batch:', NSIZH,
     &           NSIZH*8./(1024.*1024.)
         END IF
         IF (NPASS.GT.1) THEN
C
C     Initialize scalar integral buffer file and arrays
C     Use the space reserved for GMAT now as buffer storage
C
            LGFIL = LUTRA2
            NGBFSZ = (IRAT*NSIZG)/((1+IRAT)*NPASS)
            CALL MEMGET('INTE',KALGREC,NPASS,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KALGBUF,NPASS,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KAIGBUF,NPASS*NGBFSZ,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KARGBUF,NPASS*NGBFSZ,WORK,KFREE,LFREE)
            KGMAT  = KALGREC
            KLGREC = KALGREC - KGMAT + 1
            KLGBUF = KALGBUF - KGMAT + 1
            KIGBUF = KAIGBUF - KGMAT + 1
            KRGBUF = KARGBUF - KGMAT + 1
         ELSE
C
C     Write directly into the array GMAT
C
            CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KGMAT),NSIZG)
            KLGBUF  = KFREE
            KALGREC = KFREE
         ENDIF
C
C     Calculate current distribution
C
         NODV = NASHT.EQ.0
         NOPV = NASHT.LT.2
         CALL TIMEC('START ',TIMSTR,TIMEND)
         CALL CALDIS(I2TYP,WORK(KGMAT),INDX,WORK(KDXAB),NODV,NOPV,
     &        GABRAO,DRIJ,DINTSKP,SCRTRA,
     &        WORK(KFREE),LFREE,IPRINT)
         CALL TIMEC('CALDIS',TIMSTR,TIMEND) 
         TIM1(1) = TIM1(1) + TIMEND - TIMSTR
C
C     De-allocate buffers and allocate GMAT if necessary.
C
         IF (NPASS.GT.1) THEN
            CALL MEMREL('SHLTRA.1',WORK,1,KALGBUF,KFREE,LFREE)
            CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KGMAT),NSIZG)
         ENDIF
C
C     Do the transformation of the first index pair and sort the
C     half-transformed integrals to memory
C
         CALL TIMEC('START ',TIMSTR,TIMEND)
         CALL MS4IN1E(WORK,KFREE,LFREE,IPRINT,WORK(KALGREC),ICS,ICF,
     &        NSTR(1,0,1),NSTR(1,0,2),NDMOQR,ICMOQR,WORK(KIJPASS),
     &        WORK(KQ(1,1)),WORK(KQ(1,2)),WORK(KGMAT),NSIZG,WORK(KHMAT)) 
         CALL TIMEC('MS4IN1',TIMSTR,TIMEND) 
         TIM1(2) = TIM1(2) + TIMEND - TIMSTR
         CALL MEMREL('SHLTRA.2',WORK,1,KGMAT,KFREE,LFREE)
         IF (NPASS.GT.1) CALL DELGBUF(LGFIL)
C
C     3/4 transformation
C
         CALL TIMEC('START ',TIMSTR,TIMEND)
         CALL TQTRNS(WORK,KFREE,LFREE,IPRINT,IC,IHM,INDX,WORK(KDXAB),
     &        WORK(KHMAT),ITQMT,ICMOQR(1,3),NDMOQR(1,1,3),WORK(KTQMAT),
     &        NSTR(1,0,3),WORK(KQ(1,3)))
         CALL TIMEC('TQTRNS',TIMSTR,TIMEND) 
         TIM1(3) = TIM1(3) + TIMEND - TIMSTR
         CALL MEMREL('SHLTRA',WORK,1,KFRSV2,KFREE,LFREE)
      ENDDO
C
C     End of loop over S shell index
C
 100  CONTINUE
      CALL MEMREL('SHLTRA',WORK,1,KFRSV1,KFREE,LFREE)
      CALL QEXIT('SHLTRA')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck NANTSY */
      SUBROUTINE NANTSY(WORK,KFREE,LFREE,IPRINT,NSTR,NSTRT,NCLASS,LSIZ,
     &     K4QMAT,KANT)
C
C     Naive anti-symmetrize (extracted from TRDR23)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "symmet.h"
      DIMENSION NSTR(2,0:2,4),NSTRT(4),WORK(*)
C
C
C     FQMAT stored as (ij,iz1,kl,iz2,Pk,Grs,Pkl) and the order
C     of ij (and kl) is gg,uu,ug,gu or as (i,j,iz1,k,l,iz2,Pk,Grs,Pkl)
C
C     F.ex. for Pkl=g, iz2 loops faster than loop over gg/uu for kl,
C     but slower than kl within gg/uu
C
      IOFF3 = 0
      DO IREPKL = 1,NFSYM
         CALL MEMGET('REAL',KFTMAT,2*NFPCK12(IREPKL),WORK,KFREE,LFREE)
         LKR = 0
         DO IREPL = 1, NFSYM
            IREPK = MOD(IREPL+IREPKL,2) + 1
C     Offset for k and l symmetry (i.e. for Pkl=g, offset for uu)
            IOFF2 = IFPCK34(IREPK,IREPL)*NFPCK12(IREPKL)*NZ*NZ
            DO L = 1, NSTR(IREPL,0,4)
               LKR = LKR + 1
               DO K = 1, NSTR(IREPK,0,3)
                  IF (IREPK.EQ.1) THEN
                     KKR = K
                  ELSE
                     KKR = NSTR(1,0,3) + K
                  ENDIF
C     
C     Real and complex classes
C     
C     Offset for K and L
                  IOFF1 = NFPCK12(IREPKL)*NZ*
     &                 (NSTR(IREPK,0,3)*(L-1)+(K-1))
                  DO JCLASS = 1, 4
                     CALL DZERO(WORK(KFTMAT),2*NFPCK12(IREPKL))
                     ICNT1 = 0
                     DO IREPRS = 0, NBSYM-1
                        IF (IREPKL.EQ.JBTOF(IREPRS,1)) THEN
C     Counter for size of IREPKL=1 contribution and offset for
C     IREPKL
                           IOFF6 = ICNT1
                           ICNT1=ICNT1
     &                          +NFPCK12(IREPKL)*NFPCK34(IREPKL)*NZ*NZ
C     Offset for Grs
                           IRPCN1 = 0
                           DO IZ2 = 1, NZ
C     Offset for iz2 
                              IOFF4 = NFPCK12(IREPKL)*NZ*IRPCN1*
     &                             NSTR(IREPK,0,3)*NSTR(IREPL,0,4)
C
                              IRPCN1 = IRPCN1 + 1
                              IQ2 = IPQTOQ(IZ2,IREPRS)
                              IRPCN2 = 0
                              DO IZ1 = 1, NZ
C     Offest for iz1
                                 IOFF5 = NFPCK12(IREPKL)*IRPCN2
                                 IRPCN2 = IRPCN2 + 1
                                 IQ1 = IPQTOQ(IZ1,IREPRS)
                                 IFAC=INTCL(JCLASS,IQ1,IQ2)
                                 IOFF=IOFF1+IOFF2+IOFF3
                                 IOFF=IOFF+IOFF4+IOFF5+IOFF6
C
                                 IF (IFAC.NE.0) THEN
                                    FACTR=IFAC
                                    CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                   WORK(K4QMAT+IOFF),1,
     &                                   WORK(KFTMAT),1)
                                 ENDIF
C 
C     Non-real groups
C
                                 IF (NZ.GT.1) THEN
                                    IFAC=INTCL(JCLASS+4,IQ1,IQ2)
                                    IF (IFAC.NE.0) THEN
                                       FACTR=IFAC
                                       CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                      WORK(K4QMAT+IOFF),1,
     &                                      WORK(KFTMAT+
     &                                      NFPCK12(IREPKL)),1)
                                    ENDIF
                                 ENDIF
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDDO
                     CALL MP2FI1(JCLASS,IREPKL,NSTR,KKR,LKR,
     &                    WORK(KFTMAT),WORK(KANT),
     &                    NFPCK12(IREPKL),NSTRT(1),NSTRT(2),
     &                    NSTRT(3),NSTRT(4),LSIZ)
                  ENDDO
C     
C     Quaternion groups
C     
                  IF(NZ.EQ.4) THEN
                     DO JCLASS = 9, 12
                        CALL DZERO(WORK(KFTMAT),2*NFPCK12(IREPKL))
                        ICNT1 = 0
                        DO IREPRS = 0, NBSYM-1
                           IF (IREPKL.EQ.JBTOF(IREPRS,1)) THEN
C     Counter for size of IREPKL=1 contribution and offset for
C     IREPKL
                              IOFF6 = ICNT1
                              ICNT1=ICNT1
     &                             +NFPCK12(IREPKL)
     &                             *NFPCK34(IREPKL)*NZ*NZ
C     Offset for Grs
                              IRPCN1 = 0
                              DO IZ2 = 1, NZ
C     Offset for iz2 
                                 IOFF4 = NFPCK12(IREPKL)*NZ*IRPCN1*
     &                                NSTR(IREPK,0,3)*NSTR(IREPL,0,4)
C     
                                 IRPCN1 = IRPCN1 + 1
                                 IQ2 = IPQTOQ(IZ2,IREPRS)
                                 IRPCN2 = 0
                                 DO IZ1 = 1, NZ
C     Offest for iz1
                                    IOFF5 = NFPCK12(IREPKL)*IRPCN2
                                    IRPCN2 = IRPCN2 + 1
                                    IQ1 = IPQTOQ(IZ1,IREPRS)
                                    IFAC=INTCL(JCLASS,IQ1,IQ2)
                                    IOFF=IOFF1+IOFF2+IOFF3
                                    IOFF=IOFF+IOFF4+IOFF5+IOFF6
C
                                    IF (IFAC.NE.0) THEN
                                       FACTR=IFAC
                                       CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                      WORK(K4QMAT+IOFF),1,
     &                                      WORK(KFTMAT),1)
                                    ENDIF
                                    IFAC=INTCL(JCLASS+4,IQ1,IQ2)
                                    IF (IFAC.NE.0) THEN
                                       FACTR=IFAC
                                       CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                      WORK(K4QMAT+IOFF),1,
     &                                      WORK(KFTMAT
     &                                      +NFPCK12(IREPKL)),1)
                                    ENDIF
                                 ENDDO
                              ENDDO
                           ENDIF
                        ENDDO
                        CALL MP2FI1(JCLASS,IREPKL,NSTR,KKR,LKR,
     &                       WORK(KFTMAT),WORK(KANT),
     &                       NFPCK12(IREPKL),NSTRT(1),NSTRT(2),
     &                       NSTRT(3),NSTRT(4),LSIZ)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CALL MEMREL('NANTSY',WORK,1,KFTMAT,KFREE,LFREE)
         IOFF3=ICNT1
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TQTRNS */
      SUBROUTINE TQTRNS (WORK,KFREE,LFREE,IPRINT,IC,IHM,
     &     INDX,INDXAB,HMAT,ITQMT,ICMOQR,NDMOQR,
     &     TQMAT,NSTR3,Q3)
C
C     Written by Jon Laerdahl, June 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "symmet.h"
C
      PARAMETER (D0 = 0.0D00, D1=1.0D00)
      INTEGER   NSTR3(2)
      DIMENSION NDMOQR(2,2),ICMOQR(2)
      DIMENSION INDX(3,*),INDXAB(NINSHA,NINSHB,5)
      DIMENSION ITQMT(0:7,0:7,2)
      DIMENSION IHM(0:7)
      DIMENSION HMAT(*),TQMAT(*),Q3(*),QCOF(4)
      DIMENSION WORK(*)
      DIMENSION IND(2)
C
#include "ibtfun.h"
C     
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
      CALL QENTER('TQTRNS')
C
      KFRSAV = KFREE
      NIJMAX = 0
C     
C     TQMAT is stored as (ij,iz1,k,s,iz2,ic,Gs,Grs).
C     For real groups this is identical to
C     (ij,iz1,k,iz2,s,ic,Gs,Grs) which is generated
C     directly in QAXPY. For non-real groups we first
C     generate (ij,iz1,k,iz2* which is accumulated in TQMAT.
C     The correct format for TQMAT is needed in LQTRNS
C     for QGEMM: (n,s,iz)*(s,l,iz) = (n,l,iz)
C
      DO IFRP=1,NFSYM
         IF(NFPCK12(IFRP)*NZ**2.GT.NIJMAX) 
     &        NIJMAX = NFPCK12(IFRP)*NZ**2
      ENDDO
      CALL MEMGET('REAL',KT3MAT,NIJMAX,WORK,KFREE,LFREE)
C
      DO20 INDAR = 1, NINSHA
         DO 10 INDBR = 1, NINSHB
            INDAB  = INDXAB(INDAR,INDBR,1)
            IF (INDAB.EQ.0) GOTO 10
            IREPAB = INDXAB(INDAR,INDBR,3)
C     Pointer AP	
            IPAB   = INDXAB(INDAR,INDBR,4)
            IREPIJ = JBTOF(IREPAB,1)
            NIJNZ1 = NFPCK12(IREPIJ)*NZ
            NIJNZ2 = NFPCK12(IREPIJ)*NZ**2
C     IOFF for HMAT
            IOFF = (IHM(IREPAB)+(IPAB-1)*NFPCK12(IREPIJ))*NZ+1
            CALL IUNPCK(INDAB,2,IND)
            INDA = IND(1)
            INDB = IND(2)
            IREPA = INDX(2,INDA)
            IREPB = INDX(2,INDB)
C     This is actual index R and S in this irrep and IC.
            INDA1 = INDA - ICOS(IREPA+1,IC)
            INDB1 = INDB - ICOS(IREPB+1,IC)
C
C     Loop over K-transformations. 
C     Accumulate first C_rk(rs), then C_sk(rs).
C     Don't accumulate twice if functions A.EQ.B
C
            KLST = 2
            IF (INDA.EQ.INDB) KLST = 1
            INDR = INDA1
            INDS = INDB1
            IREPR = IREPA
            IREPS = IREPB
            DO KTYP = 1,KLST
               IREPK = JBTOF(IREPR,IC)
               NK = NSTR3(IREPK)
               NS = NBBAS(IREPS,IC)
               NIJKS = NK*NS*NIJNZ1
C     
C     Offset for TQMAT.
C     
               JOFF1 = ITQMT(IREPAB,IREPS,IC) 
     &              + (INDS-1)*NK*NIJNZ1+1
C     
C     Offset for Coeffs
C     
               KOFF1 = ICMOQR(IREPK) + IBBAS(IREPR,IC)
     &               - IBAS(IREPK) + INDR - 1
C
C     Dimensions of Q3
C
               NRQ1 = NDMOQR(1,IREPK)
               NCQ1 = NDMOQR(2,IREPK)
               IF (IC.EQ.1) THEN
                  IREPR2 = IREPR
               ELSE
                  IREPR2 = IBTXOR(IXYZ,IREPR)
               ENDIF
C     
C     Start looping over K index
C
               DO K = 1,NK
C     
                  JOFF2 = JOFF1 + (K-1)*NIJNZ1
                  KOFF = KOFF1 + (K-1)*NRQ1
                  DO I = 1,NZ
                     QCOF(I) = Q3(KOFF+(I-1)*NRQ1*NCQ1)
                  ENDDO
C     
C     Non-real groups
C
                  IF (NZ.GT.1) THEN
                     CALL DZERO(WORK(KT3MAT),NIJNZ2)
                     CALL QAXPY(NIJNZ1,'C','N',QCOF,IPQTOQ(1,IREPR2),
     &                    NZ,HMAT(IOFF),NIJNZ1,1,IPQTOQ(1,0),1,
     &                    WORK(KT3MAT),NIJNZ1,1,IPQTOQ(1,IREPR2),
     &                    NZ)
C     
                     DO IZ = 1,NZ
                        JOFF = JOFF2 + (IZ-1)*NIJKS
                        CALL DAXPY(NIJNZ1,D1,
     &                       WORK(KT3MAT+NIJNZ1*(IZ-1)),
     &                       1,TQMAT(JOFF),1)
                     ENDDO
C
C     Real groups
C
                  ELSE
                     CALL QAXPY(NIJNZ1,'C','N',QCOF,IPQTOQ(1,IREPR2),
     &                    NZ,HMAT(IOFF),NIJNZ1,1,IPQTOQ(1,0),NZ,
     &                    TQMAT(JOFF2),NIJNZ1,1,IPQTOQ(1,IREPR2),NZ)
                  ENDIF
               ENDDO
C     
C     Accumulate C_sk(rs)
C     
               INDR = INDB1
               INDS = INDA1
               IREPR = IREPB
               IREPS = IREPA
            ENDDO
 10      CONTINUE
 20   CONTINUE
      CALL MEMREL('TQTRNS',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TQTRNS')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck LQTRNS */
      SUBROUTINE LQTRNS (WORK,KFREE,LFREE,IPRINT,
     &     ITQMT,ICMOQR,NDMOQR,
     &     TQMAT,NSTR3,NSTR4,Q4,FQMAT,ICS,ICF)
C
C     Written by Jon Laerdahl, June 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "dcbbas.h"
#include "dcbtra.h"
#include "symmet.h"
C
      PARAMETER (D0 = 0.0D00, D1=1.0D00)
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION NDMOQR(2,2),ICMOQR(2)
      DIMENSION ITQMT(0:7,0:7,2)
      DIMENSION TQMAT(*),Q4(*),FQMAT(*)
      DIMENSION WORK(*)
C
#include "ibtfun.h"
C
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
      CALL QENTER('LQTRNS')
C
      IOFF=1
      DO IREPIJ = 1,NFSYM
         NIJ = NFPCK12(IREPIJ)
         NKL = NFPCK34(IREPIJ)
         DO IREPRS = 0, NBSYM-1
            IF (IREPIJ.EQ.JBTOF(IREPRS,1)) THEN
               DO IREPS = 0, NBSYM-1
                  DO IC = ICS,ICF
                     IREPR = IBTXOR(IREPRS,IREPS)
                     IREPK = JBTOF(IREPR,IC)
                     IREPL = JBTOF(IREPS,IC)
                     NK = NSTR3(IREPK)
                     NL = NSTR4(IREPL)
                     N1 = NIJ*NK*NZ**2
                     N1 = NIJ*NK*NZ
                     N2 = NBBAS(IREPS,IC)
                     N3 = NL
                     IF (NBBAS(IREPR,IC).EQ.0) GOTO 100
C     
C     Offset for TQMAT.
C     
                     JOFF = ITQMT(IREPRS,IREPS,IC) + 1
C     
C     Offset for FQMAT.
C     
                     IOFF1 = IOFF + 
     &                    IFPCK34(IREPK,IREPL)*NFPCK12(IREPIJ)*NZ**2
C
C     Offset for Coeffs
C     
                     KOFF = ICMOQR(IREPL) + IBBAS(IREPS,IC) 
     &                    - IBAS(IREPL)
C     
C     Dimensions of Q4
C     
                     NRQ1 = NDMOQR(1,IREPL)
                     NCQ1 = NDMOQR(2,IREPL)
C
                     IF (IC.EQ.1) THEN
                        IREPR2 = IREPR
                     ELSE
                        IREPR2 = IBTXOR(IXYZ,IREPR)
                     ENDIF
C
                     IF (IC.EQ.1) THEN
                        IREPS2 = IREPS
                     ELSE
                        IREPS2 = IBTXOR(IXYZ,IREPS)
                     ENDIF
C     
                     IF (N1*N2.GT.0) THEN
C
CTROND  check this.... TODO
                        CALL QGEMM(N1,N3,N2,D1,
     &                       'N','N',IPQTOQ(1,IREPR2),TQMAT(JOFF),
     &                       N1,N2,NZ,
     &                       'N','N',IPQTOQ(1,IREPS2),Q4(KOFF),
     &                       NRQ1,NCQ1,NZ,
     &                       D1,IPQTOQ(1,IREPRS),FQMAT(IOFF1),
     &                       N1,N3,NZ)
                     ENDIF
 100                 CONTINUE
                  ENDDO
               ENDDO
               IOFF = IOFF + NIJ*NKL*NZ*NZ
            ENDIF
         ENDDO
      ENDDO
      CALL QEXIT('LQTRNS')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MP2FI1 */
      SUBROUTINE MP2FI1(ICLASS,IREPKL,NSTR,KKR,LKR,FMAT,SPMAT,
     &     NFM,NI,NJ,NK,NL,LS)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
      DIMENSION NSTR(2,0:2,4)
      DIMENSION FMAT(NFM,2)
      DIMENSION SPMAT(NI,NJ,NK,NL,LS)
C
c       WRITE(LUPRI,'(A,4I4)') 'MP2FI1:', NI,NJ,NK,NL
      M = 0
      JJ = 0
      DO IREPJ = 1, NFSYM
         IREPI = MOD(IREPJ+IREPKL,2) + 1
         DO J = 1, NSTR(IREPJ,0,2)
            JJ = JJ + 1
            DO I = 1, NSTR(IREPI,0,1)
               IF (IREPI.EQ.1) THEN
                  II = I
               ELSE
                  II = I + NSTR(1,0,1)
               ENDIF
               M = M + 1
c               WRITE(LUPRI,'(I4,A,4I4,2X,F20.10)') 
c     &              ICLASS,'##',II,JJ,KKR,LKR,FMAT(M,1)
               SPMAT(II,JJ,KKR,LKR,ICLASS)  =FMAT(M,1)
               SPMAT(II,JJ,KKR,LKR,ICLASS+4)=FMAT(M,2)
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRIINT */
      SUBROUTINE PRIINT(TEXT,NUMINT,VALINT,THRESH,NODE)
C
C     Print out all integrals above THRESH
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) TEXT
      DIMENSION VALINT(NUMINT)
      IINT = 0
      IF (NODE .EQ. 0) THEN
         WRITE(LUPRI,'(/1X,2A/A,I10)') TEXT, ' on master',
     &        ' Number of integrals: ',NUMINT
      ELSE IF (NODE .GT. 0) THEN
         WRITE(LUPRI,'(/1X,2A,I4/A,I10)') TEXT, ' on node #', NODE,
     &        ' Number of integrals: ',NUMINT
      ELSE
         WRITE(LUPRI,'(/1X,A/A,I10)') TEXT,' Number of integrals: ',
     &        NUMINT
      END IF
      DO I = 1, NUMINT
         IF (ABS(VALINT(I)) .GT. THRESH) THEN
            WRITE(LUPRI,*) I, VALINT(I)
            IINT = IINT + 1
         END IF
      END DO
      IF (IINT .EQ. 0) THEN
         WRITE(LUPRI,'(A)') ' No integrals above threshold.'
      ELSE
         WRITE(LUPRI,'(A,I10/)')
     &        ' Number of integrals above threshold: ', IINT
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TIMEC */
      SUBROUTINE TIMEC(TEXT,TIMSTR,TIMEND)
#include "implicit.h"
#include "priunit.h"
      CHARACTER TEXT*6
C
      IF (TEXT .EQ. 'START ') THEN
         TIMSTR = SECOND()
      ELSE
         TIMEND = SECOND()
         TIME   = TIMEND - TIMSTR
C        CALL TIMTXT(' Time used in '//TEXT//' is',TIME, LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TSKCHK */
      SUBROUTINE TSKCHK(MAXTSK,LENTSK,LWORK,NODES,NOCC,IPRINT)
C
C     Written by Jon Laerdahl, June 1997
C     Parallel polish by T.Saue Sep 18 1998
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbmpt.h"
C
C     Maximum number of active occupied spinors in each
C     task is either set in input (IJTSK) or
C     on the basis of avaliable memory.
C
      WRITE(LUPRI,'(A,I10/)') 
     &     'Memory available for Direct MP2 is (Words):',
     &     LWORK
      IF (IJTSK.GT.0) THEN
        LENTSK = MIN(IJTSK,NOCC)
        IF (IPRINT.GE.0) THEN
           WRITE(LUPRI,'(A,/A,I10/)') 
     &          'Maximum number of occupied spinors I and J in batch',
     &          'set in input:', IJTSK
        ENDIF
      ELSE
        IF(NODES.GT.1) THEN
C....total number of tasks divided by number of nodes
           MAXTSK = NOCC*(NOCC+1)/2
           LENTSK = 1
           DO IJTSK = 1,NOCC
             NUM = INT(NOCC/IJTSK)
             IF(MOD(NOCC,IJTSK).NE.0) NUM = NUM + 1
             NTSK = NUM*(NUM+1)/2
             IF(NTSK.LT.NODES) GOTO 10
             IF(NTSK.NE.MAXTSK) THEN
               MAXTSK = NTSK
               LENTSK = IJTSK
             ENDIF
           ENDDO
  10       CONTINUE
           IF(NODES.GT.MAXTSK) THEN
             WRITE(LUPRI,'(A)') 
     &        '** WARNING: More nodes than tasks !'
           ENDIF
        ELSE
          LENTSK = NOCC
        ENDIF 
        WRITE(LUPRI,'(A,/A,I6,A/)') 
     &     'Maximum number of occupied spinors in batch not set'//
     &     ' in input.', 
     &     'Testing ALL spinors :',NOCC,' in batch.'
        IF (IPRINT.GE.1) THEN
          WRITE(LUPRI,'(A/A/A/A/A/)')
     &          '==========================================',
     &          '=== If not enough memory is avaliable: ===',
     &          '===   Reduce memory requirements by    ===',
     &          '===  setting .IJTSK to a lower value   ===',
     &          '=========================================='
        ENDIF
      ENDIF
      IF ((MOD(NOCC,LENTSK)).EQ.0) THEN
         NUM = NOCC/LENTSK
         IF (IPRINT.GE.1) 
     &        WRITE(LUPRI,'(/A)') 'All Batches have same Size.'
      ELSE
         NUM = NOCC/LENTSK + 1
         IF (IPRINT.GE.1) 
     &        WRITE(LUPRI,'(/A)') 'Last Batch has reduced Size.'
      ENDIF
      MAXTSK = NUM*(NUM+1)/2
      IF (IPRINT.GE.1) THEN
         WRITE(LUPRI,'(/A,I5)')'Number of tasks    :', MAXTSK
         WRITE(LUPRI,'(A,I5)') 'Number in Max Batch:', LENTSK
         CALL FLSHFO(LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MS4IN1E*/
      SUBROUTINE MS4IN1E(WORK,KFREE,LFREE,IPRINT,LGREC,ICS,ICF,
     &                   NSTR1,NSTR2,NDMOQR,ICMOQR,IJPASS,Q1,Q2,
     &                   GMAT,NSIZG,HMAT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher May 1997.
C
C     PURPOSE : Driver for the first index transformation.
C               Loop over batches for a given shell combination and
C               write the half-transformed integrals to memory
C
C     Input :
C
C     - IPRINT       Print flag
C     - NPASS        Number of passes through the scalar integrals
C     - ICS          First class of integrals 1 : (LL|XX), 2 : (SS|XX)
C     - ICL          Last class of integrals
C     - NSTR1        Number of active spinors for index 1
C     - NSTR2        Number of active spinors for index 2
C     - NFPCK12      Number of spinor pairs for each compound symmetry
C     - IFPCK12      Pointers to spinor pairs
C     - NDMOQR       Dimensions of the coefficient array
C     - ICMOQR       Pointers to coefficients
C     - IJPASS       Number of blocks for each boson symmetry
C     - Q1           Coefficients for index 1
C     - Q2           Coefficients for index 2
C     - GMAT         Symmetry packed scalar integrals
C     - NSIZG        Size of GMAT array
C     - MTHM         Number of scalar integral pairs that contribute
C                    to a given fermion ircop (multiplied with NZ)
C     - ITHM         Start address of a boson representation on the
C                    file of half-transformed integrals.
C     - INDX         For each boson function the boson irrep and the
C                    position in the shell
C
C     Output :
C
C     -              Updated values of buffer blocks
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
      INTEGER   NSTR1(2),NSTR2(2)
      DIMENSION WORK(*)
      DIMENSION Q1(*),Q2(*),GMAT(NSIZG)
      DIMENSION IJPASS(0:7,NPASS)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2)
      DIMENSION HMAT(*)
C
C     For the buffered input
C
      DIMENSION LGREC(NPASS)
C
C     For the buffered output
C
      DIMENSION IHM(0:7)
C
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbibt.h"
C
      CALL QENTER('MS4IN1E')
      KFRSAV = KFREE
C
      IHOFF = 1
      DO IPASS = 1, NPASS
C
C        Get the scalar integrals that are processed in this pass
C        The integral are already in GMAT when NPASS = 1
C
          IF (NPASS.GT.1) THEN
             CALL DZERO(GMAT,NSIZG)
             IREC = LGREC(IPASS)
             CALL READGBF(LGFIL,IREC,NGBFSZ,GMAT)
          ENDIF
C
C        Calculate the size of the H-matrix for this pass
C
         NSIZH = 0
         DO IREPAB = 0, NBSYM-1
            IHM(IREPAB) = NSIZH
            IREPIJ = JBTOF(IREPAB,1)
            NSIZH = NSIZH + IJPASS(IREPAB,IPASS)*NFPCK12(IREPIJ)
         ENDDO
C
C        Do first step of 4-index transformation :
C        Transform first pair of indices
C
         CALL MS4IN1 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,
     &                .FALSE.,NSTR1,NSTR2,NDMOQR,
     &                ICMOQR,IJPASS(0,IPASS),Q1,Q2,
     &                GMAT,HMAT(IHOFF))
         IHOFF = IHOFF + NSIZH*NZ
C
C        We can also rewrite code and call the 3/4 transformation here.
C        This reduces the memory requirements for HMAT.
C
      ENDDO
C
      CALL QEXIT('MS4IN1E')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRIEIN */
      SUBROUTINE PRIEIN(E1,E2,E3,E4,N1,N2,N3,N4,NSTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
      DIMENSION E1(N1),E2(N2),E3(N3),E4(N4)
      DIMENSION NSTR(2,0:2,4)
      NX = 0
      DO IFRP = 1, NFSYM
         WRITE(LUPRI,'(A,I1)') 'Index 1, Irrep ',IFRP
         WRITE(LUPRI,'(5D15.8)') 
     &        (E1(NX+IXX), IXX=1,NSTR(IFRP,0,1))
         NX = NX + NSTR(IFRP,0,1)
      ENDDO
      NX = 0
      DO IFRP = 1, NFSYM
         WRITE(LUPRI,'(A,I1)') 'Index 2, Irrep ',IFRP
         WRITE(LUPRI,'(5D15.8)') 
     &        (E2(NX+IXX), IXX=1,NSTR(IFRP,0,2))
         NX = NX + NSTR(IFRP,0,2)
      ENDDO
      NX = 0
      DO IFRP = 1, NFSYM
         WRITE(LUPRI,'(A,I1)') 'Index 3, Irrep ',IFRP
         WRITE(LUPRI,'(5D15.8)') 
     &        (E3(NX+IXX), IXX=1,NSTR(IFRP,0,3))
         NX = NX + NSTR(IFRP,0,3)
      ENDDO
      NX = 0
      DO IFRP = 1, NFSYM
         WRITE(LUPRI,'(A,I1)') 'Index 4, Irrep ',IFRP
         WRITE(LUPRI,'(5D15.8)') 
     &        (E4(NX+IXX), IXX=1,NSTR(IFRP,0,4))
         NX = NX + NSTR(IFRP,0,4)
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CONMP2 */
      SUBROUTINE CONMP2(EMP2,SPMAT,
     &                  E1,E2,E3,E4,N1,N2,N3,N4,
     &                  NSTR,ISPNR,JSPNR,IPRINT)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
#include "dgroup.h"
#include "dcbmpt.h"
      LOGICAL IBEQJB
      DIMENSION NSTR(2,0:2,4)
      DIMENSION E1(N1),E2(N2),E3(N3),E4(N4)
      DIMENSION SPMAT(N1,N2,N3,N4,16)
      DIMENSION ICLDAT(16)
      DATA ICLDAT/1,1,2,4,3,3,4,2,9,9,10,12,11,11,12,10/
C
      IF (IPRINT.LT.10) GOTO 123
      DO ICL = 1,4
         WRITE(LUPRI,'(/A,I5/)') 'CLASS:', ICL
         ICL2 = ICL + 4
         DO II = 1,N1
            DO IA = 1,N2
               DO IJ = 1,N3
                  DO IB = 1,N4
                     WRITE(LUPRI,'(4I4,2F20.10)')
     &                    II,IA,IJ,IB,SPMAT(II,IA,IJ,IB,ICL),
     &                    SPMAT(II,IA,IJ,IB,ICL2)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO ICL = 9,12
         WRITE(LUPRI,'(/A,I5/)') 'CLASS:', ICL
         ICL2 = ICL + 4
         DO II = 1,N1
            DO IA = 1,N2
               DO IJ = 1,N3
                  DO IB = 1,N4
                     WRITE(LUPRI,'(4I4,2F20.10)')
     &                    II,IA,IJ,IB,SPMAT(II,IA,IJ,IB,ICL),
     &                    SPMAT(II,IA,IJ,IB,ICL2)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
 123  CONTINUE
      IBEQJB = ISPNR.EQ.JSPNR
      IF (IPRINT.GE.2) THEN
         WRITE(LUPRI,'(//A/A/A/)') '====================',
     &        '  MP2 Eigenvalues', '===================='
         CALL PRIEIN(E1,E2,E3,E4,N1,N2,N3,N4,NSTR)
      ENDIF
      IF (TESTLS .OR .TESTSL) THEN
         WRITE(LUPRI,'(/A/A//A/A/A/A/)')
     &        '==============================================',
     &        '               W A R N I N G !',
     &        '  (LL|SS)- and (SS|LL)-ONLY tests should NOT',
     &        '    be run as a parallel job or loop over',
     &        '     more than a single batch of I and J!',
     &        '=============================================='
      ENDIF
C     
C     Real and complex groups, Class 1..8
C     
      EC = D0
      DO ICL = 0,3
         ICL1= ICL*2+1
         ICL2= ICL1+1
         IF (TESTLS .OR. TESTSL) THEN
            CALL ENERC1(E,E1,E2,E3,E4,N1,N2,N3,N4,
     &           SPMAT(1,1,1,1,ICLDAT(ICL1)),
     &           SPMAT(1,1,1,1,ICLDAT(ICL1)+4),
     &           SPMAT(1,1,1,1,ICLDAT(ICL2)),
     &           SPMAT(1,1,1,1,ICLDAT(ICL2)+4),IBEQJB,IPRINT)
         ELSE
            CALL ENERC2(E,E1,E2,E3,E4,N1,N2,N3,N4,
     &           SPMAT(1,1,1,1,ICLDAT(ICL1)),
     &           SPMAT(1,1,1,1,ICLDAT(ICL1)+4),
     &           SPMAT(1,1,1,1,ICLDAT(ICL2)),
     &           SPMAT(1,1,1,1,ICLDAT(ICL2)+4),IBEQJB,IPRINT)
         ENDIF
         EC=EC+E
      ENDDO
C
C     Quaternion groups, Class 9..16
C
      IF (NZ.EQ.4) THEN
         DO ICL = 4,7
            ICL1= ICL*2+1
            ICL2= ICL1+1
            IF (TESTLS .OR. TESTSL) THEN
               CALL ENERC1(E,E1,E2,E3,E4,N1,N2,N3,N4,
     &              SPMAT(1,1,1,1,ICLDAT(ICL1)),
     &              SPMAT(1,1,1,1,ICLDAT(ICL1)+4),
     &              SPMAT(1,1,1,1,ICLDAT(ICL2)),
     &              SPMAT(1,1,1,1,ICLDAT(ICL2)+4),IBEQJB,IPRINT)
            ELSE
               CALL ENERC2(E,E1,E2,E3,E4,N1,N2,N3,N4,
     &              SPMAT(1,1,1,1,ICLDAT(ICL1)),
     &              SPMAT(1,1,1,1,ICLDAT(ICL1)+4),
     &              SPMAT(1,1,1,1,ICLDAT(ICL2)),
     &              SPMAT(1,1,1,1,ICLDAT(ICL2)+4),IBEQJB,IPRINT)
            ENDIF
            EC=EC+E
         ENDDO
      ENDIF
      EMP2 = EMP2 + EC
      IF(IPRINT.GE.1) THEN
         WRITE(LUPRI,'(/A,D20.10)') 'Contribution, MP2', EC
         WRITE(LUPRI,'(A,D20.10/)') 'Accumulated, MP2 ', EMP2
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ENERC2 */
      SUBROUTINE ENERC2(E,EIG1,EIG2,EIG3,EIG4,N1,N2,N3,N4,
     & FMR1,FMI1,FMR2,FMI2,IBEQJB,IPRINT)
C
C     Written by Jon K. Laerdahl 1995
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      INTEGER I,A,J,B
      LOGICAL IBEQJB
      DIMENSION
     &     FMR1(N1,N2,N3,N4),FMI1(N1,N2,N3,N4),
     &     FMR2(N1,N2,N3,N4),FMI2(N1,N2,N3,N4)
      DIMENSION EIG1(N1),EIG2(N2),EIG3(N3),EIG4(N4)
C
C     IBEQJB is true if i and j are in same block
C     a and b should always be in same block!
C     -> NE2 = NE4 AND same range too! 
C
      E = 0.D0
      DO 100 I=1,N1
         DO 110 A=1,N2
            IF (IBEQJB) THEN
               JTOP = I
            ELSE
               JTOP = N3
            ENDIF
            DO 120 J=1,JTOP
               DO 130 B=1,A
c                  WRITE(LUPRI,'(/A,4I4)') 'aa', I,A,J,B
                  AR = FMR1(I,A,J,B)-FMR2(I,B,J,A)
                  AI = FMI1(I,A,J,B)-FMI2(I,B,J,A)
c                  WRITE(LUPRI,'(6F15.8)') 
c    &                 FMR1(I,A,J,B),FMR2(I,B,J,A),
c    &                 FMI1(I,A,J,B),FMI2(I,B,J,A),AR,AI
                  EN = AR**2 + AI**2
                  DIV = EIG1(I)+EIG3(J)-EIG2(A)
     &  -EIG4(B)
cc                  IF (DABS(EN).GT.1D-15)
c          WRITE(LUPRI,'(6F15.8)') 
c     & DIV,EIG1(I),EIG1(J),EIG2(A),EIG2(B),EN
                  DELIJ = 0.d00
                  IF ((I.EQ.J).AND.IBEQJB) DELIJ = 1.d00
                  IF (A.EQ.B) THEN
                     DELAB = 1.d00
                  ELSE
                     DELAB = 0.d00
                  ENDIF
                  E = E + (2.d00-DELIJ)*(2.d00-DELAB)*EN/(2.d00*DIV)
 130           CONTINUE
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
      IF(IPRINT.GE.1) 
     &     WRITE(LUPRI,'(A,D25.16)')  'MP2-energy, 1/4 matrix    :',E
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ENERC1 */
      SUBROUTINE ENERC1(E,EIG1,EIG2,EIG3,EIG4,N1,N2,N3,N4,
     &     FMR1,FMI1,FMR2,FMI2,IBEQJB,IPRINT)
C     
C     Written by Jon K. Laerdahl 1995
C     
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      INTEGER I,A,J,B
      LOGICAL IBEQJB
      DIMENSION
     &     FMR1(N1,N2,N3,N4),FMI1(N1,N2,N3,N4),
     &     FMR2(N1,N2,N3,N4),FMI2(N1,N2,N3,N4)
      DIMENSION EIG1(N1),EIG2(N2),EIG3(N3),EIG4(N4)
C     
      E = 0.D0
      DO 100 I=1,N1
         DO 110 A=1,N2
            DO 120 J=1,N3
               DO 130 B=1,N4
                  AR = FMR1(I,A,J,B)-FMR2(I,B,J,A)
                  AI = FMI1(I,A,J,B)-FMI2(I,B,J,A)
                  EN = AR**2 + AI**2
                  DIV = EIG1(I)+EIG3(J)-EIG2(A)-EIG4(B)
                  E = E + EN/(2.d00*DIV)
 130           CONTINUE
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
      IF(IPRINT.GE.1) 
     &     WRITE(LUPRI,'(A,D25.16)')  'MP2-energy, 1/4 matrix    :',E
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TSKARR */
      SUBROUTINE TSKARR(ISPNR,JSPNR,LENTSK,NSTRT,
     &               NSTR,NSTRIJ,ICMOQR,ICMOIJ,
     &               IPRINT)
C***********************************************************************
C
C     Set-up routine of defining arrays for given IJ-batch
C
C     Written by T.Saue Sep 9, based on previous code
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
C     Global variables
      DIMENSION NSTR(2,0:2,4),  ICMOQR(2,4)
      DIMENSION NSTRIJ(2,0:2,4),ICMOIJ(2,4)
      DIMENSION NSTRT(4)
C     Local variables
      DIMENSION NSTRS(2,4),NOCC(4)
C
C     Initialization
C
      CALL ICOPY( 8,-1,0,NSTRS,1)
      CALL ICOPY(24, 0,0,NSTRIJ,1)
C
      NOCC(1) = (ISPNR-1)*LENTSK + 1
      NOCC(2) = MIN(ISPNR*LENTSK,NSTRT(1))
      NOCC(3) = (JSPNR-1)*LENTSK + 1
      NOCC(4) = MIN(JSPNR*LENTSK,NSTRT(3))
      DO IN = 0,1
        IND1=IN*2+1
        IND2=IN*2+2
        LSTIR1 = NSTR(1,0,IND1)
C     
C       Whole batch in irrep 1
C     
        IF((NOCC(IND1).LE.LSTIR1).AND.
     &     (NOCC(IND2).LE.LSTIR1)) THEN
          NSTRS(1,IND1)    = NOCC(IND1)
          NSTRIJ(1,0,IND1) = NOCC(IND2)-NOCC(IND1) + 1
C     
C       Batch contributions from both irreps
C     
        ELSEIF((NOCC(IND1).LE.LSTIR1).AND.
     &         (NOCC(IND2).GT.LSTIR1)) THEN
          NSTRIJ(1,0,IND1) = LSTIR1-NOCC(IND1)+1
          NSTRIJ(2,0,IND1) = NOCC(IND2)-LSTIR1
          NSTRS(1,IND1)    = NOCC(IND1)
          NSTRS(2,IND1)    = 1
C     
C       Whole batch in irrep 2
C         
        ELSE
          NSTRS(2,IND1)    = NOCC(IND1) - LSTIR1
          NSTRIJ(2,0,IND1) = NOCC(IND2) - NOCC(IND1) + 1
        ENDIF
        NSTRIJ(1,1,IND1) = NSTRIJ(1,0,IND1)
        NSTRIJ(1,2,IND1) = 0
        NSTRIJ(2,1,IND1) = NSTRIJ(2,0,IND1)
        NSTRIJ(2,2,IND1) = 0
      ENDDO
C     
      CALL ICOPY (6,NSTR(1,0,2),1,NSTRIJ(1,0,2),1)
      CALL ICOPY (6,NSTR(1,0,4),1,NSTRIJ(1,0,4),1)
      NSTRS(1,2) = 1
      NSTRS(1,4) = 1
      NSTRS(2,2) = 1
      NSTRS(2,4) = 1
C     
C     Give relative offset for the coefficients
C     
      DO I = 1, 4
        DO IFRP = 1, NFSYM
          ICMOIJ(IFRP,I) = ICMOQR(IFRP,I)
     &              + (NSTRS(IFRP,I)-1)*NFBAS(IFRP,0)
        ENDDO
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.1) THEN
        WRITE(LUPRI,'(/A)') 'Task-info 2'
        WRITE(LUPRI,'(A,3I5)') 'I, J, LENTSK: ',ISPNR,JSPNR,LENTSK
        WRITE(LUPRI,'(A,2I5)') 'Occupied index I (first/last) :',
     &        NOCC(1),NOCC(2)
        WRITE(LUPRI,'(A,2I5)') 'Occupied index J (first/last) :',
     &        NOCC(3),NOCC(4)
        DO IFRP = 1,NFSYM
          WRITE(LUPRI,'(A,A3)') '** Fermion ircop ',FREP(IFRP)
          DO K = 1,4
            WRITE(LUPRI,'(3X,A,I3,2(A,I5))')
     &       '* Index:',K,' - offset :',NSTRS(IFRP,K),
     &                 ' - vectors:',NSTRIJ(IFRP,0,K)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MP2HDR */
      SUBROUTINE MP2HDR(LUNIT,MXBTCH,MXTASK,NODES)
#include "implicit.h"
      WRITE(LUNIT,'(58A1/A,12X,A,13X,A/58A1/A/A)')
     &   ('=', I=1,58),
     &   '===', 'Progress of MP2 Calculation','===',
     &   ('=', I=1,58),
     &   ' Current   Tasks    Tasks    Slaves',
     &   '  batch    sent    received released         Time'
      IF (NODES .LT. 0) THEN
         WRITE(LUNIT,'(2(A,I4,A,3X)/58A1)')
     &        ' [', MXBTCH, ']', '[', MXTASK, ']',('-', I=1,58)
      ELSE
         WRITE(LUNIT,'(4(A,I4,A,3X)/58A1)')
     &        ' [', MXBTCH, ']', '[', MXTASK, ']', '[', MXTASK, ']',
     &        '[', NODES, ']',
     &        ('-', I=1,58)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MKNSTT */
      SUBROUTINE MKNSTT(NST,NSTT)
#include "implicit.h"
#include "dgroup.h"
      DIMENSION NST(2,0:2,4), NSTT(4)
      DO I = 1, 4
        NSTT(I) = 0
        DO IFRP = 1, NFSYM
          NSTT(I) = NSTT(I) + NST(IFRP,0,I)
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PCKINI */
      SUBROUTINE PCKINI(NST,ITQMT,NTQMT,N4QMT,IPRINT)
C***********************************************************************
C
C     Setup packing of integrals and determine memory for
C     3/4 and 4/4 transformed matrix
C
C***********************************************************************
#include "implicit.h"
#include "dcbtra.h"
#include "dgroup.h"
      DIMENSION NST(2,0:2,4), ITQMT(0:7,0:7,2)
      LOGICAL TRIAN(2)
C
      TRIAN(1) = ISAME(2) .EQ. ISAME(1)
      TRIAN(2) = .FALSE.
      CALL PCK2IN(NST,TRIAN,IPRINT)
      CALL PCK3IN(ITQMT,NTQMT,NST,IPRINT)
      N4QMT = 0
      DO IREPIJ = 1, NFSYM
         DO IREPRS = 0, NBSYM-1
            IF(IREPIJ.EQ.JBTOF(IREPRS,1)) THEN
               NIJ = NFPCK12(IREPIJ) ! Number of ij pairs 
               NKL = NFPCK34(IREPIJ) ! Number of kl pairs 
               N4QMT = N4QMT + NIJ*NKL*NZ**2
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck T44INI */
      SUBROUTINE T44INI(INTFLG,ICS2,ICF2)
C***********************************************************************
C
C     Setup routine for 4/4 transformation
C
C***********************************************************************
#include "implicit.h"
#include "dcbmpt.h"
      LOGICAL L2TYP(0:3)
C
      ICS2 = 1
      ICF2 = 2
      L2TYP(1) = MOD(INTFLG,2) .EQ. 1
      L2TYP(2) = MOD((INTFLG/2),2) .EQ. 1
      L2TYP(3) = MOD(((INTFLG/2)/2),2) .EQ. 1
      IF ((.NOT. L2TYP(1)) .AND. (.NOT. L2TYP(2))) ICS2 = 2
      IF ((.NOT. L2TYP(2)) .AND. (.NOT. L2TYP(3))) ICF2 = 1
      IF (TESTLS) ICS2 = 2
      IF (TESTSL) ICF2 = 1
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck LSZINI */
      SUBROUTINE LSZINI(NSTT,LSIZ)
C***********************************************************************
C
C     Determine size of final anti-symmetrized integrals
C
C***********************************************************************
#include "implicit.h"
#include "dgroup.h"
      DIMENSION NSTT(4)
C
      IF (NZ .EQ. 1) THEN
         NCLASS = 4
      ELSE IF (NZ .EQ. 2) THEN
         NCLASS = 8
      ELSE
         NCLASS = 16
      ENDIF
      CALL DEFCL(NCLASS)
      LSIZ = 16
      IF (NZ .LT. 4) LSIZ = 8
      DO I = 1, 4
         LSIZ = LSIZ*NSTT(I)
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C Start of parallel section : activitate compilation only when necessary
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TR2MST */
      SUBROUTINE TR2MST(ERGMP2,WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,
     &     KINDX,KGAB,KDRIJ,INTFLG,ICMOIJ,NSTR,NSTRT,NSTRIJ,NSTTIJ,
     &     ANTIS,LMP2,NDSLST,IDSLST,TIMDIS,DINTSKP)
C***********************************************************************
C
C     Master driver for new parallel MP2 scheme.
C
C     Vebjorn Bakken, Summer 2003
C
C***********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
#include "maxorb.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "infpar.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbmp2.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION KQ(2,4), KE(2,2), KIBE(2,4), ICMOIJ(2,4),
     &     NSTR(2,0:2,4), NSTRT(4), NSTRIJ(2,0:2,4), NSTTIJ(4),
     &     WORK(*), IDSLST(6,NDSLST), TIMDIS(NDSLST),DINTSKP(8),
     &     ITQMT(0:7,0:7,2), IDSNFO(6), IRELEA(6), IBTNFO(2)
      CHARACTER*(24) TMSTMP
      DATA IRELEA /-1, 0, 0, 0, 0, 0/
      LOGICAL ANTIS, LMP2, TRIAN(2), MP2FIN
C
      CALL QENTER('TR2MST')
C
C     Initialize
C        Get hold of the slaves and initialize them
C        (ITASK = 2 for parallel integral transformation)
C     
      CALL DIRAC_PARCTL( MOLTRA_PAR )
      KFRSAV = KFREE
      TMMP2 = 0.D0
      CALL MKNSTT(NSTR,NSTRT)
      CALL TRAPARI(WORK,KFREE,LFREE,KQ,KE,KIBE,
     &     KINDX,KGAB,KDRIJ,INTFLG,NSTR,ANTIS,LMP2)
C     
C     Task set-up
C        Set up tasks with different batches of spinors if requested.
C        In this parallel scheme the number of tasks is not related
C        to the number of nodes in any way.
C     
      CALL TSKCHK(MAXTSK,LENTSK,LFREE,1,NSTRT(1),IPRINT)
C ======================================================================
C     Outer loop over batches of occupied spinors
C ======================================================================
      DO IJBTCH = 1, MAXTSK
         CALL UNPKIJ(IJBTCH,ISPNR,JSPNR)
         CALL TSKARR(ISPNR,JSPNR,LENTSK,NSTRT,
     &        NSTR,NSTRIJ,ICMOQR,ICMOIJ,IPRINT)
         CALL MKNSTT(NSTRIJ,NSTTIJ)
C
C     Allocation
C        Determine and allocate memory for 4/4 transformed matrix,
C        WORK(K4QMAT), and accumulation matrix, WORK(K4QACC). The former
C        will be used to receive integrals from the slaves, the latter
C        to accumulate them.
C
         CALL PCKINI(NSTRIJ,ITQMT,NTQMT,N4QMT,IPRINT)
         IF (IPRINT .GE. 1) WRITE (LUPRI,1011) 2*N4QMT*8./(1024.*1024.)
         CALL MEMGET('REAL',K4QMAT,N4QMT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',K4QACC,N4QMT,WORK,KFREE,LFREE)
         CALL DZERO(WORK(K4QACC),N4QMT)
C
C     Tasks set-up
C        Set up the shells to be calculated and send out
C        batch information (batch number and length of batch)
C
         CALL DZERO(TIMDIS,NDSLST)
         CALL SHLLST(NDSLST,IDSLST,INTFLG,WORK(KINDX),
     &               LFREE,IPRINT,SRTSHL)
         IF ((IPRINT .EQ. 0) .AND. (IJBTCH .EQ. 1))
     &        CALL MP2HDR(LUPRI,MAXTSK,NDSLST,NODES)
         IBTNFO(1) = IJBTCH
         IBTNFO(2) = LENTSK
         call interface_mpi_BCAST(IBTNFO,2,MPARID,
     &                  global_communicator)
         IF (CMMNFO) WRITE(LUPRI,'(/A,I4,/A,I4)')
     &      ' <COMMINFO> Master broadcasted batch number: ', IBTNFO(1),
     &      ' <COMMINFO>             and length of batch: ', IBTNFO(2)
C
C     Loop over tasks until all tasks have been sent and received,
C     and all slaves have been released (reset for each batch)
C
         ITSKSD = 0
         ITSKRV = 0
         ISLVRL = 0
C ----------------------------------------------------------------------
C     Inner loop over tasks
C ----------------------------------------------------------------------
 100     CONTINUE
C
C     Master is idle until it receives a signal from one of the slaves
C
         IF (CMMNFO) WRITE(LUPRI,'(A/A)') ' <COMMINFO>',
     &        ' <COMMINFO> Master is listening...'
         CALL FLSHFO(LUPRI)
#if defined (MPE)
         CALL MPE_LOG_EVENT(3, 0, "start wait")
         call interface_mpi_RECV(NODINT,1,df_MPI_ANY_SOURCE,90,
     &        global_communicator,ISTAT)
         CALL MPE_LOG_EVENT(4, 0, "end wait")
#else
         call interface_mpi_RECV(NODINT,1,df_MPI_ANY_SOURCE,90,
     &        global_communicator,ISTAT)
#endif
         NODE = ISTAT(df_MPI_SOURCE)
         IF (CMMNFO) WRITE(LUPRI,'(2(A,I4))')
     &        ' <COMMINFO> Master received ', NODINT,
     &        ' from node #', NODE
C
C     Slave is idle and there's more tasks available:
C        Send next task to slave
C
         IF ((NODINT .EQ. 0) .AND. (ITSKSD .LT. NDSLST)) THEN
            ITSKSD = ITSKSD + 1
            DO I = 1, 6
               IDSNFO(I) = IDSLST(I,ITSKSD)
            END DO
#if defined (MPE)
            CALL MPE_LOG_EVENT(5, 0, "start sndtsk")
            call interface_mpi_SEND(IDSNFO,6,NODE,92,
     &           global_communicator)
            CALL MPE_LOG_EVENT(6, 0, "end sndtsk")
#else
            call interface_mpi_SEND(IDSNFO,6,NODE,92,
     &           global_communicator)
#endif
            IF (CMMNFO) WRITE(LUPRI,*) '<COMMINFO> Master sent task ',
     &           ITSKSD, ' to node #', NODE
C
C     Slave is idle but there are no more tasks available:
C        Release slave by sending negative shell number
C
         ELSE IF ((NODINT .EQ. 0) .AND. (ITSKSD .GE. NDSLST)) THEN
            call interface_mpi_SEND(IRELEA,6,NODE,92,
     &           global_communicator)
            IF (CMMNFO) WRITE(LUPRI,'(A,I4)')
     &           ' <COMMINFO> Master released node #', NODE
            ISLVRL = ISLVRL + 1
C
C     Slave has integrals to send back
C
         ELSE
#if defined (MPE)
            CALL MPE_LOG_EVENT(7, 0, "start rcvint")
            call interface_mpi_RECV(WORK(K4QMAT),N4QMT,
     &           NODE,94,global_communicator)
            CALL MPE_LOG_EVENT(8, 0, "end rcvint")
#else
            call interface_mpi_RECV(WORK(K4QMAT),N4QMT,
     &           NODE,94,global_communicator)
#endif
            IF (CMMNFO) WRITE(LUPRI,'(2A,I4)')
     &           ' <COMMINFO> Master received ',
     &           'fully transformed integrals from node #', NODE
Csknecht: see if this can be replaced by MPI_REDUCE and MPI_IN_PLACE
            DO I = 1, N4QMT
               WORK(K4QACC+I-1) = WORK(K4QACC+I-1) + WORK(K4QMAT+I-1)
            END DO
            ITSKRV = ITSKRV + NODINT
         END IF
C
C     Print out information to indicate progress
C
         IF (CMMNFO) THEN
            WRITE(LUPRI,*) '<COMMINFO>'
            WRITE(LUPRI,*) '<COMMINFO> ### STATUS ### '
            IF (MAXTSK .GT. 1)
     &           WRITE(LUPRI,*) '<COMMINFO> Batch number   : ',
     &           IJBTCH,   '/', MAXTSK
            WRITE(LUPRI,*) '<COMMINFO> Tasks sent     : ',
     &           ITSKSD,   '/', NDSLST
            WRITE(LUPRI,*) '<COMMINFO> Tasks received : ',
     &           ITSKRV, '/', NDSLST
            WRITE(LUPRI,*) '<COMMINFO> Slaves released: ',
     &           ISLVRL, '/', NODES
            CALL FLSHFO(LUPRI)
         ELSE
            IF (IPRINT .GE. 1) THEN
               IF (MAXTSK .GT. 1) THEN
                  WRITE(LUPRI,'(4(A,I3,A,I3),A)')
     &              'Batch ', IJBTCH, '/', MAXTSK,
     &              ' -=- Sent ', ITSKSD, '/', NDSLST,
     &              ' tasks -=- Received ', ITSKRV, '/', NDSLST,
     &              ' tasks -=- Released ', ISLVRL, '/', NODES,' slaves'
               ELSE
                  WRITE(LUPRI,'(3(A,I3,A,I3),A)')
     &              'Sent ', ITSKSD, '/', NDSLST,
     &              ' tasks -=- Received ', ITSKRV, '/', NDSLST,
     &              ' tasks -=- Released ', ISLVRL, '/', NODES,' slaves'
               END IF
            ELSE IF (IPRINT .EQ. 0) THEN
               CALL GTINFO(TMSTMP)
               WRITE(LUPRI,'(I6,3I9,5X,A19)') IJBTCH, ITSKSD, ITSKRV,
     &              ISLVRL, TMSTMP(1:19)
            END IF
            CALL FLSHFO(LUPRI)
         END IF
C
C     Check if current batch is finished
C
         IF ((ITSKSD .LT. NDSLST) .OR. (ISLVRL .LT. NODES)
     &        .OR. (ITSKRV .LT. NDSLST)) GOTO 100
C ----------------------------------------------------------------------
C     End inner loop over tasks
C ----------------------------------------------------------------------
         IF (IPRINT .EQ. 0)
     &        WRITE(LUPRI,'(58A1)') ('-', I=1,58)
C
C     Print accumulated fully transformed integrals if requested
C     (above threshold TPRI44)
C
         IF (TPRI44 .GT. 0.D0)
     &        CALL PRIINT('The accumulated fully transformed integrals',
     &        N4QMT, WORK(K4QACC), TPRI44, 0)
C
C     Calculate contribution to MP2 energy for this complete batch
C
#if defined (MPE)
         CALL MPE_LOG_EVENT(9, 0, "start calc")
#endif
         CALL TIMEC('START ',TIMSTR,TIMEND)
         CALL LSZINI(NSTTIJ,LSIZ)
C
C     Allocate memory and perform naive anti-symmetrize of final integrals
C
         CALL MEMGET('REAL',KANT,LSIZ,WORK,KFREE,LFREE)
         IF (IPRINT .GE. 1) WRITE (LUPRI,1021) LSIZ*8./(1024.*1024.)
         CALL DZERO(WORK(KANT),LSIZ)
         CALL NANTSY(WORK,KFREE,LFREE,IPRINT,NSTRIJ,NSTTIJ,NCLASS,LSIZ,
     &        K4QACC,KANT)
C
C     Calculate MP2 contribution
C
         I1 = KE(1,1) + (ISPNR-1)*LENTSK
         I2 = KE(1,2)
         I3 = KE(1,1) + (JSPNR-1)*LENTSK
         I4 = KE(1,2)
         CALL CONMP2(ERGMP2,WORK(KANT),
     &        WORK(I1),WORK(I2),WORK(I3),WORK(I4),
     &        NSTTIJ(1),NSTTIJ(2),NSTTIJ(3),NSTTIJ(4),
     &        NSTRIJ,ISPNR,JSPNR,IPRINT)
         CALL TIMEC('CONMP2',TIMSTR,TIMEND) 
         TMMP2 = TMMP2 + TIMEND - TIMSTR
#if defined (MPE)
         CALL MPE_LOG_EVENT(10, 0, "end calc")
#endif
CTROND
      WRITE(6,*) 'After2 CONMP2...'
C
C     Release memory after current batch
C
         CALL MEMREL('TR2MST',WORK,1,KFRSAV,KFREE,LFREE)
C
C     Slaves want to know if there are more batches
C
         MP2FIN = .FALSE.
         IF (IJBTCH .EQ. MAXTSK) MP2FIN = .TRUE.
         call interface_mpi_bcast_l0(MP2FIN,1,MPARID,
     &                                   global_communicator)
         IF (CMMNFO) WRITE(LUPRI,'(A,L5)')
     &        ' <COMMINFO> Master broadcasted MP2FIN =', MP2FIN
      END DO
C ======================================================================
C     End outer loop over batches
C ======================================================================

C
C     Collect timing information from all slaves before releasing them
C
      KFRSAV = KFREE
      NTIM = 7*(NODES+1)
      CALL MEMGET('REAL',KTIM,NTIM,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KTIM),NTIM)
      CALL MEMGET('REAL',KBUF,8,WORK,KFREE,LFREE)
CTROND
      WRITE(6,*) 'After3 CONMP2...',NODES
      DO INOD = 1,NODES
         call interface_mpi_RECV(WORK(KTIM+INOD*7),7,
     &        INOD,96,global_communicator)
         IF (CMMNFO) WRITE(LUPRI,*) '<COMMINFO> Master received '//
     &        'timing information from node #', INOD
          call interface_mpi_RECV(WORK(KBUF),8,
     &          INOD,40,global_communicator)
          CALL DAXPY(8,D1,WORK(KBUF),1,DINTSKP,1)
      END DO
CTROND
      WRITE(6,*) 'After4 CONMP2...'
      CALL DIRAC_PARCTL( EXIT_NODEMENU )
#if defined (VAR_PFS)
C
C       close global file (opened in TRAPARI)
C
        CALL GLOBAL_FILE_CLOSE_READ
#endif
C
C     Print timing info
C
      WRITE(LUPRI,'(/A)') '** Master:'
      CALL TIMTXT('Time used in MP2CAL is: ',TMMP2,LUPRI)
      DO INOD = 1,NODES
         WRITE(LUPRI,'(A,I4)') '** Node :',INOD
         CALL TIMTXT('Time used in CALDIS is: ',
     &                                      WORK(KTIM+INOD*7), LUPRI)
         CALL TIMTXT('Time used in MS4IN1 is: ',
     &                                      WORK(KTIM+INOD*7+1), LUPRI)
         CALL TIMTXT('Time used in TQTRNS is: ',
     &                                      WORK(KTIM+INOD*7+2), LUPRI)
         CALL TIMTXT('Time used in LQTRNS is: ',
     &                                      WORK(KTIM+INOD*7+3), LUPRI)
      END DO
C
CTROND
      WRITE(6,*) 'After5 CONMP2...'
      CALL MEMREL('TR2MST',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TR2MST')
      RETURN
C
 1011 FORMAT (/
     &     ' Master: Storing all fully-transformed integrals to memory',
     &     /' Memory requirements ',F10.3,' Megabytes')
 1021 FORMAT (/' Master: Storing all final integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TR2NOD */
      SUBROUTINE TR2NOD(WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,
     &     KINDX,GABRAO,DRIJ,INTFLG,NSTR,ANTIS,LMP2,
     &     DINTSKP,TIM1)
C***********************************************************************
C
C     Slave driver for new parallel MP2 scheme.
C
C     Vebjorn Bakken, Summer 2003
C
C***********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "maxorb.h"
#include "infpar.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbmp2.h"
#include "dcbbas.h"
#include "dgroup.h"
      LOGICAL ANTIS, LMP2, L2TYP(0:3), TRIAN(2), MP2FIN
      DIMENSION KQ(2,4), KE(2,2), NSTR(2,0:2,4), NSTRT(4), KIBE(2,4),
     &     WORK(*), TIM1(7), DINTSKP(8), IDSNFO(6), IBTNFO(2),
     &     ITQMT(0:7,0:7,2),DRIJ(*),GABRAO(*)
      DIMENSION NSTRIJ(2,0:2,4), NSTTIJ(4), ICMOIJ(2,4)
      DIMENSION ISAMEIJ(4)
C
C     Initialization
C
      KFRSAV = KFREE
      CALL MKNSTT(NSTR,NSTRT)
      CALL DZERO(TIM1,7)
C ======================================================================
C     Loop until told to stop by master
C     There's an outer loop over batches of occupied spinors,
C     and an inner loop over tasks (shells) for a given IJ-batch
C ======================================================================
 100  CONTINUE
      call interface_mpi_BCAST(IBTNFO,2,MPARID,global_communicator)
      IJBTCH = IBTNFO(1)
      LENTSK = IBTNFO(2)
      IF (CMMNFO) WRITE(LUPRI,'(3(A,I4),A)')
     &     ' <COMMINFO> Node #', MYTID,' received batchinfo (', 
     &     IJBTCH,', ', LENTSK, ') from master'
C
C     Initialize for this batch
C
      INTGTH = 0
      CALL UNPKIJ(IJBTCH,ISPNR,JSPNR)
      CALL TSKARR(ISPNR,JSPNR,LENTSK,NSTRT,
     &     NSTR,NSTRIJ,ICMOQR,ICMOIJ,IPRINT)
      CALL MKNSTT(NSTRIJ,NSTTIJ)
C
C     Define the packing of the 2-index and 3-index transformed integrals
C     and determine memory for 3/4 and 4/4 transformation matrices
C
      CALL PCKINI(NSTRIJ,ITQMT,NTQMT,N4QMT,IPRINT)
C
C     Get screening matrices
C
      IF (SCRTRA .GT. D0) THEN
         DO I = 1, 4
            ISAMEIJ(I) = I
         ENDDO
         CALL PR4SC2(KDRIJ,WORK,KFREE,LFREE,GABRAO,
     &        WORK(KQ(1,1)),WORK(KQ(1,2)),
     &        WORK(KQ(1,3)),WORK(KQ(1,4)),
     &        NDMOQR,ICMOIJ,NSTRIJ,ISAMEIJ,IPRINT)
      ELSE
        KDRIJ = KFREE
      ENDIF
C
C     Allocate memory for accumulation of 3/4 transformed integrals
C
      IF (IPRINT .GE. 1) WRITE (LUPRI,1002) NTQMT*8./(1024.*1024.)
      IF (LFREE.LT.NTQMT) THEN
         WRITE(LUPRI,'(/A,I4,A//A,/A,/A)') 
     &       '!!! ERROR on node ', MYTID, ' !!!',
     &       '### Not enough memory for 3-index transf. matrix.   ###',
     &       '### Reduce memory requirements by setting           ###',
     &       '### .IJTSK to a lower value!                        ###'
         CALL QUIT('TR2NOD: Not enough memory for 3-index matrix')
      ENDIF
      CALL MEMGET('REAL',KTQMAT,NTQMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KTQMAT),NTQMT)
C ----------------------------------------------------------------------
C     Inner loop over tasks
C ----------------------------------------------------------------------
 120  CONTINUE
C
C     Sending 0 to master signals that slave is idle
C
      IF (CMMNFO) WRITE(LUPRI,'(A,I4)')
     &     ' <COMMINFO> Sending 0 from node #', MYTID
#if defined (MPE)
      CALL MPE_LOG_EVENT(3, 0, "start wait")
      call interface_mpi_SEND(0,1,MPARID,90,global_communicator)
      CALL MPE_LOG_EVENT(4, 0, "end wait")
C
      CALL MPE_LOG_EVENT(5, 0, "start rcvtsk")
      call interface_mpi_RECV(IDSNFO,6,MPARID,92,
     &     global_communicator)
      CALL MPE_LOG_EVENT(6, 0, "end rcvtsk")
#else
      call interface_mpi_SEND(0,1,MPARID,90,global_communicator)
      call interface_mpi_RECV(IDSNFO,6,MPARID,92,
     &     global_communicator)
#endif
C
C     IDSNFO(1) < 0 signals no more available tasks --> exit inner loop
C
      IF (IDSNFO(1) .LT. 0) GOTO 200
      IF (CMMNFO) WRITE(LUPRI,'(A,I4,A,5I4,I8)') ' <COMMINFO> Node #',
     &     MYTID, ' received task:', IDSNFO
C
C     Update 3/4 transformed matrix with contribution from current shell
C
#if defined (MPE)
      CALL MPE_LOG_EVENT(9, 0, "start calc")
#endif
      CALL SHLTRA(IDSNFO,WORK,KFREE,LFREE,IPRINT,KQ,
     &        WORK(KINDX),INTFLG,NSTRIJ,NSTTIJ,DINTSKP,GABRAO,
     &        WORK(KDRIJ),NDMOQR,ICMOIJ,ITQMT,NTQMT,KTQMAT,TIM1)
      INTGTH = INTGTH + 1
C
C     Print integrals if requested (above threshold TPRI34):
C
      IF (TPRI34 .GT. 0.D0)
     &     CALL PRIINT('The accumulated 3/4 transformed integrals',
     &     NTQMT, WORK(KTQMAT), TPRI34, MYTID)
#if defined (MPE)
      CALL MPE_LOG_EVENT(10, 0, "end calc")
#endif
      GOTO 120
C ----------------------------------------------------------------------
C     End inner loop over tasks
C ----------------------------------------------------------------------
 200  CONTINUE
C     
C     4/4 transformation
C
      IF (IPRINT .GE. 1) WRITE (LUPRI,1012) N4QMT*8./(1024.*1024.)
      CALL MEMGET('REAL',K4QMAT,N4QMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(K4QMAT),N4QMT)
      CALL T44INI(INTFLG,ICS2,ICF2)
      CALL TIMEC('START ',TIMSTR,TIMEND)
#if defined (MPE)
      CALL MPE_LOG_EVENT(9, 0, "start calc")
#endif
      CALL LQTRNS(WORK,KFREE,LFREE,IPRINT,
     &     ITQMT,ICMOQR(1,4),NDMOQR(1,1,4),
     &     WORK(KTQMAT),NSTRIJ(1,0,3),NSTRIJ(1,0,4),
     &     WORK(KQ(1,4)),WORK(K4QMAT),ICS2,ICF2)
#if defined (MPE)
      CALL MPE_LOG_EVENT(10, 0, "end calc")
#endif
      CALL TIMEC('LQTRNS',TIMSTR,TIMEND) 
      TIM1(4) = TIM1(4) + TIMEND - TIMSTR
C
C     Print fully transformed integrals if requested
C     (above threshold TPRI44):
C
      IF (TPRI44 .GT. 0.D0)
     &     CALL PRIINT('The fully transformed integrals',
     &     N4QMT, WORK(K4QMAT), TPRI44, MYTID)
C
C     Send integrals as soon as master is ready to accept them
C
      IF (CMMNFO) WRITE(LUPRI,'(2(A,I4))') ' <COMMINFO> Sending ',
     &     INTGTH, ' from node #', MYTID
#if defined (MPE)
      CALL MPE_LOG_EVENT(3, 0, "start wait")
      call interface_mpi_SEND(INTGTH,1,
     &     MPARID,90,global_communicator)
      CALL MPE_LOG_EVENT(4, 0, "end wait")
C
      CALL MPE_LOG_EVENT(7, 0, "start sndint")
      call interface_mpi_SEND(WORK(K4QMAT),N4QMT,
     &     MPARID,94,global_communicator)
      CALL MPE_LOG_EVENT(8, 0, "end sndint")
#else
      call interface_mpi_SEND(INTGTH,1,
     &     MPARID,90,global_communicator)
      call interface_mpi_SEND(WORK(K4QMAT),N4QMT,
     &     MPARID,94,global_communicator)
#endif
      IF (CMMNFO) WRITE(LUPRI,'(A,I4,A,I4)') ' <COMMINFO> Node #',
     &     MYTID, ' sent fully transformed integrals for batch ', IJBTCH
C
C     If master sends MP2FIN = .TRUE., there are no more batches
C
      call interface_mpi_bcast_l0(MP2FIN,1,MPARID,
     &                                global_communicator)
      IF (CMMNFO) WRITE(LUPRI,'(A,I4,A,L5,A)') ' <COMMINFO> Node #',
     &     MYTID, ' received MP2FIN =', MP2FIN, ' from master'
      IF (.NOT. MP2FIN) THEN
         CALL MEMREL('TR2NOD',WORK,1,KFRSAV,KFREE,LFREE)
         GOTO 100
      END IF
C ======================================================================
C     End outer loop over batches
C ======================================================================

C
C     Send timing information when all batches are done, then exit
C
      IF (CMMNFO) WRITE(LUPRI,*)
     &     '<COMMINFO> Sending timing information from node #', MYTID
      call interface_mpi_SEND(TIM1,7,
     &     MPARID,96,global_communicator)
      call interface_mpi_SEND(DINTSKP,8,
     &          MPARID,40,global_communicator)
      RETURN
C     
 1002 FORMAT (/'Node #', I4, 
     &     ': Storing all 3-index transformed integrals to memory',
     &     /' Memory requirements ',F10.3,' Megabytes')
 1012 FORMAT (/'Node #', I4, 
     &     ': Storing all fully-transformed integrals to memory',
     &     /' Memory requirements ',F10.3,' Megabytes')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TR2PAR */
      SUBROUTINE TR2PAR(WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,
     &        KINDX,GABRAO,DRIJ,INTFLG,NSTR,ANTIS,LMP2,
     &        EMP2,DINTSKP,TIM1)
C***********************************************************************
C
C     Slave driver for scheme 2.
C
C     Written by T.Saue Sep 17 1998
C
C***********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "maxorb.h"
#include "infpar.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
#include "dgroup.h"
      LOGICAL ANTIS,LMP2,TRIAN(2)
      DIMENSION KQ(2,4),KE(2,2),NSTR(2,0:2,4),NSTRT(4),KIBE(2,4)
      DIMENSION NSTRIJ(2,0:2,4),ICMOIJ(2,4)
      DIMENSION WORK(*),TIM1(7),DINTSKP(8),GABRAO(*),DRIJ(*)
C
C     Initialization
C
      IF (LMP2) THEN
         EMP2 = D0
         CALL DZERO(TIM1,7)
      ENDIF
      DO IND = 1,4
        NSTRT(IND) = 0
        DO IFRP = 1,NFSYM
          NSTRT(IND) = NSTRT(IND) + NSTR(IFRP,0,IND)
        ENDDO
      ENDDO
C
C     Wait for tasks until we receive the end message
C
  100 CONTINUE
         CALL SLEEP(1)
         ITEST = 0
         call interface_mpi_SEND(ITEST,1,MPARID,20,
     &                 global_communicator)
C        ITEST .eq. 0: code for "send me new task"
#if defined (MPE)
         CALL MPE_LOG_EVENT(5, 0, "start sndtsk")
#endif
         call interface_mpi_RECV(IJTSK,1,MPARID,30,
     &                 global_communicator)
#if defined (MPE)
         CALL MPE_LOG_EVENT(6, 0, "end sndtsk")
#endif
         IF (IJTSK.EQ.-1) GOTO 200
         call interface_mpi_RECV(LENTSK,1,MPARID,30,
     &                 global_communicator)
         CALL UNPKIJ(IJTSK,ISPNR,JSPNR)
         CALL TSKARR(ISPNR,JSPNR,LENTSK,NSTRT,
     &               NSTR,NSTRIJ,ICMOQR,ICMOIJ,
     &               IPRINT)
#if defined (MPE)
         CALL MPE_LOG_EVENT(9, 0, "start calc")
#endif
         CALL TRDR23(WORK,KFREE,LFREE,IPRINT,
     &        EMP2,ISPNR,JSPNR,LENTSK,
     &        INTFLG,NDMOQR,ICMOIJ,NSTRIJ,
     &        ANTIS,LMP2,WORK(KINDX),KQ,KE,
     &        DINTSKP,WORK(KGAB),TIM1)
#if defined (MPE)
         CALL MPE_LOG_EVENT(10, 0, "end calc")
#endif
         GOTO 100
  200 CONTINUE
C
C     send whatever results we have
C
#if defined (MPE)
      CALL MPE_LOG_EVENT(7, 0, "start sndint")
#endif
      CALL TRARES(WORK,KFREE,LFREE,DINTSKP,LMP2,EMP2,TIM1,DUM)
#if defined (MPE)
      CALL MPE_LOG_EVENT(8, 0, "end sndint")
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C& End of parallel section : activitate compilation only when necessary
#endif
