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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RDINTS(IPRNT_LOC,NSP,KR,INDK,INDL,CBUF,IBUF,NVT,
     &                  BUF1,BUF2,
     &                  VOOOO,VOOOOA,VVOOO,VVOOOA,VVVOO,VVVOOA,
     &                  VVOVO,VVOVOA,VVOVV,VVOVVA,VVVVV,VVVVVA,
     &                  EPS,FOO,FVO,FVV,BREIT_LOC,SRTALL,INTERFACE_LOC)
C
      use interface_to_mpi
      use spinor_indexing
      implicit none
C
C---------------Description--------------------------------------------
C
C     Reads 1-electron and 2-electron integrals
C     Constructs Fock matrix elements
C     This version reads from files created by TMOONE & ROTRAN or DIRAC
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "complex.inc"
#include "results.inc"
#include "files.inc"
#include "ccpar.inc"
#include "inpt.inc"
#include "waio.h"
#include "dcbgen.h"
C
C---------------Calling variables--------------------------------------
C
      INTEGER*8 NVT
      INTEGER NSP,iprnt_loc
      INTEGER KR(-NSP:NSP),INDK(*),INDL(*)
      INTEGER IBUF(7,*)
      REAL*8 FOO(*),FVO(*),FVV(*)
      REAL*8 VOOOO(RCW,*),VOOOOA(RCW,*),VVOOO(RCW,*),VVOOOA(RCW,*)
      REAL*8 VVVOO(RCW,*),VVVOOA(RCW,*),VVOVO(RCW,*),VVOVOA(RCW,*)
      REAL*8 VVOVV(RCW,*),VVOVVA(RCW,*),VVVVV(RCW,*),VVVVVA(RCW,*)
      REAL*8 BUF1(*),BUF2(*),CBUF(RCW,*)
      REAL*8 EPS(*), ORIG_TOTSCF
      LOGICAL BREIT_LOC,SRTALL,TOBE
      CHARACTER*10 INTERFACE_LOC
C
C---------------Local variables--------------------------------------
C
      REAL*8 E0,E1
      CHARACTER*10 DATEX,TIMEX*8
      LOGICAL FOCKSP,PERMSYM
      INTEGER*8 N5F,N5L,N6F,N6L,IND,ICOUNT
      CHARACTER*10 FNODE,LFNAM
      integer dummy,i,idum,idummy,ind1,ind2,inze,ipass,irep,itrs,itrsgn
      integer ityp,jtrs,mreader,n5,n6,nkr,nkr2,nonzr,npass,nsec,nsec5
      integer nsec6,nze,INTTYP
C
C---------------Executable code--------------------------------------
C
C
      IF (SRTALL) THEN
      WRITE (IW,1004) 'OOOO',NV1,'VOOO',NV2,'VVOO',NV3,'VOVO',NV4,
     &                'VOVV',NV5,'VVVV',NV6
      ELSE
      WRITE (IW,1004) 'OOOO',NV1,'VOOO',NV2,'VVOO',NV3,'VOVO',NV4
      ENDIF
C
      CALL FLSHFO (IW)
      FOCKSP = EQNS.EQ.'FOCKSP'
      IF (INTERFACE_LOC(1:6).EQ.'DIRAC6'.AND..NOT.SPFR) THEN
         PERMSYM = .FALSE.
      ELSE
         PERMSYM = .TRUE.
      END IF
C
C  Check buffer size for passing of VOVV and VVVV integrals
C
      IF (MOD(NVT*IRW*RCW,NWORDS_ON_RECORD8).NE.0)
     &    CALL QUIT("INVALID BUFFER SIZE IN RDINTS")
C
C  INITIALIZE READING OF COULOMB INTEGRALS
C
      IF (MYPROC.EQ.0) THEN
         FNODE = 'MDCINT    '
      ELSE
         FNODE = LFNAM('MDCINT')
      ENDIF
      OPEN(MDINT,FILE=FNODE,FORM='UNFORMATTED')
      REWIND(MDINT)
      IF (INTERFACE_LOC(1:7).EQ.'MOLFDIR') THEN
         IF (MYPROC.EQ.MASTER) THEN
C           The oldest MOLFDIR file format
            READ (MDINT,ERR=10000,END=10) DATEX,TIMEX,NKR,
     &       (IDUM,I=1,4*NKR),(KR(I),KR(-I),I=1,NKR)
            GO TO 11
C           The newer MOLFDIR file format = identical to older DIRAC
   10       READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,
     &       (KR(I),KR(-I),I=1,NKR)
   11       CONTINUE
         ENDIF
      ELSEIF (INTERFACE_LOC(1:5).EQ.'DIRAC') THEN
         IF (MYPROC.EQ.MASTER.OR.INTERFACE_LOC(6:6).EQ.'6') THEN
         READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,
     &                                    (KR(I),KR(-I),I=1,NKR)
         ENDIF
         MREADER = MASTER
      ELSE
         WRITE (IW,'(//1X,3A)') 'Interface ',INTERFACE_LOC,
     & ' not defined'
         CALL QUIT('Adapt RDINTS')
      ENDIF
C
C     In the older schemes, the slaves have no knowledge of nkr and kr, broadcast
C     and check for consistency afterwards
C
#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         call interface_mpi_BCAST(NKR,1,
     &                  MASTER,global_communicator)
         call interface_mpi_BCAST(KR,2*NSP+1,
     &                  MASTER,global_communicator)
      END IF
#endif
      IF (2*NKR.NE.NSP)
     &   CALL QUIT("INCONSISTENT MRCONEE AND MDCINT FILES")
      ITRSGN = 1
C
C  Initialize writing to sorted integral files
C
      CALL WAIO_OPEN(ITAPT+0)
      CALL WAIO_OPEN(ITAPT+1)
      CALL WAIO_OPEN(ITAPT+2)
      CALL WAIO_OPEN(ITAPT+3)
      CALL WAIO_OPEN(ITAPT+4)
C
C  We may have already sorted these integrals when restarting, in that case jump ahead
C
      IF (IMSTAT(1).EQ.3) GOTO 150
C
C     RE-ENTRY POINT IN CASE THE GAUNT INTERACTION HAS TO BE ADDED
C
    1 ICOUNT = 0
C
C     GET THE OOOO, VOOO, VVOO AND VOVO INTEGRALS
C
      CALL XCOPY (NV1,A0,0,VOOOO,1)
      CALL XCOPY (NV2,A0,0,VVOOO,1)
      CALL XCOPY (NV3,A0,0,VVVOO,1)
      CALL XCOPY (NV4,A0,0,VVOVO,1)
      CALL XCOPY (NV1,A0,0,VOOOOA,1)
      CALL XCOPY (NV2,A0,0,VVOOOA,1)
      CALL XCOPY (NV3,A0,0,VVVOOA,1)
      CALL XCOPY (NV4,A0,0,VVOVOA,1)
C
      CALL DAYTIME (DATEX,TIMEX)
      WRITE(IW,*)
      WRITE(IW,*) 
     &     'Start sorting of integral classes at '//datex//' '//timex
      WRITE(IW,*)
C
  100 CALL GET2INTS (MDINT,KR,NSP,ITRSGN,1,INDK,INDL,
     &               NONZR,NZE,CBUF,IBUF,MREADER,PERMSYM)
      IF (NONZR.LT.0) THEN
C        We ran out of integrals on this node
         IF (MREADER.LT.NMPROC-1.AND.INTERFACE_LOC(6:6).EQ.'6') THEN
C            There will be more integrals to read on the next node
             MREADER = MREADER + 1
             GOTO 100
         ELSE
C            This was the last node, we are done processing the integrals, goto the antisymmetrization step
             GOTO 110
         ENDIF
      ELSE
         ICOUNT = ICOUNT + NONZR
      ENDIF
      DO 102 INZE = 1, NZE
            ITYP = IBUF(1,INZE)
            IREP = IBUF(2,INZE) 
            IND1 = IBUF(3,INZE)
            IND2 = IBUF(4,INZE)
            ITRS = IBUF(5,INZE)
            JTRS = IBUF(6,INZE)
            I    = IBUF(7,INZE)
            IF (ITYP.EQ.1) THEN
               IND = IOOOOTT(IREP) + (IND1-1)*NOOT(IREP) + IND2
               VOOOO(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VOOOO(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.-1) THEN         
               IND = IOOOOTT(IREP) + (IND1-1)*NOOT(IREP) + IND2
               VOOOOA(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VOOOOA(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.2) THEN       
               IND = IVOOOT(IREP) + (IND1-1)*NVO(IREP) + IND2
               VVOOO(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVOOO(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.-2) THEN
               IND = IVOOOT(IREP) + (IND1-1)*NVO(IREP) + IND2
               VVOOOA(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVOOOA(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.4) THEN
               IND = IVVOOTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
               VVVOO(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVVOO(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.-4) THEN
               IND = IVVOOTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
               VVVOOA(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVVOOA(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.6) THEN
               IND = IVOVO(IREP) + (IND1-1)*NVO(IREP) + IND2
               VVOVO(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVOVO(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.-6) THEN
               IND = IVOVO(IREP) + (IND1-1)*NVO(IREP) + IND2
               VVOVOA(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVOVOA(2,IND) = CBUF(2,I) * JTRS
            ENDIF
  102 CONTINUE
      GOTO 100
C
C     ADD THE DIRECT AND EXCHANGE INTEGRALS
C
 110  CALL XAXPY (NV1,-A1,VOOOOA,1,VOOOO,1)
      CALL XAXPY (NV2,-A1,VVOOOA,1,VVOOO,1)
      CALL XAXPY (NV3,-A1,VVVOOA,1,VVVOO,1)
      CALL XAXPY (NV4,-A1,VVOVOA,1,VVOVO,1)
C
C     ADD BREIT AND COULOMB INTEGRALS
C
      IF (BREIT_LOC.AND.ITRSGN.LT.0) THEN
         CALL RREAD(ITAPT+1,VOOOOA,NV1*RCW*IRW,1)
         CALL RREAD(ITAPT+2,VVOOOA,NV2*RCW*IRW,1)
         CALL RREAD(ITAPT+3,VVVOOA,NV3*RCW*IRW,1)
         CALL RREAD(ITAPT+4,VVOVOA,NV4*RCW*IRW,1)
         CALL XAXPY(NV1,A1,VOOOOA,1,VOOOO,1)
         CALL XAXPY(NV2,A1,VVOOOA,1,VVOOO,1)
         CALL XAXPY(NV3,A1,VVVOOA,1,VVVOO,1)
         CALL XAXPY(NV4,A1,VVOVOA,1,VVOVO,1)
      ENDIF
      CALL RWRIT(ITAPT+1,VOOOO,NV1*RCW*IRW,1)
      CALL RWRIT(ITAPT+2,VVOOO,NV2*RCW*IRW,1)
      CALL RWRIT(ITAPT+3,VVVOO,NV3*RCW*IRW,1)
      CALL RWRIT(ITAPT+4,VVOVO,NV4*RCW*IRW,1)
C
C     This is a restart point : save the files and update status.
C
  150 CONTINUE
      CALL WAIO_SAVE(ITAPT+1)
      CALL WAIO_SAVE(ITAPT+2)
      CALL WAIO_SAVE(ITAPT+3)
      CALL WAIO_SAVE(ITAPT+4)

      IF (IMSTAT(1).EQ.3) THEN
C        call relcc_restart_check('OOOO to VOVO sorted integrals ',IW)
         WRITE(IW,*)
         WRITE(IW,*) 'Restart: using sorted OOOO to VOVO integrals'
      ELSE
         CALL DAYTIME (DATEX,TIMEX)
         WRITE(IW,*)
         WRITE(IW,*) 
     &         'Sorting of first 4 classes done at '//datex//' '//timex
      ENDIF
      CALL FLSHFO (IW)
C
      IMSTAT(1) = 3
      IF(EQNS.NE.'RELADC') CALL CCDUMP
C
C Now get the VOVV integrals
C
      CALL WAIO_OPEN(ITAPT+5)
C
C  We may have already sorted these integrals when restarting, in that case jump ahead
C
      IF (IMSTAT(2).EQ.3) GOTO 450
C
      IF (.NOT.SRTALL) GOTO 600
C
C     Determine # of passes necessary for treating VOVV integrals
C     on a specific node
C
      NPASS = (NV5PAR-1) / NVT + 1
      NSEC5 = 1
      N5L = 0
      WRITE (IW,1005) NPASS,'VOVV'
      WRITE(IW,*) 
      CALL FLSHFO (IW)
C
C   start of IPASS loop for the VOVV
C ----------------------------
      DO IPASS = 1, NPASS
C ----------------------------

      N5F = N5L + 1
      N5L = MIN(N5L+NVT,NV5PAR)
      N5 = N5L - N5F + 1
      NSEC = N5*IRW*RCW / NWORDS_ON_RECORD
      IF (MOD(N5*IRW*RCW,NWORDS_ON_RECORD).NE.0) NSEC = NSEC + 1
C
      IF (MYPROC.EQ.MASTER.OR.INTERFACE_LOC(6:6).EQ.'6') THEN
         CLOSE (MDINT,STATUS='KEEP')
         IF (BREIT_LOC.AND.ITRSGN.EQ.-1) THEN
             OPEN(MDINT,FILE='MDBINT',FORM='UNFORMATTED')
         ELSE
             OPEN(MDINT,FILE=FNODE,FORM='UNFORMATTED')
         ENDIF
         READ (MDINT)
         MREADER = MASTER
      ENDIF
C
      CALL XCOPY (N5,A0,0,VVOVV,1)
      CALL XCOPY (N5,A0,0,VVOVVA,1)
C
  200 CALL GET2INTS (MDINT,KR,NSP,ITRSGN,14,INDK,INDL,
     &               NONZR,NZE,CBUF,IBUF,MREADER,PERMSYM)
      IF (NONZR.LT.0) THEN
         IF (MREADER.LT.NMPROC-1.AND.INTERFACE_LOC(6:6).EQ.'6') THEN
             MREADER = MREADER + 1
             GOTO 200
         ELSE
             GOTO 210
         ENDIF
      ENDIF
      DO 202 INZE = 1, NZE
            ITYP = IBUF(1,INZE)
            IREP = IBUF(2,INZE) 
            IND1 = IBUF(3,INZE)
            IND2 = IBUF(4,INZE)
            ITRS = IBUF(5,INZE)
            JTRS = IBUF(6,INZE)
            I    = IBUF(7,INZE)
            IF (ITYP.EQ.14) THEN       

C             Check absolute address and decide if this integral goes
C             to this node

               IND = IVOVVT(IREP) + (IND1-1)*NVO(IREP) + IND2
               IF (IND.GT.IDIST(3,2,IREP).AND.
     &           IND.LE.IDIST(4,2,IREP)) THEN

C                now check relative offset if we need it in
C                the current pass

                 IND = IND -IDIST(3,2,IREP) + IDIST(5,2,IREP)
                 IF (IND.GE.N5F.AND.IND.LE.N5L) THEN
                    IND = IND - N5F + 1
                    VVOVV(1,IND) = CBUF(1,I) * ITRS
                    IF (CARITH) VVOVV(2,IND) = CBUF(2,I) * JTRS
                 ENDIF
               ENDIF
            ELSEIF (ITYP.EQ.-14) THEN
               IND = IVOVVT(IREP) + (IND1-1)*NVO(IREP) + IND2

C             Check absolute address and decide if this integral goes
C             to this node

               IF (IND.GT.IDIST(3,2,IREP).AND.
     &           IND.LE.IDIST(4,2,IREP)) THEN

                 IND = IND -IDIST(3,2,IREP) + IDIST(5,2,IREP)
                 IF (IND.GE.N5F.AND.IND.LE.N5L) THEN
                    IND = IND - N5F + 1
                    VVOVVA(1,IND) = CBUF(1,I) * ITRS
                    IF (CARITH) VVOVVA(2,IND) = CBUF(2,I) * JTRS
                 ENDIF
               ENDIF
            ENDIF
  202 CONTINUE
      GOTO 200
C
C     ADD THE DIRECT AND EXCHANGE INTEGRALS
C
  210 CALL XAXPY (N5,-A1,VVOVVA,1,VVOVV,1)
C
C     ADD BREIT AND COULOMB INTEGRALS
C
      IF (BREIT_LOC.AND.ITRSGN.LT.0) THEN
         CALL RREAD(ITAPT+5,VVOVVA,N5*RCW*IRW,NSEC5)
         CALL XAXPY(N5,A1,VVOVVA,1,VVOVV,1)
      ENDIF
      CALL RWRIT(ITAPT+5,VVOVV,N5*RCW*IRW,NSEC5)
      NSEC5 = NSEC5 + NSEC
C
C  end of IPASS LOOP
C ----------------------
      CALL DAYTIME (DATEX,TIMEX)
      WRITE(IW,'(A,I6,A)') ' Pass',ipass,' ended at '//datex//' '//timex
      CALL FLSHFO (IW)
      ENDDO
C ----------------------
C
C     This is a restart point : save the files and update status.
C
  450 CONTINUE
      CALL WAIO_SAVE(ITAPT+5)

      IF (IMSTAT(2).EQ.3) THEN
C        call relcc_restart_check('VOVV sorted integrals        ',IW)
         WRITE(IW,*)
         WRITE(IW,*) 'Restart: using sorted VOVV integrals'
      ELSE
         CALL DAYTIME (DATEX,TIMEX)
         WRITE(IW,*)
         WRITE(IW,*) 'VOVV sorting done at '//datex//' '//timex
      ENDIF
      CALL FLSHFO (IW)

      IMSTAT(2) = 3
      IF(EQNS.NE.'RELADC') CALL CCDUMP

C
C     GET THE VVVV INTEGRALS
C
      CALL WAIO_OPEN(ITAPT+6)
C
      IF (IMSTAT(3).EQ.3) GOTO 550 
C
      NPASS = (NV6PAR-1) / NVT + 1
      NSEC6 = 1
      N6L = 0
      WRITE (IW,1005) NPASS,'VVVV'
      WRITE(IW,*) 
      CALL FLSHFO (IW)
C
C    start of IPASS loop for the VVVV
C --------------------------
      DO IPASS = 1, NPASS
C --------------------------

      N6F = N6L + 1
      N6L = MIN(N6L+NVT,NV6PAR)
      N6 = N6L - N6F + 1
      NSEC = N6*IRW*RCW / NWORDS_ON_RECORD
      IF (MOD(N6*IRW*RCW,NWORDS_ON_RECORD).NE.0) NSEC = NSEC + 1
C
      IF (MYPROC.EQ.MASTER.OR.INTERFACE_LOC(6:6).EQ.'6') THEN
         CLOSE (MDINT,STATUS='KEEP')
         IF (BREIT_LOC.AND.ITRSGN.EQ.-1) THEN
             OPEN(MDINT,FILE='MDBINT',FORM='UNFORMATTED')
         ELSE
             OPEN(MDINT,FILE=FNODE,FORM='UNFORMATTED')
         ENDIF
         READ (MDINT)
         MREADER = MASTER
      ENDIF
C
      CALL XCOPY (N6,A0,0,VVVVV,1)
      CALL XCOPY (N6,A0,0,VVVVVA,1)
C
  300 CALL GET2INTS (MDINT,KR,NSP,ITRSGN,16,INDK,INDL,
     &               NONZR,NZE,CBUF,IBUF,MREADER,PERMSYM)
      IF (NONZR.LT.0) THEN
         IF (MREADER.LT.NMPROC-1.AND.INTERFACE_LOC(6:6).EQ.'6') THEN
             MREADER = MREADER + 1
             GOTO 300
         ELSE
             GOTO 310
         ENDIF
      ENDIF
      DO 302 INZE = 1, NZE
            ITYP = IBUF(1,INZE)
            IREP = IBUF(2,INZE) 
            IND1 = IBUF(3,INZE)
            IND2 = IBUF(4,INZE)
            ITRS = IBUF(5,INZE)
            JTRS = IBUF(6,INZE)
            I    = IBUF(7,INZE)
            IF (ITYP.EQ.16) THEN       
C
C              The integral has the right type : check the absolute
C              address to see if it belongs to this node. 
C
               IND = IVVVVTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
               IF (IND.GT.IDIST(3,1,IREP).AND.
     &             IND.LE.IDIST(4,1,IREP)) THEN
C
C                 Now check the relative offset, whether we need it in
C                 the current pass.
C
                  IND = IND - IDIST(3,1,IREP) + IDIST(5,1,IREP)
                  IF (IND.GE.N6F.AND.IND.LE.N6L) THEN
                     IND = IND - N6F + 1
                     VVVVV(1,IND) = CBUF(1,I) * ITRS
                     IF (CARITH) VVVVV(2,IND) = CBUF(2,I) * JTRS
                  ENDIF
               ENDIF
            ELSEIF (ITYP.EQ.-16) THEN
              IND = IVVVVTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
              IF (IND.GT.IDIST(3,1,IREP).AND.
     &             IND.LE.IDIST(4,1,IREP)) THEN
C
C                 Now check the relative offset, whether we need it in
C                 the current pass.
C
                  IND = IND - IDIST(3,1,IREP) + IDIST(5,1,IREP)
                  IF (IND.GE.N6F.AND.IND.LE.N6L) THEN
                     IND = IND - N6F + 1
                     VVVVVA(1,IND) = CBUF(1,I) * ITRS
                     IF (CARITH) VVVVVA(2,IND) = CBUF(2,I) * JTRS
                  ENDIF
               ENDIF
            ENDIF
  302 CONTINUE
      GOTO 300
C
C     ADD THE DIRECT AND EXCHANGE INTEGRALS
C
  310 CALL XAXPY (N6,-A1,VVVVVA,1,VVVVV,1)
C
C     ADD BREIT AND COULOMB INTEGRALS
C
      IF (BREIT_LOC.AND.ITRSGN.LT.0) THEN
         CALL RREAD(ITAPT+6,VVVVVA,N6*RCW*IRW,NSEC6)
         CALL XAXPY(N6,A1,VVVVVA,1,VVVVV,1)
      ENDIF
      CALL RWRIT(ITAPT+6,VVVVV,N6*RCW*IRW,NSEC6)
      NSEC6 = NSEC6 + NSEC
C
C     end of IPASS loop for the VVVV
C --------------------------
      CALL DAYTIME (DATEX,TIMEX)
      WRITE(IW,'(A,I6,A)') ' Pass',ipass,' ended at '//datex//' '//timex
      CALL FLSHFO (IW)
      ENDDO
C --------------------------
C
  500 CLOSE(MDINT)
      CALL DAYTIME (DATEX,TIMEX)
      write(iw,*)
      WRITE(IW,*) 'VVVV sorting done at '//datex//' '//timex
      CALL FLSHFO (IW)  
C
C
C  NOW PASS THE BREIT INTEGRAL FILE IF DESIRED
C  Breit is only used with MOLFDIR in this version
C
      IF (BREIT_LOC.AND.ITRSGN.EQ.1) THEN
         WRITE (IW,1000) 'Coulomb',DATEX,TIMEX,ICOUNT
         OPEN(MDINT,FILE='MDBINT',FORM='UNFORMATTED')
         REWIND(MDINT)
         IF (MYPROC.EQ.MASTER) THEN
            READ (MDINT,ERR=10001,END=10001) DATEX,TIMEX,NKR2
            IF (NKR2.NE.NKR) GOTO 10003
         ENDIF
         ITRSGN = -1
         GOTO 1
      ELSEIF (BREIT_LOC) THEN
         WRITE (IW,1000) 'Breit',DATEX,TIMEX,ICOUNT
      ELSE
         WRITE (IW,1000) 'Coulomb',DATEX,TIMEX,ICOUNT
      ENDIF
C
C     This is a restart point : save the files and update status.
C
  550 CONTINUE
      CALL WAIO_SAVE(ITAPT+6)
      IF (IMSTAT(3).EQ.3) THEN
C        call relcc_restart_check('VVVV sorted integrals         ',IW)
         WRITE(IW,*)
         WRITE(IW,*) 'Restart: using sorted VVVV integrals'
      ENDIF

      IMSTAT(3) = 3
      IF(EQNS.NE.'RELADC') CALL CCDUMP
C
  600 WRITE (IW,1006)
C
      CALL GETFM (NSP,EPS,FOO,FVO,FVV,E0,E1)
C
C     This is a restart point : save the files and update status.
C
      CALL WAIO_SAVE(ITAPT+0)
      IMSTAT(4) = 3
      IF(EQNS.NE.'RELADC') CALL CCDUMP
C
      ESCF = E0 + E1 
C
C     SK - a first step towards a consistent treatment of 
C     NORECMP == .TRUE. + USEOE == .TRUE. 
C     when we run with the MOLMF(3) hamiltonian
      IF(NORECMP .and. USEOE) THEN
        WRITE(IW,'(/A,/A,/A,/A)') 
     & ' ***********************************************************',
     & ' * Info from RDINTS: No recomputation of SCF energy.       *',
     & ' * Value is taken from the file DFCOEF instead.            *',
     & ' ***********************************************************'
        ORIG_TOTSCF = 0.0D0
        IF( MYPROC .eq. MASTER )THEN
          INQUIRE(FILE='DFCOEF',EXIST=TOBE)
          IF(.NOT.TOBE) THEN
              WRITE(IW,'(A)')
     &          'RDINTS: Coefficient file DFCOEF not found !'
              CALL QUIT('RDINTS: Coefficient file DFCOEF not found !')
          ENDIF
          CALL OPNFIL(LUCOEF,'DFCOEF','OLD','RDINTS')
          CALL REACMO(LUCOEF,'DFCOEF',DUMMY,DUMMY,IDUMMY,ORIG_TOTSCF,1)
        ENDIF
#if defined VAR_MPI
        call interface_mpi_BCAST(ORIG_TOTSCF,1,MASTER,
     &                 global_communicator)
#endif
        ESCF = ORIG_TOTSCF
        WRITE(IW,1013) ESCF
      ELSE 
        WRITE(IW,1003) ECORE,E0,E1,ESCF,ESCF+ECORE
        ESCF = ESCF + ECORE
      END IF
      ETOT = ESCF
C
 1000 FORMAT (/' Reading ',A7,' integrals :'/' File date :',T18,A10
     &/' File time :',T20,A8,/' # of integrals',I20)
 1003 FORMAT (//' Nuclear repulsion + core energy :',T40,F25.15
     &/' Zero order electronic energy :',T40,F25.15 
     &/' First order electronic energy :',T40,F25.15 
     &/' Electronic energy :',T40,F25.15 
     &/' SCF energy :',T40,F25.15) 
 1004 FORMAT (/' Expanding and sorting integrals to unique types :'/
     & (' Type ',A4,' :',I16,' integrals'))
 1005 FORMAT (/' Need ',I6,' passes to sort ',A4,' integrals')
 1006 FORMAT (/' Finished sorting of integrals')
 1013 FORMAT (//' SCF energy :',T40,F25.15) 
      RETURN
10000 CALL QUIT('ERROR READING HEADER MDCINT')
10001 CALL QUIT('ERROR READING HEADER MDBINT')
10003 CALL QUIT('MDCINT AND MDBINT FILES ARE INCOMPATIBLE')
10010 CALL QUIT('ERROR READING INTEGRALS FROM MDCINT')
10011 CALL QUIT('ERROR READING INTEGRALS FROM MDBINT')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine relcc_restart_check(message,file)
         integer, intent(in) :: file
         character(30) :: message 

         write (file,'(A)') ' '
         WRITE (file,'(A)') ' Restart: using prior '//message(1:28)
         write (file,'(A)') ' '
         write (file,'(A)') '     Warning: No consistency checks yet'
         write (file,'(A)') ' '
         write (file,'(A)') '      Please verify the number of '
     &                    //'active spinors and the size of the model '
     &                    //'space'
         write (file,'(A)') '      (in FSCC calculations) has not'
     &                    //' changed between this restart and the '
     &                    //'prior run'

      end subroutine
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADDPROP (PROP_NAME,FF_STRENGTH,NSP,FOCKR,FOCKI)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Reads property integral files and adds the specified operator
C     to the Fock matrix (after multiplying it with the specified
C     field strength). This allows for treatment of operators that
C     break the Kramers but not the point group symmetry.
C     This version reads from the file created by PRTRAN
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
C
C---------------Calling variables--------------------------------------
C
      INTEGER NSP,NSPQ
      REAL*8 FOCKR(NSP,NSP),FOCKI(NSP,NSP)
      REAL*8 FF_STRENGTH(2), F_ANTI_HERMITIAN
      CHARACTER*8 PROP_NAME
C
C---------------Local variables--------------------------------------
C
      REAL*8 AH_CHECK
      PARAMETER (AH_CHECK=1.D-5)
      REAL*8 PROPR(NSP,NSP),PROPI(NSP,NSP)
      CHARACTER*32 ACHAR
      integer i,j,inop
C
C---------------Executable code--------------------------------------
C
      OPEN (MDPROP,FILE='MDPROP',FORM='UNFORMATTED')
      INOP = 0
    1 READ (MDPROP,ERR=10,END=11) ACHAR
      IF (ACHAR(1:8).NE.'********'.OR.ACHAR(25:32).NE.PROP_NAME) GOTO 1
      WRITE (IW,1000) PROP_NAME,ACHAR(9:16),ACHAR(17:24),FF_STRENGTH
      READ (MDPROP) ((PROPR(I,J),PROPI(I,J),I=1,NSP),J=1,NSP)
      CLOSE (MDPROP,STATUS='KEEP')
      GOTO 12
   10 INOP = 1
      GOTO 12
   11 INOP = 2
   12 CONTINUE
C
C     Error exit if the integrals could not be read
C
      IF (INOP.EQ.1) GOTO 101
      IF (INOP.EQ.2) GOTO 102
C
C     Multiply with the real part of the FF strength parameter
      NSPQ = NSP * NSP
      CALL DAXPY (NSPQ,FF_STRENGTH(1),PROPR,1,FOCKR,1)
      CALL DAXPY (NSPQ,FF_STRENGTH(1),PROPI,1,FOCKI,1)
C     Multiply with the imaginary part of the FF strength parameter
      CALL DAXPY (NSPQ,-FF_STRENGTH(2),PROPI,1,FOCKR,1)
      CALL DAXPY (NSPQ,FF_STRENGTH(2),PROPR,1,FOCKI,1)
C
C     Check the hermiticity of the new Fock matrix
C
      DO J = 1, NSP
         DO I = 1, NSP
            F_ANTI_HERMITIAN = FOCKR(I,J) - FOCKR(J,I)
            IF (ABS(F_ANTI_HERMITIAN).GT.AH_CHECK)
     &         CALL QUIT ("Perturbed (real)Fock matrix is nonhermitian")
            F_ANTI_HERMITIAN = FOCKI(I,J) + FOCKI(J,I)
            IF (ABS(F_ANTI_HERMITIAN).GT.AH_CHECK)
     &         CALL QUIT ("Perturbed (imag)Fock matrix is nonhermitian")
         ENDDO
      ENDDO
      RETURN
C
  101 WRITE (6,*) ' Error reading property ',PROP_NAME,' on file MDPROP'
      CALL QUIT(' Error reading property integrals')
  102 WRITE (6,*) ' Property ',PROP_NAME,' not found on file MDPROP'
      CALL QUIT(' Property integrals missing')
 1000 FORMAT (/' Read integral type ',A8,' created ',
     & A8,' storage info : ',A8
     &/' Added this to the Fock matrix with field strength ',2G16.8)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RDPROP(IPRNT,NAMEA,NSP,RDDIAG,PROP,IREPSPI,APHASE,
     &                  AVO,AOV,AOO,AVV)
C
      use interface_to_mpi
      implicit none
C
C---------------Description--------------------------------------------
C
C     Reads property integral files for operator A
C     If (RDDIAG) it reads also the diagonal blocks (OO and VV)
C     This version reads from the file created by PRTRAN
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
C
C---------------Calling variables--------------------------------------
C
      INTEGER NSP,IREPSPI(NSP,MXREP,2),iprnt
      REAL*8 AVO(*),AOV(*),AOO(*),AVV(*)
      REAL*8 PROP(2,NSP,NSP)
      REAL*8 APHASE(2)
      CHARACTER*8 NAMEA
      LOGICAL RDDIAG,PHASE
C
C---------------Local variables--------------------------------------
C
      CHARACTER*32 ACHAR
      REAL*8 AVOI,AOVI,AMAXR,AMAXI,ANORM
      integer a,aa,ab,ai,airep,arep,b,bb,brep,i,ii,iimax,ij,inop,irep
      integer j,jj,jrep,aamax
C
C---------------Executable code--------------------------------------
C
C
C     Read symmetry information (master node only)
C
      IF (MYPROC.EQ.MASTER) THEN
C
      OPEN (MDPROP,FILE='MDPROP',FORM='UNFORMATTED')
      INOP = 0
    1 READ (MDPROP,ERR=10,END=11) ACHAR
      IF (ACHAR(1:8).NE.'********'.OR.ACHAR(25:32).NE.NAMEA) GOTO 1
      IF (IPRNT.GE.1) WRITE (IW,1000) NAMEA,ACHAR(9:16),ACHAR(17:24)
      READ (MDPROP) PROP
      CLOSE (MDPROP,STATUS='KEEP')
      GOTO 12
   10 INOP = 1
      GOTO 12
   11 INOP = 2
   12 CONTINUE
C
      ENDIF
C
C     Distribute the information read from MDPROP
C
#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         call interface_mpi_BCAST(INOP,1,MASTER,
     &                  global_communicator)
         call interface_mpi_BCAST(PROP,2*NSP*NSP,MASTER,
     &                  global_communicator)
      END IF
#endif
C
C     Error exit if the integrals could not be read
C
      IF (INOP.EQ.1) GOTO 101
      IF (INOP.EQ.2) GOTO 102
C
      APHASE(1) = D1
      APHASE(2) = D0
      PHASE = .FALSE.
C
    2 AMAXI = D0
      AAMAX = 0
      IIMAX = 0
      AI = 1
      DO AIREP = 1, NREP
         DO IREP = 1, NREP
            AREP = MULTB (IREP,AIREP+NREP,2)
            AREP = MULTB (AREP,NREP+1,2)
            DO I = 1, NO(IREP)
               II = IREPSPI(I,IREP,1)
               DO A = 1, NV(AREP)
                  AA = IREPSPI(A,AREP,2)
                  IF (CARITH) THEN
                      AVO(RCW*AI-1) = PROP(1,AA,II)
                      AVO(RCW*AI)   = PROP(2,AA,II)
                      AOV(RCW*AI-1) = PROP(1,II,AA)
                      AOV(RCW*AI)   = PROP(2,II,AA)
                   ELSE
                      AVO(AI) = PROP(1,AA,II)*APHASE(1) -
     &                          PROP(2,AA,II)*APHASE(2)
                      AOV(AI) = PROP(1,II,AA)*APHASE(1) -
     &                          PROP(2,II,AA)*APHASE(2)
                      AVOI = PROP(1,AA,II)*APHASE(2) +
     &                       PROP(2,AA,II)*APHASE(1)
                      AOVI = PROP(1,II,AA)*APHASE(2) +
     &                       PROP(2,II,AA)*APHASE(1)
                      IF (ABS(AVOI).GT.ABS(AMAXI)) THEN
                         AMAXR = AVO(AI) 
                         AMAXI = AVOI
                         AAMAX = AA
                         IIMAX = II
                      ENDIF
                      IF (ABS(AOVI).GT.ABS(AMAXI)) THEN
                         AMAXR = AOV(AI) 
                         AMAXI = AOVI 
                         AAMAX = II
                         IIMAX = AA
                      ENDIF
                   ENDIF
                   AI = AI + 1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      IF (.NOT.CARITH.AND.ABS(AMAXI).GT.ACCUR) THEN
         IF (PHASE) THEN
            WRITE (IW,1010) AA,II,AMAXI
            CALL QUIT('USE COMPLEX ARITHMETICS')
         ELSE
            APHASE(1) = AMAXR
            APHASE(2) = AMAXI 
            ANORM = SQRT(AMAXR*AMAXR+AMAXI*AMAXI)
            APHASE(1) = APHASE(1) / ANORM
            APHASE(2) = APHASE(2) / ANORM
            PHASE = .TRUE.
            IF (IPRNT.GE.1) WRITE (6,1020) NAMEA,APHASE
            GOTO 2
         ENDIF
      ENDIF
C
      IF (.NOT.RDDIAG) RETURN
C
      IJ = 1
      DO JREP = 1, NREP
         IREP = JREP
         DO J = 1, NO(JREP)
            JJ = IREPSPI(J,JREP,1)
            DO I = 1, NO(IREP)
               II = IREPSPI(I,IREP,1)
               IF (CARITH) THEN
                  AOO(RCW*IJ-1) = PROP(1,II,JJ)
                  AOO(RCW*IJ)   = PROP(2,II,JJ)
               ELSE
                  AOO(IJ) = PROP(1,II,JJ)*APHASE(1) -
     &                      PROP(2,II,JJ)*APHASE(2)
               ENDIF
               IJ = IJ + 1
            ENDDO
         ENDDO
      ENDDO
C
      AB = 1
      DO BREP = 1, NREP
         AREP = BREP
         DO B = 1, NV(BREP)
            BB = IREPSPI(B,BREP,2)
            DO A = 1, NV(AREP)
               AA = IREPSPI(A,AREP,2)
               IF (CARITH) THEN
                  AVV(RCW*AB-1) = PROP(1,AA,BB)
                  AVV(RCW*AB)   = PROP(2,AA,BB)
               ELSE
                  AVV(AB) = PROP(1,AA,BB)*APHASE(1) -
     &                      PROP(2,AA,BB)*APHASE(2)
               ENDIF
               AB = AB + 1
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
  101 WRITE (6,*) ' Error reading property ',NAMEA,' on file MDPROP'
      CALL QUIT(' Error reading property integrals')
  102 WRITE (6,*) ' Property ',NAMEA,' not found on file MDPROP'
      CALL QUIT(' Property integrals missing')
 1000 FORMAT (/' Read integral type ',A8,' created ',
     & A8,' storage info : ',A8)
 1010 FORMAT (/' Largest imaginary part of matrix element',2I5,F10.2)
 1020 FORMAT (/' Property ',A8,' scaled with phase factor ',2F10.4)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET2INTS (MDINT,KR,NSP,ITRSGN,IATYP,INDK,INDL,
     &                     NONZR,NZE,CBUF,IBUF,MREADER,PERMSYM)
C
      use interface_to_mpi
      use spinor_indexing
      implicit none
C
C---------------Description--------------------------------------------
C
C     Reads block of 2-electron integrals from MOLFDIR file MDINT
C     Expands Kramers reduced list into full list.
C     Reorder from (ij|kl) to <ik||jl>
C     Expanded (16 permutations) indices are written in array IBUF :
C     IBUF(1,INZE) : Type (OOOO, VOOO, VVOO, VOVO, VVVO, VVVV)
C     IBUF(2,INZE) : Number of combined irrep
C     IBUF(3,INZE) : Index of ket pair
C     IBUF(4,INZE) : Index of bra pair
C     IBUF(5,INZE) : Sign under time-reversion (real part)
C     IBUF(6,INZE) : Sign under time-reversion and conjugations 
C     IBUF(7,INZE) : Pointer to integral values array 
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "ccpar.inc"
#include "symm.inc"
! This dependency on a DIRAC common block should go but otherwise NZ is not known.
#include "dgroup.h"
C
C---------------Calling variables--------------------------------------
C
      INTEGER KR(-NSP:NSP),INDK(*),INDL(*),MREADER
      REAL*8 CBUF(RCW,*)
      INTEGER IBUF(7,*)
      INTEGER MDINT,NSP,ITRSGN,IATYP,NONZR,NZE
      LOGICAL PERMSYM
C
C---------------Local variables--------------------------------------
C
      INTEGER IKR,JKR,KKR,LKR
      integer   INZ,l1,j1

      REAL*8 RDUM
      LOGICAL FOCKSP
      CHARACTER*10 NEXTFIL
      integer i,i0,i1,icase,icc,ijsg,ijsign,ikind,iktyp,ipi,itr1,itr2
      integer ityp,j,j0,jcase,jlind,jltyp,k,k0,k1,kiind,kityp,klsg,l,l0
      integer itr, itrsg_r,itrsg_i
C
C---------------Executable code--------------------------------------
C
      FOCKSP = EQNS.EQ.'FOCKSP'
C
C     Depending on the MOLTRA transformation algorithm the integrals may 
C     be distributed over the nodes or only on the master
C     They will always need to be dsitributed at this point
C 

   1  CONTINUE
      IF (MYPROC.EQ.MREADER) THEN
C
C     A number of different formats are possible. 
      IF (CARITH) THEN 
       IF (NZ.EQ.1) THEN
         READ (MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,
     &        (INDK(INZ),INDL(INZ),INZ=1,NONZR),
     &        (CBUF(1,INZ),INZ=1,NONZR) 
CMI      ... zero imaginary part of two-el. integrals
         if(rcw > 1)then
           DO INZ=1,NONZR
            CBUF(2,INZ) = 0.0D0
           END DO
         end if
       ELSE
C     For complex groups we read the real and imaginary part.
         IF (SPFR) THEN
          READ (MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,
     &         (INDK(INZ),INDL(INZ),INZ=1,NONZR),
     &         (CBUF(1,INZ),INZ=1,NONZR) 
         ELSE
          READ (MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,
     &         (INDK(INZ),INDL(INZ),INZ=1,NONZR),
     &         (CBUF(1,INZ),CBUF(2,INZ),INZ=1,NONZR) 
         ENDIF
       ENDIF

      ELSEIF (INTERFACE(1:7).NE.'MOLFDIR') THEN
C     For real groups we read only the real part in DIRAC format.
         READ (MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,
     &        (INDK(INZ),INDL(INZ),INZ=1,NONZR),
     &        (CBUF(1,INZ),INZ=1,NONZR) 
      ELSEIF (INTERFACE(1:7).EQ.'MOLFDIR') THEN
C     For real groups we skip the zero imaginary part in MOLFDIR format.
         READ (MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,
     &        (INDK(INZ),INDL(INZ),INZ=1,NONZR),
     &        (CBUF(1,INZ),RDUM,INZ=1,NONZR) 
      ENDIF

      IF (IKR.EQ.0) THEN
         IF (JKR.EQ.0) THEN
            NONZR = -1
         ELSE
            READ (MDINT) NEXTFIL
            CLOSE (MDINT,STATUS='KEEP')
            OPEN(MDINT,FILE=NEXTFIL,FORM='UNFORMATTED')
            GOTO 1
         ENDIF
      ENDIF

      ENDIF
C     End of reader branch, now all nodes are active again

#if defined (VAR_MPI)
      IF (NMPROC.GT.1 ) THEN
        call interface_mpi_BCAST(IKR,  1,MREADER,global_communicator)
        call interface_mpi_BCAST(JKR,  1,MREADER,global_communicator)
        call interface_mpi_BCAST(NONZR,1,MREADER,global_communicator)
        IF (NONZR.GT.0) THEN
          call interface_mpi_bcast_i1_work_f77(INDK,NONZR,
     &                   MREADER,global_communicator)
          call interface_mpi_bcast_i1_work_f77(INDL,NONZR,
     &                   MREADER,global_communicator)
          call interface_mpi_bcast_r2_work_f77(CBUF,RCW,NONZR,
     &                   MREADER,global_communicator)
        ENDIF
      ENDIF
#endif
      IF (NONZR.EQ.-1) RETURN

      IJSIGN = SIGN(1,IKR*JKR) * ITRSGN
C
C
      IF(.NOT.PERMSYM.AND..NOT.FOCKSP) THEN
C
C     We only need to account for time-reversal symmetry
C
      NZE = 0
      DO INZ = 1, NONZR
         KKR = INDK(INZ)
         LKR = INDL(INZ)
         I = KR(IKR)
         J = KR(JKR)
         K = KR(KKR)
         L = KR(LKR)
         ITRSG_R = 1
         ITRSG_I = 1
         DO ITR = 1, 2
            JLIND = PAIR_INDEX(J,L,4,1)
            IKIND = PAIR_INDEX(I,K,4,1)
            IF (JLIND.NE.0.AND.I.NE.K) THEN
C              Canonical integral : Direct or exchange
               NZE = NZE + 1
               IF (IKIND.NE.0) THEN
C                 Contribution to direct part
                  ITYP = (PAIR_INDEX(J,L,1,1)-1)*4
     &                 + PAIR_INDEX(I,K,1,1)
                  IBUF(1,NZE) = ITYP
                  IBUF(2,NZE) = PAIR_INDEX(J,L,3,1)
                  IBUF(3,NZE) = PAIR_INDEX(J,L,4,1)
                  IBUF(4,NZE) = PAIR_INDEX(I,K,4,1)
                  IBUF(5,NZE) = ITRSG_R
                  IBUF(6,NZE) = ITRSG_I
                  IBUF(7,NZE) = INZ
C                 Skip unnecessary integrals
                  IF (IATYP.NE.ITYP.AND.ITYP.GT.6) NZE=NZE-1
               ELSE
C                 Contribution to exchange part
                  ITYP = (PAIR_INDEX(J,L,1,1)-1)*4
     &                 + PAIR_INDEX(K,I,1,1)
                  IBUF(1,NZE) = - ITYP
                  IBUF(2,NZE) = PAIR_INDEX(J,L,3,1)
                  IBUF(3,NZE) = PAIR_INDEX(J,L,4,1)
                  IBUF(4,NZE) = PAIR_INDEX(K,I,4,1)
                  IBUF(5,NZE) = ITRSG_R
                  IBUF(6,NZE) = ITRSG_I
                  IBUF(7,NZE) = INZ
C                 Skip unnecessary integrals
                  IF (IATYP.NE.ITYP.AND.ITYP.GT.6) NZE=NZE-1
               ENDIF
            ENDIF
            I = KR(-IKR)
            J = KR(-JKR)
            K = KR(-KKR)
            L = KR(-LKR)
            ! Calculate sign for applying with time-reveral, no need to specifally
            ! account for Coulomb / Gaunt (ITRSGN=+/-1) because we operate on all 4 indices
            ITRSG_R = SIGN(1,IKR*JKR*KKR*LKR)
            ITRSG_I = - ITRSG_R
         END DO
      END DO

      ELSEIF(.NOT.PERMSYM.AND.FOCKSP) THEN
C     We are doing a Fock space calculation in which integrals containing active
C     spinors need to be stored at multiple places (as occ and as virt)
      NZE = 0
      DO INZ = 1, NONZR
         KKR = INDK(INZ)
         LKR = INDL(INZ)
         I = KR(IKR)
         J = KR(JKR)
         K = KR(KKR)
         L = KR(LKR)
         ITRSG_R = 1
         ITRSG_I = 1
         DO ITR = 1, 2
C
C           One integral can contribute to a maximum of
C           16 classes in Fock space CC calculations
C
            DO ICASE = 1, 4
               DO JCASE = 1, 4
                  JLTYP = PAIR_INDEX(J,L,1,ICASE)
                  IKTYP = PAIR_INDEX(I,K,1,JCASE)
                  KITYP = PAIR_INDEX(K,I,1,JCASE)
                  JLIND = PAIR_INDEX(J,L,4,JLTYP)
                  IKIND = PAIR_INDEX(I,K,4,IKTYP)
                  KIIND = PAIR_INDEX(K,I,4,KITYP)
                  IF (JLIND.NE.0.AND.IKIND.NE.0) THEN
C                   Canonical integral : direct
                    NZE = NZE + 1
                    ITYP = (JLTYP-1)*4 + IKTYP
                    IBUF(1,NZE) = ITYP
                    IBUF(2,NZE) = PAIR_INDEX(J,L,3,JLTYP)
                    IBUF(3,NZE) = JLIND
                    IBUF(4,NZE) = IKIND
                    IBUF(5,NZE) = ITRSG_R
                    IBUF(6,NZE) = ITRSG_I
                    IBUF(7,NZE) = INZ
C                   Skip unnecessary integrals
                    IF (IATYP.NE.ITYP.AND.ITYP.GT.6) NZE=NZE-1
                  ENDIF

                  IF (JLIND.NE.0.AND.KIIND.NE.0) THEN
C                   Canonical integral : exchange
                    NZE = NZE + 1
                    ITYP = (JLTYP-1)*4 + KITYP
                    IBUF(1,NZE) = - ITYP
                    IBUF(2,NZE) = PAIR_INDEX(J,L,3,JLTYP)
                    IBUF(3,NZE) = JLIND
                    IBUF(4,NZE) = KIIND
                    IBUF(5,NZE) = ITRSG_R
                    IBUF(6,NZE) = ITRSG_I
                    IBUF(7,NZE) = INZ
C                   Skip unnecessary integrals
                    IF (IATYP.NE.ITYP.AND.ITYP.GT.6) NZE=NZE-1
                  ENDIF
               ENDDO
            ENDDO
            I = KR(-IKR)
            J = KR(-JKR)
            K = KR(-KKR)
            L = KR(-LKR)
            ! Calculate sign for applying with time-reveral, no need to specifally
            ! account for Coulomb / Gaunt (ITRSGN=+/-1) because we operate on all 4 indices
            ITRSG_R = SIGN(1,IKR*JKR*KKR*LKR)
            ITRSG_I = - ITRSG_R
         END DO
      END DO

      ELSEIF (PERMSYM.AND.FOCKSP) THEN

C     Expand to the full list in which integrals may be needed at several places
C     the most complicated (and time-consuming) sort
      NZE = 0 
      DO 50 INZ = 1, NONZR
         KKR = INDK(INZ)
         LKR = INDL(INZ)
         I0 = KR(IKR)
         J0 = KR(JKR)
         IJSG = 1
         DO 40 ITR1 = 1, 2
            K0 = KR(KKR)
            L0 = KR(LKR)
            KLSG = 1
            DO 30 ITR2 = 1, 2
               I1 = I0
               J1 = J0
               K1 = K0
               L1 = L0
               DO 20 IPI = 1, 2
                  I = I1
                  J = J1
                  K = K1
                  L = L1
                  DO 10 ICC = 1, -1, -2
C
C                     One integral can contribute to a maximum of
C                     16 classes in Fock space CC calculations
C
                      DO ICASE = 1, 4
                      DO JCASE = 1, 4
                        JLTYP = PAIR_INDEX(J,L,1,ICASE)
                        IKTYP = PAIR_INDEX(I,K,1,JCASE)
                        KITYP = PAIR_INDEX(K,I,1,JCASE)
                        JLIND = PAIR_INDEX(J,L,4,JLTYP)
                        IKIND = PAIR_INDEX(I,K,4,IKTYP)
                        KIIND = PAIR_INDEX(K,I,4,KITYP)
                        IF (JLIND.NE.0.AND.IKIND.NE.0) THEN
C                         Canonical integral : direct
                          NZE = NZE + 1
                          ITYP = (JLTYP-1)*4 + IKTYP
                          IBUF(1,NZE) = ITYP
                          IBUF(2,NZE) = PAIR_INDEX(J,L,3,JLTYP)
                          IBUF(3,NZE) = JLIND
                          IBUF(4,NZE) = IKIND
                          IBUF(5,NZE) = IJSG * KLSG
                          IBUF(6,NZE) = ICC * IJSG * KLSG
                          IBUF(7,NZE) = INZ
C                         Skip unnecessary integrals
                          IF (IATYP.NE.ITYP.AND.ITYP.GT.6) NZE=NZE-1
                        ENDIF
C
                        IF (JLIND.NE.0.AND.KIIND.NE.0) THEN
C                         Canonical integral : exchange
                          NZE = NZE + 1
                          ITYP = (JLTYP-1)*4 + KITYP
                          IBUF(1,NZE) = - ITYP
                          IBUF(2,NZE) = PAIR_INDEX(J,L,3,JLTYP)
                          IBUF(3,NZE) = JLIND
                          IBUF(4,NZE) = KIIND
                          IBUF(5,NZE) = IJSG * KLSG
                          IBUF(6,NZE) = ICC * IJSG * KLSG
                          IBUF(7,NZE) = INZ
C                         Skip unnecessary integrals
                          IF (IATYP.NE.ITYP.AND.ITYP.GT.6) NZE=NZE-1
                        ENDIF
                      ENDDO
                      ENDDO
                     J = I1
                     I = J1
                     L = K1
                     K = L1
   10             CONTINUE
                  I1 = K0
                  J1 = L0
                  K1 = I0
                  L1 = J0
   20          CONTINUE
               K0 = KR(-LKR)
               L0 = KR(-KKR)
               KLSG = SIGN(1,KKR*LKR) * ITRSGN
   30       CONTINUE
            I0 = KR(-JKR)
            J0 = KR(-IKR)
            IJSG = IJSIGN
   40    CONTINUE
   50 CONTINUE
      ELSEIF (PERMSYM.AND..NOT.FOCKSP) THEN
C
C     integrals have permutation symmetry, we need to expand
C     to the full list
C
      NZE = 0 
      DO 150 INZ = 1, NONZR
         KKR = INDK(INZ)
         LKR = INDL(INZ)
         I0 = KR(IKR)
         J0 = KR(JKR)
         IJSG = 1
         DO 140 ITR1 = 1, 2
            K0 = KR(KKR)
            L0 = KR(LKR)
            KLSG = 1
            DO 130 ITR2 = 1, 2
               I1 = I0
               J1 = J0
               K1 = K0
               L1 = L0
               DO 120 IPI = 1, 2
                  I = I1
                  J = J1
                  K = K1
                  L = L1
                  DO 110 ICC = 1, -1, -2
                       JLIND = PAIR_INDEX(J,L,4,1)
                       IKIND = PAIR_INDEX(I,K,4,1)
                       IF (JLIND.NE.0.AND.I.NE.K) THEN
C                         Canonical integral : Direct or exchange
                          NZE = NZE + 1
                          IF (IKIND.NE.0) THEN
C                            Contribution to direct part
                             ITYP = (PAIR_INDEX(J,L,1,1)-1)*4
     &                            + PAIR_INDEX(I,K,1,1)
                             IBUF(1,NZE) = ITYP
                             IBUF(2,NZE) = PAIR_INDEX(J,L,3,1)
                             IBUF(3,NZE) = PAIR_INDEX(J,L,4,1)
                             IBUF(4,NZE) = PAIR_INDEX(I,K,4,1)
                             IBUF(5,NZE) = IJSG * KLSG
                             IBUF(6,NZE) = ICC * IJSG * KLSG
                             IBUF(7,NZE) = INZ
C                            Skip unnecessary integrals
                             IF (IATYP.NE.ITYP.AND.ITYP.GT.6) NZE=NZE-1
                          ELSE
C                            Contribution to exchange part
                             ITYP = (PAIR_INDEX(J,L,1,1)-1)*4
     &                            + PAIR_INDEX(K,I,1,1)
                             IBUF(1,NZE) = - ITYP
                             IBUF(2,NZE) = PAIR_INDEX(J,L,3,1)
                             IBUF(3,NZE) = PAIR_INDEX(J,L,4,1)
                             IBUF(4,NZE) = PAIR_INDEX(K,I,4,1)
                             IBUF(5,NZE) = IJSG * KLSG
                             IBUF(6,NZE) = ICC * IJSG * KLSG
                             IBUF(7,NZE) = INZ
C                            Skip unnecessary integrals
                             IF (IATYP.NE.ITYP.AND.ITYP.GT.6) NZE=NZE-1
                          ENDIF
                       ENDIF
                     J = I1
                     I = J1
                     L = K1
                     K = L1
  110             CONTINUE
                  I1 = K0
                  J1 = L0
                  K1 = I0
                  L1 = J0
  120          CONTINUE
               K0 = KR(-LKR)
               L0 = KR(-KKR)
               KLSG = SIGN(1,KKR*LKR) * ITRSGN
  130       CONTINUE
            I0 = KR(-JKR)
            J0 = KR(-IKR)
            IJSG = IJSIGN
  140    CONTINUE
  150 CONTINUE
      ENDIF
C
      RETURN
10010 WRITE(6,'(/A/A,L10,2A)')
     & 'END OF FILE branch taken when reading 2-e integrals',
     & 'CARITH =',CARITH,', INTERFACE =',INTERFACE
      CALL QUIT('END OF FILE READING 2-E INTEGRALS')
10020 WRITE(6,'(/A/A,L10,2A)')
     & 'ERROR branch taken when reading 2-e integrals',
     & 'CARITH =',CARITH,', INTERFACE =',INTERFACE
      CALL QUIT('ERROR READING 2-E INTEGRALS')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETFM (NSP,EPS,FOO,FVO,FVV,E0,E1)
C
      use interface_to_mpi
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generates or converts Fock matrix. Calculates E0 and E1
C
C---------------Routines called----------------------------------------
C
C     MAKEFM, MPI_BCAST
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 E0,E1,EPS(*)
      REAL*8 FOO(*),FVO(*),FVV(*)
      INTEGER NSP
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "complex.inc"
#include "results.inc"
#include "files.inc"
#include "ccpar.inc"
#include "ihm.inc"
C---------------Local variables--------------------------------------
C
      LOGICAL FOCKSP
      integer i,ii,ij,ind,ioper,irep,j
C
C---------------Executable code--------------------------------------
C
      FOCKSP = EQNS.EQ.'FOCKSP'
C
C     We excute all of this only on master and synchronize at the end
C
      IF (MYPROC.EQ.MASTER) THEN
C
      IF (NORECMP) THEN
C        If the user specified norecmp = no recompute
C        we assume that the Fock matrix is diagonal and given by the orbital energies
         CALL XCOPY (NFOO,A0,0,FOO,1)
         CALL XCOPY (NFVO,A0,0,FVO,1)
         CALL XCOPY (NFVV,A0,0,FVV,1)
C        Add orbital energies to the diagonal
         CALL ADD_EPS_TO_FM (FOCKSP,EPS,FOO,FVO,FVV)
      ELSE
C        calculation of Fock matrix from scratch, effective 1e contribution
         CALL INIT_1E_FM (FOCKSP,NSP,FOO,FVO,FVV)
C        Complete the Fock matrix by adding the two-electron contributions
         CALL ADD_2E_TO_FM (FOO,FVO,FVV)
      END IF
C
C     We calculate the 1-determinant energy, which can be compared
C     to the value read from the MRCONEE file as a check on the correctness
C     of the index transformations.
C
      CALL COMPUTE_HF_ENERGY (FOO,E0,E1)
C
C     Write human-readable output (mostly concerning checks of the unperturbed Fock matrix)
C
      CALL WRITE_FM_OUTPUT (E0,E1,EPS,FOO,FVV)
C
C     If desired: add a finite field. The thus perturbed Fock matrix will be used
C     in the following so we need to recompute the energy for consistency
C
      IF (ADD_FINITE_FIELD) THEN
         CALL ADD_PERTURBATION (FOCKSP,NSP,FOO,FVO,FVV)
         CALL COMPUTE_HF_ENERGY (FOO,E0,E1)
      END IF

C
C     We normally we take the diagonal of the Fock matrix as zeroth order Hamiltonian
C     this can be modified to using the original orbital energies
C     (note that there is no difference between these two options if norecmp is true)

      IF (.NOT.USEOE)  CALL COPY_DIAG_TO_EPS (EPS,FOO,FVV)

      ENDIF ! end of serial block executed only on master

C     Provide the Fock matrix to the other nodes as well
#if defined (VAR_MPI)
      IF (NMPROC.GT.1) THEN
         call interface_mpi_bcast_r1_work_f77(FOO,NFOO*RCW,MASTER,
     &                  global_communicator)
         call interface_mpi_bcast_r1_work_f77(FVO,NFVO*RCW,MASTER,
     &                  global_communicator)
         call interface_mpi_bcast_r1_work_f77(FVV,NFVV*RCW,MASTER,
     &                  global_communicator)
         call interface_mpi_bcast_r1_work_f77(EPS,NSP,MASTER,
     &                  global_communicator)
      END IF
#endif

C     Write Fock matrix to file for restart and post-processing purposes
      CALL FMTOFILE (LTR,FVO,FOO,FVV)

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WRITE_FM_OUTPUT (E0,E1,EPS,FOO,FVV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Writes human-readable output regarding Fock matrix construction.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 E0,E1,EPS(*)
      REAL*8 FOO(*),FVV(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "complex.inc"
#include "results.inc"
#include "files.inc"
#include "ihm.inc"
C---------------Local variables--------------------------------------
C
      REAL*8 EPSDIF
      LOGICAL FOCKSP
      CHARACTER*1 INPI
      CHARACTER*2 ORBTYP
      integer i,ii,irep,j
C
C---------------Executable code--------------------------------------
C
      FOCKSP = EQNS.EQ.'FOCKSP'
C
      WRITE (IW,1000) E0+E1
C
      WRITE (IW,1001)
      II = 1
      I  = 1
      DO IREP = 1, NREP
         DO J = 1, NO(IREP)
            ORBTYP = 'O '
            IF (FOCKSP) THEN
               IF (J.LE.NIO(IREP)) THEN
                  ORBTYP = 'Oi'
               ELSEIF (J.LE.NIO(IREP)+NAO(IREP)) THEN
                  ORBTYP = 'Oa'
               ELSE 
                  ORBTYP = 'Va'
               ENDIF
            ENDIF
            EPSDIF = ABS(EPS(I)-FOO(II))
            IF (ABS(FOO(II)).GT.1.D0) EPSDIF = EPSDIF / ABS(FOO(II))

            IF (EPSDIF.GT.1.D-7.OR.IPRNT.GE.1) THEN
CMI    ... addition for IH-FSCC
             IF (FOCKSP.AND.DOIH) THEN
               INPI=" "
CMI            IF (IPIORB(I).EQ.1) INPI="*"
               IF (IPIORB(I).EQ.1.AND.
     &            (ORBTYP.EQ.'Va'.OR.ORBTYP.EQ.'Oa')) INPI="*"
               write(IW,'(2x,A2,I4,I5,2X,A4,4X,2F16.10,1X,A)') 
     &          ORBTYP,J,I,REPNA(IREP),EPS(I),FOO(II),INPI 
             ELSE
               WRITE (IW,1002) ORBTYP,J,I,REPNA(IREP),EPS(I),FOO(II)
             ENDIF
            ENDIF

            I  = I  + 1
            II = II + (NO(IREP) + 1) * RCW
         ENDDO
         II = II - NO(IREP) * RCW
      ENDDO
      II = 1
      DO IREP = 1, NREP
         DO J = 1, NV(IREP)
            ORBTYP = 'V '
            IF (FOCKSP) THEN
               IF (J.LE.NAV(IREP)) THEN
                  ORBTYP = 'Va'
               ELSEIF (J.LE.NAV(IREP)+NIV(IREP)) THEN
                  ORBTYP = 'Vi'
               ELSE 
                  ORBTYP = 'Oa'
               ENDIF
            ENDIF
            EPSDIF = ABS(EPS(I)-FVV(II))
            IF (ABS(FVV(II)).GT.1.D0) EPSDIF = EPSDIF / ABS(FVV(II))

            IF (EPSDIF.GT.1.D-7.OR.IPRNT.GE.2) THEN
             IF (FOCKSP.AND.DOIH) THEN
               INPI=" "
               IF (IPIORB(I).EQ.1.AND.
     &            (ORBTYP.EQ.'Va'.OR.ORBTYP.EQ.'Oa')) INPI="*"
CMI            IF (IPIORB(I).EQ.1) INPI="*"
               write(IW,'(2x,A2,I4,I5,2X,A4,4X,2F16.10,1X,A)') 
     &          ORBTYP,J,I,REPNA(IREP),EPS(I),FVV(II),INPI 
             ELSE
               WRITE (IW,1002) ORBTYP,J,I,REPNA(IREP),EPS(I),FVV(II)
             ENDIF
            ENDIF

            I  = I  + 1
            II = II + (NV(IREP) + 1) * RCW
         ENDDO
         II = II - NV(IREP) * RCW
      ENDDO
C
      IF (IPRNT.GE.0.AND.(FOCKSP.AND.DOIH)) THEN
        WRITE (IW,'(/1X,A)')
     & 'Spinors in the IH Fock space Pi are marked with "*"'

        WRITE(IW,'(1X,A,F20.7,A,F20.7,A)')
     &  'Energy range for the occupied Pi space is: <',
     &  ER_IH(1,1),';',ER_IH(2,1),'>'

        WRITE(IW,'(1X,A,F20.7,A,F20.7,A)')
     &  'Energy range for the virtual Pi space is : <',
     &  ER_IH(1,2),';',ER_IH(2,2),'>'
      ENDIF

      IF (USEOE) THEN
         WRITE (IW,1003)
      ELSE
         WRITE (IW,1004)
      ENDIF
      WRITE (IW,1005)
C
 1000 FORMAT (//' Single determinant electronic energy :',T40,F25.15)
 1001 FORMAT (//' Checking the orbital energies, the program computes',
     & ' the diagonal elements of the'/' reconstructed Fock matrix.',
     & ' Differences with the reference orbital energies'/
     & ' are given if above a treshold or if iprnt > 1'
     &//' Spinor   Abelian Rep.',9X,'Energy   Recalc. Energy')
 1002 FORMAT (2x,A2,I4,I5,2X,A4,4X,2F16.10)
 1003 FORMAT(/' The original energies (left column) are used in',
     &' perturbation expressions.')
 1004 FORMAT(/' The diagonal elements of the recomputed Fock matrix',
     &' (right column) are used in perturbation expressions.')
 1005 FORMAT(/' Use the perturbative values (MP2, CCSD[T]/(T)/-T)',
     &' with care, especially',/' in  open shell calculations',
     &' because the orbitals need not always be'/
     &' semi-canonical as was assumed in the derivation',
     &' of the expressions.'/' The missing terms may be important !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INIT_1E_FM (FOCKSP,NSP,FOO,FVO,FVV)
      use spinor_indexing
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Initialize Fock matrix with effective 1-electron matrix elements 
C     that are read from MRCONEE
C
C---------------Routines called----------------------------------------
C
C     None
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      LOGICAL FOCKSP
      REAL*8 FOO(*),FVO(*),FVV(*)
      INTEGER NSP
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
C---------------Local variables--------------------------------------
C
      integer i,icase,ii,ij,ijoo,ijvo,ind,irep,ityp,j,ncase,nspq
      real*8, allocatable :: buf1(:),buf2(:)
C
C---------------Executable code--------------------------------------
C
      NSPQ = NSP*NSP
C
      allocate (buf1(nspq),buf2(nspq)) ! scratch arrays for unsorted hamiltonian matrix elements
C
C  READ ONE ELECTRON INTEGRALS
C
      OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED')
      READ (MRCONEE)
      READ (MRCONEE)
      READ (MRCONEE)
      READ (MRCONEE)
      READ (MRCONEE)
      READ (MRCONEE) (BUF1(IJ),BUF2(IJ),IJ=1,NSPQ)
      CLOSE (MRCONEE)

C  ORDER AND PUT IN FOCK MATRIX
C
      IF (FOCKSP) THEN
         NCASE = 4
      ELSE
         NCASE = 1
      ENDIF
C
      DO ICASE = 1, NCASE
         IJ = 0
         DO J = 1, NSP
            DO I = 1, NSP
               IJ = IJ + 1
               ITYP = PAIR_INDEX(I,J,1,ICASE)
               IND = PAIR_INDEX(I,J,2,ICASE)
               IF (FOCKSP) IND = PAIR_INDEX(I,J,2,ITYP)
               IF (IND.NE.0) THEN
               IF (ITYP.EQ.1) THEN
                  IF (CARITH) THEN
                     FOO(RCW*IND-1) = BUF1(IJ)
                     FOO(RCW*IND)   = BUF2(IJ)
                  ELSE
                     FOO(IND) = BUF1(IJ)
                  ENDIF
               ELSEIF (ITYP.EQ.2) THEN
                  IF (CARITH) THEN
                     FVO(RCW*IND-1) = BUF1(IJ)
                     FVO(RCW*IND)   = BUF2(IJ)
                  ELSE
                     FVO(IND) = BUF1(IJ)
                  ENDIF
               ELSEIF (ITYP.EQ.4) THEN
                  IF (CARITH) THEN
                     FVV(RCW*IND-1) = BUF1(IJ)
                     FVV(RCW*IND)   = BUF2(IJ)
                  ELSE
                     FVV(IND) = BUF1(IJ)
                  ENDIF
               ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      deallocate (buf1,buf2)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADD_PERTURBATION (FOCKSP,NSP,FOO,FVO,FVV)
      use spinor_indexing
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Adds property operators multiplied by finite field strength to Fock matrix
C
C---------------Routines called----------------------------------------
C
C     MAKEFM, MPI_BCAST
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FOO(*),FVO(*),FVV(*)
      INTEGER NSP
      LOGICAL FOCKSP
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "inpt.inc"
#include "complex.inc"
#include "files.inc"
C---------------Local variables--------------------------------------
C
      REAL*8 FF_STRENGTH(2)
      CHARACTER*8 PROP_NAME
      integer i,icase,ii,ij,ijoo,ijvo,ind,ioper,irep,ityp,j,ncase,nspq
      real*8, allocatable :: buf1(:),buf2(:)
C
C---------------Executable code--------------------------------------
C
      NSPQ = NSP*NSP
C
      allocate (buf1(nspq),buf2(nspq)) ! scratch arrays for unsorted matrix elements
      buf1 = 0.D0
      buf2 = 0.D0
C
C     We add one or more property operators here, each scaled with a field strength

      WRITE(IW,"(/,1X,A,I2,A)") "Going to add ",NFFOPER,
     &    " finite field operator(s) to the unrestricted Fock matrix." 
      DO IOPER=1,NFFOPER
         FF_STRENGTH(1) = FF_PROP_STRENGTHS(1,IOPER)
         FF_STRENGTH(2) = FF_PROP_STRENGTHS(2,IOPER)
         PROP_NAME = FF_PROP_NAMES(IOPER)
         CALL ADDPROP (PROP_NAME,FF_STRENGTH,NSP,BUF1,BUF2)
      ENDDO
C
C  
C  ORDER AND PUT IN FOCK MATRIX
C
      IF (FOCKSP) THEN
         NCASE = 4
      ELSE
         NCASE = 1
      ENDIF
C
      DO ICASE = 1, NCASE
         IJ = 0
         DO J = 1, NSP
            DO I = 1, NSP
               IJ = IJ + 1
               ITYP = PAIR_INDEX(I,J,1,ICASE)
               IND = PAIR_INDEX(I,J,2,ICASE)
               IF (FOCKSP) IND = PAIR_INDEX(I,J,2,ITYP)
               IF (IND.NE.0) THEN
               IF (ITYP.EQ.1) THEN
                  IF (CARITH) THEN
                     FOO(RCW*IND-1) = FOO(RCW*IND-1) + BUF1(IJ)
                     FOO(RCW*IND)   = FOO(RCW*IND)   + BUF2(IJ)
                  ELSE
                     FOO(IND) = FOO(IND) + BUF1(IJ)
                  ENDIF
               ELSEIF (ITYP.EQ.2) THEN
                  IF (CARITH) THEN
                     FVO(RCW*IND-1) = FVO(RCW*IND-1) + BUF1(IJ)
                     FVO(RCW*IND)   = FVO(RCW*IND)   + BUF2(IJ)
                  ELSE
                     FVO(IND) = FVO(IND) + BUF1(IJ)
                  ENDIF
               ELSEIF (ITYP.EQ.4) THEN
                  IF (CARITH) THEN
                     FVV(RCW*IND-1) = FVV(RCW*IND-1) + BUF1(IJ)
                     FVV(RCW*IND)   = FVV(RCW*IND)   + BUF2(IJ)
                  ELSE
                     FVV(IND) = FVV(IND) + BUF1(IJ)
                  ENDIF
               ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      deallocate (buf1,buf2)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADD_2E_TO_FM (FOO,FVO,FVV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generates 2e-part of Fock matrix directly from 2e-integrals
C
C---------------Routines called----------------------------------------
C
C     GETOOOO, GETVOOO, GETVOVO, SRT1TT4, SRT1ST4, SRT1SS4
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FOO(*),FVO(*),FVV(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "files.inc"
C
C---------------Local variables--------------------------------------
C
      real*8, allocatable :: buf1(:),buf2(:)
      integer ii,irep,k,kk,krep
      integer*8 nvtdim1, nvtdim2
C
C---------------Executable code--------------------------------------
C

C     Calculate size needed for buffer array 1: should be able to hold all unsorted integrals of types 1, 2, and 4
      nvtdim1 = INT(NV1,8)
      nvtdim1 = MAX(nvtdim1,INT(NV2,8))
      nvtdim1 = MAX(nvtdim1,INT(NV4,8))
C     Calculate size needed for buffer array 2: should be able to hold all sorted integrals of types 1, 2, and 4
C     for which the product irrep is totally symmetric. This is usually a smaller number, but due to less permutation
C     symmetry it may become larger than nvtdim1 in case of no point group symmetry.
      nvtdim2 = INT(JOOOO(2),8)
      nvtdim2 = MAX(nvtdim2,INT(JVOOO(2),8))
      nvtdim2 = MAX(nvtdim2,INT(JVVOO(2),8))

      allocate(buf1(rcw*nvtdim1))
      allocate(buf2(rcw*nvtdim2))

      CALL GETOOOO (BUF1)
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NO,NO,NO,NO,MOO,JOOOO,JJOO,JJOO,
     &              BUF1,BUF2)
      KK = 1
      DO KREP = 1, NREP
         DO K = 1, NO(KREP)
C---------------------------------------------------------------------
C  F(IJ) = F(IJ) + V(IJ,KK)
C---------------------------------------------------------------------
            IF (K.LE.NE(KREP)) CALL XAXPY(NFOO,A1,BUF2(KK),1,FOO,1)
            KK = KK + (NO(KREP)+1) * NFOO * RCW
         ENDDO
         KK = KK - NO(KREP) * NFOO * RCW
      ENDDO
C
      CALL GETVOOO (BUF1)
      CALL SRT1ST4 (NREP,MULTB,LTR,LFA,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &              BUF1,BUF2)
      KK = 1
      DO KREP = 1, NREP
         DO K = 1, NO(KREP)
C---------------------------------------------------------------------
C  F(AI) = F(AI) + V(AI,KK)
C---------------------------------------------------------------------
            IF (K.LE.NE(KREP)) CALL XAXPY(NFVO,A1,BUF2(KK),1,FVO,1)
            KK = KK + (NO(KREP)+1) * NFVO * RCW
         ENDDO
         KK = KK - NO(KREP) * NFVO * RCW
      ENDDO
C
      CALL GETVOVO (BUF1)
      CALL SRT1SS4 (NREP,MULTB,LTR,LFA,NV,NO,NV,NO,MVV,JVVOO,JJVV,JJOO,
     &              BUF1,BUF2)
      KK = 1
      DO KREP = 1, NREP
         DO K = 1, NO(KREP)
C---------------------------------------------------------------------
C  F(AB) = F(AB) + V(AB,KK)
C---------------------------------------------------------------------
            IF (K.LE.NE(KREP)) CALL XAXPY(NFVV,A1,BUF2(KK),1,FVV,1)
            KK = KK + (NO(KREP)+1) * NFVV * RCW
         ENDDO
         KK = KK - NO(KREP) * NFVV * RCW
      ENDDO
C
      deallocate(buf1)
      deallocate(buf2)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE COMPUTE_HF_ENERGY (FOO,E0,E1)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Evaluate HF energy contributions 
C
C---------------Routines called----------------------------------------
C
C     None
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 E0,E1
      REAL*8 FOO(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CSUM(2)
      real*8, allocatable :: buf1(:),buf2(:)
      integer ii,irep,k,kk,krep
      integer*8 nvtdim1, nvtdim2
C
C---------------Executable code--------------------------------------
C

C     Calculate size needed for buffer array 1: should be able to hold all unsorted integrals of type 1
      nvtdim1 = INT(NV1,8)
C     Calculate size needed for buffer array 2: should be able to hold all sorted integrals of type 1
C     for which the product irrep is totally symmetric. This is usually a smaller number, but due to less permutation
C     symmetry it may become larger than nvtdim1 in case of no point group symmetry.
      nvtdim2 = INT(JOOOO(2),8)

      allocate(buf1(rcw*nvtdim1))
      allocate(buf2(rcw*nvtdim2))

      E0 = AR0
      E1 = AR0
C
C---------------------------------------------------------------------
C  E0 = E0 + FOO(I,I)
C---------------------------------------------------------------------
      II = 1
      DO IREP = 1, NREP
         CALL XSUM(CSUM,NE(IREP),FOO(II),NO(IREP)+1)
         E0 = E0 + CSUM(1)
         II = II + NO(IREP) * NO(IREP) * RCW
      ENDDO
C
      CALL GETOOOO (BUF1)
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NO,NO,NO,NO,MOO,JOOOO,JJOO,JJOO,
     &              BUF1,BUF2)
      KK = 1
      DO KREP = 1, NREP
         DO K = 1, NO(KREP)
C---------------------------------------------------------------------
C  E1 = E1 - V(II,KK)
C---------------------------------------------------------------------
            II = 0
            DO IREP = 1, NREP
               IF (K.LE.NE(KREP) .AND. NE(IREP).GT.0) THEN
                  CALL XSUM(CSUM,NE(IREP),BUF2(II+KK),NO(IREP)+1)
                  E1 = E1 - CSUM(1)
               ENDIF
               II = II + NO(IREP) * NO(IREP) * RCW
            ENDDO
            KK = KK + (NO(KREP)+1) * NFOO * RCW
         ENDDO
         KK = KK - NO(KREP) * NFOO * RCW
      ENDDO
      E1 = E1 * DP5
C
      deallocate(buf1)
      deallocate(buf2)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADD_EPS_TO_FM (FOCKSP,EPS,FOO,FVO,FVV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Adds orbital energies (EPS) at diagonal positions in Fock matrix
C
C---------------Routines called----------------------------------------
C
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      LOGICAL FOCKSP
      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "files.inc"
C
C---------------Local variables--------------------------------------
C
      integer i,ii,irep,j,k,kk,krep,ijoo,ijvo
C
C---------------Executable code--------------------------------------
C
C     We use the diagonal matrix of orbital energies and do not recalculate
C     the Fock matrix. This will only work if the orbitals are
C     generated with the same fock operator as used here (in general
C     only for closed shell systems) !
C
C     Note that we do not initialize but add the orbital energies to FOO and FVV
C     This is to allow for potential addition of a finite field operator elsewhere
C
      I = 1
      II = 1
      DO IREP = 1, NREP
         DO J = 1, NO(IREP)
            FOO(II) = FOO(II) + EPS(I)
            I = I + 1
            II = II + (NO(IREP) + 1) * RCW
         ENDDO
         II = II - NO(IREP) * RCW
      ENDDO
C
      II = 1
      DO IREP = 1, NREP
         DO J = 1, NV(IREP)
            FVV(II) = FVV(II) + EPS(I)
            I = I + 1
            II = II + (NV(IREP) + 1) * RCW
         ENDDO
         II = II - NV(IREP) * RCW
      ENDDO
C
      IF (FOCKSP) THEN
C     
C     We also need to take the oa,oa and va,va elements in FVO into account
C     Copy them from FOO
C
         DO IREP = 1, NREP
C           The (oa,oa) elements
            IJOO = IOO(IREP)
            IJVO = IVO(IREP)
            IJOO = (IJOO + NO(IREP) * NIO(IREP) + NIO(IREP))*RCW + 1
            IJVO = (IJVO + NV(IREP) * NIO(IREP) + NSV(IREP))*RCW + 1
            DO J = NIO(IREP)+1, NSO(IREP)
               CALL XCOPY (NAO(IREP),FOO(IJOO),1,FVO(IJVO),1)
               IJOO = IJOO + NO(IREP) * RCW
               IJVO = IJVO + NV(IREP) * RCW
            ENDDO
C           The (va,va) elements
            IJOO = IOO(IREP)
            IJVO = IVO(IREP)
            IJOO = (IJOO + NO(IREP) * NSO(IREP) + NSO(IREP))*RCW + 1
            IJVO = (IJVO + NV(IREP) * NSO(IREP)            )*RCW + 1
            DO J = NSO(IREP)+1, NO(IREP)
               CALL XCOPY (NAV(IREP),FOO(IJOO),1,FVO(IJVO),1)
               IJOO = IJOO + NO(IREP) * RCW
               IJVO = IJVO + NV(IREP) * RCW
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE COPY_DIAG_TO_EPS (EPS,FOO,FVV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Copy Fock matrix diagonal into EPS array
C
C---------------Routines called----------------------------------------
C
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*),FOO(*),FVV(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "files.inc"
C
C---------------Local variables--------------------------------------
C
      integer i,ii,irep,j
C
C---------------Executable code--------------------------------------
C
C
      I = 1
      II = 1
      DO IREP = 1, NREP
         DO J = 1, NO(IREP)
            EPS(I) = FOO(II)
            I = I + 1
            II = II + (NO(IREP) + 1) * RCW
         ENDDO
         II = II - NO(IREP) * RCW
      ENDDO
C
      II = 1
      DO IREP = 1, NREP
         DO J = 1, NV(IREP)
            EPS(I) = FVV(II)
            I = I + 1
            II = II + (NV(IREP) + 1) * RCW
         ENDDO
         II = II - NV(IREP) * RCW
      ENDDO

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE REMOVE_F0 (EPS,FVO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Remove zeroth order Hamiltonian from Fock matrix by
C     subtracting orbital energies (EPS) at diagonal positions
C
C---------------Routines called----------------------------------------
C
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*),FVO(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "files.inc"
C
C---------------Local variables--------------------------------------
C
      integer i,ii,irep,j,k,kk
C
C---------------Executable code--------------------------------------
C
C
C     We substract the diagonal in the oa,oa and va,va elements of FVO
C
      I  = 0 ! index for 1D orbital energy (= zeroth order Fock matrix)
      DO IREP = 1, NREP
         DO J = 1, NO(IREP)
            I = I  + 1
            IF (J.LE.NIO(IREP)) CYCLE ! skip inactive occupied
            ! determine position of diagonal element in FVO
            IF (J.LE.NIO(IREP)+NAO(IREP)) THEN
              K = J-NIO(IREP)+NAV(IREP)+NIV(IREP) ! the oa,oa case
            ELSE
              K = J-NIO(IREP)-NAO(IREP)           ! the va,va case
            END IF
            II = (IVO(IREP) + NV(IREP)*(J-1) + K-1)*RCW + 1
            FVO(II) = FVO(II) - EPS(I)
         ENDDO
      ENDDO

      RETURN
      END
