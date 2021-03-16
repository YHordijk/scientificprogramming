      SUBROUTINE DIISX (NDIMT1,NDIMT2,T1N,T2N,IT,
     &                  CC,BB,ITC,BB2,RMS)
C
      IMPLICIT NONE
C
C---------------Description--------------------------------------------
C
C     DIIS extrapolation of the t1 and t2 amplitudes
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher (based on a routine by T.J. Lee)
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T2N(*),T1N(*)
      REAL*8 RMS
      INTEGER NDIMT1, NDIMT2, IT, ITC  
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "diis.inc"
#include "files.inc"
#include "waio.h"
#include "ccpar.inc"
C

      REAL*8 CC(MAXDIM),BB(MAXDIM+1,MAXDIM+2),BB2(MAXDIM+1,MAXDIM+2)
C---------------Local variables--------------------------------------
C
      REAL*8 XM,XMAX,XEN,XFAC,XADD,DET
      REAL*8 SUM2,SSUM,DSUM,SDOT,DDOT
      COMPLEX*16 CFAC
      real*8, allocatable :: T1O(:),T2O(:)
      integer it1,it2,jcode,m,n,sect1,sect2,intowp
C
C---------------Executable code--------------------------------------
C
      IF (MYPROC.NE.MASTER) RETURN
      allocate (T1O(RCW*NDIMT1))
      allocate (T2O(RCW*NDIMT2))

       T1O = 0.0d0
       T2O = 0.0d0
C
C     Handle special cases of initialization or finalization first
C
      IF (IT.EQ.-1) THEN
C
C       Initialization of the DIIS procedure;
C       compute record lengths, open files, write first set of amplitudes
C       This information is stored in the common block diis.inc
C
        NST1 = INTOWP(RCW*NDIMT1)/ NWORDS_ON_RECORD
        IF(NST1*NWORDS_ON_RECORD.NE.INTOWP(RCW*NDIMT1)) NST1 = NST1 + 1
        NST2 = INTOWP(RCW*NDIMT2)/ NWORDS_ON_RECORD
        IF(NST2*NWORDS_ON_RECORD.NE.INTOWP(RCW*NDIMT2)) NST2 = NST2 + 1
        CALL WAIO_OPEN(ITAMPL)
        CALL WAIO_OPEN(ITAMPD)
        SECT1 = MAXDIM*(NST1+NST2) + 1
        SECT2 = SECT1 + NST1
C       Place first input amplitudes at record maxdim+1 of ITAMPL
        CALL RWRIT(ITAMPL,T1N,INTOWP(NDIMT1*RCW),SECT1)
        CALL RWRIT(ITAMPL,T2N,INTOWP(NDIMT2*RCW),SECT2)
        IT = 0
        ITC = 0
        RMS = 1.0
C       Return to caller
        deallocate (T1O)
        deallocate (T2O)
        RETURN
      ENDIF
C
      IF (IT.EQ.-2) THEN
C
C        Stop DIIS and clean up files
C
        JCODE = 4
        CALL WAIO_CLOSE(ITAMPL,JCODE)
        CALL WAIO_CLOSE(ITAMPD,JCODE)
        RETURN
      ENDIF
C
C     Normal procedure : start by reading the old amplitudes that were used to produce T1N and T2N
C
      SECT1 = MAXDIM*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RREAD(ITAMPL,T1O,INTOWP(NDIMT1*RCW),SECT1)
      CALL RREAD(ITAMPL,T2O,INTOWP(NDIMT2*RCW),SECT2)
C
      IT = IT + 1
      IF(IT.GT.MAXDIM) IT = MAXDIM
      ITC = ITC + 1
      IF(ITC.GT.MAXDIM) ITC = 1
C
C     Calculate the difference between the input and output amplitudes
C     This is our error vector. Below we will consider the real and
C     imaginary differences as independent so that we can work with
C     a real error vector that is twice as long in complex as in real
C     cases.
C
      CALL XAXPY (NDIMT1,-A1,T1N,1,T1O,1)
      CALL XAXPY (NDIMT2,-A1,T2N,1,T2O,1)
      CALL XSCAL (NDIMT1,-A1,T1O,1)
      CALL XSCAL (NDIMT2,-A1,T2O,1)
C
      IF (CARITH) THEN
         SUM2 = DDOT (RCW*NDIMT1,T1O,1,T1O,1)
         SUM2 = SUM2 + DDOT (RCW*NDIMT2,T2O,1,T2O,1)
      ELSE
         SUM2 = DDOT (NDIMT1,T1O,1,T1O,1)
         SUM2 = SUM2 + DDOT (NDIMT2,T2O,1,T2O,1)
      ENDIF
C
      RMS = SQRT(SUM2/(NDIMT1+NDIMT2))
C
      SECT1 = (ITC-1)*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RWRIT(ITAMPD,T1O,INTOWP(NDIMT1*RCW),SECT1)
      CALL RWRIT(ITAMPL,T1N,INTOWP(NDIMT1*RCW),SECT1)
      CALL RWRIT(ITAMPD,T2O,INTOWP(NDIMT2*RCW),SECT2)
      CALL RWRIT(ITAMPL,T2N,INTOWP(NDIMT2*RCW),SECT2)
C
      DO 30 N = 1,IT
C
      SECT1 = (N-1)*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RREAD(ITAMPD,T1N,INTOWP(NDIMT1*RCW),SECT1)
      CALL RREAD(ITAMPD,T2N,INTOWP(NDIMT2*RCW),SECT2)
C
      IF (CARITH) THEN
         SUM2 = DDOT (RCW*NDIMT1,T1N,1,T1O,1)
         SUM2 = SUM2 + DDOT (RCW*NDIMT2,T2N,1,T2O,1)
      ELSE
         SUM2 = DDOT (NDIMT1,T1N,1,T1O,1)
         SUM2 = SUM2 + DDOT (NDIMT2,T2N,1,T2O,1)
      ENDIF
C
      BB2(N,ITC) = SUM2
      BB2(ITC,N) = BB2(N,ITC)
C
   30 CONTINUE
C
C  TRANSFER INTO THE BB ARRAY
C
      DO 60 N = 1,IT
      DO 70 M = 1,N-1
      BB(N,M) = BB2(N,M)
      BB(M,N) = BB2(M,N)
  70  CONTINUE
      BB(N,N) = BB2(N,N)
  60  CONTINUE
C
C  FIND THE MAXIMUM AND SCALE
C
      XM =  ABS(BB(1,1))
      DO 80 N = 1,IT
      DO 90 M = 1,N
      XFAC =  ABS(BB(N,M))
      XMAX =  DMAX1(XM,XFAC)
      XM = XMAX
   90 CONTINUE
   80 CONTINUE
C
      XM = AR1/XM
      DO 100 N = 1,IT
      DO 110 M = 1,IT
      BB(N,M) = BB(N,M)*XM
  110 CONTINUE
  100 CONTINUE
C
      IT1 = IT+1
      IT2 = IT+2
      DO 120 N = 1,IT
      BB(N,IT1) = -AR1
      BB(IT1,N) = -AR1
  120 CONTINUE
C
      BB(IT1,IT1) = AR0
      DO 130 N = 1,IT
      BB(N,IT2) = AR0
  130 CONTINUE
      BB(IT1,IT2) = -AR1
C
      CALL FLIN(BB,MAXDIM+1,IT1,1,DET)
      XADD = AR0
      DO 140 N = 1,IT
      CC(N) = BB(N,IT2)
      XADD = XADD + CC(N)
  140 CONTINUE
C
      CALL XCOPY(NDIMT1,A0,0,T1N,1)
      CALL XCOPY(NDIMT2,A0,0,T2N,1)
C
      DO 150 N = 1,IT
C
      SECT1 = (N-1)*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RREAD(ITAMPL,T1O,INTOWP(NDIMT1*RCW),SECT1)
      CALL RREAD(ITAMPL,T2O,INTOWP(NDIMT2*RCW),SECT2)
C
      CFAC = DCMPLX(CC(N),AR0)
      CALL XAXPY (NDIMT1,CFAC,T1O,1,T1N,1)
      CALL XAXPY (NDIMT2,CFAC,T2O,1,T2N,1)
C
  150 CONTINUE
C
C     Write extrapolated amplitudes (which become input in the next iteration)
C     at the end of the amplitude file.
C
      SECT1 = MAXDIM*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RWRIT(ITAMPL,T1N,INTOWP(NDIMT1*RCW),SECT1)
      CALL RWRIT(ITAMPL,T2N,INTOWP(NDIMT2*RCW),SECT2)
      deallocate (T1O)
      deallocate (T2O)
C
      RETURN
      END
