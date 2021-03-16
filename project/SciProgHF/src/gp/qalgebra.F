!define QDIAG_TIMING
!define QDIAG_OMP_INFO
C  /* Deck iqpack */
      SUBROUTINE IQPACK(IQPA,NZA,IQPB,NZB,IQPC,NZC)
C***********************************************************************
C
C     IQPA and IQPB are pointers to quaternion units 1,i,j and k
C     IQPACK considers what quaternion units span the quaternion product
C     A * B and makes packs these in pointer IQPC as generated.
C
C     Written by T.Saue - June 27 1996
C     Last revision: June 27 1996 - tsaue
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
C
      LOGICAL LBUF(4)
      DIMENSION IQPA(NZA),IQPB(NZB),IQPC(*)
C
C     Initialize LBUF
C     ===============
C
      LBUF(1) = .FALSE.
      LBUF(2) = .FALSE.
      LBUF(3) = .FALSE.
      LBUF(4) = .FALSE.
C
C     Determine quaternion units that contribute to A*B
C     =================================================
C
      NZC = 0
      DO IZA = 1,NZA
        IQA = IQPA(IZA)
        DO IZB = 1,NZB
          IQB = IQPB(IZB)
          IQC = IQMULT(IQA,IQB,1)
          IF(.NOT.LBUF(IQC)) THEN
            NZC = NZC + 1
            IQPC(NZC) = IQC
            LBUF(IQC) = .TRUE.
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qgemm */
      SUBROUTINE QGEMM(M,N,K,ALPHA,FA,TA,IQPA,A,LRA,LCA,NZA,
     &                             FB,TB,IQPB,B,LRB,LCB,NZB,
     &                              BETA,IQPC,C,LRC,LCC,NZC)
C*****************************************************************************
C
C     Performs matrix-matrix operations on general matrices:
C
C     C := alpha * op(A) * op(B) + beta * C
C
C     where alpha and beta are real constants. A,B and C
C     are matrices that may be real(NZ = 1), complex(NZ = 2)
C     or quaternion (NZ = 4). Calls BLAS routine DGEMM.
C
C     LRX and LCX are leading rows and columns of matrix X.
C     op(A) is an M by K matrix, op(B) is a K by N matrix
C     and C an M by N matrix.
C
C     op(X) is determined by the corresponding character variable
C     FX:
C       FX = 'N'        Normal matrix
C       FX = 'T'        Transpose..
C       FX = 'C'        Complex conjugate
C       FX = 'H'        Hermitian conjugate
C     and TX
C       TX = 'N'        No transform    a+bj
C       TX = 'I'        i-transform   -i(a+bj)i = a - bj
C       TX = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TX = 'K'        k-transform   -k(a+bj)k = a*-b*j
C
C     AFUL and BFUL indicates on input what components of
C     matrices A and B are non-zero, whereas CFUL on output
C     indicates what components of C are non-zero.
C
C     Written by T.Saue, University of Oslo, Nov. 1994
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
C
C     Global variables
C
      DIMENSION A(LRA*LCA,NZA),B(LRB*LCB,NZB),C(LRC*LCC,NZC),
     &          IQPA(NZA),IQPB(NZB),IQPC(NZC)
      CHARACTER  FA*1,FB*1,TA*1,TB*1,OA*1,OB*1
C
C     Local variables
C
      DIMENSION IPB(4)
      LOGICAL   LBUF(4)      
C
C
C     Initialize LBUF
C     ===============
C
      LBUF(1:4) = .FALSE.
C
C     Make pointer IPB
C     ================
C
      CALL IZERO(IPB,4)
      DO IZB = 1,NZB
        IPB(IQPB(IZB)) = IZB
      ENDDO
C
C     MATRIX A
C

C     Determine FX:
C     =============
C     Normal matrix...
      IF    (FA.EQ.'N') THEN
        OA  = 'N'
        IFA = 1
C     Transpose...
      ELSEIF(FA.EQ.'T') THEN
        OA  = 'T'
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        OA  = 'N'
        IFA = 2
C     Hermitian conjugate...
      ELSEIF(FA.EQ.'H') THEN
        OA  = 'T'
        IFA = 2
      ELSE
        CALL QUIT('QGEMM:Unknown FA '//FA//' of matrix A!')
      ENDIF
C     Determine TX:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QGEMM:Unknown TA '//TA//' of matrix A!')
      ENDIF
C
C     MATRIX B
C
C     Determine FX:
C     =============
C     Normal matrix...
      IF    (FB.EQ.'N') THEN
        OB  = 'N'
        IFB = 1
C     Transpose...
      ELSEIF(FB.EQ.'T') THEN
        OB  = 'T'
        IFB = 1
C     Complex conjugate....
      ELSEIF(FB.EQ.'C') THEN
        OB  = 'N'
        IFB = 2
C     Hermitian conjugate...
      ELSEIF(FB.EQ.'H') THEN
        OB  = 'T'
        IFB = 2
      ELSE
        CALL QUIT('QGEMM:Unknown format '//FA//' of matrix B!')
      ENDIF
C     Determine TX:
C     =============
C     No transform
      IF    (TB.EQ.'N') THEN
        ITB = 1
C     i-transform
      ELSEIF(TB.EQ.'I') THEN
        ITB = 2
C     j-transform
      ELSEIF(TB.EQ.'J') THEN
        ITB = 3
C     k-transform
      ELSEIF(TB.EQ.'K') THEN
        ITB = 4
      ELSE
        CALL QUIT('QGEMM:Unknown TB '//TB//' of matrix B!')
      ENDIF
C
      DO 10 IZC = 1,NZC
        IQC = IQPC(IZC)
        DO 20 IZA = 1,NZA
          IQA = IQPA(IZA)
          IQB = IQMULT(IQA,IQC,1)
          IZB = IPB(IQB)
          IF(IZB.EQ.0) GOTO 20
          FACA = ALPHA*(IQSIGN(IQA,IFA,ITA)*IQSIGN(IQB,IFB,ITB)
     &                *IQPHASE(IQA,IQB,1))
          IF(LBUF(IZC)) THEN
            FACB = D1
          ELSE
            FACB = BETA
          ENDIF
          CALL DGEMM(OA,OB,M,N,K,FACA,A(1,IZA),LRA,B(1,IZB),LRB,
     &               FACB,C(1,IZC),LRC)
          LBUF(IZC) = .TRUE.
   20   CONTINUE
   10 CONTINUE
C
C     Check that matrix C is filled
C
      DO IZC = 1,NZC
      IF(.NOT.LBUF(IZC)) 
     &    WRITE(LUPRI,'(A,I5)') '* WARNING(QGEMM): Not filled: ',IZC
      ENDDO
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qgeold */
      SUBROUTINE QGEOLD(M,N,K,ALPHA,FA,TA,IQPA,A,LRA,LCA,NZA,
     &                             FB,TB,IQPB,B,LRB,LCB,NZB,
     &                              BETA,IQPC,C,LRC,LCC,NZC)
C*****************************************************************************
C
C     Performs matrix-matrix operations on general matrices:
C
C     C := alpha * op(A) * op(B) + beta * C
C
C     where alpha and beta are real constants. A,B and C
C     are matrices that may be real(NZ = 1), complex(NZ = 2)
C     or quaternion (NZ = 4). Calls BLAS routine DGEMM.
C
C     LRX and LCX are leading rows and columns of matrix X.
C     op(A) is an M by K matrix, op(B) is a K by N matrix
C     and C an M by N matrix.
C
C     op(X) is determined by the corresponding character variable
C     FX:
C       FX = 'N'        Normal matrix
C       FX = 'T'        Transpose..
C       FX = 'C'        Complex conjugate
C       FX = 'H'        Hermitian conjugate
C     and TX
C       TX = 'N'        No transform    a+bj
C       TX = 'I'        i-transform   -i(a+bj)i = a - bj
C       TX = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TX = 'K'        k-transform   -k(a+bj)k = a*-b*j
C
C     AFUL and BFUL indicates on input what components of
C     matrices A and B are non-zero, whereas CFUL on output
C     indicates what components of C are non-zero.
C
C     Written by T.Saue, University of Oslo, Nov. 1994
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
C
C     Global variables
C
      DIMENSION A(LRA*LCA,NZA),B(LRB*LCB,NZB),C(LRC*LCC,NZC),
     &          IQPA(NZA),IQPB(NZB),IQPC(NZC)
      CHARACTER  FA*1,FB*1,TA*1,TB*1,OA*1,OB*1
C
C     Local variables
C
      LOGICAL   LBUF(4)      
      DIMENSION IPC(4)
C
C
C
C     Initialize LBUF
C     ===============
C
      LBUF(1:4) = .FALSE.
C
C     Make pointer IPC
C     ================
C
      CALL IZERO(IPC,4)
      DO IZC = 1,NZC
        IPC(IQPC(IZC)) = IZC
      ENDDO
      IF(NZA.GT.NZC) CALL QUIT('QGEMM:NZA.GT.NZC')
      IF(NZB.GT.NZC) CALL QUIT('QGEMM:NZB.GT.NZC')
C
C     MATRIX A
C

C     Determine FX:
C     =============
C     Normal matrix...
      IF    (FA.EQ.'N') THEN
        OA  = 'N'
        IFA = 1
C     Transpose...
      ELSEIF(FA.EQ.'T') THEN
        OA  = 'T'
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        OA  = 'N'
        IFA = 2
C     Hermitian conjugate...
      ELSEIF(FA.EQ.'H') THEN
        OA  = 'T'
        IFA = 2
      ELSE
        CALL QUIT('QGEMM:Unknown FA '//FA//' of matrix A!')
      ENDIF
C     Determine TX:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QGEMM:Unknown TA '//TA//' of matrix A!')
      ENDIF
C
C     MATRIX B
C
C     Determine FX:
C     =============
C     Normal matrix...
      IF    (FB.EQ.'N') THEN
        OB  = 'N'
        IFB = 1
C     Transpose...
      ELSEIF(FB.EQ.'T') THEN
        OB  = 'T'
        IFB = 1
C     Complex conjugate....
      ELSEIF(FB.EQ.'C') THEN
        OB  = 'N'
        IFB = 2
C     Hermitian conjugate...
      ELSEIF(FB.EQ.'H') THEN
        OB  = 'T'
        IFB = 2
      ELSE
        CALL QUIT('QGEMM:Unknown format '//FA//' of matrix B!')
      ENDIF
C     Determine TX:
C     =============
C     No transform
      IF    (TB.EQ.'N') THEN
        ITB = 1
C     i-transform
      ELSEIF(TB.EQ.'I') THEN
        ITB = 2
C     j-transform
      ELSEIF(TB.EQ.'J') THEN
        ITB = 3
C     k-transform
      ELSEIF(TB.EQ.'K') THEN
        ITB = 4
      ELSE
        CALL QUIT('QGEMM:Unknown TB '//TB//' of matrix B!')
      ENDIF
C
      DO 10 IZA = 1,NZA
        IQA = IQPA(IZA)
        DO 20 IZB = 1,NZB
          IQB = IQPB(IZB)
          IQC = IQMULT(IQA,IQB,1)
          IZC = IPC(IQC)
          IF(IZC.EQ.0) THEN
            WRITE(LUPRI,'(A)') 'QGEMM: Error in IPC !'
            WRITE(LUPRI,'(A,4I5)') 'IQPC: ',(IQPC(I),I=1,NZC)
            WRITE(LUPRI,'(A,4I5)') 'IPC : ',(IPC(I), I=1,4)
            CALL QUIT('QGEMM: Error in IPC !')
          ENDIF
          FACA = ALPHA*(IQSIGN(IQA,IFA,ITA)*IQSIGN(IQB,IFB,ITB)
     &                *IQPHASE(IQA,IQB,1))
          IF(LBUF(IZC)) THEN
            FACB = D1
          ELSE
            FACB = BETA
          ENDIF
          CALL DGEMM(OA,OB,M,N,K,FACA,A(1,IZA),LRA,B(1,IZB),LRB,
     &               FACB,C(1,IZC),LRC)
          LBUF(IZC) = .TRUE.
   20   CONTINUE
   10 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qdiag */
      SUBROUTINE QDIAG(NZ,N,A,LRA,LCA,EIG,MATZ,VEC,LRV,LCV,
     &                 WORK,LWORK,IERR)
C***********************************************************************
C
C       This is a driver routine for the diagonalization of
C         * real symmetric matrices                     NZ = 1
C         * complex Hermitian matrices                  NZ = 2
C         * quaternion Hermitian matrices               NZ = 4
C
C       INPUT:
C
C       A       matrix to be diagonalized
C
C       Control parameter MATZ:
C               MATZ = 0        Only eigenvalues desired.
C               MATZ = 1        Eigenvalues and eigenvectors
C
C       OUTPUT:
C
C       EIG     eigenvalues in ascending order.
C       VEC     eigenvectors when MATZ=1.
C
C       Error parameter IERR:
C               IERR = 0        Normal completion
C               IERR.NE.0       Erroneous completion
C
C       TEMPORARY STORAGE:      FV1,FV2,TAU
C       Written by T.Saue November 1994 - Odense
C       Last revision : Nov 15 1994 - tsaue
C***********************************************************************
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(LRA,LCA,NZ),EIG(N),VEC(LRV,LCV,NZ),WORK(LWORK)
      CHARACTER SECTID*12,CPUTID*12,WALLTID*12
      CALL QENTER('QDIAG')
#include "memint.h"
      IF(MATZ.EQ.1.AND.(NZ.EQ.1.OR.NZ.EQ.2).AND.(LRA.NE.LRV)) THEN
        CALL QUIT('QDIAG: LRA.ne.LRV')
      ENDIF

#ifdef LAPACK_QDIAG
C     Use the fast DSYEVR for real case
      IF (NZ.EQ.1.AND.MATZ.EQ.1) THEN
            DDUM = 0.0D0
            IDUM = 0
            WORKDUM = 0.0D0
            IWORKDUM = 0
C     eigval and eigvec
c     first let lapack estimate memory use
CMI ...  I have to initialize it for intel mkl-lapack
            INFO = 0
            IERR = 0
            CALL DSYEVR('V','A','U',N,A,LRA,DDUM,DDUM,IDUM,IDUM,0.0D0,
     &           M,EIG,VEC,LRV,IDUM,WORKDUM,-1,IWORKDUM,-1,INFO)
            LLWORK = NINT(WORKDUM)
            LILWORK = IWORKDUM
            CALL MEMGET('REAL',KLWORK, LLWORK, WORK, KFREE, LFREE)
            CALL MEMGET('INTE',KLIWORK, LILWORK, WORK, KFREE, LFREE)
            CALL MEMGET('INTE',KISUPPZ, 2*N, WORK, KFREE, LFREE)
            CALL DSYEVR('V','A','U',N,A,LRA,DDUM,DDUM,IDUM,IDUM,0.0D0,
     &           M,EIG,VEC,LRV,WORK(KISUPPZ),WORK(KLWORK),LLWORK,
     &           WORK(KLIWORK),LILWORK,IERR)
C     DSYEVR does not guarantee a sorting order, but I have never seen it
C     do anything other than ascending, so the sorting should be a noop.
            CALL ORDER3(VEC,EIG,LRV,N,N,1)
            CALL QEXIT('QDIAG')
            RETURN
      ENDIF
#endif
C
C     Memory allocation
C
      CALL MEMGET('REAL',KFV1,N   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFV2,N   ,WORK,KFREE,LFREE)
      IF(NZ.GT.1) THEN
        NTAU = N*NZ
        CALL MEMGET('REAL',KTAU,NTAU,WORK,KFREE,LFREE)
      ENDIF
C
C       Householder transform to real symmetric tridiagonal matrix
C       ==========================================================
C
      CALL GETTIM(CPU1,WALL1)
      IF (NZ.EQ.4) THEN
        CALL  QTRIDI_omp(A,N,LRA,LCA,EIG,
     &               WORK(KFV1),WORK(KFV2),WORK(KTAU))
        IF(MATZ.NE.0) GO TO 10
      ELSEIF(NZ.EQ.2) THEN
        CALL  HTRIDI_omp(A,N,LRA,LCA,EIG,
     &               WORK(KFV1),WORK(KFV2),WORK(KTAU))
        IF(MATZ.NE.0) GO TO 10
      ELSEIF(NZ.EQ.1) THEN
        IF(MATZ.EQ.0) THEN
          CALL  TRED1(LRA,N,A(1,1,1),EIG,WORK(KFV1),WORK(KFV2))
        ELSE
          CALL  TRED2(LRA,N,A(1,1,1),EIG,WORK(KFV1),VEC(1,1,1))
          GOTO 20
        ENDIF
      ELSE
        CALL QUIT('QDIAG: Illegal value of NZ !')
      ENDIF
C
C       FIND EIGENVALUES ONLY
C       =====================
C
      CALL GETTIM(CPU2,WALL2)
      CALL TQLRAT_omp(N,EIG,WORK(KFV2),IERR)
      CALL GETTIM(CPU3,WALL3)
      GO TO 30
C
C       FIND BOTH EIGENVALUES AND EIGENVECTORS
C       ======================================
C
   10 CONTINUE
      CALL DUNIT2(VEC,N,LRV,LCV,NZ)
C
C       Find eigenvalues/vectors of real symmetric
C       tridiagonal matrix
C
   20 CONTINUE
      CALL GETTIM(CPU2,WALL2)
      CALL TQL2_omp(LRV,N,EIG,WORK(KFV1),VEC(1,1,1),IERR)
      CALL GETTIM(CPU3,WALL3)
      IF (IERR.NE.0) GO TO 30
C
C       Backtransform to find eigenvectors of
C       quaternion Hermitian matrix
C
      IF    (NZ.EQ.4) THEN
        CALL  QTRIBK_omp(A,N,LRA,LCA,WORK(KTAU),VEC,N,LRV,LCV)
      ELSEIF(NZ.EQ.2) THEN
        CALL  HTRIBK_omp(A,N,LRA,LCA,WORK(KTAU),VEC,N,LRV,LCV)
      ENDIF
   30 CONTINUE

#ifdef QDIAG_TIMING
         WRITE(LUPRI,'(/A,I2,2X,I0)') 'QDIAG timings. NZ,N =',NZ,N
         CALL GETTIM(CPU4,WALL4)
         CPUTID = SECTID(CPU2-CPU1)
         WALLTID = SECTID(WALL2-WALL1)
         WRITE(LUPRI,9020) CPUTID,WALLTID
         CPUTID = SECTID(CPU3-CPU2)
         WALLTID = SECTID(WALL3-WALL2)
         WRITE(LUPRI,9021) CPUTID,WALLTID
         CPUTID = SECTID(CPU4-CPU3)
         WALLTID = SECTID(WALL4-WALL3)
         WRITE(LUPRI,9022) CPUTID,WALLTID
         flush(LUPRI)

 9020 FORMAT(/'QDIAG: QTRIDI CPU (WALL) time: ',A12,' (',A12,')')
 9021 FORMAT( 'QDIAG: TQL2   CPU (WALL) time: ',A12,' (',A12,')')
 9022 FORMAT( 'QDIAG: QTRIBK CPU (WALL) time: ',A12,' (',A12,')')
#endif
C
C     Memory deallocation
C
      CALL MEMREL('QDIAG.QDIAG',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('QDIAG')
      RETURN
C
      END

      SUBROUTINE QDIAG90(NZ,N,A,LRA,LCA,EIG,MATZ,VEC,LRV,LCV,IERR)
C     Like QDIAG but uses f90 memory allocation internally.
      use memory_allocator
#include "implicit.h"
      real(8), allocatable :: work(:)

      lenid = 3 ! lenid used in MEMGET in QDIAG
      lwork = 4*lenid + 2*N
      if(nz.gt.1) lwork = lwork + N*NZ

#ifdef LAPACK_QDIAG
      if(matz.eq.1)then
        DDUM = 0.0D0
        IDUM = 0
        WORKDUM  = 0.0D0
        IWORKDUM = 0
c       first let lapack estimate memory use
        INFO = 0
        CALL DSYEVR('V','A','U',N,A,LRA,DDUM,DDUM,IDUM,IDUM,0.0D0,
     &              M,EIG,VEC,LRV,IDUM,WORKDUM,-1,IWORKDUM,-1,INFO)
        LLWORK  = NINT(WORKDUM)
        LILWORK = IWORKDUM
        lwork = LLWORK + LILWORK + 2*N + 6*lenid
      end if
#endif
      if(lwork.le.0)then
        call quit(' *** error in QDIAG90: buffer space allocation not
     & possible for length <= 0.***')
      end if
      call alloc(work,lwork, id='QDIAG90-buffer-space')
      work = 0
      call qdiag(NZ,N,A,LRA,LCA,EIG,MATZ,VEC,LRV,LCV,work,lwork,IERR)
      call dealloc(work, id='QDIAG90-buffer-space')
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qmgath */
      SUBROUTINE QMGATH(AIN,NRW1,NCL1,AOUT,NRW2,NCL2,NZ,
     &                  NROW,NCOL,IROW,ICOL)
C***********************************************************************
C
C     Gather operations for generally quaternion matrices
C
C     Written by T.Saue November 1994
C     Last revision: Nov 15 1994
C***********************************************************************
#include "implicit.h"
      DIMENSION AIN(NRW1,NCL1,NZ),AOUT(NRW2,NCL2,NZ)
      DIMENSION IROW(*),ICOL(*)
C
      DO 10 IZ = 1,NZ
        DO 20 J = 1, NCOL
          JJ = ICOL(J)
          DO 30 I = 1, NROW
           AOUT(I,J,IZ) = AIN(IROW(I),JJ,IZ)
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QTRIDI_omp(A,N,LRA,LCA,D,E,E2,TAU)
C
      IMPLICIT NONE
      INTEGER ::  N, LRA, LCA
      REAL*8  ::  A(LRA,LCA,4),D(N),E(N),E2(N),TAU(4,N)
C***********************************************************************
C

C       Trond Saue, University of Tromsoe, Norway , December 1993
C
C       Rewritten by Hans Joergen Aa. Jensen, December 2017.
C       - optimized loop structure (the inner loops were good for algol (or C), but bad in Fortran)
C       - some additional restructuring for performance
C       - implemented openMP parallelization
C       (The old EISPACK routines went through a matrix A(i,j) with j in
C       the inner loop, good for algol (or C), but bad in Fortran.)
C
C       This subroutine is a quaternion analogue of the
C       EISPACK routine HTRIDI. QTRIDI reduces a quaternion
C       Hermitian matrix to a real symmetric tridiagonal matrix.
C       The unitary transformation is done in two steps:
C
C       1. Reduction to quaternion Hermitian tridiagonal
C       matrix using (N-2) Householder matrices.
C
C       2.Reduction to real symmetric tridiagonal matrix
C       using a quaternion diagonal matrix TAU.
C
C       Householder matrices has the general form:
C
C       P = 1 - (uu(+))/H       ;       H = u(+)u/2
C
C       INPUT:
C
C       [AR+iAI+jAJ+kAK]        is a quaternion Hermitian matrix of
C                               order N and with leading dimension LRA.
C
C       OUTPUT:
C
C       The lower triangle of A contains the row vectors u
C       from the (n-2) Householder matrices. The strict upper
C       triangle and the diagonal of A is unaltered.
C
C       D(i)    diagonal elements T(i,i) of the real tridiagonal matrix
C       E(i)    off-diagonal elements T(i+1,i) of the real tridiagonal matrix
C       E2      squares of elements in E
C       TAU     the quaternion diagonal of TAU
C     ------------------------------------------------------------------
C
      INTEGER  ::  I, IM1, J, JM1, K, L
      REAL*8   ::  TR, TI, TJ, TK, FR, FI, FJ, FK, GR, GI, GJ, GK, FRR
      REAL*8   ::  H, HH, D_I, FHH, FR_INV, ASCALE, ASCALE_INV
      REAL*8,  PARAMETER :: D0 = 0.0D0, D1 = 1.0D0

#ifdef QDIAG_OMP_INFO
      integer   :: omp_get_num_threads, omp_get_thread_num
#endif

      TAU(1,N) = D1
      TAU(2,N) = D0
      TAU(3,N) = D0
      TAU(4,N) = D0
C
C     Store diagonal elements of A temporarily in D
C
      DO I = 1, N
         D(I) = A(I,I,1)
      END DO
C
C    Reduce lower triangle of A
C
      DO 300 I = N, 2, -1
         IM1 = I - 1

         FR = A(IM1,I,1)**2 + A(IM1,I,2)**2
     &      + A(IM1,I,3)**2 + A(IM1,I,4)**2

         HH = D0
         DO K = 1, I-2
            HH = HH+A(K,I,1)*A(K,I,1)+A(K,I,2)*A(K,I,2)
     +             +A(K,I,3)*A(K,I,3)+A(K,I,4)*A(K,I,4)
         END DO

         H  = HH + FR ! norm2 of A(1:i-1,i)

         IF (H .EQ. D0) THEN ! matrix is block diagonal, starting new Householder for next block
            TAU(1,IM1) = D1
            TAU(2,IM1) = D0
            TAU(3,IM1) = D0
            TAU(4,IM1) = D0
            E(I)     = D0
            E2(I)    = D0
            D_I        = D(I)
            D(I)       = A(I,I,1) ! modified diagonal
            A(I,I,1)   = D_I      ! original diagonal
            A(I,I,2)   = D0       ! = HH
            GO TO 300
         END IF

         FR    = SQRT(FR)
         GR    = SQRT(H)
         E(I)  = GR
         E2(I) = H
C
C     Form next diagonal element of matrix T: TAU(I+1) = + TAU(I)*A(I+1,I)/F
C     where F is the absolute value of A(I+1,I)
C
C     TAU(I+1) is stored in TR, TI, TJ, TK because TAU is used
C     for other numbers below. It is set at the end of the DO I loop.
C
         IF (FR.EQ.D0) THEN
           ! special case for |A(I+1,I)| = 0
           ! We choose positive real value: U(I+1) = GR,
           ! in this case TAU(I+1) = -TAU(I)
            TR  =  -TAU(1,I)
            TI  =  -TAU(2,I)
            TJ  =  -TAU(3,I)
            TK  =  -TAU(4,I)
            A(IM1,I,1)  = GR
            GO TO 170
         END IF

         FR_INV = D1/FR
         TR = -FR_INV*( TAU(1,I)*A(IM1,I,1)+TAU(2,I)*A(IM1,I,2)
     +                 +TAU(3,I)*A(IM1,I,3)+TAU(4,I)*A(IM1,I,4))
         TI = -FR_INV*(-TAU(1,I)*A(IM1,I,2)+TAU(2,I)*A(IM1,I,1)
     +                 -TAU(3,I)*A(IM1,I,4)+TAU(4,I)*A(IM1,I,3))
         TJ = -FR_INV*(-TAU(1,I)*A(IM1,I,3)+TAU(2,I)*A(IM1,I,4)
     +                 +TAU(3,I)*A(IM1,I,1)-TAU(4,I)*A(IM1,I,2))
         TK = -FR_INV*(-TAU(1,I)*A(IM1,I,4)-TAU(2,I)*A(IM1,I,3)
     +                 +TAU(3,I)*A(IM1,I,2)+TAU(4,I)*A(IM1,I,1))

C
C     Form element IM1 in vector U
C        U(IM1) = A(IM1,I) + A(IM1,I)[SQRT(H)/QVAL(A(IM1,I))]
C
         IF (IM1 .EQ. 1) THEN ! no Householder transformation for last element
            TR = -TR
            TI = -TI
            TJ = -TJ
            TK = -TK
            GO TO 270
         END IF
         ! IF (HH .EQ. 0.0D0) GO TO 270 ! no reason to do housholder if already tridiagonal

         H          =  H  + FR    *GR ! eq. (D.15)
         GR         =  D1 + FR_INV*GR ! eq. (D.15)
         A(IM1,I,1) =  GR * A(IM1,I,1)
         A(IM1,I,2) =  GR * A(IM1,I,2)
         A(IM1,I,3) =  GR * A(IM1,I,3)
         A(IM1,I,4) =  GR * A(IM1,I,4)

  170    CONTINUE ! goto 170 if FR.eq.0.0d0 above
C
         FR = D0
!$OMP PARALLEL IF(IM1 .GT. 10) DEFAULT(PRIVATE)
!$OMP& SHARED(I,IM1,H,FR,A,TAU,N,FHH)
#ifdef QDIAG_OMP_INFO
         if (I.eq.N) 
     &   print *,'QTRIDI_omp, openMP threads',omp_get_num_threads()
     &   ,  omp_get_thread_num()
#endif
!$OMP DO SCHEDULE(GUIDED) REDUCTION(+:FR)
         DO J = 1,IM1
C
C     Form element j of A*U
C       (A*U)(j) is stored in (GR,GI,GJ,GK)
C       (A*U)(j) = A(J,K) * U(K) = A*(K,J) * A(K,I)
C       as A(K,I) contains U
C
            GR = D0
            GI = D0
            GJ = D0
            GK = D0
            DO K = 1,IM1 ! use full A(K,J), i.e. do not assume only OK in upper half
               GR = GR+A(K,J,1)*A(K,I,1)+A(K,J,2)*A(K,I,2)
     +                +A(K,J,3)*A(K,I,3)+A(K,J,4)*A(K,I,4)
               GI = GI+A(K,J,1)*A(K,I,2)-A(K,J,2)*A(K,I,1)
     +                -A(K,J,3)*A(K,I,4)+A(K,J,4)*A(K,I,3)
               GJ = GJ+A(K,J,1)*A(K,I,3)+A(K,J,2)*A(K,I,4)
     +                -A(K,J,3)*A(K,I,1)-A(K,J,4)*A(K,I,2)
               GK = GK+A(K,J,1)*A(K,I,4)-A(K,J,2)*A(K,I,3)
     +                +A(K,J,3)*A(K,I,2)-A(K,J,4)*A(K,I,1)
            END DO ! K
C
C
C     Form element j of P
C        P(j) = (A*U)(j)/H
C             - stored in (TAU(1,j),TAU(2,j),TAU(3,j),TAU(4,j))
C
            TAU(1,J) = GR/H ! eq. (D.32)
            TAU(2,J) = GI/H
            TAU(3,J) = GJ/H
            TAU(4,J) = GK/H
            FR       = FR+A(J,I,1)*TAU(1,J)+A(J,I,2)*TAU(2,J) ! +=u^H(j) * p(j)
     +                   +A(J,I,3)*TAU(3,J)+A(J,I,4)*TAU(4,J)
         END DO ! J = 1,IM1
!$OMP END DO
C
C     Evaluate FHH = SUM[j](U*(j)P(j )/(2H))
C
!$OMP SINGLE
         FHH=FR/(H+H) ! K in eq. (D.35)
!$OMP END SINGLE
C
C     Form reduced A
C
C     Computational formula: A'(k,j) = A(k,j)-Q(k)U*(j)-U(k)Q*(j)
C        where Q(j) = P(j) - FHH*U(j)
C     Q(k) is stored in TAU(1:4,k), k = 1:i-1
!$OMP DO SCHEDULE(GUIDED)
         DO J = 1,IM1
            TAU(1,J) =  TAU(1,J) - FHH*A(J,I,1)
            TAU(2,J) =  TAU(2,J) - FHH*A(J,I,2)
            TAU(3,J) =  TAU(3,J) - FHH*A(J,I,3)
            TAU(4,J) =  TAU(4,J) - FHH*A(J,I,4)
         END DO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED)
         DO J = 1,IM1
C
C     Q(k) is stored in TAU(1:4,k)
C     U(k) is stored in (AR(k,i),AI(k,i),AJ(k,i),AK(k,i))
C
C     U*(j) is stored in (FR,FI,FJ,FK)
C     Q*(j) is stored in (GR,GI,GJ,GK)
C
            FRR      =  A(J,I,1)
            FI       = -A(J,I,2)
            FJ       = -A(J,I,3)
            FK       = -A(J,I,4)
            GR       =  TAU(1,J)
            GI       = -TAU(2,J)
            GJ       = -TAU(3,J)
            GK       = -TAU(4,J)
C
C     Computational formula: A'(k,j) = A(k,j) - Q(k)U*(j) - U(k)Q*(j)
C            = A(k,j) - TAU(k).F - A*(k,i).G
C
            DO K = 1,IM1
              A(K,J,1) = A(K,J,1)
     +            -FRR*TAU(1,K)+FI*TAU(2,K)+FJ*TAU(3,K)+FK*TAU(4,K)
     +            -GR*A(K,I,1)+GI*A(K,I,2)+GJ*A(K,I,3)+GK*A(K,I,4)
              A(K,J,2) = A(K,J,2)
     +            -FRR*TAU(2,K)-FI*TAU(1,K)+FJ*TAU(4,K)-FK*TAU(3,K)
     +            -GR*A(K,I,2)-GI*A(K,I,1)+GJ*A(K,I,4)-GK*A(K,I,3)
              A(K,J,3) = A(K,J,3)
     +            -FRR*TAU(3,K)-FI*TAU(4,K)-FJ*TAU(1,K)+FK*TAU(2,K)
     +            -GR*A(K,I,3)-GI*A(K,I,4)-GJ*A(K,I,1)+GK*A(K,I,2)
              A(K,J,4) = A(K,J,4)
     +            -FRR*TAU(4,K)+FI*TAU(3,K)-FJ*TAU(2,K)-FK*TAU(1,K)
     +            -GR*A(K,I,4)+GI*A(K,I,3)-GJ*A(K,I,2)-GK*A(K,I,1)
            END DO ! K
         END DO ! J
!$OMP END DO
!$OMP END PARALLEL
C
  270    CONTINUE
C
         TAU(1,IM1) = TR
         TAU(2,IM1) = TI
         TAU(3,IM1) = TJ
         TAU(4,IM1) = TK

         D_I     = D(I)
         D(I)    = A(I,I,1)
         A(I,I,1) = D_I
         A(I,I,2) = H
  300 CONTINUE ! DO 300 I = 1,N

         E(1)    = D0
         E2(1)   = D0
         D_I     = D(1)
         D(1)    = A(1,1,1)
         A(1,1,1) = D_I
         A(1,1,2) = 0.0D0
C
      RETURN
      END
      SUBROUTINE QTRIBK_omp(A,N,LRA,LCA, TAU, Z,M,LRZ,LCZ)
C
      IMPLICIT NONE
      INTEGER  :: N, LRA, LCA, M, LRZ,LCZ
      REAL*8   :: A(LRA,LCA,4),TAU(4,N),Z(LRZ,LCZ,4)
C***********************************************************************
C
C       Trond Saue, University of Tromsoe, Norway , December 1993
C
C       Rewritten by Hans Joergen Aa. Jensen, December 2017.
C       - optimized loop structure (the inner loops were good for algol (or C), but bad in Fortran)
C       - some additional restructuring for performance
C       - implemented openMP parallelization
C
C       This subroutine is a quaternion analogue of the
C       EISPACK routine HTRIBK. It forms the eigenvectors
C       of a quaternion Hermitian matrix by backtransforming
C       those of the corresponding real symmetric matrix
C       determined by QTRIDI.
C
C       INPUT:
C
C       [AR+iAI+jAJ+kAK]        is a quaternion Hermitian matrix of
C                               order N and with leading dimension LRA.
C                               The lower triangle contains the
C                               u-vectors of the (N-2) Householder
C                               matrices used in the reduction of
C                               A to a quaternion Hermitian tridiagonal
C                               matrix. The diagonal and upper triangle
C                               of the original matrix is preserved.
C       TAU                     TAU is the diagonal quaternion elements
C                               of the diagonal matrix transforming
C                               the quaternion Hermitian tridiagonal
C                               matrix to a real tridiagonal matrix.
C       OUTPUT:
C
C       [ZR+iZI+jZJ+kZK]        eigenvectors
C
C     Note that the last component of each returned vector
C     is real and that vector Euclidean norms are preserved.
C
C     ------------------------------------------------------------------
C
      INTEGER   :: I, J, K
      REAL*8    :: H, SR, SI, SJ, SK
#ifdef QDIAG_OMP_INFO
      integer   :: omp_get_num_threads, omp_get_thread_num
#endif

C
C     Transform the eigenvector of the real symmetric tridiagonal matrix
C     to those of the quaternion Hermitian tridiagonal matrix
C
C      Z = TAU* . Z
C
!$OMP PARALLEL DEFAULT(NONE) SHARED(M,N,TAU,A,Z)
!$OMP& PRIVATE(I,J,K,H,SR,SI,SJ,SK)
#ifdef QDIAG_OMP_INFO
      print '(A,3I10)','QTRIBK_omp - openMP threads, my rank',
     &   omp_get_num_threads(), omp_get_thread_num(), M
#endif
!$OMP DO SCHEDULE(GUIDED)
      DO J = 1, M
#ifdef QDIAG_OMP_INFO
!     if (MOD(J, 99) .eq. 1) then
!        print *, 'J, rank',J, omp_get_thread_num()
!     end if
#endif

         DO K = 1, N
            Z(K,J,4) = -TAU(4,K)*Z(K,J,1)
            Z(K,J,3) = -TAU(3,K)*Z(K,J,1)
            Z(K,J,2) = -TAU(2,K)*Z(K,J,1)
            Z(K,J,1) =  TAU(1,K)*Z(K,J,1)
         END DO
C
C     Recover and apply the Householder matrices
C
         DO I = 3,N
            H  = A(I,I,2)
         IF (H .EQ. 0.0D0) CYCLE ! see QTRIDI_omp
C
         ! s(j) = sum(k) u*(k) Z(k,j) / H
            SR = 0.0D0
            SI = 0.0D0
            SJ = 0.0D0
            SK = 0.0D0
            DO 110 K = 1,I-1
               SR = SR+A(K,I,1)*Z(K,J,1)+A(K,I,2)*Z(K,J,2)
     +                +A(K,I,3)*Z(K,J,3)+A(K,I,4)*Z(K,J,4)
               SI = SI+A(K,I,1)*Z(K,J,2)-A(K,I,2)*Z(K,J,1)
     +                -A(K,I,3)*Z(K,J,4)+A(K,I,4)*Z(K,J,3)
               SJ = SJ+A(K,I,1)*Z(K,J,3)+A(K,I,2)*Z(K,J,4)
     +                -A(K,I,3)*Z(K,J,1)-A(K,I,4)*Z(K,J,2)
               SK = SK+A(K,I,1)*Z(K,J,4)-A(K,I,2)*Z(K,J,3)
     +                +A(K,I,3)*Z(K,J,2)-A(K,I,4)*Z(K,J,1)
  110       CONTINUE
C
            SR = SR/H
            SI = SI/H
            SJ = SJ/H
            SK = SK/H
C
C     Do second part
         ! Z(k,j) = Z(k,j) - u(k) . s(j)
C
            DO 120 K = 1,I-1
               Z(K,J,1) = Z(K,J,1)-SR*A(K,I,1)+SI*A(K,I,2)
     +                            +SJ*A(K,I,3)+SK*A(K,I,4)
               Z(K,J,2) = Z(K,J,2)-SR*A(K,I,2)-SI*A(K,I,1)
     +                            +SJ*A(K,I,4)-SK*A(K,I,3)
               Z(K,J,3) = Z(K,J,3)-SR*A(K,I,3)-SI*A(K,I,4)
     +                            -SJ*A(K,I,1)+SK*A(K,I,2)
               Z(K,J,4) = Z(K,J,4)-SR*A(K,I,4)+SI*A(K,I,3)
     +                            -SJ*A(K,I,2)-SK*A(K,I,1)
  120       CONTINUE
         END DO ! I = 3,N
      END DO ! J = 1,M
!$OMP END DO
!$OMP END PARALLEL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HTRIDI_omp(A,N,LRA,LCA,D,E,E2,TAU)
C
      IMPLICIT NONE
      INTEGER ::  N, LRA, LCA
      REAL*8  ::  A(LRA,LCA,2),D(N),E(N),E2(N),TAU(2,N)
C***********************************************************************
C
C       Written by Hans Joergen Aa. Jensen, October 2018.
C
C       This subroutine is an optimized version of EISPACK routine HTRIDI.
C       HTRIDI reduces a complex Hermitian matrix to a real symmetric tridiagonal matrix.
C       It was created by stripping the previously optimized QTRIDI_omp for the j and k imaginary parts.
C       Optimizations compared to HTRIDI:
C       - optimized loop structure (the inner loops were good for algol (or C), but bad in Fortran)
C       - some additional restructuring for performance
C       - implemented openMP parallelization
C       (The old EISPACK routines went through a matrix A(i,j) with j in
C       the inner loop, good for algol (or C), but bad in Fortran.)
C
C       1. Reduction to complex Hermitian tridiagonal
C       matrix using (N-2) Householder matrices.
C
C       2.Reduction to real symmetric tridiagonal matrix
C       using a complex diagonal matrix TAU.
C
C       Householder matrices has the general form:
C
C       P = 1 - (uu(+))/H       ;       H = u(+)u/2
C
C       INPUT:
C
C       [AR+iAI]        is a complex Hermitian matrix of
C                       order N and with leading dimension LRA.
C
C       OUTPUT:
C
C       The A matrix is destroyed and contains the row vectors u
C       from the (n-2) Householder matrices.
C
C       D(i)    diagonal elements T(i,i) of the real tridiagonal matrix
C       E(i)    off-diagonal elements T(i+1,i) of the real tridiagonal matrix
C       E2      squares of elements in E
C       TAU     the complex diagonal of TAU
C     ------------------------------------------------------------------
C
      INTEGER  ::  I, IM1, J, JM1, K, L
      REAL*8   ::  TR, TI, TJ, TK, FR, FI, FJ, FK, GR, GI, GJ, GK, FRR
      REAL*8   ::  H, HH, D_I, FHH, FR_INV, ASCALE, ASCALE_INV
      REAL*8,  PARAMETER :: D0 = 0.0D0, D1 = 1.0D0

#ifdef QDIAG_OMP_INFO
      integer   :: omp_get_num_threads, omp_get_thread_num
#endif

      TAU(1,N) = D1
      TAU(2,N) = D0
C
C     Store diagonal elements of A temporarily in D
C
      DO I = 1, N
         D(I) = A(I,I,1)
      END DO
C
C    Reduce lower triangle of A
C
      DO 300 I = N, 2, -1
         IM1 = I - 1

         FR = A(IM1,I,1)**2 + A(IM1,I,2)**2

         HH = D0
         DO K = 1, I-2
            HH = HH+A(K,I,1)*A(K,I,1)+A(K,I,2)*A(K,I,2)
         END DO

         H  = HH + FR ! norm2 of A(1:i-1,i)

         IF (H .EQ. D0) THEN ! matrix is block diagonal, starting new Householder for next block
            TAU(1,IM1) = D1
            TAU(2,IM1) = D0
            E(I)     = D0
            E2(I)    = D0
            D_I        = D(I)
            D(I)       = A(I,I,1) ! modified diagonal
            A(I,I,1)   = D_I      ! original diagonal
            A(I,I,2)   = D0       ! = HH
            GO TO 300
         END IF

         FR    = SQRT(FR)
         GR    = SQRT(H)
         E(I)  = GR
         E2(I) = H
C
C     Form next diagonal element of matrix T: TAU(I+1) = + TAU(I)*A(I+1,I)/F
C     where F is the absolute value of A(I+1,I)
C
C     TAU(I+1) is stored in TR, TI, TJ, TK because TAU is used
C     for other numbers below. It is set at the end of the DO I loop.
C
         IF (FR.EQ.D0) THEN
           ! special case for |A(I+1,I)| = 0
           ! We choose positive real value: U(I+1) = GR,
           ! in this case TAU(I+1) = -TAU(I)
            TR  =  -TAU(1,I)
            TI  =  -TAU(2,I)
            A(IM1,I,1)  = GR
            GO TO 170
         END IF

         FR_INV = D1/FR
         TR = -FR_INV*( TAU(1,I)*A(IM1,I,1)+TAU(2,I)*A(IM1,I,2) )
         TI = -FR_INV*(-TAU(1,I)*A(IM1,I,2)+TAU(2,I)*A(IM1,I,1) )

C
C     Form element IM1 in vector U
C        U(IM1) = A(IM1,I) + A(IM1,I)[SQRT(H)/QVAL(A(IM1,I))]
C
         IF (IM1 .EQ. 1) THEN ! no Householder transformation for last element
            TR = -TR
            TI = -TI
            GO TO 270
         END IF
         ! IF (HH .EQ. 0.0D0) GO TO 270 ! no reason to do housholder if already tridiagonal

         H          =  H  + FR    *GR ! eq. (D.15)
         GR         =  D1 + FR_INV*GR ! eq. (D.15)
         A(IM1,I,1) =  GR * A(IM1,I,1)
         A(IM1,I,2) =  GR * A(IM1,I,2)

  170    CONTINUE ! goto 170 if FR.eq.0.0d0 above
C
         FR = D0
!$OMP PARALLEL IF(IM1 .GT. 10) DEFAULT(PRIVATE)
!$OMP& SHARED(I,IM1,H,FR,A,TAU,N,FHH)
#ifdef QDIAG_OMP_INFO
         if (I.eq.N) 
     &   print *,'HTRIDI_omp, openMP threads',omp_get_num_threads()
     &   ,  omp_get_thread_num()
#endif
!$OMP DO SCHEDULE(GUIDED) REDUCTION(+:FR)
         DO J = 1,IM1
C
C     Form element j of A*U
C       (A*U)(j) is stored in (GR,GI,GJ,GK)
C       (A*U)(j) = A(J,K) * U(K) = A*(K,J) * A(K,I)
C       as A(K,I) contains U
C
            GR = D0
            GI = D0
            DO K = 1,IM1 ! use full A(K,J), i.e. do not assume only OK in upper half
               GR = GR+A(K,J,1)*A(K,I,1)+A(K,J,2)*A(K,I,2)
               GI = GI+A(K,J,1)*A(K,I,2)-A(K,J,2)*A(K,I,1)
            END DO ! K
C
C
C     Form element j of P
C        P(j) = (A*U)(j)/H
C             - stored in (TAU(1,j),TAU(2,j),TAU(3,j),TAU(4,j))
C
            TAU(1,J) = GR/H ! eq. (D.32)
            TAU(2,J) = GI/H
            FR       = FR+A(J,I,1)*TAU(1,J)+A(J,I,2)*TAU(2,J) ! +=u^H(j) * p(j)
         END DO ! J = 1,IM1
!$OMP END DO
C
C     Evaluate FHH = SUM[j](U*(j)P(j )/(2H))
C
!$OMP SINGLE
         FHH=FR/(H+H) ! K in eq. (D.35)
!$OMP END SINGLE
C
C     Form reduced A
C
C     Computational formula: A'(k,j) = A(k,j)-Q(k)U*(j)-U(k)Q*(j)
C        where Q(j) = P(j) - FHH*U(j)
C     Q(k) is stored in TAU(1:4,k), k = 1:i-1
!$OMP DO SCHEDULE(GUIDED)
         DO J = 1,IM1
            TAU(1,J) =  TAU(1,J) - FHH*A(J,I,1)
            TAU(2,J) =  TAU(2,J) - FHH*A(J,I,2)
         END DO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED)
         DO J = 1,IM1
C
C     Q(k) is stored in TAU(1:4,k)
C     U(k) is stored in (AR(k,i),AI(k,i),AJ(k,i),AK(k,i))
C
C     U*(j) is stored in (FR,FI,FJ,FK)
C     Q*(j) is stored in (GR,GI,GJ,GK)
C
            FRR      =  A(J,I,1)
            FI       = -A(J,I,2)
            GR       =  TAU(1,J)
            GI       = -TAU(2,J)
C
C     Computational formula: A'(k,j) = A(k,j) - Q(k)U*(j) - U(k)Q*(j)
C            = A(k,j) - TAU(k).F - A*(k,i).G
C
            DO K = 1,IM1
              A(K,J,1) = A(K,J,1)
     +            -FRR*TAU(1,K)+FI*TAU(2,K)
     +            -GR*A(K,I,1)+GI*A(K,I,2)
              A(K,J,2) = A(K,J,2)
     +            -FRR*TAU(2,K)-FI*TAU(1,K)
     +            -GR*A(K,I,2)-GI*A(K,I,1)
            END DO ! K
         END DO ! J
!$OMP END DO
!$OMP END PARALLEL
C
  270    CONTINUE
C
         TAU(1,IM1) = TR
         TAU(2,IM1) = TI

         D_I      = D(I)
         D(I)     = A(I,I,1)
         A(I,I,1) = D_I
         A(I,I,2) = H
  300 CONTINUE ! DO 300 I = 1,N

         E(1)     = D0
         E2(1)    = D0
         D_I      = D(1)
         D(1)     = A(1,1,1)
         A(1,1,1) = D_I
         A(1,1,2) = D0
C
      RETURN
      END
      SUBROUTINE HTRIBK_omp(A,N,LRA,LCA, TAU, Z,M,LRZ,LCZ)
C
      IMPLICIT NONE
      INTEGER  :: N, LRA, LCA, M, LRZ,LCZ
      REAL*8   :: A(LRA,LCA,2),TAU(2,N),Z(LRZ,LCZ,2)
C***********************************************************************
C
C
C       Written by Hans Joergen Aa. Jensen, October 2018.
C
C       This subroutine is an optimized version of EISPACK routine HTRIBK
C       It forms the eigenvectors of a complex Hermitian matrix by backtransforming
C       those of the corresponding real symmetric matrix determined by HTRIDI_omp.
C
C       It was created by stripping the previously optimized QTRIBK_omp for the j and k imaginary parts.
C       Optimizations compared to HTRIBK:
C       - optimized loop structure (the inner loops were good for algol (or C), but bad in Fortran)
C       - some additional restructuring for performance
C       - implemented openMP parallelization
C       (The old EISPACK routines went through a matrix A(i,j) with j in
C       the inner loop, good for algol (or C), but bad in Fortran.)
C
C       INPUT:
C
C       [AR+iAJ]          is a complex Hermitian matrix of
C                         order N and with leading dimension LRA.
C                         It contains the u-vectors of the (N-2) Householder
C                         matrices used in the reduction of A
C                         to a complex Hermitian tridiagonal matrix.
C       TAU               TAU is the diagonal complex elements
C                         of the diagonal matrix transforming
C                         the complex Hermitian tridiagonal
C                         matrix to a real tridiagonal matrix.
C       OUTPUT:
C
C       [ZR+iZI]          eigenvectors
C
C     Note that the last component of each returned vector
C     is real and that vector Euclidean norms are preserved.
C
C     ------------------------------------------------------------------
C
      INTEGER   :: I, J, K
      REAL*8    :: H, SR, SI, SJ, SK
#ifdef QDIAG_OMP_INFO
      integer   :: omp_get_num_threads, omp_get_thread_num
#endif

C
C     Transform the eigenvector of the real symmetric tridiagonal matrix
C     to those of the complex Hermitian tridiagonal matrix
C
C      Z = TAU* . Z
C
!$OMP PARALLEL DEFAULT(NONE) SHARED(M,N,TAU,A,Z)
!$OMP& PRIVATE(I,J,K,H,SR,SI,SJ,SK)
#ifdef QDIAG_OMP_INFO
      print '(A,3I10)','HTRIBK_omp - openMP threads, my rank',
     &   omp_get_num_threads(), omp_get_thread_num(), M
#endif
!$OMP DO SCHEDULE(GUIDED)
      DO J = 1, M

         DO K = 1, N
            Z(K,J,2) = -TAU(2,K)*Z(K,J,1)
            Z(K,J,1) =  TAU(1,K)*Z(K,J,1)
         END DO
C
C     Recover and apply the Householder matrices
C
         DO I = 3,N
            H  = A(I,I,2)
         IF (H .EQ. 0.0D0) CYCLE ! see HTRIDI_omp
C
         ! s(j) = sum(k) u*(k) Z(k,j) / H
            SR = 0.0D0
            SI = 0.0D0
            DO 110 K = 1,I-1
               SR = SR+A(K,I,1)*Z(K,J,1)+A(K,I,2)*Z(K,J,2)
               SI = SI+A(K,I,1)*Z(K,J,2)-A(K,I,2)*Z(K,J,1)
  110       CONTINUE
C
            SR = SR/H
            SI = SI/H
C
C     Do second part
         ! Z(k,j) = Z(k,j) - u(k) . s(j)
C
            DO 120 K = 1,I-1
               Z(K,J,1) = Z(K,J,1)-SR*A(K,I,1)+SI*A(K,I,2)
               Z(K,J,2) = Z(K,J,2)-SR*A(K,I,2)-SI*A(K,I,1)
  120       CONTINUE
         END DO ! I = 3,N
      END DO ! J = 1,M
!$OMP END DO
!$OMP END PARALLEL
C
      RETURN
      END
      SUBROUTINE TQL2_omp(NM,N,D,E,Z,IERR)
C
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL*8  D(N),E(N),Z(NM,N)
      REAL*8  C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1.
C
C        E HAS BEEN DESTROYED.
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0D0
      TST1 = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0
         H = DABS(D(L)) + DABS(E(L))
         IF (TST1 .LT. H) TST1 = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + DABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = SQRT(1.0D0 + P*P)
         D(L) = E(L) / (P + DSIGN(R,P))
         D(L1) = E(L) * (P + DSIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = SQRT(P**2 + E(I)**2)
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + DABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE TQLRAT_omp(N,D,E2,IERR)
C
      IMPLICIT NONE
      INTEGER  I,J,L,M,N,II,L1,MML,IERR
      REAL*8   D(N),E2(N)
      REAL*8   B,C,F,G,H,P,R,S,T
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E2 HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0D0
      T = 0.0D0
      E2(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
         H = DABS(D(L)) + DSQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = 1.D-14*T ! was EPSLON(T)
         C = B * B
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = SQRT(1.0D0 + P**2)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0D0) G = B
            H = G * P / R
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qopvec */
      SUBROUTINE QOPVEC(FA,TA,IQPA,A,N,NZA)
C***********************************************************************
C
C     Modify quaternion vector as follows:
C
C     FA:
C       FA = 'N'        Normal array
C       FA = 'C'        Complex conjugate
C     and TA
C       TA = 'N'        No transform    a+bj
C       TA = 'I'        i-transform   -i(a+bj)i = a - bj
C       TA = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TA = 'K'        k-transform   -k(a+bj)k = a*-b*j
C
C     Written by T.Saue July 1996
C     Last revision: July 14 196 - tsaue
C     
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
      CHARACTER FA*1, TA*1
      DIMENSION A(N,NZA),IQPA(*)
C
C
C     Determine FA:
C     =============
C     Normal ARRAY...
      IF    (FA.EQ.'N') THEN
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        IFA = 2
      ELSE
        CALL QUIT('QVECOP:Unknown FA '//FA//' of matrix A!')
      ENDIF
C
C     Determine TA:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QVECOP:Unknown TA '//TA//' of matrix A!')
      ENDIF
C
      DO IZA = 1,NZA
        IQA = IQPA(IZA)
        FACA = IQSIGN(IQA,IFA,ITA)
        IF(FACA.NE.D1) CALL DSCAL(N,FACA,A(1,IZA),1)
      ENDDO      
C
      RETURN
      END
! end of qalgebra.F
