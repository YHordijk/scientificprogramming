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
      SUBROUTINE XAXPY (N,A,X,INCX,Y,INCY)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SAXPY, DAXPY, CAXPY, ZAXPY
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX,INCY
      REAL*8 A,X(*),Y(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CALL ZAXPY (N,A,X,INCX,Y,INCY)
C64B     CALL CAXPY (N,A,X,INCX,Y,INCY)
      ELSE
         CALL DAXPY (N,A,X,INCX,Y,INCY)
C64B     CALL SAXPY (N,A,X,INCX,Y,INCY)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOPY (N,X,INCX,Y,INCY)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SCOPY, DCOPY, CCOPY, ZCOPY
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX,INCY
      REAL*8 X(*),Y(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CALL ZCOPY (N,X,INCX,Y,INCY)
C64B     CALL CCOPY (N,X,INCX,Y,INCY)
      ELSE
         CALL DCOPY (N,X,INCX,Y,INCY)
C64B     CALL SCOPY (N,X,INCX,Y,INCY)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     &                  BETA,C,LDC)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SGEMM, DGEMM, CGEMM, ZGEMM
C
C         C:= ALPHA*op(A)*op(B) + BETA*C
C
C     where 
C        op(A) has dimension M x K
C        op(B) has dimension K x N
C        C     has dimension M x N
C
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*1 TRANSA,TRANSB
      INTEGER M,N,K,LDA,LDB,LDC
      REAL*8 ALPHA,BETA,A(*),B(*),C(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*1 TRANSAL,TRANSBL
C
C---------------Executable code--------------------------------------
C
      IF (M.EQ.0.OR.N.EQ.0.OR.K.EQ.0) RETURN
C
      IF (CARITH) THEN
         CALL ZGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     &               BETA,C,LDC)
C64B     CALL CGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
      ELSE
         TRANSAL = TRANSA
         TRANSBL = TRANSB
         IF (TRANSA.EQ.'C') TRANSAL = 'T'
         IF (TRANSB.EQ.'C') TRANSBL = 'T'
         IF (TRANSA.EQ.'c') TRANSAL = 't'
         IF (TRANSB.EQ.'c') TRANSBL = 't'
         CALL DGEMM (TRANSAL,TRANSBL,M,N,K,ALPHA,A,LDA,B,LDB,
     &               BETA,C,LDC)
C64B     CALL SGEMM (TRANSAL,TRANSBL,M,N,K,ALPHA,A,LDA,B,LDB,
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XGEMV (TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SGEMV, DGEMV, CGEMV, ZGEMV
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*1 TRANS
      INTEGER M,N,LDA,INCX,INCY
      REAL*8 ALPHA,BETA,A(*),X(*),Y(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*1 TRANSL
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CALL ZGEMV (TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
C64B     CALL CGEMV (TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      ELSE
         TRANSL = TRANS
         IF (TRANS.EQ.'C') TRANSL = 'T'
         IF (TRANS.EQ.'c') TRANSL = 't'
         CALL DGEMV (TRANSL,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
C64B     CALL SGEMV (TRANSL,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XSCAL (N,A,X,INCX)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SSCAL, DSCAL, CSCAL, ZSCAL
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX
      REAL*8 A,X(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CALL ZSCAL (N,A,X,INCX)
C64B     CALL CSCAL (N,A,X,INCX)
      ELSE
         CALL DSCAL (N,A,X,INCX)
C64B     CALL SSCAL (N,A,X,INCX)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XSUM (RESULT,N,X,INCX)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SSUM, DSUM, CSUM, ZSUM
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX
      REAL*8 X(*)
      REAL*8 RESULT(*),DSUM
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         RESULT(1) = DSUM (N,X(1),2*INCX) 
         RESULT(2) = DSUM (N,X(2),2*INCX) 
      ELSE
         RESULT(1) = DSUM (N,X,INCX) 
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#if defined (SYS_AIX)  || defined (SYS_IRIX) || defined (SYS_LINUX) \
 || defined (SYS_HPUX) || defined (SYS_T3D)  || defined (SYS_SUN) || defined (SYS_DARWIN)
C
      COMPLEX*16 FUNCTION ZSUM (N,X,INCX)
      COMPLEX*16 X(*),SUM
      SUM = DCMPLX(0.D0,0.D0)
      II = 1
      DO I = 1, N
         SUM = SUM + X(II)
         II = II + INCX
      ENDDO
      ZSUM = SUM
      RETURN
      END
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XTRANSP (M,N,AR,BR,AC,BC)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generic transpose routine.
C     Copies matrix A(M,N) into matrix B(N,M)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER M,N
      REAL*8 AR(M,N),BR(N,M)
      COMPLEX*16 AC(M,N),BC(N,M)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER I,J
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         DO I = 1, M
            DO J = 1, N
               BC(J,I) = AC(I,J)
            ENDDO
         ENDDO
      ELSE
         DO I = 1, M
            DO J = 1, N
               BR(J,I) = AR(I,J)
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XDOTC (RESULT,N,X,INCX,Y,INCY)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines DOTC and DOTU
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C     Note : The routine was written due to problems with the BLAS 
C            ZDOTC and ZDOTU functions, be careful if reintroducing them
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C     Date : 25-9-2003
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX,INCY
      REAL*8 X(*),Y(*)
      REAL*8 DDOT, RESULT(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         RESULT(1) = DDOT (N,X(1),2*INCX,Y(1),2*INCY) 
     &             + DDOT (N,X(2),2*INCX,Y(2),2*INCY)
         RESULT(2) = DDOT (N,X(1),2*INCX,Y(2),2*INCY) 
     &             - DDOT (N,X(2),2*INCX,Y(1),2*INCY)
      ELSE
         RESULT(1) = DDOT (N,X,INCX,Y,INCY)
      ENDIF
C
      RETURN
C
      ENTRY XDOTU (RESULT,N,X,INCX,Y,INCY)
C
      IF (CARITH) THEN
         RESULT(1) = DDOT (N,X(1),2*INCX,Y(1),2*INCY) 
     &             - DDOT (N,X(2),2*INCX,Y(2),2*INCY)
         RESULT(2) = DDOT (N,X(1),2*INCX,Y(2),2*INCY) 
     &             + DDOT (N,X(2),2*INCX,Y(1),2*INCY)
      ELSE
         RESULT(1) = DDOT (N,X,INCX,Y,INCY)
      ENDIF
C
      RETURN
C
      END
