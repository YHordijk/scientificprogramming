! mxmpack.F - Two old Cray library routines from the 1980s still called
!
      SUBROUTINE MXM(A,NAR,B,NAC,C,NBC)
! For the old Cray MXM library routine use optimized DGEMM
      implicit none
      INTEGER  NAR, NAC, NBC
      REAL*8   A(NAR,NAC), B(NAC,NBC), C(NAR,NBC)

      CALL DGEMM('N','N',NAR,NBC,NAC,1.0D0,A,NAR,B,NAC,0.0D0,C,NAR)
      RETURN
      END

      SUBROUTINE MXMA(A,ICA,IRA,B,ICB,IRB,C,ICC,IRC,NROW,NSUM,NCOL)
!
! General row and column spacing in A, B, and C.
! Can be used for e.g. C = A(transposed) * B
! or using a subblock of e.g. A.

#include "implicit.h"
      REAL*8  A(*),B(*),C(*)
      DO 100 IROW=0,NROW-1
      DO 100 ICOL=0,NCOL-1
         CSUM=0.0D0
         DO 110 ISUM=0,NSUM-1
            CSUM=CSUM +
     &          A(1+IROW*ICA+ISUM*IRA) *
     &          B(1+ISUM*ICB+ICOL*IRB)
110      CONTINUE
         INDEXC=1+ICC*IROW+IRC*ICOL
         C(INDEXC)=CSUM
100   CONTINUE
      RETURN
      END
! description:
!     CALL MXMA (A, NA, IAD, B, NB, IBD, C, NC, ICD, NAR, NAC, NBC)
!     where the arguments have the following meanings:
!     A first matrix of the product
!     NA spacing between column elements of A
!     IAD spacing between row elements of A
!     B second matrix of the product
!     NB spacing between column elements of B
!     IBD spacing between row elements of B
!     C output matrix
!     NC spacing between column elements of C
!     ICD spacing between row elements of C.
!     NAR number of rows in the first operand and result
!     NAC number of columns in the first operand and number of rows in the second
!     NBC number of columns in the second operand and result
