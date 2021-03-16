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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* CUBIC_SPLINE_SETUP */
      SUBROUTINE CUBIC_SPLINE_SETUP(N,XA,YA,D2YA)
C=======================================================================
C     you provide:
C       N
C       XA = (x_1, x_2, ..., x_N)
C       YA = (y_1, y_2, ..., y_N)
C
C     routine gives back:
C       D2YA = (d2y_1, d2y_2, ..., d2y_N)
C       second derivatives at x_1, x_2, ..., x_N
C
C     natural boundary condition:
C       d2y_1 = d2y_N = 0
C
C     this routine is typically called once
C     and CUBIC_SPLINE_INTERPOLATION with the here obtained
C     D2YA can then be called as often as needed
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                           last revision: december 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 9999)

      DIMENSION XA(NMAX),
     &          YA(NMAX),
     &          D2YA(NMAX),
     &          TEMP(NMAX)

      D2YA(1) = D0
      D2YA(N) = D0

      TEMP(1) = D0

      DO I = 2,N-1
        A = (XA(I) - XA(I-1))/(XA(I+1) - XA(I-1))
        B = A*D2YA(I-1) + D2

        D2YA(I) = (A - D1)/B

        TEMP(I) = (YA(I+1) - YA(I))/(XA(I+1) - XA(I))
        TEMP(I) = TEMP(I) - (YA(I) - YA(I-1))/(XA(I) - XA(I-1))
        TEMP(I) = 6.0D0*TEMP(I)/(XA(I+1) - XA(I-1)) - A*TEMP(I-1)
        TEMP(I) = TEMP(I)/B
      ENDDO

      DO K = N-1,1,-1
        D2YA(K) = D2YA(K)*D2YA(K+1) + TEMP(K)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* CUBIC_SPLINE_INTERPOLATION */
      SUBROUTINE CUBIC_SPLINE_INTERPOLATION(N,XA,YA,D2YA,X,Y)
C=======================================================================
C     you provide:
C         N
C         XA = (x_1, x_2, ..., x_N)
C         YA = (y_1, y_2, ..., y_N)
C       D2YA = (d2y_1, d2y_2, ..., d2y_N)
C         X
C
C     routine gives back:
C         Y(X)
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                           last revision: december 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 9999)

      DIMENSION XA(NMAX),
     &          YA(NMAX),
     &          D2YA(NMAX)

      IS = 1
      IE = N

 1    CONTINUE
      IF(IE - IS .GT. 1) THEN
          K = (IE + IS)/2
        IF(XA(K) .GT. X) THEN
          IE = K
        ELSE
          IS = K
        ENDIF
        GOTO 1
      ENDIF

      H = XA(IE) - XA(IS)
      A = (XA(IE) - X)/H
      B = (X - XA(IS))/H
      Y = A*YA(IS) + B*YA(IE)
     &  + ((A**3 - A)*D2YA(IS) + (B**3 - B)*D2YA(IE))*(H**2)/6.0D0

      RETURN
      END
