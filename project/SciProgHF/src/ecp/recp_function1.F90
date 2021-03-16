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

MODULE RECP_FUNCTION1
CONTAINS

SUBROUTINE RECP_FUNCTION_FACAB (binom,An,Bn,crda,crdb,xab)
  IMPLICIT NONE
! variables : global
  INTEGER :: An, Bn
! dimension binom(*), crda(*), crdb(*), xab(*)
  REAL(8) :: binom,crda,crdb,xab
  DIMENSION binom(*), crda(*), crdb(*), xab(*)
! variables : local
  INTEGER :: I,Ai,Bi,SUM_A,SUM_B

  do I = 1, (An+Bn-1) 
     xab(I) = 0.d0
  enddo
 
  SUM_A=(An*(An-1))/2  ! sum from 1 to An
  SUM_B=(Bn*(Bn-1))/2
  do Ai=1,An
    do Bi=1,Bn
      xab((Ai-1)+Bi) = xab((Ai-1)+Bi)        &
           + binom(SUM_A+Ai)*crda((An+1)-Ai) &
           * binom(SUM_B+Bi)*crdb((Bn+1)-Bi)
    enddo
  enddo
END SUBROUTINE RECP_FUNCTION_FACAB


SUBROUTINE RECP_FN1_DFAC(ndfac)
! compute double factorials.
  USE RECP_IPT
  IMPLICIT NONE
  INTEGER ndfac,I
  REAL(8) TMP

  IPT_SH_DFAC(1)=1.0D0
  IPT_SH_DFAC(2)=1.0D0
  TMP=0.0D0
  DO I=1,ndfac-2
     TMP=TMP+1.0D0
     IPT_SH_DFAC(I+2) = TMP * IPT_SH_DFAC(I)
  ENDDO
END SUBROUTINE RECP_FN1_DFAC


SUBROUTINE RECP_SETZERO_I1(A,D1)
  IMPLICIT NONE
  INTEGER :: A(*)
  INTEGER :: I,D1
  DO I = 1,D1
     A(I) = 0
  ENDDO
END SUBROUTINE RECP_SETZERO_I1


SUBROUTINE RECP_SETZERO_I2(A,D1,D2)
  IMPLICIT NONE
  INTEGER :: I,J,D1,D2,A(D1,D2)
  DO I = 1,D1
     DO J = 1,D2
        A(I,J) = 0
     ENDDO
  ENDDO
END SUBROUTINE RECP_SETZERO_I2


SUBROUTINE RECP_SETZERO_I3(A,D1,D2,D3)
  IMPLICIT NONE
  INTEGER :: I,J,K,D1,D2,D3,A(D1,D2,D3)
  DO I = 1,D1
     DO J = 1,D2
        DO K = 1,D3
           A(I,J,K) = 0
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE RECP_SETZERO_I3


SUBROUTINE RECP_SETZERO_R1(A,D1)
  IMPLICIT NONE
  REAL(8) :: A(:)
  INTEGER :: I,D1
  DO I = 1,D1
     A(I) = 0.0d0
  ENDDO
END SUBROUTINE RECP_SETZERO_R1


END MODULE RECP_FUNCTION1
