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
C
      SUBROUTINE ZZERO (LV,VECTOR)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
      REAL*8 VECTOR (LV)
C
      DO 100 I = 1, LV
      VECTOR (I) = 0.D0
  100 CONTINUE
C
      RETURN
      END
C     ========================================
      SUBROUTINE Q2JOB(ICODE,STRING,IPAR,IVAL)
C     ========================================
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     THIS ROUTINE (QUIT TO JOB) FINISHES THE CALLING PROGRAM AFTER
C     ISSUING AN APPROPRIATE ERROR-MESSAGE. THE FOLLOWING OPTIONS
C     ARE AVAILABLE :
C     1---END OF FILE ENCOUNTERED
C         STRING---DEFENITION OF THE RECORD EXPECTED
C     2---USER ERROR (AS INCONSISTENT INPUT)
C         STRING---INDICATION OF THE ERROR
C     3---DIMENSION BOUND EXCEEDED BECAUSE PROGRAM PARAMETERS ARE
C         TOO SMALL TO RUN WITH THE INPUT-DATA
C         STRING---INDICATION OF SECTION WHERE ERROR WAS DETECTED
C         IPAR  ---NAME OF THE PARAMETER TO BE CHANGED
C         IVAL  ---INDICATION OF THE LEAST VALUE TO ASSIGN TO PARAMETER
C     4---DIMENSION BOUND EXCEEDED BECAUSE OF INTERNAL PROGRAM ERRORS
C         STRING---INDICATION OF THE LOCATION WHERE ERROR WAS DETECTED
C     5---NON MATCHING PARAMETERS
C         STRING---PLACE ERROR OCCURRED
C         IPAR  ---PARAMETER IN ERROR
C         IVAL  ---PREVIOUS VALUE FOR PARAMETER
C     6---ANY PROGRAM ERROR
C         STRING---INDICATION/LOCATION OF ERROR
C
      CHARACTER STRING*(*),AST*4,IPAR*6
      DATA AST/' ** '/
      PRINT 100,AST
      GOTO (10,20,30,40,50,60) ICODE
   10 PRINT 400,AST,'>>>>>EOF ENCOUNTERED : ',STRING,' EXPECTED <<<<<'
      PRINT 200,AST
      CALL QUIT('--CHECK INPUT DECK FOR MISSING DATA--')
   20 PRINT 200,AST,STRING,AST
      PRINT 200,AST
      CALL QUIT('** USER ERROR--PROGRAM ABORTED **')
   30 PRINT 200,AST,'DIMENSION BOUND EXCEEDED :',AST
      WRITE(6,*) 'Change ',IPAR,' to: ',IVAL
      PRINT 200,AST,'ERROR OCCURRED AT ',STRING,AST
      PRINT 200,AST,'PARAMETER TO BE CHANGED : ',IPAR,AST
      PRINT 300,AST,'TO AT LEAST ',IVAL,AST
      PRINT 200,AST
      CALL QUIT('-- ENLARGE PARAMETER  --')
   40 PRINT 200,AST,'DIMENSION BOUND EXCEEDED :',AST
      PRINT 200,AST//'ERROR OCCURRED AT ',STRING,AST
      PRINT 200,AST,'PROGRAM ERROR',AST
      PRINT 200,AST
      CALL QUIT('-- CORRECT PROGRAM --')
   50 PRINT 400,AST,'NON MATCHING PARAMETER AT ',STRING,AST
      PRINT 300,AST,'PARAMETER '//IPAR//' SHOULD BE AT LEAST',IVAL,AST
      PRINT 200,AST
      CALL QUIT('--NON MATCHING PARAMETERS--')
   60 PRINT 200,AST
      PRINT 200,AST,STRING,AST
      PRINT 200,AST,'PROGRAM ERROR',AST
      PRINT 200,AST
      CALL QUIT('--PROGRAM ERROR--')
  100 FORMAT(/A)
  200 FORMAT(2A,T55,A)
  300 FORMAT(2A,I10,T55,A)
  400 FORMAT(4A)
      END
C
      SUBROUTINE CRDIAG (NM,N,AR,AI,W,ZR,ZI,FV1,FV2,FM1,IERR)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     DIAGONALIZE COMPLEX,HERMITIAN MATRIX BY CALLING APPROPRIATE LIBRARY
C     Checks for reality of matrix and calls real diagonalizer if possible.
C
C     1) CALL EISPACK LIBRARY
C
C     AR, AI  : MATRIX TO BE DIAGONALIZED ON INPUT; EIGENVECTORS ON OUTPUT
C     ZR, ZI  : WORKING SPACE
C     W       : EIGENVALUES
C     FV1,FM1 : WORKING SPACE
C
      DIMENSION AR(NM,NM),AI(NM,NM),W(NM),ZR(NM,NM),ZI(NM,NM)
      DIMENSION FV1(NM*NM),FM1(NM*NM*2),FV2(NM*NM)

C
      DO I = 1, N
      DO J = 1, N
         IF (ABS(AI(J,I)).GT.1.D-12) GOTO 11
      ENDDO
      ENDDO
      CALL RS (NM,N,AR,W,1,AI,FV1,FM1,IERR)
      DO I = 1, N
      DO J = 1, N
         AR(J,I) = AI(J,I)
         AI(J,I) = 0.D0
      ENDDO
      ENDDO
      RETURN
   11 CONTINUE
      CALL CH (NM,N,AR,AI,W,1,ZR,ZI,FV1,FV2,FM1,IERR)
      DO 10 I = 1, N
      DO 10 J = 1, N
         AR(I,J) = ZR(I,J)
         AI(I,J) = ZI(I,J)
   10 CONTINUE
C
      RETURN
      END
