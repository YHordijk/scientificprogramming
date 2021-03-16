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

MODULE RECP_OUTPUT
CONTAINS

! ---------------------------------------------------------------

SUBROUTINE OUTPUT_TITLE
! OUTPUT : print title
  IMPLICIT NONE
#include "inc_print.h"
  WRITE(RECP_OUT,'(/)')
  WRITE(RECP_OUT,'(10X,A)')'-----------------------------------------------------------'
  WRITE(RECP_OUT,'(10X,A)')'Relativistic Effective Core Potential Integral routine'
  WRITE(RECP_OUT,'(10X,A)')'This routine is based on ARGOS integral'
  WRITE(RECP_OUT,'(10X,A)')'with the permission of R. M. Pitzer (Ohio State University)'
  WRITE(RECP_OUT,'(10X,A)')'and maintained by Y. S. Lee and Y. C. Park (KAIST)' 
  WRITE(RECP_OUT,'(10X,A)')''
  WRITE(RECP_OUT,'(10X,A)')'Electronic mail:'
  WRITE(RECP_OUT,'(10X,A)')'YoonSupLee@kaist.ac.kr'
  WRITE(RECP_OUT,'(10X,A)')'youngc@kaist.ac.kr'
  WRITE(RECP_OUT,'(10X,A)')''
  WRITE(RECP_OUT,'(10X,A)')'Reference:'
  WRITE(RECP_OUT,'(10X,A)')'Bull.Korean.Chem.Soc. v.33, p.803 (2012)'
  WRITE(RECP_OUT,'(10X,A)')'-----------------------------------------------------------'
  WRITE(RECP_OUT,'(/)')
END SUBROUTINE OUTPUT_TITLE

! ---------------------------------------------------------------

SUBROUTINE OUTPUT_LOC(RECP_LOC,RECP_STATUS)
!  OUTPUT : print location
   IMPLICIT NONE
#include "inc_print.h"
   CHARACTER(LEN=*) :: RECP_LOC
   CHARACTER*1  :: RECP_STATUS
   IF (RECP_DBG.GE.5) THEN
      IF(RECP_STATUS .EQ. 'E') WRITE(RECP_OUT,*) ' - RECP : Entering.... ', RECP_LOC
      IF(RECP_STATUS .EQ. 'X') WRITE(RECP_OUT,*) ' - RECP : Exiting..... ', RECP_LOC
      IF(RECP_STATUS .EQ. 'H') WRITE(RECP_OUT,*) ' - RECP : Here........ ', RECP_LOC
   ENDIF
END SUBROUTINE OUTPUT_LOC

! ---------------------------------------------------------------

SUBROUTINE OUTPUT_ERROR(RECP_LOC,ERR_MESSAGE,ERR_VALUE,ERR_TYPE)
  IMPLICIT NONE
#include "inc_print.h"
  CHARACTER(LEN=*) :: RECP_LOC, ERR_MESSAGE
  INTEGER          :: ERR_VALUE, ERR_TYPE
  IF (ERR_TYPE.EQ.3) THEN    
     WRITE(RECP_OUT,*) RECP_LOC,':',ERR_MESSAGE,':',ERR_VALUE
     STOP
  ENDIF
END SUBROUTINE OUTPUT_ERROR

! ---------------------------------------------------------------

SUBROUTINE OUTPUT_DBG_TITLE(DBG_TITLE)
  IMPLICIT NONE
#include "inc_print.h"
  CHARACTER(LEN=*) DBG_TITLE
  INTEGER  I
  IF (RECP_DBG.GE.2) THEN
     WRITE(RECP_OUT,'(/,3X,A)') DBG_TITLE
     WRITE(RECP_OUT,'(3X,100A)') ('-',I=1,LEN(DBG_TITLE))      
  ENDIF
END SUBROUTINE OUTPUT_DBG_TITLE

! ---------------------------------------------------------------

SUBROUTINE OUTPUT_DBG_CHAR1(TITLE0,DBG_LEVEL)
  IMPLICIT NONE
#include "inc_print.h"
  CHARACTER(LEN=*) TITLE0
  INTEGER  DBG_LEVEL
  IF (RECP_DBG.GE.DBG_LEVEL) THEN
     WRITE(RECP_OUT,'(/,A,A)') ' * ',TITLE0
  ENDIF
END SUBROUTINE OUTPUT_DBG_CHAR1

! ---------------------------------------------------------------

SUBROUTINE OUTPUT_DBG_IVAL(DBG_TITLE,VAL,VAL_TYPE,DBG_LEVEL)
  IMPLICIT NONE
#include "inc_print.h"
  CHARACTER(LEN=*) DBG_TITLE,VAL_TYPE
  INTEGER  VAL,DBG_LEVEL
  IF (RECP_DBG.GE.DBG_LEVEL) THEN
     IF (VAL_TYPE.EQ.'I3') THEN
        WRITE(RECP_OUT,'(A,A,A,I3)') &
           ' * ',DBG_TITLE,' : ',VAL
     ELSE
        WRITE(RECP_OUT,'(A,A)') ' * Wrong output type of variable ',DBG_TITLE
     ENDIF
  ENDIF 
END SUBROUTINE OUTPUT_DBG_IVAL

! ---------------------------------------------------------------

SUBROUTINE RECP_EXIT_IPT(RECP_LOC)
   IMPLICIT NONE
#include "inc_print.h"
   CHARACTER(LEN=*) :: RECP_LOC
   WRITE(RECP_OUT,*) ' - RECP : Error in allocating variable', RECP_LOC
   STOP
END SUBROUTINE RECP_EXIT_IPT

! ---------------------------------------------------------------
END MODULE RECP_OUTPUT
