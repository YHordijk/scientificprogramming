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
#if defined (VAR_MPI2)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C======================================================================C
C                                                                      C
C    Macros for MPI-2 calls that make calls to MPI-2 routines          C
C    easier, more transparent, and independent from a MPI-1            C 
C    only library.                                                     C
C                                                                      C
C    This collection is inspired by gp/mpimacro.F (DALTON package)     C
C    collected by Paal Dahle Jan-2001                                  C
C                                                                      C
C    Stefan Knecht - October 2007                                      C
C                                                                      C
C======================================================================C
C
C     /* deck mpixget */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MPIXGET( RBUF, JCOUNT, NTYPE_R, ITARGET, IDISPL, 
     &                    JCOUNT_T, NTYPE_T, MY_WIN)
C***********************************************************************
C
C     Communicate a scalar/vector/matrix via a remote memory access (RMA) 
C     routine provided by MPI-2. Data are put from the target memory  
C     to the origin.
C
C     OUTPUT: 
C            RBUF:      updated origin buffer with JCOUNT elements of
C                       type NTYPE_R
C     INPUT: 
C            MY_WIN:    memory window on ITARGET (must be initialized!)
C                       accessed at displacement IDISPL with JCOUNT_T 
C                       elements of type NTYPE_T
C
C     Written by S. Knecht - Oct 17 2007
C     Last revision : 
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "priunit.h"
#include "mpif.h"
#include "mxdim_mpi2.h"
      INTEGER ISTAT(MPI_STATUS_SIZE)
      DIMENSION RBUF(*), MY_WIN(*)
      INTEGER NTYPE_R, NTYPE_T, DATATYPE_OUT, DATATYPE_IN
      INTEGER JCOUNT, JCOUNT_T
      INTEGER(KIND=MPI_ADDRESS_KIND) IDISPL
C
      IF( NTYPE_R .eq. 1 ) THEN
        DATATYPE_OUT = my_MPI_INTEGER
      ELSE IF ( NTYPE_R .eq. 2 ) THEN
        DATATYPE_OUT = MPI_REAL8
      ELSE
        WRITE(LUPRI,*)' MPIXGET: output datatype "',NTYPE_R,
     &                  '" does not exist!' 
        CALL Abend2('ERROR in MPIXGET:: Nonexisting NTYPE_R')
      END IF
C
      IF( NTYPE_T .eq. 1 ) THEN
        DATATYPE_IN = my_MPI_INTEGER
      ELSE IF ( NTYPE_T .eq. 2 ) THEN
        DATATYPE_IN = MPI_REAL8
      ELSE
        WRITE(LUPRI,*)' MPIXGET: output datatype "',NTYPE_T,
     &                  '" does not exist!' 
        CALL Abend2('ERROR in MPIXGET:: Nonexisting NTYPE_T')
      END IF
C
C     lock window (MPI_LOCK_SHARED mode)
C
      CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,ITARGET,MPI_MODE_NOCHECK,
     &                  MY_WIN,IERR)
C
C     transfer data     
C
      CALL MPI_GET(RBUF,JCOUNT,DATATYPE_OUT,ITARGET,IDISPL,JCOUNT_T,
     &              DATATYPE_IN, MY_WIN, IERR)
C
C     unlock
C
      CALL MPI_WIN_UNLOCK(ITARGET,MY_WIN,IERR)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     /* deck mpixaccum */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MPIXACCUM( RBUF, JCOUNT, NTYPE_R, ITARGET, IDISPL, 
     &                     JCOUNT_T, NTYPE_T, MY_WIN)
C***********************************************************************
C
C     accumulate a scalar/vector/matrix via a remote memory access (RMA) 
C     routine provided by MPI-2. Data are put from the origin memory  
C     to the target memory.
C
C     OUTPUT: 
C            MY_WIN::   updated target memory window with JCOUNT_T 
C                       elements of type NTYPE_T
C     INPUT: 
C            RBUF   :   memory buffer on origin with JCOUNT elements of 
C                       type NTYPE_R.
C
C     allowed OPERATION: MPI_SUM
C
C     Written by S. Knecht - Nov 28 2007
C     Last revision : 
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
#include "priunit.h"
#include "mxdim_mpi2.h"
      INTEGER ISTAT(MPI_STATUS_SIZE)
      DIMENSION RBUF(*), MY_WIN(*)
      INTEGER NTYPE_R, NTYPE_T, DATATYPE_OUT, DATATYPE_IN
      INTEGER JCOUNT, JCOUNT_T
      INTEGER(KIND=MPI_ADDRESS_KIND) IDISPL
C
      IF( NTYPE_T .eq. 1 ) THEN
        DATATYPE_OUT = my_MPI_INTEGER
      ELSE IF ( NTYPE_T .eq. 2 ) THEN
        DATATYPE_OUT = MPI_DOUBLE_PRECISION
      ELSE
        WRITE(LUPRI,*)' MPIXACCUM: output datatype "',NTYPE_T,
     &                  '" does not exist!' 
        CALL Abend2('ERROR in MPIXACCUM: Nonexisting NTYPE_T')
      END IF
C
      IF( NTYPE_R .eq. 1 ) THEN
        DATATYPE_IN = my_MPI_INTEGER
      ELSE IF ( NTYPE_R .eq. 2 ) THEN
        DATATYPE_IN = MPI_DOUBLE_PRECISION
      ELSE
        WRITE(LUPRI,*)' MPIXACCUM: output datatype "',NTYPE_R,
     &                  '" does not exist!' 
        CALL Abend2('ERROR in MPIXACCUM: Nonexisting NTYPE_R')
      END IF
C
C     lock window (MPI_LOCK_SHARED mode)
C
      CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,ITARGET,MPI_MODE_NOCHECK,
     &                  MY_WIN,IERR)
C
      CALL MPI_ACCUMULATE( RBUF, JCOUNT, DATATYPE_IN, ITARGET, IDISPL,
     &                    JCOUNT_T, DATATYPE_OUT, MPI_SUM, MY_WIN, IERR)
C
C     unlock
C
      CALL MPI_WIN_UNLOCK(ITARGET,MY_WIN,IERR)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     /* deck mpixmem_alloc */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MPIXMEM_ALLOC( MY_POINT, NTYPE, NELEMENT, 
     &                          MY_INFO, PROV_INFO )
C***********************************************************************
C
C     allocate memory by MPI_ALLOC_MEM that is in particular suitable
C     for memory window handling (lock + unlock). 
C
C     INPUT: 
C            memory pointer MY_POINT (KIND=MPI_ADDRESS_KIND)
C            number of elements NELEMENT of type NTYPE
C     SAFETY NOTE: NELEMENT must be an INTEGER*8 (KIND=MPI_OFFSET_KIND)
C            
C            NTYPE == 1: REAL*8
C            NTYPE == 2: my_MPI_INTEGER
C
C     OUTPUT: 
C            allocated memory
C                       
C     Written by S. Knecht - Oct 19 2007
C     Last revision : 
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
#include "priunit.h"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
      INTEGER(kind=MPI_ADDRESS_KIND) MY_MM_NEEDS, MY_MM_NEEDSSCR
      INTEGER(kind=MPI_ADDRESS_KIND) MY_POINT
      INTEGER(kind=MPI_OFFSET_KIND) NELEMENT
      LOGICAL PROV_INFO
C
C     how much memory?
C
      MY_MM_NEEDS    = 0 
      MY_MM_NEEDSSCR = 0
C
C     do we provide any information (in particular for 
C     MPI implementations of NEC, IBM, ...)
C
      IF( PROV_INFO )THEN
        IS_INFO = MY_INFO
      ELSE
        IS_INFO = MPI_INFO_NULL
      END IF
C
C     allocate memory type specific
C
      IF( NTYPE .eq. 1 ) THEN
        MY_MM_NEEDSSCR = ISIZE_DP
        MY_MM_NEEDS    = MY_MM_NEEDSSCR * NELEMENT
      ELSE IF ( NTYPE .eq. 2 )THEN
        MY_MM_NEEDSSCR = ISIZE_my_int
        MY_MM_NEEDS    = MY_MM_NEEDSSCR * NELEMENT
      ELSE IF ( NTYPE .eq. 3 )THEN
        MY_MM_NEEDSSCR = ISIZE_int8
        MY_MM_NEEDS    = MY_MM_NEEDSSCR * NELEMENT
      ELSE IF ( NTYPE .eq. 4 )THEN
        MY_MM_NEEDSSCR = ISIZE_log
        MY_MM_NEEDS    = MY_MM_NEEDSSCR * NELEMENT
      END IF
C
csk     WRITE(LUPRI,*) ' allocation of MY_MM_NEEDS',MY_MM_NEEDS
C
      CALL MPI_ALLOC_MEM(MY_MM_NEEDS,IS_INFO,MY_POINT,IERR)
C
C     WRITE(LUPRI,*) ' allocation successful'
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     /* deck mpixmem_free */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MPIXMEM_FREE(RBUF)
C***********************************************************************
C
C     free memory buffer RBUF
C
C     INPUT: 
C            allocated memory buffer RBUF
C            
C     OUTPUT: 
C                       
C     Written by S. Knecht - Oct 19 2007
C     Last revision : 
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
C
      CALL MPI_FREE_MEM(RBUF,IERR)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     /* deck mpixcalc_typesz */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MPIXCALC_TYPESZ()
C***********************************************************************
C
C     calculate various MPI type sizes (in byte) - store on common block
C     TYPE_SZMPI2 (typesz_mpi2.h)
C
C     Written by S. Knecht - Oct 19 2007
C     Last revision : 
C
C
C***********************************************************************
      use interface_to_mpi
      IMPLICIT REAL*8 (A-H,O-Z)
#include "priunit.h"
#include "infpar.h"
#include "typesz_mpi2.h"
C
      call interface_MPI_TYPE_EXTENT(df_MPI_REAL8,ISIZE_DP)
      call interface_MPI_TYPE_EXTENT(df_my_MPI_INTEGER,ISIZE_my_int)
      call interface_MPI_TYPE_EXTENT(df_MPI_INTEGER8,ISIZE_int8)
      call interface_MPI_TYPE_EXTENT(df_MPI_INTEGER,ISIZE_intdef)
      call interface_MPI_TYPE_EXTENT(df_MPI_LOGICAL,ISIZE_log)

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     /* deck mpixwindow_crt */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MPIXWINDOW_CRT(RBUF,NELEMENT,NTYPE,IWIN_COMM,MY_WIN)
C***********************************************************************
C
C     open memory window to be used in one-sided MPI communication
C
C     INPUT:
C            array RBUF (should be allocated by mpixmem_alloc to
C            assure lock functionality)
C            number of elements NELEMENT of type NTYPE
C
C            NTYPE == 1: REAL*8
C            NTYPE == 2: my_MPI_INTEGER
C
C
C     OUTPUT:
C            new memory window handle MY_WIN shared by all
C            processes in calling communicator IWIN_COMM (group).
C            extension of memory window on each process may be
C            unsymmetric.
C
C     Written by S. Knecht - Oct 19 2007
C     Last revision :
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
#include "priunit.h"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
      INTEGER(kind=MPI_ADDRESS_KIND) IBUF_LEN, IBUF_LENSCR
      INTEGER(kind=MPI_ADDRESS_KIND) NELEMENT
C
      IBUF_LENSCR = 0
      IBUF_LEN    = 0
C
C     open memory window on each process shared by IWIN_COMM
C     ISIZE_dp resp. ISIZE_int is scaling unit in memory window
C     --> scale by REAL*8 resp. size_of_int
C
      IF( NTYPE .eq. 1 ) THEN
C
        IBUF_LENSCR = ISIZE_dp * NELEMENT
        IBUF_LEN    = IBUF_LENSCR
C
csk        WRITE(LUPRI,*) ' MPIXWINDOW_CRT: IBUF_LEN and ISIZE_dp',
csk     &                                   IBUF_LEN, ISIZE_dp
        CALL MPI_WIN_CREATE(RBUF,IBUF_LEN,ISIZE_dp,MPI_INFO_NULL,
     &                      IWIN_COMM,MY_WIN,IERR)
C
      ELSE IF ( NTYPE .eq. 2 )THEN
C
        IBUF_LENSCR = ISIZE_int * NELEMENT
        IBUF_LEN    = IBUF_LENSCR
csk     WRITE(LUPRI,*) ' MPIXWINDOW_CRT: IBUF_LEN and ISIZE_int',
csk  &                                   IBUF_LEN, ISIZE_int
C
        CALL MPI_WIN_CREATE(RBUF,IBUF_LEN,ISIZE_int,MPI_INFO_NULL,
     &                      IWIN_COMM,MY_WIN,IERR)
C
      END IF
C
csk      WRITE(LUPRI,*) ' window opened: MY_WIN', MYID
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     /* deck mpixwindow_free */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MPIXWINDOW_FREE(MY_WIN)
C***********************************************************************
C
C     close memory window used in one-sided MPI communication
C
C     INPUT: 
C            memory window handle MY_WIN shared by all processes in
C            calling communicator IWIN_COMM (group).
C
C     OUTPUT: 
C                       
C     Written by S. Knecht - Oct 19 2007
C     Last revision : 
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
#include "priunit.h"
C
C     close memory window on each process shared by IWIN_COMM
C
      CALL MPI_WIN_FREE(MY_WIN,IERR)
C
C     WRITE(LUPRI,*) ' window closed', MYID
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     /* deck mpixwingroup_info */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MPIXWINGROUP_INFO(MY_WIN,INFO_GROUP)
C***********************************************************************
C
C     get info on group of processes sharing a memory window MY_WIN
C
C     INPUT: 
C            shared memory window handle MY_WIN 
C
C     OUTPUT: 
C            duplicate (INFO_GROUP) of initial group handle of processes
C            sharing the memory window MY_WIN
C                       
C     Written by S. Knecht - Oct 22 2007
C     Last revision : 
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
#include "priunit.h"
C
C     get group information
C
      CALL MPI_WIN_GET_GROUP(MY_WIN, INFO_GROUP, IERR)
C
C     WRITE(LUPRI,*) ' info gathered', MYID
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     /* deck mpixwin_attr_info */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MPIXWIN_ATTR_INFO(MY_WIN,INUM)
C***********************************************************************
C
C     get cached attributes on shared memory window
C
C     INPUT: 
C            shared memory window handle MY_WIN 
C            attribute number: INUM     == 1: window base address 
C                                       == 2: window size (in bytes)
C                                       == 3: displacement unit (in
C                                             memory window)
C
C     OUTPUT: 
C                       
C     Written by S. Knecht - Oct 22 2007
C     Last revision : 
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
#include "priunit.h"
      INTEGER IKEY_VAL
      INTEGER(KIND=MPI_ADDRESS_KIND) IWIN_ATTR
      LOGICAL VALUE_PROV
C
      VALUE_PROV = .FALSE.
      IKEY_VAL  = 0
      IWIN_ATTR = 0
C
C     which information is needed?
C
      IF( INUM .eq. 1 ) THEN
        IKEY_VAL = MPI_WIN_BASE
      ELSE IF( INUM .eq. 2 ) THEN
        IKEY_VAL = MPI_WIN_SIZE
      ELSE IF ( INUM .eq. 3 ) THEN
        IKEY_VAL = MPI_WIN_DISP_UNIT
      ELSE
         WRITE(LUPRI,*) 'invalid memory window attribute requested',
     &                   INUM
         GOTO 100
      END IF
C
      CALL MPI_WIN_GET_ATTR(MY_WIN,IKEY_VAL,IWIN_ATTR,VALUE_PROV,IERR)
C
C     print information
C
      IF( INUM .eq. 1 .and. VALUE_PROV ) THEN 
        WRITE(LUPRI,*) ' window base address at', IWIN_ATTR 
      ELSE IF( INUM .eq. 2 .and. VALUE_PROV ) THEN 
        WRITE(LUPRI,*) ' size of window is (in bytes)', IWIN_ATTR 
      ELSE IF( INUM .eq. 3 .and. VALUE_PROV ) THEN 
        WRITE(LUPRI,*) ' displacement unit for window is', IWIN_ATTR 
      END IF 
C
 100  CONTINUE
C
      END

#else
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DUMMY_MPIXXX2
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#endif /* (VAR_MPI2)*/
