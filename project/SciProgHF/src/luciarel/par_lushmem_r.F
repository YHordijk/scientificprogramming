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
***********************************************************************

      SUBROUTINE CALC_WIN_OFFSET(IT_TTPL,IT_TTOL,NDIM,IMULT)
C
C     Written by  S. Knecht         - November 28 2007
C
C**********************************************************************
C
C     calculate CPU-specific window offsets for totat 
C     T block window MY_T_WIN
C     
C
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "krmc_shmem.h"
#include "parluci.h"
      DIMENSION IT_TTPL(*) 
      INTEGER(KIND=MPI_OFFSET_KIND) IT_TTOL(*), nelm_cc
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_PROC, ISLEN_ACT_BLK
      IOFF_PROC     = 0
      ISLEN_ACT_BLK = 0
C
C     calculate window offset for each CPU forming communicator
C     MYNEW_COMM_SM
C
csk      WRITE(LUWRT,*) ' CALC_WIN_OFFSET start:'
csk      WRITE(LUWRT,*) ' NEWCOMM_PROC_SM      :',NEWCOMM_PROC_SM
csk      WRITE(LUWRT,*) ' IMULT                :',IMULT          
csk      WRITE(LUWRT,*) ' NINT_TOTAL           :',NINT_TOTAL     
csk      WRITE(LUWRT,*) ' NDIM                 :',NDIM            
C
      DO IPROC = 1, NEWCOMM_PROC_SM
        IOFF_PROC = 0
        DO IRUN = 1, IMULT
          DO IBLK = 1, NDIM
csk            WRITE(LUWRT,*) 'request: IBLK + (NDIM*(IRUN -1))',
csk     &      IBLK + (NDIM*(IRUN - 1)), IT_TTPL( IBLK + (NDIM*(IRUN - 1)))
            IF( IT_TTPL( IBLK + (NDIM*(IRUN - 1)) ) .eq. IPROC - 1) THEN
              ISLEN_ACT_BLK = NELM_CC(IBLK,NDIM,NINT_TOTAL)
              IT_TTOL( IBLK + (NDIM*(IRUN -1)) ) = IOFF_PROC
              IOFF_PROC     = IOFF_PROC + ISLEN_ACT_BLK
csk              WRITE(LUWRT,*) ' IBLK, IOFF_PROC, ISLEN_ACT_BLK:',
csk     &                         IBLK, IOFF_PROC, ISLEN_ACT_BLK  
            END IF 
          END DO
        END DO
      END DO
C
      END
***********************************************************************

      SUBROUTINE GET_MEM_T( MY_T_INFO )
C
C     Written by  S. Knecht         - November 17 2007
C
C**********************************************************************
C
C     routine for memory allocation of the integral/density block T 
C     to use in parallel shared-memory MPI LUCIAREL runs
C
C     pointer has to be defined in calling routine as 
C     POINTER (MY_T_PTR, T_BUFF(IDUMMY_BUFF_SZ))
C
C     MY_T_PTR   stored on /MPI_MEMPOINTER/ (krmc_shmem.h)
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "krmc_shmem.h"
#include "parluci.h"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
      LOGICAL CINFO
      CINFO = .FALSE.
C
C     allocate T memory
C
      IF( MY_T_INFO .ne. 0 ) CINFO = .TRUE.
C
      CALL MPIXMEM_ALLOC( MY_T_PTR, 1, MY_T_LEN , MY_T_INFO, CINFO )
C
      END
***********************************************************************

      SUBROUTINE GET_MEM_T_D( MY_T_INFO )
C
C     Written by  S. Knecht         - February 11 2008
C
C**********************************************************************
C
C     memory allocation of the pure diagonal integral block T
C     to use in parallel shared-memory MPI LUCIAREL runs
C
C     pointer has to be defined in calling routine as 
C     POINTER (MY_T_D_PTR, T_BUFF_D(IDUMMY_BUFF_SZ))
C
C     LEN_T_D_BUFF stored on /IBUFF_LENGTH/   (krmc_shmem.h)
C     MY_T_D_PTR   stored on /MPI_MEMPOINTER/ (krmc_shmem.h)
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
#include "krmc_shmem.h"
      LOGICAL CINFO
      CINFO = .FALSE.
C
C     allocate T_DIAG memory
C
      IF( MY_T_INFO .ne. 0 ) CINFO = .TRUE.
C
      CALL MPIXMEM_ALLOC(MY_T_D_PTR, 1, LEN_T_D_BUFF, MY_T_INFO, CINFO)
C
      END
***********************************************************************

      SUBROUTINE GET_MEM_T_SCR_CC( MY_T_INFO )
C
C     Written by  S. Knecht         - November 28 2007
C
C**********************************************************************
C
C     routine for memory allocation of the scratch integral/density 
C     block T 
C     to use in parallel shared-memory MPI LUCIAREL runs
C
C     pointer has to be defined in calling routine as 
C     POINTER (MY_T_SCR_PTR, T_SCR_CC(IDUMMY_BUFF_SZ))
C
C     IS_LENGTH_TT   stored on /LUCIPAR_I8/   (parluci.h)
C     MY_T_SCR_PTR   stored on /MPI_MEMPOINTER/ (krmc_shmem.h)
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
#include "krmc_shmem.h"
      LOGICAL CINFO
      CINFO = .FALSE.
C
C     allocate T memory
C
      IF( MY_T_INFO .ne. 0 ) CINFO = .TRUE.
C
      CALL MPIXMEM_ALLOC( MY_T_SCR_PTR,1,IS_LENGTH_TT,MY_T_INFO,CINFO)
C
      END
***********************************************************************

      SUBROUTINE GET_MEM_XINT_SCR( MY_XINT_INFO )
C
C     Written by  S. Knecht         - December 03 2007
C
C**********************************************************************
C
C     routine for memory allocation of the scratch integral 
C     block XINT_SCR 
C     to use in parallel shared-memory MPI LUCIAREL runs
C
C     pointer has to be defined in calling routine as 
C     POINTER (MY_XINT_SCR_PTR, XINT_SCR(IDUMMY_BUFF_SZ))
C
C     LEN_ALL_INT       stored on /LUCIPAR_I8/   (parluci.h)
C     MY_XINT_SCR_PTR   stored on /MPI_MEMPOINTER/ (krmc_shmem.h)
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
#include "krmc_shmem.h"
      LOGICAL CINFO
      CINFO = .FALSE.
C
C     allocate XINT_SCR memory
C
      IF( MY_XINT_INFO .ne. 0 ) CINFO = .TRUE.
C
      CALL MPIXMEM_ALLOC( MY_XINT_SCR_PTR, 1, LEN_ALL_INT, MY_XINT_INFO,
     &                    CINFO)
C
      END
***********************************************************************

      SUBROUTINE GET_T_WIN( T_BUFF, MY_T_COMM, MY_LEN_T )
C
C     Written by  S. Knecht         - November 17 2007
C
C**********************************************************************
C
C     routine for window creation of the integral/density block T to 
C     use in parallel shared-memory MPI LUCIAREL runs
C
C     T_BUFF is allocated using MPI_ALLOC_MEM.
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
#include "krmc_shmem.h"
      INTEGER(KIND=MPI_ADDRESS_KIND) MY_LEN_T
C
C     open T window
C
csk      write(luwrt,*) ' GET_T_WIN: MY_LEN_T, ISIZE_dp, MY_T_COMM',
csk     &                            MY_LEN_T, ISIZE_dp, MY_T_COMM
      CALL MPIXWINDOW_CRT( T_BUFF, MY_LEN_T, 1, MY_T_COMM, MY_T_WIN)
C
      END
***********************************************************************

      SUBROUTINE GET_CVEC_WIN( C_VEC, MY_C_COMM, MY_LEN_C_VEC )
C
C     Written by  S. Knecht         - February 25 2008
C                                     back in Odense :)
C
C**********************************************************************
C
C     routine for window creation of the vector block C_VEC to 
C     use in parallel shared-memory MPI LUCIAREL runs
C
C     C_VEC is allocated using MPI_ALLOC_MEM.
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
#include "krmc_shmem.h"
      INTEGER(KIND=MPI_ADDRESS_KIND) MY_LEN_C_VEC
C
C     open C_VEC window
C
csk   write(luwrt,*) ' GET_CVEC_WIN: MY_LEN_C_VEC,ISIZE_dp,MY_C_COMM',
csk  &                               MY_LEN_C_VEC,ISIZE_dp,MY_C_COMM
      CALL MPIXWINDOW_CRT( C_VEC, MY_LEN_C_VEC, 1, MY_C_COMM, MY_C_WIN)
C
      END
***********************************************************************

      SUBROUTINE REL_MEM_T_WIN( T_BUFF )
C
C     Written by  S. Knecht         - November 19 2007
C
C**********************************************************************
C
C     deallocate and and destroy memory window for integral/density 
C     block T. 
C
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "krmc_shmem.h"
#include "parluci.h"
C
C     destroy T memory window
C
      CALL MPIXWINDOW_FREE( MY_T_WIN )
C
C     release memory
C
      CALL MPIXMEM_FREE(T_BUFF)
C
      END
***********************************************************************

      SUBROUTINE REL_MEM_C_WIN( C_VEC , C_VEC_2)
C
C     Written by  S. Knecht         - February 25 2008
C                                     back in Odense :)
C
C**********************************************************************
C
C     deallocate and and destroy memory window for vector block 
C     C_VEC resp. release memory for C_VEC_2.
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "krmc_shmem.h"
#include "parluci.h"
C
C     destroy C memory window
C
      CALL MPIXWINDOW_FREE( MY_C_WIN )
C
C     release memory
C
      CALL MPIXMEM_FREE( C_VEC )
      CALL MPIXMEM_FREE( C_VEC_2 )
C
      END
***********************************************************************

      SUBROUTINE REL_MEM_T_SCR( T_SCR_CC )
C
C     Written by  S. Knecht         - November 28 2007
C
C**********************************************************************
C
C     deallocate integral/density scratch block T_SCR_CC. 
C
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "krmc_shmem.h"
#include "parluci.h"
C
C     release memory
C
      CALL MPIXMEM_FREE( T_SCR_CC )
C
      END
***********************************************************************

      SUBROUTINE REL_MEM_XINT_SCR( XINT_SCR )
C
C     Written by  S. Knecht         - November 28 2007
C
C**********************************************************************
C
C     deallocate integral scratch block XINT_SCR. 
C
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#include "krmc_shmem.h"
#include "parluci.h"
C
C     release memory
C
      CALL MPIXMEM_FREE( XINT_SCR )
C
      END
#else
***********************************************************************

      SUBROUTINE PAR_LUSHMEM_REL
*     dummy routine for normal compilation
      END
#endif

