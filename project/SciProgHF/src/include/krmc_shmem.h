!
! FILE    : krmc_shmem.h
!
!                MPI memory pointer 
!     ======================================
!
      integer, parameter :: my_MPI_ADDRESS_KIND = 8
      integer, parameter :: my_MPI_OFFSET_KIND  = 8
#if defined (VAR_MPI2) 
#if ! defined (SYS_AIX) && ! defined (VAR_XLF)
      INTEGER(KIND=my_MPI_ADDRESS_KIND) MY_T_PTR,MY_T_D_PTR,
     &                               MY_T_SCR_PTR,
     &                               MY_TTPL_PTR, MY_TCTR_PTR,
     &                               MY_TRED_PTR, MY_RCCTOS_PTR,
     &                               MY_IAPROC_PTR, MY_TTOL_PTR,
     &                               MY_CHECK_F_PTR, MY_XT_PTR,
     &                               MY_VEC1_PTR, MY_VEC2_PTR,
     &                               MY_BWEIGHT_PTR, MY_DL_PTR,
     &                               MY_T_PTRXXX 
#endif
!
      COMMON/MPI_MEMPOINTER1/  MY_T_SCR_PTR
      COMMON/MPI_MEMPOINTER2/  MY_T_PTR
      COMMON/MPI_MEMPOINTER3/  MY_T_D_PTR
      COMMON/MPI_MEMPOINTER4/  MY_TTPL_PTR
      COMMON/MPI_MEMPOINTER5/  MY_TTOL_PTR
      COMMON/MPI_MEMPOINTER6/  MY_CHECK_F_PTR
      COMMON/MPI_MEMPOINTER7/  MY_XT_PTR
      COMMON/MPI_MEMPOINTER8/  MY_VEC1_PTR
      COMMON/MPI_MEMPOINTER9/  MY_VEC2_PTR
      COMMON/MPI_MEMPOINTER10/ MY_BWEIGHT_PTR
      COMMON/MPI_MEMPOINTER11/ MY_BPROC_PTR
      COMMON/MPI_MEMPOINTER12/ MY_DL_PTR
      COMMON/MPI_MEMPOINTER13/ MY_T_PTRXXX
!
!             MPI memory window handler 
!     ======================================
!
      INTEGER MY_T_WIN
      INTEGER MY_C_WIN
!
      COMMON/MPI_WIN_HAND/ MY_T_WIN, MY_C_WIN
!
!             array length
!     ======================================
!
      INTEGER(KIND=my_MPI_OFFSET_KIND) LEN_T_BUFF, LEN_T_D_BUFF
      INTEGER(KIND=my_MPI_OFFSET_KIND) MY_TTPL_LEN, MY_TTOL_LEN
      INTEGER(KIND=my_MPI_OFFSET_KIND) MY_CHECK_F_LEN
      INTEGER(KIND=my_MPI_OFFSET_KIND) MY_XT_BUFF_LEN, MY_DL_TMP_LEN
      INTEGER(KIND=my_MPI_OFFSET_KIND) MY_VEC1_BUFF_LEN,MY_VEC2_BUFF_LEN
      INTEGER(KIND=my_MPI_OFFSET_KIND) MY_BWEIGHT_LEN, MY_BPROC_LEN
!
      COMMON/IBUFF_LENGTH/ LEN_T_BUFF, LEN_T_D_BUFF, MY_TTPL_LEN,
     &                     MY_TTOL_LEN, MY_CHECK_F_LEN, MY_XT_BUFF_LEN,
     &                     MY_VEC1_BUFF_LEN, MY_VEC2_BUFF_LEN, 
     &                     MY_BWEIGHT_LEN, MY_BPROC_LEN, MY_DL_TMP_LEN
!
#endif
