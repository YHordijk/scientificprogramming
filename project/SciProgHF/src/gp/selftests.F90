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

!
! Collection of self-tests of various functionalities used in DIRAC
! called at the beginning of dirac.x run by main.F90
!
! Tests involve:
!
! - Integer typecasting in lapack
! - MPI compatibility
! - Fortran - C/C++ data interoperability
!
! Written by : Ulf Eckstrom, Stefan Knecht,  Radovan Bast, Miro Ilias
!

      subroutine selftest_lapack
!     Call Lapack with some small tests to see if things look ok
!     written by Ulf Eckstrom
      implicit none
      integer N, N2, N3, lwork, liwork, info, idum, M, idum2,idum3,idum4
#ifdef INT_STAR8
      integer*4 N_4, N2_4, N3_4, lwork_4, liwork_4, info_4, idum_4, M_4, idum2_4,idum3_4,idum4_4, mone_4, mone2_4
#endif
      double precision A,DDUM,DDUM2,EIG,VEC,workdum
      N = 11
      N2 = N
      N3 = N
!     Query dsyevr work workspace size
!     If we use integer*8 for info and dsyevr expects integer*4
!     only parts of info will changed, and we can detect this later.
      INFO = -1
      liwork = -1
      CALL DSYEVR('V','A','U',N,A,N2,DDUM,DDUM2,IDUM,IDUM2,0.0D0,IDUM3,EIG,VEC,N3,IDUM4,WORKDUM,-1,liwork,-1,INFO)
      if (info.ne.0.or.liwork.lt.0) then
#ifdef INT_STAR8
!     integer*8 call failed, try with integer*4
      N_4 = 11
      N2_4 = N_4
      N3_4 = N_4
      INFO_4 = -1
      liwork_4 = -1
      mone_4 = -1
      mone2_4 = -1
      CALL DSYEVR('V','A','U',N_4,A,N2_4,DDUM,DDUM2,IDUM_4,        &
          IDUM2_4,0.0D0, IDUM3_4,EIG,VEC,N3_4,IDUM4_4,WORKDUM,     &
          mone_4,liwork_4,mone2_4,INFO_4)
      if (info_4.eq.0.and.liwork_4.gt.0) then
         print *,'WARNING: Self test failed!'
         print *,'Call to LAPACK(DSYEVR) failed with integer*8, but'
         print *,'worked with integer*4 arguments. However, DIRAC is'
         print *,'configured to use integer*8 (64 bit integers).'
         print *,' '
         print *,'You can fix this mismatch by linking to a integer*8'
         print *,'LAPACK library, or recompile DIRAC using integer*4.'
      else
         print *,'ERRROR: Lapack integer*8 test failed with info =',INFO,'LIWORK=',liwork
         print *,'ERRROR: Lapack integer*4 test failed with info =',INFO_4,'LIWORK=',liwork_4
      endif
#else
         print *,'ERRROR: Lapack test failed with info =',INFO,'LIWORK=',liwork
#endif
         print *,'Cannot continue without a working LAPACK, quitting.'
         call quit('LAPACK integer*4/8 test failed')
      endif
      print *,'LAPACK integer*4/8 selftest passed'
      end
!---------------------------------------------------------------------------------------------
!
!     MPI compatibility test
!
!---------------------------------------------------------------------------------------------
#ifdef VAR_MPI
      subroutine selftest_mpi
!     Call MPI routines with some small tests to see if things look ok
      use interface_to_mpi
      use integer_model_in_mpi
      implicit none
#include "infpar.h"
      INTEGER :: MPIERR, IREQ, IMESS,i
 !Miro: df_MPI_STATUS_SIZE is defined somwhere in interface_to_mpi ?
      INTEGER :: ISTAT(df_MPI_STATUS_SIZE)
      integer :: big_int, big_int_recv, big_int_arr(3), big_int_recv_arr(3)
      integer :: my_MPI_master=0
      integer :: my_MPI_slave=1
      integer :: my_MPI_rank, num_MPI_proc

      big_int=huge(i) ! the biggest signed integer*4/8 number, 2^31-1 or 2^63-1, resp.
      big_int_arr=big_int
      big_int_recv=0
      big_int_recv_arr=0

!     check number of MPI processes
      CALL interface_mpi_comm_size(global_communicator,num_MPI_proc)
!     if num_MPI_proc=1 only MASTER else we have SLAVES

!     get rank of processor
      CALL interface_mpi_comm_rank(global_communicator,my_MPI_rank)

      if (my_MPI_rank.eq.my_MPI_master) then
!     send dummy message (from MASTER to MASTER...)
!
      CALL interface_MPI_ISEND(MYTID,1,0,25,global_communicator,IREQ)
      CALL interface_MPI_PROBE(df_MPI_ANY_SOURCE,df_MPI_ANY_TAG, global_communicator,ISTAT)

      IF( ISTAT(df_MPI_SOURCE) .ne. MPARID )THEN
         print *,' *** ERROR: MPI self test failed! *** '
         print *,' sender of message should have ID 0'
         print *,' but has ID ',ISTAT(df_MPI_SOURCE)
         print *,' '
         print *,' Make sure that you are using the right mpirun for'
         print *,' your particular MPI library version.'
         print *,' '
#ifdef INT_STAR8
         print *,' DIRAC is configured to use integer*8 ' &
              //'(64 bit integers).'
         print *,' You need to link a MPI library that is compiled'
         print *,' using integer*8, make sure that this is the case.'
#else
         print *,' DIRAC is configured to use integer*4 '   &
              //' (32 bit integers).'
         print *,' You need to link a MPI library that is compiled'
         print *,' using integer*4, make sure that this is the case.'
#endif
         print *,' '
         print *,' Cannot continue without a working MPI library.'
         print *,' Quitting.'
         CALL QUIT(' *** ERROR: MPI self test failed ***')
      END IF

!
!     ... it looks ok, receive dummy message
!
      CALL interface_MPI_IRECV(IMESS,1,0,df_MPI_ANY_TAG,global_communicator,IREQ)
      CALL interface_MPI_WAIT(IREQ,ISTAT)

!  .. Miro: check single "big" integer*4/8 variable MPI sending and receiving
      CALL interface_MPI_ISEND(big_int,1,0,25,global_communicator,IREQ)
      CALL interface_MPI_PROBE(df_MPI_ANY_SOURCE,df_MPI_ANY_TAG, global_communicator,ISTAT)
      CALL interface_MPI_RECV(big_int_recv,1,0,df_MPI_ANY_TAG,global_communicator)

      if (big_int.ne.big_int_recv) then
        print *,'selftest_mpi: big_int, big_int_recv=',big_int, big_int_recv
        call quit('selftest_mpi:  probles with your MPI-integer setting ! ')
      endif

      endif
!  .. Miro: check array integer*4/8 variable MPI sending and receiving
!     only if we have a master and at least one slave

      if (num_MPI_proc.gt.1) then
        if (my_MPI_rank.eq.my_MPI_master) then
          !only master sending
          CALL interface_MPI_SEND(big_int_arr,3,my_MPI_slave,26,global_communicator)
        endif
        if (my_MPI_rank.eq.my_MPI_slave) then
          !only slave receiving
          !CALL interface_MPI_PROBE(df_MPI_ANY_SOURCE,df_MPI_ANY_TAG, global_communicator,ISTAT)
          CALL interface_MPI_RECV(big_int_recv_arr,3,my_MPI_master,df_MPI_ANY_TAG,global_communicator)
          if (big_int_arr(1).ne.big_int_recv_arr(1).or.  &
              big_int_arr(2).ne.big_int_recv_arr(2).or.  &
              big_int_arr(3).ne.big_int_recv_arr(3)) then
            print *,'selftest_mpi:    big_int_arr  =',big_int_arr
            print *,'selftest_mpi: big_int_recv_arr=',big_int_recv_arr
            call quit('selftest_mpi:  problems with your MPI-integer array setting ! ')
          endif
        endif
      endif

      if (my_MPI_rank.eq.my_MPI_master) then
      ! Finally, flush out the "all OK" message
      write(*,'(1X,A,I1,A)') 'MPI selftest passed with MPI_INTEGER of the size ',get_integer_model_in_mpi(),' bytes'
      endif

      END
#endif

subroutine selftest_fortran_c_cxx_interoperability()
!------------------------------------------------------------------------------------
!
! Test and tutorial example of Fortran - C/C++ interoperability.
!
! The Fortran-to-C/C++ interfacing is based solely upon the ISO_C_BINDING scheme.
!
! For the ISO_C_BINDING description, see  for example:
!         https://stackoverflow.com/tags/fortran-iso-c-binding/info
!         https://www-n.oca.eu/pichon/inter_C.pdf
!
! Passing Fortran Integer*4/Integer*8 data to equivalent ISO_C_BINDING
! c_int/c_long integer data types must be carried out as the ISO_C_BINDING
! integer data types have constant length.
!
! Developers can add more data types (incuding derived) for passing checks,
! depending on how heavily they employ C/C++ routines
! called from inside DIRAC Fortran sources.
!
! Written by Miro Ilias, Nov. 2015 after consulting with Roberto, Radovan, Ulf
!
!---------------------------------------------------------------------------------
  use iso_c_binding
  integer(c_int)    :: int_c_int, int_c_int_saved, string_len
  integer(c_long)   :: int_c_long, int_c_long_saved
  integer(c_long), target   :: int_c_long_target
  integer(c_long), allocatable, target :: int_c_long_array(:),int_c_long_array_saved(:)
  real(c_double), target    :: real_c_double, real_c_double_saved
  real(c_double), allocatable, target :: real_c_double_array(:),real_c_double_array_saved(:)

  type(c_ptr)       :: int_c_long_ptr, real_c_double_ptr
  integer(c_int)    :: arr_size = 3

  integer     :: i,fortran_integer  ! affected with Fortran integer*8 flag !
  character(len=20), target :: string_c_char,string_c_char_saved
  logical :: verify_int_long_double, verify_arrays, verify_chars, verify_sizeofs

  integer(kind=4) :: big_i4_integer=+2147483647_4 !  the biggest signed integer*4 number, 2^31-1
  integer(kind=8) :: big_i8_integer=+9223372036854775807_8 ! the biggest signed integer*8 number, 2^63-1

  interface
     subroutine selftest_fortran_c_interoperability(arr_size, int_c_int, int_c_long, int_c_long_array_ptr, &
                            real_c_double, real_c_double_array_ptr,string_c_char_c_ptr, string_len, &
                                        int_c_long_ptr, real_c_double_ptr) bind(c)
       use iso_c_binding
       type(c_ptr),value :: string_c_char_c_ptr
       integer(c_int)    :: int_c_int, arr_size, string_len   ! always 4-bytes
       integer(c_long)   :: int_c_long  ! always 8-bytes
       real(c_double)    :: real_c_double ! always 8-bytes
       type(c_ptr),value :: real_c_double_array_ptr
       type(c_ptr),value :: int_c_long_array_ptr
       type(c_ptr),value :: int_c_long_ptr, real_c_double_ptr
     end subroutine
  end interface

! assign numbers you are sending to C/C++ routine and bringing them back
!
  real_c_double=-1.797693D+308
  real_c_double_saved=real_c_double

! It's programmer's responsibility to ensure proper Fortran-to-C/C++ integer passing
! depending on the Fortran integer size (integer*4 or integer*8) !
#if defined (INT_STAR8)
  fortran_integer=+big_i8_integer
! we have Fortran integer*8 of the same length as int_c_long
  int_c_long=fortran_integer
  verify_sizeofs = sizeof(fortran_integer) == sizeof(int_c_long)
  int_c_int = +big_i4_integer
#else
  fortran_integer=+big_i4_integer
! we have Fortran integer*4 of the same length as int_c_int
  int_c_int=fortran_integer
  verify_sizeofs = sizeof(fortran_integer) == sizeof(int_c_int)
  int_c_long=+big_i8_integer
#endif
  int_c_int_saved=int_c_int
  int_c_long_saved=int_c_long

  allocate(real_c_double_array(arr_size)); allocate(real_c_double_array_saved(arr_size))
  real_c_double_array(1)=-real_c_double
  real_c_double_array(2)=+real_c_double
  real_c_double_array(3)=-real_c_double

  allocate(int_c_long_array(arr_size)); allocate(int_c_long_array_saved(arr_size))
  int_c_long_array(1)=+big_i8_integer
  int_c_long_array(2)=-big_i8_integer
  int_c_long_array(3)=+big_i8_integer

  do i=1,arr_size
    real_c_double_array_saved(i)=real_c_double_array(i)
    int_c_long_array_saved(i)=int_c_long_array(i)
  enddo

  string_c_char="Hello! 12345678+-@^/"
  string_len=len(string_c_char)
  do i=1,string_len
     string_c_char_saved(i:i)=string_c_char(i:i)
  enddo

  int_c_long_target=int_c_long
  int_c_long_ptr = c_loc(int_c_long_target)
  real_c_double_ptr = c_loc(real_c_double)

  call selftest_fortran_c_interoperability(arr_size, int_c_int, int_c_long, c_loc(int_c_long_array), &
  real_c_double, c_loc(real_c_double_array), c_loc(string_c_char(1:1)), string_len, int_c_long_ptr,real_c_double_ptr)


  ! verify that data are intact after returning from C routine
  verify_int_long_double = (int_c_int==int_c_int_saved) .and. (int_c_long==int_c_long_saved) .and. &
                            ( real_c_double==real_c_double_saved ) .and. verify_sizeofs
  verify_arrays=.true.
  do i=1,arr_size
     if (.not.(real_c_double_array_saved(i)==real_c_double_array(i) .and. &
               int_c_long_array_saved(i)==int_c_long_array(i))) verify_arrays=.false.
  enddo
  verify_chars=.true.
  do i=1, len(string_c_char)
      if (.not.(string_c_char(i:i) == string_c_char_saved(i:i) )) verify_chars=.false.
  enddo

  if (.not.(verify_int_long_double .and. verify_arrays .and. verify_chars)) then
     print *,"int_c_int,int_c_int_saved:",int_c_int,int_c_int_saved
     print *,"int_c_long,int_c_long_saved:",int_c_long,int_c_long_saved
     print *,"real_c_double,real_c_double_saved:",real_c_double,real_c_double_saved
     print *,"WARNING: Selftest of ISO_C_BINDING Fortran - C/C++ interoperability FAILED !"
  else
     print *,"Selftest of ISO_C_BINDING Fortran - C/C++ interoperability PASSED"
  endif

end subroutine
