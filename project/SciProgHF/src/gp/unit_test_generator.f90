!================================================================================!
!                                                                                !
!  Program:      ReSpect, Dirac                                                  !
!                                                                                !
!  Module:       UNIT_TEST_GENERATOR                                             !
!                                                                                !
!  Description:  automatically generates unit tests                              !
!                                                                                !
!  Dependencies: string_manipulations                                            !
!                                                                                !
!  Author:       Stanislav Komorovsky                                            !
!                                                                                !
!  Date:         Jan. - July 2012 (Toulouse, Tromso)                             !
!                                                                                !
!================================================================================!
!                                                                                !
!  EXAMPLE how to use unit_test_generator                                        !
!      (for more examples see unit tests of unit_test_generator)                 !
!      (      source_path/RS/gp/gp_test/unit_test_generator    )                 !
!                                                                                !
!  module some_module                                                            !
!    use unit_test_generator                                                     !
!    ...                                                                         !
!  contains                                                                      !
!    subroutine some_routine(inp1,inp2,out3,...)                                 !
!    implicit none                                                               !
!    ...                                                                         !
!  #if defined CREATE_UNIT_TESTS                                                 !
!    call unit_test__constructor()                                               !
!    call unit_test__set_suffix(var1)                                            !
!    call unit_test__set_routine_name('some_routine')                            !
!    call unit_test__add_module('some_module')                                   !
!    call unit_test__add_inp_variable(inp1  ,'inp1  ')                           !
!    call unit_test__add_inp_variable(inp2  ,'inp2  ')                           !
!    ...                                                                         !
!  #endif                                                                        !
!    ...                                                                         !
!    (body of the routine)                                                       !
!    ...                                                                         !
!  #if defined CREATE_UNIT_TESTS                                                 !
!    call unit_test__add_out_variable(out3  ,'out3  ')                           !
!    call unit_test__add_out_variable(out4  ,'out4  ')                           !
!    call unit_test__generate_test()                                             !
!  #endif                                                                        !
!    end subroutine                                                              !
!  end module                                                                    !
!                                                                                !
!--------------------------------------------------------------------------------!
!                                                                                !
!  Create test program                                                           !
!  $ setup --unit-tests=on ...                                                   !
!  $ cd build_dir                                                                !
!  $ make                                                                        !
!  Run some reference job                                                        !
!  In scratch dir you will find testing program with name test.F90_inp1          !
!  Copy it to directory source_dir/...../generalname_test/name-of-the-test       !
!    -->  "_test" is mandatory                                                   !
!  Then setup new compilation with --tests=on option, compile and run unit tests !
!  $ setup --tests=on ...                                                        !
!  $ cd build_dir                                                                !
!  $ make                                                                        !
!  $ make test                                                                   !
!                                                                                !
!================================================================================!

MODULE UNIT_TEST_GENERATOR

use string_manipulations, i2c  => string__i2c  &
                        , l2c  => string__l2c
implicit none

          ! Array used in this calls must have assumed size (*):
public :: unit_test__add_inp_int_array     ! routine for adding integer input
public :: unit_test__add_inp_real_array    ! routine for adding real input
public :: unit_test__add_out_int_array     ! routine for adding integer output
public :: unit_test__add_out_real_array    ! routine for adding real output

          ! Array used in this calls must have explicit (n) or assumed shape (:):
public :: unit_test__add_inp_variable      ! general routine for adding input variables
public :: unit_test__add_out_variable      ! general routine for adding output variables
public :: unit_test__add_opt_variable      ! general routine for adding optional variables

public :: unit_test__set_routine_name      ! set name of driver routine
public :: unit_test__add_input_file        ! add name of input file
public :: unit_test__add_output_file       ! add name of output file
public :: unit_test__add_module            ! add module name
public :: unit_test__set_keys              ! set keywords
public :: unit_test__redefine_variable     ! redefine variable parameters
public :: unit_test__set_suffix            ! set suffix for program name
public :: unit_test__get_output_name       ! return name of the output variable
public :: unit_test__generate_test         ! generate unit test
public :: unit_test__constructor           ! initialize unit test generator
public :: unit_test__destructor            ! deinitialize unit test generator

private

interface unit_test__add_inp_variable
  module procedure add_inp_string_variable
  module procedure add_inp_logical_variable
  module procedure add_inp_integer_variable
  module procedure add_inp_integer_1d_array
  module procedure add_inp_integer_2d_array
  module procedure add_inp_real_variable
  module procedure add_inp_real_1d_array
  module procedure add_inp_real_2d_array
  module procedure add_inp_real_3d_array
  module procedure add_inp_real_4d_array
  module procedure add_inp_real_5d_array
  module procedure add_inp_complex_variable
  module procedure add_inp_complex_1d_array
  module procedure add_inp_complex_2d_array
  module procedure add_inp_complex_3d_array
  module procedure add_inp_complex_4d_array
end interface unit_test__add_inp_variable

interface unit_test__add_out_variable
  ! compare     --> not present ==> compare two arrays
  ! check       --> not present ==> 'total'   string will be passed to unit_testing module
  ! threshlods  --> not present ==> '1.0d-15' string will be passed to unit_testing module
  module procedure add_out_integer_variable
  module procedure add_out_integer_1d_array
  module procedure add_out_integer_2d_array
  module procedure add_out_real_variable
  module procedure add_out_real_1d_array
  module procedure add_out_real_2d_array
  module procedure add_out_real_3d_array
  module procedure add_out_real_4d_array
  module procedure add_out_real_5d_array
  module procedure add_out_complex_variable
  module procedure add_out_complex_1d_array
  module procedure add_out_complex_2d_array
  module procedure add_out_complex_3d_array
  module procedure add_out_complex_4d_array
end interface unit_test__add_out_variable

interface unit_test__set_suffix
  module procedure set_suffix_string
  module procedure set_suffix_integer
  module procedure set_suffix_logical
end interface unit_test__set_suffix

private :: errall, erropen
private :: add_data_routine, find_name, is_equal, add_compare_info
private :: write_data, write_integer_array, write_real_array, write_complex_array
private :: associate_status, gets, append, change, get_size, get_char, init, free

interface is_equal
  module procedure is_equal_type_type
  module procedure is_equal_string_type
end interface is_equal

interface write_data
  module procedure write_string_variable
  module procedure write_logical_variable
  module procedure write_integer_variable
  module procedure write_real_variable
  module procedure write_complex_variable
end interface write_data

integer, parameter, private :: iunit = 60  ! unit used for storing some data and test.F90
integer, parameter, private :: max_nr_of_records = 100, max_length_of_string = 70

type, private :: list
  type( string ) :: records(max_nr_of_records)
  integer :: nr_records
end type list

type( list ), private, save :: input_files, output_files, driver, test_name
type( list ), private, save :: var_name, var_type, var_dim
type( list ), private, save :: var_info, var_value, var_routine
type( list ), private, save :: modules, suffix
type( list ), private, save :: compare_technique, compare_info, thresholds, threshold_type
type( list ), private, save :: compare_name, compare_type, compare_dim

type, private :: keys
  logical :: mod_call = .true.       ! driver routine is part of f90 module
  logical :: filter   = .true.       ! use filter for storing of output arrays
end type keys

type( keys ), private, save :: key

logical, private, save :: iam_initialized = .False.  ! remember if the module is initialized
character(len=max_length_of_string), save :: current_routine = ''

CONTAINS

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__add_inp_int_array(ivar,vname,ndim,work)

  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ndim, ivar(*)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: ier
  integer, allocatable :: itmp(:)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      allocate( itmp(ndim), stat=ier ); call errall(ier,'unit_test__add_inp_int_array','itmp')
      itmp(1:ndim) = ivar(1:ndim)

      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'unit_test__add_inp_int_array',trim(gets(input_files,1)))

      write(iunit,'(a)')    vname
      write(iunit,'(6i12)') itmp

      close(iunit,status='keep')

      deallocate(itmp)
    else
      call write_integer_array(ivar,ndim,vname)
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('integer, allocatable, dimension(:)',var_type)
  call append('(' // trim(i2c(ndim)) // ')',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine unit_test__add_inp_int_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__add_inp_real_array(dvar,vname,ndim,work)

  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ndim
  real(8),          intent(in) :: dvar(*)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: ier
  real(8), allocatable :: dtmp(:)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      allocate( dtmp(ndim), stat=ier ); call errall(ier,'unit_test__add_inp_real_array','dtmp')
      dtmp(1:ndim) = dvar(1:ndim)

      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'unit_test__add_inp_real_array',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(3D24.16)') dtmp

      close(iunit,status='keep')

      deallocate(dtmp)
    else
      call write_real_array(dvar,ndim,vname,'input')
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('real(8), allocatable, dimension(:)',var_type)
  call append('(' // trim(i2c(ndim)) // ')',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine unit_test__add_inp_real_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__add_out_int_array(vname,ndim)

  ! External variables
  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ndim

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('integer, allocatable, dimension(:)',var_type)
  call append('(' // trim(i2c(ndim)) // ')',var_dim)
  call append('output',var_info)
  call append('own_name',var_value)

  end subroutine unit_test__add_out_int_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__add_out_real_array(vname,ndim)

  ! External variables
  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ndim

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('real(8), allocatable, dimension(:)',var_type)
  call append('(' // trim(i2c(ndim)) // ')',var_dim)
  call append('output',var_info)
  call append('own_name',var_value)

  end subroutine unit_test__add_out_real_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_string_variable(svar,vname,work)

  character(len=*), intent(in) :: vname
  character(len=*), intent(in) :: svar
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: ier, n

  n = len(svar)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_string_variable',trim(gets(input_files,1)))

      write(iunit,'(a)') vname
      write(iunit,'(a)') svar

      close(iunit,status='keep')
    else
      call write_data(svar,vname)
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('character(len=' // trim(i2c(n)) // ')',var_type)
  call append('0',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_logical_variable(lvar,vname,work)

  character(len=*), intent(in) :: vname
  logical,          intent(in) :: lvar
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: ier

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_logical_variable',trim(gets(input_files,1)))

      write(iunit,'(a)') vname
      write(iunit,'(l1)') lvar

      close(iunit,status='keep')
    else
      call write_data(lvar,vname)
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('logical',var_type)
  call append('0',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_integer_variable(ivar,vname,work)

  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ivar
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: ier

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_integer_variable',trim(gets(input_files,1)))

      write(iunit,'(a)')   vname
      write(iunit,'(i12)') ivar

      close(iunit,status='keep')
    else
      call write_data(ivar,vname)
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('integer',var_type)
  call append('0',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_integer_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_integer_1d_array(ivar,vname,work)

  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ivar(:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, ier

  n1 = size(ivar,1)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_integer_1d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')    vname
      write(iunit,'(6i12)') ivar

      close(iunit,status='keep')
    else
      call write_integer_array(ivar,n1,vname)
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('integer, allocatable, dimension(:)',var_type)
  call append( '(' // trim(i2c(n1)) // ')' , var_dim )
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_integer_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_integer_2d_array(ivar,vname,work)

  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ivar(:,:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, n2, ier

  n1 = size(ivar,1)
  n2 = size(ivar,2)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_integer_2d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')    vname
      write(iunit,'(6i12)') ivar

      close(iunit,status='keep')
    else
      call write_integer_array(ivar,n1*n2,vname)
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('integer, allocatable, dimension(:,:)',var_type)
  call append( '(' // trim(i2c(n1)) // ',' // trim(i2c(n2)) // ')' , var_dim )
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_integer_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_variable(dvar,vname,work)

  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: ier

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_real_variable',trim(gets(input_files,1)))

      write(iunit,'(a)')      vname
      write(iunit,'(D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_data(dvar,vname)
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('real(8)',var_type)
  call append('0',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_real_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_1d_array(dvar,vname,work)

  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar(:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, ier

  n1 = size(dvar,1)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_real_1d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(3D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_real_array(dvar,n1,vname,'input')
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('real(8), allocatable, dimension(:)',var_type)
  call append( '(' // trim(i2c(n1)) // ')' , var_dim )
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_real_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_2d_array(dvar,vname,work)

  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar(:,:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, n2, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_real_2d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(3D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_real_array(dvar,n1*n2,vname,'input')
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('real(8), allocatable, dimension(:,:)',var_type)
  call append('('//trim(i2c(n1))//','//trim(i2c(n2))//')',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_real_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_3d_array(dvar,vname,work)

  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar(:,:,:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, n2, n3, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_real_3d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(3D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_real_array(dvar,n1*n2*n3,vname,'input')
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('real(8), allocatable, dimension(:,:,:)',var_type)
  call append('('//trim(i2c(n1))//','//trim(i2c(n2))//','//trim(i2c(n3))//')',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_real_3d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_4d_array(dvar,vname,work)

  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar(:,:,:,:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, n2, n3, n4, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_real_4d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(3D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_real_array(dvar,n1*n2*n3*n4,vname,'input')
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('real(8), allocatable, dimension(:,:,:,:)',var_type)
  call append('('//trim(i2c(n1))//','//trim(i2c(n2))//','//trim(i2c(n3))// &
              ','//trim(i2c(n4))//')',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_real_4d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_5d_array(dvar,vname,work)

  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar(:,:,:,:,:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, n2, n3, n4, n5, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)
  n5 = size(dvar,5)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_real_5d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(3D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_real_array(dvar,n1*n2*n3*n4*n5,vname,'input')
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('real(8), allocatable, dimension(:,:,:,:,:)',var_type)
  call append('('//trim(i2c(n1))//','//trim(i2c(n2))//','//trim(i2c(n3))// &
              ','//trim(i2c(n4))//','//trim(i2c(n5))//')',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_real_5d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_complex_variable(dvar,vname,work)

  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: ier

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_complex_variable',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(2D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_data(dvar,vname)
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('complex(8)',var_type)
  call append('0',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_complex_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_complex_1d_array(dvar,vname,work)

  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar(:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, ier

  n1 = size(dvar,1)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_complex_1d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(4D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_complex_array(dvar,n1,vname,'input')
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('complex(8), allocatable, dimension(:)',var_type)
  call append( '(' // trim(i2c(n1)) // ')' , var_dim )
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_complex_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_complex_2d_array(dvar,vname,work)

  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar(:,:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, n2, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_complex_2d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(4D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_complex_array(dvar,n1*n2,vname,'input')
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('complex(8), allocatable, dimension(:,:)',var_type)
  call append('('//trim(i2c(n1))//','//trim(i2c(n2))//')',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_complex_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_complex_3d_array(dvar,vname,work)

  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar(:,:,:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, n2, n3, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_complex_3d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(4D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_complex_array(dvar,n1*n2*n3,vname,'input')
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('complex(8), allocatable, dimension(:,:,:)',var_type)
  call append('('//trim(i2c(n1))//','//trim(i2c(n2))//','//trim(i2c(n3))//')',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_complex_3d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_complex_4d_array(dvar,vname,work)

  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar(:,:,:,:)
  logical, intent(in), optional :: work

  logical :: work_local
  integer :: n1, n2, n3, n4, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)

  work_local = .False.
  if( present(work) ) work_local = work
  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if( .not.work_local )then
    if(associate_status(input_files))then
      open(iunit,file=trim(gets(input_files,1)),status='unknown',form='formatted', &
           access='sequential',position='append',iostat=ier)
      call erropen(ier,'add_inp_complex_4d_array',trim(gets(input_files,1)))

      write(iunit,'(a)')       vname
      write(iunit,'(4D24.16)') dvar

      close(iunit,status='keep')
    else
      call write_complex_array(dvar,n1*n2*n3*n4,vname,'input')
    endif
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append('complex(8), allocatable, dimension(:,:,:,:)',var_type)
  call append('('//trim(i2c(n1))//','//trim(i2c(n2))//','//trim(i2c(n3))//','//trim(i2c(n4))//')',var_dim)
  if(     work_local) call append('work',var_info)
  if(.not.work_local) call append('input',var_info)
  call append('own_name',var_value)

  end subroutine add_inp_complex_4d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_string_variable(svar,vname)

  ! External variables
  character(len=*), intent(in) :: vname
  character(len=*), intent(in) :: svar

  ! Internal variables
  integer :: ier

  open(iunit,file=trim(adjustl(vname)),status='unknown',form='formatted', &
       access='sequential',iostat=ier)
  call erropen(ier,'write_string_variable',vname)
  rewind iunit

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', trim(adjustl(vname)), &
                                                 '(', trim(adjustl(vname)), ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a,a)')      '  character(len=*) :: ',trim(adjustl(vname))
  write(iunit,'(a)')
  write(iunit,'(a,a,a,a,a)')'  ', trim(adjustl(vname)), " = '", svar, "'"
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_logical_variable(lvar,vname)

  ! External variables
  character(len=*), intent(in) :: vname
  logical,          intent(in) :: lvar

  ! Internal variables
  integer :: ier

  open(iunit,file=trim(adjustl(vname)),status='unknown',form='formatted', &
       access='sequential',iostat=ier)
  call erropen(ier,'write_logical_variable',vname)
  rewind iunit

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', trim(adjustl(vname)), &
                                                 '(', trim(adjustl(vname)), ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a,a)')      '  logical :: ',trim(adjustl(vname))
  write(iunit,'(a)')
  if(lvar)then
    write(iunit,'(a,a,a)')  '  ', trim(adjustl(vname)), ' = .True.'
  else
    write(iunit,'(a,a,a)')  '  ', trim(adjustl(vname)), ' = .False.'
  endif
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_integer_variable(ivar,vname)

  ! External variables
  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ivar

  ! Internal variables
  integer :: ier

  open(iunit,file=trim(adjustl(vname)),status='unknown',form='formatted', &
       access='sequential',iostat=ier)
  call erropen(ier,'write_integer_variable',vname)
  rewind iunit

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', trim(adjustl(vname)), &
                                                 '(', trim(adjustl(vname)), ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a,a)')      '  integer :: ',trim(adjustl(vname))
  write(iunit,'(a)')
  write(iunit,'(2x,a,a,i0)') trim(adjustl(vname)), ' = ', ivar
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine write_integer_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_real_variable(dvar,vname)

  ! External variables
  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar

  ! Internal variables
  integer :: ier

  open(iunit,file=trim(adjustl(vname)),status='unknown',form='formatted', &
       access='sequential',iostat=ier)
  call erropen(ier,'write_real_variable',vname)
  rewind iunit

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', trim(adjustl(vname)), &
                                                 '(', trim(adjustl(vname)), ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a,a)')      '  real(8) :: ',trim(adjustl(vname))
  write(iunit,'(a)')
  write(iunit,'(2x,a,a,e24.16e3,a)') trim(adjustl(vname)), ' = ', dvar, '_8'
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine write_real_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_complex_variable(dvar,vname)

  ! External variables
  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar

  ! Internal variables
  integer :: ier

  open(iunit,file=trim(adjustl(vname)),status='unknown',form='formatted', &
       access='sequential',iostat=ier)
  call erropen(ier,'write_complex_variable',vname)
  rewind iunit

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', trim(adjustl(vname)), &
                                                 '(', trim(adjustl(vname)), ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a,a)')      '  complex(8) :: ',trim(adjustl(vname))
  write(iunit,'(a)')
  write(iunit,'(2x,a,a,e24.16e3,a,e24.16e3,a)') &
               trim(adjustl(vname)), ' = (', dreal(dvar), '_8,', dimag(dvar), '_8)'
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine write_complex_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_integer_array(ivar,ndim,vname)

  ! External variables
  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ndim, ivar(*)

  ! Internal variables
  integer :: i, ier

  open(iunit,file=trim(adjustl(vname)),status='unknown',form='formatted', &
       access='sequential',iostat=ier)
  call erropen(ier,'write_integer_array',vname)
  rewind iunit

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', trim(adjustl(vname)), &
                                                 '(', trim(adjustl(vname)), ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a)')        '  integer :: i'
  write(iunit,'(a,a)')      '  integer, dimension(*) :: ',trim(adjustl(vname))
  write(iunit,'(a)')
  write(iunit,'(a,i0)')     '  do i=1,',ndim
  write(iunit,'(a,a,a)')    '    ',trim(adjustl(vname)),'(i) = 0'
  write(iunit,'(a)')        '  enddo'
  write(iunit,'(a)')
  do i=1,ndim
    if(ivar(i) /= 0) write(iunit,'(2x,a,a,i0,a,i0)') trim(adjustl(vname)),'(',i,')=',ivar(i)
  enddo
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_real_array(dvar,ndim,vname,iotype)

  character(len=*), intent(in) :: vname, iotype
  integer,          intent(in) :: ndim
  real(8),          intent(in) :: dvar(*)

  character(len=max_length_of_string) :: str_tmp
  integer :: i, ier
  real(8) :: thr

  ! set threshold "thr"
  select case( iotype )

    case( 'input' )
      thr = 0.0d0

    case( 'output' )
      thr = 0.0d0
      if( key%filter )then
        if( is_equal( "'total'",    threshold_type, threshold_type%nr_records) .or. &
            is_equal( "'absolute'", threshold_type, threshold_type%nr_records) )then
          str_tmp = trim(gets(thresholds,thresholds%nr_records))
          read(str_tmp, *, iostat=ier) thr
          if( ier /= 0 )then
            write(6,*)
            write(6,*) ' ERROR while reading real variable from string. ier =',ier
            write(6,*) '       (write_real_array routine, unit_test_generator module)'
            write(6,*)
            stop
          endif
        endif
      endif

    case default
      write(6,*)
      write(6,*) ' ERROR in write_real_array routine in unit_test_generator module'
      write(6,*)
      stop

  end select

  ! use 10 times smaller threshold just to be sure
  thr = thr*0.1d0

  ! write real numbers to file over threshold "thr"
  open(iunit,file=trim(adjustl(vname)),status='unknown',form='formatted', &
       access='sequential',iostat=ier)
  call erropen(ier,'write_real_array',vname)
  rewind iunit

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', trim(adjustl(vname)), &
                                                 '(', trim(adjustl(vname)), ')'
  if( thr /= 0.0d0) &
    write(iunit,'(a,es8.1,a)') '  ! numbers with abs value below ',thr,' were not saved'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a)')        '  integer :: i'
  write(iunit,'(a,a)')      '  real(8), dimension(*) :: ',trim(adjustl(vname))
  write(iunit,'(a)')
  write(iunit,'(a,i0)')     '  do i=1,',ndim
  write(iunit,'(a,a,a)')    '    ',trim(adjustl(vname)),'(i) = 0.0d0'
  write(iunit,'(a)')        '  enddo'
  write(iunit,'(a)')
  do i=1,ndim
    if( dabs(dvar(i)) > thr ) write(iunit,'(2x,a,a,i0,a,e24.16e3,a)') &
                                    trim(adjustl(vname)), '(', i, ')=', dvar(i), '_8'
  enddo
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_complex_array(dvar,ndim,vname,iotype)

  character(len=*), intent(in) :: vname, iotype
  integer,          intent(in) :: ndim
  complex(8),       intent(in) :: dvar(*)

  character(len=max_length_of_string) :: str_tmp
  integer :: i, ier
  real(8) :: tmp, thr

  ! set threshold "thr"
  select case( iotype )

    case( 'input' )
      thr = 0.0d0

    case( 'output' )
      thr = 0.0d0
      if( key%filter )then
        if( is_equal( "'total'",    threshold_type, threshold_type%nr_records) .or. &
            is_equal( "'absolute'", threshold_type, threshold_type%nr_records) )then
          str_tmp = trim(gets(thresholds,thresholds%nr_records))
          read(str_tmp, *, iostat=ier) thr
          if( ier /= 0 )then
            write(6,*)
            write(6,*) ' ERROR while reading real variable from string. ier =',ier
            write(6,*) '       (write_real_array routine, unit_test_generator module)'
            write(6,*)
            stop
          endif
        endif
      endif

    case default
      write(6,*)
      write(6,*) ' ERROR in write_real_array routine in unit_test_generator module'
      write(6,*)
      stop

  end select

  ! use 10 times smaller threshold just to be sure
  thr = thr*0.1d0

  open(iunit,file=trim(adjustl(vname)),status='unknown',form='formatted', &
       access='sequential',iostat=ier)
  call erropen(ier,'write_complex_array',vname)
  rewind iunit

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', trim(adjustl(vname)), &
                                                 '(', trim(adjustl(vname)), ')'
  if( thr /= 0.0d0) &
    write(iunit,'(a,es8.1,a)') '  ! numbers with abs value below ',thr,' were not saved'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a)')        '  integer :: i'
  write(iunit,'(a,a)')      '  complex(8), dimension(*) :: ',trim(adjustl(vname))
  write(iunit,'(a)')
  write(iunit,'(a,i0)')     '  do i=1,',ndim
  write(iunit,'(a,a,a)')    '    ',trim(adjustl(vname)),'(i) = (0.0d0,0.0d0)'
  write(iunit,'(a)')        '  enddo'
  write(iunit,'(a)')
  do i=1,ndim
    tmp = dreal(dvar(i))*dreal(dvar(i)) + dimag(dvar(i))*dimag(dvar(i))
    tmp = dsqrt(tmp)
    if(tmp > thr) write(iunit,'(2x,a,a,i0,a,e24.16e3,a,e24.16e3,a)') &
      trim(adjustl(vname)), '(', i, ')=(', dreal(dvar(i)), '_8,', dimag(dvar(i)), '_8)'
  enddo
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_data_routine(iw,vname)

  ! External variables
  character(len=*), intent(in) :: vname
  integer,          intent(in) :: iw

  ! Internal variables
  integer :: i, io, ier
  character(len=250) :: line

  io = iw + 1
  open(io,file=vname,status='old',form='formatted',access='sequential',iostat=ier)
  call erropen(ier,'add_data_routine',vname)
  rewind io

  do

  read(io,'(a)',iostat=ier) line
  if(ier == -1) exit
  write(iw,'(a)') trim(line)

  enddo

  close(io,status='delete')

  end subroutine add_data_routine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  integer function find_name(string) result(ipos)
    ! find position of string in var_name list
    ! if there is no such string return 0
    character(len=*), intent(in) :: string
    integer :: i
    ipos = 0
    do i=1,get_size(var_name)
      if(string == trim(gets(var_name,i)))then
        ipos = i
        exit
      endif
    enddo
  end function

  logical function is_equal_type_type(type_list1,ir1,type_list2,ir2)
    ! compare string with record "ir" in type_list
    integer, intent(in) :: ir1, ir2
    type( list ), intent(in) :: type_list1, type_list2
    logical :: lresult
    lresult = string__compare_strings(type_list1%records(ir1),type_list2%records(ir2))
    is_equal_type_type = lresult
  end function

  logical function is_equal_string_type(string,type_list,ir)
    ! compare string with record "ir" in type_list
    integer, intent(in) :: ir
    type( list ), intent(in) :: type_list
    character(len=*), intent(in) :: string
    is_equal_string_type = string__compare_strings(string,type_list%records(ir))
  end function

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_compare_info(irecord,ndim,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  integer,          intent(in) :: irecord, ndim

  call append(compare,compare_info)
  select case(compare)
    case('number-number')
      call append('  call compare_2numbers(', compare_technique)
      call append(trim(adjustl(gets(var_name,irecord))), compare_name)
      call append(trim(gets(var_dim ,irecord)), compare_dim )
      call append(trim(gets(var_type,irecord)), compare_type)

    case('array-array')
      call append('  call compare_2vectors(' // trim(i2c(ndim)), compare_technique)
      call append(trim(adjustl(gets(var_name,irecord))), compare_name)
      call append(trim(gets(var_dim ,irecord)), compare_dim )
      call append(trim(gets(var_type,irecord)), compare_type)

    case('file-array')
      call append('  call compare_vector', compare_technique)
      call append(trim(adjustl(gets(var_name,irecord))), compare_name)
      call append(trim(gets(var_dim ,irecord)), compare_dim )
      call append(trim(gets(var_type,irecord)), compare_type)
      if(.not.present(file))then
        write(6,*) 'ERROR: You forgot to define "file" name unit_test__add_out_variable call'
        stop
      endif
      ! add file (with suffix) to the list of output files (output_files)
      if(associate_status(suffix))then
        call append(file // '_' // trim(gets(suffix,1)), output_files)
      else
        call append(file, output_files)
      endif

    case('file-file')
      call append('  call compare_files', compare_technique)
      call append('0', compare_name)
      call append('0', compare_dim )
      call append('0', compare_type)
      if(.not.present(file))then
        write(6,*) 'ERROR: You forgot to define "file" name in unit_test__add_out_variable call'
        stop
      endif
      ! add file (with suffix) to the list of output files (output_files)
      if(associate_status(suffix))then
        call append(file // '_' // trim(gets(suffix,1)), output_files)
      else
        call append(file, output_files)
      endif

    case default
      write(6,*) 'ERROR: Wrong optional parameter "compare" in unit_test__add_out_variable call'
      write(6,*) '       Use string "file-array" or "file-file"'
      stop
  end select

  ! set threshold
  if(present(threshold))then
    call append(threshold, thresholds)
  else
    call append('1.0d-15', thresholds)
  endif

  ! set what will be compared against threshold
  if(present(check))then
    call append("'" // check // "'", threshold_type)
  else
    call append("'total'", threshold_type)
  endif

  end subroutine add_compare_info

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_integer_variable(ivar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ivar

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('integer',var_type)
    call append('0',var_dim)
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),0,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),0,'array-array',threshold,file,check)
    call write_data(ivar,trim(adjustl(vname)) // '_ref')
    write(6,*) ' STOP: output integer data are not yet implemented'
    stop
  endif

  end subroutine add_out_integer_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_integer_1d_array(ivar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ivar(:)

  ! Internal variables
  integer :: n1

  n1 = size(ivar,1)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('integer, allocatable, dimension(:)',var_type)
    call append( '(' // trim(i2c(n1)) // ')' , var_dim )
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1,'array-array',threshold,file,check)
    call write_integer_array(ivar,n1,trim(adjustl(vname)) // '_ref')
    write(6,*) ' STOP: output integer data are not yet implemented'
    stop
  endif

  end subroutine add_out_integer_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_integer_2d_array(ivar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  integer,          intent(in) :: ivar(:,:)

  ! Internal variables
  integer :: n1, n2

  n1 = size(ivar,1)
  n2 = size(ivar,2)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('integer, allocatable, dimension(:,:)',var_type)
    call append( '(' // trim(i2c(n1)) // ',' // trim(i2c(n2)) // ')' , var_dim )
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1*n2,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1*n2,'array-array',threshold,file,check)
    call write_integer_array(ivar,n1*n2,trim(adjustl(vname)) // '_ref')
    write(6,*) ' STOP: output integer data are not yet implemented'
    stop
  endif

  end subroutine add_out_integer_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_variable(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('real(8)',var_type)
    call append('0',var_dim)
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),0,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),0,'number-number',threshold,file,check)
    call write_data(dvar,trim(adjustl(vname)) // '_ref')
  endif

  end subroutine add_out_real_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_1d_array(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar(:)

  ! Internal variables
  integer :: n1

  n1 = size(dvar,1)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('real(8), allocatable, dimension(:)',var_type)
    call append( '(' // trim(i2c(n1)) // ')' , var_dim )
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1,'array-array',threshold,file,check)
    call write_real_array(dvar,n1,trim(adjustl(vname)) // '_ref','output')
  endif

  end subroutine add_out_real_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_2d_array(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar(:,:)

  ! Internal variables
  integer :: n1, n2

  n1 = size(dvar,1)
  n2 = size(dvar,2)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('real(8), allocatable, dimension(:,:)',var_type)
    call append('('//trim(i2c(n1))//','//trim(i2c(n2))//')',var_dim)
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1*n2,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1*n2,'array-array',threshold,file,check)
    call write_real_array(dvar,n1*n2,trim(adjustl(vname)) // '_ref','output')
  endif

  end subroutine add_out_real_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_3d_array(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar(:,:,:)

  ! Internal variables
  integer :: n1, n2, n3

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('real(8), allocatable, dimension(:,:,:)',var_type)
    call append('('//trim(i2c(n1))//','//trim(i2c(n2))//','//trim(i2c(n3))//')',var_dim)
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1*n2*n3,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1*n2*n3,'array-array',threshold,file,check)
    call write_real_array(dvar,n1*n2*n3,trim(adjustl(vname)) // '_ref','output')
  endif

  end subroutine add_out_real_3d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_4d_array(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar(:,:,:,:)

  ! Internal variables
  integer :: n1, n2, n3, n4

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('real(8), allocatable, dimension(:,:,:,:)',var_type)
    call append('('//trim(i2c(n1))//','//trim(i2c(n2))//','//trim(i2c(n3))// &
                ','//trim(i2c(n4))//')',var_dim)
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1*n2*n3*n4,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1*n2*n3*n4,'array-array',threshold,file,check)
    call write_real_array(dvar,n1*n2*n3*n4,trim(adjustl(vname)) // '_ref','output')
  endif

  end subroutine add_out_real_4d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_5d_array(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  real(8),          intent(in) :: dvar(:,:,:,:,:)

  ! Internal variables
  integer :: n1, n2, n3, n4, n5

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)
  n5 = size(dvar,5)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('real(8), allocatable, dimension(:,:,:,:,:)',var_type)
    call append('('//trim(i2c(n1))//','//trim(i2c(n2))//','//trim(i2c(n3))// &
                ','//trim(i2c(n4))//','//trim(i2c(n5))//')',var_dim)
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1*n2*n3*n4*n5,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1*n2*n3*n4*n5,'array-array',threshold,file,check)
    call write_real_array(dvar,n1*n2*n3*n4*n5,trim(adjustl(vname)) // '_ref','output')
  endif

  end subroutine add_out_real_5d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_complex_variable(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('complex(8)',var_type)
    call append('0',var_dim)
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),0,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),0,'number-number',threshold,file,check)
    call write_data(dvar,trim(adjustl(vname)) // '_ref')
  endif

  end subroutine add_out_complex_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_complex_1d_array(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar(:)

  ! Internal variables
  integer :: n1

  n1 = size(dvar,1)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('complex(8), allocatable, dimension(:)',var_type)
    call append( '(' // trim(i2c(n1)) // ')' , var_dim )
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1,'array-array',threshold,file,check)
    call write_complex_array(dvar,n1,trim(adjustl(vname)) // '_ref','output')
  endif

  end subroutine add_out_complex_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_complex_2d_array(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar(:,:)

  ! Internal variables
  integer :: n1, n2

  n1 = size(dvar,1)
  n2 = size(dvar,2)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('complex(8), allocatable, dimension(:,:)',var_type)
    call append('('//trim(i2c(n1))//','//trim(i2c(n2))//')',var_dim)
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1*n2,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1*n2,'array-array',threshold,file,check)
    call write_complex_array(dvar,n1*n2,trim(adjustl(vname)) // '_ref','output')
  endif

  end subroutine add_out_complex_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_complex_3d_array(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar(:,:,:)

  ! Internal variables
  integer :: n1, n2, n3

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('complex(8), allocatable, dimension(:,:,:)',var_type)
    call append('('//trim(i2c(n1))//','//trim(i2c(n2))//','//trim(i2c(n3))//')',var_dim)
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1*n2*n3,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1*n2*n3,'array-array',threshold,file,check)
    call write_complex_array(dvar,n1*n2*n3,trim(adjustl(vname)) // '_ref','output')
  endif

  end subroutine add_out_complex_3d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_complex_4d_array(dvar,vname,compare,threshold,file,check)

  ! External variables
  character(len=*), intent(in), optional :: compare, threshold, file, check
  character(len=*), intent(in) :: vname
  complex(8),       intent(in) :: dvar(:,:,:,:)

  ! Internal variables
  integer :: n1, n2, n3, n4

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call append(vname,var_name)
    call append(trim(adjustl(current_routine)),var_routine)
    call append('complex(8), allocatable, dimension(:,:,:,:)',var_type)
    call append('('//trim(i2c(n1))//','//trim(i2c(n2))//','//trim(i2c(n3))//','//trim(i2c(n4))//')',var_dim)
    call append('output',var_info)
    call append('own_name',var_value)
  else
    if( is_equal('work',var_info,find_name(vname)) )then
      call change('output',var_info,find_name(vname))
    else
      call change('inout',var_info,find_name(vname))
    endif
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(find_name(vname),n1*n2*n3*n4,compare,threshold,file,check)
  else
    call add_compare_info(find_name(vname),n1*n2*n3*n4,'array-array',threshold,file,check)
    call write_complex_array(dvar,n1*n2*n3*n4,trim(adjustl(vname)) // '_ref','output')
  endif

  end subroutine add_out_complex_4d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__add_opt_variable(vname,vvalue,vtype)
  character(len=*), intent(in) :: vname, vvalue, vtype

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call append(vname,var_name)
  call append(trim(adjustl(current_routine)),var_routine)
  call append(vtype,var_type)
  call append('0',var_dim)
  call append('optional',var_info)
  call append(vvalue,var_value)

  end subroutine unit_test__add_opt_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__set_routine_name(drv_name)
    character(len=*), intent(in) :: drv_name
    ! add name of driver routine to the list
    call append(trim(adjustl(drv_name)),driver)
    current_routine = trim(adjustl(drv_name))
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__add_input_file(file_name)
  character(len=*), intent(in) :: file_name

  ! add file_name (with suffix) to the list of input files (input_files)
  if(associate_status(suffix))then
    call append(file_name // '_' // trim(gets(suffix,1)), input_files)
  else
    call append(file_name, input_files)
  endif

  end subroutine unit_test__add_input_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__add_output_file(file_name)
  character(len=*), intent(in) :: file_name

  ! add file_name (with suffix) to the list of output files (output_files)
  if(associate_status(suffix))then
    call append(file_name // '_' // trim(gets(suffix,1)), output_files)
  else
    call append(file_name, output_files)
  endif

  end subroutine unit_test__add_output_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__add_module(string)
  character(len=*), intent(in) :: string

  ! add module to the test
  call append(string,modules)

  end subroutine unit_test__add_module

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__set_keys(mod_call,filter)
    logical, intent(in), optional :: mod_call, filter
    ! set keys
    if(present(mod_call)) key%mod_call = mod_call
    if(present(filter  )) key%filter   = filter
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__redefine_variable(variable_type,variable_size,nvar)

  character(len=*), intent(in), optional :: variable_type, variable_size
  integer,          intent(in), optional :: nvar
  integer :: i

  ! position of the string
  i = get_size(var_name)
  if(present(nvar)) i = nvar

  ! reset
  if(present(variable_type)) call change(variable_type,var_type,i)
  if(present(variable_size)) call change(variable_size,var_dim,i)

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine set_suffix_string(string)
    character(len=*), intent(in) :: string
    ! add string to the name of the test source file
    if(associate_status(suffix))then
      call change(trim(gets(suffix,1)) // '_' // string, suffix, 1)
    else
      call append(string, suffix)
    endif
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine set_suffix_integer(ivar)
    integer, intent(in) :: ivar
    ! add integer to the name of the test source file
    if(associate_status(suffix))then
      call change(trim(gets(suffix,1)) // '_' // trim(i2c(ivar)), suffix, 1)
    else
      call append(trim(i2c(ivar)), suffix)
    endif
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine set_suffix_logical(lvar)
    logical, intent(in) :: lvar
    ! add integer to the name of the test source file
    if(associate_status(suffix))then
      call change(trim(gets(suffix,1)) // '_' // l2c(lvar), suffix, 1)
    else
      call append(l2c(lvar), suffix)
    endif
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  function unit_test__get_output_name(nout) result(out_name)
    integer, intent(in) :: nout
    character(len=max_length_of_string) :: out_name
    integer :: i
    if( max_length_of_string < get_size(output_files,nout) )then
      write(6,*) ' ERROR: Maximum allowed length of string exceeded in ', &
                             'unit_test_generator module (unit_test__get_output_name function)'
      write(6,*) '        max_string_length = ',max_length_of_string
      write(6,*) '        string length     = ',get_size(output_files,nout)
      write(0,*) ' ERROR in unit_test_generator module (unit_test__get_output_name function)'
      stop
    endif
    do i=1,max_length_of_string
      out_name(i:i) = ' '
    enddo
    do i=1,get_size(output_files,nout)
      out_name(i:i) = get_char(output_files,nout,i)
    enddo
  end function

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine errall(ier,chsub,charr)
    integer          :: ier
    character(len=*) :: chsub, charr
    if(ier /= 0)then
      write(6,*)
      write(6,'(a,a,a,a)')' Error in allocation of ',charr,' in subroutine ',chsub
      write(6,'(a,i4)'   )' iostat =',ier
      stop 'ERROR in allocation (see output for more information)'
    endif
  end subroutine errall

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine erropen(ier,chsub,chfile)
    integer          :: ier
    character(len=*) :: chsub, chfile
    if(ier /= 0)then
      write(6,*)
      write(6,'(a,a,a,a)')' ERROR while opening ',chfile,' in subroutine ',chsub
      write(6,'(a,i4)'   )' iostat =',ier
      stop 'ERROR in opening file (see output for more information)'
    endif
  end subroutine erropen

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__generate_test()
  ! generate unit test program with name test.F90_sufix

  ! Internal variables
  character(len=100) :: tmpf
  character(len=1)   :: coma
  integer :: i, j, ier, itmp1, itmp2
  logical :: file_exist, lll

  ! check if the module is initialized
  if(.not. iam_initialized)then
    write(6,'(/a/)') ' ERROR: Module unit_test__generator was not initialized!'
    stop
  endif

  ! set name of the test source file
  if(associate_status(suffix))then
    call append('test.F90_' // trim(gets(suffix,1)), test_name)
  else
    call append('test.F90', test_name)
  endif

  inquire(file=trim(gets(test_name,1)), exist=file_exist)
  if(file_exist)then
    write(6,'(a)')
    write(6,'(a,a,a)') ' WARNING: file ',trim(gets(test_name,1)),' already exist'
    write(6,'(a)')     ' WARNING:   => it will be overwritten'
    write(6,'(a)')
  endif

  open(iunit,file=trim(gets(test_name,1)),status='unknown',form='formatted', &
       position='rewind',iostat=ier)
  if( ier /= 0 )then
    write(6,'(a)')
    write(6,'(a,a,a)')' ERROR while opening ',trim(gets(test_name,1)), &
                      ' file in unit_test__generate_test routine'
    write(6,'(a)'    )' Please, check its existence, format and length'
    write(6,'(a,i4)' )' iostat =',ier
    write(6,'(a)')
    stop ' ERROR while opening some file (see output file for more information)'
  endif

  ! header
  write(iunit,'(a)')'program test'
  write(iunit,'(a)')
  write(iunit,'(a)')'  use unit_testing'
  do i=1,get_size(modules)
    write(iunit,'(a,100a)')'  use ', trim(gets(modules,i))
  enddo
  write(iunit,'(a)')'  implicit none'
  write(iunit,'(a)')
  write(iunit,'(a)')'  integer :: ier'

  ! ---------------------------------------------------------------------------- !
  ! declaration of variables
  do i=1,get_size(var_name)
    write(tmpf,'(a,i0,a,i0,a)') '(2x,a', get_size(var_type,i), ',a,a', get_size(var_name,i), ')'
    write(iunit,tmpf) gets(var_type,i), ' :: ', gets(var_name,i)
  enddo
  ! declare reference variables
  do i=1,get_size(compare_info)
    if(trim(gets(compare_info,i)) == 'array-array' .or. &
       trim(gets(compare_info,i)) == 'number-number')then
      itmp1 = get_size(compare_type,i)
      itmp2 = get_size(compare_name,i) + 4
      write(tmpf,'(a,i0,a,i0,a)') '(2x,a', itmp1, ',a,a', itmp2, ')'
      write(iunit,tmpf) gets(compare_type,i), ' :: ', trim(gets(compare_name,i)) // '_ref'
    endif
  enddo

  ! ---------------------------------------------------------------------------- !
  ! allocate arrays
  write(iunit,'(a)')
  do i=1,get_size(var_name)
    if( get_char(var_dim,i,1) /= '0')then
      itmp1 = get_size(var_name,i)
      itmp2 = get_size(var_dim,i)
      write(tmpf,'(a,i0,a,i0,a,i0,a)') '(a,a', itmp1, ',a', itmp2, ',a,a', itmp1, ',a)'
      write(iunit,tmpf) '  allocate( ', gets(var_name,i), gets(var_dim,i), &
                        " , stat=ier ); call errall(ier,'test','", gets(var_name,i), "')"
    endif
  enddo
  ! allocate reference array if we want to use array-array comparison
  do i=1,get_size(compare_info)
    if(trim(gets(compare_info,i)) == 'array-array')then
      itmp1 = get_size(compare_name,i) + 4
      itmp2 = get_size(compare_dim,i)
      write(tmpf,'(a,i0,a,i0,a,i0,a)') '(a,a', itmp1, ',a', itmp2, ',a,a', itmp1, ',a)'
      write(iunit,tmpf) '  allocate( ', trim(gets(compare_name,i)) // '_ref',           &
                        gets(compare_dim,i), " , stat=ier ); call errall(ier,'test','", &
                        trim(gets(compare_name,i)) // '_ref', "')"
    endif
  enddo

  ! ---------------------------------------------------------------------------- !
  ! open input file (if data are in the file)
  if(associate_status(input_files))then
    if(get_size(input_files) > 1)then
      write(6,*) 'ERROR: For now the program works only with one input file.'
      stop
    endif
    write(iunit,'(a)')
    write(iunit,'(a,70a,a)') "  open(60,file='", trim(gets(input_files,1)), &
                             "',status='old',form='formatted',access='sequential',iostat=ier)"
    write(iunit,'(a)')       '  if(ier /= 0)then'
    write(iunit,'(a,70a,a)') "    write(6,*) 'test: ERROR while creating ", &
                             trim(gets(input_files,1)), " file'"
    write(iunit,'(a)')       '    stop'
    write(iunit,'(a)')       '  endif'
  endif

  ! ---------------------------------------------------------------------------- !
  ! get input data
  write(iunit,'(a)')
  do i=1,get_size(var_name)
    if( trim(gets(var_info,i)) == 'input' .or. trim(gets(var_info,i)) == 'inout' )then
      if(associate_status(input_files))then
        write(tmpf,'(a,i0,a)') '(a,a', get_size(var_name,i), ')'
        write(iunit,'(a)')'  read(60,*)'
        write(iunit,tmpf) '  read(60,*) ', gets(var_name,i)
      else
        itmp1 = get_size(var_name,i)
        write(tmpf,'(a,i0,a,i0,a)') '(a,a', itmp1, ',a,a', itmp1, ',a)'
        write(iunit,tmpf) '  call get_data_', gets(var_name,i), '(', gets(var_name,i), ')'
      endif
    endif
  enddo
  ! close input file
  if(associate_status(input_files))then
    write(iunit,'(a)')
    write(iunit,'(a)')"  close(60,status='keep')"
  endif

  ! ---------------------------------------------------------------------------- !
  ! set optional parameters
  write(iunit,'(a)')
  do i=1,get_size(var_name)
    if( trim(gets(var_info,i)) == 'optional' )then
      write(iunit,'(2x,70a)') trim(gets(var_name,i)), ' = ', trim(gets(var_value,i))
    endif
  enddo

  ! ---------------------------------------------------------------------------- !
  ! call driver or initialization routines
  do i=1,get_size(driver)
    write(iunit,'(a)')
    write(iunit,'(a,70a,a)')"  call ", trim(gets(driver,i)), '( &'
    coma = ' '
    if(key%mod_call)then
      do j=1,get_size(var_name)
        if( .not.is_equal(var_routine,j,driver,i) ) cycle
        itmp1 = get_size(var_name,j)
        write(tmpf,'(a,i0,a,i0,a)') '(9x,a,a', itmp1, ',a,a', itmp1, ',a)'
        write(iunit,tmpf) coma, gets(var_name,j), ' = ', gets(var_name,j), ' &'
        coma = ','
      enddo
    else
      do j=1,get_size(var_name)
        if( .not.is_equal(var_routine,j,driver,i) ) cycle
        write(tmpf,'(a,i0,a)') '(9x,a,a', get_size(var_name,j), ',a)'
        write(iunit,tmpf) coma, gets(var_name,j), ' &'
        coma = ','
      enddo
    endif
    write(iunit,'(9x,a)') ')'
  enddo

  ! ---------------------------------------------------------------------------- !
  ! get output data
  do i=1,get_size(compare_info)
    if(trim(gets(compare_info,i)) == 'array-array' .or. &
       trim(gets(compare_info,i)) == 'number-number')then
      write(iunit,'(a)')
      write(iunit,'(a,70a,a,70a,a)') '  call get_data_', trim(gets(compare_name,i)), '_ref(', &
                                                         trim(gets(compare_name,i)), '_ref)'
    endif
  enddo
  ! compare results
  write(iunit,'(a)')
  do i=1,get_size(compare_technique)
    if(trim(gets(compare_info,i)) == 'file-file')then
      write(iunit,'(70a,a,70a,a,70a,a,70a,a,70a,a,70a,a)') &
           trim(gets(compare_technique,i)), "(60,'", trim(gets(output_files,i)), "','", &
           trim(gets(output_files,i)), "_ref',", trim(gets(thresholds,i)), ',check=',   &
           trim(gets(threshold_type,i)), ",comparing='",trim(gets(output_files,i)), "')"
    elseif(trim(gets(compare_info,i)) == 'file-array')then
      write(iunit,'(70a,a,70a,a,70a,a,70a,a,70a,a,70a,a)') &
           trim(gets(compare_technique,i)), "(60,'", trim(gets(output_files,i)), "',", &
           trim(gets(compare_name,i)), ",", trim(gets(thresholds,i)), ',check=',       &
           trim(gets(threshold_type,i)), ",comparing='",trim(gets(compare_name,i)), "')"
    elseif(trim(gets(compare_info,i)) == 'array-array')then
      write(iunit,'(70a,a,70a,a,70a,a,70a,a,70a,a,70a,a)') &
           trim(gets(compare_technique,i)), ',', trim(gets(compare_name,i)), '_ref,', &
           trim(gets(compare_name,i)), ',', trim(gets(thresholds,i)), ',check=',      &
           trim(gets(threshold_type,i)), ",comparing='",trim(gets(compare_name,i)), "')"
    elseif(trim(gets(compare_info,i)) == 'number-number')then
      write(iunit,'(70a,70a,a,70a,a,70a,a,70a,a,70a,a)') &
           trim(gets(compare_technique,i)), trim(gets(compare_name,i)), '_ref,', &
           trim(gets(compare_name,i)), ',', trim(gets(thresholds,i)), ',check=', &
           trim(gets(threshold_type,i)), ",comparing='",trim(gets(compare_name,i)), "')"
    endif
  enddo

  ! ---------------------------------------------------------------------------- !
  ! deallocation
  write(iunit,'(a)')
  do i=1,get_size(var_name)
    if( get_char(var_dim,i,1) /= '0')then
      write(tmpf,'(a,i0,a)') '(a,a', get_size(var_name,i), ',a)'
      write(iunit,tmpf) "  deallocate( ", gets(var_name,i), " )"
    endif
  enddo
  ! deallocate reference arrays
  do i=1,get_size(compare_info)
    if(trim(gets(compare_info,i)) == 'array-array')then
      write(iunit,'(a,70a,a)')"  deallocate( ", trim(gets(compare_name,i)), "_ref )"
    endif
  enddo
  ! ---------------------------------------------------------------------------- !
  ! test finished successfully
  write(iunit,'(a)')
  write(iunit,'(a)')"  if(all_tests_passed()) write(6,*) ' test result: =====  OK  ====='"

  ! ---------------------------------------------------------------------------- !
  ! write error routine
  write(iunit,'(a)')
  write(iunit,'(a)')"contains"
  write(iunit,'(a)')
  write(iunit,'(a)')"  subroutine errall(ier,chsub,charr)"
  write(iunit,'(a)')"  implicit none"
  write(iunit,'(a)')
  write(iunit,'(a)')"  integer          :: ier"
  write(iunit,'(a)')"  character(len=*) :: chsub, charr"
  write(iunit,'(a)')
  write(iunit,'(a)')"  if(ier /= 0)then"
  write(iunit,'(a)')"    write(6,*)"
  write(iunit,'(a)')"    write(6,'(a,a,a,a)')' Error in allocation of ',charr,' in subroutine ',chsub"
  write(iunit,'(a)')"    write(6,'(a,i4)'   )' iostat =',ier"
  write(iunit,'(a)')"    stop 'ERROR in allocation (see output form more information)'"
  write(iunit,'(a)')"  endif"
  write(iunit,'(a)')
  write(iunit,'(a)')"  end subroutine"

  ! ---------------------------------------------------------------------------- !
  ! add input data routines
  if(.not. associate_status(input_files))then
    do i=1,get_size(var_name)
      if( trim(gets(var_info,i)) == 'input' .or. trim(gets(var_info,i)) == 'inout' )then
        write(iunit,'(a)')
        call add_data_routine(iunit,trim(gets(var_name,i)))
      endif
    enddo
  endif
  ! add output data routines
  do i=1,get_size(compare_info)
    if(trim(gets(compare_info,i)) == 'array-array' .or. &
       trim(gets(compare_info,i)) == 'number-number')then
      write(iunit,'(a)')
      call add_data_routine(iunit,trim(gets(compare_name,i)) // '_ref')
    endif
  enddo
  ! end of the test program
  write(iunit,'(a)')"end program"

  ! ---------------------------------------------------------------------------- !
  ! finilize creation of unit test
  close(iunit,status='keep')
  call unit_test__destructor()

  end subroutine unit_test__generate_test

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unit_test__constructor()
    ! set maximum length of string for string_manipulation module
    call string__set_max_length(max_length_of_string)
    ! set output unit for string_manipulation module
    call string__set_iunit(6)
    ! initialize unit_test_generator
    ! remember that the module was already initialized
    if(iam_initialized)then
      write(6,'(a)')
      write(6,'(a)')'  ERROR: Module unit_test_generator is already initialized!'
      write(6,'(a)')'         Use different enviroment variable for this unit test.'
      write(6,'(a)')'         This stop will prevent simultaneous creating of unit tests.'
      write(6,'(a)')
      stop
    endif
    iam_initialized = .True.
    ! some compilers will not do it --> can cause problems in append routine
    call init( input_files ); call init( output_files )
    call init( driver ); call init( test_name )
    call init( var_name ); call init( var_type ); call init( var_dim )
    call init( var_info ); call init( var_value ); call init( var_routine )
    call init( modules ); call init( suffix )
    call init( compare_technique ); call init( compare_info ); call init( thresholds )
    call init( threshold_type )
    call init( compare_name ); call init( compare_type ); call init( compare_dim )
  end subroutine

  subroutine unit_test__destructor()
    ! job is done
    iam_initialized = .False.
    ! deallocate all global pointers
    call free( input_files ); call free( output_files )
    call free( driver ); call free( test_name )
    call free( var_name ); call free( var_type ); call free( var_dim )
    call free( var_info ); call free( var_value ); call free( var_routine )
    call free( modules ); call free( suffix )
    call free( compare_technique ); call free( compare_info ); call free( thresholds )
    call free( threshold_type )
    call free( compare_name ); call free( compare_type ); call free( compare_dim )
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!                      SET OF OPERATIONS ON TYPE LIST                            !
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  pure function associate_status(type_list) result(l)
    type( list ), intent(in) :: type_list
    logical :: l
    l = .True.
    if(type_list%nr_records == 0) l = .False.
  end function

  function gets(type_list,ir) result(c)
    type( list ), intent(in) :: type_list
    integer,      intent(in) :: ir
    character(len=max_length_of_string) :: c
    integer :: i
    if( max_length_of_string < get_size(type_list,ir) )then
      write(6,*) ' ERROR: Maximum allowed length of string exceeded in ', &
                             'unit_test_generator module (gets function)'
      write(6,*) '        max_string_length = ',max_length_of_string
      write(6,*) '        string length     = ',get_size(type_list,ir)
      write(0,*) ' ERROR in unit_test_generator module (gets function)'
      stop
    endif
    do i=1,max_length_of_string
      c(i:i) = ' '
    enddo
    c = trim(string__gets(type_list%records(ir)))
  end function

  subroutine append(add_str,type_list)
    type( list ),     intent(inout) :: type_list
    character(len=*), intent(in)    :: add_str
    type_list%nr_records = type_list%nr_records + 1
    if(type_list%nr_records > max_nr_of_records)then
      write(6,*)
      write(6,*) ' ERROR in unit_test_generator.f90'
      write(6,*) '       Maximum number of records exceeded! '
      write(6,*) '       Enlarge max_nr_of_records in unit_test_generator.f90 module.'
      write(6,*) '       You tried add string --> ',add_str
      write(6,*)
      stop
    endif
    call string__append(add_str,type_list%records(type_list%nr_records))
  end subroutine

  subroutine change(add_str,type_list,ir)
    type( list ),     intent(inout) :: type_list
    character(len=*), intent(in)    :: add_str
    integer,          intent(in)    :: ir
    call string__append(add_str,type_list%records(ir))
  end subroutine

  pure function get_size(type_list,ir) result(i)
    type( list ), intent(in) :: type_list
    integer, intent(in), optional :: ir
    integer :: i
    if(     present(ir)) i = string__get_size(type_list%records(ir))
    if(.not.present(ir)) i = type_list%nr_records
  end function

  function get_char(type_list,ir,ic) result(c)
    type( list ), intent(in) :: type_list
    integer,      intent(in) :: ir, ic
    character :: c
    c = string__get_char(type_list%records(ir),ic)
  end function

  subroutine init(type_list)
    type( list ), intent(inout) :: type_list
    integer :: i
    do i=1,max_nr_of_records
      call string__init(type_list%records(i))
    enddo
    type_list%nr_records = 0
  end subroutine

  subroutine free(type_list)
    type( list ), intent(inout) :: type_list
    integer :: i
    do i=1,max_nr_of_records
      call string__free(type_list%records(i))
    enddo
    type_list%nr_records = 0
  end subroutine

END MODULE UNIT_TEST_GENERATOR
