!================================================================================!
!                                                                                !
!  Program:      ReSpect                                                         !
!                                                                                !
!  Module:       UNIT_TESTING                                                    !
!                                                                                !
!  Description:  - compare complex/real arrays or text files                     !
!                - write some statistics on output                               !
!                - this module was designed to be used by unit_test_generator    !
!                  module, but it can be used also standalone                    !
!                                                                                !
!  Dependencies: None                                                            !
!                                                                                !
!  Author:       Stanislav Komorovsky                                            !
!                                                                                !
!  Date:         May - July 2012 (Toulouse, Tromso), ...                         !
!                                                                                !
!================================================================================!

MODULE UNIT_TESTING

  implicit none

  public all_tests_passed    ! logical function; return True if all tests passed
  public compare_2numbers    ! compare two real numbers
  public compare_vector      ! compare two vectors;  input: from parameters and from file
  public compare_2vectors    ! compare two vectors;  input: both vectors are from parameters
  public compare_files       ! compare two sets of data from two files
  public compare_text_files  ! compare two text files
  public unpackme            ! unpack vector to matrix

  private

  interface compare_2numbers
    module procedure compare_real
    module procedure compare_complex
  end interface compare_2numbers

  interface compare_2vectors
    module procedure compare_1d_real_vectors
    module procedure compare_2d_real_vectors
    module procedure compare_3d_real_vectors
    module procedure compare_4d_real_vectors
    module procedure compare_5d_real_vectors
    module procedure compare_1d_complex_vectors
    module procedure compare_2d_complex_vectors
    module procedure compare_3d_complex_vectors
    module procedure compare_4d_complex_vectors
  end interface compare_2vectors

  integer, private, parameter :: file_unit = 137
  logical, private            :: file_exists
  real(8), private            :: diff_threshold = 1.0d-50

  private :: init_info, print_err_info, compare_real_vectors, compare_complex_vectors

  logical, save :: test_status = .True.  ! remember status of tests
  
  type, private :: info
    integer :: ival     ! position of val1 and val2 numbers
    real(8) :: val1     ! variable from first set
    real(8) :: val2     ! variable from second set
    real(8) :: reld     ! relative difference of val1 and val2
    real(8) :: absd     ! absolute difference of val1 and val2
    real(8) :: val1_img ! variable from first set  (complex part)
    real(8) :: val2_img ! variable from second set (complex part)
    logical :: complx   ! if true _img contain imaginary part of complex number
  end type info

  type, private :: maxe
    logical :: complx
    real(8) :: val1     = 0.0d0
    real(8) :: val2     = 0.0d0
    real(8) :: val1_img = 0.0d0
    real(8) :: val2_img = 0.0d0
  end type maxe

contains

   subroutine init_info(type_info)
     type(info) :: type_info
     type_info%ival = 0
     type_info%val1 = 0.0d0
     type_info%val2 = 0.0d0
     type_info%reld = 0.0d0
     type_info%absd = 0.0d0
     type_info%val1_img = 0.0d0
     type_info%val2_img = 0.0d0
     type_info%complx   = .False.
   end subroutine

   function all_tests_passed()
     logical :: all_tests_passed
     all_tests_passed = test_status
     test_status = .True.  ! restart test_status
   end function

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine compare_real(value1, value2, thresh, check, comparing)

  character(len=*), intent(in), optional :: check, comparing
  real(8), intent(in) :: thresh, value1, value2

  logical :: diff_found
  real(8) :: diff_abs, diff_rel
  character(4) :: check_local

  ! initialize variables
  check_local = 'tota'
  if(present(check))then
    if(check(1:4) == 'abso') check_local = 'abso'
    if(check(1:4) == 'rela') check_local = 'rela'
    if(check(1:4) == 'tota') check_local = 'tota'
  endif

  diff_abs = dabs(value1 - value2)
  if( dabs(value1 + value2) > diff_threshold )then
    diff_rel = 2.0d0 * dabs(value1 - value2) / dabs(value1 + value2)
  else
    diff_rel = 0.0d0
  endif

  write(6,'(a)') ' ================================================================= '
  if(     present(comparing))  write(6,'(a,a)') ' Compare two real numbers: ', comparing
  if(.not.present(comparing))  write(6,'(a)')   ' Compare two real numbers '
  write(6,'(a)') '   a --> reference  number'
  write(6,'(a)') '   b --> calculated number'
  write(6,'(a,es10.2)') ' Threshold  = ',thresh
  write(6,'(a)')
  write(6,'(a)') ' ----------------------------------------------------------------- '
  write(6,'(a,d12.4)')  ' Absolute difference:     max|a-b|    = ',diff_abs
  write(6,'(a,d12.4)')  ' Relative difference: |2*(a-b)/(a+b)| = ',diff_rel
  write(6,'(a,d25.16)') ' a = ',value1
  write(6,'(a,d25.16)') ' b = ',value2
  write(6,'(a)')

  ! check result
  diff_found = .False.
  if(check_local == 'abso')then
    write(6,'(a)') ' Threshold is applied on the absolute difference.'
    if(diff_abs > thresh) diff_found = .True.
  endif
  if(check_local == 'rela')then
    write(6,'(a)') ' Threshold is applied on the relative difference.'
    if(diff_rel > thresh) diff_found = .True.
  endif
  if(check_local == 'tota')then
    if(diff_abs > thresh .and. diff_rel > thresh) diff_found = .True.
  endif

  if(diff_found)then
    test_status = .False.
    write(6,'(a)')
    write(6,'(a)') ' ERROR in compare_real routine: Examined numbers differ'
    write(0,'(a)') ' ERROR in compare_real routine: Examined numbers differ'
  endif

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine compare_complex(value1, value2, thresh, check, comparing)

  character(len=*), intent(in), optional :: check, comparing
  real(8),    intent(in) :: thresh
  complex(8), intent(in) :: value1, value2

  logical :: diff_found
  real(8) :: diff_abs, diff_rel, tmp
  character(len=4) :: check_local

  ! initialize variables
  check_local = 'tota'
  if(present(check))then
    if(check(1:4) == 'abso') check_local = 'abso'
    if(check(1:4) == 'rela') check_local = 'rela'
    if(check(1:4) == 'tota') check_local = 'tota'
  endif

  ! absolute difference
  diff_abs = (dreal(value1)-dreal(value2)) * (dreal(value1)-dreal(value2)) &
            +(dimag(value1)-dimag(value2)) * (dimag(value1)-dimag(value2))
  diff_abs = dsqrt(diff_abs)

  ! relative difference
  tmp = (dreal(value1)-dreal(value2)) * (dreal(value1)-dreal(value2)) &
      + (dimag(value1)-dimag(value2)) * (dimag(value1)-dimag(value2))
  diff_rel = dsqrt(tmp)
  tmp = (dreal(value1)+dreal(value2)) * (dreal(value1)+dreal(value2)) &
      + (dimag(value1)+dimag(value2)) * (dimag(value1)+dimag(value2))

  if( dabs(tmp) > diff_threshold )then
    diff_rel = 2.0d0 * diff_rel / dsqrt(tmp)
  else
    diff_rel = 0.0d0
  endif

  write(6,'(a)') ' ================================================================= '
  if(     present(comparing))  write(6,'(a,a)') ' Compare two complex numbers: ', comparing
  if(.not.present(comparing))  write(6,'(a)')   ' Compare two complex numbers '
  write(6,'(a)') '   a --> reference  number'
  write(6,'(a)') '   b --> calculated number'
  write(6,'(a,es10.2)') ' Threshold  = ',thresh
  write(6,'(a)')
  write(6,'(a)') ' ----------------------------------------------------------------- '
  write(6,'(a,d12.4)')  ' Absolute difference:     max|a-b|    = ',diff_abs
  write(6,'(a,d12.4)')  ' Relative difference: |2*(a-b)/(a+b)| = ',diff_rel
  write(6,'(a,2d25.16)') ' a = ',dreal(value1),dimag(value1)
  write(6,'(a,2d25.16)') ' b = ',dreal(value2),dimag(value2)
  write(6,'(a)')

  ! check result
  diff_found = .False.
  if(check_local == 'abso')then
    write(6,'(a)') ' Threshold is applied on the absolute difference.'
    if(diff_abs > thresh) diff_found = .True.
  endif
  if(check_local == 'rela')then
    write(6,'(a)') ' Threshold is applied on the relative difference.'
    if(diff_rel > thresh) diff_found = .True.
  endif
  if(check_local == 'tota')then
    if(diff_abs > thresh .and. diff_rel > thresh) diff_found = .True.
  endif

  if(diff_found)then
    test_status = .False.
    write(6,'(a)')
    write(6,'(a)') ' ERROR in compare_complex routine: Examined numbers differ'
    write(0,'(a)') ' ERROR in compare_complex routine: Examined numbers differ'
  endif

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine compare_vector(iunit, file_name, vector, thresh, check, comparing)
  ! compare two vectors
  ! one vector is read from file "file_name" the other is in "vector" array
  ! stop if there are two elements which maximum relative difference is bigger than threshold

  character(len=*), intent(in) :: file_name
  integer,          intent(in) :: iunit
  real(8),          intent(in) :: thresh, vector(*)
  character(len=*), intent(in), optional :: check, comparing

  integer :: ier, ic
  logical :: diff_found
  real(8) :: diff_abs, diff_rel, check_val, sdev, dnum
  type( info ) :: absolute, relative, total_rel, total_abs
  type( maxe ) :: maxelement
  character(4) :: check_local

  ! open file_name
  open(iunit,file=file_name,status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'compare_vector: ERROR while opening file ',file_name
    write(0,*) 'compare_vector: ERROR while opening file ',file_name
    stop
  endif
  rewind iunit

  ! initialize variables
  maxelement%complx = .False.
  check_local = 'tota'
  if(present(check))then
    if(check(1:4) == 'stan') check_local = 'stan'
    if(check(1:4) == 'abso') check_local = 'abso'
    if(check(1:4) == 'rela') check_local = 'rela'
    if(check(1:4) == 'tota') check_local = 'tota'
  endif
  call init_info(absolute)
  call init_info(relative)
  call init_info(total_rel)
  call init_info(total_abs)
  diff_found = .False.
  sdev = 0.0d0

  ! loop over all elements in the file
  do

    ic = ic + 1

    ! read data from file
    read(iunit,*,iostat=ier) dnum
    if(ier == -1)exit
    if(ier /=  0)then
      write(6,*) 'compare_vector: ERROR reading file ',file_name
      write(0,*) 'compare_vector: ERROR reading file ',file_name
      write(6,*) 'ier: ',ier
      stop
    endif

    ! absolute difference, relative difference and standard deviation
    diff_abs = dabs( dnum-vector(ic) )
    if( dabs(dnum+vector(ic)) > diff_threshold )then
      diff_rel = 2.0d0 * dabs( dnum-vector(ic) ) / dabs( dnum+vector(ic) )
    else
      diff_rel = 0.0d0
    endif
    sdev  = sdev + diff_abs*diff_abs

    ! find maximum values
    if( dabs(vector(ic)) > maxelement%val1 ) maxelement%val1 = vector(ic)
    if( dabs(dnum)       > maxelement%val2 ) maxelement%val2 = dnum
    ! store info about values with the maximum relative error
    if( diff_rel > relative%reld )then
      relative%ival = ic
      relative%val1 = dnum
      relative%val2 = vector(ic)
      relative%reld = diff_rel
      relative%absd = diff_abs
    endif
    ! store info about values with the maximum absolute error
    if( diff_abs > absolute%absd )then
      absolute%ival = ic
      absolute%val1 = dnum
      absolute%val2 = vector(ic)
      absolute%reld = diff_rel
      absolute%absd = diff_abs
    endif
    ! by default both relative and absolute errors must be under threshold for test to fail
    if( (diff_rel > thresh) .and. (diff_abs > thresh) )then
      if( diff_rel > total_rel%reld )then
       total_rel%ival = ic
       total_rel%val1 = dnum
       total_rel%val2 = vector(ic)
       total_rel%reld = diff_rel
       total_rel%absd = diff_abs
      endif
      if( diff_abs > total_abs%absd )then
       total_abs%ival = ic
       total_abs%val1 = dnum
       total_abs%val2 = vector(ic)
       total_abs%reld = diff_rel
       total_abs%absd = diff_abs
      endif
      diff_found = .True.
    endif

  enddo

  ! close file_name
  close(iunit,status='keep')

  ! print info
  ic = ic - 1
  sdev = dsqrt(sdev/dfloat(ic))
  call print_err_info(ic,thresh,absolute,relative,total_rel,total_abs,sdev, &
                      diff_found,comparing,maxelement)

  ! check result
  if(check_local /= 'tota')then
    if(check_local == 'stan')then
      write(6,'(a)') ' Threshold is applied on the standard deviation.'
      check_val = sdev
    endif
    if(check_local == 'abso')then
      write(6,'(a)') ' Threshold is applied on the absolute difference.'
      check_val = absolute%absd
    endif
    if(check_local == 'rela')then
      write(6,'(a)') ' Threshold is applied on the relative difference.'
      check_val = relative%reld
    endif
    diff_found = .False.
    if(check_val > thresh) diff_found = .True.
  endif
  if(diff_found)then
    test_status = .False.
    write(6,'(a)')
    write(6,'(a)') ' ERROR in compare_vector routine: Examined file and vector differ'
    write(0,'(a)') ' ERROR in compare_vector routine: Examined file and vector differ'
  endif

  end subroutine compare_vector

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine compare_1d_real_vectors(ndim, vector1, vector2, thresh, check, print_level, &
                                     comparing)
    integer, intent(in) :: ndim
    real(8), intent(in) :: thresh, vector1(:), vector2(:)
    character(len=*), intent(in), optional :: check, comparing
    integer,          intent(in), optional :: print_level
    call compare_real_vectors(ndim, vector1, vector2, thresh, check, print_level, comparing)
  end subroutine

  subroutine compare_2d_real_vectors(ndim, vector1, vector2, thresh, check, print_level, &
                                     comparing)
    integer, intent(in) :: ndim
    real(8), intent(in) :: thresh, vector1(:,:), vector2(:,:)
    character(len=*), intent(in), optional :: check, comparing
    integer,          intent(in), optional :: print_level
    call compare_real_vectors(ndim, vector1, vector2, thresh, check, print_level, comparing)
  end subroutine

  subroutine compare_3d_real_vectors(ndim, vector1, vector2, thresh, check, print_level, &
                                     comparing)
    integer, intent(in) :: ndim
    real(8), intent(in) :: thresh, vector1(:,:,:), vector2(:,:,:)
    character(len=*), intent(in), optional :: check, comparing
    integer,          intent(in), optional :: print_level
    call compare_real_vectors(ndim, vector1, vector2, thresh, check, print_level, comparing)
  end subroutine

  subroutine compare_4d_real_vectors(ndim, vector1, vector2, thresh, check, print_level, &
                                     comparing)
    integer, intent(in) :: ndim
    real(8), intent(in) :: thresh, vector1(:,:,:,:), vector2(:,:,:,:)
    character(len=*), intent(in), optional :: check, comparing
    integer,          intent(in), optional :: print_level
    call compare_real_vectors(ndim, vector1, vector2, thresh, check, print_level, comparing)
  end subroutine

  subroutine compare_5d_real_vectors(ndim, vector1, vector2, thresh, check, print_level, &
                                     comparing)
    integer, intent(in) :: ndim
    real(8), intent(in) :: thresh, vector1(:,:,:,:,:), vector2(:,:,:,:,:)
    character(len=*), intent(in), optional :: check, comparing
    integer,          intent(in), optional :: print_level
    call compare_real_vectors(ndim, vector1, vector2, thresh, check, print_level, comparing)
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine compare_real_vectors(ndim, vector1, vector2, thresh, check, print_level, &
                                  comparing)
  ! compare two real vectors
  ! data is stored in "vector1" and "vector2" arrays (ndim is number of elements to be compared)
  ! stop if there are two elements which maximum relative difference is bigger than threshold

  integer, intent(in) :: ndim
  real(8), intent(in) :: thresh, vector1(*), vector2(*)
  character(len=*), intent(in), optional :: check, comparing
  integer,          intent(in), optional :: print_level

  integer :: ic
  logical :: diff_found
  real(8) :: diff_abs, diff_rel, check_val, sdev
  type( info ) :: absolute, relative, total_rel, total_abs
  type( maxe ) :: maxelement
  character(4) :: check_local

  ! initialize variables
  check_local = 'tota'
  if(present(check))then
    if(check(1:4) == 'stan') check_local = 'stan'
    if(check(1:4) == 'abso') check_local = 'abso'
    if(check(1:4) == 'rela') check_local = 'rela'
    if(check(1:4) == 'tota') check_local = 'tota'
  endif
  call init_info(absolute)
  call init_info(relative)
  call init_info(total_rel)
  call init_info(total_abs)
  diff_found = .False.
  sdev = 0.0d0

  ! loop over elements
  do ic=1,ndim

    ! absolute difference, relative difference and standard deviation
    diff_abs = dabs( vector1(ic)-vector2(ic) )
    if( dabs(vector1(ic)+vector2(ic)) > diff_threshold )then
      diff_rel = 2.0d0 * dabs( vector1(ic)-vector2(ic) ) / dabs( vector1(ic)+vector2(ic) )
    else
      diff_rel = 0.0d0
    endif
    sdev  = sdev + diff_abs*diff_abs

    ! store info about values with the maximum relative error
    if( diff_rel > relative%reld )then
      relative%ival = ic
      relative%val1 = vector1(ic)
      relative%val2 = vector2(ic)
      relative%reld = diff_rel
      relative%absd = diff_abs
    endif
    ! store info about values with the maximum absolute error
    if( diff_abs > absolute%absd )then
      absolute%ival = ic
      absolute%val1 = vector1(ic)
      absolute%val2 = vector2(ic)
      absolute%reld = diff_rel
      absolute%absd = diff_abs
    endif
    ! by default both relative and absolute errors must be under threshold for test to fail
    if( (diff_rel > thresh) .and. (diff_abs > thresh) )then
      if( diff_rel > total_rel%reld )then
       total_rel%ival = ic
       total_rel%val1 = vector1(ic)
       total_rel%val2 = vector2(ic)
       total_rel%reld = diff_rel
       total_rel%absd = diff_abs
      endif
      if( diff_abs > total_abs%absd )then
       total_abs%ival = ic
       total_abs%val1 = vector1(ic)
       total_abs%val2 = vector2(ic)
       total_abs%reld = diff_rel
       total_abs%absd = diff_abs
      endif
      diff_found = .True.
    endif

  enddo

  if(present(print_level))then
    write(6,'(/a)')' --------------------------------'
    write(6,'(/a)')'       reference set      |     calculated set      |  abs_diff  |  rel_diff'
    do ic=1,ndim
      diff_abs = dabs( vector1(ic)-vector2(ic) )
      diff_rel = 2.0d0 * dabs( vector1(ic)-vector2(ic) ) / dabs( vector1(ic)+vector2(ic) )
      if( diff_abs > absolute%absd/10**print_level) &
        write(6,'(d25.16,1x,d25.16,2d13.4)') vector1(ic), vector2(ic), diff_abs, diff_rel
    enddo
  endif

  ! if there are no differences get maximum values
  if( (absolute%ival==0) .and. (relative%ival==0) )then
    maxelement%complx = .False.
    maxelement%val1   = maxval(vector1(1:ndim))
    maxelement%val2   = maxval(vector2(1:ndim))
  endif

  ! print info
  sdev = dsqrt(sdev/dfloat(ndim))
  call print_err_info(ndim,thresh,absolute,relative,total_rel,total_abs,sdev, &
                      diff_found,comparing,maxelement)

  ! check result
  if(check_local /= 'tota')then
    if(check_local == 'stan')then
      write(6,'(a)') ' Threshold is applied on the standard deviation.'
      check_val = sdev
    endif
    if(check_local == 'abso')then
      write(6,'(a)') ' Threshold is applied on the absolute difference.'
      check_val = absolute%absd
    endif
    if(check_local == 'rela')then
      write(6,'(a)') ' Threshold is applied on the relative difference.'
      check_val = relative%reld
    endif
    diff_found = .False.
    if(check_val > thresh) diff_found = .True.
  endif
  if(diff_found)then
    test_status = .False.
    write(6,'(a)')
    write(6,'(a)') ' ERROR in compare_real_vectors routine: Examined arrays differ too much'
    write(0,'(a)') ' ERROR in compare_real_vectors routine: Examined arrays differ too much'
  endif

  end subroutine compare_real_vectors

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine compare_1d_complex_vectors(ndim, vector1, vector2, thresh, check, comparing)
    integer,    intent(in) :: ndim
    real(8),    intent(in) :: thresh
    complex(8), intent(in) :: vector1(:), vector2(:)
    character(len=*), intent(in), optional :: check, comparing
    call compare_complex_vectors(ndim, vector1, vector2, thresh, check, comparing)
  end subroutine

  subroutine compare_2d_complex_vectors(ndim, vector1, vector2, thresh, check, comparing)
    integer,    intent(in) :: ndim
    real(8),    intent(in) :: thresh
    complex(8), intent(in) :: vector1(:,:), vector2(:,:)
    character(len=*), intent(in), optional :: check, comparing
    call compare_complex_vectors(ndim, vector1, vector2, thresh, check, comparing)
  end subroutine

  subroutine compare_3d_complex_vectors(ndim, vector1, vector2, thresh, check, comparing)
    integer,    intent(in) :: ndim
    real(8),    intent(in) :: thresh
    complex(8), intent(in) :: vector1(:,:,:), vector2(:,:,:)
    character(len=*), intent(in), optional :: check, comparing
    call compare_complex_vectors(ndim, vector1, vector2, thresh, check, comparing)
  end subroutine

  subroutine compare_4d_complex_vectors(ndim, vector1, vector2, thresh, check, comparing)
    integer,    intent(in) :: ndim
    real(8),    intent(in) :: thresh
    complex(8), intent(in) :: vector1(:,:,:,:), vector2(:,:,:,:)
    character(len=*), intent(in), optional :: check, comparing
    call compare_complex_vectors(ndim, vector1, vector2, thresh, check, comparing)
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine compare_complex_vectors(ndim, vector1, vector2, thresh, check, comparing)
  ! compare two complex vectors
  ! data is stored in "vector1" and "vector2" arrays (ndim is number of elements to be compared)
  ! stop if there are two elements which maximum relative difference is bigger than threshold

  integer,    intent(in) :: ndim
  real(8),    intent(in) :: thresh
  complex(8), intent(in) :: vector1(*), vector2(*)
  character(len=*), intent(in), optional :: check, comparing

  integer :: ic
  logical :: diff_found
  real(8) :: sdev, diff_abs, diff_rel, tmp, check_val
  type( info ) :: absolute, relative, total_rel, total_abs
  type( maxe ) :: maxelement
  character(4) :: check_local

  ! initialize variables
  check_local = 'tota'
  if(present(check))then
    if(check(1:4) == 'stan') check_local = 'stan'
    if(check(1:4) == 'abso') check_local = 'abso'
    if(check(1:4) == 'rela') check_local = 'rela'
    if(check(1:4) == 'tota') check_local = 'tota'
  endif
  call init_info(absolute)
  call init_info(relative)
  call init_info(total_rel)
  call init_info(total_abs)
  absolute%complx  = .True.
  relative%complx  = .True.
  total_rel%complx = .True.
  total_abs%complx = .True.
  diff_found = .False.
  sdev = 0.0d0

  ! loop over elements
  do ic=1,ndim

    ! absolute difference
    diff_abs=(dreal(vector1(ic))-dreal(vector2(ic))) * (dreal(vector1(ic))-dreal(vector2(ic))) &
            +(dimag(vector1(ic))-dimag(vector2(ic))) * (dimag(vector1(ic))-dimag(vector2(ic)))
    diff_abs=dsqrt(diff_abs)
    ! standard deviation
    sdev = sdev + diff_abs*diff_abs

    ! relative difference
    tmp = (dreal(vector1(ic))-dreal(vector2(ic))) * (dreal(vector1(ic))-dreal(vector2(ic))) &
        + (dimag(vector1(ic))-dimag(vector2(ic))) * (dimag(vector1(ic))-dimag(vector2(ic)))
    diff_rel = dsqrt(tmp)
    tmp = (dreal(vector1(ic))+dreal(vector2(ic))) * (dreal(vector1(ic))+dreal(vector2(ic))) &
        + (dimag(vector1(ic))+dimag(vector2(ic))) * (dimag(vector1(ic))+dimag(vector2(ic)))
    if( dabs(tmp) > diff_threshold )then
      diff_rel = 2.0d0 * diff_rel / dsqrt(tmp)
    else
      diff_rel = 0.0d0
    endif

    ! store info about values with the maximum relative error
    if( diff_rel > relative%reld )then
      relative%ival = ic
      relative%val1 = dreal(vector1(ic))
      relative%val2 = dreal(vector2(ic))
      relative%val1_img = dimag(vector1(ic))
      relative%val2_img = dimag(vector2(ic))
      relative%reld = diff_rel
      relative%absd = diff_abs
    endif
    ! store info about values with the maximum absolute error
    if( diff_abs > absolute%absd )then
      absolute%ival = ic
      absolute%val1 = dreal(vector1(ic))
      absolute%val2 = dreal(vector2(ic))
      absolute%val1_img = dimag(vector1(ic))
      absolute%val2_img = dimag(vector2(ic))
      absolute%reld = diff_rel
      absolute%absd = diff_abs
    endif
    ! by default both relative and absolute errors must be under threshold for test to fail
    if( (diff_rel > thresh) .and. (diff_abs > thresh) )then
      if( diff_rel > total_rel%reld )then
       total_rel%ival = ic
       total_rel%val1 = dreal(vector1(ic))
       total_rel%val2 = dreal(vector2(ic))
       total_rel%val1_img = dimag(vector1(ic))
       total_rel%val2_img = dimag(vector2(ic))
       total_rel%reld = diff_rel
       total_rel%absd = diff_abs
      endif
      if( diff_abs > total_abs%absd )then
       total_abs%ival = ic
       total_abs%val1 = dreal(vector1(ic))
       total_abs%val2 = dreal(vector2(ic))
       total_abs%val1_img = dimag(vector1(ic))
       total_abs%val2_img = dimag(vector2(ic))
       total_abs%reld = diff_rel
       total_abs%absd = diff_abs
      endif
      diff_found = .True.
    endif

  enddo

  ! if there are no differences get maximum values
  if( (absolute%ival==0) .and. (relative%ival==0) )then
    maxelement%complx = .True.
    maxelement%val1     = maxval(dreal(vector1(1:ndim)))
    maxelement%val2     = maxval(dreal(vector2(1:ndim)))
    maxelement%val1_img = maxval(dimag(vector1(1:ndim)))
    maxelement%val2_img = maxval(dimag(vector2(1:ndim)))
  endif

  ! print info
  sdev = dsqrt(sdev/dfloat(ndim))
  call print_err_info(ndim,thresh,absolute,relative,total_rel,total_abs,sdev, &
                      diff_found,comparing,maxelement)

  ! check result
  if(check_local /= 'tota')then
    if(check_local == 'stan')then
      write(6,'(a)') ' Threshold is applied on the standard deviation.'
      check_val = sdev
    endif
    if(check_local == 'abso')then
      write(6,'(a)') ' Threshold is applied on the absolute difference.'
      check_val = absolute%absd
    endif
    if(check_local == 'rela')then
      write(6,'(a)') ' Threshold is applied on the relative difference.'
      check_val = relative%reld
    endif
    diff_found = .False.
    if(check_val > thresh) diff_found = .True.
  endif
  if(diff_found)then
    test_status = .False.
    write(6,'(a)')
    write(6,'(a)')' ERROR in compare_complex_vectors routine: Examined arrays differ too much'
    write(0,'(a)')' ERROR in compare_complex_vectors routine: Examined arrays differ too much'
  endif

  end subroutine compare_complex_vectors

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine compare_files(iunit, file1, file2, thresh, check, comparing)
  ! compare real data
  ! data is stored in files "file1" and "file2"
  ! stop if there are two elements which maximum relative difference is bigger than threshold

  character(len=*), intent(in) :: file1, file2
  integer,          intent(in) :: iunit
  real(8),          intent(in) :: thresh
  character(len=*), intent(in), optional :: check, comparing

  integer :: iu1, iu2, ier, ier1, ier2, ic
  logical :: diff_found
  real(8) :: dnum1, dnum2, sdev, diff_abs, diff_rel, check_val
  type( info ) :: absolute, relative, total_rel, total_abs
  type( maxe ) :: maxelement
  character(4) :: check_local

  iu1 = iunit
  iu2 = iunit + 1

  ! open file1
  open(iu1,file=file1,status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'compare_files: ERROR while opening file ',file1
    write(0,*) 'compare_files: ERROR while opening file ',file1
    stop
  endif
  rewind iu1

  ! open file2
  open(iu2,file=file2,status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'compare_files: ERROR while opening file ',file2
    write(0,*) 'compare_files: ERROR while opening file ',file2
    stop
  endif
  rewind iu2

  ! initialize variables
  maxelement%complx = .False.
  check_local = 'tota'
  if(present(check))then
    if(check(1:4) == 'stan') check_local = 'stan'
    if(check(1:4) == 'abso') check_local = 'abso'
    if(check(1:4) == 'rela') check_local = 'rela'
    if(check(1:4) == 'tota') check_local = 'tota'
  endif
  call init_info(absolute)
  call init_info(relative)
  call init_info(total_rel)
  call init_info(total_abs)
  diff_found = .False.
  sdev = 0.0d0
  ic         = 0

  ! loop over all elements in files
  do

    ic = ic + 1

    ! read data from both files
    read(iu1,*,iostat=ier1) dnum1
    read(iu2,*,iostat=ier2) dnum2
    if( (ier1 == -1) .and. (ier2 == -1) ) exit
    if( (ier1 == -1) .or.  (ier2 == -1) )then
      write(6,*) ' ERROR in compare_files: Examined files have different length '
      write(0,*) ' ERROR in compare_files: Examined files have different length '
      stop
    endif
    if( (ier1 /=  0) .or. (ier2 /=  0) )then
      write(6,*) ' ERROR in compare_files: reading one of the files ',file1,' ',file2
      write(0,*) ' ERROR in compare_files: reading one of the files ',file1,' ',file2
      write(6,*) ' ier1, ier2: ',ier1, ier2
      stop
    endif

    ! absolute difference, relative difference and standard deviation
    diff_abs = dabs( dnum1-dnum2 )
    if( dabs(dnum1+dnum2) > diff_threshold )then
      diff_rel = 2.0d0 * dabs( dnum1-dnum2 ) / dabs( dnum1+dnum2 )
    else
      diff_rel = 0.0d0
    endif
    sdev = sdev + diff_abs*diff_abs

    ! find maximum values
    if( dabs(dnum1) > maxelement%val1 ) maxelement%val1 = dnum1
    if( dabs(dnum2) > maxelement%val2 ) maxelement%val2 = dnum2
    ! store info about values with the maximum relative error
    if( diff_rel > relative%reld )then
      relative%ival = ic
      relative%val1 = dnum2  ! reference file
      relative%val2 = dnum1
      relative%reld = diff_rel
      relative%absd = diff_abs
    endif
    ! store info about values with the maximum absolute error
    if( diff_abs > absolute%absd )then
      absolute%ival = ic
      absolute%val1 = dnum2  ! reference file
      absolute%val2 = dnum1
      absolute%reld = diff_rel
      absolute%absd = diff_abs
    endif
    ! by default both relative and absolute errors must be under threshold for test to fail
    if( (diff_rel > thresh) .and. (diff_abs > thresh) )then
      if( diff_rel > total_rel%reld )then
       total_rel%ival = ic
       total_rel%val1 = dnum2  ! reference file
       total_rel%val2 = dnum1
       total_rel%reld = diff_rel
       total_rel%absd = diff_abs
      endif
      if( diff_abs > total_abs%absd )then
       total_abs%ival = ic
       total_abs%val1 = dnum2  ! reference file
       total_abs%val2 = dnum1
       total_abs%reld = diff_rel
       total_abs%absd = diff_abs
      endif
      diff_found = .True.
    endif

  enddo

  ! close file1
  close(iu1,status='keep')

  ! close file2
  close(iu2,status='keep')

  ! print info
  ic = ic - 1
  sdev = dsqrt(sdev/dfloat(ic))
  call print_err_info(ic,thresh,absolute,relative,total_rel,total_abs,sdev, &
                      diff_found,comparing,maxelement)

  ! check result
  if(check_local /= 'tota')then
    if(check_local == 'stan')then
      write(6,'(a)') ' Threshold is applied on the standard deviation.'
      check_val = sdev
    endif
    if(check_local == 'abso')then
      write(6,'(a)') ' Threshold is applied on the absolute difference.'
      check_val = absolute%absd
    endif
    if(check_local == 'rela')then
      write(6,'(a)') ' Threshold is applied on the relative difference.'
      check_val = relative%reld
    endif
    diff_found = .False.
    if(check_val > thresh) diff_found = .True.
  endif
  if(diff_found)then
    test_status = .False.
    write(6,'(a)')
    write(6,'(a)') ' ERROR in compare_files routine: Examined files differ'
    write(0,'(a)') ' ERROR in compare_files routine: Examined files differ'
  endif

  end subroutine compare_files

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine compare_text_files(iunit, file1, file2)
  ! compare two text files; they must be identical

  character(*), intent(in) :: file1, file2
  integer,      intent(in) :: iunit

  character(250) :: line1, line2
  integer :: ic, iu1, iu2, ier, ier1, ier2
  logical :: all_is_good

  iu1 = iunit
  iu2 = iunit + 1
  all_is_good = .True.

  ! open file1
  open(iu1,file=file1,status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'compare_text_files: ERROR while opening file ',file1
    write(0,*) 'compare_text_files: ERROR while opening file ',file1
    stop
  endif
  rewind iu1

  ! open file2
  open(iu2,file=file2,status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'compare_text_files: ERROR while opening file ',file2
    write(0,*) 'compare_text_files: ERROR while opening file ',file2
    stop
  endif
  rewind iu2

  ! compare two files
  ic = 0
  do

    ic = ic + 1
    read(iu1,'(a)',iostat=ier1) line1
    read(iu2,'(a)',iostat=ier2) line2
    if( (ier1 == -1) .and. (ier2 == -1) ) exit
    if( (ier1 == -1) .or.  (ier2 == -1) )then
      write(6,*) ' ERROR: Examined files have different length '
      write(0,*) ' ERROR: Examined files have different length '
      stop
    endif
    if( (ier1 /=  0) .or. (ier2 /=  0) )then
      write(6,*) 'compare_text_files: ERROR reading one of the files ',file1,' ',file2
      write(0,*) 'compare_text_files: ERROR reading one of the files ',file1,' ',file2
      write(6,*) 'ier1, ier2: ',ier1, ier2
      stop
    endif

    if( line1 /= line2 .and. all_is_good )then
      all_is_good = .False.
      test_status = .False.
      write(6,*) ' ERROR in compare_text_files routine'
      write(6,*) '       files (',file1,' and ',file2,') differ on line ',ic
      write(6,'(a)') trim(line1)
      write(6,'(a)') trim(line2)
      write(0,*) ' ERROR: files differ (see output for more details) '
    endif

  enddo

  ! close file1
  close(iu1,status='keep')

  ! close file2
  close(iu2,status='keep')

  end subroutine compare_text_files

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine print_err_info(ndim,thresh,absolute,relative,total_rel,total_abs,sdev, &
                            diff_found,comparing,maxelement)

  character(len=*), intent(in), optional :: comparing
  logical, intent(in) :: diff_found
  integer, intent(in) :: ndim
  real(8), intent(in) :: thresh, sdev
  type( info ), intent(in) :: absolute, relative, total_rel, total_abs
  type( maxe ) :: maxelement

  write(6,'(a)') ' ================================================================= '
  if(     present(comparing)) write(6,'(a,a)') ' Compare two sets of numbers: ', comparing
  if(.not.present(comparing)) write(6,'(a)')   ' Compare two sets of numbers '
  write(6,'(a)') '   a --> reference  set'
  write(6,'(a)') '   b --> calculated set'
  write(6,'(a,es10.2)') ' Threshold  = ',thresh
  write(6,'(a,i0)')     ' Number of compared elements  = ',ndim
  write(6,'(a)')
  write(6,'(a)') ' ----------------------------------------------------------------- '

  if( (absolute%ival==0) .and. (relative%ival==0) )then

    write(6,'(a)') '                 There were found no differences.'
    write(6,'(a)')
    write(6,'(a)') ' Maximum element of given array: '
    if(maxelement%complx)then
      write(6,'(a,2d25.16)' ) ' a = ',maxelement%val1,maxelement%val1_img
      write(6,'(a,2d25.16)' ) ' b = ',maxelement%val2,maxelement%val2_img
    else
      write(6,'(a,d25.16)' ) ' a = ',maxelement%val1
      write(6,'(a,d25.16)' ) ' b = ',maxelement%val2
    endif
    write(6,'(a)')

  else

    write(6,'(a,d12.4)')  ' The largest absolute difference: max|a-b|  = ',absolute%absd
    write(6,'(a,d12.4)')  ' Their relative difference: |2*(a-b)/(a+b)| = ',absolute%reld
    if(absolute%complx)then
      write(6,'(a,2d25.16)') ' a = ',absolute%val1,absolute%val1_img
      write(6,'(a,2d25.16)') ' b = ',absolute%val2,absolute%val2_img
    else
      write(6,'(a,d25.16)') ' a = ',absolute%val1
      write(6,'(a,d25.16)') ' b = ',absolute%val2
    endif
    write(6,'(a,i0)')     ' Position of these elements = ',absolute%ival
    write(6,'(a,d12.4)')  ' Standard deviation: sqrt( (a-b)^2 / n ) = ',sdev
    write(6,'(a)')
    write(6,'(a)') ' ----------------------------------------------------------------- '
    write(6,'(a,d12.4)')  ' The largest relative difference: |2*max|a-b|/(a+b)| = ',relative%reld
    write(6,'(a,d12.4)')  ' Their absolute difference:       max|a-b|           = ',relative%absd
    if(relative%complx)then
      write(6,'(a,2d25.16)') ' a = ',relative%val1,relative%val1_img
      write(6,'(a,2d25.16)') ' b = ',relative%val2,relative%val2_img
    else
      write(6,'(a,d25.16)') ' a = ',relative%val1
      write(6,'(a,d25.16)') ' b = ',relative%val2
    endif
    write(6,'(a,i0)')     ' Position of these elements = ',relative%ival
    write(6,'(a)')
    if(diff_found)then
      write(6,'(a)') ' ----------------------------------------------------------------- '
      write(6,'(a)') ' In analized sets are present numbers which have both '
      write(6,'(a)') '   relative and absolute difference bigger than given threshold.'
      write(6,'(a)')
      write(6,'(a)') ' Number with the bigest absolute error:'
      write(6,'(a,d12.4)')  ' Relative difference: |2*max|a-b|/(a+b)| = ',total_abs%reld
      write(6,'(a,d12.4)')  ' Absolute difference:       max|a-b|     = ',total_abs%absd
      if(total_abs%complx)then
        write(6,'(a,2d25.16)') ' a = ',total_abs%val1,total_abs%val1_img
        write(6,'(a,2d25.16)') ' b = ',total_abs%val2,total_abs%val2_img
      else
        write(6,'(a,d25.16)') ' a = ',total_abs%val1
        write(6,'(a,d25.16)') ' b = ',total_abs%val2
      endif
      write(6,'(a,i0)')     ' Position of these elements = ',total_abs%ival
      write(6,'(a)')
      write(6,'(a)') ' Number with the bigest relative error:'
      write(6,'(a,d12.4)')  ' Relative difference: |2*max|a-b|/(a+b)| = ',total_rel%reld
      write(6,'(a,d12.4)')  ' Absolute difference:       max|a-b|     = ',total_rel%absd
      if(total_abs%complx)then
        write(6,'(a,2d25.16)') ' a = ',total_rel%val1,total_rel%val1_img
        write(6,'(a,2d25.16)') ' b = ',total_rel%val2,total_rel%val2_img
      else
        write(6,'(a,d25.16)') ' a = ',total_rel%val1
        write(6,'(a,d25.16)') ' b = ',total_rel%val2
      endif
      write(6,'(a,i0)')     ' Position of these elements = ',total_rel%ival
      write(6,'(a)')
    else
      write(6,'(a)') ' There are no numbers which have both relative and '
      write(6,'(a)') '   absolute difference over given threshold. '
      write(6,'(a)')
    endif

  endif

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine unpackme(n,vec,mat)

  integer, intent( in) :: n
  real(8), intent( in) :: vec(:,:)
  real(8), intent(out) :: mat(:,:)

  integer :: i, j, ij, m1, m2, v1, v2

  m1 = size(mat,1)
  m2 = size(mat,2)
  v1 = size(vec,1)
  v2 = size(vec,2)

  if( m1*m2 /=  v1 )then
    write(6,'(/a)'   ) ' ERROR: Wrong dimensions in unpackme.'
    write(6,'(a,i0)' ) '        m1 = ',m1
    write(6,'(a,i0)' ) '        m2 = ',m2
    write(6,'(a,i0/)') '        v1 = ',v1
    stop
  endif

  if( v2 < n )then
    write(6,'(/a)'    ) ' ERROR: Vector array do not have enough components.'
    write(6,'(a,i0)'  ) '        n  = ',n
    write(6,'(a,i0/)' ) '        v2 = ',v2
    stop
  endif

  ij = 0
  do j=1,m2
    do i=1,m1
      ij = ij + 1
      mat(i,j) = vec(ij,n)
    enddo
  enddo

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

END MODULE UNIT_TESTING
