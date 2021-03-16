!================================================================================!
!                                                                                !
!  Program:      ReSpect, Dirac                                                  !
!                                                                                !
!  Module:       STRING_MANIPULATIONS                                            !
!                                                                                !
!  Description:  This module contains operations on derived data type "string".  !
!                It is used in unit_test_generator and file_reader module.       !
!                                                                                !
!  Dependencies: None                                                            !
!                                                                                !
!  Author:       Stanislav Komorovsky                                            !
!                                                                                !
!  Date:         Jan. - July 2012 (Toulouse, Tromso)                             !
!                                                                                !
!================================================================================!

MODULE STRING_MANIPULATIONS

implicit none

type, public :: string
  character, pointer :: str(:)
end type string

public string__set_max_length   ! set maximum alowed number of characters per string
public string__set_iunit        ! set unit for output
public string__count_words      ! integer function: count number of words in one line
public string__compare_strings  ! logical function: compare two strings
public string__associate_status ! get associate status of type string
public string__gets             ! string function: converts type string to the character string
public string__i2c              ! string function: converts integer variable to the string
public string__l2c              ! string function: converts logical variable to the string
public string__get_size         ! integer function: returns length of %str in type string
public string__get_char         ! character function: returns character %str(i) in type string
public string__set_char         ! set ith character of %str array in type string
public string__print            ! print %str from type string
public string__split            ! split line (type string) to words; return word(s)
public string__append           ! add string to the end of list of strings
public string__free             ! free (deallocate) pointer or variable type string
public string__init             ! initialize (nullify) pointer or variable type string

interface string__compare_strings
  module procedure compare_strings_type_type
  module procedure compare_strings_string_type
end interface

private

private :: errall, csize, add_char

integer, private, save :: max_string_length = 70, iunit = 6

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine string__set_max_length(n)
    integer :: n
    max_string_length = n
  end subroutine

  subroutine string__set_iunit(n)
    integer :: n
    iunit = n
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  function string__count_words(line,separator) result(nwords)
  ! split line into words, where "separator" contains character which divides words

  type(string), intent(in) :: line
  character, intent(in), optional :: separator

  character :: sep
  integer   :: i, nwords

  ! initialize
  nwords = 0

  ! set separator
  sep = ' '
  if(present(separator)) sep = separator

  ! split line
  i = 0
  do
    i = i + 1
    if(i > string__get_size(line))exit
    if(string__get_char(line,i) /= sep)then
      do
        i = i + 1
        if(i > string__get_size(line))exit
        if(string__get_char(line,i) == sep)exit
      enddo
      nwords = nwords + 1
    endif
  enddo

  end function

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  function compare_strings_type_type(string1,string2) result(match)
    ! compare two type string variables and return logical value of the comparison
    type(string), intent(in) :: string1, string2
    logical :: match
    integer :: i
    match = .True.
    do i=1,min(size(string1%str),size(string2%str))
      if(string1%str(i) /= string2%str(i)) match = .False.
    enddo
    if(size(string1%str) /= size(string2%str)) match = .False.
  end function

  function compare_strings_string_type(string1,string2) result(match)
    ! compare character string with type string variables
    ! return logical value of the comparison
    character(len=*), intent(in) :: string1
    type( string ),   intent(in) :: string2
    integer :: i
    logical :: match
    match = .True.
    do i=1,min(len(string1),size(string2%str))
      if(string1(i:i) /= string2%str(i)) match = .False.
    enddo
    if(len(string1) /= size(string2%str)) match = .False.
  end function

  pure function string__associate_status(type_string) result(l)
    type(string), intent(in) :: type_string
    logical :: l
    l = associated(type_string%str)
  end function

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  function string__gets(type_string) result(c)
    type(string), intent(in) :: type_string
    character(len=max_string_length) :: c
    integer :: i
    if( max_string_length < string__get_size(type_string) )then
      write(iunit,*) ' ERROR: Maximum allowed length of string exceeded in ', &
                             'string_manipulation module (string__gets function)'
      write(iunit,*) '        max_string_length = ',max_string_length
      write(iunit,*) '        string length     = ',string__get_size(type_string)
      write(0,*) ' ERROR in string_manipulation module (string__gets function)'
      stop
    endif
    do i=1,max_string_length
      c(i:i) = ' '
    enddo
    do i=1,string__get_size(type_string)
      c(i:i)=type_string%str(i)
    enddo
  end function

  function string__i2c(i) result(c)
    integer, intent(in) :: i
    character(len=max_string_length) :: c
    integer :: j
    if( max_string_length < csize(i) )then
      write(iunit,*) ' ERROR: Maximum allowed length of string exceeded in ', &
                             'string_manipulation module (string__i2c function)'
      write(iunit,*) '        max_string_length = ',max_string_length
      write(iunit,*) '        string length     = ',csize(i)
      write(0,*) ' ERROR in string_manipulation module (string__i2c function)'
      stop
    endif
    do j=1,max_string_length
      c(j:j) = ' '
    enddo
    write(c,'(i0)')i
  end function

  function string__l2c(l) result(c)
    logical, intent(in) :: l
    character(len=1) :: c
    if(l) c='T'
    if(.not.l) c='F'
  end function

  pure function csize(i) result(is)
    integer, intent(in) :: i
    integer :: is
    if(i == 0) is = 1
    if(i >  0) is = floor(log10(real( i,kind(1.0d0)))) + 1
    if(i <  0) is = floor(log10(real(-i,kind(1.0d0)))) + 2
  end function

  pure function string__get_size(type_string) result(i)
    type(string), intent(in) :: type_string
    integer :: i
    i = size(type_string%str)
  end function

  pure function string__get_char(type_string,i) result(c)
    type(string), intent(in) :: type_string
    integer, intent(in) :: i
    character :: c
    c = type_string%str(i)
  end function

  subroutine string__set_char(type_string,i,c)
    type(string), intent(inout) :: type_string
    character, intent(in) :: c
    integer, intent(in) :: i
    type_string%str(i) = c
  end subroutine

  subroutine string__print(char_string)
    type(string), intent(in) :: char_string
    integer :: i
    character(len=20) :: for
    for = ' '
    write(for,'(a,i0,a)') '(',string__get_size(char_string),'a)'
    write(iunit,for) (char_string%str(i),i=1,string__get_size(char_string))
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_char(word,chvar)
  ! add character to the end of word%str array
  ! n   --> dimension of input array
  ! n+1 --> dimension of output array

  character, intent(in) :: chvar
  character, intent(inout), pointer :: word(:)

  character, pointer :: tmp_str(:)
  integer :: i, n, ier

  ! initialize
  nullify(tmp_str)

  if(associated(word))then
    n = size(word)   ! store original size of word
    tmp_str => word  ! store temporarily word
  else
    n = 0
  endif

  ! reallocate word with greater dimension
  nullify( word )
  allocate( word(n+1), stat = ier ); call errall(ier,'add_char','word')

  ! restore old values and add new one
  do i=1,n
    word(i) = tmp_str(i)
  enddo
  word(n+1) = chvar

  ! nullify temporary pointer
  if(associated(tmp_str))then
    deallocate(tmp_str)
    nullify(tmp_str)
  endif

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine string__split(line,word,word_number,word_found)
  ! extract n-th word from line (n=word_number) and store it in the word variable
  ! if word_found is pressent use it to return if there exists n-th (word_number) word

  logical, intent(out), optional :: word_found
  integer, intent(in) :: word_number
  type( string ), intent(in   ) :: line
  type( string ), intent(inout) :: word

  logical :: reading_word
  integer :: ier, i, icount, nwfirst, nwlast, nwsize

  ! initialize
  nwfirst = 0
  nwlast  = 0
  icount  = 0
  reading_word = .False.
  if(present(word_found)) word_found = .True.
  if(associated(word%str)) call string__free(word)

  ! get size of the word
  do i=1,size(line%str)
    if( (line%str(i) /= ' ') .and. (reading_word .eqv. .False.) )then
      icount = icount + 1
      reading_word = .True.
      if(icount == word_number) nwfirst = i
    endif
    if( (line%str(i) == ' ') .and. (reading_word .eqv. .True.)  )then
      reading_word = .False.
      if(icount == word_number) nwlast = i
    endif
    if( nwfirst /= 0 .and. nwlast == 0 .and. i == size(line%str) ) nwlast = i + 1
  enddo
  nwsize = nwlast - nwfirst

  ! stop if n-th word does not exist
  if( nwlast == 0 .and. nwfirst == 0 )then
    if(present(word_found))then
      word_found = .False.
    else
      write(iunit,'(a)'   )
      write(iunit,'(a,i0)')'  ERROR while trying to read non-existing word '
      write(iunit,'(a,$)' )'  line     -->  '; call string__print(line)
      write(iunit,'(a,i0)')'  position -->  ',word_number
      write(iunit,'(a)'   )'  There is no word on that position!'
      write(iunit,'(a)'   )
      stop
    endif
  endif

  ! allocate word character array
  nullify( word%str )
  allocate( word%str( nwsize ), stat = ier ); call errall(ier,'string__split','word%str')

  ! copy string to newly allocated str array
  do i=1,nwsize
    word%str(i) = line%str(nwfirst+i-1)
  enddo

  end subroutine string__split

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine string__append(string_name,type_string)
  ! copy "string_name" to the "type_string" array

  character(len=*),   intent(in )   :: string_name
  type( string ), intent(inout) :: type_string

  integer :: i, ier, nstr

  if(associated(type_string%str)) deallocate(type_string%str)

  ! allocate str character array
  nstr = len(string_name)
  nullify( type_string%str )
  allocate( type_string%str(nstr), stat = ier ); call errall(ier,'string__append','type_string%str')

  ! copy string to newly allocated str array
  do i=1,nstr
    type_string%str(i) = string_name(i:i)
  enddo

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine string__free(string_value)
    type( string ), intent(inout) :: string_value
    if(associated(string_value%str))then
      deallocate( string_value%str )
      nullify( string_value%str )
    endif
  end subroutine

  subroutine string__init(string_value)
    type( string ), intent(out) :: string_value
    nullify(string_value%str)
  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine errall(ier,chsub,charr)
    integer      :: ier
    character(len=*) :: chsub, charr
    if(ier /= 0)then
      write(iunit,*)
      write(iunit,'(a,a,a,a)')' Error in allocation of ',charr,' in subroutine ',chsub
      write(iunit,'(a,i4)'   )' iostat =',ier
      stop 'ERROR in allocation (see output for more information)'
    endif
  end subroutine errall

END MODULE STRING_MANIPULATIONS
