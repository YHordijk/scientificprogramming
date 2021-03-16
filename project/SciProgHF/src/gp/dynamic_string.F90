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

module mystring
!--------------------------------------------------------------------
!
!   Module for simple dynamic string.
!
!  Example:  
!             use mystring
!             type(string) :: s,s1
!             call make_string(s,"this is my string")
!             call copy_string(s,s1)
!             call print_string(s1) or print get_string(s1)
!             call cancel_string(s); call cancel_string(s1) 
! 
!--------------------------------------------------------------------
 type string
  private
   character*1, dimension (:), pointer :: strp  => NULL()
   logical :: is_alloc = .false. ! this default does not work for older F90 compilers 
   integer :: str_length = 0 ! length of the allocated string
 end type string
 integer, private :: iprint_string = 0  ! control print-out variable
 contains 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine make_string(self, text)
!-----------------------------------------------------------------------
!       create (allocate) the string
!-----------------------------------------------------------------------
implicit none
#include "priunit.h"
 type(string), dimension(:), intent(inout) :: self
 !type(string), dimension(1), intent(inout) :: self
 character(*), intent(in) :: text
 integer :: i, alloc_stat

 if (.not.self(1)%is_alloc) then
   allocate( self(1)%strp(len(text)), stat = alloc_stat )
   ! ... store the actual size
   self(1)%str_length = size(self(1)%strp)
   if (self(1)%str_length.ne.len(text)) then
     write(LUPRI,*) 'new string lengths:',len(text),self(1)%str_length
     call FLSHFO(LUPRI)
     call QUIT('make_string: error in string length !')
   endif
   if ( alloc_stat > 0 ) then
     write(LUPRI,*) 'alloc_stat erro !'
     call FLSHFO(LUPRI)
     call QUIT('mystring/make_string: error in string allocation !')
     stop 
   endif

  ! copy the text-string into allocated array, pointed by strp(1)
  do i = 1, len(text)
    self(1)%strp(i) = text(i:i)
  enddo
  self(1)%is_alloc=.true.
 else
     write(LUPRI,*) 'mystring/make_string: this string '// &
     'is already allocated, see: ',self(1)%strp
 endif
 end subroutine make_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine cancel_string(self)
!----------------------------------------------------------------
!   cancels (deallocates) the string
!----------------------------------------------------------------
implicit none
#include "priunit.h"
 type(string), dimension(:), intent(inout) :: self
 integer :: dealloc_stat

 if ( self(1)%is_alloc) then
     deallocate(self(1)%strp, stat = dealloc_stat)
     if (dealloc_stat /= 0) then
       write(LUPRI,*) 'dellocation error !, dealloc_stat=',dealloc_stat
       write(LUPRI,*) ' length of the delloc. string was:',self(1)%str_length
       call FLSHFO(LUPRI)
       call QUIT('mystring/cancel_string: deallocation error !')
     endif
     nullify( self(1)%strp )
     self(1)%is_alloc=.false.
     self(1)%str_length = 0
 endif
 end subroutine cancel_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_string(self)
!----------------------------------------------------------------------------
!    prints out the whole string (into LUPRI channel)
!----------------------------------------------------------------------------
implicit none
#include "priunit.h"
type(string), dimension(:), intent(in) :: self
  if (self(1)%is_alloc) then
    write(LUPRI,*) (self(1)%strp)
  endif
  call FLSHFO(LUPRI)
end subroutine print_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 function is_string(self) result(is)
!---------------------------------------------------------------------------
! Returns true/false depending on whether the string was allocated/deallocated
!----------------------------------------------------------------------------
implicit none
   type(string), dimension(:), intent(in) :: self
   logical (kind=1) :: is
     is = self(1)%is_alloc
 end function is_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_string(self) result(str_out)
!-------------------------------------------------------------------
!
!   returns the string - pointer to the character array
!
!-------------------------------------------------------------------
  implicit none
  type(string), dimension(:), intent(in) :: self
  character*1, dimension (:), pointer :: str_out  
  integer :: i

  ! first clean returning pointer
  if ( associated(str_out) ) then
   deallocate(str_out); nullify(str_out)
  endif

  if (self(1)%is_alloc) then
   allocate( str_out(self(1)%str_length) )
   do i = 1, self(1)%str_length
    str_out(i:i) = self(1)%strp(i)
   enddo
  else
    str_out=' '
  endif
end function get_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine copy_strings(source,dest)
!-----------------------------------------------------------------------
!
!  copy allocated "source" -> "dest"; allocate "dest" first
!
!-----------------------------------------------------------------------
implicit none
#include "priunit.h"
type(string), dimension(1), intent(inout) :: source, dest
integer :: ierr,i
if (.not.source(1)%is_alloc) then
   call Quit('copy_strings: "source" is NOT allocated !')
endif
if (associated(dest(1)%strp)) then
  call cancel_string(dest)
endif
allocate( dest(1)%strp(source(1)%str_length), stat =  ierr )
if (ierr /= 0) then
  call Quit('copy_strings: allocation error of "dest(1)%strp" !')
endif
dest(1)%is_alloc   = .true.
dest(1)%str_length = source(1)%str_length
!...copy the text
do i = 1, dest(1)%str_length
 dest(1)%strp(i) = source(1)%strp(i)
enddo
end subroutine copy_strings

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_no_string(self)
!------------------------------------------------------------------------
!    sets logical flag indicating allocated string to .FALSE.
!------------------------------------------------------------------------
implicit none
 type(string), dimension(:), intent(inout) :: self
 self(1)%is_alloc = .FALSE.
 self(1)%str_length = 0
end subroutine set_no_string

end module mystring


! replacement of mystring module; mystring shall be canceled later (on TODO) 
module dyn_string
 type string
  private
   character*1, dimension (:), pointer :: strp  => NULL()
   logical :: is_alloc = .false. ! this default does not work for older F90 compilers 
   integer :: str_length = 0 ! length of the allocated string
 end type string
 integer, private :: iprint_string = 0  ! control print-out variable
 contains 

subroutine create_string(self, text)
implicit none
#include "priunit.h"
 type(string), intent(inout) :: self
 character(*), intent(in) :: text
 integer :: i, alloc_stat

 if (.not.self%is_alloc) then
   allocate(self%strp(len(text)), stat = alloc_stat)
   ! ... store the actual size
   self%str_length = size(self%strp)
   if (self%str_length.ne.len(text)) then
     write(LUPRI,*) 'new string lengths:',len(text),self%str_length
     call FLSHFO(LUPRI)
     call QUIT('make_string: error in string length !')
   endif
   if ( alloc_stat > 0 ) then
     write(LUPRI,*) 'alloc_stat erro !'
     call FLSHFO(LUPRI)
     call QUIT('mystring/make_string: error in string allocation !')
     stop 
   endif

  ! copy the text-string into allocated array, pointed by strp(1)
  do i = 1, len(text)
    self%strp(i) = text(i:i)
  enddo
  self%is_alloc=.true.
 else
    write(LUPRI,*) 'mystring/make_string: this string '// &
    'was already allocated, here: ',self%strp
 endif
end subroutine 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_string(self) result(str_out)
!-------------------------------------------------------------------
!   returns the string - pointer to the character array
!-------------------------------------------------------------------
  implicit none
  type(string), intent(in) :: self
  character*1, dimension (:), pointer :: str_out  
  integer :: i

  ! first clean returning pointer
  if ( associated(str_out) ) then
     nullify(str_out)
  endif

  if (self%is_alloc) then
   allocate( str_out(self%str_length) )
   do i = 1, self%str_length
    str_out(i:i) = self%strp(i)
   enddo
  else
    str_out=' '
  endif
end function get_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine cancel_string(self)
!----------------------------------------------------------------
!   cancels (deallocates) the string
!----------------------------------------------------------------
implicit none
#include "priunit.h"
 type(string), intent(inout) :: self
 integer :: dealloc_stat

 if ( self%is_alloc) then
     deallocate(self%strp, stat = dealloc_stat)
     if (dealloc_stat /= 0) then
       write(LUPRI,*) 'dellocation error !, dealloc_stat=',dealloc_stat
       write(LUPRI,*) ' length of the delloc. string was:',self%str_length
       call FLSHFO(LUPRI)
       call QUIT('mystring/cancel_string: deallocation error !')
     endif
     nullify( self%strp )
     self%is_alloc=.false.
     self%str_length = 0
 endif
end subroutine cancel_string

end module dyn_string
