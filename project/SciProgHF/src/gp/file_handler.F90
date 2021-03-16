!==============================================================================!
module file_handler
!------------------------------------------------------------------------------!
!
! module for I/O 
!
!  (1) managing unit numbers of open files
!
!  (2) open files in sequential or direct access mode
!
!
!
!  BHP, summer 2015
!
!------------------------------------------------------------------------------!

 implicit none

 private

 ! maximum number of open files
 integer, parameter :: iunend  = 99
 
 ! smallest possible unit number
 integer, parameter :: iunstr = 10 

 ! maximum file length
 integer, parameter :: mxfilen = 256

 ! multiplication factor for direct-access files
 integer, parameter :: mulrec = 8

 ! default file name if unitialized
 character(len=mxfilen), parameter :: filuninit = " --- uninitialized file name ---"

 ! names of files
 character(len=mxfilen) :: files(iunstr:iunend) = filuninit

 public :: mulrec, mxfilen, getunit, releaseunit, open_seq, prt_open, open_da, &
           remove_file

 contains

!==============================================================================!
integer function getunit(filnam)

 implicit none

#include "stdunit.h"

! parameter:
 character(len=*), parameter :: chrout = 'getunit>'

! input:
 character(len=*),intent(in) :: filnam

! local:
 integer :: lunit, ipos
 character(len=mxfilen) :: filtmp

 filtmp = adjustl(filnam)

 do lunit = iunstr,iunend
  ipos = index(files(lunit),filuninit)
  if (ipos.gt.0) then
   files(lunit) = trim(filtmp)
   exit
  end if
 end do

 if (lunit.gt.iunend) call quit(chrout//"Too many files open. Increase iunend!")

 getunit = lunit

end function getunit
!==============================================================================!

!==============================================================================!
integer function open_seq(filnam,cform,caction)
   
 implicit none

#include "stdunit.h"

! parameter:
 character(len=*), parameter :: chrout = 'open_seq>'

! input:
 character(len=*),intent(in) :: filnam
 character(len=*),intent(in),optional :: cform, caction

! local:
 logical :: lopen
 integer :: lunit, ipos, ios
 character(len=16) :: caccess
 character(len=11) :: cform_loc
 character(len=9)  :: caction_loc
 character(len=mxfilen) :: filtmp

 filtmp = adjustl(filnam)
 lunit = getunit(trim(filtmp))

 ! check if file is already open
 inquire(file=trim(filtmp),opened=lopen,access=caccess,form=cform_loc,action=caction_loc,iostat=ios)
 if (ios.ne.0) then
  call prt_open()
  call quit('Cannot inquire file '//trim(filtmp)//'!')
 end if

 if (lopen) then
  ipos = index(trim(caccess),'sequential')
  if (ipos.eq.0) then
   call prt_open()
   call quit(trim(filtmp)//' is not a sequential file!')
  end if
  
  if (present(cform)) then
   ipos = index(trim(cform),trim(cform_loc))
   if (ipos.eq.0) then
    call prt_open()
    call quit('Format inconsistent in open '//trim(filtmp)//'!')
   end if
  end if

  if (present(caction)) then
   ipos = index(trim(caction),trim(caction_loc))
   if (ipos.eq.0) then
    call prt_open()
    call quit('Action inconsistent in open '//trim(filtmp)//'!')
   end if
  end if

  rewind(unit=lunit,iostat=ios)
  if (ios.ne.0) then
   call prt_open()
   call quit('Cannot rewind file '//trim(filtmp)//'!')
  end if

 else

  if (present(cform)) then
   cform_loc = trim(cform)
  else
   cform_loc = "UNFORMATTED"
  end if

  if (present(caction)) then
   caction_loc = trim(caction)
  else
   caction_loc = "READWRITE"
  end if

  
  if (caction_loc(1:9).eq."READ     ") then
   open(unit=lunit,file=trim(filtmp),status='old',access='sequential',form=cform_loc,action=caction_loc,iostat=ios)
  else
   open(unit=lunit,file=trim(filtmp),access='sequential',form=cform_loc,action=caction_loc, iostat=ios)
  end if

  if (ios.ne.0) then
   call prt_open()
   call quit('Cannot open file '//trim(filtmp)//'!')
  end if
 end if

 open_seq = lunit

end function open_seq
!==============================================================================!

!==============================================================================!
subroutine releaseunit(lunit,cstatus)

 implicit none

#include "stdunit.h"

! 
 character(len=*), parameter :: chrdbg = "releaseunit>"

! input:
 integer,intent(in) :: lunit
 character(len=*),intent(in),optional :: cstatus

! local:
 integer :: ipos, ios
 logical :: lopen
 character(len=mxfilen) :: filtmp
 character(len=6) :: cstatus_loc

 ! check range of unit number
 if (lunit.lt.iunstr .or. lunit.gt.iunend) then
  write(6,*) chrdbg,"lunit:",lunit,iunstr,iunend
  call prt_open()
  call quit("unit number out of range!")
 end if

 ipos = index(files(lunit),filuninit)
 if (ipos.gt.0) then
  call prt_open()
  call quit("unit number is not assigned to any file!")
 end if

 inquire(unit=lunit,opened=lopen,name=filtmp,iostat=ios)
 if (ios.ne.0) then
  call prt_open()
  call quit('Cannot inquire unit number!')
 end if

! ipos = index(files(lunit),trim(filtmp))
! ipos = max(ipos,index(trim(filtmp),files(lunit)))
! if (ipos.eq.0) then
!  call prt_open()
!  write(6,*) 'debug:', files(lunit),trim(filtmp),lunit
!  call quit("unit number assigned to file is inconsistent with file record!")
! end if
!
 ! check if status is available and correct
 if (present(cstatus)) then
  select case (trim(cstatus))
   case ('keep','KEEP')
    cstatus_loc(1:6) = 'keep  '
   case ('delete','DELETE')
    cstatus_loc(1:6) = 'delete'
   case default
    call quit("Unknown status for closing unit!")
  end select
 else
  cstatus_loc(1:6) = 'keep  '
 end if

 if (lopen) then
  close(unit=lunit,iostat=ios,status=cstatus_loc)
  if (ios.ne.0) then
   call prt_open()
   call quit('Cannot close file!')
  end if
 end if

 files(lunit) = filuninit

end subroutine releaseunit
!==============================================================================!
  
!==============================================================================!
subroutine prt_open()

 implicit none

#include "stdunit.h"

! local:
 integer :: lunit, ipos

 write(istdout,'(/,a)') '--------------------------------------------'
 write(istdout,  '(a)') '  unit num.            file name            '
 write(istdout,  '(a)') '--------------------------------------------'
 do lunit = iunstr, iunend
  ipos = index(files(lunit),filuninit)
  if (ipos.eq.0) then
   write(istdout,'(4x,i4,4x,a32)') lunit,files(lunit)
  end if
 end do
 write(istdout,'(a,/)') '--------------------------------------------'

end subroutine prt_open
!==============================================================================!
!==============================================================================!
integer function open_da(filnam,lenrec,cstat)
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

 implicit none

#include "stdunit.h"

! parameter:
 character(len=*), parameter :: chrout = 'open_da>'

! input
 character(len=*), intent(in) :: filnam
 character(len=*), intent(in),optional :: cstat
 integer, intent(in) :: lenrec

! local
 integer :: lunit, junit, ios
 logical :: lopen
 character(len=mxfilen) :: filtmp
 character(len=6) :: cstatus_loc

 filtmp = adjustl(filnam)

 if (lenrec.lt.1) call quit('Record length should be >= 0!')

! ! check if status is available and correct
! if (present(cstatus)) then
!  select case (trim(cstatus))
!   case ('keep','KEEP')
!    cstatus_loc(1:6) = 'keep  '
!   case ('delete','DELETE')
!    cstatus_loc(1:6) = 'delete'
!   case default
!    call quit("Unknown status for closing unit!")
!  end select
! else
!  cstatus_loc(1:6) = 'keep  '
! end if

 lunit = getunit(trim(filtmp))

 ! check if file is already open
!$OMP CRITICAL
 inquire(file=trim(filtmp),opened=lopen,iostat=ios)
 if (ios.ne.0) then
  call prt_open()
  call quit('Cannot inquire open status of file '//trim(filtmp)//'!')
 end if

 if (lopen) then
  inquire(file=trim(filtmp),number=junit,iostat=ios)
  if (ios.ne.0) then
   call prt_open()
   call quit('Cannot inquire unit number of file '//trim(filtmp)//'!')
  end if

  if (lunit.ne.junit) then
   call quit('Unit number inconsistent!')
  end if

  rewind(unit=lunit,iostat=ios)
  if (ios.ne.0) then
   call prt_open()
   call quit('Cannot rewind file '//trim(filtmp)//'!')
  end if
 else
  open(unit=lunit,file=trim(filtmp),status='unknown',form='unformatted',&
       access='direct',recl=lenrec,iostat=ios)

  if (ios.ne.0) then
   call prt_open()
   call quit('Cannot open DA file '//trim(filtmp)//'!')
  end if
 end if
!$OMP END CRITICAL

 open_da = lunit

end function open_da
!==============================================================================!

!==============================================================================!
subroutine remove_file(filnam)
   
 implicit none

#include "stdunit.h"

! parameter:
 character(len=*), parameter :: chrout = 'remove_file>'

! input:
 character(len=*),intent(in) :: filnam

! local:
 logical :: lopen
 integer :: lunit, ios
 character(len=16) :: caccess
 character(len=mxfilen) :: filtmp

 filtmp = adjustl(filnam)
 lunit = getunit(trim(filtmp))

 ! check if file is already open
 inquire(file=trim(filtmp),access=caccess,opened=lopen,iostat=ios)
 if (ios.ne.0) then
  call prt_open()
  call quit('Cannot inquire file '//trim(filtmp)//'!')
 end if

 if (.not.lopen) then
  open(unit=lunit,file=trim(filtmp),status='old',iostat=ios)

  if (ios.ne.0) then
   call prt_open()
   call quit('Cannot open file '//trim(filtmp)//'!')
  end if
 end if

 close(unit=lunit,iostat=ios,status='delete')
 if (ios.ne.0) then
  call prt_open()
  call quit('Cannot close file!')
 end if

 files(lunit) = filuninit

end subroutine remove_file
!==============================================================================!

end module file_handler
!==============================================================================!
