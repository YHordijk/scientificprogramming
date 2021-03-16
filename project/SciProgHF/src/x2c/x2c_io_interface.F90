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
!
! file-I/O routines used in the X2C routines.
!
! written by sknecht april 2010
!
module x2c_fio

  implicit none

  public x2c_write
  public x2c_read

  private

contains

!**********************************************************************
  subroutine x2c_write(flabel,fmat,nlen,funit)
!**********************************************************************
     real(8),            intent(in)     :: fmat(*)
     integer,            intent(in)     :: nlen
     integer,            intent(in)     :: funit
     character (len=12), intent(in)     :: flabel
!----------------------------------------------------------------------
     logical                            :: fndlab12
     character (len=800)                :: funit_name
     integer                            :: lfunit_name
!**********************************************************************

       if(nlen >= 0)then
         rewind funit                                                    
         if(fndlab12('EOFLABEL-x2c',funit))then
           backspace(funit)
           call newlab12(flabel,funit,6)
           call writt(funit,nlen,fmat)
           call newlab12('EOFLABEL-x2c',funit,6)
         else
           inquire(unit=funit, name=funit_name)
           lfunit_name = len_trim(funit_name)
           print '(/2x,a,a,a)', '*** error in x2c_write: end-of-file label is missing in the file ',            &
                              funit_name(1:lfunit_name),'. The program therefore stops. ***'                            
           call quit('*** error in x2c_write: end-of-file label is missing. ***')
         end if
       else 
         call newlab12('EOFLABEL-x2c',funit,6)
       end if

  end subroutine x2c_write

!**********************************************************************
  subroutine x2c_read(flabel,fmat,nlen,funit)
!**********************************************************************
     real(8),            intent(inout)  :: fmat(*)
     integer,            intent(in)     :: nlen
     integer,            intent(in)     :: funit
     character (len=12), intent(in)     :: flabel
!----------------------------------------------------------------------
     logical                            :: fndlab12
     character (len=800)                :: funit_name
     character (len=1000)               :: error_stream
     integer                            :: lfunit_name
     integer                            :: lerror_stream
!**********************************************************************
 
       rewind(funit)
       if(nlen > 0)then
         if(fndlab12(flabel,funit))then
           call readt(funit,nlen,fmat)
         else
           inquire(unit=funit, name=funit_name)
           lfunit_name = len_trim(funit_name)
           print '(a,a,a,a)', '   *** error in x2c_read: data with label ',flabel,' is not present in the file ',  &
                                  funit_name(1:lfunit_name)
           write(error_stream,'(a,a,a,a)') & 
           '   *** error in x2c_read: data with label ',flabel,' is not present in the file ',funit_name(1:lfunit_name)
           lerror_stream = len_trim(error_stream)
           call quit(error_stream(1:lerror_stream))
         end if
       else
         print '(2x,a)', ' *** warning in x2c_read: attempt to read data array which is supposed to be zero. ***'
!        call quit('bla bla')
       end if
       
  end subroutine x2c_read
!**********************************************************************

end module x2c_fio
