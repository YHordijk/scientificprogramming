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

! This module contains functions to communicate with the operating system.
! ulfek march 2008
module os_utils

       private

       public get_environment_integer
    
       interface get_environment_integer
          module procedure get_environment_integer_i4
          module procedure get_environment_integer_i8
       end interface get_environment_integer

       contains
! Gets an environment variable and converts it to an integer.
! Return default if the variable is not defined or conversion fails
function get_environment_integer_i4(name, default)
  implicit none
  integer(kind=4) :: default, get_environment_integer_i4, length, status, i, ios
  character(len=*) :: name
  character(len=128) :: value
#ifndef FORTRAN2003  
  call getenv(name, value)
  if (len_trim(value) == 0) then
     status = 1
  else  
     status = 0 
  endif
#else
  call get_environment_variable(name, value, length, status)
#endif  
  if (status==0) then
     read (value,*,iostat=ios) i
     if (ios/=0) then
        print *,'ERROR reading environment variable ',name
        print *,'ERROR cannot convert to 32 bit integer, value = ',value
        call quit('Invalid environment variable value')
        get_environment_integer_i4 = default
     else
        get_environment_integer_i4 = i
     endif
  else
        get_environment_integer_i4 = default
  endif
end function get_environment_integer_i4


function get_environment_integer_i8(name, default)
  implicit none
  integer(kind=8) :: default, get_environment_integer_i8, i
  integer(kind=4) :: length, status, ios
  character(len=*) :: name
  character(len=128) :: value
#ifndef FORTRAN2003  
  call getenv(name, value)
  if (len_trim(value) == 0) then
     status = 1
  else
     status = 0
  endif
#else
  call get_environment_variable(name, value, length, status)
#endif  
  if (status==0) then
     read (value,*,iostat=ios) i
     if (ios/=0) then
        print *,'ERROR reading environment variable ',name
        print *,'ERROR cannot convert to 64 bit integer, value = ',value
        call quit('Invalid environment variable value')
        get_environment_integer_i8 = default
     else
        get_environment_integer_i8 = i
     endif
  else
        get_environment_integer_i8 = default
  endif
end function get_environment_integer_i8


end module os_utils
