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

! this file by ulfek march 2008

program dirac_program

  use os_utils
  use dft_cfg
  use xmlout
  use interface_ao
  use interface_mo
  use memory_allocator
  use machine_parameters
  use dirac_cfg

use iso_c_binding

#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#ifdef VAR_MPI
#include "mpif.h"
#endif
#endif

interface
   subroutine print_info() bind(c)
   end subroutine
end interface

#include "priunit.h"
#include "infpar.h"
#include "dcbgen.h"
#include "dcbpsi.h"


  real(kind=8), allocatable :: memory(:)

  integer :: iparcal, status, length, ios, idummy, ierr
  integer :: WORK_memsize
  integer(kind=8) :: maxsize = -1
  integer(kind=8) :: iscr_mem
  integer :: high_water_mark

! initialize common variables (parallel and sequential run)
  mytid   = 0
  mparid  = 0
  numnod  = 0
  parcal  = .false.
  iparex  = 0
! detect machine parameters
  call get_machine_parameters()
#ifdef VAR_MPI
! check for sequential run with parallel executable 
! iparex = 1 for parallel run mode with MPI executable (default)
  iparex = get_environment_integer('DIRPAR',1)
!
! initialize MPI environment and wake up slaves -
! iff numnod > 0 set parcal == .true.
  if (iparex .gt. 0) call mpixinit()
  if (parcal) then 
    call stopwatch(1,numnod+1)
  else
    call stopwatch(1,numnod)
  endif
#else
! check for parallel run with sequential executable 
! this may happen if the user has set --mpi in the .diracrc file
! we should then quit and report an error message
  iparex = get_environment_integer('DIRPAR',0)
  if (iparex .ne. 0) then
     print*, ' This is a serial executable, do not set --mpi in pam or .diracrc '
     call quit (' Attempted parallel run with serial executable')
  end if
  call stopwatch(1,numnod)
#endif
  iparcal = 0
  if (parcal) iparcal = 1

  if (mytid == 0) then
     iscr_mem     = get_environment_integer('DIRWRK',INSTALL_WRKMEM)
     WORK_memsize = iscr_mem
     if (WORK_memsize .ne. iscr_mem) then
        print *,'Problem reading DIRWRK:',WORK_memsize,iscr_mem
        call quit('Invalid entry in environment variable DIRWRK')
     end if
     maxsize = get_environment_integer('DIRMAX',maxsize)

  else
!    iscr_mem     = get_environment_integer('DIRNOD',INSTALL_WRKMEM) 
!    SK - Sep 08: DIRNOD is passed to slaves in mpixinit and 
!                 stored on common block GNCPAM (dcbgen.h)
     if( CDIRNOD(1:LENDND) .eq. ' ' )then
       WORK_memsize = INSTALL_WRKMEM
     else
       read (CDIRNOD(1:LENDND),*,iostat=ios) iscr_mem
       if(ios/=0)then
         print *,'Warning: Cannot parse environment variable DIRNOD'
         print *,'         Will use the default value:',INSTALL_WRKMEM
         WORK_memsize = INSTALL_WRKMEM
       else
         WORK_memsize = iscr_mem
       endif
     endif
     maxsize = get_environment_integer('DIRMAX',maxsize)
  endif

  call allocator_init()
  call allocator_set_max(maxsize)
  call allocator_get_max_words(maxsize)

#ifdef VAR_MPI
  if (parcal) then
    call sync_nodes(int(WORK_memsize,8),0)
    call sync_nodes(maxsize,2)
    call legacy_lwork_set(WORK_memsize)
  else
#else
  if (.not. parcal) then
#endif
    call legacy_lwork_set(WORK_memsize)
    print '(a,i0,a,f10.2,a,f10.3,a,$)', 'DIRAC serial starts by allocating ',WORK_memsize, &
    ' words (',DFLOAT(WORK_memsize)/DFLOAT(2**17),' MB - ', DFLOAT(WORK_memsize)/DFLOAT(2**27),' GB)' ,' of memory'
    if (maxsize.gt.0) then
       print '(a,i0,a,f10.2,a,f10.3,a)', '    out of the allowed maximum of ',  &
          maxsize,' words (',DFLOAT(maxsize)/DFLOAT(2**17),' MB - ',DFLOAT(maxsize)/DFLOAT(2**27),' GB)'
    else
       print '(a)', 'DIRAC serial has no limitations in place for the amount of dynamically allocated memory'
    endif
    print '(/a/)','Note: maximum allocatable memory for serial run can be set by pam --aw/--ag'
  end if
  if (mytid == 0) then
     ! test if memory is available by allocating and deallocating
     call alloc(memory,WORK_memsize, id="test allocation of work array in DIRAC main program")
     call dealloc(memory)
     call printtitle
     call printlogo
     call printsubtitle
     call flush(lupri)
     ! print out git_info & build_info by the executable
     call print_info()
     call flush(lupri)
! miro: do the self-tests after printing dirac system info
#ifndef NO_SELFTEST
     call selftest_lapack()
     call selftest_fortran_c_cxx_interoperability()
#ifdef VAR_MPI
!    checks MPI library compatibility 
     endif
     if ( parcal ) call selftest_mpi()
     if (mytid == 0) then
#endif
#endif
! debug     print*,'++++ calling dirac with',iparcal,mytid,mparid,numnod
     call dirac(iparcal,mytid,mparid,numnod)
  else
#ifdef VAR_MPI
     ! test if memory is available by allocating and deallocating
     call alloc(memory,WORK_memsize, id="test allocation of work array at main for slave")
     call dealloc(memory)
! debug     print*,'++++ calling dirnod with',iparcal,mytid,mparid,numnod
     call dirnod(iparcal,mytid,mparid,numnod)
#else
     call quit('Internal error in main.F90. Quitting.')
#endif
  endif
#ifdef VAR_MPI
  if ( iparcal .eq. 1 ) then
    call sync_nodes(idummy,1)
  endif
#endif
  if (mytid == 0) then
     write (lupri, '(/a)') '*****************************************************'
     write (lupri, '(a)')  '********** E N D   of   D I R A C  output  **********'
     write (lupri, '(a/)') '*****************************************************'
     CALL TSTAMP(' ',LUPRI)
     write (lupri,'()')
     flush(lupri)
     call stopwatch(0,mytid)
     flush(lupri)
  endif
#ifdef VAR_MPI
  if (iparex .gt. 0) call mpixfinalize()
#endif

 call allocator_report(mytid,0)
  if (mytid == 0) then
      write (lupri, '(A,F8.2,A)') ' MEMGET high-water mark: ',(high_water_mark()*8)*2.0D0**(-20),' MB'
      write (lupri, '(/a)' ) '*****************************************************'
  endif

   call allocator_cleanup()

   call interface_ao_clear()
   call interface_mo_clear()

end program
