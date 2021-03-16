!      Copyright (c) 2018 by the authors of DIRAC.
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

! $Id$
!
!
! $Id$
!
! memory allocation subsystem
!
! initially conceived to be used in the dirac and dalton packages
! as a replacement to the memget/memrel routines and related
!
! design/discussions    : ulf ekstrom, andre gomes, jetze sikkema
! inital implementation : andre gomes
!                         vu amsterdam, winter/spring 2008
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test program for memory allocator for dirac/dalton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 program alloc_test
   call andre_alloc_test
#if defined INT_STAR8
   call test_xcopy_with_big_array
#endif
 end program alloc_test

 subroutine andre_alloc_test
!miro: adapted for current version of memory_allocator
    use memory_allocator
    use os_utils

    integer, parameter :: kreal = 8, kcomplex = 8
    integer, parameter :: charArrayLength = 132
    integer                :: v2size, v1size, p_mem2,p_mem3
    integer                :: maxsize, minimal = 20, verbose = 0
    double precision, allocatable, target :: vdatadb(:), vdatarl(:)
    real(kind=kreal), allocatable :: matrix(:,:)
    real(kind=kreal), pointer             :: p_mem, p_mem42(:)
    complex(kind=kreal), allocatable :: cmatrix(:,:)
    character(len=charArrayLength) :: kbytesAvailable
    logical :: test_huge_alloc=.False.

 ! maxsize to be extracted from environmental variable
 maxsize = get_environment_integer('DIRMAX',maxsize)
 print *,'variable from DIRMAX, maxsize=',maxsize
         !maxsize   = 1500
         !maxsize   = 15000 ! MI: increase
         minimal   = 16
         nstep     = 10
         direction = 1

         v2size = 10

!        write (*,*) 'memory address 0a: ',c_loc(vdatarl)
!        write (*,*) 'memory address 0b: ',c_loc(vdatadb)

 call allocator_init
 call allocator_set_max(maxsize)

!         call allocator_setGroupName('main')

         call alloc(vdatadb,v2size,id="vdatadb")
         call alloc(vdatarl,v2size)
         call allocp(p_mem42,v2size)
         call alloc(matrix,v2size,v2size,id="real matrix")
         call alloc(cmatrix,v2size,v2size,id="complx matrix")

         p_mem => vdatadb(0)

         call get_memory_address(vdatadb, p_mem2)
         write (*,*) 'memory address 2: [',p_mem2,']'
         call get_memory_address(vdatarl, p_mem3)
         write (*,*) 'memory address 3: [',p_mem3,']'

         forall (i=1:v2size, j=1:v2size) matrix(i,j) = 0.001

         do i=1,v2size
            do j=1,v2size 
               write(*,*) 'm(',i,',',j,') =',matrix(i,j)
            enddo
         enddo

!        write (*,*) 'memory address 4: ',c_loc(vdatarl)
!        write (*,*) 'memory address 5: ',c_loc(vdatadb)

         call compare_memory(p_mem2,p_mem3)

         write (*,*) 'outside rbla'
         do i=1,v2size
            vdatarl(i) = 0.002 
            write (*,*) 'prep data ',vdatarl(i)
         enddo

         write (*,*) 'before rbla'
         call rbla(vdatarl,v2size)
         write (*,*) 'before dbla'
         call dbla(vdatarl,vdatadb,3)

         do i=1,v2size
            write (*,*) 'd(',i,')=',vdatadb(i),' r(',i,')=',vdatarl(i)
            write (*,*) 'kind(d) =',kind(vdatadb(i))
            write (*,*) 'kind(r) =',kind(vdatarl(i))
         enddo

         call dealloc(vdatadb)
         call dealloc(vdatarl,id="vdatarl")

         call get_memory_address(vdatadb, p_mem2)
         write (*,*) 'memory address 2: [',p_mem2,']'
         call get_memory_address(vdatarl, p_mem3)
         write (*,*) 'memory address 3: [',p_mem3,']'
         call compare_memory(p_mem2,p_mem3)
!        write (*,*) 'memory address 4: ',c_loc(vdatarl)
!        write (*,*) 'memory address 5: ',c_loc(vdatadb)
         write (*,*) 'is 4 allocated? ',allocated(vdatarl)
         write (*,*) 'is 5 allocated? ',allocated(vdatadb)
         write (*,*) 'is p_mem42 allocated? ',associated(p_mem42)

      !   call allocator_PrintStats()

         call get_memory_info(kbytesAvailable,charArrayLength)
         ! MI: how to find out the amount of free memory ?
         write(*,*) 'before large allocations:', & 
                   kbytesAvailable,' kbytes available'
         !call catfile() 

         write (*,*) ' '
         write (*,*) 'testing alloc_maxsize, no previous allocations'
         write (*,*) '   minimal mwords required   :',minimal
         write (*,*) '   maximum mwords start      :',maxsize
         !call alloc_maxsize(maxsize,minimal,verbose) 
         call allocator_set_max(maxsize)
         write (*,*) '   maximum mwords end        :',maxsize
         write (*,*) ' '

         v2size = 1000
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size,id="datadb30")
         write (*,*) 'is it allocated? ',allocated(vdatadb)
 
         call foo()

         write (*,*) ' '
         write (*,*) 'testing alloc_maxsize, allocation of',(v2size/(1204*1024*8)),'mwords active' 
         write (*,*) '   minimal mwords required   :',minimal
         write (*,*) '   maximum mwords start      :',maxsize
         !call alloc_maxsize(maxsize,minimal,verbose) 
         call allocator_set_max(maxsize)
         write (*,*) '   maximum mwords end        :',maxsize
         write (*,*) ' '

         call dealloc(vdatadb)
          
         v2size = 10000
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size,id="v2size_10000")

         maxsize = 0
         write (*,*) ' '
         write (*,*) 'testing alloc_maxsize, allocation of',(v2size/(1204*1024*8)),'mwords active' 
         write (*,*) '   minimal mwords required   :',minimal
         write (*,*) '   maximum mwords start      :',maxsize
         !call alloc_maxsize(maxsize,minimal,verbose)
         call allocator_set_max(maxsize)
         write (*,*) '   maximum mwords end        :',maxsize
         write (*,*) ' '


         write (*,*) 'is it allocated? ',allocated(vdatadb)
         call dealloc(vdatadb) ! mi: 2 times, it should continue

         print *,'before 2. large allocations...available: ', &
                   kbytesAvailable,' kbytes'
          
         v2size = 100000
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size)

         maxsize = 0
         write (*,*) ' '
         write (*,*) 'testing alloc_maxsize, allocation of',(v2size/(1204*1024*8)),'mwords active' 
         write (*,*) '   minimal mwords required   :',minimal
         write (*,*) '   number of steps attepmted :',nstep
         write (*,*) '   direction of search       :',direction
         write (*,*) '   maximum mwords start      :',maxsize
         !call alloc_maxsize(maxsize,minimal,verbose)
         call allocator_set_max(maxsize)
         write (*,*) '   maximum mwords end        :',maxsize
         write (*,*) ' '

         write (*,*) 'is it allocated? ',allocated(vdatadb)
         call dealloc(vdatadb)
         call dealloc(vdatadb)
          
         v2size = 1000000
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size,id="reallylarge")

         maxsize = 0
         write (*,*) ' '
         write (*,*) 'testing alloc_maxsize, allocation of',(v2size/(1204*1024*8)),'mwords active' 
         write (*,*) '   minimal mwords required   :',minimal
         write (*,*) '   number of steps attepmted :',nstep
         write (*,*) '   direction of search       :',direction
         write (*,*) '   maximum mwords start      :',maxsize
         !call alloc_maxsize(maxsize,minimal,verbose)
         call allocator_set_max(maxsize)
         write (*,*) '   maximum mwords end        :',maxsize
         write (*,*) ' '
         write (*,*) 'is it allocated? ',allocated(vdatadb)
         call dealloc(vdatadb)

  if (test_huge_alloc) then
         v2size = huge(maxsize) !  value of 9223372036854775807, not allocatable
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size)
         maxsize = 0
         write (*,*) ' '
         write (*,*) 'testing alloc_maxsize, allocation of',(v2size/(1204*1024*8)),'mwords active' 
         write (*,*) '   minimal mwords required   :',minimal
         write (*,*) '   number of steps attepmted :',nstep
         write (*,*) '   direction of search       :',direction
         write (*,*) '   maximum mwords start      :',maxsize
         !call alloc_maxsize(maxsize,minimal,verbose)
         call allocator_set_max(maxsize)
         write (*,*) '   maximum mwords end        :',maxsize
         write (*,*) ' '
         write (*,*) 'is it allocated? ',allocated(vdatadb)
         call dealloc(vdatadb)
  endif

 call allocator_cleanup()

 end subroutine andre_alloc_test


      subroutine compare_memory(location_1, location_2)
           integer           :: location_1, location_2 

           write (*,*) 'a=',location_1,'  b=',location_2
           if (location_1 .eq. location_2) then
              write (*,*) 'addresses are the same'
           else
              write (*,*) 'addresses are different'
           endif

      end subroutine compare_memory

      subroutine rbla(data, size)
           character(len=1)  :: dummy
           integer           :: size
           real              :: data(*)

!           call allocator_setGroupName('rbla')
           write (*,*) 'inside rbla 1'
           do i=1,size
              write (*,*) ' data ',data(i)
           enddo
           write (*,*) 'inside rbla 2'
           do i=1,size
              data(i) = 2.1 
              write (*,*) ' new data ',data(i)
           enddo
           write (*,*) 'leaving rbla'
!           call get_traceback_info (printPrettyTrace,dummy,1)
!           call get_traceback_info (printAllocCaller,dummy,1)
      end subroutine rbla

      subroutine dbla(data1, data2, size)
           integer           :: size
           double precision  :: data2(*)
           real              :: data1(*)

           write (*,*) 'inside dbla'
           do i=1,size
              data2(i) = data1(i) / 1000
           enddo
           write (*,*) 'leaving dbla'
      end subroutine dbla

      subroutine foo
         character bla*1
         write (*,*) 'first the pretty print'
         CALL GET_TRACEBACK_INFO (3,BLA,1)
         write (*,*) 'then the basic print'
         CALL GET_TRACEBACK_INFO (4,BLA,1)
         write (*,*) 'and we are done'
      end subroutine foo
     
#if defined INT_STAR8
subroutine test_xcopy_with_big_array
!-------------------------------------------------------
! Miro Ilias routine for testing some alloc stuff
!-------------------------------------------------------
use memory_allocator
use os_utils
 implicit none
 integer, parameter :: kreal    = 8
 real (kind=kreal), pointer ::     a_v4(:)
 real (kind=kreal), allocatable, target :: a_vt(:)
!---------------------
! COMMON /COMPI/ RCW
! INTEGER RCW
!---------------------
#include "../src/relccsd/complex.inc"
 integer :: NVTMIN 
 integer(kind=8) :: maxsize = -1

 RCW = 1
 NVTMIN = 2147483948 ! ! larger than 2^31-1=2147483648
 print *,'NVTMIN=',NVTMIN,' RCW=',RCW
 print *,'a_VT of size NVTMIN*RCW=',NVTMIN*RCW

 maxsize = get_environment_integer('DIRMAX',maxsize)
 print *,'evironmental variable, DIRMAX=',maxsize
 call allocator_init
 call allocator_set_max(maxsize)

 !allocate(a_VT(NVTMIN*RCW))
 call alloc(a_VT ,NVTMIN*RCW, id="vt" )
 
 a_V4  => a_VT
 call my_rdints(a_V4)

 call allocator_cleanup()

 end subroutine test_xcopy_with_big_array

 subroutine my_rdints(VVVVV)
implicit none
 real*8 :: A0r=0.0d0
 complex*16 :: A0=(0.0d0,0.0d0)
 integer :: N = 2147483650 ! larger than 2^31-1
!---------------------
! COMMON /COMPI/ RCW
! INTEGER RCW
!---------------------
#include "../src/relccsd/complex.inc"
 REAL*8 :: VVVVV(RCW,*)
 print *,'rdints: befor call dcopy RCW=',RCW
 print *,' big N=',N
 print *,' A0r=',A0r
 print *,' A0=',A0
 call xcopy(N,A0r,0,VVVVV,1)
 print *,' after xcopy with real A0r, big N'
 call xcopy(N,A0,0,VVVVV,1)
 print *,' after xcopy with complex A0, big N'

 return
 end subroutine my_rdints
#endif
 
