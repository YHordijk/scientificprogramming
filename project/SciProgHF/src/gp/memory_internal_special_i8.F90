!*
!*
!* Copyright (c) 2008-2013 Andre Severo Pereira Gomes <andre.gomes@univ-lille1.fr>
!* All rights reserved.
!*
!* Redistribution and use in source and binary forms, with or without
!* modification, are permitted provided that the following conditions
!* are met:
!* 1. Redistributions of source code must retain the above copyright
!*    notice, this list of conditions and the following disclaimer.
!* 2. Redistributions in binary form must reproduce the above copyright
!*    notice, this list of conditions and the following disclaimer in the
!*    documentation and/or other materials provided with the distribution.
!*
!* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
!* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
!* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
!* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
!* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
!* SUCH DAMAGE.
!*
!*

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      module allocator_internal_A_1d_fixed_i8_indexing_R
      
         use allocator_parameters
         use allocator_cfg
         use allocator_track_if
         use allocator_eh
         use allocator_internal_init_R

         implicit none
         
         private
 
         public alloc_A_normal_1d_fixed_i8_indexing_R
         public alloc_A_sliced_1d_fixed_i8_indexing_R

         contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! allocation routines 
!

!
! one-dimension variable
!
            subroutine alloc_backend_A_1d_fixed_i8_indexing_R(data,&
     &           start_index,end_index,id)
               integer                               :: status
               integer(kind=8)                       :: size
               integer(kind=8)                       :: start_index
               integer(kind=8)                       :: end_index
               character(len=*)                      :: id
               integer(kind=klongint)                :: memAddr
               real(kind=kreal), allocatable, target   :: data(:)
              
               size = end_index - start_index + 1
               
               if ( size .eq. 0 ) then 
!
! calls with same start and index are dummy calls
! we go back to the caller after doing nothing
!
                  return
!
               else if (size.lt.0) then
!
! negative size calls are calls we can't satisfy (for instance, 
! we want to set a buffer of a given size but there's nearly no 
! memory left, or the size variable has under/overflowed
!
                  call allocator_errorHandler(negativeSizeCall,size, &
     &                                        kind(data),var_id=id)
!
               endif
!
! now check if the size we want to allocate will still give less 
! memory usage than what has been established as the ceiling 
!
               if ( .not. allocator_withinRange(size,kind(data)) ) then
                  call allocator_errorHandler(reqMemoryOverLimit,size, &
     &                                           kind(data),var_id=id)
               endif


               if (.not.varIsAliveA(data)) then
                  allocate(data(start_index:end_index),stat=status)

                  if ( status.ne.0 ) then
                     call allocator_errorHandler(allocateFailure,size, &
     &                                           kind(data),var_id=id)
                  else
!                    call initData_R(data,start_index,end_index)
                     call get_memory_address(data, memAddr)
                     call allocator_track(size,kind(data),memAddr,id)
                  endif
               else
                 call allocator_errorHandler(alreadyAllocated,var_id=id)
               endif

            end subroutine alloc_backend_A_1d_fixed_i8_indexing_R


!
! wrapper for the allocation routine
! handling arrays with start indexes equal to one
!
            subroutine     alloc_A_normal_1d_fixed_i8_indexing_R(data,size,id)
               integer(kind=8)                       :: size, start_index
               character(len=*), optional            :: id
               character(len=charArrayLength)        :: local_id
               real(kind=kreal), allocatable, target   :: data(:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'unnamed variable'
               endif
               start_index = 1
               call alloc_backend_A_1d_fixed_i8_indexing_R(data,start_index,size,local_id) 

            end subroutine alloc_A_normal_1d_fixed_i8_indexing_R


!
! wrapper for the allocation routine
! handling arrays with start indexes not equal to one
!
            subroutine     alloc_A_sliced_1d_fixed_i8_indexing_R(data,start_index,end_index,id)
               integer(kind=8)                       :: start_index, end_index
               character(len=*), optional            :: id
               character(len=charArrayLength)        :: local_id
               real(kind=kreal), allocatable, target   :: data(:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'unnamed variable'
               endif

               call alloc_backend_A_1d_fixed_i8_indexing_R(data,start_index,end_index,local_id)

            end subroutine alloc_A_sliced_1d_fixed_i8_indexing_R

      end module allocator_internal_A_1d_fixed_i8_indexing_R
!*
!*
!* Copyright (c) 2008-2013 Andre Severo Pereira Gomes <andre.gomes@univ-lille1.fr>
!* All rights reserved.
!*
!* Redistribution and use in source and binary forms, with or without
!* modification, are permitted provided that the following conditions
!* are met:
!* 1. Redistributions of source code must retain the above copyright
!*    notice, this list of conditions and the following disclaimer.
!* 2. Redistributions in binary form must reproduce the above copyright
!*    notice, this list of conditions and the following disclaimer in the
!*    documentation and/or other materials provided with the distribution.
!*
!* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
!* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
!* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
!* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
!* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
!* SUCH DAMAGE.
!*
!*


   
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      module allocator_internal_P_1d_fixed_i8_indexing_R
      
         use allocator_parameters
         use allocator_cfg
         use allocator_track_if
         use allocator_eh
         use allocator_internal_init_R

         implicit none
         
         private
 
         public alloc_P_normal_1d_fixed_i8_indexing_R
         public alloc_P_sliced_1d_fixed_i8_indexing_R

         contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! allocation routines 
!

!
! one-dimension variable
!
            subroutine alloc_backend_P_1d_fixed_i8_indexing_R(data,start_index,end_index,id)
               integer                               :: status
               integer(kind=8)                       :: size
               integer(kind=8)                       :: start_index, end_index
               character(len=*)                      :: id
               integer(kind=klongint)                :: memAddr
               real(kind=kreal), pointer   :: data(:)
              
               size = end_index - start_index + 1
               
               if ( size .eq. 0 ) then 
!
! calls with same start and index are dummy calls
! we go back to the caller after doing nothing
!
                  return
!
               else if (size.lt.0) then
!
! negative size calls are calls we can't satisfy (for instance, 
! we want to set a buffer of a given size but there's nearly no 
! memory left, or the size variable has under/overflowed
!
                  call allocator_errorHandler(negativeSizeCall,size, &
     &                                        kind(data),var_id=id)
!
               endif
!
! now check if the size we want to allocate will still give less 
! memory usage than what has been established as the ceiling 
!
               if ( .not. allocator_withinRange(size,kind(data)) ) then
                     call allocator_errorHandler(reqMemoryOverLimit,size, &
     &                                           kind(data),var_id=id)
               endif


               if (.not.varIsAliveP(data)) then
                  allocate(data(start_index:end_index),stat=status)

                  if ( status.ne.0 ) then
                     call allocator_errorHandler(allocateFailure,size, &
     &                                           kind(data),var_id=id)
                  else
!                    call initData_R(data,start_index,end_index)
                     call get_memory_address(data, memAddr)
                     call allocator_track(size,kind(data),memAddr,id)
                  endif
               else
                  call allocator_errorHandler(alreadyAllocated,var_id=id)
               endif

            end subroutine alloc_backend_P_1d_fixed_i8_indexing_R


!
! wrapper for the allocation routine
! handling arrays with start indexes equal to one
!
            subroutine     alloc_P_normal_1d_fixed_i8_indexing_R(data,size,id)
               integer(kind=8)                       :: size, start_index
               character(len=*), optional            :: id
               character(len=charArrayLength)        :: local_id
               real(kind=kreal), pointer   :: data(:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'unnamed variable'
               endif
               start_index = 1 
               call alloc_backend_P_1d_fixed_i8_indexing_R(data,start_index,size,local_id) 

            end subroutine alloc_P_normal_1d_fixed_i8_indexing_R


!
! wrapper for the allocation routine
! handling arrays with start indexes not equal to one
!
            subroutine     alloc_P_sliced_1d_fixed_i8_indexing_R(data,start_index,end_index,id)
               integer(kind=8)                       :: start_index, end_index
               character(len=*), optional            :: id
               character(len=charArrayLength)        :: local_id
               real(kind=kreal), pointer   :: data(:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'unnamed variable'
               endif

               call alloc_backend_P_1d_fixed_i8_indexing_R(data,start_index,end_index,local_id)

            end subroutine alloc_P_sliced_1d_fixed_i8_indexing_R

      end module allocator_internal_P_1d_fixed_i8_indexing_R

