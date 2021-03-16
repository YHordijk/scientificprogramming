module nice_module

   ! it is good to list only what you will use
   ! this keeps things nicely separated
   use some_module, only: some_routine

   ! always use this
   ! it will be valid module-wide
   implicit none

   public nice_module_init
   public nice_module_finalize
   public nice_module_do_something

   ! everything else is private
   private

   ! if false the module will refuse to be accessed
   logical :: is_initialized = .false.

   ! this array is visible to everything below "contains"
   ! but not visible outside
   ! it is not a good idea to make it visible to outside
   ! rather access through get/set functions that first check whether
   ! the array is allocated
   real(8), allocatable :: some_array(:)

contains

   subroutine nice_module_init(n)
      integer, intent(in) :: n
      allocate(some_array(n))
      is_initialized = .true.
   end subroutine

   subroutine nice_module_finalize()
      ! this routine should reset everything set by this module
      ! and deallocate all allocated memory
      deallocate(some_array)
      is_initialized = .false.
   end subroutine

   subroutine check_if_module_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access nice_module'
         print *, '       but this module is not initialized'
         stop 1
      end if
   end subroutine

   subroutine nice_module_do_something()
      ! first check whether the module is initialized
      call check_if_module_is_initialized()
      ! only now continue with real work
      ! ...
   end subroutine

end module
