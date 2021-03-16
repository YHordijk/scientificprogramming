module debug_mem_alloc

   implicit none

   public allocate_deallocate_gb

   private

   type gb
      real(8), pointer :: array(:)
   end type

contains

   subroutine allocate_deallocate_gb(nr_gb)

     integer, intent(in) :: nr_gb
     type(gb)            :: t(nr_gb)
     integer             :: i

     do i = 1, nr_gb
         allocate(t(i)%array(131*1000000))
         t(i)%array = 137.0d0
     end do
     print *, 'successfully allocated', nr_gb, 'GB'

     do i = 1, nr_gb
         deallocate(t(i)%array)
     end do
     print *, 'successfully deallocated', nr_gb, 'GB'

   end subroutine

end module
