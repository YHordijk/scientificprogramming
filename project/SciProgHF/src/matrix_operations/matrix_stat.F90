module matrix_stat

   implicit none
   
   public add_allocation
   public remove_allocation
   public print_nr_allocations
   public matrix_stat_reset
   
   private
   
   integer :: nr_allocations      = 0
   integer :: nr_allocations_peak = 0

contains

   subroutine matrix_stat_reset()
      nr_allocations      = 0
      nr_allocations_peak = 0
   end subroutine

   subroutine add_allocation()
      nr_allocations = nr_allocations + 1
      if (nr_allocations > nr_allocations_peak) then
         nr_allocations_peak = nr_allocations
      end if
   end subroutine
 
   subroutine remove_allocation()
      nr_allocations = nr_allocations - 1
   end subroutine
 
   subroutine print_nr_allocations(text)
      character(*), intent(in) :: text
      write(*, *) 'debug nr allocations      ', nr_allocations,      text
      write(*, *) 'debug nr allocations peak ', nr_allocations_peak, text
   end subroutine

end module
