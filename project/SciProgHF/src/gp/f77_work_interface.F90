module f77_work_interface

!  if you can, you should introduce allocations using allocate/deallocate or alloc/dealloc
!  but sometimes you have no other choice than to use work (because you don't know how much a routine will consume)
!  this module is for the situation where you want to use the work array but it is not available 
!  in your routine/module

!  example usage:

!  use f77_work_interface
!  ...
!  call plug_f77_work(work, lwork) !here work and lwork are available, lwork is used, not set
!  call routine()
!  ...
!  call unplug_f77_work()
!
!     subroutine routine() !has no work array
!        call another_routine()
!     end subroutine
!
!     subroutine another_routine() !has no work array
!        use f77_work_interface
!        ... 
!        real(8), pointer :: work(:)
!        integer          :: lwork
!        ...
!        call select_f77_work(work, lwork) !lwork is set, not used
!        call old_routine(var1, var2, work, lwork)
!        call unselect_f77_work(work)
!     end subroutine

   implicit none

   public plug_f77_work
   public unplug_f77_work
   public select_f77_work
   public unselect_f77_work

   private

   real(8), pointer :: f77_work(:)
   integer, public  :: len_f77_work = 0
   integer          :: work_count   = 0

contains

   subroutine plug_f77_work(work, lwork)

!     --------------------------------------------------------------------------
      real(8), target     :: work(:)
      integer, intent(in) :: lwork
!     --------------------------------------------------------------------------

      len_f77_work = lwork
      f77_work     => work

   end subroutine

   subroutine unplug_f77_work()

      nullify(f77_work)
      len_f77_work = 0

   end subroutine

   subroutine select_f77_work(work, lwork)

!     --------------------------------------------------------------------------
      real(8), pointer     :: work(:)
      integer, intent(out) :: lwork
!     --------------------------------------------------------------------------

      if (work_count > 1) then
         call quit('work array is already associated')
      end if

      work => f77_work
      work_count = work_count + 1
      lwork = len_f77_work

   end subroutine

   subroutine unselect_f77_work(work)

!     --------------------------------------------------------------------------
      real(8), pointer :: work(:)
!     --------------------------------------------------------------------------

      nullify(work)
      work_count = work_count - 1

   end subroutine

end module
