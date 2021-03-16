module conf_parameters

   implicit none

   public allocate_conf_array
   public destroy_conf_array

   private

   real(8), allocatable, target :: conf_array(:)
!  integer, allocatable, target :: orbital_rotation_indices_pp(:)
!  integer, allocatable, target :: orbital_rotation_indices_pn(:)

   logical :: is_initialized = .false.

contains

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'ERROR: you try to access conf_parameters'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   subroutine allocate_conf_array(n)
      integer, intent(in) :: n

      ! routine is called often ...
      if(allocated(conf_array)) return

      allocate(conf_array(n))
      is_initialized = .true.
   end subroutine

   subroutine destroy_conf_array()
      if(allocated(conf_array)) deallocate(conf_array)
      is_initialized = .false.
   end subroutine

end module
