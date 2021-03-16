
module interface_file_io

   implicit none

   public matrix_to_file
   public matrix_from_file

   private

   integer, parameter :: file_unit = 137
   logical            :: file_exists

contains

   subroutine matrix_to_file(file_name, l, matrix)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: file_name
      integer,      intent(in) :: l
      real(8)                  :: matrix(*)
!     --------------------------------------------------------------------------
      integer                  :: i
!     --------------------------------------------------------------------------

      inquire(file = file_name, exist = file_exists)
      if (file_exists) then
         open(file_unit,             &
              file   = file_name,    &
              status = 'unknown',    &
              form   = 'formatted',  &
              access = 'sequential')
         close(file_unit, status = 'delete')
      end if
      open(file_unit,             &
           file   = file_name,    &
           status = 'new',        &
           form   = 'formatted',  &
           access = 'sequential')
      rewind(file_unit)

      write(file_unit, *) l
      do i = 1, l
         write(file_unit, *) matrix(i)
      end do

      close(file_unit, status = 'keep')

   end subroutine

   subroutine matrix_from_file(file_name, l, matrix)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: file_name
      integer,      intent(in) :: l
      real(8)                  :: matrix(:, :)
!     --------------------------------------------------------------------------
      integer                  :: i, j, n
!     --------------------------------------------------------------------------

      inquire(file = file_name, exist = file_exists)
      if (file_exists) then
         open(file_unit,             &
              file   = file_name,    &
              status = 'unknown',    &
              form   = 'formatted',  &
              access = 'sequential')
      else
         print *, 'error: file '//file_name//' not found'
         stop
      end if
      rewind(file_unit)

      read(file_unit, *) n
      if (n /= l*l) then
         print *, 'error: matrix dimensions do not match in file '//file_name
         stop
      end if

      do i = 1, l
         do j = 1, l
            read(file_unit, *) matrix(j, i)
         end do
      end do

      close(file_unit, status = 'keep')

   end subroutine

end module
