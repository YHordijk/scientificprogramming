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

module interface_grid

   implicit none

   public num_grid_write
   public num_grid_read

   private

contains

   subroutine num_grid_write(          &
                             x,        &
                             y,        &
                             z,        &
                             w,        &
                             lunit,    &
                             nr_points &
                            )

!     --------------------------------------------------------------------------
      real(8), intent(in) :: x(*)
      real(8), intent(in) :: y(*)
      real(8), intent(in) :: z(*)
      real(8), intent(in) :: w(*)
      integer, intent(in) :: lunit
      integer, intent(in) :: nr_points
!     --------------------------------------------------------------------------
      integer             :: i
!     --------------------------------------------------------------------------

      write(lunit, *) nr_points
      if (nr_points > 0) then
         do i = 1, nr_points
            write(lunit, '(4d20.12)') x(i), y(i), z(i), w(i)
         enddo
      end if

   end subroutine

   subroutine num_grid_read(          &
                            x,        &
                            y,        &
                            z,        &
                            w,        &
                            lunit,    &
                            nr_points &
                           )

!     --------------------------------------------------------------------------
      real(8), intent(out) :: x(*)
      real(8), intent(out) :: y(*)
      real(8), intent(out) :: z(*)
      real(8), intent(out) :: w(*)
      integer, intent(in)  :: lunit
      integer, intent(in)  :: nr_points
!     --------------------------------------------------------------------------
      integer              :: i
!     --------------------------------------------------------------------------

      do i = 1, nr_points
         read(lunit, *) x(i), y(i), z(i), w(i)
      enddo

   end subroutine

end module
