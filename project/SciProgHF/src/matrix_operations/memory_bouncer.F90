module memory_bouncer

   use matrix_defop_old
   use matrix_stat

   implicit none

   public test_available_memory

   private

   integer, parameter :: max_nr_matrices = 500

contains

   subroutine test_available_memory(nr_matrices, &
                                    nr_ao,       &
                                    nz)

!     --------------------------------------------------------------------------
      integer, intent(in) :: nr_matrices
      integer, intent(in) :: nr_ao
      integer, intent(in) :: nz
!     --------------------------------------------------------------------------
      type(matrix)        :: M(nr_matrices)
      integer             :: i, mb
!     --------------------------------------------------------------------------

      mb = nr_matrices*nz*nr_ao*nr_ao/128000
      write(*, *) 'this calculation will consume approximately Mb = ', mb

      write(*, *) 'will now attempt to allocate the memory'
      do i = 1, nr_matrices
         call init_mat(M(i), nr_ao, nr_ao)
      end do
      write(*, *) 'successfully allocated the memory'

      do i = 1, nr_matrices
         M(i) = 0
      end do
      write(*, *) 'successfully deallocated the memory'
      call matrix_stat_reset()

   end subroutine

end module
