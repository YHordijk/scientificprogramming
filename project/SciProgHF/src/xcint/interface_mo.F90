
module interface_mo

   use file_units

   implicit none

   public interface_mo_read
   public interface_mo_clear

   logical, public :: interface_mo_set

   integer, public :: nr_electrons_inactive
   integer, public :: nr_electrons_active
   integer, public :: nr_electrons
   integer, public :: nr_mo
   integer, public :: nr_mo_gerade
   integer, public :: nr_mo_gerade_negative_secondary
   integer, public :: nr_mo_gerade_positive_inactive
   integer, public :: nr_mo_gerade_positive_active
   integer, public :: nr_mo_gerade_positive_secondary
   integer, public :: nr_mo_ungerade
   integer, public :: nr_mo_ungerade_negative_secondary
   integer, public :: nr_mo_ungerade_positive_inactive
   integer, public :: nr_mo_ungerade_positive_active
   integer, public :: nr_mo_ungerade_positive_secondary
   integer, public :: nr_open_shells
   logical, public :: mo_eigenvalues_available
   logical, public :: mo_coef_available

   real(8), allocatable, public :: mo_eigenvalues(:)
   real(8), allocatable, public :: mo_coef(:)
   real(8), allocatable, public :: df_open(:)

   private

contains

   subroutine interface_mo_clear()

      interface_mo_set                  = .false.

      nr_electrons_inactive             = 0
      nr_electrons_active               = 0
      nr_electrons                      = 0
      nr_mo                             = 0
      nr_mo_gerade                      = 0
      nr_mo_gerade_negative_secondary   = 0
      nr_mo_gerade_positive_inactive    = 0
      nr_mo_gerade_positive_active      = 0
      nr_mo_gerade_positive_secondary   = 0
      nr_mo_ungerade                    = 0
      nr_mo_ungerade_negative_secondary = 0
      nr_mo_ungerade_positive_inactive  = 0
      nr_mo_ungerade_positive_active    = 0
      nr_mo_ungerade_positive_secondary = 0
      nr_open_shells                    = 0
      mo_eigenvalues_available          = .false.
      mo_coef_available                 = .false.

      if (allocated(df_open        )) deallocate(df_open        )
      if (allocated(mo_coef        )) deallocate(mo_coef        )
      if (allocated(mo_eigenvalues )) deallocate(mo_eigenvalues )

   end subroutine

   subroutine interface_mo_read()

      use xc_mpi

!     --------------------------------------------------------------------------
      character(120) :: line
      integer        :: i, idummy, n, size_mo_coef
!     --------------------------------------------------------------------------

      call interface_mo_clear()

      if (xc_mpi_is_master()) then
         inquire(file = interface_file_name_mo, exist = interface_file_exists)
         if (.not. interface_file_exists) then
            print *, 'interface_mo file not found'
            stop
         end if
         open(interface_file_unit,             &
              file   = interface_file_name_mo, &
              status = 'unknown',              &
              form   = 'formatted',            &
              access = 'sequential')
         rewind(interface_file_unit)

         do while (.true.)
            read(interface_file_unit, '(a120)', end=1) line

            if (word_contains(line, 'nr_electrons_inactive_total')) then
               read(interface_file_unit, *) nr_electrons_inactive
            end if

            if (word_contains(line, 'nr_electrons_active_total')) then
               read(interface_file_unit, *) nr_electrons_active
            end if

            if (word_contains(line, 'nr_mo_gerade_negative_secondary')) then
               read(interface_file_unit, *) nr_mo_gerade_negative_secondary
            end if

            if (word_contains(line, 'nr_mo_gerade_positive_inactive')) then
               read(interface_file_unit, *) nr_mo_gerade_positive_inactive
            end if

            if (word_contains(line, 'nr_mo_gerade_positive_active')) then
               read(interface_file_unit, *) nr_mo_gerade_positive_active
            end if

            if (word_contains(line, 'nr_mo_gerade_positive_secondary')) then
               read(interface_file_unit, *) nr_mo_gerade_positive_secondary
            end if

            if (word_contains(line, 'nr_mo_ungerade_negative_secondary')) then
               read(interface_file_unit, *) nr_mo_ungerade_negative_secondary
            end if

            if (word_contains(line, 'nr_mo_ungerade_positive_inactive')) then
               read(interface_file_unit, *) nr_mo_ungerade_positive_inactive
            end if

            if (word_contains(line, 'nr_mo_ungerade_positive_active')) then
               read(interface_file_unit, *) nr_mo_ungerade_positive_active
            end if

            if (word_contains(line, 'nr_mo_ungerade_positive_secondary')) then
               read(interface_file_unit, *) nr_mo_ungerade_positive_secondary
            end if

            if (word_contains(line, 'nr_open_shells')) then
               read(interface_file_unit, *) nr_open_shells
               if (nr_open_shells > 0) then
                  allocate(df_open(nr_open_shells))
                  df_open = 0.0d0
                  do i = 1, nr_open_shells
                     read(interface_file_unit, *) idummy, df_open(i)
                  end do
               end if
            end if

         end do
1        continue

         nr_electrons = nr_electrons_inactive + nr_electrons_active
         nr_mo_gerade = nr_mo_gerade_negative_secondary &
                      + nr_mo_gerade_positive_inactive  &
                      + nr_mo_gerade_positive_active    &
                      + nr_mo_gerade_positive_secondary
         nr_mo_ungerade = nr_mo_ungerade_negative_secondary &
                        + nr_mo_ungerade_positive_inactive  &
                        + nr_mo_ungerade_positive_active    &
                        + nr_mo_ungerade_positive_secondary
         nr_mo = nr_mo_gerade + nr_mo_ungerade

         rewind(interface_file_unit)
         do while (.true.)

            read(interface_file_unit, '(a120)', end=2) line

            if (word_contains(line, 'mo_eigenvalues')) then
               allocate(mo_eigenvalues(nr_mo))
               do i = 1, nr_mo
                  read(interface_file_unit, *) idummy, mo_eigenvalues(i)
               end do
               mo_eigenvalues_available = .true.
            end if

         end do
2        continue
         close(interface_file_unit, status = 'keep')

         inquire(file = interface_file_name_mo_coef, exist = interface_file_exists)
         if (interface_file_exists) then
            open(interface_file_unit,                  &
                 file   = interface_file_name_mo_coef, &
                 status = 'unknown',         &
                 form   = 'formatted',       &
                 access = 'sequential')
            rewind(interface_file_unit)
            read(interface_file_unit, *) n
            allocate(mo_coef(n))
            do i = 1, n
               read(interface_file_unit, *) mo_coef(i)
            end do
            close(interface_file_unit, status = 'keep')
            mo_coef_available = .true.
         end if

         interface_mo_set = .true.
      end if

#ifdef VAR_MPI
      if (xc_mpi_get_nr_proc() > 1) then

         call xc_mpi_bcast(interface_mo_set                  )

         call xc_mpi_bcast(nr_electrons_inactive             )
         call xc_mpi_bcast(nr_electrons_active               )
         call xc_mpi_bcast(nr_electrons                      )
         call xc_mpi_bcast(nr_mo                             )
         call xc_mpi_bcast(nr_mo_gerade                      )
         call xc_mpi_bcast(nr_mo_gerade_negative_secondary   )
         call xc_mpi_bcast(nr_mo_gerade_positive_inactive    )
         call xc_mpi_bcast(nr_mo_gerade_positive_active    )
         call xc_mpi_bcast(nr_mo_gerade_positive_secondary   )
         call xc_mpi_bcast(nr_mo_ungerade                    )
         call xc_mpi_bcast(nr_mo_ungerade_negative_secondary )
         call xc_mpi_bcast(nr_mo_ungerade_positive_inactive  )
         call xc_mpi_bcast(nr_mo_ungerade_positive_active  )
         call xc_mpi_bcast(nr_mo_ungerade_positive_secondary )
         call xc_mpi_bcast(nr_open_shells                        )
         call xc_mpi_bcast(mo_eigenvalues_available          )
         call xc_mpi_bcast(mo_coef_available                 )

         if (nr_open_shells > 0) then
            if (.not. xc_mpi_is_master()) then
               allocate(df_open(nr_open_shells))
            end if
            call xc_mpi_bcast(df_open)
         end if

         if (mo_coef_available) then
            if (xc_mpi_is_master()) then
               size_mo_coef = size(mo_coef)
            end if
            call xc_mpi_bcast(size_mo_coef)
            if (.not. xc_mpi_is_master()) then
               allocate(mo_coef(size_mo_coef))
            end if
            call xc_mpi_bcast(mo_coef)
         end if

         if (mo_eigenvalues_available) then
            if (.not. xc_mpi_is_master()) then
               allocate(mo_eigenvalues(nr_mo))
            end if
            call xc_mpi_bcast(mo_eigenvalues)
         end if
      end if
#endif

   end subroutine

   function word_contains(word, substring)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: word
      character(*), intent(in) :: substring
!     --------------------------------------------------------------------------
      logical                  :: word_contains
!     --------------------------------------------------------------------------

      if (index(word, substring) > 0) then
         word_contains = .true.
      else
         word_contains = .false.
      end if

   end function

end module
