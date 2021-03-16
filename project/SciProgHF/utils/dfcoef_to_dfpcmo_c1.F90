
! author: Radovan Bast (based on Ulf's dfcoef_read.F)

!         convert DFCOEF to DFPCMO
!         only C1
!         only 4-component DFCOEF
!         use at own risk

program dfcoef_read

   implicit none

   integer, parameter   :: io = 1
   integer              :: nr_arg
   character(80)        :: arg
   character(74)        :: title
   real(8), allocatable :: coef(:)
   real(8), allocatable :: eigs(:)
   real(8)              :: total_energy
   integer              :: nr_mo_negative
   integer              :: nr_mo_positive
   integer              :: nr_mo
   integer              :: nr_ao
   integer              :: idummy
   integer              :: i

   nr_arg = iargc()

   if (nr_arg /= 1) then
      print *, 'usage: dfcoef_to_dfpcmo_c1 DFCOEF'
      stop 1
   end if

   call getarg(1, arg)

   open(io, file=arg, status='unknown', form='unformatted')
   rewind io

   read(io) title, idummy, nr_mo_negative, nr_mo_positive, nr_ao, total_energy

   nr_mo = nr_mo_negative + nr_mo_positive

   allocate(coef(nr_ao*nr_mo*4))
   allocate(eigs(nr_mo))

   read(io) coef
   read(io) eigs

   write(*, '(a74)'    ) title
   write(*, '(i4,3i10)') 1, nr_mo_negative, nr_mo_positive, nr_ao
   write(*, '(d24.16)' ) total_energy
   write(*, '(6f22.16)') (coef(i), i = 1, size(coef))
   write(*, '(6d22.12)') (eigs(i), i = 1, size(eigs))
   write(*, '(66i2)'   ) (0,       i = 1, size(eigs))

   deallocate(coef)
   deallocate(eigs)

   close(io, status='keep')

end program
