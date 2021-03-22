program HartreeFock

! Stand-alone version of skeleton Hartree-Fock program
! Only needs the compilation and linking of the driver module 
! and the auxilliary files in this directory.

         use datatypes
         use mobasis_hartree_fock
         
         implicit none
         type(one_el_t)          :: oneint
         complex(8), allocatable :: twoint(:,:,:,:)

         call read_mo_integrals (oneint, twoint)
         call hartree_fock_driver(oneint,twoint)

end program HartreeFock
