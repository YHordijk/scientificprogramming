module mobasis_hartree_fock

!This module contains an in-core Hartree-Fock implementation

        use datatypes
        use diagonalization

        implicit none
        private
        public hartree_fock_driver, read_mo_integrals 

       contains

        subroutine hartree_fock_driver(oneint,twoint)

         implicit none

         type(one_el_t) :: oneint
         complex(8)     :: twoint(:,:,:,:)
         
         print*,"Hello World: let's diagonalize the one-electron matrix testtesttesttesttes!"
 
         call diagonalize_complex_matrix (oneint%h_core)

         print*,"That was interesting, but now I'm leaving"

        end subroutine hartree_fock_driver

        subroutine read_mo_integrals (oneint, twoint)

         implicit none
       
         integer, parameter :: max_orbs=255 ! maximum number or orbitals for which data can be read
         integer, parameter :: funit=15     ! fileunit
         integer            :: norb, nelec, orbsym(max_orbs), isym, iuhf

         type(one_el_t)          :: oneint
         complex(8), allocatable :: twoint(:,:,:,:)
         real(8)                 :: data(2)
         integer                 :: labels(4)
         namelist / FCI / norb, nelec, orbsym, isym, iuhf  

         open (funit,file='FCIDUMP',form='formatted')
         read (funit,FCI)
         allocate (twoint(NORB,NORB,NORB,NORB))
         allocate (oneint%h_core(norb,norb))

         oneint%n_spinor = norb
         labels(1) = -1
         do while (labels(1) /= 0)
            read (funit,*) data,labels
            if (labels(3) /= 0) then 
               twoint(labels(1),labels(2),labels(3),labels(4)) = dcmplx(data(1),data(2)) ! 2e integrals have 4 labels
            else if (labels(2) /= 0) then 
               oneint%h_core(labels(1),labels(2)) = dcmplx(data(1),data(2))  ! 1e integrals have 2 labels
            else if (labels(1) /= 0) then 
               print*, "orbital energy ",labels(1),data(1)  ! orbital energies have 1 label ( print-out)
            else
               oneint%e_core = data(1)
            end if
         end do

        end subroutine read_mo_integrals

        subroutine diagonalize_complex_matrix (the_matrix)

!       Illustrates the diagonalization of a complex Hermitian matrix
!       Note that this routines requires the matrix in a different
!       format, with the real/imaginary dimension as the last instead
!       of the first.

        complex(8), intent(in) :: the_matrix(:,:)
        real(8), allocatable :: reordered_matrix(:,:,:)
        real(8), allocatable :: eigenvalues(:), eigenvectors(:,:,:)
       
        integer ierr, i, n
        
        n = size(the_matrix,1)
        if (size(the_matrix,2) /= n) stop 'Matrix has to be square'

        ! allocate the required memory
        allocate( reordered_matrix(n,n,2) )
        allocate( eigenvectors(n,n,2) )
        allocate( eigenvalues(n) )

        reordered_matrix(:,:,1) = real(the_matrix)
        reordered_matrix(:,:,2) = aimag(the_matrix)
        call qdiag(2,n,reordered_matrix,n,n,eigenvalues,1,eigenvectors,n,n,ierr)

        if (ierr.ne.0) then
            stop 'qdiag failed '
        endif
        print *,"the eigenvalues of the complex hermitian matrix are:"
        do i=1,n
           print *,i,eigenvalues(i)
        enddo

        ! free the memory
        deallocate( reordered_matrix )
        deallocate( eigenvectors )
        deallocate( eigenvalues )

end subroutine diagonalize_complex_matrix

end module mobasis_hartree_fock
