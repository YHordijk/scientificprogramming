module mobasis_hartree_fock

!This module contains an in-core Hartree-Fock implementation

        use exacorr_datatypes

        implicit none
        private
        public hartree_fock_driver

       contains

        subroutine hartree_fock_driver(oneint,twoint)

         implicit none

         type(one_el_t) :: oneint
         complex(8)     :: twoint(:,:,:,:)
         
         print*,"Hello World: let's diagonalize the one-electron matrix !"
 
         call diagonalize_complex_matrix (oneint%h_core)

         print*,"That was interesting, but now I'm leaving"

        end subroutine hartree_fock_driver

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
        if (size(the_matrix,2) /= n) call quit ('Matrix has to be square')

        ! allocate the required memory
        allocate( reordered_matrix(n,n,2) )
        allocate( eigenvectors(n,n,2) )
        allocate( eigenvalues(n) )

        reordered_matrix(:,:,1) = real(the_matrix)
        reordered_matrix(:,:,2) = aimag(the_matrix)
        call qdiag90(2,n,reordered_matrix,n,n,eigenvalues,1,eigenvectors,n,n,ierr)

        if (ierr.ne.0) then
            call quit('qdiag90 failed ')
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

