module mobasis_hartree_fock

    !This module contains an in-core Hartree-Fock implementation

    !STARTING NOW FROM UHF
    
        use exacorr_datatypes
    
        implicit none
        private
        public hartree_fock_driver
    
    
        type :: settings 
            logical :: ReadFromFile = .false.
            logical :: UHF = .false. !whether to do UHF
            integer :: Nel = -1 !number of electrons
            logical :: usemixing = .true.
            real(8) :: mixing = 0.1
            integer :: maxiter = 20
        end type settings
    
        contains
    
        !main subroutine which is called from dirac
        subroutine hartree_fock_driver(oneint,twoint)
    
            implicit none
    
            type(one_el_t) :: oneint
            complex(8)     :: twoint(:,:,:,:)
            integer        :: nspin, i, j, p, q, r, s
            real(8), allocatable     :: D(:,:,:), F(:,:,:), Dn(:,:,:)
            real(8), allocatable :: eigenvalues(:), eigenvectors(:,:,:)
            real(8) :: thresh = 1e-9, eps1, eps2, energy
            logical :: converged = .false.
            character(256) :: setpath, cwd
            type(settings) :: set
    
            !path to settings file
            call getcwd(cwd)
            setpath = trim(cwd) // '/settings'
            call read_settings(setpath, set)
    
            !if we are doing UHF we double the number of spinorbitals
            nspin = oneint%n_spinor / 2
            if(set%UHF) nspin = nspin*2

            print*, nspin
    
            allocate(D(nspin,nspin,2), Dn(nspin,nspin,2), F(nspin,nspin,2))
            allocate(eigenvalues(nspin), eigenvectors(nspin,nspin,2))
    
            !initialize D as eye
            D = 0
            do i = 1, nspin
                D(i,i,1) = 1
            enddo
            call calculate_energy(oneint%h_core, oneint%e_core, twoint, D, energy, set)
    
            call print_data(D, eigenvalues, energy)
    
            i = 1
            do while(i<set%maxiter .and..not.converged)
                eps1 = eps2
                !get fock matrix first
                call calculate_Fock(oneint%h_core, D, twoint, F, set)
                !diagonalize F to get eigenvectors
                call diagonalize_complex_matrix (F, eigenvalues, eigenvectors)
                !calculate convergence criterium using density from previous iteration
                call calculate_eps(F, D, eps2)
                !check for convergence:
                if(abs(eps1-eps2)<thresh .and. i>1) then
                    ! converged = .true.
                endif
                !if not converged calculate the new density
                call calculate_dens(eigenvectors, Dn, set)
                call calculate_energy(oneint%h_core, oneint%e_core, twoint, D, energy, set)
                
                !if we are using mixing then we are going to mix the density matrices
                if(set%usemixing) then
                    D = (1 - set%mixing) * D + set%mixing * Dn
                else
                    D = Dn
                endif
                call print_data(D, eigenvalues, energy)
                i = i + 1
            enddo
    
        end subroutine hartree_fock_driver
    
    
        subroutine read_settings(f, set)
            character(256), intent(in) :: f
            type(settings), intent(out) :: set
            integer :: iostat = 0
            logical :: fexist
            
            character(64) :: var 
    
            inquire(file=f, exist=fexist)
            
            if (.not.fexist) return
    
            print*, "Found settings file"
            set%ReadFromFile = .true.
    
            open(20, file=f)
    
            !settings file should look like follows:
            !UHF 
            !logical
            !Nel 
            !integer
    
            read(20,*)
            read(20,*) set%UHF
            read(20,*)
            read(20,*) set%Nel
            read(20,*)
            read(20,*) set%usemixing
            read(20,*)
            read(20,*) set%mixing
            read(20,*)
            read(20,*) set%maxiter
    
            print *, "Use UHF: ", set%UHF 
            print *, "Number of electrons: ", set%Nel
            print *, "Use mixing: ", set%usemixing
            print *, "Mixing strenght: ", set%mixing
            print *, "Max iteration: ", set%maxiter
    
        end subroutine read_settings   
    
    
        !subroutine that prints values for logging and use in python script later on
        subroutine print_data(D, eigenvalues, energy)
            real(8), intent(in) :: D(:,:,:), eigenvalues(:), energy
            integer :: nspin, j
    
            nspin = size(eigenvalues)
            !log energies
            print*, "Energy: ", energy
            print*, "begin eigenvalues"
            do j = 1, nspin
                print*, eigenvalues(j)
            enddo
            print*, "end eigenvalues"
            !log the densities for visualization using python and debugging
            !real part of dens
            print*, "begin density"
            do j = 1, nspin
                print*, D(j,:,1)
            enddo
            print*, "end density"
            !imag part of dens
            print*, "begin density"
            do j = 1, nspin
                print*, D(j,:,2)
            enddo
            print*, "end density"
        end subroutine print_data
    
    
        !subroutine for calculating convergence as described in project description
        subroutine calculate_eps(F, D, eps)
            real(8), intent(in) :: D(:,:,:), F(:,:,:)
            real(8), allocatable :: comm2(:,:,:)
            integer :: n, p, q
            real(8), intent(out) :: eps
    
    
            n = size(F, 1)
            allocate(comm2(n,n,2))
    
            !calculate real and imaginary parts of the commutator
            comm2(:,:,1) = matmul(F(:,:,1),D(:,:,1)) - matmul(D(:,:,1),F(:,:,1))
            comm2(:,:,2) = matmul(F(:,:,2),D(:,:,2)) - matmul(D(:,:,2),F(:,:,2))
    
            eps = 0.0
            do p = 1, n 
                do q = 1, n
                    !square commutator
                    eps = comm2(p,q,1)*comm2(p,q,1) + comm2(p,q,2)*comm2(p,q,2)
                enddo
            enddo
    
            eps = sqrt(eps)
        end subroutine calculate_eps
    
    
        !subroutine for calculating density matrix
        subroutine calculate_dens(C, D, set)
            real(8), intent(in) :: C(:,:,:)
            real(8), intent(out) :: D(:,:,:)
            integer :: p, q, n, i, n_occ
            type(settings) :: set  
    
            !here n is the number of spinors, this is already set previously for UHF
            !we need to change the number of occupied spinors however.
            !In UHF case we have that n_occ = n_elec
            !In RHF case we have that n_occ = floor(n_elec/2)
            n = size(D, 1)
    
            if(set%UHF) then
                n_occ = set%nel
            else
                n_occ = set%nel / 2
            endif

            !iterate over all spinors p and q and over occupied spinors i
            do p = 1, n
                do q = 1, n
                    do i = 1, n_occ
                        !calculate real part of density
                        D(p,q,1) = C(p,i,1) * C(q,i,1) + C(p,i,2) * C(q,i,2)
                        !and imaginary part of density
                        D(p,q,2) = C(p,i,1) * C(q,i,2) - C(p,i,2) * C(q,i,1)
                    enddo
                enddo
            enddo
    
        end subroutine calculate_dens
    
    
        !subroutine for calculating fock matrix
        subroutine calculate_Fock(h, D, g, F, set)
            !we have that size(D,1) == size(F,1) and size(h,1) == size(g,1)
            complex(8), intent(in) :: h(:,:), g(:,:,:,:)
            real(8), intent(in) :: D(:,:,:)
            real(8), intent(out) :: F(:,:,:)
            integer :: n, p, q, r, s, n_occ
            type(settings) :: set !settings object containing a few parameters for the calculations
            
            !n is the number of spinors in D, so twice as large as h in UHF
            !n was already set correctly earlier (correct for UHF and RHF)
            n = size(D, 1)
            !get number of occupied spinors
            if(set%UHF) then
                n_occ = set%nel
            else
                n_occ = set%nel / 2
            endif
        
            !copy the core_hamiltonian into the Fock matrix
            do p = 1, n
                do q = 1, n
                    F(p,q,1) = real(h(p,q))
                    F(p,q,2) = aimag(h(p,q))
                enddo
            enddo

            !we iterate p and q over all spinorbitals
            do p = 1, n
                do q = 1, n
                    !and r and s over occupied spinorbitals
                    do r = 1, n_occ
                        do s = 1, n_occ
                            !calculate fock matrix for real and imag parts
                            !for UHF decompg is twice as small so we must correct the indices for that
                            print *, set%UHF, p, q, r, s

                            F(p,q,1) = F(p,q,1) + real(g(p,r,q,s))*D(r,s,1) - aimag(g(p,r,q,s))*D(r,s,2)
                            F(p,q,2) = F(p,q,2) + aimag(g(p,r,q,s))*D(r,s,2) + real(g(p,r,q,s))*D(r,s,2)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine calculate_Fock 
    
    
        subroutine calculate_energy(h, Enuc, g, D, energy, set)
            complex(8), intent(in) :: h(:,:), g(:,:,:,:)
            real(8), intent(in) :: D(:,:,:), Enuc
            real(8), intent(out) :: energy
            real(8) :: energy1, energy2
            complex(8), allocatable :: CD(:,:)
            type(settings) :: set
    
            integer :: p, q, r, s, n, n_occ
    
            n = size(D, 1)


            if(set%UHF) then
                n_occ = set%nel
            else
                n_occ = set%nel / 2
            endif
    
            allocate(CD(n,n))
            CD = CMPLX(D(:,:,1), D(:,:,2))

            !calcualte the first term
            energy1 = 0
            do p = 1, n
                do q = 1, n
                    energy1 = energy1 + h(p,q) * CD(p,q)
                enddo
            enddo
    
    
            !calculate second term
            energy2 = 0
            do p = 1, n
                do q = 1, n
                    do r = 1, n_occ
                        do s = 1, n_occ
                            ! if(set%UHF) then
                            !     energy2 = energy2 + g((p-1)/2+1,(r-1)/2+1,(q-1)/2+1,(s-1)/2+1) * CD(p,q) * CD(r,s)  
                            ! else    
                            !     energy2 = energy2 + g(p,r,q,s) * CD(p,q) * CD(r,s) 
                            ! endif
                            energy2 = energy2 + g(p,r,q,s) * CD(p,q) * CD(r,s) 
                        enddo
                    enddo
                enddo
            enddo
    
            energy = energy1 + 0.5d0 * energy2 + Enuc
    
        end subroutine calculate_energy
    
    
        !subroutine for diagonalizing a complex matrix
        subroutine diagonalize_complex_matrix (the_matrix, eigenvalues, eigenvectors)
    
        !       Illustrates the diagonalization of a complex Hermitian matrix
        !       Note that this routines requires the matrix in a different
        !       format, with the real/imaginary dimension as the last instead
        !       of the first.
    
            real(8), intent(in) :: the_matrix(:,:,:)
            real(8), allocatable :: reordered_matrix(:,:,:)
            real(8), intent(out) :: eigenvalues(:), eigenvectors(:,:,:)
    
            integer ierr, i, n
            
            n = size(the_matrix,1)
            if (size(the_matrix,2) /= n) call quit ('Matrix has to be square')
    
            call qdiag90(2,n,the_matrix,n,n,eigenvalues,1,eigenvectors,n,n,ierr)
    
            if (ierr.ne.0) then
            call quit('qdiag90 failed ')
            endif
    
        end subroutine diagonalize_complex_matrix
    
    end module mobasis_hartree_fock