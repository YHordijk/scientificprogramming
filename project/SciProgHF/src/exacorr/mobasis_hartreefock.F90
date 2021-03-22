module mobasis_hartree_fock

!This module contains an in-core Hartree-Fock implementation

    use exacorr_datatypes

    implicit none
    private
    public hartree_fock_driver


    type :: settings 
        logical :: ReadFromFile = .false.
        logical :: UHF = .false. !whether to do UHF
        integer :: Nel = -1 !number of electrons
    end type settings

    contains

    !main subroutine which is called from dirac
    subroutine hartree_fock_driver(oneint,twoint)

        implicit none

        type(one_el_t) :: oneint
        complex(8)     :: twoint(:,:,:,:)
        integer        :: nspin, i, j, p, q, r, s, maxiter=20
        real(8), allocatable     :: D(:,:,:), F(:,:,:)
        real(8), allocatable :: eigenvalues(:), eigenvectors(:,:,:)
        real(8) :: thresh = 1e-9, eps1, eps2
        logical :: converged = .false.
        character(256) :: setpath
        type(settings) :: set

        !path to settings file
        setpath = "/mnt/d/Users/Yuman/Desktop/Programmeren/Fortran/scientificprogra&
        &mming/project/SciProgHF/runs/settings"
        call read_settings(setpath, set)


        nspin = oneint%n_spinor

        !if we are doing UHF we double the number of spinors
        if(set%UHF) then
            nspin = nspin*2
            !and half the values in the core hamiltonian and 2e integrals
            oneint%h_core = oneint%h_core/2
            twoint = twoint/2
        endif

        allocate(D(nspin,nspin,2), F(nspin,nspin,2))
        allocate(eigenvalues(nspin), eigenvectors(nspin,nspin,2))

        !initialize D as eye
        D = 0
        do i = 1, nspin
            D(i,i,1) = 1
        enddo

        call print_data(D, eigenvalues)

        i = 1
        do while(i<maxiter .and..not.converged)
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
            call calculate_dens(eigenvectors, D, set)
            call print_data(D, eigenvalues)
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

        print *, "Use UHF: ", set%UHF 
        print *, "Number of electrons: ", set%Nel

    end subroutine read_settings   

    !subroutine that prints values for logging and use in python script later on
    subroutine print_data(D, eigenvalues)
        real(8), intent(in) :: D(:,:,:), eigenvalues(:)
        integer :: nspin, j

        nspin = size(eigenvalues)
        !log energies
        print*, "begin energies"
        do j = 1, nspin
            print*, eigenvalues(j)
        enddo
        print*, "end energies"
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
        real(8), allocatable :: decompg(:,:,:,:,:)
        integer :: n, p, q, r, s, n_occ
        type(settings) :: set
        
        !n is the number of spinors in D, so twice as large as h in UHF
        !n was already set correctly earlier (correct for UHF and RHF)
        n = size(D, 1)
        !get number of occupied spinors
        if(set%UHF) then
            n_occ = set%nel
        else
            n_occ = set%nel / 2
        endif
        
        !initialize the fock matrix as the decomposed core hamiltonian
        !in UHF case we have to set it differently, namely we must make a matrix:

        !F \in \mathbb(R)^{2n} \times \mathbb(R)^{2n}:
        !h(1,1) h(1,1) h(1,2) h(1,2) ... h(1,n) h(1,n)
        !h(1,1) h(1,1) h(1,2) h(1,2) ... h(1,n) h(1,n)
        !h(2,1) h(2,1) h(2,2) h(2,2) ... h(2,n) h(2,n)
        !h(2,1) h(2,1) h(2,2) h(2,2) ... h(2,n) h(2,n)
        ! ...    ...    ...    ...   ...  ...    ...
        !h(n,1) hn2,1) h(n,2) h(n,2) ... hn1,n) h(n,n)
        !h(n,1) hn2,1) h(n,2) h(n,2) ... hn1,n) h(n,n)
        
        !so that each element in h forms a 2x2 subarray in F

        !we use indexing (p-1)/2+1 which gives:
        !p      (p-1)/2+1
        !1      1
        !2      1
        !3      2
        !4      2
        !...
        !n-1    n/2
        !n      n/2

        if(set%UHF) then
            do p = 1, n
                do q = 1, n
                    F(p,q,1) = real(h( (p-1)/2+1, (q-1)/2+1) )
                    F(p,q,2) = aimag(h( (p-1)/2+1, (q-1)/2+1) )
                enddo
            enddo
        else
        !in the case of RHF we can simply copy the core hamiltonian
            F(:,:,1) = real(h)
            F(:,:,2) = aimag(h)
        endif


        !we explicitely decompose the complex matrix g into real parts
        !since it is easier to work with it for now
        allocate(decompg(n,n,n,n,2))
        decompg(:,:,:,:,1) = real(g)
        decompg(:,:,:,:,2) = aimag(g)


        !we iterate p and q over all spinors
        do p = 1, n
            do q = 1, n
                do r = 1, n_occ
                    do s = 1, n_occ
                        !calculate fock matrix for real and imag parts
                        !for UHF decompg is twice as small so we must correct the indices for that

                        if(set%UHF) then
                            F(p,q,1) = F(p,q,1) + decompg((p-1)/2+1,(r-1)/2+1,(q-1)/2+1,(s-1)/2+1,1)*D(r,s,1)
                            F(p,q,2) = F(p,q,2) + decompg((p-1)/2+1,(r-1)/2+1,(q-1)/2+1,(s-1)/2+1,2)*D(r,s,2)
                        else    
                            F(p,q,1) = F(p,q,1) + decompg(p,r,q,s,1)*D(r,s,1)
                            F(p,q,2) = F(p,q,2) + decompg(p,r,q,s,2)*D(r,s,2)
                        endif
                    enddo
                enddo
            enddo
        enddo

        deallocate(decompg)
    end subroutine calculate_Fock 

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

