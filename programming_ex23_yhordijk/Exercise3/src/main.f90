program eigenproblem

    use Diagonalization
    use arrayfunc

    ! Solve the eigenproblem for a symmetric matrix
    integer :: natoms, con1, con2, i
    real*8, allocatable :: H(:,:) !huckel matrix
    real*8, allocatable :: eigvec(:,:), eigval(:), homo(:), lumo(:) !eigenvectors and values
    real :: alpha, beta

    character(100) :: fin

    fin = "src\input\benzene.in"

    print *, ' ' 
    print *, "Reading: ", fin
    print *, ' ' 
    open(unit=15, file=fin)
    alpha = -11.
    beta = -0.7

    read(15, *) natoms
    allocate(H(natoms, natoms))
    allocate(eigvec(natoms, natoms))
    allocate(eigval(natoms))
    allocate(Homo(natoms))
    allocate(Lumo(natoms))

    !initialize huckel matrix as zero except diagonals as alpha
    do i = 1, natoms
        do j = 1, natoms
            if (i .eq. j) then
                H(i,j) = alpha
            else
                H(i,j) = 0
            end if 
        end do 
    end do 

    !continue reading the file and set all given connections as beta
    do
        read(15, *) con1, con2
        !check if we hit eof
        if(con1 < 0) exit

        H(con1, con2) = beta
        H(con2, con1) = beta
    end do

    !print out the huckel matrix
    print *, "Generated Huckel matrix:"
    call arrayprint(H)

    call diagonalize(H, eigvec, eigval)

    print *, ' ' 
    print *, "Diagonalization complete"
    print *, "  Eigenvalues (eV):"
    call vectorprint(eigval)
    print *, "  Eigenvectors (as column vectors) with elements as coefficients of p-orbital basis:"
    call arrayprint(eigvec)
     
    print *, ' ' 
    HOMO = getcolumn(eigvec, int(natoms/2))
    print *, "HOMO with energy", eigval(int(natoms/2)), "eV"
    print *, "  Coefficients in p-orbital basis:"
    call vectorprint(HOMO)

    LUMO = getcolumn(eigvec, int(natoms/2)+1)
    print *, "LUMO with energy", eigval(int(natoms/2)+1), "eV"
    print *, "  Coefficients in p-orbital basis:"
    call vectorprint(LUMO)

end program eigenproblem

