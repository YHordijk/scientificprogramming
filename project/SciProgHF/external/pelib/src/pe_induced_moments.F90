!
!   Polarizable Embedding (PE) library
!   Copyright (C) 2013, 2014 The PE library developers. See the CONTRIBUTORS file
!                            in the top-level directory of this distribution.
!
!   This file is part of the PE library.
!
!   The PE library is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as
!   published by the Free Software Foundation, either version 3 of the
!   License, or (at your option) any later version.
!
!   The PE library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the PE library. If not, see <http://www.gnu.org/licenses/>.
!
!   Contact information:
!
!   Jogvan Magnus Haugaard Olsen
!   E-mail: foeroyingur@gmail.com
!
module pe_induced_moments

    use pe_precision

    implicit none

    private

    public :: induced_moments, response_matrix

contains

!------------------------------------------------------------------------------

subroutine induced_moments(Mkinds, Fks)

    use pe_variables

    real(dp), dimension(:,:), intent(out), optional :: Mkinds
    real(dp), dimension(:,:), intent(in), optional :: Fks

    integer :: i, j, k

    if (any(lpol)) then
        if (.not. present(Mkinds) .and. .not. present(Fks)) then
            stop 'ERROR: missing in-/output pe_variables'
        end if
    else
        stop 'ERROR: pe_polar should be .false. here'
    end if

    if (pe_iter) then
        call iterative_solver(Mkinds, Fks)
    else
        if (myid == master) then
            call direct_solver(Mkinds=Mkinds, Fks=Fks)
        end if
    end if

end subroutine induced_moments

!------------------------------------------------------------------------------

subroutine iterative_solver(Mkinds, Fks)

    use pe_mpi
    use pe_variables
    use pe_utils, only: openfile
    use pe_blas_interfaces, only: nrm2, spmv
    use pe_multipole_tools, only: Tk_tensor, Tk_damped_tensor

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fks

    integer :: lu, iter
    integer :: i, j, k, l, m, n
    logical :: exclude, lexist, restart, converged
    real(dp) :: alpha_i, alpha_j, norm, redthr
    real(dp), parameter :: i3 = 1.0_dp / 3.0_dp
    real(dp), dimension(:), allocatable :: T2, Rij, Ftmp, M1tmp

    if (myid == master) then
        if (pe_redthr .and. fock .and. .not. pe_restart) then
            redthr = 10.0_dp**(- nint(log10(thriter)) - 2 * scf_cycle - redlvl)
            if (redthr * thriter > thriter) then
                write(luout,'(4x,a,es8.1)') 'INFO: using reduced threshold to&
                                            & determine induced dipole moments: ', redthr * thriter
            else
                redthr = 1.0_dp
            end if
        else
            redthr = 1.0_dp
        end if
    end if

    if (myid == master) then
        inquire(file='pe_induced_moments.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif

    if (lexist .and. (fock .or. energy) .and. (scf_cycle > 1 .or. pe_restart) .and. .not. pe_nomb) then
        restart = .true.
    else
        restart = .false.
    end if

    if (restart) then
        if (myid == master) then
            call openfile('pe_induced_moments.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Mkinds
            close(lu)
        end if
    else
        Mkinds = 0.0_dp
    end if

    allocate(T2(6), Rij(3), Ftmp(3), M1tmp(3))
    do n = 1, size(Mkinds, 2)
        if (.not. restart) then
#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + poldists(i-1)
                end do
                call mpi_scatterv(Fks(1,n), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
            else if (myid /= master) then
                call mpi_scatterv(dummy_real, poldists, displs, rmpi, Fks(1,n), poldists(myid), rmpi, master, comm, ierr)
            end if
#endif

            l = 1
            do i = site_start, site_finish
                if (zeroalphas(i)) cycle
                call spmv(P11s(:,i), Fks(l:l+2,n), Mkinds(l:l+2,n), 'L')
                l = l + 3
            end do

#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + poldists(i-1)
                end do
                call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Mkinds(1,n), poldists, displs, rmpi, master, comm, ierr)
            else if (myid /= master) then
                call mpi_gatherv(Mkinds(1,n), poldists(myid), rmpi, dummy_real, poldists, displs, rmpi, master, comm, ierr)
            end if
#endif
        end if

        if (pe_nomb) cycle

#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + poldists(i-1)
            end do
            call mpi_scatterv(Mkinds(1,n), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
        else if (myid /= master) then
            call mpi_scatterv(dummy_real, poldists, displs, rmpi, Mkinds(1,n), poldists(myid), rmpi, master, comm, ierr)
        end if
#endif

        iter = 1
        do
            norm = 0.0_dp
            l = 1
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                if (pe_ind_damp) then
                    alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
                end if
                m = 1
                Ftmp = 0.0_dp
                do j = site_start, site_finish
                    if (zeroalphas(j)) cycle
                    exclude = .false.
                    do k = 1, lexlst
                        if (exclists(k,i) == exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (i == j .or. exclude) then
                        m = m + 3
                        cycle
                    end if
                    Rij = Rs(:,j) - Rs(:,i)
                    if (pe_ind_damp) then
                        alpha_j = (P11s(1,j) + P11s(4,j) + P11s(6,j)) * i3
                        call Tk_damped_tensor(T2, Rij, alpha_i, alpha_j, ind_damp)
                    else
                        call Tk_tensor(T2, Rij)
                    end if
                    call spmv(T2, Mkinds(m:m+2,n), Ftmp, 'L', 1.0_dp, 1.0_dp)
                    m = m + 3
                end do

#if defined(VAR_MPI)
                if (myid == master .and. nprocs > 1) then
                    call mpi_reduce(mpi_in_place, Ftmp(1), 3, rmpi, mpi_sum, master, comm, ierr)
                else if (myid /= master) then
                    call mpi_reduce(Ftmp(1), dummy_real, 3, rmpi, mpi_sum, master, comm, ierr)
                end if
#endif

                if (myid == master) then
                    M1tmp = Mkinds(l:l+2,n)
                    Ftmp = Ftmp + Fks(l:l+2,n)
                    call spmv(P11s(:,i), Ftmp, Mkinds(l:l+2,n), 'L')
                    M1tmp = Mkinds(l:l+2,n) - M1tmp
                    norm = norm + nrm2(M1tmp)
                end if

#if defined(VAR_MPI)
                if (myid == master .and. nprocs > 1) then
                    displs(0) = 0
                    do j = 1, nprocs
                        displs(j) = displs(j-1) + poldists(j-1)
                    end do
                    call mpi_scatterv(Mkinds(1,n), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
                else if (myid /= master) then
                    call mpi_scatterv(dummy_real, poldists, displs, rmpi, Mkinds(1,n), poldists(myid), rmpi, master, comm, ierr)
                end if
#endif
                l = l + 3
            end do

            if (myid == master) then
                if (norm < redthr * thriter) then
                    if (pe_verbose .and. .not. energy) then
                        write(luout, '(6x,a,i2,a)') 'Induced dipole moments converged in ', iter, ' iterations.'
                    end if
                    converged = .true.
                else if (iter > 50) then
                    write(luout, *) 'ERROR: could not converge induced dipole moments'
                    write(luout, *) 'Current norm: ', norm
                    write(luout, *) 'Target: ', redthr * thriter
                    stop 'ERROR: could not converge induced dipole moments'
                else
                    converged = .false.
                    iter = iter + 1
                end if
            end if

#if defined(VAR_MPI)
            if (nprocs > 1) then
                call mpi_bcast(converged, 1, lmpi, master, comm, ierr)
            end if
#endif
            if (converged) exit
        end do
    end do

    if (fock) then
        if (myid == master) then
            call openfile('pe_induced_moments.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Mkinds
            close(lu)
        end if
    end if

    deallocate(T2, Rij, Ftmp, M1tmp)

end subroutine iterative_solver

!------------------------------------------------------------------------------

subroutine direct_solver(Mkinds, Fks)

    use pe_variables
    use pe_utils, only: openfile
    use pe_lapack_interfaces, only: lansp, pptrf, sptrf, pptrs, sptrs, ppcon, spcon

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fks

    logical :: lexist
    integer :: nelems, lu, info
    integer, dimension(:), allocatable :: ipiv
    real(dp) :: anorm, rcond
    real(dp), dimension(:), allocatable :: B
    real(dp), dimension(:,:), allocatable :: F
    real(dp) :: eps_fac

    nelems = 3 * npols

    allocate(B(nelems*(nelems+1)/2))
    allocate(F(nelems,ndens))

    inquire(file='pe_response_matrix.bin', exist=lexist)
    if (lexist .and. ((scf_cycle > 1) .or. pe_restart)) then
        call openfile('pe_response_matrix.bin', lu, 'old', 'unformatted')
        rewind(lu)
        if (chol) then
            read(lu) B
        else
            allocate(ipiv(nelems))
            read(lu) B, ipiv
        end if
        close(lu)
    else
        call response_matrix(B)
        if (pe_debug) then
            anorm = lansp('1', B, 'L')
            write(luout, '(/4x,a,f15.8)') '1-norm of response matrix B: ', anorm
        end if

        if (chol) then
            call pptrf(B, 'L', info)
            if (info /= 0) then
                write(luout,*) 'INFO: Cholesky factorization of classical&
                               & response matrix failed.'
                write(luout,*) 'INFO: Cholesky factorization can be disabled with&
                               & NOCHOL under .DIRECT.'
                stop 'ERROR: cannot create classical response matrix.'
            end if
        else
            allocate(ipiv(nelems))
            call sptrf(B, 'L', ipiv, info)
            if (info /= 0) then
                stop 'ERROR: cannot create classical response matrix.'
            end if
        end if

        if (pe_debug) then
            if (chol) then
                call ppcon(B, anorm, rcond, 'L')
                write(luout, '(4x,a,f15.8/)') 'Condition number of response matrix B: ', 1.0_dp / rcond
            else
                call spcon(B, ipiv, anorm, rcond, 'L')
                write(luout, '(4x,a,f15.8/)') 'Condition number of response matrix B: ', 1.0_dp / rcond
            end if
        end if

        call openfile('pe_response_matrix.bin', lu, 'unknown', 'unformatted')
        rewind(lu)
        if (chol) then
            write(lu) B
        else
            write(lu) B, ipiv
        end if
        close(lu)

    end if

    F = Fks

    if (chol) then
        call pptrs(B, F, 'L', info)
        if (info /= 0) then
            write(luout,*) 'INFO: Cholesky solver failed.'
            write(luout,*) 'INFO: Cholesky can be disabled with NOCHOL under .DIRECT.'
            stop 'ERROR: cannot solve for induced moments.'
        end if
        deallocate(B)
    else
        call sptrs(B, F, ipiv, 'L', info)
        if (info /= 0) then
            stop 'ERROR: cannot solve for induced moments.'
        end if
        deallocate(B, ipiv)
    end if

    Mkinds = F

    if (fock) then
        call openfile('pe_induced_moments.bin', lu, 'unknown', 'unformatted')
        rewind(lu)
        write(lu) Mkinds
        close(lu)
    end if

end subroutine direct_solver

!------------------------------------------------------------------------------

subroutine response_matrix(B)

    use pe_variables
    use pe_constants
    use pe_multipole_tools, only: damping_coefficient
    use pe_blas_interfaces, only: nrm2
    use pe_lapack_interfaces, only: pptrf, sptrf, pptri, sptri

! TODO: Cutoff radius

    real(dp), dimension(:), intent(out) :: B

    logical :: exclude
    integer :: info
    integer :: i, j, k, l, m
    integer, dimension(3) :: ipiv
    real(dp), parameter :: fourpi = 4.0_dp * pi
    real(dp), parameter :: i3 = 1.0_dp / 3.0_dp
    real(dp), parameter :: i6 = 1.0_dp / 6.0_dp
    real(dp) :: dism0
    real(dp) :: fE = 1.0_dp
    real(dp) :: fT = 1.0_dp
    real(dp) :: v, alpha_i, alpha_j
    real(dp) :: R, R3, R5, T
    real(dp), dimension(3) :: Rij
    real(dp), dimension(6) :: P11inv
    real(dp) :: eps_fac

    B = 0.0_dp

    m = 0
    do i = 1, nsites
        if (zeroalphas(i)) cycle
        P11inv = P11s(:,i)
        call pptrf(P11inv, 'L', info)
        if (info /= 0) then
            P11inv = P11s(:,i)
            call sptrf(P11inv, 'L', ipiv, info)
            if (info /= 0) then
                stop 'ERROR: could not factorize polarizability.'
            else if (chol) then
                chol = .false.
            end if
            call sptri(P11inv, ipiv, 'L')
        else
            call pptri(P11inv, 'L')
        end if
        if (pe_ind_damp) then
            alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
        end if
        do l = 3, 1, - 1
            do j = i, nsites
                if (zeroalphas(j)) cycle
                if (j == i) then
                    if (l == 3) then
                        do k = 1, l
                            B(m+k) = P11inv(k)
                        end do
                    else if (l == 2) then
                        do k = 1, l
                            B(m+k) = P11inv(3+k)
                        end do
                    else if (l == 1) then
                        do k = 1, l
                            B(m+k) = P11inv(5+k)
                        end do
                    end if
                    m = m + l
                else
                    if (pe_nomb) then
                        m = m + 3
                        cycle
                    end if
                    exclude = .false.
                    do k = 1, lexlst
                        if (exclists(k,i) == exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (exclude) then
                        m = m + 3
                        cycle
                    end if
                    Rij = Rs(:,j) - Rs(:,i)
                    R = nrm2(Rij)
                    R3 = R**3
                    R5 = R**5
! TODO: cutoff radius
!                        if (R > cutoff) then
!                            m = m + 3
!                            cycle
!                        end if
                    ! Thole damping
                    ! JPC A 102 (1998) 2399 and Mol. Sim. 32 (2006) 471
                    ! v = a * u, where a = 2.1304 (default) and u = R / (alpha_i * alpha_j)**(1/6)
                    ! fE = 1-(v^2/2+v+1)*exp(-v)
                    ! fT = 1-(v^3/6+v^2/2+v+1)*exp(-v)
                    if (pe_ind_damp) then
                        alpha_j = (P11s(1,j) + P11s(4,j) + P11s(6,j)) * i3
                        call damping_coefficient(Rij, alpha_i, alpha_j, ind_damp, v)
                        fE = 1.0_dp - (0.5_dp * v**2 + v + 1.0_dp) * exp(- v)
                        fT = fE - i6 * v**3 * exp(- v)
                    end if
                    if (l == 3) then
                        do k = 1, 3
                            T = 3.0_dp * Rij(1) * Rij(k) * fT / R5
                            if (k == 1) T = T - fE / R3
                            B(m+k) = - T
                        end do
                    else if (l == 2) then
                        do k = 1, 3
                            T = 3.0_dp * Rij(2) * Rij(k) * fT / R5
                            if (k == 2) T = T - fE / R3
                            B(m+k) = - T
                        end do
                    else if (l == 1) then
                        do k = 1, 3
                            T = 3.0_dp * Rij(3) * Rij(k) * fT / R5
                            if (k == 3) T = T - fE / R3
                            B(m+k) = - T
                        end do
                    end if
                    m = m + 3
                end if ! i /= j
            end do  ! do j = i, nsites
        end do  ! do l = 3, 1, - 1
    end do  ! do i = 1, nsites

    if (pe_debug) then
        do i = 1, (3 * npols) * (3 * npols + 1) / 2
            write (luout,*) 'Response matrix(i)',i, B(i)
        end do
    end if

end subroutine response_matrix

!------------------------------------------------------------------------------

end module pe_induced_moments
